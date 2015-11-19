(ns bcbio.prioritize.known
  "Prioritize a set of calls with a set of known regions of interest"
  (:import  [htsjdk.variant.vcf VCFInfoHeaderLine VCFHeaderLineType VCFHeaderLineCount])
  (:require [bcbio.prioritize.utils :as utils]
            [bcbio.run.clhelp :as clhelp]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.variation.variantcontext :as vc]
            [clojure.algo.generic.functor :refer [fmap]]
            [clojure.data.csv :as csv]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.set :as cset]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]))

(def ^{:doc "We avoid flattening all intermediate genes and look at ends only for larger events."}
  LARGE_EVENT_SIZE 50000)

(defn- intersect
  "Intersect two BED files returning a tsv of overlaps."
  [a-file b-file base-name out-dir]
  (let [out-file (str (io/file out-dir (str base-name "-intersect.tsv")))
        cat-cmd (if (.endsWith a-file ".gz") "zcat" "cat")]
    (itx/run-cmd out-file
                 "~{cat-cmd} ~{a-file} | bedtools intersect -a - -b ~{b-file} -wa -wb "
                 "> ~{out-file}")))

(defmulti get-start-end
  "Parse start/end coordinates from flat bedtools outputs"
  (fn [hit]
    (cond
      (> (count hit) 12) :vcf
      :else :bed)))

(defmethod get-start-end :bed
  [hit]
  (let [[contig start end] (take 3 hit)]
    [contig (Integer/parseInt start) contig (Integer/parseInt end)]))

(defmethod get-start-end :vcf
  [hit]
  (let [contig (str (first hit))
        start (Integer/parseInt (second hit))
        info (nth hit 7)
        end-info-tag (first (filter #(.startsWith % "END=") (string/split info #";")))
        end (when end-info-tag (Integer/parseInt (last (string/split end-info-tag #"="))))]
    (when-not (nil? end)
      [contig start contig end])))

(defn- hit-near-end?
  [[contig1 pos1 contig2 pos2] hit]
  (let [buffer 5000
        [hit-contig hit-start-str hit-end-str] (->> hit (take-last 4) drop-last)
        hit-start (Integer/parseInt hit-start-str)
        hit-end (Integer/parseInt hit-end-str)]
    (letfn [(in-end? [contig pos]
              (and (= hit-contig contig)
                   (>= pos (- hit-start buffer))
                   (<= pos (+ hit-end buffer))))]
      (or (in-end? contig1 pos1) (in-end? contig2 pos2)))))

(defn- limit-hits
  "For long SV events, limit hits to those that overlap endpoints."
  [hits]
  (let [[contig1 start contig2 end] (get-start-end (first hits))]
    (if (or (nil? contig1)
            (and (= contig1 contig2) (< (- end start) LARGE_EVENT_SIZE)))
      hits
      (filter (partial hit-near-end? [contig1 start contig2 end]) hits))))

(defn- combine-hits
  "Combine a set of binned hits into a short descriptive string about a call"
  [hits-plus-coords]
  (letfn [(parse-hit [hit-txt]
            ;; Parse hit, handling EDN and plain text cases
            (if (.startsWith hit-txt "{")
              (edn/read-string hit-txt)
              hit-txt))
          (prep-hit [hit]
            ;; Handle raw inputs where we convert a string into the name
            (if (every? nil? (map #(get hit %) [:origin :name]))
              {:name (set (string/split hit #","))}
              hit))
          (merge-hit [coll hit]
            {:name (cset/union (get coll :name #{}) (:name hit))
             :origin (cset/union (get coll :support #{}) (set (or (get-in hit [:support :origin])
                                                                  (map :origin (:support hit)))))})
          (merge->str [coll]
            (let [origin (string/join "," (sort (:origin coll)))
                  names (string/join "," (sort (:name coll)))]
              (->> [origin names]
                   (remove empty?)
                   (string/join ":"))))]
    (->> hits-plus-coords
         limit-hits
         (map #(parse-hit (last %)))
         (map prep-hit)
         (reduce merge-hit {})
         merge->str)))

(defn- summarize-matches [[coords hits]]
  (conj (vec coords) (combine-hits hits)))

(defmulti summarize
  "Handle summarization of multiple output formats"
  (fn [in-file _ out-file]
    (case (last (fsp/split-ext+ out-file))
      (".bed.gz" ".bed") :bed
      (".vcf.gz" ".vcf") :vcf)))

(defmethod summarize :bed
  ^{:doc "Summarize an intersected bedtools TSV file into a prioritized bed."}
  [in-file _ out-file-orig]
  (let [out-file (string/replace out-file-orig ".gz" "")]
    (itx/with-tx-file [tx-out-file out-file]
      (with-open [rdr (io/reader in-file)
                  wtr (io/writer out-file)]
        (as-> rdr $
          (csv/read-csv $ :separator \tab)
          (group-by (partial take 3) $)
          (map summarize-matches $)
          (sort-by (fn [[c s & xs]] [c (Integer/parseInt s)]) $)
          (csv/write-csv wtr $ :separator \tab))))
    out-file))

(defn- parse-intersects
  "Retrieve a lookup map of VCF coordinates to prioritization hits"
  [in-file]
  (with-open [rdr (io/reader in-file)]
    (as-> rdr $
      (csv/read-csv $ :separator \tab)
      (group-by (partial take 4) $)
      (fmap combine-hits $))))

(defn- vc-add-hit
  "Add hit information to a variant context if it passes."
  [hits vc]
  (when-let [hit (get hits [(:chr vc) (str (:start vc)) (:id vc) (.getBaseString (:ref-allele vc))])]
    (vc/vc-add-attr (:vc vc) "KNOWN" hit)))

(defmethod summarize :vcf
  ^{:doc "Summarize intersected bedtools TSV into an output VCF."}
  [in-file orig-file out-file]
  (let [hits (parse-intersects in-file)]
    (with-open [rdr (io/reader in-file)
                vcf-iter (vc/get-vcf-iterator orig-file)]
      (vc/write-vcf-w-template orig-file {:out out-file}
                               (->> (vc/parse-vcf vcf-iter)
                                    (map (partial vc-add-hit hits))
                                    (remove nil?))
                               :new-md #{(VCFInfoHeaderLine. "KNOWN" VCFHeaderLineCount/UNBOUNDED
                                                             VCFHeaderLineType/String
                                                             "Known gene associations from existing databases")})))
  out-file)

(defn prioritize
  [input-file known-file out-file]
  (itx/with-named-tempdir [work-dir (str (fsp/file-root out-file) "-work")]
    (-> input-file
        (utils/bgzip-index work-dir)
        (intersect known-file "i2k" work-dir)
        (summarize input-file out-file)
        utils/bgzip-index))
  (when (.endsWith out-file ".bed.gz")
    (fsp/remove-path (string/replace out-file ".bed.gz" ".bed")))
  out-file)

(defn- usage [options-summary]
  (->> ["Prioritize a set of calls with based on binned regions of interest"
        ""
        "Usage: bcbio-prioritize known [options] -i input -k known -o output"
        ""
        "  input:  File of calls to prioritize (bed or bed.gz)"
        "  known:  Prepared file of known regions to prioritize on"
        "  output: Output with prioritized calls(a bed.gz file)"
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [opt-spec [["-o" "--output OUTPUT" "Output file to write priority regions with annotations (bed.gz)"]
                  ["-i" "--input INPUT" "Input file of calls for prioritization (bed.gz or vcf.gz)"]
                  ["-k" "--known KNOWN" "Organized set of known regions from 'create' (bed.gz or vcf.gz)"]
                  ["-h" "--help"]]
        {:keys [options arguments errors summary]} (parse-opts args opt-spec)
        missing (clhelp/check-missing options #{:output :input :known})]
    (cond
      (:help options) (clhelp/exit 0 (usage summary))
      errors (clhelp/exit 1 (clhelp/error-msg errors))
      (not (empty? missing)) (clhelp/exit 1 (str (clhelp/error-msg missing) \newline \newline (usage summary)))
      :else (prioritize (:input options) (:known options) (:output options)))))
