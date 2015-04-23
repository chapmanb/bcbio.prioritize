(ns bcbio.prioritize.known
  "Prioritize a set of calls with a set of known regions of interest"
  (:require [bcbio.prioritize.utils :as utils]
            [bcbio.run.clhelp :as clhelp]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.data.csv :as csv]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.set :as cset]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]))

(defn- intersect
  "Intersect two BED files returning a tsv of overlaps."
  [a-file b-file base-name out-dir]
  (let [out-file (str (io/file out-dir (str base-name "-intersect.tsv")))]
    (itx/run-cmd out-file
                 "bedtools intersect -a ~{a-file} -b ~{b-file} -wa -wb "
                 "> ~{out-file}")))

(defn- combine-hits
  "Combine a set of binned hits into a short descriptive string about a call"
  [hits-plus-coords]
  (letfn [(merge-hit [coll hit]
            {:name (conj (get coll :name #{}) (:name hit))
             :origin (cset/union (get coll :support #{}) (set (map :origin (:support hit))))})
          (merge->str [coll]
            (format "%s:%s" (string/join "," (sort (:origin coll)))
                    (string/join "," (sort (:name coll)))))]
    (->> hits-plus-coords
         (map #(edn/read-string (last %)))
         (reduce merge-hit {})
         merge->str)))

(defn- summarize-matches [[coords hits]]
  (conj (vec coords) (combine-hits hits)))

(defn- summarize
  "Summarize an intersected bedtools TSV file into a prioritized bed"
  [in-file out-file-orig]
  (let [out-file (string/replace out-file-orig ".bed.gz" ".bed")]
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

(defn prioritize
  [input-file known-file out-file]
  (itx/with-named-tempdir [work-dir (str (fsp/file-root out-file) "-work")]
    (-> input-file
        (utils/bgzip-index work-dir)
        (intersect known-file "i2k" work-dir)
        (summarize out-file)
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
                  ["-i" "--input INPUT" "Input file of calls for prioritization (bed.gz "]
                  ["-k" "--known KNOWN" "Organized set of known regions from 'create' (bed.gz)"]
                  ["-h" "--help"]]
        {:keys [options arguments errors summary]} (parse-opts args opt-spec)
        missing (clhelp/check-missing options #{:output :input :known})]
    (cond
      (:help options) (clhelp/exit 0 (usage summary))
      errors (clhelp/exit 1 (clhelp/error-msg errors))
      (not (empty? missing)) (clhelp/exit 1 (str (clhelp/error-msg missing) \newline \newline (usage summary)))
      :else (prioritize (:input options) (:known options) (:output options)))))
