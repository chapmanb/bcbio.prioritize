(ns bcbio.prioritize.create
  "Create database of priority regions based on genes and existing biological evidence"
  (:import [htsjdk.samtools.util BlockCompressedOutputStream BlockCompressedInputStream]
           [htsjdk.tribble AbstractFeatureReader]
           [htsjdk.variant.vcf VCFCodec]
           [htsjdk.tribble.readers TabixReader LineIteratorImpl LineReaderUtil]
           [htsjdk.tribble.bed BEDCodec])
  (:require [bcbio.run.clhelp :as clhelp]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.prioritize.provider.civic :as civic]
            [bcbio.prioritize.provider.clinvar :as clinvar]
            [bcbio.prioritize.provider.intogen :as intogen]
            [bcbio.prioritize.provider.oncomine :as oncomine]
            [bcbio.prioritize.utils :as utils]
            [clojure.set :as cset]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [gavagai.core :as gavagai]))

(defmulti hit->rec
  "Convert a VCF hit into an annotation specific record"
  (fn [f hit]
    (cond
      (.contains (string/lower-case f) "clinvar") :clinvar
      (.contains (string/lower-case f) "cosmic") :cosmic
      (.contains (string/lower-case f) "oncomine") :oncomine
      (.contains (string/lower-case f) "clinvar") :clinvar)))

(defmethod hit->rec :cosmic
  [_ hit]
  {:origin "cosmic"
   :id (.getID hit)})

(defmethod hit->rec :oncomine
  [_ hit]
  (oncomine/vc->rec hit))

(defmethod hit->rec :clinvar
  [_ hit]
  (clinvar/vc->rec hit))

(defmethod hit->rec :civic
  [_ hit]
  (civic/hit->rec hit))

(defmethod hit->rec :default
  [f hit]
  (throw (Exception. (str "Need to implement method to convert input file to records: " f))))

(defn- find-hits-one
  "Find elements from a single VCF input in a region."
  [known]
  (fn [cur-bin]
    (with-open [reader (TabixReader. known)
                liter (LineIteratorImpl. (LineReaderUtil/fromBufferedStream
                                          (BlockCompressedInputStream. (io/file known))))]
      (let [q (.query reader (:chr cur-bin) (:start cur-bin) (:end cur-bin))
            codec (doto (case (last (fsp/split-ext+ known))
                          ".vcf.gz" (VCFCodec.)
                          ".bed.gz" (BEDCodec.))
                    (.readActualHeader liter))]
        (->> (take-while (complement nil?) (repeatedly #(.next q)))
             (map #(.decode codec %))
             (map (partial hit->rec known))
             vec)))))

(defn- bin->known-bed
  "Convert a bin region into a BED compatible output file with output information."
  [get-known-fns]
  (fn [cur-bin]
    (when-let [support (seq (mapcat #(% cur-bin) get-known-fns))]
      (assoc cur-bin :name
             (pr-str {:support (->> support
                                    (map #(reduce-kv (fn [m k v] (assoc m k (if (set? v) v #{v}))) {} %))
                                    (apply merge-with cset/union))
                      :name (set (string/split (:name cur-bin) #","))})))))

(defmulti get-known
  "Create a function to return known supporting information based on inputs."
  (fn [known-file work-dir]
    (letfn [(is-intogen? [known-file]
              (.contains known-file "intogen_cancer_drivers"))
            (is-oncomine? [known-file]
              (.contains (string/lower-case known-file) "oncomine"))
            (is-civic? [known-file]
              (.contains (string/lower-case known-file) "civic"))
            (is-bed? [known-file]
              (or (.endsWith known-file ".bed.gz") (.endsWith known-file ".bed")))
            (is-vcf? [known-file]
              (or (.endsWith known-file ".vcf.gz") (.endsWith known-file ".vcf")))]
      (cond
        (and (is-vcf? known-file) (is-oncomine? known-file)) :oncomine
        (is-vcf? known-file) :vcf
        (is-bed? known-file) :bed
        (and (is-bed? known-file) (is-civic? known-file)) :bed
        (is-intogen? known-file) :intogen))))

(defmethod get-known :vcf
  ^{:doc "Retrieve from an input VCF files, handling COSMIC."}
  [known-file work-dir]
  (find-hits-one (utils/bgzip-index known-file work-dir)))

(defmethod get-known :oncomine
  ^{:doc "Retrieve from Oncomine VCF file"}
  [known-file work-dir]
  (-> known-file
      (oncomine/clean-file work-dir)
      (utils/bgzip-index work-dir)
      find-hits-one))

(defmethod get-known :intogen
  ^{:doc "Retrieve from a directory download from IntoGen."}
  [known-file _]
  (intogen/get-known known-file))

(defmethod get-known :bed
  ^{:doc "Parse and return a standard BED file, handling CIViC preparation."}
  [known-file work-dir]
  (find-hits-one (utils/bgzip-index known-file work-dir)))

(defmethod get-known :default
  [known-file _]
  (throw (Exception. (str "Do not know how to prepare input file: " known-file))))

(defn from-known
  "Create a new database grouped by bins with information from the known input files."
  [bin-file known-files out-file]
  (itx/with-named-tempdir [work-dir (str (fsp/file-root out-file) "-work")]
    (itx/with-tx-file [tx-out-file out-file]
      (let [prep-bin (-> bin-file
                         (utils/merge-bed work-dir)
                         (utils/bgzip-index work-dir))
            get-known-fn (bin->known-bed (map #(get-known % work-dir) known-files))]
      (with-open [bin-reader (AbstractFeatureReader/getFeatureReader prep-bin (BEDCodec.) false)
                  wtr (io/writer (BlockCompressedOutputStream. tx-out-file))]
        (gavagai/with-translator (gavagai/register-converters
                                  [["htsjdk.tribble.bed.FullBEDFeature" :only [:chr :start :end :name]]])
          (doseq [cur-bin (->> (.iterator bin-reader)
                               (map gavagai/translate)
                               (map get-known-fn)
                               (remove nil?))]
            (.write wtr (str (string/join "\t" (map cur-bin [:chr :start :end :name]))
                             "\n"))))))))
  (utils/bgzip-index out-file))

(defn- usage [options-summary]
  (->> ["Create file of priority regions based on gene and domain regions with biological evidence:"
        ""
        "Usage: bcbio-prioritize create [options] -o output -b summary-bins -k known-input"
        ""
        "  output   : Output to use in subsequent prioritization (a bed.gz file)"
        "  summary-bins: BED file to guide binning of known variations. Can be transcript"
        "                or domain based"
        "  known-input : Known variations to apply to summary bins. Can be BED or VCF file, and"
        "                specified multiple times"
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [opt-spec [["-o" "--output OUTPUT" "Output file to write results to"]
                  ["-b" "--bins BINS" "BED file to bin regions to evaluate in: transcripts or domains"]
                  ["-k" "--known KNOWN" "File of known changes to associate. Can specify multiple times"
                   :assoc-fn (fn [m k v] (update-in m [k] (fnil conj #{}) v))]
                  ["-h" "--help"]]
        {:keys [options arguments errors summary]} (parse-opts args opt-spec)
        missing (clhelp/check-missing options #{:output :bins :known})]
    (cond
      (:help options) (clhelp/exit 0 (usage summary))
      errors (clhelp/exit 1 (clhelp/error-msg errors))
      (not (empty? missing)) (clhelp/exit 1 (str (clhelp/error-msg missing) \newline \newline (usage summary)))
      :else (from-known (:bins options) (:known options) (:output options)))))
