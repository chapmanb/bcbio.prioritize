(ns bcbio.prioritize.create
  "Create database of priority regions based on genes and existing biological evidence"
  (:import [htsjdk.samtools.util BlockCompressedOutputStream]
           [htsjdk.tribble AbstractFeatureReader]
           [htsjdk.tribble.readers TabixReader]
           [htsjdk.tribble.bed BEDCodec])
  (:require [bcbio.run.clhelp :as clhelp]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [bcbio.prioritize.utils :as utils]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [gavagai.core :as gavagai]))

(defmulti hit->rec
  "Convert a hit into an annotation specific record"
  (fn [f hit]
    (cond
      (.contains (string/lower-case f) "cosmic") :cosmic)))

(defmethod hit->rec :cosmic
  [_ hit]
  {:origin "cosmic"
   :id (nth (string/split hit #"\t") 2)})

(defn- find-hits-one
  "Find elements from a single input in a region"
  [cur-bin known]
  (with-open [reader (TabixReader. known)]
    (let [q (.query reader (:chr cur-bin) (:start cur-bin) (:end cur-bin))]
      (vec (map (partial hit->rec known) (take-while (complement nil?) (repeatedly #(.next q))))))))

(defn- find-hits-many
  "Retrieve hits found in the known inputs files by region."
  [cur-bin knowns]
  (mapcat (partial find-hits-one cur-bin) knowns))

(defn- bin->known-bed
  "Convert a bin region into a BED compatible output file with output information."
  [prep-known cur-bin]
  (when-let [support (seq (find-hits-many cur-bin prep-known))]
    (assoc cur-bin :name
           (pr-str {:support support
                    :name (:name cur-bin)}))))

(defn from-known
  "Create a new database grouped by bins with information from the known input files."
  [bin-file known-files out-file]
  (itx/with-named-tempdir [work-dir (str (fsp/file-root out-file) "-work")]
    (itx/with-tx-file [tx-out-file out-file]
      (let [prep-bin (utils/bgzip-index bin-file work-dir)
            prep-known (map #(utils/bgzip-index % work-dir) known-files)]
        (with-open [bin-reader (AbstractFeatureReader/getFeatureReader prep-bin (BEDCodec.) false)
                    wtr (io/writer (BlockCompressedOutputStream. tx-out-file))]
          (gavagai/with-translator (gavagai/register-converters
                                    [["htsjdk.tribble.bed.FullBEDFeature" :only [:chr :start :end :name]]])
            (doseq [cur-bin (->> (.iterator bin-reader)
                                 (map gavagai/translate)
                                 (map (partial bin->known-bed prep-known))
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
      missing (clhelp/exit 1 (str (clhelp/error-msg missing) \newline \newline (usage summary)))
      :else (from-known (:bins options) (:known options) (:output options)))))
