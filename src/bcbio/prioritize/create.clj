(ns bcbio.prioritize.create
  "Create database of priority regions based on genes and existing biological evidence"
  (:require [bcbio.run.clhelp :as clhelp]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]))

(defn from-known
  "Create a new database grouped by bins with information from the known input files."
  [bin-file known-files out-file])

(defn- usage [options-summary]
  (->> ["Create a database of priority regions based on gene and domain regions with biological evidence:"
        ""
        "Usage: bcbio-prioritize createdb [options] -o output-db -b summary-bins -k known-input"
        ""
        "  output-db   : Output database to use in subsequent prioritization"
        "  summary-bins: BED file to guide binning of known variations. Can be transcript"
        "                or domain based"
        "  known-input : Known variations to apply to summary bins. Can be BED or VCF files and"
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
