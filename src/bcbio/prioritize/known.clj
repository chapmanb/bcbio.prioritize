(ns bcbio.prioritize.known
  "Prioritize a set of calls with a set of known regions of interest"
  (:require [bcbio.run.clhelp :as clhelp]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]))

(defn prioritize
  [input-file known-file out-file]
  (println input-file known-file out-file))

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
      missing (clhelp/exit 1 (str (clhelp/error-msg missing) \newline \newline (usage summary)))
      :else (prioritize (:input options) (:known options) (:output options)))))
