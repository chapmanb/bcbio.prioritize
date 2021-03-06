(ns bcbio.prioritize.main
  (:require [bcbio.prioritize.create :as create]
            [bcbio.prioritize.known :as known]
            [bcbio.prioritize.missing :as missing]
            [bcbio.prioritize.provider.civic :as civic]
            [clojure.java.io :as io]
            [taoensso.timbre :as timbre])
  (:gen-class))

(defn version
  []
  (let [version (with-open [reader (-> "META-INF/maven/bcbio.prioritize/bcbio.prioritize/pom.properties"
                                       io/resource
                                       io/reader)]
                  (-> (doto (java.util.Properties.)
                        (.load reader))
                      (.getProperty "version")))]
    (println "bcbio.prioritize" version)))

(def ^{:private true} progs
  {:create {:main create/-main
            :doc "Create file of priority regions based on genes and existing biological evidence"}
   :known {:main known/-main
           :doc "Prioritize a set of calls with based on binned regions of interest"}
   :missing {:main missing/-main
             :doc "Identify regions with missing coverage based on regions of interest"}
   :create-civic {:main civic/-main
                  :doc "Create file of priority regions from current CIViC database"}
   :version {:main version
             :doc "Print version"}})

(defn -main [& args]
  (if-let [to-run (get progs (keyword (first args)))]
    (do
      (try
        (apply (:main to-run) (rest args))
        (catch Exception e
          (timbre/error e)
          (shutdown-agents)
          (System/exit 1)))
      (shutdown-agents)
      (System/exit 0))
    (do
      (println "Prioritize small variants, structural variants and coverage based on biological inputs\n")
      (println "Commands:")
      (doseq [k (sort (keys progs))]
        (println (format "%-15s %s" (name k) (-> progs k :doc))))
      (System/exit 1))))
