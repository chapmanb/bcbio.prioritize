(ns bcbio.prioritize.main
  (:require [clojure.java.io :as io]
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
  {:version {:main version
             :doc "Print version"}})

(defn -main [& args]
  (if-let [to-run (get progs (keyword (first args)))]
    (try
      (apply (:main to-run) (rest args))
      (catch Exception e
        (timbre/error e)
        (shutdown-agents)
        (System/exit 1))
      (finally
        (shutdown-agents)
        (System/exit 0)))
    (do
      (println "Prioritize small variants, structural variants and coverage based on biological inputs\n")
      (println "Commands:")
      (doseq [k (sort (keys progs))]
        (println (format "%-15s %s" (name k) (-> progs k :doc))))
      (System/exit 1))))
