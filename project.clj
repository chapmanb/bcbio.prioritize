(defproject bcbio.prioritize "0.0.1-SNAPSHOT"
  :description "Prioritize small variants, structural variants and coverage based on biological inputs"
  :url "https://github.com/chapmanb/bcbio.prioritize"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [com.taoensso/timbre "3.4.0"]
                 [bcbio.run "0.0.1"]]
  :profiles {:uberjar {:aot [bcbio.prioritize.main]}}
  :main ^:skip-aot bcbio.prioritize.main)
