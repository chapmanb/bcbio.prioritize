(defproject bcbio.prioritize "0.0.3"
  :description "Prioritize small variants, structural variants and coverage based on biological inputs"
  :url "https://github.com/chapmanb/bcbio.prioritize"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.7.0"]
                 [org.clojure/algo.generic "0.1.2"]
                 [org.clojure/data.csv "0.1.2"]
                 [semantic-csv "0.1.0"]
                 [org.clojure/tools.cli "0.3.1"]
                 [clj-biosequence "0.2.6"]
                 [clj-http "1.1.1"]
                 [clj-time "0.11.0"]
                 [slingshot "0.12.2"]
                 [com.taoensso/timbre "3.4.0"]
                 [gavagai "0.3.2"]
                 [bcbio.variation.recall "0.1.4"]
                 [com.github.samtools/htsjdk "1.140"]
                 [bcbio.run "0.0.5"]]
  :profiles {:uberjar {:aot [bcbio.prioritize.main]}}
  :main ^:skip-aot bcbio.prioritize.main)
