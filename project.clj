(defproject bcbio.prioritize "0.0.9-SNAPSHOT"
  :description "Prioritize small variants, structural variants and coverage based on biological inputs"
  :url "https://github.com/chapmanb/bcbio.prioritize"
  :license {:name "MIT" :url "http://www.opensource.org/licenses/mit-license.html"}
  :dependencies [[org.clojure/clojure "1.10.0"]
                 [org.clojure/algo.generic "0.1.2"]
                 [org.clojure/data.csv "0.1.2"]
                 [semantic-csv "0.1.0"]
                 [org.clojure/tools.cli "0.3.1"]
                 [clj-http "3.9.1"]
                 [cheshire "5.8.1"]
                 [clj-time "0.15.1"]
                 [slingshot "0.12.2"]
                 [gavagai "0.3.2"]
                 [bcbio.variation.recall "0.2.0"]
                 [com.github.samtools/htsjdk "1.140"]
                 [bcbio.run "0.0.6"]
                 [dk.ative/docjure "1.13.0"]
                 [com.datomic/datomic-free "0.9.5697"]
                 [io.rkn/conformity "0.5.1" :exclusions [com.datomic/datomic-free]]
                 [hodur/engine "0.1.5"]
                 [hodur/datomic-schema "0.1.0"]]
  :profiles {:uberjar {:aot [bcbio.prioritize.main]}}
  :main ^:skip-aot bcbio.prioritize.main)
