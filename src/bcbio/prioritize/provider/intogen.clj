(ns bcbio.prioritize.provider.intogen
  "Prioritization inputs using IntoGen http://www.intogen.org/downloads
   Tested on release 2014.12"
  (:require [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [semantic-csv.core :as sc]))

(defn- read-driver-file
  "Create list of driver mutations along with type and role."
  [in-file]
  (with-open [rdr (io/reader in-file)]
    (->> (csv/read-csv rdr :separator \tab)
         (sc/remove-comments)
         (sc/mappify)
         (map #(dissoc % :OncodriveROLE_prob))
         doall)))

(defn- read-tumor-type
  "Create map of gene symbols to list of tumor types"
  [in-file]
  (letfn [(type-by-gene [coll cur]
            (assoc coll (:geneHGNCsymbol cur)
                   (conj (get coll (:geneHGNCsymbol cur) #{}) (:Tumor_type cur))))]
    (with-open [rdr (io/reader in-file)]
      (->> (csv/read-csv rdr :separator \tab)
           (sc/remove-comments)
           (sc/mappify)
           (reduce type-by-gene {})
           doall))))

(defn get-known
  "Retrieve known IntoGen drivers mapped by gene name."
  [intogen-dir]
  (let [tumor-types (read-tumor-type (str (io/file intogen-dir "Mutational_drivers_per_tumor_type.tsv")))
        drivers (->> (str (io/file intogen-dir "Drivers_type_role.tsv"))
                     read-driver-file
                     (map #(assoc % :Tumor_type (get tumor-types (:geneHGNCsymbol %) #{})))
                     (map (juxt :geneHGNCsymbol identity))
                     (into {})
                     (#(assoc % :origin #{"oncomine"})))]
    (fn [hit]
      (remove nil? (map #(get drivers %) (string/split (:name hit) #","))))))
