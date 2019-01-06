(ns bcbio.prioritize.db
  "Initialize and query a database of priors for variant prioritization"
  (:require [clojure.java.io :as io]
            [hodur-engine.core :as hodur]
            [hodur-datomic-schema.core :as hodur-datomic]
            [io.rkn.conformity :as c]
            [datomic.api :as d]))

(defn init
  "Transactionally initialize database with schema and seed data"
  []
  (let [v-schema (-> "schema/variants.edn"
                     io/resource
                     hodur/init-path
                     hodur-datomic/schema)
        uri "datomic:mem://variants"]
    (d/create-database uri)
    (let [conn (d/connect uri)]
      (c/ensure-conforms conn {:variants {:txes [v-schema]}} [:variants])
      conn)))

(defn- prepare-evidence-data
  "Prepare data input for supplied evidence"
  [e]
  {:db/id (d/tempid :db.part/user)
   :evidence/origin (:origin e)
   :evidence/diseases (or (:diseases e) [])
   :evidence/drugs (or (:drugs e) [])
   :evidence/variant-type (:variant e)})

(defn- prepare-location-data
  "Prepare data input for map with location"
  [l]
  {:db/id (d/tempid :db.part/user)
   :location/build (:assembly_name l)
   :location/contig (:seq_region_name l)
   :location/start (:start l)
   :location/end (:end l)})

(defn get-gene
  "Retrieve an existing gene from the database, allowing adding additional evidence."
  [db g-in]
  (ffirst (d/q '[:find (pull ?g [*])
                 :in $ ?name
                 :where [?g :gene/name ?name]]
               db g-in)))

(defn- prepare-gene-data
  "Prepare gene, updating existing with additional evidence or creating new."
  [g-in loc-fn db]
  (let [evs (map prepare-evidence-data (:evidence g-in))]
    (if-let [g (get-gene db g-in)]
      (concat [(update-in g [:gene/evidence] #(concat % evs))] evs)
      (let [locs (map prepare-location-data (loc-fn (:name g-in)))
            g {:db/id (d/tempid :db.part/user)
               :gene/name (:name g-in)
               :gene/location (map :db/id locs)
               :gene/evidence (map :db/id evs)}]
        (concat [g] locs evs)))))

(defn load-data
  "Load gene evidence data into prepared database."
  [conn gs loc-fn]
  (let [db (d/db conn)]
    (c/ensure-conforms conn {:d {:txes [(mapcat #(prepare-gene-data % loc-fn db) gs)]}} [:d])
    db))
