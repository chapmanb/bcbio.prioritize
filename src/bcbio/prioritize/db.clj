(ns bcbio.prioritize.db
  "Initialize and query a database of priors for variant prioritization"
  (:require [clojure.java.io :as io]
            [mount.core :refer [defstate]]
            [config.core :refer [load-env]]
            [hodur-engine.core :as hodur]
            [hodur-datomic-schema.core :as hodur-datomic]
            [io.rkn.conformity :as c]
            [datomic.api :as d]))

(defstate env
  "Environmental variables, defined by config.edn in resources"
  :start (load-env))

(defn- init
  "Transactionally initialize database with schema"
  []
  (let [v-schema (-> "schema/variants.edn"
                     io/resource
                     hodur/init-path
                     hodur-datomic/schema)
        uri (get-in env [:datomic :uri])]
    (d/create-database uri)
    (let [conn (d/connect uri)]
      (c/ensure-conforms conn {:variants {:txes [v-schema]}} [:variants])
      conn)))

(defn- close [conn]
  (let [uri (get-in env [:datomic :uri])]
    (.release conn)
    ;(d/delete-database uri)
    ))

(defstate conn
  "Create a database connection to configured datomic db."
  :start (init)
  :stop (close conn))

(defn get-sample
  "Retrieve an existing sample from the database by name, or create a new one"
  [db name]
  (if-let [s (ffirst (d/q '[:find (pull ?s [*])
                            :in $ ?name
                            :where [?s :sample/name ?name]]
                          db name))]
    s
    {:db/id (d/tempid :db.part/user)
     :sample/name name}))

(defn- prepare-evidence-data
  "Prepare data input for supplied evidence"
  [e g db]
  (cond-> {:db/id (d/tempid :db.part/user)
           :evidence/origin (:origin e)
           :evidence/variant-type (:variant e)
           :evidence/gene (:db/id g)
           :evidence/diseases (or (:diseases e) [])
           :evidence/drugs (or (:drugs e) [])}
    (:sample e) (assoc :evidence/sample (get-sample db (:sample e)))))

(defn- prepare-location-data
  "Prepare data input for map with location"
  [l]
  {:db/id (d/tempid :db.part/user)
   :location/build (:assembly_name l)
   :location/contig (:seq_region_name l)
   :location/start (:start l)
   :location/end (:end l)})

(defn get-gene
  "Retrieve an existing gene from the database"
  [db g-in]
  (ffirst (d/q '[:find (pull ?g [*])
                 :in $ ?name
                 :where [?g :gene/name ?name]]
               db g-in)))

(defn get-evidence
  [db vtype]
  (first (d/q '[:find (pull ?e [*])
                :in $ ?type
                :where [?e :evidence/variant-type ?type]]
              db vtype)))

(defn get-gene-with-location
  "Retrieve a gene, potentially existing and add location information if missing."
  [db g-in loc-fn]
  (if-let [g (get-gene db g-in)]
    {:gene g :locs []}
    (let [locs (map prepare-location-data (loc-fn g-in))
          g {:db/id (d/tempid :db.part/user)
             :gene/name g-in
             :gene/location (map :db/id locs)}]
      {:gene g :locs locs})))

(defn- prepare-evidence
  "Prepare gene, updating existing with additional evidence or creating new."
  [g-in loc-fn db]
  (let [g-info (get-gene-with-location db (:name g-in) loc-fn)
        evs (map #(prepare-evidence-data % (:gene g-info) db) (:evidence g-in))]
    (if (:locs g-info)
      (concat evs [(:gene g-info)] (:locs g-info))
      evs)))

(defn load-data
  "Load gene evidence data into prepared database."
  [conn gs loc-fn]
  (c/ensure-conforms conn {:d {:txes [(mapcat #(prepare-evidence % loc-fn (d/db conn)) gs)]}} [:d]))
