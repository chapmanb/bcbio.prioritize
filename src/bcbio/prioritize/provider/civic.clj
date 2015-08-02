(ns bcbio.prioritize.provider.civic
  "Interaction with the CIViC API for retrieving clinical mutations
   https://civic.genome.wustl.edu/#/api-documentation"
  (:require [clj-http.client :as client]))

;; ## Interaction with prioritzation

(defn hit->rec
  "Parse out BED name field for passing into build of binned priority targets."
  [hit]
  (println hit))

;; ## Preparing CIViC BED file

(defn- url->json*
  [base url]
  (-> (str base "/" url)
      (client/get  {:as :json
                    :insecure? true})
      :body))

(defmulti url->json
  "Dispatch to multiple locations parsing JSON output"
  (fn [loc & _]
    loc))

(defmethod url->json :civic
  [_ url]
  (url->json* "https://civic.genome.wustl.edu" url))

(defmethod url->json :ensembl-GRCh37
  [_ url]
  (url->json* "http://grch37.rest.ensembl.org" url))

(defmethod url->json :ensembl
  [_ url]
  (url->json* "http://rest.ensembl.org" url))

(defn gene-info
  "Fetch gene information from Ensembl for given display name from CIViC"
  [species build gene-name]
  (let [loc (if (= build :GRCh37) :ensembl-GRCh37 :ensembl)
        species-url (get {:human "homo_sapiens"} species species)]
    (url->json loc (format "lookup/symbol/%s/%s?content-type=application/json" species-url gene-name))))

(defn gene
  "Retrieve a single gene from CIViC"
  [civic-id]
  (url->json :civic (format "api/genes/%s" civic-id)))

(defn genes
  "Retrieve all genes available in CIViC"
  []
  (url->json :civic "api/genes"))

(defn variant-summary
  "Organize variants into a summary for a CIViC gene."
  [g]
  {:origin "civic"
   :variant-groups (map :name (:variant_groups g))
   :url (format "https://civic.genome.wustl.edu//#/events/genes/%s/summary" (:id g))
   :drugs (->> (map :id (:variants g))
               (mapcat #(url->json :civic (format "/api/variants/%s/evidence_items" %)))
               (mapcat :drugs)
               (map :name)
               set)})

(defn- gene->build
  "Retrieve human genome build information from CIViC gene object.
   Defaults to GRCh37 if not found."
  [g]
  (keyword (or (-> g :variant_groups first :variants first :coordinates :reference_build)
               "GRCh37")))

(defn gene->bed
  "Convert gene into a BED-ready output with metadata encoded in the name field."
  [g]
  (let [g-info (merge (gene-info :human (gene->build g) (:name g))
                      (variant-summary g))]
    (-> (select-keys g-info [:assembly_name :seq_region_name :start :end])
        (assoc :name (select-keys g-info [:url :variant-groups :origin :drugs :display_name])))))
