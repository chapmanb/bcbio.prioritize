(ns bcbio.prioritize.provider.civic
  "Interaction with the CIViC API for retrieving clinical mutations
   https://civic.genome.wustl.edu/#/api-documentation"
  (:import [htsjdk.samtools.util BlockCompressedOutputStream])
  (:require [bcbio.prioritize.utils :as utils]
            [bcbio.run.clhelp :as clhelp]
            [bcbio.run.itx :as itx]
            [clj-http.client :as client]
            [clj-time.core]
            [clj-time.format]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [clojure.tools.cli :refer [parse-opts]]
            [slingshot.slingshot :refer [try+]]))

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
  (try+
   (url->json :civic (format "api/genes/%s" civic-id))
   (catch [:status 404] {})))

(defn genes
  "Retrieve all genes available in CIViC"
  []
  (url->json :civic "api/genes"))

(defn genes-hack
  "Retrieve genes via a hack on IDs, until /api/genes returns all items"
  []
  (->> (range 200)
       (map gene)
       (remove nil?)))

(defn variant-summary
  "Organize variants into a summary for a CIViC gene."
  [g]
  {:origin #{"civic"}
   :variant-groups (map :name (:variant_groups g))
   :url (format "https://civic.genome.wustl.edu//#/events/genes/%s/summary" (:id g))
   :drugs (->> (map :id (:variants g))
               (mapcat #(url->json :civic (format "/api/variants/%s/evidence_items" %)))
               (mapcat :drugs)
               (map :name)
               (remove (partial = "N/A"))
               set)})

(defn- gene->build
  "Retrieve human genome build information from CIViC gene object.
   Defaults to GRCh37 if not found."
  [g]
  (keyword (or (-> g :variant_groups first :variants first :coordinates :reference_build)
               "GRCh37")))

(defn gene->bed
  "Convert gene into a BED-ready output with metadata encoded in the name field."
  [options g]
  (let [g-info (merge (gene-info :human (:build options) (:name g))
                      (variant-summary g))]
    (-> (select-keys g-info [:assembly_name :seq_region_name :start :end])
        (assoc :name {:support (select-keys g-info [:url :variant-groups :origin :drugs])
                      :name #{(:display_name g-info)}}))))

(defn- gene-w-no-info?
  "Check if an input gene has no information"
  [g]
  (and (empty? (:variants g))
       (empty? (:variant_groups g))
       (empty? (:sources g))
       (nil? (:clinical_description g))))

(defn curdb->bed
  "Dump the current CIViC database to a BED file for prioritization."
  [options]
  (let [out-file (format "civic-%s-%s.bed.gz" (:build options)
                         (clj-time.format/unparse (clj-time.format/formatters :year-month-day)
                                                  (clj-time.core/now)))]
    (itx/with-tx-file [tx-out-file out-file]
      (with-open [wtr (io/writer (BlockCompressedOutputStream. (io/file tx-out-file)))]
        (doseq [g (->> (genes-hack)
                       (remove gene-w-no-info?)
                       (map (partial gene->bed options))
                       (sort-by (juxt :seq_region_name :start)))]
          (.write wtr (str (string/join "\t" (map (update-in g [:name] pr-str)
                                                  [:seq_region_name :start :end :name]))
                           "\n")))))
    (utils/bgzip-index out-file)))

(defn- usage [options-summary]
  (->> ["Create file of priority regions with genes from current CIViC database (https://civic.genome.wustl.edu)"
        ""
        "Usage: bcbio-prioritize create-civic"
        ""
        "Options:"
        options-summary]
       (string/join \newline)))

(defn -main [& args]
  (let [opt-spec [["-h" "--help"]
                  ["-b" "--build BUILD" "Genome build to use (defaults to GRCh37)" :default "GRCh37"]]
        {:keys [options summary errors]} (parse-opts args opt-spec)
        missing (clhelp/check-missing options #{})]
    (cond
      (:help options) (clhelp/exit 0 (usage summary))
      errors (clhelp/exit 1 (clhelp/error-msg errors))
      (not (empty? missing)) (clhelp/exit 1 (str (clhelp/error-msg missing) \newline \newline (usage summary)))
      :else (curdb->bed options))))
