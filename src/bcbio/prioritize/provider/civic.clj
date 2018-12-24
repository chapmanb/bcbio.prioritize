(ns bcbio.prioritize.provider.civic
  "Interaction with the CIViC API for retrieving clinical mutations
   https://griffithlab.github.io/civic-api-docs/"
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
  (url->json* "https://civicdb.org" url))

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
    (try+
     (url->json loc (format "lookup/symbol/%s/%s?content-type=application/json" species-url gene-name))
     (catch [:status 400] {}))))

(defn gene
  "Retrieve a single gene from CIViC"
  [civic-id]
  (try+
   (url->json :civic (format "api/genes/%s" civic-id))
   (catch [:status 404] {})))

(defn genes
  "Retrieve all genes available in CIViC"
  []
  (:records (url->json :civic "api/genes?count=500")))

(defn genes-hack
  "Retrieve genes via a hack on IDs, until /api/genes returns all items"
  []
  (->> (range 200)
       (map gene)
       (remove nil?)))

(defn variant-summary
  "Organize variants into a summary for a CIViC gene."
  [g]
  (let [evidence-groups (->> (map :id (:variants g))
                             (map #(url->json :civic (format "/api/variants/%s" %)))
                             (mapcat :evidence_items))]
    {:origin #{"civic"}
     :url (format "https://civicdb.org//#/events/genes/%s/summary" (:id g))
     :diseases (->> (map :disease evidence-groups)
                    (map :name)
                    (remove (partial = "N/A"))
                    set)
     :drugs (->> (mapcat :drugs evidence-groups)
                 (map :name)
                 (remove (partial = "N/A"))
                 set)}))

(defn- gene->build
  "Retrieve human genome build information from CIViC gene object.
   Defaults to GRCh37 if not found."
  [g]
  (keyword (or (-> g :variant_groups first :variants first :coordinates :reference_build)
               "GRCh37")))

(defn gene->bed
  "Convert gene into a BED-ready output with metadata encoded in the name field."
  [options g]
  (println (:name g))
  (let [g-info (merge (gene-info :human (:build options) (:name g))
                      (variant-summary g))]
    (when (:assembly_name g-info)
      (-> (select-keys g-info [:assembly_name :seq_region_name :start :end])
          (assoc :name {:support (select-keys g-info [:url :diseases :origin :drugs])
                        :name #{(:display_name g-info)}})))))

(defn- gene-w-no-info?
  "Check if an input gene has no information"
  [g]
  (and (empty? (:variants g))
       (nil? (:clinical_description g))))

(defn curdb->bed
  "Dump the current CIViC database to a BED file for prioritization."
  [options]
  (let [out-file (format "civic-%s-%s.bed.gz" (:build options)
                         (clj-time.format/unparse (clj-time.format/formatters :year-month-day)
                                                  (clj-time.core/now)))]
    (itx/with-tx-file [tx-out-file out-file]
      (with-open [wtr (io/writer (BlockCompressedOutputStream. (io/file tx-out-file)))]
        (doseq [g (->> (genes)
                       (remove gene-w-no-info?)
                       (map (partial gene->bed options))
                       (sort-by (juxt :seq_region_name :start)))]
          (.write wtr (str (string/join "\t" (map (update-in g [:name] pr-str)
                                                  [:seq_region_name :start :end :name]))
                           "\n")))))
    (utils/bgzip-index out-file)))

(defn- usage [options-summary]
  (->> ["Create file of priority regions with genes from current CIViC database (https://civicdb.org)"
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
