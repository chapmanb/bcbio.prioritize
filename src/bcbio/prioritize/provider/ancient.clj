(ns bcbio.prioritize.provider.ancient
  "Genes found in pharmGKB based on selection from ancient DNA sequencing"
  (:import [htsjdk.samtools.util BlockCompressedOutputStream])
  (:require [bcbio.prioritize.provider.civic :as civic]
            [bcbio.prioritize.db :as db]
            [bcbio.run.itx :as itx]
            [clj-time.core]
            [clj-time.format]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [dk.ative.docjure.spreadsheet :as ss]))

(defn- excel->genes*
  "Raw dictionary of evidence keyed by genes based on excel input."
  [excel-file]
  (letfn [(collect-by-gene [out-coll cur]
            (reduce (fn [coll k]
                      (assoc-in coll [(:gene cur) k]
                                (conj (get-in coll [(:gene cur) k] #{}) (get cur k))))
                    out-coll
                    (remove #{:gene} (keys cur))))]
    (reduce collect-by-gene {}
            (concat
             (->> (ss/load-workbook-from-file excel-file)
                  (ss/select-sheet "Disease")
                  (ss/select-columns {:B :gene :O :disease})
                  (take 250))
             (->> (ss/load-workbook-from-file excel-file)
                  (ss/select-sheet "Chemical")
                  (ss/select-columns {:B :gene :P :drug}))))))

(defn excel->genes
  "Cleaned list of genes based on excel input."
  [excel-file]
  (->> (excel->genes* excel-file)
       (map (fn [[k v]] (assoc v :gene k)))
       (remove #(= (:gene %) "Symbol"))))

(defn excel->bed
  "Create a summarized BED file of genes of interest from excel summary"
  [in-file]
  (let [options {:build :GRCh38 :genome :human}
        out-file (format "ancient-%s-%s.bed.gz" (name (:build options))
                         (clj-time.format/unparse (clj-time.format/formatters :year-month-day)
                                                  (clj-time.core/now)))]
    (itx/with-tx-file [tx-out-file out-file]
      (with-open [wtr (io/writer (BlockCompressedOutputStream. (io/file tx-out-file)))]
        (doseq [g (->> (excel->genes in-file)
                       (map #(assoc (civic/gene-info (:genome options) (:build options) (:gene %)) :name %))
                       (sort-by (fn [x] [(Integer/parseInt (:seq_region_name x)) (:start x)]))
                       (map #(assoc % :seq_region_name (if (contains? #{:GRCh38 :hg19} (:build options))
                                                         (format "chr%s" (:seq_region_name %))
                                                         (:seq_region_name %)))))]
          (.write wtr (str (string/join "\t" (map (update-in g [:name] pr-str)
                                                  [:seq_region_name :start :end :name]))
                           "\n")))) )))

(defn- get-locations
  "Retrieve locations based on gene name for builds of interest"
  [options name]
  (map #(civic/gene-info (:genome options) % name) (:builds options)))

(defn- reorder-evidence
  "Reorder evidence and gene names for database storage"
  [orig]
  {:name (:gene orig)
   :evidence [{:origin "ancient" :variant "sweep"
               :diseases (:disease orig) :drugs (:drug orig)}]})

(defn excel->db
  "Load gene information from excel summary into queryable database."
  [in-file]
  (let [options {:builds [:GRCh37 :GRCh38] :genome :human}
        conn (db/init)]
    (db/load-data conn (map reorder-evidence (excel->genes in-file)) (partial get-locations options))))
