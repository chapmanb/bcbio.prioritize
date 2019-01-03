(ns bcbio.prioritize.provider.ancient
  "Genes found in pharmGKB based on selection from ancient DNA sequencing"
  (:import [htsjdk.samtools.util BlockCompressedOutputStream])
  (:require [bcbio.prioritize.provider.civic :as civic]
            [bcbio.run.itx :as itx]
            [clj-time.core]
            [clj-time.format]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [dk.ative.docjure.spreadsheet :as ss]))

(defn excel->genes
  "Create a summarized BED file of genes of interest from excel summary"
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
                       (map (fn [[k v]] (assoc v :gene k)))
                       (remove #(= (:gene %) "Symbol"))
                       (map #(assoc (civic/gene-info (:genome options) (:build options) (:gene %)) :name %))
                       (sort-by (fn [x] [(Integer/parseInt (:seq_region_name x)) (:start x)]))
                       (map #(assoc % :seq_region_name (if (contains? #{:GRCh38 :hg19} (:build options))
                                                         (format "chr%s" (:seq_region_name %))
                                                         (:seq_region_name %)))))]
          (.write wtr (str (string/join "\t" (map (update-in g [:name] pr-str)
                                                  [:seq_region_name :start :end :name]))
                           "\n")))) )))
