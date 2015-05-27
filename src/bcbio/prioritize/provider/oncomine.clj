(ns bcbio.prioritize.provider.oncomine
  "Parse VCF exported from oncomine: https://www.oncomine.org"
  (:import [htsjdk.samtools.util BlockCompressedOutputStream BlockCompressedInputStream])
  (:require [bcbio.prioritize.utils :as utils]
            [bcbio.run.fsp :as fsp]
            [bcbio.run.itx :as itx]
            [clojure.java.io :as io]
            [clojure.string :as string]))

(defn- passing-line?
  "Check if a line is a header or has correct REF and ALT alleles"
  [line]
  (or (.startsWith line "#")
      (let [parts (string/split line #"\t")
            ref (nth parts 3)
            alt (nth parts 4)]
        (and (not= "." ref)
             (not= "." alt)
             (not= ref alt)))))

(defn clean-file
  "Clean a VCF file, removing variants with missing or identical REF and ALTs."
  [in-file work-dir]
  (let [out-file (fsp/add-file-part in-file "clean" work-dir)]
    (when-not (itx/up-to-date? out-file in-file)
      (itx/with-tx-file [tx-out-file out-file]
        (with-open [rdr (io/reader (BlockCompressedInputStream. (io/file in-file)))
                    wtr (io/writer (BlockCompressedOutputStream. (io/file tx-out-file)))]
          (doseq [line (filter passing-line? (line-seq rdr))]
            (.write wtr (str line "\n"))))))
    (utils/bgzip-index out-file work-dir)))

(defn vc->rec
  "Convert variant context into a map of descriptive key pairs."
  [vc]
  {:origin "oncomine"
   :id (.getID vc)
   :mutation-class (.getAttribute vc "om_MutClass" "")
   :patient (.getAttribute vc "om_PATIENT")
   :cancer (.getAttribute vc "om_Cancer")})
