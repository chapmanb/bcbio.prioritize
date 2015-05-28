(ns bcbio.prioritize.utils
  (:require [bcbio.run.itx :as itx]
            [bcbio.run.fsp :as fsp]
            [clojure.java.io :as io]
            [clojure.core.strint :refer [<<]]
            [me.raynes.fs :as fs]))

;; bgzip and tabix support

(defn- tabix-index
  "Tabix index input an input file inside a transactional directory."
  [bgzip-file ftype]
  (let [tabix-file (str bgzip-file ".tbi")]
    (when (or (itx/needs-run? tabix-file) (not (itx/up-to-date? tabix-file bgzip-file)))
      (itx/with-tx-file [tx-tabix-file tabix-file]
        (let [tx-bgzip-file (fsp/file-root tx-tabix-file)
              full-bgzip-file (str (fs/file bgzip-file))
              tmp-dir (str (fs/parent tx-bgzip-file))]
          (itx/check-run (<< "ln -s ~{full-bgzip-file} ~{tx-bgzip-file}"))
          (itx/check-run (<< "bcftools tabix -p ~{ftype} ~{tx-bgzip-file}")))))
    bgzip-file))

(defn- bgzip-file
  [in-file out-file]
  (when-not (itx/up-to-date? out-file in-file)
    (fsp/remove-path out-file))
  (itx/run-cmd out-file "bgzip -c ~{in-file} > ~{out-file}"))

(defmulti do-bgzip
  "Worker functionality for bgzip and indexing"
  (fn [in-file out-dir]
    (last (fsp/split-ext+ in-file))))

(defmethod do-bgzip".bed"
  [in-file out-dir]
  (let [out-file (str (fsp/file-root in-file out-dir) ".bed.gz")]
    [(bgzip-file in-file out-file) "bed"]))

(defmethod do-bgzip ".bed.gz"
  [in-file out-dir]
  (let [out-file (str (fsp/file-root in-file out-dir) ".bed.gz")]
    (when-not (or (= in-file out-file) (fs/exists? out-file))
      (fs/sym-link out-file in-file))
    [out-file "bed"]))

(defmethod do-bgzip ".vcf"
  [in-file out-dir]
  (let [out-file (str (fsp/file-root in-file out-dir) ".vcf.gz")]
    [(bgzip-file in-file out-file) "vcf"]))

(defmethod do-bgzip ".vcf.gz"
  [in-file out-dir]
  (let [out-file (str (fsp/file-root in-file out-dir) ".vcf.gz")]
    (when-not (or (= in-file out-file) (fs/exists? out-file))
      (fs/sym-link out-file in-file))
    [out-file "vcf"]))

(defn bgzip-index
  "Prepare an input file, handling multiple input types."
  ([in-file out-dir]
   (apply tabix-index (do-bgzip in-file out-dir)))
  ([in-file]
   (bgzip-index in-file (fs/parent in-file))))

;; BED support

(defn merge-bed
  [in-file work-dir]
  (let [out-file (str (fsp/file-root in-file work-dir) "-merged.bed.gz")
        cat-cmd (if (.endsWith in-file "bed.gz") "zcat" "cat")]
    (itx/run-cmd out-file
                 "~{cat-cmd} ~{in-file} | bedtools merge -c 4 -o distinct -i - |"
                 "bgzip -c > ~{out-file}")
    (bgzip-index out-file work-dir)))
