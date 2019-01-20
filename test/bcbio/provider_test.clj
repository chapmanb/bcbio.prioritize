(ns bcbio.provider-test
  (:require [bcbio.prioritize.provider.civic :as civic]
            [bcbio.prioritize.provider.ancient :as ancient]
            [bcbio.prioritize.db :as pdb]
            [bcbio.prioritize.create :as create]
            [clojure.java.io :as io]
            [clojure.test :refer :all]
            [datomic.api :as d]
            [mount.core :as mount]
            [me.raynes.fs :as fs]))

(def dirs {:data "test/data" :work "test/work"})
(def data-files
  {:ref (str (io/file (:data dirs) "ref" "GRCh37-transcripts.bed"))
   :ancient (str (io/file (:data dirs) "prioritize" "ancient_prioritize.xlsx"))})

(deftest ^:db db-test
  (testing "Loading and retrieval of ancient prioritization into datomic"
    (mount/start)
    (ancient/excel->db (:ancient data-files))
    (let [db (d/db pdb/conn)
          g (pdb/get-gene db "CDK4")
          es (pdb/get-evidence db "sweep")]
      (println g)
      (println (first es))
      (is (= "CDK4" (:gene/name g)))
      (is (= "ancient" (:evidence/origin (first es)))))
    (mount/stop)))

;; (deftest civic-test
;;   (testing "Retrieval of genes from CIViC"
;;     (let [g (civic/gene "20")]
;;       (clojure.pprint/pprint (civic/gene->bed g))
;;       (is (= "ERBB2" (:name g)))
;;       (is (= ["HER2 Activating"] (vec (map :name (:variant_groups g))))))))

;; (deftest intogen-test
;;   (testing "Preparation of bin files from IntoGen"
;;     (let [dirname (str (fs/expand-home (io/file "~" "tmp" "intogen" "intogen_cancer_drivers-2014.12")))
;;           out-file (str (fs/file (:work dirs) "binned-intogen.bed.gz"))]
;;       (is (= out-file
;;              (create/from-known (:ref data-files) [dirname] out-file))))))
