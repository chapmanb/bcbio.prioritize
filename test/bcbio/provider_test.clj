(ns bcbio.provider-test
  (:require [bcbio.prioritize.provider.civic :as civic]
            [bcbio.prioritize.create :as create]
            [clojure.java.io :as io]
            [clojure.test :refer :all]
            [me.raynes.fs :as fs]))

(def dirs {:data "test/data" :work "test/work"})
(def data-files
  {:ref (str (io/file (:data dirs) "ref" "GRCh37-transcripts.bed"))})

(deftest civic-test
  (testing "Retrieval of genes from CIViC"
    (let [g (civic/gene "20")]
      (clojure.pprint/pprint (civic/gene->bed g))
      (is (= "ERBB2" (:name g)))
      (is (= ["HER2 Activating"] (vec (map :name (:variant_groups g))))))))

;; (deftest intogen-test
;;   (testing "Preparation of bin files from IntoGen"
;;     (let [dirname (str (fs/expand-home (io/file "~" "tmp" "intogen" "intogen_cancer_drivers-2014.12")))
;;           out-file (str (fs/file (:work dirs) "binned-intogen.bed.gz"))]
;;       (is (= out-file
;;              (create/from-known (:ref data-files) [dirname] out-file))))))
