(ns bcbio.prioritize-test
  (:require [clojure.test :refer :all]
            [clojure.java.io :as io]
            [bcbio.prioritize.create :as create]
            [me.raynes.fs :as fs]))

(def dirs {:data "test/data" :work "test/work"})
(def data-files
  {:ref (str (io/file (:data dirs) "ref" "GRCh37-transcripts.bed"))
   :known (str (io/file (:data dirs) "annotation" "cosmic-v68-GRCh37.vcf.gz"))})

(deftest create-test
  (testing "Creation of binned items for prioritization based on known inputs."
    (let [out-file (str (fs/file (:work dirs) "binned-known.bed.gz"))]
      (is (= out-file
             (create/from-known (:ref data-files) #{(:known data-files)} out-file))))))
