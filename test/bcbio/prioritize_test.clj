(ns bcbio.prioritize-test
  (:require [clojure.test :refer :all]
            [clojure.java.io :as io]
            [bcbio.prioritize.create :as create]
            [bcbio.prioritize.known :as known]
            [bcbio.run.fsp :as fsp]
            [me.raynes.fs :as fs]))

(def dirs {:data "test/data" :work "test/work"})
(def data-files
  {:ref (str (io/file (:data dirs) "ref" "GRCh37-transcripts.bed"))
   :known (str (io/file (:data dirs) "annotation" "cosmic-v68-GRCh37.vcf.gz"))
   :call (str (io/file (:data dirs) "calls" "ex1.bed.gz"))})

(deftest known-test
  (fsp/remove-path (:work dirs))
  (let [bin-file (str (fs/file (:work dirs) "binned-known.bed.gz"))]
    (testing "Creation of binned items for prioritization based on known inputs."
      (is (= bin-file
             (create/from-known (:ref data-files) #{(:known data-files)} bin-file))))
    (testing "Prioritization of BED file calls given binned inputs"
      (let [priority-file (str (fs/file (:work dirs) "ex1-priority.bed.gz"))]
        (is (= priority-file
               (known/prioritize (:call data-files) bin-file priority-file)))))))
