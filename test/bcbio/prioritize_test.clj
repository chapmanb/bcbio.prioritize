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
   :cosmic (str (io/file (:data dirs) "annotation" "cosmic-v68-GRCh37.vcf.gz"))
   :oncomine (str (io/file (:data dirs) "annotation" "oncomine-GRCh37.vcf.gz"))
   :clinvar (str (io/file (:data dirs) "annotation" "clinvar-GRCh37.vcf.gz"))
   :simple-known (str (io/file (:data dirs) "annotation" "binned-simple.bed"))
   :simple-known-bed6 (str (io/file (:data dirs) "annotation" "binned-simple2.bed"))
   :call (str (io/file (:data dirs) "calls" "ex1.bed.gz"))
   :call-vcf (str (io/file (:data dirs) "calls" "ex2.vcf.gz"))
   })

(deftest known-test
  (fsp/remove-path (:work dirs))
  (let [bin-file (str (fs/file (:work dirs) "binned-cosmic.bed.gz"))
        onc-file (str (fs/file (:work dirs) "binned-oncomine.bed.gz"))
        clinvar-file (str (fs/file (:work dirs) "binned-clinvar.bed.gz"))]
    (testing "Cosmic: Creation of binned items for prioritization based on known inputs."
      (is (= bin-file
             (create/from-known (:ref data-files) #{(:cosmic data-files)} bin-file))))
    (testing "Prioritization of BED file calls given binned inputs"
      (let [priority-file (str (fs/file (:work dirs) "ex1-priority.bed.gz"))]
        (is (= priority-file
               (known/prioritize (:call data-files) bin-file priority-file)))))
    (testing "Prioritization of VCF file calls given binned inputs"
      (let [priority-file (str (fs/file (:work dirs) "ex2-priority-prebinned.vcf.gz"))]
        (is (= priority-file
               (known/prioritize (:call-vcf data-files) (:simple-known data-files) priority-file)))))
    (testing "Prioritize VCF file given a BED6 file"
      (let [priority-file (str (fs/file (:work dirs) "ex2-priority-prebinned-bed6.vcf.gz"))]
        (is (= priority-file
               (known/prioritize (:call-vcf data-files) (:simple-known-bed6 data-files) priority-file)))))
    (testing "Prioritize VCF file given a prepped known input file"
      (let [priority-file (str (fs/file (:work dirs) "ex2-priority-known.vcf.gz"))]
        (is (= priority-file
               (known/prioritize (:call-vcf data-files) bin-file priority-file)))))
    (testing "Oncomine: creation of binned items for prioritization."
      (is (= onc-file
             (create/from-known (:ref data-files) #{(:oncomine data-files)} onc-file))))
    (testing "ClinVar: creation of binned items for prioritization."
      (is (= clinvar-file
             (create/from-known (:ref data-files) #{(:clinvar data-files)} clinvar-file))))))
