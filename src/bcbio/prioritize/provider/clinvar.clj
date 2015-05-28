(ns bcbio.prioritize.provider.clinvar
  "Prepare prioritization of bins based on evidence from ClinVar: http://www.clinvar.com/"
  )

(defn vc->rec
  "Convert variant context into a map of descriptive key pairs."
  [vc]
  (let [clnsig {"0" "Uncertain significance"
                "1" "not provided"
                "2" "Benign"
                "3" "Likely benign"
                "4" "Likely pathogenic"
                "5" "Pathogenic"
                "6" "drug response"
                "7" "histocompatibility"
                "255" "other"}
        clnorigin {"0" "unknown"
                   "1" "germline"
                   "2" "somatic"
                   "4" "inherited"
                   "8" "paternal"
                   "16" "maternal"
                   "32" "de-novo"
                   "64" "biparental"
                   "128" "uniparental"
                   "256" "not-tested"
                   "512" "tested-inconclusive"
                   "1073741824" "other"}
        clnrevstat {"mult" "Classified by multiple submitters"
                    "single" "Classified by single submitter"
                    "not" "Not classified by submitter"
                    "exp" "Reviewed by expert panel"
                    "prof" "Reviewed by professional society"}]
    {:origin "clinvar"
     :id (.getID vc)
     :clinical-origin (get clnorigin (.getAttribute vc "CLNORIGIN" "0") "unknown")
     :significance (get clnsig (.getAttribute vc "CLNSIG" "1") "not provided")
     :review-status (get clnrevstat (.getAttribute vc "CLNREVSTAT" "not") "Not classified by submitter")
     :disease-name (.getAttribute vc "CLNDBN" "")
     :disease-db (.getAttribute vc "CLNDSDB" "")}))
