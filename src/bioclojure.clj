;; To generate the uberjar:
;; - lein deps
;; - lein compile
;; - lein uberjar
(ns bioclojure
  (:use [incanter core io stats charts]
        [bioclojure.vcf]
        [bioclojure.fastq]
        [bioclojure.nucleotide-db]
        [bioclojure.translation])
  (:gen-class)
)

(def help-string
  "The bioclojure command takes a list of arguments. The first one is always the
  command that should be run, for example: 'java -jar bioclojure.jar say-hi'. There
  might be additional arguments that will act as the input for the command that
  is selected. Try:
  * java -jar bioclojure.jar say-hi
  * java -jar bioclojure.jar help
  * java -jar bioclojure.jar vcf2tsv ./data/sample.vcf ./data/sample.tsv")

(defn -main
  [command & args]
  (case command
        "help" (println help-string)
        "vcf2tsv" (vcf2tsv (first args))
        "say-hi" (println "Hi there! This actually works...")
        (println "Unknown command")))
