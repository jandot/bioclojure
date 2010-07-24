(ns bioclojure
  (:use [bioclojure.vcf])
)

(vcf2tsv "../data/sample.vcf" "../data/sample.tsv")
