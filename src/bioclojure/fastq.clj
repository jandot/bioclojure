(comment "
  FASTQ data is held in simple map structure:
    {:id \"abcd\"
     :seq \"AACCGGTT\"
     :qual \"27NN4Jk?\"}
")

(ns bioclojure.fastq
  (:use [clojure.contrib.io :only (reader)])
  (:require [clojure.string :as str])
)

(defn fastqfile2map
  "Reads a FASTQ stream and returns lazy list of maps"
  [filename]
  (let [s (partition 4 (line-seq (reader filename)))
        m1 (map #(zipmap [:id :seq :dummy :qual] %) s)]
      (map #(dissoc % :dummy) m1)))

(defn fastq2fasta
  "Given a map, returns the FASTA representation"
  [m]
  (str/join "\n" [(str/replace-first (:id m) "@" ">") (:seq m)]))

(defn fastqfile2fasta
  "Reads a FASTQ file and prints out FASTA format"
  [filename]
  (let [list-of-maps (fastqfile2map filename)
        list-of-fasta-strings (map fastq2fasta list-of-maps)]
      (str/join "\n" list-of-fasta-strings)))