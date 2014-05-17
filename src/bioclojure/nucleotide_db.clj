(comment "
  FASTA data is held in a map structure:
    {:description \"first line in the FASTA format\"
     :sequence \"AACCGGTT\"}
")

(ns bioclojure.nucleotide-db
  (:use [clojure.contrib.io :only (reader)]
        [clojure.contrib.http.agent :only (string http-agent) :as http :exclude (bytes)]
        [clojure.contrib.str-utils :only (re-split) :as str-utils])
  (:require [clojure.string :as str])
)

(defn entrez-nucleotide
  "Take a string accession ID and download a fasta file from NCBI and parse fasta format. The function returns a map with two keys, :description and :sequence."
  [accessionID]

  ;request nucleotide record based on accession ID from caller
(def string-request (string (http/http-agent (str "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=" accessionID "&rettype=fasta&retmode=text"))))

  ;split record into meaningful data
  (def request-description (first (str-utils/re-split #"\n" string-request 2)))
  (def request-sequence (str/replace  (first (rest (str-utils/re-split #"\n" string-request 2)))
                                     #"\n" ""))

  ;return map
  {:description request-description
    :sequence request-sequence})
