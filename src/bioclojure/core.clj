(ns bioclojure.core
  (:use [clojure.contrib.io :only (reader with-out-writer)])
  (:require [clojure.string :as str])
)

(defn is-comment?
  "Checks if argument is a comment (i.e. starts with a '#').
   Returns: boolean."
  [line]
  (= \# (first line)))

(defn is-file-header?
  "Checks if argument is part of meta-information (i.e. starts with '##').
   Returns: boolean."
  [line]
  (= [\# \#] (take 2 (str line))))

(defn header
  "Returns header for file (i.e. all lines at top that start with '#')
   Returns: sequence containing header lines"
  [filename]
  (take-while is-comment? (line-seq (reader filename))))

(defn data-lines
  "Returns data lines in file (i.e. all lines that do not start with '#')
  Returns: sequence containing data lines"
  [filename]
  (drop-while is-comment? (line-seq (reader filename))))

(defn parsed-data-lines
  "Extract data elements from VCF file.
  Returns: sequence of sequences
  Example: (parsed-data-line \"example.vcf\"
           ; => ((\"20\" \"14370\" \"rs6054257\" \"G\" \"A\" \"29\" \"0\"
           ;      \"NS=58;DP=258;AF=0.786;DB;H2\"
           ;      \"GT:GQ:DP:HQ\" \"0|0:48:1:51,51\" \"1|0:48:8:51,51\" \"1/1:43:5\")
           ;     (\"20\" \"13330\" \".\" \"T\" \"A\" \"3\"
           ;      \"q10\" \"NS=55;DP=202;AF=0.024\"
           ;      \"GT:GQ:DP:HQ\" \"0|0:49:3:58,50\" \"0|1:3:5:65,3\" \"0/0:41:3\")
           ;     (\"20\" \"1110696\" \"rs6040355\" \"A\" \"G,T\" \"67\" \"0\"
           ;      \"NS=55;DP=276;AF=0.421,0.579;AA=T;DB\" \"GT:GQ:DP:HQ\" \"1|2:21:6:23,27\" \"2|1:2:0:18,2\" \"2/2:35:4\")"
  [filename]
  (map #(str/split % #"\t") (data-lines filename)))

(defn meta-information
  "Returns header for file (i.e. all lines at top that start with '##')
  Returns: sequence containing meta-information of VCF file"
  [filename]
  (take-while is-file-header? (header filename)))

(defn column-header
  "Returns column header line of file
  Returns: string containing column header line of VCF file"
  [filename]
  (str/replace-first (first (take 1 (drop-while is-file-header? (header filename)))) "#" ""))

(defn column-names
  "Get column names in VCF file
  Returns: sequence containing VCF column names"
  [filename]
  (str/split (str/replace (column-header filename) ":" "_") #"\t"))

(defn make-tag-value
  "Adds a '=1' to a tag that does not have a value"
  [s]
  (str/replace s #"^([^=]+)$" #(str (% 1) "=true")))

(defn create-map-for-info
  "Takes a string representing the INFO column for one variation and returns a 
  map of tag-value pairs. Tags that do not have a value (e.g. H2) are assigned
  a value of '1'. The values are all stored as strings, whether or not they 
  actually should represent integers or floats.
  Returns: map
  Example: (create-map-for-info \"NS=58;DP=258;AF=0.786;DB;H2\")
           ; => {\"NS\" \"58\", \"DP\" \"258\", \"AF\" \"0.786\", \"DB\" \"1\", \"H2\" \"1\"}"
  [line]
  (let [fields (str/split line #";")]
    (apply hash-map (clojure.core/flatten (map #(str/split % #"=") (map #(make-tag-value %) fields))))))

(defn extract-data
  [filename]
  (map #(apply hash-map (interleave (column-names filename) %)) (parsed-data-lines filename)))

(defn all-info-data
  "Extract all data from the INFO column.
  Returns: sequence of strings"
  [data]
  (map #(get % "INFO") data))

(defn all-info-tags
  "Extracts the unique tags that are present in the INFO column. Note: only
  tags that have a value (e.g. *not* H3).
  Returns: sorted sequence"
  [data]
  (sort (set (flatten (map #(keys %) (map #(create-map-for-info %) (all-info-data data)))))))

(defn all-format-data
  "Extract all data from the FORMAT column.
  Returns: sequence of strings"
  [data]
  (map  #(get % "FORMAT") data))

(defn all-format-tags
  "Extracts the unique tags that are present in the FORMAT column
  Returns: sorted sequence"
  [data]
  (sort (set (flatten (map #(str/split % #":") (all-format-data data))))))

(defn info-header
  "Create the header for the INFO columns: a sorted list of 'INFO-' concatenated to the tag.
  Returns: list"
  [data]
  (map #(str "INFO-" %) (all-info-tags data)))

(defn extract-info-value
  "Extracts the value for a given tag from an INFO string. Returns an empty string if
  that tag is not present in the INFO field.
  Example: (extract-info-value \"DP=17;CQ=INTRONIC;AB=0.75\" \"DP\") ; => \"17\"
  Returns: string"
  [string tag]
  (get (create-map-for-info string) tag "empty"))

(defn sample-names
  "Return sorted list of sample names"
  [filename]
  (sort (drop 9 (column-names filename))))

(defn sample-header
  "Create the header for the sample columns: a sorted list of sample-names concatenated
  to the FORMAT strings: \"NA00001-DP NA00001-GT NA00002-DP NA00002-GT\"
  Returns: list"
  [filename data]
  (for [s (sample-names filename) t (all-format-tags data)] (str s "-" t)))

(defn get-line-part-info
  "Create the part of the output line that concerns the INFO field"
  [m ait]
  (map #(extract-info-value (get m "INFO") %) ait))

(defn get-line-part-sample
  "Create the part of the output line that concerns a single sample"
  [sample m aft]
  (let [sample-data (str/split (get m sample) #":")
        sample-tags (str/split (get m "FORMAT") #":")
        sample-map (apply hash-map (interleave sample-tags sample-data))]
    (map #(get sample-map % "empty") aft)))

(defn get-line-part-all-samples
  "Create the part of the output line that concerns all samples"
  [m aft sn]
  (flatten (map #(get-line-part-sample % m aft) sn)))

(defn get-line
  "Create a complete data output line"
  [m cn ait aft sn]
  (let [common-fields (map #(get m %) cn)
        info-fields (get-line-part-info m ait)
        sample-fields (get-line-part-all-samples m aft sn)]
    (flatten (conj sample-fields info-fields common-fields))))

(defn get-all-lines
  "Create all data output lines"
  [filename data]
  (let [cn (take 7 (column-names filename))
        ait (all-info-tags data)
        aft (all-format-tags data)
        sn (sample-names filename)]
    (map #(get-line % cn ait aft sn) data)))

(defn vcf2tsv
  "Convert a VCF file to real tab-delimited format"
  [input-file output-file]
  (let [data (extract-data input-file)
        common-fields (take 7 (column-names input-file))]
    (with-out-writer output-file
      (println (str/join "\n" (meta-information input-file)))
      (println (str/join "\t" (flatten (conj (sample-header input-file data) (all-info-tags data) common-fields))))
      (println (str/join "\n" (map #(str/join "\t" %) (get-all-lines input-file data)))))))

;;;;;;;;;;;;;;;;;;

;(vcf2tsv "../../data/sample.vcf" "../../data/sample.tsv")
