(use '(incanter core io charts stats))
(use '[clojure.contrib.duck-streams :only (read-lines reader)])
(use '[clojure.contrib.str-utils])
(use '[clojure.contrib.str-utils2 :only (lower-case split)])

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

(defn replace-item
  "Returns a list with the n-th item of l replaced by v.
  Example: (replace-item [1 2 3 4] 2 7) ; => [1 2 7 4]
  Returns: sequence"
  [l n v]
  (concat (take n l) (list v) (drop (inc n) l)))

(defn element-as-float
  "Casts nth element in sequence from string to float
  Returns: sequence"
  [lst n]
  (replace-item lst n (Float. (nth lst n))))

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
  (map #(re-split #"\t" %) (data-lines filename)))

(defn meta-information
  "Returns header for file (i.e. all lines at top that start with '##')
  Returns: sequence containing meta-information of VCF file"
  [filename]
  (take-while is-file-header? (line-seq (reader filename))))

(defn column-header
  "Returns column header line of file
  Returns: string containing column header line of VCF file"
  [filename]
  (re-sub #"#" "" (first (take 1 (drop-while is-file-header? (line-seq (reader filename)))))))

(defn column-names
  "Get column names in VCF file
  Returns: sequence containing VCF column names"
  [filename]
  (lazy-seq (map #(keyword %) (map #(lower-case %) (re-split #"\t" (re-gsub #":" "_" (column-header filename)))))))

(defn create-map-for-info
  "Takes a string representing the INFO column for one variation and returns a 
  map of tag-value pairs. In it current implementation any elements in the INFO
  column that are *not* tag-value (e.g. 'H3') are discarded. The values are also
  all stored as strings, whether or not they actually should represent integers
  or floats.
  Returns: map
  Example: (create-map-for-info \"NS=58;DP=258;AF=0.786;DB;H2\")
           ; => {\"NS\" \"58\", \"DP\" \"258\", \"AF\" \"0.786\"}"
  [line]
  (let [fields (split line #";")]
    (apply hash-map (flatten (filter #(= (count %) 2) (map #(split % #"=") fields))))))

(defn read-vcf
  "Read VCF file into incanter Dataset"
  [filename]
  (dataset (column-names filename) (map #(element-as-float % 5) (parsed-data-lines filename))))

;;;;;;;;;;;;;;;;;;

; For example: load a VCF file and print the sequence depth for all SNPs
(def a (read-vcf "./data/sample.vcf"))

(def all-info (map :info (:rows a)))

(println (map #(get (create-map-for-info %) "DP") all-info))

