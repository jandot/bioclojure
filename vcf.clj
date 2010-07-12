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

(defn qual-as-float
  "Casts quality score column (col 6) from string to float
  Returns: sequence"
  [split-line]
  (replace-item split-line 5 (Float. (nth split-line 5))))

(defn parsed-data-lines
  ""
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
  (dataset (column-names filename) (map #(qual-as-float %) (parsed-data-lines filename))))

;;;;;;;;;;;;;;;;;;

(def a (read-vcf "./data/sample.vcf"))

(def all-info (map :info (:rows a)))

(map #(get (create-map-for-info %) "DP") all-info)

