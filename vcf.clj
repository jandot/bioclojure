(use '(incanter core io charts stats))
(use '[clojure.contrib.duck-streams :only (read-lines reader)])
(use '[clojure.contrib.str-utils])
(use '[clojure.contrib.str-utils2 :only (lower-case split)])

(defn is-comment?
	"Checks if argument is a comment (i.e. starts with a '#')"
	[line]
	(= \# (first line)))

(defn is-file-header?
	"Checks if argument is part of meta-information (i.e. starts with '##')"
	[line]
	(= [\# \#] (take 2 (str line))))

(defn header
	"Returns header for file (i.e. all lines at top that start with '#')"
	[filename]
	(take-while is-comment? (line-seq (reader filename))))

(defn data-lines
	"Returns data lines in file (i.e. all lines that do not start with '#')"
	[filename]
	(drop-while is-comment? (line-seq (reader filename))))

(defn replace-item
	"Returns a list with the n-th item of l replaced by v."
	[l n v]
	(concat (take n l) (list v) (drop (inc n) l)))

(defn qual-as-float
	"Casts quality score column (col 6) from string to float"
	[split-line]
	(replace-item split-line 5 (Float. (nth split-line 5))))

(defn parsed-data-lines
	""
	[filename]
	(map #(re-split #"\t" %) (data-lines filename)))

(defn meta-information
	"Returns header for file (i.e. all lines at top that start with '##')"
	[filename]
	(take-while is-file-header? (line-seq (reader filename))))

(defn column-header
	"Returns column header line of file"
	[filename]
	(re-sub #"#" "" (first (take 1 (drop-while is-file-header? (line-seq (reader filename)))))))

(defn column-names
	[filename]
	(lazy-seq (map #(keyword %) (map #(lower-case %) (re-split #"\t" (re-gsub #":" "_" (column-header filename)))))))

(defn read-vcf
	"Read VCF file into incanter Dataset"
	[filename]
	(dataset (column-names filename) (map #(qual-as-float %) (parsed-data-lines filename))))

(defn create-map-for-info
	[line]
	(let [fields (split line #";")]
          (apply hash-map (flatten (filter #(= (count %) 2) (map #(split % #"=") fields))))))

;;;;;;;;;;;;;;;;;;

(def a (read-vcf "./data/sample.vcf"))

(def all-info (map :info (:rows a)))

(map #(get (create-map-for-info %) "DP") all-info)

