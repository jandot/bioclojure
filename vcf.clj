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

;;;;;;;;;;;;;;;;;;

;(use '[clojure.set])
;(def w (sel (read-vcf "w.vcf") :cols [:chrom :pos :qual]))
;(def w_ids (map #(str (:chrom %) "_" (:pos %)) (:rows w)))
;(def b (sel (read-vcf "b.vcf") :cols [:chrom :pos :qual]))
;(def b_ids (map #(str (:chrom %) "_" (:pos %)) (:rows b)))
;(def c (conj-rows (:rows b) (take 2 (:rows w))))
;(def c_ids (map #(str (:chrom %) "_" (:pos %)) (:rows c)))
;(def i (intersection (set w_ids) (set c_ids)))
;(count i)
;(def d (col-names (conj-cols c c_ids) [:qual1 :pos1 :chrom1 :id]))
;(def x (col-names (conj-cols w w_ids) [:chrom2 :pos2 :qual2 :id]))
;(def y (dataset [:id :chrom1 :pos1 :qual1 :chrom2 :pos2 :qual2] (data-merge :id (:rows d) (:rows x))))
;(defn data-merge
;  "Merge two data sets on the given key"
;  [merge-key a b]
;  (let [indexed-b (zipmap (map merge-key b) b)]
;    (map #(into % (indexed-b (merge-key %))) (filter #(contains? indexed-b (merge-key %)) a))))
;
;(def v (sel (read-vcf "v.vcf") :cols [:chrom :pos :qual]))
;(def v_ids (map #(str (:chrom %) "_" (:pos %)) (:rows v)))
;(def v (col-names (conj-cols v_ids v) [:chrom-pos :chrom :pos :qual]))
;(def golden (col-names (read-dataset "golden_standard_na12878_filtered_CTR_and_exons.tsv" :delim \tab) [:id-golden :is-present]))
;(def v-golden ($join [:id-golden :chrom-pos] golden v))

;(def a (sel (read-vcf "a.vcf") :cols [:chrom :pos :qual]))
;(def a_ids (map #(str (:chrom %) "_" (:pos %)) (:rows a)))
;(count (intersection (set a_ids) (set v_ids)))
;(def a (col-names (conj-cols a a_ids) [:chrom1 :pos1 :qual1 :id]))
;(def v (col-names (conj-cols v v_ids) [:chrom2 :pos2 :qual2 :id]))
;
;(def i (dataset [:id :chrom1 :pos1 :qual1 :chrom2 :pos2 :qual2] (data-merge :id (:rows a) (:rows v))))
;(save i "intersection.tsv")
;(view (scatter-plot :qual1 :qual2 :data i))
