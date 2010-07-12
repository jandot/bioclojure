(use '[clojure.set])
(def w (sel (read-vcf "w.vcf") :cols [:chrom :pos :qual]))
(def w_ids (map #(str (:chrom %) "_" (:pos %)) (:rows w)))
(def b (sel (read-vcf "b.vcf") :cols [:chrom :pos :qual]))
(def b_ids (map #(str (:chrom %) "_" (:pos %)) (:rows b)))
(def c (conj-rows (:rows b) (take 2 (:rows w))))
(def c_ids (map #(str (:chrom %) "_" (:pos %)) (:rows c)))
(def i (intersection (set w_ids) (set c_ids)))
(count i)
(def d (col-names (conj-cols c c_ids) [:qual1 :pos1 :chrom1 :id]))
(def x (col-names (conj-cols w w_ids) [:chrom2 :pos2 :qual2 :id]))
(def y (dataset [:id :chrom1 :pos1 :qual1 :chrom2 :pos2 :qual2] (data-merge :id (:rows d) (:rows x))))
(defn data-merge
  "Merge two data sets on the given key"
  [merge-key a b]
  (let [indexed-b (zipmap (map merge-key b) b)]
    (map #(into % (indexed-b (merge-key %))) (filter #(contains? indexed-b (merge-key %)) a))))

(def v (sel (read-vcf "v.vcf") :cols [:chrom :pos :qual]))
(def v_ids (map #(str (:chrom %) "_" (:pos %)) (:rows v)))
(def v (col-names (conj-cols v_ids v) [:chrom-pos :chrom :pos :qual]))
(def golden (col-names (read-dataset "golden_standard_na12878_filtered_CTR_and_exons.tsv" :delim \tab) [:id-golden :is-present]))
(def v-golden ($join [:id-golden :chrom-pos] golden v))

(def a (sel (read-vcf "a.vcf") :cols [:chrom :pos :qual]))
(def a_ids (map #(str (:chrom %) "_" (:pos %)) (:rows a)))
(count (intersection (set a_ids) (set v_ids)))
(def a (col-names (conj-cols a a_ids) [:chrom1 :pos1 :qual1 :id]))
(def v (col-names (conj-cols v v_ids) [:chrom2 :pos2 :qual2 :id]))

(def i (dataset [:id :chrom1 :pos1 :qual1 :chrom2 :pos2 :qual2] (data-merge :id (:rows a) (:rows v))))
(save i "intersection.tsv")
(view (scatter-plot :qual1 :qual2 :data i))
