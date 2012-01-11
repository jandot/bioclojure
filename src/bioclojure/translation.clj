(ns bioclojure.translation
  (:require [clojure.string :as str])
  )

(defn six-frame
  [DNA-sequence]
  (def #^Integer size (count DNA-sequence))
  (println "sequence length:" (Integer/toString size))
  (println (str "Forward: " DNA-sequence))
  (let [revcomp (dna-complement (str/reverse DNA-sequence))]
    (println (str "Reverse-complement: " revcomp))
    (println "Reading frames:")
    (println (str "frame 1: " DNA-sequence))
    (println (str "frame 2: " (subs DNA-sequence 1)))
    (println (str "frame 3: " (subs DNA-sequence 2)))
    (println (str "frame 4: " revcomp))
    (println (str "frame 5: " (subs revcomp 1)))
    (println (str "frame 6: " (subs revcomp 2)))))

(defn dna-complement
  [DNA-sequence]
  (let [seq-dictionary
        {\A \T
         \T \A
         \G \C
         \C \G}]
    (apply str (replace seq-dictionary DNA-sequence))))
