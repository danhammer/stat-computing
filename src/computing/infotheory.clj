(ns computing.infotheory
  (:use [clojure.math.numeric-tower :only (sqrt floor abs expt)])
  (:require [incanter.core :as i]))

(defn- ordinal-idx
  "Returns a sequence of indices that rank the values of the supplied
  time series in ascending order.  If there are equal values, the
  lexicographic ordering kicks in and the order at which the values
  appear is used to order the indices.

  Example:
    (take 5 ndvi) => (6217 8599 7074 8437 8471)
    (ordinal-idx (take 5 ndvi)) => (0 2 3 4 1)"
  [sub-ts]
  (let [indexed-vector (map-indexed vector sub-ts)]
    (map first
         (sort-by second indexed-vector))))

(defn- permutation-count
  "Returns a map of the ordinal sequences of length `D` and their
  count; note that the offset is fixed at 1"
  [D ts]
  (let [subs (partition D 1 ts)]
    (frequencies (map ordinal-idx subs))))

(defn- log-fn [x]
  (* x (i/log2 x)))

(defn- to-freq
  "Returns the normalized frequency of the supplied column"
  [coll]
  (let [total (reduce + coll)]
    (map #(/ % total) coll)))

(defn permutation-entropy
  "Normalixed permutation entropy based on the Shannon entropy
  distribution"
  [D ts]
  (let [pi-seq (to-freq (vals (permutation-count D ts)))
        scale-by (* -1 (/ 1 (i/log2 (i/factorial D))))]
    (* scale-by
       (reduce + (map log-fn pi-seq)))))

(defn kl-entropy
  "Returns the normalized Kullback-Leibler entropy (KLE) information
  measure, which quantifies the distance between the orginal pattern
  probability distribution of `ts` and the uniform distribution"
  [D ts]
  (- 1 (permutation-entropy D ts)))

(defn critically-bad?
  [bad all-types full]
  (let [total-nil (count (positions nil? full))]
    (if (or (nil? bad) (neg? total-nil))
      (false? all-types))))


(defn break-series
  [len & {:keys [avg1 avg2] :or {avg1 1 avg2 0.5}}]
  {:pre [(even? len)]}
  (let [half-len (/ len 2)]
    (concat (s/sample-normal half-len :mean avg1 :sd 0.05)
            (s/sample-normal half-len :mean avg2 :sd 0.05))))

(defn std-series
  [len & {:keys [avg] :or {avg 1}}]
  (s/sample-normal len :mean avg :sd 0.05))

(defn boot-hist
  [D B len f]
  (let [series (map f (repeat B len))]
    (apply merge-with +
           (map (partial permutation-count D) series))))

(defn plot-hist [f]
  (i/view (c/histogram (boot-hist 5 1000 4000 f) :nbins 120 :density true)))

(defn chi-sq [dist1 dist2]
  (let [norm-diff (fn [x y] (double (/ (i/sq (- x y)) y)))]
    (reduce + (map norm-diff dist1 dist2))))
