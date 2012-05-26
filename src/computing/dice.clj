(ns computing.dice
  (:use [cascalog.api]
        [clojure.contrib.combinatorics :only [subsets]]
        [computing.data :only [ndvi]])
  (:require [incanter.core :as i]
            [incanter.io :as io]
            [incanter.stats :as s]
            [incanter.optimize :as optimize]
            [incanter.charts :as c]))

;; Maximize entropy subject to constraints, specifically that the
;; average equal the true average and the sum of the individual
;; probabilities is equal to 1

(defn prob-die [beta i]
  (let [exp-fn (fn [b x] (i/exp (* -1 x b)))]
    (i/div (exp-fn beta i)
           (reduce + (map (partial exp-fn beta) (range 1 7))))))

(defn tot-sum [beta]
  (let [iseq (range 1 7)
        prob-seq (map (partial prob-die beta) iseq)]
    (reduce + (map * iseq prob-seq))))

(defn true-beta [true-avg]
  (let [beta-range (range -2 2 0.001)
        avg-seq (map tot-sum beta-range)
        diff-seq (map i/abs (map #(- % true-avg) avg-seq))]
    (nth beta-range (.indexOf diff-seq (reduce min diff-seq)))))

(defn prob-seq [true-avg]
  (let [beta (true-beta true-avg)]
    (map (partial prob-die beta) (range 1 7))))

;; http://www.bestinclass.dk/index.clj/2010/01/hadoop-feeding-reddit-to-hadoop.html
