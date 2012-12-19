(ns computing.matching
  "Functions for a particular example of matching estimation."
  (:use [clojure.contrib.math :only (abs)])
  (:require [incanter.core :as i]
            [incanter.stats :as s]
            [incanter.charts :as c]))

(defn weighted-sum
  [weight-vec element-vec]
  (reduce + (map * weight-vec element-vec)))

(defn treated-diff
  [yvec-untreated wi-vec yi]
  (- yi (weighted-sum wi-vec yvec-untreated)))

(defn normalize-row
  [row]
  (let [sum (reduce + row)]
    (map #(/ % sum) row)))

(defn matching-estimator
  "w-vec is a vector of vectors.  Each vector i is a series of weights
  for the control units.  Nt x Nc different weights (matrix, or vector
  of vectors)."
  [w-mat yvec-treated yvec-untreated]
  (let [normalized-w (map normalize-row w-mat)]
    (/ (reduce + (map (partial treated-diff yvec-untreated) normalized-w yvec-treated))
       (count yvec-treated))))

(defn indicator-fn
  [])
