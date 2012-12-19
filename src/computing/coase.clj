(ns computing.coase
  (:use [incanter.charts]
        [incanter.stats :only (euclidean-distance sample-uniform)])
  (:require [incanter.core :as i]))

(defn wtp-consumer [])

(defn marginal-benefits
  [p]
  (+ 5 (* -0.5 p)))

(defn marginal-damages
  [p]
  (+ 2 (* 0.1 p)))


