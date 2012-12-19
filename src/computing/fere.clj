(ns computing.fere
  (:use [incanter.charts]
        [incanter.stats :only (euclidean-distance sample-uniform sample-normal)])
  (:require [incanter.core :as i]))

(defn expand [n coll]
  (flatten (map (partial repeat n) coll)))

(defn create-data [N T]
  (let [idx (range N)
        data-idx (expand T idx)
        x (sample-normal (* N T))
        ones (repeat (* N T) 1)
        id-rand (map * data-idx (sample-normal (* N T)))
        y (map + x id-rand)
        X (i/bind-columns ones data-idx x)]
    [y X]))
