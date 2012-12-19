(ns computing.dice
  (:require [incanter.core :as i]
            [incanter.optimize :as optimize]
            [incanter.charts :as c]
            [incanter.stats :as s]))

;; Given the observed average from rolling a six-sided die many times,
;; what are the probabilities associated with each side of the die?
;; For a fair die, all sides are equally likely, such that the
;; probability is 1/6 for each side.  If we roll the fair die
;; thousands of times, the distribution of observed averages will be
;; centered around 3.5.  Suppose, however, we observe an average of
;; 4.5.  There are an infinite number of probability sequences that
;; would yield a distribution centered around 4.5.  We find the
;; distribution that is centered around the empirical average AND that
;; maximizes the entropy around the empirical average, subject to the
;; constraint that the probabilities sum to one.  This is a
;; prototypical example in information theory.

;; There are two functions that are surfaced:
;; Example:
;;   (prob-seq 4.5) => (0.054 0.079 0.114 0.165 0.240 0.347)
;; Example:
;;   (graph-probs 4.5) => ... bar chart of above probabilities ...

(defn- prob-die [beta i]
  (let [exp-fn (fn [b x] (i/exp (* -1 x b)))]
    (i/div (exp-fn beta i)
           (reduce + (map (partial exp-fn beta) (range 1 7))))))

(defn- tot-sum
  [beta]
  (let [iseq (range 1 7)
        prob-seq (map (partial prob-die beta) iseq)]
    (reduce + (map * iseq prob-seq))))

(defn- beta-param
  "Accepts the empirical average and finds the value where the
  function is equal to zero using the Newton-Raphson method."
  [empirical-avg] {:pre [(> empirical-avg 1) (< empirical-avg 6)]}
  (let [diff-fn (fn [t] (- (tot-sum t) empirical-avg))]
    (loop [x0 0 iter 100 diff 100]
      (if (or (zero? iter) (< diff 1e-10))
        x0
        (let [x1  (- x0 (/ (diff-fn x0) ((optimize/derivative diff-fn) x0)))]
          (recur x1 (dec iter) (i/abs (- x1 x0))))))))

(defn prob-seq
  "Accepts the average of the simulated die rolls, and returns the
  sequence of probabilities that maximizes the entropy subject to (1)
  the average is equal to the empirical average, and (2) the
  probabilities sum to one."
  [empirical-avg]
  (let [beta (beta-param empirical-avg)]
    (map (partial prob-die beta) (range 1 7))))

(defn graph-probs
  "Plots the probabilities associated with each face of the dice,
  given the observed average from simulation."
  [empirical-avg]
  (doto (c/bar-chart (range 1 7) (prob-seq empirical-avg))
    (c/set-y-label "probabilities")
    (c/set-x-label "die faces")
    (c/set-title (str "observed average: " empirical-avg))))

(defn dyn-graph
  "plot the probabilitiy distribution with a slider to adjust the
  empirical average"
  []
  (let [x (range 1 7)]
    (def pdf-chart (doto (c/scatter-plot)
                     (c/set-y-label "probabilities")
                     (c/set-x-label "die faces")))
    (i/view pdf-chart) 
    (c/sliders [mean (range 1.01 5.99 0.01)]
               (i/set-data pdf-chart [x  (prob-seq mean)]))))
