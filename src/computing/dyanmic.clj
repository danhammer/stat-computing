(ns computing.dynamic
  (:use [incanter.charts]
        [incanter.stats :only (euclidean-distance)])
  (:require [incanter.core :as i]))

(def iteration-spec
  (let [Rgrid (range 0.1 100 0.1)]
    {:thresh 0.0001
     :grid Rgrid
     :init-vals (repeat (count Rgrid) 0)}))

(def util-spec
  "A map that defines the utility and resource transition functions."
  {:eta 2
   :eps 0.01
   :alpha 1
   :beta 0.9
   :gamma 0
   :k 0
   :delta 0})

(defn consumption-vector
  "Returns a vector of /positive/ consumption possibilities for a
  given resource endowment `r`"
  [{:keys [eta eps alpha beta k gamma delta]} r] 
 (let [val (+ (i/pow r gamma) (* r (- 1 delta)))]
    (filter pos? (map #(- val %) (:grid iteration-spec)))))

(defn util
  "Accepts the parameter map `params` and returns the value of the
  utility function at the value of consumption `c`"
  [{:keys [eta eps alpha beta k gamma delta]} c]
  (let [eta-sub (- 1 eta)]
    (* (/ alpha eta-sub)
       (- (i/pow (+ c eps) eta-sub) k))))

(defn best-value
  "Accepts the parameter map that defines the utility function and
  discount rate; returns the value of the best possible consumption
  allocation.  Relies on the fact that the value function is
  monotonically increasing in resource consumpion (nonsatiation).
  That is, we can collect the maximum value function instead of
  keeping track of the index of the maximum consumption, and then
  calculating the value."
  [util-params old-vals r]
  (let [util-vec (map (partial util util-params)
                      (consumption-vector util-params r))
        bellman-fn (fn [x y] (+ x (* (:beta util-params) y)))]
    (reduce max (map bellman-fn util-vec old-vals))))

(defn iter-val
  "Accepts a parameterization of the utility function and the
  iteration process, and returns an incanter plot instance of the
  convergence process to the optimal resource consumption path.

    Example usage: (iter-val iteration-spec util-spec)"
  [iteration-params util-params]
  (let [graph (xy-plot)
        {:keys [thresh grid init-vals]} iteration-params]
    (loop [old-val init-vals
           diff 100]
      (if (< thresh diff)
        (let [new-val (map (partial best-value util-params old-val) grid)]
          (do (add-lines graph grid old-val)
              (recur new-val (euclidean-distance new-val old-val))))
        (-> graph
            (set-title "Convergence of value functions")
            (set-x-label "Resource")
            (set-y-label "Value"))))))
