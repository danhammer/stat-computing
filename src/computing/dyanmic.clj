(ns computing.dynamic
  (:use [incanter.charts]
        [incanter.stats :only (euclidean-distance)])
  (:require [incanter.core :as i]))

(def iteration-spec
  (let [Rgrid (range 0.1 8 0.1)]
    {:thresh 0.0001
     :grid Rgrid
     :init-vals (repeat (count Rgrid) 0)}))

(def util-spec
  "A map that defines the utility and resource transition functions."
  {:nu 2
   :eps 0.1
   :alpha 1
   :beta 0.95
   :cons 1
   :apprec 0.3
   :delta 0.1})

(defn consumption-vector
  "Returns a vector of /positive/ consumption possibilities for a
  given resource endowment `r`"
  [{:keys [nu eps alpha beta cons apprec delta]} r]
  (let [val (+ (i/pow r apprec) (* r (- 1 delta)))]
    (filter pos? (map #(- val %) Rgrid))))

(defn util
  "Accepts the parameter map `params` and returns the value of the
  utility function at the value of consumption `c`"
  [{:keys [nu eps alpha beta cons a]} c]
  (let [nu-sub (- 1 nu)]
    (* (/ alpha nu-sub)
       (- (i/pow (+ c eps) nu-sub) cons))))

(defn best-consumption
  "Accepts the "
  [param-map old-vals r]
  (let [util-vec (map (partial util param-map)
                      (consumption-vector param-map r))
        disc-val (map #(* (:beta param-map) %) old-vals)
        last-vec (map + util-vec disc-val)]
    (reduce max last-vec)))

(defn iter-val
  [iteration-params util-params n]
  (let [plot1 (xy-plot)
        {:keys [thresh grid init-vals]} iteration-params]
    (loop [old-val init-vals
           diff 100]
      (if (< thresh diff)
        (let [new-val (map (partial best-consumption util-params old-val) grid)]
          (do (add-lines plot1 grid old-val)
              (recur new-val (euclidean-distance new-val old-val))))
        plot1))))
