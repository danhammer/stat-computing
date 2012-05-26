(ns computing.econometrics
  (:use [cascalog.api]
        [clojure.math.numeric-tower :only (expt)]
        [clojure-csv.core]
        [computing.data :only [ndvi]])
  (:require [incanter.core :as i]
            [incanter.io :as io]
            [incanter.stats :as s]
            [incanter.charts :as c]))

;; Monte Carlo estimation of mean and standard error of the mean
;; page 154-155 of Statistical Computing for R

(defn demean
  "returns a vector with each element, less the mean of the vector"
  [coll]
  (let [m (s/mean coll)]
    (map #(- % m) coll)))

(defn rnorm-diff
  "returns a lazy sequence of the absolute difference between two
  standard normal distributions"
  []
  (i/abs (apply - (s/sample-normal 2))))

(defn mc-est
  "returns the monte carlo estimate of the difference based on `m`
  iterations between the random variables drawn from standard normal
  distributions"
  [m]
  (let [rands (repeatedly rnorm-diff)]
    (s/mean (take m rands))))

(defn mc-stderr
  "calculates the MC estimate of the std. error of the difference
  between two RVs distributed std. normal, based on `m` iterations"
  [m]
  (let [rands (repeatedly rnorm-diff)
        r-vec (take m rands)]
    (/ (i/sqrt (i/sum-of-squares (demean r-vec))) m)))

;; Feasible Generalized Least Squares (FGLS)

(defn linear-coefs
  "returns the vector of coefficients "
  [y x & {:keys [intercept] :or {intercept true}}]
  (let [_x (if intercept (i/bind-columns (repeat (i/nrow x) 1) x) x)
        xt (i/trans _x)
        xtx (i/mmult xt _x)]
    (i/mmult (i/solve xtx) xt y)))

(defn linear-fitted
  "fitted values from ols model"
  [y x & {:keys [intercept] :or {intercept true}}]
  (let [_x (if intercept (i/bind-columns (repeat (i/nrow x) 1) x) x)
        coefs (linear-coefs y _x :intercept false)]
    (i/mmult _x coefs)))

(defn linear-residuals
  "returns the residuals from a linear model; cribbed from
  incanter.stats linear model"
  [y x & {:keys [intercept] :or {intercept true}}]
  (let [_x (if intercept (i/bind-columns (repeat (i/nrow x) 1) x) x)
        coefs (linear-coefs y _x :intercept false)]
    (i/minus y (i/mmult _x coefs))))

(defn cofactor-mat
  "returns a cofactor matrix with intercept using built-in incanter
  load functions"
  [dataset & col-names]
  (let [my-func (fn [x] (i/$ x dataset))]
    (i/trans (i/matrix (map my-func col-names)))))

(defn head
  "returns the first n rows of a matrix, like the head function in R"
  [mat & n]
  (if (nil? n)
    (i/sel mat :rows (range 6))
    (i/sel mat :rows (apply range n))))

(defn mmult-col
  "element-wise multiplication of each column in a matrix `mat` by a
  column `v`"
  [mat v]
  (let [ks (range (i/ncol mat))
        mult-fn (fn [k] (i/mult (i/sel mat :cols k) v))]
    (apply i/bind-columns (map mult-fn ks))))

(defn fgls
  "returns the coefficients for feasible generalized least squares."
  [y x & {:keys [intercept] :or {intercept true}}]
  (let [_x (if intercept (i/bind-columns (repeat (i/nrow x) 1) x) x)
        resid (linear-residuals y _x :intercept false)
        resid-beta (linear-coefs (i/sq resid) _x :intercept false)
        weight (i/sqrt (i/pow (i/mmult _x resid-beta) -1))]
    (linear-coefs (i/mult y weight)
                  (mmult-col _x weight)
                  :intercept false)))

(def fgls-coeffs
  "returns the coefficients of FGLS estimation on the Card data set."
  (let [data (io/read-dataset "data/NLS80.csv" :header true)
        y (i/matrix (i/$ :lwage my-data))
        X (cofactor-mat my-data :exper :tenure :married :south :urban :black :educ)]
    (fgls y X)))

;; Newey-west

(defn S-naught
  [y x & {:keys [intercept] :or {intercept true}}]
  (let [_x (if intercept (i/bind-columns (repeat (i/nrow x) 1) x) x)
        resid (linear-residuals y _x :intercept false)
        xwt (mmult-col _x resid)]
    (i/div (i/mmult (i/trans xwt) xwt)
           (i/nrow xwt))))

;; Estimate \rho

(defn X2 [ts]
  (let [ones (repeat (count ts) 1)]
    (i/bind-columns ones (s/sample-uniform (count ts)))))

(def resid (:residuals (s/linear-model ndvi (X2 ndvi) :intercept false)))

(def rho-1
  (let [T (count resid)
        et (rest resid)
        et-1 (take (dec T) resid)]
    (:coefs (s/linear-model et et-1 :intercept false))))

(def rho-2
  (let [my-func (fn [[x y]] (* x y))]
    (/ (reduce + (map my-func (partition 2 1 resid)))
       (i/sum-of-squares (take (dec (count resid)) resid)))))

(defn rho-lag
  [coll lag]
  (let [my-func (fn [[& args]] (* (first args) (last args)))
        offset  (dec (- (count coll) lag))]
    (/ (reduce + (map my-func (partition (inc lag) 1 coll)))
       (i/sum-of-squares (take offset coll)))))

(defn box-pierce-stat
  [resid-series max-lag]
  (let [T (count resid-series)]
    (* (reduce + (map (comp i/sq (partial rho-lag resid-series)) (range 1 max-lag)))
       T)))


;; (i/view (c/xy-plot (range (count Yt)) Yt))
;; (i/view (c/xy-plot (range 1 100) (map (partial rho-lag resid) (range 1 100))))
;; (i/view (c/xy-plot (range 1 100) (map (partial box-pierce-stat
;; resid) (range 1 100))))


