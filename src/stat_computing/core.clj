(ns stat-computing.core
  (:use [cascalog.api]
        [clojure.math.numeric-tower :only (expt)]
        [clojure-csv.core])
  (:require [incanter.core :as i]
            [incanter.io :as io]
            [incanter.stats :as s]
            [incanter.charts :as c]))

;; Monte Carlo estimation of mean and standard error of the mean
;; page 154-155 of Statistical Computing for R

(defn demean [coll]
  (let [m (s/mean coll)]
    (map #(- % m) coll)))

(defn rnorm-diff []
  (i/abs (apply - (s/sample-normal 2))))

(defn mc-est [m]
  (let [rands (repeatedly rnorm-diff)]
    (s/mean (take m rands))))

(defn mc-stderr [m]
  (let [rands (repeatedly rnorm-diff)
        r-vec (take m rands)]
    (/ (i/sqrt (i/sum-of-squares (demean r-vec))) m)))

;; Feasible Generalized Least Squares (FGLS)

(defn linear-coefs
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

(def path "data/NLS80.csv")
(def my-data (io/read-dataset path :header true))
(def y (i/matrix (i/$ :lwage my-data)))
(def X (cofactor-mat my-data :exper :tenure :married :south :urban :black :educ))

(defn fgls
  "returns the coefficients for feasible generalized least squares."
  [y x & {:keys [intercept] :or {intercept true}}]
  (let [_x (if intercept (i/bind-columns (repeat (i/nrow x) 1) x) x)
        resid (linear-residuals y _x :intercept false)
        resid-beta (linear-coefs (i/sq resid) _x :intercept false)
        weight (i/sqrt (i/pow (i/mmult _x resid-beta) -1))]
    (linear-coefs (i/mult y weight) (mmult-col _x weight) :intercept false)))


;; Time-series; and introduction to time
(def ndvi
  [6217 8599 7074 8437 8471 8285 8342 9035 8356 7612 8538 6439 8232 7277 8588 7651 8824 4981 7251 7332 6179 4618 6506 8188 8320 6262 8094 8129 6773 7230 6417 6791 6285 6013 5786 8020 7588 6423 5734 6522 6481 7924 8067 7328 4249 8490 8591 4472 6335 8706 8076 8376 8861 8183 8712 6426 8314 8441 6643 7673 5193 8813 7902 7275 4480 7004 5691 5630 7540 8610 8981 5181 8947 8681 9072 8931 8879 8770 8702 6578 9027 8846 8530 6927 9128 5984 6133 4775 6707 6707 3392 7081 5806 6580 9108 5748 6784 8520 8597 9130 7585 6531 6768 7249 4992 4048 7988 8088 7418 4082 8056 2715 1899 8579 8852 8896 3010 8063 7985 8377 5503 8139 8672 8319 5995 8252 8835 8593 8909 6817 8488 7206 8561 8549 4261 5659 5924 8601 7302 2610 7610 7416 8978 8704 8528 8236 5400 6372 8387 9279 9175 8652 4637 4167 5624 5707 5404 4369 8607 2557 7840 9053 9502 8350 5512 8692 8274 4387 8192 8341 8042 6401 6284 7568 8354 7423 7064 4733 8441 5717 6456 4626 8160 7142 7135 5727 6847 8186 8179 8377 6998 6936 6722 5768 8552 6355 7360 7270 6069 3242 4972 5147 3720 6407 3887 6666 4915 5367 7383 5451 7442 7432 7807 7115 8622 7946 8488 7488 5482 4718 8206 8280 8822 8530 7810 8141 3207 5628 7737 7662 8606 8226 6252 8267 8668 8808 4407 8330 7473 8432 9038 7924 8581 9452 8347 5745 8741 5246 7643 6559 7669 3969 7377 4248 3973 8415 8031 8867 8481 4967 8804 8598 7935 8080 5842 8896 8789 7931 5019 8511 5378])

;; Estimate \rho

(def X
  (let [ones (repeat (count ndvi) 1)]
    (i/bind-columns ones (s/sample-uniform (count ndvi)))))

(def resid (:residuals (s/linear-model ndvi X :intercept false)))

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

;; (i/view (c/xy-plot (range 1 36) (map (partial rho-lag resid) (range
;; 1 36))))

(defn box-pierce-stat
  [resid-series max-lag]
  (let [T (count resid-series)]
    (/ (reduce + (map (partial rho-lag resid-series) (range 1 max-lag)))
       T)))

