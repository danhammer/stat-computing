(ns computing.kernel
  "Functions for a particular example of kernel density smoothing
toward kernel regression analysis.  Not suitable for general kernel
smoothing."
  (:use [clojure.contrib.math :only (abs)]
        [forma.classify.logistic])
  (:require [incanter.core :as i]
            [incanter.stats :as s]
            [incanter.charts :as c]))

(defn tri-kern
  "Returns the value of the triangular kernel function with bandwidth
  `h` for a reference value `xs` that corresponds to a value at a
  regular interval on the x-axis and the objective point `xi` which
  corresponds to a point in the histogram that is the objective of the
  kernel density smoothing.  If 5 parameters are supplied, then
  returns the value of the bivariate kernel density."
  ([h xs xi]
     (let [abs-u (abs (/ (- xs xi) h))
           m (if (< abs-u 1) 1 0)]
       (* m (- 1 abs-u))))
  ([h xs ys xi yi]
     (let [x-diff (abs (/ (- xs xi) h))
           y-diff (abs (/ (- ys yi) h))
           m (if (and (< x-diff 1) (< y-diff 1)) 1 0)]
       (* m (- 1 x-diff) (- 1 y-diff)))))

;; weights are summing to one

(defn unif-kern
  "Returns the value of the uniform kernel function with bandwidth `h`
  for a reference value `xs` that corresponds to a value at a regular
  interval on the x-axis and the objective point `xi` which
  corresponds to a point in the histogram that is the objective of the
  kernel density smoothing."
  ([h xs xi]
     (let [abs-u (abs (- xs xi))
           m (if (< abs-u h) 1 0)]
       (* m abs-u)))
  ([h xs xi ys yi]
     (let [abs-u (abs (- xs xi))
           m (if (< abs-u h) 1 0)]
       (* m abs-u))))

(defn kern-dens
  "Accepts a kernel smoothing function, a bandwidth length `h`, a
  vector of hits that generate the histogram `xi-vec` and a single
  point on the x-axis `xs`.  Returns the value of the kernel density
  smoothing at `xs`. ONLY one bandwidth for both x and y is allowed."
  ([f h xi-vec xs]
     (/ (reduce + (map (partial f h xs) xi-vec))
        (* h (count xi-vec))))
  ([f h yi-vec xi-vec xs]
     {:pre [(= (count yi-vec) (count xi-vec))]}
     (/ (reduce + (map * yi-vec (map (partial f h xs) xi-vec)))
        (* h (count xi-vec)))))

(defn bivariate-kexp
  "Conditional expectation of y, given xs.  Avoid divide by zero
  error."
  [f h yi-vec xi-vec xs]
  (let [kern-x (kern-dens f h xi-vec xs)]
    (if (zero? kern-x)
      nil
      (/ (kern-dens f h yi-vec xi-vec xs) kern-x))))

(defn bivariate-kreg
  "Returns a vector of conditional expectations of the dependent
  variable in `yi-vec`, conditional on each of the domain values
  `xs-vec`, knowing the independent variable in `xi-vec`

   Example usage:
     (def xs (range -3 3 0.1))
     (def ys (range -3 3 0.1))
     (def x (s/sample-normal 100))
     (def y (map + x (s/sample-normal 100 :sd 0.1)))
     (bivariate-kreg tri-kern 1 y x xs)"
  [f h yi-vec xi-vec xs-vec]
  (map (partial bivariate-kexp f h yi-vec xi-vec) xs-vec))

(defn plot-kreg
  "Note that there are possible nil values, which are interpolated by
  incanter (flat interpolation for line graphs).

    Example usage: (plot-kreg tri-kern 0.2)"
  [f h]
  (let [xs (range -3 3 0.1)
        ys (range -3 3 0.1)
        x (s/sample-normal 100)
        y (map + x (s/sample-normal 100 :sd 0.2))]
    (prn (bivariate-kreg f h y x xs))
    (-> (c/scatter-plot x y)
        (c/add-lines xs (bivariate-kreg f h y x xs))
        (i/view))))

(defn plot-kdens
  "Returns a view of the kernel density smoothing for the supplied
  kernel smoothing funciton and bandwidth for a sample of 1000 points
  from a standard normal distriubtion."
  [f h]
  (let [x (range -3 3 0.1)
        x-rand (s/sample-normal 1000)
        fx (map (partial kern-dens f h x-rand) x)]
    (i/view (c/xy-plot x fx))))

(defn min-idx
  "Returns the index of the minimum value of the supplied collection.
  If there is more than one minimum value, then the first index will
  be returned."
  [coll]
  (.indexOf coll (reduce min coll)))

(defn min-val
  "Returns the value in the domain `x-coll` that corresponds to the
  minimum value in the range `y-coll`"
  [x-coll y-coll]
  (nth x-coll (min-idx y-coll)))

(defn jack-kdens
  "Returns the scalar value of the jackknife-like procedure toward an
  estimate of the ISE(h) value at a given bandwidth h.  The `rand-vec`
  parameter corresponds to the vector of values that generate the
  histogram."
  [f h rand-vec]
  (let [fx (map (partial kern-dens f h rand-vec) rand-vec)]
    (* 2 (s/mean fx))))

(defn sq-kdens
  "Returns the scalar value of the sum of squared values of the kernel
  density values."
  [f h x-vec rand-vec]
  (let [fx (map (partial kern-dens f h rand-vec) x-vec)]
    (reduce + (map i/sq fx))))

(defn cv
  "Accepts a functional form `f`, a bandwidth `h`, and a standard
  deviation `sd` for the random sample from the normal 
 distribution (with std. deviation). Returns the
  Cross-Validation (CV) value of the estimator for the integrated mean
  squared error (IMSE)."
  [f sd h]
  (let [x (range -3 3 0.1)
        xr (s/sample-normal 1000 :sd sd)]
    (- (sq-kdens f h x xr)
       (jack-kdens f h xr))))

(defn optimal-bw
  "Returns the value of the optimal bandwidth, given a kernel
  functional and the standard deviation of the normal distribution
  from which the points are sampled."
  [f sd]
  (let [h-vec (range 0.3 0.45 0.01)
        v (map (partial cv f sd) h-vec)]
    (min-val h-vec v)))

(defn bw-vec
  "Used for simulation.  Accepts a kernel functional and the number of
  simulations to run `B` in order to bootstrap the mean and standard
  deviation of the optimal kernel bandwidth."
  [f B]
  (let [ones (repeat B 1)]
    (map (partial optimal-bw f) ones)))

(defn plot-cv [f]
  (let [h-vec (range 0.35 0.45 0.005)
        v (map (partial cv f) h-vec)]
    (i/view (c/xy-plot h-vec v))))

;; Propensity score example; selection on observables k = 3, three
;; observables, one treatment, and intercept (5 factors to generate y
;; variable)

(defn cofactor-mat
  [n k]
  (let [rand-vec (s/sample-normal (* n k))]
    (i/matrix rand-vec k)))

(defn make-binary [thresh x]
  (if (> x thresh) 1 0))

(defn treatment-vec
  [n X]
  (let [true-beta (i/matrix [1 2 3 2])
        e (s/sample-normal n :sd (i/sqrt 0.5))
        with-int (i/bind-columns X (i/matrix (repeat n 1)))]
    (i/matrix
     (map (partial make-binary 1.5)
          (i/to-vect (i/plus e (i/mmult with-int true-beta)))))))

(defn outcome-vec [X D]
  (let [ones (i/matrix (repeat (i/nrow X) 1))
        full-mat (i/bind-columns ones X D)
        true-beta (i/matrix [1 1 1 1 3])]
    (i/mmult full-mat true-beta)))

(defn propensity-scores
  [D X]
  (let [with-int (i/bind-columns X (i/matrix (repeat (i/nrow X) 1)))
        labels (to-double-rowmat (i/to-vect D))
        features (to-double-matrix (vec (i/to-vect with-int)))
        beta (logistic-beta-vector labels features 1e-8 1e-6 250)
        predict (fn [b x] (vec (.toArray (logistic-prob (to-double-rowmat b) (to-double-rowmat x)))))]
    (flatten (map (partial predict beta) (i/to-vect with-int)))))

(defn treat-this [y d p]
  (/ (* y d) p))

(defn treat-wt [d p]
  (/ d p))

(defn notreat-this [y d p]
  (/ (* y (- 1 d)) (- 1 p)))

(defn notreat-wt [d p]
  (/  (- 1 d) (- 1 p)))

(defn estimates [Y D P]
  (let [A (reduce + (map treat-this Y D P))
        B (reduce + (map treat-wt D P))
        C (reduce + (map notreat-this Y D P))
        D (reduce + (map notreat-wt D P))]
    (prn B)
    (prn D)
    (- (/ A B) (/ C D))))

(defn get-estimate [n]
  (let [X (cofactor-mat n 3)
        D (treatment-vec n X)
        Y (outcome-vec X D)
        P (propensity-scores D X)]
    (estimates Y D P)))
