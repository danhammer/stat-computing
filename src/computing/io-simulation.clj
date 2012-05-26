(ns computing.io-simulation
  (:use [cascalog.api]
        [clojure.math.numeric-tower :only (expt)])
  (:require [incanter.core :as i]
            [incanter.io :as io]
            [incanter.charts :as c]
            [incanter.stats :as s]))

(defn sample-bernoulli [n p]
  (map #(if (> % p) 1 0) (s/sample-uniform n)))

(defn expected-profit
  "returns expected profit, given high and low profit states along
  with the probability of the high state."
  [high low p]
  (+ (* p high) (* (- 1 p) low)))

(defn prob-seq
  "A sequence of length B of probabilities, each equivalent to finding
  the average of a sequence of bernoulli trials of length n with
  probability p of success."
  [B n p]
  (let [sd (i/sqrt (/ (* p (- 1 p)) n))]
    (s/sample-normal B :mean p :sd sd)))

(defn crit-discount-rate 
  "find the critical distcount rate.  Any rate above the critical
  value implies that the firm will not cheat; otherwise, the firm will
  cheat."
  [high-profit low-profit p]
  (/ high-profit
     (+ (* 2 (expected-profit high-profit low-profit p)) high-profit)))

(defn incentive-compat
  "incentive compatability (the firm will collude) when the return
  value is 0.  If the return value is 1, then the firm cheats."
  [d high-profit low-profit p]
  (if (> d (crit-discount-rate high-profit low-profit p)) 0 1))

(defn cheat-probability
  "monte-carlo integration of the area under the probability density
  function of observed probabilities that is to the left of the true
  discount value of the firm `d`.  The story is that a firm faces a
  `high-profit` and `low-profit` state with true probability `true-p`
  of the high profit state [and likewise (1 - `true-prob`) of the
  low-profit state]. The firm does not know `true-p` and can only
  infer its value from the history of states of length `n`.  By
  monte-carlo integration, we can find the number of times that the
  firm cheats -- indicated by a 1 as output from `incentive-compat`"
  [d high-profit low-profit B n true-p]
  (s/mean
   (map (partial incentive-compat d high-profit low-profit)
        (prob-seq B n true-p))))

;; Suppose my discount rate is 0.4.  At this discount rate, the
;; critical value of the probability of the high-profit state is 0.5.
;; Anything above 0.5, and the firm will not cheat.  The expected
;; profits in the next period are too high to risk cheating.  Suppose
;; that the TRUE probability of the high-profit state (p or alpha in
;; the problem set) is 0.58.  If this was known, then the firm will
;; NEVER cheat. The value will ALWAYS be 0:

;; (incentive-compat 0.4 0.58) => 0

;; However, suppose that the industries learn the probability of the
;; high-profit state, based on observed history.  The standard
;; deviation of a Bernoulli series is inversely related to the length
;; of that history -- the longer the history, the more tight the
;; estimate of the probability. Monte Carlo integration to find out
;; what the area uner the tail of the spread of the probabilities.

;; (cheat-probability 0.4 20 10 100000 100 0.58)

;; New industries will cheat more under these assumptions that older
;; industries.

(defn my-function1 [x]
  (s/pdf-normal x :mean 0.58 :sd (i/sqrt 0.02346)))

(defn my-function2 [x]
  (s/pdf-normal x :mean 0.58 :sd (i/sqrt 0.002346)))

(def myplot (c/function-plot my-function1 0.2 1))
(c/add-function myplot my-function2 0.2 1)
(c/add-lines myplot (repeat 100 0.5) (range 0 8.5))
(c/add-lines myplot (repeat 100 0.5) (range 0 8.5))
(c/add-lines myplot (repeat 100 0.5) (range 0 8.5))
;; (i/save myplot "~/Dropbox/berkeley/coursework/ARE202/mc-test/pdf.png")


;; NOW assume we know the states of the world; high and low, with
;; known probabilities.  The expected profit is therefore fixed.
;; Note that T = 1/lambda.  T is the time in between collusion events, and
;; lambda is the rate parameter in the exponential distribution.  The
;; maximum likelihood estimate is distributed N(lambda, lambda^2/n)

(defn incentive-compat2
  [d high-profit low-profit p l]
  (let [T (/ 1 l)
        split-profit (/ high-profit 2)
        ep (expected-profit high-profit low-profit p)
        A (/ (* d ep) (- 1 d))
        discount-cond (* A (- 1 (i/pow d T)))]
    (if (> discount-cond split-profit) 0 1)))

;; (incentive-compat2 0.4 20 10 0.51 0.1) => 0

(defn prob-seq-pois
  "A sequence of length B of probabilities, each equivalent to finding
  the average of a sequence of bernoulli trials of length n with
  probability p of success."
  [B n l]
  (let [sd (/ l (i/sqrt n))]
    (s/sample-normal B :mean l :sd sd)))

(defn cheat-probability2
  [d high-profit low-profit p B n true-l]
  (s/mean
   (map (partial incentive-compat2 d high-profit low-profit p)
        (prob-seq-pois B n true-l))))

(cheat-probability2 0.4 20 10 0.501 10000 2 0.08)
(cheat-probability2 0.4 20 10 0.501 10000 10 0.09)


;; Actually depends a lot on the spread of profits; difference high
;; and low.

(cheat-probability2 0.35 12 9 0.75 1000 100 0.2)

(defn inequality-cond
  [d high-profit low-profit p l]
  (let [T (/ 1 l)
        ep (expected-profit high-profit low-profit p)
        A (/ (* d ep) (- 1 d))]
    (* A (- 1 (i/pow d T)))))

(defn empirical-pdf
  [d high-profit low-profit p B n true-l]
  (let [cheat-cond (map (partial inequality-cond d high-profit low-profit p)
                        (prob-seq-pois B n true-l))
        rescale-cond (map #(- % 6) cheat-cond)]
    (doto (c/histogram rescale-cond :density true :nbins 30 :x-label "Cheat condition")
      (c/add-lines (repeat 100 0) (range 0 27))
      i/view)))

;; (empirical-pdf 0.35 12 9 0.75 100000 100 0.2)
;; (incentive-compat2 0.35 12 9 0.75 0.2) => 0 ;; should NEVER cheat,
;; but instead cheats 10% of the time when n=100:
;; (cheat-probability2 0.35 12 9 0.75 100000 100 0.2) => 0.10173
;; (cheat-probability2 0.35 12 9 0.75 100000 10 0.2) => 0.34513

(defn random-mat
  "creates a random matrix from the supplied probability distribution
  function of dimension `row` by `col`, with options.

  Example:
    (random-mat s/sample-uniform 10 10)
    (random-mat s/sample-uniform 10 10 :min 1 :max 10)"
  [pdf row col & opt]
  (i/matrix (partition col (apply pdf (* row col) opt))))

(def iv-projection-matrix
  (let [X (random-mat s/sample-uniform 100 3)
        Z (random-mat s/sample-normal  100 3 :mean 5)
        zt (i/trans Z)
        ztzi (i/solve (i/mmult zt Z))]
    (i/mmult Z ztzi zt)))

;; Bertrand (1) There are at least two firms producing homogeneous
;; (undifferentiated) products; (2) Firms do not cooperate; (3) Firms
;; compete by setting prices simultaneously; (4) Consumers buy
;; everything from a firm with a lower price. If all firms charge the
;; same price, consumers randomly select among them.

;; symmetric firms, homogenous goods
;; q1 = a - bp_1 + dbp_2
;; q2 = a - bp_2 + dbp_1

(def param-map {:a 5 :b 10 :c 2 :d 0.5})

;; (let [[a b c d] (map param-map [:a :b :c :d])]
;;   (prn (+ a b)))

(defn linear-util [y t a g]
  (+ (* y (- 1 t))
     (* a g)))

(defn gov-budget [t y r theta] (reduce + theta (* t y) r))

(defn convex-cost [r] (i/pow r 2))
(defn concave-revenue [r] (i/pow r 0.5))

;; Informed consumers, externality on uniformed consumers

(defn quality-condition
  [c-high c-low alpha]
  (let [inv-alpha (- 1 alpha)]
    (/ (+ c-high (* inv-alpha c-low)) alpha)))

(defn choose-high
  [price c-high c-low alpha]
  (let [quality-cond (quality-condition c-high c-low alpha)]
    (if (> price quality-cond) 1 0)))

;; (choose-high 16 5 2 0.4) => 1 (qual-cond = 15.5)

(defn avg-quality
  [B price c-high c-low alpha-mean alpha-sd]
  (let [partial-qual (partial choose-high price c-high c-low)
        alpha-seq (s/sample-normal B :mean alpha-mean :sd alpha-sd)]
    (s/mean (map partial-qual alpha-seq))))

(defn graph-qual
  []
  (let [qual-fn (partial avg-quality 10000 16 5 2 0.4)
        x-seq (range 0.0001 1 0.0001)]
    (i/view (c/xy-plot x-seq (map qual-fn x-seq)))))

