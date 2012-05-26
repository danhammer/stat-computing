(ns computing.test.dice
  (:use [computing.dice]
        [clojure.test]
        [clojure.contrib.generic.math-functions :only [approx=]]))

(deftest dice-problem
  "the probabilities associated with each face should all be equal to
  1/6 when the observed average is 3.5 -- the average of a fair die."
  (let [prob-one (first (prob-seq 3.5))]
    (is (approx= prob-one (/ 1 6) 0.001))))
