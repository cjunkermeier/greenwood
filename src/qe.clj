(ns qe
 (:require [clojure.core.matrix :as cmtx]
   [filters])
 (:use     clojure.core.matrix.operators))



(defn cubic-fcc
[a]
 (* 0.5 a [[-1 0 1][0 1 1][-1 1 0]]))



(defn triclinic
  ""
  [a b c alpha beta gamma]
  (vector [a 0.0 0.0]
   [(* b (cmtx/cos gamma)), (* b (cmtx/sin gamma)), 0.0]
   [(* c (cmtx/cos beta)) , (/ (* c (- (cmtx/cos alpha) (* (cmtx/cos beta) (cmtx/cos gamma)))) (cmtx/sin gamma)),
    (/ (* c (cmtx/sqrt (+ 1  (* 2 (cmtx/cos alpha) (cmtx/cos beta) (cmtx/cos gamma))
                      (* -1 (cmtx/cos alpha) (cmtx/cos alpha))
                      (* -1 (cmtx/cos beta) (cmtx/cos beta))
                      (* -1 (cmtx/cos gamma) (cmtx/cos gamma))))) (cmtx/sin gamma))]))



(defn QE-triclinic
  [celldm1 celldm2 celldm3 celldm4 celldm4 celldm5 celldm6]
  (triclinic celldm1 (* celldm2 celldm1) (* celldm3 celldm1) celldm4 celldm5 celldm6))



(defn find-band-gap
  "Used with quantum espresso"
  [bands fermienergy]
  (let [above-below (filters/multi-filter [#(> % fermienergy) #(< % fermienergy)] (flatten bands))]
    (- (apply min (first above-below)) (apply max (second above-below)))))
