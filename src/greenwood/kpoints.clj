(ns greenwood.kpoints
  (:require [greenwood.utils :as utils]
            [greenwood.empirical-data :as ed])
  (:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix)
  (:use clojure.core.matrix.operators))




(defn generate-middle-kpoints [kstart n kend]
  "This generates a col of kpoints starting a kstart followed by n kpoints
incrementally moving linearly toward kend, but the col does not contain kend."
  (let [changevec (map #(/ % (inc n))(- kend kstart))]
    (take (inc n) (iterate #(+ 0.0 changevec %) kstart))))


(defn generate-band-kpoints [col-kpoints col-integers]
  "This computes the kpoints needed to produce a band structure.  All of the
kpoints are weighted equally. col-kpoints is a list of the high symmetry kpoints
in the order that is desired.  col-integers is a list of the number of kpoints
in between each high symmetry kpoint.  The output will be a list of 4-tuples,
with the first being the first high symmetry kpoint and the last being the last
high symmetry kpoint.

Usage: (generate-band-kpoints [Gamma K M Gamma] [10 5 8])"
  (if (= (count col-kpoints) (inc (count col-integers)))
  (let [factor (/ 1. (+ (count col-kpoints) (reduce + col-integers)))
        endpoints (partition 2 1 col-kpoints)
        kpoint-elements (inc (count (first col-kpoints)))]
    (partition kpoint-elements
    (flatten
      (interleave
        (concat
          (utils/flatten-n 1
            (map #(generate-middle-kpoints (first %1) %2 (second %1)) endpoints col-integers)) (vector (last col-kpoints)))
        (repeat factor)))))
    (println "error in calling generate-band-kpoints: (count col-kpoints) should equal (inc (count col-integers))")))





(defn MP-grid
  "Produces a Monkhorst-Pack grid of kpoints in the Brillouin-zone.

Usage:
  4x4x4 grid centered at the origin:
    (MP-grid 4 4 4) or
    (MP-grid 4 4 4 [0.0 0.0 0.0])
  4x4x4 grid centered at [0.1 0.1 0.1]:
    (MP-grid 4 4 4 [0.1 0.1 0.1])

  To include the gamma point use odd numbers:
    (MP-grid 5 7 9)"
   ([n1 n2 n3]
    (MP-grid n1 n2 n3 [0.0 0.0 0.0]))
  ([n1 n2 n3 P]
       (for [a (take n1 (iterate inc 1))
             b (take n2 (iterate inc 1))
             c (take n3 (iterate inc 1))]
        (+ P [(double (/ (- (* 2 a) n1 1) (* 2 n1)))
         (double (/ (- (* 2 b) n2 1) (* 2 n2)))
         (double (/ (- (* 2 c) n3 1) (* 2 n3)))]))))




(defn recipcal-lvs
  "Computes the recipcal lattice vectors associated to a real space lattice."
  [lvs]
  (vec (* ed/tau (cross (second lvs) (last lvs)) (/ 1 (dot (first lvs) (cross (second lvs) (last lvs)))))
       (* ed/tau (cross (last lvs) (first lvs)) (/ 1 (dot (first lvs) (cross (second lvs) (last lvs)))))
       (* ed/tau (cross (first lvs) (second lvs)) (/ 1 (dot (first lvs) (cross (second lvs) (last lvs)))))))





;(def Gamma [0 0 0])
;(def K [1/3 1/3 0])
;(def M [0.5     0.0     0.0])
;(generate-band-kpoints [Gamma K M Gamma] [10 5 8])
