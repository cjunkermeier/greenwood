
(ns poresizes
(:require [greenwood.empirical-data :as ed]
          [greenwood.math :as gmath]
          [clojure.core.reducers :as r])
(:refer-clojure :exclude [* - + == /])
(:use clojure.core.matrix
clojure.core.matrix.operators
clojure.math.combinatorics
 ultra-csv.core))

(defn MP-sphere-centers
  "Produces an evenly distributed set of points (Monkhorst-Pack like)
   in real space.
  Usage: (MP-sphere-centers 0.5 0.5 0.5 -3 -3 -3 3 3 3)."
  [stepsize-x stepsize-y stepsize-z min-x min-y min-z max-x max-y max-z]
  (clojure.math.combinatorics/cartesian-product
           (range min-x (+ max-x (* stepsize-x 0.5)) stepsize-x)
           (range min-y (+ max-y (* stepsize-y 0.5)) stepsize-y)
           (range min-z (+ max-z (* stepsize-z 0.5)) stepsize-z)))

(defn create-spheres
  "This function creates multiple spheres around the center of
   all given points. It can be used to determine how large of
   a sphere can fit inside the pore."
  [points radius-min radius-max stepsize]
  (map #(hash-map :coordinates (first %) :radius (second %))
      (clojure.math.combinatorics/cartesian-product
       points
       (range radius-min (+ radius-max (* stepsize 0.5)) stepsize))))

(defn sphere-overlapping-mol?
  "This checks to make sure that a sphere is not too
  close to any atom in the mol. If the sphere does not
  overlap with any of the atoms then this function returns
  'false'.
  sphere = (hash-map :coordinates [x y z] :radius r)"
  [mol sphere]
   (not (not-any? false?
        (map #(< (+ ((comp ed/atomic-radius ed/atomic-numbers :species) %)
                      (:radius sphere))
                    (distance (:coordinates sphere) (:coordinates %))) mol))))

(defn determine-pore-size
  "Given a mole and a set of spheres this creates a seq of
   hash-maps that gives the positions and radii of spheres
   that don't overlap with any of the atoms in the mol.
   The seq can be written out to a csv file.

   This is a good choice is you want to show the size of sphere at each point."
  [mol spheres]
    (let [max-by-coord (fn [c] (as-> (group-by :coordinates c) cc
                                 (map #(hash-map :x (first %1)
                                                 :y (second  %1)
                                                 :z (last  %1)
                                                 :radius  %2)
                    (map first cc)
                    (map (comp #(apply max %) #(map :radius %) second) cc) )))]
  (->> (r/filter #(not (sphere-overlapping-mol? mol %)) spheres)
       (into [])
       (max-by-coord))))


(defn determine-pore-size-max
  "Given a mole and a set of spheres this determines the sphere with then
   largest radii that doesn't overlap with any of the atoms in the mol."
  [mol spheres]
    (let [max-by-coord (fn [c] (as-> (group-by :coordinates c) cc
                                 (map #(hash-map :x (first %1)
                                                 :y (second  %1)
                                                 :z (last  %1)
                                                 :radius  %2)
                    (map first cc)
                    (map (comp #(apply max %) #(map :radius %) second) cc) )))]
  (->> (r/filter #(not (sphere-overlapping-mol? mol %)) spheres)
       (into [])
       (max-by-coord)
       (sort-by :radius > )
       (first ))))


(defn random-points
    [boxpoint1 boxpoint2 npoints]
    (let [f (fn [x] (gmath/random-point boxpoint1 boxpoint2))]
    (map f
       (range npoints))))

(defn max-pore-size
   [mol rmin rmax stepsize boxpoint1 boxpoint2 npoints repeats]
   (loop [previous-sphere (hash-map :radius rmin)
          n repeats
          t 0]
      (if (or (< n 0) (= t 100))  previous-sphere
          (let [m (determine-pore-size-max mol
                    (create-spheres
                      (random-points boxpoint1 boxpoint2 npoints)
                      (:radius previous-sphere) rmax stepsize))]
          (cond
            (or (nil? m) (>= (:radius previous-sphere) (:radius m)))  (recur previous-sphere (dec n) (inc t))
            (< (:radius previous-sphere) (:radius m)) (recur m repeats (inc t)))))))








;
