(ns greenwood.neighbors-species
  (:require [clojure.core.reducers :as r]
            [greenwood.basics :as basic]
            [greenwood.empirical-data :as ed]
            [greenwood.mol :as jmol]
            [greenwood.math :as jmath]
            [greenwood.utils :as utils]
            [clojure.math.combinatorics :as mathcomb])
  (:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix)
  (:use clojure.core.matrix.operators))


;(use 'greenwood.testsystems :reload)

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Using user defined maximum radius
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn atom-neighbors-maxdis
   [mol atomm  max-distance]
   (let [ap (:coordinates atomm)]
    (filter (complement nil?)
            (map #(if (> max-distance (distance (:coordinates %) ap) 0.1)
          (basic/neigh-struct (:pos %) (:species %) (distance (:coordinates %) ap) (:coordinates %)))
         mol))))



#_(defn atom-neighbors-maxdis
   [mol atomm  max-distance]
   (let [ap (:coordinates atomm)]
     (loop [a (first mol)
            mmol (rest mol)
            n []]
       (if (nil? a)
         n
         (if (> max-distance (distance (:coordinates a) ap) 0.1)
             (recur (first mmol) (rest mmol)
                (into n (basic/neigh-struct (:pos a) (:species a) (distance (:coordinates a) ap) (:coordinates %))))
  (recur (first mmol) (rest mmol) n))))))





(defn neighbors-maxdis
    "This computes the neighbors of all of the atoms in one time step."
  [mol max-distance]
  (doall (map #(assoc-in % [:neigh] (atom-neighbors-maxdis mol % max-distance)) mol)))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Using user defined maximum radius -- PERIODIC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





#_(defn atom-neighbors-maxdis-periodic
   [mol atomm  max-distance projectors]
   (let [ap (:coordinates atomm)
         nei (fn [x] (map #(if (> max-distance (distance (+ x (:coordinates %)) ap) 0.1)
          (basic/neigh-struct (:pos %) (:species %) (distance (+ x (:coordinates %)) ap)))
           mol))]))



(defn atom-neighbors-maxdis-projector
   [mol atomm  max-distance projector]
   (let [ap (:coordinates atomm)]
    (filter (complement nil?)
            (map #(if (> max-distance (distance (+ (:coordinates %) projector) ap) 0.1)
          (basic/neigh-struct (:pos %) (:species %) (distance (:coordinates %) ap) (+ (:coordinates %) projector)))
         mol))))



(defn atom-neighbors-maxdis-periodic
   [mol atomm  max-distance projectors]
   (loop [p (first projectors)
          pp (rest projectors)
          neig []]
     (if (nil? p)
       neig
       (recur (first pp) (rest pp)
     (concat neig (atom-neighbors-maxdis-projector mol atomm  max-distance p))))))





(defn neighbors-maxdis-periodic
    "This computes the neighbors of all of the atoms in one time step."
  [mol max-distance projectors]
  (doall (map #(assoc-in % [:neigh] (atom-neighbors-maxdis-periodic mol % max-distance projectors)) mol)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Using User defined interactions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(def interactions [(hash-map :species "O" :inter [(hash-map :species "H" :min 1.0 :max 1.8)])
                   (hash-map :species "H" :inter [(hash-map :species "F" :min 1.9 :max 2.5)])])


;(defn )






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Using covalent radius
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn atom-neighbors
   [mol atomm]
   (let [fr (comp ed/atomic-radius ed/atomic-numbers :species)
         ar (fr atomm)
         ap (:coordinates atomm)]
    (map #(if (> (* 1.4 (+ (fr %) ar)) (distance (:coordinates %) ap) 0.1)
          (basic/neigh-struct (:pos %) (:species %) (distance (:coordinates %) ap) (:coordinates %) ))
         mol)))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Old method of computing neighbors-but still great for running other analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn- distances- [atomm]
  "mol is one time step of the xyz file. atomm is some atom."
  (let [a (:coordinates atomm)]
  (comp (map :coordinates) (map (partial distance a)))))


(defn distances [mol atomm]
  "mol is one time step of the xyz file. atomm is some atom."
  (sequence (distances- atomm) mol))



(defn neighbors-distances [mol atomm min-distance max-distance]
  "In this case mol is one time step of the xyz file.
If you run this and you get the error message:
#<CompilerException java.lang.IllegalArgumentException: No value supplied for key: clojure.lang.LazySeq@0
Then you need to apply xyz-parser/atom-pos to the mol."
  (let [neigh-atoms (filter (comp #(and (> max-distance %) (< min-distance %))
                              #(distance (:coordinates %) (:coordinates atomm))))
        neighv (map #(basic/neigh-struct (:pos %) (:species %) (distance (:coordinates %) (:coordinates atomm)) (:coordinates %) ))]
(sequence (comp neigh-atoms neighv) mol)))





(defn overlapping-atoms
  "This checks to make sure that no two atoms in a particular mol are too close to each other. The
output is a set of sets, where each of the subsets contains the atoms that are
too close."
  [mol max-distance]
  (letfn [(listset [x]
            (let [neighs (:npos (neighbors-distances mol x -0.1 max-distance))]
              (if ((comp not empty?) neighs)
                (map #(set [(:pos x) %]) neighs))))]
    (remove nil? (set (utils/flatten-n 1 (doall (pmap listset mol)))))))




(defn remove-overlapping
  "Use with care."
  [mol max-distance]
  (->> (overlapping-atoms mol max-distance)
       (map (comp first sort vec) )
       (flatten )
       (#(jmol/mol-filter-not-vec :pos % mol))))




(defn nearest-atom-point
  "Returns the atom in mol closest to point-vec.
Usage:  (neartest-atom-point graphene [0 0 0])"
  [mol point-vec]
  (let [distancevec (map (comp length #(- point-vec (:coordinates %)))  mol)
        minval (apply min distancevec)
        pos (utils/positions #{minval} distancevec)]
    (jmol/mol-nth mol (first pos))))




(defn neighbor-order
  "This allows you to specify the neighbor order from the atom specified by atom-num.
In this case atom-num is the value of :pos.  This does not take into account
periodic boundaries.  In order for this to work you will have had to run neighbors
on mol first."
  ([mol atom-num]
    (neighbor-order mol atom-num (* 2 (count mol))))
  ([mol atom-num maxnum]
  (letfn [(group [done intermediate whatsleft i]
            (let [a (vec (set (flatten (map (comp (partial map :npos) :neigh) intermediate))))]
            (if (or (= i (inc maxnum))(not (seq intermediate)))
              (flatten done)
             (group
               (concat done intermediate)
               (-> (jmol/mol-filter-vec :pos a whatsleft)
                 (flatten )
                 (jmol/update-mol-name map? i))
                (jmol/mol-filter-not-vec :pos a whatsleft)
               (inc i)))))]
       (group [] (jmol/update-mol-name (jmol/mol-filter {:pos atom-num} mol) map? 0) (jmol/mol-filter-not {:pos atom-num} mol) 1))))





(defn atoms-near-line
  "Filters a mol so that it gives atoms along a line (or within 0.1 Angstroms)."
  [mol pt1 pt2 distance]
    (jmol/mol-filter {:coordinates #(> distance (jmath/point-line-distance  pt1 pt2 %))} mol))



(defn count-bonds
  "This counts the total number of bonds in a system.  You should have ran nieghbors on the mol before running this."
  [mol]
  (->> mol
    (map :neigh)
    (map count)
    (map (partial * 0.5))
     (reduce +)))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Angles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




(defn atom-angles-trios
  "Helper function to create a list of triplets, where each triplet is the list of atoms involved in a three atom angle."
       [atomm]
  (as-> (:neigh atomm) x
       (map :npos x)
       (mathcomb/combinations x 2)
       (mapv #(vector (first %) (:pos atomm) (second %)) x)))





(defn atom-angles-trios-coords
  "Helper function to create a list of triplets, where each triplet is the list of coordinates involved in a three atom angle."
       [atomm]
  (as-> (:neigh atomm) x
       (map :ncoord x)
       (mathcomb/combinations x 2)
       (mapv #(vector (first %) (:coordinates atomm) (second %)) x)))




(defn atom-angles
  "Computes all of the angles between atomm and it's neighboring atoms, where atomm is at the pivot of the angle.
  The user must have run neighbors on mol before using this function."
  [mol atomm]
    (if (>= (count (:neigh atomm)) 2)
        (map #(basic/angle-struct %1 (jmath/three-point-angle (first %2)
                                    (second %2)
                                    (last %2))) (atom-angles-trios atomm) (atom-angles-trios-coords atomm))
        nil))





(defn angles
  "Use this to compute the angles between the nearest neighbors of all of the
atoms."
  [mol]
  (map #(assoc-in % [:angles] (atom-angles mol %)) mol))



(defn angles-stand-alone
  "Use this to compute the angles between the nearest neighbors of all of the
atoms."
  [mol]
  (map #(atom-angles mol %) mol))























