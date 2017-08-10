(ns greenwood.neighbors
  (:require [clojure.core.reducers :as r]
            [greenwood.basics :as basic]
            [greenwood.mol :as gmol]
            [greenwood.math :as gmath]
            [greenwood.utils :as utils]
            [clojure.set :as cset]
            [clojure.math.combinatorics :as mathcomb]
            [clojure.core.matrix :as cmat]
            [clojure.core.matrix.operators :as cmato])
  (:refer-clojure :exclude [* - + == /]))

;We can turn off the vectorz-clj implementation of core.matrix by commenting out this line.
;(set-current-implementation :vectorz)

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)



(defn- distances- [atomm]
  "mol is one time step of the xyz file. atomm is some atom."
  (let [a (:coordinates atomm)]
  (comp (map :coordinates) (map (partial cmat/distance a)))))


(defn distances [mol atomm]
  "mol is one time step of the xyz file. atomm is some atom."
  (sequence (distances- atomm) mol))


;(time (dotimes [_ 1000] (distances C8F (first C8F))))
;"Elapsed time: 144.838 msecs" - persistent vector
;(time (dotimes [_ 1000] (distances C8F (first C8F))))
;"Elapsed time: 112.447 msecs" - core.matrix/matrix


#_(defn distances-periodic
  "mol is one time step of the xyz file. atomm is some atom."
   [mol atomm projectors]
  (let [a (:coordinates atomm)]
    (mapv
  #(mapv (comp (partial distance a) (partial map + %) :coordinates) mol)
     projectors)))


(defn distances-non-atomm [mol atomm]
  "This computes the distances between atomm and all of the other atoms in mol.  Thus the length of the resulting seq will be of length (dec (count mol)).
  mol is one time step of the xyz file. atomm is some atom."
  (let [a (:coordinates atomm)]
  (doall (mapv (comp (partial cmat/distance a) :coordinates) (gmol/mol-filter-not {:pos (:pos atomm)} mol)))))



(defn all-distances
  [mol]
  (doall (mapv (partial distances-non-atomm mol) mol)))







(defn neighbors-distances
  "In this case mol is one time step of the xyz file.
If you run this and you get the error message:
#<CompilerException java.lang.IllegalArgumentException: No value supplied for key: clojure.lang.LazySeq@0
Then you need to apply xyz-parser/atom-pos to the mol."
   [mol atomm min-distance max-distance]
  (let [neigh-atoms (filter (comp #(and (> max-distance %) (< min-distance %))
                              #(cmat/distance (:coordinates %) (:coordinates atomm))) )
        neighv (map #(basic/neigh-struct (:pos %) (:species %) (cmat/distance (:coordinates %) (:coordinates atomm)) (:coordinates %)))]
(sequence (comp neigh-atoms neighv) mol)))









#_(defn neighbors-distances-periodic
  "In this case mol is one time step of the xyz file.
If you run this and you get the error message:
#<CompilerException java.lang.IllegalArgumentException: No value supplied for key: clojure.lang.LazySeq@0
Then you need to apply xyz-parser/atom-pos to the mol."
   [mol atomm min-distance max-distance projectors]
  (let [neigh-atoms (flatten (doall (filter seq (map (fn [x] (filter (comp #(and (> max-distance %) (< min-distance %))
                              #(euclidean (map + x (:coordinates %)) (:coordinates atomm))) mol)) projectors))))
        list-positions (map :pos neigh-atoms)
        list-species (map :species neigh-atoms)
        list-distances (flatten (doall (filter seq (map (fn [x] (filter #(and (> max-distance %) (< min-distance %))
                              (map #(euclidean (map + x (:coordinates %)) (:coordinates atomm)) mol))) projectors))))]
    (neigh-struct list-positions list-species list-distances)))



(defn neighbors
    "This computes the neighbors of all of the atoms in one time step."
  [mol min-distance max-distance]
  (doall (map #(assoc-in % [:neigh] (neighbors-distances mol % min-distance max-distance)) mol )))








#_(defn neighbors-periodic
  [mol min-distance max-distance projectors]
  (doall (pmap #(assoc-in % [:neigh] (neighbors-distances-periodic mol % min-distance max-distance projectors)) mol)))





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
  (as-> (overlapping-atoms mol max-distance) x
       (map (comp first sort vec) x)
       (flatten x)
       (gmol/mol-filter-not-vec :pos x mol)))





(defn nearest-atom-point
  "Returns the atom in mol closest to point-vec.
Usage:  (neartest-atom-point graphene [0 0 0])"
  [mol point-vec]
  (let [distancevec (map (comp cmat/length #(cmato/- point-vec (:coordinates %)))  mol)
        minval (apply min distancevec)
        pos (utils/positions #{minval} distancevec)]
    (gmol/mol-nth mol (first pos))))






(defn atom-angles-trios
  "Helper function to create a list of triplets, where each triplet is the list of atoms involved in a three atom angle."
       [atomm]
  (as-> (:neigh atomm) x
       (map :npos x)
       (mathcomb/combinations x 2)
       (mapv #(vector (first %) (:pos atomm) (second %)) x)))




(defn atom-angles
  "Computes all of the angles between atomm and it's neighboring atoms, where atomm is at the pivot of the angle."
  [mol atomm]
    (if (>= (count (:neigh atomm)) 2)
        (map #(basic/angle-struct % (gmath/three-point-angle (:coordinates (gmol/mol-nth mol (first %)))
                                   (:coordinates (gmol/mol-nth mol (second %)))
                                   (:coordinates (gmol/mol-nth mol (last %)))))  (atom-angles-trios atomm))
        nil))






(defn angles
  "Use this to compute the angles between the nearest neighbors of all of the
atoms."
  [mol]
  (map #(assoc-in % [:angles] (atom-angles mol %)) mol))



(defn angles-standalone
  "Use this to compute the angles between the nearest neighbors of all of the
atoms."
  [mol]
  (map #(atom-angles mol %) mol))




(defn neighbor-order
  "This allows you to specify the neighbor order from the atom specified by atom-num.
In this case atom-num is the value of :pos.  This does not take into account
periodic boundaries.  In order for this to work you will have had to run neighbors
on mol first.  This is Manhattan distance used in statistics."
  ([mol atom-num]
    (neighbor-order mol atom-num (cmato/* 2 (count mol))))
  ([mol atom-num maxnum]
  (letfn [(group [done intermediate whatsleft i]
            (let [a (vec (set (flatten (map (comp (partial map :npos) :neigh) intermediate))))]
            (if (or (= i (inc maxnum))(not (seq intermediate)))
              (flatten done)
             (group
               (concat done intermediate)
               (-> (gmol/mol-filter-vec :pos a whatsleft)
                 (flatten )
                 (gmol/update-mol-name map? i))
                (gmol/mol-filter-not-vec :pos a whatsleft)
               (inc i)))))]
       (group [] (gmol/update-mol-name (gmol/mol-filter {:pos atom-num} mol) map? 0) (gmol/mol-filter-not {:pos atom-num} mol) 1))))




(defn atoms-near-line
  "Filters a mol so that it gives atoms along a line (or within 0.1 Angstroms)."
  [mol pt1 pt2 distance]
    (gmol/mol-filter {:coordinates #(> distance (gmath/point-line-distance  pt1 pt2 %))} mol))



(defn count-bonds
  "This counts the total number of bonds in a system.  You should have ran nieghbors on the mol before running this."
  [mol]
  (->> mol
    (map :neigh)
    (map count)
    (map (partial cmato/* 0.5))
     (reduce cmato/+)))





(defn distances-distribution
  [mol max-dist]
  (loop [a (first mol)
         m (rest mol)
         d []]
    (if (empty? m)
      (concat d (filter #(> max-dist %) (distances m a)))
      (recur (first m) (rest m) (concat d (filter #(> max-dist %) (distances m a)))))))





(defn atom-atom-dipole-moment
  "This computes the dipole moment of a bond, it assumes that
   you have charge of each atom."
  [atom1 atom2]
  (cmato/* (cmato/- (:charge atom1) (:charge atom2))
    (cmato/- (:coordinates atom1) (:coordinates atom2))))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn simple-df
  [min-distance max-distance atom1 atom2]
  (#(and (> max-distance %) (< min-distance %)) (cmat/distance (:coordinates atom1) (:coordinates atom2))))


(defrecord dist-rec [^int pos1 ^double charg1 ^int pos2 ^double charg2  ^doubles dist ])
(defn dist-struct [pos1 charg1 pos2 charg2 dist]
  (->dist-rec pos1 charg1 pos2 charg2 dist))




(defn atom-distances
  "molna is the part of the mol listed after atomm."
  [df atomm molna]
    (loop [b (first molna)
         m (rest molna)
         v (transient (vector ))]
    (cond (nil? b) (persistent! v)
          (false? (df atomm b)) (recur (first m) (rest m) v)
          (true? (df atomm b)) (recur (first m) (rest m)
                                  (conj! v (dist-struct (:pos atomm) (:charge atomm)
                                                    (:pos b) (:charge b)
                                                    (cmato/- (:coordinates atomm) (:coordinates b))))))))


(defn- computation-projectors-
"This is the same function as is found in greenwood.supercell, I brought this here
  and made it private so that this name space wouldn't depend on a name space that depends on it."
  ([la1 ma2 na3]
    (mathcomb/cartesian-product (concat (range (inc la1)) (range (cmato/- la1) 0))
      (concat (range (inc ma2)) (range (cmato/- ma2) 0))
      (concat (range (inc na3)) (range (cmato/- na3) 0))))

  ([lvs la1 ma2 na3]
    (map #(cmato/+ (cmato/* (first %) (first lvs))
            (cmato/* (second %) (second lvs))
            (cmato/* (last %) (last lvs)))
      (computation-projectors- la1 ma2 na3))))



(defn s-distances
  "molna is the part of the mol listed after atomm."
  ([df mol]
  (loop [a (first mol)
         m (rest mol)
         V (transient (vector ))]
    (if (empty? m) (flatten (persistent! V))
          (let [n (atom-distances df a m)]
            (if (empty?  n)
              (recur (first m) (rest m) V)
              (recur (first m) (rest m) (conj! V n)))))))
  ([df mol lvs]
  (let [p (flatten (map #(gmol/shift % mol)
                        (cset/difference (set (computation-projectors- lvs  1 1 0))
                                         (set [[0.0 0.0 0.0]]))))]
   (loop [a (first mol)
         m (rest mol)
         V (transient [])]
    (if (empty? m) (flatten (persistent! V))
          (let [n (atom-distances df a (flatten [m p]))]
            (if (empty?  n)
              (recur (first m) (rest m) V)
              (recur (first m) (rest m) (conj! V n)))))))))





#_(defn s-distances
  "molna is the part of the mol listed after atomm."
  ([df mol]
  (loop [a (first mol)
         m (rest mol)
         V (transient (vector ))]
    (if (empty? m) (flatten (persistent! V))
          (let [n (atom-distances df a m)]
            (if (empty?  n)
              (recur (first m) (rest m) V)
              (recur (first m) (rest m) (conj! V n)))))))
  ([df mol lvs]
  (let [p (flatten (map #(gmol/shift % mol)
                        (cset/difference (set (computation-projectors- lvs  1 1 0))
                                         (set [[0.0 0.0 0.0]]))))]
   (loop [a (first mol)
         m (rest mol)
         V (transient [])]
    (if (empty? m) (flatten (persistent! V))
          (let [n (atom-distances df a (flatten [mol p]))]
            (if (empty?  n)
              (recur (first m) (rest m) V)
              (recur (first m) (rest m) (conj! V n)))))))))








(defn dipole-moment
  ([df mol]
  (->> (s-distances df mol)
      (map #(cmato/* (cmato/- (:charg1 %) (:charg2 %))  (:dist %)))
     (reduce cmato/+ [0.0 0.0 0.0] )))
  ([df mol lvs]
  (->> (s-distances df mol lvs)
      (map #(cmato/* (cmato/- (:charg1 %) (:charg2 %))  (:dist %)))
     (reduce cmato/+ [0.0 0.0 0.0] ))))






#_(require '[greenwood.xyz :as xyz])
#_(def graphene (->> (xyz/foldable-chunks "/Users/chadjunkermeier/Desktop/graphene.xyz" )
     (r/map (partial xyz/parse-xmolout-readable 5))
       (into [])
       (last )))

#_(def dffff (partial simple-df 0.1 1.08))

#_(s-distances dffff (:mol graphene))

#_(dipole-moment  dffff  (:mol graphene) (:lvs graphene))

#_(dipole-moment  dffff  (:mol graphene) )



;(s-distances dffff (:mol graphene) (:lvs graphene))

