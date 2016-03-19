(ns greenwood.solution
  (:require [greenwood.empirical-data :as ed]
            [greenwood.math :as jmath]
            [greenwood.mol :as jmol]
            [greenwood.supercell :as sc]
            [greenwood.xyz :as xyz])
  (:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix)
  (:use clojure.core.matrix.operators))



(def ^:const closest-radius 5.2)


(defn rand-placement
  "Places mol at a random point within the box described by the lvs.
  This does not check to see if the mol is completely within the box."
  [lvs mol]
  (->> mol
       (jmol/mol-center)
       (jmol/random-rotate-mol )
       (jmol/shift (mmul lvs (jmath/rand-point-fractional)) )))








(defn overlapping-mols?
  "This checks to make sure that none of the atoms in mol1 are too
  close to any atom in mol2. A true value means that the mols overlap."
  [min-dist mol1 mol2]
  (let [f :coordinates]
 (not (not-any? false?
    (map
     (fn [a]
       (not-any? true?
        (map #(> (abs min-dist) (length (map - (f a) (f %))))  mol1))) mol2)))))



(defn- excluded-region-
  "This determines if a point p is within, without, or on some volume, surface, or line.
  Usage: (excluded-region- [1 0 2 1 0 2 1 0 2 > 4] [2 2 0]) => true"
  [v p]
  (let [[ax xn xe ay yn ye az zn ze equ nmbr] v]
   (#(equ (+ (* ax (pow (- (first %) xn) xe)) (* ay (pow (- (second %) yn) ye)) (* az (pow (- (last %) zn) ze))) nmbr) p)))



(defn within?
 "This determines if a point p is within, without,
  or on some volume(s), surface(s), line(s),
  or some combination of such.  The volumes, surfaces,
  and lines are defined by the col, v,
  which can either be a col of the form [ax xe ay ye az ze equ nmbr]
  or a col were each element is a [ax xe ay ye az ze equ nmbr].
  A few examples will clearify this.

  The cols are used to define the parameters used in the equation
  ax x^{xe} + ay y^{ye} + az z^{ze} equ nmber,
  where equ is one of { <=, <, =, >=, > }.

  Usage: (within? [1 0 2 1 0 2 1 0 2 > 4] [2 2 0]) => true
  Usage: (within? [[1 0 2 1 0 2 1 0 2 > 4][1 0 2 1 0 2 1 0 2 < 9]] [2 2 0]) => true"
  [v p]
  (if (number? (first v))
       (excluded-region- v p)
     (every? identity (map #(excluded-region- % p) v))))





(defn sol-density
  "Gives the density of a solution in units of kg/liter.

  previous mol is some structure(s) that is(are) already within
  the bounding box defined by lvs (i.e. a surface).
  sol-mol contains all of the atoms/molecules that make us the solution."
  [lvs previous-mol sol-mol]
  (let [excld-vol nil
        prev-vol (->> previous-mol
                      (map :species)
                      (map ed/atomic-numbers)
                      (map ed/atomic-radius)
                      (map #(* 4/3 ed/pi (pow %)))
                      (reduce + 0.0))
        mass (->> sol-mol
                      (flatten)
                      (map :species)
                      (map ed/atomic-numbers)
                      (map ed/atomic-mass)
                      (reduce + 0.0))
        conversion-factor 1.6605387831627259]
  (str "The effective density of the solution in the cell is "
                (* conversion-factor (/ mass (- (jmath/lvs-volume lvs) prev-vol))) " g/ml.")))






(defn density->nmol
 "Determines how many molecules you need to get the desired density.
 desired-density should be in g/ml.
 This currently assumes that LVS is for a rectangular cell."
 [vmols lvs desired-density]
 (let [f #(reduce + (map (comp ed/atomic-mass ed/atomic-numbers :species) %))
       mass (reduce + (map f vmols))
       conversion-factor 1.6605387831627259]
       (floor (* desired-density (jmath/lvs-volume lvs) (/ 1 (* conversion-factor mass))))))


(defn- mol-radius-
  "Helper function to determine automate the determination of the
  minimum allowed spacing between mols. It might need to be changed.
  What would probably work better is something based on bond order."
  [mol]
   (+ 2.1 (apply max (map #(distance [0 0 0] (:coordinates %))(jmol/mol-center mol)))))



(defn solution
  "Use this to create a solution, contained within the box produced by
  the lattice vectors, lvs.  The solution is made up of one type of mol.

  already-there-mol is in case there is a mol in the box already (ie. a surface).

  lvs is the lattice vectors of the box.

  the box is assumed to be located with one of its corners at the origin and three
  of its edges on the pos-x-, pos-y-, pos-z-axes: offset should be a 3-tupple that
  moves the corner off of the origin.

  n-new-mol are the number of new mols to be placed in the box. In some cases
  the function will not be able to accommodate all of the new mols in the box
  and it will only privide as many as it can; a standard out will specify how
  many new mols were placed in the box.

  new-mol is the mol that you want replicated and randomly placed (random position
  and random rotation) within the box.

  (solution [] lvs-ts [0 0 0] 400 h2o-ts)"
  [already-there-mol lvs offset n-new-mol new-mol]
  (let [cell (sc/define-cell lvs offset)
         f (partial sc/within-cell? cell)
        spacing (mol-radius- new-mol)
        mols (map (partial rand-placement lvs) (repeat new-mol))]
    (loop [new (first mols)
           rst (rest mols)
           sol already-there-mol
           m 0
           n 0]
      (if (or (== m n-new-mol) (== n (* 100 n-new-mol)))
       (do (println "Number of molecules in solution = " m) (jmol/mol-center sol))
        (recur (first rst)
               (rest rst)
          (if (and (every? (comp f :coordinates) new)
               (false? (overlapping-mols? spacing sol new)))
               (concat sol (jmol/update-mol-name new map? m))
               sol)
          (if (and (every? (comp f :coordinates) new)
              (false? (overlapping-mols? spacing sol new)))
              (inc m)
              m)
       (inc n))))))




(defn multi-solution
  "This works the same as solution but with new-mol now being a seq of mols.

  To specify a ratio of mols populate the seq new-mol with the relative number
  of copies that you want.  For example, if you want a ratio of two h2o-ts
  molecules to three DMMP molecules you would set new-mol to be
  [h2o-ts h2o-ts DMMP DMMP DMMP].

  Currently this doesn't test to see if it has kept the ratio of mols that you
  asked for.

  Usage: (multi-solution [] lvs-ts [0 0 0] 400 [water water DMMP DMMP DMMP])"
  [already-there-mol lvs offset n-new-mol new-mol]
  (let [cell (sc/define-cell lvs offset)
         f (partial sc/within-cell? cell)
        mols (map (partial rand-placement lvs) (cycle new-mol))]
    (loop [new (first mols)
           rst (rest mols)
           sol already-there-mol
           m 0
           n 0]
      (if (or (== m n-new-mol) (== n (* 100 n-new-mol)))
       (do (println "Number of molecules in solution = " m) (xyz/atom-pos sol))
        (recur (first rst)
               (rest rst)
          (if (and (every? (comp f :coordinates) new)
               (false? (overlapping-mols? (mol-radius- new) sol new)))
               (concat sol (jmol/update-mol-name new map?  m))
               sol)
          (if (and (every? (comp f :coordinates) new)
              (false? (overlapping-mols? (mol-radius- new) sol new)))
              (inc m)
              m)
       (inc n))))))





(defn solution-without
  "Use this to create a solution, contained within the box produced by
  the lattice vectors, lvs.  The solution is made up of one type of mol.

  already-there-mol is in case there is a mol in the box already (ie. a surface).

  lvs is the lattice vectors of the box.

  the box is assumed to be located with one of its corners at the origin and three
  of its edges on the pos-x-, pos-y-, pos-z-axes: offset should be a 3-tupple that
  moves the corner off of the origin.

  n-new-mol are the number of new mols to be placed in the box. In some cases
  the function will not be able to accommodate all of the new mols in the box
  and it will only privide as many as it can; a standard out will specify how
  many new mols were placed in the box.

  new-mol is the mol that you want replicated and randomly placed (random position
  and random rotation) within the box.

  (solution [] within-fn lvs-ts [0 0 0] 400 h2o-ts)"
  [already-there-mol within-fn lvs offset n-new-mol new-mol]
  (let [cell (sc/define-cell lvs offset)
         f (partial sc/within-cell? cell)
        spacing (mol-radius- new-mol)
        mols (map (partial rand-placement lvs) (repeat new-mol))]
    (loop [new (first mols)
           rst (rest mols)
           sol already-there-mol
           m 0
           n 0]
      (if (or (== m n-new-mol) (== n (* 100 n-new-mol)))
       (do (println "Number of molecules in solution = " m) (xyz/atom-pos sol))
        (recur (first rst)
               (rest rst)
          (if (and (every? (comp f :coordinates) new)
               (false? (overlapping-mols? spacing sol new))
                   (every? false? (map (comp within-fn :coordinates) new)))
               (concat sol (jmol/update-mol-name new map? m))
               sol)
          (if (and (every? (comp f :coordinates) new)
              (false? (overlapping-mols? spacing sol new))
                   (every? false? (map (comp within-fn :coordinates) new)))
              (inc m)
              m)
       (inc n))))))










(defn solution
  "Use this to create a solution, contained within the box produced by
  the lattice vectors, lvs.  The solution is made up of one type of mol.

  already-there-mol is in case there is a mol in the box already (ie. a surface).

  lvs is the lattice vectors of the box.

  the box is assumed to be located with one of its corners at the origin and three
  of its edges on the pos-x-, pos-y-, pos-z-axes: offset should be a 3-tupple that
  moves the corner off of the origin.

  n-new-mol are the number of new mols to be placed in the box. In some cases
  the function will not be able to accommodate all of the new mols in the box
  and it will only privide as many as it can; a standard out will specify how
  many new mols were placed in the box.

  new-mol is the mol that you want replicated and randomly placed (random position
  and random rotation) within the box.

  (solution [] lvs-ts [0 0 0] 400 h2o-ts)"
  [already-there-mol lvs offset n-new-mol new-mol]
  (let [f (partial sc/within-cell?? lvs offset)
        spacing (mol-radius- new-mol)
        mols (map (partial rand-placement lvs) (repeat new-mol))]
    (loop [new (first mols)
           rst (rest mols)
           sol already-there-mol
           m 0
           n 0]
      (if (or (== m n-new-mol) (== n (* 100 n-new-mol)))
       (do (println "Number of molecules in solution = " m) (xyz/atom-pos sol))
        (recur (first rst)
               (rest rst)
          (if (and (every? (comp f :coordinates) new)
               (false? (overlapping-mols? spacing sol new)))
               (concat sol (jmol/update-mol-name new map? m))
               sol)
          (if (and (every? (comp f :coordinates) new)
              (false? (overlapping-mols? spacing sol new)))
              (inc m)
              m)
       (inc n))))))



