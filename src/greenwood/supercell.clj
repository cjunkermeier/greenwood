(ns greenwood.supercell
  (:require [greenwood.basics :as b]
            [greenwood.contrib-math :as cmath]
            [greenwood.utils :as utils]
            [greenwood.math :as gmath]
            [greenwood.mol :as gmol]
            [greenwood.neighbors :as gngh]
            [greenwood.xyz :as xyz]
            [clojure.math.combinatorics :as combine]
            [clojure.string :as st])
  (:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix)
  (:use clojure.core.matrix.operators))





(defn cartesian->fractional
    "This function is designed to take atoms in a fractional (or internal, or primitive)
coordiante system and transform them into the real space cartesian coordinate
system.

mol is the vector of atoms
lvs is the lattice vectors."
  [mol lvs]
  (gmol/update-mol :coordinates #(mmul ((comp transpose inverse) lvs)  %) mol))




(defn fractional->cartesian
  "This function is designed to take atoms in a cell (or internal, or primitive)
coordiante system and transform them into the real space cartesian coordinate
system.
lvs are the lattice vectors of the primitive unit cell."
   [mol lvs]
    (gmol/update-mol  :coordinates
            #(mmul (transpose lvs) %) mol))



(defn scale-lat-vec
  "Rescales lattice vectors to use for interpolating structures."
  [lvs factor]
  (letfn [(scale-vec [x] (* factor x))]
    (map scale-vec lvs)))


(defn linear-interpolate-supercell
  "Use this to linear interpolate the cartesian coordinates of atoms in mol."
  [mol initial-lvs final-lvs]
  (as-> mol x
   (cartesian->fractional x initial-lvs)
   (fractional->cartesian x final-lvs)))



(defn cell-projectors
  "When only the miller indicies are given, [la1 ma2 na3], this computes cell
projectors for use with a set of atoms coordinates in the crystal (internal)
coordinates.  Thus, position = a A1 + b A2 + c A3. Where A1, A2, and A3 are
lattice vectors.

To define a supercell, with atoms in those coordinates, you only have to add
integers to the first, second or utils/third indices.

latoms is the number of unit cells in the A1 direction
matoms is the number of unit cells in the A2 direction
natoms is the number of unit cells in the A3 direction

Usage: (cell-projectors 1 1 1) => ([0 0 0])
(supercell-projectors 3 1 1) => ([0 0 0] [1 0 0] [2 0 0])

If the lattice vectors and the Miller indicies, [a1 a2 a3 la1 ma2 na3] are given
this produces cell projectors in cartesian coordinates (or whatever coordinate
system that you are using).

Usage: (cell-projectors [[2.06 -3.568 0] [2.06 3.568 0] [0 0 6.840]] 2 2 2) =>
((0.0 0.0 0.0)
  (0.0 0.0 6.84)
  (2.06 3.568 0.0)
  (2.06 3.568 6.84)
  (2.06 -3.568 0.0)
  (2.06 -3.568 6.84)
  (4.12 0.0 0.0)
  (4.12 0.0 6.84))


If the lattice parameters and the Miller indicies, [a b c alpha beta gamma  la1 ma2 na3],
are given this produces cell projectors in cartesian coordinates (or whatever coordinate
system that you are using) where alpha beta and gamma are in units of degrees.

Usage: (cell-projectors 2 2 2 (/ Math/PI 2) (/ Math/PI 2) (/ Math/PI 2) 2 2 2)"
  ([la1 ma2 na3]
    (combine/cartesian-product (range la1) (range ma2) (range na3)))

  ([lvs la1 ma2 na3]
    (map #(+ (* (first %) (first lvs))
            (* (second %) (second lvs))
            (* (last %) (last lvs)))
      (combine/cartesian-product (range la1) (range ma2) (range na3)))))




(defn computation-projectors
"This is slightly different from cell-projectors. This places boxes all round the central system.
Since the central system is the only one that is important from time-step to
time-step, it will be written first, all others will follow."
  ([la1 ma2 na3]
    (combine/cartesian-product (concat (range (inc la1)) (range (- la1) 0))
      (concat (range (inc ma2)) (range (- ma2) 0))
      (concat (range (inc na3)) (range (- na3) 0))))

  ([lvs la1 ma2 na3]
    (map #(+ (* (first %) (first lvs))
            (* (second %) (second lvs))
            (* (last %) (last lvs)))
      (computation-projectors la1 ma2 na3))))




(defn create-supercell
  "This will create a supercell whether the atoms (and projectors) are in real
space cartesian coordinates or are in crystal (internal) coordinates.  More useful
for real space coordinates.

This also could be used to make a bigger supercell out of a supercell."
  [mol projectors]
  (flatten (map #(gmol/shift % mol) projectors)))





(defn supercell
  "This will create a supercell whether the atoms (and projectors) are in real
space cartesian coordinates or are in crystal (internal) coordinates.  More useful
for real space coordinates.

This also could be used to make a bigger supercell out of a supercell."
  [mol lvs la1 ma2 na3]
  (b/unitcell  [(* la1 (first lvs)) (* ma2 (second lvs)) (* na3 (last lvs))]
             (xyz/atom-pos (flatten (map #(gmol/shift % mol) (cell-projectors lvs la1 ma2 na3))))))




(defn computation-supercell
  "This will create a supercell whether the atoms (and projectors) are in real
space cartesian coordinates or are in crystal (internal) coordinates.  More useful
for real space coordinates.

This also could be used to make a bigger supercell out of a supercell."
  [mol lvs la1 ma2 na3]
  (b/unitcell [(* (inc (* 2 la1)) (first lvs)) (* (inc (* 2 ma2)) (second lvs)) (* (inc (* 2 na3)) (last lvs))]
            (flatten (map #(gmol/shift % mol) (computation-projectors lvs la1 ma2 na3)))))






(defn create-filtered-supercell
 "filters atoms out during the creation process, thus using less memory"
 [mol projectors f]
 (flatten (map (comp f #(gmol/shift % mol)) projectors)))






(defn filtered-supercell
  "This will create a supercell whether the atoms (and projectors) are in real
space cartesian coordinates or are in crystal (internal) coordinates.  More useful
for real space coordinates.

This also could be used to make a bigger supercell out of a supercell."
  [mol lvs la1 ma2 na3 f]
  (b/unitcell  [(* la1 (first lvs)) (* ma2 (second lvs)) (* na3 (last lvs))]
             (flatten (map (comp f #(gmol/shift % mol)) (cell-projectors lvs la1 ma2 na3)))))






(defn filtered-computation-supercell
  "This will create a supercell whether the atoms (and projectors) are in real
space cartesian coordinates or are in crystal (internal) coordinates.  More useful
for real space coordinates.

This also could be used to make a bigger supercell out of a supercell."
  [mol lvs la1 ma2 na3 f]
  (b/unitcell [(* (inc (* 2 la1)) (first lvs)) (* (inc (* 2 ma2)) (second lvs)) (* (inc (* 2 na3)) (last lvs))]
            (flatten (map (comp f #(gmol/shift % mol)) (computation-projectors lvs la1 ma2 na3)))))














(defn cartesian-supercell
  "This will create a supercell in real space cartesian coordinates.  The :coordinates
of the input mol are expected to be in crystal (fractional) coordinates.  The
output is a vector, where the first element is the supercell lattice vectors and
the second element is the mol of the atoms in the supercell.

This also could be used to make a bigger supercell out of a supercell.

Usage: (create-cartesian-supercell graphene [2.461 0 0] [1.2305 2.13129 0] [0.0 0.0 0.0] 4 4 1)
Usage: (create-cartesian-supercell graphene 2.461 2.461 0 (/ Math/PI 2) (/ Math/PI 2) (/ Math/PI 3) 4 4 1)"
  [mol lvs la1 ma2 na3]
  (let [atoms (fractional->cartesian mol lvs)
        projs (cell-projectors lvs la1 ma2 na3)]
    (hash-map :lvs [(* la1 (first lvs))
             (* ma2 (second lvs))
             (* na3 (last lvs))]
              :mol (create-supercell atoms projs))))






(defn patchwork-supercell-helper-
  "This function is designed to do the work of breaking the patchwork supercell
into it's component parts, and then delete the offending atoms from the portion
that the user said not to keep.

keep-atoms-from is the string of either name-one or of name-two"
  [mol overlap keep-atoms-from]
  (let [definitely-keep (gmol/mol-filter :name #{keep-atoms-from} mol)
        possibly-keep (gmol/mol-filter-not :name #{keep-atoms-from} mol)
        flattened (flatten (map vec overlap))]
    (concat definitely-keep
      (filter (fn [x](not (eval (concat '[or] (map  #(= (:pos x) %) flattened)))))  possibly-keep))))





(defn patchwork-supercell
  "This is designed to take two supercells and piece them together in a kind of
patchwork quilt type of thing.  Thus, start off with two supercells that are the
same size, and use rules to replace some of the atoms in mol-1 with the atoms in
mol-2.
Usage: (patchwork-supercell supercell-graphene supercell-CF #(< 10 (first %)))
Usage: (patchwork-supercell supercell-graphene supercell-CF #(and (< 10 (first %))(> 20 (first %))))"
  ([mol-1 mol-2 rules]
   (patchwork-supercell mol-1 mol-2 rules 3))
  ([mol-1 mol-2 rules keep-atoms-from]
(let [name-one (str "name-one-" (rand))
      name-two (str "name-two-" (rand))
      quilt-1 (future (gmol/update-mol-name (gmol/mol-filter-not {:coordinates rules} mol-1) map? name-one))
      quilt-2 (future (gmol/update-mol-name (gmol/mol-filter {:coordinates rules} mol-2) map? name-two))
      atoms (concat @quilt-1 @quilt-2)
      overlap (gngh/overlapping-atoms atoms 0.08)]
    (gmol/update-mol-name (cond
      (and (= keep-atoms-from 1) (set? (first overlap)))
         (xyz/atom-pos (patchwork-supercell-helper- atoms overlap name-one))
      (and (= keep-atoms-from 2) (set? (first overlap)))
        (xyz/atom-pos (patchwork-supercell-helper- atoms overlap name-two))
      (and (> keep-atoms-from 2)(set? (first overlap)))
    (do (println "You have atoms that are too close together.")
     (println overlap)
    (xyz/atom-pos atoms))
      :else
      (xyz/atom-pos atoms))   #{name-one name-two})  )))






(defn drop-symmetrically-redundant-points
  "Given a set of points as col of vectors, and also given a col of computation-projectors
this will delete one any of the points that symmetricaly redundant.  This function
automatically scans for and deletes the origin projectors, [0 0 0 ...], else it would
give false positives.
Usage: (drop-symmetry-redundant-points [[0.5 0.5 0.5] [0.75 0.75 0.75] [1.5 1.5 1.5]] [[0 0 0][1 1 1][-1 -1 -1]])
  => [(0.75 0.75 0.75) (1.5 1.5 1.5)]"
  [points projectors]
  (let [new-projectors (filter #(not (gmath/origin? %)) projectors) ;gets rid of [0 0 0] projector
        new-points (map #(map gmath/round-decimal %) points)] ;may not be needed
  (letfn [(test-points [in-points out-points]
            (if (empty? in-points)
              out-points
              (let [fpoint (first in-points)
                    compare-to (utils/flatten-n 1 (map (fn [x](map #(map + x %) (rest in-points))) new-projectors))]
                (recur (rest in-points)
                  (if (some true? (map #(gmath/vectors-equal? % fpoint) compare-to))
                    out-points
                    (conj out-points fpoint))))))]
    (test-points new-points []))))





(defn drop-overlapping-sc-atoms
 "This is a helper function to get rid of atoms thatare at the same cartesian coordinates even though they are associated with different projectors.
  Usage: (drop-overlapping-sc-atoms (:mol rolled) (:lvs rolled) 0 1 0)"
  [mol lvs nx ny nz]
  (let [proatoms (create-supercell mol (computation-projectors lvs nx ny nz))
        overlap (map second (gngh/overlapping-atoms proatoms 0.4))]
    (gmol/mol-filter-not-vec :pos overlap mol)))








(defn cartesian->create-surface
  "atoms is a unit cell in cartesian coordinates.
origin is defined in the cartesian coordinates.
n1 n2 n3 are integers defining how big a supercell to create (this is the cell
                                                               from which the surface will be chisled).
lvs are the lattice vectors defining the unit cell.
C1, C3 are vectors defined to be the new x- and z-directions.
x y z are the height, width, and depth of the surface from -x/2 to x/2, -y/2 to
y/2, and -z/2 to z/2.
epsilon is a parameter determining how close to a boundary an atom can be and be
considered to be ON the boundary.

Usage: (cartesian->create-surface graphene [0 0 0] 50 50 0
         [(a-one 1.421) (a-three 1.421) [0 0 10]](a-one 1.421)[0 0 10]
         (* 4 24.61) (* 5 2.461))

Usage:  "
  ([atoms origin lvs C1 C3 x]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x 0 0])
          y 100
          z 100
          epsilon 0.01]
      (cartesian->create-surface atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x y 0])
          z 100
          epsilon 0.01]
      (cartesian->create-surface atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y z]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x y z])
          epsilon 0.01]
      (cartesian->create-surface atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin n1 n2 n3 lvs C1 C3 x y z epsilon]
    (let [f #(mmul (gmath/rotate-vec-to-axis C1 :x) %)
          alpha (Math/asin (/ (dot C3 [0 0 1]) (length C3)))
          [xx yy zz] (* -0.5 [x y z])
          [xxx yyy zzz] (* 0.5 [x y z])
          projectors (computation-projectors lvs n1 n2 n3)
          box? #(and (gmath/tolerated-gte (first %) xx epsilon)
                  (gmath/tolerated-lt (first %) xxx epsilon)
                  (gmath/tolerated-gte (second %) yy epsilon)
                  (gmath/tolerated-lt (second %) yyy epsilon)
                  (gmath/tolerated-gte (last %) zz epsilon)
                  (gmath/tolerated-lt (last %) zzz epsilon))]
      (b/unitcell [[x 0.0 0.0] [0.0 y 0.0] [0.0 0.0 z]]
      ((comp
         #(gmol/mol-filter {:coordinates box?} %)
         #(gmol/rotate-mol % [0 0 0] [0 0 1] alpha)
         #(gmol/update-mol  :coordinates f %)
         #(gmol/shift (* -1 origin) %))
        (create-supercell atoms projectors))))))






(defn create-surface-cartesian
  "atoms is a unit cell in fractional coordinates.
origin is defined in the fractional coordinates.
n1 n2 n3 are integers defining how big a supercell to create
      (this is the cell from which the surface will be chisled).
lvs are the lattice vectors defining the unit cell.
C1, C3 are vectors defined to be the new x- and z-directions.
nn1 nn2 nn3 are integers the determine how many lattice vectors to move over
(from the origin) before starting the chisiling processes.
x y z are the height, width, and depth of the surface from -x/2 to x/2, -y/2 to
y/2, and -z/2 to z/2.

Usage: (create-surface-cartesian graphene [0 0 0] 50 50 1
         (a-one 1.421) (a-three 1.421) [0 0 10](a-one 1.421)[0 0 10]
         10 10 0 (* 4 24.61) (* 5 2.461) 1)"
  ([atoms origin lvs C1 C3 x]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x 0 0])
          y 100
          z 100
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x y 0])
          z 100
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y z]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x y z])
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin n1 n2 n3 lvs C1 C3 x y z epsilon]
    (let [[A1 A2 A3] lvs
          f #(mmul (gmath/rotate-vec-to-axis C1 :x) %)
          alpha (Math/asin (/ (dot C3 [0 0 1]) (length C3)))
          cart-origin (+ (* (first origin) A1)
                        (* (second origin) A2)
                        (* (last origin) A3))
          projectors (computation-projectors n1 n2 n3)
          [xx yy zz] (* -0.5 [x y z])
          [xxx yyy zzz] (* 0.5 [x y z])
          box? #(and (gmath/tolerated-gte (first %) xx epsilon)
                  (gmath/tolerated-lt (first %) xxx epsilon)
                  (gmath/tolerated-gte (second %) yy epsilon)
                  (gmath/tolerated-lt (second %) yyy epsilon)
                  (gmath/tolerated-gte (last %) zz epsilon)
                  (gmath/tolerated-lt (last %) zzz epsilon))]
      (b/unitcell [[x 0.0 0.0] [0.0 y 0.0] [0.0 0.0 z]]
      ((comp
         #(gmol/mol-filter {:coordinates box?} %)
         #(gmol/rotate-mol % [0 0 0] [0 0 1] alpha)
         #(gmol/update-mol :coordinates f %)
         #(gmol/shift (* -1 cart-origin) %)
         #(fractional->cartesian % lvs))
        (create-supercell atoms projectors))))))












(defn create-surface-cartesian
  "atoms is a unit cell in fractional coordinates.
origin is defined in the fractional coordinates.
n1 n2 n3 are integers defining how big a supercell to create
      (this is the cell from which the surface will be chisled).
lvs are the lattice vectors defining the unit cell.
C1, C3 are vectors defined to be the new x- and z-directions.
nn1 nn2 nn3 are integers the determine how many lattice vectors to move over
(from the origin) before starting the chisiling processes.
x y z are the height, width, and depth of the surface from -x/2 to x/2, -y/2 to
y/2, and -z/2 to z/2.

Usage: (create-surface-cartesian graphene [0 0 0] 50 50 1
         (a-one 1.421) (a-three 1.421) [0 0 10](a-one 1.421)[0 0 10]
         10 10 0 (* 4 24.61) (* 5 2.461) 1)"
  ([atoms origin lvs C1 C3 x]
    (let [[nn1 nn2 nn3] (map int (* 0.6 [x 0 0]))
          y 100
          z 100
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y]
    (let [[nn1 nn2 nn3] (map int (* 3.9 [x y 0]))
          z 100
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y z]
    (let [[nn1 nn2 nn3] (map int (* 0.6 [x y z]))
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin n1 n2 n3 lvs C1 C3 x y z epsilon]
    (let [[A1 A2 A3] lvs
          f #(mmul (gmath/rotate-vec-to-axis C1 :x) %)
          alpha (Math/asin (/ (dot C3 [0 0 1]) (length C3)))
          cart-origin (+ (* (first origin) A1)
                        (* (second origin) A2)
                        (* (last origin) A3))
          projectors (computation-projectors n1 n2 n3)
          [xx yy zz] (* -0.5 [x y z])
          [xxx yyy zzz] (* 0.5 [x y z])
          box? #(and (gmath/tolerated-gte (first %) xx epsilon)
                  (gmath/tolerated-lt (first %) xxx epsilon)
                  (gmath/tolerated-gte (second %) yy epsilon)
                  (gmath/tolerated-lt (second %) yyy epsilon)
                  (gmath/tolerated-gte (last %) zz epsilon)
                  (gmath/tolerated-lt (last %) zzz epsilon))]
      (b/unitcell [[x 0.0 0.0] [0.0 y 0.0] [0.0 0.0 z]]
      ((comp
         #(gmol/mol-filter {:coordinates box?} %)
         #(gmol/rotate-mol % [0 0 0] [0 0 1] alpha)
         #(gmol/update-mol :coordinates f %)
         #(gmol/shift (* -1 cart-origin) %)
         #(fractional->cartesian % lvs))
        (create-supercell atoms projectors))))))





















(defn define-cell
  "This defines the bounding box (parallelepiped) of a unit cell.
  The output is a vector where the first two elements are opposite
  corners of the bounding box and the last three elements are the
  unit gmath/normal-vectors of the sides of the parallelepiped.  The
  normal vectors point into the box for those sides (planes) that
  have the point [0 0 0] in them."
  ([lvs]
   (define-cell lvs [0 0 0]))
  ([lvs offset]
   (let [[l1 l2 l3] lvs]
     [offset
        (+ l1 l2 l3 offset)
        (gmath/normal-vector l2 l3)
        (gmath/normal-vector l3 l1)
        (gmath/normal-vector l1 l2)])))




(defn within-cell?
  "This determines if point lies within the parallelepiped defined by cell.
  cell is the output of define-cell.  This only really works for cubic cells."
  [cell point]
  (let [m (- point (first cell))
        n (- point (second cell))
        p? (partial <= 0 )
        n? (partial >= 0 )]
    (and (p? (dot m (nth cell 2)))
         (p? (dot m (nth cell 3)))
         (p? (dot m (nth cell 4)))
         (n? (dot n (nth cell 2)))
         (n? (dot n (nth cell 3)))
         (n? (dot n (nth cell 4))))))




(defn within-cell??
  "This determines if a point, P, lies within the parallelepiped defined by lvs.
  cell is the output of define-cell.

  This approach is about the property of parallelopipeds since they are such
  nice regular objects:
  The three sides emanating from a vertex will form the basis of the 3D vector
  space. If you choose them as positive x, y, z axes, then any point can be
  expressed as a linear combination of them.  A point lies inside iff all the
  coefficients are positive and are less than lengths of the corresponding sides.
  http://www.mathworks.com/matlabcentral/newsreader/view_thread/292833"
  [lvs V P]
  (let [[A1 A2 A3] (first lvs)
        [B1 B2 B3] (second lvs)
        [C1 C2 C3] (last lvs)
        [V1 V2 V3] V
        [P1 P2 P3] P
        a (/ (+ (* -1 B3 C2 P1) (* B2 C3 P1) (* B3 C1 P2) (* -1 B1 C3 P2) (* -1 B2 C1 P3) (* B1 C2 P3) (* B3 C2 V1) (* -1 B2 C3 V1) (* -1 B3 C1 V2) (* B1 C3 V2) (* B2 C1 V3) (* -1 B1 C2 V3))
           (+ (* -1 A3 B2 C1)  (* A2 B3 C1)  (* A3 B1 C2) (* -1 A1 B3 C2) (* -1 A2 B1 C3) (* A1 B2 C3)))
        b (/ (+ (* A3 C2 P1) (* -1 A2 C3 P1) (* -1 A3 C1 P2) (* A1 C3 P2) (* A2 C1 P3) (* -1 A1 C2 P3) (* -1 A3 C2 V1) (* A2 C3 V1) (* A3 C1 V2) (* -1 A1 C3 V2) (* -1 A2 C1 V3) (* A1 C2 V3))
           (+ (* -1 A3 B2 C1) (* A2 B3 C1) (* A3 B1 C2) (* -1 A1 B3 C2) (* -1 A2 B1 C3) (* A1 B2 C3)))
        c (/ (+ (* -1 A3 B2 P1) (* A2 B3 P1) (* A3 B1 P2) (* -1 A1 B3 P2) (* -1 A2 B1 P3) (* A1 B2 P3) (* A3 B2 V1) (* -1 A2 B3 V1) (* -1 A3 B1 V2) (* A1 B3 V2) (* A2 B1 V3) (* -1 A1 B2 V3))
           (+ (* -1 A3 B2 C1) (* A2 B3 C1) (* A3 B1 C2) (* -1 A1 B3 C2) (* -1 A2 B1 C3) (* A1 B2 C3)))]
        (every? #(and (gmath/tolerated-gte % 0.0) (gmath/tolerated-gte 1.0 %))  [a b c])))




(defn create-surface-cartesian-filtered
  "atoms is a unit cell in fractional coordinates.
origin is defined in the fractional coordinates.
n1 n2 n3 are integers defining how big a supercell to create
      (this is the cell from which the surface will be chisled).
lvs are the lattice vectors defining the unit cell.
C1, C3 are vectors defined to be the new x- and z-directions.
nn1 nn2 nn3 are integers the determine how many lattice vectors to move over
(from the origin) before starting the chisiling processes.
x y z are the height, width, and depth of the surface from -x/2 to x/2, -y/2 to
y/2, and -z/2 to z/2.

Usage: (create-surface-cartesiann graphene [0 0 0] 50 50 1
         (a-one 1.421) (a-three 1.421) [0 0 10](a-one 1.421)[0 0 10]
         10 10 0 (* 4 24.61) (* 5 2.461) 1)"
  ([atoms origin lvs C1 C3 x]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x 0 0])
          y 100
          z 100
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x y 0])
          z 100
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin lvs C1 C3 x y z]
    (let [[nn1 nn2 nn3] (map #(int (* 0.6 %)) [x y z])
          epsilon 0.01]
      (create-surface-cartesian atoms origin nn1 nn2 nn3 lvs C1 C3 x y z epsilon)))

  ([atoms origin n1 n2 n3 lvs C1 C3 x y z epsilon]
    (let [[A1 A2 A3] lvs
          f #(mmul (gmath/rotate-vec-to-axis C1 :x) %)
          alpha (Math/asin (/ (dot C3 [0 0 1]) (length C3)))
          cart-origin (+ (* (first origin) A1)
                        (* (second origin) A2)
                        (* (last origin) A3))
          projectors (computation-projectors n1 n2 n3)
          [xx yy zz] (* -0.5 [x y z])
          [xxx yyy zzz] (* 0.5 [x y z])
          box? #(and (gmath/tolerated-gte (first %) xx epsilon)
                  (gmath/tolerated-lt (first %) xxx epsilon)
                  (gmath/tolerated-gte (second %) yy epsilon)
                  (gmath/tolerated-lt (second %) yyy epsilon)
                  (gmath/tolerated-gte (last %) zz epsilon)
                  (gmath/tolerated-lt (last %) zzz epsilon))]
        (create-filtered-supercell atoms projectors
                          (comp
         (partial gmol/mol-filter {:coordinates box?})
         #(gmol/rotate-mol % [0 0 0] [0 0 1] alpha)
         #(gmol/update-mol :coordinates f %)
         #(gmol/shift (* -1 cart-origin) %)
         #(fractional->cartesian % lvs))))))



(defn shift->scell
  "This will automatically gmol/shift the mol so that it will fit inside the surface
unit cell that GULP uses, assuming, of course, that the user gave GULP the correct
lattice parameters for mol."
  [cell-x-size cell-y-size mol]
  (let [extrema (gmol/min-max-coordinates mol)
        mol-size (map #(Math/abs (reduce - %)) (partition 2 extrema))
        x-shift (- cell-x-size (first mol-size))
        y-shift (- cell-y-size (second mol-size))]
    (if (or (neg? x-shift) (neg? y-shift))
       (throw (Exception. "The mol's size is larger than cell-x/y-size."))
    ((comp
    #(gmol/shift (* -1 [(first extrema) (nth extrema 2) 0]) %)
      #(gmol/shift [(/ x-shift 2) (/ y-shift 2) 0] %))
      mol))))


(defn define-shift->scell
  "Like shift->scell this function will determine the shift needed to move
  a mol to place it within a gulp cell, but instead of moving the mol it
  will write out the amount it should be shifted."
  [cell-x-size cell-y-size mol]
  (let [extrema (gmol/min-max-coordinates mol)
        mol-size (map #(Math/abs (reduce - %)) (partition 2 extrema))
        x-shift (- cell-x-size (first mol-size))
        y-shift (- cell-y-size (second mol-size))]
    (if (or (neg? x-shift) (neg? y-shift))
       (throw (Exception. "The mol's size is larger than cell-x/y-size."))
   (+ (* -1 [(first extrema) (nth extrema 2) 0])
       [(/ x-shift 2) (/ y-shift 2) 0] ))))





(defn multiples
  "Given the lengths of lattice vectors, in the same direction, from two different unit
  cells this helper function determines how many of each unit cells are needed in that direction
  for the supercells to be the same size.
  a and b are real numbers"
  ([a b]
  (let [f (comp count second #(st/split % #"\.") str double)
        t (apply * (take (apply max (map f  [a b])) (repeat 10)))
        abc (map (comp int (partial * t)) [a b])]
    (map #(/ (apply cmath/lcm abc)  %) abc)))
  ([a b & more]
  (let [f (comp count second #(st/split % #"\.") str double)
        t (apply * (take (apply max (map f (flatten [a b more]))) (repeat 10)))
        abc (map (comp int (partial * t)) (flatten [a b more]))]
    (map #(/ (apply cmath/lcm abc)  %) abc))))




(defn incompatible-layers
  "This is used to "
  ([lvs1 lvs2]
    (let [lvss  [lvs1 lvs2]]
       (map #(hash-map :x %1 :y %2 :z %3)
           (apply multiples (map (comp gmath/round-decimal length first) lvss))
           (apply multiples (map (comp gmath/round-decimal length second) lvss))
           (apply multiples (map (comp gmath/round-decimal length last) lvss)))))
  ([lvs1 lvs2 & more]
    (let [lvss  [lvs1 lvs2 more]]
       (map #(hash-map :x %1 :y %2 :z %3)
           (apply multiples (map (comp gmath/round-decimal length first) lvss))
           (apply multiples (map (comp gmath/round-decimal length second) lvss))
           (apply multiples (map (comp gmath/round-decimal length last) lvss))))))































