(ns greenwood.mol
  (:require ;[clojure.core.reducers :as r]
            [greenwood.utils :as utils]
            [greenwood.math :as jmath]
            [greenwood.empirical-data :as ed]
            [clojure.set :as cset])
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :as cmatrix])
  (:use clojure.core.matrix.operators))



(defn- mol-filter-pred-
""
[keyvals]
  (let [f #(cond
          (set? %2)
            (comp (partial cset/subset? %2) %1)
          (fn? %2)
            (comp %2 %1)
          :else
            (comp (partial = %2) %1))]
    (apply utils/every-p? (doall (map (partial apply f ) keyvals)))))







(defn mol-filter
  "Filters through a col of maps (mol) selecting the elements of
the mol that has keywords (kw) that point to the value.  If the value of the
element is supposed to be within a set, place it within a set (i.e. #{someval})
when calling mol-filter. Val may also be a predicate, which will act on the
elements associated with the kw.

Returns a transducer when no colmap is provided.


Usage:  Suppose a val, bas, is defined as
(def bas
[{:species 8, :coordinates [0 0 0], :name #{:cheesehead}}
 {:species 1, :coordinates [0 0 0.96], :name nil}
 {:species 1, :coordinates [0 0.96 0], :name nil}])

Then (mol-filter {:species 1} bas) returns
({:species 1, :coordinates [0 0 0.96]}
  {:species 1, :coordinates [0 0.96 0]})

Usage:  (mol-filter {:name #{:cheesehead}} bas)
returns ({:species 8, :coordinates [0 0 0], :name #{:cheesehead}})

Usage:  (mol-filter {:coordinates #(> 0.5 (last %))} bas)
returns ({:species 1, :coordinates [0 0 0.96], :name nil}).

Usage: (mol-filter {:species 1 :coordinates #(zero? (second %))} bas)
returns ({:species 1, :coordinates [0 0 0.96], :name nil})"
  ([keyvals]
   (let [pred (mol-filter-pred- keyvals)]
    (filter pred)))
  ([keyvals mol]
   (let [pred (mol-filter-pred- keyvals)]
     (filter pred mol))))



(defn mol-filter-not
  "Filters through a col of maps (mol) selecting the elements of
the mol that has keywords (kw) that DO NOT point to the value.  If the value of the
element is supposed to be within a set, place it within a set (i.e. #{someval})
when calling mol-filter.  Val may also be a predicate, which will act on the
elements associated with the kw.

Usage:  Suppose a val, bas, is defined such that
user=>bas
[{:species 8, :coordinates [0 0 0], :name #{:cheesehead}}
 {:species 1, :coordinates [0 0 0.96], :name nil}
 {:species 1, :coordinates [0 0.96 0], :name nil}]

Then (mol-filter-not :species 1 bas) returns
({:species 8, :coordinates [0 0 0], :name #{:cheesehead}})

Usage: (mol-filter-not :name #{:cheesehead} bas)
returns ({:species 1, :coordinates [0 0 0.96]}
  {:species 1, :coordinates [0 0.96 0]})

Usage: (mol-filter-not :coordinates #(> 0.5 (last %)) bas)
returns ({:species 8, :coordinates [0 0 0], :name #{:cheesehead}}
          {:species 1, :coordinates [0 0.96 0], :name nil})"
  ([keyvals]
   (let [pred (mol-filter-pred- keyvals)]
    (filter (complement pred))))
  ([keyvals mol]
   (let [pred (mol-filter-pred- keyvals)]
     (filter (complement pred) mol))))



(defn mol-filter-vec
  "Filters through a col of maps (mol) selecting the elements of
the mol that has keywords (kw) that point to the value.  If the value of the
element is supposed to be within a set, place it within a set (i.e. #{someval})
when calling mol-filter. Val may also be a predicate, which will act on the
elements associated with the kw.

Usage (mol-filter-vec :pos [0 2] bas)
returns ({:species 8, :coordinates [0 0 0], :name #{:cheesehead}}{:species 1, :coordinates [0 0.96 0], :name nil :pos 2})."
  ([kw values]
   (filter (apply utils/some-p? (map #(mol-filter-pred- {kw %}) values))))
  ([kw values mol]
   (filter (apply utils/some-p? (map #(mol-filter-pred- {kw %}) values)) mol)))



(defn mol-filter-not-vec
  "Filters through a col of maps (mol) selecting the elements of
the mol that has keywords (kw) that point to the value.  If the value of the
element is supposed to be within a set, place it within a set (i.e. #{someval})
when calling mol-filter. Val may also be a predicate, which will act on the
elements associated with the kw.

Usage (mol-filter-vec :pos [0 2] bas)
returns ({:species 8, :coordinates [0 0 0], :name #{:cheesehead}}{:species 1, :coordinates [0 0.96 0], :name nil :pos 2})."
  ([kw values]
   (filter (apply utils/not-any-p? (map #(mol-filter-pred- {kw %}) values))))
  ([kw values mol]
   (filter (apply utils/not-any-p? (map #(mol-filter-pred- {kw %}) values)) mol)))







(defn mol-nth
  "Returns the atom at the :pos index. Returns nil if atom is not found."
  [mol ^long index ]
  (loop [i (first mol)
         m (rest mol)]
    (cond
      (== (:pos i) index )
        i
      (empty? m)
        nil
     :else
      (recur (first m) (rest m)))))





(defn take-mol-by-pos
  "Given a col containing a listing of the pos to keep, this will select create
a new mol from the input mol."
  [mol col]
  (flatten (map #(mol-nth mol %) col)))





(defn min-max-coordinates
  "Sometimes you need to know what the bounding values of your coordinate system
are."
  [mol]
  (let [x (future (sort (map (comp first :coordinates) mol)))
        y (future (sort (map (comp second :coordinates) mol)))
        z (future (sort (map (comp last :coordinates) mol)))]
  (vector (first @x) (last @x) (first @y) (last @y) (first @z) (last @z))))



(defn update-mol
 "Maps f across the key of mol, returning an 'updated' mol."
  ([key f]
  (map #(update-in % [key] f)))
  ([key f mol]
  (map #(update-in % [key] f) mol)))



(defn col->mol
  "Replaces the value of key with the elements of col.  There needs
to be the same number of elements in col and mol."
  ([key col]
    (map #(assoc-in %2 [key] %1) col))
  ([key col mol]
    (map #(assoc-in %2 [key] %1) col mol)))



(defn col->mol+
  "Adds another element to the value of key with the elements of col.
There needs to be the same number of elements in col and mol."
  ([key col]
    (map #(assoc-in %2 [key] %1) col))
  ([key col mol]
    (map #(assoc-in %2 [key] %1) col mol)))







#_(defn dfdf
  ""
  [keypred keyf mol]
  (map #(if ((mol-filter-pred- keypred) %)
          ((fn [x](update-in x [(first keyf)] (second keyf))) %)
          %) mol))




(defn dfdf
  "Given a mol dfdf associates into the key of keyval the val whenever the predicate of keypred is true.
  Usage: (dfdf {:neigh empty?} {:species 'B'} mol)"
  ([keypred keyval ]
  (map #(if ((mol-filter-pred- keypred) %)
          (assoc-in % [(first %)] (second %))
          %) ))
  ([keypred keyval mol]
  (map #(if ((mol-filter-pred- keypred) %)
          (assoc-in % [(first keyval)] (second keyval))
          %) mol)))




(defn dfdf
  "Given a mol dfdf associates into the key of keyval the val whenever the predicate of keypred is true.
  Usage: (dfdf {:neigh empty?} mol)"

  ([keypred mol]
  (map #(if ((mol-filter-pred- keypred) %)
          (assoc-in % [:species] "N")
          %) mol)))













(defn find-assoc-in
  "Given a vector of maps (or a vector of records),mapvec, this will assoc-in a value,
  update-val, associated with update-key into a map if the value test-val of the key test-key.
  Usage: (def aa [{:a 0 :b 0 :c 0} {:a 50 :b 0 :c 0}])
              (def bb  [50 20])
              (find-assoc-in [:a :c] bb aa) => ({:a 0 :c 0 :b 0} {:a 50 :c 20 :b 0})"
  ([[test-key update-key] [test-val update-val] mapvec]
  (reduce
   (fn [idx v] (if (= (get-in v [idx test-key]) test-val)
                   (assoc-in v [idx update-key] update-val)
                 v))
   (iterate inc 0)
   (vec mapvec))))







;This needs to be turned into a transducer
#_(defn find-assoc-in
  "Given a vector of maps (or a vector of records),mapvec, this will assoc-in a value,
  update-val, associated with update-key into a map if the value test-val of the key test-key.
  Usage: (def aa [{:a 0 :b 0 :c 0} {:a 50 :b 0 :c 0}])
              (def bb  [50 20])
              (find-assoc-in [:a :c] bb aa) => ({:a 0 :c 0 :b 0} {:a 50 :c 20 :b 0})"
  [[test-key update-key] [test-val update-val] mapvec]
  (reduce
   (fn [v idx] (if (= (get-in v [idx test-key]) test-val)
                   (assoc-in v [idx update-key] update-val)
                 v))
   (vec mapvec)
   (range (count mapvec))))


;Does not transduce
(defn find-assoc-in
  "Given a vector of maps (or a vector of records),mapvec, this will assoc-in a value,
  update-val, associated with update-key into a map if the value test-val of the key test-key.
  Usage: (def aa [{:a 0 :b 0 :c 0} {:a 50 :b 0 :c 0}])
              (def bb  [50 20])
              (find-assoc-in [:a :c] bb aa) => ({:a 0 :c 0 :b 0} {:a 50 :c 20 :b 0})"
  ([[test-key update-key] [test-val update-val] mapvec]
  (reduce
   (fn [idx v] (if (= (get-in v [idx test-key]) test-val)
                   (assoc-in v [idx update-key] update-val)
                 v))
   (range (count mapvec))
   (vec mapvec))))









(defn find-assoc-in-vec
  "Given a vector of maps (or a vector of records),mapvec, this will assoc-in a value,
  update-val, associated with update-key into a map if the value test-val of the key test-key.
  Usage: (def aa [{:a 10 :b 0 :c 0} {:a 50 :b 0 :c 0}])
              (def bb  [[50 20][10 4]])
              (find-assoc-in [:a :c] bb aa) => ({:a 4 :c 0 :b 0} {:a 50 :c 20 :b 0})"
  [[test-key update-key] test-update-vals mapvec]
  (loop [f (first test-update-vals)
         r (rest test-update-vals)
         m mapvec]
    (if (empty? r)
      (find-assoc-in [test-key update-key] f m)
      (recur (first r) (rest r) (find-assoc-in [test-key update-key] f m)))))










(defn- flatten-set
  "This is like flatten but works on a set, thus (flatten-set #{:a :b #{:c :d}})
returns #{:a :b :c :d}."
  [x]
  (set (flatten (map #(if (set? %) (vec %) %) x))))

;Does not transduce
(defn update-mol-name
  "This will add a value (be it a number, string, or keyword) to :name.
The assumption is that you will only want to add values to :name based on some
predicate, but that you might want to erase a value or values based on the
circumstances.  Thus there are three ways to delete the values.

The value associated with the key name will be a set.  Thus, the same
atom in mol can have several names attached to it.  To add a name value you do
the following:
Usage: (update-mol-name some-mol some-pred some-name)

If you want to delete all of the name values of all of the atoms, you do the this:
Usage: (update-mol-name some-mol)

If you want to delete a particular name from all of the atoms (it doesn't matter
                                                                if a particular atom actually has that name), you should do the following:
Usage (update-mol-name some-mol some-name)

Finally, if you want to remove a name based on some predicate you do this:
Usage (update-mol-name some-mol some-pred some-name :remove)"
  ([mol]
    (map #(assoc-in % [:name] nil) mol))
  ([mol removed-name]
    (map #(assoc-in % [:name] (cset/difference (:name %) (if (set? removed-name)
                                                      removed-name
                                                      #{removed-name}))) mol))
  ([mol pred added-name]
    (map #(if (pred %) (if (:name %)
                         (assoc-in % [:name] (flatten-set (conj (:name %) added-name)))
                         (assoc-in % [:name] (flatten-set #{added-name})))
            (assoc-in % [:name] (:name %))) mol))
  ([mol pred removed-name k]
            (map #(if (pred %)
                    (assoc-in % [:name] (cset/difference (:name %) (if (set? removed-name)
                                                                removed-name
                                                                #{removed-name})))
                    (assoc-in % [:name] (:name %))) mol)))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




(defn coord-average
  "Averages the x,y,z coordinates of coords independently,
returning the averages as a new coordinate."
  [mol]
  (map jmath/average (utils/transpose (map :coordinates mol))))



(defn shift
  "Displaces all atoms in mol by the vector pnt."
  ([pnt]
   (update-mol :coordinates #(+ % pnt)))
  ([pnt mol]
    (update-mol :coordinates #(+ % pnt)  mol)))



;Does not transduce
(defn shift-to
  "Displaces all atoms in mol such that atom n is at point P."
  [P n mol]
  (let [shift-vec (->> mol
                    (mol-filter {:pos n})
                    (first)
                    (:coordinates)
                    (- P))]
    (shift shift-vec  mol)))




(defn mol-center
  "Returns a molecule whose center of the atm positions is (0, 0, 0)"
  [mol]
  (shift (- (coord-average mol)) mol))



(defn rotate-mol
  "pnt2 and pnt1 are both 3-tuples that define the axis.
angle is the amount of ration this is a scalar value."
   [mol pt1 pt2 angle]
  (let [f #(jmath/the-rotation-function % pt1 pt2 angle)]
           (update-mol :coordinates f mol)))


(defn random-rotate-mol
  "Usage: (random-rotate-mol DMMP) "
  ([mol]
         (let [pnt1 (jmath/random-point [0 0 0] [1 1 1])
               pnt2 (jmath/random-point  [2 2 2] [4 4 4])]
           (random-rotate-mol mol pnt1 pnt2)))
  ([mol pnt1 pnt2]
    (let [angle (* (rand 2) Math/PI)]
      (rotate-mol mol pnt1 pnt2 angle)))
  ([mol pnt1 pnt2 angle]
    (let [pnt3 (first
                (drop-while
                       (fn [x](> (jmath/find-angle (- pnt2 pnt1) (map - x pnt1))
                                angle))
                  (repeatedly #(jmath/random-point [2 2 2] [4 4 4]))))]
      (random-rotate-mol pnt1 pnt3))))


(defn apply-coord-transform-matrix
  "This will apply the coordiante transformation matrix, mat, to the coordinates
of the atoms in mol, leaving the mol in the transformed frame.  The transformation
could be a rotation, translation, mirror, etc."
  ([mat mol]
  (update-mol :coordinates #(cmatrix/mmul mat %) mol))
  ([mat]
  (update-mol :coordinates #(cmatrix/mmul mat %))))


(defn replace-mol-coordinates
  "This is designed to replace all of the coordinates in a mol with coordinates
from another mol. In order for this to make sense the atoms have to be in the
same order.

I expect this to be useful when using clojure to run fireball/gulp/etc and you
want to read the output coordinates back in, but leave the other fields in the
struct unchanged."
  [replace-into-mol replace-from-mol]
  (map #(assoc-in %1 [:coordinates] (:coordinates %2)) replace-into-mol replace-from-mol))


(defn coords->zero
  "Many times in the creation of mols we end up with positions that have a
coordinate that should be zero, but is represented by something stupid
like: 3.3306690738754696E-16.  I hate that.  This function gets rid of those
nonzero zeros.  This is also of importance when this is used to make input for
legacy code."
  ([mol]
  (update-mol :coordinates #(map jmath/setzero %) mol))
  ([]
  (update-mol :coordinates #(map jmath/setzero %))))





(defn- center-of-position-
  [mol]
  (map #(/ % (count mol)) (reduce #(map + %1 %2) (map :coordinates mol))))





(defn- pivot-axis-
  "I came across a use case where my original rotate-sub-mol needed some help,
if the center-of-position of the submol is along the line connecting the pivot
and pntb; in that case it still worked, but not in the way that I wanted it.  I
am adding this helper function to try to make it work the way that I want."
  [cop pivot pntb]
  (cond
    (= (first pivot) (first pntb))
    (cmatrix/cross (- pivot pntb) (+ pivot [1 0 0]))
    (= (second pivot) (second pntb))
    (cmatrix/cross (- pivot pntb) (+ pivot [0 1 0]))
    (= (last pivot) (last pntb))
    (do
    (cmatrix/cross (- pivot pntb) (+ pivot [0 0 1]))
      (println [cop pivot pntb])
     (println (cmatrix/cross (- pivot pntb) (+ pivot [0 0 1]))) )
    :else
    (cmatrix/cross (- pntb pivot) (- cop pivot))))




(defn rotate-sub-mol
  "This is used to rotate part of a mol with respect to the rest of the mol.
This works a lot like how GuassView works.  After defining the mol and submol you
specify three atoms that will be involved.  The first atom is atom in an atom in
the submol, this atom will be moving.  You then specify an atom that will be the pivot.
The last atom is used to define the direction of the movement.  A positive angle,
specifies that it should rotate away from the third atom, a negetive angle means
the sub mol will rotate closer."
  [mol submol pivot pntb angle]
  (let [mmol (cset/difference (set mol) (set submol))
        cop (center-of-position- submol)
        V (if (:coordinates pivot)
            (:coordinates pivot)
            pivot)
        P (if (:coordinates pntb)
            (:coordinates pntb)
            pntb)
        X (if (zero? (jmath/setzero (cmatrix/dot (- P V) (- cop V))))
            (do (println "chad")(pivot-axis- cop V P))
            (do (println "junkermeier")(cmatrix/cross (- P V) (- cop V))))]
    (concat mmol (rotate-mol submol V (+ V X) (- angle)))))



(defn mol-find-dihedral
  [mol atoms]
  (let [[a b c d] atoms
        f #(- (:coordinates (mol-nth mol %2)) (:coordinates (mol-nth mol %1)))]
    (jmath/dihedral (f b a) (f b c) (f c d))))



(defn mol-angle
  [mol atms]
  (let [[a P c] atms]
        (jmath/three-point-angle (:coordinates (mol-nth mol a))
                    (:coordinates (mol-nth mol P))
                    (:coordinates (mol-nth mol c)))))


(defn mol-vector
  "This is a little helper function that I wrote to find the vector between
  two atoms.  I wrote this because I am tired of always having to write
  the code over and over again.  It finds the vector that points from atm1 to atm2."
  [mol atm2 atm1]
  (- (:coordinates (mol-nth mol atm2))
         (:coordinates (mol-nth mol atm1))))


(defn arrange-mol
   "This is used to rotate a mol such that the vector defined by (map - atm2 atm1)
   is parallel to the vector direction. It also places atm1 at the origin.

   Usage: (arrange-mol h2o-ts 0 1 [0 0 1])"
   [mol atm1 atm2 direction]
   (shift-to
    (apply-coord-transform-matrix mol
     (jmath/align-vectors (mol-vector mol atm2 atm1) direction)) atm1 [0 0 0]))




(defn- NEB-helper-
  "Adds another element to the value of key with the elements of col.
There needs to be the same number of elements in col and mol."
[mol  col]
  (let [coll (map #(+ %1 %2) (map :coordinates mol) col)]
    (map #(assoc-in %1 [:coordinates] %2) mol coll)))



(defn NEB-guess
  "Produces intermediate images for an Nudged Elastic band calculation."
  [mol1 mol2 n]
  (let [inv (/ 1 n)
        avec (map #(- (:coordinates %1)(:coordinates %2)) mol1 mol2)
        bvec (fn [x] (map #(map (partial * x inv) %) avec))]
    (pmap #(NEB-helper- mol2 (bvec %)) (reverse (take (inc n) (iterate inc 0))))))



(defn resize-bond
  "This moves atm2 such that the bond length between atm1 and atm2 is lngth."
  [mol atm1 atm2 lngth]
  (let [a (map (partial * lngth) (cmatrix/normalise (mol-vector mol atm1 atm2)))
        b (+ (:coordinates (mol-nth mol atm2)) a)]
    (find-assoc-in [:pos :coordinates] [atm1 b] mol)))




(defn center-of-mass
  "Returns the center of mass of mol."
  [mol]
  (/
     (->> mol
          (map #(* ((comp ed/atomic-mass ed/atomic-numbers :species) %) (:coordinates %)))
          (reduce + [0.0 0.0 0.0]))
    (->> mol
         (map :species)
         (map ed/atomic-numbers)
         (map ed/atomic-mass)
         (reduce + 0.0))))




(defn update-mol-name
  "This will add a value (be it a number, string, or keyword) to :name.
The assumption is that you will only want to add values to :name based on some
predicate, but that you might want to erase a value or values based on the
circumstances.  Thus there are three ways to delete the values.

The value associated with the key name will be a set.  Thus, the same
atom in mol can have several names attached to it.  To add a name value you do
the following:
Usage: (update-mol-name some-mol some-pred some-name)

If you want to delete all of the name values of all of the atoms, you do the this:
Usage: (update-mol-name some-mol)

If you want to delete a particular name from all of the atoms (it doesn't matter
                                                                if a particular atom actually has that name), you should do the following:
Usage (update-mol-name some-mol some-name)

Finally, if you want to remove a name based on some predicate you do this:
Usage (update-mol-name some-mol some-pred some-name :remove)"
  ([mol]
    (map #(assoc-in % [:name] nil) mol))
  ([mol removed-name]
    (map #(assoc-in % [:name] (cset/difference (:name %) (if (set? removed-name)
                                                      removed-name
                                                      #{removed-name}))) mol))
  ([mol pred added-name]
    (map #(if (pred %) (if (:name %)
                         (assoc-in % [:name] (flatten-set (conj (:name %) added-name)))
                         (assoc-in % [:name] (flatten-set #{added-name})))
            (assoc-in % [:name] (:name %))) mol))
  ([mol pred removed-name k]
            (map #(if (pred %)
                    (assoc-in % [:name] (cset/difference (:name %) (if (set? removed-name)
                                                                removed-name
                                                                #{removed-name})))
                    (assoc-in % [:name] (:name %))) mol)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





;(use 'greenwood.xyz)
;(require '[clojure.core.reducers :as r])

#_(def g (second (->> (foldable-chunks "/Users/chadjunkermeier/Desktop/graphene.xyz" [[0 77] [77 82]])
(r/map (partial drop 2))
     (r/map xyz-iota->atoms)
     (into []))))


#_(->> g
(shift [0 0 10])
)









