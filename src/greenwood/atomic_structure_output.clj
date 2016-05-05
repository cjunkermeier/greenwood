(ns greenwood.atomic-structure-output
  (:require
   [clojure.math.combinatorics :as mathcomb]
   [clojure.string :as strng]
   [greenwood.empirical-data :as empdata]
   [greenwood.mol :as jmol]
   [greenwood.math :as jmath]
   [clojure.pprint :as pp]
   [greenwood.utils :as utils]))



(defn ^:private ^:static
  reduce1
  ([f coll]
    (let [s (seq coll)]
      (if s
        (reduce1 f (first s) (next s))
        (f))))
  ([f val coll]
    (let [s (seq coll)]
      (if s
        (if (chunked-seq? s)
          (recur f
            (.reduce (chunk-first s) f val)
            (chunk-next s))
          (recur f (f val (first s)) (next s)))
        val))))

(comment "This is in no way a complete specification of writing output in either
the PDB or CrystalMaker format, it is only the parts that I needed in a
particular project.")



(defn commutative-cartesian-product
  "This produces a set of elements, similar to the cartesian product on two
instances of the col, where the ordering of the items in the elements don't
matter, ie (1 2) = (2 1).
Usage: (commutative-cartesian-product [1 2]) => ((1 1) (2 1) (2 2))"
  [col]
  (letfn [(filterer [unfiltered filtered]
            (if (empty? unfiltered)
              filtered
              (let [sorted-equals #(not (= (sort (first unfiltered)) (sort %)))]
                (recur (filter sorted-equals unfiltered) (cons (first unfiltered) filtered)))))]
    (filterer (mathcomb/cartesian-product col col) [])))



(defn get-atoms [atoms convert-to-nZ]
  "get atoms is used to turn a mol into a col of cols where the subcols are nothing
more than vectors of the species and the coordinates values.  In this function
we also test to see if the mol has values for the :pos key, if all of the atoms
do then we sort the mol by :pos."
  (let [sorted-atoms (if (not-any? (comp nil? :pos) atoms)
                       (sort-by :pos atoms)
                       atoms)
    f1 (comp #(if convert-to-nZ (empdata/atomic-numbers %) %) :species)
	f2 (comp #(concat (map float (take 3 %)) (drop 3 %)) :coordinates)]
    (map #(cons (f1 %) (f2 %)) sorted-atoms)))



(defn- check-charge-
  [mol]
  (if ((comp nil? :charges first) mol)
    (let [one (repeat 1)]
      (jmol/col->mol mol :charges one))
    mol))


(defn species-charge-
  "This will do it's best to find the correct Z value for :species.
For example, if {:species \"C\"} the Z value should be 6 for C^{12}."
  [species]
  (if (integer? species)
    species
    (empdata/atomic-numbers species)))



(defn write-cmtx-ATOM
  "The listing of atoms must be in the asymmetric unit cell."
  [mol]
  (let [species (set (map #(:species %) mol))
        blanked-names (jmol/update-mol-name mol)
        properly-named (sort-by :pos (mapcat (fn [x] (map #(assoc-in %1 [:name] (str x %2))
                                                       (jmol/mol-filter :species x blanked-names) (iterate inc 1))) species))
        f #(strng/join "     " [(:species %) (:name %) (first (:coordinates %)) (second (:coordinates %)) (last (:coordinates %))])]
    (str "ATOM" utils/endline (strng/join utils/endline (map f properly-named)) utils/endline)))



(defn write-cmtx-lat-params
  "Listing of lattice parameters and space group in Crystal Maker Format."
  [params space-group]
  (str "CELL "
    (strng/join "     " (take 3 params)) "     "
    (strng/join "     " (map #(* % (/ 180 Math/PI)) (drop 3 params)))
    utils/endline "SPGR " space-group utils/endline))



(defn write-cmtx-BOND
  "The cmtx data file expects a listing of the atomic species that may have bonds
between them.  This function is going to do the dumbest thing possible and just
write out all combinations of atom species."
  [mol]
  (let [bonding-pairs (commutative-cartesian-product (set (map #(:species %) mol)))]
    (str (strng/join utils/endline (map #(str "BOND     " (strng/join "     " %)) bonding-pairs)) utils/endline)))



(defn write-cmtx
  "This is the one to call."
  [mol lat-params space-group]
  (str (write-cmtx-lat-params lat-params space-group)
    (write-cmtx-BOND mol)
    "XYZR  -1.5 1.5 -1.5 1.5 -1.5 1.5" utils/endline
    (write-cmtx-ATOM mol)))



(defn write-bas [atoms]
  "Returns atoms as a string in fireball's bas format."
  (str (count atoms)  utils/endline (utils/inter-cat-tree [utils/endline "   "] (get-atoms atoms true)) utils/endline))



(comment "ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf")

(defn write-pdb-lat
  "The Hermann-Mauguin space group symbol is given without parenthesis, e.g.,
P 43 21 2. Please note that the screw axis is described as a two digit number.
The full International Table’s Hermann-Mauguin symbol is used, e.g., P 1 21 1
instead of P 21.  For a rhombohedral space group in the hexagonal setting, the
lattice type symbol used is H."
  [params space-group]
  (let [[a b c al be ga] params
        alpha (* al (/ 180 Math/PI))
        beta (* be (/ 180 Math/PI))
        gamma (* ga (/ 180 Math/PI))]
    (str (pp/cl-format nil "CRYST1~9,3F~9,3F~9,3F~7,2F~7,2F~7,2F~11:<~A~>"
           a b c alpha beta gamma space-group) utils/endline)))



(defn write-pdb-atom-name
    ""
    [atom-from-mol]
    (let [names (filter string? (:name atom-from-mol))]
      (if (empty? names)
        "nil"
        (apply str (take 5 (first names))))))



(defn write-bfactor-
       ""
       [atom-from-mol]
       (if (nil? (:charge atom-from-mol))
         0
        ((comp :total :charge) atom-from-mol)))





(defn write-pdb-HETATM [mol]
  (let [mmol (check-charge- mol)
        f #(pp/cl-format nil "HETATM~5D ~4@<~A~> ~3:<~A~>  ~4D    ~8,3F~8,3F~8,3F~6,2F~6,3F          ~2:<~A~>"
             (inc (:pos %))
             (:species %)
             "nil"
             1
             (first (:coordinates %))
             (second (:coordinates %))
             (last (:coordinates %))
             1.0
             (write-bfactor- %)
             (species-charge- (:species %)))]
    (str (strng/join utils/endline (map f mmol)) utils/endline)))



(defn write-pdb-HETATM [mol]
  (let [mmol (check-charge- mol)
        f #(pp/cl-format nil "HETATM~5D ~4@<~A~> ~3:<~A~>  ~4D    ~8,3F~8,3F~8,3F~6,2F~6,3F          ~2:<~A~>"
             (inc (:pos %))
             (:species %)
             "nil"
             1
             (first (:coordinates %))
             (second (:coordinates %))
             (last (:coordinates %))
             1.0
             (write-bfactor- %)
             (:species %))]
    (str (strng/join utils/endline (map f mmol)) utils/endline)))





(defn write-pdb-HETATM [mol]
  (let [mmol (check-charge- mol)
        f #(pp/cl-format nil "HETATM~5D ~4@<~A~> ~3:<~A~>  ~4D    ~8,3F~8,3F~8,3F~6,2F~6,2F          ~2:<~A~>"
             (inc (:pos %))
             (:species %)
             "nil"
             1
             (first (:coordinates %))
             (second (:coordinates %))
             (last (:coordinates %))
             1.0
             (write-bfactor- %)
             (:species %))]
    (str (strng/join utils/endline (map f mol)) utils/endline)))


(defn get-npos
  "This is a helper function that will create a seq of the :npos of the atoms that are a neighbor to atomm"
  [atomm]
  (loop [n (first (:neigh atomm))
         nn (rest (:neigh atomm))
         m []]
    (if (nil? n)
      m
      (recur (first nn) (rest nn) (conj m (:npos n))))))




(defn write-pdb-connect
  [mol]
  (let [un-named (fn [x y]
                   (cond
                     (= 4 (count y))
                     (let [[a b c d] (map inc y)]
                       (pp/cl-format nil "CONECT~5D~5D~5D~5D~5D" (inc x) a b c d))
                     (= 3 (count y))
                     (let [[a b c] (map inc y)]
                       (pp/cl-format nil "CONECT~5D~5D~5D~5D" (inc x) a b c))
                     (= 2 (count y))
                     (let [[a b] (map inc y)]
                       (pp/cl-format nil "CONECT~5D~5D~5D" (inc x) a b))
                     (= 1 (count y))
                     (let [[a] (map inc y)]
                       (pp/cl-format nil "CONECT~5D~5D" (inc x) a))
                     (= 0 (count y))
                     (pp/cl-format nil "CONECT~5D" (inc x))))]
    (str (strng/join utils/endline (map #(un-named %1 %2) (iterate inc 0) (map get-npos mol))) utils/endline)))



(defn write-pdb
  ([mol lat-param space-group]
   (if (nil? ((comp :neigh first) mol))
                 (str "HEADER" utils/endline
                  "AUTHOR    GENERATED IN JMD" utils/endline
                 (write-pdb-lat lat-param space-group)
                 (write-pdb-HETATM mol)
                 "END")
    (str "HEADER" utils/endline
      "AUTHOR    GENERATED IN JMD" utils/endline
      (write-pdb-lat lat-param space-group)
      (write-pdb-HETATM mol)
      (write-pdb-connect mol)
      "END")))
  ([mol lat-vec]
    (if (nil? ((comp :neigh first) mol))
                 (str "HEADER" utils/endline
                  "AUTHOR    GENERATED IN JMD" utils/endline
      (write-pdb-lat (jmath/lvs->parameters lat-vec) "P 1")
                 (write-pdb-HETATM mol)
                 "END")
    (str "HEADER" utils/endline
      "AUTHOR    GENERATED IN JMD" utils/endline
      (write-pdb-lat (jmath/lvs->parameters lat-vec) "P 1")
      (write-pdb-HETATM mol)
      (write-pdb-connect mol)
      "END"))))




(defn write-xyz [timesteps]
  "Returns a seq of atoms as a string in xyz format.  This version of write-xyz
allows for writing a whole bunch of time steps (arranged in a col of cols of maps)
or a single time step which is a col of maps.

Usage:  Suppose (def test (xyz-str->atoms 'C 0 0 0 \n C 0.3333 0.6667 0')) then
the following would both give the same result (write-xyz test) => '2\n\n C 0 0 0 \n C 0.3333 0.6667 0'
and (write-xyz (vector test)) => '2\n\n C 0 0 0 \n C 0.3333 0.6667 0'."
  (if (:species (first timesteps))
    (str (count timesteps)  utils/endline utils/endline (utils/inter-cat-tree [utils/endline "   "] (get-atoms timesteps false)) utils/endline)
    (apply str (map #(str (count %) utils/endline utils/endline (utils/inter-cat-tree [utils/endline "   "] (get-atoms % false)) utils/endline) timesteps))))



(defn write-just-xyz [timesteps]
  "Returns a seq of atoms as a string in xyz format, minus the number of atoms
being tacked to the front.

Usage:  Suppose (def test (xyz-str->atoms 'C 0 0 0 \n C 0.3333 0.6667 0')"
  (str (utils/inter-cat-tree [utils/endline "   "] (get-atoms timesteps false))))



(defn write-xyz-limitdec
  "Returns a seq of atoms as a string in xyz format.  It works just like write-xyz
                  but with the caviot that it limits the number of decimal points that it will
                  print.  The user can also supply the decimal to which it should be rounded."
  ([timesteps]
    (write-xyz-limitdec timesteps 0.000001))
  ([timesteps  decimal]
    (if (:species (first timesteps))
      (write-xyz (jmol/update-mol timesteps :coordinates (fn [y] (map (fn [x](jmath/round-decimal x 0.000001)) y))))
      (write-xyz
        (map #(jmol/update-mol % :coordinates (fn [y] (map (fn [x](jmath/round-decimal x 0.000001)) y))) timesteps)))))



(defn keys->xyz
  "At times we want to be able to add other information onto the end of each line
                  to make a new, if temporary, file format.  For example you
                  might want a line that looks like:
                  species  xpos ypos zpos totalcharge.
                  Once you run this then apply write-xyz to it as you would any
                  other mol.
                  Usage: (keys->xyz g (comp :total :charge))
                  Usage (doesn't work yet): (keys->xyz g (comp :total :charge) :pos) "
  ([mol k]
    (map #(assoc-in % [:coordinates] (flatten (conj (:coordinates %) (k %)))) mol)))



(defn write-xyz-cshell
  ([timestep]
    "Returns a seq of atoms as a string in xyz format.  THIS ONLY WORKS ON ONE
                  TIMESTEP.  On some machines there is an abnormally small maximum number of
                  characters that can be in a word, if this is the case, give it a maximum number of lines."
    (str (utils/inter-cat-tree ["\\\n" "\t"] (get-atoms timestep false))))
  ([timestep max-lines name-string]
    (strng/join "\n\n" (map #(str "set " name-string  "=\""
                         (utils/inter-cat-tree ["\\\n" "\t"] (get-atoms %2 false)) "\"")
                   (iterate inc 1) (partition-all max-lines timestep)))))


(defn write-lvs-cshell
  [lvs name]
  (str "set LVS" name " = \""
       (utils/inter-cat-tree ["\\\n" "\t"] lvs) "\"" ))



(defn write-xyz-qe-NEB [timesteps]
  ""
  (let [num ((comp dec dec count) timesteps)
        images (flatten ["FIRST_IMAGE" (take num (repeat "INTERMEDIATE_IMAGE")) "LAST_IMAGE"])]
      (apply str (map #(str %2 utils/endline "ATOMIC_POSITIONS angstrom" utils/endline
                            (utils/inter-cat-tree [utils/endline "   "] (get-atoms %1 false)) utils/endline) timesteps images))))



(defn write-dftbplus-xyz
  ([mol]
  (let [s (set (map :species mol))
        lc (apply merge (map #(hash-map %1 %2) s (iterate inc 1)))] ;lc is used to do list comprehension.
    (strng/join utils/endline
     (concat [(str (count mol) " C")
      (strng/join " " s)]
             (map #(str %1 " " (lc (:species %2)) " " (strng/join " " (:coordinates %2))) (iterate inc 1) mol)))))
        ([mol lvs]
  (let [s (set (map :species mol))
        lc (apply merge (map #(hash-map %1 %2) s (iterate inc 1)))] ;lc is used to do list comprehension.
    (str
    (strng/join utils/endline
     (concat [(str (count mol) " S")
      (strng/join " " s)]
             (map #(str %1 " " (lc (:species %2)) " " (strng/join " " (:coordinates %2))) (iterate inc 1) mol)))
    utils/endline
     "0.00 0.00 0.00"
     utils/endline
     (utils/inter-cat-tree [utils/endline " "] lvs)))))




(defn write-dftbplus-xyz-ricardo
    [mol lvs]
  (let [s (set (map :species mol))
        lc (apply merge (map #(hash-map %1 %2) s (iterate inc 1)))] ;lc is used to do list comprehension.
    (strng/join utils/endline
[" Geometry = {"
               (str "TypeNames = { "\" (strng/join "\" \"" s) "\" }")
               "TypesAndCoordinates [Angstrom] = {"
             (strng/join utils/endline (map #(str (lc (:species %)) " " (strng/join " " (:coordinates %))) mol))
           " }
 Periodic = Yes
 LatticeVectors [Angstrom] = {"
             (utils/inter-cat-tree [utils/endline " "] lvs)
"  }
}"])))




(defn lammps-datain
  [mol lvs]
  (let [species (set (map :species mol))
        lc (apply merge (map #(hash-map %1 %2) species (iterate inc 1)))]
          (str
            "Input created by CEJ\n\n"
            (str (count mol) " atoms\n")
            (str (count species) " atom types\n")
            utils/endline
            (str 0 " " (ffirst lvs) " xlo xhi\n")
            (str 0 " " ((comp second second) lvs) " ylo yhi\n")
            (str 0 " " ((comp last last) lvs) " zlo zhi\n")
            utils/endline
            "Masses\n\n"
            (strng/join utils/endline (map #(str (lc %) " " ((comp empdata/atomic-mass empdata/atomic-numbers) %) "   #" %) species))
            utils/endline utils/endline
            "Atoms\n"
            utils/endline
               (strng/join utils/endline
           (map #(str %1 " "
                      ((comp lc :species) %2) " 0 "
                      (strng/join " " (:coordinates %2))) (iterate inc 1) mol)))))



(defn lammps-NEB-final
  [mol]
  (strng/join utils/endline (map #(str %1 " " (strng/join " " (:coordinates %2))) (iterate inc 1) mol)))



(defn write-jaguar-xyz
  [mol]
  (let [species (map #(str (:species %1) %2) mol (iterate inc 1))]
    (write-just-xyz (jmol/col->mol mol :species species))))



(defn write-lammps-trajectory
  [timesteps]
    (let [species (set (map :species (first timesteps)))
          lc (apply merge (map #(hash-map %1 %2) species (iterate inc 1)))]
      (letfn [(tm [timestep natoms]
              (strng/join utils/endline
                    ["ITEM: TIMESTEP"
                     timestep
                     "ITEM: NUMBER OF ATOMS"
                     natoms
                     "ITEM: BOX BOUNDS pp pp pp
-1000 1000
-1000 1000
-1000 1000
ITEM: ATOMS id type x y z vx v_trstress"]))
          (m [mol]
             (strng/join utils/endline (map #(str %1 " " ((comp lc :species) %2) " "(strng/join " " (:coordinates %2)) " " (:charge %2) " 0.0")
                  (iterate inc 1) mol)))]
    (strng/join utils/endline (map #(str (tm %1 (count %2)) utils/endline (m %2)) (iterate inc 1) timesteps)))))





(defn write-xyz-string [timesteps strings]
  "Returns a seq of atoms as a string in xyz format.  This version of write-xyz
allows for writing a whole bunch of time steps (arranged in a col of cols of maps)
or a single time step which is a col of maps.

Usage:  Suppose (def test (xyz-str->atoms 'C 0 0 0 \n C 0.3333 0.6667 0')) then
the following would both give the same result (write-xyz-string test) => '2\n\n C 0 0 0 \n C 0.3333 0.6667 0'
and (write-xyz (vector test)) => '2\n\n C 0 0 0 \n C 0.3333 0.6667 0'."
  (if (:species (first timesteps))
    (str (count timesteps)  utils/endline strings utils/endline (utils/inter-cat-tree [utils/endline "   "] (get-atoms timesteps false)) utils/endline)
    (apply str (map #(str (count %1) utils/endline %2 utils/endline (utils/inter-cat-tree [utils/endline "   "] (get-atoms %1 false)) utils/endline) timesteps strings))))





(defn write-pureMD-lat
  "takes params and outputs a string that is used in the PuReMD input format"
  [params]
  (let [[a b c alpha beta gamma] params]
    (str (pp/cl-format nil "BOXGEO~9,3F~9,3F~9,3F~7,2F~7,2F~7,2F"
           a b c alpha beta gamma) utils/endline)))


(defn pureMD-atoms [mol]
  "get atoms is used to turn a mol into a col of cols where the subcols are nothing
more than vectors of the species and the coordinates values.  In this function
we also test to see if the mol has values for the :pos key, if all of the atoms
do then we sort the mol by :pos."
  (let [sorted-atoms (if (not-any? (comp nil? :pos) mol)
                       (sort-by :pos mol)
                       mol)
    f1 (comp #(if false (empdata/atomic-numbers %) %) :species)
	f2 (comp #(if (nil? %) "atm" %) :name)
  f3 (comp #(strng/join "  " %) :coordinates)]
    (strng/join "\n" (map #(str %1 "  " (f1 %2) "  " (f2 %2) "  " (f3 %2) ) (iterate inc 1)  sorted-atoms))))



(defn pureMD-input
  ""
  ([mol lvs]
    (str
     (write-pureMD-lat (jmath/lvs->parameters lvs))
    (count mol)
     "\n"
     (pureMD-atoms mol)))
  ([mol a b c alpha beta gamma]
    (str
     (write-pureMD-lat [a b c alpha beta gamma])
    (count mol)
     "\n"
     (pureMD-atoms mol))))











