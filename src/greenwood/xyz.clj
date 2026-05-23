(ns greenwood.xyz
  "Utilities for determining chunk position in XYZ files and processing
  chunks using reducers.

  Uses iota, which uses mmap()."
  ;(:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.reducers :as r]
            [clojure.string :as strng]
            [greenwood.basics :as basic]
            [greenwood.utils :as utils]
            [greenwood.mol :as gmol]
            [greenwood.math :as gmath]
            [clojure.core.matrix :as cmat]
            iota
            [clojure.core.matrix.operators :as cmato]))


(defn atom-pos
  "This will associate file positions to the :pos of each atom (starting with zero).
Usage: (atoms-pos mol)
Usage: (atoms-pos mol  5)"
    ([mol]
    (map (fn [x y] (assoc-in x [:pos] y)) mol (iterate inc 0)))
    ([mol start]
     (map (fn [x y] (assoc-in x [:pos] y)) mol (iterate inc start))))




(defn- step-chunk-starts
  "chad"
  ([] [(vector-of :int) 0 :find-start])
  ([[v i state]] (conj v i))
  ([[v i state] ^String l]
    (case state
      :find-start (if (re-matches #"[ \t]*\d+[ \t]*" l)
                    [(conj v i) (inc i) :comment]
                    [v (inc i) :find-start])
      :comment [v (inc i) :find-start])))


(defn chunk-starts
  "Returns `[chunk-start-line ... total-num-of-lines]` from a reducible
  `lines` of XYZ file lines."
  [lines]
  (->> lines
       (r/reduce step-chunk-starts)
       (step-chunk-starts)))


(defn chunk-ranges
  "Returns `[[start-idx0 end-idx0] [end-idx0 start-idx1] ...]`.
  Same as `(partition 2 1 start+total)` but non-lazy.

  Designed to produce arguments for subvec to grab groups of items at a time."
  [start+total]
  (->> (partition 2 1 start+total)
       (reduce (fn [ranges [s e]] (conj! ranges (vector-of :int s e)))
         (transient []))
       (persistent!)))


(defn index-xyz*
  "Return a vec of vecs, each of which is a range of indexes in coll which
  constitutes a single XYZ chunk. Retrieve the chunks with
  `(map #(apply subvec coll %) (index-xyz* coll))`.
  This "
  [coll]
  (->> coll
       chunk-starts
       chunk-ranges))


(defn index-xyz
  "Returns the chunk index of an XYZ file.  This can be used to cache the
  index or even write it to disk for iota to read later.  This would be useful
  for very large files.  It would allow for copying the xyz file and index file to
  several machines and processing the file in parallel."
  [xyzfile]
  (index-xyz* (iota/seq xyzfile)))



(defn natoms-index
  "Returns the chunk index of an XYZ file assuming that each time step has
  the same number of atoms.  The user will also need to determine the number
  of lines in the xyz file before running.  This will greatly speed up the
  reading of the the xyz-file.
  Usage: (natoms-index 3 10) => ((0 5) (5 10))"
  [natoms nlines]
  (partition 2 1 (range 0 (inc nlines) (cmato/+ 2 natoms))))





(defn reax-index-timesteps
  [start stop howoften natoms]
 (partition 2 1 (map int (range (* (/ start howoften) (+ 2 natoms))
  (* ((comp inc inc)  (/ stop howoften)) (+ 2 natoms)) (+ 2 natoms)))))





(defn foldable-chunks*
  "Return a foldable collection of the chunks in coll."
  ([coll]
    (foldable-chunks* coll (index-xyz* coll)))
  ([coll index]
    (r/map (fn [[s e]] (subvec coll s e)) index)))



(defn foldable-chunks
  "Return a foldable collection chunks in an XYZ file.
  index-xyz can be used to precompute an index of the
  chunk ranges.  Which can then be stored for later
  consumation of the xyzfile."
  ([xyzfile]
    (foldable-chunks* (iota/vec xyzfile)))
  ([xyzfile index]
    (foldable-chunks* (iota/vec xyzfile) index)))










(defn take-nth-foldable-chunks
  "Return a foldable collection chunks in an XYZ file.
  index-xyz can be used to precompute an index of the
  chunk ranges.  Which can then be stored for later
  consumation of the xyzfile.

  Returns a lazy seq of every nth item in index."
  ([xyzfile n]
    (foldable-chunks* (iota/vec xyzfile) (take-nth n (index-xyz xyzfile))))
  ([xyzfile index n]
    (foldable-chunks* (iota/vec xyzfile) (take-nth n index))))


(defn take-foldable-chunks
  "Return a foldable collection chunks in an XYZ file.
  index-xyz can be used to precompute an index of the
  chunk ranges.  Which can then be stored for later
  consumation of the xyzfile.

  Returns a lazy seq of every nth item in index."
  ([xyzfile n]
    (foldable-chunks* (iota/vec xyzfile) (take n (index-xyz xyzfile))))
  ([xyzfile index n]
    (foldable-chunks* (iota/vec xyzfile) (take n index))))








(comment
  "Silly example: count all the chunks."
  (->> (foldable-chunks "myfile.xyz")
       (r/map (constantly 1))
       (r/fold cmato/+))
  "Same."
  (->> (foldable-chunks "myfile.xyz")
       (r/fold cmato/+ (fn ([] 0) ([x _] (inc x)))))
  "Get comment line of each chunk."
  (->> (foldable-chunks "myfile.xyz")
       (r/map (fn [atom-count comment & atoms] comment))
       (r/foldcat)))




(defn xyz-iota->atoms
  "This will parse a string into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.

Thus if: (def test 'C 0 0 0 \n C 0.3333 0.6667 0')
then the usage would be (xyz-str->atoms test)."
 ([lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix :double-array  (map read-string (take 3 (rest %)))) nil nil nil nil y)
                (strng/split (strng/triml x) #"\s+")))
 lines (iterate inc 0)))
 ([charge-column lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix :double-array  (map read-string (take 3 (rest %))))  (double (read-string (nth % charge-column))) nil nil nil y)
                (strng/split (strng/triml x) #"\s+")))
     lines (iterate inc 0))))


(defn xyz-iota->atoms-readable
  "This will parse a string into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.

Thus if: (def test 'C 0 0 0 \n C 0.3333 0.6667 0')
then the usage would be (xyz-str->atoms test)."
  ([lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix (map read-string (take 3 (rest %)))) nil nil nil nil y)
                (strng/split (strng/triml x) #"\s+")))
 lines (iterate inc 0)))
 ([charge-column lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix (map read-string (take 3 (rest %))))  (double (read-string (nth % charge-column))) nil nil nil y)
                (strng/split (strng/triml x) #"\s+")))
     lines (iterate inc 0))))



(defn xyz-str->atoms
  "This will parse a string into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.

Thus if: (def test 'C 0 0 0 \n C 0.3333 0.6667 0')
then the usage would be (xyz-str->atoms test)."
  ([string]
  (let [lines (strng/split-lines string)]
    (xyz-iota->atoms lines)))
  ([charge-column string]
    (let [lines (strng/split-lines string)]
      (xyz-iota->atoms charge-column lines))))


(defn xyz-str->atoms-readable
  "This will parse a string into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.

Thus if: (def test 'C 0 0 0 \n C 0.3333 0.6667 0')
then the usage would be (xyz-str->atoms-readable test)."
  ([string]
  (let [lines (strng/split-lines string)]
    (xyz-iota->atoms-readable lines)))
  ([charge-column string]
  (let [lines (strng/split-lines string)]
     (xyz-iota->atoms-readable charge-column lines))))








(defn xyz-reax-iota->atoms
  "This will parse a set of strings into the atoms struct.
The assumption is that the set of lines are from a reaxff xmolout file, where ixmolo was
set such that the molecule number was printed out in the 5th column of every line.
Note that the string should start with the first atom, not with the number of atoms
in the system.  Also, this assumes that there is a newline character between atoms."
  ([lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix :double-array  (map read-string (take 3 (rest %))))  nil nil nil (read-string (nth % 4)) y)
                (strng/split (strng/triml x) #"\s+")))
     lines (iterate inc 0)))
  ([charge-column lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix :double-array  (map read-string (take 3 (rest %))))  (read-string (nth % charge-column)) nil nil (read-string (nth % 4)) y)
                (strng/split (strng/triml x) #"\s+")))
     lines (iterate inc 0))))


(defn xyz-reax-iota->atoms-readable
  "This will parse a set of strings into the atoms struct.
The assumption is that the set of lines are from a reaxff xmolout file, where ixmolo was
set such that the molecule number was printed out in the 5th column of every line.
Note that the string should start with the first atom, not with the number of atoms
in the system.  Also, this assumes that there is a newline character between atoms."
  ([lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix  (map read-string (take 3 (rest %))))  nil nil nil (read-string (nth % 4)) y)
                (strng/split (strng/triml x) #"\s+")))
     lines (iterate inc 0)))
  ([charge-column lines]
    (mapv (fn [x y] (#(basic/new-atom (.intern (first %)) (cmat/matrix  (map read-string (take 3 (rest %))))  (read-string (nth % charge-column)) nil nil (read-string (nth % 4)) y)
                (strng/split (strng/triml x) #"\s+")))
     lines (iterate inc 0))))




(defn parse-xyz
  "This reads the whole xyz-file into memory, thus this should be used only for small files.
  Since there are a number of files types that are based on the xyz-file system where they
  include additional data (ie. charge, velocity, strain) we are allowing for charge to also
  be read in; but in order to do that you will have to specify which column the charge is
  placed in.

This produces a col of cols, where each of the sub-cols is a time step.

Usage: (second (parse-xyz PATH)), where PATH is a string containing the path
to some xyz file.

  Usage: (second (parse-xyz PATH 5)), where PATH is a string containing the path
to some xyz file, and the charge is given in the 5 column.  The column number starts
  counting from one."
  ([filename]
  (->> (foldable-chunks filename)
       (r/map (partial drop 2))
       (r/map xyz-iota->atoms-readable)
       (into [])))
  ([filename charge-column]
  (->> (foldable-chunks filename)
       (r/map (partial drop 2))
       (r/map (partial xyz-iota->atoms-readable charge-column))
       (into []))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  xmolout
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- reaxff-cell-params-lvs-
  "Cell parameters should be banned!  But since they are not, and Reaxff uses them,
  and like all other codes that use them, has its own particular method for defining
  their directions"
  [a b c alpha beta gamma]
  (let [alph (cmato/* 0.0174532925199433 alpha)
        bet (cmato/* 0.0174532925199433 beta)
        gamm (cmato/* 0.0174532925199433 gamma)
        cosphi (cmato// (cmato/- (cmat/cos gamm) (cmato/* (cmat/cos alph) (cmat/cos bet))) (cmat/sin alph)  (cmat/sin bet))
        cphi (if (> cosphi 1) 1 cosphi)
        sinphi  (cmat/pow (cmato/- 1 (cmato/* cphi cphi)) 0.5)]
    [[(cmato/* a (cmat/sin bet) sinphi) (cmato/* a (cmat/sin bet) cosphi) (cmato/* a (cmat/cos bet))]
     [0  (cmato/* b (cmat/sin alph)) (cmato/* b (cmat/cos alph))]
     [0 0 c]]))


(defn parse-xmolout
  ([lines]
  (let [x (as-> lines x
                (second x)
                (strng/split x #"[ X]+"))]
  (basic/system (first x)
          (try (read-string (second x)) (catch Exception e (str "Reported time-step is not a number.  Check xmolout file for \"********\"." )))
          (apply reaxff-cell-params-lvs- (map read-string (drop 3 x)))
          (xyz-reax-iota->atoms (drop 2 lines)))))
  ([charge-column lines]
  (let [x (as-> lines x
                (second x)
                (strng/split x #"[ X]+"))]
  (basic/system (first x)
          (try (read-string (second x)) (catch Exception e (str "Reported time-step is not a number.  Check xmolout file for \"********\"." )))
          (apply reaxff-cell-params-lvs- (map read-string (drop 3 x)))
          (xyz-reax-iota->atoms charge-column (drop 2 lines))))))



(defn parse-xmolout-readable
  ([lines]
  (let [x (as-> lines x
                (second x)
                (strng/split x #"[ X]+"))]
  (basic/system (first x)
          (try (read-string (second x)) (catch Exception e (str "Reported time-step is not a number.  Check xmolout file for \"********\"." )))
          (apply reaxff-cell-params-lvs- (map read-string (drop 3 x)))
          (xyz-reax-iota->atoms-readable (drop 2 lines)))))
  ([charge-column lines]
  (let [x (as-> lines x
                (second x)
                (strng/split x #"[ X]+"))]
  (basic/system (first x)
          (try (read-string (second x)) (catch Exception e (str "Reported time-step is not a number.  Check xmolout file for \"********\"." )))
          (apply reaxff-cell-params-lvs- (map read-string (drop 3 x)))
          (xyz-reax-iota->atoms-readable charge-column (drop 2 lines))))))



(defn parse-xmoloutt
  ([lines]
  (let [x (as-> lines x
                (second x)
                (strng/split x #"[ X]+"))]

          (map read-string (drop 3 x)))))






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#_(defn get-atoms [atoms]
  "get atoms is used to turn a mol into a col of cols where the subcols are nothing
more than vectors of the species and the coordinates values.  In this function
we also test to see if the mol has values for the :pos key, if all of the atoms
do then we sort the mol by :pos."
  (let [sorted-atoms (if (not-any? (comp nil? :pos) atoms)
                       (sort-by :pos atoms)
                       atoms)
	f2 (comp #(concat (map float (take 3 %)) (drop 3 %)) :coordinates)]
    (map #(cons (:species %) (f2 %)) sorted-atoms)))




#_(defn write-xyz [mol]
  "Returns a seq of atoms as a string in xyz format.  This version of write-xyz
allows for writing a whole bunch of time steps (arranged in a col of cols of maps)
or a single time step which is a col of maps.

Usage:  Suppose (def test (xyz-str->atoms 'C 0 0 0 \n C 0.3333 0.6667 0')) then
the following would both give the same result (write-xyz test) => '2\n\n C 0 0 0 \n C 0.3333 0.6667 0'."
(str (count mol)  "\n\n" (utils/inter-cat-tree ["\n" "   "] (get-atoms mol)) "\n"))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



(defn parse-geo-HETATM-
  "This is a helper function that will be used in parse-geo that will parse the species
  and positions of all the atoms assocated with one BIOGRF/XTLGRF record in a geo file."
  [x]
  (->> x
        (strng/split-lines )
        (utils/grep #"HETATM" )
        (map (comp
               #(basic/new-atom (second %) (map read-string (take 3 (rest (rest %)))) nil nil nil nil (read-string (first %)))
               rest
               #(strng/split % #"\s+")) )))



(defn- parse-geo-CONECT-
  "This is a helper function that will be used in parse-geo that will parse the species
  and positions of all the atoms assocated with one BIOGRF/XTLGRF record in a geo file."
  [x]
  (let [l (->> x
       (strng/split-lines )
       (utils/grep #"CONECT" )
       (rest ))]
    (if (empty? l)
      (repeat nil)
       (map (comp
              #(basic/neigh-struct (map (comp dec read-string) %) nil nil)
               (partial drop 2)
              #(strng/split % #"\s+")) l))))




(defn parse-geo
"This is currently not full featured, and will only read in the positions of the atoms."
[filename]
  (map #(gmol/col->mol (parse-geo-HETATM- %) :neigh (parse-geo-CONECT- %))
    (utils/lazy-chunk-file filename #"BIOGRF|XTLGRF")))



(defn parse-geo
"This is currently not full featured, and will only read in the positions of the atoms."
[filename]
  (map parse-geo-HETATM-
    (utils/lazy-chunk-file filename #"BIOGRF|XTLGRF")))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn parse-QE-pw
  "This will produce a col of atomic coordinates from the output file of Quantum
Espresso's PW.x program.  As far as I can tell this only works with relaxation results."
  [filename]
  (map (comp
         xyz-str->atoms
         #(strng/join utils/endline %)
         #(utils/grep #"\s*[A-Z][a-zA-Z0-9]*([ \t]+)(-*\d+\.\d*[DE]?-?\d*[ \t]*){3}" %)
         strng/split-lines)
    (rest (utils/lazy-chunk-file filename #"ATOMIC_POSITIONS \([alat|bohr|angstrom|crystal]*\)\n"))))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   MDL Molfile  (.mol)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn- mdl-mol->atoms-
  "This will parse a string into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.

"
  [lines]
  (mapv (fn [x y] (#(basic/new-atom (.intern (nth % 3)) (cmat/matrix  (map read-string (take 3 %))) nil  nil nil nil y)
                (strng/split (strng/trim x) #"\s+")))
                  lines (iterate inc 0)))


(defn- mdl-mol->bonds-
  "This will parse a string into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.

"
  [atoms neighlines]
  (let [fneigh #(vector (first %) (basic/neigh-struct (second %)
                             (:species (gmol/mol-nth atoms (second %)))
                             (cmat/length (gmol/mol-vector atoms (first %) (second %)))
                             (gmol/mol-vector atoms (first %) (second %))))]
        (->> neighlines
            (r/map (fn [x](strng/split (strng/triml x) #"\s+")))
            (r/map (partial take 2))
            (r/map (partial map (comp dec read-string)))
            (r/map #(vector % (reverse %)))
            (into [])
            (utils/flatten-n 1)
            (set )
            (vec )
            (map fneigh))))







(defn parse-mdl-mol
  "This function parses the positions and neighbor information of atoms
in a MDL *.mol file.  The same basic information is contained in *.sdf files.
Assuming that the *.mol file information is at the beginning of the *.sdf this
function will parse *.sdf files."
  [filename]
  (let [lines (flatten (utils/clean-parse-empty filename))
        natoms (as-> lines x
                      (nth x 3)
                      (strng/split x #"[ X]+")
                      (first x)
                      (read-string x))
        nbonds (as-> lines x
                      (nth x 3)
                      (strng/split x #"[ X]+")
                      (second x)
                      (read-string x))
      atomss (mdl-mol->atoms- (take natoms (drop 4 lines)))]
      (gmol/find-assoc-in-vec [:pos :neigh]
      (map #(vector (first %) (map second (second %)))
        (group-by first
        (mdl-mol->bonds- atomss (take nbonds (drop (+ 4 natoms) lines)))))
      atomss)))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  POSCAR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn parse-simple-poscar
  "Poscars is a file format used by VASP to create specify the
   lattice coordinates and atomic positions of a system.  This "
  [filename]
  (let [lines (flatten (utils/clean-parse-empty filename))
        lvss (as-> lines x
                    (second x)
                    (read-string x))
        lvs  (as-> lines x
                    (take 5 x)
                    (take-last 3 x)
                    (strng/join " " x)
                    (strng/split x #"\s")
                    (map read-string x)
                    (cmato/* lvss x)
                    (partition-all 3 3 x))
        natoms (as-> lines x
                    (nth x 6)
                    (read-string x))
        species (as-> lines x
                    (nth x 5)
                    (strng/split x #"\s")
                    (map read-string x)
                    (map #(take %1 (repeat (str "S" %2))) x (iterate inc 0))
                    (flatten x))
        dircar (as-> lines x
                    (nth x 7))]
        (basic/system (first lines) 0  lvs
        (cond
             (= dircar "Cartesian")
             (as-> (drop 8 lines) x
                   (take natoms x)
                   (strng/join " " x)
                   (strng/split x #"\s")
                   (map read-string x)
                   (partition-all 3 3 x)
                   (map #(basic/new-atom %1 %2 nil nil nil nil %3) species x (iterate inc 0)))
             (= dircar "Direct")
             (as-> (drop 8 lines) x
                   (take natoms x)
                   (strng/join " " x)
                   (strng/split x #"\s")
                   (map read-string x)
                   (partition-all 3 3 x)
                   (map #(cmato/+ (cmato/* (first %) (first lvs))
                                 (cmato/* (second %) (second lvs))
                                 (cmato/* (last %) (last lvs))) x)
                   (map #(basic/new-atom %1 %2 nil nil nil nil %3) species x (iterate inc 0)))))))


(defn parse-simple-poscar-s
  "Poscars is a file format used by VASP to create specify the
   lattice coordinates and atomic positions of a system.  This "
  [filename]
  (let [lines (flatten (utils/clean-parse-empty filename))
        lvss (as-> lines x
                    (second x)
                    (read-string x))
        lvs  (as-> lines x
                    (take 5 x)
                    (take-last 3 x)
                    (strng/join " " x)
                    (strng/split x #"\s")
                    (map read-string x)
                    (cmato/* lvss x)
                    (partition-all 3 3 x))
        natoms (as-> lines x
                    (nth x 6)
                    (read-string x))
        species (repeat "C")
        dircar (as-> lines x
                    (nth x 7))]
        (utils/serialize
        (basic/system (first lines) 0 lvs
        (cond
             (= dircar "Cartesian")
             (as-> (drop 8 lines) x
                   (take natoms x)
                   (strng/join " " x)
                   (strng/split x #"\s")
                   (map read-string x)
                   (partition-all 3 3 x)
                   (map #(basic/new-atom %1 %2 nil nil nil nil %3) species x (iterate inc 0)))
             (= dircar "Direct")
             (as-> (drop 8 lines) x
                   (take natoms x)
                   (strng/join " " x)
                   (strng/split x #"\s")
                   (map read-string x)
                   (partition-all 3 3 x)
                   (map #(cmato/+ (cmato/* (first %) (first lvs))
                                 (cmato/* (second %) (second lvs))
                                 (cmato/* (last %) (last lvs))) x)
                   (map #(basic/new-atom %1 %2 nil nil nil nil %3) species x (iterate inc 0))))))))







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  LAMMPS Trajectory File Parser (.lammpstrj)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn- parse-lammpstrj-box-bounds
  "Parse the BOX BOUNDS line to extract box bounds.
  Format: ITEM: BOX BOUNDS pp pp pp
  Followed by lines like: 4.0833781541581260e-01 1.1529322184584935e+01"
  [lines]
  (let [bounds-lines (take 3 lines)
        parsed (mapv (fn [line]
                       (let [parts (strng/split (strng/trim line) #"\s+")]
                         (if (>= (count parts) 2)
                           [(read-string (first parts)) (read-string (second parts))]
                           nil)))
                     bounds-lines)]
    {:xlo (nth (nth parsed 0) 0)
     :xhi (nth (nth parsed 0) 1)
     :ylo (nth (nth parsed 1) 0)
     :yhi (nth (nth parsed 1) 1)
     :zlo (nth (nth parsed 2) 0)
     :zhi (nth (nth parsed 2) 1)}))




(defn- box-bounds->lvs
  "Convert LAMMPS BOX BOUNDS to lattice vectors.
  LAMMPS BOX BOUNDS format:
  ITEM: BOX BOUNDS pp pp pp
  xlo xhi
  ylo yhi
  zlo zhi
  
  Returns lattice vectors as a vector of 3 coordinate vectors.
  
  THIS ONLY WORKS FOR RECTANGULAR BOXES."
  [boxbounds]
  [[(- (:xhi boxbounds) (:xlo boxbounds)) 0.0 0.0]
   [0.0 (- (:yhi boxbounds) (:ylo boxbounds)) 0.0]
   [0.0 0.0 (- (:zhi boxbounds) (:zlo boxbounds))]])



(defn- lammps-atom-line->atom
  "Convert a LAMMPS atom line to an atom structure.
  Dynamically handles the column header to support various output formats."
  [line column-map pos-idx]
  (let [parts (strng/split (strng/trim line) #"\s+")
        element-idx (get column-map "element")
        x-idx (get column-map "x")
        y-idx (get column-map "y")
        z-idx (get column-map "z")
        q-idx (get column-map "q")
        bonds-idx (get column-map "c_bonds_count")]
    (when (and (not-any? nil? [element-idx x-idx y-idx z-idx])
               (>= (count parts) (max element-idx x-idx y-idx z-idx)))
      (let [element (when element-idx (nth parts element-idx))
            x (when x-idx (read-string (nth parts x-idx)))
            y (when y-idx (read-string (nth parts y-idx)))
            z (when z-idx (read-string (nth parts z-idx)))
            q (when q-idx (try (read-string (nth parts q-idx)) (catch Exception _ nil)))
            bonds (when bonds-idx (try (read-string (nth parts bonds-idx)) (catch Exception _ nil)))]
        (when element
          (basic/new-atom (.intern element)
                         (cmat/matrix [x y z])
                         q
                         nil
                         nil
                         bonds
                         pos-idx))))))


(defn- parse-lammpstrj-timestep
  "Parse a single LAMMPS timestep from a sequence of lines.
  Lines should be the chunk for one timestep."
  [lines]
  (try
    (let [lines-vec (vec (map strng/trim lines))
          timestep (when (and (> (count lines-vec) 1)
                              (= "ITEM: TIMESTEP" (nth lines-vec 0)))
                     (try (read-string (nth lines-vec 1)) (catch Exception _ nil)))
          natoms-line-idx (first (utils/positions #(= % "ITEM: NUMBER OF ATOMS") lines-vec))
          natoms (when natoms-line-idx
                   (try (read-string (nth lines-vec (+ natoms-line-idx 1))) (catch Exception _ nil)))
          box-line-idx (first (utils/positions #(= % "ITEM: BOX BOUNDS pp pp pp") lines-vec))
          box-bounds (when box-line-idx
                       (parse-lammpstrj-box-bounds (drop (+ box-line-idx 1) lines-vec)))
          atoms-line-idx (first (utils/positions #(strng/starts-with? % "ITEM: ATOMS") lines-vec))
          atoms-header (when atoms-line-idx (nth lines-vec atoms-line-idx))
          column-names (when atoms-header
                         (let [header-parts (strng/split atoms-header #"\s+")
                               cols (drop 2 header-parts)]
                           (zipmap cols (range (count cols)))))
          atom-lines (when (and atoms-line-idx natoms)
                       (take natoms (drop (+ atoms-line-idx 1) lines-vec)))]
      (when (and timestep natoms)
        {:timestep timestep
         :natoms natoms
         :lvs (box-bounds->lvs box-bounds)
         :column-names column-names
         :mol (filterv some? (mapv (fn [line idx] (lammps-atom-line->atom line column-names idx))
                                     atom-lines
                                     (range natoms)))}))
    (catch Exception e
      (println "Error parsing timestep:" e)
      nil)))









(defn parse-lammpstrj
  "Parse a LAMMPS trajectory file and return a vector of timesteps.
  Each timestep contains: :timestep, :natoms, :box-bounds, :column-names, and :atoms

  This function should only be used to parse small lammpstrj files.
	
  Usage: (parse-lammpstrj \"/path/to/file.lammpstrj\")"
  [filename]
  (let [lines (with-open [rdr (clojure.java.io/reader filename)]
                (into [] (line-seq rdr)))
        timestep-indices (into [] (utils/positions #(= % "ITEM: TIMESTEP") lines))
        timestep-ranges (partition 2 1 (conj (vec timestep-indices) (count lines)))]
    (filterv some? (mapv (fn [[start end]]
                           (parse-lammpstrj-timestep (subvec lines start end)))
                         timestep-ranges))))






(defn lammpstrj-chunks
  "Create foldable chunks directly from a LAMMPS trajectory file.
  Returns a foldable collection of timestep lines.

 This function should be used to parse large lammpstrj files.
  
  Usage: (lammpstrj-chunks filename start stop taken)
  taken means take every nth recorded time step
  "
  [filename start stop taken]
  (let [lines (vec (iota/seq filename))
        timestep-indices (into [] (utils/positions #(re-find #"ITEM: TIMESTEP" %) lines))
        ;; Filter timesteps in the desired range
        filtered-indices (filterv (fn [idx]
                                    (when (< (inc idx) (count lines))
                                      (let [ts-num (try (read-string (strng/trim (nth lines (inc idx))))
                                                        (catch Exception _ -1))]
                                        (and (>= ts-num start) (<= ts-num stop)))))
                                  timestep-indices)
        ;; Take every taken recorded time step, (i.e., take every nth recorded time step)
        sampled-indices (take-nth taken filtered-indices)
        ;; Create ranges for each timestep to the next ITEM: TIMESTEP
        timestep-ranges (map (fn [idx]
                               (let [next-idx (first (drop-while #(<= % idx) timestep-indices))]
                                 (if next-idx
                                   [idx next-idx]
                                   [idx (count lines)])))
                             sampled-indices)]
    (r/map (fn [[start end]]
             (subvec lines start end))
           timestep-ranges)))




(defn parse-lammpstrj-for-foldable
  "Parse a LAMMPS timestep chunk for use with foldable-chunks.
  Returns a map compatible with downstream processing like parse-xmolout.

  This function should be used to parse large lammpstrj files.

  Similar interface to parse-xmolout but for LAMMPS files."
  [lines]
  (let [parsed (parse-lammpstrj-timestep lines)]
    (when parsed
      {:timestep (:timestep parsed)
       :mol (:mol parsed)
       :lvs (:lvs parsed)})))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



#_(->> (foldable-chunks "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C2F1256ud400MPa/N3/xmolout")
(r/map (partial drop 2))
     (r/map xyz-iota->atoms)
       (r/map #(shift [0 0 -200] %) )
       (r/map out/write-xyz)
     (r/map #(utils/append-file "/Users/chadjunkermeier/Desktop/graphene2.xyz" %))
     (into []))





;(require '[greenwood.atomic-structure-output :as out])
#_(->> (foldable-chunks "/Users/chadjunkermeier/Desktop/graphene.xyz" )
(r/map (partial drop 2))
     (r/map xyz-iota->atoms)
     (r/map out/write-xyz)
     (r/map #(utils/append-file "/Users/chadjunkermeier/Desktop/graphene2.xyz" %))
     (into []))

;(parse-xyz "/Users/chadjunkermeier/Dropbox/Fgraphene-NEB/Results/NEB/FF/F0F12/F0F12.xyz")

#_(->> (foldable-chunks "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/all300MPa/N10/xmolout" (reax-index-timesteps 1532000 1532100 100 (* 2 576)))
     (r/map parse-xmolout)
     (r/map :mol)
     (r/map (partial dfdf {:name 2} ))
     (r/map write-xyz)
     (r/map #(utils/append-file "/Users/chadjunkermeier/Desktop/graphene.xyz" %))
       (into []))

#_(->> (take-foldable-chunks "/Users/chadjunkermeier/Desktop/xmolout" 1)
(r/map (partial drop 2))
     (r/map xyz-iota->atoms)
     (into []))


;this will parse every other record in /Users/junky/3netC10.lammpstrj
#_(->> (xyz/lammpstrj-chunks "/Users/junky/3netC10.lammpstrj" 0 300 2)
     (r/map xyz/parse-lammpstrj-for-foldable)
     (r/map :mol)
     (r/map out/write-xyz)
     (r/map #(utils/append-file "/Users/junky/Desktop/graphene.xyz" %))
     (into []))

;this will parse every record in /Users/junky/3netC10.lammpstrj
#_(->> (xyz/lammpstrj-chunks "/Users/junky/3netC10.lammpstrj" 0 300 1)
     (r/map xyz/parse-lammpstrj-for-foldable)
     (r/map :mol)
     (r/map out/write-xyz)
     (r/map #(utils/append-file "/Users/junky/Desktop/graphene.xyz" %))
     (into []))

;(def netc (xyz/parse-lammpstrj "/Users/junky/3netC10.lammpstrj"))

