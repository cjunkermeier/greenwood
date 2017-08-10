(ns reaxff
  (:require
          [clojure.string :as cstrng]
          [clojure.pprint :as cpp]
          [clojure.set :as cset]
          [greenwood.basics :as bas]
          [greenwood.empirical-data :as ed]
          [greenwood.math :as jmath]
          [greenwood.mol :as jmol]
          [greenwood.utils :as utils]
          [greenwood.xyz :as xyz])
  (:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix)
  (:use clojure.core.matrix.operators))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  WRITING INPUT GEOMETRY RECORDS FOR REAXFF
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn write-reac-HETATM [mol]
  "Defines atom type and atom position. In this order, the ATOM INFO consists of
  the atom number, the atom type, the x, y and z-coordinates of the atom in Å,
  the force field type (same as atom type for ReaxFF), two switches not used
  by ReaxFF and the atom partial charge. ReaxFF does not use these partial charges.

  transduce.xyz/atom-pos must be applied to mol in order for this function to work.

  ReaxFF always uses Cartesian (not fractional) coordinates to define atom positions.
  This uses the fortran format (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)"
  (let [f #(cpp/cl-format nil "~6A~1@T~5D~1@T~5A~1@T~3A~1@T~1A~1@T~5A~10,5,0,'*F~10,5,0,'*F~10,5,0,'*F~1@T~5A~3D~2D~1@T~8,5,0,'*F"
             "HETATM"
             (inc (:pos %))
             (:species %)
              "" "" ""
             (first (:coordinates %))
             (second (:coordinates %))
             (utils/third (:coordinates %))
             (:species %)
             1
             1
             0.00)]
    (str (cstrng/join utils/endline (map f (sort-by :pos mol))))))




(defn write-bond-restraint
  "(a6,2i4,f8.4,f8.2,f8.5,f10.7)"
  ([mol a1 a2]
   (write-bond-restraint mol a1 a2 500 1))
 ([mol a1 a2 distance]
 (let [s "~15A~4D~4D~8,4,0,'*F~8,2,0,'*F~8,5,0,'*F~10,7,0,'*F"
      dist (distance (:coordinates (nth mol a1)) (:coordinates (nth mol a2)))]
      (cpp/cl-format nil s "BOND RESTRAINT" (inc a1) (inc a2) distance  500 1  0.0000000  0  0)))
 ([mol a1 a2 R1 R2]
 (let [s "~15A~4D~4D~8,4,0,'*F~8,2,0,'*F~8,5,0,'*F~10,7,0,'*F"
      dist (distance (:coordinates (nth mol a1)) (:coordinates (nth mol a2)))]
      (cpp/cl-format nil s "BOND RESTRAINT" (inc a1) (inc a2) dist  R1  R2  0.0000000  0  0)))
 ([mol a1 a2 R1 R2 distance-add]
 (let [s "~15A~4D~4D~8,4,0,'*F~8,2,0,'*F~8,5,0,'*F~10,7,0,'*F"
      dist (distance (:coordinates (nth mol a1)) (:coordinates (nth mol a2)))]
      (cpp/cl-format nil s "BOND RESTRAINT" (inc a1) (inc a2) (+ dist distance-add)  R1  R2  0.0000000  0  0))))



(defn atom-atom-bond-restraint
  "In this case you give the two records of the atoms that you want to have a bond restraint between.
  (a6,2i4,f8.4,f8.2,f8.5,f10.7)"
 [a1 a2]
 (let [s "~15A~4D~4D~8,4,0,'*F~8,2,0,'*F~8,5,0,'*F~10,7,0,'*F"
      dist (distance (:coordinates a1) (:coordinates a2))]
      (cpp/cl-format nil s "BOND RESTRAINT" (inc (:pos a1)) (inc (:pos a2)) dist  500.00  0.5000  0.0000000  0  0)))


(defn write-angle-restraint
  "FORMAT ANGLE RESTRAINT (a16,3i4,2f8.2,f8.4,f9.6)"
 [mol a1 a2 a3]
 (let [s "~16A~4D~4D~4D~8,2,0,'*F~8,2,0,'*F~8,4,0,'*F~9,6,0,'*F"
      angle (jmath/round-decimal 0.001 (ed/radians->degrees (jmath/three-point-angle (:coordinates (jmol/mol-nth mol a1)) (:coordinates (jmol/mol-nth mol a2))
                         (:coordinates (jmol/mol-nth mol a3)))))]
      (cpp/cl-format nil s "ANGLE RESTRAINT" (inc a1) (inc a2) (inc a3) angle  200.00 1.00000 0.0000)))

(defn write-angle-restraint-angle
  "FORMAT ANGLE RESTRAINT (a16,3i4,2f8.2,f8.4,f9.6)"
 ([a1 a2 a3 angle]
 (let [s "~16A~4D~4D~4D~8,2,0,'*F~8,2,0,'*F~8,4,0,'*F~9,6,0,'*F"]
      (cpp/cl-format nil s "ANGLE RESTRAINT" (inc a1) (inc a2) (inc a3) angle  200.00 1.00000 0.0000)))
 ([a1 a2 a3 angle strength]
 (let [s "~16A~4D~4D~4D~8,2,0,'*F~8,2,0,'*F~8,4,0,'*F~9,6,0,'*F"]
      (cpp/cl-format nil s "ANGLE RESTRAINT" (inc a1) (inc a2) (inc a3) angle  strength 1.00000 0.0000))))






(defn write-torsion-restraint
  "FORMAT TORSION RESTRAINT (18a,4i4,2f8.2,f8.4,f9.6)"
 [mol a1 a2 a3 a4]
 (let [s "~18A~4D~4D~4D~4D~8,2,0,'*F~8,2,0,'*F~8,4,0,'*F~9,6,0,'*F"
      angle (jmath/round-decimal 0.001 (ed/radians->degrees (jmol/mol-find-dihedral mol [a1 a2 a3 a4])))]
      (cpp/cl-format nil s "TORSION RESTRAINT " (inc a1) (inc a2) (inc a3) (inc a4) angle  5000.00 5.00000 0.0000)))


(defn write-torsion-restraint-angle
  "FORMAT TORSION RESTRAINT (18a,4i4,2f8.2,f8.4,f9.6)"
 [mol a1 a2 a3 a4 angle]
 (let [s "~18A~4D~4D~4D~4D~8,2,0,'*F~8,2,0,'*F~8,4,0,'*F~9,6,0,'*F"]
      (cpp/cl-format nil s "TORSION RESTRAINT " (inc a1) (inc a2) (inc a3) (inc a4) angle  5000.00 5.00000 0.0000)))

(defn write-FIXATOMS
    "Formats FIXATOMS call in geo files"
  [a1 a2]
 (let [s "~8A~6D~6D"]
  (cpp/cl-format nil s "FIXATOMS" (inc a1) (inc a2))))

(defn write-descrp
  "System description. This description can be used in trainset.in to define a force field training set.
  Fortran format: (a7,a40)"
  [string]
  (let [s "~7A~40A"]
    (cpp/cl-format nil s "DESCRP" (apply str (take 40 string)))))

(defn write-remark
  "Remarks. Multiple REMARK lines are allowed.
  Fortran format: (a7,a40)"
  [string]
  (let [s "~7A~40A"]
    (cpp/cl-format nil s "REMARK" (apply str (take 40 string)))))

(defn write-crystx
  "Defines cell lengths (in Å) and cell angles (in degrees) for periodic system.
  Fortran format: (8a,6f11.5)"
  [a b c alpha beta gamma]
  (let [s "~8A~11,5,0,'*F~11,5,0,'*F~11,5,0,'*F~11,5,0,'*F~11,5,0,'*F~11,5,0,'*F"]
    (cpp/cl-format nil s "CRYSTX" a b c alpha beta gamma)))

(defn write-ff-charges
  [mol name-str]
  (let [f (fn [x y] (-  ((comp bigdec  (partial apply unchecked-add-int) ed/atomic-charge ed/atomic-numbers :species) y) (bigdec (x y))))
        data (if ((comp :total :charge first) mol)
                 (map (partial f (comp :total :charge)) mol)
                 (map (partial f :charge) mol))]
    (join utils/endline (map #(str name-str " 0.1 " %1 " " %2) (iterate inc 1) data))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  READING IN REAXFF OUTPUT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn parse-ixmolo1
  "This will utils/parse a str into the atoms struct.  Note that the string should start
with the first atom, not with the number of atoms in the system.  Also, this
assumes that there is a newline character between atoms.  Further, this assumes
  that the user set ixmolo=1.

Thus if: (def test 'C 0 0 0 \n C 0.3333 0.6667 0 1')
then the usage would be (xyz-str->atoms test)."
  [string]
  (let [lines  (cstrng/split-lines string)]
    (map (comp #(bas/new-atom (.intern (first %1)) ((comp (partial map read-string) (partial take 3) next) %1) nil nil nil (read-string (nth %1 4)) %2)
                   #(cstrng/split % #"\s+")
                    #(cstrng/trim %1)) lines (iterate inc 0))))

(defn parse-ixmolo345
  [string]
  (let [lines (cstrng/split-lines string)]
    (map (comp
           #(bas/new-atom (.intern (first %1))
                      ((comp (partial map read-string) (partial take 3) next) %1)
                      (read-string (last %1))
                      nil nil
                      (read-string (nth %1 4))
                      %2)
           #(cstrng/split %1 #"\s+")
           #(cstrng/trim %1))
          lines (iterate inc 0))))


(defn parse-ixmolo345
  [string]
  (let [lines (cstrng/split-lines string)]
    (map (comp
           #(bas/new-atom (.intern (first %))
                      ((comp (partial map read-string) (partial take 3) next) %)
                      (read-string (last %))
                      nil nil
                      (read-string (nth % 4))
                      nil)
           #(cstrng/split % #"\s+")
           cstrng/trim )
          lines )))






#_(defn parse-xmolout
  "This lazily utils/parses an xmolout-file.  It also includes the xyz/atom-pos functionality.
This produces a col of cols, where each of the sub-cols is a time step. This assumes
  that the user set ixmolo=1, if they didn't use xyz-utils/parser/utils/parse-xyz.

Usage: (second (utils/parse-xmolout PATH 3)), where PATH is a string containing the path
to some xmolout file, and 3 is the value that ixmolo is set to in the reaxff control file."
  ([filename ixmolo]
  (let [f ({1 utils/parse-ixmolo1 3 utils/parse-ixmolo345 4 utils/parse-ixmolo345 5 utils/parse-ixmolo345} ixmolo)]
  (mapv f
    (utils/lazy-chunk-file filename #"\n[ \t]*\d+[ \t]*\n.*\n|^[ \t]*\d+[ \t]*\n.*\n"))))
  ([filename ixmolo steps-f]
  (let [f ({1 utils/parse-ixmolo1 3 utils/parse-ixmolo345 4 utils/parse-ixmolo345 5 utils/parse-ixmolo345} ixmolo)]
  (mapv f
   (steps-f (utils/lazy-chunk-file filename #"\n[ \t]*\d+[ \t]*\n.*\n|^[ \t]*\d+[ \t]*\n.*\n"))))))





(defn parse-xmolout-params
  [filename]
  (as-> filename x
       (utils/clean-parse x)
       (second )
      (cstrng/split x #"\s+")
      (drop 3 x)
       (map read-string x)))















(defrecord fort71
   [Iter Nmol Nmol2 Epot Ekin Etot T Eaver-block Eaver-tot Taver Tmax Press sdev-Epot sdev-Eaver Tset Timestep RMSG totaltime])

(defn parse-71
  [folder]
  (utils/parse (str folder "/fort.71")
    (comp #(apply ->fort71 %) (partial map read-string) #(cstrng/split % #"\s+")) 1))


(defn analyze-74
  [filename]
  (let [raw (utils/clean-parse filename)
        averageIT (map (comp read-string #(nth % 5) #(cstrng/split % #"\s+")) raw)]
    (* 1. (jmath/average averageIT))))


(defn analyze-99
  [filename]
  (let [raw (utils/clean-parse filename)
        f (fn [x] (reduce + (map (comp read-string first #(take-last 2 %) #(cstrng/split % #"\s+")) x)))]
    (hash-map :a (f (utils/grep  #" a:" raw))
              :charge (f (utils/grep #"Charge atom:" raw))
              :bond (f (utils/grep  #"Bond distance:" raw))
              :angle (f (utils/grep #"Valence angle:" raw))
              :energy (f (utils/grep #"Energy " raw)))))


(defn analyze-99-percent-diff
  [filename]
  (let [raw (utils/clean-parse filename)
        g #(/ (abs (- (first %) (second %)) )  (abs (second %) ))
        f (fn [x] (reduce + (map (comp
                                  g
                                  (partial map read-string)
                                  (partial take 2)
                                  (partial take-last 5)
                                  #(cstrng/split % #"\s+")) x)))]
        (hash-map :a (f (utils/grep  #" a:" raw))
              :charge (f (utils/grep #"Charge atom:" raw))
              :bond (f (utils/grep  #"Bond distance:" raw))
              :angle (f (utils/grep #"Valence angle:" raw))
              :energy (f (utils/grep #"Energy " raw)))))




#_(defn analyze-99-charge
  [filename]
  (let [raw (utils/clean-parse filename)
        g #(/ (expt (- (first %) (second %)) 2)  (expt (second %) 2))
        f (fn [x]  (map (comp
                                  g
                                  (partial map read-string)
                                  (partial take 2)
                                  (partial take-last 5)
                                  #(cstrng/split % #"\s+")) x))
        h (fn [x]  (map (comp read-string first #(take-last 2 %) #(cstrng/split % #"\s+")) x))]
        (map #(hash-map
              :percentdiff %1    :error %2)   (f (utils/grep #"Charge atom:" raw))   (h (utils/grep #"Charge atom:" raw)))))


(defn parse-mep
  [filename n]
  (let [f (->> filename
      utils/clean-parse-empty
     first
      utils/ls-to-ll
      utils/transpose
      last
      (map ed/kcalpermol->eV))]
    (map #(- % (nth f n)) f)))


(defn NEB-barriers-heights
  [en-col]
  (let [a (ed/kcalpermol->eV (#(-  (apply max %) (first %)) en-col))
        b (ed/kcalpermol->eV (#(-  (apply max %) (last %)) en-col))]
      (hash-map :-> a :<- b )))


(defn reaxff-cell-params-lvs
  "Cell parameters should be banned!  But since they are not, and Reaxff uses them,
  and like all other codes that use them, has its own particular method for defining
  their directions"
  [a b c alpha beta gamma]
  (let [alph (ed/degrees->radians alpha)
        bet (ed/degrees->radians beta)
        gamm (ed/degrees->radians gamma)
        cosphi (/ (- (cos gamm) (* (cos alph) (cos bet))) (sin alph)  (sin bet))
        cphi (if (> cosphi 1) 1 cosphi)
        sinphi  (pow (- 1 (* cphi cphi)) 0.5)]
    [[(* a (sin bet) sinphi) (* a (sin bet) cosphi) (* a (cos bet))]
     [0  (* b (sin alph)) (* b (cos alph))]
     [0 0 c]]))


(reaxff-cell-params-lvs 5.21200    5.21200   15.00000   90.00000   90.00000  120.00000)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      Beginning of definitions used for parsing molfra.out files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn  molfra-names
  "produces a list of all of the molecules given in molfra.out"
[filename]
(with-open [rdr (clojure.java.io/reader  filename)]
  (let [a (utils/grep  #"\d+\s+\d+\sx\s+[1-9a-zA-Z]*" (line-seq rdr))]
    (cset/difference (set
     (map (comp  #(nth % 4) #(cstrng/split % #"\s+")) a)) #{"x" "atoms:" "Total"}))))


(defn molfra-iteration
  "produces a list of all of the iterations for which there is output in molfra.out
  This is needed because I often find that reaxff writes the output of some of
  the iterations twice."
  [filename]
  (with-open [rdr (clojure.java.io/reader  filename)]
  (let [a (utils/grep  #"\d+\s+\d+\sx\s+[1-9a-zA-Z]*" (line-seq rdr))]
     (sort (set (map (comp read-string first #(cstrng/split % #"\s+") cstrng/trim) a))))))

(defn molfra-iteration
  "produces a list of all of the iterations for which there is output in molfra.out
  This is needed because I often find that reaxff writes the output of some of
  the iterations twice."
  [filename]
  (with-open [rdr (clojure.java.io/reader  filename)]
  (let [a (utils/grep  #"\d+\s+\d+\sx\s+[1-9a-zA-Z]*" (line-seq rdr))]
     (sort (set (map #(try ((comp read-string first (fn [x] (cstrng/split x #"\s+")) cstrng/trim) %) (catch java.lang.ClassCastException e (str %))) a))))))


(defn molfra-iteration
  "produces a list of all of the iterations for which there is output in molfra.out
  This is needed because I often find that reaxff writes the output of some of
  the iterations twice."
  [filename]
  (with-open [rdr (clojure.java.io/reader  filename)
                    wrtr (clojure.java.io/writer  "/Users/chadjunkermeier/Desktop/bbbblah.txt" :append false)]
  (let [a (utils/grep  #"\d+\s+\d+\sx\s+[1-9a-zA-Z]*" (line-seq rdr))]
     (sort (set (map #(do (.write wrtr (str % "\n"))
                             ((comp read-string first (fn [x] (cstrng/split x #"\s+")) cstrng/trim) %)) a))))))



(defn  molfra-molecule
[filename molname]
(with-open [rdr (clojure.java.io/reader  filename)]
 (let [a (sort (set (utils/grep (re-pattern (str "\\d+\\s+\\d+\\s+x\\s+" molname "\\s+")) (line-seq rdr))))]
    (map (comp #(map read-string %) #(vector (first %) (second %)) #(cstrng/split % #"\s+") cstrng/trim) a))))


(defn gg
  [k v mv]
  (let [f (partial jmol/find-assoc-in [:iteration k])]
  (loop [a (first v)
         b (rest v)
         c  (vec mv)]
    (if (empty? b)
      (f a c)
      (recur (first b) (rest b) (f a c))))))


(defn molfra-analysis
  "This is used to utils/parse the molfra*****.out files produced by reaxff.  There are
  several possible problems with the formatting produced by reaxff and by Adri's
  molfra_annal fortran programs.  This function is designed to compensate for
  those problems and correctly produce what would be expected from molfra_annal.

  Usage (molfra-analysis '/Users/junky/Dropbox/PennState-data/BATTERY/EC_interface/200K/amorphous_EC200/molfra_ig.out')"
  [filename]
    (let [kw (molfra-names filename)
          its  (molfra-iteration filename)
          hm (zipmap (concat [:iteration] (map keyword kw)) (repeat 0))
          colhm (map (partial assoc hm :iteration ) its)]
      (loop [x colhm
             y (first kw)
             z (rest kw)]
        (if (empty? z)
         (gg (keyword y) (molfra-molecule filename y) x)
          (recur (gg (keyword y) (molfra-molecule filename y) x) (first z) (rest z))))))


(defn molfra-analysis-selective
  "This is used to utils/parse the molfra*****.out files produced by reaxff.  There are
  several possible problems with the formatting produced by reaxff and by Adri's
  molfra_annal fortran programs.  This function is designed to compensate for
  those problems and correctly produce what would be expected from molfra_annal.

  Usage (molfra-analysis-selective '/Users/junky/Dropbox/PennState-data/BATTERY/EC_interface/200K/amorphous_EC200/molfra_ig.out' 'HF')"
  [filename molecules-keys]
    (let [its  (molfra-iteration filename)
          hm (zipmap (concat [:iteration] (map keyword molecules-keys)) (repeat 0))
          colhm (map (partial assoc hm :iteration ) its)]
      (loop [x colhm
             y (first molecules-keys)
             z (rest molecules-keys)]
        (if (empty? z)
         (gg (keyword y) (molfra-molecule filename y) x)
          (recur (gg (keyword y) (molfra-molecule filename y) x) (first z) (rest z))))))




(defn molfra-data
  "Sometimes you only want a couple of the molecules from a
  file and sometimes you only want one molecule from a bunch of files,
  this allows you to run molfra-molecule for only the data you want and
  then collect them in one molfra data set.

  Usage (def PATH '/Users/junky/Dropbox/PennState-data/Fgraphene/TEMP_INCREASE_periodic/')
        (def C8F (molfra-molecule (str PATH '/Temp_increase_C8F/molfra.out') 'F'))
        (def C4F (molfra-molecule (str PATH '/Temp_increase_C4F/molfra.out') 'F'))
        (def it (take 1000000/50 (iterate (partial + 50 ) 0)))
        (molfra-data it [:C8F :C4F] [C8F C4F])"
  [its kw molfra-molecule-data]
    (let [hm (zipmap (concat [:iteration] kw) (repeat 0))
          colhm (map (partial assoc hm :iteration ) its)]
      (loop [x colhm
             y (first kw)
             z (rest kw)
             yy (first molfra-molecule-data)
             zz (rest molfra-molecule-data)]
        (if (empty? z)
         (gg (keyword y) yy x)
          (recur (gg (keyword y) yy x) (first z) (rest z)(first zz) (rest zz))))))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;      End of definitions used for parsing molfra.out files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;











