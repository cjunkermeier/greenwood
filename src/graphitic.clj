(ns graphitic
  (:use
   greenwood.supercell
   greenwood.xyz
   )
  (:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix.operators
              clojure.core.matrix)
  (:require [greenwood.solution :as gsol]
            [greenwood.empirical-data :as ed]
            [greenwood.mol :as gmol]
            [clojure.set :as cset]
            [clojure.core.matrix :as cmat]
            [greenwood.neighbors :as gneigh]
           [greenwood.math :as gmath]
            [greenwood.contrib-math :as gcmath]
            [greenwood.basics :as basic]
            [greenwood.utils :as gutils]
            [greenwood.xyz :as xyz]
))




(comment   "The atomic gutils/positions, lattice vectors, reciprical lattice vectors all came from
Peres et al. 'Scanning Tunneling Microscopy currents on locally disordered graphene'
PRB 79 155442.

The honeycomb lattice has a unit cell represented in Fig. 1 by the
vectors a1 and a2, such that |a1| = |a2| = a, with a ≃ 2.461.
Where a0 = a/√3 ≃ 1.421 is the carbon-carbon distance.

Carbon-hydrogen bonds have a bond length of about 1.09 Å (1.09 × 10-10 m) and a
bond energy of about 413 kJ/mol.

The carbon–fluorine bond length is typically about 1.35 Angstrom.  The carbon–fluorine
bond length varies by several hundredths of an angstrom depending on the hybridization
of the carbon atom and the presence of other substituents on the carbon or even in atoms
farther away. These fluctuations can be used as indication of subtle hybridization
changes and stereoelectronic interactions.


the B-N bond length in a sheet is 1.45Å")


(defn a-one [a]
  "a is the a lattice constant"
  [(* a 0.5 (cmat/sqrt 3)), (* 0.5 a ), 0])

(defn a-two [a]
  "a is the a lattice constant"
  [(* a 0.5 (cmat/sqrt 3)), (* -0.5 a ), 0])

(defn a-three [a]
  "This is for use in quantum espresso."
  (- (a-two a) (a-one a)))


(defn graphene-primitive-unit-cell
  "The honeycomb lattice has a unit cell represented in Fig. 1 by the
vectors a1 and a2, such that |a1| = |a2| = a, with a ≃ 2.461.
Where a0 = a/√3 ≃ 1.421 is the carbon-carbon distance.
  Usage: (graphene-primitive-unit-cell 'C' 'C' 2.461)"
  [C1 C2 a]
  (hash-map :lvs [(a-one a) (a-two a) [0 0 30]]
   :mol [(basic/new-atom (.intern C1) (* 1/3 (+ (a-one a) (a-two a))) nil nil nil nil 0)
	(basic/new-atom (.intern C2) (* 2/3 (+ (a-one a) (a-two a))) nil nil nil nil 1)]))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; I used these in working with Quantum Espresso
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn QE-to-xyz
  "Used with my 2x2 Fgraphene paper."
  [mol alat nxn]
(let [a (ed/Bohr->Angstrom (* alat nxn))
      c (ed/Bohr->Angstrom 46.5)]
(atom-pos (fractional->cartesian mol [[a, 0.0, 0.0] [(* -0.5 a) (* 0.5 a (cmat/sqrt 3)) 0.0] [0.0 0.0 c]]))))



(defn QE-to-lvs
    "Used with my 2x2 Fgraphene paper."
  [alat nxn]
(let [a (ed/Bohr->Angstrom (* alat nxn))
c (ed/Bohr->Angstrom 46.5)
A-one [a, 0.0, 0.0]
A-two [(* -0.5 a) (* 0.5 a (greenwood.contrib-math/sqrt 3)) 0.0]
A-three [0.0 0.0 c]]
[A-one A-two A-three]))


(defn graphene-QE-unit-cell
  "The honeycomb lattice has a unit cell represented in Fig. 1 by the
vectors a1 and a2, such that |a1| = |a2| = a, with a ≃ 2.461.
Where a0 = a/√3 ≃ 1.421 is the carbon-carbon distance.
  Usage (graphene-primitive-unit-cell 'C' 'C' 2.461)"
  [C1 C2 a]
  (hash-map :lvs (QE-to-lvs a 1)
   :mol [(basic/new-atom (.intern C1) [0 0 0] nil nil nil nil 0)
   (basic/new-atom (.intern C2) (+ (* 1/3 (a-one a)) (* 2/3 (a-three a))) nil nil nil nil 1)]))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; rectangular armchair unit cell
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn armchair-uc
  "a is the C-C distance"
  [C1 C2 C3 C4 a]
  (let [c (cmat/cos (/ ed/pi 3))
        s (cmat/sin (/ ed/pi 3))]
  (vector
    (basic/new-atom C1 [(* 0.5 a c) (* 0.5 a s) 0] nil nil nil nil 1)
	(basic/new-atom C2 [(+ a (* 0.5 a c)) (* 0.5 a s) 0] nil nil nil nil 1)
	(basic/new-atom C3 [(+ a (* 1.5 a c)) (* 1.5 a s) 0] nil nil nil nil 1)
	(basic/new-atom C4 [(+ a a (* 1.5 a c)) (* 1.5 a s) 0] nil nil nil nil 1))))


(defn armchair-lvs
  [a]
  (vector
    [(* 2 a (+ 1 (cmat/cos (/ ed/pi 3)))) 0 0]
	  [0 (* 2 a (cmat/sin (/ ed/pi 3))) 0]
	  [0 0 20]))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; rectangular zigzag unit cell
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn zigzag-uc
  [C1 C2 C3 C4 a]
  (let [c (cmat/cos (/ ed/pi 6))
        s (cmat/sin (/ ed/pi 6))]
  (vector
    (basic/new-atom C1 [0 (* 0.5 a) 0] nil nil nil nil 1)
	(basic/new-atom C2 [0 (* 5/2 a) 0] nil nil nil nil 1)
	(basic/new-atom C3 [(* a c) (* a (+ 0.5 s)) 0] nil nil nil nil 1)
	(basic/new-atom C4 [(* a c) (* a (+ 1.5 s)) 0] nil nil nil nil 1))))



(defn zigzag-lvs
  [a]
  (vector
    [(* 2 a (cmat/cos (/ ed/pi 6))) 0 0]
	  [0 (* 2 a (+ 1 (cmat/sin (/ ed/pi 6)))) 0]
	  [0 0 30]))


;;;;;;;;;;;;;;;;; removing the line of C of the highest y ;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn remove-highest-y
  "The above supercells only make structures with an even number of atoms in the
  y-direction.  This allows there to be an odd number."
  [mol]
  (let [f #(> (- (nth (gmol/min-max-coordinates mol) 3) 0.2) (second %))]
    (gmol/mol-filter {:coordinates f} mol)))


;;;;;;;;;;;;;;;;; Figuring out how to place H-atoms ;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn two-point-unit-vec
         "Creates a unit vector that points from p-one to p-two in cartesian coord."
  [p-one p-two]
         (cmat/normalise (- p-two p-one)))



(defn find-edge-C
  "This determines if an atom is an edge atom."
  [mol]
  (filter #(>= 2 ((comp count  :neigh) %)) mol))


(defn find-edge-C-periodic
  "This determines if an atom is an edge atom."
  [mol lvs]
  (as->  (computation-supercell mol lvs 1 1 0) x
         (:mol x)
         (atom-pos x)
         (gneigh/neighbors x 0.1 1.8)
         (take (count mol) x)
  (filter #(>= 2 ((comp count :neigh) %)) x)))



(defn add-bonding-H
  "This works for binding H to graphene."
  [mol atomm]
  (let [numneigh ((comp count :neigh) atomm)
        cone (:coordinates atomm)
        ctwo (:coordinates (nth mol ((comp :npos first :neigh) atomm)))
        CH-length 1.09]
    (cond
      (= 1 numneigh)
      (let [unit-2-rotate (cmat/normalise (- ctwo cone))
            newz (+ cone [0 0 1])]
            (vector (basic/new-atom "H" (map + cone (map (partial * CH-length) (gmath/zrotation (* ed/pi 2/3) unit-2-rotate))) nil nil nil nil -1)
                      (basic/new-atom "H" (map + cone (map (partial * CH-length) (gmath/zrotation (* ed/pi -2/3) unit-2-rotate))) nil nil nil nil -1)))
      (= 2 numneigh)
      (let [cthree (:coordinates (nth mol ((comp :npos second :neigh) atomm)))]
         [(basic/new-atom "H" (+ cone (* CH-length (two-point-unit-vec (gmath/point-line-intersection ctwo cthree cone) cone))) nil nil nil nil -1)])
     (< 2 numneigh)
      (println "Atom number " (:pos atomm) "(at coordinates " (:coordinates atomm) ") has more than 2 gneigh/neighbors"))))



(defn add-H
  ([mol]
   (concat mol (flatten (map #(add-bonding-H mol %) (find-edge-C mol)))))
  ([mol lvs]
   (concat mol (flatten (map #(add-bonding-H (:mol (computation-supercell mol lvs 1 1 0)) %) (find-edge-C-periodic mol lvs))))))



;;;;;;;;;;;;;;;;;;;;;;;;; NANORIBBONS ;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-armchair-nanoribbon
"This will make a carbon nanoribbon that runs in the x-direction.
  The edges are armchair and have and have H atoms bonded appropriately.

CC-dist is the bond length between neighboring C atoms,
nxcells is the number of unit cells in the direction of the ribbon,
nycells is the number of unit cells in the direction perpendicular to the ribbon,
rtl is a boolean that determins if the C atoms with the largest y-values are discarded

Usage: (make-armchair-nanoribbon 'C 'C 'C 'C 1.421 2 2 true)
Usage: (make-armchair-nanoribbon 'C 'C 'C 'C 1.421 2 2 nil)"
[C1 C2 C3 C4 CC-dist nxcells nycells rtl]
(let [uclvs (armchair-lvs CC-dist)
         sc (if (true? rtl)
                (remove-highest-y (:mol (supercell (armchair-uc C1 C2 C3 C4 CC-dist)  uclvs nxcells nycells 1)))
                (:mol (supercell (armchair-uc C1 C2 C3 C4 CC-dist) uclvs nxcells nycells 1)))
      lvs [(* nxcells (first uclvs)) (+ [0.0 30.0 0.0] (* nycells (second uclvs))) [0.0 0.0 30.0]]]
  (basic/unitcell lvs
  (as-> sc x
    (atom-pos x)
    (gneigh/neighbors x 0.1 1.8)
    (add-H x lvs)))))







(defn make-zigzag-nanoribbon
"This will make a carbon nanoribbon that runs in the x-direction.
  The edges are zigzag and have and have H atoms bonded appropriately.

  CC-dist is the bond length between neighboring C atoms,
nxcells is the number of unit cells in the direction of the ribbon,
nycells is the number of unit cells in the direction perpendicular to the ribbon,
rtl is a boolean that determins if the C atoms with the largest y-values are discarded

Usage: (make-zigzag-nanoribbon 'C 'C 'C 'C 1.421 2 2 true)
Usage: (make-zigzag-nanoribbon 'C 'C 'C 'C 1.421 2 2 nil)"
[C1 C2 C3 C4 CC-dist nxcells nycells rtl]
(let [uclvs (zigzag-lvs CC-dist)
         sc (if (true? rtl)
                (remove-highest-y (:mol (supercell (zigzag-uc C1 C2 C3 C4 CC-dist)  uclvs nxcells nycells 1)))
                (:mol (supercell (zigzag-uc C1 C2 C3 C4 CC-dist) uclvs nxcells nycells 1)))
      lvs [(* nxcells (first uclvs)) (+ [0.0 30.0 0.0] (* nycells (second uclvs))) [0.0 0.0 30.0]]]
  (basic/unitcell lvs
  (as-> sc x
    (atom-pos x)
    (gneigh/neighbors x 0.1 1.8)
    (add-H x lvs)))))







;;;;;;;;;;;;;;;;;;;;;;;;; graphene supercells ;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-armchair-graphene-supercell
"This will make a graphene supercell that runs in the x-direction.
  The y-direction edges are armchair.

CC-dist is the bond length between neighboring C atoms,
nxcells is the number of unit cells in the direction of the ribbon,
nycells is the number of unit cells in the direction perpendicular to the ribbon,
rtl is a boolean that determins if the C atoms with the largest y-values are discarded"
[C1 C2 C3 C4 CC-dist nxcells nycells rtl]
(let [uclvs (armchair-lvs CC-dist)
         sc (if (true? rtl)
                (remove-highest-y (:mol (supercell (armchair-uc C1 C2 C3 C4 CC-dist)  uclvs nxcells nycells 1)))
                (:mol (supercell (armchair-uc C1 C2 C3 C4 CC-dist) uclvs nxcells nycells 1)))]
  (hash-map :lvs [(map * (first uclvs) (repeat nxcells))
    (map * (second uclvs) (repeat nycells))
    (last uclvs)]
            :mol (atom-pos sc))))





(defn make-zigzag-graphene-supercell
"This will make a graphene supercell that runs in the x-direction.
  The y-direction edges are zigzag.

  CC-dist is the bond length between neighboring C atoms,
nxcells is the number of unit cells in the direction of the ribbon,
nycells is the number of unit cells in the direction perpendicular to the ribbon,
rtl is a boolean that determins if the C atoms with the largest y-values are discarded"
[C1 C2 C3 C4 CC-dist nxcells nycells rtl]
(let [uclvs (zigzag-lvs CC-dist)
         sc (if (true? rtl)
               ((comp remove-highest-y remove-highest-y) (:mol (supercell (zigzag-uc C1 C2 C3 C4 CC-dist)  uclvs nxcells nycells 1)))
               (:mol (supercell (zigzag-uc C1 C2 C3 C4 CC-dist)  uclvs nxcells nycells 1)))]
  (hash-map :lvs [(map * (first uclvs) (repeat nxcells))
    (map * (second uclvs) (repeat nycells))
    (last uclvs)]
            :mol (atom-pos sc))))



;;;;;;;;;;;;;;;;;;;;;;;;; Other 2-D Allotropes ;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn gamma-graphyne
  "a is the C-C bond length in the rings
   b is the C-C bond length between a ring carbon and a carbon involved in the triple bonds
   c is the C-C bond length of the triple bonds
   n is the number of triple bonds between rings

(citation: J. Chern. Phys. 87 (11))

 Usage: (gamma-graphyne 'C 'C 1.42 1.4 1.22 1)"
  [species1 species2 a b c n]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        f1 #(if (odd? %) s1 s2)
        f2 #(if (odd? %) s2 s1)
        Auv [1 0 0]
        Buv [0.5 (* (cmat/sqrt 3.0) 0.5) 0]
        abc (+ a a (* (inc n) b) (* n c))
        A (* abc Auv)
        B (* abc Buv)
        U1 (+ A B)
        U2 (cmat/normalise (- B A))
        atm1 (+ (* 0.5 U1) (* a Buv))
        atm2 (+ (* 0.5 U1) (* a Auv))
        atm3 (+ (* 0.5 U1) (* -1 a U2))
        atm4 (+ (* 0.5 U1) (* -1 a Buv))
        atm5 (+ (* 0.5 U1) (* -1 a Auv))
        atm6 (+ (* 0.5 U1) (* a U2))
        f #(reduce + (take % (cycle [b c])))
        atm7 #(+ atm1 (* (f %) Buv))
        atm8 #(+ atm2 (* (f %) Auv))
        atm9 #(+ atm3 (* -1 (f %) U2))
        atm10 #(+ atm4 (* -1 (f %) Buv))
        atm11 #(+ atm5 (* -1 (f %) Auv))
        atm12 #(+ atm6 (* (f %) U2))]
    (hash-map :lvs [A B [0 0 30]]
     :mol (flatten [[(basic/new-atom s1 atm1 nil nil nil nil 1)
      (basic/new-atom s2 atm2 nil nil nil nil 2)
      (basic/new-atom s1 atm3 nil nil nil nil 3)
      (basic/new-atom s2 atm4 nil nil nil nil 4)
      (basic/new-atom s1 atm5 nil nil nil nil 5)
      (basic/new-atom s2 atm6 nil nil nil nil 6)]
     (map #(vector
            (basic/new-atom (f2 %) (atm7 %) nil nil nil nil 7)
            (basic/new-atom (f1 %) (atm8 %) nil nil nil nil 8)
            (basic/new-atom (f2 %) (atm9 %) nil nil nil nil 9)
            (basic/new-atom (f1 %) (atm10 %) nil nil nil nil 10)
            (basic/new-atom (f2 %) (atm11 %) nil nil nil nil 11)
            (basic/new-atom (f1 %) (atm12 %) nil nil nil nil 12)) (range 1 (inc n)))]))))








(defn beta-graphyne
  "a is the C-C bond length in the rings
   b is the C-C bond length between a ring carbon and a carbon involved in the triple bonds
   c is the C-C bond length of the triple bonds
(citation: J. Chern. Phys. 87 (11))

  Usage: (beta-graphyne 'C 'C 1.43 1.4 1.2)"
  [species1 species2 a b c]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        a1uv [1 0 0]
        a2uv [-0.5 (* (cmat/sqrt 3.0) 0.5) 0]
        abc (+ a a a a b c c)]
        (hash-map :lvs [[abc 0 0] (* abc a2uv) [0 30 0]]
        :mol [(basic/new-atom s1 [(+ a a c) 0 0] nil nil nil nil 0)
         (basic/new-atom s2 (+ [(+ a a c) 0 0] (* a a2uv)) nil nil nil nil 1)
         (basic/new-atom s1 (+ [(+ a a c) 0 0] (* (+ a c) a2uv)) nil nil nil nil 2)
         (basic/new-atom s2 (+ [(+ a a c) 0 0] (* (+ a a c) a2uv)) nil nil nil nil 3)
         (basic/new-atom s1 (+ [(+ a c) 0 0] (* (+ a a c) a2uv)) nil nil nil nil 4)
         (basic/new-atom s2 (+ [a 0 0] (* (+ a a c) a2uv)) nil nil nil nil 5)
         (basic/new-atom s1 (* (+ a a c) a2uv) nil nil nil nil 6) ;atm7
         (basic/new-atom s2 (+ [(- 0 a a c) 0 0] (* (+ a c) (cmat/normalise (+ a1uv a2uv)))) nil nil nil nil 7)
         (basic/new-atom s1 (+ [(- 0 a a c) 0 0] (* a (cmat/normalise (+ a1uv a2uv)))) nil nil nil nil 8)
         (basic/new-atom s2 [(- 0 a a c) 0 0] nil nil nil nil 9)
         (basic/new-atom s1 (+ [(- 0 a a c) 0 0] (* (- a) a2uv)) nil nil nil nil 10)
         (basic/new-atom s2 (+ [(- 0 a a c) 0 0] (* (- 0 a c) a2uv)) nil nil nil nil 11)
         (basic/new-atom s1 (+ [(- 0 a a c) 0 0] (* (- 0 a a c) a2uv)) nil nil nil nil 12) ;13
         (basic/new-atom s2 (+ [(- 0  a c) 0 0] (* (- 0 a a c) a2uv)) nil nil nil nil 13)
         (basic/new-atom s1 (+ [(- a) 0 0] (* (- 0 a a c) a2uv)) nil nil nil nil 14)
         (basic/new-atom s2 (* (- 0 a a c) a2uv) nil nil nil nil 15) ;16
         (basic/new-atom s1 (+ (map * (repeat (- 0 a a c)) a2uv) (* a (cmat/normalise (+ a1uv a2uv)))) nil nil nil nil 16)
         (basic/new-atom s2 (+ (* (- 0 a a c) a2uv) (* (+ a c) (cmat/normalise (+ a1uv a2uv)))) nil nil nil nil 17)])))





(defn alpha-graphyne
  "a is the C-C bond length in the rings
   b is the C-C bond length between a ring carbon and a carbon involved in the triple bonds
   c is the C-C bond length of the triple bonds
(citation: J. Chern. Phys. 87 (11))

  Usage: (alpha-graphyne Si 1.417055  1.4 1.244)"
  [species1 a b c]
  (let [s1 (.intern species1)
        a1uv [(* (cmat/sqrt 3.0) 0.5) -0.5 0]
        a2uv [(* (cmat/sqrt 3.0) 0.5)  0.5 0]
        abc (+ a a a a b c c)
        lvs [(* abc a1uv) (* abc a2uv) [0 30 0]]]
        (hash-map :lvs lvs
          [(basic/new-atom s1 (* (* 3 b) (cmat/normalise (+ a1uv a2uv))) nil nil nil nil 0)
           (basic/new-atom s1 (* (* 4 b) (cmat/normalise (+ a1uv a2uv))) nil nil nil nil 1)
           (basic/new-atom s1 (* (* 5 b) (cmat/normalise (+ a1uv a2uv))) nil nil nil nil 2)
           (basic/new-atom s1 (* (* 6 b) (cmat/normalise (+ a1uv a2uv))) nil nil nil nil 3)])))






(defn biphenylene-carbon
  "a is the C-C bond length in the rings
   b is the C-C bond length between a ring carbon and a carbon involved in the triple bonds
   c is the C-C bond length of the triple bonds

  The origin in this puc is set at the center of the four member ring, thus many of the
  points fall outside of the parallelepiped described by the lattice vectors if their vertex is at the origin.

  This is called graphenylene in The Open Organic Chemistry Journal, 2011, 5, (Suppl 1-M8) 117-126.

  Usage: (biphenylene-carbon 'B 'N 1.485 1.476 1.365)"
  [species1 species2 a b c]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        buv (a-two 1.0)
        auv (a-one 1.0)
        aa (+ b c (* 2 a (cmat/cos (/ ed/pi 6))) (* 2 c (cmat/cos (/ ed/pi 3))))
        apb (+ (* aa auv) (* aa buv))
        apb_mag (cmat/normalise  apb)
        apb_uv (cmat/normalise apb)
        bma (- (* aa buv) (* aa auv))
        bma_mag (cmat/normalise  bma)
        bma_uv (cmat/normalise bma)
        atm1 (+ (* (* 0.5 (- apb_mag a)) apb_uv)(* (* 0.5 (- b)) bma_uv))
        atm2 (+ (* (* 0.5 (+ apb_mag a)) apb_uv)(* (* 0.5 (- b)) bma_uv))
        atm3 (+ (* (* 0.5 (- apb_mag a)) apb_uv)(* (* 0.5 b) bma_uv))
        atm4 (+ (* (* 0.5 (+ apb_mag a)) apb_uv)(* (* 0.5 b) bma_uv))
        atm5 (+ atm4 (* c buv))
        atm8 (+ atm2 (* c auv))
        atm6 (+ atm5 (* b auv))
        atm7 (+ atm8 (* b buv))
        atm9 (- atm1 (* c buv))
        atm10 (- atm9 (* b auv))
        atm12 (- atm3 (* c auv))
        atm11 (- atm12 (* b buv))]
        (hash-map :lvs  [(* aa auv) (* aa buv) [0 0 30]]
          :mol [(basic/new-atom s2 atm1 nil nil nil nil 0)
           (basic/new-atom s1 atm2 nil nil nil nil 1)
           (basic/new-atom s1 atm3 nil nil nil nil 2)
           (basic/new-atom s2 atm4 nil nil nil nil 3)
           (basic/new-atom s1 atm5 nil nil nil nil 4)
           (basic/new-atom s2 atm6 nil nil nil nil 5)
           (basic/new-atom s1 atm7 nil nil nil nil 6)
           (basic/new-atom s2 atm8 nil nil nil nil 7)
           (basic/new-atom s1 atm9 nil nil nil nil 8)
           (basic/new-atom s2 atm10 nil nil nil nil 9)
           (basic/new-atom s1 atm11 nil nil nil nil 10)
           (basic/new-atom s2 atm12 nil nil nil nil 11)])))




(defn- four-six-carbophene-odd
  "help function.  There is probably a fix that means we don't
   need the helper functions but I don't have time."
  ;[species1 species2 species3 a b c d units]
  [species1 species2  a b c  units]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        buv (a-two 1.0)
        auv (a-one 1.0)
        f1 #(if (odd? %) s1 s2)
        f2 #(if (odd? %) s2 s1)
        aa (* (dec units) (+ b c (* 2 a (cmat/cos (/ ed/pi 6))) (* 2 c (cmat/cos (/ ed/pi 3)))))
        apb (+ (* aa auv) (* aa buv))
        apb_mag (cmat/normalise  apb)
        apb_uv (cmat/normalise apb)
        bma (- (* aa buv) (* aa auv))
        bma_mag (cmat/normalise  bma)
        bma_uv (cmat/normalise bma)
        atm1 (+ (* (* 0.5 (- apb_mag a)) apb_uv)(* (* 0.5 (- b)) bma_uv))
        atm2 (+ (* (* 0.5 (+ apb_mag a)) apb_uv)(* (* 0.5 (- b)) bma_uv))
        atm3 (+ (* (* 0.5 (- apb_mag a)) apb_uv)(* (* 0.5 b) bma_uv))
        atm4 (+ (* (* 0.5 (+ apb_mag a)) apb_uv)(* (* 0.5 b) bma_uv))
        atm5 (+ atm4 (* c buv))
        atm8 (+ atm2 (* c auv))
        atm6 (+ atm5 (* b auv))
        atm7 (+ atm8 (* b buv))
        atm9 (- atm1 (* c buv))
        atm10 (- atm9 (* b auv))
        atm12 (- atm3 (* c auv))
        atm11 (- atm12 (* b buv))
        up  (+ a (* 2 c (cmat/cos (/ ed/pi 6))))
        xup  (* 2 c (cmat/cos (/ ed/pi 6)))
        vxup   [0.5000000000000001 0.8660254037844386 0]
        vbxup   [-0.5000000000000001 0.8660254037844386 0]
        vxdown [0.5000000000000001 -0.8660254037844386 0]
        vbxdown [-0.5000000000000001 -0.8660254037844386 0]
        lvs [(* aa auv) (* aa buv) [0 0 30]]
        mol (atom-pos (flatten
          [(basic/new-atom s2 atm1 nil nil nil nil 0)
           (basic/new-atom s1 atm2 nil nil nil nil 1)
           (basic/new-atom s1 atm3 nil nil nil nil 2)
           (basic/new-atom s2 atm4 nil nil nil nil 3)
           (basic/new-atom s1 atm5 nil nil nil nil 4)
           (basic/new-atom s2 atm6 nil nil nil nil 5)
           (basic/new-atom s1 atm7 nil nil nil nil 6)
           (basic/new-atom s2 atm8 nil nil nil nil 7)
           (basic/new-atom s1 atm9 nil nil nil nil 8)
           (basic/new-atom s2 atm10 nil nil nil nil 9)
           (basic/new-atom s1 atm11 nil nil nil nil 10)
           (basic/new-atom s2 atm12 nil nil nil nil 11)
     (map #(vector
            (basic/new-atom (f2 %) (+ atm2 (* % [up 0 0])) nil nil nil nil 7)
            (basic/new-atom (f1 %) (+ atm4 (* % [up 0 0])) nil nil nil nil 8)
            (basic/new-atom (f2 %) (+ atm5 (* % [up 0 0])) nil nil nil nil 9)
            (basic/new-atom (f1 %) (+ atm6 (* % [up 0 0])) nil nil nil nil 10)
            (basic/new-atom (f2 %) (+ atm7 (* % [up 0 0])) nil nil nil nil 11)
            (basic/new-atom (f1 %) (+ atm8 (* % [up 0 0])) nil nil nil nil 12)) (range 1 (dec units)))
(map #(vector
       (basic/new-atom (f1 %) (+ atm2 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 7)
       (basic/new-atom (f2 %) (+ atm4 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 8)
       (basic/new-atom (f1 %) (+ atm5 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 9)
       (basic/new-atom (f2 %) (+ atm6 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 10)
       (basic/new-atom (f1 %) (+ atm7 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 11)
       (basic/new-atom (f2 %) (+ atm8 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 12)) (range 1 (int (Math/ceil (/ units 2)))))
(map #(vector
       (basic/new-atom (f1 %) (+ atm2 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 7) ;only when units is even
       (basic/new-atom (f2 %) (+ atm4 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 8) ;only when units is even
       (basic/new-atom (f1 %) (+ atm5 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 9) ;only when units is even
       (basic/new-atom (f2 %) (+ atm6 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 10) ;only when units is even
       (basic/new-atom (f1 %) (+ atm7 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 11) ;only when units is even
       (basic/new-atom (f2 %) (+ atm8 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 12)) (range 1 (int (Math/floor (/ units 2)))))
(map #(vector
       (basic/new-atom (f1 %) (+ atm1  (* % up vbxup)) nil nil nil nil 7)
       (basic/new-atom (f2 %) (+ atm3  (* % up vbxup)) nil nil nil nil 8)
       (basic/new-atom (f2 %) (+ atm9  (* % up vbxup)) nil nil nil nil 9)
       (basic/new-atom (f1 %) (+ atm10  (* % up vbxup)) nil nil nil nil 10)
       (basic/new-atom (f2 %) (+ atm11  (* % up vbxup)) nil nil nil nil 11)
       (basic/new-atom (f1 %) (+ atm12  (* % up vbxup)) nil nil nil nil 12)) (range 1 (int (Math/ceil (/ units 2)))))
(map #(vector
       (basic/new-atom (f1 %) (+ atm1  (* % up vbxdown)) nil nil nil nil 7)
       (basic/new-atom (f2 %) (+ atm3  (* % up vbxdown)) nil nil nil nil 8)
       (basic/new-atom (f2 %) (+ atm9  (* % up vbxdown)) nil nil nil nil 9)
       (basic/new-atom (f1 %) (+ atm10 (* % up vbxdown)) nil nil nil nil 10)
       (basic/new-atom (f2 %) (+ atm11 (* % up vbxdown)) nil nil nil nil 11)
       (basic/new-atom (f1 %) (+ atm12 (* % up vbxdown)) nil nil nil nil 12)) (range 1 (int (Math/floor (/ units 2)))))
]))]
        (hash-map :lvs lvs
                  :mol (gmol/shift (* 0.385 (+ (first lvs) (second lvs)))
                        (add-H mol lvs)))))



(defn- four-six-carbophene-even
  "help function.  There is probably a fix that means we don't
   need the helper functions but I don't have time."
  ;[species1 species2 species3 a b c d units]
  [species1 species2  a b c  units]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        buv (a-two 1.0)
        auv (a-one 1.0)
        f1 #(if (odd? %) s1 s2)
        f2 #(if (odd? %) s2 s1)
        aa (* (dec units) (+ b c (* 2 a (cmat/cos (/ ed/pi 6))) (* 2 c (cmat/cos (/ ed/pi 3)))))
        apb (+ (* aa auv) (* aa buv))
        apb_mag (cmat/normalise  apb)
        apb_uv (cmat/normalise apb)
        bma (- (* aa buv) (* aa auv))
        bma_mag (cmat/normalise  bma)
        bma_uv (cmat/normalise bma)
        atm1 (+ (* (* 0.5 (- apb_mag a)) apb_uv)(* (* 0.5 (- b)) bma_uv))
        atm2 (+ (* (* 0.5 (+ apb_mag a)) apb_uv)(* (* 0.5 (- b)) bma_uv))
        atm3 (+ (* (* 0.5 (- apb_mag a)) apb_uv)(* (* 0.5 b) bma_uv))
        atm4 (+ (* (* 0.5 (+ apb_mag a)) apb_uv)(* (* 0.5 b) bma_uv))
        atm5 (+ atm4 (* c buv))
        atm8 (+ atm2 (* c auv))
        atm6 (+ atm5 (* b auv))
        atm7 (+ atm8 (* b buv))
        atm9 (- atm1 (* c buv))
        atm10 (- atm9 (* b auv))
        atm12 (- atm3 (* c auv))
        atm11 (- atm12 (* b buv))
        up  (+ a (* 2 c (cmat/cos (/ ed/pi 6))))
        xup  (* 2 c (cmat/cos (/ ed/pi 6)))
        vxup   [0.5000000000000001 0.8660254037844386 0]
        vbxup   [-0.5000000000000001 0.8660254037844386 0]
        vxdown [0.5000000000000001 -0.8660254037844386 0]
        vbxdown [-0.5000000000000001 -0.8660254037844386 0]
        lvs [(* aa auv) (* aa buv) [0 0 30]]
        mol (atom-pos (flatten
          [(basic/new-atom s2 atm1 nil nil nil nil 0)
           (basic/new-atom s1 atm2 nil nil nil nil 1)
           (basic/new-atom s1 atm3 nil nil nil nil 2)
           (basic/new-atom s2 atm4 nil nil nil nil 3)
           (basic/new-atom s1 atm5 nil nil nil nil 4)
           (basic/new-atom s2 atm6 nil nil nil nil 5)
           (basic/new-atom s1 atm7 nil nil nil nil 6)
           (basic/new-atom s2 atm8 nil nil nil nil 7)
           (basic/new-atom s1 atm9 nil nil nil nil 8)
           (basic/new-atom s2 atm10 nil nil nil nil 9)
           (basic/new-atom s1 atm11 nil nil nil nil 10)
           (basic/new-atom s2 atm12 nil nil nil nil 11)
     (map #(vector
            (basic/new-atom (f2 %) (+ atm2 (* % [up 0 0])) nil nil nil nil 7)
            (basic/new-atom (f1 %) (+ atm4 (* % [up 0 0])) nil nil nil nil 8)
            (basic/new-atom (f2 %) (+ atm5 (* % [up 0 0])) nil nil nil nil 9)
            (basic/new-atom (f1 %) (+ atm6 (* % [up 0 0])) nil nil nil nil 10)
            (basic/new-atom (f2 %) (+ atm7 (* % [up 0 0])) nil nil nil nil 11)
            (basic/new-atom (f1 %) (+ atm8 (* % [up 0 0])) nil nil nil nil 12)) (range 1 (dec units)))
(map #(vector
       (basic/new-atom (f2 %) (+ atm2 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 7)
       (basic/new-atom (f1 %) (+ atm4 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 8)
       (basic/new-atom (f2 %) (+ atm5 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 9)
       (basic/new-atom (f1 %) (+ atm6 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 10)
       (basic/new-atom (f2 %) (+ atm7 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 11)
       (basic/new-atom (f1 %) (+ atm8 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxup)) nil nil nil nil 12)) (range 1 (int (Math/ceil (/ units 2)))))
(map #(vector
       (basic/new-atom (f2 %) (+ atm2 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 7) ;only when units is even
       (basic/new-atom (f1 %) (+ atm4 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 8) ;only when units is even
       (basic/new-atom (f2 %) (+ atm5 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 9) ;only when units is even
       (basic/new-atom (f1 %) (+ atm6 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 10) ;only when units is even
       (basic/new-atom (f2 %) (+ atm7 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 11) ;only when units is even
       (basic/new-atom (f1 %) (+ atm8 (* (- units 3) [up 0 0]) [up 0 0] (* % up vxdown)) nil nil nil nil 12)) (range 1 (int (Math/floor (/ units 2)))))
(map #(vector
       (basic/new-atom (f1 %) (+ atm1  (* % up vbxup)) nil nil nil nil 7)
       (basic/new-atom (f2 %) (+ atm3  (* % up vbxup)) nil nil nil nil 8)
       (basic/new-atom (f2 %) (+ atm9  (* % up vbxup)) nil nil nil nil 9)
       (basic/new-atom (f1 %) (+ atm10  (* % up vbxup)) nil nil nil nil 10)
       (basic/new-atom (f2 %) (+ atm11  (* % up vbxup)) nil nil nil nil 11)
       (basic/new-atom (f1 %) (+ atm12  (* % up vbxup)) nil nil nil nil 12)) (range 1 (int (Math/ceil (/ units 2)))))
(map #(vector
       (basic/new-atom (f1 %) (+ atm1  (* % up vbxdown)) nil nil nil nil 7)
       (basic/new-atom (f2 %) (+ atm3  (* % up vbxdown)) nil nil nil nil 8)
       (basic/new-atom (f2 %) (+ atm9  (* % up vbxdown)) nil nil nil nil 9)
       (basic/new-atom (f1 %) (+ atm10 (* % up vbxdown)) nil nil nil nil 10)
       (basic/new-atom (f2 %) (+ atm11 (* % up vbxdown)) nil nil nil nil 11)
       (basic/new-atom (f1 %) (+ atm12 (* % up vbxdown)) nil nil nil nil 12)) (range 1 (int (Math/floor (/ units 2)))))
]))]
        (hash-map :lvs lvs
                  :mol (gmol/shift (* 0.385 (+ (first lvs) (second lvs)))
                        (add-H mol lvs)))))



(defn four-six-carbophene
  "a is the C-C bond length in the rings
   b is the C-C bond length between a ring carbon and a carbon involved in the triple bonds
   c is the C-C bond length of the triple bonds

  The origin in this puc is set at the center of the four member ring, thus many of the
  points fall outside of the parallelepiped described by the lattice vectors if their vertex is at the origin.

  This is called graphenylene in The Open Organic Chemistry Journal, 2011, 5, (Suppl 1-M8) 117-126.

  Usage: (four-six-carbophene 'B 'N  1.485 1.476 1.365  3)"
  ;[species1 species2 species3 a b c d units]
  [species1 species2  a b c  units]
  (if (odd? units) (four-six-carbophene-odd species1 species2 a b c units)
                   (four-six-carbophene-even species1 species2 a b c units)))



(defn junkermeier
  "a is the graphene lattice constant

  Usage: (junkermeier 'B 'N 1.485 'H 1.09)"
  [species1 species2 a species3 b]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        s3 (.intern species3)
        buv (a-two 1.0)
        auv (a-one 1.0)
        A (* 3 a auv)
        B (* 3 a buv)
        apb (+ A B)
        apb_uv (cmat/normalise apb)
        bma (- B A)
        bma_uv (cmat/normalise bma)
        hx (* b (cmat/cos (/ ed/pi 3)))
        hy (* b (cmat/sin (/ ed/pi 3)))
        h1 (* b (cmat/normalise (+ (* 0.45 apb_uv) (* -0.45 bma_uv) (* 0.1 [0 0 1]))))
        h2 (* b (cmat/normalise (+ (* 0.45 apb_uv) (* 0.45 bma_uv) (* -0.1 [0 0 1]))))
        h3 (* b (cmat/normalise (+ (* -0.45 apb_uv) (* -0.45 bma_uv) (* -0.1 [0 0 1]))))
        h4 (* b (cmat/normalise (+ (* -0.45 apb_uv) (* 0.45 bma_uv) (* 0.1 [0 0 1]))))
        ]
        (hash-map :lvs  [A B [0 0 30]]
           :mol [(basic/new-atom s1 (- (* 0.5 apb) (* 0.5 a apb_uv)) nil nil nil nil 0)
            (basic/new-atom s2 (+ (* 0.5 apb) (* 0.5 a apb_uv)) nil nil nil nil 1)
            (basic/new-atom s2 (+ (* 0.5 A) (* 0.5 a apb_uv)) nil nil nil nil 2)
            ;(basic/new-atom s3 (+ (* 0.5 A) (* 0.5 a apb_uv) (* hx apb_uv) (* -1 hy bma_uv)) nil nil nil nil nil)
            (basic/new-atom s3 (+ (* 0.5 A) (* 0.5 a apb_uv) h1) nil nil nil nil 3)

            (basic/new-atom s2 (+ (* 0.5 B) (* 0.5 a apb_uv)) nil nil nil nil 4)
            (basic/new-atom s3 (+ (* 0.5 B) (* 0.5 a apb_uv) h2) nil nil nil nil 5)

            (basic/new-atom s1 (+ A (* 0.5 B) (* -0.5 a apb_uv)) nil nil nil nil 6)
            (basic/new-atom s3 (+ A (* 0.5 B) (* -0.5 a apb_uv) h3) nil nil nil nil 7)

            (basic/new-atom s1 (+ B (* 0.5 A) (* -0.5 a apb_uv)) nil nil nil nil 8)
            (basic/new-atom s3 (+ B (* 0.5 A) (* -0.5 a apb_uv) h4) nil nil nil nil 9)
            ])))



;(def sclbcaa (shift (* -0.5 (+ ((comp first :lvs) sclbc) ((comp second :lvs) sclbc)))  (:mol sclbc))


(defn cyclic-3naphthylene
  "a is the C-C distance

  The idea for this comes from: The Open Organic Chemistry Journal, 2011, 5, (Suppl 1-M8) 117-126
https://benthamopen.com/ABSTRACT/TOOCJ-5-117

  Usage: (cyclic-3naphthylene 'B 'N 1.485)"
  [species1 species2 a]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        buv (a-two 1.0)
        auv (a-one 1.0)
        csn (cmat/cos (/ ed/pi 6))
        sn  (cmat/sin (/ ed/pi 6))
        A (* (+ (* 6 csn a) a) auv);fix
        B (* (+ (* 6 csn a) a) buv);fix
        U1 (+ A B)
        U1_uv (cmat/normalise U1)
        U2 (- A B)
        U2_uv (cmat/normalise U2)]
        (hash-map :lvs  [A B [0 0 30]]
           :mol [(basic/new-atom s1 (+ (* 0.5 U1) (* 0.5 a U1_uv) (* 0.5 a U2_uv)) nil nil nil nil 0)
            (basic/new-atom s2 (+ (* 0.5 U1) (* 0.5 a U1_uv) (* -0.5 a U2_uv)) nil nil nil nil 1)
            (basic/new-atom s2 (+ (* 0.5 U1) (* -0.5 a U1_uv) (* 0.5 a U2_uv)) nil nil nil nil 2)
            (basic/new-atom s1 (+ (* 0.5 U1) (* -0.5 a U1_uv) (* -0.5 a U2_uv)) nil nil nil nil 3)

            (basic/new-atom s2 (+ (* 0.5 U1) (* 0.5 a U1_uv) (* (+ 0.5 csn csn) a U2_uv)) nil nil nil nil 4)
            (basic/new-atom s1 (+ (* 0.5 U1) (* -0.5 a U1_uv) (* (+ 0.5 csn csn) a U2_uv)) nil nil nil nil 5)
            (basic/new-atom s1 (+ (* 0.5 U1) (* a (+ 0.5 sn) U1_uv) (* (+ 0.5 csn) a U2_uv)) nil nil nil nil 6)
            (basic/new-atom s2 (+ (* 0.5 U1) (* -1 a (+ 0.5 sn) U1_uv) (* (+ 0.5 csn) a U2_uv)) nil nil nil nil 7)

            (basic/new-atom s2 (+ (* 0.5 U1) (* 0.5 a U1_uv) (* -1 (+ 0.5 csn csn) a U2_uv)) nil nil nil nil 8)
            (basic/new-atom s1 (+ (* 0.5 U1) (* -0.5 a U1_uv) (* -1 (+ 0.5 csn csn) a U2_uv)) nil nil nil nil 9)
            (basic/new-atom s1 (+ (* 0.5 U1) (* a (+ 0.5 sn) U1_uv) (* -1 (+ 0.5 csn) a U2_uv)) nil nil nil nil 10)
            (basic/new-atom s2 (+ (* 0.5 U1) (* -1 a (+ 0.5 sn) U1_uv) (* -1 (+ 0.5 csn) a U2_uv)) nil nil nil nil 11)

            (basic/new-atom s2 (* a U1_uv)  nil nil nil nil 12)
            (basic/new-atom s1 (* 2 a U1_uv) nil nil nil nil 13)
            (basic/new-atom s1 (+ (* a U1_uv) (* 2 csn a auv)) nil nil nil nil 14)
            (basic/new-atom s2 (+ (* a U1_uv) (* 2 csn a auv) (* a auv)) nil nil nil nil 15)
            (basic/new-atom s1 (+ (* a U1_uv) (* 2 csn a buv)) nil nil nil nil 16)
            (basic/new-atom s2 (+ (* a U1_uv) (* 2 csn a buv) (* a buv)) nil nil nil nil 17)

            (basic/new-atom s2 (+ U1 (* -1 a U1_uv))  nil nil nil nil 18)
            (basic/new-atom s1 (+ U1 (* -2 a U1_uv)) nil nil nil nil 19)
            (basic/new-atom s1 (+ U1 (* -1 a U1_uv) (* -2 csn a auv)) nil nil nil nil 20)
            (basic/new-atom s2 (+ U1 (* -1 a U1_uv) (* -2 csn a auv) (* -1 a auv)) nil nil nil nil 21)
            (basic/new-atom s1 (+ U1 (* -1 a U1_uv) (* -2 csn a buv)) nil nil nil nil 22)
            (basic/new-atom s2 (+ U1 (* -1 a U1_uv) (* -2 csn a buv) (* -1 a buv)) nil nil nil nil 23)])))






(defn foureighteight-archimedean-tile
  "Found in The Open Organic Chemistry Journal, 2011, 5, (Suppl 1-M8) 117-126.
 Usage: (foureighteight-archimedean-tile 'C 1.485 1.476)

  THIS DOESN'T HAVE THE CORRECT LVS STRUCTURE FOR MAKING TUBES."
  [species1 a b]
  (let [s1 (.intern species1)
        cp4 (cmat/cos (/ ed/pi 4))
        A [(+ b (* 2 a cp4)) 0 0]
        B [0 (+ b (* 2 a cp4)) 0]]
    (hash-map :lvs [A B [0 0 30]]
     :mol [(basic/new-atom s1 (+ (* 0.5 A) (* 0.5 b (cmat/normalise B))) nil nil nil nil 0)
      (basic/new-atom s1 (+ (* 0.5 B) (* (+ (* 0.5 b) (* 2 a cp4)) (cmat/normalise A))) nil nil nil nil 1)
      (basic/new-atom s1 (+ (* 0.5 A) (* (+ (* 0.5 b) (* 2 a cp4)) (cmat/normalise B))) nil nil nil nil 2)
      (basic/new-atom s1 (+ (* 0.5 b (cmat/normalise A)) (* 0.5 B)) nil nil nil nil 3)])))




(defn kagome
  "a is the C-C bond length in the rings

  I don't know if this one will work.

  Usage: (kagome species1 1.42 1.4)"
  [species1 c d]
  (let [s1 (.intern species1)
        Buv [1 0 0]
        Auv [0.5 (* (cmat/sqrt 3.0) 0.5) 0]
        aa (+ c c c c d d)
        A (* aa Auv)
        B (* aa Buv)
        ApB (+ A B)
        AmB (cmat/normalise (- A B))]
    (hash-map :lvs [A B [0 0 30]]
     :mol [(basic/new-atom s1 (* 0.5 ApB) nil nil nil nil 0)
      (basic/new-atom s1 (- (* 0.5 ApB) (* c (cmat/normalise A))) nil nil nil nil 1)
      (basic/new-atom s1 (- (* 0.5 ApB) (* (+ c d) (cmat/normalise A))) nil nil nil nil 2)
      (basic/new-atom s1 (* 0.5 B) nil nil nil nil 3)
      (basic/new-atom s1 (- (* 0.5 B) (* c AmB)) nil nil nil nil 4)
      (basic/new-atom s1 (- (* 0.5 B) (* (+ c d) AmB)) nil nil nil nil 5)
      (basic/new-atom s1 (* 0.5 A) nil nil nil nil 6)
      (basic/new-atom s1 (+ (* 0.5 A) (* c (cmat/normalise B))) nil nil nil nil 7)
      (basic/new-atom s1 (+ (* 0.5 A) (* (+ c d) (cmat/normalise B))) nil nil nil nil 8)
      (basic/new-atom s1 (+ (* 0.5 ApB) (* c (cmat/normalise B))) nil nil nil nil 9)
      (basic/new-atom s1 (+ (* 0.5 ApB) (* (+ c d) (cmat/normalise B))) nil nil nil nil 10)
      (basic/new-atom s1 (+ (* 0.5 ApB) (* c (cmat/normalise A))) nil nil nil nil 11)
      (basic/new-atom s1 (+ (* 0.5 ApB) (* (+ c d) (cmat/normalise A))) nil nil nil nil 12)
      (basic/new-atom s1 (+ (* 0.5 A) B (* -1 c AmB)) nil nil nil nil 13)
      (basic/new-atom s1 (+ (* 0.5 A) B (* -1 (+ c d) AmB)) nil nil nil nil 14)])))






(defn porous-graphene
  "a is the C-C bond length in the rings
  b is the C-H bond length.

  I don't know if this one will work.

  Usage: (porous-graphene 'C 'C 'H 1.421 1.421 1.09 1.09)
  Usage: (porous-graphene 'B 'N 'H 1.506590 1.460214 1.219986 1.047633)"
  [species1 species2 species3 a b h1 h2]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        s3 (.intern species3)
        buv (a-two 1.0)
        auv (a-one 1.0)
        b-one (cmat/cross [0 0 1] auv)
        b-two (cmat/cross [0 0 1] buv)
        a3  (a-three 1.0)
        A (+ (* 2 (a-one (* (cmat/sqrt 3) b)))  (a-one (* (cmat/sqrt 3) a)))
        B (+ (* 2 (a-two (* (cmat/sqrt 3) b)))  (a-two (* (cmat/sqrt 3) a)))
        apb (+ A B)
        apb_uv (cmat/normalise apb)
        p1 (+ (* 0.5 apb) (* 0.5 a apb_uv))
        p2 (+ (* 0.5 B) (* 0.5 a b-two) (* b apb_uv))
        p3 (+ (* 0.5 B) (* (+ (* 0.5 a) b b) b-two) (* -1 b apb_uv))
        c  (cmat/normalise (- p3 p2))
        d  (cmat/normalise (- (+ p1 (* -1 a c) )  (+ p1 (* a c) (* 2 a apb_uv))))
        ee (cmat/normalise  (- p2 (+ p1 (* -1 a c) (* a apb_uv))))
        ff (cmat/normalise (- (+ p1 (* a c) (* a apb_uv))   (+ p1 (* -1 a c) (* a apb_uv)) ))
        ]
        (hash-map :lvs  [A B [0 0 30]]
                :mol [(basic/new-atom s1 (- (* 0.5 apb) (* 0.5 a apb_uv)) nil nil nil nil 0)
                 (basic/new-atom s2 p1 nil nil nil nil 1)
                 (basic/new-atom s1 (+ (* 0.5 apb) (* (+ (* 0.5 a) b b) apb_uv)) nil nil nil nil 2)
                 (basic/new-atom s2 (- (* 0.5 apb) (* (+ (* 0.5 a) b b) apb_uv)) nil nil nil nil 3)
                 (basic/new-atom s3 (+ (* 0.5 apb) (* (+ (* 0.5 a) b b h2) apb_uv)) nil nil nil nil 4)
                 (basic/new-atom s3 (- (* 0.5 apb) (* (+ (* 0.5 a) b b h1) apb_uv)) nil nil nil nil 5)
                 (basic/new-atom s1 (+ (* 0.5 B) (* 0.5 a b-two)) nil nil nil nil 6)
                 (basic/new-atom s2 (+ (* 0.5 B) (* (+ (* 0.5 a) b b) b-two)) nil nil nil nil 7)
                 (basic/new-atom s3 (+ (* 0.5 B) (* (+ (* 0.5 a) b b h1) b-two)) nil nil nil nil 8)
                 (basic/new-atom s2 p2 nil nil nil nil 9)
                 (basic/new-atom s1 p3 nil nil nil nil 10)
                 (basic/new-atom s3 (- p2 (* h1 c)) nil nil nil nil 11)
                 (basic/new-atom s1 (+ p1 (* b c) (* b apb_uv)) nil nil nil nil 12)
                 (basic/new-atom s2 (+ p1 (* b c) (* 2 b apb_uv)) nil nil nil nil 13)
                 (basic/new-atom s2 (+ p1 (* -1 b c) (* b apb_uv)) nil nil nil nil 14)
                 (basic/new-atom s1 (+ p1 (* -1 b c)) nil nil nil nil 15)
                 (basic/new-atom s3 (+ p1 (* -1 b c) (* h2 d)) nil nil nil nil 16)
                 (basic/new-atom s3 (+ (+ p1 (* b c) (* b apb_uv)) (* h2 ff)) nil nil nil nil 17)])))





(defn triazine-based-graphitic-carbon-nitride
  "From DOI: 10.1002/anie.201402191
  a is the C-N bond length
  Usage: (triazine-based-graphitic-carbon-nitride 'N 'C 1.455)"
  [species1 species2 a]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        A (* 4 a (cmat/cos (/ ed/pi 6)) (a-one 1.0))
        B (* 4 a (cmat/cos (/ ed/pi 6)) (a-two 1.0))
        apb (+ A B)
        apb_uv (cmat/normalise apb)]
        (hash-map :lvs  [A B [0 0 30]]
            :mol [(basic/new-atom s1 [0 0 0] nil nil nil nil 0)
                 (basic/new-atom s1 (* 0.5 A) nil nil nil nil 1)
                 (basic/new-atom s1 A nil nil nil nil 2)
                 (basic/new-atom s1 (* 0.5 B) nil nil nil nil 3)
                 (basic/new-atom s1 B nil nil nil nil 4)
                 (basic/new-atom s1 (* 0.5 apb) nil nil nil nil 5)
                 (basic/new-atom s2 (* a apb_uv) nil nil nil nil 6)
                 (basic/new-atom s2 (+ (* 0.5 A) (* a apb_uv)) nil nil nil nil 7)
                 (basic/new-atom s2 (+ (* 0.5 B) (* a apb_uv)) nil nil nil nil 8)])))



(defn Octafunctionalized-Biphenylenes
  "Precurser synthesized in: 'Octafunctionalized Biphenylenes: Molecular Precursors for Isomeric Graphene Nanostructures'
  -by Florian Schlütter, Tomohiko Nishiuchi, Volker Enkelmann, and Klaus Müllen
   DOI: 10.1002/ange.2013093246767

  Usage (Octafunctionalized-Biphenylenes 'C 1.537659 1.426094 1.409636  1.469558) "
  [species1 a b c d]
  (let [aprime d
        bprime c
        cprime b
        dprime a
        s1 (.intern species1)
        A [(+ dprime (* 2 bprime (cmat/cos (/ ed/pi 6)))) 0 0]
        B [0 (+ aprime cprime (* 2 bprime (cmat/sin (/ ed/pi 6)))) 0]
        atm1 (+ (* 0.5 A) (* 0.5 aprime [0 1 0]))
        atm2 (+ (* 0.5 dprime [1 0 0]) (* 0.5 B) (* -0.5 cprime [0 1 0]))
        atm3 (+ atm2 (* 2 bprime (cmat/cos (/ ed/pi 6)) [1 0 0]))
        atm4 (+ atm2 (* cprime [0 1 0]))
        atm5 (+ atm3 (* cprime [0 1 0]))
        atm6 (+ atm1 (* (+ cprime (* 2 bprime (cmat/sin (/ ed/pi 6)))) [0 1 0]))]
    (hash-map :lvs  [A B [0 0 30]]
            :mol [(basic/new-atom s1 atm1 nil nil nil nil 1)
                  (basic/new-atom s1 atm2 nil nil nil nil 2)
                  (basic/new-atom s1 atm3 nil nil nil nil 3)
                  (basic/new-atom s1 atm4 nil nil nil nil 4)
                  (basic/new-atom s1 atm5 nil nil nil nil 5)
                  (basic/new-atom s1 atm6 nil nil nil nil 6)])))






(defn Inorganic-Octafunctionalized-Biphenylenes
  "Precurser synthesized in: 'Octafunctionalized Biphenylenes: Molecular Precursors for Isomeric Graphene Nanostructures'
  -by Florian Schlütter, Tomohiko Nishiuchi, Volker Enkelmann, and Klaus Müllen
   DOI: 10.1002/ange.2013093246767

  Usage (Octafunctionalized-Biphenylenes 'B 'N 1.4995052857142857  1.465623 1.49851875  1.51432575) "
  [species1 species2 a b c d]
  (let [aprime d
        bprime c
        cprime b
        dprime a
        s1 (.intern species1)
        s2 (.intern species2)
        A [(+ dprime (* 2 bprime (cmat/cos (/ ed/pi 6)))) 0 0]
        B [0 (+ aprime cprime (* 2 bprime (cmat/sin (/ ed/pi 6)))) 0]
        atm1 (+ (* 0.5 A) (* 0.5 aprime [0 1 0]))
        atm2 (+ (* 0.5 dprime [1 0 0]) (* 0.5 B) (* -0.5 cprime [0 1 0]))
        atm3 (+ atm2 (* 2 bprime (cmat/cos (/ ed/pi 6)) [1 0 0]))
        atm4 (+ atm2 (* cprime [0 1 0]))
        atm5 (+ atm3 (* cprime [0 1 0]))
        atm6 (+ atm1 (* (+ cprime (* 2 bprime (cmat/sin (/ ed/pi 6)))) [0 1 0]))]
    (hash-map :lvs  [(* 2.0 A) B [0 0 30]]
            :mol [(basic/new-atom s1 atm1 nil nil nil nil 1)
                  (basic/new-atom s2 atm2 nil nil nil nil 2)
                  (basic/new-atom s2 atm3 nil nil nil nil 3)
                  (basic/new-atom s1 atm4 nil nil nil nil 4)
                  (basic/new-atom s1 atm5 nil nil nil nil 5)
                  (basic/new-atom s2 atm6 nil nil nil nil 6)
                  (basic/new-atom s2 (+ A atm1) nil nil nil nil 7)
                  (basic/new-atom s1 (+ A atm2) nil nil nil nil 8)
                  (basic/new-atom s1 (+ A atm3) nil nil nil nil 9)
                  (basic/new-atom s2 (+ A atm4) nil nil nil nil 10)
                  (basic/new-atom s2 (+ A atm5) nil nil nil nil 11)
                  (basic/new-atom s1 (+ A atm6) nil nil nil nil 12)])))






(defn Octafunctionalized-Biphenylenes-type2-rectagularsc
  "Precurser synthesized in: 'Octafunctionalized Biphenylenes: Molecular Precursors for Isomeric Graphene Nanostructures'
  -by Florian Schlütter, Tomohiko Nishiuchi, Volker Enkelmann, and Klaus Müllen
   DOI: 10.1002/ange.201309324

  Usage: (Octafunctionalized-Biphenylenes-type2 'C 'C 1.539807 1.415841  1.394072  1.450921 1.436087 1.441179 1.443518)"
  [species1 species2 a b c d e f g]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        cs (cmat/cos (/ ed/pi 3))
        sn (cmat/sin (/ ed/pi 3))
        A [(+ b (* 2 c cs) d) 0.0 0.0]
        B [0.0 (+ a a (* 2 (+ c c e e g) sn)) 0.0]
        atm1 [(* 0.5 b) (* 0.5 a) 0.0]
        atm13 (+ A [(* -0.5 b) (* 0.5 a) 0.0])
        atm2 (+ atm1 [(* c cs) (* c sn) 0.0])
        atm14 (+ atm2 [d 0.0 0.0])
        atm3 (+ atm2 [(* -1.0 e cs) (* e sn) 0.0])
        atm15 (+ atm14 [(* e cs) (* e sn) 0.0])
        atm4 (+ atm3 [(* g cs) (* g sn) 0.0])
        atm16 (+ atm4 [f 0.0 0.0])
        atm5 (+ atm4 [(* -1.0 e cs) (* e sn) 0.0])
        atm17 (+ atm16 [(* e cs) (* e sn) 0.0])
        atm6 (+ atm5 [(* c cs) (* c sn) 0.0])
        atm18 (+ atm6 [b 0.0 0.0])
        atm7 (+ atm6 [0.0 a 0.0])
        atm19 (+ atm7 [b 0.0  0.0])
        atm8 (+ atm7 [(* -1.0 c cs) (* c sn) 0.0])
        atm9 (+ atm8 [(* e cs) (* e sn) 0.0])
        atm10 (+ atm9 [(* -1.0 g cs) (* g sn) 0.0])
        atm21 (+ atm9 [f 0.0 0.0])
        atm22 (+ atm21 [(* g cs) (* g sn) 0.0])
        atm20 (+ atm19 [(* c cs) (* c sn) 0.0])
        atm11 (+ atm10 [(* e cs) (* e sn) 0.0])
        atm23 (+ atm11 [d 0.0 0.0])
        atm24 (+ atm23 [(* c cs) (* c sn) 0.0])
        atm12 (+ atm11 [(* -1.0 c cs) (* c sn) 0.0])]
     (hash-map :lvs  [A B [0 0 30]]
            :mol [(basic/new-atom s1 atm1 nil nil nil nil 1)
                  (basic/new-atom s2 atm2 nil nil nil nil 2)
                  (basic/new-atom s1 atm3 nil nil nil nil 3)
                  (basic/new-atom s2 atm4 nil nil nil nil 4)
                  (basic/new-atom s1 atm5 nil nil nil nil 5)
                  (basic/new-atom s2 atm6 nil nil nil nil 6)
                  (basic/new-atom s1 atm7 nil nil nil nil 7)
                  (basic/new-atom s2 atm8 nil nil nil nil 8)
                  (basic/new-atom s1 atm9 nil nil nil nil 9)
                  (basic/new-atom s2 atm10 nil nil nil nil 10)
                  (basic/new-atom s1 atm11 nil nil nil nil 11)
                  (basic/new-atom s2 atm12 nil nil nil nil 12)
                  (basic/new-atom s2 atm13 nil nil nil nil 13)
                  (basic/new-atom s1 atm14 nil nil nil nil 14)
                  (basic/new-atom s2 atm15 nil nil nil nil 15)
                  (basic/new-atom s1 atm16 nil nil nil nil 16)
                  (basic/new-atom s2 atm17 nil nil nil nil 17)
                  (basic/new-atom s1 atm18 nil nil nil nil 18)
                  (basic/new-atom s2 atm19 nil nil nil nil 19)
                  (basic/new-atom s1 atm20 nil nil nil nil 20)
                  (basic/new-atom s2 atm21 nil nil nil nil 21)
                  (basic/new-atom s1 atm22 nil nil nil nil 22)
                  (basic/new-atom s2 atm23 nil nil nil nil 23)
                  (basic/new-atom s1 atm24 nil nil nil nil 24)])))





(defn Octafunctionalized-Biphenylenes-type2
  "Precurser synthesized in: 'Octafunctionalized Biphenylenes: Molecular Precursors for Isomeric Graphene Nanostructures'
  -by Florian Schlütter, Tomohiko Nishiuchi, Volker Enkelmann, and Klaus Müllen
   DOI: 10.1002/ange.201309324

  Usage: (Octafunctionalized-Biphenylenes-type2 'C 'C 1.539807 1.415841  1.394072  1.450921 1.436087 1.441179 1.443518)"
  [species1 species2 a b c d e f g]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        cs (cmat/cos (/ ed/pi 3))
        sn (cmat/sin (/ ed/pi 3))
        A [(+ b (* 2 c cs) d) 0.0 0.0]
        B [0.0 (+ a a (* 2 (+ c c e e g) sn)) 0.0]
        Aprime [(* 0.5 (+ b (* 2 c cs) d)) (* -0.5 (+ a a (* 2 (+ c c e e g) sn))) 0.0]
        Bprime [(* 0.5 (+ b (* 2 c cs) d)) (* 0.5 (+ a a (* 2 (+ c c e e g) sn))) 0.0]
        atm1 [(* 0.5 b) (* 0.5 a) 0.0]
        atm13 (+ A [(* -0.5 b) (* 0.5 a) 0.0])
        atm2 (+ atm1 [(* c cs) (* c sn) 0.0])
        atm14 (+ atm2 [d 0.0 0.0])
        atm3 (+ atm2 [(* -1.0 e cs) (* e sn) 0.0])
        atm15 (+ atm14 [(* e cs) (* e sn) 0.0])
        atm4 (+ atm3 [(* g cs) (* g sn) 0.0])
        atm16 (+ atm4 [f 0.0 0.0])
        atm5 (+ atm4 [(* -1.0 e cs) (* e sn) 0.0])
        atm17 (+ atm16 [(* e cs) (* e sn) 0.0])
        atm6 (+ atm5 [(* c cs) (* c sn) 0.0])
        atm18 (+ atm6 [b 0.0 0.0])
        atm7 (+ atm6 [0.0 a 0.0])
        atm19 (+ atm7 [b 0.0  0.0])
        atm8 (+ atm7 [(* -1.0 c cs) (* c sn) 0.0])
        atm9 (+ atm8 [(* e cs) (* e sn) 0.0])
        atm10 (+ atm9 [(* -1.0 g cs) (* g sn) 0.0])
        atm21 (+ atm9 [f 0.0 0.0])
        atm22 (+ atm21 [(* g cs) (* g sn) 0.0])
        atm20 (+ atm19 [(* c cs) (* c sn) 0.0])
        atm11 (+ atm10 [(* e cs) (* e sn) 0.0])
        atm23 (+ atm11 [d 0.0 0.0])
        atm24 (+ atm23 [(* c cs) (* c sn) 0.0])
        atm12 (+ atm11 [(* -1.0 c cs) (* c sn) 0.0])]
     (hash-map :lvs  [Aprime Bprime [0 0 30]]
            :mol [;(basic/new-atom s1 atm1 nil nil nil nil 1)
                  ;(basic/new-atom s2 atm2 nil nil nil nil 2)
                  ;(basic/new-atom s1 atm3 nil nil nil nil 3)
                  (basic/new-atom s2 atm4 nil nil nil nil 4)
                  (basic/new-atom s1 atm5 nil nil nil nil 5)
                  (basic/new-atom s2 atm6 nil nil nil nil 6)
                  (basic/new-atom s1 atm7 nil nil nil nil 7)
                  (basic/new-atom s2 atm8 nil nil nil nil 8)
                  (basic/new-atom s1 atm9 nil nil nil nil 9)
                  ;(basic/new-atom s2 atm10 nil nil nil nil 10)
                  ;(basic/new-atom s1 atm11 nil nil nil nil 11)
                  ;(basic/new-atom s2 atm12 nil nil nil nil 12)
                  ;(basic/new-atom s2 atm13 nil nil nil nil 13)
                  ;(basic/new-atom s1 atm14 nil nil nil nil 14)
                  ;(basic/new-atom s2 atm15 nil nil nil nil 15)
                  (basic/new-atom s1 atm16 nil nil nil nil 16)
                  (basic/new-atom s2 atm17 nil nil nil nil 17)
                  (basic/new-atom s1 atm18 nil nil nil nil 18)
                  (basic/new-atom s2 atm19 nil nil nil nil 19)
                  (basic/new-atom s1 atm20 nil nil nil nil 20)
                  (basic/new-atom s2 atm21 nil nil nil nil 21)
                  ;(basic/new-atom s1 atm22 nil nil nil nil 22)
                  ;(basic/new-atom s2 atm23 nil nil nil nil 23)
                  ;(basic/new-atom s1 atm24 nil nil nil nil 24)
])))






(defn Octafunctionalized-Biphenylenes-type2-Cpdos
  "Precurser synthesized in: 'Octafunctionalized Biphenylenes: Molecular Precursors for Isomeric Graphene Nanostructures'
  -by Florian Schlütter, Tomohiko Nishiuchi, Volker Enkelmann, and Klaus Müllen
   DOI: 10.1002/ange.201309324

  I used this function for the special purpose of defining certain carbon atoms differently for doing PDOS calculations in DFTB+.

  Usage: (Octafunctionalized-Biphenylenes-type2-Cpdos 'C 'Cc 1.539807 1.415841  1.394072  1.450921 1.436087 1.441179 1.443518)"
  [species1 species2 a b c d e f g]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        cs (cmat/cos (/ ed/pi 3))
        sn (cmat/sin (/ ed/pi 3))
        A [(+ b (* 2 c cs) d) 0.0 0.0]
        B [0.0 (+ a a (* 2 (+ c c e e g) sn)) 0.0]
        Aprime [(* 0.5 (+ b (* 2 c cs) d)) (* -0.5 (+ a a (* 2 (+ c c e e g) sn))) 0.0]
        Bprime [(* 0.5 (+ b (* 2 c cs) d)) (* 0.5 (+ a a (* 2 (+ c c e e g) sn))) 0.0]
        atm1 [(* 0.5 b) (* 0.5 a) 0.0]
        atm13 (+ A [(* -0.5 b) (* 0.5 a) 0.0])
        atm2 (+ atm1 [(* c cs) (* c sn) 0.0])
        atm14 (+ atm2 [d 0.0 0.0])
        atm3 (+ atm2 [(* -1.0 e cs) (* e sn) 0.0])
        atm15 (+ atm14 [(* e cs) (* e sn) 0.0])
        atm4 (+ atm3 [(* g cs) (* g sn) 0.0])
        atm16 (+ atm4 [f 0.0 0.0])
        atm5 (+ atm4 [(* -1.0 e cs) (* e sn) 0.0])
        atm17 (+ atm16 [(* e cs) (* e sn) 0.0])
        atm6 (+ atm5 [(* c cs) (* c sn) 0.0])
        atm18 (+ atm6 [b 0.0 0.0])
        atm7 (+ atm6 [0.0 a 0.0])
        atm19 (+ atm7 [b 0.0  0.0])
        atm8 (+ atm7 [(* -1.0 c cs) (* c sn) 0.0])
        atm9 (+ atm8 [(* e cs) (* e sn) 0.0])
        atm10 (+ atm9 [(* -1.0 g cs) (* g sn) 0.0])
        atm21 (+ atm9 [f 0.0 0.0])
        atm22 (+ atm21 [(* g cs) (* g sn) 0.0])
        atm20 (+ atm19 [(* c cs) (* c sn) 0.0])
        atm11 (+ atm10 [(* e cs) (* e sn) 0.0])
        atm23 (+ atm11 [d 0.0 0.0])
        atm24 (+ atm23 [(* c cs) (* c sn) 0.0])
        atm12 (+ atm11 [(* -1.0 c cs) (* c sn) 0.0])]
     (hash-map :lvs  [Aprime Bprime [0 0 30]]
            :mol [;(basic/new-atom s1 atm1 nil nil nil nil 1)
                  ;(basic/new-atom s2 atm2 nil nil nil nil 2)
                  ;(basic/new-atom s1 atm3 nil nil nil nil 3)
                  (basic/new-atom s2 atm4 nil nil nil nil 4)
                  (basic/new-atom s2 atm5 nil nil nil nil 5)
                  (basic/new-atom s1 atm6 nil nil nil nil 6)
                  (basic/new-atom s1 atm7 nil nil nil nil 7)
                  (basic/new-atom s2 atm8 nil nil nil nil 8)
                  (basic/new-atom s2 atm9 nil nil nil nil 9)
                  ;(basic/new-atom s2 atm10 nil nil nil nil 10)
                  ;(basic/new-atom s1 atm11 nil nil nil nil 11)
                  ;(basic/new-atom s2 atm12 nil nil nil nil 12)
                  ;(basic/new-atom s2 atm13 nil nil nil nil 13)
                  ;(basic/new-atom s1 atm14 nil nil nil nil 14)
                  ;(basic/new-atom s2 atm15 nil nil nil nil 15)
                  (basic/new-atom s2 atm16 nil nil nil nil 16)
                  (basic/new-atom s2 atm17 nil nil nil nil 17)
                  (basic/new-atom s1 atm18 nil nil nil nil 18)
                  (basic/new-atom s1 atm19 nil nil nil nil 19)
                  (basic/new-atom s2 atm20 nil nil nil nil 20)
                  (basic/new-atom s2 atm21 nil nil nil nil 21)
                  ;(basic/new-atom s1 atm22 nil nil nil nil 22)
                  ;(basic/new-atom s2 atm23 nil nil nil nil 23)
                  ;(basic/new-atom s1 atm24 nil nil nil nil 24)
])))









(defn Octafunctionalized-Biphenylenes-type3
  "From: 'Two dimensional Dirac carbon allotropes from graphene'
  -by Li-Chun Xu, et al.
   DOI: 10.1039/c3nr04463g

This is a rectangular cell structure of the primitive cell structure
given in the net-W function below.


  Usage: (Octafunctionalized-Biphenylenes-type3 'C 'C 1.539807 1.415841 1.394072 1.450921 1.436087)"
  [species1 species2 a b c d e]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        cs (cmat/cos (/ ed/pi 3))
        sn (cmat/sin (/ ed/pi 3))
        A [(+ b (* 2 c cs) e) 0.0 0.0]
        B [0.0 (+ a a (* 2 (+ c c d) sn)) 0.0]
        atm1  [(* -0.5 b) (*  0.5 a) 0.0]
        atm2  [(*  0.5 b) (*  0.5 a) 0.0]
        atm3  [(* -0.5 b) (* -0.5 a) 0.0]
        atm4  [(*  0.5 b) (* -0.5 a) 0.0]
        atm5  (+ atm1  [(* -1 cs c) (*    sn c) 0.0])
        atm6  (+ atm2  [(*    cs c) (*    sn c) 0.0])
        atm7  (+ atm5  [(*    cs d) (*    sn d) 0.0])
        atm8  (+ atm6  [(* -1 cs d) (*    sn d) 0.0])
        atm9  (+ atm7  [(* -1 cs c) (*    sn c) 0.0])
        atm10 (+ atm8  [(*    cs c) (*    sn c) 0.0])
        atm11 (+ atm3  [(* -1 cs c) (* -1 sn c) 0.0])
        atm12 (+ atm4  [(*    cs c) (* -1 sn c) 0.0])
        atm13 (+ atm11 [(*    cs d) (* -1 sn d) 0.0])
        atm14 (+ atm12 [(* -1 cs d) (* -1 sn d) 0.0])
        atm15 (+ atm13 [(* -1 cs c) (* -1 sn c) 0.0])
        atm16 (+ atm14 [(*    cs c) (* -1 sn c) 0.0])]
     (hash-map :lvs  [A B [0 0 30]]
            :mol (gmol/shift (* -0.5 (+ A B))
                 [(basic/new-atom s1 atm1 nil nil nil nil 0)
                  (basic/new-atom s2 atm2 nil nil nil nil 1)
                  (basic/new-atom s2 atm3 nil nil nil nil 2)
                  (basic/new-atom s1 atm4 nil nil nil nil 3)
                  (basic/new-atom s2 atm5 nil nil nil nil 4)
                  (basic/new-atom s1 atm6 nil nil nil nil 5)
                  (basic/new-atom s1 atm7 nil nil nil nil 6)
                  (basic/new-atom s2 atm8 nil nil nil nil 7)
                  (basic/new-atom s2 atm9 nil nil nil nil 8)
                  (basic/new-atom s1 atm10 nil nil nil nil 9)
                  (basic/new-atom s1 atm11 nil nil nil nil 10)
                  (basic/new-atom s2 atm12 nil nil nil nil 11)
                  (basic/new-atom s2 atm13 nil nil nil nil 12)
                  (basic/new-atom s1 atm14 nil nil nil nil 13)
                  (basic/new-atom s1 atm15 nil nil nil nil 14)
                  (basic/new-atom s2 atm16 nil nil nil nil 15)]))))



(defn net-W
  "From: 'Prediction of a new two-dimensional metallic carbon allotrope'
  -by Xin-Quan Wang, et al.
   Phys. Chem. Chem. Phys., 2013, 15:2024


  3-D Spacegroup: Cmmm (number 65)

  The origin in this puc is set at the center of the four member ring, thus many of the
  points fall outside of the parallelepiped described by the lattice vectors if their vertex is at the origin.

  Usage: (net-W 'C 'C 1.51695375 1.43764092 1.395370978 1.456896124 1.45869914)"
  [species1 species2 a b c d e]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        cs (cmat/cos (/ ed/pi 3))
        sn (cmat/sin (/ ed/pi 3))
        A [(+ b (* 2 c cs) e) 0.0 0.0]
        B [0.0 (+ a a (* 2 (+ c c d) sn)) 0.0]

        Aprime (- (* 0.5 A) (* 0.5 B))
        Bprime (- (+ (* 0.5 A) B) (* 0.5 B))

        atm1  [(* -0.5 b) (*  0.5 a) 0.0]
        atm2  [(*  0.5 b) (*  0.5 a) 0.0]
        atm3  [(* -0.5 b) (* -0.5 a) 0.0]
        atm4  [(*  0.5 b) (* -0.5 a) 0.0]
        atm5  (+ atm1  [(* -1 cs c) (*    sn c) 0.0])
        atm6  (+ atm2  [(*    cs c) (*    sn c) 0.0])
        atm7  (+ atm5  [(*    cs d) (*    sn d) 0.0])
        atm8  (+ atm6  [(* -1 cs d) (*    sn d) 0.0])
        atm9  (+ atm7  [(* -1 cs c) (*    sn c) 0.0])
        atm10 (+ atm8  [(*    cs c) (*    sn c) 0.0])
        atm11 (+ atm3  [(* -1 cs c) (* -1 sn c) 0.0])
        atm12 (+ atm4  [(*    cs c) (* -1 sn c) 0.0])
        atm13 (+ atm11 [(*    cs d) (* -1 sn d) 0.0])
        atm14 (+ atm12 [(* -1 cs d) (* -1 sn d) 0.0])
        atm15 (+ atm13 [(* -1 cs c) (* -1 sn c) 0.0])
        atm16 (+ atm14 [(*    cs c) (* -1 sn c) 0.0])]
     (hash-map :lvs  [Aprime Bprime [0 0 30]]
            :mol
                 [(basic/new-atom s1 atm1 nil nil nil nil 0)
                  (basic/new-atom s2 atm2 nil nil nil nil 1)
                  (basic/new-atom s2 atm3 nil nil nil nil 2)
                  (basic/new-atom s1 atm4 nil nil nil nil 3)
                  (basic/new-atom s1 atm7 nil nil nil nil 6)
                  (basic/new-atom s2 atm8 nil nil nil nil 7)
                  (basic/new-atom s2 atm13 nil nil nil nil 12)
                  (basic/new-atom s1 atm14 nil nil nil nil 13)])))











(defn net-W-Cpdos
  "From: 'Prediction of a new two-dimensional metallic carbon allotrope'
  -by Xin-Quan Wang, et al.
   Phys. Chem. Chem. Phys., 2013, 15:2024


  3-D Spacegroup: Cmmm (number 65)

  The origin in this puc is set at the center of the four member ring, thus many of the
  points fall outside of the parallelepiped described by the lattice vectors if their vertex is at the origin.

  Usage: (net-W-Cpdos 'Cc 'C 1.51695375 1.43764092 1.395370978 1.456896124 1.45869914)"
  [species1 species2 a b c d e]
  (let [s1 (.intern species1)
        s2 (.intern species2)
        cs (cmat/cos (/ ed/pi 3))
        sn (cmat/sin (/ ed/pi 3))
        A [(+ b (* 2 c cs) e) 0.0 0.0]
        B [0.0 (+ a a (* 2 (+ c c d) sn)) 0.0]

        Aprime (- (* 0.5 A) (* 0.5 B))
        Bprime (- (+ (* 0.5 A) B) (* 0.5 B))

        atm1  [(* -0.5 b) (*  0.5 a) 0.0]
        atm2  [(*  0.5 b) (*  0.5 a) 0.0]
        atm3  [(* -0.5 b) (* -0.5 a) 0.0]
        atm4  [(*  0.5 b) (* -0.5 a) 0.0]
        atm5  (+ atm1  [(* -1 cs c) (*    sn c) 0.0])
        atm6  (+ atm2  [(*    cs c) (*    sn c) 0.0])
        atm7  (+ atm5  [(*    cs d) (*    sn d) 0.0])
        atm8  (+ atm6  [(* -1 cs d) (*    sn d) 0.0])
        atm9  (+ atm7  [(* -1 cs c) (*    sn c) 0.0])
        atm10 (+ atm8  [(*    cs c) (*    sn c) 0.0])
        atm11 (+ atm3  [(* -1 cs c) (* -1 sn c) 0.0])
        atm12 (+ atm4  [(*    cs c) (* -1 sn c) 0.0])
        atm13 (+ atm11 [(*    cs d) (* -1 sn d) 0.0])
        atm14 (+ atm12 [(* -1 cs d) (* -1 sn d) 0.0])
        atm15 (+ atm13 [(* -1 cs c) (* -1 sn c) 0.0])
        atm16 (+ atm14 [(*    cs c) (* -1 sn c) 0.0])]
     (hash-map :lvs  [Aprime Bprime [0 0 30]]
            :mol
                 [(basic/new-atom s1 atm1 nil nil nil nil 0)
                  (basic/new-atom s1 atm2 nil nil nil nil 1)
                  (basic/new-atom s1 atm3 nil nil nil nil 2)
                  (basic/new-atom s1 atm4 nil nil nil nil 3)
                  (basic/new-atom s2 atm7 nil nil nil nil 4)
                  (basic/new-atom s2 atm8 nil nil nil nil 5)
                  (basic/new-atom s2 atm13 nil nil nil nil 6)
                  (basic/new-atom s2 atm14 nil nil nil nil 7)])))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;; Structures as reported in the literature ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(comment "Above I have endevored to produce functions that will allow the user to create graphitic
systems that are tailored to their needs, using the proper unit cells.  In this section I am
including structures that were published elsewhere.  This will both ease the user in comparing their
results with published results and will allow me to quickly add lots of graphitic structures that
researchers might find useful, but that I don't necessarily want to program in.")



#_(def ZhenhaiWang2015
  "Structures Phagraphene: A Low-energy Graphene Allotrope composed of 5-6-7 Carbon Rings with Distorted Dirac Cones.")



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;; Nanotubes ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn chiral-vector
  "This defines the chiral vector of a nanotube.  If you use the primitive unit cell lattice vectors of graphene, this will give you a typical carbon nanotube.  But you can use the unit vectors of other graphitic structures to obtrain different types of tubes."
  [n m puc]
  (let [[a1 a2 z] (:lvs puc)]
  (+ (* n a1) (* m a2))))


(defn metallic?
  "Used to determine if a carbon nanotube is metallic or semiconducting."
  [n m]
  (if (integer? (/ (- n m) 3))
    true
    false))


(defn nanotube-circumference
  [n m puc]
   (cmat/length (chiral-vector n m puc)))


(defn nanotube-radius
  "Computes the radius of a n,m-nanotube."
  [n m puc]
  (/  (nanotube-circumference n m puc) ed/tau))


(defn d_R
  "Used in defining T-vector."
  [n m]
  (gcmath/gcd (+ n n m) (+ n m m)))


(defn natoms-nanotube
  "computes the number of atoms in a the unit cell of a nanotube"
  [n m puc]
  (let [natom-puc (count (:mol puc))]
    (/ (* 2 natom-puc (+ (* n n) (* m m)(* n m))) (d_R n m))))



(defn T-vector
  "Vector perpendicular to the chiral vector, but still in the graphene plane.  This is also the direction of the tube axis."
  [n m puc]
  (let [[a1 a2 z] (:lvs puc)]
  (+ (* (+ n m m) (/ 1 (d_R n m)) a1) (* (+ n n m) (/ -1 (d_R n m)) a2))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;The method of Carter White, doi = {10.1103/PhysRevB.47.5485} to create a simple nanotube.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn primitive-cell-onto-tube
  "puc stands for the primitive unit cell-this includeds the lvs as the first element and the mol as the second element."
  [n m puc]
  (let [[a1 a2 z] (:lvs puc)
        mol (:mol puc)
        r (nanotube-radius n m puc)
        C (chiral-vector n m puc)
        T (cmat/normalise (T-vector n m puc))
       rot-angle (fn [p] (* 2 ed/pi (/ 1 (cmat/length C))(/ 1 (cmat/length C)) (cmat/dot p C)))
       rot (fn [p] (gmath/the-rotation-function [0 0 (+ r (last p))] [0 0 0] T (rot-angle p)))
       translation (fn [p] (* (/ (cmat/length (cmat/cross p C)) (cmat/length C)) T))]
    (do (map #(comp (partial println rot-angle) :coordinates) mol)
    (gmol/update-mol :coordinates #(+ (rot %) (translation %))  mol))))



(defn create-ring
 ""
  [n m puc]
  (let [T (cmat/normalise (T-vector n m puc))
        d (gcmath/gcd n m)
        alpha2 (/ ed/tau d)
        mol (primitive-cell-onto-tube n m puc)]
    (flatten (map #(gmol/rotate-mol mol [0 0 0] T (* alpha2 %)) (range 1 (inc d))))))





(defn findd-h1h2
  "This is used as part of defining the screw."
  [n m]
  (if (or (zero? n) (zero? m))
    [1 1]
  (let [h1 (fn [h2] (/ (- (* n h2) (gcmath/gcd n m)) m))]
  (loop [x (iterate inc 1)]
    (if (integer? (h1 (first x)))
      [(h1 (first x)) (first x)]
      (recur (rest x)))))))



(defn screw
    [n m puc]
  (let [[a1 a2 z] (:lvs puc)
        C (chiral-vector n m puc)
        T (cmat/normalise (T-vector n m puc))
        H (#(+ (* (first %) a1) (* (second %) a2)) (findd-h1h2 n m))
       rot-angle (* ed/tau (/ 1 (cmat/length C))(/ 1 (cmat/length C)) (cmat/dot H C))
       rot (fn [p] (gmol/rotate-mol p [0 0 0] T rot-angle))
       translation (* (/ (cmat/length (cmat/cross H C)) (cmat/length C)) T)]
    [translation  #(gmol/shift (rot %) translation)]))




#_(defn create-nanotube-fast
  "This uses the method of Carter White, doi = {10.1103/PhysRevB.47.5485} to
  create a simple nanotube.  This is the best method to use for creating simple nanotubes.

  When moving from pre-release greenwood to greenwood this function broke."
  [n m puc]
  (let [[a1 a2 z] (:lvs puc)
        h (create-ring n m puc)
        s (screw n m puc)
        T (T-vector n m puc)
        scrw #(apply comp (take % (repeat (second s))))
        mm (/ (natoms-nanotube n m puc)(count h))]
    (do (println (natoms-nanotube n m puc))
    (hash-map :lvs
    [[1000 0 0] [0 1000 0] [0 0 (cmat/length T)]]
  :mol (->> (flatten [h (map #((scrw %) h) (range 1 mm))])
    (gmol/apply-coord-transform-matrix (gmath/rotate-vec-to-axis T :z))
      (gmol/mol-center ))))))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Creating a nanotube by producing a flat sheet of graphene (or the like) and rolling it up.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn same-point?
  "This determines if two atoms in a mol, atom1 and atom2, lay on the same crystalagraphic
  point, assuming that those atoms are only a linear displacement away from each other.
  It does this by checking to see if the atoms are the same species and then if they are,
  checking to see if the same vectors may be used to get from one of the atoms to all of
  it's gneigh/neighbors."
  [mol atom1 atom2]
  (if (= (:species atom1) (:species atom2))
  (let [a1n (map #(- (:coordinates atom1) (:coordinates %)) (gmol/take-mol-by-pos mol ((comp (fn [x] (map :npos x)) :neigh) atom1)))
        a2n (map #(- (:coordinates atom2) (:coordinates %)) (gmol/take-mol-by-pos mol ((comp (fn [x] (map :npos x)) :neigh) atom2)))]
    (every? (fn [x] (some #(gmath/vectors-equal? % x 1.0E-4) a1n)) a2n))
    false))




(defn nearest-similar-crystal-point
  [mol atom1]
  (let [a (filter (partial same-point? mol atom1) mol)
        b (gmol/mol-filter-not {:pos (:pos atom1)} a)
        distancevec (map #(cmat/distance (:coordinates atom1) (:coordinates %))  b)
        pos (gutils/positions #{(apply min distancevec)} distancevec)]
    (nth b (first pos))))




(defn unrolled-nanotube
  "Creating a flat sheet of graphene (or the like) which can be rolled up to make a
  nanotube. n and m are the normal n and m values used in defining a carbon nanotube.
  puc is the primitive unit cell (in units of Angstroms) of the structure that is to
  be rolled into a tube.  The puc is a hash-map of :lvs and :mol keyword/values.
  It is assumed that the unit cell is parallelogram of the shape used in the graphene
  primitive unit cell.

  Note that this produces a mol that is larger in the y-direction than what the
  lvs says, this is to make sure that it works for making nanotubes.  If you wish
  to use this in creating a 2-D supercell then you will need to use them
  unrolled->parallelopipedsc function.

  (n,0) zigzag nanotube, (n,n) armchair nanotube, (n,m) chiral nanotube."
  [n m puc]
  (let [chiral (chiral-vector n m puc)
        LC (cmat/length chiral)
        T (T-vector n m puc)
        LT (cmat/length T)
        rot-mat (gmath/rotate-vec-to-axis chiral :x)
        lvs [(* 1.001 chiral) (* 5.0 T) [0 0 20]]
        cell (define-cell lvs [0 0 -10])]
       (basic/unitcell [[LC 0 0] [0 (* 3 LT) 0] [0 0 LC]]
                              (as-> (computation-supercell (:mol puc)  (:lvs puc) 40 40 0) x
                                    (:mol x)
                                    (gmol/mol-filter {:coordinates (partial within-cell? cell)} x)
                                    (atom-pos x)
                                    (gmol/apply-coord-transform-matrix rot-mat x)))))




(defn unrolled->parallelopipedsc
"Creates a supercell of a graphitic material where the size and orientation of the supercell
uses the n and m values used in nanotubes.  n must be greater than or equal to m."
 [n m puc]
 (let [b (unrolled-nanotube n m puc)
       blvs [(first (:lvs b)) [0 (/ (second (second (:lvs b))) 3.0) 0] (last (:lvs b))]]
   (basic/unitcell blvs
   (xyz/atom-pos (gmol/mol-filter {:coordinates (partial within-cell?? blvs [0 0 0])} (:mol b))))))









(defn create-nanotube
  "This is the function to use when you want to create a more complicated nanotube
  (ie. Fgraphene nanotubes, nanotubes made out of biphenylene-carbon).  puc is the
  primitive unit cell (in units of Angstroms) of the structure that is to be rolled
  into a tube.  The puc is a hash-map of :lvs and :mol keyword/values.  It is assumed
  that the unit cell is parallelogram of the shape used in the graphene primitive unit cell.
  (n,0) zigzag nanotube, (n,n) armchair nanotube, (n,m) chiral nanotube."
  [n m puc]
  (let [unrolled (unrolled-nanotube n m puc)
        lvs (:lvs unrolled)
        C (first lvs)
        T (second lvs)
        r (/ (cmat/length C) ed/tau)
       rot-angle (fn [p] (/ (first p) r))
       rot (fn [p] (+ [0 (second p) 0] (gmath/the-rotation-function [0 0 (- r (last p))] [0 0 0] T (rot-angle p))))
       tube (-> (gmol/update-mol :coordinates rot (:mol unrolled))
                (atom-pos )
                (gneigh/remove-overlapping 0.4)
                (gneigh/neighbors 0.1 1.6))
        a (map #(gneigh/nearest-atom-point tube (gmath/mat-vect-mult (gmath/card-rot-mat (* ed/pi %) :y) [0 (* 0.25 (second T)) r] )) [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.8 1.0])
        dis (- (reduce max (map #(cmat/length (- (:coordinates %)
                      (:coordinates (nearest-similar-crystal-point tube %)))) a))  0.1)
        p-tube (gmol/mol-filter {:coordinates #(and (>= (second %) 0) (<= (second %) (cmat/abs dis)))} tube)
        tubelvs [C [0 dis 0] (last lvs)]]
    (hash-map :lvs tubelvs :mol (atom-pos (drop-overlapping-sc-atoms p-tube tubelvs 0 1 0)))))




#_(defn create-spiral-multiwalled-nanotube
  "It has been suggested [26] that the presence of dislocation-like defects in the
  scroll type nanotubes is responsible for the transition from the scroll type to
  the nested type multiwalled carbon nanotube.

  Currently a work in progress."
    [n m puc nx]
  (let [unrolled (unrolled-nanotube n m puc)
        lvs (:lvs unrolled)
        C (first lvs)
        T (second lvs)
        r (/ (cmat/length C) ed/tau)
       rot-angle (fn [p] (/ (first p) r))
       rot (fn [p] (+ [0 (second p) 0] (gmath/the-rotation-function [0 0 (- r (last p))] [0 0 0] T (rot-angle p))))
       tube (-> (gmol/update-mol (:mol unrolled) :coordinates rot)
                (atom-pos )
                (gneigh/remove-overlapping 0.4)
                (gneigh/neighbors 0.1 1.6))
        a (gneigh/nearest-atom-point tube [0 (* 0.25 (second T)) r])
        dis (cmat/length (- (:coordinates a)
                      (:coordinates (nearest-similar-crystal-point tube a))))
        p-unrolled (create-supercell
                    (gmol/mol-filter {:coordinates #(and (>= (second %) 0) (<= (second %) (cmat/abs dis)))} (:mol unrolled))
                    (cell-projectors [C [0 dis 0] (last lvs)] nx 1 1))]

    ))





#_(defn supercell->nanotube
"This is the function to use when you have a rectangular supercell and want to
    turn it into a nanotube.  puc is the structure that is to be rolled into a
    tube; it is a hash-map of :lvs and :mol keyword/values.

  Currently it hasn't been tested."
[puc]
(let [lvs (:lvs puc)
      C (first lvs)
      T (second lvs)
      r (/ (cmat/length C) ed/tau)
     rot-angle (fn [p] (/ (first p) r))
     rot (fn [p] (+ [0 (second p) 0] (gmath/the-rotation-function [0 0 (- r (last p))] [0 0 0] T (rot-angle p))))
     tube (gmol/update-mol (:mol puc) :coordinates rot)]
  (hash-map :lvs [C [0 dis 0] (last lvs)] :mol tube)))


;;;;;;;;;;;;;;;;;;;;;;;;; turbostratic stacking ;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(defn mn->theta
"I found this in 'Bending modes, elastic constants and mechanical stability of graphitic systems' by G. Savini, Y.J. Dappe, S. O ̈ berg, J.-C. Charlier, M.I. Katsnelson, A. Fasolino.
For this to work properly they state the condition n>m must hold."
  [n m]
  (cond (> n m)
  (acos
  (/
  (+ (* 2 n n) (* 2 n m) (* -1 m m))
  (+ (* 2 n n) (* 2 n m) (* 2 m m))))
  :else nil))






#_(defn turbostratic-stacking
    "This assumes a cell structure similar to graphene."
  [n m h puc1 puc2]
  (let [chiral (chiral-vector n m puc1)
        theta (mn->theta n m)
        rot-mat (gmath/rotate-vec-to-axis chiral :x)
        lvs [(* (cmat/length chiral) (cmat/normalise (first (:lvs puc1))))
                        (* (cmat/length chiral) (cmat/normalise (second (:lvs puc1))))
                        [0 0 10]]
        avec (* -0.5 (+ ((comp first :lvs) puc1) ((comp second :lvs) puc1)))]
    (hash-map
       :puc1
       (basic/unitcell lvs
                       (as-> (create-supercell (:mol puc1) (computation-projectors (:lvs puc1) 40 40 0)) x
                                    (gmol/apply-coord-transform-matrix rot-mat x)
                                    (gmol/rotate-mol x [0.0 0.0 0.0] [0.0 0.0 1.0] (ed/degrees->radians -30))
                                    (gmol/shift [0.0 0.0 0.1] x)
                                    (gmol/mol-filter {:coordinates (partial within-cell?? lvs avec)} x)))
       :puc2
       (basic/unitcell lvs
                        (as-> (create-supercell (:mol puc2) (computation-projectors (:lvs puc1) 40 40 0)) x
                                    (gmol/apply-coord-transform-matrix rot-mat x)
                                    (gmol/rotate-mol x [0.0 0.0 0.0] [0.0 0.0 1.0] (* -1 (mn->theta n m)))
                                    (gmol/rotate-mol x [0.0 0.0 0.0] [0.0 0.0 1.0] (ed/degrees->radians -30))
                                    (gmol/shift [0.0 0.0 (+ 0.1 h)] x)
                                    (gmol/mol-filter {:coordinates (partial within-cell?? lvs avec)} x))))))


;(def b (turbostratic-stacking 3 2 3.0 a ap))

;(def c  (supercell (:mol b) (:lvs b) 3 3 1))

;(spit "/Users/chadjunkermeier/Desktop/graphene.xyz" (out/write-xyz (:mol c)))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   Adding defects to a graphene sheet
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn top-site-adsorption
  "gmol is the mol of the graphitic sheet
  admol is the mol of the adsorbate
  nsites is a coll of integers, :pos, to which the adsorbates are bound.
  height is the length of the adsorbate C bond.
  Usage: (def ggg (second (make-zigzag-graphene-supercell 'C' 'C' 'C' 'C' 1.421 2 2 false)))
         (top-site-adsorption ggg (xyz-str->atoms 'F 0 0 0') [1 2] 1.44)"
  [grmol admol nsites height]
  (let [ngrmol (count grmol)
        added (flatten (doall (map #(as-> (:coordinates %) x
                      (+ [0 0 height] x)
                       (gmol/shift-to x 0 admol))
                        (gmol/take-mol-by-pos grmol nsites))))
        a-mol (gmol/col->mol  :pos (doall (take (count added) (iterate inc ngrmol))) added)]
  (atom-pos (flatten [grmol a-mol]))))







(defn random-up-down
  [h]
  ({0 h 1 (- h)} (rand-int 2)))


(defn random-topsite-adsorption
  "This is designed to work with graphene, but it should work with all surfaces
  assuming the normal to the surface is in the z-direction, you pass in only the
  surface atoms into grmol, and you set one-two-key to :one.
  Usage: (random-topsite-adsorption graphene (xyz-str->atoms 'F 0 0 0') 1/4 1.4 :two)"
  [grmol admol percentage height one-two-key]
  (let [nC (count grmol)
        nf (int (* percentage nC))
        fpos (take nf (shuffle (range nC)))
        added (flatten (doall (map #(as-> (:coordinates %) x
                      (+ [0 0 ({:one height :two (random-up-down height)} one-two-key)] x)
                       (gmol/shift-to admol 0 x)) (gmol/take-mol-by-pos grmol fpos))))
        a-mol (gmol/col->mol added :pos (doall (take (count added) (iterate inc nC))))]
  (gneigh/neighbors (flatten [grmol a-mol]) 0.1 1.8)))




(defn pre-optimize-adatom
 "This is a poor mans optimization.  Really only designed to work with graphene.
 It determines which side of the graphene plane the adatom is on and then moves the
 adatom and the C atom in that direction.
You must have run gneigh/neighbors on the mol before using this.  It does not check to see
if you have, and if you haven't it won't optimize, it just spits the mol back out.
  YOU MUST HAVE THE GRAPHENE SHEET SET AT z=0."
 [mol species v]
 (let [func #(map :npos (gmol/mol-filter {:nspecies "C"} (flatten (map :neigh %))))
         upf (gmol/update-mol :coordinates #(+ v %) (gmol/mol-filter {:species species :coordinates #(pos? (last %))} mol))
       downf (gmol/update-mol :coordinates #(- % v) (gmol/mol-filter {:species species :coordinates #(neg? (last %))} mol))
         upC (gmol/update-mol :coordinates #(+ % v) (gmol/mol-filter-vec :pos (func upf) mol))
       downC (gmol/update-mol :coordinates #(- % v) (gmol/mol-filter-vec :pos (func downf) mol))
       moved (concat upf upC downf downC)]
   (concat (gmol/mol-filter-not-vec :pos (map :pos moved) mol) moved)))




(defn pre-optimize-neighboring-adatoms
  [mol species]
  (let [f #(some (partial = "F") (:nspecies %))
        F-to-move (gneigh/neighbors (gmol/mol-filter {:species species :neigh f} mol) 0.1 1.6)
        pairs (mapv #(vector (:pos %) ((comp first :npos :neigh) %)) F-to-move)]
    (concat (gmol/mol-filter-not {:species species :neigh f} mol)
            (flatten (map #(gmol/take-mol-by-pos (gmol/resize-bond F-to-move (first %) (second %) (+ 0.4 (cmat/length (gmol/mol-vector F-to-move (first %) (second %))))) [(first %)]) pairs)))))







(defn bond-centered-adsorption
  "This will place an adsorbate 'above' a bond gmath/midpoint (e.g. an epoxy group over
the bond connecting two C atoms in graphene). As is, this is incomplete.  I would
like to change this so that it uses a routine to find the normal to the surface;
which would allow for irregular surfaces.
  bond is a 2-tuple of the atoms' :pos"
  [surface-mol bond adsorbate-mol ads-1 height]
  (let [pos (apply gmath/midpoint (map :coordinates (gmol/mol-filter-vec :pos bond surface-mol)))
        normal [0 0 1]
        a-mol (atom-pos adsorbate-mol)
        ad-mol (gmol/shift-to (+ pos (* (+ 0.0 height) normal))  ads-1  a-mol )
        moved (gmol/update-mol :coordinates #(+ [0 0 0.0] %) (gmol/mol-filter-vec :pos bond surface-mol))]
    (atom-pos (sort-by :pos (flatten [(gmol/mol-filter-not-vec :pos bond surface-mol) moved (atom-pos ad-mol (count surface-mol))])))))



#_(defn bond-centered-adsorption
  "This will place an adsorbate 'above' a bond gmath/midpoint (e.g. an epoxy group over
the bond connecting two C atoms in graphene). As is, this is incomplete.  I would
like to change this so that it uses a routine to find the normal to the surface;
which would allow for irregular surfaces.
  bond is a 2-tuple of the atoms' :pos"
  [surface-mol bond adsorbate-mol ads-1 height]
  (let [pos (apply gmath/midpoint (map :coordinates (gmol/mol-filter-vec :pos bond surface-mol)))
        normal [0 0 1]
        a-mol (atom-pos adsorbate-mol )
        ad-mol (gmol/shift-to (+ pos (* (+ 0.0 height) normal))  ads-1  a-mol )
        moved (gmol/update-mol :coordinates #(+ [0 0 0.0] %) (gmol/mol-filter-vec :pos bond surface-mol))]
    (flatten [surface-mol ad-mol ])))












#_(defn H-terminated-substitional-defect
  "In order to use this, gneigh/neighbors needs to have been run on mol previous to
using this command."
  [mol atom num-H]
  (let [C-H-length 1.09
        C-C-length 1.421
        neighboring-C (flatten (map #(gmol/mol-filter :species "C" (gmol/take-mol-by-pos mol [%]))
                        (flatten (map (comp :npos :neigh)
                          (gmol/take-mol-by-pos mol [(:pos atom)])))))
        C-C-bonds (map #(map - (:coordinates atom) (:coordinates %)) neighboring-C)
        Hbonds1 (#(map + (*
                    (* C-H-length (/ (cmat/length %1) C-C-length))
                    (cmat/normalise %1)) (:coordinates %2))
                 (first C-C-bonds) (first neighboring-C))
        Hbonds2 (#(map + (*
                    (* 0.8660254037844387 C-H-length (/ (cmat/length %1) C-C-length))
                    (cmat/normalise %1))
                    (:coordinates %2)
                    (*
                    (* 0.5 C-H-length (/ (cmat/length %1) C-C-length))
                    [0 0 1]))
                 (second C-C-bonds) (second neighboring-C))
        Hbonds3 (#(map + (*
                    (* 0.8660254037844387 C-H-length (/ (cmat/length %1) C-C-length))
                    (cmat/normalise %1))
                    (:coordinates %2)
                    (*
                    (* 0.5 C-H-length (/ (cmat/length %1) C-C-length))
                    [0 0 -1]))
                 (last C-C-bonds) (last neighboring-C))]
    (cond
      (== num-H 1)
      (concat (gmol/mol-filter-not :pos (:pos atom) mol)
        (basic/new-atom "H" Hbonds1 nil nil nil nil nil))
      (== num-H 2)
      (concat (gmol/mol-filter-not :pos (:pos atom) mol)
        (map #(basic/new-atom "H" % nil nil nil nil nil) (vector Hbonds2 Hbonds3)))
      (== num-H 3)
      (concat (gmol/mol-filter-not :pos (:pos atom) mol)
        (map #(basic/new-atom "H" % nil nil nil nil nil) (vector Hbonds1 Hbonds2 Hbonds3))))))





(defn- rand-GO-carbons-
  [mol col]
  (loop [c1 (rand-nth col)
         [c2 c3] ((comp (partial take 2) :npos :neigh first) (gmol/take-mol-by-pos mol [c1]))
         i 1]
        (if (cset/subset? (set [c1 c2 c3]) (set col))  [c1 c2 c3]
          (if (= i 150) (do (println "rand-GO Failed") false)
          (recur (rand-nth col)
                 ((comp (partial take 2) :npos :neigh first) (gmol/take-mol-by-pos mol [c1]))
                 (inc i))))))



(defn random-GO
  "This creates graphene oxide with random placements of the atoms."
[mol proportion]
(let [nC (count mol)]
  (loop [moll mol
         m (map :pos mol)
         i 1]
    (let [a (rand-GO-carbons- mol m)
          h (random-up-down 1)]
    (if (or (false? a) (> (* 2 i) (* nC proportion))) moll
      (recur
       (-> moll
           (bond-centered-adsorption [(first a)(second a)] (xyz-str->atoms "O 0 0 0") 0 (* h 1.2))
           (top-site-adsorption (gneigh/neighbors (xyz-str->atoms "O 0 0 0\nH 0.96 0 0") 0.2 1.4) [(last a)] (* h -1.44)))
       (vec (cset/difference (set m) (set a)))
       (inc i)))))))







(defn tri-layer-F-graphene
  "This is to make tri-layer graphene.  This allows for the top and bottom
layers to have different coverages.  Currently, I assume that
all three layers have the same a lattice constant."
  [top middle bottom distance]
  (concat
    (gmol/shift [0 0 distance] top)
    middle
      (gmol/shift [0 0 (- distance)] bottom)))








#_(defn Weighted-HOMO-LUMO-kinetic-stability
  "This is based on the the paper Theor Chem Acc (1999) 102:134±138,  DOI 10.1007/s002149800m93.
  Usage: (Weighted-HOMO-LUMO-kinetic-stability fullerene -0.1 -2.2)
    I think there description leaves something to be desired and this won't work without a better understanding of what is m_{HOMO} and m_{LUMO}."
  [mol HOMO-LUMO-diff]
  (* (gneigh/count-bonds mol) HOMO-LUMO-diff))




(defn graphene-resonator-frequency
  "This computes the frequency of a graphene resonator with a gate voltage of zero.
The form for this equation comes from Nature Nanotechnology 4, 861 - 867 (2009),
http://www.nature.com/nnano/journal/v4/n12/full/nnano.2009.267.html.
This needs the mol and lvs to determine the mass density.  In this we will
actually calculate the f/f-ref, where L,w,To are respectively equal to L-ref,
w-ref, To-ref, and Te = 0."
  [mol-ref lvs-ref mol lvs]
  (let [mass  #(reduce + (map (comp ed/atomic-mass ed/atomic-numbers :species) %))]
    (sqrt (/
      (/ (mass mol-ref) (gmath/lvs-volume lvs-ref))
      (/ (mass mol) (gmath/lvs-volume lvs))))))




;;;;;;;;;;;;;;;;;;;;;;;;; Patterning of adsorbates ;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn middle-half?
  [min-max vectors]
  (let [middle (gmath/midpoint [(first min-max)][(second min-max)])
        lower-quartarian (first (gmath/midpoint [(first min-max)] middle))
        higher-quartarian (first (gmath/midpoint middle [(second min-max)]))]
    (and
      (> (first vectors) lower-quartarian)
      (< (first vectors) higher-quartarian))))


(defn middle-half-cheat?
  [min-max cheat-distanceL cheat-distanceR vectors]
  (let [middle (gmath/midpoint [(first min-max)][(second min-max)])
        lower-quartarian (first (gmath/midpoint [(first min-max)] middle))
        higher-quartarian (first (gmath/midpoint middle [(second min-max)]))]
    (and
      (> (first vectors) (+ lower-quartarian cheat-distanceL))
      (< (first vectors) (- higher-quartarian cheat-distanceR)))))

(defn middle-third?
  [min-max vectors]
  (let [thir (/ (- (second min-max) (first min-max)) 3.)]
    (and
      (> (first vectors) (+ (first min-max) thir))
      (< (first vectors) (+ (first min-max) thir thir)))))



(defn four-square?
    "This function assumes that the graphene sheet is centered at the origin.
  Usage: (patchwork-supercell graphene-sc CO-sc #(four-square? % ))"
  [vectors]
    (or
      (and (< (first vectors) 0)
        (< (second vectors) 0))
      (and (>= (first vectors) 0)
        (>= (second vectors) 0))))


(defn big-circle?
  "Circle centered at zero with a radius close to the size of the cell."
  [min-max vectors]
  (let [min (first (sort (map #(abs (reduce - %))(partition 2 (drop-last 2 min-max)))))]
    (<= (length vectors) (* min 0.5))))

(defn half-width-circle? [min-max vectors]
  (let [min (first (sort (map #(abs (reduce - %))(partition 2 (drop-last 2 min-max)))))]
    (< (length vectors) (* 0.25 min))))


(defn quarter-width-circle? [min-max vectors]
  (let [min (first (sort (map #(abs (reduce - %))(partition 2 (drop-last 2 min-max)))))]
    (< (length vectors) (* 0.125 min))))
