(ns greenwood.empirical-data
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :as cmat]
            [clojure.core.matrix.operators :as cmato]))



(def ^:const Avogadro-number (cmato/* 6.02214179 (cmat/pow  10. 23)))
(def  ^:const kB 0.0000861734315);Boltzmann constant units of eV/K
(def ^:const  pi 3.141592653589793)
(def ^:const tau 6.283185307179586)

(defn Angstrom->Bohr [x] (cmato/* x 1.889682673))
(defn Bohr->Angstrom [x] (cmato/* x 0.529189379))



(defn Ry->eV [x] (cmato/* x 13.605698066))
(defn eV->Ry [x] (cmato// x 13.605698066))
(defn kelvin->eV [T] (cmato/* kB T))

(Ry->eV (cmato/- -3331.31268534  -3331.321586419))

(defn kcalpermol->eV [x] (cmato/* x 0.04336))
(defn eV->kcalpermol [x] (cmato// x 0.04336))

(defn Ry->kcalpermol [x] (eV->kcalpermol (Ry->eV x)))


(defn kJpermol->eV [x] (cmato/* x 0.010153795839444585))
(defn eV->kJpermol [x] (cmato// x 0.010153795839444585))

(defn J->eV [x] (cmato/* x 6.24150974E18))

(defn eV->J [x] (cmato// x 6.24150974E18))

(defn eV->K [x] (cmato/* x 11604.5))

(defn eA->D [x] (cmato/* x 0.20819434)) ;Debye (D), unit used in dipole moment

(defn radians->degrees [x]
  (let [d (cmato/* x 57.2957795130823)]
    (cond (> d 360.0)
           (cmato/- d 360.0)
          (< d -360.0)
          (cmato/+ d 360.0)
          :else d)))
(defn degrees->radians [x] (cmato/* x 0.0174532925199433))

(defn Angstrom->meters [x] (cmato/* x 1.0E-10))

(defn Hartree->kcalpermol [x] (cmato/* x 627.5095))
(defn Hartree->eV [x] (cmato/* x 27.2116))
(defn Hz->inverse_cm [x] (cmato/* x 29979245800.0))

(eV->kcalpermol 9.17);Hydroxyl
(eV->kcalpermol 4.68);Epoxy
(eV->kcalpermol (kJpermol->eV 750)); H—F
(eV->kcalpermol (kJpermol->eV 190)); O—F

;(eV->kcalpermol (kelvin->eV 1000))


;;;;;;;;;;;;;;;;;;;;   PERIODIC VALUES EMPIRICAL DATA   ;;;;;;;;;;;;;;;;;;;;

(def periodic-atomic-symbols ["H", "He", "Li", "Be", "B", "C", "N", "O",
            "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
            "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
            "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
            "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
            "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
            "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
            "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
            "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
                        "Uun", "Uuu", "Uub","", "Uuq", "", "Uuh", "", "Uuo"])

(def periodic-atomic-radius
  #^{:doc "J.C. Slater (1964). Atomic Radii in Crystals. J. Chem.
           Phys. 41: 3199. doi:10.1063/1.1725697"}
  [0.35 nil 1.45 1.05 0.85 0.7 0.65 0.6 0.5 nil 1.8 1.5 1.25 1.1 1. 1.
   1. 0.71 2.2 1.8 1.6 1.4 1.35 1.4 1.4 1.4 1.35 1.35 1.35 1.35 1.3 1.25
   1.15 1.15 1.15 nil 2.35 2. 1.8 1.55 1.45 1.45 1.35 1.3 1.35 1.4 1.6
   1.55 1.55 1.45 1.45 1.4 1.4 nil 2.6 2.15 1.95 1.85 1.85 1.85 1.85 1.85
   1.85 1.8 1.75 1.75 1.75 1.75 1.75 1.75 1.75 1.55 1.45 1.35 1.35 1.3 1.35
   1.35 1.35 1.5 1.9 1.8 1.6 1.9 nil nil nil 2.15 1.95 1.8 1.8 1.75 1.75
   1.75 1.75 nil nil nil nil nil nil nil nil nil nil nil nil nil nil nil
   nil nil nil nil nil nil])

(def periodic-atomic-mass [1.00794 4.002602 6.941 9.012182 10.811 12.0107
                           14.0067 15.9994 18.9984032 20.1797 22.98976928
                           24.305 26.9815386 28.0855 30.973762 32.065 35.453
                           39.948 39.0983 40.078 44.955912 47.867 50.9415
                           51.9961 54.938045 55.845 58.6934 58.933195 63.546
                           65.38 69.723 72.64 74.9216 78.96 79.904 83.798
                           85.4678 87.62 88.90585 91.224 92.90638 95.96 98
                           101.07 102.9055 106.42 107.8682 112.411 114.818
                           118.71 121.76 127.6 126.90447 131.293 132.9054519
                           137.327 138.90547 140.116 140.90765 144.242 145
                           150.36 151.964 157.25 158.92535 162.5 164.93032
                           167.259 168.93421 173.054 174.9668 178.49 180.94788
                           183.84 186.207 190.23 192.217 195.084 196.966569
                           200.59 204.3833 207.2 208.9804 210 210 222 223 226
                           227 232.03806 231.03588 238.02891 237 244 243 247
                           247 251 252 257 258 259 262 261 262 266 264 267 268
                           271 272 285 284 289 288 292 294])



(def periodic-neutral-charges [[1] [2] [1] [2] [2 1] [2 2] [2 3] [2 4] [2 5] [2 6]])







;;;;;;;;;;;;;;;;;;;;   BINDING VALUES PERIODIC VALUES   ;;;;;;;;;;;;;;;;;;;;


#_(def ^:dynamic *periodic-atomic-symbols* (zipmap (range 1 (count periodic-atomic-symbols)) periodic-atomic-symbols))
#_(def ^:dynamic *periodic-atomic-radii* (zipmap (range 1 (count periodic-atomic-radius)) periodic-atomic-radius))
#_(def ^:dynamic *periodic-atomic-mass* (zipmap (range 1 (count periodic-atomic-mass)) periodic-atomic-mass))
#_(def ^:dynamic *periodic-charges* (zipmap (range 1 (count periodic-neutral-charges)) periodic-neutral-charges))


(def periodic-atomic-symbols (zipmap (range 1 (count periodic-atomic-symbols)) periodic-atomic-symbols))
(def periodic-atomic-radii (zipmap (range 1 (count periodic-atomic-radius)) periodic-atomic-radius))
(def periodic-atomic-mass (zipmap (range 1 (count periodic-atomic-mass)) periodic-atomic-mass))
(def periodic-charges (zipmap (range 1 (count periodic-neutral-charges)) periodic-neutral-charges))


#_(defn set-atomic-symbols [hashmap]
  "This can be used to override the values hardcoded into this ns.  The user
should supply a hashmap with the key values desired.  It is expected that this
could be used to increase the efficiency of a calculation, or to define pseudo-atoms
that have special traits."
  (def ^:dynamic *periodic-atomic-symbols* hashmap))

#_(defn set-atomic-radii [hashmap]
  "This can be used to override the values hardcoded into this ns.  The user
should supply a hashmap with the key values desired.  It is expected that this
could be used to increase the efficiency of a calculation, or to define pseudo-atoms
that have special traits."
  (def ^:dynamic *periodic-atomic-radii* hashmap))

#_(defn set-atomic-mass [hashmap]
  "This can be used to override the values hardcoded into this ns.  The user
should supply a hashmap with the key values desired.  It is expected that this
could be used to increase the efficiency of a calculation, or to define pseudo-atoms
that have special traits."
  (def ^:dynamic *periodic-atomic-mass* hashmap))


#_(defn set-atomic-charges [hashmap]
  "This can be used to override the values hardcoded into this ns.  The user
should supply a hashmap with the key values desired.  It is expected that this
could be used to increase the efficiency of a calculation, or to define pseudo-atoms
that have special traits."
  (def ^:dynamic *periodic-charges* hashmap))




;;;;;;;;;;;;;;;;;;;;   USING PERIODIC VALUES DATA   ;;;;;;;;;;;;;;;;;;;;

(defn isnil? [value string]
  (if (empty? value)
    (println string)
    (value)))

(defn atomic-symbols [nZ]
  "nZ is an atomic number.  The output will be a string of the corresponding
atomic symbol.

Usage: (atomic-symbols 1)  => H"
    (periodic-atomic-symbols nZ))



(defn atomic-numbers [Xx]
    "Xx is a string containing an atomic symbol.  The output will be a
integer of the corresponding atomic number.

Usage: (periodic-atomic-symbols H) => 1"
    ((zipmap (vals periodic-atomic-symbols) (keys periodic-atomic-symbols)) Xx))



(defn atomic-radius [nZ]
  "Gives the atomic radius for any element. nZ is an atomic number. The output
will be a real value corresponding to the atomic radius of nZ.

Atomic radius is given in units of Angtroms.

Usage: (atomic-radius 5) => 1.17"
(periodic-atomic-radii nZ))




(defn atomic-mass [nZ]
  "Gives the atomic mass for any element.  nZ is an atomic number. The output
will be a real value corresponding to the atomic radius of nZ.

Atomic mass is given in units of a.u.

Usage: (atomic-weight 5) => 10.811"
(periodic-atomic-mass nZ))


(defn atomic-charge [nZ]
  (periodic-charges nZ))



