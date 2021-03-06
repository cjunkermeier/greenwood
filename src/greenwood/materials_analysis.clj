(ns greenwood.materials-analysis
    (:refer-clojure :exclude [* - + == /])
  (:use
        [filters :only [multi-filter]]
        [clojure.core.matrix]
        [clojure.core.matrix.operators])
  (:require
   [greenwood.empirical-data :as ged]
   [greenwood.math :as gmath]
   [greenwood.mol :as gmol]
   [greenwood.neighbors :as gn]
   [greenwood.utils  :as gutil]
        [incanter.core :as inccore ]
[incanter.optimize :as incopt]))




(defn Birch-Murnaghan-eqn [theta V]
  "Equation for the Birch-Murnaghan equation of state.  This was derived for a 3-D
cubic crystal structures but has been used for many other structures.
For explanation of physics:
http://www.fhi-berlin.mpg.de/th/Meetings/DFT-workshop-Berlin2009/scientificprogram.html
http://www.fhi-berlin.mpg.de/th/Meetings/DFT-workshop-Berlin2009/Talks/OnlinePublication/0624-4T_20090701-1_Yoon_2perA4_-_aims_workshop_problem_set.pdf
http://www.fhi-berlin.mpg.de/th/Meetings/DFT-workshop-Berlin2009/Talks/OnlinePublication/0630-4T_20090630-1_FH_-_phonon_tutorial_script_web.pdf

Units of Results:
Enot = the input unit of energy (probably eV)
Vnot = the input unit of volume (when I did this for graphene it was Angstroms)
Bnot = Enot/Vnot
Bnot_prime = unitless"
  (let [[Enot Vnot Bnot Bnot-prime] theta]
    (inccore/plus Enot
      (inccore/mult (inccore/div (inccore/mult Bnot V) Bnot-prime)(inccore/plus 1 (inccore/div (inccore/pow (inccore/div Vnot V) Bnot-prime) (inccore/minus Bnot-prime 1))))
      (inccore/div (inccore/mult Bnot Vnot)(inccore/minus 1 Bnot-prime)))))


(defn Birch-Murnaghan-EOS
  "This is used to compute the cohesive properties of a material using the
Birch-Murnaghan equation. The input is a col of lattice volumes (lattice constants),
vol-col, for the material along with the total energy, en-col, associated with
vol-col.  The output is a col containing, in order, the equilibrium total energy,
equilibrium lattice volume, bulk modulus, and the derivative of the bulk modulus
with respect to pressure.  This can be used with just a col of lattice constants
in place of lattice volumes, but you will no longer obtain the bulk modulus.

If the output ends up being (NaN NaN NaN NaN) you could have exchanged the positioning
of the volume and energy vectors in the call of this function.

The user while the values of bnot and bnotprime will generally work, I have found that
  at times the user may need to try a few different values to get a result that converges."
  ([vol-col en-col]
  (let [min-en (apply min en-col)
        where (first (gutil/positions #(= min-en %) en-col))
        min-lat (nth vol-col where)]
   (#(hash-map :Enot (first (:coefs %)) :Vnot (second (:coefs %))
               :Bnot (nth (:coefs %) 2) :Bnot-prime (nth (:coefs %) 3)
               :rss (:rss %))
     (incopt/non-linear-model Birch-Murnaghan-eqn en-col vol-col [min-en min-lat 4 2]))))
  ([vol-col en-col bnot bnotprime]
  (let [min-en (apply min en-col)
        where (first (gutil/positions #(= min-en %) en-col))
        min-lat (nth vol-col where)]
   (#(hash-map :Enot (first (:coefs %)) :Vnot (second (:coefs %))
               :Bnot (nth (:coefs %) 2) :Bnot-prime (nth (:coefs %) 3)
               :rss (:rss %))
     (incopt/non-linear-model Birch-Murnaghan-eqn en-col vol-col [min-en min-lat bnot bnotprime])))))


#_(def g (transpose [[4.55  -139.8214545264]
[4.6  -139.8389697975]
[4.625  -139.8442017475]
[4.65  -139.8471905565]
[4.66  -139.8477729292]
[4.67  -139.8480217128]
[4.675  -139.8480189493]
[4.68  -139.8479320836]
[4.69  -139.8475103507]
[4.700  -139.8467658516]
[4.725  -139.8435146113]
[4.750  -139.8383433422]
[4.775  -139.8313303156]
[4.800  -139.8225552678]]))
#_(Birch-Murnaghan-EOS h en  4 2)


(defn polynomial-eqn
"A simple polynomial equation."
  [theta v]
  (let [[v0 c0 c1 c2 c3 c4] theta]
              (inccore/plus c0
              (inccore/mult c1 (- v v0))
              (inccore/mult c2 (- v v0) (- v v0))
              (inccore/mult c3 (- v v0) (- v v0) (- v v0))
              (inccore/mult c4 (- v v0) (- v v0) (- v v0) (- v v0)))))




(defn polynomial-EOS
  "This is used to compute the cohesive properties of a material using a
polynomial equation. The input is a col of lattice volumes (lattice constants),
vol-col, for the material along with the total energy, en-col, associated with
vol-col.  The output is a col containing, in order, the equilibrium total energy,
equilibrium lattice volume, bulk modulus, and the derivative of the bulk modulus
with respect to pressure.  This can be used with just a col of lattice constants
in place of lattice volumes, but you will no longer obtain the bulk modulus.

If the output ends up being (NaN NaN NaN NaN) you could have exchanged the potitioning
of the volume and energy vectors in the call of this function."
   ([vol-col en-col]
    (polynomial-EOS vol-col en-col [0.00001 1.0 1 1]))
([vol-col en-col theta]
 (let [[c1 c2 c3 c4] theta
        min-en (apply min en-col)
     where (first (gutil/positions #(= min-en %) en-col))
     min-v (nth vol-col where)]
 (#(hash-map :Enot (second (:coefs %)) :Vnot (first (:coefs %))
            :rss (:rss %))
  (incopt/non-linear-model polynomial-eqn en-col vol-col [min-v min-en c1 c2 c3 c4])))))



(defn young-modulus
  "This formulation is taken from the paper:
  'Equilibrium configuration and continuum elastic properties of finite sized graphene'
  by C D Reddy, S Rajendran, and K M Liew.  Appearing in Nanotechnology 17 (2006) 864–870
  (http://iopscience.iop.org/0957-4484/17/3/042).
  Specifically, we use equations 7 and 11 from this paper."
   [volumes energies]
    (let [min-en (apply min energies)
        where (first (gutil/positions #(= min-en %) energies))
        min-v (nth volumes where)
        model (incopt/non-linear-model polynomial-eqn energies volumes [min-v min-en 0.00001 1.0 1 1])]
      (hash-map :ym (/ (nth (:coefs model) 2) (first (:coefs model))) :rss (:rss model) )))


(defn young-modulus2
  "This formulation is taken from the paper:
  'Equilibrium configuration and continuum elastic properties of finite sized graphene'
  by C D Reddy, S Rajendran, and K M Liew.  Appearing in Nanotechnology 17 (2006) 864–870
  (http://iopscience.iop.org/0957-4484/17/3/042).
  Specifically, we use equations 7 and 11 from this paper."
   ([strains energies]
    (let [min-en (apply min energies)
        where (first (gutil/positions #(= min-en %) energies))
        min-s (nth strains where)
        model (incopt/non-linear-model polynomial-eqn energies strains [min-s min-en 0.00001 1.0 1 1])]
      (hash-map :ym (/ (nth (:coefs model) 2) (first (:coefs model))) :rss (:rss model) )))
   ([strains energies theta]
    (let [[c1 c2 c3 c4] theta
        min-en (apply min energies)
        where (first (gutil/positions #(= min-en %) energies))
        min-s (nth strains where)
        model (incopt/non-linear-model polynomial-eqn energies strains [min-s min-en c1 c2 c3 c4])]
      (hash-map :ym (/ (nth (:coefs model) 2) (first (:coefs model))) :rss (:rss model) ))))








(defn local-lindemann-index
  "The Lindemann index effectively measures the degree of translation
  freedom of a molecule and in this sense is similar to the mean square
  displacement (MSD). The Lindemann index measures the average displacement
  of molecules relative to their neighbors, whereas the MSD measures the
  average displacement of a molecule relative to its own position at
  an earlier time. So in principle, just as the magnitude of the MSD can
  distinguish between a solid and liquid, or a liquid and gas, the Lindemann
  index should be able to do the same.

  The boiling point depends on the pressure, but at a given pressure, the
  magnitude of the Lindemann index should still be able to distinguish
  between 'a translationally free - small volume' liquid phase and a
  'translationally free - large volume' gas phase.
  The advantage of using the Lindemann index is that it does not depend as
  much on the size of the simulation system. For a simulation with periodic
  boundary conditions in all three directions, the Lindemann index and MSD
  give similar information, but in the case of a nanoparticle simulation, or
  a simulation of a confined system, the MSD will depend on system dimension
  and not give a clear picture of the molecular motion.
  --Saman Alavi · University of British Columbia - Vancouver
  The Lindemann index[1] is a simple measure of thermally driven disorder in
  atoms or molecules. In condensed matter physics a departure from linearity in
  the behaviour of the global Lindemann index or an increase above a threshold
  value related to the spacing between atoms (or micelles, particles, globules,
  etc.) is often taken as the indication that a solid-liquid phase transition
  has taken place.  THIS FORMULATION IS USED FOR NON-PERIODIC SYSTEMS.
  This function outputs a seq containing the lindemann-index for each atom in
  the system."
  [timesteps]
  (let [ntimesteps (/ 1. (count timesteps))
        Nminus1 (/ 1. (dec (count (first timesteps))))
        distances (pmap gn/all-distances timesteps)
        d (reduce + distances)
        d2 (reduce + (map #(* % %) distances))
        f #(/ (sqrt (- (* ntimesteps %1) (* ntimesteps ntimesteps %2 %2))) (* ntimesteps %2))]
           (* Nminus1 (reduce + (transpose (map f d2 d))))))


(defn rolling-average-lindemann-index
  "The Lindemann index[1] is a simple measure of thermally driven disorder in atoms or molecules.
  In condensed matter physics a departure from linearity in the behaviour of the global Lindemann
  index or an increase above a threshold value related to the spacing between atoms (or micelles,
  particles, globules, etc.) is often taken as the indication that a solid-liquid phase transition
  has taken place.  This function outputs a seq containing the average of the local-lindemann-index
  for n time steps in the system.
  This uses the partition function to break up all of the timesteps into small chunks so that you
  can analyze the progression of the lindemann index over the course of the calculation."
  [timesteps n m]
  (map (comp #(/ (reduce + %) (count %)) local-lindemann-index) (partition n m timesteps)))




(defn POAV1-Pyramidalization-Angle
  "The π-orbital axis vector (POAV) analysis provides a complete description
  of the electronic structure of nonplanar conjugated organic molecules, and
  the current interest in fullerenes and carbon nanotubes has led to widespread
  applica- tion of the POAV method. As may be inferred from the name, the
  method is based on vector algebra, and the equations necessary to solve for
  the various quantities of interest (pyramidalization angles, dihedral angles,
  hybridizations, and resonance integrals) are best handled with a PC and the
  computer program POAV3.11 The only quantities necessary for this analysis are
  the atomic coordinates of the atoms in the molecule or fragment. Thus, to
  solve for the pyramidalization angle and hybridization at a single nonplanar
  conjugated carbon merely requires the atomic coordinates of the conjugated atom
  and its three attached atoms (denoted 1, 2, and 3). --Haddon, R. C., 'Comment on
  the Relationship of the Pyramidalization Angle at a Conjugated Carbon
  Atom to the sigma Bond Angles.' J. Phys. Chem. A 2001, 105, 4164-4165.

  Usage: (POAV1-Pyramidalization-Angle coronene 1 2 3 4)

  "
  [mol conjugatedC C1 C2 C3]
  (let [B11 (cos (gmol/mol-angle mol [C1 conjugatedC C2]))
        B12 (cos (gmol/mol-angle mol [C2 conjugatedC C3]))
        B13 (cos (gmol/mol-angle mol [C3 conjugatedC C1]))
        B14 (sin (gmol/mol-angle mol [C1 conjugatedC C2]))
        B20 (/ (- B12 (* B11 B13)) B14)
        B21 (sqrt (- 1 (* B13 B13) (* B20 B20)))
        B23 (* -1 B14 B21)
        B24 (* B14 B21 B14 B21)
        B25 (* (- 1 B11) (- 1 B11) B21 B21)
        B26 (** (- (* (- 1 B11) B20) (* B14 (- 1 B13))) 2)
        B30 (/ B23 (sqrt (+ B24 B25 B26)))      ;cos(theta_sigma.pi)
        B31 (/ (- B23) (sqrt (+ B24 B25 B26)))  ;cos(theta_sigma.pi)
        B35 (ged/radians->degrees (acos B30))   ;theta_sigma.pi
        B36 (ged/radians->degrees (acos B31))]  ;theta_sigma.pi
      (- B35 90)))
