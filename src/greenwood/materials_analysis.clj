(ns JMD.materials-analysis
    (:refer-clojure :exclude [* - + == /])
  (:use
        [filters :only [multi-filter]]
        [clojure.core.matrix]
        [clojure.core.matrix.operators])
  (:require
   [greenwood.empirical-data :as ged]
   [greenwood.math :as gmath]
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


(defn Birch-Murnaghan-cohesive-properties [vol-col en-col]
  "This is used to compute the cohesive properties of a material using the
Birch-Murnaghan equation. The input is a col of lattice volumes (lattice constants),
vol-col, for the material along with the total energy, en-col, associated with
vol-col.  The output is a col containing, in order, the equilibrium total energy,
equilibrium lattice volume, bulk modulus, and the derivative of the bulk modulus
with respect to pressure.  This can be used with just a col of lattice constants
in place of lattice volumes, but you will no longer obtain the bulk modulus.

If the output ends up being (NaN NaN NaN NaN) you could have exchanged the positioning
of the volume and energy vectors in the call of this function."
  (let [min-en (apply min en-col)
        where (first (gutil/positions #(= min-en %) en-col))
        min-lat (nth vol-col where)]
    (:coefs(incopt/non-linear-model Birch-Murnaghan-eqn en-col vol-col [min-en min-lat 4 2]))))



#_(defn Birch-Murnaghan-cohesive-properties [vol-col en-col]
  "This is used to compute the cohesive properties of a material using the
Birch-Murnaghan equation. The input is a col of lattice volumes (lattice constants),
vol-col, for the material along with the total energy, en-col, associated with
vol-col.  The output is a col containing, in order, the equilibrium total energy,
equilibrium lattice volume, bulk modulus, and the derivative of the bulk modulus
with respect to pressure.  This can be used with just a col of lattice constants
in place of lattice volumes, but you will no longer obtain the bulk modulus.

If the output ends up being (NaN NaN NaN NaN) you could have exchanged the positioning
of the volume and energy vectors in the call of this function."
  (let [min-en (apply min en-col)
        where (first (gutil/positions #(= min-en %) en-col))
        min-lat (nth vol-col where)]
    (:coefs(incopt/non-linear-model Birch-Murnaghan-eqn en-col vol-col [min-en min-lat 11 4] :newton-raphson true))))



(defn graphene-resonator-frequency
  "This computes the frequency of a graphene resonator with a gate voltage of zero.
The form for this equation comes from Nature Nanotechnology 4, 861 - 867 (2009),
http://www.nature.com/nnano/journal/v4/n12/full/nnano.2009.267.html.
This needs the mol and lvs to determine the mass density.  In this we will
actually calculate the f/f-ref, where L,w,To are respectively equal to L-ref,
w-ref, To-ref, and Te = 0."
  [mol-ref lvs-ref mol lvs]
  (let [mass  #(reduce + (map (comp ged/atomic-mass ged/atomic-numbers :species) %))]
    (sqrt (/
      (/ (mass mol-ref) (gmath/lvs-volume lvs-ref))
      (/ (mass mol) (gmath/lvs-volume lvs))))))







(defn find-band-gap
  "Used with quantum espresso"
  [bands fermienergy]
  (let [above-below (multi-filter [#(> % fermienergy) #(< % fermienergy)] (flatten bands))]
    (- (apply min (first above-below)) (apply max (second above-below)))))





; Third order Birch-M. EOS http://en.wikipedia.org/wiki/Birch–Murnaghan_equation_of_state
;http://gilgamesh.cheme.cmu.edu/doc/software/jacapo/appendices/appendix-eos.html
;https://wiki.fysik.dtu.dk/ase/ase/utils.html
#_(defmacro polynomial-eqn
"it is assumed that t is the wines root of the volume of the cell."
  [theta t]
  (let [[c0 c1 c2 c3 c4 c5 c6 t] theta]
(inccore/plus c0 (inccore/mult c1 t) (inccore/mult c2 t t)
   (inccore/mult c3 t t t)
    (inccore/mult c4 t t t t)
    (inccore/mult c5 t t t t t)
    (inccore/mult c6 t t t t t))))



#_(defn polynomial-cohesive-properties [vol-col en-col]
  "This is used to compute the cohesive properties of a material using a
polynomial equation. The input is a col of lattice volumes (lattice constants),
vol-col, for the material along with the total energy, en-col, associated with
vol-col.  The output is a col containing, in order, the equilibrium total energy,
equilibrium lattice volume, bulk modulus, and the derivative of the bulk modulus
with respect to pressure.  This can be used with just a col of lattice constants
in place of lattice volumes, but you will no longer obtain the bulk modulus.

If the output ends up being (NaN NaN NaN NaN) you could have exchanged the potitioning
of the volume and energy vectors in the call of this function."
  (let [min-en (apply min en-col)
        where (first (gutil/positions #(= min-en %) en-col))
        min-lat (nth vol-col where)]
    (:coefs(incopt/non-linear-model polynomial-eqn en-col vol-col [min-en min-lat 4 2]))))





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



