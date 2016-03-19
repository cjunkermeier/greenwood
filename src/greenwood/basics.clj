(ns greenwood.basics)



(defrecord atoms ^{:doc "Must use (import 'JMD.mol.atoms) in any ns that you use this in"}
   [^String species ^doubles coordinates charge neigh angles name ^int pos])


(defrecord neigh-rec [^int npos ^String nspecies ^double ndistance ^doubles ncoord])
(defrecord charge-rec  [^doubles orbitals ^double total])

(defn new-atom [ species coordinates charge neigh angles name  pos]
  (->atoms species coordinates charge neigh angles name pos))


(defn neigh-struct [npos nspecies ndistance ncoord]
  (->neigh-rec  npos  nspecies ndistance ncoord))


(defn charge-struct [orbitals total]
  (->charge-rec orbitals total))


(defrecord cell ^{:doc "Must use (import 'JMD.crystals.unitcell) in any ns that you use this in"}
   [lvs mol])


(defn unitcell [lvs mol]
  (->cell lvs mol))


(defrecord angle-rec [apos radians])

(defn angle-struct [apos  radians]
  (->angle-rec apos radians))


(defrecord system-rec [^String name ^int time-step lvs mol])
(defn system [name time-step lvs mol]
  (->system-rec name time-step lvs mol))




