(ns greenwood.core
  (:gen-class  :main true))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Jejoon's count OH and H2O
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;(use 'greenwood.xyz :reload)
;(use 'greenwood.atomic-structure-output :reload)
;(use 'greenwood.neighbors :reload)
;(require '[greenwood.utils :as utils])
;(use 'greenwood.mol :reload)
;(require '[clojure.core.reducers :as r])
;(require '[clojure.java.io :as io])

;(:refer-clojure :exclude [* - + == /])
;(use 'clojure.core.matrix)
;(use 'clojure.core.matrix.operators)
;(use 'greenwood.basics :reload)



#_(defn count-OH-H2O
  [mol]
  (->> mol
      (mol-filter {:species "O", :neigh #(> 3 (count %) 0)})
      (map :neigh)
       (map #(map :nspecies %))
       (map #(sort (into ["O"] %)) )
       (map #(clojure.string/join "" %) )
       (frequencies)
       (#(str  (% "HO") "  "  (% "HHO") "\n"))))




#_(defn get-JEJOON-OH-H2O-data
  [in-file indx out-file]
  (append-file out-file "HO  HHO\n")
  (->> (foldable-chunks in-file indx)
         (r/map (partial drop 2))
         (r/map xyz-iota->atoms )
          (r/map #(neighbors % 0.2 1.0))
       (r/map count-OH-H2O)
       (r/map #(append-file out-file %))
       (into [])))





#_(defn -main
  "I don't do a whole lot."
  ([infile  outfile]
   (let [indx (index-xyz infile)]
    (get-JEJOON-OH-H2O-data infile indx outfile)))
  ([infile  outfile natoms  nlines]
   (let [indx (natoms-index natoms  nlines)]
    (get-JEJOON-OH-H2O-data infile indx outfile))))




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Determine adsorption time step, temperature, strain, and movie
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(use 'greenwood.xyz)
;(use 'desorption)

#_(defn -main
  "This will output a file with the set of maps."
   [path]
   (let [NNN (N-folders path)]
     (map #(append-file (str path "/DESORPTION.dat") (str (within-N %) "\n"))  NNN)))


















