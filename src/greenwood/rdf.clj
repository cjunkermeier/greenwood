(ns greenwood.rdf
  (:require [greenwood.utils :as gu]
            [greenwood.neighbors :as gn]
            [clojure.string :as cstrng]
       ))

(defn hist-bin-edges [rmin rmax nbins]
 "This defines a col of the points that act as delimiters of the histogram bins.
Usage: (hist-bin-edges 0 10 10) => (0 1 2 3 4 5 6 7 8 9 10)"
 (let [step-size (/ (- rmax rmin) nbins)]
   (doall (map #(+ rmin (* step-size %)) (range (inc nbins))))))


(defn bin-points [edges lcr]
  "This can be used to define the points at which the binned data is to be plotted.
By setting lcr to :r the point will be on the right hand edge of the bins.  :c
stands for center of the bin, and :l the LHS.

edges is a 1-D col listing where the bin min and max are.

This computes the mid-point of each bin in case the user wants to use some weird
bin system that has different sizes of bins.

Usage: (bin-pointss (hist-bin-edges 0 3 3) :c) +> (1/2 3/2 5/2)"
  (condp = lcr
    :l (drop-last edges)
    :c (map #(+ (/ (- %2 %1) 2.0) %1) edges (next edges))
    :r (drop 1 edges)))



(defn binned-data [sorted-data edges]
  "sorted-data is a col of values, edges is a col of histogram bin edges.

Usage: (binned-data [1 2 3 4 5 6 7 8 9 10] [1.5 2.5 6.5 7.6])
=> ((2) (3 4 5 6) (7))"
  (map
    (fn [x]
      (take-while #(> (second x) %) (drop-while #(> (first x) %) sorted-data)))
    (partition 2 1 edges)))




(defn- pair-distribution-binning-
  "This computes the radial pair distribution function centered on one
particular atom, called atomm.  It is called by pair-distribution."
  [mol atomm edges]
  (let [dist (sort (gn/distances mol atomm))]
    (map count (binned-data dist edges))))



(defn sph-shell-volume [bmin bmax]
  "This defines the volume of a spherical shell.  This can be used to define the dV of a histogram bin
in the radial-pair-distribution function or even the volume of the whole thing."
  (* 4.18879020479 (- (* bmax bmax bmax) (* bmin bmin bmin))))



(defn rdf
  "This computes the pair distribution function or radial distribution function of a mol."
  [edges mol]
  (let [V (sph-shell-volume (first edges)(last edges))
        dV  (map sph-shell-volume edges (next edges))
        in-holes (map #(reduce + %)(gu/transpose (pmap #(doall (pair-distribution-binning- mol % edges)) mol)))
        tot-num (reduce + in-holes)
        tot-density (/ tot-num V)]
    (map #(hash-map :r %1 :y %2) (bin-points edges :c) (map #(/ %1 (* %2 tot-density)) in-holes dV))))




(defn distribution->plot-txt
  "This is the function to use if you are going to write out and plot (in your favorite program) the distribution due to one time step."
  [distribution]
  (gu/inter-cat-tree [gu/endline ", "] (map #(vector (:r %) (:y %)) distribution)))


(defn distribution->append-file
  "This is the function to use if you are going to write out and plot (in your favorite program) the distribution due to one time step."
  [filename distribution]
  (gu/append-file filename  (cstrng/join ", " (map :y distribution))))



;(require '[clojure.core.reducers :as r] '[greenwood.xyz :as xyz])
#_(->> (foldable-chunks "/Users/chadjunkermeier/Desktop/graphene.xyz" )
     (r/map (partial drop 2))
     (r/map xyz-iota->atoms)
     (r/map (partial rdf (hist-bin-edges 1 10 18)))
     (r/map #(distribution->append-file "/Users/chadjunkermeier/Desktop/graphene.txt" %))
     (into []))


#_(inccore/with-data (inccore/to-dataset (rdf (hist-bin-edges 1 10 36) graphene))
    (inccore/view (inccharts/bar-chart :r :y
                     :title "CO2 Uptake"

                     :x-label "Grass Types" :y-label "Uptake"
                    :legend true)))



