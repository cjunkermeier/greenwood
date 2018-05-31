
(ns find-nodes
  (:require [greenwood.utils :as utils]
            [greenwood.mol :as gmol]
            [ubergraph.core :as uber]
            [ubergraph.alg :as alg]
            [ubergraph.protocols :as prots])
(:refer-clojure :exclude [* - + == /])
  (:use clojure.core.matrix)
  (:use clojure.core.matrix.operators))
(use '[incanter.core :only [to-dataset to-list]] 'incanter.excel)
(require 'ultra-csv.core )


(defn cell-repeat
  "This creates a repeating cell (like a supercell of atoms in atomic coordinates)
of the displacement in A-B space.  This is used when the path we want to determine
goes outside of the unit cell."
  [la1 ma2 data]
    (cond (and (zero? la1) (zero? ma2)) data
           (and (pos? la1) (zero? ma2)) (concat data (gmol/update-mol :A #(+ % la1) data))
           (and (zero? la1) (pos? ma2)) (concat data (gmol/update-mol :B #(+ % ma2) data))
           (and (pos? la1) (pos? ma2)) (concat data (gmol/update-mol :A #(+ % la1) data)
                                                    (gmol/update-mol :B #(+ % ma2) data)
                                        (gmol/update-mol :B #(+ % ma2) (gmol/update-mol :A #(+ % la1) data)))))



(defn parse-graph-data
  [file la1 ma2]
  (->>  (utils/parse-table-hash file)
        (cell-repeat la1 ma2)
        (sort-by (juxt :A :B) )
      (map #(assoc %2 :pos %1) (iterate inc 0))))


(defn vertices-coords
  [hm]
  (map #(assoc % :coordinates [(:A %) (:B %)]) hm))


(defn create-ubernodes
  [hm]
  (map  #(vector ((comp keyword str) (:pos %)) %) hm))


(defn endpoint-pos
  ""
  [hm A B]
 ((comp keyword str :pos first)
  (gmol/mol-filter {:A A, :B B} hm)))



(defn- find-edges-
  "This finds the cost of moving th"
  [hs dmin dmax vertex]
  (loop [a (first hs)
         b (rest hs)
         edges '()]
    (cond (empty? a)
             edges
          (>= dmax (distance (:coordinates a) (:coordinates vertex)) dmin)
            (recur (first b) (rest b)
             (conj edges
             (vector ((comp keyword str :pos)  vertex) ((comp keyword str :pos) a)
                       (abs (- (:Energy a) (:Energy vertex)))))) ;cost given by slope
          :else
              (recur (first b) (rest b) edges))))






(defn create-uberedges
    [hm dmin dmax]
  (let [h (vertices-coords hm)]
  (utils/flatten-n 1 (map (partial find-edges- h dmin dmax) h))))


(defn create-ubergraph
  [hm dmin dmax]
  (apply
      (partial uber/add-directed-edges (apply uber/multigraph (create-ubernodes hm)))
   (create-uberedges hm dmin dmax)))


(defn path-nodes
"This will transform the path into a col of hash-maps that can be used to
transform the data into a format that is useful for writing it out into a
csv file."
  [g p]
    (let [edges (map (comp read-string name)
                     (flatten [(map uber/src (prots/edges-in-path p))
                               (uber/dest (last (prots/edges-in-path p)))]))
          nodes (uber/nodes g)]
    (->> (map (partial uber/attrs g) nodes)
         (gmol/mol-filter-vec :pos  edges )
         (map #(dissoc % :pos) ))))





(defn energy-surface-paths
  [file la1 ma2 dmin dmax A1 B1 A2 B2]
  (let [d (parse-graph-data file la1 ma2)
        p1 (endpoint-pos d A1 B1)
        p2 (endpoint-pos d A2 B2)
        ug (create-ubergraph d dmin dmax)]
(path-nodes ug
(alg/shortest-path ug {:start-node p1, :end-node p2, :cost-attr :weight}))))










;(save-xls (to-dataset bilayer-graphene-path1) "/Users/junky/Documents/PennState/bilayer_BPC/bilayer_graphene/bilayer_graphene_path1.xls")



(def graphite-CC "/Users/junky/Documents/PennState/bilayer_BPC/bilayer_graphene/new_graphite_graph.txt")
(def bulk-graphite-path1 (energy-surface-paths graphite-CC 0 0 0.0001  0.015 0.33 0.33 0.67 0.67))
(ultra-csv.core/write-csv!  "/Users/junky/Documents/PennState/bilayer_BPC/bilayer_graphene/bulk_graphite_path1.csv" bulk-graphite-path1)





