(use 'greenwood.xyz :reload)
(use 'greenwood.math :reload)
(require '[greenwood.mol :as gmol])
(use 'greenwood.neighbors-species :reload)
(:refer-clojure :exclude [* - + == /])
(use 'clojure.core.matrix)
(use 'clojure.core.matrix.operators)
(require '[clojure.core.reducers :as r])
(use 'ultra-csv.core)
(require '[clojure.string :as cstr])



(defn total_and_average_bond_order
  [b puc]
  (let [Fbondorder (as-> (:mol puc) x
                     (gmol/mol-filter {:species "F"} x)
                     (map :charge x)
                      (average x))
        adsorbedCbondorder (as-> (:mol puc) x
                             (neighbors-maxdis x  b)
                             (gmol/mol-filter {:species "C", :neigh #(some (fn [y](= "F" y)) (map :nspecies %))} x)
                             (map :charge x)
                             [ x (average x)])
        notadsorbedCbondorder (as-> (:mol puc) x
                             (neighbors-maxdis x  b)
                             (gmol/mol-filter {:species "C", :neigh #(not-any? (fn [y](= "F" y)) (map :nspecies %))} x)
                              (if (empty? x) [() "NA"]
                                   (let [y (map :charge x)]
                                     [y (average y)])))
        Name (first (cstr/split (:name puc) #"N"))
        N (read-string (second (cstr/split (:name puc) #"N")))]
    (hash-map   :SYSTEM Name :N N
                    :mean.Fbo Fbondorder
                    :mean.aCbo (second adsorbedCbondorder)
                    :mean.naCbo (second notadsorbedCbondorder)
                    :mean.Cbo (average (flatten [(first adsorbedCbondorder)(first notadsorbedCbondorder)])))))


(defn total_and_average_bond_order
  [b puc]
  (let [notadsorbedCbondorder (as-> (:mol puc) x
                             (neighbors-maxdis x  b)
                             (gmol/mol-filter {:species "C", :neigh #(not-any? (fn [y](= "F" y)) (map :nspecies %))} x)
                              (map :pos x))]
    notadsorbedCbondorder))








#_(->> (foldable-chunks "/Users/chadjunkermeier/Desktop/Cdisplacement_energy/xmolout")
     (r/map (partial parse-xmolout 5))
     (r/map (partial total_and_average_bond_order 2.6))
     (into [])
     (take-nth 2)
     (ultra-csv.core/write-csv! "/Users/chadjunkermeier/Desktop/Cdisplacement_energy/bond-order.csv"))

#_(def ind (->> (index-xyz "/Users/chadjunkermeier/Desktop/Cdisplacement_energy/xmolout")
     (take-last 4 )
     (take-nth 2)))

#_(->> (foldable-chunks "/Users/chadjunkermeier/Desktop/Cdisplacement_energy/xmolout" ind)
     (r/map (partial parse-xmolout 5))
     (r/map (partial total_and_average_bond_order 2.6))
     (into [])
)

(* 1180300 1154)

(take-nth 2 [0 1 2 3])
