
(ns desorption)

(require '[clojure.string :as cstrng])
(use 'greenwood.xyz :reload)
(use 'greenwood.atomic-structure-output :reload)
;(use 'greenwood.neighbors :reload)
(use 'greenwood.neighbors-species :reload)
(require '[greenwood.utils :as utils])
(use 'greenwood.mol :reload)
(require '[clojure.core.reducers :as r])
(require '[clojure.java.io :as io])

(:refer-clojure :exclude [* - + == /])
(require '[clojure.core.matrix :as cmat])
(use 'clojure.core.matrix.operators)
(use 'greenwood.basics :reload)
(require '[me.raynes.fs :as fs])
(use 'reaxff)
(use 'ultra-csv.core)
(use 'greenwood.math :reload)


(require '[clojure.string :as strng])





(defn desorption?
  "Determines if a fluorine atom has desorbed from graphene."
  [mol]
  (let [cmol (mol-filter {:species "C"} mol)
        fmol (mol-filter {:species "F"} mol)
        d #(distances cmol %)]
    (loop [ff (d (first fmol))
           ffmol (rest fmol)]
      (cond
        (every? #(>  % 4.0) ff) true
        (empty? ffmol) 1
        :else (recur (d (first ffmol)) (rest ffmol))))))





(defn determine-desorption-time
  [filename]
  (->> (foldable-chunks filename)
         (r/map (partial drop 2))
         (r/map (partial xyz-iota->atoms ))
         (r/map desorption?)
         (r/take-while #(= 1 %))
       (r/reduce + 0 )
       (#(* % 100) )))





(defn determine-desorption-time
  [filename rearrange-ts]
    (let [natoms (-> (iota/seq filename)
              (first)
              (read-string))
          nlines (->> (iota/seq filename)
                     (r/map (fn [x](let [a x] 1)) )
                     (r/fold +))
        total-indx (natoms-index natoms nlines)
          rts (+ -100 (if (re-matches #"[0-9]+50" (str rearrange-ts))
                          (- rearrange-ts 50)
                          rearrange-ts))
                rearrang-indx (count (reax-index-timesteps 0 rts 100 natoms))
        indx (drop  rearrang-indx total-indx)]
  (->> (foldable-chunks filename indx)
         (r/map (partial drop 2))
         (r/map (partial xyz-iota->atoms ))
         (r/map desorption?)
         (r/take-while #(= 1 %))
       (r/reduce + rearrang-indx )
       (#(* % 100) ))))








;(determine-desorption-time "/Users/chadjunkermeier/Desktop/graphene.xyz" 2)


(defn Desorption-strain
  [mol]
    (->> (mol-filter {:species "C"}  mol)
         (map :charge )
         (map cmat/abs)
         (reduce +  0.0 )))





(defn movement-time
  "movement could be rearrangement or desorption."
  [filename]
  (let [n (-> (iota/seq filename)
              (first)
              (read-string))
        timestep #(* 50  %)
        ]
    (->> (iota/seq filename)
              (map read-string)
              (take-while #(= n %))
           (count )
           (timestep))))



(defn which-F?
  "Determines if a fluorine atom has desorbed from graphene."
  [mol]
  (let [cmol (mol-filter {:species "C"} mol)
        fmol (mol-filter {:species "F"} mol)
        d #(distances cmol %)]
    (loop [ff (first fmol)
           ffmol (rest fmol)]
      (cond
        (every? #(>  % 4.0) (d ff)) (:pos ff)
        :else (recur (first ffmol) (rest ffmol))))))




(defn get-mol
  [path timestep]
  (let [natoms (-> (iota/seq (str path "/xmolout"))
              (first)
              (read-string))
        indx (reax-index-timesteps timestep timestep 100 natoms)]
    (->> (foldable-chunks (str path "/xmolout") indx)
         (r/map (partial drop 2))
         (r/map (partial xyz-iota->atoms 5))
         (into [])
         (first ))))


(defn N-folders
  "This determines which N[1-13] folders are included in the results."
  [path]
  (utils/grep #"N[1-9]" (fs/list-dir path )))


(defn within-N
  [path]
  (let [Rtime (movement-time (str path "/percent_coverage.txt"))
        Dtime (determine-desorption-time (str path "/xmolout") Rtime) ]
    (hash-map :system (first (take-last 2 (strng/split (str path) #"/")))
              :N (last  (strng/split (str path) #"/"))
              :rearrange Rtime
              :desorption Dtime
              :F (which-F?  (get-mol (str path) Dtime))
              :strain (Desorption-strain (get-mol path (- Dtime 100))))))






(defn get-desorption-data
  "This will output a file with the set of maps."
   [NNN]
     (map #(utils/append-file (str % "/DESORPTION.dat") (str (within-N %) "\n"))  NNN))




(defn- write-BIOGRF-
[name puc]
(cstrng/join "\n"
["XTLGRF 200"
(str "DESCRP " name  )
 "RUTYPE SINGLE POINT"
(write-crystx ((comp ffirst :lvs) puc) ((comp second second :lvs) puc)  ((comp last last :lvs) puc) 90 90 90)
(write-reac-HETATM (:mol puc))
"END"
 "\n"]))






(defn create-singpoint-geo-records
  "We use this to create the geo files."
  [name path]
    (let [rearrange-ts (movement-time (str path "/percent_coverage.txt"))
          natoms (-> (iota/seq (str path "/xmolout"))
              (first)
              (read-string))
          rts (+ -100 (if (re-matches #"[0-9]+50" (str rearrange-ts))
                          (- rearrange-ts 50)
                          (- rearrange-ts 100)))
          indx (reax-index-timesteps rts rts 100 natoms)]
  (->> (foldable-chunks (str path "/xmolout") indx)
         (r/map parse-xmolout)
         (into [])
         (first)
         (write-BIOGRF- name))))





(defn total_and_average_bond_order
  [puc]
  (let [Fbondorder (as-> (:mol puc) x
                     (mol-filter {:species "F"} x)
                     (map :charge x)
                      (average x))
        adsorbedCbondorder (as-> (:mol puc) x
                             (neighbors-maxdis x  1.8)
                             (mol-filter {:species "C", :neigh #(some (fn [y](= "F" y)) (map :nspecies %))} x)
                             (map :charge x)
                             [ x (average x)])
        notadsorbedCbondorder (as-> (:mol puc) x
                             (neighbors-maxdis x  1.8)
                             (mol-filter {:species "C", :neigh #(not-any? (fn [y](= "F" y)) (map :nspecies %))} x)
                              (if (empty? x) [() "NA"]
                                   (let [y (map :charge x)]
                                     [y (average y)])))
 ]
    (hash-map   :mean.Fbo Fbondorder
                    :mean.aCbo (second adsorbedCbondorder)
                    :mean.naCbo (second notadsorbedCbondorder)
                    :mean.Cbo (average (flatten [(first adsorbedCbondorder)(first notadsorbedCbondorder)])))))








#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C32F0MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C32F0MPaN" %) (str "/Volumes/HAWAII/DESORPTION/RAMP/SC32F/N" % )) [1 2 3 4 5 6 7 8  10 11 12 13])))


#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C32F100MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C32F100MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C32F100Mpa/N" % )) [2 3 4 5 6 7 8 9 10 11 12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C32F300MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C32F300MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C32F300Mpa/N" % )) [1 2 3 4 5 6 7 8 9 10 11 12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C32F500MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C32F500MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C32F500Mpa/N" % )) [1  3 4 5 6 7 8 9  12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C32F700MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C32F700MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C32F700Mpa/N" % )) [1 2 3 4 5 6 7 8 9 10 11 12 13])))




;;;;;
#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C4F0MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C4F0MPaN" %) (str "/Volumes/HAWAII/DESORPTION/RAMP/SincreaseC4F/N" % )) [1 2 3 4  6 7 8 9 10 11 12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C4F100MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C4F100MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C4F100MPa/N" % )) [1 2 3 4 5 6 7 8 9 10 11 12 13])))


#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C4F300MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C4F300MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C4F300MPa/N" % )) [1 2 3 4 5 6 7 8 9 10 11 12 13])))


#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C4F500MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C4F500MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C4F500MPa/N" % )) [1 2 3 4 5 6  8 9 10 11 12 13])))




#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C2F1256ud0MPa.geo"
 (cstrng/join ""
 (map #(create-singpoint-geo-records (str "C2F0MPaN" %) (str "/Volumes/HAWAII/DESORPTION/RAMP/C2F1256ud/N" % )) [1 2 3 4 5 6 7 8 9 10 11 12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C2F1256ud100MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C2F100MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C2F1256ud100f/N" % )) [1 2 3 4 5 6 7 8 9 10  ])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/C2F1256ud300MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "C2F300MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/C2F1256ud300MPa/N" % )) [1 2 3 4 5 6 7 8 9 10  12 13])))







#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/all0MPa.geo"
 (cstrng/join ""
 (map #(create-singpoint-geo-records (str "CF0MPaN" %) (str "/Volumes/HAWAII/DESORPTION/RAMP/Sall/N" % )) [1 2 3 4 5 6 7  12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/all100MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "CF100MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/all100MPa/N" %)) [2 5 6 7 8 9  12 13])))

#_(spit "/Users/chadjunkermeier/Desktop/SINGLEPOINTGEOS/all300MPa.geo"
 (cstrng/join "\n"
 (map #(create-singpoint-geo-records (str "CF300MPaN" %) (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/all300MPa/N" % )) [2  4  6 7 10 11 12 ])))







#_(->> (foldable-chunks "/Users/chadjunkermeier/Desktop/xmolout.xyz" )
(r/map (partial drop 2))
     (r/map xyz-iota->atoms)
     (r/map #(shift [0 0 -200] %))
     (r/map write-xyz)
    (r/map #(append-file "/Users/chadjunkermeier/Desktop/graphene.xyz" %))
     (into []))


;(def g (hash-map :desorption 1519300, :system "all300MPa", :F 662, :strain 29941.37983999998, :N "N2", :rearrange 1212950))

(defn parse-energy-log
  [path hm]
  (utils/clean-parse (re-pattern (:rearrange hm)) (str path (:N hm))))



(defn parse-energy-log
  [ hm]
  (let [rts (if (re-matches #"[0-9]+50" (str (:rearrange hm)))
                          (- (:rearrange hm) 50)
                          (- (:rearrange hm) 100))]
  (->> (iota/seq (str "/Volumes/HAWAII/DESORPTION/NEGATIVEPRESSURE/" (:system hm) "/" (:N hm) "/energylog"))
       (r/filter (partial re-find (re-pattern (str rts))) )
       (into [])
       (first )
       (str (:system hm) "  " (:N hm) "  "))))



;(parse-energy-log  g)


;(get-desorption-data ["/Volumes/Untitled/C4F_05_14"])

