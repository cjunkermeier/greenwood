(require '[greenwood.xyz :as xyz] '[greenwood.atomic-structure-output :as out])
(require '[greenwood.mol :as gmol])
(require '[greenwood.supercell :as gsc])
(require '[graphitic :as g])
(require '[clojure.string :as cstr])
(require '[greenwood.basics :as b])

(use 'graphitic :reload)


;Here are the list of the 2x2 super cell configurations.
;	    		4
;	             \
;                 \
;   0              5
;	 \            /
;     \          /
;      1--------6
;     /          \
;    /            \
;   2              7
;    \
;     \
;      3

(def OHOreo [
  [[[0 1]] [2]]   ;done
  [[[0 1] [2 3]] [5 6]]
  [[[0 1] [4 5]] [2 6]]
  [[[0 1] [5 6]] [2 7]]
  [[[0 1] [6 7]] [2 5]]])




(def names ["O01OH2" "O0123OH56" "O0145OH26" "O0156OH27" "O0167OH25"])







(defn make-structure
[O OH name]
(let [g (g/graphene-QE-unit-cell "C" "C" 2.461)
        C-sc (gsc/supercell (:mol g) (:lvs g)  2 2 1)
      GO  (loop [sc (xyz/atom-pos (:mol C-sc))
             b (first O)
             bleft (rest O)]
            (if (nil? b)
                (as-> sc x
                      (gmol/shift [0 0 5] x)
                      (assoc C-sc :mol  x))
                (recur
                 (g/bond-centered-adsorption sc b [(b/new-atom "O" [0.0 0.0 0.0] nil nil nil nil 999)] 0 1.22)
                (first bleft)
                (rest bleft))))
    GOH (assoc GO :mol (g/top-site-adsorption (:mol GO) [(b/new-atom "O" [0.0 0.0 0.0] nil nil nil nil 0)(b/new-atom "H" [0.0 0.0 1.02] nil nil nil nil 999)] OH 1.47))
    ]
  (str
  (out/write-xyz-cshell (:mol GOH)  5000 name)
   "\n\n"
  (out/write-lvs-cshell (:lvs GOH) name)
   "\n\n")))








(spit "/Users/chadjunkermeier/Desktop/sandra_epoxy.xyz"
      (cstr/join "\n" (map  #(make-structure (first %1)(second %1) %2) OHOreo names)))






