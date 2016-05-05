(require '[greenwood.xyz :as xyz] '[greenwood.atomic-structure-output :as out])
(require '[greenwood.mol :as gmol])
(require '[greenwood.math :as gmath])
(require '[greenwood.supercell :as gsc])
(require '[graphitic :as g])
(require '[clojure.string :as cstr])
(require '[greenwood.empirical-data :as ged])
(use 'reaxff)

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











 (defn get-data
 [a]
 (let [
 mols (as-> (xyz/xyz-str->atoms  "C        0.000000000   0.000000295   0.176677042
C        0.000000002   1.465801589  -0.171048541
C       -1.213648240   2.166503355  -0.171039869
C       -1.213648923   3.567902925  -0.171048672
C        2.483070600  -0.032200705  -0.171039745
C        2.483070600   1.433601280   0.176681727
C        1.213648240   2.166503355  -0.171039869
C        1.213648923   3.567902925  -0.171048672
F        0.000000000  -0.000000222   1.638647596
F        2.483070600   1.433601698   1.638644283") x
          (gsc/supercell x  [[4.966141117175168	0	0]
[-2.483070558587584	4.300804366252128	0]
[0	0	24.6073061235]] 2  2  1))]
  mols))




(def n (get-data 5.313061365159999))




(spit "/Users/chadjunkermeier/Desktop/graphene.xyz"
       (write-xyz  n ))















