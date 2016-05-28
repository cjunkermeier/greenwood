(require '[greenwood.xyz :as xyz] '[greenwood.atomic-structure-output :as out])
(require '[greenwood.mol :as gmol]'[greenwood.math :as gmath])
(require '[greenwood.supercell :as gsc]'[greenwood.empirical-data :as ed])
(require '[graphitic :as g])
(require '[clojure.string :as cstr])
(require '[greenwood.basics :as b])
(require '[greenwood.neighbors :as gneigh])
(use 'ultra-csv.core)


(def twotwoF1 "C        0.000000000  -0.019879553  -0.090245310
C        0.000000002   1.427513666   0.211046501
C       -1.253478390   2.151208580  -0.090245196
C       -1.240177126   3.566523916  -0.105519069
C        2.472525600  -0.000000453  -0.077745281
C        2.472525600   1.432033515  -0.105519082
C        1.253478390   2.151208580  -0.090245196
C        1.240177126   3.566523916  -0.105519069
F        0.000000002   1.427513137   1.756184276")

(def LVStwotwoF1  [[4.94505096047861	0	0]
[-2.472525480239305	4.282539754783114	0]
[0	0	24.6073061235]])


(def twotwoF05 "C        0.000000000   0.000000295   0.175940229
C        0.000000002   1.465975006  -0.170451016
C       -1.213498056   2.166590064  -0.170442344
C       -1.213498739   3.567816216  -0.170451147
C        2.483070600  -0.032374122  -0.170442220
C        2.483070600   1.433601280   0.175944914
C        1.213498056   2.166590064  -0.170442344
C        1.213498739   3.567816216  -0.170451147
F        0.000000000  -0.000000222   1.637591834
F        2.483070600   1.433601698   1.637588521")

(def LVStwotwoF05 [[4.966141117175168	0	0]
[-2.483070558587584	4.300804366252128	0]
[0	0	24.6073061235]])







(def C8F3u014 "C       -0.003270478  -0.094092911   0.190908204
C       -0.003206177   1.486516452   0.183158777
C       -1.276518200   2.177501941  -0.289604443
C       -1.252534071   3.571648814  -0.221888761
C        2.508330207  -0.012800502   0.095837268
C        2.508326457   1.461641860  -0.289951997
C        1.270048794   2.177479885  -0.289848640
C        1.246054124   3.571661315  -0.222149059
F       -0.003530246  -0.464936570   1.611147860
F       -0.002987043   1.837129216   1.574358028
F        2.508424729   0.163559051   1.564610293")

(def LVSC8F3u014 [[5.023102708049032	0	0]
[-2.511551354024516	4.350134550988869	0]
[0	0	24.6073061235]])




(def C2F1256 "C       -0.000000003  -0.019166891  -0.300529286
C        0.000000000   1.364862614   0.403355085
C       -1.255775687   2.260249634   0.403355036
C       -1.255775678   3.644279075  -0.300529341
C        2.511551400  -0.019166967  -0.300529426
C        2.511551400   1.364862586   0.403354848
C        1.255775687   2.260249634   0.403355036
C        1.255775678   3.644279075  -0.300529341
F        0.000000000   0.817473204   1.720243946
F       -1.255775767   2.807639031   1.720243827
F        2.511551400   0.817472939   1.720243600
F        1.255775767   2.807639031   1.720243827")

(def LVSC2F1256  [[5.023102708049032	0	0]
[-2.511551354024516	4.350134550988869	0]
[0	0	24.6073061235]])




(defn structure-characteristics
[second-g lvs adsorb-species name]
(let [CCC (xyz/xyz-str->atoms second-g)
      CC (gmol/mol-filter {:species "C"} CCC)
      CC-sc (as-> (gsc/create-supercell CC (gsc/computation-projectors lvs 1 1 0)) x
                  (xyz/atom-pos x)
                  (gneigh/neighbors x 0.2 1.9))
     A (as-> (gneigh/angles CC-sc) x
             (take 8 x)
             (map (comp  #(map :radians %) :angles) x)
             (flatten x)
             (map ed/radians->degrees x)
             (gmath/simple-stats x))
     B (as-> CC-sc x
             (take 8 x)
             (map (comp  #(map :ndistance %) :neigh) x)
             (flatten x)
             (gmath/simple-stats x))
     CFB (as-> CCC x
             (gsc/create-supercell x (gsc/computation-projectors lvs 1 1 0))
             (xyz/atom-pos x)
             (gneigh/neighbors x 0.2 1.9)
             (take (count CCC) x)
             (gmol/mol-filter {:species  adsorb-species} x)
             (map (comp #(map :ndistance %) :neigh) x)
             (flatten x)
             (gmath/simple-stats x))]
    (hash-map :b-mean (:mean B) :b-var (:var B) :a-mean (:mean A) :a-var (:var A) :Ab-mean (:mean CFB) :Ab-var (:var CFB) :name name )))




(ultra-csv.core/write-csv!  "/Users/chadjunkermeier/Desktop/Fgraphene.csv" [
(structure-characteristics twotwoF1 LVStwotwoF1  "F" "C8F")
(structure-characteristics twotwoF05 LVStwotwoF05  "F" "C4F")
(structure-characteristics C8F3u014 LVSC8F3u014  "F" "C8F3")
(structure-characteristics C2F1256 LVSC2F1256  "F"  "C2F")])







