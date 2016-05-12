(require '[greenwood.xyz :as xyz] '[greenwood.atomic-structure-output :as out])
(require '[greenwood.mol :as gmol])
(require '[greenwood.math :as gmath])
(require '[greenwood.supercell :as gsc])
(require '[graphitic :as g])
(require '[clojure.string :as cstr])
(require '[greenwood.empirical-data :as ged])
(use 'reaxff)


(def twotwoF1 (xyz/xyz-str->atoms "C       -0.072156043   0.026875558   0.321729674
C       -0.076023351   1.615815774   0.312268832
C       -1.344371317   2.288794863  -0.205319164
C       -1.322463763   3.684856708  -0.142421093
C        2.424729138   0.106661783   0.221532047
C        2.424842294   1.577469908  -0.187579834
C        1.193466751   2.293180372  -0.187520108
C        1.172637254   3.685299977  -0.122914876
O       -0.071280800  -0.383336599   1.785059949
H        0.896969102  -0.257217842   2.0372136065"))
(def twotwoF1_lvs (QE-to-lvs 4.6846611209931845 2))














