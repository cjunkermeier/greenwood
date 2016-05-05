# Introduction to greenwood

TODO: write [great documentation](http://jacobian.org/writing/what-to-write/)





;(require 'reaxff)
(require '[greenwood.solution :as gsol])
(require '[greenwood.xyz :as gxyz])
(require '[greenwood.atomic-structure-output :as gout])
(require '[clojure.string :as cstr])
(require '[greenwood.utils :as gutils])
;(require '[greenwood.empirical-data :as ged])
;(require '[greenwood.neighbors :as gn])



(def DCAH (gxyz/xyz-str->atoms-readable "N   12.32170  13.00473  11.78595
C   13.29902  13.45354  12.74665
N   14.10215  13.81901  13.53892
C   11.73234  11.80858  12.05561
N   11.18361  10.64494  12.34951
H   10.63120   9.81316  12.43579"))





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;gas phase DCAH, rhp=0.25, n=35
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(def v 25.0)
(def lvs [[v 0 0] [0 v 0] [0 0 v]])
;(gsol/sol-density lvs [] (vector (take 35 (repeat DCAH))))
(def gas-phase-single-component-DCAH (gsol/solution [] lvs [0 0 0] 35 DCAH))
(spit "/Users/chadjunkermeier/Desktop/gas-phase-single-component-DCAH.xyz" (gout/write-xyz gas-phase-single-component-DCAH))


