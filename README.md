# greenwood
A library for creating molecular models and processing molecular dynamics simulations.


Better documentation is fourthcoming.

##Example of how to create a box of 35 DCAH molecules at a density of 0.25 g/ml.

(require '[greenwood.solution :as gsol])
(require '[greenwood.xyz :as gxyz])
(require '[greenwood.atomic-structure-output :as gout])
(require '[clojure.string :as cstr])
(require '[greenwood.utils :as gutils])



(def DCAH (gxyz/xyz-str->atoms-readable "N   12.32170  13.00473  11.78595
C   13.29902  13.45354  12.74665
N   14.10215  13.81901  13.53892
C   11.73234  11.80858  12.05561
N   11.18361  10.64494  12.34951
H   10.63120   9.81316  12.43579"))

(def v 25.0)
(def lvs [[v 0 0] [0 v 0] [0 0 v]])
;(gsol/sol-density lvs [] (vector (take 35 (repeat DCAH))))
(def gas-phase-single-component-DCAH (gsol/solution [] lvs [0 0 0] 35 DCAH))
(spit "/Users/chadjunkermeier/Desktop/gas-phase-single-component-DCAH.xyz" (gout/write-xyz gas-phase-single-component-DCAH))


