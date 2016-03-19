
(use 'greenwood.math)
(use 'greenwood.xyz 'greenwood.atomic-structure-output )
(use 'greenwood.contrib-math :reload)
(use 'greenwood.mol :reload)
(use 'greenwood.supercell :reload)
(use 'greenwood.empirical-data)
(use 'graphitic :reload)
(use 'greenwood.neighbors)
(use '[clojure.string :only [join]])
(use 'reaxff)
;'





(def graphene "C	0.0	0.0	0
C	0.16666667	0.33333334	0
C	0.0	0.5	0
C	0.16666667	0.8333333	0
C	0.5	0.0	0
C	0.6666667	0.33333334	0
C	0.5	0.5	0
C	0.6666667	0.8333333	0")


(def C8F "C       -0.002324590  -0.004649180  -0.003953565
C        0.166666667   0.333333333   0.008185541
C       -0.002324590   0.502324590  -0.003953565
C        0.165337100   0.832668550  -0.004628464
C        0.500000000   0.000000000  -0.003569322
C        0.667331450   0.334662900  -0.004628464
C        0.504649180   0.502324590  -0.003953565
C        0.667331450   0.832668550  -0.004628464
F        0.166666667   0.333333333   0.071129871")
(def C8F_a 4.672288557475574)

(def C4F_05 "C        0.000000000   0.000000000   0.006587937
C        0.170383763   0.340767525  -0.007521498
C        0.007434192   0.503717096  -0.007521498
C        0.170383763   0.829616237  -0.007521498
C        0.496282904  -0.007434192  -0.007521498
C        0.666666667   0.333333333   0.006587937
C        0.496282904   0.503717096  -0.007521498
C        0.659232475   0.829616237  -0.007521498
F        0.000000000   0.000000000   0.065976559
F        0.666666667   0.333333333   0.065976559")
(def C4F_05_a 4.692215409311123)



(def C2F1256ud "C	-0.003445883	-0.006891767	-0.001616076
C	0.16874869	0.33749738	-0.0113299
C	-0.002081999	0.495836	0.011329786
C	0.1701125	0.840225	0.001616171
C	0.4965541	-0.006891805	-0.001616058
C	0.6687487	0.33749732	-0.011329875
C	0.49791798	0.495836	0.011329786
C	0.6701125	0.840225	0.001616171
F	0.1937216	0.3874432	-0.06832481
F	-0.027054835	0.44589034	0.0683248
F	0.6937215	0.38744298	-0.06832479
F	0.47294518	0.44589034	0.0683248")
(def C2F1256ud_a (* 0.5 9.475301725563062))



(def all "C        0.000000000   0.000000000   0.009868572
C        0.166666791   0.333333582  -0.009868544
F        0.000000000   0.000000000   0.066205232
F        0.166666285   0.333332570  -0.066205098
C        0.000000302   0.500000151   0.009868551
C        0.166666791   0.833333179  -0.009868544
F        0.000000387   0.500000194   0.066204873
F        0.166666285   0.833333685  -0.066205098
C        0.499999849  -0.000000302   0.009868551
C        0.666666700   0.333333340  -0.009868637
F        0.499999806  -0.000000387   0.066204873
F        0.666666700   0.333333340  -0.066204510
C        0.499999849   0.500000151   0.009868551
C        0.666666458   0.833333179  -0.009868544
F        0.499999806   0.500000194   0.066204873
F        0.666667470   0.833333685  -0.066205098")
(def all_a 4.924182376003504)







(defn multilayer-F-graphene
  "This is to make tri-layer graphene.  This allows for the top and bottom
layers to have different coverages.  Currently, I assume that
all three layers have the same a lattice constant."
  [top graphene ngraphene distance]
  (flatten [
    (shift top [0 0 (* ngraphene distance)])
    (map #(shift graphene [0 0 (* % distance)]) (range ngraphene))]))





(defn write-BIOGRF
[mol name xlat ylat]
(join "\n"
["XTLGRF 200"
(str "DESCRP " name )
(write-crystx (ffirst lvs) xlat ylat 400 90.00 90.00 90.00)
(write-reac-HETATM mol)
"END"]))


(defn make-structure
[second-g glat size ngraphene job-name]
(let [nxn 2
alat (Bohr->Angstrom (* nxn glat))
xlat (* alat (/ 3 (Math/sqrt 3)))
lvs [(a-one alat) (a-three alat) (last (QE-to-lvs glat nxn))]
CF (xyz-str->atoms second-g)
CF-sc (shift (atom-pos (create-surface-cartesian CF [0 0 0] lvs (a-one alat) [0 0 10] (* size nxn xlat) (* size nxn alat))) [0 0 10])
C (xyz-str->atoms graphene)
C-sc (shift (atom-pos (create-surface-cartesian C [0 0 0] lvs (a-one alat) [0 0 10] (* size nxn xlat) (* size nxn alat))) [0 0 10])]
(write-BIOGRF  (shift (shift->scell (* size nxn xlat) (* size nxn alat) (atom-pos (multilayer-F-graphene CF-sc C-sc ngraphene 3.2))) [0 0 150])  job-name (* size nxn xlat) (* size nxn alat))))



(defn Fgraphene-layer-multiples
[second-g Fglat glat size]
(let [nxn 2
alat (Bohr->Angstrom (* nxn glat))
xlat (* alat (/ 3 (Math/sqrt 3)))
lvs [(a-one alat) (a-three alat) (last (QE-to-lvs glat nxn))]
Falat (Bohr->Angstrom (* nxn Fglat))
Fxlat (* Falat (/ 3 (Math/sqrt 3)))
Flvs [(a-one Falat) (a-three Falat) (last (QE-to-lvs Fglat nxn))]
CF (xyz-str->atoms second-g)
CF-sc (create-surface-cartesian CF [0 0 0] lvs (a-one Falat) [0 0 10] (* size nxn Fxlat) (* size nxn Falat))
C (xyz-str->atoms graphene)
C-sc (create-surface-cartesian C [0 0 0] lvs (a-one alat) [0 0 10] (* size nxn xlat) (* size nxn alat))]
[(:lvs CF-sc) (:lvs C-sc)(incompatible-layers (:lvs CF-sc) (:lvs C-sc))]))










(def gmols [ C4F_05] )
(def gas [ C4F_05_a ])
(def gnames [ "C4F_05" ])


(Fgraphene-layer-multiples C4F_05 C4F_05_a (Angstrom->Bohr 2.461) 0.5)





