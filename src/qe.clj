(ns qe
 (require [clojure.core.matrix :as cmtx])
 (use     clojure.core.matrix.operators))



(defn cubic-fcc
[a]
 (* 0.5 a [[-1 0 1][0 1 1][-1 1 0]]))




(defn triclinic
  [a b c alpha beta gamma]
  [[a 0.0 0.0]
   [(* b cos(gamma)), (* b sin(gamma)), 0.0]
   [(* c cos(beta)) , (/ (* c (- cos(alpha) (* cos(beta) cos(gamma)))) sin(gamma)),
    (/ (* c (sqrt (+ 1  (* 2 cos(alpha) cos(beta) cos(gamma))
                      (* -1 cos(alpha) cos(alpha))
                      (* -1 cos(beta) cos(beta))
                      (* -1 cos(gamma) cos(gamma))))) sin(gamma))]])


(defn QE-triclinic
  [celldm1 celldm2 celldm3 celldm4 celldm4 celldm5 celldm6]
  (triclinic celldm1 (* celldm2 celldm1) (* celldm3 celldm1) celldm4 celldm5 celldm6))





