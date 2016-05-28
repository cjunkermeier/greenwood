(use 'greenwood.utils)
(use 'greenwood.materials-analysis)


(def EOS (transpose (partition 2 [4.7    -157.3809507681 
4.8    -157.4430301543 
4.9    -157.4759926435 
5.0    -157.4806714261  
5.1    -157.4606814228 
5.2    -157.4191637762])))


(println (second (Birch-Murnaghan-cohesive-properties (first EOS) (second EOS )))) 

