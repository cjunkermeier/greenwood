 (ns greenwood.math
  (:require ;[clojure.core.reducers :as r]
               ;[clojure.set :as cset]
            [greenwood.contrib-math :as cm]
            [clojure.core.matrix :as cmat]
            [clojure.core.matrix.operators :as cmato]
            [clojure.core.matrix.stats :as cms]
            [clojure.core.matrix.random :as cmr])
  (:refer-clojure :exclude [* - + == /]))



;(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(def ^:private ^:const  pi 3.141592653589793)

(def ^:private ^:const tau 6.283185307179586)

(defn degrees
  "Transforms the angle a from radians to degrees."
  ^double [a]
  (cmato/* 0.31830988618379  180.  a))

(defn radians
  "Transforms the angle a from degrees to radians."
 ^double  [a]
  (cmato/* pi  0.555555555555556  a))

(defn scalar-times-vector
  "multiplies each element of a vector by the value of scalar."
  [scalar vect]
  (cmato/* scalar vect))

(defn dot-product
  ^double [x y]
  (cmat/dot x y))


(defn cross-product
  "pnt1 and pnt2 are 3-tuples."
  [pnt1 pnt2]
  (cmat/cross pnt1 pnt2))




(defn magnitude
  "Calculates the sqrt of the sum of the squares of coords"
  [coords]
  (cmat/length coords))



(defn unit-vec [vect]
  "Creates a unit vector in the direction of the input vect."
  (cmat/normalise vect))

(defn euclidean [x y]
  "Returns the euclidean distance between two coordinates, represented as seqs of numbers."
  (cmat/distance x y))


(defn mat-vect-mult [mat vect]
  "Performs matrix vector multplication.
   mat must be a seq of row vectors, not column vectors"
  (map #(reduce cmato/+ (map cmato/* % vect)) mat))


(defn mat-mat-mult [mat1 mat2]
  "Performs matrix matrix multiplication.
   mat1 and mat2 are seqs of row vectors, not column vectors."
  (let [t (cmat/transpose mat2)]
    (map #(mat-vect-mult t %) mat1)))



(defn average
  "computes the average value of a col"
  ^double [v]
  (cms/mean v))





(defn median
  "computes the median value of a col"
  [v]
  (let [c (cmato// (count v) 2)]
  (if (integer? c)
    ((comp cms/mean (partial take 2) (partial drop (dec (cmat/floor c))) sort) v)
    ((comp first (partial drop (dec c)) sort) v))))


(defn quantiles
"Computes the median and upper and lower quartiles a seq of values.  Outputs the a hash-map."
  [v]
  (let [m (median v)]
  (hash-map :upper  (median (filter  #(> % m) v))  :lower  (median (filter  #(< % m) v)) :median m)))



(defn sum-sqrs [v]
  (cms/sum-of-squares v))

(defn p-sd
"Computes the population standard deviation of a vector of values.
  If you collect data from all members of a population or set, you apply the population standard deviation."
  [v]
 (cmat/sqrt (cmato// (sum-sqrs (cmato/- v (cms/mean v))) (count v))))


(defn s-sd
"Computes the sample standard deviation of a vector of values.
  If you collect data from all members of a population or set, you apply the population standard deviation."
  [v]
  (if (cmato/== (cmato/- (apply max (map cmato/- v))) (apply min v))
    0
  (cms/sd v)))


(defn simple-stats
"Computes the mean, median, standard deviation, and upper and lower quartiles a seq of values.  Outputs the a hash-map."
 [v]
  (if (cmato/== (cmato/- (apply max (map cmato/- v))) (apply min v))
    (hash-map :upper nil :lower nil :median (median v) :mean (cms/mean v) :var 0)
    (assoc (quantiles v) :mean (cms/mean v) :var (cms/variance v))))





(defn midpoint
  "Returns the midpoint between points x and y."
  [x y]
   (cmato/+ x (cmato/* (cmat/length (cmato/- y x)) 0.5 (cmat/normalise (cmato/- y x)))))




(defn abs-angle
  "Given an angle, ang, this will transform it to being a positive value in
the y>0 region of the 2 pi sphere."
  ^double [^double ang]
  (let [a (mod (cmat/abs ang) tau)]
  (cond
    (= pi a) pi
    (< pi a) (cmato/- tau a)
    (> pi a) a)))




(defn normal-vector
  "Given two vectors that lie on the plane, and are not parallel, this will
find the vector that is normal to the plane.
Usage:  (normal-vector [1 0 0] [0 1 0]) => [0.0 0.0 1.0]"
  [v1 v2]
  (let [descrimant (cmat/sqrt (cmato/- (cmato/* (cmat/length v1) (cmat/length v2)) (cmat/abs (cmat/dot v1 v2))))]
    (cmato/* (cmato// 1.0 descrimant) (cmat/cross v1 v2))))



(defn normalize-sum [col]
  "Divides each element of col by the sum of col, such that the new sum is equal to 1."
  (let [sum (reduce cmato/+ col)]
    (if (= 0.0 sum)
      (let [c (count col)]
        (take c (repeat (cmato// 1.0 c))))
      (cmato/* (cmato// 1.0 sum) col))))



(defn find-angle
  "Returns the angle between the vectors a and b."
  [a b]
  (cmat/acos (cmato// (cmat/dot a b) (cmat/length a) (cmat/length b))))




(defn three-point-angle
  "This computes the angle described the three points a P b, where P is
  the pivot between a and b."
  ^double [a P b]
  (let [A (cmato/- a P)
        B (cmato/- b P)]
    (find-angle A B)))



(defn bisect-angle-midpoint
  "Given three points that define the angle A-P-C; this finds the unit vector
that points from P to the midpoint of the line connecting A and C."
  [A P C]
(cmat/normalise (cmato/- (midpoint A C) P)))



(defn point-line-distance
  "This comes from mathworld.wolfram.com.  This finds the shortest distance
between a point and a line."
  ^double [end1 end2 point]
  (cmato//
    (cmat/length (cmat/cross (cmato/- point end1) (cmato/- point end2)))
    (cmat/length (cmato/- end1 end2))))



(defn point-line-intersection
  "Returns the point on a line (defined by the segment end1 and end2)
   which is closest to another point. all args are 3-tuples of nums"
   [end1 end2 point]
  (let [l (unit-vec (cmato/- end2 end1))
	      h (cmato/- point end1)
	      theta (find-angle l h)]
    (cmato/+ end1 (cmato/*  (cmat/length h) (cmat/cos theta) l))))



(defn point-plane-distance
  "In a three-dimensional space one way of defining a plane is by specifying a point and a normal vector to the plane.
Let pl-point be the position vector of some known point in the plane, let n be a nonzero vector normal to the plane
and p be the point off of the plane.
  The resultant value is positive if p is on the same side of the plane as the normal vector and negative if it is on the opposite side."
  [n pl-point p]
  (cmat/dot n (cmato/- p pl-point)))



(defmacro tolerance? [value tol]
  "This function is a two-tuple predicate.
   It is TRUE if abs(val) < abs(tol), FALSE otherwise.

Usage: (tolerance? 2 1.0E-8) => false
       (tolerance? 1.0E-10 1.0E-8) => true"
   `(< (cmat/abs ~value) (cmat/abs ~tol)))



(defn setzero  [ value]
  "This sets a value less than 1.0E-8 to zero."
  (if (tolerance? value 1.0E-8) 0 value))


(defn exclude-zero
  "Usage: (exclude-zero [0.00000000001 1 2 3 4 5 6]) => (1 2 3 4 5 6)"
  [col]
  (filter (comp not zero? setzero) col))


(defn origin?
  [vec1]
  (every? zero? (map setzero vec1)))


(defn round-decimal
  "This rounds a real number at 1E-6."
  ([^double real-num]
  (-> real-num
    (cmato/* 1E6M)
    (cmat/round )
    (cmato/* 1E-6M)
    (double )))
  ([^double val ^double real-num ]
  (-> real-num
    (cmato// val)
    (cmat/round )
    (cmato/* (bigdec val))
    (double ))))




(defn lvs->parameters
  "Converts a unit cell specified in lattice vectors to lattice parameters.
   Assumes the a vector is aligned with the cartesian x axis, and the b vector is in the cartesian xy plane."
  [lattice]
  (let [[a b c] lattice]
    (concat (map cmat/length lattice)
	    (map find-angle [b c a] [c a b]))))


(defn lvs-volume
  "Given a set of three unit vectors this computes the volume of the unit cell."
  [lvs]
  (cmat/dot (first lvs) (cmat/cross (second lvs) (last lvs))))



(defn vectors-equal?
"Determines if two vectors are equal, meaning that they have the same number of
elements and that the difference between each corresponding element is less than
1.0E-6.
Usage: (vectors-equal? [1 1 1][1 1 1]) -> true
Usage: (vectors-equal? [1 1 1][1 1 1 1]) -> false
Usage: (vectors-equal? [1 1 1][1 1 0.9999]) -> false
Usage: (vectors-equal? [1 1 1][1 1 0.9999999]) -> true"
  ([vec1 vec2]
   (vectors-equal? vec1 vec2 1.0E-6))
  ([vec1 vec2 tol]
  (let [num-args (= (count vec1) (count vec2))
        differences (cmato/- vec1 vec2)]
    (every? true? (cons num-args (map #(tolerance? % tol) differences))))))





(defn- rot-arb-axis-not-parallel-z [point pnt1 pnt2 angle]
  "This algorithm assumes that the axis of rotation is not parallel to the z-axis.

We must give the axis an orientation so that positive and negative angles of
rotation are defined. If the axis of rotation is given by two points
pnt1 = (a,b,c) and pnt2 = (d,e,f), then the (oriented) vector of rotation can be
given by (u,v,w) = (d-a,e-b,f-c).

This was derived at http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/

point is a 3-tuple containing the point you want to rotate around the axis that
is defined by the line (map - pnt2 pnt1), where pnt2 and pnt1 are both 3-tuples.
angle is the amount of ration this is a scalar value.

Usage: (rot-arb-axis-not-parallel-z [0 0 1] [1 0 0] [2 0 0] (/ Math/PI 2)) => [0.0 1.0 6.123233995736766E-17]"
  (let [axis (cmat/normalise (cmato/- pnt2 pnt1))
        x (first point)
        y (second point)
        z (nth point 2)
        a (first pnt1)
        b (second pnt1)
        c (nth pnt1 2)
        u (first axis)
        v (second axis)
        w (nth axis 2)
        A (sum-sqrs axis)
        B (sum-sqrs [v w])
        G (sum-sqrs [u w])
        F (sum-sqrs [u v])
        magaxis (cmat/length axis)
        rmpa (cmat/dot point axis)
        ca (cmat/cos angle)
        sa (cmat/sin angle)]
 (vector
    (cmato/+ (cmato/* (cmato/- (cmato/* a B) (cmato/* u (cmato/- (cmato/* b v) (cmato/* -1 c w) rmpa))) (cmato/- 1 ca))
             (cmato/* x ca)
             (cmato/* (cmato/+ (cmato/* -1 c v) (cmato/* b w) (cmato/* -1 w y)  (cmato/* v z)) sa))
    (cmato/+ (cmato/* (cmato/- (cmato/* b G) (cmato/* v (cmato/- (cmato/* a u) (cmato/* -1 c w) rmpa))) (cmato/- 1 ca))
             (cmato/* y ca)
             (cmato/* (cmato/+ (cmato/* c u) (cmato/* -1 a w) (cmato/* w x)  (cmato/* -1 u z)) sa))
    (cmato/+ (cmato/* (cmato/- (cmato/* c F) (cmato/* w (cmato/- (cmato/* a u) (cmato/* -1 b v) rmpa))) (cmato/- 1 ca))
             (cmato/* z ca)
             (cmato/* (cmato/+ (cmato/* -1 b u) (cmato/* a v) (cmato/* -1 v x)  (cmato/* u y)) sa)))))




(defn xrotation [theta inputvec]
  (let [s (cmat/sin theta)
        c (cmat/cos theta)
        [i1 i2 i3] inputvec]
    [i1,
     (cmato/+ (cmato/* i2 c) (cmato/* i3 s)),
     (cmato/- (cmato/* i3 c) (cmato/* i2 s))]))

(defn yrotation [theta inputvec]
  (let [s (cmat/sin theta)
        c (cmat/cos theta)
        [i1 i2 i3] inputvec]
    [(cmato/- (cmato/* i1 c) (cmato/* i3 s)),
     i2,
     (cmato/+ (cmato/* i3 c) (cmato/* i1  s))]))

(defn zrotation [theta inputvec]
  (let [s (cmat/sin theta)
        c (cmat/cos theta)
        [i1 i2 i3] inputvec]
    [(cmato/-  (cmato/* i1  c) (cmato/* i2 s)),
     (cmato/+ (cmato/* i1 s) (cmato/* i2 c)),
     i3]))


(defn the-rotation-function [point pnt1 pnt2 angle]
  "This is really designed for rotating a bunch of vectors (ie atom coordinates)
  where the thing that is being rotated can't just be assumed to be centered at the origin.

  If you have something like a unit vector that can be moved to the origin without changing
  it's meaning then assume that it is at the origin and place pnt1 at the origin.

point is a 3-tuple containing the point you want to rotate around the axis that
is defined by the line (map - pnt2 pnt1), where pnt2 and pnt1 are both 3-tuples
that define the axis.
angle is the amount of ration this is a scalar value."
  (let [axis (cmato/- pnt2 pnt1)
        output (cond
                 (and (tolerance? (second axis) 1.0E-6) (tolerance? (last axis) 1.0E-6))
                 (cmato/+ pnt1 (xrotation angle (cmato/- point pnt1)))
                 (and (tolerance? (first axis) 1.0E-6) (tolerance? (last axis) 1.0E-6))
                 (cmato/+ pnt1 (yrotation angle (cmato/- point pnt1)))
                 (and (tolerance? (first axis) 1.0E-6) (tolerance? (second axis) 1.0E-6))
                 (cmato/+ pnt1 (zrotation angle (cmato/- point pnt1)))
                 :else
                 (rot-arb-axis-not-parallel-z point pnt1 pnt2 angle))]
    (map #(setzero %) output)))




(defn- rotate-vec-2-xaxis [u v w]
  "This produces a rotation matrix that will rotate a vector onto the x-axis.
This function is only a helper function and is desigened to only be used within
rotate-vec-to-axis within this namespace."
  (let [vw (cmat/length [v w])
        uvw (cmat/length [u v w])]
    (cond
      (and (tolerance? w 1.0E-6) (tolerance? v 1.0E-6))
      [[1.0 0.0 0.0][0.0 1.0 0.0][0.0 0.0 1.0]]
    :else [[(cmato// u uvw)   (cmato// v uvw) (cmato// w uvw)]
   [(cmato// vw uvw -1.0) (cmato// (cmato/* u v) vw uvw) (cmato// (cmato/* u w) vw uvw)]
   [0.0 (cmato// w vw -1.0)  (cmato// v vw)]])))


(defn- rotate-vec-2-yaxis [u v w]
  "This produces a rotation matrix that will rotate a vector onto the y-axis.
This function is only a helper function and is desigened to only be used within
rotate-vec-to-axis within this namespace."
  (let [uw (cmat/length [u w])
        uvw (cmat/length [u v w])]
    (cond
      (and (tolerance? w 1.0E-6) (tolerance? v 1.0E-6))
      [[1.0 0.0 0.0][0.0 1.0 0.0][0.0 0.0 1.0]]
    :else [[(cmato// u uw)  0.0 (cmato// w uw -1.0)]
   [(cmato// u uvw) (cmato// v uvw) (cmato// w uvw)]
   [(cmato// (cmato/* u v) uw uvw) (cmato// uw uvw -1.0)  (cmato// (cmato/* w v) uw uvw)]])))


(defn- rotate-vec-2-zaxis [u v w]
  "This produces a rotation matrix that will rotate a vector onto the z-axis.
This function is only a helper function and is desigened to only be used within
rotate-vec-to-axis within this namespace."
  (let [uv (cmat/length [u v])
        uvw (cmat/length [u v w])]
    (cond
      (and (tolerance? u 1.0E-6) (tolerance? v 1.0E-6))
      [[1.0 0.0 0.0][0.0 1.0 0.0][0.0 0.0 1.0]]
    :else [[(cmato// (cmato/* u w) uv uvw)   (cmato// (cmato/* w v) uv uvw)   (cmato// uv uvw -1.0)]
   [(cmato// v uv -1.0) (cmato// u uv) 0.0]
   [(cmato// u uvw) (cmato// v uvw)  (cmato// w uvw)]])))




(defn rotate-vec-to-axis [somevec axis]
  "This produces a rotation matrix that will rotate somevec, where somevec is a
3-tuple, onto the x-, y-, or z-axis.  The choice of which axis is made by
specifying axis, where if the user wants to rotate somevec onto the x-axis then
they should use the key :x.

This matrix was derived from R_{xz} and R_{xz2z} found at
  http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/"
  (let [x (first somevec)
        y (second somevec)
        z (last somevec)]
  (cond
    (= axis :x) (rotate-vec-2-xaxis x y z)
    (= axis :y) (rotate-vec-2-yaxis x y z)
    (= axis :z) (rotate-vec-2-zaxis x y z))))





(defn reflection-matrix
  "Reflects a vector through a plane.  Currently, this only allows for
reflections through the xy-, xz-, yz-planes.  In a future version this will allow
for greater flexability for choosing a plane.  To choose the xy-plane use the
keyword :xy, similarly for the xz- and yz-planes.

Usage (reflection-matrix :yz) => [[-1 0 0][0 0 0][0 0 0]]"
  [plane]
  (cond
    (= plane :xy) [[1.0 0.0 0.0][0.0 1.0 0.0][0.0 0.0 -1.0]]
    (= plane :yz) [[-1.0 0.0 0.0][0.0 1.0 0.0][0.0 0.0 1.0]]
    (= plane :xz) [[1.0 0.0 0.0][0.0 -1.0 0.0][0.0 0.0 1.0]]))




(defn linear-interpolation
  "Suppose I have a range of numbers [26.43534 234.24] and a value of 118.12 that
lies within that range, and I want to do a rescale the range so that it is [0 255]
and then determine what my value is within that range.  The technical terminology
for this is linear-interpolation.

Usage: (linear-interpolation [26.43534 234.24] [0 255] 118.12) => 112.50752653958772"
  [[initial-min initial-max] [target-min target-max] value]
  (let [scale (cmato//  (cmato/- target-max target-min) (cmato/- initial-max initial-min))]
    (cmato/+ target-min (cmato/* (cmato/- value initial-min) scale))))



(defn linear-interpolation-definite
  "Suppose I have a range of numbers [26.43534 234.24] and a value of 118.12 that
lies within that range, and I want to do a rescale the range so that it is [0 255]
and then determine what my value is within that range.  The technical terminology
for this is linear-interpolation.

Usage: (linear-interpolation-definite [0 5] [0 1] 6) => 1"
  [[initial-min initial-max] [target-min target-max] value]
  (let [scale (cmato//  (cmato/- target-max target-min) (cmato/- initial-max initial-min))
        v (cond
                      (< value initial-min) initial-min
                      (> value initial-max) initial-max
                      :else  value)]
    (cmato/+ target-min (cmato/* (cmato/- v initial-min) scale))))




(defn card-rot-mat
  "This is used to define a rotation matrix around one of the cardinal axes.
Usage: (card-rot-mat pi :z) =>   "
  [angle axis]
   (let [s (cmat/sin angle)
         c (cmat/cos angle)]
     (cond
       (= axis :x) [[1.0 0.0 0.0][0.0 c s][0.0 (cmato/- s) c]]
       (= axis :y) [[c 0.0 (cmato/* -1 s)][0.0 1.0 0.0][s 0.0 c]]
       (= axis :z) [[c (cmato/- s) 0.0][s c 0.0][0.0 0.0 1.0]])))


(defn dihedral
"If we think of a clock with b1 and b3 as the hands of the clock, with b2 as the
pivot, this computes the angle between b1 and b3."
  [b1 b2 b3]
  (-> (cm/atan2
    (cmat/dot (cmato/* (cmat/length b2) b1) (cmat/cross b2 b3))
    (cmat/dot (cmat/cross b1 b2) (cmat/cross b2 b3)))
    (cmato/+ pi)))




(defmacro tolerated-eq
  "Computers often lie to us about the precision of a number, mostly due to the
damn idiots who implemented the IEEE floating point specification."
  ([x y] `(tolerated-eq ~x ~y 0.00001))
  ([x y epsilon]
    `(<= (abs (cmato/- ~x ~y)) ~epsilon)))



(defmacro tolerated-lt
  "Computers often lie to us about the precision of a number. They will often
give us a real number like 12.305 as 12.304999999999998,
thus (< 12.304999999999998 12.305) = true.  But since 12.30499999999998 was
supposed to be 12.305 the inequality should have been false.  This function is
designed to correct this failure.

Usage: (tolerated-lt 12.304999999999998 12.305) => false
Usage: (tolerated-lt 12.304999999999998 12.306) => true"
  ([a b]
    `(tolerated-lt ~a ~b 10E-6))
  ([a b epsilon]
      `(and (< ~a ~b) (not (tolerance? (cmato/- ~a ~b) ~epsilon)))))

(defmacro tolerated-gt
  "Computers often lie to us about the precision of a number. They will often
give us a real number like 12.305 as 12.304999999999998,
thus (> 12.305 12.304999999999998) = true.  But since 12.30499999999998 was
supposed to be 12.305 the inequality should have been false.  This function is
designed to correct this failure.

Usage: (tolerated-gt 12.305 12.304999999999998) => false
Usage: (tolerated-gt 12.306 12.304999999999998) => true"
  ([a b]
    `(tolerated-gt ~a ~b 10E-6))
  ([a b epsilon]
      `(and (> ~a ~b) (not (tolerance? (cmato/- ~a ~b) ~epsilon)))))

(defmacro tolerated-lte
  "Computers often lie to us about the precision of a number. They will often
give us a real number like 12.305 as 12.304999999999998,
thus (< 12.304999999999998 12.305) = true.  But since 12.30499999999998 was
supposed to be 12.305 the inequality should have been false.  This function is
designed to correct this failure.

Usage: (tolerated-lte 12.304999999999998 12.305) => false
Usage: (tolerated-lte 12.304999999999998 12.306) => true"
  ([a b]
    ~(tolerated-lte ~a ~b 10E-6))
  ([a b epsilon]
      (or
        `(and (< ~a ~b) (not (tolerance? (cmato/- ~a ~b) ~epsilon)))
        (or (< ~a ~b) (tolerance? (cmato/- ~a ~b) ~epsilon)))))

(defmacro tolerated-gte
  "Computers often lie to us about the precision of a number. They will often
give us a real number like 12.305 as 12.304999999999998,
thus (< 12.304999999999998 12.305) = true.  But since 12.30499999999998 was
supposed to be 12.305 the inequality should have been false.  This function is
designed to correct this failure.

Usage: (tolerated-gte 12.304999999999998 12.305) => true
Usage: (tolerated-gte 12.304999999999998 12.306) => true"
  ([a b]
    `(tolerated-gte ~a ~b 10E-6))
  ([a b epsilon]
      `(or
        (and (> ~a ~b) (not (tolerance? (cmato/- ~a ~b) ~epsilon)))
        (or (> ~a ~b) (tolerance? (cmato/- ~a ~b) ~epsilon)))))



(defn random-point [mins maxes]
  "Generates a random point inside the specified boundaries."
  (map #(cmato/+ (cmato/* (cmato/- %2 %1) (rand)) %1)
       mins maxes))

(defn rand-point-fractional
  "Returns a random point [0..1] in fractional coordinates.
  Useage: (rand-point)"
   []
  (take 3 (repeatedly rand)))



(defn rand-displacement
  "Displaces a point, pnt, by a random amount where that amount is drawn
  from a normal distribution with a standard deviation of sd.

  Usage: (rand-displacement :xy 0.1 [1 0 0]) => (1.0498892020635808 -0.11030772423707187 0)
  Usage: (rand-displacement :xyz 0.1 [1 0 0]) => (1.0056703785638446 -0.0034408339313103174 0.002919629292373202)"
  [tgl sd pnt]
  (cond
    (= tgl :xyz) (the-rotation-function (cmato/+ pnt (cmat/abs [(cmato/* (first (cmr/sample-normal 1)) sd) 0.0 0.0])) pnt (rand-point-fractional) (rand tau))
    (= tgl :xy) (the-rotation-function (cmato/+ pnt (cmat/abs [(cmato/* (first (cmr/sample-normal 1)) sd) 0.0 0.0])) pnt (cmato/+ pnt [0.0 0.0 1.0]) (rand tau))
    (= tgl :z) (cmato/+ pnt (cmat/abs [0.0 0.0 (cmato/* (first (cmr/sample-normal 1)) sd)]))
    (= tgl :x) (cmato/+ pnt (cmat/abs [(cmato/* (first (cmr/sample-normal 1)) sd) 0.0 0.0]))
    (= tgl :y) (cmato/+ pnt (cmat/abs [0.0 (cmato/* (first (cmr/sample-normal 1)) sd)  0.0]))))




(defn intersecting-spheres?
  "Determines if a sphere of radius r1 intersections a sphere or radius 2
  when they are centered at p1 and p2, resptectively.

  Usage: (intersecting-spheres? 1 2 [0 0 0] [0 0 2.5]) => true
  Usage: (intersecting-spheres? 1 2 [0 0 0] [0 0 3.5]) => false"
  [r1 r2 p1 p2]
  (> (cmato/+ r1 r2) (cmat/distance p1 p2)))




(defn normalize-mag [lst]
  "Divides each element of lst by the magnitude of lst, such that the new magnitude is equal to 1."
  (let [mag (cmat/length lst)]
    (if (= 0.0 mag)
      (recur (take (count lst) (repeat 1)))
      (cmato/* (cmato// 1 mag) lst))))



(defn get-rotation-matrix [vectr angle]
  "Returns the Rodrigues' rotation matrix associated with the vector and angle.
           | 00  01  02 |
     rmat: | 10  11  12 |
           | 20  21  22 |"
  (let [[a b c] (normalize-mag vectr)
	      x (float a)
        y (float b)
        z (float c)
        one (int 1);primitive locals = 10X faster!
	sina (cmat/sin angle) cosa (cmat/cos angle) dcosa (cmato/- one cosa)
	r00 (cmato/+ (cmato/* x x) (cmato/* (cmato/- one (cmato/* x x)) cosa))
	r01 (cmato/- (cmato/* x y dcosa) (cmato/* z sina))
	r02 (cmato/+ (cmato/* x z dcosa) (cmato/* y sina))
	r10 (cmato/+ (cmato/* x y dcosa) (cmato/* z sina))
	r11 (cmato/+ (cmato/* y y) (cmato/* (cmato/- one (cmato/* y y)) cosa))
	r12 (cmato/- (cmato/* y z dcosa) (cmato/* x sina))
	r20 (cmato/- (cmato/* x z dcosa) (cmato/* y sina))
	r21 (cmato/+ (cmato/* y z dcosa) (cmato/* x sina))
	r22 (cmato/+ (cmato/* z z) (cmato/* (cmato/- one (cmato/* z z)) cosa))]
    [[r00 r01 r02] [r10 r11 r12] [r20 r21 r22]]))


(defn align-vectors
  "Rotates an arbitrary vector A so that it points in the direction of
  some other arbitrary vector B. This produces a rotation matrix."
  [A B]
  (let [angle1 (find-angle A B)
	      vector1 (cmat/cross A B)]
	       (get-rotation-matrix vector1 angle1)))



(defn common-factorize
  [v-reals]
  (let [f #(if (= (cmat/floor %) %) 1 ((comp denominator rationalize) %))
        fa (map f v-reals)
        m (apply max (map f v-reals))]
        (map (comp int (partial cmato/* m)) v-reals)))













