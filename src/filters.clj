(ns filters)

(defn indexed
  "Returns a lazy sequence of [index, item] pairs, where items come
  from 's' and indexes count up from zero.
  
  by Stuart Sierra, http://stuartsierra.com/ from clojure.contrib.seq
  
  (indexed '(a b c d))  =>  ([0 a] [1 b] [2 c] [3 d])"
  [s]
  (map vector (iterate inc 0) s))

(defn group-filter
  "Given a col A, this creates a col, B, of cols [b1,b2,...,bn], where the bi
  are made up of elements of A that satisfy the filter and are neighboring elements.

  Usage: (group-filter even? [1 3 4 6 7 8 9 10 11]) => ((4 6)(8)(10))"
  [pred col]
  (letfn [(group [picked whatsleft]
            (if (not (seq whatsleft))
              picked
                (if (seq (take-while pred whatsleft))
                  (recur (conj picked (take-while pred whatsleft)) (drop-while pred whatsleft))
                  (recur picked (drop-while (comp not pred) whatsleft)))))]
    (group [] col)))



(defn multi-filter
 "Given a list of predicates and a col of on which to act, this will produce a
col of cols where the first element is all of the members of the input col that
satisfy the first predicate.  The succeeding elements of the output are all of the
members of the input coll that satisfy the corresponding elements in pred-col.
This means that an element in coll can satisfy more than one element of pred-col
and thus be listed more than once in the output.

Usage: (multi-filter [even? odd? zero?] [0 1 2 3 4 5]) => ((0 2 4) (1 3 5) (0))"
  [pred-col coll]
 (pmap #(filter % coll) pred-col))




(defn map-filter
  "(map-filter :k [{:a :b :c :d :k :found} {:a :b :c :d} {:a :b :c :d :k :other}])
 =>  (:found :other)"
  [f coll]
    (filter identity (map f coll)))


(defn index-filter
  "Similar to Clojure’s filter but that returns the indices instead of the
matches themselves. (Programming Clojure, Chapter 2.6, Where is my for loop?
on page 72)"
  [pred coll]
  (when pred
    (for [[idx elt] (indexed coll) :when (pred elt)] idx)))


(defn filter-repeats
  "from npt11tpn
Returns a vector of the elements of a sequence that are repeated n times.
Usage: (filter-repeats '(2 3 1 4 2 2 4 3 5 ) 2) => [3 4]"
  [n coll]
  (let [counts (reduce (fn [m k] (assoc m k (inc (m k 0))))
                       {} coll)]
    (for [[k v] counts :when (= v n)] k)))


#_(defn multi-filter [filters coll]
  (let [c (count filters)
        ignore (Object.)
        merged-seq
          (mapcat
            (fn [k]
              (map
                (fn [p e]
                  (if (p e) e ignore))
                filters (repeat c k)))
            coll)]
    (map
      (fn [i]
        (remove #(= % ignore)
          (take-nth c
            (drop i
              merged-seq))))
      (range c))))


#_(defn multi-filter [filters coll]
  (let [c (count filters)
        ignore (Object.)]
    (map
      (fn [i]
        (remove #(= % ignore)
          (take-nth c
            (drop i
              (mapcat #(map (fn [p e] (if (p e) e ignore)) filters (repeat c
%)) coll)))))
      (range c))))



(comment "There is an interesting discussion of some filtering functions that
would be a good fit for this library.
http://groups.google.com/group/clojure/browse_thread/thread/fda5a5b3ccf9678f")




(defn between-elements
  "Given a col and two predicates this will find all of the elements between the two .
Usage: (between-elements [1 2 3 4 5 6] #(= % 2) #(= % 5)) => (3 4)
Usage: (between-elements [1 2 3 4 5 6] #(< % 3) #(> % 4)) => (3 4)
Usage: (between-elements [1 2 3 4 5 6] #(< % 3) #(> % 2)) => (3 4 5 6)"
[col f1 f2]
((comp first #(partition-by f2 %) last #(partition-by f1 %)) col))


(comment "http://groups.google.com/group/clojure/browse_thread/thread/99b3d792b1d34b56")

