(ns greenwood.utils
  (:require [clojure.string :as strng]
            [clojure.set :as cset]
            [clojure.java.io :as jio]))

;(set *warn-on-reflection* true)

(def endline (System/getProperty "line.separator"))

(def slash (System/getProperty "file.separator"))

(def home (System/getProperty "user.home"))

(def tmp (System/getProperty "java.io.tmpdir"))

(defn nanotime [] (System/nanoTime))

(def real-num #"-*\d+\.\d*[DE]?-?\d*")

(def integer-num #"-*\d+")

;(defn empty-str [#^String s]
;  (= "" (trim s)))


(defn dirname
  "Return directory name of path.\n\t(dirname \"a/b/c\") -> \"/a/b\""
  [path]
  (.getParent (jio/as-file path)))



(defn append-file
  ""
  [file string]
  (do (with-open [w (clojure.java.io/writer  file :append true)]
    (.write w string))
    1))


(defn flatten-n [n coll]
  "Like flatten, but only goes n levels deep."
  (if (= n 0)
    coll
    (recur (dec n) (apply concat (map #(if (sequential? %) % (list %)) coll)))))



(defn filter-user-code [trace]
  (let [java #{"LazySeq.java" "Cons.java" "Compiler.java" "RT.java" "Var.java" "AFn.java" "RestFn.java" "Thread.java"}
	clj #{"core.clj" "basic.clj" "NO_SOURCE_FILE"}
	code (clojure.set/union java clj)]
    (remove #(code (.getFileName %)) trace)))



(defn get-name-number [e]
  (let [s (filter-user-code (.getStackTrace e))]
    (map #(list (.getFileName %) (.getLineNumber %)) s)))



(defn find-cause [exception]
  (let [c (.getCause exception)]
    (if c (distinct (concat (get-name-number exception) (find-cause c)))
      (get-name-number exception))))



(defn serialize
  "This will produce a string from the data in a way that preserves
its type when read in later."
  [data]
  (binding [*print-dup* true] (println-str data)))



(defn transpose [x]
  (if x (apply map list x) (throw (RuntimeException. "Can't transpose a nil matrix."))))



(defn grep
  "Filters elements of coll by a regular expression.  The String
  representation (with str) of each element is tested with re-find."
  [re coll]
  (filter (fn [x] (re-find re (str x))) coll))



(defmacro dbg
  "Prints out intermediate values in functions.

Here's an example function that we might want to debug:
  (defn pythag [ x y ] (* (* x x) (* y y)))

And here is a version enhanced to print out its thought process as it runs
  (defn pythag [ x y ]  (dbg (* (dbg (* x x)) (dbg (* y y)))))

  (pythag 4 5)
;; (* x x) = 16
;; (* y y) = 25
;; (* (dbg (* x x)) (dbg (* y y))) = 400
;; 400"
  [x] `(let [x# ~x] (println '~x "=" x#) x#))




(defn clean-parse [file]
"reading in, chomping, triming a data file."
  (with-open [rdr (jio/reader file)]
     (remove empty? (into [] (line-seq rdr)))))




(defn lazy-chunk-file [filename regex]
  "This outputs a coll of strings, each string containing all of the characters
in the file in between two instances of the regex.

Usage: (lazy-chunk-file filename regex)"
  (let [scanner (java.util.Scanner. (jio/file filename))
        get-next (fn get-next []
                   (if (.hasNext scanner) (lazy-seq (cons (.next scanner) (get-next)))))]
    (.useDelimiter scanner regex)
    (lazy-seq (get-next))))



(defn clean-parse-empty [file]
"reading in, chomping, triming a data file.  This creates a col of cols where
the sub cols are created at each empty line.  This is nice for using with data
files where blank lines are often significant to the data structure."
  (->> (lazy-chunk-file file #"\n(\s*\n)+")
	(map #(strng/split % #"\n") )
       (map (fn [x](map strng/trim x)) )))



(defn ls-to-ll [listlist]
  "Takes a list of strings and turns it into a list of lists.
will be used in turning strings that contain values into lists
that have those values in them. Usage:
list of strings (\"48 -6 -6 1\" \"48 -7 -4 -1\")
becomes
list of lists ((48 -6 -6 1) (48 -7 -4 -1))"
  (map #(map read-string (seq (strng/split % #"\s+"))) listlist))



(defn parse-table
  "Parses a table where the first line  of the file is made up of the headings for each column.  This creates a hash-map
  where each column of the table is an element of the hash-map."
  [filename]
  (let [data (clean-parse filename)
        headings (map keyword (strng/split (strng/replace (first data) #":" "") #"\s+"))]
    (apply hash-map (interleave headings ((comp transpose ls-to-ll rest) data)))))

(defn parse-table-new
  "Parses a table where the first line  of the file is made up of the headings for each column.  This creates a hash-map
  where each column of the table is an element of the hash-map."
  [filename]
  (let [data (clean-parse filename)
        headings (map keyword (strng/split (strng/replace (first data) #":" "") #"\s+|\t+|\^I"))]
    (apply hash-map (interleave headings ((comp transpose ls-to-ll rest) data)))))



(defn parse-table-hash
  "Parses a table where the first line of the file is made up of the headings for each column.  This creates
  a seq of hash-maps."
  [filename]
  (let [data (clean-parse filename)
        headings (map keyword (strng/split (strng/replace (first data) #":" "") #"\s+"))]
     (mapv (comp (partial apply hash-map) (partial interleave headings)) ((comp ls-to-ll rest) data))))


(defn parse-table-hash-new
  "Parses a table where the first line of the file is made up of the headings for each column.  This creates
  a seq of hash-maps."
  [filename]
  (let [data (clean-parse filename)
        headings (map keyword (strng/split (strng/replace (first data) #":" "") #"\s+|\t+|\^I"))]
     (mapv (comp (partial apply hash-map) (partial interleave headings)) ((comp ls-to-ll rest) data))))



(defn indexed
  "Returns a lazy sequence of [index, item] pairs, where items come
  from 's' and indexes count up from zero.

  (indexed '(a b c d))  =>  ([0 a] [1 b] [2 c] [3 d])"
  [s]
  (map vector (iterate inc 0) s))



(defn positions
  "Returns a lazy sequence containing the positions at which pred
   is true for items in coll."
  [pred coll]
  (for [[idx elt] (indexed coll) :when (pred elt)] idx))



(defn inter-cat-tree [items tree]
  "applies join recursively to tree"
  (if (nil? (next items))
    (strng/join (first items) tree)
    (strng/join (first items)
      (map #(inter-cat-tree (next items) %)
        tree))))



(defn parse
  "Use this to lazily parse a datafile.
  f is some function that is used upon each line.
  n is the number of lines that is dropped from the beginning of the file."
  [filename f n]
  (with-open [rdr (jio/reader filename)]
     (doall (map (comp f strng/trim) (drop n (line-seq  rdr))))))


(defn parse-grep
  "Use this to lazily parse a datafile."
  [filename re]
  (with-open [rdr (jio/reader filename)]
     (filter (partial re-find re ) (line-seq  rdr))))


(defn
 ^{:doc "Same as (first (nnext x))"
   :arglists '([x])}
 third [x]
  (first (nnext x)))



(defn txt->LaTeX-table
  [filename]
 (strng/join " \\\\ \n"  (map #(clojure.string/replace % "\t" " & ") (clean-parse filename))))



(defn hm->LaTeX-table
  "converts a vector of hash maps into a LaTeX style table."
  [hm]
  (inter-cat-tree [" \\\\ \n" " & "]
          (conj  (map vals hm)  (keys (first hm)))))



(defn ignore-sequential
  "This is used with clojure.walk/prewalk or with clojure.walk/postwalk.
Usage:  (prewalk (ignore-sequential inc) [1 [2] 3]) => [2 [3] 4]"
  [f]
    #(if (sequential? %) % (f %)))

(defn drop-n [n col]
  "Works like nth, but instead of keeping only the nth term, this keeps
everything but the nth term."
  (concat (take n col) (drop (inc n) col)))

(defn disjoint?
  "This should be part of the clojure.set namespace, but isn't."
  [x y]
  (empty? (cset/intersection x y)))

(defn create-keyval-vectormap
  "Given a set of keys and a seq of maps, ie. ({} {2 1} {1 1, 2 1}),
  this will create a new seq of maps where each of the new maps has
  all of the keys.  It will then assoc in the values from the
  original seq and set the values of the remaining keys to zero.

  There does not need to be a key-val pair for nil (zero),
  because nil is implied by an absence of the other values.

  Usage: "
  [all-k in-vector-map]
  (let [hm (zipmap (filter (complement nil?) all-k) (repeat 0))]
    (map #(merge hm %1) in-vector-map )))


(defn vec-maps->llist
  "Turns a col of maps (i.e. a mol) into a linked list where
  the first list is the set of keys.
  Useage: (vec-maps->llist ({2 0, 1 0} {2 1, 1 0} {2 1, 1 1})) => ((2 1) (0 0) (1 0) (1 1))"
  [x]
  (let [k (keys (first x))]
    (concat [k] (map vals x))))


#_(ns pretzel.combine
  {:doc "Combine predicates
   Copyright (C) 2011
   Joost Diepenmaat, Zeekat Softwareontwikkeling joost@zeekat.nl - http://joost.zeekat.nl/"})

(defn combine-p
  "Combine predicates using collection function f.
Returns a new predicate that tests the collection of predicates
against the given value(s) using f."
  [f predicates]
  (fn [& values]
    (f #(apply % values) predicates)))

(defn every-p?
  "Create a new predicate which is true if every predicate is true"
  [& predicates]
  (combine-p every? predicates))

(defn some-p?
  "Create a new predicate returning the first true value of predicates"
  [& predicates]
  (combine-p some predicates))

(defn not-any-p?
  "Create a new predicate returning false if all predicates return true"
  [& predicates]
  (combine-p not-any? predicates))

(defn not-every-p?
  "Create a new predicate returning true if any predicate returns false"
  [& predicates]
  (combine-p not-every? predicates))

