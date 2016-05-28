(ns fortran-format
 ^{:doc "Many legacy (and new) scientific and engineering programs are written in fortran.  In
 many of these programs datafiles that need to be read in at runtime must follow a rigid
 formatting structure designed by the developer.  Fortran has a precise set of rules that may
 be used in formatted input and output.  This code is not exhaustive in its ability to follow
 all of the rules and intricies of the more exotic formatting capabilities of Fortran, but this
 should handle most use cases."
   :author "Chad E. Junkermeier, PHD."}
 (:use [clojure.pprint :only [cl-format]]
      [clojure.string :only [join split replace-first]]))



(defn nil-repeat
 "Fortran can include repeating of formatting.  This is indicated by an integer before the form indicator."
 [i string]
 (if (empty? i) string (join "" (repeat (read-string i) string))))



(defn int-frmt
 "This produces the correct format cl-format description of a
 fortran integer format description.  Where the integer can be
 any of the fortran integer types: i, b, o, or z which
 correspond to the cl-format types: D, B, O, and X.

 The general format for this is Iw.m where w is the number of
 positions to be used and m is the minimum number of positions
 to be used.  Currently this doesn't implement m correctly, it assumes m=w."
 [s I-re]
 (let [[a b] (split s I-re)
       [c d] (split b #"\.")
       II (cond
          (= (str I-re) "[Ii]") "D"
          (= (str I-re) "[Zz]") "X"
          :else (str I-re))]
 (cond
   (nil? d) (nil-repeat a (str "~" c  II))
   (integer? (read-string d))(nil-repeat a (str "~" c ",'0" II)))))




(defn dec-frmt
 "This produces the correct cl-format description from a fortran
 decimal format description."
 [s]
 (let [[a b] (split s #"[fF]")
       [c d] (split b #"\.")]
   (cond
    (nil? d) (nil-repeat a (str "~" c "F" ))
    (integer? (read-string d))(nil-repeat a (str "~" c "," d ",0,'*F")))))


(defn space-frmt
 "This produces the correct cl-format description from a fortran
 space format description."
 [s]
 (if-let [x  (first (split s #"[Xx]"))]
   (str "~" (read-string x) "@T")
   (str "~@T")))


(defn char-frmt
 "This produces the correct cl-format description from a fortran
 character format description."
 [s]
 (let [[a b] (split s #"[aA]")]
   (nil-repeat a (str "~" b "A"))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn logic-frmt
 "This produces the correct cl-format description from a fortran
 Logic format description."
   [s]
(let [[a b] (split s #"[lL]")]
   (nil-repeat a (str "~" b "A"))))

(defn expo-frmt
  "Transforms fortran exponential format to cl-format exponential."
  [s]
  (let [[a b] (split s #"[eE]")
       [c d] (split b #"\.")]
   (cond
    (nil? d) (nil-repeat a (str "~" c "E" ))
    (integer? (read-string d))(nil-repeat a (str "~" c "," d ",0,'*E")))))


(defn tab-frmt
      [s]
 (let [[a b] (split s #"[tT]")]
   (nil-repeat a (str "~" b "A"))))

(defn tab-frmt
      [s]
 (let [[a b] (split s #"[tT]")]
   (nil-repeat a (str "~" b "A"))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defn- replace-parens-
      "Helper function that actually does the work or replacing a set of parens.
      The function replace-parens counts the number of sets of parens and then
      calls this function in a loop until all parens are replaced."
 [s]
 (let [re #"[0-9]*\(.*?\)"
       m (re-find  re s)
     ss (split m #"[\(\)]")
      sss (if (empty? (first ss))
                  "1"
                  (first ss))]
 (replace-first s re (join "," (repeat (read-string sss) (second ss))))))


(defn replace-parens
 "Replace each set of parens with the number preseding the left paren.
 Usage: (replace-parens 's3(tri)n3(g)') => 'stritritringgg'"
[s]
(let [n (count (re-seq #"[0-9]*\(.*?\)" s))]
  (loop [string s
         i n]
    (if (zero? i)
      string
  (recur (replace-parens- string) (dec i))))))



(defn convert-frmt
 "This is the main function to use.  Given a string containing the fortran
 formatting statement it will produce the equivalent statement to be used in cl-format.

 Usage: (convert-frmt \"15x,2i4,f8.4,f8.2,f8.5,f10.7\") => \"~15@T~4D~4D~8,4,0,'*F~8,2,0,'*F~8,5,0,'*F~10,7,0,'*F\""
 [s]
 (let [f #(cond
           (re-find #"[Ff]" %) (dec-frmt %)
           (re-find #"[xX]" %) (space-frmt %)
           (re-find #"[Aa]" %) (char-frmt %)
           (re-find #"[eE]" %) (expo-frmt %)
           (re-find #"[iI]" %) (int-frmt % #"[Ii]")
           (re-find #"[bB]" %) (int-frmt % #"[bB]")
           (re-find #"[oO]" %) (int-frmt % #"[oO]")
           (re-find #"[zZ]" %) (int-frmt % #"[zZ]")
           (re-find #"[Ll]" %) (logic-frmt %))]
   (join "" (map f (split (replace-parens s) #",")))))




