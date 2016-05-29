(defproject greenwood "0.1.0"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [ [org.clojure/clojure "1.8.0"]
                 [org.clojure/math.combinatorics "0.1.1"]
                 [net.mikera/core.matrix "0.44.0"]
                 [net.mikera/vectorz-clj "0.29.0"]
    [incanter/incanter "1.9.0"]
                 [iota "1.1.2"]
                  [criterium "0.4.3"]
                  [me.raynes/fs "1.4.6"]
                  [ultra-csv "0.2.0"]
                  [random-seed "1.0.0"]]
  :main  greenwood.core
  :target-path "target/%s"
  :jvm-opts ^:replace ["-server" "-Xmx4G"])

