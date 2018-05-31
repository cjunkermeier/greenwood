(defproject greenwood "0.1.0"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [ [org.clojure/clojure "1.9.0"]
                  [org.clojure/math.combinatorics "0.1.1"]
                  [net.mikera/core.matrix "0.62.0"]
                  [net.mikera/vectorz-clj "0.47.0"]
                  [ubergraph "0.5.0"]
                  [incanter "1.5.5"]
                  [iota "1.1.3"]
                  [criterium "0.4.3"]
                  [me.raynes/fs "1.4.6"]
                  [ultra-csv "0.2.1"]
                  [random-seed "1.0.0"]
                  [proto-repl "0.3.1"]]
  :main  greenwood.core
  :target-path "target/%s"
  :jvm-opts ^:replace ["-server" "-Xmx4G"])
