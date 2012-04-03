(defproject stat-computing "1.0.0-SNAPSHOT"
  :description "Applied econometrics and statistics in Clojure"
  :source-path "src/clj"
  :resources-path "resources"
  :repositories {"conjars" "http://conjars.org/repo/"}
  :marginalia {:javascript ["mathjax/MathJax.js"]}
  :javac-options {:debug "true" :fork "true"}
  :dependencies [[org.clojure/clojure "1.3.0"]
                 [org.clojure/math.numeric-tower "0.0.1"]
                 [incanter/incanter-core "1.3.0-SNAPSHOT"]
                 [incanter/incanter-io "1.3.0-SNAPSHOT"]
                 [incanter/incanter-charts "1.3.0-SNAPSHOT"]
                 [clj-time "0.3.4"]
                 [clojure-csv/clojure-csv "2.0.0-alpha1"]
                 [cascalog "1.9.0-wip"]
                 [cascalog-checkpoint "0.1.1"]
                 [backtype/dfs-datastores "1.1.0"]
                 [backtype/dfs-datastores-cascading "1.1.1"]]
  :dev-dependencies [[org.apache.hadoop/hadoop-core "0.20.2-dev"]
                     [midje-cascalog "0.4.0"]])
