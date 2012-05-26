(defproject stat-computing "1.0.0-SNAPSHOT"
  :description "Applied econometrics and statistics in Clojure"
  :resources-path "resources"
  :repositories {"conjars" "http://conjars.org/repo/"}
  :marginalia {:javascript ["mathjax/MathJax.js"]}
  :javac-options {:debug "true" :fork "true"}
  :plugins [[lein-swank "1.4.4"]]
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
                 [backtype/dfs-datastores-cascading "1.1.1"]
                 [lein-swank "1.4.4"]
                 [clj-http "0.4.1"]
                 [cheshire "4.0.0"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [org.clojure/data.json "0.1.2"]
                 [org.apache.mahout/mahout-core "0.5"]
                 [org.apache.mahout/mahout-math "0.5"]
                 [org.apache.mahout/mahout-utils "0.5"]]
  :dev-dependencies [[org.apache.hadoop/hadoop-core "0.20.2-dev"]
                     [midje-cascalog "0.4.0"]])

