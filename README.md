# stat-computing

A project to describe simulations in economics and information theory.  Based on course notes by Max Auffhammer and Sofia Berto Villas-Boas, as well as _An Information Theoretic Approach to Econometrics_ by George Judge and Ron Mittelhammer.

# Let's get started #

To get started, you'll need to install a few tools, but it's painless.

* stat-computing (this project)
* Leiningen (Build tool for clojure, located [on github](https://github.com/technomancy/leiningen))
* Plugins

## forma-clj

Fire up your command line and:


```bash
    git clone https://github.com/danhammer/stat-computing.git
    cd stat-computing
```

## Leiningen

Next install Leiningen, the build tool for Clojure. These instructions are copied from the Leiningen README:

* [Download this script](https://raw.github.com/technomancy/leiningen/stable/bin/lein) which is named `lein`
* Place it on your path so that you can execute it. (I like to use `~/bin`)
* Set it to be executable. (`chmod 755 ~/bin/lein`)

## Plugins

Finally, install the plugins using the `lein` command. This part's easy!

```bash
lein plugin install swank-clojure "1.4.0-SNAPSHOT"
lein plugin install lein-marginalia "0.6.1"
lein plugin install lein-midje "1.0.7"
```

And then, just run `lein deps` to download the dependencies, and run `lein deps` a second time to install them.

And you are DONE. As a sanity check, try compiling via `lein compile`.

## License

Distributed under the Eclipse Public License, the same as Clojure.
