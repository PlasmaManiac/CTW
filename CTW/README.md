# The Context Tree Weighting Method - A Python Implementation

The context tree weighting method is universal data compression algorithm. This project is a Python implementation of the algorithm for the purpose of analyzing DNA sequences. This project was developed for the Bachelor End Project of Luke Martin a electrical engineering student at TU/e/


## Getting Started

You will need to install these packages to run the project on your own Python development environment.

### Prerequisites

* **BioPython** - link[http://biopython.org/]


## Documentation

An explanation of each function and object within the project

### Node(self, sym, prnt, tree, lvl):

* 'sym' is the symbol that is traversed from the parent to reach this node.

* `prnt` is parent node to this  node (None if it is the root)

* `tree` the tree in which this node exists

* `lvl` the level of the tree in which this node exists

This object also has some number of stored information

* `self.children` a dictionary containing the children nodes of this node. None if there is not one.
* `self.counts` a dictionary that keeps track of the amount of times this node has seen a specific symbol
* `self.probability` stores the KT - Estimator probability of that node.
* `self.weighted_probability` stores the weighted probability of the  node.

### update_probility(self, symbol)
