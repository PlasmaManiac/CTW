# The Context Tree Weighting Method - A Python Implementation

The context tree weighting method is universal data compression algorithm. This project is a Python implementation of the algorithm for the purpose of analyzing DNA sequences. This project was developed for the Bachelor End Project of Luke Martin a electrical engineering student at TU/e/


## Getting Started

You will need to install these packages to run the project on your own Python development environment.

### Prerequisites

* **BioPython** - link[http://biopython.org/]


## Documentation

An explanation of each function and object within the project

### Node(self, sym, prnt, tree, lvl):

* `sym` is the symbol that is traversed from the parent to reach this node.

* `prnt` is parent node to this  node (None if it is the root)

* `tree` the tree in which this node exists

* `lvl` the level of the tree in which this node exists

This object also has some number of stored information

* `self.children` a dictionary containing the children nodes of this node. None if there is not one.
* `self.counts` a dictionary that keeps track of the amount of times this node has seen a specific symbol
* `self.probability` stores the KT - Estimator probability of that node.
* `self.weighted_probability` stores the weighted probability of the  node.

### update_probility_log(self, symbol)
This function will update both the  `self.probability` and `self.weighted_probability` of the node based on the inputted `symbol`. This function does not output anything.

### print(self, tab):
Will print information about the node and its children. The `tab` input determines what symbol is used as indents between levels of the tree

Example Print out:
(using tab as a space)
```
C Prob: -30.108163475962257
 A Prob: -9.281181752544054
  A Prob: -5.415037499278844
  G Prob: -2.415037499278844
  C Prob: -1.0
```
### Tree(self, name, depth)

* `self` the name given to this tree object
* `depth` the maximum depth of the tree
