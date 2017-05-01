# The Context Tree Weighting Method - A Python Implementation

The context tree weighting method is universal data compression algorithm. This project is a Python implementation of the algorithm for the purpose of analyzing DNA sequences. This project was developed for the Bachelor End Project of Luke Martin a electrical engineering student at TU/e/


## Getting Started

You will need to install these packages to run the project on your own Python development environment.

### Prerequisites

* **BioPython** - [link](http://biopython.org/)


## Documentation

An explanation of each function and object within the project

### Node(self, sym, prnt, tree, lvl):

This is the primary object used in the program. These node objects are what populate the context tree.

* `sym` is the symbol that is traversed from the parent to reach this node. "root" if the node is at the root of the tree

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
(with `tab` as a space character)
```
C Prob: -30.108163475962257
 A Prob: -9.281181752544054
  A Prob: -5.415037499278844
  G Prob: -2.415037499278844
  C Prob: -1.0
```
### Tree(self, name, depth)

This object is a context-tree for a sequence of information. It is made up of `Nodes` that point to each other in a typical tree structure.

* `self` - the name given to this tree object
* `depth` - the maximum depth of the tree
The tree also creates its root node when created as follows:
```
      self.root = Node("root", None, self, 0)
```
### add_data_point(self, symbol, context)
This function will add a symbol and context pair into the tree. This is done by calling the `traverse_context` function and then passing the dirty nodes into the `reweight_tree` function

### read_file(self, file)
Will load the `seq` object from the given file and generate the symbol/context pairs and add them to tree using `add_data_point` function

### reweight_tree(self, dirtynodes, symbol)
* `dirtynodes` a list of nodes that need to be updated for the tree to accurately reflect its latest state.
* `symbol` the symbol that was added that caused these nodes to be dirty, needed to properly update the probability of each node

Will reweight all of the probabilities of the dirty nodes strored in `dirtynodes` with regards to the given `symbol`. This is done by calling the `update_probility_log` function on each dirty node

### traverse_context(self, symbol, context)
* `symbol` the symbol that has the given context
* `context` a list of the preceding symbols

Given a `symbol` and `context` pair this function will traverse the tree, iterating the appropriate counter of each Node it reaches. The function will then return a list of the affected nodes. This list represent the "dirty" nodes, meaning the probability stored at that node doesn't represent the current counter values.

### print(self, tab)
Will call the `Node.print` function on its root and therefore printing the whole tree.

### commands(self, command)
Used when interfacing from the console. Takes in the `command` and calls the appropriate functions.


## Console Commands:

### Tree-Manager Commands
```
Commands:
a - Add a new tree
o - Open an existing tree
s - List existing trees
d - Delete a tree
c - show commands

```
### a - Add a new tree:
  Will prompt you too name and create a new tree with a inputted depth

### o - Open an existing tree
  Open a tree, allowing for specific control of that tree

### Tree Commands
```
Commands:
a - Add a node
r - Add random node
p - Print tree
l - Load file into tree
c - Show commands
x - Exit current tree
```
