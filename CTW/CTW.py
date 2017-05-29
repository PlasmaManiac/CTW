import math
import random
import sys
import gzip
import cProfile
import pstats

from tkinter.filedialog import askopenfilename
from Bio import SeqIO
from collections import defaultdict

def print_tree_commands():
    print("Commands:")
    print("a - Add a node")
    print("r - Add random node")
    print("p - Print tree")
    print("l - Load file into tree")
    print("c - Show commands")
    print("x - Exit current tree")


def print_commands():
    print("Commands:")
    print("a - Add a new tree")
    print("o - Open an existing tree")
    print("s - List existing trees")
    print("d - Delete a tree")
    print("c - show commands")


def add_logs(log1, log2):
    """
    :param log1: Log1 is larger than log2
    :param log2: Second log value
    :return: will return the addition of the two logs whilst only in the log domain
    """
    ans = log1 + math.log(1 + pow(2, (log2 - log1)), 2)
    return ans


def int_to_ACGT(int_list):
    """Helper Function - Allows for the randomly generated sequence to be converted into chars"""
    ACGT_list = []
    for x in int_list:
        if x == 1:
            ACGT_list.append("A")
        elif x == 2:
            ACGT_list.append("C")
        elif x == 3:
            ACGT_list.append("G")
        elif x == 4:
            ACGT_list.append("T")

    return ACGT_list


def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 3, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    sys.stdout.write('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix))
    # Print New Line on Complete
    if iteration == total:
        print()

def num_to_ACGT(x):
    if x == 1:
        return "A"
    elif x == 2:
        return "C"
    elif x == 3:
        return "G"
    elif x == 4:
        return "T"


def print_trees(trees):
    print("\n Listing all current trees: ")
    print("\t Name:", "\t Depth", "\t\t| Root Prob.|\t\t Counts")
    print("\t ____________________________________________________")
    for key in trees:
        print("\t", key, "\t \t", trees[key].depth,
              "\t\t", trees[key].root.weighted_probability,
              "\t\t",  trees[key].root.counts)


def inverted_complement(symbols):
    reverse = [None] * len(symbols)

    for i, sym in enumerate(symbols):
        if sym == "A":
            reverse[i] = "T"
        elif sym == "T":
            reverse[i] = "A"
        elif sym == "G":
            reverse[i] = "C"
        elif sym == "C":
            reverse[i] = "G"

    reverse.reverse()
    return reverse


class Node:
    def __init__(self, sym, prnt, tree, lvl):
        # symbol sequence of that node
        self.symbol = sym

        # the tree in which this Node exists
        self.tree = tree

        # Parent Node
        self.parent = prnt
        # The edges of a node, corresponding to a new node.
        # eg: A => The A branch from current node
        self.children = {"A": None, "C": None, "G": None, "T": None}

        # Counts
        self.counts = {"A": 0, "C": 0, "G": 0, "T": 0}

        # Estimated probability given the KT Estimator
        self.estimated_probability = 0

        # Estimated probability given the KT Estimator
        # With an alpha of 0.05
        self.estimated_probability_alpha = 0

        # Used to store the weighted probability of the node
        self.weighted_probability = 0

        # The level of the node in the tree
        # level = 0 is the root
        self.level = lvl

        self.weighted_probabilities = dict()

        for x in range(self.level, self.tree.depth + 1):
            self.weighted_probabilities[x] = 0

    def get_weighted_prob(self, depth, maxdepth):

        if depth is 0:
            if maxdepth > 11:
                return self.estimated_probability_alpha
            else:
                return self.estimated_probability
        else:
            children_prob = 0
            for key, child in self.children.items():
                if child is not None:
                    children_prob += child.get_weighted_prob(depth - 1, maxdepth)

        if maxdepth < 12:
            return children_prob + math.log(1 + pow(2, (self.estimated_probability - children_prob)), 2) + math.log(1/2, 2)
        else:
            return children_prob + math.log(1 + pow(2, (self.estimated_probability_alpha - children_prob)), 2) + math.log(1/2, 2)


    def update_probability_log(self, symbol):
        """When called will calculate the probability of the node"""

        total_count = 0

        for key, count in self.counts.items():
            total_count += count

        self.tree.fixed_depth_probabilities[self.level] -= self.estimated_probability
        self.tree.fixed_depth_probabilities_alpha[self.level] -= self.estimated_probability_alpha
        # Base 2 log
        # CHANGED: self.tree.alpha -> 0.5
        self.estimated_probability += math.log((self.counts[symbol] + 0.5)
                                               / (total_count + 4*0.5), 2)

        self.estimated_probability_alpha += math.log((self.counts[symbol] + 0.01)
                                               / (total_count + 4*0.01), 2)

        self.tree.fixed_depth_probabilities[self.level] += self.estimated_probability
        self.tree.fixed_depth_probabilities_alpha[self.level] += self.estimated_probability_alpha

        children_prob = 0

        for key,value in self.weighted_probabilities.items():
            if key is self.level:
                if key > 11:
                    self.weighted_probabilities[key] = self.estimated_probability_alpha
                else:
                    self.weighted_probabilities[key] = self.estimated_probability
            else:
                child_weight = 0
                for child_key, child in self.children.items():
                    if child is not None:
                        child_weight += child.weighted_probabilities[key]
                if key > 11:
                    self.weighted_probabilities[key] = add_logs(child_weight, self.estimated_probability_alpha) + math.log(1 / 2, 2)
                else:
                    self.weighted_probabilities[key] = add_logs(child_weight, self.estimated_probability) + math.log(1 / 2, 2)


        if self.level is not self.tree.depth:
            for key, child in self.children.items():
                if child is not None:
                    children_prob += child.weighted_probability
            # log(1/2,2) is equivalent to multiplying it all 1/2
            self.weighted_probability = (
                children_prob
                + math.log(1 + pow(2, (self.estimated_probability - children_prob)), 2)
                + math.log(1 / 2, 2))
        else:
            self.weighted_probability = self.estimated_probability

    def update_probability(self, symbol):

        total_count = 0
        for key, count in self.counts.items():
            total_count += count

        self.estimated_probability *= (self.counts[symbol] + 0.5) / (total_count + 2)

        children_prob = 1
        if self.level is not self.tree.depth:
            for key, child in self.children.items():
                if child is not None:
                    children_prob *= child.probability
            self.weighted_probability = (self.estimated_probability + children_prob) / 2
        else:
            self.weighted_probability = self.estimated_probability

    def print(self, tab):
        # Todo: Print self, based on the depth of the node and  children
        print(tab * (self.level - 1), self.symbol, "Prob: " + str(self.estimated_probability))

        for key, child in self.children.items():
            if child is not None:
                child.print(tab)
        return


class Tree:
    def __init__(self, name, depth, best=False):
        self.name = name
        self.depth = depth
        self.best = best

        self.root = Node("root", None, self, 0)

        self.depthTable = defaultdict(list)
        self.depthTable[0] = [self.root]
        self.fixed_depth_probabilities = dict()
        self.fixed_depth_probabilities_alpha = dict()
        for x in range (0, depth + 1):
            self.fixed_depth_probabilities[x] = 0
            self.fixed_depth_probabilities_alpha[x] = 0


        if depth > 11:
            self.alpha = 0.05
        else:
            self.alpha = 0.5

    def calculate_markov_prob(self, depth):
        probability_sum = 0

        for node in self.depthTable[int(depth)]:
            node_prob = 0
            if int(depth) > 11:
                node_prob = node.estimated_probability_alpha
            else:
                node_prob = node.estimated_probability

            if probability_sum is 0:
                probability_sum = node_prob
            else:
                # probability_sum = add_logs(probability_sum, node_prob)
                probability_sum += node_prob

        return probability_sum

    def best_tree(self):

        # TODO: Finish this
        best_tree = Tree("Best Tree", self.depth, True)

        # First check if the root itself is not best, IE: Just KT
        if self.root.weighted_probability >= self.root.estimated_probability:
            # Effectively making a new tree,
            best_tree.root.counts = self.root.counts.copy()
            best_tree.root.children = self.root.children.copy()

            for key, child in best_tree.root.children.items():
                if child.weighted_probability >= child.estimated_probability:
                    print("SADKALSJD")
        return best_tree

    def read_weighted_probability(self, depth):
        """
        Will calculate and return the weighted probability of a CTW with given depth
        :param depth: The depth of the tree you want to read
        :return: A probability in log base 2 of the CTW of given depth
        """

        # TODO: Add in a conisderation for how the depth may affect ALPHA IE: Depths > 11 -> alpha = 0.05 DONE

        if depth > self.depth:
            return None
        if depth < 0:
            return None

        children_prob = 0

        for key, child in self.root.children.items():
            if child is not None:
                children_prob += child.get_weighted_prob(depth - 1, depth)

        if depth < 11 and depth is not 0:
            return children_prob + math.log(1 + pow(2, (self.root.estimated_probability - children_prob)), 2) \
                   + math.log(1 / 2, 2)
        elif depth is 0:
            return self.root.estimated_probability
        else:
            return children_prob + math.log(1 + pow(2, (self.root.estimated_probability_alpha - children_prob)), 2) \
                   + math.log(1 / 2, 2)

    def add_data_point(self, symbol, context):
        """Will add the data point to tree."""
        dirty_nodes = self.traverse_context(symbol, context)
        self.reweight_tree(dirty_nodes, symbol)
        self.add_inverted_point(symbol, context)

    def add_inverted_point(self, symbol, context):
        other_context = list(context)
        other_context.insert(0, symbol)
        # print(other_context);

        reverse = inverted_complement(other_context)
        new_symbol = reverse.pop(0)

        # print(reverse, "Symbol: ", new_symbol)

        dirty_nodes = self.traverse_context(new_symbol, reverse)
        self.reweight_tree(dirty_nodes, new_symbol)

    def read_file(self, file):

        file_seq = file.seq

        print(file_seq)

        context = []

        for index, symbol in enumerate(file_seq):
            if len(context) is not self.depth:
                #  Add to the context
                context.insert(0, symbol)
            else:
                #  Add to the tree
                # print("context: ", context)
                # print("symbol: ", symbol)

                self.add_data_point(symbol, context)

                context.pop()
                context.insert(0, symbol)

    def reweight_tree(self, dirtynodes, symbol):
        self.root.counts[symbol] += 1

        for node in dirtynodes:
            node.update_probability_log(symbol)

        self.root.update_probability_log(symbol)

    def traverse_context(self, symbol, context):
        # Start at the root
        current_node = self.root
        full_path = ""
        dirty_nodes = []
        depth = 1
        # For Each symbol in the context, traverse down from the root.
        for sym in context:

            full_path += sym
            if current_node.children[sym] is not None:
                current_node = current_node.children[sym]
            else:
                current_node.children[sym] = Node(sym, current_node, current_node.tree, depth)
                # self.depthTable[depth].append(current_node.children[sym])
                current_node = current_node.children[sym]
            depth += 1
            dirty_nodes.insert(0, current_node)
            current_node.counts[symbol] += 1

        return dirty_nodes

    def print(self, tab):
        print("Root: ", self.root.weighted_probability, " ", self.root.counts)

        for key, child in self.root.children.items():
            if child is not None:
                child.print(tab)

    def commands(self, command):
        # Adds a random sequence to the tree
        if command == "r" or command == "R":
            print("    Random:")
            context = []
            for x in range(0, self.depth):
                context.append(random.randint(1, 4))
            context = int_to_ACGT(context)
            symbol = num_to_ACGT(random.randint(1, 4))
            print(context, "  symbol: ", symbol)
            self.add_data_point(symbol, context)

        # Print out the tree
        elif command == "p" or command == "P":
            print("     Print:")
            self.print(input('     Enter your tab char: '))

        elif command == "l" or command == "L":
            print("Loading a file: ")
            # TODO: Add in a way to manage the records, so to load a specific record

            filename = askopenfilename()
            print(filename)

            file_extension = filename.split(".")[-1]
            iterator = None

            if file_extension == "gbk":
                iterator = SeqIO.parse(filename, "genbank")
            elif file_extension == "gz":
                handle = gzip.open(filename, "rt")
                iterator = SeqIO.parse(handle, "genbank")

            record = next(iterator)
            print(record.seq)

            self.read_file(record)

            iterator.close()

        elif command == "c" or command == "C":
            print_tree_commands()


class Competition:
    def __init__(self, max_depth, blocksize=200):
        self.max_depth = max_depth
        self.trees = dict()
        self.tree = Tree("MD:" + str(max_depth), max_depth)
        self.blocksize = blocksize
        #  Keeps track of the Pw and Pm of when each block was added
        self.block_Pe = defaultdict(list)

    def read_block(self, block):
        context = []
        previous_root_prob = dict()
        previous_root_prob_markov = dict()

        root_prob = dict()
        markov_prob = dict()

        for x in range(1, self.max_depth + 1):
            #  previous_root_prob[x] = self.tree.read_weighted_probability(x)
            previous_root_prob[x] = self.tree.root.weighted_probabilities[x]
            #previous_root_prob_markov[x] = self.tree.calculate_markov_prob(x)
            previous_root_prob_markov[x] = self.tree.fixed_depth_probabilities[x]

            self.block_Pe[x].append((previous_root_prob[x],previous_root_prob_markov[x]))

        for index, symbol in enumerate(block):
            if len(context) is not self.max_depth:
                #  Add to the context
                context.insert(0, symbol)
            else:
                #  Add to the tree
                #  print("context: ", context)
                #  print("symbol: ", symbol)

                # Add to all to trees
                self.add_to_all(symbol, context)

                context.pop()
                context.insert(0, symbol)

        difference = dict()
        difference_markov = dict()
        for x in range(1, self.max_depth + 1):
            # root_prob[x] = self.tree.read_weighted_probability(x)
            root_prob[x] = self.tree.root.weighted_probabilities[x]
            #markov_prob[x] = self.tree.calculate_markov_prob(x)
            markov_prob[x] = self.tree.fixed_depth_probabilities[x]

            difference[x] = root_prob[x] - previous_root_prob[x]
            difference_markov[x] = markov_prob[x] - previous_root_prob_markov[x]

        # print( "\n",previous_root_prob, "\n")
        # print(root_prob, "\n")
        # print(difference, "\n")
        # print("MAX: ", max(difference))
        # print("MIN: ", min(difference))

        return max(difference, key=difference.get), difference[max(difference, key=difference.get)], max(difference_markov, key=difference_markov.get), difference_markov[max(difference_markov, key=difference_markov.get)]

    def add_to_all(self, symbol, context_full):
        #  sys.stdout.write("\rAdding symbol: %s With context: %s " % (symbol, context_full))
        self.tree.add_data_point(symbol, context_full)


    def read_file(self, file_handle):

        records = dict()

        for r in SeqIO.parse(file_handle, "gb"):
            records[r.id] = r.seq

        for key, value in records.items():
            print("ID: ",key)
        print("\n")
        file = input("What record do you want to open?")
        file_seq = records[file]
        #  Just a randomly chosen one
        #file_seq = records["NW_004929430.1"]

        block = list()
        blockseq = list()
        symcounter = 0
        max_symbols = 40000

        for index, symbol in enumerate(file_seq):

            if symbol is not "N":

                if symcounter > max_symbols:
                    break

                block.append(symbol)
                symcounter += 1

                if len(block) is self.blocksize:
                    printProgressBar(symcounter, max_symbols)
                    blockseq.append(self.read_block(block))
                    del block[:]

        print('\n',"Length of block sequence" , len(blockseq))
        print_to_file(blockseq)


class Codon:
    def __init__(self, max_depth, type = "normal"):

        # TODO: Make this change based on compettive or not
        if type is "competitive":
            print("TODO: SOMETHIGN ")

        # Creates the three seperate trees that will be used for
        self.trees = {1 : Tree("Codon1",max_depth),
                      2 : Tree("Codon2",max_depth),
                      3 : Tree("Codon3",max_depth)}

        # Keeps track of the next tree that will be used
        self.current_tree = 1
        self.max_depth = max_depth

        def read_block(block):
            if len(block) % 3 is not 0:
                print("Block is not a multiple of three, this could lead to some errors")

            context = []

            for index, symbol in enumerate(block):
                if len(context) is not self.max_depth:
                    #  Add to the context
                    context.insert(0, symbol)
                else:
                    self.trees[self.current_tree].add_data_point(symbol, context)
                    self.current_tree += 1
                    if self.current_tree is 3:
                        self.current_tree = 1

                    context.pop()
                    context.insert(0, symbol)


def print_to_file(output, file1="CTW.txt", file2="markov.txt"):
    with open(file1, 'w') as f:
        out = ""
        for item in output:
            out += str(item[0]) + ", " + str(item[1]) + "; \n"
        f.write(out)

    if len(output[0]) < 4:
        return

    with open(file2, 'w') as f:
        out = ""
        for item in output:
            out += str(item[2]) + ", " + str(item[3]) + "; \n"
        f.write(out)


def main():

    # TODO: Investigate why the probabilities of the return values are wrong
    # TODO: The ones return by the two methods of calculating the weighted probability are not correct
    # TODO: Try loading it into a tree with fixed depth and check there what is correct
    # TODO: do this for each tree

    # TODO: Figure out how ALPHA is affecting the probabilities, since its having a strong effect on
    # TODO: the overall probability, also send email asking this question

    filename = "tree.profile"
    pr = cProfile.Profile()
    pr.enable()

    # Allows for user input to interact with the functions above
    test_tree = Tree("Test5", 5)
    test_tree2 = Tree("Test16", 16)

    trees = dict()
    test = Competition(16)
    codon_test = Codon(10)

    block = []


    handle = gzip.open("hs_alt_CHM1_1.1_chr22.gbk.gz", "rt")
    test.read_file(handle)


    print_trees(test.trees)
    current_tree = test.tree

    trees[test_tree.name] = test_tree
    trees[test.tree.name] = test.tree
    trees[test_tree2.name] = test_tree2

    # current_tree = current_tree.best_tree()
    current_tree = test.tree

    # Table of commands
    print_commands()

    print_trees(trees)

    pr.disable()
    pr.dump_stats(filename)
    stats = pstats.Stats(filename)
    # stats.strip_dirs().sort_stats(-1).print_stats()
    # stats.sort_stats('time').print_stats(10)

    while 1:

        if current_tree is None:

            command = input('Enter your command: ')
            if command is "a":
                tree_name = input("Name of the new tree: ")
                depth_name = input("Depth:")
                trees[tree_name] = Tree(tree_name, int(depth_name))
                print_trees(trees)

            elif command is "o":
                name = input("\tWhich tree do you want to open?\n\t")
                current_tree = trees[name]
                print_tree_commands()

            elif command is "s":
                print_trees(trees)

            elif command is "d":
                name = input("\tWhich tree do you want to remove?\n\t")
                trees.pop(name, None)
                print_trees(trees)

            elif command is "c":
                print_commands()

        else:
            print("\n Currently in tree: ", current_tree.name)
            command = input("\t Enter your command: ")

            if command is "b":
                depth = input("what max depth do you want to calculate for? \n")
                if depth.isdigit():
                    for x in range(0, int(depth) + 1):
                        print(str(x) + ": ", current_tree.read_weighted_probability(int(x)))
            if command is "t":
                depth = input("what depth do you want to access? \n")
                if depth.isdigit():
                    print(current_tree.depthTable[int(depth)])

            if command is "m":
                depth = input("what depth do you want to calculate for? \n")
                if depth.isdigit():
                    print(current_tree.calculate_markov_prob(depth))

            if command is "e":
                for key, value in test.block_Pe.items():
                    print_to_file(value, "outputs/ " + str(key) + "Pw and Markov.mat")

            if command is "q":
                print("Weighted Probabilities  stored in root: (a = 0.5)")
                for key, value in current_tree.root.weighted_probabilities.items():
                    print("\t", key, ": ", value)

            if command is "z":
                print("Fixed Depth Probabilities stored in Tree: Alpha = 0.5")
                for key, value in current_tree.fixed_depth_probabilities.items():
                    print("\t", key, ": ", value)

                print("Fixed Depth Probabilities stored in Tree: Alpha = 0.05")
                for key, value in current_tree.fixed_depth_probabilities_alpha.items():
                    print("\t", key, ": ", value)

            if command is "x" or command is "X":
                current_tree = None
                print_trees(trees)

            else:
                current_tree.commands(command)

if __name__ == "__main__": main()



