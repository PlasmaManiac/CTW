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
    """
    Helper Function - Prints the options available when within a tree
    """
    print("Commands:")
    print("a - Add a node")
    print("r - Add random node")
    print("p - Print tree")
    print("l - Load file into tree")
    print("c - Show commands")
    print("x - Exit current tree")

def print_commands():
    """
    Helper function - Prints the options when not-in a tree
    :return:
    """
    print("Commands:")
    print("a - Add a new tree")
    print("o - Open an existing tree")
    print("s - List existing trees")
    print("d - Delete a tree")
    print("c - show commands")


def add_logs(log1, log2):
    """
    Help function - Adds two log values under the assumption that log1 >= log2
    :param log1: Log1 is larger than log2
    :param log2: Second log value
    :return: will return the addition of the two logs
    """
    return log1 + math.log(1 + pow(2, (log2 - log1)), 2)



def int_to_ACGT(int_list):
    """
    Helper Function - Allows for the randomly generated sequence to be converted into chars
    converts as follows:
        1 - A
        2 - C
        3 - G
        4 - T

    :param int_list: A list of ints of either 1,2,3 or 4
    :return will return a list of a A,C,G or T characters converted from int_list
    """
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

def num_to_ACGT(x):
    """
    Convert a single value to ACGT
    :param x: an int of 1,2,3 or 4
    :return: A char
    """
    if x == 1:
        return "A"
    elif x == 2:
        return "C"
    elif x == 3:
        return "G"
    elif x == 4:
        return "T"

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



def print_trees(trees):
    """
    Provides a print out of relevent information
    :param trees: A list of trees
    """
    print("\n Listing all current trees: ")
    print("\t Name:", "\t Depth", "\t\t| Root Prob.|\t\t Counts")
    print("\t ____________________________________________________")
    for key, value1 in trees.items():
        if isinstance(value1, Codon):
            print("\t",key, ":")
            for tree,value2 in value1.trees.items():
                print("\t\t",key, tree, "\t \t", value2.depth,
                      "\t\t", value2.root.weighted_probability,
                      "\t\t", value2.root.counts)
        else:
            print("\t", key, "\t \t", trees[key].depth,
            "\t\t", trees[key].root.weighted_probability,
            "\t\t",  trees[key].root.counts)


def inverted_complement(symbols):
    """
    This function will take a set of symbols, and return the inverted reverse of them

    :param symbols: A set of symbols that will be converted into their
    :return: The reversed and inverted sequence
    """
    #  Initialize an empty list
    reverse = [None] * len(symbols)

    #  for all the elements in symbols, convert and insert into reverse
    for i, sym in enumerate(symbols):
        if sym == "A":
            reverse[i] = "T"
        elif sym == "T":
            reverse[i] = "A"
        elif sym == "G":
            reverse[i] = "C"
        elif sym == "C":
            reverse[i] = "G"

    #  reverse the list
    reverse.reverse()
    #  return the list
    return reverse


class Node:
    def __init__(self, sym, prnt, tree, lvl):
        # symbol used to access that node from parent
        # So the A-child of a node, A is stored here
        # Solely used for printing the tree
        self.symbol = sym

        # the tree in which this Node exists
        self.tree = tree

        # Parent Node
        self.parent = prnt
        # The edges of a node, corresponding to a new node.
        # eg: A => The A branch from current node
        self.children = {"A": None, "C": None, "G": None, "T": None}

        # Counts of each symbol seen by that node
        self.counts = {"A": 0, "C": 0, "G": 0, "T": 0}

        # Estimated probability given the KT Estimator
        # with an alpha of 0.5
        self.estimated_probability = 0

        # Estimated probability given the KT Estimator
        # With an alpha of 0.05
        self.estimated_probability_alpha = 0

        # The weighted probability of the node
        self.weighted_probability = 0

        # The level of the node in the tree
        # level = 0 is the root
        self.level = lvl

        # A dictionary containing all the possible weighted probabilities that could use this node
        # Eg: a Node at level 2 of a tree with depth 4, can be part of a weighted probabilty of depth 2,3 or 4.
        # So it has entries for 2 3 and 4
        self.weighted_probabilities = dict()

        # Adds the correct keys to dictionary for later access
        # Initialize the weighted_probabilities to 0

        for x in range(self.level, self.tree.depth + 1):
            self.weighted_probabilities[x] = 0


    def update_probability_log(self, symbol):
        """
        Given a symbol, update the probability stored in the node
        :param symbol: A char of the symbol that has just been added to the tree
        :return:
        """

        # Initialize a value to store the counts of all symbols stored
        total_count = 0

        # For each count value, add it to the total_count
        for key, count in self.counts.items():
            total_count += count

        # By removing the estimated probability from the current summation
        # It is possible to then change it, and re-add it to the probability
        self.tree.fixed_depth_probabilities[self.level] -= self.estimated_probability
        self.tree.fixed_depth_probabilities_alpha[self.level] -= self.estimated_probability_alpha

        # Calculation of the estimated probability with alpha = 0.5
        # Ths is the KT-Estimater as described for standard CTW usage.
        # The '-1' is to account for the fact that the symbol counts have already been incremented
        self.estimated_probability += math.log((self.counts[symbol] - 1 + 0.5)
                                               / (total_count - 1 + 4*0.5), 2)

        # Calculation of the estimated probability (PINHO)
        # This is the same KT-estimater, expect with an alpha of 0.05 instead of 0.5
        self.estimated_probability_alpha += math.log((self.counts[symbol] - 1 + 0.05)
                                               / (total_count - 1 + 4*0.05), 2)

        # Re-adding to the fixed depth proabilities to account for the change
        # These are the probabilities as described in the Pinho paper
        self.tree.fixed_depth_probabilities[self.level] += self.estimated_probability
        self.tree.fixed_depth_probabilities_alpha[self.level] += self.estimated_probability_alpha

        # stores the probability of the children to this node
        children_prob = 0

        # Begin by going through all of the weighted probabilities
        # These are the probabilities as described for CTW
        for key,value in self.weighted_probabilities.items():
            # If the level of the weighted probability of is the same as the level of the node:
            # only return the estimated value.
            # For example when calculating for a max depth of 8, the Pw of nodes at depth 8 is
            # only the estimated probability.
            if key is self.level:
                # If the level is greater than 11
                if key > 11:
                    # If the tree is supposed to use the 0.05 alpha value
                    # then store the 0.05 value
                    if self.tree.alpha:
                        self.weighted_probabilities[key] = self.estimated_probability_alpha
                    else:
                        self.weighted_probabilities[key] = self.estimated_probability
                else:
                    # below 12, nominal estimated probability
                    self.weighted_probabilities[key] = self.estimated_probability
            else:
                # So now we know that the children to this node should be taken into consideration
                # when calculating the weighted probability of this node

                # Stores the weighted probability of the children such that it can be used
                # to update the weighted probability at this node
                child_weight = 0

                # Begins by looping through all the children of this node
                for child_key, child in self.children.items():
                    # Check to see if the child exists
                    if child is not None:
                        # If it exist, add to child_weight
                        child_weight += child.weighted_probabilities[key]
                if key > 11:
                    # check to see if it is needed to use the 0.05 alpha value
                    if self.tree.alpha:
                        self.weighted_probabilities[key] = add_logs(child_weight, self.estimated_probability_alpha) + math.log(1 / 2, 2)
                    else:
                        self.weighted_probabilities[key] = add_logs(child_weight, self.estimated_probability) + math.log(1 / 2, 2)
                else:
                    # NOTE: the math.log(1/2,2) is equivalent to dividing the resulting sum by 2
                    #       since pw = (pe + pw_children)/2 the division by 2 is necessary
                    self.weighted_probabilities[key] = add_logs(child_weight, self.estimated_probability) + math.log(1 / 2, 2)

        # This portion of code is used to update the singular
        # weighted probability, this was used before an update was made to consider
        # weighted probabilities of different levels
        # These values are not used but is a hold over from a previous
        # version of the implementation

        # First check to see if the depth is equal to level
        # Will affect the value used in the weighted prob.
        if self.level is not self.tree.depth:
            # Loop through all the children
            for key, child in self.children.items():
                # if the child exists..
                if child is not None:
                    # add to child prob
                    children_prob += child.weighted_probability
            # Once all children have been added
            # log(1/2,2) is equivalent to multiplying it all 1/2
            self.weighted_probability = (
                children_prob
                + math.log(1 + pow(2, (self.estimated_probability - children_prob)), 2)
                + math.log(1 / 2, 2))
        else:
            # If at the max-depth only use the estimated prob.
            # because there are no children available.
            self.weighted_probability = self.estimated_probability

    def print(self, tab):
        """
        Will printout the tree in it current state.
        :param tab: The character used as the space to show depth
        :return:
        """
        # Todo: Print self, based on the depth of the node and  children
        print(tab * (self.level - 1), self.symbol, "Prob: " + str(self.estimated_probability))

        for key, child in self.children.items():
            if child is not None:
                child.print(tab)
        return


class Tree:
    def __init__(self, name, depth, alpha = False):
        # The name of tree, used to identify it
        self.name = name
        # The maximum depth of the tree
        self.depth = depth
        # A boolean as to whether or not to use the
        # 0.05 alpha value for depths >11
        self.alpha = alpha

        # Create the root node of the tree
        # has special symbol: root
        self.root = Node("root", None, self, 0)

        # Table that holds all the nodes per depth
        # Is a dictionary of lists
        self.depthTable = defaultdict(list)
        self.depthTable[0] = [self.root]

        # This is used to store the changes in Pe as blocks are loaded in
        self.block_Pe = defaultdict(list)

        # A dictionary to store the fixed depth probilities
        self.fixed_depth_probabilities = dict()
        self.fixed_depth_probabilities_alpha = dict()

        for x in range (0, depth + 1):
            self.fixed_depth_probabilities[x] = 0
            self.fixed_depth_probabilities_alpha[x] = 0

    def read_block(self, block, inverted, competitive, previous_context):
        """
        Takes a block of symbols, generate symbol/context pairs and adds to the context
        :param block: a set of symbols that will be loaded into the tree
        :param inverted: A character that indicates whether or not to also add the inverted complement
        :param competitive: A boolean that indicates whether or not the depths be considered competitively
        :return: 4 values:
                    1 - Depth of max difference of the CTW
                    2 - The value of that difference (CTW)
                    3 - Depth of max difference of the fixed depth
                    4 - The value of the difference (fixed depth)
        """
        # Initialize a list used for the context
        context = previous_context
        # A dictionary used to store the previous root values
        previous_root_prob = dict()
        # A dictionary that stores the previous fixed depth prob
        previous_root_prob_Fixed_Depth = dict()

        # These two dicts will stores the value after the change
        root_prob = dict()
        Fixed_Depth_prob = dict()

        # start by going through all the depths
        for x in range(1, self.depth + 1):
            # Store the previous probabilities
            previous_root_prob[x] = self.root.weighted_probabilities[x]
            previous_root_prob_Fixed_Depth[x] = self.fixed_depth_probabilities[x]

            # stores the probabilities such that they can be accessed later
            self.block_Pe[x].append((previous_root_prob[x],previous_root_prob_Fixed_Depth[x]))

        # Now begin to loop through the block
        for index, symbol in enumerate(block):
            # If there is not enough symbols to go to full depth:
            if len(context) is not self.depth:
                #  Add to the context
                context.insert(0, symbol)
            else:
                #  Add to the tree
                self.add_data_point(symbol, context, inverted)

                # Remove a symbol, the last one added
                context.pop()
                # Insert a new one at the back (ie: the most recent one)
                context.insert(0, symbol)

        # Initialize dicts that stores the differences
        difference = dict()
        difference_Fixed_Depth = dict()

        # Loop over all the depths
        for x in range(1, self.depth + 1):
            # Store the updated probabilities
            root_prob[x] = self.root.weighted_probabilities[x]
            Fixed_Depth_prob[x] = self.fixed_depth_probabilities[x]

            # Calculate the differences
            difference[x] = root_prob[x] - previous_root_prob[x]
            difference_Fixed_Depth[x] = Fixed_Depth_prob[x] - previous_root_prob_Fixed_Depth[x]

        # If the output should be competitively chosen
        if competitive:

            # Return the values with the max differnce between before and after
            # MAX because the probability is in the log domain
            return max(difference, key=difference.get), \
                   difference[max(difference, key=difference.get)],\
                   max(difference_Fixed_Depth, key=difference_Fixed_Depth.get),\
                   difference_Fixed_Depth[max(difference_Fixed_Depth, key=difference_Fixed_Depth.get)],\
                   context
        else:
            # Otherwise return the value at the max depth
            return self.depth, \
                   difference[self.depth], \
                   max(difference_Fixed_Depth, key=difference_Fixed_Depth.get), \
                   difference_Fixed_Depth[max(difference_Fixed_Depth, key=difference_Fixed_Depth.get)],\
                   context

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

    def add_data_point(self, symbol, context, inverted):
        """
        Will take a symbol and context and add them to the tree
        also takes a char that determines if it should add the inverted complement
        :param symbol: The symbol to be added
        :param context: The context to that symbol
        :param inverted: A char or either N or Y that determines if the inverted
                        complement should also be added
        :return returns the change in probability from max depth
        """
        # Store the previous prob:
        previous_prob = self.root.weighted_probabilities[self.depth]

        # Store the dirty nodes return by traversing the tree
        dirty_nodes = self.traverse_context(symbol, context)
        # Re-weight the dirty nodes
        self.reweight_tree(dirty_nodes, symbol)

        # Add the inverted symbol/context pair if needed
        if inverted is 'Y':
            self.add_inverted_point(symbol, context)
        return self.root.weighted_probabilities[self.depth] - previous_prob

    def add_inverted_point(self, symbol, context):
        """
        Will take symbol and context, transform them into the inverted complement
        version and then add them to the tree
        :param symbol: The symbol to be inverted
        :param context: The context to be inverted
        """

        # Copy the context
        other_context = list(context)
        # Adds the symbol to the context
        other_context.insert(0, symbol)
        # Invert the context
        inverted = inverted_complement(other_context)
        # Take the first symbol as the new symbol to be added
        new_symbol = inverted.pop(0)

        # store the dirty nodes
        dirty_nodes = self.traverse_context(new_symbol, inverted)
        #reweight the dirty ndoes
        self.reweight_tree(dirty_nodes, new_symbol)

    def read_file(self):
        """
        Will open a file that is specified by the user and give them options
        on how to load it into the tree
        :return:
        """
        # Opens an interface to select a file

        p = cProfile.Profile()
        p.enable()
        filename = askopenfilename()

        # Shows the user the file-name
        print(filename)

        # gets the file extension
        file_extension = filename.split(".")[-1]
        # initialize a dict for storing the records
        records = dict()

        #If compressed, uncompress
        if file_extension == "gz":
            filename = gzip.open(filename, "rt")

        # initialize variables to find longest record
        longest = 0
        longest_id = None

        # checks for file-extension types
        if file_extension == "fa":
            for seq_record in SeqIO.parse(filename, "fasta"):
                if(len(seq_record) > longest):
                    longest_id = seq_record.id
                    longest = len(seq_record)
                print(seq_record.id, " ", len(seq_record.seq))
                records[seq_record.id] = seq_record.seq

        else:
            for seq_record in SeqIO.parse(filename, "gb"):
                if(len(seq_record) > longest):
                    longest_id = seq_record.id
                    longest = len(seq_record)
                records[seq_record.id] = seq_record.seq
                print(seq_record.id, " ", len(seq_record.seq))

        print("\n")
        file = input("What record do you want to open? [Input \"L\" for longest seqeunce]")
        if file == 'L':
            file_seq = records[longest_id]
            file = longest_id
        else:
            file_seq = records[file]

        print("Opened: " + file)
        print("Length: " + str(len(file_seq)))

        self.load_seq(file_seq, file)
        p.disable()
        p.dump_stats("CTW.profile")

    def load_seq(self,file_seq, file):
        """
        A function to load in a given file sequence into the tree
        :param file_seq: The sequence that will be added
        :param file: The file's name
        :return:
        """

        context = []
        read_type = input("How do you want to load the information: \n B : Blocks \n S : Sequentially \n")
        inverted = input("Do you want to add inverted complements?[Y/N]")
        competitive = input("Do you to load it competitively?[Y/N]")

        comp = 'Y' == competitive

        # Checks for read type
        if read_type is "B":
            block = list()
            block_seq = list()
            symbol_counter = 0

            max_symbols = int(input("How many symbols do you want to load?"))
            block_size = int(input("What block-size?"))

            # Generate a filename based on the parameters inputted by the user
            output_filename = file + "_S" + str(max_symbols) + "_BS" + str(block_size) + "_d" + str(self.depth) + "_Inv" + inverted + "_Comp" + str(competitive) + ".txt"

            previous_context = []
            for index, symbol in enumerate(file_seq):

                if symbol is not "N":

                    # Check to see if the desired amount of symbols has been loaded
                    if symbol_counter > max_symbols:
                        break

                    # Add to the block
                    block.append(symbol)
                    symbol_counter += 1

                    # If enough symbols in block, add it to the tree
                    if len(block) == block_size:
                        printProgressBar(symbol_counter, max_symbols)
                        block_seq.append(self.read_block(block, inverted, comp,previous_context))
                        previous_context = block_seq[len(block_seq) - 1][4]
                        del block[:]

            print('\n', "Length of block sequence", len(block_seq))

            # Print the results to a file
            print_to_file(block_seq, "data_graphs/data/CTW/CTW_" + output_filename, "data_graphs/data/Fixed_Depth/Fixed_Depth_" + output_filename)

        elif read_type is "S":
            output_file = []
            total_length = len(file_seq)
            for index, symbol in enumerate(file_seq):
                printProgressBar(index, total_length)
                if symbol is not 'N':
                    if len(context) is not self.depth:
                        #  Add to the context
                        context.insert(0, symbol)
                    else:
                        output_file.append((self.depth, self.add_data_point(symbol, context, inverted)))

                        context.pop()
                        context.insert(0, symbol)
            print_to_file(output_file, "data_graphs/data/sequentially/Sequential_"+ file + "_depth"+str(self.depth)+"_Inv" + str(inverted) + ".txt")


    def reweight_tree(self, dirtynodes, symbol):
        """
        Will re-weight all of the nodes affected by the addition of a new symbol/context pair

        :param dirtynodes: A list of all the nodes that have been changed by the recent symbol/context pair
        :param symbol: The symbol that was just added ( A C G or T)
        :return: N/a        """

        # For each node in the dirtnodes
        for node in dirtynodes:
            # Update the log probability
            node.update_probability_log(symbol)
        # And then update the probability of the root

    def traverse_context(self, symbol, context):
        """
        Will go through the tree, follow the path of the given context
        Will add nodes if they are not already present
        :param symbol: The new symbol that is being added to the tree
        :param context: The context of that symbol in the form of a list where
                        context[0] is the symbol directly before the symbol being added
        :return: A list of all the nodes either visited or created
        """

        # Start at the root
        current_node = self.root

        # A list to store all the nodes that have been accessed
        dirty_nodes = []

        # increment the symbol count of the current_node
        current_node.counts[symbol] += 1

        # Add the current node to the list of dirty nodes
        dirty_nodes.insert(0, current_node)

        # A counter to keep track of the depth
        depth = 1
        # For Each symbol in the context, traverse down from the root.
        # Begin by looping through the context
        for sym in context:


            # Check if the child to the current node exists
            if current_node.children[sym] is not None:
                # move to that node
                current_node = current_node.children[sym]
            else:
                # create a new node
                current_node.children[sym] = Node(sym, current_node, current_node.tree, depth)
                # self.depthTable[depth].append(current_node.children[sym])
                # Was commented out because was adding a lot of computation time
                # Move to the newly created node
                current_node = current_node.children[sym]
            # Increment depth
            depth += 1

            # increment the symbol count of the current_node
            current_node.counts[symbol] += 1

            # Add the current node to the list of dirty nodes
            dirty_nodes.insert(0, current_node)

        return dirty_nodes

    def print(self, tab):
        """
        Prints a representation of the Tree
        :param tab: the character used for spacing
        :return:
        """
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
            print("Loading file: ")
            self.read_file()

        elif command == "c" or command == "C":
            print_tree_commands()

        if command is "t":
            depth = input("what depth do you want to access? \n")
            if depth.isdigit():
                print(self.depthTable[int(depth)])

        if command is "e":
            for key, value in self.block_Pe.items():
                print_to_file(value, "outputs/ " + str(key) + "Pw and Fixed_Depth.mat")

        if command is "q":
            print("Weighted Probabilities  stored in root: (a = 0.5)")
            for key, value in self.root.weighted_probabilities.items():
                print("\t", key, ": ", value)

        if command is "z":
            print("Fixed Depth Probabilities stored in Tree: Alpha = 0.5")
            for key, value in self.fixed_depth_probabilities.items():
                print("\t", key, ": ", value)

            print("Fixed Depth Probabilities stored in Tree: Alpha = 0.05")
            for key, value in self.fixed_depth_probabilities_alpha.items():
                print("\t", key, ": ", value)


class Codon(Tree):
    """
    An extenstion of the Tree object
    """
    def __init__(self, name, depth, alpha = False):

        # The name of codon tree, used to identify it
        self.name = name
        # The maximum depth of the trees
        self.depth = depth

        # This is used to store the changes in Pe as blocks are loaded in
        self.block_Pe = defaultdict(list)

        # A dictionary to store the fixed depth probilities
        self.fixed_depth_probabilities = dict()
        self.fixed_depth_probabilities_alpha = dict()

        # Initialized the dictionaries
        for x in range(0, depth + 1):
            self.fixed_depth_probabilities[x] = 0
            self.fixed_depth_probabilities_alpha[x] = 0

        # Creates the three seperate trees that will be used to
        # encode the inputted symbols
        self.trees = {1 : Tree("Codon1",depth,alpha),
                      2 : Tree("Codon2",depth,alpha),
                      3 : Tree("Codon3",depth,alpha)}

        # Keeps track of the next tree that will be used
        self.current_tree = 1

        # The currently considered root
        self.root = self.trees[self.current_tree].root

    def read_block(self, block, inverted, competitive, previous_context):
        """
        Takes a block of symbols, generate symbol/context pairs and adds to the context
        :param block: a set of symbols that will be loaded into the tree
        :param inverted: A character that indicates whether or not to also add the inverted complement
        :param competitive: A boolean that indicates whether or not the depths be considered competitively
        :return: 4 values:
                    1 - Depth of the tree
                    2 - The value of that difference
        """

        context = previous_context
        previous_pw = 0

        # Add together the weighted probability of all three trees
        for keys, values in self.trees.items():
            previous_pw += values.root.weighted_probabilities[self.depth]

        # Loop through block
        for index, symbol in enumerate(block):
            # If not enough symbols in context
            if len(context) is not self.depth:
                #  Add to the context
                context.insert(0, symbol)
            else:
                # Access the current tree and add the data point
                self.trees[self.current_tree].add_data_point(symbol, context,inverted)
                #Increment the current tree
                self.current_tree += 1

                # If > 3 change to 1
                if self.current_tree is 4:
                    self.current_tree = 1

                # change the root
                self.root = self.trees[self.current_tree].root

                context.pop()
                context.insert(0, symbol)

        new_probability = 0

        for keys, values in self.trees.items():
            new_probability += values.root.weighted_probabilities[self.depth]

        return self.depth, new_probability - previous_pw, context

    def load_seq(self, file_seq, file):
        context = []
        read_type = input("How do you want to load the information: \n B : Blocks \n S : Sequentially \n")
        inverted = input("Do you want to add inverted complements?[Y/N]")

        block = list()
        block_seq = list()
        symbol_counter = 0

        if read_type is 'B':
            max_symbols = int(input("How many symbols do you want to load?"))
            block_size = int(input("What block-size?"))

            output_filename = file + "_MaxS" + str(max_symbols) + "_d" + str(self.depth) + "_bs" + str(block_size) +"_Inv" + inverted + ".txt"

            previous_context = []
            for index, symbol in enumerate(file_seq):

                if symbol is not "N":

                    if symbol_counter > max_symbols:
                        break

                    block.append(symbol)
                    symbol_counter += 1

                    if len(block) == int(block_size):
                        printProgressBar(symbol_counter, max_symbols)
                        block_seq.append(self.read_block(block, inverted, False, previous_context))
                        previous_context = block_seq[len(block_seq) - 1][2]
                        del block[:]

            print('\n', "Length of block sequence", len(block_seq))

            print_to_file(block_seq, "data_graphs/data/Codon_CTW_" + output_filename, "data_graphs/data/Codon_Fixed_Depth_" + output_filename)

        elif read_type is 'S':
            output_file = []
            total_length = len(file_seq)
            for index, symbol in enumerate(file_seq):
                printProgressBar(index, total_length)
                if symbol is not 'N':
                    if len(context) is not self.depth:
                        #  Add to the context
                        context.insert(0, symbol)
                    else:
                        output_file.append((self.depth, self.add_data_point(symbol, context, inverted)))
                        context.pop()
                        context.insert(0, symbol)

                        # Increment the current tree
                        self.current_tree += 1

                        # If > 3 change to 1
                        if self.current_tree is 4:
                            self.current_tree = 1

                        # change the root
                        self.root = self.trees[self.current_tree].root

            print_to_file(output_file, "data_graphs/data/sequentially/Codon_Sequential_" + file + "_depth" + str(
                self.depth) + "_Inv" + str(inverted) + ".txt")


def print_to_file(output, file1="CTW.txt", file2="Fixed_Depth.txt"):
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

    test_tree = Tree("Test7", 7)
    test_tree2 = Tree("Test16", 16, True)

    trees = dict()

    test = Tree("test", 16)
    #test_output = list()

    #print_trees(test.trees)
    #current_tree = test.tree

    trees[test_tree.name] = test_tree
    #trees[test.tree.name] = test.tree
    trees[test_tree2.name] = test_tree2

    # current_tree = current_tree.best_tree()
    current_tree = test_tree2

    # Table of commands
    print_commands()

    print_trees(trees)
    while 1:

        if current_tree is None:

            command = input('Enter your command: ')
            if command is "a":
                tree_name = input("Name of the new tree: ")
                depth = input("Depth:")
                type = input("Type? \n[N - Normal] \n[C- Codon]")
                if type is "N" or type is "n":
                    trees[tree_name] = Tree(tree_name, int(depth))
                elif type is 'C' or type is 'c':
                    trees[tree_name] = Codon(tree_name,int(depth))
                print_trees(trees)

            elif command is "o":
                name = input("\tWhich tree do you want to open?\n\t")
                current_tree = trees[name]
                print_tree_commands()

            elif command is "s":
                print_trees(trees)

            elif command is "d":
                name = input("\tWhich tree do you want to remove?\n\t")
                del trees[name]
                trees.pop(name, None)
                print_trees(trees)

            elif command is "c":
                print_commands()

        else:
            print("\n Currently in tree: ", current_tree.name)
            command = input("\t Enter your command: ")

            if command is "x" or command is "X":
                current_tree = None
                print_trees(trees)

            else:
                current_tree.commands(command)

if __name__ == "__main__": main()



