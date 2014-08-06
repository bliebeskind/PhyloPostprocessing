#! /usr/bin/env python

## Adapted from BioPython Cookbook 6.2.1
# Takes a fasta alignment as input and replaces names with seq0, seq1...
# and converts to phylip.

Usage = """
tree_map.py <newick_tree> <mapping> 
"""

from Bio import Phylo
import sys                
import copy

infile = sys.argv[1]
infilecore = infile.split('.')[0]
mapping = sys.argv[2]

def get_mapping(infile):  
        """Opens file with mapping info and returns dictionary"""
        with open(infile) as map:
                my_map = eval(map.read().strip('\n'))
                return my_map   
                
def get_key(dict, value):
        """Returns key from a dictionary given a value"""
        return [key for key, val in dict.iteritems() if val == value][0]
                                             
tree = Phylo.read(infile, "newick")
newtree = copy.deepcopy(tree) ## Copy input tree

#replace terminal names with their keys (seq0, seq1 ...) from input
#dictionary
for clade in newtree.get_terminals():
        clade.name = get_key(get_mapping(mapping), clade.name.strip("'"))
                
Phylo.write([newtree], infilecore + "_map.nhx", "newick")
                                                                     
