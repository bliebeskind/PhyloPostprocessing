#! /usr/bin/env python

## Takes a mapped tree and dictionary file from phylip_map.py as 
# input and replaces the mapped taxon names (seq0, seq1..) with their
# original names.

Usage = """
tree_unmap.py <newick_tree> <mapping> 
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
                
def get_value(dict, key):
        """Returns dictionary value for key"""
        return dict[key]
               
tree = Phylo.read(infile, "nexus")
newtree = copy.deepcopy(tree) ## Copy input tree

#replace terminal names (keys: seq0...) with their values from input
#dictionary.
for clade in newtree.get_terminals():
        clade.name = get_value(get_mapping(mapping), clade.name.strip("'"))
                
Phylo.write([newtree], infilecore + "_unmap.nex", "nexus")
                                                                     
