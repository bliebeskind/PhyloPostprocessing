#! /usr/bin/env python


from Bio import Phylo
import sys                
import copy

## This is hard coded - not desirable.
exceptions = {'sb': 'Pleurobrachia', 'ML': 'Mnemiopsis','c0': 'Acropora','GE': 'Porites'}

def get_mapping(map_file):
	'''Opens file with mapping info and returns dictionary'''
	with open(map_file) as mapping:
		my_map = eval(mapping.read().strip('\n'))
		return my_map
		
def add_labels(infile,map_file,format="newick"):
	'''Read in gene tree and dictionary file to map genes to species. 
	Add NHX comment section with species information to tips.'''
	infilecore = infile.split('.')[0]
	tree = Phylo.read(infile, format)
	newtree = copy.deepcopy(tree) ## Copy input tree
	mapping = get_mapping(map_file)
	for clade in newtree.get_terminals():
		if clade.name.strip("'")[:2] in exceptions:
			key = clade.name.strip("'")[:2]
			species = exceptions[key]
		else:
			key = clade.name.strip("'")[:4]
			species = mapping[key]
		clade.name = str(clade.name) + "[&&NHX:S=%s]" % species
	Phylo.write([newtree], infilecore + "_tagged.nhx", "newick")
	
if __name__ == '__main__':
	Usage = """
		For adding species tags to a tree in the NHX format:
		[&&NHX:S=species]

		Will use a text dictionary as an infile to map gene tree
		names to species information:

		add_nhx_tags.py <newick_tree> <mapping> <format>
		"""
	try:
		infile = sys.argv[1]
		mapping = sys.argv[2]
		format = sys.argv[3]
	except IndexError:
		print Usage
	add_labels(infile,mapping,format)
