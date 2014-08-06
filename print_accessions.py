#! /usr/bin/env python

## Should print accession numbers, database info, taxon name in spreadsheet
## given a tree of taxon names

## ! Depends on tree_unmap.py ! ##

from Bio import Phylo
from tree_unmap import get_mapping, get_value
import sys, csv, re

def get_mapping(infile):  
	"""Opens file with mapping info and returns dictionary"""
	with open(infile) as map:
		my_map = eval(map.read().strip('\n'))
		return my_map   
                
def get_value(dict, key):
	"""Returns dictionary value for key"""
	return dict[key]

def taxon_names(treefile, format="nexus"):
	'''Return list of terminal names given input tree file and format.
	Treefile should have just one tree in it.'''
	tree = Phylo.read(treefile, format)
	taxa = []
	for clade in tree.get_terminals():
		taxon = re.sub('_D\d','',clade.name) #strip domain names
		if taxon not in taxa: #Add only one instance of 4-domain types
			taxa.append(taxon)
	return sorted(taxa)
	
def get_info(taxa, mapping):
	'''Yield tuple of taxon and accession'''
	D = get_mapping(mapping)
	for taxon in taxa:
		try:
			accession = get_value(D, taxon)
			yield taxon, accession
		except KeyError:
			sys.stderr.write("Taxon: %s, not in mapping\n" % taxon)
			
if __name__ == '__main__':
	try:
		TREEFILE = sys.argv[1]
		MAPPING = sys.argv[2]
		csvout = csv.writer(sys.stdout)
		taxa = taxon_names(TREEFILE)
		for tup in get_info(taxa, MAPPING):
			csvout.writerow(tup)
	except IndexError:
		sys.stderr.write("Index Error: ", detail)
