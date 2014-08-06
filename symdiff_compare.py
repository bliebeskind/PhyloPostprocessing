#! /usr/bin/env python
#Adapted from Dendropy tutorial: http://packages.python.org/DendroPy/tutorial/trees.html

##### Note: both input tree files must be in nexus format and have the	#####
##### same taxon list.				      								#####

import dendropy,sys
from dendropy import tree_source_iter
from dendropy import treecalc


def get_mle_tree(mle_tree):
	taxa = dendropy.TaxonSet()
	mle_tree = dendropy.Tree.get_from_path(mle_tree, 'nexus', taxon_set=taxa)
	return mle_tree, taxa
	
def print_distances(tree_list,mle,uniq_flag=False):
	mle_tree, taxa = get_mle_tree(mle)
	distances = []
	uniq_trees = dendropy.TreeList()
	count = 1
	for t in tree_source_iter(stream=open(tree_list, 'rU'),schema='nexus',taxon_set=taxa):
		dist = treecalc.symmetric_difference(mle_tree, t)
		print "Distance between MLE tree and tree %i: %i" % (count,dist)
		distances.append(dist)
		count +=1
		if uniq_flag and dist > 0:
			uniq_trees.append(t)
	print("Mean symmetric distance between MLE and tree list: %d" \
		% float(sum(distances)/len(distances)))
	return uniq_trees, len(uniq_trees)
        
if __name__ == '__main__':
	Usage = '''
symdist.py <MLE tree> <TreeList> (optional args -->) <--U> <outfile>
## '--U' will print unique topologies to 'outfile' ##\n'''
	mletree = sys.argv[1]
	treelist = sys.argv[2]
	if len(sys.argv) > 3 and sys.argv[3] == '--U':
		out = sys.argv[4]
		utrees, num = print_distances(treelist, mletree, uniq_flag=True)
		dendropy.TreeList.write_to_path(utrees, out,'newick', suppress_edge_lengths=True)
		sys.stderr.write("Found %i unique trees\n" % num)
	else:
		print_distances(treelist, mletree)
	
