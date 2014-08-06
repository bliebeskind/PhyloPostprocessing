#! /usr/bin/env python

import dendropy,sys
from dendropy import treecalc, TreeList, tree_source_iter

def tree_iter(infile,format,burnin,taxonset=None):
    '''Returns an iterator over a multi-tree file.'''
    return tree_source_iter(
    stream=open(infile, 'rU'),
    schema=format,
    taxon_set=taxonset,
    tree_offset=burnin)

def unique_trees(tree_list,mcmc_trees,format,burnin=0,taxonset=None):
    '''Takes a list and a Mr. Bayes mcmc sample as input.  Returns
    a list of non-redundant tree topologies using symmetric difference, 
    and the number of redundant topologies in the sample.'''
    redundant_count = 0
    for tree in tree_iter(mcmc_trees,format,burnin,taxonset):
    	for ut in tree_list:
    	    sd = treecalc.symmetric_difference(tree,ut)
            #print sd ## error check
            if sd == 0:
            	redundant_count +=1
                break
        else:
            tree_list.append(tree)
    return tree_list, redundant_count


if __name__ == '__main__':
    #inputs#
    mle_tree = raw_input("File with Maximum Likelihood tree: ")
    mcmc_trees = raw_input("File with MCMC trees: ")
    burnin = int(raw_input("Burnin: "))
    outfile = raw_input("Name of outfile: ")
    
    uts = [] #list of unique topologies
    taxa = dendropy.TaxonSet() #initialize TaxonSet object
    mle_tree = dendropy.Tree.get_from_path(mle_tree, 'nexus', taxon_set=taxa)
    uts.append(mle_tree) #MLE tree is the first topology in unique list
    
    uts, redundant_count = unique_trees(uts,mcmc_trees,'nexus',burnin,taxonset=taxa)
    print "\nNumber of redundant trees: %d" % redundant_count
    print "Number of unique trees: %d\n" % len(uts)
    unique_tree_list = TreeList(uts)
    unique_tree_list.write_to_path(outfile,'newick',suppress_edge_lengths=True)
	
    	    
