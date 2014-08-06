#! /usr/bin/env python

## Currently works except for dendropy's behavior of quoting annotations for
## distal branches.  Can't figure out how to suppress this behavior
## Could do a post processing hack like this to remove quotes:
## sed "s/'//g" my_annotated.tre >> my_clean_annotated.tre

## Annotates branches of trees to be used as input for PAML branch models. This 
## script puts these annotations there by reading in a file whose syntax tells
## the program where to annotate.  The syntax of the infile is as follows:

##	[Outgroup]
##	taxon1, taxon2, #1, c
##	taxon3, taxon4, $2, b

## The two taxon names are used to get a most recent common ancestor (mrca).
## The third column is the annotation, and the fourth tells the program
## whether to put the annotation on the whole clade (minus the mrca) = "c"
## or whether to put the annotation just on the mrca = "b"
## The code above puts a '#1' on all edges subtended by the common ancestor
## of taxon1 and taxon2, but not this ancestral branch. And it puts a $2 on
## the mrca of taxon3 and taxon4, but not the descendent branches.

## Always double check that these annotations 
## have been put in the right place by opening the tree file in FigTree.

import sys, dendropy

def parse_infile(taxafile):
	'''Take a file of taxon names and annotations symbols, separated by
	commas, and returns a list of tuples.'''
	taxa_list = []
	with open(taxafile) as f:
		line = f.readline()
		while line:
			line = line.strip()
			while line.startswith('['):
				start = line.find('[') + 1
				outgroup = line[start: line.find(']')]
				line = f.readline()
			line = line.split(',')
			assert len(line) == 4, \
			"Format problem in taxa file: each line should have just two"\
			" taxon names and a branch lable each separated by a comma"
			#dedropy has a 'bug' that replaces a "_" in a taxon name with a space.
			taxa_list.append(tuple([taxon.strip().replace("_"," ") for taxon in line]))
			line = f.readline()
	return taxa_list, outgroup
	## should have some output showing what was read in?
	
def annotate_children(tree, taxa, label):
	'''Given a tree and two taxa, do a postorder trace from mrca of both 
	taxa, labeling each edge except mrca with label. Return tree'''
	anc = tree.mrca(taxon_labels=taxa)
	for child in anc.postorder_iter():
		child.label = "%s" % label
	anc.label = ''
	return tree

def root_on_outgroup(tree, instring):
	'''Takes dendropy tree object, and string of comma separated taxa. Roots
	tree on mrca of specified taxa if more than one in string, otherwise roots 
	on edge associated with taxon.  Returns tree object'''
	string = instring.replace(' ','') #remove spaces
	outgroups = string.split(',')
	assert len(outgroups) > 0, "No outgroups specified"
	if len(outgroups) > 1: # is there more than one taxon specified?
		anc = tree.mrca(taxon_labels=outgroups)
	elif len(outgroups) == 1:
		anc = tree.find_node_with_taxon_label(outgroups[0])
	if anc.parent_node == None:
			tree.is_rooted = True
	else:
		tree.reroot_at_edge(anc.edge, update_splits=True)
	return tree
	
def annotate_tree(tree, taxafile):
	'''Given a dendropy tree obeject ...
	number the branches that are ancestral to the two taxa with $0, $1 etc.
	Outputs a tree (default name and format are "annotated_tree.tre", nexus).
	Requires dendropy and is made to used with taxa_list().'''
	tree.is_rooted = True
	#count = 0
	taxalist, outgroup = parse_infile(taxafile)
	tree = root_on_outgroup(tree, outgroup)
	for t in taxalist:
		try:
			taxa = t[:2]
			label, clade = t[2], t[3]
			if clade.lower() == 'b':
				anc = tree.mrca(taxon_labels=t[:2])
				anc.label = "%s" % label
				#count +=1
			elif clade.lower() == 'c':
				tree = annotate_children(tree, taxa, label)
				#count +=1
		except KeyError, detail:
			sys.stderr.write("%s: %s\n" % (str(detail), str(t)))
	return tree
	
	
if __name__ == '__main__':
	Usage = "\nUsage: annotate_branches.py <newick treefile> <taxa file>" \
	" <Optional: outfile>\n"\
	"### Taxa file should have just two comma separated taxa per line ###\n"
	try:
		treefile = sys.argv[1]
		tree = dendropy.Tree.get_from_path(treefile,'newick',preserve_underscores=True)
		taxafile = sys.argv[2]
		if len(sys.argv) > 3:
			outfile = sys.argv[3]
			annotate_tree(tree, taxafile)
			tree.write_to_path(outfile,'newick',preserve_spaces=True, suppress_leaf_node_labels=False, node_label_element_separator=' ')
		else:
			annotate_tree(tree, taxafile)
			tree.write_to_path('annotated.nhx','newick',preserve_spaces=True, suppress_leaf_node_labels=False, node_label_element_separator=' ')
	except IndexError:
		print Usage	
