#! /usr/bin/env python

import dendropy

exceptions = {'sb': 'Ctenophora', 'ML': 'Ctenophora','c0': 'Cnidaria','GE': 'Cnidaria'}

def get_mapping(map_file):
	'''Opens file with mapping info and returns dictionary'''
	with open(map_file) as mapping:
		my_map = eval(mapping.read().strip('\n'))
		return my_map
		
def annotate_branches(tree_file,map_file,format='newick'):
	tree = dendropy.Tree.get_from_path(
		tree_file,format,preserve_underscores=True,extract_comment_metadata=True)
	mapping = get_mapping(map_file)
	for i in tree.postorder_node_iter():
		if i.is_leaf():
			try:
				group = mapping[str(i.taxon)[:4]]
			except KeyError:
				group = exceptions[str(i.taxon)[:2]]
			i.annotations.add_new(name='taxon',value="%s" % group)
		else: # is internal
			func = lambda x,y: x if x == y else False
			child_groups = (
				c.annotations.get_value('taxon')
				for c in i.child_nodes()) # generator of child taxa
			group = reduce(func,child_groups) # returns group if all children match, else False
			if group:
				i.annotations.add_new(name='taxon',value="%s" % group)
	return tree
				
	
	
