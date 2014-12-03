#! /usr/bin/env python

import dendropy

exceptions = {'sb': 'Ctenophora',
	'ML': 'Ctenophora',
	'c0': 'Cnidaria',
	'GE': 'Cnidaria'}

def get_mapping(map_file):
	'''Opens file with mapping info and returns dictionary'''
	with open(map_file) as mapping:
		my_map = eval(mapping.read().strip('\n'))
		return my_map
		
def annotate_branches(tree_file,taxon_mapping,attribute_mapping=None,format='newick'):
	tree = dendropy.Tree.get_from_path(
		tree_file,format,preserve_underscores=True,extract_comment_metadata=True)
	for i in tree.postorder_node_iter():
		if i.is_leaf():
			try:
				group = taxon_mapping[str(i.taxon)[:4]]
			except KeyError:
				group = exceptions[str(i.taxon)[:2]]
			i.annotations.add_new(name='taxon',value="%s" % group)
			if attribute_mapping:
				for n,v in attribute_mapping[group].iteritems():
					i.annotations.add_new(name=n,value=v)
		else: # is internal
			func = lambda x,y: x if x == y else False
			child_groups = (
				c.annotations.get_value('taxon')
				for c in i.child_nodes()) # generator of child taxa
			group = reduce(func,child_groups) # group if all children match, else False
			if group: # if all children are same taxon
				i.annotations.add_new(name='taxon',value="%s" % group)
				if attribute_mapping:
					for n,v in attribute_mapping[group].iteritems():
						i.annotations.add_new(name=n,value=v)
	return tree
				
	
	
