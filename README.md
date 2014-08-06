PhyloPostprocessing
========

Collection of little scripts for doing stuff with phylogenetic trees.

Script | Useful for you? | Description |
------ | ----------------- | ----------- |
add_nhx_tags | no | add nhx style tags to nodes mapping a gene tree to species tree labels. Good for Notung, but is hardcoded for my data.
annotate_tree | yes | add PAML tags to tree nodes for branch models
print_accessions | no | for getting info on taxa in tree from upstream data
symdiff_compare | yes | for comparing topology of one tree (usually ML tree) with several (often MCMC trees). adapted from the [dendropy tutorial](http://packages.python.org/DendroPy/tutorial/trees.html)
symmdiff_table | yes | print table of all-by-all symmetric differences for a set of trees. Careful, not wise for lots of trees. Used in [Liebeskind et al. 2012](http://mbe.oxfordjournals.org/content/29/12/3613.full?sid=41b3a144-0227-4b66-9a30-97931870ca91)
tree_map | no | replace names in tree with seq1,seq2 etc...
unique_trees | maybe | Get unique topologies from an MCMC sample