#! /usr/bin/env python

import sys, dendropy

def scratch(trees):
	'''Take nexus file of trees and return list with list[0]
	being tree0 - treeN, list[1:] being symmetric distances
	between tree and all trees --> treeN.'''
	ts = dendropy.TreeList.get_from_path(trees,'nexus')
	L = []
	for i in range(len(ts)):
		L.append(tuple(["tree%i" % i]+\
			[dendropy.treecalc.symmetric_difference(ts[i], t) for t in ts[:i+1]]))
	return L
	
def symdiffs(infile):
	'''Call scratch(infile) and print table of symdiff values
	plus an average.'''
	vals = []
	L = scratch(infile)
	for row in L:
		print ''.join(str(i).rjust(7) for i in row)
		for symdiff in row[1:]:
			vals.append(symdiff)
	#print vals
	AVG = sum(vals)/float(len(vals))
	return AVG
	
if __name__ == '__main__':
	AVG = symdiffs(sys.argv[1])
	print "\nAverage symmetric difference: %i" % AVG
