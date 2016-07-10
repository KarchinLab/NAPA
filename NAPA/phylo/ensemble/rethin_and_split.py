#!/usr/bin/python

# This script takes an input MrBayes treefile and creates an output treefile that
# contains a subset of the input trees sampled at regular intervals.   If a tree number
# is a multiple of the the specified interval, the tree is kept; otherwise it is discarded.
#
# This script currently does no error checking other than checking the command line and
# verifying file existence.  It is assumed that the input file is formatted correctly.
# Dave Swofford 20 April 2008
# modified to include burninSteps by Violeta Beleva Guthrie
# January 2014

import sys
from Bio import Phylo
def main():

 	if len(sys.argv) == 5:
  		infile_name = sys.argv[1]
 		thinning = int(sys.argv[2])
		burninIter = int(sys.argv[3]) # the burnin in number of iterations
		outputBase = sys.argv[4] # the output directory and base file name
		inputValid = thinning != 0
		print infile_name, thinning, burninIter
	else:
		inputValid = False
		print sys.argv

 	if (not inputValid):
 		scriptName = sys.argv[0].split('/')
 		scriptName = scriptName[len(scriptName) - 1]
 		print "Usage: %s inputtreefile thinningfrequency burninIter outputBase" % (scriptName)
 		return 1

	try:
		infile = open(infile_name, "rb")
	except:
		print "Could not open input treefile"
		return None

	header=""
	keepIter=False
	for line in infile:
		if (line.find("tree gen") < 0): header += line.replace('-','_')
		else:
			num = int(line.split("tree gen.")[-1].split(" = ")[0])
			keepIter = True if(num > burninIter and (num % thinning) == 0) else False
		if (keepIter):
			# write tree to individual file
			outFileName='/'.join([outputBase,outputBase])+'.'+str(num)+'.nex.tree'
			outfile = open(outFileName, 'wb')
			outfile.write(header+line+'end;\n')
			outfile.close()
			Phylo.convert(outFileName, 'nexus', outFileName.replace('.nex.','.nwk.'), 'newick')
			print "CONVERTED:", outFileName, "to NEWICK format"
	return 0

if __name__ == '__main__': main()
