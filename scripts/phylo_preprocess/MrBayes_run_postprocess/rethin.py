#!/usr/bin/python

# This script takes an input MrBayes treefile and creates an output treefile that
# contains a subset of the input trees sampled at regular intervals.   If a tree number
# is a multiple of the the specified interval, the tree is kept; otherwise it is discarded.
#
# This script currently does no error checking other than checking the command line and
# verifying file existence.  It is assumed that the input file is formatted correctly.
#
# Dave Swofford 20 April 2008
# modified to include burninSteps by Violeta Beleva Guthrie
# January 2014

import sys
def main():

 	if len(sys.argv) == 4:
  		infile_name = sys.argv[1]
 		thinning = int(sys.argv[2])
		burninIter = int(sys.argv[3]) # the burnin in number of iterations
		inputValid = thinning != 0
		print infile_name, thinning, burninIter
	else:
		inputValid = False
		print sys.argv

 	if (not inputValid):
 		scriptName = sys.argv[0].split('/')
 		scriptName = scriptName[len(scriptName) - 1]
 		print "Usage: %s inputtreefile thinningfrequency burninIter > outputtreefile" % (scriptName)
 		return 1

	try:
		infile = open(infile_name, "r")
	except:
		print "Could not open input treefile"
		return None

	for line in infile:
		if (line.find("tree gen") >= 0):
			num = int(line.split("tree gen.")[-1].split(" = ")[0])
			keepLine = True if(num > burninIter and (num % thinning) == 0) else False
		else:
			keepLine = True
		if (keepLine):
			print line,

	return 0

if __name__ == '__main__': main()
