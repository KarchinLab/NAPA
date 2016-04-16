from cogent import LoadTree, LoadSeqs, DNA
from cogent.evolve import models
from cogent.evolve.models import GTR, JC69

def optimize_params(treeFile, sequenceFile):
	#treeFileDir ='/'.join( treeFile.split('/')[0:-1] ) + '/'
	#treeFileName = '.'.join( treeFile.split('/')[-1].split('.')[0:-1] )
	#treeFileExt = treeFile.split('/')[-1].split('.')[-1]

	tree = LoadTree(treeFile)
	tree.root().Name='root'
	aln = LoadSeqs(sequenceFile,moltype=DNA,aligned=True)
    
    # substitution model definition
    # see /Library/Python/2.7/site-packages/cogent/evolve/models.py
	#subMod = GTR(with_rate=True, distribution='gamma') 
	subMod = JC69() # for testing purposes
	print subMod
	
	# likelihood function based on substitution model
	lf = subMod.makeLikelihoodFunction(tree)
	lf.setAlignment(aln)
	print lf # print model
	
	
	# optimization step
	# use simmulated annealing first and then the more rigorous local Powell
	# use checkpointing in case of a crash
	print "Optimizing LF"
	lf.optimise(global_tolerance=2.0,tolerance=1e-6,max_restarts=5)

	# lf statistics
	tables = lf.getStatistics(with_motif_probs=True, with_titles=True)
	for table in tables: print table

	# display optimized log-likelihood
	lnL=lf.getLogLikelihood()
	print "Optimized ln(L):", lnL

    # Does not work for some reason
	#annot_tree = lf.getAnnotatedTree() # tree annotated by ancestor
	#annot_tree.writeToFile(treeFileDir+treeFileName+".annot."+treeFileExt)

	return lf

def ancestralGivenLF(lf):
    """ get ancestral sequences given optimized likelihood function """
    # get the most likely ancestral sequences
    ancestors = lf.likelyAncestralSeqs()
    return ancestors

def ancestralProbs(lf, seqId):
	ancestral_probs = lf.reconstructAncestralSeqs()
	if seqId in ancestral_probs: 
		print ancestral_probs[seqId]
	return ancestral_probs[seqId]

def main():
	treeFile="test.nwk.tree"

	sequenceFile="test.fasta" # make sure tip names match

	lf = optimize_params(treeFile, sequenceFile)
	ancestors = ancestralGivenLF(lf)
	print type(ancestors)
	print str(ancestors)
	outanc = treeFile+".anc.seq"
	oa = open(outanc,'wb')
	oa.write(str(ancestors))
		
if __name__ == "__main__": main() 
    
