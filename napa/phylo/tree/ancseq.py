import argparse
from cogent import LoadTree, LoadSeqs, DNA
from cogent.evolve import models
from cogent.evolve.models import GTR, JC69 # these are DNA evolution only

#Current implementation lacks options for protein vs. nucleotide models etc.
#http://pycogent.org/cookbook/using_likelihood_to_perform_evolutionary_analyses.html

def argParser():
        parser = argparse.ArgumentParser()
        parser.add_argument('-t', dest="treeFile", 
                            help="Phylogenetic tree in newick format: *.nwk.tree",
                            required=True)
        parser.add_argument('-s', dest="sequenceFile",
                            help="Sequence file in fasta format (aligned).",
                            required=True)
        parser.add_argument('-a', dest="ancSeqFile", 
                            help="Ancestral sequences are written to this file.",
                            required=True)        
        args = parser.parse_args()
        return args

def optimize_params(treeFileDir, treeFileName, sequenceFile):
        '''Create an ancestral likelihood model in pycogent and
        optimize it. Currently implemented for nucleotides.'''

	tree = LoadTree(treeFile)
	tree.root().Name='root'
	aln = LoadSeqs(sequenceFile,moltype=DNA,aligned=True)
    
        # substitution model definition
        # see /Library/Python/2.7/site-packages/cogent/evolve/models.py
	subMod = GTR(with_rate=True, distribution='gamma') 
	#subMod = JC69() # for testing purposes
	print subMod
	
	# likelihood function based on substitution model
	lf = subMod.makeLikelihoodFunction(tree)
	lf.setAlignment(aln)
	print lf # print model
	
	# optimization step
	# uses simmulated annealing first and then the more rigorous local Powell
	
	print "Optimizing Likelihood Function"
	checkpointLF = treeFileName+".anc.LF.chk" # use checkpointing in case of a crash
	lf.optimise(filename=checkpointLF, global_tolerance=2.0, tolerance=1e-6,
                    max_restarts=5) # should allow user to choose params 

	# likelihood function statistics
	tables = lf.getStatistics(with_motif_probs=True, with_titles=True)
	for table in tables: print table

	# display optimized log-likelihood
	lnL=lf.getLogLikelihood()
	print "Optimized ln(L):", lnL

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
        args=argParser()
        
	treeFile=args.treeFile 

	treeFileDir ='/'.join( treeFile.split('/')[0:-1] ) + '/'
	treeFileName = treeFile.split('/')[-1].split('.')

	sequenceFile=args.sequenceFile #"PSE4-TEM-SHV_noLen.cds.fasta" # make sure tip names match

	lf = optimize_params(treeFileDir, treeFileName, sequenceFile)
	ancestors = ancestralGivenLF(lf)
	outanc = args.treeFile + ".anc.seq"
	
        with open(args.ancSeqFile,'wb') as f:
                f.write(str(ancestors))
		
if __name__ == "__main__": main() 
    
