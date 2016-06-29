import argparse
import numpy as np
from os import listdir
from os.path import isfile, join
from UndirectedPairClass import *
from datetime import datetime
import sys
import itertools
from collections import defaultdict

getInt = lambda s: int(re.sub('\D','',s))

posList = [ pos for pos in range(3, 33+1)] +\
          [pos for pos in range (33, 300+1) if pos not in [239,253] ]

temSeq="MSIQHFRVALIPFFAAFCLPVFAHPETLVK-VKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLL"
temSeq+="CGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELT"
temSeq+="AFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLR"
temSeq+="SALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW-----"

pos2wt = dict(zip(posList,temSeq))
subName = "TEM"

def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

def parseArgs():
        parser = argparse.ArgumentParser()
        parser.add_argument('-d', dest = "runDir", required=True, 
                            help = "Directory with nwk trees and ancestral sequences")
        parser.add_argument('-r', dest = 'maxRuns', required=True, type=int,
                            help = "Maximum number of independent runs to include.")
        parser.add_argument('-m', dest = "maxTrees", type=int,
                            help = "The maximum number of trees to analyze."+\
                            " Default = all.")
        parser.add_argument('-t', dest = "pvalThresh", required=True, type=float,
                            help = "The p-value threshold to apply to Fisher's exact")
        parser.add_argument('-dt', dest="distThresh", required=True,type=float,
                            help = "Distance threshold for phylo nodes.")
        parser.add_argument('-o', dest = "outPath", required=True,
                            help = "Path for undirected pairs tree summary.")
        args = parser.parse_args()
        return args

def get_allUndirPairsRun(runDir, maxTrees, undirPairs2counts, negPairs2counts, 
                         distThresh, pvalThresh):
    stderr_write(["\tAdding counts from:", runDir.split('/')[-1]])
    suffix = str(distThresh)+'.undirPairsPheno.txt'
    undirFiles=sorted([ runDir+'/'+f for f in listdir(runDir) \
                       if f.strip().endswith(suffix)])
    if maxTrees<1: maxTrees = len(undirFiles) 
    for fileIdx in range(maxTrees):
        undirPairFile = undirFiles[fileIdx]
        with open(undirPairFile, 'rb') as epF:
            for line in epF:
                if line.startswith('_'):continue
                [mut0, mut1] = line.strip().split('\t')[0:2]
                if getInt(mut0) < 23 or getInt(mut1) < 23:
                    continue
                [pvalLess, pvalMore] = [eval(p) for p in line.strip().split('\t')[4:6]]
                if pvalMore <  pvalThresh:
                    undirPairs2counts[(mut0,mut1)] += 1
                elif pvalLess < pvalThresh:
                    negPairs2counts[(mut0,mut1)] += 1
    return undirPairs2counts,negPairs2counts
                           
def get_allUndirPairs(runDir, maxRuns, maxTrees, outPath, distThresh, pvalThresh):
    undirPairs2counts, negPairs2counts = defaultdict(int), defaultdict(int)
    for runNum in range(1,maxRuns+1):
        stderr_write(["Collecting pair counts from run#", runNum])
        runRunDir = runDir+'/'+str(runNum)
        undirPairs2counts, negPairs2counts = get_allUndirPairsRun(runRunDir, 
                                            maxTrees,undirPairs2counts,
                                            negPairs2counts,
                                            distThresh, pvalThresh)
        
    outFile = outPath.rstrip('/')+'/'+str(datetime.now()).split()[0]+\
                 '.undirPhenoPair'+str(pvalThresh)+'.dist'+str(distThresh)+'.txt'
    with open(outFile, 'wb') as of:
        of.write('\t'.join(['_Mut1_','_Mut2_','treeCount'])+'\n')
        for pair in undirPairs2counts:
            of.write('\t'.join(list(pair)+[str(undirPairs2counts[pair])])+'\n')

    outFile = outPath.rstrip('/')+'/'+str(datetime.now()).split()[0]+\
                 '.negPhenoPair'+str(pvalThresh)+'.dist'+str(distThresh)+'.txt'
    with open(outFile, 'wb') as of:
        of.write('\t'.join(['_Mut1_','_Mut2_','treeCount'])+'\n')
        for pair in negPairs2counts:
            of.write('\t'.join(list(pair)+[str(negPairs2counts[pair])])+'\n')

def main():
    args = parseArgs()
    get_allUndirPairs(args.runDir, args.maxRuns, args.maxTrees, args.outPath, 
                      args.distThresh, args.pvalThresh)

if __name__ == "__main__": main()

        
