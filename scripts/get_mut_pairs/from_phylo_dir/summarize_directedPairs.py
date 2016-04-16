import argparse
import numpy as np
from os import listdir
from os.path import isfile, join
from DirectedPairClass import *
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
        parser.add_argument('-s', dest = "distThresh", required=True, type=float,
                            help = "The distance threshold")
        parser.add_argument('-o', dest = "outPath", required=True,
                            help = "Path for directed pairs tree summary.")
        args = parser.parse_args()
        return args

class MutPair(object):
    def __init__(self, mut1, mut2):
        self.mut1, self.mut2 = mut1, mut2
        self.pair = (mut1, mut2)
        self.allowed12 = 0 # count of trees allowing directed pair mut1,mut2
        self.disallowed12 = 0
        self.allowed21 = 0
        self.disallowed21 = 0
        self.directionality = "unknown"
    
    def __repr__(self):
        return '\t'.join([self.mut1, self.mut2, str(self.allowed12),
                          str(self.disallowed12), str(self.allowed21), 
                          str(self.disallowed21), self.directionality])

def get_dirPairs(posKeys,negKeys):
    posPairList = []
    for mutPair in posKeys+negKeys:
        if getInt(mutPair[0]) > getInt(mutPair[1]):
            posPairList.append((mutPair[1],mutPair[0]))
        else:
            posPairList.append(mutPair)
    return set(posPairList)


def write_summary(outFile, posDirPairs2counts, negDirPairs2counts, pvalThresh):
    mutPairObjects = []
    mutPairSet = get_dirPairs(posDirPairs2counts.keys(),
                               negDirPairs2counts.keys())
    for pair in mutPairSet:
        mut1,mut2 = pair[0],pair[1]
        mutPairObject = MutPair(mut1,mut2)
        if pair in posDirPairs2counts:
            mutPairObject.allowed12 += posDirPairs2counts[pair]
        if pair in negDirPairs2counts:
            mutPairObject.disallowed12 += negDirPairs2counts[pair]
        # check for reversed pair
        if (mut2, mut1) in posDirPairs2counts:
            mutPairObject.allowed21 += posDirPairs2counts[(mut2,mut1)]
        if (mut2,mut1) in negDirPairs2counts:
            mutPairObject.disallowed21 += negDirPairs2counts[(mut2,mut1)]
        mutPairObjects.append(mutPairObject)
    
    mutPairObjects.sort(key=lambda x: (getInt(x.mut1),getInt(x.mut2)))
    with open(outFile, 'wb') as of:
        of.write('\t'.join(['mut1','mut2','incr12','decr12','incr21','decr21'])\
                 +'\n')
        for m in mutPairObjects:
            of.write(str(m)+'\n')

def get_allDirPairsRun(runDir, maxTrees, posDirPairs2counts, negDirPairs2counts,
                       distThresh,pvalThresh):
    stderr_write(["\tAdding counts from:", runDir])

    dirFiles=sorted([ runDir+'/'+f for f in listdir(runDir) \
                       if f.strip().endswith(str(distThresh)+'.dirPairsPheno.txt')])
    stderr_write(["Number of files from run(s)", len(dirFiles)])
    if maxTrees<1: maxTrees = len(dirFiles) 
    for fileIdx in range(maxTrees):
        dirPairFile = dirFiles[fileIdx]
        with open(dirPairFile, 'rb') as epF:
            for line in epF:
                if line.startswith('_'):continue
                [mut1, mut2] = line.strip().split('\t')[0:2]
                if getInt(mut1) < 23 or getInt(mut2) < 23:
                    continue
                [pValLessFDR, pValMoreFDR] = line.strip().split('\t')[-2:]
                if float(pValLessFDR) < pvalThresh: 
                    negDirPairs2counts[(mut1,mut2)]+= 1 
                if float(pValMoreFDR) < pvalThresh: 
                    posDirPairs2counts[(mut1,mut2)]+= 1 
    return posDirPairs2counts, negDirPairs2counts

                           
def get_allDirPairs(runDir, maxRuns, maxTrees, outPath, distThresh,pvalThresh):
    posDirPairs2counts, negDirPairs2counts = defaultdict(int),defaultdict(int)
    for runNum in range(1,maxRuns+1):
        stderr_write(["Collecting pair counts from run#", runNum])
        runRunDir = runDir+'/'+str(runNum)
        posDirPairs2counts, negDirPairs2counts = \
                                get_allDirPairsRun(runRunDir, maxTrees, 
                                                   posDirPairs2counts, 
                                                   negDirPairs2counts,
                                                   distThresh,pvalThresh)
        
    outFile = outPath.rstrip('/')+'/'+str(datetime.now()).split()[0]+\
                 '.directPhenoPair_p'+str(pvalThresh)+'_d'+str(distThresh)+\
                 '.txt'
    write_summary(outFile, posDirPairs2counts, negDirPairs2counts, 
                  pvalThresh)

def main():
    args = parseArgs()
    get_allDirPairs(args.runDir, args.maxRuns, args.maxTrees, args.outPath, 
                    args.distThresh,args.pvalThresh)
    
if __name__ == "__main__":
    main()

        
