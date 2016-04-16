import argparse
from os import listdir
from os.path import isfile, join
from DirectedPairClass import *
from datetime import datetime
import sys
import itertools
from collections import defaultdict

posList = [ pos for pos in range(3, 33+1)] +\
          [pos for pos in range (33, 300+1) if pos not in [239,253] ]

temSeq="MSIQHFRVALIPFFAAFCLPVFAHPETLVK-VKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLL"
temSeq+="CGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELT"
temSeq+="AFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLR"
temSeq+="SALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW-----"

pos2wt = dict(zip(posList,temSeq))
subName = "TEM"

tip2pheno = {}
phenoClades = defaultdict(list)
phenotypeList = ["2b","2be","2ber"]

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
        parser.add_argument('-g', dest = "globalFasta", required=True,
                            help = "Path of fasta file  wih the MSA.")
        parser.add_argument('-p', dest = "leafPheno", required=True,
                            help="File where each known sequence's phenotype is stored")
        parser.add_argument('-o', dest = "outPath", required=True,
                            help = "Path for directed pairs tree summary.")
        parser.add_argument('-t', dest="distThresh", required=True,type=float,
                            help = "Distance threshold for phylo nodes.")
        args = parser.parse_args()
        return args


def get_internalNodes(nwkFile,fastaFile, gfStr):

    with open(fastaFile, 'rb') as fF: fString = fF.read()
    fStringOut = '\n'.join([gfStr.strip(),fString.replace('>edge.','>edge')])

    # use regexp to infer internal node names sequentially (as in PyCogent)
    with open(nwkFile, 'rb') as nF: nFstring = nF.read()
    cnt = itertools.count()
    outNwkStr = re.sub(r'\)\:', lambda x: ')edge{}:'.format(next(cnt)), nFstring)
    outNwkStr = re.sub(r'(.*)edge(.*)\:',r'\1root:',outNwkStr)
    
    fastaOutFile = fastaFile.replace('anc.prot.fasta','all.prot.fasta')
    with open(fastaOutFile, 'wb') as fOut:
        fOut.write(fStringOut)

    nwkOutFile = nwkFile.replace('nwk.tree','int_nwk.tree')
    with open(nwkOutFile, 'wb') as nOut:
        nOut.write(outNwkStr)

    return fastaOutFile, nwkOutFile

def get_dirPairs(nwkFile, fastaFile, outPath, distThresh):

    dirPairSet = DirectedPairSet(nwkFile, fastaFile, posList, pos2wt, 
                                 tip2pheno, phenoClades, distThresh)
    for dirPair in dirPairSet.dirPairs:
        dirPair.get_FisherPvals()

    stderr_write(["\t\tNumber of dir (mut) pairs", len(dirPairSet.dirPairs)])

    outDirFile = nwkFile.split('/')[-1]
    for ext in ['.tree','.nwk', '.int_nwk', '.int_subPheno','.txt','..','..']:
        outDirFile = outDirFile.replace(ext,'')
    outDirFile = outPath.rstrip('/')+'/' + outDirFile + '.'+str(distThresh)
    outDirFile += '.dirPairsPheno.txt'

    stderr_write(["\tWrite dirPairs to file",outDirFile.split('/')[-1],
                  datetime.now()])
    with open(outDirFile, 'wb') as eF:
        eF.write(str(dirPairSet))

    return dirPairSet

def get_allDirpairsRun(runDir,globalFasta, maxTrees, outPath, distThresh):
    with open(globalFasta, 'rb') as gF: msaStr = gF.read()
    treeFiles=sorted([ runDir+'/'+f for f in listdir(runDir) \
                       if f.strip().endswith('.nwk.tree')])
    ancProtFastas = sorted([ runDir+'/'+f for f in listdir(runDir) \
                             if f.strip().endswith('.nwk.tree.anc.prot.fasta')])

    if len(treeFiles) != len(ancProtFastas):
        print "Error: Tree and sequence file lists of different lengths."
        print len(treeFiles), len(ancProtFastas)

    if maxTrees<1: maxTrees = len(treeFiles) 
    for fileIdx in range(maxTrees):
        nwkFile = treeFiles[fileIdx]
        fastaFile = ancProtFastas[fileIdx]
        fastaFileOut, nwkFileOut = get_internalNodes(nwkFile,fastaFile,msaStr)
        #nwkFileOut = nwkFileOut.replace('int_nwk.','int_subPheno.')
        stderr_write(["\nFASTA/nwk:",fastaFileOut.split('/')[-1], '\t', 
                      nwkFileOut.split('/')[-1]])
        dirPairSet = get_dirPairs(nwkFileOut, fastaFileOut, outPath, distThresh)

def get_allDirPairs(runDir, maxRuns, globalFasta, maxTrees, tipPhenoFile, outPath,
                    distThresh):
    global tip2pheno
    global phenoClades

    phenoDict = {'1':'2br','2':'2b','3':'2ber','4':'2be'}
    with open(tipPhenoFile,'rb') as pf:
        lineRecs = [line.strip().split() \
                    for line in pf.read().replace('"','').split('\n')]
        for rec in lineRecs:
            if len(rec)>1:
                if rec[1].strip() not in phenoDict: continue
                tip2pheno[rec[0]]=phenoDict[rec[1].strip()]

    for tip in tip2pheno:
        if tip2pheno[tip] in ["2be","2ber"]:
            phenoClades[tip2pheno[tip]].append(tip)

    runNum = maxRuns
    runRunDir = runDir+'/'+'run'+str(runNum)
    get_allDirpairsRun(runRunDir,globalFasta, maxTrees, outPath, distThresh)

def main():
    args = parseArgs()
    distThresh = args.distThresh
    get_allDirPairs(args.runDir, args.maxRuns, args.globalFasta, 
                    args.maxTrees, args.leafPheno, args.outPath, distThresh)
    
if __name__ == "__main__": main()

        
