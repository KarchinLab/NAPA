import argparse
from os import listdir
from os.path import isfile, join
from UndirectedPairClass import *
from datetime import datetime
import sys
import itertools
from collections import defaultdict
from cogent import LoadTree, LoadSeqs, DNA

posList = [ pos for pos in range(3, 33+1)] +\
          [pos for pos in range (33, 300+1) if pos not in [239,253] ]

temSeq="MSIQHFRVALIPFFAAFCLPVFAHPETLVK-VKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLL"
temSeq+="CGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELT"
temSeq+="AFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLR"
temSeq+="SALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW-----"

pos2wt = dict(zip(posList,temSeq))
subName = "TEM"
tip2pheno = {}
phenoClades = []
phenotypeList = ["2b","2be","2ber"]




def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

def parseArgs():
        parser = argparse.ArgumentParser()
        parser.add_argument('-c', dest = "clinFile", required=True,
                            help = "Path of fasta files  wih the MSA, "+\
                            "separated by ;.")
        parser.add_argument('-l', dest = "labFile", required=True, default="None",
                            help = "Path of fasta files  wih the MSA, "+\
                            "separated by ;.")
        parser.add_argument('-p', dest = "leafPheno", required=True,
                            help = "Path of file with names of leafnodes"+\
                            " and their phenotype.")
        parser.add_argument('-o', dest = "outPath", required=True,
                            help = "Path for undirected alignment pairs summary.")
        args = parser.parse_args()
        return args

def get_undirPairs(fastaFile, labFile,tipPhenoFile, outPath):
    phenoDict = {'1':'2br','2':'2b','3':'2ber','4':'2be'}
    with open(tipPhenoFile,'rb') as pf:
        lineRecs = [line.strip().split() \
                    for line in pf.read().replace('"','').split('\n')]
        for rec in lineRecs:
            if len(rec)>1:
                if rec[1].strip() not in phenoDict: continue
		tip2pheno[rec[0]]=phenoDict[rec[1].strip()]

    #stderr_write([tip2pheno])

    # First the network from clinical sequences
    undirPairSetClin = UndirPairSet([fastaFile], posList, pos2wt,
                                  tip2pheno)
    stderr_write(["\t\tNumber of undir (mut) pairs", len(undirPairSetClin.undirPairs)])
    for undirPair in undirPairSetClin.undirPairs:
        undirPair.get_FisherPvals()


    outUndirFile = outPath.rstrip('/')+'/TEM_undirPairsPhenoClin.txt'
    with open(outUndirFile, 'wb') as eF:
        eF.write(str(undirPairSetClin))
        
    outUndirFile = outPath.rstrip('/')+'/TEM_undirPairsPhenoClin_Ji.txt'
    with open(outUndirFile, 'wb') as eF:
        eF.write(undirPairSetClin.str_jaccard())

    #Then the network that also includes laboratory evolved sequences
    undirPairSet = UndirPairSet([fastaFile,labFile], posList, pos2wt,
                                  tip2pheno)
    stderr_write(["\t\tNumber of undir (mut) pairs", len(undirPairSet.undirPairs)])
    for undirPair in undirPairSet.undirPairs:
        undirPair.get_FisherPvals()


    outUndirFile = outPath.rstrip('/')+'/TEM_undirPairsPhenoClinDirevo.txt'
    
    stderr_write(["\tWrite undirPairs to file",outUndirFile.split('/')[-1],
                  datetime.now()])
    with open(outUndirFile, 'wb') as eF:
        eF.write(str(undirPairSet))

    outUndirFile = outPath.rstrip('/')+'/TEM_undirPairsPhenoClinDirevo_Ji.txt'
    with open(outUndirFile, 'wb') as eF:
        eF.write(undirPairSet.str_jaccard())
        

def main():
    args = parseArgs()
    get_undirPairs(args.clinFile, args.labFile, args.leafPheno, args.outPath) 
    
if __name__ == "__main__":
    main()
