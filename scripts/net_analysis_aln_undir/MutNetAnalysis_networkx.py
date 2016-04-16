import networkx as nx
import argparse
import sys
import os.path
from collections import defaultdict
from itertools import islice
import math
import re
from kpath import *
from kpath import kpath_probPath_centrality
from scipy.misc import comb

getInt = lambda s: int(re.sub('\D','',s))

def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def parseArgs():
    parser = argparse.ArgumentParser( \
    description = 'Characterize mutation graph/vertices/edges')
    
    parser.add_argument('-d', dest = "dirPairFileName", type=str, 
                        required=True,
                        help = 'Path to where weighted directed mutation'+\
                        ' pairs are stored.')
    parser.add_argument('-f', dest = "fileIsNet", type=bool, default=False, 
                        help = 'Whether or not the file is a ready-to-read'+\
                        ' tab-delimited network')
    parser.add_argument('-l', dest = "pathLength", type=int, default=3, 
                        help = 'Length of central paths to consider, in number'+\
                        ' of nodes.')
    args = parser.parse_args()
    return args

class MutGraphNX(object):
    def __init__(self,g=None,dirNetFile=None,fileIsNet=None):
        self.dirNetFile = dirNetFile
        print self.dirNetFile

        if dirNetFile!=None and os.path.isfile(dirNetFile):
            stderr_write([fileIsNet])
            if fileIsNet:
                print dirNetFile
                with open(dirNetFile, 'rb') as dnf:
                    lines = dnf.readlines()[1:]
            else:
                lines = self.mutPairs_to_dirNet()       
            self.g=nx.parse_edgelist(lines, data=(('weight',float),), 
                                     create_using=nx.Graph())
            for u,v,d in self.g.edges(data=True):
                d['invWeight'] = 1./d['weight']
            print self.g
        
        elif g!=None:
            self.g = g
        else:
            print "WARNING: Network File not found:\n", dirNetFile
            print "Generating empty graph!"
            self.g = nx.DiGraph()
        

    def mutPairs_to_dirNet(self):
        lineList = []
        with open(self.dirNetFile, 'rb') as dnf:
            for line in dnf:
                if line.lower().startswith('_mut'): continue
                recs = line.strip().split('\t')
                #stderr_write(recs)
                [mut1,mut2,pval] = [recs[0],recs[1], float(recs[-1])]
                assoc = 1.0-pval
                if assoc >0:
                    lineList.append('\t'.join([mut1,mut2,str(assoc)]))
            outNetFile = self.dirNetFile.replace('.txt','_nx.txt')
            stderr_write([outNetFile])
            with open(outNetFile, 'wb') as onf:
                onf.write('\n'.join(lineList))
        return lineList

    def get_path_cent(self,pathLength):
        kappa = None
        #kappa = (pathLength+2)*2
        self.kPath_pathCent = kpath_path_centrality(self.g, k=kappa,
                                                        seed=12345,alpha=0.1, 
                                                        weight='invWeight',
                                                        pathLen=pathLength)

        self.kPath_pathProb = kpath_probPath_centrality(self.g, k=kappa,
                                                        seed=12345,alpha=0.1, 
                                                        weight='invWeight',
                                                        pathLen=pathLength)
            

    def print_path_cent(self,pathLength):
        path_outFile = self.dirNetFile.replace('.txt','_nx_pathBetw_'+str(pathLength)+'.txt')
        with open(path_outFile, 'wb') as pof:
            pof.write('\t'.join(["Path","kPathCentrality"])+'\n')
            for path in sorted(self.kPath_pathCent.keys()):
                pof.write(path+'\t'+str(self.kPath_pathCent[path])+'\n')
        
        pathProb_outFile = self.dirNetFile.replace('.txt',
                                            '_nx_pathProb_'+str(pathLength)+'.txt')
        with open(pathProb_outFile, 'wb') as pof:
            pof.write('\t'.join(["Path","kPathProbCent"])+'\n')
            for path in self.kPath_pathCent:
                pof.write(path+'\t'+str(self.kPath_pathProb[path])+'\n')
        

    def print_node_properties(self,pathLength):
        node_outFile = self.dirNetFile.replace('.txt','_nx_node_props.txt')
        path_outFile = self.dirNetFile.replace('.txt',
                                            '_nx_path_'+str(pathLength)+'.txt')
        pathProb_outFile = self.dirNetFile.replace('.txt',
                                            '_nx_pathProb_'+str(pathLength)+'.txt')
            
        with open(node_outFile, 'wb') as nof:
            nof.write('\t'.join(["Node","kPathCentrality"])+'\n')
            for node in self.kPath_nodeCent:
                nof.write(node+'\t'+str(self.kPath_nodeCent[node])+'\n')

        with open(path_outFile, 'wb') as pof:
            pof.write('\t'.join(["Path","kPathCentrality"])+'\n')
            for path in self.kPath_pathCent:
                pof.write(path+'\t'+str(self.kPath_pathCent[path])+'\n')

        with open(pathProb_outFile, 'wb') as pof:
            pof.write('\t'.join(["Path","kPathProbCent"])+'\n')
            for path in self.kPath_pathCent:
                pof.write(path+'\t'+str(self.kPath_pathProb[path])+'\n')

        if not nx.is_directed(self.g):
            print "\nNode current flow closeness:\n", self.nodeCurrentFlowCloseness
            print "\nNode current flow betweenness:\n", self.nodeCurrentFlowBetweenness

            
                                
                
def main():
    #python MutNetAnalysis_networkx.py -d $netFileName -f True -l 3

    args = parseArgs()
    mg = MutGraphNX(dirNetFile=args.dirPairFileName,fileIsNet=False)
    stderr_write(["Obtaining kpath path centralities"])

    mg.get_path_cent(args.pathLength)
    mg.print_path_cent(args.pathLength)


if __name__ == "__main__": main()
