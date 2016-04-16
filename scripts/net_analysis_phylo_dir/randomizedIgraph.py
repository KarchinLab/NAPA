from igraph import *
import argparse
import os.path
from collections import defaultdict, Counter
from itertools import islice
import math
import re
import copy
import random

getInt = lambda s: int(re.sub('\D','',s))
stderr_write = lambda myList: sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

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
    parser.add_argument('-e', dest = "numEdgeSwaps", type=int, default=10000, 
                            help = 'How many times to swap the edges.')
    parser.add_argument('-n', dest = "numRandNets", type=int, default=500,
                            help = "How many random networks to generate.")
    args = parser.parse_args()
    return args

class MutGraph(object):
    def __init__(self,g=None,dirNetFile=None,fileIsNet=None,fromGraph=None):
        if fromGraph == True and g!=None:
            self.g = copy.deepcopy(g)

        elif dirNetFile!=None and os.path.isfile(dirNetFile):
            if fileIsNet:
                self.g=Graph.Read_Ncol(dirNetFile,directed=True,names=True)
            else:
                self.mutPairs_to_dirNet(dirNetFile)
            vtxToDelete = [v.index for v in self.g.vs if "_" in v["name"]]
            self.g.delete_vertices(vtxToDelete)
        else:
            stderr_write([ "Generating empty graph!"])
            g.Graph(directed=True) 
                        
        self.vertexNames = [v['name'] for v in self.g.vs]
        self.edgeNames = [self.g.vs[ei[0]]['name']+'\t'+self.g.vs[ei[1]]['name'] \
                          for ei in self.g.get_edgelist()]
        self.dirNetFile = dirNetFile
        self.update_degrees()
        
    def deepcopy(self):
        return MutGraph(g=copy.deepcopy(self.g), fromGraph=True)
        
    def write_edgeList(self, outFile):
        with open(outFile,'wb') as of:
            for ei,e in enumerate(self.g.get_edgelist()): 
                of.write('\t'.join([self.g.vs[e[0]]["name"],self.g.vs[e[1]]["name"],
                                    str(int(self.g.es[ei]["weight"]))])+'\n')
    def write_edgeList_PICKL(self,outFile):
        self.g.write_pickle(fname=outFile)

    def mutPairs_to_dirNet(self,dirNetFile):
        tupleList = []
        with open(dirNetFile, 'rb') as dnf:
            for line in dnf:
                if line.startswith('mut'): continue
                recs = line.split('\t')
                [mut1,mut2,incr12,incr21] = [recs[0],recs[1],recs[2],recs[4]]
                if int(incr12)>=1: tupleList.append((mut1,mut2,float(incr12))) 
            if int(incr21)>=1: tupleList.append((mut2,mut1,float(incr21))) 
            
        self.g = Graph.TupleList(edges=tupleList, directed=True,edge_attrs=("weight"))
        outNetFile = dirNetFile.replace('.txt','_net.txt')
        with open(outNetFile, 'wb') as onf:
            for t in tupleList:
                onf.write('\t'.join(t[0:-1])+'\t'+str(t[-1])+'\n')

    def update_degrees(self):
        self.vertexWeightedDegreesOut = self.g.strength(vertices=self.g.vs,mode=OUT,
                                                        loops=False,weights=self.g.es['weight'])  
        self.vertexWeightedDegreesIn =  self.g.strength(vertices=self.g.vs,mode=IN,
                                                        loops=False,weights=self.g.es['weight'])
        self.outDegreeNodes = [v for i,v in enumerate(self.g.vs) \
                               if self.vertexWeightedDegreesOut[i] > 0]
        self.node2outNeighb = {}
        for v in self.g.vs:
            self.node2outNeighb[v] = set([self.g.vs[ni] for ni in self.g.neighbors(v,mode=OUT)])

    def dirSwap(self):
        [v1,v2] = random.sample(self.outDegreeNodes,2)
        v1v2neighb = self.node2outNeighb[v1] & self.node2outNeighb[v2]
        v1neighb = self.node2outNeighb[v1] - set([v2])
        v2neighb = self.node2outNeighb[v2] - set([v1])
        if len(v1neighb) and len(v2neighb):
            v1n, v2n = random.choice(list(v1neighb)),random.choice(list(v2neighb))
            
            v1nEdgeWeight=self.g.es[self.g.get_eid(v1.index,v1n.index)]["weight"]
            v2nEdgeWeight=self.g.es[self.g.get_eid(v2.index,v2n.index)]["weight"]
            minWeight = min([v1nEdgeWeight,v2nEdgeWeight])
            subWeight = random.randint(1,minWeight) if minWeight>1 else 1
            
            if subWeight < self.g.es[self.g.get_eid(v1.index,v1n.index)]["weight"]:
                self.g.es[self.g.get_eid(v1.index,v1n.index)]["weight"] -= subWeight
            else:
                self.g.delete_edges(self.g.get_eid(v1.index,v1n.index))

            if subWeight < self.g.es[self.g.get_eid(v2.index,v2n.index)]["weight"]:
                self.g.es[self.g.get_eid(v2.index,v2n.index)]["weight"] -= subWeight
            else:
                self.g.delete_edges(self.g.get_eid(v2.index,v2n.index))
                
            if self.g.get_eid(v1.index,v2n.index,error=False) > -1:
                self.g.es[self.g.get_eid(v1.index,v2n.index)]['weight']+=subWeight 
            else:
                self.g.add_edge(v1, v2n, weight = subWeight) # actual swap

            if self.g.get_eid(v2.index,v1n.index,error=False) > -1: 
                self.g.es[self.g.get_eid(v2.index,v1n.index)]['weight']+=subWeight
            else:
                self.g.add_edge(v2, v1n, weight = subWeight) # actual swap
            
            self.update_degrees()

    def randomizeNet(self,nSwaps):
        for n in range(nSwaps):
            self.dirSwap()
                                
def make_randomSet(mg, setSize,nSwaps, dirPairFileName):
    for i in range(setSize):
        rmg = mg.deepcopy()
        rmg.randomizeNet(nSwaps)
        if mg.vertexWeightedDegreesIn != rmg.vertexWeightedDegreesIn:
            stderr_write([ "original IN degrees\n\t", mg.vertexWeightedDegreesIn])
            stderr_write(["random", i+1, "IN degrees\n\t", rmg.vertexWeightedDegreesIn])
            stderr_write(["original OUT degrees\n\t", mg.vertexWeightedDegreesOut])
            stderr_write(["random", i+1, "OUT degrees\n\t", rmg.vertexWeightedDegreesOut])
        outF = dirPairFileName.replace(".txt",".rand"+str(i+1)+".swp"+str(nSwaps)+".pyPickle")
        outF = '/'.join(outF.split('/')[0:-1])+'/randomNets/'+outF.split('/')[-1]
        #rmg.write_edgeList(outF.replace('.pyPickle','txt'))
        rmg.write_edgeList_PICKL(outF)
    
def main():
    # randomizedIgraph.py -d $netFileName -e 1 -n 1 -f True
    args = parseArgs()
    stderr_write([args])
    mg = MutGraph(dirNetFile=args.dirPairFileName, fileIsNet=args.fileIsNet)
    make_randomSet(mg, args.numRandNets, args.numEdgeSwaps, 
                   args.dirPairFileName)

if __name__ == "__main__": main()
