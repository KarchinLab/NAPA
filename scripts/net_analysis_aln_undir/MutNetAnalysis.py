from igraph import *
import argparse
import os.path
import pandas as pd
from collections import defaultdict, Counter
from itertools import islice
import math
import re
from scipy.misc import comb

getInt = lambda s: int(re.sub('\D','',s))

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
        parser.add_argument('-d', dest = "epiPairFileName", type=str, 
                            required=True,
                            help = 'Path to where weighted mutation'+\
                            ' pairs are stored.')
        parser.add_argument('-f', dest = "fileIsNet", type=int, default=0, 
                            help = 'Whether or not the file is a ready-to-read'+\
                            ' tab-delimited network')
        parser.add_argument('-g', dest = "fromGraph", type=int, default=0,
                            help = "Read a (pickled) igraph?")
        parser.add_argument('-l', dest = "pathLength", type=int, default=3, 
                            help = 'Length of central paths to consider, in number'+\
                            ' of nodes.')
        args = parser.parse_args()
        return args

class MutSubPath(object):
    def __init__(self, subPathStr, source, target, stPathNum, stPathNumSP):
        self.subPathStr = subPathStr
        self.length = len(subPathStr.split('_'))
        self.sourceTarget2frac = {'_'.join([source,target]):float(stPathNumSP)/stPathNum}
    
    def __repr__(self):
        return '\t'.join([self.subPathStr,str(self.centrality)])

    def add_stPair(self, source, target, stPathNum, stPathNumSP):
        self.sourceTarget2frac['_'.join([source,target])] = float(stPathNumSP)/stPathNum

    def calc_betweenness(self, n=3):
        scale = float(comb(n-2,self.length))/(comb(n-self.length,2)) 
        self.centrality = sum(self.sourceTarget2frac.values()) * scale

class MutGraph(object):
    def __init__(self,g=None,epiNetFile=None,fileIsNet=None,fromGraph=None):
        if fromGraph == True and epiNetFile==None:
            if g!=None:
                self.g = copy.deepcopy(g)
            else:
                stderr_write(["No Igraph supplied, generating empty!"])
                self.g = Graph(directed=False)
        elif fromGraph == True and epiNetFile !=None:
            if os.path.isfile(epiNetFile):
                self.g = Graph.Read_Pickle(epiNetFile)
            else:
                stderr_write(["Pickled igraph d.n.e, generating empty!"]) 
                self.g = Graph(directed=False)
        elif fromGraph == False and epiNetFile!=None and os.path.isfile(epiNetFile):
            if fileIsNet:
                #self.g=Graph.Read_Ncol(epiNetFile,directed=False,names=True)
                self.readNet(epiNetFile)
            else:
                self.mutPairs_to_epiNet(epiNetFile)
            vtxToDelete = [v.index for v in self.g.vs if "_" in v["name"]]
            self.g.delete_vertices(vtxToDelete)
        else:
            stderr_write([ "No valid inputs. Generating empty graph!"])
            self.g = Graph(directed=True) 
              
        self.g.es['invWeight'] = [5040.0/w for w in self.g.es['weight']]
        self.vertexNames = [v['name'] for v in self.g.vs]
        self.edgeNames = [self.g.vs[ei[0]]['name']+'\t'+self.g.vs[ei[1]]['name'] 
                          for ei in self.g.get_edgelist()]
            
        self.modularity = None
        self.epiNetFile = epiNetFile

    def readNet(self, dirNetFile):
        tupleList=[]
        with open(dirNetFile,'rb') as dnf:
            for line in dnf:
                if line.startswith('mut') or line.startswith("_"): 
                    continue
                recs = line.split('\t')
                if len(recs) < 3: continue
                [mut1,mut2,weight] = [recs[0],recs[1],recs[2]]
                tupleList.append((mut1,mut2,float(weight)))
        self.g = Graph.TupleList(edges=tupleList, directed=False,
                                 edge_attrs=("weight"))

    def mutPairs_to_dirNet(self):
        lineList = []
        with open(self.dirNetFile, 'rb') as dnf:
            for line in dnf:
                if line.lower().strip('_').startswith('mut'): continue
                recs = line.strip().split('\t')
                stderr_write(recs)
                [mut1,mut2,pval] = [recs[0],recs[1], float(recs[-1])]
                assoc = 1.0-pval
                if assoc >0:
                    lineList.append('\t'.join([mut1,mut2,str(assoc)]))
            outNetFile = self.dirNetFile.replace('.txt','_nx.txt')
            stderr_write([outNetFile])
            with open(outNetFile, 'wb') as onf:
                onf.write('\n'.join(lineList))
        return lineList

    def mutPairs_to_epiNet(self,epiNetFile):
        tupleList = []
        with open(epiNetFile, 'rb') as dnf:
            for line in dnf:
                if line.startswith('mut'): continue
                recs = line.split('\t')
                [mut1,mut2,incr12,incr21] = [recs[0],recs[1],recs[2],
                                             recs[4]]
                if int(incr12)>=1: tupleList.append((mut1,mut2,float(incr12))) 
                if int(incr21)>=1: tupleList.append((mut2,mut1,float(incr21))) 
            
        self.g = Graph.TupleList(edges=tupleList, directed=False,
                                 edge_attrs=("weight"))
        outNetFile = epiNetFile.replace('.txt','_net.txt')
        with open(outNetFile, 'wb') as onf:
            for t in tupleList:
                onf.write('\t'.join(t[0:-1])+'\t'+str(t[-1])+'\n')

    def get_node_properties(self):
        self.get_node_degree()
        self.get_node_closeness()
        self.get_node_betweenness()
        self.get_node_multilevel_communities()

    def get_edge_properties(self):
        self.get_edge_betweenness()
        self.get_edge_degree()
        
    def get_graph_properties(self):
        if not self.modularity:
            self.get_node_multilevel_communities()
            self.degAssort = self.g.assortativity_degree(directed=False)

    def get_node_degree(self):
        self.vertexUnweightedDegreesOut = self.g.degree(vertices=self.g.vs,mode=OUT,
                                                        loops=False)  
        self.vertexUnweightedDegreesIn =  self.g.degree(vertices=self.g.vs,mode=IN,
                                                        loops=False)  
        self.vertexUnweightedTotalDegree =self.g.degree(vertices=self.g.vs, mode=ALL,
                                                        loops=False)
        self.vertexWeightedDegreesOut = self.g.strength(vertices=self.g.vs,mode=OUT,
                                                        loops=False,
                                                        weights=self.g.es['weight'])  
        self.vertexWeightedDegreesIn =  self.g.strength(vertices=self.g.vs,mode=IN,
                                                        loops=False,
                                                        weights=self.g.es['weight'])  
        self.vertexWeightedTotalDegree = self.g.strength(vertices=self.g.vs, mode=ALL,
                                                         loops=False,
                                                         weights=self.g.es['weight'])
    def get_node_closeness(self):
        self.vertexOutCloseness = self.g.closeness(mode="OUT",weights='invWeight')
        self.vertexInCloseness  = self.g.closeness(mode="IN", weights='invWeight')
        self.vertexTotCloseness = self.g.closeness(mode="ALL",weights="invWeight")
        
    def get_node_betweenness(self):
        self.nodeUndirBetweeness = self.g.betweenness(directed=False,
                                                      weights=self.g.es['invWeight'])
    def get_node_multilevel_communities(self):
        '''VD Blondel, J-L Guillaume, R Lambiotte and E Lefebvre: Fast 
        unfolding of community hierarchies in large networks, J Stat Mech 
        P10008 (2008), http://arxiv.org/abs/0803.0476'''
        self.g_undir = self.g.copy()
        self.communities = self.g_undir.community_multilevel(\
                            weights=self.g_undir.es['weight'],return_levels=False)
        self.nodeCommunities = self.communities.membership
                
        self.modularity = self.g_undir.modularity(\
                                        membership=self.communities.membership,
                                            weights=self.g_undir.es['weight'])

    def get_edge_betweenness(self):
        self.edgeBetweenness = self.g.edge_betweenness(directed=True,
                                                               weights='invWeight')
    def get_edge_degree(self):
        [sourceIndices, targetIndices] = zip(*self.g.get_edgelist()) # bad for large graphs?
        self.edgeWeightedDegreesOut = [self.vertexWeightedDegreesOut[tI] \
                                       for tI in targetIndices]
        self.edgeUnweightedDegreesOut = [self.vertexUnweightedDegreesOut[tI] \
                                         for tI in targetIndices]
        self.edgeWeightedDegreesIn = [self.vertexWeightedDegreesIn[sI] \
                                      for sI in sourceIndices]
        self.edgeUnweightedDegreesIn = [self.vertexUnweightedDegreesIn[sI] \
                                        for sI in sourceIndices]
                        
    def print_node_properties(self):
        outFile = self.epiNetFile.replace('.txt','_net_nodeProps.txt')
        with open(outFile, 'wb') as of:
            of.write('\t'.join(["Node","uwOutDeg","uwInDeg","uwTotDegree",
                                "wOutDegree","wInDeg","wTotDeg",
                                "wOutCloseness","wInCloseness","wTotCloseness",
                                "wDirBetweeness", "wUndirBetweenness","community"])+'\n')
        for ni in range(len(self.g.vs)):
            degStr = '\t'.join(str(d) for d in [self.vertexUnweightedDegreesOut[ni],
                                                self.vertexUnweightedDegreesIn[ni],
                                                self.vertexUnweightedTotalDegree[ni],
                                                self.vertexWeightedDegreesOut[ni],
                                                self.vertexWeightedDegreesIn[ni],
                                                self.vertexWeightedTotalDegree[ni]])
                                
            centStr = '\t'.join(str(c) for c in [self.vertexOutCloseness[ni],
                                                 self.vertexInCloseness[ni],
                                                 self.vertexTotCloseness[ni],
                                                 self.nodeDirBetweenness[ni],
                                                 self.nodeUndirBetweeness[ni]])
                                
            communityStr = str(self.nodeCommunities[ni])
                                
            of.write('\t'.join([self.vertexNames[ni],degStr,centStr,communityStr])+'\n')

    def print_node_communities(self):
        outFile = self.epiNetFile.replace('.txt','_net_nodeCommunities.txt')
        with open(outFile, 'wb') as of:
            of.write('\t'.join(["Node","community"])+'\n')
            for ni in range(len(self.g.vs)):                                
                communityStr = str(self.nodeCommunities[ni])
                of.write('\t'.join([self.vertexNames[ni],communityStr])+'\n')

    def print_edge_properties(self):
        outFile = self.epiNetFile.replace('.txt','_net_edgeProps.txt')
        with open(outFile,'w') as of:
            of.write('\t'.join(["Node1","Node2","edgeWeight","invWeight",
                                "uwDegOut","uwDegIn","uwTotDeg","wDegOut","wDegIn","wTotDeg",
                                "wDirBetweenness"])+'\n')
            for ei in range(len(self.g.es)):
                degStr='\t'.join(str(d) for d in [self.edgeUnweightedDegreesOut[ei],
                                                  self.edgeUnweightedDegreesIn[ei],
                                                  self.edgeUnweightedDegreesOut[ei]+\
                                                  self.edgeUnweightedDegreesIn[ei],
                                                  self.edgeWeightedDegreesOut[ei],
                                                  self.edgeWeightedDegreesIn[ei],
                                                  self.edgeWeightedDegreesOut[ei]+\
                                                  self.edgeWeightedDegreesIn[ei]])
            of.write('\t'.join([self.edgeNames[ei],str(self.g.es['weight'][ei]),
                                str(self.g.es['invWeight'][ei]),degStr,
                                str(self.edgeBetweenness[ei])])+'\n')
        
    def print_graph_properties(self):
        outFile = self.epiNetFile.replace('.txt','_netProps.txt')
        with open(outFile,'w') as of:
            of.write('\t'.join(["Graph Modularity:",str(self.modularity)])+'\n')
            of.write('\t'.join(["Graph Vertex Degree Assortativity:",
                                    str(self.degAssort)])+'\n')

def epi_path_betweenness(mg,pathLen):
    subPath2obj = {}

    for sourceI, source in enumerate(mg.g.vs):
        if sourceI == len(mg.g.vs)-1: break
        reachableTargetI = list(set(range(sourceI+1, len(mg.g.vs))) &\
                                set(mg.g.neighborhood(vertices=source,order=len(mg.g.vs), 
                                                      mode='all')))
        for targetI in reachableTargetI:
            if targetI == sourceI: continue
            target = mg.g.vs[targetI]
            #print source['name'], target['name']
            sourceTargetSPs = mg.g.get_all_shortest_paths(v=source, to=target,
                                                          weights=mg.g.es['invWeight'])
            #print source['name'], target['name'], len(sourceTargetSPs), sourceTargetSPs
            
            stSubPaths = []
            for shortPath in sourceTargetSPs:
                #slideRight = [w for w in window(shortPath[1:-1],n=pathLen)]
                #slideLeft = [w for w in window(shortPath[-2:0:-1],n=pathLen)]
                #stSubPaths += list(set(slideRight) | set(slideLeft))
                stSubPaths += [w for w in window(shortPath[1:-1],n=pathLen)]

            subPath2count = Counter(stSubPaths)
            for subPath in subPath2count:
                subPathStrList = [str(mg.g.vs[v]['name']) for v in subPath]
                subPathPosList = [getInt(m) for m in subPathStrList]
                if len(set(subPathPosList)) < len(subPathPosList): continue
                
                subPathStr = '_'.join(subPathStrList)
                if subPathStr in subPath2obj:
                    subPath2obj[subPathStr].add_stPair(source['name'],target['name'],
                                                       len(sourceTargetSPs),
                                                       subPath2count[subPath])
                else:
                    subPath2obj[subPathStr]=MutSubPath(subPathStr,
                                                       source['name'],target['name'],
                                                       len(sourceTargetSPs),
                                                       subPath2count[subPath])
            
    for subPath in subPath2obj:
        subPath2obj[subPath].calc_betweenness(n=len(mg.g.vs))

    return subPath2obj

def print_epi_betw_path(subPath2obj,outFile):
    with open(outFile,'w') as of:
        for subPath in sorted(subPath2obj.keys()): 
            of.write(str(subPath2obj[subPath])+'\n')
                                
def main():
    # MutNetAnalysis.py -d $netFileName -l 3
    args = parseArgs()
    mg = MutGraph(epiNetFile=args.epiPairFileName,fileIsNet=bool(args.fileIsNet),
                  fromGraph=bool(args.fromGraph))

    #mg.get_node_properties()
    #mg.get_edge_properties()
    mg.get_graph_properties()
    #mg.print_node_properties()
    #mg.print_edge_properties()
    mg.print_node_communities()

    pathLength = args.pathLength
    betweennessFile = args.epiPairFileName.replace('txt','spBetw.'+str(pathLength)+'.txt')
    p2obj = epi_path_betweenness(mg,pathLength)
    print_epi_betw_path(p2obj,betweennessFile)

if __name__ == "__main__": main()
