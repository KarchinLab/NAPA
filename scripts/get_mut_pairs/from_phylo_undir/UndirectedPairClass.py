#! /usr/bin/env python 
import sys, re
from collections import defaultdict
from datetime import datetime
from ete2 import *
from scipy import stats
from itertools import *
#from rpy2.robjects.packages import importr
#from rpy2.robjects.vectors import FloatVector
#rstats = importr('stats')

def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

allowed_positions = set(range(300)) - set(range(23+1))
flatten = lambda x: [subitem for item in x for subitem in item]
listPairs = lambda x: [(x[i1],x[i2]) for i1 in range(len(x)-1) \
                       for i2 in range(i1+1,len(x))]
getInt = lambda s: int(re.sub('\D','',s))
sortByDigits = lambda x: sorted(x, key = lambda s: getInt(s))
sortListbyOtherList = lambda x,y: [x for (y,x) in sorted(zip(Y,X),
                                                key=lambda pair: pair[0])]
# allowed transitions between phenotypes
phenTransitions = {'2b':['2be','2ber'],'2be':['2be'],'2br':['2be','2ber'], 
                   '2ber':['2be','2ber']} # don't lose e, don;t gain just r


class PhyloEdge(object):
    def __init__(self,nodePair=None,phenoDict={},parent=None,
                 children=[],posList=[],pos2wt={}):
        self.nodePair = nodePair
        self.phenoDict = phenoDict.copy()
        #stderr_write([[n.name for n in nodePair if n!=None],'phenoDict', phenoDict])
        nodePair0name = "ANC" if not nodePair[0] else nodePair[0].name
        self.name = nodePair0name +'__'+nodePair[1].name
        
        self.pos2wt = pos2wt
        self.posList = posList
        self.get_mutations()
        self.positions = [getInt(m) for m in self.mutations]
        self.parent = parent # parent edge = another node pair
        self.children = children
        if not(len(self.children)): self.get_node_children()
        self.precursors = [] 
        self.followers = []
        self.orderedPrecursors = []

    def build_EdgeTree(self):
        self.children = [child for child in self.children \
                         if child.parent.name==self.name]

        stderr_write(["\tBuilding edge followers", datetime.now()])
        self.build_edge_followers() 
        stderr_write(["\tBuilding edge precursors", datetime.now()])
        self.build_edge_precursors()
        
    def __iter__(self):
        for pEdge in chain(*imap(iter, self.children)):
            yield pEdge
        yield self
            
    def __repr__(self):
        toPrint = [self.name, 'mutations: '+';'.join(self.mutations), 
                   'children: '+';'.join([c.name for c in self.children]), 
                   'followers: '+';'.join([f.name for f in self.followers]),
                   'precursors: '+';'.join([p.name for p in self.precursors])]
        if len(self.orderedPrecursors):
            toPrint.append('orderedPrec: '+';'.join([op.name \
                                            for op in self.orderedPrecursors \
                                            if op])) 
        return '\n\t'.join(toPrint)

    def get_mutations(self):
        ''' Get mutations connecting two phylogeny nodes (sequences)
        nodePair: an instance of the PhyloEdge class
        seqId2seq: dictionary of sequence ids vs their AA sequences
        '''
        standardAA="ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids
        nodeSeq = self.nodePair[1].sequence # phyloTree linked to MSA
        ancSeq = self.nodePair[0].sequence if self.nodePair[0] else nodeSeq
        self.mutations =  [ancSeq[i]+str(self.posList[i])+nodeSeq[i] \
                for i in range(len(nodeSeq)) if ancSeq[i]!=nodeSeq[i] \
                and ancSeq[i]==self.pos2wt[self.posList[i]] \
                           and nodeSeq[i] in standardAA \
                           and self.posList[i] > 23]

    def get_mutPositions(self,posList,pos2wt):
        standardAA="ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids  
        ancSeq = self.nodePair[0].sequence # phyloTree must be linked to alignmnent
        nodeSeq = self.nodePair[1].sequence
        return [posList[i] for i in range(len(nodeSeq)) \
                if ancSeq[i]!=nodeSeq[i] and \
                ancSeq[i]==pos2wt[posList[i]] and \
                nodeSeq[i] in standardAA and self.posList[i] > 23]
    
    def isLeafEdge(self):
        return self.nodePair[1].is_leaf()

    def check_dist(self, other, distThresh):
        return self.nodePair[0].get_distance(other.nodePair[1]) <= distThresh

    def check_pheno(self):
        phen0 = self.phenoDict[self.nodePair[0].name] \
                 if self.nodePair[0] and self.nodePair[0].name in self.phenoDict \
                    else 'unknown'
        phen1 = self.phenoDict[self.nodePair[1].name] \
                 if self.nodePair[1].name in self.phenoDict else 'unknown'
        #stderr_write([self.nodePair[0].name, self.nodePair[1].name, phen0,phen1])
        phenoCheck = len(self.mutations)>0 and  phen0 in phenTransitions \
            and phen1 in phenTransitions[phen0]
        return phenoCheck

    def check_pos(self):
        if not len(self.positions): return False
        return set(self.positions) <= set(allowed_positions)

    def add_child(self,edgeNode):
        if edgeNode.name not in [c.name for c in self.children]:
            self.children.append(edgeNode)
            edgeNode.parent = self

    def get_node_children(self):
        self.children = []
        phyloNode = self.nodePair[1]
        nodeChildren = phyloNode.get_children()
        for childNode in nodeChildren:
            child = PhyloEdge(nodePair=[phyloNode,childNode],
                              phenoDict=self.phenoDict,
                              parent=self,posList=self.posList,
                              pos2wt=self.pos2wt)
            self.add_child(child)

    def build_edge_followers(self):
        for child in self.children:
            child.build_edge_followers()
            self.followers += [child] + child.followers
        self.followers = list(set(self.followers))

    def build_edge_precursors(self):
        for child in self.children:
            child.precursors = list(set(self.precursors + [self]))
            child.build_edge_precursors()
   
    def get_common_precursors(self, other):
        commonPrec = []
        if self.parent==None or other.parent==None:
            commonPrec = [] # case of one of the edges being root
        elif self == other:
            commonPrec = [self, self.parent]
        elif self in other.precursors:
            commonPrec = [self.parent, self]
        elif other in self.precursors:
            commonPrec = [other.parent, other]
        else: # they are in different lineages: look for MRCA
            commonPrec = list(set(self.precursors) & set(other.precursors))
            if not len(commonPrec):
                print "No common precursors for", self.name, other.name
            elif len(commonPrec) > 1:
                commonPrec = [self.orderedPrecursors[len(commonPrec)-1]]
        commonPrec = list(set([p for p in commonPrec if p.check_pheno()]))
        return commonPrec

    def get_intermediatePrec(self, precursorList):
        intermediates = [self] if self.check_pheno else []

        if not len(precursorList):
            intermediates += [p for p in self.orderedPrecursors \
                              if p.check_pheno()]
            return intermediates

        firstPrec = self.orderedPrecursors.index(\
                                        next((p for p in self.orderedPrecursors \
                                              if p in precursorList), 
                                             self.orderedPrecursors[-1]))
        return intermediates + [im for im in self.orderedPrecursors[:firstPrec] \
                                if im.check_pheno()]

    def assign_ordered_precursors(self):
        if not len(self.precursors): 
            self.orderedPrecursors = []
            return
        nodePairDict = {}
        for p in self.precursors:
            pName = p.nodePair[0].name if p.nodePair[0] else 'ANC'
            npName = p.nodePair[1].name if p.nodePair[1] else 'ANC'
            nodePairDict[pName] = npName

        mostAncNodeName = list(set(nodePairDict.keys()) - \
                           set(nodePairDict.values()))
        if not(len(mostAncNodeName)): print "No ancestor found!", self.name
        if len(mostAncNodeName) > 1: print "Found multiple ancestors!"
        
        currAncestorName = mostAncNodeName[0]
        maxLen = len(set(nodePairDict.keys()+nodePairDict.values()))

        for i in range(maxLen-1):
            nextAncestorName = nodePairDict[currAncestorName]
            nextPrecursor = [p for p in self.precursors if p.name == \
                             currAncestorName+'__'+nextAncestorName]
            self.orderedPrecursors.append(nextPrecursor[0])
            currAncestorName = nextAncestorName

        if self.orderedPrecursors[-1] != self.parent or \
           len(self.orderedPrecursors)!= len(self.precursors):
            print "Precursors not properly ordered/found:", self.name
            print "\tnodePairDict", nodePairDict
            print "\tordered", [p.name for p in self.orderedPrecursors if p!=None]
            print "\tunordered", [p.name for p in self.precursors]

    def get_newick(self,nwkStr='',rootParent=None):
        if self.parent == rootParent:
            if not len(self.children): 
                return nwkStr+str(self.name)
            childNames = [c.name for c in self.children]
            childNames.sort()
            childrenStr = ','.join([str(c) for c in childNames])
            if str(self.name) not in re.findall(r"[\w']+",nwkStr):
                nwkStr+='('+childrenStr+')' + str(self.name)
                for child in self.children:
                    nwkStr = child.get_newick(nwkStr)
                else:
                    for child in self.children:
                        nwkStr = child.get_newick(nwkStr)
        
        if not len(self.children): return nwkStr

        possibleStrings = ['('+str(self.name)+')']
        possibleStrings += ['('+str(self.name)+',']
        possibleStrings += [','+str(self.name)+',']
        possibleStrings += [','+str(self.name)+')']
        childNames = sorted([c.name for c in self.children])
        childrenStr = ','.join([str(c) for c in childNames])
        for pStr in possibleStrings:
            if pStr in nwkStr:
                nwkStr = nwkStr.replace(pStr, pStr[0]+'('+childrenStr+')'+ \
                                        str(self.name)+pStr[-1])
                break
            ordered_children = sorted(self.children, key=lambda c: c.name)
            for child in ordered_children:
                nwkStr = child.get_newick(nwkStr)
        return nwkStr


#*******************************************************************
class MutationPair(object):
    '''
    Pair of mutations from tree.
    '''
   # __slots__ = ['mutationPair','commonPrecursors',
   #              'mut1precursors','mut2precursors','otherPrecursors',
   #              'commonFollowers',
   #              'mut1followers','mut2followers','otherFollowers',
   #              'commonLineage', 
   #              'mut1lineage','mut2lineage','otherLineage']
    
    def __init__(self,phyloEdgePairs=set(),mutationPair=None,
                 allEdges=set()):
        # List of pairs of phyloEdgeObjects
        if not len(phyloEdgePairs):
            raise ValueError("Set of PhyloEdge pairs must be"+\
                             " nonempty!")
        else:
            self.phyloEdgePairs = set(phyloEdgePairs)
        # A unique pair of AA mutations
        if mutationPair == None:
            print "WARNING: No mutation pair given. Set to"+\
                " first pair of mutations for first "+\
                " PhyloEdge pair!"
            self.mutationPair = (phyloEdgePairs[0][0].get_mutations()[0],
                                 phyloEdgePairs[0][1].get_mutations()[1])
        else:
            self.mutationPair = mutationPair

        self.allEdges = allEdges

    def assign_lineageDistribution(self):
        commonLineage, mut0lineage, mut1lineage = [],[],[]
        for pep in self.phyloEdgePairs:
            if pep[0]==pep[1]:
                commonLineage += [pep[0].parent] + [pep[0]] + pep[0].followers 
            else:
                commonPrec = pep[0].get_common_precursors(pep[1])
                commonFollow = list(set([pep[0]]+pep[0].followers) & \
                                    set([pep[1]]+pep[1].followers))
                commonLineage += commonPrec + commonFollow

                mut0prec = pep[0].get_intermediatePrec(commonPrec)
                mut1prec = pep[1].get_intermediatePrec(commonPrec)

                mut0lineage += list(set(mut0prec+pep[0].followers) - \
                                    set(commonPrec+commonFollow))
                mut1lineage += list(set(mut1prec+pep[1].followers) - \
                                    set(commonPrec+commonFollow))
        
        self.commonLineage = set([p for p in commonLineage if p.mutations])  
        self.mut0lineage = set([p for p in mut0lineage if p.mutations]) - \
                           self.commonLineage
        self.mut1lineage = set([p for p in mut1lineage if p.mutations]) - \
                           self.commonLineage

        assocLineage = self.commonLineage|self.mut0lineage|self.mut1lineage
        self.otherLineage = self.allEdges  - assocLineage

    
    def print_lineages(self):
        print "\tcommonLin:", [p.name for p in self.commonLineage]
        print "\tmut0Lin:", [p.name for p in self.mut0lineage]
        print "\tmut1Lin:", [p.name for p in self.mut1lineage]
        print "\totherLin:", [p.name for p in self.otherLineage]
        print 

    def get_contingencyTable(self):
        common = len(self.commonLineage)
        only0 = len(self.mut0lineage) 
        only1 = len(self.mut1lineage)
        other = len(self.otherLineage)
        self.contingencyTable = [[common, only0],[only1,other]]

    def get_contingencyTable_edges(self):
        occur01 = len(self.phyloEdgePairs)
        occur0 = len(self.mutationPair[0].phyloEdgePairs - self.phyloEdgePairs) 
        occur1 = len(self.mutationPair[1].phyloEdgePairs - self.phyloEdgePairs)
        occurOther = len(self.allPhyloEdgePairs\
                        -self.mutationPair[0].phyloEdgePairs \
                        -self.mutationPair[1].phyloEdgePairs)

        self.contingencyTable = [[occur01, occur0],[occur1,occurOther]]
        
    def get_FisherPvals(self):
        self.get_contingencyTable()
        self.oddsRatioLess, self.pValLess = \
                        stats.fisher_exact(self.contingencyTable,
                                           alternative='less')

        self.oddsRatioMore, self.pValMore = \
                        stats.fisher_exact(self.contingencyTable,
                                           alternative='greater')
class Mutation(object):
    def __init__(self,phyloEdgePairs=set(),mutationStr=None):
        self.mutationStr = mutationStr
        self.phyloEdgePairs = set(phyloEdgePairs)
    
    def add_phyloEdgePair(self,phyloEdgePair):
        self.phyloEdgePairs.add(phyloEdgePair) 
    
    def add_phyloEdgePairs(self,phyloEdgePairs):
        self.phyloEdgePairs.update(phyloEdgePairs)
    def __repr__(self):
        return self.mutationStr
#******************************************************************* 
class UndirPairSet(object):
    def __init__(self, nwkFile, fastaFile,posList,pos2wt,tip2pheno,phenoClades,
                 distThresh):
        self.nwkFile = nwkFile
        self.fastaFile = fastaFile
        self.posList = posList
        self.pos2wt = pos2wt
        self.node2pheno = tip2pheno.copy()
        self.allPhyloEdgePairs = set()
        # read the tree from newick(nwk) file and link nodes
        # to their corresponding fasta sequences
        stderr_write(["Reading in:\n", self.fastaFile, '\n', self.nwkFile])
        self.phyloTree = PhyloTree(newick=self.nwkFile,format=1)
        self.phyloTree.link_to_alignment(alignment=self.fastaFile,alg_format='fasta')

        self.phenoNameNodes = filter(lambda phyloNode: \
                                     phyloNode.name in phenoClades["2be"],
                                     self.phyloTree.get_leaves())

        self.get_phenoMuts()
        self.ancNode = self.phyloTree.get_common_ancestor(self.phenoNameNodes)
        stderr_write(['\tMRCA 2be',self.ancNode.name, '\t'+str(datetime.now())])



        intPhenoFile =  self.nwkFile.replace("int_nwk.tree","internalStates.txt")
        intPhenoFile =  intPhenoFile.replace("int_subPheno.tree","internalStates.txt")
        stderr_write(["intPhenoFile",intPhenoFile])

        with open(intPhenoFile,'rb') as pf:
            lineRecs = [line.strip().split() \
                        for line in pf.read().replace('"','').split('\n')]
            for recs in lineRecs:
                if len(recs)>1:
                    self.node2pheno[recs[0]] = recs[1]
        stderr_write(["Node 2 pheno number", len(self.node2pheno)])

        # The folling generates a linegraph object of the whole tree
        # it is a tree rooted in the edge between ancNode and its ancestor
        self.ancEdge = PhyloEdge(nodePair=[self.ancNode.up,self.ancNode],
                                 phenoDict = self.node2pheno.copy(),
                                 parent=None,children=[],posList=self.posList,
                                 pos2wt=self.pos2wt)
        self.ancEdge.build_EdgeTree()
        stderr_write(["Built edge tree!"])
        #stderr_write([self.ancEdge.get_newick()])

        # Now get all possible undirected pairs
        self.undirPairs = []
        stderr_write(["\tGetting all undirected pairs"])
        self.get_undirPairs(distThresh)
        stderr_write(["\tFinished extracting pairs:",self.nwkFile.split('/')[-1],
                      datetime.now()])

    def get_phenoMuts(self):
        standardAA="ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids

        mutations = []
        for node in self.phenoNameNodes:
            nodeSeq = node.sequence
            ancSeq = [self.pos2wt[p] for p in self.posList if p in self.pos2wt]
            mutations +=  [ancSeq[i]+str(self.posList[i])+nodeSeq[i] \
                        for i in range(len(nodeSeq)) if ancSeq[i]!=nodeSeq[i]\
                           and nodeSeq[i] in standardAA \
                           and self.posList[i]>23]
        self.phenoMutations = set(mutations)
        stderr_write(["phenotype mutations", self.phenoMutations])

    def get_undirPairs(self, distThresh):
        mutPair2phyloEdgePairs = defaultdict(list) 
        mut2phyloEdgePairs = defaultdict(list)
        allEdges = [self.ancEdge]+[fe for fe in self.ancEdge.followers]
        stderr_write(["\t\tNumber of tree edges:", 
                      len(allEdges)])
        
        stderr_write(["\t\tAssigning ordered precurors for all edges"])
        for edge in allEdges:
           edge.assign_ordered_precursors()
        
       
        for startEdge in self.ancEdge.followers:
            if startEdge.isLeafEdge(): continue
            if not startEdge.check_pheno(): continue

            if not startEdge.check_dist(startEdge,distThresh): continue
            

            startEdgeSortedMutations = sortByDigits(list(\
                                                    set(startEdge.mutations) & \
                                                         self.phenoMutations))
            if not len(startEdgeSortedMutations): continue

            #stderr_write([startEdge.name, startEdge.mutations])
            if len(startEdgeSortedMutations)>1:
                startEdgeMutationPairs = listPairs(startEdgeSortedMutations)
                # get mutations within same edge for all ancEdge+followers
                for mutPair in startEdgeMutationPairs:
                    mutPair2phyloEdgePairs[mutPair].append((startEdge,startEdge))
                    mut2phyloEdgePairs[mutPair[0]].append((startEdge,startEdge))
                    mut2phyloEdgePairs[mutPair[1]].append((startEdge,startEdge))
            
            # get mutation pairs for all start,follow edge combinations
            for followEdge in startEdge.followers:
                #if followEdge.isLeafEdge(): continue
                if not followEdge.check_pheno(): continue
                if not startEdge.check_dist(followEdge,distThresh): continue
               
                
                followEdgeSortedMutations = list((set(followEdge.mutations) &\
                                             self.phenoMutations)-\
                                            set(startEdge.mutations))
                
                if not len(followEdgeSortedMutations): continue
        
                for startMut in startEdgeSortedMutations:
                    for followMut in followEdgeSortedMutations:
                        if getInt(startMut) < getInt(followMut):
                            mutPair2phyloEdgePairs[(startMut,followMut)].\
                            append((startEdge,followEdge))
                            mut2phyloEdgePairs[startMut].\
                                append((startEdge,followEdge))
                            mut2phyloEdgePairs[followMut].\
                                append((startEdge,followEdge))

                        elif getInt(startMut) > getInt(followMut):
                            mutPair2phyloEdgePairs[(followMut,startMut)].\
                                append((followEdge,startEdge))
                            mut2phyloEdgePairs[startMut].\
                                append((followEdge,startEdge))
                            mut2phyloEdgePairs[followMut].\
                                append((followEdge,startEdge))


        stderr_write(["Mutations and pheno muts", mut2phyloEdgePairs.keys()])
        allMutations = set(mut2phyloEdgePairs.keys()) & set(self.phenoMutations)
        stderr_write(["\t\tNumber of AA mutations:",len(allMutations)])
        
        allAAmutEdges = set([edge for edge in self.ancEdge.followers \
                             if edge.check_pheno()])
        stderr_write(["\t\tNumber of edges with non-silent mutations for phenotype:",
                      len(allAAmutEdges)])
 
        self.mut2obj={}
        for mutPair in mutPair2phyloEdgePairs:
            if mutPair[0] not in self.mut2obj:
                self.mut2obj[mutPair[0]] = Mutation(mutationStr=mutPair[0])
            if mutPair[1] not in self.mut2obj:
                self.mut2obj[mutPair[1]] = Mutation(mutationStr=mutPair[1])
            
            self.mut2obj[mutPair[0]].add_phyloEdgePairs(set(\
                                                mutPair2phyloEdgePairs[mutPair]))
            self.mut2obj[mutPair[1]].add_phyloEdgePairs(set(\
                                                mutPair2phyloEdgePairs[mutPair]))

            undirMutPair = MutationPair(phyloEdgePairs = list(set(\
                                        mutPair2phyloEdgePairs[mutPair])),
                                        mutationPair = (self.mut2obj[mutPair[0]],
                                                        self.mut2obj[mutPair[1]]),
                                        allEdges=allAAmutEdges)
            undirMutPair.assign_lineageDistribution()
            self.undirPairs.append(undirMutPair)

        stderr_write(["\t\tPairwise lineage dsitributions assigned."])
        stderr_write(["\t\tBuild contingencies."])
        for mutPair in self.undirPairs:
            mutPair.get_contingencyTable()
                        

    def __repr__(self):
        outStr='\t'.join(['_Mut1_','_Mut2_','edgePairList',
                          'contTable(common_1_2_other)',
                          "FisherPvalLess","FisherPvalMore"])+'\n'
        
        for undirPairIdx in range(len(self.undirPairs)):
            undirPair = self.undirPairs[undirPairIdx]
            
            if len(undirPair.contingencyTable)<2:
                continue
            
            if len(flatten(undirPair.contingencyTable))<4:
                continue
            if undirPair.contingencyTable[0][1]==0:
                if undirPair.contingencyTable[1][0]==0:
                    if undirPair.contingencyTable[0][0]<=2:
                        continue
           
            toPrint = [str(undirPair.mutationPair[0]), 
                       str(undirPair.mutationPair[1])]
            toPrint.append(str(len(undirPair.phyloEdgePairs)))
            toPrint.append(str(undirPair.contingencyTable))
            toPrint.append(str(undirPair.pValLess))
            toPrint.append(str(undirPair.pValMore))
            outStr+= '\t'.join(toPrint)+'\n'

        return outStr
    
                                         



        
        
    
