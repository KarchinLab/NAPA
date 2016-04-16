#! /usr/bin/env python 
import sys, re
from collections import defaultdict
from datetime import datetime
from ete2 import *
from scipy import stats
from itertools import *
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
rstats = importr('stats')

def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

allowedPositions = set(range(300)) - set(range(23+1))        

flatten = lambda x: [subitem for item in x for subitem in item]
listPairs = lambda x: [(x[i1],x[i2]) for i1 in range(len(x)-1) \
                       for i2 in range(i1+1,len(x))]
getInt = lambda s: int(re.sub('\D','',s))
sortByDigits = lambda x: sorted(x, key = lambda s: getInt(s))
sortListbyOtherList = lambda x,y: [x for (y,x) in sorted(zip(Y,X),
                                                key=lambda pair: pair[0])]

standardAA="ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids

phenTransitions = {'2b':['2b','2be','2ber'],'2be':['2be'],'2br':['2b','2be','2ber'], 
                   '2ber':['2be','2ber']} # don't lose e, don't just gain r
allSequences = set()

def fasta_iter(fasta_name):
    """given a fasta file. yield tuples of header, sequence"""
    if fasta_name.endswith("gz"): fh = gzip.open(fasta_name)
    else: fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip() # drop the '>'
        seq = "".join(s.strip() for s in faiter.next()) 
        #join all sequence lines into one string
        yield header, seq

def read_fasta(FAinputFile):
    fasta_sequences = fasta_iter(FAinputFile)
    id2seq = {}
    for fs in fasta_sequences:
        id2seq[fs[0].replace(' ','|')]=fs[1]
    return id2seq

#*******************************************************************
class MutationPair(object):
    '''
    Pair of mutations from sequence.
    '''
    def __init__(self,sequences=set(),mutationPair=None):
        # List of pairs of phyloEdgeObjects
        if not len(sequences):
            raise ValueError("Set of PhyloEdge pairs must be"+\
                             " nonempty!")
        else:
            self.sequences = set(sequences)
        # A unique pair of AA mutations
        if mutationPair == None:
            print "WARNING: No mutation pair given. Set to"+\
                " first pair of mutations for first "+\
                " sequence!"
            self.mutationPair = (Mutation(sequences = self.sequences,
                                          mutationStr=\
                                          sequences[0].get_mutations()[0]),
                                 Mutation(sequences = self.sequences,
                                          mutationStr=\
                                          sequences[0].get_mutations()[1]))
        else:
            self.mutationPair = mutationPair

    def add_sequence(self,sequence):
        self.sequences.update(sequence) 

    def get_contingencyTable(self):
        occur01 = len(self.sequences)
        occur0 = len(self.mutationPair[0].sequences - self.sequences) 
        occur1 = len(self.mutationPair[1].sequences - self.sequences)
        occurOther = len(allSequences\
                        -self.mutationPair[0].sequences \
                        -self.mutationPair[1].sequences)

        self.contingencyTable = [[occur01, occur0],[occur1,occurOther]]
        self.jaccard = (float(occur01)-1./len(allSequences)) / (float(occur01+occur0+occur1))

    def get_FisherPvals(self):
        self.get_contingencyTable()
        self.oddsRatioLess, self.pValLess = \
                        stats.fisher_exact(self.contingencyTable,
                                           alternative='less')

        self.oddsRatioMore, self.pValMore = \
                        stats.fisher_exact(self.contingencyTable,
                                           alternative='greater')

class Mutation(object):
    def __init__(self,sequences=set(),mutationStr=None):
        self.mutationStr = mutationStr
        self.sequences = set(sequences)
    
    def add_sequence(self,sequence):
        self.sequences.add(sequence) 
    
    def add_sequences(self,sequences):
        self.sequences.update(sequences)
    def __repr__(self):
        return self.mutationStr
#******************************************************************* 
class UndirPairSet(object):
    def __init__(self, fastaFileList,posList,pos2wt,tip2pheno):
        self.fastaFileList = fastaFileList
        self.posList = posList
        self.pos2wt = pos2wt
        self.node2pheno = tip2pheno.copy()
        self.phenoNameNodes = [node for node in tip2pheno if tip2pheno[node] in ['2be']]
        print self.phenoNameNodes
        # read fasta file
        self.node2sequence = {}
        for ff in self.fastaFileList:
            node2sequence = read_fasta(ff)
            self.node2sequence.update(node2sequence)
            stderr_write([len(self.node2sequence), ff])
        self.node2sequenceNew = {}
        for node in self.node2sequence:
            if node not in tip2pheno: continue
            if tip2pheno[node] not in ['2be', '2ber']: continue
            self.node2sequenceNew[node] = self.node2sequence[node]
        self.node2sequence = self.node2sequenceNew


        # Get mutations for each node
        self.get_pheno_muts() #only count mutations associated with phenotype
        self.get_muts()

        # Now get all possible undirected pairs
        self.undirPairs = []
        self.get_undirPairs()

    def get_pheno_muts(self):
        standardAA="ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids

        self.phenoMutations = []
        for node in self.node2sequence: #mutations in both 2be, and 2ber
            if node not in self.phenoNameNodes: continue
            nodeSeq = self.node2sequence[node]
            ancSeq = [self.pos2wt[p] for p in self.posList if p in self.pos2wt]
            for i in range(len(nodeSeq)):
                if ancSeq[i]!=nodeSeq[i]:
                    if nodeSeq[i] in standardAA:
                        if self.posList[i] in allowedPositions:
                            self.phenoMutations.append(ancSeq[i]+str(self.posList[i])+\
                                                       nodeSeq[i])
    

    def get_muts(self):
        standardAA="ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids

        self.node2mutations = {}
        for node in self.node2sequence:
            nodeSeq = self.node2sequence[node]
            ancSeq = [self.pos2wt[p] for p in self.posList if p in self.pos2wt]
            self.node2mutations[node] = []
            for i in range(len(nodeSeq)):
                if ancSeq[i]!=nodeSeq[i]:
                    if nodeSeq[i] in standardAA:
                        if self.posList[i] in allowedPositions:
                            mut = ancSeq[i]+str(self.posList[i])+nodeSeq[i]
                            if mut not in self.phenoMutations: continue
                            self.node2mutations[node].append(mut)

    def get_undirPairs(self):
        global allSequences
        mutPair2sequences = defaultdict(list) 
        mut2sequences = defaultdict(list)
        for node in self.node2mutations:
                nodeMutationPairs = listPairs(self.node2mutations[node])
                for mutPair in nodeMutationPairs:
                    mutPair2sequences[mutPair].append(node)
                    mut2sequences[mutPair[0]].append(node)
                    mut2sequences[mutPair[1]].append(node)
        
        allMutations = set(mut2sequences.keys()) 
        stderr_write(["\t\tNumber of AA mutations:",len(allMutations)])
    
        self.mut2obj={}
        for mutPair in mutPair2sequences:
            allSequences.update(set(mutPair2sequences[mutPair]))
            if mutPair[0] not in self.mut2obj:
                self.mut2obj[mutPair[0]] = Mutation(mutationStr=mutPair[0])
            if mutPair[1] not in self.mut2obj:
                self.mut2obj[mutPair[1]] = Mutation(mutationStr=mutPair[1])
            
            self.mut2obj[mutPair[0]].add_sequences(set(\
                                                mutPair2sequences[mutPair]))
            self.mut2obj[mutPair[1]].add_sequences(set(\
                                                mutPair2sequences[mutPair]))

            undirMutPair = MutationPair(sequences = list(set(\
                                        mutPair2sequences[mutPair])),
                                        mutationPair = (self.mut2obj[mutPair[0]],
                                                        self.mut2obj[mutPair[1]]))
            self.undirPairs.append(undirMutPair)

        
        stderr_write(["\t\tNumber of sequences considered",
                      len(allSequences)])
        stderr_write(sorted(list(allSequences)))
        for mutPair in self.undirPairs:
            mutPair.get_contingencyTable()
                        

    def __repr__(self):
        outStr='\t'.join(['_Mut1_','_Mut2_','edgePairList',
                          'contTable(common_1_2_other)',
                          "FisherPvalLess","FisherPvalMore"])+'\n'
        
        for undirPairIdx in range(len(self.undirPairs)):
            undirPair = self.undirPairs[undirPairIdx]
            #if undirPair.contingencyTable[0][0] < 3: continue
            if (undirPair.contingencyTable[0][1] < 1) and \
                    (undirPair.contingencyTable[1][0] < 1):
                continue
            toPrint = [str(undirPair.mutationPair[0]), 
                       str(undirPair.mutationPair[1])]
            toPrint.append(str(len(undirPair.sequences)))
            toPrint.append(str(undirPair.contingencyTable))
            toPrint.append(str(undirPair.pValLess))
            toPrint.append(str(undirPair.pValMore))
            outStr+= '\t'.join(toPrint)+'\n'

        return outStr

    def str_jaccard(self):
        outStr='\t'.join(['_Mut1_','_Mut2_','edgePairList',
                          'contTable(common_1_2_other)',
                          'Fisher pval R', "OneMinusJaccard"])+'\n'
        
        for undirPairIdx in range(len(self.undirPairs)):
            undirPair = self.undirPairs[undirPairIdx]
            #if undirPair.contingencyTable[0][0] < 3: continue
            if (undirPair.contingencyTable[0][1] < 1) and \
                    (undirPair.contingencyTable[1][0] < 1):
                continue
            toPrint = [str(undirPair.mutationPair[0]), str(undirPair.mutationPair[1])]
            toPrint.append(str(len(undirPair.sequences)))
            toPrint.append(str(undirPair.contingencyTable))
            toPrint.append(str(undirPair.pValMore))
            toPrint.append(str(1.0-undirPair.jaccard))
            outStr+= '\t'.join(toPrint)+'\n'

        return outStr
                                         



        
        
    
