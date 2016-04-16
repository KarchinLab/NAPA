from igraph import *
import networkx as nx
import argparse
import os.path
from os import listdir
from collections import defaultdict
import math
import re
import scipy, scipy.stats
from scipy.stats import rankdata
import pickle
import numpy as np
from  MutNetAnalysis_networkx import *
from  MutNetAnalysis import *
from combinePathPredictions import *

flatten = lambda x: [subitem for item in x for subitem in item]
getInt = lambda s: int(re.sub('\D','',s))
rankDict = lambda l: dict( zip(l, rankdata([-e for e in l])) )
test=False
maxNets, maxPath = 5, 6

def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest = "randomNetDir", type=str, required=True)
    parser.add_argument('-n', dest = "numSwaps", type=int, required=True)
    parser.add_argument('-p', dest = "pathLength", type=int, required=True)
    args = parser.parse_args()
    return args

def read_randIgraph(randDir,numSwaps):
    randFiles=sorted([ randDir+'/'+f for f in listdir(randDir) \
                       if f.strip().endswith('swp'+str(numSwaps)+'.pyPickle')])
    if test:
        return [MutGraph(dirNetFile=rf,fromGraph=True) for rf in randFiles[0:maxNets]]
    return [MutGraph(dirNetFile=rf,fromGraph=True) for rf in randFiles]


def convert_igraph2NX(ig):
    #nxg = nx.DiGraph(ig.get_edgeList()) #will need node names
    # add weights
    dod = {}
    for e in ig.get_edgelist():
        source, target = ig.vs[e[0]]['name'],ig.vs[e[1]]['name']
        eObj = ig.es[ig.get_eid(e[0],e[1])]
        weight = eObj['weight']
        invWeight = eObj['invWeight']
        if source not in dod:
            dod[source] = {target:{'weight':weight, 'invWeight':invWeight}}
        elif target not in dod[source]:
                dod[source][target]={'weight':weight,'invWeight':invWeight}
    nxg = nx.DiGraph(dod)       
    return nxg

def get_allSPbetweenCents(mgList,pathLength,minBetw):
    p2allBetw = {}
    for i,mg in enumerate(mgList):
        p2mp = dir_path_betweenness(mg,pathLength)
        p2betw = dict([(p,p2mp[p].centrality) for p in p2mp])
        for p in p2betw:
            if p not in p2allBetw:
               p2allBetw[p]=[minBetw for mi in range(len(mgList))]
            p2allBetw[p][i]=p2betw[p]
    return p2allBetw

def get_RWbetweenCent(g,pathLength):
    if test:
        return kpath_path_centrality(g, k=maxPath,seed=12345,alpha=0.2,
                                     weight='invWeight',pathLen=pathLength)
    return kpath_path_centrality(g, k=(pathLength+2)*2,seed=12345,alpha=0.2,
                                 weight='invWeight',pathLen=pathLength)

def get_RWprobCent(g,pathLength):
    if test:
        return kpath_probPath_centrality(g, k=maxPath,seed=12345,alpha=0.2,
                                     weight='invWeight',pathLen=pathLength)
    return kpath_probPath_centrality(g, k=(pathLength+2)*2,seed=12345,alpha=0.2,
                                 weight='invWeight',pathLen=pathLength)


def get_allRWbetweenCents(ngList,pathLength,minBetw):
    print "---Getting RW betweenness."
    p2allCent={}
    for i,ng in enumerate(ngList):
        p2cent = get_RWbetweenCent(ng,pathLength)
        for p in p2cent:
            if p not in p2allCent:
                p2allCent[p]=[minBetw for mi in range(len(ngList))]
            p2allCent[p][i]=p2cent[p]
    return p2allCent

def get_allRWprobCents(ngList,pathLength,minProb):
    print "---Getting RW probable."
    p2allCent={}
    for i,ng in enumerate(ngList):
        p2cent = get_RWprobCent(ng,pathLength)
        for p in p2cent:
            if p not in p2allCent:
                p2allCent[p]=[minProb for mi in range(len(ngList))]
            p2allCent[p][i]=p2cent[p]
    return p2allCent


def get_path2rpval(path2cents,path,cent):
    return sum(np.array(path2cents[path]) >= cent)/float(len(path2cents[path]))
    

def main():
    args = parseArgs()
    defaultBetween=0.0
    mgs = read_randIgraph(args.randomNetDir,args.numSwaps)
    stderr_write(["Read in", len(mgs), "random nets"])

    sp_p2cent = get_allSPbetweenCents(mgs,args.pathLength,defaultBetween)
    stderr_write(["Got betweenness centralities for",len(sp_p2cent),"paths"])
    pickle.dump(sp_p2cent, open("SPbetweenCent.pathLen"+str(args.pathLength)+\
                                ".numNets"+str(len(mgs))+".numSwp"+\
                                str(args.numSwaps)+'.pyPickle','wb'))

    networkxRandNets = [convert_igraph2NX(mgr.g) for mgr in mgs]
    stderr_write(["Converted",str(len(mgs)),"random nets to networkX."])
    rw_p2cent = get_allRWbetweenCents(networkxRandNets,args.pathLength, 
                                       defaultBetween)
    stderr_write(["Got RW betweenness centralities for",len(rw_p2cent),"paths"])
    pickle.dump(rw_p2cent, open("RWbetweenCent.pathLen"+str(args.pathLength)+\
                                ".numNets"+str(len(mgs))+".numSwp"+\
                                str(args.numSwaps)+'.pyPickle','wb'))

    rwp_p2cent = get_allRWprobCents(networkxRandNets,args.pathLength, 
                                       defaultBetween)
    stderr_write(["Got RW prob centralities for",len(rwp_p2cent),"paths"])
    pickle.dump(rwp_p2cent, open("RWprobCent.pathLen"+str(args.pathLength)+\
                                ".numNets"+str(len(mgs))+".numSwp"+\
                                str(args.numSwaps)+'.pyPickle','wb'))

    if test:
        print "compare centr"
        for p in  sorted(list(set(sp_p2cent.keys()) & set(rwp_p2cent.keys())))[0:10]:
            print p, [b for b in sp_p2cent[p] if b > defaultBetween],\
                [b for b in rwp_p2cent[p] if b> defaultBetween]


if __name__ == "__main__": main()
