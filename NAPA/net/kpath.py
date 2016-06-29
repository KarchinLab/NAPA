"""
k-path centrality
"""

import math, random, re
import networkx as nx
from collections import defaultdict
from networkx.utils import weighted_choice
from itertools import islice
from scipy.misc import comb

__all__ = ['kpath_node_centrality', 'kpath_path_centrality']

getInt = lambda s: int(re.sub('\D','',s))
getDiffs=lambda k,pl:sum([max(ki-pl,0.) for ki in range(1,k)])

#getScale = lambda k,n,pathLen: comb(n,pathLen)#/comb(n-pathLen,2)
getScale = lambda k,n,pathLen: 20.*comb(n,pathLen)/float(comb(n-pathLen,2))
#getScale = lambda k,n,pathLen: 10000./getDiffs(k,pathLen)
#getScale = lambda k,n,pathLen: 1000.*((k-pathLen+1.)*n)/comb(n-pathLen,2)
#getScale = lambda k,n,pathLen: ((k-pathLen+1.)*n) /comb(n-pathLen,2)
#getScale = lambda k,n,pathLen:  ( (k-pathLen+1.) * comb(n,pathLen) )

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

def kpath_node_centrality(G, k=None, alpha=0.2, weight=None, seed=None):
    if G.is_multigraph():
        raise nx.NetworkXError("Not implemented for multigraphs.")
    n = G.number_of_nodes()
    if k is None:
        k = math.log(n + G.number_of_edges())
    if k < 3: k=3

    T = 2 * k**2 * n**(1-2*alpha) * math.log(n)
    #print "kpath length:",k,'; T:',T

    random.seed(seed)
    centrality = dict.fromkeys(G, 0.0) 
    for i in range(int(T+1)):
        stack = []
        s = random.choice(G.nodes()) # choose source node
        l = random.randint(3,k) # choose a random path length
        stack.append(s)
        
        for j in range(l): # fill out a path of length l
            nbrs = dict((nbr,1.0/d.get(weight,1.0)) 
                        for nbr,d in G[s].items() if nbr not in stack)
            if not nbrs: break
            v = weighted_choice(nbrs)
            stack.append(v)
            s = v # set the current path source (current node) to v

        if len(stack) < 3: continue

        for node in stack[1:-1]:
            centrality[node] += 1

    for node in centrality:
        centrality[node] *= ((k-1.) * n) / T

    return centrality

def kpath_path_centrality(G, k=None, alpha=0.2, weight=None, seed=None, pathLen=1):

    if G.is_multigraph():
        raise nx.NetworkXError("Not implemented for multigraphs.")
    n = G.number_of_nodes()

    if k is None: k = int(math.log(n + G.number_of_edges()))
    if k < pathLen+2: k = pathLen+2
    T = 2. * k**2 * n**(1-2*alpha) * math.log(n)
    #print "kpath length:",k,'; T:',T
    
    random.seed(seed)
    centrality = defaultdict(float)
    for i in range(int(T+1)):
        stack = []

        s = random.choice(G.nodes()) # choose source node
        l = random.randint(pathLen+2,k) # choose a random path length
        stack.append(s)
        
        for j in range(l): # fill out a path of length l
            nbrs = dict((nbr,1.0/d.get(weight,1.0)) 
                        for nbr,d in G[s].items() if nbr not in stack)
            if not nbrs: break
            v = weighted_choice(nbrs)
            stack.append(v)
            s = v # set the current path source (current node) to v

        if len(stack) < pathLen+2: continue

        subPaths = [wn for wn in window(stack[1:-1],n=pathLen)]
        for subPath in subPaths:
            subPathStrList = [str(node) for node in subPath] 
            subPathPosList = [getInt(mut) for mut in subPathStrList]
            if pathLen == 1:
                centrality['_'.join(subPathStrList)] += 2
                continue
            if len(set(subPathPosList)) < len(subPathPosList):
                continue #repeated positions in path!
            centrality['_'.join(reversed(subPathStrList))] += 1
            centrality['_'.join(subPathStrList)] += 1

    for path in centrality:
        scale =  getScale(k,n,pathLen)
        centrality[path] *= (scale / T)*0.5

    return centrality

def kpath_probPath_centrality(G, k=None, alpha=0.2, weight=None, seed=None, pathLen=1):

    if G.is_multigraph():
        raise nx.NetworkXError("Not implemented for multigraphs.")
    n = G.number_of_nodes()

    if k is None: k = int(math.log(n + G.number_of_edges()))
    if k <= pathLen+1: k = pathLen+2
    T = 2 * k**2 * n**(1-2*alpha) * math.log(n)
    print "kpath length:",k,'; T:',T
    
    random.seed(seed)
    centrality = defaultdict(float)
    for i in range(int(T+1)):
        stack = []
        
        allOutDegrees = {}
        for v in G:
            weightSum = 0
            for nbr,d in G[v].items():
                weightSum += d.get(weight,1.0)
            allOutDegrees[v] = weightSum

        s = weighted_choice(allOutDegrees) # choose source node by degree
        l = random.randint(pathLen,k) # choose a random path length
        stack.append(s)
    
        for j in range(l): # fill out a path of length l
            nbrs = dict((nbr,1.0/d.get(weight,1.0)) 
                        for nbr,d in G[s].items() if nbr not in stack)
            if not nbrs: break
            v = weighted_choice(nbrs)
            stack.append(v)
            s = v # set the current path source (current node) to v

        if len(stack) < pathLen: continue

        subPaths = [wn for wn in window(stack,n=pathLen)]
        for subPath in subPaths:
            subPathStrList = [str(node) for node in subPath] 
            subPathPosList = [getInt(mut) for mut in subPathStrList]

            if pathLen == 1:
                centrality['_'.join(subPathStrList)] += 2
                continue

            if len(set(subPathPosList)) < len(subPathPosList):
                continue
            centrality['_'.join(reversed(subPathStrList))] += 1
            centrality['_'.join(subPathStrList)] += 1

    for path in centrality:
        scale =  getScale(k,n,pathLen)
        centrality[path] *= (scale / T)*0.5

    return centrality

