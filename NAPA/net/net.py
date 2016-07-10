import argparse
import sys
import os.path
from collections import defaultdict
from itertools import islice
import math
import random

import networkx as nx
from networkx.utils import weighted_choice

from NAPA.utils.general import *


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


def file_line_iterator(file_name):
    with open(file_name, 'rb') as f:
        for line in f: 
            yield line.strip()


def file_line_list(file_name):
    with open(file_name, 'rb') as f:
        return f.readlines()


class MutNet(object):

    '''
    Mutation network analysis methods
    '''

    def __init__(self, net_file = None, net_type = None):
        '''
        Networks based on tab-delimited files
        source\ttarget\tweight
        weight: weights from alignment/phylo networks
        net_type: 'directed' or 'undir'ected
        '''
        self.create_net(net_file, net_type)
        self.invert_weights()
        
    def create_net(self, net_file, net_type):
        '''
        Create networkX graph object (directed or undirected)
        '''

        print net_file
        if net_file == None or not os.path.isfile(net_file):
            raise MutationNetworkCreationException('Failed to create network'+ \
                                                   'Network file does not exist.')
        if net_type == 'undir':
            self.g = nx.parse_edgelist(file_line_list(net_file), 
                                       data=(('weight',float),), 
                                       create_using=nx.Graph())
            self.directed = False
        else:
            self.g = nx.parse_edgelist(file_line_list(net_file), 
                                       data=(('weight',float),), 
                                       create_using=nx.DiGraph())
            self.directed = True

    def invert_weights(self):
        '''
        Weight inversion to convert weights of association to distances between
        linked nodes. The higher the assoc., the lower the distance.
        Used in all shortest paths and path length methods.
        '''
        for u, v, d in self.g.edges(data=True):
            d['invWeight'] = 1. / d['weight']        

    def get_node_centralities(self, cent_list = ['degree', 'closeness', 
                                                 'betweenness','kpath']):
        '''
        Combined standard and extended cetntrality metrics for 
        single nodes. Needs directed versions of standard centralities
        already implemented in NetworkX.
        '''

        self.node_cent = {}

        if 'degree' in cent_list:
            self.node_cent['degree'] = nx.degree_centrality(self.g)
            # Need to add in/out-degree for directed network

        # Note: all (shortest) path-length based centralities use the inverse weight
        if 'closeness' in cent_list:
            self.node_cent['closeness'] = \
                nx.closeness_centrality(self.g, distance = 'invWeight', normalized = True)

        if 'betweenness' in cent_list:
            self.node_cent['betweenness'] = \
                nx.betweenness_centrality(weight = 'invWeight', normalized = True)

        if 'kpath' in cent_list:
            #see default options for path_kpath_centrality
            self.node_cent['k_path'] = self.path_kpath_centrality(alpha=0.)

            # get 0 k-path centrality nodes
        for node in self.g.nodes():
            if str(node) not in self.node_cent['k_path']:
                self.node_cent['k_path'][str(node)] = 0.

    def write_node_centralities(self, output_file):
        '''
        Not tested: write centralities for individual nodes in the network.
        '''
        with open(output_file, 'wb') as f:
            sorted_cent_names = sorted(self.node_cent.keys())
            f.write('\t'.join(['Node']+[cent for cent in sorted_cent_names]) + '\n')

            for node in self.g.nodes:
                f.write(node+'\t'+ '\t'.join([str(self.node_cent[cent][node]) \
                                              for cent in sorted_cent_names]) + '\n')

            
    def get_path_between_path_cent(self, path_node_length = 2,
                                          cent_list = ['shortest_path', 'k_path']):
        '''
        Get shortest and k-path betweenness centralities for all 
        paths consisting of path_node_length number of nodes.
        '''
        self.path_cent = {}
        
        if 'shortest_path' in cent_list:
            self.path_shortest_path_b_cent(path_len = path_node_length)

        if 'k_path' in cent_list:            
            self.path_kpath_centrality(alpha = 0., path_len = path_node_length)

            
    def write_path_betw_path_cent(self, output_file):
        with open(output_file, 'wb') as f:
            sorted_cent_names = sorted(self.path_cent.keys())
            f.write('\t'.join(['Path']+[cent+'_cent' for cent in sorted_cent_names]) + '\n')
            
            path_set = set()
            print self.path_cent
            for path_dict in self.path_cent.values():
                path_set |= set(path_dict.keys())

            for path in path_set:
                f.write(path + '\t' + \
                        '\t'.join([str(self.path_cent[cent][path]) if path in self.path_cent[cent] \
                                   else str(0.) for cent in sorted_cent_names]) + '\n')

    def get_gn_communities(self):
        self.components = nx.girvan_newman(self.g, weight='weight')


    def write_communities(self, output_file):
        '''
        Not tested: check that mutations are output rather than
        node numbers from networkX.
        '''
        with open(output_file, 'wb') as f:
            f.write('Node\tcommunity_label\n')

            for comp_i, component in enumerate(self.components):
                for node_name in component:
                    f.write('%s\t%d\n'%(node_name, comp_i + 1)) 
                                            
        

    def path_shortest_path_b_cent(self, path_len = 2):
        '''
        Extension of the single-node shortest path (sp) betweenness
        to a set of contiguous nodes in the network
        can be further optimized
        '''
        num_all_paths = 0
        self.path_cent['shortest_path'] = defaultdict(float)

        for source in self.g.nodes():
            for target in self.g.nodes():
                if source == target: 
                    continue #suboptimal, better iteration strategy?
                st_paths = [p for p in nx.all_shortest_paths(self.g, source, target)]
                num_all_paths += len(st_paths)
                for st_path in st_paths:
                    sub_paths = [wn for wn in window(st_path[1:-1], n = path_len)]
                    for sub_path in sub_paths:
                        sub_path_str_list = [str(node) for node in sub_path]

                        if path_len == 1:
                            self.path_cent['shortest_path']['_'.join(sub_path_str_list)] += 2
                            continue

                        if not self.directed:
                            self.path_cent['shortest_path']['_'.join(reversed(sub_path_str_list))] += 1
                            self.path_cent['shortest_path']['_'.join(sub_path_str_list)] += 1
                        else:
                            self.path_cent['shortest_path']['_'.join(sub_path_str_list)] += 2


        for path in self.path_cent['shortest_path']:
            self.path_cent['shortest_path'][path] *= (0.5 / num_all_paths)


    def path_kpath_centrality(self, k = None, alpha = 0.2, 
                              weight = 'invWeight', seed = 123456, path_len = 1):
        '''
        Adapted from Alakahoon et al., by extending from single node to path centrality
        and published code
        '''

        self.path_cent['k_path'] = defaultdict(float)

        if self.g.is_multigraph():
            raise nx.NetworkXError("Not implemented for multigraphs.")

        n = self.g.number_of_nodes()

        if k is None: 
            k = int(math.log(n + self.g.number_of_edges()))

        if k < path_len + 4: 
            k = min(path_len + 4, len(self.g.nodes()))

        T = 2. * k**2 * n**(1-2*alpha) * math.log(n)
        print "kpath length:",k,'; T:',T, 'path length (nodes)', path_len

        random.seed(seed)
        self.path_cent['k_path'] = defaultdict(float)

        for i in range(int(T+1)):
            st_path = []
            s = random.choice(self.g.nodes()) # choose source node
            l = random.randint(path_len + 2, k) # choose a random path length
            st_path.append(s)

            for j in range(l): # fill out a path of length l
                nbrs = {nbr: 1./d.get(weight,1.0) for nbr,d in self.g[s].items() \
                        if nbr not in st_path} #neighbors, weight needs to be inverted here too
                if not nbrs: break
                v = weighted_choice(nbrs)
                st_path.append(v)
                s = v # set the current path source (current node) to v

            if len(st_path) < path_len + 4: continue

            # sub_paths exclude end-points
            sub_paths = [wn for wn in window(st_path[1:-1], n = path_len)]
            for sub_path in sub_paths:
                sub_path_str_list = [str(node) for node in sub_path] 

                if path_len == 1:
                    self.path_cent['k_path']['_'.join(sub_path_str_list)] += 2.
                    continue

                if not self.directed:
                    self.path_cent['k_path']['_'.join(reversed(sub_path_str_list))] += 1.
                    self.path_cent['k_path']['_'.join(sub_path_str_list)] += 1.
                else:
                    self.path_cent['k_path']['_'.join(sub_path_str_list)] += 2.


        for path in self.path_cent['k_path']:
            self.path_cent['k_path'][path] *= (0.5 / T)

