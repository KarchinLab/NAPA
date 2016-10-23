import os
from collections import defaultdict
from copy import deepcopy

import math
import random

import networkx as nx
from networkx.utils import weighted_choice

from napa.utils.serials import *
from napa.utils.io import *


class MutNet(object):

    '''
    Mutation network analysis methods
    '''

    def __init__(self, net = None, net_file = None, 
                 net_type = None):
        '''
        Networks based on tab-delimited files
        source\ttarget\tweight
        weight: weights from alignment/phylo nets
        net_type: 'directed' or 'undir'ected
        '''
        if net != None:
            self.from_nx_net(net, net_type)
        else:
            self.read_net(net_file, net_type)

        # inverted edges net for some centralities
        self.get_reverse_net()
        
        # inverted weights for centralities 
        #requiring distance
        self.get_inv_w()
        
    #------------------------------------------------------#
    def from_nx_net(self, net, net_type):
        try:
            self.g = deepcopy(net)
            self.directed = nx.is_directed(net)
        except:
            raise Exception('NetworkX format error.')

    #------------------------------------------------------#
    def get_reverse_net(self):
        if self.directed:
            self.rg = nx.reverse(self.g) 
        else:    
            self.rg = deepcopy(self.g)
   
    #------------------------------------------------------#
    def read_net(self, net_file, net_type):
        '''
        Read from file
        Create networkX graph object 
        (directed or undirected)
        '''
        if net_file == None or not os.path.isfile(net_file):
            stderr_write(['WARN: Missing net file. '+ \
                          'Empty network created.'])
            if net_type == 'undir':
                self.g = nx.Graph()
                self.directed = False
            else:
                self.g = nx.DiGraph()
                self.directed = True
            return

        if net_type == 'undir':
            self.g = \
            nx.parse_edgelist(file_line_iterator(net_file), 
                              data=(('weight',float),), 
                              create_using=nx.Graph())
            self.directed = False
        
        else:
            self.g = \
            nx.parse_edgelist(file_line_iterator(net_file), 
                              data=(('weight',float),), 
                              create_using=nx.DiGraph())

            self.directed = True

    #------------------------------------------------------#
    def get_inv_w(self):
        '''
        Weight inversion to convert weights of association to
        distances between linked nodes. 
        The higher the association, the lower the distance.
        Used in all shortest paths and path length methods.
        '''

        if not len(self.g): return
        if not len(nx.get_edge_attributes(self.g, 
                                          'weight').values()):
            return
        max_weight = float(max([d['weight'] for u, v, d in \
                                self.g.edges(data=True)]))
        min_weight = float(min([d['weight'] for u, v, d in \
                                self.g.edges(data=True)]))

        nx.set_edge_attributes(self.g, 'weight',
            {(u, v):d['weight']/max_weight \
             for u, v, d in self.g.edges(data=True)})

        nx.set_edge_attributes(self.g, 'invWeight',
            {(u, v):(min_weight/max_weight)/d['weight'] \
             for u, v, d in self.g.edges(data=True)})


        nx.set_edge_attributes(self.rg, 'weight',
            {(u, v):d['weight']/max_weight \
             for u, v, d in self.rg.edges(data=True)})

        nx.set_edge_attributes(self.rg, 'invWeight',
            {(u, v):(min_weight/max_weight)/d['weight'] \
             for u, v, d in self.rg.edges(data=True)})
            
    #------------------------------------------------------#
    def normalize(self, node_att_dict = {}, norm_factor = 1.):
        norm_factor = 1. if norm_factor == 0. else norm_factor
        d = defaultdict(float)
        for k in node_att_dict:
            d[k] = float(node_att_dict[k]) / norm_factor
        return d

    #------------------------------------------------------#
    def get_centralities(self, path_len, cent_rank_type,
                         cent_list, outfile):
        if path_len == 1:
            self.get_node_centralities(cent_list)
            
            write_table_str(outfile, '',
                self.str_node_centralities(header = True))

        elif path_len > 1:
            if 'abs' in cent_rank_type:
                self.get_path_between_path_cent(\
                        path_node_length = path_len)
                write_table_str(outfile, '',
                    self.str_path_betw_path_cent())

            elif 'rel' in cent_rank_type:
                header = 'node.or.path' + \
                 (path_len - 1) * '\t' + \
                 '\t'.join(cent_list) + \
                 '\tnum.net.nodes\n'
                
                out_str = \
                self.get_rel_cent(cent_list = cent_list, 
                                  path_len = path_len)
                write_table_str(outfile, header, out_str)

    #------------------------------------------------------#
    def get_node_centralities(self, cent_list):
        '''
        Combined standard (local and shortest-path)
        and kpath centrality metrics for 
        single nodes. 
        '''

        if 'loc.in.deg' in cent_list:
            norm_factor = self.g.number_of_nodes() - 1.
            # number of in neighbors 
            # (disreguard link weights)
            # normalize each node by remaining number 
            # of nodes in network
            nx.set_node_attributes(self.g, 'loc.in.deg', 
                self.normalize(self.g.in_degree(weight=None),
                                             norm_factor))

        norm_factor = self.g.number_of_nodes() - 1.

        if 'loc.out.deg' in cent_list:
            nx.set_node_attributes(self.g, 'loc.out.deg', 
                self.normalize(self.g.out_degree(weight=None),
                               norm_factor))
        
        if 'loc.deg' in cent_list:
            nx.set_node_attributes(self.g, 'loc.deg', 
                self.normalize(\
                    self.g.degree(weight=None),
                               norm_factor))

        if 'loc.strength' in cent_list:
            nx.set_node_attributes(self.g, 'loc.strength', 
                self.normalize(\
                    self.g.degree(weight = 'weight'),
                               norm_factor))

        if 'loc.in.strength' in cent_list:
            nx.set_node_attributes(self.g, 'loc.in.strength', 
                self.normalize(\
                    self.g.in_degree(weight = 'weight'),
                               norm_factor))

        if 'loc.out.strength' in cent_list:
            nx.set_node_attributes(self.g, 
                'loc.out.strength', 
                self.normalize(\
                    self.g.out_degree(weight = 'weight'),
                                            norm_factor))
        #------------------------------------------------------#
        # Note: all (shortest) path-length based centralities 
        # use the inverse weight -- 
        # this is a proxy for distance, rather than assoc.,
        # s.t. stronger association weight = shorter dist.
        if 'glob.close' in cent_list:
            nx.set_node_attributes(self.g, 'glob.close', 
                nx.closeness_centrality(self.g, 
                                        distance = 'invWeight', 
                                        normalized = True))

        if 'glob.info' in cent_list:
            nx.set_node_attributes(self.g, 'glob.info',
                nx.current_flow_closeness_centrality(self.g, 
                                        weight = 'weight'))

        if 'glob.eigen' in cent_list:
            try:
                ecent = nx.eigenvector_centrality(self.rg, 
                    max_iter = 500, tol = 1e-05,
                    weight = 'weight')

            except nx.exception.NetworkXError:
                ecent = {n:0. for n in self.g.nodes()}
            
            nx.set_node_attributes(self.g, 'glob.eigen',
                                   ecent)
                                       

        if 'glob.pagerank' in cent_list:
            nx.set_node_attributes(self.g, 'glob.pagerank',
                nx.pagerank(self.rg, weight = 'weight'))

        # betweenness is not distance based in networkx
        if 'glob.betw' in cent_list:
            nx.set_node_attributes(self.g, 'glob.betw',
                nx.betweenness_centrality(self.g, 
                    weight = 'weight', normalized = True))

        if 'glob.kpath' in cent_list:
            node_kpath_cents = \
            self.path_kpath_centrality(alpha = 0.)

            for node in self.g.nodes():
                if str(node) in node_kpath_cents:
                    self.g.node[node]['glob.kpath'] = \
                    node_kpath_cents[str(node)]
                else: 
                    self.g.node[node]['glob.kpath'] = 0.

  
    #------------------------------------------------------#
    def str_node_centralities(self, header = False, 
                              prefix = '', suffix = '',
                              nodes = None, 
                              sorted_cent_names = None):
        '''
        Write centralities for (individual) nodes 
        in the network.
        '''
        
        if sorted_cent_names == None:

            sorted_cent_names = \
            flatten([self.g.node[node].keys() for node in \
                     self.g.nodes_iter()])
            sorted_cent_names = \
            sorted(list(set(sorted_cent_names)))
        
        out_list = []
        if header:
            header = \
            'node.or.path\t' + (prefix + suffix) + '\t' + \
            '\t'.join(sorted_cent_names)+'\tnum.net.nodes'
            
            out_list += [header]

        if nodes != None:
            for node in set(nodes) & set(self.g.nodes()):
                out_list +=  \
                ['%s\t%s%s\t'%(prefix, str(node), suffix) + \
                 '\t'.join([\
                    '{:.4f}'.format(self.g.node[node][cent]) \
                    for cent in sorted_cent_names]) + \
                 '\t' + str(self.g.number_of_nodes())]

        else:
            for node in self.g.nodes():
                out_list += \
                [prefix + str(node) + suffix + '\t'+ \
                 '\t'.join(['{:.4f}'.format(\
                    self.g.node[node][cent]) \
                    for cent in sorted_cent_names]) + '\t' + \
                 str(self.g.number_of_nodes())]

        if not len(out_list): return ''

        return '\n'.join(out_list) + '\n'


    #------------------------------------------------------#
    def get_desc_net_node_cent(self, source = None, 
                               exclude_nodes = [],
                               cent_list = None):

        nodes_in_graph =  nx.descendants(self.g, source)
        nodes_in_graph |= set([source])
        nodes_in_graph -=  set(exclude_nodes)

        if not len(nodes_in_graph): 
            desc_net = MutNet()
        
        desc_net = \
        MutNet(net = self.g.subgraph(nodes_in_graph))
        
        desc_net.get_node_centralities(cent_list = cent_list)
        return desc_net

    #------------------------------------------------------#
    def get_rel_cent(self, cent_list = None, prev_nodes = [], 
                     nodes = None, path_len = 3):
        
        if len(prev_nodes) >= path_len: return ''
        
        if nodes != None and not len(nodes): nodes = None

        
        prefix = '\t'.join(prev_nodes)
        suffix = (path_len - len(prev_nodes) - 1)*'\t----'
        
        self.get_node_centralities(cent_list = cent_list)
        out_str = \
        self.str_node_centralities(prefix = prefix,
                                   suffix = suffix,
                                   nodes = nodes)

        for node in self.g:

            if node in prev_nodes: continue
            desc_net = \
            self.get_desc_net_node_cent(source = node,
                            exclude_nodes = prev_nodes,
                                cent_list = cent_list)

            if desc_net.g.size() <= 3: continue

            out_str += \
            desc_net.get_rel_cent(cent_list = cent_list, 
                prev_nodes = prev_nodes + [node],
                nodes = list(set(self.g.neighbors(node)) - \
                             set(prev_nodes)),
                path_len = path_len)

        return out_str

    #------------------------------------------------------#
    def get_path_between_path_cent(self, path_node_length = 2,
                    cent_list = ['glob.betw', 'glob.kpath']):
        '''
        Get shortest and k-path betweenness centralities for 
        all paths consisting of path_node_length number of nodes.
        Here path_len refers to number of nodes rather than 
        edges in path.
        '''
        self.path_cent = {}
        
        if 'glob.betw' in cent_list:
            self.path_shortest_path_b_cent(\
                                path_len = path_node_length)

        if 'glob.kpath' in cent_list:            
            self.path_kpath_centrality(alpha = 0., 
                        path_len = path_node_length)

    #------------------------------------------------------#
    def str_path_betw_path_cent(self):

        sorted_cent_names = sorted(self.path_cent.keys())
        out_str =  '\t'.join(['Path'] + \
            [cent for cent in sorted_cent_names]) + '\n'
            
        path_set = set()
        for path_dict in self.path_cent.values():
            path_set |= set(path_dict.keys())

        for path in path_set:
            out_str += path + '\t' + \
                '\t'.join([str(self.path_cent[cent][path]) \
                           if path in self.path_cent[cent] \
                           else str(0.) for cent in \
                           sorted_cent_names]) + '\n'
        return out_str

    #------------------------------------------------------# 
    def get_gn_communities(self): 
        self.components =  girvan_newman(self.g, 
                                         weight='weight')
        
    #------------------------------------------------------#
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
                                            
    #------------------------------------------------------#
    def all_sps(self, source, target):
        '''
        No error all_shortest_paths between source and target
        '''
        try:
            return \
            [p for p in nx.all_shortest_paths(self.g, 
                        source, target, weight = 'invWeight')]

        except nx.NetworkXNoPath:
            return []

    #------------------------------------------------------#
    def path_shortest_path_b_cent(self, path_len = 2):
        '''
        Extension of the single-node shortest path (sp)
        betweenness to a set of contiguous nodes 
        in the network (can be further optimized)
        '''
        num_all_paths = 0
        self.path_cent['shortest_path'] = defaultdict(float)

        for source in self.g.nodes():

            for target in self.g.nodes():
                if source == target: 
                    continue #better iteration strategy?

                st_paths = self.all_sps(source, target)    
                num_all_paths += len(st_paths)
                for st_path in st_paths:
                    sub_paths = \
                    [wn for wn in window(st_path[1:-1], 
                                         n = path_len)]

                    for sub_path in sub_paths:
                        sub_path_str_list = \
                        [str(node) for node in sub_path]

                        if path_len == 1:
                            self.path_cent['shortest_path'][\
                                '_'.join(sub_path_str_list)] \
                                += 2
                            continue

                        if not self.directed:
                            self.path_cent['shortest_path'][\
                            '_'.join(\
                                reversed(sub_path_str_list))] \
                                += 1

                            self.path_cent['shortest_path'][\
                            '_'.join(sub_path_str_list)] += 1

                        else:
                            self.path_cent['shortest_path'][\
                                '_'.join(sub_path_str_list)] \
                                += 2

        for path in self.path_cent['shortest_path']:
            self.path_cent['shortest_path'][path] *= \
                                    (0.5 / num_all_paths)


    #------------------------------------------------------#
    def path_kpath_centrality(self, k = None, alpha = 0.2, 
        weight = 'invWeight', seed = 123456, path_len = 1):
        '''
        Adapted from Alakahoon et al. and published code
        by extending from single node to path centrality
        '''

        kpath_path_cent = defaultdict(float)

        if self.g.is_multigraph():
            raise nx.NetworkXError("Not implemented" + \
                                   "for multigraphs.")

        n = self.g.number_of_nodes()
        ne = self.g.number_of_edges()

        if n <= path_len + 4 or ne < path_len + 4:
            # network is too small for calculating
            # path centr. of this length
            return {}
            
        if k is None: 
            k = int(math.log(n + ne))

        if k < path_len + 4: k = min(path_len + 4, n)
            
        T = 2. * k**2 * n**(1-2*alpha) * math.log(n)
        
        random.seed(seed)

        for i in range(int(T+1)):
            st_path = []
            
            # choose source node
            s = random.choice(self.g.nodes())

            # choose a random path length
            l = random.randint(path_len + 2, k) 
            st_path.append(s)

             # fill out a path of length l
            for j in range(l):
                # invert distance again - 
                # get association
                nbrs = {nbr: 1./d.get(weight, 1.0) \
                        for nbr,d in self.g[s].items() \
                        if nbr not in st_path} 

                if not nbrs: break

                v = weighted_choice(nbrs)
                st_path.append(v)

                # set the current path source 
                # (current node) to v
                s = v 

            if len(st_path) < path_len + 4: continue

            # sub_paths exclude end-points
            sub_paths = [wn for wn in window(st_path[1:-1], 
                                             n = path_len)]

            for sub_path in sub_paths:
                sub_path_str_list = \
                [str(node) for node in sub_path] 

                if path_len == 1:
                    kpath_path_cent[\
                        '_'.join(sub_path_str_list)] += 2.
                    continue

                if not self.directed:
                    kpath_path_cent['_'.join(\
                        reversed(sub_path_str_list))] += 1.

                    kpath_path_cent[\
                        '_'.join(sub_path_str_list)] += 1.

                else:
                    kpath_path_cent[\
                        '_'.join(sub_path_str_list)] += 2.


        for path in kpath_path_cent:
            kpath_path_cent[path] *= (0.5 / T)

        return kpath_path_cent
