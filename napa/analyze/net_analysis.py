#!/usr/bin/env python
from collections import defaultdict

from napa.utils.config import Config
from napa.utils.serials import *
from napa.utils.io import *

from napa.net.net import *

class NetAnalysisInput(object):
    '''
    Prepare all inputs for phylogeny-based 
    network reconstruction.
    Builds on alignment inputs as provided for 
    alingment based network and adds 
    phylogeny specific parameters
    '''
    def __init__(self, config):
        self.__dict__.update(vars(config))
            
        if self.calculate_centralities:
            self.get_path_len_list()
            self.cent_list = defaultdict(list)
            for pl in self.path_len:
                self.get_cent_list(pl)
        
            self.cent_files = \
                [self.net_file.replace('.net.txt', 
                    '.'.join(['.cent_' + self.cent_type,
                              'path_len' + str(pl),
                              self.cent_rank_type[0:3],
                              'txt'])) \
                 for pl in self.path_len]

            stderr_write(['\nCentralities will be calculated',
                          'and printed to:\n\t' + \
                          '\n\t'.join(self.cent_files)])

        if self.cluster_nodes:
        
            self.node_clust_file = \
                self.net_file.replace('.net.txt', 
                                      '.communities.txt')

            stderr_write(['\nNode communities will be',
                          'printed to:\n',
                          self.node_clust_file])
            
    #--------------------------------------------------------#
    def get_path_len_list(self):

        if self.path_len == None:
            self.path_len = [1]
            stderr_write(['\nPath length not provided',
                'for path or rel. path centr.',
                '\nSetting it to default 1.'])
            
        if 'rel' in self.cent_rank_type.lower():
            path_len = [pl for pl in self.path_len if pl > 1]
            if not path_len or path_len != self.path_len:
                if not len(path_len):
                    self.path_len = [2]
                else:
                    self.path_len = path_len
                
                stderr_write(['\nWARNING: Path lengths',
                    '< 2 not defined for relative centralities.',
                    'Calculating centralities for path lengths:'] +\
                    self.path_len)

    #---------------------------------------------------------#
    def get_cent_list(self, path_len):
        '''
        Get all centralities that can be calculated based 
        on applicability to directed or undirected graphs.
        '''
       

        if 'loc' in self.cent_type or \
           'both' in self.cent_type:
            if self.edge_type.startswith('dir'):
                # Directed local centralities
                self.cent_list[path_len] +=  \
                ['loc.in.deg', 'loc.out.deg', 
                 'loc.in.strength', 
                 'loc.out.strength']
            else:
                # Undirected local centralities
                self.cent_list[path_len] += \
                ['loc.deg', 'loc.strength']

        if 'glob' in self.cent_type or \
           'both' in self.cent_type:
            self.cent_list[path_len] +=  \
                ['glob.close', 'glob.eigen', 
                 'glob.pagerank', 
                 'glob.betw', 'glob.kpath']

        if path_len > 1 and 'abs' in \
           self.cent_rank_type:

            self.cent_list[path_len] =  \
            ['glob.betw', 'glob.kpath']
            
                
        stderr_write(['\nWill calculate', 
                      self.cent_rank_type,
                      'centralities, for', 
                      self.cent_type,
                      'centrality types:\n'] + \
                     self.cent_list[path_len])

        if not len(self.cent_list[path_len]):
            stderr_write(['\nInvalid centrality',
                          'type provided',
                          'calculating undirected degree',
                          'centrality only.'])
            self.cent_list[path_len] = ['loc.deg']
            

#=========================================================#
def run_net_analysis(config):
    '''
    Processes network input information.
    Generates and writes network to text file(s).
    '''
    inp = NetAnalysisInput(config)
    mn = MutNet(net_file = inp.net_file, 
                net_type = inp.edge_type)

    if to_bool(inp.calculate_centralities):
        mn.get_centralities(\
            path_len_list = inp.path_len,
            cent_rank_type = inp.cent_rank_type,
            path_len_cent_list = inp.cent_list,
            outfiles = inp.cent_files)

    if to_bool(inp.cluster_nodes):
        mn.get_node_clusters(inp.node_clust_file)
    
