#!/usr/bin/env python
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

        if self.path_len == None:
            self.path_len = 1
            stderr_write(['Path length not provided',
                'for path or rel. path centr.',
                '\nSetting it to default 1.'])
            
        if self.path_len == 1 and \
           'rel' in self.cent_rank_type:
            self.cent_rank_type = 'absolute'
            
            stderr_write(['Centrality calculation',
                'was reset to absolute, because',
                'relative centralities not defined for',
                'paths of length 1.'])

        if self.calculate_centralities:
            self.get_cent_list()
        
            self.cent_file = \
                self.net_file.replace('.net.txt', 
                    '.'.join(['.cent_' + self.cent_type,
                    'path_len' + str(self.path_len),
                    self.cent_rank_type[0:3],'txt']))

            stderr_write(['Centralities will be calculated',
                          'and printed to:\n',
                          self.cent_file])
                             

    #---------------------------------------------------------#
    def get_cent_list(self):
        '''
        Get all centralities that can be calculated based 
        on applicability to directed or undirected graphs.
        '''
        self.cent_list = []

        if 'loc' in self.cent_type or \
           'both' in self.cent_type:
            if self.edge_type.startswith('dir'):
                # Directed local centralities
                self.cent_list +=  \
                ['loc.in.deg', 'loc.out.deg', 
                 'loc.in.strength', 
                 'loc.out.strength']
            else:
                # Undirected local centralities
                self.cent_list += \
                ['loc.deg', 'loc.strength']

        if 'glob' in self.cent_type or \
           'both' in self.cent_type:
            self.cent_list +=  \
                ['glob.close', 'glob.eigen', 
                 'glob.pagerank', 
                 'glob.betw', 'glob.kpath']

        if not len(self.cent_list):
            stderr_write(['Invalid centrality type provided',
                          'calculating undirected degree',
                          'centrality only.'])
            self.cent_list = ['loc.deg']
            

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
            path_len = inp.path_len,
            cent_rank_type = inp.cent_rank_type,
            cent_list = inp.cent_list,
            outfile = inp.cent_file)
                        
            
    elif cent_type == 'path':
        mn.get_path_between_path_cent(\
            path_node_length = args.pathLength)
        mn.write_path_betw_path_cent(args.centOutFile)

    elif 'rel' in cent_type: 
        if cent_type == 'rel':
            # all relative centralities
            centralities =  mn.get_default_cent_list()

        elif cent_type == 'relloc':
            # local centrality types only: 
            # degree, clustering...
            centralities =  \
            [c for c in mn.get_default_cent_list() \
             if 'loc' in c]

        elif cent_type == 'relglob':
            # global centrality types only: closeness, 
            # eigenvalue, betweenness...
            centralities =  \
            [c for c in mn.get_default_cent_list() \
             if 'glob' in c]


        header = 'node.or.path' + \
                 (args.pathLength - 1) * '\t' + \
                 '\t'.join(centralities) + \
                 '\tnum.net.nodes\n'

        out_str = mn.get_rel_cent(\
            cent_list = centralities, 
            path_len = args.pathLength)

        with open(args.centOutFile, 'wb') as cf:
            cf.write(header + out_str)


