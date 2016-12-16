#!/usr/bin/env python
import os
import datetime
import yaml
from napa.utils.serials import *
from napa.utils.io import *

class Config():

    def __init__(self, path, command_type):
        '''
        Generate run configuration object from 
        user-provided config file in yaml format
        '''
        with open(path, 'r') as f:
            config = yaml.load(f)
        
        self.__dict__.update(config)
        self.set_run_type(command_type)
        print vars(self)
        self.check_input()

    #---------------------------------------------------------#
    def set_run_type(self, command_type):
        self.build_network = False
        self.analyze_network = False
        
        if 'build' in command_type:
            self.build_network = True
        
        if 'analy' in command_type:
            self.analyze_network = True

    #---------------------------------------------------------#
    def check_input(self):
        '''
        Checks input for consistently defined parameters.
        '''

        # Check main command and parameter options
        self.check_param_list(['net_type', 'edge_type'],
                              'Required command parameters')

        # Missing directories
        self.check_param_list(['working_dir', 'data_dir', 
                               'results_dir'], 
                              'Required Directories')
        self.set_net_file_base()
        self.set_net_file_name()

        # Check network reconstruction required parameters
        if self.build_network:
            self.build_net_input()

        # Check network analysis required parameters
        if self.analyze_network:
            self.analyze_net_input()

    #---------------------------------------------------------#
    def set_net_file_base(self):
        if not hasattr(self, 'net_file_base') and \
           not hasattr(self, 'net_file'):

            stderr_write(['No network name/base given.'])
            self.net_file_base = \
                str(datetime.datetime.now()).replace(' ',
                                        '_').replace(':','-')

            if hasattr(self, 'aln_fasta_file'):
                netname = os.path.basename(self.aln_fasta_file)
                netname = os.path.splitext(netname)[0]

            self.net_file_base += '_' + netname

            stderr_write(['Set network base name to:',
                          self.net_file_base])

    #---------------------------------------------------------#
    def set_net_file_name(self):
        if self.net_file == None or not \
           to_bool(self.net_file):

            suffix = '-' + self.net_type
            method_str = '-m_'
            if hasattr(self, 'method'):
                method_str += self.method
            else:
                method_str += 'na'
            suffix += method_str 

            thresh_str = '-t_'
            if hasattr(self, 'thresh'):
                thresh_str += str(self.thresh)
            else:
                thresh_str += 'na'
            suffix += thresh_str 

            if self.net_type == 'aln':
                occur_str = '-mo_'
                if hasattr(self, 'min_co_occur'):
                    occur_str += str(self.min_co_occur)
                suffix += occur_str 

            if self.net_type == 'phylo':

                edge_type_str = '-et_'
                if hasattr(self, 'edge_type'):
                    edge_type_str += self.edge_type
                suffix += edge_type_str

                dist_thresh_str = '-d_'
                if hasattr(self, 'dist_thresh'):
                    dist_thresh_str += \
                    str(self.dist_thresh)

                suffix += dist_thresh_str
                
                tree_support_str = '-mintrees_'
                if hasattr(self, 'tree_support_thresh'):
                    tree_support_str += \
                    str(self.tree_support_thresh)
                
                suffix += tree_support_str
            
            self.net_file = os.path.join(self.working_dir, 
                self.results_dir, self.net_file_base + \
                suffix + '.net.txt')

    #---------------------------------------------------------#
    def check_param_list(self, param_list, error_prefix = ''):
        '''
        Checks that all parameters in param_list are present
        and non-empty and first prints all missing parameters
        then throws a ValueError
        '''
        required_missing = False

        for param in param_list:
            if not hasattr(self, param):
                stderr_write(['Required parameter missing',
                              'from config:', param])
                required_missing = True
            elif self.__dict__[param] == None:
                stderr_write(['Empty value for required',
                              'parameter:', param])
                required_missing = True

        if required_missing:
            raise ValueError(error_prefix + " ERROR: " + \
                             "Missing or empty parameters!")
    #---------------------------------------------------------#
    def build_net_input(self):
        '''
        Checks and prepares input for building network 
        from alignment or phylogeny
        '''

        # REQUIRED for network construction 
        self.check_param_list(['method', 'thresh', 
                               'output_edge_weight_table', 
                               'wild_type', 'aln_fasta_file'], 
                              'Net Construct Input')
        
        if to_bool(self.output_edge_weight_table):
            self.net_table_file = \
            self.net_file.replace('.net.txt', '.table.txt')
                                  
        
        # REQUIRED: FASTA format alignment for alignment 
        # and phylo networks
        self.aln_fasta_file = \
            os.path.join(self.working_dir, self.data_dir,
                         self.aln_fasta_file)
        if not os.path.isfile(self.aln_fasta_file):
            raise ValueError('Please provide a ' + \
                             'valid FASTA file (path)!\n' + \
                             self.aln_fasta_file + \
                             'does not exist!')


        # OPTIONAL: List of positions/column numbers in alignment
        # Full path to position list file if provided
        if hasattr(self, 'pos_list_file') and \
           self.pos_list_file is not None:
            self.pos_list_file = \
                os.path.join(self.working_dir, self.data_dir,
                             self.pos_list_file)
        else:
            self.pos_list_file = ''

        # OPTIONAL: File with list of protein residue 
        # positions considered for network (column format)
        # (remove residues from signaling peptide, and/or
        # other regions not studied).
        if  hasattr(self, 'pos_subset_file') \
            and self.pos_subset_file is not None:
            self.pos_subset_file = \
                os.path.join(self.working_dir, self.data_dir,
                             self.pos_subset_file)
        else:
            self.pos_subset_file = ''

        # OPTIONAL: Ranges of protein residues to include
        # Can be used instead of pos_subset_file
        if not hasattr(self, 'protein_ranges') or not \
           to_bool(self.protein_ranges):
            self.protein_ranges = []
        


        # OPTIONAL: When using protein functional selection
        # for alignment or phylo leaf sequences, we need 
        # sequence functional annotation and 
        # list of selected functions.
        if hasattr(self, 'prot_func_file') and \
           hasattr(self, 'sel_prot_func_file'):
            if self.sel_prot_func_file != None and \
               self.prot_func_file == None:

                raise ValueError(' '.join([\
                    'Protein functional assignments needed',
                    'when selecting a subset by function.']))

            if self.sel_prot_func_file != None:
                self.sel_prot_func_file = \
                    os.path.join(self.working_dir,
                                 self.data_dir,
                                 self.sel_prot_func_file)
            else:
                self.sel_prot_func_file = ''

            if len(self.sel_prot_func_file): 
                if not os.path.isfile(self.sel_prot_func_file):
                    raise ValueError(' '.join([\
                        'If you are using the selected',
                        'protein function option',
                        'please provide a valid',
                        'selected prot func. file (path).']))
                
                self.prot_func_file = \
                    os.path.join(self.working_dir, self.data_dir,
                                 self.prot_func_file)
            
                if not os.path.isfile(self.prot_func_file):
                    raise ValueError(\
                        'Please provide a valid protein ' + \
                        'functional annotation file (path).')

        # REQUIRED when constructing the phylo network
        if 'phylo' in self.net_type:
            self.check_param_list(\
                ['dist_thresh', 'tree_support_thresh',
                 'tree_file_prefix_list',
                 'tree_nwk_file_suffix',
                 'tree_internal_node_seq_suffix'], 
                 'PhyloNet Construct Input')
            self.build_phylo_net_input()

    #---------------------------------------------------------#
    def build_phylo_net_input(self):
        '''
        Check and assign input dirs/params/opts for
        reconstructing phylogeny-based network
        '''
        # Functional transitions between phylogeny sequences 
        # allowed (evolution of specific function)
        self.func_transitions = {}
        if hasattr(self, 'func_transitions_file'):
            if self.func_transitions_file != None:
                self.func_transitions_file = \
                    os.path.join(self.working_dir, 
                                 self.data_dir,
                                 self.func_transitions_file)

                if not hasattr(self,
                        'tree_internal_node_state_suffix'):
                    raise ValueError(' '.join([\
                        'Internal node functions needed when',
                        'constraining by',
                        'functional transitions.']))
        
        if not len(self.tree_file_prefix_list):
           raise ValueError(' '.join([\
                'PhyloNet Build Input: Invalid newick tree'
                'directory \nand/or wrong tree name prefix!']))
    
        self.dist_thresh = get_float_from_str(self.dist_thresh, 
                                              -1.)
        if self.dist_thresh < 0.:
            stderr_write(['WARNING: Distance threshold',
                          'was not given.',
                          '\nNo distance constraints for pairs',
                          'of mutations will be applied.'])
        else:
            stderr_write(['Phylogeny distance threshold set to',
                          self.dist_thresh])

        # A mutation pair has to be significantly associated 
        # in at least this many ensemble trees, default = 2
        self.tree_support_thresh = \
                get_int_from_str(self.tree_support_thresh, 2)

    #---------------------------------------------------------#
    def analyze_net_input(self):
        if not hasattr(self, 'net_file'):
             raise ValueError(\
                    'Must provide network file when ' + \
                    'analyzing network!')

        if not os.path.isfile(self.net_file):
            raise ValueError(\
                    'Network file does not exist ' + \
                    'or is not a regular file:\n' + \
                    self.net_file)

        self.check_param_list(\
            ['cluster_nodes', 'calculate_centralities'], 
            'Required for netw. analysis')
        
        if not hasattr(self, 'cent_type') or \
           self.cent_type == None:
            self.cent_type = 'loc'

    #---------------------------------------------------------#
    def print_input_summary(self):
        stderr_write(['\nINITIAL INPUT PARAMETER SUMMARY:'])

        if to_bool(self.build_network):
            stderr_write(['Network will be reconstructed.'])
                
        if to_bool(self.analyze_network):
            stderr_write(['Network will be analyzed.'])

        stderr_write(['Network is', self.net_type, 'based.'])

        stderr_write(['Tab-delimited network file path:\n' + \
                      self.net_file])
        
        

        
