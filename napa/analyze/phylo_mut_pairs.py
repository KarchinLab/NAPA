from collections import defaultdict

from napa.utils.io import *
from napa.utils.serials import *

from napa.mutpair.phylo import * 
from napa.analyze.aln_mut_pairs import *

class PhyloNetInput(AlnNetInput):
    '''
    Prepare all inputs for phylogeny-based 
    network reconstruction.
    Builds on alignment inputs as provided for 
    alingment based network and adds 
    phylogeny specific parameters
    '''
    def __init__(self, config):
        # (Leaf sequence alignment inputs)
        super(PhyloNetInput, self).__init__(config)

        # Phylogeny specific inputs
        self.func_transitions = {}
        if len(self.func_transitions_file):
            self.func_transitions = \
                parse_keyval_dict(self.func_transitions_file)
            stderr_write(['Using', len(self.func_transitions),
                          'functional transitions'])
        else:
            stderr_write(['Functional transitions between',
                          'internal node functions when',
                          'traversing branches',
                          'not constrained.'])
            
        self.get_tree_files()

    def get_tree_files(self):

        self.tree_files =  defaultdict(list)
        self.int_seq_files = defaultdict(list)
        self.int_func_files = defaultdict(list)

        for prefix in self.tree_file_prefix_list:
            prefix = os.path.join(self.working_dir,
                                  self.data_dir, prefix)

            self.tree_files[prefix] = \
                list_display_files('Nwk trees with internal '+ \
                                   'nodes', prefix, 
                                self.tree_nwk_file_suffix)

            self.int_seq_files[prefix] = list_display_files(\
                'Internal node sequences', prefix, 
                self.tree_internal_node_seq_suffix)

            if hasattr(self,'tree_internal_node_state_suffix'):
                if self.tree_internal_node_state_suffix != None:
                    self.int_func_files[prefix] = \
                        list_display_files('Reconstructed ' + \
                            'functions on internalmnodes', 
                                           prefix, 
                            self.tree_internal_node_state_suffix)

            if len(self.int_seq_files[prefix]) != \
               len(self.tree_files[prefix]):
                raise ValueError(\
                    ' '.join(['PhyloNet Build Input:',
                              'Check that each Newick tree file',
                              'has a corresponding',
                              'internal states file!']))

            if len(self.int_func_files[prefix]):
                if len(self.int_func_files[prefix]) != \
                   len(self.tree_files[prefix]):
                    raise ValueError(\
                        ' '.join(['PhyloNet Build Input:',
                                  'Check that each Newick tree',
                                  'file has a corresponding', 
                                  'internal states file!']))


def run_phylo_mut_pairs(config):
    '''
    Processes network input information.
    Generates and writes network to text file(s).
    '''
    inp = PhyloNetInput(config)

    tree_ensemble = TreeEnsembleMutPairs(inp)
    mut_pair_dicts = get_mut_pair_dicts(tree_ensemble)
    tree_ensemble.assign_mut_pair_weights(mut_pair_dicts)
    phylo_mut_pairs = tree_ensemble
    
    # Output network in tab-delimited format
    # Source\tTarget\tWeight
    stderr_write(['\nWriting network to file:\n' + inp.net_file]) 
    phylo_mut_pairs.write_network_to_file(inp.net_file)

    # outputs additional information about edge weights
    # to a separate table
    if to_bool(inp.output_edge_weight_table):
        stderr_write(['\nWriting network extended network',
                      'table to file:\n' + inp.net_table_file]) 
        phylo_mut_pairs.write_table_to_file(inp.net_table_file)
 
    # Prints basic stats for reconstructed network
    phylo_mut_pairs.print_stats(aln_pos = inp.pos_list,
                                 pos_subset = inp.pos_subset)


