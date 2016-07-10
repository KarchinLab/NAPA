import argparse
import sys
from os import listdir
from os.path import isfile, join
import itertools
from collections import defaultdict

from NAPA.utils.general import *
from NAPA.mutpair.phylo import * # Netw. construct from alignment


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', dest = 'alnFASTA', required = True, 
                        help = 'Path of multiple sequence alignment fasta file.')

    parser.add_argument('-pl', dest = 'posList', required = True, 
                        help = 'File containing protein residue positions ' + \
                        'corresponding to each column in the alignment.')

    parser.add_argument('-d', dest = 'treeDistThresh', type = float, required = True, 
                        help = 'Maximum distance between mutations on the tree. ' + \
                        'Mutations further appart not linked in network.')


    parser.add_argument('-wt', dest = 'wtId', required = True,  
                        help = 'Id of wildtype or ancestral sequence, ' + \
                        'assuming it is already in alignment. Can also provide '+ \
                        'separate FASTA with wildtype sequence.')

    parser.add_argument('-ps', dest = 'posSubset', required = True, 
                        help = 'File containing subset of aln/protein positions' + \
                        'to be considered in network')

    parser.add_argument('-af', dest = 'protFunc', required = True, 
                        help = 'File containing known protein functional annot. ' + \
                        'of sequences in alignment.')

    parser.add_argument('-ft', dest = 'funcTransitions', required = True, 
                        help = 'String/File containing allowed transitions between ' + \
                        'protein functions in the evolution of a final function.')

    parser.add_argument('-sf', dest = 'selProtFuncs', required = True, 
                        help = 'Path of file containing the selected final function.')


    parser.add_argument('-tl', dest = 'treeFileList', required = True,
                        help = 'File with a list of tree files from the ensemble.')
        
    parser.add_argument('-is', dest = 'intSeqFileList', required = True,
                        help = 'File with a list of fasta files where the internal' + \
                        'node sequences are stored.')
        
    parser.add_argument('-if', dest = 'intFuncFileList', required = True,
                        help = 'File with a list of fasta files where the internal' + \
                        'node functions are stored.')

    parser.add_argument('-mt', dest = 'mutPairType', type = str, required = True,
                        help = 'The type of network: dir(ected) or undir(ected).')
    
    parser.add_argument('-pt', dest = 'pvalThresh', type = float, required = True,
                        help = 'The maximum p-value for Fisher\'s exact test.' + \
                        'Used to extract significantly associated mut. pairs.')
    

    parser.add_argument('-nf', dest = 'netFile', required = True, 
                        help = 'Outputfile for result network.')

    args = parser.parse_args()

    return args


class PhyloNetInput(object):
    '''
    Prepare input files for alignment-based network reconstruction.
    '''
    def __init__(self, args):

        # List of positions/column numbers in alignment
        # If not provided, positions start with 1,2..alignment_length
        self.pos_list = [int(rec) for rec in parse_column(args.posList)]

        # List of protein residue positions considered for network.
        # (remove residues from signaling peptide, or other domains not studied).
        self.pos_subset = [int(rec) for rec in parse_column(args.posSubset)]

        # Functional transitions between phylogeny sequences allowed
        # (evolution of specific function)
        self.prot_func_transitions = parse_keyval_dict(args.funcTransitions)

        # List of absolute paths for phylogenetic trees in ensemble
        self.tree_files = parse_column(args.treeFileList)

        # List of internal node sequence fasta files for each tree
        self.int_seq_files = parse_column(args.intSeqFileList)

        # List of internal node function files for each tree
        self.int_func_files = parse_column(args.intFuncFileList)

        # check we have required inputs for each tree
        if len(self.tree_files) != len(self.int_seq_files):
            stderr_write(['Unequal number of tree files and ',
                          'internal node files.',
                          'Using minimum number of files of both.'])

            min_trees = min(len(self.tree_files), len(self.int_seq_files))
            self.tree_files = self.tree_files[0:min_trees]
            self.int_seq_files = self.int_seq_files[0:min_trees]


def main():
    args = parse_args()
    inp = PhyloNetInput(args)
    for i, (tf, isf, iff) in enumerate(zip(inp.tree_files, 
                                           inp.int_seq_files, inp.int_func_files)):

        stderr_write(['\n','Extracting mutation pairs for tree file:\n',tf])
        # Set of mutation pairs occurring along a single tree
        phylo_mut_pair_set = PhyloMutPairSet(leaf_fasta_file = args.alnFASTA,
                                             int_node_fasta_file = isf,
                                             aln_pos = inp.pos_list,
                                             wt_id = args.wtId,
                                             pos_subset = inp.pos_subset,
                                             leaf_seqid_to_prot_func_file = args.protFunc,
                                             int_seqid_to_prot_func_file = iff,
                                             prot_func_transitions_file = args.funcTransitions,
                                             sel_prot_func_file = args.selProtFuncs,
                                             dist_thresh = args.treeDistThresh,
                                             tree_nwk_file = tf,
                                             mut_pair_type = args.mutPairType)
        # Combined sets of mutations from all trees 
        # The number of trees in ensemble with a significant right-tailed p-value
        # are the weights output for each mutation pair
        if i == 0:
            cmb_mut_pair_set = CombinedPhyloMutPairs(phylo_mut_pair_set, args.pvalThresh)
        else:
            cmb_mut_pair_set.add_mut_pairs(phylo_mut_pair_set)

        
    with open(args.netFile, 'wb') as f:
        f.write(str(cmb_mut_pair_set))

if __name__ == '__main__': main()        
