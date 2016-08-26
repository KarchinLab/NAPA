import argparse
import sys
import itertools
import sys

#print sys.path

from napa.utils.general import *
from napa.mutpair.aln import * # Netw. construct from alignment


def parse_args():
        parser = argparse.ArgumentParser()

        parser.add_argument('-a', dest = 'alnFASTA', required = True, 
                            help = 'Path of multiple sequence alignment fasta file.')

        parser.add_argument('-p', dest = 'posList', required = True, 
                            help = 'File containing protein residue positions ' + \
                            'corresponding to each column in the alignment.')

        parser.add_argument('-wt', dest = 'wtId', required = True,  
                            help = 'Id of wildtype or ancestral sequence, ' + \
                            'assuming it is already in alignment. Can also provide '+ \
                            'separate FASTA with wildtype sequence.')

        parser.add_argument('-ps', dest = 'posSubset', required = True, 
                            help = 'File containing protein positions' + \
                            'to be considered in network')

        parser.add_argument('-pf', dest = 'protFunc', required = True, 
                            help = 'File containing known protein functional annot. ' + \
                            'of sequences in alignment.')

        parser.add_argument('-sf', dest = 'selProtFunc', required = True, 
                            help = 'String/File containing protein funcs' + \
                            'to be included in net.')

        parser.add_argument('-nf', dest = 'netFile', required = True, 
                            help = 'Outputfile for result network.')

        args = parser.parse_args()

        return args

class AlnNetInput(object):
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


                # prot. func. for each seq. in aln with known func. annot.
                self.seqid_to_prot_func = parse_keyval_dict(args.protFunc)

                # subset of protein selected functions considered 
                self.sel_prot_func_list = parse_column(args.selProtFunc)


def main():
    args = parse_args()

    # Construct network based on alignment
    inp = AlnNetInput(args)
    alnMutPairSet = AlnMutPairSet(aln_fasta_file = args.alnFASTA,
                                  aln_pos = inp.pos_list,
                                  wt_id = args.wtId,
                                  pos_subset = inp.pos_subset, 
                                  seqid_to_prot_func = inp.seqid_to_prot_func,
                                  sel_prot_func = inp.sel_prot_func_list)
        
    with open(args.netFile, 'wb') as f:
        f.write(alnMutPairSet.write_jaccard_weights_network(min_co_occur = 2))
    
if __name__ == '__main__': main()


