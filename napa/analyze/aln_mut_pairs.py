import argparse
import sys
import itertools
import sys

from napa.seq.parse import *
from napa.utils.general import * #General parsing
from napa.mutpair.aln import * # Network construct from alignment


def parse_args():
        parser = argparse.ArgumentParser()

        parser.add_argument('-a', dest = 'alnFASTA', required = True, 
                            help = 'Path of multiple sequence alignment' +  
                            'fasta file.')
        
        parser.add_argument('-wt', dest = 'wtId', required = True,  
                            help = 'Id of wildtype or ancestral seq., ' + \
                            'assuming it is already in alignment. ' + \
                            'Can also provide a separate FASTA path ' + \
                            'with wildtype sequence.')

        parser.add_argument('-p', dest = 'posList', required = False, 
                            default = '',
                            help = 'OPTIONAL: File containing protein' + \
                            'residue positions corresponding to each ' + \
                            'column in the alignment.')

        parser.add_argument('-ps', dest = 'posSubset', required = False, 
                            default = '',
                            help = 'OPTIONAL: File containing '+ \
                            'protein positions' + \
                            'to be considered in network. ' + \
                            'By default all positions are considered.')

        parser.add_argument('-pf', dest = 'protFunc', required = False, 
                            default = '',
                            help = 'OPTIONAL: File containing ' + 
                            'known protein ' + \
                            'functional annotation of ' + \
                            'sequences in alignment.')

        parser.add_argument('-sf', dest = 'selProtFunc', required = False, 
                            default = '',
                            help = 'OPTIONAL: String/File containing ' + \
                            'protein funcstions' + \
                            'to be included in network.')

        parser.add_argument('-nf', dest = 'netFile', required = True, 
                            help = 'Outputfile for result network.')

        parser.add_argument('-tf', dest = 'tableFile', required = False,
                            default = '',
                            help = 'Outputfile for weights table.')

        args = parser.parse_args()


        return args

class AlnNetInput(object):
        '''
        Prepare input files for alignment-based network reconstruction.
        '''
        def __init__(self, args):
                '''
                List of positions/column numbers in alignment
                If not provided, positions start with 1,2..alignment_length
                '''
                if os.path.isfile(args.wtId):
                        self.wt_dict = fasta_to_dict(args.wtId)
                        if len(self.wt_dict):
                                self.wt_id = self.wt_dict.keys()[0]
                                self.wt_seq = self.wt_dict[self.wt_id]
                else:
                        self.wt_id = args.wtId
                        self.wt_seq = ''

                if len(args.posList) == 0:
                        self.pos_list = []
                else:
                        self.pos_list = [int(rec) for rec \
                                         in parse_column(args.posList)]


                # List of protein residue positions considered for network.
                # Allows removal of residues from signaling peptide, 
                # other domains/regions not studied.
                if len(args.posSubset) == 0:
                        self.pos_subset = []
                else:
                        self.pos_subset = [int(rec) for rec \
                                           in parse_column(args.posSubset)]


                # Prot. func. for each seq. in aln with known func. annot.
                if len(args.protFunc) == 0:
                        self.seqid_to_prot_func = {'function':'default'}
                else:
                        self.seqid_to_prot_func = \
                                        parse_keyval_dict(args.protFunc)

                # subset of protein selected functions considered 
                if len(args.selProtFunc) == 0:
                        self.sel_prot_func_list = ['default']
                else:
                        self.sel_prot_func_list = \
                                        parse_column(args.selProtFunc)


def main():
    args = parse_args()

    # Construct network based on alignment
    inp = AlnNetInput(args)
    alnMutPairSet = AlnMutPairSet(aln_fasta_file = args.alnFASTA,
                                  aln_pos = inp.pos_list,
                                  wt_id = inp.wt_id,
                                  wt_seq_str = inp.wt_seq, 
                                  pos_subset = inp.pos_subset, 
                                  seqid_to_prot_func = \
                                  inp.seqid_to_prot_func,
                                  sel_prot_func = inp.sel_prot_func_list)
        
    with open(args.netFile, 'wb') as f:
        f.write(alnMutPairSet.write_jaccard_weights_network(\
                                                min_co_occur = 2))
    

    if len(args.tableFile):
        with open(args.tableFile, 'wb') as f:
                f.write(alnMutPairSet.write_jaccard_weights_table(\
                                                min_co_occur = 2))
        
    
if __name__ == '__main__': main()


