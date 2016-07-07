#! /usr/bin/env python 

from collections import defaultdict
from scipy import stats

from NAPA.utils.general import *
from NAPA.seq.BioSeq import *


class AlnMut(object):

    def __init__(self, seqs=set(), mut_str = None):
        self.mut_str = mut_str
        self.seqs = set(seqs)
    
    def add_seq(self,seq):
        self.seqs.add(seq) 
    
    def add_seqs(self,seqs):
        self.seqs.update(seqs)

    def __repr__(self):
        return self.mut_str


class AlnMutPair(object):
    '''
    Pair of muts from seq alignment.
    '''
    def __init__(self, seqs = set(), mut_pair = None):

        # Unique set of seqs containing mut
        if not len(seqs):
            raise ValueError("No seqs associated with mut pair.")
        else:
            self.seqs = set(seqs)

        if not mut_pair:
            raise ValueError("No muts given for mut pair")
        else:
            # Mutations in mut_pair should be AlnMut instances
            self.mut_pair = mut_pair 

    def __repr__(self):
        return '%s\t%s'%(str(self.mut_pair[0]), str(self.mut_pair[1])) 

    def add_seq(self, seq):
        self.seqs.update(seq) 


    def get_contingency_table(self, all_seqs):
        ''' Contingency table of occurrence counts of muts
        individually or together. 
        Occurrence count = number of seqs containing mut.
        '''

        occur01 = len(self.seqs) #occurrence count of both muts

        # Number of seqs containing first or second mut, not both 
        occur0 = len(self.mut_pair[0].seqs - self.seqs) 
        occur1 = len(self.mut_pair[1].seqs - self.seqs)

        # Number of seqs containing neither mut
        occur_other = len(all_seqs \
                        - self.mut_pair[0].seqs \
                        - self.mut_pair[1].seqs)

        self.contingency_table = [[occur01, occur0], [occur1, occur_other]]
        self.num_seqs = len(all_seqs)


    def get_jaccard_weight(self, min_co_occur = 0):
        ''' Jaccard weight with 
        1. minimum requirement for coocurrence count
        2. if no minimum is set, introduce correction for pairs 
        of muts occurring only once, and only together in one seq.'''
        [[occur01, occur0], [occur1, occur_other]] = self.contingency_table
        if min_co_occur == 0:
            self.jaccard = (float(occur01) - 1./self.num_seqs) / \
                            (float(occur01 + occur0 + occur1))            
        elif occur01 < min_co_occur:
            self.jaccard = 0.
        else:
            self.jaccard = (float(occur01) - 1./self.num_seqs) / \
                            (float(occur01 + occur0 + occur1))


    def get_fisher_pval_weight(self):
        ''' Unlike Jaccard, the smaller the p-value, 
        the greater the association '''
        self.odds_ratio_less, self.p_val_less = \
                        stats.fisher_exact(self.contingency_table, 
                                           alternative = 'less')

        self.odds_ratio_more, self.p_val_more = \
                        stats.fisher_exact(self.contingency_table, 
                                           alternative = 'greater')



class AlnMutPairSet(object):
    '''
    Set of mut pairs and associated metrics from an MSA.
    '''
    def __init__(self, aln_fasta_file = '', aln_pos = [], 
                 wt_seq_str = '', wt_id = '', pos_subset = [],
                 seqid_to_prot_func = {'':''}, 
                 seqid_to_prot_func_file = '',
                 sel_prot_func = [''], 
                 sel_prot_func_file = ''):

        # Alignment from input fasta file
        # Do not select sequences by function yet
        # (that would remove the WT/reference sequence from set
        self.aln = BioSeqAln(aln_fasta_file = aln_fasta_file,
                             aln_pos = aln_pos, 
                             pos_subset = pos_subset,
                             annot_key = 'function',
                             seqid_to_annot = seqid_to_prot_func, 
                             seqid_to_annot_file = seqid_to_prot_func_file,
                             sel_annot_key = '', 
                             sel_annot_list = [''],
                             sel_annot_file = '')

        # Get mutations from ancestral/wt sequence 
        # for each sequence in alignment
        # option to focus on subset of aln. positions
        self.get_wt_seq(wt_id, wt_seq_str)
        self.wt_seq.extract_pos(pos_subset)


        #Pick out sequences with selected function of interest
        #(Alignment may consist of sequences with different annotations)
        self.aln.subset_annot_seq_dict(sel_annot_key = 'function',
                                       sel_annot_list = sel_prot_func,
                                       sel_annot_file = sel_prot_func_file)

        self.aln.get_seq_muts(self.wt_seq)
        # stderr_write([self.aln_func.seqid_to_mut])

        # Get all possible mut pairs in alignment, 
        #given selected function
        self.get_aln_mut_pairs()
        self.print_stats()

    
    def get_wt_seq(self, wt_id, wt_seq_str):
        '''
        Obtain WT/ancestral sequence to which all other sequences in 
        the alignment should be compared. Can provide id+sequence,
        or the id of a sequence already in the alignment.
        '''
        if len(wt_id) and len(wt_seq_str) == self.aln.length:
            self.wt_seq = BioSeq(seq_id = wt_id, seq_str = wt_seq_str,
                                 seq_type = 'Protein', 
                                 seq_pos_list = self.aln_pos)
        elif len(wt_id) and len(wt_seq_str) != self.aln.length:
            if wt_id in self.aln.seqid_to_seq:
                stderr_write(['AlnMutPairSet WARNING:',
                              'Sequence length does not match alignment,',
                              'but sequence id in alignment. Replacing',
                              'with corresp. alignment sequence.'])
                self.wt_seq = self.aln.seqid_to_seq[wt_id]
            else:
                stderr_write(['Invalid  input for WT sequence.',
                              'Please enter valid id already in alignment.',
                               'or a new sequence.'])
                              
        elif not len(wt_id) and len(wt_seq_str) == self.aln.length:
            self.wt_seq = BioSeq(seq_id = 'Wild_type', 
                                 seq_str = wt_seq_str,
                                 seq_type = 'Protein', 
                                 seq_pos_list = aln_pos)

        else:
            stderr_write(['Invalid  input for WT sequence.',
                          'Please enter valid id already in alignment.',
                          'or a new sequence.'])
        
    def get_aln_mut_pairs(self):
        ''' Extract pairs of mutations occurring in at least one seq.'''

        self.mut_str_to_obj, self.mut_pair_to_obj = {}, {}
        for seqid in self.aln.seqid_to_mut:
            for mut_str in self.aln.seqid_to_mut[seqid]:
                if mut_str not in self.mut_str_to_obj:
                    self.mut_str_to_obj[mut_str] = \
                            AlnMut(seqs = set([seqid]), 
                                   mut_str = mut_str)
                else:
                    self.mut_str_to_obj[mut_str].add_seq(seqid)

            for mut_str_pair in list_pairs(self.aln.seqid_to_mut[seqid]):
                if mut_str_pair not in self.mut_pair_to_obj:
                    self.mut_pair_to_obj[mut_str_pair] = \
                        AlnMutPair(seqs = set([seqid]), 
                                   mut_pair = (self.mut_str_to_obj[mut_str_pair[0]],
                                               self.mut_str_to_obj[mut_str_pair[1]]))
                else:
                    self.mut_pair_to_obj[mut_str_pair].add_seq(seqid)

        # get contingency table for mutation pair stats
        for mut_pair in self.mut_pair_to_obj:
            self.mut_pair_to_obj[mut_pair].get_contingency_table(all_seqs = \
                                           set(self.aln.seqid_to_mut.keys()))
                        


    def write_jaccard_weights_table(self, min_co_occur = 0):
        out_str = '\t'.join(['_Mut1_', '_Mut2_', 'Num_sequences_co-occur', 
                             'Contingency_Tab(common_1_2_other)',
                             'Jaccard_index']) + '\n'
        

        for mut_str_pair in self.mut_pair_to_obj:
            aln_mut_pair = self.mut_pair_to_obj[mut_str_pair]
            aln_mut_pair.get_jaccard_weight(min_co_occur = min_co_occur)

            if aln_mut_pair.jaccard == 0.:
                continue

            to_print = [str(aln_mut_pair)]
            to_print.append(str(len(aln_mut_pair.seqs)))
            to_print.append(str(aln_mut_pair.contingency_table))
            to_print.append(str(aln_mut_pair.jaccard))

            out_str +=  '\t'.join(to_print) + '\n'

        return out_str
                                         


    def __repr__(self):
        '''
        Writes table of mutation pairs and statistics.
        '''
        out_str = '\t'.join(['_Mut1', '_Mut2', 'numSharedSeqs', 
                             'contTab(common_1_2_other)']) + '\n'
        for aln_mut_pair in self.aln_mut_pairs:
            to_print = [str(aln_mut_pair)]
            to_print.append(str(len(aln_mut_pair.seqs)))
            to_print.append(str(aln_mut_pair.contingency_table))
            out_str +=  '\t'.join(to_print) + '\n'

        return out_str


    def print_stats(self):

        stderr_write(['Alignment depth and length:', self.aln.depth, self.aln.length])

        stderr_write(["Total positions considered:", len(self.aln.aln_pos)])
                     
        stderr_write(["Number of functional mutation pairs", 
                      len(self.mut_pair_to_obj)])


