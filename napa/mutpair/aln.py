#! /usr/bin/env python 

from collections import defaultdict
from scipy import stats

from napa.utils.serials import * 
from napa.utils.io import *

from napa.seq.bioseq import *

def add_seq_id_to_set(seq_id, curr_set):
    '''
    expands current set with a (str of) single 
    sequence id or or a (listlike of) multiple
    ids
    '''
    if not isinstance(seq_id, basestring):
        curr_set.update(set(seq_id))
    else:
        curr_set.add(seq_id)
    return curr_set

# Mutation pairs are not allowed to share
# same numerical position
list_mut_pairs = lambda mut_list: \
                 [pair for pair in list_pairs(mut_list) \
                  if get_int_substring(str(pair[0])) != \
                  get_int_substring(str(pair[1]))]

class AlnMut(object):

    def __init__(self, seqs=set(), mut_str = None):
        self.mut_str = mut_str
        self.seqs = add_seq_id_to_set(seqs, set())
    
    def add_seq(self,seq):
        self.seqs = add_seq_id_to_set(seq, self.seqs)

    def __repr__(self):
        return self.mut_str


class AlnMutPair(object):
    '''
    Pair of mutations from seq alignment.
    '''
    def __init__(self, seqs = set(), mut_pair = None):

        # Unique set of seqs containing mutation
        if not len(seqs):
            raise ValueError(\
                "No sequences or empty sequence IDs "+ \
                    "associated with mut pair.")
        else:
            if not isinstance(seqs, basestring):
                self.seqs = set(seqs)
            else:
                self.seqs = set([seqs])

        if not mut_pair:
            raise ValueError('No mutations given ' + \
                             'for mutation pair')
        else:
            # Mutations in mut_pair should be AlnMut instances
            self.mut_pair = mut_pair 

    def __repr__(self):
        return '%s\t%s'%(self.mut_pair[0], self.mut_pair[1])

    def add_seq(self,seq):
        self.seqs = add_seq_id_to_set(seq, self.seqs)

    def get_contingency_table(self, all_seqs):
        ''' 
        Contingency table of occurrence counts of
        mutations, occurring individually and together. 
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

        self.contingency_table = [[occur01, occur0], 
                                  [occur1, occur_other]]
        self.num_seqs = len(all_seqs)
        

    def get_jaccard_weight(self, min_co_occur = 2):
        ''' 
        OBSOLETE 
        Jaccard weight with 
        1. minimum requirement for coocurrence count
        2. if no minimum is set, introduce correction for pairs 
        of muts occurring only once, and only together in one seq.
        '''
        [[num11, num10], [num01, num00]] = self.contingency_table

        if num11 < min_co_occur:
            self.jaccard = 0.
            return

        # Penalizes for small co-ocurrence num11
        epsilon = 1./self.num_seqs

        self.jaccard = max(0.,
            (num11 - epsilon) / float(num11 + num01 + num10))


    def get_mod_jaccard_weight(self, min_co_occur = 2):
        '''
        Modified Jaccard score based on sequence occurrence counts
        corrected for low mutation co-occurrence counts
        '''
        [[num11, num10], [num01, num00]] = self.contingency_table
        if num11 < min_co_occur:
            self.jaccard = 0.
            return

        # Penalizes for small co-ocurrence num11
        epsilon = float(num00)/(num11 + num10 + num01 + num00)
        self.jaccard = max(0.,
            (num11 - epsilon) / float(num11 + num01 + num10))


    def get_fisher_pval_weight(self, min_co_occur = 2):
        ''' 
        Unlike Jaccard, the smaller the p-value, 
        the greater the association '''
        if self.contingency_table[0][0] < min_co_occur:

            [[num11, num10], [num01, num00]] = \
                                self.contingency_table

            self.odds_ratio_more = float(num11 * num00)
            self.odds_ratio_more /= num10 * num01
            self.odds_ratio_less = float(num10 * num01)
            self.odds_ratio_less /= num11 * num00
            self.pval_less = 1.
            self.pval_more = 1.
            return

        self.odds_ratio, self.pval_less = \
        stats.fisher_exact(self.contingency_table, 
                           alternative = 'less')

        self.odds_ratio, self.pval_more = \
        stats.fisher_exact(self.contingency_table, 
                           alternative = 'greater')



class AlnMutPairSet(object):
    '''
    Set of mut pairs and associated metrics from an MSA.
    '''
    def __init__(self, inp):
        self.__dict__.update(vars(inp))

        # Pick out sequences with selected function of interest
        # (Alignment may consist of sequences with different 
        # annotations)
        self.aln.subset_annot_seq_dict(\
            sel_annot_key = 'function', 
            sel_annot_list = self.sel_prot_func)

        self.aln.get_seq_muts(self.wt_seq)

        # Print mutations from WT for each sequence in aln
        self.print_aln_mut_list()

        # Get all possible mut pairs in alignment, 
        #given selected function
        self.get_aln_mut_pairs()


    def __repr__(self):
        out_str = ''
        for mp in self.mut_pair_to_obj:
            
            mut_pair = self.mut_pair_to_obj[mp]
    
            if 'jaccard' in self.method.lower(): 
                if mut_pair.jaccard < self.thresh:
                    continue

                out_str += \
                '%s\t%.6f\n' % (mut_pair, 
                                mut_pair.jaccard)

            elif 'fisher' in self.method.lower(): 
                if mut_pair.pval_more > self.thresh:
                    continue

                out_str += \
                '%s\t%.6f\n' % (mut_pair, 
                                mut_pair.pval_more)

        return out_str


    def get_aln_mut_pairs(self):
        ''' 
        Extract pairs of mutations occurring 
        in at least one seq.
        '''

        self.mut_str_to_obj = {} 
        self.mut_pair_to_obj = {}

        for seqid in self.aln.seqid_to_mut:
            
            for m in self.aln.seqid_to_mut[seqid]:
                if m not in self.mut_str_to_obj:
                    self.mut_str_to_obj[m] = \
                        AlnMut(seqs = seqid, mut_str = m)
                else:
                    self.mut_str_to_obj[m].add_seq(seqid)

            for mp in \
                list_mut_pairs(self.aln.seqid_to_mut[seqid]):
                if mp not in self.mut_pair_to_obj:
                    self.mut_pair_to_obj[mp] = \
                    AlnMutPair(seqs = seqid, 
                               mut_pair = \
                            (self.mut_str_to_obj[mp[0]],
                             self.mut_str_to_obj[mp[1]]))
                else:
                    self.mut_pair_to_obj[mp].add_seq(seqid)

        # get contingency tables and network weights 
        # for pairs of co-occurring mutations
        all_seqs = set(self.aln.seqid_to_mut.keys())

        for mp in self.mut_pair_to_obj:
           
            mut_pair = self.mut_pair_to_obj[mp]
            
            mut_pair.get_contingency_table(all_seqs = \
                                           all_seqs)
            
            if 'mod_jaccard' in self.method.lower():
                mut_pair.get_mod_jaccard_weight(\
                                self.min_co_occur)

            elif 'jaccard' in self.method.lower():
                mut_pair.get_jaccard_weight(\
                                self.min_co_occur)
            
            elif 'fisher' in self.method.lower():
                mut_pair.get_fisher_pval_weight(\
                                self.min_co_occur)
            
    def print_aln_mut_list(self):
        outfile = ''
        if hasattr(self, 'print_seq_muts'):
            if to_bool(self.print_seq_muts):
                if hasattr(self, 'aln'):
                    outfile = \
                    '.'.join(self.aln_fasta_file.split('.')[0:-1])
                    outfile = outfile + '.mutation-list.txt'
                    header = 'sequence_id\tmutations_from_'
                    header += header + self.wt_id.replace(' ','-')

        if len(outfile):            
            stderr_write(['Printing alignment mutations from', 
                          self.wt_id, 'to file', outfile])
            
            write_keyval_dlist(self.aln.seqid_to_mut,
                               outfile, header, ';')

    def write_network_to_file(self, file_path):
        '''
        Writes network to file with path file_path
        Columns: source_node  target_node  weight
        '''
        with open(file_path, 'wb') as f:
            f.write(str(self))      
                    

    def write_table_to_file(self, file_path):
        '''
        Prints a table of all co-occurring mutations.
        Columns: Each mutation in pair, number of co-occurrences,
        contingency table (for Fisher's exact) and Jaccard Index.
        '''
        with open(file_path, 'wb') as f:
            f.write('\t'.join(\
                ['Mut1', 'Mut2', 'Weight' + str(self.method),
                 'Contingency_Table', 'Co-occur_Count', 
                 'weight_' + str(self.method)]) + '\n')
                    
            for mp in self.mut_pair_to_obj:
                mut_pair = self.mut_pair_to_obj[mp]

                if 'jaccard' in self.method.lower():
                    weight = mut_pair.jaccard

                elif 'fisher' in self.method.lower():
                    weight = mut_pair.pval_more
                
                f.write('%s\t%.6f\t%s\t%d\t%s\n' % \
                        (mut_pair, weight, 
                         mut_pair.contingency_table,
                         len(mut_pair.seqs),
                         ';'.join(sorted(list(mut_pair.seqs)))))
     

    def print_stats(self, aln_pos = [], pos_subset = []):

        stderr_write(['\nALIGNMENT NETWORK PROPERTIES:'])

        length = max(len(aln_pos), self.aln.length)

        stderr_write(['Alignment depth (num. rows)',
                      str(self.aln.depth) + \
                      '\nAlignment length (num. columns):', 
                      length])

        stderr_write([len(self.aln.aln_pos), 
                      "positions considered for network."])

        stderr_write([len(self.mut_pair_to_obj), 
                    "nonzero weight mutation pairs in network.\n"])

