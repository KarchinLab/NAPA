from NAPA.utils.general import *
from NAPA.seq.parse import *

class AlnMut(object):
    def __init__(self, seqs=set(), mut_str = None):
        self.mut_str = mut_str
        self.seqs = set(seqs)
    
    def add_sequence(self,sequence):
        self.seqs.add(sequence) 
    
    def add_seqs(self,seqs):
        self.seqs.update(seqs)
    def __repr__(self):
        return self.mut_str


class AlnMutPair(object):
    '''
    Pair of muts from sequence alignment.
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

    def add_sequence(self, sequence):
        self.seqs.update(sequence) 


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
        of muts occurring only once, and only together in one sequence.'''
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
    ''' Set of mut pairs and associated metrics from an MSA. '''

    def __init__(self, fasta_file_list = [], pos_list = [], 
                 pos_to_wt = {}, pos_subset = [],
                 seqid_to_prot_func = {}, 
                 prot_func_list =[]):

        self.fasta_file_list = fasta_file_list
        self.pos_list = pos_list
        self.pos_to_wt = pos_to_wt
        self.pos_subset = pos_list if not len(pos_subset) else pos_subset
        self.seqid_to_prot_func = seqid_to_prot_func.copy()
        self.prot_func_list = prot_func_list #list of selected functions

        # Get seqs having a protein function of interest
        self.get_prot_func_seqs()
        # Get muts within seqs
        self.get_prot_func_muts()

        # Get all possible mut pairs in alignment, 
        #given selected function
        self.aln_mut_pairs = []
        self.get_aln_mut_pairs()


    def __repr__(self):
        ''' Writes table of mutation pairs and statistics '''

        out_str = '\t'.join(['_Mut1', '_Mut2', 'numSharedSeqs', 
                             'contTab(common_1_2_other)']) + '\n'
        
        for aln_mut_pair in self.aln_mut_pairs:
            to_print = [str(aln_mut_pair)]
            to_print.append(str(len(aln_mut_pair.seqs)))
            to_print.append(str(aln_mut_pair.contingency_table))
            out_str +=  '\t'.join(to_print) + '\n'

        return out_str


    def get_prot_func_seqs(self):
        '''Extract seqs based on associated prot. function'''

        # read in fasta seqs as {id:seq} dictionary
        self.seqid_to_sequence = {}
        for ff in self.fasta_file_list:
            seqid_to_sequence = fasta_to_dict(ff)
            for seqid in seqid_to_sequence:
                if not len(prot_func_list): #any function allowed by default
                    self.seqid_to_sequence[seqid] = seqid_to_sequence[seqid]
                else:
                    if seqid_to_prot_func[seqid] not in self.prot_func_list: 
                        continue
                    self.seqid_to_sequence[seqid] = seqid_to_sequence[seqid]

    
    def get_prot_func_muts(self):
        '''Extract muts in seqs with prot. function of interest.'''
        
        standard_aa = "ACDEFGHIKLMNPQRSTVWY" # list of standard amino acids
        # Always need reference sequence / common ancestor
        anc_seq = [self.pos_to_wt[p] for p in self.pos_list \
                   if p in self.pos_to_wt]

        self.seqid_to_muts = {}
        for seqid in self.seqid_to_sequence:
            seq = self.seqid_to_sequence[seqid]
            self.seqid_to_muts[seqid] = []
            for i in range(len(seq)):
                if seqid_seq[i] not in standard_aa: continue
                if anc_seq[i]! = seqid_seq[i]:
                    pos = self.pos_list[i]
                    if pos not in self.pos_subset: continue
                    self.seqid_to_muts[seqid].append(anc_seq[i] + \
                                                     pos + seqid_seq[i])


    def get_aln_mut_pairs(self):
        ''' Extract pairs of mutations occurring in at least one seq.'''
        mut_pair_to_seqs = defaultdict(list) 
        mut_to_seqs = defaultdict(list)
        all_seqs = set()
        for seqid in self.seqid_to_muts:
                seqid_mut_pairs = list_pairs(self.seqid_to_muts[seqid])
                for mut_pair in seqid_mut_pairs:
                    mut_pair_to_seqs[mut_pair].append(seqid)
                    mut_to_seqs[mut_pair[0]].append(seqid)
                    mut_to_seqs[mut_pair[1]].append(seqid)
                    all_seqs.update(seqid)

        all_muts = set(mut_to_seqs.keys()) 
    
        self.mut_to_obj = {}

        for mut_pair in mut_pair_to_seqs:
            if mut_pair[0] not in self.mut_to_obj:
                self.mut_to_obj[mut_pair[0]] = AlnMut(mut_str = mut_pair[0])
            if mut_pair[1] not in self.mut_to_obj:
                self.mut_to_obj[mut_pair[1]] = AlnMut(mut_str = mut_pair[1])

            self.mut_to_obj[mut_pair[0]].add_seqs(set(mut_pair_to_seqs[mut_pair]))
            self.mut_to_obj[mut_pair[1]].add_seqs(set(mut_pair_to_seqs[mut_pair]))

            aln_mut_pair = AlnMutPair(seqs = list(set(mut_pair_to_seqs[mut_pair])), 
                                      mut_pair = (self.mut_to_obj[mut_pair[0]], 
                                                  self.mut_to_obj[mut_pair[1]]))
            self.aln_mut_pairs.append(aln_mut_pair)

        # get contingency table for mutation pair stats
        for mut_pair in self.aln_mut_pairs:
            mut_pair.get_contingency_table(all_seqs)
                        


    def jaccard_table(self, min_co_occur = 0):
        out_str = '\t'.join(['_Mut1_', '_Mut2_', 'numMuts', 
                             'contTable(common_1_2_other)','Jaccard']) + '\n'
        

        for aln_mut_pair in self.aln_mut_pairs:
            aln_mut_pair.get_jaccard_weight(min_co_occur)

            to_print = [str(aln_mut_pair)]
            to_print.append(str(len(aln_mut_pair.seqs)))
            to_print.append(str(aln_mut_pair.contingency_table))
            to_print.append(str(aln_mut_pair.jaccard))
            out_str +=  '\t'.join(toPrint) + '\n'

        return out_str
                                         



