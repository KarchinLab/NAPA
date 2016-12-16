import os

# General parsing/container manipulation
from napa.utils.serials import * 
from napa.utils.io import *

# Configuration and input checking
from napa.utils.config import Config

# Bio sequence manipulation
from napa.seq.parse import *
from napa.seq.bioseq import *

# Network construction from alignment
from napa.mutpair.aln import *


class AlnNetInput(object):
    '''
    Prepare all inputs for alignment-based 
    network reconstruction.
    '''
    def __init__(self, config):
        self.__dict__.update(vars(config))

        print
        print '\n'.join(sorted(['%s\t%s'%(k,v) for (k,v) \
                                in vars(self).iteritems()]))

        # REQUIRED: Reads wild-type sequence or
        # provides the ID of this sequence in the alignment
        if os.path.isfile(self.wild_type):
            self.wt_dict = fasta_to_dict(self.wild_type)
            if len(self.wt_dict):
                self.wt_id = self.wt_dict.keys()[0]
                self.wt_seq_str = self.wt_dict[self.wt_id]
        else:
            self.wt_id = self.wild_type
            self.wt_seq_str = ''

        # OPTIONAL: Allows for custom numbering of positions
        # pos_list length must match alignment width
        if len(self.pos_list_file) == 0:
            self.pos_list = []
        else:
            self.pos_list = \
            [int(p) for p in \
             parse_column(self.pos_list_file)]

        # OPTIONAL: List of protein residue positions 
        # considered for network.
        # Allows removal of residues from signaling peptide, 
        # other domains/regions not studied.
        if len(self.pos_subset_file) == 0:
            self.pos_subset = []
        else:
            self.pos_subset = \
            [get_int_substring(p) for p in \
             parse_column(self.pos_subset_file)]


        # OPTIONAL: Ranges of protein residues to include
        # Can be used instead of pos_subset_file
        if len(self.protein_ranges):
            if len(self.pos_subset):
                stderr_write(['WARNING:Position subset already',
                              'defined in file',
                              self.pos_subset_file,
                              'Protein range list ignored.'])
            else:
                self.pos_subset = set()
                stderr_write([str(self.protein_ranges)])
                for range_str in self.protein_ranges:
                    recs = range_str.strip().split('-')
                    if len(recs) != 2:
                        raise ValueError('Invalid protein' + \
                                'range specified in config:' + \
                                str(recs))
                    begin = get_int_substring(recs[0])
                    end = get_int_substring(recs[1])
                    self.pos_subset |= set(range(begin, 
                                             end + 1))
                self.pos_subset = sorted(list(self.pos_subset))
                
        # OPTIONAL: Threshold for including links in network
        # For p-values from Fisher's Exact Test, this is the 
        # maximum p-value; 
        # For Jaccard indices, this is the minimum Jaccard 
        # index weight
        if not hasattr(self, 'thresh'):
            self.thresh = 0.05
            stderr_write(['WARNING: No threshold provided',
                          'for network method', self.method,
                          '\nDefault set to:', self.thresh])
        elif not to_bool(self.thresh) and not \
             self.thresh > 0: 
            self.thresh = 0.05
            stderr_write(['WARNING: No threshold provided',
                          'for network method', self.method,
                          '\nDefault set to:', self.thresh])

        if not hasattr(self, 'min_co_occur'):
            self.min_co_occur = 1
            stderr_write(['WARNING: No minimum co-occurrence',
                          'count provided for network method', 
                          self.method,
                          '\nDefault set to:', self.min_co_occur])

        elif not to_bool(self.min_co_occur):
            self.min_co_occur = 1
            stderr_write(['WARNING: No minimum co-occurrence',
                          'count provided for network method', 
                          self.method,
                          '\nDefault set to:', self.min_co_occur])

        # OPTIONAL: Selection by functional annotation
        # Prot. func. for each seq. in aln with known func. annot.
        if not hasattr(self, 'prot_func_file'):
            self.seqid_to_prot_func = {'function':'default'}
        elif not to_bool(self.prot_func_file):
            self.seqid_to_prot_func = {'function':'default'}
        else:
            self.seqid_to_prot_func = \
                parse_keyval_dict(self.prot_func_file)

        # OPTIONAL: Subset of protein selected functions 
        # considered 
        if not hasattr(self, 'sel_prot_func_file'):
            self.sel_prot_func = ['default']
        elif not to_bool(self.sel_prot_func_file):
            self.sel_prot_func = ['default']
        else:
            self.sel_prot_func = \
            parse_column(self.sel_prot_func_file)

        # Get protein sequence alignment object
        self.get_alignment()
        
        # Get wt sequence id and string
        self.get_wt_seq()




    def print_input_summary(self):
        stderr_write(['Building network from alignment.'])
        stderr_write(['Alignment used:\n', 
                      self.aln_fasta_file])
        

    def get_alignment(self):
        '''
        Alignment from input fasta file
        Do not select sequences by function yet
        (that would remove the WT/reference sequence from set
        '''
        self.aln = \
            BioSeqAln(aln_fasta_file = self.aln_fasta_file,
                      aln_pos = self.pos_list, 
                      pos_subset = self.pos_subset,
                      annot_key = 'function',
                      seqid_to_annot = self.seqid_to_prot_func, 
                      sel_annot_list = self.sel_prot_func)


    def get_wt_seq(self):
        '''
        Obtain WT/ancestral sequence to which all other 
        sequences in the alignment should be compared. 
        Can provide id+sequence, or the id of a sequence 
        already in the alignment.
        If sequence is not in alignment, sequence 
        still needs to match in positions (be aligned) 
        to the multiple sequence alignment.
        No checking or alignment is done.
        '''
        self.wt_id = self.wt_id if len(self.wt_id.strip()) \
                else 'Wild_type'
        
        if len(self.wt_seq_str) == self.aln.length:
            self.wt_seq = \
                BioSeq(seq_id = self.wt_id, 
                       seq_str = self.wt_seq_str,
                       seq_type = 'Protein', 
                       seq_pos_list = self.aln_pos)
        else:
            #make sure wt sequence is trimmed to 
            #the subset of positions (if any)
            #before rejecting due to length mismatch
            if len(self.pos_subset) and \
               self.aln.length == len(self.pos_subset) and \
               len(self.wt_seq_str) >= len(self.pos_subset):
                self.wt_seq = \
                    BioSeq(seq_id = self.wt_id, 
                           seq_str = self.wt_seq_str,
                           seq_type = 'Protein', 
                           seq_pos_list = self.aln.aln_pos, 
                           pos_subset = self.pos_subset)
            
            elif self.wt_id in self.aln.seqid_to_seq:
                if self.wt_seq_str.strip() != '':
                    stderr_write([\
                        'WARNING (AlnMutPairSet):',
                        'Input wild type sequence length',
                        len(self.wt_seq_str), 'does not match', 
                        'alignment length,', self.aln.length,
                        '\nbut sequence id is in alignment.', 
                        '\nReplacing with corresponding',
                        'alignment sequence.'])
                self.wt_seq = self.aln.seqid_to_seq[self.wt_id]
            else:
                raise ValueError(\
                    ' '.join(['ERROR: Invalid  input for',
                              'WT sequence.',
                              '\nPlease enter a valid ID',
                              'already in alignment, '
                              '\n OR a new ID with sequence', 
                              'that aligns to the alignment',
                              '(or aln. column subset chosen).\n',
                              self.wt_id, 'not in\n', 
                              str(sorted(self.aln.seqid_to_seq))]))



def run_aln_mut_pairs(config):
    '''
    Processes network input information.
    Generates and writes network to text file(s).
    '''
    # Obtains alignment, wt sequence, parameters 
    # for alignment based network reconstruction
    inp = AlnNetInput(config)

    # Calculates edge weights for co-occuring mutation pairs
    aln_mut_pair_set = AlnMutPairSet(inp)



    # Output network in tab-delimited format
    # Source\tTarget\tWeight
    stderr_write(['\nWriting network to file:\n' + inp.net_file]) 
    aln_mut_pair_set.write_network_to_file(inp.net_file)

    # outputs additional information about edge weights
    # to a separate table
    if to_bool(inp.output_edge_weight_table):
        stderr_write(['\nWriting network extended network',
                      'table to file:\n' + inp.net_table_file]) 
        aln_mut_pair_set.write_table_to_file(inp.net_table_file)
 
    # Prints basic stats for reconstructed network
    aln_mut_pair_set.print_stats(aln_pos = inp.pos_list,
                                 pos_subset = inp.pos_subset)
