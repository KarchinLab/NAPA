from collections import OrderedDict
import copy
import re

from napa.utils.serials import * 
from napa.utils.io import *

from napa.seq.parse import *
from napa.seq.format import *

class BioSeq(object):

    def __init__(self, seq_id = 'seqid', 
                 seq_str = '', 
                 seq_type = '', 
                 seq_pos_list = [],
                 pos_subset = [],
                 seq_annot = {'':['']}):

        self.seq_id = str(seq_id)
        self.seq_str = str(seq_str)
        self.length = len(self.seq_str)

        if len(seq_type): 
            self.seq_type = seq_type 
        else:
            self.determine_seq_type(seq_str)

        self.assign_standard_chars()        
        self.get_seq_pos(seq_pos_list, pos_subset)
        self.seq_annot = seq_annot


    def __repr__(self):
        return '%s\t%s\t%s\t%s' % (self.seq_id, self.seq_str, 
                                   self.seq_type, self.seq_annot)
    def copy(self):
        return copy.deepcopy(self)

    def add_annot(self, annot_dict):
        '''
        Add sequence annotation in the form of dictionary
        key == property or feature, 
        value == value of that feature for this sequence
        '''
        self.seq_annot.update(annot_dict)


    def get_seq_pos(self, seq_pos_list, pos_subset):
        '''
        Number alignment columns. 
        If no column numbers provided, start at 1.
        '''
        self.length = len(self.seq_str)

        if not len(seq_pos_list):
            self.seq_pos_list = range(1, len(self.seq_str)+1)

        elif len(seq_pos_list) == self.length:
            self.seq_pos_list = seq_pos_list

        elif len(seq_pos_list) < self.length:
            if len(pos_subset) == self.length:
                self.seq_pos_list = pos_subset
            elif len(pos_subset) > self.length:
                self.seq_pos_list = [pos_subset[i] for i in \
                                     range(self.length)]
            else:
                raise ValueError(' '.join(['Could not assign',
                                'sequence positions:', 
                                '\nSequence pos list length',
                                str(len(seq_pos_list)), 
                                '\ndoes not match alignment',
                                'length,', str(self.length), 
                                '\nor size of chosen subset', 
                                str(len(pos_subset))]))

        elif len(seq_pos_list) > self.length:
            if len(pos_subset) == self.length:
                self.seq_pos_list = pos_subset
            else:
                self.seq_pos_list = [seq_pos_list[i] \
                                     for i in range(self.length)]

        if len(pos_subset) \
                    and len(pos_subset) < len(self.seq_pos_list):
            self.extract_pos(pos_subset)



    def extract_pos(self, pos_subset):
        ''' 
        If pos_subset is not empty
        update sequence and position lists to that subset
        '''
        if len(pos_subset) > 0 and len(pos_subset) < self.length:
            seq_str = ''
            for posi, pos in enumerate(self.seq_pos_list):
                if pos in pos_subset:
                    if posi < len(self.seq_str):
                        seq_str += self.seq_str[posi]
            self.seq_str = seq_str
            self.seq_pos_list = copy.deepcopy(pos_subset)

        self.length = len(self.seq_str)


    def determine_seq_type(self, seq_str):
        '''
        Based on characters,
        find if sequence is DNA/RNA/Protein
        '''
        seq_char_set =  set(list(seq_str))
        if 'O' in seq_char_set or 'J' in seq_char_set:
            self.seq_type = 'unknown'
        elif seq_char_set <= set('ACDEFGHIKLMNPQRSTVWY.*-BZX'):
            self.seq_type = 'Protein'
        elif seq_char_set <= set('ATGC'+'RYKMSW'+'BDHVN.-'):
            self.seq_type = 'DNA'
        elif seq_char_set <= set('AUGC'+'IYRWSKM'+'BDHVN.-'):
            self.seq_type = 'RNA'
        else:
            self.seq_type = 'unknown'
        assert (self.seq_type in ['DNA','RNA','Protein']), \
            'Could not determine sequence type.' + \
            'Please enter standard Protein/DNA/RNA seq.'


    def assign_standard_chars(self):
        seq_type_to_chars = {'DNA':'ATGC', 'RNA':'AUGC', 
                             'Protein':'ACDEFGHIKLMNPQRSTVWY'}
        self.standard_chars = seq_type_to_chars[self.seq_type]


    def get_substitutions(self, other):
        '''
        Finds substitutions only (not insertions or deletions)
        between two sequences in alignment, considering 
        standard and unambiguous sequence characters only.
        (Needs proper error handling).
        '''

        if self.length != other.length or \
           self.seq_pos_list != other.seq_pos_list:
            stderr_write(['Need to align sequences',
                          'before comparing'])
            stderr_write([self.seq_id, 'length:', self.length, '\n',
                          other.seq_id, 'length:', other.length])
            raise ValueError('Trying to extract mutations' + \
                             'between unaligned sequences!') 

        if self.seq_type != other.seq_type:
            stderr_write(['Sequences have different',
                          'sequence types!\n',
                          'Please make sure ',
                          'they are both DNA/RNA/Protein.'])
            stderr_write(['Cannot compare sequences with IDs:', 
                          self.seq_id, 'and', other.seq_id,
                          'have types', self.seq_type, 'and',
                          other.seq_type, 'respectively.'])
            return 

        # only compare same sequence types and the 
        # standard, unambiguous characters
        mut_list = []
        for i in range(self.length):
            if self.seq_str[i] != other.seq_str[i]:
                if self.seq_str[i] in self.standard_chars and \
                   other.seq_str[i] in other.standard_chars:
                    pos = str(self.seq_pos_list[i])
                    mut_list.append(self.seq_str[i] + pos + \
                                    other.seq_str[i])

        return mut_list



class BioSeqAln(object):
    
    def __init__(self, **kwargs):
        self.seq_type = kwargs.get('seq_type', 'Protein')
        self.seq_annot = kwargs.get('seq_annot', 
                                    {'function':'default'}) 
        self.aln_pos = kwargs.get('aln_pos', [])
        self.orig_aln_pos = copy.deepcopy(self.aln_pos)
        self.pos_subset = kwargs.get('pos_subset', [])

        if len(kwargs.get('seqid_to_seq', '')):
            self.seqid_to_seq = \
                    copy.deepcopy(kwargs.get('seqid_to_seq', {}))
        
        elif len(kwargs.get('aln_fasta_file', '')):
            self.seqs_from_fasta_file(kwargs.get('aln_fasta_file', 
                                                 ''))

        
        if len(kwargs.get('annot_key', '')) and \
           (len(kwargs.get('seqid_to_annot_file', '')) or \
            len(kwargs.get('seqid_to_annot', {}))):
            self.annot_seqs(kwargs.get('annot_key', ''),
                            kwargs.get('seqid_to_annot_file', ''),
                            kwargs.get('seqid_to_annot', {}))
        

        if len(kwargs.get('sel_annot_key', '')) and \
            (kwargs.get('sel_annot_list', ['']) != [''] or \
            len(kwargs.get('sel_annot_file', ''))):
            self.subset_annot_seq_dict(\
                sel_annot_key = kwargs.get('sel_annot_key', ''),
                sel_annot_file = kwargs.get('sel_annot_file', ''),
                sel_annot_list = kwargs.get('sel_annot_list', ['']))


        if len(self.pos_subset):
            self.subset_aln_pos(self.pos_subset)

        if not len(self.aln_pos):
            stderr_write(['WARNING (BioSeqAln):',
                         'Default alignment positions set.'])
            self.aln_pos = self.seqid_to_seq[\
                            self.seqid_to_seq.keys()[0]].seq_pos_list

        self.depth = len(self.seqid_to_seq)
        self.length = len(self.aln_pos)

    def copy(self):
        return copy.deepcopy(self)

    def seqs_from_fasta_file(self, aln_fasta_file):
        self.seqid_to_seq = OrderedDict()
        self.add_sequences_from_file(aln_fasta_file)


    def add_sequences_from_file(self, aln_fasta_file):
        # To new sequences, assign 
        # positions in alignment before subsetting
        # (truncating columns)
        # Then apply the same subsetting
        aln_pos = self.orig_aln_pos

        fasta_recs = fasta_iter(aln_fasta_file)
    
        for r in fasta_recs:

            seq_id, seq_str = r[0], r[1]
            seq_id = re.sub('[!@#$.]|','', seq_id)
        
            if seq_id in self.seqid_to_seq:
                stderr_write(['Sequence ID', seq_id, 
                              'from FASTA file:\n',
                              aln_fasta_file,
                              'already in alignment.'
                              '\nSequence will not be added!'])
                continue

            self.seqid_to_seq[seq_id] = \
                    BioSeq(seq_id = seq_id, seq_str = seq_str,
                           seq_type = self.seq_type,
                           seq_pos_list = aln_pos,
                           pos_subset = self.pos_subset,
                           seq_annot = copy.deepcopy(self.seq_annot))



    def annot_seqs(self, annot_key = 'function', 
                   seqid_to_annot_file = '',
                   seqid_to_annot = {}):
        '''
        Add sequence function class annotation to 
        sequences in alignment.
        '''
        if len(seqid_to_annot_file):
            aln_seqid_to_annot = \
                            parse_keyval_dict(seqid_to_annot_file)
        elif len(seqid_to_annot):
            aln_seqid_to_annot = seqid_to_annot
        else:
            stderr_write(['WARNING (BioSeqAln):',
                          'No functions assigned'])
            return 

        for seqid in aln_seqid_to_annot:
            seqid_fasta = re.sub('[!@#$.]|','', seqid)
            if seqid_fasta in self.seqid_to_seq:
                self.seqid_to_seq[seqid_fasta].add_annot(\
                        {annot_key: aln_seqid_to_annot[seqid]})


    def subset_aln_pos(self, pos_subset = []):
        if not len(pos_subset):
            return # No subsetting of columns
        for seqid in self.seqid_to_seq:
            self.seqid_to_seq[seqid].extract_pos(pos_subset)
        self.aln_pos = pos_subset
        self.length = len(pos_subset)


    def seqids_with_annot(self, annot_key, annot_list):
        seqid_annot = []
        for seqid in self.seqid_to_seq:
            if annot_key not in self.seqid_to_seq[seqid].seq_annot:
                continue
            if self.seqid_to_seq[seqid].seq_annot[annot_key] \
               in annot_list:
                seqid_annot.append(seqid)
        return seqid_annot

                                                         
    def subset_annot_seq_dict(self,sel_annot_key = 'function', 
                              sel_annot_file = '', 
                              sel_annot_list = ['default']):
                                        
        '''
        Subset alignment sequences by a sequence id 
        using a sequence dictionary based on sequence annotation
        e.g. only sequences having a given function.
        '''
        if sel_annot_list != ['default'] and len(sel_annot_file):
            sel_annot_list = parse_column(sel_annot_file)

        seqid_to_seq_annot = OrderedDict()
        for seqid in self.seqids_with_annot(sel_annot_key, 
                                            sel_annot_list):
            seqid_to_seq_annot[seqid] = \
                            copy.deepcopy(self.seqid_to_seq[seqid])
        
        self.seqid_to_seq.clear()
        self.seqid_to_seq = copy.deepcopy(seqid_to_seq_annot)



    def get_seq_muts(self, wt_seq = ''):
        '''
        Obtain mutations between each alignment sequence and
        a reference sequence.
        '''
        self.seqid_to_mut = OrderedDict()
        for seqid in self.seqid_to_seq:
            self.seqid_to_mut[seqid] = \
                wt_seq.get_substitutions(self.seqid_to_seq[seqid])


    def get_seq_muts_id(self, wt_id = ''):
        '''
        Obtain set of mutations in each sequence 
        when a WT/reference sequence ID is given.
        '''       
        try:
            assert wt_id in self.seqid_to_seq
            wt_seq = self.seqid_to_seq[wt_id]
            self.get_seq_muts(wt_seq)
        
        except AssertionError:
            stderr_write(['BioSeqAln ERROR: ',
                          'Sequence id not in alignment. ',
                          'Cannot extract mutations!'])



    def get_seq_muts_id_seq(self, wt_id = '', wt_seq_str = ''):
        '''
        Obtain set of mutations in each sequence when a WT 
        sequence or sequence ID is given.
        '''
        if len(wt_id) and len(wt_seq_str):
            if wt_id in self.seqid_to_seq:
                if self.seqid_to_seq[wt_id].seq_str == wt_seq_str:
                    self.get_muts_id(wt_id)
                else:
                    stderr_write(['BioSeqAln WARNING:',
                                  'WT present in alignment but',
                                  'with different id.'])
                    wt_seq = BioSeq(wt_id, wt_seq_str, self.aln_pos)
                    self.get_seq_muts(wt_seq)
            else:
                wt_seq = BioSeq(wt_id, wt_seq_str, self.aln_pos)
                self.get_seq_muts(wt_seq)

        elif not len(wt_id) and len(wt_seq_str):
            wt_id = 'WildType'
            wt_seq = BioSeq(wt_id, wt_seq_str, self.aln_pos)
            self.get_seq_muts(wt_seq)

        elif len(wt_id) and not len(wt_seq_str):
            self.get_seq_muts_id(wt_id)


    def __repr__(self):
        out_str = ''
        for seqid in self.seqid_to_seq:
            out_str += write_wrap_fasta_seq(seqid + '___' + \
                '__'.join(['_'.join(ann_pair) for ann_pair in \
                        self.seqid_to_seq[seqid].seq_annot.items() \
                           if '' not in ann_pair and \
                           'default' not in ann_pair]), 
                        self.seqid_to_seq[seqid].seq_str)
        return out_str



    def write_muts(self, out_file, pos_subset):
        with open(out_file, 'wb') as of:
            for seqid in self.seqid_to_mut:
                if len(pos_subset):
                    of.write('\t'.join(seqid, 
                    ','.join([mut \
                              for mut in self.seqid_to_mut[seqid] \
                              if mut in pos_subset]))+'\n')
                else:
                    of.write('\t'.join(seqid, 
                            ','.join(self.seqid_to_mut[seqid]))+'\n')

