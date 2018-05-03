# General dependencies
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import pandas as pd
import sys

# General parsing/container manipulation
from napa.utils.serials import * 
from napa.utils.io import *

## Bio sequence manipulation
#from napa.seq.parse import *
#from napa.seq.bioseq import *

# Network construction from alignment
from napa.pospair.align import AlnPosPairSet


class PosNetInput(object):
    '''
    Prepare all inputs for alignment-based 
    position network reconstruction.
    '''
    def __init__(self, config):
        #inherits same attributes from Config class
        self.__dict__.update(vars(config)) 
        # OPTIONAL: Threshold for including links in network
        # For p-values from Fisher's Exact Test, this is the 
        # maximum p-value; 
        # For Jaccard indices, this is the minimum Jaccard 
        # index weight
        if not hasattr(self, 'thresh'):
            self.thresh = 1
            stderr_write(['WARNING: No threshold provided',
                          'for network method', self.method,
                          '\nDefault set to:', self.thresh])
        elif not to_bool(self.thresh) and not \
             self.thresh > 0: 
            self.thresh = 1
            stderr_write(['WARNING: No threshold provided',
                          'for network method', self.method,
                          '\nDefault set to:', self.thresh])
        #has not been implemented, but would be useful to control min_co_occur
#        if not hasattr(self, 'min_co_occur'):
#            self.min_co_occur = 2
#            stderr_write(['WARNING: No minimum co-occurrence',
#                          'count provided for network method', 
#                          self.method,
#                          '\nDefault set to:', self.min_co_occur])
#
#        elif not to_bool(self.min_co_occur):
#            self.min_co_occur = 2
#            stderr_write(['WARNING: No minimum co-occurrence',
#                          'count provided for network method', 
#                          self.method,
#                          '\nDefault set to:', self.min_co_occur])

        # default dictionary of patient IDs and Biopython SeqRecords for patients 
        self.patient_records = self.make_patient_records()
        
        # sets sequence length based on alignment sequences provided to yaml
        self.seq_length = self.check_seq_lengths()
        
    
    def check_seq_lengths(self):
        ''' Check that the length in all sequences for provided alignment file
            are all the same. If they are, returns the length of all of them
        '''
        seq_lengths = [len(record.seq) for _k, record_list in self.patient_records.items()
                       for record in record_list]
        assert max(seq_lengths) == min(seq_lengths)
        return min(seq_lengths)
    
    def make_patient_records(self):
        ''' Create default dictionary with key as unique patient id from study
            and list of Biopython SeqRecords corresponding to patient id
        '''
        patient_records = defaultdict(list)
        with open(self.aln_fasta_file, 'rU') as handle:
            for index, record in enumerate(SeqIO.parse(handle, 'fasta')):
                try:
                    '''
                    Assumes that the alignment FASTA headers hold this form
                    >968_12446_...
                    where "968_12446" is a unique identifier for a patient
                    '''
                    a = record.name.split('_')   
                    study_patient_id = a[0]+'_'+a[1]
                    patient_records[study_patient_id].append(record)        
                except IndexError:
                    print >> sys.stderr, 'Error with fasta header: %s' % record.name
        return patient_records   

def assign_drug_type(PosNetInput, 
                     nNRTI = set(['ETR', 'NVP', 'RPV', 'EFV', 'NNRTI', 'DLV', 'AAPA']),
                     NRTI = set(['C',
                                 'ABC',
                                 'TAF',
                                 '3TC',
                                 'ZDV',
                                 'DTG',
                                 'FTC',
                                 'TDF',
                                 'DDI',
                                 'EVG',
                                 'D4T',
                                 'DDC',
                                 'NRTI',
                                 'AZT'])):
    '''
    Classifies patient with drug treatment program by parsing from FASTA headers
    of form 
    >~~~~~~__DrugName1_DrugName2_DrugName3
        
    and returns a dictionary where keys are unique patient ID strings
    and the values are a list of drug types used for each SeqRecord of each patient
    '''
    patient_drug_assignment = defaultdict(list)
    for patient_id, record_list in PosNetInput.patient_records.items():
        record_assignments = []
        for i,record in enumerate(record_list):
            try:
                record_drugs = record.name.split('__')[1].split('_')
                record_drugs = [d.upper() for d in record_drugs]

                record_drugs = set(record_drugs) - set(['NONE'])
                NRTI_intersect = NRTI & record_drugs
                nNRTI_intersect = nNRTI & record_drugs
                
                if record_drugs == NRTI_intersect:
                    record_assignments.append('NRTI')
                elif record_drugs == nNRTI_intersect:
                    record_assignments.append('nNRTI')
                elif record_drugs == (NRTI_intersect | nNRTI_intersect):
                    record_assignments.append('Both')
                elif not (NRTI_intersect | nNRTI_intersect):
                    record_assignments.append('None')
                else:
                    record_assignments.append('N/A')
                    print >> sys.stderr, \
                        record.name + ' thrown out due to non (n)NRTI drug present.'
            except IndexError:
                record_assignments.append('N/A')
                print >> sys.stderr,\
                    'Could not split fasta header: %s' % (record.name)
                continue
        patient_drug_assignment[patient_id] = record_assignments
    #print(patient_drug_assignment)    
    return patient_drug_assignment
    
def assign_patient_drugs(PosNetInput, 
                     nNRTI = set(['ETR', 'NVP', 
                                  'RPV', 'EFV', 
                                  'NNRTI', 'DLV', 'AAPA']),
                     NRTI = set(['C','ABC',
                                 'TAF','3TC',
                                 'ZDV','DTG',
                                 'FTC','TDF',
                                 'DDI','EVG',
                                 'D4T','DDC',
                                 'NRTI','AZT']) ):
    '''
    Classifies patient with drug treatment program by parsing from FASTA headers
    of form 
        >~~etc~~~~__DrugName1_DrugName2_DrugName3
    and returns a dictionary where keys are unique patient ID strings
    and the values are a list of drugs used for each SeqRecord for each patient
    '''
    patient_drug_assignment = dict()
    for patient_id, record_list in PosNetInput.patient_records.items():
        record_assignments = []
        for i,record in enumerate(record_list):
            try:
                record_drugs = record.name.split('__')[1].split('_')
                record_drugs = [d.upper() for d in record_drugs]
                '''
                if len(record_drugs) == 1:
                    if record_drugs[0] == 'None':
                        record_assignments.append('None')
                        continue
                '''
                
                # some headers had drugs listed followed by None for some reason
                record_drugs = set(record_drugs) - set(['NONE'])
                record_assignments.append(list(record_drugs))
            except IndexError:
                record_assignments.append([])
                print >> sys.stderr,\
                    'Could not identify drugs used for %s' % (record.name)
                continue
        patient_drug_assignment[patient_id] = record_assignments
        #patient_drug_assignment[patient_id] = set(record_assignments)
    return patient_drug_assignment
    #print(patient_drug_assignment)

def write_drug_set_to_fasta(patient_drug_dict, PosNetInput):
    '''
    Writes  patient drug classification calculated by assign_drug_type to fastas
    '''
    with open('./NRTI_aligned_big_aids.fasta', 'w') as NRTI,\
    open('./nNRTI_aligned_big_aids.fasta', 'w') as nNRTI,\
    open('./both_aligned_big_aids.fasta', 'w') as both,\
    open('./none_aligned_big_aids.fasta', 'w') as none:
        thrown_out = 0
        NRTI_patients = 0
        NRTI_seqs = 0
        nNRTI_patients = 0
        nNRTI_seqs = 0
        both_patients = 0
        both_seqs = 0
        none_patients = 0
        none_seqs = 0
        for patient_id, record_list in patient_drug_dict.items():
            record_set = set(record_list[1:]) #skips position 0 because know's its None
            if record_set == set(['NRTI']):
                SeqIO.write(PosNetInput.patient_records[patient_id], NRTI, 'fasta')
                NRTI_patients += 1
                NRTI_seqs += len(PosNetInput.patient_records[patient_id])
            elif record_set == set(['nNRTI']):
                SeqIO.write(PosNetInput.patient_records[patient_id], nNRTI, 'fasta')
                nNRTI_patients += 1
                nNRTI_seqs += len(PosNetInput.patient_records[patient_id])
            elif record_set == set(['Both']):
                SeqIO.write(PosNetInput.patient_records[patient_id], both, 'fasta')
                both_patients += 1
                both_seqs += len(PosNetInput.patient_records[patient_id])
            elif record_set == set(['None']):
                SeqIO.write(PosNetInput.patient_records[patient_id], none, 'fasta')
                none_patients += 1
                none_seqs += len(PosNetInput.patient_records[patient_id])
            else:
                thrown_out += 1
    print >> sys.stderr, '%i patients thrown out' % thrown_out
    print >> sys.stderr, '%i patients using NRTI drugs, %i sequences'\
             % (NRTI_patients, NRTI_seqs)
    print >> sys.stderr, '%i patients using nNRTI drugs, %i sequences'\
             % (nNRTI_patients, nNRTI_seqs)
    print >> sys.stderr, '%i patients using both nNRTI and NRTI drugs, %i sequences'\
             % (both_patients, both_seqs)
    print >> sys.stderr, '%i patients not using drugs, %i sequences'\
             % (none_patients, none_seqs)

def calc_drug_patient_count(AlnPosPairSet, patient_drug_dict, 
                            list_of_drugs = ['ETR', 'NVP', 
                                             'RPV', 'EFV', 
                                             'NNRTI', 'DLV', 
                                             'AAPA','UNKNOWN', 
                                             'C','ABC',
                                             'TAF','3TC',
                                             'ZDV','DTG',
                                             'FTC','TDF',
                                             'DDI','EVG',
                                             'D4T','DDC',
                                             'NRTI','AZT'],
                            outfile=''):
    ''' 
    Creates dataframe with mutation positions as indices, and drugs at columns.
    Cell(i,j) would be the # of patients that had mutation at a position i and
    also took drug j
    '''
    node_drugs_df = pd.DataFrame(index=list(range(1, AlnPosPairSet.seq_length+1)),
                                 columns=list_of_drugs)
    node_drugs_df.fillna(0, inplace=True)
    patient_mutation_dict = AlnPosPairSet.patient_mutations
    
    # mutatations are represented like this example
    # pos   1 2 3 4 5 6
    # t=0   P I S P I E
    # t=1   K I S P I E
    #        
    # mut_table = 
    #     [[0 0 0 0 0 0]
    #      [1 0 0 0 0 0]]
    for patient_id, mut_table in patient_mutation_dict.items():
        drug_records = patient_drug_dict[patient_id] #list of drugs used by patient
        assert len(drug_records) == len(mut_table)
        #assert len(mut_table[(i,)]) == AlnPosPairSet.seq_length
        for i,drug_list in enumerate(drug_records):
            positions = np.where(mut_table[(i,)] == 1)[0] #positions NOT zero-indexed
            if list(positions) and list(drug_list): #if neither is empty
                node_drugs_df.loc[positions, drug_list] += 1
    # remove a columns and rows that only have zeros in them
    node_drugs_df = node_drugs_df.loc[(node_drugs_df.sum(axis=1) != 0), (node_drugs_df.sum(axis=0) != 0)]
    if outfile:
        node_drugs_df.to_csv(outfile, index=True, header=True)
    return node_drugs_df
    

def run_aln_pos_pairs(config):
    '''
    Processes patient-specific position mutation network input information.
    Generates and writes network to text file(s).
    '''
    # Obtains alignment file and parameters 
    # for alignment based network reconstruction
    inp = PosNetInput(config)
    print('Calculated pospairs input for %i patients' % len(inp.patient_records))    
    
    # #Uncomment line below if want to do drug type subsetting for patients
    #write_drug_set_to_fasta(assign_drug_type(inp), inp)
    

    # Calculates edge weights for co-occuring mutation pairs
    aln_mut_pair_set = AlnPosPairSet(inp)
    print('Edge weights calculated')
    
        
    # #Uncomment line below if want to export matrix of patient counts for
    # #mutation position vs. drug administered to file
    #calc_drug_patient_count(aln_mut_pair_set, assign_patient_drugs(inp), 
    #    outfile=aln_mut_pair_set.aln_fasta_file+'_node_drug_counts.csv')
    
    # Output network in tab-delimited format
    # Source\tTarget\tWeight
    aln_mut_pair_set.write_network_to_file(inp.net_file)
    print('Writing network to %s' % inp.net_file)
    