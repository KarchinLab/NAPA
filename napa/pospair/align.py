# -*- coding: utf-8 -*-

from collections import defaultdict
import numpy as np

from napa.utils.serials import * 
from napa.utils.io import *

class AlnPosPairSet(object):
    '''
    Set of mut pairs and associated metrics from an MSA.
    '''
    def __init__(self, inp):
        #inherits same attributes from PosNetInp class
        self.__dict__.update(vars(inp))
        
        #each entry consists of a {(source_node, target_node) : edge_weight_val}
        self.edge_weights = defaultdict(int)
        
        #key is unique patient ID string, value is table of mutations 
        #see update_edge_weights(..)
        self.patient_mutations = dict()
        
        #loops through each patient and their longitudinal sequences
        for patient_id, record_list in self.patient_records.items():
            if record_list:
                wt_record = record_list[0]
                
                #update edge weight dict
                mutation_table = self.update_edge_weights(record_list, wt_record)
                
                #add to mutation table
                self.patient_mutations[patient_id] = mutation_table
    
    def __repr__(self):
        out_str = ''
        for (e1,e2), raw_count in self.edge_weights.items():
            if raw_count >= self.thresh:
                out_str += \
                '%s\t%.6f\n' % (str(e1)+'\t'+str(e2), 
                                raw_count)

        return out_str
    
    def update_edge_weights(self, record_list, wt_record, 
                            ignore_codes=['X','-','B','Z','J','O','U']):
        '''
        Returns matrix of positions that mutated in a list of SeqRecords
        ignoring mutations from/to X, J, O, U, B, Z and gaps (-)
        
        Assumes record_list is ordered temporally (oldest time to newest)
        For future probably best to check index in fasta header
        '''
        
        # What is a mutation table?
        #
        # mutatations are represented like this example
        # pos   1 2 3 4 5 6
        # t=0   P I S P I E
        # t=1   K I S P I E
        #        
        # mut_table = 
        #     [[0 0 0 0 0 0]
        #      [1 0 0 0 0 0]]
        
        mut_table = np.zeros((len(record_list), self.seq_length))
        previous_mutations = []
        previous_record = []
        for i,record in enumerate(record_list):
            if i==0:
                previous_record = wt_record #wildtype to compare first sequence to
            else:
                previous_record = record_list[i-1]
                
            # compare wildtype sequence with record sequence and output changed positions
            mutated_positions = [j for j in xrange(self.seq_length) 
                                 if str(previous_record.seq)[j] != str(record.seq)[j]]
           
            # remove any mutation positions that contain blacklisted codes
            # add 1 to preserve biological position indexing 
            # (don't want python zero indexing)
            mutated_positions = [p+1 for p in mutated_positions 
                                 if str(record.seq)[p] not in ignore_codes 
                                    and str(previous_record.seq)[p] not in ignore_codes]
            #update mutation table
            mut_table[(i,mutated_positions)] += 1
            
            #print(previous_mutations)
            #print(i,mutated_positions)
            
            if mutated_positions:
                # update edge weights for mutations in same time point
                for e1 in mutated_positions:
                    for e2 in list(mutated_positions):
                        if e1 != e2:
                            self.edge_weights[(e1,e2)] += 1
                            #print('%f -> %f incremented to %f' % 
                            #      (e1,e2, self.edge_weights[(e1,e2)]))
                
                # update edge weights for mutations at current time point 
                # w/ mutations at previous time point
                if previous_mutations:
                    for i in previous_mutations:
                        for j in mutated_positions:
                            if i != j: # condition that ensures no self loops
                                self.edge_weights[(i,j)] += 1
                                #print('%f -> %f incremented to %f' % 
                                #      (i,j, self.edge_weights[(i,j)]))
#                            elif i==j:
#                                print('Self-loop at %i for %s' % (i, record.name))
#                                print('Protein at previous time: %s' 
#                                        % previous_record.seq[i-1])
#                                print('Protein mutated at current time: %s' 
#                                        % record.seq[i-1])
#                                print('\n')
                                
            previous_mutations = mutated_positions
        return mut_table
    
    def write_network_to_file(self, file_path):
        '''
        Writes network to file with path file_path
        Columns: source_node  target_node  weight
        '''
        with open(file_path, 'wb') as f:
            f.write(str(self))    
        