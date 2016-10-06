import os
from itertools import groupby
import re


def fasta_iter(fasta_name):
    '''
    Given a fasta file. yield tuples of header,  sequence
    '''
    with open(fasta_name, 'rb') as fh:
        faiter = (x[1] for x in groupby(fh,  
                                        lambda line: line[0] == ">"))
    
        for header in faiter:
            header = header.next()[1:].strip() # drop the '>'
            #join all sequence lines into one string
            seq = "".join(s.strip() for s in faiter.next()) 
            yield header,  seq

def fasta_to_dict(fasta_input_name):
    ''' Reads: Fasta-formatted sequence file.
    Returns dictionary: {fasta_id: sequence_string}
    '''
    if os.path.isfile(fasta_input_name):
        fasta_sequences = fasta_iter(fasta_input_name)
        id_sequence_dict = {}
        for fs in fasta_sequences:
            id_sequence_dict[fs[0]]=fs[1]
        return id_sequence_dict
    else:
        return {}
