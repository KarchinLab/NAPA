#! /usr/bin/env python 
'''Custom functions for writing sequence files. 
Should be replaced with optimized inbuilt functions from PyCogent(GPL).'''

def write_wrap_fasta_seq(seqId,  seqStr,  wrap_at = 80,  
                         endline = '\n'):
    ''' 
    Print FASTA sequence and wrap at desired character count 
    '''
    s = []
    s.append('>%s%s' % (seqId.lstrip('>'),  endline))
    exploded_seq = list(seqStr)
    wrap_points = range(0, len(exploded_seq), wrap_at)
    wrap_points.reverse()
    for i in wrap_points[:-1]:
        exploded_seq.insert(i,  endline)
    s = s + exploded_seq + [endline]
    return ''.join(s)


def write_wrap_fasta_dict(id_sequence_dict={}, 
                          fasta_output_name='fasta_dict.fa',
                          wrap_at=80, sort=False):
    ''' 
    Prints FASTA sequences from id vs. sequence dictionary 
    allows for sorting fasta records alphabetically by fasta id.
    '''
    with open(fasta_output_name,  'w') as f:
        if not sort:
            for seq_id in id_sequence_dict:
                f.write(write_wrap_fasta_seq(seq_id, 
                                         id_sequence_dict[seq_id], 
                                         wrap_at))
        else:
            for seq_id in sorted(id_sequence_dict.keys()):
                f.write(write_wrap_fasta_seq(seq_id, 
                                         id_sequence_dict[seq_id], 
                                         wrap_at))
                                        
            
