import argparse
import sys
import itertools

from NAPA.mutpairs.UndirectedPairClass import *


def parseArgs():
        parser = argparse.ArgumentParser()
        parser.add_argument('-a', dest = 'alnFASTA', required = True,
                            help = 'Path of multiple sequence alignment fasta file.')
        parser.add_argument('-fi', dest = 'seqId', required = True, default = 'None',
                            help = 'IDs of sequences in alignment studied. Default is all sequences.')
        parser.add_argument('-af', dest = 'alnFunc', required = True,
                            help = 'File containing known functions of sequences in alignment.')
        parser.add_argument('-fn', dest = 'functNum', required = True,
                            help = 'File containing numerical value to name of function.')
        parser.add_argument('-f', dest = 'function', required = True,
                            help = 'String/File containing functions to be included in net.')
        parser.add_argument('-pa', dest = 'pos2aa', required = True,
                            help = 'File containing protein position vs. wt residue corresponding to aln.')
        parser.add_argument('-p', dest = 'positionList', required = True,
                            help = 'File containing protein positions to be considered in network')
        parser.add_argument('-nf', dest = 'netFile', required = True,
                            help = 'Outputfile for result network.')
        args = parser.parse_args()

        return args

def stderr_write(myList):
    sys.stderr.write(' '.join([str(el) for el in myList])+'\n')


def parse_key_val(file_name):
    out_dict = {}
    with open(file_name, 'w') as f:
        for line in f:
            out_dict[line.strip().split('\t')[0]] = line.strip().split('\t')[1].strip()
    return out_dict

def main():
    args = parseArgs()
    num_function = parse_key_val(args.funcNum) # numerical function representation
    leafnode_func_num = parse_key_val(leafnode_func_file) # num. func. for each leafnode in aln
    
    # leafnode to function name (R APE ACE only accepts numerical arguments for anc.reconst)
    leafnode_func = {leaf:num_function[num] for leaf, num in leafnode_func_num.iteritems()}
    
    # list of positions considered for network.
    # (Remove residues from signaling peptide, or other domains not studied).
    with open(args.positionList,'r') as f:
        pos_list = [int(line.strip()) for line in f] 

    # position and corresponding wildtype residue at position
    pos_wt = parse_key_val(args.pos2aa)

    # Construct network based on alignment
    undirMutationPairs = UndirPairSet([args.alnFASTA], pos_list, 
                                      pos_wt, leafnode_function)

    stderr_write(['\t\tNumber of undir (mut) pairs', 
                  len(undirMutationPairs.undirPairs)])
    for mutPair in undirMutationPairs.undirPairs:
        mutPair.get_FisherPvals()
        
    with open(args.netFile, 'wb') as f:
        f.write(undirMutationPairsClin.str_jaccard())
    
if __name__  ==  '__main__':
    main()
