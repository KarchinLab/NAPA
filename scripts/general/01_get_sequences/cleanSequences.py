import argparse
import sequenceManipulation as seqManip

def parseArguments():
    """                                                                                                                    
    Builds the argument parser object                                                                                      
    """
    parser = argparse.ArgumentParser(prog="cleanSequences.py",
                                     description="Remove non-unique sequences from a fasta sequence alignment.",
                                     usage="cleanSequences.py -i inFasta -o outFasta")
    parser.add_argument("-i", dest="inFasta", required=True,
                        help="input fasta file")
    parser.add_argument("-o", dest="outFasta", required=True,
                        help="output fasta file")
    args = parser.parse_args()
    return args


def main():
    args = parseArguments()
    FAinputFile = args.inFasta
    FAoutputFile= args.outFasta
    id2seq = seqManip.get_FASTA_records(FAinputFile) # read fasta recs into dictionary
    replacePairs = [['/','_'],['|','_'],['-',''],['TEM00','TEM'],['TEM0','TEM']]
    truncStart,truncEnd = 0,6
    cleanId2seq = seqManip.fasta_sequence_cleaner(id2seq,replacePairs,truncStart,truncEnd) # clean non-unique seqs, ids truncated at 10 chars
    seqManip.print_FASTA_sequences_ordered(cleanId2seq,FAoutputFile) # print sequences

if __name__ == "__main__": main()
