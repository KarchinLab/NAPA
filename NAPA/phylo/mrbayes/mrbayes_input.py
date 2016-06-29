# Author: V Beleva Guthrie
# Created: 05/28/2013
# Last Modified: 05/28/2013 by: VBG

import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped

MBstr1 = '\n\t'.join(['begin mrbayes;','log start replace;','set autoclose = no nowarn=yes;'])
MBstr2 = '\n\t'.join(['\ncharset 1st_pos = 1-.\\3;','charset 2nd_pos = 2-.\\3;','charset 3rd_pos = 3-.\\3;',
                      'partition by_codon = 3:1st_pos,2nd_pos,3rd_pos;','set partition = by_codon;',
                      'lset applyto = (all) nst=6 rates = gamma;',
                      'prset applyto = (all); [For JC add opotion statefreqpr = fixed(equal]',
                      'unlink revmat=(all) shape=(all) pinvar=(all) statefreq=(all) tratio= (all);',
                      'mcmc ngen=5000000 nchains=8 temp=0.04 printfreq=1000 samplefreq=100 savebrlens=yes;',
                      'sumt burnin=1250 contype = halfcompat;','log stop;'])
MBstr3 = 'end;\n'

def parseArguments():
    """
    Builds the argument parser object
    """
    parser = argparse.ArgumentParser(prog="fasta2nexus.py", 
                                     description="Convert a fasta alignment to nexus format.",
                                     usage="fasta2nexus.py -b baseFileName -s alnSuffix -o outGroup")
    parser.add_argument("-b", dest="baseFileName", required=True,
                        help="The base file name (minus extension) of the aln and tree files")
    parser.add_argument("-s", dest="alnSuffix",required=True,
                default="fasta", help="the suffix of the aligned fasta file [fasta]")
    parser.add_argument("-o", dest="outGroup", required=True,
                        help="The taxon which will be used as an outgroup")
    args = parser.parse_args()
    return args

def convert_seq_format(inName, outName, inForm, outForm):
    input_handle = open(inName, "r")
    output_handle = open(outName, "w")
 
    sequences = SeqIO.parse(input_handle, inForm, alphabet=Gapped(IUPAC.unambiguous_dna))
    count = SeqIO.write(sequences, output_handle, outForm)
 
    output_handle.close()
    input_handle.close()
 
    print "Converted %i records" % count

def nexus_for_MrBayes(nexBaseFilename, nexSuffix, mbSuffix):
    linestring = open(nexBaseFilename+'.'+nexSuffix, 'r').read()
    with open(nexBaseFilename+'.'+mbSuffix, 'w') as mrBayesNexFile:
        mrBayesNexFile.write(linestring.replace('\'', ''))

def mrBayesBatch(nexBaseFilename,nexSuffix,mbSuffix,outGroup):
    with open(nexBaseFilename+'.'+mbSuffix, 'w') as mrBayesBatchFile:
        execute = "execute "+nexBaseFilename.split('/')[-1]+'.bay.'+nexSuffix+';'
        outgroup = "outgroup "+outGroup+';'
        mrBayesBatchFile.write('\n\t'.join([MBstr1,execute,outgroup,MBstr2])+'\n'+MBstr3)

def main():
    args = parseArguments()
    inFileName = args.baseFileName + '.' + args.alnSuffix
    outFileName = args.baseFileName + '.' + "nexus"
    convert_seq_format(inFileName,outFileName,"fasta","nexus")
    nexus_for_MrBayes(args.baseFileName, "nexus", "bay.nexus")
    mrBayesBatch(args.baseFileName, "nexus", "bay.batch", args.outGroup)
    
if __name__ == "__main__": main()
