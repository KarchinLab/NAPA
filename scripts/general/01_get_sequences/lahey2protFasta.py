# Author: V Beleva Guthrie
# Created: 04/09/2013
# Last Modified: 04/09/2013
# Modified by: VBG
#
# Reads in the latest TEM table from Lahey, which has been manipulated, such that:
# The first two rows are:
#       Lactamase GenBank_Accession_Number Phenotype 6 21 35 ...
#       TEM-001   J01749;AF397068          2b        Q  L  D ...
# It generates a set of FASTA sequences corresponding to the TEM mutants
# Needs: Option and argument parsing
# Check for correctness
 
import csv
import pprint
import sequenceManipulation as seqManip

tem1seq = "MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQE"
tem1seq += "QLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPN"
tem1seq += "DERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVI"
tem1seq += "YTTGSQATMDERNRQIAEIGASLIKHW"

def readLahey(laheyFile):
    """
    Reads the Lahey file into a dictionary of 
    key = TEM Id for a given TEM sequence
    value = mutations in that TEM sequence 
            mutations are in the form:
            {temPosition:(wildtypeResidue,mutatedResidue)} 
    """
    col2temPos = {} #col(key) the column in the csv file; temPos(val):the TEM position for that column
    temPos2wt = {} # TEM position to wild-type residue
    tem2mut = {} #holds tem sequence id vs mutations(value) in that sequence(key)
    with open(laheyFile, 'rb') as csvFile:
        laheyReader = csv.reader(csvFile, delimiter=',')#an iterator through the lines
        row = laheyReader.next()
        col2temPos = dict([ (i,row[i]) for i in range(3,len(row)) ])
        print col2temPos
        row = laheyReader.next()
        col2wt = dict([ (i,row[i]) for i in range(3,len(row)) ])
        for row in laheyReader: # row is a list of words(strings)
            col2mut = dict([ (i+3,e) for i,e in enumerate(row[3:]) ])
            dictOfMutations = {} # {temPosition:(wildtypeResidue,mutatedResidue)}
            for col in col2mut:
                mutString = ""
                if col in col2temPos and col in col2wt:
                    pos = int(col2temPos[col]) if len(col2temPos[col]) else 0
                    wtRes = col2wt[col]
                    mut = col2mut[col]
                    if pos and len(wtRes) and len(mut):
                        dictOfMutations[pos] = (wtRes,mut)
            tem2mut[row[0].strip()] = dictOfMutations.copy()
    return tem2mut

def temPositions():
    """ Produces a dictionary of 
    Ambler numbering of TEM protein sequence positions
    vs. the actual position in the TEM sequence string
    """
    startPos = 3
    endPos = 290
    skippedPos = [239,253]
    temPos = [] # contains all Ambler TEM positions
    for pos in range(startPos,endPos+1):
        if pos not in skippedPos: temPos.append(pos)
    ambler2temPos = dict(zip(temPos,range(len(tem1seq))))# Check This!
    return ambler2temPos

def lahey2sequenceDict(tem2mut,ambler2temPos):
    # initialize mutant sequences to wild-type
    temMut2seq = dict.fromkeys(tem2mut.keys(), tem1seq)
    temMut2seq["TEM-001"] = tem1seq # add the TEM-1 sequence
    # then introduce mutations
    for tem in tem2mut:
        for pos in tem2mut[tem]:
            currSequence = temMut2seq[tem]
            currSeqList = list(currSequence)
            currSeqList[ambler2temPos[pos]]=tem2mut[tem][pos][1]
            temMut2seq[tem] = ''.join(currSeqList)# CHECK THIS!
    return temMut2seq

def main():
    laheyFile = "2be_sequs.csv"
    temFastaFile = "2be_sequs.fa"
    tems2mutations = readLahey(laheyFile)
    print "TEM2mut"
    pprint.pprint(tems2mutations)
    ambler2temPos = temPositions()
    print "AMBLER2TEMposition:"
    print ambler2temPos
    tems2sequences = lahey2sequenceDict(tems2mutations,ambler2temPos)
    print tems2sequences
    seqManip.print_FASTA_sequences_ordered(tems2sequences,"TEM-001",temFastaFile)    

if __name__ == "__main__": main()            
