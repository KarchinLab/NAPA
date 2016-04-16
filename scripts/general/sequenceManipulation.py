from itertools import groupby
import re
import numpy as np

def fasta_iter(fasta_name):
    """given a fasta file. yield tuples of header, sequence"""
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip() # drop the '>'
        seq = "".join(s.strip() for s in faiter.next()) #join all sequence lines into one string
        yield header, seq

def entry2seqStr(seqId, seqStr, wrap_at=80, endline='\n'):
    """ print FASTA sequence and wrap at desired character count"""
    s = []
    s.append('>%s%s' % (seqId.lstrip('>'), endline))
    exploded_seq = list(seqStr)
    wrap_points = range(0,len(exploded_seq),wrap_at)
    wrap_points.reverse()
    for i in wrap_points[:-1]:
        exploded_seq.insert(i, endline)
    s = s + exploded_seq + [endline]
    return ''.join(s)

def get_FASTA_records(FAinputFile):
    fasta_sequences = fasta_iter(FAinputFile)
    id2seq = {}
    for fs in fasta_sequences:
        id2seq[fs[0]]=fs[1]
    return id2seq

def print_FASTA_sequences(id2seq,mainID,FAoutputFile):
    """ prints fasta sequences in an ortholog set"""
    with open(FAoutputFile, 'w') as faOut:
        faOut.write(entry2seqStr(mainID, id2seq[mainID]))
        for otherID in id2seq.keys():
            if otherID != mainID:
                faOut.write(entry2seqStr(otherID,id2seq[otherID]))

def print_FASTA_orthologs_ordered(id2seq,mainID,FAoutputFile):
    """ prints fasta sequences in an ortholog set"""
    with open(FAoutputFile, 'w') as faOut:
        faOut.write(entry2seqStr(mainID, id2seq[mainID]))
        for otherID in sorted(id2seq.keys()):
            if otherID != mainID:
                faOut.write(entry2seqStr(otherID,id2seq[otherID]))

def print_FASTA_sequences_ordered(id2seq,FAoutputFile):
    with open(FAoutputFile, 'w') as faOut:
        for seqID in sorted(id2seq.keys()): faOut.write(entry2seqStr(seqID,id2seq[seqID]))

def rmGaps(seq):
    chars2remove = ['-','*','\n']
    for ch in chars2remove: seq = seq.replace(ch, '')
    return seq

def percentSequenceIdentity(seq1,seq2):
    """ computes the sequence identity of two sequences without aligning"""
    seq1 = rmGaps(seq1)
    seq2 = rmGaps(seq2)
    if seq1 == seq2: return 100.0
    if len(seq1) < len(seq2):
        if seq1 in seq2: return 100.0*float(len(seq1))/float(len(seq2))
        identCnt=0
        for i in range(len(seq1)):
            if seq1[i]==seq2[i]: identCnt += 1
        return float(identCnt)/float(len(seq2))
    if len(seq2) < len(seq1):
        if seq2 in seq1: return 100.0*float(len(seq2))/float(len(seq1))
        identCnt=0
        for i in range(len(seq2)):
            if seq2[i]==seq1[i]: identCnt += 1
        return float(identCnt)/float(len(seq1))

def fasta_sequence_cleaner(seqDict,replacePairs,truncPosStart,truncPosEnd):
    """Cleans a set of fasta sequences from an MSA:
    * Removes duplicate sequences with different sequence IDs (uses dictionary inversion)
    * Removes some punctuation characters from sequence IDs (uses replaceDict)
    * Makes IDs no longer than 10 chars (useful for PHYLIP formats)
    * Truncates by truncatePosStart,truncatePosEnd positions at start or end of each sequence
    """
    invSeqDict = {}
    uniqSeqDict = {}
    print "Starting with", len(seqDict), "sequences"

    #when two different IDs have same sequence, the ID alphabetically smaller ID is chosen
    for seqId in sorted(seqDict.keys(), reverse=True):
        newSeqId = str(seqId)
        for s in range(len(replacePairs)):
            newSeqId = newSeqId.replace(replacePairs[s][0], replacePairs[s][1])
        truncSeq = seqDict[seqId][truncPosStart:-truncPosEnd]
        invSeqDict[truncSeq] = newSeqId[0:10]

    uniqSeqDict = {seqId:seq for seq,seqId in invSeqDict.iteritems()}

    print "Ending with", len(uniqSeqDict), "sequences"
    return uniqSeqDict
        
def seqLenDiff(seq1, seq2):
    """ chooses longer sequence from two sequences from same organism 
    and protein but from different databases"""
    seq1 = rmGaps(seq1)
    seq2 = rmGaps(seq2)
    return len(seq1) - len(seq2)

def mutLists2matrix(mutLists=[], posList=[]):
    if not len(mutLists): return []
    if not len(posList):
        posList = sorted(set([int(re.sub("[^0-9]", "",mut)) \
                              for mutList in mutLists for mut in mutList]))
    pos2posNum = dict(zip(posList,range(len(posList))))
    outMutMatrix = np.zeros((len(mutLists),len(posList)))
    for i in range(len(mutLists)):
        for mut in mutLists[i]:
            outMutMatrix[i,pos2posNum[int(re.sub("[^0-9]", "",mut))]]=1
    return outMutMatrix, posList
