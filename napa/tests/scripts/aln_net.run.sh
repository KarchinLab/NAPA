projDir=$1
projDir=$projDir/napa
dataDir=$projDir/tests/data
resultsDir=$projDir/tests/results
analyzeDir=$projDir/analyze

#------------Input files------------------
# Required: Input MSA in FASTA format
alnFASTA=$dataDir/classA.protein.fasta 

# Optional: Sequential numbering of positions
posList=$dataDir/TEM-pos2aa.txt  

# Optional: Subset of positions
posSubset=$dataDir/TEM-posSubset.txt

# Optional: Functional assignment for each sequence
protFunc=$dataDir/TEM-functions.txt

# Optional: List of selected functions assoc. with sequences
selProtFunc=$dataDir/sel-funcs.txt

# Required output file path: the network in tab-delimited format
netFile=$resultsDir/TEM-aln-net.txt

#Alignment network construction and analysis
#Jaccard index weighted undirected pairs
python $analyzeDir/aln_mut_pairs.py \
    -a  $alnFASTA \
    -p  $posList \
    -wt "TEM_1" \
    -ps $posSubset \
    -pf $protFunc \
    -sf $selProtFunc \
    -nf $netFile
