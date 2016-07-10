projDir=$1
projDir=$projDir/NAPA
dataDir=$projDir/tests/data
resultsDir=$projDir/tests/results
scriptsDir=$projDir/tests/scripts

#Input file
alnFASTA=$dataDir/classA.protein.fasta
#fastaIDs=$dataDir/TEM-clin-ids.txt
posList=$dataDir/TEM-pos2aa.txt
posSubset=$dataDir/TEM-posSubset.txt
protFunc=$dataDir/TEM-functions.txt
selProtFunc=$dataDir/sel-funcs.txt

#Output file: the network in tab-delimited format
netFile=$resultsDir/TEM-aln-net.txt

#Alignment network construction and analysis
#Jaccard index weighted undirected pairs
python $scriptsDir/aln_mut_pairs.py \
    -a  $alnFASTA \
    -p  $posList \
    -wt "TEM_1" \
    -ps $posSubset \
    -pf $protFunc \
    -sf $selProtFunc \
    -nf $netFile
