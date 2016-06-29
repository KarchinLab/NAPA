projDir=./
dataDir=$projDir/data
resultsDir=$projDir/results
scriptsDir=$projDir/scripts
alnDir=$dataDir/MSA

msaFASTA=$dataDir/classA.codingDNA.fasta
fastaIDs=$dataDir/TEM-clin-ids.txt
msaFunctions=$dataDir/TEM-functions-num.txt
pos2aa=$dataDir/TEM-pos2aa.txt
posList=$dataDir/TEM-posList.txt

#Alignment network construction and analysis
#Jaccard index weighted undirected pairs
python msa_mut_pairs.py \
    -a $msaFASTA \
    -fi $fastaIDs \
    -af $msaFunctions \
    -fn  $funcNum \
    -f $func \
    -pa  $pos2aa \
    -p $posListm \
    -nf $netFile
