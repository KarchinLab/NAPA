projDir=/home/vbeleva/phylogeny/pairs
runDir=$projDir/data
workDir=$projDir/alignmentPairs
clinFasta=$runDir/PSE4-TEM-SHV_noLen.cds.prot.fasta
labFasta=$runDir/tem_direvo.fasta
outPath=$workDir/results/
mkdir -p $outPath
phenoFile=$workDir/TEM_pheno_clinDirEvo.txt

nohup python get_undirectedPairs_single.py -c $clinFasta -l $labFasta -p $phenoFile -o $outPath >& getAlignmentPairs.oe &


