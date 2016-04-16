projDir=/home/vbeleva/phylogeny/pairs
runDir=$projDir/data/old
workDir=$projDir/undirectedPairs
globalFasta=$runDir/PSE4-TEM-SHV_noLen.cds.prot.fasta
maxTrees=$2
maxRuns=$1
echo $maxRuns $maxTrees
outPath=$workDir/results/$maxRuns
mkdir -p $outPath
phenoFile=$workDir/TEM_pheno.txt

distThreshs=( 0.05 )
for d in "${distThreshs[@]}"
do
    nohup python get_undirectedPairs_single.py -d $runDir -r $maxRuns -m $maxTrees \
	-g $globalFasta -p $phenoFile -t $d -o $outPath >& getEpiPairs.$maxRuns.$maxTrees.$d.oe &
done

