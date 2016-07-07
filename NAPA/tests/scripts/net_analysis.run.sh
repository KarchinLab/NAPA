projDir=$1
projDir=$projDir/NAPA
dataDir=$projDir/tests/data
resultsDir=$projDir/tests/results
scriptsDir=$projDir/tests/scripts


#Input arguments
mutPairType='undir'
#mutPairType='dir'
pathLen=3


#Input file
netBase=TEM-phylo-net
netFile=$resultsDir/$netBase.$mutPairType.txt

#Output file(s)
pathCentOutFile=$resultsDir/$netBase.$mutPairType.path_betw.txt

#Network analysis
python $scriptsDir/net_analysis.py \
    -nf $netFile \
    -mt $mutPairType \
    -pl $pathLen \
    -pc $pathCentOutFile 



