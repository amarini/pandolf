#!/bin/bash 
# usage: ./makeFileLists.sh [dirname] [analyzerType=PhotonJet] [recoType="pf"] [jetAlgo="akt5"]

maindir=$1
if [ $# -eq 1 ]
then
analyzerType="PhotonJet"
else
analyzerType=$2
fi

if [ $# -lt 3 ]
then
recoType="pf"
else
recoType=$3
fi

if [ $# -lt 4 ]
then
jetAlgo="akt5"
else
jetAlgo=$4
fi

echo "scanning $maindir"

for i in `ls $maindir`
do 
rm -f files_${analyzerType}_2ndLevel_${i}_${recoType}${jetAlgo}.txt
for j in `ls $maindir/$i/${analyzerType}*${recoType}${jetAlgo}*`
do
echo $j >> files_${analyzerType}_2ndLevel_${i}_${recoType}${jetAlgo}.txt
#echo $j
#echo $maindir/$i/$j 
done
done
