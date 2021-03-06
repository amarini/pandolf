#! /bin/bash

myrand=$1
mass=$2
OUTDIR=datacards
echo "Starting HiggsCombination with seed=$myrand at $( date +%c ) on $hostname."

startdir=$( pwd )

#set CMSSW environment
RELDIR=/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/

#algo="MarkovChainMC"
#algo="HybridNew"
algo="ProfileLikelihood"
hint="ProfileLikelihood" # before the algo method, run the hint method for restricting integration field
label="lvqq"
ntoys=100
WORKDIR=${RELDIR}/UserCode/pandolf/HWWlvjj/HWWlvjjAnalyzer/FIT/${OUTDIR}/${mass}/
datacard="CMS_hwwlvqq_${mass}_2channels.txt" 
OUTDIR="combine_${label}_${algo}_"$( basename $datacard .txt )

cd $RELDIR
export SCRAM_ARCH=slc5_amd64_gcc434
#cmsenv
eval `scramv1 runtime -sh`
cd $curdir

TMPDIR="/tmp/$(whoami)"
mkdir ${TMPDIR}/combine_${myrand}
cd $TMPDIR/combine_${myrand}
cp $WORKDIR/*input*root .
echo "I am in $( pwd ) (it should be: $TMPDIR/combine_${myrand} )"
echo


if [ ! -d ${WORKDIR}/$OUTDIR/ ]
    then
    mkdir ${WORKDIR}/$OUTDIR/
fi
echo "Datacard: $datacard"
# if algo=HybridNew
#combine -M $algo -n $label -m 400 -s $myrand -t $ntoys -U  -d $WORKDIR/$datacard --freq --singlePoint 1

#if algo="MarkovChainMC"
#expected
#combine -M $algo -n $label -m $mass  -s $myrand -d $WORKDIR/$datacard  -H $hint  -t $ntoys -U 
#observed
#combine -M $algo -n $label -m $mass  -s $myrand -d $WORKDIR/$datacard  -H $hint  -U

#if algo="MarkovChainMC" and function RooCB
#expected
#combine -M $algo -n $label -m $mass  -s $myrand -d $WORKDIR/$datacard  -H $hint  -t $ntoys -U -L /cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/HWWlvjj/HWWlvjjAnalyzer/FIT/PDFs/RooCB_cc.so 
#combine -M $algo -n $label -m $mass  -s $myrand -d $WORKDIR/$datacard  -H $hint  -t $ntoys -U 
#observed
echo "combine -M $algo -n $label -m $mass  -s $myrand -d $WORKDIR/$datacard  -H $hint  -U "
combine -M $algo -n $label -m $mass  -s $myrand -d $WORKDIR/$datacard  -H $hint  -U 


#if algo="ProfileLikelihood"
#combine -M $algo   -n $label -m $mass  -s $myrand  -d $WORKDIR/$datacard -U  -t $ntoys

echo "List of files in $( pwd ):"
ls -lh

mv $TMPDIR/combine_${myrand}/higgsCombine${label}*.${myrand}.root  ${WORKDIR}/$OUTDIR/
###mv $TMPDIR/combine_${myrand}/log_combine_${label}_${mass}.${myrand}.out  ${WORKDIR}/$OUTDIR/
