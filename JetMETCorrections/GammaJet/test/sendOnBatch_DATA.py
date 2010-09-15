#! /usr/bin/env python
import os
import sys
import time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 5:
    print "usage sendOnBatch.py dataset recoType jetAlgo njobs"
    sys.exit(1)
dataset = sys.argv[1]
recoType = sys.argv[2]
jetAlgo = sys.argv[3]
inputlist = "files_"+dataset+".txt"
#inputlist = "files_EG_Run2010A_PromptReco_v4.txt"
#inputlist = "files_MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1.txt"
#settingfile = "config/RSZZsettings.txt"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = sys.argv[5]
queue = "8nh"
#ijobmax = 40
ijobmax = int(sys.argv[4])
application = "do2ndLevel_PhotonJet_batch"
# to write on the cmst3 cluster disks
################################################
# to write on local disks
################################################
#diskoutputdir = "/cmsrm/pc21_2/pandolf/DATA/EG/Run2010A-PromptReco-v4"
diskoutputdir = "/cmsrm/pc21_2/pandolf/DATA/EG/Run2010A-PromptReco-v4_runs143337-144114"
diskoutputmain = diskoutputdir
# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir -p "+dataset+"_"+recoType+jetAlgo)
os.system("mkdir -p "+dataset+"_"+recoType+jetAlgo+"/log/")
os.system("mkdir -p "+dataset+"_"+recoType+jetAlgo+"/input/")
os.system("mkdir -p "+dataset+"_"+recoType+jetAlgo+"/src/")

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = pwd+"/"+dataset+"_"+recoType+jetAlgo+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    if ijob != (ijobmax-1):
        for line in range(filesperjob):
            ntpfile = input.readline() 
            inputfile.write(ntpfile)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            inputfile.write(ntpfile)
            ntpfile = input.readline()
            continue
    inputfile.close()

    # prepare the script to run
    outputname = dataset+"_"+recoType+jetAlgo+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp '+pwd+'/Cert_143337-144114_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt $WORKDIR\n')
    outputfile.write('cp '+pwd+'/csv_runs143337_144114.txt $WORKDIR\n')
    #outputfile.write('cp '+pwd+'/Cert_132440-143336_7TeV_StreamExpress_Collisions10_CMSSWConfig_v2.txt $WORKDIR\n')
    #outputfile.write('cp '+pwd+'/csvfile_upto143336.txt $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    outputfile.write(pwd+'/'+application+" "+dataset+" "+recoType+" "+jetAlgo+" "+inputfilename+" _"+str(ijob)+"\n")
    outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm21:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dataset+"_"+recoType+jetAlgo+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dataset+"_"+recoType+jetAlgo+"/log/"+dataset+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
    ijob = ijob+1
    time.sleep(5.)
    continue
