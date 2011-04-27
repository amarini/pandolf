#! /usr/bin/env python
import os
import sys
import re
import time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 6:
    print "usage sendOnBatch.py PDname SDname recoType jetAlgo files_per_job"
    print "example : sendOnBatch.py QCD_Spring10 Pt80 pf akt5 10"
    print "example : sendOnBatch.py EG Run2010A-PromptReco-v4 calo kt6 5"
    sys.exit(1)
PDname = sys.argv[1]
SDname = sys.argv[2]
dataset_name = PDname + "_" + SDname
dataset_path = PDname + "/" + SDname
print "path: " + dataset_path + "  name: " + dataset_name
recoType = sys.argv[3]
jetAlgo = sys.argv[4]
inputlist = "files_"+dataset_name+".txt"
#inputlist = "files_EG_Run2010A_PromptReco_v4.txt"
#inputlist = "files_MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1.txt"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = sys.argv[5]
queue = "cmsshort"
#ijobmax = 40
ijobmax = int(sys.argv[5])
application = "do2ndLevel_PhotonJet_batch"
# to write on the cmst3 cluster disks
################################################
# to write on local disks
################################################
#diskoutputdir = "/cmsrm/pc21_2/pandolf/DATA/EG/Run2010A-PromptReco-v4"
if PDname=="EG" or PDname=="Photon":
  diskoutputdir = "/cmshome/pandolf/STORE/DATA/" + dataset_path
else:
  diskoutputdir = "/cmshome/pandolf/STORE/MC/Spring11/" + dataset_path
diskoutputmain = diskoutputdir
# prepare job to write on the cmst3 cluster disks
################################################
dir = dataset_name+"_"+recoType+jetAlgo
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")
os.system("mkdir -p "+dir+"/res/")

if diskoutputdir != "none": 
    os.system("mkdir -p "+diskoutputmain)


#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
#numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
input = open(inputlist)
inputfiles = input.readlines()
######################################
ijob=0

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+dir+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)

    inputfile.close()


    # prepare the script to run
    outputname = dataset_name+"_"+recoType+jetAlgo+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
    outputfile.write('cd /cmshome/pandolf/CMSSW_4_1_5/src ; eval `scramv1 runtime -sh` ; cd -\n')
    #    outputfile.write('cd '+pwd)
    #outputfile.write('cd $WORKDIR\n')
    outputfile.write("mkdir -p "+dir+"/res/job_"+str(ijob)+"\n")
    outputfile.write("cp Cert* "+dir+"/res/job_"+str(ijob)+"\n")
    #outputfile.write("cp "+application+" "+dir+"/res/job_"+str(ijob)+"\n")
    outputfile.write("cd "+dir+"/res/job_"+str(ijob)+"\n")
    outputfile.write(pwd+'/'+application+" "+dataset_name+" "+recoType+" "+jetAlgo+" "+inputfilename+" _"+str(ijob)+"\n") 
    # select this for GENJETS ntuples:
    #outputfile.write(pwd+'/'+application+" "+dataset_name+" "+recoType+" "+jetAlgo+" "+inputfilename+" _"+str(ijob)+" true\n")
    outputfile.write('mv *.root '+diskoutputmain+'\n') 
    #outputfile.write('ls *.root | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} cmsrm-ui02:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+dir+"/log/"+dataset_name+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dir+"/log/"+dataset_name+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+dataset_name+"_"+str(ijob))
    ijob = ijob+1
    #time.sleep(2.)
    continue
