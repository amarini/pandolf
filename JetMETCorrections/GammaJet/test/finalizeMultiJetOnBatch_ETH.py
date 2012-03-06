#!/usr/bin/env python
import os
import sys
import re
import time

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print "usage finalizeMultiJetOnBatch.py dataset nBlocks analyzerType=\"MultiJet\""
    sys.exit(1)
dataset = sys.argv[1]
dataset_path = dataset
nBlocks = int(sys.argv[2])
analyzerType = "MultiJet"
if len(sys.argv)==4 :
  analyzerType = sys.argv[3]

queue = "short.q"

#if analyzerType=="MultiJet_DIJET":
#  analyzerType = "MultiJet"
#  flags


#application = "finalize_" + analyzerType
application = "finalize_MultiJet"


match_Fall11 = re.search( r'Fall11', dataset, re.M|re.I)
isDATA = re.search( r'Run2011', dataset, re.M|re.I)

if isDATA:
  diskoutputdir = "/shome/pandolf/FinalizedTrees/DATA/" + dataset_path
else:
  if match_Fall11:
    diskoutputdir = "/shome/pandolf/FinalizedTrees/MC/Fall11/" + dataset_path
  else:
    diskoutputdir = "/shome/pandolf/FinalizedTrees/MC/Summer11/" + dataset_path

diskoutputmain = diskoutputdir


dir = analyzerType + "_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")

if diskoutputdir != "none": 
    os.system("mkdir -p "+diskoutputmain)


pwd = os.environ['PWD']


ijob=0

while (ijob<nBlocks):

  print( str(ijob) + " / " + str(nBlocks) +"\n\n" )
  WORKDIR=os.getenv("PWD")
  outputname = dir+"/src/submit_"+str(ijob)+".src"
  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
  outputfile.write('source /swshare/cms/cmsset_default.sh\n')
  outputfile.write('cd /shome/pandolf/CMSSW_4_2_8/src/ ; eval `scramv1 runtime -sh` ; cd -\n')
  outputfile.write('cd '+WORKDIR+'\n')
  outputfile.write(pwd+'/'+application+" "+dataset+" "+str(ijob)+" "+ str(nBlocks) + "\n") 
  # select this for GENJETS ntuples:
  #outputfile.write(pwd+'/'+application+" "+dataset+" "+recoType+" "+jetAlgo+" "+inputfilename+" _"+str(ijob)+" true\n")
  outputfile.write('mv ' + analyzerType + "_" + dataset + "_" + str(ijob) + '.root '+diskoutputmain+'\n') 
  outputfile.close()
  #print "chmod +x "+outputname
  os.popen("chmod +x "+outputname)
  #os.system("echo qsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log "+pwd+"/"+outputname)
  #out = os.popen("echo qsub -q "+queue+" -o "+dir+"/log/"+dataset+"_"+str(ijob)+".log "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob))
  os.system("ls -lah "+pwd+"/"+outputname)
  command = "qsub -q "+queue+" -o "+pwd+"/"+dir+"/log/"+dataset+"_"+str(ijob)+".log  -e "+pwd+"/"+dir+"/log/"+dataset+"_"+str(ijob)+".err "+pwd+"/"+outputname+" -copyInput="+dataset+"_"+str(ijob)
  print command
  os.system(command)
  ijob = ijob+1
