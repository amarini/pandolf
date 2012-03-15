#!/usr/bin/env python
import os
import sys
import re
import time

if len(sys.argv) != 3 :
    print "usage finalizeQGStudiesOnBatch.py dataset nBlocks"
    sys.exit(1)
dataset = sys.argv[1]
dataset_path = dataset
nBlocks = int(sys.argv[2])
analyzerType = "QGStudies"

queue = "8nh"



#application = "finalize_" + analyzerType
application = "finalize_QGStudies"


match_Fall11 = re.search( r'Fall11', dataset, re.M|re.I)
isDATA = re.search( r'Run2011', dataset, re.M|re.I)

if isDATA:
  diskoutputdir = "/cmsrm/pc24_2/pandolf/FinalizedTrees/DATA/" + dataset_path
else:
  if match_Fall11:
    diskoutputdir = "/cmsrm/pc24_2/pandolf/FinalizedTrees/MC/Fall11/" + dataset_path
  else:
    diskoutputdir = "/cmsrm/pc24_2/pandolf/FinalizedTrees/MC/Summer11/" + dataset_path

diskoutputmain = diskoutputdir


dir = "batch_" + analyzerType + "_" + dataset
os.system("mkdir -p "+dir)
os.system("mkdir -p "+dir+"/log/")
os.system("mkdir -p "+dir+"/input/")
os.system("mkdir -p "+dir+"/src/")

if diskoutputdir != "none": 
    os.system("mkdir -p "+diskoutputmain)

os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm25 mkdir -p "+diskoutputmain)

pwd = os.environ['PWD']


ijob=0

while (ijob<nBlocks):

  print( str(ijob) + " / " + str(nBlocks) +"\n\n" )
  WORKDIR=os.getenv("PWD")
  outputname = dir+"/src/submit_"+str(ijob)+".src"
  outputfile = open(outputname,'w')
  outputfile.write('#!/bin/bash\n')
  outputfile.write('export STAGE_HOST=castorcms\n')
  outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
  outputfile.write('cd /afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_3_patch5/src/ ; eval `scramv1 runtime -sh` ; cd -\n')
  outputfile.write('cd '+WORKDIR+'\n')
  outputfile.write(pwd+'/'+application+" "+dataset+" " + str(ijob)+" "+ str(nBlocks) + "\n") 
  # select this for GENJETS ntuples:
  #outputfile.write(pwd+'/'+application+" "+dataset+" "+recoType+" "+jetAlgo+" "+inputfilename+" _"+str(ijob)+" true\n")
  outputfile.write('mv /tmp/pandolf/' + analyzerType + "_" + dataset + "_pfakt5_" + str(ijob) + '.root '+diskoutputmain+'\n') 
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
