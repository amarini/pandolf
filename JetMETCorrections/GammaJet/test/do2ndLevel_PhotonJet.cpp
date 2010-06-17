#include "TreeAnalyzer_PhotonJet.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON, bool useGenJets);



int main(int argc, char* argv[]) { 
//do2ndLevel(std::string dataset, std::string recoType, std::string jetAlgo) {

  if( argc != 4 ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJet [dataset] [recoType] [jetAlgo]" << std::endl;
    exit(23);
  }

  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);

  std::string name;
  bool useGenJets = true;
    
    
  if( dataset=="PhotonJet_Summer09" || dataset=="QCD_Summer09") {

    name = dataset + "_Pt15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt30";
    //doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt170";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt300";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt470";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt800";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt1400";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);

  } else if( dataset=="PhotonJet_Spring10" ) {
    
    name = dataset + "_Pt0to15";
//    doSingleLoop(name, recoType, jetAlgo, (bool)false, (bool)false);
    name = dataset + "_Pt15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, (bool)false);
    name = dataset + "_Pt30";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, (bool)false);
    name = dataset + "_Pt80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, (bool)false);
    
  } else if( dataset=="QCD_Spring10" ) {
    
    name = dataset + "_Pt0to15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt15to20";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt20to30";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt30to50";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt30to50";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt50to80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt170";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt300";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt470";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt800";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);
    name = dataset + "_Pt1400";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, useGenJets);

  } else if( dataset=="MinimumBias_Commissioning10-GOODCOLL-v9" 
          || dataset=="MinimumBias_Commissioning10-Apr20Skim_GOODCOLL-v1" 
          || dataset=="MinimumBias_Commissioning10_May6thPDSkim2_SD_EG" 
          || dataset=="MinimumBias_Commissioning10_SD_EG-v9"
          || dataset=="EG_Run2010A-PromptReco-v1"
          || dataset=="EG_Run2010A-PromptReco-v2"
          ) {

    doSingleLoop(dataset, recoType, jetAlgo, (bool)true, (bool)false);

  } else { //generic case:

    doSingleLoop(dataset, recoType, jetAlgo, (bool)false, (bool)false);

  }

}


void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON, bool useGenJets) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_PhotonJet* t = new TreeAnalyzer_PhotonJet(name.c_str(), recoType, jetAlgo, useGenJets);
  if( useJSON ) {
    //t->ReadJSONFile("json_L1_HLT_PIX_STRIP_ECAL_HCAL_132440_135735.txt");
    t->ReadJSONFile("Cert_132440-136100_7TeV_StreamExpress_Collisions10_CMSSWConfig_PRIVATE.txt");
    t->ReadCSVFile("lumi_by_LS_132440_136100.csv");
  }
  t->Loop();
  delete t;

}
