#include "TreeAnalyzer_PhotonJet.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON, std::string flags, bool useGenJets);



int main(int argc, char* argv[]) { 
//do2ndLevel(std::string dataset, std::string recoType, std::string jetAlgo) {

  if( argc != 4 &&  argc != 5  ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJet [dataset] [recoType] [jetAlgo] [useGenJets=false]" << std::endl;
    exit(23);
  }

  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);

  std::string name;
  bool useGenJets = false;
  if( argc == 5 ) {
    std::string useGenJets_str(argv[4]);
    if( useGenJets_str=="true" )
      useGenJets = true;
  }
    
  if( dataset=="PhotonJet_Summer09" || dataset=="QCD_Summer09") {

    name = dataset + "_Pt15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt30";
    //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt170";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt300";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt470";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt800";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt1400";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);

  } else if( dataset=="PhotonJet_Spring10" ) {
    
    name = dataset + "_Pt5to15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", (bool)false);
    name = dataset + "_Pt15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt30";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt170";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    
  } else if( dataset=="PhotonJet_Summer1036X" ) {
    
    name = dataset + "_Pt5to15";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt15to20";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt20to30";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt30to50";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt50to80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt80to120";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt120to170";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt170to300";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    
  } else if( dataset=="QCD_Spring10" ) {
    
  //name = dataset + "_Pt5to15";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt15to20";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt20to30";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt30to50";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt30to50";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt50to80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt80";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt170";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
    name = dataset + "_Pt300";
    doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt470";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt800";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);
  //name = dataset + "_Pt1400";
  //doSingleLoop(name, recoType, jetAlgo, (bool)false, "", useGenJets);

  } else if( dataset=="MinimumBias_Commissioning10-GOODCOLL-v9" 
          || dataset=="MinimumBias_Commissioning10-Apr20Skim_GOODCOLL-v1" 
          || dataset=="MinimumBias_Commissioning10_May6thPDSkim2_SD_EG" 
          || dataset=="MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1"
          || dataset=="MinimumBias_Commissioning10_SD_EG-v9"
          || dataset=="EG_Run2010A-PromptReco-v1"
          || dataset=="EG_Run2010A-PromptReco-v2"
          || dataset=="EG_Run2010A-PromptReco-v4"
          ) {

    doSingleLoop(dataset, recoType, jetAlgo, (bool)true, "", useGenJets);

  } else { //generic case:

    doSingleLoop(dataset, recoType, jetAlgo, (bool)false, "", useGenJets);

  }

}


void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON, std::string flags, bool useGenJets) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_PhotonJet* t = new TreeAnalyzer_PhotonJet(name.c_str(), recoType, jetAlgo, flags, useGenJets);


  if( useGenJets ) {
    t->SetFlags("_GENJETS");
  }

  if( useJSON ) {
    t->ReadJSONFile("Cert_132440-139459_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt");
    t->ReadCSVFile("lumi_by_LS_upto138939.csv");
  }

  t->LoadInput();
  t->Loop();
  delete t;

}
