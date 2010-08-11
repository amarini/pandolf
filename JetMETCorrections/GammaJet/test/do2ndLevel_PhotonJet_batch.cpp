#include "TreeAnalyzer_PhotonJet.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string fileName, std::string name, std::string recoType, std::string jetAlgo, std::string flags, bool useJSON, bool useGenJets);



int main(int argc, char* argv[]) { 

  if( argc != 6 && argc != 7 ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJet_batch [dataset] [recoType] [jetAlgo] [inputFileList] [flags] [useGenJets=false]" << std::endl;
    exit(23);
  }


  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);
  std::string inputFileList(argv[4]);
  std::string flags(argv[5]);
  bool useGenJets = false;

  std::cout << "argv[1]: " << argv[1] << std::endl;
  std::cout << "argv[2]: " << argv[2] << std::endl;
  std::cout << "argv[3]: " << argv[3] << std::endl;
  std::cout << "argv[4]: " << argv[4] << std::endl;
  std::cout << "argv[5]: " << argv[5] << std::endl;

  if( argc == 7 ) {
    std::string useGenJets_str(argv[6]);
    if( useGenJets_str=="true" ) useGenJets = true;
  }
  
  if( dataset=="EG_Run2010A-PromptReco-v1" ||
      dataset=="EG_Run2010A-PromptReco-v2" ||
      dataset=="EG_Run2010A-PromptReco-v4" ||
      dataset=="EG_Run2010A_PromptReco_v4" ||
      dataset=="EG_Run2010A_PromptReco_v4_80MeV_MET" ||
      dataset=="EG_Run2010A_PromptReco_v4_139347_139375" ||
      dataset=="EG_Run2010A_PromptReco_v4_139376_139459" ||
      dataset=="EG_Run2010A_Jun14thReReco_v1" ||
      dataset=="EG_Run2010A_Jun14thReReco_v1_80MeV_MET" ||
      dataset=="EG_Run2010A_Jul6thReReco_v1" ||
      dataset=="EG_Run2010A_Jul15thReReco_v1" ||
      dataset=="EG_Run2010A_Jul26thReReco_v1" ||
      dataset=="MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1" || 
      dataset=="MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1_80MeV"  ||
      dataset=="MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1_80MeV_MET" ) {
    doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags, (bool)true, (bool)false);
  } else {
    doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags, (bool)false, useGenJets);
  }


}


void doSingleLoop(std::string fileName, std::string name, std::string recoType, std::string jetAlgo, std::string flags, bool useJSON, bool useGenJets) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_PhotonJet* t = new TreeAnalyzer_PhotonJet(name.c_str(), recoType, jetAlgo, flags, useGenJets);
  t->LoadInputFromFile(fileName);
  if( useJSON ) {
    t->ReadJSONFile("Cert_132440-140399_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt");
    t->ReadCSVFile("lumi_by_LS_132440_140401.csv");
  }
  t->Loop();
  delete t;

}
