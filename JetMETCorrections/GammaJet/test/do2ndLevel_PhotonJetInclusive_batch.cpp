#include "TreeAnalyzer_PhotonJetInclusive.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include "TRegexp.h"
#include "TString.h"



void doSingleLoop(std::string fileName, std::string name, std::string recoType, std::string jetAlgo, std::string flags, bool useJSON, bool useGenJets);



int main(int argc, char* argv[]) { 

  if( argc != 4 && argc != 5 ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJetInclusive_batch [dataset] [inputFileList] [flags] [useGenJets=false]" << std::endl;
    exit(23);
  }


  std::string dataset(argv[1]);
  std::string recoType = "pf";
  std::string jetAlgo = "akt5";
  std::string inputFileList(argv[2]);
  std::string flags(argv[3]);
  bool useGenJets = false;

  if( argc == 5 ) {
    std::string useGenJets_str(argv[4]);
    if( useGenJets_str=="true" ) useGenJets = true;
  }

  TRegexp run2010("Run2010");
  TRegexp minbias("MinimumBias_Commissioning10");
  TRegexp nov4("Nov4ReReco");
  TRegexp dec22("Dec22ReReco");
  TString dataset_str(dataset);
  
  if( dataset_str.Contains(run2010) || dataset_str.Contains(minbias) || dataset_str.Contains(nov4) || dataset_str.Contains(dec22) ) { // then it's data
    doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags, (bool)true, (bool)false);
  } else {
    doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags, (bool)false, useGenJets);
  }


}


void doSingleLoop(std::string fileName, std::string name, std::string recoType, std::string jetAlgo, std::string flags, bool useJSON, bool useGenJets) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_PhotonJetInclusive* t = new TreeAnalyzer_PhotonJetInclusive(name.c_str(), recoType, jetAlgo, flags, useGenJets);
  t->LoadInputFromFile(fileName);
  if( useJSON ) {
    t->ReadJSONFile("Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_CMSSWConfig.txt");
  }
  t->Loop();
  delete t;

}
