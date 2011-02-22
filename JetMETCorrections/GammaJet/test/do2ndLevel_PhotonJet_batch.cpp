#include "TreeAnalyzer_PhotonJet.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include "TRegexp.h"
#include "TString.h"



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
  TreeAnalyzer_PhotonJet* t = new TreeAnalyzer_PhotonJet(name.c_str(), recoType, jetAlgo, flags, useGenJets);
  t->LoadInputFromFile(fileName);
  if( useJSON ) {
    t->ReadJSONFile("Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_CMSSWConfig.txt");
  //t->ReadCSVFile("csv_runs143337_144114.txt");
  //t->ReadJSONFile("Cert_132440-143336_7TeV_StreamExpress_Collisions10_CMSSWConfig_v2.txt");
  //t->ReadCSVFile("csvfile_upto143336.txt");
  }
  t->Loop();
  delete t;

}
