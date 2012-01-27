#include "TreeAnalyzer_MultiJet.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include "TRegexp.h"
#include "TString.h"



void doSingleLoop(std::string fileName, std::string name, std::string recoType, std::string jetAlgo, std::string flags, bool useJSON);



int main(int argc, char* argv[]) { 

  if( argc != 4 ) {
    std::cout << "USAGE: ./do2ndLevel_MultiJet_batch [dataset] [inputFileList] [flags]" << std::endl;
    exit(23);
  }


  std::string dataset(argv[1]);
  std::string recoType = "pf";
  std::string jetAlgo = "akt5";
  std::string inputFileList(argv[2]);
  std::string flags(argv[3]);

  TRegexp run2010("Run2010");
  TRegexp run2011("Run2011");
  TRegexp minbias("MinimumBias_Commissioning10");
  TRegexp nov4("Nov4ReReco");
  TRegexp dec22("Dec22ReReco");
  TString dataset_str(dataset);
  
  if( dataset_str.Contains(run2010) || dataset_str.Contains(run2011) ||dataset_str.Contains(minbias) || dataset_str.Contains(nov4) || dataset_str.Contains(dec22) ) { // then it's data
    doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags, (bool)true);
  } else {
    doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags, (bool)false);
  }


}


void doSingleLoop(std::string fileName, std::string name, std::string recoType, std::string jetAlgo, std::string flags, bool useJSON) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_MultiJet* t = new TreeAnalyzer_MultiJet(name.c_str(), recoType, jetAlgo, flags);
  t->LoadInputFromFile(fileName);
  if( useJSON ) {
//    t->ReadJSONFile("Cert_136033-149442_7TeV_Nov4ReReco_Collisions10_CMSSWConfig.txt");
//    t->ReadJSONFile("Cert_160404-171116_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
    t->ReadJSONFile("Cert_160404-180252_7TeV_PromptReco_Collisions11_CMSSWConfig.txt");
  //t->ReadCSVFile("csv_runs143337_144114.txt");
  //t->ReadJSONFile("Cert_132440-143336_7TeV_StreamExpress_Collisions10_CMSSWConfig_v2.txt");
  //t->ReadCSVFile("csvfile_upto143336.txt");
  }
  t->Loop();
  delete t;

}
