#include "TreeAnalyzer_Inclusive.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string inputFileList, std::string dataset, std::string recoType, std::string jetAlgo, std::string flags);



int main(int argc, char* argv[]) { 

  if( argc != 6 ) {
    std::cout << "USAGE: ./do2ndLevel_Inclusive_batch [dataset] [recoType] [jetAlgo] [inputFileList] [flags]" << std::endl;
    exit(23);
  }


  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);
  std::string inputFileList(argv[4]);
  std::string flags(argv[5]);
  
  doSingleLoop(inputFileList, dataset, recoType, jetAlgo, flags);

}


void doSingleLoop(std::string inputFileList, std::string dataset, std::string recoType, std::string jetAlgo, std::string flags) {

  std::cout << "Starting: " << dataset << std::endl;
  TreeAnalyzer_Inclusive* t = new TreeAnalyzer_Inclusive(dataset, recoType, jetAlgo, flags);
  t->LoadInputFromFile( inputFileList );
  t->ReadJSONFile("Cert_132440-137028_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt");
  t->ReadCSVFile("lumi_by_LS_upto137028.csv");
  t->Loop();
  delete t;

}
