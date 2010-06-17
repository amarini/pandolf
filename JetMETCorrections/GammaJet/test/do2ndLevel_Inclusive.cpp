#include "TreeAnalyzer_Inclusive.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo);



int main(int argc, char* argv[]) { 
//do2ndLevel(std::string dataset, std::string recoType, std::string jetAlgo) {

  if( argc != 4 ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJet [dataset] [recoType] [jetAlgo]" << std::endl;
    exit(23);
  }

  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);


  doSingleLoop(dataset, recoType, jetAlgo);


}


void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_Inclusive* t = new TreeAnalyzer_Inclusive(name, recoType, jetAlgo);
  t->ReadJSONFile("json_L1_HLT_PIX_STRIP_ECAL_HCAL_132440_134725.txt");
  t->ReadCSVFile("lumi_by_LS_all.csv");
  t->Loop();
  delete t;

}
