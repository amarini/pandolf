#include "TreeAnalyzer_HZZlljj.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON);



int main(int argc, char* argv[]) { 
//do2ndLevel(std::string dataset, std::string recoType, std::string jetAlgo) {

  if( argc != 4 ) {
    std::cout << "USAGE: ./do2ndLevel_HZZlljj [dataset] [recoType] [jetAlgo]" << std::endl;
    exit(23);
  }

  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);

  std::string name;
    
  doSingleLoop(dataset, recoType, jetAlgo, (bool)false);


}


void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON ) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_HZZlljj* t = new TreeAnalyzer_HZZlljj(name.c_str(), recoType, jetAlgo);
  if( useJSON ) {
    //t->ReadJSONFile("json_L1_HLT_PIX_STRIP_ECAL_HCAL_132440_135735.txt");
    t->ReadJSONFile("Cert_132440-135735_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt");
    t->ReadCSVFile("lumi_by_LS_132440_135735.csv");
  }
  t->Loop();
  delete t;

}
