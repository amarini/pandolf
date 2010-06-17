#include "TreeAnalyzer_DiJet.h"
#include <iostream>
#include <string>
#include <cstdlib>



void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON=false);



int main(int argc, char* argv[]) { 
//do2ndLevel(std::string dataset, std::string recoType, std::string jetAlgo) {

  if( argc != 4 ) {
    std::cout << "USAGE: ./do2ndLevel_PhotonJet [dataset] [recoType] [jetAlgo]" << std::endl;
    exit(23);
  }

  std::string dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);

  if( dataset=="QCD_Spring10" ) {

    doSingleLoop("QCD_Spring10_Pt0to15", recoType, jetAlgo);
    doSingleLoop("QCD_Spring10_Pt15", recoType, jetAlgo);
    doSingleLoop("QCD_Spring10_Pt30", recoType, jetAlgo);
    doSingleLoop("QCD_Spring10_Pt80", recoType, jetAlgo);

  } else if( dataset=="MinimumBias_Commissioning10-GOODCOLL-v9" ) {

    doSingleLoop( dataset, recoType, jetAlgo, (bool)true );

  } else {

    doSingleLoop(dataset, recoType, jetAlgo);

  }


}


void doSingleLoop(std::string name, std::string recoType, std::string jetAlgo, bool useJSON) {

  std::cout << "Starting: " << name << std::endl;
  TreeAnalyzer_DiJet* t = new TreeAnalyzer_DiJet(name, recoType, jetAlgo);
  if( useJSON ) {
    t->ReadJSONFile("Cert_132440-135735_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt");
    t->ReadCSVFile("lumi_by_LS_132440_136100.csv");
  }
  t->Loop();
  delete t;

}
