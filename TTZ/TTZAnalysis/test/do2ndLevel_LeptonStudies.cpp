#include "Ntp1Analyzer_LeptonStudies.h"
#include <stdlib.h>

#include "TString.h"



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_LeptonStudies [dataset] [inputfile=\"\"] [flags=\"\"]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_LeptonStudies* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_LeptonStudies(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_LeptonStudies(dataset, flags);
  }

  TString dataset_tstr(dataset);


  if( argc==2 ) {
    std::string fileName = "files_" + dataset + ".txt";
    na->LoadInputFromFile(fileName.c_str());
    //na->LoadInput();
  } else {
    std::string inputfile(argv[2]);
    na->LoadInputFromFile(inputfile.c_str());
  }

  na->Loop();

  delete na;

  return 0;

}
