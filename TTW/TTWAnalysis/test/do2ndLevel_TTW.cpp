#include "Ntp1Analyzer_TTW.h"
#include <stdlib.h>



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_TTW [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_TTW* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_TTW(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_TTW(dataset, flags);
  }


 // na->AddRequiredTrigger( "HLT_DoubleMu7" );
 // na->AddRequiredTrigger( "HLT_Mu13_Mu8" );
 // na->AddRequiredTrigger( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL" );


  if( argc==2 ) {
    na->LoadInput();
  } else {
    std::string inputfile(argv[2]);
    na->LoadInputFromFile(inputfile.c_str());
  }

  na->Loop();

  delete na;

  return 0;

}
