#include "Ntp1Analyzer_TTZ.h"
#include <stdlib.h>

#include "TString.h"



int main( int argc, char* argv[]) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./do2ndLevel_TTZ [dataset] [inputfile=""] [flags=""]" << std::endl;
    exit(31);
  }

  std::string dataset(argv[1]);

  Ntp1Analyzer_TTZ* na;

  if( argc<4 ) {
    na = new Ntp1Analyzer_TTZ(dataset);
  } else {
    std::string flags(argv[3]);
    na = new Ntp1Analyzer_TTZ(dataset, flags);
  }

  TString dataset_tstr(dataset);

  if( dataset_tstr.Contains("Run2011") || dataset_tstr.Contains("Run2012") ) { //is data

    if( dataset_tstr.BeginsWith("DoubleMu") ) {

      na->AddRequiredTrigger( "HLT_DoubleMu7_v" );
      na->AddRequiredTrigger( "HLT_Mu13_Mu8_v" );
      na->AddRequiredTrigger( "HLT_Mu17_Mu8_v" );

    } else if( dataset_tstr.BeginsWith("DoubleElectron") ) {

      na->AddRequiredTrigger( "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v" );
      na->AddRequiredTrigger( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v" );
      na->AddRequiredTrigger( "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v" );
 
    } else if( dataset_tstr.BeginsWith("MuEG") ) {

      na->AddRequiredTrigger( "HLT_Mu17_Ele8_CaloIdL_v" );
      na->AddRequiredTrigger( "HLT_Mu8_Ele17_CaloIdL_v" );
      na->AddRequiredTrigger( "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v" );

    } else if( dataset_tstr.BeginsWith("SingleMu") ) {

      na->AddRequiredTrigger( "HLT_IsoMu24_v" );

    }

  }  //if is data



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
