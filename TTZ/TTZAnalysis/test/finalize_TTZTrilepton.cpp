
#include "Ntp1Finalizer_TTZTrilepton.h"
#include "TMath.h"
#include <iostream>








int main( int argc, char* argv[] ) {

  if( argc!=3 && argc!=4 && argc!=5 && argc!=6 ) {
    std::cout << "USAGE: ./finalize_TTZTrilepton [dataset] [selectionType] [bTaggerType=\"SSVHE\"] [PUType=\"HR11_73pb\"] [leptType=\"ALL\"]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);

  std::string bTaggerType="SSVHE";
  if( argc==4 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }

  std::string PUType="HR11_73pb";
  if( argc==5 ) {
    std::string PUType_str(argv[4]);
    PUType = PUType_str;
  }

  std::string leptType="ALL";
  if( argc==6 ) {
    std::string leptType_str(argv[5]);
    leptType = leptType_str;
  }



  Ntp1Finalizer_TTZTrilepton* nf = new Ntp1Finalizer_TTZTrilepton( dataset, selectionType, bTaggerType, PUType, leptType );
  nf->set_inputAnalyzerType("TTZ");


  if( dataset=="DATA_Run2011_FULL" ) {
   
    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");
    nf->addFile("DoubleElectron_Run2011B_v2");
    nf->addFile("MuEG_Run2011A_FULL");
    nf->addFile("MuEG_Run2011B");

  } else if( dataset=="VV_Summer11" ) {

    nf->addFile("WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1");
    nf->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");

  } else {
  
    nf->addFile( dataset );

  }

  nf->finalize();


  return 0;

}


