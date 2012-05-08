
#include "Ntp1Finalizer_TTZTrilepton.h"
#include "TMath.h"
#include <iostream>



void runOnAllBackgrounds( Ntp1Finalizer_TTZTrilepton* nf );
void runOnSingleBackground( Ntp1Finalizer_TTZTrilepton* nf, const std::string& dataset );




int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./finalize_TTZTrilepton [selectionType] [bTaggerType=\"TCHE\"] [leptType=\"ALL\"]" <<std::endl;
    return 13;
  }


  std::string selectionType(argv[1]);

  std::string bTaggerType="TCHE";
  if( argc==3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }


  std::string leptType="ALL";
  if( argc==4 ) {
    std::string leptType_str(argv[3]);
    leptType = leptType_str;
  }



  Ntp1Finalizer_TTZTrilepton* nf = new Ntp1Finalizer_TTZTrilepton( "syst", selectionType, bTaggerType, leptType );
  nf->set_inputAnalyzerType("TTZ");

  std::cout << std::endl << "++++++++++++++++++++++" << std::endl;
  std::cout << "+++ BTag Syst +1 Sigma" << std::endl;
  std::cout << "++++++++++++++++++++++" << std::endl << std::endl;
  nf->set_btagSyst(1);

  runOnAllBackgrounds( nf );


  std::cout << std::endl << "++++++++++++++++++++++" << std::endl;
  std::cout << "+++ BTag Syst -1 Sigma" << std::endl;
  std::cout << "++++++++++++++++++++++" << std::endl << std::endl;
  nf->set_btagSyst(-1);

  runOnAllBackgrounds( nf );

  nf->set_btagSyst(0);


  std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
  std::cout << "+++ JES Syst +1 Sigma" << std::endl;
  std::cout << "+++++++++++++++++++++" << std::endl << std::endl;
  nf->set_jes(1);

  runOnAllBackgrounds( nf );

  std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
  std::cout << "+++ JES Syst -1 Sigma" << std::endl;
  std::cout << "+++++++++++++++++++++" << std::endl << std::endl;
  nf->set_jes(-1);

  runOnAllBackgrounds( nf );

  return 0;

}




void runOnAllBackgrounds( Ntp1Finalizer_TTZTrilepton* nf ) {

  
  runOnSingleBackground( nf, "TTJ_Fall11_highstat" );
  runOnSingleBackground( nf, "DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11" );
  runOnSingleBackground( nf, "VV_Summer11" );

}


void runOnSingleBackground( Ntp1Finalizer_TTZTrilepton* nf, const std::string& dataset ) {

  Ntp1Finalizer_TTZTrilepton* nf_syst = new Ntp1Finalizer_TTZTrilepton( *nf );

  nf_syst->addFile( dataset );

  nf_syst->finalize();

  delete nf;

} 
