
#include "Ntp1Finalizer_TTZDilepton.h"
#include "TMath.h"
#include <iostream>



double delta_phi(double phi1, double phi2) {

  double dphi = fabs(phi1 - phi2);
  return (dphi <= TMath::Pi())? dphi : TMath::TwoPi() - dphi;
}


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}





int main( int argc, char* argv[] ) {

  if( argc!=3 && argc!=4 && argc!=5 && argc!=6 ) {
    std::cout << "USAGE: ./finalize_TTZDilepton [dataset] [selectionType] [bTaggerType=\"SSVHE\"] [PUType=\"HR11_73pb\"] [leptType=\"ALL\"]" <<std::endl;
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



  Ntp1Finalizer_TTZDilepton* nf = new Ntp1Finalizer_TTZDilepton( dataset, selectionType, bTaggerType, PUType, leptType );
  nf->set_inputAnalyzerType("TTZ");


  if( dataset=="DATA_HR11_v2" ) {
   
    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("SingleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");
    nf->addFile("DoubleElectron_Run2011B_v2");

  } else {
  
    nf->addFile( dataset );

  }

  nf->finalize();


  return 0;

}


