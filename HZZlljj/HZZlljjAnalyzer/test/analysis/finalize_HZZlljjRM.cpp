
#include "Ntp1Finalizer_HZZlljjRM.h"
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

  if( argc!=3 && argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./finalize_HZZlljjRM [dataset] [selectionType] [PUType=\"HR11\"] [leptType=\"ALL\"]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);

  std::string PUType="HR11_73pb";
  if( argc>3 ) {
    std::string PUType_str(argv[3]);
    PUType = PUType_str;
  }

  std::string leptType="ALL";
  if( argc>4 ) {
    std::string leptType_str(argv[4]);
    leptType = leptType_str;
  }



  Ntp1Finalizer_HZZlljjRM* nf = new Ntp1Finalizer_HZZlljjRM( dataset, selectionType, PUType, leptType );
  nf->set_inputAnalyzerType("HZZlljj");

  if( dataset=="DATA_EG_37X" ) {

    nf->addFile( "EG_Run2010A_Jul15thReReco_v1" );
    nf->addFile( "EG_Run2010A_Jul26thReReco_v1" );

  } else if( dataset=="Electron_Run2010B" ) {

    nf->addFile( "Electron_Run2010B_PromptReco_v2_runs146240_146733_prod2" );
    nf->addFile( "Electron_Run2010B_PromptReco_v2_runs146733_147111" );

  } else if( dataset=="Run2010B_runs146240_146733" ) {

    nf->addFile( "Electron_Run2010B_PromptReco_v2_runs146240_146733_prod2" );
    nf->addFile( "MU_Run2010B_PromptReco_v2_runs146240_146733" );

  } else if( dataset=="Mu_upto147589" ) {

    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146734_147589");
    nf->addFile("Mu_Run2010A-Sep17ReReco_v2_runs135821_144114");
    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146240_146733");

  } else if( dataset=="EleMu_Nov4ReReco_PU" ) {

    nf->addFile("Mu_Nov4ReReco_PU");
    nf->addFile("Electron_Nov4ReReco_PU");

  } else if( dataset=="DATA_Run2011A" ) {

    nf->addFile("DoubleMu_Run2011A"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A");

  } else if( dataset=="DATA_Run2011A_v2_Sub2" ) {

    nf->addFile("DoubleMu_Run2011A_v2_Sub2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_v2_Sub2");

  } else if( dataset=="DATA_175pb" ) {

    nf->addFile("DoubleMu_1fb"); //first muons! important!
    nf->addFile("DoubleElectron_1fb");


  } else if( dataset=="DATA_1fb" ) {

    nf->addFile("DoubleMu_1fb"); //first muons! important!
    nf->addFile("DoubleElectron_1fb");

  } else if( dataset=="DATA_EPS_FINAL" ) {

    nf->addFile("DoubleMu_EPS_FINAL"); //first muons! important!
    nf->addFile("DoubleElectron_EPS_FINAL");

  } else if( dataset=="DATA_EPS_FINAL_FULL" ) {

    nf->addFile("DoubleMu_EPS_FINAL_FULL"); //first muons! important!
    nf->addFile("SingleMu_EPS_FINAL_FULL"); //first muons! important!
    nf->addFile("DoubleElectron_EPS_FINAL_FULL");

  } else if( dataset=="DATA_EPS_FINAL_FULL_plusSingleMu" ) {

    nf->addFile("DoubleMu_EPS_FINAL_FULL"); //first muons! important!
    nf->addFile("SingleMu_EPS_FINAL_FULL", "passed_HLT_IsoMu4 && !passed_HLT_DoubleMu7 && !passed_HLT_Mu13_Mu8" ); //first muons! important!
    nf->addFile("DoubleElectron_EPS_FINAL_FULL");

  } else if( dataset=="DATA_LP11" ) {

    nf->addFile("DoubleMu_EPS_FINAL_FULL"); //first muons! important!
    nf->addFile("SingleMu_EPS_FINAL_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Aug05ReReco"); //first muons! important!
    nf->addFile("SingleMu_Aug05ReReco"); //first muons! important!
    nf->addFile("DoubleMu_PromptReco_v6"); //first muons! important!
    nf->addFile("SingleMu_PromptReco_v6"); //first muons! important!
    nf->addFile("DoubleElectron_EPS_FINAL_FULL");
    nf->addFile("DoubleElectron_Aug05ReReco");
    nf->addFile("DoubleElectron_PromptReco_v6");

  } else if( dataset=="DATA_Run2011A_FULL" ) {

    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("SingleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");

  } else if( dataset=="DATA_Run2011B_v1" ) {

    nf->addFile("DoubleMu_Run2011B_v1"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_v1"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011B_v1");

  } else if( dataset=="DATA_Run2011B_Oct21JSON" ) {

    nf->addFile("DoubleMu_Run2011B_Oct21JSON"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_Oct21JSON"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011B_Oct21JSON");

  } else if( dataset=="DATA_HR11" ) {

    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Run2011B_v1"); //first muons! important!
    nf->addFile("SingleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_v1"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");
    nf->addFile("DoubleElectron_Run2011B_v1");

  } else if( dataset=="DATA_HR11_v2" ) {

    nf->addFile("DoubleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("DoubleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("SingleMu_Run2011A_FULL"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011A_FULL");
    nf->addFile("DoubleElectron_Run2011B_v2");

  } else if( dataset=="DATA_Run2011B_v2" ) {

    nf->addFile("DoubleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("SingleMu_Run2011B_v2"); //first muons! important!
    nf->addFile("DoubleElectron_Run2011B_v2");

  } else if( dataset=="DATA_EPS_FINAL_plusSingleMu" ) {

    nf->addFile("DoubleMu_EPS_FINAL"); //first muons! important!
    nf->addFile("SingleMu_EPS_FINAL_FULL"); //first muons! important!
    nf->addFile("DoubleElectron_EPS_FINAL");

  } else if( dataset=="DATA_EPS" ) {

    nf->addFile("DoubleMu_EPS"); //first muons! important!
    nf->addFile("DoubleElectron_EPS");

  } else if( dataset=="DATA_postEPS" ) {

    nf->addFile("DoubleMu_postEPS"); //first muons! important!
    nf->addFile("DoubleElectron_postEPS");

  } else if( dataset=="ZJets_alpgen_TuneZ2_Fall10" ) {

    nf->addFile( "Z0Jets_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z1Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_3" );
    nf->addFile( "Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_3" );

  } else if( dataset=="ZJets_alpgen_TuneZ2_Fall10_v2" ) {

    nf->addFile( "Z0Jets_TuneZ2_7TeV-alpgen-tauola" );
    nf->addFile( "Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z1Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_v2" );
    nf->addFile( "Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_v2" );

  } else if( dataset=="ZJets_alpgen_TuneZ2_Spring11" ) {

    nf->addFile( "Z0Jets_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3" );
    nf->addFile( "Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_4" );
    nf->addFile( "Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3" );
    nf->addFile( "Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3" );
    nf->addFile( "Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );

  } else if( dataset=="ZJets_alpgen_TuneZ2_Spring11_v2" ) {

    nf->addFile( "Z0Jets_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z1Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z1Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z1Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z2Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z2Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z2Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z2Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z3Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z3Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z3Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z3Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z4Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z4Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z4Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z4Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z5Jets_ptZ-0to100_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z5Jets_ptZ-100to300_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z5Jets_ptZ-300to800_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "Z5Jets_ptZ-800to1600_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );

  } else if( dataset=="ZBB_alpgen_TuneZ2_Spring11_v2" ) {

    nf->addFile( "ZBB0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "ZBB1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "ZBB2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "ZBB3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );

  } else if( dataset=="ZCC_alpgen_TuneZ2_Spring11_v2" ) {

    nf->addFile( "ZCC0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "ZCC1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "ZCC2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "ZCC3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );

  } else if( dataset=="ZBB_alpgen_TuneZ2_Spring11" ) {

    nf->addFile( "ZBB0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3" );
    nf->addFile( "ZBB1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "ZBB2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "ZBB3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="ZCC_alpgen_TuneZ2_Spring11" ) {

    nf->addFile( "ZCC0JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "ZCC1JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "ZCC2JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1" );
    nf->addFile( "ZCC3JetsToLNu_TuneZ2_7TeV-alpgen-tauola_Spring11-PU_S1_START311_V1G1-v1_3" );

  } else if( dataset=="VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" ) {

    nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" );
    nf->addFile( "WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" );
    nf->addFile( "WWtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10" );

  } else if( dataset=="VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2" ) {

    nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "WWtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );

  } else if( dataset=="VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_OLDPU" ) {

    nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_OLDPU" );
    nf->addFile( "WWtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_OLDPU" );
    nf->addFile( "WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_OLDPU" );

  } else if( dataset=="VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1" ) {

    //nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1" );
    nf->addFile( "ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1" );
    nf->addFile( "WZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1" );
    nf->addFile( "WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1" );

  } else if( dataset=="TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2" ) {

    nf->addFile( "TTTo2L2Nu2B_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );
    nf->addFile( "TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2" );

  } else if( dataset=="TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11_v2_OLDPU" ) {

    nf->addFile( "TTTo2L2Nu2B_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_OLDPU" );
    nf->addFile( "TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_Spring11_v2_OLDPU" );

  } else if( dataset=="TT_TW_TuneZ2_7TeV-pythia6-tauola_Spring11" ) {

    //nf->addFile( "TT_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_3" );
    nf->addFile( "TTTo2L2Nu2B_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1" ) {

    nf->addFile( "TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1" ) ;
    nf->addFile( "T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1" );
    nf->addFile( "Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1" );

  } else if( dataset=="all" ) {


    if( leptType=="ELE" ) {
      argv[1] = new char[strlen("EG_upto146724")];
      strcpy( argv[1], "EG_upto146724");
      main( argc, argv );
    } else if( leptType=="MU") {
      argv[1] = new char[strlen("Mu_upto147589")];
      strcpy( argv[1], "Mu_upto147589");
      main( argc, argv );
      //finalize( "Mu_upto147589", leptType);
    } else if( leptType=="ALL" ) {
      argv[1] = new char[strlen("ELEMUCombined")];
      strcpy( argv[1], "ELEMUCombined");
      main( argc, argv );
    }


    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M130")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M130");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M150")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M150");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M200")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M200");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M300")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M300");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M400")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M400");
    main( argc, argv );
    argv[1] = new char[strlen("HZZ_qqll_gluonfusion_M500")];
    strcpy( argv[1], "HZZ_qqll_gluonfusion_M500");
    main( argc, argv );
    argv[1] = new char[strlen("TTbar_2l_Spring10")];
    strcpy( argv[1], "TTbar_2l_Spring10");
    main( argc, argv );
    argv[1] = new char[strlen("ZZ_Spring10")];
    strcpy( argv[1], "ZZ_Spring10");
    main( argc, argv );
    argv[1] = new char[strlen("ZJets_madgraph")];
    strcpy( argv[1], "ZJets_madgraph");
    main( argc, argv );
    argv[1] = new char[strlen("ZJets_alpgen")];
    strcpy( argv[1], "ZJets_alpgen");
    main( argc, argv );
    //return 0;
    exit(1111);

  } else {
  
    nf->addFile( dataset );

  }

  nf->finalize();


  return 0;

}


