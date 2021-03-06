#include "Ntp1Finalizer_HWWlvjj.h"
#include "TMath.h"
#include <iostream>


float delta_phi(float phi1, float phi2) {

  float dphi = fabs(phi1 - phi2);
  float sgn = (phi1 >= phi2 ? +1. : -1.);
  return sgn * (dphi <= TMath::Pi() ? dphi : TMath::TwoPi() - dphi);
}





int main( int argc, char* argv[] ) {

  if( /*argc!=4 &&*/ argc!=5 ) {
    std::cout << "USAGE: ./finalize_HWWlvjj [dataset] [selectionType] [HiggsMass] [leptType=\"ALL\"]"  <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);

  float HiggsMass;
  if( argc==5 || argc==4 ) {
    HiggsMass = atoi(argv[3]);
  }

  std::string leptType="ALL";
  if( argc==4 ) {
    std::string leptType_str(argv[4]);
    leptType = leptType_str;
  }

  Ntp1Finalizer_HWWlvjj* nf = new Ntp1Finalizer_HWWlvjj( dataset, selectionType, HiggsMass );



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

  } else if( dataset=="ELEMUCombined" ) {

    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146734_147589");
    nf->addFile("Mu_Run2010A-Sep17ReReco_v2_runs135821_144114");
    nf->addFile("Mu_Run2010B_PromptReco_v2_runs146240_146733");
    nf->addFile("EG_upto146724");

  } else if( dataset=="EleMu_Nov4ReReco_PU" ) {

    nf->addFile("Electron_Nov4ReReco_PU");
    nf->addFile("Mu_Nov4ReReco_PU");

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


  } else if( dataset=="VVtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2" ) {

    //nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "WWtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1" );

  } else if( dataset=="VVtoAnything_TuneZ2_7TeV-pythia6-tauola" ) {

    nf->addFile( "WWtoAnything_TuneZ2_7TeV-pythia6-tauola" );
    nf->addFile( "WZtoAnything_TuneZ2_7TeV-pythia6-tauola" );
    nf->addFile( "ZZtoAnything_TuneZ2_7TeV-pythia6-tauola" );

  } else if( dataset=="DY_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2" ) {

    nf->addFile( "DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "DYToEE_M-20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2" );

  } else if( dataset=="DY_TuneZ2_7TeV-pythia6" ) {

    nf->addFile( "DYToEE_M-10To20_TuneZ2_7TeV-pythia6" );
    nf->addFile( "DYToEE_M-20_TuneZ2_7TeV-pythia6" );
    nf->addFile( "DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6" );
    nf->addFile( "DYToMuMu_M-20_TuneZ2_7TeV-pythia6" );

  } else if( dataset=="TToBLNu_TuneZ2_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2" ) {

    nf->addFile( "TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2" );
    nf->addFile( "TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2" );

  } else if( dataset=="TToBLNu_TuneZ2_7TeV-madgraph" ) {

     nf->addFile( "TToBLNu_TuneZ2_s-channel_7TeV-madgraph" );
     nf->addFile( "TToBLNu_TuneZ2_t-channel_7TeV-madgraph" ); 
     nf->addFile( "TToBLNu_TuneZ2_tW-channel_7TeV-madgraph" );

  } else if( dataset=="QCD_EMEnriched_TuneZ2_7TeV-pythia6_3" ) {

    nf->addFile( "QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6_3" );
    nf->addFile( "QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia6_3" );
    nf->addFile( "QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6_3" );

  } else if( dataset=="QCD_BCtoE_TuneZ2_7TeV-pythia6" ) {

    nf->addFile( "QCD_Pt-20to30_BCtoE_TuneZ2_7TeV-pythia6" );
    nf->addFile( "QCD_Pt-30to80_BCtoE_TuneZ2_7TeV-pythia6" );
    nf->addFile( "QCD_Pt-80to170_BCtoE_TuneZ2_7TeV-pythia6" );

  } else if( dataset=="GJet_TuneZ2_7TeV-alpgen" ) {

    nf->addFile( "G1Jet_Pt-120to180_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G1Jet_Pt-180to240_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G1Jet_Pt-20to60_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G1Jet_Pt-300to5000_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G1Jet_Pt-60to120_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G2Jets_Pt-120to180_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G2Jets_Pt-180to240_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G2Jets_Pt-20to60_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G2Jets_Pt-240to300_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G2Jets_Pt-300to5000_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G2Jets_Pt-60to120_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G3Jets_Pt-120to180_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G3Jets_Pt-180to240_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G3Jets_Pt-20to60_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G3Jets_Pt-240to300_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G3Jets_Pt-300to5000_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G3Jets_Pt-60to120_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G4Jets_Pt-120to180_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G4Jets_Pt-180to240_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G4Jets_Pt-20to60_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G4Jets_Pt-240to300_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G4Jets_Pt-300to5000_TuneZ2_7TeV-alpgen" );
    nf->addFile( "G4Jets_Pt-60to120_TuneZ2_7TeV-alpgen" );

  } else if ( dataset == "DATA_6july" ) {

     nf->addFile("SingleMu_6july");
     nf->addFile("SingleElectron_6july");

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


