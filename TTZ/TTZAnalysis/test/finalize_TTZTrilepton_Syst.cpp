
#include "Ntp1Finalizer_TTZTrilepton.h"
#include "TMath.h"
#include <iostream>



void runOnAllDatasets( Ntp1Finalizer_TTZTrilepton* nf );
void runOnSingleDataset( Ntp1Finalizer_TTZTrilepton* nf, const std::string& dataset );
void haddBGFiles( const std::string& selection, const std::string& flag );
std::string getFileName( const std::string& dataset, const std::string& selection, const std::string& suffix );




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

  // run once with no syst
  runOnAllDatasets( nf );
  haddBGFiles( selectionType, "" );


  std::cout << std::endl << "++++++++++++++++++++++" << std::endl;
  std::cout << "+++ BTag Syst +1 Sigma" << std::endl;
  std::cout << "++++++++++++++++++++++" << std::endl << std::endl;
  nf->set_btagSyst(1);

  runOnAllDatasets( nf );
  haddBGFiles( selectionType, "BTagUP" );


  std::cout << std::endl << "++++++++++++++++++++++" << std::endl;
  std::cout << "+++ BTag Syst -1 Sigma" << std::endl;
  std::cout << "++++++++++++++++++++++" << std::endl << std::endl;
  nf->set_btagSyst(-1);

  runOnAllDatasets( nf );
  haddBGFiles( selectionType, "BTagDOWN" );

  nf->set_btagSyst(0);


  std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
  std::cout << "+++ JES Syst +1 Sigma" << std::endl;
  std::cout << "+++++++++++++++++++++" << std::endl << std::endl;
  nf->set_jes(1);

  runOnAllDatasets( nf );
  haddBGFiles( selectionType, "JESUP" );

  std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
  std::cout << "+++ JES Syst -1 Sigma" << std::endl;
  std::cout << "+++++++++++++++++++++" << std::endl << std::endl;
  nf->set_jes(-1);

  runOnAllDatasets( nf );
  haddBGFiles( selectionType, "JESDOWN" );

  nf->set_jes(0);

  std::cout << std::endl << "+++++++++++++++++++++" << std::endl;
  std::cout << "+++ JER Syst +1 Sigma" << std::endl;
  std::cout << "+++++++++++++++++++++" << std::endl << std::endl;
  nf->set_jer(true);

  runOnAllDatasets( nf );
  haddBGFiles( selectionType, "JERUP" );

  return 0;

}




void runOnAllDatasets( Ntp1Finalizer_TTZTrilepton* nf ) {

  
  runOnSingleDataset( nf, "TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi" );
  runOnSingleDataset( nf, "VV_Summer11" );
  runOnSingleDataset( nf, "TTJ_Fall11_highstat" );
  runOnSingleDataset( nf, "DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11" );
  //runOnSingleDataset( nf, "BG" );

}


void runOnSingleDataset( Ntp1Finalizer_TTZTrilepton* nf, const std::string& dataset ) {

  nf->clearTree();

  Ntp1Finalizer_TTZTrilepton* nf_syst = new Ntp1Finalizer_TTZTrilepton( *nf );
  nf_syst->set_dataset(dataset);

  if( dataset=="VV_Summer11" ) {
  
    nf_syst->addFile("WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1");
    nf_syst->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf_syst->addFile("WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
  
  } else if( dataset=="BG" ) {
  
    nf_syst->addFile("WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1");
    nf_syst->addFile("ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf_syst->addFile("WW_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1");
    nf_syst->addFile("TTJ_Fall11_highstat");
    nf_syst->addFile("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11");
  
  } else {

    nf_syst->addFile( dataset );
 
  }

  nf_syst->finalize();

} 




void haddBGFiles( const std::string& selection, const std::string& flag ) {

  std::string suffix = flag;
  if( suffix!="" ) suffix = "_" + suffix;

  std::string hadd_command = "hadd -f ";

  hadd_command += getFileName( "BG", selection, suffix );
  hadd_command += " ";
  hadd_command += getFileName( "VV_Summer11", selection, suffix );
  hadd_command += " ";
  hadd_command += getFileName( "TTJ_Fall11_highstat", selection, suffix );
  hadd_command += " ";
  hadd_command += getFileName( "DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11", selection, suffix );
  hadd_command += " ";

  std::cout << "-> Hadding BG files:" << std::endl;
  std::cout << hadd_command << std::endl;
  system(hadd_command.c_str());

}


std::string getFileName( const std::string& dataset, const std::string& selection, const std::string& suffix ) {

  std::string fileName = "TTZTrilepton_" + dataset + "_" + selection + "_TCHE_ALL" + suffix + ".root";

  return fileName;

}
