#include <stdlib.h>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {

  if( argc != 7 && argc!=8 ) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_BG_dataset] [mc_SIGNAL_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')] [flags=""]" << std::endl;
    exit(23);
  }

  std::string data_dataset( argv[1] );
  std::string mc_bg( argv[2] );
  std::string mc_photonjet( argv[3] );
  std::string recoType( argv[4] );
  std::string jetAlgo( argv[5] );
  std::string norm( argv[6] );
  if( norm!="LUMI" && norm!="SHAPE" ) {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  std::string algoType;
  if( recoType=="calo" )
    algoType = jetAlgo;
  else
    algoType = recoType+jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5") algoType = "jptak5";
  

  Float_t etamax = 3.;

  std::string flags = "";
  if( argc==8 ) {
    std::string flags_str(argv[7]);
    flags = flags_str;
  }

  bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  DrawBase* db = new DrawBase("PhotonJet", recoType, jetAlgo);
  db->set_pdf_aussi((bool)false);

  char outputdir_char[200];
  if( flags!="" ) {
    sprintf( outputdir_char, "PhotonJetPlots_%s_vs_%s_plus_%s_%s_%s_%s", data_dataset.c_str(), mc_bg.c_str(), mc_photonjet.c_str(), algoType.c_str(), norm.c_str(), flags.c_str());
  } else {
    sprintf( outputdir_char, "PhotonJetPlots_%s_vs_%s_plus_%s_%s_%s", data_dataset.c_str(), mc_bg.c_str(), mc_photonjet.c_str(), algoType.c_str(), norm.c_str());
  }
  std::string outputdir_str(outputdir_char);

  std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();
  if( ptPhot_binning[0]==15. ) outputdir_str = outputdir_str + "_Pt15";

  db->set_outputdir(outputdir_str);

  char dataFileName[150];
  if( flags != "" ) {
    sprintf( dataFileName, "PhotonJet_%s_%s_%s.root", data_dataset.c_str(), algoType.c_str(), flags.c_str());
  } else {
    sprintf( dataFileName, "PhotonJet_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  }
  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->set_dataFile( dataFile );

  char mcBGFile[150];
  if( flags!="" ) {
    sprintf( mcBGFile, "PhotonJet_%s_%s_%s.root", mc_bg.c_str(), algoType.c_str(), flags.c_str());
  } else {
    sprintf( mcBGFile, "PhotonJet_%s_%s.root", mc_bg.c_str(), algoType.c_str());
  }
  TFile* mcMinBiasFile = TFile::Open(mcBGFile);
  std::cout << "Opened mc file '" << mcBGFile << "'." << std::endl;

  db->set_mcFile( mcMinBiasFile );

  char mcPhotonJetFileName[150];
  if( flags!="" && flags!="GENJETS" ) {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s_%s.root", mc_photonjet.c_str(), algoType.c_str(), flags.c_str());
  } else {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName);
  std::cout << "Opened mc file '" << mcPhotonJetFileName << "'." << std::endl;

  db->set_mcFile2( mcPhotonJetFile );

  if( norm=="LUMI" )
    db->set_crossSectionNormalization();
  else
    db->set_shapeNormalization();


  bool log = true;

  std::string MCbgName = "QCD MC";
  std::string MCsignalName = "#gamma+jet MC";

  db->set_mcName( MCbgName );
  db->set_mcName2( MCsignalName );

  db->drawHisto_2bkg( "ptPhot", "", "clusterOK_isolated",1, log);
  db->drawHisto_2bkg( "ptPhot", "", "loose",1, log);
  db->drawHisto_2bkg( "ptPhot", "", "medium",1, log);
  db->drawHisto_onlyData( "ptPhot", "", "clusterOK_isolated",1, log);
  db->drawHisto_2bkg( "phiPhot", "", "clusterOK_isolated",1);
  db->drawHisto_2bkg( "etaPhot", "", "clusterOK_isolated",1);
  db->drawHisto_2bkg( "deltaPhi", "", "Nm1",1);
  db->drawHisto_2bkg( "deltaPhi", "", "clusterOK",1);

  db->drawHisto_2bkg( "ptSecondJetRel", "", "Nm1",1, log);
  db->drawHisto_2bkg( "hcalIsoPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "hcalIsoEnergyPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "ecalIsoPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "ecalIsoEnergyPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "nTrkIsoPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "ptTrkIsoPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "clusterMajPhotReco", "", "Nm1",1, log);
  db->drawHisto_2bkg( "clusterMinPhotReco", "", "Nm1",1, log);

  db->drawHisto_2bkg( "clusterMajPhotReco", "", "",1, log);
  db->drawHisto_2bkg( "clusterMinPhotReco", "", "",1, log);

  db->drawHisto_2bkg( "ptSecondJetRel", "", "clusterOK",1);
  db->drawHisto_2bkg( "hcalIsoPhotReco", "", "clusterOK",1);
  db->drawHisto_2bkg( "hcalIsoEnergyPhotReco", "", "clusterOK",1);
  db->drawHisto_2bkg( "ecalIsoPhotReco", "", "clusterOK",1);
  db->drawHisto_2bkg( "ecalIsoEnergyPhotReco", "", "clusterOK",1);
  db->drawHisto_2bkg( "nTrkIsoPhotReco", "", "clusterOK",1);
  db->drawHisto_2bkg( "ptTrkIsoPhotReco", "", "clusterOK",1);
  db->drawHisto_2bkg( "clusterMajPhotReco", "", "isolated",1);
  db->drawHisto_2bkg( "clusterMinPhotReco", "", "isolated",1);

  db->drawHisto_2bkg( "ptSecondJetRel", "", "clusterOK_isolated",1);
  db->drawHisto_2bkg( "deltaPhi", "", "clusterOK_isolated",1);
  db->drawHisto_onlyData( "deltaPhi", "", "clusterOK_isolated",1);

  db->drawHisto_2bkg( "response", "", "loose",1);
  db->drawHisto_2bkg( "response", "", "",1);
  db->drawHisto_2bkg( "response", "", "clusterOK",1);
  if( recoType=="pf" ) {
    db->drawHisto_2bkg( "responseMPF", "", "",1);
    db->drawHisto_2bkg( "responseMPF", "", "loose",1);
  }

  delete db;
  db = 0;

  return 0;

}  


