#include <stdlib.h>
#include "DrawBase.h"
#include "fitTools.h"





int main(int argc, char* argv[]) {

  if( argc != 7 && argc!=8 ) {
    std::cout << "USAGE: ./drawPhotonJet [(string)data_dataset] [(string)mc_BG_dataset] [(string)mc_SIGNAL_dataset] [(string)recoType] [(string)jetAlgo] [(string)norm ('LUMI' or 'SHAPE')] [(bool)useGenJets=false]" << std::endl;
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

  bool useGenJets = false;
  if( argc==8 ) {
    std::string useGenJets_str(argv[7]);
    if( useGenJets_str=="true") useGenJets=true;
  }

  bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  DrawBase* db = new DrawBase("PhotonJet", recoType, jetAlgo);
  db->set_pdf_aussi((bool)false);

  char outputdir_char[200];
  if( useGenJets ) {
    sprintf( outputdir_char, "PhotonJetPlots_%s_vs_%s_GENJETS_plus_%s_%s_%s", data_dataset.c_str(), mc_bg.c_str(), mc_photonjet.c_str(), algoType.c_str(), norm.c_str());
  } else {
    sprintf( outputdir_char, "PhotonJetPlots_%s_vs_%s_plus_%s_%s_%s", data_dataset.c_str(), mc_bg.c_str(), mc_photonjet.c_str(), algoType.c_str(), norm.c_str());
  }
  std::string outputdir_str(outputdir_char);

  std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();
  if( ptPhot_binning[0]==15. ) outputdir_str = outputdir_str + "_Pt15";

  db->set_outputdir(outputdir_str);

  char dataFileName[150];
  sprintf( dataFileName, "PhotonJet_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->set_dataFile( dataFile );

  char mcMinBiasFileName[150];
  if( useGenJets ) {
    sprintf( mcMinBiasFileName, "PhotonJet_%s_%s_GENJETS.root", mc_bg.c_str(), algoType.c_str());
  } else {
    sprintf( mcMinBiasFileName, "PhotonJet_%s_%s.root", mc_bg.c_str(), algoType.c_str());
  }
  TFile* mcMinBiasFile = TFile::Open(mcMinBiasFileName);
  std::cout << "Opened mc file '" << mcMinBiasFileName << "'." << std::endl;

  db->set_mcFile( mcMinBiasFile );

  char mcPhotonJetFileName[150];
//if( MCassoc ) {
//  sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s_MCassoc.root", mc_photonjet.c_str(), algoType.c_str());
//} else {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
//}
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

  db->drawHisto_2bkg( "ptPhot", MCbgName, MCsignalName, "", "clusterOK_isolated",1, log);
  db->drawHisto_onlyData( "ptPhot", "", "clusterOK_isolated",1, log);
  db->drawHisto_2bkg( "phiPhot", MCbgName, MCsignalName, "", "clusterOK_isolated",1);
  db->drawHisto_2bkg( "etaPhot", MCbgName, MCsignalName, "", "clusterOK_isolated",1);
  db->drawHisto_2bkg( "deltaPhi", MCbgName, MCsignalName, "", "Nm1",1);
  db->drawHisto_2bkg( "deltaPhi", MCbgName, MCsignalName, "", "clusterOK",1);

  db->drawHisto_2bkg( "ptSecondJetRel", MCbgName, MCsignalName, "", "Nm1",1, log);
  db->drawHisto_2bkg( "hcalIsoPhotReco", MCbgName, MCsignalName, "", "Nm1",1, log);
  db->drawHisto_2bkg( "ecalIsoPhotReco", MCbgName, MCsignalName, "", "Nm1",1, log);
  db->drawHisto_2bkg( "nTrkIsoPhotReco", MCbgName, MCsignalName, "", "Nm1",1, log);
  db->drawHisto_2bkg( "ptTrkIsoPhotReco", MCbgName, MCsignalName, "", "Nm1",1, log);
  db->drawHisto_2bkg( "clusterMajPhotReco", MCbgName, MCsignalName, "", "Nm1",1, log);
  db->drawHisto_2bkg( "clusterMinPhotReco", MCbgName, MCsignalName, "", "Nm1",1, log);

  db->drawHisto_2bkg( "clusterMajPhotReco", MCbgName, MCsignalName, "", "",1, log);
  db->drawHisto_2bkg( "clusterMinPhotReco", MCbgName, MCsignalName, "", "",1, log);

  db->drawHisto_2bkg( "ptSecondJetRel", MCbgName, MCsignalName, "", "clusterOK",1);
  db->drawHisto_2bkg( "hcalIsoPhotReco", MCbgName, MCsignalName, "", "clusterOK",1);
  db->drawHisto_2bkg( "ecalIsoPhotReco", MCbgName, MCsignalName, "", "clusterOK",1);
  db->drawHisto_2bkg( "nTrkIsoPhotReco", MCbgName, MCsignalName, "", "clusterOK",1);
  db->drawHisto_2bkg( "ptTrkIsoPhotReco", MCbgName, MCsignalName, "", "clusterOK",1);
  db->drawHisto_2bkg( "clusterMajPhotReco", MCbgName, MCsignalName, "", "isolated",1);
  db->drawHisto_2bkg( "clusterMinPhotReco", MCbgName, MCsignalName, "", "isolated",1);

  db->drawHisto_2bkg( "ptSecondJetRel", MCbgName, MCsignalName, "", "clusterOK_isolated",1);
  db->drawHisto_2bkg( "deltaPhi", MCbgName, MCsignalName, "", "clusterOK_isolated",1);
  db->drawHisto_onlyData( "deltaPhi", "", "clusterOK_isolated",1);

  db->drawHisto_2bkg( "response", MCbgName, MCsignalName, "", "loose",1);
  db->drawHisto_2bkg( "response", MCbgName, MCsignalName, "", "",1);
  db->drawHisto_2bkg( "response", MCbgName, MCsignalName, "", "clusterOK",1);
  if( recoType=="pf" )
    db->drawHisto_2bkg( "responseMPF", MCbgName, MCsignalName, "", "",1);

  delete db;
  db = 0;

  return 0;

}  


