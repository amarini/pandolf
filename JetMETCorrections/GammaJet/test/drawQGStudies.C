#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"


bool useMCassoc_ = false;
bool ONEVTX = false;



int main(int argc, char* argv[]) {

  if( argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawQGStudies [LUMI_SHAPE] [photonID=\"medium\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset = "Photon_Run2011A-May10ReReco-v1";
  std::string mc_photonjet = "G_Summer11";
  //std::string mc_photonjet = "G_Pt-80to120_Tune23_7TeV_herwigpp_Summer11-PU_S3_START42_V11-v2";
  //std::string mc_QCD = mc_photonjet;
  std::string mc_QCD = "QCD_EMEnriched_Spring11";
  std::string recoType = "pf";
  std::string jetAlgo = "akt5";
  std::string norm( argv[1] );
  if( norm!="LUMI" && norm!="SHAPE" ) {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  std::string photonID="medium";
  if( argc==3 ) {
    std::string photonID_tmp(argv[2]);
    photonID = photonID_tmp;
  }

  std::string algoType;
  if( recoType=="calo" )
    algoType = jetAlgo;
  else
    algoType = recoType+jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5") algoType = "jptak5";
  


  bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  DrawBase* db = new DrawBase("QGStudies", recoType, jetAlgo);

  char dataFileName[150];
  if( photonID=="medium" )
    sprintf( dataFileName, "QGStudies_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  else
    sprintf( dataFileName, "QGStudies_%s_%s_%s.root", data_dataset.c_str(), algoType.c_str(), photonID.c_str());
  TFile* dataFile = TFile::Open(dataFileName);

  db->add_dataFile( dataFile, data_dataset );

  char mcPhotonJetFileName[150];
  if( photonID=="medium" )
    sprintf( mcPhotonJetFileName, "QGStudies_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
  else
    sprintf( mcPhotonJetFileName, "QGStudies_%s_%s_%s.root", mc_photonjet.c_str(), algoType.c_str(), photonID.c_str());
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName);

  db->add_mcFile( mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);


  char mcQCDFile_char[150];
  if( photonID=="medium" )
    sprintf( mcQCDFile_char, "QGStudies_%s_%s.root", mc_QCD.c_str(), algoType.c_str());
  else
    sprintf( mcQCDFile_char, "QGStudies_%s_%s_%s.root", mc_QCD.c_str(), algoType.c_str(), photonID.c_str());

  if( mc_QCD!=mc_photonjet ) {
    TFile* mcQCDFile = TFile::Open(mcQCDFile_char);
    db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);
  }



  if( norm=="LUMI" )
    db->set_lumiNormalization(191.);
  else
    db->set_shapeNormalization();

  db->set_lumi(191.);

  db->set_outputdir();


  bool log = true;

  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "GeV/c", "Events", log);
  db->drawHisto( "nvertex_PUW", "Number of Reconstructed Vertexes", "GeV/c", "Events", log);

  db->drawHisto( "ptJetReco", "Jet Transverse Momentum", "GeV/c", "Events", log);
  db->set_rebin(5);
  db->drawHisto( "ptPhot", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto( "ptPhot_3050", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto( "ptPhot_5080", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto( "ptPhot_80120", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->set_rebin(1);

  db->set_legendTitle("30 < p_{T} < 50 GeV");
  db->drawHisto( "nChargedJetReco_3050", "Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJetReco_3050", "Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "ptDJetReco_3050", "Jet p_{T} Distribution", "", "Events", log);
  db->set_rebin(5);
  db->drawHisto( "QGLikelihoodJetReco_3050", "Jet Q-G Likelihood", "", "Events", log);
  db->set_rebin(1);

  db->drawHisto( "nChargedJetReco_antibtag_3050", "Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJetReco_antibtag_3050", "Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "ptDJetReco_antibtag_3050", "Jet p_{T} Distribution", "", "Events", log);
  db->set_rebin(5);
  db->drawHisto( "QGLikelihoodJetReco_antibtag_3050", "Jet Q-G Likelihood", "", "Events", log);
  db->set_rebin(1);

  db->set_legendTitle("50 < p_{T} < 80 GeV");
  db->drawHisto( "nChargedJetReco_5080", "Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJetReco_5080", "Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "ptDJetReco_5080", "Jet p_{T} Distribution", "", "Events", log);
  db->set_rebin(2);
  db->drawHisto( "QGLikelihoodJetReco_5080", "Jet Q-G Likelihood", "", "Events", log);
  db->set_rebin(1);

  db->drawHisto( "nChargedJetReco_antibtag_5080", "Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJetReco_antibtag_5080", "Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "ptDJetReco_antibtag_5080", "Jet p_{T} Distribution", "", "Events", log);
  db->set_rebin(2);
  db->drawHisto( "QGLikelihoodJetReco_antibtag_5080", "Jet Q-G Likelihood", "", "Events", log);
  db->set_rebin(1);

  db->set_legendTitle("80 < p_{T} < 120 GeV");
  db->drawHisto( "nChargedJetReco_80120", "Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJetReco_80120", "Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "ptDJetReco_80120", "Jet p_{T} Distribution", "", "Events", log);
  db->set_rebin(2);
  db->drawHisto( "QGLikelihoodJetReco_80120", "Jet Q-G Likelihood", "", "Events", log);
  db->set_rebin(1);

  db->drawHisto( "nChargedJetReco_antibtag_80120", "Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJetReco_antibtag_80120", "Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "ptDJetReco_antibtag_80120", "Jet p_{T} Distribution", "", "Events", log);
  db->set_rebin(2);
  db->drawHisto( "QGLikelihoodJetReco_antibtag_80120", "Jet Q-G Likelihood", "", "Events", log);
  db->set_rebin(1);


  delete db;
  db = 0;

  return 0;

}  


