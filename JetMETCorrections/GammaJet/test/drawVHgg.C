#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if( argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawVHgg [LUMI_SHAPE] [photonID=\"medium\"]" << std::endl;
    exit(23);
  }

  



  DrawBase* db = new DrawBase("VHgg");


  std::string mcSignalFileName = "VHgg_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root";
  TFile* mcSignalFile = TFile::Open(mcSignalFileName.c_str());
  db->add_mcFile( mcSignalFile, "signal", "VH (120)", 46);

  std::string mcBGFileName = "VHgg_DiPhotonBoxBorn";
  TFile* mcBGFile = TFile::Open(mcBGFileName.c_str());
  db->add_mcFile( mcBGFile, "BG", "SM DiPhoton", 38);

  db->set_shapeNormalization();

  db->set_outputdir();


  bool log = true;

  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "", "Events", log);
  db->drawHisto( "nvertex_PU", "Number of Reconstructed Vertexes", "", "Events", log);

  db->drawHisto( "ptPhot1", "Lead Photon p_{T}", "GeV", "Events", log);
  db->drawHisto( "etaPhot1", "Lead Photon #eta", "", "Events", log);

  db->drawHisto( "ptPhot2", "Sublead Photon p_{T}", "GeV", "Events", log);
  db->drawHisto( "etaPhot2", "Sublead Photon #eta", "", "Events", log);

  db->drawHisto( "mgg", "DiPhoton Invariant Mass", "GeV", "Events", log);



  db->drawHisto( "ptJet1", "Lead Jet p_{T}", "GeV", "Events", log);
  db->drawHisto( "etaJet1", "Lead Jet #eta", "", "Events", log);
  db->drawHisto( "QGLikelihoodJet1", "Lead Jet QG Likelihood", "", "Events", log);

  db->drawHisto( "ptJet2", "Sublead Jet p_{T}", "GeV", "Events", log);
  db->drawHisto( "etaJet2", "Sublead Jet #eta", "", "Events", log);
  db->drawHisto( "QGLikelihoodJet2", "Sublead Jet QG Likelihood", "", "Events", log);

  db->drawHisto( "mjj", "DiJet Invariant Mass", "GeV", "Events", log);
  db->drawHisto( "QGLikelihoodProd", "QG Likelihood Product", "", "Events", log);



  delete db;
  db = 0;

  return 0;

}  




