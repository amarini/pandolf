#include <stdlib.h>
#include "DrawBase.h"
#include "fitTools.h"




int main(int argc, char* argv[]) {

  if( argc != 5 ) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [norm ('LUMI' or 'SHAPE')]" << std::endl;
    exit(23);
  }

  std::string data_dataset( argv[1] );
  std::string mc_photonjet( argv[2] );
  std::string mc_QCD( argv[3] );
  std::string norm( argv[4] );
  if( norm!="LUMI" && norm!="SHAPE" ) {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }



  DrawBase* db = new DrawBase("NP", "pf", "akt5");
  db->set_pdf_aussi((bool)false);
  //db->set_isCMSArticle( (bool)true );



  char dataFileName[150];
  sprintf( dataFileName, "NP_%s.root", data_dataset.c_str());
  TFile* dataFile = TFile::Open(dataFileName);

  db->add_dataFile( dataFile, data_dataset );

  char mcPhotonJetFileName[150];
  sprintf( mcPhotonJetFileName, "NP_%s.root", mc_photonjet.c_str());
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName);

  db->add_mcFile( mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);


  char mcQCDFile_char[150];
  sprintf( mcQCDFile_char, "NP_%s.root", mc_QCD.c_str());
  TFile* mcQCDFile = TFile::Open(mcQCDFile_char);

  if( mc_QCD!=mc_photonjet )
    db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);



  if( norm=="LUMI" )
    db->set_lumiNormalization(36.);
  else
    db->set_shapeNormalization();

  db->set_lumi(36.);

  db->set_outputdir();


  bool log = true;


  db->drawHisto("etaPhot", "Photon Pseudorapidity", "", "Events");
  db->drawHisto("phiPhot", "Photon Azimuth", "", "Events");

  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 20 GeV/c)", "", "Events", log);

  db->set_rebin(2);
  db->drawHisto("ptPhot", "Photon Transverse Momentum", "GeV/c", "Events", log );
  db->set_legendTitle("nJets = 0");
  db->drawHisto("ptPhot_0jet", "Photon Transverse Momentum", "GeV/c", "Events", log);

  db->set_legendTitle("nJets = 1");
  db->drawHisto("ptPhot_1jet", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_1jet", "Jet Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(5);
  db->drawHisto("mGammaJet_1jet", "Gamma-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtGammaJet_1jet", "Gamma-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(2);
  db->drawHisto("met_1jet", "Missing Transverse Energy", "GeV", "Events", log);

  db->set_legendTitle("nJets = 2");
  db->drawHisto("ptPhot_2jet", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_2jet", "Jet Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(5);
  db->drawHisto("mGammaJet_2jet", "Gamma-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtGammaJet_2jet", "Gamma-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(2);
  db->drawHisto("met_2jet", "Missing Transverse Energy", "GeV", "Events", log);

  db->set_legendTitle("nJets = 3");
  db->drawHisto("ptPhot_3jet", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_3jet", "Jet Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(5);
  db->drawHisto("mGammaJet_3jet", "Gamma-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtGammaJet_3jet", "Gamma-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(2);
  db->drawHisto("met_3jet", "Missing Transverse Energy", "GeV", "Events", log);

  db->set_legendTitle("nJets = 4");
  db->drawHisto("ptPhot_4jet", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto("mJet_4jet", "Jet Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(5);
  db->drawHisto("mGammaJet_4jet", "Gamma-Jet Invariant Mass", "GeV/c^{2}", "Events", log);
  db->drawHisto("mtGammaJet_4jet", "Gamma-Jet Transverse Mass", "GeV/c^{2}", "Events", log);
  db->set_rebin(2);
  db->drawHisto("met_4jet", "Missing Transverse Energy", "GeV", "Events", log);


  delete db;
  db = 0;

  return 0;

}  


