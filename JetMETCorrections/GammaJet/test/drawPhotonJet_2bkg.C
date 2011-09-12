#include <stdlib.h>
#include "DrawBase.h"
#include "fitTools.h"


bool useMCassoc_ = false;
bool ONEVTX = false;



int main(int argc, char* argv[]) {

  if( argc != 7 && argc!=8 ) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')] [flags=\"\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset( argv[1] );
  std::string mc_photonjet( argv[2] );
  std::string mc_QCD( argv[3] );
  //std::string mc_Wenu( argv[4] );
  std::string recoType( argv[4] );
  std::string jetAlgo( argv[5] );
  std::string norm( argv[6] );
  if( norm!="LUMI" && norm!="SHAPE" ) {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }
  std::string flags = "";
  if( argc==8 ) {
    std::string flags_str(argv[7]);
    flags = flags_str;
  }


  std::string algoType;
  if( recoType=="calo" )
    algoType = jetAlgo;
  else
    algoType = recoType+jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5") algoType = "jptak5";
  

  Float_t etamax = 3.;

  bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  DrawBase* db = new DrawBase("PhotonJet", recoType, jetAlgo);
  db->set_pdf_aussi((bool)false);
//if( useMCassoc_ ) {
//  if( flags=="" ) { flags = "MCassoc"; }
//  else { flags = "MCassoc_" + flags; }
//}
  db->set_flags( flags );

  db->set_isCMSArticle( (bool)true );

std::cout << "flags set." << std::endl;
//char outputdir_char[200];
//std::string dirflags(flags);
//if( dirflags!="" ) {
//  sprintf( outputdir_char, "PhotonJetPlots_%s_vs_%s_plus_%s_%s_%s_%s", data_dataset.c_str(), mc_photonjet.c_str(), mc_QCD.c_str(),algoType.c_str(), norm.c_str(), dirflags.c_str());
//} else {
//  sprintf( outputdir_char, "PhotonJetPlots_%s_vs_%s_plus_%s_%s_%s", data_dataset.c_str(), mc_photonjet.c_str(), mc_QCD.c_str(), algoType.c_str(), norm.c_str());
//}
//std::string outputdir_str(outputdir_char);

//std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();
//if( ptPhot_binning[0]==15. ) outputdir_str = outputdir_str + "_Pt15";

  std::string ONEVTX_text = (ONEVTX) ? "_ONEVTX":"";

  char dataFileName[150];
  if( flags != "" ) {
    sprintf( dataFileName, "PhotonJet_%s_%s%s_%s.root", data_dataset.c_str(), algoType.c_str(), ONEVTX_text.c_str(), flags.c_str());
  } else {
    sprintf( dataFileName, "PhotonJet_%s_%s%s.root", data_dataset.c_str(), algoType.c_str(), ONEVTX_text.c_str());
  }
  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->add_dataFile( dataFile, data_dataset );

  char mcPhotonJetFileName[150];
  std::string signalflags(flags);
  if( useMCassoc_ ) signalflags = flags;
  if( signalflags!="" && flags!="GENJETS" ) {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s_%s.root", mc_photonjet.c_str(), algoType.c_str(), signalflags.c_str());
  } else {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName);
  std::cout << "Opened mc file '" << mcPhotonJetFileName << "'." << std::endl;

  db->add_mcFile( mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);


  char mcQCDFile_char[150];
  std::string bgflags(flags);
  if( useMCassoc_ && mc_QCD==mc_photonjet ) bgflags = signalflags;
  if( bgflags!="" ) {
    sprintf( mcQCDFile_char, "PhotonJet_%s_%s_%s.root", mc_QCD.c_str(), algoType.c_str(), bgflags.c_str());
  } else {
    sprintf( mcQCDFile_char, "PhotonJet_%s_%s.root", mc_QCD.c_str(), algoType.c_str());
  }
  TFile* mcQCDFile = TFile::Open(mcQCDFile_char);
  std::cout << "Opened mc file '" << mcQCDFile_char << "'." << std::endl;

  if( mc_QCD!=mc_photonjet )
    db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);


//char mcWenuFile_char[150];
//if( flags!="" ) {
//  sprintf( mcWenuFile_char, "PhotonJet_%s_%s_%s.root", mc_Wenu.c_str(), algoType.c_str(), flags.c_str());
//} else {
//  sprintf( mcWenuFile_char, "PhotonJet_%s_%s.root", mc_Wenu.c_str(), algoType.c_str());
//}
//TFile* mcWenuFile = TFile::Open(mcWenuFile_char);
//std::cout << "Opened mc file '" << mcWenuFile_char << "'." << std::endl;

//db->add_mcFile( mcWenuFile, "W #rightarrow e #nu MC", 8 );


  if( norm=="LUMI" ) 
    db->set_lumiNormalization(1000.);
  else
    db->set_shapeNormalization();

  db->set_lumi(191.);

  if( ONEVTX ) {
    flags += "_ONEVTX";
    db->set_flags(flags);
  }
  db->set_outputdir();

std::cout << "outputdir set." << std::endl;

  bool log = true;

  db->drawHisto( "ptPhot", "Photon Transverse Momentum", "GeV", "Events", log, 1, "passedID");
  db->drawHisto( "ptPhot", "Photon Transverse Momentum", "GeV", "Events", log);

  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "", "Events", log);
  db->drawHisto( "nvertexPU", "Number of Reconstructed Vertexes", "", "Events", log);
//db->drawHisto( "nvertex_passedID", "Number of Reconstructed Vertexes", "", "Events", log);
//db->drawHisto( "nvertex_passedID_pt50", "Number of Reconstructed Vertexes", "", "Events", log);

//db->drawHisto( "ptSecondJetRel_passedID", "Second Jet p_{T} / Photon p_{T}", "", "Events", log);
//db->drawHisto( "ptSecondJetRel_passedID_pt50", "Second Jet p_{T} / Photon p_{T}", "", "Events", log);
//db->drawHisto( "deltaPhi_2ndJet_passedID", "Photon - 2ndJet Azimuth Difference", "rad", "Events", log);
//db->drawHisto( "deltaPhi_2ndJet_passedID_pt50", "Photon - 2ndJet Azimuth Difference", "rad", "Events", log);

//db->drawHisto( "ptPhot", "", "loose", "", 1, log);
//db->drawHisto( "ptPhot", "", "medium", "", 1, log);
////db->drawHisto_onlyData( "ptPhot", "", "clusterOK_isolated", "", 1, log);
//db->drawHisto( "phiPhot", "", "clusterOK_isolated", "", 1);
//db->drawHisto( "etaPhot", "", "clusterOK_isolated", "", 1);
////db->drawHisto( "deltaPhi", "", "Nm1", "", 1);
//db->drawHisto( "deltaPhi", "", "clusterOK", "", 1);

//db->drawHisto( "deltaPhi_2ndJet", "", "loose", "", 1);

//db->drawHisto( "ptSecondJetRel", "", "Nm1", "", 1, log);
//db->drawHisto( "hcalIsoPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "hcalIsoEnergyPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "ecalIsoPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "ecalIsoEnergyPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "nTrkIsoPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "ptTrkIsoPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "clusterMajPhotReco", "", "Nm1", "", 1, log);
//db->drawHisto( "clusterMinPhotReco", "", "Nm1", "", 1, log);

//db->drawHisto( "clusterMajPhotReco", "", "", "", 1, log);
//db->drawHisto( "clusterMinPhotReco", "", "", "", 1, log);

//db->drawHisto( "ptSecondJetRel", "", "clusterOK", "", 1);
//db->drawHisto( "hcalIsoPhotReco", "", "clusterOK", "", 1);
//db->drawHisto( "hcalIsoEnergyPhotReco", "", "clusterOK", "", 1);
//db->drawHisto( "ecalIsoPhotReco", "", "clusterOK", "", 1);
//db->drawHisto( "ecalIsoEnergyPhotReco", "", "clusterOK", "", 1);
//db->drawHisto( "nTrkIsoPhotReco", "", "clusterOK", "", 1);
//db->drawHisto( "ptTrkIsoPhotReco", "", "clusterOK", "", 1);
//db->drawHisto( "clusterMajPhotReco", "", "isolated", "", 1);
//db->drawHisto( "clusterMinPhotReco", "", "isolated", "", 1);

//db->drawHisto( "ptSecondJetRel", "", "clusterOK_isolated", "", 1);
//db->drawHisto( "deltaPhi", "", "clusterOK_isolated", "", 1);
////db->drawHisto_onlyData( "deltaPhi", "", "clusterOK_isolated", "", 1);

//db->drawHisto( "response_loose", "Response", "", "Events", log);
  std::string metName;
  if( recoType=="pf" ) metName = "PF Missing E_{T}";
  if( recoType=="calo" ) metName = "PF Missing E_{T}";
  if( recoType=="jpt" ) metName = "PF Missing E_{T}";
   
  db->drawHisto( "met", metName, "GeV", "Events", log, 1, "");
  db->drawHisto( "met_noSmaj", (std::string)(metName+" (no cut on S_{maj})"), "GeV", "Events", log, 1, "");

  db->set_rebin(5);
  db->set_legendTitle("|#eta| < 1.1");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta011");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta011", (bool)true);
  db->drawHisto( "responseMPF", "MPF Response", "", "Events", log, 1, "eta011");
//db->set_legendTitle("|#eta| < 1.3");
//db->drawHisto( "response", "p_{T}/p_{T}^{#gamma} Response", "", "Events", log, 1, "eta013");
  db->set_legendTitle("|#eta^{#gamma}| < 1.3");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta013");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta013", (bool)true);
  db->drawHisto( "responseMPF", "MPF Response", "", "Events", log, 1, "eta013");
  db->set_legendTitle("1.5 < |#eta| < 2.4");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta1524");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta1524", (bool)true);
  db->drawHisto( "responseMPF", "MPF Response", "", "Events", log, 1, "eta1524");
  db->set_legendTitle("2.4 < |#eta| < 3");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta243");
  db->drawHisto( "response", "Balancing Response", "", "Events", log, 1, "eta243", (bool)true);
  db->drawHisto( "responseMPF", "MPF Response", "", "Events", log, 1, "eta243");
//db->set_legendTitle("3 < |#eta| < 5");
//db->drawHisto( "response", "p_{T}/p_{T}^{#gamma} Response", "", "Events", log, 1, "eta35");
//db->drawHisto( "responseMPF", "MPF Response", "", "Events", log, 1, "eta35");

  delete db;
  db = 0;

  return 0;

}  


