#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



bool dijet_selection=true;


void drawHistoWithQuarkGluonComponents( DrawBase* db, const std::string& treeName, const std::string& varName, const std::string& axisName, const std::string& units="", const std::string& instanceName="Events", bool log=false, float ptMin=0., float ptMax=10000. );
std::pair<TH1D,TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db );
//void drawSignalPtMix( std::pair<TH1D*,TH1D*> h1pair_3050, std::pair<TH1D*,TH1D*> h1pair_5080, std::pair<TH1D*,TH1D*> h1pair_80120, int mass, float frac_3050, float frac_5080, float frac_80120, DrawBase* db );


int main(int argc, char* argv[]) {

  if( argc!=1 && argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawMultiJetQG [data_dataset=\"HT_Run2011B-PromptReco-v1_HLT\"] [normalization=\"SHAPE\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset = "HT_Run2011_FULL";
  if( argc>1 ) {
    std::string dataset_tmp(argv[1]);
    data_dataset = dataset_tmp;
  }

  //std::string data_dataset = "DATA_Run2011A_1fb";
  //std::string mc_QCD = "QCD_TuneZ2_HT-500To1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1";
  std::string mc_QCD = "QCD_HT_Summer11";
  std::string mc_PhotonJet = "G_Summer11";


  std::string recoType = "pf";
  std::string jetAlgo = "akt5";
  std::string norm = "SHAPE";
  if( argc==3 ) {
    std::string norm_tmp(argv[2]);
    norm = norm_tmp;
  }
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
  



  DrawBase* db = new DrawBase("MultiJetQG", recoType, jetAlgo);

  char dataFileName[150];
  if( dijet_selection )
    sprintf( dataFileName, "MultiJet_%s_DIJET.root", data_dataset.c_str());
  else
    sprintf( dataFileName, "MultiJet_%s.root", data_dataset.c_str());
  TFile* dataFile = TFile::Open(dataFileName);

  db->add_dataFile( dataFile, data_dataset );



  char mcQCDFile_char[150];
  if( dijet_selection )
    sprintf( mcQCDFile_char, "MultiJet_%s_DIJET.root", mc_QCD.c_str());
  else
    sprintf( mcQCDFile_char, "MultiJet_%s.root", mc_QCD.c_str());

  TFile* mcQCDFile = TFile::Open(mcQCDFile_char);
  db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);


  //char mcPhotonJetFile_char[150];
  //sprintf( mcPhotonJetFile_char, "MultiJet_%s.root", mc_PhotonJet.c_str());

  //TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFile_char);
  //db->add_mcFile( mcPhotonJetFile, mc_PhotonJet, "#gamma+Jet MC", 46);



  if( norm=="LUMI" ) {
    db->set_lumiNormalization(4600.);
    std::cout << "-> Lumi normalization." << std::endl;
  } else {
    std::cout << "-> Shape normalization." << std::endl;
    db->set_shapeNormalization();
    db->set_lumi(4600.);
  }





  db->set_outputdir();


  db->set_lumiOnRightSide();


/*
  bool log = true;

  db->set_yAxisMaxScaleLog( 1000. );
  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "", "Events", log);
  db->drawHisto( "nvertexPU", "Number of Reconstructed Vertexes", "", "Events", log);

  db->drawHisto( "rhoPF", "Particle Flow Energy Density", "GeV", "Events", log);
  db->drawHisto( "rhoPFPU", "Particle Flow Energy Density", "GeV", "Events", log);

  db->set_xAxisRange( 650., 1300. );
  db->drawHisto( "ht_akt5", "CaloJet H_{T}", "GeV", "Events", log);
  db->set_xAxisRange();
  db->drawHisto( "htmet_akt5", "CaloJet H_{T} + CaloME_{T}", "GeV", "Events", log);
  db->drawHisto( "sumpt_pfakt5", "Sump_{T}^{4}", "GeV", "Events", log);

  db->set_rebin(5);
  db->drawHisto( "ptJet0", "First Jet p_{T}", "GeV", "Events", log);
  db->drawHisto( "ptJet1", "Second Jet p_{T}", "GeV", "Events", log);
  db->drawHisto( "ptJet2", "Third Jet p_{T}", "GeV", "Events", log);
  db->drawHisto( "ptJet3", "Fourth Jet p_{T}", "GeV", "Events", log);

  db->set_rebin();
  db->drawHisto( "ptDJet0", "First Jet p_{T}D",  "", "Events", log);
  db->drawHisto( "ptDJet1", "Second Jet p_{T}D", "", "Events", log);
  db->drawHisto( "ptDJet2", "Third Jet p_{T}D",  "", "Events", log);
  db->drawHisto( "ptDJet3", "Fourth Jet p_{T}D", "", "Events", log);

  db->set_xAxisMax(60.);
  db->drawHisto( "nChargedJet0", "First Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nChargedJet1", "Second Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nChargedJet2", "Third Jet Charged Multiplicity", "", "Events", log);
  db->drawHisto( "nChargedJet3", "Fourth Jet Charged Multiplicity", "", "Events", log);

  db->drawHisto( "nNeutralJet0", "First Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJet1", "Second Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJet2", "Third Jet Neutral Multiplicity", "", "Events", log);
  db->drawHisto( "nNeutralJet3", "Fourth Jet Neutral Multiplicity", "", "Events", log);
  db->set_xAxisMax();

*/

  std::string selection_nopt;
  std::string selection_noptW;
  std::string selection_pt0;
  std::string selection_pt1;



/*
  // HLT_HT150
  selection_nopt = "ptJet2<0.2*(ptJet0+ptJet1)/2. && (passed_HT150||run<5) && ht_akt5>160.";
  selection_noptW = "eventWeight*(" + selection_nopt + ")";
  selection_pt0 = "eventWeight*(" + selection_nopt + " && ptJet0 > 50. && ptJet0 < 100.)";
  selection_pt1 = "eventWeight*(" + selection_nopt + " && ptJet1 > 50. && ptJet1 < 100.)";

  db->drawHisto_fromTree( "tree_passedEvents", "nvertex", selection_noptW.c_str(), 30, 0.5, 30.5, "nvertex_HT150", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "rhoPF", selection_noptW.c_str(), 50, 0., 20., "rhoPF_HT150", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "ht_akt5", selection_noptW.c_str(), 50, 150., 500., "ht_akt5_HT150", "Calo H_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet0", selection_noptW.c_str(), 50, 0., 250., "ptJet0_HT150", "Lead Jet p_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet1", selection_noptW.c_str(), 50, 0., 250., "ptJet1_HT150", "Sublead Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "50 < p_{T} < 100 GeV");
  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_HT150", "Lead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_HT150", "Lead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_HT150", "Lead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_HT150", "Lead Jet Q-G LD", "", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_HT150", "Sublead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_HT150", "Sublead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_HT150", "Sublead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_HT150", "Sublead Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");

  // HLT_HT250
  selection_nopt = "ptJet2<0.2*(ptJet0+ptJet1)/2. && (passed_HT250||run<5) && ht_akt5>265.";
  selection_noptW = "eventWeight*(" + selection_nopt + ")";
  selection_pt0 = "eventWeight*(" + selection_nopt + " && ptJet0 > 100. && ptJet0 < 150.)";
  selection_pt1 = "eventWeight*(" + selection_nopt + " && ptJet1 > 100. && ptJet1 < 150.)";

  db->drawHisto_fromTree( "tree_passedEvents", "nvertex", selection_noptW.c_str(), 30, 0.5, 30.5, "nvertex_HT250", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "rhoPF", selection_noptW.c_str(), 50, 0., 20., "rhoPF_HT250", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "ht_akt5", selection_noptW.c_str(), 50, 150., 550., "ht_akt5_HT250", "Calo H_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet0", selection_noptW.c_str(), 50, 0., 250., "ptJet0_HT250", "Lead Jet p_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet1", selection_noptW.c_str(), 50, 0., 250., "ptJet1_HT250", "Sublead Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "100 < p_{T} < 150 GeV");
  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_HT250", "Lead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_HT250", "Lead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_HT250", "Lead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_HT250", "Lead Jet Q-G LD", "", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_HT250", "Sublead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_HT250", "Sublead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_HT250", "Sublead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_HT250", "Sublead Jet Q-G LD", "", "Events", true);


  // HLT_HT350
  selection_nopt = "ptJet2<0.2*(ptJet0+ptJet1)/2. && (passed_HT350||run<5) && ht_akt5>365.";
  selection_noptW = "eventWeight*(" + selection_nopt + ")";
  selection_pt0 = "eventWeight*(" + selection_nopt + " && ptJet0 > 150. && ptJet0 < 200.)";
  selection_pt1 = "eventWeight*(" + selection_nopt + " && ptJet1 > 150. && ptJet1 < 200.)";

  db->drawHisto_fromTree( "tree_passedEvents", "nvertex", selection_noptW.c_str(), 30, 0.5, 30.5, "nvertex_HT350", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "rhoPF", selection_noptW.c_str(), 50, 0., 20., "rhoPF_HT350", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "ht_akt5", selection_noptW.c_str(), 50, 150., 600., "ht_akt5_HT350", "Calo H_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet0", selection_noptW.c_str(), 50, 0., 250., "ptJet0_HT350", "Lead Jet p_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet1", selection_noptW.c_str(), 50, 0., 250., "ptJet1_HT350", "Sublead Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "150 < p_{T} < 200 GeV");
  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_HT350", "Lead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_HT350", "Lead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_HT350", "Lead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_HT350", "Lead Jet Q-G LD", "", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_HT350", "Sublead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_HT350", "Sublead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_HT350", "Sublead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_HT350", "Sublead Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");

  // HLT_HT400
  selection_nopt = "ptJet2<0.2*(ptJet0+ptJet1)/2. && (passed_HT400||run<5) && ht_akt5>420.";
  selection_noptW = "eventWeight*(" + selection_nopt + ")";
  selection_pt0 = "eventWeight*(" + selection_nopt + " && ptJet0 > 200. && ptJet0 < 250.)";
  selection_pt1 = "eventWeight*(" + selection_nopt + " && ptJet1 > 200. && ptJet1 < 250.)";

  db->drawHisto_fromTree( "tree_passedEvents", "nvertex", selection_noptW.c_str(), 30, 0.5, 30.5, "nvertex_HT400", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "rhoPF", selection_noptW.c_str(), 50, 0., 20., "rhoPF_HT400", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "ht_akt5", selection_noptW.c_str(), 50, 150., 700., "ht_akt5_HT400", "Calo H_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet0", selection_noptW.c_str(), 50, 0., 250., "ptJet0_HT400", "Lead Jet p_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet1", selection_noptW.c_str(), 50, 0., 250., "ptJet1_HT400", "Sublead Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "200 < p_{T} < 300 GeV");
  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_HT400", "Lead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_HT400", "Lead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_HT400", "Lead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_HT400", "Lead Jet Q-G LD", "", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_HT400", "Sublead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_HT400", "Sublead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_HT400", "Sublead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_HT400", "Sublead Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");

*/
  // HLT_HT500
  selection_nopt = "ptJet2<0.2*(ptJet0+ptJet1)/2. && (passed_HT500||run<5) && ht_akt5>525.";
  selection_noptW = "eventWeight*(" + selection_nopt + ")";
  selection_pt0 = "eventWeight*(" + selection_nopt + " && ptJet0 > 250. && ptJet0 < 300. && QGLikelihoodJet0>0. && QGLikelihoodJet0<1.)";
  selection_pt1 = "eventWeight*(" + selection_nopt + " && ptJet1 > 250. && ptJet1 < 300. && QGLikelihoodJet1>0. && QGLikelihoodJet1<1.)";

  db->drawHisto_fromTree( "tree_passedEvents", "nvertex", selection_noptW.c_str(), 30, 0.5, 30.5, "nvertex_HT500", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "rhoPF", selection_noptW.c_str(), 50, 0., 20., "rhoPF_HT500", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "ht_akt5", selection_noptW.c_str(), 50, 200., 1000., "ht_akt5_HT500", "Calo H_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet0", selection_noptW.c_str(), 50, 0., 250., "ptJet0_HT500", "Lead Jet p_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet1", selection_noptW.c_str(), 50, 0., 250., "ptJet1_HT500", "Sublead Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "250 < p_{T} < 300 GeV");
  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_HT500", "Lead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_HT500", "Lead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_HT500", "Lead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_HT500", "Lead Jet Q-G LD", "", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_HT500", "Sublead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_HT500", "Sublead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_HT500", "Sublead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_HT500", "Sublead Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");

  // HLT_HT600
  selection_nopt = "ptJet2<0.2*(ptJet0+ptJet1)/2. && (passed_HT600||run<5) && ht_akt5>640.";
  selection_noptW = "eventWeight*(" + selection_nopt + ")";
  selection_pt0 = "eventWeight*(" + selection_nopt + " && ptJet0 > 300. && ptJet0 < 400.&& QGLikelihoodJet0>0. && QGLikelihoodJet0<1.)";
  selection_pt1 = "eventWeight*(" + selection_nopt + " && ptJet1 > 300. && ptJet1 < 400.&& QGLikelihoodJet1>0. && QGLikelihoodJet1<1.)";

  db->drawHisto_fromTree( "tree_passedEvents", "nvertex", selection_noptW.c_str(), 30, 0.5, 30.5, "nvertex_HT600", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "rhoPF", selection_noptW.c_str(), 50, 0., 20., "rhoPF_HT600", "Particle Flow Energy Density (#rho)", "GeV", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "ht_akt5", selection_noptW.c_str(), 50, 300., 1200., "ht_akt5_HT600", "Calo H_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet0", selection_noptW.c_str(), 50, 0., 250., "ptJet0_HT600", "Lead Jet p_{T}", "GeV", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptJet1", selection_noptW.c_str(), 50, 0., 250., "ptJet1_HT600", "Sublead Jet p_{T}", "GeV", "Events", true);

  db->set_legendTitle( "300 < p_{T} < 400 GeV");
  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nChargedJet0_HT600", "Lead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet0", selection_pt0.c_str(), 51, -0.5, 50.5, "nNeutralJet0_HT600", "Lead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet0", selection_pt0.c_str(),      50, 0., 1.0001, "ptDJet0_HT600", "Lead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet0", selection_pt0.c_str(),      50, 0., 1.0001, "QGLikelihoodJet0_HT600", "Lead Jet Q-G LD", "", "Events", true);

  db->drawHisto_fromTree( "tree_passedEvents", "nChargedJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nChargedJet1_HT600", "Sublead Jet Charged Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "nNeutralJet1", selection_pt1.c_str(), 51, -0.5, 50.5, "nNeutralJet1_HT600", "Sublead Jet Neutral Multiplicity", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "ptDJet1", selection_pt1.c_str(),      50, 0., 1.0001, "ptDJet1_HT600", "Sublead Jet p_{T}D", "", "Events", true);
  db->drawHisto_fromTree( "tree_passedEvents", "QGLikelihoodJet1", selection_pt1.c_str(),      50, 0., 1.0001, "QGLikelihoodJet1_HT600", "Sublead Jet Q-G LD", "", "Events", true);

  db->set_legendTitle("");

/*
  db->set_rebin(2);
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet0", "First Jet Q-G Likelihood", "", "Events");
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet1", "Second Jet Q-G Likelihood", "", "Events");
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet2", "Third Jet Q-G Likelihood", "", "Events");
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet3", "Fourth Jet Q-G Likelihood", "", "Events");

  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet0", "First Jet Q-G Likelihood", "", "Events", false, 100., 200.);
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet0", "First Jet Q-G Likelihood", "", "Events", false, 200., 300.);
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet0", "First Jet Q-G Likelihood", "", "Events", false, 300., 400.);

  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet1", "Second Jet Q-G Likelihood", "", "Events", false, 100., 200.);
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet1", "Second Jet Q-G Likelihood", "", "Events", false, 200., 300.);
  drawHistoWithQuarkGluonComponents( db, "tree_passedEvents", "QGLikelihoodJet1", "Second Jet Q-G Likelihood", "", "Events", false, 300., 400.);
*/

  delete db;
  db = 0;

  return 0;

}  



void drawHistoWithQuarkGluonComponents( DrawBase* db, const std::string& treeName, const std::string& varName, const std::string& axisName, const std::string& units, const std::string& instanceName, bool log, float ptMin, float ptMax ) {

  db->drawHisto( varName, axisName, units, instanceName, log );

  TString varName_tstr(varName);

  int jetNumber=-1;
  if( varName_tstr.Contains("Jet0") )
    jetNumber=0;
  else if( varName_tstr.Contains("Jet1") )
    jetNumber=1;
  else if( varName_tstr.Contains("Jet2") )
    jetNumber=2;
  else if( varName_tstr.Contains("Jet3") )
    jetNumber=3;

  if( jetNumber<0 ) {
    std::cout << "There must be a problem, this is not possible." << std::endl;
    exit(333);
  }

  //TH1D* h1_data = db->get_lastHistos_data()[0];
  std::vector<TH1D*> lastHistos = db->get_lastHistos_mc();

  int nBins = lastHistos[0]->GetNbinsX();
  float xMin = lastHistos[0]->GetXaxis()->GetXmin();
  float xMax = lastHistos[0]->GetXaxis()->GetXmax();

  TH1D* h1_data = new TH1D( "data", "", nBins, xMin, xMax );
  TH1D* h1_all = new TH1D( "all", "", nBins, xMin, xMax );
  TH1D* h1_quark = new TH1D( "quark", "", nBins, xMin, xMax );
  TH1D* h1_gluon = new TH1D( "gluon", "", nBins, xMin, xMax );
  TH1D* h1_b = new TH1D( "b", "", nBins, xMin, xMax );


  //char allCondition[400];
  //sprintf( allCondition, "eventWeight*(pdgIdPartJet%d>-10 && ptJet%d>%f && ptJet%d<%f && ptDJet%d<0.9)", jetNumber, jetNumber, ptMin, jetNumber, ptMax, jetNumber );
  //char quarkCondition[400];
  //sprintf( quarkCondition, "eventWeight*(abs(pdgIdPartJet%d)<5 && ptJet%d>%f && ptJet%d<%f && ptDJet%d<0.9)", jetNumber, jetNumber, ptMin, jetNumber, ptMax, jetNumber );
  //char gluonCondition[400];
  //sprintf( gluonCondition, "eventWeight*(pdgIdPartJet%d==21 && ptJet%d>%f && ptJet%d<%f && ptDJet%d<0.9)", jetNumber, jetNumber, ptMin, jetNumber, ptMax, jetNumber );

  char commonCondition[200];
  //sprintf( commonCondition, "ptJet%d>%f && ptJet%d<%f", jetNumber, ptMin, jetNumber, ptMax ); 
  sprintf( commonCondition, "ptJet%d>%f && ptJet%d<%f && ptDJet%d<0.9", jetNumber, ptMin, jetNumber, ptMax, jetNumber ); 

  char allCondition[400];
  sprintf( allCondition, "eventWeight*(%s)", commonCondition );
  char quarkCondition[400];
  sprintf( quarkCondition, "eventWeight*(abs(pdgIdPartJet%d)<5 && %s)", jetNumber, commonCondition );
  char gluonCondition[400];
  sprintf( gluonCondition, "eventWeight*(pdgIdPartJet%d==21 && %s)", jetNumber, commonCondition );
  char bCondition[400];
  sprintf( bCondition, "eventWeight*(abs(pdgIdPartJet%d)==5 && %s)", jetNumber, commonCondition );

  TTree* treeDATA = (TTree*)(db->get_dataFile(0).file->Get(treeName.c_str()));
  treeDATA->Project( "data", varName.c_str(), commonCondition );


  TTree* treeMC = (TTree*)(db->get_mcFile(0).file->Get(treeName.c_str()));
  treeMC->Project( "all",   varName.c_str(), allCondition );
  treeMC->Project( "quark", varName.c_str(), quarkCondition );
  treeMC->Project( "gluon", varName.c_str(), gluonCondition );
  treeMC->Project( "b", varName.c_str(), bCondition );

  float data_int = h1_data->Integral();
  float mc_int = h1_all->Integral();
  float scaleFactor = data_int/mc_int;

  float quark_fraction = h1_quark->Integral()/mc_int;
  float gluon_fraction = h1_gluon->Integral()/mc_int;
  float b_fraction = h1_b->Integral()/mc_int;
  float other_fraction = 1.-quark_fraction-gluon_fraction-b_fraction;

  char quarkText[300];
  sprintf( quarkText, "udsc (%.1f%%)", 100.*quark_fraction );
  char gluonText[300];
  sprintf( gluonText, "Gluons (%.1f%%)", 100.*gluon_fraction );
  char bText[300];
  sprintf( bText, "b (%.1f%%)", 100.*b_fraction );
  char otherText[300];
  sprintf( otherText, "Undefined (%.1f%%)", 100.*other_fraction );

  
  TLegend* legend;
  if( ptMin !=0. && ptMax != 10000. ) {
    char legendTitle[150];
    sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
    legend = new TLegend( 0.5, 0.55, 0.88, 0.9, legendTitle );
  } else {
    legend = new TLegend( 0.5, 0.6, 0.88, 0.9 );
  }
  legend->SetFillColor( kWhite );
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_data, "Data", "p" );
  legend->AddEntry( h1_gluon, gluonText, "F" );
  legend->AddEntry( h1_quark, quarkText, "F" );
  legend->AddEntry( h1_b, bText, "F" );
  legend->AddEntry( h1_all, otherText, "F" );

//h1_all->Rebin( 2 );
//h1_gluon->Rebin( 2 );
//h1_quark->Rebin( 2 );
//h1_b->Rebin( 2 );
//h1_data->Rebin( 2 );
  
  h1_all->Scale( scaleFactor );
  h1_gluon->Scale( scaleFactor );
  h1_quark->Scale( scaleFactor );
  h1_b->Scale( scaleFactor );
  
  h1_data->SetMarkerStyle( 20 );
  h1_data->SetMarkerSize( 1. );
  h1_all->SetFillColor( kGray );
  h1_gluon->SetFillColor( 46 );
  h1_quark->SetFillColor( 38 );
  h1_b->SetFillColor( kYellow );

  THStack* stack = new THStack();
  stack->Add(h1_gluon );
  stack->Add(h1_quark);
  stack->Add(h1_b);

  float dataMax = h1_data->GetMaximum();
  float mcMax = h1_all->GetMaximum();
  float yMax = (dataMax>mcMax) ? dataMax : mcMax;
  yMax *= 1.3;


  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  legend->Draw("same");
  h1_all->Draw("same");
  stack->Draw("histo same");
  h1_data->Draw("p same");
  sqrtLabel->Draw("Same");

  gPad->RedrawAxis();

  //std::string canvasName = db->get_outputdir() + "/" + varName + "_components.eps";

  char canvasName[200];
  if( ptMin !=0. && ptMax != 10000. )
    sprintf( canvasName, "%s/%s_pt%.0f_%.0f_components.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );
  else
    sprintf( canvasName, "%s/%s_components.eps", db->get_outputdir().c_str(), varName.c_str() );

  c1->SaveAs(canvasName);

  
  delete c1;
  delete h2_axes;

  delete h1_data;
  delete h1_all;
  delete h1_quark;
  delete h1_gluon;
  delete h1_b;
  
}





std::pair<TH1D, TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db ) {

  TFile* fileMC_photonjet = db->get_mcFile(0).file;
  TFile* fileMC_qcd = db->get_mcFile(1).file;

  if( fileMC_qcd==0 ) {
    std::cout << "Didn't find QCD file. Exiting." << std::endl;
    exit(199);
  }

  TFile* file_data = db->get_dataFile(0).file;

  char histoName[200];
  sprintf( histoName, "%s_%d%d", varName.c_str(), ptMin, ptMax);

  TH1D* h1_photonjet = (TH1D*)fileMC_photonjet->Get(histoName);
  TH1D* h1_qcd = (TH1D*)fileMC_qcd->Get(histoName);
  TH1D* h1_data = (TH1D*)file_data->Get(histoName);

  int nBins = h1_qcd->GetXaxis()->GetNbins();


  char histoName_quark[200];
  sprintf( histoName_quark, "%s_quark_noPhotID_%d%d", varName.c_str(), ptMin, ptMax);
  char histoName_gluon[200];
  sprintf( histoName_gluon, "%s_gluon_noPhotID_%d%d", varName.c_str(), ptMin, ptMax);

  //TH1D* h1_quark = (TH1D*)fileMC_qcd->Get(histoName_quark);
  //TH1D* h1_gluon = (TH1D*)fileMC_qcd->Get(histoName_gluon);
  TH1D* h1_quark = (TH1D*)fileMC_photonjet->Get(histoName_quark);
  TH1D* h1_gluon = (TH1D*)fileMC_photonjet->Get(histoName_gluon);

  // same area:
  float quark_integral = h1_quark->Integral(1, nBins);
  float gluon_integral = h1_gluon->Integral(1, nBins);
  h1_quark->Scale(1./quark_integral);
  h1_gluon->Scale(1./gluon_integral);


  char histoName_quarkFraction[200];
  sprintf( histoName_quarkFraction, "quarkFraction_%d%d", ptMin, ptMax );
  //sprintf( histoName_quarkFraction, "quarkFraction_antibtag_%d%d", ptMin, ptMax );

  TH1D* h1_quarkFraction_qcd = (TH1D*)fileMC_qcd->Get(histoName_quarkFraction);
  TH1D* h1_quarkFraction_photonjet = (TH1D*)fileMC_photonjet->Get(histoName_quarkFraction);

  std::cout <<  h1_quarkFraction_qcd << std::endl;
  std::cout <<  h1_quarkFraction_photonjet << std::endl;

  float quarkFraction_qcd = h1_quarkFraction_qcd->GetBinContent(1);
  float quarkFraction_photonjet = h1_quarkFraction_photonjet->GetBinContent(1);

  float qcd_integral = h1_qcd->Integral(1, nBins);
  float photonjet_integral = h1_photonjet->Integral(1, nBins);
  float data_integral = h1_data->Integral(1, nBins);

  TH1D* h1_qcd_fromFractions = new TH1D(*h1_qcd);

  for( unsigned iBin=1; iBin<nBins+1; ++iBin ) {

    float thisBinValue = quarkFraction_qcd*h1_quark->GetBinContent(iBin) + (1.-quarkFraction_qcd)*h1_gluon->GetBinContent(iBin);
    h1_qcd_fromFractions->SetBinContent( iBin, thisBinValue );

  } //for bins


  float qcd_fromFraction_integral = h1_qcd_fromFractions->Integral(1, nBins);
  h1_qcd_fromFractions->Scale( 1./qcd_fromFraction_integral );
  h1_photonjet->Scale( 1./photonjet_integral );

  float ratio = photonjet_integral/qcd_integral;
  h1_qcd_fromFractions->Scale( data_integral/(ratio + 1.) );
  h1_photonjet->Scale( data_integral*ratio/(ratio + 1.) );

  h1_photonjet->SetFillColor( 46 );
  h1_qcd_fromFractions->SetFillColor( 38 );
  h1_qcd->SetFillColor( kRed );
  h1_qcd->SetLineColor( kRed );
  h1_qcd->SetFillStyle( 3004 );

  h1_data->SetMarkerStyle(20);
  h1_data->SetMarkerSize(1.2);

  int rebin = 2;
  if( ptMin==30 ) rebin = 5;

  h1_quark->Rebin(rebin);
  h1_gluon->Rebin(rebin);
  h1_photonjet->Rebin(rebin);
  h1_qcd->Rebin(rebin);
  h1_qcd_fromFractions->Rebin(rebin);
  h1_data->Rebin(rebin);

  nBins = h1_quark->GetXaxis()->GetNbins();

  THStack* mcstack = new THStack();
  mcstack->Add( h1_qcd_fromFractions );
  mcstack->Add( h1_photonjet );

  float data_max = h1_data->GetMaximum();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();


  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., 1.5*data_max);
  h2_axes->SetXTitle("Q-G Likelihood");
  h2_axes->SetYTitle("Events");


  TLegend* legend_qcd = new TLegend( 0.2, 0.75, 0.5, 0.9 );
  legend_qcd->SetFillColor(0);
  legend_qcd->SetTextSize(0.04);
  legend_qcd->AddEntry( h1_qcd, "QCD MC", "F" );
  legend_qcd->AddEntry( h1_qcd_fromFractions, "QCD From Comp", "F" );

  
  //h2_axes->Draw(); 

  h1_qcd_fromFractions->DrawNormalized("histo");
  h1_qcd->DrawNormalized("histo same");

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  legend_qcd->Draw("same");

  gPad->RedrawAxis();

  char canvasName_qcd[500];
  sprintf( canvasName_qcd, "%s/compare_qcd_%s_%d%d.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );
  
  c1->SaveAs( canvasName_qcd );
  
  c1->Clear();

  char legendTitle[200];
  sprintf( legendTitle, "%d < p_{T} < %d GeV", ptMin , ptMax );
  TLegend* legend = new TLegend( 0.2, 0.65, 0.5, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_data, "Data", "P" );
  legend->AddEntry( h1_photonjet, "#gamma+Jet MC", "F" );
  legend->AddEntry( h1_qcd_fromFractions, "QCD MC", "F" );
  

  h2_axes->Draw(); 

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  legend->Draw("same");

  mcstack->Draw("histo same");

  h1_data->Draw("E same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/%s_%d%d_fromComponents.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );

  c1->SaveAs(canvasName);


  // and now subtract gluon component both from data and MC


  TH1D* h1_photonjet_quarkOnly = new TH1D(*h1_quark);
  h1_photonjet_quarkOnly->SetName("photonjet_quarkOnly");
  float photonjet_integral_quarkOnly = h1_photonjet_quarkOnly->Integral(1, nBins);
  h1_photonjet_quarkOnly->Scale( 1./photonjet_integral_quarkOnly );
  h1_photonjet_quarkOnly->Scale( quarkFraction_photonjet*h1_photonjet->Integral(1,nBins) );

  TH1D* h1_qcd_quarkOnly = new TH1D(*h1_quark);
  h1_qcd_quarkOnly->SetName("qcd_quarkOnly");
  float qcd_integral_quarkOnly = h1_qcd_quarkOnly->Integral(1, nBins);
  h1_qcd_quarkOnly->Scale( 1./qcd_integral_quarkOnly );
  h1_qcd_quarkOnly->Scale( quarkFraction_qcd*h1_qcd_fromFractions->Integral(1,nBins) );

  TH1D* h1_photonjet_gluonOnly = new TH1D(*h1_gluon);
  h1_photonjet_gluonOnly->SetName("photonjet_gluonOnly");
  float photonjet_integral_gluonOnly = h1_photonjet_gluonOnly->Integral(1, nBins);
  h1_photonjet_gluonOnly->Scale( 1./photonjet_integral_gluonOnly );
  h1_photonjet_gluonOnly->Scale( (1.-quarkFraction_photonjet)*h1_photonjet->Integral(1,nBins) );

  TH1D* h1_qcd_gluonOnly = new TH1D(*h1_gluon);
  h1_qcd_gluonOnly->SetName("qcd_gluonOnly");
  float qcd_integral_gluonOnly = h1_qcd_gluonOnly->Integral(1, nBins);
  h1_qcd_gluonOnly->Scale( 1./qcd_integral_gluonOnly );
  h1_qcd_gluonOnly->Scale( (1.-quarkFraction_qcd)*h1_qcd_fromFractions->Integral(1,nBins) );

//TH1D* h1_qcd_quark = new TH1D(*h1_qcd);
//h1_qcd_quark->SetName("qcd_quark");
//h1_qcd_quark->Add( h1_quark, -quarkFraction_qcd );
////h1_qcd_quark->Add( h1_quark, -quarkFraction_qcd*(h1_qcd_quark->Integral(1, nBins)) );
//TH1D* h1_qcd_gluon = new TH1D(*h1_qcd);
//h1_qcd_gluon->SetName("qcd_gluon");
//h1_qcd_gluon->Add( h1_gluon, -(1.-quarkFraction_qcd) );
////h1_qcd_gluon->Add( h1_gluon, -(1.-quarkFraction_qcd)*(h1_qcd_gluon->Integral(1, nBins)) );


  TH1D* h1_mc_gluon_only = new TH1D(*h1_photonjet_gluonOnly);
  h1_mc_gluon_only->SetName("mc_gluon_only");
  h1_mc_gluon_only->Add( h1_qcd_gluonOnly );

  h1_data->SetName("data");


  // subtract from mc:
//  TH1D* h1_mc_quark_only = new TH1D("mc_quark_only", "", nBins, 0., 1.0001);
  TH1D* h1_mc_quark_only = new TH1D(*h1_photonjet);
  h1_mc_quark_only->Add(h1_qcd_fromFractions);
//  h1_mc_quark_only = (mcstack->GetHistogram());
  h1_mc_quark_only->SetName("mc_quark_only");
  h1_mc_quark_only->Add( h1_mc_gluon_only, -1. );

  // subtract from data:
  TH1D* h1_data_quark_only = new TH1D(*h1_data);
  h1_data_quark_only->SetName("data_quark_only");
  h1_data_quark_only->Add( h1_mc_gluon_only, -1. );


  //h1_mc_quark_only->Scale( h1_data_quark_only->Integral(1,nBins)/h1_mc_quark_only->Integral(1,nBins) );

TFile* file_prova = TFile::Open("prova.root", "recreate");
file_prova->cd();
h1_qcd_quarkOnly->Write();
h1_qcd_gluonOnly->Write();
h1_photonjet_quarkOnly->Write();
h1_photonjet_gluonOnly->Write();
  h1_mc_gluon_only->Write();
  h1_mc_quark_only->Write();
h1_data->Write();
h1_data_quark_only->Write();
file_prova->Close();
  h1_mc_quark_only->SetFillColor(46);

  h1_data_quark_only->SetMarkerStyle(20);
  h1_data_quark_only->SetMarkerSize(1.1);


  c1->Clear();

  h2_axes->Draw();

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  TLegend* legend_quark = new TLegend(0.2, 0.65, 0.5, 0.9, legendTitle );
  legend_quark->SetFillColor(0);
  legend_quark->SetTextSize(0.04);
  legend_quark->AddEntry( h1_data_quark_only, "Data (gluon subtracted)", "P" );
  legend_quark->AddEntry( h1_mc_quark_only, "MC (gluon subtracted)", "F" );
  legend_quark->Draw("same");

  h1_mc_quark_only->Draw("histo same"); 
  h1_data_quark_only->Draw("E same"); 

  gPad->RedrawAxis();

  char canvasName_quark[500];
  sprintf( canvasName_quark, "%s/%s_%d%d_quarkOnly.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );

  c1->SaveAs( canvasName_quark );

  
  //check efficiencies:

  int iBin_cut = h1_data_quark_only->FindBin(0.2);

  float eff_mc   = h1_mc_quark_only->Integral(iBin_cut+1, nBins)/h1_mc_quark_only->Integral(1, nBins);
  float eff_data = h1_data_quark_only->Integral(iBin_cut+1, nBins)/h1_data_quark_only->Integral(1, nBins);

  float effErr_mc   = sqrt( eff_mc*(1.-eff_mc)/h1_mc_quark_only->GetEntries());
  float effErr_data = sqrt( eff_data*(1.-eff_data)/h1_data_quark_only->GetEntries());

  std::cout << std::endl << "*** pt bin: " << ptMin << "-" << ptMax << std::endl;
  std::cout << "requiring QGLikelihood>0.2" << std::endl;
  std::cout << "eff(mc): " << eff_mc << " +/- " << effErr_mc << std::endl;
  std::cout << "eff(data): " << eff_data << " +/- " << effErr_data << std::endl;
  std::cout << std::endl;

  std::pair< TH1D, TH1D > returnPair;
  returnPair.first = *h1_mc_quark_only;
  returnPair.second = *h1_data_quark_only;

  return returnPair;

}



//void drawSignalPtMix( std::pair<TH1D*,TH1D*> h1pair_3050, std::pair<TH1D*,TH1D*> h1pair_5080, std::pair<TH1D*,TH1D*> h1pair_80120, int mass, float frac_3050, float frac_5080, float frac_80120, DrawBase* db ) {
//
//
//
//
//}
