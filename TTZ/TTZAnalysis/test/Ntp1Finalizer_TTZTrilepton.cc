#include "Ntp1Finalizer_TTZTriplepton.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "CommonTools/fitTools.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"
#include "KinematicFit/DiJetKinFitter.h"

#include "PUWeight.h"




bool USE_MC_MASS=false;

int DEBUG_EVENTNUMBER = 98901397;






// constructor:

Ntp1Finalizer_TTZTriplepton::Ntp1Finalizer_TTZTriplepton( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, const std::string& PUType, const std::string& leptType ) : Ntp1Finalizer( "TTZTriplepton", dataset, leptType ) {

  if( leptType!="ALL" && leptType!="MU" && leptType!="ELE" ) {
    std::cout << "Lept type '" << leptType << "' currently not supported. Exiting." << std::endl;
    exit(9177);
  }

  if( bTaggerType!="SSVHE" && bTaggerType!="TCHE" ) {
    std::cout << "b-Tagger type '" << bTaggerType << "' currently not supported. Exiting." << std::endl;
    exit(9179);
  }

  bTaggerType_ = bTaggerType;
  leptType_ = leptType;
  PUType_ = PUType;

  setSelectionType(selectionType);

}




void Ntp1Finalizer_TTZTriplepton::finalize() {

  //if( outFile_==0 ) this->createOutputFile();
  
  Int_t run;
  tree_->SetBranchAddress("run", &run);
  tree_->GetEntry(0);
  bool isMC = (run < 10);
  std::string fullFlags = selectionType_ + "_" + bTaggerType_;
  if( isMC ) fullFlags = fullFlags + "_PU" + PUType_;
  fullFlags = fullFlags + "_" + leptType_;
  this->set_flags(fullFlags); //this is for the outfile name
  this->createOutputFile();



  TTree* tree_passedEvents = new TTree("tree_passedEvents", "Unbinned data for statistical treatment");

  TH1D* h1_nCounter = new TH1D("nCounter", "", 1, 0., 1.);
  h1_nCounter->Sumw2();
  TH1D* h1_nCounterW = new TH1D("nCounterW", "", 1, 0., 1.);
  h1_nCounterW->Sumw2();
  TH1D* h1_nCounterPU = new TH1D("nCounterPU", "", 1, 0., 1.);
  h1_nCounterPU->Sumw2();




  TH1D* h1_nvertex = new TH1D("nvertex", "", 36, -0.5, 35.5);
  h1_nvertex->Sumw2();
  TH1D* h1_nvertex_PUW = new TH1D("nvertex_PUW", "", 36, -0.5, 35.5);
  h1_nvertex_PUW->Sumw2();
  TH1D* h1_nvertex_PUW_ave = new TH1D("nvertex_PUW_ave", "", 36, -0.5, 35.5);
  h1_nvertex_PUW_ave->Sumw2();

  TH1D* h1_pfMet = new TH1D("pfMet", "", 500, 0., 500.);
  h1_pfMet->Sumw2();

  TH1D* h1_metSignificance= new TH1D("metSignificance", "", 80, 0., 40.);
  h1_metSignificance->Sumw2();

  TH1D* h1_mEtSig= new TH1D("mEtSig", "", 60, 0., 15.);
  h1_mEtSig->Sumw2();


  TH1D* h1_rhoPF_presel = new TH1D("rhoPF_presel", "", 50, 0., 20.);
  h1_rhoPF_presel->Sumw2();
  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., 20.);
  h1_rhoPF->Sumw2();


  TH1D* h1_ptLeptZ1 = new TH1D("ptLeptZ1", "", 500, 20., 520.);
  h1_ptLeptZ1->Sumw2();
  TH1D* h1_ptLeptZ2 = new TH1D("ptLeptZ2", "", 200, 20., 220.);
  h1_ptLeptZ2->Sumw2();
  TH1D* h1_etaLeptZ1 = new TH1D("etaLeptZ1", "", 50, -2.5, 2.5);
  h1_etaLeptZ1->Sumw2();
  TH1D* h1_etaLeptZ2 = new TH1D("etaLeptZ2", "", 50, -2.5, 2.5);
  h1_etaLeptZ2->Sumw2();

  TH1D* h1_ptLept3 = new TH1D("ptLept3", "", 200, 20., 220.);
  h1_ptLept3->Sumw2();
  TH1D* h1_etaLept3 = new TH1D("etaLept3", "", 50, -2.5, 2.5);
  h1_etaLept3->Sumw2();


  TH1D* h1_nJets_presel = new TH1D("nJets_presel", "", 7, 1.5, 8.5);
  h1_nJets_presel->Sumw2();


  TH1D* h1_bTagJet1B = new TH1D("bTagJet1B", "", 420, 0., 20.);
  h1_bTagJet1B->Sumw2();
  TH1D* h1_bTagJet2B = new TH1D("bTagJet2B", "", 420, -1., 20.);
  h1_bTagJet2B->Sumw2();




  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 100, -2.4, 2.4);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 100, -2.4, 2.4);
  h1_etaJet2->Sumw2();






  Int_t nPU;
  tree_->SetBranchAddress("nPU", &nPU);
  Int_t nvertex;
  tree_->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree_->SetBranchAddress("rhoPF", &rhoPF);
  Int_t LS;
  tree_->SetBranchAddress("LS", &LS);
  unsigned int event;
  tree_->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree_->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeightPU;
  tree_->SetBranchAddress("eventWeightPU", &eventWeightPU);
  Float_t eventWeightPU_ave;
  tree_->SetBranchAddress("eventWeightPU_ave", &eventWeightPU_ave);
  Float_t eventWeight_Zee;
  tree_->SetBranchAddress("eventWeight_Zee", &eventWeight_Zee);
  Float_t eventWeight_Zmm;
  tree_->SetBranchAddress("eventWeight_Zmm", &eventWeight_Zmm);

  Float_t ptHat;
  tree_->SetBranchAddress("ptHat", &ptHat);

  Float_t pfMet;
  tree_->SetBranchAddress("epfMet", &pfMet);
  Float_t metSignificance;
  tree_->SetBranchAddress("metSignificance", &metSignificance);
  Float_t mEtSig;
  tree_->SetBranchAddress("mEtSig", &mEtSig);
  Float_t phiMet;
  tree_->SetBranchAddress("phipfMet", &phiMet);


  int leptType;
  tree_->SetBranchAddress("leptType", &leptType);

  Float_t eLeptZ1;
  tree_->SetBranchAddress("eLeptZ1", &eLeptZ1);
  Float_t ptLeptZ1;
  tree_->SetBranchAddress("ptLeptZ1", &ptLeptZ1);
  Float_t etaLeptZ1;
  tree_->SetBranchAddress("etaLeptZ1", &etaLeptZ1);
  Float_t phiLeptZ1;
  tree_->SetBranchAddress("phiLeptZ1", &phiLeptZ1);
  Int_t chargeLeptZ1;
  tree_->SetBranchAddress("chargeLeptZ1", &chargeLeptZ1);

  Float_t eLeptZ2;
  tree_->SetBranchAddress("eLeptZ2", &eLeptZ2);
  Float_t ptLeptZ2;
  tree_->SetBranchAddress("ptLeptZ2", &ptLeptZ2);
  Float_t etaLeptZ2;
  tree_->SetBranchAddress("etaLeptZ2", &etaLeptZ2);
  Float_t phiLeptZ2;
  tree_->SetBranchAddress("phiLeptZ2", &phiLeptZ2);
  Int_t chargeLeptZ2;
  tree_->SetBranchAddress("chargeLeptZ2", &chargeLeptZ2);


  Int_t nLept;
  tree_->SetBranchAddress("nLept", &nLept);
  Int_t leptTypeLept[10];
  tree_->SetBranchAddress("leptTypeLept", leptTypeLept);
  Float_t eLept[10];
  tree_->SetBranchAddress("eLept", eLept);
  Float_t ptLept[10];
  tree_->SetBranchAddress("ptLept", ptLept);
  Float_t etaLept[10];
  tree_->SetBranchAddress("etaLept", etaLept);
  Float_t phiLept[10];
  tree_->SetBranchAddress("phiLept", phiLept);
  Int_t chargeLept[10];
  tree_->SetBranchAddress("chargeLept", chargeLept);


  Int_t nJets;
  tree_->SetBranchAddress("nJets", &nJets);

  Float_t eJet[50];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[50];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t etaJet[50];
  tree_->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[50];
  tree_->SetBranchAddress("phiJet", phiJet);
  Float_t eJetGen[50];
  tree_->SetBranchAddress("eJetGen", eJetGen);
  Float_t ptJetGen[50];
  tree_->SetBranchAddress("ptJetGen", ptJetGen);
  Float_t etaJetGen[50];
  tree_->SetBranchAddress("etaJetGen", etaJetGen);
  Float_t phiJetGen[50];
  tree_->SetBranchAddress("phiJetGen", phiJetGen);
  Float_t eChargedHadronsJet[50];
  tree_->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t rmsCandJet[50];
  tree_->SetBranchAddress("rmsCandJet", rmsCandJet);
  Float_t ptDJet[50];
  tree_->SetBranchAddress("ptDJet", ptDJet);
  Int_t nChargedJet[50];
  tree_->SetBranchAddress("nChargedJet", nChargedJet);
  Int_t nNeutralJet[50];
  tree_->SetBranchAddress("nNeutralJet", nNeutralJet);
  Float_t eMuonsJet[50];
  tree_->SetBranchAddress("eMuonsJet", eMuonsJet);
  Float_t eElectronsJet[50];
  tree_->SetBranchAddress("eElectronsJet", eElectronsJet);
  Float_t trackCountingHighEffBJetTagJet[50];
  tree_->SetBranchAddress("trackCountingHighEffBJetTagJet", trackCountingHighEffBJetTagJet);
  Float_t trackCountingHighPurBJetTagJet[50];
  tree_->SetBranchAddress("trackCountingHighPurBJetTagJet", trackCountingHighPurBJetTagJet);
  Float_t simpleSecondaryVertexHighEffBJetTagJet[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagJet", simpleSecondaryVertexHighEffBJetTagJet);
  Float_t simpleSecondaryVertexHighPurBJetTagJet[50];
  tree_->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagJet", simpleSecondaryVertexHighPurBJetTagJet);
  Float_t jetBProbabilityBJetTagJet[50];
  tree_->SetBranchAddress("jetBProbabilityBJetTagJet", jetBProbabilityBJetTagJet);
  Float_t jetProbabilityBJetTagJet[50];
  tree_->SetBranchAddress("jetProbabilityBJetTagJet", jetProbabilityBJetTagJet);



  Int_t nPart;
  tree_->SetBranchAddress("nPart", &nPart);
  Float_t ePart[20];
  tree_->SetBranchAddress("ePart", ePart);
  Float_t ptPart[20];
  tree_->SetBranchAddress("ptPart", ptPart);
  Float_t etaPart[20];
  tree_->SetBranchAddress("etaPart", etaPart);
  Float_t phiPart[20];
  tree_->SetBranchAddress("phiPart", phiPart);
  Int_t pdgIdPart[20];
  tree_->SetBranchAddress("pdgIdPart", pdgIdPart);


  // HLT:
  Bool_t passed_HLT_DoubleMu6;
  tree_->SetBranchAddress("passed_HLT_DoubleMu6", &passed_HLT_DoubleMu6);
  Bool_t passed_HLT_DoubleMu7;
  tree_->SetBranchAddress("passed_HLT_DoubleMu7", &passed_HLT_DoubleMu7);
  Bool_t passed_HLT_Mu13_Mu8;
  tree_->SetBranchAddress("passed_HLT_Mu13_Mu8", &passed_HLT_Mu13_Mu8);
  Bool_t passed_HLT_IsoMu17;
  tree_->SetBranchAddress("passed_HLT_IsoMu17", &passed_HLT_IsoMu17);
  Bool_t passed_HLT_IsoMu24;
  tree_->SetBranchAddress("passed_HLT_IsoMu24", &passed_HLT_IsoMu24);
  Bool_t passed_HLT_Mu8_Jet40;
  tree_->SetBranchAddress("passed_HLT_Mu8_Jet40", &passed_HLT_Mu8_Jet40);
  Bool_t passed_HLT_L2DoubleMu23_NoVertex;
  tree_->SetBranchAddress("passed_HLT_L2DoubleMu23_NoVertex", &passed_HLT_L2DoubleMu23_NoVertex);
  Bool_t passed_HLT_L2DoubleMu30_NoVertex;
  tree_->SetBranchAddress("passed_HLT_L2DoubleMu30_NoVertex", &passed_HLT_L2DoubleMu30_NoVertex);
  Bool_t passed_HLT_TripleMu5;
  tree_->SetBranchAddress("passed_HLT_TripleMu5", &passed_HLT_TripleMu5);

  Bool_t passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
  tree_->SetBranchAddress("passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", &passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  Bool_t passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  tree_->SetBranchAddress("passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL);
  





  int nEntries = tree_->GetEntries();
  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;


  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
  float Zmass = 91.1876;


  std::string puType = "Spring11_Flat10";
  std::string puType_ave = "Spring11_Flat10";
  TString dataset_tstr(dataset_);
  if( dataset_tstr.Contains("Summer11") && dataset_tstr.Contains("PU_S4") ) {
    puType = "Summer11_S4";
    puType_ave = "Summer11_S4_ave";
  } else if( dataset_tstr.Contains("Fall11") ) {
    puType = "Fall11";
  }
  PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
  PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
  std::string puFileName;
  if( PUType_=="HR11" || PUType_=="HR11_v3")
    puFileName = "Pileup_DATA_up_to_178479_SiXie.root";
    //puFileName = "Pileup_DATA_up_to_178078.root";
  else if( PUType_=="Run2011A" )
    puFileName = "Pileup_DATA_up_to_173692.root";
  else if( PUType_=="Run2011A_73pb" )
    puFileName = "all2011A.pileup_v2_73mb.root";
  else if( PUType_=="Run2011B" )
    puFileName = "Pileup_DATA_Run2011B.root";
    //puFileName = "Pileup_DATA_173692_to_178078.root";
  else if( PUType_=="Run2011B_73pb" )
    puFileName = "all2011B.pileup_v2_73mb.root";
  else if( PUType_=="HR11_73pb" || PUType_=="HR11_73pb_DY" )
    puFileName = "all2011AB.pileup_v2_73mb.root";
  else if( PUType_!="HR11_v2" ) {
    std::cout << "-> Unknown PU Type: '" << PUType_ << "'. Will use HR11 default." << std::endl;
    puFileName = "Pileup_DATA_up_to_178078.root";
  }


  if( PUType_!="HR11_v2" ) {
    std::cout << std::endl << "-> Using data pileup file: " << puFileName << std::endl;
    TFile* filePU = TFile::Open(puFileName.c_str());
    TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
    fPUWeight->SetDataHistogram(h1_nPU_data);
    fPUWeight_ave->SetDataHistogram(h1_nPU_data);
  } else {  // HR11_v2: 4.6fb-1 = 2.1 (A) + 2.5 (B)
    TFile* filePU_RunA = TFile::Open("all2011A.pileup_v2_73mb.root");
    TFile* filePU_RunB = TFile::Open("all2011B.pileup_v2_73mb.root");
    TH1F* h1_PURunA = (TH1F*)filePU_RunA->Get("pileup");
    TH1F* h1_PURunB = (TH1F*)filePU_RunB->Get("pileup");
    h1_PURunA->Scale(2.1/h1_PURunA->Integral());
    h1_PURunB->Scale(2.5/h1_PURunB->Integral());
    TH1F* h1_PU_weightedAverage = new TH1F(*h1_PURunA);
    h1_PU_weightedAverage->Add(h1_PURunB);
    fPUWeight->SetDataHistogram(h1_PU_weightedAverage);
    fPUWeight_ave->SetDataHistogram(h1_PU_weightedAverage);
  }
    
     


  if( PUType_=="HR11_73pb_DY" ) {
    TFile* filePUMC = TFile::Open("generatedpileup_Zjets_MADGRAPH_AOD423.root");
    TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("GenLevelInfoModule/npileup");
    std::cout << "-> Switching MC PU file to: generatedpileup_Zjets_MADGRAPH_AOD423.root" << std::endl;
    fPUWeight->SetMCHistogram(h1_nPU_mc);
  } else if( dataset_tstr.Contains("Summer11") && dataset_tstr.Contains("PU_S4") && PUType_!="HR11_v3" ) {
    TFile* filePUMC = TFile::Open("Pileup_MC_Summer11_S4.root");
    TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("hNPU");
    std::cout << "-> Switching MC PU file to: Pileup_MC_Summer11_S4.root" << std::endl;
    fPUWeight->SetMCHistogram(h1_nPU_mc);
  } else if( dataset_tstr.Contains("Fall11") ) {
    TFile* filePUMC = TFile::Open("s6MCPileUp.root");
    TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("pileup");
    std::cout << "-> Switching MC PU file to: s6MCPileUp.root" << std::endl;
    fPUWeight->SetMCHistogram(h1_nPU_mc);
  }



  int maxBTag_found = -1;
  float mZll;
  float ptLeptZ1_t, ptLeptZ2_t, etaLeptZ1_t, etaLeptZ2_t;
  float ptLept3_t, etaLept3_t;
  float ptJetB1_t, ptJetB2_t, etaJetB1_t, etaJetB2_t;
  float ptJet3_t, ptJet4_t, etaJet3_t, etaJet4_t;
  float HLTSF;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "leptType", &leptType, "leptType/I" );
  tree_passedEvents->Branch( "ptLeptZ1", &ptLeptZ1_t, "ptLeptZ1_t/F" );
  tree_passedEvents->Branch( "ptLeptZ2", &ptLeptZ2_t, "ptLeptZ2_t/F" );
  tree_passedEvents->Branch( "ptLept3", &ptLept3_t, "ptLept3_t/F" );
  tree_passedEvents->Branch( "etaLeptZ1", &etaLeptZ1_t, "etaLeptZ1_t/F" );
  tree_passedEvents->Branch( "etaLeptZ2", &etaLeptZ2_t, "etaLeptZ2_t/F" );
  tree_passedEvents->Branch( "etaLept3", &etaLept3_t, "etaLept3_t/F" );
  tree_passedEvents->Branch( "ptJetB1", &ptJetB1_t, "ptJetB1_t/F" );
  tree_passedEvents->Branch( "ptJetB2", &ptJetB2_t, "ptJetB2_t/F" );
  tree_passedEvents->Branch( "ptJet3", &ptJet3_t, "ptJet3_t/F" );
  tree_passedEvents->Branch( "ptJet4", &ptJet4_t, "ptJet4_t/F" );
  tree_passedEvents->Branch( "etaJetB1", &etaJetB1_t, "etaJetB1_t/F" );
  tree_passedEvents->Branch( "etaJetB2", &etaJetB2_t, "etaJetB2_t/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1_t, "etaJet1_t/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2_t, "etaJet2_t/F" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "HLTSF", &HLTSF, "HLTSF/F" );
  tree_passedEvents->Branch( "PUWeight", &eventWeightPU, "eventWeightPU/F" );




ofstream ofs("run_event.txt");




  std::cout << std::endl << std::endl;
  std::cout << "+++ BEGINNING ANALYSIS LOOP" << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  if( isMC ) std::cout << "----> PU REWEIGHING: " << PUType_ << std::endl;
  std::cout << std::endl << std::endl;



  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 20000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    if( leptType_!="ALL" ) {
      if( leptType_=="ELE" && leptType==0 ) continue;
      if( leptType_=="MU" && leptType==1 ) continue;
    }




    h1_nvertex->Fill(nvertex, eventWeight);

    if( isMC ) {

//    // scale factor for double mu triggers:
//    if( leptType==0 ) {

//      float effDouble1_Run2011A = getMuonHLTSF_DoubleTrigger( ptLept1, etaLept1, "Run2011A" );
//      float effDouble2_Run2011A = getMuonHLTSF_DoubleTrigger( ptLept2, etaLept2, "Run2011A" );

//      float effDouble1_Run2011B = getMuonHLTSF_DoubleTrigger( ptLept1, etaLept1, "Run2011B" );
//      float effDouble2_Run2011B = getMuonHLTSF_DoubleTrigger( ptLept2, etaLept2, "Run2011B" );

//      float effSingle1_Run2011A1 = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011A1");
//      float effSingle2_Run2011A1 = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011A1");

//      float effSingle1_Run2011A2 = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011A2");
//      float effSingle2_Run2011A2 = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011A2");

//      float effSingle1_Run2011A3 = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011A3");
//      float effSingle2_Run2011A3 = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011A3");

//      float effSingle1_Run2011B = getMuonHLTSF_SingleTrigger( ptLept1, etaLept1, "Run2011B");
//      float effSingle2_Run2011B = getMuonHLTSF_SingleTrigger( ptLept2, etaLept2, "Run2011B");


//      float HLTSF_Run2011A1 = getEventHLTSF( effSingle1_Run2011A1, effSingle2_Run2011A1, effDouble1_Run2011A, effDouble2_Run2011A );
//      float HLTSF_Run2011A2 = getEventHLTSF( effSingle1_Run2011A2, effSingle2_Run2011A2, effDouble1_Run2011A, effDouble2_Run2011A );
//      float HLTSF_Run2011A3 = getEventHLTSF( effSingle1_Run2011A3, effSingle2_Run2011A3, effDouble1_Run2011A, effDouble2_Run2011A );
//      float HLTSF_Run2011B  = getEventHLTSF( effSingle1_Run2011B, effSingle2_Run2011B, effDouble1_Run2011B, effDouble2_Run2011B );


//      // weighted average over full run (weighted with lumi):
//      // LP11:
//      //HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 478.*HLTSF_Run2011A3)/(217.+920.+478.);
//      if( PUType_=="Run2011A" || PUType_=="Run2011A_73pb" )
//        HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 1000.*HLTSF_Run2011A3)/(217.+920.+1000.);
//      else if( PUType_=="HR11" ) 
//        HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 1000.*HLTSF_Run2011A3 + 2100.*HLTSF_Run2011B)/(217.+920.+1000.+2100.);
//      else if( PUType_=="HR11_v2" || PUType_=="HR11_73pb" )
//        HLTSF = (217.*HLTSF_Run2011A1 + 920.*HLTSF_Run2011A2 + 1000.*HLTSF_Run2011A3 + 2500.*HLTSF_Run2011B)/(217.+920.+1000.+2500.);

//      eventWeight *= HLTSF;

//    } else { //electrons

//      HLTSF = 1.;

//    }


      eventWeight *= fPUWeight->GetWeight(nPU);

    } // if is MC



    h1_nvertex_PUW->Fill(nvertex, eventWeight);


    if( !isMC ) { 

      // remove duplicate events:

      std::map<int, std::map<int, std::vector<int> > >::iterator it;

      it = run_lumi_ev_map.find(run);


      if( it==run_lumi_ev_map.end() ) {

        std::vector<int> events;
        events.push_back(event);
        std::map<int, std::vector<int> > lumi_ev_map;
        lumi_ev_map.insert( std::pair<int,std::vector<int> >(LS, events));
        run_lumi_ev_map.insert( std::pair<int, std::map<int, std::vector<int> > > (run, lumi_ev_map) );

      } else { //run exists, look for LS


        std::map<int, std::vector<int> >::iterator it_LS;
        it_LS = it->second.find( LS );

        if( it_LS==(it->second.end())  ) {

          std::vector<int> events;
          events.push_back(event);
          it->second.insert( std::pair<int, std::vector<int> > (LS, events) );

        } else { //LS exists, look for event

          std::vector<int>::iterator ev;
          for( ev=it_LS->second.begin(); ev!=it_LS->second.end(); ++ev )
            if( *ev==event ) break;


          if( ev==it_LS->second.end() ) {

            it_LS->second.push_back(event);

          } else {

            std::cout << "DISCARDING DUPLICATE EVENT!! Run: " << run << " LS: " << LS << " event: " << event << std::endl;

            continue;

          }
        }
      }


      h1_run->Fill( run, eventWeight );

    
    } //if is not mc



    // this is trilepton channel: require at least one other lepton:
    if( nLept<1 ) continue;


    nEvents_presel += eventWeight;


    h1_rhoPF_presel->Fill( rhoPF, eventWeight);


    h1_pfMet->Fill( pfMet, eventWeight );
    h1_metSignificance->Fill( metSignificance, eventWeight );



    TLorentzVector leptZ1, leptZ2;
    leptZ1.SetPtEtaPhiE( ptLeptZ1, etaLeptZ1, phiLeptZ1, eLeptZ1 );
    leptZ2.SetPtEtaPhiE( ptLeptZ2, etaLeptZ2, phiLeptZ2, eLeptZ2 );

    TLorentzVector otherLept;
    leptZ1.SetPtEtaPhiE( ptLept[0], etaLept[0], phiLept[0], eLept[0] );

    TLorentzVector diLepton = lept1+lept2;


    h1_ptLeptZ1->Fill( leptZ1.Pt(), eventWeight );
    h1_ptLeptZ2->Fill( leptZ2.Pt(), eventWeight );
    h1_etaLeptZ1->Fill( leptZ1.Eta(), eventWeight );
    h1_etaLeptZ2->Fill( leptZ2.Eta(), eventWeight );

    h1_deltaRllZ->Fill( leptZ2.DeltaR(&leptZ2), eventWeight );

    h1_ptZll->Fill( diLepton.Pt(), eventWeight );
    h1_etaZll->Fill( diLepton.Eta(), eventWeight );
    h1_mZll->Fill( diLepton.M(), eventWeight );



    TLorentzVector neutrino;
    neutrino.SetPtEtaPhiE( pfMet, 0., phiMet, pfMet );

    TLorentzVector W = otherLept + neutrino;

    h1_mTW->Fill( W.Mt() , eventWeight );

    TLorentzVector lZ = otherLept + diLepton;
    TLorentzVector lZ_plusMet = lZ + neutrino;

    float mT_lZmet = ( sqrt( lZ.Pt()*lZ.Pt() + lZ.M()*lZ.M() ) + pfMet )*( sqrt( lZ.Pt()*lZ.Pt() + lZ.M()*lZ.M() ) + pfMet )  -  lZ_plusMet.Pt()*lZ_plusMet.Pt();
    mT_lZmet = sqrt(mT_lZmet);

    h1_mT_lZmet->Fill( mT_lZmet, eventWeight );




    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "** LOG FOR RUN: " << run << "   EVENT: " << DEBUG_EVENTNUMBER << std::endl << std::endl;
      std::cout << "leptType: " << leptType << std::endl; 
      std::cout << "leptZ1.Pt(): " << leptZ1.Pt() << " leptZ1.Eta(): " << leptZ1.Eta() << std::endl;
      std::cout << "leptZ2.Pt(): " << leptZ2.Pt() << " leptZ2.Eta(): " << leptZ2.Eta() << std::endl;
      std::cout << "diLepton.M(): " << diLepton.M() << std::endl;
    }


    // ----------------------------
    // KINEMATIC SELECTION: LEPTONS
    // ----------------------------

    if( leptZ1.Pt() < ptLeptZ1_thresh_ ) continue;
    if( leptZ2.Pt() < ptLeptZ2_thresh_ ) continue;
    if( fabs(leptZ1.Eta()) > etaLeptZ1_thresh_ ) continue;
    if( fabs(leptZ2.Eta()) > etaLeptZ2_thresh_ ) continue;
    if( diLepton.M() < mZll_threshLo_ || diLepton.M() > mZll_threshHi_ ) continue;



    h1_nJets->Fill( nJets , eventWeight );

    if( nJets<4 ) continue;

    AnalysisJet jet1B, jet2B, jet3, jet4;

    // default: order by pt
    jet1B.SetPtEtaPhiE( ptJet[0], etaJet[0], phiJet[0], eJet[0]);
    jet2B.SetPtEtaPhiE( ptJet[1], etaJet[1], phiJet[1], eJet[1]);
    jet3.SetPtEtaPhiE( ptJet[2], etaJet[2], phiJet[2], eJet[2]);
    jet4.SetPtEtaPhiE( ptJet[3], etaJet[3], phiJet[3], eJet[3]);




    float bestBtag=-9999.;

    for( unsigned iJet=0; iJet<nJet; ++iJet) {

      AnalysisJet thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      thisJet.rmsCand = rmsCandJet[iJet];
      thisJet.ptD = ptDJet[iJet];
      thisJet.nCharged = nChargedJet[iJet];
      thisJet.nNeutral = nNeutralJet[iJet];
      thisJet.muonEnergyFraction = eMuonsJet[iJet]/thisJet.Energy();
      thisJet.electronEnergyFraction = eElectronsJet[iJet]/thisJet.Energy();

      thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet[iJet];
      thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet[iJet];
      thisJet.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet[iJet];
      thisJet.jetProbabilityBJetTag               = jetProbabilityBJetTagJet[iJet];

      thisJet.ptGen = ptJetGen[iJet];
      thisJet.etaGen = etaJetGen[iJet];
      thisJet.phiGen = phiJetGen[iJet];
      thisJet.eGen = eJetGen[iJet];

      //match to parton:
      int partFlavor=0;
      float deltaRmin=999.;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = thisJet.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin ) {
          partFlavor = pdgIdPart[iPart];
          deltaRmin = thisDeltaR;
        }
      }
      thisJet.pdgIdPart = partFlavor;


      float thisBtag;
      if( bTaggerType_=="TCHE" )
        thisBtag = thisJet.trackCountingHighEffBJetTag;
      else if( bTaggerType_=="SSVHE" ) 
        thisBtag = thisJet.simpleSecondaryVertexHighEffBJetTag;
;

      if( thisBtag > bestTag ) {
        bestTag = thisBtag;
        // slide them all:
        jet4 = jet3;
        jet3 = jet2B;
        jet2B = jet1B;
        jet1B = thisJet;
      }

    } // for jets



    if( bTaggerType_=="TCHE" ) {
      h1_bTagJet1B->Fill( jet1B.trackCountingHighEffBJetTag, eventWeight );
      h1_bTagJet2B->Fill( jet2B.trackCountingHighEffBJetTag, eventWeight );
    } else if( bTaggerType_=="SSVHE" ) {
      h1_bTagJet1B->Fill( jet1B.simpleSecondaryVertexHighEffBJetTag, eventWeight );
      h1_bTagJet2B->Fill( jet2B.simpleSecondaryVertexHighEffBJetTag, eventWeight );
    }

    h1_ssvheJet_best->Fill( jet_bestSSVHE.simpleSecondaryVertexHighEffBJetTag, eventWeight );
    h1_ssvheJet_2ndbest->Fill( jet_2ndbestSSVHE.simpleSecondaryVertexHighEffBJetTag, eventWeight );

    h1_deltaRbb_tche->Fill( jet_bestTCHE.DeltaR(&jet_2ndbestTCHE), eventWeight );
    h1_deltaRbb_ssvhe->Fill( jet_bestSSVHE.DeltaR(&jet_2ndbestSSVHE), eventWeight );


    




    // -------------------------
    // KINEMATIC SELECTION: JETS
    // -------------------------

    if( jet1.Pt() < ptJet1_thresh_ ) continue;
    if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet pt OK" << std::endl;
    if( jet2.Pt() < ptJet2_thresh_ ) continue;
    if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet pt OK" << std::endl;
    if( fabs(jet1.Eta()) > etaJet1_thresh_ ) continue;
    if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet eta OK" << std::endl;
    if( fabs(jet2.Eta()) > etaJet2_thresh_ ) continue;
    if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet eta OK" << std::endl;






    ptLeptZ1_t = ptLeptZ1;
    ptLeptZ2_t = ptLeptZ2;
    ptLept3_t = ptLept3;
    etaLeptZ1_t = etaLeptZ1;
    etaLeptZ2_t = etaLeptZ2;
    etaLept3_t = etaLept3;

    ptJetB1_t = jetB1.Pt();
    ptJetB2_t = jetB2.Pt();
    etaJetB1_t = jetB1.Pt();
    etaJetB2_t = jetB2.Eta();

    bTagJetB1_t = 

    ptJet3_t = jet3.Pt();
    ptJet4_t = jet4.Pt();
    etaJet3_t = jet3.Pt();
    etaJet4_t = jet4.Eta();


    isSignalRegion = (Zjj_nokinfit.M()>=75. && Zjj_nokinfit.M()<=105.);
    isSidebands = !isSignalRegion &&  Zjj_nokinfit.M()>60. && Zjj_nokinfit.M()<130.;

    // and fill tree:
    tree_passedEvents->Fill();



    // fill mZjj plots for all events:

    h1_mZjj->Fill( Zjj_nokinfit.M(), eventWeight);
 
    if( leptType==0 ) {
      h1_mZjj_MU->Fill( Zjj_nokinfit.M(), eventWeight);
      if( maxBTag_found>=0 ) h1_mZjj_nogluetag_MU->Fill( Zjj_nokinfit.M(), eventWeight);
    }
    if( leptType==1 ) {
      h1_mZjj_ELE->Fill( Zjj_nokinfit.M(), eventWeight);
      if( maxBTag_found>=0 ) h1_mZjj_nogluetag_ELE->Fill( Zjj_nokinfit.M(), eventWeight);
    }
    if( mZZ>150. && mZZ<250. ) h1_mZjj_loMass->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( mZZ>250. && mZZ<400. ) h1_mZjj_medMass->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( mZZ>400. ) h1_mZjj_hiMass->Fill( Zjj_nokinfit.M(), eventWeight);
    if( maxBTag_found==0 ) h1_mZjj_0btag->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( maxBTag_found==1 ) h1_mZjj_1btag->Fill( Zjj_nokinfit.M(), eventWeight);
    else if( maxBTag_found==2 ) h1_mZjj_2btag->Fill( Zjj_nokinfit.M(), eventWeight);

    h2_mZjj_vs_mZZ->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    if( maxBTag_found==0 ) h2_mZjj_vs_mZZ_0btag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==1 ) h2_mZjj_vs_mZZ_1btag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==2 ) h2_mZjj_vs_mZZ_2btag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==-1 ) h2_mZjj_vs_mZZ_gluetag->Fill( ZZ_nokinfit.M(), Zjj_nokinfit.M() );

    h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );
    if( maxBTag_found==0 )       h2_mZjj_vs_mZZ_kinfit_0btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==1 )  h2_mZjj_vs_mZZ_kinfit_1btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==2 )  h2_mZjj_vs_mZZ_kinfit_2btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
    else if( maxBTag_found==-1 ) h2_mZjj_vs_mZZ_kinfit_gluetag->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );


    //if( Zjj_nokinfit.M()<60. || Zjj_nokinfit.M()>130. ) continue;


    
    if( isSidebands ) {
 
      //fill sideband plots:

      h2_mZjj_vs_mZZ->Fill( mZZ, Zjj_nokinfit.M(), eventWeight );
      h2_mZjj_vs_mZZ_kinfit->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M(), eventWeight );

      h1_helicityLD_sidebands->Fill( helicityLD_selected, eventWeight );

      if( maxBTag_found==0 ) {
        h2_mZjj_vs_mZZ_0btag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_0btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_0btag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_0btag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_0btag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      } else if( maxBTag_found==1 ) {
        h2_mZjj_vs_mZZ_1btag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_1btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_1btag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_1btag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_1btag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      } else if( maxBTag_found==2 ) {
        h2_mZjj_vs_mZZ_2btag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_2btag->Fill(   ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_2btag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_2btag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_2btag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      } else if( maxBTag_found==-1 ) {
        h2_mZjj_vs_mZZ_gluetag->Fill( mZZ, Zjj_nokinfit.M() );
        h2_mZjj_vs_mZZ_kinfit_gluetag->Fill( ZZ_kinfit.M(), Zjj_nokinfit.M() );
        h1_mZZ_kinfit_hiMass_sidebands_gluetag->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==0 ) h1_mZZ_kinfit_hiMass_sidebands_gluetag_MU->Fill( ZZ_kinfit.M(), eventWeight );
        if( leptType==1 ) h1_mZZ_kinfit_hiMass_sidebands_gluetag_ELE->Fill( ZZ_kinfit.M(), eventWeight );
      }

      

    } else if( isSignalRegion ) {
    


      if( maxBTag_found==0 ) nEvents_presel_mZll_mZjj_0btag += eventWeight;
      if( maxBTag_found==1 ) nEvents_presel_mZll_mZjj_1btag += eventWeight;
      if( maxBTag_found==2 ) nEvents_presel_mZll_mZjj_2btag += eventWeight;
  
  
  
      if( helicityLD_selected < 0. ) 
        std::cout << "helicityLD_selected is less than 0!!! THIS IS NOT POSSIBLE!!" << std::endl;
  
  
        if( maxBTag_found>=0 ) ofs << run << " " << LS << " " << event << std::endl;
  
  
  
  
  
      //if( Zjj_nokinfit.M() < mZjj_threshLo_ || Zjj_nokinfit.M() > mZjj_threshHi_ ) continue;
  
    
  
      // percentages which define the cut and count windows:
      float mZZ_minPerc = 0.94;  //settled on -6/+10%
      float mZZ_maxPerc = 1.1;
  
      if( ZZ_kinfit.M() > 250.*mZZ_minPerc && ZZ_kinfit.M() < 250.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_250  += eventWeight;
          nEventsPassed_0btag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_250_MU  += eventWeight;
            nEventsPassed_0btag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_250_ELE  += eventWeight;
            nEventsPassed_0btag_250_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_250 += eventWeight;
          nEventsPassed_1btag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_250_MU  += eventWeight;
            nEventsPassed_1btag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_250_ELE  += eventWeight;
            nEventsPassed_1btag_250_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_250 += eventWeight;
          nEventsPassed_2btag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_250_MU  += eventWeight;
            nEventsPassed_2btag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_250_ELE  += eventWeight;
            nEventsPassed_2btag_250_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_250 += eventWeight;
          nEventsPassed_gluetag_250++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_250_MU  += eventWeight;
            nEventsPassed_gluetag_250_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_250_ELE  += eventWeight;
            nEventsPassed_gluetag_250_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 300.*mZZ_minPerc && ZZ_kinfit.M() < 300.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_300  += eventWeight;
          nEventsPassed_0btag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_300_MU  += eventWeight;
            nEventsPassed_0btag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_300_ELE  += eventWeight;
            nEventsPassed_0btag_300_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_300 += eventWeight;
          nEventsPassed_1btag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_300_MU  += eventWeight;
            nEventsPassed_1btag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_300_ELE  += eventWeight;
            nEventsPassed_1btag_300_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_300 += eventWeight;
          nEventsPassed_2btag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_300_MU  += eventWeight;
            nEventsPassed_2btag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_300_ELE  += eventWeight;
            nEventsPassed_2btag_300_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_300 += eventWeight;
          nEventsPassed_gluetag_300++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_300_MU  += eventWeight;
            nEventsPassed_gluetag_300_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_300_ELE  += eventWeight;
            nEventsPassed_gluetag_300_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 350.*mZZ_minPerc && ZZ_kinfit.M() < 350.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_350  += eventWeight;
          nEventsPassed_0btag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_350_MU  += eventWeight;
            nEventsPassed_0btag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_350_ELE  += eventWeight;
            nEventsPassed_0btag_350_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_350 += eventWeight;
          nEventsPassed_1btag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_350_MU  += eventWeight;
            nEventsPassed_1btag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_350_ELE  += eventWeight;
            nEventsPassed_1btag_350_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_350 += eventWeight;
          nEventsPassed_2btag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_350_MU  += eventWeight;
            nEventsPassed_2btag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_350_ELE  += eventWeight;
            nEventsPassed_2btag_350_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_350 += eventWeight;
          nEventsPassed_gluetag_350++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_350_MU  += eventWeight;
            nEventsPassed_gluetag_350_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_350_ELE  += eventWeight;
            nEventsPassed_gluetag_350_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 400.*mZZ_minPerc && ZZ_kinfit.M() < 400.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_400  += eventWeight;
          nEventsPassed_0btag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_400_MU  += eventWeight;
            nEventsPassed_0btag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_400_ELE  += eventWeight;
            nEventsPassed_0btag_400_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_400 += eventWeight;
          nEventsPassed_1btag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_400_MU  += eventWeight;
            nEventsPassed_1btag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_400_ELE  += eventWeight;
            nEventsPassed_1btag_400_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_400 += eventWeight;
          nEventsPassed_2btag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_400_MU  += eventWeight;
            nEventsPassed_2btag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_400_ELE  += eventWeight;
            nEventsPassed_2btag_400_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_400 += eventWeight;
          nEventsPassed_gluetag_400++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_400_MU  += eventWeight;
            nEventsPassed_gluetag_400_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_400_ELE  += eventWeight;
            nEventsPassed_gluetag_400_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 450.*mZZ_minPerc && ZZ_kinfit.M() < 450.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_450  += eventWeight;
          nEventsPassed_0btag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_450_MU  += eventWeight;
            nEventsPassed_0btag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_450_ELE  += eventWeight;
            nEventsPassed_0btag_450_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_450 += eventWeight;
          nEventsPassed_1btag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_450_MU  += eventWeight;
            nEventsPassed_1btag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_450_ELE  += eventWeight;
            nEventsPassed_1btag_450_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_450 += eventWeight;
          nEventsPassed_2btag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_450_MU  += eventWeight;
            nEventsPassed_2btag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_450_ELE  += eventWeight;
            nEventsPassed_2btag_450_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_450 += eventWeight;
          nEventsPassed_gluetag_450++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_450_MU  += eventWeight;
            nEventsPassed_gluetag_450_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_450_ELE  += eventWeight;
            nEventsPassed_gluetag_450_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 500.*mZZ_minPerc && ZZ_kinfit.M() < 500.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_500  += eventWeight;
          nEventsPassed_0btag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_500_MU  += eventWeight;
            nEventsPassed_0btag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_500_ELE  += eventWeight;
            nEventsPassed_0btag_500_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_500 += eventWeight;
          nEventsPassed_1btag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_500_MU  += eventWeight;
            nEventsPassed_1btag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_500_ELE  += eventWeight;
            nEventsPassed_1btag_500_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_500 += eventWeight;
          nEventsPassed_2btag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_500_MU  += eventWeight;
            nEventsPassed_2btag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_500_ELE  += eventWeight;
            nEventsPassed_2btag_500_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_500 += eventWeight;
          nEventsPassed_gluetag_500++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_500_MU  += eventWeight;
            nEventsPassed_gluetag_500_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_500_ELE  += eventWeight;
            nEventsPassed_gluetag_500_ELE++;
          }
        }
      } 
      if( ZZ_kinfit.M() > 600.*mZZ_minPerc && ZZ_kinfit.M() < 600.*mZZ_maxPerc ) {
        if( maxBTag_found==0 ) {
          nEventsPassed_fb_0btag_600  += eventWeight;
          nEventsPassed_0btag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_0btag_600_MU  += eventWeight;
            nEventsPassed_0btag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_0btag_600_ELE  += eventWeight;
            nEventsPassed_0btag_600_ELE++;
          }
        } else if( maxBTag_found==1 ) {
          nEventsPassed_fb_1btag_600 += eventWeight;
          nEventsPassed_1btag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_1btag_600_MU  += eventWeight;
            nEventsPassed_1btag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_1btag_600_ELE  += eventWeight;
            nEventsPassed_1btag_600_ELE++;
          }
        } else if( maxBTag_found==2 ) {
          nEventsPassed_fb_2btag_600 += eventWeight;
          nEventsPassed_2btag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_2btag_600_MU  += eventWeight;
            nEventsPassed_2btag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_2btag_600_ELE  += eventWeight;
            nEventsPassed_2btag_600_ELE++;
          }
        } else if( maxBTag_found==-1 ) {
          nEventsPassed_fb_gluetag_600 += eventWeight;
          nEventsPassed_gluetag_600++;
          if( leptType==0 ) {
            nEventsPassed_fb_gluetag_600_MU  += eventWeight;
            nEventsPassed_gluetag_600_MU++;
          } else if( leptType==1 ) {
            nEventsPassed_fb_gluetag_600_ELE  += eventWeight;
            nEventsPassed_gluetag_600_ELE++;
          }
        }
      } 
  
  
  
      // match to partons:
      TLorentzVector matchedPart1, matchedPart2;
      float bestDeltaRPart1=999.;
      float bestDeltaRPart2=999.;
      for( unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( abs(pdgIdPart[iPart])>6 ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        if( jet1_selected.DeltaR(thisPart) < bestDeltaRPart1 ) {
          bestDeltaRPart1 = jet1_selected.DeltaR(thisPart);
          matchedPart1 = thisPart;
        }
        if( jet2_selected.DeltaR(thisPart) < bestDeltaRPart2 ) {
          bestDeltaRPart2 = jet2_selected.DeltaR(thisPart);
          matchedPart2 = thisPart;
        }
      }
      float ptReso1_before = (isMC) ? ( jet1_nokinfit.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
      float ptReso2_before = (isMC) ? ( jet2_nokinfit.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
      h1_ptResoJet1_beforeKin->Fill( ptReso1_before, eventWeight );
      h1_ptResoJet2_beforeKin->Fill( ptReso2_before, eventWeight );
  
      float ptReso1_after = (isMC) ? ( jet1_selected.Pt()-matchedPart1.Pt() )/matchedPart1.Pt() : 0.;
      float ptReso2_after = (isMC) ? ( jet2_selected.Pt()-matchedPart2.Pt() )/matchedPart2.Pt() : 0.;
      h1_ptResoJet1_afterKin->Fill( ptReso1_after, eventWeight );
      h1_ptResoJet2_afterKin->Fill( ptReso2_after, eventWeight );
  
  
  
      TLorentzVector matchedZ;
      float bestDeltaRZ=999.;
      for( unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( pdgIdPart[iPart]!=23 ) continue;
        TLorentzVector thisZ;
        thisZ.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        if( Zjj_kinfit.DeltaR(thisZ) < bestDeltaRZ ) {
          bestDeltaRZ = Zjj_kinfit.DeltaR(thisZ);
          matchedZ = thisZ;
        }
      }
  
      bool eventIsMatched = bestDeltaRZ<0.2;
  
      float ptZreso_before = (isMC) ? (Zjj_nokinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
      float ptZreso_after  = (isMC) ? (Zjj_kinfit.Pt()-matchedZ.Pt())/matchedZ.Pt() : 0.;
      h1_ptZreso_beforeKin->Fill( ptZreso_before, eventWeight);
      h1_ptZreso_afterKin->Fill( ptZreso_after, eventWeight);
  
  
      TLorentzVector matchedH;
      bool foundH=false;
      for( unsigned iPart=0; iPart<nPart && !foundH; ++iPart ) {
        if( pdgIdPart[iPart]!=25 && pdgIdPart[iPart]!=39 ) continue;
        matchedH.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        foundH=true;
      }
  
      if( foundH ) {
        float mHreso_beforeKin = (ZZ_nokinfit.M()-matchedH.M())/matchedH.M();
        float mHreso_afterKin = (ZZ_kinfit.M()-matchedH.M())/matchedH.M();
        h1_mHreso_beforeKin->Fill(mHreso_beforeKin, eventWeight);
        h1_mHreso_afterKin->Fill(mHreso_afterKin, eventWeight);
      }
  
      //compare kinfit to Zmass constraint
      TLorentzVector Zjj_constr;
      Zjj_constr.SetXYZM( Zjj_nokinfit.Px(), Zjj_nokinfit.Py(), Zjj_nokinfit.Pz(), Zmass);
  
      TLorentzVector ZZ_constr = diLepton + Zjj_constr;
  
  
      float chiSquareProb = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());
      h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
      h1_kinfit_chiSquareProb->Fill( chiSquareProb, eventWeight ); 
  
  
  
  
  
  
      h1_pfMet->Fill( pfMet, eventWeight );
      h1_pfMetOverMZZ->Fill( pfMet/ZZ_kinfit.M(), eventWeight );
      h1_metSignificance->Fill( metSignificance, eventWeight );
      h1_mEtSig->Fill( mEtSig, eventWeight );
      if( maxBTag_found==2 ) {
        h1_pfMet_2btag->Fill( pfMet, eventWeight );
        h1_pfMetOverMZZ_2btag->Fill( pfMet/ZZ_kinfit.M(), eventWeight );
        h1_metSignificance_2btag->Fill( metSignificance, eventWeight );
        h1_mEtSig_2btag->Fill( mEtSig, eventWeight );
      }
  
      h1_rhoPF->Fill( rhoPF, eventWeight );
  
      h2_helicityLD_vs_mZZ->Fill( ZZ_kinfit.M(), helicityLD_selected, eventWeight );

      h1_tcheJet->Fill( jet1_selected.trackCountingHighEffBJetTag, eventWeight );
      h1_tcheJet->Fill( jet2_selected.trackCountingHighEffBJetTag, eventWeight );

      h1_ptJetRecoil->Fill( jetRecoil_selected.Pt(), eventWeight );
      h1_ptHiggs->Fill( ZZ_kinfit.Pt(), eventWeight );

      if( jetRecoil_selected.Pt()>0. && fabs(jetRecoil_selected.Eta())<2.4 ) {
        float QGLikelihoodJetRecoil = qglikeli->computeQGLikelihoodPU( jetRecoil_selected.Pt(), rhoPF, jetRecoil_selected.nCharged, jetRecoil_selected.nNeutral, jetRecoil_selected.ptD, -1. );
        h1_QGLikelihoodJetRecoil->Fill( QGLikelihoodJetRecoil, eventWeight);
        h1_QGLikelihoodProdRecoil->Fill( QGLikelihoodJetRecoil*jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight);
        if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
          h1_QGLikelihoodJetRecoil_MW400->Fill( QGLikelihoodJetRecoil, eventWeight);
          h1_QGLikelihoodProdRecoil_MW400->Fill( QGLikelihoodJetRecoil*jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight);
        }
      } else {
        h1_QGLikelihoodJetRecoil->Fill( -1., eventWeight);
        h1_QGLikelihoodProdRecoil->Fill( -1., eventWeight );
        if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
          h1_QGLikelihoodJetRecoil_MW400->Fill( -1., eventWeight);
          h1_QGLikelihoodProdRecoil_MW400->Fill( -1., eventWeight);
        }
      }
     
  
      if( jet1_selected.Pt()>jet2_selected.Pt() ) {
        h1_ptJet1->Fill( jet1_selected.Pt(), eventWeight );
        h1_ptJet2->Fill( jet2_selected.Pt(), eventWeight );
        h1_ptJet1_prekin->Fill( jet1_nokinfit.Pt(), eventWeight );
        h1_ptJet2_prekin->Fill( jet2_nokinfit.Pt(), eventWeight );
        h1_etaJet1->Fill( jet1_selected.Eta(), eventWeight );
        h1_etaJet2->Fill( jet2_selected.Eta(), eventWeight );
        h1_tcheJet1->Fill( jet1_selected.trackCountingHighEffBJetTag, eventWeight );
        h1_tcheJet2->Fill( jet2_selected.trackCountingHighEffBJetTag, eventWeight );
      } else {
        h1_ptJet1->Fill( jet2_selected.Pt(), eventWeight );
        h1_ptJet2->Fill( jet1_selected.Pt(), eventWeight );
        h1_ptJet1_prekin->Fill( jet2_nokinfit.Pt(), eventWeight );
        h1_ptJet2_prekin->Fill( jet1_nokinfit.Pt(), eventWeight );
        h1_etaJet1->Fill( jet2_selected.Eta(), eventWeight );
        h1_etaJet2->Fill( jet1_selected.Eta(), eventWeight );
        h1_tcheJet1->Fill( jet2_selected.trackCountingHighEffBJetTag, eventWeight );
        h1_tcheJet2->Fill( jet1_selected.trackCountingHighEffBJetTag, eventWeight );
      }
      h1_eMuonsJet1->Fill( jet1_selected.muonEnergyFraction, eventWeight );
      h1_eMuonsJet2->Fill( jet2_selected.muonEnergyFraction, eventWeight );
      h1_eElectronsJet1->Fill( jet1_selected.electronEnergyFraction, eventWeight );
      h1_eElectronsJet2->Fill( jet2_selected.electronEnergyFraction, eventWeight );
  
      h1_ptLept1->Fill( lept1.Pt(), eventWeight );
      h1_ptLept2->Fill( lept2.Pt(), eventWeight );
      h1_deltaRjj->Fill( jet1_selected.DeltaR(jet2_selected), eventWeight);
      h1_deltaRjj_prekin->Fill( jet1_nokinfit.DeltaR(jet2_nokinfit), eventWeight);
      h1_ptZll->Fill( diLepton.Pt(), eventWeight);
      h1_ptZjj->Fill( Zjj_kinfit.Pt(), eventWeight);
      if( leptType==0 )
        h1_mZmumu->Fill( diLepton.M(), eventWeight );
      else
        h1_mZee->Fill( diLepton.M(), eventWeight );
      h1_mZll->Fill( diLepton.M(), eventWeight);
  
  
      // fill QG plots only for the 0- and glue-tag category:
      if( maxBTag_found<=0 ) {
  
        h1_nChargedJet1->Fill(jet1_selected.nCharged, eventWeight);
        h1_nNeutralJet1->Fill(jet1_selected.nNeutral, eventWeight);
        h1_ptDJet1->Fill(jet1_selected.ptD, eventWeight);
      
        h1_nChargedJet2->Fill(jet2_selected.nCharged, eventWeight);
        h1_nNeutralJet2->Fill(jet2_selected.nNeutral, eventWeight);
        h1_ptDJet2->Fill(jet2_selected.ptD, eventWeight);
      
        h1_QGLikelihoodJet1->Fill( jet1_selected.QGLikelihood, eventWeight );
        h1_QGLikelihoodJet2->Fill( jet2_selected.QGLikelihood, eventWeight );
        h1_QGLikelihoodProd->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
  
        h1_QGLikelihoodNoPUJet1->Fill( jet1_selected.QGLikelihoodNoPU, eventWeight );
        h1_QGLikelihoodNoPUJet2->Fill( jet2_selected.QGLikelihoodNoPU, eventWeight );
        h1_QGLikelihoodNoPUProd->Fill( jet1_selected.QGLikelihoodNoPU*jet2_selected.QGLikelihoodNoPU, eventWeight );
  
        if( jet1_selected.Pt()>100. && jet1_selected.Pt()<123. ) h1_QGLikelihood_100_123->Fill( jet1_selected.QGLikelihood, eventWeight );
        if( jet2_selected.Pt()>100. && jet2_selected.Pt()<123. ) h1_QGLikelihood_100_123->Fill( jet2_selected.QGLikelihood, eventWeight );
        if( jet1_selected.Pt()>66. && jet1_selected.Pt()<81. ) h1_QGLikelihood_66_81->Fill( jet1_selected.QGLikelihood, eventWeight );
        if( jet2_selected.Pt()>66. && jet2_selected.Pt()<81. ) h1_QGLikelihood_66_81->Fill( jet2_selected.QGLikelihood, eventWeight );
  
  
        if( ZZ_kinfit.M()>0.94*300. && ZZ_kinfit.M()<1.1*300. ) {
          h1_QGLikelihoodJet1_MW300->Fill( jet1_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodJet2_MW300->Fill( jet2_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodProd_MW300->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        }
        if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
          h1_QGLikelihoodJet1_MW400->Fill( jet1_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodJet2_MW400->Fill( jet2_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodProd_MW400->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        }
        if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
          h1_QGLikelihoodJet1_MW500->Fill( jet1_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodJet2_MW500->Fill( jet2_selected.QGLikelihood, eventWeight );
          h1_QGLikelihoodProd_MW500->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
        }
  
  
        if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood < 0.1 ) {
          h1_mZZ_kinfit_hiMass_loQG->Fill(ZZ_kinfit.M(), eventWeight);
        } else {
          h1_mZZ_kinfit_hiMass_hiQG->Fill(ZZ_kinfit.M(), eventWeight);
        }
  
      } // if 0/glue tags
  
      h1_mZZ_nokinfit_hiMass_all->Fill( ZZ_nokinfit.M(), eventWeight);
      h1_mZZ_ZjjMassConstr_hiMass->Fill(ZZ_constr.M(), eventWeight);
      h1_mZZ_kinfit_hiMass_all->Fill( ZZ_kinfit.M(), eventWeight);
      if( nvertex>5 ) h1_mZZ_kinfit_hiMass_hiPU->Fill( ZZ_kinfit.M(), eventWeight);
      else h1_mZZ_kinfit_hiMass_loPU->Fill( ZZ_kinfit.M(), eventWeight);
      if( maxBTag_found>=0 ) h1_mZZ_kinfit_hiMass_nogluetag->Fill( ZZ_kinfit.M(), eventWeight);
      if( maxBTag_found==0 ) {
        h1_mZZ_kinfit_hiMass_0btag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_0btag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_0btag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      } else if( maxBTag_found==1 ) {
        h1_mZZ_kinfit_hiMass_1btag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_1btag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_1btag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      } else if( maxBTag_found==2 ) {
        h1_mZZ_kinfit_hiMass_2btag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_2btag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_2btag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      } else if( maxBTag_found==-1 ) {
        h1_mZZ_kinfit_hiMass_gluetag->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==0 )  h1_mZZ_kinfit_hiMass_gluetag_MU->Fill( ZZ_kinfit.M(), eventWeight);
        if( leptType==1 )  h1_mZZ_kinfit_hiMass_gluetag_ELE->Fill( ZZ_kinfit.M(), eventWeight);
      }
  
      h1_deltaRZmatching->Fill( bestDeltaRZ, eventWeight );
      if( maxBTag_found==0 && eventIsMatched ) h1_mZZ_kinfit_hiMass_0btag_matched->Fill( ZZ_kinfit.M(), eventWeight);
  
  
      h1_deltaRZZ->Fill(Zjj_nokinfit.DeltaR(diLepton), eventWeight);
  
      h1_ptZZ->Fill( ZZ_nokinfit.Pt(), eventWeight );
      h1_ptZZ_kinfit->Fill( ZZ_kinfit.Pt(), eventWeight );
      h1_etaZZ->Fill( ZZ_nokinfit.Eta(), eventWeight );
      h1_etaZZ_kinfit->Fill( ZZ_kinfit.Eta(), eventWeight );
  
      h1_helicityLD->Fill( helicityLD_selected, eventWeight );
      if( maxBTag_found>=0 ) h1_helicityLD_nogluetag->Fill( helicityLD_selected, eventWeight );
      h1_helicityLD_nokinfit->Fill( helicityLD_nokinfit_selected, eventWeight );
  
      h1_cosThetaStar->Fill(hangles_selected.helCosThetaStar, eventWeight);
      h1_cosTheta1->Fill(hangles_selected.helCosTheta1, eventWeight);
      h1_cosTheta2->Fill(hangles_selected.helCosTheta2, eventWeight);
      h1_phi->Fill(hangles_selected.helPhi, eventWeight);
      h1_phi1->Fill(hangles_selected.helPhi1, eventWeight);
  
  
  
      int partFlavor1=0;
      float deltaRmin1=999.;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( abs(pdgIdPart[iPart])>6 && pdgIdPart[iPart]!=21 ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet1_selected.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin1 ) {
          partFlavor1 = pdgIdPart[iPart];
          deltaRmin1 = thisDeltaR;
        }
      }
      h1_deltaR_part1->Fill(deltaRmin1, eventWeight);
      h1_partFlavorJet1->Fill( partFlavor1, eventWeight );
      if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) h1_partFlavorJet1_MW400->Fill( partFlavor1, eventWeight );
      if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) h1_partFlavorJet1_MW500->Fill( partFlavor1, eventWeight );
      if( deltaRmin1<0.5 ) h1_partFlavorJet1_matched->Fill( partFlavor1, eventWeight );
      else h1_partFlavorJet1_notmatched->Fill( partFlavor1, eventWeight );
  
      if( maxBTag_found==0 ) h1_partFlavorJet1_0btag->Fill( partFlavor1, eventWeight );
      else if( maxBTag_found==1 ) h1_partFlavorJet1_1btag->Fill( partFlavor1, eventWeight );
      else if( maxBTag_found==2 ) h1_partFlavorJet1_2btag->Fill( partFlavor1, eventWeight );
      else if( maxBTag_found==-1 ) h1_partFlavorJet1_gluetag->Fill( partFlavor1, eventWeight );
      jet1_selected.pdgIdPart = partFlavor1;
      if( partFlavor1==21 ) h1_QGLikelihoodJet1_gluonMatched->Fill( jet1_selected.QGLikelihood, eventWeight );
      else if( abs(partFlavor1)<5 ) h1_QGLikelihoodJet1_quarkMatched->Fill( jet1_selected.QGLikelihood, eventWeight );
  
      float deltaRmin2=999.;
      int partFlavor2=0;
      for(unsigned iPart=0; iPart<nPart; ++iPart ) {
        if( abs(pdgIdPart[iPart])>6 && pdgIdPart[iPart]!=21 ) continue;
        TLorentzVector thisPart;
        thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
        float thisDeltaR = jet2_selected.DeltaR(thisPart);
        if( thisDeltaR<deltaRmin2 ) {
          partFlavor2 = pdgIdPart[iPart];
          deltaRmin2 = thisDeltaR;
        }
      }
      h1_deltaR_part2->Fill(deltaRmin2, eventWeight);
      h1_partFlavorJet2->Fill( partFlavor2, eventWeight );
      if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) h1_partFlavorJet2_MW400->Fill( partFlavor2, eventWeight );
      if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) h1_partFlavorJet2_MW500->Fill( partFlavor2, eventWeight );
      if( deltaRmin2<0.5 ) h1_partFlavorJet2_matched->Fill( partFlavor2, eventWeight );
      else h1_partFlavorJet2_notmatched->Fill( partFlavor2, eventWeight );
  
      if( maxBTag_found==0 ) h1_partFlavorJet2_0btag->Fill( partFlavor2, eventWeight );
      else if( maxBTag_found==1 ) h1_partFlavorJet2_1btag->Fill( partFlavor2, eventWeight );
      else if( maxBTag_found==2 ) h1_partFlavorJet2_2btag->Fill( partFlavor2, eventWeight );
      else if( maxBTag_found==-1 ) h1_partFlavorJet2_gluetag->Fill( partFlavor2, eventWeight );
      jet2_selected.pdgIdPart = partFlavor2;
      if( partFlavor2==21 ) h1_QGLikelihoodJet2_gluonMatched->Fill( jet2_selected.QGLikelihood, eventWeight );
      else if( abs(partFlavor2)<5 ) h1_QGLikelihoodJet2_quarkMatched->Fill( jet2_selected.QGLikelihood, eventWeight );
  
  
      if( maxBTag_found<=0 ) {

        if( partFlavor1==21 || partFlavor2==21 ) {
          h1_QGLikelihoodProd_oneGluon->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          h1_mZZ_kinfit_hiMass_oneGluon->Fill( mZZ, eventWeight );
        }
        if( partFlavor1==21 && partFlavor2==21 ) {
          h1_QGLikelihoodProd_twoGluon->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          h1_mZZ_kinfit_hiMass_twoGluon->Fill( mZZ, eventWeight );
        }
        if( abs(partFlavor1)<5 && abs(partFlavor2)<5 ) {
          h1_QGLikelihoodProd_twoQuark->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          h1_mZZ_kinfit_hiMass_twoQuark->Fill( mZZ, eventWeight );
          if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
            h1_QGLikelihoodProd_twoQuark_MW400->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          }
          if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
            h1_QGLikelihoodProd_twoQuark_MW500->Fill( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood, eventWeight );
          }
          if( maxBTag_found<=0 ) {
            h1_nEvents_partFlavor_beforeQG->Fill( 0., eventWeight );
            if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEvents_partFlavor_afterQG->Fill( 0., eventWeight );
            if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
              h1_nEventsMW400_partFlavor_beforeQG->Fill( 0., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW400_partFlavor_afterQG->Fill( 0., eventWeight );
            }
            if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
              h1_nEventsMW500_partFlavor_beforeQG->Fill( 0., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW500_partFlavor_afterQG->Fill( 0., eventWeight );
            }
          }
        } else { 
          if( maxBTag_found<=0 ) {
            h1_nEvents_partFlavor_beforeQG->Fill( 1., eventWeight );
            if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEvents_partFlavor_afterQG->Fill( 1., eventWeight );
            if( ZZ_kinfit.M()>0.94*400. && ZZ_kinfit.M()<1.1*400. ) {
              h1_nEventsMW400_partFlavor_beforeQG->Fill( 1., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW400_partFlavor_afterQG->Fill( 1., eventWeight );
            }
            if( ZZ_kinfit.M()>0.94*500. && ZZ_kinfit.M()<1.1*500. ) {
              h1_nEventsMW500_partFlavor_beforeQG->Fill( 1., eventWeight );
              if( jet1_selected.QGLikelihood*jet2_selected.QGLikelihood>0.1 ) h1_nEventsMW500_partFlavor_afterQG->Fill( 1., eventWeight );
            }
          }
        }

      } //if 0/glue tag
  

    } //if signal region
  
  } //for entries


  std::cout << "nEvents_presel: " << nEvents_presel << std::endl;
  std::cout << "nEvents_presel_mZll: " << nEvents_presel_mZll << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj: " << nEvents_presel_mZll_mZjj << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj_0btag: " << nEvents_presel_mZll_mZjj_0btag << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj_1btag: " << nEvents_presel_mZll_mZjj_1btag << std::endl;
  std::cout << "nEvents_presel_mZll_mZjj_2btag: " << nEvents_presel_mZll_mZjj_2btag << std::endl;

  h1_nCounter->SetBinContent(1, nCounter_);
  h1_nCounterW->SetBinContent(1, nCounterW_);
  h1_nCounterPU->SetBinContent(1, nCounterPU_);

  float eff_gluetag_250 = nEventsPassed_fb_gluetag_250/nCounterW_;
  float eff_0btag_250 = nEventsPassed_fb_0btag_250/nCounterW_;
  float eff_1btag_250 = nEventsPassed_fb_1btag_250/nCounterW_;
  float eff_2btag_250 = nEventsPassed_fb_2btag_250/nCounterW_;

  float eff_gluetag_250_ELE = nEventsPassed_fb_gluetag_250_ELE/nCounterW_;
  float eff_0btag_250_ELE = nEventsPassed_fb_0btag_250_ELE/nCounterW_;
  float eff_1btag_250_ELE = nEventsPassed_fb_1btag_250_ELE/nCounterW_;
  float eff_2btag_250_ELE = nEventsPassed_fb_2btag_250_ELE/nCounterW_;

  float eff_gluetag_250_MU = nEventsPassed_fb_gluetag_250_MU/nCounterW_;
  float eff_0btag_250_MU = nEventsPassed_fb_0btag_250_MU/nCounterW_;
  float eff_1btag_250_MU = nEventsPassed_fb_1btag_250_MU/nCounterW_;
  float eff_2btag_250_MU = nEventsPassed_fb_2btag_250_MU/nCounterW_;

  float eff_gluetag_300 = nEventsPassed_fb_gluetag_300/nCounterW_;
  float eff_0btag_300 = nEventsPassed_fb_0btag_300/nCounterW_;
  float eff_1btag_300 = nEventsPassed_fb_1btag_300/nCounterW_;
  float eff_2btag_300 = nEventsPassed_fb_2btag_300/nCounterW_;

  float eff_gluetag_300_ELE = nEventsPassed_fb_gluetag_300_ELE/nCounterW_;
  float eff_0btag_300_ELE = nEventsPassed_fb_0btag_300_ELE/nCounterW_;
  float eff_1btag_300_ELE = nEventsPassed_fb_1btag_300_ELE/nCounterW_;
  float eff_2btag_300_ELE = nEventsPassed_fb_2btag_300_ELE/nCounterW_;

  float eff_gluetag_300_MU = nEventsPassed_fb_gluetag_300_MU/nCounterW_;
  float eff_0btag_300_MU = nEventsPassed_fb_0btag_300_MU/nCounterW_;
  float eff_1btag_300_MU = nEventsPassed_fb_1btag_300_MU/nCounterW_;
  float eff_2btag_300_MU = nEventsPassed_fb_2btag_300_MU/nCounterW_;

  float eff_gluetag_350 = nEventsPassed_fb_gluetag_350/nCounterW_;
  float eff_0btag_350 = nEventsPassed_fb_0btag_350/nCounterW_;
  float eff_1btag_350 = nEventsPassed_fb_1btag_350/nCounterW_;
  float eff_2btag_350 = nEventsPassed_fb_2btag_350/nCounterW_;

  float eff_gluetag_350_ELE = nEventsPassed_fb_gluetag_350_ELE/nCounterW_;
  float eff_0btag_350_ELE = nEventsPassed_fb_0btag_350_ELE/nCounterW_;
  float eff_1btag_350_ELE = nEventsPassed_fb_1btag_350_ELE/nCounterW_;
  float eff_2btag_350_ELE = nEventsPassed_fb_2btag_350_ELE/nCounterW_;

  float eff_gluetag_350_MU = nEventsPassed_fb_gluetag_350_MU/nCounterW_;
  float eff_0btag_350_MU = nEventsPassed_fb_0btag_350_MU/nCounterW_;
  float eff_1btag_350_MU = nEventsPassed_fb_1btag_350_MU/nCounterW_;
  float eff_2btag_350_MU = nEventsPassed_fb_2btag_350_MU/nCounterW_;

  float eff_gluetag_400 = nEventsPassed_fb_gluetag_400/nCounterW_;
  float eff_0btag_400 = nEventsPassed_fb_0btag_400/nCounterW_;
  float eff_1btag_400 = nEventsPassed_fb_1btag_400/nCounterW_;
  float eff_2btag_400 = nEventsPassed_fb_2btag_400/nCounterW_;

  float eff_gluetag_400_ELE = nEventsPassed_fb_gluetag_400_ELE/nCounterW_;
  float eff_0btag_400_ELE = nEventsPassed_fb_0btag_400_ELE/nCounterW_;
  float eff_1btag_400_ELE = nEventsPassed_fb_1btag_400_ELE/nCounterW_;
  float eff_2btag_400_ELE = nEventsPassed_fb_2btag_400_ELE/nCounterW_;

  float eff_gluetag_400_MU = nEventsPassed_fb_gluetag_400_MU/nCounterW_;
  float eff_0btag_400_MU = nEventsPassed_fb_0btag_400_MU/nCounterW_;
  float eff_1btag_400_MU = nEventsPassed_fb_1btag_400_MU/nCounterW_;
  float eff_2btag_400_MU = nEventsPassed_fb_2btag_400_MU/nCounterW_;

  float eff_gluetag_450 = nEventsPassed_fb_gluetag_450/nCounterW_;
  float eff_0btag_450 = nEventsPassed_fb_0btag_450/nCounterW_;
  float eff_1btag_450 = nEventsPassed_fb_1btag_450/nCounterW_;
  float eff_2btag_450 = nEventsPassed_fb_2btag_450/nCounterW_;

  float eff_gluetag_450_ELE = nEventsPassed_fb_gluetag_450_ELE/nCounterW_;
  float eff_0btag_450_ELE = nEventsPassed_fb_0btag_450_ELE/nCounterW_;
  float eff_1btag_450_ELE = nEventsPassed_fb_1btag_450_ELE/nCounterW_;
  float eff_2btag_450_ELE = nEventsPassed_fb_2btag_450_ELE/nCounterW_;

  float eff_gluetag_450_MU = nEventsPassed_fb_gluetag_450_MU/nCounterW_;
  float eff_0btag_450_MU = nEventsPassed_fb_0btag_450_MU/nCounterW_;
  float eff_1btag_450_MU = nEventsPassed_fb_1btag_450_MU/nCounterW_;
  float eff_2btag_450_MU = nEventsPassed_fb_2btag_450_MU/nCounterW_;

  float eff_gluetag_500 = nEventsPassed_fb_gluetag_500/nCounterW_;
  float eff_0btag_500 = nEventsPassed_fb_0btag_500/nCounterW_;
  float eff_1btag_500 = nEventsPassed_fb_1btag_500/nCounterW_;
  float eff_2btag_500 = nEventsPassed_fb_2btag_500/nCounterW_;

  float eff_gluetag_500_ELE = nEventsPassed_fb_gluetag_500_ELE/nCounterW_;
  float eff_0btag_500_ELE = nEventsPassed_fb_0btag_500_ELE/nCounterW_;
  float eff_1btag_500_ELE = nEventsPassed_fb_1btag_500_ELE/nCounterW_;
  float eff_2btag_500_ELE = nEventsPassed_fb_2btag_500_ELE/nCounterW_;

  float eff_gluetag_500_MU = nEventsPassed_fb_gluetag_500_MU/nCounterW_;
  float eff_0btag_500_MU = nEventsPassed_fb_0btag_500_MU/nCounterW_;
  float eff_1btag_500_MU = nEventsPassed_fb_1btag_500_MU/nCounterW_;
  float eff_2btag_500_MU = nEventsPassed_fb_2btag_500_MU/nCounterW_;

  float eff_gluetag_600 = nEventsPassed_fb_gluetag_600/nCounterW_;
  float eff_0btag_600 = nEventsPassed_fb_0btag_600/nCounterW_;
  float eff_1btag_600 = nEventsPassed_fb_1btag_600/nCounterW_;
  float eff_2btag_600 = nEventsPassed_fb_2btag_600/nCounterW_;

  float eff_gluetag_600_ELE = nEventsPassed_fb_gluetag_600_ELE/nCounterW_;
  float eff_0btag_600_ELE = nEventsPassed_fb_0btag_600_ELE/nCounterW_;
  float eff_1btag_600_ELE = nEventsPassed_fb_1btag_600_ELE/nCounterW_;
  float eff_2btag_600_ELE = nEventsPassed_fb_2btag_600_ELE/nCounterW_;

  float eff_gluetag_600_MU = nEventsPassed_fb_gluetag_600_MU/nCounterW_;
  float eff_0btag_600_MU = nEventsPassed_fb_0btag_600_MU/nCounterW_;
  float eff_1btag_600_MU = nEventsPassed_fb_1btag_600_MU/nCounterW_;
  float eff_2btag_600_MU = nEventsPassed_fb_2btag_600_MU/nCounterW_;


  h1_nEvents_fb_gluetag_250->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_250);
  h1_nEvents_fb_0btag_250->SetBinContent(1,1000.*nEventsPassed_fb_0btag_250);
  h1_nEvents_fb_1btag_250->SetBinContent(1,1000.*nEventsPassed_fb_1btag_250);
  h1_nEvents_fb_2btag_250->SetBinContent(1,1000.*nEventsPassed_fb_2btag_250);

  h1_nEvents_fb_gluetag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_250_ELE);
  h1_nEvents_fb_0btag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_250_ELE);
  h1_nEvents_fb_1btag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_250_ELE);
  h1_nEvents_fb_2btag_250_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_250_ELE);

  h1_nEvents_fb_gluetag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_250_MU);
  h1_nEvents_fb_0btag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_250_MU);
  h1_nEvents_fb_1btag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_250_MU);
  h1_nEvents_fb_2btag_250_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_250_MU);

  h1_eff_gluetag_250->SetBinContent(1,eff_gluetag_250);
  h1_eff_0btag_250->SetBinContent(1,eff_0btag_250);
  h1_eff_1btag_250->SetBinContent(1,eff_1btag_250);
  h1_eff_2btag_250->SetBinContent(1,eff_2btag_250);

  h1_eff_gluetag_250_ELE->SetBinContent(1,eff_gluetag_250_ELE);
  h1_eff_0btag_250_ELE->SetBinContent(1,eff_0btag_250_ELE);
  h1_eff_1btag_250_ELE->SetBinContent(1,eff_1btag_250_ELE);
  h1_eff_2btag_250_ELE->SetBinContent(1,eff_2btag_250_ELE);

  h1_eff_gluetag_250_MU->SetBinContent(1,eff_gluetag_250_MU);
  h1_eff_0btag_250_MU->SetBinContent(1,eff_0btag_250_MU);
  h1_eff_1btag_250_MU->SetBinContent(1,eff_1btag_250_MU);
  h1_eff_2btag_250_MU->SetBinContent(1,eff_2btag_250_MU);


  h1_nEvents_fb_gluetag_300->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_300);
  h1_nEvents_fb_0btag_300->SetBinContent(1,1000.*nEventsPassed_fb_0btag_300);
  h1_nEvents_fb_1btag_300->SetBinContent(1,1000.*nEventsPassed_fb_1btag_300);
  h1_nEvents_fb_2btag_300->SetBinContent(1,1000.*nEventsPassed_fb_2btag_300);

  h1_nEvents_fb_gluetag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_300_ELE);
  h1_nEvents_fb_0btag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_300_ELE);
  h1_nEvents_fb_1btag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_300_ELE);
  h1_nEvents_fb_2btag_300_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_300_ELE);

  h1_nEvents_fb_gluetag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_300_MU);
  h1_nEvents_fb_0btag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_300_MU);
  h1_nEvents_fb_1btag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_300_MU);
  h1_nEvents_fb_2btag_300_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_300_MU);

  h1_eff_gluetag_300->SetBinContent(1,eff_gluetag_300);
  h1_eff_0btag_300->SetBinContent(1,eff_0btag_300);
  h1_eff_1btag_300->SetBinContent(1,eff_1btag_300);
  h1_eff_2btag_300->SetBinContent(1,eff_2btag_300);

  h1_eff_gluetag_300_ELE->SetBinContent(1,eff_gluetag_300_ELE);
  h1_eff_0btag_300_ELE->SetBinContent(1,eff_0btag_300_ELE);
  h1_eff_1btag_300_ELE->SetBinContent(1,eff_1btag_300_ELE);
  h1_eff_2btag_300_ELE->SetBinContent(1,eff_2btag_300_ELE);

  h1_eff_gluetag_300_MU->SetBinContent(1,eff_gluetag_300_MU);
  h1_eff_0btag_300_MU->SetBinContent(1,eff_0btag_300_MU);
  h1_eff_1btag_300_MU->SetBinContent(1,eff_1btag_300_MU);
  h1_eff_2btag_300_MU->SetBinContent(1,eff_2btag_300_MU);


  h1_nEvents_fb_gluetag_350->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_350);
  h1_nEvents_fb_0btag_350->SetBinContent(1,1000.*nEventsPassed_fb_0btag_350);
  h1_nEvents_fb_1btag_350->SetBinContent(1,1000.*nEventsPassed_fb_1btag_350);
  h1_nEvents_fb_2btag_350->SetBinContent(1,1000.*nEventsPassed_fb_2btag_350);

  h1_nEvents_fb_gluetag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_350_ELE);
  h1_nEvents_fb_0btag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_350_ELE);
  h1_nEvents_fb_1btag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_350_ELE);
  h1_nEvents_fb_2btag_350_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_350_ELE);

  h1_nEvents_fb_gluetag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_350_MU);
  h1_nEvents_fb_0btag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_350_MU);
  h1_nEvents_fb_1btag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_350_MU);
  h1_nEvents_fb_2btag_350_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_350_MU);

  h1_eff_gluetag_350->SetBinContent(1,eff_gluetag_350);
  h1_eff_0btag_350->SetBinContent(1,eff_0btag_350);
  h1_eff_1btag_350->SetBinContent(1,eff_1btag_350);
  h1_eff_2btag_350->SetBinContent(1,eff_2btag_350);

  h1_eff_gluetag_350_ELE->SetBinContent(1,eff_gluetag_350_ELE);
  h1_eff_0btag_350_ELE->SetBinContent(1,eff_0btag_350_ELE);
  h1_eff_1btag_350_ELE->SetBinContent(1,eff_1btag_350_ELE);
  h1_eff_2btag_350_ELE->SetBinContent(1,eff_2btag_350_ELE);

  h1_eff_gluetag_350_MU->SetBinContent(1,eff_gluetag_350_MU);
  h1_eff_0btag_350_MU->SetBinContent(1,eff_0btag_350_MU);
  h1_eff_1btag_350_MU->SetBinContent(1,eff_1btag_350_MU);
  h1_eff_2btag_350_MU->SetBinContent(1,eff_2btag_350_MU);


  h1_nEvents_fb_gluetag_400->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_400);
  h1_nEvents_fb_0btag_400->SetBinContent(1,1000.*nEventsPassed_fb_0btag_400);
  h1_nEvents_fb_1btag_400->SetBinContent(1,1000.*nEventsPassed_fb_1btag_400);
  h1_nEvents_fb_2btag_400->SetBinContent(1,1000.*nEventsPassed_fb_2btag_400);

  h1_nEvents_fb_gluetag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_400_ELE);
  h1_nEvents_fb_0btag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_400_ELE);
  h1_nEvents_fb_1btag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_400_ELE);
  h1_nEvents_fb_2btag_400_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_400_ELE);

  h1_nEvents_fb_gluetag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_400_MU);
  h1_nEvents_fb_0btag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_400_MU);
  h1_nEvents_fb_1btag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_400_MU);
  h1_nEvents_fb_2btag_400_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_400_MU);

  h1_eff_gluetag_400->SetBinContent(1,eff_gluetag_400);
  h1_eff_0btag_400->SetBinContent(1,eff_0btag_400);
  h1_eff_1btag_400->SetBinContent(1,eff_1btag_400);
  h1_eff_2btag_400->SetBinContent(1,eff_2btag_400);

  h1_eff_gluetag_400_ELE->SetBinContent(1,eff_gluetag_400_ELE);
  h1_eff_0btag_400_ELE->SetBinContent(1,eff_0btag_400_ELE);
  h1_eff_1btag_400_ELE->SetBinContent(1,eff_1btag_400_ELE);
  h1_eff_2btag_400_ELE->SetBinContent(1,eff_2btag_400_ELE);

  h1_eff_gluetag_400_MU->SetBinContent(1,eff_gluetag_400_MU);
  h1_eff_0btag_400_MU->SetBinContent(1,eff_0btag_400_MU);
  h1_eff_1btag_400_MU->SetBinContent(1,eff_1btag_400_MU);
  h1_eff_2btag_400_MU->SetBinContent(1,eff_2btag_400_MU);


  h1_nEvents_fb_gluetag_450->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_450);
  h1_nEvents_fb_0btag_450->SetBinContent(1,1000.*nEventsPassed_fb_0btag_450);
  h1_nEvents_fb_1btag_450->SetBinContent(1,1000.*nEventsPassed_fb_1btag_450);
  h1_nEvents_fb_2btag_450->SetBinContent(1,1000.*nEventsPassed_fb_2btag_450);

  h1_nEvents_fb_gluetag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_450_ELE);
  h1_nEvents_fb_0btag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_450_ELE);
  h1_nEvents_fb_1btag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_450_ELE);
  h1_nEvents_fb_2btag_450_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_450_ELE);

  h1_nEvents_fb_gluetag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_450_MU);
  h1_nEvents_fb_0btag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_450_MU);
  h1_nEvents_fb_1btag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_450_MU);
  h1_nEvents_fb_2btag_450_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_450_MU);

  h1_eff_gluetag_450->SetBinContent(1,eff_gluetag_450);
  h1_eff_0btag_450->SetBinContent(1,eff_0btag_450);
  h1_eff_1btag_450->SetBinContent(1,eff_1btag_450);
  h1_eff_2btag_450->SetBinContent(1,eff_2btag_450);

  h1_eff_gluetag_450_ELE->SetBinContent(1,eff_gluetag_450_ELE);
  h1_eff_0btag_450_ELE->SetBinContent(1,eff_0btag_450_ELE);
  h1_eff_1btag_450_ELE->SetBinContent(1,eff_1btag_450_ELE);
  h1_eff_2btag_450_ELE->SetBinContent(1,eff_2btag_450_ELE);

  h1_eff_gluetag_450_MU->SetBinContent(1,eff_gluetag_450_MU);
  h1_eff_0btag_450_MU->SetBinContent(1,eff_0btag_450_MU);
  h1_eff_1btag_450_MU->SetBinContent(1,eff_1btag_450_MU);
  h1_eff_2btag_450_MU->SetBinContent(1,eff_2btag_450_MU);


  h1_nEvents_fb_gluetag_500->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_500);
  h1_nEvents_fb_0btag_500->SetBinContent(1,1000.*nEventsPassed_fb_0btag_500);
  h1_nEvents_fb_1btag_500->SetBinContent(1,1000.*nEventsPassed_fb_1btag_500);
  h1_nEvents_fb_2btag_500->SetBinContent(1,1000.*nEventsPassed_fb_2btag_500);

  h1_nEvents_fb_gluetag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_500_ELE);
  h1_nEvents_fb_0btag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_500_ELE);
  h1_nEvents_fb_1btag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_500_ELE);
  h1_nEvents_fb_2btag_500_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_500_ELE);

  h1_nEvents_fb_gluetag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_500_MU);
  h1_nEvents_fb_0btag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_500_MU);
  h1_nEvents_fb_1btag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_500_MU);
  h1_nEvents_fb_2btag_500_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_500_MU);

  h1_eff_gluetag_500->SetBinContent(1,eff_gluetag_500);
  h1_eff_0btag_500->SetBinContent(1,eff_0btag_500);
  h1_eff_1btag_500->SetBinContent(1,eff_1btag_500);
  h1_eff_2btag_500->SetBinContent(1,eff_2btag_500);

  h1_eff_gluetag_500_ELE->SetBinContent(1,eff_gluetag_500_ELE);
  h1_eff_0btag_500_ELE->SetBinContent(1,eff_0btag_500_ELE);
  h1_eff_1btag_500_ELE->SetBinContent(1,eff_1btag_500_ELE);
  h1_eff_2btag_500_ELE->SetBinContent(1,eff_2btag_500_ELE);

  h1_eff_gluetag_500_MU->SetBinContent(1,eff_gluetag_500_MU);
  h1_eff_0btag_500_MU->SetBinContent(1,eff_0btag_500_MU);
  h1_eff_1btag_500_MU->SetBinContent(1,eff_1btag_500_MU);
  h1_eff_2btag_500_MU->SetBinContent(1,eff_2btag_500_MU);


  h1_nEvents_fb_gluetag_600->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_600);
  h1_nEvents_fb_0btag_600->SetBinContent(1,1000.*nEventsPassed_fb_0btag_600);
  h1_nEvents_fb_1btag_600->SetBinContent(1,1000.*nEventsPassed_fb_1btag_600);
  h1_nEvents_fb_2btag_600->SetBinContent(1,1000.*nEventsPassed_fb_2btag_600);

  h1_nEvents_fb_gluetag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_600_ELE);
  h1_nEvents_fb_0btag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_0btag_600_ELE);
  h1_nEvents_fb_1btag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_1btag_600_ELE);
  h1_nEvents_fb_2btag_600_ELE->SetBinContent(1,1000.*nEventsPassed_fb_2btag_600_ELE);

  h1_nEvents_fb_gluetag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_gluetag_600_MU);
  h1_nEvents_fb_0btag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_0btag_600_MU);
  h1_nEvents_fb_1btag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_1btag_600_MU);
  h1_nEvents_fb_2btag_600_MU->SetBinContent(1,1000.*nEventsPassed_fb_2btag_600_MU);

  h1_eff_gluetag_600->SetBinContent(1,eff_gluetag_600);
  h1_eff_0btag_600->SetBinContent(1,eff_0btag_600);
  h1_eff_1btag_600->SetBinContent(1,eff_1btag_600);
  h1_eff_2btag_600->SetBinContent(1,eff_2btag_600);

  h1_eff_gluetag_600_ELE->SetBinContent(1,eff_gluetag_600_ELE);
  h1_eff_0btag_600_ELE->SetBinContent(1,eff_0btag_600_ELE);
  h1_eff_1btag_600_ELE->SetBinContent(1,eff_1btag_600_ELE);
  h1_eff_2btag_600_ELE->SetBinContent(1,eff_2btag_600_ELE);

  h1_eff_gluetag_600_MU->SetBinContent(1,eff_gluetag_600_MU);
  h1_eff_0btag_600_MU->SetBinContent(1,eff_0btag_600_MU);
  h1_eff_1btag_600_MU->SetBinContent(1,eff_1btag_600_MU);
  h1_eff_2btag_600_MU->SetBinContent(1,eff_2btag_600_MU);



  std::cout << std::endl << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  if( isMC ) std::cout << "----> PU REWEIGHING: " << PUType_ << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << "----> 250 GeV (235-275): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_250 << " ev/fb-1  (" << nEventsPassed_0btag_250 << " events)" << " Efficiency: " << 100.*eff_0btag_250 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_250 << " ev/fb-1  (" << nEventsPassed_1btag_250 << " events)" << " Efficiency: " << 100.*eff_1btag_250 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_250 << " ev/fb-1  (" << nEventsPassed_2btag_250 << " events)" << " Efficiency: " << 100.*eff_2btag_250 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 300 GeV (282-330): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_300 << " ev/fb-1  (" << nEventsPassed_0btag_300 << " events)" << " Efficiency: " << 100.*eff_0btag_300 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_300 << " ev/fb-1  (" << nEventsPassed_1btag_300 << " events)" << " Efficiency: " << 100.*eff_1btag_300 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_300 << " ev/fb-1  (" << nEventsPassed_2btag_300 << " events)" << " Efficiency: " << 100.*eff_2btag_300 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 350 GeV (329-385): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_350 << " ev/fb-1  (" << nEventsPassed_0btag_350 << " events)" << " Efficiency: " << 100.*eff_0btag_350 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_350 << " ev/fb-1  (" << nEventsPassed_1btag_350 << " events)" << " Efficiency: " << 100.*eff_1btag_350 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_350 << " ev/fb-1  (" << nEventsPassed_2btag_350 << " events)" << " Efficiency: " << 100.*eff_2btag_350 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 400 GeV (376-440): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_400 << " ev/fb-1  (" << nEventsPassed_0btag_400 << " events)" << " Efficiency: " << 100.*eff_0btag_400 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_400 << " ev/fb-1  (" << nEventsPassed_1btag_400 << " events)" << " Efficiency: " << 100.*eff_1btag_400 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_400 << " ev/fb-1  (" << nEventsPassed_2btag_400 << " events)" << " Efficiency: " << 100.*eff_2btag_400 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 450 GeV (423-495): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_450 << " ev/fb-1  (" << nEventsPassed_0btag_450 << " events)" << " Efficiency: " << 100.*eff_0btag_450 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_450 << " ev/fb-1  (" << nEventsPassed_1btag_450 << " events)" << " Efficiency: " << 100.*eff_1btag_450 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_450 << " ev/fb-1  (" << nEventsPassed_2btag_450 << " events)" << " Efficiency: " << 100.*eff_2btag_450 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 500 GeV (470-550): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_500 << " ev/fb-1  (" << nEventsPassed_0btag_500 << " events)" << " Efficiency: " << 100.*eff_0btag_500 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_500 << " ev/fb-1  (" << nEventsPassed_1btag_500 << " events)" << " Efficiency: " << 100.*eff_1btag_500 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_500 << " ev/fb-1  (" << nEventsPassed_2btag_500 << " events)" << " Efficiency: " << 100.*eff_2btag_500 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << "----> 600 GeV (564-660): " << std::endl;
  std::cout << "            0 btag: " << 1000.*nEventsPassed_fb_0btag_600 << " ev/fb-1  (" << nEventsPassed_0btag_600 << " events)" << " Efficiency: " << 100.*eff_0btag_600 << "%" << std::endl;
  std::cout << "            1 btag: " << 1000.*nEventsPassed_fb_1btag_600 << " ev/fb-1  (" << nEventsPassed_1btag_600 << " events)" << " Efficiency: " << 100.*eff_1btag_600 << "%" << std::endl;
  std::cout << "            2 btag: " << 1000.*nEventsPassed_fb_2btag_600 << " ev/fb-1  (" << nEventsPassed_2btag_600 << " events)" << " Efficiency: " << 100.*eff_2btag_600 << "%" << std::endl;
  std::cout << "----------------------------------------------------------------" << std::endl;
  std::cout << std::endl;



  outFile_->cd();

  tree_passedEvents->Write();

  h1_nCounter->Write();
  h1_nCounterW->Write();
  h1_nCounterPU->Write();

  h1_nEventsCategories_presel->Write();

  h1_nEvents_fb_gluetag_250->Write();
  h1_nEvents_fb_0btag_250->Write();
  h1_nEvents_fb_1btag_250->Write();
  h1_nEvents_fb_2btag_250->Write();

  h1_nEvents_fb_gluetag_250_ELE->Write();
  h1_nEvents_fb_0btag_250_ELE->Write();
  h1_nEvents_fb_1btag_250_ELE->Write();
  h1_nEvents_fb_2btag_250_ELE->Write();

  h1_nEvents_fb_gluetag_250_MU->Write();
  h1_nEvents_fb_0btag_250_MU->Write();
  h1_nEvents_fb_1btag_250_MU->Write();
  h1_nEvents_fb_2btag_250_MU->Write();

  h1_eff_gluetag_250->Write();
  h1_eff_0btag_250->Write();
  h1_eff_1btag_250->Write();
  h1_eff_2btag_250->Write();

  h1_eff_gluetag_250_ELE->Write();
  h1_eff_0btag_250_ELE->Write();
  h1_eff_1btag_250_ELE->Write();
  h1_eff_2btag_250_ELE->Write();

  h1_eff_gluetag_250_MU->Write();
  h1_eff_0btag_250_MU->Write();
  h1_eff_1btag_250_MU->Write();
  h1_eff_2btag_250_MU->Write();

  h1_nEvents_fb_gluetag_300->Write();
  h1_nEvents_fb_0btag_300->Write();
  h1_nEvents_fb_1btag_300->Write();
  h1_nEvents_fb_2btag_300->Write();

  h1_nEvents_fb_gluetag_300_ELE->Write();
  h1_nEvents_fb_0btag_300_ELE->Write();
  h1_nEvents_fb_1btag_300_ELE->Write();
  h1_nEvents_fb_2btag_300_ELE->Write();

  h1_nEvents_fb_gluetag_300_MU->Write();
  h1_nEvents_fb_0btag_300_MU->Write();
  h1_nEvents_fb_1btag_300_MU->Write();
  h1_nEvents_fb_2btag_300_MU->Write();

  h1_eff_gluetag_300->Write();
  h1_eff_0btag_300->Write();
  h1_eff_1btag_300->Write();
  h1_eff_2btag_300->Write();

  h1_eff_gluetag_300_ELE->Write();
  h1_eff_0btag_300_ELE->Write();
  h1_eff_1btag_300_ELE->Write();
  h1_eff_2btag_300_ELE->Write();

  h1_eff_gluetag_300_MU->Write();
  h1_eff_0btag_300_MU->Write();
  h1_eff_1btag_300_MU->Write();
  h1_eff_2btag_300_MU->Write();

  h1_nEvents_fb_gluetag_350->Write();
  h1_nEvents_fb_0btag_350->Write();
  h1_nEvents_fb_1btag_350->Write();
  h1_nEvents_fb_2btag_350->Write();

  h1_nEvents_fb_gluetag_350_ELE->Write();
  h1_nEvents_fb_0btag_350_ELE->Write();
  h1_nEvents_fb_1btag_350_ELE->Write();
  h1_nEvents_fb_2btag_350_ELE->Write();

  h1_nEvents_fb_gluetag_350_MU->Write();
  h1_nEvents_fb_0btag_350_MU->Write();
  h1_nEvents_fb_1btag_350_MU->Write();
  h1_nEvents_fb_2btag_350_MU->Write();

  h1_eff_gluetag_350->Write();
  h1_eff_0btag_350->Write();
  h1_eff_1btag_350->Write();
  h1_eff_2btag_350->Write();

  h1_eff_gluetag_350_ELE->Write();
  h1_eff_0btag_350_ELE->Write();
  h1_eff_1btag_350_ELE->Write();
  h1_eff_2btag_350_ELE->Write();

  h1_eff_gluetag_350_MU->Write();
  h1_eff_0btag_350_MU->Write();
  h1_eff_1btag_350_MU->Write();
  h1_eff_2btag_350_MU->Write();

  h1_nEvents_fb_gluetag_400->Write();
  h1_nEvents_fb_0btag_400->Write();
  h1_nEvents_fb_1btag_400->Write();
  h1_nEvents_fb_2btag_400->Write();

  h1_nEvents_fb_gluetag_400_ELE->Write();
  h1_nEvents_fb_0btag_400_ELE->Write();
  h1_nEvents_fb_1btag_400_ELE->Write();
  h1_nEvents_fb_2btag_400_ELE->Write();

  h1_nEvents_fb_gluetag_400_MU->Write();
  h1_nEvents_fb_0btag_400_MU->Write();
  h1_nEvents_fb_1btag_400_MU->Write();
  h1_nEvents_fb_2btag_400_MU->Write();

  h1_eff_gluetag_400->Write();
  h1_eff_0btag_400->Write();
  h1_eff_1btag_400->Write();
  h1_eff_2btag_400->Write();

  h1_eff_gluetag_400_ELE->Write();
  h1_eff_0btag_400_ELE->Write();
  h1_eff_1btag_400_ELE->Write();
  h1_eff_2btag_400_ELE->Write();

  h1_eff_gluetag_400_MU->Write();
  h1_eff_0btag_400_MU->Write();
  h1_eff_1btag_400_MU->Write();
  h1_eff_2btag_400_MU->Write();

  h1_nEvents_fb_gluetag_450->Write();
  h1_nEvents_fb_0btag_450->Write();
  h1_nEvents_fb_1btag_450->Write();
  h1_nEvents_fb_2btag_450->Write();

  h1_nEvents_fb_gluetag_450_ELE->Write();
  h1_nEvents_fb_0btag_450_ELE->Write();
  h1_nEvents_fb_1btag_450_ELE->Write();
  h1_nEvents_fb_2btag_450_ELE->Write();

  h1_nEvents_fb_gluetag_450_MU->Write();
  h1_nEvents_fb_0btag_450_MU->Write();
  h1_nEvents_fb_1btag_450_MU->Write();
  h1_nEvents_fb_2btag_450_MU->Write();

  h1_eff_gluetag_450->Write();
  h1_eff_0btag_450->Write();
  h1_eff_1btag_450->Write();
  h1_eff_2btag_450->Write();

  h1_eff_gluetag_450_ELE->Write();
  h1_eff_0btag_450_ELE->Write();
  h1_eff_1btag_450_ELE->Write();
  h1_eff_2btag_450_ELE->Write();

  h1_eff_gluetag_450_MU->Write();
  h1_eff_0btag_450_MU->Write();
  h1_eff_1btag_450_MU->Write();
  h1_eff_2btag_450_MU->Write();

  h1_nEvents_fb_gluetag_500->Write();
  h1_nEvents_fb_0btag_500->Write();
  h1_nEvents_fb_1btag_500->Write();
  h1_nEvents_fb_2btag_500->Write();

  h1_nEvents_fb_gluetag_500_ELE->Write();
  h1_nEvents_fb_0btag_500_ELE->Write();
  h1_nEvents_fb_1btag_500_ELE->Write();
  h1_nEvents_fb_2btag_500_ELE->Write();

  h1_nEvents_fb_gluetag_500_MU->Write();
  h1_nEvents_fb_0btag_500_MU->Write();
  h1_nEvents_fb_1btag_500_MU->Write();
  h1_nEvents_fb_2btag_500_MU->Write();

  h1_eff_gluetag_500->Write();
  h1_eff_0btag_500->Write();
  h1_eff_1btag_500->Write();
  h1_eff_2btag_500->Write();

  h1_eff_gluetag_500_ELE->Write();
  h1_eff_0btag_500_ELE->Write();
  h1_eff_1btag_500_ELE->Write();
  h1_eff_2btag_500_ELE->Write();

  h1_eff_gluetag_500_MU->Write();
  h1_eff_0btag_500_MU->Write();
  h1_eff_1btag_500_MU->Write();
  h1_eff_2btag_500_MU->Write();

  h1_nEvents_fb_gluetag_600->Write();
  h1_nEvents_fb_0btag_600->Write();
  h1_nEvents_fb_1btag_600->Write();
  h1_nEvents_fb_2btag_600->Write();

  h1_nEvents_fb_gluetag_600_ELE->Write();
  h1_nEvents_fb_0btag_600_ELE->Write();
  h1_nEvents_fb_1btag_600_ELE->Write();
  h1_nEvents_fb_2btag_600_ELE->Write();

  h1_nEvents_fb_gluetag_600_MU->Write();
  h1_nEvents_fb_0btag_600_MU->Write();
  h1_nEvents_fb_1btag_600_MU->Write();
  h1_nEvents_fb_2btag_600_MU->Write();

  h1_eff_gluetag_600->Write();
  h1_eff_0btag_600->Write();
  h1_eff_1btag_600->Write();
  h1_eff_2btag_600->Write();

  h1_eff_gluetag_600_ELE->Write();
  h1_eff_0btag_600_ELE->Write();
  h1_eff_1btag_600_ELE->Write();
  h1_eff_2btag_600_ELE->Write();

  h1_eff_gluetag_600_MU->Write();
  h1_eff_0btag_600_MU->Write();
  h1_eff_1btag_600_MU->Write();
  h1_eff_2btag_600_MU->Write();


  h1_run->Write();

  h1_nvertex->Write();
  h1_nvertex_PUW->Write();
  h1_nvertex_PUW_ave->Write();

  h1_rhoPF_presel->Write();
  h1_rhoPF->Write();
  
  h1_pfMet->Write();
  h1_metSignificance->Write();
  h1_metSignificance_2btag->Write();
  h1_mEtSig->Write();
  h1_mEtSig_2btag->Write();
  h1_pfMet_2btag->Write();
  h1_pfMetOverMZZ->Write();
  h1_pfMetOverMZZ_2btag->Write();

  h1_nCandidates->Write();

  h1_ptJet_all_presel->Write();
  h1_etaJet_all_presel->Write();
  h1_nJets_presel->Write();
  h1_nPairs_presel->Write();

  h1_ptResoJet1_beforeKin->Write();
  h1_ptResoJet2_beforeKin->Write();
  h1_ptResoJet1_afterKin->Write();
  h1_ptResoJet2_afterKin->Write();

  h1_ptZreso_beforeKin->Write();
  h1_ptZreso_afterKin->Write();

  h1_mHreso_beforeKin->Write();
  h1_mHreso_afterKin->Write();

  h1_deltaRll_presel->Write();
  h1_deltaRjj_all_presel->Write();

  h1_ptLept1->Write();
  h1_ptLept2->Write();

  h1_ptLept1_presel->Write();
  h1_ptLept2_presel->Write();
  h1_etaLept1_presel->Write();
  h1_etaLept2_presel->Write();

  h1_mZll->Write();
  h1_mZll_presel->Write();
  h1_mZmumu->Write();
  h1_mZmumu_presel->Write();
  h1_mZmumu_presel_0jets->Write();
  h1_mZee->Write();
  h1_mZee_presel->Write();
  h1_mZee_presel_0jets->Write();

  h1_mZjjMC->Write();

  h1_mZjj->Write();
  h1_mZjj_0btag->Write();
  h1_mZjj_1btag->Write();
  h1_mZjj_2btag->Write();
  h1_mZjj_loMass->Write();
  h1_mZjj_medMass->Write();
  h1_mZjj_hiMass->Write();
  h1_mZjj_MU->Write();
  h1_mZjj_ELE->Write();
  h1_mZjj_nogluetag->Write();
  h1_mZjj_nogluetag_MU->Write();
  h1_mZjj_nogluetag_ELE->Write();
  h1_mZjj_loChiSquareProb->Write();
  h1_mZjj_hiChiSquareProb->Write();
  h1_mZjj_all_presel->Write();

  h1_ptZll_presel->Write();
  h1_ptZjj_all_presel->Write();

  h1_deltaRjj_sidebands->Write();
  h1_deltaRjj_prekin_sidebands->Write();
  h1_deltaRjj->Write();
  h1_deltaRjj_prekin->Write();

  h1_ptZjj->Write();
  h1_ptZll->Write();

  h1_cosThetaStar->Write();
  h1_cosTheta1->Write();
  h1_cosTheta2->Write();
  h1_phi->Write();
  h1_phi1->Write();

  h1_cosThetaStar_kinfit->Write();
  h1_cosTheta1_kinfit->Write();
  h1_cosTheta2_kinfit->Write();
  h1_phi_kinfit->Write();
  h1_phi1_kinfit->Write();

  h1_kinfit_chiSquare->Write();
  h1_kinfit_chiSquareProb->Write();
  
  h1_helicityLD->Write();
  h1_helicityLD_nogluetag->Write();
  h1_helicityLD_nokinfit->Write();
  h1_helicityLD_sidebands->Write();
  h1_helicityLD_MW200->Write();
  h1_helicityLD_MW250->Write();
  h1_helicityLD_MW300->Write();
  h1_helicityLD_MW400->Write();
  h1_helicityLD_MW500->Write();
  h1_helicityLD_kinfit->Write();

  h2_helicityLD_vs_mZZ->Write();

  h1_mZZ_hiChiSquareProb->Write();
  h1_mZZ_loChiSquareProb->Write();
  h1_mZZ_mZjj_cut->Write();
  h1_mZZ_mZjj_notcut->Write();
  h1_mZZ_ZjjMassConstr_hiMass->Write();
  h1_mZZ_nokinfit_hiMass_all->Write();
  h1_mZZ_kinfit_hiMass_all->Write();
  h1_mZZ_kinfit_hiMass_hiPU->Write();
  h1_mZZ_kinfit_hiMass_loPU->Write();
  h1_mZZ_kinfit_hiMass_nogluetag->Write();
  h1_mZZ_kinfit_hiMass_gluetag->Write();
  h1_mZZ_kinfit_hiMass_gluetag_ELE->Write();
  h1_mZZ_kinfit_hiMass_gluetag_MU->Write();
  h1_mZZ_kinfit_hiMass_0btag->Write();
  h1_mZZ_kinfit_hiMass_0btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_0btag_MU->Write();
//h1_mZZ_nokinfit_hiMass_0btag->Write();
//h1_mZZ_nokinfit_hiMass_0btag_ELE->Write();
//h1_mZZ_nokinfit_hiMass_0btag_MU->Write();
  h1_mZZ_kinfit_hiMass_1btag->Write();
  h1_mZZ_kinfit_hiMass_1btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_1btag_MU->Write();
//h1_mZZ_nokinfit_hiMass_1btag->Write();
//h1_mZZ_nokinfit_hiMass_1btag_ELE->Write();
//h1_mZZ_nokinfit_hiMass_1btag_MU->Write();
  h1_mZZ_kinfit_hiMass_2btag->Write();
  h1_mZZ_kinfit_hiMass_2btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_2btag_MU->Write();
//h1_mZZ_nokinfit_hiMass_2btag->Write();
//h1_mZZ_nokinfit_hiMass_2btag_ELE->Write();
//h1_mZZ_nokinfit_hiMass_2btag_MU->Write();
  h1_mZZ_kinfit_hiMass_hiQG->Write();
  h1_mZZ_kinfit_hiMass_loQG->Write();

  h1_mZZ_kinfit_hiMass_sidebands_gluetag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_gluetag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_gluetag_MU->Write();
  h1_mZZ_kinfit_hiMass_sidebands_0btag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_0btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_0btag_MU->Write();
  h1_mZZ_kinfit_hiMass_sidebands_1btag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_1btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_1btag_MU->Write();
  h1_mZZ_kinfit_hiMass_sidebands_2btag->Write();
  h1_mZZ_kinfit_hiMass_sidebands_2btag_ELE->Write();
  h1_mZZ_kinfit_hiMass_sidebands_2btag_MU->Write();

  h1_deltaRZmatching->Write();
  h1_mZZ_kinfit_hiMass_0btag_matched->Write();

  h1_ptZZ->Write();
  h1_ptZZ_kinfit->Write();
  h1_etaZZ->Write();
  h1_etaZZ_kinfit->Write();

  h1_ptJetRecoil->Write();
  h1_ptHiggs->Write();

  h1_deltaR_part1->Write();
  h1_ptJet1->Write();
  h1_ptJet1_prekin->Write();
  h1_etaJet1->Write();
  h1_eElectronsJet1->Write();
  h1_eMuonsJet1->Write();
  h1_partFlavorJet1->Write();
  h1_partFlavorJet1_MW400->Write();
  h1_partFlavorJet1_MW500->Write();
  h1_partFlavorJet1_matched->Write();
  h1_partFlavorJet1_notmatched->Write();
  h1_partFlavorJet1_0btag->Write();
  h1_partFlavorJet1_1btag->Write();
  h1_partFlavorJet1_2btag->Write();
  h1_partFlavorJet1_gluetag->Write();
  h1_QGLikelihoodJet1_gluonMatched->Write();
  h1_QGLikelihoodJet1_quarkMatched->Write();
  h1_tcheJet1->Write();

  h1_tcheJet->Write();

  h1_QGLikelihoodJetRecoil->Write();
  h1_QGLikelihoodJetRecoil_MW400->Write();
  h1_QGLikelihoodProdRecoil->Write();
  h1_QGLikelihoodProdRecoil_MW400->Write();

  h1_deltaR_part2->Write();
  h1_ptJet2->Write();
  h1_ptJet2_prekin->Write();
  h1_etaJet2->Write();
  h1_eElectronsJet2->Write();
  h1_eMuonsJet2->Write();
  h1_partFlavorJet2->Write();
  h1_partFlavorJet2_MW400->Write();
  h1_partFlavorJet2_MW500->Write();
  h1_partFlavorJet2_matched->Write();
  h1_partFlavorJet2_notmatched->Write();
  h1_partFlavorJet2_0btag->Write();
  h1_partFlavorJet2_1btag->Write();
  h1_partFlavorJet2_2btag->Write();
  h1_partFlavorJet2_gluetag->Write();
  h1_QGLikelihoodJet2_gluonMatched->Write();
  h1_QGLikelihoodJet2_quarkMatched->Write();
  h1_tcheJet2->Write();

  h1_QGLikelihoodProd_oneGluon->Write();
  h1_mZZ_kinfit_hiMass_oneGluon->Write();
  h1_QGLikelihoodProd_twoGluon->Write();
  h1_mZZ_kinfit_hiMass_twoGluon->Write();
  h1_QGLikelihoodProd_twoQuark->Write();
  h1_QGLikelihoodProd_twoQuark_MW400->Write();
  h1_QGLikelihoodProd_twoQuark_MW500->Write();
  h1_mZZ_kinfit_hiMass_twoQuark->Write();

  h1_nEvents_partFlavor_beforeQG->Write();
  h1_nEvents_partFlavor_afterQG->Write();

  h1_nEventsMW400_partFlavor_beforeQG->Write();
  h1_nEventsMW400_partFlavor_afterQG->Write();

  h1_nEventsMW500_partFlavor_beforeQG->Write();
  h1_nEventsMW500_partFlavor_afterQG->Write();


  h1_deltaRZZ->Write();

//h1_mZZ_MCassoc->Write();
//h1_mZZ_MCassoc_ZjjMassConstr->Write();
//h1_mZZ_MCassoc_kinfit->Write();
//h1_mZZ_MCassoc_kinfit_cands->Write();

  h2_mZjj_vs_mZZ->Write();
  h2_mZjj_vs_mZZ_0btag->Write();
  h2_mZjj_vs_mZZ_1btag->Write();
  h2_mZjj_vs_mZZ_2btag->Write();
  h2_mZjj_vs_mZZ_gluetag->Write();

  h2_mZjj_vs_mZZ_kinfit->Write();
  h2_mZjj_vs_mZZ_kinfit_0btag->Write();
  h2_mZjj_vs_mZZ_kinfit_1btag->Write();
  h2_mZjj_vs_mZZ_kinfit_2btag->Write();
  h2_mZjj_vs_mZZ_kinfit_gluetag->Write();


  h1_deltaE_ch->Write();
  h1_deltaE_gamma->Write();
  h1_deltaE_nh->Write();
  h1_deltaEta_ch->Write();
  h1_deltaEta_gamma->Write();
  h1_deltaEta_nh->Write();
  h1_deltaPhi_ch->Write();
  h1_deltaPhi_gamma->Write();
  h1_deltaPhi_nh->Write();
  h1_deltaPt_ch->Write();
  h1_deltaPt_gamma->Write();
  h1_deltaPt_nh->Write();

  h1_ptDJet1->Write();
  h1_ptDJet2->Write();

  h1_nChargedJet1->Write();
  h1_nNeutralJet1->Write();
  h1_ptDJet1->Write();

  h1_nChargedJet2->Write();
  h1_nNeutralJet2->Write();
  h1_ptDJet2->Write();

  h1_QGLikelihoodJet1->Write();
  h1_QGLikelihoodJet2->Write();
  h1_QGLikelihoodProd->Write();

  h1_QGLikelihoodNoPUJet1->Write();
  h1_QGLikelihoodNoPUJet2->Write();
  h1_QGLikelihoodNoPUProd->Write();

  h1_QGLikelihoodJet1_MW300->Write();
  h1_QGLikelihoodJet2_MW300->Write();
  h1_QGLikelihoodProd_MW300->Write();

  h1_QGLikelihoodJet1_MW400->Write();
  h1_QGLikelihoodJet2_MW400->Write();
  h1_QGLikelihoodProd_MW400->Write();

  h1_QGLikelihoodJet1_MW500->Write();
  h1_QGLikelihoodJet2_MW500->Write();
  h1_QGLikelihoodProd_MW500->Write();


  h1_QGLikelihood_100_123->Write();
  h1_QGLikelihood_66_81->Write();
//outFile_->mkdir("QGbins");
//outFile_->cd("QGbins");

//for( unsigned iPtBin=0; iPtBin<nPtBins; ++iPtBin ) {

//  vh1_rmsCandJet1[iPtBin]->Write();
//  vh1_ptDJet1[iPtBin]->Write();
//  vh1_nChargedJet1[iPtBin]->Write();
//  vh1_nNeutralJet1[iPtBin]->Write();
//  vh1_QGLikelihoodJet1[iPtBin]->Write();

//  vh1_rmsCandJet2[iPtBin]->Write();
//  vh1_ptDJet2[iPtBin]->Write();
//  vh1_nChargedJet2[iPtBin]->Write();
//  vh1_nNeutralJet2[iPtBin]->Write();
//  vh1_QGLikelihoodJet2[iPtBin]->Write();

//}


  outFile_->Close();


} // finalize()



void Ntp1Finalizer_TTZTriplepton::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="presel" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.;
    helicityLD_intercept_1btags_ = 0.;
    helicityLD_intercept_2btags_ = 0.;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 9999999999.;
    use_looseBTags_ = true;

  } else if( selectionType_=="presel_invMZll" ) {

    ptLept1_thresh_ = 10.;
    ptLept2_thresh_ = 10.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = true;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.;
    helicityLD_intercept_1btags_ = 0.;
    helicityLD_intercept_2btags_ = 0.;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 9999999999.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_v1" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00124;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.1433;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 9999999999.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_v2" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00025;
    helicityLD_slope_1btags_ = 0.000656;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.55;
    helicityLD_intercept_1btags_ = 0.302;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_noBTagCat" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00025;
    helicityLD_slope_1btags_ = 0.00025;;
    helicityLD_slope_2btags_ = 0.00025;
    helicityLD_intercept_0btags_ = 0.55;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.55;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_mediumBTags_v1" ) { //"option 6"

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00054;
    helicityLD_slope_1btags_ = 0.00098;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.428;
    helicityLD_intercept_1btags_ = 0.156;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.5;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = false;

  } else if( selectionType_=="optLD_looseBTags_noQG" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00124;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.1433;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_v2_noQG" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.00025;
    helicityLD_slope_1btags_ = 0.000656;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.55;
    helicityLD_intercept_1btags_ = 0.302;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="optLD_looseBTags_fix400OLD" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    //helicityLD_intercept_0btags_ = 0.639;
    helicityLD_intercept_0btags_ = 0.65;
    helicityLD_intercept_1btags_ = 0.55;
    helicityLD_intercept_2btags_ = 0.5;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else if( selectionType_=="noCutLD" ) {

    ptLept1_thresh_ = 40.;
    ptLept2_thresh_ = 20.;
    etaLept1_thresh_ = 3.;
    etaLept2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    invert_mZll_ = false;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;
    mZjj_threshLo_ = 75.;
    mZjj_threshHi_ = 105.;
    helicityLD_slope_0btags_ = 0.;
    helicityLD_slope_1btags_ = 0.;
    helicityLD_slope_2btags_ = 0.;
    helicityLD_intercept_0btags_ = 0.;
    helicityLD_intercept_1btags_ = 0.;
    helicityLD_intercept_2btags_ = 0.;
    helicityLD_minThresh_0btags_ = 0.;
    helicityLD_minThresh_1btags_ = 0.;
    helicityLD_minThresh_2btags_ = 0.;
    helicityLD_maxThresh_0btags_ = 1.;
    helicityLD_maxThresh_1btags_ = 1.;
    helicityLD_maxThresh_2btags_ = 1.;
    QGLikelihoodProd_thresh_ = 0.1;
    metSignificance_thresh_ = 10.;
    use_looseBTags_ = true;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType



