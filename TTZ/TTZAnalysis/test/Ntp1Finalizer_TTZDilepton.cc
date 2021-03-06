#include "Ntp1Finalizer_TTZDilepton.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"
#include "KinematicFit/DiJetKinFitter.h"

#include "PUWeight.h"




bool USE_MC_MASS=false;

int DEBUG_EVENTNUMBER = 98901397;






// constructor:

Ntp1Finalizer_TTZDilepton::Ntp1Finalizer_TTZDilepton( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, const std::string& PUType, const std::string& leptType ) : Ntp1Finalizer( "TTZDilepton", dataset, leptType ) {

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




void Ntp1Finalizer_TTZDilepton::finalize() {

  //if( outFile_==0 ) this->createOutputFile();
  
  Int_t run;
  tree_->SetBranchAddress("run", &run);
  tree_->GetEntry(0);
  bool isMC = (run < 160000);
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


  TH1D* h1_deltaRZll = new TH1D("deltaRZll", "", 500, 0., 5.);
  h1_deltaRZll->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 400., 0., 400.);
  h1_ptZll->Sumw2();
  TH1D* h1_etaZll = new TH1D("etaZll", "", 200, -5., 5.);
  h1_etaZll->Sumw2();
  TH1D* h1_mZll = new TH1D("mZll", "", 240, 40., 160.);
  h1_mZll->Sumw2();



  TH1D* h1_nJets = new TH1D("nJets", "", 7, 5.5, 12.5);
  h1_nJets->Sumw2();


  TH1D* h1_bTagJet1 = new TH1D("bTagJet1", "", 420, -1., 20.);
  h1_bTagJet1->Sumw2();
  TH1D* h1_bTagJet2 = new TH1D("bTagJet2", "", 420, -1., 20.);
  h1_bTagJet2->Sumw2();


  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 400, 20., 420.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 400, 20., 420.);
  h1_ptJet2->Sumw2();
  TH1D* h1_ptJet3 = new TH1D("ptJet3", "", 400, 20., 420.);
  h1_ptJet3->Sumw2();
  TH1D* h1_ptJet4 = new TH1D("ptJet4", "", 400, 20., 420.);
  h1_ptJet4->Sumw2();
  TH1D* h1_ptJet5 = new TH1D("ptJet5", "", 400, 20., 420.);
  h1_ptJet5->Sumw2();
  TH1D* h1_ptJet6 = new TH1D("ptJet6", "", 400, 20., 420.);
  h1_ptJet6->Sumw2();

  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 200, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 200, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_etaJet3 = new TH1D("etaJet3", "", 200, -5., 5.);
  h1_etaJet3->Sumw2();
  TH1D* h1_etaJet4 = new TH1D("etaJet4", "", 200, -5., 5.);
  h1_etaJet4->Sumw2();
  TH1D* h1_etaJet5 = new TH1D("etaJet5", "", 200, -5., 5.);
  h1_etaJet5->Sumw2();
  TH1D* h1_etaJet6 = new TH1D("etaJet6", "", 200, -5., 5.);
  h1_etaJet6->Sumw2();

  TH1D* h1_partFlavorJet1 = new TH1D("partFlavorJet1", "", 38, -15.5, 22.5);
  h1_partFlavorJet1->Sumw2();
  TH1D* h1_partFlavorJet2 = new TH1D("partFlavorJet2", "", 38, -15.5, 22.5);
  h1_partFlavorJet2->Sumw2();
  TH1D* h1_partFlavorJet3 = new TH1D("partFlavorJet3", "", 38, -15.5, 22.5);
  h1_partFlavorJet3->Sumw2();
  TH1D* h1_partFlavorJet4 = new TH1D("partFlavorJet4", "", 38, -15.5, 22.5);
  h1_partFlavorJet4->Sumw2();
  TH1D* h1_partFlavorJet5 = new TH1D("partFlavorJet5", "", 38, -15.5, 22.5);
  h1_partFlavorJet5->Sumw2();
  TH1D* h1_partFlavorJet6 = new TH1D("partFlavorJet6", "", 38, -15.5, 22.5);
  h1_partFlavorJet6->Sumw2();


  TH1D* h1_deltaRbb = new TH1D("deltaRbb", "", 500, 0., 5.);
  h1_deltaRbb->Sumw2();
  
  TH1D* h1_mb1jj = new TH1D("mb1jj", "", 500., 0., 500.);
  h1_mb1jj->Sumw2();
  TH1D* h1_mb2jj = new TH1D("mb2jj", "", 500., 0., 500.);
  h1_mb2jj->Sumw2();
  TH1D* h1_mbjj_best = new TH1D("mbjj_best", "", 500., 0., 500.);
  h1_mbjj_best->Sumw2();
  
  TH1D* h1_mb1jjZ = new TH1D("mb1jjZ", "", 500., 0., 500.);
  h1_mb1jjZ->Sumw2();
  TH1D* h1_mb2jjZ = new TH1D("mb2jjZ", "", 500., 0., 500.);
  h1_mb2jjZ->Sumw2();
  TH1D* h1_mbjjZ_best = new TH1D("mbjjZ_best", "", 500., 0., 500.);
  h1_mbjjZ_best->Sumw2();




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
  float tmass = 172.9;


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
//} else if( dataset_tstr.Contains("Fall11") ) {
//  TFile* filePUMC = TFile::Open("s6MCPileUp.root");
//  TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("pileup");
//  std::cout << "-> Switching MC PU file to: s6MCPileUp.root" << std::endl;
//  fPUWeight->SetMCHistogram(h1_nPU_mc);
  }



  float mZll_t, ptZll_t;
  float ptLeptZ1_t, ptLeptZ2_t, etaLeptZ1_t, etaLeptZ2_t;
  float ptJet1_t, ptJet2_t, ptJet3_t, ptJet4_t, ptJet5_t, ptJet6_t;
  float etaJet1_t, etaJet2_t, etaJet3_t, etaJet4_t, etaJet5_t, etaJet6_t;
  float bTagJet1_t, bTagJet2_t;
  float HLTSF;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "leptType", &leptType, "leptType/I" );
  tree_passedEvents->Branch( "ptLeptZ1", &ptLeptZ1_t, "ptLeptZ1_t/F" );
  tree_passedEvents->Branch( "ptLeptZ2", &ptLeptZ2_t, "ptLeptZ2_t/F" );
  tree_passedEvents->Branch( "etaLeptZ1", &etaLeptZ1_t, "etaLeptZ1_t/F" );
  tree_passedEvents->Branch( "etaLeptZ2", &etaLeptZ2_t, "etaLeptZ2_t/F" );
  tree_passedEvents->Branch( "mZll", &mZll_t, "mZll_t/F" );
  tree_passedEvents->Branch( "ptZll", &ptZll_t, "ptZll_t/F" );
  tree_passedEvents->Branch( "ptJet1", &ptJet1_t, "ptJet1_t/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2_t, "ptJet2_t/F" );
  tree_passedEvents->Branch( "bTagJet1", &bTagJet1_t, "bTagJet1_t/F" );
  tree_passedEvents->Branch( "bTagJet2", &bTagJet2_t, "bTagJet2_t/F" );
  tree_passedEvents->Branch( "ptJet3", &ptJet3_t, "ptJet3_t/F" );
  tree_passedEvents->Branch( "ptJet4", &ptJet4_t, "ptJet4_t/F" );
  tree_passedEvents->Branch( "ptJet5", &ptJet5_t, "ptJet5_t/F" );
  tree_passedEvents->Branch( "ptJet6", &ptJet6_t, "ptJet6_t/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1_t, "etaJet1_t/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2_t, "etaJet2_t/F" );
  tree_passedEvents->Branch( "etaJet3", &etaJet3_t, "etaJet3_t/F" );
  tree_passedEvents->Branch( "etaJet4", &etaJet4_t, "etaJet4_t/F" );
  tree_passedEvents->Branch( "etaJet5", &etaJet5_t, "etaJet5_t/F" );
  tree_passedEvents->Branch( "etaJet6", &etaJet6_t, "etaJet6_t/F" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "HLTSF", &HLTSF, "HLTSF/F" );
  tree_passedEvents->Branch( "PUWeight", &eventWeightPU, "eventWeightPU/F" );




ofstream ofs("run_event.txt");




  std::cout << std::endl << std::endl;
  std::cout << "+++ BEGINNING ANALYSIS LOOP" << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  std::cout << "----> B-TAGGER: " << bTaggerType_ << std::endl;
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


      //eventWeight *= fPUWeight->GetWeight(nPU);

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


    } //if is not mc



    // this is dilepton channel: no other lepton in the event 
    if( nLept!=0 ) continue;


    h1_rhoPF_presel->Fill( rhoPF, eventWeight);


    h1_pfMet->Fill( pfMet, eventWeight );
    h1_metSignificance->Fill( metSignificance, eventWeight );



    TLorentzVector leptZ1, leptZ2;
    leptZ1.SetPtEtaPhiE( ptLeptZ1, etaLeptZ1, phiLeptZ1, eLeptZ1 );
    leptZ2.SetPtEtaPhiE( ptLeptZ2, etaLeptZ2, phiLeptZ2, eLeptZ2 );


    TLorentzVector diLepton = leptZ1+leptZ2;


    h1_ptLeptZ1->Fill( leptZ1.Pt(), eventWeight );
    h1_ptLeptZ2->Fill( leptZ2.Pt(), eventWeight );
    h1_etaLeptZ1->Fill( leptZ1.Eta(), eventWeight );
    h1_etaLeptZ2->Fill( leptZ2.Eta(), eventWeight );

    h1_deltaRZll->Fill( leptZ2.DeltaR(leptZ2), eventWeight );

    h1_ptZll->Fill( diLepton.Pt(), eventWeight );
    h1_etaZll->Fill( diLepton.Eta(), eventWeight );
    h1_mZll->Fill( diLepton.M(), eventWeight );







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




    if( nJets<6 ) continue;

    h1_nJets->Fill( nJets , eventWeight );

    std::vector<AnalysisJet> jets;

  //// default: order by pt
  //for( unsigned int iJet=0; iJet<6; ++iJet ) {
  //  AnalysisJet* newJet = new AnalysisJet();
  //  newJet->SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);
  //  jets.push_back(*newJet);
  //}


//std::cout << "check jets: " << std::endl;
//for( unsigned iSelectedJet=0; iSelectedJet<6; ++iSelectedJet )
//std::cout << iSelectedJet << " " << jets[iSelectedJet].Pt() << std::endl;



    float bestBtag=-9999.;
    float bestBtag2=-9999.;
    int i_bestBtag=-1;
    int i_bestBtag2=-1;

    for( unsigned iJet=0; iJet<nJets; ++iJet) {

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


      if( thisBtag > bestBtag ) {

        bestBtag2 = bestBtag;
        bestBtag = thisBtag;
        i_bestBtag2 = i_bestBtag;
        i_bestBtag = iJet;

        if( jets.size()==0 ) {
          AnalysisJet* newJet = new AnalysisJet(thisJet);
          jets.push_back(*newJet);
        } else if( jets.size()==1 ) {
          AnalysisJet oldJet = jets[0];
          AnalysisJet* newJet = new AnalysisJet(thisJet);
          jets.push_back(*newJet);
          jets[0] = *newJet;
          jets[1] = oldJet;
        } else if( jets.size()==2 ) {
          jets[1] = jets[0];
          jets[0] = thisJet;
        }

      } else if( thisBtag > bestBtag2 ) { //means that at least one jet was found

        bestBtag2 = thisBtag;
        i_bestBtag2 = iJet;
        
        if( jets.size()==1 ) {
          AnalysisJet* newJet = new AnalysisJet(thisJet);
          jets.push_back(*newJet);
        } else if( jets.size()==2 ) {
          jets[1] = thisJet;
        }
  
      }

    } // for jets



    // now add other jets ordered in pt:
    for( unsigned iJet=0; iJet<nJets; ++iJet) {

      if( iJet==i_bestBtag || iJet==i_bestBtag2 ) continue;

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


      AnalysisJet* newJet = new AnalysisJet(thisJet);
      jets.push_back(*newJet);

    } //for additional jets
 
        
 // std::cout << "new order: " << std::endl;
 // for( unsigned iSelectedJet=0; iSelectedJet<6; ++iSelectedJet )
 // std::cout << iSelectedJet << " " << jets[iSelectedJet].Pt()  << std::endl;
 // std::cout << std::endl << std::endl << std::endl;
      //// slide them all:
      //for( unsigned iSelectedJet=0; iSelectedJet<5; ++iSelectedJet )
      //  jets[5-iSelectedJet] = jets[4-iSelectedJet];
      //jets[0] = thisJet;

    h1_ptJet1->Fill( jets[0].Pt(), eventWeight );
    h1_ptJet2->Fill( jets[1].Pt(), eventWeight );
    h1_ptJet3->Fill( jets[2].Pt(), eventWeight );
    h1_ptJet4->Fill( jets[3].Pt(), eventWeight );
    h1_ptJet5->Fill( jets[4].Pt(), eventWeight );
    h1_ptJet6->Fill( jets[5].Pt(), eventWeight );

    h1_etaJet1->Fill( jets[0].Eta(), eventWeight );
    h1_etaJet2->Fill( jets[1].Eta(), eventWeight );
    h1_etaJet3->Fill( jets[2].Eta(), eventWeight );
    h1_etaJet4->Fill( jets[3].Eta(), eventWeight );
    h1_etaJet5->Fill( jets[4].Eta(), eventWeight );
    h1_etaJet6->Fill( jets[5].Eta(), eventWeight );


    h1_partFlavorJet1->Fill( jets[0].pdgIdPart, eventWeight );
    h1_partFlavorJet2->Fill( jets[1].pdgIdPart, eventWeight );
    h1_partFlavorJet3->Fill( jets[2].pdgIdPart, eventWeight );
    h1_partFlavorJet4->Fill( jets[3].pdgIdPart, eventWeight );
    h1_partFlavorJet5->Fill( jets[4].pdgIdPart, eventWeight );
    h1_partFlavorJet6->Fill( jets[5].pdgIdPart, eventWeight );


    float bTaggerJet1, bTaggerJet2;
    if( bTaggerType_=="TCHE" ) {
      bTaggerJet1 = jets[0].trackCountingHighEffBJetTag;
      bTaggerJet2 = jets[1].trackCountingHighEffBJetTag;
    } else if( bTaggerType_=="SSVHE" ) {
      bTaggerJet1 = jets[0].simpleSecondaryVertexHighEffBJetTag;
      bTaggerJet2 = jets[1].simpleSecondaryVertexHighEffBJetTag;
    }


    h1_bTagJet1->Fill( bTaggerJet1, eventWeight );
    h1_bTagJet2->Fill( bTaggerJet1, eventWeight );



    h1_deltaRbb->Fill( jets[0].DeltaR(jets[1]), eventWeight );

/*
    TLorentzVector b1jj = jetB1 + jet3 + jet4;
    TLorentzVector b2jj = jetB2 + jet3 + jet4;

    TLorentzVector b1jjZ = jetB1 + jet3 + jet4 + diLepton;
    TLorentzVector b2jjZ = jetB2 + jet3 + jet4 + diLepton;

    h1_mb1jj->Fill( b1jj.M(), eventWeight );
    h1_mb2jj->Fill( b2jj.M(), eventWeight );

    if( fabs(b1jj.M()-tmass) < fabs(b2jj.M()-tmass) )
      h1_mbjj_best->Fill( b1jj.M(), eventWeight );
    else
      h1_mbjj_best->Fill( b2jj.M(), eventWeight );


    h1_mb1jjZ->Fill( b1jjZ.M(), eventWeight );
    h1_mb2jjZ->Fill( b2jjZ.M(), eventWeight );

    if( fabs(b1jjZ.M()-tmass) < fabs(b2jjZ.M()-tmass) )
      h1_mbjjZ_best->Fill( b1jjZ.M(), eventWeight );
    else
      h1_mbjjZ_best->Fill( b2jjZ.M(), eventWeight );



    TLorentzVector b1W = jetB1 + lept3 + met;
    TLorentzVector b2W = jetB2 + lept3 + met;

    TLorentzVector b1WZ = jetB1 + diLepton + lept3 + met;
    TLorentzVector b2WZ = jetB2 + diLepton + lept3 + met;

    h1_mTb1W->Fill( b1W.Mt(), eventWeight );
    h1_mTb2W->Fill( b2W.Mt(), eventWeight );

    if( fabs(b1W.Mt()-tmass) < fabs(b2W.Mt()-tmass) )
      h1_mTbW_best->Fill( b1W.Mt(), eventWeight );
    else
      h1_mTbW_best->Fill( b2W.Mt(), eventWeight );

    
    h1_mTb1WZ->Fill( b1WZ.Mt(), eventWeight );
    h1_mTb2WZ->Fill( b2WZ.Mt(), eventWeight );

    if( fabs(b1WZ.Mt()-tmass) < fabs(b2WZ.Mt()-tmass) )
      h1_mTbWZ_best->Fill( b1WZ.Mt(), eventWeight );
    else
      h1_mTbWZ_best->Fill( b2WZ.Mt(), eventWeight );
    
*/




//  // -------------------------
//  // KINEMATIC SELECTION: JETS
//  // -------------------------

//  if( jet1.Pt() < ptJet1_thresh_ ) continue;
//  if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet pt OK" << std::endl;
//  if( jet2.Pt() < ptJet2_thresh_ ) continue;
//  if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet pt OK" << std::endl;
//  if( fabs(jet1.Eta()) > etaJet1_thresh_ ) continue;
//  if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet eta OK" << std::endl;
//  if( fabs(jet2.Eta()) > etaJet2_thresh_ ) continue;
//  if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet eta OK" << std::endl;




    ptZll_t = diLepton.Pt(); 
    mZll_t = diLepton.M(); 

    ptLeptZ1_t = leptZ1.Pt();
    ptLeptZ2_t = leptZ2.Pt();
    etaLeptZ1_t = leptZ1.Eta();
    etaLeptZ2_t = leptZ2.Eta();

    ptJet1_t = jets[0].Pt();
    ptJet2_t = jets[1].Pt();
    ptJet3_t = jets[2].Pt();
    ptJet4_t = jets[3].Pt();
    ptJet5_t = jets[4].Pt();
    ptJet6_t = jets[5].Pt();

    bTagJet1_t = bTaggerJet1;
    bTagJet2_t = bTaggerJet2;

    etaJet1_t = jets[0].Eta();
    etaJet2_t = jets[1].Eta();
    etaJet3_t = jets[2].Eta();
    etaJet4_t = jets[3].Eta();
    etaJet5_t = jets[4].Eta();
    etaJet6_t = jets[5].Eta();



    // and fill tree:
    tree_passedEvents->Fill();

  

  
  } //for entries



  h1_nCounter->SetBinContent(1, nCounter_);
  h1_nCounterW->SetBinContent(1, nCounterW_);
  h1_nCounterPU->SetBinContent(1, nCounterPU_);




  // write all stuff in files:

  outFile_->cd();

  tree_passedEvents->Write();

  h1_nCounter->Write();
  h1_nCounterW->Write();
  h1_nCounterPU->Write();



  h1_nvertex->Write();
  h1_nvertex_PUW->Write();
  h1_nvertex_PUW_ave->Write();

  h1_pfMet->Write();

  h1_metSignificance->Write();


  h1_rhoPF_presel->Write();
  h1_rhoPF->Write();


  h1_ptLeptZ1->Write();
  h1_ptLeptZ2->Write();
  h1_etaLeptZ1->Write();
  h1_etaLeptZ2->Write();


  h1_deltaRZll->Write();

  h1_ptZll->Write();
  h1_etaZll->Write();
  h1_mZll->Write();



  h1_nJets->Write();


  h1_bTagJet1->Write();
  h1_bTagJet2->Write();


  h1_ptJet1->Write();
  h1_ptJet2->Write();
  h1_ptJet3->Write();
  h1_ptJet4->Write();

  h1_etaJet1->Write();
  h1_etaJet2->Write();
  h1_etaJet3->Write();
  h1_etaJet4->Write();

  h1_partFlavorJet1->Write();
  h1_partFlavorJet2->Write();
  h1_partFlavorJet3->Write();
  h1_partFlavorJet4->Write();


  h1_deltaRbb->Write();
  
  h1_mb1jj->Write();
  h1_mb2jj->Write();
  h1_mbjj_best->Write();
  
  h1_mb1jjZ->Write();
  h1_mb2jjZ->Write();
  h1_mbjjZ_best->Write();



  outFile_->Close();


} // finalize()



void Ntp1Finalizer_TTZDilepton::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="presel" ) {

    ptLeptZ1_thresh_ = 10.;
    ptLeptZ2_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    ptJet1_thresh_ = 30.;
    ptJet2_thresh_ = 30.;
    ptJet3_thresh_ = 30.;
    ptJet4_thresh_ = 30.;
    etaJet1_thresh_ = 2.4;
    etaJet2_thresh_ = 2.4;
    etaJet3_thresh_ = 2.4;
    etaJet4_thresh_ = 2.4;
    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType



