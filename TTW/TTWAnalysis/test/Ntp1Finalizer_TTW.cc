#include "Ntp1Finalizer_TTW.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TGraphAsymmErrors.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"

#include "PUWeight.h"




bool USE_MC_MASS=false;

int DEBUG_EVENTNUMBER = 675033372;


std::vector<int> getJets( const std::string& choice, std::vector< AnalysisJet > jets );
std::vector<int> getSingleJetPair( const std::string& choice, std::vector<AnalysisJet> jets, std::vector<int> *vetoIndexes=0 );
bool matchedToVeto( int iJet, std::vector<int> *vetoIndexes );


bool isMatchedPair( AnalysisJet *jet1, AnalysisJet *jet2 );
bool isMatchedToPart( AnalysisJet jet );
bool isSignalJet( AnalysisJet jet );



// constructor:

Ntp1Finalizer_TTW::Ntp1Finalizer_TTW( const std::string& dataset, const std::string& selectionType, const std::string& jetChoice, const std::string& bTaggerType ) : Ntp1Finalizer( "TTW", dataset ) {

  if( bTaggerType!="SSVHE" && bTaggerType!="TCHE"&& bTaggerType!="JP" ) {
    std::cout << "b-Tagger type '" << bTaggerType << "' currently not supported. Exiting." << std::endl;
    exit(9179);
  }

  jetChoice_ = jetChoice;

  bTaggerType_ = bTaggerType;

  setSelectionType(selectionType);

}




void Ntp1Finalizer_TTW::finalize() {

  //if( outFile_==0 ) this->createOutputFile();
  
  Int_t run;
  tree_->SetBranchAddress("run", &run);
  tree_->GetEntry(0);
  bool isMC = (run < 160000);
  std::string fullFlags = selectionType_;
  fullFlags = fullFlags + "_" + jetChoice_;
  if( bTaggerType_!="TCHE" ) fullFlags = fullFlags + "_" + bTaggerType_;
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

  TH1D* h1_rhoPF_presel = new TH1D("rhoPF_presel", "", 50, 0., 20.);
  h1_rhoPF_presel->Sumw2();
  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., 20.);
  h1_rhoPF->Sumw2();


  TH1D* h1_leptType = new TH1D("leptType", "", 3, -0.5, 2.5 );
  h1_leptType->Sumw2();
  //h1_leptType->GetXaxis()->SetLabelSize(0.1);
  h1_leptType->GetXaxis()->SetBinLabel(1, "#mu#mu");
  h1_leptType->GetXaxis()->SetBinLabel(2, "ee");
  h1_leptType->GetXaxis()->SetBinLabel(3, "e#mu");

  TH1D* h1_pfMet = new TH1D("pfMet", "", 500, 0., 500.);
  h1_pfMet->Sumw2();

  TH1D* h1_nJets = new TH1D("nJets", "", 11, -0.5, 10.5);
  h1_nJets->Sumw2();
  TH1D* h1_nBJets = new TH1D("nBJets", "", 11, -0.5, 10.5);
  h1_nBJets->Sumw2();
  TH1D* h1_nSignalJets = new TH1D("nSignalJets", "", 11, -0.5, 10.5);
  h1_nSignalJets->Sumw2();
  

  TH1D* h1_mjj0 = new TH1D("mjj0", "", 100, 0., 1000.);
  h1_mjj0->Sumw2();
  TH1D* h1_mjj1 = new TH1D("mjj1", "", 100, 0., 1000.);
  h1_mjj1->Sumw2();

  TH1D* h1_mjjSum = new TH1D("mjjSum", "", 100, 0., 1000.);
  h1_mjjSum->Sumw2();
  TH1D* h1_mjjDiff = new TH1D("mjjDiff", "", 100, 0., 1000.);
  h1_mjjDiff->Sumw2();

  

  TH1D* h1_ptLept1 = new TH1D("ptLept1", "", 500, 20., 520.);
  h1_ptLept1->Sumw2();
  TH1D* h1_ptLept2 = new TH1D("ptLept2", "", 200, 20., 220.);
  h1_ptLept2->Sumw2();
  TH1D* h1_etaLept1 = new TH1D("etaLept1", "", 50, -2.5, 2.5);
  h1_etaLept1->Sumw2();
  TH1D* h1_etaLept2 = new TH1D("etaLept2", "", 50, -2.5, 2.5);
  h1_etaLept2->Sumw2();
  TH1D* h1_combinedIsoRelLept1 = new TH1D("combinedIsoRelLept1", "", 30, 0., 0.15);
  h1_combinedIsoRelLept1->Sumw2();
  TH1D* h1_combinedIsoRelLept2 = new TH1D("combinedIsoRelLept2", "", 30, 0., 0.15);
  h1_combinedIsoRelLept2->Sumw2();

  TH1D* h1_passedCombinedIsoRel015 = new TH1D("passedCombinedIsoRel015", "", 20, 0., 20.);
  TH1D* h1_passedCombinedIsoRel005 = new TH1D("passedCombinedIsoRel005", "", 20, 0., 20.);

  TH1D* h1_etaMu = new TH1D("etaMu", "", 50, -2.5, 2.5);
  h1_etaMu->Sumw2();
  TH1D* h1_etaEle = new TH1D("etaEle", "", 50, -2.5, 2.5);
  h1_etaEle->Sumw2();


  TH1D* h1_deltaRll = new TH1D("deltaRll", "", 500, 0., 5.);
  h1_deltaRll->Sumw2();

  TH1D* h1_mT1 = new TH1D("mT1", "", 200, 0., 200.);
  h1_mT1->Sumw2();
  TH1D* h1_mT2 = new TH1D("mT2", "", 200, 0., 200.);
  h1_mT2->Sumw2();
  TH1D* h1_mTMax = new TH1D("mTMax", "", 400, 0., 400.);
  h1_mTMax->Sumw2();



  TH1D* h1_bTagJet1 = new TH1D("bTagJet1", "", 420, -1., 20.);
  h1_bTagJet1->Sumw2();
  TH1D* h1_bTagJet2 = new TH1D("bTagJet2", "", 420, -1., 20.);
  h1_bTagJet2->Sumw2();

  TH1D* h1_kinfit_chiSquare = new TH1D("kinfit_chiSquare", "", 60, 0., 10.);
  h1_kinfit_chiSquare->Sumw2();
  TH1D* h1_kinfit_chiSquareProb = new TH1D("kinfit_chiSquareProb", "", 60, 0., 1.0001);
  h1_kinfit_chiSquareProb->Sumw2();


  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 400, 20., 420.);
  h1_ptJet1->Sumw2();
  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 400, 20., 420.);
  h1_ptJet2->Sumw2();
  TH1D* h1_ptJet3 = new TH1D("ptJet3", "", 400, 20., 420.);
  h1_ptJet3->Sumw2();
  TH1D* h1_ptJet4 = new TH1D("ptJet4", "", 400, 20., 420.);
  h1_ptJet4->Sumw2();

  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 200, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 200, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_etaJet3 = new TH1D("etaJet3", "", 200, -5., 5.);
  h1_etaJet3->Sumw2();
  TH1D* h1_etaJet4 = new TH1D("etaJet4", "", 200, -5., 5.);
  h1_etaJet4->Sumw2();

  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 50, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 50, 0., 1.0001);
  h1_QGLikelihoodJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet3 = new TH1D("QGLikelihoodJet3", "", 50, 0., 1.0001);
  h1_QGLikelihoodJet3->Sumw2();
  TH1D* h1_QGLikelihoodJet4 = new TH1D("QGLikelihoodJet4", "", 50, 0., 1.0001);
  h1_QGLikelihoodJet4->Sumw2();

  TH1D* h1_partFlavorJet1 = new TH1D("partFlavorJet1", "", 38, -15.5, 22.5);
  h1_partFlavorJet1->Sumw2();
  TH1D* h1_partFlavorJet2 = new TH1D("partFlavorJet2", "", 38, -15.5, 22.5);
  h1_partFlavorJet2->Sumw2();
  TH1D* h1_partFlavorJet3 = new TH1D("partFlavorJet3", "", 38, -15.5, 22.5);
  h1_partFlavorJet3->Sumw2();
  TH1D* h1_partFlavorJet4 = new TH1D("partFlavorJet4", "", 38, -15.5, 22.5);
  h1_partFlavorJet4->Sumw2();


  TH1D* h1_deltaRbb = new TH1D("deltaRbb", "", 200, 0., 5.);
  h1_deltaRbb->Sumw2();
  TH1D* h1_deltaRqq = new TH1D("deltaRqq", "", 200, 0., 5.);
  h1_deltaRqq->Sumw2();
  
  TH1D* h1_deltaR_b_lept_max = new TH1D("deltaR_b_lept_max", "", 200, 0., 5.);
  h1_deltaR_b_lept_max->Sumw2();
  TH1D* h1_deltaR_b_lept_min = new TH1D("deltaR_b_lept_min", "", 200, 0., 5.);
  h1_deltaR_b_lept_min->Sumw2();
  
  TH1D* h1_mbb = new TH1D("mbb", "", 500, 0., 1000.);
  h1_mbb->Sumw2();
  TH1D* h1_mqq = new TH1D("mqq", "", 500, 0., 1000.);
  h1_mqq->Sumw2();
  
  TH1D* h1_ptbb = new TH1D("ptbb", "", 500, 0., 1000.);
  h1_ptbb->Sumw2();
  TH1D* h1_ptqq = new TH1D("ptqq", "", 500, 0., 1000.);
  h1_ptqq->Sumw2();
  



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

  bool is_ttHWWEvent;
  tree_->SetBranchAddress("is_ttHWWEvent", &is_ttHWWEvent);
  bool is_ttHBBEvent;
  tree_->SetBranchAddress("is_ttHBBEvent", &is_ttHBBEvent);
  bool is_ttHZZEvent;
  tree_->SetBranchAddress("is_ttHZZEvent", &is_ttHZZEvent);
  bool is_ttHTauTauEvent;
  tree_->SetBranchAddress("is_ttHTauTauEvent", &is_ttHTauTauEvent);
  bool is_ttHWW4QEvent;
  tree_->SetBranchAddress("is_ttHWW4QEvent", &is_ttHWW4QEvent);
  bool is_ttHWW4QEvent_reco;
  tree_->SetBranchAddress("is_ttHWW4QEvent_reco", &is_ttHWW4QEvent_reco);
  bool is_ssdlEvent;
  tree_->SetBranchAddress("is_ssdlEvent", &is_ssdlEvent);
  bool is_ssdlEvent_reco;
  tree_->SetBranchAddress("is_ssdlEvent_reco", &is_ssdlEvent_reco);
  bool is_ssdl_ttHWW4QEvent;
  tree_->SetBranchAddress("is_ssdl_ttHWW4QEvent", &is_ssdl_ttHWW4QEvent);
  bool is_ssdl_ttHWW4QEvent_reco;
  tree_->SetBranchAddress("is_ssdl_ttHWW4QEvent_reco", &is_ssdl_ttHWW4QEvent);

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


  Float_t eLept1;
  tree_->SetBranchAddress("eLept1", &eLept1);
  Float_t ptLept1;
  tree_->SetBranchAddress("ptLept1", &ptLept1);
  Float_t etaLept1;
  tree_->SetBranchAddress("etaLept1", &etaLept1);
  Float_t phiLept1;
  tree_->SetBranchAddress("phiLept1", &phiLept1);
  Int_t chargeLept1;
  tree_->SetBranchAddress("chargeLept1", &chargeLept1);
  Int_t leptTypeLept1;
  tree_->SetBranchAddress("leptTypeLept1", &leptTypeLept1);
  Float_t combinedIsoRelLept1;
  tree_->SetBranchAddress("combinedIsoRelLept1", &combinedIsoRelLept1);

  Float_t eLept2;
  tree_->SetBranchAddress("eLept2", &eLept2);
  Float_t ptLept2;
  tree_->SetBranchAddress("ptLept2", &ptLept2);
  Float_t etaLept2;
  tree_->SetBranchAddress("etaLept2", &etaLept2);
  Float_t phiLept2;
  tree_->SetBranchAddress("phiLept2", &phiLept2);
  Int_t chargeLept2;
  tree_->SetBranchAddress("chargeLept2", &chargeLept2);
  Int_t leptTypeLept2;
  tree_->SetBranchAddress("leptTypeLept2", &leptTypeLept2);
  Float_t combinedIsoRelLept2;
  tree_->SetBranchAddress("combinedIsoRelLept2", &combinedIsoRelLept2);


  //Int_t nLept;
  //tree_->SetBranchAddress("nLept", &nLept);
  //Int_t leptTypeLept[10];
  //tree_->SetBranchAddress("leptTypeLept", leptTypeLept);
  //Float_t eLept[10];
  //tree_->SetBranchAddress("eLept", eLept);
  //Float_t ptLept[10];
  //tree_->SetBranchAddress("ptLept", ptLept);
  //Float_t etaLept[10];
  //tree_->SetBranchAddress("etaLept", etaLept);
  //Float_t phiLept[10];
  //tree_->SetBranchAddress("phiLept", phiLept);
  //Int_t chargeLept[10];
  //tree_->SetBranchAddress("chargeLept", chargeLept);


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
  Float_t eGenJet[50];
  tree_->SetBranchAddress("eGenJet", eGenJet);
  Float_t ptGenJet[50];
  tree_->SetBranchAddress("ptGenJet", ptGenJet);
  Float_t etaGenJet[50];
  tree_->SetBranchAddress("etaGenJet", etaGenJet);
  Float_t phiGenJet[50];
  tree_->SetBranchAddress("phiGenJet", phiGenJet);
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
  Float_t QGLikelihoodJet[50];
  tree_->SetBranchAddress("QGLikelihoodJet", QGLikelihoodJet);
  Float_t ePartJet[50];
  tree_->SetBranchAddress("ePartJet", ePartJet);
  Float_t ptPartJet[50];
  tree_->SetBranchAddress("ptPartJet", ptPartJet);
  Float_t etaPartJet[50];
  tree_->SetBranchAddress("etaPartJet", etaPartJet);
  Float_t phiPartJet[50];
  tree_->SetBranchAddress("phiPartJet", phiPartJet);
  Int_t pdgIdPartJet[50];
  tree_->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);
  Int_t pdgIdPartMomJet[50];
  tree_->SetBranchAddress("pdgIdPartMomJet", pdgIdPartMomJet);
  Int_t pdgIdPartMomMomJet[50];
  tree_->SetBranchAddress("pdgIdPartMomMomJet", pdgIdPartMomMomJet);



  Int_t nPart;
  tree_->SetBranchAddress("nPart", &nPart);
  Float_t ePart[100];
  tree_->SetBranchAddress("ePart", ePart);
  Float_t ptPart[100];
  tree_->SetBranchAddress("ptPart", ptPart);
  Float_t etaPart[100];
  tree_->SetBranchAddress("etaPart", etaPart);
  Float_t phiPart[100];
  tree_->SetBranchAddress("phiPart", phiPart);
  Int_t pdgIdPart[100];
  tree_->SetBranchAddress("pdgIdPart", pdgIdPart);
  Int_t statusPart[100];
  tree_->SetBranchAddress("statusPart", statusPart);


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
  float tmass = 172.9;
  float Wmass = 80.399;


//std::string puType = "Spring11_Flat10";
//std::string puType_ave = "Spring11_Flat10";
//TString dataset_tstr(dataset_);
//if( dataset_tstr.Contains("Summer11") && dataset_tstr.Contains("PU_S4") ) {
//  puType = "Summer11_S4";
//  puType_ave = "Summer11_S4_ave";
//} else if( dataset_tstr.Contains("Fall11") ) {
//  puType = "Fall11";
//}
//PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
//PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
//std::string puFileName;
//puFileName = "all2011AB.pileup_v2_73mb.root";
//std::cout << std::endl << "-> Using data pileup file: " << puFileName << std::endl;
//TFile* filePU = TFile::Open(puFileName.c_str());
//TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
//fPUWeight->SetDataHistogram(h1_nPU_data);
//fPUWeight_ave->SetDataHistogram(h1_nPU_data);
    
//TFile* filePUMC = TFile::Open("Pileup_MC_Summer11_S4.root");
//TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("hNPU");
//std::cout << "-> Switching MC PU file to: Pileup_MC_Summer11_S4.root" << std::endl;
//fPUWeight->SetMCHistogram(h1_nPU_mc);



  int leptType;
  float HLTSF;
  int nbjets, nbjetsmed, nnbjets, nnbjetsmed;
  float ptLept1_t, ptLept2_t, etaLept1_t, etaLept2_t;
  float ptJet1_t, ptJet2_t, ptJet3_t, ptJet4_t;
  float ptbjet_t[20], etabjet_t[20];
  float etaJet1_t, etaJet2_t, etaJet3_t, etaJet4_t;
  float QGLikelihoodJet1_t, QGLikelihoodJet2_t, QGLikelihoodJet3_t, QGLikelihoodJet4_t;
  float bTagJet1_t, bTagJet2_t;
  float mll, deltaRll;
  float mjj1_t, mjj2_t;
  float ht;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "HLTSF", &HLTSF, "HLTSF/F" );

  tree_passedEvents->Branch( "leptType", &leptType, "leptType/I" );
  tree_passedEvents->Branch( "ptLept1", &ptLept1_t, "ptLept1_t/F" );
  tree_passedEvents->Branch( "ptLept2", &ptLept2_t, "ptLept2_t/F" );
  tree_passedEvents->Branch( "etaLept1", &etaLept1_t, "etaLept1_t/F" );
  tree_passedEvents->Branch( "etaLept2", &etaLept2_t, "etaLept2_t/F" );
  tree_passedEvents->Branch( "mll", &mll, "mll/F" );
  tree_passedEvents->Branch( "deltaRll", &deltaRll, "deltaRll/F" );
  tree_passedEvents->Branch( "pfMet", &pfMet, "pfMet/F" );
  tree_passedEvents->Branch( "ht", &ht, "ht/F" );
  tree_passedEvents->Branch( "nbjetsmed", &nbjetsmed, "nbjetsmed/I" );
  tree_passedEvents->Branch( "nbjets", &nbjets, "nbjets/I" );
  tree_passedEvents->Branch( "ptbjet", ptbjet_t, "ptbjet_t[nbjets]/F" );
  tree_passedEvents->Branch( "etabjet", etabjet_t, "etabjet_t[nbjets]/F" );

  tree_passedEvents->Branch( "nnbjets", &nnbjets, "nnbjets/I" ); // non btagged jets
  tree_passedEvents->Branch( "nnbjetsmed", &nnbjetsmed, "nnbjetsmed/I" ); // non btagged jets

  tree_passedEvents->Branch( "ptJet1", &ptJet1_t, "ptJet1_t/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2_t, "ptJet2_t/F" );
  tree_passedEvents->Branch( "ptJet3", &ptJet3_t, "ptJet3_t/F" );
  tree_passedEvents->Branch( "ptJet4", &ptJet4_t, "ptJet4_t/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1_t, "etaJet1_t/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2_t, "etaJet2_t/F" );
  tree_passedEvents->Branch( "etaJet3", &etaJet3_t, "etaJet3_t/F" );
  tree_passedEvents->Branch( "etaJet4", &etaJet4_t, "etaJet4_t/F" );

  tree_passedEvents->Branch( "mjj1", &mjj1_t, "mjj1_t/F" );
  tree_passedEvents->Branch( "mjj2", &mjj2_t, "mjj2_t/F" );

  tree_passedEvents->Branch( "is_ttHWWEvent", &is_ttHWWEvent, "is_ttHWWEvent/O" );
  tree_passedEvents->Branch( "is_ttHBBEvent", &is_ttHBBEvent, "is_ttHBBEvent/O" );
  tree_passedEvents->Branch( "is_ttHTauTauEvent", &is_ttHTauTauEvent, "is_ttHTauTauEvent/O" );
  tree_passedEvents->Branch( "is_ttHZZEvent", &is_ttHZZEvent, "is_ttHZZEvent/O" );
  tree_passedEvents->Branch( "is_ttHWW4QEvent_reco", &is_ttHWW4QEvent_reco, "is_ttHWW4QEvent_reco/O" );




ofstream ofs("run_event.txt");




  std::cout << std::endl << std::endl;
  std::cout << "+++ BEGINNING ANALYSIS LOOP" << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  std::cout << "----> JET CHOICE: " << jetChoice_ << std::endl;
  std::cout << "----> B-TAGGER: " << bTaggerType_ << std::endl;
  std::cout << std::endl << std::endl;


  int nPairsTot0 = 0;
  int nPairsTot1 = 0;
  int foundCorrectJetPair0 = 0;
  int foundCorrectJetPair1 = 0;

  int nPairsTot0_matched = 0;
  int nPairsTot1_matched = 0;
  int foundCorrectJetPair0_matched = 0;
  int foundCorrectJetPair1_matched = 0;


  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 20000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;



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



//  // this is dilepton channel: no other lepton in the event 
//  if( nLept!=0 ) continue;



    h1_rhoPF_presel->Fill( rhoPF, eventWeight);


    if( pfMet < pfMet_thresh_ ) continue;

    if( chargeLept1!=chargeLept2 ) continue;


    TLorentzVector Lept1, Lept2;
    Lept1.SetPtEtaPhiE( ptLept1, etaLept1, phiLept1, eLept1 );
    Lept2.SetPtEtaPhiE( ptLept2, etaLept2, phiLept2, eLept2 );


    TLorentzVector diLepton = Lept1+Lept2;

    if( leptTypeLept1==0 && leptTypeLept2==0 ) 
      leptType=0;
    else if( leptTypeLept1==1 && leptTypeLept2==1 )
      leptType=1;
    else 
      leptType=2;
  







    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "** LOG FOR RUN: " << run << "   EVENT: " << DEBUG_EVENTNUMBER << std::endl << std::endl;
      std::cout << "leptType: " << leptType << std::endl; 
      std::cout << "Lept1.Pt(): " << Lept1.Pt() << " Lept1.Eta(): " << Lept1.Eta() << std::endl;
      std::cout << "Lept2.Pt(): " << Lept2.Pt() << " Lept2.Eta(): " << Lept2.Eta() << std::endl;
      std::cout << "combinedIsoRelLept1: " << combinedIsoRelLept1 << " combinedIsoRelLept2: " << combinedIsoRelLept2 << std::endl;
      std::cout << "diLepton.M(): " << diLepton.M() << std::endl;
    }


    // ----------------------------
    // KINEMATIC SELECTION: LEPTONS
    // ----------------------------

    mll = diLepton.M();

    if( Lept1.Pt() < ptLept1_thresh_ ) continue;
    if( Lept2.Pt() < ptLept2_thresh_ ) continue;
    if( fabs(Lept1.Eta()) > etaLept1_thresh_ ) continue;
    if( fabs(Lept2.Eta()) > etaLept2_thresh_ ) continue;
    if( combinedIsoRelLept1 > 0.05 ) continue;
    if( combinedIsoRelLept2 > 0.05 ) continue;
    if( mll < mll_thresh_ ) continue;


    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "-> LEPTONS PASSED CUTS" << std::endl;
    }

    
    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "-> TOTAL NUMBER OF JETS IN EVENT: " << nJets << std::endl;
    }
    
    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "-> FIRST: LOOK FOR BEST B-TAGGED ONES" << std::endl;
    }
    

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



    std::vector<AnalysisJet> btaggedJets;
    std::vector<AnalysisJet> nonBtaggedJets;

    std::vector<AnalysisJet> btaggedMedJets;
    std::vector<AnalysisJet> nonBtaggedMedJets;

    int nSignalJets=0;
    ht = 0.;

    for( unsigned iJet=0; iJet<nJets; ++iJet) {

      AnalysisJet thisJet;
      thisJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      thisJet.rmsCand = rmsCandJet[iJet];
      thisJet.ptD = ptDJet[iJet];
      thisJet.nCharged = nChargedJet[iJet];
      thisJet.nNeutral = nNeutralJet[iJet];
      thisJet.eMuons= eMuonsJet[iJet]/thisJet.Energy();
      thisJet.eElectrons = eElectronsJet[iJet]/thisJet.Energy();

      thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet[iJet];
      thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet[iJet];
      thisJet.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet[iJet];
      thisJet.jetProbabilityBJetTag               = jetProbabilityBJetTagJet[iJet];

      thisJet.QGLikelihood               = QGLikelihoodJet[iJet];

      thisJet.ptGen = ptGenJet[iJet];
      thisJet.etaGen = etaGenJet[iJet];
      thisJet.phiGen = phiGenJet[iJet];
      thisJet.eGen = eGenJet[iJet];

      thisJet.ptPart = ptPartJet[iJet];
      thisJet.etaPart = etaPartJet[iJet];
      thisJet.phiPart = phiPartJet[iJet];
      thisJet.ePart = ePartJet[iJet];
      thisJet.pdgIdPart = pdgIdPartJet[iJet];
      thisJet.pdgIdPartMom = pdgIdPartMomJet[iJet];
      thisJet.pdgIdPartMomMom = pdgIdPartMomMomJet[iJet];



      if( isSignalJet(thisJet) ) nSignalJets++;

      if( event==DEBUG_EVENTNUMBER ) {
        std::cout << std::endl << "Jet N." << iJet << " Pt: " << ptJet[iJet] << " \tEta: " << etaJet[iJet];
      }

      if( thisJet.Pt()<ptJet_thresh_ ) {
        if( event==DEBUG_EVENTNUMBER ) std::cout << "\t (didn't pass pt cut)";
        continue;
      }
      if( fabs(thisJet.Eta())>etaJet_thresh_ ) {
        if( event==DEBUG_EVENTNUMBER ) std::cout << "\t (didn't pass eta cut)";
        continue;
      }

      ht += thisJet.Pt();



      if( thisJet.jetProbabilityBJetTag > 0.275 ) { // is btagged loose (JP)
        btaggedJets.push_back(thisJet);
      } else {
        nonBtaggedJets.push_back(thisJet);
      }

      if( thisJet.jetProbabilityBJetTag > 0.545 ) { // is btagged medium (JP)
        btaggedMedJets.push_back(thisJet);
      } else {
        nonBtaggedMedJets.push_back(thisJet);
      }



      ////match to parton:
      //int partFlavor=0;
      //float deltaRmin=999.;
      //for(unsigned iPart=0; iPart<nPart; ++iPart ) {
      //  TLorentzVector thisPart;
      //  if( !( fabs(pdgIdPart[iPart])<7 || pdgIdPart[iPart]==21) ) continue;
      //  thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
      //  float thisDeltaR = thisJet.DeltaR(thisPart);
      //  if( thisDeltaR<deltaRmin ) {
      //    partFlavor = pdgIdPart[iPart];
      //    deltaRmin = thisDeltaR;
      //  }
      //}
      //thisJet.pdgIdPart = partFlavor;


      //float thisBtag;
      //if( bTaggerType_=="TCHE" )
      //  thisBtag = thisJet.trackCountingHighEffBJetTag;
      //else if( bTaggerType_=="SSVHE" ) 
      //  thisBtag = thisJet.simpleSecondaryVertexHighEffBJetTag;

      //if( event==DEBUG_EVENTNUMBER ) {
      //  std::cout << " \tbtag: " << thisBtag;
      //}


      //if( thisBtag > bestBtag ) {

      //  if( event==DEBUG_EVENTNUMBER ) {
      //    std::cout << "\t <--- THIS IS NEW LEADING BTAG JET IN EVENT";
      //  }
      //  
      //  bestBtag2 = bestBtag;
      //  bestBtag = thisBtag;
      //  i_bestBtag2 = i_bestBtag;
      //  i_bestBtag = iJet;

      //  if( jets.size()==0 ) {
      //    AnalysisJet* newJet = new AnalysisJet(thisJet);
      //    jets.push_back(*newJet);
      //  } else if( jets.size()==1 ) {
      //    AnalysisJet oldJet = jets[0];
      //    AnalysisJet* newJet = new AnalysisJet(thisJet);
      //    jets.push_back(*newJet);
      //    jets[0] = *newJet;
      //    jets[1] = oldJet;
      //  } else if( jets.size()==2 ) {
      //    jets[1] = jets[0];
      //    jets[0] = thisJet;
      //  }

      //} else if( thisBtag > bestBtag2 ) { //means that at least one jet was found

      //  bestBtag2 = thisBtag;
      //  i_bestBtag2 = iJet;

      //  if( event==DEBUG_EVENTNUMBER ) {
      //    std::cout << "\t <--- THIS IS NEW SUBLEADING BTAG JET IN EVENT";
      //  }
      //  
      //  if( jets.size()==1 ) {
      //    AnalysisJet* newJet = new AnalysisJet(thisJet);
      //    jets.push_back(*newJet);
      //  } else if( jets.size()==2 ) {
      //    jets[1] = thisJet;
      //  }
  
      //}
      
    } // for jets


    if( is_ttHWWEvent ) 
      h1_nSignalJets->Fill( nSignalJets, eventWeight );

    
    if( event==DEBUG_EVENTNUMBER ) 
      std::cout << std::endl << std::endl << std::endl << "-> Event HT: " << ht << "   (thresh is " << ht_thresh_ << ")" << std::endl;


    nbjetsmed = btaggedMedJets.size();
    nnbjetsmed = nonBtaggedMedJets.size();

    nbjets = btaggedJets.size();
    nnbjets = nonBtaggedJets.size();

    int njets  = nbjets + nnbjets;

    if( nbjets < nbjets_thresh_ ) continue;
    if( njets  < njets_thresh_ ) continue;
    
    if( ht < ht_thresh_ ) continue;


    ptJet1_t = -1.;
    ptJet2_t = -1.;
    ptJet3_t = -1.;
    ptJet4_t = -1.;

    etaJet1_t = -20.;
    etaJet2_t = -20.;
    etaJet3_t = -20.;
    etaJet4_t = -20.;

    mjj1_t = -1.;
    mjj2_t = -1.;


    // if at least 4 jets, try to find the two hadronic W's
    if( nonBtaggedJets.size() >= 2 ) {

      //std::vector<int> indexJets = getJets( jetChoice_, nonBtaggedJets );

      std::vector<int> firstPair  = getSingleJetPair( jetChoice_, nonBtaggedJets );

      AnalysisJet dijetpair0_0 = nonBtaggedJets[firstPair[0]];
      AnalysisJet dijetpair0_1 = nonBtaggedJets[firstPair[1]];

      bool correctJetPair0 = isMatchedPair( &dijetpair0_0, &dijetpair0_1 );

      TLorentzVector dijet0 = dijetpair0_0 + dijetpair0_1;

      ptJet1_t = dijetpair0_0.Pt();
      ptJet2_t = dijetpair0_1.Pt();

      etaJet1_t = dijetpair0_0.Eta();
      etaJet2_t = dijetpair0_1.Eta();

      mjj1_t = dijet0.M();


      if( is_ttHWWEvent ) nPairsTot0_matched += 1;
      if( is_ttHWWEvent && correctJetPair0 ) foundCorrectJetPair0_matched += 1;

      nPairsTot0 += 1;
      if( correctJetPair0 ) foundCorrectJetPair0 += 1;

      // fill histos
      h1_mjj0->Fill( dijet0.M(), eventWeight );

      if( nonBtaggedJets.size() >= 4 ) {

        std::vector<int> secondPair = getSingleJetPair( jetChoice_, nonBtaggedJets, &firstPair );

        AnalysisJet dijetpair1_0 = nonBtaggedJets[secondPair[0]];
        AnalysisJet dijetpair1_1 = nonBtaggedJets[secondPair[1]];

        bool correctJetPair1 = isMatchedPair( &dijetpair1_0, &dijetpair1_1 );

        //if( is_ssdl_ttHWW4QEvent_reco ) nPairsTot_matched += 1;
        //if( is_ssdl_ttHWW4QEvent_reco && correctJetPair0 ) foundCorrectJetPair0_matched += 1;
        //if( is_ssdl_ttHWW4QEvent_reco && correctJetPair1 ) foundCorrectJetPair1_matched += 1;

        if( is_ttHWWEvent ) nPairsTot1_matched += 1;
        if( is_ttHWWEvent && correctJetPair1 ) foundCorrectJetPair1_matched += 1;

        nPairsTot1 += 1;
        if( correctJetPair1 ) foundCorrectJetPair1 += 1;

        TLorentzVector dijet1 = dijetpair1_0 + dijetpair1_1;

        ptJet3_t = dijetpair1_0.Pt();
        ptJet4_t = dijetpair1_1.Pt();

        etaJet3_t = dijetpair1_0.Eta();
        etaJet4_t = dijetpair1_1.Eta();

        mjj2_t = dijet1.M();


        // fill histos
        h1_mjj1->Fill( dijet1.M(), eventWeight );
  
        h1_mjjDiff->Fill( dijet0.M()-dijet1.M(), eventWeight );
        h1_mjjSum->Fill( dijet0.M()+dijet1.M(), eventWeight );

      } // if 4 jets

    }  // if 2 jets


    h1_nJets->Fill( njets, eventWeight );
    h1_nBJets->Fill( nbjets, eventWeight );


 // if( event==DEBUG_EVENTNUMBER ) {
 //   std::cout << std::endl << std::endl << std::endl << "-> So here are two best-btag jets: " << std::endl;
 //   std::cout << "[1]  Pt: " << jets[0].Pt() << " \tEta: " << jets[0].Eta() << " btag: " << bestBtag << std::endl;
 //   std::cout << "[2]  Pt: " << jets[1].Pt() << " \tEta: " << jets[1].Eta() << " btag: " << bestBtag2 << std::endl;
 //   std::cout << std::endl << "NOW ADDING OTHER JETS (ORDERED IN PT)" << std::endl;
 // }





//  if( event==DEBUG_EVENTNUMBER ) {
//    std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
//    std::cout << "-> JETS PASSING REQUIREMENTS: " << std::endl;
//    for( unsigned ijet=0; ijet<jets.size(); ++ijet ) 
//      std::cout << ijet << "  Pt: " << jets[ijet].Pt() << " \tEta: " << jets[ijet].Eta() << std::endl;
//  }
    
 
        


    // -------------------------
    // KINEMATIC SELECTION: JETS
    // -------------------------
  
    //if( jet1.Pt() < ptJet1_thresh_ ) continue;
    //if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet pt OK" << std::endl;
    //if( jet2.Pt() < ptJet2_thresh_ ) continue;
    //if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet pt OK" << std::endl;
    //if( fabs(jet1.Eta()) > etaJet1_thresh_ ) continue;
    //if( event==DEBUG_EVENTNUMBER ) std::cout << "first jet eta OK" << std::endl;
    //if( fabs(jet2.Eta()) > etaJet2_thresh_ ) continue;
    //if( event==DEBUG_EVENTNUMBER ) std::cout << "second jet eta OK" << std::endl;




//  // -------------
//  // KINEMATIC FIT
//  // -------------
//
//  std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(jets[2], jets[3]);

//  float chiSquareProb = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());
//  h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
//  h1_kinfit_chiSquareProb->Fill( chiSquareProb, eventWeight );


//  
//  float bTaggerJet1, bTaggerJet2;
//  if( bTaggerType_=="TCHE" ) {
//    bTaggerJet1 = jets[0].trackCountingHighEffBJetTag;
//    bTaggerJet2 = jets[1].trackCountingHighEffBJetTag;
//  } else if( bTaggerType_=="SSVHE" ) {
//    bTaggerJet1 = jets[0].simpleSecondaryVertexHighEffBJetTag;
//    bTaggerJet2 = jets[1].simpleSecondaryVertexHighEffBJetTag;
//  }

//  if( event==DEBUG_EVENTNUMBER ) {
//    std::cout << std::endl << "-- FOUND JETS: " << std::endl;
//    std::cout << "jet1: \tpt:" << jets[0].Pt() << " \teta: " << jets[0].Eta() << "\tbtag: " << bTaggerJet1 << std::endl;
//    std::cout << "jet2: \tpt:" << jets[1].Pt() << " \teta: " << jets[1].Eta() << "\tbtag: " << bTaggerJet2 << std::endl;
//    std::cout << "jet3: \tpt:" << jets[2].Pt() << " \teta: " << jets[2].Eta() << std::endl;
//    std::cout << "jet4: \tpt:" << jets[3].Pt() << " \teta: " << jets[3].Eta() << std::endl;
//  }


//  h1_bTagJet1->Fill( bTaggerJet1, eventWeight );
//  h1_bTagJet2->Fill( bTaggerJet2, eventWeight );

//  if( !(this->passedBTag( bTaggerJet1, bTaggerJet2, bTaggerType_ ) ) ) continue;
//  //if( bTaggerJet1<btag_threshold && bTaggerJet2<btag_threshold ) continue;


//  h1_leptType->Fill( leptType, eventWeight );

//  h1_ptLept1->Fill( Lept1.Pt(), eventWeight );
//  h1_ptLept2->Fill( Lept2.Pt(), eventWeight );
//  h1_etaLept1->Fill( Lept1.Eta(), eventWeight );
//  h1_etaLept2->Fill( Lept2.Eta(), eventWeight );
//  if( leptTypeLept1==0 ) {
//    h1_etaMu->Fill( Lept1.Eta(), eventWeight );
//  } else {
//    h1_etaEle->Fill( Lept1.Eta(), eventWeight );
//  }
//  if( leptTypeLept2==0 ) {
//    h1_etaMu->Fill( Lept2.Eta(), eventWeight );
//  } else {
//    h1_etaEle->Fill( Lept2.Eta(), eventWeight );
//  }
//  h1_combinedIsoRelLept1->Fill( combinedIsoRelLept1, eventWeight );
//  h1_combinedIsoRelLept2->Fill( combinedIsoRelLept2, eventWeight );

//  // these ones with no weight (for eff vs. PU):
//  h1_passedCombinedIsoRel015->Fill( rhoPF );
//  if( combinedIsoRelLept1<0.05 && combinedIsoRelLept2<0.05 ) 
//    h1_passedCombinedIsoRel005->Fill( rhoPF );

//  deltaRll = Lept1.DeltaR(Lept2);
//  h1_deltaRll->Fill( deltaRll, eventWeight );


//  h1_nJets->Fill( jets.size(), eventWeight );
//  h1_pfMet->Fill( pfMet, eventWeight );
//  h1_metSignificance->Fill( metSignificance, eventWeight );

//  TLorentzVector v_met;
//  v_met.SetPtEtaPhiE( pfMet, 0., phiMet, pfMet );

//  float deltaPhiLept1_met = Lept1.DeltaPhi(v_met);
//  float deltaPhiLept2_met = Lept2.DeltaPhi(v_met);

//  float mT1 = sqrt( Lept1.Pt()*pfMet*(1.-cos(deltaPhiLept1_met)) );
//  float mT2 = sqrt( Lept2.Pt()*pfMet*(1.-cos(deltaPhiLept2_met)) );

//  h1_mT1->Fill( mT1, eventWeight );
//  h1_mT2->Fill( mT2, eventWeight );
//  h1_mTMax->Fill( TMath::Max(mT1, mT2), eventWeight );

//  h1_ptJet1->Fill( jets[0].Pt(), eventWeight );
//  h1_ptJet2->Fill( jets[1].Pt(), eventWeight );
//  h1_ptJet3->Fill( jets[2].Pt(), eventWeight );
//  h1_ptJet4->Fill( jets[3].Pt(), eventWeight );

//  h1_etaJet1->Fill( jets[0].Eta(), eventWeight );
//  h1_etaJet2->Fill( jets[1].Eta(), eventWeight );
//  h1_etaJet3->Fill( jets[2].Eta(), eventWeight );
//  h1_etaJet4->Fill( jets[3].Eta(), eventWeight );

//  h1_QGLikelihoodJet1->Fill( jets[0].QGLikelihood, eventWeight );
//  h1_QGLikelihoodJet2->Fill( jets[1].QGLikelihood, eventWeight );
//  h1_QGLikelihoodJet3->Fill( jets[2].QGLikelihood, eventWeight );
//  h1_QGLikelihoodJet4->Fill( jets[3].QGLikelihood, eventWeight );

//  h1_partFlavorJet1->Fill( jets[0].pdgIdPart, eventWeight );
//  h1_partFlavorJet2->Fill( jets[1].pdgIdPart, eventWeight );
//  h1_partFlavorJet3->Fill( jets[2].pdgIdPart, eventWeight );
//  h1_partFlavorJet4->Fill( jets[3].pdgIdPart, eventWeight );


//  TLorentzVector diJet_bb = jets[0] + jets[1];
//  TLorentzVector diJet_qq = jets[2] + jets[3];

//  h1_mbb->Fill( diJet_bb.M(), eventWeight );
//  h1_mqq->Fill( diJet_qq.M(), eventWeight );

//  h1_ptbb->Fill( diJet_bb.Pt(), eventWeight );
//  h1_ptqq->Fill( diJet_qq.Pt(), eventWeight );

//  h1_deltaRbb->Fill( jets[0].DeltaR(jets[1]), eventWeight );
//  h1_deltaRqq->Fill( jets[2].DeltaR(jets[3]), eventWeight );


//  float deltaR_b0_lept1 = jets[0].DeltaR(Lept1);
//  float deltaR_b0_lept2 = jets[0].DeltaR(Lept2);
//  float deltaR_b1_lept1 = jets[1].DeltaR(Lept1);
//  float deltaR_b1_lept2 = jets[1].DeltaR(Lept2);

//  float deltaR_b0_lept_max, deltaR_b0_lept_min;
//  if( deltaR_b0_lept1>deltaR_b0_lept2 ) {
//    deltaR_b0_lept_max = deltaR_b0_lept1;
//    deltaR_b0_lept_min = deltaR_b0_lept2;
//  } else {
//    deltaR_b0_lept_min = deltaR_b0_lept1;
//    deltaR_b0_lept_max = deltaR_b0_lept2;
//  }

//  float deltaR_b1_lept_max, deltaR_b1_lept_min;
//  if( deltaR_b1_lept1>deltaR_b1_lept2 ) {
//    deltaR_b1_lept_max = deltaR_b1_lept1;
//    deltaR_b1_lept_min = deltaR_b1_lept2;
//  } else {
//    deltaR_b1_lept_min = deltaR_b1_lept1;
//    deltaR_b1_lept_max = deltaR_b1_lept2;
//  }

//  float deltaR_b_lept_max = (deltaR_b1_lept_max>deltaR_b0_lept_max) ? deltaR_b1_lept_max : deltaR_b0_lept_max;
//  float deltaR_b_lept_min = (deltaR_b1_lept_min>deltaR_b0_lept_min) ? deltaR_b1_lept_min : deltaR_b0_lept_min;

//  h1_deltaR_b_lept_max->Fill(deltaR_b_lept_max, eventWeight);
//  h1_deltaR_b_lept_min->Fill(deltaR_b_lept_min, eventWeight);


    ptLept1_t = Lept1.Pt();
    ptLept2_t = Lept2.Pt();
    etaLept1_t = Lept1.Eta();
    etaLept2_t = Lept2.Eta();

//  QGLikelihoodJet1_t = jets[0].QGLikelihood;
//  QGLikelihoodJet2_t = jets[1].QGLikelihood;
//  QGLikelihoodJet3_t = jets[2].QGLikelihood;
//  QGLikelihoodJet4_t = jets[3].QGLikelihood;



    // and fill tree:
    tree_passedEvents->Fill();

  

  
  } //for entries



  std::cout << std::endl << std::endl;
  std::cout << "Jet choice: " << jetChoice_ << std::endl;
  std::cout << "Found correct pair 0 in " << foundCorrectJetPair0 << "/" << nPairsTot0 << " = " << (float)100.*foundCorrectJetPair0/nPairsTot0 << "% of cases." << std::endl;
  std::cout << "Found correct pair 1 in " << foundCorrectJetPair1 << "/" << nPairsTot1 << " = " << (float)100.*foundCorrectJetPair1/nPairsTot1 << "% of cases." << std::endl;
  std::cout << "**MATCHED** Found correct pair 0 in " << foundCorrectJetPair0_matched << "/" << nPairsTot0_matched << " = " << (float)100.*foundCorrectJetPair0_matched/nPairsTot0_matched << "% of cases." << std::endl;
  std::cout << "**MATCHED** Found correct pair 1 in " << foundCorrectJetPair1_matched << "/" << nPairsTot1_matched << " = " << (float)100.*foundCorrectJetPair1_matched/nPairsTot1_matched << "% of cases." << std::endl;



  TGraphAsymmErrors* gr_eff_combinedIsoRel_vs_rho = new TGraphAsymmErrors(0);
  gr_eff_combinedIsoRel_vs_rho->SetName("eff_combinedIsoRel_vs_rho");
  gr_eff_combinedIsoRel_vs_rho->BayesDivide( h1_passedCombinedIsoRel005, h1_passedCombinedIsoRel015 );


  h1_nCounter->SetBinContent(1, nCounter_);
  h1_nCounterW->SetBinContent(1, nCounterW_);
  h1_nCounterPU->SetBinContent(1, nCounterPU_);




  // write all stuff in files:

  outFile_->cd();

  tree_passedEvents->Write();

  gr_eff_combinedIsoRel_vs_rho->Write();

  h1_nCounter->Write();
  h1_nCounterW->Write();
  h1_nCounterPU->Write();



  h1_nvertex->Write();
  h1_nvertex_PUW->Write();
  h1_nvertex_PUW_ave->Write();

  h1_leptType->Write();

  h1_nSignalJets->Write();

  h1_nJets->Write();
  h1_nBJets->Write();

  h1_mjj0->Write();
  h1_mjj1->Write();

  h1_mjjDiff->Write();
  h1_mjjSum->Write();

//h1_pfMet->Write();

//h1_mT1->Write();
//h1_mT2->Write();
//h1_mTMax->Write();

//h1_metSignificance->Write();


//h1_rhoPF_presel->Write();
//h1_rhoPF->Write();

//

//h1_ptLept1->Write();
//h1_ptLept2->Write();
//h1_etaLept1->Write();
//h1_etaLept2->Write();
//h1_combinedIsoRelLept1->Write();
//h1_combinedIsoRelLept2->Write();

//h1_passedCombinedIsoRel015->Write();
//h1_passedCombinedIsoRel005->Write();

//h1_etaMu->Write();
//h1_etaEle->Write();


//h1_deltaRll->Write();



//h1_nJets->Write();


//h1_bTagJet1->Write();
//h1_bTagJet2->Write();

//h1_kinfit_chiSquare->Write();
//h1_kinfit_chiSquareProb->Write();

//h1_ptJet1->Write();
//h1_ptJet2->Write();
//h1_ptJet3->Write();
//h1_ptJet4->Write();

//h1_etaJet1->Write();
//h1_etaJet2->Write();
//h1_etaJet3->Write();
//h1_etaJet4->Write();

//h1_QGLikelihoodJet1->Write();
//h1_QGLikelihoodJet2->Write();
//h1_QGLikelihoodJet3->Write();
//h1_QGLikelihoodJet4->Write();

//h1_partFlavorJet1->Write();
//h1_partFlavorJet2->Write();
//h1_partFlavorJet3->Write();
//h1_partFlavorJet4->Write();


//h1_deltaRbb->Write();
//h1_deltaRqq->Write();

//h1_deltaR_b_lept_max->Write();
//h1_deltaR_b_lept_min->Write();

//h1_mbb->Write();
//h1_mqq->Write();
//
//h1_ptbb->Write();
//h1_ptqq->Write();
  

  outFile_->Close();


} // finalize()



void Ntp1Finalizer_TTW::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;


  ptLept1_thresh_ = 20.;
  ptLept2_thresh_ = 20.;
  etaLept1_thresh_ = 3.;
  etaLept2_thresh_ = 3.;

  mll_thresh_ = 8.;

  pfMet_thresh_ = 0.;
  ht_thresh_ = 0.;

  nbjets_thresh_ = 0;
  nbjetsmed_thresh_ = 0;
  njets_thresh_ = 0;

  ptJet_thresh_ = 20.;
  etaJet_thresh_ = 2.4;


  if( selectionType_=="presel" ) {

    // nothing to change

  } else if( selectionType_=="sel1" ) {

    ptLept1_thresh_ = 55.;
    ptLept2_thresh_ = 30.;

    njets_thresh_ = 3;
    nbjetsmed_thresh_ = 1;
    ht_thresh_ = 100.;

  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType




//bool Ntp1Finalizer_TTW::passedBTag( float btag1, float btag2, const std::string& btagger ) {
//
//  float loose_thresh = -99999.;
//  float med_thresh = -99999.;
//
//  if( btagger=="TCHE" ) {
//    loose_thresh = 1.7;
//    med_thresh = 3.3;
//  } else if( btagger=="JP" ) {
//    loose_thresh = 
//    med_thresh = 3.3;
//  } else if( btagger=="SSVHE" ) {
//    med_thresh = 2.0;
//  }
//
//
//  bool returnBool = false;
//
//  if( btagSelectionType_=="looseloose" ) {
//
//    returnBool = ( btag1>loose_thresh && btag2>loose_thresh );  
//
//  } else if( btagSelectionType_=="loosemed" ) {
//
//    returnBool = ( (btag1>med_thresh && btag2>loose_thresh) || (btag1>loose_thresh && btag2>med_thresh) );  
//
//  } else if( btagSelectionType_=="medmed" ) {
//
//    returnBool = ( btag1>med_thresh && btag2>med_thresh );  
// 
//  } else if( btagSelectionType_=="loosenone" ) {
//
//    returnBool = ( btag1>loose_thresh || btag2>loose_thresh );  
// 
//  }
//
//  return returnBool;
//
//}



bool isMatchedPair( AnalysisJet *jet1, AnalysisJet *jet2 ) {


  // jet 1
  if( !isMatchedToPart(*jet1) ) return false;
//std::cout << "passed deltaR" << std::endl;
//std::cout << "jet 1 part pdg id: " << jet1->pdgIdPart << std::endl;
  if( abs(jet1->pdgIdPart) > 4 ) return false; // only udsc
//std::cout << "passed jet1 part pdg id" << std::endl;
//std::cout << "jet 1 part mom pdg id: " << jet1->pdgIdPartMom << std::endl;
  if( abs(jet1->pdgIdPartMom) != 24 ) return false; // only W->jj
//std::cout << "passed jet1 part mom pdg id" << std::endl;

  // jet 2
  if( !isMatchedToPart(*jet2) ) return false;
//std::cout << "jet 2 part pdg id: " << jet2->pdgIdPart << std::endl;
  if( abs(jet2->pdgIdPart) > 4 ) return false; // only udsc
//std::cout << "passed jet2 part pdg id" << std::endl;
//std::cout << "jet 2 part mom pdg id: " << jet2->pdgIdPartMom << std::endl;

  // compare moms (for info on your mom, see http://i.imgur.com/SFfXU.gif)
  //std::cout << "same mom? " << jet1->pdgIdPartMom << " " <<  jet2->pdgIdPartMom << std::endl;
  if( jet1->pdgIdPartMom != jet2->pdgIdPartMom ) return false; // same W
  //std::cout << "same mommom? " << jet1->pdgIdPartMomMom << " " <<  jet2->pdgIdPartMomMom << std::endl;
  if( jet1->pdgIdPartMomMom != jet2->pdgIdPartMomMom ) return false; // same W

  return true;

}



bool isMatchedToPart( AnalysisJet jet ) {

  TLorentzVector part;
  part.SetPtEtaPhiE( jet.ptPart, jet.etaPart, jet.phiPart, jet.ePart );

  float deltaR = jet.DeltaR(part);

  return (deltaR<0.5);

}




std::vector<int> getJets( const std::string& choice, std::vector< AnalysisJet > jets ) {


  std::vector<int> firstPair  = getSingleJetPair( choice, jets );
  std::vector<int> secondPair = getSingleJetPair( choice, jets, &firstPair );

  std::vector<int> returnIndexes;
  returnIndexes.push_back(firstPair[0]);
  returnIndexes.push_back(firstPair[1]);
  returnIndexes.push_back(secondPair[0]);
  returnIndexes.push_back(secondPair[1]);

  return returnIndexes;

}



std::vector<int> getSingleJetPair( const std::string& choice, std::vector<AnalysisJet> jets, std::vector<int> *vetoIndexes ) {

  
  std::vector<int> jetIndexes;
  
  if( choice == "leading" ) {

    for( unsigned iJet=0; iJet<jets.size() && jetIndexes.size()<2; ++iJet ) {

      if( matchedToVeto( iJet, vetoIndexes ) ) continue;

      jetIndexes.push_back(iJet);

    }
  
  } else if( choice == "closestToLead" ) {
  

    float deltaRmin = 999.;
    int indexJet0 = -1;
    int indexJet1 = -1;
  
    // start from hardest, take closest
    for( unsigned iJet=0; iJet<jets.size(); ++iJet ) {

      if( matchedToVeto( iJet, vetoIndexes ) ) continue;
  
      if( indexJet0==-1 ) {

        indexJet0 = iJet; //hardest non-vetoed jet

      } else { //find closest

        float thisDeltaR = jets[iJet].DeltaR(jets[indexJet0]);
    
        if( thisDeltaR<deltaRmin ) {
          deltaRmin = thisDeltaR;
          indexJet1 = iJet;
        }
    
      } //for jets i - looking for closest

    }

    jetIndexes.push_back(indexJet0);
    jetIndexes.push_back(indexJet1);

  
  } else if( choice == "closest" ) {
  
    float deltaRmin = 999.;
    int indexJet0 = -1;
    int indexJet1 = -1;
  
    for( unsigned iJet=0; iJet<jets.size(); ++iJet ) {

      if( matchedToVeto( iJet, vetoIndexes ) ) continue;

      for( unsigned jJet=iJet+1; jJet<jets.size(); ++jJet ) {
  
        if( matchedToVeto( jJet, vetoIndexes ) ) continue;

        float thisDeltaR = jets[iJet].DeltaR(jets[jJet]);
  
        if( thisDeltaR<deltaRmin ) {
          deltaRmin = thisDeltaR;
          indexJet0 = iJet;
          indexJet1 = jJet;
        }
  
      } //for jets j - looking for closest
    } //for jets i - looking for closest

    jetIndexes.push_back(indexJet0);
    jetIndexes.push_back(indexJet1);

  
  } else if( choice == "bestWmass" ) {
  
    float closestMass = 10000.;
    float wMass = 80.4;
  
    int indexJet0 = -1;
    int indexJet1 = -1;
  
    for( unsigned iJet=0; iJet<jets.size(); ++iJet ) {

      if( matchedToVeto( iJet, vetoIndexes ) ) continue;

      for( unsigned jJet=iJet+1; jJet<jets.size(); ++jJet ) {
  
        if( matchedToVeto( iJet, vetoIndexes ) ) continue;

        TLorentzVector thisDiJet = jets[iJet] + jets[jJet];
        float thisMass = thisDiJet.M();
  
        if( fabs(thisMass-wMass)<fabs(closestMass-wMass) ) {
          closestMass = thisMass;
          indexJet0 = iJet;
          indexJet1 = jJet;
        }
  
      } //for jets j - looking for closest
    } //for jets i - looking for closest
  
    jetIndexes.push_back(indexJet0);
    jetIndexes.push_back(indexJet1);


  } else {
  
    std::cout << "Unkown jet choice criterium: '" << choice << "'. Exiting." << std::endl;
    exit(3939);
  
  } // jet choice
   
  
  
  return jetIndexes;

}




bool matchedToVeto( int iJet, std::vector<int> *vetoIndexes ) {

  if( vetoIndexes==0 ) return false;

  for( unsigned iVeto=0; iVeto<vetoIndexes->size(); ++iVeto ) 
    if( iJet==vetoIndexes->at(iVeto) ) return true;

  return false;

}



bool isSignalJet( AnalysisJet jet ) {

  if( !isMatchedToPart(jet) ) return false;
  if( abs(jet.pdgIdPart) > 4 ) return false; // only udsc
  if( abs(jet.pdgIdPartMom) != 24 ) return false; // only W->jj

  return true;

}

