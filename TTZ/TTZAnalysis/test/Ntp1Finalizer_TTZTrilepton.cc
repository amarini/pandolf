#include "Ntp1Finalizer_TTZTrilepton.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TRandom3.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"
#include "HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.h"
#include "KinematicFit/DiJetKinFitter.h"
#include "BTagSFUtil/BTagSFUtil.h"

#include "PUWeight.h"




bool USE_MC_MASS=false;

int DEBUG_EVENTNUMBER = -1;


struct ValueAndError {

  float val;
  float err;

};


ValueAndError getMuonHLTSF_DoubleTrigger( float pt, float eta, const std::string& runPeriod );
ValueAndError getMuonHLTSF_SingleTrigger( float pt, float eta, const std::string& runPeriod );
ValueAndError getEventHLTSF( ValueAndError effSingle1, ValueAndError effSingle2, ValueAndError effDouble1, ValueAndError effDouble2 );
ValueAndError getMuonRecoSF( float pt, float eta );
ValueAndError getMuonIsoSF( float pt, float eta );
ValueAndError getElectronRecoSF( float pt, float eta );
ValueAndError getElectronIsoSF( float pt, float eta );
float getJERSF( float eta );


// constructor:

Ntp1Finalizer_TTZTrilepton::Ntp1Finalizer_TTZTrilepton( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, const std::string& leptType ) : Ntp1Finalizer( "TTZTrilepton", dataset, leptType ) {

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

  leptSyst_ = 0;
  jes_ = 0;
  btagSyst_ = 0;
  jer_ = false;

  setSelectionType(selectionType);

}




void Ntp1Finalizer_TTZTrilepton::finalize( ) {

  //if( outFile_==0 ) this->createOutputFile();
  
  Int_t run;
  tree_->SetBranchAddress("run", &run);
  tree_->GetEntry(0);
  bool isMC = (run < 160000);

  std::string fullFlags = selectionType_ + "_" + bTaggerType_;
  fullFlags = fullFlags + "_" + leptType_;
  if( jes_ == 1 ) fullFlags = fullFlags + "_JESUP";
  else if( jes_ == -1 ) fullFlags = fullFlags + "_JESDOWN";
  else if( jes_ != 0 ) {
    char fullFlags_char[100];
    if( jes_ < 0 )
      sprintf( fullFlags_char, "%s_JESDOWN%d", fullFlags.c_str(), jes_ );
    else
      sprintf( fullFlags_char, "%s_JESUP%d", fullFlags.c_str(), jes_ );
    std::string fullFlags_tmp(fullFlags_char);
    fullFlags = fullFlags_tmp;
  }
  if( btagSyst_ == 1 ) fullFlags = fullFlags + "_BTagUP";
  else if( btagSyst_ == -1 ) fullFlags = fullFlags + "_BTagDOWN";
  else if( btagSyst_ != 0 ) {
    char fullFlags_char[100];
    if( btagSyst_ < 0 )
      sprintf( fullFlags_char, "%s_BTagDOWN%d", fullFlags.c_str(), btagSyst_ );
    else
      sprintf( fullFlags_char, "%s_BTagUP%d", fullFlags.c_str(), btagSyst_ );
    std::string fullFlags_tmp(fullFlags_char);
    fullFlags = fullFlags_tmp;
  }
  if( leptSyst_ == 1 ) fullFlags = fullFlags + "_LeptUP";
  else if( leptSyst_ == -1 ) fullFlags = fullFlags + "_LeptDOWN";
  else if( leptSyst_ != 0 ) {
    char fullFlags_char[100];
    if( leptSyst_ < 0 )
      sprintf( fullFlags_char, "%s_LeptDOWN%d", fullFlags.c_str(), leptSyst_ );
    else
      sprintf( fullFlags_char, "%s_LeptUP%d", fullFlags.c_str(), leptSyst_ );
    std::string fullFlags_tmp(fullFlags_char);
    fullFlags = fullFlags_tmp;
  }
  if( jer_ ) fullFlags = fullFlags + "_JERUP";
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

  TH1D* h1_pfMet_presel = new TH1D("pfMet_presel", "", 500, 0., 500.);
  h1_pfMet_presel->Sumw2();
  TH1D* h1_pfMet = new TH1D("pfMet", "", 500, 0., 500.);
  h1_pfMet->Sumw2();
  TH1D* h1_ht = new TH1D("ht", "", 1000, 0., 1000.);
  h1_ht->Sumw2();

  TH1D* h1_metSignificance= new TH1D("metSignificance", "", 80, 0., 40.);
  h1_metSignificance->Sumw2();


  TH1D* h1_rhoPF_noPUW = new TH1D("rhoPF_noPUW", "", 50, 0., 30.);
  h1_rhoPF_noPUW->Sumw2();
  TH1D* h1_rhoPF_prepresel = new TH1D("rhoPF_prepresel", "", 50, 0., 30.);
  h1_rhoPF_prepresel->Sumw2();
  TH1D* h1_rhoPF_presel = new TH1D("rhoPF_presel", "", 50, 0., 30.);
  h1_rhoPF_presel->Sumw2();
  TH1D* h1_rhoPF = new TH1D("rhoPF", "", 50, 0., 30.);
  h1_rhoPF->Sumw2();


  TH1D* h1_ptLeptZ1 = new TH1D("ptLeptZ1", "", 500, 20., 520.);
  h1_ptLeptZ1->Sumw2();
  TH1D* h1_ptLeptZ2 = new TH1D("ptLeptZ2", "", 200, 20., 220.);
  h1_ptLeptZ2->Sumw2();
  TH1D* h1_etaLeptZ1 = new TH1D("etaLeptZ1", "", 50, -2.5, 2.5);
  h1_etaLeptZ1->Sumw2();
  TH1D* h1_etaLeptZ2 = new TH1D("etaLeptZ2", "", 50, -2.5, 2.5);
  h1_etaLeptZ2->Sumw2();
  TH1D* h1_combinedIsoRelLeptZ1 = new TH1D("combinedIsoRelLeptZ1", "", 100, 0., 1.);
  h1_combinedIsoRelLeptZ1->Sumw2();
  TH1D* h1_combinedIsoRelLeptZ2 = new TH1D("combinedIsoRelLeptZ2", "", 100, 0., 1.);
  h1_combinedIsoRelLeptZ2->Sumw2();

  TH1D* h1_ptLept3_presel = new TH1D("ptLept3_presel", "", 200, 0., 200.);
  h1_ptLept3_presel->Sumw2();
  TH1D* h1_etaLept3_presel = new TH1D("etaLept3_presel", "", 50, -2.5, 2.5);
  h1_etaLept3_presel->Sumw2();
  TH1D* h1_leptTypeLept3_presel = new TH1D("leptTypeLept3_presel", "", 2, -0.5, 1.5);
  h1_leptTypeLept3_presel->Sumw2();
  TH1D* h1_combinedIsoRelLept3_presel = new TH1D("combinedIsoRelLept3_presel", "", 100, 0., 1.);
  h1_combinedIsoRelLept3_presel->Sumw2();

  TH1D* h1_ptLept3 = new TH1D("ptLept3", "", 200, 10., 220.);
  h1_ptLept3->Sumw2();
  TH1D* h1_etaLept3 = new TH1D("etaLept3", "", 50, -2.5, 2.5);
  h1_etaLept3->Sumw2();

  TH1D* h1_deltaRZll = new TH1D("deltaRZll", "", 500, 0., 5.);
  h1_deltaRZll->Sumw2();

  TH1D* h1_ptZll = new TH1D("ptZll", "", 400., 0., 400.);
  h1_ptZll->Sumw2();
  TH1D* h1_etaZll = new TH1D("etaZll", "", 200, -5., 5.);
  h1_etaZll->Sumw2();
  TH1D* h1_mZll_prepresel = new TH1D("mZll_prepresel", "", 220, 50., 160.);
  h1_mZll_prepresel->Sumw2();
  TH1D* h1_mZll_prepresel_ELE = new TH1D("mZll_prepresel_ELE", "", 220, 50., 160.);
  h1_mZll_prepresel_ELE->Sumw2();
  TH1D* h1_mZll_prepresel_MU = new TH1D("mZll_prepresel_MU", "", 220, 50., 160.);
  h1_mZll_prepresel_MU->Sumw2();
  TH1D* h1_mZll_prepresel_antibtag = new TH1D("mZll_prepresel_antibtag", "", 220, 50., 160.);
  h1_mZll_prepresel_antibtag->Sumw2();
  TH1D* h1_mZll_OF_prepresel = new TH1D("mZll_OF_prepresel", "", 220, 50., 160.);
  h1_mZll_OF_prepresel->Sumw2();
  TH1D* h1_mZll_OF2_prepresel = new TH1D("mZll_OF2_prepresel", "", 220, 50., 160.);
  h1_mZll_OF2_prepresel->Sumw2();
  TH1D* h1_mZll_OF3_prepresel = new TH1D("mZll_OF3_prepresel", "", 220, 50., 160.);
  h1_mZll_OF3_prepresel->Sumw2();
  TH1D* h1_mZll_presel = new TH1D("mZll_presel", "", 220, 50., 160.);
  h1_mZll_presel->Sumw2();
  TH1D* h1_mZll_presel_antibtag = new TH1D("mZll_presel_antibtag", "", 220, 50., 160.);
  h1_mZll_presel_antibtag->Sumw2();
  TH1D* h1_mZll = new TH1D("mZll", "", 220, 50., 160.);
  h1_mZll->Sumw2();


  TH1D* h1_nJets_prepresel = new TH1D("nJets_prepresel", "", 11, -0.5, 10.5);
  h1_nJets_prepresel->Sumw2();
  TH1D* h1_nJets_presel = new TH1D("nJets_presel", "", 11, -0.5, 10.5);
  h1_nJets_presel->Sumw2();
  TH1D* h1_nJets = new TH1D("nJets", "", 7, 3.5, 10.5);
  h1_nJets->Sumw2();

  TH1D* h1_nBJets_loose_presel = new TH1D("nBJets_loose_presel", "", 11, -0.5, 10.5);
  h1_nBJets_loose_presel->Sumw2();
  TH1D* h1_nBJets_loose = new TH1D("nBJets_loose", "", 7, 3.5, 10.5);
  h1_nBJets_loose->Sumw2();

  TH1D* h1_nBJets_medium_presel = new TH1D("nBJets_medium_presel", "", 11, -0.5, 10.5);
  h1_nBJets_medium_presel->Sumw2();
  TH1D* h1_nBJets_medium = new TH1D("nBJets_medium", "", 7, 3.5, 10.5);
  h1_nBJets_medium->Sumw2();


  TH1D* h1_bTagJetB1 = new TH1D("bTagJetB1", "", 420, -1., 20.);
  h1_bTagJetB1->Sumw2();
  TH1D* h1_bTagJetB2 = new TH1D("bTagJetB2", "", 420, -1., 20.);
  h1_bTagJetB2->Sumw2();


  TH1D* h1_ptJetB1 = new TH1D("ptJetB1", "", 400, 20., 420.);
  h1_ptJetB1->Sumw2();
  TH1D* h1_ptJetB2 = new TH1D("ptJetB2", "", 400, 20., 420.);
  h1_ptJetB2->Sumw2();
  TH1D* h1_ptJet3 = new TH1D("ptJet3", "", 400, 20., 420.);
  h1_ptJet3->Sumw2();
  TH1D* h1_ptJet4 = new TH1D("ptJet4", "", 400, 20., 420.);
  h1_ptJet4->Sumw2();

  TH1D* h1_etaJetB1 = new TH1D("etaJetB1", "", 200, -5., 5.);
  h1_etaJetB1->Sumw2();
  TH1D* h1_etaJetB2 = new TH1D("etaJetB2", "", 200, -5., 5.);
  h1_etaJetB2->Sumw2();
  TH1D* h1_etaJet3 = new TH1D("etaJet3", "", 200, -5., 5.);
  h1_etaJet3->Sumw2();
  TH1D* h1_etaJet4 = new TH1D("etaJet4", "", 200, -5., 5.);
  h1_etaJet4->Sumw2();

  TH1D* h1_partFlavorJetB1 = new TH1D("partFlavorJetB1", "", 38, -15.5, 22.5);
  h1_partFlavorJetB1->Sumw2();
  TH1D* h1_partFlavorJetB2 = new TH1D("partFlavorJetB2", "", 38, -15.5, 22.5);
  h1_partFlavorJetB2->Sumw2();
  TH1D* h1_partFlavorJet3 = new TH1D("partFlavorJet3", "", 38, -15.5, 22.5);
  h1_partFlavorJet3->Sumw2();
  TH1D* h1_partFlavorJet4 = new TH1D("partFlavorJet4", "", 38, -15.5, 22.5);
  h1_partFlavorJet4->Sumw2();


  TH1D* h1_mTW = new TH1D("mTW", "", 300, 0., 300.);
  h1_mTW->Sumw2();
  TH1D* h1_mT_lZmet = new TH1D("mT_lZmet", "", 1000, 0., 1000.);
  h1_mT_lZmet->Sumw2();

  TH1D* h1_deltaRbb = new TH1D("deltaRbb", "", 500, 0., 5.);
  h1_deltaRbb->Sumw2();
  
  TH1D* h1_mb1jj = new TH1D("mb1jj", "", 1000., 0., 1000.);
  h1_mb1jj->Sumw2();
  TH1D* h1_mb2jj = new TH1D("mb2jj", "", 1000., 0., 1000.);
  h1_mb2jj->Sumw2();
  TH1D* h1_mbjj_best = new TH1D("mbjj_best", "", 1000., 0., 1000.);
  h1_mbjj_best->Sumw2();
  TH1D* h1_mbjj_max = new TH1D("mbjj_max", "", 1000., 0., 1000.);
  h1_mbjj_max->Sumw2();
  
  TH1D* h1_mb1jjZ = new TH1D("mb1jjZ", "", 1000., 0., 1000.);
  h1_mb1jjZ->Sumw2();
  TH1D* h1_mb2jjZ = new TH1D("mb2jjZ", "", 1000., 0., 1000.);
  h1_mb2jjZ->Sumw2();
  TH1D* h1_mbjjZ_best = new TH1D("mbjjZ_best", "", 1000., 0., 1000.);
  h1_mbjjZ_best->Sumw2();
  TH1D* h1_mbjjZ_max = new TH1D("mbjjZ_max", "", 1000., 0., 1000.);
  h1_mbjjZ_max->Sumw2();

  TH1D* h1_mTb1W = new TH1D("mTb1W", "", 1000, 0., 1000.);
  h1_mTb1W->Sumw2();
  TH1D* h1_mTb2W = new TH1D("mTb2W", "", 1000, 0., 1000.);
  h1_mTb2W->Sumw2();
  TH1D* h1_mTbW_best = new TH1D("mTbW_best", "", 1000, 0., 1000.);
  h1_mTbW_best->Sumw2();
  TH1D* h1_mTbW_max = new TH1D("mTbW_max", "", 1000, 0., 1000.);
  h1_mTbW_max->Sumw2();

  TH1D* h1_mTb1WZ = new TH1D("mTb1WZ", "", 1000, 0., 1000.);
  h1_mTb1WZ->Sumw2();
  TH1D* h1_mTb2WZ = new TH1D("mTb2WZ", "", 1000, 0., 1000.);
  h1_mTb2WZ->Sumw2();
  TH1D* h1_mTbWZ_best = new TH1D("mTbWZ_best", "", 1000, 0., 1000.);
  h1_mTbWZ_best->Sumw2();
  TH1D* h1_mTbWZ_max = new TH1D("mTbWZ_max", "", 1000, 0., 1000.);
  h1_mTbWZ_max->Sumw2();


  TH1D* h1_channelYields = new TH1D("channelYields", "", 5, 0., 5.);
  h1_channelYields->GetXaxis()->SetLabelSize(0.085);
  h1_channelYields->GetXaxis()->SetBinLabel(1, "(ee)e");
  h1_channelYields->GetXaxis()->SetBinLabel(2, "(ee)#mu");
  h1_channelYields->GetXaxis()->SetBinLabel(3, "(#mu#mu)e");
  h1_channelYields->GetXaxis()->SetBinLabel(4, "(#mu#mu)#mu");
  h1_channelYields->GetXaxis()->SetBinLabel(5, "Total" );
  h1_channelYields->Sumw2();



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
  //Float_t eventWeight_Zee;
  //tree_->SetBranchAddress("eventWeight_Zee", &eventWeight_Zee);
  //Float_t eventWeight_Zmm;
  //tree_->SetBranchAddress("eventWeight_Zmm", &eventWeight_Zmm);

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
  Float_t combinedIsoRelLeptZ1;
  tree_->SetBranchAddress("combinedIsoRelLeptZ1", &combinedIsoRelLeptZ1);

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
  Float_t combinedIsoRelLeptZ2;
  tree_->SetBranchAddress("combinedIsoRelLeptZ2", &combinedIsoRelLeptZ2);


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
  Float_t combinedIsoRelLept[10];
  tree_->SetBranchAddress("combinedIsoRelLept", combinedIsoRelLept);


  Int_t nJets;
  tree_->SetBranchAddress("nJets", &nJets);

  Float_t eJet[50];
  tree_->SetBranchAddress("eJet", eJet);
  Float_t ptJet[50];
  tree_->SetBranchAddress("ptJet", ptJet);
  Float_t ptUncertJet[50];
  tree_->SetBranchAddress("ptUncertJet", ptUncertJet);
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
  Int_t pdgIdPartJet[50];
  tree_->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);
  Float_t etaPartJet[50];
  tree_->SetBranchAddress("etaPartJet", etaPartJet);
  Float_t phiPartJet[50];
  tree_->SetBranchAddress("phiPartJet", phiPartJet);



  //Int_t nPart;
  //tree_->SetBranchAddress("nPart", &nPart);
  //Float_t ePart[20];
  //tree_->SetBranchAddress("ePart", ePart);
  //Float_t ptPart[20];
  //tree_->SetBranchAddress("ptPart", ptPart);
  //Float_t etaPart[20];
  //tree_->SetBranchAddress("etaPart", etaPart);
  //Float_t phiPart[20];
  //tree_->SetBranchAddress("phiPart", phiPart);
  //Int_t pdgIdPart[20];
  //tree_->SetBranchAddress("pdgIdPart", pdgIdPart);


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


  BTagSFUtil* btsfutil = new BTagSFUtil(bTaggerType_, 13);
  std::string meanminmax;
  if( btagSyst_==0 )       meanminmax = "mean";
  else if( btagSyst_==1 )  meanminmax = "max";
  else if( btagSyst_==-1 ) meanminmax = "min";
  else {
    std::cout << "Only allowed values for btagSyst are 0, 1, -1. Exiting." << std::endl;
    exit(55);
  }

  
  //QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
 
  float Zmass = 91.1876;
  float tmass = 172.9;


  std::string puType = "Summer11_S4";
  TString dataset_tstr(dataset_);
  if( dataset_tstr.Contains("Fall11") ) {
    puType = "Fall11";
    //puType = "Fall11Truth";
  }


  PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
  //PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
  std::string puFileName;
  //puFileName = "all2011AB.pileup_v2_73mb.root";
  puFileName = "Pileup_DATA_up_to_178479_SiXie.root";
  //puFileName = "PileupTruth_v2.root";
  //puFileName = "PileupObs_v2.root";
  //puFileName = "FullData_178078.root"; 


  std::cout << std::endl << "-> Using data pileup file: " << puFileName << std::endl;
  TFile* filePU = TFile::Open(puFileName.c_str());
  TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
  fPUWeight->SetDataHistogram(h1_nPU_data);


  if( dataset_tstr.Contains("spadhi") )
    fPUWeight->SetMCHistogram(h1_nPU_gen_);
  

  //if( dataset_tstr.BeginsWith("TTW") || dataset_tstr.BeginsWith("TTZ") ) {
  //  TFile* filePUMC = TFile::Open("RareSM_Sanjay-v1.root");
  //  TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("pileup");
  //  std::cout << "-> Switching MC PU file to: RareSM_Sanjay-v1.root" << std::endl;
  //  fPUWeight->SetMCHistogram(h1_nPU_mc);
  //}

  //if( dataset_tstr.BeginsWith("TTJ") && dataset_tstr.Contains("Fall11") ) {
  //  TFile* filePUMC = TFile::Open("TTJets_Fall11-S6.root");
  //  TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("pileup");
  //  std::cout << "-> Switching MC PU file to: TTJets_Fall11-S6.root" << std::endl;
  //  fPUWeight->SetMCHistogram(h1_nPU_mc);
  //}

  //if( dataset_ == "DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1" ) {
  //  TFile* filePUMC = TFile::Open("generatedpileup_Zjets_MADGRAPH_AOD423.root");
  //  TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("GenLevelInfoModule/npileup");
  //  std::cout << "-> Switching MC PU file to: generatedpileup_Zjets_MADGRAPH_AOD423.root" << std::endl;
  //  fPUWeight->SetMCHistogram(h1_nPU_mc);
  //}
  
  //TFile* filePUMC = TFile::Open("PUFile_DYJets_Summer11.root");
  //TH1F* h1_nPU_mc = (TH1F*)filePUMC->Get("pileup");
  //std::cout << "-> Switching MC PU file to: PUFile_DYJets_Summer11.root" << std::endl;
  //fPUWeight->SetMCHistogram(h1_nPU_mc);


  float mZll_t, ptZll_t;
  float ptLeptZ1_t, ptLeptZ2_t, etaLeptZ1_t, etaLeptZ2_t;
  float ptLept3_t, etaLept3_t;
  float ptJetB1_t, ptJetB2_t, etaJetB1_t, etaJetB2_t;
  //float bTagJetB1_t, bTagJetB2_t;
  //float ptJet3_t, ptJet4_t, etaJet3_t, etaJet4_t;
  float HLTSF, leptonSF;
  int leptType3;
  bool isMZllSignalRegion;
  bool passed_btag;
  int njets;
  int nBjets_loose;
  int nBjets_medium;
  float ht, mt;
  float mT_lZmet;
  float mb1jj;
  float mb2jj;
  float mbjj_best;
  float mbjj_max;
  float mb1jjZ;
  float mb2jjZ;
  float mbjjZ_best;
  float mbjjZ_max;
  float mTb1W;
  float mTb2W;
  float mTbW_best;
  float mTbW_max;
  float mTb1WZ;
  float mTb2WZ;
  float mTbWZ_best;
  float mTbWZ_max;

  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "LS", &LS, "LS/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "HLTSF", &HLTSF, "HLTSF/F" );
  tree_passedEvents->Branch( "leptonSF", &leptonSF, "leptonSF/F" );
  tree_passedEvents->Branch( "PUWeight", &eventWeightPU, "eventWeightPU/F" );
  tree_passedEvents->Branch( "pfMet", &pfMet, "pfMet/F" );
  tree_passedEvents->Branch( "leptType", &leptType, "leptType/I" );
  tree_passedEvents->Branch( "leptType3", &leptType3, "leptType3/I" );
  tree_passedEvents->Branch( "ptLeptZ1", &ptLeptZ1_t, "ptLeptZ1_t/F" );
  tree_passedEvents->Branch( "ptLeptZ2", &ptLeptZ2_t, "ptLeptZ2_t/F" );
  tree_passedEvents->Branch( "ptLept3", &ptLept3_t, "ptLept3_t/F" );
  tree_passedEvents->Branch( "etaLeptZ1", &etaLeptZ1_t, "etaLeptZ1_t/F" );
  tree_passedEvents->Branch( "etaLeptZ2", &etaLeptZ2_t, "etaLeptZ2_t/F" );
  tree_passedEvents->Branch( "ptZll", &ptZll_t, "ptZll_t/F" );
  tree_passedEvents->Branch( "mZll", &mZll_t, "mZll_t/F" );
  tree_passedEvents->Branch( "etaLept3", &etaLept3_t, "etaLept3_t/F" );
  tree_passedEvents->Branch( "ptJetB1", &ptJetB1_t, "ptJetB1_t/F" );
  //tree_passedEvents->Branch( "ptJetB2", &ptJetB2_t, "ptJetB2_t/F" );
  //tree_passedEvents->Branch( "bTagJetB1", &bTagJetB1_t, "bTagJetB1_t/F" );
  //tree_passedEvents->Branch( "bTagJetB2", &bTagJetB2_t, "bTagJetB2_t/F" );
  //tree_passedEvents->Branch( "ptJet3", &ptJet3_t, "ptJet3_t/F" );
  //tree_passedEvents->Branch( "ptJet4", &ptJet4_t, "ptJet4_t/F" );
  //tree_passedEvents->Branch( "etaJetB1", &etaJetB1_t, "etaJetB1_t/F" );
  //tree_passedEvents->Branch( "etaJetB2", &etaJetB2_t, "etaJetB2_t/F" );
  //tree_passedEvents->Branch( "etaJet3", &etaJet3_t, "etaJet3_t/F" );
  //tree_passedEvents->Branch( "etaJet4", &etaJet4_t, "etaJet4_t/F" );
  tree_passedEvents->Branch( "ht", &ht, "ht/F" );
  tree_passedEvents->Branch( "mt", &mt, "mt/F" );
  tree_passedEvents->Branch( "isMZllSignalRegion", &isMZllSignalRegion, "isMZllSignalRegion/O" );
  tree_passedEvents->Branch( "passed_btag", &passed_btag, "passed_btag/O" );
  tree_passedEvents->Branch( "njets", &njets, "njets/I" );
  tree_passedEvents->Branch( "nBjets_loose", &nBjets_loose, "nBjets_loose/I" );
  tree_passedEvents->Branch( "nBjets_medium", &nBjets_medium, "nBjets_medium/I" );
  tree_passedEvents->Branch("mT_lZmet"  , &mT_lZmet  , "mT_lZmet/F");
  tree_passedEvents->Branch("mb1jj"     , &mb1jj     , "mb1jj/F");
  tree_passedEvents->Branch("mb2jj"     , &mb2jj     , "mb2jj/F");
  tree_passedEvents->Branch("mbjj_best" , &mbjj_best , "mbjj_best/F");
  tree_passedEvents->Branch("mbjj_max" , &mbjj_max , "mbjj_max/F");
  tree_passedEvents->Branch("mb1jjZ"    , &mb1jjZ    , "mb1jjZ/F");
  tree_passedEvents->Branch("mb2jjZ"    , &mb2jjZ    , "mb2jjZ/F");
  tree_passedEvents->Branch("mbjjZ_best", &mbjjZ_best, "mbjjZ_best/F");
  tree_passedEvents->Branch("mbjjZ_max", &mbjjZ_max, "mbjjZ_max/F");
  tree_passedEvents->Branch("mTb1W"     , &mTb1W     , "mTb1W/F");
  tree_passedEvents->Branch("mTb2W"     , &mTb2W     , "mTb2W/F");
  tree_passedEvents->Branch("mTbW_best" , &mTbW_best , "mTbW_best/F");
  tree_passedEvents->Branch("mTbW_max" , &mTbW_max , "mTbW_max/F");
  tree_passedEvents->Branch("mTb1WZ"    , &mTb1WZ    , "mTb1WZ/F");
  tree_passedEvents->Branch("mTb2WZ"    , &mTb2WZ    , "mTb2WZ/F");
  tree_passedEvents->Branch("mTbWZ_best", &mTbWZ_best, "mTbWZ_best/F");
  tree_passedEvents->Branch("mTbWZ_max", &mTbWZ_max, "mTbWZ_max/F");
      



ofstream ofs("run_event.txt");


  TRandom3* rand = new TRandom3(1313);


  std::cout << std::endl << std::endl;
  std::cout << "+++ BEGINNING ANALYSIS LOOP" << std::endl;
  std::cout << "----> DATASET: " << dataset_ << std::endl;
  std::cout << "----> SELECTION: " << selectionType_ << std::endl;
  std::cout << "----> B-TAGGER: " << bTaggerType_ << std::endl;
  std::cout << "----> LEPT TYPE: " << leptType_ << std::endl;
  if( jes_!= 0 )
    std::cout << "----> JES SYST: " << jes_ << " sigma" << std::endl;
  if( btagSyst_!= 0 )
    std::cout << "----> BTag SYST: " << btagSyst_ << " sigma" << std::endl;
  if( leptSyst_!= 0 )
    std::cout << "----> LEPT SYST: " << leptSyst_ << " sigma" << std::endl;
  if( jer_ )
    std::cout << "----> JER SYST ON" << std::endl;
  std::cout << std::endl << std::endl;



  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree_->GetEntry(iEntry);

//std::cout << "new event" << std::endl;

    if( eventWeight <= 0. ) eventWeight = 1.;

    if( leptType_!="ALL" ) {
      if( leptType_=="ELE" && leptType==0 ) continue;
      if( leptType_=="MU" && leptType==1 ) continue;
    }

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

   

    

    h1_nvertex->Fill(nvertex, eventWeight);
    h1_rhoPF_noPUW->Fill( rhoPF, eventWeight);

    if( isMC ) {

      HLTSF = 1.;
      float HLTSF_err = 0.;

      float leptonRecoSF = 1.;
      float leptonIsoSF  = 1.;

      float leptonRecoSF_err = 0.;
      float leptonIsoSF_err  = 0.;

      // scale factor for double mu triggers:
      if( leptType==0 ) {

        ValueAndError effDouble1_Run2011A = getMuonHLTSF_DoubleTrigger( ptLeptZ1, etaLeptZ1, "Run2011A" );
        ValueAndError effDouble2_Run2011A = getMuonHLTSF_DoubleTrigger( ptLeptZ2, etaLeptZ2, "Run2011A" );

        ValueAndError effDouble1_Run2011B = getMuonHLTSF_DoubleTrigger( ptLeptZ1, etaLeptZ1, "Run2011B" );
        ValueAndError effDouble2_Run2011B = getMuonHLTSF_DoubleTrigger( ptLeptZ2, etaLeptZ2, "Run2011B" );

        ValueAndError effSingle1_Run2011A1 = getMuonHLTSF_SingleTrigger( ptLeptZ1, etaLeptZ1, "Run2011A1");
        ValueAndError effSingle2_Run2011A1 = getMuonHLTSF_SingleTrigger( ptLeptZ2, etaLeptZ2, "Run2011A1");

        ValueAndError effSingle1_Run2011A2 = getMuonHLTSF_SingleTrigger( ptLeptZ1, etaLeptZ1, "Run2011A2");
        ValueAndError effSingle2_Run2011A2 = getMuonHLTSF_SingleTrigger( ptLeptZ2, etaLeptZ2, "Run2011A2");

        ValueAndError effSingle1_Run2011A3 = getMuonHLTSF_SingleTrigger( ptLeptZ1, etaLeptZ1, "Run2011A3");
        ValueAndError effSingle2_Run2011A3 = getMuonHLTSF_SingleTrigger( ptLeptZ2, etaLeptZ2, "Run2011A3");

        ValueAndError effSingle1_Run2011B = getMuonHLTSF_SingleTrigger( ptLeptZ1, etaLeptZ1, "Run2011B");
        ValueAndError effSingle2_Run2011B = getMuonHLTSF_SingleTrigger( ptLeptZ2, etaLeptZ2, "Run2011B");


        ValueAndError HLTSF_Run2011A1 = getEventHLTSF( effSingle1_Run2011A1, effSingle2_Run2011A1, effDouble1_Run2011A, effDouble2_Run2011A );
        ValueAndError HLTSF_Run2011A2 = getEventHLTSF( effSingle1_Run2011A2, effSingle2_Run2011A2, effDouble1_Run2011A, effDouble2_Run2011A );
        ValueAndError HLTSF_Run2011A3 = getEventHLTSF( effSingle1_Run2011A3, effSingle2_Run2011A3, effDouble1_Run2011A, effDouble2_Run2011A );
        ValueAndError HLTSF_Run2011B  = getEventHLTSF( effSingle1_Run2011B, effSingle2_Run2011B, effDouble1_Run2011B, effDouble2_Run2011B );


        // weighted average over full run (weighted with lumi):
        HLTSF = (217.*HLTSF_Run2011A1.val + 920.*HLTSF_Run2011A2.val + 1000.*HLTSF_Run2011A3.val + 2100.*HLTSF_Run2011B.val)/(217.+920.+1000.+2100.);
        HLTSF_err = (217.*217.*HLTSF_Run2011A1.err*HLTSF_Run2011A1.err + 
                     920.*920.*HLTSF_Run2011A2.err*HLTSF_Run2011A2.err + 
                     1000.*1000.*HLTSF_Run2011A3.err*HLTSF_Run2011A3.err +
                     2100.*2100.*HLTSF_Run2011B.err*HLTSF_Run2011B.err);
        HLTSF_err = HLTSF_err / ( (217.+920.+1000.+2100.)*(217.+920.+1000.+2100.) );
        HLTSF_err = sqrt( HLTSF_err );


        ValueAndError muonRecoSF1 = getMuonRecoSF( ptLeptZ1, etaLeptZ1 );
        ValueAndError muonRecoSF2 = getMuonRecoSF( ptLeptZ2, etaLeptZ2 );
        float muonRecoSF = muonRecoSF1.val*muonRecoSF2.val;
        float muonRecoSF_err = muonRecoSF1.val*muonRecoSF1.val*muonRecoSF2.err*muonRecoSF2.err +
                               muonRecoSF2.val*muonRecoSF2.val*muonRecoSF1.err*muonRecoSF1.err;
        muonRecoSF_err = sqrt( muonRecoSF_err );

        leptonRecoSF = muonRecoSF;
        leptonRecoSF_err = muonRecoSF_err;


        ValueAndError muonIsoSF1 = getMuonIsoSF( ptLeptZ1, etaLeptZ1 );
        ValueAndError muonIsoSF2 = getMuonIsoSF( ptLeptZ2, etaLeptZ2 );
        float muonIsoSF = muonIsoSF1.val*muonIsoSF2.val;
        float muonIsoSF_err = muonIsoSF1.val*muonIsoSF1.val*muonIsoSF2.err*muonIsoSF2.err +
                              muonIsoSF2.val*muonIsoSF2.val*muonIsoSF1.err*muonIsoSF1.err;
        muonIsoSF_err = sqrt( muonIsoSF_err );

        leptonIsoSF = muonIsoSF;
        leptonIsoSF_err = muonIsoSF_err;


      } else if( leptType==1 ) { //electrons

        HLTSF = 1.;
        HLTSF_err = 0.01;

        ValueAndError electronRecoSF1 = getElectronRecoSF( ptLeptZ1, etaLeptZ1 );
        ValueAndError electronRecoSF2 = getElectronRecoSF( ptLeptZ2, etaLeptZ2 );
        float electronRecoSF = electronRecoSF1.val*electronRecoSF2.val;
        float electronRecoSF_err = electronRecoSF1.val*electronRecoSF1.val*electronRecoSF2.err*electronRecoSF2.err +
                               electronRecoSF2.val*electronRecoSF2.val*electronRecoSF1.err*electronRecoSF1.err;
        electronRecoSF_err = sqrt( electronRecoSF_err );

        leptonRecoSF = electronRecoSF;
        leptonRecoSF_err = electronRecoSF_err;


        ValueAndError electronIsoSF1 = getElectronIsoSF( ptLeptZ1, etaLeptZ1 );
        ValueAndError electronIsoSF2 = getElectronIsoSF( ptLeptZ2, etaLeptZ2 );
        float electronIsoSF = electronIsoSF1.val*electronIsoSF2.val;
        float electronIsoSF_err = electronIsoSF1.val*electronIsoSF1.val*electronIsoSF2.err*electronIsoSF2.err +
                              electronIsoSF2.val*electronIsoSF2.val*electronIsoSF1.err*electronIsoSF1.err;
        electronIsoSF_err = sqrt( electronIsoSF_err );

        leptonIsoSF = electronIsoSF;
        leptonIsoSF_err = electronIsoSF_err;

      } else if( leptType==2 ) { // e+mu-

        //float effSingle_Run2011A  = ( chargeLeptZ1<0 ) ? getMuonHLTSF_DoubleTrigger( ptLeptZ1, etaLeptZ1, "Run2011A") : getMuonHLTSF_DoubleTrigger( ptLeptZ2, etaLeptZ2, "Run2011A");
        //float effSingle_Run2011B  = ( chargeLeptZ1<0 ) ? getMuonHLTSF_DoubleTrigger( ptLeptZ1, etaLeptZ1, "Run2011B") : getMuonHLTSF_DoubleTrigger( ptLeptZ2, etaLeptZ2, "Run2011B");

        //HLTSF = ( (217.+920.+1000.)*effSingle_Run2011A + 2100.*effSingle_Run2011B)/(217.+920.+1000.+2100.);
        HLTSF = 1.; // not needed as signal doesnt have OF leptons and ttbar is scaled
        HLTSF_err = 0.; // not needed as signal doesnt have OF leptons and ttbar is scaled

        leptonRecoSF = 1.;
        leptonRecoSF_err = 0.;

        leptonIsoSF = 1.;
        leptonIsoSF_err = 0.;

      } else if( leptType==3 ) { // e-mu+

        //float effSingle_Run2011A = ( chargeLeptZ1>0 ) ? getMuonHLTSF_DoubleTrigger( ptLeptZ1, etaLeptZ1, "Run2011A") : getMuonHLTSF_DoubleTrigger( ptLeptZ2, etaLeptZ2, "Run2011A");
        //float effSingle_Run2011B = ( chargeLeptZ1>0 ) ? getMuonHLTSF_DoubleTrigger( ptLeptZ1, etaLeptZ1, "Run2011B") : getMuonHLTSF_DoubleTrigger( ptLeptZ2, etaLeptZ2, "Run2011B");

        //HLTSF = ( (217.+920.+1000.)*effSingle_Run2011A + 2100.*effSingle_Run2011B)/(217.+920.+1000.+2100.);
        HLTSF = 1.; // not needed as signal doesnt have OF leptons and ttbar is scaled
        HLTSF_err = 0.; // not needed as signal doesnt have OF leptons and ttbar is scaled

        leptonRecoSF = 1.;
        leptonRecoSF_err = 0.;

        leptonIsoSF = 1.;
        leptonIsoSF_err = 0.;

      }

      leptonSF = leptonRecoSF*leptonIsoSF*HLTSF;
      float leptonSF_err = leptonRecoSF*leptonRecoSF*leptonIsoSF*leptonIsoSF*HLTSF_err*HLTSF_err + 
                           leptonRecoSF*leptonRecoSF*leptonIsoSF_err*leptonIsoSF_err*HLTSF*HLTSF + 
                           leptonRecoSF_err*leptonRecoSF_err*leptonIsoSF*leptonIsoSF*HLTSF*HLTSF ;
      leptonSF_err = sqrt( leptonSF_err );

      leptonSF = leptonSF + (float)leptSyst_*leptonSF_err;

      eventWeight *= leptonSF;
      eventWeight *= fPUWeight->GetWeight(nPU);


    } // if is MC



//    // first: count them
//    ht = 0.;
//    njets=0;
//    for( unsigned iJet=0; iJet<nJets; ++iJet) {
//
//      if( ptJet[iJet] < 20. ) continue;
//      if( fabs(etaJet[iJet]) > etaJet_thresh_ ) continue;
//
//      ht += ptJet[iJet];
//
//      njets++;
//
//    }
//if( njets<3 ) continue;


    // define jets:
    
    int i_jetB1=-1;
    int i_jetB2=-1;


    float bestBtag=-9999.;
    float bestBtag2=-9999.;
    ht = 0.;
    njets=0;
    nBjets_loose = 0;
    nBjets_medium  = 0;

    AnalysisJet jetB1, jetB2;
    jetB1.SetPtEtaPhiE( 0., 0., 0., 0. );
    jetB2.SetPtEtaPhiE( 0., 0., 0., 0. );
  
    TLorentzVector pfMet_vector;
    pfMet_vector.SetPtEtaPhiE( pfMet, 0., phiMet, pfMet );

    std::vector< TLorentzVector > selectedJets;


    for( unsigned iJet=0; iJet<nJets; ++iJet) {

      // JES syst:
      float ptJet_corr = ptJet[iJet] + (float)jes_*ptUncertJet[iJet]*ptJet[iJet];

      // JER syst:
      if( jer_ ) {
        AnalysisJet thisJet_tmp;
        thisJet_tmp.SetPtEtaPhiE( ptJet_corr, etaJet[iJet], phiJet[iJet], eJet[iJet]);
        AnalysisJet thisGenJet_tmp;
        thisGenJet_tmp.SetPtEtaPhiE( ptGenJet[iJet], etaGenJet[iJet], phiGenJet[iJet], eGenJet[iJet]);
        float jer_SF = getJERSF( thisJet_tmp.Eta() );
        bool matched = (thisGenJet_tmp.Pt()>3.) ? thisJet_tmp.DeltaR(thisGenJet_tmp)<0.5 : false;
        if( matched )
          ptJet_corr = TMath::Max(0., thisGenJet_tmp.Pt()+jer_SF*(thisJet_tmp.Pt()-thisGenJet_tmp.Pt() ) );
        else
          ptJet_corr *= rand->Gaus( 1., fabs(1.-jer_SF) );
      }
      

      if( ptJet_corr < ptJet_thresh_ ) continue;
      if( fabs(etaJet[iJet]) > etaJet_thresh_ ) continue;

      ht += ptJet_corr;

      AnalysisJet thisJet;
      thisJet.SetPtEtaPhiE( ptJet_corr, etaJet[iJet], phiJet[iJet], eJet[iJet]);

      // MET syst:
      AnalysisJet thisJet_uncorr;
      thisJet_uncorr.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet]);

      pfMet_vector -= thisJet_uncorr;
      pfMet_vector += thisJet;
      

      thisJet.rmsCand = rmsCandJet[iJet];
      thisJet.ptD = ptDJet[iJet];
      thisJet.nCharged = nChargedJet[iJet];
      thisJet.nNeutral = nNeutralJet[iJet];
      thisJet.eMuons = eMuonsJet[iJet]/thisJet.Energy();
      thisJet.eElectrons = eElectronsJet[iJet]/thisJet.Energy();

      thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet[iJet];
      thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet[iJet];
      thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet[iJet];
      thisJet.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet[iJet];
      thisJet.jetProbabilityBJetTag               = jetProbabilityBJetTagJet[iJet];

      thisJet.ptGen = ptGenJet[iJet];
      thisJet.etaGen = etaGenJet[iJet];
      thisJet.phiGen = phiGenJet[iJet];
      thisJet.eGen = eGenJet[iJet];

      ////match to parton:
      //int partFlavor=0;
      //float deltaRmin=999.;
      //for(unsigned iPart=0; iPart<nPart; ++iPart ) {
      //  TLorentzVector thisPart;
      //  thisPart.SetPtEtaPhiE( ptPart[iPart], etaPart[iPart], phiPart[iPart], ePart[iPart] );
      //  float thisDeltaR = thisJet.DeltaR(thisPart);
      //  if( thisDeltaR<deltaRmin ) {
      //    partFlavor = pdgIdPart[iPart];
      //    deltaRmin = thisDeltaR;
      //  }
      //}

      if( isMC ) {
        TLorentzVector parton;
        parton.SetPtEtaPhiE( thisJet.Pt(), etaPartJet[iJet], phiPartJet[iJet], thisJet.Energy() );
        float deltaR = parton.DeltaR( thisJet );
        thisJet.pdgIdPart = (deltaR<0.5) ? pdgIdPartJet[iJet] : 21; //needed only for btag SF's
      } else {
        thisJet.pdgIdPart = 0;
      }


      float thisBtag;
      if( bTaggerType_=="TCHE" )
        thisBtag = thisJet.trackCountingHighEffBJetTag;
      else if( bTaggerType_=="SSVHE" ) 
        thisBtag = thisJet.simpleSecondaryVertexHighEffBJetTag;

      bool isBtagged_loose = ( thisBtag > this->get_btagThresh("loose") );
      bool isBtagged_medium = ( thisBtag > this->get_btagThresh("medium") );

      // take into account btag scale factors
      btsfutil->modifyBTagsWithSF_fast(isBtagged_loose, isBtagged_medium, thisJet.Pt(), thisJet.Eta(), thisJet.pdgIdPart, meanminmax );
//      btsfutil->modifyBTagsWithSF_fast(isBtagged_loose, isBtagged_medium, thisJet.Pt(), thisJet.Eta(), thisJet.pdgIdPart, "min" );

      njets += 1;
      if( isBtagged_loose ) nBjets_loose += 1;
      if( isBtagged_medium ) nBjets_medium += 1;

      if( thisBtag > bestBtag ) {
        bestBtag2 = bestBtag;
        bestBtag = thisBtag;
        jetB2 = jetB1;
        jetB1 = thisJet;
        i_jetB2 = i_jetB1;
        i_jetB1 = iJet;
      } else if( thisBtag > bestBtag2 ) {
        bestBtag2 = thisBtag;
        i_jetB2 = iJet;
        jetB2 = thisJet;
      }

      selectedJets.push_back(thisJet);


    } // for jets


    if( njets < njets_thresh_ ) continue;
    //if( nBjets_loose < nBjets_loose_thresh_ ) continue;
    //if( nBjets_medium < nBjets_medium_thresh_ ) continue;

    passed_btag = (nBjets_loose >= nBjets_loose_thresh_) && ( nBjets_medium >= nBjets_medium_thresh_ );

    if( jetB1.Pt()==0. || jetB2.Pt()==0. ) {
      std::cout << "JetB1/B2 are not defined. There must be a problem." << std::endl;
      exit(33);
    }

    AnalysisJet jet3, jet4;
    jet3.SetPtEtaPhiE( 0., 0., 0., 0. );
    jet4.SetPtEtaPhiE( 0., 0., 0., 0. );


    // takes care of JES syst:
    pfMet = pfMet_vector.Pt();
    
    // now add other jets ordered in pt:
    int istep=0;
    for( unsigned iJet=0; iJet<selectedJets.size(); ++iJet) {

      AnalysisJet thisJet = selectedJets[iJet];

      if( thisJet==jetB1 || thisJet==jetB2 ) continue;
      //if( iJet==i_jetB1 || iJet==i_jetB2 ) continue;
 
      //if( ptJet_corr < ptJet_thresh_ ) continue;
      //if( fabs(etaJet[iJet]) > etaJet_thresh_ ) continue;

      //AnalysisJet thisJet;
      //thisJet.SetPtEtaPhiE( ptJet_corr, etaJet[iJet], phiJet[iJet], eJet[iJet]);

      //thisJet.rmsCand = rmsCandJet[iJet];
      //thisJet.ptD = ptDJet[iJet];
      //thisJet.nCharged = nChargedJet[iJet];
      //thisJet.nNeutral = nNeutralJet[iJet];
      //thisJet.eMuons = eMuonsJet[iJet]/thisJet.Energy();
      //thisJet.eElectrons = eElectronsJet[iJet]/thisJet.Energy();

      //thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagJet[iJet];
      //thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagJet[iJet];
      //thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagJet[iJet];
      //thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagJet[iJet];
      //thisJet.jetBProbabilityBJetTag              = jetBProbabilityBJetTagJet[iJet];
      //thisJet.jetProbabilityBJetTag               = jetProbabilityBJetTagJet[iJet];

      //thisJet.ptGen = ptGenJet[iJet];
      //thisJet.etaGen = etaGenJet[iJet];
      //thisJet.phiGen = phiGenJet[iJet];
      //thisJet.eGen = eGenJet[iJet];

      //if( isMC ) {
      //  TLorentzVector parton;
      //  parton.SetPtEtaPhiE( thisJet.Pt(), etaPartJet[iJet], phiPartJet[iJet], thisJet.Energy() );
      //  float deltaR = parton.DeltaR( thisJet );
      //  thisJet.pdgIdPart = (deltaR<0.5) ? pdgIdPartJet[iJet] : 21; //needed only for btag SF's
      //} else {
      //  thisJet.pdgIdPart = 0;
      //}


      if( istep==0 ) {
        jet3 = thisJet;
        istep++;
      } else if( istep==1 ) {
        jet4 = thisJet;
        istep++;
      } else {
        break;
      }

    } //for additional jets




    
    TLorentzVector leptZ1, leptZ2;
    leptZ1.SetPtEtaPhiE( ptLeptZ1, etaLeptZ1, phiLeptZ1, eLeptZ1 );
    leptZ2.SetPtEtaPhiE( ptLeptZ2, etaLeptZ2, phiLeptZ2, eLeptZ2 );

    TLorentzVector diLepton = leptZ1+leptZ2;

    if( event==DEBUG_EVENTNUMBER ) {
      std::cout << std::endl << std::endl << "----------------------------------" << std::endl;
      std::cout << "** LOG FOR RUN: " << run << "   EVENT: " << DEBUG_EVENTNUMBER << std::endl << std::endl;
      std::cout << "eventWeight: " << eventWeight << std::endl;
      std::cout << "leptType: " << leptType << std::endl; 
      std::cout << "leptZ1.Pt(): " << leptZ1.Pt() << " leptZ1.Eta(): " << leptZ1.Eta() << std::endl;
      std::cout << "leptZ2.Pt(): " << leptZ2.Pt() << " leptZ2.Eta(): " << leptZ2.Eta() << std::endl;
      std::cout << "diLepton.M(): " << diLepton.M() << std::endl;
    }


    if( diLepton.M()<50. ) continue; // gen cut in DY sample




    if( leptType<=1 ) {
      h1_mZll_prepresel->Fill( diLepton.M(), eventWeight );
      if( leptType==0 )
        h1_mZll_prepresel_MU->Fill( diLepton.M(), eventWeight );
      else
        h1_mZll_prepresel_ELE->Fill( diLepton.M(), eventWeight );
    } else {  //opposite flavour leptons: ttbar control region

      if( event==DEBUG_EVENTNUMBER ) std::cout << "Filling OF plot and continuing." << std::endl;
      h1_mZll_OF_prepresel->Fill( diLepton.M(), eventWeight );
      if( leptType==2 )
        h1_mZll_OF2_prepresel->Fill( diLepton.M(), eventWeight );
      if( leptType==3 )
        h1_mZll_OF3_prepresel->Fill( diLepton.M(), eventWeight );
      continue;
    }



    // btag free region: Z+jets and WZ control region
    //if( nBjets_loose == 0 ) {
    if( nBjets_medium == 0 ) {
    //if( !passed_btag ) {
      h1_mZll_prepresel_antibtag->Fill( diLepton.M(), eventWeight );
      if( nLept>0 )
        h1_mZll_presel_antibtag->Fill( diLepton.M(), eventWeight );
    }

    



    // fill some histos before requiring third lepton:
    h1_nvertex_PUW->Fill(nvertex, eventWeight);

    h1_rhoPF_prepresel->Fill( rhoPF, eventWeight);
    h1_nJets_prepresel->Fill( njets, eventWeight );



    // this is the trilepton channel: require at least one other lepton:
    if( nLept<1 ) continue;



    h1_nJets_presel->Fill( njets, eventWeight );
    h1_nBJets_loose_presel->Fill( nBjets_loose, eventWeight );
    h1_nBJets_medium_presel->Fill( nBjets_medium, eventWeight );
    h1_pfMet_presel->Fill( pfMet, eventWeight);
    h1_rhoPF_presel->Fill( rhoPF, eventWeight);


    TLorentzVector lept3;
    lept3.SetPtEtaPhiE( ptLept[0], etaLept[0], phiLept[0], eLept[0] );

    h1_ptLept3_presel->Fill( lept3.Pt(), eventWeight );
    h1_etaLept3_presel->Fill( lept3.Pt(), eventWeight );
    h1_leptTypeLept3_presel->Fill( leptTypeLept[0], eventWeight );
    h1_combinedIsoRelLept3_presel->Fill( combinedIsoRelLept[0], eventWeight );


    if( lept3.Pt() < ptLept3_thresh_ ) continue;
    if( combinedIsoRelLept[0] > combinedIsoRelLept3_thresh_ ) continue;


 


    // ----------------------------
    // KINEMATIC SELECTION: LEPTONS
    // ----------------------------

    h1_mZll_presel->Fill( diLepton.M(), eventWeight );


    if( leptZ1.Pt() < ptLeptZ1_thresh_ ) continue;
    if( leptZ2.Pt() < ptLeptZ2_thresh_ ) continue;
    if( fabs(leptZ1.Eta()) > etaLeptZ1_thresh_ ) continue;
    if( fabs(leptZ2.Eta()) > etaLeptZ2_thresh_ ) continue;
    if( diLepton.Pt() < ptZll_thresh_ ) continue;
    isMZllSignalRegion=true;
    if( diLepton.M() < mZll_threshLo_ || diLepton.M() > mZll_threshHi_ ) isMZllSignalRegion=false;



      
    if( pfMet < met_thresh_ ) continue;  
      
      

      
      
    // ------------
    // AND NOW JETS
    // ------------

    //std::cout <<  jetB1.Pt() << std::endl;
    //std::cout <<  jetB2.Pt() << std::endl;
    //std::cout << jet3.Pt() << std::endl;
    //std::cout << jet4.Pt() << std::endl;
    //std::cout << std::endl << std::endl << std::endl;

   
    h1_nJets->Fill( njets , eventWeight );
    h1_nBJets_loose->Fill( nBjets_loose , eventWeight );
    h1_nBJets_medium->Fill( nBjets_medium , eventWeight );



    // define all kinematic variables:
    
    TLorentzVector met;
    met.SetPtEtaPhiE( pfMet, 0., phiMet, pfMet );

    TLorentzVector W = lept3 + met;
    // this is for tree_passedEvents:
    mt = W.Mt();


    h1_mTW->Fill( W.Mt() , eventWeight );

    TLorentzVector lZ = lept3 + diLepton;
    TLorentzVector lZ_plusMet = lZ + met;

    mT_lZmet = ( sqrt( lZ.Pt()*lZ.Pt() + lZ.M()*lZ.M() ) + pfMet )*( sqrt( lZ.Pt()*lZ.Pt() + lZ.M()*lZ.M() ) + pfMet )  -  lZ_plusMet.Pt()*lZ_plusMet.Pt();
    mT_lZmet = sqrt(mT_lZmet);

    ptJetB1_t = jetB1.Pt();
    ptJetB2_t = jetB2.Pt();

    TLorentzVector b1jj = jetB1 + jet3 + jet4;
    TLorentzVector b2jj = jetB2 + jet3 + jet4;

    TLorentzVector b1jjZ = jetB1 + jet3 + jet4 + diLepton;
    TLorentzVector b2jjZ = jetB2 + jet3 + jet4 + diLepton;

    mb1jj = b1jj.M();
    mb2jj = b2jj.M();

    if( fabs(mb1jj-tmass) < fabs(mb2jj-tmass) )
      mbjj_best = mb1jj;
    else
      mbjj_best = mb2jj;

    mbjj_max = TMath::Max(mb1jj,mb2jj);

    mb1jjZ = b1jjZ.M();
    mb2jjZ = b2jjZ.M();

    if( fabs(mb1jjZ-tmass) < fabs(mb2jjZ-tmass) )
      mbjjZ_best = mb1jjZ;
    else
      mbjjZ_best = mb2jjZ;

    mbjjZ_max = TMath::Max(mb1jjZ,mb2jjZ);


    TLorentzVector b1W = jetB1 + lept3 + met;
    TLorentzVector b2W = jetB2 + lept3 + met;

    TLorentzVector b1WZ = jetB1 + diLepton + lept3 + met;
    TLorentzVector b2WZ = jetB2 + diLepton + lept3 + met;

    mTb1W = b1W.Mt();
    mTb2W = b2W.Mt();

    if( fabs(mTb1W-tmass) < fabs(mTb2W-tmass) )
      mTbW_best = mTb1W;
    else
      mTbW_best = mTb2W;
    
    mTbW_max = TMath::Max(mTb1W,mTb2W);

    mTb1WZ = b1WZ.Mt();
    mTb2WZ = b2WZ.Mt();

    if( fabs(mTb1WZ-tmass) < fabs(mTb2WZ-tmass) )
      mTbWZ_best = mTb1WZ;
    else
      mTbWZ_best = mTb2WZ;

    mTbWZ_max = TMath::Max(mTb1WZ,mTb2WZ);



    if( ht<ht_thresh_ ) continue;
    if( mbjjZ_best<mbjjZ_best_thresh_ ) continue;



    if( passed_btag ) {

      // fill histograms:
      
      h1_mZll->Fill( diLepton.M(), eventWeight );


      if( isMZllSignalRegion ) {

        h1_pfMet->Fill( pfMet, eventWeight );
        h1_ht->Fill( ht, eventWeight );
        h1_metSignificance->Fill( metSignificance, eventWeight );

        h1_mT_lZmet->Fill( mT_lZmet, eventWeight );


        h1_ptLeptZ1->Fill( leptZ1.Pt(), eventWeight );
        h1_ptLeptZ2->Fill( leptZ2.Pt(), eventWeight );
        h1_etaLeptZ1->Fill( leptZ1.Eta(), eventWeight );
        h1_etaLeptZ2->Fill( leptZ2.Eta(), eventWeight );

        if( nLept>0 ) 
          h1_ptLept3->Fill( lept3.Pt(), eventWeight );

        h1_deltaRZll->Fill( leptZ2.DeltaR(leptZ2), eventWeight );

        h1_ptZll->Fill( diLepton.Pt(), eventWeight );
        h1_etaZll->Fill( diLepton.Eta(), eventWeight );


        h1_ptJetB1->Fill( ptJetB1_t, eventWeight );
        h1_ptJetB2->Fill( ptJetB2_t, eventWeight );
        h1_ptJet3->Fill( jet3.Pt(), eventWeight );
        h1_ptJet4->Fill( jet4.Pt(), eventWeight );

        h1_etaJetB1->Fill( jetB1.Eta(), eventWeight );
        h1_etaJetB2->Fill( jetB2.Eta(), eventWeight );
        h1_etaJet3->Fill( jet3.Eta(), eventWeight );
        if( jet4.Pt()>0. )
          h1_etaJet4->Fill( jet4.Eta(), eventWeight );


        h1_partFlavorJetB1->Fill( jetB1.pdgIdPart, eventWeight );
        h1_partFlavorJetB2->Fill( jetB2.pdgIdPart, eventWeight );
        h1_partFlavorJet3->Fill( jet3.pdgIdPart, eventWeight );
        h1_partFlavorJet4->Fill( jet4.pdgIdPart, eventWeight );


        float bTaggerJetB1, bTaggerJetB2;
        if( bTaggerType_=="TCHE" ) {
          bTaggerJetB1 = jetB1.trackCountingHighEffBJetTag;
          bTaggerJetB2 = jetB2.trackCountingHighEffBJetTag;
        } else if( bTaggerType_=="SSVHE" ) {
          bTaggerJetB1 = jetB1.simpleSecondaryVertexHighEffBJetTag;
          bTaggerJetB2 = jetB2.simpleSecondaryVertexHighEffBJetTag;
        }


        h1_bTagJetB1->Fill( bTaggerJetB1, eventWeight );
        h1_bTagJetB2->Fill( bTaggerJetB2, eventWeight );


        h1_deltaRbb->Fill( jetB1.DeltaR(jetB2), eventWeight );

        h1_mb1jj->Fill( mb1jj, eventWeight );
        h1_mb2jj->Fill( mb2jj, eventWeight );

        h1_mbjj_best->Fill( mbjj_best, eventWeight );
        h1_mbjj_max->Fill( mbjj_max, eventWeight );


        h1_mb1jjZ->Fill( mb1jjZ, eventWeight );
        h1_mb2jjZ->Fill( mb2jjZ, eventWeight );


        h1_mbjjZ_best->Fill( mbjjZ_best, eventWeight );
        h1_mbjjZ_max->Fill( mbjjZ_max, eventWeight );


        h1_mTb1W->Fill( mTb1W, eventWeight );
        h1_mTb2W->Fill( mTb2W, eventWeight );

        h1_mTbW_best->Fill( mTbW_best, eventWeight );
        h1_mTbW_max->Fill( mTbW_max, eventWeight );

        
        h1_mTb1WZ->Fill( mTb1WZ, eventWeight );
        h1_mTb2WZ->Fill( mTb2WZ, eventWeight );

        h1_mTbWZ_best->Fill( mTbWZ_best, eventWeight );
        h1_mTbWZ_max->Fill( mTbWZ_max, eventWeight );


        int leptType_forYield=-1;
        if( leptType==1 && leptTypeLept[0]==1 ) leptType_forYield = 0;
        if( leptType==1 && leptTypeLept[0]==0 ) leptType_forYield = 1;
        if( leptType==0 && leptTypeLept[0]==1 ) leptType_forYield = 2;
        if( leptType==0 && leptTypeLept[0]==0 ) leptType_forYield = 3;
        if( leptType_forYield<0 ) {
          std::cout << "There must be a problem, leptType_forYield can't be negative" << std::endl;
          exit(55);
        }
        h1_channelYields->Fill( leptType_forYield, eventWeight );
        h1_channelYields->Fill( 4, eventWeight ); //this is the total

      
      } // is mzll signal region

    }


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





    leptType3 = leptTypeLept[0];

    ptLeptZ1_t = leptZ1.Pt();
    ptLeptZ2_t = leptZ2.Pt();
    ptLept3_t = lept3.Pt();
    etaLeptZ1_t = leptZ1.Eta();
    etaLeptZ2_t = leptZ2.Eta();
    etaLept3_t = lept3.Eta();

    ptZll_t = diLepton.Pt();
    mZll_t = diLepton.M();

//  ptJetB1_t = jetB1.Pt();
//  ptJetB2_t = jetB2.Pt();
//  etaJetB1_t = jetB1.Eta();
//  etaJetB2_t = jetB2.Eta();

//  //bTagJetB1_t = bTaggerJetB1;
//  //bTagJetB2_t = bTaggerJetB2;

//  ptJet3_t = jet3.Pt();
//  ptJet4_t = jet4.Pt();
//  etaJet3_t = jet3.Eta();
//  etaJet4_t = jet4.Eta();



    // and fill tree (remember this includes the mZll sidebands):
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
  h1_pfMet_presel->Write();
  h1_ht->Write();

  h1_metSignificance->Write();


  h1_rhoPF_noPUW->Write();
  h1_rhoPF_prepresel->Write();
  h1_rhoPF_presel->Write();
  h1_rhoPF->Write();


  h1_ptLeptZ1->Write();
  h1_ptLeptZ2->Write();
  h1_etaLeptZ1->Write();
  h1_etaLeptZ2->Write();

  h1_ptLept3_presel->Write();
  h1_etaLept3_presel->Write();
  h1_leptTypeLept3_presel->Write();
  h1_combinedIsoRelLept3_presel->Write();

  h1_ptLept3->Write();
  h1_etaLept3->Write();

  h1_deltaRZll->Write();

  h1_ptZll->Write();
  h1_etaZll->Write();
  h1_mZll_prepresel->Write();
  h1_mZll_prepresel_MU->Write();
  h1_mZll_prepresel_ELE->Write();
  h1_mZll_prepresel_antibtag->Write();
  h1_mZll_OF_prepresel->Write();
  h1_mZll_OF2_prepresel->Write();
  h1_mZll_OF3_prepresel->Write();
  h1_mZll_presel->Write();
  h1_mZll_presel_antibtag->Write();
  h1_mZll->Write();


  h1_mTW->Write();
  h1_mT_lZmet->Write();

  h1_nJets_prepresel->Write();
  h1_nJets_presel->Write();
  h1_nJets->Write();

  h1_nBJets_loose_presel->Write();
  h1_nBJets_loose->Write();

  h1_nBJets_medium_presel->Write();
  h1_nBJets_medium->Write();


  h1_bTagJetB1->Write();
  h1_bTagJetB2->Write();


  h1_ptJetB1->Write();
  h1_ptJetB2->Write();
  h1_ptJet3->Write();
  h1_ptJet4->Write();

  h1_etaJetB1->Write();
  h1_etaJetB2->Write();
  h1_etaJet3->Write();
  h1_etaJet4->Write();

  h1_partFlavorJetB1->Write();
  h1_partFlavorJetB2->Write();
  h1_partFlavorJet3->Write();
  h1_partFlavorJet4->Write();


  h1_deltaRbb->Write();
  
  h1_mb1jj->Write();
  h1_mb2jj->Write();
  h1_mbjj_best->Write();
  h1_mbjj_max->Write();
  
  h1_mb1jjZ->Write();
  h1_mb2jjZ->Write();
  h1_mbjjZ_best->Write();
  h1_mbjjZ_max->Write();

  h1_mTb1W->Write();
  h1_mTb2W->Write();
  h1_mTbW_best->Write();
  h1_mTbW_max->Write();

  h1_mTb1WZ->Write();
  h1_mTb2WZ->Write();
  h1_mTbWZ_best->Write();
  h1_mTbWZ_max->Write();

  h1_channelYields->Write();

  outFile_->Close();


} // finalize()



void Ntp1Finalizer_TTZTrilepton::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;

  if( selectionType_=="presel" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 0;

    mZll_threshLo_ = 70.;
    mZll_threshHi_ = 110.;

    ptZll_thresh_ = 0.;

    met_thresh_ = 0.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 0.;
 
  } else if( selectionType_=="preselplus" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 0.;

    met_thresh_ = 30.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 0.;
 
  } else if( selectionType_=="sel1" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 20.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 4;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 0.;

    met_thresh_ = 30.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 0.;

  } else if( selectionType_=="sel2" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 20.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 4;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 0.;

    met_thresh_ = 0.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 0.;

  } else if( selectionType_=="sel3" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 20.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 0.;

    met_thresh_ = 0.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 0.;

  } else if( selectionType_=="optsel1" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 23.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 28.;

    met_thresh_ = 30.;
    ht_thresh_ = 184.;
    mbjjZ_best_thresh_ = 0.;

  } else if( selectionType_=="optsel2" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 30.;

    met_thresh_ = 30.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 200.;

  } else if( selectionType_=="optsel2_noHT" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 0;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 30.;

    met_thresh_ = 30.;
    ht_thresh_ = 0.;
    mbjjZ_best_thresh_ = 0.;

  } else if( selectionType_=="optsel3" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 2;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 35.;

    met_thresh_ = 0.;
    ht_thresh_ = 120.;
    mbjjZ_best_thresh_ = 0.;

  } else if( selectionType_=="presel_2LB_1MB" ) {

    ptLeptZ1_thresh_ = 20.;
    ptLeptZ2_thresh_ = 20.;
    ptLept3_thresh_ = 10.;
    etaLeptZ1_thresh_ = 3.;
    etaLeptZ2_thresh_ = 3.;
    etaLept3_thresh_ = 3.;

    combinedIsoRelLept3_thresh_ = 1.;

    ptJet_thresh_ = 20.;
    etaJet_thresh_ = 2.4;

    njets_thresh_ = 3;
    nBjets_loose_thresh_ = 2;
    nBjets_medium_thresh_ = 1;

    mZll_threshLo_ = 81.;
    mZll_threshHi_ = 101.;

    ptZll_thresh_ = 35.;

    met_thresh_ = 0.;
    ht_thresh_ = 120.;
    mbjjZ_best_thresh_ = 0.;


  } else {

    std::cout << "Unknown selection type '" << selectionType << "'. Exiting." << std::endl;
    exit(1112);

  }

  
} //setSelectionType





float Ntp1Finalizer_TTZTrilepton::get_btagThresh( const std::string& btag_OP_ ) {

  if( btag_OP_ == "none" ) return -9999.;

  else if( btag_OP_ == "loose" ) {

    if( bTaggerType_=="TCHE" ) return  1.7;
    else if( bTaggerType_=="SSVHE" ) return  -9999.;

  } else if( btag_OP_ == "medium" ) {

    if( bTaggerType_=="TCHE" ) return  3.3;
    else if( bTaggerType_=="SSVHE" ) return  1.74;

  } else {

    std::cout << "Didn't find btag OP. Returning -9999." << std::endl;
  
  }

  return -9999.;

}




ValueAndError getMuonHLTSF_DoubleTrigger( float pt, float eta, const std::string& runPeriod ) {

  float hltsf = 0.;
  float hltsf_err = 0.;

  // these numbers taken from AN2011-399-v4
  if( runPeriod=="Run2011A" ) {

    if( fabs(eta)<0.8 ) {
      hltsf = 0.975;
      hltsf_err = 0.004;
    } else if( fabs(eta)<2.1 ) {
      hltsf = 0.950;
      hltsf_err = 0.005;
    } else  {
      hltsf = 0.910;
      hltsf_err = 0.01;
    }

  } else if( runPeriod=="Run2011B" ) {

    if( fabs(eta)<0.8 ) {
      if( pt<40. ) {
        hltsf = 0.977;
        hltsf_err = 0.001;
      } else {
        hltsf = 0.975;
        hltsf_err = 0.001;
      }
    } else if( fabs(eta)<2.1 ){
      if( pt<40. ) {
        hltsf = 0.955;
        hltsf_err = 0.002;
      } else {
        hltsf = 0.955;
        hltsf_err = 0.001;
      }
    } else { // eta 2.1 -> 2.4
      if( pt<40. ) {
        hltsf = 0.89;
        hltsf_err = 0.007;
      } else {
        hltsf = 0.90;
        hltsf_err = 0.006;
      }
    }

  } else {

    std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;

  }

  ValueAndError ve_hlt;
  ve_hlt.val = hltsf;
  ve_hlt.err = hltsf_err;

  return ve_hlt;

}




ValueAndError getMuonHLTSF_SingleTrigger( float pt, float eta, const std::string& runPeriod ) {

  if( pt<25. ) { 
    ValueAndError ve_hlt;
    ve_hlt.val = 0.;
    ve_hlt.err = 0.;
    return ve_hlt;
  }

  float hltsf = 0.;
  float hltsf_err = 0.;

  if( runPeriod=="Run2011A1" ) { //up to may10 technical stop

    if( fabs(eta)<0.8 ) {
      hltsf = 0.896;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.807;
      hltsf_err = 0.001;
    } else {
      hltsf = 0.608;
      hltsf_err = 0.001;
    }

  } else if( runPeriod=="Run2011A2" ) { //from may10 to EPS

    if( fabs(eta)<0.8 ) {
      hltsf = 0.895;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.838;
      hltsf_err = 0.001;
    } else {
      hltsf = 0.738;
      hltsf_err = 0.001;
    }

  } else if( runPeriod=="Run2011A3" ) { //from EPS to end of Run2011A

    if( fabs(eta)<0.8 ) {
      hltsf = 0.890;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.809;
      hltsf_err = 0.001;
    } else {
      hltsf = 0.493;
      hltsf_err = 0.001;
    }

  } else if( runPeriod=="Run2011B" ) {

    if( fabs(eta)<0.8 ) {
      hltsf = 0.87;
      hltsf_err = 0.001;
    } else if( fabs(eta)<2.1 ){
      hltsf = 0.79;
      hltsf_err = 0.001;
    } else  { //using HLT_IsoMu24_eta2p1
      hltsf = 0.;
      hltsf_err = 0.;
    }

  } else {

    std::cout << "WARNING! Unknown run period: " << runPeriod << "! Returning HLTSF=0." << std::endl;

  }

  ValueAndError ve_hlt;
  ve_hlt.val = hltsf;
  ve_hlt.err = hltsf_err;

  return ve_hlt;


}


ValueAndError getEventHLTSF( ValueAndError effSingle1, ValueAndError effSingle2, ValueAndError effDouble1, ValueAndError effDouble2 ) {

  float HLTSF = effDouble1.val * effDouble2.val +
                effSingle2.val * (1. - effDouble2.val ) +
                effSingle1.val * (1. - effDouble1.val );

  float HLTSF_err = effDouble1.err * effDouble1.err* effDouble2.val * effDouble2.val +
                    effDouble1.val * effDouble1.val* effDouble2.err * effDouble2.err +
                    effSingle2.err * effSingle2.err * (1. - effDouble2.val ) * (1. - effDouble2.val ) +
                    effSingle2.val * effSingle2.val * effDouble2.err * effDouble2.err +
                    effSingle1.err * effSingle1.err * (1. - effDouble1.val ) * (1. - effDouble1.val ) +
                    effSingle1.val * effSingle1.val * effDouble1.err * effDouble1.err;

  HLTSF_err = sqrt( HLTSF_err );

  ValueAndError ve_hlt;
  ve_hlt.val = HLTSF;
  ve_hlt.err = HLTSF_err;

  return ve_hlt;

}



ValueAndError getMuonRecoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( fabs(eta)<1.2 ) {
    recoSF = 0.996;
    recoSF_err = 0.001;
  } else {
    recoSF = 0.986;
    recoSF_err = 0.001;
  }

  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}


ValueAndError getElectronRecoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( fabs(eta)<0.8 ) {
    recoSF = 0.999;
    recoSF_err = 0.005;
  } else if( fabs(eta)<1.44 ) {
    recoSF = 0.964;
    recoSF_err = 0.003;
  } else if( fabs(eta)<1.57 ) {
    recoSF = 0.99;
    recoSF_err = 0.04;
  } else if( fabs(eta)<2.0 ) {
    recoSF = 0.992;
    recoSF_err = 0.006;
  } else {
    recoSF = 1.001;
    recoSF_err = 0.006;
  }

  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}



ValueAndError getMuonIsoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( pt<40. ) {

    if( fabs(eta)<0.9 ) {
      recoSF = 0.987;
      recoSF_err = 0.006;
    } else {
      recoSF = 0.995;
      recoSF_err = 0.005;
    }
   
  } else {

    if( fabs(eta)<0.9 ) {
      recoSF = 0.994;
      recoSF_err = 0.002;
    } else {
      recoSF = 0.996;
      recoSF_err = 0.002;
    }
   
  }


  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}


ValueAndError getElectronIsoSF( float pt, float eta ) {

  float recoSF;
  float recoSF_err;

  if( pt<40. ) {

    if( fabs(eta)<1.5 ) {
      recoSF = 0.988;
      recoSF_err = 0.006;
    } else {
      recoSF = 0.998;
      recoSF_err = 0.011;
    }
   
  } else {

    if( fabs(eta)<1.5 ) {
      recoSF = 0.988;
      recoSF_err = 0.003;
    } else {
      recoSF = 1.016;
      recoSF_err = 0.064;
    }
   
  }


  ValueAndError ve_reco;
  ve_reco.val = recoSF;
  ve_reco.err = recoSF_err;

  return ve_reco;

}





float getJERSF( float eta ) {

  // scale factors taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
 
  float SF=-1.;

  if( fabs(eta)<0.5 ) {
    SF = 1.052;
  } else if( fabs(eta)<1.1 ) {
    SF = 1.057;
  } else if( fabs(eta)<1.7 ) {
    SF = 1.096;
  } else if( fabs(eta)<2.3 ) {
    SF = 1.134;
  } else {
    SF = 1.288;
  }

  if( SF<0. ) {
    std::cout << "JERSF is negative, this shouldn't be possible." << std::endl;
    exit(1123);
  }

  return SF;

}


