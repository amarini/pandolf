#include <iostream>
#include <vector>
#include <cmath>

#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "TQObject.h"
//#include "RQ_OBJECT.h"

#include "fitTools.C"


struct photonidegcuts {

  float hovereiso;
  float hcaliso_rel;
  float hcaliso_abs;
  float ecaliso_rel;
  float ecaliso_abs;
  float trackiso_rel;
  float trackiso_abs;
  float setaetaEB;
  float setaetaEE;

};


class AnalysisPhoton : public TLorentzVector {

//RQ_OBJECT("AnalysisPhoton")
ClassDef (AnalysisPhoton, 0) 

public:

  AnalysisPhoton( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {};


  Float_t hcalIso;
  Float_t ecalIso;
  Int_t   nTrkIso;
  Float_t ptTrkIso;
  Float_t clusterMaj;
  Float_t clusterMin;
  Int_t hasPixelSeed;

  Float_t pid_hlwTrack;
  Float_t pid_twrHCAL;
  Float_t pid_HoverE;
  Float_t pid_jurECAL;
  Float_t pid_etawid;

};


bool DEBUG_ = false;

TChain* tree;


void addInput(const std::string& dataset);

bool passedPhotonID(const AnalysisPhoton & photon, photonidegcuts const& pid, vector<bool> *vpass);

std::vector<TH1D*> allocateHistoVector( int vectorSize, const std::string& name, int nBins, float xMin, float xMax );
float getThresh_constEff( TTree* tree, float eff, int iBinRho, float threshAbs);



void finalize(const std::string& dataset, std::string photonID="medium" ) {



  tree = new TChain("jetTree");



  std::string infileName, treeName;




  if( dataset=="G_TuneZ2_7TeV_pythia6" ) {

    addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
  //addInput( "G_Pt_800to1400_TuneZ2_7TeV_pythia6" );
  //addInput( "G_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="G_TuneZ2_7TeV_pythia6_CORR" ) {

    addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6_CORR" );

  } else if( dataset=="QCD_TuneZ2_7TeV_pythia6" ) {

    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="G_QCD_TuneZ2_7TeV_pythia6" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="DATA_EG_35X" ) {

    addInput( "MinimumBias_Commissioning10_May6thPDSkim2_SD_EG" );
    addInput( "MinimumBias_Commissioning10_SD_EG-v9" );
    addInput( "EG_Run2010A-PromptReco-v1" );
    addInput( "EG_Run2010A-PromptReco-v2" );

  } else if( dataset=="DATA_EG" ) {

    addInput( "MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1" );
    addInput( "EG_Run2010A-Jun14thReReco_v1" );
    addInput( "EG_Run2010A-PromptReco-v4" );
    addInput( "EG_Run2010A_PromptReco_v4_139347_139375" );
    addInput( "EG_Run2010A_PromptReco_v4_139376_139459" );

  } else if( dataset=="DATA_EG_37X" ) {

    addInput( "EG_Run2010A_Jul15thReReco_v1" );
    addInput( "EG_Run2010A_Jul26thReReco_v1" );

  } else if( dataset=="EG_3pb" ) {

    addInput( "EG_Run2010A-PromptReco-v4" );
    addInput( "EG_Run2010A-PromptReco-v4_runs143337-144114" );

  } else if( dataset=="DATA_Nov4ReReco" ) {

    addInput( "EG_Run2010A-Nov4ReReco_v1" );
    addInput( "Photon_Run2010B_Nov4ReReco" );

  } else {

    addInput( dataset );

  }




  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Float_t ptMax = 100.;
  int nBins_photBinning = ptPhot_binning.size()-1;
  Double_t ptPhot_binning_array[100];
  for( unsigned i=0; i<nBins_photBinning; ++i)
    ptPhot_binning_array[i] = ptPhot_binning[i];



  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Float_t eMet;
  tree->SetBranchAddress("epfMet", &eMet);
  Float_t phiMet;
  tree->SetBranchAddress("phipfMet", &phiMet);


  Float_t rhoCalo;
  tree->SetBranchAddress("rhoCalo", &rhoCalo);
  Float_t rhoPF;
  tree->SetBranchAddress("rhoPF", &rhoPF);

  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t ePhotReco;
  tree->SetBranchAddress("ePhotReco", &ePhotReco);
  Float_t ptPhotReco;
  tree->SetBranchAddress("ptPhotReco", &ptPhotReco);
  Float_t etaPhotReco;
  tree->SetBranchAddress("etaPhotReco", &etaPhotReco);
  Float_t phiPhotReco;
  tree->SetBranchAddress("phiPhotReco", &phiPhotReco);

  Float_t ePhotGen;
  tree->SetBranchAddress("ePhotGen", &ePhotGen);
  Float_t ptPhotGen;
  tree->SetBranchAddress("ptPhotGen", &ptPhotGen);
  Float_t etaPhotGen;
  tree->SetBranchAddress("etaPhotGen", &etaPhotGen);
  Float_t phiPhotGen;
  tree->SetBranchAddress("phiPhotGen", &phiPhotGen);

  Float_t hcalIsoPhotReco;
  tree->SetBranchAddress("hcalIsoPhotReco", &hcalIsoPhotReco);
  Float_t ecalIsoPhotReco;
  tree->SetBranchAddress("ecalIsoPhotReco", &ecalIsoPhotReco);
  Int_t nTrkIsoPhotReco;
  tree->SetBranchAddress("nTrkIsoPhotReco", &nTrkIsoPhotReco);
  Float_t ptTrkIsoPhotReco;
  tree->SetBranchAddress("ptTrkIsoPhotReco", &ptTrkIsoPhotReco);
  Float_t clusterMajPhotReco;
  tree->SetBranchAddress("clusterMajPhotReco", &clusterMajPhotReco);
  Float_t clusterMinPhotReco;
  tree->SetBranchAddress("clusterMinPhotReco", &clusterMinPhotReco);
  Float_t hasPixelSeedPhotReco;
  tree->SetBranchAddress("hasPixelSeedPhotReco", &hasPixelSeedPhotReco);

  Float_t pid_towerHCALPhotReco;
  tree->SetBranchAddress("pid_towerHCALPhotReco", &pid_towerHCALPhotReco);
  Float_t pid_HoverEPhotReco;
  tree->SetBranchAddress("pid_HoverEPhotReco", &pid_HoverEPhotReco);
  Float_t pid_jurECALPhotReco;
  tree->SetBranchAddress("pid_jurECALPhotReco", &pid_jurECALPhotReco);
  Float_t pid_hlwTrackPhotReco;
  tree->SetBranchAddress("pid_hlwTrackPhotReco", &pid_hlwTrackPhotReco);
  Float_t pid_sIEtaIEtaPhotReco;
  tree->SetBranchAddress("pid_sIEtaIEtaPhotReco", &pid_sIEtaIEtaPhotReco);

  Bool_t matchedToMC;
  tree->SetBranchAddress("matchedToMC", &matchedToMC);

  Float_t eJetReco;
  tree->SetBranchAddress("eJetReco", &eJetReco);
  Float_t ptJetReco;
  tree->SetBranchAddress("ptJetReco", &ptJetReco);
  Float_t ptCorrJetReco;
  tree->SetBranchAddress("ptCorrJetReco", &ptCorrJetReco);
  Float_t etaJetReco;
  tree->SetBranchAddress("etaJetReco", &etaJetReco);
  Float_t phiJetReco;
  tree->SetBranchAddress("phiJetReco", &phiJetReco);
  Float_t eTracksReco;
  tree->SetBranchAddress("eTracksReco", &eTracksReco);

  Float_t eJetGen;
  tree->SetBranchAddress("eJetGen", &eJetGen);
  Float_t ptJetGen;
  tree->SetBranchAddress("ptJetGen", &ptJetGen);
  Float_t etaJetGen;
  tree->SetBranchAddress("etaJetGen", &etaJetGen);
  Float_t phiJetGen;
  tree->SetBranchAddress("phiJetGen", &phiJetGen);
  Float_t eTracksGen;
  tree->SetBranchAddress("eTracksGen", &eTracksGen);

  Float_t ptPart;
  tree->SetBranchAddress("ptPart", &ptPart);
  Float_t etaPart;
  tree->SetBranchAddress("etaPart", &etaPart);
  Float_t phiPart;
  tree->SetBranchAddress("phiPart", &phiPart);
  Int_t pdgIdPart;
  tree->SetBranchAddress("pdgIdPart", &pdgIdPart);

  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t ptCorr2ndJetReco;
  tree->SetBranchAddress("ptCorr2ndJetReco", &ptCorr2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);

  Float_t pt2ndJetGen;
  tree->SetBranchAddress("pt2ndJetGen", &pt2ndJetGen);
  Float_t eta2ndJetGen;
  tree->SetBranchAddress("eta2ndJetGen", &eta2ndJetGen);
  Float_t phi2ndJetGen;
  tree->SetBranchAddress("phi2ndJetGen", &phi2ndJetGen);

  Float_t ptSecondaryJetsReco;
  tree->SetBranchAddress("ptSecondaryJetsReco", &ptSecondaryJetsReco);
  Float_t ptSecondaryJetsGen;
  tree->SetBranchAddress("ptSecondaryJetsGen", &ptSecondaryJetsGen);

  Bool_t passed_Photon10;
  tree->SetBranchAddress("passed_Photon10", &passed_Photon10);
  Bool_t passed_Photon15;
  tree->SetBranchAddress("passed_Photon15", &passed_Photon15);
  Bool_t passed_Photon20;
  tree->SetBranchAddress("passed_Photon20", &passed_Photon20);
  Bool_t passed_Photon25;
  tree->SetBranchAddress("passed_Photon25", &passed_Photon25);
  Bool_t passed_Photon30;
  tree->SetBranchAddress("passed_Photon30", &passed_Photon30);
  Bool_t passed_Photon35;
  tree->SetBranchAddress("passed_Photon35", &passed_Photon35);


  float pid_hcalIso;
  float pid_ecalIso;
  float pid_hlwTrackIso;

  int bin_rhoPF;
  int bin_rhoCalo;

  int nBinsRho = 10;
  float rhoCalo_max = 6.;
  float rhoPF_max = 16.;


  TProfile* hp_hcalIsoMean_vs_rhoCalo = new TProfile("hcalIsoMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_hcalIsoMean_vs_rhoCalo->Sumw2();
  TProfile* hp_hcalIsoPIDMean_vs_rhoCalo = new TProfile("hcalIsoPIDMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_hcalIsoPIDMean_vs_rhoCalo->Sumw2();
  TProfile* hp_hcalIsoMean_vs_rhoPF = new TProfile("hcalIsoMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_hcalIsoMean_vs_rhoPF->Sumw2();
  TProfile* hp_hcalIsoPIDMean_vs_rhoPF = new TProfile("hcalIsoPIDMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_hcalIsoPIDMean_vs_rhoPF->Sumw2();

  TProfile* hp_ecalIsoMean_vs_rhoCalo = new TProfile("ecalIsoMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_ecalIsoMean_vs_rhoCalo->Sumw2();
  TProfile* hp_ecalIsoPIDMean_vs_rhoCalo = new TProfile("ecalIsoPIDMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_ecalIsoPIDMean_vs_rhoCalo->Sumw2();
  TProfile* hp_ecalIsoMean_vs_rhoPF = new TProfile("ecalIsoMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_ecalIsoMean_vs_rhoPF->Sumw2();
  TProfile* hp_ecalIsoPIDMean_vs_rhoPF = new TProfile("ecalIsoPIDMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_ecalIsoPIDMean_vs_rhoPF->Sumw2();

  TProfile* hp_ptTrkIsoMean_vs_rhoCalo = new TProfile("ptTrkIsoMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_ptTrkIsoMean_vs_rhoCalo->Sumw2();
  TProfile* hp_hlwTrackIsoMean_vs_rhoCalo = new TProfile("pid_hlwTrackMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_hlwTrackIsoMean_vs_rhoCalo->Sumw2();
  TProfile* hp_ptTrkIsoMean_vs_rhoPF = new TProfile("ptTrkIsoMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_ptTrkIsoMean_vs_rhoPF->Sumw2();
  TProfile* hp_hlwTrackIsoMean_vs_rhoPF = new TProfile("pid_hlwTrackMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_hlwTrackIsoMean_vs_rhoPF->Sumw2();


  int nBinsDistrib = 200;

  std::vector<TH1D*> vh1_hcalIsoPhotReco_rhoPF  = allocateHistoVector( nBinsRho, "hcalIsoPhotReco_rhoPF", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_hcalIsoPID_rhoPF       = allocateHistoVector( nBinsRho, "hcalIsoPID_rhoPF", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_ecalIsoPhotReco_rhoPF  = allocateHistoVector( nBinsRho, "ecalIsoPhotReco_rhoPF", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_ecalIsoPID_rhoPF       = allocateHistoVector( nBinsRho, "ecalIsoPID_rhoPF", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_ptTrkIsoPhotReco_rhoPF = allocateHistoVector( nBinsRho, "ptTrkIsoPhotReco_rhoPF", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_hlwTrackIso_rhoPF      = allocateHistoVector( nBinsRho, "hlwTrackIso_rhoPF", nBinsDistrib, 0., 1.);

  std::vector<TH1D*> vh1_hcalIsoPhotReco_rhoCalo  = allocateHistoVector( nBinsRho, "hcalIsoPhotReco_rhoCalo", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_hcalIsoPID_rhoCalo       = allocateHistoVector( nBinsRho, "hcalIsoPID_rhoCalo", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_ecalIsoPhotReco_rhoCalo  = allocateHistoVector( nBinsRho, "ecalIsoPhotReco_rhoCalo", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_ecalIsoPID_rhoCalo       = allocateHistoVector( nBinsRho, "ecalIsoPID_rhoCalo", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_ptTrkIsoPhotReco_rhoCalo = allocateHistoVector( nBinsRho, "ptTrkIsoPhotReco_rhoCalo", nBinsDistrib, 0., 1.);
  std::vector<TH1D*> vh1_hlwTrackIso_rhoCalo      = allocateHistoVector( nBinsRho, "hlwTrackIso_rhoCalo", nBinsDistrib, 0., 1.);


  TTree* tree_hcalIsoPhotReco_rhoPF  = new TTree( "tree_hcalIsoPhotReco_rhoPF", "" );
  TTree* tree_hcalIsoPID_rhoPF       = new TTree( "tree_hcalIsoPID_rhoPF", "" );
  TTree* tree_ecalIsoPhotReco_rhoPF  = new TTree( "tree_ecalIsoPhotReco_rhoPF", "" );
  TTree* tree_ecalIsoPID_rhoPF       = new TTree( "tree_ecalIsoPID_rhoPF", "" );
  TTree* tree_ptTrkIsoPhotReco_rhoPF = new TTree( "tree_ptTrkIsoPhotReco_rhoPF", "" );
  TTree* tree_hlwTrackIso_rhoPF      = new TTree( "tree_hlwTrackIso_rhoPF", "" );

  TTree* tree_hcalIsoPhotReco_rhoCalo  = new TTree( "tree_hcalIsoPhotReco_rhoCalo", "" );
  TTree* tree_hcalIsoPID_rhoCalo       = new TTree( "tree_hcalIsoPID_rhoCalo", "" );
  TTree* tree_ecalIsoPhotReco_rhoCalo  = new TTree( "tree_ecalIsoPhotReco_rhoCalo", "" );
  TTree* tree_ecalIsoPID_rhoCalo       = new TTree( "tree_ecalIsoPID_rhoCalo", "" );
  TTree* tree_ptTrkIsoPhotReco_rhoCalo = new TTree( "tree_ptTrkIsoPhotReco_rhoCalo", "" );
  TTree* tree_hlwTrackIso_rhoCalo      = new TTree( "tree_hlwTrackIso_rhoCalo", "" );

  tree_hcalIsoPhotReco_rhoPF->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_hcalIsoPID_rhoPF->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_ecalIsoPhotReco_rhoPF->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_ecalIsoPID_rhoPF->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_ptTrkIsoPhotReco_rhoPF->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_hlwTrackIso_rhoPF->Branch( "eventWeight", &eventWeight, "eventWeight/F");

  tree_hcalIsoPhotReco_rhoCalo->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_hcalIsoPID_rhoCalo->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_ecalIsoPhotReco_rhoCalo->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_ecalIsoPID_rhoCalo->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_ptTrkIsoPhotReco_rhoCalo->Branch( "eventWeight", &eventWeight, "eventWeight/F");
  tree_hlwTrackIso_rhoCalo->Branch( "eventWeight", &eventWeight, "eventWeight/F");

  tree_hcalIsoPhotReco_rhoPF->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F"); 
  tree_hcalIsoPID_rhoPF->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F");
  tree_ecalIsoPhotReco_rhoPF->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F"); 
  tree_ecalIsoPID_rhoPF->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F");
  tree_ptTrkIsoPhotReco_rhoPF->Branch( "ePhotReco", &ptPhotReco, "ptPhotReco/F"); //ptPhot!!
  tree_hlwTrackIso_rhoPF->Branch( "ePhotReco", &ptPhotReco, "ptPhotReco/F"); //ptPhot!!

  tree_hcalIsoPhotReco_rhoCalo->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F"); 
  tree_hcalIsoPID_rhoCalo->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F");
  tree_ecalIsoPhotReco_rhoCalo->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F"); 
  tree_ecalIsoPID_rhoCalo->Branch( "ePhotReco", &ePhotReco, "ePhotReco/F");
  tree_ptTrkIsoPhotReco_rhoCalo->Branch( "ePhotReco", &ptPhotReco, "ptPhotReco/F"); //ptPhot!!
  tree_hlwTrackIso_rhoCalo->Branch( "ePhotReco", &ptPhotReco, "ptPhotReco/F"); //ptPhot!!

  tree_hcalIsoPhotReco_rhoPF->Branch( "bin_rho", &bin_rhoPF, "bin_rhoPF/I");
  tree_hcalIsoPID_rhoPF->Branch( "bin_rho", &bin_rhoPF, "bin_rhoPF/I");
  tree_ecalIsoPhotReco_rhoPF->Branch( "bin_rho", &bin_rhoPF, "bin_rhoPF/I");
  tree_ecalIsoPID_rhoPF->Branch( "bin_rho", &bin_rhoPF, "bin_rhoPF/I");
  tree_ptTrkIsoPhotReco_rhoPF->Branch( "bin_rho", &bin_rhoPF, "bin_rhoPF/I");
  tree_hlwTrackIso_rhoPF->Branch( "bin_rho", &bin_rhoPF, "bin_rhoPF/I");

  tree_hcalIsoPhotReco_rhoCalo->Branch( "bin_rho", &bin_rhoCalo, "bin_rhoCalo/I");
  tree_hcalIsoPID_rhoCalo->Branch( "bin_rho", &bin_rhoCalo, "bin_rhoCalo/I");
  tree_ecalIsoPhotReco_rhoCalo->Branch( "bin_rho", &bin_rhoCalo, "bin_rhoCalo/I");
  tree_ecalIsoPID_rhoCalo->Branch( "bin_rho", &bin_rhoCalo, "bin_rhoCalo/I");
  tree_ptTrkIsoPhotReco_rhoCalo->Branch( "bin_rho", &bin_rhoCalo, "bin_rhoCalo/I");
  tree_hlwTrackIso_rhoCalo->Branch( "bin_rho", &bin_rhoCalo, "bin_rhoCalo/I");

  tree_hcalIsoPhotReco_rhoPF->Branch( "isoVar", &hcalIsoPhotReco, "hcalIsoPhotReco/F");
  tree_hcalIsoPID_rhoPF->Branch( "isoVar", &pid_hcalIso, "pid_hcalIso/F");
  tree_ecalIsoPhotReco_rhoPF->Branch( "isoVar", &ecalIsoPhotReco, "ecalIsoPhotReco/F");
  tree_ecalIsoPID_rhoPF->Branch( "isoVar", &pid_ecalIso, "pid_ecalIso/F");
  tree_ptTrkIsoPhotReco_rhoPF->Branch( "isoVar", &ptTrkIsoPhotReco, "ptTrkIsoPhotReco/F");
  tree_hlwTrackIso_rhoPF->Branch( "isoVar", &pid_hlwTrackIso, "pid_hlwTrackIso/F");

  tree_hcalIsoPhotReco_rhoCalo->Branch( "isoVar", &hcalIsoPhotReco, "hcalIsoPhotReco/F");
  tree_hcalIsoPID_rhoCalo->Branch( "isoVar", &pid_hcalIso, "pid_hcalIso/F");
  tree_ecalIsoPhotReco_rhoCalo->Branch( "isoVar", &ecalIsoPhotReco, "ecalIsoPhotReco/F");
  tree_ecalIsoPID_rhoCalo->Branch( "isoVar", &pid_ecalIso, "pid_ecalIso/F");
  tree_ptTrkIsoPhotReco_rhoCalo->Branch( "isoVar", &ptTrkIsoPhotReco, "ptTrkIsoPhotReco/F");
  tree_hlwTrackIso_rhoCalo->Branch( "isoVar", &pid_hlwTrackIso, "pid_hlwTrackIso/F");




  float denom_hcalIso = 0.;
  float denom_pid_hcalIso = 0.;
  float denom_ecalIso = 0.;
  float denom_pid_ecalIso = 0.;
  float denom_ptTrkIso = 0.;
  float denom_pid_hlwTrackIso = 0.;

  float num_hcalIso = 0.;
  float num_pid_hcalIso = 0.;
  float num_ecalIso = 0.;
  float num_pid_ecalIso = 0.;
  float num_ptTrkIso = 0.;
  float num_pid_hlwTrackIso = 0.;


  photonidegcuts hggtightid;
  hggtightid.hovereiso=           0.02;
  hggtightid.hcaliso_rel=         0.0025;
  hggtightid.hcaliso_abs=         2.;
  hggtightid.ecaliso_rel=         0.006;
  hggtightid.ecaliso_abs=         2.;
  hggtightid.trackiso_rel=        0.001;
  hggtightid.trackiso_abs=        1.5;
  hggtightid.setaetaEB=           0.010;
  hggtightid.setaetaEE=           0.028;

  float hcalIsoPhotReco_threshRel = 0.05;
  float hcalIsoPhotReco_threshAbs = 2.4;
  float ecalIsoPhotReco_threshRel = 0.05;
  float ecalIsoPhotReco_threshAbs = 3.;
  float ptTrkIsoPhotReco_threshRel = 0.10;
    
//float pid_hcalIso_threshRel = hggtightid.hcaliso_rel;
//float pid_hcalIso_threshAbs = hggtightid.hcaliso_abs;
//float pid_ecalIso_threshRel = hggtightid.ecaliso_rel;
//float pid_ecalIso_threshAbs = hggtightid.ecaliso_abs;
//float pid_hlwTrackIso_threshRel = hggtightid.trackiso_rel;
//float pid_hlwTrackIso_threshAbs = hggtightid.trackiso_abs;
    

  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;


    AnalysisPhoton photon;
    photon.SetPtEtaPhiE( ptPhotReco, etaPhotReco, phiPhotReco, ePhotReco );
    photon.hcalIso = hcalIsoPhotReco;
    photon.ecalIso = ecalIsoPhotReco;
    photon.nTrkIso = nTrkIsoPhotReco;
    photon.ptTrkIso = ptTrkIsoPhotReco;
    photon.clusterMaj = clusterMajPhotReco;
    photon.clusterMin = clusterMinPhotReco;
    photon.hasPixelSeed = hasPixelSeedPhotReco;
   
    photon.pid_hlwTrack = pid_hlwTrackPhotReco;
    photon.pid_twrHCAL = pid_towerHCALPhotReco;
    photon.pid_HoverE = pid_HoverEPhotReco;
    photon.pid_jurECAL = pid_jurECALPhotReco;
    photon.pid_etawid = pid_sIEtaIEtaPhotReco;




    if( ptPhotReco<25. ) continue;
    //if( fabs(etaPhotReco)>1.3 ) continue;
    if( clusterMinPhotReco<0.15 ) continue; //protection vs EB spikes

    if( clusterMinPhotReco>0.3 ) continue;
    if( clusterMajPhotReco>0.35 ) continue;



    // match photon to MC:
    TLorentzVector photonGen;
    photonGen.SetPtEtaPhiE( ptPhotGen, etaPhotGen, phiPhotGen, ePhotGen );

    bool isMatched = photon.DeltaR(photonGen)<0.5;

    if( !isMatched ) continue;
    

    //------------
    // PHOTON ID:
    //------------

    std::vector<bool> idpasseg(5);

    bool passedID = passedPhotonID(photon, hggtightid, &idpasseg);
    passedID = ( passedID && !hasPixelSeedPhotReco );

    //if( !passedID ) continue; 

    pid_hcalIso = photon.pid_HoverE + photon.pid_twrHCAL / photon.Pt();
    pid_ecalIso = photon.pid_jurECAL*cosh(photon.Eta()) / photon.Energy();
    pid_hlwTrackIso = photon.pid_hlwTrack / photon.Pt();
  
    hp_hcalIsoMean_vs_rhoCalo->Fill( rhoCalo, hcalIsoPhotReco, eventWeight );
    hp_hcalIsoPIDMean_vs_rhoCalo->Fill( rhoCalo, pid_hcalIso, eventWeight );
    hp_hcalIsoMean_vs_rhoPF->Fill( rhoPF, hcalIsoPhotReco, eventWeight );
    hp_hcalIsoPIDMean_vs_rhoPF->Fill( rhoPF, pid_hcalIso, eventWeight );

    hp_ecalIsoMean_vs_rhoCalo->Fill( rhoCalo, ecalIsoPhotReco, eventWeight );
    hp_ecalIsoPIDMean_vs_rhoCalo->Fill( rhoCalo, pid_ecalIso, eventWeight );
    hp_ecalIsoMean_vs_rhoPF->Fill( rhoPF, ecalIsoPhotReco, eventWeight );
    hp_ecalIsoPIDMean_vs_rhoPF->Fill( rhoPF, pid_ecalIso, eventWeight );

    hp_ptTrkIsoMean_vs_rhoCalo->Fill( rhoCalo, ptTrkIsoPhotReco, eventWeight );
    hp_hlwTrackIsoMean_vs_rhoCalo->Fill( rhoCalo, pid_hlwTrackIso, eventWeight );
    hp_ptTrkIsoMean_vs_rhoPF->Fill( rhoPF, ptTrkIsoPhotReco, eventWeight );
    hp_hlwTrackIsoMean_vs_rhoPF->Fill( rhoPF, pid_hlwTrackIso, eventWeight );

    bin_rhoPF   = hp_hcalIsoMean_vs_rhoPF->FindFixBin(rhoPF);
    bin_rhoCalo = hp_hcalIsoMean_vs_rhoCalo->FindFixBin(rhoCalo);

    // arrays start from zero:
    bin_rhoPF--;
    bin_rhoCalo--;

    if( bin_rhoPF >= 0  && bin_rhoPF < nBinsRho ) {
      tree_hcalIsoPhotReco_rhoPF->Fill();
      tree_hcalIsoPID_rhoPF->Fill();
      tree_ecalIsoPhotReco_rhoPF->Fill();
      tree_ecalIsoPID_rhoPF->Fill();
      tree_ptTrkIsoPhotReco_rhoPF->Fill();
      tree_hlwTrackIso_rhoPF->Fill();

      vh1_hcalIsoPhotReco_rhoPF[bin_rhoPF]->Fill( hcalIsoPhotReco, eventWeight );
      vh1_hcalIsoPID_rhoPF[bin_rhoPF]->Fill( pid_hcalIso, eventWeight );
      vh1_ecalIsoPhotReco_rhoPF[bin_rhoPF]->Fill( ecalIsoPhotReco, eventWeight );
      vh1_ecalIsoPID_rhoPF[bin_rhoPF]->Fill( pid_ecalIso, eventWeight );
      vh1_ptTrkIsoPhotReco_rhoPF[bin_rhoPF]->Fill( ptTrkIsoPhotReco, eventWeight );
      vh1_hlwTrackIso_rhoPF[bin_rhoPF]->Fill( pid_hlwTrackIso, eventWeight );
    }

    if( bin_rhoCalo >= 0  && bin_rhoCalo < nBinsRho ) {
      tree_hcalIsoPhotReco_rhoCalo->Fill();
      tree_hcalIsoPID_rhoCalo->Fill();
      tree_ecalIsoPhotReco_rhoCalo->Fill();
      tree_ecalIsoPID_rhoCalo->Fill();
      tree_ptTrkIsoPhotReco_rhoCalo->Fill();
      tree_hlwTrackIso_rhoCalo->Fill();

      vh1_hcalIsoPhotReco_rhoCalo[bin_rhoCalo]->Fill( hcalIsoPhotReco, eventWeight );
      vh1_hcalIsoPID_rhoCalo[bin_rhoCalo]->Fill( pid_hcalIso, eventWeight );
      vh1_ecalIsoPhotReco_rhoCalo[bin_rhoCalo]->Fill( ecalIsoPhotReco, eventWeight );
      vh1_ecalIsoPID_rhoCalo[bin_rhoCalo]->Fill( pid_ecalIso, eventWeight );
      vh1_ptTrkIsoPhotReco_rhoCalo[bin_rhoCalo]->Fill( ptTrkIsoPhotReco, eventWeight );
      vh1_hlwTrackIso_rhoCalo[bin_rhoCalo]->Fill( pid_hlwTrackIso, eventWeight );
    }


    // fill histograms used for working point evaluation
    // only for the no PU case
    if( nvertex==1 ) {

      denom_hcalIso += eventWeight;
      denom_pid_hcalIso += eventWeight;
      denom_ecalIso += eventWeight;
      denom_pid_ecalIso += eventWeight;
      denom_ptTrkIso += eventWeight;
      denom_pid_hlwTrackIso += eventWeight;

    
      bool passed_hcalIsoPhotReco = ( hcalIsoPhotReco < hcalIsoPhotReco_threshRel ) || ( hcalIsoPhotReco*photon.Energy() < hcalIsoPhotReco_threshAbs );
      bool passed_pid_hcalIso = ( pid_hcalIso < hcalIsoPhotReco_threshRel ) || ( pid_hcalIso*photon.Energy() < hcalIsoPhotReco_threshAbs );
      bool passed_ecalIsoPhotReco = ( ecalIsoPhotReco < ecalIsoPhotReco_threshRel ) || ( ecalIsoPhotReco*photon.Energy() < ecalIsoPhotReco_threshAbs );
      bool passed_pid_ecalIso = ( pid_ecalIso < ecalIsoPhotReco_threshRel ) || ( pid_ecalIso*photon.Energy() < ecalIsoPhotReco_threshAbs );
      bool passed_ptTrkIsoPhotReco = ( ptTrkIsoPhotReco < ptTrkIsoPhotReco_threshRel );
      bool passed_pid_hlwTrackIso = ( pid_hlwTrackIso < ptTrkIsoPhotReco_threshRel );

    //bool passed_hcalIsoPhotReco = ( hcalIsoPhotReco < hcalIsoPhotReco_threshRel ) || ( hcalIsoPhotReco*photon.Energy() < hcalIsoPhotReco_threshAbs );
    //bool passed_pid_hcalIso = ( pid_hcalIso < (pid_hcalIso_threshRel*photon.Pt()+pid_hcalIso_threshAbs) );
    //bool passed_ecalIsoPhotReco = ( ecalIsoPhotReco < ecalIsoPhotReco_threshRel ) || ( ecalIsoPhotReco*photon.Energy() < ecalIsoPhotReco_threshAbs );
    //bool passed_pid_ecalIso = ( pid_ecalIso < (pid_ecalIso_threshRel*photon.Pt()+pid_ecalIso_threshAbs) );
    //bool passed_ptTrkIsoPhotReco = ( ptTrkIsoPhotReco < ptTrkIsoPhotReco_threshRel );
    //bool passed_pid_hlwTrackIso = ( pid_hlwTrackIso < (pid_hlwTrackIso_threshRel*photon.Pt()+pid_hlwTrackIso_threshAbs) );

      if( passed_hcalIsoPhotReco )       num_hcalIso += eventWeight;
      if( passed_pid_hcalIso )           num_pid_hcalIso += eventWeight;
      if( passed_ecalIsoPhotReco )       num_ecalIso += eventWeight;
      if( passed_pid_ecalIso )           num_pid_ecalIso += eventWeight;
      if( passed_ptTrkIsoPhotReco )      num_ptTrkIso += eventWeight;
      if( passed_pid_hlwTrackIso )       num_pid_hlwTrackIso += eventWeight;

    }


  } //for entries


  // compute working point efficiencies:

  float eff_hcalIso     = num_hcalIso / denom_hcalIso;
  float eff_hcalIsoPID  = num_pid_hcalIso / denom_pid_hcalIso;
  float eff_ecalIso     = num_ecalIso / denom_ecalIso;
  float eff_ecalIsoPID  = num_pid_ecalIso / denom_pid_ecalIso;
  float eff_ptTrkIso    = num_ptTrkIso / denom_ptTrkIso;
  float eff_hlwTrackIso = num_pid_hlwTrackIso / denom_pid_hlwTrackIso;

  TH1D* h1_hcalIsoPhotRecoThresh_vs_rhoPF = new TH1D("hcalIsoPhotRecoThresh_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  TH1D* h1_hcalIsoPIDThresh_vs_rhoPF = new TH1D("hcalIsoPIDThresh_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  TH1D* h1_ecalIsoPhotRecoThresh_vs_rhoPF = new TH1D("ecalIsoPhotRecoThresh_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  TH1D* h1_ecalIsoPIDThresh_vs_rhoPF = new TH1D("ecalIsoPIDThresh_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  TH1D* h1_ptTrkIsoPhotRecoThresh_vs_rhoPF = new TH1D("ptTrkIsoPhotRecoThresh_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  TH1D* h1_hlwTrackIsoThresh_vs_rhoPF = new TH1D("hlwTrackIsoThresh_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);

  TH1D* h1_hcalIsoPhotRecoThresh_vs_rhoCalo = new TH1D("hcalIsoPhotRecoThresh_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  TH1D* h1_hcalIsoPIDThresh_vs_rhoCalo = new TH1D("hcalIsoPIDThresh_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  TH1D* h1_ecalIsoPhotRecoThresh_vs_rhoCalo = new TH1D("ecalIsoPhotRecoThresh_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  TH1D* h1_ecalIsoPIDThresh_vs_rhoCalo = new TH1D("ecalIsoPIDThresh_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  TH1D* h1_ptTrkIsoPhotRecoThresh_vs_rhoCalo = new TH1D("ptTrkIsoPhotRecoThresh_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  TH1D* h1_hlwTrackIsoThresh_vs_rhoCalo = new TH1D("hlwTrackIsoThresh_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);


  // find const efficiency thresholds:

  for( unsigned iBinRho=0; iBinRho<nBinsRho; ++iBinRho ) {

    float thresh_hcalIsoPF = getThresh_constEff( tree_hcalIsoPhotReco_rhoPF, eff_hcalIso, iBinRho, hcalIsoPhotReco_threshAbs );
    h1_hcalIsoPhotRecoThresh_vs_rhoPF->SetBinContent( iBinRho+1, thresh_hcalIsoPF );
    float thresh_hcalIsoPIDPF = getThresh_constEff( tree_hcalIsoPID_rhoPF, eff_hcalIsoPID, iBinRho, hcalIsoPhotReco_threshAbs );
    h1_hcalIsoPIDThresh_vs_rhoPF->SetBinContent( iBinRho+1, thresh_hcalIsoPIDPF );

    float thresh_ecalIsoPF = getThresh_constEff( tree_ecalIsoPhotReco_rhoPF, eff_ecalIso, iBinRho, ecalIsoPhotReco_threshAbs );
    h1_ecalIsoPhotRecoThresh_vs_rhoPF->SetBinContent( iBinRho+1, thresh_ecalIsoPF );
    float thresh_ecalIsoPIDPF = getThresh_constEff( tree_ecalIsoPID_rhoPF, eff_ecalIsoPID, iBinRho, ecalIsoPhotReco_threshAbs );
    h1_ecalIsoPIDThresh_vs_rhoPF->SetBinContent( iBinRho+1, thresh_ecalIsoPIDPF );

    float thresh_ptTrkIsoPF = getThresh_constEff( tree_ptTrkIsoPhotReco_rhoPF, eff_ptTrkIso, iBinRho, 0. );
    h1_ptTrkIsoPhotRecoThresh_vs_rhoPF->SetBinContent( iBinRho+1, thresh_ptTrkIsoPF );
    float thresh_hlwTrackIsoPF = getThresh_constEff( tree_hlwTrackIso_rhoPF, eff_hlwTrackIso, iBinRho, 0. );
    h1_hlwTrackIsoThresh_vs_rhoPF->SetBinContent( iBinRho+1, thresh_hlwTrackIsoPF );


    float thresh_hcalIsoCalo = getThresh_constEff( tree_hcalIsoPhotReco_rhoCalo, eff_hcalIso, iBinRho, hcalIsoPhotReco_threshAbs );
    h1_hcalIsoPhotRecoThresh_vs_rhoCalo->SetBinContent( iBinRho+1, thresh_hcalIsoCalo );
    float thresh_hcalIsoPIDCalo = getThresh_constEff( tree_hcalIsoPID_rhoCalo, eff_hcalIsoPID, iBinRho, hcalIsoPhotReco_threshAbs );
    h1_hcalIsoPIDThresh_vs_rhoCalo->SetBinContent( iBinRho+1, thresh_hcalIsoPIDCalo );

    float thresh_ecalIsoCalo = getThresh_constEff( tree_ecalIsoPhotReco_rhoCalo, eff_ecalIso, iBinRho, ecalIsoPhotReco_threshAbs );
    h1_ecalIsoPhotRecoThresh_vs_rhoCalo->SetBinContent( iBinRho+1, thresh_ecalIsoCalo );
    float thresh_ecalIsoPIDCalo = getThresh_constEff( tree_ecalIsoPID_rhoCalo, eff_ecalIsoPID, iBinRho, ecalIsoPhotReco_threshAbs );
    h1_ecalIsoPIDThresh_vs_rhoCalo->SetBinContent( iBinRho+1, thresh_ecalIsoPIDCalo );

    float thresh_ptTrkIsoCalo = getThresh_constEff( tree_ptTrkIsoPhotReco_rhoCalo, eff_ptTrkIso, iBinRho, 0. );
    h1_ptTrkIsoPhotRecoThresh_vs_rhoCalo->SetBinContent( iBinRho+1, thresh_ptTrkIsoCalo );
    float thresh_hlwTrackIsoCalo = getThresh_constEff( tree_hlwTrackIso_rhoCalo, eff_hlwTrackIso, iBinRho, 0. );
    h1_hlwTrackIsoThresh_vs_rhoCalo->SetBinContent( iBinRho+1, thresh_hlwTrackIsoCalo );


  }


  // FIT TRENDS VS. RHO:

  TF1* linePF = new TF1("linePF", "[0] + [1]*x", 0., rhoPF_max);
  std::cout << std::endl << std::endl << "### HCAL vs rhoPF: " << std::endl;
  hp_hcalIsoMean_vs_rhoPF->Fit(linePF, "R");
  std::cout << std::endl << std::endl << "### HCAL PID vs rhoPF: " << std::endl;
  hp_hcalIsoPIDMean_vs_rhoPF->Fit(linePF, "R");
  //std::cout << std::endl << std::endl << "### ECAL vs rhoPF: " << std::endl;
  //hp_ecalIsoMean_vs_rhoPF->Fit(linePF, "R");
  std::cout << std::endl << std::endl << "### ECAL PID vs rhoPF: " << std::endl;
  hp_ecalIsoPIDMean_vs_rhoPF->Fit(linePF, "R");

  TF1* lineCalo = new TF1("lineCalo", "[0] + [1]*x", 0., rhoCalo_max);
  std::cout << std::endl << std::endl << "### HCAL vs rhoCalo: " << std::endl;
  hp_hcalIsoMean_vs_rhoCalo->Fit(lineCalo, "R");
  std::cout << std::endl << std::endl << "### HCAL PID vs rhoCalo: " << std::endl;
  hp_hcalIsoPIDMean_vs_rhoCalo->Fit(lineCalo, "R");
  //std::cout << std::endl << std::endl << "### ECAL vs rhoCalo: " << std::endl;
  //hp_ecalIsoMean_vs_rhoCalo->Fit(lineCalo, "R");
  std::cout << std::endl << std::endl << "### ECAL PID vs rhoCalo: " << std::endl;
  hp_ecalIsoPIDMean_vs_rhoCalo->Fit(lineCalo, "R");



  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaRhoStudies_"+dataset;
  else {
   if(dataset!="") outfileName = "RhoStudies_"+dataset;
   else outfileName = "PhotonJet";
  }

  outfileName = outfileName;
  if( photonID!="medium" ) outfileName = outfileName + "_" + photonID;
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();

  hp_hcalIsoMean_vs_rhoCalo->Write();
  hp_hcalIsoPIDMean_vs_rhoCalo->Write();
  hp_hcalIsoMean_vs_rhoPF->Write();
  hp_hcalIsoPIDMean_vs_rhoPF->Write();

  hp_ecalIsoMean_vs_rhoCalo->Write();
  hp_ecalIsoPIDMean_vs_rhoCalo->Write();
  hp_ecalIsoMean_vs_rhoPF->Write();
  hp_ecalIsoPIDMean_vs_rhoPF->Write();

  hp_ptTrkIsoMean_vs_rhoCalo->Write();
  hp_hlwTrackIsoMean_vs_rhoCalo->Write();
  hp_ptTrkIsoMean_vs_rhoPF->Write();
  hp_hlwTrackIsoMean_vs_rhoPF->Write();

  h1_hcalIsoPhotRecoThresh_vs_rhoPF->Write();
  h1_hcalIsoPIDThresh_vs_rhoPF->Write();
  h1_ecalIsoPhotRecoThresh_vs_rhoPF->Write();
  h1_ecalIsoPIDThresh_vs_rhoPF->Write();
  h1_ptTrkIsoPhotRecoThresh_vs_rhoPF->Write();
  h1_hlwTrackIsoThresh_vs_rhoPF->Write();

  h1_hcalIsoPhotRecoThresh_vs_rhoCalo->Write();
  h1_hcalIsoPIDThresh_vs_rhoCalo->Write();
  h1_ecalIsoPhotRecoThresh_vs_rhoCalo->Write();
  h1_ecalIsoPIDThresh_vs_rhoCalo->Write();
  h1_ptTrkIsoPhotRecoThresh_vs_rhoCalo->Write();
  h1_hlwTrackIsoThresh_vs_rhoCalo->Write();

  outFile->Close();

  delete tree;
  tree = 0;


}


void addInput( const std::string& dataset ) {


  // opening from filelist now deprecated (files have to be merged with merge_and_setWeights
  //std::string infileName = "files_PhotonJet_2ndLevel_" + dataset+"_" + ALGOTYPE_ +".txt";


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "PhotonJet_2ndLevelTreeW_" + dataset + "_pfakt5";
    infileName += ".root";
    TFile* infile = TFile::Open(infileName.c_str(), "read");
    if( infile==0 ) {
      std::cout << "Didn't find file '" << infileName << "'. Did you forget to finalize (i.e. the \"W\")?" << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(77);
    }
    std::string treeName = infileName +"/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infile->Close();


} //addinput


  

bool passedPhotonID(const AnalysisPhoton& photon, photonidegcuts const& pid, vector<bool> *vpass) {

 // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
 bool ptiso =   (photon.pid_hlwTrack < photon.Pt() * pid.trackiso_rel + pid.trackiso_abs);
 bool ecaliso = (photon.pid_jurECAL  < photon.Pt() * pid.ecaliso_rel + pid.ecaliso_abs);
 bool hcaliso = (photon.pid_twrHCAL  < photon.Pt() * pid.hcaliso_rel + pid.hcaliso_abs);
 bool hoveiso = (photon.pid_HoverE < pid.hovereiso);
 bool setaeta =  photon.pid_etawid < pid.setaetaEB;
 if(TMath::Abs(photon.Eta()) > 1.44) {
   setaeta = photon.pid_etawid < pid.setaetaEE;
 }

 if (vpass) {
   //assert((*vpass).size()==7);
   if((*vpass).size()!=5) { cout << "major failure! (*vpass).size()!=7.. die!" << endl; exit(0) ; }
   (*vpass)[0] = ptiso;
   (*vpass)[1] = ecaliso;
   (*vpass)[2] = hcaliso;
   (*vpass)[3] = hoveiso;
   (*vpass)[4] = setaeta;
 }

 return (ptiso && hcaliso && ecaliso && hoveiso && setaeta);
}



std::vector<TH1D*> allocateHistoVector( int vectorSize, const std::string& name, int nBins, float xMin, float xMax ) {

  std::vector<TH1D*> returnVector;

  for( unsigned iHisto=0; iHisto<vectorSize; ++iHisto ) {

    char iName[500];
    sprintf( iName, "%s_%d", name.c_str(), iHisto );
    
    TH1D* h1_new = new TH1D( iName, "", nBins, xMin, xMax );
    h1_new->Sumw2();

    returnVector.push_back( h1_new );

  }

  return returnVector;

}



float getThresh_constEff( TTree* ttree, float eff, int iBinRho, float threshAbs) {

  std::cout << "Scanning: " << tree->GetName() << std::endl;

  float foundThresh = 0.;

  char denomCondition[500];
  sprintf( denomCondition, "bin_rho == %d", iBinRho );
  float nEntries = ttree->GetEntries(denomCondition);


  for( unsigned i=0; i<100; ++i) {

    float threshRel = (float)i/100.;

    char numCondition[800];
    sprintf( numCondition, "bin_rho == %d  && ( isoVar < %f  || isoVar*ePhotReco < %f )", iBinRho, threshRel, threshAbs );
    float thisEff = ttree->GetEntries(numCondition)/nEntries;

    if( thisEff > eff ) {
      foundThresh = threshRel;
      break;
    }

  } //for

  return foundThresh;

}



