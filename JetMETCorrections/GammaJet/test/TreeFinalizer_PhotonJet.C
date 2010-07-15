#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "fitTools.C"




bool DEBUG_ = false;
bool BINNINGFINO_ = false;
bool ONEVTX_ = false;
bool NO2ndJETABS = true;
bool TIGHTDELTAPHI_ = false;
std::string RECOTYPE_;
std::string ALGOTYPE_;
std::string suffix;

TChain* tree;
Double_t totalLumi=0.;


void addInput(const std::string& dataset);
std::vector<TH1F*> getResponseHistos(const std::string& name);



void finalize(const std::string& dataset, std::string recoType, std::string jetAlgo="akt5", float secondJetThreshold=0.5, bool useGenJets=false, bool MCassoc=false) {

  bool noJetSelection = ( secondJetThreshold < 0. );

  tree = new TChain("jetTree");

  RECOTYPE_ = recoType;

  ALGOTYPE_ = (recoType=="calo") ? jetAlgo : recoType+jetAlgo;
  if( recoType=="jpt"&&jetAlgo=="akt5" ) ALGOTYPE_="jptak5"; 

  std::string infileName, treeName;

  suffix = "_"+ALGOTYPE_;
  if( useGenJets ) suffix = suffix + "_GENJETS";



  if( dataset=="PhotonJet_Summer09") {

    addInput( "PhotonJet_Summer09_Pt15" );
    addInput( "PhotonJet_Summer09_Pt30" );
    addInput( "PhotonJet_Summer09_Pt80" );
    addInput( "PhotonJet_Summer09_Pt170" );
    addInput( "PhotonJet_Summer09_Pt300" );
    addInput( "PhotonJet_Summer09_Pt470" );
    addInput( "PhotonJet_Summer09_Pt800" );
    addInput( "PhotonJet_Summer09_Pt1400" );

  } else if( dataset=="QCD_Summer09" ) {

    addInput( "QCD_Summer09_Pt15" );
    addInput( "QCD_Summer09_Pt30" );
    addInput( "QCD_Summer09_Pt80" );
    addInput( "QCD_Summer09_Pt170" );
    addInput( "QCD_Summer09_Pt300" );
    addInput( "QCD_Summer09_Pt470" );
    addInput( "QCD_Summer09_Pt800" );
    addInput( "QCD_Summer09_Pt1400" );

  } else if( dataset=="PhotonJetQCD_Summer09" ) {

    addInput( "PhotonJet_Summer09_Pt15" );
    addInput( "PhotonJet_Summer09_Pt30" );
    addInput( "PhotonJet_Summer09_Pt80" );
    addInput( "PhotonJet_Summer09_Pt170" );
    addInput( "PhotonJet_Summer09_Pt300" );
    addInput( "PhotonJet_Summer09_Pt470" );
    addInput( "PhotonJet_Summer09_Pt800" );
    addInput( "PhotonJet_Summer09_Pt1400" );
    addInput( "QCD_Summer09_Pt15" );
    addInput( "QCD_Summer09_Pt30" );
    addInput( "QCD_Summer09_Pt80" );
    addInput( "QCD_Summer09_Pt170" );
    addInput( "QCD_Summer09_Pt300" );
    addInput( "QCD_Summer09_Pt470" );
    addInput( "QCD_Summer09_Pt800" );
    addInput( "QCD_Summer09_Pt1400" );

  } else if( dataset=="PhotonJet_Spring10" ) {

    addInput( "PhotonJet_Spring10_Pt5to15" );
    addInput( "PhotonJet_Spring10_Pt15" );
    addInput( "PhotonJet_Spring10_Pt30" );
    addInput( "PhotonJet_Spring10_Pt80" );
    addInput( "PhotonJet_Spring10_Pt170" );

  } else if( dataset=="PhotonJet_Summer1036X" ) {

    addInput( "PhotonJet_Summer1036X_Pt5to15" );
    addInput( "PhotonJet_Summer1036X_Pt15to20" );
    addInput( "PhotonJet_Summer1036X_Pt20to30" );
    addInput( "PhotonJet_Summer1036X_Pt30to50" );
    addInput( "PhotonJet_Summer1036X_Pt50to80" );
    addInput( "PhotonJet_Summer1036X_Pt80to120" );
    addInput( "PhotonJet_Summer1036X_Pt120to170" );
    addInput( "PhotonJet_Summer1036X_Pt170to300" );

  } else if( dataset=="QCD_Spring10" ) {

    //addInput( "QCD_Spring10_Pt5to15" );
    addInput( "QCD_Spring10_Pt15to20" );
    addInput( "QCD_Spring10_Pt20to30" );
    addInput( "QCD_Spring10_Pt30to50" );
    addInput( "QCD_Spring10_Pt50to80" );
    addInput( "QCD_Spring10_Pt80" );
    addInput( "QCD_Spring10_Pt170" );
    //addInput( "QCD_Spring10_Pt300" );

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

  } else if( dataset=="DATA_EG_NEW" ) {

    addInput( "MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1_80MeV_MET" );
    addInput( "EG_Run2010A_Jun14thReReco_v1_80MeV_MET" );
    addInput( "EG_Run2010A_PromptReco_v4_80MeV_MET" );

  } else if( dataset=="DATA_EG_1" ) {

    addInput( "MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1" );
    addInput( "EG_Run2010A-Jun14thReReco_v1" );

  } else if( dataset=="DATA_EG_2" ) {

    addInput( "EG_Run2010A-PromptReco-v4" );
    addInput( "EG_Run2010A_PromptReco_v4_139347_139375" );
    addInput( "EG_Run2010A_PromptReco_v4_139376_139459" );

  } else {
  
    addInput( dataset );

  }



  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Float_t ptMax = 100.;
  int nBins_photBinning = ptPhot_binning.size()-1;
  Double_t ptPhot_binning_array[100];
  for( unsigned i=0; i<nBins_photBinning; ++i)
    ptPhot_binning_array[i] = ptPhot_binning[i];

  TH1D* h1_ptPhot = new TH1D("ptPhot", "", nBins_photBinning-1, ptPhot_binning_array);
  //TH1D* h1_ptPhot = new TH1D("ptPhot", "", 10, ptPhot_binning[0], ptMax);
  h1_ptPhot->Sumw2();
  TH1D* h1_ptJetReco = new TH1D("ptJetReco", "", 10, 0., ptMax);
  h1_ptJetReco->Sumw2();
  TH1D* h1_pt2ndJetReco = new TH1D("pt2ndJetReco", "", 10, 5., ptMax);
  h1_pt2ndJetReco->Sumw2();
  TH1D* h1_etaPhot = new TH1D("etaPhot", "", 15, -1.3, 1.3);
  h1_etaPhot->Sumw2();
  TH1D* h1_phiPhot = new TH1D("phiPhot", "", 15, -3.1416, 3.1416);
  h1_phiPhot->Sumw2();
  TH1D* h1_deltaPhi_Nm1 = new TH1D("deltaPhi_Nm1", "", 15, 2.1, 3.1416);
  h1_deltaPhi_Nm1->Sumw2();
  TH1D* h1_ptSecondJetRel_Nm1 = new TH1D("ptSecondJetRel_Nm1", "", 15, 0., 1.5);
  h1_ptSecondJetRel_Nm1->Sumw2();

  TH1D* h1_deltaPhi_2ndJet_medium = new TH1D("deltaPhi_2ndJet_medium", "", 15, 0., 3.1416 );
  h1_deltaPhi_2ndJet_medium->Sumw2();
  TH1D* h1_deltaPhi_2ndJet_loose = new TH1D("deltaPhi_2ndJet_loose", "", 15, 0., 3.1416 );
  h1_deltaPhi_2ndJet_loose->Sumw2();

  TH1D* h1_ptPhot_loose = new TH1D("ptPhot_loose", "", nBins_photBinning-1, ptPhot_binning_array);
  //TH1D* h1_ptPhot_loose = new TH1D("ptPhot_loose", "", 10, ptPhot_binning[0], ptMax);
  h1_ptPhot_loose->Sumw2();
  TH1D* h1_etaPhot_loose = new TH1D("etaPhot_loose", "", 15, -1.3, 1.3);
  h1_etaPhot_loose->Sumw2();
  TH1D* h1_phiPhot_loose = new TH1D("phiPhot_loose", "", 15, -3.1416, 3.1416);
  h1_phiPhot_loose->Sumw2();

  TH1D* h1_met_loose = new TH1D("met_loose", "", 100., 0., 500.);
  h1_met_loose->Sumw2();
  TH1D* h1_met_medium = new TH1D("met_medium", "", 100., 0., 500.);
  h1_met_medium->Sumw2();

  TH1D* h1_deltaPhi_phot_met_loose = new TH1D("deltaPhi_phot_met_loose", "", 15, -3.1416, 3.1416);
  h1_deltaPhi_phot_met_loose->Sumw2();
  TH1D* h1_deltaPhi_phot_met_medium = new TH1D("deltaPhi_phot_met_medium", "", 15, -3.1416, 3.1416);
  h1_deltaPhi_phot_met_medium->Sumw2();

  //TH1D* h1_ptPhot_medium = new TH1D("ptPhot_medium", "", 10, ptPhot_binning[0], ptMax);
  TH1D* h1_ptPhot_medium = new TH1D("ptPhot_medium", "", nBins_photBinning-1, ptPhot_binning_array);
  h1_ptPhot_medium->Sumw2();

  //TH1D* h1_ptPhot_clusterOK = new TH1D("ptPhot_clusterOK", "", 10, ptPhot_binning[0], ptMax);
  TH1D* h1_ptPhot_clusterOK = new TH1D("ptPhot_clusterOK", "", nBins_photBinning-1, ptPhot_binning_array);
  h1_ptPhot_clusterOK->Sumw2();
  TH1D* h1_etaPhot_clusterOK = new TH1D("etaPhot_clusterOK", "", 15, -1.3, 1.3);
  h1_etaPhot_clusterOK->Sumw2();
  TH1D* h1_phiPhot_clusterOK = new TH1D("phiPhot_clusterOK", "", 15, -3.1416, 3.1416);
  h1_phiPhot_clusterOK->Sumw2();
  TH1D* h1_deltaPhi_clusterOK = new TH1D("deltaPhi_clusterOK", "", 15, 2.1, 3.1416);
  h1_deltaPhi_clusterOK->Sumw2();

  //TH1D* h1_ptPhot_clusterOK_isolated = new TH1D("ptPhot_clusterOK_isolated", "", 28, ptPhot_binning[0], ptMax);
  TH1D* h1_ptPhot_clusterOK_isolated = new TH1D("ptPhot_clusterOK_isolated", "", nBins_photBinning-1, ptPhot_binning_array);
  h1_ptPhot_clusterOK_isolated->Sumw2();
  TH1D* h1_etaPhot_clusterOK_isolated = new TH1D("etaPhot_clusterOK_isolated", "", 15, -1.3, 1.3);
  h1_etaPhot_clusterOK_isolated->Sumw2();
  TH1D* h1_phiPhot_clusterOK_isolated = new TH1D("phiPhot_clusterOK_isolated", "", 15, -3.1416, 3.1416);
  h1_phiPhot_clusterOK_isolated->Sumw2();

  TH1D* h1_hcalIsoPhotReco_Nm1 = new TH1D("hcalIsoPhotReco_Nm1", "", 10, 0., 0.5);
  h1_hcalIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_hcalIsoEnergyPhotReco_Nm1 = new TH1D("hcalIsoEnergyPhotReco_Nm1", "", 10, 0., 15.);
  h1_hcalIsoEnergyPhotReco_Nm1->Sumw2();
  TH1D* h1_ecalIsoPhotReco_Nm1 = new TH1D("ecalIsoPhotReco_Nm1", "", 20, 0., 1.);
  h1_ecalIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_ecalIsoEnergyPhotReco_Nm1 = new TH1D("ecalIsoEnergyPhotReco_Nm1", "", 20, 0., 30.);
  h1_ecalIsoEnergyPhotReco_Nm1->Sumw2();
  TH1D* h1_ptTrkIsoPhotReco_Nm1 = new TH1D("ptTrkIsoPhotReco_Nm1", "", 20, 0., 1.);
  h1_ptTrkIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_nTrkIsoPhotReco_Nm1 = new TH1D("nTrkIsoPhotReco_Nm1", "", 11, -0.5, 10.5);
  h1_nTrkIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_clusterMajPhotReco_Nm1 = new TH1D("clusterMajPhotReco_Nm1", "", 50, 0., 1.5);
  h1_clusterMajPhotReco_Nm1->Sumw2();
  TH1D* h1_clusterMinPhotReco_Nm1 = new TH1D("clusterMinPhotReco_Nm1", "", 50, 0., 0.6);
  h1_clusterMinPhotReco_Nm1->Sumw2();

  TH1D* h1_hcalIsoPhotReco_clusterOK = new TH1D("hcalIsoPhotReco_clusterOK", "", 10, 0., 0.5);
  h1_hcalIsoPhotReco_clusterOK->Sumw2();
  TH1D* h1_hcalIsoEnergyPhotReco_clusterOK = new TH1D("hcalIsoEnergyPhotReco_clusterOK", "", 10, 0., 15.);
  h1_hcalIsoEnergyPhotReco_clusterOK->Sumw2();
  TH1D* h1_ecalIsoPhotReco_clusterOK = new TH1D("ecalIsoPhotReco_clusterOK", "", 20, 0., 1.);
  h1_ecalIsoPhotReco_clusterOK->Sumw2();
  TH1D* h1_ecalIsoEnergyPhotReco_clusterOK = new TH1D("ecalIsoEnergyPhotReco_clusterOK", "", 20, 0., 30.);
  h1_ecalIsoEnergyPhotReco_clusterOK->Sumw2();
  TH1D* h1_ptTrkIsoPhotReco_clusterOK = new TH1D("ptTrkIsoPhotReco_clusterOK", "", 20, 0., 1.);
  h1_ptTrkIsoPhotReco_clusterOK->Sumw2();
  TH1D* h1_nTrkIsoPhotReco_clusterOK = new TH1D("nTrkIsoPhotReco_clusterOK", "", 11, -0.5, 10.5);
  h1_nTrkIsoPhotReco_clusterOK->Sumw2();

  TH1D* h1_clusterMajPhotReco_isolated = new TH1D("clusterMajPhotReco_isolated", "", 30, 0., 1.5);
  h1_clusterMajPhotReco_isolated->Sumw2();
  TH1D* h1_clusterMinPhotReco_isolated = new TH1D("clusterMinPhotReco_isolated", "", 20, 0., 0.6);
  h1_clusterMinPhotReco_isolated->Sumw2();

  TH1D* h1_clusterMajPhotReco= new TH1D("clusterMajPhotReco", "", 30, 0., 1.5);
  h1_clusterMajPhotReco->Sumw2();
  TH1D* h1_clusterMinPhotReco= new TH1D("clusterMinPhotReco", "", 20, 0., 0.6);
  h1_clusterMinPhotReco->Sumw2();

  TH1D* h1_deltaPhi_clusterOK_isolated = new TH1D("deltaPhi_clusterOK_isolated", "", 15, 2.1, 3.1416);
  h1_deltaPhi_clusterOK_isolated->Sumw2();
  TH1D* h1_ptSecondJetRel_clusterOK_isolated = new TH1D("ptSecondJetRel_clusterOK_isolated", "", 15, 0., 1.5);
  h1_ptSecondJetRel_clusterOK_isolated->Sumw2();
//TH1D* h1_response = new TH1D("response", "", 15, 0., 2.);
//h1_response->Sumw2();

//std::vector<TH1F*> h1_ptPhotMean;
//for( unsigned i=0; i<(ptPhot_binning.size()-1); ++i ) {
//  Float_t ptMin = ptPhot_binning[i];
//  Float_t ptMax = ptPhot_binning[i+1];
//  char histoName[100];
//  sprintf( histoName, "ptPhot_%.0f_%.0f", ptMin, ptMax );
//  TH1F* tmp = new Th1F(histoName, "", 


  Double_t ptPhotBinning_array[200]; //ugly! no more than 200 pt bins supported
  for( unsigned i=0; i<ptPhot_binning.size(); ++i )
    ptPhotBinning_array[i] = ptPhot_binning[i];


  TProfile* hp_ptJetGenMean = new TProfile("ptJetGenMean", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  TProfile* hp_ptPhotMean = new TProfile("ptPhotMean", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  TProfile* hp_ptPhotMean_loose = new TProfile("ptPhotMean_loose", "", ptPhot_binning.size()-1, ptPhotBinning_array);

  std::vector<TH1F*> h1_response              = getResponseHistos("response");
  std::vector<TH1F*> h1_responseGEN           = getResponseHistos("responseGEN");
  std::vector<TH1F*> h1_responseMPF           = getResponseHistos("responseMPF");

  std::vector<TH1F*> h1_response_clusterOK    = getResponseHistos("response_clusterOK");
  std::vector<TH1F*> h1_responseGEN_clusterOK = getResponseHistos("responseGEN_clusterOK");

  std::vector<TH1F*> h1_response_loose        = getResponseHistos("response_loose");
  std::vector<TH1F*> h1_responseGEN_loose     = getResponseHistos("responseGEN_loose");
  std::vector<TH1F*> h1_responseMPF_loose     = getResponseHistos("responseMPF_loose");



  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeight_loose;
  if( useGenJets )
    tree->SetBranchAddress("eventWeight_loose", &eventWeight_loose);
  Float_t eventWeight_medium;
  if( useGenJets )
    tree->SetBranchAddress("eventWeight_medium", &eventWeight_medium);

  Float_t eMet;
  Float_t phiMet;
  if( recoType == "pf" ) {
    tree->SetBranchAddress("epfMet", &eMet);
    tree->SetBranchAddress("phipfMet", &phiMet);
  } else if( recoType == "calo" ) {
    tree->SetBranchAddress("eMet", &eMet);
    tree->SetBranchAddress("phiMet", &phiMet);
  } else if( recoType == "jpt" ) {
    tree->SetBranchAddress("etcMet", &eMet);
    tree->SetBranchAddress("phitcMet", &phiMet);
  }

//Bool_t isIsolated_hcal_loose;
//tree->SetBranchAddress("isIsolated_hcal_loose", &isIsolated_hcal_loose);
//Bool_t isIsolated_ecal_loose;
//tree->SetBranchAddress("isIsolated_ecal_loose", &isIsolated_ecal_loose);
//Bool_t isIsolated_ptTracks_loose;
//tree->SetBranchAddress("isIsolated_ptTracks_loose", &isIsolated_ptTracks_loose);
//Bool_t isIsolated_nTracks_loose;
//tree->SetBranchAddress("isIsolated_nTracks_loose", &isIsolated_nTracks_loose);
//Bool_t clusterMajOK_loose;
//tree->SetBranchAddress("clusterMajOK_loose", &clusterMajOK_loose);
//Bool_t clusterMinOK_loose;
//tree->SetBranchAddress("clusterMinOK_loose", &clusterMinOK_loose);

//Bool_t isIsolated_hcal_medium;
//tree->SetBranchAddress("isIsolated_hcal_medium", &isIsolated_hcal_medium);
//Bool_t isIsolated_ecal_medium;
//tree->SetBranchAddress("isIsolated_ecal_medium", &isIsolated_ecal_medium);
//Bool_t isIsolated_ptTracks_medium;
//tree->SetBranchAddress("isIsolated_ptTracks_medium", &isIsolated_ptTracks_medium);
//Bool_t isIsolated_nTracks_medium;
//tree->SetBranchAddress("isIsolated_nTracks_medium", &isIsolated_nTracks_medium);
//Bool_t clusterMajOK_medium;
//tree->SetBranchAddress("clusterMajOK_medium", &clusterMajOK_medium);
//Bool_t clusterMinOK_medium;
//tree->SetBranchAddress("clusterMinOK_medium", &clusterMinOK_medium);

//Bool_t passedPhotonID_loose;
//tree->SetBranchAddress("passedPhotonID_loose", &passedPhotonID_loose);
//Bool_t passedPhotonID_medium;
//tree->SetBranchAddress("passedPhotonID_medium", &passedPhotonID_medium);

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

  Float_t eJetReco;
  tree->SetBranchAddress("eJetReco", &eJetReco);
  Float_t ptJetReco;
  tree->SetBranchAddress("ptJetReco", &ptJetReco);
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





  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;


    if( ONEVTX_ && dataset!="QCD_Spring10" ) {
      if( nvertex>1 ) continue;
    }


    bool jetInBarrel = (fabs(etaJetReco)<1.3);

    if( ptPhotReco<ptPhot_binning[0] ) continue;
    if( fabs(etaPhotReco)>1.3 ) continue;
    if( clusterMinPhotReco<0.15 && !MCassoc && !useGenJets ) continue; //protection vs EB spikes


    //first find correct photon pt bin:
    int theBin = hp_ptPhotMean->FindBin( ptPhotReco );
    theBin -= 1; //because arrays start from 0


    //leading jet and photon back2back in transverse plane
    bool back2back = true;
    Float_t deltaPhi_jet = fabs(fitTools::delta_phi(phiPhotReco, phiJetReco));
    Float_t pi = TMath::Pi();
    float deltaPhiThreshold = (TIGHTDELTAPHI_) ? 0.2 : 1.;
    if( fabs(deltaPhi_jet) < (pi - 1.) ) back2back = false; //loose back to back for now


    Float_t deltaPhi_2ndJet = fabs(fitTools::delta_phi(phiPhotReco, phi2ndJetReco));



    bool secondJetOK;
    if( NO2ndJETABS )
      secondJetOK = ( pt2ndJetReco < secondJetThreshold*ptPhotReco );
    else
      secondJetOK = ( pt2ndJetReco < secondJetThreshold*ptPhotReco || pt2ndJetReco < 5. );


    // do them by hand just to be sure:
    Bool_t isIsolated_hcal_loose = ( hcalIsoPhotReco<0.1 || hcalIsoPhotReco*ePhotReco<4. );
    Bool_t isIsolated_ecal_loose = ( ecalIsoPhotReco<0.1  || ecalIsoPhotReco*ePhotReco<4.5 );
    Bool_t isIsolated_ptTracks_loose = ( ptTrkIsoPhotReco<0.2 );
    Bool_t isIsolated_nTracks_loose = (nTrkIsoPhotReco < 5 );
    Bool_t clusterMajOK_loose = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
    Bool_t clusterMinOK_loose = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    Bool_t isIsolated_hcal_medium = ( hcalIsoPhotReco<0.05 || hcalIsoPhotReco*ePhotReco<2.4 );
    Bool_t isIsolated_ecal_medium = ( ecalIsoPhotReco<0.05  || ecalIsoPhotReco*ePhotReco<3. );
    Bool_t isIsolated_ptTracks_medium = ( ptTrkIsoPhotReco<0.1 );
    Bool_t isIsolated_nTracks_medium = (nTrkIsoPhotReco < 3 );
    Bool_t clusterMajOK_medium = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
    Bool_t clusterMinOK_medium = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

  ////temporary fix:
    isIsolated_ecal_loose  = ( ecalIsoPhotReco<0.1  || ecalIsoPhotReco*ePhotReco<4.5 );
    isIsolated_ecal_medium = ( ecalIsoPhotReco<0.05 || ecalIsoPhotReco*ePhotReco<3.  );
  //isIsolated_ecal_loose  = ( ecalIsoPhotReco<0.1  || ecalIsoPhotReco*ePhotReco<3. );
  //isIsolated_ecal_medium = ( ecalIsoPhotReco<0.05 || ecalIsoPhotReco*ePhotReco<1.7  );
  

    //before selection fill N-1 isolation plots (no event topology for isolation variables):
    if(                           isIsolated_ecal_medium  && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_hcalIsoPhotReco_Nm1->Fill( hcalIsoPhotReco, eventWeight);
    if(                           isIsolated_ecal_medium  && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_hcalIsoEnergyPhotReco_Nm1->Fill( hcalIsoPhotReco*ePhotReco, eventWeight);
    if( isIsolated_hcal_medium                           && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_ecalIsoPhotReco_Nm1->Fill( ecalIsoPhotReco, eventWeight);
    if( isIsolated_hcal_medium                           && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_ecalIsoEnergyPhotReco_Nm1->Fill( ecalIsoPhotReco*ePhotReco, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium                                && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_ptTrkIsoPhotReco_Nm1->Fill( ptTrkIsoPhotReco, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium  && isIsolated_ptTracks_medium                              && clusterMajOK_medium && clusterMinOK_medium  ) h1_nTrkIsoPhotReco_Nm1->Fill( nTrkIsoPhotReco, eventWeight);
    //no cluster cuts on cluster N-1's:
    if( isIsolated_hcal_medium && isIsolated_ecal_medium  && isIsolated_ptTracks_medium && isIsolated_nTracks_medium   ) h1_clusterMajPhotReco_Nm1->Fill( clusterMajPhotReco, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium  && isIsolated_ptTracks_medium && isIsolated_nTracks_medium   ) h1_clusterMinPhotReco_Nm1->Fill( clusterMinPhotReco, eventWeight);
    // yes topology for topology variables:
    if( isIsolated_hcal_medium && isIsolated_ecal_medium  && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium              && secondJetOK && jetInBarrel) h1_deltaPhi_Nm1->Fill( deltaPhi_jet, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium  && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium && back2back                && jetInBarrel) h1_ptSecondJetRel_Nm1->Fill( pt2ndJetReco/ptPhotReco, eventWeight);


    bool isIsolated_loose = (isIsolated_hcal_loose && isIsolated_ecal_loose && isIsolated_ptTracks_loose && isIsolated_nTracks_loose);
    bool isIsolated_medium = (isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium);
    bool clusterShapeOK_medium = (clusterMajOK_medium && clusterMinOK_medium );



    //////////////////////////////////////////////
    /////      CLUSTER SHAPE ONLY SELECTION 
    //////////////////////////////////////////////


    if( clusterShapeOK_medium ) {
      h1_hcalIsoPhotReco_clusterOK->Fill( hcalIsoPhotReco, eventWeight);
      h1_hcalIsoEnergyPhotReco_clusterOK->Fill( hcalIsoPhotReco*ePhotReco, eventWeight);
      h1_ecalIsoPhotReco_clusterOK->Fill( ecalIsoPhotReco, eventWeight);
      h1_ecalIsoEnergyPhotReco_clusterOK->Fill( ecalIsoPhotReco*ePhotReco, eventWeight);
      h1_ptTrkIsoPhotReco_clusterOK->Fill( ptTrkIsoPhotReco, eventWeight);
      h1_nTrkIsoPhotReco_clusterOK->Fill( nTrkIsoPhotReco, eventWeight);
      h1_deltaPhi_clusterOK->Fill( deltaPhi_jet, eventWeight);

      if( isIsolated_medium ) {
        h1_deltaPhi_clusterOK_isolated->Fill( deltaPhi_jet, eventWeight);
        h1_ptSecondJetRel_clusterOK_isolated->Fill( pt2ndJetReco/ptPhotReco, eventWeight);
        h1_phiPhot_clusterOK_isolated->Fill( phiPhotReco, eventWeight );
        h1_etaPhot_clusterOK_isolated->Fill( etaPhotReco, eventWeight );
        h1_ptPhot_clusterOK_isolated->Fill( ptPhotReco, eventWeight );
      }

      if( (back2back && secondJetOK && jetInBarrel) || noJetSelection ) {
        h1_ptPhot_clusterOK->Fill( ptPhotReco, eventWeight );
        h1_phiPhot_clusterOK->Fill( phiPhotReco, eventWeight );
        h1_etaPhot_clusterOK->Fill( etaPhotReco, eventWeight );
        h1_response_clusterOK[theBin]->Fill( ptJetReco/ptPhotReco, eventWeight );
        h1_responseGEN_clusterOK[theBin]->Fill( ptJetReco/ptJetGen, eventWeight );
      } //if back2back
    }

    if( isIsolated_medium ) {
      h1_clusterMajPhotReco_isolated->Fill( clusterMajPhotReco, eventWeight );
      h1_clusterMinPhotReco_isolated->Fill( clusterMinPhotReco, eventWeight );
    }

    h1_clusterMajPhotReco->Fill( clusterMajPhotReco, eventWeight );
    h1_clusterMinPhotReco->Fill( clusterMinPhotReco, eventWeight );


    bool photonOK_medium = (isIsolated_medium && clusterShapeOK_medium) || MCassoc || useGenJets;
    bool photonOK_loose  = (isIsolated_loose  && clusterShapeOK_medium) || MCassoc || useGenJets;

    // compute mpf :
    Float_t phi_Phot_Met = fitTools::delta_phi( phiPhotReco, phiMet );
    Float_t mpfResponse = 1. + eMet*ptPhotReco*cos( phi_Phot_Met ) / (ptPhotReco*ptPhotReco);

    

  
    //////////////////////////////////////////////
    /////      EVENT SELECTION: MEDIUM ID
    //////////////////////////////////////////////

    bool passedMedium_no2ndJet = photonOK_medium && ( (back2back && jetInBarrel) || noJetSelection );
    bool passedMedium_FULL     = passedMedium_no2ndJet && (secondJetOK || noJetSelection);


    if( passedMedium_no2ndJet ) {

      Float_t correctWeight = (useGenJets) ? eventWeight_medium : eventWeight;

    //// no cut in 2nd jet to produce gen response plot to have more stat
    //h1_responseGEN[theBin]->Fill( ptJetReco/ptJetGen, correctWeight );
      
      if( passedMedium_FULL ) {

        h1_ptPhot_medium->Fill( ptPhotReco, correctWeight );

        h1_deltaPhi_2ndJet_medium->Fill( deltaPhi_2ndJet, correctWeight );
        
        int theBinGEN = hp_ptPhotMean->FindBin( ptJetGen );
        theBinGEN -= 1; //because arrays start from 0
      

        //h1_responseGEN[theBin]->Fill( ptJetReco/ptJetGen, correctWeight );
        
        if( ptJetGen>ptPhot_binning[0] ) {
          h1_responseGEN[theBinGEN]->Fill( ptJetReco/ptJetGen, eventWeight );
          hp_ptJetGenMean->Fill( ptJetGen, ptJetGen, eventWeight );
        }

        h1_phiPhot->Fill( phiPhotReco, correctWeight );
        h1_etaPhot->Fill( etaPhotReco, correctWeight );
        h1_ptPhot->Fill( ptPhotReco, correctWeight );
        h1_ptJetReco->Fill( ptJetReco, correctWeight );
        h1_pt2ndJetReco->Fill( pt2ndJetReco, correctWeight );
        hp_ptPhotMean->Fill( ptPhotReco, ptPhotReco, correctWeight );
        h1_response[theBin]->Fill( ptJetReco/ptPhotReco, correctWeight );

   //   std::cout << "eMet: " << eMet << "\tptPhot: " << ptPhotReco << "\tphiMet: " << phiMet << "\tphiPhot: " << phiPhotReco << "\tdeltaPhi: " << phi_Phot_Met << "\tmpf: " << mpfResponse << std::endl;
//        if( recoType=="pf" ) {
          h1_deltaPhi_phot_met_medium->Fill( phi_Phot_Met, correctWeight );
          h1_met_medium->Fill( eMet, correctWeight );
          h1_responseMPF[theBin]->Fill( mpfResponse, correctWeight );
//        } //if pf
          
      } //if second jet ok

    } 


    //////////////////////////////////////////////
    /////      EVENT SELECTION: LOOSE ID
    //////////////////////////////////////////////

    bool passedLoose_no2ndJet = photonOK_loose && ( (back2back && jetInBarrel) || noJetSelection );
    bool passedLoose_FULL     = passedLoose_no2ndJet && (secondJetOK || noJetSelection);


    if( passedLoose_no2ndJet ) {

      Float_t correctWeight = (useGenJets) ? eventWeight_loose : eventWeight;
      
      
    
      if( passedLoose_FULL ) {

        h1_responseGEN_loose[theBin]->Fill( ptJetReco/ptJetGen, correctWeight );

        h1_ptPhot_loose->Fill( ptPhotReco, correctWeight );
        h1_phiPhot_loose->Fill( phiPhotReco, correctWeight );
        h1_etaPhot_loose->Fill( etaPhotReco, correctWeight );

        h1_deltaPhi_2ndJet_loose->Fill( deltaPhi_2ndJet, correctWeight );

        hp_ptPhotMean_loose->Fill( ptPhotReco, ptPhotReco, correctWeight );
        h1_response_loose[theBin]->Fill( ptJetReco/ptPhotReco, correctWeight );

//        if( recoType=="pf" ) {
          h1_deltaPhi_phot_met_loose->Fill( phi_Phot_Met, correctWeight );
          h1_met_loose->Fill( eMet, correctWeight );
          h1_responseMPF_loose[theBin]->Fill( mpfResponse, correctWeight );
//        } //if pf
      } //if second jet ok

    } //if loose


  } //for entries



  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaPhotonJet_"+dataset;
  else {
   if(dataset!="") outfileName = "PhotonJet_"+dataset;
   else outfileName = "PhotonJet";
  }

  outfileName = outfileName + suffix;
  if( noJetSelection ) outfileName = outfileName + "_NOJETSEL";
  if( MCassoc ) outfileName = outfileName + "_MCassoc";
  if( useGenJets ) outfileName = outfileName + "_GENJETS";
  if( BINNINGFINO_ ) outfileName = outfileName + "_BINNINGFINO";
  if( ONEVTX_ ) outfileName = outfileName + "_ONEVTX";
  if( TIGHTDELTAPHI_ ) outfileName = outfileName + "_TIGHTDPHI";
  if( !noJetSelection && secondJetThreshold!=0.5 ) {
    std::string R = ( NO2ndJETABS ) ? "R" : "";
    char outfileName_char[300];
    sprintf( outfileName_char, "%s_2ndJet%d%s", outfileName.c_str(), (int)(100.*secondJetThreshold), R.c_str());
    std::string outfileName_str( outfileName_char );
    outfileName = outfileName_str;
  }
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h1_totalLumi->Write();

  h1_ptPhot->Write();
  h1_phiPhot->Write();
  h1_etaPhot->Write();
  h1_ptJetReco->Write();
  h1_pt2ndJetReco->Write();

  h1_ptPhot_loose->Write();
  h1_phiPhot_loose->Write();
  h1_etaPhot_loose->Write();

  h1_ptPhot_medium->Write();

  h1_ptPhot_clusterOK->Write();
  h1_phiPhot_clusterOK->Write();
  h1_etaPhot_clusterOK->Write();

  h1_ptPhot_clusterOK_isolated->Write();
  h1_phiPhot_clusterOK_isolated->Write();
  h1_etaPhot_clusterOK_isolated->Write();
  
  h1_met_loose->Write();
  h1_met_medium->Write();

  h1_deltaPhi_phot_met_loose->Write();
  h1_deltaPhi_phot_met_medium->Write();

  h1_deltaPhi_Nm1->Write();
  h1_ptSecondJetRel_Nm1->Write();

  h1_deltaPhi_clusterOK_isolated->Write();
  h1_ptSecondJetRel_clusterOK_isolated->Write();

  h1_deltaPhi_2ndJet_medium->Write();
  h1_deltaPhi_2ndJet_loose->Write();

  h1_hcalIsoPhotReco_Nm1->Write();
  h1_hcalIsoEnergyPhotReco_Nm1->Write();
  h1_ecalIsoPhotReco_Nm1->Write();
  h1_ecalIsoEnergyPhotReco_Nm1->Write();
  h1_ptTrkIsoPhotReco_Nm1->Write();
  h1_nTrkIsoPhotReco_Nm1->Write();
  h1_clusterMajPhotReco_Nm1->Write();
  h1_clusterMinPhotReco_Nm1->Write();

  h1_hcalIsoPhotReco_clusterOK->Write();
  h1_hcalIsoEnergyPhotReco_clusterOK->Write();
  h1_ecalIsoPhotReco_clusterOK->Write();
  h1_ecalIsoEnergyPhotReco_clusterOK->Write();
  h1_ptTrkIsoPhotReco_clusterOK->Write();
  h1_nTrkIsoPhotReco_clusterOK->Write();
  h1_deltaPhi_clusterOK->Write();

  h1_clusterMajPhotReco_isolated->Write();
  h1_clusterMinPhotReco_isolated->Write();

  h1_clusterMajPhotReco->Write();
  h1_clusterMinPhotReco->Write();

  hp_ptJetGenMean->Write();
  hp_ptPhotMean->Write();
  hp_ptPhotMean_loose->Write();

  for( unsigned i=0; i<h1_response.size(); ++i ) {
    h1_response[i]->Write();
    h1_responseGEN[i]->Write();
    h1_response_loose[i]->Write();
    h1_responseGEN_loose[i]->Write();
    h1_response_clusterOK[i]->Write();
    h1_responseGEN_clusterOK[i]->Write();

//    if( recoType=="pf" ) {
      h1_responseMPF[i]->Write();
      h1_responseMPF_loose[i]->Write();
//    }
  }


  outFile->Close();

  delete h1_hcalIsoPhotReco_Nm1;
  h1_hcalIsoPhotReco_Nm1 = 0;
  delete h1_hcalIsoEnergyPhotReco_Nm1;
  h1_hcalIsoEnergyPhotReco_Nm1 = 0;
  delete h1_ecalIsoPhotReco_Nm1;
  h1_ecalIsoPhotReco_Nm1 = 0;
  delete h1_ecalIsoEnergyPhotReco_Nm1;
  h1_ecalIsoEnergyPhotReco_Nm1 = 0;
  delete h1_ptTrkIsoPhotReco_Nm1;
  h1_ptTrkIsoPhotReco_Nm1 = 0;
  delete h1_nTrkIsoPhotReco_Nm1;
  h1_nTrkIsoPhotReco_Nm1 = 0;
  delete h1_clusterMajPhotReco_Nm1;
  h1_clusterMajPhotReco_Nm1 = 0;
  delete h1_clusterMinPhotReco_Nm1;
  h1_clusterMinPhotReco_Nm1 = 0;

  delete h1_ptJetReco;
  h1_ptJetReco = 0;
  delete h1_pt2ndJetReco;
  h1_pt2ndJetReco = 0;
  delete h1_ptPhot;
  h1_ptPhot = 0;
  delete h1_etaPhot;
  h1_etaPhot = 0;
  delete h1_phiPhot;
  h1_phiPhot = 0;
  delete h1_deltaPhi_Nm1;
  h1_deltaPhi_Nm1 = 0;
  delete h1_ptSecondJetRel_Nm1;
  h1_ptSecondJetRel_Nm1 = 0;

  delete h1_ptPhot_loose;
  h1_ptPhot_loose = 0;
  delete h1_etaPhot_loose;
  h1_etaPhot_loose = 0;
  delete h1_phiPhot_loose;
  h1_phiPhot_loose = 0;

  delete h1_ptPhot_medium;
  h1_ptPhot_medium = 0;

  delete h1_hcalIsoPhotReco_clusterOK;
  h1_hcalIsoPhotReco_clusterOK = 0;
  delete h1_hcalIsoEnergyPhotReco_clusterOK;
  h1_hcalIsoEnergyPhotReco_clusterOK = 0;
  delete h1_ecalIsoPhotReco_clusterOK;
  h1_ecalIsoPhotReco_clusterOK = 0;
  delete h1_ecalIsoEnergyPhotReco_clusterOK;
  h1_ecalIsoEnergyPhotReco_clusterOK = 0;
  delete h1_ptTrkIsoPhotReco_clusterOK;
  h1_ptTrkIsoPhotReco_clusterOK = 0;
  delete h1_nTrkIsoPhotReco_clusterOK;
  h1_nTrkIsoPhotReco_clusterOK = 0;

  delete h1_clusterMajPhotReco_isolated;
  h1_clusterMajPhotReco_isolated = 0;
  delete h1_clusterMinPhotReco_isolated;
  h1_clusterMinPhotReco_isolated = 0;
//delete h1_response;
//h1_response = 0;

  delete h1_totalLumi;
  h1_totalLumi = 0;
  delete h1_deltaPhi_2ndJet_medium;
  h1_deltaPhi_2ndJet_medium= 0;
  delete h1_deltaPhi_2ndJet_loose;
  h1_deltaPhi_2ndJet_loose= 0;
  delete h1_met_loose;
  h1_met_loose= 0;
  delete h1_met_medium;
  h1_met_medium;
  delete h1_deltaPhi_phot_met_loose;
  h1_deltaPhi_phot_met_loose = 0;
  delete h1_deltaPhi_phot_met_medium;
  h1_deltaPhi_phot_met_medium = 0;
  delete h1_ptPhot_clusterOK;
  h1_ptPhot_clusterOK = 0;
  delete h1_etaPhot_clusterOK;
  h1_etaPhot_clusterOK = 0;
  delete h1_phiPhot_clusterOK;
  h1_phiPhot_clusterOK = 0;
  delete h1_deltaPhi_clusterOK;
  h1_deltaPhi_clusterOK = 0;
  delete h1_ptPhot_clusterOK_isolated;
  h1_ptPhot_clusterOK_isolated = 0;
  delete h1_etaPhot_clusterOK_isolated;
  h1_etaPhot_clusterOK_isolated = 0;
  delete h1_phiPhot_clusterOK_isolated;
  h1_phiPhot_clusterOK_isolated = 0;
  delete h1_clusterMajPhotReco;
  h1_clusterMajPhotReco = 0;
  delete h1_clusterMinPhotReco;
  h1_clusterMinPhotReco = 0;
  delete h1_deltaPhi_clusterOK_isolated;
  h1_deltaPhi_clusterOK_isolated = 0;
  delete h1_ptSecondJetRel_clusterOK_isolated;
  h1_ptSecondJetRel_clusterOK_isolated = 0;
  delete hp_ptJetGenMean;
  hp_ptJetGenMean = 0;
  delete hp_ptPhotMean;
  hp_ptPhotMean = 0;
  delete hp_ptPhotMean_loose;
  hp_ptPhotMean_loose = 0;

  for( unsigned i=0; i< h1_response.size(); ++i) {
    delete h1_response[i];
    h1_response[i]=0;
    delete h1_responseGEN[i];
    h1_responseGEN[i]=0;
    delete h1_responseMPF[i];
    h1_responseMPF[i]=0;

    delete h1_response_clusterOK[i];
    h1_response_clusterOK[i]=0;
    delete h1_responseGEN_clusterOK[i];
    h1_responseGEN_clusterOK[i]=0;

    delete h1_response_loose[i];
    h1_response_loose[i]=0;
    delete h1_responseGEN_loose[i];
    h1_responseGEN_loose[i]=0;
    delete h1_responseMPF_loose[i];
    h1_responseMPF_loose[i]=0;
  }

  delete tree;
  tree = 0;

  totalLumi = 0.;

}


void addInput( const std::string& dataset ) {

  std::string infileName = "files_PhotonJet_2ndLevel_" + dataset+"_" + ALGOTYPE_ +".txt";
  TH1F* h1_lumi;


  //open from file.txt:
  FILE* iff = fopen(infileName.c_str(),"r");
  if(iff == 0) {
    std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    infileName = "PhotonJet_2ndLevelTree_" + dataset + suffix + ".root";
    std::string treeName = infileName +"/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    TFile* infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
      std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    infile->Close();

  } else {

    char singleLine[500];

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/jetTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
      }
      infile->Close();

    }
    fclose(iff);

  }

} //addinput


  
std::vector<TH1F*> getResponseHistos(const std::string& name) {

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();
  std::vector<TH1F*> returnVector;

  for( unsigned i=0; i<(ptPhot_binning.size()-1); ++i ) {
    char histoName[100];
    sprintf( histoName, "%s_ptPhot_%.0f_%.0f", name.c_str(), ptPhot_binning[i], ptPhot_binning[i+1]);
    int nbins = (BINNINGFINO_) ? 1000 : 14;
    float xmin = (BINNINGFINO_) ? -5. : 0.1;
    float xmax = (BINNINGFINO_) ? 10. : 2.;
    //TH1F* newHisto = new TH1F(histoName, "", 15, 0., 2.);
    TH1F* newHisto = new TH1F(histoName, "", nbins, xmin, xmax);
    newHisto->Sumw2();
    returnVector.push_back(newHisto);
  }

  return returnVector;

}

