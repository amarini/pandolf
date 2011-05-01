
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "fitTools.C"




bool DEBUG_ = false;
bool MCassoc_ = false;
bool BINNINGFINO_ = false;
bool ONEVTX_ = false;
bool ONEVTX_2232_ = false;
bool ONEVTX_2247_ = false;
bool NO2ndJETABS = false;
bool ADD12_ = false;
bool HEMISPHERE_ = false;
float DELTAPHI_ = 1.;
std::string PARTTYPE_;
std::string suffix;

TChain* tree;
Double_t totalLumi=0.;


void addInput(const std::string& dataset, bool genjets=false);




void finalize(const std::string& dataset, std::string photonID="medium", float secondJetThreshold=0.1, bool useGenJets=false, std::string partType="") {

//if( dataset=="PhotonJet_Summer1036X" )
//  MCassoc_ = true;

  //float secondJetThreshold=0.1;

  bool noJetSelection = ( secondJetThreshold < 0. );

  tree = new TChain("jetTree");


  PARTTYPE_ = partType;

  std::string infileName, treeName;

  if( photonID!="medium" && photonID!="loose" && photonID!="clusterOK") {
    std::cout << "Photon ID: '" << photonID << "' currently not implemented. Exiting." << std::endl;
    exit(11);
  }

  if( useGenJets ) suffix = suffix + "_GENJETS";



  if( dataset=="G_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" ) {

    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );

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

  } else if( dataset=="QCD_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" ) {

    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );

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

  } else if( dataset=="G_QCD_TuneZ2_7TeV_pythia6_GENJETS" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6", (bool)true );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6", (bool)true );

  } else if( dataset=="G_QCD_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data_BX156_final_38" );
    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_Fall10_ProbDist_2010Data" );

  } else if( dataset=="DATA_Nov4ReReco" ) {

    addInput( "EG_Run2010A-Nov4ReReco_v1" );
    addInput( "Photon_Run2010B_Nov4ReReco" );

  } else if( dataset=="DATA_Nov4ReReco_L1Offset" ) {

    addInput( "EG_Run2010A-Nov4ReReco_L1Offset" );
    addInput( "Photon_Run2010B-Nov4ReReco_L1Offset" );

  } else if( dataset=="DATA_Nov4ReReco_L1Offset_ALL" ) {

    addInput( "EG_Run2010A-Nov4ReReco_L1Offset_ALL" );
    addInput( "Photon_Run2010B-Nov4ReReco_L1Offset" );

  } else if( dataset=="DATA_Nov4ReReco_final" ) {

    addInput( "EG_Run2010A-Nov4ReReco_final" );
    addInput( "Photon_Run2010B-Nov4ReReco_final" );

  } else {

    addInput( dataset );

  }



  //std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);

  TH1F* h1_nvertex = new TH1F("nvertex", "", 11, -0.5, 10.5);
  h1_nvertex->Sumw2();
  TH1F* h1_nvertex_passedID = new TH1F("nvertex_passedID", "", 11, -0.5, 10.5);
  h1_nvertex_passedID->Sumw2();
  TH1F* h1_nvertex_passedID_pt50 = new TH1F("nvertex_passedID_pt50", "", 11, -0.5, 10.5);
  h1_nvertex_passedID_pt50->Sumw2();

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Float_t ptMax = 100.;
  int nBins_photBinning = ptPhot_binning.size()-1;
  Double_t ptPhot_binning_array[100];
  for( unsigned i=0; i<nBins_photBinning; ++i)
    ptPhot_binning_array[i] = ptPhot_binning[i];




  TH1D* h1_ptPhot= new TH1D("ptPhot", "", 200, 50., 250.);
  h1_ptPhot->Sumw2();
  TH1D* h1_etaPhot= new TH1D("etaPhot", "", 15, -1.3, 1.3);
  h1_etaPhot->Sumw2();
  TH1D* h1_phiPhot= new TH1D("phiPhot", "", 15, -3.1416, 3.1416);
  h1_phiPhot->Sumw2();

  TH1D* h1_nJets= new TH1D("nJets", "", 7, -0.5, 6.5);
  h1_nJets->Sumw2();

  TH1D* h1_ptPhot_0jet = new TH1D("ptPhot_0jet", "", 200, 50., 250.);
  h1_ptPhot_0jet->Sumw2();
  TH1D* h1_ptPhot_1jet = new TH1D("ptPhot_1jet", "", 200, 50., 250.);
  h1_ptPhot_1jet->Sumw2();
  TH1D* h1_ptPhot_2jet = new TH1D("ptPhot_2jet", "", 200, 50., 250.);
  h1_ptPhot_2jet->Sumw2();
  TH1D* h1_ptPhot_3jet = new TH1D("ptPhot_3jet", "", 200, 50., 250.);
  h1_ptPhot_3jet->Sumw2();
  TH1D* h1_ptPhot_4jet = new TH1D("ptPhot_4jet", "", 200, 50., 250.);
  h1_ptPhot_4jet->Sumw2();

  TH1D* h1_mJet_1jet = new TH1D("mJet_1jet", "", 150, 0., 150.);
  h1_mJet_1jet->Sumw2();
  TH1D* h1_mJet_2jet = new TH1D("mJet_2jet", "", 150, 0., 150.);
  h1_mJet_2jet->Sumw2();
  TH1D* h1_mJet_3jet = new TH1D("mJet_3jet", "", 150, 0., 150.);
  h1_mJet_3jet->Sumw2();
  TH1D* h1_mJet_4jet = new TH1D("mJet_4jet", "", 150, 0., 150.);
  h1_mJet_4jet->Sumw2();
  
  TH1D* h1_mGammaJet_1jet = new TH1D("mGammaJet_1jet", "", 500, 0., 500.);
  h1_mGammaJet_1jet->Sumw2();
  TH1D* h1_mGammaJet_2jet = new TH1D("mGammaJet_2jet", "", 500, 0., 500.);
  h1_mGammaJet_2jet->Sumw2();
  TH1D* h1_mGammaJet_3jet = new TH1D("mGammaJet_3jet", "", 500, 0., 500.);
  h1_mGammaJet_3jet->Sumw2();
  TH1D* h1_mGammaJet_4jet = new TH1D("mGammaJet_4jet", "", 500, 0., 500.);
  h1_mGammaJet_4jet->Sumw2();
  
  TH1D* h1_mGammaJet_eta2_1jet = new TH1D("mGammaJet_eta2_1jet", "", 500, 0., 500.);
  h1_mGammaJet_eta2_1jet->Sumw2();
  TH1D* h1_mGammaJet_eta2_2jet = new TH1D("mGammaJet_eta2_2jet", "", 500, 0., 500.);
  h1_mGammaJet_eta2_2jet->Sumw2();
  TH1D* h1_mGammaJet_eta2_3jet = new TH1D("mGammaJet_eta2_3jet", "", 500, 0., 500.);
  h1_mGammaJet_eta2_3jet->Sumw2();
  TH1D* h1_mGammaJet_eta2_4jet = new TH1D("mGammaJet_eta2_4jet", "", 500, 0., 500.);
  h1_mGammaJet_eta2_4jet->Sumw2();
  
  TH1D* h1_mtGammaJet_1jet = new TH1D("mtGammaJet_1jet", "", 500, 0., 500.);
  h1_mtGammaJet_1jet->Sumw2();
  TH1D* h1_mtGammaJet_2jet = new TH1D("mtGammaJet_2jet", "", 500, 0., 500.);
  h1_mtGammaJet_2jet->Sumw2();
  TH1D* h1_mtGammaJet_3jet = new TH1D("mtGammaJet_3jet", "", 500, 0., 500.);
  h1_mtGammaJet_3jet->Sumw2();
  TH1D* h1_mtGammaJet_4jet = new TH1D("mtGammaJet_4jet", "", 500, 0., 500.);
  h1_mtGammaJet_4jet->Sumw2();
  
  TH1D* h1_met_1jet = new TH1D("met_1jet", "", 250, 0., 250.);
  h1_met_1jet->Sumw2();
  TH1D* h1_met_2jet = new TH1D("met_2jet", "", 250, 0., 250.);
  h1_met_2jet->Sumw2();
  TH1D* h1_met_3jet = new TH1D("met_3jet", "", 250, 0., 250.);
  h1_met_3jet->Sumw2();
  TH1D* h1_met_4jet = new TH1D("met_4jet", "", 250, 0., 250.);
  h1_met_4jet->Sumw2();






  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeight_genjets;
  if( useGenJets ) {
    if( photonID=="medium" )
      tree->SetBranchAddress("eventWeight_medium", &eventWeight_genjets);
    else if( photonID=="loose" )
      tree->SetBranchAddress("eventWeight_loose", &eventWeight_genjets);
    else {
      std::cout << "Photon ID: '" << photonID << "' currently not implemented for GENJET technology. Exiting." << std::endl;
      exit(12);
    }
  }
//Float_t eventWeight_medium;
//if( useGenJets )
//  tree->SetBranchAddress("eventWeight_medium", &eventWeight_medium);

  Float_t eMet;
  Float_t phiMet;
    tree->SetBranchAddress("epfMet", &eMet);
    tree->SetBranchAddress("phipfMet", &phiMet);


  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t ePhot;
  tree->SetBranchAddress("ePhotReco", &ePhot);
  Float_t ptPhot;
  tree->SetBranchAddress("ptPhotReco", &ptPhot);
  Float_t etaPhot;
  tree->SetBranchAddress("etaPhotReco", &etaPhot);
  Float_t phiPhot;
  tree->SetBranchAddress("phiPhotReco", &phiPhot);

  Float_t ePhotGen;
  tree->SetBranchAddress("ePhotGen", &ePhotGen);
  Float_t ptPhotGen;
  tree->SetBranchAddress("ptPhotGen", &ptPhotGen);
  Float_t etaPhotGen;
  tree->SetBranchAddress("etaPhotGen", &etaPhotGen);
  Float_t phiPhotGen;
  tree->SetBranchAddress("phiPhotGen", &phiPhotGen);

  Float_t hcalIsoPhot;
  tree->SetBranchAddress("hcalIsoPhotReco", &hcalIsoPhot);
  Float_t ecalIsoPhot;
  tree->SetBranchAddress("ecalIsoPhotReco", &ecalIsoPhot);
  Int_t nTrkIsoPhot;
  tree->SetBranchAddress("nTrkIsoPhotReco", &nTrkIsoPhot);
  Float_t ptTrkIsoPhot;
  tree->SetBranchAddress("ptTrkIsoPhotReco", &ptTrkIsoPhot);
  Float_t clusterMajPhot;
  tree->SetBranchAddress("clusterMajPhotReco", &clusterMajPhot);
  Float_t clusterMinPhot;
  tree->SetBranchAddress("clusterMinPhotReco", &clusterMinPhot);
  Int_t hasPixelSeedPhot;
  tree->SetBranchAddress("hasPixelSeedPhotReco", &hasPixelSeedPhot);

  Bool_t matchedToMC;
  tree->SetBranchAddress("matchedToMC", &matchedToMC);

  Int_t nJet;
  tree->SetBranchAddress("nJet", &nJet);

  Float_t eJet[20];
  tree->SetBranchAddress("eJet", eJet);
  Float_t ptJet[20];
  tree->SetBranchAddress("ptJet", ptJet);
  Float_t ptCorrJet[20];
  tree->SetBranchAddress("ptCorrJet", ptCorrJet);
  Float_t etaJet[20];
  tree->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[20];
  tree->SetBranchAddress("phiJet", phiJet);

  Float_t eChargedHadronsJet[20];
  tree->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t eNeutralHadronsJet[20];
  tree->SetBranchAddress("eNeutralHadronsJet", eNeutralHadronsJet);
  Float_t ePhotonsJet[20];
  tree->SetBranchAddress("ePhotonsJet", ePhotonsJet);
  Float_t eElectronsJet[20];
  tree->SetBranchAddress("eElectronsJet", eElectronsJet);
  Float_t eMuonsJet[20];
  tree->SetBranchAddress("eMuonsJet", eMuonsJet);

  Int_t nChargedHadronsJet[20];
  tree->SetBranchAddress("nChargedHadronsJet", nChargedHadronsJet);
  Int_t nNeutralHadronsJet[20];
  tree->SetBranchAddress("nNeutralHadronsJet", nNeutralHadronsJet);
  Int_t nPhotonsJet[20];
  tree->SetBranchAddress("nPhotonsJet", nPhotonsJet);
  Int_t nElectronsJet[20];
  tree->SetBranchAddress("nElectronsJet", nElectronsJet);
  Int_t nMuonsJet[20];
  tree->SetBranchAddress("nMuonsJet", nMuonsJet);


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





  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;
    Float_t correctWeight = (useGenJets) ? eventWeight_genjets : eventWeight;

    h1_nvertex->Fill( nvertex, correctWeight);

    if( ONEVTX_ && dataset!="QCD_Spring10" ) {
      if( nvertex>1 ) continue;
    }



    if( ptPhot<50. ) continue;
    if( fabs(etaPhot)>1.3 ) continue;
    if( clusterMinPhot<0.15 && !MCassoc_ && !useGenJets ) continue; //protection vs EB spikes



    // do them by hand just to be sure:
    Bool_t isIsolated_hcal;
    Bool_t isIsolated_ecal;
    Bool_t isIsolated_ptTracks;
    Bool_t isIsolated_nTracks;
    Bool_t clusterMajOK;
    Bool_t clusterMinOK;
    Bool_t pixelSeedOK;

    if( photonID=="medium" ) {

      isIsolated_hcal = ( hcalIsoPhot<0.05 || hcalIsoPhot*ePhot<2.4 );
      isIsolated_ecal = ( ecalIsoPhot<0.05  || ecalIsoPhot*ePhot<3. );
      isIsolated_ptTracks = ( ptTrkIsoPhot<0.1 );
      isIsolated_nTracks = (nTrkIsoPhot < 3 );
      clusterMajOK = ( clusterMajPhot>0.15 && clusterMajPhot<0.35 );
      clusterMinOK = ( clusterMinPhot>0.15 && clusterMinPhot<0.3 );
      pixelSeedOK = !hasPixelSeedPhot;
      if( dataset=="G_TuneZ2_7TeV_pythia6_CORR" || dataset=="QCD_TuneZ2_7TeV_pythia6" ) pixelSeedOK = true;

    } else if( photonID=="loose" ) {

      isIsolated_hcal = ( hcalIsoPhot<0.1 || hcalIsoPhot*ePhot<4. );
      isIsolated_ecal = ( ecalIsoPhot<0.1  || ecalIsoPhot*ePhot<4.5 );
      isIsolated_ptTracks = ( ptTrkIsoPhot<0.2 );
      isIsolated_nTracks = (nTrkIsoPhot < 5 );
      clusterMajOK = ( clusterMajPhot>0.15 && clusterMajPhot<0.35 );
      clusterMinOK = ( clusterMinPhot>0.15 && clusterMinPhot<0.3 );

    } else if( photonID=="clusterOK" ) {

      isIsolated_hcal = true;
      isIsolated_ecal = true;
      isIsolated_ptTracks = true;
      isIsolated_nTracks = true;
      clusterMajOK = ( clusterMajPhot>0.15 && clusterMajPhot<0.35 );
      clusterMinOK = ( clusterMinPhot>0.15 && clusterMinPhot<0.3 );

    }

  



    bool isIsolated = (isIsolated_hcal && isIsolated_ecal && isIsolated_ptTracks && isIsolated_nTracks);
    bool clusterShapeOK = (clusterMajOK && clusterMinOK );

    if( MCassoc_ && matchedToMC ) {
      isIsolated = true;
      clusterShapeOK= true;
    }

  
    //////////////////////////////////////////////
    /////      EVENT SELECTION
    //////////////////////////////////////////////


    bool photonOK = (isIsolated && clusterShapeOK && pixelSeedOK) || MCassoc_ || useGenJets;

    TLorentzVector photon;
    photon.SetPtEtaPhiE( ptPhot, etaPhot, phiPhot, ePhot );

    std::vector<int> jetIndexes;

    if( photonOK ) {

      h1_ptPhot->Fill( photon.Pt(), eventWeight );
      h1_etaPhot->Fill( photon.Eta(), eventWeight );
      h1_phiPhot->Fill( photon.Phi(), eventWeight );

      for( unsigned iJet=0; iJet<nJet; ++iJet ) {

        if( iJet >=10 ) continue;

        float corr = ptCorrJet[iJet]/ptJet[iJet];
        TLorentzVector thisJet;
        thisJet.SetPtEtaPhiE( ptCorrJet[iJet], etaJet[iJet], phiJet[iJet], corr*eJet[iJet] );

        if( thisJet.DeltaR(photon) > 0.5 && thisJet.Pt() > 20. ) {
          jetIndexes.push_back(iJet);
        }

      } //for jets

      //if( jetIndexes.size()==0 ) continue;

      h1_nJets->Fill( jetIndexes.size(), eventWeight );

      if( jetIndexes.size()==0 ) {
        h1_ptPhot_0jet->Fill( photon.Pt(), eventWeight );
      } else if( jetIndexes.size()==1 ) {
        h1_ptPhot_1jet->Fill( photon.Pt(), eventWeight ); 
        for( unsigned iJet=0; iJet<jetIndexes.size(); ++iJet ) {
          int thisJetIndex = jetIndexes[iJet];
          float corr = ptCorrJet[thisJetIndex]/ptJet[thisJetIndex];
          TLorentzVector thisJet;
          thisJet.SetPtEtaPhiE( ptCorrJet[thisJetIndex], etaJet[thisJetIndex], phiJet[thisJetIndex], corr*eJet[thisJetIndex] );
          h1_mJet_1jet->Fill( thisJet.M(), eventWeight );
          TLorentzVector gammaJet = photon+thisJet;
          h1_mGammaJet_1jet->Fill( gammaJet.M(), eventWeight );
          if( fabs(thisJet.Eta())>2. ) h1_mGammaJet_eta2_1jet->Fill( gammaJet.M(), eventWeight );
          h1_mtGammaJet_1jet->Fill( gammaJet.Mt(), eventWeight );
          h1_met_1jet->Fill( eMet, eventWeight );
        }
      } else if( jetIndexes.size()==2 ) {
        h1_ptPhot_2jet->Fill( photon.Pt(), eventWeight ); 
        for( unsigned iJet=0; iJet<jetIndexes.size(); ++iJet ) {
          int thisJetIndex = jetIndexes[iJet];
          float corr = ptCorrJet[thisJetIndex]/ptJet[thisJetIndex];
          TLorentzVector thisJet;
          thisJet.SetPtEtaPhiE( ptCorrJet[thisJetIndex], etaJet[thisJetIndex], phiJet[thisJetIndex], corr*eJet[thisJetIndex] );
          h1_mJet_2jet->Fill( thisJet.M(), eventWeight );
          TLorentzVector gammaJet = photon+thisJet;
          h1_mGammaJet_2jet->Fill( gammaJet.M(), eventWeight );
          if( fabs(thisJet.Eta())>2. ) h1_mGammaJet_eta2_2jet->Fill( gammaJet.M(), eventWeight );
          h1_mtGammaJet_2jet->Fill( gammaJet.Mt(), eventWeight );
          h1_met_2jet->Fill( eMet, eventWeight );
        }
      } else if( jetIndexes.size()==3 ) {
        h1_ptPhot_3jet->Fill( photon.Pt(), eventWeight ); 
        for( unsigned iJet=0; iJet<jetIndexes.size(); ++iJet ) {
          int thisJetIndex = jetIndexes[iJet];
          float corr = ptCorrJet[thisJetIndex]/ptJet[thisJetIndex];
          TLorentzVector thisJet;
          thisJet.SetPtEtaPhiE( ptCorrJet[thisJetIndex], etaJet[thisJetIndex], phiJet[thisJetIndex], corr*eJet[thisJetIndex] );
          h1_mJet_3jet->Fill( thisJet.M(), eventWeight );
          TLorentzVector gammaJet = photon+thisJet;
          h1_mGammaJet_3jet->Fill( gammaJet.M(), eventWeight );
          if( fabs(thisJet.Eta())>2. ) h1_mGammaJet_eta2_3jet->Fill( gammaJet.M(), eventWeight );
          h1_mtGammaJet_3jet->Fill( gammaJet.Mt(), eventWeight );
          h1_met_3jet->Fill( eMet, eventWeight );
        }
      } else if( jetIndexes.size()==4 ) {
        h1_ptPhot_4jet->Fill( photon.Pt(), eventWeight ); 
        for( unsigned iJet=0; iJet<jetIndexes.size(); ++iJet ) {
          int thisJetIndex = jetIndexes[iJet];
          float corr = ptCorrJet[thisJetIndex]/ptJet[thisJetIndex];
          TLorentzVector thisJet;
          thisJet.SetPtEtaPhiE( ptCorrJet[thisJetIndex], etaJet[thisJetIndex], phiJet[thisJetIndex], corr*eJet[thisJetIndex] );
          h1_mJet_4jet->Fill( thisJet.M(), eventWeight );
          TLorentzVector gammaJet = photon+thisJet;
          h1_mGammaJet_4jet->Fill( gammaJet.M(), eventWeight );
          if( fabs(thisJet.Eta())>2. ) h1_mGammaJet_eta2_4jet->Fill( gammaJet.M(), eventWeight );
          h1_mtGammaJet_4jet->Fill( gammaJet.Mt(), eventWeight );
          h1_met_4jet->Fill( eMet, eventWeight );
        }
      }

    } //if photonOK


  } //for entries



  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaPhotonJet_"+dataset;
  else {
   if(dataset!="") outfileName = "NP_"+dataset;
   else outfileName = "NP";
  }

  outfileName = outfileName + suffix;
  if( noJetSelection ) outfileName = outfileName + "_NOJETSEL";
  if( photonID!="medium" ) outfileName = outfileName + "_" + photonID;
  if( MCassoc_ ) outfileName = outfileName + "_MCassoc";
  //if( useGenJets ) outfileName = outfileName + "_GENJETS";
  if( PARTTYPE_!="" ) outfileName = outfileName + "_" + PARTTYPE_;
  if( BINNINGFINO_ ) outfileName = outfileName + "_BINNINGFINO";
  if( ONEVTX_ ) outfileName = outfileName + "_ONEVTX";
  if( ONEVTX_2232_ ) outfileName = outfileName + "_ONEVTX2232";
  if( ONEVTX_2247_ ) outfileName = outfileName + "_ONEVTX2247";
  if( DELTAPHI_ != 1. ) {
    char outfileName_char[300];
    sprintf( outfileName_char, "%s_DPHI%d", outfileName.c_str(), (int)(10.*DELTAPHI_));
    std::string outfileName_str( outfileName_char );
    outfileName = outfileName_str;
  }
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h1_totalLumi->Write();
  h1_nvertex->Write();
  h1_nvertex_passedID->Write();
  h1_nvertex_passedID_pt50->Write();

  h1_ptPhot->Write();
  h1_etaPhot->Write();
  h1_phiPhot->Write();

  h1_nJets->Write();

  h1_ptPhot_0jet->Write();
  h1_ptPhot_1jet->Write();
  h1_ptPhot_2jet->Write();
  h1_ptPhot_3jet->Write();
  h1_ptPhot_4jet->Write();

  h1_mJet_1jet->Write();
  h1_mJet_2jet->Write();
  h1_mJet_3jet->Write();
  h1_mJet_4jet->Write();
  
  h1_mGammaJet_1jet->Write();
  h1_mGammaJet_2jet->Write();
  h1_mGammaJet_3jet->Write();
  h1_mGammaJet_4jet->Write();
  
  h1_mGammaJet_eta2_1jet->Write();
  h1_mGammaJet_eta2_2jet->Write();
  h1_mGammaJet_eta2_3jet->Write();
  h1_mGammaJet_eta2_4jet->Write();
  
  h1_mtGammaJet_1jet->Write();
  h1_mtGammaJet_2jet->Write();
  h1_mtGammaJet_3jet->Write();
  h1_mtGammaJet_4jet->Write();
  
  h1_met_1jet->Write();
  h1_met_2jet->Write();
  h1_met_3jet->Write();
  h1_met_4jet->Write();
  
  outFile->Close();

  delete tree;
  tree = 0;

  totalLumi = 0.;

}


void addInput( const std::string& dataset, bool genjets ) {


  // opening from filelist now deprecated (files have to be merged with merge_and_setWeights
  TH1F* h1_lumi;


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "NP_2ndLevelTreeW_" + dataset + suffix;
    if( genjets ) infileName += "_GENJETS";
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
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
      std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    infile->Close();

  /*} else {

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

  } */

} //addinput


  
