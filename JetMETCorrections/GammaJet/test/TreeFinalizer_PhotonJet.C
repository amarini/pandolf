
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
bool MCassoc_ = false;
bool BINNINGFINO_ = false;
bool ONEVTX_ = false;
bool ONEVTX_2232_ = false;
bool ONEVTX_2247_ = false;
bool NO2ndJETABS = false;
bool ADD12_ = false;
bool HEMISPHERE_ = false;
float DELTAPHI_ = 1.;
std::string RECOTYPE_;
std::string ALGOTYPE_;
std::string PARTTYPE_;
std::string suffix;

TChain* tree;
Double_t totalLumi=0.;


void addInput(const std::string& dataset, bool genjets=false);
std::vector<TH1F*> getResponseHistos(const std::string& name);
std::vector< std::vector< TH1D* > > getExtrapHistoVector(const std::string& name, const std::string& abs_rel, Int_t nPoints, bool isL2L3=false, bool is_pt=false);
void deleteExtrapHistoVector(std::vector< std::vector< TH1D* > > histoVector, int nPoints);




void finalize(const std::string& dataset, std::string recoType, std::string jetAlgo="akt5", std::string photonID="medium", bool useGenJets=false, std::string partType="") {

//if( dataset=="PhotonJet_Summer1036X" )
//  MCassoc_ = true;

  float secondJetThreshold=0.1;

  bool noJetSelection = ( secondJetThreshold < 0. );

  tree = new TChain("jetTree");

  RECOTYPE_ = recoType;

  ALGOTYPE_ = (recoType=="calo") ? jetAlgo : recoType+jetAlgo;
  if( recoType=="jpt"&&jetAlgo=="akt5" ) ALGOTYPE_="jptak5"; 

  PARTTYPE_ = partType;

  std::string infileName, treeName;

  if( photonID!="medium" && photonID!="loose" && photonID!="clusterOK") {
    std::cout << "Photon ID: '" << photonID << "' currently not implemented. Exiting." << std::endl;
    exit(11);
  }

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
    addInput( "PhotonJet_Summer1036X_Pt300to500" );
    addInput( "PhotonJet_Summer1036X_Pt500toInf" );

  } else if( dataset=="QCD_Spring10" ) {

    //addInput( "QCD_Spring10_Pt5to15" );
    addInput( "QCD_Spring10_Pt15to20" );
    addInput( "QCD_Spring10_Pt20to30" );
    addInput( "QCD_Spring10_Pt30to50" );
    addInput( "QCD_Spring10_Pt50to80" );
    addInput( "QCD_Spring10_Pt80" );
    addInput( "QCD_Spring10_Pt170" );
    addInput( "QCD_Spring10_Pt300" );
  //addInput( "QCD_Spring10_Pt470" );
  //addInput( "QCD_Spring10_Pt800" );

  } else if( dataset=="PhotonJet_Summer1036X_QCD_Spring10" ) { //this option is the one to be used in case of extrapolation

    addInput( "PhotonJet_Summer1036X_Pt5to15" );
    addInput( "PhotonJet_Summer1036X_Pt15to20" );
    addInput( "PhotonJet_Summer1036X_Pt20to30" );
    addInput( "PhotonJet_Summer1036X_Pt30to50" );
    addInput( "PhotonJet_Summer1036X_Pt50to80" );
    addInput( "PhotonJet_Summer1036X_Pt80to120" );
    addInput( "PhotonJet_Summer1036X_Pt120to170" );
    addInput( "PhotonJet_Summer1036X_Pt170to300" );
    addInput( "QCD_Spring10_Pt15to20" );
    addInput( "QCD_Spring10_Pt20to30" );
    addInput( "QCD_Spring10_Pt30to50" );
    addInput( "QCD_Spring10_Pt50to80" );
    addInput( "QCD_Spring10_Pt80" );
    addInput( "QCD_Spring10_Pt170" );
    addInput( "QCD_Spring10_Pt300" );

  } else if( dataset=="G_TuneZ2_7TeV_pythia6" ) {

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



  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Float_t ptMax = 100.;
  int nBins_photBinning = ptPhot_binning.size()-1;
  Double_t ptPhot_binning_array[100];
  for( unsigned i=0; i<nBins_photBinning; ++i)
    ptPhot_binning_array[i] = ptPhot_binning[i];

  TH1D* h1_ptJetReco = new TH1D("ptJetReco", "", 10, 0., ptMax);
  h1_ptJetReco->Sumw2();
  TH1D* h1_pt2ndJetReco = new TH1D("pt2ndJetReco", "", 10, 5., ptMax);
  h1_pt2ndJetReco->Sumw2();
  TH1D* h1_etaPhot = new TH1D("etaPhot", "", 15, -1.3, 1.3);
  h1_etaPhot->Sumw2();
  TH1D* h1_phiPhot = new TH1D("phiPhot", "", 15, -3.1416, 3.1416);
  h1_phiPhot->Sumw2();
  TH1D* h1_deltaPhi_Nm1 = new TH1D("deltaPhi_Nm1", "", 15, 3.1416/2., 3.1416);
  h1_deltaPhi_Nm1->Sumw2();
  TH1D* h1_ptSecondJetRel_Nm1 = new TH1D("ptSecondJetRel_Nm1", "", 15, 0., 1.5);
  h1_ptSecondJetRel_Nm1->Sumw2();

  TH1D* h1_deltaPhi_2ndJet = new TH1D("deltaPhi_2ndJet", "", 15, 0., 3.1416 );
  h1_deltaPhi_2ndJet->Sumw2();
  TH1D* h1_met = new TH1D("met", "", 50, 0., 300.);
  h1_met->Sumw2();
  TH1D* h1_met_noSmaj = new TH1D("met_noSmaj", "", 50, 0., 300.);
  h1_met_noSmaj->Sumw2();

  TH1D* h1_deltaPhi_phot_met = new TH1D("deltaPhi_phot_met", "", 15, -3.1416, 3.1416);
  h1_deltaPhi_phot_met->Sumw2();

  TH1D* h1_ptPhot = new TH1D("ptPhot", "", nBins_photBinning-1, ptPhot_binning_array);
  h1_ptPhot->Sumw2();
  TH1D* h1_ptPhot_noSmaj = new TH1D("ptPhot_noSmaj", "", nBins_photBinning-1, ptPhot_binning_array);
  h1_ptPhot_noSmaj->Sumw2();

  TH1D* h1_ptPhot_passedID = new TH1D("ptPhot_passedID", "", nBins_photBinning-1, ptPhot_binning_array);
  h1_ptPhot_passedID->Sumw2();
  TH1D* h1_etaPhot_passedID = new TH1D("etaPhot_passedID", "", 15, -1.3, 1.3);
  h1_etaPhot_passedID->Sumw2();
  TH1D* h1_phiPhot_passedID = new TH1D("phiPhot_passedID", "", 15, -3.1416, 3.1416);
  h1_phiPhot_passedID->Sumw2();

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

//TH1D* h1_hcalIsoPhotReco_clusterOK = new TH1D("hcalIsoPhotReco_clusterOK", "", 10, 0., 0.5);
//h1_hcalIsoPhotReco_clusterOK->Sumw2();
//TH1D* h1_hcalIsoEnergyPhotReco_clusterOK = new TH1D("hcalIsoEnergyPhotReco_clusterOK", "", 10, 0., 15.);
//h1_hcalIsoEnergyPhotReco_clusterOK->Sumw2();
//TH1D* h1_ecalIsoPhotReco_clusterOK = new TH1D("ecalIsoPhotReco_clusterOK", "", 20, 0., 1.);
//h1_ecalIsoPhotReco_clusterOK->Sumw2();
//TH1D* h1_ecalIsoEnergyPhotReco_clusterOK = new TH1D("ecalIsoEnergyPhotReco_clusterOK", "", 20, 0., 30.);
//h1_ecalIsoEnergyPhotReco_clusterOK->Sumw2();
//TH1D* h1_ptTrkIsoPhotReco_clusterOK = new TH1D("ptTrkIsoPhotReco_clusterOK", "", 20, 0., 1.);
//h1_ptTrkIsoPhotReco_clusterOK->Sumw2();
//TH1D* h1_nTrkIsoPhotReco_clusterOK = new TH1D("nTrkIsoPhotReco_clusterOK", "", 11, -0.5, 10.5);
//h1_nTrkIsoPhotReco_clusterOK->Sumw2();

//TH1D* h1_clusterMajPhotReco_isolated = new TH1D("clusterMajPhotReco_isolated", "", 30, 0., 1.5);
//h1_clusterMajPhotReco_isolated->Sumw2();
//TH1D* h1_clusterMinPhotReco_isolated = new TH1D("clusterMinPhotReco_isolated", "", 20, 0., 0.6);
//h1_clusterMinPhotReco_isolated->Sumw2();

  TH1D* h1_clusterMajPhotReco= new TH1D("clusterMajPhotReco", "", 30, 0., 1.5);
  h1_clusterMajPhotReco->Sumw2();
  TH1D* h1_clusterMinPhotReco= new TH1D("clusterMinPhotReco", "", 20, 0., 0.6);
  h1_clusterMinPhotReco->Sumw2();

  TH1D* h1_deltaPhi_passedID = new TH1D("deltaPhi_passedID", "", 15, 3.1416/2., 3.1416);
  h1_deltaPhi_passedID->Sumw2();
  TH1D* h1_ptSecondJetRel_passedID = new TH1D("ptSecondJetRel_passedID", "", 15, 0., 1.5);
  h1_ptSecondJetRel_passedID->Sumw2();


  Double_t ptPhotBinning_array[200]; //ugly! no more than 200 pt bins supported
  for( unsigned i=0; i<ptPhot_binning.size(); ++i )
    ptPhotBinning_array[i] = ptPhot_binning[i];


//TProfile* hp_ptJetGenMean = new TProfile("ptJetGenMean", "", ptPhot_binning.size()-1, ptPhotBinning_array);
//TProfile* hp_ptJetGenMean_loose = new TProfile("ptJetGenMean_loose", "", ptPhot_binning.size()-1, ptPhotBinning_array);
//TProfile* hp_ptJetGenMean_clusterOK = new TProfile("ptJetGenMean_clusterOK", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  TProfile* hp_ptPhotMean = new TProfile("prof_ptPhotMean", "", ptPhot_binning.size()-1, ptPhotBinning_array);
//TH2D* hp_ptPhotMean_no2ndJet = new TH2D("ptPhotMean_no2ndJet", "", ptPhot_binning.size()-1, ptPhotBinning_array, (int)floor(ptPhot_binning[ptPhot_binning.size()-1]), 0., ptPhot_binning[ptPhot_binning.size()-1]);
//TProfile* hp_ptPhotMean_loose = new TProfile("ptPhotMean_loose", "", ptPhot_binning.size()-1, ptPhotBinning_array);
//TProfile* hp_ptPhotMean_clusterOK = new TProfile("ptPhotMean_clusterOK", "", ptPhot_binning.size()-1, ptPhotBinning_array);

  TH2D* h2_ptPhotMean_eta011 = new TH2D("ptPhotMean_eta011", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_eta011->Sumw2();
  TH2D* h2_ptPhotMean_eta1524 = new TH2D("ptPhotMean_eta1524", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_eta1524->Sumw2();
  TH2D* h2_ptPhotMean_eta243 = new TH2D("ptPhotMean_eta243", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_eta243->Sumw2();
  TH2D* h2_ptPhotMean_eta35 = new TH2D("ptPhotMean_eta35", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_eta35->Sumw2();

  TH2D* h2_ptJetGenMean_eta011 = new TH2D("ptJetGenMean_eta011", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptJetGenMean_eta011->Sumw2();
  TH2D* h2_ptJetGenMean_eta1524 = new TH2D("ptJetGenMean_eta1524", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptJetGenMean_eta1524->Sumw2();
  TH2D* h2_ptJetGenMean_eta243 = new TH2D("ptJetGenMean_eta243", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptJetGenMean_eta243->Sumw2();
  TH2D* h2_ptJetGenMean_eta35 = new TH2D("ptJetGenMean_eta35", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptJetGenMean_eta35->Sumw2();

  TH2D* h2_ptPhotMean_no2ndJet_eta011 = new TH2D("ptPhotMean_no2ndJet_eta011", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_no2ndJet_eta011->Sumw2();
  TH2D* h2_ptPhotMean_no2ndJet_eta1524 = new TH2D("ptPhotMean_no2ndJet_eta1524", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_no2ndJet_eta1524->Sumw2();
  TH2D* h2_ptPhotMean_no2ndJet_eta243 = new TH2D("ptPhotMean_no2ndJet_eta243", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_no2ndJet_eta243->Sumw2();
  TH2D* h2_ptPhotMean_no2ndJet_eta35 = new TH2D("ptPhotMean_no2ndJet_eta35", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
  h2_ptPhotMean_no2ndJet_eta35->Sumw2();

//TH2D* h2_ptJetGenMean_loose = new TH2D("ptJetGenMean_loose", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
//h2_ptJetGenMean_loose->Sumw2();
//TH2D* h2_ptJetGenMean_clusterOK = new TH2D("ptJetGenMean_clusterOK", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
//h2_ptJetGenMean_clusterOK->Sumw2();
//TH2D* h2_ptPhotMean_loose = new TH2D("ptPhotMean_loose", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
//h2_ptPhotMean_loose->Sumw2();
//TH2D* h2_ptPhotMean_clusterOK = new TH2D("ptPhotMean_clusterOK", "", ptPhot_binning.size()-1, ptPhotBinning_array, 3000, 0., 3000.);
//h2_ptPhotMean_clusterOK->Sumw2();

  std::vector<TH1F*> h1_response_eta011              = getResponseHistos("response_eta011");
  std::vector<TH1F*> h1_responseMPF_eta011           = getResponseHistos("responseMPF_eta011");
  std::vector<TH1F*> h1_response_L2L3_eta011         = getResponseHistos("response_L2L3_eta011");
  std::vector<TH1F*> h1_responseMPF_L2L3_eta011      = getResponseHistos("responseMPF_L2L3_eta011");
  std::vector<TH1F*> h1_responseGEN_eta011           = getResponseHistos("responseGEN_eta011");
  std::vector<TH1F*> h1_responseGEN_L2L3_eta011      = getResponseHistos("responseGEN_L2L3_eta011");

  std::vector<TH1F*> h1_response_eta1524              = getResponseHistos("response_eta1524");
  std::vector<TH1F*> h1_responseMPF_eta1524           = getResponseHistos("responseMPF_eta1524");
  std::vector<TH1F*> h1_response_L2L3_eta1524         = getResponseHistos("response_L2L3_eta1524");
  std::vector<TH1F*> h1_responseMPF_L2L3_eta1524      = getResponseHistos("responseMPF_L2L3_eta1524");
  std::vector<TH1F*> h1_responseGEN_eta1524           = getResponseHistos("responseGEN_eta1524");
  std::vector<TH1F*> h1_responseGEN_L2L3_eta1524      = getResponseHistos("responseGEN_L2L3_eta1524");

  std::vector<TH1F*> h1_response_eta243              = getResponseHistos("response_eta243");
  std::vector<TH1F*> h1_responseMPF_eta243           = getResponseHistos("responseMPF_eta243");
  std::vector<TH1F*> h1_response_L2L3_eta243         = getResponseHistos("response_L2L3_eta243");
  std::vector<TH1F*> h1_responseMPF_L2L3_eta243      = getResponseHistos("responseMPF_L2L3_eta243");
  std::vector<TH1F*> h1_responseGEN_eta243           = getResponseHistos("responseGEN_eta243");
  std::vector<TH1F*> h1_responseGEN_L2L3_eta243      = getResponseHistos("responseGEN_L2L3_eta243");

  std::vector<TH1F*> h1_response_eta35              = getResponseHistos("response_eta35");
  std::vector<TH1F*> h1_responseMPF_eta35           = getResponseHistos("responseMPF_eta35");
  std::vector<TH1F*> h1_response_L2L3_eta35         = getResponseHistos("response_L2L3_eta35");
  std::vector<TH1F*> h1_responseMPF_L2L3_eta35      = getResponseHistos("responseMPF_L2L3_eta35");
  std::vector<TH1F*> h1_responseGEN_eta35           = getResponseHistos("responseGEN_eta35");
  std::vector<TH1F*> h1_responseGEN_L2L3_eta35      = getResponseHistos("responseGEN_L2L3_eta35");

//std::vector<TH1F*> h1_response_genpart_eta011     = getResponseHistos("response_genpart_eta011");
//std::vector<TH1F*> h1_response_partgamma_eta011   = getResponseHistos("response_partgamma_eta011");
//std::vector<TH1F*> h1_response_photgamma_eta011   = getResponseHistos("response_photgamma_eta011");

//std::vector<TH1F*> h1_response_clusterOK    = getResponseHistos("response_clusterOK");
//std::vector<TH1F*> h1_responseGEN_clusterOK = getResponseHistos("responseGEN_clusterOK");
//std::vector<TH1F*> h1_responsePART_clusterOK= getResponseHistos("responsePART_clusterOK");

//std::vector<TH1F*> h1_response_loose        = getResponseHistos("response_loose");
//std::vector<TH1F*> h1_response_loose_L2L3   = getResponseHistos("response_loose_L2L3");
//std::vector<TH1F*> h1_responseMPF_loose     = getResponseHistos("responseMPF_loose");
//std::vector<TH1F*> h1_responseMPF_loose_L2L3= getResponseHistos("responseMPF_loose_L2L3");
//std::vector<TH1F*> h1_responseGEN_loose     = getResponseHistos("responseGEN_loose");
//std::vector<TH1F*> h1_responsePART_loose    = getResponseHistos("responsePART_loose");

  std::vector<TH2D*> corr_genPart_partGamma;
  for( unsigned i=0; i<(ptPhot_binning.size()-1); ++i ) {
    char histoName[200];
    sprintf( histoName, "corr_genPart_partGamma_ptPhot_%.0f_%.0f", ptPhot_binning[i], ptPhot_binning[i+1]);
    TH2D* newHisto = new TH2D(histoName, "", 100, 0., 2., 100, 0., 2.);
    newHisto->Sumw2();
    corr_genPart_partGamma.push_back(newHisto);
  }

  // and now allocate histograms for 2nd jet pt extrapolation:
  bool isL2L3 = true;
  int nPoints = 5;
  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_recoRel_eta011 = getExtrapHistoVector("RecoPhot_vs_RecoRel_eta011", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_recoRel_eta011  = getExtrapHistoVector("RecoGen_vs_RecoRel_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_recoRel_eta011  = getExtrapHistoVector("GenPhot_vs_RecoRel_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoRelMean_eta011 = getExtrapHistoVector("pt2ndJetRecoRelMean_eta011", "rel", nPoints, !isL2L3, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhotL2L3_vs_recoRel_eta011 = getExtrapHistoVector("RecoPhotL2L3_vs_RecoRel_eta011", "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_RecoGenL2L3_vs_recoRel_eta011  = getExtrapHistoVector("RecoGenL2L3_vs_RecoRel_eta011",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_GenPhotL2L3_vs_recoRel_eta011  = getExtrapHistoVector("GenPhotL2L3_vs_RecoRel_eta011",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_GenPartL2L3_vs_recoRel_eta011  = getExtrapHistoVector("GenPartL2L3_vs_RecoRel_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenGammaL2L3_vs_recoRel_eta011  = getExtrapHistoVector("GenGammaL2L3_vs_RecoRel_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_PartGammaL2L3_vs_recoRel_eta011  = getExtrapHistoVector("PartGammaL2L3_vs_RecoRel_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_PhotGammaL2L3_vs_recoRel_eta011  = getExtrapHistoVector("PhotGammaL2L3_vs_RecoRel_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetL2L3RecoRelMean_eta011 = getExtrapHistoVector("pt2ndJetL2L3RecoRelMean_eta011", "rel", nPoints, isL2L3, (bool)true);
  std::vector< std::vector< TH1D* > > h1_RecoPhotL2L3_vs_recoRelRaw_eta011 = getExtrapHistoVector("RecoPhotL2L3_vs_RecoRelRaw_eta011", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGenL2L3_vs_recoRelRaw_eta011  = getExtrapHistoVector("RecoGenL2L3_vs_RecoRelRaw_eta011",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhotL2L3_vs_recoRelRaw_eta011  = getExtrapHistoVector("GenPhotL2L3_vs_RecoRelRaw_eta011",  "rel", nPoints);

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_recoRel_eta1524 = getExtrapHistoVector("RecoPhot_vs_RecoRel_eta1524", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_recoRel_eta1524  = getExtrapHistoVector("RecoGen_vs_RecoRel_eta1524",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_recoRel_eta1524  = getExtrapHistoVector("GenPhot_vs_RecoRel_eta1524",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoRelMean_eta1524 = getExtrapHistoVector("pt2ndJetRecoRelMean_eta1524", "rel", nPoints, !isL2L3, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhotL2L3_vs_recoRel_eta1524 = getExtrapHistoVector("RecoPhotL2L3_vs_RecoRel_eta1524", "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_RecoGenL2L3_vs_recoRel_eta1524  = getExtrapHistoVector("RecoGenL2L3_vs_RecoRel_eta1524",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_GenPhotL2L3_vs_recoRel_eta1524  = getExtrapHistoVector("GenPhotL2L3_vs_RecoRel_eta1524",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetL2L3RecoRelMean_eta1524 = getExtrapHistoVector("pt2ndJetL2L3RecoRelMean_eta1524", "rel", nPoints, isL2L3, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_recoRel_eta243 = getExtrapHistoVector("RecoPhot_vs_RecoRel_eta243", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_recoRel_eta243  = getExtrapHistoVector("RecoGen_vs_RecoRel_eta243",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_recoRel_eta243  = getExtrapHistoVector("GenPhot_vs_RecoRel_eta243",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoRelMean_eta243 = getExtrapHistoVector("pt2ndJetRecoRelMean_eta243", "rel", nPoints, !isL2L3, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhotL2L3_vs_recoRel_eta243 = getExtrapHistoVector("RecoPhotL2L3_vs_RecoRel_eta243", "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_RecoGenL2L3_vs_recoRel_eta243  = getExtrapHistoVector("RecoGenL2L3_vs_RecoRel_eta243",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_GenPhotL2L3_vs_recoRel_eta243  = getExtrapHistoVector("GenPhotL2L3_vs_RecoRel_eta243",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetL2L3RecoRelMean_eta243 = getExtrapHistoVector("pt2ndJetL2L3RecoRelMean_eta243", "rel", nPoints, isL2L3, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_recoRel_eta35 = getExtrapHistoVector("RecoPhot_vs_RecoRel_eta35", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_recoRel_eta35  = getExtrapHistoVector("RecoGen_vs_RecoRel_eta35",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_recoRel_eta35  = getExtrapHistoVector("GenPhot_vs_RecoRel_eta35",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoRelMean_eta35 = getExtrapHistoVector("pt2ndJetRecoRelMean_eta35", "rel", nPoints, !isL2L3, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhotL2L3_vs_recoRel_eta35 = getExtrapHistoVector("RecoPhotL2L3_vs_RecoRel_eta35", "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_RecoGenL2L3_vs_recoRel_eta35  = getExtrapHistoVector("RecoGenL2L3_vs_RecoRel_eta35",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_GenPhotL2L3_vs_recoRel_eta35  = getExtrapHistoVector("GenPhotL2L3_vs_RecoRel_eta35",  "rel", nPoints, isL2L3);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetL2L3RecoRelMean_eta35 = getExtrapHistoVector("pt2ndJetL2L3RecoRelMean_eta35", "rel", nPoints, isL2L3, (bool)true);


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
//Bool_t passed_Photon40;
//tree->SetBranchAddress("passed_Photon40", &passed_Photon40);
//Bool_t passed_Photon50;
//tree->SetBranchAddress("passed_Photon50", &passed_Photon50);
//Bool_t passed_Photon60;
//tree->SetBranchAddress("passed_Photon60", &passed_Photon60);
//Bool_t passed_Photon70;
//tree->SetBranchAddress("passed_Photon70", &passed_Photon70);





  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;
    Float_t correctWeight = (useGenJets) ? eventWeight_genjets : eventWeight;


    if( ONEVTX_ && dataset!="QCD_Spring10" ) {
      if( nvertex>1 ) continue;
    }


    bool jetInBarrel = (fabs(etaJetReco)<1.1);

    if( ptPhotReco<ptPhot_binning[0] ) continue;
    if( fabs(etaPhotReco)>1.3 ) continue;
    if( clusterMinPhotReco<0.15 && !MCassoc_ && !useGenJets ) continue; //protection vs EB spikes


    // trigger selection, only on data, to solve bias in response:
    bool isMC = run<5;
    if( !isMC ) {

      if( ptPhotReco<22. ) {
        if( !passed_Photon15 && !passed_Photon10 ) continue;
      } else if( ptPhotReco<32. ) {
        if( !passed_Photon20 && !passed_Photon25 ) continue;
      } else if( ptPhotReco<47. ) {
        if( !passed_Photon30 /*&& !passed_Photon35*/ ) continue;
      }

      if( ONEVTX_2232_ ) if( ptPhotReco>22. && ptPhotReco<32. && nvertex>1 ) continue;
      if( ONEVTX_2247_ ) if( ptPhotReco>22. && ptPhotReco<47. && nvertex>1 ) continue;

    //// select only Run2010A for first bins to avoid trigger prescale related issues:
    //if( ptPhotReco < 47. && run>145000. ) continue;

    }

    //first find correct photon pt bin:
    int theBin = hp_ptPhotMean->FindBin( ptPhotReco );
    theBin -= 1; //because arrays start from 0

    int theBinGEN = hp_ptPhotMean->FindBin( ptJetGen );
    theBinGEN -= 1; //because arrays start from 0


    //leading jet and photon back2back in transverse plane
    bool back2back = true;
    Float_t deltaPhi_jet = fabs(fitTools::delta_phi(phiPhotReco, phiJetReco));
    Float_t pi = TMath::Pi();
    float deltaPhiThreshold = DELTAPHI_;
    if( fabs(deltaPhi_jet) < (pi - deltaPhiThreshold) ) back2back = false; //loose back to back for now


    Float_t deltaPhi_2ndJet = fabs(fitTools::delta_phi(phiPhotReco, phi2ndJetReco));

    
    // parton flavour selection:
    bool partFlavorOK = true;
    if( PARTTYPE_=="QUARK" ) {
      partFlavorOK = ( (fabs(pdgIdPart)==1) || (fabs(pdgIdPart)==2) || (fabs(pdgIdPart)==3) );
    } else if( PARTTYPE_=="GLUON" ) {
      partFlavorOK = ( (fabs(pdgIdPart)==9) || (fabs(pdgIdPart)==21) );
    } else if( PARTTYPE_!="" ) {
      std::cout << " --> WARNING!!! Parton type '" << PARTTYPE_ << "' not implemented yet. No selection will be applied." << std::endl;
    }

    if( !partFlavorOK ) continue;


    bool secondJetOK;
    if( NO2ndJETABS )
      secondJetOK = ( pt2ndJetReco < secondJetThreshold*ptPhotReco );
    else
      //secondJetOK = ( ptCorr2ndJetReco < secondJetThreshold*ptPhotReco || pt2ndJetReco < 5. );
      secondJetOK = ( pt2ndJetReco < secondJetThreshold*ptPhotReco || pt2ndJetReco < 5. );


    // do them by hand just to be sure:
    Bool_t isIsolated_hcal;
    Bool_t isIsolated_ecal;
    Bool_t isIsolated_ptTracks;
    Bool_t isIsolated_nTracks;
    Bool_t clusterMajOK;
    Bool_t clusterMinOK;

    if( photonID=="medium" ) {

      isIsolated_hcal = ( hcalIsoPhotReco<0.05 || hcalIsoPhotReco*ePhotReco<2.4 );
      isIsolated_ecal = ( ecalIsoPhotReco<0.05  || ecalIsoPhotReco*ePhotReco<3. );
      isIsolated_ptTracks = ( ptTrkIsoPhotReco<0.1 );
      isIsolated_nTracks = (nTrkIsoPhotReco < 3 );
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    } else if( photonID=="loose" ) {

      isIsolated_hcal = ( hcalIsoPhotReco<0.1 || hcalIsoPhotReco*ePhotReco<4. );
      isIsolated_ecal = ( ecalIsoPhotReco<0.1  || ecalIsoPhotReco*ePhotReco<4.5 );
      isIsolated_ptTracks = ( ptTrkIsoPhotReco<0.2 );
      isIsolated_nTracks = (nTrkIsoPhotReco < 5 );
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    } else if( photonID=="clusterOK" ) {

      isIsolated_hcal = true;
      isIsolated_ecal = true;
      isIsolated_ptTracks = true;
      isIsolated_nTracks = true;
      clusterMajOK = ( clusterMajPhotReco>0.15 && clusterMajPhotReco<0.35 );
      clusterMinOK = ( clusterMinPhotReco>0.15 && clusterMinPhotReco<0.3 );

    }

  

    //before selection fill N-1 isolation plots (no event topology for isolation variables):
    if(                    isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_hcalIsoPhotReco_Nm1->Fill( hcalIsoPhotReco, correctWeight);
    if(                    isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_hcalIsoEnergyPhotReco_Nm1->Fill( hcalIsoPhotReco*ePhotReco, correctWeight);
    if( isIsolated_hcal                     && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_ecalIsoPhotReco_Nm1->Fill( ecalIsoPhotReco, correctWeight);
    if( isIsolated_hcal                     && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_ecalIsoEnergyPhotReco_Nm1->Fill( ecalIsoPhotReco*ePhotReco, correctWeight);
    if( isIsolated_hcal && isIsolated_ecal                         && isIsolated_nTracks && clusterMajOK && clusterMinOK  ) h1_ptTrkIsoPhotReco_Nm1->Fill( ptTrkIsoPhotReco, correctWeight);
    if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks                       && clusterMajOK && clusterMinOK  ) h1_nTrkIsoPhotReco_Nm1->Fill( nTrkIsoPhotReco, correctWeight);
    //no cluster cuts on cluster N-1's:
    if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks                                  ) h1_clusterMajPhotReco_Nm1->Fill( clusterMajPhotReco, correctWeight);
    if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks                                  ) h1_clusterMinPhotReco_Nm1->Fill( clusterMinPhotReco, correctWeight);
    // yes topology for topology variables:
    if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK              && secondJetOK && jetInBarrel) h1_deltaPhi_Nm1->Fill( deltaPhi_jet, correctWeight);
    if( isIsolated_hcal && isIsolated_ecal  && isIsolated_ptTracks && isIsolated_nTracks && clusterMajOK && clusterMinOK && back2back                && jetInBarrel) h1_ptSecondJetRel_Nm1->Fill( pt2ndJetReco/ptPhotReco, correctWeight);


//    bool isIsolated_loose = (isIsolated_hcal_loose && isIsolated_ecal_loose && isIsolated_ptTracks_loose && isIsolated_nTracks_loose);
    bool isIsolated = (isIsolated_hcal && isIsolated_ecal && isIsolated_ptTracks && isIsolated_nTracks);
    bool clusterShapeOK = (clusterMajOK && clusterMinOK );

    if( MCassoc_ && matchedToMC ) {
      isIsolated = true;
      //isIsolated_medium = true;
      clusterShapeOK= true;
    }

  
    //////////////////////////////////////////////
    /////      EVENT SELECTION
    //////////////////////////////////////////////

    // compute mpf :
    Float_t phi_Phot_Met = fitTools::delta_phi( phiPhotReco, phiMet );
    Float_t mpfResponse = 1. + eMet*ptPhotReco*cos( phi_Phot_Met ) / (ptPhotReco*ptPhotReco);

    //bool passedMedium_no2ndJet = photonOK_medium && ( (back2back && jetInBarrel) || noJetSelection );
    bool photonOK = (isIsolated && clusterShapeOK) || MCassoc_ || useGenJets;
    bool passedID_no2ndJet = photonOK && ( back2back || noJetSelection );
    bool passedID_FULL     = passedID_no2ndJet && (secondJetOK || noJetSelection);
    bool passedID_noSmaj     = isIsolated && clusterMinOK && (secondJetOK || noJetSelection);

    if( passedID_noSmaj ) {
      h1_met_noSmaj->Fill( eMet, correctWeight );
      h1_ptPhot_noSmaj->Fill( ptPhotReco, correctWeight );
    }


    if( passedID_no2ndJet ) {

      h1_deltaPhi_passedID->Fill( deltaPhi_jet, correctWeight);
      h1_ptSecondJetRel_passedID->Fill( pt2ndJetReco/ptPhotReco, correctWeight);
      h1_phiPhot_passedID->Fill( phiPhotReco, correctWeight );
      h1_etaPhot_passedID->Fill( etaPhotReco, correctWeight );
      h1_ptPhot_passedID->Fill( ptPhotReco, correctWeight );


      // --------------------------------------
      // BEGIN  extrapolation to pt(2ndjet)->0:
      // --------------------------------------
      if( HEMISPHERE_ )
        if( fabs(deltaPhi_2ndJet)< (3.14159/2.) ) continue;
   
        
      Float_t minPerc = h1_pt2ndJetRecoRelMean_eta011[theBin][0]->GetXaxis()->GetXmin();
      Float_t percStep = h1_pt2ndJetRecoRelMean_eta011[theBin][0]->GetXaxis()->GetXmax()  - minPerc;
      Float_t minPercL2L3 = h1_pt2ndJetL2L3RecoRelMean_eta011[theBin][0]->GetXaxis()->GetXmin();
      Float_t percStepL2L3 = h1_pt2ndJetL2L3RecoRelMean_eta011[theBin][0]->GetXaxis()->GetXmax()  - minPercL2L3;
      Double_t pt2ndJetRecoRel = 100.*pt2ndJetReco/ptPhotReco; //in percentage
      int iRecoRel = (int)floor((pt2ndJetRecoRel-minPerc)/percStep);
      Double_t ptCorr2ndJetRecoRel = 100.*ptCorr2ndJetReco/ptPhotReco; //in percentage
      int iRecoRelL2L3 = (int)floor((ptCorr2ndJetRecoRel-minPercL2L3)/percStepL2L3);
      //int iRecoRelL2L3 = iRecoRel;
    
      Float_t r_RecoPhot = ptJetReco/ptPhotReco;
      Float_t r_RecoGen  = ptJetReco/ptJetGen;
      Float_t r_RecoL2L3Phot = ptCorrJetReco/ptPhotReco;
      Float_t r_RecoL2L3Gen  = ptCorrJetReco/ptJetGen;
      Float_t r_GenPhot  = ptJetGen/ptPhotReco;
      Float_t r_GenPart  = ptJetGen/ptPart;
      Float_t r_GenGamma  = ptJetGen/ptPhotGen;
      Float_t r_PartGamma  = ptPart/ptPhotGen;
      Float_t r_PhotGamma  = ptPhotReco/ptPhotGen;

      float etaThresh = 1.1;

      if( fabs(etaJetReco) < etaThresh ) {

        h2_ptPhotMean_no2ndJet_eta011->Fill( ptPhotReco, ptPhotReco, correctWeight );

        if( (iRecoRel>=0)&&(iRecoRel<5) ) {
          if(r_RecoPhot!=0.) h1_RecoPhot_vs_recoRel_eta011[theBin][iRecoRel]->Fill(r_RecoPhot, correctWeight);
          if(r_RecoGen!=0.)  h1_RecoGen_vs_recoRel_eta011[theBin][iRecoRel]->Fill(r_RecoGen, correctWeight);
          if(r_GenPhot!=0.)  h1_GenPhot_vs_recoRel_eta011[theBin][iRecoRel]->Fill(r_GenPhot, correctWeight);
          if(r_RecoL2L3Phot!=0.) h1_RecoPhotL2L3_vs_recoRelRaw_eta011[theBin][iRecoRel]->Fill(r_RecoL2L3Phot, correctWeight);
          if(r_RecoL2L3Gen!=0.)  h1_RecoGenL2L3_vs_recoRelRaw_eta011[theBin][iRecoRel]->Fill(r_RecoL2L3Gen, correctWeight);
          if(r_GenPhot!=0.)      h1_GenPhotL2L3_vs_recoRelRaw_eta011[theBin][iRecoRel]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetRecoRelMean_eta011[theBin][iRecoRel]->Fill(pt2ndJetRecoRel, correctWeight);
        } 

        if( (iRecoRelL2L3>=0)&&(iRecoRelL2L3<5) ) {
          if(r_RecoL2L3Phot!=0.) h1_RecoPhotL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Phot, correctWeight);
          if(r_RecoL2L3Gen!=0.)  h1_RecoGenL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Gen, correctWeight);
          if(r_GenPhot!=0.)      h1_GenPhotL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_GenPhot, correctWeight);
          if(r_GenPart!=0.)      h1_GenPartL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_GenPart, correctWeight);
          if(r_GenGamma!=0.)     h1_GenGammaL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_GenGamma, correctWeight);
          if(r_PartGamma!=0.)    h1_PartGammaL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_PartGamma, correctWeight);
          if(r_PhotGamma!=0.)    h1_PhotGammaL2L3_vs_recoRel_eta011[theBin][iRecoRelL2L3]->Fill(r_PhotGamma, correctWeight);
          h1_pt2ndJetL2L3RecoRelMean_eta011[theBin][iRecoRelL2L3]->Fill(ptCorr2ndJetRecoRel, correctWeight);
        } 

      } else if( fabs(etaJetReco)>1.5 && fabs(etaJetReco)<2. ) {

        h2_ptPhotMean_no2ndJet_eta1524->Fill( ptPhotReco, ptPhotReco, correctWeight );

        if( (iRecoRel>=0)&&(iRecoRel<5) ) {
          if(r_RecoPhot!=0.) h1_RecoPhot_vs_recoRel_eta1524[theBin][iRecoRel]->Fill(r_RecoPhot, correctWeight);
          if(r_RecoGen!=0.)  h1_RecoGen_vs_recoRel_eta1524[theBin][iRecoRel]->Fill(r_RecoGen, correctWeight);
          if(r_GenPhot!=0.)  h1_GenPhot_vs_recoRel_eta1524[theBin][iRecoRel]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetRecoRelMean_eta1524[theBin][iRecoRel]->Fill(pt2ndJetRecoRel, correctWeight);
        } 

        if( (iRecoRelL2L3>=0)&&(iRecoRelL2L3<5) ) {
          if(r_RecoL2L3Phot!=0.) h1_RecoPhotL2L3_vs_recoRel_eta1524[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Phot, correctWeight);
          if(r_RecoL2L3Gen!=0.)  h1_RecoGenL2L3_vs_recoRel_eta1524[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Gen, correctWeight);
          if(r_GenPhot!=0.)      h1_GenPhotL2L3_vs_recoRel_eta1524[theBin][iRecoRelL2L3]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetL2L3RecoRelMean_eta1524[theBin][iRecoRelL2L3]->Fill(ptCorr2ndJetRecoRel, correctWeight);
        } 

      } else if( fabs(etaJetReco)<3. ) {

        h2_ptPhotMean_no2ndJet_eta243->Fill( ptPhotReco, ptPhotReco, correctWeight );

        if( (iRecoRel>=0)&&(iRecoRel<5) ) {
          if(r_RecoPhot!=0.) h1_RecoPhot_vs_recoRel_eta243[theBin][iRecoRel]->Fill(r_RecoPhot, correctWeight);
          if(r_RecoGen!=0.)  h1_RecoGen_vs_recoRel_eta243[theBin][iRecoRel]->Fill(r_RecoGen, correctWeight);
          if(r_GenPhot!=0.)  h1_GenPhot_vs_recoRel_eta243[theBin][iRecoRel]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetRecoRelMean_eta243[theBin][iRecoRel]->Fill(pt2ndJetRecoRel, correctWeight);
        } 

        if( (iRecoRelL2L3>=0)&&(iRecoRelL2L3<5) ) {
          if(r_RecoL2L3Phot!=0.) h1_RecoPhotL2L3_vs_recoRel_eta243[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Phot, correctWeight);
          if(r_RecoL2L3Gen!=0.)  h1_RecoGenL2L3_vs_recoRel_eta243[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Gen, correctWeight);
          if(r_GenPhot!=0.)      h1_GenPhotL2L3_vs_recoRel_eta243[theBin][iRecoRelL2L3]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetL2L3RecoRelMean_eta243[theBin][iRecoRelL2L3]->Fill(ptCorr2ndJetRecoRel, correctWeight);
        } 

      } else if( fabs(etaJetReco)<5. ) {

        h2_ptPhotMean_no2ndJet_eta35->Fill( ptPhotReco, ptPhotReco, correctWeight );

        if( (iRecoRel>=0)&&(iRecoRel<5) ) {
          if(r_RecoPhot!=0.) h1_RecoPhot_vs_recoRel_eta35[theBin][iRecoRel]->Fill(r_RecoPhot, correctWeight);
          if(r_RecoGen!=0.)  h1_RecoGen_vs_recoRel_eta35[theBin][iRecoRel]->Fill(r_RecoGen, correctWeight);
          if(r_GenPhot!=0.)  h1_GenPhot_vs_recoRel_eta35[theBin][iRecoRel]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetRecoRelMean_eta35[theBin][iRecoRel]->Fill(pt2ndJetRecoRel, correctWeight);
        } 

        if( (iRecoRelL2L3>=0)&&(iRecoRelL2L3<5) ) {
          if(r_RecoL2L3Phot!=0.) h1_RecoPhotL2L3_vs_recoRel_eta35[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Phot, correctWeight);
          if(r_RecoL2L3Gen!=0.)  h1_RecoGenL2L3_vs_recoRel_eta35[theBin][iRecoRelL2L3]->Fill(r_RecoL2L3Gen, correctWeight);
          if(r_GenPhot!=0.)      h1_GenPhotL2L3_vs_recoRel_eta35[theBin][iRecoRelL2L3]->Fill(r_GenPhot, correctWeight);
          h1_pt2ndJetL2L3RecoRelMean_eta35[theBin][iRecoRelL2L3]->Fill(ptCorr2ndJetRecoRel, correctWeight);
        } 

      }

      // --------------------------------------
      // END  extrapolation to pt(2ndjet)->0:
      // --------------------------------------


      
      if( passedID_FULL ) {

        h1_ptPhot->Fill( ptPhotReco, correctWeight );

        h1_deltaPhi_2ndJet->Fill( deltaPhi_2ndJet, correctWeight );
        

        h1_phiPhot->Fill( phiPhotReco, correctWeight );
        h1_etaPhot->Fill( etaPhotReco, correctWeight );
        h1_ptPhot->Fill( ptPhotReco, correctWeight );
        h1_ptJetReco->Fill( ptJetReco, correctWeight );
        h1_pt2ndJetReco->Fill( pt2ndJetReco, correctWeight );

        h1_deltaPhi_phot_met->Fill( phi_Phot_Met, correctWeight );
        h1_met->Fill( eMet, correctWeight );

        if( fabs(etaJetReco)<etaThresh ) {
          h2_ptPhotMean_eta011->Fill( ptPhotReco, ptPhotReco, correctWeight );
          h1_response_eta011[theBin]->Fill( r_RecoPhot, correctWeight );
          h1_response_L2L3_eta011[theBin]->Fill( r_RecoL2L3Phot, correctWeight );
          h1_responseMPF_eta011[theBin]->Fill( mpfResponse, correctWeight );
          //h1_responseMPF_L2L3_eta011[theBin]->Fill( mpfResponseL2L3, correctWeight );
          if( ptJetGen>ptPhot_binning[0] ) {
            h1_responseGEN_eta011[theBinGEN]->Fill( r_RecoGen, correctWeight );
            h1_responseGEN_L2L3_eta011[theBinGEN]->Fill( r_RecoL2L3Gen, correctWeight );
            h2_ptJetGenMean_eta011->Fill( ptJetGen, ptJetGen, correctWeight );
          corr_genPart_partGamma[theBinGEN]->Fill( r_PartGamma, r_GenPart, correctWeight );
          //h1_response_genpart_eta011[theBinGEN]->Fill( ptJetGen/ptPart, correctWeight );
          //h1_response_partgamma_eta011[theBinGEN]->Fill( ptPart/ptPhotGen, correctWeight );
          //h1_response_photgamma_eta011[theBinGEN]->Fill( ptPhotReco/ptPhotGen, correctWeight );
          }
        } else if( fabs(etaJetReco)>1.5 && fabs(etaJetReco)<2.4 ) {
          h2_ptPhotMean_eta1524->Fill( ptPhotReco, ptPhotReco, correctWeight );
          h1_response_eta1524[theBin]->Fill( r_RecoPhot, correctWeight );
          h1_response_L2L3_eta1524[theBin]->Fill( r_RecoL2L3Phot, correctWeight );
          h1_responseMPF_eta1524[theBin]->Fill( mpfResponse, correctWeight );
          //h1_responseMPF_L2L3_eta1524[theBin]->Fill( mpfResponseL2L3, correctWeight );
          if( ptJetGen>ptPhot_binning[0] ) {
            h1_responseGEN_eta1524[theBinGEN]->Fill( r_RecoGen, correctWeight );
            h1_responseGEN_L2L3_eta1524[theBinGEN]->Fill( r_RecoL2L3Gen, correctWeight );
            h2_ptJetGenMean_eta1524->Fill( ptJetGen, ptJetGen, correctWeight );
          }
        } else if( fabs(etaJetReco)>2.4 && fabs(etaJetReco)<3. ) {
          h2_ptPhotMean_eta243->Fill( ptPhotReco, ptPhotReco, correctWeight );
          h1_response_eta243[theBin]->Fill( r_RecoPhot, correctWeight );
          h1_response_L2L3_eta243[theBin]->Fill( r_RecoL2L3Phot, correctWeight );
          h1_responseMPF_eta243[theBin]->Fill( mpfResponse, correctWeight );
          //h1_responseMPF_L2L3_eta243[theBin]->Fill( mpfResponseL2L3, correctWeight );
          if( ptJetGen>ptPhot_binning[0] ) {
            h1_responseGEN_eta243[theBinGEN]->Fill( r_RecoGen, correctWeight );
            h1_responseGEN_L2L3_eta243[theBinGEN]->Fill( r_RecoL2L3Gen, correctWeight );
            h2_ptJetGenMean_eta243->Fill( ptJetGen, ptJetGen, correctWeight );
          }
        } else if( fabs(etaJetReco)<5. ) {
          h2_ptPhotMean_eta35->Fill( ptPhotReco, ptPhotReco, correctWeight );
          h1_response_eta35[theBin]->Fill( r_RecoPhot, correctWeight );
          h1_response_L2L3_eta35[theBin]->Fill( r_RecoL2L3Phot, correctWeight );
          h1_responseMPF_eta35[theBin]->Fill( mpfResponse, correctWeight );
          //h1_responseMPF_L2L3_eta35[theBin]->Fill( mpfResponseL2L3, correctWeight );
          if( ptJetGen>ptPhot_binning[0] ) {
            h1_responseGEN_eta35[theBinGEN]->Fill( r_RecoGen, correctWeight );
            h1_responseGEN_L2L3_eta35[theBinGEN]->Fill( r_RecoL2L3Gen, correctWeight );
            h2_ptJetGenMean_eta35->Fill( ptJetGen, ptJetGen, correctWeight );
          }
        }
          
/*
      // --------------------------------------
      // BEGIN  extrapolation to pt(2ndjet)->0:
      // --------------------------------------
      h2_ptPhotMean_no2ndJet->Fill( ptPhotReco, ptPhotReco, correctWeight );
      Float_t minPerc = h1_pt2ndJetRecoRelMean[theBin][0]->GetXaxis()->GetXmin();
      Float_t percStep = h1_pt2ndJetRecoRelMean[theBin][0]->GetXaxis()->GetXmax()  - minPerc;
      Double_t pt2ndJetRecoRel = 100.*pt2ndJetReco/ptPhotReco; //in percentage
      int iRecoRel = (int)floor((pt2ndJetRecoRel-minPerc)/percStep);
    
      Float_t r_RecoPhot = ptJetReco/ptPhotReco;
      Float_t r_RecoGen  = ptJetReco/ptJetGen;
      Float_t r_GenPhot  = ptJetGen/ptPhotReco;

      if( (iRecoRel>=0)&&(iRecoRel<5) ) {
        if(r_RecoPhot!=0.) h1_RecoPhot_vs_recoRel[theBin][iRecoRel]->Fill(r_RecoPhot, correctWeight);
        if(r_RecoGen!=0.)  h1_RecoGen_vs_recoRel[theBin][iRecoRel]->Fill(r_RecoGen, correctWeight);
        if(r_GenPhot!=0.)  h1_GenPhot_vs_recoRel[theBin][iRecoRel]->Fill(r_GenPhot, correctWeight);
        h1_pt2ndJetRecoRelMean[theBin][iRecoRel]->Fill(pt2ndJetRecoRel, correctWeight);
      } 
      // --------------------------------------
      // END  extrapolation to pt(2ndjet)->0:
      // --------------------------------------


      
      if( passedMedium_FULL ) {

        h1_ptPhot_medium->Fill( ptPhotReco, correctWeight );

        h1_deltaPhi_2ndJet_medium->Fill( deltaPhi_2ndJet, correctWeight );
        
        if( ptJetGen>ptPhot_binning[0] ) {
          h1_responseGEN[theBinGEN]->Fill( ptJetReco/ptJetGen, correctWeight );
          h1_responsePART[theBinGEN]->Fill( ptJetGen/ptPart, correctWeight );
          h2_ptJetGenMean->Fill( ptJetGen, ptJetGen, correctWeight );
        }

        h1_phiPhot->Fill( phiPhotReco, correctWeight );
        h1_etaPhot->Fill( etaPhotReco, correctWeight );
        h1_ptPhot->Fill( ptPhotReco, correctWeight );
        h1_ptJetReco->Fill( ptJetReco, correctWeight );
        h1_pt2ndJetReco->Fill( pt2ndJetReco, correctWeight );
        h2_ptPhotMean->Fill( ptPhotReco, ptPhotReco, correctWeight );
        h1_response[theBin]->Fill( response, correctWeight );

        h1_deltaPhi_phot_met_medium->Fill( phi_Phot_Met, correctWeight );
        h1_met_medium->Fill( eMet, correctWeight );
        h1_responseMPF[theBin]->Fill( mpfResponse, correctWeight );
          
*/
      } //if second jet ok

    } 


/*
    //////////////////////////////////////////////
    /////      EVENT SELECTION: LOOSE ID
    //////////////////////////////////////////////

    bool passedLoose_no2ndJet = photonOK_loose && ( (back2back && jetInBarrel) || noJetSelection );
    bool passedLoose_FULL     = passedLoose_no2ndJet && (secondJetOK || noJetSelection);


    if( passedLoose_no2ndJet ) {

      Float_t correctWeight = (useGenJets) ? eventWeight_loose : eventWeight;
      
      
    
      if( passedLoose_FULL ) {

        if( ptJetGen>ptPhot_binning[0] ) {
          h1_responseGEN_loose[theBinGEN]->Fill( ptJetReco/ptJetGen, correctWeight );
          h1_responsePART_loose[theBinGEN]->Fill( ptJetGen/ptPart, correctWeight );
          h2_ptJetGenMean_loose->Fill( ptJetGen, ptJetGen, correctWeight );
        }

        h1_ptPhot_loose->Fill( ptPhotReco, correctWeight );
        h1_phiPhot_loose->Fill( phiPhotReco, correctWeight );
        h1_etaPhot_loose->Fill( etaPhotReco, correctWeight );

        h1_deltaPhi_2ndJet_loose->Fill( deltaPhi_2ndJet, correctWeight );

        h2_ptPhotMean_loose->Fill( ptPhotReco, ptPhotReco, correctWeight );
        h1_response_loose[theBin]->Fill( response, correctWeight );

        h1_deltaPhi_phot_met_loose->Fill( phi_Phot_Met, correctWeight );
        h1_met_loose->Fill( eMet, correctWeight );
        h1_responseMPF_loose[theBin]->Fill( mpfResponse, correctWeight );

      } //if second jet ok

    } //if loose

*/

  } //for entries



  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaPhotonJet_"+dataset;
  else {
   if(dataset!="") outfileName = "PhotonJet_"+dataset;
   else outfileName = "PhotonJet";
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
  if( ADD12_ ) outfileName = outfileName + "_ADD12";
  //if( !noJetSelection && secondJetThreshold!=0.5 ) {
  if( !noJetSelection ) {
    std::string R = ( NO2ndJETABS ) ? "R" : "";
    char outfileName_char[300];
    sprintf( outfileName_char, "%s_2ndJet%d%s", outfileName.c_str(), (int)(100.*secondJetThreshold), R.c_str());
    std::string outfileName_str( outfileName_char );
    outfileName = outfileName_str;
  }
  if( HEMISPHERE_ ) outfileName += "_HEMI";
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h1_totalLumi->Write();

  h1_ptPhot->Write();
  h1_phiPhot->Write();
  h1_etaPhot->Write();
  h1_ptJetReco->Write();
  h1_pt2ndJetReco->Write();

  h1_ptPhot_noSmaj->Write();

  h1_ptPhot_passedID->Write();
  h1_phiPhot_passedID->Write();
  h1_etaPhot_passedID->Write();
  
  h1_met->Write();
  h1_met_noSmaj->Write();

  h1_deltaPhi_phot_met->Write();

  h1_deltaPhi_Nm1->Write();
  h1_ptSecondJetRel_Nm1->Write();

  h1_deltaPhi_passedID->Write();
  h1_ptSecondJetRel_passedID->Write();

  h1_deltaPhi_2ndJet->Write();

  h1_hcalIsoPhotReco_Nm1->Write();
  h1_hcalIsoEnergyPhotReco_Nm1->Write();
  h1_ecalIsoPhotReco_Nm1->Write();
  h1_ecalIsoEnergyPhotReco_Nm1->Write();
  h1_ptTrkIsoPhotReco_Nm1->Write();
  h1_nTrkIsoPhotReco_Nm1->Write();
  h1_clusterMajPhotReco_Nm1->Write();
  h1_clusterMinPhotReco_Nm1->Write();

  h1_clusterMajPhotReco->Write();
  h1_clusterMinPhotReco->Write();

  h2_ptPhotMean_eta011->Write();
  h2_ptPhotMean_eta1524->Write();
  h2_ptPhotMean_eta243->Write();
  h2_ptPhotMean_eta35->Write();

  h2_ptJetGenMean_eta011->Write();
  h2_ptJetGenMean_eta1524->Write();
  h2_ptJetGenMean_eta243->Write();
  h2_ptJetGenMean_eta35->Write();

  h2_ptPhotMean_no2ndJet_eta011->Write();
  h2_ptPhotMean_no2ndJet_eta1524->Write();
  h2_ptPhotMean_no2ndJet_eta243->Write();
  h2_ptPhotMean_no2ndJet_eta35->Write();


  for( unsigned i=0; i<h1_response_eta011.size(); ++i ) {

    outFile->cd();

    h1_response_eta011[i]->Write();
    h1_response_L2L3_eta011[i]->Write();
    h1_responseGEN_eta011[i]->Write();
    h1_responseGEN_L2L3_eta011[i]->Write();
    h1_responseMPF_eta011[i]->Write();
    h1_responseMPF_L2L3_eta011[i]->Write();

    h1_response_eta1524[i]->Write();
    h1_response_L2L3_eta1524[i]->Write();
    h1_responseGEN_eta1524[i]->Write();
    h1_responseGEN_L2L3_eta1524[i]->Write();
    h1_responseMPF_eta1524[i]->Write();
    h1_responseMPF_L2L3_eta1524[i]->Write();

    h1_response_eta243[i]->Write();
    h1_response_L2L3_eta243[i]->Write();
    h1_responseGEN_eta243[i]->Write();
    h1_responseGEN_L2L3_eta243[i]->Write();
    h1_responseMPF_eta243[i]->Write();
    h1_responseMPF_L2L3_eta243[i]->Write();

    h1_response_eta35[i]->Write();
    h1_response_L2L3_eta35[i]->Write();
    h1_responseGEN_eta35[i]->Write();
    h1_responseGEN_L2L3_eta35[i]->Write();
    h1_responseMPF_eta35[i]->Write();
    h1_responseMPF_L2L3_eta35[i]->Write();

    corr_genPart_partGamma[i]->Write();

    TF1* gaussian = new TF1("gaussian", "gaus");
    Float_t nSigma = 2.5;

    char dirName[40];
    sprintf(dirName, "extrap_ptBin_%d_%d", (int)ptPhot_binning[i], (int)ptPhot_binning[i+1]);
    outFile->mkdir(dirName);
    outFile->cd(dirName);

    for( int ip=0; ip<nPoints; ++ip) {
      
      fitTools::fitProjection(h1_RecoPhot_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhot_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_RecoGen_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_GenPhot_vs_recoRel_eta011[i][ip]->Write();
      h1_pt2ndJetRecoRelMean_eta011[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhotL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhotL2L3_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGenL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_RecoGenL2L3_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhotL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_GenPhotL2L3_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_GenPartL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_GenPartL2L3_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_GenGammaL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_GenGammaL2L3_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_PartGammaL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_PartGammaL2L3_vs_recoRel_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_PhotGammaL2L3_vs_recoRel_eta011[i][ip], gaussian, nSigma, "QR");
      h1_PhotGammaL2L3_vs_recoRel_eta011[i][ip]->Write();
      h1_pt2ndJetL2L3RecoRelMean_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_RecoPhotL2L3_vs_recoRelRaw_eta011[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhotL2L3_vs_recoRelRaw_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGenL2L3_vs_recoRelRaw_eta011[i][ip], gaussian, nSigma, "QR");
      h1_RecoGenL2L3_vs_recoRelRaw_eta011[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhotL2L3_vs_recoRelRaw_eta011[i][ip], gaussian, nSigma, "QR");
      h1_GenPhotL2L3_vs_recoRelRaw_eta011[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhot_vs_recoRel_eta1524[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhot_vs_recoRel_eta1524[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_recoRel_eta1524[i][ip], gaussian, nSigma, "QR");
      h1_RecoGen_vs_recoRel_eta1524[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_recoRel_eta1524[i][ip], gaussian, nSigma, "QR");
      h1_GenPhot_vs_recoRel_eta1524[i][ip]->Write();
      h1_pt2ndJetRecoRelMean_eta1524[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhotL2L3_vs_recoRel_eta1524[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhotL2L3_vs_recoRel_eta1524[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGenL2L3_vs_recoRel_eta1524[i][ip], gaussian, nSigma, "QR");
      h1_RecoGenL2L3_vs_recoRel_eta1524[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhotL2L3_vs_recoRel_eta1524[i][ip], gaussian, nSigma, "QR");
      h1_GenPhotL2L3_vs_recoRel_eta1524[i][ip]->Write();
      h1_pt2ndJetL2L3RecoRelMean_eta1524[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhot_vs_recoRel_eta243[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhot_vs_recoRel_eta243[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_recoRel_eta243[i][ip], gaussian, nSigma, "QR");
      h1_RecoGen_vs_recoRel_eta243[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_recoRel_eta243[i][ip], gaussian, nSigma, "QR");
      h1_GenPhot_vs_recoRel_eta243[i][ip]->Write();
      h1_pt2ndJetRecoRelMean_eta243[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhotL2L3_vs_recoRel_eta243[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhotL2L3_vs_recoRel_eta243[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGenL2L3_vs_recoRel_eta243[i][ip], gaussian, nSigma, "QR");
      h1_RecoGenL2L3_vs_recoRel_eta243[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhotL2L3_vs_recoRel_eta243[i][ip], gaussian, nSigma, "QR");
      h1_GenPhotL2L3_vs_recoRel_eta243[i][ip]->Write();
      h1_pt2ndJetL2L3RecoRelMean_eta243[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhot_vs_recoRel_eta35[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhot_vs_recoRel_eta35[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_recoRel_eta35[i][ip], gaussian, nSigma, "QR");
      h1_RecoGen_vs_recoRel_eta35[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_recoRel_eta35[i][ip], gaussian, nSigma, "QR");
      h1_GenPhot_vs_recoRel_eta35[i][ip]->Write();
      h1_pt2ndJetRecoRelMean_eta35[i][ip]->Write();
     
      fitTools::fitProjection(h1_RecoPhotL2L3_vs_recoRel_eta35[i][ip], gaussian, nSigma, "QR");
      h1_RecoPhotL2L3_vs_recoRel_eta35[i][ip]->Write();
      fitTools::fitProjection(h1_RecoGenL2L3_vs_recoRel_eta35[i][ip], gaussian, nSigma, "QR");
      h1_RecoGenL2L3_vs_recoRel_eta35[i][ip]->Write();
      fitTools::fitProjection(h1_GenPhotL2L3_vs_recoRel_eta35[i][ip], gaussian, nSigma, "QR");
      h1_GenPhotL2L3_vs_recoRel_eta35[i][ip]->Write();
      h1_pt2ndJetL2L3RecoRelMean_eta35[i][ip]->Write();
     
    } //for npoints
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

  delete h1_totalLumi;
  h1_totalLumi = 0;
  delete h1_deltaPhi_2ndJet;
  h1_deltaPhi_2ndJet= 0;
  delete h1_ptPhot_passedID;
  h1_ptPhot_passedID = 0;
  delete h1_etaPhot_passedID;
  h1_etaPhot_passedID = 0;
  delete h1_phiPhot_passedID;
  h1_phiPhot_passedID = 0;
  delete h1_deltaPhi_passedID;
  h1_deltaPhi_passedID = 0;
  delete h1_ptSecondJetRel_passedID;
  h1_ptSecondJetRel_passedID = 0;

  delete h2_ptPhotMean_eta011;
  h2_ptPhotMean_eta011 = 0;
  delete h2_ptPhotMean_eta1524;
  h2_ptPhotMean_eta1524 = 0;
  delete h2_ptPhotMean_eta243;
  h2_ptPhotMean_eta243 = 0;
  delete h2_ptPhotMean_eta35;
  h2_ptPhotMean_eta35 = 0;

  delete h2_ptJetGenMean_eta011;
  h2_ptJetGenMean_eta011 = 0;
  delete h2_ptJetGenMean_eta1524;
  h2_ptJetGenMean_eta1524 = 0;
  delete h2_ptJetGenMean_eta243;
  h2_ptJetGenMean_eta243 = 0;
  delete h2_ptJetGenMean_eta35;
  h2_ptJetGenMean_eta35 = 0;

  delete h2_ptPhotMean_no2ndJet_eta011;
  h2_ptPhotMean_no2ndJet_eta011 = 0;
  delete h2_ptPhotMean_no2ndJet_eta1524;
  h2_ptPhotMean_no2ndJet_eta1524 = 0;
  delete h2_ptPhotMean_no2ndJet_eta243;
  h2_ptPhotMean_no2ndJet_eta243 = 0;
  delete h2_ptPhotMean_no2ndJet_eta35;
  h2_ptPhotMean_no2ndJet_eta35 = 0;

  for( unsigned i=0; i< h1_response_eta011.size(); ++i) {

    delete h1_response_eta011[i];
    h1_response_eta011[i]=0;
    delete h1_response_L2L3_eta011[i];
    h1_response_L2L3_eta011[i]=0;
    delete h1_responseGEN_eta011[i];
    h1_responseGEN_eta011[i]=0;
    delete h1_responseGEN_L2L3_eta011[i];
    h1_responseGEN_L2L3_eta011[i]=0;
    delete h1_responseMPF_eta011[i];
    h1_responseMPF_eta011[i]=0;
    delete h1_responseMPF_L2L3_eta011[i];
    h1_responseMPF_L2L3_eta011[i]=0;

    delete h1_response_eta1524[i];
    h1_response_eta1524[i]=0;
    delete h1_response_L2L3_eta1524[i];
    h1_response_L2L3_eta1524[i]=0;
    delete h1_responseGEN_eta1524[i];
    h1_responseGEN_eta1524[i]=0;
    delete h1_responseGEN_L2L3_eta1524[i];
    h1_responseGEN_L2L3_eta1524[i]=0;
    delete h1_responseMPF_eta1524[i];
    h1_responseMPF_eta1524[i]=0;
    delete h1_responseMPF_L2L3_eta1524[i];
    h1_responseMPF_L2L3_eta1524[i]=0;

    delete h1_response_eta243[i];
    h1_response_eta243[i]=0;
    delete h1_response_L2L3_eta243[i];
    h1_response_L2L3_eta243[i]=0;
    delete h1_responseGEN_eta243[i];
    h1_responseGEN_eta243[i]=0;
    delete h1_responseGEN_L2L3_eta243[i];
    h1_responseGEN_L2L3_eta243[i]=0;
    delete h1_responseMPF_eta243[i];
    h1_responseMPF_eta243[i]=0;
    delete h1_responseMPF_L2L3_eta243[i];
    h1_responseMPF_L2L3_eta243[i]=0;

    delete h1_response_eta35[i];
    h1_response_eta35[i]=0;
    delete h1_response_L2L3_eta35[i];
    h1_response_L2L3_eta35[i]=0;
    delete h1_responseGEN_eta35[i];
    h1_responseGEN_eta35[i]=0;
    delete h1_responseGEN_L2L3_eta35[i];
    h1_responseGEN_L2L3_eta35[i]=0;
    delete h1_responseMPF_eta35[i];
    h1_responseMPF_eta35[i]=0;
    delete h1_responseMPF_L2L3_eta35[i];
    h1_responseMPF_L2L3_eta35[i]=0;

  }


  deleteExtrapHistoVector(h1_RecoPhot_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_RecoGen_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_GenPhot_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetRecoRelMean_eta011, nPoints);

  deleteExtrapHistoVector(h1_RecoPhotL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_RecoGenL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_GenPhotL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_GenPartL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_GenGammaL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_PartGammaL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_PhotGammaL2L3_vs_recoRel_eta011, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetL2L3RecoRelMean_eta011, nPoints);
  deleteExtrapHistoVector(h1_RecoPhotL2L3_vs_recoRelRaw_eta011, nPoints);
  deleteExtrapHistoVector(h1_RecoGenL2L3_vs_recoRelRaw_eta011, nPoints);
  deleteExtrapHistoVector(h1_GenPhotL2L3_vs_recoRelRaw_eta011, nPoints);

  deleteExtrapHistoVector(h1_RecoPhot_vs_recoRel_eta1524, nPoints);
  deleteExtrapHistoVector(h1_RecoGen_vs_recoRel_eta1524, nPoints);
  deleteExtrapHistoVector(h1_GenPhot_vs_recoRel_eta1524, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetRecoRelMean_eta1524, nPoints);

  deleteExtrapHistoVector(h1_RecoPhotL2L3_vs_recoRel_eta1524, nPoints);
  deleteExtrapHistoVector(h1_RecoGenL2L3_vs_recoRel_eta1524, nPoints);
  deleteExtrapHistoVector(h1_GenPhotL2L3_vs_recoRel_eta1524, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetL2L3RecoRelMean_eta1524, nPoints);

  deleteExtrapHistoVector(h1_RecoPhot_vs_recoRel_eta243, nPoints);
  deleteExtrapHistoVector(h1_RecoGen_vs_recoRel_eta243, nPoints);
  deleteExtrapHistoVector(h1_GenPhot_vs_recoRel_eta243, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetRecoRelMean_eta243, nPoints);

  deleteExtrapHistoVector(h1_RecoPhotL2L3_vs_recoRel_eta243, nPoints);
  deleteExtrapHistoVector(h1_RecoGenL2L3_vs_recoRel_eta243, nPoints);
  deleteExtrapHistoVector(h1_GenPhotL2L3_vs_recoRel_eta243, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetL2L3RecoRelMean_eta243, nPoints);

  deleteExtrapHistoVector(h1_RecoPhot_vs_recoRel_eta35, nPoints);
  deleteExtrapHistoVector(h1_RecoGen_vs_recoRel_eta35, nPoints);
  deleteExtrapHistoVector(h1_GenPhot_vs_recoRel_eta35, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetRecoRelMean_eta35, nPoints);

  deleteExtrapHistoVector(h1_RecoPhotL2L3_vs_recoRel_eta35, nPoints);
  deleteExtrapHistoVector(h1_RecoGenL2L3_vs_recoRel_eta35, nPoints);
  deleteExtrapHistoVector(h1_GenPhotL2L3_vs_recoRel_eta35, nPoints);
  deleteExtrapHistoVector(h1_pt2ndJetL2L3RecoRelMean_eta35, nPoints);

  delete tree;
  tree = 0;

  totalLumi = 0.;

}


void addInput( const std::string& dataset, bool genjets ) {


  // opening from filelist now deprecated (files have to be merged with merge_and_setWeights
  //std::string infileName = "files_PhotonJet_2ndLevel_" + dataset+"_" + ALGOTYPE_ +".txt";
  TH1F* h1_lumi;


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "PhotonJet_2ndLevelTreeW_" + dataset + suffix;
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


  
std::vector<TH1F*> getResponseHistos(const std::string& name) {

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();
  std::vector<TH1F*> returnVector;

  for( unsigned i=0; i<(ptPhot_binning.size()-1); ++i ) {
    char histoName[100];
    sprintf( histoName, "%s_ptPhot_%.0f_%.0f", name.c_str(), ptPhot_binning[i], ptPhot_binning[i+1]);
    int nbins = (BINNINGFINO_) ? 1000 : 140;
    float xmin = (BINNINGFINO_) ? -5. : 0.1;
    float xmax = (BINNINGFINO_) ? 10. : 2.;
    //TH1F* newHisto = new TH1F(histoName, "", 15, 0., 2.);
    TH1F* newHisto = new TH1F(histoName, "", nbins, xmin, xmax);
    newHisto->Sumw2();
    returnVector.push_back(newHisto);
  }

  return returnVector;

}


std::vector< std::vector< TH1D* > > getExtrapHistoVector(const std::string& name, const std::string& abs_rel, Int_t nPoints, bool isL2L3, bool is_pt) {


  std::vector< std::vector< TH1D* > > returnVector;

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();


  for( int i_ptBin=0; i_ptBin<(ptPhot_binning.size()-1); ++i_ptBin ) {

    std::vector< TH1D* > histoVector;
    returnVector.push_back( histoVector );

    Float_t minPt;
    Float_t ptStep;
    if( abs_rel=="abs" ) {
      minPt = 5.;
      ptStep = 5.;
    } else { //minpt and step are in percent:

      if( RECOTYPE_=="pf"||RECOTYPE_=="jpt" ) {
        if( ptPhot_binning[i_ptBin]<=80.) {
          minPt = 8.;
          ptStep = 3.;
        } else if( ptPhot_binning[i_ptBin]<=350.) {
          minPt = 5.;
          ptStep = 3.;
        } else {
          minPt = 2.;
          ptStep = 2.;
        }
        if( isL2L3 && ptPhot_binning[i_ptBin]<=40. ) minPt+=2.*ptStep;
      } else { //is calo
        if( ptPhot_binning[i_ptBin]<=80.) {
          minPt = 8.;
          ptStep = 1.5;
        } else if( ptPhot_binning[i_ptBin]<=350.) {
          minPt = 6.;
          ptStep = 2.;
        } else {
          minPt = 2.;
          ptStep = 2.;
        }
        if( isL2L3 && ptPhot_binning[i_ptBin]<=40. ) minPt+=2.*ptStep;
      } //if calo-pf
    }
 
    double iPt = minPt;
    Int_t nBins = 50;

    for( int i=0; i<nPoints; ++i) {

      Double_t xMin = (is_pt) ? iPt : 0.;
      Double_t xMax = (is_pt) ? (iPt+ptStep) : 2.;

      TH1D* tmp;

      std::string prefix = (is_pt) ? "" : "r_";
      char histName[100];
      sprintf(histName, "%s%s_%d_%d", prefix.c_str(), name.c_str(), i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, xMin, xMax); 
      tmp->Sumw2();
      returnVector[i_ptBin].push_back(tmp);

      iPt+=ptStep;

    } //for nPoints

  } //for ptBins


  return returnVector;

}
  

void deleteExtrapHistoVector(std::vector< std::vector< TH1D* > > histoVector, int nPoints) {

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  for( int i_ptBin=0; i_ptBin<(ptPhot_binning.size()-1); ++i_ptBin ) {

    for( int i=0; i<nPoints; ++i) {

      delete histoVector[i_ptBin][i];
      histoVector[i_ptBin][i] = 0;

    }
  }

}
