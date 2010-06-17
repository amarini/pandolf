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

#include "fitTools.C"


bool NOWEIGHTS_ = true;
bool DEBUG_ = false;

class Jet {
public:
  Float_t eventWeight;

  Float_t eReco;
  Float_t ptReco;
  Float_t etaReco;
  Float_t phiReco;
  Float_t eGen;
  Float_t ptGen;
  Float_t etaGen;
  Float_t phiGen;
  Float_t Rch;
  Float_t Rgamma;
  Float_t Rnh;
  Float_t Rmu;
  Float_t Rhfhad;
  Float_t Rhfem;
  Float_t Re;
  Float_t Ptch;
  Float_t Ptgamma;
  Float_t Ptnh;
  Float_t Ptmu;
  Float_t Pthf;
  Float_t Pte;
  Int_t Nch;
  Int_t Ngamma;
  Int_t Nnh;
  Int_t Nmu;
  Int_t Nhf;
  Int_t Ne;

  Float_t RchGen;
  Float_t RgammaGen;
  Float_t RnhGen;
  Float_t RmuGen;
  Float_t ReGen;
};



void fillFile(const std::string& outFileName, TTree* tree, std::string algoType, bool draw=false);



bool LONGSCALE_ = false;


void finalize(const std::string& dataset) {


  if( (dataset=="RelValQCDFlat15_3000")||(dataset=="PJ_long") )
    LONGSCALE_=true;

  std::string infile_dataset = dataset;
  if( dataset=="PJ_long" )
    infile_dataset = "PJ";

  TChain* chain = new TChain("jetTree");
  
  std::string infileName;

  if( dataset=="QCDptBins") {

    chain->Add("2ndLevelOutputFile_QCD_227_0_15.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_15_20.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_20_30.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_30_50.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_50_80.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_80_120.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_120_170.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_170_230.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_230_300.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_300_380.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_380_470.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_470_600.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_600_800.root/jetTree");
    chain->Add("2ndLevelOutputFile_QCD_227_800_1000.root/jetTree");

  } else if( dataset=="PhotonJet_Summer09") {

    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt15.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt15.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt30.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt30.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt80.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt80.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt170.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt170.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt300.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt300.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt470.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt470.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt800.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt800.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;
    chain->Add("2ndLevelOutputFile_PhotonJet_Summer09_Pt1400.root/jetTree");
    std::cout << "Added 2ndLevelOutputFile_PhotonJet_Summer09_Pt1400.root/jetTree. Chain has " << chain->GetEntries() << " entries." << std::endl;

  } else {

    if( dataset=="" )
      infileName = "2ndLevelOutputFile.root/jetTree";
    else
      infileName = "2ndLevelOutputFile_"+infile_dataset+".root/jetTree";

    std::cout << "Opening file: " << infileName << std::endl;
    chain->Add(infileName.c_str());

  }
  
  std::cout << "Done adding." << std::endl;

  bool draw=true;

  std::string fileName;
  if( DEBUG_ )
    fileName = "finalOutputFile_prova";
  else {
    if( dataset=="" ) fileName = "finalOutputFile";
    else fileName = "finalOutputFile_"+dataset;
  }

  if( NOWEIGHTS_ ) fileName = fileName + "_NOWEIGHTS";

  fileName = fileName + ".root";

  TFile* outfile = new TFile(fileName.c_str(), "RECREATE");
  outfile->Close();


//  std::cout << std::endl << std::endl << "**** starting calo anti-kt 5" << std::endl;
//  fillFile(fileName, chain, "caloakt5", draw);
//  std::cout << std::endl << std::endl << "**** starting calo SIS cone 5" << std::endl;
//  fillFile(fileName, chain, "caloSISCone5");
//  std::cout << std::endl << std::endl << "**** starting calo kt4" << std::endl;
//  fillFile(fileName, chain, "calokt4");
  std::cout << std::endl << std::endl << "**** starting PF anti-kt 5" << std::endl;
  fillFile(fileName, chain, "PFakt5", draw);

  outfile->Write();
  outfile->Close();
  delete outfile;
  outfile = 0;


}  //finalize




void fillFile(const std::string& outFileName, TTree* tree, std::string algoType, bool draw) {

  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  int nReco;
  std::string nReco_name = "nJetReco_"+algoType;
  tree->SetBranchAddress(nReco_name.c_str(), &nReco);
  Float_t ptReco[40];
  std::string ptReco_name = "ptJetReco_"+algoType;
  tree->SetBranchAddress(ptReco_name.c_str(), ptReco);
  Float_t etaReco[40];
  std::string etaReco_name = "etaJetReco_"+algoType;
  tree->SetBranchAddress(etaReco_name.c_str(), etaReco);
  Float_t phiReco[40];
  std::string phiReco_name = "phiJetReco_"+algoType;
  tree->SetBranchAddress(phiReco_name.c_str(), phiReco);

  //PF only variables:
  Float_t eReco[40];
  std::string eReco_name = "eJetReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eReco_name.c_str(), eReco);
  Int_t nTracksReco[40];
  std::string nTracksReco_name = "nTracksReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(nTracksReco_name.c_str(), nTracksReco);
  Int_t nPhotonsReco[40];
  std::string nPhotonsReco_name = "nPhotonsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(nPhotonsReco_name.c_str(), nPhotonsReco);
  Float_t ptTracksReco[40];
  std::string ptTracksReco_name = "ptTracksReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(ptTracksReco_name.c_str(), ptTracksReco);
  Float_t ptPhotonsReco[40];
  std::string ptPhotonsReco_name = "ptPhotonsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(ptPhotonsReco_name.c_str(), ptPhotonsReco);
  Float_t eTracksReco[40];
  std::string eTracksReco_name = "eTracksReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eTracksReco_name.c_str(), eTracksReco);
  Float_t ePhotonsReco[40];
  std::string ePhotonsReco_name = "ePhotonsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(ePhotonsReco_name.c_str(), ePhotonsReco);
  Float_t eNeutralHadronsReco[40];
  std::string eNeutralHadronsReco_name = "eNeutralHadronsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eNeutralHadronsReco_name.c_str(), eNeutralHadronsReco);
  Float_t ptNeutralHadronsReco[40];
  std::string ptNeutralHadronsReco_name = "ptNeutralHadronsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(ptNeutralHadronsReco_name.c_str(), ptNeutralHadronsReco);
  Int_t nNeutralHadronsReco[40];
  std::string nNeutralHadronsReco_name = "nNeutralHadronsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(nNeutralHadronsReco_name.c_str(), nNeutralHadronsReco);
  Float_t eMuonsReco[40];
  std::string eMuonsReco_name = "eMuonsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eMuonsReco_name.c_str(), eMuonsReco);
  Float_t eElectronsReco[40];
  std::string eElectronsReco_name = "eElectronsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eElectronsReco_name.c_str(), eElectronsReco);
  Float_t eHFhadronsReco[40];
  std::string eHFhadronsReco_name = "eHFHadronsReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eHFhadronsReco_name.c_str(), eHFhadronsReco);
  Float_t eHFEMReco[40];
  std::string eHFEMReco_name = "eHFEMReco_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eHFEMReco_name.c_str(), eHFEMReco);

  Float_t eGen[40];
  std::string eGen_name = "eJetGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eGen_name.c_str(), eGen);
  Float_t eTracksGen[40];
  std::string eTracksGen_name = "eTracksGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eTracksGen_name.c_str(), eTracksGen);
  Int_t nTracksGen[40];
  std::string nTracksGen_name = "nTracksGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(nTracksGen_name.c_str(), nTracksGen);
  Float_t ePhotonsGen[40];
  std::string ePhotonsGen_name = "ePhotonsGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(ePhotonsGen_name.c_str(), ePhotonsGen);
  Int_t nPhotonsGen[40];
  std::string nPhotonsGen_name = "nPhotonsGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(nPhotonsGen_name.c_str(), nPhotonsGen);
  Float_t eNeutralHadronsGen[40];
  std::string eNeutralHadronsGen_name = "eNeutralHadronsGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eNeutralHadronsGen_name.c_str(), eNeutralHadronsGen);
  Int_t nNeutralHadronsGen[40];
  std::string nNeutralHadronsGen_name = "nNeutralHadronsGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(nNeutralHadronsGen_name.c_str(), nNeutralHadronsGen);
  Float_t eMuonsGen[40];
  std::string eMuonsGen_name = "eMuonsGen_"+algoType;
  if( algoType=="PFakt5")
    tree->SetBranchAddress(eMuonsGen_name.c_str(), eMuonsGen);

  Float_t ptGen[40];
  std::string ptGen_name = "ptJetGen_"+algoType;
  tree->SetBranchAddress(ptGen_name.c_str(), ptGen);
  Float_t etaGen[40];
  std::string etaGen_name = "etaJetGen_"+algoType;
  tree->SetBranchAddress(etaGen_name.c_str(), etaGen);
  Float_t phiGen[40];
  std::string phiGen_name = "phiJetGen_"+algoType;
  tree->SetBranchAddress(phiGen_name.c_str(), phiGen);


  std::string histoName;

  histoName = "deltaR_barrelEndcap_"+algoType;
  TH1F* h1_deltaR_barrelEndcap = new TH1F(histoName.c_str(), "Matching #DeltaR Distribution (15<p_{T}^{GEN}<300 GeV/c, 0<|#eta_{GEN}|<3)", 100, 0., 1.);
  h1_deltaR_barrelEndcap->SetXTitle("#DeltaR");
  h1_deltaR_barrelEndcap->SetYTitle("a.u.");

  histoName = "deltaR_forward_"+algoType;
  TH1F* h1_deltaR_forward = new TH1F(histoName.c_str(), "Matching #DeltaR Distribution (15<p_{T}^{GEN}<300 GeV/c, 3<|#eta_{GEN}|<5)", 100, 0., 1.);
  h1_deltaR_forward->SetXTitle("#DeltaR");
  h1_deltaR_forward->SetYTitle("a.u.");

  histoName = "deltaR_veryForward_"+algoType;
  TH1F* h1_deltaR_veryForward = new TH1F(histoName.c_str(), "Matching #DeltaR Distribution (15<p_{T}^{GEN}<300 GeV/c, 4<|#eta_{GEN}|<5)", 100, 0., 1.);
  h1_deltaR_veryForward->SetXTitle("#DeltaR");
  h1_deltaR_veryForward->SetYTitle("a.u.");

  


  //for pT spectrums, variable binning to avoid huge error bars at high pT:
  bool from10GeV = true;

  int nBins = 24;
  if( from10GeV) nBins = 26;

  Double_t Lower[nBins];

  Int_t i=0;
  if( from10GeV) {
    Lower[i++] = 10.;
    Lower[i++] = 15.;
  }
  Lower[i++] = 20.;
  Lower[i++] = 25.;
  Lower[i++] = 30.;
  Lower[i++] = 36.;
  Lower[i++] = 42.;
  Lower[i++] = 48.;
  Lower[i++] = 55.;
  Lower[i++] = 62.;
  Lower[i++] = 70.;
  Lower[i++] = 78.;
  Lower[i++] = 87.;
  Lower[i++] = 97.;
  Lower[i++] = 108.;
  Lower[i++] = 121.;
  Lower[i++] = 135.;
  Lower[i++] = 152.;
  Lower[i++] = 172.;
  Lower[i++] = 197.;
  Lower[i++] = 231.;
  Lower[i++] = 283.;
  Lower[i++] = 397.;
  Lower[i++] = 600.;
  Lower[i++] = 1000.;
  Lower[i++] = 2000.;




  // histograms booking for all algorithms:
  histoName = "deltaEta_vs_pt_barrel_"+algoType;
  TH2D* h2_deltaEta_vs_pt_barrel = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_barrel->Sumw2();

  histoName = "deltaPhi_vs_pt_barrel_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_barrel = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_barrel->Sumw2();

  histoName = "deltaEta_vs_pt_endcap_"+algoType;
  TH2D* h2_deltaEta_vs_pt_endcap = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_endcap->Sumw2();

  histoName = "deltaPhi_vs_pt_endcap_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_endcap = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_endcap->Sumw2();

  histoName = "deltaEta_vs_pt_forward_"+algoType;
  TH2D* h2_deltaEta_vs_pt_forward = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaEta_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_forward->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_forward->Sumw2();

  histoName = "deltaPhi_vs_pt_forward_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_forward = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_forward->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_forward->Sumw2();


  histoName = "response_vs_pt_barrel_"+algoType;
  TH2D* h2_response_vs_pt_barrel = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_barrel->Sumw2();

  histoName = "response_vs_pt_centralBarrel_"+algoType;
  TH2D* h2_response_vs_pt_centralBarrel = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_centralBarrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_centralBarrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_centralBarrel->Sumw2();

  histoName = "response_vs_pt_endcap_"+algoType;
  TH2D* h2_response_vs_pt_endcap = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_endcap->Sumw2();

  histoName = "response_vs_pt_eta0_20_"+algoType;
  TH2D* h2_response_vs_pt_eta0_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta0_20->Sumw2();

  histoName = "response_vs_pt_eta0_3_"+algoType;
  TH2D* h2_response_vs_pt_eta0_3 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_3->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_3->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta0_3->Sumw2();

  histoName = "response_vs_pt_eta14_20_"+algoType;
  TH2D* h2_response_vs_pt_eta14_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta14_20->Sumw2();

  histoName = "response_vs_pt_eta26_32_"+algoType;
  TH2D* h2_response_vs_pt_eta26_32 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta26_32->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta26_32->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta26_32->Sumw2();

  histoName = "response_vs_pt_forward_"+algoType;
  TH2D* h2_response_vs_pt_forward = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_forward->Sumw2();

  histoName = "response_vs_eta_"+algoType;
  TH2D* h2_response_vs_eta = new TH2D(histoName.c_str(), "Jet Response vs. True #eta", 100, -5.5, 5.5, 120, 0., 2.);
  h2_response_vs_eta->SetXTitle("#eta");
  h2_response_vs_eta->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_eta->Sumw2();

  histoName = "ptGenMean_vs_pt_barrel_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_barrel = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel->Sumw2();

  histoName = "ptGenMean_vs_pt_centralBarrel_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_centralBarrel = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_centralBarrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_endcap = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap->Sumw2();


  histoName = "ptGenMean_vs_pt_eta0_20_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta0_20 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20->Sumw2();

  histoName = "ptGenMean_vs_pt_eta0_3_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta0_3 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_3->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_3->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_3->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta14_20 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20->Sumw2();

  histoName = "ptGenMean_vs_pt_eta26_32_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta26_32 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta26_32->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta26_32->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta26_32->Sumw2();



  //additional histograms for PF only:


  //--------------------- Rch --------------------------------------------

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_barrel->Sumw2();

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch0->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_barrel_Rch0->Sumw2();

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch10_30->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_barrel_Rch10_30->Sumw2();

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch10_30->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_barrel_Rch10_30->Sumw2();

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch50_70->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_barrel_Rch50_70->Sumw2();

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch50_70->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_barrel_Rch50_70->Sumw2();

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch90_100->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_barrel_Rch90_100->Sumw2();

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch90_100->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_barrel_Rch90_100->Sumw2();

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch0->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_endcap_Rch0->Sumw2();

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch0->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_endcap_Rch0->Sumw2();

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch10_30->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_endcap_Rch10_30->Sumw2();

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch10_30->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_endcap_Rch10_30->Sumw2();

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch50_70->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_endcap_Rch50_70->Sumw2();

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch50_70->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_endcap_Rch50_70->Sumw2();

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch90_100->SetYTitle("#Delta#eta");
  h2_deltaEta_vs_pt_endcap_Rch90_100->Sumw2();

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch90_100->SetYTitle("#Delta#phi");
  h2_deltaPhi_vs_pt_endcap_Rch90_100->Sumw2();


  histoName = "response_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_barrel_Rch0->Sumw2();

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch0_30";
  TH2D* h2_response_vs_pt_barrel_Rch0_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch0_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch0_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_barrel_Rch0_30->Sumw2();

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_barrel_Rch10_30->Sumw2();

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_barrel_Rch50_70->Sumw2();

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_barrel_Rch90_100->Sumw2();

  histoName = "response_vs_pt_centralBarrel_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_centralBarrel_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_centralBarrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_centralBarrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_centralBarrel_Rch0->Sumw2();

  histoName = "response_vs_pt_centralBarrel_"+algoType+"_Rch0_30";
  TH2D* h2_response_vs_pt_centralBarrel_Rch0_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_centralBarrel_Rch0_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_centralBarrel_Rch0_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_centralBarrel_Rch0_30->Sumw2();

  histoName = "response_vs_pt_centralBarrel_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_centralBarrel_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_centralBarrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_centralBarrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_centralBarrel_Rch10_30->Sumw2();

  histoName = "response_vs_pt_centralBarrel_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_centralBarrel_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_centralBarrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_centralBarrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_centralBarrel_Rch50_70->Sumw2();

  histoName = "response_vs_pt_centralBarrel_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_centralBarrel_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_centralBarrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_centralBarrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_centralBarrel_Rch90_100->Sumw2();

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_endcap_Rch0->Sumw2();

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch0_30";
  TH2D* h2_response_vs_pt_endcap_Rch0_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch0_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch0_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_endcap_Rch0_30->Sumw2();

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_endcap_Rch10_30->Sumw2();

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_endcap_Rch50_70->Sumw2();

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_endcap_Rch90_100->Sumw2();


  histoName = "response_vs_pt_eta0_20_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_eta0_20_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_20_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_20_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta0_20_Rch0->Sumw2();

  histoName = "response_vs_pt_eta0_20_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_eta0_20_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_20_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_20_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta0_20_Rch10_30->Sumw2();

  histoName = "response_vs_pt_eta0_20_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_eta0_20_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_20_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_20_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta0_20_Rch50_70->Sumw2();

  histoName = "response_vs_pt_eta0_20_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_eta0_20_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_20_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_20_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta0_20_Rch90_100->Sumw2();

  histoName = "response_vs_pt_eta14_20_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_eta14_20_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_20_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_20_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta14_20_Rch0->Sumw2();

  histoName = "response_vs_pt_eta14_20_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_eta14_20_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_20_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_20_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta14_20_Rch10_30->Sumw2();

  histoName = "response_vs_pt_eta14_20_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_eta14_20_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_20_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_20_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta14_20_Rch50_70->Sumw2();

  histoName = "response_vs_pt_eta14_20_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_eta14_20_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_20_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_20_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta14_20_Rch90_100->Sumw2();

  histoName = "response_vs_pt_eta26_32_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_eta26_32_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta26_32_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta26_32_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h2_response_vs_pt_eta26_32_Rch0->Sumw2();


  histoName = "ptGenMean_vs_pt_barrel_"+algoType+"Rch0";
  TH2D* h2_ptGenMean_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch0->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch0->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType+"Rch0";
  TH2D* h2_ptGenMean_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch0->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch0->Sumw2();

  histoName = "ptGenMean_vs_pt_eta0_20_"+algoType+"Rch0";
  TH2D* h2_ptGenMean_vs_pt_eta0_20_Rch0 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_20_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch0->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch0->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_"+algoType+"Rch0";
  TH2D* h2_ptGenMean_vs_pt_eta14_20_Rch0 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch0->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch0->Sumw2();

  histoName = "ptGenMean_vs_pt_eta26_32_"+algoType+"Rch0";
  TH2D* h2_ptGenMean_vs_pt_eta26_32_Rch0 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta26_32_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta26_32_Rch0->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta26_32_Rch0->Sumw2();


  histoName = "ptGenMean_vs_pt_barrel_"+algoType+"Rch0_30";
  TH2D* h2_ptGenMean_vs_pt_barrel_Rch0_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel_Rch0_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch0_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch0_30->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType+"Rch0_30";
  TH2D* h2_ptGenMean_vs_pt_endcap_Rch0_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap_Rch0_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch0_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch0_30->Sumw2();


  histoName = "ptGenMean_vs_pt_barrel_"+algoType+"Rch10_30";
  TH2D* h2_ptGenMean_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch10_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch10_30->Sumw2();

  histoName = "ptGenMean_vs_pt_centralBarrel_"+algoType+"Rch10_30";
  TH2D* h2_ptGenMean_vs_pt_centralBarrel_Rch10_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_centralBarrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel_Rch10_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel_Rch10_30->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType+"Rch10_30";
  TH2D* h2_ptGenMean_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch10_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch10_30->Sumw2();

  histoName = "ptGenMean_vs_pt_eta0_20_"+algoType+"Rch10_30";
  TH2D* h2_ptGenMean_vs_pt_eta0_20_Rch10_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_20_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch10_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch10_30->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_"+algoType+"Rch10_30";
  TH2D* h2_ptGenMean_vs_pt_eta14_20_Rch10_30 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch10_30->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch10_30->Sumw2();


  histoName = "ptGenMean_vs_pt_barrel_"+algoType+"Rch50_70";
  TH2D* h2_ptGenMean_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch50_70->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch50_70->Sumw2();

  histoName = "ptGenMean_vs_pt_centralBarrel_"+algoType+"Rch50_70";
  TH2D* h2_ptGenMean_vs_pt_centralBarrel_Rch50_70 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_centralBarrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel_Rch50_70->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel_Rch50_70->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType+"Rch50_70";
  TH2D* h2_ptGenMean_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch50_70->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch50_70->Sumw2();

  histoName = "ptGenMean_vs_pt_eta0_20_"+algoType+"Rch50_70";
  TH2D* h2_ptGenMean_vs_pt_eta0_20_Rch50_70 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_20_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch50_70->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch50_70->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_"+algoType+"Rch50_70";
  TH2D* h2_ptGenMean_vs_pt_eta14_20_Rch50_70 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch50_70->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch50_70->Sumw2();


  histoName = "ptGenMean_vs_pt_barrel_"+algoType+"Rch90_100";
  TH2D* h2_ptGenMean_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch90_100->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_Rch90_100->Sumw2();

  histoName = "ptGenMean_vs_pt_centralBarrel_"+algoType+"Rch90_100";
  TH2D* h2_ptGenMean_vs_pt_centralBarrel_Rch90_100 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_centralBarrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel_Rch90_100->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_centralBarrel_Rch90_100->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType+"Rch90_100";
  TH2D* h2_ptGenMean_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch90_100->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_Rch90_100->Sumw2();

  histoName = "ptGenMean_vs_pt_eta0_20_"+algoType+"Rch90_100";
  TH2D* h2_ptGenMean_vs_pt_eta0_20_Rch90_100 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_20_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch90_100->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_20_Rch90_100->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_"+algoType+"Rch90_100";
  TH2D* h2_ptGenMean_vs_pt_eta14_20_Rch90_100 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch90_100->SetYTitle("p_{T}^{GEN} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_Rch90_100->Sumw2();




  histoName = "Rch_"+algoType;
  TH1F* h1_Rch = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_Rch->SetXTitle("R_{ch}");
  h1_Rch->SetYTitle("entries");

  histoName = "Rch_barrel_"+algoType;
  TH1F* h1_Rch_barrel = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_Rch_barrel->SetXTitle("R_{ch}");
  h1_Rch_barrel->SetYTitle("entries");

  histoName = "Rch_barrel_50_80_"+algoType;
  TH1F* h1_Rch_barrel_50_80 = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_Rch_barrel_50_80->SetXTitle("R_{ch}");
  h1_Rch_barrel_50_80->SetYTitle("entries");

  histoName = "Rch_eta14_20_"+algoType;
  TH1F* h1_Rch_eta14_20 = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_Rch_eta14_20->SetXTitle("R_{ch}");
  h1_Rch_eta14_20->SetYTitle("entries");

  histoName = "eta_Rch0_"+algoType;
  TH1F* h1_eta_Rch0 = new TH1F(histoName.c_str(), "Eta distribution for R_{ch}=0", 100, -5.5, 5.5);
  h1_eta_Rch0->SetXTitle("#eta");
  h1_eta_Rch0->SetYTitle("entries");

  histoName = "eta_Rgamma0_"+algoType;
  TH1F* h1_eta_Rgamma0 = new TH1F(histoName.c_str(), "Eta distribution for R_{gamma}=0", 100, -5.5, 5.5);
  h1_eta_Rgamma0->SetXTitle("#eta");
  h1_eta_Rgamma0->SetYTitle("entries");

  histoName = "Rch_vs_eta_"+algoType;
  TProfile* hp_Rch_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rch_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rch_vs_eta->SetYTitle("entries");

  histoName = "Rch_vs_eta_80_90_"+algoType;
  TProfile* hp_Rch_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rch_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Rch_vs_eta_80_90->SetYTitle("entries");

  histoName = "Rch_vs_pt_"+algoType;
  TProfile* hp_Rch_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rch_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt->SetYTitle("entries");

  histoName = "Rch_vs_pt_barrel_"+algoType;
  TProfile* hp_Rch_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rch_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt_barrel->SetYTitle("entries");

  histoName = "N_vs_pt_eta0_20_"+algoType;
  TProfile* hp_N_vs_pt_eta0_20 = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_N_vs_pt_eta0_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_N_vs_pt_eta0_20->SetYTitle("entries");

  histoName = "N_vs_pt_barrel_"+algoType;
  TProfile* hp_N_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_N_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_N_vs_pt_barrel->SetYTitle("entries");

  histoName = "N_vs_pt_eta14_20_"+algoType;
  TProfile* hp_N_vs_pt_eta14_20 = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_N_vs_pt_eta14_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_N_vs_pt_eta14_20->SetYTitle("entries");

  histoName = "Nch_vs_pt_"+algoType;
  TProfile* hp_Nch_vs_pt = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Nch_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Nch_vs_pt->SetYTitle("entries");

  histoName = "Nch_vs_pt_barrel_"+algoType;
  TProfile* hp_Nch_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Nch_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Nch_vs_pt_barrel->SetYTitle("entries");

  histoName = "Ptch_vs_pt_"+algoType;
  TProfile* hp_Ptch_vs_pt = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 5000.);
  hp_Ptch_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ptch_vs_pt->SetYTitle("entries");

  histoName = "Ptch_vs_pt_barrel_"+algoType;
  TProfile* hp_Ptch_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{ch} vs. pt", nBins-1, Lower, 0., 5000.);
  hp_Ptch_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ptch_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rch_vs_pt_endcap_"+algoType;
  TProfile* hp_Rch_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rch_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt_endcap->SetYTitle("entries");


  histoName = "Rgamma_"+algoType;
  TH1F* h1_Rgamma = new TH1F(histoName.c_str(), "Jet Photon Energy Fraction", 100, 0., 1.1);
  h1_Rgamma->SetXTitle("R_{#gamma}");
  h1_Rgamma->SetYTitle("entries");

  histoName = "Rgamma_vs_pt_"+algoType;
  TProfile* hp_Rgamma_vs_pt = new TProfile(histoName.c_str(), "R_{#gamma} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rgamma_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rgamma_vs_pt->SetYTitle("entries");

  histoName = "Rgamma_vs_pt_barrel_"+algoType;
  TProfile* hp_Rgamma_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{#gamma} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rgamma_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rgamma_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rgamma_vs_pt_endcap_"+algoType;
  TProfile* hp_Rgamma_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{#gamma} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rgamma_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rgamma_vs_pt_endcap->SetYTitle("entries");

  histoName = "Ngamma_vs_pt_"+algoType;
  TProfile* hp_Ngamma_vs_pt = new TProfile(histoName.c_str(), "N_{#gamma} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Ngamma_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ngamma_vs_pt->SetYTitle("entries");

  histoName = "Ngamma_vs_pt_barrel_"+algoType;
  TProfile* hp_Ngamma_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{#gamma} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Ngamma_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ngamma_vs_pt_barrel->SetYTitle("entries");

  histoName = "Ptgamma_vs_pt_"+algoType;
  TProfile* hp_Ptgamma_vs_pt = new TProfile(histoName.c_str(), "N_{#gamma} vs. pt", nBins-1, Lower, 0., 5000.);
  hp_Ptgamma_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ptgamma_vs_pt->SetYTitle("entries");

  histoName = "Ptgamma_vs_pt_barrel_"+algoType;
  TProfile* hp_Ptgamma_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{#gamma} vs. pt", nBins-1, Lower, 0., 5000.);
  hp_Ptgamma_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ptgamma_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rgamma_vs_eta_"+algoType;
  TProfile* hp_Rgamma_vs_eta = new TProfile(histoName.c_str(), "R_{#gamma} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rgamma_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rgamma_vs_eta->SetYTitle("entries");

  histoName = "Rgamma_vs_eta_80_90_"+algoType;
  TProfile* hp_Rgamma_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{#gamma} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rgamma_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Rgamma_vs_eta_80_90->SetYTitle("entries");

  histoName = "Rgamma_vs_Rch_"+algoType;
  TProfile* hp_Rgamma_vs_Rch = new TProfile(histoName.c_str(), "R_{#gamma} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rgamma_vs_Rch->SetXTitle("R_{ch}");
  hp_Rgamma_vs_Rch->SetYTitle("R_{#gamma}");

  histoName = "Rgamma_vs_Rch_barrel_"+algoType;
  TProfile* hp_Rgamma_vs_Rch_barrel = new TProfile(histoName.c_str(), "R_{#gamma} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rgamma_vs_Rch_barrel->SetXTitle("R_{ch}");
  hp_Rgamma_vs_Rch_barrel->SetYTitle("R_{#gamma}");


  histoName = "Rnh_"+algoType;
  TH1F* h1_Rnh = new TH1F(histoName.c_str(), "Jet Neutral Hadron Energy Fraction", 100, 0., 1.1);
  h1_Rnh->SetXTitle("R_{#gamma}");
  h1_Rnh->SetYTitle("entries");

  histoName = "Rnh_vs_pt_"+algoType;
  TProfile* hp_Rnh_vs_pt = new TProfile(histoName.c_str(), "R_{nh} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rnh_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rnh_vs_pt->SetYTitle("entries");

  histoName = "Rnh_vs_pt_barrel_"+algoType;
  TProfile* hp_Rnh_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{nh} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rnh_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rnh_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rnh_vs_pt_endcap_"+algoType;
  TProfile* hp_Rnh_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{nh} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rnh_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rnh_vs_pt_endcap->SetYTitle("entries");

  histoName = "Nnh_vs_pt_"+algoType;
  TProfile* hp_Nnh_vs_pt = new TProfile(histoName.c_str(), "N_{nh} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Nnh_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Nnh_vs_pt->SetYTitle("entries");

  histoName = "Nnh_vs_pt_barrel_"+algoType;
  TProfile* hp_Nnh_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{nh} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Nnh_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Nnh_vs_pt_barrel->SetYTitle("entries");

  histoName = "Ptnh_vs_pt_"+algoType;
  TProfile* hp_Ptnh_vs_pt = new TProfile(histoName.c_str(), "N_{nh} vs. pt", nBins-1, Lower, 0., 5000.);
  hp_Ptnh_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ptnh_vs_pt->SetYTitle("entries");

  histoName = "Ptnh_vs_pt_barrel_"+algoType;
  TProfile* hp_Ptnh_vs_pt_barrel = new TProfile(histoName.c_str(), "N_{nh} vs. pt", nBins-1, Lower, 0., 5000.);
  hp_Ptnh_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Ptnh_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rnh_vs_eta_"+algoType;
  TProfile* hp_Rnh_vs_eta= new TProfile(histoName.c_str(), "R_{nh} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rnh_vs_eta->SetYTitle("entries");

  histoName = "Rnh_vs_eta_80_90_"+algoType;
  TProfile* hp_Rnh_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{nh} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Rnh_vs_eta_80_90->SetYTitle("entries");

  histoName = "Rnh_vs_Rch_"+algoType;
  TProfile* hp_Rnh_vs_Rch = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rnh_vs_Rch->SetXTitle("R_{ch}");
  hp_Rnh_vs_Rch->SetYTitle("R_{nh}");

  histoName = "Rnh_vs_Rch_barrel_"+algoType;
  TProfile* hp_Rnh_vs_Rch_barrel = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rnh_vs_Rch_barrel->SetXTitle("R_{ch}");
  hp_Rnh_vs_Rch_barrel->SetYTitle("R_{nh}");


 
  histoName = "Re_"+algoType;
  TH1F* h1_Re = new TH1F(histoName.c_str(), "Jet Electron Energy Fraction", 100, 0., 1.1);
  h1_Re->SetXTitle("R_{e}");
  h1_Re->SetYTitle("entries");

  histoName = "Re_vs_pt_"+algoType;
  TProfile* hp_Re_vs_pt = new TProfile(histoName.c_str(), "R_{e} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Re_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Re_vs_pt->SetYTitle("entries");

  histoName = "Re_vs_pt_barrel_"+algoType;
  TProfile* hp_Re_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{e} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Re_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Re_vs_pt_barrel->SetYTitle("entries");

  histoName = "Re_vs_pt_endcap_"+algoType;
  TProfile* hp_Re_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{e} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Re_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Re_vs_pt_endcap->SetYTitle("entries");

  histoName = "Re_vs_eta_"+algoType;
  TProfile* hp_Re_vs_eta= new TProfile(histoName.c_str(), "R_{e} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Re_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Re_vs_eta->SetYTitle("entries");

  histoName = "Re_vs_eta_80_90_"+algoType;
  TProfile* hp_Re_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{e} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Re_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Re_vs_eta_80_90->SetYTitle("entries");

  histoName = "Re_vs_Rch_"+algoType;
  TProfile* hp_Re_vs_Rch = new TProfile(histoName.c_str(), "R_{e} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Re_vs_Rch->SetXTitle("R_{ch}");
  hp_Re_vs_Rch->SetYTitle("R_{e}");

 
  histoName = "Rmu_"+algoType;
  TH1F* h1_Rmu = new TH1F(histoName.c_str(), "Jet Muon Energy Fraction", 100, 0., 1.1);
  h1_Rmu->SetXTitle("R_{#gamma}");
  h1_Rmu->SetYTitle("entries");

  histoName = "Rmu_vs_pt_"+algoType;
  TProfile* hp_Rmu_vs_pt = new TProfile(histoName.c_str(), "R_{#mu} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rmu_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rmu_vs_pt->SetYTitle("entries");

  histoName = "Rmu_vs_pt_barrel_"+algoType;
  TProfile* hp_Rmu_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{#mu} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rmu_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rmu_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rmu_vs_pt_endcap_"+algoType;
  TProfile* hp_Rmu_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{#mu} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rmu_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rmu_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rmu_vs_eta_"+algoType;
  TProfile* hp_Rmu_vs_eta= new TProfile(histoName.c_str(), "R_{#mu} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rmu_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rmu_vs_eta->SetYTitle("entries");

  histoName = "Rmu_vs_eta_80_90_"+algoType;
  TProfile* hp_Rmu_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{#mu} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rmu_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Rmu_vs_eta_80_90->SetYTitle("entries");

  histoName = "Rmu_vs_Rch_"+algoType;
  TProfile* hp_Rmu_vs_Rch = new TProfile(histoName.c_str(), "R_{#mu} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rmu_vs_Rch->SetXTitle("R_{ch}");
  hp_Rmu_vs_Rch->SetYTitle("R_{#mu}");


  histoName = "Rhfhad_"+algoType;
  TH1F* h1_Rhfhad = new TH1F(histoName.c_str(), "Jet HF Hadron Energy Fraction", 100, 0., 1.1);
  h1_Rhfhad->SetXTitle("R_{HF}");
  h1_Rhfhad->SetYTitle("entries");

  histoName = "Rhfhad_vs_pt_"+algoType;
  TProfile* hp_Rhfhad_vs_pt = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rhfhad_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhfhad_vs_pt->SetYTitle("entries");

  histoName = "Rhfhad_vs_pt_endcap_"+algoType;
  TProfile* hp_Rhfhad_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rhfhad_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhfhad_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rhfhad_vs_eta_"+algoType;
  TProfile* hp_Rhfhad_vs_eta= new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhfhad_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rhfhad_vs_eta->SetYTitle("entries");

  histoName = "Rhfhad_vs_eta_80_90_"+algoType;
  TProfile* hp_Rhfhad_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhfhad_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Rhfhad_vs_eta_80_90->SetYTitle("entries");

  histoName = "Rhfem_"+algoType;
  TH1F* h1_Rhfem = new TH1F(histoName.c_str(), "Jet HF EM Energy Fraction", 100, 0., 1.1);
  h1_Rhfem->SetXTitle("R_{HF}");
  h1_Rhfem->SetYTitle("entries");

  histoName = "Rhfem_vs_pt_"+algoType;
  TProfile* hp_Rhfem_vs_pt = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rhfem_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhfem_vs_pt->SetYTitle("entries");

  histoName = "Rhfem_vs_pt_endcap_"+algoType;
  TProfile* hp_Rhfem_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rhfem_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhfem_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rhfem_vs_eta_"+algoType;
  TProfile* hp_Rhfem_vs_eta= new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhfem_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rhfem_vs_eta->SetYTitle("entries");

  histoName = "Rhfem_vs_eta_80_90_"+algoType;
  TProfile* hp_Rhfem_vs_eta_80_90 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhfem_vs_eta_80_90->SetXTitle("#eta^{GEN}");
  hp_Rhfem_vs_eta_80_90->SetYTitle("entries");


  // ----------------  R_x GEN plots:
  histoName = "RchGen_"+algoType;
  TH1F* h1_RchGen = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_RchGen->SetXTitle("R_{ch}");
  h1_RchGen->SetYTitle("entries");

  histoName = "RchGen_barrel_"+algoType;
  TH1F* h1_RchGen_barrel = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_RchGen_barrel->SetXTitle("R_{ch}");
  h1_RchGen_barrel->SetYTitle("entries");

  histoName = "RchGen_barrel_50_80_"+algoType;
  TH1F* h1_RchGen_barrel_50_80 = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_RchGen_barrel_50_80->SetXTitle("R_{ch}");
  h1_RchGen_barrel_50_80->SetYTitle("entries");

  histoName = "RchGen_endcap_"+algoType;
  TH1F* h1_RchGen_endcap = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_RchGen_endcap->SetXTitle("R_{ch}");
  h1_RchGen_endcap->SetYTitle("entries");

  histoName = "RchGen_vs_eta_Rch0_"+algoType;
  TProfile* hp_RchGen_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RchGen_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_RchGen_vs_eta_Rch0->SetYTitle("entries");

  histoName = "RchGen_vs_eta_"+algoType;
  TProfile* hp_RchGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RchGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RchGen_vs_eta->SetYTitle("entries");

  histoName = "RchGen_vs_pt_"+algoType;
  TProfile* hp_RchGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RchGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt->SetYTitle("entries");

  histoName = "RchGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RchGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RchGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RchGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RchGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RchGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RgammaGen_"+algoType;
  TH1F* h1_RgammaGen = new TH1F(histoName.c_str(), "Jet Photon Energy Fraction", 100, 0., 1.1);
  h1_RgammaGen->SetXTitle("R_{#gamma}");
  h1_RgammaGen->SetYTitle("entries");

  histoName = "RgammaGen_vs_eta_Rch0_"+algoType;
  TProfile* hp_RgammaGen_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RgammaGen_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_RgammaGen_vs_eta_Rch0->SetYTitle("entries");

  histoName = "RgammaGen_vs_eta_"+algoType;
  TProfile* hp_RgammaGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RgammaGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RgammaGen_vs_eta->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_"+algoType;
  TProfile* hp_RgammaGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RgammaGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RgammaGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RgammaGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RgammaGen_vs_RchGen_"+algoType;
  TProfile* hp_RgammaGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{#gamma} vs. R_{ch}", 100, 0., 1.1, 0., 1.1);
  hp_RgammaGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RgammaGen_vs_RchGen->SetYTitle("R_{#gamma}");

  histoName = "RnhGen_"+algoType;
  TH1F* h1_RnhGen = new TH1F(histoName.c_str(), "Jet Neutral Hadron Energy Fraction", 100, 0., 1.1);
  h1_RnhGen->SetXTitle("R_{#gamma}");
  h1_RnhGen->SetYTitle("entries");

  histoName = "RnhGen_vs_eta_Rch0_"+algoType;
  TProfile* hp_RnhGen_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RnhGen_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_RnhGen_vs_eta_Rch0->SetYTitle("entries");

  histoName = "RnhGen_vs_eta_"+algoType;
  TProfile* hp_RnhGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RnhGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RnhGen_vs_eta->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_"+algoType;
  TProfile* hp_RnhGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RnhGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RnhGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RnhGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RnhGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RnhGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RnhGen_vs_RchGen_"+algoType;
  TProfile* hp_RnhGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1.1);
  hp_RnhGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RnhGen_vs_RchGen->SetYTitle("R_{nh}");

  histoName = "RmuGen_"+algoType;
  TH1F* h1_RmuGen = new TH1F(histoName.c_str(), "Jet Neutral Hadron Energy Fraction", 100, 0., 1.1);
  h1_RmuGen->SetXTitle("R_{#gamma}");
  h1_RmuGen->SetYTitle("entries");

  histoName = "RmuGen_vs_eta_"+algoType;
  TProfile* hp_RmuGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1.1);
  hp_RmuGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RmuGen_vs_eta->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_"+algoType;
  TProfile* hp_RmuGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RmuGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RmuGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RmuGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RmuGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1.1);
  hp_RmuGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RmuGen_vs_RchGen_"+algoType;
  TProfile* hp_RmuGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1.1);
  hp_RmuGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RmuGen_vs_RchGen->SetYTitle("R_{nh}");


  // E_x resolutions:
  histoName = "Ech_vs_pt_barrel_"+algoType;
  TH2D* h2_Ech_vs_pt_barrel = new TH2D(histoName.c_str(), "Ech reco over Ech gen",  nBins-1, Lower, 400, 0., 4.);
  h2_Ech_vs_pt_barrel->Sumw2();

  histoName = "Egammanh_vs_pt_barrel_"+algoType;
  TH2D* h2_Egammanh_vs_pt_barrel = new TH2D(histoName.c_str(), "Egammanh reco over Egammanh gen",  nBins-1, Lower, 400, 0., 4.);
  h2_Egammanh_vs_pt_barrel->Sumw2();

  histoName = "Ech_vs_pt_endcap_"+algoType;
  TH2D* h2_Ech_vs_pt_endcap = new TH2D(histoName.c_str(), "Ech reco over Ech gen",  nBins-1, Lower, 400, 0., 4.);
  h2_Ech_vs_pt_endcap->Sumw2();

  histoName = "Egammanh_vs_pt_endcap_"+algoType;
  TH2D* h2_Egammanh_vs_pt_endcap = new TH2D(histoName.c_str(), "Egammanh reco over Egammanh gen",  nBins-1, Lower, 400, 0., 4.);
  h2_Egammanh_vs_pt_endcap->Sumw2();


  std::vector< TH2D* > vec_h2_Egammanh_vs_Ech_barrel;
  for( int i=0; i<nBins; ++i) {
    char h2_name[70];
    sprintf(h2_name, "Egammanh_vs_Ech_barrel_bin%d", (i+1));
    TH2D* h2 = new TH2D(h2_name, "", 400, 0., 4., 400, 0., 4.);
    h2->Sumw2();
    vec_h2_Egammanh_vs_Ech_barrel.push_back(h2);
  }

  std::vector< TH2D* > vec_h2_Egammanh_vs_Ech_endcap;
  for( int i=0; i<nBins; ++i) {
    char h2_name[70];
    sprintf(h2_name, "Egammanh_vs_Ech_endcap_bin%d", (i+1));
    TH2D* h2 = new TH2D(h2_name, "", 400, 0., 4., 400, 0., 4.);
    h2->Sumw2();
    vec_h2_Egammanh_vs_Ech_endcap.push_back(h2);
  }






  int nEntries = tree->GetEntries();

  if( DEBUG_ ) nEntries = 100000.;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( NOWEIGHTS_ )
      eventWeight = 1.;

    Jet firstJet, secondJet;

    firstJet.ptReco = 0.;
    secondJet.ptReco = 0.;

    for(int iRecoJet = 0; iRecoJet<nReco; ++iRecoJet) {

      if( iRecoJet >= 40 ) continue;

      Jet iJet;

      iJet.eventWeight = (eventWeight<=0.) ? 1. : eventWeight;

      iJet.etaReco = etaReco[iRecoJet];
      iJet.phiReco = phiReco[iRecoJet];
      iJet.ptReco = ptReco[iRecoJet];
      iJet.eReco = eReco[iRecoJet];

      iJet.etaGen = etaGen[iRecoJet];
      iJet.phiGen = phiGen[iRecoJet];
      iJet.ptGen = ptGen[iRecoJet];
      iJet.eGen = eGen[iRecoJet];


      //response_onlyCh = 0.;

      if( algoType == "PFakt5") {
        iJet.Rch = eTracksReco[iRecoJet];
        iJet.Rgamma = ePhotonsReco[iRecoJet];
        iJet.Rnh = eNeutralHadronsReco[iRecoJet];
        iJet.Rmu = eMuonsReco[iRecoJet];
        iJet.Re = eElectronsReco[iRecoJet];
        iJet.Rhfhad = eHFhadronsReco[iRecoJet];
        iJet.Rhfem = eHFEMReco[iRecoJet];

        iJet.Ptch = ptTracksReco[iRecoJet];
        iJet.Ptgamma = ptPhotonsReco[iRecoJet];
        iJet.Ptnh = ptNeutralHadronsReco[iRecoJet];

        iJet.Nch = nTracksReco[iRecoJet];
        iJet.Ngamma = nPhotonsReco[iRecoJet];
        iJet.Nnh = nNeutralHadronsReco[iRecoJet];

        //response_onlyCh =ptTracksReco[iRecoJet]/iPtGen;

        iJet.RchGen = eTracksGen[iRecoJet];
        iJet.RgammaGen = ePhotonsGen[iRecoJet];
        iJet.RnhGen = eNeutralHadronsGen[iRecoJet];
        iJet.RmuGen = eMuonsGen[iRecoJet];

      }

      Float_t deltaEta = iJet.etaReco-iJet.etaGen;
      Float_t deltaPhi = iJet.phiReco-iJet.phiGen;
      Float_t pi = 3.14159;
      if(deltaPhi < -pi) deltaPhi += 2.*pi;
      if(deltaPhi > pi) deltaPhi -= 2.*pi;

      Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

      bool matched = false;

      if( algoType=="PFakt5" ) {
        if( deltaR<0.1 ) matched = true;
      } else {
        if( deltaR<0.3 ) matched = true;
      }

      if( (iJet.ptGen>15.)&&(iJet.ptGen<300.) ) {
        if( (fabs(iJet.etaGen)>0.)&&(fabs(iJet.etaGen)<3.) ) h1_deltaR_barrelEndcap->Fill(deltaR);
        if( (fabs(iJet.etaGen)>3.)&&(fabs(iJet.etaGen)<5.) ) h1_deltaR_forward->Fill(deltaR);
        if( (fabs(iJet.etaGen)>4.)&&(fabs(iJet.etaGen)<5.) ) h1_deltaR_veryForward->Fill(deltaR);
      } //if

      if( matched ) {

        if( iJet.ptReco > firstJet.ptReco ) {
          secondJet = firstJet;
          firstJet = iJet;
        } else if( iJet.ptReco > secondJet.ptReco ) {
          secondJet = iJet;
        }

      } //if matched

    } // for jets

    std::vector<Jet> leadingJets;
    if( firstJet.ptReco>5. ) {
      leadingJets.push_back(firstJet);
      if( secondJet.ptReco>5. ) {
        leadingJets.push_back(secondJet);
      }
    }


    for(std::vector<Jet>::const_iterator iJet = leadingJets.begin(); iJet!=leadingJets.end(); iJet++) {

      Float_t eventWeight = iJet->eventWeight;

      Float_t iEtaReco = iJet->etaReco;
      Float_t iPhiReco = iJet->phiReco;
      Float_t iEReco = iJet->eReco;
      Float_t iPtReco = iJet->ptReco;

      Float_t iEtaGen = iJet->etaGen;
      Float_t iPhiGen = iJet->phiGen;
      Float_t iEGen = iJet->eGen;
      Float_t iPtGen = iJet->ptGen;

      Int_t iPtGenBin = -1;
      for( unsigned iBin=0; iBin<(nBins-1); ++iBin) {
        if( iPtGen>Lower[iBin] && iPtGen<Lower[iBin+1] )
          iPtGenBin = iBin;
      }


      Float_t response = iPtReco/iPtGen;

      Float_t Rch = iJet->Rch;
      Float_t Rgamma = iJet->Rgamma;
      Float_t Rnh = iJet->Rnh;
      Float_t Rmu = iJet->Rmu;
      Float_t Rhfhad = iJet->Rhfhad;
      Float_t Rhfem = iJet->Rhfem;
      Float_t Re = iJet->Re;

      Float_t Ptch = iJet->Ptch;
      Float_t Ptgamma = iJet->Ptgamma;
      Float_t Ptnh = iJet->Ptnh;

      Int_t Nch = iJet->Nch;
      Int_t Ngamma = iJet->Ngamma;
      Int_t Nnh = iJet->Nnh;
      Int_t N = Nch+Ngamma+Nnh;

      Float_t RchGen = iJet->RchGen;
      Float_t RgammaGen = iJet->RgammaGen;
      Float_t RnhGen = iJet->RnhGen;
      Float_t RmuGen = iJet->RmuGen;
      Float_t ReGen = iJet->ReGen;

      Float_t Ech = Rch*iEReco;
      Float_t Egamma = Rgamma*iEReco;
      Float_t Enh = Rnh*iEReco;
      Float_t Egammanh = Egamma+Enh;

      Float_t EchGen = RchGen*iEGen;
      Float_t EgammaGen = RgammaGen*iEGen;
      Float_t EnhGen = RnhGen*iEGen;
      Float_t EgammanhGen = EgammaGen+EnhGen;


      if( iPtGen>20. )
        h2_response_vs_eta->Fill(iEtaGen, response, eventWeight);

      if( Rch==0. ) {
        hp_RchGen_vs_eta_Rch0->Fill(iEtaGen, RchGen, eventWeight);
        hp_RgammaGen_vs_eta_Rch0->Fill(iEtaGen, RgammaGen, eventWeight);
        hp_RnhGen_vs_eta_Rch0->Fill(iEtaGen, RnhGen, eventWeight);
      }

      Float_t deltaEta = fabs(iEtaReco) - fabs(iEtaGen);
      Float_t deltaPhi = iPhiReco-iPhiGen;
      Float_t pi = 3.14159;
      if(deltaPhi < -pi) deltaPhi += 2.*pi;
      if(deltaPhi > pi) deltaPhi -= 2.*pi;

      if( fabs(iEtaGen)<3. ) {
        h2_response_vs_pt_eta0_3->Fill(iPtGen, response, eventWeight);
        h2_ptGenMean_vs_pt_eta0_3->Fill(iPtGen, iPtGen, eventWeight);
      }

      //central barrel:
      if( fabs(iEtaGen)<1. ) {
        h2_response_vs_pt_centralBarrel->Fill(iPtGen, response, eventWeight);
        h2_ptGenMean_vs_pt_centralBarrel->Fill(iPtGen, iPtGen, eventWeight);
        if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
          h2_response_vs_pt_centralBarrel_Rch10_30->Fill(iPtGen, response, eventWeight);
          h2_ptGenMean_vs_pt_centralBarrel_Rch10_30->Fill(iPtGen, iPtGen, eventWeight);
        } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
          h2_response_vs_pt_centralBarrel_Rch50_70->Fill(iPtGen, response, eventWeight);
          h2_ptGenMean_vs_pt_centralBarrel_Rch50_70->Fill(iPtGen, iPtGen, eventWeight);
        } else if( (Rch >= 0.9) && (Rch<=1.) )  {
          h2_response_vs_pt_centralBarrel_Rch90_100->Fill(iPtGen, response, eventWeight);
          h2_ptGenMean_vs_pt_centralBarrel_Rch90_100->Fill(iPtGen, iPtGen, eventWeight);
        }
      }

      //full-eta PFlow plots:
      if( algoType=="PFakt5" ) {
        hp_Rch_vs_eta->Fill(iEtaGen, Rch, eventWeight);
        hp_Rgamma_vs_eta->Fill(iEtaGen, Rgamma, eventWeight);
        hp_Rnh_vs_eta->Fill(iEtaGen, Rnh, eventWeight);
        hp_Rmu_vs_eta->Fill(iEtaGen, Rmu, eventWeight);
        hp_Re_vs_eta->Fill(iEtaGen, Re, eventWeight);
        hp_Rhfhad_vs_eta->Fill(iEtaGen, Rhfhad, eventWeight);
        hp_Rhfem_vs_eta->Fill(iEtaGen, Rhfem, eventWeight);
        hp_RchGen_vs_eta->Fill(iEtaGen, RchGen, eventWeight);
        hp_RgammaGen_vs_eta->Fill(iEtaGen, RgammaGen, eventWeight);
        hp_RnhGen_vs_eta->Fill(iEtaGen, RnhGen, eventWeight);
        hp_RmuGen_vs_eta->Fill(iEtaGen, RmuGen, eventWeight);
        if( iPtGen>80. && iPtGen<90. ) {
          hp_Rch_vs_eta_80_90->Fill(iEtaGen, Rch, eventWeight);
          hp_Rgamma_vs_eta_80_90->Fill(iEtaGen, Rgamma, eventWeight);
          hp_Rnh_vs_eta_80_90->Fill(iEtaGen, Rnh, eventWeight);
          hp_Rmu_vs_eta_80_90->Fill(iEtaGen, Rmu, eventWeight);
          hp_Re_vs_eta_80_90->Fill(iEtaGen, Re, eventWeight);
          hp_Rhfhad_vs_eta_80_90->Fill(iEtaGen, Rhfhad, eventWeight);
          hp_Rhfem_vs_eta_80_90->Fill(iEtaGen, Rhfem, eventWeight);
        } //if iptgen
      } //if algotype
   
      //barrel:
      if( fabs(iEtaGen)<1.4 ) {
        h2_deltaEta_vs_pt_barrel->Fill(iPtGen, deltaEta, eventWeight);
        h2_deltaPhi_vs_pt_barrel->Fill(iPtGen, deltaPhi, eventWeight);
        h2_response_vs_pt_barrel->Fill(iPtGen, response, eventWeight);
        h2_ptGenMean_vs_pt_barrel->Fill(iPtGen, iPtGen, eventWeight);
        h2_response_vs_pt_eta0_20->Fill(iPtGen, response, eventWeight);
        h2_ptGenMean_vs_pt_eta0_20->Fill(iPtGen, iPtGen, eventWeight);

        if( algoType=="PFakt5" ) {

          //h2_response_vs_pt_barrel_onlyCh->Fill(iPtGen, response_onlyCh);

          hp_Rgamma_vs_Rch_barrel->Fill(Rch, Rgamma, eventWeight);
          hp_Rnh_vs_Rch_barrel->Fill(Rch, Rnh, eventWeight);

          hp_Rgamma_vs_Rch->Fill(Rch, Rgamma, eventWeight);
          hp_Rnh_vs_Rch->Fill(Rnh, Rgamma, eventWeight);
          hp_Rmu_vs_Rch->Fill(Rmu, Rgamma, eventWeight);

          h1_Rch->Fill(Rch, eventWeight);
          h1_Rch_barrel->Fill(Rch, eventWeight);
          if( iPtGen>50. && iPtGen<80. ) h1_Rch_barrel_50_80->Fill(Rch, eventWeight);
          hp_Rch_vs_pt->Fill(iPtGen, Rch, eventWeight);
          hp_Ptch_vs_pt->Fill(iPtGen, Ptch, eventWeight);
          hp_Nch_vs_pt->Fill(iPtGen, Nch, eventWeight);

          h1_Rgamma->Fill(Rgamma, eventWeight);
          hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma, eventWeight);
          hp_Ptgamma_vs_pt->Fill(iPtGen, Ptgamma, eventWeight);
          hp_Ngamma_vs_pt->Fill(iPtGen, Ngamma, eventWeight);

          h1_Rnh->Fill(Rnh, eventWeight);
          hp_Rnh_vs_pt->Fill(iPtGen, Rnh, eventWeight);
          hp_Ptnh_vs_pt->Fill(iPtGen, Ptnh, eventWeight);
          hp_Nnh_vs_pt->Fill(iPtGen, Nnh, eventWeight);

          h1_Rmu->Fill(Rmu, eventWeight);
          hp_Rmu_vs_pt->Fill(iPtGen, Rmu, eventWeight);

          h1_Re->Fill(Re, eventWeight);
          hp_Re_vs_pt->Fill(iPtGen, Re, eventWeight);

          h1_Rhfhad->Fill(Rhfhad, eventWeight);
          hp_Rhfhad_vs_pt->Fill(iPtGen, Rhfhad, eventWeight);

          h1_Rhfem->Fill(Rhfem, eventWeight);
          hp_Rhfem_vs_pt->Fill(iPtGen, Rhfem, eventWeight);

          hp_N_vs_pt_eta0_20->Fill(iPtGen, N, eventWeight);
          hp_N_vs_pt_barrel->Fill(iPtGen, N, eventWeight);

          hp_Rch_vs_pt_barrel->Fill(iPtGen, Rch, eventWeight);
          hp_Ptch_vs_pt_barrel->Fill(iPtGen, Ptch, eventWeight);
          hp_Nch_vs_pt_barrel->Fill(iPtGen, Nch, eventWeight);

          hp_Rgamma_vs_pt_barrel->Fill(iPtGen, Rgamma, eventWeight);
          hp_Ptgamma_vs_pt_barrel->Fill(iPtGen, Ptgamma, eventWeight);
          hp_Ngamma_vs_pt_barrel->Fill(iPtGen, Ngamma, eventWeight);

          hp_Rnh_vs_pt_barrel->Fill(iPtGen, Rnh, eventWeight);
          hp_Ptnh_vs_pt_barrel->Fill(iPtGen, Ptnh, eventWeight);
          hp_Nnh_vs_pt_barrel->Fill(iPtGen, Nnh, eventWeight);

          hp_Rmu_vs_pt_barrel->Fill(iPtGen, Rmu, eventWeight);

          hp_RgammaGen_vs_RchGen->Fill(RchGen, RgammaGen, eventWeight);
          hp_RnhGen_vs_RchGen->Fill(RnhGen, RgammaGen, eventWeight);
          hp_RmuGen_vs_RchGen->Fill(RmuGen, RgammaGen, eventWeight);

          h1_RchGen->Fill(RchGen, eventWeight);
          h1_RchGen_barrel->Fill(RchGen, eventWeight);
          if( iPtGen>50. && iPtGen<80. ) h1_RchGen_barrel_50_80->Fill(RchGen, eventWeight);
          hp_RchGen_vs_pt->Fill(iPtGen, RchGen, eventWeight);

          h1_RgammaGen->Fill(RgammaGen, eventWeight);
          hp_RgammaGen_vs_pt->Fill(iPtGen, RgammaGen, eventWeight);

          h1_RnhGen->Fill(RnhGen, eventWeight);
          hp_RnhGen_vs_pt->Fill(iPtGen, RnhGen, eventWeight);

          h1_RmuGen->Fill(RmuGen, eventWeight);
          hp_RmuGen_vs_pt->Fill(iPtGen, RmuGen, eventWeight);

          hp_RchGen_vs_pt_barrel->Fill(iPtGen, RchGen, eventWeight);
          hp_RgammaGen_vs_pt_barrel->Fill(iPtGen, RgammaGen, eventWeight);
          hp_RnhGen_vs_pt_barrel->Fill(iPtGen, RnhGen, eventWeight);
          hp_RmuGen_vs_pt_barrel->Fill(iPtGen, RmuGen, eventWeight);

          if( EchGen!=0. ) h2_Ech_vs_pt_barrel->Fill( iPtGen, Ech/EchGen, eventWeight);
          if( EgammanhGen!=0. ) h2_Egammanh_vs_pt_barrel->Fill( iPtGen, Egammanh/EgammanhGen, eventWeight);
          
          if( EchGen!=0. && EgammanhGen!=0. && iPtGenBin!=-1 ) 
            vec_h2_Egammanh_vs_Ech_barrel[iPtGenBin]->Fill( Ech/EchGen, Egammanh/EgammanhGen, eventWeight);

          if( Rch==0. ) { 
            h1_eta_Rch0->Fill( iEtaGen , eventWeight);
            h2_deltaEta_vs_pt_barrel_Rch0->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_barrel_Rch0->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_barrel_Rch0->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_barrel_Rch0->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2.5 ) {
              h2_response_vs_pt_eta0_20_Rch0->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch0->Fill(iPtGen, iPtGen, eventWeight);
            }
          } else if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
            h2_deltaEta_vs_pt_barrel_Rch10_30->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_barrel_Rch10_30->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_barrel_Rch10_30->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_barrel_Rch10_30->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2.5 ) {
              h2_response_vs_pt_eta0_20_Rch10_30->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch10_30->Fill(iPtGen, iPtGen, eventWeight);
            }
          } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
            h2_deltaEta_vs_pt_barrel_Rch50_70->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_barrel_Rch50_70->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_barrel_Rch50_70->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_barrel_Rch50_70->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2.5 ) {
              h2_response_vs_pt_eta0_20_Rch50_70->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch50_70->Fill(iPtGen, iPtGen, eventWeight);
            }
          } else if( (Rch >= 0.9) && (Rch<=1.) )  {
            h2_deltaEta_vs_pt_barrel_Rch90_100->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_barrel_Rch90_100->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_barrel_Rch90_100->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_barrel_Rch90_100->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2.5 ) {
              h2_response_vs_pt_eta0_20_Rch90_100->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch90_100->Fill(iPtGen, iPtGen, eventWeight);
            }
          }

        } //if PF 

      } //if barrel

      //endcaps:
      if( (fabs(iEtaGen)>=1.4) && (fabs(iEtaGen)<3.) ) {
        h2_deltaEta_vs_pt_endcap->Fill(iPtGen, deltaEta, eventWeight);
        h2_deltaPhi_vs_pt_endcap->Fill(iPtGen, deltaPhi, eventWeight);
        h2_response_vs_pt_endcap->Fill(iPtGen, response, eventWeight);
        h2_ptGenMean_vs_pt_endcap->Fill(iPtGen, iPtGen, eventWeight);

        if(fabs(iEtaGen)<2.) {//tracker fiducial region
          h2_response_vs_pt_eta0_20->Fill(iPtGen, response, eventWeight);
          h2_ptGenMean_vs_pt_eta0_20->Fill(iPtGen, iPtGen, eventWeight);
          h2_response_vs_pt_eta14_20->Fill(iPtGen, response, eventWeight);
          h2_ptGenMean_vs_pt_eta14_20->Fill(iPtGen, iPtGen, eventWeight);
        }

        if( algoType=="PFakt5" ) {

          if(fabs(iEtaGen)<2.) {//tracker fiducial region
            hp_N_vs_pt_eta0_20->Fill(iPtGen, N, eventWeight);
            hp_N_vs_pt_eta14_20->Fill(iPtGen, N, eventWeight);

            hp_Rgamma_vs_Rch->Fill(Rch, Rgamma, eventWeight);
            hp_Rnh_vs_Rch->Fill(Rnh, Rgamma, eventWeight);
            hp_Rmu_vs_Rch->Fill(Rmu, Rgamma, eventWeight);
          
            h1_Rch_eta14_20->Fill(Rch, eventWeight);

            hp_Rch_vs_pt->Fill(iPtGen, Rch, eventWeight);
            hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma, eventWeight);
            hp_Rnh_vs_pt->Fill(iPtGen, Rnh, eventWeight);
            hp_Re_vs_pt->Fill(iPtGen, Re, eventWeight);
            hp_Rmu_vs_pt->Fill(iPtGen, Rmu, eventWeight);
           // hp_Rhf_vs_pt->Fill(iPtGen, Rhf);

            hp_Rch_vs_pt_endcap->Fill(iPtGen, Rch, eventWeight);
            hp_Rgamma_vs_pt_endcap->Fill(iPtGen, Rgamma, eventWeight);
            hp_Rnh_vs_pt_endcap->Fill(iPtGen, Rnh, eventWeight);
           // hp_Rhf_vs_pt_endcap->Fill(iPtGen, Rhf);
            hp_Rmu_vs_pt_endcap->Fill(iPtGen, Rmu, eventWeight);
          } //tracker region

          h1_Rch->Fill(Rch, eventWeight);
          h1_Rgamma->Fill(Rgamma, eventWeight);
          h1_Rnh->Fill(Rnh, eventWeight);
          h1_Re->Fill(Re, eventWeight);
          h1_Rmu->Fill(Rmu, eventWeight);
          h1_Rhfhad->Fill(Rhfhad, eventWeight);
          h1_Rhfem->Fill(Rhfem, eventWeight);

          hp_RgammaGen_vs_RchGen->Fill(RchGen, RgammaGen, eventWeight);
          hp_RnhGen_vs_RchGen->Fill(RnhGen, RgammaGen, eventWeight);
          hp_RmuGen_vs_RchGen->Fill(RmuGen, RgammaGen, eventWeight);

          h1_RchGen->Fill(RchGen, eventWeight);
          h1_RchGen_endcap->Fill(RchGen, eventWeight);
          h1_RgammaGen->Fill(RgammaGen, eventWeight);
          h1_RnhGen->Fill(RnhGen, eventWeight);
          h1_RmuGen->Fill(RmuGen, eventWeight);

          hp_RchGen_vs_pt->Fill(iPtGen, RchGen, eventWeight);
          hp_RgammaGen_vs_pt->Fill(iPtGen, RgammaGen, eventWeight);
          hp_RnhGen_vs_pt->Fill(iPtGen, RnhGen, eventWeight);
          hp_RmuGen_vs_pt->Fill(iPtGen, RmuGen, eventWeight);

          hp_RchGen_vs_pt_endcap->Fill(iPtGen, RchGen, eventWeight);
          hp_RgammaGen_vs_pt_endcap->Fill(iPtGen, RgammaGen, eventWeight);
          hp_RnhGen_vs_pt_endcap->Fill(iPtGen, RnhGen, eventWeight);
          hp_RmuGen_vs_pt_endcap->Fill(iPtGen, RmuGen, eventWeight);

          if( EchGen!=0. ) h2_Ech_vs_pt_endcap->Fill( iPtGen, Ech/EchGen, eventWeight);
          if( EgammanhGen!=0. ) h2_Egammanh_vs_pt_endcap->Fill( iPtGen, Egammanh/EgammanhGen, eventWeight);


          if( Rch==0. ) { 
            h1_eta_Rch0->Fill( iEtaGen , eventWeight);
            h2_deltaEta_vs_pt_endcap_Rch0->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_endcap_Rch0->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_endcap_Rch0->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_endcap_Rch0->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2. ) {
              h2_response_vs_pt_eta0_20_Rch0->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch0->Fill(iPtGen, iPtGen, eventWeight);
              h2_response_vs_pt_eta14_20_Rch0->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta14_20_Rch0->Fill(iPtGen, iPtGen, eventWeight);
            }
          } else if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
            h2_deltaEta_vs_pt_endcap_Rch10_30->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_endcap_Rch10_30->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_endcap_Rch10_30->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_endcap_Rch10_30->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2. ) {
              h2_response_vs_pt_eta0_20_Rch10_30->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch10_30->Fill(iPtGen, iPtGen, eventWeight);
              h2_response_vs_pt_eta14_20_Rch10_30->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta14_20_Rch10_30->Fill(iPtGen, iPtGen, eventWeight);
            }
          } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
            h2_deltaEta_vs_pt_endcap_Rch50_70->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_endcap_Rch50_70->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_endcap_Rch50_70->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_endcap_Rch50_70->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2. ) {
              h2_response_vs_pt_eta0_20_Rch50_70->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch50_70->Fill(iPtGen, iPtGen, eventWeight);
              h2_response_vs_pt_eta14_20_Rch50_70->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta14_20_Rch50_70->Fill(iPtGen, iPtGen, eventWeight);
            }
          } else if( (Rch >= 0.9) && (Rch <=1.) )  {
            h2_deltaEta_vs_pt_endcap_Rch90_100->Fill(iPtGen, deltaEta, eventWeight);
            h2_deltaPhi_vs_pt_endcap_Rch90_100->Fill(iPtGen, deltaPhi, eventWeight);
            h2_response_vs_pt_endcap_Rch90_100->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_endcap_Rch90_100->Fill(iPtGen, iPtGen, eventWeight);
            if( fabs(iEtaGen)<2. ) {
              h2_response_vs_pt_eta0_20_Rch90_100->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta0_20_Rch90_100->Fill(iPtGen, iPtGen, eventWeight);
              h2_response_vs_pt_eta14_20_Rch90_100->Fill(iPtGen, response, eventWeight);
              h2_ptGenMean_vs_pt_eta14_20_Rch90_100->Fill(iPtGen, iPtGen, eventWeight);
            }
          }

        } //if PF 

      } //if endcaps

      if( (fabs(iEtaGen)>2.6)&&(fabs(iEtaGen)<3.2) ) {
        h2_response_vs_pt_eta26_32->Fill(iPtGen, response, eventWeight);
        h2_ptGenMean_vs_pt_eta26_32->Fill(iPtGen, iPtGen, eventWeight);

        if( algoType=="PFakt5" ) {
          if( Rch==0. ) {
            h2_response_vs_pt_eta26_32_Rch0->Fill(iPtGen, response, eventWeight);
            h2_ptGenMean_vs_pt_eta26_32_Rch0->Fill(iPtGen, iPtGen, eventWeight);
          }
        }
      }

      if( (fabs(iEtaGen)>=3.) && (fabs(iEtaGen)<5.) ) {
        h2_deltaEta_vs_pt_forward->Fill(iPtGen, deltaEta, eventWeight);
        h2_deltaPhi_vs_pt_forward->Fill(iPtGen, deltaPhi, eventWeight);
        h2_response_vs_pt_forward->Fill(iPtGen, response, eventWeight);

        if( algoType=="PFakt5" ) {

          h1_Rch->Fill(Rch, eventWeight);
          h1_Rgamma->Fill(Rgamma, eventWeight);
          h1_Rnh->Fill(Rnh, eventWeight);
          h1_Re->Fill(Re, eventWeight);
          h1_Rmu->Fill(Rmu, eventWeight);
          h1_Rhfhad->Fill(Rhfhad, eventWeight);
          h1_Rhfem->Fill(Rhfem, eventWeight);

          hp_RgammaGen_vs_RchGen->Fill(RchGen, RgammaGen, eventWeight);
          hp_RnhGen_vs_RchGen->Fill(RnhGen, RgammaGen, eventWeight);
          hp_RmuGen_vs_RchGen->Fill(RmuGen, RgammaGen, eventWeight);

          h1_RchGen->Fill(RchGen, eventWeight);
          hp_RchGen_vs_pt->Fill(iPtGen, RchGen, eventWeight);

          h1_RgammaGen->Fill(RgammaGen, eventWeight);
          hp_RgammaGen_vs_pt->Fill(iPtGen, RgammaGen, eventWeight);

          h1_RnhGen->Fill(RnhGen, eventWeight);
          hp_RnhGen_vs_pt->Fill(iPtGen, RnhGen, eventWeight);

          h1_RmuGen->Fill(RmuGen, eventWeight);
          hp_RmuGen_vs_pt->Fill(iPtGen, RmuGen, eventWeight);


        } //if PF


      } //if forward


    } //for recoJets


  } //for entries


  
  if( draw ) {
    std::string name;
    name = algoType + "_barrel";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel, h2_ptGenMean_vs_pt_barrel, "response", "pt", "barrel", "", algoType, outFileName, name);
    name = algoType + "_centralBarrel";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_centralBarrel, h2_ptGenMean_vs_pt_centralBarrel, "response", "pt", "centralBarrel", "", algoType, outFileName, name);
    name = algoType + "_endcap";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap, h2_ptGenMean_vs_pt_endcap, "response", "pt", "endcap", "", algoType, outFileName, name);
    name = algoType + "_eta0_20";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20, h2_ptGenMean_vs_pt_eta0_20, "response", "pt", "eta0_20", "", algoType, outFileName, name);
    name = algoType + "_eta14_20";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20, h2_ptGenMean_vs_pt_eta14_20, "response", "pt", "eta14_20", "", algoType, outFileName, name);
    name = algoType + "_eta26_32";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta26_32, h2_ptGenMean_vs_pt_eta26_32, "response", "pt", "eta26_32", "", algoType, outFileName, name);
    name = algoType + "_eta0_3";
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_3, h2_ptGenMean_vs_pt_eta0_3, "response", "pt", "eta0_3", "", algoType, outFileName, name);
  } else {
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel, h2_ptGenMean_vs_pt_barrel, "response", "pt", "barrel", "", algoType, outFileName);
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_centralBarrel, h2_ptGenMean_vs_pt_centralBarrel, "response", "pt", "centralBarrel", "", algoType, outFileName);
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap, h2_ptGenMean_vs_pt_endcap, "response", "pt", "endcap", "", algoType, outFileName);
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20, h2_ptGenMean_vs_pt_eta0_20,"response",  "pt", "eta0_20", "", algoType, outFileName);
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20, h2_ptGenMean_vs_pt_eta14_20, "response", "pt", "eta14_20", "", algoType, outFileName);
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta26_32, h2_ptGenMean_vs_pt_eta26_32, "response", "pt", "eta26_32", "", algoType, outFileName);
    fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_3, h2_ptGenMean_vs_pt_eta0_3, "response", "pt", "eta0_3", "", algoType, outFileName);
  }


  if( algoType=="PFakt5" ) {
    if( draw ) {
      std::string name;
      name = "Ech_"+algoType+"_barrel";
      fitTools::fitDistribution_TGraph(h2_Ech_vs_pt_barrel, h2_ptGenMean_vs_pt_barrel, "Ech", "pt", "barrel", "", algoType, outFileName, name, 1.2);
      name = "Ech_"+algoType+"_endcap";
      fitTools::fitDistribution_TGraph(h2_Ech_vs_pt_endcap, h2_ptGenMean_vs_pt_endcap, "Ech", "pt", "endcap", "", algoType, outFileName, name, 1.2);
      name = "Egammanh_"+algoType+"_barrel";
      fitTools::fitDistribution_TGraph(h2_Egammanh_vs_pt_barrel, h2_ptGenMean_vs_pt_barrel, "Egammanh", "pt", "barrel", "", algoType, outFileName, name, 1.2);
      name = "Egammanh_"+algoType+"_endcap";
      fitTools::fitDistribution_TGraph(h2_Egammanh_vs_pt_endcap, h2_ptGenMean_vs_pt_endcap, "Egammanh", "pt", "endcap", "", algoType, outFileName, name, 1.2);

      name = "response_" + algoType + "_barrel_Rch0";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch0, h2_ptGenMean_vs_pt_barrel_Rch0, "response", "pt", "barrel", "Rch0", algoType, outFileName, name);
      name = "response_" + algoType + "_barrel_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch10_30, h2_ptGenMean_vs_pt_barrel_Rch10_30, "response", "pt", "barrel", "Rch10_30", algoType, outFileName, name);
      name = "response_" + algoType + "_barrel_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch50_70, h2_ptGenMean_vs_pt_barrel_Rch50_70, "response", "pt", "barrel", "Rch50_70", algoType, outFileName, name);
      name = "response_" + algoType + "_barrel_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch90_100, h2_ptGenMean_vs_pt_barrel_Rch90_100, "response", "pt", "barrel", "Rch90_100", algoType, outFileName, name);

      name = "response_" + algoType + "_centralBarrel_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_centralBarrel_Rch10_30, h2_ptGenMean_vs_pt_centralBarrel_Rch10_30, "response", "pt", "centralBarrel", "Rch10_30", algoType, outFileName, name);
      name = "response_" + algoType + "_centralBarrel_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_centralBarrel_Rch50_70, h2_ptGenMean_vs_pt_centralBarrel_Rch50_70, "response", "pt", "centralBarrel", "Rch50_70", algoType, outFileName, name);
      name = "response_" + algoType + "_centralBarrel_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_centralBarrel_Rch90_100, h2_ptGenMean_vs_pt_centralBarrel_Rch90_100, "response", "pt", "centralBarrel", "Rch90_100", algoType, outFileName, name);

      name = "response_" + algoType + "_endcap_Rch0";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch0, h2_ptGenMean_vs_pt_endcap_Rch0, "response", "pt", "endcap", "Rch0", algoType, outFileName, name);
      name = "response_" + algoType + "_endcap_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch10_30, h2_ptGenMean_vs_pt_endcap_Rch10_30, "response", "pt", "endcap", "Rch10_30", algoType, outFileName, name);
      name = "response_" + algoType + "_endcap_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch50_70, h2_ptGenMean_vs_pt_endcap_Rch50_70, "response", "pt", "endcap", "Rch50_70", algoType, outFileName, name);
      name = "response_" + algoType + "_endcap_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch90_100, h2_ptGenMean_vs_pt_endcap_Rch90_100, "response", "pt", "endcap", "Rch90_100", algoType, outFileName, name);

      name = "response_" + algoType + "_eta0_20_Rch0";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch0, h2_ptGenMean_vs_pt_eta0_20_Rch0, "response", "pt", "eta0_20", "Rch0", algoType, outFileName, name);
      name = "response_" + algoType + "_eta0_20_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch10_30, h2_ptGenMean_vs_pt_eta0_20_Rch10_30, "response", "pt", "eta0_20", "Rch10_30", algoType, outFileName, name);
      name = "response_" + algoType + "_eta0_20_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch50_70, h2_ptGenMean_vs_pt_eta0_20_Rch50_70, "response", "pt", "eta0_20", "Rch50_70", algoType, outFileName, name);
      name = "response_" + algoType + "_eta0_20_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch90_100, h2_ptGenMean_vs_pt_eta0_20_Rch90_100, "response", "pt", "eta0_20", "Rch90_100", algoType, outFileName, name);

      name = "response_" + algoType + "_eta14_20_Rch0";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch0, h2_ptGenMean_vs_pt_eta14_20_Rch0,"response",  "pt", "eta14_20", "Rch0", algoType, outFileName, name);
      name = "response_" + algoType + "_eta14_20_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch10_30, h2_ptGenMean_vs_pt_eta14_20_Rch10_30, "response", "pt", "eta14_20", "Rch10_30", algoType, outFileName, name);
      name = "response_" + algoType + "_eta14_20_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch50_70, h2_ptGenMean_vs_pt_eta14_20_Rch50_70, "response", "pt", "eta14_20", "Rch50_70", algoType, outFileName, name);
      name = "response_" + algoType + "_eta14_20_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch90_100, h2_ptGenMean_vs_pt_eta14_20_Rch90_100, "response", "pt", "eta14_20", "Rch90_100", algoType, outFileName, name);

      name = "response_" + algoType + "_eta26_32_Rch0";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta26_32_Rch0, h2_ptGenMean_vs_pt_eta26_32_Rch0, "response", "pt", "eta26_32", "Rch0", algoType, outFileName, name);

      name = "deltaEta_" + algoType + "_barrel_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch10_30, h2_ptGenMean_vs_pt_barrel_Rch10_30, "deltaEta", "pt", "barrel", "Rch10_30", algoType, outFileName, name);
      name = "deltaEta_" + algoType + "_barrel_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch50_70, h2_ptGenMean_vs_pt_barrel_Rch50_70, "deltaEta", "pt", "barrel", "Rch50_70", algoType, outFileName, name);
      name = "deltaEta_" + algoType + "_barrel_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch90_100, h2_ptGenMean_vs_pt_barrel_Rch90_100, "deltaEta", "pt", "barrel", "Rch90_100", algoType, outFileName, name);

      name = "deltaPhi_" + algoType + "_barrel_Rch10_30";
      fitTools::fitDistribution_TGraph(h2_deltaPhi_vs_pt_barrel_Rch10_30, h2_ptGenMean_vs_pt_barrel_Rch10_30, "deltaPhi", "pt", "barrel", "Rch10_30", algoType, outFileName, name);
      name = "deltaPhi_" + algoType + "_barrel_Rch50_70";
      fitTools::fitDistribution_TGraph(h2_deltaPhi_vs_pt_barrel_Rch50_70, h2_ptGenMean_vs_pt_barrel_Rch50_70, "deltaPhi", "pt", "barrel", "Rch50_70", algoType, outFileName, name);
      name = "deltaPhi_" + algoType + "_barrel_Rch90_100";
      fitTools::fitDistribution_TGraph(h2_deltaPhi_vs_pt_barrel_Rch90_100, h2_ptGenMean_vs_pt_barrel_Rch90_100, "deltaPhi", "pt", "barrel", "Rch90_100", algoType, outFileName, name);

    } else {

      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch0, h2_ptGenMean_vs_pt_barrel_Rch0, "response", "pt", "barrel", "Rch0", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch10_30, h2_ptGenMean_vs_pt_barrel_Rch10_30, "response", "pt", "barrel", "Rch10_30", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch50_70, h2_ptGenMean_vs_pt_barrel_Rch50_70, "response", "pt", "barrel", "Rch50_70", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch90_100, h2_ptGenMean_vs_pt_barrel_Rch90_100, "response", "pt", "barrel", "Rch90_100", algoType, outFileName);

      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch0, h2_ptGenMean_vs_pt_endcap_Rch0, "response", "pt", "endcap", "Rch0", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch10_30, h2_ptGenMean_vs_pt_endcap_Rch10_30, "response", "pt", "endcap", "Rch10_30", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch50_70, h2_ptGenMean_vs_pt_endcap_Rch50_70, "response", "pt", "endcap", "Rch50_70", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_endcap_Rch90_100, h2_ptGenMean_vs_pt_endcap_Rch90_100, "response", "pt", "endcap", "Rch90_100", algoType, outFileName);

      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch0, h2_ptGenMean_vs_pt_eta0_20_Rch0, "response", "pt", "eta0_20", "Rch0", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch10_30, h2_ptGenMean_vs_pt_eta0_20_Rch10_30, "response", "pt", "eta0_20", "Rch10_30", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch50_70, h2_ptGenMean_vs_pt_eta0_20_Rch50_70, "response", "pt", "eta0_20", "Rch50_70", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta0_20_Rch90_100, h2_ptGenMean_vs_pt_eta0_20_Rch90_100, "response", "pt", "eta0_20", "Rch90_100", algoType, outFileName);

      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch0, h2_ptGenMean_vs_pt_eta14_20_Rch0, "response", "pt", "eta14_20", "Rch0", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch10_30, h2_ptGenMean_vs_pt_eta14_20_Rch10_30,"response",  "pt", "eta14_20", "Rch10_30", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch50_70, h2_ptGenMean_vs_pt_eta14_20_Rch50_70, "response", "pt", "eta14_20", "Rch50_70", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_eta14_20_Rch90_100, h2_ptGenMean_vs_pt_eta14_20_Rch90_100, "response", "pt", "eta14_20", "Rch90_100", algoType, outFileName);

      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch10_30, h2_ptGenMean_vs_pt_barrel_Rch10_30, "response", "pt", "barrel", "Rch10_30", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch50_70, h2_ptGenMean_vs_pt_barrel_Rch50_70, "sigmaEta", "pt", "barrel", "Rch50_70", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_response_vs_pt_barrel_Rch90_100, h2_ptGenMean_vs_pt_barrel_Rch90_100, "sigmaEta", "pt", "barrel", "Rch90_100", algoType, outFileName);

      fitTools::fitDistribution_TGraph(h2_deltaPhi_vs_pt_barrel_Rch10_30, h2_ptGenMean_vs_pt_barrel_Rch10_30, "sigmaPhi", "pt", "barrel", "Rch10_30", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_deltaPhi_vs_pt_barrel_Rch50_70, h2_ptGenMean_vs_pt_barrel_Rch50_70, "sigmaPhi", "pt", "barrel", "Rch50_70", algoType, outFileName);
      fitTools::fitDistribution_TGraph(h2_deltaPhi_vs_pt_barrel_Rch90_100, h2_ptGenMean_vs_pt_barrel_Rch90_100, "sigmaPhi", "pt", "barrel", "Rch90_100", algoType, outFileName);
    }

  }





  histoName = "response_vs_eta_FIT_"+algoType;
  TH1F* h1_response_vs_eta_FIT = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5.5, 5.5);
  h1_response_vs_eta_FIT->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_FIT->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_eta_FIT->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_eta_FIT_"+algoType;
  TH1F* h1_resolution_vs_eta_FIT = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet #eta", 100, -5.5, 5.5);
  h1_resolution_vs_eta_FIT->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_FIT->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_FIT->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_eta_MEAN_"+algoType;
  TH1F* h1_response_vs_eta_MEAN = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5.5, 5.5);
  h1_response_vs_eta_MEAN->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_MEAN->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_response_vs_eta_MEAN->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_eta_RMS_"+algoType;
  TH1F* h1_resolution_vs_eta_RMS = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet #eta", 100, -5.5, 5.5);
  h1_resolution_vs_eta_RMS->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_RMS->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_RMS->GetYaxis()->SetRangeUser(0., 0.5);


  fitTools::fillProfile(h1_response_vs_eta_FIT, h1_resolution_vs_eta_FIT, h1_response_vs_eta_MEAN, h1_resolution_vs_eta_RMS, h2_response_vs_eta);
  

 
  delete h2_response_vs_pt_barrel_Rch0;
  h2_response_vs_pt_barrel_Rch0 = 0;
  delete h2_response_vs_pt_endcap_Rch0;
  h2_response_vs_pt_endcap_Rch0 = 0;
  delete h2_response_vs_pt_barrel_Rch10_30;
  h2_response_vs_pt_barrel_Rch10_30 = 0;
  delete h2_response_vs_pt_endcap_Rch10_30;
  h2_response_vs_pt_endcap_Rch10_30 = 0;
  delete h2_response_vs_pt_barrel_Rch50_70;
  h2_response_vs_pt_barrel_Rch50_70 = 0;
  delete h2_response_vs_pt_endcap_Rch50_70;
  h2_response_vs_pt_endcap_Rch50_70 = 0;
  delete h2_response_vs_pt_barrel_Rch90_100;
  h2_response_vs_pt_barrel_Rch90_100 = 0;
  delete h2_response_vs_pt_endcap_Rch90_100;
  h2_response_vs_pt_endcap_Rch90_100 = 0;

  delete h2_deltaEta_vs_pt_barrel;
  h2_deltaEta_vs_pt_barrel = 0;
  delete h2_deltaEta_vs_pt_endcap;
  h2_deltaEta_vs_pt_endcap = 0;
  delete h2_deltaEta_vs_pt_forward;
  h2_deltaEta_vs_pt_forward = 0;

  delete h2_deltaEta_vs_pt_barrel_Rch0;
  h2_deltaEta_vs_pt_barrel_Rch0 = 0;
  delete h2_deltaEta_vs_pt_endcap_Rch0;
  h2_deltaEta_vs_pt_endcap_Rch0 = 0;
  delete h2_deltaEta_vs_pt_barrel_Rch10_30;
  h2_deltaEta_vs_pt_barrel_Rch10_30 = 0;
  delete h2_deltaEta_vs_pt_endcap_Rch10_30;
  h2_deltaEta_vs_pt_endcap_Rch10_30 = 0;
  delete h2_deltaEta_vs_pt_barrel_Rch50_70;
  h2_deltaEta_vs_pt_barrel_Rch50_70 = 0;
  delete h2_deltaEta_vs_pt_endcap_Rch50_70;
  h2_deltaEta_vs_pt_endcap_Rch50_70 = 0;
  delete h2_deltaEta_vs_pt_barrel_Rch90_100;
  h2_deltaEta_vs_pt_barrel_Rch90_100 = 0;
  delete h2_deltaEta_vs_pt_endcap_Rch90_100;
  h2_deltaEta_vs_pt_endcap_Rch90_100 = 0;

  delete h2_deltaPhi_vs_pt_barrel;
  h2_deltaPhi_vs_pt_barrel = 0;
  delete h2_deltaPhi_vs_pt_endcap;
  h2_deltaPhi_vs_pt_endcap = 0;
  delete h2_deltaPhi_vs_pt_forward;
  h2_deltaPhi_vs_pt_forward = 0;

  delete h2_deltaPhi_vs_pt_barrel_Rch0;
  h2_deltaPhi_vs_pt_barrel_Rch0 = 0;
  delete h2_deltaPhi_vs_pt_endcap_Rch0;
  h2_deltaPhi_vs_pt_endcap_Rch0 = 0;
  delete h2_deltaPhi_vs_pt_barrel_Rch10_30;
  h2_deltaPhi_vs_pt_barrel_Rch10_30 = 0;
  delete h2_deltaPhi_vs_pt_endcap_Rch10_30;
  h2_deltaPhi_vs_pt_endcap_Rch10_30 = 0;
  delete h2_deltaPhi_vs_pt_barrel_Rch50_70;
  h2_deltaPhi_vs_pt_barrel_Rch50_70 = 0;
  delete h2_deltaPhi_vs_pt_endcap_Rch50_70;
  h2_deltaPhi_vs_pt_endcap_Rch50_70 = 0;
  delete h2_deltaPhi_vs_pt_barrel_Rch90_100;
  h2_deltaPhi_vs_pt_barrel_Rch90_100 = 0;
  delete h2_deltaPhi_vs_pt_endcap_Rch90_100;
  h2_deltaPhi_vs_pt_endcap_Rch90_100 = 0;



  TFile* outFile = TFile::Open(outFileName.c_str(), "UPDATE");
  outFile->cd();


  for( unsigned ih2=0; ih2<nBins; ++ih2 )
    vec_h2_Egammanh_vs_Ech_barrel[ih2]->Write();


  h2_response_vs_pt_barrel->Write();
  h2_response_vs_pt_endcap->Write();
  h2_response_vs_pt_forward->Write();

  h2_ptGenMean_vs_pt_barrel->Write();
  h2_ptGenMean_vs_pt_endcap->Write();

  h1_deltaR_barrelEndcap->Write();
  h1_deltaR_forward->Write();
  h1_deltaR_veryForward->Write();

  h1_response_vs_eta_MEAN->Write();
  h1_response_vs_eta_FIT->Write();
  h1_resolution_vs_eta_RMS->Write();
  h1_resolution_vs_eta_FIT->Write();

  if( algoType == "PFakt5" ) { 

    hp_Rgamma_vs_Rch_barrel->Write();
    hp_Rnh_vs_Rch_barrel->Write();

    hp_Rgamma_vs_Rch->Write();
    hp_Rnh_vs_Rch->Write();
    hp_Rmu_vs_Rch->Write();

    h1_Rch->Write();
    h1_Rch_barrel->Write();
    h1_Rch_barrel_50_80->Write();
    h1_Rch_eta14_20->Write();
    hp_Rch_vs_eta->Write();
    hp_Rch_vs_eta_80_90->Write();
    hp_Rch_vs_pt->Write();

    h1_Rgamma->Write();
    hp_Rgamma_vs_eta->Write();
    hp_Rgamma_vs_eta_80_90->Write();
    hp_Rgamma_vs_pt->Write();

    h1_Rnh->Write();
    hp_Rnh_vs_eta->Write();
    hp_Rnh_vs_eta_80_90->Write();
    hp_Rnh_vs_pt->Write();

    h1_Re->Write();
    hp_Re_vs_eta->Write();
    hp_Re_vs_eta_80_90->Write();
    hp_Re_vs_pt->Write();

    h1_Rmu->Write();
    hp_Rmu_vs_eta->Write();
    hp_Rmu_vs_eta_80_90->Write();
    hp_Rmu_vs_pt->Write();

    h1_Rhfhad->Write();
    hp_Rhfhad_vs_eta->Write();
    hp_Rhfhad_vs_eta_80_90->Write();
    hp_Rhfhad_vs_pt->Write();

    h1_Rhfem->Write();
    hp_Rhfem_vs_eta->Write();
    hp_Rhfem_vs_eta_80_90->Write();
    hp_Rhfem_vs_pt->Write();

    hp_Ptch_vs_pt->Write();
    hp_Ptgamma_vs_pt->Write();
    hp_Ptnh_vs_pt->Write();

    hp_N_vs_pt_eta0_20->Write();
    hp_N_vs_pt_barrel->Write();
    hp_N_vs_pt_eta14_20->Write();

    hp_Nch_vs_pt->Write();
    hp_Ngamma_vs_pt->Write();
    hp_Nnh_vs_pt->Write();

    hp_Rch_vs_pt_barrel->Write();
    hp_Rgamma_vs_pt_barrel->Write();
    hp_Rnh_vs_pt_barrel->Write();
    hp_Rmu_vs_pt_barrel->Write();

    hp_Ptch_vs_pt_barrel->Write();
    hp_Ptgamma_vs_pt_barrel->Write();
    hp_Ptnh_vs_pt_barrel->Write();

    hp_Nch_vs_pt_barrel->Write();
    hp_Ngamma_vs_pt_barrel->Write();
    hp_Nnh_vs_pt_barrel->Write();

    hp_Rch_vs_pt_endcap->Write();
    hp_Rgamma_vs_pt_endcap->Write();
    hp_Rnh_vs_pt_endcap->Write();
    hp_Rhfhad_vs_pt_endcap->Write();
    hp_Rhfem_vs_pt_endcap->Write();
    hp_Rmu_vs_pt_endcap->Write();

    hp_RgammaGen_vs_RchGen->Write();
    hp_RnhGen_vs_RchGen->Write();
    hp_RmuGen_vs_RchGen->Write();

    h1_RchGen->Write();
    h1_RchGen_barrel->Write();
    h1_RchGen_barrel_50_80->Write();
    h1_RchGen_endcap->Write();
    hp_RchGen_vs_eta->Write();
    hp_RchGen_vs_eta_Rch0->Write();
    hp_RchGen_vs_pt->Write();

    h1_RgammaGen->Write();
    hp_RgammaGen_vs_eta->Write();
    hp_RgammaGen_vs_eta_Rch0->Write();
    hp_RgammaGen_vs_pt->Write();

    h1_RnhGen->Write();
    hp_RnhGen_vs_eta->Write();
    hp_RnhGen_vs_eta_Rch0->Write();
    hp_RnhGen_vs_pt->Write();

    h1_RmuGen->Write();
    hp_RmuGen_vs_eta->Write();
    hp_RmuGen_vs_pt->Write();

    hp_RchGen_vs_pt_barrel->Write();
    hp_RgammaGen_vs_pt_barrel->Write();
    hp_RnhGen_vs_pt_barrel->Write();
    hp_RmuGen_vs_pt_barrel->Write();

    hp_RchGen_vs_pt_endcap->Write();
    hp_RgammaGen_vs_pt_endcap->Write();
    hp_RnhGen_vs_pt_endcap->Write();
    hp_RmuGen_vs_pt_endcap->Write();

  } //if PF


  outFile->Write();
  outFile->Close();

} //fillFile


