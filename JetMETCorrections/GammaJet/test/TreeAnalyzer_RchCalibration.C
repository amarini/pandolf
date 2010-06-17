#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>

#include "fitTools.C"


bool DOFITS_ = true;
bool NOWEIGHTS_ = true;


class FitResults {
 public:
  Float_t p0;
  Float_t p1;
  Float_t p2;
  Float_t p3;
  Float_t p4;
  Float_t p5;
};


std::vector<FitResults> fitTheHellOutOfThem(std::vector<TH2D*> vec_h2, std::vector<TH2D*> vec_h2_Rch, std::string name, const std::string& dataset);
std::vector< TH2D* > bookResponseHistoVector(const std::string& name, Int_t nBins);
std::vector< TH2D* > bookRchTH2Vector(const std::string& name, Int_t nBins);
std::vector< TH1F* > bookRchTH1Vector(const std::string& name, Int_t nBins);


bool DEBUG_ = false;
Int_t _nBins_Rch=15;



void RchCalibration(std::string dataset) {


  std::string infile_dataset = dataset;

  TChain* tree = new TChain("jetTree");
  
  std::string infileName;

  if( dataset=="PhotonJet") {

    tree->Add("PhotonJet_2ndLevelTree_PhotonJet15.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet30.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet80.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet170.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet300.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet470.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet800.root/jetTree");

  } else if( dataset=="PhotonJet_Summer09") {

    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt15.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt30.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt80.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt170.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt300.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt470.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt800.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet_Summer09_Pt1400.root/jetTree");

  } else {

    if( dataset=="" )
      infileName = "PhotonJet_2ndLevelTree.root/jetTree";
    else
      infileName = "PhotonJet_2ndLevelTree_"+infile_dataset+".root/jetTree";

    std::cout << "Opening file: " << infileName << std::endl;
    tree->Add(infileName.c_str());

  }


  std::string fileName;
  if( DEBUG_ )
    fileName = "provaRchCalibration.root";
  else {
    if( dataset=="" ) fileName = "RchCalibration";
    else fileName = "RchCalibration_"+dataset;
  }

  if( NOWEIGHTS_ ) fileName = fileName + "_NOWEIGHTS";
  else fileName = fileName + "_CHI";

  fileName = fileName + ".root";

  TFile* file = TFile::Open(fileName.c_str(), "RECREATE");
  file->Close();

 
  TH1F* h1_deltaR_phot_2ndJet = new TH1F("deltaR_phot_2ndJet", "", 100, 0., 1.);


  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

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
  Float_t ePhotonsReco;
  tree->SetBranchAddress("ePhotonsReco", &ePhotonsReco);
  Float_t eNeutralHadronsReco;
  tree->SetBranchAddress("eNeutralHadronsReco", &eNeutralHadronsReco);

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
  Float_t ePhotonsGen;
  tree->SetBranchAddress("ePhotonsGen", &ePhotonsGen);
  Float_t eNeutralHadronsGen;
  tree->SetBranchAddress("eNeutralHadronsGen", &eNeutralHadronsGen);

  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);

  
  std::string histoName;
  


  //for pT spectrums, variable binning to avoid huge error bars at high pT:
  int nBins = 26;

  Double_t Lower[nBins];

  Int_t i=0;
  Lower[i++] = 10.;
  Lower[i++] = 15.;
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

/*
  Double_t Lower_Rch[_nBins_Rch];

  for( int i_Rch=0; i_Rch<_nBins_Rch; ++i_Rch )
    Lower_Rch[i_Rch]= (double)i_Rch/(double)_nBins_Rch;
*/


  // histogram booking:
  std::vector< TH2D* > vec_h2_response_vs_Rch_barrel    = bookResponseHistoVector("response_vs_Rch_barrel", nBins);
  std::vector< TH2D* > vec_h2_response_vs_Rch_endcap    = bookResponseHistoVector("response_vs_Rch_endcap", nBins);
  std::vector< TH2D* > vec_h2_response_vs_Rch_eta0_3    = bookResponseHistoVector("response_vs_Rch_eta0_3", nBins);
  std::vector< TH2D* > vec_h2_response_vs_Rch_eta14_20  = bookResponseHistoVector("response_vs_Rch_eta14_20", nBins);
  std::vector< TH2D* > vec_h2_response_vs_Rch_eta155_3  = bookResponseHistoVector("response_vs_Rch_eta155_3", nBins);
  std::vector< TH2D* > vec_h2_responseGEN_vs_Rch_barrel = bookResponseHistoVector("responseGEN_vs_Rch_barrel", nBins);
  std::vector< TH2D* > vec_h2_responseGEN_vs_Rch_endcap = bookResponseHistoVector("responseGEN_vs_Rch_endcap", nBins);
  std::vector< TH2D* > vec_h2_responseGEN_vs_Rch_eta14_20 = bookResponseHistoVector("responseGEN_vs_Rch_eta14_20", nBins);
  std::vector< TH2D* > vec_h2_responseGEN_vs_RchGEN_barrel = bookResponseHistoVector("responseGEN_vs_RchGEN_barrel", nBins);
  std::vector< TH2D* > vec_h2_responseGEN_vs_RchGEN_endcap = bookResponseHistoVector("responseGEN_vs_RchGEN_endcap", nBins);

  std::vector< TH2D* > vec_h2_Rch_barrel =       bookRchTH2Vector("h2_Rch_barrel", nBins);
  std::vector< TH2D* > vec_h2_Rch_endcap =       bookRchTH2Vector("h2_Rch_endcap", nBins);
  std::vector< TH2D* > vec_h2_Rch_eta14_20 =     bookRchTH2Vector("h2_Rch_eta14_20", nBins);
  std::vector< TH2D* > vec_h2_Rch_eta0_3 =       bookRchTH2Vector("h2_Rch_eta0_3", nBins);
  std::vector< TH2D* > vec_h2_Rch_eta155_3 =     bookRchTH2Vector("h2_Rch_eta155_3", nBins);
  std::vector< TH2D* > vec_h2_Rch_barrelGEN =    bookRchTH2Vector("h2_Rch_barrelGEN", nBins);
  std::vector< TH2D* > vec_h2_Rch_endcapGEN =    bookRchTH2Vector("h2_Rch_endcapGEN", nBins);
  std::vector< TH2D* > vec_h2_Rch_eta14_20GEN =  bookRchTH2Vector("h2_Rch_eta14_20GEN", nBins);
  std::vector< TH2D* > vec_h2_RchGEN_barrelGEN = bookRchTH2Vector("h2_RchGEN_barrelGEN", nBins);
  std::vector< TH2D* > vec_h2_RchGEN_endcapGEN = bookRchTH2Vector("h2_RchGEN_endcapGEN", nBins);

  //Rch distributions for all pt bins:
  std::vector< TH1F* > vec_h1_RchGEN_barrelGEN = bookRchTH1Vector("RchGEN_barrelGEN", nBins);
  std::vector< TH1F* > vec_h1_RchGEN_endcapGEN = bookRchTH1Vector("RchGEN_endcapGEN", nBins);


  std::string algoType = "PFakt5";

  histoName = "response_vs_ptPhot_barrel_PRE_"+algoType;
  TH2D* h2_response_vs_ptPhot_barrel_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_barrel_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_barrel_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_barrel_PRE->Sumw2();

  histoName = "response_vs_ptPhot_endcap_PRE_"+algoType;
  TH2D* h2_response_vs_ptPhot_endcap_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_endcap_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_endcap_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_endcap_PRE->Sumw2();

  histoName = "responseGEN_vs_ptGen_barrel_PRE_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_barrel_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_barrel_PRE->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_barrel_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_barrel_PRE->Sumw2();

  histoName = "responseGEN_vs_ptGen_endcap_PRE_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_endcap_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_endcap_PRE->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_endcap_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_endcap_PRE->Sumw2();

  histoName = "responseGEN_vs_ptGen_eta14_20_PRE_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_eta14_20_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_eta14_20_PRE->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_eta14_20_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_eta14_20_PRE->Sumw2();

  histoName = "ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE_"+algoType;
  TH2D* h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE->Sumw2();

  histoName = "ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE_"+algoType;
  TH2D* h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE->Sumw2();

  histoName = "response_vs_ptPhot_barrel_"+algoType;
  TH2D* h2_response_vs_ptPhot_barrel = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_barrel->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_barrel->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_barrel->Sumw2();

  histoName = "response_vs_ptPhot_eta14_20_PRE_"+algoType;
  TH2D* h2_response_vs_ptPhot_eta14_20_PRE = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{#gamma}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_eta14_20_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_eta14_20_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_eta14_20_PRE->Sumw2();

  histoName = "response_vs_ptPhot_eta14_20_"+algoType;
  TH2D* h2_response_vs_ptPhot_eta14_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{#gamma}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_eta14_20->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_eta14_20->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_eta14_20->Sumw2();

  histoName = "response_vs_ptPhot_endcap_"+algoType;
  TH2D* h2_response_vs_ptPhot_endcap = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{#gamma}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_endcap->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_endcap->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_endcap->Sumw2();

  histoName = "response_vs_ptPhot_eta0_3_"+algoType;
  TH2D* h2_response_vs_ptPhot_eta0_3 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{#gamma}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot_eta0_3->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_response_vs_ptPhot_eta0_3->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_response_vs_ptPhot_eta0_3->Sumw2();

  histoName = "responseGEN_vs_ptGen_barrel_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_barrel= new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_barrel->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_barrel->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_barrel->Sumw2();

  histoName = "responseGEN_vs_ptGen_eta14_20_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_eta14_20= new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_eta14_20->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_eta14_20->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_eta14_20->Sumw2();

  histoName = "responseGEN_vs_ptGen_endcap_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_endcap= new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_endcap->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_endcap->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_endcap->Sumw2();

  histoName = "responseGEN_vs_ptGen_eta0_3_"+algoType;
  TH2D* h2_responseGEN_vs_ptGen_eta0_3= new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_responseGEN_vs_ptGen_eta0_3->SetXTitle("p_{T}^{gen} [GeV/c]");
  h2_responseGEN_vs_ptGen_eta0_3->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");
  h2_responseGEN_vs_ptGen_eta0_3->Sumw2();


  histoName = "response_vs_etaPhot_PRE_"+algoType;
  TProfile* hp_response_vs_etaPhot_PRE = new TProfile(histoName.c_str(), "Jet Response vs. True Jet #eta", 50, -5.5, 5.5, 0., 2.);
  hp_response_vs_etaPhot_PRE->SetXTitle("#eta^{#gamma}");
  hp_response_vs_etaPhot_PRE->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");

  histoName = "response_vs_etaPhot_"+algoType;
  TProfile* hp_response_vs_etaPhot= new TProfile(histoName.c_str(), "Jet Response vs. True Jet #eta", 50, -5.5, 5.5, 0., 2.);
  hp_response_vs_etaPhot->SetXTitle("#eta^{#gamma}");
  hp_response_vs_etaPhot->SetYTitle("p_{T}^{jet}/p_{T}^{#gamma}");


  histoName = "ptGenMean_vs_pt_barrel_PRE_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_barrel_PRE= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel_PRE->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_PRE->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel_PRE->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_PRE_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_endcap_PRE= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap_PRE->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_PRE->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap_PRE->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_PRE_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta14_20_PRE= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20_PRE->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_PRE->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20_PRE->Sumw2();

  histoName = "ptGenMean_vs_pt_barrel_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_barrel= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_barrel->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_barrel->Sumw2();

  histoName = "ptGenMean_vs_pt_eta14_20_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta14_20= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta14_20->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_eta14_20->Sumw2();

  histoName = "ptGenMean_vs_pt_endcap_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_endcap= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_endcap->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_endcap->Sumw2();

  histoName = "ptGenMean_vs_pt_eta0_3_"+algoType;
  TH2D* h2_ptGenMean_vs_pt_eta0_3= new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptGenMean_vs_pt_eta0_3->SetXTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_3->SetYTitle("p_{T}^{#gen} [GeV/c]");
  h2_ptGenMean_vs_pt_eta0_3->Sumw2();

  histoName = "ptPhotMean_vs_pt_barrel_PRE_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_barrel_PRE = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_barrel_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_barrel_PRE->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_barrel_PRE->Sumw2();

  histoName = "ptPhotMean_vs_pt_endcap_PRE_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_endcap_PRE = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_endcap_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_endcap_PRE->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_endcap_PRE->Sumw2();

  histoName = "ptPhotMean_vs_pt_barrel_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_barrel = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_barrel->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_barrel->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_barrel->Sumw2();

  histoName = "ptPhotMean_vs_pt_eta14_20_PRE_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_eta14_20_PRE = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_eta14_20_PRE->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta14_20_PRE->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta14_20_PRE->Sumw2();

  histoName = "ptPhotMean_vs_pt_eta14_20_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_eta14_20 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_eta14_20->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta14_20->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta14_20->Sumw2();

  histoName = "ptPhotMean_vs_pt_endcap_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_endcap = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_endcap->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_endcap->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_endcap->Sumw2();

  histoName = "ptPhotMean_vs_pt_eta0_3_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_eta0_3 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_eta0_3->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta0_3->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta0_3->Sumw2();

  histoName = "ptPhotMean_vs_pt_eta155_3_"+algoType;
  TH2D* h2_ptPhotMean_vs_pt_eta155_3 = new TH2D(histoName.c_str(), "", nBins-1, Lower, 2000, 0., 2000.);
  h2_ptPhotMean_vs_pt_eta155_3->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta155_3->SetYTitle("p_{T}^{#gamma} [GeV/c]");
  h2_ptPhotMean_vs_pt_eta155_3->Sumw2();

 

  Float_t minPhotPt = 15.;
  Float_t max2ndJetPt = 10.;

  Int_t nTotalEvents=0;
  Float_t totalEvents=0.;
  Int_t nGoodEvents=0;
  Float_t goodEvents=0.;
 
  Int_t nEvents_ptPhot = 0;
  Int_t nEvents_ptPhot_etaPhot = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal_nTrk = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd = 0;
  Int_t nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd_phi = 0;

  Int_t nEvents_illegalWeight=0;

  int nEntries = tree->GetEntries();
//nEntries = 100000;


  //Masta Loop

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( NOWEIGHTS_ ) eventWeight=1.;

    if( eventWeight<=0. ) {
      eventWeight = 1.;
      nEvents_illegalWeight+=1;
    }

    Float_t Rch = eTracksReco;


    Float_t deltaEta_recoGen = etaJetReco - etaJetGen;
    Float_t deltaPhi_recoGen = phiJetReco - phiJetGen;
    Float_t deltaR_recoGen = sqrt( deltaEta_recoGen*deltaEta_recoGen + deltaPhi_recoGen*deltaPhi_recoGen );

    Float_t responseGEN = ptJetReco/ptJetGen;
    Float_t RchGEN = eTracksGen;

    if( deltaR_recoGen<0.1 ) {

      if( fabs(etaJetGen)<1.4 ) {
        h2_ptGenMean_vs_pt_barrel_PRE->Fill(ptJetGen, ptJetGen, eventWeight);
        h2_responseGEN_vs_ptGen_barrel_PRE->Fill(ptJetGen, responseGEN, eventWeight);
      } else if( fabs(etaJetGen)<3. ) {
        h2_ptGenMean_vs_pt_endcap_PRE->Fill(ptJetGen, ptJetGen, eventWeight);
        h2_responseGEN_vs_ptGen_endcap_PRE->Fill(ptJetGen, responseGEN, eventWeight);
      }

      if( (fabs(etaJetGen)>1.4)&&(fabs(etaJetGen)<2.) ) {
        h2_ptGenMean_vs_pt_eta14_20_PRE->Fill(ptJetGen, ptJetGen, eventWeight);
        h2_responseGEN_vs_ptGen_eta14_20_PRE->Fill(ptJetGen, responseGEN, eventWeight);
      }

      Int_t iBinGEN=-1;
      for(int i=0; i<(nBins-1); ++i) {
        if( (ptJetGen<=Lower[i+1])&&(ptJetGen>Lower[i]) )
          iBinGEN=i;
      }
      if( iBinGEN==-1 ) continue;
    
      if( fabs(etaJetGen)<1.4 ) {
        vec_h2_responseGEN_vs_Rch_barrel[iBinGEN]->Fill(Rch, responseGEN, eventWeight);
        vec_h2_Rch_barrelGEN[iBinGEN]->Fill(Rch, Rch, eventWeight);
        vec_h2_responseGEN_vs_RchGEN_barrel[iBinGEN]->Fill(RchGEN, responseGEN, eventWeight);
        vec_h2_RchGEN_barrelGEN[iBinGEN]->Fill(RchGEN, RchGEN, eventWeight);
        vec_h1_RchGEN_barrelGEN[iBinGEN]->Fill(RchGEN, eventWeight);
        h2_ptGenMean_vs_pt_barrel->Fill(ptJetGen, ptJetGen, eventWeight);
      }
    
      if( (fabs(etaJetGen)<2.) && (fabs(etaJetGen)>1.4) ) {
        vec_h2_responseGEN_vs_Rch_eta14_20[iBinGEN]->Fill(Rch, responseGEN, eventWeight);
        vec_h2_Rch_eta14_20GEN[iBinGEN]->Fill(Rch, Rch, eventWeight);
        h2_ptGenMean_vs_pt_eta14_20->Fill(ptJetGen, ptJetGen, eventWeight);
      }

      if( (fabs(etaJetGen)<3.) && (fabs(etaJetGen)>1.4) ) {
        vec_h2_responseGEN_vs_Rch_endcap[iBinGEN]->Fill(Rch, responseGEN, eventWeight);
        vec_h2_Rch_endcapGEN[iBinGEN]->Fill(Rch, Rch, eventWeight);
        vec_h2_responseGEN_vs_RchGEN_endcap[iBinGEN]->Fill(RchGEN, responseGEN, eventWeight);
        vec_h2_RchGEN_endcapGEN[iBinGEN]->Fill(RchGEN, RchGEN, eventWeight);
        vec_h1_RchGEN_endcapGEN[iBinGEN]->Fill(RchGEN, eventWeight);
        h2_ptGenMean_vs_pt_endcap->Fill(ptJetGen, ptJetGen, eventWeight);
      }

    } //if deltaR

    totalEvents += eventWeight;
    nTotalEvents += 1;

    if( ptPhotReco<minPhotPt ) continue;
    nEvents_ptPhot += 1;

    if( fabs(etaPhotReco)>1.3 ) continue;
    nEvents_ptPhot_etaPhot += 1;


    //photonID
    if( hcalIsoPhotReco >= 0.053) continue;
    nEvents_ptPhot_etaPhot_hcal += 1;
    if( ecalIsoPhotReco >= 0.05) continue;
    nEvents_ptPhot_etaPhot_hcal_ecal += 1;
    if( nTrkIsoPhotReco >= 3) continue;
    nEvents_ptPhot_etaPhot_hcal_ecal_nTrk += 1;
    if( ptTrkIsoPhotReco >= 0.1) continue;
    nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk += 1;
    if( clusterMajPhotReco >= 0.3) continue;
    nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj += 1;
    if( clusterMinPhotReco >= 0.25) continue;
    nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min += 1;

    //2nd jet ID
    Float_t deltaEta_2ndJet = etaPhotReco - eta2ndJetReco;
    Float_t deltaPhi_2ndJet = phiPhotReco - phi2ndJetReco;
    Float_t deltaR_2ndJet = sqrt( deltaEta_2ndJet*deltaEta_2ndJet + deltaPhi_2ndJet*deltaPhi_2ndJet );
    if( (pt2ndJetReco>0.1*ptPhotReco)&&(pt2ndJetReco>max2ndJetPt)) continue;
    h1_deltaR_phot_2ndJet->Fill(deltaR_2ndJet, eventWeight);
    nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd += 1;

    //leading jet and photon back2back in transverse plane
    Float_t deltaPhi_jet = phiPhotReco - phiJetReco;
    Float_t pi = 3.14159;
    if( fabs(deltaPhi_jet) < (pi - 0.2) ) continue;
    nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd_phi += 1;
    

    goodEvents+=eventWeight;
    nGoodEvents+=1;

    Float_t response = ptJetReco/ptPhotReco;


    hp_response_vs_etaPhot_PRE->Fill(etaJetReco, response, eventWeight);



    if( fabs(etaJetReco)<1.4 ) {
      h2_ptPhotMean_vs_pt_barrel_PRE->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_barrel_PRE->Fill(ptPhotReco, response, eventWeight);
      h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE->Fill(ptPhotReco, ptJetGen/ptPhotReco, eventWeight);
    } else if( fabs(etaJetReco)<3. ) {
      h2_ptPhotMean_vs_pt_endcap_PRE->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_endcap_PRE->Fill(ptPhotReco, response, eventWeight);
      h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE->Fill(ptPhotReco, ptJetGen/ptPhotReco, eventWeight);
    }

    if( (fabs(etaJetReco)<2.) && (fabs(etaJetReco)>1.4) ) {
      h2_ptPhotMean_vs_pt_eta14_20_PRE->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_eta14_20_PRE->Fill(ptPhotReco, response, eventWeight);
    }

    Int_t iBin=-1;
    for(int i=0; i<(nBins-1); ++i) {
      if( (ptPhotReco<=Lower[i+1])&&(ptPhotReco>Lower[i]) )
        iBin=i;
    }
    if( iBin==-1 ) continue;



    if( fabs(etaJetReco)<1.4 ) {
      vec_h2_response_vs_Rch_barrel[iBin]->Fill(Rch, response, eventWeight);
      vec_h2_Rch_barrel[iBin]->Fill(Rch, Rch, eventWeight);
      h2_ptPhotMean_vs_pt_barrel->Fill(ptPhotReco, ptPhotReco, eventWeight);
    }

    if( (fabs(etaJetReco)<2.) && (fabs(etaJetReco)>1.4) ) {
      vec_h2_response_vs_Rch_eta14_20[iBin]->Fill(Rch, response, eventWeight);
      vec_h2_Rch_eta14_20[iBin]->Fill(Rch, Rch, eventWeight);
      h2_ptPhotMean_vs_pt_eta14_20->Fill(ptPhotReco, ptPhotReco, eventWeight);
    }

    if( (fabs(etaJetReco)<3.) && (fabs(etaJetReco)>1.4) ) {
      vec_h2_response_vs_Rch_endcap[iBin]->Fill(Rch, response, eventWeight);
      vec_h2_Rch_endcap[iBin]->Fill(Rch, Rch, eventWeight);
      h2_ptPhotMean_vs_pt_endcap->Fill(ptPhotReco, ptPhotReco, eventWeight);
    }

    if( (fabs(etaJetReco)<3.) && (fabs(etaJetReco)>1.55) ) {
      vec_h2_response_vs_Rch_eta155_3[iBin]->Fill(Rch, response, eventWeight);
      vec_h2_Rch_eta155_3[iBin]->Fill(Rch, Rch, eventWeight);
      h2_ptPhotMean_vs_pt_eta155_3->Fill(ptPhotReco, ptPhotReco, eventWeight);
    }




  } //for entries

  if( nEvents_illegalWeight > 0 ) {
    std::cout << std::endl << "---->>> WARNING!!!!!!" << std::endl;
    std::cout << std::endl << "---->>> WARNING!!!!!!" << std::endl;
    std::cout << "---->>> " << nEvents_illegalWeight << " events had a weight <= 0!!!" << std::endl;
  }

std::cout << std::endl << std::endl << "@@@@@ N TOTAL EVENTS: " << nTotalEvents << std::endl;
std::cout << std::endl << std::endl << "@@@@@ N GOOD EVENTS: " << nGoodEvents << std::endl;

  std::cout << "nEvents_ptPhot: " <<                                              nEvents_ptPhot << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot: " <<                                      nEvents_ptPhot_etaPhot << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal: " <<                                 nEvents_ptPhot_etaPhot_hcal << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal: " <<                            nEvents_ptPhot_etaPhot_hcal_ecal << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal_nTrk: " <<                       nEvents_ptPhot_etaPhot_hcal_ecal_nTrk << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk: " <<                 nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj: " <<             nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min: " <<         nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd: " <<     nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd << std::endl;
  std::cout << "nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd_phi: " << nEvents_ptPhot_etaPhot_hcal_ecal_nTrk_ptTrk_maj_min_2nd_phi << std::endl;

  fitTools::fitDistribution_TGraph(h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE, h2_ptPhotMean_vs_pt_barrel, "response", "pt", "barrel", "GENoverPHOT_PRE", algoType, fileName, "barrel_GENoverPHOT_PRE");
  fitTools::fitDistribution_TGraph(h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE, h2_ptPhotMean_vs_pt_endcap,  "response","pt", "endcap", "GENoverPHOT_PRE", algoType, fileName, "endcap_GENoverPHOT_PRE");

  fitTools::fitDistribution_TGraph(h2_ptJetGenOverPtPhot_vs_ptPhot_barrel_PRE, h2_ptGenMean_vs_pt_barrel,  "response","pt", "barrel", "GENoverPHOT_GEN", algoType, fileName, "barrel_GENoverPHOT_GEN");
  fitTools::fitDistribution_TGraph(h2_ptJetGenOverPtPhot_vs_ptPhot_endcap_PRE, h2_ptGenMean_vs_pt_endcap,  "response","pt", "endcap", "GENoverPHOT_GEN", algoType, fileName, "endcap_GENoverPHOT_GEN");

  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_barrel_PRE, h2_ptPhotMean_vs_pt_barrel, "response", "pt", "barrel", "PRE", algoType, fileName, "Rch_barrel_RECO_PRE");
  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_endcap_PRE, h2_ptPhotMean_vs_pt_endcap, "response", "pt", "endcap", "PRE", algoType, fileName, "Rch_endcap_RECO_PRE");
  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_eta14_20_PRE, h2_ptPhotMean_vs_pt_eta14_20, "response", "pt", "eta14_20", "PRE", algoType, fileName, "Rch_eta14_20_RECO_PRE");

  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_barrel_PRE, h2_ptGenMean_vs_pt_barrel_PRE, "response", "pt", "barrel", "GEN_PRE", algoType, fileName, "Rch_barrel_GEN_PRE");
  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_endcap_PRE, h2_ptGenMean_vs_pt_endcap_PRE, "response", "pt", "endcap", "GEN_PRE", algoType, fileName, "Rch_endcap_GEN_PRE");
  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_eta14_20_PRE, h2_ptGenMean_vs_pt_eta14_20_PRE, "response", "pt", "eta14_20", "GEN_PRE", algoType, fileName, "Rch_eta14_20_GEN_PRE");

  std::vector<FitResults> fitResults_barrel   = fitTheHellOutOfThem(vec_h2_response_vs_Rch_barrel, vec_h2_Rch_barrel, "barrel", dataset);
  std::vector<FitResults> fitResults_eta14_20 = fitTheHellOutOfThem(vec_h2_response_vs_Rch_eta14_20, vec_h2_Rch_eta14_20, "eta14_20", dataset);
  std::vector<FitResults> fitResults_endcap   = fitTheHellOutOfThem(vec_h2_response_vs_Rch_endcap, vec_h2_Rch_endcap, "endcap", dataset);
  std::vector<FitResults> fitResults_eta155_3 = fitTheHellOutOfThem(vec_h2_response_vs_Rch_eta155_3, vec_h2_Rch_eta155_3, "eta155_3", dataset);

  std::vector<FitResults> fitResults_barrelGEN = fitTheHellOutOfThem(vec_h2_responseGEN_vs_Rch_barrel, vec_h2_Rch_barrelGEN, "barrelGEN", dataset);
  std::vector<FitResults> fitResults_endcapGEN = fitTheHellOutOfThem(vec_h2_responseGEN_vs_Rch_endcap, vec_h2_Rch_endcapGEN, "endcapGEN", dataset);
  std::vector<FitResults> fitResults_eta14_20GEN = fitTheHellOutOfThem(vec_h2_responseGEN_vs_Rch_eta14_20, vec_h2_Rch_eta14_20GEN, "eta14_20GEN", dataset);

  std::vector<FitResults> fitResultsGEN_barrelGEN = fitTheHellOutOfThem(vec_h2_responseGEN_vs_RchGEN_barrel, vec_h2_RchGEN_barrelGEN, "barrelGENGEN", dataset);
  std::vector<FitResults> fitResultsGEN_endcapGEN = fitTheHellOutOfThem(vec_h2_responseGEN_vs_RchGEN_endcap, vec_h2_RchGEN_endcapGEN, "endcapGENGEN", dataset);


  TFile* outfile = TFile::Open(fileName.c_str(), "update");
  outfile->cd();
  h1_deltaR_phot_2ndJet->Write();
  for( int i = 0; i<nBins; ++i) {
    vec_h2_response_vs_Rch_barrel[i]->Write();
    vec_h2_response_vs_Rch_eta14_20[i]->Write();
    vec_h2_response_vs_Rch_endcap[i]->Write();
    vec_h2_response_vs_Rch_eta155_3[i]->Write();

    vec_h2_responseGEN_vs_Rch_barrel[i]->Write();
    vec_h2_responseGEN_vs_Rch_endcap[i]->Write();
    vec_h2_responseGEN_vs_Rch_eta14_20[i]->Write();

    vec_h2_Rch_barrel[i]->Write();
    vec_h2_Rch_eta14_20[i]->Write();
    vec_h2_Rch_endcap[i]->Write();
    vec_h2_Rch_eta155_3[i]->Write();

    vec_h1_RchGEN_barrelGEN[i]->Write();
    vec_h1_RchGEN_endcapGEN[i]->Write();
    
  }
  outfile->Close();
  outfile->Write();

//return;

  //re-loop to apply corrections:

  std::cout << std::endl << "Applying corrections." << std::endl;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( NOWEIGHTS_ ) eventWeight=1.;

    Float_t Rch = eTracksReco;


    Float_t deltaEta_recoGen = etaJetReco - etaJetGen;
    Float_t deltaPhi_recoGen = phiJetReco - phiJetGen;
    Float_t deltaR_recoGen = sqrt( deltaEta_recoGen*deltaEta_recoGen + deltaPhi_recoGen*deltaPhi_recoGen );

    if( deltaR_recoGen < 0.1 ) {

      Int_t iBinGEN=-1;
      for(int i=0; i<(nBins-1); ++i) {
        if( (ptJetGen<=Lower[i+1])&&(ptJetGen>Lower[i]) )
          iBinGEN=i;
      }
      if( iBinGEN==-1 ) continue;


      Float_t etaMaxCorrect = 3.;

      FitResults paramsGEN;
      if( fabs(etaJetGen)<1.4 )
        paramsGEN = fitResults_barrelGEN[iBinGEN];
      else if( fabs(etaJetGen)<etaMaxCorrect )
        paramsGEN = fitResults_eta14_20GEN[iBinGEN];
      
      Float_t responseGEN = ptJetReco/ptJetGen;

      Float_t correctionGEN;
      if( fabs(etaJetGen)<etaMaxCorrect ) {
        correctionGEN = paramsGEN.p0 + paramsGEN.p1*Rch + paramsGEN.p2*Rch*Rch + paramsGEN.p3*Rch*Rch*Rch;
        correctionGEN += paramsGEN.p4*Rch*Rch*Rch*Rch + paramsGEN.p5*Rch*Rch*Rch*Rch*Rch;
      } else {
        correctionGEN = 1.;
      }

      if( correctionGEN > 0. )
        responseGEN /= correctionGEN;

      if( fabs(etaJetGen)<1.4 ) {
        h2_responseGEN_vs_ptGen_barrel->Fill(ptJetGen, responseGEN, eventWeight);
      }
    
      if( (fabs(etaJetGen)<2.)&&(fabs(etaJetGen)>1.4) ) {
        h2_responseGEN_vs_ptGen_eta14_20->Fill(ptJetGen, responseGEN, eventWeight);
      }
    
      if( (fabs(etaJetGen)<3.)&&(fabs(etaJetGen)>1.4) ) {
        h2_responseGEN_vs_ptGen_endcap->Fill(ptJetGen, responseGEN, eventWeight);
      }
    
      if( fabs(etaJetGen)<3. ) {
        h2_responseGEN_vs_ptGen_eta0_3->Fill(ptJetGen, responseGEN, eventWeight);
      }

    } //if deltaR

    if( ptPhotReco<minPhotPt ) continue;
    if( fabs(etaPhotReco)>1.3 ) continue;

    //photonID
    if( hcalIsoPhotReco >= 0.053) continue;
    if( ecalIsoPhotReco >= 0.05) continue;
    if( nTrkIsoPhotReco >= 3) continue;
    if( ptTrkIsoPhotReco >= 0.1) continue;
    if( clusterMajPhotReco >= 0.3) continue;
    if( clusterMinPhotReco >= 0.25) continue;

    //2nd jet ID
    if( (pt2ndJetReco>0.1*ptPhotReco)&&(pt2ndJetReco>max2ndJetPt)) continue;

    //leading jet and photon back2back in transverse plane
    Float_t deltaPhi_phot_jet = phiPhotReco - phiJetReco;
    Float_t pi = 3.14159;
    if( fabs(deltaPhi_phot_jet) < pi - 0.2 ) continue;
    

    Float_t response = ptJetReco/ptPhotReco;

    Int_t iBin=-1;
    for(int i=0; i<(nBins-1); ++i) {
      if( (ptPhotReco<=Lower[i+1])&&(ptPhotReco>Lower[i]) )
        iBin=i;
    }
    if( iBin==-1 ) continue;

    FitResults params;
    if( fabs(etaJetReco)<1.4 )
      params = fitResults_barrel[iBin];
    else if( fabs(etaJetReco)<3. )
      params = fitResults_endcap[iBin];



    Float_t correction;
    if( fabs(etaJetReco)<3. ) {
      correction = params.p0 + params.p1*Rch + params.p2*Rch*Rch + params.p3*Rch*Rch*Rch;
      correction += params.p4*Rch*Rch*Rch*Rch + params.p5*Rch*Rch*Rch*Rch*Rch;
    } else {
      correction = 1.;
    }

    if( correction > 0. )
      response /= correction;

    hp_response_vs_etaPhot->Fill( etaJetReco, response, eventWeight );



    if( fabs(etaJetReco)<1.4 ) {
      //h2_ptPhotMean_vs_pt_barrel->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_barrel->Fill(ptPhotReco, response, eventWeight);
    }

    if( (fabs(etaJetReco)<2.)&&(fabs(etaJetReco)>1.4) ) {
      //h2_ptPhotMean_vs_pt_eta14_20->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_eta14_20->Fill(ptPhotReco, response, eventWeight);
    }

    if( (fabs(etaJetReco)<3.)&&(fabs(etaJetReco)>1.4) ) {
      //h2_ptPhotMean_vs_pt_endcap->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_endcap->Fill(ptPhotReco, response, eventWeight);
    }

    if( fabs(etaJetReco)<3. ) {
      //h2_ptPhotMean_vs_pt_eta0_3->Fill(ptPhotReco, ptPhotReco, eventWeight);
      h2_response_vs_ptPhot_eta0_3->Fill(ptPhotReco, response, eventWeight);
    }


      
  } //for entries


  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_barrel, h2_ptGenMean_vs_pt_barrel, "response","pt", "barrel", "GEN", algoType, fileName, "Rch_barrel_GEN_POST");
  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_eta14_20, h2_ptGenMean_vs_pt_eta14_20, "response", "pt", "eta14_20", "GEN", algoType, fileName, "Rch_eta14_20_GEN_POST");
  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_endcap, h2_ptGenMean_vs_pt_endcap, "response", "pt", "endcap", "GEN", algoType, fileName, "Rch_endcap_GEN_POST");
  fitTools::fitDistribution_TGraph(h2_responseGEN_vs_ptGen_eta0_3, h2_ptGenMean_vs_pt_eta0_3, "response", "pt", "eta0_3", "GEN", algoType, fileName, "Rch_eta0_3_GEN_POST");

  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_barrel, h2_ptPhotMean_vs_pt_barrel, "response", "pt", "barrel", "", algoType, fileName, "Rch_barrel_RECO_POST");
  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_eta14_20, h2_ptPhotMean_vs_pt_eta14_20, "response", "pt", "eta14_20", "", algoType, fileName, "Rch_eta14_20_RECO_POST");
  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_endcap, h2_ptPhotMean_vs_pt_endcap, "response", "pt", "endcap", "", algoType, fileName, "Rch_endcap_RECO_POST");
  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot_eta0_3, h2_ptPhotMean_vs_pt_eta0_3, "response", "pt", "eta0_3", "", algoType, fileName, "Rch_eta0_3_RECO_POST");


  outfile = TFile::Open(fileName.c_str(), "update");
  outfile->cd();
  hp_response_vs_etaPhot_PRE->Write();
  hp_response_vs_etaPhot->Write();
  h2_ptPhotMean_vs_pt_barrel->Write();
  h2_ptPhotMean_vs_pt_endcap->Write();
  h2_ptPhotMean_vs_pt_eta14_20->Write();
  h2_ptPhotMean_vs_pt_eta0_3->Write();

  h2_ptGenMean_vs_pt_barrel->Write();
  h2_ptGenMean_vs_pt_endcap->Write();
  h2_ptGenMean_vs_pt_eta14_20->Write();
  h2_ptGenMean_vs_pt_eta0_3->Write();

  outfile->Close();



} //fillFile



std::vector<FitResults> fitTheHellOutOfThem(std::vector<TH2D*> vec_h2, std::vector<TH2D*> vec_h2_Rch, std::string name, const std::string& dataset) {


  std::string dirName = "Projections_RchCalib_" + dataset;
  if( NOWEIGHTS_ )
    dirName = dirName + "_NOWEIGHTS";
  else
    dirName = dirName + "_CHI";
  dirName = dirName + "/";
  std::string fileName = dirName + name + ".root";

  if( name!="" ) std::cout << "Saving in : " << fileName << std::endl;

  TFile* projectionFile;
  if( name!= "" ) {
    projectionFile = TFile::Open(fileName.c_str(), "RECREATE");
    projectionFile->cd();
  }

  std::vector<FitResults> returnData;

  Int_t nBins = vec_h2.size();

  //for(std::vector<TH2D*>::const_iterator iH2 = vec_h2.begin(); iH2!=vec_h2.end(); ++iH2) {
  for(int iBin=0; iBin<nBins; ++iBin) {

    //char histName[30];
    //sprintf(histName, "fit_%d", i);
    //TH1F* h1_fit = new TH1F(histName, "Response vs. R_{ch}", _nBins_Rch, 0., 1.);
    //h1_fit->SetXTitle("R_{ch}");
    //h1_fit->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

    char dirName[30];
    sprintf(dirName, "ptbin_%d", iBin);
    projectionFile->mkdir(dirName);
    projectionFile->cd(dirName);

    Float_t x[_nBins_Rch];
    Float_t x_err[_nBins_Rch];
    Float_t response[_nBins_Rch];
    Float_t response_err[_nBins_Rch];
    Float_t resolution[_nBins_Rch];
    Float_t resolution_err[_nBins_Rch];

    for(int iBin_Rch=0; iBin_Rch<_nBins_Rch; ++iBin_Rch) {

      //first: get the x value:
      char proj_Rch_name[50];
      sprintf(proj_Rch_name, "Rch_Rchbin_%d", iBin_Rch+1);
      TH1D* h1_projection_Rch = vec_h2_Rch[iBin]->ProjectionY(proj_Rch_name, iBin_Rch+1, iBin_Rch+1);

      x[iBin_Rch] = h1_projection_Rch->GetMean();
      x_err[iBin_Rch] = h1_projection_Rch->GetRMS();

      //second: get the y value:
      char proj_name[50];
      sprintf(proj_name, "response_Rchbin_%d", iBin_Rch+1);
      TH1D* h1_projection = vec_h2[iBin]->ProjectionY(proj_name, iBin_Rch+1, iBin_Rch+1);

      Float_t mean = h1_projection->GetMean();
      Float_t rms = h1_projection->GetRMS();
      Float_t rms_err = h1_projection->GetMeanError();
      Float_t mean_err = h1_projection->GetRMSError();
      //Float_t rms_err = rms / sqrt( (Float_t) h1_projection->GetEntries() );
      //Float_t mean_err = rms_err;
  
      Int_t entries = h1_projection->GetEntries();

      if( entries>0 ) {
        if( !DOFITS_ ) {
          //h1_fit->SetBinContent(iBin_Rch, mean);
          response[iBin_Rch] = mean;
          response_err[iBin_Rch] = mean_err;
          //h1_fit->SetBinError(iBin_Rch, rms);

          resolution[iBin_Rch] = rms / mean;
          resolution_err[iBin_Rch] = rms_err;

        } else {
      
          TF1* gaussian_LL = new TF1("gaussian_LL", "gaus");
          gaussian_LL->SetLineColor(kRed);
          fitTools::fitProjection(h1_projection, gaussian_LL, 2., "RQLL");
      
          TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
          gaussian_chi->SetLineColor(kBlue);
          fitTools::fitProjection(h1_projection, gaussian_chi, 2., "RQ+");
          
          Float_t mu;
          Float_t sigma;
          Float_t sigma_err;

          if( NOWEIGHTS_ ) {
            mu         = gaussian_LL->GetParameter(1);
            sigma      = gaussian_LL->GetParameter(2);
            sigma_err  = gaussian_LL->GetParError(2);
          } else {
            mu         = gaussian_chi->GetParameter(1);
            sigma      = gaussian_chi->GetParameter(2);
            sigma_err  = gaussian_chi->GetParError(2);
          }

          response[iBin_Rch] = mu;
          response_err[iBin_Rch] = sigma/sqrt(entries);
      
          if( NOWEIGHTS_ ) {
            resolution[iBin_Rch] = (mean!=0.) ? rms/mean : -1.;
            resolution_err[iBin_Rch] = (mean!=0.) ? sqrt( rms_err*rms_err/(mean*mean) + mean_err*mean_err*rms*rms/(mean*mean*mean*mean) ) : 0.;
          } else {
            resolution[iBin_Rch] = sigma / mu;
            resolution_err[iBin_Rch] = sigma_err;
          }

        }
      } //if entries

      h1_projection->Write();
      //h1_projection_Rch->Write();

    } //for bins

    char graphName[30];
    sprintf(graphName, "fitResponse_%d", iBin);
    TGraphErrors* gr_response_vs_Rch = new TGraphErrors(_nBins_Rch, x, response, x_err, response_err);
    gr_response_vs_Rch->SetName(graphName);

    TF1* func = new TF1("func", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x", 0., 1.);
    func->SetParameter(0, response[1]);
    Float_t deltaY = response[5] - response[4];
    Float_t deltaX = 10.;
    func->SetParameter(1, deltaY/deltaX);
    func->SetParameter(2, 0.);
    func->SetParameter(3, 0.);
    func->SetParameter(4, 0.);
    func->SetParameter(5, 0.);

    gr_response_vs_Rch->Fit(func, "RQ");

    sprintf(graphName, "fitResolution_%d", iBin);
    TGraphErrors* gr_resolution_vs_Rch = new TGraphErrors(_nBins_Rch, x, resolution, x_err, resolution_err);
    gr_resolution_vs_Rch->SetName(graphName);

    if( name!="") {
      gr_response_vs_Rch->Write();
      gr_resolution_vs_Rch->Write();
    }

    FitResults fitResults;

    fitResults.p0 = func->GetParameter(0);
    fitResults.p1 = func->GetParameter(1);
    fitResults.p2 = func->GetParameter(2);
    fitResults.p3 = func->GetParameter(3);
    fitResults.p4 = func->GetParameter(4);
    fitResults.p5 = func->GetParameter(5);

    delete func;
    delete gr_response_vs_Rch;
    delete gr_resolution_vs_Rch;

    returnData.push_back(fitResults);

  } //for TH2Ds

  projectionFile->Write();
  projectionFile->Close();

  return returnData;

}


std::vector< TH2D* > bookResponseHistoVector(const std::string& name, Int_t nBins) {

  std::vector< TH2D* > vec_h2;

  for( int i=0; i<nBins; ++i) {
    char h2_name[70];
    sprintf(h2_name, "%s_bin%d", name.c_str(), (i+1));
    TH2D* h2 = new TH2D(h2_name, "", _nBins_Rch, 0., 1.001, 25, 0., 2.);
    h2->Sumw2();
    vec_h2.push_back(h2);
  }

  return vec_h2;

}


std::vector< TH2D* > bookRchTH2Vector(const std::string& name, Int_t nBins) {

  std::vector< TH2D* > vec_h2;

  for( int i=0; i<nBins; ++i) {
    char h2_name[70];
    sprintf(h2_name, "%s_bin%d", name.c_str(), (i+1));
    TH2D* h2 = new TH2D(h2_name, "", _nBins_Rch, 0., 1.001, 10*_nBins_Rch, 0., 1.001);
    h2->Sumw2();
    vec_h2.push_back(h2);
  }

  return vec_h2;

}


std::vector< TH1F* > bookRchTH1Vector(const std::string& name, Int_t nBins) {

  std::vector< TH1F* > vec_h1;

  for( int i=0; i<nBins; ++i) {
    char h1_name[70];
    sprintf(h1_name, "%s_bin%d", name.c_str(), (i+1));
    TH1F* h1 = new TH1F(h1_name, "", 100, 0., 1.1);
    h1->Sumw2();
    vec_h1.push_back(h1);
  }

  return vec_h1;

}
