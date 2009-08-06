#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>

#include "fitTools.C"


class Jet {
public:
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
  Float_t Rhf;
  Float_t Re;
  Float_t RchGen;
  Float_t RgammaGen;
  Float_t RnhGen;
  Float_t RmuGen;
  Float_t ReGen;
};



void fillFile(TFile* file, TTree* tree, std::string algoType, bool draw=false);


bool DEBUG_ = false;
bool DEBUG_PROVA_ = false;

bool LONGSCALE_ = false;


void finalize() {


  TChain* chain = new TChain("jetTree");

  int n_files = 49;
  if( DEBUG_ ) n_files = 1;

  if( DEBUG_PROVA_ ) {
    chain->Add("prova2ndLevel.root/jetTree");
  } else {
    for(int i=1; i<(n_files+1); ++i) {

      char treeName[70];
      sprintf(treeName, "2ndLevelOutputFiles/2ndLevelOutputFile_%d.root/jetTree", i);
      chain->Add(treeName);

    } //for
  } //else
  
  bool draw=true;

  std::string fileName;
  if( DEBUG_ )
    fileName = "provaOutputFile.root";
  else 
    fileName = "finalOutputFile.root";
  

  TFile* outfile = new TFile(fileName.c_str(), "RECREATE");
  outfile->cd();

  std::cout << std::endl << std::endl << "**** starting calo it cone 5" << std::endl;
  fillFile(outfile, chain, "caloItCone5", draw);
  std::cout << std::endl << std::endl << "**** starting calo SIS cone 5" << std::endl;
  fillFile(outfile, chain, "caloSISCone5");
  std::cout << std::endl << std::endl << "**** starting calo kt4" << std::endl;
  fillFile(outfile, chain, "calokt4");
  std::cout << std::endl << std::endl << "**** starting PF it cone 5" << std::endl;
  fillFile(outfile, chain, "PFItCone5", draw);

  outfile->Write();
  outfile->Close();
  delete outfile;
  outfile = 0;


}  //finalize




void fillFile(TFile* file, TTree* tree, std::string algoType, bool draw) {

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
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eReco_name.c_str(), eReco);
  Float_t ptTracksReco[40];
  std::string ptTracksReco_name = "ptTracksReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ptTracksReco_name.c_str(), ptTracksReco);
  Float_t ptPhotonsReco[40];
  std::string ptPhotonsReco_name = "ptPhotonsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ptPhotonsReco_name.c_str(), ptPhotonsReco);
  Float_t eTracksReco[40];
  std::string eTracksReco_name = "eTracksReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eTracksReco_name.c_str(), eTracksReco);
  Float_t ePhotonsReco[40];
  std::string ePhotonsReco_name = "ePhotonsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ePhotonsReco_name.c_str(), ePhotonsReco);
  Float_t eNeutralHadronsReco[40];
  std::string eNeutralHadronsReco_name = "eNeutralHadronsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eNeutralHadronsReco_name.c_str(), eNeutralHadronsReco);
  Float_t eMuonsReco[40];
  std::string eMuonsReco_name = "eMuonsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eMuonsReco_name.c_str(), eMuonsReco);
  Float_t eElectronsReco[40];
  std::string eElectronsReco_name = "eElectronsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eElectronsReco_name.c_str(), eElectronsReco);
  Float_t eHFhadronsReco[40];
  std::string eHFhadronsReco_name = "eHFhadronsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eHFhadronsReco_name.c_str(), eHFhadronsReco);

  Float_t eGen[40];
  std::string eGen_name = "eJetGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eGen_name.c_str(), eGen);
  Float_t eTracksGen[40];
  std::string eTracksGen_name = "eTracksGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eTracksGen_name.c_str(), eTracksGen);
  Float_t ePhotonsGen[40];
  std::string ePhotonsGen_name = "ePhotonsGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ePhotonsGen_name.c_str(), ePhotonsGen);
  Float_t eNeutralHadronsGen[40];
  std::string eNeutralHadronsGen_name = "eNeutralHadronsGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eNeutralHadronsGen_name.c_str(), eNeutralHadronsGen);
  Float_t eMuonsGen[40];
  std::string eMuonsGen_name = "eMuonsGen_"+algoType;
  if( algoType=="PFItCone5")
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
  int nBins;

  if( LONGSCALE_ )
    nBins = 23;
  else
    nBins = 20;


  Double_t Lower[nBins];

  for(int iLower=0; iLower<11; ++iLower)
    Lower[iLower] = iLower*10.;
  
  Lower[11] = 120.;
  Lower[12] = 140.;
  Lower[13] = 160.;
  Lower[14] = 180.;
  Lower[15] = 200.;
  Lower[16] = 250.;
  Lower[17] = 300.;
  Lower[18] = 400.;
  Lower[19] = 500.;

  if( LONGSCALE_ ) {
  Lower[20] = 700.;
    Lower[21] = 1000.;
    Lower[22] = 1500.;
  }



  // histograms booking for all algorithms:
  histoName = "deltaEta_vs_pt_barrel_"+algoType;
  TH2D* h2_deltaEta_vs_pt_barrel = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_barrel = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType;
  TH2D* h2_deltaEta_vs_pt_endcap = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_endcap = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_forward_"+algoType;
  TH2D* h2_deltaEta_vs_pt_forward = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaEta_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_forward->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_forward_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_forward = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_forward->SetYTitle("#Delta#phi");


  histoName = "response_vs_pt_barrel_"+algoType;
  TH2D* h2_response_vs_pt_barrel = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType;
  TH2D* h2_response_vs_pt_endcap = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_forward_"+algoType;
  TH2D* h2_response_vs_pt_forward = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_eta_"+algoType;
  TH2D* h2_response_vs_eta = new TH2D(histoName.c_str(), "Jet Response vs. True #eta", 100, -5., 5., 120, 0., 2.);
  h2_response_vs_eta->SetXTitle("#eta");
  h2_response_vs_eta->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");



  //additional histograms for PF only:
//  histoName = "response_vs_pt_barrel_onlyCh_"+algoType;
//  TH2D* h2_response_vs_pt_barrel_onlyCh = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower, 120, 0., 2.);
//  h2_response_vs_pt_barrel_onlyCh->SetXTitle("p_{T}^{GEN} [GeV/c]");
//  h2_response_vs_pt_barrel_onlyCh->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_Rch_"+algoType;
  TH2D* h2_response_vs_Rch = new TH2D(histoName.c_str(), "Jet Response vs. R_{ch}", 50, 0., 1., 120, 0., 2.);
  h2_response_vs_Rch->SetXTitle("R_{ch}");
  h2_response_vs_Rch->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_Rgamma_"+algoType;
  TH2D* h2_response_vs_Rgamma = new TH2D(histoName.c_str(), "Jet Response vs. R_{#gamma}", 50, 0., 1., 120, 0., 2.);
  h2_response_vs_Rgamma->SetXTitle("R_{#gamma}");
  h2_response_vs_Rgamma->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");


  //--------------------- Rch --------------------------------------------

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch10_30->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch10_30->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch50_70->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch50_70->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch90_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch90_100->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch0->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch10_30->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch10_30->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch50_70->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch50_70->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch90_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch90_100->SetYTitle("#Delta#phi");


  //---------------------- Rgamma
/*
  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma5_20 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rgamma5_20->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma5_20 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma5_20->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma30_60 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rgamma30_60->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma30_60 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma30_60->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma60_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rgamma60_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma60_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma60_100->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma0->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma5_20 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma5_20->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma5_20 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma5_20->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma30_60 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma30_60->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma30_60 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma30_60->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma60_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma60_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma60_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma60_100->SetYTitle("#Delta#phi");
*/


  histoName = "response_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");



  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH2D* h2_response_vs_pt_barrel_Rgamma0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, R_{gamma}=0%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH2D* h2_response_vs_pt_barrel_Rgamma5_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 5<R_{gamma}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH2D* h2_response_vs_pt_barrel_Rgamma30_60 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH2D* h2_response_vs_pt_barrel_Rgamma60_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH2D* h2_response_vs_pt_endcap_Rgamma0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, R_{gamma}=0%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH2D* h2_response_vs_pt_endcap_Rgamma5_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 5<R_{gamma}<20%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH2D* h2_response_vs_pt_endcap_Rgamma30_60 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH2D* h2_response_vs_pt_endcap_Rgamma60_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");



  histoName = "Rch_"+algoType;
  TH1F* h1_Rch = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1.1);
  h1_Rch->SetXTitle("R_{ch}");
  h1_Rch->SetYTitle("entries");

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

  histoName = "Rch_vs_pt_"+algoType;
  TProfile* hp_Rch_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rch_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt->SetYTitle("entries");

  histoName = "Rch_vs_pt_barrel_"+algoType;
  TProfile* hp_Rch_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rch_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt_barrel->SetYTitle("entries");

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

  histoName = "Rgamma_vs_eta_"+algoType;
  TProfile* hp_Rgamma_vs_eta = new TProfile(histoName.c_str(), "R_{#gamma} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rgamma_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rgamma_vs_eta->SetYTitle("entries");

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

  histoName = "Rnh_vs_eta_"+algoType;
  TProfile* hp_Rnh_vs_eta = new TProfile(histoName.c_str(), "R_{nh} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rnh_vs_eta->SetYTitle("entries");

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
  TProfile* hp_Re_vs_eta = new TProfile(histoName.c_str(), "R_{e} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Re_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Re_vs_eta->SetYTitle("entries");

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
  TProfile* hp_Rmu_vs_eta = new TProfile(histoName.c_str(), "R_{#mu} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rmu_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rmu_vs_eta->SetYTitle("entries");

  histoName = "Rmu_vs_Rch_"+algoType;
  TProfile* hp_Rmu_vs_Rch = new TProfile(histoName.c_str(), "R_{#mu} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rmu_vs_Rch->SetXTitle("R_{ch}");
  hp_Rmu_vs_Rch->SetYTitle("R_{#mu}");


  histoName = "Rhf_"+algoType;
  TH1F* h1_Rhf = new TH1F(histoName.c_str(), "Jet HF Hadron Energy Fraction", 100, 0., 1.1);
  h1_Rhf->SetXTitle("R_{HF}");
  h1_Rhf->SetYTitle("entries");

  histoName = "Rhf_vs_pt_"+algoType;
  TProfile* hp_Rhf_vs_pt = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rhf_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhf_vs_pt->SetYTitle("entries");

  histoName = "Rhf_vs_pt_endcap_"+algoType;
  TProfile* hp_Rhf_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_Rhf_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhf_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rhf_vs_eta_"+algoType;
  TProfile* hp_Rhf_vs_eta = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhf_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rhf_vs_eta->SetYTitle("entries");


  // ----------------  R_x GEN plots:
  histoName = "RchGen_"+algoType;
  TH1F* h1_RchGen = new TH1F(histoName.c_str(), "Jet Charged p_{T} Fraction", 100, 0., 1000.);
  h1_RchGen->SetXTitle("R_{ch}");
  h1_RchGen->SetYTitle("entries");

  histoName = "RchGen_vs_eta_Rch0_"+algoType;
  TProfile* hp_RchGen_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RchGen_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_RchGen_vs_eta_Rch0->SetYTitle("entries");

  histoName = "RchGen_vs_eta_"+algoType;
  TProfile* hp_RchGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RchGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RchGen_vs_eta->SetYTitle("entries");

  histoName = "RchGen_vs_pt_"+algoType;
  TProfile* hp_RchGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RchGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt->SetYTitle("entries");

  histoName = "RchGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RchGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RchGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RchGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RchGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RchGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RgammaGen_"+algoType;
  TH1F* h1_RgammaGen = new TH1F(histoName.c_str(), "Jet Photon Energy Fraction", 100, 0., 1000.);
  h1_RgammaGen->SetXTitle("R_{#gamma}");
  h1_RgammaGen->SetYTitle("entries");

  histoName = "RgammaGen_vs_eta_Rch0_"+algoType;
  TProfile* hp_RgammaGen_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RgammaGen_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_RgammaGen_vs_eta_Rch0->SetYTitle("entries");

  histoName = "RgammaGen_vs_eta_"+algoType;
  TProfile* hp_RgammaGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RgammaGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RgammaGen_vs_eta->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_"+algoType;
  TProfile* hp_RgammaGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RgammaGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RgammaGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RgammaGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RgammaGen_vs_RchGen_"+algoType;
  TProfile* hp_RgammaGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{#gamma} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_RgammaGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RgammaGen_vs_RchGen->SetYTitle("R_{#gamma}");

  histoName = "RnhGen_"+algoType;
  TH1F* h1_RnhGen = new TH1F(histoName.c_str(), "Jet Neutral Hadron Energy Fraction", 100, 0., 1000.);
  h1_RnhGen->SetXTitle("R_{#gamma}");
  h1_RnhGen->SetYTitle("entries");

  histoName = "RnhGen_vs_eta_Rch0_"+algoType;
  TProfile* hp_RnhGen_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RnhGen_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_RnhGen_vs_eta_Rch0->SetYTitle("entries");

  histoName = "RnhGen_vs_eta_"+algoType;
  TProfile* hp_RnhGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RnhGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RnhGen_vs_eta->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_"+algoType;
  TProfile* hp_RnhGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RnhGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RnhGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RnhGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RnhGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RnhGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RnhGen_vs_RchGen_"+algoType;
  TProfile* hp_RnhGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_RnhGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RnhGen_vs_RchGen->SetYTitle("R_{nh}");

  histoName = "RmuGen_"+algoType;
  TH1F* h1_RmuGen = new TH1F(histoName.c_str(), "Jet Neutral Hadron Energy Fraction", 100, 0., 1000.);
  h1_RmuGen->SetXTitle("R_{#gamma}");
  h1_RmuGen->SetYTitle("entries");

  histoName = "RmuGen_vs_eta_"+algoType;
  TProfile* hp_RmuGen_vs_eta = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RmuGen_vs_eta->SetXTitle("#eta^{GEN}");
  hp_RmuGen_vs_eta->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_"+algoType;
  TProfile* hp_RmuGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RmuGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RmuGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RmuGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RmuGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-1, Lower, 0., 1000.);
  hp_RmuGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RmuGen_vs_RchGen_"+algoType;
  TProfile* hp_RmuGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_RmuGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RmuGen_vs_RchGen->SetYTitle("R_{nh}");


  // E_x resolutions:
  histoName = "Ech_resolution_"+algoType;
  TH1F* h1_Ech_resolution = new TH1F(histoName.c_str(), "R_{ch} Resolution", 100, -1.1, 1.1);
  h1_Ech_resolution->SetXTitle("(R_{ch}^{RECO}-R_{ch}^{GEN})/R_{ch}^{GEN} [%]");
  h1_Ech_resolution->SetYTitle("entries");

  histoName = "Egamma_resolution_"+algoType;
  TH1F* h1_Egamma_resolution = new TH1F(histoName.c_str(), "R_{#gamma} Resolution", 100, -1.1, 1.1);
  h1_Egamma_resolution->SetXTitle("(R_{#gamma}^{RECO}-R_{#gamma}^{GEN})/R_{#gamma}^{GEN} [%]");
  h1_Egamma_resolution->SetYTitle("entries");

  histoName = "Enh_resolution_"+algoType;
  TH1F* h1_Enh_resolution = new TH1F(histoName.c_str(), "R_{nh} Resolution", 100, -1.1, 1.1);
  h1_Enh_resolution->SetXTitle("(R_{nh}^{RECO}-R_{nh}^{GEN})/R_{nh}^{GEN} [%]");
  h1_Enh_resolution->SetYTitle("entries");

  histoName = "Egammanh_resolution_"+algoType;
  TH1F* h1_Egammanh_resolution = new TH1F(histoName.c_str(), "R_{nh} Resolution", 100, -1.1, 1.1);
  h1_Egammanh_resolution->SetXTitle("(R_{nh}^{RECO}-R_{nh}^{GEN})/R_{nh}^{GEN} [%]");
  h1_Egammanh_resolution->SetYTitle("entries");

  histoName = "Ee_resolution_"+algoType;
  TH1F* h1_Ee_resolution = new TH1F(histoName.c_str(), "E_{e} Resolution", 100, -1.1, 1.1);
  h1_Ee_resolution->SetXTitle("(E_{e}^{RECO}-E_{e}^{GEN})/E_{e}^{GEN} [%]");
  h1_Ee_resolution->SetYTitle("entries");


  histoName = "Ech_resolution_50_80_"+algoType;
  TH1F* h1_Ech_resolution_50_80 = new TH1F(histoName.c_str(), "R_{ch} Resolution", 100, -1.1, 1.1);
  h1_Ech_resolution_50_80->SetXTitle("(R_{ch}^{RECO}-R_{ch}^{GEN})/R_{ch}^{GEN} [%]");
  h1_Ech_resolution_50_80->SetYTitle("entries");

  histoName = "Egamma_resolution_50_80_"+algoType;
  TH1F* h1_Egamma_resolution_50_80 = new TH1F(histoName.c_str(), "R_{#gamma} Resolution", 100, -1.1, 1.1);
  h1_Egamma_resolution_50_80->SetXTitle("(R_{#gamma}^{RECO}-R_{#gamma}^{GEN})/R_{#gamma}^{GEN} [%]");
  h1_Egamma_resolution_50_80->SetYTitle("entries");

  histoName = "Enh_resolution_50_80_"+algoType;
  TH1F* h1_Enh_resolution_50_80 = new TH1F(histoName.c_str(), "R_{nh} Resolution", 100, -1.1, 1.1);
  h1_Enh_resolution_50_80->SetXTitle("(R_{nh}^{RECO}-R_{nh}^{GEN})/R_{nh}^{GEN} [%]");
  h1_Enh_resolution_50_80->SetYTitle("entries");

  histoName = "Egammanh_resolution_50_80_"+algoType;
  TH1F* h1_Egammanh_resolution_50_80 = new TH1F(histoName.c_str(), "R_{nh} Resolution", 100, -1.1, 1.1);
  h1_Egammanh_resolution_50_80->SetXTitle("(R_{nh}^{RECO}-R_{nh}^{GEN})/R_{nh}^{GEN} [%]");
  h1_Egammanh_resolution_50_80->SetYTitle("entries");

  histoName = "Ee_resolution_50_80_"+algoType;
  TH1F* h1_Ee_resolution_50_80 = new TH1F(histoName.c_str(), "E_{e} Resolution", 100, -1.1, 1.1);
  h1_Ee_resolution_50_80->SetXTitle("(E_{e}^{RECO}-E_{e}^{GEN})/E_{e}^{GEN} [%]");
  h1_Ee_resolution_50_80->SetYTitle("entries");





  int nEntries = tree->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    Jet firstJet, secondJet;

    firstJet.ptReco = 0.;
    secondJet.ptReco = 0.;

    for(int iRecoJet = 0; iRecoJet<nReco; ++iRecoJet) {

      if( iRecoJet >= 40 ) continue;

      Jet iJet;

      iJet.etaReco = etaReco[iRecoJet];
      iJet.phiReco = phiReco[iRecoJet];
      iJet.ptReco = ptReco[iRecoJet];
      iJet.eReco = eReco[iRecoJet];

      iJet.etaGen = etaGen[iRecoJet];
      iJet.phiGen = phiGen[iRecoJet];
      iJet.ptGen = ptGen[iRecoJet];
      iJet.eGen = eGen[iRecoJet];


      //response_onlyCh = 0.;

      if( algoType == "PFItCone5") {
        iJet.Rch = eTracksReco[iRecoJet]/iJet.eReco;
        iJet.Rgamma = ePhotonsReco[iRecoJet]/iJet.eReco;
        iJet.Rnh = eNeutralHadronsReco[iRecoJet]/iJet.eReco;
        iJet.Rmu = eMuonsReco[iRecoJet]/iJet.eReco;
        iJet.Re = eElectronsReco[iRecoJet]/iJet.eReco;
        iJet.Rhf = eHFhadronsReco[iRecoJet]/iJet.eReco;

        //response_onlyCh =ptTracksReco[iRecoJet]/iPtGen;

        iJet.RchGen = eTracksGen[iRecoJet]/iJet.eGen;
        iJet.RgammaGen = ePhotonsGen[iRecoJet]/iJet.eGen;
        iJet.RnhGen = eNeutralHadronsGen[iRecoJet]/iJet.eGen;
        iJet.RmuGen = eMuonsGen[iRecoJet]/iJet.eGen;

      }

      Float_t deltaEta = iJet.etaReco-iJet.etaGen;
      Float_t deltaPhi = iJet.phiReco-iJet.phiGen;
      Float_t pi = 3.14159;
      if(deltaPhi < -pi) deltaPhi += 2.*pi;
      if(deltaPhi > pi) deltaPhi -= 2.*pi;

      Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

      bool matched = false;

      if( algoType=="PFItCone5" ) {
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

      Float_t iEtaReco = iJet->etaReco;
      Float_t iPhiReco = iJet->phiReco;
      Float_t iEReco = iJet->eReco;
      Float_t iPtReco = iJet->ptReco;

      Float_t iEtaGen = iJet->etaGen;
      Float_t iPhiGen = iJet->phiGen;
      Float_t iEGen = iJet->eGen;
      Float_t iPtGen = iJet->ptGen;

      Float_t response = iPtReco/iPtGen;

      Float_t Rch = iJet->Rch;
      Float_t Rgamma = iJet->Rgamma;
      Float_t Rnh = iJet->Rnh;
      Float_t Rmu = iJet->Rmu;
      Float_t Rhf = iJet->Rhf;
      Float_t Re = iJet->Re;

      Float_t RchGen = iJet->RchGen;
      Float_t RgammaGen = iJet->RgammaGen;
      Float_t RnhGen = iJet->RnhGen;
      Float_t RmuGen = iJet->RmuGen;
      Float_t ReGen = iJet->ReGen;

      Float_t Ech = Rch*iPtReco;
      Float_t Egamma = Rgamma*iPtReco;
      Float_t Enh = Rnh*iPtReco;
      Float_t Egammanh = Egamma+Enh;

      Float_t EchGen = RchGen*iPtGen;
      Float_t EgammaGen = RgammaGen*iPtGen;
      Float_t EnhGen = RnhGen*iPtGen;
      Float_t EgammanhGen = EgammaGen+EnhGen;

      if( (algoType=="PFItCone5")&&(fabs(iEtaGen)<2.5)&&(iPtGen>15.) ) {
      
        h1_Ech_resolution->Fill( (Ech-EchGen)/EchGen );
        h1_Egamma_resolution->Fill( (Egamma-EgammaGen)/EgammaGen );
        h1_Enh_resolution->Fill( (Enh-EnhGen)/EnhGen );
        h1_Egammanh_resolution->Fill( (Egammanh-EgammanhGen)/EgammanhGen );

        if( (iPtGen>50.)&&(iPtGen<80.) ) {
          h1_Ech_resolution_50_80->Fill( (Ech-EchGen)/EchGen );
          h1_Egamma_resolution_50_80->Fill( (Egamma-EgammaGen)/EgammaGen );
          h1_Enh_resolution_50_80->Fill( (Enh-EnhGen)/EnhGen );
          h1_Egammanh_resolution_50_80->Fill( (Egammanh-EgammanhGen)/EgammanhGen );
        }

      }

      h2_response_vs_eta->Fill(iEtaGen, response);

      if( Rch==0. ) {
        hp_RchGen_vs_eta_Rch0->Fill(iEtaGen, RchGen);
        hp_RgammaGen_vs_eta_Rch0->Fill(iEtaGen, RgammaGen);
        hp_RnhGen_vs_eta_Rch0->Fill(iEtaGen, RnhGen);
      }

      Float_t deltaEta = fabs(iEtaReco) - fabs(iEtaGen);
      Float_t deltaPhi = iPhiReco-iPhiGen;
      Float_t pi = 3.14159;
      if(deltaPhi < -pi) deltaPhi += 2.*pi;
      if(deltaPhi > pi) deltaPhi -= 2.*pi;

   
      if( fabs(iEtaGen)<1.4 ) {
        h2_deltaEta_vs_pt_barrel->Fill(iPtGen, deltaEta);
        h2_deltaPhi_vs_pt_barrel->Fill(iPtGen, deltaPhi);
        h2_response_vs_pt_barrel->Fill(iPtGen, response);

        if( algoType=="PFItCone5" ) {

          //h2_response_vs_pt_barrel_onlyCh->Fill(iPtGen, response_onlyCh);

          hp_Rgamma_vs_Rch_barrel->Fill(Rch, Rgamma);
          hp_Rnh_vs_Rch_barrel->Fill(Rch, Rnh);

          hp_Rgamma_vs_Rch->Fill(Rch, Rgamma);
          hp_Rnh_vs_Rch->Fill(Rnh, Rgamma);
          hp_Rmu_vs_Rch->Fill(Rmu, Rgamma);

          h1_Rch->Fill(Rch);
          hp_Rch_vs_eta->Fill(iEtaGen, Rch);
          hp_Rch_vs_pt->Fill(iPtGen, Rch);

          h1_Rgamma->Fill(Rgamma);
          hp_Rgamma_vs_eta->Fill(iEtaGen, Rgamma);
          hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma);

          h1_Rnh->Fill(Rnh);
          hp_Rnh_vs_eta->Fill(iEtaGen, Rnh);
          hp_Rnh_vs_pt->Fill(iPtGen, Rnh);

          h1_Rmu->Fill(Rmu);
          hp_Rmu_vs_eta->Fill(iEtaGen, Rmu);
          hp_Rmu_vs_pt->Fill(iPtGen, Rmu);

          h1_Re->Fill(Re);
          hp_Re_vs_eta->Fill(iEtaGen, Re);
          hp_Re_vs_pt->Fill(iPtGen, Re);

          h1_Rhf->Fill(Rhf);
          hp_Rhf_vs_eta->Fill(iEtaGen, Rhf);
          hp_Rhf_vs_pt->Fill(iPtGen, Rhf);

          hp_Rch_vs_pt_barrel->Fill(iPtGen, Rch);
          hp_Rgamma_vs_pt_barrel->Fill(iPtGen, Rgamma);
          hp_Rnh_vs_pt_barrel->Fill(iPtGen, Rnh);
          hp_Rmu_vs_pt_barrel->Fill(iPtGen, Rmu);

          hp_RgammaGen_vs_RchGen->Fill(RchGen, RgammaGen);
          hp_RnhGen_vs_RchGen->Fill(RnhGen, RgammaGen);
          hp_RmuGen_vs_RchGen->Fill(RmuGen, RgammaGen);

          h1_RchGen->Fill(RchGen);
          hp_RchGen_vs_eta->Fill(iEtaGen, RchGen);
          hp_RchGen_vs_pt->Fill(iPtGen, RchGen);

          h1_RgammaGen->Fill(RgammaGen);
          hp_RgammaGen_vs_eta->Fill(iEtaGen, RgammaGen);
          hp_RgammaGen_vs_pt->Fill(iPtGen, RgammaGen);

          h1_RnhGen->Fill(RnhGen);
          hp_RnhGen_vs_eta->Fill(iEtaGen, RnhGen);
          hp_RnhGen_vs_pt->Fill(iPtGen, RnhGen);

          h1_RmuGen->Fill(RmuGen);
          hp_RmuGen_vs_eta->Fill(iEtaGen, RmuGen);
          hp_RmuGen_vs_pt->Fill(iPtGen, RmuGen);

          hp_RchGen_vs_pt_barrel->Fill(iPtGen, RchGen);
          hp_RgammaGen_vs_pt_barrel->Fill(iPtGen, RgammaGen);
          hp_RnhGen_vs_pt_barrel->Fill(iPtGen, RnhGen);
          hp_RmuGen_vs_pt_barrel->Fill(iPtGen, RmuGen);

          if( Rch==0. ) { 
            h1_eta_Rch0->Fill( iEtaGen );
            h2_deltaEta_vs_pt_barrel_Rch0->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rch0->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rch0->Fill(iPtGen, response);
          } else if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
            h2_deltaEta_vs_pt_barrel_Rch10_30->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rch10_30->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rch10_30->Fill(iPtGen, response);
          } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
            h2_deltaEta_vs_pt_barrel_Rch50_70->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rch50_70->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rch50_70->Fill(iPtGen, response);
          } else if( (Rch >= 0.9) && (Rch<=1.) )  {
            h2_deltaEta_vs_pt_barrel_Rch90_100->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rch90_100->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rch90_100->Fill(iPtGen, response);
          }

/*          if( Rgamma==0. ) { 
            h1_eta_Rgamma0->Fill( iEtaGen );
            h2_deltaEta_vs_pt_barrel_Rgamma0->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rgamma0->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rgamma0->Fill(iPtGen, response);
          } else if( (Rgamma >= 0.05) &&( Rgamma <= 0.2 ) ) {
            h2_deltaEta_vs_pt_barrel_Rgamma5_20->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rgamma5_20->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rgamma5_20->Fill(iPtGen, response);
          } else if( (Rgamma >= 0.3) &&( Rgamma < 0.6 ) ) {
            h2_deltaEta_vs_pt_barrel_Rgamma30_60->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rgamma30_60->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rgamma30_60->Fill(iPtGen, response);
          } else if( (Rgamma >= 0.6) && (Rgamma<=1.))  {
            h2_deltaEta_vs_pt_barrel_Rgamma60_100->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_barrel_Rgamma60_100->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_barrel_Rgamma60_100->Fill(iPtGen, response);
          }
*/

          if( iPtGen>70. ) {
            h2_response_vs_Rch->Fill(Rch, response);
            h2_response_vs_Rgamma->Fill(Rgamma, response);
          }

        } //if PF 

      } //if barrel

      if( (fabs(iEtaGen)>1.4) && (fabs(iEtaGen)<3.) ) {
        h2_deltaEta_vs_pt_endcap->Fill(iPtGen, deltaEta);
        h2_deltaPhi_vs_pt_endcap->Fill(iPtGen, deltaPhi);
        h2_response_vs_pt_endcap->Fill(iPtGen, response);

        if( algoType=="PFItCone5" ) {

          if(fabs(iEtaGen)<2.5) {//tracker region
            hp_Rgamma_vs_Rch->Fill(Rch, Rgamma);
            hp_Rnh_vs_Rch->Fill(Rnh, Rgamma);
            hp_Rmu_vs_Rch->Fill(Rmu, Rgamma);
          }

          h1_Rch->Fill(Rch);
          hp_Rch_vs_eta->Fill(iEtaGen, Rch);
          hp_Rch_vs_pt->Fill(iPtGen, Rch);

          h1_Rgamma->Fill(Rgamma);
          hp_Rgamma_vs_eta->Fill(iEtaGen, Rgamma);
          hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma);

          h1_Rnh->Fill(Rnh);
          hp_Rnh_vs_eta->Fill(iEtaGen, Rnh);
          hp_Rnh_vs_pt->Fill(iPtGen, Rnh);

          h1_Re->Fill(Re);
          hp_Re_vs_eta->Fill(iEtaGen, Re);
          hp_Re_vs_pt->Fill(iPtGen, Re);

          h1_Rmu->Fill(Rmu);
          hp_Rmu_vs_eta->Fill(iEtaGen, Rmu);
          hp_Rmu_vs_pt->Fill(iPtGen, Rmu);

          h1_Rhf->Fill(Rhf);
          hp_Rhf_vs_eta->Fill(iEtaGen, Rhf);
          hp_Rhf_vs_pt->Fill(iPtGen, Rhf);

          hp_Rch_vs_pt_endcap->Fill(iPtGen, Rch);
          hp_Rgamma_vs_pt_endcap->Fill(iPtGen, Rgamma);
          hp_Rnh_vs_pt_endcap->Fill(iPtGen, Rnh);
          hp_Rhf_vs_pt_endcap->Fill(iPtGen, Rhf);
          hp_Rmu_vs_pt_endcap->Fill(iPtGen, Rmu);

          hp_RgammaGen_vs_RchGen->Fill(RchGen, RgammaGen);
          hp_RnhGen_vs_RchGen->Fill(RnhGen, RgammaGen);
          hp_RmuGen_vs_RchGen->Fill(RmuGen, RgammaGen);

          h1_RchGen->Fill(RchGen);
          hp_RchGen_vs_eta->Fill(iEtaGen, RchGen);
          hp_RchGen_vs_pt->Fill(iPtGen, RchGen);

          h1_RgammaGen->Fill(RgammaGen);
          hp_RgammaGen_vs_eta->Fill(iEtaGen, RgammaGen);
          hp_RgammaGen_vs_pt->Fill(iPtGen, RgammaGen);

          h1_RnhGen->Fill(RnhGen);
          hp_RnhGen_vs_eta->Fill(iEtaGen, RnhGen);
          hp_RnhGen_vs_pt->Fill(iPtGen, RnhGen);

          h1_RmuGen->Fill(RmuGen);
          hp_RmuGen_vs_eta->Fill(iEtaGen, RmuGen);
          hp_RmuGen_vs_pt->Fill(iPtGen, RmuGen);

          hp_RchGen_vs_pt_endcap->Fill(iPtGen, RchGen);
          hp_RgammaGen_vs_pt_endcap->Fill(iPtGen, RgammaGen);
          hp_RnhGen_vs_pt_endcap->Fill(iPtGen, RnhGen);
          hp_RmuGen_vs_pt_endcap->Fill(iPtGen, RmuGen);



          if( Rch==0. ) { 
            h1_eta_Rch0->Fill( iEtaGen );
            h2_deltaEta_vs_pt_endcap_Rch0->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rch0->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rch0->Fill(iPtGen, response);
          } else if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
            h2_deltaEta_vs_pt_endcap_Rch10_30->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rch10_30->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rch10_30->Fill(iPtGen, response);
          } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
            h2_deltaEta_vs_pt_endcap_Rch50_70->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rch50_70->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rch50_70->Fill(iPtGen, response);
          } else if( (Rch >= 0.9) && (Rch <=1.) )  {
            h2_deltaEta_vs_pt_endcap_Rch90_100->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rch90_100->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rch90_100->Fill(iPtGen, response);
          }

/*          if( Rgamma==0. ) { 
            h1_eta_Rgamma0->Fill( iEtaGen );
            h2_deltaEta_vs_pt_endcap_Rgamma0->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rgamma0->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rgamma0->Fill(iPtGen, response);
          } else if( (Rgamma >= 0.05) &&( Rgamma <= 0.2 ) ) {
            h2_deltaEta_vs_pt_endcap_Rgamma5_20->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rgamma5_20->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rgamma5_20->Fill(iPtGen, response);
          } else if( (Rgamma >= 0.3) &&( Rgamma < 0.6 ) ) {
            h2_deltaEta_vs_pt_endcap_Rgamma30_60->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rgamma30_60->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rgamma30_60->Fill(iPtGen, response);
          } else if( (Rgamma >= 0.6) && (Rgamma<=1.) )  {
            h2_deltaEta_vs_pt_endcap_Rgamma60_100->Fill(iPtGen, deltaEta);
            h2_deltaPhi_vs_pt_endcap_Rgamma60_100->Fill(iPtGen, deltaPhi);
            h2_response_vs_pt_endcap_Rgamma60_100->Fill(iPtGen, response);
          }
*/
          if( iPtGen>70. ) {
            h2_response_vs_Rch->Fill(Rch, response);
            h2_response_vs_Rgamma->Fill(Rgamma, response);
          }

        } //if PF 

      } //if endcaps

      if( (fabs(iEtaGen)>3.) && (fabs(iEtaGen)<5.) ) {
        h2_deltaEta_vs_pt_forward->Fill(iPtGen, deltaEta);
        h2_deltaPhi_vs_pt_forward->Fill(iPtGen, deltaPhi);
        h2_response_vs_pt_forward->Fill(iPtGen, response);

        if( algoType=="PFItCone5" ) {

          h1_Rch->Fill(Rch);
          hp_Rch_vs_eta->Fill(iEtaGen, Rch);
          hp_Rch_vs_pt->Fill(iPtGen, Rch);

          h1_Rgamma->Fill(Rgamma);
          hp_Rgamma_vs_eta->Fill(iEtaGen, Rgamma);
          hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma);

          h1_Rnh->Fill(Rnh);
          hp_Rnh_vs_eta->Fill(iEtaGen, Rnh);
          hp_Rnh_vs_pt->Fill(iPtGen, Rnh);

          h1_Re->Fill(Re);
          hp_Re_vs_eta->Fill(iEtaGen, Re);
          hp_Re_vs_pt->Fill(iPtGen, Re);

          h1_Rmu->Fill(Rmu);
          hp_Rmu_vs_eta->Fill(iEtaGen, Rmu);
          hp_Rmu_vs_pt->Fill(iPtGen, Rmu);

          h1_Rhf->Fill(Rhf);
          hp_Rhf_vs_eta->Fill(iEtaGen, Rhf);
          hp_Rhf_vs_pt->Fill(iPtGen, Rhf);

          hp_RgammaGen_vs_RchGen->Fill(RchGen, RgammaGen);
          hp_RnhGen_vs_RchGen->Fill(RnhGen, RgammaGen);
          hp_RmuGen_vs_RchGen->Fill(RmuGen, RgammaGen);

          h1_RchGen->Fill(RchGen);
          hp_RchGen_vs_eta->Fill(iEtaGen, RchGen);
          hp_RchGen_vs_pt->Fill(iPtGen, RchGen);

          h1_RgammaGen->Fill(RgammaGen);
          hp_RgammaGen_vs_eta->Fill(iEtaGen, RgammaGen);
          hp_RgammaGen_vs_pt->Fill(iPtGen, RgammaGen);

          h1_RnhGen->Fill(RnhGen);
          hp_RnhGen_vs_eta->Fill(iEtaGen, RnhGen);
          hp_RnhGen_vs_pt->Fill(iPtGen, RnhGen);

          h1_RmuGen->Fill(RmuGen);
          hp_RmuGen_vs_eta->Fill(iEtaGen, RmuGen);
          hp_RmuGen_vs_pt->Fill(iPtGen, RmuGen);


        } //if PF


      } //if forward


    } //for recoJets


  } //for entries

 
  histoName = "response_vs_pt_FIT_barrel_"+algoType;
  TH1F* h1_response_vs_pt_FIT_barrel = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_endcap_"+algoType;
  TH1F* h1_response_vs_pt_FIT_endcap = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_forward_"+algoType;
  TH1F* h1_response_vs_pt_FIT_forward = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-1, Lower);
  h1_response_vs_pt_FIT_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_forward->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_pt_FIT_barrel_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_endcap_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_endcap = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_endcap->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_forward_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_forward = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_forward->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_forward->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_barrel = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_endcap_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_endcap = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_forward_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_forward = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_forward->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_pt_RMS_barrel_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_endcap_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_endcap = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_forward_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_forward = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_forward->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_forward->GetYaxis()->SetRangeUser(0., 0.5);

  if(draw) {
    std::string name;
    name = algoType + "_barrel";
    fillProfile(h1_response_vs_pt_FIT_barrel, h1_resolution_vs_pt_FIT_barrel, h1_response_vs_pt_MEAN_barrel, h1_resolution_vs_pt_RMS_barrel, h2_response_vs_pt_barrel, name);
    name = algoType + "_endcap";
    fillProfile(h1_response_vs_pt_FIT_endcap, h1_resolution_vs_pt_FIT_endcap, h1_response_vs_pt_MEAN_endcap, h1_resolution_vs_pt_RMS_endcap, h2_response_vs_pt_endcap, name);
    name = algoType + "_forward";
    fillProfile(h1_response_vs_pt_FIT_forward, h1_resolution_vs_pt_FIT_forward, h1_response_vs_pt_MEAN_forward, h1_resolution_vs_pt_RMS_forward, h2_response_vs_pt_forward, name);
  } else { 
    fillProfile(h1_response_vs_pt_FIT_barrel, h1_resolution_vs_pt_FIT_barrel, h1_response_vs_pt_MEAN_barrel, h1_resolution_vs_pt_RMS_barrel, h2_response_vs_pt_barrel);
    fillProfile(h1_response_vs_pt_FIT_endcap, h1_resolution_vs_pt_FIT_endcap, h1_response_vs_pt_MEAN_endcap, h1_resolution_vs_pt_RMS_endcap, h2_response_vs_pt_endcap);
    fillProfile(h1_response_vs_pt_FIT_forward, h1_resolution_vs_pt_FIT_forward, h1_response_vs_pt_MEAN_forward, h1_resolution_vs_pt_RMS_forward, h2_response_vs_pt_forward);
  }

  histoName = "response_vs_eta_FIT_"+algoType;
  TH1F* h1_response_vs_eta_FIT = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5., 5.);
  h1_response_vs_eta_FIT->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_FIT->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_eta_FIT->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_eta_FIT_"+algoType;
  TH1F* h1_resolution_vs_eta_FIT = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet #eta", 100, -5., 5.);
  h1_resolution_vs_eta_FIT->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_FIT->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_FIT->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_eta_MEAN_"+algoType;
  TH1F* h1_response_vs_eta_MEAN = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5., 5.);
  h1_response_vs_eta_MEAN->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_MEAN->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_response_vs_eta_MEAN->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_eta_RMS_"+algoType;
  TH1F* h1_resolution_vs_eta_RMS = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet #eta", 100, -5., 5.);
  h1_resolution_vs_eta_RMS->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_RMS->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_RMS->GetYaxis()->SetRangeUser(0., 0.5);


  fillProfile(h1_response_vs_eta_FIT, h1_resolution_vs_eta_FIT, h1_response_vs_eta_MEAN, h1_resolution_vs_eta_RMS, h2_response_vs_eta);
  

  //PF only histograms: 
/*
  histoName = "response_vs_pt_FIT_barrel_onlyCh_"+algoType;
  TH1F* h1_response_vs_pt_FIT_barrel_onlyCh = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_onlyCh->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_onlyCh->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_onlyCh->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_barrel_onlyCh_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_barrel_onlyCh = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_onlyCh->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_onlyCh->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_onlyCh->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "resolution_vs_pt_FIT_barrel_onlyCh_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_barrel_onlyCh = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_onlyCh->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_onlyCh->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel_onlyCh->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_barrel_onlyCh_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_barrel_onlyCh = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_onlyCh->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_onlyCh->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_onlyCh->GetYaxis()->SetRangeUser(0., 0.5);
*/ 
  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch10_30->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch50_70->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch90_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch10_30->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch50_70->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch90_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma5_20->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma30_60->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma60_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma5_20->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma30_60->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-1, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma60_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);




  histoName = "response_vs_Rch_FIT_"+algoType;
  TH1F* h1_response_vs_Rch_FIT = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{ch}", 50, 0., 1.);
  h1_response_vs_Rch_FIT->SetXTitle("R_{ch}");
  h1_response_vs_Rch_FIT->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rch_FIT->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_Rgamma_FIT_"+algoType;
  TH1F* h1_response_vs_Rgamma_FIT = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{gamma}", 50, 0., 1.);
  h1_response_vs_Rgamma_FIT->SetXTitle("R_{#gamma}");
  h1_response_vs_Rgamma_FIT->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rgamma_FIT->GetYaxis()->SetRangeUser(0.2, 1.);


  histoName = "response_vs_Rch_MEAN_"+algoType;
  TH1F* h1_response_vs_Rch_MEAN = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{ch}", 50, 0., 1.);
  h1_response_vs_Rch_MEAN->SetXTitle("R_{ch}");
  h1_response_vs_Rch_MEAN->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rch_MEAN->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_Rgamma_MEAN_"+algoType;
  TH1F* h1_response_vs_Rgamma_MEAN = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{gamma}", 50, 0., 1.);
  h1_response_vs_Rgamma_MEAN->SetXTitle("R_{#gamma}");
  h1_response_vs_Rgamma_MEAN->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rgamma_MEAN->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_Rch_FIT_"+algoType;
  TH1F* h1_resolution_vs_Rch_FIT = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{ch}", 50, 0., 1.);
  h1_resolution_vs_Rch_FIT->SetXTitle("R_{ch}");
  h1_resolution_vs_Rch_FIT->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rch_FIT->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_Rgamma_FIT_"+algoType;
  TH1F* h1_resolution_vs_Rgamma_FIT = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{#gamma}", 50, 0., 1.);
  h1_resolution_vs_Rgamma_FIT->SetXTitle("R_{#gamma}");
  h1_resolution_vs_Rgamma_FIT->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rgamma_FIT->GetYaxis()->SetRangeUser(0., 0.5);


  histoName = "resolution_vs_Rch_RMS_"+algoType;
  TH1F* h1_resolution_vs_Rch_RMS = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{ch}", 50, 0., 1.);
  h1_resolution_vs_Rch_RMS->SetXTitle("R_{ch}");
  h1_resolution_vs_Rch_RMS->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rch_RMS->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_Rgamma_RMS_"+algoType;
  TH1F* h1_resolution_vs_Rgamma_RMS = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{#gamma}", 50, 0., 1.);
  h1_resolution_vs_Rgamma_RMS->SetXTitle("R_{#gamma}");
  h1_resolution_vs_Rgamma_RMS->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rgamma_RMS->GetYaxis()->SetRangeUser(0., 0.5);


  if( algoType == "PFItCone5" ) { 
std::cerr << "sto per fare Rch etc" << std::endl;
    if(draw) {

      fillProfile(h1_response_vs_pt_FIT_barrel_Rch0, h1_resolution_vs_pt_FIT_barrel_Rch0, h1_response_vs_pt_MEAN_barrel_Rch0, h1_resolution_vs_pt_RMS_barrel_Rch0, h2_response_vs_pt_barrel_Rch0, "Rch0_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch0, h1_resolution_vs_pt_FIT_endcap_Rch0, h1_response_vs_pt_MEAN_endcap_Rch0, h1_resolution_vs_pt_RMS_endcap_Rch0, h2_response_vs_pt_endcap_Rch0, "Rch0_endcap");
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch10_30, h1_resolution_vs_pt_FIT_barrel_Rch10_30, h1_response_vs_pt_MEAN_barrel_Rch10_30, h1_resolution_vs_pt_RMS_barrel_Rch10_30, h2_response_vs_pt_barrel_Rch10_30, "Rch10_30_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch10_30, h1_resolution_vs_pt_FIT_endcap_Rch10_30, h1_response_vs_pt_MEAN_endcap_Rch10_30, h1_resolution_vs_pt_RMS_endcap_Rch10_30, h2_response_vs_pt_endcap_Rch10_30, "Rch10_30_endcap");
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch50_70, h1_resolution_vs_pt_FIT_barrel_Rch50_70, h1_response_vs_pt_MEAN_barrel_Rch50_70, h1_resolution_vs_pt_RMS_barrel_Rch50_70, h2_response_vs_pt_barrel_Rch50_70, "Rch50_70_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch50_70, h1_resolution_vs_pt_FIT_endcap_Rch50_70, h1_response_vs_pt_MEAN_endcap_Rch50_70, h1_resolution_vs_pt_RMS_endcap_Rch50_70, h2_response_vs_pt_endcap_Rch50_70, "Rch50_70_endcap");
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch90_100, h1_resolution_vs_pt_FIT_barrel_Rch90_100, h1_response_vs_pt_MEAN_barrel_Rch90_100, h1_resolution_vs_pt_RMS_barrel_Rch90_100, h2_response_vs_pt_barrel_Rch90_100, "Rch90_100_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch90_100, h1_resolution_vs_pt_FIT_endcap_Rch90_100, h1_response_vs_pt_MEAN_endcap_Rch90_100, h1_resolution_vs_pt_RMS_endcap_Rch90_100, h2_response_vs_pt_endcap_Rch90_100, "Rch90_100_endcap");
    } else {

//      fillProfile(h1_response_vs_pt_FIT_barrel_onlyCh, h1_resolution_vs_pt_FIT_barrel_onlyCh, h1_response_vs_pt_MEAN_barrel_onlyCh, h1_resolution_vs_pt_RMS_barrel_onlyCh, h2_response_vs_pt_barrel_onlyCh);

      fillProfile(h1_response_vs_pt_FIT_barrel_Rch0, h1_resolution_vs_pt_FIT_barrel_Rch0, h1_response_vs_pt_MEAN_barrel_Rch0, h1_resolution_vs_pt_RMS_barrel_Rch0, h2_response_vs_pt_barrel_Rch0);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch0, h1_resolution_vs_pt_FIT_endcap_Rch0, h1_response_vs_pt_MEAN_endcap_Rch0, h1_resolution_vs_pt_RMS_endcap_Rch0, h2_response_vs_pt_endcap_Rch0);
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch10_30, h1_resolution_vs_pt_FIT_barrel_Rch10_30, h1_response_vs_pt_MEAN_barrel_Rch10_30, h1_resolution_vs_pt_RMS_barrel_Rch10_30, h2_response_vs_pt_barrel_Rch10_30);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch10_30, h1_resolution_vs_pt_FIT_endcap_Rch10_30, h1_response_vs_pt_MEAN_endcap_Rch10_30, h1_resolution_vs_pt_RMS_endcap_Rch10_30, h2_response_vs_pt_endcap_Rch10_30);
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch50_70, h1_resolution_vs_pt_FIT_barrel_Rch50_70, h1_response_vs_pt_MEAN_barrel_Rch50_70, h1_resolution_vs_pt_RMS_barrel_Rch50_70, h2_response_vs_pt_barrel_Rch50_70);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch50_70, h1_resolution_vs_pt_FIT_endcap_Rch50_70, h1_response_vs_pt_MEAN_endcap_Rch50_70, h1_resolution_vs_pt_RMS_endcap_Rch50_70, h2_response_vs_pt_endcap_Rch50_70);
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch90_100, h1_resolution_vs_pt_FIT_barrel_Rch90_100, h1_response_vs_pt_MEAN_barrel_Rch90_100, h1_resolution_vs_pt_RMS_barrel_Rch90_100, h2_response_vs_pt_barrel_Rch90_100);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch90_100, h1_resolution_vs_pt_FIT_endcap_Rch90_100, h1_response_vs_pt_MEAN_endcap_Rch90_100, h1_resolution_vs_pt_RMS_endcap_Rch90_100, h2_response_vs_pt_endcap_Rch90_100);
    }

/*    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma0, h1_resolution_vs_pt_FIT_barrel_Rgamma0, h1_response_vs_pt_MEAN_barrel_Rgamma0, h1_resolution_vs_pt_RMS_barrel_Rgamma0, h2_response_vs_pt_barrel_Rgamma0);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma0, h1_resolution_vs_pt_FIT_endcap_Rgamma0, h1_response_vs_pt_MEAN_endcap_Rgamma0, h1_resolution_vs_pt_RMS_endcap_Rgamma0, h2_response_vs_pt_endcap_Rgamma0);
    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma5_20, h1_resolution_vs_pt_FIT_barrel_Rgamma5_20, h1_response_vs_pt_MEAN_barrel_Rgamma5_20, h1_resolution_vs_pt_RMS_barrel_Rgamma5_20, h2_response_vs_pt_barrel_Rgamma5_20);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma5_20, h1_resolution_vs_pt_FIT_endcap_Rgamma5_20, h1_response_vs_pt_MEAN_endcap_Rgamma5_20, h1_resolution_vs_pt_RMS_endcap_Rgamma5_20, h2_response_vs_pt_endcap_Rgamma5_20);
    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma30_60, h1_resolution_vs_pt_FIT_barrel_Rgamma30_60, h1_response_vs_pt_MEAN_barrel_Rgamma30_60, h1_resolution_vs_pt_RMS_barrel_Rgamma30_60, h2_response_vs_pt_barrel_Rgamma30_60);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma30_60, h1_resolution_vs_pt_FIT_endcap_Rgamma30_60, h1_response_vs_pt_MEAN_endcap_Rgamma30_60, h1_resolution_vs_pt_RMS_endcap_Rgamma30_60, h2_response_vs_pt_endcap_Rgamma30_60);
    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma60_100, h1_resolution_vs_pt_FIT_barrel_Rgamma60_100, h1_response_vs_pt_MEAN_barrel_Rgamma60_100, h1_resolution_vs_pt_RMS_barrel_Rgamma60_100, h2_response_vs_pt_barrel_Rgamma60_100);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma60_100, h1_resolution_vs_pt_FIT_endcap_Rgamma60_100, h1_response_vs_pt_MEAN_endcap_Rgamma60_100, h1_resolution_vs_pt_RMS_endcap_Rgamma60_100, h2_response_vs_pt_endcap_Rgamma60_100);
*/
    fillProfile(h1_response_vs_Rch_FIT, h1_resolution_vs_Rch_FIT, h1_response_vs_Rch_MEAN, h1_resolution_vs_Rch_RMS, h2_response_vs_Rch);
    fillProfile(h1_response_vs_Rgamma_FIT, h1_resolution_vs_Rgamma_FIT, h1_response_vs_Rgamma_MEAN, h1_resolution_vs_Rgamma_RMS, h2_response_vs_Rgamma);
  }
  
  histoName = "sigmaEta_vs_pt_barrel_"+algoType;
  TH1F* h1_sigmaEta_vs_pt_barrel = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType;
  TH1F* h1_sigmaEta_vs_pt_endcap = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_forward_"+algoType;
  TH1F* h1_sigmaEta_vs_pt_forward = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_forward->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_forward->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_barrel_"+algoType;
  TH1F* h1_sigmaPhi_vs_pt_barrel = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_barrel->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType;
  TH1F* h1_sigmaPhi_vs_pt_endcap = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_endcap->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_forward_"+algoType;
  TH1F* h1_sigmaPhi_vs_pt_forward = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_forward->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_forward->GetYaxis()->SetRangeUser(0., 0.2);

std::cerr << "sto per fare posres" << std::endl;
  fillPositionResolution(h1_sigmaEta_vs_pt_barrel, h1_sigmaPhi_vs_pt_barrel, h2_deltaEta_vs_pt_barrel, h2_deltaPhi_vs_pt_barrel);
  fillPositionResolution(h1_sigmaEta_vs_pt_endcap, h1_sigmaPhi_vs_pt_endcap, h2_deltaEta_vs_pt_endcap, h2_deltaPhi_vs_pt_endcap);
  fillPositionResolution(h1_sigmaEta_vs_pt_forward, h1_sigmaPhi_vs_pt_forward, h2_deltaEta_vs_pt_forward, h2_deltaPhi_vs_pt_forward);

std::cerr << "fatto posres" << std::endl;




 

  // **************        PF only  


  //--------------------------------- Rch ------------------------------------
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch0";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch0";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
  
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch0";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch0";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
  
  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);


  //-------------------- Rgamma --------------------------------------------
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-1, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);



 
  if( algoType=="PFItCone5" ) {

std::cerr << "sto per fare pos res" << std::endl;
    //Rch
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rch0, h1_sigmaPhi_vs_pt_barrel_Rch0, h2_deltaEta_vs_pt_barrel_Rch0, h2_deltaPhi_vs_pt_barrel_Rch0);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rch0, h1_sigmaPhi_vs_pt_endcap_Rch0, h2_deltaEta_vs_pt_endcap_Rch0, h2_deltaPhi_vs_pt_endcap_Rch0);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rch10_30, h1_sigmaPhi_vs_pt_barrel_Rch10_30, h2_deltaEta_vs_pt_barrel_Rch10_30, h2_deltaPhi_vs_pt_barrel_Rch10_30);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rch10_30, h1_sigmaPhi_vs_pt_endcap_Rch10_30, h2_deltaEta_vs_pt_endcap_Rch10_30, h2_deltaPhi_vs_pt_endcap_Rch10_30);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rch50_70, h1_sigmaPhi_vs_pt_barrel_Rch50_70, h2_deltaEta_vs_pt_barrel_Rch50_70, h2_deltaPhi_vs_pt_barrel_Rch50_70);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rch50_70, h1_sigmaPhi_vs_pt_endcap_Rch50_70, h2_deltaEta_vs_pt_endcap_Rch50_70, h2_deltaPhi_vs_pt_endcap_Rch50_70);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rch90_100, h1_sigmaPhi_vs_pt_barrel_Rch90_100, h2_deltaEta_vs_pt_barrel_Rch90_100, h2_deltaPhi_vs_pt_barrel_Rch90_100);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rch90_100, h1_sigmaPhi_vs_pt_endcap_Rch90_100, h2_deltaEta_vs_pt_endcap_Rch90_100, h2_deltaPhi_vs_pt_endcap_Rch90_100);
   
    //Rgamma
/*    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma0,      h1_sigmaPhi_vs_pt_barrel_Rgamma0,      h2_deltaEta_vs_pt_barrel_Rgamma0,      h2_deltaPhi_vs_pt_barrel_Rgamma0);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma0,      h1_sigmaPhi_vs_pt_endcap_Rgamma0,      h2_deltaEta_vs_pt_endcap_Rgamma0,      h2_deltaPhi_vs_pt_endcap_Rgamma0);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma5_20,  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20,  h2_deltaEta_vs_pt_barrel_Rgamma5_20,  h2_deltaPhi_vs_pt_barrel_Rgamma5_20);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma5_20,  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20,  h2_deltaEta_vs_pt_endcap_Rgamma5_20,  h2_deltaPhi_vs_pt_endcap_Rgamma5_20);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma30_60,  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60,  h2_deltaEta_vs_pt_barrel_Rgamma30_60,  h2_deltaPhi_vs_pt_barrel_Rgamma30_60);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma30_60,  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60,  h2_deltaEta_vs_pt_endcap_Rgamma30_60,  h2_deltaPhi_vs_pt_endcap_Rgamma30_60);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma60_100, h1_sigmaPhi_vs_pt_barrel_Rgamma60_100, h2_deltaEta_vs_pt_barrel_Rgamma60_100, h2_deltaPhi_vs_pt_barrel_Rgamma60_100);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma60_100, h1_sigmaPhi_vs_pt_endcap_Rgamma60_100, h2_deltaEta_vs_pt_endcap_Rgamma60_100, h2_deltaPhi_vs_pt_endcap_Rgamma60_100);
*/  }

 
  delete h2_response_vs_pt_barrel;
  h2_response_vs_pt_barrel = 0;
  delete h2_response_vs_pt_endcap;
  h2_response_vs_pt_endcap = 0;
  delete h2_response_vs_pt_forward;
  h2_response_vs_pt_forward = 0;

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


  file->cd();


  h1_Ech_resolution->Write();
  h1_Egamma_resolution->Write();
  h1_Enh_resolution->Write();
  h1_Egammanh_resolution->Write();

  h1_Ech_resolution_50_80->Write();
  h1_Egamma_resolution_50_80->Write();
  h1_Enh_resolution_50_80->Write();
  h1_Egammanh_resolution_50_80->Write();

  h1_deltaR_barrelEndcap->Write();
  h1_deltaR_forward->Write();
  h1_deltaR_veryForward->Write();

  h1_response_vs_eta_MEAN->Write();
  h1_response_vs_eta_FIT->Write();
  h1_resolution_vs_eta_RMS->Write();
  h1_resolution_vs_eta_FIT->Write();

  h1_response_vs_pt_MEAN_barrel->Write();
  h1_response_vs_pt_MEAN_endcap->Write();
  h1_response_vs_pt_MEAN_forward->Write();

  h1_resolution_vs_pt_RMS_barrel->Write();
  h1_resolution_vs_pt_RMS_endcap->Write();
  h1_resolution_vs_pt_RMS_forward->Write();

  h1_response_vs_pt_FIT_barrel->Write();
  h1_response_vs_pt_FIT_endcap->Write();
  h1_response_vs_pt_FIT_forward->Write();

  h1_resolution_vs_pt_FIT_barrel->Write();
  h1_resolution_vs_pt_FIT_endcap->Write();
  h1_resolution_vs_pt_FIT_forward->Write();

  h1_sigmaEta_vs_pt_barrel->Write();
  h1_sigmaEta_vs_pt_endcap->Write();
  h1_sigmaEta_vs_pt_forward->Write();

  h1_sigmaPhi_vs_pt_barrel->Write();
  h1_sigmaPhi_vs_pt_endcap->Write();
  h1_sigmaPhi_vs_pt_forward->Write();

  if( algoType == "PFItCone5" ) { 

//    h1_response_vs_pt_MEAN_barrel_onlyCh->Write();
//    h1_response_vs_pt_FIT_barrel_onlyCh->Write();
//    h1_resolution_vs_pt_RMS_barrel_onlyCh->Write();
//    h1_resolution_vs_pt_FIT_barrel_onlyCh->Write();

    hp_Rgamma_vs_Rch_barrel->Write();
    hp_Rnh_vs_Rch_barrel->Write();

    hp_Rgamma_vs_Rch->Write();
    hp_Rnh_vs_Rch->Write();
    hp_Rmu_vs_Rch->Write();

    h1_Rch->Write();
    hp_Rch_vs_eta->Write();
    hp_Rch_vs_pt->Write();

    h1_Rgamma->Write();
    hp_Rgamma_vs_eta->Write();
    hp_Rgamma_vs_pt->Write();

    h1_Rnh->Write();
    hp_Rnh_vs_eta->Write();
    hp_Rnh_vs_pt->Write();

    h1_Re->Write();
    hp_Re_vs_eta->Write();
    hp_Re_vs_pt->Write();

    h1_Rmu->Write();
    hp_Rmu_vs_eta->Write();
    hp_Rmu_vs_pt->Write();

    h1_Rhf->Write();
    hp_Rhf_vs_eta->Write();
    hp_Rhf_vs_pt->Write();

    hp_Rch_vs_pt_barrel->Write();
    hp_Rgamma_vs_pt_barrel->Write();
    hp_Rnh_vs_pt_barrel->Write();
    hp_Rmu_vs_pt_barrel->Write();

    hp_Rch_vs_pt_endcap->Write();
    hp_Rgamma_vs_pt_endcap->Write();
    hp_Rnh_vs_pt_endcap->Write();
    hp_Rhf_vs_pt_endcap->Write();
    hp_Rmu_vs_pt_endcap->Write();

    hp_RgammaGen_vs_RchGen->Write();
    hp_RnhGen_vs_RchGen->Write();
    hp_RmuGen_vs_RchGen->Write();

    h1_RchGen->Write();
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

    h1_response_vs_Rch_FIT->Write();
    h1_response_vs_Rch_MEAN->Write();
    h1_response_vs_Rgamma_FIT->Write();
    h1_response_vs_Rgamma_MEAN->Write();

    h1_resolution_vs_Rch_FIT->Write();
    h1_resolution_vs_Rch_RMS->Write();
    h1_resolution_vs_Rgamma_FIT->Write();
    h1_resolution_vs_Rgamma_RMS->Write();

    h1_response_vs_pt_MEAN_barrel_Rch0->Write();
    h1_response_vs_pt_MEAN_endcap_Rch0->Write();
    h1_response_vs_pt_MEAN_barrel_Rch10_30->Write();
    h1_response_vs_pt_MEAN_endcap_Rch10_30->Write();
    h1_response_vs_pt_MEAN_barrel_Rch50_70->Write();
    h1_response_vs_pt_MEAN_endcap_Rch50_70->Write();
    h1_response_vs_pt_MEAN_barrel_Rch90_100->Write();
    h1_response_vs_pt_MEAN_endcap_Rch90_100->Write();

    h1_resolution_vs_pt_RMS_barrel_Rch0->Write();
    h1_resolution_vs_pt_RMS_endcap_Rch0->Write();
    h1_resolution_vs_pt_RMS_barrel_Rch10_30->Write();
    h1_resolution_vs_pt_RMS_endcap_Rch10_30->Write();
    h1_resolution_vs_pt_RMS_barrel_Rch50_70->Write();
    h1_resolution_vs_pt_RMS_endcap_Rch50_70->Write();
    h1_resolution_vs_pt_RMS_barrel_Rch90_100->Write();
    h1_resolution_vs_pt_RMS_endcap_Rch90_100->Write();

    h1_response_vs_pt_FIT_barrel_Rch0->Write();
    h1_response_vs_pt_FIT_endcap_Rch0->Write();
    h1_response_vs_pt_FIT_barrel_Rch10_30->Write();
    h1_response_vs_pt_FIT_endcap_Rch10_30->Write();
    h1_response_vs_pt_FIT_barrel_Rch50_70->Write();
    h1_response_vs_pt_FIT_endcap_Rch50_70->Write();
    h1_response_vs_pt_FIT_barrel_Rch90_100->Write();
    h1_response_vs_pt_FIT_endcap_Rch90_100->Write();

    h1_resolution_vs_pt_FIT_barrel_Rch0->Write();
    h1_resolution_vs_pt_FIT_endcap_Rch0->Write();
    h1_resolution_vs_pt_FIT_barrel_Rch10_30->Write();
    h1_resolution_vs_pt_FIT_endcap_Rch10_30->Write();
    h1_resolution_vs_pt_FIT_barrel_Rch50_70->Write();
    h1_resolution_vs_pt_FIT_endcap_Rch50_70->Write();
    h1_resolution_vs_pt_FIT_barrel_Rch90_100->Write();
    h1_resolution_vs_pt_FIT_endcap_Rch90_100->Write();

    h1_response_vs_pt_FIT_barrel_Rgamma0->Write();
    h1_response_vs_pt_FIT_endcap_Rgamma0->Write();
    h1_response_vs_pt_FIT_barrel_Rgamma5_20->Write();
    h1_response_vs_pt_FIT_endcap_Rgamma5_20->Write();
    h1_response_vs_pt_FIT_barrel_Rgamma30_60->Write();
    h1_response_vs_pt_FIT_endcap_Rgamma30_60->Write();
    h1_response_vs_pt_FIT_barrel_Rgamma60_100->Write();
    h1_response_vs_pt_FIT_endcap_Rgamma60_100->Write();

    h1_resolution_vs_pt_FIT_barrel_Rgamma0->Write();
    h1_resolution_vs_pt_FIT_endcap_Rgamma0->Write();
    h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->Write();
    h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->Write();
    h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->Write();
    h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->Write();
    h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->Write();
    h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->Write();


    h1_sigmaEta_vs_pt_barrel_Rch0->Write();
    h1_sigmaEta_vs_pt_endcap_Rch0->Write();
    h1_sigmaEta_vs_pt_barrel_Rch10_30->Write();
    h1_sigmaEta_vs_pt_endcap_Rch10_30->Write();
    h1_sigmaEta_vs_pt_barrel_Rch50_70->Write();
    h1_sigmaEta_vs_pt_endcap_Rch50_70->Write();
    h1_sigmaEta_vs_pt_barrel_Rch90_100->Write();
    h1_sigmaEta_vs_pt_endcap_Rch90_100->Write();

    h1_sigmaPhi_vs_pt_barrel_Rch0->Write();
    h1_sigmaPhi_vs_pt_endcap_Rch0->Write();
    h1_sigmaPhi_vs_pt_barrel_Rch10_30->Write();
    h1_sigmaPhi_vs_pt_endcap_Rch10_30->Write();
    h1_sigmaPhi_vs_pt_barrel_Rch50_70->Write();
    h1_sigmaPhi_vs_pt_endcap_Rch50_70->Write();
    h1_sigmaPhi_vs_pt_barrel_Rch90_100->Write();
    h1_sigmaPhi_vs_pt_endcap_Rch90_100->Write();

    h1_sigmaEta_vs_pt_barrel_Rgamma0->Write();
    h1_sigmaEta_vs_pt_endcap_Rgamma0->Write();
    h1_sigmaEta_vs_pt_barrel_Rgamma5_20->Write();
    h1_sigmaEta_vs_pt_endcap_Rgamma5_20->Write();
    h1_sigmaEta_vs_pt_barrel_Rgamma30_60->Write();
    h1_sigmaEta_vs_pt_endcap_Rgamma30_60->Write();
    h1_sigmaEta_vs_pt_barrel_Rgamma60_100->Write();
    h1_sigmaEta_vs_pt_endcap_Rgamma60_100->Write();

    h1_sigmaPhi_vs_pt_barrel_Rgamma0->Write();
    h1_sigmaPhi_vs_pt_endcap_Rgamma0->Write();
    h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->Write();
    h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->Write();
    h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->Write();
    h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->Write();
    h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->Write();
    h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->Write();

  } //if PF


/*
  if( h1_response_vs_pt_FIT_barrel!=0 ) {
    delete h1_response_vs_pt_FIT_barrel;
    h1_response_vs_pt_FIT_barrel = 0;
  }
  if( h1_response_vs_pt_FIT_endcap!=0 ) {
    delete h1_response_vs_pt_FIT_endcap;
    h1_response_vs_pt_FIT_endcap = 0;
  }
  if( h1_response_vs_pt_FIT_forward!=0 ) {
    delete h1_response_vs_pt_FIT_forward;
    h1_response_vs_pt_FIT_forward = 0;
  }
*/
} //fillFile


