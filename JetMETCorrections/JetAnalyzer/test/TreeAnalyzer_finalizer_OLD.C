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


void fillFile(TFile* file, TTree* tree, std::string algoType, bool draw=false);
void fillProfile(TH1F* h1_response_FIT, TH1F* h1_resolution_FIT, TH1F* h1_response_MEAN, TH1F* h1_resolution_RMS, TH2D* h2, std::string name="");
void fillPositionResolution(TH1F* h1_sigmaEta, TH1F* h1_sigmaPhi, TH2D* h2_deltaEta, TH2D* h2_deltaPhi);
void fitProjection(TH1D* h1_projection, TF1* gaussian);


bool DEBUG_ = false;
bool DEBUG_PROVA_ = false;


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
  //std::cout << std::endl << std::endl << "**** starting calo SIS cone 5" << std::endl;
  //fillFile(outfile, chain, "caloSISCone5");
  //std::cout << std::endl << std::endl << "**** starting calo kt4" << std::endl;
  //fillFile(outfile, chain, "calokt4");
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
  int nBins = 21;
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
  Lower[20] = 600.;


  // histograms booking for all algorithms:
  histoName = "deltaEta_vs_pt_barrel_"+algoType;
  TH2D* h2_deltaEta_vs_pt_barrel = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_barrel = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType;
  TH2D* h2_deltaEta_vs_pt_endcap = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_endcap = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_forward_"+algoType;
  TH2D* h2_deltaEta_vs_pt_forward = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaEta_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_forward->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_forward_"+algoType;
  TH2D* h2_deltaPhi_vs_pt_forward = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_forward->SetYTitle("#Delta#phi");


  histoName = "response_vs_pt_barrel_"+algoType;
  TH2D* h2_response_vs_pt_barrel = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType;
  TH2D* h2_response_vs_pt_endcap = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_forward_"+algoType;
  TH2D* h2_response_vs_pt_forward = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta0_25"+algoType;
  TH2D* h2_response_vs_pt_eta0_25 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_25->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta14_25"+algoType;
  TH2D* h2_response_vs_pt_eta14_25 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_25->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_eta_"+algoType;
  TH2D* h2_response_vs_eta = new TH2D(histoName.c_str(), "Jet Response vs. True #eta", 100, -5., 5., 120, 0., 2.);
  h2_response_vs_eta->SetXTitle("#eta");
  h2_response_vs_eta->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");



  //additional histograms for PF only:
  histoName = "response_vs_pt_onlyCh_barrel_"+algoType;
  TH2D* h2_response_vs_pt_onlyCh_barrel = new TH2D(histoName.c_str(), "", nBins-2, Lower, 120, 0., 2.);
  
  histoName = "response_vs_pt_onlyCh_endcap_"+algoType;
  TH2D* h2_response_vs_pt_onlyCh_endcap = new TH2D(histoName.c_str(), "", nBins-2, Lower, 120, 0., 2.);
  


  histoName = "response_vs_Rch_"+algoType;
  TH2D* h2_response_vs_Rch_barrel = new TH2D(histoName.c_str(), "Jet Response vs. R_{ch}", 50, 0., 1., 120, 0., 2.);
  h2_response_vs_Rch_barrel->SetXTitle("R_{ch}");
  h2_response_vs_Rch_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_Rgamma_"+algoType;
  TH2D* h2_response_vs_Rgamma_barrel = new TH2D(histoName.c_str(), "Jet Response vs. R_{#gamma}", 50, 0., 1., 120, 0., 2.);
  h2_response_vs_Rgamma_barrel->SetXTitle("R_{#gamma}");
  h2_response_vs_Rgamma_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_Rnh_"+algoType;
  TH2D* h2_response_vs_Rnh_barrel = new TH2D(histoName.c_str(), "Jet Response vs. R_{nh}", 50, 0., 1., 120, 0., 2.);
  h2_response_vs_Rnh_barrel->SetXTitle("R_{nh}");
  h2_response_vs_Rnh_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");


  //--------------------- Rch --------------------------------------------

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch10_30->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch10_30->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch50_70->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch50_70->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_deltaEta_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rch90_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rch90_100->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch0->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch10_30->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch10_30->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch50_70->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch50_70->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_deltaEta_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rch90_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rch90_100->SetYTitle("#Delta#phi");


  //---------------------- Rgamma

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma5_20 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rgamma5_20->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma5_20 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma5_20->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma30_60 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rgamma30_60->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma30_60 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma30_60->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaEta_vs_pt_barrel_Rgamma60_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_Rgamma60_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaPhi_vs_pt_barrel_Rgamma60_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_Rgamma60_100->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma0 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma0->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma0 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma0->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma5_20 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma5_20->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma5_20 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma5_20->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma30_60 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma30_60->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma30_60 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma30_60->SetYTitle("#Delta#phi");

  histoName = "deltaEta_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaEta_vs_pt_endcap_Rgamma60_100 = new TH2D(histoName.c_str(), "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_Rgamma60_100->SetYTitle("#Delta#eta");

  histoName = "deltaPhi_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH2D* h2_deltaPhi_vs_pt_endcap_Rgamma60_100 = new TH2D(histoName.c_str(), "#phi^{RECO}-#phi^{GEN} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower,  100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_Rgamma60_100->SetYTitle("#Delta#phi");



  histoName = "response_vs_pt_barrel_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_barrel_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_barrel_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_barrel_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_barrel_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_endcap_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_endcap_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_endcap_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_endcap_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");


  //------ phillip plots
  histoName = "response_vs_pt_eta0_25_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_eta0_25_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta0_25_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_eta0_25_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta0_25_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_eta0_25_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta0_25_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_eta0_25_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta0_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta0_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta14_25_"+algoType+"_Rch0";
  TH2D* h2_response_vs_pt_eta14_25_Rch0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta14_25_"+algoType+"_Rch10_30";
  TH2D* h2_response_vs_pt_eta14_25_Rch10_30 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta14_25_"+algoType+"_Rch50_70";
  TH2D* h2_response_vs_pt_eta14_25_Rch50_70 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_eta14_25_"+algoType+"_Rch90_100";
  TH2D* h2_response_vs_pt_eta14_25_Rch90_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_eta14_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_eta14_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  //------


  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH2D* h2_response_vs_pt_barrel_Rgamma0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, R_{gamma}=0%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH2D* h2_response_vs_pt_barrel_Rgamma5_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 5<R_{gamma}<20%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH2D* h2_response_vs_pt_barrel_Rgamma30_60 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH2D* h2_response_vs_pt_barrel_Rgamma60_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH2D* h2_response_vs_pt_endcap_Rgamma0 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, R_{gamma}=0%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH2D* h2_response_vs_pt_endcap_Rgamma5_20 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 5<R_{gamma}<20%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH2D* h2_response_vs_pt_endcap_Rgamma30_60 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower, 120, 0., 2.);
  h2_response_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH2D* h2_response_vs_pt_endcap_Rgamma60_100 = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower, 120, 0., 2.);
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

  histoName = "Rch_vs_eta_50_80_"+algoType;
  TProfile* hp_Rch_vs_eta_50_80 = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rch_vs_eta_50_80->SetXTitle("#eta^{GEN}");
  hp_Rch_vs_eta_50_80->SetYTitle("entries");

  histoName = "Rch_vs_pt_"+algoType;
  TProfile* hp_Rch_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rch_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt->SetYTitle("entries");

  histoName = "Rch_vs_pt_barrel_"+algoType;
  TProfile* hp_Rch_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rch_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rch_vs_pt_endcap_"+algoType;
  TProfile* hp_Rch_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rch_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rch_vs_pt_endcap->SetYTitle("entries");


  histoName = "Rgamma_"+algoType;
  TH1F* h1_Rgamma = new TH1F(histoName.c_str(), "Jet Photon Energy Fraction", 100, 0., 1.1);
  h1_Rgamma->SetXTitle("R_{#gamma}");
  h1_Rgamma->SetYTitle("entries");

  histoName = "Rgamma_vs_pt_"+algoType;
  TProfile* hp_Rgamma_vs_pt = new TProfile(histoName.c_str(), "R_{#gamma} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rgamma_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rgamma_vs_pt->SetYTitle("entries");

  histoName = "Rgamma_vs_pt_barrel_"+algoType;
  TProfile* hp_Rgamma_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{#gamma} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rgamma_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rgamma_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rgamma_vs_pt_endcap_"+algoType;
  TProfile* hp_Rgamma_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{#gamma} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rgamma_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rgamma_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rgamma_vs_eta_"+algoType;
  TProfile* hp_Rgamma_vs_eta = new TProfile(histoName.c_str(), "R_{#gamma} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rgamma_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rgamma_vs_eta->SetYTitle("entries");

  histoName = "Rgamma_vs_eta_50_80_"+algoType;
  TProfile* hp_Rgamma_vs_eta_50_80 = new TProfile(histoName.c_str(), "R_{#gamma} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rgamma_vs_eta_50_80->SetXTitle("#eta^{GEN}");
  hp_Rgamma_vs_eta_50_80->SetYTitle("entries");

  histoName = "Rgamma_vs_eta_Rch0_"+algoType;
  TProfile* hp_Rgamma_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rgamma_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_Rgamma_vs_eta_Rch0->SetYTitle("entries");

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
  TProfile* hp_Rnh_vs_pt = new TProfile(histoName.c_str(), "R_{nh} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rnh_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rnh_vs_pt->SetYTitle("entries");

  histoName = "Rnh_vs_pt_barrel_"+algoType;
  TProfile* hp_Rnh_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{nh} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rnh_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rnh_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rnh_vs_pt_endcap_"+algoType;
  TProfile* hp_Rnh_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{nh} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rnh_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rnh_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rnh_vs_eta_"+algoType;
  TProfile* hp_Rnh_vs_eta = new TProfile(histoName.c_str(), "R_{nh} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rnh_vs_eta->SetYTitle("entries");

  histoName = "Rnh_vs_eta_50_80_"+algoType;
  TProfile* hp_Rnh_vs_eta_50_80 = new TProfile(histoName.c_str(), "R_{nh} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_vs_eta_50_80->SetXTitle("#eta^{GEN}");
  hp_Rnh_vs_eta_50_80->SetYTitle("entries");

  histoName = "Rnh_vs_eta_Rch0_"+algoType;
  TProfile* hp_Rnh_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_Rnh_vs_eta_Rch0->SetYTitle("entries");

  histoName = "Rnh_vs_Rch_"+algoType;
  TProfile* hp_Rnh_vs_Rch = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rnh_vs_Rch->SetXTitle("R_{ch}");
  hp_Rnh_vs_Rch->SetYTitle("R_{nh}");

  histoName = "Rnh_vs_Rch_barrel_"+algoType;
  TProfile* hp_Rnh_vs_Rch_barrel = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rnh_vs_Rch_barrel->SetXTitle("R_{ch}");
  hp_Rnh_vs_Rch_barrel->SetYTitle("R_{nh}");


 
  histoName = "Rmu_"+algoType;
  TH1F* h1_Rmu = new TH1F(histoName.c_str(), "Jet Muon Energy Fraction", 100, 0., 1.1);
  h1_Rmu->SetXTitle("R_{#gamma}");
  h1_Rmu->SetYTitle("entries");

  histoName = "Rmu_vs_pt_"+algoType;
  TProfile* hp_Rmu_vs_pt = new TProfile(histoName.c_str(), "R_{#mu} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rmu_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rmu_vs_pt->SetYTitle("entries");

  histoName = "Rmu_vs_pt_barrel_"+algoType;
  TProfile* hp_Rmu_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{#mu} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rmu_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rmu_vs_pt_barrel->SetYTitle("entries");

  histoName = "Rmu_vs_pt_endcap_"+algoType;
  TProfile* hp_Rmu_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{#mu} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rmu_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rmu_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rmu_vs_eta_"+algoType;
  TProfile* hp_Rmu_vs_eta = new TProfile(histoName.c_str(), "R_{#mu} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rmu_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rmu_vs_eta->SetYTitle("entries");

  histoName = "Rmu_vs_eta_50_80_"+algoType;
  TProfile* hp_Rmu_vs_eta_50_80 = new TProfile(histoName.c_str(), "R_{#mu} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rmu_vs_eta_50_80->SetXTitle("#eta^{GEN}");
  hp_Rmu_vs_eta_50_80->SetYTitle("entries");

  histoName = "Rmu_vs_eta_Rch0_"+algoType;
  TProfile* hp_Rmu_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rmu_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_Rmu_vs_eta_Rch0->SetYTitle("entries");

  histoName = "Rmu_vs_Rch_"+algoType;
  TProfile* hp_Rmu_vs_Rch = new TProfile(histoName.c_str(), "R_{#mu} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_Rmu_vs_Rch->SetXTitle("R_{ch}");
  hp_Rmu_vs_Rch->SetYTitle("R_{#mu}");


  histoName = "Rhf_"+algoType;
  TH1F* h1_Rhf = new TH1F(histoName.c_str(), "Jet HF Hadron Energy Fraction", 100, 0., 1.1);
  h1_Rhf->SetXTitle("R_{HF}");
  h1_Rhf->SetYTitle("entries");

  histoName = "Rhf_vs_pt_"+algoType;
  TProfile* hp_Rhf_vs_pt = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rhf_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhf_vs_pt->SetYTitle("entries");

  histoName = "Rhf_vs_pt_endcap_"+algoType;
  TProfile* hp_Rhf_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{hf} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_Rhf_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_Rhf_vs_pt_endcap->SetYTitle("entries");

  histoName = "Rhf_vs_eta_"+algoType;
  TProfile* hp_Rhf_vs_eta = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhf_vs_eta->SetXTitle("#eta^{GEN}");
  hp_Rhf_vs_eta->SetYTitle("entries");

  histoName = "Rhf_vs_eta_50_80_"+algoType;
  TProfile* hp_Rhf_vs_eta_50_80 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhf_vs_eta_50_80->SetXTitle("#eta^{GEN}");
  hp_Rhf_vs_eta_50_80->SetYTitle("entries");

  histoName = "Rhf_vs_eta_Rch0_"+algoType;
  TProfile* hp_Rhf_vs_eta_Rch0 = new TProfile(histoName.c_str(), "R_{hf} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_Rhf_vs_eta_Rch0->SetXTitle("#eta^{GEN}");
  hp_Rhf_vs_eta_Rch0->SetYTitle("entries");


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

  histoName = "RchGen_vs_pt_Rch0_"+algoType;
  TProfile* hp_RchGen_vs_pt_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RchGen_vs_pt_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt_Rch0->SetYTitle("entries");

  histoName = "RchGen_vs_pt_"+algoType;
  TProfile* hp_RchGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RchGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt->SetYTitle("entries");

  histoName = "RchGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RchGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RchGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RchGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RchGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RchGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
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

  histoName = "RgammaGen_vs_pt_Rch0_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RgammaGen_vs_pt_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt_Rch0->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_"+algoType;
  TProfile* hp_RgammaGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RgammaGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RgammaGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RgammaGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RgammaGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RgammaGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
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

  histoName = "RnhGen_vs_pt_Rch0_"+algoType;
  TProfile* hp_RnhGen_vs_pt_Rch0 = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RnhGen_vs_pt_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt_Rch0->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_"+algoType;
  TProfile* hp_RnhGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RnhGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RnhGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RnhGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RnhGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RnhGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RnhGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
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
  TProfile* hp_RmuGen_vs_pt = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RmuGen_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_barrel_"+algoType;
  TProfile* hp_RmuGen_vs_pt_barrel = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RmuGen_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt_barrel->SetYTitle("entries");

  histoName = "RmuGen_vs_pt_endcap_"+algoType;
  TProfile* hp_RmuGen_vs_pt_endcap = new TProfile(histoName.c_str(), "R_{ch} vs. pt", nBins-2, Lower, 0., 1000.);
  hp_RmuGen_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  hp_RmuGen_vs_pt_endcap->SetYTitle("entries");

  histoName = "RmuGen_vs_RchGen_"+algoType;
  TProfile* hp_RmuGen_vs_RchGen = new TProfile(histoName.c_str(), "R_{nh} vs. R_{ch}", 100, 0., 1.1, 0., 1000.);
  hp_RmuGen_vs_RchGen->SetXTitle("R_{ch}");
  hp_RmuGen_vs_RchGen->SetYTitle("R_{nh}");



  int nEntries = tree->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    for(int iRecoJet = 0; iRecoJet<nReco; ++iRecoJet) {

      if( iRecoJet >= 40 ) continue;

      Float_t iEtaReco = etaReco[iRecoJet];
      Float_t iPhiReco = phiReco[iRecoJet];
      Float_t iEtaGen = etaGen[iRecoJet];
      Float_t iPhiGen = phiGen[iRecoJet];
      Float_t iThetaReco = 2*atan(exp(-iEtaReco));

      Float_t iPtReco = ptReco[iRecoJet];
      Float_t iEReco = eReco[iRecoJet];
      Float_t iEtReco = iEReco*sin(iThetaReco);
      Float_t iPtGen = ptGen[iRecoJet];
      Float_t iEGen = eGen[iRecoJet];
      Float_t response = iPtReco/iPtGen;


      Float_t Rch=0.;
      Float_t Rgamma=0.;
      Float_t Rnh=0.;
      Float_t Rmu=0.;
      Float_t Rhf=0.;

      Float_t RchGen=0.;
      Float_t RgammaGen=0.;
      Float_t RnhGen=0.;
      Float_t RmuGen=0.;

      if( algoType == "PFItCone5") {

        Rch = eTracksReco[iRecoJet]/iEReco;
        Rgamma = ePhotonsReco[iRecoJet]/iEReco;
        Rnh = eNeutralHadronsReco[iRecoJet]/iEReco;
        Rmu = eMuonsReco[iRecoJet]/iEReco;
        Rhf = eHFhadronsReco[iRecoJet]/iEReco;

        RchGen = eTracksGen[iRecoJet]/iEGen;
        RgammaGen = ePhotonsGen[iRecoJet]/iEGen;
        RnhGen = eNeutralHadronsGen[iRecoJet]/iEGen;
        RmuGen = eMuonsGen[iRecoJet]/iEGen;

      }

      bool fillPFhistos=true;
/*
      if( (Rch>1.)||(Rgamma>1.)||(Rnh>1.)||(Rmu>1.)||(Rhf>1.) )
        fillPFhistos=false;
*/
      Float_t deltaEta = iEtaReco-iEtaGen;
      Float_t deltaPhi = iPhiReco-iPhiGen;
      Float_t pi = 3.14159;
      if(deltaPhi < -pi) deltaPhi += 2.*pi;
      if(deltaPhi > pi) deltaPhi -= 2.*pi;

      Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

      if( (iPtGen>15.)&&(iPtGen<300.) ) {
        if( (fabs(iEtaGen)>0.)&&(fabs(iEtaGen)<3.) ) h1_deltaR_barrelEndcap->Fill(deltaR);
        if( (fabs(iEtaGen)>3.)&&(fabs(iEtaGen)<5.) ) h1_deltaR_forward->Fill(deltaR);
        if( (fabs(iEtaGen)>4.)&&(fabs(iEtaGen)<5.) ) h1_deltaR_veryForward->Fill(deltaR);
      } //if

      bool matched = false;

      if( algoType=="PFItCone5") {
        if( deltaR<0.1 ) matched = true;
      } else {
        if( deltaR<0.3 ) matched = true;
      }

      if( matched ) {
    
        h2_response_vs_eta->Fill(iEtaGen, response);

        if( algoType=="PFItCone5" ) {

          hp_Rch_vs_eta->Fill(iEtaGen, Rch);
          hp_Rgamma_vs_eta->Fill(iEtaGen, Rgamma);
          hp_Rnh_vs_eta->Fill(iEtaGen, Rnh);
          hp_Rmu_vs_eta->Fill(iEtaGen, Rmu);
          hp_Rhf_vs_eta->Fill(iEtaGen, Rhf);

          if( (iPtGen>50.) && (iPtGen<80.) ) {
            hp_Rch_vs_eta_50_80->Fill(iEtaGen, Rch);
            hp_Rgamma_vs_eta_50_80->Fill(iEtaGen, Rgamma);
            hp_Rnh_vs_eta_50_80->Fill(iEtaGen, Rnh);
            hp_Rmu_vs_eta_50_80->Fill(iEtaGen, Rmu);
            hp_Rhf_vs_eta_50_80->Fill(iEtaGen, Rhf);
          }
          


          if( Rch==0. ) {
            hp_Rgamma_vs_eta_Rch0->Fill(iEtaGen, Rgamma);
            hp_Rnh_vs_eta_Rch0->Fill(iEtaGen, Rnh);
            hp_Rmu_vs_eta_Rch0->Fill(iEtaGen, Rmu);
            hp_Rhf_vs_eta_Rch0->Fill(iEtaGen, Rhf);

            hp_RchGen_vs_eta_Rch0->Fill(iEtaGen, RchGen);
            hp_RgammaGen_vs_eta_Rch0->Fill(iEtaGen, RgammaGen);
            hp_RnhGen_vs_eta_Rch0->Fill(iEtaGen, RnhGen);

            if( fabs(iEtaGen)<1.4 ) {
              hp_RchGen_vs_pt_Rch0->Fill(iPtGen, RchGen);
              hp_RgammaGen_vs_pt_Rch0->Fill(iPtGen, RgammaGen);
              hp_RnhGen_vs_pt_Rch0->Fill(iPtGen, RnhGen);
            }
          }
        }


        deltaEta = fabs(iEtaReco) - fabs(iEtaGen);
    
        if( fabs(iEtaGen)<1.4 ) {
          h2_deltaEta_vs_pt_barrel->Fill(iPtGen, deltaEta);
          h2_deltaPhi_vs_pt_barrel->Fill(iPtGen, deltaPhi);
          h2_response_vs_pt_barrel->Fill(iPtGen, response);
          h2_response_vs_pt_eta0_25->Fill(iPtGen, response);

          if( (algoType=="PFItCone5")&&(fillPFhistos) ) {

            h2_response_vs_pt_onlyCh_barrel->Fill(iPtGen, Rch);

            hp_Rgamma_vs_Rch_barrel->Fill(Rch, Rgamma);
            hp_Rnh_vs_Rch_barrel->Fill(Rch, Rnh);

            hp_Rgamma_vs_Rch->Fill(Rch, Rgamma);
            hp_Rnh_vs_Rch->Fill(Rnh, Rgamma);
            hp_Rmu_vs_Rch->Fill(Rmu, Rgamma);

            h1_Rch->Fill(Rch);
            hp_Rch_vs_pt->Fill(iPtGen, Rch);

            h1_Rgamma->Fill(Rgamma);
            hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma);

            h1_Rnh->Fill(Rnh);
            hp_Rnh_vs_pt->Fill(iPtGen, Rnh);

            h1_Rmu->Fill(Rmu);
            hp_Rmu_vs_pt->Fill(iPtGen, Rmu);

            h1_Rhf->Fill(Rhf);
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
              h2_response_vs_pt_eta0_25_Rch0->Fill(iPtGen, response);
            } else if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
              h2_deltaEta_vs_pt_barrel_Rch10_30->Fill(iPtGen, deltaEta);
              h2_deltaPhi_vs_pt_barrel_Rch10_30->Fill(iPtGen, deltaPhi);
              h2_response_vs_pt_barrel_Rch10_30->Fill(iPtGen, response);
              h2_response_vs_pt_eta0_25_Rch10_30->Fill(iPtGen, response);
            } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
              h2_deltaEta_vs_pt_barrel_Rch50_70->Fill(iPtGen, deltaEta);
              h2_deltaPhi_vs_pt_barrel_Rch50_70->Fill(iPtGen, deltaPhi);
              h2_response_vs_pt_barrel_Rch50_70->Fill(iPtGen, response);
              h2_response_vs_pt_eta0_25_Rch50_70->Fill(iPtGen, response);
            } else if( (Rch >= 0.9) && (Rch<=1.) )  {
              h2_deltaEta_vs_pt_barrel_Rch90_100->Fill(iPtGen, deltaEta);
              h2_deltaPhi_vs_pt_barrel_Rch90_100->Fill(iPtGen, deltaPhi);
              h2_response_vs_pt_barrel_Rch90_100->Fill(iPtGen, response);
              h2_response_vs_pt_eta0_25_Rch90_100->Fill(iPtGen, response);
            }

            if( Rgamma==0. ) { 
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


            if( (iPtGen>30.)&&(iPtGen<40.) ) {
              h2_response_vs_Rch_barrel->Fill(Rch, response);
              h2_response_vs_Rgamma_barrel->Fill(Rgamma, response);
              h2_response_vs_Rnh_barrel->Fill(Rnh, response);
            }

          } //if PF 

        } //if barrel

        if( (fabs(iEtaGen)>1.4) && (fabs(iEtaGen)<3.) ) {
          h2_deltaEta_vs_pt_endcap->Fill(iPtGen, deltaEta);
          h2_deltaPhi_vs_pt_endcap->Fill(iPtGen, deltaPhi);
          h2_response_vs_pt_endcap->Fill(iPtGen, response);
          if( fabs(iEtaGen)<2.5 ) {
            h2_response_vs_pt_eta0_25->Fill(iPtGen, response);
            h2_response_vs_pt_eta14_25->Fill(iPtGen, response);
          }
         

          if( (algoType=="PFItCone5")&&(fillPFhistos) ) {

            if(fabs(iEtaGen)<2.5) {//fiducial tracker region
              hp_Rgamma_vs_Rch->Fill(Rch, Rgamma);
              hp_Rnh_vs_Rch->Fill(Rnh, Rgamma);
              hp_Rmu_vs_Rch->Fill(Rmu, Rgamma);
            }

            h2_response_vs_pt_onlyCh_endcap->Fill(iPtGen, Rch);

            h1_Rch->Fill(Rch);
            hp_Rch_vs_pt->Fill(iPtGen, Rch);

            h1_Rgamma->Fill(Rgamma);
            hp_Rgamma_vs_pt->Fill(iPtGen, Rgamma);

            h1_Rnh->Fill(Rnh);
            hp_Rnh_vs_pt->Fill(iPtGen, Rnh);

            h1_Rmu->Fill(Rmu);
            hp_Rmu_vs_pt->Fill(iPtGen, Rmu);

            h1_Rhf->Fill(Rhf);
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
              if( fabs(iEtaGen)<2.5 ) {
                 h2_response_vs_pt_eta0_25_Rch0->Fill(iPtGen, response);
                h2_response_vs_pt_eta14_25_Rch0->Fill(iPtGen, response);
              }
            } else if( (Rch >= 0.1) &&( Rch <= 0.3 ) ) {
              h2_deltaEta_vs_pt_endcap_Rch10_30->Fill(iPtGen, deltaEta);
              h2_deltaPhi_vs_pt_endcap_Rch10_30->Fill(iPtGen, deltaPhi);
              h2_response_vs_pt_endcap_Rch10_30->Fill(iPtGen, response);
              if( fabs(iEtaGen)<2.5 ) {
                 h2_response_vs_pt_eta0_25_Rch10_30->Fill(iPtGen, response);
                h2_response_vs_pt_eta14_25_Rch10_30->Fill(iPtGen, response);
              }
            } else if( (Rch >= 0.5) &&( Rch <= 0.7 ) ) {
              h2_deltaEta_vs_pt_endcap_Rch50_70->Fill(iPtGen, deltaEta);
              h2_deltaPhi_vs_pt_endcap_Rch50_70->Fill(iPtGen, deltaPhi);
              h2_response_vs_pt_endcap_Rch50_70->Fill(iPtGen, response);
              if( fabs(iEtaGen)<2.5 ) {
                 h2_response_vs_pt_eta0_25_Rch50_70->Fill(iPtGen, response);
                h2_response_vs_pt_eta14_25_Rch50_70->Fill(iPtGen, response);
              }
            } else if( (Rch >= 0.9) && (Rch <=1.) )  {
              h2_deltaEta_vs_pt_endcap_Rch90_100->Fill(iPtGen, deltaEta);
              h2_deltaPhi_vs_pt_endcap_Rch90_100->Fill(iPtGen, deltaPhi);
              h2_response_vs_pt_endcap_Rch90_100->Fill(iPtGen, response);
              if( fabs(iEtaGen)<2.5 ) {
                 h2_response_vs_pt_eta0_25_Rch90_100->Fill(iPtGen, response);
                h2_response_vs_pt_eta14_25_Rch90_100->Fill(iPtGen, response);
              }
            }

            if( Rgamma==0. ) { 
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


          } //if PF 

        } //if endcaps

        if( (fabs(iEtaGen)>3.) && (fabs(iEtaGen)<5.) ) {
          h2_deltaEta_vs_pt_forward->Fill(iPtGen, deltaEta);
          h2_deltaPhi_vs_pt_forward->Fill(iPtGen, deltaPhi);
          h2_response_vs_pt_forward->Fill(iPtGen, response);

          if( (algoType=="PFItCone5")&&(fillPFhistos) ) {

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

    } //if matched

  } //for entries

 
  histoName = "response_vs_pt_FIT_barrel_"+algoType;
  TH1F* h1_response_vs_pt_FIT_barrel = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_endcap_"+algoType;
  TH1F* h1_response_vs_pt_FIT_endcap = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_forward_"+algoType;
  TH1F* h1_response_vs_pt_FIT_forward = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_FIT_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_forward->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_eta0_25_"+algoType;
  TH1F* h1_response_vs_pt_FIT_eta0_25 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta0_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta0_25->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta0_25->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_eta14_25_"+algoType;
  TH1F* h1_response_vs_pt_FIT_eta14_25 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta14_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta14_25->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta14_25->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_pt_FIT_barrel_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_endcap_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_endcap = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_endcap->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_forward_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_forward = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_forward->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_forward->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta0_25_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_eta0_25 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta0_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta0_25->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_eta0_25->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta14_25_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_eta14_25 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta14_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta14_25->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_eta14_25->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_barrel = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_endcap_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_endcap = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_forward_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_forward = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_forward->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_forward->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_eta0_25_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_eta0_25 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta0_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta0_25->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta0_25->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_eta14_25_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_eta14_25 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta14_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta14_25->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta14_25->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_pt_RMS_barrel_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_endcap_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_endcap = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_forward_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_forward = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_forward->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_forward->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta0_25_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_eta0_25 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta0_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta0_25->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_eta0_25->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta14_25_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_eta14_25 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta14_25->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta14_25->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_eta14_25->GetYaxis()->SetRangeUser(0., 0.5);


  if(draw) {
    std::string name;
    name = algoType+"_barrel";
    fillProfile(h1_response_vs_pt_FIT_barrel, h1_resolution_vs_pt_FIT_barrel, h1_response_vs_pt_MEAN_barrel, h1_resolution_vs_pt_RMS_barrel, h2_response_vs_pt_barrel, name);
    name = algoType+"_endcap";
    fillProfile(h1_response_vs_pt_FIT_endcap, h1_resolution_vs_pt_FIT_endcap, h1_response_vs_pt_MEAN_endcap, h1_resolution_vs_pt_RMS_endcap, h2_response_vs_pt_endcap, name);
    name = algoType+"_forward";
    fillProfile(h1_response_vs_pt_FIT_forward, h1_resolution_vs_pt_FIT_forward, h1_response_vs_pt_MEAN_forward, h1_resolution_vs_pt_RMS_forward, h2_response_vs_pt_forward, name);
    name = algoType+"_eta0_25";
    fillProfile(h1_response_vs_pt_FIT_eta0_25, h1_resolution_vs_pt_FIT_eta0_25, h1_response_vs_pt_MEAN_eta0_25, h1_resolution_vs_pt_RMS_eta0_25, h2_response_vs_pt_eta0_25, name);
    name = algoType+"_eta14_25";
    fillProfile(h1_response_vs_pt_FIT_eta14_25, h1_resolution_vs_pt_FIT_eta14_25, h1_response_vs_pt_MEAN_eta14_25, h1_resolution_vs_pt_RMS_eta14_25, h2_response_vs_pt_eta14_25, name);
  } else { 
    fillProfile(h1_response_vs_pt_FIT_barrel, h1_resolution_vs_pt_FIT_barrel, h1_response_vs_pt_MEAN_barrel, h1_resolution_vs_pt_RMS_barrel, h2_response_vs_pt_barrel);
    fillProfile(h1_response_vs_pt_FIT_endcap, h1_resolution_vs_pt_FIT_endcap, h1_response_vs_pt_MEAN_endcap, h1_resolution_vs_pt_RMS_endcap, h2_response_vs_pt_endcap);
    fillProfile(h1_response_vs_pt_FIT_forward, h1_resolution_vs_pt_FIT_forward, h1_response_vs_pt_MEAN_forward, h1_resolution_vs_pt_RMS_forward, h2_response_vs_pt_forward);
    fillProfile(h1_response_vs_pt_FIT_eta0_25, h1_resolution_vs_pt_FIT_eta0_25, h1_response_vs_pt_MEAN_eta0_25, h1_resolution_vs_pt_RMS_eta0_25, h2_response_vs_pt_eta0_25);
    fillProfile(h1_response_vs_pt_FIT_eta14_25, h1_resolution_vs_pt_FIT_eta14_25, h1_response_vs_pt_MEAN_eta14_25, h1_resolution_vs_pt_RMS_eta14_25, h2_response_vs_pt_eta14_25);
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
  histoName = "response_vs_pt_onlyCh_MEAN_barrel_"+algoType;
  TH1F* h1_response_vs_pt_onlyCh_MEAN_barrel = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_onlyCh_MEAN_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_onlyCh_MEAN_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_onlyCh_MEAN_barrel->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_onlyCh_FIT_barrel_"+algoType;
  TH1F* h1_response_vs_pt_onlyCh_FIT_barrel = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_onlyCh_FIT_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_onlyCh_FIT_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_onlyCh_FIT_barrel->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "resolution_vs_pt_onlyCh_RMS_barrel_"+algoType;
  TH1F* h1_resolution_vs_pt_onlyCh_RMS_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_onlyCh_RMS_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_onlyCh_RMS_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_onlyCh_RMS_barrel->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_onlyCh_FIT_barrel_"+algoType;
  TH1F* h1_resolution_vs_pt_onlyCh_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_onlyCh_FIT_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_onlyCh_FIT_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_onlyCh_FIT_barrel->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "response_vs_pt_onlyCh_MEAN_endcap_"+algoType;
  TH1F* h1_response_vs_pt_onlyCh_MEAN_endcap = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_onlyCh_MEAN_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_onlyCh_MEAN_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_onlyCh_MEAN_endcap->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_onlyCh_FIT_endcap_"+algoType;
  TH1F* h1_response_vs_pt_onlyCh_FIT_endcap = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_onlyCh_FIT_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_onlyCh_FIT_endcap->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_onlyCh_FIT_endcap->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "resolution_vs_pt_onlyCh_RMS_endcap_"+algoType;
  TH1F* h1_resolution_vs_pt_onlyCh_RMS_endcap = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_onlyCh_RMS_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_onlyCh_RMS_endcap->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_onlyCh_RMS_endcap->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_onlyCh_FIT_endcap_"+algoType;
  TH1F* h1_resolution_vs_pt_onlyCh_FIT_endcap = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_onlyCh_FIT_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_onlyCh_FIT_endcap->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_onlyCh_FIT_endcap->GetYaxis()->SetRangeUser(0., 0.5);
 

  //Rch: response vs pt FIT

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_FIT_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_FIT_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_eta0_25_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_FIT_eta0_25_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta0_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta0_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta0_25_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_eta0_25_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_FIT_eta0_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta0_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta0_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta0_25_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_eta0_25_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_FIT_eta0_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta0_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta0_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta0_25_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_eta0_25_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_FIT_eta0_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta0_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta0_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta0_25_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_eta14_25_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_FIT_eta14_25_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta14_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta14_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta14_25_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_eta14_25_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_FIT_eta14_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta14_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta14_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta14_25_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_FIT_eta14_25_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_FIT_eta14_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta14_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta14_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta14_25_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_eta14_25_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_FIT_eta14_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_eta14_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_eta14_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_eta14_25_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 

  // Rch: resolution vs pt FIT

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta0_25_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_FIT_eta0_25_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta0_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta0_25_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_FIT_eta0_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta0_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta0_25_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_FIT_eta0_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta0_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta0_25_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_FIT_eta0_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta0_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta0_25_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta14_25_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_FIT_eta14_25_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta14_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta14_25_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_FIT_eta14_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta14_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta14_25_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_FIT_eta14_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta14_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_eta14_25_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_FIT_eta14_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_eta14_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_eta14_25_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);


  // Rch: response vs pt mean

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_eta0_25_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_MEAN_eta0_25_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta0_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta0_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta0_25_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_eta0_25_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_MEAN_eta0_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta0_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta0_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta0_25_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_eta0_25_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_MEAN_eta0_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta0_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta0_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta0_25_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_eta0_25_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_MEAN_eta0_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta0_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta0_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta0_25_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_eta14_25_"+algoType+"_Rch0";
  TH1F* h1_response_vs_pt_MEAN_eta14_25_Rch0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta14_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta14_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta14_25_Rch0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_eta14_25_"+algoType+"_Rch10_30";
  TH1F* h1_response_vs_pt_MEAN_eta14_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta14_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta14_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta14_25_Rch10_30->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "response_vs_pt_MEAN_eta14_25_"+algoType+"_Rch50_70";
  TH1F* h1_response_vs_pt_MEAN_eta14_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta14_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta14_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta14_25_Rch50_70->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_eta14_25_"+algoType+"_Rch90_100";
  TH1F* h1_response_vs_pt_MEAN_eta14_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_eta14_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_eta14_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_eta14_25_Rch90_100->GetYaxis()->SetRangeUser(0.2, 1.);
 

  //Rch: resolution vs pt RMS

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch10_30->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch50_70->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rch90_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch10_30->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch50_70->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rch90_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta0_25_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_RMS_eta0_25_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta0_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta0_25_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_RMS_eta0_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta0_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta0_25_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_RMS_eta0_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta0_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta0_25_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_RMS_eta0_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta0_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta0_25_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta14_25_"+algoType+"_Rch0";
  TH1F* h1_resolution_vs_pt_RMS_eta14_25_Rch0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta14_25_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta14_25_"+algoType+"_Rch10_30";
  TH1F* h1_resolution_vs_pt_RMS_eta14_25_Rch10_30 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta14_25_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch10_30->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch10_30->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta14_25_"+algoType+"_Rch50_70";
  TH1F* h1_resolution_vs_pt_RMS_eta14_25_Rch50_70 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta14_25_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch50_70->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch50_70->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_eta14_25_"+algoType+"_Rch90_100";
  TH1F* h1_resolution_vs_pt_RMS_eta14_25_Rch90_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_eta14_25_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch90_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_RMS_eta14_25_Rch90_100->GetYaxis()->SetRangeUser(0., 0.5);



  // Rgamma

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_FIT_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_FIT_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_FIT_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_FIT_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_FIT_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN} [%]");
  h1_resolution_vs_pt_FIT_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_MEAN_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma0->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma0->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma5_20->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma30_60->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_pt_MEAN_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_response_vs_pt_MEAN_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_Rgamma60_100->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma5_20->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma30_60->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_RMS_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma60_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{#gamma}<20%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma0->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma5_20->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{#gamma}<50%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma30_60->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_resolution_vs_pt_RMS_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{#gamma}<100%)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma60_100->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.5);




  histoName = "response_vs_Rch_FIT_barrel_"+algoType;
  TH1F* h1_response_vs_Rch_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{ch}", 50, 0., 1.);
  h1_response_vs_Rch_FIT_barrel->SetXTitle("R_{ch}");
  h1_response_vs_Rch_FIT_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rch_FIT_barrel->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_Rgamma_FIT_barrel_"+algoType;
  TH1F* h1_response_vs_Rgamma_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{gamma}", 50, 0., 1.);
  h1_response_vs_Rgamma_FIT_barrel->SetXTitle("R_{#gamma}");
  h1_response_vs_Rgamma_FIT_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rgamma_FIT_barrel->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_Rnh_FIT_barrel_"+algoType;
  TH1F* h1_response_vs_Rnh_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{nh}", 50, 0., 1.);
  h1_response_vs_Rnh_FIT_barrel->SetXTitle("R_{nh}");
  h1_response_vs_Rnh_FIT_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rnh_FIT_barrel->GetYaxis()->SetRangeUser(0.2, 1.);


  histoName = "response_vs_Rch_MEAN_barrel_"+algoType;
  TH1F* h1_response_vs_Rch_MEAN_barrel = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{ch}", 50, 0., 1.);
  h1_response_vs_Rch_MEAN_barrel->SetXTitle("R_{ch}");
  h1_response_vs_Rch_MEAN_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rch_MEAN_barrel->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_Rgamma_MEAN_barrel_"+algoType;
  TH1F* h1_response_vs_Rgamma_MEAN_barrel = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{gamma}", 50, 0., 1.);
  h1_response_vs_Rgamma_MEAN_barrel->SetXTitle("R_{#gamma}");
  h1_response_vs_Rgamma_MEAN_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rgamma_MEAN_barrel->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "response_vs_Rnh_MEAN_barrel_"+algoType;
  TH1F* h1_response_vs_Rnh_MEAN_barrel = new TH1F(histoName.c_str(), "Jet Energy Response vs. R_{nh}", 50, 0., 1.);
  h1_response_vs_Rnh_MEAN_barrel->SetXTitle("R_{nh}");
  h1_response_vs_Rnh_MEAN_barrel->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_Rnh_MEAN_barrel->GetYaxis()->SetRangeUser(0.2, 1.);

  histoName = "resolution_vs_Rch_FIT_barrel_"+algoType;
  TH1F* h1_resolution_vs_Rch_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{ch}", 50, 0., 1.);
  h1_resolution_vs_Rch_FIT_barrel->SetXTitle("R_{ch}");
  h1_resolution_vs_Rch_FIT_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rch_FIT_barrel->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_Rgamma_FIT_barrel_"+algoType;
  TH1F* h1_resolution_vs_Rgamma_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{#gamma}", 50, 0., 1.);
  h1_resolution_vs_Rgamma_FIT_barrel->SetXTitle("R_{#gamma}");
  h1_resolution_vs_Rgamma_FIT_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rgamma_FIT_barrel->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_Rnh_FIT_barrel_"+algoType;
  TH1F* h1_resolution_vs_Rnh_FIT_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{nh}", 50, 0., 1.);
  h1_resolution_vs_Rnh_FIT_barrel->SetXTitle("R_{nh}");
  h1_resolution_vs_Rnh_FIT_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rnh_FIT_barrel->GetYaxis()->SetRangeUser(0., 0.5);


  histoName = "resolution_vs_Rch_RMS_barrel_"+algoType;
  TH1F* h1_resolution_vs_Rch_RMS_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{ch}", 50, 0., 1.);
  h1_resolution_vs_Rch_RMS_barrel->SetXTitle("R_{ch}");
  h1_resolution_vs_Rch_RMS_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rch_RMS_barrel->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_Rgamma_RMS_barrel_"+algoType;
  TH1F* h1_resolution_vs_Rgamma_RMS_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{#gamma}", 50, 0., 1.);
  h1_resolution_vs_Rgamma_RMS_barrel->SetXTitle("R_{#gamma}");
  h1_resolution_vs_Rgamma_RMS_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rgamma_RMS_barrel->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_Rnh_RMS_barrel_"+algoType;
  TH1F* h1_resolution_vs_Rnh_RMS_barrel = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. R_{nh}", 50, 0., 1.);
  h1_resolution_vs_Rnh_RMS_barrel->SetXTitle("R_{nh}");
  h1_resolution_vs_Rnh_RMS_barrel->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_Rnh_RMS_barrel->GetYaxis()->SetRangeUser(0., 0.5);


  if( algoType == "PFItCone5" ) { 
    if(draw) {
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch0, h1_resolution_vs_pt_FIT_barrel_Rch0, h1_response_vs_pt_MEAN_barrel_Rch0, h1_resolution_vs_pt_RMS_barrel_Rch0, h2_response_vs_pt_barrel_Rch0, "Rch0_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch0, h1_resolution_vs_pt_FIT_endcap_Rch0, h1_response_vs_pt_MEAN_endcap_Rch0, h1_resolution_vs_pt_RMS_endcap_Rch0, h2_response_vs_pt_endcap_Rch0, "Rch0_endcap");
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch10_30, h1_resolution_vs_pt_FIT_barrel_Rch10_30, h1_response_vs_pt_MEAN_barrel_Rch10_30, h1_resolution_vs_pt_RMS_barrel_Rch10_30, h2_response_vs_pt_barrel_Rch10_30, "Rch10_30_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch10_30, h1_resolution_vs_pt_FIT_endcap_Rch10_30, h1_response_vs_pt_MEAN_endcap_Rch10_30, h1_resolution_vs_pt_RMS_endcap_Rch10_30, h2_response_vs_pt_endcap_Rch10_30, "Rch10_30_endcap");
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch50_70, h1_resolution_vs_pt_FIT_barrel_Rch50_70, h1_response_vs_pt_MEAN_barrel_Rch50_70, h1_resolution_vs_pt_RMS_barrel_Rch50_70, h2_response_vs_pt_barrel_Rch50_70, "Rch50_70_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch50_70, h1_resolution_vs_pt_FIT_endcap_Rch50_70, h1_response_vs_pt_MEAN_endcap_Rch50_70, h1_resolution_vs_pt_RMS_endcap_Rch50_70, h2_response_vs_pt_endcap_Rch50_70, "Rch50_70_endcap");
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch90_100, h1_resolution_vs_pt_FIT_barrel_Rch90_100, h1_response_vs_pt_MEAN_barrel_Rch90_100, h1_resolution_vs_pt_RMS_barrel_Rch90_100, h2_response_vs_pt_barrel_Rch90_100, "Rch90_100_barrel");
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch90_100, h1_resolution_vs_pt_FIT_endcap_Rch90_100, h1_response_vs_pt_MEAN_endcap_Rch90_100, h1_resolution_vs_pt_RMS_endcap_Rch90_100, h2_response_vs_pt_endcap_Rch90_100, "Rch90_100_endcap");

      fillProfile(h1_response_vs_pt_FIT_eta0_25_Rch0, h1_resolution_vs_pt_FIT_eta0_25_Rch0, h1_response_vs_pt_MEAN_eta0_25_Rch0, h1_resolution_vs_pt_RMS_eta0_25_Rch0, h2_response_vs_pt_eta0_25_Rch0, "Rch0_eta0_25");
      fillProfile(h1_response_vs_pt_FIT_eta14_25_Rch0, h1_resolution_vs_pt_FIT_eta14_25_Rch0, h1_response_vs_pt_MEAN_eta14_25_Rch0, h1_resolution_vs_pt_RMS_eta14_25_Rch0, h2_response_vs_pt_eta14_25_Rch0, "Rch0_eta14_25");
      fillProfile(h1_response_vs_pt_FIT_eta0_25_Rch10_30, h1_resolution_vs_pt_FIT_eta0_25_Rch10_30, h1_response_vs_pt_MEAN_eta0_25_Rch10_30, h1_resolution_vs_pt_RMS_eta0_25_Rch10_30, h2_response_vs_pt_eta0_25_Rch10_30, "Rch10_30_eta0_25");
      fillProfile(h1_response_vs_pt_FIT_eta14_25_Rch10_30, h1_resolution_vs_pt_FIT_eta14_25_Rch10_30, h1_response_vs_pt_MEAN_eta14_25_Rch10_30, h1_resolution_vs_pt_RMS_eta14_25_Rch10_30, h2_response_vs_pt_eta14_25_Rch10_30, "Rch10_30_eta14_25");
      fillProfile(h1_response_vs_pt_FIT_eta0_25_Rch50_70, h1_resolution_vs_pt_FIT_eta0_25_Rch50_70, h1_response_vs_pt_MEAN_eta0_25_Rch50_70, h1_resolution_vs_pt_RMS_eta0_25_Rch50_70, h2_response_vs_pt_eta0_25_Rch50_70, "Rch50_70_eta0_25");
      fillProfile(h1_response_vs_pt_FIT_eta14_25_Rch50_70, h1_resolution_vs_pt_FIT_eta14_25_Rch50_70, h1_response_vs_pt_MEAN_eta14_25_Rch50_70, h1_resolution_vs_pt_RMS_eta14_25_Rch50_70, h2_response_vs_pt_eta14_25_Rch50_70, "Rch50_70_eta14_25");
      fillProfile(h1_response_vs_pt_FIT_eta0_25_Rch90_100, h1_resolution_vs_pt_FIT_eta0_25_Rch90_100, h1_response_vs_pt_MEAN_eta0_25_Rch90_100, h1_resolution_vs_pt_RMS_eta0_25_Rch90_100, h2_response_vs_pt_eta0_25_Rch90_100, "Rch90_100_eta0_25");
      fillProfile(h1_response_vs_pt_FIT_eta14_25_Rch90_100, h1_resolution_vs_pt_FIT_eta14_25_Rch90_100, h1_response_vs_pt_MEAN_eta14_25_Rch90_100, h1_resolution_vs_pt_RMS_eta14_25_Rch90_100, h2_response_vs_pt_eta14_25_Rch90_100, "Rch90_100_eta14_25");
    } else {
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch0, h1_resolution_vs_pt_FIT_barrel_Rch0, h1_response_vs_pt_MEAN_barrel_Rch0, h1_resolution_vs_pt_RMS_barrel_Rch0, h2_response_vs_pt_barrel_Rch0);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch0, h1_resolution_vs_pt_FIT_endcap_Rch0, h1_response_vs_pt_MEAN_endcap_Rch0, h1_resolution_vs_pt_RMS_endcap_Rch0, h2_response_vs_pt_endcap_Rch0);
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch10_30, h1_resolution_vs_pt_FIT_barrel_Rch10_30, h1_response_vs_pt_MEAN_barrel_Rch10_30, h1_resolution_vs_pt_RMS_barrel_Rch10_30, h2_response_vs_pt_barrel_Rch10_30);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch10_30, h1_resolution_vs_pt_FIT_endcap_Rch10_30, h1_response_vs_pt_MEAN_endcap_Rch10_30, h1_resolution_vs_pt_RMS_endcap_Rch10_30, h2_response_vs_pt_endcap_Rch10_30);
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch50_70, h1_resolution_vs_pt_FIT_barrel_Rch50_70, h1_response_vs_pt_MEAN_barrel_Rch50_70, h1_resolution_vs_pt_RMS_barrel_Rch50_70, h2_response_vs_pt_barrel_Rch50_70);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch50_70, h1_resolution_vs_pt_FIT_endcap_Rch50_70, h1_response_vs_pt_MEAN_endcap_Rch50_70, h1_resolution_vs_pt_RMS_endcap_Rch50_70, h2_response_vs_pt_endcap_Rch50_70);
      fillProfile(h1_response_vs_pt_FIT_barrel_Rch90_100, h1_resolution_vs_pt_FIT_barrel_Rch90_100, h1_response_vs_pt_MEAN_barrel_Rch90_100, h1_resolution_vs_pt_RMS_barrel_Rch90_100, h2_response_vs_pt_barrel_Rch90_100);
      fillProfile(h1_response_vs_pt_FIT_endcap_Rch90_100, h1_resolution_vs_pt_FIT_endcap_Rch90_100, h1_response_vs_pt_MEAN_endcap_Rch90_100, h1_resolution_vs_pt_RMS_endcap_Rch90_100, h2_response_vs_pt_endcap_Rch90_100);
    }

    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma0, h1_resolution_vs_pt_FIT_barrel_Rgamma0, h1_response_vs_pt_MEAN_barrel_Rgamma0, h1_resolution_vs_pt_RMS_barrel_Rgamma0, h2_response_vs_pt_barrel_Rgamma0);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma0, h1_resolution_vs_pt_FIT_endcap_Rgamma0, h1_response_vs_pt_MEAN_endcap_Rgamma0, h1_resolution_vs_pt_RMS_endcap_Rgamma0, h2_response_vs_pt_endcap_Rgamma0);
    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma5_20, h1_resolution_vs_pt_FIT_barrel_Rgamma5_20, h1_response_vs_pt_MEAN_barrel_Rgamma5_20, h1_resolution_vs_pt_RMS_barrel_Rgamma5_20, h2_response_vs_pt_barrel_Rgamma5_20);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma5_20, h1_resolution_vs_pt_FIT_endcap_Rgamma5_20, h1_response_vs_pt_MEAN_endcap_Rgamma5_20, h1_resolution_vs_pt_RMS_endcap_Rgamma5_20, h2_response_vs_pt_endcap_Rgamma5_20);
    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma30_60, h1_resolution_vs_pt_FIT_barrel_Rgamma30_60, h1_response_vs_pt_MEAN_barrel_Rgamma30_60, h1_resolution_vs_pt_RMS_barrel_Rgamma30_60, h2_response_vs_pt_barrel_Rgamma30_60);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma30_60, h1_resolution_vs_pt_FIT_endcap_Rgamma30_60, h1_response_vs_pt_MEAN_endcap_Rgamma30_60, h1_resolution_vs_pt_RMS_endcap_Rgamma30_60, h2_response_vs_pt_endcap_Rgamma30_60);
    fillProfile(h1_response_vs_pt_FIT_barrel_Rgamma60_100, h1_resolution_vs_pt_FIT_barrel_Rgamma60_100, h1_response_vs_pt_MEAN_barrel_Rgamma60_100, h1_resolution_vs_pt_RMS_barrel_Rgamma60_100, h2_response_vs_pt_barrel_Rgamma60_100);
    fillProfile(h1_response_vs_pt_FIT_endcap_Rgamma60_100, h1_resolution_vs_pt_FIT_endcap_Rgamma60_100, h1_response_vs_pt_MEAN_endcap_Rgamma60_100, h1_resolution_vs_pt_RMS_endcap_Rgamma60_100, h2_response_vs_pt_endcap_Rgamma60_100);

    fillProfile(h1_response_vs_Rch_FIT_barrel, h1_resolution_vs_Rch_FIT_barrel, h1_response_vs_Rch_MEAN_barrel, h1_resolution_vs_Rch_RMS_barrel, h2_response_vs_Rch_barrel);
    fillProfile(h1_response_vs_Rgamma_FIT_barrel, h1_resolution_vs_Rgamma_FIT_barrel, h1_response_vs_Rgamma_MEAN_barrel, h1_resolution_vs_Rgamma_RMS_barrel, h2_response_vs_Rgamma_barrel);
    fillProfile(h1_response_vs_Rnh_FIT_barrel, h1_resolution_vs_Rnh_FIT_barrel, h1_response_vs_Rnh_MEAN_barrel, h1_resolution_vs_Rnh_RMS_barrel, h2_response_vs_Rnh_barrel);

    fillProfile(h1_response_vs_pt_onlyCh_FIT_barrel, h1_resolution_vs_pt_onlyCh_FIT_barrel, h1_response_vs_pt_onlyCh_MEAN_barrel, h1_resolution_vs_pt_onlyCh_RMS_barrel, h2_response_vs_pt_onlyCh_barrel);
    fillProfile(h1_response_vs_pt_onlyCh_FIT_endcap, h1_resolution_vs_pt_onlyCh_FIT_endcap, h1_response_vs_pt_onlyCh_MEAN_endcap, h1_resolution_vs_pt_onlyCh_RMS_endcap, h2_response_vs_pt_onlyCh_endcap);

  }
  
  histoName = "sigmaEta_vs_pt_barrel_"+algoType;
  TH1F* h1_sigmaEta_vs_pt_barrel = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType;
  TH1F* h1_sigmaEta_vs_pt_endcap = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_forward_"+algoType;
  TH1F* h1_sigmaEta_vs_pt_forward = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_forward->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_forward->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_barrel_"+algoType;
  TH1F* h1_sigmaPhi_vs_pt_barrel = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_barrel->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType;
  TH1F* h1_sigmaPhi_vs_pt_endcap = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_endcap->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_forward_"+algoType;
  TH1F* h1_sigmaPhi_vs_pt_forward = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_forward->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_forward->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_forward->GetYaxis()->SetRangeUser(0., 0.2);

  fillPositionResolution(h1_sigmaEta_vs_pt_barrel, h1_sigmaPhi_vs_pt_barrel, h2_deltaEta_vs_pt_barrel, h2_deltaPhi_vs_pt_barrel);
  fillPositionResolution(h1_sigmaEta_vs_pt_endcap, h1_sigmaPhi_vs_pt_endcap, h2_deltaEta_vs_pt_endcap, h2_deltaPhi_vs_pt_endcap);
  fillPositionResolution(h1_sigmaEta_vs_pt_forward, h1_sigmaPhi_vs_pt_forward, h2_deltaEta_vs_pt_forward, h2_deltaPhi_vs_pt_forward);





 

  // **************        PF only  


  //--------------------------------- Rch ------------------------------------
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch0";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch0";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
  
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch0";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch0";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{ch}<20%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch0->GetYaxis()->SetRangeUser(0., 0.2);
  
  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch10_30";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch10_30 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch10_30->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch10_30->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch10_30->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch50_70";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch50_70 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{ch}<50%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch50_70->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch50_70->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch50_70->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);
 
  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rch90_100";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rch90_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{ch}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rch90_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rch90_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rch90_100->GetYaxis()->SetRangeUser(0., 0.2);


  //-------------------- Rgamma --------------------------------------------
  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaEta_vs_pt_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaEta_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaEta_vs_pt_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #eta Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 0<R_{gamma}<20%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma0";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma0 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 0<R_{gamma}<20%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma0->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma0->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma0->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma5_20";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma5_20 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 30<R_{gamma}<50%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma30_60";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma30_60 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 30<R_{gamma}<50%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_barrel_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaPhi_vs_pt_barrel_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_barrel_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);

  histoName = "sigmaPhi_vs_pt_endcap_"+algoType+"_Rgamma60_100";
  TH1F* h1_sigmaPhi_vs_pt_endcap_Rgamma60_100 = new TH1F(histoName.c_str(), "Jet Position #phi Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0, 60<R_{gamma}<100%)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaPhi_vs_pt_endcap_Rgamma60_100->GetYaxis()->SetRangeUser(0., 0.2);



 
  if( algoType=="PFItCone5" ) {
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
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma0,      h1_sigmaPhi_vs_pt_barrel_Rgamma0,      h2_deltaEta_vs_pt_barrel_Rgamma0,      h2_deltaPhi_vs_pt_barrel_Rgamma0);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma0,      h1_sigmaPhi_vs_pt_endcap_Rgamma0,      h2_deltaEta_vs_pt_endcap_Rgamma0,      h2_deltaPhi_vs_pt_endcap_Rgamma0);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma5_20,  h1_sigmaPhi_vs_pt_barrel_Rgamma5_20,  h2_deltaEta_vs_pt_barrel_Rgamma5_20,  h2_deltaPhi_vs_pt_barrel_Rgamma5_20);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma5_20,  h1_sigmaPhi_vs_pt_endcap_Rgamma5_20,  h2_deltaEta_vs_pt_endcap_Rgamma5_20,  h2_deltaPhi_vs_pt_endcap_Rgamma5_20);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma30_60,  h1_sigmaPhi_vs_pt_barrel_Rgamma30_60,  h2_deltaEta_vs_pt_barrel_Rgamma30_60,  h2_deltaPhi_vs_pt_barrel_Rgamma30_60);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma30_60,  h1_sigmaPhi_vs_pt_endcap_Rgamma30_60,  h2_deltaEta_vs_pt_endcap_Rgamma30_60,  h2_deltaPhi_vs_pt_endcap_Rgamma30_60);
    fillPositionResolution(h1_sigmaEta_vs_pt_barrel_Rgamma60_100, h1_sigmaPhi_vs_pt_barrel_Rgamma60_100, h2_deltaEta_vs_pt_barrel_Rgamma60_100, h2_deltaPhi_vs_pt_barrel_Rgamma60_100);
    fillPositionResolution(h1_sigmaEta_vs_pt_endcap_Rgamma60_100, h1_sigmaPhi_vs_pt_endcap_Rgamma60_100, h2_deltaEta_vs_pt_endcap_Rgamma60_100, h2_deltaPhi_vs_pt_endcap_Rgamma60_100);
  }

 
  delete h2_response_vs_pt_barrel;
  h2_response_vs_pt_barrel = 0;
  delete h2_response_vs_pt_endcap;
  h2_response_vs_pt_endcap = 0;
  delete h2_response_vs_pt_forward;
  h2_response_vs_pt_forward = 0;

  delete h2_response_vs_pt_eta0_25;
  h2_response_vs_pt_eta0_25 = 0;
  delete h2_response_vs_pt_eta14_25;
  h2_response_vs_pt_eta14_25 = 0;

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

  delete h2_response_vs_pt_eta0_25_Rch0;
  h2_response_vs_pt_eta0_25_Rch0 = 0;
  delete h2_response_vs_pt_eta14_25_Rch0;
  h2_response_vs_pt_eta14_25_Rch0 = 0;
  delete h2_response_vs_pt_eta0_25_Rch10_30;
  h2_response_vs_pt_eta0_25_Rch10_30 = 0;
  delete h2_response_vs_pt_eta14_25_Rch10_30;
  h2_response_vs_pt_eta14_25_Rch10_30 = 0;
  delete h2_response_vs_pt_eta0_25_Rch50_70;
  h2_response_vs_pt_eta0_25_Rch50_70 = 0;
  delete h2_response_vs_pt_eta14_25_Rch50_70;
  h2_response_vs_pt_eta14_25_Rch50_70 = 0;
  delete h2_response_vs_pt_eta0_25_Rch90_100;
  h2_response_vs_pt_eta0_25_Rch90_100 = 0;
  delete h2_response_vs_pt_eta14_25_Rch90_100;
  h2_response_vs_pt_eta14_25_Rch90_100 = 0;

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


  h1_deltaR_barrelEndcap->Write();
  h1_deltaR_forward->Write();
  h1_deltaR_veryForward->Write();

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

  h1_response_vs_pt_MEAN_eta0_25->Write();
  h1_response_vs_pt_MEAN_eta14_25->Write();

  h1_resolution_vs_pt_RMS_eta0_25->Write();
  h1_resolution_vs_pt_RMS_eta14_25->Write();

  h1_response_vs_pt_FIT_eta0_25->Write();
  h1_response_vs_pt_FIT_eta14_25->Write();

  h1_resolution_vs_pt_FIT_eta0_25->Write();
  h1_resolution_vs_pt_FIT_eta14_25->Write();

  h1_sigmaEta_vs_pt_barrel->Write();
  h1_sigmaEta_vs_pt_endcap->Write();
  h1_sigmaEta_vs_pt_forward->Write();

  h1_sigmaPhi_vs_pt_barrel->Write();
  h1_sigmaPhi_vs_pt_endcap->Write();
  h1_sigmaPhi_vs_pt_forward->Write();

  if( algoType == "PFItCone5" ) { 

    h1_response_vs_pt_onlyCh_MEAN_barrel->Write();
    h1_response_vs_pt_onlyCh_MEAN_endcap->Write();

    h1_response_vs_pt_onlyCh_FIT_barrel->Write();
    h1_response_vs_pt_onlyCh_FIT_endcap->Write();

    h1_resolution_vs_pt_onlyCh_RMS_barrel->Write();
    h1_resolution_vs_pt_onlyCh_RMS_endcap->Write();

    h1_resolution_vs_pt_onlyCh_FIT_barrel->Write();
    h1_resolution_vs_pt_onlyCh_FIT_endcap->Write();

    hp_Rgamma_vs_Rch_barrel->Write();
    hp_Rnh_vs_Rch_barrel->Write();

    hp_Rgamma_vs_Rch->Write();
    hp_Rnh_vs_Rch->Write();
    hp_Rmu_vs_Rch->Write();

    h1_Rch->Write();
    hp_Rch_vs_eta->Write();
    hp_Rch_vs_eta_50_80->Write();
    hp_Rch_vs_pt->Write();

    h1_Rgamma->Write();
    hp_Rgamma_vs_eta->Write();
    hp_Rgamma_vs_eta_50_80->Write();
    hp_Rgamma_vs_eta_Rch0->Write();
    hp_Rgamma_vs_pt->Write();

    h1_Rnh->Write();
    hp_Rnh_vs_eta->Write();
    hp_Rnh_vs_eta_50_80->Write();
    hp_Rnh_vs_eta_Rch0->Write();
    hp_Rnh_vs_pt->Write();

    h1_Rmu->Write();
    hp_Rmu_vs_eta->Write();
    hp_Rmu_vs_eta_50_80->Write();
    hp_Rmu_vs_eta_Rch0->Write();
    hp_Rmu_vs_pt->Write();

    h1_Rhf->Write();
    hp_Rhf_vs_eta->Write();
    hp_Rhf_vs_eta_50_80->Write();
    hp_Rhf_vs_eta_Rch0->Write();
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
    hp_RchGen_vs_pt_Rch0->Write();

    h1_RgammaGen->Write();
    hp_RgammaGen_vs_eta->Write();
    hp_RgammaGen_vs_eta_Rch0->Write();
    hp_RgammaGen_vs_pt->Write();
    hp_RgammaGen_vs_pt_Rch0->Write();

    h1_RnhGen->Write();
    hp_RnhGen_vs_eta->Write();
    hp_RnhGen_vs_eta_Rch0->Write();
    hp_RnhGen_vs_pt->Write();
    hp_RnhGen_vs_pt_Rch0->Write();

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

    h1_response_vs_Rch_FIT_barrel->Write();
    h1_response_vs_Rch_MEAN_barrel->Write();
    h1_response_vs_Rgamma_FIT_barrel->Write();
    h1_response_vs_Rgamma_MEAN_barrel->Write();
    h1_response_vs_Rnh_FIT_barrel->Write();
    h1_response_vs_Rnh_MEAN_barrel->Write();

    h1_resolution_vs_Rch_FIT_barrel->Write();
    h1_resolution_vs_Rch_RMS_barrel->Write();
    h1_resolution_vs_Rgamma_FIT_barrel->Write();
    h1_resolution_vs_Rgamma_RMS_barrel->Write();
    h1_resolution_vs_Rnh_FIT_barrel->Write();
    h1_resolution_vs_Rnh_RMS_barrel->Write();

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


    h1_response_vs_pt_FIT_eta0_25_Rch0->Write();
    h1_response_vs_pt_FIT_eta14_25_Rch0->Write();
    h1_response_vs_pt_FIT_eta0_25_Rch10_30->Write();
    h1_response_vs_pt_FIT_eta14_25_Rch10_30->Write();
    h1_response_vs_pt_FIT_eta0_25_Rch50_70->Write();
    h1_response_vs_pt_FIT_eta14_25_Rch50_70->Write();
    h1_response_vs_pt_FIT_eta0_25_Rch90_100->Write();
    h1_response_vs_pt_FIT_eta14_25_Rch90_100->Write();

    h1_resolution_vs_pt_FIT_eta0_25_Rch0->Write();
    h1_resolution_vs_pt_FIT_eta14_25_Rch0->Write();
    h1_resolution_vs_pt_FIT_eta0_25_Rch10_30->Write();
    h1_resolution_vs_pt_FIT_eta14_25_Rch10_30->Write();
    h1_resolution_vs_pt_FIT_eta0_25_Rch50_70->Write();
    h1_resolution_vs_pt_FIT_eta14_25_Rch50_70->Write();
    h1_resolution_vs_pt_FIT_eta0_25_Rch90_100->Write();
    h1_resolution_vs_pt_FIT_eta14_25_Rch90_100->Write();

    h1_response_vs_pt_MEAN_eta0_25_Rch0->Write();
    h1_response_vs_pt_MEAN_eta14_25_Rch0->Write();
    h1_response_vs_pt_MEAN_eta0_25_Rch10_30->Write();
    h1_response_vs_pt_MEAN_eta14_25_Rch10_30->Write();
    h1_response_vs_pt_MEAN_eta0_25_Rch50_70->Write();
    h1_response_vs_pt_MEAN_eta14_25_Rch50_70->Write();
    h1_response_vs_pt_MEAN_eta0_25_Rch90_100->Write();
    h1_response_vs_pt_MEAN_eta14_25_Rch90_100->Write();

    h1_resolution_vs_pt_RMS_eta0_25_Rch0->Write();
    h1_resolution_vs_pt_RMS_eta14_25_Rch0->Write();
    h1_resolution_vs_pt_RMS_eta0_25_Rch10_30->Write();
    h1_resolution_vs_pt_RMS_eta14_25_Rch10_30->Write();
    h1_resolution_vs_pt_RMS_eta0_25_Rch50_70->Write();
    h1_resolution_vs_pt_RMS_eta14_25_Rch50_70->Write();
    h1_resolution_vs_pt_RMS_eta0_25_Rch90_100->Write();
    h1_resolution_vs_pt_RMS_eta14_25_Rch90_100->Write();


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


void fillProfile(TH1F* h1_response_FIT, TH1F* h1_resolution_FIT, TH1F* h1_response_MEAN, TH1F* h1_resolution_RMS, TH2D* h2, std::string name) {

  std::string fileName = "Projections/"+name+".root";
  TFile* projectionFile;
  if( name!= "" ) {
    projectionFile = TFile::Open(fileName.c_str(), "RECREATE");
    projectionFile->cd();
  }

  for(int iBin=1; iBin<(h2->GetNbinsX()+1); ++iBin) {
   
    char histName[50];
    sprintf(histName, "projection_%d",iBin);
    TH1D* h1_projection = h2->ProjectionY(histName, iBin, iBin); 

    TF1* gaussian = new TF1("gaussian", "gaus");

    fitProjection(h1_projection, gaussian);

    if( name!="" ) {
        h1_projection->Write();
    }

    Float_t mu = gaussian->GetParameter(1);
    Float_t mu_err = gaussian->GetParError(1);
    h1_response_FIT->SetBinContent(iBin, mu);
    h1_response_FIT->SetBinError(iBin, mu_err);

    Float_t sigma = gaussian->GetParameter(2);
    Float_t resolution = (mu!=0.) ? sigma/mu : -1.;
    h1_resolution_FIT->SetBinContent(iBin, resolution);

    Float_t sigma_err = gaussian->GetParError(2);
    Float_t res_err = (mu!=0.) ? sqrt( sigma_err*sigma_err/(mu*mu) + mu_err*mu_err*sigma*sigma/(mu*mu*mu*mu) ) : 0.;
    h1_resolution_FIT->SetBinError(iBin, res_err);

    delete gaussian;
    gaussian = 0;

    Float_t n = h1_projection->GetEntries();
    Float_t mean = h1_projection->GetMean();
    Float_t mean_err = (n!=0) ? h1_projection->GetRMS()/sqrt(n) : 0.;
    h1_response_MEAN->SetBinContent(iBin, mean);
    h1_response_MEAN->SetBinError(iBin, mean_err);

    Float_t rms = h1_projection->GetRMS();
    Float_t rms_err = (n!=0) ? h1_projection->GetRMS()/sqrt(n) : 0.;
    resolution = (mean!=0.) ? rms/mean : -1.;
    res_err = (mean!=0.) ? sqrt( rms_err*rms_err/(mean*mean) + mean_err*mean_err*rms*rms/(mean*mean*mean*mean) ) : 0.;
    if( resolution != 0. ) {
      h1_resolution_RMS->SetBinContent(iBin, resolution);
      h1_resolution_RMS->SetBinError(iBin, res_err);
    }

    h1_projection = 0;

  } //for bins

  if(name!="") {
    projectionFile->Write();
    projectionFile->Close();
    delete projectionFile;
  }
  projectionFile = 0;

} //fill profile


void fillPositionResolution(TH1F* h1_sigmaEta, TH1F* h1_sigmaPhi, TH2D* h2_deltaEta, TH2D* h2_deltaPhi) {

  for(int iBin=1; iBin<(h2_deltaEta->GetNbinsX()+1); ++iBin) {

    TH1D* h1_projection = h2_deltaEta->ProjectionY("projectiony", iBin, iBin);

    TF1* gaussian = new TF1("gaussian", "gaus");
  
    fitProjection(h1_projection, gaussian);

    h1_sigmaEta->SetBinContent(iBin, gaussian->GetParameter(2));
    h1_sigmaEta->SetBinError(iBin, gaussian->GetParError(2));
  
    delete gaussian;
    gaussian = 0;

    h1_projection = 0;

  } //for bins eta

  for(int iBin=1; iBin<(h2_deltaPhi->GetNbinsX()+1); ++iBin) {

    TH1D* h1_projection = h2_deltaPhi->ProjectionY("projectiony", iBin, iBin);
  
    TF1* gaussian = new TF1("gaussian", "gaus");
  
    fitProjection(h1_projection, gaussian);

    h1_sigmaPhi->SetBinContent(iBin, gaussian->GetParameter(2));
    h1_sigmaPhi->SetBinError(iBin, gaussian->GetParError(2));
  
    delete gaussian;
    gaussian = 0;

    h1_projection = 0;

  } //for bins

} //fill position resolution


void fitProjection(TH1D* h1_projection, TF1* gaussian) {

   Float_t histMean = h1_projection->GetMean();
   Float_t histRMS = h1_projection->GetRMS();
 
   gaussian->SetParameter(0, h1_projection->GetMaximum());
   gaussian->SetParameter(1, histMean);
   gaussian->SetParameter(2, histRMS);
   if( histRMS==0. ) return;

   gaussian->SetParLimits(1, 0., 2*histMean);

   Float_t nSigma = 5.;

   Float_t lowerBound = histMean - nSigma*histRMS;
   Float_t upperBound = histMean + nSigma*histRMS;
   gaussian->SetRange( lowerBound, upperBound );
 
   h1_projection->Fit(gaussian, "R");

   int n_iter = 3;

   for(int i=0; i<n_iter; ++i) {

     Float_t lowerBound = gaussian->GetParameter(1) - nSigma*gaussian->GetParameter(2);
     Float_t upperBound = gaussian->GetParameter(1) + nSigma*gaussian->GetParameter(2);
     gaussian->SetRange( lowerBound, upperBound );
 
     h1_projection->Fit(gaussian, "R");
     
   }

}
