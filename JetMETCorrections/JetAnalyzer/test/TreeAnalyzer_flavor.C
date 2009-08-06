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


void fillFile(TFile* file, TTree* tree, std::string algoType, bool draw=false);


void flavor() {

  TChain* chain = new TChain("jetTree");

  int nFiles = 2;

  for(int i=1; i<(nFiles+1); ++i) {

    char treeName[70];
    sprintf(treeName, "2ndLevelOutputFiles/2ndLevelOutputFile_%d.root/jetTree", i);
    chain->Add(treeName);

  } //for


  bool draw=true;

  TFile* outfile = new TFile("flavorOutputFile.root", "RECREATE");
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
  tree->SetBranchAddress(eGen_name.c_str(), eGen);
  Float_t ptGen[40];
  std::string ptGen_name = "ptJetGen_"+algoType;
  tree->SetBranchAddress(ptGen_name.c_str(), ptGen);
  Float_t etaGen[40];
  std::string etaGen_name = "etaJetGen_"+algoType;
  tree->SetBranchAddress(etaGen_name.c_str(), etaGen);
  Float_t phiGen[40];
  std::string phiGen_name = "phiJetGen_"+algoType;
  tree->SetBranchAddress(phiGen_name.c_str(), phiGen);
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
 
  Int_t pdgIdPart[40];
  std::string pdgIdPart_name = "pdgIdPart_"+algoType;
  tree->SetBranchAddress(pdgIdPart_name.c_str(), pdgIdPart);
  Float_t ptPart[40];
  std::string ptPart_name = "ptPart_"+algoType;
  tree->SetBranchAddress(ptPart_name.c_str(), ptPart);
  Float_t etaPart[40];
  std::string etaPart_name = "etaPart_"+algoType;
  tree->SetBranchAddress(etaPart_name.c_str(), etaPart);
  Float_t phiPart[40];
  std::string phiPart_name = "phiPart_"+algoType;
  tree->SetBranchAddress(phiPart_name.c_str(), phiPart);



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


  std::string histoName;

  histoName = "response_vs_pt_light_"+algoType;
  TH2D* h2_response_vs_pt_light= new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T}", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_light_"+algoType;
  TH2D* h2_response_vs_pt_barrel_light = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_barrel_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_light_"+algoType;
  TH2D* h2_response_vs_pt_endcap_light = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_endcap_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_forward_light_"+algoType;
  TH2D* h2_response_vs_pt_forward_light = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_forward_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_eta_light_"+algoType;
  TH2D* h2_response_vs_eta_light = new TH2D(histoName.c_str(), "Jet Response vs. True #eta", 100, -5.5, 5.5, 100, 0., 1.5);
  h2_response_vs_eta_light->SetXTitle("#eta");
  h2_response_vs_eta_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");



  histoName = "response_vs_pt_gluon_"+algoType;
  TH2D* h2_response_vs_pt_gluon= new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T}", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_barrel_gluon_"+algoType;
  TH2D* h2_response_vs_pt_barrel_gluon = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_barrel_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_endcap_gluon_"+algoType;
  TH2D* h2_response_vs_pt_endcap_gluon = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_endcap_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_pt_forward_gluon_"+algoType;
  TH2D* h2_response_vs_pt_forward_gluon = new TH2D(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_forward_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  histoName = "response_vs_eta_gluon_"+algoType;
  TH2D* h2_response_vs_eta_gluon = new TH2D(histoName.c_str(), "Jet Response vs. True #eta", 100, -5.5, 5.5, 100, 0., 1.5);
  h2_response_vs_eta_gluon->SetXTitle("#eta");
  h2_response_vs_eta_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");


  histoName = "RchGen_vs_eta_Rch0_gluon_"+algoType;
  TProfile* hp_RchGen_vs_eta_Rch0_gluon = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RchGen_vs_eta_Rch0_gluon->SetXTitle("#eta^{GEN}");
  hp_RchGen_vs_eta_Rch0_gluon->SetYTitle("entries");

  histoName = "RgammaGen_vs_eta_Rch0_gluon_"+algoType;
  TProfile* hp_RgammaGen_vs_eta_Rch0_gluon = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RgammaGen_vs_eta_Rch0_gluon->SetXTitle("#eta^{GEN}");
  hp_RgammaGen_vs_eta_Rch0_gluon->SetYTitle("entries");

  histoName = "RnhGen_vs_eta_Rch0_gluon_"+algoType;
  TProfile* hp_RnhGen_vs_eta_Rch0_gluon = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RnhGen_vs_eta_Rch0_gluon->SetXTitle("#eta^{GEN}");
  hp_RnhGen_vs_eta_Rch0_gluon->SetYTitle("entries");

  histoName = "RchGen_vs_eta_Rch0_light_"+algoType;
  TProfile* hp_RchGen_vs_eta_Rch0_light = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RchGen_vs_eta_Rch0_light->SetXTitle("#eta^{GEN}");
  hp_RchGen_vs_eta_Rch0_light->SetYTitle("entries");

  histoName = "RgammaGen_vs_eta_Rch0_light_"+algoType;
  TProfile* hp_RgammaGen_vs_eta_Rch0_light = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RgammaGen_vs_eta_Rch0_light->SetXTitle("#eta^{GEN}");
  hp_RgammaGen_vs_eta_Rch0_light->SetYTitle("entries");

  histoName = "RnhGen_vs_eta_Rch0_light_"+algoType;
  TProfile* hp_RnhGen_vs_eta_Rch0_light = new TProfile(histoName.c_str(), "R_{ch} vs. eta", 100, -5.5, 5.5, 0., 1000.);
  hp_RnhGen_vs_eta_Rch0_light->SetXTitle("#eta^{GEN}");
  hp_RnhGen_vs_eta_Rch0_light->SetYTitle("entries");


  int nEntries = tree->GetEntries();

  //for(int iEntry=0; iEntry<1000000; ++iEntry) {
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

      Float_t Rch = eTracksReco[iRecoJet]/iEReco;
      Float_t Rgamma = ePhotonsReco[iRecoJet]/iEReco;
      Float_t Rnh = eNeutralHadronsReco[iRecoJet]/iEReco;
      Float_t Rmu = eMuonsReco[iRecoJet]/iEReco;
      Float_t Rhf = eHFhadronsReco[iRecoJet]/iEReco;

      Float_t RchGen = eTracksGen[iRecoJet]/iEGen;
      Float_t RgammaGen = ePhotonsGen[iRecoJet]/iEGen;
      Float_t RnhGen = eNeutralHadronsGen[iRecoJet]/iEGen;

      Float_t deltaEta = iEtaReco-iEtaGen;
      Float_t deltaPhi = iPhiReco-iPhiGen;
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

      if( matched ) {
   
        if( pdgIdPart[iRecoJet] == 21) { //gluon

          if( algoType=="PFItCone5" ) {
            if( Rch==0. ) {
              hp_RchGen_vs_eta_Rch0_gluon->Fill(iEtaGen, RchGen);
              hp_RgammaGen_vs_eta_Rch0_gluon->Fill(iEtaGen, RgammaGen);
              hp_RnhGen_vs_eta_Rch0_gluon->Fill(iEtaGen, RnhGen); 
            }
          }

          if( iPtGen>15. )
            h2_response_vs_eta_gluon->Fill(iEtaGen, response);

          h2_response_vs_pt_gluon->Fill(iPtGen, response);

          if( fabs(iEtaGen)<1.4 ) 
            h2_response_vs_pt_barrel_gluon->Fill(iPtGen, response);
          else if( fabs(iEtaGen)<3. ) 
            h2_response_vs_pt_endcap_gluon->Fill(iPtGen, response);
          else if( fabs(iEtaGen)<5. ) 
            h2_response_vs_pt_forward_gluon->Fill(iPtGen, response);

        } //if gluon

        if( (fabs(pdgIdPart[iRecoJet]) == 1) ||
            (fabs(pdgIdPart[iRecoJet]) == 2) ||
            (fabs(pdgIdPart[iRecoJet]) == 3) ) { //light quark

          if( algoType=="PFItCone5" ) {
            if( Rch==0. ) {
              hp_RchGen_vs_eta_Rch0_light->Fill(iEtaGen, RchGen);
              hp_RgammaGen_vs_eta_Rch0_light->Fill(iEtaGen, RgammaGen);
              hp_RnhGen_vs_eta_Rch0_light->Fill(iEtaGen, RnhGen); 
            }
          }


          if( iPtGen>15. )
            h2_response_vs_eta_light->Fill(iEtaGen, response);

          h2_response_vs_pt_light->Fill(iPtGen, response);

          if( fabs(iEtaGen)<1.4 ) 
            h2_response_vs_pt_barrel_light->Fill(iPtGen, response);
          else if( fabs(iEtaGen)<3. ) 
            h2_response_vs_pt_endcap_light->Fill(iPtGen, response);
          else if( fabs(iEtaGen)<5. ) 
            h2_response_vs_pt_forward_light->Fill(iPtGen, response);

        } //if light quark



      } //if matched

    } //for reco jets

  } //for entries

 
  // ----- gluon plots:

  histoName = "response_vs_eta_FIT_gluon_"+algoType;
  TH1F* h1_response_vs_eta_FIT_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5.5, 5.5);
  h1_response_vs_eta_FIT_gluon->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_FIT_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_eta_FIT_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_gluon_"+algoType;
  TH1F* h1_response_vs_pt_FIT_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_FIT_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_barrel_gluon_"+algoType;
  TH1F* h1_response_vs_pt_FIT_barrel_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_endcap_gluon_"+algoType;
  TH1F* h1_response_vs_pt_FIT_endcap_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_forward_gluon_"+algoType;
  TH1F* h1_response_vs_pt_FIT_forward_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_FIT_forward_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_forward_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_forward_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_eta_FIT_gluon_"+algoType;
  TH1F* h1_resolution_vs_eta_FIT_gluon = new TH1F(histoName.c_str(), "Jet Resolution vs. True Jet #eta", 100, -5.5, 5.5);
  h1_resolution_vs_eta_FIT_gluon->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_FIT_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_FIT_gluon->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_gluon->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_barrel_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_barrel_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel_gluon->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_endcap_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_endcap_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_endcap_gluon->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_forward_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_forward_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_forward_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_forward_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_forward_gluon->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_eta_MEAN_gluon_"+algoType;
  TH1F* h1_response_vs_eta_MEAN_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5.5, 5.5);
  h1_response_vs_eta_MEAN_gluon->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_MEAN_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_eta_MEAN_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_gluon_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_MEAN_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_barrel_gluon_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_barrel_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_endcap_gluon_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_endcap_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_forward_gluon_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_forward_gluon = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_forward_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_forward_gluon->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_forward_gluon->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_eta_RMS_gluon_"+algoType;
  TH1F* h1_resolution_vs_eta_RMS_gluon = new TH1F(histoName.c_str(), "Jet Resolution vs. True Jet #eta", 100, -5.5, 5.5);
  h1_resolution_vs_eta_RMS_gluon->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_RMS_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_RMS_gluon->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_gluon->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_barrel_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_gluon->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_endcap_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_endcap_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_gluon->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_forward_gluon_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_forward_gluon = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_forward_gluon->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_forward_gluon->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_forward_gluon->GetYaxis()->SetRangeUser(0., 0.5);



  // ----- light quark plots:

  histoName = "response_vs_eta_FIT_light_"+algoType;
  TH1F* h1_response_vs_eta_FIT_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5.5, 5.5);
  h1_response_vs_eta_FIT_light->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_FIT_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_eta_FIT_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_light_"+algoType;
  TH1F* h1_response_vs_pt_FIT_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_FIT_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_barrel_light_"+algoType;
  TH1F* h1_response_vs_pt_FIT_barrel_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_endcap_light_"+algoType;
  TH1F* h1_response_vs_pt_FIT_endcap_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_FIT_forward_light_"+algoType;
  TH1F* h1_response_vs_pt_FIT_forward_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_FIT_forward_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_forward_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_forward_light->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_eta_FIT_light_"+algoType;
  TH1F* h1_resolution_vs_eta_FIT_light = new TH1F(histoName.c_str(), "Jet Resolution vs. True Jet #eta", 100, -5.5, 5.5);
  h1_resolution_vs_eta_FIT_light->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_FIT_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_FIT_light->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_FIT_light_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_light->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_barrel_light_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_barrel_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel_light->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_endcap_light_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_endcap_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_endcap_light->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_FIT_forward_light_"+algoType;
  TH1F* h1_resolution_vs_pt_FIT_forward_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_forward_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_forward_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_forward_light->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "response_vs_eta_MEAN_light_"+algoType;
  TH1F* h1_response_vs_eta_MEAN_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet #eta", 100, -5.5, 5.5);
  h1_response_vs_eta_MEAN_light->SetXTitle("#eta^{GEN}");
  h1_response_vs_eta_MEAN_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_eta_MEAN_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_light_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_MEAN_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_barrel_light_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_barrel_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_endcap_light_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_endcap_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_light->GetYaxis()->SetRangeUser(0.2, 1.);
  
  histoName = "response_vs_pt_MEAN_forward_light_"+algoType;
  TH1F* h1_response_vs_pt_MEAN_forward_light = new TH1F(histoName.c_str(), "Jet Response vs. True Jet p_{T} (3<|#eta_{GEN}|<5)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_forward_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_forward_light->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_forward_light->GetYaxis()->SetRangeUser(0.2, 1.);
 
  histoName = "resolution_vs_eta_RMS_light_"+algoType;
  TH1F* h1_resolution_vs_eta_RMS_light = new TH1F(histoName.c_str(), "Jet Resolution vs. True Jet #eta", 100, -5.5, 5.5);
  h1_resolution_vs_eta_RMS_light->SetXTitle("#eta^{GEN}");
  h1_resolution_vs_eta_RMS_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_eta_RMS_light->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_light_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_light->GetYaxis()->SetRangeUser(0., 0.5);

  histoName = "resolution_vs_pt_RMS_barrel_light_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_barrel_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_light->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_endcap_light_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_endcap_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_light->GetYaxis()->SetRangeUser(0., 0.5);
 
  histoName = "resolution_vs_pt_RMS_forward_light_"+algoType;
  TH1F* h1_resolution_vs_pt_RMS_forward_light = new TH1F(histoName.c_str(), "Jet Energy Resolution vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_forward_light->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_forward_light->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_forward_light->GetYaxis()->SetRangeUser(0., 0.5);


  std::string name = "";

  if(draw) {

    fillProfile(h1_response_vs_eta_FIT_gluon, h1_resolution_vs_eta_FIT_gluon, h1_response_vs_eta_MEAN_gluon, h1_resolution_vs_eta_RMS_gluon, h2_response_vs_eta_gluon);
    fillProfile(h1_response_vs_eta_FIT_light, h1_resolution_vs_eta_FIT_light, h1_response_vs_eta_MEAN_light, h1_resolution_vs_eta_RMS_light, h2_response_vs_eta_light);

    char name[40];
    sprintf(name, "%s_gluon", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_gluon, h1_resolution_vs_pt_FIT_gluon, h1_response_vs_pt_MEAN_gluon, h1_resolution_vs_pt_RMS_gluon, h2_response_vs_pt_gluon, name);
    sprintf(name, "%s_barrel_gluon", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_barrel_gluon, h1_resolution_vs_pt_FIT_barrel_gluon, h1_response_vs_pt_MEAN_barrel_gluon, h1_resolution_vs_pt_RMS_barrel_gluon, h2_response_vs_pt_barrel_gluon, name);
    sprintf(name, "%s_endcap_gluon", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_endcap_gluon, h1_resolution_vs_pt_FIT_endcap_gluon, h1_response_vs_pt_MEAN_endcap_gluon, h1_resolution_vs_pt_RMS_endcap_gluon, h2_response_vs_pt_endcap_gluon, name);
    sprintf(name, "%s_forward_gluon", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_forward_gluon, h1_resolution_vs_pt_FIT_forward_gluon, h1_response_vs_pt_MEAN_forward_gluon, h1_resolution_vs_pt_RMS_forward_gluon, h2_response_vs_pt_forward_gluon, name);

    sprintf(name, "%s_light", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_light, h1_resolution_vs_pt_FIT_light, h1_response_vs_pt_MEAN_light, h1_resolution_vs_pt_RMS_light, h2_response_vs_pt_light, name);
    sprintf(name, "%s_barrel_light", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_barrel_light, h1_resolution_vs_pt_FIT_barrel_light, h1_response_vs_pt_MEAN_barrel_light, h1_resolution_vs_pt_RMS_barrel_light, h2_response_vs_pt_barrel_light, name);
    sprintf(name, "%s_endcap_light", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_endcap_light, h1_resolution_vs_pt_FIT_endcap_light, h1_response_vs_pt_MEAN_endcap_light, h1_resolution_vs_pt_RMS_endcap_light, h2_response_vs_pt_endcap_light, name);
    sprintf(name, "%s_forward_light", algoType.c_str());
    fillProfile(h1_response_vs_pt_FIT_forward_light, h1_resolution_vs_pt_FIT_forward_light, h1_response_vs_pt_MEAN_forward_light, h1_resolution_vs_pt_RMS_forward_light, h2_response_vs_pt_forward_light, name);

  } else {

    fillProfile(h1_response_vs_eta_FIT_gluon, h1_resolution_vs_eta_FIT_gluon, h1_response_vs_eta_MEAN_gluon, h1_resolution_vs_eta_RMS_gluon, h2_response_vs_eta_gluon);
    fillProfile(h1_response_vs_pt_FIT_gluon, h1_resolution_vs_pt_FIT_gluon, h1_response_vs_pt_MEAN_gluon, h1_resolution_vs_pt_RMS_gluon, h2_response_vs_pt_gluon);
    fillProfile(h1_response_vs_pt_FIT_barrel_gluon, h1_resolution_vs_pt_FIT_barrel_gluon, h1_response_vs_pt_MEAN_barrel_gluon, h1_resolution_vs_pt_RMS_barrel_gluon, h2_response_vs_pt_barrel_gluon);
    fillProfile(h1_response_vs_pt_FIT_endcap_gluon, h1_resolution_vs_pt_FIT_endcap_gluon, h1_response_vs_pt_MEAN_endcap_gluon, h1_resolution_vs_pt_RMS_endcap_gluon, h2_response_vs_pt_endcap_gluon);
    fillProfile(h1_response_vs_pt_FIT_forward_gluon, h1_resolution_vs_pt_FIT_forward_gluon, h1_response_vs_pt_MEAN_forward_gluon, h1_resolution_vs_pt_RMS_forward_gluon, h2_response_vs_pt_forward_gluon);

    fillProfile(h1_response_vs_eta_FIT_light, h1_resolution_vs_eta_FIT_light, h1_response_vs_eta_MEAN_light, h1_resolution_vs_eta_RMS_light, h2_response_vs_eta_light);
    fillProfile(h1_response_vs_pt_FIT_light, h1_resolution_vs_pt_FIT_light, h1_response_vs_pt_MEAN_light, h1_resolution_vs_pt_RMS_light, h2_response_vs_pt_light);
    fillProfile(h1_response_vs_pt_FIT_barrel_light, h1_resolution_vs_pt_FIT_barrel_light, h1_response_vs_pt_MEAN_barrel_light, h1_resolution_vs_pt_RMS_barrel_light, h2_response_vs_pt_barrel_light);
    fillProfile(h1_response_vs_pt_FIT_endcap_light, h1_resolution_vs_pt_FIT_endcap_light, h1_response_vs_pt_MEAN_endcap_light, h1_resolution_vs_pt_RMS_endcap_light, h2_response_vs_pt_endcap_light);
    fillProfile(h1_response_vs_pt_FIT_forward_light, h1_resolution_vs_pt_FIT_forward_light, h1_response_vs_pt_MEAN_forward_light, h1_resolution_vs_pt_RMS_forward_light, h2_response_vs_pt_forward_light);

  }

  delete h2_response_vs_pt_gluon;
  h2_response_vs_pt_gluon = 0;
  delete h2_response_vs_pt_barrel_gluon;
  h2_response_vs_pt_barrel_gluon = 0;
  delete h2_response_vs_pt_endcap_gluon;
  h2_response_vs_pt_endcap_gluon = 0;
  delete h2_response_vs_pt_forward_gluon;
  h2_response_vs_pt_forward_gluon = 0;

  delete h2_response_vs_pt_light;
  h2_response_vs_pt_light = 0;
  delete h2_response_vs_pt_barrel_light;
  h2_response_vs_pt_barrel_light = 0;
  delete h2_response_vs_pt_endcap_light;
  h2_response_vs_pt_endcap_light = 0;
  delete h2_response_vs_pt_forward_light;
  h2_response_vs_pt_forward_light = 0;


  file->cd();


  h1_response_vs_eta_MEAN_gluon->Write();
  h1_response_vs_pt_MEAN_gluon->Write();
  h1_response_vs_pt_MEAN_barrel_gluon->Write();
  h1_response_vs_pt_MEAN_endcap_gluon->Write();
  h1_response_vs_pt_MEAN_forward_gluon->Write();

  h1_resolution_vs_eta_RMS_gluon->Write();
  h1_resolution_vs_pt_RMS_gluon->Write();
  h1_resolution_vs_pt_RMS_barrel_gluon->Write();
  h1_resolution_vs_pt_RMS_endcap_gluon->Write();
  h1_resolution_vs_pt_RMS_forward_gluon->Write();

  h1_response_vs_eta_FIT_gluon->Write();
  h1_response_vs_pt_FIT_gluon->Write();
  h1_response_vs_pt_FIT_barrel_gluon->Write();
  h1_response_vs_pt_FIT_endcap_gluon->Write();
  h1_response_vs_pt_FIT_forward_gluon->Write();

  h1_resolution_vs_eta_FIT_gluon->Write();
  h1_resolution_vs_pt_FIT_gluon->Write();
  h1_resolution_vs_pt_FIT_barrel_gluon->Write();
  h1_resolution_vs_pt_FIT_endcap_gluon->Write();
  h1_resolution_vs_pt_FIT_forward_gluon->Write();


  h1_response_vs_eta_MEAN_light->Write();
  h1_response_vs_pt_MEAN_light->Write();
  h1_response_vs_pt_MEAN_barrel_light->Write();
  h1_response_vs_pt_MEAN_endcap_light->Write();
  h1_response_vs_pt_MEAN_forward_light->Write();

  h1_resolution_vs_eta_RMS_light->Write();
  h1_resolution_vs_pt_RMS_light->Write();
  h1_resolution_vs_pt_RMS_barrel_light->Write();
  h1_resolution_vs_pt_RMS_endcap_light->Write();
  h1_resolution_vs_pt_RMS_forward_light->Write();

  h1_response_vs_eta_FIT_light->Write();
  h1_response_vs_pt_FIT_light->Write();
  h1_response_vs_pt_FIT_barrel_light->Write();
  h1_response_vs_pt_FIT_endcap_light->Write();
  h1_response_vs_pt_FIT_forward_light->Write();

  h1_resolution_vs_eta_FIT_light->Write();
  h1_resolution_vs_pt_FIT_light->Write();
  h1_resolution_vs_pt_FIT_barrel_light->Write();
  h1_resolution_vs_pt_FIT_endcap_light->Write();
  h1_resolution_vs_pt_FIT_forward_light->Write();

  if( algoType=="PFItCone5" ) {

    hp_RchGen_vs_eta_Rch0_gluon->Write();
    hp_RgammaGen_vs_eta_Rch0_gluon->Write();
    hp_RnhGen_vs_eta_Rch0_gluon->Write();

    hp_RchGen_vs_eta_Rch0_light->Write();
    hp_RgammaGen_vs_eta_Rch0_light->Write();
    hp_RnhGen_vs_eta_Rch0_light->Write();

  }


} //fillFile

