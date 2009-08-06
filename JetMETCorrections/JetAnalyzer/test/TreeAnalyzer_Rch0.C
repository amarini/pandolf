#include <TH2D.h>
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



void Rch0() {

  TChain* chain = new TChain("jetTree");

  for(int i=1; i<50; ++i) {

    char treeName[70];
    sprintf(treeName, "2ndLevelOutputFiles/2ndLevelOutputFile_%d.root/jetTree", i);
    chain->Add(treeName);

  } //for


  TFile* outFile = new TFile("Rch0.root", "RECREATE");

  Int_t nJetReco_calo;
  chain->SetBranchAddress("nJetReco_caloItCone5", &nJetReco_calo);
  Float_t etaReco_calo[40];
  chain->SetBranchAddress("etaJetReco_caloItCone5", etaReco_calo);
  Float_t ptReco_calo[40];
  chain->SetBranchAddress("ptJetReco_caloItCone5", ptReco_calo);
  Float_t eReco_calo[40];
  chain->SetBranchAddress("eJetReco_caloItCone5", eReco_calo);
  Float_t phiReco_calo[40];
  chain->SetBranchAddress("phiJetReco_caloItCone5", phiReco_calo);

  // gen jets matched to calo it Cone 0.5:
  Float_t etaGen_calo[40];
  chain->SetBranchAddress("etaJetGen_caloItCone5", etaGen_calo);
  Float_t ptGen_calo[40];
  chain->SetBranchAddress("ptJetGen_caloItCone5", ptGen_calo);
  Float_t eGen_calo[40];
  chain->SetBranchAddress("eJetGen_caloItCone5", eGen_calo);
  Float_t phiGen_calo[40];
  chain->SetBranchAddress("phiJetGen_caloItCone5", phiGen_calo);


  Int_t nJetReco_PF;
  chain->SetBranchAddress("nJetReco_PFItCone5", &nJetReco_PF);
  Float_t etaReco_PF[40];
  chain->SetBranchAddress("etaJetReco_PFItCone5", etaReco_PF);
  Float_t ptReco_PF[40];
  chain->SetBranchAddress("ptJetReco_PFItCone5", ptReco_PF);
  Float_t eReco_PF[40];
  chain->SetBranchAddress("eJetReco_PFItCone5", eReco_PF);
  Float_t phiReco_PF[40];
  chain->SetBranchAddress("phiJetReco_PFItCone5", phiReco_PF);
  Float_t eTracksReco_PF[40];
  chain->SetBranchAddress("eTracksReco_PFItCone5", eTracksReco_PF);

  // gen jets matched to PF it Cone 0.5:
  Float_t etaGen_PF[40];
  chain->SetBranchAddress("etaJetGen_PFItCone5", etaGen_PF);
  Float_t ptGen_PF[40];
  chain->SetBranchAddress("ptJetGen_PFItCone5", ptGen_PF);
  Float_t eGen_PF[40];
  chain->SetBranchAddress("eJetGen_PFItCone5", eGen_PF);
  Float_t phiGen_PF[40];
  chain->SetBranchAddress("phiJetGen_PFItCone5", phiGen_PF);


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


  //----------------- CALO histograms

  TH2D* h2_deltaEta_vs_pt_calo = new TH2D("deltaEta_vs_pt_calo", "#eta^{RECO}-#eta^{GEN}", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_calo->SetYTitle("#Delta#eta");

  TH2D* h2_deltaEta_vs_pt_barrel_calo = new TH2D("deltaEta_vs_pt_barrel_calo", "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_calo->SetYTitle("#Delta#eta");

  TH2D* h2_deltaEta_vs_pt_endcap_calo = new TH2D("deltaEta_vs_pt_endcap_calo", "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_calo->SetYTitle("#Delta#eta");

  TH2D* h2_deltaEta_vs_pt_forward_calo = new TH2D("deltaEta_vs_pt_forward_calo", "#eta^{RECO}-#eta^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_forward_calo->SetYTitle("#Delta#eta");

  TH2D* h2_deltaPhi_vs_pt_calo = new TH2D("deltaPhi_vs_pt_calo", "#phi^{RECO}-#phi^{GEN}", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_calo->SetYTitle("#Delta#phi");

  TH2D* h2_deltaPhi_vs_pt_barrel_calo = new TH2D("deltaPhi_vs_pt_barrel_calo", "#phi^{RECO}-#phi^{GEN} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_calo->SetYTitle("#Delta#phi");

  TH2D* h2_deltaPhi_vs_pt_endcap_calo = new TH2D("deltaPhi_vs_pt_endcap_calo", "#phi^{RECO}-#phi^{GEN} (1.4<|#phi_{GEN}|<3.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_calo->SetYTitle("#Delta#phi");

  TH2D* h2_deltaPhi_vs_pt_forward_calo = new TH2D("deltaPhi_vs_pt_forward_calo", "#phi^{RECO}-#phi^{GEN} (3.0<|#phi_{GEN}|<5.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_forward_calo->SetYTitle("#Delta#phi");

  TH2D* h2_response_vs_pt_calo = new TH2D("response_vs_pt_calo", "Jet Response vs. True Jet p_{T}", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_response_vs_pt_barrel_calo = new TH2D("response_vs_pt_barrel_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_response_vs_pt_endcap_calo = new TH2D("response_vs_pt_endcap_calo", "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_response_vs_pt_forward_calo = new TH2D("response_vs_pt_forward_calo", "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_responseE_vs_E_calo = new TH2D("responseE_vs_E_calo", "Jet Response vs. True Jet p_{T}", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_calo->SetYTitle("E^{RECO}/E^{GEN}");

  TH2D* h2_responseE_vs_E_barrel_calo = new TH2D("responseE_vs_E_barrel_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_barrel_calo->SetYTitle("E^{RECO}/E^{GEN}");

  TH2D* h2_responseE_vs_E_endcap_calo = new TH2D("responseE_vs_E_endcap_calo", "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_endcap_calo->SetYTitle("E^{RECO}/E^{GEN}");

  TH2D* h2_responseE_vs_E_forward_calo = new TH2D("responseE_vs_E_forward_calo", "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_forward_calo->SetYTitle("E^{RECO}/E^{GEN}");


  //----------------- PF histograms

  TH2D* h2_deltaEta_vs_pt_PF = new TH2D("deltaEta_vs_pt_PF", "#eta^{RECO}-#eta^{GEN}", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_PF->SetYTitle("#Delta#eta");

  TH2D* h2_deltaEta_vs_pt_barrel_PF = new TH2D("deltaEta_vs_pt_barrel_PF", "#eta^{RECO}-#eta^{GEN} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_barrel_PF->SetYTitle("#Delta#eta");

  TH2D* h2_deltaEta_vs_pt_endcap_PF = new TH2D("deltaEta_vs_pt_endcap_PF", "#eta^{RECO}-#eta^{GEN} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_endcap_PF->SetYTitle("#Delta#eta");

  TH2D* h2_deltaEta_vs_pt_forward_PF = new TH2D("deltaEta_vs_pt_forward_PF", "#eta^{RECO}-#eta^{GEN} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaEta_vs_pt_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEta_vs_pt_forward_PF->SetYTitle("#Delta#eta");

  TH2D* h2_deltaPhi_vs_pt_PF = new TH2D("deltaPhi_vs_pt_PF", "#phi^{RECO}-#phi^{GEN}", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_PF->SetYTitle("#Delta#phi");

  TH2D* h2_deltaPhi_vs_pt_barrel_PF = new TH2D("deltaPhi_vs_pt_barrel_PF", "#phi^{RECO}-#phi^{GEN} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_barrel_PF->SetYTitle("#Delta#phi");

  TH2D* h2_deltaPhi_vs_pt_endcap_PF = new TH2D("deltaPhi_vs_pt_endcap_PF", "#phi^{RECO}-#phi^{GEN} (1.4<|#phi_{GEN}|<3.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_endcap_PF->SetYTitle("#Delta#phi");

  TH2D* h2_deltaPhi_vs_pt_forward_PF = new TH2D("deltaPhi_vs_pt_forward_PF", "#phi^{RECO}-#phi^{GEN} (3.0<|#phi_{GEN}|<5.0)", nBins-2, Lower, 100, -0.8, 0.8);
  h2_deltaPhi_vs_pt_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaPhi_vs_pt_forward_PF->SetYTitle("#Delta#phi");

  TH2D* h2_response_vs_pt_PF = new TH2D("response_vs_pt_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_response_vs_pt_barrel_PF = new TH2D("response_vs_pt_barrel_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_barrel_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_response_vs_pt_endcap_PF = new TH2D("response_vs_pt_endcap_PF", "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_endcap_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_response_vs_pt_forward_PF = new TH2D("response_vs_pt_forward_PF", "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_response_vs_pt_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_response_vs_pt_forward_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  TH2D* h2_responseE_vs_E_PF = new TH2D("responseE_vs_E_PF", "Jet Response vs. True Jet p_{T}", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_PF->SetYTitle("E^{RECO}/E^{GEN}");

  TH2D* h2_responseE_vs_E_barrel_PF = new TH2D("responseE_vs_E_barrel_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_barrel_PF->SetYTitle("E^{RECO}/E^{GEN}");

  TH2D* h2_responseE_vs_E_endcap_PF = new TH2D("responseE_vs_E_endcap_PF", "Jet Response vs. True Jet p_{T} (1.4<|#eta_{GEN}|<3.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_endcap_PF->SetYTitle("E^{RECO}/E^{GEN}");

  TH2D* h2_responseE_vs_E_forward_PF = new TH2D("responseE_vs_E_forward_PF", "Jet Response vs. True Jet p_{T} (3.0<|#eta_{GEN}|<5.0)", nBins-2, Lower, 100, 0., 1.5);
  h2_responseE_vs_E_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_responseE_vs_E_forward_PF->SetYTitle("E^{RECO}/E^{GEN}");




  TH1F* h1_deltaR_PF2calo = new TH1F("deltaR_PF2calo", "", 100, 0., 1.);



  int nEntries = chain->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    chain->GetEntry(iEntry);

    for(int iRecoJet_PF = 0; iRecoJet_PF<nJetReco_PF; ++iRecoJet_PF) {

      if( iRecoJet_PF >= 40 ) continue;

      Float_t iEtaGen_PF = etaGen_PF[iRecoJet_PF];
      Float_t iPhiGen_PF = phiGen_PF[iRecoJet_PF];
      Float_t iPtGen_PF = ptGen_PF[iRecoJet_PF];
      Float_t iEGen_PF = eGen_PF[iRecoJet_PF];

      Float_t iEtaReco_PF = etaReco_PF[iRecoJet_PF];
      Float_t iPhiReco_PF = phiReco_PF[iRecoJet_PF];
      Float_t iThetaReco_PF = 2*atan(exp(-iEtaReco_PF));
      Float_t iPtReco_PF = ptReco_PF[iRecoJet_PF];
      Float_t iEReco_PF = eReco_PF[iRecoJet_PF];
      Float_t iEtReco_PF = iEReco_PF*sin(iThetaReco_PF);

      Float_t Rch = eTracksReco_PF[iRecoJet_PF]/iEReco_PF;

      if( Rch>0. ) continue;

      Float_t deltaEta_PF = iEtaReco_PF-iEtaGen_PF;
      Float_t deltaPhi_PF = iPhiReco_PF-iPhiGen_PF;
      Float_t pi = 3.14159;
      if(deltaPhi_PF < -pi) deltaPhi_PF += 2.*pi;
      if(deltaPhi_PF > pi) deltaPhi_PF -= 2.*pi;

      Float_t deltaR_PF2gen = sqrt( deltaEta_PF*deltaEta_PF + deltaPhi_PF*deltaPhi_PF );

      if( deltaR_PF2gen <0.1 ) { //matched
    
        Float_t response_PF = iPtReco_PF/iPtGen_PF;
        Float_t responseE_PF = iEReco_PF/iEGen_PF;


        // PF jet matched and Rch=0. lets look for the corresponding caloJet:
        for(unsigned iRecoJet_calo=0; iRecoJet_calo<nJetReco_calo; ++iRecoJet_calo) {

          Float_t iEtaGen_calo = etaGen_calo[iRecoJet_calo];
          Float_t iPhiGen_calo = phiGen_calo[iRecoJet_calo];
          Float_t iPtGen_calo = ptGen_calo[iRecoJet_calo];
          Float_t iEGen_calo = eGen_calo[iRecoJet_calo];

          Float_t iEtaReco_calo = etaReco_calo[iRecoJet_calo];
          Float_t iPhiReco_calo = phiReco_calo[iRecoJet_calo];
          Float_t iThetaReco_calo = 2*atan(exp(-iEtaReco_calo));
          Float_t iPtReco_calo = ptReco_calo[iRecoJet_calo];
          Float_t iEReco_calo = eReco_calo[iRecoJet_calo];

          Float_t deltaEta_calo = iEtaReco_calo-iEtaGen_calo;
          Float_t deltaPhi_calo = iPhiReco_calo-iPhiGen_calo;

          Float_t deltaEta_PF2calo = iEtaReco_PF-iEtaReco_calo;
          Float_t deltaPhi_PF2calo = iPhiReco_PF-iPhiReco_calo;

          Float_t deltaR_PF2calo = sqrt( deltaEta_PF2calo*deltaEta_PF2calo + deltaPhi_PF2calo*deltaPhi_PF2calo );

          if( deltaR_PF2calo<0.1 ) { //matched

            Float_t response_calo = iPtReco_calo/iPtGen_calo; 
            Float_t responseE_calo = iEReco_calo/iEGen_calo; 

            if( iEtaGen_PF != iEtaGen_calo ) continue; //should be same gen Jet

            h2_response_vs_pt_PF->Fill(iPtGen_PF, response_PF);
            h2_responseE_vs_E_PF->Fill(iEGen_PF, responseE_PF);
            h2_deltaEta_vs_pt_PF->Fill(iPtGen_PF, deltaEta_PF);
            h2_deltaPhi_vs_pt_PF->Fill(iPtGen_PF, deltaPhi_PF);

            h2_response_vs_pt_calo->Fill(iPtGen_calo, response_calo);
            h2_responseE_vs_E_calo->Fill(iEGen_calo, responseE_calo);
            h2_deltaEta_vs_pt_calo->Fill(iPtGen_calo, deltaEta_calo);
            h2_deltaPhi_vs_pt_calo->Fill(iPtGen_calo, deltaPhi_calo);


            if( fabs(iEtaGen_PF)<1.4 ) {
              h2_response_vs_pt_barrel_PF->Fill(iPtGen_PF, response_PF);
              h2_responseE_vs_E_barrel_PF->Fill(iEGen_PF, responseE_PF);
              h2_deltaEta_vs_pt_barrel_PF->Fill(iPtGen_PF, deltaEta_PF);
              h2_deltaPhi_vs_pt_barrel_PF->Fill(iPtGen_PF, deltaPhi_PF);

              h2_response_vs_pt_barrel_calo->Fill(iPtGen_calo, response_calo);
              h2_responseE_vs_E_barrel_calo->Fill(iEGen_calo, responseE_calo);
              h2_deltaEta_vs_pt_barrel_calo->Fill(iPtGen_calo, deltaEta_calo);
              h2_deltaPhi_vs_pt_barrel_calo->Fill(iPtGen_calo, deltaPhi_calo);

            } //if barrel

            if( (fabs(iEtaGen_PF)>1.4)&&(iEtaGen_PF<3.) ) { //endcaps

              h2_response_vs_pt_endcap_PF->Fill(iPtGen_PF, response_PF);
              h2_responseE_vs_E_endcap_PF->Fill(iEGen_PF, responseE_PF);
              h2_deltaEta_vs_pt_endcap_PF->Fill(iPtGen_PF, deltaEta_PF);
              h2_deltaPhi_vs_pt_endcap_PF->Fill(iPtGen_PF, deltaPhi_PF);

              h2_response_vs_pt_endcap_calo->Fill(iPtGen_calo, response_calo);
              h2_responseE_vs_E_endcap_calo->Fill(iEGen_calo, responseE_calo);
              h2_deltaEta_vs_pt_endcap_calo->Fill(iPtGen_calo, deltaEta_calo);
              h2_deltaPhi_vs_pt_endcap_calo->Fill(iPtGen_calo, deltaPhi_calo);

            } //if endcaps
             
            if( (fabs(iEtaGen_PF)>3.)&&(iEtaGen_PF<5.) ) { //forward

              h2_response_vs_pt_forward_PF->Fill(iPtGen_PF, response_PF);
              h2_deltaEta_vs_pt_forward_PF->Fill(iPtGen_PF, deltaEta_PF);
              h2_deltaPhi_vs_pt_forward_PF->Fill(iPtGen_PF, deltaPhi_PF);

              h2_response_vs_pt_forward_calo->Fill(iPtGen_calo, response_calo);
              h2_deltaEta_vs_pt_forward_calo->Fill(iPtGen_calo, deltaEta_calo);
              h2_deltaPhi_vs_pt_forward_calo->Fill(iPtGen_calo, deltaPhi_calo);

            } //if forward
   

          } //if matched

        } //for calo jets

      }  //if matched

    } //for PF jets

  } //for entries


  // -------- calo jet final plots:

  TH1F* h1_response_vs_pt_FIT_calo = new TH1F("response_vs_pt_FIT_calo", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_FIT_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_calo->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_response_vs_pt_FIT_barrel_calo = new TH1F("response_vs_pt_FIT_barrel_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_calo->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_response_vs_pt_FIT_endcap_calo = new TH1F("response_vs_pt_FIT_endcap_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_FIT_forward_calo = new TH1F("response_vs_pt_FIT_forward_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_forward_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_forward_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_MEAN_calo = new TH1F("response_vs_pt_MEAN_calo", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_MEAN_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_MEAN_barrel_calo = new TH1F("response_vs_pt_MEAN_barrel_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_calo->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_response_vs_pt_MEAN_endcap_calo = new TH1F("response_vs_pt_MEAN_endcap_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_MEAN_forward_calo = new TH1F("response_vs_pt_MEAN_forward_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_forward_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_forward_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_resolution_vs_pt_FIT_calo = new TH1F("resolution_vs_pt_FIT_calo", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_FIT_barrel_calo = new TH1F("resolution_vs_pt_FIT_barrel_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_FIT_endcap_calo = new TH1F("resolution_vs_pt_FIT_endcap_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_endcap_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_FIT_forward_calo = new TH1F("resolution_vs_pt_FIT_forward_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_forward_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_forward_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_calo = new TH1F("resolution_vs_pt_RMS_calo", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_barrel_calo = new TH1F("resolution_vs_pt_RMS_barrel_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_endcap_calo = new TH1F("resolution_vs_pt_RMS_endcap_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_forward_calo = new TH1F("resolution_vs_pt_RMS_forward_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_forward_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_forward_calo->GetYaxis()->SetRangeUser(0., 0.5);


  TH1F* h1_responseE_vs_E_FIT_calo = new TH1F("responseE_vs_E_FIT_calo", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_responseE_vs_E_FIT_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_calo->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_responseE_vs_E_FIT_barrel_calo = new TH1F("responseE_vs_E_FIT_barrel_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_FIT_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_barrel_calo->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_responseE_vs_E_FIT_endcap_calo = new TH1F("responseE_vs_E_FIT_endcap_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_FIT_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_endcap_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_endcap_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_responseE_vs_E_FIT_forward_calo = new TH1F("responseE_vs_E_FIT_forward_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_FIT_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_forward_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_forward_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_responseE_vs_E_MEAN_calo = new TH1F("responseE_vs_E_MEAN_calo", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_responseE_vs_E_MEAN_barrel_calo = new TH1F("responseE_vs_E_MEAN_barrel_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_barrel_calo->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_responseE_vs_E_MEAN_endcap_calo = new TH1F("responseE_vs_E_MEAN_endcap_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_endcap_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_endcap_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_responseE_vs_E_MEAN_forward_calo = new TH1F("responseE_vs_E_MEAN_forward_calo", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_forward_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_forward_calo->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_resolutionE_vs_E_FIT_calo = new TH1F("resolutionE_vs_E_FIT_calo", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_FIT_barrel_calo = new TH1F("resolutionE_vs_E_FIT_barrel_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_barrel_calo->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_barrel_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_FIT_endcap_calo = new TH1F("resolutionE_vs_E_FIT_endcap_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_endcap_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_endcap_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_FIT_forward_calo = new TH1F("resolutionE_vs_E_FIT_forward_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_forward_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_forward_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_calo = new TH1F("resolutionE_vs_E_RMS_calo", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_barrel_calo = new TH1F("resolutionE_vs_E_RMS_barrel_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_barrel_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_barrel_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_endcap_calo = new TH1F("resolutionE_vs_E_RMS_endcap_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_endcap_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_endcap_calo->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_forward_calo = new TH1F("resolutionE_vs_E_RMS_forward_calo", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_forward_calo->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_forward_calo->GetYaxis()->SetRangeUser(0., 0.5);



  // --------- PF plots:
  TH1F* h1_response_vs_pt_FIT_PF = new TH1F("response_vs_pt_FIT_PF", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_FIT_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_FIT_barrel_PF = new TH1F("response_vs_pt_FIT_barrel_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_barrel_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_barrel_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_FIT_endcap_PF = new TH1F("response_vs_pt_FIT_endcap_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_endcap_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_endcap_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_FIT_forward_PF = new TH1F("response_vs_pt_FIT_forward_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_FIT_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_FIT_forward_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_FIT_forward_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_MEAN_PF = new TH1F("response_vs_pt_MEAN_PF", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_response_vs_pt_MEAN_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_MEAN_barrel_PF = new TH1F("response_vs_pt_MEAN_barrel_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_barrel_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_barrel_PF->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_response_vs_pt_MEAN_endcap_PF = new TH1F("response_vs_pt_MEAN_endcap_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_endcap_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_endcap_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_response_vs_pt_MEAN_forward_PF = new TH1F("response_vs_pt_MEAN_forward_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_response_vs_pt_MEAN_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_response_vs_pt_MEAN_forward_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_response_vs_pt_MEAN_forward_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_resolution_vs_pt_FIT_PF = new TH1F("resolution_vs_pt_FIT_PF", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_FIT_barrel_PF = new TH1F("resolution_vs_pt_FIT_barrel_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_barrel_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_barrel_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_FIT_endcap_PF = new TH1F("resolution_vs_pt_FIT_endcap_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_endcap_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_endcap_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_FIT_forward_PF = new TH1F("resolution_vs_pt_FIT_forward_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_FIT_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_FIT_forward_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_FIT_forward_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_PF = new TH1F("resolution_vs_pt_RMS_PF", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_barrel_PF = new TH1F("resolution_vs_pt_RMS_barrel_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_barrel_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_barrel_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_endcap_PF = new TH1F("resolution_vs_pt_RMS_endcap_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_endcap_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_endcap_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolution_vs_pt_RMS_forward_PF = new TH1F("resolution_vs_pt_RMS_forward_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolution_vs_pt_RMS_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolution_vs_pt_RMS_forward_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolution_vs_pt_RMS_forward_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_responseE_vs_E_FIT_PF = new TH1F("responseE_vs_E_FIT_PF", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_responseE_vs_E_FIT_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_PF->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_responseE_vs_E_FIT_barrel_PF = new TH1F("responseE_vs_E_FIT_barrel_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_FIT_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_barrel_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_barrel_PF->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_responseE_vs_E_FIT_endcap_PF = new TH1F("responseE_vs_E_FIT_endcap_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_FIT_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_FIT_endcap_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_FIT_endcap_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_responseE_vs_E_MEAN_PF = new TH1F("responseE_vs_E_MEAN_PF", "Jet Response vs. True Jet p_{T}", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_responseE_vs_E_MEAN_barrel_PF = new TH1F("responseE_vs_E_MEAN_barrel_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_barrel_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_barrel_PF->GetYaxis()->SetRangeUser(0.2, 1.);
  
  TH1F* h1_responseE_vs_E_MEAN_endcap_PF = new TH1F("responseE_vs_E_MEAN_endcap_PF", "Jet Response vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_responseE_vs_E_MEAN_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_responseE_vs_E_MEAN_endcap_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");
  h1_responseE_vs_E_MEAN_endcap_PF->GetYaxis()->SetRangeUser(0.2, 1.);

  TH1F* h1_resolutionE_vs_E_FIT_PF = new TH1F("resolutionE_vs_E_FIT_PF", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_FIT_barrel_PF = new TH1F("resolutionE_vs_E_FIT_barrel_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_barrel_PF->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_barrel_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_FIT_endcap_PF = new TH1F("resolutionE_vs_E_FIT_endcap_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_FIT_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_FIT_endcap_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_FIT_endcap_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_PF = new TH1F("resolutionE_vs_E_RMS_PF", "Jet resolution vs. True Jet p_{T}", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_barrel_PF = new TH1F("resolutionE_vs_E_RMS_barrel_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_barrel_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_barrel_PF->GetYaxis()->SetRangeUser(0., 0.5);

  TH1F* h1_resolutionE_vs_E_RMS_endcap_PF = new TH1F("resolutionE_vs_E_RMS_endcap_PF", "Jet resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_resolutionE_vs_E_RMS_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_resolutionE_vs_E_RMS_endcap_PF->SetYTitle("#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]");
  h1_resolutionE_vs_E_RMS_endcap_PF->GetYaxis()->SetRangeUser(0., 0.5);



  fillProfile(h1_response_vs_pt_FIT_calo, h1_resolution_vs_pt_FIT_calo, h1_response_vs_pt_MEAN_calo, h1_resolution_vs_pt_RMS_calo, h2_response_vs_pt_calo);
  fillProfile(h1_response_vs_pt_FIT_PF, h1_resolution_vs_pt_FIT_PF, h1_response_vs_pt_MEAN_PF, h1_resolution_vs_pt_RMS_PF, h2_response_vs_pt_PF);

  fillProfile(h1_response_vs_pt_FIT_barrel_calo, h1_resolution_vs_pt_FIT_barrel_calo, h1_response_vs_pt_MEAN_barrel_calo, h1_resolution_vs_pt_RMS_barrel_calo, h2_response_vs_pt_barrel_calo);
  fillProfile(h1_response_vs_pt_FIT_endcap_calo, h1_resolution_vs_pt_FIT_endcap_calo, h1_response_vs_pt_MEAN_endcap_calo, h1_resolution_vs_pt_RMS_endcap_calo, h2_response_vs_pt_endcap_calo);
  fillProfile(h1_response_vs_pt_FIT_forward_calo, h1_resolution_vs_pt_FIT_forward_calo, h1_response_vs_pt_MEAN_forward_calo, h1_resolution_vs_pt_RMS_forward_calo, h2_response_vs_pt_forward_calo);

  fillProfile(h1_response_vs_pt_FIT_barrel_PF, h1_resolution_vs_pt_FIT_barrel_PF, h1_response_vs_pt_MEAN_barrel_PF, h1_resolution_vs_pt_RMS_barrel_PF, h2_response_vs_pt_barrel_PF);
  fillProfile(h1_response_vs_pt_FIT_endcap_PF, h1_resolution_vs_pt_FIT_endcap_PF, h1_response_vs_pt_MEAN_endcap_PF, h1_resolution_vs_pt_RMS_endcap_PF, h2_response_vs_pt_endcap_PF);
  fillProfile(h1_response_vs_pt_FIT_forward_PF, h1_resolution_vs_pt_FIT_forward_PF, h1_response_vs_pt_MEAN_forward_PF, h1_resolution_vs_pt_RMS_forward_PF, h2_response_vs_pt_forward_PF);


  fillProfile(h1_responseE_vs_E_FIT_calo, h1_resolutionE_vs_E_FIT_calo, h1_responseE_vs_E_MEAN_calo, h1_resolutionE_vs_E_RMS_calo, h2_responseE_vs_E_calo);
  fillProfile(h1_responseE_vs_E_FIT_PF, h1_resolutionE_vs_E_FIT_PF, h1_responseE_vs_E_MEAN_PF, h1_resolutionE_vs_E_RMS_PF, h2_responseE_vs_E_PF);

  fillProfile(h1_responseE_vs_E_FIT_barrel_calo, h1_resolutionE_vs_E_FIT_barrel_calo, h1_responseE_vs_E_MEAN_barrel_calo, h1_resolutionE_vs_E_RMS_barrel_calo, h2_responseE_vs_E_barrel_calo);
  fillProfile(h1_responseE_vs_E_FIT_endcap_calo, h1_resolutionE_vs_E_FIT_endcap_calo, h1_responseE_vs_E_MEAN_endcap_calo, h1_resolutionE_vs_E_RMS_endcap_calo, h2_responseE_vs_E_endcap_calo);

  fillProfile(h1_responseE_vs_E_FIT_barrel_PF, h1_resolutionE_vs_E_FIT_barrel_PF, h1_responseE_vs_E_MEAN_barrel_PF, h1_resolutionE_vs_E_RMS_barrel_PF, h2_responseE_vs_E_barrel_PF);
  fillProfile(h1_responseE_vs_E_FIT_endcap_PF, h1_resolutionE_vs_E_FIT_endcap_PF, h1_responseE_vs_E_MEAN_endcap_PF, h1_resolutionE_vs_E_RMS_endcap_PF, h2_responseE_vs_E_endcap_PF);

 
  // ######## Position resolution


  // -------- calo jet plots:
  
  TH1F* h1_sigmaEta_vs_pt_calo = new TH1F("sigmaEta_vs_pt_calo", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_calo->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_calo->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaEta_vs_pt_barrel_calo = new TH1F("sigmaEta_vs_pt_barrel_calo", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_calo->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_calo->GetYaxis()->SetRangeUser(0., 0.2);
 
  TH1F* h1_sigmaEta_vs_pt_endcap_calo = new TH1F("sigmaEta_vs_pt_endcap_calo", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_calo->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_calo->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaEta_vs_pt_forward_calo = new TH1F("sigmaEta_vs_pt_forward_calo", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_forward_calo->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_forward_calo->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_calo = new TH1F("sigmaPhi_vs_pt_calo", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_calo->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_calo->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_barrel_calo = new TH1F("sigmaPhi_vs_pt_barrel_calo", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_calo->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_barrel_calo->GetYaxis()->SetRangeUser(0., 0.2);
 
  TH1F* h1_sigmaPhi_vs_pt_endcap_calo = new TH1F("sigmaPhi_vs_pt_endcap_calo", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_calo->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_endcap_calo->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_forward_calo = new TH1F("sigmaPhi_vs_pt_forward_calo", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_forward_calo->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_forward_calo->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_forward_calo->GetYaxis()->SetRangeUser(0., 0.2);


  // -------- PF jet plots:

  TH1F* h1_sigmaEta_vs_pt_PF = new TH1F("sigmaEta_vs_pt_PF", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_PF->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaEta_vs_pt_barrel_PF = new TH1F("sigmaEta_vs_pt_barrel_PF", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_barrel_PF->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_barrel_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaEta_vs_pt_endcap_PF = new TH1F("sigmaEta_vs_pt_endcap_PF", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_endcap_PF->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_endcap_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaEta_vs_pt_forward_PF = new TH1F("sigmaEta_vs_pt_forward_PF", "Jet Position #eta Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaEta_vs_pt_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaEta_vs_pt_forward_PF->SetYTitle("#sigma(#Delta#eta)");
  h1_sigmaEta_vs_pt_forward_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_PF = new TH1F("sigmaPhi_vs_pt_PF", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_PF->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_barrel_PF = new TH1F("sigmaPhi_vs_pt_barrel_PF", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_barrel_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_barrel_PF->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_barrel_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_endcap_PF = new TH1F("sigmaPhi_vs_pt_endcap_PF", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_endcap_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_endcap_PF->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_endcap_PF->GetYaxis()->SetRangeUser(0., 0.2);

  TH1F* h1_sigmaPhi_vs_pt_forward_PF = new TH1F("sigmaPhi_vs_pt_forward_PF", "Jet Position #phi Resolution vs. True Jet p_{T} (0.0<|#phi_{GEN}|<1.4)", nBins-2, Lower);
  h1_sigmaPhi_vs_pt_forward_PF->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_sigmaPhi_vs_pt_forward_PF->SetYTitle("#sigma(#Delta#phi)");
  h1_sigmaPhi_vs_pt_forward_PF->GetYaxis()->SetRangeUser(0., 0.2);


  fillPositionResolution(h1_sigmaEta_vs_pt_calo, h1_sigmaPhi_vs_pt_calo, h2_deltaEta_vs_pt_calo, h2_deltaPhi_vs_pt_calo);
  fillPositionResolution(h1_sigmaEta_vs_pt_PF, h1_sigmaPhi_vs_pt_PF, h2_deltaEta_vs_pt_PF, h2_deltaPhi_vs_pt_PF);

  fillPositionResolution(h1_sigmaEta_vs_pt_barrel_calo, h1_sigmaPhi_vs_pt_barrel_calo, h2_deltaEta_vs_pt_barrel_calo, h2_deltaPhi_vs_pt_barrel_calo);
  fillPositionResolution(h1_sigmaEta_vs_pt_endcap_calo, h1_sigmaPhi_vs_pt_endcap_calo, h2_deltaEta_vs_pt_endcap_calo, h2_deltaPhi_vs_pt_endcap_calo);
  fillPositionResolution(h1_sigmaEta_vs_pt_forward_calo, h1_sigmaPhi_vs_pt_forward_calo, h2_deltaEta_vs_pt_forward_calo, h2_deltaPhi_vs_pt_forward_calo);

  fillPositionResolution(h1_sigmaEta_vs_pt_barrel_PF, h1_sigmaPhi_vs_pt_barrel_PF, h2_deltaEta_vs_pt_barrel_PF, h2_deltaPhi_vs_pt_barrel_PF);
  fillPositionResolution(h1_sigmaEta_vs_pt_endcap_PF, h1_sigmaPhi_vs_pt_endcap_PF, h2_deltaEta_vs_pt_endcap_PF, h2_deltaPhi_vs_pt_endcap_PF);
  fillPositionResolution(h1_sigmaEta_vs_pt_forward_PF, h1_sigmaPhi_vs_pt_forward_PF, h2_deltaEta_vs_pt_forward_PF, h2_deltaPhi_vs_pt_forward_PF);


 
  delete h2_response_vs_pt_calo;
  h2_response_vs_pt_calo = 0;
  delete h2_response_vs_pt_barrel_calo;
  h2_response_vs_pt_barrel_calo = 0;
  delete h2_response_vs_pt_endcap_calo;
  h2_response_vs_pt_endcap_calo = 0;
  delete h2_response_vs_pt_forward_calo;
  h2_response_vs_pt_forward_calo = 0;

  delete h2_response_vs_pt_PF;
  h2_response_vs_pt_PF = 0;
  delete h2_response_vs_pt_barrel_PF;
  h2_response_vs_pt_barrel_PF = 0;
  delete h2_response_vs_pt_endcap_PF;
  h2_response_vs_pt_endcap_PF = 0;
  delete h2_response_vs_pt_forward_PF;
  h2_response_vs_pt_forward_PF = 0;

  delete h2_deltaEta_vs_pt_calo;
  h2_deltaEta_vs_pt_calo = 0;
  delete h2_deltaEta_vs_pt_barrel_calo;
  h2_deltaEta_vs_pt_barrel_calo = 0;
  delete h2_deltaEta_vs_pt_endcap_calo;
  h2_deltaEta_vs_pt_endcap_calo = 0;
  delete h2_deltaEta_vs_pt_forward_calo;
  h2_deltaEta_vs_pt_forward_calo = 0;

  delete h2_deltaPhi_vs_pt_calo;
  h2_deltaPhi_vs_pt_calo = 0;
  delete h2_deltaPhi_vs_pt_barrel_calo;
  h2_deltaPhi_vs_pt_barrel_calo = 0;
  delete h2_deltaPhi_vs_pt_endcap_calo;
  h2_deltaPhi_vs_pt_endcap_calo = 0;
  delete h2_deltaPhi_vs_pt_forward_calo;
  h2_deltaPhi_vs_pt_forward_calo = 0;

  delete h2_deltaEta_vs_pt_PF;
  h2_deltaEta_vs_pt_PF = 0;
  delete h2_deltaEta_vs_pt_barrel_PF;
  h2_deltaEta_vs_pt_barrel_PF = 0;
  delete h2_deltaEta_vs_pt_endcap_PF;
  h2_deltaEta_vs_pt_endcap_PF = 0;
  delete h2_deltaEta_vs_pt_forward_PF;
  h2_deltaEta_vs_pt_forward_PF = 0;

  delete h2_deltaPhi_vs_pt_PF;
  h2_deltaPhi_vs_pt_PF = 0;
  delete h2_deltaPhi_vs_pt_barrel_PF;
  h2_deltaPhi_vs_pt_barrel_PF = 0;
  delete h2_deltaPhi_vs_pt_endcap_PF;
  h2_deltaPhi_vs_pt_endcap_PF = 0;
  delete h2_deltaPhi_vs_pt_forward_PF;
  h2_deltaPhi_vs_pt_forward_PF = 0;


  outFile->cd();


  h1_response_vs_pt_MEAN_calo->Write();
  h1_response_vs_pt_MEAN_barrel_calo->Write();
  h1_response_vs_pt_MEAN_endcap_calo->Write();
  h1_response_vs_pt_MEAN_forward_calo->Write();

  h1_resolution_vs_pt_RMS_calo->Write();
  h1_resolution_vs_pt_RMS_barrel_calo->Write();
  h1_resolution_vs_pt_RMS_endcap_calo->Write();
  h1_resolution_vs_pt_RMS_forward_calo->Write();

  h1_response_vs_pt_FIT_calo->Write();
  h1_response_vs_pt_FIT_barrel_calo->Write();
  h1_response_vs_pt_FIT_endcap_calo->Write();
  h1_response_vs_pt_FIT_forward_calo->Write();

  h1_resolution_vs_pt_FIT_calo->Write();
  h1_resolution_vs_pt_FIT_barrel_calo->Write();
  h1_resolution_vs_pt_FIT_endcap_calo->Write();
  h1_resolution_vs_pt_FIT_forward_calo->Write();

  h1_responseE_vs_E_MEAN_calo->Write();
  h1_responseE_vs_E_MEAN_barrel_calo->Write();
  h1_responseE_vs_E_MEAN_endcap_calo->Write();
  h1_responseE_vs_E_MEAN_forward_calo->Write();

  h1_resolutionE_vs_E_RMS_calo->Write();
  h1_resolutionE_vs_E_RMS_barrel_calo->Write();
  h1_resolutionE_vs_E_RMS_endcap_calo->Write();
  h1_resolutionE_vs_E_RMS_forward_calo->Write();

  h1_responseE_vs_E_FIT_calo->Write();
  h1_responseE_vs_E_FIT_barrel_calo->Write();
  h1_responseE_vs_E_FIT_endcap_calo->Write();
  h1_responseE_vs_E_FIT_forward_calo->Write();

  h1_resolutionE_vs_E_FIT_calo->Write();
  h1_resolutionE_vs_E_FIT_barrel_calo->Write();
  h1_resolutionE_vs_E_FIT_endcap_calo->Write();
  h1_resolutionE_vs_E_FIT_forward_calo->Write();

  h1_sigmaEta_vs_pt_calo->Write();
  h1_sigmaEta_vs_pt_barrel_calo->Write();
  h1_sigmaEta_vs_pt_endcap_calo->Write();
  h1_sigmaEta_vs_pt_forward_calo->Write();

  h1_sigmaPhi_vs_pt_calo->Write();
  h1_sigmaPhi_vs_pt_barrel_calo->Write();
  h1_sigmaPhi_vs_pt_endcap_calo->Write();
  h1_sigmaPhi_vs_pt_forward_calo->Write();

  h1_response_vs_pt_MEAN_PF->Write();
  h1_response_vs_pt_MEAN_barrel_PF->Write();
  h1_response_vs_pt_MEAN_endcap_PF->Write();
  h1_response_vs_pt_MEAN_forward_PF->Write();

  h1_resolution_vs_pt_RMS_PF->Write();
  h1_resolution_vs_pt_RMS_barrel_PF->Write();
  h1_resolution_vs_pt_RMS_endcap_PF->Write();
  h1_resolution_vs_pt_RMS_forward_PF->Write();

  h1_response_vs_pt_FIT_PF->Write();
  h1_response_vs_pt_FIT_barrel_PF->Write();
  h1_response_vs_pt_FIT_endcap_PF->Write();
  h1_response_vs_pt_FIT_forward_PF->Write();

  h1_resolution_vs_pt_FIT_PF->Write();
  h1_resolution_vs_pt_FIT_barrel_PF->Write();
  h1_resolution_vs_pt_FIT_endcap_PF->Write();
  h1_resolution_vs_pt_FIT_forward_PF->Write();

  h1_sigmaEta_vs_pt_PF->Write();
  h1_sigmaEta_vs_pt_barrel_PF->Write();
  h1_sigmaEta_vs_pt_endcap_PF->Write();
  h1_sigmaEta_vs_pt_forward_PF->Write();

  h1_sigmaPhi_vs_pt_PF->Write();
  h1_sigmaPhi_vs_pt_barrel_PF->Write();
  h1_sigmaPhi_vs_pt_endcap_PF->Write();
  h1_sigmaPhi_vs_pt_forward_PF->Write();

  
  outFile->Write();
  outFile->Close();


}//Rch0

