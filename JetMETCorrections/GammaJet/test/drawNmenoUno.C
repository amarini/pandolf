#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TChain.h>
#include <iostream>


void fillEfficiencyPlot( TH1F* h1_num, TH1F* h1_denom, TH1F* h1_eff );


void drawNmenoUno() {

  TStyle *simpleStyle = new TStyle("simpleStyle","");
  simpleStyle->SetCanvasColor(0);
  simpleStyle->SetFrameFillColor(0);
  simpleStyle->SetStatColor(0);
  simpleStyle->SetOptStat(0);
  simpleStyle->SetTitleFillColor(0);
  simpleStyle->SetCanvasBorderMode(0);
  simpleStyle->SetPadBorderMode(0);
  simpleStyle->SetFrameBorderMode(0);
  simpleStyle->cd();

  TChain* chain = new TChain("jetTree");

  //chain->Add("2ndLevelOutputFile_PhotonJet15.root");
  //chain->Add("2ndLevelOutputFile_PhotonJet30.root");
  //chain->Add("2ndLevelOutputFile_PhotonJet80.root");
  chain->Add("2ndLevelOutputFile_PhotonJet170.root");
  //chain->Add("2ndLevelOutputFile_PhotonJet300.root");
  //chain->Add("2ndLevelOutputFile_PhotonJet470.root");
  //chain->Add("2ndLevelOutputFile_PhotonJet800.root");

  Float_t eventWeight;
  chain->SetBranchAddress("eventWeight", &eventWeight);

  Float_t hcalIso;
  chain->SetBranchAddress("hcalIsoPhotReco", &hcalIso);

  Float_t ecalIso;
  chain->SetBranchAddress("ecalIsoPhotReco", &ecalIso);

  Int_t nTrkIso;
  chain->SetBranchAddress("nTrkIsoPhotReco", &nTrkIso);

  Float_t ptTrkIso;
  chain->SetBranchAddress("ptTrkIsoPhotReco", &ptTrkIso);

  Float_t clusterMaj;
  chain->SetBranchAddress("clusterMajPhotReco", &clusterMaj);

  Float_t clusterMin;
  chain->SetBranchAddress("clusterMinPhotReco", &clusterMin);

  Float_t ptPhot;
  chain->SetBranchAddress("ptPhotReco", &ptPhot);

  Float_t etaPhot;
  chain->SetBranchAddress("etaPhotReco", &etaPhot);

  Float_t phiPhot;
  chain->SetBranchAddress("phiPhotReco", &phiPhot);

  Float_t pt2ndJet;
  chain->SetBranchAddress("pt2ndJetReco", &pt2ndJet);

  Float_t phi2ndJet;
  chain->SetBranchAddress("phi2ndJetReco", &phi2ndJet);

  Float_t eta2ndJet;
  chain->SetBranchAddress("eta2ndJetReco", &eta2ndJet);

  Float_t phiJet;
  chain->SetBranchAddress("phiJetReco", &phiJet);


  TH1F* h1_hcalIso = new TH1F("hcalIso", "", 100, -0.2, 1.4 );
  TH1F* h1_ecalIso = new TH1F("ecalIso", "", 100, -0.2, 1.4 );
  TH1F* h1_nTrkIso = new TH1F("nTrkIso", "", 26, -0.5, 25.5 );
  TH1F* h1_ptTrkIso = new TH1F("ptTrkIso", "", 100, 0., 5. );
  TH1F* h1_clusterMaj = new TH1F("clusterMaj", "", 100, 0., 2. );
  TH1F* h1_clusterMin = new TH1F("clusterMin", "", 100, 0., 1.4 );
  //TH1F* h1_ptPhot = new TH1F("ptPhot", "", 100, 0., 70. );
  //TH1F* h1_etaPhot = new TH1F("etaPhot", "", 100, -3.5, 3.5 );
  TH1F* h1_pt2ndJet = new TH1F("pt2ndJet", "", 100, 0., 5.);
  TH1F* h1_deltaPhi = new TH1F("deltaPhi", "", 100, 0., 3.14159);

  TH1F* h1_deltaR = new TH1F("deltaR", "", 100, 0., 1.);


  //efficiency plots:
  Int_t nBins_pt = 15;
  Float_t ptMin = 0.;
  Float_t ptMax = 700.;
  TH1F* h1_ptPhot_denom = new TH1F("ptPhot_denom", "", nBins_pt, ptMin, ptMax);
  h1_ptPhot_denom->SetMarkerStyle(20);
  TH1F* h1_ptPhot_nTrkIso_abs = new TH1F("ptPhot_nTrkIso_abs", "", nBins_pt, ptMin, ptMax);
  h1_ptPhot_nTrkIso_abs->SetMarkerStyle(20);
  TH1F* h1_ptPhot_nTrkIso_cascade = new TH1F("ptPhot_nTrkIso_cascade", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_ptTrkIso_abs = new TH1F("ptPhot_ptTrkIso_abs", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_ptTrkIso_cascade = new TH1F("ptPhot_ptTrkIso_cascade", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_hcalIso_abs = new TH1F("ptPhot_hcalIso_abs", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_hcalIso_cascade = new TH1F("ptPhot_hcalIso_cascade", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_clusterMaj_abs = new TH1F("ptPhot_clusterMaj_abs", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_clusterMaj_cascade = new TH1F("ptPhot_clusterMaj_cascade", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_clusterMin_abs = new TH1F("ptPhot_clusterMin_abs", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_clusterMin_cascade = new TH1F("ptPhot_clusterMin_cascade", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_ecalIso_abs = new TH1F("ptPhot_ecalIso_abs", "", nBins_pt, ptMin, ptMax);
  TH1F* h1_ptPhot_ecalIso_cascade = new TH1F("ptPhot_ecalIso_cascade", "", nBins_pt, ptMin, ptMax);
  
  TH1F* h1_eff_nTrkIso_abs = new TH1F("eff_nTrkIso_abs", "", nBins_pt, ptMin, ptMax);
  h1_eff_nTrkIso_abs->SetMarkerStyle(20);
  h1_eff_nTrkIso_abs->SetMarkerColor(kRed);
  TH1F* h1_eff_nTrkIso_cascade = new TH1F("eff_nTrkIso_cascade", "", nBins_pt, ptMin, ptMax);
  h1_eff_nTrkIso_cascade->SetMarkerStyle(20);
  h1_eff_nTrkIso_cascade->SetMarkerColor(kRed);
  TH1F* h1_eff_ptTrkIso_abs = new TH1F("eff_ptTrkIso_abs", "", nBins_pt, ptMin, ptMax);
  h1_eff_ptTrkIso_abs->SetMarkerStyle(21);
  h1_eff_ptTrkIso_abs->SetMarkerColor(kBlue);
  TH1F* h1_eff_ptTrkIso_cascade = new TH1F("eff_ptTrkIso_cascade", "", nBins_pt, ptMin, ptMax);
  h1_eff_ptTrkIso_cascade->SetMarkerStyle(21);
  h1_eff_ptTrkIso_cascade->SetMarkerColor(kBlue);
  TH1F* h1_eff_hcalIso_abs = new TH1F("eff_hcalIso_abs", "", nBins_pt, ptMin, ptMax);
  h1_eff_hcalIso_abs->SetMarkerStyle(22);
  h1_eff_hcalIso_abs->SetMarkerColor(kGreen);
  TH1F* h1_eff_hcalIso_cascade = new TH1F("eff_hcalIso_cascade", "", nBins_pt, ptMin, ptMax);
  h1_eff_hcalIso_cascade->SetMarkerStyle(22);
  h1_eff_hcalIso_cascade->SetMarkerColor(kGreen);
  TH1F* h1_eff_clusterMaj_abs = new TH1F("eff_clusterMaj_abs", "", nBins_pt, ptMin, ptMax);
  h1_eff_clusterMaj_abs->SetMarkerStyle(23);
  h1_eff_clusterMaj_abs->SetMarkerColor(kMagenta);
  TH1F* h1_eff_clusterMaj_cascade = new TH1F("eff_clusterMaj_cascade", "", nBins_pt, ptMin, ptMax);
  h1_eff_clusterMaj_cascade->SetMarkerStyle(23);
  h1_eff_clusterMaj_cascade->SetMarkerColor(kMagenta);
  TH1F* h1_eff_clusterMin_abs = new TH1F("eff_clusterMin_abs", "", nBins_pt, ptMin, ptMax);
  h1_eff_clusterMin_abs->SetMarkerStyle(24);
  h1_eff_clusterMin_abs->SetMarkerColor(kCyan);
  TH1F* h1_eff_clusterMin_cascade = new TH1F("eff_clusterMin_cascade", "", nBins_pt, ptMin, ptMax);
  h1_eff_clusterMin_cascade->SetMarkerStyle(24);
  h1_eff_clusterMin_cascade->SetMarkerColor(kCyan);
  TH1F* h1_eff_ecalIso_abs = new TH1F("eff_ecalIso_abs", "", nBins_pt, ptMin, ptMax);
  h1_eff_ecalIso_abs->SetMarkerStyle(25);
  h1_eff_ecalIso_abs->SetMarkerColor(kYellow);
  TH1F* h1_eff_ecalIso_cascade = new TH1F("eff_ecalIso_cascade", "", nBins_pt, ptMin, ptMax);
  h1_eff_ecalIso_cascade->SetMarkerStyle(25);
  h1_eff_ecalIso_cascade->SetMarkerColor(kYellow);
  


  Int_t nEntries = chain->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    chain->GetEntry(iEntry);

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << "/" << nEntries << std::endl;

    if( ptPhot < 25. ) continue;
    if( fabs(etaPhot)>1.3 ) continue;


    bool b_hcalIso = true;
    bool b_ecalIso = true;
    bool b_nTrkIso = true;
    bool b_ptTrkIso = true;
    bool b_clusterMaj= true;
    bool b_clusterMin = true;
    bool b_pt2ndJet = true;
    bool b_deltaPhi = true;

    h1_ptPhot_denom->Fill(ptPhot, eventWeight);

    if( nTrkIso>= 3) {
      b_nTrkIso = false;
      h1_ptPhot_nTrkIso_abs->Fill(ptPhot, eventWeight);
    }
    if( ptTrkIso>= 0.1) {
      b_ptTrkIso = false;
      h1_ptPhot_ptTrkIso_abs->Fill(ptPhot, eventWeight);
    }
    if( hcalIso>= 0.053) {
      b_hcalIso = false;
      h1_ptPhot_hcalIso_abs->Fill(ptPhot, eventWeight);
    }
    if( clusterMaj>= 0.3) {
      b_clusterMaj = false;
      h1_ptPhot_clusterMaj_abs->Fill(ptPhot, eventWeight);
    }
    if( clusterMin>= 0.25) {
      b_clusterMin = false;
      h1_ptPhot_clusterMin_abs->Fill(ptPhot, eventWeight);
    }
    if( ecalIso>= 0.05) {
      b_ecalIso = false;
      h1_ptPhot_ecalIso_abs->Fill(ptPhot, eventWeight);
    }

    if( b_nTrkIso && b_ptTrkIso )
      h1_ptPhot_ptTrkIso_cascade->Fill(ptPhot, eventWeight);

    if( b_nTrkIso && b_ptTrkIso && b_hcalIso )
      h1_ptPhot_hcalIso_cascade->Fill(ptPhot, eventWeight);

    if( b_nTrkIso && b_ptTrkIso && b_hcalIso && b_clusterMaj)
      h1_ptPhot_clusterMaj_cascade->Fill(ptPhot, eventWeight);

    if( b_nTrkIso && b_ptTrkIso && b_hcalIso && b_clusterMaj && b_clusterMin )
      h1_ptPhot_clusterMin_cascade->Fill(ptPhot, eventWeight);

    if( b_nTrkIso && b_ptTrkIso && b_hcalIso && b_clusterMaj && b_clusterMin && b_ecalIso )
      h1_ptPhot_ecalIso_cascade->Fill(ptPhot, eventWeight);



    if( (pt2ndJet>0.1*ptPhot)&&(pt2ndJet>10.)) b_pt2ndJet = false;
    
    //leading jet and photon back2back in transverse plane
    Float_t deltaPhi_jet = phiPhot- phiJet;
    Float_t pi = 3.14159;
    if( fabs(deltaPhi_jet) < (pi - 0.2) ) b_deltaPhi = false;

    Float_t deltaEta = etaPhot - eta2ndJet;
    Float_t deltaPhi = phiPhot - phi2ndJet;
    Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

    if( b_ecalIso && b_nTrkIso && b_ptTrkIso && b_clusterMaj && b_clusterMin ) 
      h1_hcalIso->Fill(hcalIso);

    if( b_hcalIso && b_nTrkIso && b_ptTrkIso && b_clusterMaj && b_clusterMin ) 
      h1_ecalIso->Fill(ecalIso);

    if( b_hcalIso && b_ecalIso && b_ptTrkIso && b_clusterMaj && b_clusterMin ) 
      h1_nTrkIso->Fill(nTrkIso);

    if( b_hcalIso && b_ecalIso && nTrkIso && b_clusterMaj && b_clusterMin) 
      h1_ptTrkIso->Fill(ptTrkIso);

    if( b_hcalIso && b_ecalIso && nTrkIso && b_ptTrkIso && b_clusterMin) 
      h1_clusterMaj->Fill(clusterMaj);

    if( b_hcalIso && b_ecalIso && nTrkIso && b_ptTrkIso && b_clusterMaj) 
      h1_clusterMin->Fill(clusterMin);

    if( b_hcalIso && b_ecalIso && b_nTrkIso && b_ptTrkIso && b_clusterMaj && b_clusterMin && b_pt2ndJet) {
      h1_deltaPhi->Fill( deltaPhi_jet, eventWeight );
    }

  } //for entries


  fillEfficiencyPlot( h1_ptPhot_nTrkIso_abs, h1_ptPhot_denom, h1_eff_nTrkIso_abs );

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h1_deltaPhi->Draw();
  c1->SaveAs("deltaPhi.eps");

  h1_ptPhot_denom->Draw("P");
  c1->SaveAs("ptPhot_denom.eps");

  h1_ptPhot_nTrkIso_abs->Draw("P");
  c1->SaveAs("ptPhot_nTrkIso.eps");

  h1_eff_nTrkIso_abs->Draw("P");
  c1->SaveAs("eff_nTrkIso.eps");


  c1->SetLogy();

  h1_hcalIso->Draw();
  c1->SaveAs("hcalIso.eps");

  h1_ecalIso->Draw();
  c1->SaveAs("ecalIso.eps");

  h1_nTrkIso->Draw();
  c1->SaveAs("nTrkIso.eps");

  h1_ptTrkIso->Draw();
  c1->SaveAs("ptTrkIso.eps");

  h1_clusterMaj->Draw();
  c1->SaveAs("clusterMaj.eps");

  h1_clusterMin->Draw();
  c1->SaveAs("clusterMin.eps");

}


void fillEfficiencyPlot( TH1F* h1_num, TH1F* h1_denom, TH1F* h1_eff ) {

  Int_t nBins = h1_denom->GetNbinsX();

  for( unsigned int iBin = 1; iBin<(nBins+1); ++iBin) {

    Float_t denom = h1_denom->GetBinContent(iBin);
    Float_t num = h1_num->GetBinContent(iBin);

    Float_t eff = ( denom > 0. ) ? num/denom : 0.;

    h1_eff->SetBinContent( iBin, eff );

  }

}
