#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TChain.h>
#include <iostream>




void drawPtSpectrum() {

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

  chain->Add("2ndLevelOutputFile_PhotonJet15.root");
  chain->Add("2ndLevelOutputFile_PhotonJet30.root");
  chain->Add("2ndLevelOutputFile_PhotonJet80.root");
  chain->Add("2ndLevelOutputFile_PhotonJet170.root");
  chain->Add("2ndLevelOutputFile_PhotonJet300.root");
  chain->Add("2ndLevelOutputFile_PhotonJet470.root");
  chain->Add("2ndLevelOutputFile_PhotonJet800.root");

  Float_t eventWeight;
  chain->SetBranchAddress("eventWeight", &eventWeight);

  Float_t ptHat;
  chain->SetBranchAddress("ptHat", &ptHat);

  TH1F* h1_ptSpectrum = new TH1F("ptSpectrum", "", 100, 0., 1000.);
  h1_ptSpectrum->SetXTitle("#hat{p}_{T} [GeV/c]");
  h1_ptSpectrum->SetMarkerStyle(20);
  h1_ptSpectrum->SetMarkerSize(1.4);



  Int_t nEntries = chain->GetEntries();


  for( Int_t iEntry = 0; iEntry<nEntries; ++iEntry) {
 
    if( (iEntry%100000) == 0 ) std::cout << "Entry: " << iEntry << "/" << nEntries << std::endl;
    chain->GetEntry(iEntry);

    h1_ptSpectrum->Fill(ptHat, eventWeight);

  }

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogy();
  h1_ptSpectrum->Draw("P");
  c1->SaveAs("ptSpectrum.eps");

}
