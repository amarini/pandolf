#include "TCanvas.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"

void drawTracks() {


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

  TFile* dataFile = TFile::Open("Runs_123592_123596_123615_123732/output_1.root");
  TTree* dataTree = (TTree*)dataFile->Get("myanalysis/pippo");
  TFile* mcFile = TFile::Open("MC_MinBias_900GeV/output_12.root");
  TTree* mcTree = (TTree*)dataFile->Get("myanalysis/pippo");

std::cout << "1" << std::endl;
  Int_t nEventTracks_data;
  dataTree->SetBranchAddress("nEventTracks", &nEventTracks_data);
  Int_t nChargedHadrons_data;
  dataTree->SetBranchAddress("nChargedHadrons", &nChargedHadrons_data);

std::cout << "2" << std::endl;
  Int_t nEventTracks_mc;
  mcTree->SetBranchAddress("nEventTracks", &nEventTracks_mc);
  Int_t nChargedHadrons_mc;
  mcTree->SetBranchAddress("nChargedHadrons", &nChargedHadrons_mc);

std::cout << "3" << std::endl;
  TH1F* h1_data = new TH1F("data", "", 100, 0., 100.);
  h1_data->SetMarkerStyle(20);

  TH1F* h1_mc = new TH1F("mc", "", 100, 0., 100.);
  h1_mc->SetFillColor(kGray);

std::cout << "4" << std::endl;
  for( int iEntry=0; iEntry<dataTree->GetEntries(); ++iEntry) {
    dataTree->GetEntry(iEntry);
    Float_t trackRatio = (nChargedHadrons_data>0) ? (Float_t)(nEventTracks_data/nChargedHadrons_data) : 0.;
    h1_data->Fill( trackRatio );
  }

std::cout << "5" << std::endl;
  for( int iEntry=0; iEntry<mcTree->GetEntries(); ++iEntry) {
    Float_t trackRatio = (nChargedHadrons_mc>0) ? (Float_t)(nEventTracks_mc/nChargedHadrons_mc) : 0.;
    h1_mc->Fill( trackRatio );
  }

std::cout << "6" << std::endl;
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h1_data->Draw("E");
  h1_mc->Draw("histo same");
  c1->SaveAs("tracks.eps");
 
}

