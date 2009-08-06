#include "TProfile.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"


void Rx_vs_Rch() {

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

  TFile* file = TFile::Open("finalOutputFile_226_SAVE.root");

  TProfile* hp_Rgamma_vs_Rch = (TProfile*)file->Get("Rgamma_vs_Rch_PFItCone5");
  hp_Rgamma_vs_Rch->SetMarkerColor(kGreen);
  hp_Rgamma_vs_Rch->SetMarkerStyle(20);
  hp_Rgamma_vs_Rch->SetTitle("");
  hp_Rgamma_vs_Rch->GetYaxis()->SetTitle("");

  TProfile* hp_Rnh_vs_Rch = (TProfile*)file->Get("Rnh_vs_Rch_PFItCone5");
  hp_Rnh_vs_Rch->SetMarkerColor(kBlue);
  hp_Rnh_vs_Rch->SetMarkerStyle(21);
  hp_Rnh_vs_Rch->SetTitle("");
  hp_Rnh_vs_Rch->GetYaxis()->SetTitle("");


  Int_t nBins = hp_Rnh_vs_Rch->GetXaxis()->GetNbins();
  Float_t xMin = hp_Rnh_vs_Rch->GetXaxis()->GetXmin();
  Float_t xMax = hp_Rnh_vs_Rch->GetXaxis()->GetXmax();
  Float_t yMin = hp_Rnh_vs_Rch->GetYaxis()->GetXmin();
  Float_t yMax = hp_Rnh_vs_Rch->GetYaxis()->GetXmax();

  TH1F* hp_divide = new TH1F("divide", "", nBins, xMin, xMax, yMin, yMax);

  for(int iBin=1; iBin<nBins+1; ++iBin)
    hp_divide->SetBinContent(iBin, hp_Rnh_vs_Rch->GetBinContent(iBin)/hp_Rgamma_vs_Rch->GetBinContent(iBin));

  TLegend* legend = new TLegend(0.7, 0.7, 0.88, 0.85);
  legend->SetFillColor(kWhite);

  legend->AddEntry(hp_Rgamma_vs_Rch, "R_{#gamma}", "P");
  legend->AddEntry(hp_Rnh_vs_Rch, "R_{nh}", "P");


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  //hp_Rgamma_vs_Rch->Draw("");
  //hp_Rnh_vs_Rch->Draw("same");
  hp_divide->Draw("same");
  legend->Draw("Psame");
  c1->SaveAs("Rx_vs_Rch.eps");

}
