#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <string>



void drawRx_vs_Rch() {

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


  //TFile* file = TFile::Open("finalOutputFile_226_SAVE.root");
  TFile* file = TFile::Open("provaOutputFile.root");

  TH1F* h1_Rgamma = (TH1F*)file->Get("Rgamma_vs_Rch_barrel_PFItCone5");
  h1_Rgamma->SetTitle("");
  h1_Rgamma->GetXaxis()->SetTitleOffset(1.1);
  h1_Rgamma->GetXaxis()->SetTitleSize(1.2);
  h1_Rgamma->GetXaxis()->SetTitle("R_{ch}");
  h1_Rgamma->GetYaxis()->SetTitleOffset(1.2);
  h1_Rgamma->GetYaxis()->SetTitle("");
  h1_Rgamma->SetMarkerSize(1.2);
  h1_Rgamma->SetMarkerStyle(20);
  h1_Rgamma->SetMarkerColor(kGreen);

  TH1F* h1_Rnh = (TH1F*)file->Get("Rnh_vs_Rch_barrel_PFItCone5");
  h1_Rnh->SetTitle("");
  h1_Rnh->GetXaxis()->SetTitleSize(1.5);
  h1_Rnh->GetXaxis()->SetTitleOffset(1.1);
  h1_Rnh->GetYaxis()->SetTitleOffset(1.2);
  h1_Rnh->GetYaxis()->SetTitle("");
  h1_Rnh->SetMarkerSize(1.2);
  h1_Rnh->SetMarkerStyle(21);
  h1_Rnh->SetMarkerColor(kBlue);


  TLegend* l1 = new TLegend(0.7, 0.7, 0.85, 0.8);
  l1->AddEntry(h1_Rgamma, "R_{#gamma}", "P");
  l1->SetFillColor(kWhite);

  TLegend* l2 = new TLegend(0.2, 0.2, 0.35, 0.3);
  l2->AddEntry(h1_Rnh, "R_{nh}", "P");
  l2->SetFillColor(kWhite);


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h1_Rgamma->Draw("P");
  //h1_Rnh->Draw("Psame");
  l1->Draw("same");
  //l2->Draw("same");

  c1->SaveAs("Rx_vs_Rch.eps");


  TH1F* h1_RgammaGen = (TH1F*)file->Get("RgammaGen_vs_RchGen_PFItCone5");
  h1_RgammaGen->SetTitle("");
  h1_RgammaGen->GetXaxis()->SetTitleOffset(1.1);
  h1_RgammaGen->GetXaxis()->SetTitleSize(1.5);
  h1_RgammaGen->GetYaxis()->SetTitleOffset(1.2);
  h1_RgammaGen->GetYaxis()->SetTitle("");
  h1_RgammaGen->SetMarkerSize(1.2);
  h1_RgammaGen->SetMarkerStyle(20);
  h1_RgammaGen->SetMarkerColor(kGreen);

  TH1F* h1_RnhGen = (TH1F*)file->Get("RnhGen_vs_RchGen_PFItCone5");
  h1_RnhGen->SetTitle("");
  h1_RnhGen->GetXaxis()->SetTitleSize(1.5);
  h1_RnhGen->GetXaxis()->SetTitleOffset(1.1);
  h1_RnhGen->GetYaxis()->SetTitleOffset(1.2);
  h1_RnhGen->GetYaxis()->SetTitle("");
  h1_RnhGen->SetMarkerSize(1.2);
  h1_RnhGen->SetMarkerStyle(21);
  h1_RnhGen->SetMarkerColor(kBlue);


  TLegend* l1 = new TLegend(0.7, 0.7, 0.85, 0.8);
  l1->AddEntry(h1_RgammaGen, "R_{#gamma}", "P");
  l1->SetFillColor(kWhite);

  TLegend* l2 = new TLegend(0.2, 0.2, 0.35, 0.3);
  l2->AddEntry(h1_RnhGen, "R_{nh}", "P");
  l2->SetFillColor(kWhite);


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h1_RgammaGen->Draw("P");
  //h1_RnhGen->Draw("Psame");
  l1->Draw("same");
  //l2->Draw("same");

  c1->SaveAs("Rx_vs_RchGen.eps");

}
  
