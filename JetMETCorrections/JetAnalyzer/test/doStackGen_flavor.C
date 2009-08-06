#include "TCanvas.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include <string>



void doStackGen() {

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

  TFile* inFile = TFile::Open("flavorOutputFile.root");



  TH1D* h1_RchGen_vs_eta_gluon = (TH1D*)inFile->Get("RchGen_vs_eta_Rch0_gluon_PFItCone5");
  TH1D* h1_RgammaGen_vs_eta_gluon = (TH1D*)inFile->Get("RgammaGen_vs_eta_Rch0_gluon_PFItCone5");
  TH1D* h1_RnhGen_vs_eta_gluon = (TH1D*)inFile->Get("RnhGen_vs_eta_Rch0_gluon_PFItCone5");
//  TH1D* h1_RmuGen_vs_eta = (TH1D*)inFile->Get("RmuGen_vs_eta_PFItCone5");

  int nBinsX_RchGen = h1_RchGen_vs_eta_gluon->GetNbinsX();
  int nBinsX_RgammaGen = h1_RgammaGen_vs_eta_gluon->GetNbinsX();
  int nBinsX_RnhGen = h1_RnhGen_vs_eta_gluon->GetNbinsX();
//  int nBinsX_RmuGen = h1_RmuGen_vs_eta->GetNbinsX();

  double minX_RchGen = h1_RchGen_vs_eta_gluon->GetXaxis()->GetXmin();
  double maxX_RchGen = h1_RchGen_vs_eta_gluon->GetXaxis()->GetXmax();

  double minX_RgammaGen = h1_RgammaGen_vs_eta_gluon->GetXaxis()->GetXmin();
  double maxX_RgammaGen = h1_RgammaGen_vs_eta_gluon->GetXaxis()->GetXmax();

  double minX_RnhGen = h1_RnhGen_vs_eta_gluon->GetXaxis()->GetXmin();
  double maxX_RnhGen = h1_RnhGen_vs_eta_gluon->GetXaxis()->GetXmax();

//  double minX_RmuGen = h1_RmuGen_vs_eta->GetXaxis()->GetXmin();
//  double maxX_RmuGen = h1_RmuGen_vs_eta->GetXaxis()->GetXmax();

  //TH1D* h1_RchGen_vs_eta_gluon_stack = new TH1D("RchGen_vs_eta_gluon_stack", "", nBinsX_RchGen, minX_RchGen, maxX_RchGen);
  TH1D* h1_RchGen_vs_eta_gluon_stack = new TH1D("RchGen_vs_eta_gluon_stack", "", 100, -5.5, 5.5);
  h1_RchGen_vs_eta_gluon_stack->SetFillColor(kRed);
  h1_RchGen_vs_eta_gluon_stack->SetXTitle("#eta^{GEN}");
 
  //TH1D* h1_RgammaGen_vs_eta_gluon_stack = new TH1D("RgammaGen_vs_eta_gluon_stack", "", nBinsX_RgammaGen, minX_RgammaGen, maxX_RgammaGen);
  TH1D* h1_RgammaGen_vs_eta_gluon_stack = new TH1D("RgammaGen_vs_eta_gluon_stack", "", 100, -5.5, 5.5);
  h1_RgammaGen_vs_eta_gluon_stack->SetFillColor(kGreen);
  h1_RgammaGen_vs_eta_gluon_stack->SetXTitle("#eta^{GEN}");

  //TH1D* h1_RnhGen_vs_eta_gluon_stack = new TH1D("RnhGen_vs_eta_gluon_stack", "", nBinsX_RnhGen, minX_RnhGen, maxX_RnhGen);
  TH1D* h1_RnhGen_vs_eta_gluon_stack = new TH1D("RnhGen_vs_eta_gluon_stack", "", 100, -5.5, 5.5);
  h1_RnhGen_vs_eta_gluon_stack->SetFillColor(kBlue);
  h1_RnhGen_vs_eta_gluon_stack->SetXTitle("#eta^{GEN}");

//  TH1D* h1_RmuGen_vs_eta_stack = new TH1D("RmuGen_vs_eta", "", nBinsX_RmuGen, minX_RmuGen, maxX_RmuGen);
//  h1_RmuGen_vs_eta_stack->SetFillColor(kYellow);
 // h1_RmuGen_vs_eta_stack->SetXTitle("#eta^{GEN}");


  TH1D* h1_RchGen_vs_eta_light = (TH1D*)inFile->Get("RchGen_vs_eta_Rch0_light_PFItCone5");
  TH1D* h1_RgammaGen_vs_eta_light = (TH1D*)inFile->Get("RgammaGen_vs_eta_Rch0_light_PFItCone5");
  TH1D* h1_RnhGen_vs_eta_light = (TH1D*)inFile->Get("RnhGen_vs_eta_Rch0_light_PFItCone5");
//  TH1D* h1_RmuGen_vs_eta = (TH1D*)inFile->Get("RmuGen_vs_eta_PFItCone5");

  nBinsX_RchGen = h1_RchGen_vs_eta_light->GetNbinsX();
  nBinsX_RgammaGen = h1_RgammaGen_vs_eta_light->GetNbinsX();
  nBinsX_RnhGen = h1_RnhGen_vs_eta_light->GetNbinsX();
//  int nBinsX_RmuGen = h1_RmuGen_vs_eta->GetNbinsX();

  minX_RchGen = h1_RchGen_vs_eta_light->GetXaxis()->GetXmin();
  maxX_RchGen = h1_RchGen_vs_eta_light->GetXaxis()->GetXmax();

  minX_RgammaGen = h1_RgammaGen_vs_eta_light->GetXaxis()->GetXmin();
  maxX_RgammaGen = h1_RgammaGen_vs_eta_light->GetXaxis()->GetXmax();

  minX_RnhGen = h1_RnhGen_vs_eta_light->GetXaxis()->GetXmin();
  maxX_RnhGen = h1_RnhGen_vs_eta_light->GetXaxis()->GetXmax();

//  double minX_RmuGen = h1_RmuGen_vs_eta->GetXaxis()->GetXmin();
//  double maxX_RmuGen = h1_RmuGen_vs_eta->GetXaxis()->GetXmax();

  TH1D* h1_RchGen_vs_eta_light_stack = new TH1D("RchGen_vs_eta_light_stack", "", 100, -5.5, 5.5);
  h1_RchGen_vs_eta_light_stack->SetFillColor(kRed);
  h1_RchGen_vs_eta_light_stack->SetXTitle("#eta^{GEN}");
 
  TH1D* h1_RgammaGen_vs_eta_light_stack = new TH1D("RgammaGen_vs_eta_light_stack", "", 100, -5.5, 5.5);
  h1_RgammaGen_vs_eta_light_stack->SetFillColor(kGreen);
  h1_RgammaGen_vs_eta_light_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_RnhGen_vs_eta_light_stack = new TH1D("RnhGen_vs_eta_light_stack", "", 100, -5.5, 5.5);
  h1_RnhGen_vs_eta_light_stack->SetFillColor(kBlue);
  h1_RnhGen_vs_eta_light_stack->SetXTitle("#eta^{GEN}");


  for(int i=1; i<101; ++i) {
    h1_RchGen_vs_eta_gluon_stack->SetBinContent( i, h1_RchGen_vs_eta_gluon->GetBinContent(i) );
    h1_RgammaGen_vs_eta_gluon_stack->SetBinContent( i, h1_RgammaGen_vs_eta_gluon->GetBinContent(i) );
    h1_RnhGen_vs_eta_gluon_stack->SetBinContent( i, h1_RnhGen_vs_eta_gluon->GetBinContent(i) );

    h1_RchGen_vs_eta_light_stack->SetBinContent( i, h1_RchGen_vs_eta_light->GetBinContent(i) );
    h1_RgammaGen_vs_eta_light_stack->SetBinContent( i, h1_RgammaGen_vs_eta_light->GetBinContent(i) );
    h1_RnhGen_vs_eta_light_stack->SetBinContent( i, h1_RnhGen_vs_eta_light->GetBinContent(i) );
  }


  THStack* stack_gluon = new THStack("stack_gluon", "");
  stack_gluon->Add(h1_RchGen_vs_eta_gluon_stack);
  stack_gluon->Add(h1_RgammaGen_vs_eta_gluon_stack);
  stack_gluon->Add(h1_RnhGen_vs_eta_gluon_stack);
//  stack->Add(h1_RmuGen_vs_eta_stack);
  stack_gluon->Draw();
  stack_gluon->GetXaxis()->SetTitle("#eta^{GEN}");
  stack_gluon->GetYaxis()->SetTitle("% of Jet Energy");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  stack_gluon->Draw();
  c1->SaveAs("Plots/genstack_vs_eta_Rch0_gluon.eps");


  THStack* stack_light = new THStack("stack_light", "");
  stack_light->Add(h1_RchGen_vs_eta_light_stack);
  stack_light->Add(h1_RgammaGen_vs_eta_light_stack);
  stack_light->Add(h1_RnhGen_vs_eta_light_stack);
//  stack->Add(h1_RmuGen_vs_eta_stack);
  stack_light->Draw();
  stack_light->GetXaxis()->SetTitle("#eta^{GEN}");
  stack_light->GetYaxis()->SetTitle("% of Jet Energy");

  c1->cd();
  stack_light->Draw();
  c1->SaveAs("Plots/genstack_vs_eta_Rch0_light.eps");



} 

