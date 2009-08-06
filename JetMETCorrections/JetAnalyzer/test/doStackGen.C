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

  TFile* inFile = TFile::Open("RxGenOutputFile.root");



  TH1D* h1_RchGen_vs_eta = (TH1D*)inFile->Get("RchGen_vs_eta");
  TH1D* h1_RgammaGen_vs_eta = (TH1D*)inFile->Get("RgammaGen_vs_eta");
  TH1D* h1_RnhGen_vs_eta = (TH1D*)inFile->Get("RnhGen_vs_eta");
//  TH1D* h1_RmuGen_vs_eta = (TH1D*)inFile->Get("RmuGen_vs_eta_PFItCone5");

  int nBinsX_RchGen = h1_RchGen_vs_eta->GetNbinsX();
  int nBinsX_RgammaGen = h1_RgammaGen_vs_eta->GetNbinsX();
  int nBinsX_RnhGen = h1_RnhGen_vs_eta->GetNbinsX();
//  int nBinsX_RmuGen = h1_RmuGen_vs_eta->GetNbinsX();

  double minX_RchGen = h1_RchGen_vs_eta->GetXaxis()->GetXmin();
  double maxX_RchGen = h1_RchGen_vs_eta->GetXaxis()->GetXmax();

  double minX_RgammaGen = h1_RgammaGen_vs_eta->GetXaxis()->GetXmin();
  double maxX_RgammaGen = h1_RgammaGen_vs_eta->GetXaxis()->GetXmax();

  double minX_RnhGen = h1_RnhGen_vs_eta->GetXaxis()->GetXmin();
  double maxX_RnhGen = h1_RnhGen_vs_eta->GetXaxis()->GetXmax();

//  double minX_RmuGen = h1_RmuGen_vs_eta->GetXaxis()->GetXmin();
//  double maxX_RmuGen = h1_RmuGen_vs_eta->GetXaxis()->GetXmax();

  TH1D* h1_RchGen_vs_eta_stack = new TH1D("RchGen_vs_eta_stack", "", nBinsX_RchGen, minX_RchGen, maxX_RchGen);
  h1_RchGen_vs_eta_stack->SetFillColor(kRed);
  h1_RchGen_vs_eta_stack->SetXTitle("#eta^{GEN}");
 
  TH1D* h1_RgammaGen_vs_eta_stack = new TH1D("RgammaGen_vs_eta_stack", "", nBinsX_RgammaGen, minX_RgammaGen, maxX_RgammaGen);
  h1_RgammaGen_vs_eta_stack->SetFillColor(kGreen);
  h1_RgammaGen_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_RnhGen_vs_eta_stack = new TH1D("RnhGen_vs_eta_stack", "", nBinsX_RnhGen, minX_RnhGen, maxX_RnhGen);
  h1_RnhGen_vs_eta_stack->SetFillColor(kBlue);
  h1_RnhGen_vs_eta_stack->SetXTitle("#eta^{GEN}");

//  TH1D* h1_RmuGen_vs_eta_stack = new TH1D("RmuGen_vs_eta", "", nBinsX_RmuGen, minX_RmuGen, maxX_RmuGen);
//  h1_RmuGen_vs_eta_stack->SetFillColor(kYellow);
 // h1_RmuGen_vs_eta_stack->SetXTitle("#eta^{GEN}");

  for(int i=1; i<101; ++i) {
    h1_RchGen_vs_eta_stack->SetBinContent( i, h1_RchGen_vs_eta->GetBinContent(i) );
    h1_RgammaGen_vs_eta_stack->SetBinContent( i, h1_RgammaGen_vs_eta->GetBinContent(i) );
    h1_RnhGen_vs_eta_stack->SetBinContent( i, h1_RnhGen_vs_eta->GetBinContent(i) );
//    h1_RmuGen_vs_eta_stack->SetBinContent( i, h1_RmuGen_vs_eta->GetBinContent(i) );
  }


  THStack* stack = new THStack("stack", "");
  stack->Add(h1_RchGen_vs_eta_stack);
  stack->Add(h1_RgammaGen_vs_eta_stack);
  stack->Add(h1_RnhGen_vs_eta_stack);
//  stack->Add(h1_RmuGen_vs_eta_stack);
  stack->Draw();
  stack->GetXaxis()->SetTitle("#eta^{GEN}");
  stack->GetYaxis()->SetTitle("% of Jet Energy");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  stack->Draw();
  c1->SaveAs("Plots/genstack_vs_eta.eps");

  TH1D* h1_Rnh_neutron_vs_eta = (TH1D*)inFile->Get("Rnh_neutron_vs_eta");
  TH1D* h1_Rnh_K0L_vs_eta = (TH1D*)inFile->Get("Rnh_K0L_vs_eta");
  TH1D* h1_Rnh_K0S_vs_eta = (TH1D*)inFile->Get("Rnh_K0S_vs_eta");
  TH1D* h1_Rnh_lambda_vs_eta = (TH1D*)inFile->Get("Rnh_lambda_vs_eta");
  TH1D* h1_Rnh_csi_vs_eta = (TH1D*)inFile->Get("Rnh_csi_vs_eta");

  int nBinsX_Rnh_neutron = h1_Rnh_neutron_vs_eta->GetNbinsX();
  int nBinsX_Rnh_K0L = h1_Rnh_K0L_vs_eta->GetNbinsX();
  int nBinsX_Rnh_K0S = h1_Rnh_K0S_vs_eta->GetNbinsX();
  int nBinsX_Rnh_lambda = h1_Rnh_lambda_vs_eta->GetNbinsX();
  int nBinsX_Rnh_csi = h1_Rnh_csi_vs_eta->GetNbinsX();

  double minX_Rnh_neutron = h1_Rnh_neutron_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh_neutron = h1_Rnh_neutron_vs_eta->GetXaxis()->GetXmax();

  double minX_Rnh_K0L = h1_Rnh_K0L_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh_K0L = h1_Rnh_K0L_vs_eta->GetXaxis()->GetXmax();

  double minX_Rnh_K0S = h1_Rnh_K0S_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh_K0S = h1_Rnh_K0S_vs_eta->GetXaxis()->GetXmax();

  double minX_Rnh_lambda = h1_Rnh_lambda_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh_lambda = h1_Rnh_lambda_vs_eta->GetXaxis()->GetXmax();

  double minX_Rnh_csi = h1_Rnh_csi_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh_csi = h1_Rnh_csi_vs_eta->GetXaxis()->GetXmax();


  TH1D* h1_Rnh_neutron_vs_eta_stack = new TH1D("Rnh_neutron_vs_eta_stack", "", nBinsX_Rnh_neutron, minX_Rnh_neutron, maxX_Rnh_neutron);
  h1_Rnh_neutron_vs_eta_stack->SetFillColor(kBlue+4);
  h1_Rnh_neutron_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rnh_K0L_vs_eta_stack = new TH1D("Rnh_K0L_vs_eta_stack", "", nBinsX_Rnh_K0L, minX_Rnh_K0L, maxX_Rnh_K0L);
  h1_Rnh_K0L_vs_eta_stack->SetFillColor(kBlue);
  h1_Rnh_K0L_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rnh_K0S_vs_eta_stack = new TH1D("Rnh_K0S_vs_eta_stack", "", nBinsX_Rnh_K0S, minX_Rnh_K0S, maxX_Rnh_K0S);
  h1_Rnh_K0S_vs_eta_stack->SetFillColor(kCyan);
  h1_Rnh_K0S_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rnh_lambda_vs_eta_stack = new TH1D("Rnh_lambda_vs_eta_stack", "", nBinsX_Rnh_lambda, minX_Rnh_lambda, maxX_Rnh_lambda);
  h1_Rnh_lambda_vs_eta_stack->SetFillColor(kCyan+2);
  h1_Rnh_lambda_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rnh_csi_vs_eta_stack = new TH1D("Rnh_csi_vs_eta_stack", "", nBinsX_Rnh_csi, minX_Rnh_csi, maxX_Rnh_csi);
  h1_Rnh_csi_vs_eta_stack->SetFillColor(kCyan+4);
  h1_Rnh_csi_vs_eta_stack->SetXTitle("#eta^{GEN}");

  for(int i=1; i<101; ++i) {
    h1_Rnh_neutron_vs_eta_stack->SetBinContent( i, h1_Rnh_neutron_vs_eta->GetBinContent(i) );
    h1_Rnh_K0L_vs_eta_stack->SetBinContent( i, h1_Rnh_K0L_vs_eta->GetBinContent(i) );
    h1_Rnh_K0S_vs_eta_stack->SetBinContent( i, h1_Rnh_K0S_vs_eta->GetBinContent(i) );
    h1_Rnh_lambda_vs_eta_stack->SetBinContent( i, h1_Rnh_lambda_vs_eta->GetBinContent(i) );
    h1_Rnh_csi_vs_eta_stack->SetBinContent( i, h1_Rnh_csi_vs_eta->GetBinContent(i) );
  }


  THStack* stackNH = new THStack("stackNH", "");
  stackNH->Add(h1_Rnh_neutron_vs_eta_stack);
  stackNH->Add(h1_Rnh_K0L_vs_eta_stack);
  stackNH->Add(h1_Rnh_K0S_vs_eta_stack);
  stackNH->Add(h1_Rnh_lambda_vs_eta_stack);
  stackNH->Add(h1_Rnh_csi_vs_eta_stack);
  stackNH->Draw();
  stackNH->GetXaxis()->SetTitle("#eta^{GEN}");
  stackNH->GetYaxis()->SetTitle("% of Jet Neutral Hadron Energy");

  TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
  c3->cd();
  stackNH->Draw();
  c3->SaveAs("Plots/genstackNH_vs_eta.eps");

/*
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


  TH1D* h1_RchGen_vs_pt = (TH1D*)inFile->Get("RchGen_vs_pt_PFItCone5");
  TH1D* h1_RgammaGen_vs_pt = (TH1D*)inFile->Get("RgammaGen_vs_pt_PFItCone5");
  TH1D* h1_RnhGen_vs_pt = (TH1D*)inFile->Get("RnhGen_vs_pt_PFItCone5");
//  TH1D* h1_RmuGen_vs_pt = (TH1D*)inFile->Get("RmuGen_vs_pt_PFItCone5");

  TH1D* h1_RchGen_vs_pt_stack = new TH1D("RchGen_vs_pt", "", nBins-2, Lower);
  h1_RchGen_vs_pt_stack->SetFillColor(kRed);
  h1_RchGen_vs_pt_stack->SetXTitle("p_{T}^{GEN}");
 
  TH1D* h1_RgammaGen_vs_pt_stack = new TH1D("RgammaGen_vs_pt", "", nBins-2, Lower);
  h1_RgammaGen_vs_pt_stack->SetFillColor(kGreen);
  h1_RgammaGen_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_RnhGen_vs_pt_stack = new TH1D("RnhGen_vs_pt", "", nBins-2, Lower);
  h1_RnhGen_vs_pt_stack->SetFillColor(kBlue);
  h1_RnhGen_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

//  TH1D* h1_RmuGen_vs_pt_stack = new TH1D("RmuGen_vs_pt", "", nBins-2, Lower);
//  h1_RmuGen_vs_pt_stack->SetFillColor(kYellow);
//  h1_RmuGen_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  for(int i=1; i<21; ++i) {
    h1_RchGen_vs_pt_stack->SetBinContent( i, h1_RchGen_vs_pt->GetBinContent(i) );
    h1_RgammaGen_vs_pt_stack->SetBinContent( i, h1_RgammaGen_vs_pt->GetBinContent(i) );
    h1_RnhGen_vs_pt_stack->SetBinContent( i, h1_RnhGen_vs_pt->GetBinContent(i) );
//    h1_RmuGen_vs_pt_stack->SetBinContent( i, h1_RmuGen_vs_pt->GetBinContent(i) );
  }


  THStack* stack2 = new THStack("stack2", "");
  stack2->Add(h1_RchGen_vs_pt_stack);
  stack2->Add(h1_RgammaGen_vs_pt_stack);
  stack2->Add(h1_RnhGen_vs_pt_stack);
//  stack2->Add(h1_RmuGen_vs_pt_stack);
  stack2->Draw();
  stack2->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  stack2->GetYaxis()->SetTitle("% of Jet Energy");

  TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();
  stack2->Draw();
  c2->SaveAs("Plots/genstack_vs_pt.eps");


  TH1D* h1_RchGen_vs_pt_barrel = (TH1D*)inFile->Get("RchGen_vs_pt_barrel_PFItCone5");
  TH1D* h1_RgammaGen_vs_pt_barrel = (TH1D*)inFile->Get("RgammaGen_vs_pt_barrel_PFItCone5");
  TH1D* h1_RnhGen_vs_pt_barrel = (TH1D*)inFile->Get("RnhGen_vs_pt_barrel_PFItCone5");
//  TH1D* h1_RmuGen_vs_pt_barrel = (TH1D*)inFile->Get("RmuGen_vs_pt_barrel_PFItCone5");

  TH1D* h1_RchGen_vs_pt_barrel_stack = new TH1D("RchGen_vs_pt_barrel", "", nBins-2, Lower);
  h1_RchGen_vs_pt_barrel_stack->SetFillColor(kRed);
  h1_RchGen_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");
 
  TH1D* h1_RgammaGen_vs_pt_barrel_stack = new TH1D("RgammaGen_vs_pt_barrel", "", nBins-2, Lower);
  h1_RgammaGen_vs_pt_barrel_stack->SetFillColor(kGreen);
  h1_RgammaGen_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_RnhGen_vs_pt_barrel_stack = new TH1D("RnhGen_vs_pt_barrel", "", nBins-2, Lower);
  h1_RnhGen_vs_pt_barrel_stack->SetFillColor(kBlue);
  h1_RnhGen_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");

//  TH1D* h1_RmuGen_vs_pt_barrel_stack = new TH1D("RmuGen_vs_pt_barrel", "", nBins-2, Lower);
//  h1_RmuGen_vs_pt_barrel_stack->SetFillColor(kYellow);
//  h1_RmuGen_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");

  for(int i=1; i<21; ++i) {
    h1_RchGen_vs_pt_barrel_stack->SetBinContent( i, h1_RchGen_vs_pt_barrel->GetBinContent(i) );
    h1_RgammaGen_vs_pt_barrel_stack->SetBinContent( i, h1_RgammaGen_vs_pt_barrel->GetBinContent(i) );
    h1_RnhGen_vs_pt_barrel_stack->SetBinContent( i, h1_RnhGen_vs_pt_barrel->GetBinContent(i) );
 //   h1_RmuGen_vs_pt_barrel_stack->SetBinContent( i, h1_RmuGen_vs_pt_barrel->GetBinContent(i) );
  }


  THStack* stack3 = new THStack("stack3", "");
  stack3->Add(h1_RchGen_vs_pt_barrel_stack);
  stack3->Add(h1_RgammaGen_vs_pt_barrel_stack);
  stack3->Add(h1_RnhGen_vs_pt_barrel_stack);
//  stack3->Add(h1_RmuGen_vs_pt_barrel_stack);
  stack3->Draw();
  stack3->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  stack3->GetYaxis()->SetTitle("% of Jet Energy");

  c2->cd();
  stack3->Draw();
  c2->SaveAs("Plots/genstack_vs_pt_barrel.eps");


  TH1D* h1_RchGen_vs_pt_endcap = (TH1D*)inFile->Get("RchGen_vs_pt_endcap_PFItCone5");
  TH1D* h1_RgammaGen_vs_pt_endcap = (TH1D*)inFile->Get("RgammaGen_vs_pt_endcap_PFItCone5");
  TH1D* h1_RnhGen_vs_pt_endcap = (TH1D*)inFile->Get("RnhGen_vs_pt_endcap_PFItCone5");
//  TH1D* h1_RmuGen_vs_pt_endcap = (TH1D*)inFile->Get("RmuGen_vs_pt_endcap_PFItCone5");

  TH1D* h1_RchGen_vs_pt_endcap_stack = new TH1D("RchGen_vs_pt_endcap", "", nBins-2, Lower);
  h1_RchGen_vs_pt_endcap_stack->SetFillColor(kRed);
  h1_RchGen_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_RgammaGen_vs_pt_endcap_stack = new TH1D("RgammaGen_vs_pt_endcap", "", nBins-2, Lower);
  h1_RgammaGen_vs_pt_endcap_stack->SetFillColor(kGreen);
  h1_RgammaGen_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_RnhGen_vs_pt_endcap_stack = new TH1D("RnhGen_vs_pt_endcap", "", nBins-2, Lower);
  h1_RnhGen_vs_pt_endcap_stack->SetFillColor(kBlue);
  h1_RnhGen_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

//  TH1D* h1_RmuGen_vs_pt_endcap_stack = new TH1D("RmuGen_vs_pt_endcap", "", nBins-2, Lower);
//  h1_RmuGen_vs_pt_endcap_stack->SetFillColor(kYellow);
//  h1_RmuGen_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  for(int i=1; i<21; ++i) {
    h1_RchGen_vs_pt_endcap_stack->SetBinContent( i, h1_RchGen_vs_pt_endcap->GetBinContent(i) );
    h1_RgammaGen_vs_pt_endcap_stack->SetBinContent( i, h1_RgammaGen_vs_pt_endcap->GetBinContent(i) );
    h1_RnhGen_vs_pt_endcap_stack->SetBinContent( i, h1_RnhGen_vs_pt_endcap->GetBinContent(i) );
    //h1_RmuGen_vs_pt_endcap_stack->SetBinContent( i, h1_RmuGen_vs_pt_endcap->GetBinContent(i) );
  }


  THStack* stack4 = new THStack("stack4", "");
  stack4->Add(h1_RchGen_vs_pt_endcap_stack);
  stack4->Add(h1_RgammaGen_vs_pt_endcap_stack);
  stack4->Add(h1_RnhGen_vs_pt_endcap_stack);
//  stack4->Add(h1_RmuGen_vs_pt_endcap_stack);
  stack4->Draw();
  stack4->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  stack4->GetYaxis()->SetTitle("% of Jet Energy");

  c2->cd();
  stack4->Draw();
  c2->SaveAs("Plots/genstack_vs_pt_endcap.eps");
*/
} 

