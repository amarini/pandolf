#include "TCanvas.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include <string>



void doStack() {

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

  TFile* inFile = TFile::Open("finalOutputFile.root");



  TH1D* h1_Rch_vs_eta = (TH1D*)inFile->Get("Rch_vs_eta_PFItCone5");
  TH1D* h1_Rgamma_vs_eta = (TH1D*)inFile->Get("Rgamma_vs_eta_PFItCone5");
  TH1D* h1_Rnh_vs_eta = (TH1D*)inFile->Get("Rnh_vs_eta_PFItCone5");
  TH1D* h1_Rmu_vs_eta = (TH1D*)inFile->Get("Rmu_vs_eta_PFItCone5");
  TH1D* h1_Rhf_vs_eta = (TH1D*)inFile->Get("Rhf_vs_eta_PFItCone5");

  int nBinsX_Rch = h1_Rch_vs_eta->GetNbinsX();
  int nBinsX_Rgamma = h1_Rgamma_vs_eta->GetNbinsX();
  int nBinsX_Rnh = h1_Rnh_vs_eta->GetNbinsX();
  int nBinsX_Rmu = h1_Rmu_vs_eta->GetNbinsX();
  int nBinsX_Rhf = h1_Rhf_vs_eta->GetNbinsX();

  double minX_Rch = h1_Rch_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rch = h1_Rch_vs_eta->GetXaxis()->GetXmax();

  double minX_Rgamma = h1_Rgamma_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rgamma = h1_Rgamma_vs_eta->GetXaxis()->GetXmax();

  double minX_Rnh = h1_Rnh_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh = h1_Rnh_vs_eta->GetXaxis()->GetXmax();

  double minX_Rmu = h1_Rmu_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rmu = h1_Rmu_vs_eta->GetXaxis()->GetXmax();

  double minX_Rhf = h1_Rhf_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rhf = h1_Rhf_vs_eta->GetXaxis()->GetXmax();

  TH1D* h1_Rch_vs_eta_stack = new TH1D("Rch_vs_eta", "", nBinsX_Rch, minX_Rch, maxX_Rch);
  h1_Rch_vs_eta_stack->SetFillColor(kRed);
  h1_Rch_vs_eta_stack->SetXTitle("#eta^{GEN}");
 
  TH1D* h1_Rgamma_vs_eta_stack = new TH1D("Rgamma_vs_eta", "", nBinsX_Rgamma, minX_Rgamma, maxX_Rgamma);
  h1_Rgamma_vs_eta_stack->SetFillColor(kGreen);
  h1_Rgamma_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rnh_vs_eta_stack = new TH1D("Rnh_vs_eta", "", nBinsX_Rnh, minX_Rnh, maxX_Rnh);
  h1_Rnh_vs_eta_stack->SetFillColor(kBlue);
  h1_Rnh_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rmu_vs_eta_stack = new TH1D("Rmu_vs_eta", "", nBinsX_Rmu, minX_Rmu, maxX_Rmu);
  h1_Rmu_vs_eta_stack->SetFillColor(kYellow);
  h1_Rmu_vs_eta_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rhf_vs_eta_stack = new TH1D("Rhf_vs_eta", "", nBinsX_Rhf, minX_Rhf, maxX_Rhf);
  h1_Rhf_vs_eta_stack->SetFillColor(kMagenta);
  h1_Rhf_vs_eta_stack->SetXTitle("#eta^{GEN}");

  for(int i=1; i<101; ++i) {
    h1_Rch_vs_eta_stack->SetBinContent( i, h1_Rch_vs_eta->GetBinContent(i) );
    h1_Rgamma_vs_eta_stack->SetBinContent( i, h1_Rgamma_vs_eta->GetBinContent(i) );
    h1_Rnh_vs_eta_stack->SetBinContent( i, h1_Rnh_vs_eta->GetBinContent(i) );
    h1_Rmu_vs_eta_stack->SetBinContent( i, h1_Rmu_vs_eta->GetBinContent(i) );
    h1_Rhf_vs_eta_stack->SetBinContent( i, h1_Rhf_vs_eta->GetBinContent(i) );
  }


  THStack* stack = new THStack("stack", "");
  stack->Add(h1_Rch_vs_eta_stack);
  stack->Add(h1_Rgamma_vs_eta_stack);
  stack->Add(h1_Rnh_vs_eta_stack);
  stack->Add(h1_Rhf_vs_eta_stack);
  stack->Add(h1_Rmu_vs_eta_stack);
  stack->Draw();
  stack->GetXaxis()->SetTitle("#eta^{GEN}");
  stack->GetYaxis()->SetTitle("% of Jet Energy");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  stack->Draw();
  c1->SaveAs("Plots/stack_vs_eta.eps");


  TH1D* h1_Rgamma_vs_eta_Rch0 = (TH1D*)inFile->Get("Rgamma_vs_eta_Rch0PFItCone5");
  TH1D* h1_Rnh_vs_eta_Rch0 = (TH1D*)inFile->Get("Rnh_vs_eta_Rch0PFItCone5");
  TH1D* h1_Rmu_vs_eta_Rch0 = (TH1D*)inFile->Get("Rmu_vs_eta_Rch0PFItCone5");
  TH1D* h1_Rhf_vs_eta_Rch0 = (TH1D*)inFile->Get("Rhf_vs_eta_Rch0PFItCone5");

  TH1D* h1_Rch_vs_eta_Rch0_stack = new TH1D("Rch_vs_eta_Rch0", "", 100, -5.5, 5.5);
  h1_Rch_vs_eta_Rch0_stack->SetFillColor(kRed);
  h1_Rch_vs_eta_Rch0_stack->SetXTitle("#eta^{GEN}");
 
  TH1D* h1_Rgamma_vs_eta_Rch0_stack = new TH1D("Rgamma_vs_eta_Rch0", "", 100, -5.5, 5.5);
  h1_Rgamma_vs_eta_Rch0_stack->SetFillColor(kGreen);
  h1_Rgamma_vs_eta_Rch0_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rnh_vs_eta_Rch0_stack = new TH1D("Rnh_vs_eta_Rch0", "", 100, -5.5, 5.5);
  h1_Rnh_vs_eta_Rch0_stack->SetFillColor(kBlue);
  h1_Rnh_vs_eta_Rch0_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rmu_vs_eta_Rch0_stack = new TH1D("Rmu_vs_eta_Rch0", "", 100, -5.5, 5.5);
  h1_Rmu_vs_eta_Rch0_stack->SetFillColor(kYellow);
  h1_Rmu_vs_eta_Rch0_stack->SetXTitle("#eta^{GEN}");

  TH1D* h1_Rhf_vs_eta_Rch0_stack = new TH1D("Rhf_vs_eta_Rch0", "", 100, -5.5, 5.5);
  h1_Rhf_vs_eta_Rch0_stack->SetFillColor(kMagenta);
  h1_Rhf_vs_eta_Rch0_stack->SetXTitle("#eta^{GEN}");

  for(int i=1; i<101; ++i) {
    h1_Rgamma_vs_eta_Rch0_stack->SetBinContent( i, h1_Rgamma_vs_eta_Rch0->GetBinContent(i) );
    h1_Rnh_vs_eta_Rch0_stack->SetBinContent( i, h1_Rnh_vs_eta_Rch0->GetBinContent(i) );
    h1_Rmu_vs_eta_Rch0_stack->SetBinContent( i, h1_Rmu_vs_eta_Rch0->GetBinContent(i) );
    h1_Rhf_vs_eta_Rch0_stack->SetBinContent( i, h1_Rhf_vs_eta_Rch0->GetBinContent(i) );
  }


  THStack* stack_Rch0 = new THStack("stack_Rch0", "");
  stack_Rch0->Add(h1_Rgamma_vs_eta_Rch0_stack);
  stack_Rch0->Add(h1_Rnh_vs_eta_Rch0_stack);
  stack_Rch0->Add(h1_Rhf_vs_eta_Rch0_stack);
  stack_Rch0->Add(h1_Rmu_vs_eta_Rch0_stack);
  stack_Rch0->Draw();
  stack_Rch0->GetXaxis()->SetTitle("#eta^{GEN}");
  stack_Rch0->GetYaxis()->SetTitle("% of Jet Energy");

  c1->cd();
  stack_Rch0->Draw();
  c1->SaveAs("Plots/stack_vs_eta_Rch0.eps");


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


  TH1D* h1_Rch_vs_pt = (TH1D*)inFile->Get("Rch_vs_pt_PFItCone5");
  TH1D* h1_Rgamma_vs_pt = (TH1D*)inFile->Get("Rgamma_vs_pt_PFItCone5");
  TH1D* h1_Rnh_vs_pt = (TH1D*)inFile->Get("Rnh_vs_pt_PFItCone5");
  TH1D* h1_Rhf_vs_pt = (TH1D*)inFile->Get("Rhf_vs_pt_PFItCone5");
  TH1D* h1_Rmu_vs_pt = (TH1D*)inFile->Get("Rmu_vs_pt_PFItCone5");

  TH1D* h1_Rch_vs_pt_stack = new TH1D("Rch_vs_pt", "", nBins-2, Lower);
  h1_Rch_vs_pt_stack->SetFillColor(kRed);
  h1_Rch_vs_pt_stack->SetXTitle("p_{T}^{GEN}");
 
  TH1D* h1_Rgamma_vs_pt_stack = new TH1D("Rgamma_vs_pt", "", nBins-2, Lower);
  h1_Rgamma_vs_pt_stack->SetFillColor(kGreen);
  h1_Rgamma_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rnh_vs_pt_stack = new TH1D("Rnh_vs_pt", "", nBins-2, Lower);
  h1_Rnh_vs_pt_stack->SetFillColor(kBlue);
  h1_Rnh_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rhf_vs_pt_stack = new TH1D("Rhf_vs_pt", "", nBins-2, Lower);
  h1_Rhf_vs_pt_stack->SetFillColor(kMagenta);
  h1_Rhf_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rmu_vs_pt_stack = new TH1D("Rmu_vs_pt", "", nBins-2, Lower);
  h1_Rmu_vs_pt_stack->SetFillColor(kYellow);
  h1_Rmu_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  for(int i=1; i<21; ++i) {
    h1_Rch_vs_pt_stack->SetBinContent( i, h1_Rch_vs_pt->GetBinContent(i) );
    h1_Rgamma_vs_pt_stack->SetBinContent( i, h1_Rgamma_vs_pt->GetBinContent(i) );
    h1_Rnh_vs_pt_stack->SetBinContent( i, h1_Rnh_vs_pt->GetBinContent(i) );
    h1_Rhf_vs_pt_stack->SetBinContent( i, h1_Rhf_vs_pt->GetBinContent(i) );
    h1_Rmu_vs_pt_stack->SetBinContent( i, h1_Rmu_vs_pt->GetBinContent(i) );
  }


  THStack* stack2 = new THStack("stack2", "");
  stack2->Add(h1_Rch_vs_pt_stack);
  stack2->Add(h1_Rgamma_vs_pt_stack);
  stack2->Add(h1_Rnh_vs_pt_stack);
  stack2->Add(h1_Rhf_vs_pt_stack);
  stack2->Add(h1_Rmu_vs_pt_stack);
  stack2->Draw();
  stack2->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  stack2->GetYaxis()->SetTitle("% of Jet Energy");

  TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();
  stack2->Draw();
  c2->SaveAs("Plots/stack_vs_pt.eps");


  TH1D* h1_Rch_vs_pt_barrel = (TH1D*)inFile->Get("Rch_vs_pt_barrel_PFItCone5");
  TH1D* h1_Rgamma_vs_pt_barrel = (TH1D*)inFile->Get("Rgamma_vs_pt_barrel_PFItCone5");
  TH1D* h1_Rnh_vs_pt_barrel = (TH1D*)inFile->Get("Rnh_vs_pt_barrel_PFItCone5");
  TH1D* h1_Rmu_vs_pt_barrel = (TH1D*)inFile->Get("Rmu_vs_pt_barrel_PFItCone5");

  TH1D* h1_Rch_vs_pt_barrel_stack = new TH1D("Rch_vs_pt_barrel", "", nBins-2, Lower);
  h1_Rch_vs_pt_barrel_stack->SetFillColor(kRed);
  h1_Rch_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");
 
  TH1D* h1_Rgamma_vs_pt_barrel_stack = new TH1D("Rgamma_vs_pt_barrel", "", nBins-2, Lower);
  h1_Rgamma_vs_pt_barrel_stack->SetFillColor(kGreen);
  h1_Rgamma_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rnh_vs_pt_barrel_stack = new TH1D("Rnh_vs_pt_barrel", "", nBins-2, Lower);
  h1_Rnh_vs_pt_barrel_stack->SetFillColor(kBlue);
  h1_Rnh_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rmu_vs_pt_barrel_stack = new TH1D("Rmu_vs_pt_barrel", "", nBins-2, Lower);
  h1_Rmu_vs_pt_barrel_stack->SetFillColor(kYellow);
  h1_Rmu_vs_pt_barrel_stack->SetXTitle("p_{T}^{GEN}");

  for(int i=1; i<21; ++i) {
    h1_Rch_vs_pt_barrel_stack->SetBinContent( i, h1_Rch_vs_pt_barrel->GetBinContent(i) );
    h1_Rgamma_vs_pt_barrel_stack->SetBinContent( i, h1_Rgamma_vs_pt_barrel->GetBinContent(i) );
    h1_Rnh_vs_pt_barrel_stack->SetBinContent( i, h1_Rnh_vs_pt_barrel->GetBinContent(i) );
    h1_Rmu_vs_pt_barrel_stack->SetBinContent( i, h1_Rmu_vs_pt_barrel->GetBinContent(i) );
  }


  THStack* stack3 = new THStack("stack3", "");
  stack3->Add(h1_Rch_vs_pt_barrel_stack);
  stack3->Add(h1_Rgamma_vs_pt_barrel_stack);
  stack3->Add(h1_Rnh_vs_pt_barrel_stack);
  stack3->Add(h1_Rmu_vs_pt_barrel_stack);
  stack3->Draw();
  stack3->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  stack3->GetYaxis()->SetTitle("% of Jet Energy");

  c2->cd();
  stack3->Draw();
  c2->SaveAs("Plots/stack_vs_pt_barrel.eps");


  TH1D* h1_Rch_vs_pt_endcap = (TH1D*)inFile->Get("Rch_vs_pt_endcap_PFItCone5");
  TH1D* h1_Rgamma_vs_pt_endcap = (TH1D*)inFile->Get("Rgamma_vs_pt_endcap_PFItCone5");
  TH1D* h1_Rnh_vs_pt_endcap = (TH1D*)inFile->Get("Rnh_vs_pt_endcap_PFItCone5");
  TH1D* h1_Rhf_vs_pt_endcap = (TH1D*)inFile->Get("Rhf_vs_pt_endcap_PFItCone5");
  TH1D* h1_Rmu_vs_pt_endcap = (TH1D*)inFile->Get("Rmu_vs_pt_endcap_PFItCone5");

  TH1D* h1_Rch_vs_pt_endcap_stack = new TH1D("Rch_vs_pt_endcap", "", nBins-2, Lower);
  h1_Rch_vs_pt_endcap_stack->SetFillColor(kRed);
  h1_Rch_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rgamma_vs_pt_endcap_stack = new TH1D("Rgamma_vs_pt_endcap", "", nBins-2, Lower);
  h1_Rgamma_vs_pt_endcap_stack->SetFillColor(kGreen);
  h1_Rgamma_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rhf_vs_pt_endcap_stack = new TH1D("Rhf_vs_pt_endcap", "", nBins-2, Lower);
  h1_Rhf_vs_pt_endcap_stack->SetFillColor(kMagenta);
  h1_Rhf_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rnh_vs_pt_endcap_stack = new TH1D("Rnh_vs_pt_endcap", "", nBins-2, Lower);
  h1_Rnh_vs_pt_endcap_stack->SetFillColor(kBlue);
  h1_Rnh_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  TH1D* h1_Rmu_vs_pt_endcap_stack = new TH1D("Rmu_vs_pt_endcap", "", nBins-2, Lower);
  h1_Rmu_vs_pt_endcap_stack->SetFillColor(kYellow);
  h1_Rmu_vs_pt_endcap_stack->SetXTitle("p_{T}^{GEN}");

  for(int i=1; i<21; ++i) {
    h1_Rch_vs_pt_endcap_stack->SetBinContent( i, h1_Rch_vs_pt_endcap->GetBinContent(i) );
    h1_Rgamma_vs_pt_endcap_stack->SetBinContent( i, h1_Rgamma_vs_pt_endcap->GetBinContent(i) );
    h1_Rnh_vs_pt_endcap_stack->SetBinContent( i, h1_Rnh_vs_pt_endcap->GetBinContent(i) );
    h1_Rhf_vs_pt_endcap_stack->SetBinContent( i, h1_Rhf_vs_pt_endcap->GetBinContent(i) );
    h1_Rmu_vs_pt_endcap_stack->SetBinContent( i, h1_Rmu_vs_pt_endcap->GetBinContent(i) );
  }


  THStack* stack4 = new THStack("stack4", "");
  stack4->Add(h1_Rch_vs_pt_endcap_stack);
  stack4->Add(h1_Rgamma_vs_pt_endcap_stack);
  stack4->Add(h1_Rnh_vs_pt_endcap_stack);
  stack4->Add(h1_Rhf_vs_pt_endcap_stack);
  stack4->Add(h1_Rmu_vs_pt_endcap_stack);
  stack4->Draw();
  stack4->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  stack4->GetYaxis()->SetTitle("% of Jet Energy");

  c2->cd();
  stack4->Draw();
  c2->SaveAs("Plots/stack_vs_pt_endcap.eps");

} 

