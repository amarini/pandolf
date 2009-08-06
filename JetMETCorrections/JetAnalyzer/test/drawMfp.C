#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"





void drawMfp() {

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

  TH1F* total_K0S = new TH1F("total_K0S", "", 100, 0., 250.);
  TH1F* total_lambda = new TH1F("total_lambda", "", 100, 0., 250.);


    
  TFile* file = TFile::Open("QCD_2_2_6_NH/output_3.root");
  TH1F* h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  TH1F* h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_4.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_6.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_9.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_10.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_13.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_17.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_19.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_20.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);

  file = TFile::Open("QCD_2_2_6_NH/output_21.root");
  h1_K0S = (TH1F*)file->Get("myanalysis/mfpK0S_itCone5");
  h1_lambda = (TH1F*)file->Get("myanalysis/mfpLambda_itCone5");
  total_K0S->Add(h1_K0S);
  total_lambda->Add(h1_lambda);


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  TLegend* legend = new TLegend(0.2, 0.68, 0.3, 0.81);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.05);

  legend->AddEntry(total_K0S, "K_{S}^{0}", "P");
  legend->AddEntry(total_lambda, "#Lambda", "P");

  total_K0S->SetMarkerStyle(21);
  total_K0S->SetMarkerColor(kBlue);
  total_K0S->SetMarkerSize(1.2);
  total_K0S->SetXTitle("#lambda_{T} [cm]");
  total_K0S->GetXaxis()->SetRangeUser(0., 170.);
  total_K0S->DrawNormalized("P");

  total_lambda->SetMarkerStyle(24);
  total_lambda->SetMarkerColor(kRed);
  total_lambda->SetMarkerSize(1.2);
  total_lambda->SetXTitle("#lambda_{T} [cm]");
  total_lambda->GetXaxis()->SetRangeUser(0., 170.);
  total_lambda->DrawNormalized("Psame");

  legend->Draw("same");

  c1->SaveAs("mfp.eps");

}
