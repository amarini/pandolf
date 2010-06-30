#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.C"


class Result_Error {
 public:
  Double_t result;
  Double_t error;
};



Result_Error getSinglePtBin(int iPtBin);



void drawCompare() {

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

  int nBins = 23;

  Double_t Lower[nBins];
  Lower[0] = 20.;
  Lower[1] = 25.;
  Lower[2] = 30.;
  Lower[3] = 36.;
  Lower[4] = 42.;
  Lower[5] = 48.;
  Lower[6] = 55.;
  Lower[7] = 62.;
  Lower[8] = 70.;
  Lower[9] = 78.;
  Lower[10] = 87.;
  Lower[11] = 97.;
  Lower[12] = 108.;
  Lower[13] = 121.;
  Lower[14] = 135.;
  Lower[15] = 152.;
  Lower[16] = 172.;
  Lower[17] = 197.;
  Lower[18] = 231.;
  Lower[19] = 283.;
  Lower[20] = 397.;
  Lower[21] = 600.;
  Lower[22] = 1000.;


  Double_t res[nBins];
  Double_t res_err[nBins];

  res[0] = -1.;
  res_err[0] = 0.;

  for(int i=1; i<(nBins-1); ++i) {
    Result_Error resolution = getSinglePtBin(i+1);
    res[i] = resolution.result;
    res_err[i] = resolution.error;
  }

  TFile* file = TFile::Open("RchCalibration_PhotonJet.root");

  TGraphErrors* gr_response_GEN = (TGraphErrors*)file->Get("gr_response_vs_pt_FIT_barrel_PFItCone5_GEN_PRE");
  gr_response_GEN->SetMarkerStyle(29);
  gr_response_GEN->SetMarkerSize(1.5);
  gr_response_GEN->SetMarkerColor(kRed);

  TGraphErrors* gr_response_PHOT = (TGraphErrors*)file->Get("gr_response_vs_pt_FIT_barrel_PFItCone5_PRE");
  gr_response_PHOT->SetMarkerStyle(27);
  gr_response_PHOT->SetMarkerSize(2.);
  gr_response_PHOT->SetMarkerColor(kBlue);
  gr_response_PHOT->SetMarkerSize(1.3);

  TGraphErrors* gr_resolution_GEN = (TGraphErrors*)file->Get("gr_resolution_vs_pt_FIT_barrel_PFItCone5_GEN_PRE");
  gr_resolution_GEN->SetMarkerStyle(29);
  gr_resolution_GEN->SetMarkerSize(1.5);
  gr_resolution_GEN->SetMarkerColor(kRed);

  TGraphErrors* gr_resolution_PHOT = (TGraphErrors*)file->Get("gr_resolution_vs_pt_FIT_barrel_PFItCone5_PRE");
  gr_resolution_PHOT->SetMarkerStyle(27);
  gr_resolution_PHOT->SetMarkerSize(2.);
  gr_resolution_PHOT->SetMarkerColor(kBlue);
  gr_resolution_PHOT->SetMarkerSize(1.3);


  Double_t *x;
  x = gr_resolution_PHOT->GetX();
  Double_t *x_err;
  x_err = gr_resolution_PHOT->GetEX();

  for( int i=0; i<nBins; ++i)
    std::cout << "x[" << i << "]: " << x[i] << "res[" << i << "]: " << res[i] << std::endl;
 
  TGraphErrors* gr_resolution_2ndJet = new TGraphErrors(nBins, x, res, x_err, res_err);
  gr_resolution_2ndJet->SetMarkerStyle(20);    
  gr_resolution_2ndJet->SetMarkerColor(kRed);    

  TH2D* h2_aux = new TH2D("aux", "", 10, 10., 600., 10, 0., 0.5);
  h2_aux->SetXTitle("p_{T}^{#gamma/gen} [GeV/c]");
  h2_aux->SetYTitle("resolution [%]");
  h2_aux->GetYaxis()->SetTitleOffset(1.2);
  h2_aux->GetXaxis()->SetTitleOffset(1.1);

  TLegend* legend = new TLegend(0.6, 0.6, 0.8, 0.8);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(gr_resolution_PHOT, "#gamma + jet", "P");
  legend->AddEntry(gr_resolution_GEN, "MC", "P");
  //legend->AddEntry(gr_resolution_2ndJet, "p_{T}^{2ndJet} #rightarrow 0", "P");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h2_aux->Draw();
//  gr_resolution_2ndJet->Draw("P");
  gr_resolution_GEN->Draw("Psame");
  gr_resolution_PHOT->Draw("Psame");
  legend->Draw("same");

  c1->SaveAs("resolution_compare.eps");

  TLegend* legend2 = new TLegend(0.6, 0.2, 0.8, 0.4);
  legend2->SetFillColor(kWhite);
  legend2->SetFillStyle(0);
  legend2->SetTextSize(0.05);
  legend2->AddEntry(gr_response_PHOT, "#gamma + jet", "P");
  legend2->AddEntry(gr_response_GEN, "MC", "P");

  TH2D* h2_aux2 = new TH2D("aux2", "", 10, 10., 600., 10, 0.7, 1.);
  h2_aux2->SetXTitle("p_{T}^{#gamma/gen} [GeV/c]");
  h2_aux2->SetYTitle("response");
  h2_aux2->GetYaxis()->SetTitleOffset(1.2);
  h2_aux2->GetXaxis()->SetTitleOffset(1.1);
  
  h2_aux2->Draw();
  gr_response_GEN->Draw("Psame");
  gr_response_PHOT->Draw("Psame");
  legend2->Draw("same");

  c1->SaveAs("response_compare.eps");

} //drawPlots




Result_Error getSinglePtBin(int iPtBin) {

  TFile* f_5 =  TFile::Open("Projections_pt2ndJet_5/barrel_PRE.root");
  TFile* f_7 =  TFile::Open("Projections_pt2ndJet_7/barrel_PRE.root");
  TFile* f_10 = TFile::Open("Projections_pt2ndJet_10/barrel_PRE.root");
  TFile* f_12 = TFile::Open("Projections_pt2ndJet_12/barrel_PRE.root");
  TFile* f_15 = TFile::Open("Projections_pt2ndJet_15/barrel_PRE.root");
  //TFile* f_20 = TFile::Open("Projections_pt2ndJet_20/barrel_PRE.root");

  char histName[40];
  sprintf(histName, "projection_%d", iPtBin);

  TH1F* h1_5 = (TH1F*)f_5->Get(histName);
  TH1F* h1_7 = (TH1F*)f_7->Get(histName);
  TH1F* h1_10 = (TH1F*)f_10->Get(histName);
  TH1F* h1_12 = (TH1F*)f_12->Get(histName);
  TH1F* h1_15 = (TH1F*)f_15->Get(histName);
  //TH1F* h1_20 = (TH1F*)f_20->Get(histName);

  Int_t nPoints = 5;

  Float_t x[nPoints];
  Float_t x_err[nPoints];
  Float_t y[nPoints];
  Float_t y_err[nPoints];
  
  x[0] = 5.;
  x[1] = 7.;
  x[2] = 10.;
  x[3] = 12.;
  x[4] = 15.;
  
  x_err[0] = 0.;
  x_err[1] = 0.;
  x_err[2] = 0.;
  x_err[3] = 0.;
  x_err[4] = 0.;
  
  Float_t mu_5 =  h1_5->GetFunction("gaussian")->GetParameter(1);
  Float_t mu_7 =  h1_7->GetFunction("gaussian")->GetParameter(1);
  Float_t mu_10 = h1_10->GetFunction("gaussian")->GetParameter(1);
  Float_t mu_12 = h1_12->GetFunction("gaussian")->GetParameter(1);
  Float_t mu_15 = h1_15->GetFunction("gaussian")->GetParameter(1);

  Float_t mu_err_5 =  h1_5->GetFunction("gaussian")->GetParError(1);
  Float_t mu_err_7 =  h1_7->GetFunction("gaussian")->GetParError(1);
  Float_t mu_err_10 = h1_10->GetFunction("gaussian")->GetParError(1);
  Float_t mu_err_12 = h1_12->GetFunction("gaussian")->GetParError(1);
  Float_t mu_err_15 = h1_15->GetFunction("gaussian")->GetParError(1);

  Float_t sigma_5 =  h1_5->GetFunction("gaussian")->GetParameter(2);
  Float_t sigma_7 =  h1_7->GetFunction("gaussian")->GetParameter(2);
  Float_t sigma_10 = h1_10->GetFunction("gaussian")->GetParameter(2);
  Float_t sigma_12 = h1_12->GetFunction("gaussian")->GetParameter(2);
  Float_t sigma_15 = h1_15->GetFunction("gaussian")->GetParameter(2);

  Float_t sigma_err_5 =  h1_5->GetFunction("gaussian")->GetParError(2);
  Float_t sigma_err_7 =  h1_7->GetFunction("gaussian")->GetParError(2);
  Float_t sigma_err_10 = h1_10->GetFunction("gaussian")->GetParError(2);
  Float_t sigma_err_12 = h1_12->GetFunction("gaussian")->GetParError(2);
  Float_t sigma_err_15 = h1_15->GetFunction("gaussian")->GetParError(2);

  y[0] = sigma_5/mu_5;
  y[1] = sigma_7/mu_7;
  y[2] = sigma_10/mu_10;
  y[3] = sigma_12/mu_12;
  y[4] = sigma_15/mu_15;

  y_err[0] = sqrt( sigma_err_5*sigma_err_5/(mu_5*mu_5) + sigma_5*sigma_5*mu_err_5*mu_err_5/(mu_5*mu_5*mu_5*mu_5) );
  y_err[1] = sqrt( sigma_err_7*sigma_err_7/(mu_7*mu_7) + sigma_7*sigma_7*mu_err_7*mu_err_7/(mu_7*mu_7*mu_7*mu_7) );
  y_err[2] = sqrt( sigma_err_10*sigma_err_10/(mu_10*mu_10) + sigma_10*sigma_10*mu_err_10*mu_err_10/(mu_10*mu_10*mu_10*mu_10) );
  y_err[3] = sqrt( sigma_err_12*sigma_err_12/(mu_12*mu_12) + sigma_12*sigma_12*mu_err_12*mu_err_12/(mu_12*mu_12*mu_12*mu_12) );
  y_err[4] = sqrt( sigma_err_15*sigma_err_15/(mu_15*mu_15) + sigma_15*sigma_15*mu_err_15*mu_err_15/(mu_15*mu_15*mu_15*mu_15) );

  TGraphErrors* gr = new TGraphErrors(nPoints-1, x, y, x_err, y_err);

  Float_t deltaX = x[1]-x[0];
  Float_t deltaY = y[1]-y[0];
  Float_t m = deltaY/deltaX;
  TF1* line = new TF1("line", "[0] + [1]*x");
  line->SetRange(0., 12.);
  line->SetParameter(0, 0.1);
  line->SetParameter(1, m);
  gr->Fit(line, "RQ");

  TH2F* h2_aux = new TH2F("aux", "", 10, 0., 15., 10, 0., 0.5);
  
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h2_aux->Draw();
  gr->Draw("*same");
  char canvasName[60];
  sprintf(canvasName, "Fits_2ndJet/fit%d.eps", iPtBin);
  c1->SaveAs(canvasName);

  Result_Error r;
  r.result = line->GetParameter(0);
  r.error = line->GetParError(0);

  f_5->Close();
  f_7->Close();
  f_10->Close();
  f_12->Close();
  f_15->Close();

  delete c1;

  return r;

}


