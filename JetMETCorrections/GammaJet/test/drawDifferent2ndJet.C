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



void different2ndJet(int iPtBin) {

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


