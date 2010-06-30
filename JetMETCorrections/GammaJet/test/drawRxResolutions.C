#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <string>
#include <iostream>


std::string PICTURE_FORMAT = "eps";


void drawResolution(std::string dataset,std::string varName, bool logScale=false);
void drawResolution_ptBin(std::string dataset,std::string varName, Int_t minPt, Int_t maxPt, bool logScale);

void drawRxResolutions(const std::string& dataset) {

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


  bool logScale=true;
  drawResolution(dataset, "Ech");
  drawResolution(dataset, "Egamma");
  drawResolution(dataset, "Egammanh");
  //drawResolution(dataset, "Ee");
  drawResolution(dataset, "Enh", logScale);

}

void drawResolution( std::string dataset,std::string varName, bool logScale) {

  drawResolution_ptBin( dataset, varName, 15, 20, logScale);
  drawResolution_ptBin( dataset, varName, 62, 70, logScale);
  drawResolution_ptBin( dataset, varName, 121, 135, logScale);
  drawResolution_ptBin( dataset, varName, 397, 600, logScale);

}

void drawResolution_ptBin(std::string dataset,std::string varName, Int_t minPt, Int_t maxPt, bool logScale) {

  std::string fileName = "finalOutputFile_"+dataset+".root";
  TFile* file = TFile::Open(fileName.c_str());

  std::string xTitle;

  char histName[150];
  sprintf(histName, "%s_resolution_barrel_%d_%d_PFakt5", varName.c_str(), minPt, maxPt);
  TH1F* h1= (TH1F*)file->Get(histName);
  h1->SetTitle("");
  h1->GetXaxis()->SetTitleOffset(1.1);
  //h1->GetXaxis()->SetTitleSize(1.5);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitle("Normalized Distribution");
  h1->SetLineWidth(1.5);
  if( varName=="Ech") {
    h1->SetFillColor(kRed);
    xTitle = "(E_{ch}^{RECO} - E_{ch}^{GEN}) / E_{ch}^{GEN}";
  }
  if( varName=="Egamma") {
    h1->SetFillColor(kGreen);
    xTitle = "(E_{#gamma}^{RECO} - E_{#gamma}^{GEN}) / E_{#gamma}^{GEN}";
  }
  if( varName=="Enh") {
    h1->SetFillColor(kBlue);
    xTitle = "(E_{nh}^{RECO} - E_{nh}^{GEN}) / E_{nh}^{GEN}";
  }
  if( varName=="Egammanh") {
    h1->SetFillColor(kCyan);
    xTitle = "(E_{#gamma+nh}^{RECO} - E_{#gamma+nh}^{GEN}) / E_{#gamma+nh}^{GEN}";
  }

  h1->GetXaxis()->SetTitle(xTitle.c_str());

  Float_t textSize = 0.04;

  char meanLabelText[100];
  sprintf(meanLabelText, "mean = %.3f", h1->GetMean());
  TPaveText* meanLabel = new TPaveText(0.15, 0.8, 0.3, 0.85, "brNDC");
  meanLabel->SetFillColor(kWhite);
  meanLabel->SetTextSize(textSize);
  meanLabel->AddText(meanLabelText);

  char rmsLabelText[100];
  sprintf(rmsLabelText, "RMS = %.3f", h1->GetRMS());
  TPaveText* rmsLabel = new TPaveText(0.15, 0.7, 0.3, 0.78, "brNDC");
  rmsLabel->SetFillColor(kWhite);
  rmsLabel->SetTextSize(textSize);
  rmsLabel->AddText(rmsLabelText);

  char ptLabelText[100];
  sprintf(ptLabelText, "%d < p_{T}^{GEN} < %d GeV/c", minPt, maxPt);
  TPaveText* ptLabel = new TPaveText(0.67, 0.8, 0.8, 0.88, "brNDC");
  ptLabel->SetFillColor(kWhite);
  ptLabel->SetTextSize(textSize);
  ptLabel->AddText(ptLabelText);

  TPaveText* etaLabel = new TPaveText(0.67, 0.7, 0.8, 0.78, "brNDC");
  etaLabel->SetFillColor(kWhite);
  etaLabel->SetTextSize(textSize);
  etaLabel->AddText("|#eta_{GEN}| < 1.4");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetBottomMargin(0.105);
  if(logScale)
    c1->SetLogy();
  h1->DrawNormalized();
  meanLabel->Draw("same");
  rmsLabel->Draw("same");
  ptLabel->Draw("same");
  etaLabel->Draw("same");
  std::string plotName(histName);
  plotName = "Plots_"+dataset+"/"+plotName+"."+PICTURE_FORMAT;
  c1->SaveAs(plotName.c_str());

  delete c1;

}
