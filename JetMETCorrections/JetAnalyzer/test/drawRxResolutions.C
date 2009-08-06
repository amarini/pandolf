#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <string>
#include <iostream>


void drawResolution(TFile* file,std::string varName, bool logScale=false);

void drawRxResolutions() {

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


  TFile* file = TFile::Open("finalOutputFile.root");

  bool logScale=true;
  drawResolution(file, "Ech");
  drawResolution(file, "Egamma");
  drawResolution(file, "Egammanh");
  drawResolution(file, "Enh", logScale);

}


void drawResolution(TFile* file,std::string varName, bool logScale) {

  std::string histName;

  histName = varName+"_resolution_PFItCone5";
  TH1F* h1= (TH1F*)file->Get(histName.c_str());
  h1->SetTitle("");
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleSize(1.5);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitle("a.u.");
  h1->SetLineWidth(1.5);
  if( varName=="Ech")
    h1->SetFillColor(kRed);
  if( varName=="Egamma")
    h1->SetFillColor(kGreen);
  if( varName=="Enh")
    h1->SetFillColor(kBlue);
  if( varName=="Egammanh")
    h1->SetFillColor(kCyan);

  std::cout << varName << "\t RMS: " << h1->GetRMS() << std::endl;
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  if(logScale)
    c1->SetLogy();
  h1->DrawNormalized();
  std::string plotName = "Plots/"+histName+".eps";
  c1->SaveAs(plotName.c_str());

  histName = varName+"_resolution_50_80_PFItCone5";
  h1= (TH1F*)file->Get(histName.c_str());
  h1->SetTitle("");
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetXaxis()->SetTitleSize(1.5);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitle("a.u.");
  h1->SetLineWidth(1.5);
  if( varName=="Ech")
    h1->SetFillColor(kRed);
  if( varName=="Egamma")
    h1->SetFillColor(kGreen);
  if( varName=="Enh")
    h1->SetFillColor(kBlue);
  if( varName=="Egammanh")
    h1->SetFillColor(kCyan);

  std::cout << varName << " 50_80 RMS: " << h1->GetRMS() << std::endl;
  c1->cd();
  if(logScale)
    c1->SetLogy();
  h1->DrawNormalized();
  plotName = "Plots/"+histName+".eps";
  c1->SaveAs(plotName.c_str());

  delete c1;

}

 
