#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <string>
#include <iostream>


void compare_variable(std::string varName, std::string etaRegion, int legendQuadrant);

void compare_fit() {

  compare_variable("response_vs_pt", "barrel", 4);
  compare_variable("response_vs_pt", "endcap", 4);
  compare_variable("resolution_vs_pt", "barrel", 1);
  compare_variable("resolution_vs_pt", "endcap", 1);


}


void compare_variable(std::string varName, std::string etaRegion, int legendQuadrant) {


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

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();


  std::string var_fit = varName + "_FIT_" + etaRegion+"_PFItCone5";
  std::string var_nofit;
  if( varName == "resolution_vs_pt" ) 
    var_nofit = varName + "_RMS_"+etaRegion+"_PFItCone5";
  else
    var_nofit = varName + "_MEAN_"+etaRegion+"_PFItCone5";
  


  TFile* file= TFile::Open("finalOutputFile_226.root");



  TH1F* h1_fit = (TH1F*)file->Get(var_fit.c_str());
  h1_fit->SetMarkerStyle(20);
  h1_fit->SetMarkerSize(1.5);
  h1_fit->SetMarkerColor(kBlue);
  h1_fit->SetTitle("");
  h1_fit->GetXaxis()->SetTitleOffset(1.1);
  h1_fit->GetYaxis()->SetTitleOffset(1.2);
  h1_fit->Draw("P");

  TH1F* h1_nofit = (TH1F*)file->Get(var_nofit.c_str());
  h1_nofit->SetMarkerStyle(24);
  h1_nofit->SetMarkerSize(1.5);
  h1_nofit->SetMarkerColor(kBlue);
  h1_nofit->SetTitle("");
  h1_nofit->GetXaxis()->SetTitleOffset(1.1);
  h1_nofit->GetYaxis()->SetTitleOffset(1.2);
  h1_nofit->Draw("Psame");


  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.3, 0.5, 0.68, 0.85);
  if( legendQuadrant==1 )  legend = new TLegend(0.62, 0.6, 0.88, 0.8);
  if( legendQuadrant==2 )  legend = new TLegend(0.15, 0.6, 0.4, 0.8);
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.2, 0.4, 0.4);
  if( legendQuadrant==4 )  legend = new TLegend(0.62, 0.2, 0.88, 0.4);

  legend->SetTextSize(0.04);
  legend->AddEntry(h1_fit, "Gaussian Fit", "P");
  if( varName == "resolution_vs_pt" )
    legend->AddEntry(h1_nofit, "RMS", "P");
  else
    legend->AddEntry(h1_nofit, "Mean", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->Draw("same");

  std::string plotName = "compareFIT_"+varName+"_"+etaRegion+".eps";
  c1->SaveAs(plotName.c_str());

  delete c1;

}
