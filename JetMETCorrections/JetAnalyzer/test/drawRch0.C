#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include <string>


void draw(TFile* file, TFile* controlFile, std::string name, int legendQuadrant);


void drawRch0() {

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


  TFile* file = TFile::Open("Rch0.root");
  TFile* controlFile = TFile::Open("finalOutputFile.root");

  draw(file, controlFile, "response_vs_pt_MEAN_barrel", 4);
  draw(file, controlFile, "response_vs_pt_MEAN_endcap", 4);
  draw(file, controlFile, "response_vs_pt_MEAN_forward", 4);
  draw(file, controlFile, "resolution_vs_pt_RMS_barrel", 1);
  draw(file, controlFile, "resolution_vs_pt_RMS_endcap", 1);
  draw(file, controlFile, "resolution_vs_pt_RMS_forward", 1);

}


void draw(TFile* file, TFile* controlFile, std::string name, int legendQuadrant) {

  std::string histName_PF = name+"_PF";
  TH1F* h1_PF = (TH1F*)file->Get(histName_PF.c_str());
  h1_PF->SetTitle("");
  h1_PF->GetXaxis()->SetTitleOffset(1.1);
  h1_PF->GetYaxis()->SetTitleOffset(1.1);
  h1_PF->SetMarkerSize(1.2);
  h1_PF->SetMarkerStyle(20);
  h1_PF->SetMarkerColor(kGreen);

  std::string histName_calo = name+"_calo";
  TH1F* h1_calo = (TH1F*)file->Get(histName_calo.c_str());
  h1_calo->SetTitle("");
  h1_calo->GetXaxis()->SetTitleOffset(1.1);
  h1_calo->GetYaxis()->SetTitleOffset(1.1);
  h1_calo->SetMarkerSize(1.2);
  h1_calo->SetMarkerStyle(21);
  h1_calo->SetMarkerColor(kRed);

  std::string histName_control = name+"_caloItCone5";
  TH1F* h1_control = (TH1F*)controlFile->Get(histName_control.c_str());
  h1_control->SetTitle("");
  h1_control->GetXaxis()->SetTitleOffset(1.1);
  h1_control->GetYaxis()->SetTitleOffset(1.1);
  h1_control->SetLineWidth(2.);
  h1_control->SetMarkerSize(0.);
  h1_control->SetLineColor(kBlack);

  Int_t nBins = h1_PF->GetNbinsX();

  for(int iBin=1; iBin<(nBins+1); ++iBin) {
    if( h1_PF->GetBinContent(iBin)<0.0000001 ) h1_PF->SetBinContent(iBin, -1.);
    if( h1_calo->GetBinContent(iBin)<0.0000001 ) h1_calo->SetBinContent(iBin, -1.);
    if( h1_control->GetBinContent(iBin)<0.0000001 ) h1_control->SetBinContent(iBin, -1.);
  }

  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.68, 0.85);
  if( legendQuadrant==1 )  legend = new TLegend(0.48, 0.52, 0.88, 0.85);
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.58, 0.85);
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.2, 0.55, 0.5);
  if( legendQuadrant==4 )  legend = new TLegend(0.48, 0.2, 0.88, 0.5);

  legend->SetTextSize(0.05);
  legend->AddEntry(h1_PF, "PFlow, R_{ch}=0", "P");
  legend->AddEntry(h1_calo, "Calo, \"R_{ch}=0\"", "P");
  legend->AddEntry(h1_control, "Calo, all events", "L"); 
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd(0);
  h1_PF->Draw("p");
  h1_calo->Draw("psame");
  h1_control->Draw("lsame");
  legend->Draw("same");
  std::string saveName = "Rch0_vs_calo_"+name+".eps";
  c1->SaveAs(saveName.c_str());
  
}  
