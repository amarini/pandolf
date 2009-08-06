#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"



void draw(TFile* file, std::string name, int legendQuadrant);

void drawFlavor() {

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

  TFile* file = TFile::Open("flavorOutputFile.root");


  draw(file, "response_vs_pt_MEAN", 4);
  draw(file, "response_vs_pt_MEAN_barrel", 4);
  draw(file, "response_vs_pt_MEAN_endcap", 4);
  draw(file, "response_vs_pt_MEAN_forward", 4);

  draw(file, "resolution_vs_pt_RMS", 1);
  draw(file, "resolution_vs_pt_RMS_barrel", 1);
  draw(file, "resolution_vs_pt_RMS_endcap", 1);
  draw(file, "resolution_vs_pt_RMS_forward", 1);

}


void draw(TFile* file, std::string name, int legendQuadrant) {

  std::string histName;

  histName = name + "_gluon_PFItCone5";
  TH1F* h1_gluon_PF = (TH1F*)file->Get(histName.c_str());
  h1_gluon_PF->SetMarkerSize(1.5);
  h1_gluon_PF->SetMarkerStyle(20);
  h1_gluon_PF->SetMarkerColor(kRed);
  h1_gluon_PF->SetTitle("");
  h1_gluon_PF->GetXaxis()->SetTitleOffset(1.1);
  h1_gluon_PF->GetYaxis()->SetTitleOffset(1.2);

  histName = name + "_light_PFItCone5";
  TH1F* h1_light_PF = (TH1F*)file->Get(histName.c_str());
  h1_light_PF->SetMarkerSize(1.5);
  h1_light_PF->SetMarkerStyle(21);
  h1_light_PF->SetMarkerColor(kBlue);
  h1_light_PF->SetTitle("");
  h1_light_PF->GetXaxis()->SetTitleOffset(1.1);
  h1_light_PF->GetYaxis()->SetTitleOffset(1.2);


  histName = name + "_gluon_caloItCone5";
  TH1F* h1_gluon_calo = (TH1F*)file->Get(histName.c_str());
  h1_gluon_calo->SetMarkerSize(1.5);
  h1_gluon_calo->SetMarkerStyle(24);
  h1_gluon_calo->SetMarkerColor(kRed);
  h1_gluon_calo->SetTitle("");
  h1_gluon_calo->GetXaxis()->SetTitleOffset(1.1);
  h1_gluon_calo->GetYaxis()->SetTitleOffset(1.2);

  histName = name + "_light_caloItCone5";
  TH1F* h1_light_calo = (TH1F*)file->Get(histName.c_str());
  h1_light_calo->SetMarkerSize(1.5);
  h1_light_calo->SetMarkerStyle(25);
  h1_light_calo->SetMarkerColor(kBlue);
  h1_light_calo->SetTitle("");
  h1_light_calo->GetXaxis()->SetTitleOffset(1.1);
  h1_light_calo->GetYaxis()->SetTitleOffset(1.2);


  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.26, 0.52, 0.74, 0.85);
  if( legendQuadrant==1 )  legend = new TLegend(0.41, 0.52, 0.88, 0.85);
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.65, 0.85);
  if( legendQuadrant==3 )  legend = new TLegend(0.18, 0.2, 0.65, 0.5);
  if( legendQuadrant==4 )  legend = new TLegend(0.41, 0.2, 0.88, 0.5);

  legend->SetTextSize(0.05);
  legend->AddEntry(h1_gluon_PF, "PFJet (gluon)", "P");
  legend->AddEntry(h1_light_PF, "PFJet (light quark)", "P");
  legend->AddEntry(h1_gluon_calo, "CaloJet (gluon)", "P"); 
  legend->AddEntry(h1_light_calo, "CaloJet (light quark)", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h1_gluon_PF->Draw("P");
  h1_light_PF->Draw("Psame");
  h1_gluon_calo->Draw("Psame");
  h1_light_calo->Draw("Psame");
  legend->Draw("same");
  std::string plotName = "Flavor/"+name+".eps";
  c1->SaveAs(plotName.c_str());

}
