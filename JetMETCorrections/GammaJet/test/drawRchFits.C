#include "TFile.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TPaveText.h"


void draw(const std::string& dataset, const std::string& etaRegion, const std::string& recoGen);


void drawRchFits(std::string dataset) {

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


  //draw(dataset, "barrel", "RECO");
  //draw(dataset, "endcap", "RECO");
  //draw(dataset, "eta14_25", "RECO");
  //draw(dataset, "eta155_3", "RECO");

  draw(dataset, "barrel", "GEN");
  draw(dataset, "endcap", "GEN");
//  draw(dataset, "eta14_25", "GEN");
//  draw(dataset, "eta155_3", "GEN");

  draw(dataset, "barrel", "GENGEN");
  draw(dataset, "endcap", "GENGEN");
}


void draw(const std::string& dataset, const std::string& etaRegion, const std::string& recoGen) {

  int nBins = 26;

  Double_t Lower[nBins];

  Int_t i=0;
  Lower[i++] = 10.;
  Lower[i++] = 15.;
  Lower[i++] = 20.;
  Lower[i++] = 25.;
  Lower[i++] = 30.;
  Lower[i++] = 36.;
  Lower[i++] = 42.;
  Lower[i++] = 48.;
  Lower[i++] = 55.;
  Lower[i++] = 62.;
  Lower[i++] = 70.;
  Lower[i++] = 78.;
  Lower[i++] = 87.;
  Lower[i++] = 97.;
  Lower[i++] = 108.;
  Lower[i++] = 121.;
  Lower[i++] = 135.;
  Lower[i++] = 152.;
  Lower[i++] = 172.;
  Lower[i++] = 197.;
  Lower[i++] = 231.;
  Lower[i++] = 283.;
  Lower[i++] = 397.;
  Lower[i++] = 600.;
  Lower[i++] = 1000.;
  Lower[i++] = 2000.;


  //pt bins to graph:
  Int_t bin1 = 2;
  Int_t bin2 = 10;
  Int_t bin3 = 15;
  Int_t bin4 = 23;

  std::string fileName = etaRegion;
  if( recoGen!="RECO" )
    fileName = fileName+recoGen;
  std::string filePath = "Projections_RchCalib_"+dataset+"/"+fileName+".root";
  TFile* file = TFile::Open(filePath.c_str());

  Float_t yMin = (recoGen=="RECO") ? 0.6 : 0.8;
  Float_t yMax = 1.05;

  std::string xLabel = (recoGen=="GENGEN") ? "R_{ch}^{GEN}" : "R_{ch}";
  std::string yLabel = "response";

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1., 10, yMin, yMax);
  h2_axes->SetXTitle(xLabel.c_str());
  h2_axes->SetYTitle(yLabel.c_str());
  
  char fitName[100];

  sprintf(fitName, "ptbin_%d/fitResponse_%d", bin1, bin1);
  TGraphErrors* gr_fit1 = (TGraphErrors*)file->Get(fitName);
  gr_fit1->SetMarkerStyle(20);
  gr_fit1->SetMarkerColor(kGreen);
  gr_fit1->SetMarkerSize(0.8);
  gr_fit1->SetTitle("");
  gr_fit1->GetYaxis()->SetRangeUser(yMin, yMax);
  gr_fit1->GetFunction("func")->SetLineColor(kGreen);
  gr_fit1->GetFunction("func")->SetLineWidth(0.7);
  gr_fit1->GetYaxis()->SetTitleOffset(1.2);
  
  sprintf(fitName, "ptbin_%d/fitResponse_%d", bin2, bin2);
  TGraphErrors* gr_fit2 = (TGraphErrors*)file->Get(fitName);
  gr_fit2->SetMarkerStyle(21);
  gr_fit2->SetMarkerColor(kBlue);
  gr_fit2->SetMarkerSize(0.8);
  gr_fit2->SetTitle("");
  gr_fit2->GetYaxis()->SetRangeUser(yMin, yMax);
  gr_fit2->GetFunction("func")->SetLineColor(kBlue);
  gr_fit2->GetFunction("func")->SetLineWidth(0.7);
  gr_fit2->GetYaxis()->SetTitleOffset(1.2);

  sprintf(fitName, "ptbin_%d/fitResponse_%d", bin3, bin3);
  TGraphErrors* gr_fit3 = (TGraphErrors*)file->Get(fitName);
  gr_fit3->SetMarkerStyle(22);
  gr_fit3->SetMarkerColor(kBlack);
  gr_fit3->SetMarkerSize(0.8);
  gr_fit3->SetTitle("");
  gr_fit3->GetYaxis()->SetRangeUser(yMin, yMax);
  gr_fit3->GetFunction("func")->SetLineColor(kBlack);
  gr_fit3->GetFunction("func")->SetLineWidth(0.7);
  gr_fit3->GetYaxis()->SetTitleOffset(1.2);

  sprintf(fitName, "ptbin_%d/fitResponse_%d", bin4, bin4);
  TGraphErrors* gr_fit4 = (TGraphErrors*)file->Get(fitName);
  gr_fit4->SetMarkerStyle(23);
  gr_fit4->SetMarkerColor(kRed);
  gr_fit4->SetMarkerSize(0.8);
  gr_fit4->SetTitle("");
  gr_fit4->GetYaxis()->SetRangeUser(yMin, yMax);
  gr_fit4->GetFunction("func")->SetLineColor(kRed);
  gr_fit4->GetFunction("func")->SetLineWidth(0.7);
  gr_fit4->GetYaxis()->SetTitleOffset(1.2);

  
  std::string suffix = (recoGen=="RECO") ? "PHOT" : "GEN";
  std::string legendTitle;
  if( etaRegion=="barrel" )
    legendTitle = "   |#eta^{" + suffix + "}| < 1.4";
  if( etaRegion=="endcap" )
    legendTitle = "   1.4 < |#eta^{" + suffix + "}| < 3.0";
  if( etaRegion=="eta14_25" )
    legendTitle = "   1.4 < |#eta^{" + suffix + "}| < 2.5";
  if( etaRegion=="eta155_3" )
    legendTitle = "   1.55 < |#eta^{" + suffix + "}| < 3.0";

  TLegend* legend = new TLegend(0.53, 0.13, 0.88, 0.4, legendTitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.033);

  char legendName1[100];
  char legendName2[100];
  char legendName3[100];
  char legendName4[100];

  sprintf( legendName1, "%d < p_{T}^{%s} < %d GeV/c", (int)Lower[bin1-1], suffix.c_str(), (int)Lower[bin1]);
  sprintf( legendName2, "%d < p_{T}^{%s} < %d GeV/c", (int)Lower[bin2-1], suffix.c_str(), (int)Lower[bin2]);
  sprintf( legendName3, "%d < p_{T}^{%s} < %d GeV/c", (int)Lower[bin3-1], suffix.c_str(), (int)Lower[bin3]);
  sprintf( legendName4, "%d < p_{T}^{%s} < %d GeV/c", (int)Lower[bin4-1], suffix.c_str(), (int)Lower[bin4]);

  legend->AddEntry(gr_fit1, legendName1, "P");
  legend->AddEntry(gr_fit2, legendName2, "P");
  if( etaRegion!="barrel" )
    legend->AddEntry(gr_fit3, legendName3, "P");
  legend->AddEntry(gr_fit4, legendName4, "P");
  
  TPaveText* label = new TPaveText(0.22, 0.85, 0.25, 0.88, "brNDC");
  label->SetFillColor(kWhite);
  label->SetTextSize(0.04);
  label->AddText("Anti-kt 0.5");
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  h2_axes->Draw("");
  gr_fit1->Draw("Psame");
  gr_fit2->Draw("Psame");
  if( etaRegion!="barrel" )
    gr_fit3->Draw("Psame");
  gr_fit4->Draw("Psame");
  legend->Draw("same");
  label->Draw("same");

  std::string canvasName = "PlotsRchCalib_"+dataset+"/fits_"+fileName+".eps";
  c1->SaveAs(canvasName.c_str());
  delete c1;

}
