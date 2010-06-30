#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.C"


std::string PICTURE_FORMAT = "eps";


void draw(const std::string& dataset, const std::string& name, const std::string& etaRegion, const std::string& algoType, int legendQuadrant, const std::string& suffix);
void draw_histo(const std::string& dataset, const std::string& name, const std::string& algoType, int legendQuadrant);


void drawTGraphs(std::string dataset) {

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

  std::string algoType = "akt5";

  draw(dataset, "response_vs_pt_FIT", "barrel", algoType, 4, "RECO");
  //draw(dataset, "response_vs_pt_FIT", "eta14_25", algoType, 4, "RECO");
  draw(dataset, "response_vs_pt_FIT", "endcap", algoType, 4, "RECO");
  //draw(dataset, "response_vs_pt_FIT", "eta0_3", algoType, 4, "RECO");
  draw(dataset, "response_vs_pt_MEAN", "barrel", algoType, 4, "RECO");
  //draw(dataset, "response_vs_pt_MEAN", "eta14_25", algoType, 4, "RECO");
  draw(dataset, "response_vs_pt_MEAN", "endcap", algoType, 4, "RECO");
  //draw(dataset, "response_vs_pt_MEAN", "eta0_3", algoType, 4, "RECO");

  draw(dataset, "response_vs_pt_FIT", "barrel", algoType, 4, "GEN");
  draw(dataset, "response_vs_pt_FIT", "endcap", algoType, 4, "GEN");
  draw(dataset, "response_vs_pt_FIT", "eta14_20", algoType, 4, "GEN");
  draw(dataset, "response_vs_pt_MEAN", "barrel", algoType, 4, "GEN");
  draw(dataset, "response_vs_pt_MEAN", "endcap", algoType, 4, "GEN");
  draw(dataset, "response_vs_pt_MEAN", "eta14_20", algoType, 4, "GEN");

  draw(dataset, "responseRes_vs_pt_FIT", "barrel", algoType, 1, "RECO");
  draw(dataset, "responseRes_vs_pt_FIT", "eta14_20", algoType, 1, "RECO");
  draw(dataset, "responseRes_vs_pt_FIT", "endcap", algoType, 1, "RECO");
//draw(dataset, "responseRes_vs_pt_FIT", "eta0_3", algoType, 1, "RECO");
  draw(dataset, "responseRes_vs_pt_RMS", "barrel", algoType, 1, "RECO");
  draw(dataset, "responseRes_vs_pt_RMS", "eta14_20", algoType, 1, "RECO");
  draw(dataset, "responseRes_vs_pt_RMS", "endcap", algoType, 1, "RECO");
//draw(dataset, "responseRes_vs_pt_RMS", "eta0_3", algoType, 1, "RECO");

  draw(dataset, "responseRes_vs_pt_FIT", "barrel", algoType, 1, "GEN");
  draw(dataset, "responseRes_vs_pt_FIT", "endcap", algoType, 1, "GEN");
  draw(dataset, "responseRes_vs_pt_FIT", "eta14_20", algoType, 1, "GEN");
  draw(dataset, "responseRes_vs_pt_RMS", "barrel", algoType, 1, "GEN");
  draw(dataset, "responseRes_vs_pt_RMS", "endcap", algoType, 1, "GEN");
  draw(dataset, "responseRes_vs_pt_RMS", "eta14_20", algoType, 1, "GEN");

std::cout << "e fin qui tutto bene." << std::endl;
  draw_histo(dataset, "response_vs_eta", algoType, 5);

} //drawPlots




void draw(const std::string& dataset, const std::string& name, const std::string& etaRegion, const std::string& algoType, int legendQuadrant, const std::string& suffix) {


  std::string fileName= "RchCalibration_"+dataset+".root";
  TFile* file= TFile::Open(fileName.c_str());

  Float_t yMin, yMax;
  std::string yName;
  std::string ptlabel = (suffix=="RECO") ? "#gamma" : "GEN"; 
  std::string etalabel = (suffix=="RECO") ? "RECO" : "GEN"; 

  if( (name=="response_vs_pt_MEAN")||(name=="response_vs_pt_FIT") ) {
    yMin = (suffix=="GEN") ? 0.8 : 0.6;
    yMax = 1.05;
    yName = "p_{T}^{RECO}/p_{T}^{" + ptlabel + "}";
  } else {
    yMin = 0.;
    yMax = (suffix=="GEN") ? 0.3 : 0.4;
    //yName = "#sigma(p_{T}^{RECO}/p_{T}^{" + ptlabel + "}) [%]";
    yName = "resolution [%]";
  }

  std::string xName = "p_{T}^{" + ptlabel + "} [GeV/c]";

  TH2F* h2 = new TH2F("h1_aux", "", 10, 10., 2300., 10, yMin, yMax);
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->GetXaxis()->SetTitleOffset(1.1);
  h2->SetXTitle(xName.c_str());
  h2->SetYTitle(yName.c_str());

  std::string gr_name;

  gr_name = "gr_" + name + "_"+etaRegion+"_PF" + algoType;
  if( suffix!="RECO" )
    gr_name = gr_name + "_" + suffix;
  gr_name += "_PRE";
  TGraph* gr_preCorr = (TGraph*)file->Get(gr_name.c_str());
  gr_preCorr->SetMarkerStyle(21);
  gr_preCorr->SetMarkerSize(0.8);
  gr_preCorr->SetLineColor(kBlack);

  gr_name = "gr_" + name + "_"+etaRegion+"_PF" + algoType;
  if( suffix!="RECO" )
    gr_name = gr_name + "_" + suffix;
  TGraph* gr_postCorr = (TGraph*)file->Get(gr_name.c_str());
  gr_postCorr->SetMarkerStyle(20);
  gr_postCorr->SetMarkerSize(1.5);
  gr_postCorr->SetMarkerColor(kRed);


  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  std::string legendTitle;
  if( etaRegion=="barrel" )
    legendTitle = "   |#eta^{" + etalabel + "}| < 1.4";
  if( etaRegion=="endcap" )
    legendTitle = "   1.4 < |#eta^{" + etalabel + "}| < 3";
  if( etaRegion=="eta0_25" )
    legendTitle = "   |#eta^{" + etalabel + "}|< 2.5";
  if( etaRegion=="eta14_20" )
    legendTitle = "   1.4 < |#eta^{" + etalabel + "}|< 2.0";
  if( etaRegion=="eta14_25" )
    legendTitle = "   1.4 < |#eta^{" + etalabel + "}|< 2.5";
  if( etaRegion=="eta0_3" )
    legendTitle =       "   |#eta^{" + etalabel + "}| < 3";
   

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.72, 0.85, legendTitle.c_str());
  if( legendQuadrant==1 )  legend = new TLegend(0.44, 0.52, 0.88, 0.85, legendTitle.c_str());
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.62, 0.85, legendTitle.c_str());
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.14, 0.6 , 0.5, legendTitle.c_str());
  if( legendQuadrant==4 )  legend = new TLegend(0.46, 0.14, 0.88, 0.45, legendTitle.c_str());

  legend->SetTextSize(0.045);
  legend->AddEntry(gr_preCorr, "Before Correction", "P");
  legend->AddEntry(gr_postCorr, "After Correction", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  std::string algoType_text;
  if( algoType=="akt5" )
    algoType_text = "Anti-kt 0.5";
  else if( algoType=="ite5" )
    algoType_text = "It. Cone 0.5";

  TPaveText* label;
  if(legendQuadrant==4) {
    label = new TPaveText(0.13, 0.14, 0.35, 0.17, "brNDC");
    label->SetFillColor(kWhite);
    label->SetTextSize(0.04);
    label->AddText(algoType_text.c_str());
  } else {
    label = new TPaveText(0.13, 0.85, 0.35, 0.88, "brNDC");
    label->SetFillColor(kWhite);
    label->SetTextSize(0.04);
    label->AddText(algoType_text.c_str());
  }


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogx();
  h2->Draw();
  gr_postCorr->Draw("P"); 
  gr_preCorr->Draw("Psame"); 
  legend->Draw("same");
  label->Draw("same");

  std::string canvasName = "PlotsRchCalib_"+dataset+"/gr_RchCalib_"+name+"_"+etaRegion+"_"+suffix+"."+PICTURE_FORMAT;
  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h2;

}



void draw_histo(const std::string& dataset, const std::string& name, const std::string& algoType, int legendQuadrant) {


  std::string fileName= "RchCalibration_"+dataset+".root";
  TFile* file= TFile::Open(fileName.c_str());

  std::string hist_name;

  hist_name =  name + "_PRE_PF" + algoType;
  TH1F* hist_preCorr = (TH1F*)file->Get(hist_name.c_str());
  hist_preCorr->SetMarkerStyle(21);
  hist_preCorr->SetMarkerSize(0.8);
  hist_preCorr->SetLineColor(kBlack);
  hist_preCorr->SetTitle("");
  hist_preCorr->GetYaxis()->SetRangeUser(0.2, 1.15);

  hist_name =  name + "_PF" + algoType;
  TH1F* hist_postCorr = (TH1F*)file->Get(hist_name.c_str());
  hist_postCorr->SetMarkerStyle(20);
  hist_postCorr->SetMarkerSize(1.5);
  hist_postCorr->SetMarkerColor(kRed);
  hist_postCorr->SetTitle("");
  hist_postCorr->GetYaxis()->SetRangeUser(0.2, 1.15);



  std::string legendTitle = "   p_{T}^{GEN} > 20 GeV/c";
   

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.72, 0.85, legendTitle.c_str());
  if( legendQuadrant==1 )  legend = new TLegend(0.44, 0.52, 0.88, 0.85, legendTitle.c_str());
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.62, 0.85, legendTitle.c_str());
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.14, 0.6 , 0.5, legendTitle.c_str());
  if( legendQuadrant==4 )  legend = new TLegend(0.46, 0.14, 0.88, 0.45, legendTitle.c_str());
  if( legendQuadrant==5 )  legend = new TLegend(0.3, 0.14, 0.7, 0.43, legendTitle.c_str());

  legend->SetTextSize(0.045);
  legend->AddEntry(hist_preCorr, "Before Correction", "P");
  legend->AddEntry(hist_postCorr, "After Correction", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  std::string algoType_text;
  if( algoType=="akt5" )
    algoType_text = "Anti-kt 0.5";
  else if( algoType=="ite5" )
    algoType_text = "It. Cone 0.5";

  TPaveText* label;
  if(legendQuadrant==4) {
    label = new TPaveText(0.13, 0.14, 0.35, 0.17, "brNDC");
    label->SetFillColor(kWhite);
    label->SetTextSize(0.04);
    label->AddText(algoType_text.c_str());
  } else {
    label = new TPaveText(0.13, 0.85, 0.35, 0.88, "brNDC");
    label->SetFillColor(kWhite);
    label->SetTextSize(0.04);
    label->AddText(algoType_text.c_str());
  }


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  hist_postCorr->Draw("P"); 
  hist_preCorr->Draw("Psame"); 
  legend->Draw("same");
  label->Draw("same");

  std::string canvasName = "PlotsRchCalib_"+dataset+"/"+name+".eps";
  c1->SaveAs(canvasName.c_str());

  delete c1;

}


