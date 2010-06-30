#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.C"


std::string PICTURE_FORMAT = "png";

void draw_Rch(const std::string& dataset, const std::string& name, const std::string& etaRegion, int legendQuadrant);
void drawResponse_PFvsCalo(const std::string& dataset, const std::string& name, const std::string& etaRegion, int legendQuadrant);
void drawResponse_allAlgos(const std::string& dataset, const std::string& name, const std::string& etaRegion, int legendQuadrant);


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



  /*drawResponse_allAlgos(dataset, "response_vs_pt_FIT", "barrel", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_MEAN", "barrel", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_FIT", "endcap", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_MEAN", "endcap", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_FIT", "eta0_25", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_MEAN", "eta0_25", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_FIT", "eta14_25", 4);
  drawResponse_allAlgos(dataset, "response_vs_pt_MEAN", "eta14_25", 4);

  drawResponse_allAlgos(dataset, "resolution_vs_pt_FIT", "barrel", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_RMS", "barrel", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_FIT", "endcap", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_RMS", "endcap", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_FIT", "eta0_25", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_RMS", "eta0_25", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_FIT", "eta14_25", 1);
  drawResponse_allAlgos(dataset, "resolution_vs_pt_RMS", "eta14_25", 1);
*/
  draw_Rch(dataset, "response_vs_pt_FIT",  "barrel",        3);
  draw_Rch(dataset, "response_vs_pt_MEAN", "barrel",        3);
  draw_Rch(dataset, "response_vs_pt_FIT",  "centralBarrel", 3);
  draw_Rch(dataset, "response_vs_pt_MEAN", "centralBarrel", 3);
  draw_Rch(dataset, "response_vs_pt_FIT",  "endcap",        3);
  draw_Rch(dataset, "response_vs_pt_MEAN", "endcap",        3);
  draw_Rch(dataset, "response_vs_pt_FIT",  "eta14_20",      3);
  draw_Rch(dataset, "response_vs_pt_MEAN", "eta14_20",      3);
//  draw_Rch(dataset, "response_vs_pt_FIT", "eta0_25", 4);
//  draw_Rch(dataset, "response_vs_pt_MEAN", "eta0_25", 4);
//  draw_Rch(dataset, "response_vs_pt_FIT", "eta14_25", 4);
//  draw_Rch(dataset, "response_vs_pt_MEAN", "eta14_25", 4);

  draw_Rch(dataset, "responseRes_vs_pt_FIT", "barrel",        0);
  draw_Rch(dataset, "responseRes_vs_pt_RMS", "barrel",        0);
  draw_Rch(dataset, "responseRes_vs_pt_FIT", "centralBarrel", 0);
  draw_Rch(dataset, "responseRes_vs_pt_RMS", "centralBarrel", 0);
  draw_Rch(dataset, "responseRes_vs_pt_FIT", "endcap",        0);
  draw_Rch(dataset, "responseRes_vs_pt_RMS", "endcap",        0);
  draw_Rch(dataset, "responseRes_vs_pt_FIT", "eta14_20",      0);
  draw_Rch(dataset, "responseRes_vs_pt_RMS", "eta14_20",      0);

  draw_Rch(dataset, "deltaEtaRes_vs_pt_FIT", "barrel",        0);
  draw_Rch(dataset, "deltaEtaRes_vs_pt_RMS", "barrel",        0);
  draw_Rch(dataset, "deltaPhiRes_vs_pt_FIT", "barrel",        0);
  draw_Rch(dataset, "deltaPhiRes_vs_pt_RMS", "barrel",        0);

  draw_Rch(dataset, "deltaEtaRes_vs_pt_FIT", "endcap",        0);
  draw_Rch(dataset, "deltaEtaRes_vs_pt_RMS", "endcap",        0);
  draw_Rch(dataset, "deltaPhiRes_vs_pt_FIT", "endcap",        0);
  draw_Rch(dataset, "deltaPhiRes_vs_pt_RMS", "endcap",        0);

//  drawResponse_Rch(dataset, "resolution_vs_pt_FIT", "eta0_25", 1);
//  drawResponse_Rch(dataset, "resolution_vs_pt_RMS", "eta0_25", 1);
//  drawResponse_Rch(dataset, "resolution_vs_pt_FIT", "eta14_25", 1);
//  drawResponse_Rch(dataset, "resolution_vs_pt_RMS", "eta14_25", 1);

//  drawResponse_PFvsCalo(dataset, "response_vs_pt_FIT", "barrel", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_MEAN", "barrel", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_FIT", "endcap", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_MEAN", "endcap", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_FIT", "eta0_25", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_MEAN", "eta0_25", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_FIT", "eta14_25", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_MEAN", "eta14_25", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_FIT", "eta26_32", 4);
//  drawResponse_PFvsCalo(dataset, "response_vs_pt_MEAN", "eta26_32", 4);

} //drawPlots





void draw_Rch(const std::string& dataset, const std::string& name, const std::string& etaRegion, int legendQuadrant) {

  std::string fileName = "finalOutputFile_"+dataset+".root";
  TFile* file = TFile::Open(fileName.c_str());

  Float_t yMin, yMax;
  std::string yName;

  if( name=="response_vs_pt_MEAN" || name=="response_vs_pt_FIT" ) {
    yMin = 0.7;
    yMax = 1.04;
    yName = "p_{T}^{RECO}/p_{T}^{GEN}";
  } else {
    yMin = 0.;
    yMax = 0.35;
    //yName = "#sigma(p_{T}^{RECO}/p_{T}^{GEN}) [%]";
    yName = "resolution [%]";
  }
  
  if( name=="deltaPhiRes_vs_pt_RMS" || name=="deltaPhiRes_vs_pt_FIT" ) {
    yMin = 0.;
    yMax = 0.06;
    yName = "#phi resolution";
  }

  if( name=="deltaEtaRes_vs_pt_RMS" || name=="deltaEtaRes_vs_pt_FIT" ) {
    yMin = 0.;
    yMax = 0.1;
    yName = "#eta resolution";
  }

  TH2F* h2 = new TH2F("h1_aux", "", 10, 10., 2300., 10, yMin, yMax);
  h2->GetYaxis()->SetTitleOffset(1.2);
  h2->GetXaxis()->SetTitleOffset(1.1);
  h2->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2->SetYTitle(yName.c_str());

  std::string gr_name;

/*  gr_name = "gr_" + name + "_"+etaRegion+"_caloakt5";
  TProfile* gr_calo = (TProfile*)file->Get(gr_name.c_str());
  gr_calo->SetLineWidth(2);
  gr_calo->SetMarkerSize(1.5);
  gr_calo->SetLineColor(kBlack);

  gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5_Rch0";
  TGraph* gr_Rch0 = (TGraph*)file->Get(gr_name.c_str());
  gr_Rch0->SetMarkerStyle(20);
  gr_Rch0->SetMarkerSize(1.5);
  gr_Rch0->SetMarkerColor(kGreen);
*/

  gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5_Rch10_30";
  TGraph* gr_Rch10_30 = (TGraph*)file->Get(gr_name.c_str());
  gr_Rch10_30->SetMarkerStyle(21);
  gr_Rch10_30->SetMarkerSize(1.5);
  gr_Rch10_30->SetMarkerColor(kRed);

  gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5_Rch50_70";
  TGraph* gr_Rch50_70 = (TGraph*)file->Get(gr_name.c_str());
  gr_Rch50_70->SetMarkerStyle(22);
  gr_Rch50_70->SetMarkerColor(kMagenta);
  gr_Rch50_70->SetMarkerSize(1.5);

  gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5_Rch90_100";
  TGraph* gr_Rch90_100 = (TGraph*)file->Get(gr_name.c_str());
  gr_Rch90_100->SetMarkerStyle(29);
  gr_Rch90_100->SetMarkerSize(1.5);
  gr_Rch90_100->SetMarkerColor(kBlue);


  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  std::string legendTitle;
  if( etaRegion=="barrel" )
    legendTitle = "   |#eta_{GEN}| < 1.4";
  if( etaRegion=="endcap" )
    legendTitle = "   1.4 < |#eta_{GEN}| < 3.0";
  if( etaRegion=="eta0_25" )
    legendTitle = "   0.0 < |#eta_{GEN}| < 2.5";
  if( etaRegion=="eta14_25" )
    legendTitle = "   1.4 < |#eta_{GEN}| < 2.5";
  if( etaRegion=="eta14_20" )
    legendTitle = "   1.4 < |#eta_{GEN}| < 2.0";
  if( etaRegion=="centralBarrel" )
    legendTitle = "   |#eta_{GEN}| < 1.0";
   

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.32, 0.62, 0.68, 0.88, legendTitle.c_str());
  if( legendQuadrant==1 )  legend = new TLegend(0.52, 0.6, 0.88, 0.88, legendTitle.c_str());
  if( legendQuadrant==2 )  legend = new TLegend(0.2, 0.62, 0.57, 0.88, legendTitle.c_str());
  if( legendQuadrant==3 )  legend = new TLegend(0.17, 0.14, 0.53, 0.43, legendTitle.c_str());
  if( legendQuadrant==4 )  legend = new TLegend(0.46, 0.14, 0.88, 0.45, legendTitle.c_str());

  legend->SetTextSize(0.04);
  //legend->AddEntry(gr_calo, "Calo Jets", "P");
  //legend->AddEntry(gr_Rch0, "PF Jets, R_{ch}=0%", "P");
  legend->AddEntry(gr_Rch10_30, "PF Jets, 10<R_{ch}<30%", "P");
  legend->AddEntry(gr_Rch50_70, "PF Jets, 50<R_{ch}<70%", "P");
  legend->AddEntry(gr_Rch90_100, "PF Jets, R_{ch}>90%", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  TPaveText* label;

  if( (legendQuadrant == 1) || (legendQuadrant==0) )
    label = new TPaveText(0.15, 0.14, 0.3, 0.17, "brNDC");
  else
    label = new TPaveText(0.7, 0.14, 0.85, 0.17, "brNDC");
  label->SetFillColor(kWhite);
  label->SetTextSize(0.04);
  label->AddText("Anti-kt 0.5");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogx();
  h2->Draw();
  //gr_calo->Draw("*same"); 
  //gr_Rch0->Draw("Psame"); 
  gr_Rch10_30->Draw("Psame"); 
  gr_Rch50_70->Draw("Psame"); 
  gr_Rch90_100->Draw("Psame"); 
  legend->Draw("same");
  label->Draw("same");

  std::string canvasName = "Plots_"+dataset+"/gr_Rch_"+name+"_"+etaRegion+"."+PICTURE_FORMAT;
  //std::string canvasName = "Philipp/Rch_"+name+"_"+etaRegion+".pdf";
  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h2;

}


void drawResponse_PFvsCalo(const std::string& dataset, const std::string& name, const std::string& etaRegion, int legendQuadrant) {

  std::string fileName = "finalOutputFile_"+dataset+".root";
  TFile* file = TFile::Open(fileName.c_str());

  TH2F* h2 = new TH2F("h1_aux", "", 10, 20., 600., 10, 0.2, 1.05);
  h2->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  std::string gr_name;

  gr_name = "gr_" + name + "_"+etaRegion+"_caloItCone5";
  TProfile* gr_calo = (TProfile*)file->Get(gr_name.c_str());
  gr_calo->SetLineWidth(2);
  gr_calo->SetMarkerSize(1.5);
  gr_calo->SetLineColor(kBlack);

  if( etaRegion=="eta26_30" )
    gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5_Rch0";
  else 
    gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5";
  TGraph* gr_PF = (TGraph*)file->Get(gr_name.c_str());
  gr_PF->SetMarkerStyle(20);
  gr_PF->SetMarkerSize(1.5);
  gr_PF->SetMarkerColor(kGreen);

  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  std::string legendTitle;
  if( etaRegion=="barrel" )
    legendTitle = "   0.0<|#eta_{GEN}|<1.4";
  if( etaRegion=="endcap" )
    legendTitle = "   1.4<|#eta_{GEN}|<3";
  if( etaRegion=="eta0_25" )
    legendTitle = "   0.0<|#eta_{GEN}|<2.5";
  if( etaRegion=="eta14_25" )
    legendTitle = "   1.4<|#eta_{GEN}|<2.5";
  if( etaRegion=="eta26_32" )
    legendTitle = "   2.6<|#eta_{GEN}|<3.2";
   

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.72, 0.85, legendTitle.c_str());
  if( legendQuadrant==1 )  legend = new TLegend(0.44, 0.52, 0.88, 0.85, legendTitle.c_str());
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.62, 0.85, legendTitle.c_str());
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.14, 0.6 , 0.5, legendTitle.c_str());
  if( legendQuadrant==4 )  legend = new TLegend(0.52, 0.14, 0.88, 0.38, legendTitle.c_str());

  legend->SetTextSize(0.045);
  legend->AddEntry(gr_calo, "Calo Jets", "P");
  if( etaRegion=="eta26_30" )
    legend->AddEntry(gr_PF, "PF Jets, R_{ch}=0%", "P");
  else
    legend->AddEntry(gr_PF, "PF Jets", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  TPaveText* label;

  if(legendQuadrant==4) {
    label = new TPaveText(0.13, 0.14, 0.35, 0.17, "brNDC");
    label->SetFillColor(kWhite);
    label->SetTextSize(0.04);
    label->AddText("Iterative Cone 0.5");
  }

  TPaveText* label2;

  if(legendQuadrant==4) {
    label2 = new TPaveText(0.13, 0.18, 0.25, 0.21, "brNDC");
    //label2 = new TPaveText(0.7, 0.53, 0.85, 0.55, "brNDC");
    label2->SetFillColor(kWhite);
    label2->SetTextSize(0.04);
 }


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h2->Draw();
  gr_calo->Draw("*"); 
  gr_PF->Draw("Psame"); 
  legend->Draw("same");
  label->Draw("same");
  label2->Draw("same");

  std::string canvasName = "Plots_"+dataset+"/gr_PFvsCalo_"+name+"_"+etaRegion+"."+PICTURE_FORMAT;
  //std::string canvasName = "Philipp/Rch_"+name+"_"+etaRegion+".pdf";
  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h2;

}


void drawResponse_allAlgos(const std::string& dataset, const std::string& name, const std::string& etaRegion, int legendQuadrant) {

  std::string fileName = "finalOutputFile_"+dataset+".root";
  TFile* file = TFile::Open(fileName.c_str());


  Float_t yMin, yMax;

  if( name=="response_vs_pt_MEAN" || name=="response_vs_pt_FIT" ) {
    yMin = 0.2;
    yMax = 1.05;
  }
  if( name=="resolution_vs_pt_RMS" || name=="resolution_vs_pt_FIT" ) {
    yMin = 0.;
    yMax = 0.5;
  }
  
  TH2F* h2 = new TH2F("h1_aux", "", 10, 20., 600., 10, yMin, yMax);
  h2->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2->SetYTitle("p_{T}^{RECO}/p_{T}^{GEN}");

  std::string gr_name;

  gr_name = "gr_" + name + "_"+etaRegion+"_caloItCone5";
  TGraph* gr_caloItCone5 = (TGraph*)file->Get(gr_name.c_str());
  gr_caloItCone5->SetMarkerStyle(20);
  gr_caloItCone5->SetMarkerSize(1.5);
  gr_caloItCone5->SetMarkerColor(kRed);

  gr_name = "gr_" + name + "_"+etaRegion+"_caloSISCone5";
  TGraph* gr_caloSISCone5 = (TGraph*)file->Get(gr_name.c_str());
  gr_caloSISCone5->SetMarkerStyle(21);
  gr_caloSISCone5->SetMarkerSize(1.5);
  gr_caloSISCone5->SetMarkerColor(kGreen);

  gr_name = "gr_" + name + "_"+etaRegion+"_calokt4";
  TGraph* gr_calokt4 = (TGraph*)file->Get(gr_name.c_str());
  gr_calokt4->SetMarkerStyle(22);
  gr_calokt4->SetMarkerColor(kBlue);
  gr_calokt4->SetMarkerSize(1.5);

  gr_name = "gr_" + name + "_"+etaRegion+"_PFakt5";
  TGraph* gr_PFakt5 = (TGraph*)file->Get(gr_name.c_str());
  gr_PFakt5->SetMarkerStyle(29);
  gr_PFakt5->SetMarkerSize(1.5);
  gr_PFakt5->SetMarkerColor(kMagenta);


  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.68, 0.85);
  if( legendQuadrant==1 )  legend = new TLegend(0.48, 0.52, 0.88, 0.85);
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.58, 0.85);
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.2, 0.55, 0.5);
  if( legendQuadrant==4 )  legend = new TLegend(0.48, 0.2, 0.88, 0.5);

  legend->SetTextSize(0.05);
  legend->AddEntry(gr_caloItCone5, "Iterative Cone 0.5", "P");
  legend->AddEntry(gr_caloSISCone5, "SIS Cone 0.5", "P");
  legend->AddEntry(gr_calokt4, "k_{T} D=0.4", "P");
  legend->AddEntry(gr_PFakt5, "PF anti-kt 0.5", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h2->Draw();
  gr_caloItCone5->Draw("Psame"); 
  gr_caloSISCone5->Draw("Psame"); 
  gr_calokt4->Draw("Psame"); 
  gr_PFakt5->Draw("Psame"); 

  std::string canvasName = "Plots_"+dataset+"/gr_"+name+"_"+etaRegion+"."+PICTURE_FORMAT;
  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h2;

}

