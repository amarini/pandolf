#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include <iostream>


bool DEBUG_ = false;


void draw(const std::string& dataset, const std::string& name, int legendQuadrant=1, bool normalized=false);
void draw_Rch(const std::string& dataset, const std::string& name, int legendQuadrant, bool normalized=false);
void draw_Rgamma(const std::string& dataset, const std::string& name, int legendQuadrant, bool normalized=false);
void draw_variable(const std::string& dataset, const std::string& name, bool normalized=false, bool setLogY=false);
void draw_variable2D(const std::string& dataset, const std::string& name);


void drawPlots(const std::string& dataset="") {

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


  bool drawNormalized=true;
  bool setLogY=true;

/*  draw(dataset, "deltaR_barrelEndcap", 1, drawNormalized);
  draw(dataset, "deltaR_forward", 1, drawNormalized);

  draw(dataset, "response_vs_pt_FIT_barrel", 4);
  draw(dataset, "response_vs_pt_FIT_endcap", 4);
  draw(dataset, "response_vs_pt_FIT_forward", 4);

  draw(dataset, "response_vs_pt_MEAN_barrel", 4);
  draw(dataset, "response_vs_pt_MEAN_endcap", 4);
  draw(dataset, "response_vs_pt_MEAN_forward", 4);

  draw(dataset, "response_vs_eta_FIT", 4);
  draw(dataset, "resolution_vs_eta_FIT", 0);

  draw(dataset, "response_vs_eta_MEAN", 4);
  draw(dataset, "resolution_vs_eta_RMS", 0);

  draw(dataset, "resolution_vs_pt_FIT_barrel", 1);
  draw(dataset, "resolution_vs_pt_FIT_endcap", 1);
  draw(dataset, "resolution_vs_pt_FIT_forward", 1);

  draw(dataset, "resolution_vs_pt_RMS_barrel", 1);
  draw(dataset, "resolution_vs_pt_RMS_endcap", 1);
  draw(dataset, "resolution_vs_pt_RMS_forward", 1);

  draw(dataset, "sigmaEta_vs_pt_barrel", 1);
  draw(dataset, "sigmaEta_vs_pt_endcap", 1);
  draw(dataset, "sigmaEta_vs_pt_forward", 1);

  draw(dataset, "sigmaPhi_vs_pt_barrel", 1);
  draw(dataset, "sigmaPhi_vs_pt_endcap", 1);
  draw(dataset, "sigmaPhi_vs_pt_forward", 1);
*/
  draw_Rch(dataset, "response_vs_pt_FIT_barrel", 4);
  draw_Rch(dataset, "response_vs_pt_FIT_endcap", 4);
  draw_Rch(dataset, "response_vs_pt_MEAN_barrel", 4);
  draw_Rch(dataset, "response_vs_pt_MEAN_endcap", 4);

  //draw_Rch(dataset, "response_vs_pt_FIT_eta0_25", 4);
  //draw_Rch(dataset, "response_vs_pt_FIT_eta14_25", 4);
  //draw_Rch(dataset, "response_vs_pt_MEAN_eta0_25", 4);
  //draw_Rch(dataset, "response_vs_pt_MEAN_eta14_25", 4);

  draw_Rch(dataset, "resolution_vs_pt_FIT_barrel", 1);
  draw_Rch(dataset, "resolution_vs_pt_FIT_endcap", 1);
  draw_Rch(dataset, "resolution_vs_pt_RMS_barrel", 1);
  draw_Rch(dataset, "resolution_vs_pt_RMS_endcap", 1);

  //draw_Rch(dataset, "resolution_vs_pt_FIT_eta0_25", 1);
  //draw_Rch(dataset, "resolution_vs_pt_FIT_eta14_25", 1);
  //draw_Rch(dataset, "resolution_vs_pt_RMS_eta0_25", 1);
  //draw_Rch(dataset, "resolution_vs_pt_RMS_eta14_25", 1);

  draw_Rch(dataset, "sigmaEta_vs_pt_barrel", 1);
  draw_Rch(dataset, "sigmaEta_vs_pt_endcap", 1);
  draw_Rch(dataset, "sigmaPhi_vs_pt_barrel", 1);
  draw_Rch(dataset, "sigmaPhi_vs_pt_endcap", 1);

  draw_Rgamma(dataset, "response_vs_pt_FIT_barrel", 4);
  draw_Rgamma(dataset, "response_vs_pt_FIT_endcap", 4);
  draw_Rgamma(dataset, "response_vs_pt_MEAN_barrel", 4);
  draw_Rgamma(dataset, "response_vs_pt_MEAN_endcap", 4);

  draw_Rgamma(dataset, "resolution_vs_pt_FIT_barrel", 1);
  draw_Rgamma(dataset, "resolution_vs_pt_FIT_endcap", 1);
  draw_Rgamma(dataset, "resolution_vs_pt_RMS_barrel", 1);
  draw_Rgamma(dataset, "resolution_vs_pt_RMS_endcap", 1);

  draw_Rgamma(dataset, "sigmaEta_vs_pt_barrel", 1);
  draw_Rgamma(dataset, "sigmaEta_vs_pt_endcap", 1);
  draw_Rgamma(dataset, "sigmaPhi_vs_pt_barrel", 1);
  draw_Rgamma(dataset, "sigmaPhi_vs_pt_endcap", 1);

  draw_variable(dataset, "Rch_PFItCone5", drawNormalized, setLogY);
  draw_variable(dataset, "Rgamma_PFItCone5", drawNormalized);
  draw_variable(dataset, "eta_Rch0_PFItCone5", drawNormalized);
  draw_variable2D(dataset, "Rch_vs_eta_PFItCone5");
  draw_variable2D(dataset, "Rgamma_vs_Rch_PFItCone5");
  draw_variable2D(dataset, "response_vs_Rch_MEAN_PFItCone5");
  draw_variable2D(dataset, "resolution_vs_Rch_RMS_PFItCone5");
  draw_variable2D(dataset, "response_vs_Rgamma_MEAN_PFItCone5");
  draw_variable2D(dataset, "resolution_vs_Rgamma_RMS_PFItCone5");



} //drawPlots


void draw(const std::string& dataset, const std::string& name, int legendQuadrant, bool normalized) {

  std::string infileName;
  if( dataset=="") infileName = "finalOutputFile.root";
  else infileName = "finalOutputFile_"+dataset+".root";

  if( DEBUG_ ) infileName = "provaOutputFile.root";
std::cout << "Opening " << infileName << std::endl;

  TFile* file = TFile::Open(infileName.c_str());


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

  std::string histoName;

  histoName = name+"_caloItCone5";
  TH1F* h1_itCone5 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_caloSISCone5";
  TH1F* h1_sisCone5 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_calokt4";
  TH1F* h1_kt4 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5";
  TH1F* h1_pFlow = (TH1F*)file->Get(histoName.c_str());

  h1_itCone5->SetMarkerStyle(21);
  h1_itCone5->SetMarkerSize(1.5);
  h1_itCone5->SetMarkerColor(kRed);
  h1_itCone5->SetTitle("");
  h1_itCone5->GetXaxis()->SetTitleOffset(1.1);
  h1_itCone5->GetYaxis()->SetTitleOffset(1.2);
  //h1_itCone5->GetYaxis()->SetTitleSize(0.06);
  //h1_itCone5->GetXaxis()->SetTitleSize(0.06);
  //h1_itCone5->GetYaxis()->SetLabelSize(0.05);
  //h1_itCone5->GetXaxis()->SetLabelSize(0.05);

  h1_sisCone5->SetMarkerStyle(22);
  h1_sisCone5->SetMarkerSize(1.5);
  h1_sisCone5->SetMarkerColor(kGreen);
  h1_sisCone5->SetTitle("");
  h1_sisCone5->GetYaxis()->SetTitleOffset(1.2);
  //h1_sisCone5->GetYaxis()->SetTitleSize(0.06);
  //h1_sisCone5->GetXaxis()->SetTitleSize(0.06);
  //h1_sisCone5->GetYaxis()->SetLabelSize(0.05);
  //h1_sisCone5->GetXaxis()->SetLabelSize(0.05);

  h1_kt4->SetMarkerStyle(20);
  h1_kt4->SetMarkerSize(1.5);
  h1_kt4->SetMarkerColor(kBlue);
  h1_kt4->SetTitle("");
  h1_kt4->GetYaxis()->SetTitleOffset(1.2);
  //h1_kt4->GetYaxis()->SetTitleSize(0.06);
  //h1_kt4->GetXaxis()->SetTitleSize(0.06);
  //h1_kt4->GetYaxis()->SetLabelSize(0.05);
  //h1_kt4->GetXaxis()->SetLabelSize(0.05);

  h1_pFlow->SetMarkerStyle(29);
  h1_pFlow->SetMarkerSize(1.5);
  h1_pFlow->SetMarkerColor(kMagenta);
  h1_pFlow->SetTitle("");
  h1_pFlow->GetYaxis()->SetTitleOffset(1.2);
  //h1_pFlow->GetYaxis()->SetTitleSize(0.06);
  //h1_pFlow->GetXaxis()->SetTitleSize(0.06);
  //h1_pFlow->GetYaxis()->SetLabelSize(0.05);
  //h1_pFlow->GetXaxis()->SetLabelSize(0.05);

  Int_t nBins = h1_pFlow->GetNbinsX();

  for(int iBin=1; iBin<(nBins+1); ++iBin) {
    if( h1_itCone5->GetBinContent(iBin)<0.0000001 ) h1_itCone5->SetBinContent(iBin, -1.);
    if( h1_sisCone5->GetBinContent(iBin)<0.0000001 ) h1_sisCone5->SetBinContent(iBin, -1.);
    if( h1_kt4->GetBinContent(iBin)<0.0000001 ) h1_kt4->SetBinContent(iBin, -1.);
    if( h1_pFlow->GetBinContent(iBin)<0.0000001 ) h1_pFlow->SetBinContent(iBin, -1.);
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
  legend->AddEntry(h1_itCone5, "Iterative Cone 0.5", "P");
  legend->AddEntry(h1_sisCone5, "SIS Cone 0.5", "P");
  legend->AddEntry(h1_kt4, "k_{T} D=0.4", "P");
  legend->AddEntry(h1_pFlow, "PF It. Cone 0.5", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  c1->cd();

  if( normalized ) {
    h1_pFlow->DrawNormalized("P");
    h1_itCone5->DrawNormalized("Psame");
    h1_sisCone5->DrawNormalized("Psame");
    h1_kt4->DrawNormalized("Psame");
  } else {
    h1_itCone5->Draw("P");
    h1_sisCone5->Draw("Psame");
    h1_kt4->Draw("Psame");
    h1_pFlow->Draw("Psame");
  } //else

  legend->Draw("same");

  std::string plotName;
  if(dataset=="") plotName = "Plots/"+name+".eps";
  else plotName = "Plots_"+dataset+"/"+name+".eps";
  c1->SaveAs(plotName.c_str());


  delete c1;

} //draw


void draw_Rch(const std::string& dataset, const std::string& name, int legendQuadrant, bool normalized) {

  std::string infileName;
  if( dataset=="") infileName = "finalOutputFile.root";
  else infileName = "finalOutputFile_"+dataset+".root";

  if( DEBUG_ ) infileName = "provaOutputFile.root";
  TFile* file = TFile::Open(infileName.c_str());

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

  std::string histoName;

  histoName = name+"_caloItCone5";
  TH1F* h1_itCone5 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rch0";
  TH1F* h1_pFlow_Rch0 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rch10_30";
  TH1F* h1_pFlow_Rch10_30 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rch50_70";
  TH1F* h1_pFlow_Rch50_70 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rch90_100";
  TH1F* h1_pFlow_Rch90_100 = (TH1F*)file->Get(histoName.c_str());

  h1_itCone5->SetLineColor(kBlack);
  h1_itCone5->SetLineWidth(2);
  h1_itCone5->SetTitle("");
  h1_itCone5->SetMarkerSize(0);
  h1_itCone5->GetXaxis()->SetTitleOffset(1.1);
  h1_itCone5->GetYaxis()->SetTitleOffset(1.2);
  //h1_itCone5->GetYaxis()->SetRangeUser(0., 1.);
  //h1_itCone5->GetYaxis()->SetTitleSize(0.06);
  //h1_itCone5->GetXaxis()->SetTitleSize(0.06);
  //h1_itCone5->GetYaxis()->SetLabelSize(0.05);
  //h1_itCone5->GetXaxis()->SetLabelSize(0.05);

  h1_pFlow_Rch0->SetMarkerStyle(20);
  h1_pFlow_Rch0->SetMarkerSize(1.2);
  h1_pFlow_Rch0->SetMarkerColor(kGreen);
  h1_pFlow_Rch0->SetTitle("");
  h1_pFlow_Rch0->GetYaxis()->SetTitleOffset(1.2);
  //h1_sisCone5->GetYaxis()->SetTitleSize(0.06);
  //h1_sisCone5->GetXaxis()->SetTitleSize(0.06);
  //h1_sisCone5->GetYaxis()->SetLabelSize(0.05);
  //h1_sisCone5->GetXaxis()->SetLabelSize(0.05);

  h1_pFlow_Rch10_30->SetMarkerStyle(21);
  h1_pFlow_Rch10_30->SetMarkerSize(1.2);
  h1_pFlow_Rch10_30->SetMarkerColor(kRed);
  h1_pFlow_Rch10_30->SetTitle("");
  h1_pFlow_Rch10_30->GetYaxis()->SetTitleOffset(1.2);

  h1_pFlow_Rch50_70->SetMarkerStyle(22);
  h1_pFlow_Rch50_70->SetMarkerSize(1.5);
  h1_pFlow_Rch50_70->SetMarkerColor(kMagenta);
  h1_pFlow_Rch50_70->SetTitle("");
  h1_pFlow_Rch50_70->GetYaxis()->SetTitleOffset(1.2);
  //h1_kt4->GetYaxis()->SetTitleSize(0.06);
  //h1_kt4->GetXaxis()->SetTitleSize(0.06);
  //h1_kt4->GetYaxis()->SetLabelSize(0.05);
  //h1_kt4->GetXaxis()->SetLabelSize(0.05);

  h1_pFlow_Rch90_100->SetMarkerStyle(29);
  h1_pFlow_Rch90_100->SetMarkerSize(1.5);
  h1_pFlow_Rch90_100->SetMarkerColor(kBlue);
  h1_pFlow_Rch90_100->SetTitle("");
  h1_pFlow_Rch90_100->GetYaxis()->SetTitleOffset(1.2);
  //h1_pFlow->GetYaxis()->SetTitleSize(0.06);
  //h1_pFlow->GetXaxis()->SetTitleSize(0.06);
  //h1_pFlow->GetYaxis()->SetLabelSize(0.05);
  //h1_pFlow->GetXaxis()->SetLabelSize(0.05);

  Int_t nBins = h1_itCone5->GetNbinsX();

  for(int iBin=1; iBin<(nBins+1); ++iBin) {
    if( h1_itCone5->GetBinContent(iBin)<0.0000001 ) h1_itCone5->SetBinContent(iBin, -1.);
    if( h1_pFlow_Rch0->GetBinContent(iBin)<0.0000001 ) h1_pFlow_Rch0->SetBinContent(iBin, -1.);
    if( h1_pFlow_Rch10_30->GetBinContent(iBin)<0.0000001 ) h1_pFlow_Rch10_30->SetBinContent(iBin, -1.);
    if( h1_pFlow_Rch50_70->GetBinContent(iBin)<0.0000001 ) h1_pFlow_Rch50_70->SetBinContent(iBin, -1.);
    if( h1_pFlow_Rch90_100->GetBinContent(iBin)<0.0000001 ) h1_pFlow_Rch90_100->SetBinContent(iBin, -1.);
  }

  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.72, 0.85);
  if( legendQuadrant==1 )  legend = new TLegend(0.44, 0.52, 0.88, 0.85);
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.62, 0.85);
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.2, 0.6, 0.5);
  if( legendQuadrant==4 )  legend = new TLegend(0.44, 0.2, 0.88, 0.5);

  legend->SetTextSize(0.05);
  legend->AddEntry(h1_itCone5, "CALO It. Cone 0.5", "L");
  legend->AddEntry(h1_pFlow_Rch0, "PFlow, R_{ch}=0%", "P");
  legend->AddEntry(h1_pFlow_Rch10_30, "PFlow, 10<R_{ch}<30%", "P");
  legend->AddEntry(h1_pFlow_Rch50_70, "PFlow, 50<R_{ch}<70%", "P");
  legend->AddEntry(h1_pFlow_Rch90_100, "PFlow, R_{ch}>90%", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  c1->cd();

  if( normalized ) {
    
    h1_pFlow_Rch0->DrawNormalized("P");
    h1_pFlow_Rch10_30->DrawNormalized("Psame");
    h1_pFlow_Rch50_70->DrawNormalized("Psame");
    h1_pFlow_Rch90_100->DrawNormalized("Psame");
    h1_itCone5->DrawNormalized("same");
  } else {
    h1_itCone5->Draw("L");
    h1_pFlow_Rch0->Draw("Psame");
    h1_pFlow_Rch10_30->Draw("Psame");
    h1_pFlow_Rch50_70->Draw("Psame");
    h1_pFlow_Rch90_100->Draw("Psame");
  } //else

  legend->Draw("same");

  std::string plotName;
  if(dataset=="") plotName = "Plots/Rch_"+name+".eps";
  else plotName = "Plots_"+dataset+"/Rch_"+name+".eps";
  c1->SaveAs(plotName.c_str());
 

  delete c1;

}




void draw_Rgamma(const std::string& dataset, const std::string& name, int legendQuadrant, bool normalized) {

  std::string infileName;
  if( dataset=="") infileName = "finalOutputFile.root";
  else infileName = "finalOutputFile_"+dataset+".root";

  TFile* file = TFile::Open(infileName.c_str());

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);

  std::string histoName;

  histoName = name+"_caloItCone5";
  TH1F* h1_itCone5 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rgamma0";
  TH1F* h1_pFlow_Rgamma0 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rgamma5_20";
  TH1F* h1_pFlow_Rgamma5_20 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rgamma30_60";
  TH1F* h1_pFlow_Rgamma30_60 = (TH1F*)file->Get(histoName.c_str());
  histoName = name+"_PFItCone5_Rgamma60_100";
  TH1F* h1_pFlow_Rgamma60_100 = (TH1F*)file->Get(histoName.c_str());

  h1_itCone5->SetLineColor(kBlack);
  h1_itCone5->SetLineWidth(2);
  h1_itCone5->SetTitle("");
  h1_itCone5->SetMarkerSize(0);
  h1_itCone5->GetXaxis()->SetTitleOffset(1.1);
  h1_itCone5->GetYaxis()->SetTitleOffset(1.2);

  h1_pFlow_Rgamma0->SetMarkerStyle(20);
  h1_pFlow_Rgamma0->SetMarkerSize(1.2);
  h1_pFlow_Rgamma0->SetMarkerColor(kGreen);
  h1_pFlow_Rgamma0->SetTitle("");
  h1_pFlow_Rgamma0->GetYaxis()->SetTitleOffset(1.2);

  h1_pFlow_Rgamma5_20->SetMarkerStyle(21);
  h1_pFlow_Rgamma5_20->SetMarkerSize(1.2);
  h1_pFlow_Rgamma5_20->SetMarkerColor(kRed);
  h1_pFlow_Rgamma5_20->SetTitle("");
  h1_pFlow_Rgamma5_20->GetYaxis()->SetTitleOffset(1.2);

  h1_pFlow_Rgamma30_60->SetMarkerStyle(22);
  h1_pFlow_Rgamma30_60->SetMarkerSize(1.5);
  h1_pFlow_Rgamma30_60->SetMarkerColor(kMagenta);
  h1_pFlow_Rgamma30_60->SetTitle("");
  h1_pFlow_Rgamma30_60->GetYaxis()->SetTitleOffset(1.2);

  h1_pFlow_Rgamma60_100->SetMarkerStyle(29);
  h1_pFlow_Rgamma60_100->SetMarkerSize(1.5);
  h1_pFlow_Rgamma60_100->SetMarkerColor(kBlue);
  h1_pFlow_Rgamma60_100->SetTitle("");
  h1_pFlow_Rgamma60_100->GetYaxis()->SetTitleOffset(1.2);


  if( (legendQuadrant!=0)&&(legendQuadrant!=1)&&(legendQuadrant!=2)&&(legendQuadrant!=3)&&(legendQuadrant!=4) )
    legendQuadrant=1;

  TLegend* legend;
  if( legendQuadrant==0 )  legend = new TLegend(0.28, 0.52, 0.72, 0.85);
  if( legendQuadrant==1 )  legend = new TLegend(0.44, 0.52, 0.88, 0.85);
  if( legendQuadrant==2 )  legend = new TLegend(0.18, 0.52, 0.62, 0.85);
  if( legendQuadrant==3 )  legend = new TLegend(0.15, 0.2, 0.6, 0.5);
  if( legendQuadrant==4 )  legend = new TLegend(0.44, 0.2, 0.88, 0.5);

  legend->SetTextSize(0.05);
  legend->AddEntry(h1_itCone5, "CALO It. Cone 0.5", "L");
  legend->AddEntry(h1_pFlow_Rgamma0, "PFlow, R_{#gamma}=0%", "P");
  legend->AddEntry(h1_pFlow_Rgamma5_20, "PFlow, 5<R_{#gamma}<20%", "P");
  legend->AddEntry(h1_pFlow_Rgamma30_60, "PFlow, 30<R_{#gamma}<60%", "P");
  legend->AddEntry(h1_pFlow_Rgamma60_100, "PFlow, R_{#gamma}>60%", "P");
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
 
  c1->cd();

  if( normalized ) {
    
    h1_pFlow_Rgamma0->DrawNormalized("P");
    h1_pFlow_Rgamma5_20->DrawNormalized("Psame");
    h1_pFlow_Rgamma30_60->DrawNormalized("Psame");
    h1_pFlow_Rgamma60_100->DrawNormalized("Psame");
    h1_itCone5->DrawNormalized("same");
  } else {
    h1_itCone5->Draw("L");
    h1_pFlow_Rgamma0->Draw("Psame");
    h1_pFlow_Rgamma5_20->Draw("Psame");
    h1_pFlow_Rgamma30_60->Draw("Psame");
    h1_pFlow_Rgamma60_100->Draw("Psame");
  } //else

  legend->Draw("same");

  std::string plotName;
  if(dataset=="") plotName = "Plots/Rgamma_"+name+".eps";
  else plotName = "Plots_"+dataset+"/Rgamma_"+name+".eps";
  c1->SaveAs(plotName.c_str());


  delete c1;

}



void draw_variable(const std::string& dataset, const std::string& varName, bool drawNormalized, bool setLogY) {

  std::string infileName;
  if( dataset=="") infileName = "finalOutputFile.root";
  else infileName = "finalOutputFile_"+dataset+".root";

  TFile* file = TFile::Open(infileName.c_str());

  TCanvas* c1 = new TCanvas("c", "c", 800, 600);

  TH1F* h1 = (TH1F*)file->Get(varName.c_str());
  h1->SetTitle("");
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetTitleOffset(1.1);
  h1->SetYTitle("a.u.");
  h1->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h1->SetFillColor(kRed);
  h1->SetFillStyle(3004);

  if( setLogY ) c1->SetLogy(1);

  c1->cd();

  if(drawNormalized)
    h1->DrawNormalized();
  else
    h1->Draw();

  
  std::string plotName;
  if(dataset=="") plotName = "Plots/"+varName+".eps";
  else plotName = "Plots_"+dataset+"/"+varName+".eps";
  c1->SaveAs(plotName.c_str());

  delete c1;

}


void draw_variable2D(const std::string& dataset, const std::string& varName) {

  std::string infileName;
  if( dataset=="") infileName = "finalOutputFile.root";
  else infileName = "finalOutputFile_"+dataset+".root";

  TFile* file = TFile::Open(infileName.c_str());

  TCanvas* c1 = new TCanvas("c", "c", 800, 600);

  TH2F* h2 = (TH2F*)file->Get(varName.c_str());
  h2->SetTitle("");
  h2->GetXaxis()->SetTitleOffset(1.1);
  h2->GetYaxis()->SetTitleOffset(1.1);
  h2->SetYTitle("a.u.");
  h2->SetMarkerStyle(20);
  h2->SetMarkerColor(kRed);
  h2->SetMarkerSize(1.5);

  c1->cd();
  h2->Draw("P");


  std::string plotName;
  if(dataset=="") plotName = "Plots/"+varName+".eps";
  else plotName = "Plots_"+dataset+"/"+varName+".eps";
  c1->SaveAs(plotName.c_str());

  delete c1;

}
