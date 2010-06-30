#include "TCanvas.h"
#include "TH1.h"
#include "TProfile.h"
#include "THStack.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include <string>
#include <iostream>

std::string PICTURE_FORMAT = "eps";

void doSingleStack_vs_eta(const std::string& recoGen, const std::string& dataset, const std::string& varY, const std::string& ptCut="");
void doSingleStack_vs_pt(const std::string& recoGen, const std::string& dataset, const std::string& varY, const std::string& etaRegion="");



void drawStack(std::string dataset) {

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


  doSingleStack_vs_eta("RECO", dataset, "R");
  //doSingleStack_vs_eta("RECO", dataset, "R", "50_80");
  doSingleStack_vs_eta("GEN", dataset, "R");
  //doSingleStack_vs_eta("GEN", dataset, "R", "50_80");

  doSingleStack_vs_pt("RECO", dataset, "R");
  doSingleStack_vs_pt("RECO", dataset, "N");
  doSingleStack_vs_pt("RECO", dataset, "Pt");

  doSingleStack_vs_pt("RECO", dataset, "R", "barrel");
  doSingleStack_vs_pt("RECO", dataset, "N", "barrel");
  doSingleStack_vs_pt("RECO", dataset, "Pt", "barrel");

//  doSingleStack_vs_pt("GEN", dataset, "R");
  doSingleStack_vs_pt("GEN", dataset, "N");
  doSingleStack_vs_pt("GEN", dataset, "Pt");

//  doSingleStack_vs_pt("GEN", dataset, "R", "barrel");
  doSingleStack_vs_pt("GEN", dataset, "N", "barrel");
  doSingleStack_vs_pt("GEN", dataset, "Pt", "barrel");

}


void doSingleStack_vs_eta(const std::string& recoGen, const std::string& dataset, const std::string& varY, const std::string& ptCut) {


  std::string fileName;

  if( recoGen== "RECO" ) 
    fileName = "finalOutputFile_"+dataset+".root";
  else if( recoGen=="GEN" )
    fileName = "RxGenOutputFile_"+dataset+".root";
  
  //ugly fix: to be removed
  if( recoGen == "GEN" )
    fileName = "RxGenOutputFile_PhotonJet_Summer09_Pt80.root";

  TFile* inFile = TFile::Open(fileName.c_str());

  std::string histoName;

  std::string suffix;

  if( recoGen=="RECO" )
    suffix = "_vs_eta";
  else if( recoGen=="GEN" )
    suffix = "Gen_vs_eta";
  if( ptCut!="" )
    suffix = suffix + "_" + ptCut;
  if( recoGen=="RECO" )
    suffix = suffix + "_PFakt5";
  
  histoName = "Rch" + suffix;
  TH1D* h1_Rch_vs_eta = (TH1D*)inFile->Get(histoName.c_str());

  histoName = "Rgamma" + suffix;
  TH1D* h1_Rgamma_vs_eta = (TH1D*)inFile->Get(histoName.c_str());

  histoName = "Rnh" + suffix;
  TH1D* h1_Rnh_vs_eta = (TH1D*)inFile->Get(histoName.c_str());

  histoName = "Rmu" + suffix;
  TH1D* h1_Rmu_vs_eta = (TH1D*)inFile->Get(histoName.c_str());


  TH1D* h1_Re_vs_eta = 0; //Re must not be put to 0. modify TreeAnalyzer_RxGen
  TH1D* h1_Rhfhad_vs_eta = 0;
  TH1D* h1_Rhfem_vs_eta = 0;
  if( recoGen=="RECO" ) {
    histoName = "Re" + suffix;
    h1_Re_vs_eta = (TH1D*)inFile->Get(histoName.c_str());
    histoName = "Rhfhad" + suffix;
    h1_Rhfhad_vs_eta = (TH1D*)inFile->Get(histoName.c_str());
    histoName = "Rhfem" + suffix;
    h1_Rhfem_vs_eta = (TH1D*)inFile->Get(histoName.c_str());
  }

  int nBinsX_Rch = h1_Rch_vs_eta->GetNbinsX();
  int nBinsX_Rgamma = h1_Rgamma_vs_eta->GetNbinsX();
  int nBinsX_Rnh = h1_Rnh_vs_eta->GetNbinsX();
  int nBinsX_Rmu = h1_Rmu_vs_eta->GetNbinsX();
  int nBinsX_Re = (h1_Re_vs_eta!=0) ? h1_Re_vs_eta->GetNbinsX() : 0;
  int nBinsX_Rhfhad = (h1_Rhfhad_vs_eta!=0) ? h1_Rhfhad_vs_eta->GetNbinsX(): 0;
  int nBinsX_Rhfem = (h1_Rhfem_vs_eta!=0) ? h1_Rhfem_vs_eta->GetNbinsX(): 0;

  double minX_Rch = h1_Rch_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rch = h1_Rch_vs_eta->GetXaxis()->GetXmax();

  double minX_Rgamma = h1_Rgamma_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rgamma = h1_Rgamma_vs_eta->GetXaxis()->GetXmax();

  double minX_Rnh = h1_Rnh_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rnh = h1_Rnh_vs_eta->GetXaxis()->GetXmax();

  double minX_Rmu = h1_Rmu_vs_eta->GetXaxis()->GetXmin();
  double maxX_Rmu = h1_Rmu_vs_eta->GetXaxis()->GetXmax();

  double minX_Re = (h1_Re_vs_eta!=0) ? h1_Re_vs_eta->GetXaxis()->GetXmin() : 0;
  double maxX_Re = (h1_Re_vs_eta!=0) ? h1_Re_vs_eta->GetXaxis()->GetXmax() : 0;

  double minX_Rhfhad = (h1_Rhfhad_vs_eta!=0) ? h1_Rhfhad_vs_eta->GetXaxis()->GetXmin():0;
  double maxX_Rhfhad = (h1_Rhfhad_vs_eta!=0) ? h1_Rhfhad_vs_eta->GetXaxis()->GetXmax():0;

  double minX_Rhfem = (h1_Rhfem_vs_eta!=0) ? h1_Rhfem_vs_eta->GetXaxis()->GetXmin():0;
  double maxX_Rhfem = (h1_Rhfem_vs_eta!=0) ? h1_Rhfem_vs_eta->GetXaxis()->GetXmax():0;

std::cout << "nBins Rch: " <<    nBinsX_Rch << std::endl;
std::cout << "nBins Rgamma: " << nBinsX_Rgamma <<std::endl;
std::cout << "nBins Rnh: " <<    nBinsX_Rnh << std::endl;
std::cout << "nBins Re: " <<     nBinsX_Re << std::endl;
std::cout << "nBins Rmu: " <<    nBinsX_Rmu << std::endl;
std::cout << "nBins Rhfhad: " << nBinsX_Rhfhad << std::endl;
std::cout << "nBins Rhfem: " <<  nBinsX_Rhfem << std::endl;

std::cout << "min Rch: " << minX_Rch << "\tmax Rch: " << maxX_Rch << std::endl;
std::cout << "min Rgamma: " << minX_Rgamma << "\tmax Rgamma: " << maxX_Rgamma << std::endl;
std::cout << "min Rnh: " << minX_Rnh << "\tmax Rnh: " << maxX_Rnh << std::endl;
std::cout << "min Re: " << minX_Re << "\tmax Re: " << maxX_Re << std::endl;
std::cout << "min Rmu: " << minX_Rmu << "\tmax Rmu: " << maxX_Rmu << std::endl;
std::cout << "min Rhfhad: " << minX_Rhfhad << "\tmax Rhfhad: " << maxX_Rhfhad << std::endl;
std::cout << "min Rhfem: " << minX_Rhfem << "\tmax Rhfem: " << maxX_Rhfem << std::endl;

  histoName = varY+"ch_vs_eta";
  TH1D* h1_Rch_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Rch, minX_Rch, maxX_Rch);
  h1_Rch_vs_eta_stack->SetFillColor(kRed);
  h1_Rch_vs_eta_stack->SetXTitle("#eta^{GEN}");
 
  histoName = varY+"gamma_vs_eta";
  TH1D* h1_Rgamma_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Rgamma, minX_Rgamma, maxX_Rgamma);
  h1_Rgamma_vs_eta_stack->SetFillColor(kBlue);
  h1_Rgamma_vs_eta_stack->SetXTitle("#eta^{GEN}");

  histoName = varY+"nh_vs_eta";
  TH1D* h1_Rnh_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Rnh, minX_Rnh, maxX_Rnh);
  h1_Rnh_vs_eta_stack->SetFillColor(kGreen+1);
  h1_Rnh_vs_eta_stack->SetXTitle("#eta^{GEN}");

  histoName = varY+"mu_vs_eta";
  TH1D* h1_Rmu_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Rmu, minX_Rmu, maxX_Rmu);
  h1_Rmu_vs_eta_stack->SetFillColor(kOrange);
  h1_Rmu_vs_eta_stack->SetXTitle("#eta^{GEN}");

  histoName = varY+"e_vs_eta";
  TH1D* h1_Re_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Re, minX_Re, maxX_Re);
  h1_Re_vs_eta_stack->SetFillColor(kCyan);
  h1_Re_vs_eta_stack->SetXTitle("#eta^{GEN}");

  histoName = varY+"hfhad_vs_eta";
  TH1D* h1_Rhfhad_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Rhfhad, minX_Rhfhad, maxX_Rhfhad);
  h1_Rhfhad_vs_eta_stack->SetFillColor(kMagenta);
  h1_Rhfhad_vs_eta_stack->SetXTitle("#eta^{GEN}");

  histoName = varY+"hfhem_vs_eta";
  TH1D* h1_Rhfem_vs_eta_stack = new TH1D(histoName.c_str(), "", nBinsX_Rhfem, minX_Rhfem, maxX_Rhfem);
  h1_Rhfem_vs_eta_stack->SetFillColor(kBlue-5);
  h1_Rhfem_vs_eta_stack->SetXTitle("#eta^{GEN}");


  int nBins = h1_Rch_vs_eta->GetNbinsX();

  for(int i=1; i<(nBins+1); ++i) {
    h1_Rch_vs_eta_stack->SetBinContent( i, h1_Rch_vs_eta->GetBinContent(i) );
    h1_Rgamma_vs_eta_stack->SetBinContent( i, h1_Rgamma_vs_eta->GetBinContent(i) );
    h1_Rnh_vs_eta_stack->SetBinContent( i, h1_Rnh_vs_eta->GetBinContent(i) );
    h1_Rmu_vs_eta_stack->SetBinContent( i, h1_Rmu_vs_eta->GetBinContent(i) );
    if( h1_Re_vs_eta!=0 )
      h1_Re_vs_eta_stack->SetBinContent( i, h1_Re_vs_eta->GetBinContent(i) );
    if( h1_Rhfhad_vs_eta!=0 )
      h1_Rhfhad_vs_eta_stack->SetBinContent( i, h1_Rhfhad_vs_eta->GetBinContent(i) );
    if( h1_Rhfem_vs_eta!=0 )
      h1_Rhfem_vs_eta_stack->SetBinContent( i, h1_Rhfem_vs_eta->GetBinContent(i) );
  }

  THStack* stack = new THStack("stack", "");
  stack->Add(h1_Rch_vs_eta_stack);
  stack->Add(h1_Rgamma_vs_eta_stack);
  stack->Add(h1_Rnh_vs_eta_stack);
  if( recoGen=="RECO" ) {
    stack->Add(h1_Rhfhad_vs_eta_stack);
    stack->Add(h1_Rhfem_vs_eta_stack);
    stack->Add(h1_Re_vs_eta_stack);
  }
  stack->Add(h1_Rmu_vs_eta_stack);
  stack->Draw();
  std::string xTitle = "#eta^{GEN}";

  std::string yTitle;
  if( varY == "R" )
    yTitle = "% of Jet Energy";
  else if( varY == "N" )
    yTitle = "Particle Multiplicity";
  else if( varY == "Pt" ) {
    if( recoGen=="RECO" )
      yTitle = "<p_{T}^{RECO}> [GeV/c]";
    else if( recoGen=="GEN" )
      yTitle = "<p_{T}^{GEN}> [GeV/c]";
  }
  stack->GetXaxis()->SetTitle(xTitle.c_str());
  stack->GetYaxis()->SetTitle(yTitle.c_str());

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  stack->Draw();
  std::string canvasName;
  if( recoGen=="RECO" )
    canvasName = "Plots_"+dataset+"/stack_" + varY +"_vs_eta";
  else if( recoGen=="GEN" )
    canvasName = "Plots_"+dataset+"/genstack_" + varY +"_vs_eta";
  if( ptCut!="" )
    canvasName = canvasName + "_" + ptCut;
  canvasName = canvasName + "." + PICTURE_FORMAT;
  c1->SaveAs(canvasName.c_str());

  delete c1;

} //doSingleStack



void doSingleStack_vs_pt(const std::string& recoGen, const std::string& dataset, const std::string& varY, const std::string& etaRegion) {

  std::string fileName;

  if( recoGen== "RECO" )
    fileName = "finalOutputFile_"+dataset+".root";
  else if( recoGen=="GEN" )
    fileName = "RxGenOutputFile_"+dataset+".root";

  TFile* inFile = TFile::Open(fileName.c_str());

  int nBins = 20;
  Double_t Lower[nBins];
    
  for(int iLower=0; iLower<11; ++iLower)
    Lower[iLower] = iLower*10.; 
  
  Lower[11] = 120.;
  Lower[12] = 140.;
  Lower[13] = 160.;
  Lower[14] = 180.;
  Lower[15] = 200.;
  Lower[16] = 250.;
  Lower[17] = 300.;
  Lower[18] = 400.;
  Lower[19] = 500.;


  std::string histoName;

  if( recoGen=="RECO" )
    histoName = varY+"ch_vs_pt";
  else if( recoGen=="GEN" )
    histoName = varY+"chGen_vs_pt";
  if( etaRegion!="")
    histoName = histoName + "_" + etaRegion;
  if( recoGen=="RECO" )
    histoName = histoName + "_PFakt5";
  //else if( recoGen=="GEN" )
  //  histoName = histoName + "_akt5";
  TH1D* h1_Rch_vs_pt = (TH1D*)inFile->Get(histoName.c_str());

  if( recoGen=="RECO" )
    histoName = varY+"gamma_vs_pt";
  else if( recoGen=="GEN" )
    histoName = varY+"gammaGen_vs_pt";
  if( etaRegion!="")
    histoName = histoName + "_" + etaRegion;
  if( recoGen=="RECO" )
    histoName = histoName + "_PFakt5";
  //else if( recoGen=="GEN" )
  //  histoName = histoName + "_akt5";
  TH1D* h1_Rgamma_vs_pt = (TH1D*)inFile->Get(histoName.c_str());

  if( recoGen=="RECO" )
    histoName = varY+"nh_vs_pt";
  else if( recoGen=="GEN" )
    histoName = varY+"nhGen_vs_pt";
  if( etaRegion!="")
    histoName = histoName + "_" + etaRegion;
  if( recoGen=="RECO" )
    histoName = histoName + "_PFakt5";
  //else if( recoGen=="GEN" )
  //  histoName = histoName + "_akt5";
  TH1D* h1_Rnh_vs_pt = (TH1D*)inFile->Get(histoName.c_str());
/*  histoName = varY+"hf_vs_pt_PFakt5";
  TH1D* h1_Rhf_vs_pt = (TH1D*)inFile->Get(histoName.c_str());
  histoName = varY+"mu_vs_pt_PFakt5";
  TH1D* h1_Rmu_vs_pt = (TH1D*)inFile->Get(histoName.c_str());
*/
  histoName = varY+"ch_vs_pt";
  TH1D* h1_Rch_vs_pt_stack = new TH1D(histoName.c_str(), "", nBins-1, Lower);
  h1_Rch_vs_pt_stack->SetFillColor(kRed);
  h1_Rch_vs_pt_stack->SetXTitle("p_{T}^{GEN}");
 
  histoName = varY+"gamma_vs_pt";
  TH1D* h1_Rgamma_vs_pt_stack = new TH1D(histoName.c_str(), "", nBins-1, Lower);
  h1_Rgamma_vs_pt_stack->SetFillColor(kGreen);
  h1_Rgamma_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  histoName = varY+"nh_vs_pt";
  TH1D* h1_Rnh_vs_pt_stack = new TH1D(histoName.c_str(), "", nBins-1, Lower);
  h1_Rnh_vs_pt_stack->SetFillColor(kBlue);
  h1_Rnh_vs_pt_stack->SetXTitle("p_{T}^{GEN}");
/*
  histoName = varY+"hf_vs_pt";
  TH1D* h1_Rhf_vs_pt_stack = new TH1D(histoName.c_str(), "", nBins-1, Lower);
  h1_Rhf_vs_pt_stack->SetFillColor(kMagenta);
  h1_Rhf_vs_pt_stack->SetXTitle("p_{T}^{GEN}");

  histoName = varY+"mu_vs_pt";
  TH1D* h1_Rmu_vs_pt_stack = new TH1D(histoName.c_str(), "", nBins-1, Lower);
  h1_Rmu_vs_pt_stack->SetFillColor(kYellow);
  h1_Rmu_vs_pt_stack->SetXTitle("p_{T}^{GEN}");
*/
  for(int i=1; i<(nBins+1); ++i) {
    h1_Rch_vs_pt_stack->SetBinContent( i, h1_Rch_vs_pt->GetBinContent(i) );
    h1_Rgamma_vs_pt_stack->SetBinContent( i, h1_Rgamma_vs_pt->GetBinContent(i) );
    h1_Rnh_vs_pt_stack->SetBinContent( i, h1_Rnh_vs_pt->GetBinContent(i) );
 //   h1_Rhf_vs_pt_stack->SetBinContent( i, h1_Rhf_vs_pt->GetBinContent(i) );
 //   h1_Rmu_vs_pt_stack->SetBinContent( i, h1_Rmu_vs_pt->GetBinContent(i) );
  }


  //h1_Rch_vs_pt_stack->GetXaxis()->SetRangeUser(20., 500.);
  //h1_Rgamma_vs_pt_stack->GetXaxis()->SetRangeUser(20., 500.);
  //h1_Rnh_vs_pt_stack->GetXaxis()->SetRangeUser(20., 500.);

  if( varY=="N" ) {
    h1_Rch_vs_pt_stack->GetYaxis()->SetRangeUser(0., 40.);
    h1_Rgamma_vs_pt_stack->GetYaxis()->SetRangeUser(0., 40.);
    h1_Rnh_vs_pt_stack->GetYaxis()->SetRangeUser(0., 40.);
//    h1_Rmu_vs_pt_stack->GetYaxis()->SetRangeUser(0., 40.);
//    if( recoGen=="RECO" )
//      h1_Rhf_vs_pt_stack->GetYaxis()->SetRangeUser(0., 40.);
  }


  THStack* stack2 = new THStack("stack2", "");
  stack2->Add(h1_Rch_vs_pt_stack);
  stack2->Add(h1_Rgamma_vs_pt_stack);
  stack2->Add(h1_Rnh_vs_pt_stack);
  //stack2->Add(h1_Rhf_vs_pt_stack);
  //stack2->Add(h1_Rmu_vs_pt_stack);

  stack2->Draw();
  stack2->GetXaxis()->SetTitle("p_{T}^{GEN} [GeV/c]");
  std::string yTitle;
  if( varY == "R" )
    yTitle = "% of Jet Energy";
  else if( varY == "N" )
    yTitle = "Particle Multiplicity";
  else if( varY == "Pt" )
    yTitle = "<p_{T}^{GEN}> [GeV/c]";
  stack2->GetYaxis()->SetTitle(yTitle.c_str());
  stack2->GetXaxis()->SetRangeUser(20., 500.);
    

  TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  c2->cd();
  stack2->Draw();
  std::string canvasName;
  if( recoGen=="RECO" )
    canvasName = "Plots_"+dataset+"/stack_" + varY;
  else if( recoGen=="GEN" )
    canvasName = "Plots_"+dataset+"/genstack_" + varY;
  canvasName = canvasName + "_vs_pt";
  if( etaRegion!="" )
    canvasName = canvasName + "_" + etaRegion;
  canvasName = canvasName + "." + PICTURE_FORMAT;
  c2->SaveAs(canvasName.c_str());

  delete c2;

}
