#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.C"


TF1* getSyst_q();
TF1* getSyst_flavor();
TF1* getResolutionCurve(std::string fileName, std::string graphName, std::string funcName);



TCanvas* drawSyst() {

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


  TF1* syst_q = getSyst_q();
  syst_q->SetLineColor(kBlack);
  syst_q->SetLineStyle(2);

  TF1* syst_flavor = getSyst_flavor();
  syst_flavor->SetLineColor(kRed);
  syst_flavor->SetLineStyle(3);

  TF1* syst_total = new TF1("syst_total", "sqrt(syst_q*syst_q + syst_flavor*syst_flavor)", 20., 1400.);
  syst_total->SetLineColor(kGray);
  syst_total->SetLineWidth(0.);
  syst_total->SetFillColor(kGray);
  syst_total->SetFillStyle(3001);

  TLegend* legend = new TLegend(0.65, 0.7, 0.88, 0.88);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->AddEntry( syst_q, "MC Model", "L");
  legend->AddEntry( syst_flavor, "Jet Flavour", "L");
  legend->AddEntry( syst_total, "Total", "F");

  TH2D* axes = new TH2D("axes", "", 10, 20., 1400., 10, 0., 20.);
  axes->GetXaxis()->SetTitleOffset(1.1);
  axes->GetYaxis()->SetTitleOffset(1.2);
  axes->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  axes->SetYTitle("Relative Uncertainty [%]");
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();
  
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogx();
  axes->Draw();
  syst_total->Draw("Chist same");
  syst_q->Draw("L same");
  syst_flavor->Draw("L same");
  legend->Draw("same");
  return c1;

}

TF1* getResolutionCurve(std::string fileName, std::string graphName, std::string funcName) {

  TFile* file = TFile::Open(fileName.c_str());
  TGraphErrors* graph = (TGraphErrors*)file->Get(graphName.c_str());
  TF1* fitFunction = new TF1(funcName.c_str(), "sqrt( [0]*[0] + [1]*[1]/sqrt(x) + [2]*[2]/x)", 20., 1400.);
  graph->Fit(fitFunction, "RQ");
  file->Close();

  return fitFunction;

}

TF1* getSyst_q() {

  bool q15 = false;
  TF1* fit_extrapReso_q = getResolutionCurve("SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_pfakt5_FIT_95.root", "gr_extrapResoRecoRel_vs_pt", "fit_extrapReso_q");
  std::string fileName; 
  if( q15 )
    fileName = "SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_pfakt5_FIT_95_15Q.root";
  else
    fileName = "SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_pfakt5_FIT_95_NOQ.root";

  TF1* fit_extrapReso_noq = getResolutionCurve(fileName, "gr_extrapResoRecoRel_vs_pt", "fit_extrapReso_noq");


  std::string syst_q_formula;
  if( q15 )
    syst_q_formula = "100.*(fit_extrapReso_q-fit_extrapReso_noq)/fit_extrapReso_q";
  else
    syst_q_formula = "100.*(fit_extrapReso_noq-fit_extrapReso_q)/fit_extrapReso_q";

  TF1* syst_q = new TF1("syst_q", syst_q_formula.c_str(), 20., 1400.);
  
  return syst_q;

}



TF1* getSyst_flavor() {

  TF1* fit_intrReso_gluon = getResolutionCurve("SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_GLUON_pfakt5_FIT_95.root", "gr_intrResoRecoRel_vs_pt", "fit_intrReso_gluon");
  TF1* fit_intrReso_quark = getResolutionCurve("SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_QUARK_pfakt5_FIT_95.root", "gr_intrResoRecoRel_vs_pt", "fit_intrReso_quark");

  TFile* file1 = TFile::Open("SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_GLUON_pfakt5_FIT_95.root");
  TGraphErrors* graph1 = (TGraphErrors*)file1->Get("gr_intrResoRecoRel_vs_pt");
  TFile* file2 = TFile::Open("SecondJet_GraphFile_PhotonJet_Summer09_MCassoc_QUARK_pfakt5_FIT_95.root");
  TGraphErrors* graph2 = (TGraphErrors*)file2->Get("gr_intrResoRecoRel_vs_pt");
  graph2->SetMarkerStyle(22);
  graph2->SetMarkerColor(kRed);

  TF1* syst_flavor_tmp = new TF1("syst_flavor_tmp", "100.*(fit_intrReso_quark-fit_intrReso_gluon)/fit_intrReso_quark", 20., 1400.);
  TF1* syst_flavor= new TF1("syst_flavor", "sqrt(syst_flavor_tmp*syst_flavor_tmp)", 20., 1400.);

/*  TH2D* axes = new TH2D("axes", "", 10, 20., 1400., 10, 0., 0.25);
  axes->GetXaxis()->SetTitleOffset(1.1);
  axes->GetYaxis()->SetTitleOffset(1.2);
  axes->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  axes->SetYTitle("Relative Uncertainty [%]");
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogx();
  axes->Draw();
  fit_intrReso_gluon->SetLineWidth(0.5);
  fit_intrReso_quark->SetLineWidth(0.5);
  fit_intrReso_gluon->Draw("same");
  fit_intrReso_quark->Draw("same");
  graph1->Draw("P same");
  graph2->Draw("P same");
  c1->SaveAs("prova.eps");
exit(1);
*/
  delete syst_flavor_tmp;

  return syst_flavor;

}


