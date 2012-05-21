#include <iostream>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TEfficiency.h"

#include "DrawBase.h"




void computeEfficiency( DrawBase* db, const std::string& varName, const std::string& axisName, TTree* treedata, TTree* treemc, float ptMin, float ptMax );
TEfficiency* getSingleEff( const std::string& varName, TTree* tree, float ptMin, float ptMax );


int main() {


  DrawBase* db = new DrawBase("eff");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(4980.);

  TFile* filedata = TFile::Open("Omog_QGStudies_Photon_Run2011_FULL.root");
  TTree* treedata = (TTree*)filedata->Get("omog");

  db->add_dataFile( filedata, "dummy" );

  TChain* treemc = new TChain("omog");
  treemc->Add("Omog_QGStudies_G_Summer11.root/omog");
  treemc->Add("Omog_QGStudies_QCD_EMEnriched_Summer11.root/omog");


  computeEfficiency( db, "rmsCandJet0", "Jet Candidate RMS", treedata, treemc, 50., 100.);
  computeEfficiency( db, "rmsCandJet0", "Jet Candidate RMS", treedata, treemc, 100., 150.);

  computeEfficiency( db, "betaStarJet0", "Jet #beta*", treedata, treemc, 50., 100.);
  computeEfficiency( db, "betaStarJet0", "Jet #beta*", treedata, treemc, 100., 150.);

  return 0;

}


void computeEfficiency( DrawBase* db, const std::string& varName, const std::string& axisName, TTree* treedata, TTree* treemc, float ptMin, float ptMax ) {

  //std::cout << "Pt bin: " << ptMin << " - " << ptMax << std::endl;
  //std::cout << "Efficiency for " << varName << " (MC): " << getSingleEff( varName, treemc, ptMin, ptMax ) << std::endl;
  //std::cout << "Efficiency for " << varName << " (data): " <<  getSingleEff( varName, treedata, ptMin, ptMax ) << std::endl;
  TEfficiency* effmc = getSingleEff( varName, treemc, ptMin, ptMax );
  TEfficiency* effdata = getSingleEff( varName, treedata, ptMin, ptMax );

  effdata->SetMarkerStyle(20);
  effdata->SetMarkerSize(1.6);

  effmc->SetMarkerStyle(24);
  effmc->SetMarkerSize(1.6);
  effmc->SetMarkerColor(46);
  effmc->SetLineColor(46);

  std::string yAxisTitle = axisName + " Efficiency";

  TH2D* axes = new TH2D("axes", "", 10, -4.99, 4.99, 10, 0.8, 1.);
  axes->SetXTitle("Jet Pseudorapidity");
  axes->SetYTitle(yAxisTitle.c_str());

  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );

  //TLegend* legend = new TLegend( 0.6, 0.2, 0.88, 0.4, axisName.c_str() );
  TLegend* legend = new TLegend( 0.22, 0.2, 0.45, 0.42, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( effdata, "Data", "P" );
  legend->AddEntry( effmc, "MC", "P" );

  TPaveText* label_sqrt = db->get_labelSqrt();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  axes->Draw();
  effmc->Draw("P same");
  effdata->Draw("P same");
  legend->Draw("same");
  label_sqrt->Draw("same");

  char canvasName[300];
  sprintf( canvasName, "eff_%s_pt%.0f_%.0f.eps", varName.c_str(), ptMin, ptMax );
  c1->SaveAs(canvasName);

  delete c1;
  delete axes;
  delete legend;

}


TEfficiency* getSingleEff( const std::string& varName, TTree* tree, float ptMin, float ptMax ) {

  // will compute eff vs eta:
  Double_t etaBins[7];
  etaBins[0] = -4.7;
  etaBins[1] = -3.;
  etaBins[2] = -2.5;
  etaBins[3] = 0;
  etaBins[4] = 2.5;
  etaBins[5] = 3.;
  etaBins[6] = 4.7;

  TH1D* h1_denom = new TH1D("denom", "", 6, etaBins );
  h1_denom->Sumw2();
  TH1D* h1_num = new TH1D("num", "", 6, etaBins );
  h1_num->Sumw2();


  std::string varselection;
  if( varName=="rmsCandJet0" ) {
    varselection = "( (abs(etaJet0)<3. && exp(-rmsCandJet0)<0.03) || (abs(etaJet0)>3. && abs(etaJet0)<4.7 && exp(-rmsCandJet0)<0.04) )";
  } else if( varName=="betaStarJet0" ) {
    varselection = "betaStarJet0 < 0.2 * log( nvertex - 0.67 )";
  }

  char ptselection[500];
  sprintf( ptselection, "ptJet0>%f && ptJet0<%f", ptMin, ptMax );

  char denomselection[700];
  sprintf( denomselection, "(%s)", ptselection );
  //sprintf( denomselection, "eventWeight*(%s)", ptselection );

  char numselection[700];
  sprintf( numselection, "(%s && %s)", varselection.c_str(), ptselection );
  //sprintf( numselection, "eventWeight*(%s && %s)", varselection.c_str(), ptselection );


  tree->Project( "denom", "etaJet0", denomselection);
  tree->Project( "num", "etaJet0", numselection);

  TEfficiency* eff_vs_eta = new TEfficiency( *h1_num, *h1_denom );
  eff_vs_eta->SetConfidenceLevel(0.683);
  eff_vs_eta->SetStatisticOption(TEfficiency::kBUniform);

  char effName[300];
  sprintf( effName, "%s_pt%.0f_%.0f", varName.c_str(), ptMin, ptMax );
  eff_vs_eta->SetName( effName );


  delete h1_denom;
  delete h1_num;

  return eff_vs_eta;

}
