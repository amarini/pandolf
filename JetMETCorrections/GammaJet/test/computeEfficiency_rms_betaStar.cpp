#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

#include "DrawBase.h"




void computeEfficiency_vs_eta( DrawBase* db, const std::string& varName, const std::string& axisName, TTree* treedata, TTree* treemc, float ptMin, float ptMax );
void computeEfficiency_vs_pt( DrawBase* db, const std::string& varName, const std::string& axisName, TTree* treedata, TTree* treemc, float etaMin=0., float etaMax=4.7 );
TGraphAsymmErrors* getSingleEff_vs_eta( const std::string& varName, TTree* tree, float ptMin, float ptMax, bool isData=false );
TGraphAsymmErrors* getSingleEff_vs_pt( const std::string& varName, TTree* tree, float etaMin=0., float etaMax=4.7, bool isData=false );


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


  computeEfficiency_vs_eta( db, "rmsCandJet0", "Jet Candidate RMS", treedata, treemc, 50., 100.);
  computeEfficiency_vs_eta( db, "rmsCandJet0", "Jet Candidate RMS", treedata, treemc, 100., 150.);

  computeEfficiency_vs_eta( db, "betaStarJet0", "Jet #beta*", treedata, treemc, 50., 100.);
  computeEfficiency_vs_eta( db, "betaStarJet0", "Jet #beta*", treedata, treemc, 100., 150.);

  computeEfficiency_vs_eta( db, "betaStar_and_rmsCandJet0", "#beta* + RMS", treedata, treemc, 50., 100.);
  computeEfficiency_vs_eta( db, "betaStar_and_rmsCandJet0", "#beta* + RMS", treedata, treemc, 100., 150.);

  computeEfficiency_vs_pt( db, "betaStarJet0", "Jet #beta*", treedata, treemc, 0., 4.7);
  computeEfficiency_vs_pt( db, "betaStarJet0", "Jet #beta*", treedata, treemc, 0., 4.7);

  return 0;

}


void computeEfficiency_vs_eta( DrawBase* db, const std::string& varName, const std::string& axisName, TTree* treedata, TTree* treemc, float ptMin, float ptMax ) {

  //std::cout << "Pt bin: " << ptMin << " - " << ptMax << std::endl;
  //std::cout << "Efficiency for " << varName << " (MC): " << getSingleEff_vs_eta( varName, treemc, ptMin, ptMax ) << std::endl;
  //std::cout << "Efficiency for " << varName << " (data): " <<  getSingleEff_vs_eta( varName, treedata, ptMin, ptMax ) << std::endl;
  TGraphAsymmErrors* effmc = getSingleEff_vs_eta( varName, treemc, ptMin, ptMax );
  TGraphAsymmErrors* effdata = getSingleEff_vs_eta( varName, treedata, ptMin, ptMax, true );

  effdata->SetMarkerStyle(20);
  effdata->SetMarkerSize(1.6);

  effmc->SetMarkerStyle(24);
  effmc->SetMarkerSize(1.6);
  effmc->SetMarkerColor(46);
  effmc->SetLineColor(46);

  std::string yAxisTitle = axisName + " Efficiency";

  float yMin = (varName=="betStar_and_rmsCandJet0") ? 0.7 : 0.8;

  TH2D* axes = new TH2D("axes", "", 10, -4.99, 4.99, 10, yMin, 1.);
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


void computeEfficiency_vs_pt( DrawBase* db, const std::string& varName, const std::string& axisName, TTree* treedata, TTree* treemc, float etaMin, float etaMax ) {

  //std::cout << "Pt bin: " << ptMin << " - " << ptMax << std::endl;
  //std::cout << "Efficiency for " << varName << " (MC): " << getSingleEff_vs_eta( varName, treemc, ptMin, ptMax ) << std::endl;
  //std::cout << "Efficiency for " << varName << " (data): " <<  getSingleEff_vs_eta( varName, treedata, ptMin, ptMax ) << std::endl;
  TGraphAsymmErrors* effmc = getSingleEff_vs_pt( varName, treemc, etaMin, etaMax );
  TGraphAsymmErrors* effdata = getSingleEff_vs_pt( varName, treedata, etaMin, etaMax, true );

  effdata->SetMarkerStyle(20);
  effdata->SetMarkerSize(1.6);

  effmc->SetMarkerStyle(24);
  effmc->SetMarkerSize(1.6);
  effmc->SetMarkerColor(46);
  effmc->SetLineColor(46);

  std::string yAxisTitle = axisName + " Efficiency";

  float yMin = (varName=="betStar_and_rmsCandJet0") ? 0.7 : 0.8;

  TH2D* axes = new TH2D("axes", "", 10, -4.99, 4.99, 10, yMin, 1.);
  axes->SetXTitle("Jet Pseudorapidity");
  axes->SetYTitle(yAxisTitle.c_str());

  TLegend* legend;

  if( etaMin!=0. && etaMax!=4.7 ) {
    char legendTitle[300];
    sprintf( legendTitle, "%.0f < |#eta| < %.0f GeV", etaMin, etaMax );
    legend = new TLegend( 0.22, 0.2, 0.45, 0.42, legendTitle );
  } else {
    legend = new TLegend( 0.22, 0.2, 0.45, 0.42 );
  }

  //TLegend* legend = new TLegend( 0.6, 0.2, 0.88, 0.4, axisName.c_str() );
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
  sprintf( canvasName, "eff_%s_eta%.0f_%.0f.eps", varName.c_str(), etaMin, etaMax );
  c1->SaveAs(canvasName);

  delete c1;
  delete axes;
  delete legend;

}


TGraphAsymmErrors* getSingleEff_vs_eta( const std::string& varName, TTree* tree, float ptMin, float ptMax, bool isData ) {

  // will compute eff vs eta:
  Double_t etaBins[7];
  etaBins[0] = -4.7;
  etaBins[1] = -3.;
  etaBins[2] = -2.5;
  etaBins[3] = 0.;
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
  } else if( varName=="betaStar_and_rmsCandJet0" ) {
    varselection = "( (abs(etaJet0)<2.5 && exp(-rmsCandJet0)<0.03 && betaStarJet0<0.2*log( nvertex-0.67 ) ) || ( abs(etaJet0)>2.5 && abs(etaJet0)<3. && exp(-rmsCandJet0)<0.03) || (abs(etaJet0)>3. && abs(etaJet0)<4.7 && exp(-rmsCandJet0)<0.04) )";
  }

  char ptselection[500];
  sprintf( ptselection, "ptJet0>%f && ptJet0<%f", ptMin, ptMax );

  

  char denomselection[700];
  //sprintf( denomselection, "(%s)", ptselection );
  sprintf( denomselection, "eventWeight*(%s)", ptselection );

  char numselection[700];
  //sprintf( numselection, "(%s && %s)", varselection.c_str(), ptselection );
  sprintf( numselection, "eventWeight*(%s && %s)", varselection.c_str(), ptselection );


  tree->Project( "denom", "etaJet0", denomselection);
  tree->Project( "num", "etaJet0", numselection);

  //// set bin content to integers:
  //for( unsigned iBin=1; iBin<h1_denom->GetNbinsX()+1; ++iBin ) {
  //  h1_denom->SetBinContent(iBin, (int)h1_denom->GetBinContent(iBin) );
  //  h1_num->SetBinContent(iBin, (int)h1_num->GetBinContent(iBin) );
  //}

  //TFile* prova = TFile::Open("prova.root", "recreate");
  //prova->cd();
  //h1_denom->Write();
  //h1_num->Write();
  //prova->Close();
  //exit(1);

  //TEfficiency* eff_vs_eta = new TEfficiency( *h1_num, *h1_denom );
  //eff_vs_eta->SetConfidenceLevel(0.683);
  //eff_vs_eta->SetStatisticOption(TEfficiency::kBUniform);

  char effName[300];
  sprintf( effName, "%s_pt%.0f_%.0f", varName.c_str(), ptMin, ptMax );
  //eff_vs_eta->SetName( effName );

  TGraphAsymmErrors* gr_eff = new TGraphAsymmErrors();
  gr_eff->SetName( effName );

  if( isData ) {

    gr_eff->Divide( h1_num, h1_denom );

  } else {

    for(unsigned ibin=1; ibin<h1_num->GetNbinsX()+1; ++ibin ) {
 
      float eff = h1_num->GetBinContent(ibin)/h1_denom->GetBinContent(ibin);
      float eff_err = sqrt( eff*(1.-eff) / h1_denom->GetBinContent(ibin) );

      float x = h1_num->GetBinCenter( ibin );
      float xlo = h1_num->GetBinLowEdge( ibin );
      float xup = h1_num->GetBinLowEdge( ibin+1 );

      int thisPoint = gr_eff->GetN();
      gr_eff->SetPoint( thisPoint, x, eff );
      gr_eff->SetPointError( thisPoint, x-xlo, xup-x, eff_err, eff_err );

    }

  }


  delete h1_denom;
  delete h1_num;

  return gr_eff;

}


TGraphAsymmErrors* getSingleEff_vs_pt( const std::string& varName, TTree* tree, float etaMin, float etaMax, bool isData ) {

  // will compute eff vs pt:
  Double_t ptBins[4];
  ptBins[0] = 30.;
  ptBins[1] = 50.;
  ptBins[2] = 100.;
  ptBins[3] = 150.;

  TH1D* h1_denom = new TH1D("denom", "", 3, ptBins );
  h1_denom->Sumw2();
  TH1D* h1_num = new TH1D("num", "", 3, ptBins );
  h1_num->Sumw2();


  std::string varselection;
  if( varName=="rmsCandJet0" ) {
    varselection = "( (abs(etaJet0)<3. && exp(-rmsCandJet0)<0.03) || (abs(etaJet0)>3. && abs(etaJet0)<4.7 && exp(-rmsCandJet0)<0.04) )";
  } else if( varName=="betaStarJet0" ) {
    varselection = "betaStarJet0 < 0.2 * log( nvertex - 0.67 )";
  } else if( varName=="betaStar_and_rmsCandJet0" ) {
    varselection = "( (abs(etaJet0)<2.5 && exp(-rmsCandJet0)<0.03 && betaStarJet0<0.2*log( nvertex-0.67 ) ) || ( abs(etaJet0)>2.5 && abs(etaJet0)<3. && exp(-rmsCandJet0)<0.03) || (abs(etaJet0)>3. && abs(etaJet0)<4.7 && exp(-rmsCandJet0)<0.04) )";
  }

  char etaselection[500];
  sprintf( etaselection, "abs(etaJet0)>%f && abs(etaJet0)<%f", etaMin, etaMax );

  

  char denomselection[700];
  //sprintf( denomselection, "(%s)", ptselection );
  sprintf( denomselection, "eventWeight*(%s)", etaselection );

  char numselection[700];
  //sprintf( numselection, "(%s && %s)", varselection.c_str(), ptselection );
  sprintf( numselection, "eventWeight*(%s && %s)", varselection.c_str(), etaselection );


  tree->Project( "denom", "ptJet0", denomselection);
  tree->Project( "num", "ptJet0", numselection);

  //// set bin content to integers:
  //for( unsigned iBin=1; iBin<h1_denom->GetNbinsX()+1; ++iBin ) {
  //  h1_denom->SetBinContent(iBin, (int)h1_denom->GetBinContent(iBin) );
  //  h1_num->SetBinContent(iBin, (int)h1_num->GetBinContent(iBin) );
  //}

  //TFile* prova = TFile::Open("prova.root", "recreate");
  //prova->cd();
  //h1_denom->Write();
  //h1_num->Write();
  //prova->Close();
  //exit(1);

  //TEfficiency* eff_vs_eta = new TEfficiency( *h1_num, *h1_denom );
  //eff_vs_eta->SetConfidenceLevel(0.683);
  //eff_vs_eta->SetStatisticOption(TEfficiency::kBUniform);

  char effName[300];
  sprintf( effName, "%s_eta%.0f_%.0f", varName.c_str(), etaMin, etaMax );
  //eff_vs_eta->SetName( effName );

  TGraphAsymmErrors* gr_eff = new TGraphAsymmErrors();
  gr_eff->SetName( effName );

  if( isData ) {

    gr_eff->Divide( h1_num, h1_denom );

  } else {

    for(unsigned ibin=1; ibin<h1_num->GetNbinsX()+1; ++ibin ) {
 
      float eff = h1_num->GetBinContent(ibin)/h1_denom->GetBinContent(ibin);
      float eff_err = sqrt( eff*(1.-eff) / h1_denom->GetBinContent(ibin) );

      float x = h1_num->GetBinCenter( ibin );
      float xlo = h1_num->GetBinLowEdge( ibin );
      float xup = h1_num->GetBinLowEdge( ibin+1 );

      int thisPoint = gr_eff->GetN();
      gr_eff->SetPoint( thisPoint, x, eff );
      gr_eff->SetPointError( thisPoint, x-xlo, xup-x, eff_err, eff_err );

    }

  }


  delete h1_denom;
  delete h1_num;

  return gr_eff;

}



