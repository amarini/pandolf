#include <iostream>
#include <cstdlib>
#include "DrawBase.h"

#include "TTree.h"
#include "TFile.h"


void drawSingleRoC( DrawBase* db, TTree* tree_sig, TTree* tree_bg, const std::string& suffix, const std::string& additionalCuts );
TGraph* getRoC( TH1D* h1_signal, TH1D* h1_bg );


int main() {

  // just to set the style:
  DrawBase* db = new DrawBase("LeptonStudies");
  db->set_lumiOnRightSide();
  db->set_shapeNormalization();

  system("mkdir -p LeptonStudiesPlots");


  TFile* file_ttZ = TFile::Open("LeptonStudies_2ndLevelTree_TTZ_prova.root");
  TFile* file_tt  = TFile::Open("LeptonStudies_2ndLevelTree_TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola.root");
  TFile* file_dy  = TFile::Open("LeptonStudies_2ndLevelTree_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root");

  TTree* tree_ttZ = (TTree*)file_ttZ->Get("reducedTree");
  TTree* tree_tt  = (TTree*)file_tt->Get("reducedTree");
  TTree* tree_dy  = (TTree*)file_dy->Get("reducedTree");


  drawSingleRoC( db, tree_ttZ, tree_tt, "tt", "" );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_iso", "pfIsoEle/ptEle < 0.15" );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_passedHLT", "passedHLTEle" );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_passedHLT_iso", "passedHLTEle && pfIsoEle/ptEle < 0.15" );

  drawSingleRoC( db, tree_ttZ, tree_dy, "dy", "" );
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_iso", "pfIsoEle/ptEle < 0.15" );
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_passedHLT", "passedHLTEle" );
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_passedHLT_iso", "passedHLTEle && pfIsoEle/ptEle < 0.15" );

  return 0;

}



void drawSingleRoC( DrawBase* db, TTree* tree_signal, TTree* tree_bg, const std::string& suffix, const std::string& additionalCuts ) {

  TPaveText* label = db->get_labelSqrt();


  TH1D* h1_effP_num_susyTight_pt = new TH1D("effP_num_susyTight_pt", "", 100, 0., 1000.);
  h1_effP_num_susyTight_pt->Sumw2();
  TH1D* h1_effP_denom_susyTight_pt = new TH1D("effP_denom_susyTight_pt", "", 100, 0., 1000.);
  h1_effP_denom_susyTight_pt->Sumw2();

  TH1D* h1_effP_num_susyLoose_pt = new TH1D("effP_num_susyLoose_pt", "", 100, 0., 1000.);
  h1_effP_num_susyLoose_pt->Sumw2();
  TH1D* h1_effP_denom_susyLoose_pt = new TH1D("effP_denom_susyLoose_pt", "", 100, 0., 1000.);
  h1_effP_denom_susyLoose_pt->Sumw2();

  TH1D* h1_effP_num_passedHLT_pt = new TH1D("effP_num_passedHLT_pt", "", 100, 0., 1000.);
  h1_effP_num_passedHLT_pt->Sumw2();
  TH1D* h1_effP_denom_passedHLT_pt = new TH1D("effP_denom_passedHLT_pt", "", 100, 0., 1000.);
  h1_effP_denom_passedHLT_pt->Sumw2();


  TH1D* h1_effNP_num_susyTight_pt = new TH1D("effNP_num_susyTight_pt", "", 100, 0., 1000.);
  h1_effNP_num_susyTight_pt->Sumw2();
  TH1D* h1_effNP_denom_susyTight_pt = new TH1D("effNP_denom_susyTight_pt", "", 100, 0., 1000.);
  h1_effNP_denom_susyTight_pt->Sumw2();

  TH1D* h1_effNP_num_susyLoose_pt = new TH1D("effNP_num_susyLoose_pt", "", 100, 0., 1000.);
  h1_effNP_num_susyLoose_pt->Sumw2();
  TH1D* h1_effNP_denom_susyLoose_pt = new TH1D("effNP_denom_susyLoose_pt", "", 100, 0., 1000.);
  h1_effNP_denom_susyLoose_pt->Sumw2();

  TH1D* h1_effNP_num_passedHLT_pt = new TH1D("effNP_num_passedHLT_pt", "", 100, 0., 1000.);
  h1_effNP_num_passedHLT_pt->Sumw2();
  TH1D* h1_effNP_denom_passedHLT_pt = new TH1D("effNP_denom_passedHLT_pt", "", 100, 0., 1000.);
  h1_effNP_denom_passedHLT_pt->Sumw2();


  std::string numCondition, denomCondition;

  // susy loose
  numCondition = "isLooseSUSYElectronEle && matchedToGenEle";
  denomCondition = "matchedToGenEle";
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }

  tree_signal->Project("effP_num_susyLoose_pt", "ptEle", numCondition.c_str() );
  tree_signal->Project("effP_denom_susyLoose_pt", "ptEle", denomCondition.c_str() );


  numCondition = "isLooseSUSYElectronEle && !matchedToGenEle";
  denomCondition = "!matchedToGenEle";
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }
  tree_bg->Project("effNP_num_susyLoose_pt", "ptEle", numCondition.c_str());
  tree_bg->Project("effNP_denom_susyLoose_pt", "ptEle", denomCondition.c_str());


  // susy Tight
  numCondition = "isTightSUSYElectronEle && matchedToGenEle";
  denomCondition = "matchedToGenEle";
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }

  tree_signal->Project("effP_num_susyTight_pt", "ptEle", numCondition.c_str() );
  tree_signal->Project("effP_denom_susyTight_pt", "ptEle", denomCondition.c_str() );


  numCondition = "isTightSUSYElectronEle && !matchedToGenEle";
  denomCondition = "!matchedToGenEle";
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }
  tree_bg->Project("effNP_num_susyTight_pt", "ptEle", numCondition.c_str());
  tree_bg->Project("effNP_denom_susyTight_pt", "ptEle", denomCondition.c_str());


  // HLT only
  numCondition = "passedHLTEle && matchedToGenEle";
  denomCondition = "matchedToGenEle";
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }

  tree_signal->Project("effP_num_passedHLT_pt", "ptEle", numCondition.c_str() );
  tree_signal->Project("effP_denom_passedHLT_pt", "ptEle", denomCondition.c_str() );


  numCondition = "passedHLTEle && !matchedToGenEle";
  denomCondition = "!matchedToGenEle";
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }
  tree_bg->Project("effNP_num_passedHLT_pt", "ptEle", numCondition.c_str());
  tree_bg->Project("effNP_denom_passedHLT_pt", "ptEle", denomCondition.c_str());







  TH1D* h1_mva_prompt = new TH1D("mva_prompt", "", 100, -1., 1.0001);
  h1_mva_prompt->Sumw2();
  TH1D* h1_mva_nonprompt = new TH1D("mva_nonprompt", "", 100, -1., 1.0001);
  h1_mva_nonprompt->Sumw2();

  std::string promptCondition = "matchedToGenEle";
  std::string nonpromptCondition = "!matchedToGenEle";
  if( additionalCuts != "" ) {
    promptCondition = promptCondition + " && " + additionalCuts;
    nonpromptCondition = nonpromptCondition + " && " + additionalCuts;
  }

  tree_signal->Project("mva_prompt", "mvaidtrigEle", promptCondition.c_str());
  tree_bg->Project("mva_nonprompt", "mvaidtrigEle", nonpromptCondition.c_str());

  TGraph* gr_mvaRoC = getRoC( h1_mva_prompt, h1_mva_nonprompt );
  gr_mvaRoC->SetMarkerSize(1.6);
  gr_mvaRoC->SetMarkerStyle(20); 
  gr_mvaRoC->SetMarkerColor(46); 


  h1_mva_prompt->SetFillColor( 46 );
  h1_mva_prompt->SetLineColor( 46 );
  h1_mva_prompt->SetLineWidth( 2 );
  h1_mva_prompt->SetFillStyle( 3004 );
  h1_mva_prompt->Rebin( 2 );
  h1_mva_prompt->SetYTitle( "Normalized to Unity" );
  h1_mva_prompt->SetXTitle( "Electron ID MVA" );

  h1_mva_nonprompt->SetFillColor( 38 );
  h1_mva_nonprompt->SetLineColor( 38 );
  h1_mva_nonprompt->SetLineWidth( 2 );
  h1_mva_nonprompt->SetFillStyle( 3005 );
  h1_mva_nonprompt->Rebin( 2 );

  TLegend* legend0 = new TLegend( 0.25, 0.62, 0.6, 0.82 );
  legend0->SetTextSize(0.04);
  legend0->SetFillColor(0);
  legend0->AddEntry( h1_mva_prompt, "Prompt Electrons", "F" );
  legend0->AddEntry( h1_mva_nonprompt, "Non-Prompt Electrons", "F" );

  // draw them:
  TCanvas* c0 = new TCanvas("c0", "", 600, 600);
  c0->cd();
  h1_mva_prompt->DrawNormalized("h ");
  h1_mva_nonprompt->DrawNormalized("h same");
  label->Draw("Same");
  legend0->Draw("same");
  gPad->RedrawAxis();

  std::string canvasName = "LeptonStudiesPlots/mvaDistrib_" + suffix + ".eps";
  c0->SaveAs(canvasName.c_str());




  float effPrompt_susyLoose = (h1_effP_num_susyLoose_pt->Integral() / h1_effP_denom_susyLoose_pt->Integral());
  float effNonPrompt_susyLoose = (h1_effNP_num_susyLoose_pt->Integral() / h1_effNP_denom_susyLoose_pt->Integral());
  
  float effPrompt_susyTight = (h1_effP_num_susyTight_pt->Integral() / h1_effP_denom_susyTight_pt->Integral());
  float effNonPrompt_susyTight = (h1_effNP_num_susyTight_pt->Integral() / h1_effNP_denom_susyTight_pt->Integral());

  float effPrompt_passedHLT = (h1_effP_num_passedHLT_pt->Integral() / h1_effP_denom_passedHLT_pt->Integral());
  float effNonPrompt_passedHLT = (h1_effNP_num_passedHLT_pt->Integral() / h1_effNP_denom_passedHLT_pt->Integral());

  
  TGraph* gr_susyLoose = new TGraph(0);
  gr_susyLoose->SetPoint(0, 1.-effNonPrompt_susyLoose, effPrompt_susyLoose);
  gr_susyLoose->SetMarkerSize(1.6);
  gr_susyLoose->SetMarkerStyle(21);
  gr_susyLoose->SetMarkerColor(38);

  TGraph* gr_susyTight = new TGraph(0);
  gr_susyTight->SetPoint(0, 1.-effNonPrompt_susyTight, effPrompt_susyTight);
  gr_susyTight->SetMarkerSize(1.6);
  gr_susyTight->SetMarkerStyle(29);
  gr_susyTight->SetMarkerColor(kRed+3);

  TGraph* gr_passedHLT = new TGraph(0);
  gr_passedHLT->SetPoint(0, 1.-effNonPrompt_passedHLT, effPrompt_passedHLT);
  gr_passedHLT->SetMarkerSize(1.6);
  gr_passedHLT->SetMarkerStyle(22);
  gr_passedHLT->SetMarkerColor(kOrange+1);


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.00001, 10, 0.7, 1.00001);
  h2_axes->SetXTitle("Non-Prompt Lepton Rejection");
  h2_axes->SetYTitle("Prompt Lepton Efficiency");


  h2_axes->Draw();
  gr_mvaRoC->Draw("P same");
  gr_susyLoose->Draw("P same");
  gr_susyTight->Draw("P same");
  gr_passedHLT->Draw("P same");

  float xMin_legend = (additionalCuts=="") ? 0.2 : 0.55;
  float xMax_legend = xMin_legend + 0.3;

  TLegend* legend = new TLegend( xMin_legend, 0.2, xMax_legend, 0.5 );
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( gr_passedHLT, "Passed HLT", "P" );
  legend->AddEntry( gr_susyLoose, "SUSY Loose", "P" );
  legend->AddEntry( gr_susyTight, "SUSY Tight", "P" );
  legend->AddEntry( gr_mvaRoC, "MVA ID", "P" );
  legend->Draw("same");
  
  label->Draw("same");
  
  gPad->RedrawAxis();

  canvasName = "LeptonStudiesPlots/RoC_" + suffix + ".eps";
  c1->SaveAs(canvasName.c_str());

  delete c0;
  delete c1;
  delete legend0;
  delete legend;
  delete h2_axes;
  delete h1_effP_num_susyTight_pt;
  delete h1_effP_denom_susyTight_pt;
  delete h1_effP_num_susyLoose_pt;
  delete h1_effP_denom_susyLoose_pt;
  delete h1_effP_num_passedHLT_pt;
  delete h1_effP_denom_passedHLT_pt;
  delete h1_effNP_num_susyTight_pt;
  delete h1_effNP_denom_susyTight_pt;
  delete h1_effNP_num_susyLoose_pt;
  delete h1_effNP_denom_susyLoose_pt;
  delete h1_effNP_num_passedHLT_pt;
  delete h1_effNP_denom_passedHLT_pt;
  delete h1_mva_prompt;
  delete h1_mva_nonprompt;

}





TGraph* getRoC( TH1D* h1_signal, TH1D* h1_bg ) {

  TGraph* graph = new TGraph();
  graph->SetName("RoC");

  int nbins = h1_signal->GetNbinsX();

  for( unsigned ibin=1; ibin<nbins+1; ++ibin ) {

    float eff_signal_num   = h1_signal->Integral( ibin, nbins );
    float eff_signal_denom = h1_signal->Integral( 1, nbins );

    float eff_bg_num   = h1_bg->Integral( ibin, nbins );
    float eff_bg_denom = h1_bg->Integral( 1, nbins );

    float eff_signal = eff_signal_num/eff_signal_denom;
    float eff_bg = eff_bg_num/eff_bg_denom;

    graph->SetPoint( ibin-1, 1.-eff_bg, eff_signal );

  } //for bins

  return graph;

}
