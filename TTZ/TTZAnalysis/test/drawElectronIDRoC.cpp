#include <iostream>
#include <cstdlib>
#include "DrawBase.h"

#include "TTree.h"
#include "TFile.h"


void drawSingleRoC( DrawBase* db, TTree* tree_sig, TTree* tree_bg, const std::string& suffix, const std::string& additionalCuts, float ptMin=10., float ptMax=5000. );
TGraph* getRoC( TH1D* h1_signal, TH1D* h1_bg );


int main() {

  // just to set the style:
  DrawBase* db = new DrawBase("LeptonStudies");
  db->set_lumiOnRightSide();
  db->set_shapeNormalization();

  system("mkdir -p LeptonStudiesPlots");


  TFile* file_ttZ = TFile::Open("LeptonStudies_2ndLevelTree_TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root");
  TFile* file_tt  = TFile::Open("LeptonStudies_2ndLevelTree_TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola.root");
  TFile* file_dy  = TFile::Open("LeptonStudies_2ndLevelTree_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root");

  db->add_mcFile(file_ttZ, "ttZ", "tt+Z", kOrange+1);
  db->add_mcFile(file_tt, "tt", "Top", 38);
  db->add_mcFile(file_dy, "dy", "DY", 46);

  db->set_outputdir("LeptonStudiesPlots");
//db->drawHisto_fromTree( "reducedTree", "ptEle", "matchedToGenEle", 50, 0., 120., "ptEle_prompt", "Electron p_{T}", "GeV");

//TH1D* h1_prompt_signal;
//std::vector<TH1D*> vh1_lastHistos_prompt = db->get_lastHistos_mc();
//for(unsigned i=0; i<vh1_lastHistos_prompt.size(); ++i ) {
//  if( db->get_mcFile(i).datasetName=="ttZ" ) {
//    h1_prompt_signal = new TH1D( *(vh1_lastHistos_prompt[i]) );
//  }
//}


//db->drawHisto_fromTree( "reducedTree", "ptEle", "!matchedToGenEle", 50, 0., 120., "ptEle_nonprompt", "Electron p_{T}", "GeV");
//TH1D* h1_nonprompt_dy;
//TH1D* h1_nonprompt_tt;
//std::vector<TH1D*> vh1_lastHistos_nonprompt = db->get_lastHistos_mc();
//for(unsigned i=0; i<vh1_lastHistos_nonprompt.size(); ++i ) {
//  if( db->get_mcFile(i).datasetName=="dy" ) {
//    h1_nonprompt_dy = new TH1D( *(vh1_lastHistos_nonprompt[i]) );
//  }
//  if( db->get_mcFile(i).datasetName=="tt" ) {
//    h1_nonprompt_tt = new TH1D( *(vh1_lastHistos_nonprompt[i]) );
//  }
//}


//TCanvas* c_pt = new TCanvas("c_pt", "", 600, 600);
//c_pt->cd();

//float xMin = h1_nonprompt_tt->GetXaxis()->GetXmin();
//float xMax = h1_nonprompt_tt->GetXaxis()->GetXmax();
//float yMin = 0.;
//float yMax = h1_nonprompt_dy->GetMaximum()*1.4;

//TH2D* h2_axes_pt = new TH2D("axes_pt", "", 10, xMin, xMax, 10, yMin, yMax);
//h2_axes_pt->SetXTitle("Electron p_{T} [GeV]");
//h2_axes_pt->SetYTitle("Normalized to Unity");


//h2_axes_pt->Draw();
//h1_nonprompt_tt->Draw("same");
//h1_nonprompt_dy->Draw("same");
//h1_prompt_signal->Draw("same");

//TLegend* legend_pt = new TLegend( 0.5, 0.6, 0.88, 0.88 );
//legend_pt->SetTextSize(0.04);
//legend_pt->SetFillColor(0);
//legend_pt->AddEntry( h1_prompt_signal, "tt+Z (prompt)", "F");
//legend_pt->AddEntry( h1_nonprompt_tt, "Top (non-prompt)", "F");
//legend_pt->AddEntry( h1_nonprompt_dy, "DY (non-prompt)", "F");

//legend_pt->Draw("same");

//TPaveText* label = db->get_labelSqrt();
//label->Draw("same");

//gPad->RedrawAxis();
//
//std::string canvasName_pt = db->get_outputdir() + "/ptEle.eps";
//c_pt->SaveAs(canvasName_pt.c_str()); 



  TTree* tree_ttZ = (TTree*)file_ttZ->Get("reducedTree");
  TTree* tree_tt  = (TTree*)file_tt->Get("reducedTree");
  TTree* tree_dy  = (TTree*)file_dy->Get("reducedTree");


  drawSingleRoC( db, tree_ttZ, tree_tt, "tt", "", 10., 20. );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_iso", "pfIsoEle/ptEle < 0.15", 10., 20. );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_passedHLT", "passedHLTEle", 10., 20. );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_passedHLT_iso", "passedHLTEle && pfIsoEle/ptEle < 0.15", 10., 20. );

  drawSingleRoC( db, tree_ttZ, tree_tt, "tt", "", 20., 5000. );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_iso", "pfIsoEle/ptEle < 0.15", 20., 5000.);
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_passedHLT", "passedHLTEle", 20., 5000. );
  drawSingleRoC( db, tree_ttZ, tree_tt, "tt_passedHLT_iso", "passedHLTEle && pfIsoEle/ptEle < 0.15", 20., 5000. );

  drawSingleRoC( db, tree_ttZ, tree_dy, "dy", "", 10., 20.);
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_iso", "pfIsoEle/ptEle < 0.15", 10., 20.);
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_passedHLT", "passedHLTEle", 10., 20.);
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_passedHLT_iso", "passedHLTEle && pfIsoEle/ptEle < 0.15", 10., 20.);

  drawSingleRoC( db, tree_ttZ, tree_dy, "dy", "", 20., 5000. );
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_iso", "pfIsoEle/ptEle < 0.15", 20., 5000.);
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_passedHLT", "passedHLTEle", 20., 5000. );
  drawSingleRoC( db, tree_ttZ, tree_dy, "dy_passedHLT_iso", "passedHLTEle && pfIsoEle/ptEle < 0.15", 20., 5000. );

  return 0;

}



void drawSingleRoC( DrawBase* db, TTree* tree_signal, TTree* tree_bg, const std::string& suffix, const std::string& additionalCuts, float ptMin, float ptMax ) {

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



  char ptCondition_char[300];
  sprintf( ptCondition_char, "ptEle > %f && ptEle < %f", ptMin, ptMax );
  std::string ptCondition(ptCondition_char);


  std::string numCondition, denomCondition;

  // susy loose
  numCondition = "isLooseSUSYElectronEle && matchedToGenEle && " + ptCondition;
  denomCondition = "matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }

  tree_signal->Project("effP_num_susyLoose_pt", "ptEle", numCondition.c_str() );
  tree_signal->Project("effP_denom_susyLoose_pt", "ptEle", denomCondition.c_str() );


  numCondition = "isLooseSUSYElectronEle && !matchedToGenEle && " + ptCondition;
  denomCondition = "!matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }
  tree_bg->Project("effNP_num_susyLoose_pt", "ptEle", numCondition.c_str());
  tree_bg->Project("effNP_denom_susyLoose_pt", "ptEle", denomCondition.c_str());


  // susy Tight
  numCondition = "isTightSUSYElectronEle && matchedToGenEle && " + ptCondition;
  denomCondition = "matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }

  tree_signal->Project("effP_num_susyTight_pt", "ptEle", numCondition.c_str() );
  tree_signal->Project("effP_denom_susyTight_pt", "ptEle", denomCondition.c_str() );


  numCondition = "isTightSUSYElectronEle && !matchedToGenEle && " + ptCondition;
  denomCondition = "!matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }
  tree_bg->Project("effNP_num_susyTight_pt", "ptEle", numCondition.c_str());
  tree_bg->Project("effNP_denom_susyTight_pt", "ptEle", denomCondition.c_str());


  // HLT only
  numCondition = "passedHLTEle && matchedToGenEle && " + ptCondition;
  denomCondition = "matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }

  tree_signal->Project("effP_num_passedHLT_pt", "ptEle", numCondition.c_str() );
  tree_signal->Project("effP_denom_passedHLT_pt", "ptEle", denomCondition.c_str() );


  numCondition = "passedHLTEle && !matchedToGenEle && " + ptCondition;
  denomCondition = "!matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    numCondition = numCondition + " && " + additionalCuts;
    denomCondition = denomCondition + " && " + additionalCuts;
  }
  tree_bg->Project("effNP_num_passedHLT_pt", "ptEle", numCondition.c_str());
  tree_bg->Project("effNP_denom_passedHLT_pt", "ptEle", denomCondition.c_str());





  int nBins = (ptMax<100.) ? 250 : 500;

  TH1D* h1_mva_prompt = new TH1D("mva_prompt", "", nBins, -1., 1.0001);
  h1_mva_prompt->Sumw2();
  TH1D* h1_mva_nonprompt = new TH1D("mva_nonprompt", "", nBins, -1., 1.0001);
  h1_mva_nonprompt->Sumw2();

  std::string promptCondition = "matchedToGenEle && " + ptCondition;
  std::string nonpromptCondition = "!matchedToGenEle && " + ptCondition;
  if( additionalCuts != "" ) {
    promptCondition = promptCondition + " && " + additionalCuts;
    nonpromptCondition = nonpromptCondition + " && " + additionalCuts;
  }



  TString additionalCuts_tstr(additionalCuts);
  bool has_passedHLT = additionalCuts_tstr.Contains("passedHLTEle");

  if( has_passedHLT ) {
    tree_signal->Project("mva_prompt", "mvaidtrigEle", promptCondition.c_str() );
    tree_bg->Project("mva_nonprompt", "mvaidtrigEle", nonpromptCondition.c_str() );
  } else {
    tree_signal->Project("mva_prompt", "mvaidnontrigEle", promptCondition.c_str() );
    tree_bg->Project("mva_nonprompt", "mvaidnontrigEle", nonpromptCondition.c_str() );
  }


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

  char canvasName[1000];
  sprintf( canvasName, "LeptonStudiesPlots/mvaDistrib_%s_pt%.0f_%.0f.eps", suffix.c_str(), ptMin , ptMax);
  c0->SaveAs(canvasName);




  //std::cout << " effPrompt_susyLoose = " << h1_effP_num_susyLoose_pt->Integral() << " / " << h1_effP_denom_susyLoose_pt->Integral() << std::endl;
  //std::cout << " effNonPrompt_susyLoose = " << h1_effNP_num_susyLoose_pt->Integral() << " / " << h1_effNP_denom_susyLoose_pt->Integral() << std::endl;

  //std::cout << " effPrompt_susyTight = " << h1_effP_num_susyTight_pt->Integral() << " / " << h1_effP_denom_susyTight_pt->Integral() << std::endl;
  //std::cout << " effNonPrompt_susyTight = " << h1_effNP_num_susyTight_pt->Integral() << " / " << h1_effNP_denom_susyTight_pt->Integral() << std::endl;

  //std::cout << " effPrompt_passedHLT = " << h1_effP_num_passedHLT_pt->Integral() << " / " << h1_effP_denom_passedHLT_pt->Integral() << std::endl;
  //std::cout << " effNonPrompt_passedHLT = " << h1_effNP_num_passedHLT_pt->Integral() << " / " << h1_effNP_denom_passedHLT_pt->Integral() << std::endl;

  
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


  float xMax = (additionalCuts=="") ? 1.0001 : 0.5;
  float yMin = (ptMax<100.) ? 0.4 : 0.7;

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, yMin, 1.00001);
  h2_axes->SetXTitle("Non-Prompt Lepton Rejection");
  h2_axes->SetYTitle("Prompt Lepton Efficiency");


  h2_axes->Draw();
  float xMin_legend = (additionalCuts=="") ? 0.2 : 0.6;
  float xMax_legend = xMin_legend + 0.3;
  float yMin_legend = (additionalCuts!="" && ptMin>15.) ? 0.55 : 0.2;
  float yMax_legend = yMin_legend + 0.3;


  TLegend* legend;

  if( ptMin != 10. || ptMax != 5000. ) {

    char legendTitle[500];
    if( ptMax==5000. )
      sprintf( legendTitle, "p_{T} > %.0f GeV", ptMin );
    else
      sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );

    legend = new TLegend( xMin_legend, yMin_legend, xMax_legend, yMax_legend, legendTitle );

  } else {

    legend = new TLegend( xMin_legend, yMax_legend, xMax_legend, yMax_legend );

  }
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  if( !has_passedHLT )
    legend->AddEntry( gr_passedHLT, "Passed HLT", "P" );
  legend->AddEntry( gr_susyLoose, "SUSY Loose", "P" );
  legend->AddEntry( gr_susyTight, "SUSY Tight", "P" );
  legend->AddEntry( gr_mvaRoC, "MVA ID", "P" );
  legend->Draw("same");
  
  label->Draw("same");

  gr_mvaRoC->Draw("P same");
  gr_susyLoose->Draw("P same");
  gr_susyTight->Draw("P same");
  if( !has_passedHLT )
    gr_passedHLT->Draw("P same");

  
  gPad->RedrawAxis();

  sprintf( canvasName, "LeptonStudiesPlots/RoC_%s_pt%.0f_%.0f.eps", suffix.c_str(), ptMin , ptMax);
  c1->SaveAs(canvasName);

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
