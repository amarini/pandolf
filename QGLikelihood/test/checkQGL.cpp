#include <cstdlib>
#include <cmath>
#include <string>
#include "DrawBase.h"
#include "QG/QGLikelihood/interface/QGLikelihoodCalculator.h"

bool Summer12=true;
std::string plotsdir = (Summer12) ? "plots_Summer12" : "plots";



void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, TTree* tree, float ptMin, float ptMax );
void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax, const std::string& labelText="" );
void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_MLP_gluon=0, TH1D* h1_MLP_quark=0, const std::string& labelText="" );



int main() {

  DrawBase* db = new DrawBase("checkQG");
  if( !Summer12 ) db->set_is7TeV(true);

  //QGLikelihoodCalculator* qglc = new QGLikelihoodCalculator("/afs/cern.ch/user/a/amarini/scratch0/CMSSW_4_2_5/src/UserCode/pandolf/QGDev/Fit/Output/Histos.root");
  QGLikelihoodCalculator* qglc;
  if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos_2012.root");
  else           qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos.root");

  TFile* file;
  //if( Summer12 ) file = TFile::Open("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_2_TREE.root");
  if( Summer12 ) file = TFile::Open("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_2_TREE.root");
  else           file = TFile::Open("/afs/cern.ch/work/a/amarini/2ndLevel/QG/QG/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_TREE.root");


  TTree* tree = (TTree*)file->Get("tree_passedEvents");

  std::string mkdircommand = "mkdir -p " + plotsdir;
  system(mkdircommand.c_str());

  drawSinglePtBin( db, qglc, tree, 20., 25. );
  drawSinglePtBin( db, qglc, tree, 25., 30. );
  drawSinglePtBin( db, qglc, tree, 30., 40. );
  drawSinglePtBin( db, qglc, tree, 40., 50. );
  drawSinglePtBin( db, qglc, tree, 50., 65. );
  drawSinglePtBin( db, qglc, tree, 65., 80. );
  drawSinglePtBin( db, qglc, tree, 80., 100. );
  drawSinglePtBin( db, qglc, tree, 150., 200. );
  drawSinglePtBin( db, qglc, tree, 200., 250. );
  drawSinglePtBin( db, qglc, tree, 300., 400. );
  drawSinglePtBin( db, qglc, tree, 500., 600. );
  drawSinglePtBin( db, qglc, tree, 800., 1000. );

  return 0;

}




void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, TTree* tree, float ptMin, float ptMax ) {

  std::cout << "-> Processing pt bin: " << ptMin << "-" << ptMax << " GeV..." << std::endl;

  bool doFwd = (ptMin<100.);


  float pt;
  tree->SetBranchAddress("ptJet0", &pt);
  float eta;
  tree->SetBranchAddress("etaJet0", &eta);
  int pdgId;
  tree->SetBranchAddress("pdgIdPartJet0", &pdgId);
  float rho;
  tree->SetBranchAddress("rhoPF", &rho);
  int nCharged;
  tree->SetBranchAddress("nChargedJet0", &nCharged);
  int nNeutral;
  tree->SetBranchAddress("nNeutralJet0", &nNeutral);
  float ptD;
  tree->SetBranchAddress("ptDJet0", &ptD);
  float ptD_QC;
  tree->SetBranchAddress("ptD_QCJet0", &ptD_QC);
  float axis2_QC;
  tree->SetBranchAddress("axis2_QCJet0", &axis2_QC);
  int nCharged_QC;
  tree->SetBranchAddress("nChg_QCJet0", &nCharged_QC);
  int nNeutral_ptCut;
  tree->SetBranchAddress("nNeutral_ptCutJet0", &nNeutral_ptCut);
  float qglMLPJet0;
  tree->SetBranchAddress("QGLMLP", &qglMLPJet0);
  float qglJet0;
  tree->SetBranchAddress("qglJet0", &qglJet0);


  TH1D* h1_qgl_old_gluon = new TH1D("qgl_old_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_old_quark = new TH1D("qgl_old_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgl_new_gluon = new TH1D("qgl_new_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_quark = new TH1D("qgl_new_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgMLP_gluon = new TH1D("qgMLP_gluon", "", 100, 0., 1.);
  TH1D* h1_qgMLP_quark = new TH1D("qgMLP_quark", "", 100, 0., 1.);


  // in the forward:
  TH1D* h1_qgl_new_F_gluon = new TH1D("qgl_new_F_gluon", "", 100, 0., 1.0001);
  TH1D* h1_qgl_new_F_quark = new TH1D("qgl_new_F_quark", "", 100, 0., 1.0001);

  TH1D* h1_qgMLP_F_gluon = new TH1D("qgMLP_F_gluon", "", 100, 0., 1.);
  TH1D* h1_qgMLP_F_quark = new TH1D("qgMLP_F_quark", "", 100, 0., 1.);


  int nentries = tree->GetEntries();

  for( unsigned int ientry=0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( pt<ptMin || pt>ptMax ) continue;
    //if( rho>22. ) continue;

    float qgl_new = qglJet0;

    if( fabs(eta)<2.5 && h1_qgl_old_gluon->GetEntries()<10000 && h1_qgl_old_quark->GetEntries()<10000 ) { //save time

      float qgl_old = qglc->computeQGLikelihoodPU( pt, rho, nCharged, nNeutral, ptD);

      if( fabs(pdgId)<5 ) {
        h1_qgl_old_quark->Fill( qgl_old );
        h1_qgl_new_quark->Fill( qgl_new );
        h1_qgMLP_quark->Fill( qglMLPJet0 );
      }
      if( pdgId==21 ) {
        h1_qgl_old_gluon->Fill( qgl_old );
        h1_qgl_new_gluon->Fill( qgl_new );
        h1_qgMLP_gluon->Fill( qglMLPJet0 );
      }

    } else if( fabs(eta)>3. ) {

      if( fabs(pdgId)<5 ) {
        h1_qgl_new_F_quark->Fill( qgl_new );
        h1_qgMLP_F_quark->Fill( qglMLPJet0 );
      }
      if( pdgId==21 ) {
        h1_qgl_new_F_gluon->Fill( qgl_new );
        h1_qgMLP_F_gluon->Fill( qglMLPJet0 );
      }

    }
    
    if( h1_qgl_old_gluon->GetEntries()>10000 
     && h1_qgl_old_quark->GetEntries()>10000 
     && ( !doFwd || (h1_qgl_new_F_quark->GetEntries()>10000
     && h1_qgl_new_F_gluon->GetEntries()>10000) ) ) break;


  }


  drawPlot( db, h1_qgl_old_gluon, h1_qgl_old_quark, "old", ptMin, ptMax, "|#eta| < 2.5" );
  drawPlot( db, h1_qgl_new_gluon, h1_qgl_new_quark, "new", ptMin, ptMax, "|#eta| < 2.5" );
  drawPlot( db, h1_qgMLP_gluon, h1_qgMLP_quark, "MLP", ptMin, ptMax, "|#eta| < 2.5" );

  drawPlot( db, h1_qgl_new_F_gluon, h1_qgl_new_F_quark, "new_F", ptMin, ptMax, "3 < |#eta| < 5" );
  drawPlot( db, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "MLP_F", ptMin, ptMax, "3 < |#eta| < 5" );

  drawRoC(db, ptMin, ptMax, "", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_old_gluon, h1_qgl_old_quark, 0, 0, "|#eta| < 2.5");
  drawRoC(db, ptMin, ptMax, "_withMLP", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_old_gluon, h1_qgl_old_quark, h1_qgMLP_gluon, h1_qgMLP_quark, "|#eta| < 2.5");
  drawRoC(db, ptMin, ptMax, "_F", h1_qgl_new_F_gluon, h1_qgl_new_F_quark, 0, 0, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "3 < |#eta| < 5");

  delete h1_qgl_old_gluon;
  delete h1_qgl_old_quark;

  delete h1_qgl_new_gluon;
  delete h1_qgl_new_quark;

  delete h1_qgMLP_gluon;
  delete h1_qgMLP_quark;

  delete h1_qgl_new_F_gluon;
  delete h1_qgl_new_F_quark;

  delete h1_qgMLP_F_gluon;
  delete h1_qgMLP_F_quark;


}





void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax, const std::string& labelText ) {


  h1_quark->Rebin(2);
  h1_gluon->Rebin(2);

  float norm_max_g = h1_gluon->GetMaximum()/h1_gluon->Integral();
  float norm_max_q = h1_quark->GetMaximum()/h1_quark->Integral();
  float hmax = (norm_max_q>norm_max_g) ? norm_max_q : norm_max_g;

  float ymax = hmax*1.2;

  bool isMLP = h1_gluon->GetXaxis()->GetXmax()<1.;

  TH2D* h2_axes = new TH2D("axes", "", 10, h1_gluon->GetXaxis()->GetXmin(), h1_gluon->GetXaxis()->GetXmax(), 10, 0., ymax);
  if( isMLP )
    h2_axes->SetXTitle("Quark-Gluon MLP Discriminator");
  else
    h2_axes->SetXTitle("Quark-Gluon Likelihood Discriminator");
  h2_axes->SetYTitle("Normalized To Unity");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_quark->SetLineColor(38);
  h1_quark->SetLineWidth(2);
  h1_quark->SetFillColor(38);
  h1_quark->SetFillStyle(3005);

  h1_gluon->SetLineColor(46);
  h1_gluon->SetLineWidth(2);
  h1_gluon->SetFillColor(46);
  h1_gluon->SetFillStyle(3004);


  h2_axes->Draw();

  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");



  float xMin_legend = (isMLP) ? 0.2 : 0.55;
  float xMax_legend = (isMLP) ? 0.5 : 0.8;

  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend( xMin_legend, 0.7, xMax_legend, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_quark, "Quark Jets", "F");
  legend->AddEntry( h1_gluon, "Gluon Jets", "F");
  legend->Draw("same");
  
  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");



  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  float xMin_label = isMLP ? 0.7 : 0.2;
  float xMax_label = isMLP ? 0.9 : 0.4;

  TPaveText* label = new TPaveText( xMin_label, 0.83, xMax_label, 0.9, "brNDC" );
  label->SetTextSize(0.04);
  label->SetFillColor(0);
  label->AddText(labelText.c_str());
  if( labelText!="" )
    label->Draw("same");


  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.eps", plotsdir.c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);
  sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.png", plotsdir.c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;

}


void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_MLP_gluon, TH1D* h1_MLP_quark, const std::string& labelText ) {


  TGraph* gr_RoC_old = new TGraph(0);
  TGraph* gr_RoC_new = new TGraph(0);
  TGraph* gr_RoC_MLP = new TGraph(0);

  int nbins = h1_new_quark->GetNbinsX();

  for( unsigned int ibin=1; ibin<nbins+1; ++ibin ) {

    float eff_q_old = -1.;
    float eff_g_old = -1.;
  
    if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
      eff_q_old = h1_old_quark->Integral( nbins-ibin, nbins )/h1_old_quark->Integral( 1, nbins );
      eff_g_old = h1_old_gluon->Integral( nbins-ibin, nbins )/h1_old_gluon->Integral( 1, nbins );
    }
  
    float eff_q_MLP = -1.;
    float eff_g_MLP = -1.;
  
    if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) { //opposite convention:
      eff_q_MLP = h1_MLP_quark->Integral( 1, ibin )/h1_MLP_quark->Integral( 1, nbins );
      eff_g_MLP = h1_MLP_gluon->Integral( 1, ibin )/h1_MLP_gluon->Integral( 1, nbins );
    }
  
    float eff_q_new = h1_new_quark->Integral( nbins-ibin, nbins )/h1_new_quark->Integral( 1, nbins );
    float eff_g_new = h1_new_gluon->Integral( nbins-ibin, nbins )/h1_new_gluon->Integral( 1, nbins );
  
    gr_RoC_new->SetPoint( ibin-1, 1.-eff_g_new, eff_q_new );

    if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
      gr_RoC_old->SetPoint( ibin-1, 1.-eff_g_old, eff_q_old );

    if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) 
      gr_RoC_MLP->SetPoint( ibin-1, 1.-eff_g_MLP, eff_q_MLP );

  }


  gr_RoC_new->SetMarkerSize(1.3);
  gr_RoC_new->SetMarkerStyle(24);
  gr_RoC_new->SetMarkerColor(kRed+3);

  if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
    gr_RoC_old->SetMarkerSize(1.3);
    gr_RoC_old->SetMarkerStyle(20);
    gr_RoC_old->SetMarkerColor(kOrange+1);
  }

  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) {
    gr_RoC_MLP->SetMarkerSize(1.3);
    gr_RoC_MLP->SetMarkerStyle(21);
    gr_RoC_MLP->SetMarkerColor(29);
  }

  TCanvas* c1 = new TCanvas("c1_roc", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes_roc", "", 10, 0., 1.0001, 10, 0., 1.0001);
  h2_axes->SetXTitle( "Gluon Jet Rejection" );
  h2_axes->SetYTitle( "Quark Jet Efficiency" );

  h2_axes->Draw();

  TLine* diag = new TLine(0., 1., 1., 0.);
  diag->Draw("same");


  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend( 0.2, 0.2, 0.45, 0.45, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  if( h1_old_quark!=0 && h1_old_gluon!=0 )
    legend->AddEntry( gr_RoC_old, "Old LD", "P");
  legend->AddEntry( gr_RoC_new, "New LD", "P");
  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 )
    legend->AddEntry( gr_RoC_MLP, "MLP", "P");
  legend->Draw("same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  TPaveText* label = new TPaveText( 0.7, 0.83, 0.9, 0.9, "brNDC" );
  label->SetTextSize(0.04);
  label->SetFillColor(0);
  label->AddText(labelText.c_str());
  if( labelText!="" )
    label->Draw("same");

  
  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) 
    gr_RoC_MLP->Draw("p same");
  if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
    gr_RoC_old->Draw("p same");
  gr_RoC_new->Draw("p same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.eps", plotsdir.c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.png", plotsdir.c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;
}
