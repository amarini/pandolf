#include <iostream>
#include <cmath>
#include <cstdlib>
#include "DrawBase.h"
#include "fitTools.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "QG/QGLikelihood/interface/QGLikelihoodCalculator.h"


bool Summer12=true;


void drawOneVariable( DrawBase* db, TTree* tree, const std::string& varName, const std::string& axisName, int nbins, float xmin, float xmax, std::string treeVar="" );
//void compareTaggers( DrawBase *db, TH1D* h1_qgl_old, TH1D* h1_qgl_newHisto, TH1D* h1_qgMLP, const std::string& name, float ptMin, float ptMax, const std::string& labelText );

void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc,  QGLikelihoodCalculator* qglc_old, TTree* tree, float ptMin, float ptMax );
void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax, const std::string& labelText, TH1D* h1_third=0, const std::string& thirdName="Charm" );
void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_MLP_gluon=0, TH1D* h1_MLP_quark=0, const std::string& labelText="", const std::string& legendName_old="Old Likelihood", const std::string& legendName_new="New Likelihood", const std::string& legendName_MLP="MLP" );

void compareTrees( DrawBase* db, TTree* tree, TTree* tree_herwig, float ptMin, float ptMax, float etaMin, float etaMax );
void compareSingleVariable( std::string varName, const std::string& axisName, int nbins, float xmin, float xmax, DrawBase* db, TTree* tree, TTree* tree_herwig, float ptMin, float ptMax, float etaMin, float etaMax, const std::string& varExpression="" );

void drawQuarkFraction_vs_pt( DrawBase* db, TTree* tree, TTree* tree_herwig, float etaMin, float etaMax );





int main() {


  TChain* tree = new TChain("reducedTree");
  
  tree->Add("/cmsrm/pc25_2/pandolf/MC/Summer12/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_finalQG_withCHS_JEC53X/QG_2ndLevelTree_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_finalQG_withCHS_JEC53X_*.root");
  

  DrawBase* db = new DrawBase("prova");

  db->set_outputdir("QGLMCPlots");
  

  //drawOneVariable( db, tree, "ptD_QCJet", "p_{T}D", 50, 0., 1.0001);
  //drawOneVariable( db, tree, "rmsCandJet", "RMSCand", 50, 0., 0.3);
  //drawOneVariable( db, tree, "axis1_QCJet", "Axis_{1}", 50, 0., 0.3);
  //drawOneVariable( db, tree, "axis2_QCJet", "Axis_{2}", 50, 0., 0.3);
  //drawOneVariable( db, tree, "pullJet",     "Pull", 50, 0., 0.03);
  //drawOneVariable( db, tree, "RJet",        "R", 50, 0., 1.0001);
  //drawOneVariable( db, tree, "nPFCandJet", "Total Multiplicity", 50, 0., 100., "nChargedJet[0]+nNeutralJet[0]");
  //drawOneVariable( db, tree, "nChargedJet", "Charged Multiplicity", 50, 0., 100.);
  //drawOneVariable( db, tree, "nNeutralJet", "Neutral Multiplicity", 50, 0., 100.);

  QGLikelihoodCalculator* qglc, *qglc_old;
  //if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/Histos_2012_NEW.root");
  //if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/ReducedHistos_prova.root");
  //if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/ReducedHisto_2012.root");
  if( Summer12 ) qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_7_patch4_testQG_HWW/src/QuarkGluonTagger/EightTeV/data/ReducedHisto_2012.root");
  else           qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos.root");

  qglc_old = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/Histos_2012_NEW.root");

  drawSinglePtBin( db, qglc, qglc_old, tree, 20., 30. );
  drawSinglePtBin( db, qglc, qglc_old, tree, 30., 40. );
  drawSinglePtBin( db, qglc, qglc_old, tree, 50., 65. );
  drawSinglePtBin( db, qglc, qglc_old, tree, 80., 100. );
  drawSinglePtBin( db, qglc, qglc_old, tree, 200., 250. );
  drawSinglePtBin( db, qglc, qglc_old, tree, 500., 600. );

  delete qglc;
  delete qglc_old;


  TChain* tree_herwig = new TChain("reducedTree");
  tree_herwig->Add("/cmsrm/pc25_2/pandolf/MC/Summer12/QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_finalQG_withCHS_JEC53X/QG_2ndLevelTree_QCD_Pt-15to3000_TuneEE3C_Flat_8TeV_herwigpp_Summer12_DR53X-PU_S10_START53_V7A-v1_finalQG_withCHS_JEC53X_*.root/reducedTree");

  drawQuarkFraction_vs_pt( db, tree, tree_herwig, 3., 5. );
  drawQuarkFraction_vs_pt( db, tree, tree_herwig, 0., 2. );

  compareTrees( db, tree, tree_herwig, 30., 40., 0., 2. );
  compareTrees( db, tree, tree_herwig, 80., 100., 0., 2. );
  compareTrees( db, tree, tree_herwig, 200., 250., 0., 2. );
  compareTrees( db, tree, tree_herwig, 500., 600., 0., 2. );

  compareTrees( db, tree, tree_herwig, 20., 30., 3., 5. );
  compareTrees( db, tree, tree_herwig, 30., 40., 3., 5. );
  compareTrees( db, tree, tree_herwig, 50., 65., 3., 5. );
  compareTrees( db, tree, tree_herwig, 80., 120., 3., 5. );


  return 0;

}




void drawOneVariable( DrawBase* db, TTree* tree, const std::string& varName, const std::string& axisName, int nbins, float xmin, float xmax, std::string treeVar ) {

  if( treeVar=="" ) treeVar = varName + "[0]";

  std::string treeSelection = "abs(etaJet[0])<2. && ptJet[0]>80. && ptJet[0]<120.";
  std::string quarkSelection = treeSelection + "&& abs(pdgIdJet[0])<4 && abs(pdgIdJet[0])>0";
  std::string gluonSelection = treeSelection + "&& pdgIdJet[0]==21";

  std::string legendTitle = "80 < p_{T} < 120 GeV";

  TH1D* h1_quark = new TH1D("quark", "", nbins, xmin, xmax );
  TH1D* h1_gluon = new TH1D("gluon", "", nbins, xmin, xmax );

  tree->Project( "quark", treeVar.c_str(), quarkSelection.c_str() );
  tree->Project( "gluon", treeVar.c_str(), gluonSelection.c_str() );

  h1_quark->SetLineWidth(2);
  h1_quark->SetLineColor(38);
  h1_quark->SetFillColor(38);
  h1_quark->SetFillStyle(3004);

  h1_gluon->SetLineWidth(2);
  h1_gluon->SetLineColor(46);
  h1_gluon->SetFillColor(46);
  h1_gluon->SetFillStyle(3005);

  float yMax_quark = h1_quark->GetMaximum()/h1_quark->Integral();
  float yMax_gluon = h1_gluon->GetMaximum()/h1_gluon->Integral();

  float yMax = (yMax_quark>yMax_gluon) ? yMax_quark : yMax_gluon;
  yMax *= 1.3;

  TPaveText* label_eta = new TPaveText(0.2, 0.85, 0.3, 0.9, "brNDC" );
  label_eta->SetTextSize(0.038);
  label_eta->SetFillColor(0);
  label_eta->AddText("|#eta| < 2");

  TPaveText* label_top = db->get_labelTop();

  TLegend* legend = new TLegend(0.58, 0.65, 0.9, 0.9, legendTitle.c_str());
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( h1_quark, "Quark Jets", "F" );
  legend->AddEntry( h1_gluon, "Gluon Jets", "F" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, 0., yMax );
  h2_axes->SetXTitle(axisName.c_str());
  h2_axes->SetYTitle("Normalized to Unity");

  h2_axes->Draw();

  legend->Draw("same");
  h1_gluon->DrawNormalized("same");
  h1_quark->DrawNormalized("same");

  label_top->Draw("same");
  label_eta->Draw("same");

  gPad->RedrawAxis();

  std::string canvasName = db->get_outputdir() + "/" + varName.c_str() + ".eps";
  c1->SaveAs(canvasName.c_str());
  std::string canvasName_str(canvasName);
  std::string command = "epstopdf " + canvasName_str;
  system( command.c_str() );

  delete c1;
  delete h2_axes;
  delete h1_quark;
  delete h1_gluon;
  delete legend;
  delete label_eta;

}



void drawSinglePtBin( DrawBase* db, QGLikelihoodCalculator* qglc, QGLikelihoodCalculator* qglc_old, TTree* tree, float ptMin, float ptMax ) {

  std::cout << "-> Processing pt bin: " << ptMin << "-" << ptMax << " GeV..." << std::endl;

  bool doFwd = (ptMin<100.);

  int njet;
  tree->SetBranchAddress("nJet", &njet);
  float pt[20];
  tree->SetBranchAddress("ptJet", pt);
  float eta[20];
  tree->SetBranchAddress("etaJet", eta);
  int pdgId[20];
  tree->SetBranchAddress("pdgIdJet", pdgId);
  float rho;
  tree->SetBranchAddress("rhoPF", &rho);
  int nCharged[20];
  tree->SetBranchAddress("nChargedJet", nCharged);
  int nNeutral[20];
  tree->SetBranchAddress("nNeutralJet", nNeutral);
  float ptD[20];
  tree->SetBranchAddress("ptDJet", ptD);
  float ptD_QC[20];
  tree->SetBranchAddress("ptD_QCJet", ptD_QC);
  float axis2_QC[20];
  tree->SetBranchAddress("axis2_QCJet", axis2_QC);
  int nCharged_QC[20];
  tree->SetBranchAddress("nChg_QCJet", nCharged_QC);
  int nNeutral_ptCut[20];
  tree->SetBranchAddress("nNeutral_ptCutJet", nNeutral_ptCut);
  float qglMLPJet[20];
  tree->SetBranchAddress("qgMLPJet", qglMLPJet);
  float qglJet[20];
  tree->SetBranchAddress("qglJet", qglJet);

  int nbins=100;

  TH1D* h1_qgl_old_gluon = new TH1D("qgl_old_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_old_quark = new TH1D("qgl_old_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_old_charm = new TH1D("qgl_old_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_old_pu = new TH1D("qgl_old_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgl_new_gluon = new TH1D("qgl_new_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_quark = new TH1D("qgl_new_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_charm = new TH1D("qgl_new_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_pu = new TH1D("qgl_new_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgl_newHisto_gluon = new TH1D("qgl_newHisto_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_quark = new TH1D("qgl_newHisto_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_charm = new TH1D("qgl_newHisto_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_pu = new TH1D("qgl_newHisto_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgMLP_gluon = new TH1D("qgMLP_gluon", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_quark = new TH1D("qgMLP_quark", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_charm = new TH1D("qgMLP_charm", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_pu = new TH1D("qgMLP_pu", "", nbins, 0., 1.);


  // in the transition:
  TH1D* h1_qgl_old_T_gluon = new TH1D("qgl_old_T_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_old_T_quark = new TH1D("qgl_old_T_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_old_T_charm = new TH1D("qgl_old_T_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_old_T_pu = new TH1D("qgl_old_T_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgl_new_T_gluon = new TH1D("qgl_new_T_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_T_quark = new TH1D("qgl_new_T_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_T_charm = new TH1D("qgl_new_T_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_T_pu = new TH1D("qgl_new_T_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgl_newHisto_T_gluon = new TH1D("qgl_newHisto_T_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_T_quark = new TH1D("qgl_newHisto_T_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_T_charm = new TH1D("qgl_newHisto_T_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_T_pu = new TH1D("qgl_newHisto_T_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgMLP_T_gluon = new TH1D("qgMLP_T_gluon", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_T_quark = new TH1D("qgMLP_T_quark", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_T_charm = new TH1D("qgMLP_T_charm", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_T_pu = new TH1D("qgMLP_T_pu", "", nbins, 0., 1.);

  // in the forward:
  TH1D* h1_qgl_new_F_gluon = new TH1D("qgl_new_F_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_F_quark = new TH1D("qgl_new_F_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_F_charm = new TH1D("qgl_new_F_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_new_F_pu = new TH1D("qgl_new_F_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgl_newHisto_F_gluon = new TH1D("qgl_newHisto_F_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_F_quark = new TH1D("qgl_newHisto_F_quark", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_F_charm = new TH1D("qgl_newHisto_F_charm", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHisto_F_pu = new TH1D("qgl_newHisto_F_pu", "", nbins, 0., 1.0001);

  TH1D* h1_qgl_newHistoNoMult_F_gluon = new TH1D("qgl_newHistoNoMult_F_gluon", "", nbins, 0., 1.0001);
  TH1D* h1_qgl_newHistoNoMult_F_quark = new TH1D("qgl_newHistoNoMult_F_quark", "", nbins, 0., 1.0001);

  TH1D* h1_qgMLP_F_gluon = new TH1D("qgMLP_F_gluon", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_F_quark = new TH1D("qgMLP_F_quark", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_F_charm = new TH1D("qgMLP_F_charm", "", nbins, 0., 1.);
  TH1D* h1_qgMLP_F_pu = new TH1D("qgMLP_F_pu", "", nbins, 0., 1.);


  int nentries = tree->GetEntries();

  for( unsigned int ientry=0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( njet==0 ) continue;

    if( pt[0]<ptMin || pt[0]>ptMax ) continue;

    if( ptD_QC[0]>0.9 ) continue; //this is to cut out anomalous (~single particle) jets

    float qgl_new = qglJet[0];

    //if( fabs(eta[0])<2. && h1_qgl_old_gluon->GetEntries()<10000 && h1_qgl_old_quark->GetEntries()<10000 ) { //save time
    if( fabs(eta[0])<2. ) {

      float qgl_old = qglc_old->computeQGLikelihoodPU( pt[0], rho, nCharged[0], nNeutral[0], ptD[0]);
      float qgl_newHisto = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);

      if( fabs(pdgId[0])<4 && fabs(pdgId[0])>0 ) {
        h1_qgl_old_quark->Fill( qgl_old );
        h1_qgl_new_quark->Fill( qgl_new );
        h1_qgl_newHisto_quark->Fill( qgl_newHisto );
        h1_qgMLP_quark->Fill( qglMLPJet[0] );
      }
      if( pdgId[0]==21 ) {
        h1_qgl_old_gluon->Fill( qgl_old );
        h1_qgl_new_gluon->Fill( qgl_new );
        h1_qgl_newHisto_gluon->Fill( qgl_newHisto );
        h1_qgMLP_gluon->Fill( qglMLPJet[0] );
      }
      if( fabs(pdgId[0])==4 ) {
        h1_qgl_old_charm->Fill( qgl_old );
        h1_qgl_new_charm->Fill( qgl_new );
        h1_qgl_newHisto_charm->Fill( qgl_newHisto );
        h1_qgMLP_charm->Fill( qglMLPJet[0] );
      }
      if( fabs(pdgId[0])==0 ) {
        h1_qgl_old_pu->Fill( qgl_old );
        h1_qgl_new_pu->Fill( qgl_new );
        h1_qgl_newHisto_pu->Fill( qgl_newHisto );
        h1_qgMLP_pu->Fill( qglMLPJet[0] );
      }


    } else if( fabs(eta[0])<2.5 ) {

      float qgl_old = qglc_old->computeQGLikelihoodPU( pt[0], rho, nCharged[0], nNeutral[0], ptD[0]);
      float qgl_newHisto = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);

      if( fabs(pdgId[0])<4 && fabs(pdgId[0])>0 ) {
        h1_qgl_old_T_quark->Fill( qgl_old );
        h1_qgl_new_T_quark->Fill( qgl_new );
        h1_qgl_newHisto_T_quark->Fill( qgl_newHisto );
        h1_qgMLP_T_quark->Fill( qglMLPJet[0] );
      }
      if( pdgId[0]==21 ) {
        h1_qgl_old_T_gluon->Fill( qgl_old );
        h1_qgl_new_T_gluon->Fill( qgl_new );
        h1_qgl_newHisto_T_gluon->Fill( qgl_newHisto );
        h1_qgMLP_T_gluon->Fill( qglMLPJet[0] );
      }
      if( fabs(pdgId[0])==4 ) {
        h1_qgl_old_T_charm->Fill( qgl_old );
        h1_qgl_new_T_charm->Fill( qgl_new );
        h1_qgl_newHisto_T_charm->Fill( qgl_newHisto );
        h1_qgMLP_T_charm->Fill( qglMLPJet[0] );
      }
      if( fabs(pdgId[0])==0 ) {
        h1_qgl_old_T_pu->Fill( qgl_old );
        h1_qgl_new_T_pu->Fill( qgl_new );
        h1_qgl_newHisto_T_pu->Fill( qgl_newHisto );
        h1_qgMLP_T_pu->Fill( qglMLPJet[0] );
      }

    } else if( fabs(eta[0])>3. ) {

      float qgl_newHisto = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);
      float qgl_newHisto_noMult = qglc->computeQGLikelihood2012( pt[0], eta[0], rho, -1, ptD_QC[0], axis2_QC[0]);

      if( fabs(pdgId[0])<4 && fabs(pdgId[0])>0 ) {
        h1_qgl_new_F_quark->Fill( qgl_new );
        h1_qgl_newHisto_F_quark->Fill( qgl_newHisto );
        h1_qgl_newHistoNoMult_F_quark->Fill( qgl_newHisto_noMult );
        h1_qgMLP_F_quark->Fill( qglMLPJet[0] );
      }
      if( pdgId[0]==21 ) {
        h1_qgl_new_F_gluon->Fill( qgl_new );
        h1_qgl_newHisto_F_gluon->Fill( qgl_newHisto );
        h1_qgl_newHistoNoMult_F_gluon->Fill( qgl_newHisto_noMult );
        h1_qgMLP_F_gluon->Fill( qglMLPJet[0] );
      }
      if( fabs(pdgId[0])==4 ) {
        h1_qgl_new_F_charm->Fill( qgl_new );
        h1_qgl_newHisto_F_charm->Fill( qgl_newHisto );
        h1_qgMLP_F_charm->Fill( qglMLPJet[0] );
      }
      if( fabs(pdgId[0])==0 ) {
        h1_qgl_new_F_pu->Fill( qgl_new );
        h1_qgl_newHisto_F_pu->Fill( qgl_newHisto );
        h1_qgMLP_F_pu->Fill( qglMLPJet[0] );
      }

    }
    
//  if( h1_qgl_old_gluon->GetEntries()>10000 
//   && h1_qgl_old_quark->GetEntries()>10000 
//   && h1_qgl_old_T_gluon->GetEntries()>10000 
//   && h1_qgl_old_T_quark->GetEntries()>10000 
//   && ( !doFwd || (h1_qgl_new_F_quark->GetEntries()>10000
//   && h1_qgl_new_F_gluon->GetEntries()>10000) ) ) break;


  }

  drawRoC(db, ptMin, ptMax, "", h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, h1_qgl_old_gluon, h1_qgl_old_quark, 0, 0, "|#eta| < 2");

  //drawRoC(db, ptMin, ptMax, "_withMLP", h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, h1_qgl_old_gluon, h1_qgl_old_quark, h1_qgMLP_gluon, h1_qgMLP_quark, "|#eta| < 2");
  //drawRoC(db, ptMin, ptMax, "_T", h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, h1_qgl_old_T_gluon, h1_qgl_old_T_quark, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "2 < |#eta| < 2.5");
  drawRoC(db, ptMin, ptMax, "_withMLP", h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, 0, 0, h1_qgMLP_gluon, h1_qgMLP_quark, "|#eta| < 2", "", "Likelihood");
  drawRoC(db, ptMin, ptMax, "_T", h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, 0, 0, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "2 < |#eta| < 2.5", "", "Likelihood");
  drawRoC(db, ptMin, ptMax, "_F", h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, 0, 0, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "3 < |#eta| < 5", "", "Likelihood");
  drawRoC(db, ptMin, ptMax, "_F_NoMult", h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, h1_qgl_newHistoNoMult_F_gluon, h1_qgl_newHistoNoMult_F_quark, 0, 0, "3 < |#eta| < 5", "Removing Multiplicity", "Likelihood");


  drawRoC(db, ptMin, ptMax, "_charm", h1_qgl_newHisto_gluon, h1_qgl_newHisto_charm, h1_qgl_old_gluon, h1_qgl_old_charm, h1_qgMLP_gluon, h1_qgMLP_charm, "|#eta| < 2");
  drawRoC(db, ptMin, ptMax, "_charm_T", h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_charm, h1_qgl_old_T_gluon, h1_qgl_old_T_charm, h1_qgMLP_T_gluon, h1_qgMLP_T_charm, "2 < |#eta| < 2.5");
  drawRoC(db, ptMin, ptMax, "_charm_F", h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_charm, 0, 0, h1_qgMLP_F_gluon, h1_qgMLP_F_charm, "3 < |#eta| < 5");

  drawRoC(db, ptMin, ptMax, "_vsHisto", h1_qgl_new_gluon, h1_qgl_new_quark, h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, 0, 0, "|#eta| < 2", "Fit-Based Likelihood", "Histo-Based Likelihood");
  drawRoC(db, ptMin, ptMax, "_T_vsHisto", h1_qgl_new_T_gluon, h1_qgl_new_T_quark, h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, 0, 0, "2 < |#eta| < 2.5", "Fit-Based Likelihood", "Histo-Based Likelihood");
  drawRoC(db, ptMin, ptMax, "_F_vsHisto", h1_qgl_new_F_gluon, h1_qgl_new_F_quark, h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, 0, 0, "3 < |#eta| < 5", "Fit-Based Likelihood", "Histo-Based Likelihood");


  std::cout << "--- UNDERFLOWS --- " << std::endl;
  std::cout << "    |eta| < 2 : " << std::endl;
  std::cout << "    OLD: " << (float)(h1_qgl_old_gluon->GetBinContent(0) + h1_qgl_old_quark->GetBinContent(0))/(h1_qgl_old_gluon->Integral(1,nbins)+h1_qgl_old_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    NEW: " << (float)(h1_qgl_new_gluon->GetBinContent(0) + h1_qgl_new_quark->GetBinContent(0))/(h1_qgl_new_gluon->Integral(1,nbins)+h1_qgl_new_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    newHisto: " << (float)(h1_qgl_newHisto_gluon->GetBinContent(0) + h1_qgl_newHisto_quark->GetBinContent(0))/(h1_qgl_newHisto_gluon->Integral(1,nbins)+h1_qgl_newHisto_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    MLP: " << (float)(h1_qgMLP_gluon->GetBinContent(0) + h1_qgMLP_quark->GetBinContent(0))/(h1_qgMLP_gluon->Integral(1,nbins)+h1_qgMLP_quark->Integral(1,nbins)) << std::endl;

  std::cout << "    2 < |eta| < 2.5 : " << std::endl;
  std::cout << "    OLD: " << (float)(h1_qgl_old_T_gluon->GetBinContent(0) + h1_qgl_old_T_quark->GetBinContent(0))/(h1_qgl_old_T_gluon->Integral(1,nbins)+h1_qgl_old_T_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    NEW: " << (float)(h1_qgl_new_T_gluon->GetBinContent(0) + h1_qgl_new_T_quark->GetBinContent(0))/(h1_qgl_new_T_gluon->Integral(1,nbins)+h1_qgl_new_T_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    newHisto: " << (float)(h1_qgl_newHisto_T_gluon->GetBinContent(0) + h1_qgl_newHisto_T_quark->GetBinContent(0))/(h1_qgl_newHisto_T_gluon->Integral(1,nbins)+h1_qgl_newHisto_T_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    MLP: " << (float)(h1_qgMLP_T_gluon->GetBinContent(0) + h1_qgMLP_T_quark->GetBinContent(0))/(h1_qgMLP_T_gluon->Integral(1,nbins)+h1_qgMLP_T_quark->Integral(1,nbins)) << std::endl;

  std::cout << "    3 < |eta| < 5 : " << std::endl;
  std::cout << "    NEW: " << (float)(h1_qgl_new_F_gluon->GetBinContent(0) + h1_qgl_new_F_quark->GetBinContent(0))/(h1_qgl_new_F_gluon->Integral(1,nbins)+h1_qgl_new_F_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    newHisto: " << (float)(h1_qgl_newHisto_F_gluon->GetBinContent(0) + h1_qgl_newHisto_F_quark->GetBinContent(0))/(h1_qgl_newHisto_F_gluon->Integral(1,nbins)+h1_qgl_newHisto_F_quark->Integral(1,nbins)) << std::endl;
  std::cout << "    MLP: " << (float)(h1_qgMLP_F_gluon->GetBinContent(0) + h1_qgMLP_F_quark->GetBinContent(0))/(h1_qgMLP_F_gluon->Integral(1,nbins)+h1_qgMLP_F_quark->Integral(1,nbins)) << std::endl;


  h1_qgl_old_gluon->Rebin(2);
  h1_qgl_old_quark->Rebin(2);
  h1_qgl_old_charm->Rebin(2);
  h1_qgl_old_pu->Rebin(2);

  h1_qgl_new_gluon->Rebin(2);
  h1_qgl_new_quark->Rebin(2);
  h1_qgl_new_charm->Rebin(2);
  h1_qgl_new_pu->Rebin(2);

  h1_qgl_newHisto_gluon->Rebin(2);
  h1_qgl_newHisto_quark->Rebin(2);
  h1_qgl_newHisto_charm->Rebin(2);
  h1_qgl_newHisto_pu->Rebin(2);

  h1_qgMLP_gluon->Rebin(2);
  h1_qgMLP_quark->Rebin(2);
  h1_qgMLP_charm->Rebin(2);
  h1_qgMLP_pu->Rebin(2);

  h1_qgl_old_T_gluon->Rebin(2);
  h1_qgl_old_T_quark->Rebin(2);
  h1_qgl_old_T_charm->Rebin(2);
  h1_qgl_old_T_pu->Rebin(2);

  h1_qgl_new_T_gluon->Rebin(2);
  h1_qgl_new_T_quark->Rebin(2);
  h1_qgl_new_T_charm->Rebin(2);
  h1_qgl_new_T_pu->Rebin(2);

  h1_qgl_newHisto_T_gluon->Rebin(2);
  h1_qgl_newHisto_T_quark->Rebin(2);
  h1_qgl_newHisto_T_charm->Rebin(2);
  h1_qgl_newHisto_T_pu->Rebin(2);

  h1_qgMLP_T_gluon->Rebin(2);
  h1_qgMLP_T_quark->Rebin(2);
  h1_qgMLP_T_charm->Rebin(2);
  h1_qgMLP_T_pu->Rebin(2);

  h1_qgl_new_F_gluon->Rebin(2);
  h1_qgl_new_F_quark->Rebin(2);
  h1_qgl_new_F_charm->Rebin(2);
  h1_qgl_new_F_pu->Rebin(2);

  h1_qgl_newHisto_F_gluon->Rebin(2);
  h1_qgl_newHisto_F_quark->Rebin(2);
  h1_qgl_newHisto_F_charm->Rebin(2);
  h1_qgl_newHisto_F_pu->Rebin(2);

  h1_qgMLP_F_gluon->Rebin(2);
  h1_qgMLP_F_quark->Rebin(2);
  h1_qgMLP_F_charm->Rebin(2);
  h1_qgMLP_F_pu->Rebin(2);


  drawPlot( db, h1_qgl_old_gluon, h1_qgl_old_quark, "old", ptMin, ptMax, "|#eta| < 2" );
  drawPlot( db, h1_qgl_new_gluon, h1_qgl_new_quark, "new", ptMin, ptMax, "|#eta| < 2" );
  drawPlot( db, h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, "newHisto", ptMin, ptMax, "|#eta| < 2" );
  drawPlot( db, h1_qgMLP_gluon, h1_qgMLP_quark, "MLP", ptMin, ptMax, "|#eta| < 2" );

  drawPlot( db, h1_qgl_old_T_gluon, h1_qgl_old_T_quark, "old_T", ptMin, ptMax, "2 < |#eta| < 2.5" );
  drawPlot( db, h1_qgl_new_T_gluon, h1_qgl_new_T_quark, "new_T", ptMin, ptMax, "2 < |#eta| < 2.5" );
  drawPlot( db, h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, "newHisto_T", ptMin, ptMax, "2 < |#eta| < 2.5" );
  drawPlot( db, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "MLP_T", ptMin, ptMax, "2 < |#eta| < 2.5" );

  drawPlot( db, h1_qgl_new_F_gluon, h1_qgl_new_F_quark, "new_F", ptMin, ptMax, "3 < |#eta| < 5" );
  drawPlot( db, h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, "newHisto_F", ptMin, ptMax, "3 < |#eta| < 5" );
  drawPlot( db, h1_qgl_newHistoNoMult_F_gluon, h1_qgl_newHistoNoMult_F_quark, "newHistoNoMult_F", ptMin, ptMax, "3 < |#eta| < 5" );
  drawPlot( db, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "MLP_F", ptMin, ptMax, "3 < |#eta| < 5" );

  drawPlot( db, h1_qgl_old_gluon, h1_qgl_old_quark, "old", ptMin, ptMax, "|#eta| < 2", h1_qgl_old_charm );
  drawPlot( db, h1_qgl_new_gluon, h1_qgl_new_quark, "new", ptMin, ptMax, "|#eta| < 2", h1_qgl_new_charm );
  drawPlot( db, h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, "newHisto", ptMin, ptMax, "|#eta| < 2", h1_qgl_newHisto_charm );
  drawPlot( db, h1_qgMLP_gluon, h1_qgMLP_quark, "MLP", ptMin, ptMax, "|#eta| < 2", h1_qgMLP_charm );

  drawPlot( db, h1_qgl_old_T_gluon, h1_qgl_old_T_quark, "old_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgl_old_T_charm );
  drawPlot( db, h1_qgl_new_T_gluon, h1_qgl_new_T_quark, "new_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgl_new_T_charm );
  drawPlot( db, h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, "newHisto_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgl_newHisto_T_charm );
  drawPlot( db, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "MLP_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgMLP_T_charm );

  drawPlot( db, h1_qgl_new_F_gluon, h1_qgl_new_F_quark, "new_F", ptMin, ptMax, "3 < |#eta| < 5", h1_qgl_new_F_charm );
  drawPlot( db, h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, "newHisto_F", ptMin, ptMax, "3 < |#eta| < 5", h1_qgl_newHisto_F_charm );
  drawPlot( db, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "MLP_F", ptMin, ptMax, "3 < |#eta| < 5", h1_qgMLP_F_charm );

  drawPlot( db, h1_qgl_old_gluon, h1_qgl_old_quark, "old", ptMin, ptMax, "|#eta| < 2", h1_qgl_old_pu, "PU" );
  drawPlot( db, h1_qgl_new_gluon, h1_qgl_new_quark, "new", ptMin, ptMax, "|#eta| < 2", h1_qgl_new_pu, "PU" );
  drawPlot( db, h1_qgl_newHisto_gluon, h1_qgl_newHisto_quark, "newHisto", ptMin, ptMax, "|#eta| < 2", h1_qgl_newHisto_pu, "PU" );
  drawPlot( db, h1_qgMLP_gluon, h1_qgMLP_quark, "MLP", ptMin, ptMax, "|#eta| < 2", h1_qgMLP_pu, "PU" );

  drawPlot( db, h1_qgl_old_T_gluon, h1_qgl_old_T_quark, "old_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgl_old_T_pu, "PU" );
  drawPlot( db, h1_qgl_new_T_gluon, h1_qgl_new_T_quark, "new_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgl_new_T_pu, "PU" );
  drawPlot( db, h1_qgl_newHisto_T_gluon, h1_qgl_newHisto_T_quark, "newHisto_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgl_newHisto_T_pu, "PU" );
  drawPlot( db, h1_qgMLP_T_gluon, h1_qgMLP_T_quark, "MLP_T", ptMin, ptMax, "2 < |#eta| < 2.5", h1_qgMLP_T_pu, "PU" );

  drawPlot( db, h1_qgl_new_F_gluon, h1_qgl_new_F_quark, "new_F", ptMin, ptMax, "3 < |#eta| < 5", h1_qgl_new_F_pu, "PU" );
  drawPlot( db, h1_qgl_newHisto_F_gluon, h1_qgl_newHisto_F_quark, "newHisto_F", ptMin, ptMax, "3 < |#eta| < 5", h1_qgl_newHisto_F_pu, "PU" );
  drawPlot( db, h1_qgMLP_F_gluon, h1_qgMLP_F_quark, "MLP_F", ptMin, ptMax, "3 < |#eta| < 5", h1_qgMLP_F_pu, "PU" );

  //compareTaggers( db, h1_qgl_old_charm, h1_qgl_newHisto_charm, h1_qgMLP_charm,       "charm", ptMin, ptMax, "|#eta| < 2" );
  //compareTaggers( db, h1_qgl_old_T_charm, h1_qgl_newHisto_T_charm, h1_qgMLP_T_charm, "charm", ptMin, ptMax, "2 < |#eta| < 2.5" );
  //compareTaggers( db, 0, h1_qgl_newHisto_F_charm, h1_qgMLP_F_charm, "charm", ptMin, ptMax, "3 < |#eta| < 5" );


  delete h1_qgl_old_gluon;
  delete h1_qgl_old_quark;
  delete h1_qgl_old_charm;
  delete h1_qgl_old_pu;

  delete h1_qgl_new_gluon;
  delete h1_qgl_new_quark;
  delete h1_qgl_new_charm;
  delete h1_qgl_new_pu;

  delete h1_qgl_newHisto_gluon;
  delete h1_qgl_newHisto_quark;
  delete h1_qgl_newHisto_charm;
  delete h1_qgl_newHisto_pu;

  delete h1_qgMLP_gluon;
  delete h1_qgMLP_quark;
  delete h1_qgMLP_charm;
  delete h1_qgMLP_pu;

  delete h1_qgl_old_T_gluon;
  delete h1_qgl_old_T_quark;
  delete h1_qgl_old_T_charm;
  delete h1_qgl_old_T_pu;

  delete h1_qgl_new_T_gluon;
  delete h1_qgl_new_T_quark;
  delete h1_qgl_new_T_charm;
  delete h1_qgl_new_T_pu;

  delete h1_qgl_newHisto_T_gluon;
  delete h1_qgl_newHisto_T_quark;
  delete h1_qgl_newHisto_T_charm;
  delete h1_qgl_newHisto_T_pu;

  delete h1_qgMLP_T_gluon;
  delete h1_qgMLP_T_quark;
  delete h1_qgMLP_T_charm;
  delete h1_qgMLP_T_pu;

  delete h1_qgl_new_F_gluon;
  delete h1_qgl_new_F_quark;
  delete h1_qgl_new_F_charm;
  delete h1_qgl_new_F_pu;

  delete h1_qgl_newHisto_F_gluon;
  delete h1_qgl_newHisto_F_quark;
  delete h1_qgl_newHisto_F_charm;
  delete h1_qgl_newHisto_F_pu;

  delete h1_qgl_newHistoNoMult_F_gluon;
  delete h1_qgl_newHistoNoMult_F_quark;

  delete h1_qgMLP_F_gluon;
  delete h1_qgMLP_F_quark;
  delete h1_qgMLP_F_charm;
  delete h1_qgMLP_F_pu;


}





void drawPlot( DrawBase* db, TH1D* h1_gluon, TH1D* h1_quark, std::string name, float ptMin, float ptMax, const std::string& labelText, TH1D* h1_third, const std::string& thirdName ) {


  float norm_max_g = h1_gluon->GetMaximum()/h1_gluon->Integral();
  float norm_max_q = h1_quark->GetMaximum()/h1_quark->Integral();
  float hmax = (norm_max_q>norm_max_g) ? norm_max_q : norm_max_g;

  float ymax = hmax*1.2;

  bool isMLP = (name=="MLP" );

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

  if( h1_third!=0 ) {
    h1_third->SetLineColor(kBlack);
    h1_third->SetLineWidth(2);
  }


  h2_axes->Draw();



  float xMin_legend = (isMLP) ? 0.2 : 0.55;
  float xMax_legend = (isMLP) ? 0.5 : 0.8;

  std::string third_legendName = thirdName + " Jets";
  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend( xMin_legend, 0.7, xMax_legend, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_quark, "Quark Jets", "F");
  legend->AddEntry( h1_gluon, "Gluon Jets", "F");
  if( h1_third!=0 ) 
    legend->AddEntry( h1_third, third_legendName.c_str(), "F");
  legend->Draw("same");
  

  h1_quark->DrawNormalized("same");
  h1_gluon->DrawNormalized("same");
  if( h1_third!=0 ) 
    h1_third->DrawNormalized("same");



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

  
  std::string thirdCanvasName = "with" + thirdName;

  char canvasName[500];
  if( h1_third!=0 ) 
    sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f_%s.eps", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax, thirdCanvasName.c_str());
  else
    sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.eps", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);
  std::string canvasName_str(canvasName);
  std::string command = "epstopdf " + canvasName_str;
  system( command.c_str() );
  if( h1_third!=0 ) 
    sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f_%s.png", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax, thirdCanvasName.c_str());
  else
    sprintf( canvasName, "%s/qgl_%s_pt%.0f_%.0f.png", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax);
  c1->SaveAs(canvasName);


  delete c1;
  delete h2_axes;
  delete legend;

}




//void compareTaggers( DrawBase *db, TH1D* h1_qgl_old, TH1D* h1_qgl_new, TH1D* h1_qgMLP, const std::string& name, float ptMin, float ptMax, const std::string& labelText ) {
//
//
//  float norm_max_mlp = h1_qgMLP->GetMaximum()/h1_qgMLP->Integral();
//  float norm_max_new = h1_qgl_new->GetMaximum()/h1_qgl_new->Integral();
//  float hmax = (norm_max_mlp>norm_max_new) ? norm_max_mlp : norm_max_new;
//
//  float ymax = hmax*1.4;
//
//
//  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., ymax);
//  h2_axes->SetXTitle("Q-G Discriminator");
//  h2_axes->SetYTitle("Normalized To Unity");
//
//  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
//  c1->cd();
//
//  if( h1_qgl_old!=0 ) {
//    h1_qgl_old->SetLineColor(38);
//    h1_qgl_old->SetLineWidth(2);
//  //h1_qgl_old->SetFillColor(38);
//  //h1_qgl_old->SetFillStyle(3005);
//  }
//
//  h1_qgl_new->SetLineColor(46);
//  h1_qgl_new->SetLineWidth(2);
//  //h1_qgl_new->SetFillColor(46);
//  //h1_qgl_new->SetFillStyle(3004);
//
//  h1_qgMLP->SetLineColor(kBlack);
//  h1_qgMLP->SetLineWidth(2);
//
//  h2_axes->Draw();
//
//  if( h1_qgl_old!=0 ) 
//    h1_qgl_old->DrawNormalized("same");
//  h1_qgl_new->DrawNormalized("same");
//  h1_qgMLP->DrawNormalized("same");
//
//
//
//  float xMin_legend = 0.55;
//  float xMax_legend = 0.8;
//
//  char legendTitle[300];
//  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
//  TLegend* legend = new TLegend( xMin_legend, 0.7, xMax_legend, 0.9, legendTitle );
//  legend->SetFillColor(0);
//  legend->SetTextSize(0.04);
//  if( h1_qgl_old!=0 ) 
//    legend->AddEntry( h1_qgl_old, "Old LD", "F");
//  legend->AddEntry( h1_qgl_new, "New LD", "F");
//  legend->AddEntry( h1_qgMLP, "MLP", "F");
//  legend->Draw("same");
//  
//
//
//
//  TPaveText* labelTop = db->get_labelTop();
//  labelTop->Draw("same");
//
//  float xMin_label = 0.2;
//  float xMax_label = 0.4;
//
//  TPaveText* label = new TPaveText( xMin_label, 0.83, xMax_label, 0.9, "brNDC" );
//  label->SetTextSize(0.04);
//  label->SetFillColor(0);
//  label->AddText(labelText.c_str());
//  if( labelText!="" )
//    label->Draw("same");
//
//
//  gPad->RedrawAxis();
//
//  char canvasName[500];
//  sprintf( canvasName, "%s/compareTaggers_%s_pt%.0f_%.0f.eps", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax);
//  c1->SaveAs(canvasName);
//  std::string canvasName_str(canvasName);
//  std::string command = "epstopdf " + canvasName_str;
//  system( command.c_str() );
//  sprintf( canvasName, "%s/compareTaggers_%s_pt%.0f_%.0f.png", db->get_outputdir().c_str(), name.c_str(), ptMin, ptMax);
//  c1->SaveAs(canvasName);
//
//
//  delete c1;
//  delete h2_axes;
//  delete legend;
//
//
//
//}



void drawRoC( DrawBase* db, float ptMin, float ptMax, const std::string& flag, TH1D* h1_new_gluon, TH1D* h1_new_quark, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_MLP_gluon, TH1D* h1_MLP_quark, const std::string& labelText, const std::string& legendName_old, const std::string& legendName_new, const std::string& legendName_MLP ) {


  TString flag_tstr(flag);
  bool isCharm = flag_tstr.Contains("charm");

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
    gr_RoC_old->SetMarkerStyle(21);
    gr_RoC_old->SetMarkerColor(29);
  }

  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 ) {
    gr_RoC_MLP->SetMarkerSize(1.3);
    gr_RoC_MLP->SetMarkerStyle(20);
    gr_RoC_MLP->SetMarkerColor(kOrange+1);
  }

  TCanvas* c1 = new TCanvas("c1_roc", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes_roc", "", 10, 0., 1.0001, 10, 0., 1.0001);
  h2_axes->SetXTitle( "Gluon Jet Rejection" );
  if( isCharm )
    h2_axes->SetYTitle( "Charm Jet Efficiency" );
  else
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
    legend->AddEntry( gr_RoC_old, legendName_old.c_str(), "P");
  legend->AddEntry( gr_RoC_new, legendName_new.c_str(), "P");
  if( h1_MLP_quark!=0 && h1_MLP_gluon!=0 )
    legend->AddEntry( gr_RoC_MLP, legendName_MLP.c_str(), "P");
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
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.eps", db->get_outputdir().c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);
  std::string canvasName_eps(canvasName);
  std::string command = "epstopdf " + canvasName_eps;
  system( command.c_str() );
  sprintf( canvasName, "%s/RoC_pt%.0f_%.0f%s.png", db->get_outputdir().c_str(), ptMin, ptMax, flag.c_str());
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;
}


void compareTrees( DrawBase* db, TTree* tree, TTree* tree_herwig, float ptMin, float ptMax, float etaMin, float etaMax ) {


  compareSingleVariable( "ptD_QCJet", "p_{T}D", 50, 0., 1.0001, db, tree, tree_herwig, ptMin, ptMax, etaMin, etaMax );
  compareSingleVariable( "axis1_QCJet", "-ln(Axis_{1})", 50, 1., 7., db, tree, tree_herwig, ptMin, ptMax, etaMin, etaMax, "-log(axis1_QCJet[0])" );
  compareSingleVariable( "axis2_QCJet", "-ln(Axis_{2})", 50, 1., 7., db, tree, tree_herwig, ptMin, ptMax, etaMin, etaMax, "-log(axis2_QCJet[0])" );
  compareSingleVariable( "nPFCand_QC_ptCut", "PFCandidate Multiplicity", 50, 0., 50., db, tree, tree_herwig, ptMin, ptMax, etaMin, etaMax, "nChg_QCJet[0] + nNeutral_ptCutJet[0]" );

  compareSingleVariable( "qgMLPJet", "Quark-Gluon MLP Discriminator", 50, 0., 1.0001, db, tree, tree_herwig, ptMin, ptMax, etaMin, etaMax );

  //compareSingleVariable( "qgl", "Quark-Gluon Likelihood Discriminator", 50, 0., 1.0001, db, tree, tree_herwig, ptMin, ptMax, etaMin, etaMax );


}


void compareSingleVariable( std::string varName, const std::string& axisName, int nbins, float xmin, float xmax, DrawBase* db, TTree* tree, TTree* tree_herwig, float ptMin, float ptMax, float etaMin, float etaMax, const std::string& varExpression ) {

  std::cout << "Herwig-Pythia comparison: " << varName << ", "<< ptMin << " < pt < " << ptMax << " GeV, " << etaMin << " < |eta| < " << etaMax << std::endl;

  TH1D* h1_pythia_quark = new TH1D("pythia_quark", "", nbins, xmin, xmax );
  TH1D* h1_pythia_gluon = new TH1D("pythia_gluon", "", nbins, xmin, xmax );
  TH1D* h1_herwig_quark = new TH1D("herwig_quark", "", nbins, xmin, xmax );
  TH1D* h1_herwig_gluon = new TH1D("herwig_gluon", "", nbins, xmin, xmax );


  if( varName != "qgl" ) {

    char selection[500];
    sprintf( selection, "ptJetGen[0]>%f && ptJetGen[0]<%f && abs(etaJet[0])>%f && abs(etaJet[0])<%f", ptMin, ptMax, etaMin, etaMax);
   
    char selection_quark[500];
    sprintf( selection_quark, "%s && abs(pdgIdJet[0])<4 && abs(pdgIdJet[0])>0", selection );
    char selection_gluon[500];
    sprintf( selection_gluon, "%s && pdgIdJet[0]==21", selection );

    std::string treeVar;
    if( varExpression!="" ) treeVar = varExpression;
    else                    treeVar = varName + "[0]";

    tree->Project("pythia_quark", treeVar.c_str(), selection_quark);
    tree->Project("pythia_gluon", treeVar.c_str(), selection_gluon);
    tree_herwig->Project("herwig_quark", treeVar.c_str(), selection_quark);
    tree_herwig->Project("herwig_gluon", treeVar.c_str(), selection_gluon);

  } else {

    QGLikelihoodCalculator *qglc_tmp = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/ReducedHisto_2012.root");

    // first pythia:
    float pt[20];
    tree->SetBranchAddress("ptJet", pt);
    float ptGen[20];
    tree->SetBranchAddress("ptJetGen", ptGen);
    float eta[20];
    tree->SetBranchAddress("etaJet", eta);
    int pdgId[20];
    tree->SetBranchAddress("pdgIdJet", pdgId);
    float rho;
    tree->SetBranchAddress("rhoPF", &rho);
    float ptD_QC[20];
    tree->SetBranchAddress("ptD_QCJet", ptD_QC);
    float axis2_QC[20];
    tree->SetBranchAddress("axis2_QCJet", axis2_QC);
    int nCharged_QC[20];
    tree->SetBranchAddress("nChg_QCJet", nCharged_QC);
    int nNeutral_ptCut[20];
    tree->SetBranchAddress("nNeutral_ptCutJet", nNeutral_ptCut);

    for( unsigned iEntry=0; iEntry<tree->GetEntries(); ++iEntry ) {

      tree->GetEntry(iEntry);

      if( ptGen[0]<ptMin ) continue;
      if( ptGen[0]>ptMax ) continue;
      if( abs(eta[0])<etaMin ) continue;
      if( abs(eta[0])>etaMax ) continue;

      float qgl_newHisto = qglc_tmp->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);
    
      if( fabs(pdgId[0])<4 && abs(pdgId[0])>0 ) {
        h1_pythia_quark->Fill( qgl_newHisto );
      }
      if( pdgId[0]==21 ) {
        h1_pythia_gluon->Fill( qgl_newHisto );
      }

      if( h1_pythia_quark->GetEntries()>10000 && h1_pythia_gluon->GetEntries()>10000 ) break;

    }
  
    // then herwig:
    tree_herwig->SetBranchAddress("ptJetGen", ptGen);
    tree_herwig->SetBranchAddress("ptJet", pt);
    tree_herwig->SetBranchAddress("etaJet", eta);
    tree_herwig->SetBranchAddress("pdgIdJet", pdgId);
    tree_herwig->SetBranchAddress("rhoPF", &rho);
    tree_herwig->SetBranchAddress("ptD_QCJet", ptD_QC);
    tree_herwig->SetBranchAddress("axis2_QCJet", axis2_QC);
    tree_herwig->SetBranchAddress("nChg_QCJet", nCharged_QC);
    tree_herwig->SetBranchAddress("nNeutral_ptCutJet", nNeutral_ptCut);


    for( unsigned iEntry=0; iEntry<tree_herwig->GetEntries(); ++iEntry ) {

      tree_herwig->GetEntry(iEntry);

      if( ptGen[0]<ptMin ) continue;
      if( ptGen[0]>ptMax ) continue;
      if( abs(eta[0])<etaMin ) continue;
      if( abs(eta[0])>etaMax ) continue;

      float qgl_newHisto = qglc_tmp->computeQGLikelihood2012( pt[0], eta[0], rho, nCharged_QC[0]+nNeutral_ptCut[0], ptD_QC[0], axis2_QC[0]);
    
      if( fabs(pdgId[0])<4 && abs(pdgId[0])>0 ) {
        h1_herwig_quark->Fill( qgl_newHisto );
      }
      if( pdgId[0]==21 ) {
        h1_herwig_gluon->Fill( qgl_newHisto );
      }

      if( h1_herwig_quark->GetEntries()>10000 && h1_herwig_gluon->GetEntries()>10000 ) break;

    }

    delete qglc_tmp;

  } // if var == qgl


  h1_pythia_quark->SetLineColor(38);
  h1_pythia_quark->SetLineWidth(2);
  h1_pythia_quark->SetFillColor(38);
  h1_pythia_quark->SetFillStyle(3005);

  h1_pythia_gluon->SetLineColor(46);
  h1_pythia_gluon->SetLineWidth(2);
  h1_pythia_gluon->SetFillColor(46);
  h1_pythia_gluon->SetFillStyle(3004);

  h1_herwig_quark->SetMarkerStyle(20);
  //h1_herwig_quark->SetMarkerSize(1.5);
  h1_herwig_quark->SetMarkerColor(kBlue+2);

  h1_herwig_gluon->SetMarkerStyle(21);
  //h1_herwig_gluon->SetMarkerSize(1.5);
  h1_herwig_gluon->SetMarkerColor(kRed+2);

  float yMax_quark = h1_pythia_quark->GetMaximum()/h1_pythia_quark->Integral();
  float yMax_gluon = h1_pythia_gluon->GetMaximum()/h1_pythia_gluon->Integral();

  float yMax = (yMax_quark>yMax_gluon) ? yMax_quark : yMax_gluon;
  yMax *= 1.5;


  TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, 0., yMax );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Normalized to Unity" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();

  char legendTitle[300];
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );
  TLegend* legend = new TLegend(0.5, 0.6, 0.9, 0.9, legendTitle);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( h1_pythia_quark, "Quarks (Pythia)", "F" );
  legend->AddEntry( h1_pythia_gluon, "Gluons (Pythia)", "F" );
  legend->AddEntry( h1_herwig_quark, "Quarks (Herwig)", "P" );
  legend->AddEntry( h1_herwig_gluon, "Gluons (Herwig)", "P" );
  legend->Draw("same");

  h1_pythia_quark->DrawNormalized("same");
  h1_pythia_gluon->DrawNormalized("same");
  h1_herwig_quark->DrawNormalized("p same");
  h1_herwig_gluon->DrawNormalized("p same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");


  char labelText[300];
  if( etaMin > 0. )
    sprintf( labelText, "%.0f < |#eta| < %.0f", etaMin, etaMax );
  else
    sprintf( labelText, "|#eta| < %.0f", etaMax );

  TPaveText* label = new TPaveText( 0.2, 0.83, 0.4, 0.9, "brNDC" );
  label->SetTextSize(0.04);
  label->SetFillColor(0);
  label->AddText(labelText);
  if( labelText!="" )
    label->Draw("same");

  gPad->RedrawAxis();

  char canvasName[1000];
  sprintf( canvasName, "%s/herwigPythia_%s_pt%.0f_%.0f_eta%.0f_%.0f.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax, 10.*etaMin, 10.*etaMax);

  c1->SaveAs(canvasName);
  std::string canvasName_str(canvasName);
  std::string command = "epstopdf " + canvasName_str;
  system( command.c_str() );


  TString varName_tstr(varName);
  if( varName_tstr.BeginsWith("qg") ) {
    if( etaMin>2.5 ) {
      drawRoC(db, ptMin, ptMax, "_pythiaHerwig_F", h1_pythia_gluon, h1_pythia_quark, h1_herwig_gluon, h1_herwig_quark, 0, 0, labelText, "Herwig++", "Pythia 6" );
    } else {
      drawRoC(db, ptMin, ptMax, "_pythiaHerwig", h1_pythia_gluon, h1_pythia_quark, h1_herwig_gluon, h1_herwig_quark, 0, 0, labelText, "Herwig++", "Pythia 6" );
    }
  }
  

  delete legend;
  delete c1;
  delete h2_axes;
  delete h1_pythia_quark;
  delete h1_pythia_gluon;
  delete h1_herwig_quark;
  delete h1_herwig_gluon;

}





void drawQuarkFraction_vs_pt( DrawBase* db, TTree* tree, TTree* tree_herwig, float etaMin, float etaMax ) {

  float xMax = (fabs(etaMin)>2.5) ? 210. : 2000.;

  Double_t bins[21];

  fitTools::getBins_int( 21, bins, 20., xMax);



  // pythia begin


  TH1D* h1_denom_pythia    = new TH1D("denom_pythia", "", 20, bins);
  TH1D* h1_quarkNum_pythia = new TH1D("quarkNum_pythia", "", 20, bins);
  TH1D* h1_gluonNum_pythia = new TH1D("gluonNum_pythia", "", 20, bins);
  TH1D* h1_pileupNum_pythia = new TH1D("pileupNum_pythia", "", 20, bins);
  TH1D* h1_undefNum_pythia = new TH1D("undefNum_pythia", "", 20, bins);



  float ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);
  int nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  int nJet;
  tree->SetBranchAddress("nJet", &nJet);
  float ptJet[20];
  tree->SetBranchAddress("ptJet", ptJet);
  float etaJet[20];
  tree->SetBranchAddress("etaJet", etaJet);
  float phiJet[20];
  tree->SetBranchAddress("phiJet", phiJet);
  float eJet[20];
  tree->SetBranchAddress("eJet", eJet);
  float betastarJet[20];
  tree->SetBranchAddress("betastarJet", betastarJet);
  int pdgIdJet[20];
  tree->SetBranchAddress("pdgIdJet", pdgIdJet);


  for( unsigned iEntry=0; iEntry<tree->GetEntries(); ++iEntry ) {

    tree->GetEntry(iEntry);

    if( ptHat < 20. || ptHat > 2000. ) continue;

 
    if( betastarJet[0]>TMath::Log(nvertex-0.67) ) continue;

    if( nJet<2 ) continue;

    TLorentzVector jet1;
    jet1.SetPtEtaPhiE( ptJet[0], etaJet[0], phiJet[0], eJet[0] );
    TLorentzVector jet2;
    jet2.SetPtEtaPhiE( ptJet[1], etaJet[1], phiJet[1], eJet[1] );

    if( fabs(jet1.DeltaPhi(jet2))< 2.5 ) continue;

    float ptAve = 0.5*(ptJet[0] + ptJet[1]);
    if( nJet>2 ) {
      if( ptJet[2]>0.3*ptAve ) continue;
    }


    if( fabs(etaJet[0])>etaMin && fabs(etaJet[0])<etaMax ) {

      h1_denom_pythia->Fill(ptHat);
      if( abs(pdgIdJet[0])<4 && abs(pdgIdJet[0])>0 ) h1_quarkNum_pythia->Fill(ptHat);
      if( abs(pdgIdJet[0])==21 )                     h1_gluonNum_pythia->Fill(ptHat);
      if( abs(pdgIdJet[0])==0  )                     h1_pileupNum_pythia->Fill(ptHat);
      if( pdgIdJet[0]<-100     )                     h1_undefNum_pythia->Fill(ptHat);

    }


    if( fabs(etaJet[1])>etaMin && fabs(etaJet[1])<etaMax ) {

      h1_denom_pythia->Fill(ptHat);
      if( abs(pdgIdJet[1])<4 && abs(pdgIdJet[1])>0 ) h1_quarkNum_pythia->Fill(ptHat);
      if( abs(pdgIdJet[1])==21 )                     h1_gluonNum_pythia->Fill(ptHat);
      if( abs(pdgIdJet[1])==0  )                     h1_pileupNum_pythia->Fill(ptHat);
      if( pdgIdJet[1]<-100     )                     h1_undefNum_pythia->Fill(ptHat);

    }

  }



  TH1D* h1_quarkFraction_pythia = new TH1D(*h1_quarkNum_pythia);
  h1_quarkFraction_pythia->Divide(h1_denom_pythia);

  TH1D* h1_gluonFraction_pythia = new TH1D(*h1_gluonNum_pythia);
  h1_gluonFraction_pythia->Divide(h1_denom_pythia);

  TH1D* h1_pileupFraction_pythia = new TH1D(*h1_pileupNum_pythia);
  h1_pileupFraction_pythia->Divide(h1_denom_pythia);

  TH1D* h1_undefFraction_pythia = new TH1D(*h1_undefNum_pythia);
  h1_undefFraction_pythia->Divide(h1_denom_pythia);

  


  // herwig begin

  TH1D* h1_denom_herwig    = new TH1D("denom_herwig", "", 20, bins);
  TH1D* h1_quarkNum_herwig = new TH1D("quarkNum_herwig", "", 20, bins);
  TH1D* h1_gluonNum_herwig = new TH1D("gluonNum_herwig", "", 20, bins);
  TH1D* h1_pileupNum_herwig = new TH1D("pileupNum_herwig", "", 20, bins);
  TH1D* h1_undefNum_herwig = new TH1D("undefNum_herwig", "", 20, bins);



  tree_herwig->SetBranchAddress("ptHat", &ptHat);
  tree_herwig->SetBranchAddress("nvertex", &nvertex);
  tree_herwig->SetBranchAddress("nJet", &nJet);
  tree_herwig->SetBranchAddress("ptJet", ptJet);
  tree_herwig->SetBranchAddress("etaJet", etaJet);
  tree_herwig->SetBranchAddress("phiJet", phiJet);
  tree_herwig->SetBranchAddress("eJet", eJet);
  tree_herwig->SetBranchAddress("betastarJet", betastarJet);
  tree_herwig->SetBranchAddress("pdgIdJet", pdgIdJet);


  for( unsigned iEntry=0; iEntry<tree_herwig->GetEntries(); ++iEntry ) {

    tree_herwig->GetEntry(iEntry);

    if( ptHat < 20. || ptHat > 2000. ) continue;

 
    if( betastarJet[0]>TMath::Log(nvertex-0.67) ) continue;

    if( nJet<2 ) continue;

    TLorentzVector jet1;
    jet1.SetPtEtaPhiE( ptJet[0], etaJet[0], phiJet[0], eJet[0] );
    TLorentzVector jet2;
    jet2.SetPtEtaPhiE( ptJet[1], etaJet[1], phiJet[1], eJet[1] );

    if( fabs(jet1.DeltaPhi(jet2))< 2.5 ) continue;

    float ptAve = 0.5*(ptJet[0] + ptJet[1]);
    if( nJet>2 ) {
      if( ptJet[2]>0.3*ptAve ) continue;
    }


    if( fabs(etaJet[0])>etaMin && fabs(etaJet[0])<etaMax ) {

      h1_denom_herwig->Fill(ptHat);
      if( abs(pdgIdJet[0])<4 && abs(pdgIdJet[0])>0 ) h1_quarkNum_herwig->Fill(ptHat);
      if( abs(pdgIdJet[0])==21 )                     h1_gluonNum_herwig->Fill(ptHat);
      if( abs(pdgIdJet[0])==0  )                     h1_pileupNum_herwig->Fill(ptHat);
      if( pdgIdJet[0]<-100     )                     h1_undefNum_herwig->Fill(ptHat);

    }


    if( fabs(etaJet[1])>etaMin && fabs(etaJet[1])<etaMax ) {

      h1_denom_herwig->Fill(ptHat);
      if( abs(pdgIdJet[1])<4 && abs(pdgIdJet[1])>0 ) h1_quarkNum_herwig->Fill(ptHat);
      if( abs(pdgIdJet[1])==21 )                     h1_gluonNum_herwig->Fill(ptHat);
      if( abs(pdgIdJet[1])==0  )                     h1_pileupNum_herwig->Fill(ptHat);
      if( pdgIdJet[1]<-100     )                     h1_undefNum_herwig->Fill(ptHat);

    }

  }



  TH1D* h1_quarkFraction_herwig = new TH1D(*h1_quarkNum_herwig);
  h1_quarkFraction_herwig->Divide(h1_denom_herwig);

  TH1D* h1_gluonFraction_herwig = new TH1D(*h1_gluonNum_herwig);
  h1_gluonFraction_herwig->Divide(h1_denom_herwig);

  TH1D* h1_pileupFraction_herwig = new TH1D(*h1_pileupNum_herwig);
  h1_pileupFraction_herwig->Divide(h1_denom_herwig);

  TH1D* h1_undefFraction_herwig = new TH1D(*h1_undefNum_herwig);
  h1_undefFraction_herwig->Divide(h1_denom_herwig);

  
  h1_quarkFraction_pythia->SetMarkerStyle(20);
  h1_quarkFraction_pythia->SetMarkerColor(38);
  h1_quarkFraction_pythia->SetMarkerSize(1.6);

  h1_gluonFraction_pythia->SetMarkerStyle(21);
  h1_gluonFraction_pythia->SetMarkerColor(46);
  h1_gluonFraction_pythia->SetMarkerSize(1.6);

  h1_pileupFraction_pythia->SetMarkerStyle(23);
  h1_pileupFraction_pythia->SetMarkerColor(29);
  h1_pileupFraction_pythia->SetMarkerSize(1.6);

  h1_undefFraction_pythia->SetMarkerStyle(22);
  h1_undefFraction_pythia->SetMarkerColor(kGray+1);
  h1_undefFraction_pythia->SetMarkerSize(1.6);

  h1_quarkFraction_herwig->SetMarkerStyle(24);
  h1_quarkFraction_herwig->SetMarkerColor(38);
  h1_quarkFraction_herwig->SetMarkerSize(1.6);

  h1_gluonFraction_herwig->SetMarkerStyle(25);
  h1_gluonFraction_herwig->SetMarkerColor(46);
  h1_gluonFraction_herwig->SetMarkerSize(1.6);

  h1_undefFraction_herwig->SetMarkerStyle(26);
  h1_undefFraction_herwig->SetMarkerColor(kGray+1);
  h1_undefFraction_herwig->SetMarkerSize(1.6);



  char legendTitle[200];
  if( etaMin==0. )
    sprintf( legendTitle, "|#eta| < %.0f", etaMax );
  else
    sprintf( legendTitle, "%.0f < |#eta| < %.0f", etaMin, etaMax );

  TLegend* legend = new TLegend(0.2, 0.68, 0.5, 0.93, legendTitle );
  legend->SetTextSize(0.038);
  legend->SetFillColor(0);
  legend->AddEntry( h1_quarkFraction_pythia, "Quark", "P" );
  legend->AddEntry( h1_gluonFraction_pythia, "Gluon", "P" );
  legend->AddEntry( h1_undefFraction_pythia, "Undefined", "P" );
  legend->AddEntry( h1_pileupFraction_pythia, "Pile Up", "P" );

  TLegend* legend_herwig = new TLegend(0.47, 0.75, 0.77, 0.93 );
  legend_herwig->SetTextSize(0.038);
  legend_herwig->SetFillColor(0);
  legend_herwig->AddEntry( h1_quarkFraction_herwig, "Quark (Herwig)", "P" );
  legend_herwig->AddEntry( h1_gluonFraction_herwig, "Gluon (Herwig)", "P" );
  legend_herwig->AddEntry( h1_undefFraction_herwig, "Undefined (Herwig)", "P" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  c1->SetLogx();

  TH2D* h2_axes = new TH2D("axes", "", 10, 20., xMax, 10, 0., 1.15);
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->SetXTitle("#hat{p}_{T} [GeV]");
  h2_axes->SetYTitle("Flavor Fraction");

  h2_axes->Draw();

  h1_gluonFraction_herwig->Draw("p same");
  h1_quarkFraction_herwig->Draw("p same");
  h1_undefFraction_herwig->Draw("p same");
  h1_gluonFraction_pythia->Draw("p same");
  h1_quarkFraction_pythia->Draw("p same");
  h1_pileupFraction_pythia->Draw("p same");
  h1_undefFraction_pythia->Draw("p same");

  legend->Draw("same");
  legend_herwig->Draw("same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  gPad->RedrawAxis();

  char canvasName[512];
  sprintf( canvasName, "%s/fractionsVsPt_eta%.0f%.0f.eps", db->get_outputdir().c_str(), etaMin, etaMax );
  c1->SaveAs( canvasName );

  delete c1;

  delete legend;
  delete legend_herwig;
  delete h2_axes;
  delete h1_gluonFraction_herwig;
  delete h1_quarkFraction_herwig;
  delete h1_undefFraction_herwig;
  delete h1_pileupFraction_herwig;

  delete h1_gluonFraction_pythia;
  delete h1_quarkFraction_pythia;
  delete h1_pileupFraction_pythia;
  delete h1_undefFraction_pythia;

  delete h1_denom_pythia;
  delete h1_denom_herwig;

  delete h1_gluonNum_herwig;
  delete h1_quarkNum_herwig;
  delete h1_undefNum_herwig;
  delete h1_pileupNum_herwig;

  delete h1_gluonNum_pythia;
  delete h1_quarkNum_pythia;
  delete h1_pileupNum_pythia;
  delete h1_undefNum_pythia;

}
