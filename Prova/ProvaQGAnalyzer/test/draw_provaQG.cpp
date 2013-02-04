#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLegend.h"
#include "DrawBase.h"


bool addCHS = true;

void drawOne( DrawBase* db, const std::string& suffix, TFile* file, TFile* file_chs=0 );
void drawRoC( DrawBase* db, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_central_gluon, TH1D* h1_central_quark, TH1D* h1_new_gluon=0, TH1D* h1_new_quark=0, const std::string& labelText="" );
void compareOneVariable( DrawBase* db, TTree* tree, TTree* tree_etaFix, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax );
void comparePDF( DrawBase* db, const std::string& name, const std::string& savename, TH1D* h1_nPFCand_quark_central, TH1D* h1_nPFCand_gluon_central, TH1D* h1_nPFCand_quark_transition_fix, TH1D* h1_nPFCand_gluon_transition_fix, const std::string& flags="" );
void drawCompare_etaFix( DrawBase* db, TFile* file, TFile* file_etaFix );



int main() {

  DrawBase* db = new DrawBase("provaQG");
  db->set_outputdir("prova");

  db->set_rebin(2);

  TFile* file = TFile::Open("provaQG.root");
  db->add_mcFile(file, "prova", "prova");

  TFile* file_chs = 0;
  if( addCHS ) file_chs = TFile::Open("provaQG_CHS.root");

  drawOne( db, "quark_pt3050", file, file_chs );
  drawOne( db, "gluon_pt3050", file, file_chs );

  drawOne( db, "quark_pt80120", file, file_chs );
  drawOne( db, "gluon_pt80120", file, file_chs );

  drawOne( db, "quark_pt200300", file, file_chs );
  drawOne( db, "gluon_pt200300", file, file_chs );

  drawOne( db, "quark_pt2050_F", file, file_chs );
  drawOne( db, "gluon_pt2050_F", file, file_chs );

  drawOne( db, "quark_pt50100_F", file, file_chs );
  drawOne( db, "gluon_pt50100_F", file, file_chs );

  TFile* file_etaFix = TFile::Open("provaQG_etaFix.root");
  drawCompare_etaFix( db, file, file_etaFix );

  return 0;

}


void drawOne( DrawBase* db, const std::string& suffix, TFile* file, TFile* file_chs ) {

  //std::vector< HistoAndName > hn;

  //HistoAndName hn_old;
  //hn_old.histoName = "qglOLD_" + suffix;
  //hn_old.legendName = "Histo based LD";
  //hn.push_back(hn_old);

  //HistoAndName hn_new;
  //hn_new.histoName = "qglNEW_" + suffix;
  //hn_new.legendName = "Fit based LD";
  //hn_new.markerStyle = 20;
  //hn.push_back(hn_new);

  //db->compareDifferentHistos( hn, suffix, "QG LD" );

  std::string histoLD_name = "qglOLD_" + suffix;
  std::string fitLD_name = "qglNEW_" + suffix;

  TH1D* h1_histoLD  = (TH1D*)file->Get( histoLD_name.c_str() );
  TH1D* h1_fitLD    = (TH1D*)file->Get( fitLD_name.c_str() );
  TH1D* h1_fitLD_CHS = 0;
  if( addCHS )
    h1_fitLD_CHS = (TH1D*)file_chs->Get( fitLD_name.c_str() );

  h1_histoLD->Rebin(2);
  h1_fitLD->Rebin(2);
  if( addCHS )
    h1_fitLD_CHS->Rebin(2);

  h1_histoLD->SetLineWidth(2);
  h1_histoLD->SetLineColor(38);
  h1_histoLD->SetFillColor(38);
  h1_histoLD->SetFillStyle(3004);

  h1_fitLD->SetMarkerSize(1.6);
  h1_fitLD->SetMarkerStyle(20);
  h1_fitLD->SetMarkerColor(46);
  h1_fitLD->SetFillColor(46);

  if( addCHS ) {
    h1_fitLD_CHS->SetLineWidth(2);
    h1_fitLD_CHS->SetMarkerSize(1.6);
    h1_fitLD_CHS->SetMarkerStyle(24);
    h1_fitLD_CHS->SetMarkerColor(kRed+3);
    h1_fitLD_CHS->SetFillColor(kRed+3);
  }

  float yMax = (!addCHS) ? h1_fitLD->GetMaximum()*1.3/h1_fitLD->Integral() : h1_fitLD_CHS->GetMaximum()*1.3/h1_fitLD_CHS->Integral();


  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., yMax);
  h2_axes->SetXTitle("Quark-Gluon LD");
  h2_axes->SetYTitle("Normalized to Unity");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TLegend* legend = new TLegend(0.35, 0.6, 0.65, 0.9);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry(h1_histoLD, "Histo-based LD", "F");
  legend->AddEntry(h1_fitLD, "Fit-based LD", "P");
  if( addCHS )
    legend->AddEntry(h1_fitLD_CHS, "Fit-based LD (CHS)", "P");

  h2_axes->Draw();
  legend->Draw("same");
  h1_histoLD->DrawNormalized("same");
  h1_fitLD->DrawNormalized("p same");
  if( addCHS ) 
    h1_fitLD_CHS->DrawNormalized("p same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");




  std::string canvasName = db->get_outputdir() + "/qgl_prova_" + suffix;
  if( addCHS ) canvasName += "_CHS";
  canvasName += ".eps";
  
  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h2_axes;
  delete legend;

}


void drawCompare_etaFix( DrawBase* db, TFile* file, TFile* file_etaFix ) {

  TTree* tree = (TTree*)file->Get("tree_out");
  TTree* tree_etaFix = (TTree*)file_etaFix->Get("tree_out");

  compareOneVariable( db, tree, tree_etaFix, "nPFCand_QC_ptCut", "PF Candidate Multiplicity", 30, -0.5, 49.5 );
  compareOneVariable( db, tree, tree_etaFix, "axis2_QC", "Jet Minor Axis", 30, 0., 5. );
  compareOneVariable( db, tree, tree_etaFix, "ptD_QC", "Jet p_{T}D", 30, 0., 1.0001 );
  compareOneVariable( db, tree, tree_etaFix, "pt", "Jet p_{T} [GeV]", 30, 20., 300.);
  compareOneVariable( db, tree, tree_etaFix, "qglNEW", "Quark-Gluon LD", 50, 0., 1.0001);

}




void compareOneVariable( DrawBase* db, TTree* tree, TTree* tree_etaFix, const std::string& varName, const std::string& axisName, int nBins, float xMin, float xMax ) {


  std::string histoName_quark_central = varName + "quark_central";
  TH1D* h1_quark_central = new TH1D(histoName_quark_central.c_str(), "", nBins, xMin, xMax);
  std::string histoName_gluon_central = varName + "gluon_central";
  TH1D* h1_gluon_central = new TH1D(histoName_gluon_central.c_str(), "", nBins, xMin, xMax);

  std::string histoName_quark_transition = varName + "quark_transition";
  TH1D* h1_quark_transition = new TH1D(histoName_quark_transition.c_str(), "", nBins, xMin, xMax);
  std::string histoName_gluon_transition = varName + "gluon_transition";
  TH1D* h1_gluon_transition = new TH1D(histoName_gluon_transition.c_str(), "", nBins, xMin, xMax);

  std::string histoName_quark_transition_fix = varName + "quark_transition_fix";
  TH1D* h1_quark_transition_fix = new TH1D(histoName_quark_transition_fix.c_str(), "", nBins, xMin, xMax);
  std::string histoName_gluon_transition_fix = varName + "gluon_transition_fix";
  TH1D* h1_gluon_transition_fix = new TH1D(histoName_gluon_transition_fix.c_str(), "", nBins, xMin, xMax);


  tree->Project(histoName_quark_central.c_str(), varName.c_str(), "pdgid>0 && pdgid<4 && pt>20. && pt<300.&& abs(eta)<2.");
  tree->Project(histoName_gluon_central.c_str(), varName.c_str(), "pdgid==21 && pt>20. && pt<300.&& abs(eta)<2.");

  tree->Project(histoName_quark_transition.c_str(), varName.c_str(), "pdgid>0 && pdgid<4 && pt>20. && pt<300.&& abs(eta)>2. && abs(eta)<2.4");
  tree->Project(histoName_gluon_transition.c_str(), varName.c_str(), "pdgid==21 && pt>20. && pt<300.&& abs(eta)>2. && abs(eta)<2.4");

  tree_etaFix->Project(histoName_quark_transition_fix.c_str(), varName.c_str(), "pdgid>0 && pdgid<4 && pt>20. && pt<300.&& abs(eta)>2. && abs(eta)<2.4");
  tree_etaFix->Project(histoName_gluon_transition_fix.c_str(), varName.c_str(), "pdgid==21 && pt>20. && pt<300.&& abs(eta)>2. && abs(eta)<2.4");

  comparePDF( db, axisName, varName, h1_quark_central, h1_gluon_central, h1_quark_transition, h1_gluon_transition, "nofix" );
  comparePDF( db, axisName, varName, h1_quark_central, h1_gluon_central, h1_quark_transition_fix, h1_gluon_transition_fix, "etafix" );

  delete h1_quark_central;
  delete h1_gluon_central;

  delete h1_quark_transition;
  delete h1_gluon_transition;

  delete h1_quark_transition_fix;
  delete h1_gluon_transition_fix;

}



void comparePDF( DrawBase* db, const std::string& name, const std::string& savename, TH1D* h1_nPFCand_quark_central, TH1D* h1_nPFCand_gluon_central, TH1D* h1_nPFCand_quark_transition_fix, TH1D* h1_nPFCand_gluon_transition_fix, const std::string& flags ) {

  h1_nPFCand_quark_central->SetLineWidth(2);
  h1_nPFCand_quark_central->SetLineColor(38);
  h1_nPFCand_quark_central->SetFillColor(38);
  h1_nPFCand_quark_central->SetFillStyle(3004);

  h1_nPFCand_gluon_central->SetLineWidth(2);
  h1_nPFCand_gluon_central->SetLineColor(46);
  h1_nPFCand_gluon_central->SetFillColor(46);
  h1_nPFCand_gluon_central->SetFillStyle(3005);

  h1_nPFCand_quark_transition_fix->SetMarkerSize(1.4);
  h1_nPFCand_quark_transition_fix->SetMarkerColor(kBlue+1);
  h1_nPFCand_quark_transition_fix->SetMarkerStyle(24);

  h1_nPFCand_gluon_transition_fix->SetMarkerSize(1.4);
  h1_nPFCand_gluon_transition_fix->SetMarkerColor(kRed+2);
  h1_nPFCand_gluon_transition_fix->SetMarkerStyle(20);

  float xmin = h1_nPFCand_quark_transition_fix->GetXaxis()->GetXmin();
  float xmax = h1_nPFCand_quark_transition_fix->GetXaxis()->GetXmax();
  float ymax = h1_nPFCand_quark_central->GetMaximum()*1.6/h1_nPFCand_quark_central->Integral();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, 0., ymax);
  h2_axes->SetXTitle( name.c_str() );
  h2_axes->SetYTitle( "Normalized to Unity" );

  h2_axes->Draw();
 
  h1_nPFCand_quark_central->DrawNormalized("same");
  h1_nPFCand_gluon_central->DrawNormalized("same");

  h1_nPFCand_quark_transition_fix->DrawNormalized("p same");
  h1_nPFCand_gluon_transition_fix->DrawNormalized("p same");

  TPaveText* label_top = db->get_labelTop();
  label_top->Draw("same");

  TLegend* legend = new TLegend( 0.55, 0.55, 0.9, 0.9 );
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( h1_nPFCand_quark_central, "Quarks (|#eta|<2)", "F");
  legend->AddEntry( h1_nPFCand_gluon_central, "Gluons (|#eta|<2)", "F");
  legend->AddEntry( h1_nPFCand_quark_transition_fix, "Quarks (2<|#eta|<2.4)", "P");
  legend->AddEntry( h1_nPFCand_gluon_transition_fix, "Gluons (2<|#eta|<2.4)", "P");
  legend->Draw("same");

  std::string canvasName = db->get_outputdir() + "/" + savename + "_transition";
  if( flags!="" ) canvasName += "_" + flags;
  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";
  c1->SaveAs(canvasName_eps.c_str());
  c1->SaveAs(canvasName_png.c_str());

  delete c1;
  delete legend;
  delete h2_axes;

}


void drawRoC( DrawBase* db, TH1D* h1_old_gluon, TH1D* h1_old_quark, TH1D* h1_central_gluon, TH1D* h1_central_quark, TH1D* h1_new_gluon, TH1D* h1_new_quark, const std::string& labelText ) {


  TGraph* gr_RoC_old = new TGraph(0);
  TGraph* gr_RoC_new = new TGraph(0);
  TGraph* gr_RoC_central = new TGraph(0);

  int nbins = h1_old_quark->GetNbinsX();

  for( unsigned int ibin=1; ibin<nbins+1; ++ibin ) {

    float eff_q_old = -1.;
    float eff_g_old = -1.;
  
    if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
      eff_q_old = h1_old_quark->Integral( nbins-ibin, nbins )/h1_old_quark->Integral( 1, nbins );
      eff_g_old = h1_old_gluon->Integral( nbins-ibin, nbins )/h1_old_gluon->Integral( 1, nbins );
    }
  
    float eff_q_central = -1.;
    float eff_g_central = -1.;
  
    if( h1_central_quark!=0 && h1_central_gluon!=0 ) { 
      eff_q_central = h1_central_quark->Integral( nbins-ibin, nbins )/h1_central_quark->Integral( 1, nbins );
      eff_g_central = h1_central_gluon->Integral( nbins-ibin, nbins )/h1_central_gluon->Integral( 1, nbins );
    }

    float eff_q_new = -1.;
    float eff_g_new = -1.;

    if( h1_new_quark!=0 && h1_new_gluon!=0 ) {
      eff_q_new = h1_new_quark->Integral( nbins-ibin, nbins )/h1_new_quark->Integral( 1, nbins );
      eff_g_new = h1_new_gluon->Integral( nbins-ibin, nbins )/h1_new_gluon->Integral( 1, nbins );
    }
  
    if( h1_new_quark!=0 && h1_new_gluon!=0 ) 
      gr_RoC_new->SetPoint( ibin-1, 1.-eff_g_new, eff_q_new );

    if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
      gr_RoC_old->SetPoint( ibin-1, 1.-eff_g_old, eff_q_old );

    if( h1_central_quark!=0 && h1_central_gluon!=0 ) 
      gr_RoC_central->SetPoint( ibin-1, 1.-eff_g_central, eff_q_central );

  }


  if( h1_new_quark!=0 && h1_new_gluon!=0 ) {
    gr_RoC_new->SetMarkerSize(1.3);
    gr_RoC_new->SetMarkerStyle(24);
    gr_RoC_new->SetMarkerColor(kRed+3);
  }

  if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
    gr_RoC_old->SetMarkerSize(1.3);
    gr_RoC_old->SetMarkerStyle(20);
    gr_RoC_old->SetMarkerColor(kOrange+1);
  }

  if( h1_central_quark!=0 && h1_central_gluon!=0 ) {
    gr_RoC_central->SetMarkerSize(1.3);
    gr_RoC_central->SetMarkerStyle(21);
    gr_RoC_central->SetMarkerColor(29);
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
  sprintf( legendTitle, "20 < p_{T} < 300 GeV" );
  TLegend* legend = new TLegend( 0.2, 0.2, 0.45, 0.45, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  if( h1_central_quark!=0 && h1_central_gluon!=0 )
    legend->AddEntry( gr_RoC_central, "|#eta|<2", "P");
  if( h1_old_quark!=0 && h1_old_gluon!=0 ) {
    if( h1_new_quark==0 )
      legend->AddEntry( gr_RoC_old, "2<|#eta|<2.5", "P");
    else
      legend->AddEntry( gr_RoC_old, "2<|#eta|<2.5 Before Correction", "P");
  }
  if( h1_new_quark!=0 && h1_new_gluon!=0 )
    legend->AddEntry( gr_RoC_new, "2<|#eta|<2.5 After Correction", "P");
  legend->Draw("same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");

  TPaveText* label = new TPaveText( 0.7, 0.83, 0.9, 0.9, "brNDC" );
  label->SetTextSize(0.04);
  label->SetFillColor(0);
  label->AddText(labelText.c_str());
  if( labelText!="" )
    label->Draw("same");

  
  if( h1_central_quark!=0 && h1_central_gluon!=0 ) 
    gr_RoC_central->Draw("p same");
  if( h1_old_quark!=0 && h1_old_gluon!=0 ) 
    gr_RoC_old->Draw("p same");
  if( h1_new_quark!=0 && h1_new_gluon!=0 ) 
    gr_RoC_new->Draw("p same");

  gPad->RedrawAxis();

  char canvasName[500];
  if( h1_new_quark!=0 && h1_new_gluon!=0 ) {
    sprintf( canvasName, "%s/RoC_etaFix.eps", db->get_outputdir().c_str());
    c1->SaveAs(canvasName);
    sprintf( canvasName, "%s/RoC_etaFix.png", db->get_outputdir().c_str());
    c1->SaveAs(canvasName);
  } else {
    sprintf( canvasName, "%s/RoC_noFix.eps", db->get_outputdir().c_str());
    c1->SaveAs(canvasName);
    sprintf( canvasName, "%s/RoC_noFix.png", db->get_outputdir().c_str());
    c1->SaveAs(canvasName);
  }

  delete c1;
  delete h2_axes;
  delete legend;
}

