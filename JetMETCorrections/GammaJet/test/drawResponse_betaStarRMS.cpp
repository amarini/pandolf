#include <iostream>
#include <cstdlib>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"

#include "DrawBase.h"

void drawResponsePlot( DrawBase* db, float ptMin, float ptMax );

int main() {


  DrawBase* db = new DrawBase("resp");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(4980.);

  TFile* file = TFile::Open("Omog_QGStudies_Photon_Run2011_FULL.root");
  db->add_dataFile( file, "data" );

  
  drawResponsePlot( db, 30., 50. );
  drawResponsePlot( db, 50., 100. );
  drawResponsePlot( db, 100., 150. );


  return 0;

}



void drawResponsePlot( DrawBase* db, float ptMin, float ptMax ) {

  TFile* file = db->get_dataFile(0).file;
  TTree* tree = (TTree*)file->Get("omog");

  char ptSelection[200];
  sprintf( ptSelection, "ptJet0 > %f && ptJet0 < %f", ptMin, ptMax );

  std::string varSelection = "( (abs(etaJet0)<2.5 && exp(-rmsCandJet0)<0.03 && betaStarJet0<0.2*log( nvertex-0.67 ) ) || ( abs(etaJet0)>2.5 && abs(etaJet0)<3. && exp(-rmsCandJet0)<0.03) || (abs(etaJet0)>3. && abs(etaJet0)<4.7 && exp(-rmsCandJet0)<0.04) )";

  char beforeSel[200];
  sprintf( beforeSel, "%s", ptSelection );

  char afterSel[500];
  sprintf( afterSel, "%s && %s", ptSelection, varSelection.c_str() );

  TH1D* h1_response_before = new TH1D( "response_before", "", 36, 0.2, 2.);
  h1_response_before->Sumw2();
  TH1D* h1_response_after = new TH1D( "response_after", "", 36, 0.2, 2.);
  h1_response_after->Sumw2();

  tree->Project( "response_before", "ptJet0/ptPhot", beforeSel );
  tree->Project( "response_after", "ptJet0/ptPhot", afterSel );

  h1_response_before->SetMarkerStyle(24);
  h1_response_before->SetMarkerSize(1.6);
  h1_response_before->SetMarkerColor(kBlack);

  h1_response_after->SetMarkerStyle(20);
  h1_response_after->SetMarkerSize(1.6);
  h1_response_after->SetMarkerColor(46);

  TPaveText* label_sqrt = db->get_labelSqrt();

  char ptSelection_text[200];
  sprintf( ptSelection_text, "%.0f < p_{T} < %.0f GeV", ptMin, ptMax );

  TLegend* legend = new TLegend( 0.6, 0.65, 0.88, 0.88, ptSelection_text );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_response_before, "All Jets", "P" );
  legend->AddEntry( h1_response_after, "Selected", "P" );

  float yMax = 1.3*h1_response_after->GetMaximum()/h1_response_after->Integral();

  TH2D* h2_axes = new TH2D("axes", "", 10, 0.2, 2., 10, 0., yMax );
  h2_axes->SetXTitle( "Photon p_{T} / Jet p_{T}");
  h2_axes->SetYTitle( "Normalized to Unity");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  h1_response_before->DrawNormalized("e same");
  h1_response_after->DrawNormalized("e same");

  legend->Draw("Same");
  label_sqrt->Draw("Same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "response_betaStar_rms_pt%.0f_%.0f.eps", ptMin, ptMax );
  c1->SaveAs(canvasName);

  delete legend;
  delete c1;
  delete h2_axes;
  delete h1_response_before;
  delete h1_response_after;

}
