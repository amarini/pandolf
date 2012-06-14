#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TArrow.h"

#include "DrawBase.h"


void drawSingleLikelihoodPlot( const std::string& suffix );


int main() {


  drawSingleLikelihoodPlot( "trilepton_4channels" );
  drawSingleLikelihoodPlot( "ssdl_3channels" );
  drawSingleLikelihoodPlot( "7channels" );

  return 0;

}



void drawSingleLikelihoodPlot( const std::string& suffix ) {

  DrawBase* db = new DrawBase(suffix);
  db->set_lumi(4980.); 

  std::string fileName = "q_mergedToys_" + suffix + ".root";

  TFile* qFile = TFile::Open(fileName.c_str());

  TTree* qTree = (TTree*)qFile->Get("q");

  float xMax = 30.;

  TH1D* h1_bg = new TH1D("bg", "", 50, 0., xMax);
  TH1D* h1_obs = new TH1D("obs", "", 50, 0., xMax);

  qTree->Project("bg", "2.*q", "type==-1");
  qTree->Project("obs", "2.*q", "type==0");

  h1_bg->SetFillStyle(3004);
  h1_bg->SetFillColor(39);
  h1_bg->SetLineColor(39);
  h1_bg->SetLineWidth(2.);

  float yMin = 0.5;
  TH2D* h2_axes = new TH2D("axes", "", 10, 0., xMax, 10, yMin, 10.*h1_bg->GetMaximum());
  h2_axes->SetXTitle("2 ln(L)");
  h2_axes->SetYTitle("Number of Toys");

  TPaveText* labelSqrt = db->get_labelSqrt();

  float xArrow = h1_obs->GetMean();

  TArrow* arrow_obs = new TArrow( xArrow, 500., xArrow, yMin, 0.05, ">" );
  arrow_obs->SetLineColor(46);
  arrow_obs->SetFillColor(46);
  arrow_obs->SetLineWidth(2.);

  std::string legendTitle;
  if( suffix=="trilepton_4channels" ) legendTitle = "Trilepton Channel";
  else if( suffix=="ssdl_3channels" ) legendTitle = "Dilepton Channel";
  else                                legendTitle = "7 Channels Combined";

  TLegend* legend = new TLegend( 0.4, 0.7, 0.72, 0.9, legendTitle.c_str());
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_bg, "Background Hypothesis", "F" );
  legend->AddEntry( arrow_obs, "Observed", "L" );

  TCanvas* c1 = new TCanvas( "c1", "", 600, 600);
  c1->cd();
  c1->SetLogy();

  h2_axes->Draw();
  h1_bg->Draw("same");
  arrow_obs->Draw();
  labelSqrt->Draw("same");
  legend->Draw("Same");


  gPad->RedrawAxis();

  std::string canvasName = "logLikelihood_" + suffix + ".eps";

  c1->SaveAs(canvasName.c_str());

  delete c1;
  delete h2_axes;
  delete legend;
  delete h1_bg;
  delete h1_obs;
  delete arrow_obs;

}

