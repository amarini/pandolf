#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.h"
#include "DrawBase.h"



int main( int argc, char* argv[]) {

  if( argc!=2 ) {
    std::cout << "USAGE: ./draw_dijet_gammajet.C [recoType]" << std::endl;
    exit(9191);
  }

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

  std::string recoType(argv[1]);

  std::string jetalgo_dijet = "ak5"+recoType;
  std::string jetalgo_gammajet;
  if( recoType=="pf" ) {
    jetalgo_gammajet = recoType+ "akt5";
  } else if( recoType=="calo" ) {
    jetalgo_gammajet = "akt5";
  } else if( recoType=="jpt" ) {
    jetalgo_gammajet = recoType+"ak5";
  } else {
    std::cout << "Unknown recoType '" << recoType<< "'. Exiting." << std::endl;
    exit(878);
  }

  DrawBase* db = new DrawBase("Dijet_Gammajet", recoType, "akt5");

  //std::string dijetFileName = "dj_resolution_17012011/data_"+jetalgo_dijet+"_final.root";
  std::string dijetFileName = "dj_resolution_10022011_NEW/data_"+jetalgo_dijet+"_final.scaled.root";
  TFile* file_dijet = TFile::Open(dijetFileName.c_str());
  std::cout << "-> Opened file '" << dijetFileName << "'." << std::endl;
  std::string gammajetFileName = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+jetalgo_gammajet+"_LUMI_2ndJet10_eta011L2L3Raw_RMS99.root";
  //std::string gammajetFileName = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+jetalgo_gammajet+"_LUMI_2ndJet10_eta011L2L3_RMS99.root";
  TFile* file_gammajet = TFile::Open(gammajetFileName.c_str());
  std::cout << "-> Opened file '" << gammajetFileName << "'." << std::endl;

  db->add_dataFile(file_dijet, "DiJets");
  db->add_dataFile(file_gammajet, "GammaJets");

  db->set_lumiNormalization(34.);
  
  TDirectory* d = (TDirectory*)file_dijet->Get(jetalgo_dijet.c_str());
  TGraphErrors* graph_dijet = (TGraphErrors*)d->Get("RelResVsRefPt_JetEta0to1.1");
  TGraphErrors* graph_gammajet = (TGraphErrors*)file_gammajet->Get("gr_DATAReso_vs_pt");
  TGraphErrors* graph_intrinsic = (TGraphErrors*)file_gammajet->Get("gr_intrReso_vs_pt");

  if( recoType=="calo" ) {
    graph_gammajet->RemovePoint(graph_gammajet->GetN()-1);
    graph_gammajet->RemovePoint(0);
    graph_gammajet->RemovePoint(0);
    graph_gammajet->RemovePoint(0);
  }

  graph_dijet->SetMarkerStyle(20);
  graph_dijet->SetMarkerSize(1.8);
  graph_dijet->SetMarkerColor(38);

  graph_gammajet->SetMarkerStyle(22);
  graph_gammajet->SetMarkerSize(1.8);
  graph_gammajet->SetMarkerColor(46);

  float xMin = 32.;
  float xMax = 600.;

  TF1* fit_intrReso = fitTools::fitResolutionGraph(graph_intrinsic, "NSCPF", "fit_intrReso", "NR", xMax, xMin);
  fit_intrReso->SetLineColor(kBlack);
  fit_intrReso->SetLineWidth(2);
  fit_intrReso->SetLineStyle(2);
  TH1D* band_intrReso = fitTools::getBand(fit_intrReso, "band_intrReso");
  band_intrReso->SetFillColor(kYellow-9);
  band_intrReso->SetLineColor(kBlack);
  band_intrReso->SetLineWidth(2);
  band_intrReso->SetLineStyle(2);
  

  TH2D* axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 0.3);
  axes->GetXaxis()->SetTitleOffset(1.1);
  axes->GetYaxis()->SetTitleOffset(1.5);
  axes->SetXTitle("Transverse Momentum [GeV/c]");
  axes->SetYTitle("Jet p_{T} Resolution");
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();


  TPaveText* label_cms = db->get_labelCMS(2);
  TPaveText* label_sqrt = db->get_labelSqrt(2);
  TPaveText* label_algo = db->get_labelAlgo(3);

  TLegend* legend = new TLegend(0.53, 0.7, 0.88, 0.88, "  |#eta| < 1.1");
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( graph_dijet, "DiJet Asymmetry", "P" );
  legend->AddEntry( graph_gammajet, "Photon+Jet", "P" );
  //legend->AddEntry( fit_intrReso, "MC Truth", "L" );
  //legend->AddEntry( band_intrReso, "MC Truth", "FL" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->SetLogx();
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  axes->Draw();
//band_intrReso->Draw("c e3 same");
//fit_intrReso->Draw("c e3 same");
  graph_dijet->Draw("p same"); 
  graph_gammajet->Draw("p same"); 
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");
  gPad->RedrawAxis();

  std::string canvName = "reso_dijet_gammajet_" + jetalgo_gammajet + ".eps";
  c1->SaveAs(canvName.c_str());

  return 0;

}
