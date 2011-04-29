#include <cstdlib>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "DrawBase.h"



void drawSinglePlot( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisTitle, const std::string& rhoType );


int main( int argc, char* argv[] )  {


  if( argc != 2 ) {
    std::cout << "USAGE: ./drawRhoStudies [dataset]" << std::endl;
    exit(13);
  }


  std::string dataset(argv[1]);

  std::string fileName = "RhoStudies_"+dataset+".root";
  TFile* file = TFile::Open(fileName.c_str());

  DrawBase* db = new DrawBase("RhoStudies");
  std::string outputdirName = "RhoStudiesPlots_"+dataset;
  db->set_outputdir(outputdirName);

  drawSinglePlot( db, file, "hcalIsoPID", "HCAL Isolation Threshold [GeV]", "PF" );
  drawSinglePlot( db, file, "ecalIsoPID", "ECAL Isolation Threshold [GeV]", "PF" );

  return 0;

}


void drawSinglePlot( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisTitle, const std::string& rhoType ) {

  std::string histoThreshName = varName+"Thresh_vs_rho"+rhoType;
  std::string histoEffName = varName+"Eff_vs_rho"+rhoType;

  TH1D* histoThresh = (TH1D*)file->Get(histoThreshName.c_str());
  TH1D* histoEff = (TH1D*)file->Get(histoEffName.c_str());

  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();


  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  // first draw efficiency:

  TH2D* h2_axes_eff = new TH2D("axes_eff", "", 10, histoThresh->GetXaxis()->GetXmin(), histoThresh->GetXaxis()->GetXmax(), 10, 0., 1.1);
  if( rhoType=="PF" )
    h2_axes_eff->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes_eff->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  h2_axes_eff->SetYTitle("Efficiency");
  h2_axes_eff->Draw();

  histoEff->SetMarkerStyle(29);
  histoEff->SetMarkerSize(1.8);
  histoEff->SetMarkerColor(38);
  
  TF1* constLine = new TF1("constLine", "[0]");
  constLine->SetLineStyle(2);
  constLine->SetLineWidth(1);

  histoEff->Fit( constLine, "Q" );
  histoEff->Draw("p same");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");

  TPaveText* label_fitResults_eff = new TPaveText(0.2, 0.22, 0.5, 0.3, "brNDC");
  label_fitResults_eff->SetFillColor(0);
  label_fitResults_eff->SetTextSize(0.035);
  char label_fitResultsText_eff[300];
  sprintf( label_fitResultsText_eff, "FIT = %.3f #pm %.3f", constLine->GetParameter(0), constLine->GetParError(0) );
  label_fitResults_eff->AddText(label_fitResultsText_eff);
  label_fitResults_eff->Draw("same");

  std::string canvasName_eff = db->get_outputdir() + "/" + varName + "Eff_vs_PF";
  std::string canvasName_eff_png = canvasName_eff + ".png";
  std::string canvasName_eff_eps = canvasName_eff + ".eps";

  c1->SaveAs(canvasName_eff_png.c_str());
  c1->SaveAs(canvasName_eff_eps.c_str());

  
  c1->Clear();
  

  // then the thresholds:

  TH2D* h2_axes = new TH2D("axes", "", 10, histoThresh->GetXaxis()->GetXmin(), histoThresh->GetXaxis()->GetXmax(), 10, 0., 1.5*histoThresh->GetMaximum() );
  if( rhoType=="PF" )
    h2_axes->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  h2_axes->SetYTitle(yAxisTitle.c_str());
  h2_axes->Draw();

  histoThresh->SetMarkerStyle(29);
  histoThresh->SetMarkerSize(1.8);
  histoThresh->SetMarkerColor(kRed+1);
  
  std::string funcName = "line"+rhoType;
  TF1* fitLine = histoThresh->GetFunction(funcName.c_str());
  fitLine->SetLineStyle(2);
  fitLine->SetLineWidth(1);

  histoThresh->Draw("Psame");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");


  TPaveText* label_fitResults = new TPaveText(0.3, 0.73, 0.6, 0.85, "brNDC");
  label_fitResults->SetFillColor(0);
  label_fitResults->SetTextSize(0.04);
  char label_fitResultsText[300];
  sprintf( label_fitResultsText, "y = %.3f x + %.3f GeV", fitLine->GetParameter(1), fitLine->GetParameter(0) );
  label_fitResults->AddText(label_fitResultsText);
  label_fitResults->Draw("same");

  std::string canvasName = db->get_outputdir() + "/" + varName + "Thresh_vs_PF";
  std::string canvasName_png = canvasName + ".png";
  std::string canvasName_eps = canvasName + ".eps";

  c1->SaveAs(canvasName_png.c_str());
  c1->SaveAs(canvasName_eps.c_str());

}
