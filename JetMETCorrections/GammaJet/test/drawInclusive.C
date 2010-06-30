#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"

std::string PICTUREFORMAT = "eps";
std::string OUTPUTDIR;
Int_t PT_THRESH;
std::string RAW_CORR;


std::string getAxisName(std::string name);
void drawHisto( TFile* dataFile, TFile* mcFile, std::string name, std::string etaRegion, Float_t scaleFactor, int legendQuadrant=1, bool log_aussi=false);
void drawProfile( TFile* dataFile, TFile* mcFile, std::string yVar, std::string xVar, int legendQuadrant=1);
void shrinkPad(double b=0.1, double l=0.1, double r=0.1, double t=0.1);



void drawInclusive( std::string dataset, std::string mc_version, std::string algoType, Int_t pt_thresh=5, std::string raw_corr="raw" ) {

  PT_THRESH = pt_thresh;
  RAW_CORR = raw_corr;

  TStyle *simpleStyle = new TStyle("simpleStyle","");
  simpleStyle->SetCanvasColor(0);
  simpleStyle->SetPadColor(0);
  simpleStyle->SetFrameFillColor(0);
  simpleStyle->SetStatColor(0);
  simpleStyle->SetOptStat(0);
  simpleStyle->SetTitleFillColor(0);
  simpleStyle->SetCanvasBorderMode(0);
  simpleStyle->SetPadBorderMode(0);
  simpleStyle->SetFrameBorderMode(0);
  simpleStyle->cd();


  char outputdir_char[200];
  if( mc_version=="" )
    sprintf( outputdir_char, "InclusivePlots_%s_%s_%dGeV%s", dataset.c_str(), algoType.c_str(), PT_THRESH, raw_corr.c_str());
  else
    sprintf( outputdir_char, "InclusivePlots_%s_%s_%s_%dGeV%s", dataset.c_str(), mc_version.c_str(), algoType.c_str(), PT_THRESH, raw_corr.c_str());
  std::string outputdir_str(outputdir_char);
  OUTPUTDIR = outputdir_str;

  //std::string dataFileName = "Inclusive_" + dataset + "_" + algoType;
  char dataFileName[150];
  sprintf( dataFileName, "Inclusive_DATA_%s_%s_%dGeV%s.root", dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str());
  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  char mcFileName[150];
  if( mc_version=="" ) {
    sprintf( mcFileName, "Inclusive_MC_%s_%s_%dGeV%s.root", dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str());
  } else {
    sprintf( mcFileName, "Inclusive_MC_%s_%s_%s_%dGeV%s.root", dataset.c_str(), mc_version.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str());
  }
  TFile* mcFile = TFile::Open(mcFileName);
  std::cout << "Opened mc file '" << mcFileName << "'." << std::endl;

  TH1F* h1_phi_data = (TH1F*)dataFile->Get("phiJet");
  TH1F* h1_phi_mc = (TH1F*)mcFile->Get("phiJet");
 

  Float_t scaleFactor = h1_phi_data->Integral()/h1_phi_mc->Integral();


  bool log = true;
  drawHisto( dataFile, mcFile, "massJet", "", scaleFactor, 1, log);

  drawHisto( dataFile, mcFile, "ptJet", "", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "ptCorrJet", "", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "phiJet", "", scaleFactor);
  drawHisto( dataFile, mcFile, "etaJet", "", scaleFactor);

  drawHisto( dataFile, mcFile, "nCandJet", "", scaleFactor, log);

  drawHisto( dataFile, mcFile, "RchJet", "", scaleFactor, 1);
  drawHisto( dataFile, mcFile, "RnhJet", "", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "RgammaJet", "", scaleFactor, 1, log);

  drawHisto( dataFile, mcFile, "EchJet", "", scaleFactor, 1);
  drawHisto( dataFile, mcFile, "EnhJet", "", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "EnhJet", "barrel", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "EnhJet", "endcap", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "EgammaJet", "", scaleFactor, 1, log);

  drawHisto( dataFile, mcFile, "EphotAveJet", "", scaleFactor, 1, log);

  drawHisto( dataFile, mcFile, "NchJet", "", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "NnhJet", "", scaleFactor, 1, log);
  drawHisto( dataFile, mcFile, "NgammaJet", "", scaleFactor, 1, log);

  drawProfile( dataFile, mcFile, "pt", "eta", 1);
  drawProfile( dataFile, mcFile, "ptCorr", "eta", 1);

  drawProfile( dataFile, mcFile, "Rch", "eta", 1);
  drawProfile( dataFile, mcFile, "Rgamma", "eta", 1);
  drawProfile( dataFile, mcFile, "Rnh", "eta", 1);

  drawProfile( dataFile, mcFile, "Ech", "eta", 1);
  drawProfile( dataFile, mcFile, "Egamma", "eta", 1);
  drawProfile( dataFile, mcFile, "Enh", "eta", 1);

  drawProfile( dataFile, mcFile, "Nch", "eta", 1);
  drawProfile( dataFile, mcFile, "Ngamma", "eta", 1);
  drawProfile( dataFile, mcFile, "Nnh", "eta", 1);

  drawProfile( dataFile, mcFile, "Rch", "ptCorr", 1);
  drawProfile( dataFile, mcFile, "Rgamma", "ptCorr", 1);
  drawProfile( dataFile, mcFile, "Rnh", "ptCorr", 1);

  drawProfile( dataFile, mcFile, "Ech", "ptCorr", 1);
  drawProfile( dataFile, mcFile, "Egamma", "ptCorr", 1);
  drawProfile( dataFile, mcFile, "Enh", "ptCorr", 1);

  drawProfile( dataFile, mcFile, "Nch", "ptCorr", 1);
  drawProfile( dataFile, mcFile, "Ngamma", "ptCorr", 1);
  drawProfile( dataFile, mcFile, "Nnh", "ptCorr", 1);

  drawProfile( dataFile, mcFile, "Rch", "pt", 1);
  drawProfile( dataFile, mcFile, "Rgamma", "pt", 1);
  drawProfile( dataFile, mcFile, "Rnh", "pt", 1);

  drawProfile( dataFile, mcFile, "Ech", "pt", 1);
  drawProfile( dataFile, mcFile, "Egamma", "pt", 1);
  drawProfile( dataFile, mcFile, "Enh", "pt", 1);

  drawProfile( dataFile, mcFile, "Nch", "pt", 1);
  drawProfile( dataFile, mcFile, "Ngamma", "pt", 1);
  drawProfile( dataFile, mcFile, "Nnh", "pt", 1);

  drawProfile( dataFile, mcFile, "Egamma", "Ech", 1);
  drawProfile( dataFile, mcFile, "Enh",    "Ech", 1);
  drawProfile( dataFile, mcFile, "Rgamma", "Rch", 1);
  drawProfile( dataFile, mcFile, "Rnh",    "Rch", 1);
  drawProfile( dataFile, mcFile, "Rnh",    "Rgamma", 1);
  drawProfile( dataFile, mcFile, "Ngamma", "Nch", 1);
  drawProfile( dataFile, mcFile, "Nnh",    "Nch", 1);

}  


void drawHisto( TFile* dataFile, TFile* mcFile, std::string name, std::string etaRegion, Float_t scaleFactor, int legendQuadrant, bool log_aussi) {

  std::string histoName = name;
  if( etaRegion!="" ) histoName = name + "_" + etaRegion;

  TH1F* dataHisto = (TH1F*)dataFile->Get(histoName.c_str());
  TH1F* mcHisto = (TH1F*)mcFile->Get(histoName.c_str());

  mcHisto->SetFillColor(38);

  dataHisto->SetMarkerStyle(20);

  Float_t yAxisMaxScale = (name=="phiJet" || name=="etaJet") ? 1.8 : 1.6;
  Float_t xMin = dataHisto->GetXaxis()->GetXmin();
  Float_t xMax = dataHisto->GetXaxis()->GetXmax();
  Float_t yMax_data = dataHisto->GetMaximum();
  Float_t yMax_mc = dataHisto->GetMaximum();
  Float_t yMax = (yMax_data>yMax_mc) ? yAxisMaxScale*yMax_data : yAxisMaxScale*yMax_mc;
  Float_t yMin = 0.;

  if( name=="NchJet" || name=="NnhJet" || name=="NgammaJet" ) xMax = 15.5;

  std::string xAxis = getAxisName(name);

  std::string yAxis = "Jets";
  if( name=="ptJet" || name=="ptCorrJet" ) {
    char yAxis_char[50];
    sprintf(yAxis_char, "Jets/(%d GeV/c)", (Int_t)dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  


  if( name=="diJetMass" ) {
    char yAxis_char[50];
    sprintf(yAxis_char, "Events/(%d GeV/c^{2})", (Int_t)dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  
  if( name=="massJet" ) {
    char yAxis_char[50];
    sprintf(yAxis_char, "Jets/(%d GeV/c^{2})", (Int_t)dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  
  if( name=="deltaPhiJet" || name=="asymmJet" ) yAxis = "Events";

  std::string etaRange = "|#eta| < 2.4";
  if( etaRegion=="barrel") etaRange = "|#eta| < 1.4";
  if( etaRegion=="endcap") etaRange = "1.4 < |#eta| < 2.4";


  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax);
  h2_axes->SetXTitle(xAxis.c_str());
  h2_axes->SetYTitle(yAxis.c_str());
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);

  if( legendQuadrant<0 || legendQuadrant>5 ) {
    std::cout << "Invalid legendQuadrant! Exiting!" << std::endl;
    exit(1);
  }

  Float_t legend_xMin;
  Float_t legend_yMin;
  Float_t legend_xMax;
  Float_t legend_yMax;

  if( legendQuadrant==1 ) {
    legend_xMin = 0.65;
    legend_yMin = 0.7;
    legend_xMax = 0.88;
    legend_yMax = 0.88;
  } else if( legendQuadrant==0 ) {
    legend_xMin = 0.47;
    legend_yMin = 0.7;
    legend_xMax = 0.7;
    legend_yMax = 0.88;
  } else if( legendQuadrant==2 ) {
    legend_xMin = 0.15;
    legend_yMin = 0.7;
    legend_xMax = 0.4;
    legend_yMax = 0.88;
  } else if( legendQuadrant==5 ) {
    legend_xMin = 0.4;
    legend_yMin = 0.15;
    legend_xMax = 0.6;
    legend_yMax = 0.25;
  }


  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax);
  legend->SetFillColor(kWhite);
  legend->AddEntry(dataHisto, "Data", "P");
  legend->AddEntry(mcHisto, "MC", "F");

  TPaveText* label_cms = new TPaveText(0.15, 0.84, 0.35, 0.88, "brNDC");
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.038);
  label_cms->SetTextFont(62);
  label_cms->AddText("CMS Preliminary");

  TPaveText* label_sqrt = new TPaveText(0.15, 0.78, 0.35, 0.82, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->AddText("#sqrt{s} = 900 GeV");

  Float_t label_algo_xMin = 0.63;
  Float_t label_algo_yMin = 0.5;
  Float_t label_algo_xMax = 0.84;
  Float_t label_algo_yMax = 0.65;

  if( name=="asymmJet" || name=="deltaPhiJet" || name=="RchJet" ) {
   label_algo_xMin = 0.15;
   label_algo_yMin = 0.55;
   label_algo_xMax = 0.36;
   label_algo_yMax = 0.7;
  }
  if( name=="phiJet" || name=="etaJet" ) {
   label_algo_xMin = 0.4;
   label_algo_yMin = 0.6;
   label_algo_xMax = 0.6;
   label_algo_yMax = 0.75;
  }

  TPaveText* label_algo = new TPaveText(label_algo_xMin, label_algo_yMin, label_algo_xMax, label_algo_yMax,  "brNDC");
  label_algo->SetFillColor(kWhite);
  label_algo->SetTextSize(0.035);
  label_algo->SetTextFont(42);
  label_algo->AddText("Anti-k_{T} R=0.5");
  if( name != "etaJet" )
    label_algo->AddText(etaRange.c_str());
  if( name != "ptJet" && name != "ptCorrJet" ) {
    char labelText[70];
    sprintf( labelText, "p_{T}^{%s} > %d GeV/c", RAW_CORR.c_str(), PT_THRESH);
    label_algo->AddText(labelText);
  }

  mcHisto->Scale(scaleFactor);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

/*
  TPad* pad_hi = new TPad("pad_hi", "", 0., 0.3, 0.99, 0.99);
  pad_hi->Draw();
  TPad* pad_lo = new TPad("pad_lo", "", 0., 0., 0.99, 0.29);
  pad_lo->Draw();
  
  pad_hi->cd();
  shrinkPad(0.001, 0.2);
*/


  h2_axes->Draw("");
  mcHisto->Draw("same");
  dataHisto->Draw("E same");
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

/*
  pad_lo->cd();
  shrinkPad(0.001, 0.2);
  
  TH1F* h1_ratio = new TH1F(*dataHisto);
  h1_ratio->SetName("ratio");
  h1_ratio->Reset();
  h1_ratio->Divide( &(*dataHisto), &(*mcHisto) );
  h1_ratio->SetMinimum(0.5);
  h1_ratio->SetMaximum(1.5);
  h1_ratio->SetMarkerStyle(20);
  h1_ratio->SetMarkerSize(2.);
  h1_ratio->SetMarkerColor(38);
  h1_ratio->SetNdivisions(504, "Y");
  h1_ratio->SetLabelSize(0.22, "X");
  h1_ratio->SetLabelSize(0.06, "Y");
  h1_ratio->SetStats(0);
  h1_ratio->SetTitle("");
  h1_ratio->Draw("P");
*/



  std::string canvasName = OUTPUTDIR + "/" + name;
  if( etaRegion!="" )
    canvasName = canvasName + "_" + etaRegion;
  std::string canvasName_eps = canvasName + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  std::string canvasName_png = canvasName + ".png";
  c1->SaveAs(canvasName_png.c_str());

  if( log_aussi ) {
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.5, 5.*yMax);
    h2_axes_log->SetXTitle(xAxis.c_str());
    h2_axes_log->SetYTitle(yAxis.c_str());
    h2_axes_log->GetXaxis()->SetTitleOffset(1.1);
    h2_axes_log->GetYaxis()->SetTitleOffset(1.2);
    c1->SetLogy();
    h2_axes_log->Draw("");
    mcHisto->Draw("same");
    dataHisto->Draw("E same");
    legend->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");
    canvasName = OUTPUTDIR + "/" + name;
    if( etaRegion!="" )
      canvasName = canvasName + "_" + etaRegion;
    canvasName = canvasName + "_log";
    canvasName_eps = canvasName + ".eps";
    c1->SaveAs(canvasName_eps.c_str());
    canvasName_png = canvasName + ".png";
    c1->SaveAs(canvasName_png.c_str());
    delete h2_axes_log;
  }
    
  delete c1;
  delete legend;
  delete h2_axes;

}


void drawProfile( TFile* dataFile, TFile* mcFile, std::string yVar, std::string xVar, int legendQuadrant) {

  std::string name = yVar + "_vs_" + xVar;

  TProfile* dataProfile = (TProfile*)dataFile->Get(name.c_str());
  TProfile* mcProfile = (TProfile*)mcFile->Get(name.c_str());

  Float_t profile_xMin = dataProfile->GetXaxis()->GetXmin();
  Float_t profile_xMax = dataProfile->GetXaxis()->GetXmax();

  mcProfile->SetFillColor(38);

  dataProfile->SetMarkerStyle(20);

  Float_t xMin = (xVar=="eta") ? -2.5 : profile_xMin;
  Float_t xMax = (xVar=="eta") ?  2.5 : profile_xMax;

  if( xVar=="Nch" ) xMax = 11.5;
  
  Float_t dataMax = dataProfile->GetMaximum();
  Float_t mcMax = mcProfile->GetMaximum();
  Float_t plotMax = (dataMax>mcMax) ? dataMax : mcMax;

  std::string xAxisName = getAxisName(xVar);
  std::string yAxisName = getAxisName(yVar);

  Float_t yAxisMaxScale = 1.5;
  if( yVar=="pt" || yVar=="ptCorr" || yVar=="Rch" || yVar=="Nch" || yVar=="Ngamma" || yVar=="Nnh" ) yAxisMaxScale=1.8;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yAxisMaxScale*plotMax);
  h2_axes->SetXTitle(xAxisName.c_str());
  h2_axes->SetYTitle(yAxisName.c_str());
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);


  if( legendQuadrant<0 || legendQuadrant>5 ) {
    std::cout << "Invalid legendQuadrant! Exiting!" << std::endl;
    exit(1);
  }

  Float_t legend_xMin;
  Float_t legend_yMin;
  Float_t legend_xMax;
  Float_t legend_yMax;

  if( legendQuadrant==1 ) {
    legend_xMin = 0.65;
    legend_yMin = 0.7;
    legend_xMax = 0.88;
    legend_yMax = 0.88;
  } else if( legendQuadrant==2 ) {
    legend_xMin = 0.15;
    legend_yMin = 0.7;
    legend_xMax = 0.4;
    legend_yMax = 0.88;
  } else if( legendQuadrant==5 ) {
    legend_xMin = 0.4;
    legend_yMin = 0.15;
    legend_xMax = 0.6;
    legend_yMax = 0.25;
  }


  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax);
  legend->SetFillColor(kWhite);
  legend->AddEntry(dataProfile, "Data", "P");
  legend->AddEntry(mcProfile, "MC", "F");

  TPaveText* label_cms = new TPaveText(0.15, 0.84, 0.35, 0.88, "brNDC");
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.038);
  label_cms->SetTextFont(62);
  label_cms->AddText("CMS Preliminary");

  TPaveText* label_sqrt = new TPaveText(0.15, 0.78, 0.35, 0.82, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->AddText("#sqrt{s} = 900 GeV");

  Float_t label_algo_xMin = 0.4;
  Float_t label_algo_yMin = 0.55;
  Float_t label_algo_xMax = 0.6;
  Float_t label_algo_yMax = 0.7;

  if( yVar=="pt" || yVar=="ptCorr" || yVar=="Rch" || yVar=="Rgamma" ) {
    label_algo_xMin = 0.4;
    label_algo_yMin = 0.65;
    label_algo_xMax = 0.6;
    label_algo_yMax = 0.8;
  }

  if( xVar=="pt" || xVar=="ptCorr" ) {
    label_algo_xMin = 0.15;
    label_algo_xMax = 0.35;
  }



  TPaveText* label_algo = new TPaveText(label_algo_xMin, label_algo_yMin, label_algo_xMax, label_algo_yMax, "brNDC");
  label_algo->SetFillColor(kWhite);
  label_algo->SetTextSize(0.035);
  label_algo->SetTextFont(42);
  label_algo->AddText("Anti-k_{T} R=0.5");
  if( xVar != "eta" )
    label_algo->AddText("|#eta| < 2.4");
  if( yVar != "pt" && yVar != "ptCorr" ) {
    char labelText[70];
    sprintf( labelText, "p_{T}^{%s} > %d GeV/c", RAW_CORR.c_str(), PT_THRESH);
    label_algo->AddText(labelText);
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  //if(name=="ptJet") c1->SetLogy();
  h2_axes->Draw("");
  mcProfile->Draw("histo same");
  dataProfile->Draw("E same");
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");
 
  std::string canvasName = OUTPUTDIR + "/" + name;
  std::string canvasName_eps = canvasName + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  std::string canvasName_png = canvasName + ".png";
  c1->SaveAs(canvasName_png.c_str());


  delete c1;
  delete legend;
  delete h2_axes;

}

std::string getAxisName(std::string name) { 
  
  std::string axisName = "";

  if( name=="ptJet"|| name=="pt" ) axisName = "Raw p_{T} [GeV/c]";
  if( name=="ptCorrJet"|| name=="ptCorr" ) axisName = "Corrected p_{T} [GeV/c]";
  if( name=="etaJet"|| name=="eta" ) axisName = "#eta";
  if( name=="phiJet"|| name=="phi" ) axisName = "#phi";
  if( name=="nCandJet" ) axisName = "N Candidates";
  if( name=="RchJet"|| name=="Rch" ) axisName = "R_{ch}";
  if( name=="RnhJet"|| name=="Rnh" ) axisName = "R_{nh}";
  if( name=="RgammaJet"|| name=="Rgamma" ) axisName = "R_{#gamma}";
  if( name=="EchJet"|| name=="Ech" ) axisName = "E_{ch} [GeV]";
  if( name=="EnhJet"|| name=="Enh" ) axisName = "E_{nh} [GeV]";
  if( name=="EgammaJet"|| name=="Egamma" ) axisName = "E_{#gamma} [GeV]";
  if( name=="asymmJet" ) axisName = "Asymmetry";
  if( name=="deltaPhiJet" ) axisName = "#Delta#Phi";
  if( name=="massJet" ) axisName = "Jet Invariant Mass [GeV/c^{2}]";
  if( name=="diJetMass" ) axisName = "DiJet Invariant Mass [GeV/c^{2}]";
  if( name=="EphotAveJet" ) axisName = "Average Photon Energy in Jet [GeV]";
  if( name=="NchJet" || name=="Nch" ) axisName = "Number of Charged Hadrons in Jet";
  if( name=="NgammaJet" || name=="Ngamma" ) axisName = "Number of Photons in Jet";
  if( name=="NnhJet" || name=="Nnh" ) axisName = "Number of Neutral Hadrons in Jet";

  return axisName;

}


void shrinkPad(double b, double l, double r, double t) {
 gPad->SetBottomMargin(b); 
 gPad->SetLeftMargin(l);
 gPad->SetRightMargin(r);
 gPad->SetTopMargin(t);
}
