#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.C"



void drawSinglePtBin(const std::string& fileName, int iPtBin, const std::string& recoGen);
void fillFromHistos( Int_t n, Float_t *y, Float_t *y_err, const std::vector< TH1D* >& h1_r);



void draw2ndJet(const std::string& fileName) {

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


  

  int nBins_pt = 23;
   
  for( int i=0; i<(nBins_pt-1); ++i) {
    drawSinglePtBin(fileName, i, "reco");
    drawSinglePtBin(fileName, i, "gen");
    drawSinglePtBin(fileName, i, "recoRel");
    drawSinglePtBin(fileName, i, "genRel");
  }

}




void drawSinglePtBin(const std::string& fileName, int iPtBin, const std::string& x_axis) {

  int nBins_pt = 23;

  Double_t Lower[nBins_pt];
  Lower[0] = 20.;
  Lower[1] = 25.;
  Lower[2] = 30.;
  Lower[3] = 36.;
  Lower[4] = 42.;
  Lower[5] = 48.;
  Lower[6] = 55.;
  Lower[7] = 62.;
  Lower[8] = 70.;
  Lower[9] = 78.;
  Lower[10] = 87.;
  Lower[11] = 97.;
  Lower[12] = 108.;
  Lower[13] = 121.;
  Lower[14] = 135.;
  Lower[15] = 152.;
  Lower[16] = 172.;
  Lower[17] = 197.;
  Lower[18] = 231.;
  Lower[19] = 283.;
  Lower[20] = 397.;
  Lower[21] = 600.;
  Lower[22] = 1000.;

  Double_t ptMin = Lower[iPtBin];
  Double_t ptMax = Lower[iPtBin+1];

  TFile* file = TFile::Open(fileName.c_str());
 
  std::vector< TH1D* > h1_recoPhot;
  std::vector< TH1D* > h1_recoGen;
  std::vector< TH1D* > h1_genPhot;
  std::vector< TH1D* > h1_pt2ndJetMean;

  Int_t nPoints = 5;


  for( int i=0; i<nPoints; ++i ) {
    char name[100];
    sprintf( name, "ptBin_%d_%d/r_recoPhot_vs_%s_%d_%d", (int)ptMin, (int)ptMax, x_axis.c_str(), iPtBin, i);
    h1_recoPhot.push_back((TH1D*)file->Get(name));
    sprintf( name, "ptBin_%d_%d/r_recoGen_vs_%s_%d_%d", (int)ptMin, (int)ptMax, x_axis.c_str(), iPtBin, i);
    h1_recoGen.push_back((TH1D*)file->Get(name));
    sprintf( name, "ptBin_%d_%d/r_genPhot_vs_%s_%d_%d", (int)ptMin, (int)ptMax, x_axis.c_str(), iPtBin, i);
    h1_genPhot.push_back((TH1D*)file->Get(name));
    if( x_axis=="reco" )
      sprintf( name, "ptBin_%d_%d/pt2ndJetRecoMean_%d_%d", (int)ptMin, (int)ptMax, iPtBin, i);
    else if ( x_axis=="gen" )
      sprintf( name, "ptBin_%d_%d/pt2ndJetGenMean_%d_%d", (int)ptMin, (int)ptMax, iPtBin, i);
    else if ( x_axis=="recoRel" )
      sprintf( name, "ptBin_%d_%d/pt2ndJetRecoRelMean_%d_%d", (int)ptMin, (int)ptMax, iPtBin, i);
    else if ( x_axis=="genRel" )
      sprintf( name, "ptBin_%d_%d/pt2ndJetGenRelMean_%d_%d", (int)ptMin, (int)ptMax, iPtBin, i);
    h1_pt2ndJetMean.push_back((TH1D*)file->Get(name));
  }

  Float_t x[nPoints];
  Float_t x_err[nPoints];

  Float_t y_recoPhot[nPoints];
  Float_t y_recoPhot_err[nPoints];

  Float_t y_recoGen[nPoints];
  Float_t y_recoGen_err[nPoints];

  Float_t y_genPhot[nPoints];
  Float_t y_genPhot_err[nPoints];

  for( int i =0; i<nPoints; ++i ) {
    x[i] = h1_pt2ndJetMean[i]->GetMean();
    x_err[i] =  h1_pt2ndJetMean[i]->GetRMS()/sqrt((Float_t)h1_pt2ndJetMean[i]->GetEntries());
  }


  fillFromHistos( nPoints, y_recoPhot, y_recoPhot_err, h1_recoPhot);
  fillFromHistos( nPoints, y_recoGen, y_recoGen_err, h1_recoGen);
  fillFromHistos( nPoints, y_genPhot, y_genPhot_err, h1_genPhot);


  TGraphErrors* gr_recoPhot = new TGraphErrors(nPoints, x, y_recoPhot, x_err, y_recoPhot_err);
  gr_recoPhot->SetMarkerStyle(20);
  gr_recoPhot->SetMarkerColor(kRed);

  TGraphErrors* gr_recoGen = new TGraphErrors(nPoints, x, y_recoGen, x_err, y_recoGen_err);
  gr_recoGen->SetMarkerStyle(21);
  gr_recoGen->SetMarkerColor(kBlue);

  TGraphErrors* gr_genPhot = new TGraphErrors(nPoints, x, y_genPhot, x_err, y_genPhot_err);
  gr_genPhot->SetMarkerStyle(22);
  gr_genPhot->SetMarkerColor(kBlack);

  TF1* line1 = new TF1("line1", "[0] + x*[1]");
  line1->SetRange(0., 30.);
  Double_t x1, x2, y1, y2;
  Int_t ii;
  ii = gr_genPhot->GetPoint(1, x1, y1);
  ii = gr_genPhot->GetPoint(2, x2, y2);
  line1->SetParameter(0, 0.);
  line1->SetParameter(1, (y2-y1)/(x2-x1));
  line1->SetLineWidth(0.5);
  gr_genPhot->Fit(line1, "RQ");

  TF1* line2 = new TF1("line2", "[0] + x*[1]");
  line2->SetRange(0., 30.);
  ii = gr_recoGen->GetPoint(1, x1, y1);
  ii = gr_recoGen->GetPoint(2, x2, y2);
  line2->SetParameter(0, y1);
  line2->SetParameter(1, (y2-y1)/(x2-x1));
  line2->SetLineWidth(0.5);
  line2->SetLineColor(kBlue);
  gr_recoGen->Fit(line2, "RQ");

  TF1* sum = new TF1("sum", "sqrt( line1*line1 + line2*line2 )");
  sum->SetRange(0., 30.);
  //sum->SetLineWidth(1.);
  sum->SetLineColor(kGray+2);

  TF1* fit_recoPhot = new TF1("fit_recoPhot", "sqrt([0]+[1]*x*x)");
  fit_recoPhot->SetRange(0., 30.);
  fit_recoPhot->SetParameter(0, y1*y1);
  fit_recoPhot->SetLineColor(kRed);
  fit_recoPhot->SetLineWidth(0.5);
  gr_recoPhot->Fit(fit_recoPhot, "RQ");

  Float_t ymax=0.5;
  if( iPtBin>12 )
    ymax = 0.25;


  std::string xTitle;
  if( x_axis=="reco" )
    xTitle = "p_{T}^{2ndJet reco} [GeV/c]";
  else if( x_axis=="gen" )
    xTitle = "p_{T}^{2ndJet gen} [GeV/c]";
  else if( x_axis=="recoRel" )
    xTitle = "p_{T}^{2ndJet reco}/p_{T}^{#gamma} [%]";
  else if( x_axis=="genRel" )
    xTitle = "p_{T}^{2ndJet gen}/p_{T}^{#gamma} [%]";

  TH2D* h2_aux = new TH2D("aux", "", 10, 0., 35., 10, 0., ymax);
  h2_aux->SetXTitle(xTitle.c_str());
  h2_aux->GetXaxis()->SetTitleOffset(1.1);
  h2_aux->SetYTitle("resolution [%]");
  h2_aux->GetYaxis()->SetTitleOffset(1.2);

  TLegend* legend = new TLegend(0.2, 0.6, 0.45, 0.85);
  legend->SetTextSize(0.04);
  legend->SetFillStyle(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(gr_recoGen, "Intrinsic", "P");
  legend->AddEntry(gr_genPhot, "FSR", "P");
  legend->AddEntry(gr_recoPhot, "#gamma + jet", "P");
  legend->AddEntry(sum, "Intr. #oplus FSR", "L");

  char labeltext[50];
  sprintf(labeltext, "%d < p_{T}^{#gamma} < %d GeV/c", (int)ptMin, (int)ptMax);  


  TPaveText* label = new TPaveText(0.67, 0.83, 0.8, 0.86, "brNDC");
  label->SetFillColor(kWhite);
  label->SetTextSize(0.035);
  label->AddText(labeltext);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h2_aux->Draw();
  sum->Draw("same");
  gr_recoPhot->Draw("Psame");
  gr_recoGen->Draw("Psame");
  gr_genPhot->Draw("Psame");
  legend->Draw("same");
  label->Draw("same");


  char canvasName[100];
  sprintf(canvasName, "secondJetPlots/ptBin_%d_%d_%s.png", (int)ptMin, (int)ptMax, x_axis.c_str());
  c1->SaveAs(canvasName);
  delete c1;

  file->Close();

}




void fillFromHistos( Int_t n, Float_t *y, Float_t *y_err, const std::vector< TH1D* >& h1_r) {


  for( int i = 0; i<n; ++i ) {

    //this ugly fix saves from empty relative pt binning plots
    if( h1_r[i]->GetEntries()==0 ) {
      y[i]=-1.;
      y_err[i]=0.;
      continue;
    }

    TF1* gaussian = h1_r[i]->GetFunction("gaussian");

    Float_t mu, sigma, mu_err, sigma_err;
    mu = gaussian->GetParameter(1);
    sigma = gaussian->GetParameter(2);

    TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
    fitTools::fitProjection(h1_r[i], gaussian_chi, "RQN");
    mu_err = gaussian_chi->GetParError(1);
    sigma_err = gaussian_chi->GetParError(2);

    y[i] = sigma/mu;
    y_err[i] = sqrt( sigma_err*sigma_err/(mu*mu) + sigma*sigma*mu_err*mu_err/(mu*mu*mu*mu) ); 
  
  }

}
