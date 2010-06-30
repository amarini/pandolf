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

std::string OUTPUTDIR_;
std::string FIT_RMS = "FIT";
std::string PICTURE_FORMAT = "eps";
bool ONLYBLACK = true;
bool NOQ = false;
Int_t LUMI = 0;
Int_t INT_PERC = 95;
int NBINS_PT = 16;
std::string RECOTYPE_;
std::string ALGOTYPE_;


//void drawSinglePtBin(TFile* file, int iPtBin, const std::string& recoGen, Float_t& extrapReso, Float_t& extrapReso_err, Float_t& trueReso, Float_t& trueReso_err, Float_t& extrapResp, Float_t& extrapResp_err, Float_t& trueResp, Float_t& trueResp_err);
void drawSinglePtBin(TFile* file, int iPtBin, const std::string& recoGen, Float_t& extrapReso, Float_t& extrapReso_err, Float_t& trueReso, Float_t& trueReso_err, Float_t& intrReso, Float_t& intrReso_err, Float_t& extrapResp, Float_t& extrapResp_err, Float_t& trueResp, Float_t& trueResp_err, Float_t& intrResp, Float_t& intrResp_err, Float_t& imbalanceResp, Float_t& imbalanceResp_err);
void getXPoints( TFile * file, const char* xHistoName, Int_t nPoints, Float_t* x, Float_t* x_err);
void getYPoints( TFile * file, const char* yHistoName, Int_t nPoints, Float_t* y_resp, Float_t* y_resp_err, Float_t* y_reso, Float_t* y_reso_err);
//void drawGraphs_vs_pt( const std::string& varY, const std::string& canvasName, TGraphErrors* gr_resoReco_cutOn2ndJet, TGraphErrors* gr_trueResoGen_vs_pt, TGraphErrors* gr_extrapReso_vs_pt=0 );
void drawGraphs_vs_pt( const std::string& varY, const std::string& canvasName, TGraphErrors* gr_resoReco_cutOn2ndJet, TGraphErrors* gr_intrReso_vs_pt, TGraphErrors* gr_extrapReso_vs_pt=0, TGraphErrors* gr_trueReso_vs_pt=0);
void drawPullHistogram( const std::string& name, Float_t* pull, Float_t* pull_err);



void drawSecondJet(std::string dataset, std::string recoType, std::string jetAlgo, const std::string& fitrms="") {

  ALGOTYPE_ = (recoType=="pf") ? recoType+jetAlgo : jetAlgo;
  RECOTYPE_ = recoType;

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

  if( fitrms != "" ) FIT_RMS = fitrms;

  char suffix_char[200];
  sprintf(suffix_char, "%s_%s_%s_%d", dataset.c_str(), ALGOTYPE_.c_str(), FIT_RMS.c_str(), INT_PERC);
  std::string suffix(suffix_char);
  if( LUMI>0 ) {
    char lumi_char[20];
    sprintf( lumi_char, "_%dPB", LUMI);
    std::string lumi_str(lumi_char);
    suffix += lumi_str;
  }
  if( NOQ ) suffix += "_05Q";
  OUTPUTDIR_ = "SecondJetPlots_" + suffix;

  std::string fileName = "SecondJet_"+dataset + "_" + ALGOTYPE_ + ".root";
  TFile* file = TFile::Open(fileName.c_str());

  TH2D* h2_ptPhotReco_vs_pt = (TH2D*)file->Get("ptPhotMean_vs_pt");


  Float_t ptPhotReco[NBINS_PT-1];
  Float_t ptPhotReco_err[NBINS_PT-1];

  Float_t extrapReso_RecoRel[NBINS_PT-1];
  Float_t extrapReso_err_RecoRel[NBINS_PT-1];
  Float_t trueReso_RecoRel[NBINS_PT-1];
  Float_t trueReso_err_RecoRel[NBINS_PT-1];
  Float_t intrReso_RecoRel[NBINS_PT-1];
  Float_t intrReso_err_RecoRel[NBINS_PT-1];
  Float_t extrapResp_RecoRel[NBINS_PT-1];
  Float_t extrapResp_err_RecoRel[NBINS_PT-1];
  Float_t trueResp_RecoRel[NBINS_PT-1];
  Float_t trueResp_err_RecoRel[NBINS_PT-1];
  Float_t intrResp_RecoRel[NBINS_PT-1];
  Float_t intrResp_err_RecoRel[NBINS_PT-1];
  Float_t imbalanceResp_RecoRel[NBINS_PT-1];
  Float_t imbalanceResp_err_RecoRel[NBINS_PT-1];
  Float_t pullResp_RecoRel[NBINS_PT-1];
  Float_t pullResp_err_RecoRel[NBINS_PT-1];
  Float_t pullReso_RecoRel[NBINS_PT-1];
  Float_t pullReso_err_RecoRel[NBINS_PT-1];
   
   

  for( int i=0; i<(NBINS_PT-1); ++i) {
    char projName[50];
    sprintf(projName, "projection_%d",i);
    TH1D* h1_proj = h2_ptPhotReco_vs_pt->ProjectionY(projName, i+1, i+1);
    ptPhotReco[i] = h1_proj->GetMean();
    ptPhotReco_err[i] = (h1_proj->GetEntries()>1.) ? h1_proj->GetRMS()/sqrt(h1_proj->GetEntries()) : h1_proj->GetRMS();
    drawSinglePtBin(file, i, "RecoRel", extrapReso_RecoRel[i], extrapReso_err_RecoRel[i], trueReso_RecoRel[i], trueReso_err_RecoRel[i], intrReso_RecoRel[i], intrReso_err_RecoRel[i], 
                                        extrapResp_RecoRel[i], extrapResp_err_RecoRel[i], trueResp_RecoRel[i], trueResp_err_RecoRel[i], intrResp_RecoRel[i], intrResp_err_RecoRel[i], imbalanceResp_RecoRel[i], imbalanceResp_err_RecoRel[i]);

    pullResp_RecoRel[i] = 100.*(extrapResp_RecoRel[i] - intrResp_RecoRel[i])/intrResp_RecoRel[i];
    pullResp_err_RecoRel[i] = 100.*sqrt( extrapResp_err_RecoRel[i]*extrapResp_err_RecoRel[i] + intrResp_RecoRel[i]*intrResp_RecoRel[i] );

    pullReso_RecoRel[i] = 100.*(extrapReso_RecoRel[i] - intrReso_RecoRel[i])/intrReso_RecoRel[i];
    pullReso_err_RecoRel[i] = 100.*sqrt( extrapReso_err_RecoRel[i]*extrapReso_err_RecoRel[i] + intrReso_RecoRel[i]*intrReso_RecoRel[i] );

    imbalanceResp_RecoRel[i] *= 100.; //in percent
    imbalanceResp_err_RecoRel[i] *= 100.;

  }

  TGraphErrors* gr_extrapRespRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, extrapResp_RecoRel, ptPhotReco_err, extrapResp_err_RecoRel); 
  gr_extrapRespRecoRel_vs_pt->SetName("gr_extrapRespRecoRel_vs_pt");
  gr_extrapRespRecoRel_vs_pt->SetMarkerStyle(25);
  gr_extrapRespRecoRel_vs_pt->SetMarkerColor(kBlack);
  gr_extrapRespRecoRel_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_intrRespRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, intrResp_RecoRel, ptPhotReco_err, intrResp_err_RecoRel); 
  gr_intrRespRecoRel_vs_pt->SetName("gr_intrRespRecoRel_vs_pt");
  gr_intrRespRecoRel_vs_pt->SetMarkerStyle(29);
  gr_intrRespRecoRel_vs_pt->SetMarkerColor(kBlue);
  gr_intrRespRecoRel_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_trueRespRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, trueResp_RecoRel, ptPhotReco_err, trueResp_err_RecoRel); 
  gr_trueRespRecoRel_vs_pt->SetName("gr_trueRespRecoRel_vs_pt");
  gr_trueRespRecoRel_vs_pt->SetMarkerStyle(23);
  gr_trueRespRecoRel_vs_pt->SetMarkerColor(kGray+2);
  gr_trueRespRecoRel_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_extrapResoRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, extrapReso_RecoRel, ptPhotReco_err, extrapReso_err_RecoRel); 
  gr_extrapResoRecoRel_vs_pt->SetName("gr_extrapResoRecoRel_vs_pt");
  gr_extrapResoRecoRel_vs_pt->SetMarkerStyle(25);
  gr_extrapResoRecoRel_vs_pt->SetMarkerColor(kBlack);
  gr_extrapResoRecoRel_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_intrResoRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, intrReso_RecoRel, ptPhotReco_err, intrReso_err_RecoRel); 
  gr_intrResoRecoRel_vs_pt->SetName("gr_intrResoRecoRel_vs_pt");
  gr_intrResoRecoRel_vs_pt->SetMarkerStyle(29);
  gr_intrResoRecoRel_vs_pt->SetMarkerColor(kBlue);
  gr_intrResoRecoRel_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_trueResoRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, trueReso_RecoRel, ptPhotReco_err, trueReso_err_RecoRel); 
  gr_trueResoRecoRel_vs_pt->SetName("gr_trueResoRecoRel_vs_pt");
  gr_trueResoRecoRel_vs_pt->SetMarkerStyle(23);
  gr_trueResoRecoRel_vs_pt->SetMarkerColor(kGray+2);
  gr_trueResoRecoRel_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_imbalanceRecoRel_vs_pt = new TGraphErrors(NBINS_PT-1, ptPhotReco, imbalanceResp_RecoRel, ptPhotReco_err, imbalanceResp_err_RecoRel); 
  gr_imbalanceRecoRel_vs_pt->SetName("gr_imbalanceRecoRel_vs_pt");
  gr_imbalanceRecoRel_vs_pt->SetMarkerStyle(20);
  gr_imbalanceRecoRel_vs_pt->SetMarkerColor(kBlack);
  gr_imbalanceRecoRel_vs_pt->SetMarkerSize(1.5);

  drawPullHistogram("pullResp", pullResp_RecoRel, pullResp_err_RecoRel);
  drawPullHistogram("pullReso", pullReso_RecoRel, pullReso_err_RecoRel);

  std::string graphFileName= "SecondJet_GraphFile_" + suffix + ".root";
  TFile* graphFile = TFile::Open(graphFileName.c_str(), "RECREATE");
  graphFile->cd();
  gr_extrapRespRecoRel_vs_pt->Write();
  gr_trueRespRecoRel_vs_pt->Write();
  gr_intrRespRecoRel_vs_pt->Write();
  gr_extrapResoRecoRel_vs_pt->Write();
  gr_trueResoRecoRel_vs_pt->Write();
  gr_intrResoRecoRel_vs_pt->Write();
  graphFile->Close();

  std::string graphName;
  std::string canvasName;

  //response:
  graphName = "gr_response_vs_pt_MEAN_barrel_" + ALGOTYPE_;
  TGraphErrors* gr_respReco_cutOn2ndJet = (TGraphErrors*)file->Get(graphName.c_str());
  gr_respReco_cutOn2ndJet->SetMarkerStyle(20);
  gr_respReco_cutOn2ndJet->SetMarkerColor(kRed);
  gr_respReco_cutOn2ndJet->SetMarkerSize(1.5);

  canvasName = OUTPUTDIR_ + "/response_vs_pt_noextrap."+PICTURE_FORMAT;
  drawGraphs_vs_pt( "Response", canvasName, gr_respReco_cutOn2ndJet, gr_intrRespRecoRel_vs_pt );
  
  canvasName = OUTPUTDIR_ + "/response_vs_pt_RECO."+PICTURE_FORMAT;
  drawGraphs_vs_pt( "Response", canvasName, gr_respReco_cutOn2ndJet, gr_intrRespRecoRel_vs_pt, gr_extrapRespRecoRel_vs_pt);


  // resolution:
  std::string resoType = (FIT_RMS=="RMS") ? "RMS" : "FIT";
  graphName = "gr_responseRes_vs_pt_"+resoType+"_barrel_" + ALGOTYPE_;
  TGraphErrors* gr_resoReco_cutOn2ndJet = (TGraphErrors*)file->Get(graphName.c_str());
  gr_resoReco_cutOn2ndJet->SetMarkerStyle(20);
  gr_resoReco_cutOn2ndJet->SetMarkerColor(kRed);
  gr_resoReco_cutOn2ndJet->SetMarkerSize(1.5);

  canvasName = OUTPUTDIR_ + "/resolution_vs_pt_noextrap."+PICTURE_FORMAT;
  drawGraphs_vs_pt( "Resolution", canvasName, gr_resoReco_cutOn2ndJet, gr_intrResoRecoRel_vs_pt );
  
  canvasName = OUTPUTDIR_ + "/resolution_vs_pt_RECO."+PICTURE_FORMAT;
  drawGraphs_vs_pt( "Resolution", canvasName, gr_resoReco_cutOn2ndJet, gr_intrResoRecoRel_vs_pt, gr_extrapResoRecoRel_vs_pt);


  TH2D* h2_axes = new TH2D("axes", "", 10, 20., 2200., 10, -10., 5.);
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);
  h2_axes->SetYTitle("Residual Imbalance [%]");
  h2_axes->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();

  canvasName = OUTPUTDIR_ + "/imbalance." + PICTURE_FORMAT;
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogx();
  c1->SetGridx();
  c1->SetGridy();
  h2_axes->Draw();
  gr_imbalanceRecoRel_vs_pt->Draw("Psame");
  c1->SaveAs(canvasName.c_str());
  


  delete h2_axes;
  delete gr_extrapRespRecoRel_vs_pt;
  delete gr_extrapResoRecoRel_vs_pt;

}




void drawSinglePtBin(TFile* file, int iPtBin, const std::string& recoGen, Float_t& extrapReso, Float_t& extrapReso_err, Float_t& trueReso, Float_t& trueReso_err, Float_t& intrReso, Float_t& intrReso_err, Float_t& extrapResp, Float_t& extrapResp_err, Float_t& trueResp, Float_t& trueResp_err, Float_t& intrResp, Float_t& intrResp_err, Float_t& imbalanceResp, Float_t& imbalanceResp_err) {

  Double_t Lower[NBINS_PT];
  fitTools::getPtBins(NBINS_PT, Lower);

  Double_t ptMin = Lower[iPtBin];
  Double_t ptMax = Lower[iPtBin+1];

 
  Int_t nPoints = 5;

  char xHistoName[90];
  sprintf( xHistoName, "ptBin_%d_%d/pt2ndJet%sMean_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);

  Float_t x[nPoints];
  Float_t x_err[nPoints];

  getXPoints( file, xHistoName, nPoints, x, x_err);


  Float_t y_resp_recoPhot[nPoints];
  Float_t y_resp_recoPhot_err[nPoints];

  Float_t y_resp_genPhot[nPoints];
  Float_t y_resp_genPhot_err[nPoints];

  Float_t y_resp_recoGen[nPoints];
  Float_t y_resp_recoGen_err[nPoints];

  Float_t y_resp_genPart[nPoints];
  Float_t y_resp_genPart_err[nPoints];

  Float_t y_resp_partPhot[nPoints];
  Float_t y_resp_partPhot_err[nPoints];

  Float_t y_reso_recoPhot[nPoints];
  Float_t y_reso_recoPhot_err[nPoints];

  Float_t y_reso_genPhot[nPoints];
  Float_t y_reso_genPhot_err[nPoints];

  Float_t y_reso_recoGen[nPoints];
  Float_t y_reso_recoGen_err[nPoints];

  Float_t y_reso_genPart[nPoints];
  Float_t y_reso_genPart_err[nPoints];

  Float_t y_reso_partPhot[nPoints];
  Float_t y_reso_partPhot_err[nPoints];

  char yHistoName[100];
  sprintf( yHistoName, "ptBin_%d_%d/r_RecoPhot_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
  getYPoints( file, yHistoName, nPoints, y_resp_recoPhot, y_resp_recoPhot_err,  y_reso_recoPhot, y_reso_recoPhot_err);
  sprintf( yHistoName, "ptBin_%d_%d/r_GenPhot_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
  getYPoints( file, yHistoName, nPoints, y_resp_genPhot, y_resp_genPhot_err,  y_reso_genPhot, y_reso_genPhot_err);
  sprintf( yHistoName, "ptBin_%d_%d/r_RecoGen_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
  getYPoints( file, yHistoName, nPoints, y_resp_recoGen, y_resp_recoGen_err,  y_reso_recoGen, y_reso_recoGen_err);
  sprintf( yHistoName, "ptBin_%d_%d/r_GenPart_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
  getYPoints( file, yHistoName, nPoints, y_resp_genPart, y_resp_genPart_err,  y_reso_genPart, y_reso_genPart_err);
  sprintf( yHistoName, "ptBin_%d_%d/r_PartPhot_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
  getYPoints( file, yHistoName, nPoints, y_resp_partPhot, y_resp_partPhot_err,  y_reso_partPhot, y_reso_partPhot_err);

 
  //draw response histograms:

  TGraphErrors* gr_resp_recoPhot = new TGraphErrors(nPoints, x, y_resp_recoPhot, x_err, y_resp_recoPhot_err);
  gr_resp_recoPhot->SetMarkerStyle(20);
  gr_resp_recoPhot->SetMarkerColor(kRed);

  TGraphErrors* gr_resp_genPhot = new TGraphErrors(nPoints, x, y_resp_genPhot, x_err, y_resp_genPhot_err);
  gr_resp_genPhot->SetMarkerStyle(22);
  gr_resp_genPhot->SetMarkerColor(kBlack);

  TGraphErrors* gr_resp_recoGen = new TGraphErrors(nPoints, x, y_resp_recoGen, x_err, y_resp_recoGen_err);
  gr_resp_recoGen->SetMarkerStyle(21);
  gr_resp_recoGen->SetMarkerColor(kBlue);

  TGraphErrors* gr_resp_genPart = new TGraphErrors(nPoints, x, y_resp_genPart, x_err, y_resp_genPart_err);
  gr_resp_genPart->SetMarkerStyle(23);
  gr_resp_genPart->SetMarkerColor(kGreen);

  TGraphErrors* gr_resp_partPhot = new TGraphErrors(nPoints, x, y_resp_partPhot, x_err, y_resp_partPhot_err);
  gr_resp_partPhot->SetMarkerStyle(22);
  gr_resp_partPhot->SetMarkerColor(kBlack);


  Float_t lastX = x[nPoints-1];
  Float_t xMax_fit = lastX+2.;
  Float_t xMax_fit_green = x[2]+2.;

  Float_t xMax_axis;
  if( lastX <=12. ) 
    xMax_axis = 15.;
  else if( lastX <= 20. )
    xMax_axis = 25.;
  else if( lastX <= 25. )
    xMax_axis = 30.;
  else if( lastX <= 35. )
    xMax_axis = 40.;
  

  std::string xTitle;
  if( recoGen=="Reco" ) {
    xTitle = "p_{T}^{2ndJet reco} [GeV/c]";
  } else if( recoGen=="Gen" ) {
    xTitle = "p_{T}^{2ndJet gen} [GeV/c]";
  } else if( recoGen=="RecoRel" ) {
    xTitle = "p_{T}^{2ndJet reco}/p_{T}^{#gamma} [%]";
  } else if( recoGen=="GenRel" ) {
    xTitle = "p_{T}^{2ndJet gen}/p_{T}^{#gamma} [%]";
  }


  std::string fitFunct_name;
//  if( ptMin <=150. ) 
    fitFunct_name = "[0] - x*x*[1]";
//  else 
//    fitFunct_name = "[0] - x*[1]";

  TF1* fit_resp_genPart = new TF1("fit_resp_genPart", fitFunct_name.c_str());
  fit_resp_genPart->SetRange(0., xMax_fit);
  fit_resp_genPart->SetLineWidth(0.5);
  fit_resp_genPart->SetLineColor(kGreen);
  gr_resp_genPart->Fit( fit_resp_genPart, "RQ");

  TF1* fit_resp_genPhot = new TF1("fit_resp_genPhot", fitFunct_name.c_str());
  fit_resp_genPhot->SetRange(0., xMax_fit);
  fit_resp_genPhot->SetLineWidth(0.5);
  fit_resp_genPhot->SetLineColor(kBlack);
  gr_resp_genPhot->Fit( fit_resp_genPhot, "RQ");

  TF1* fit_resp_recoGen = new TF1("fit_resp_recoGen", "[0]");
  fit_resp_recoGen->SetRange(0., xMax_fit);
  fit_resp_recoGen->SetLineWidth(0.5);
  fit_resp_recoGen->SetLineColor(kBlue);
  gr_resp_recoGen->Fit( fit_resp_recoGen, "RQ");

  TF1* fit_resp_partPhot = new TF1("fit_resp_partPhot", fitFunct_name.c_str());
  fit_resp_partPhot->SetRange(0., xMax_fit);
  fit_resp_partPhot->SetLineWidth(0.5);
  fit_resp_partPhot->SetLineColor(kBlack);
  gr_resp_partPhot->Fit( fit_resp_partPhot, "RQ");

  std::string total_resp_str = (ONLYBLACK) ? "fit_resp_recoGen*fit_resp_genPhot" : "fit_resp_partPhot * fit_resp_genPart * fit_resp_recoGen";
  TF1* total_resp = new TF1("total_resp", total_resp_str.c_str());
  total_resp->SetRange(0., xMax_fit);
  total_resp->SetLineColor(kGray+2);

 
  fitFunct_name = "[0]*[1] - x*x*[2]";
  TF1* fit_respParabola = new TF1("fit_respParabola", fitFunct_name.c_str());
  fit_respParabola->SetRange(0., xMax_fit);
  if( NOQ ) //to evaluate syst
    fit_respParabola->FixParameter(1, 0.5*fit_resp_genPhot->GetParameter(0));
  else
    fit_respParabola->FixParameter(1, fit_resp_genPhot->GetParameter(0));
  fit_respParabola->SetLineColor(kRed);
  fit_respParabola->SetLineWidth(0.5);
  gr_resp_recoPhot->Fit( fit_respParabola, "RQ" );


  extrapResp = fit_respParabola->GetParameter(0);
  extrapResp_err = fit_respParabola->GetParError(0);

  intrResp = fit_resp_recoGen->GetParameter(0);
  intrResp_err = fit_resp_recoGen->GetParError(0);

  imbalanceResp = fit_resp_genPhot->GetParameter(0) - 1.;
  imbalanceResp_err = fit_resp_genPhot->GetParError(0);

  Float_t y0_partPhot = fit_resp_partPhot->GetParameter(0);
  Float_t y0_genPart = fit_resp_genPart->GetParameter(0);
  Float_t y0_recoGen = fit_resp_recoGen->GetParameter(0);
  Float_t y0_genPhot = fit_resp_genPhot->GetParameter(0);
  Float_t y0_err_partPhot = fit_resp_partPhot->GetParError(0);
  Float_t y0_err_genPart = fit_resp_genPart->GetParError(0);
  Float_t y0_err_recoGen = fit_resp_recoGen->GetParError(0);
  Float_t y0_err_genPhot = fit_resp_genPhot->GetParError(0);
  
  if( ONLYBLACK ) {
    trueResp = y0_genPhot*y0_recoGen;
    trueResp_err = y0_err_genPhot*y0_err_genPhot*y0_recoGen*y0_recoGen;
    trueResp_err += y0_genPhot*y0_genPhot*y0_err_recoGen*y0_err_recoGen;
    trueResp_err = sqrt( trueResp_err );
  } else {
    trueResp = y0_partPhot*y0_genPart*y0_recoGen;
    trueResp_err = y0_err_partPhot*y0_err_partPhot*y0_genPart*y0_genPart*y0_recoGen*y0_recoGen;
    trueResp_err += y0_partPhot*y0_partPhot*y0_err_genPart*y0_err_genPart*y0_recoGen*y0_recoGen;
    trueResp_err += y0_partPhot*y0_partPhot*y0_genPart*y0_genPart*y0_err_recoGen*y0_err_recoGen;
    trueResp_err = sqrt( trueResp_err );
  }

  //Float_t yMin_axis = (Lower[iPtBin]<50.) ? 0.7 : 0.8;
  Float_t yMin_axis = 0.8;
  TH2D* h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, yMin_axis, 1.1);
  h2_axes_resp->SetXTitle(xTitle.c_str());
  h2_axes_resp->GetXaxis()->SetTitleOffset(1.1);
  h2_axes_resp->SetYTitle("Response");
  h2_axes_resp->GetYaxis()->SetTitleOffset(1.2);

  TLegend* legend_resp = new TLegend(0.65, 0.65, 0.88, 0.88);
  legend_resp->SetTextSize(0.04);
  legend_resp->SetFillStyle(0);
  legend_resp->SetFillColor(kWhite);
  legend_resp->AddEntry(gr_resp_recoGen, "Intrinsic", "P");
  if( ONLYBLACK ) {
    legend_resp->AddEntry(gr_resp_genPhot, "Imbalance", "P");
  } else {
    legend_resp->AddEntry(gr_resp_partPhot, "Imbalance", "P");
    legend_resp->AddEntry(gr_resp_genPart, "Out of cone", "P");
  }
  legend_resp->AddEntry(total_resp, "Total", "L");
  legend_resp->AddEntry(gr_resp_recoPhot, "#gamma + jet", "P");

  char labeltext[50];
  sprintf(labeltext, "%d < p_{T}^{#gamma} < %d GeV/c", (int)ptMin, (int)ptMax);  
  TPaveText* label_resp = new TPaveText(0.15, 0.15, 0.4, 0.2, "brNDC");
  label_resp->SetFillColor(kWhite);
  label_resp->SetTextSize(0.035);
  label_resp->AddText(labeltext);

  TCanvas* c1_resp = new TCanvas("c1_resp", "c1_resp", 800, 600);
  c1_resp->cd();
  h2_axes_resp->Draw();
  total_resp->Draw("same");
  gr_resp_recoGen->Draw("Psame");
  if( ONLYBLACK ) {
    gr_resp_genPhot->Draw("Psame");
  } else { 
    gr_resp_genPart->Draw("Psame");
    gr_resp_partPhot->Draw("Psame");
  }
  gr_resp_recoPhot->Draw("Psame");
  legend_resp->Draw("same");
  label_resp->Draw("same");

  char canvasName_resp[100];
  sprintf(canvasName_resp, "%s/response_ptBin_%d_%d_%s.%s", OUTPUTDIR_.c_str(), (int)ptMin, (int)ptMax, recoGen.c_str(), PICTURE_FORMAT.c_str());
  c1_resp->SaveAs(canvasName_resp);

  delete c1_resp;
  delete h2_axes_resp;
  delete legend_resp;
  delete label_resp;



  //draw resolution histograms:

  TGraphErrors* gr_reso_recoPhot = new TGraphErrors(nPoints, x, y_reso_recoPhot, x_err, y_reso_recoPhot_err);
  gr_reso_recoPhot->SetMarkerStyle(20);
  gr_reso_recoPhot->SetMarkerColor(kRed);

  TGraphErrors* gr_reso_genPhot = new TGraphErrors(nPoints, x, y_reso_genPhot, x_err, y_reso_genPhot_err);
  gr_reso_genPhot->SetMarkerStyle(22);
  gr_reso_genPhot->SetMarkerColor(kBlack);

  TGraphErrors* gr_reso_recoGen = new TGraphErrors(nPoints, x, y_reso_recoGen, x_err, y_reso_recoGen_err);
  gr_reso_recoGen->SetMarkerStyle(21);
  gr_reso_recoGen->SetMarkerColor(kBlue);

  TGraphErrors* gr_reso_genPart = new TGraphErrors(nPoints, x, y_reso_genPart, x_err, y_reso_genPart_err);
  gr_reso_genPart->SetMarkerStyle(23);
  gr_reso_genPart->SetMarkerColor(kGreen);

  TGraphErrors* gr_reso_partPhot = new TGraphErrors(nPoints, x, y_reso_partPhot, x_err, y_reso_partPhot_err);
  gr_reso_partPhot->SetMarkerStyle(22);
  gr_reso_partPhot->SetMarkerColor(kBlack);


  Double_t x1, x2, y1, y2;
  TF1* fit_reso_genPart;
  if( FIT_RMS=="RMS" ) {
    fit_reso_genPart = new TF1("fit_reso_genPart", "[0] + x*[1]");
    fit_reso_genPart->SetRange(0., xMax_fit);
  } else {
    fit_reso_genPart = new TF1("fit_reso_genPart", "[0]");
    fit_reso_genPart->SetRange(0., xMax_fit_green);
  }
  gr_reso_genPart->GetPoint(1, x1, y1);
  fit_reso_genPart->SetParameter(0, y1);
  fit_reso_genPart->SetLineWidth(0.5);
  fit_reso_genPart->SetLineColor(kGreen);
  gr_reso_genPart->Fit(fit_reso_genPart, "RQ");
  fit_reso_genPart->SetRange(0., xMax_fit);

  TF1* fit_reso_partPhot = new TF1("fit_reso_partPhot", "[0] + x*[1]");
  fit_reso_partPhot->SetRange(0., xMax_fit);
  gr_reso_partPhot->GetPoint(1, x1, y1);
  gr_reso_partPhot->GetPoint(2, x2, y2);
  fit_reso_partPhot->SetParameter(0, 0.);
  fit_reso_partPhot->SetParameter(1, (y2-y1)/(x2-x1));
  fit_reso_partPhot->SetLineWidth(0.5);
  fit_reso_partPhot->SetLineColor(kBlack);
  gr_reso_partPhot->Fit(fit_reso_partPhot, "RQ");

  TF1* fit_reso_genPhot = new TF1("fit_reso_genPhot", "[0] + x*[1]");
  fit_reso_genPhot->SetRange(0., xMax_fit);
  fit_reso_genPhot->SetLineWidth(0.5);
  fit_reso_genPhot->SetLineColor(kBlack);
  gr_reso_genPhot->Fit(fit_reso_genPhot, "RQ");

//  std::cout << " ptBin: " << ptMin << "-" << ptMax << "\t" << "slope: " << fit_reso_genPhot->GetParameter(1) << "+/-" << fit_reso_genPhot->GetParError(1)<< std::endl;

  TF1* fit_reso_recoGen = new TF1("fit_reso_recoGen", "[0]");
  fit_reso_recoGen->SetRange(0., xMax_fit);
  gr_reso_recoGen->GetPoint(1, x1, y1);
  gr_reso_recoGen->GetPoint(2, x2, y2);
  fit_reso_recoGen->SetParameter(0, y1);
  fit_reso_recoGen->SetLineWidth(0.5);
  fit_reso_recoGen->SetLineColor(kBlue);
  gr_reso_recoGen->Fit(fit_reso_recoGen, "RQ");



  TF1* sum_intrinsic = new TF1("sum_intrinsic", "sqrt( fit_reso_genPart*fit_reso_genPart + fit_reso_recoGen*fit_reso_recoGen )");
  sum_intrinsic->SetRange(0., xMax_fit);
  sum_intrinsic->SetLineColor(kViolet);
  sum_intrinsic->SetLineStyle(3);


  std::string sum_str;
  if( ONLYBLACK ) {
    sum_str = "sqrt( fit_reso_recoGen*fit_reso_recoGen + fit_reso_genPhot*fit_reso_genPhot )";
  } else {
    sum_str = "sqrt( sum_intrinsic*sum_intrinsic + fit_reso_partPhot*fit_reso_partPhot )";
  }

  TF1* sum= new TF1("sum", sum_str.c_str());
  sum->SetRange(0., xMax_fit);
  sum->SetLineColor(kGray+2);

  Double_t x1_reco, y1_reco;
  gr_reso_recoPhot->GetPoint(1, x1_reco, y1_reco);
  Double_t x2_reco, y2_reco;
  gr_reso_recoPhot->GetPoint(2, x2_reco, y2_reco);

  // genPhot: y = mx + q
  // recoGen: y = c
  Float_t c = fit_reso_recoGen->GetParameter(0);
  Float_t q = fit_reso_genPhot->GetParameter(0);
  Float_t m = fit_reso_genPhot->GetParameter(1);
  // [0] = c; [1] = q; [2] = m
  TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "sqrt([0]*[0] + [1]*[1] + 2.*[1]*[2]*x + [2]*[2]*x*x)");
  fit_extrapToZero_sqrt->SetRange(0., xMax_fit);
  fit_extrapToZero_sqrt->SetParameter(0, c);
  fit_extrapToZero_sqrt->SetParLimits(0, 0.01, 0.3);
  if( NOQ )
    fit_extrapToZero_sqrt->FixParameter(1, 0.5*q); //to evaluate syst
  else
    fit_extrapToZero_sqrt->FixParameter(1, q); //fixed
  //fit_extrapToZero_sqrt->SetParLimits(1, 0., 0.1);
  fit_extrapToZero_sqrt->SetParameter(2, m);
  fit_extrapToZero_sqrt->SetParLimits(2, 0., 0.05);
  fit_extrapToZero_sqrt->SetLineColor(kRed);
  fit_extrapToZero_sqrt->SetLineWidth(0.5);

  TF1* fit_extrapToZero_line = new TF1("fit_extrapToZero_line", "[0]+[1]*x");
  fit_extrapToZero_line->SetRange(0., xMax_fit);
  fit_extrapToZero_line->SetParameter(0, y1_reco);
  fit_extrapToZero_line->SetParameter(1, (y2_reco-y1_reco)/(x2_reco-x1_reco));
  fit_extrapToZero_line->SetParLimits(0, 0., 0.4);
  fit_extrapToZero_line->SetLineColor(kRed);
  fit_extrapToZero_line->SetLineWidth(0.5);

  gr_reso_recoPhot->Fit(fit_extrapToZero_sqrt, "RQ");
  extrapReso = fit_extrapToZero_sqrt->GetParameter(0);
  extrapReso_err = fit_extrapToZero_sqrt->GetParError(0);
  //extrapReso = sqrt(fit_extrapToZero_sqrt->GetParameter(0));
  //extrapReso_err = ( 1. / (2.*sqrt(extrapReso)) * fit_extrapToZero_sqrt->GetParError(0)); //error propag

  intrReso = fit_reso_recoGen->GetParameter(0);
  intrReso_err = fit_reso_recoGen->GetParError(0);

  trueReso = sum->Eval(0.);
  if( ONLYBLACK )
    trueReso_err = sqrt( fit_reso_genPhot->GetParError(0)*fit_reso_genPhot->GetParError(0) + fit_reso_recoGen->GetParError(0)*fit_reso_recoGen->GetParError(0) );
  else 
    trueReso_err = sqrt( fit_reso_genPart->GetParError(0)*fit_reso_genPart->GetParError(0) + fit_reso_recoGen->GetParError(0)*fit_reso_recoGen->GetParError(0) + fit_reso_partPhot->GetParError(0)*fit_reso_partPhot->GetParError(0) );



  Float_t ymax;
  if( Lower[iPtBin] < 40. ) ymax = (RECOTYPE_=="pf") ? 0.5 : 0.8;
  else if( Lower[iPtBin] < 70. ) ymax = (RECOTYPE_=="pf") ? 0.4 : 0.6;
  else ymax = 0.3;

  TH2D* h2_axes_reso = new TH2D("axes_reso", "", 10, 0., xMax_axis, 10, 0., ymax);
  h2_axes_reso->SetXTitle(xTitle.c_str());
  h2_axes_reso->GetXaxis()->SetTitleOffset(1.1);
  h2_axes_reso->SetYTitle("Resolution");
  h2_axes_reso->GetYaxis()->SetTitleOffset(1.2);

  TPaveText* label_reso = new TPaveText(0.67, 0.83, 0.8, 0.86, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labeltext);

  Float_t minLegend = (ONLYBLACK) ? 0.2 : 0.15;
  TLegend* legend_reso = new TLegend(minLegend, 0.6, 0.45, 0.85);
  legend_reso->SetTextSize(0.04);
  legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->AddEntry(gr_reso_recoGen, "Intrinsic", "P");
  if( ONLYBLACK ) {
    legend_reso->AddEntry(gr_reso_genPhot, "Imbalance", "P");
  } else {
    legend_reso->AddEntry(gr_reso_partPhot, "Imbalance", "P");
    legend_reso->AddEntry(gr_reso_genPart, "Out of Cone", "P");
  }
  if( !ONLYBLACK )
    legend_reso->AddEntry(sum_intrinsic, "Intrinsic #oplus OoC", "L");
  legend_reso->AddEntry(sum, "Total", "L");
  legend_reso->AddEntry(gr_reso_recoPhot, "#gamma + jet", "P");

  TCanvas* c1_reso = new TCanvas("c1_reso", "c1_reso", 800, 600);
  c1_reso->cd();
  h2_axes_reso->Draw();
  sum->Draw("same");
  if( !ONLYBLACK )
    sum_intrinsic->Draw("same");
  gr_reso_recoPhot->Draw("Psame");
  gr_reso_recoGen->Draw("Psame");
  if( ONLYBLACK ) {
    gr_reso_genPhot->Draw("Psame");
  } else {
    gr_reso_genPart->Draw("Psame");
    gr_reso_partPhot->Draw("Psame");
  }
  legend_reso->Draw("same");
  label_reso->Draw("same");


  char canvasName_reso[100];
  sprintf(canvasName_reso, "%s/resolution_ptBin_%d_%d_%s.%s", OUTPUTDIR_.c_str(), (int)ptMin, (int)ptMax, recoGen.c_str(), PICTURE_FORMAT.c_str());
  c1_reso->SaveAs(canvasName_reso);

  if( ptMin==197. ) {
    h2_axes_reso->Draw();
    sum->Draw("same");
    if( !ONLYBLACK )
      sum_intrinsic->Draw("same");
    gr_reso_recoGen->Draw("Psame");
    if( ONLYBLACK ) {
      gr_reso_genPhot->Draw("Psame");
    } else {
      gr_reso_genPart->Draw("Psame");
      gr_reso_partPhot->Draw("Psame");
    }
    legend_reso->Draw("same");
    label_reso->Draw("same");
    sprintf(canvasName_reso, "%s/resolution_ptBin_%d_%d_%s_NORECO.%s", OUTPUTDIR_.c_str(), (int)ptMin, (int)ptMax, recoGen.c_str(), PICTURE_FORMAT.c_str());
    c1_reso->SaveAs(canvasName_reso);
  }

  delete c1_reso;
  delete h2_axes_reso;
  delete legend_reso;
  delete label_reso;


  delete gr_resp_recoPhot;
  delete gr_resp_genPart;
  delete gr_resp_partPhot;
  delete gr_resp_recoGen;

  delete gr_reso_recoPhot;
  delete gr_reso_genPart;
  delete gr_reso_partPhot;
  delete gr_reso_recoGen;

}



void getXPoints( TFile * file, const char* xHistoName, Int_t nPoints, Float_t* x, Float_t* x_err) {


  for( int i=0; i<nPoints; ++i ) {
    char fullName[100];
    sprintf( fullName, "%s_%d", xHistoName,  i);
    TH1D* h1_pt2ndJetMean = (TH1D*)file->Get( fullName );
    x[i] = h1_pt2ndJetMean->GetMean();
    x_err[i] =  h1_pt2ndJetMean->GetRMS()/sqrt((Float_t)h1_pt2ndJetMean->GetEntries());
  }


} //getxpoints




void getYPoints( TFile * file, const char* yHistoName, Int_t nPoints, Float_t* y_resp, Float_t* y_resp_err,  Float_t* y_reso, Float_t* y_reso_err) {

  for( int i=0; i<nPoints; ++i ) {

    char fullName[100];
    sprintf( fullName, "%s_%d", yHistoName, i);

    TH1D* h1_r = (TH1D*)file->Get(fullName);

    //this ugly fix saves from empty relative pt binning plots
    if( h1_r->GetEntries()==0 ) {
      y_resp[i]=-1.;
      y_resp_err[i]=0.;
      y_reso[i]=-1.;
      y_reso_err[i]=0.;
      continue;
    }


    Float_t mean = h1_r->GetMean(); 
    Float_t rms = h1_r->GetRMS();
    Float_t mean_err = h1_r->GetMeanError();
    //Float_t rms_err = h1_r->GetRMSError();
    Float_t rms_err = (LUMI>0) ? rms/sqrt((Float_t)LUMI*h1_r->Integral()*(Float_t)INT_PERC/100.) : h1_r->GetRMSError();

    //Float_t mean, mean_err, rms, rms_err;
    //fitTools::getProjectionMeanAndRMS(h1_r, mean, mean_err, rms, rms_err, 1., 0.9);

    TF1* gaussian = new TF1("gaussian", "gaus");
    gaussian->SetLineColor(kRed);
    fitTools::fitProjection_sameArea(h1_r, gaussian, (Float_t)(INT_PERC/100.), "RQ");
    //fitTools::fitProjection(h1_r, gaussian, 2., "RQ");

    Float_t mu = gaussian->GetParameter(1);
    Float_t sigma = gaussian->GetParameter(2);
    //TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
    //fitTools::fitProjection(h1_r, gaussian_chi, 1.5, "RQN");
    Float_t mu_err = gaussian->GetParError(1);
    //Float_t sigma_err = gaussian->GetParError(2);
    Float_t sigma_err = (LUMI>0) ? sigma/sqrt((Float_t)LUMI*h1_r->Integral()*(Float_t)INT_PERC/100.) : h1_r->GetRMSError();

    if( FIT_RMS == "FIT" ) {
      y_resp[i] = mu;
      y_resp_err[i] = mu_err;
    
      y_reso[i] = sigma/mu;
      y_reso_err[i] = sqrt( sigma_err*sigma_err/(mu*mu) + sigma*sigma*mu_err*mu_err/(mu*mu*mu*mu) ); 

    } else if( FIT_RMS == "RMS" ) {

      y_resp[i] = mean;
      y_resp_err[i] = mean_err;
      
      y_reso[i] = rms/mean;
      y_reso_err[i] = sqrt( rms_err*rms_err/(mean*mean) + rms*rms*mean_err*mean_err/(mean*mean*mean*mean) ); 

    } else if( FIT_RMS == "MIX" ) {

      y_resp[i] = mean;
      y_resp_err[i] = mean_err;
      
      y_reso[i] = sigma/mean;
      y_reso_err[i] = sqrt( sigma_err*sigma_err/(mean*mean) + sigma*sigma*mean_err*mean_err/(mean*mean*mean*mean) ); 

    }
  } //for

} //getYPoints



void drawGraphs_vs_pt( const std::string& varY, const std::string& canvasName, TGraphErrors* gr_resoReco_cutOn2ndJet, TGraphErrors* gr_intrReso_vs_pt, TGraphErrors* gr_extrapReso_vs_pt, TGraphErrors* gr_trueReso_vs_pt) {

  Float_t yMin, yMax;
  
  if( varY=="Response" ) {
    yMin = (RECOTYPE_=="pf") ? 0.75 : 0.2;
    yMax = 1.05;
  } else {
    yMin = 0.03;
    yMax = (RECOTYPE_=="pf") ? 0.4 : 0.7;
  }

  TH2D* h2_axes = new TH2D("axes", "", 10, 20., 2200., 10, yMin, yMax);
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);
  h2_axes->SetYTitle(varY.c_str());
  h2_axes->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();
  
  std::string legendTitle = "   |#eta^{jet}| < 1.4";

  Float_t leg_xMin = (varY=="Response") ? 0.5 : 0.5;
  Float_t leg_yMin = (varY=="Response") ? 0.15: 0.55;
  Float_t leg_xMax = (varY=="Response") ? 0.88: 0.88;
  Float_t leg_yMax = (varY=="Response") ? 0.45 : 0.88;

  TLegend* legend = new TLegend(leg_xMin, leg_yMin, leg_xMax, leg_yMax, legendTitle.c_str());
  legend->SetFillColor(kWhite);
  //legend->SetFillStyle(1);
  //legend->SetBorderSize(0.5);
  legend->AddEntry(gr_intrReso_vs_pt, "Intrinsic", "P");
  legend->AddEntry(gr_resoReco_cutOn2ndJet, "#gamma+Jet, cut on 2^{nd} Jet", "P");
  if( gr_trueReso_vs_pt != 0 )
    legend->AddEntry(gr_trueReso_vs_pt, "MC Truth", "P");
  if( gr_extrapReso_vs_pt != 0 )
    legend->AddEntry(gr_extrapReso_vs_pt, "#gamma+Jet, p_{T}^{2nd Jet} #rightarrow 0", "P");


  Float_t lab_yMin = (varY=="Response") ? 0.8 : 0.15;
  Float_t lab_yMax = (varY=="Response") ? 0.88 : 0.2;
  TPaveText* label = new TPaveText(0.15, lab_yMin, 0.3, lab_yMax, "brNDC");
  label->SetFillColor(kWhite);
  label->SetTextSize(0.04);
  label->AddText("Anti-kt 0.5");
  
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  c1->SetLogx();
  if( varY=="Resolution" )
    c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  h2_axes->Draw();
  gr_resoReco_cutOn2ndJet->Draw("Psame");
  if( gr_trueReso_vs_pt != 0 )
    gr_trueReso_vs_pt->Draw("Psame");
  gr_intrReso_vs_pt->Draw("P same");
  if( gr_extrapReso_vs_pt != 0 )
    gr_extrapReso_vs_pt->Draw("P same");
  legend->Draw("same");
  label->Draw("same");
  c1->SaveAs(canvasName.c_str()); 
 
  delete c1;
  delete h2_axes;

} //drawGraphs_vs_pt



void drawPullHistogram( const std::string& name, Float_t* pull, Float_t* pull_err) {

  Double_t Lower[NBINS_PT-1];
  fitTools::getPtBins(NBINS_PT, Lower);

  TH1D* h1_pull = new TH1D(name.c_str(), "", NBINS_PT-1, Lower);

  for(unsigned i = 0; i<(NBINS_PT-1); ++i ) {
    h1_pull->SetBinContent(i, pull[i]);
    h1_pull->SetBinError(i, pull_err[i]);
  }


  h1_pull->SetFillColor(43);
  if( name == "pullReso" ) {
    h1_pull->GetYaxis()->SetRangeUser(-8., 8.);
    h1_pull->GetYaxis()->SetNdivisions(10, kTRUE);
  } else {
    h1_pull->GetYaxis()->SetRangeUser(-0.25, 0.25);
    h1_pull->GetYaxis()->SetNdivisions(7, kTRUE);
  }
  h1_pull->GetXaxis()->SetTitleOffset(1.1);
  h1_pull->GetYaxis()->SetTitleOffset(1.2);
  h1_pull->SetXTitle("p_{T}^{#gamma} [GeV/c]");
  h1_pull->SetYTitle("Pull [%]");
  h1_pull->GetXaxis()->SetMoreLogLabels();
  h1_pull->GetXaxis()->SetNoExponent();

  std::string canvasName;
  TCanvas* c0 = new TCanvas("c0", "c0", 800, 600);
  c0->cd();
  c0->SetLogx();
  c0->SetGridy();
  h1_pull->Draw("HIST");
  canvasName = OUTPUTDIR_ + "/" + name + "." + PICTURE_FORMAT;
  c0->SaveAs(canvasName.c_str());
  delete c0;
  delete h1_pull;

}
