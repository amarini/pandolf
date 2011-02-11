#include "TCanvas.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "fitTools.h"

#include "TMinuit.h"
#include "TRandom3.h"
#include "TMatrixD.h"


struct TwoHistos {
  TH1D* h1_resp;
  TH1D* h1_reso;
};

struct TwoFunctions {
  TF1* f1_resp;
  TF1* f1_reso;
};


bool MCASSOC_ = false;
bool NOHEMISPHERE_ = false;
std::string FITRMS_ = "FIT";
std::string JETALGO_;
std::string RECOTYPE_;



TwoHistos getSyst_q();
TwoHistos getSyst_flavor();
TwoHistos getSyst_QCD();
TwoFunctions getSyst_photonScale();
TwoHistos getStat();
TGraphErrors* getGraph(std::string fileName, std::string graphName);
//void drawTotalError( std::string resp_reso, TH1D* syst_q, TH1D* syst_flavor, TH1D* syst_QCD, TF1* syst_photonScale, TH1D* stat );
void drawTotalError( std::string resp_reso, TH1D* syst_q, TH1D* syst_flavor, TF1* syst_photonScale, TH1D* stat );
void getBinning( Double_t* xvec, Double_t xmin=20.*0.9, Double_t xmax=1400.*1.1, Int_t npx=100, bool islog=true);
//TH1D *getBand(TF1 *f, TMatrixD const& m, std::string name, bool getRelativeBand=false, int npx=100);




int main( int argc, char* argv[]) {

  if( argc != 2 ) {
    std::cout << "USAGE: ./drawSystPhotonJetExtrap [recoType]" << std::endl;
    exit(877);
  }

  std::string recoType_str(argv[1]);
  if( recoType_str=="pf" ) {
    JETALGO_ = recoType_str + "akt5";
  } else if( recoType_str=="calo" ) {
    JETALGO_ = "akt5";
  } else if( recoType_str=="jpt" ) {
    JETALGO_ = recoType_str+"ak5";
  } else {
    std::cout << "Unknown recoType '" << recoType_str << "'. Exiting." << std::endl;
    exit(878);
  }
  RECOTYPE_ = recoType_str;


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


  //syst:
  TwoHistos histos_q = getSyst_q();
  TwoHistos histos_flavor = getSyst_flavor();
  //TwoHistos histos_QCD = getSyst_QCD();

  TwoFunctions tf = getSyst_photonScale();
  TF1* syst_photonScale_resp = tf.f1_resp; 
  TF1* syst_photonScale_reso = tf.f1_reso; 

  //stat:
  TwoHistos histos_stat = getStat();

  //total:
//drawTotalError( "resp", histos_q.h1_resp, histos_flavor.h1_resp, histos_QCD.h1_resp, syst_photonScale, histos_stat.h1_resp );
//drawTotalError( "reso", histos_q.h1_reso, histos_flavor.h1_reso, histos_QCD.h1_reso, syst_photonScale, histos_stat.h1_reso );
  drawTotalError( "resp", histos_q.h1_resp, histos_flavor.h1_resp, syst_photonScale_resp, histos_stat.h1_resp );
  drawTotalError( "reso", histos_q.h1_reso, histos_flavor.h1_reso, syst_photonScale_reso, histos_stat.h1_reso );



  return 0;

}

TGraphErrors* getGraph(std::string fileName, std::string graphName) {

  TFile* file = TFile::Open(fileName.c_str());
  TGraphErrors* graph = (TGraphErrors*)file->Get(graphName.c_str());
  file->Close();

  return graph;

}



TwoHistos getSyst_q() {

  Double_t xvec[100];
  getBinning( xvec );

  TH1D* syst_q_reso = new TH1D("syst_q_reso", "", 100, xvec);
  TH1D* syst_q_resp = new TH1D("syst_q_resp", "", 100, xvec);


  bool q15 = false;
  //std::string funcType_resp = (RECOTYPE_=="pf") ? "rpf" : "powerlawL2L3";
  std::string funcType_resp = "powerlawL2L3";
  std::string funcType_reso = (RECOTYPE_=="pf") ? "NSCPF" : "NSC";
  //std::string funcType = "NSC";

  std::string graphName = "gr_extrapReso_vs_pt";

  std::string fileName_reso = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3Raw_RMS99.root";
  if( RECOTYPE_=="pf" ) fileName_reso = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3_RMS99.root";
  TGraphErrors* gr_extrapReso_q = getGraph(fileName_reso, graphName);
  /*if( RECOTYPE_!="pf" )*/ gr_extrapReso_q->RemovePoint(0);
  if( RECOTYPE_=="calo" ) gr_extrapReso_q->RemovePoint(0);
  TF1* fit_extrapReso_q= fitTools::fitResolutionGraph(gr_extrapReso_q, funcType_reso, "fit_extrapReso_q", "RN");
  fit_extrapReso_q->SetLineColor(kRed);

  std::string fileName_resp = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3Raw_RMS99.root";
  TGraphErrors* gr_extrapResp_q = getGraph(fileName_resp, "gr_extrapResp_vs_pt");
  TF1* fit_extrapResp_q= fitTools::fitResponseGraph(gr_extrapResp_q, funcType_resp, "fit_extrapResp_q", "RN");
  fit_extrapResp_q->SetLineColor(kRed);


  fileName_resp = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3Raw_RMS99_NOQ.root";
  fileName_reso = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3Raw_RMS99_NOQ.root";
  if( RECOTYPE_=="pf" ) fileName_reso = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3_RMS99_NOQ.root";


  TGraphErrors* gr_extrapReso_noq = getGraph(fileName_reso, graphName);
  /*if( RECOTYPE_!="pf" )*/ gr_extrapReso_noq->RemovePoint(0);
  if( RECOTYPE_=="calo" ) gr_extrapReso_noq->RemovePoint(0);
  TF1* fit_extrapReso_noq = fitTools::fitResolutionGraph(gr_extrapReso_noq, funcType_reso, "fit_extrapReso_noq", "RN");
  fit_extrapReso_noq->SetLineColor(kBlue);

  TGraphErrors* gr_extrapResp_noq = getGraph(fileName_resp, "gr_extrapResp_vs_pt");
  TF1* fit_extrapResp_noq = fitTools::fitResponseGraph(gr_extrapResp_noq, funcType_resp, "fit_extrapResp_noq", "RN");
  fit_extrapResp_noq->SetLineColor(kBlue);

  TH1D* band_resp_q = fitTools::getBand(fit_extrapResp_q, "band_resp_q");
  TH1D* band_resp_noq = fitTools::getBand(fit_extrapResp_noq, "band_resp_noq");

  TH1D* band_reso_q = fitTools::getBand(fit_extrapReso_q, "band_reso_q");
  TH1D* band_reso_noq = fitTools::getBand(fit_extrapReso_noq, "band_reso_noq");


  for( unsigned iBin=0; iBin<99; ++iBin ) {

    Double_t binCenter = 0.5*(xvec[iBin+1]+xvec[iBin]);

    Double_t syst_resp_value = 100.*fabs( fit_extrapResp_q->Eval(binCenter) - fit_extrapResp_noq->Eval(binCenter) )/fit_extrapResp_q->Eval(binCenter);
    syst_q_resp->SetBinContent( iBin+1, syst_resp_value );

    Double_t syst_reso_value = 100.*fabs( fit_extrapReso_q->Eval(binCenter) - fit_extrapReso_noq->Eval(binCenter) )/fit_extrapReso_q->Eval(binCenter);
    Double_t err_reso = 100.*sqrt( band_reso_q->GetBinError(iBin+1)*band_reso_q->GetBinError(iBin+1) + band_reso_noq->GetBinError(iBin+1)*band_reso_noq->GetBinError(iBin+1) );
      syst_q_reso->SetBinContent( iBin+1, syst_reso_value );
  //if( syst_reso_value > err_reso ) {
  //  syst_q_reso->SetBinContent( iBin+1, syst_reso_value );
  //} else {
  //  syst_q_reso->SetBinContent( iBin+1, err_reso );
  //}

  }

  TwoHistos th;
  th.h1_resp = syst_q_resp;
  th.h1_reso = syst_q_reso;

  TH2D* axesreso = new TH2D("axesreso", "", 10, 20., 1400., 10, 0., 0.4);
  TCanvas* cSyst_q = new TCanvas("cSyst_q", "errorGraph", 600, 600);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  axesreso->Draw();
  //band_reso_noq->Draw("c e3 same");
  //band_reso_q->Draw("c e3 same");
  fit_extrapReso_noq->Draw("same");
  fit_extrapReso_q->Draw("same");
  gr_extrapReso_q->SetMarkerStyle(20);
  gr_extrapReso_noq->SetMarkerStyle(20);
  gr_extrapReso_q->SetMarkerColor(kRed);
  gr_extrapReso_noq->SetMarkerColor(kBlue);
  gr_extrapReso_q->Draw("p same");
  gr_extrapReso_noq->Draw("p same");
  std::string canvName = "syst_reso_q_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  TH2D* axesresp = new TH2D("axesresp", "", 10, 20., 1400., 10, 0., 1.3);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  axesresp->Draw();
  band_resp_noq->Draw("c e3 same");
  band_resp_q->Draw("c e3 same");
  fit_extrapResp_noq->Draw("same");
  fit_extrapResp_q->Draw("same");
  gr_extrapResp_q->SetMarkerStyle(20);
  gr_extrapResp_noq->SetMarkerStyle(20);
  gr_extrapResp_q->SetMarkerColor(kRed);
  gr_extrapResp_noq->SetMarkerColor(kBlue);
  gr_extrapResp_q->Draw("p same");
  gr_extrapResp_noq->Draw("p same");
  canvName = "syst_resp_q_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  delete axesresp;
  delete axesreso;
  delete cSyst_q;

  return th;

}


//void getSyst_flavor(TH1D* syst_flavor_reso, TH1D* syst_flavor_resp) {
TwoHistos getSyst_flavor() {

  Double_t xvec[100];
  getBinning( xvec );

  TH1D* syst_flavor_reso = new TH1D("syst_flavor_reso", "", 100, xvec);
  TH1D* syst_flavor_resp = new TH1D("syst_flavor_resp", "", 100, xvec);


  std::string fileName_resp, fileName_reso;

  fileName_resp = "PhotonJetGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_GLUON_2ndJet10.root";
  fileName_reso = "PhotonJetGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_GLUON_2ndJet10.root";
  TGraphErrors* gr_intrReso_gluon = getGraph(fileName_reso, "GENresolution_L2L3_eta011_vs_pt");
  //TGraphErrors* gr_intrReso_gluon = getGraph(fileName_reso, "GENresolution_eta011_vs_pt");
  TF1* fit_intrReso_gluon_NSC = fitTools::fitResolutionGraph(gr_intrReso_gluon, "NSC", "fit_intrReso_gluon_NSC");
  TF1* fit_intrReso_gluon = fitTools::fitResolutionGraph(gr_intrReso_gluon, "NSCPF", "fit_intrReso_gluon");
  TGraphErrors* gr_intrResp_gluon = getGraph(fileName_resp, "GENresponse_L2L3_eta011_vs_pt");
  //TGraphErrors* gr_intrResp_gluon = getGraph(fileName_resp, "GENresponse_eta011_vs_pt");
  //TF1* fit_intrResp_gluon = fitTools::fitResponseGraph(gr_intrResp_gluon, "rpf", "fit_intrResp_gluon");
  TF1* fit_intrResp_gluon= fitTools::fitResponseGraph(gr_intrResp_gluon, "powerlaw", "fit_intrResp_gluon", "RN", 800.);

  fileName_resp = "PhotonJetGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_QUARK_2ndJet10.root";
  fileName_reso = "PhotonJetGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_QUARK_2ndJet10.root";
  TGraphErrors* gr_intrReso_quark = getGraph(fileName_reso, "GENresolution_L2L3_eta011_vs_pt");
  //TGraphErrors* gr_intrReso_quark = getGraph(fileName_reso, "GENresolution_eta011_vs_pt");
  TF1* fit_intrReso_quark_NSC = fitTools::fitResolutionGraph(gr_intrReso_quark, "NSC", "fit_intrReso_quark_NSC");
  TF1* fit_intrReso_quark = fitTools::fitResolutionGraph(gr_intrReso_quark, "NSCPF", "fit_intrReso_quark");
  TGraphErrors* gr_intrResp_quark = getGraph(fileName_resp, "GENresponse_L2L3_eta011_vs_pt");
  //TGraphErrors* gr_intrResp_quark = getGraph(fileName_resp, "GENresponse_eta011_vs_pt");
  //TF1* fit_intrResp_quark = fitTools::fitResponseGraph(gr_intrResp_quark, "rpf", "fit_intrResp_quark");
  TF1* fit_intrResp_quark= fitTools::fitResponseGraph(gr_intrResp_quark, "powerlaw", "fit_intrResp_quark", "RN", 800.);

  TH1D* band_reso_quark = fitTools::getBand(fit_intrReso_quark_NSC, "band_intrReso_quark");
  TH1D* band_resp_quark = fitTools::getBand(fit_intrResp_quark, "band_intrResp_quark");
  TH1D* band_reso_gluon = fitTools::getBand(fit_intrReso_gluon_NSC, "band_intrReso_gluon");
  TH1D* band_resp_gluon = fitTools::getBand(fit_intrResp_gluon, "band_intrResp_gluon");

//std::cout << std::endl << std::endl;
//std::cout << "AFTER: " << band_reso_quark->GetName() << std::endl;
//for( unsigned ibin=1; ibin<band_reso_quark->GetNbinsX()+1; ++ibin)
//std::cout << band_reso_quark->GetBinContent(ibin) << "   " << band_reso_quark->GetBinError(ibin) << std::endl;

  band_reso_quark->SetFillColor(kYellow-9);
  band_resp_quark->SetFillColor(kYellow-9);
  band_reso_gluon->SetFillColor(kYellow-9);
  band_resp_gluon->SetFillColor(kYellow-9);


  for( unsigned iBin=0; iBin<99; ++iBin ) {

    Double_t binCenter = 0.5*(xvec[iBin+1]+xvec[iBin]);

    Double_t denom = 0.5*(fit_intrReso_quark->Eval(binCenter)+fit_intrReso_gluon->Eval(binCenter));
    Double_t syst_reso_value = 100.*fabs( fit_intrReso_quark->Eval(binCenter) - fit_intrReso_gluon->Eval(binCenter) )/denom;
    Double_t err_reso = 100.*sqrt( band_reso_quark->GetBinError(iBin+1)*band_reso_quark->GetBinError(iBin+1) + band_reso_gluon->GetBinError(iBin+1)*band_reso_gluon->GetBinError(iBin+1) );
    if( syst_reso_value > err_reso ) {
      syst_flavor_reso->SetBinContent( iBin+1, syst_reso_value );
    } else {
      syst_flavor_reso->SetBinContent( iBin+1, err_reso );
    }

    Double_t syst_resp_value = 100.*fabs( fit_intrResp_quark->Eval(binCenter) - fit_intrResp_gluon->Eval(binCenter))/(0.5*(fit_intrResp_quark->Eval(binCenter)+fit_intrResp_gluon->Eval(binCenter)));
    syst_flavor_resp->SetBinContent( iBin+1, syst_resp_value );

  }

  TH2D* axesreso = new TH2D("axesreso", "", 10, 20., 1400., 10, 0., 0.3);
  TCanvas* cSyst_q = new TCanvas("cSyst_q", "errorGraph", 600, 600);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  cSyst_q->SetLogy();
  axesreso->Draw();
  band_reso_quark->Draw("c e3 same");
  band_reso_gluon->Draw("c e3 same");
  gr_intrReso_quark->SetMarkerStyle(20);
  gr_intrReso_gluon->SetMarkerStyle(20);
  gr_intrReso_quark->SetMarkerColor(kRed);
  gr_intrReso_gluon->SetMarkerColor(kBlue);
  gr_intrReso_quark->Draw("p same");
  gr_intrReso_gluon->Draw("p same");
  fit_intrReso_quark->Draw("same");
  fit_intrReso_gluon->Draw("same");
  std::string canvName = "syst_reso_flavor_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  TH2D* axesresp = new TH2D("axesresp", "", 10, 20., 1400., 10, 0.7, 1.3);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogy(0);
  axesresp->Draw();
  gr_intrResp_quark->SetMarkerStyle(20);
  gr_intrResp_gluon->SetMarkerStyle(20);
  gr_intrResp_quark->SetMarkerColor(kRed);
  gr_intrResp_gluon->SetMarkerColor(kBlue);
  gr_intrResp_quark->Draw("p same");
  gr_intrResp_gluon->Draw("p same");
  fit_intrResp_quark->Draw("same");
  fit_intrResp_gluon->Draw("same");
  canvName = "syst_resp_flavor_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  TwoHistos th;
  th.h1_resp = syst_flavor_resp;
  th.h1_reso = syst_flavor_reso;
 
  return th;

}



TwoHistos getSyst_QCD() {

  Double_t xvec[100];
  getBinning( xvec );

  TH1D* syst_QCD_reso = new TH1D("syst_QCD_reso", "", 100, xvec);
  TH1D* syst_QCD_resp = new TH1D("syst_QCD_resp", "", 100, xvec);


  std::string fileName_resp, fileName_reso;

  fileName_resp = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011_RMS99.root";
  fileName_reso = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011_RMS99.root";
  TGraphErrors* gr_extrapResp_s = getGraph(fileName_resp, "gr_extrapResp_vs_pt");
  TF1* fit_extrapResp_s= fitTools::fitResponseGraph(gr_extrapResp_s, "rpf", "fit_extrapResp_s");
  TGraphErrors* gr_extrapReso_s = getGraph(fileName_reso, "gr_extrapReso_vs_pt");
  TF1* fit_extrapReso_s= fitTools::fitResolutionGraph(gr_extrapReso_s, "NSCPF", "fit_extrapReso_s");

  fileName_resp = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_QCD_TuneZ2_7TeV_pythia6_"+JETALGO_+"_LUMI_2ndJet10_eta011_RMS99.root";
  fileName_reso = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_QCD_TuneZ2_7TeV_pythia6_"+JETALGO_+"_LUMI_2ndJet10_eta011_RMS99.root";
  TGraphErrors* gr_extrapResp_sb = getGraph(fileName_resp, "gr_extrapResp_vs_pt");
  TF1* fit_extrapResp_sb= fitTools::fitResponseGraph(gr_extrapResp_sb, "rpf", "fit_extrapResp_sb");
  TGraphErrors* gr_extrapReso_sb = getGraph(fileName_reso, "gr_extrapReso_vs_pt");
  TF1* fit_extrapReso_sb= fitTools::fitResolutionGraph(gr_extrapReso_sb, "NSCPF", "fit_extrapReso_sb");

  for( unsigned iBin=0; iBin<99; ++iBin ) {

    Double_t binCenter = 0.5*(xvec[iBin+1]+xvec[iBin]);

    Double_t syst_reso_value = 100.*fabs( fit_extrapReso_s->Eval(binCenter) - fit_extrapReso_sb->Eval(binCenter) )/(0.5*(fit_extrapReso_s->Eval(binCenter)+fit_extrapReso_sb->Eval(binCenter)));
    syst_QCD_reso->SetBinContent( iBin+1, syst_reso_value );

    Double_t syst_resp_value = 100.*fabs( fit_extrapResp_s->Eval(binCenter) - fit_extrapResp_sb->Eval(binCenter))/(0.5*(fit_extrapResp_s->Eval(binCenter)+fit_extrapResp_sb->Eval(binCenter)));
    syst_QCD_resp->SetBinContent( iBin+1, syst_resp_value );

  }


  TH2D* axesreso = new TH2D("axesreso", "", 10, 20., 1400., 10, 0., 0.4);
  TCanvas* cSyst_q = new TCanvas("cSyst_q", "errorGraph", 600, 600);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  axesreso->Draw();
  gr_extrapReso_s->SetMarkerStyle(20);
  gr_extrapReso_sb->SetMarkerStyle(20);
  gr_extrapReso_s->SetMarkerColor(kRed);
  gr_extrapReso_sb->SetMarkerColor(kBlue);
  gr_extrapReso_s->Draw("p same");
  gr_extrapReso_sb->Draw("p same");
  std::string canvName = "syst_resp_QCD_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  TH2D* axesresp = new TH2D("axesresp", "", 10, 20., 1400., 10, 0.7, 1.3);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  axesresp->Draw();
  gr_extrapResp_s->SetMarkerStyle(20);
  gr_extrapResp_sb->SetMarkerStyle(20);
  gr_extrapResp_s->SetMarkerColor(kRed);
  gr_extrapResp_sb->SetMarkerColor(kBlue);
  gr_extrapResp_s->Draw("p same");
  gr_extrapResp_sb->Draw("p same");
  canvName = "syst_reso_QCD_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  TwoHistos th;
  th.h1_resp = syst_QCD_resp;
  th.h1_reso = syst_QCD_reso;

  return th;

}



TwoFunctions getSyst_photonScale() {

  TF1* syst_photonScale_resp = new TF1("syst_photonScale_resp", "[0]", 20., 1400.);
  syst_photonScale_resp->SetParameter(0, 1.);

  TF1* syst_photonScale_reso = new TF1("syst_photonScale_reso", "[0]", 20., 1400.);
  syst_photonScale_reso->SetParameter(0, 2.);

  TwoFunctions tf;
  tf.f1_resp = syst_photonScale_resp; 
  tf.f1_reso = syst_photonScale_reso; 

  return tf;

}



TwoHistos getStat() {


  std::string fileName = "PhotonJetExtrapGraphs_DATA_Nov4ReReco_vs_G_TuneZ2_7TeV_pythia6_CORR_"+JETALGO_+"_LUMI_2ndJet10_eta011L2L3Raw_RMS99.root";
  TFile* fileDATA = TFile::Open(fileName.c_str());

  TGraphErrors* gr_respDATA = (TGraphErrors*)fileDATA->Get("gr_DATAResp_vs_pt");
  TGraphErrors* gr_resoDATA = (TGraphErrors*)fileDATA->Get("gr_DATAReso_vs_pt");

  TF1* fit_respDATA = fitTools::fitResponseGraph( gr_respDATA, "powerlaw", "fit_respDATA");
  TF1* fit_resoDATA = fitTools::fitResolutionGraph( gr_resoDATA, "NSC", "fit_resoDATA");

  TH1D* band_respDATA = fitTools::getBand(fit_respDATA, "band_respDATA");
  TH1D* band_resoDATA = fitTools::getBand(fit_resoDATA, "band_resoDATA");

  Double_t xvec[100];
  getBinning( xvec );

  TH1D* stat_reso = new TH1D("stat_reso", "", 100, xvec);
  TH1D* stat_resp = new TH1D("stat_resp", "", 100, xvec);


  for( unsigned iBin=0; iBin<99; ++iBin ) {

    Double_t binCenter = 0.5*(xvec[iBin+1]+xvec[iBin]);

    float resp_fit = fit_respDATA->Eval(binCenter);
    float resp_fit_err = band_respDATA->GetBinError(iBin);
    float error_resp = 100.*resp_fit_err/resp_fit;

    float reso_fit = fit_resoDATA->Eval(binCenter);
    float reso_fit_err = band_resoDATA->GetBinError(iBin);
    float error_reso = 100.*reso_fit_err/reso_fit;

    stat_resp->SetBinContent( iBin+1, error_resp );
    stat_reso->SetBinContent( iBin+1, error_reso );

  }


  TCanvas* cSyst_q = new TCanvas("cSyst_q", "errorGraph", 600, 600);
  TH2D* axesresp = new TH2D("axesresp", "", 10, 20., 1400., 10, 0.7, 1.3);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  axesresp->Draw();
  band_respDATA->SetFillColor(kYellow-9);
  band_respDATA->Draw("c e3 same");
  fit_respDATA->SetLineWidth(2);
  fit_respDATA->Draw("l same");
  std::string canvName = "stat_resp_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());

  TH2D* axesreso = new TH2D("axesreso", "", 10, 20., 1400., 10, 0., 0.4);
  cSyst_q->cd();
  cSyst_q->Clear();
  cSyst_q->SetLogx();
  axesreso->Draw();
  band_resoDATA->SetFillColor(kYellow-9);
  band_resoDATA->Draw("c e3 same");
  fit_resoDATA->SetLineWidth(2);
  fit_resoDATA->Draw("l same");
  canvName = "stat_reso_"+JETALGO_+".eps";
  cSyst_q->SaveAs(canvName.c_str());


  TwoHistos th;
  th.h1_resp = stat_resp;
  th.h1_reso = stat_reso;

  return th;

}

void getBinning( Double_t* xvec, Double_t xmin, Double_t xmax, Int_t npx, bool islog) {

  xvec[0] = xmin;
  double dx = (islog ? pow(xmax/xmin, 1./npx) : (xmax-xmin)/npx);
  for (int i = 0; i != npx; ++i) {
    xvec[i+1] = (islog ? xvec[i]*dx : xvec[i]+dx);
  }


}




//void drawTotalError( std::string resp_reso, TH1D* syst_q, TH1D* syst_flavor, TH1D* syst_QCD, TF1* syst_photonScale, TH1D* stat ) {
void drawTotalError( std::string resp_reso, TH1D* syst_q, TH1D* syst_flavor, TF1* syst_photonScale, TH1D* stat ) {

  syst_q->SetLineColor(kBlack);
  syst_q->SetLineStyle(2);
  syst_q->SetLineWidth(3);

  syst_flavor->SetLineColor(kBlack);
  syst_flavor->SetLineStyle(3);
  syst_flavor->SetLineWidth(3);

//syst_QCD->SetLineColor(8);
//syst_QCD->SetLineStyle(4);
//syst_QCD->SetLineWidth(3);

  syst_photonScale->SetLineColor(kBlue);
  syst_photonScale->SetLineStyle(1);
  syst_photonScale->SetLineWidth(3);


  Double_t xvec[100];
  getBinning( xvec, 20.*0.9, 1400*1.1, 100, (bool)true);

  TH1D* syst_total = new TH1D("syst_total", "", 100, xvec);
  for( unsigned iBin=0; iBin<100; ++iBin) {
    Double_t binCenter = 0.5*(xvec[iBin+1]+xvec[iBin]);
    Double_t squareBinValue = syst_q->GetBinContent(iBin+1)*syst_q->GetBinContent(iBin+1) + 
                        syst_flavor->GetBinContent(iBin+1)*syst_flavor->GetBinContent(iBin+1) + 
//                        syst_QCD->GetBinContent(iBin+1)*syst_QCD->GetBinContent(iBin+1) + 
                        syst_photonScale->Eval(binCenter)*syst_photonScale->Eval(binCenter);
    syst_total->SetBinContent(iBin+1, sqrt(squareBinValue));
  }
  syst_total->SetLineColor(kGray);
  syst_total->SetLineWidth(0.);
  syst_total->SetFillColor(kGray);
  syst_total->SetFillStyle(3001);

  stat->SetLineColor(kRed);
  stat->SetLineWidth(2);
  stat->SetFillColor(kRed);
  stat->SetFillStyle(3004);

  TPaveText* label_cms = new TPaveText( 0.22, 0.83, 0.42, 0.87, "brNDC" );
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.038);
  label_cms->SetTextFont(62);
  //label_cms->AddText("CMS Preliminary 2010");
  label_cms->AddText("CMS Simulation 2010");

  TPaveText* label_sqrt = new TPaveText( 0.22, 0.78, 0.42, 0.82, "brNDC" );
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  //label_sqrt->AddText("#sqrt{s} = 7 TeV, L = 34 pb^{-1}");
  label_sqrt->AddText("#sqrt{s} = 7 TeV");

  std::string legendTitleText;
  if( JETALGO_=="pfakt5" )
    legendTitleText="Anti-k_{T} 0.5 PFJets";
  else if( JETALGO_=="akt5" )
    legendTitleText="Anti-k_{T} 0.5 CaloJets";
  else if( JETALGO_=="jptak5" )
    legendTitleText="Anti-k_{T} 0.5 JPTJets";
  else 
    legendTitleText="[Unknown jet algo]";
   

  TLegend* legend = new TLegend(0.57, 0.6, 0.85, 0.88, legendTitleText.c_str());
  legend->SetFillColor(kWhite);
  //legend->SetFillStyle(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( syst_q, "MC Model", "L");
  legend->AddEntry( syst_flavor, "Jet Flavour", "L");
 // legend->AddEntry( syst_QCD, "QCD BG", "L");
  if( resp_reso=="resp" )
    legend->AddEntry( syst_photonScale, "Photon Scale", "L");
  else
    legend->AddEntry( syst_photonScale, "Photon", "L");
  legend->AddEntry( syst_total, "Total Syst.", "F");
  //legend->AddEntry( stat, "Stat. (34 pb^{-1})", "F");

  Float_t yaxis_max = (resp_reso=="resp") ? 10. : 50.;
  if( resp_reso=="resp" && RECOTYPE_=="calo" ) yaxis_max = 20.;

  TH2D* axes = new TH2D("axes", "", 10, 20., 500., 10, 0., yaxis_max);
  axes->GetXaxis()->SetTitleOffset(1.1);
  axes->GetYaxis()->SetTitleOffset(1.2);
  axes->SetXTitle("Photon p_{T} [GeV/c]");
  axes->SetYTitle("Relative Uncertainty [%]");
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();

  
  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  c1->SetLogx();
  axes->Draw();
  syst_total->Draw("Csame");
  syst_q->Draw("L same");
  syst_photonScale->Draw("L same");
  syst_flavor->Draw("L same");
  //syst_QCD->Draw("L same");
  //stat->Draw("csame");
  gPad->RedrawAxis();
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");


  std::string canvasName = "totalSyst_" + resp_reso + "_" + JETALGO_ + ".eps";
  c1->SaveAs(canvasName.c_str());

  std::string file_total_name = "totalSyst_"+resp_reso+"_" + JETALGO_ + ".root";
  TFile* file_total = TFile::Open(file_total_name.c_str(), "recreate");
  file_total->cd();
  syst_total->Write();
  syst_q->Write();
  syst_photonScale->Write();
  syst_flavor->Write();
  stat->Write();
  file_total->Close();
  

  delete syst_total;
  delete legend;
  delete c1;
  delete axes;
 
}
