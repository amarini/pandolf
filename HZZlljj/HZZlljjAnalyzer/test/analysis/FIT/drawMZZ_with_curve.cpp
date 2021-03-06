#include <stdlib.h>
#include <iostream>
#include <string>
#include "DrawBase.h"
#include "fitTools.h"

#include "RooRealVar.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "TString.h"

#include "SidebandFitter.h"


bool QUICK_ = true;
bool root_aussi_ = true;

bool withSignal_ = true;
bool use_sherpa = false;




void drawHistoWithCurve( DrawBase* db, const std::string& data_dataset, const std::string& PUType, const std::string& data_mc, const std::string& flags, int nbtags, const std::string& leptType="ALL", string suffix="" );
void drawSidebandsWithCurve( DrawBase* db, const std::string& data_dataset, const std::string& PUType, const std::string& data_mc, const std::string& flags, int nbtags, const std::string& leptType="ALL", std::string suffix="" );




int main(int argc, char* argv[]) {

  if( argc!=3 && argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./drawMZZ_with_curve [(string)data_dataset] [(int)signalScaleFactor] [(string)data_mc] [(string)flags=\"\"]" << std::endl;
    exit(23);
  }

  std::string leptType = "ALL";

  std::string data_prefix(argv[1]);
  std::string data_dataset = "DATA_" + data_prefix;

  float signalScaleFactor=3.;
  if( argc>2 ) {
    std::string signalScaleFactor_str(argv[2]);
    signalScaleFactor = (float)atoi(signalScaleFactor_str.c_str());;
  }

  std::string data_mc="MC";
  if( argc>3 ) {
    std::string data_mc_tmp(argv[3]);
    data_mc = data_mc_tmp;
  }


  std::string flags = "";
  if( argc>4 ) {
    std::string flags_tmp(argv[4]);
    flags = "_" + flags_tmp;
  }


  std::string selType = "optLD_looseBTags_v2";

  TString dataset_tstr(data_prefix);
  std::string PUType = "Run2011A";
  if( data_prefix=="HR11" )
    PUType = "HR11";
  if( data_prefix=="HR11_v2" )
    PUType = "HR11_73pb";
  if( dataset_tstr.BeginsWith("Run2011B") )
    PUType = "HR11_73pb";





  DrawBase* db = new DrawBase("HZZlljjRM");
  db->set_pdf_aussi((bool)false);

  db->set_isCMSArticle(true);
  db->set_lumiOnRightSide(true);


  std::string outputdir_str = "HZZlljjRMPlots_" + data_dataset;
  if( withSignal_ ) {
    outputdir_str += "_plusSignal";
    if( signalScaleFactor!=1. ) {
      char scaleFactorText[100];
      sprintf( scaleFactorText, "_times%.0f", signalScaleFactor );
      std::string scaleFactorText_str(scaleFactorText);
      outputdir_str += scaleFactorText_str;
    }
  }
  outputdir_str += "_" + selType + "_PU" + PUType + "_" + leptType + "_fit" + data_mc;
  if( use_sherpa ) outputdir_str += "_sherpa";
  outputdir_str += flags;
  db->set_outputdir(outputdir_str);


  std::string dataFileName = "HZZlljjRM_" + data_dataset + "_"+selType+"_"+leptType+".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db->add_dataFile( dataFile, "THEDATA" );

  std::string signalFileName = "HZZlljjRM_GluGluToHToZZTo2L2Q_M-400_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1";
  //std::string signalFileName = "HZZlljjRM_GluGluToHToZZTo2L2Q_M-230_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1";
  signalFileName += "_" + selType;
  signalFileName += "_PU" + PUType;
  signalFileName += "_" + leptType;
  signalFileName += ".root";
  TFile* signalFile = TFile::Open(signalFileName.c_str());
  if( withSignal_ ) {
    char signalLegendText[400];
    if( signalScaleFactor==1. ) 
      sprintf( signalLegendText, "H(400 GeV)" );
    else
      sprintf( signalLegendText, "H(400 GeV) #times %.0f", signalScaleFactor);
    std::string signalLegendText_str(signalLegendText);
    db->add_mcFile( signalFile, signalScaleFactor, "H400", signalLegendText_str, kYellow, 3004);
    //db->add_mcFile_superimp( signalFile, "H400", signalLegendText_str, 5., kRed+3, 2 );
  }

  std::string mcZJetsFileName;
  mcZJetsFileName = "HZZlljjRM_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_PU" + PUType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 29, 3001);
  //db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 10, 0);



  std::string mcVVFileName = "HZZlljjRM_VV_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S4_START42_V11-v1";
  mcVVFileName += "_" + selType;
  mcVVFileName += "_PU" + PUType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  db->add_mcFile( mcVVFile, "VVtoAnything_TuneZ2", "ZZ/WZ/WW", 38, 3003);


  std::string mcTTbarFileName = "HZZlljjRM_TT_TW_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_PU" + PUType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TTtW", "t#bar{t}/tW", 39, 3002);






  if( data_dataset=="DATA_Run2011A_v2_Sub2" )
    db->set_lumiNormalization(175.);
  else if( data_dataset=="DATA_1fb" )
    db->set_lumiNormalization(859.);
  else if( data_dataset=="DATA_EPS" )
    db->set_lumiNormalization(960.); 
  else if( data_dataset=="DATA_EPS_FINAL" )
    db->set_lumiNormalization(1000.); 
  else if( data_dataset=="DATA_EPS_FINAL_FULL" )
    db->set_lumiNormalization(1143.); 
  else if( data_dataset=="DATA_EPS_FINAL_FULL_plusSingleMu" )
    db->set_lumiNormalization(1143.); 
  else if( data_dataset=="DoubleElectron_Aug05ReReco" )
    db->set_lumiNormalization(227.); 
  else if( data_dataset=="DoubleMu_Aug05ReReco" )
    db->set_lumiNormalization(285.); 
  else if( data_dataset=="DATA_EPS_FINAL_plusSingleMu" )
    db->set_lumiNormalization(1143.);
  else if( data_dataset=="DATA_LP11" )
    db->set_lumiNormalization(1580.);
  else if( data_dataset=="DATA_Run2011A_FULL" )
    db->set_lumiNormalization(2100.);
  else if( data_dataset=="DATA_HR11" )
    db->set_lumiNormalization(4200.);
  else if( data_dataset=="DATA_HR11_v2" )
    db->set_lumiNormalization(4600.);
  else if( data_dataset=="DATA_Run2011B_v2" )
    db->set_lumiNormalization(2500.);




  bool log = true;


  db->set_legendTitle("Gluon- and 0 b-tag");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags<=0)", 30, 150., 750., "mZZ_g0btag", "m_{ZZ}", "GeV");
  db->set_legendTitle("0 b-tag category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", 30, 150., 750., "mZZ_0btag", "m_{ZZ}", "GeV");
  db->set_legendTitle("Gluon-tag category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==-1)", 30, 150., 750., "mZZ_gtag", "m_{ZZ}", "GeV");


  float binWidth = 20.;
  float xMin = 183.;
  float xMax = 803.;
  int nBins = (int)((xMax-xMin)/binWidth);


  // signal box plots:
  db->set_legendTitle("0 b-tag category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", nBins, xMin, xMax, "mZZ_0btag", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 0);
  drawSidebandsWithCurve( db, data_prefix, PUType, data_mc, flags, 0, "ALL");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", 14, 183., 253., "mZZ_0btag_turnOnZOOM", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 0, "ALL", "turnOnZOOM");
  drawSidebandsWithCurve( db, data_prefix, PUType, data_mc, flags, 0, "ALL", "turnOnZOOM");
  db->set_legendTitle("0 b-tag category (muons)");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0 && leptType==0)", nBins, xMin, xMax, "mZZ_0btag_mm", "m_{#mu#mujj}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 0, "MU");
  db->set_legendTitle("0 b-tag category (electrons)");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0 && leptType==1)", nBins, xMin, xMax, "mZZ_0btag_ee", "m_{eejj}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 0, "ELE");

  db->set_legendTitle("1 b-tag category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1)", nBins, xMin, xMax, "mZZ_1btag", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 1 );
  drawSidebandsWithCurve( db, data_prefix, PUType, data_mc, flags, 1, "ALL");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1)", 14, 183., 253., "mZZ_1btag_turnOnZOOM", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 1, "ALL", "turnOnZOOM");
  drawSidebandsWithCurve( db, data_prefix, PUType, data_mc, flags, 1, "ALL", "turnOnZOOM");
  db->set_legendTitle("1 b-tag category (muons)");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1 && leptType==0)", nBins, xMin, xMax, "mZZ_0btag_mm", "m_{#mu#mujj}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 1, "MU");
  db->set_legendTitle("1 b-tag category (electrons)");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1 && leptType==1)", nBins, xMin, xMax, "mZZ_0btag_ee", "m_{eejj}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 1, "ELE");

  db->set_legendTitle("2 b-tag category");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2)", nBins, xMin, xMax, "mZZ_2btag", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 2 );
  drawSidebandsWithCurve( db, data_prefix, PUType, data_mc, flags, 2, "ALL");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2)", 14, 183., 253., "mZZ_2btag_turnOnZOOM", "m_{ZZ}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 2, "ALL", "turnOnZOOM");
  drawSidebandsWithCurve( db, data_prefix, PUType, data_mc, flags, 2, "ALL", "turnOnZOOM");
  db->set_legendTitle("2 b-tag category (muons)");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2 && leptType==0)", nBins, xMin, xMax, "mZZ_0btag_mm", "m_{#mu#mujj}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 2, "MU");
  db->set_legendTitle("2 b-tag category (electrons)");
  db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2 && leptType==1)", nBins, xMin, xMax, "mZZ_0btag_ee", "m_{eejj}", "GeV");
  drawHistoWithCurve( db, data_prefix, PUType, data_mc, flags, 2, "ELE");


  //xMin = 160.;
  //xMax = 1350.;
  //nBins = (int)((xMax-xMin)/binWidth);

  //// long range (up to 1300 gev):
  //db->set_legendTitle("0 b-tag category");
  //db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==0)", nBins, xMin, xMax, "mZZ_0btag_longRange", "m_{ZZ}", "GeV");
  //drawHistoWithCurve( db, data_prefix, PUType, 0, "longRange");

  //db->set_legendTitle("1 b-tag category");
  //db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==1)", nBins, xMin, xMax, "mZZ_1btag_longRange", "m_{ZZ}", "GeV");
  //drawHistoWithCurve( db, data_prefix, PUType, 1, "longRange");

  //db->set_legendTitle("2 b-tag category");
  //db->drawHisto_fromTree("tree_passedEvents", "CMS_hzz2l2q_mZZ", "eventWeight*(mZjj>75. && mZjj<105. && nBTags==2)", nBins, xMin, xMax, "mZZ_2btag_longRange", "m_{ZZ}", "GeV");
  //drawHistoWithCurve( db, data_prefix, PUType, 2, "longRange");


//db->set_legendTitle("0 b-tag Sidebands");
//db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag", "m_{lljj}", "GeV", "Events", log);

//db->set_legendTitle("1 b-tag Sidebands");
//db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag", "m_{lljj}", "GeV", "Events", log);

//db->set_legendTitle("2 b-tag Sidebands");
//db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag", "m_{lljj}", "GeV", "Events", log);



  delete db;
  db = 0;

  return 0;

}  





void drawHistoWithCurve( DrawBase* db, const std::string& data_dataset, const std::string& PUType, const std::string& data_mc, const std::string& flags, int nbtags, const std::string& leptType, std::string suffix ) {

  if( suffix!="" ) suffix = "_" + suffix;

  TH1F::AddDirectory(kTRUE);

  // get histograms:

  std::vector< TH1D* > lastHistos_data = db->get_lastHistos_data();
  std::vector< TH1D* > lastHistos_mc   = db->get_lastHistos_mc();
  std::vector< TH1D* > lastHistos_mc_superimp   = db->get_lastHistos_mc_superimp();


  TH1D* h1_data = new TH1D(*(lastHistos_data[0]));
  float xMin = (db->get_xAxisMin()!=9999.) ? db->get_xAxisMin() : h1_data->GetXaxis()->GetXmin();
  float xMax = (db->get_xAxisMax()!=9999.) ? db->get_xAxisMax() : h1_data->GetXaxis()->GetXmax();

  // create data graph (poisson asymm errors):
  TGraphAsymmErrors* graph_data_poisson = new TGraphAsymmErrors(0);
  graph_data_poisson = fitTools::getGraphPoissonErrors(h1_data);
  graph_data_poisson->SetMarkerStyle(20);

  THStack* mc_stack = new THStack();
  for( unsigned ihisto=0; ihisto<lastHistos_mc.size(); ++ihisto ) 
    mc_stack->Add(lastHistos_mc[lastHistos_mc.size()-ihisto-1]);



  std::string sherpa_suffix = (use_sherpa) ? "_sherpa" : "";


  // open fit results file:
  char fitResultsFileName[200];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_ALL_PU%s_fit%s%s%s.root", data_dataset.c_str(), nbtags, PUType.c_str(), data_mc.c_str(), sherpa_suffix.c_str(), flags.c_str() );
  TFile* fitResultsFile = TFile::Open(fitResultsFileName);

  // get bg workspace:
  char workspaceName[200];
  sprintf( workspaceName, "fitWorkspace_%dbtag", nbtags );
  RooWorkspace* bgws = (RooWorkspace*)fitResultsFile->Get(workspaceName);
  char fitResultsName[200];
  sprintf( fitResultsName, "fitResults_%dbtag_decorr", nbtags );
  RooFitResult* bgfr = (RooFitResult*)fitResultsFile->Get(fitResultsName);

  // get mZZ variable:
  RooRealVar* CMS_hzz2l2q_mZZ = (RooRealVar*)bgws->var("CMS_hzz2l2q_mZZ");

  // get bg shape:
  RooAbsPdf* background = (RooAbsPdf*)bgws->pdf("background_decorr");


  std::string leptType_cut;
  if( leptType=="MU" ) leptType_cut = "&& leptType==0";
  if( leptType=="ELE" ) leptType_cut = "&& leptType==1";


  SidebandFitter* sf = new SidebandFitter( data_dataset, PUType, data_mc, flags );

  // get bg normalization:
  //float expBkg = sf->get_backgroundNormalization( nbtags, leptType, "DATA" );
  std::pair<Double_t, Double_t> bgNormAndError = sf->get_backgroundNormalizationAndError( nbtags, leptType, "DATA" );
  float expBkg = bgNormAndError.first;
  float expBkg_err = bgNormAndError.second;


  //TTree* treeSidebandsDATA_alphaCorr = (TTree*)fitResultsFile->Get("sidebandsDATA_alpha");
  //TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, xMin, xMax);
  //char sidebandsCut_alpha[500];
  //sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d %s)", nbtags, leptType_cut.c_str());
  //treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
  //float expBkg = h1_mZZ_sidebands_alpha->Integral();


  RooPlot *plot_MCbkg = CMS_hzz2l2q_mZZ->frame(xMin,xMax,(int)(xMax-xMin)/h1_data->GetXaxis()->GetBinWidth(1));
  background->plotOn(plot_MCbkg,RooFit::Normalization(expBkg));
  if( suffix == "_turnOnZOOM" && !QUICK_ ) {
    background->plotOn(plot_MCbkg,RooFit::VisualizeError(*bgfr,2.0,kFALSE),RooFit::FillColor(kYellow), RooFit::Normalization(expBkg));
    background->plotOn(plot_MCbkg,RooFit::VisualizeError(*bgfr,1.0,kFALSE),RooFit::FillColor(kGreen), RooFit::Normalization(expBkg));
    background->plotOn(plot_MCbkg,RooFit::Normalization(expBkg));
    background->plotOn(plot_MCbkg,RooFit::LineStyle(2),RooFit::Normalization(expBkg+expBkg_err));
    background->plotOn(plot_MCbkg,RooFit::LineStyle(2),RooFit::Normalization(expBkg-expBkg_err));
  }

  TF1* f1_bgForLegend = new TF1("bgForLegend", "[0]");
  f1_bgForLegend->SetLineColor(kBlue);
  f1_bgForLegend->SetLineWidth(3);
  

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.3*h1_data->GetMaximum());
  char yTitle[200];
  sprintf( yTitle, "Events / (%.0f GeV)", h1_data->GetXaxis()->GetBinWidth(1) );
  h2_axes->SetYTitle(yTitle);
  if( leptType=="MU" )
    h2_axes->SetXTitle("m_{#mu#mujj} [GeV]");
  else if( leptType=="ELE" )
    h2_axes->SetXTitle("m_{eejj} [GeV]");
  else
    h2_axes->SetXTitle("m_{ZZ} [GeV]");

  float legend_xMin = 0.38;
  float legend_yMax = 0.91;
  float legend_yMin = legend_yMax - 0.07*6.;
  float legend_xMax = 0.92;

  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( graph_data_poisson, "Data", "P");
  legend->AddEntry( f1_bgForLegend, "Expected background", "L");
  for( unsigned imc=0; imc<lastHistos_mc.size(); ++imc ) 
    legend->AddEntry( lastHistos_mc[imc], (db->get_mcFile(imc).legendName).c_str(), "F");

  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  if( suffix != "_turnOnZOOM" )
    legend->Draw("same");
  mc_stack->Draw("histo same");
  for( unsigned imc=0; imc<lastHistos_mc_superimp.size(); ++imc ) 
    lastHistos_mc_superimp[imc]->Draw("same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();


  
  char canvasName[1000];
  if( leptType=="ALL" )
    sprintf( canvasName, "%s/mZZ_%dbtag_withCurve%s", (db->get_outputdir()).c_str(), nbtags, suffix.c_str() );
  else 
    sprintf( canvasName, "%s/mZZ_%dbtag_withCurve%s_%s", (db->get_outputdir()).c_str(), nbtags, suffix.c_str(), leptType.c_str() );
  std::string canvasName_str(canvasName);
  std::string canvasName_eps = canvasName_str+".eps";
  c1->SaveAs(canvasName_eps.c_str());
  std::string canvasName_root = canvasName_str+".root";
  if( root_aussi_ ) {
    c1->SetFixedAspectRatio(true);
    c1->SaveAs(canvasName_root.c_str());
  }

  c1->Clear();


  legend_yMin = legend_yMax - 0.07*4.;

  TLegend* legend2 = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend2->SetTextSize(0.04);
  legend2->SetFillColor(0);
  legend2->AddEntry( graph_data_poisson, "Data", "P");
  legend2->AddEntry( f1_bgForLegend, "Expected background", "L");
  for( unsigned imc=0; imc<lastHistos_mc.size(); ++imc ) {
    if( db->get_mcFile(imc).fillColor==kYellow ) //signal only
      legend2->AddEntry( lastHistos_mc[imc], (db->get_mcFile(imc).legendName).c_str(), "F");
  }

  h2_axes->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  legend2->Draw("same");
  for( unsigned imc=0; imc<lastHistos_mc.size(); ++imc ) {
    if( db->get_mcFile(imc).fillColor==kYellow ) //signal only
      lastHistos_mc[imc]->Draw("same");
  }
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();


  
  if( leptType=="ALL" )
    sprintf( canvasName, "%s/mZZ_%dbtag_withCurve_noMC%s", (db->get_outputdir()).c_str(), nbtags, suffix.c_str() );
  else 
    sprintf( canvasName, "%s/mZZ_%dbtag_withCurve_noMC%s_%s", (db->get_outputdir()).c_str(), nbtags, suffix.c_str(), leptType.c_str() );
  std::string canvasName2_str(canvasName);
  std::string canvasName2_eps = canvasName2_str+".eps";
  c1->SaveAs(canvasName2_eps.c_str());

  c1->Clear();






  c1->SetLogy();

  TLegend* legend_log = new TLegend(0.6, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend_log->SetTextSize(0.04);
  legend_log->SetFillColor(0);
  legend_log->AddEntry( graph_data_poisson, "Data", "P");
  legend_log->AddEntry( f1_bgForLegend, "Exp. BG", "L");
  for( unsigned imc=0; imc<lastHistos_mc.size(); ++imc ) 
    legend_log->AddEntry( lastHistos_mc[imc], (db->get_mcFile(imc).legendName).c_str(), "F");



  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.05, 200.*h1_data->GetMaximum());
  h2_axes_log->SetYTitle(yTitle);
  h2_axes_log->SetXTitle("m_{ZZ} [GeV]");
  h2_axes_log->GetYaxis()->SetNoExponent();

  h2_axes_log->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  legend_log->Draw("same");
  mc_stack->Draw("histo same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();

  std::string canvasName_log_eps = canvasName_str+"_log.eps";
  c1->SaveAs(canvasName_log_eps.c_str());

  delete sf;
  delete h2_axes;
  delete h2_axes_log;
  delete c1;

}



void drawSidebandsWithCurve( DrawBase* db, const std::string& data_dataset, const std::string& PUType, const std::string& data_mc, const std::string& flags, int nbtags, const std::string& leptType, std::string suffix ) {

  if( suffix!="" ) suffix = "_" + suffix;

  TH1F::AddDirectory(kTRUE);

  // get histograms:

  std::vector< TH1D* > lastHistos_data = db->get_lastHistos_data();


  TH1D* h1_data = new TH1D(*(lastHistos_data[0]));
  float xMin = (db->get_xAxisMin()!=9999.) ? db->get_xAxisMin() : h1_data->GetXaxis()->GetXmin();
  float xMax = (db->get_xAxisMax()!=9999.) ? db->get_xAxisMax() : h1_data->GetXaxis()->GetXmax();
  int nBins = (int)(xMax-xMin)/h1_data->GetXaxis()->GetBinWidth(1);




  std::string sherpa_suffix = (use_sherpa) ? "_sherpa" : "";


  // open fit results file:
  char fitResultsFileName[200];
  sprintf( fitResultsFileName, "fitResultsFile_%s_%dbtag_ALL_PU%s_fit%s%s%s.root", data_dataset.c_str(), nbtags, PUType.c_str(), data_mc.c_str(), sherpa_suffix.c_str(), flags.c_str() );
  TFile* fitResultsFile = TFile::Open(fitResultsFileName);

  // get bg workspace:
  char workspaceName[200];
  sprintf( workspaceName, "fitWorkspace_%dbtag", nbtags );
  RooWorkspace* bgws = (RooWorkspace*)fitResultsFile->Get(workspaceName);
  char fitResultsName[200];
  sprintf( fitResultsName, "fitResults_%dbtag_decorr", nbtags );
  RooFitResult* bgfr = (RooFitResult*)fitResultsFile->Get(fitResultsName);

  // get mZZ variable:
  RooRealVar* CMS_hzz2l2q_mZZ = (RooRealVar*)bgws->var("CMS_hzz2l2q_mZZ");

  // get bg shape:
  RooAbsPdf* background = (RooAbsPdf*)bgws->pdf("background_decorr");


  // get data sidebands:
  TTree* tree_sidebands = (TTree*)fitResultsFile->Get("sidebandsDATA_alpha");
  TH1D* h1_dataSidebands = new TH1D("dataSidebands", "", nBins, xMin, xMax);
  h1_dataSidebands->Sumw2();

  char selection[900];
  if( leptType=="ALL" )
    sprintf( selection, "eventWeight_alpha*(((mZjj>60.&&mZjj<75.)||(mZjj>105.&&mZjj<130.)) && nBTags==%d && CMS_hzz2l2q_mZZ>183. && CMS_hzz2l2q_mZZ<800.)", nbtags);
  else {
    int leptType_int = SidebandFitter::convert_leptType(leptType);
    sprintf( selection, "eventWeight_alpha*(((mZjj>60.&&mZjj<75.)||(mZjj>105.&&mZjj<130.)) && nBTags==%d && leptType==%d && CMS_hzz2l2q_mZZ>183. && CMS_hzz2l2q_mZZ<800.)", nbtags, leptType_int);
  }
  tree_sidebands->Project("dataSidebands", "CMS_hzz2l2q_mZZ", selection);



  // create data graph (poisson asymm errors):
  TGraphAsymmErrors* graph_data_poisson = new TGraphAsymmErrors(0);
  graph_data_poisson = fitTools::getGraphPoissonErrors(h1_dataSidebands);
  graph_data_poisson->SetMarkerStyle(20);


  std::string leptType_cut;
  if( leptType=="MU" ) leptType_cut = "&& leptType==0";
  if( leptType=="ELE" ) leptType_cut = "&& leptType==1";


  SidebandFitter* sf = new SidebandFitter( data_dataset, PUType, data_mc, flags );

  // get bg normalization:
  //float expBkg = sf->get_backgroundNormalization( nbtags, leptType, "DATA" );
  std::pair<Double_t, Double_t> bgNormAndError = sf->get_backgroundNormalizationAndError( nbtags, leptType, "DATA" );
  float expBkg = bgNormAndError.first;
  float expBkg_err = bgNormAndError.second;


  //TTree* treeSidebandsDATA_alphaCorr = (TTree*)fitResultsFile->Get("sidebandsDATA_alpha");
  //TH1D* h1_mZZ_sidebands_alpha = new TH1D("mZZ_sidebands_alpha", "", 65, xMin, xMax);
  //char sidebandsCut_alpha[500];
  //sprintf(sidebandsCut_alpha, "eventWeight_alpha*(isSidebands && nBTags==%d %s)", nbtags, leptType_cut.c_str());
  //treeSidebandsDATA_alphaCorr->Project("mZZ_sidebands_alpha", "CMS_hzz2l2q_mZZ", sidebandsCut_alpha);
  //float expBkg = h1_mZZ_sidebands_alpha->Integral();


  RooPlot *plot_MCbkg = CMS_hzz2l2q_mZZ->frame(xMin,xMax,nBins);
  background->plotOn(plot_MCbkg,RooFit::Normalization(expBkg));
  if( !QUICK_ ) {
    background->plotOn(plot_MCbkg,RooFit::VisualizeError(*bgfr,2.0,kFALSE),RooFit::FillColor(kYellow), RooFit::Normalization(expBkg));
    background->plotOn(plot_MCbkg,RooFit::VisualizeError(*bgfr,1.0,kFALSE),RooFit::FillColor(kGreen), RooFit::Normalization(expBkg));
    background->plotOn(plot_MCbkg,RooFit::Normalization(expBkg));
  }

  TF1* f1_bgForLegend = new TF1("bgForLegend", "[0]");
  f1_bgForLegend->SetLineColor(kBlue);
  f1_bgForLegend->SetLineWidth(3);
  

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.3*h1_data->GetMaximum());
  char yTitle[200];
  sprintf( yTitle, "Events / (%.0f GeV)", h1_data->GetXaxis()->GetBinWidth(1) );
  h2_axes->SetYTitle(yTitle);
  if( leptType=="MU" )
    h2_axes->SetXTitle("m_{#mu#mujj} [GeV]");
  else if( leptType=="ELE" )
    h2_axes->SetXTitle("m_{eejj} [GeV]");
  else
    h2_axes->SetXTitle("m_{ZZ} [GeV]");

  float legend_xMin = 0.5;
  float legend_yMax = 0.85;
  float legend_yMin = legend_yMax - 0.12*2.;
  float legend_xMax = 0.92;

  char legendTitle[300];
  sprintf( legendTitle, "%d b-tag sidebands", nbtags );

  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax, legendTitle);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( graph_data_poisson, "Data", "P");
  legend->AddEntry( f1_bgForLegend, "Fit", "L");

  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();


  TCanvas* c1 = new TCanvas( "c1", "", 600, 600);
  c1->cd();

  h2_axes->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  if( suffix != "_turnOnZOOM" )
    legend->Draw("same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();


  
  char canvasName[1000];
  if( leptType=="ALL" )
    sprintf( canvasName, "%s/mZZsidebands_%dbtag_withCurve%s", (db->get_outputdir()).c_str(), nbtags, suffix.c_str() );
  else 
    sprintf( canvasName, "%s/mZZsidebands_%dbtag_withCurve%s_%s", (db->get_outputdir()).c_str(), nbtags, suffix.c_str(), leptType.c_str() );
  std::string canvasName_str(canvasName);
  std::string canvasName_eps = canvasName_str+".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->Clear();

  c1->SetLogy();


  TLegend* legend_log = new TLegend(0.6, legend_yMin, legend_xMax, legend_yMax, (db->get_legendTitle()).c_str());
  legend_log->SetTextSize(0.04);
  legend_log->SetFillColor(0);
  legend_log->AddEntry( graph_data_poisson, "Data", "P");
  legend_log->AddEntry( f1_bgForLegend, "Exp. BG", "L");



  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.05, 200.*h1_data->GetMaximum());
  h2_axes_log->SetYTitle(yTitle);
  h2_axes_log->SetXTitle("m_{ZZ} [GeV]");

  h2_axes_log->Draw();
  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  //legend_log->Draw("same");
  legend->Draw("same");
  plot_MCbkg->Draw("same");
  graph_data_poisson->Draw("P same");

  gPad->RedrawAxis();

  std::string canvasName_log_eps = canvasName_str+"_log.eps";
  c1->SaveAs(canvasName_log_eps.c_str());


  delete sf;
  delete h2_axes;
  delete h2_axes_log;
  delete c1;

}



