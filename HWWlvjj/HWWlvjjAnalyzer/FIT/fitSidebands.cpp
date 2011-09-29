#include <cstdlib>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "CommonTools/DrawBase.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace RooFit;



void fitSidebands( TTree* treeMC, TTree* treeDATA, int leptType);


int main() {


  TFile* file_DY         = TFile::Open("HWWlvjj_DY_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  TFile* file_GJet       = TFile::Open("HWWlvjj_GJet_TuneZ2_7TeV-alpgen_helicity_ALL.root");
  TFile* file_QCD_BCtoE  = TFile::Open("HWWlvjj_QCD_BCtoE_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  TFile* file_EMEnriched = TFile::Open("HWWlvjj_QCD_EMEnriched_TuneZ2_7TeV-pythia6_3_helicity_ALL.root");
  TFile* file_QCD_Mu     = TFile::Open("HWWlvjj_QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  TFile* file_TTJets     = TFile::Open("HWWlvjj_TTJets_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root");
  TFile* file_TTOBLNu    = TFile::Open("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_helicity_ALL.root");
  TFile* file_VV         = TFile::Open("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_helicity_ALL.root");
  TFile* file_WJets      = TFile::Open("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root");



  TFile* file_DATA = TFile::Open("HWWlvjj_DATA_6july_helicity_ALL.root");
  TTree* treeDATA = (TTree*)file_DATA->Get("Tree_FITUL");

  TChain* chainMC_ele = new TChain("Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_DY_TuneZ2_7TeV-pythia6_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_GJet_TuneZ2_7TeV-alpgen_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_QCD_BCtoE_TuneZ2_7TeV-pythia6_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_QCD_EMEnriched_TuneZ2_7TeV-pythia6_3_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_TTJets_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_helicity_ALL.root/Tree_FITUL");
  chainMC_ele->Add("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root/Tree_FITUL");

  TChain* chainMC_mu = new TChain("Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_DY_TuneZ2_7TeV-pythia6_helicity_ALL.root/Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_GJet_TuneZ2_7TeV-alpgen_helicity_ALL.root/Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_helicity_ALL.root/Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_TTJets_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root/Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_helicity_ALL.root/Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_helicity_ALL.root/Tree_FITUL");
  chainMC_mu->Add("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root/Tree_FITUL");


  TTree* treeDATA_mu  = treeDATA->CopyTree("leptType==0");
  TTree* treeDATA_ele = treeDATA->CopyTree("leptType==1");

  TTree* treeMC_mu =  chainMC_mu->CopyTree("leptType==0");
  TTree* treeMC_ele = chainMC_ele->CopyTree("leptType==1");

  fitSidebands( treeMC_ele, treeDATA_ele, 1);
  fitSidebands( treeMC_mu, treeDATA_mu, 0);

  return 0;

}






void fitSidebands( TTree* treeMC, TTree* treeDATA, int leptType ) {

  std::string leptType_str;
  if( leptType==0 ) leptType_str="mu";
  else if( leptType==1 ) leptType_str="ele";
  else {
    std::cout << "Unkown lepton type: '" << leptType << "'. Exiting." << std::endl;
  }

  std::string outfilename = "alphaFile_" + leptType_str + ".root";
  TFile* outfile = TFile::Open(outfilename.c_str(), "recreate");
  outfile->cd();

  char ofs_name[400];
  sprintf( ofs_name, "FitSidebands/fitresults_%s.txt", leptType_str.c_str());
  ofstream ofs(ofs_name);



  char cut_sidebands[500];
  sprintf( cut_sidebands, "( (mJJ>40. && mJJ<60.)||(mJJ>100. && mJJ<160.) )");
  char cut_signal[500];
  sprintf( cut_signal, "( mJJ>60. && mJJ<100. )");
  

  float mWW_min = 240.;
  float mWW_max = 1000.;
  float binWidth = 20.;
  int nBins = (int)(mWW_max-mWW_min)/binWidth;

  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  //RooRealVar* eventWeight_alpha = new RooRealVar("eventWeight_alpha", "event weight (alpha corrected)", 0., 2., "");
  RooRealVar* mWW = new RooRealVar("mWW", "m_{lvjj}", mWW_min, mWW_max, "GeV");
  RooRealVar* mJJ = new RooRealVar("mJJ", "mJJ", 60., 130., "GeV");

  RooFormulaVar* weight_lumi = new RooFormulaVar("weight_lumi", "@0*1000.", RooArgList(*eventWeight));

  RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",treeMC,RooArgSet(*eventWeight,*mWW,*mJJ),cut_sidebands,"eventWeight");
  RooDataSet signalMC("signalMC","signalMC",treeMC,RooArgSet(*eventWeight,*mWW,*mJJ),cut_signal,"eventWeight");

  RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",treeDATA,RooArgSet(*eventWeight,*mWW,*mJJ),cut_sidebands);
  RooDataSet signalDATA("signalDATA","signalDATA",treeDATA,RooArgSet(*eventWeight,*mWW,*mJJ),cut_signal);



  // -------------------- exponential ---------------------------
  RooRealVar a_exp("a_exp","a_exp",-0.001, -2., 1.);
  RooExponential exp("exp","exp",*mWW,a_exp);
  


  // FIRST: fit MC sidebands:

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FIRST STEP: FIT MC SIDEBANDS " << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsMC,SumW2Error(kTRUE));
  TH1D* h1_expSlope_sidebandsMC = new TH1D("expSlope_sidebandsMC", "", 1, 0., 1.);
  h1_expSlope_sidebandsMC->SetBinContent(1, a_exp.getVal());
  h1_expSlope_sidebandsMC->SetBinError(1, a_exp.getError());


  // save value/error of parameter:
  float a_exp_sidebandsMC = a_exp.getVal();
  float a_exp_sidebandsMC_error = a_exp.getError();
  ofs << "a_sidebandsMC: \t" << a_exp_sidebandsMC << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_sidebandsMC = mWW->frame();

  //plotOnFrame(plot_sidebandsMC, &sidebandsMC, nBins, mWW, &a_exp, &exp);

  sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));

  RooRealVar* a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  RooExponential* exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mWW,*a_exp_plusSigma);
  exp_plusSigma->plotOn(plot_sidebandsMC,LineColor(38),LineStyle(2));

  RooRealVar* a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  RooExponential* exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mWW,*a_exp_minusSigma);
  exp_minusSigma->plotOn(plot_sidebandsMC,LineColor(38),LineStyle(2));

  exp.plotOn(plot_sidebandsMC, LineColor(kRed));
  sidebandsMC.plotOn(plot_sidebandsMC, Binning(nBins));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  plot_sidebandsMC->Draw();

  char canvasName[500];
  sprintf( canvasName, "FitSidebands/mWW_sidebandsMC_%s", leptType_str.c_str());
  std::string* canvasName_str = new std::string(canvasName);
  std::string canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  // SECOND: fit MC signal region:

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  SECOND STEP: FIT MC SIGNAL " << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;

  c1->Clear();
  c1->SetLogy(false);


  RooFitResult *r_signalMC = exp.fitTo(signalMC,SumW2Error(kTRUE));
  TH1D* h1_expSlope_signalMC = new TH1D("expSlope_signalMC", "", 1, 0., 1.);
  h1_expSlope_signalMC->SetBinContent(1, a_exp.getVal());
  h1_expSlope_signalMC->SetBinError(1, a_exp.getError());


  // save value/error of parameter:
  float a_exp_signalMC = a_exp.getVal();
  float a_exp_signalMC_error = a_exp.getError();
  ofs << "a_signalMC: \t" << a_exp_signalMC << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_signalMC = mWW->frame();

  //plotOnFrame(plot_signalMC, &signalMC, nBins, mWW, &a_exp, &exp);

  signalMC.plotOn(plot_signalMC, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mWW,*a_exp_plusSigma);
  exp_plusSigma->plotOn(plot_signalMC,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mWW,*a_exp_minusSigma);
  exp_minusSigma->plotOn(plot_signalMC,LineColor(38),LineStyle(2));

  exp.plotOn(plot_signalMC, LineColor(kRed));
  signalMC.plotOn(plot_signalMC, Binning(nBins));


  plot_signalMC->Draw();

  sprintf( canvasName, "FitSidebands/mWW_signalMC_%s", leptType_str.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  std::string cut_sidebands_weight(cut_sidebands);
  cut_sidebands_weight = "eventWeight*"+cut_sidebands_weight;

  std::string cut_signal_weight(cut_signal);
  cut_signal_weight = "eventWeight*"+cut_signal_weight;

  TH1D* h1_sidebandsMC = new TH1D("sidebandsMC", "", nBins, mWW_min, mWW_max);
  treeMC->Project("sidebandsMC", "mWW", cut_sidebands_weight.c_str());

  TH1D* h1_signalMC = new TH1D("signalMC", "", nBins, mWW_min, mWW_max);
  treeMC->Project("signalMC", "mWW", cut_signal_weight.c_str());

  h1_expSlope_signalMC->Write();
  h1_expSlope_sidebandsMC->Write();
  h1_signalMC->Write();
  h1_sidebandsMC->Write();
  outfile->Close();


/*
  // THIRD: define alpha:

  float alpha = a_exp_sidebandsMC - a_exp_signalMC;

  // FOURTH: fit DATA sidebands:

  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  THIRD STEP: FIT DATA SIDEBANDS (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA);
  RooFitResult *r_sidebandsDATA2 = landau_exp.fitTo(sidebandsDATA);
  //RooFitResult *r_sidebandsDATA = exp.fitTo(sidebandsDATA,SumW2Error(kFALSE),InitialHesse(kTRUE),Save());

  // save value/error of parameter:
  float a_exp_sidebandsDATA = a_exp.getVal();
  float a_exp_sidebandsDATA_error = a_exp.getError();
  ofs << "a_sidebandsDATA: \t" << a_exp_sidebandsDATA << "\t+-" << a_exp.getError() << std::endl;

  RooPlot *plot_sidebandsDATA = mWW->frame();


  sidebandsDATA.plotOn(plot_sidebandsDATA, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp.getVal()+a_exp.getError(), -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mWW,*a_exp_plusSigma);
  //exp_plusSigma->plotOn(plot_sidebandsDATA,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp.getVal()-a_exp.getError(), -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mWW,*a_exp_minusSigma);
  //exp_minusSigma->plotOn(plot_sidebandsDATA,LineColor(38),LineStyle(2));

  //exp.plotOn(plot_sidebandsDATA, LineColor(kRed));
  landau_exp.plotOn(plot_sidebandsDATA, LineColor(kGreen));
  sidebandsDATA.plotOn(plot_sidebandsDATA, Binning(nBins));

  plot_sidebandsDATA->Draw();

  sprintf( canvasName, "FitSidebands/mWW_sidebandsDATA_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());


  // FIFTH: scale data sidebands fit with alpha and superimpose to signal region:


  std::cout << std::endl << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "  FOURTH STEP: FIT DATA SIGNAL (" << btagCategory << " btags)" << std::endl;
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << std::endl << std::endl;


  c1->Clear();
  c1->SetLogy(false);

  RooRealVar a_exp_data("a_exp_data","a_exp_data", -alpha+a_exp_sidebandsDATA, -1., 0.);
  a_exp_data.setConstant(kTRUE);
  RooExponential exp_data("exp_data","exp_data",*mWW,a_exp_data);

  float a_exp_data_error = sqrt( a_exp_sidebandsMC_error*a_exp_sidebandsMC_error + 
                                 a_exp_signalMC_error*a_exp_signalMC_error +
                                 a_exp_sidebandsDATA_error*a_exp_sidebandsDATA_error );


  ofs << "a_signalDATA: \t" << a_exp_data.getVal() << "\t+-" << a_exp_data_error << std::endl;

  RooFitResult *r_signalDATA = exp_data.fitTo(signalDATA);

  RooPlot *plot_signalDATA = mWW->frame();
  signalDATA.plotOn(plot_signalDATA, Binning(nBins));

  a_exp_plusSigma = new RooRealVar("a_exp_plusSigma", "a_exp_plusSigma", a_exp_data.getVal()+a_exp_data_error, -1., 0.);
  exp_plusSigma = new RooExponential("exp_plusSigma", "exp_plusSigma",*mWW,*a_exp_plusSigma);
  //exp_plusSigma->plotOn(plot_signalDATA,LineColor(38),LineStyle(2));

  a_exp_minusSigma = new RooRealVar("a_exp_minusSigma", "a_exp_minusSigma", a_exp_data.getVal()-a_exp_data_error, -1., 0.);
  exp_minusSigma = new RooExponential("exp_minusSigma", "exp_minusSigma",*mWW,*a_exp_minusSigma);
  //exp_minusSigma->plotOn(plot_signalDATA,LineColor(38),LineStyle(2));

  //exp_data.plotOn(plot_signalDATA, LineColor(kRed));
  signalDATA.plotOn(plot_signalDATA, Binning(nBins));

  plot_signalDATA->Draw();

  sprintf( canvasName, "FitSidebands/mWW_signalDATA_%dbtag_%s", btagCategory, leptType.c_str());
  canvasName_str = new std::string(canvasName);
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  c1->SetLogy();
  *canvasName_str += "_log";
  canvasName_eps = *canvasName_str + ".eps";
  c1->SaveAs(canvasName_eps.c_str());

  std::cout << std::endl << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << " -- " << btagCategory << " btags" << std::endl;
  std::cout << " a_exp_sidebandsMC: " << a_exp_sidebandsMC << std::endl;
  std::cout << " a_exp_signalMC: " << a_exp_signalMC << std::endl;
  std::cout << " a_exp_sidebandsDATA: " << a_exp_sidebandsDATA << std::endl;
  std::cout << " Alpha: " << alpha << std::endl;
  std::cout << " alpha/a_sidebandsDATA: " << alpha/a_exp_sidebandsDATA << std::endl;
  std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout << std::endl << std::endl;

  ofs.close();

*/

  delete eventWeight;
  delete mWW;
  delete mJJ;
  delete c1;

}



