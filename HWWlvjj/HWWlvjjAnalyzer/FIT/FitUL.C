#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <cstring>
#include <string>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TChain.h>
#include <TMath.h>
#include <TLegend.h>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace RooFit;

void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog);
TTree* correctTreeWithAlpha( TTree* tree, TF1* h1_alpha, const std::string& name );

void FitUL(int leptType){

  std::string leptType_str = (leptType==0) ? "mu" : "ele";

  gStyle->SetOptStat(1111111);
  gStyle->SetPalette(1);
  gStyle->SetOptFit(111110);
  gStyle->SetOptFile(1);
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.3);
  gStyle->SetMarkerColor(1);
  gStyle->SetFillColor(0);

  TCanvas* myc1 = new TCanvas("myc1", "myc1", 600, 600);
 
  TH1F* mWW_peak = new TH1F("mWW_peak","Before Correction",56, 140., 700.);
  mWW_peak->Sumw2();
  mWW_peak->SetXTitle("WW Invariant Mass [GeV]");
  TH1F* mWW_side = new TH1F("mWW_side","Before Correction",56, 140., 700.);
  mWW_side->Sumw2();
  mWW_side->SetXTitle("WW Invariant Mass [GeV]");
  TLegend *leg = new TLegend( 0.5,0.7,0.8,0.8 );
  leg->AddEntry(mWW_side,"Side-bands region","F");
  leg->AddEntry(mWW_peak,"Peak region","F");

  TH1F* mJJ_tot = new TH1F("mJJ_tot","",100, 0., 250.);
  mJJ_tot->Sumw2();


  // HISTO with Xaxis log
  double Bin[30]={1.};
  double* BinP;
  BinP=Bin;
  getBins_int( 10, BinP, 140., 700., true);
  TH1F* mWW_peak_log = new TH1F("mWW_peak_log","", 9, BinP);
  mWW_peak_log->Sumw2();
  TH1F* mWW_side_log = new TH1F("mWW_side_log","", 9, BinP);
  mWW_side_log->Sumw2();

  TChain chMC("Tree_FITUL");
  chMC.Add("HWWlvjj_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root");
  chMC.Add("HWWlvjj_DY_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  chMC.Add("HWWlvjj_VVtoAnything_TuneZ2_7TeV-pythia6-tauola_helicity_ALL.root");
  chMC.Add("HWWlvjj_TTJets_TuneZ2_7TeV-madgraph-tauola_helicity_ALL.root");
  chMC.Add("HWWlvjj_TToBLNu_TuneZ2_7TeV-madgraph_helicity_ALL.root");
  chMC.Add("HWWlvjj_QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  chMC.Add("HWWlvjj_QCD_EMEnriched_TuneZ2_7TeV-pythia6_3_helicity_ALL.root");
  chMC.Add("HWWlvjj_QCD_BCtoE_TuneZ2_7TeV-pythia6_helicity_ALL.root");
  chMC.Add("HWWlvjj_GJet_TuneZ2_7TeV-alpgen_helicity_ALL.root");
  chMC.Draw("mWW");

  float mJJ_, mWW_, eventWeight_;
  int leptType_;
  chMC.SetBranchAddress("mJJ", &mJJ_);
  chMC.SetBranchAddress("mWW", &mWW_);
  chMC.SetBranchAddress("eventWeight", &eventWeight_);
  chMC.SetBranchAddress("leptType", &leptType_);

  for( int iEntry=0; iEntry<chMC.GetEntries() ; iEntry++ ){
   chMC.GetEntry(iEntry);
  //std::cout<<iEntry<<"/"<<ch.GetEntries()<<":  "<< mJJ <<", "<< mWW <<std::endl;
   if( leptType_!=leptType ) continue;
   if( mJJ_>60. && mJJ_<100. )                          {  mWW_peak->Fill(mWW_,eventWeight_); mWW_peak_log->Fill(mWW_,eventWeight_); }
   if( (mJJ_>40. && mJJ_<60. ) || (mJJ_>100. && mJJ_<160.) ) {  mWW_side->Fill(mWW_,eventWeight_); mWW_side_log->Fill(mWW_,eventWeight_); } 
   mJJ_tot->Fill(mJJ_,eventWeight_);
   //if(mWW>600 && mWW<650 && ( mJJ>60. && mJJ<100 ) )                          { std::cout<<"IN "<<mJJ<<"  "<<eventWeight<<std::endl; }
   //if(mWW>600 && mWW<650 && ( (mJJ>40 && mJJ<60 ) || (mJJ>100 && mJJ<160) ) ) { std::cout<<"OUT"<<mJJ<<"  "<<eventWeight<<std::endl; }
  }

  // M(JJ)
  myc1->cd();
  mJJ_tot->Draw();
  myc1->SaveAs("mJJ.eps");
  // M(WW) COMPARISON
  mWW_side->SetLineColor(2);
  mWW_side->Scale( 1/mWW_side->Integral() );
  mWW_peak->Scale( 1/mWW_peak->Integral() );
  mWW_peak->Draw("HISTO");
  mWW_side->Draw("sameHISTO");
  leg->Draw("same");
  myc1->SaveAs("mWW.eps");

  // M(WW) RATIO
  TH1F* peak_side = new TH1F( *mWW_peak );
  peak_side->Divide( mWW_side );
  peak_side->Draw();
  myc1->SaveAs("peak_side.eps");
  // M(WW)_LOG peak & side
  mWW_peak_log->Draw();
  myc1->SaveAs("mWW_peak_log.eps");
  mWW_side_log->Draw();
  myc1->SaveAs("mWW_side_log.eps");
  // M(WW) RATIO xLOG

  myc1->cd();
  TH1F* peak_side_log = new TH1F( *mWW_peak_log );
  peak_side_log->Sumw2();
  peak_side_log->SetXTitle("WW Invariant Mass [GeV]");
  peak_side_log->SetYTitle("Ratio peak/side");
  peak_side_log->Divide( mWW_side_log );
  peak_side_log->Draw();
  myc1->SaveAs("peak_side_log.eps");

  // NOW FIT
  TF1 *myfit = new TF1("myfit","[0]+[1]*x+[2]*x^2+[3]*x^3", 140, 700);
  myfit->SetParName(0,"c0");
  myfit->SetParName(1,"c1");
  myfit->SetParName(2,"c2");
  myfit->SetParName(3,"c3");
  peak_side_log->Fit("myfit");
  peak_side_log->Draw();
  myc1->SaveAs("Ratio_fitted.eps");
    
  // Now I correct my distribution
  TH1F* mWW_peak_log2 = new TH1F("mWW_peak_log2","After Correction", 56, 140., 700);
  mWW_peak_log2->Sumw2();
  mWW_peak_log2->SetXTitle("WW Invariant Mass [GeV]");
  TH1F* mWW_side_log2 = new TH1F("mWW_side_log2","After Correction", 56, 140., 700) ;
  mWW_side_log2->Sumw2();
  mWW_side_log2->SetXTitle("WW Invariant Mass [GeV]");
  TLegend *leg_CORR = new TLegend( 0.5,0.7,0.8,0.8 );
  leg_CORR->AddEntry(mWW_side_log2,"Side-bands region","F");
  leg_CORR->AddEntry(mWW_peak_log2,"Peak region","F");

  for( int iEntry=0; iEntry<chMC.GetEntries() ; iEntry++ ){
   chMC.GetEntry(iEntry);
   if( leptType_!=leptType ) continue;
   if( mJJ_>60. && mJJ_<100. )                          {  mWW_peak_log2->Fill( mWW_,eventWeight_ ); }
   if( (mJJ_>40. && mJJ_<60. ) || (mJJ_>100. && mJJ_<160.) ) {  mWW_side_log2->Fill( mWW_,eventWeight_*myfit->Eval(mWW_)  );} //( 1.26421 + (-1.79074*pow(10,-3))*mWW +
                                                                          // (-1.51523*pow(10,-6))*pow(mWW,2) + (4.98491*pow(10,-9))*pow(mWW,3) ) ); }
  }
  // M(WW) COMPARISON_RATIO
  mWW_side_log2->SetLineColor(2);
  mWW_side_log2->Scale( 1/mWW_side_log2->Integral() );
  mWW_peak_log2->Scale( 1/mWW_peak_log2->Integral() );
  mWW_side_log2->Draw("HISTO");
  mWW_peak_log2->Draw("HISTOsame");
  leg_CORR->Draw("same");
  myc1->SaveAs("mWW_CORR.eps");
  
  //RATIO CORRECTED
  TH1F* peak_side_log2 = new TH1F( *mWW_peak_log2 );
  peak_side_log2->Sumw2();
  peak_side_log2->Divide( mWW_side_log2 );
  peak_side_log2->Draw();
  myc1->SaveAs("peak_side_log2.eps");


  float mWW_min = 240.;
  float mWW_max = 800.;
  float binWidth = 20.;
  int nBins = (int)(mWW_max-mWW_min)/binWidth;


  char sidebandsCut[300];
  sprintf( sidebandsCut, "mWW>%f && ((mJJ>40. && mJJ<60.)||(mJJ>100.&&mJJ<160.)) && leptType==%d", mWW_min, leptType);
  char signalCut[300];
  sprintf( signalCut, "mWW>%f && mJJ>60. && mJJ<100. && leptType==%d", mWW_min, leptType);


  std::cout << "Correcting MC: " << std::endl;
  TTree* tree_signalMC = chMC.CopyTree(signalCut);
  TTree* tree_sidebandsMC = chMC.CopyTree(sidebandsCut);
  TTree* tree_sidebandsMC_alpha = correctTreeWithAlpha( tree_sidebandsMC , myfit, "sidebandsMC_alpha");


TFile* fileprova = TFile::Open("PROVA.root", "recreate");
tree_sidebandsMC_alpha->Write();
fileprova->Close();


  RooRealVar* eventWeight = new RooRealVar("eventWeight", "event weight", 0., 2., "");
  RooRealVar* eventWeight_alpha = new RooRealVar("eventWeight_alpha", "alpha corrected weight", 0., 2., "");
  RooRealVar* mWW = new RooRealVar("mWW", "m_{lvjj}", mWW_min, mWW_max, "GeV");
  RooRealVar* mJJ = new RooRealVar("mJJ", "mJJ", 60., 130., "GeV");

  //RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",tree_sidebandsMC_alpha,RooArgSet(*mWW,*mJJ),"");
  RooDataSet sidebandsMC_alpha("sidebandsMC_alpha","sidebandsMC_alpha",RooArgSet(*eventWeight_alpha,*mWW,*mJJ),Import(*tree_sidebandsMC_alpha),WeightVar("eventWeight_alpha"));
  //RooDataSet sidebandsMC("sidebandsMC","sidebandsMC",tree_sidebandsMC,RooArgSet(*eventWeight,*mWW,*mJJ),"","eventWeight");
  RooDataSet signalMC("signalMC","signalMC",RooArgSet(*eventWeight,*mWW,*mJJ),Import(*tree_signalMC),WeightVar("eventWeight"));

  RooPlot *plot_sidebandsMC = mWW->frame();
  sidebandsMC_alpha.plotOn(plot_sidebandsMC, Binning(nBins));
  signalMC.plotOn(plot_sidebandsMC, Binning(nBins), LineColor(kRed));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  plot_sidebandsMC->Draw();
  char canvasName[300];
  sprintf(canvasName, "ClosureMC_%s.eps", leptType_str.c_str());
  c1->SaveAs(canvasName);


  std::cout << "Correcting DATA: " << std::endl;
  TChain chDATA("Tree_FITUL");
  chDATA.Add("HWWlvjj_DATA_6july_helicity_ALL.root");
  TTree* tree_signalDATA = chDATA.CopyTree(signalCut);
  tree_signalDATA->SetName("tree_signalDATA");
  TTree* tree_sidebandsDATA = chDATA.CopyTree(sidebandsCut);
  tree_sidebandsDATA->SetName("tree_sidebandsDATA");

  std::cout << "+++                          " << tree_signalDATA->GetEntries() << std::endl;
  std::cout << "+++                          " << tree_sidebandsDATA->GetEntries() << std::endl;

  //TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( tree_sidebandsDATA , myfit, "tree_sidebandsDATA_alpha");

//std::cout << tree_sidebandsDATA_alpha->GetEntries() << std::endl;

//RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*mWW,*mJJ),"");
////RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",tree_sidebandsDATA,RooArgSet(*eventWeight,*mWW,*mJJ),"","eventWeight");
//RooDataSet signalDATA("signalDATA","signalDATA",tree_signalDATA,RooArgSet(*eventWeight,*mWW,*mJJ),"","eventWeight");


  TH1D* h1_signalDATA = new TH1D("signalDATA", "", 30, 240., 840.);
  h1_signalDATA->Sumw2();
  TH1D* h1_sidebandsDATA_alpha = new TH1D("sidebandsDATA_alpha", "", 30, 240., 840.);
  h1_sidebandsDATA_alpha->Sumw2();
  tree_signalDATA->Project("signalDATA", "mWW", "eventWeight");
std::cout << "AAJAJAJAJA" << std::endl;
  //tree_sidebandsDATA_alpha->Project("sidebandsDATA_alpha", "mWW", "eventWeight_alpha");
  tree_sidebandsDATA->Project("sidebandsDATA_alpha", "mWW", "eventWeight");
std::cout << "h1_signalDATA->GetEntries(): " << h1_signalDATA->GetEntries() << std::endl;
std::cout << "h1_sidebandsDATA_alpha->GetEntries(): " << h1_sidebandsDATA_alpha->GetEntries() << std::endl;

  float integral_ratio = h1_signalDATA->Integral()/h1_sidebandsDATA_alpha->Integral();
  h1_sidebandsDATA_alpha->Scale(integral_ratio);


  //RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",tree_sidebandsDATA_alpha,RooArgSet(*mWW,*mJJ),"");
  //RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",RooArgSet(*eventWeight_alpha,*mWW,*mJJ),Import(*tree_sidebandsDATA_alpha),WeightVar("eventWeight_alpha"));

  RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",RooArgSet(*eventWeight_alpha,*mWW,*mJJ),Import(*tree_sidebandsDATA),WeightVar("eventWeight"));
  //RooDataSet sidebandsDATA_alpha("sidebandsDATA_alpha","sidebandsDATA_alpha",RooArgSet(*eventWeight_alpha,*mWW,*mJJ),Import(*tree_sidebandsDATA_alpha),WeightVar("eventWeight_alpha"));
  //RooDataSet sidebandsDATA("sidebandsDATA","sidebandsDATA",tree_sidebandsDATA,RooArgSet(*eventWeight,*mWW,*mJJ),"","eventWeight");
  RooDataSet signalDATA("signalDATA","signalDATA",RooArgSet(*eventWeight,*mWW,*mJJ),Import(*tree_signalDATA),WeightVar("eventWeight"));

  RooPlot *plot_sidebandsDATA = mWW->frame();
  sidebandsDATA_alpha.plotOn(plot_sidebandsDATA, Binning(nBins));
  signalDATA.plotOn(plot_sidebandsDATA, Binning(nBins), LineColor(kRed));


  // -------------------- exponential ---------------------------
  RooRealVar a_exp("a_exp","a_exp",-0.001, -2., 1.);
  RooExponential exp("exp","exp",*mWW,a_exp);
  RooFitResult *r_sidebandsMC = exp.fitTo(sidebandsDATA_alpha,SumW2Error(kTRUE));
  TH1D* h1_expSlope_DATA = new TH1D("expSlope_DATA", "", 1, 0., 1.);
  h1_expSlope_DATA->SetBinContent(1, a_exp.getVal());
  h1_expSlope_DATA->SetBinError(1, a_exp.getError());

 
  TF1* f1_expDATA = new TF1("expDATA", "[0]*exp([1]*x)");
  f1_expDATA->SetRange(240., 840.);
  f1_expDATA->FixParameter(1, a_exp.getVal());
  h1_sidebandsDATA_alpha->Fit(f1_expDATA, "R");

  c1->Clear();
  //plot_sidebandsDATA->Draw();
  h1_sidebandsDATA_alpha->Draw();
  h1_signalDATA->SetLineColor(kRed);
  h1_signalDATA->Draw("same");
  sprintf(canvasName, "ClosureDATA_%s.eps", leptType_str.c_str());
  c1->SaveAs(canvasName);
  c1->Clear();
  c1->SetLogy();
  h1_sidebandsDATA_alpha->Draw();
  h1_signalDATA->SetLineColor(kRed);
  h1_signalDATA->Draw("same");
  sprintf(canvasName, "ClosureDATA_%s_log.eps", leptType_str.c_str());
  c1->SaveAs(canvasName);
  

  char outfilename[300];
  sprintf(outfilename, "BackgroundFile_DATA_%s.root", leptType_str.c_str());
  TFile* outfileDATA = TFile::Open(outfilename, "recreate");
  outfileDATA->cd();
  h1_sidebandsDATA_alpha->Write();
  h1_signalDATA->Write();
  h1_expSlope_DATA->Write();
  outfileDATA->Close();



/*
  TH1D* h1_mWW_sidebandsMC_alpha = new TH1D("mWW_sidebandsMC_alpha", "", 30, 240., 840.);
  tree_sidebandsMC_alpha->Project("mWW_sidebandsMC_alpha", "mWW", "eventWeight*alpha");
  TH1D* h1_mWW_signalMC = new TH1D("mWW_signalMC", "", 30, 240., 840.);
  tree_signalMC->Project("mWW_signalMC", "mWW", "eventWeight");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  h1_mWW_signalMC->Draw();
  h1_mWW_sidebandsMC_alpha->SetLineColor(kRed);
  h1_mWW_sidebandsMC_alpha->Draw("same");
  
  TLegend* legendMC = new TLegend(0.55, 0.6, 0.88, 0.88);
  legendMC->SetFillColor(0);
  legendMC->SetTextSize(0.035);
  legendMC->AddEntry(h1_mWW_signalMC, "Signal Region MC", "F");
  legendMC->AddEntry(h1_mWW_sidebandsMC_alpha, "Rescaled Sidebands MC", "F");
  legendMC->Draw("same");

  c1->SaveAs("ClosureMC.eps");

  c1->Clear();


  TChain chDATA("Tree_FITUL");
  chDATA.Add("HWWlvjj_DATA_6july_helicity_ALL.root");
 
  std::cout << "Correcting signal (DATA): " << std::endl;
  TTree* tree_signalDATA = chDATA.CopyTree("mJJ>60. && mJJ<100.");
  TTree* tree_sidebandsDATA = chDATA.CopyTree("(mJJ>40. && mJJ<60.)||(mJJ>100.&&mJJ<160.)");
  TTree* tree_sidebandsDATA_alpha = correctTreeWithAlpha( tree_sidebandsDATA, myfit, "sidebandsDATA_alpha");

  TH1D* h1_mWW_sidebandsDATA_alpha = new TH1D("mWW_sidebandsDATA_alpha", "", 30, 240., 840.);
  tree_sidebandsDATA_alpha->Project("mWW_sidebandsDATA_alpha", "mWW", "eventWeight_alpha*(mJJ>0.)");
  TH1D* h1_mWW_signalDATA = new TH1D("mWW_signalDATA", "", 30, 240., 840.);
  tree_signalDATA->Project("mWW_signalDATA", "mWW", "eventWeight");

  h1_mWW_signalDATA->Draw();
  h1_mWW_sidebandsDATA_alpha->SetLineColor(kRed);
  h1_mWW_sidebandsDATA_alpha->Draw("same");
  
  TLegend* legendDATA = new TLegend(0.55, 0.6, 0.88, 0.88);
  legendDATA->SetFillColor(0);
  legendDATA->SetTextSize(0.035);
  legendDATA->AddEntry(h1_mWW_signalDATA, "Signal Region DATA", "F");
  legendDATA->AddEntry(h1_mWW_sidebandsDATA_alpha, "Rescaled Sidebands DATA", "F");
  legendDATA->Draw("same");

  c1->SaveAs("ClosureMC.eps");
*/

} //End


void getBins_int( int nBins_total, Double_t* Lower, Double_t xmin, Double_t xmax, bool plotLog) {

  Double_t Lower_exact;
  int nBins = nBins_total-1;
  const double dx = (plotLog) ? pow((xmax / xmin), (1. / (double)nBins)) : ((xmax - xmin) / (double)nBins);
  Lower[0] = xmin;
  Lower_exact = Lower[0];
  for (int i = 1; i != nBins; ++i) {

    if (plotLog) {
      Lower_exact *= dx;
      Lower[i] = TMath::Ceil(Lower_exact);
    } else {
      Lower[i] = TMath::Ceil(Lower[i-1] + dx);
    }

  }

  Lower[nBins] = xmax;

}


TTree* correctTreeWithAlpha( TTree* tree, TF1* f1_alpha, const std::string& name ) {

  Float_t mJJ;
  tree->SetBranchAddress( "mJJ", &mJJ );
  Float_t mWW;
  tree->SetBranchAddress( "mWW", &mWW );
  Float_t eventWeight;
  tree->SetBranchAddress( "eventWeight", &eventWeight );

  //TTree* newTree = tree->CloneTree();
  TTree* newTree = new TTree(name.c_str(), "");
//  newTree->SetName(name.c_str());

  Float_t newWeight;
  newTree->Branch( "eventWeight_alpha", &newWeight, "newWeight/F" );
  newTree->Branch( "mJJ", &mJJ, "mJJ/F" );
  newTree->Branch( "mWW", &mWW, "mWW/F" );
  newTree->Branch( "eventWeight", &eventWeight, "eventWeight/F" );

  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry( iEntry );
    if( (iEntry % 10000)==0 ) std::cout << "Entry: " << iEntry << "/" << nentries << std::endl;

    // alpha correction
    newWeight = eventWeight*f1_alpha->Eval(mWW);
    //std::cout<< eventWeight <<"  "<<newWeight<<std::endl;
    newTree->Fill();

  }

  return newTree;

}
