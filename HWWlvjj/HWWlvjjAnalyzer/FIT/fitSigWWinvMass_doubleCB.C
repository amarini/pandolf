#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

using namespace RooFit;

vector<double> readParamsFromFile(string fileName){

  vector<double> param;

  ifstream inFile(fileName.c_str());
  if(inFile.is_open()){
    string burnIn;
    string paramName;
    double paramVal;
    while(!inFile.eof()){
      inFile >> paramName >> paramVal >> burnIn >> burnIn;
      param.push_back(paramVal);
    }

  }
  return param;
}

vector<double> fitSigWWinvMass_doubleCB(int iSig=2, int fixVar=0) {
  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".L PDFs/RooRelBW_cc.so");
  gROOT->ProcessLine(".L PDFs/RooDoubleCB_cc.so");
  gROOT->ProcessLine(".L PDFs/RooFermi_cc.so");
  gSystem->Load("libFFTW");
  //gROOT->ProcessLine(".L ~/tdrstyle.C");
  //setTDRStyle();
  // --------------------- initial values -----------------

  double meanIV[7]={200.,300.,350.,400.,450.,500.,550.};
  double sigmaIV[7]={1.43,8.445,15.39,29.42,47.01,67.52,93.15};
  double effWidth[7]={10.1,13.12,18.36,31.05,48.1,68.26,93.7};
  double fitRangeLow= meanIV[iSig]-10*effWidth[iSig]<150?150:meanIV[iSig]-10*effWidth[iSig];
  double fitRangeHigh=meanIV[iSig]+10*effWidth[iSig]>800?800:meanIV[iSig]+10*effWidth[iSig];
  string plotName;
  string massName[7]={"200","300","350","400","450","500","550"};
  string cutString = "mJJ>60. && mJJ<100.";

  float mass = meanIV[iSig];
  std::cout << "+++ MASS: " << mass << std::endl;
      
  // --------------------- measurable (WW invariant mass) ----------------
  RooRealVar mWW("mWW", "ww inv mass",fitRangeLow,fitRangeHigh);
  RooRealVar eventWeight("eventWeight","eventWeight",0,10);
  RooRealVar mJJ("mJJ","mJJ",0,1000);

  cout << "declared tree vars" << endl;
  // ====================== defining signal PDF =========================

  //float CB_mean_0 = 70.6146-.697703*mass+0.00212559*mass*mass-0.00000180624*mass*mass*mass;
  float CB_mean_0 = -8.52531+5.88265e-03*mass+1.10855e-04*mass*mass;
  float CB_sigma_0 = -1.63923+8.50605e-02*mass-1.32823e-04*mass*mass; //fp
  //float CB_alpha1_0 = 1.0;
  float CB_alpha1_0 = 0.85;
  float CB_n1_0 = 3.38183-0.00421732*mass;
  float CB_alpha2_0 = 2.82669-1.21212e-02*mass+1.84626e-05*mass*mass; //fp
  float CB_n2_0 = -1.37066+0.0190719*mass-0.0000250673*mass*mass;


//for(int i=0; i<param.size(); i++){
//  cout << "param[" << i << "]: " << param.at(i) << endl;
//}

  // -------------------- fermi ------------------------
  
  RooRealVar cutOff("cutOff","cutOff",190.-32.5+65.*meanIV[iSig]/400.,0.,1000.); //param[6],0,1000);
  cutOff.setConstant(kTRUE);
  RooRealVar g("g","g",5.-12.5+25.*meanIV[iSig]/400.,0.,100.); //param[7],0,100);
  g.setConstant(kTRUE);

  RooFermi fermi("fermi","fermi",mWW,cutOff,g);

  // ------------------- fermi for high mass cutoff --------------

  RooRealVar cutOff2("cutOff2","cutOff2",700.,0.,1000.); //param[6],0,1000);
  cutOff2.setConstant(kTRUE);
  RooRealVar g2("g2","g2",-70.,-100.,0.); //param[7],0,100);
  g2.setConstant(kTRUE);

  RooFermi fermi2("fermi2","fermi2",mWW,cutOff2,g2);

  // ------------------- Relativistic BW --------------------------------
  //
 
  RooRealVar BW_mean("BW_mean", "mean",meanIV[iSig],0.,1000.);
  BW_mean.setConstant(kTRUE);
  RooRealVar BW_sigma("BW_sigma", "sigma",sigmaIV[iSig],0,200);
  BW_sigma.setConstant(kTRUE);
  RooRealVar BW_n("BW_n","n",0.,0.,1.);
  if(iSig>3 && iSig<10) BW_n.setVal(0.0); 
  else BW_n.setVal(1.0);
  BW_n.setConstant(kTRUE);

  cout << "BW_n: " << BW_n.getVal() << endl;

  RooRelBW BW("BW","Relativistic B-W",mWW,BW_mean,BW_sigma,BW_n);

  // ------------------- Crystal Ball -------------------------------
  RooRealVar CB_mean("CB_mean","param 1 of CB",CB_mean_0,0.,100.);
//CB_mean.setConstant(kTRUE);
  RooRealVar CB_sigma("CB_sigma","param 2 of CB",CB_sigma_0,0.,100.);
//CB_sigma.setConstant(kTRUE);
  RooRealVar CB_alpha1("CB_alpha1","param 3 of CB",CB_alpha1_0,0.,100.);
//CB_alpha1.setConstant(kTRUE);
  RooRealVar CB_n1("CB_n1","param 4 of CB",CB_n1_0,0.,100.);
//CB_n1.setConstant(kTRUE);
  RooRealVar CB_alpha2("CB_alpha2","param 3 of CB",CB_alpha2_0,0.,100.);
//CB_alpha2.setConstant(kTRUE);
  RooRealVar CB_n2("CB_n2","param 4 of CB",CB_n2_0,0.,100.);
//CB_n2.setConstant(kTRUE);

//if(fixVar==1) CB_mean.setConstant(kFALSE);
//if(fixVar==2) CB_sigma.setConstant(kFALSE);
//if(fixVar==3) CB_alpha1.setConstant(kFALSE);
//if(fixVar==4) CB_n1.setConstant(kFALSE);
//if(fixVar==5) CB_alpha2.setConstant(kFALSE);
//if(fixVar==6) CB_n2.setConstant(kFALSE);

  RooDoubleCB CB("CB","Crystal Ball",mWW,CB_mean,CB_sigma,CB_alpha1,CB_n1,CB_alpha2,CB_n2);
  
  //------------------------ convolution -------------------------
  //Set #bins to be used for FFT sampling to ...
  mWW.setBins(10000,"fft");

  //RooFFTConvPdf signal("signal","Rel-BW (X) CB",mWW,BW,CB);
  //signal.setBufferFraction(1.0);
  RooFFTConvPdf sig("sig","Rel-BW (X) CB",mWW,BW,CB);
  sig.setBufferFraction(1.0);
  
  RooProdPdf signal("signal","signal",RooArgSet(sig,fermi,fermi2));
  //RooProdPdf signal("signal","signal",RooArgSet(sig));
  //RooProdPdf signal("signal","signal",RooArgSet(CB));

  cout << "PDFs defined" << endl;

  // ---------------------------- get data --------------------------

  char signalFileName[800];
  sprintf( signalFileName, "HWWlvjj_GluGluToHToWWToLNuQQ_M-%.0f_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1_helicity_ALL.root", mass);
  
  TFile* f = new TFile(signalFileName);

  cout << "file loaded" << endl;
  RooDataSet data_sig("data_sig","data",(TTree*)f->Get("Tree_FITUL"),RooArgSet(mWW,mJJ,eventWeight),cutString.c_str(),"eventWeight");

  cout << "dataset loaded" << endl;
  // ---------------------- fit signal and plot -------------------
  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  
  //if(fixVar!=0)
    signal.fitTo(data_sig,SumW2Error(kTRUE));

  RooPlot* plot_sig = mWW.frame();
  data_sig.plotOn(plot_sig,DataError(RooAbsData::SumW2),Binning(100));
  signal.plotOn(plot_sig);
  plot_sig->Draw();
  
  plotName="Parameterizations_SIG_"+massName[iSig]+"GeVHiggs.eps";
  c2->SaveAs(plotName.c_str());

  string outFileName="HWWlvjj_signal_"+massName[iSig]+"GeV_parameters.txt";

  ofstream output(outFileName.c_str());
  output << BW_mean.GetName() << " " << BW_mean.getVal() << " " << BW_mean.getError() << endl;
  output << BW_sigma.GetName() << " " << BW_sigma.getVal() << " " << BW_sigma.getError() << endl;
  output << BW_n.GetName() << " " << BW_n.getVal() << " " << BW_n.getError() << endl;
  output << CB_mean.GetName() << " " << CB_mean.getVal() << " " << CB_mean.getError() << endl;
  output << CB_sigma.GetName() << " " << CB_sigma.getVal() << " " << CB_sigma.getError() << endl;
  output << CB_alpha1.GetName() << " " << CB_alpha1.getVal() << " " << CB_alpha1.getError() << endl;
  output << CB_n1.GetName() << " " << CB_n1.getVal() << " " << CB_n1.getError() << endl;
  output << CB_alpha2.GetName() << " " << CB_alpha2.getVal() << " " << CB_alpha2.getError() << endl;
  output << CB_n2.GetName() << " " << CB_n2.getVal() << " " << CB_n2.getError() << endl;
  output << cutOff.GetName() << " " << cutOff.getVal() << " " << cutOff.getError() << endl;
  output << g.GetName() << " " << g.getVal() << " " << g.getError() << endl;
  output << cutOff2.GetName() << " " << cutOff2.getVal() << " " << cutOff2.getError() << endl;
  output << g2.GetName() << " " << g2.getVal() << " " << g2.getError() << endl;
  
  output.close();

  vector<double> result;
  if(fixVar==1){
    result.push_back(CB_mean.getVal()  );
    result.push_back(CB_mean.getError());
  }
  if(fixVar==2){
    result.push_back(CB_sigma.getVal()  );
    result.push_back(CB_sigma.getError());
  }
  if(fixVar==3){
    result.push_back(CB_alpha1.getVal()  );
    result.push_back(CB_alpha1.getError());
  }
  if(fixVar==4){
    result.push_back(CB_n1.getVal()  );
    result.push_back(CB_n1.getError());
  }
  if(fixVar==5){
    result.push_back(CB_alpha2.getVal()  );
    result.push_back(CB_alpha2.getError());
  }
  if(fixVar==6){
    result.push_back(CB_n2.getVal()  );
    result.push_back(CB_n2.getError());
  }
  return result;
}

