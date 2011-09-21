#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;
using namespace RooFit;

RooFitResult* fitBkgZZinvMass_CB_DATA(double rot=0.){

  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".L PDFs/RooFermi_cc.so");
  gROOT->ProcessLine(".L PDFs/RooCB_cc.so");
  gSystem->Load("libFFTW");
  //gROOT->ProcessLine(".L ~/tdrstyle.C");
  //setTDRStyle();


  // --------------------- initial values -----------------
  string plotName;
  string cutString = "((mJJ>40 && mJJ<60.) || (mJJ>100. && mJJ<160.))";
  vector<double> paramVal;
  paramVal.push_back(160.); //cutOff   
  paramVal.push_back(5.257);	//beta	    
  paramVal.push_back(200.);	//CB_mean  
  paramVal.push_back(53.876);	//CB_wdth   54.686 --- 1.0 /fb fits
  paramVal.push_back(25.775);	//CB_n	    13.39 
  paramVal.push_back(-.8312);	//CB_alpha  -.879 
  cout << "parameters set..."<< endl;
  // --------------------- measurable (ZZ invariant mass) ----------------

  RooRealVar mWW("mWW", "WW inv mass", 140.,800.);
  RooRealVar isSB("isSB","isSB",-1.,1.);
  RooRealVar wght("wght","wght",0,100.);
  RooRealVar mJJ("mJJ","mJJ",0,150.);

  // ==================== defining bkg PDF ==========================
  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",paramVal.at(0),0,1000);
  cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",paramVal.at(1),0,50);
  beta.setConstant(kTRUE);
	     		       
  RooFermi fermi("fermi","fermi function",mWW,cutOff,beta);
  // -------------------- double gauss ---------------------------

  double a=cos(rot)*paramVal.at(5) + sin(rot)*paramVal.at(3);
  double w=-sin(rot)*paramVal.at(5) + cos(rot)*paramVal.at(3);
  cout << "alpha true: " << a << endl;
  cout << "width true: " << w << endl;
  
  RooRealVar m("m","m",paramVal.at(2),200,1000);
  m.setConstant(kTRUE);
  RooRealVar wdth("wdth","#beta",w,-1000,1000);
  //wdth.setConstant(kTRUE);
  RooRealVar n("n","n",paramVal.at(4),0,100);
  n.setConstant(kTRUE);
  RooRealVar alpha("alpha","#alpha",a,-1000,1000); 
  //alpha.setConstant(kTRUE);
  RooRealVar theta("theta","theta",rot,-3.1415,3.1415); 
  theta.setConstant(kTRUE);
  
  RooCB CB("CB","Crystal ball",mWW,m,wdth,alpha,n,theta);
  
  RooProdPdf background("background","background",RooArgSet(fermi,CB));

  // ------------------ get data --------------------------
  RooDataSet data_bkg;  
  RooDataSet *data_temp;
  TFile *file;
 
  // for reading sideband extrapolated data...
  file = new TFile("NEW_HWWlvjj_DATA_6july_helicity_ALL.root");
  data_bkg=new RooDataSet("data_bkg","data_bkg",(TTree*)file->Get("selectedEvents"),RooArgSet(mWW,isSB,mJJ,wght),cutString.c_str(),"wght");

  cout << "check" << endl;

  // --------- draw MC data -------------------

  //RooFitResult *r = background.fitTo(data_bkg,SumW2Error(kTRUE),InitialHesse(kTRUE),Minos(kTRUE),Save());  
  RooFitResult *r = background.fitTo(data_bkg,SumW2Error(kTRUE),Save());  

  RooPlot *plot_MCbkg = mWW.frame(120,800,68);
  data_bkg.plotOn(plot_MCbkg,DataError(RooAbsData::SumW2));
  background.plotOn(plot_MCbkg,VisualizeError(*r,2.0,kFALSE),FillColor(kYellow));
  background.plotOn(plot_MCbkg,VisualizeError(*r,1.0,kFALSE),FillColor(kGreen));  
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  background.plotOn(plot_MCbkg);
  data_bkg.plotOn(plot_MCbkg,DataError(RooAbsData::SumW2));
  
  plot_MCbkg->Draw();
  
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  RooRealVar* alpha_0 = (RooRealVar*)(r->floatParsFinal().find("wdth"));
  RooRealVar* beta_0 = (RooRealVar*) (r->floatParsFinal().find("alpha"));

  Double_t x1= alpha_0->getVal();
  Double_t x2= beta_0->getVal();
  Double_t s1= alpha_0->getError();// magic number goes here
  Double_t s2= beta_0->getError();// magic number goes here
  Double_t rho= r->correlation("alpha", "wdth");
  r->Print("V");
  r->correlationMatrix().Print("v");
  r->covarianceMatrix().Print("v");

//RooEllipse *contour= new RooEllipse("contour",x1,x2,s1,s2,rho,500);
//contour->SetLineWidth(2) ;
//

//RooPlot *plot = new RooPlot(*alpha_0,*beta_0,wdth.getVal()-2*wdth.getError(),wdth.getVal()+2*wdth.getError(),
//			      alpha.getVal()-2*alpha.getError(),alpha.getVal()+2*alpha.getError());
////RooPlot *plot = new RooPlot(*alpha_0,*beta_0,40,100,-1.5,.4);
//plot->addPlotable(contour);

//r->plotOn(plot,*alpha_0,*beta_0,"ME12");

//plot->Draw();

  return r;
}

