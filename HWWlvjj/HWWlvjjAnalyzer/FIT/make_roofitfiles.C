#include "TH1F.h"
#include "TCanvas.h"
#include <Riostream.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


#include "RooRealVar.h"
#include "RooGenericPdf.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
#include "RooExponential.h"

#include "PDFs/RooRodenbach.h"
#include "PDFs/RooFermi.h"
#include "PDFs/RooDoubleCB.h"
#include "PDFs/RooCB.h"
#include "PDFs/RooRelBW.h"

#include "TH1F.h"
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"

using namespace std;
//using namespace ROOT::Math;
using namespace RooFit;


//inputs: btag category, observed bkg yield (-> expected one for MC limit calc)
//        mass of Higgs, sigma of Higgs, 4 parameters of CB (depend on mass)

void make_roofitfiles(int chan, double massH, double sigmaH, double &obs_yield, double &exp_yield, vector<double> cb_pars){

  gSystem->Load("libRooFit");
  cout<<"Trying to load custom PDFS ..."<<flush<<endl;
  gSystem->Load("libFFTW");


  string str_chan="dummychan";
  if(chan==0)str_chan="ele";
  else if(chan==1)str_chan="mu";
  else cout<<"Unrecognized number of channels: "<<chan<<endl;


  //integration window
  double effWidth=sqrt(sigmaH*sigmaH+100);  // effective width used for defining window
  double fitRangeLow=99.;
  double fitRangeHigh=101.0; 

  ///expo functions/////////////////////////////////////////////////////
  /*if(massH-10*effWidth<230.0)  fitRangeLow=230.0;
  else fitRangeLow=massH-10*effWidth;

  if(massH+10*effWidth>800.0)  fitRangeHigh=800.;
  else fitRangeHigh=massH+10*effWidth;
  cout<<"----- FIT RANGE : "<<fitRangeLow<<" - "<< fitRangeHigh<<endl;*/
  /////////////////////////////////////////////////////////////////////

  //////////////////////////eps functions /////////////////////////////
  if(massH-10*effWidth<183.0)  fitRangeLow=183.0;
  else fitRangeLow=massH-10*effWidth;

  if(massH+10*effWidth>800.0)  fitRangeHigh=800.;
  else fitRangeHigh=massH+10*effWidth;

  cout<<"----- FIT RANGE : "<<fitRangeLow<<" - "<< fitRangeHigh<<endl;

  ////////////////////////////////////////////////////////////////

  std::ostringstream ossm;
  ossm<<massH; 
  string str_massH=ossm.str() ;

  // --------------------- initial values -----------------
  RooRealVar CMS_hwwlvqq_mWW("CMS_hwwlvqq_mWW", "WW inv mass",fitRangeLow,fitRangeHigh);

  // ==================== defining bkg PDF ==========================


  //exponential ////////////////////////////////////////////////////////////
  //par0 will be the overall bkgd normalization, used in the datacard
  /*string bkgp1name="CMS_hww2l2q_bkg"+str_btag+"p1";
 
  double slope_Val[3]={-0.0150692, -0.0124851, -0.0127107};

  RooRealVar slope(bkgp1name.c_str(),"slope",slope_Val[btag],-1.0,0.0);
  slope.setConstant(kTRUE);

  RooExponential background("background","Exponential background",CMS_hwwlvqq_mWW,slope);
  */

  //new shapes from Andrew///////////////////////
  /*vector<double> BKGparam;
  if(btag==0){
    BKGparam.push_back(186.405);   //fermi: cutOff
    BKGparam.push_back(5.25694);   //fermi: beta
    BKGparam.push_back(226.78);   //m
    BKGparam.push_back(44.1607);   //width
    BKGparam.push_back(254.615);   //alpha
  }else if(btag==1){
    BKGparam.push_back(184.04);   //fermi: cutOff
    BKGparam.push_back(3.225);   //fermi: beta
    BKGparam.push_back(200);   //m
    BKGparam.push_back(59.9061);   //width
    BKGparam.push_back(244.851);   //alpha
  }else if(btag==2){
    BKGparam.push_back(182.32);   //fermi: cutOff
    BKGparam.push_back(2.8492);   //fermi: beta
    BKGparam.push_back(200.002);   //m
    BKGparam.push_back(103.69);   //width
    BKGparam.push_back(646.467);   //alpha
  }
  // ------------------------ fermi ------------------------------
  RooRealVar cutOff("cutOff","position of fermi",BKGparam.at(0),0,1000);
  cutOff.setConstant(kTRUE);
  RooRealVar beta("beta","width of fermi",BKGparam.at(1),0,50);
  beta.setConstant(kTRUE);
	     		       
  RooFermi fermi("fermi","fermi function",CMS_hwwlvqq_mWW,cutOff,beta);
  // -------------------- double gauss ---------------------------
  string bkgp1name="CMS_hww2l2q_bkg"+str_btag+"p1";
  string bkgp2name="CMS_hww2l2q_bkg"+str_btag+"p2";
  string bkgp3name="CMS_hww2l2q_bkg"+str_btag+"p3";
  RooRealVar m(bkgp1name.c_str(),bkgp1name.c_str(),BKGparam.at(2),200,1000);
  m.setConstant(kTRUE);
  RooRealVar wdth(bkgp2name.c_str(),bkgp2name.c_str(),BKGparam.at(3),0,1000);
  wdth.setConstant(kTRUE);
  RooRealVar alpha(bkgp3name.c_str(),bkgp3name.c_str(),BKGparam.at(4),200,1000); 
  alpha.setConstant(kTRUE);

  RooRodenbach Rod("Rod","Rod",CMS_hwwlvqq_mWW,m,wdth,alpha);

  RooProdPdf background("background","background",RooArgSet(fermi,Rod));
  */
  //CB for EPS/////////////////////////////////////////////////////////////

  // ------------------------ fermi ------------------------------            
  vector<double> BKGparam;
     BKGparam.push_back(186.41); //cutoff
     BKGparam.push_back(5.257); //beta
     BKGparam.push_back(222.72); //mean
     BKGparam.push_back(-.116428); //width
     BKGparam.push_back(13.39); //n
     BKGparam.push_back(54.704); //alpha
     BKGparam.push_back(1.589); //theta

  //LP ones
  /*
  if(btag==0){
    BKGparam.push_back(186.41);  //cutOff 
    BKGparam.push_back(5.257);   //beta  
    BKGparam.push_back(222.72);  //CB_mean   
    BKGparam.push_back(53.876);  //CB_wdth  ---- uncorrelated  
    BKGparam.push_back(25.775);   //CB_n              
    BKGparam.push_back(.11276);   //CB_alpha ---- uncorrelated          
    BKGparam.push_back(.01752);   //theta                       
    }
  if(btag==1){
    BKGparam.push_back(184.04);  //cutOff    
    BKGparam.push_back(3.225);   //beta   
    BKGparam.push_back(166.6);   //CB_mean                                       
    BKGparam.push_back(94.7425);  //CB_wdth   --- uncorrelated     
    BKGparam.push_back(6.5670);  //CB_n                                                                                           
    BKGparam.push_back(.26330 );  //CB_alpha   --- uncorrelated  
    BKGparam.push_back(.0183  );  //theta                                         
    }
  if(btag==2){
    BKGparam.push_back(182.32);  //cutOff   
    BKGparam.push_back(2.8492);  //beta      
    BKGparam.push_back(226.23);  //CB_mean     
    BKGparam.push_back(39.41);    //CB_wdth   --- uncorrelated  
    BKGparam.push_back(6.1580);  //CB_n   
    BKGparam.push_back(-.2518);   //CB_alpha   --- uncorrelated   
    BKGparam.push_back(.0125);    //theta                                        
    }
  */

  // -------------------- fermi ---------------------------
  RooRealVar cutOff_BKG("cutOff_BKG","position of fermi",BKGparam.at(0),0,1000);
  cutOff_BKG.setConstant(kTRUE);
  RooRealVar beta_BKG("beta_BKG","width of fermi",BKGparam.at(1),0,50);
  beta_BKG.setConstant(kTRUE);
	     		       
  RooFermi fermi_BKG("fermi_BKG","fermi function",CMS_hwwlvqq_mWW,cutOff_BKG,beta_BKG);
 // -------------------- double gauss ---------------------------
  //par0 will be the overall bkgd normalization, used in the datacard
  string bkgp1name="CMS_hwwlvqq_bkg_p1"; //m
  string bkgp2name="CMS_hwwlvqq_bkg_p2"; //width
  string bkgp3name="CMS_hwwlvqq_bkg_p3"; //n
  string bkgp4name="CMS_hwwlvqq_bkg_p4"; //alpha
  string bkgp5name="CMS_hwwlvqq_bkg_p5"; //theta (rotation)

  RooRealVar m(bkgp1name.c_str(),bkgp1name.c_str(),BKGparam.at(2),100.,1000.);
  m.setConstant(kTRUE);
  RooRealVar wdth(bkgp2name.c_str(),bkgp2name.c_str(),BKGparam.at(3),0,1000);
  wdth.setConstant(kTRUE);
  RooRealVar n(bkgp3name.c_str(),bkgp3name.c_str(),BKGparam.at(4),0.,1001.);//2.75833,0,1000);
  n.setConstant(kTRUE);
  RooRealVar alpha(bkgp4name.c_str(),bkgp4name.c_str(),BKGparam.at(5),-100,100);  //0,100);  //,-100,0);
  alpha.setConstant(kTRUE);
  RooRealVar theta(bkgp5name.c_str(),bkgp5name.c_str(),BKGparam.at(6),-3.1415,3.1415); 
  theta.setConstant(kTRUE);

  RooCB CB_BKG("CB_BKG","Crystal ball",CMS_hwwlvqq_mWW,m,wdth,alpha,n, theta);
  RooProdPdf background("background","background",RooArgSet(fermi_BKG,CB_BKG));
  
  ///////////////////////////////////////////////////////////////////////////////////

  ////Fill dataset with REAL DATA 

  RooRealVar eventWeight("eventWeight","eventWeight",0,100.);
  RooRealVar mJJ("mJJ","mJJ",0,160.);
  RooRealVar leptType("leptType","lepton type",-1,2);

  string lept_sel= chan==0 ? "leptType==0" :"leptType==1" ;//opposite convention btw Francesco and me
  string tree_sel= "mJJ>60.0 && mJJ<100.0 && "+lept_sel;
  stringstream ossmww1;
  ossmww1 << float(fitRangeLow);
  string mwwcut="CMS_hwwlvqq_mWW>"+ossmww1.str(); 
  stringstream ossmww2;
  ossmww2 << float(fitRangeHigh);
  mwwcut+="&&CMS_hwwlvqq_mWW<"+ossmww2.str();
  cout<<"$$$$$$ TEMP SEL:  "<<mwwcut.c_str()<<"  $$$$$$$$$$$$$$$$$$$$$$ "<<fitRangeLow<<" - "<< fitRangeHigh<<endl;
  tree_sel+=" && "+mwwcut;
 
 
  /* TFile *dfile = new TFile("fileout-999invpb.root");
  //RooArgList arg1(mww);
  RooFormulaVar cut1("mycut1",tree_sel.c_str(),RooArgList(mww,nBTags,mZjj,leptType));

  RooDataSet *data_b=new RooDataSet("data_bkg","data_bkg",(TTree*)dfile->Get("tree_passedEvents"),
				    RooArgSet(mww,nBTags,mZjj,leptType),cut1,"eventWeight");
  obs_yield=double(data_b->numEntries());
  //RooDataSet *data_b = background.generate(x,int(obs_yield));
  cout<<"\nBTAG "<<btag<<"   OBS_YIELDS: "<<obs_yield<<" ->   "<<int(obs_yield)<<endl<<endl;*/

  TFile* file = new TFile("./convertedTree_LP_20110811.root");
  RooFormulaVar cut1("mycut1",tree_sel.c_str(),RooArgList(CMS_hwwlvqq_mWW,mJJ,leptType));
  RooDataSet *dataset_obs_orig=new RooDataSet("dataset_obs_orig","dataset_obs_orig",(TTree*)file->Get("tree_passedEvents"),
					      RooArgSet(CMS_hwwlvqq_mWW,mJJ,leptType),
					      cut1,"eventWeight");

  obs_yield=double(dataset_obs_orig->numEntries());

  RooArgSet *newMwwargset= new RooArgSet(CMS_hwwlvqq_mWW);
  RooDataSet *dataset_obs=(RooDataSet*) dataset_obs_orig->reduce(*newMwwargset);
  dataset_obs->SetName("dataset_obs");
  cout<<"Dataset entries: ORIG "<< dataset_obs_orig->sumEntries()<< "   NEW "<<dataset_obs->sumEntries()<<endl;
  // ----------------------------------------------

  // ====================== defining signal PDF =========================

  vector<double> param;
    param.push_back(70.6146-.697703*massH+0.00212559*massH*massH-0.00000180624*massH*massH*massH);
    param.push_back(-5.967+0.05885*massH-0.00006977*massH*massH);
    param.push_back(1.0);
    param.push_back(3.38183-0.00421732*massH);
    param.push_back(1.0);
    param.push_back(-1.37066+0.0190719*massH-0.0000250673*massH*massH);

    for(int i=0; i<param.size(); i++){
    cout << "param[" << i << "]: " << param.at(i) << endl;
  }

  // -------------------- fermi ------------------------
  
  RooRealVar cutOff_SIG("cutOff_SIG","cutOff",190-32.5+65*massH/400,0,1000); 
  cutOff_SIG.setConstant(kTRUE);  
  RooRealVar g_SIG("g_SIG","g",5-12.5+25*massH/400,0,100); 
  g_SIG.setConstant(kTRUE);

  RooFermi fermi_SIG("fermi_SIG","fermi",CMS_hwwlvqq_mWW,cutOff_SIG,g_SIG);

  // ------------------- fermi for high mass cutoff --------------

  RooRealVar cutOff2_SIG("cutOff2_SIG","cutOff2",700,0,1000);
  cutOff2_SIG.setConstant(kTRUE);
  RooRealVar g2_SIG("g2_SIG","g2",-70.0,-100.0,0.0);
  g2_SIG.setConstant(kTRUE);

  RooFermi fermi2_SIG("fermi2_SIG","fermi2",CMS_hwwlvqq_mWW,cutOff2_SIG,g2_SIG);

  // ------------------- Relativistic BW --------------------------------
  //
 
  RooRealVar BW_mean("BW_mean", "mean",massH,0,1000);
  BW_mean.setConstant(kTRUE);
  RooRealVar BW_sigma("BW_sigma", "sigma",sigmaH,0,200);
  BW_sigma.setConstant(kTRUE);
  RooRealVar BW_n("BW_n","n",0.,0.,1.);
  BW_n.setConstant(kTRUE);

  RooRelBW BW("BW","Relativistic B-W",CMS_hwwlvqq_mWW,BW_mean,BW_sigma,BW_n);

  // ------------------- Crystal Ball -------------------------------
  string sigp1name="CMS_hwwlvqq_sig_p1"; //m
  string sigp2name="CMS_hwwlvqq_sig_p2"; //width
   RooRealVar CB_mean(sigp1name.c_str(),sigp1name.c_str(),param[0],0.,100.);
  CB_mean.setConstant(kTRUE);
  RooRealVar CB_sigma(sigp2name.c_str(),sigp2name.c_str(),param[1],0.,100.);
  CB_sigma.setConstant(kTRUE);
  RooRealVar CB_alpha1("CB_alpha1","param 3 of CB",param[2],0.,100.);
  CB_alpha1.setConstant(kTRUE);
  RooRealVar CB_n1("CB_n1","param 4 of CB",param[3],0.,100.);
  CB_n1.setConstant(kTRUE);
  RooRealVar CB_alpha2("CB_alpha2","param 3 of CB",param[4],0.,100.);
  CB_alpha2.setConstant(kTRUE);
  RooRealVar CB_n2("CB_n2","param 4 of CB",param[5],0.,100.);
  CB_n2.setConstant(kTRUE);

  RooDoubleCB CB_SIG("CB_SIG","Crystal Ball",CMS_hwwlvqq_mWW,CB_mean,CB_sigma,CB_alpha1,CB_n1,CB_alpha2,CB_n2);
  //------------------------ convolution -------------------------
  CMS_hwwlvqq_mWW.setBins(10000,"fft");

  RooFFTConvPdf sig("sig","Rel-BW (X) CB",CMS_hwwlvqq_mWW,BW,CB_SIG);
  sig.setBufferFraction(1.0);
  
  RooProdPdf signal("signal","signal",RooArgSet(sig,fermi_SIG,fermi2_SIG));

  //RooProdPdf signal_ggH(signal, "ggH_prodPDF");
  //RooProdPdf signal_VBF(signal, "qqH_prodPDF");
 

  //--- write everything into the workspace -------  

  RooWorkspace* w = new RooWorkspace("w","w");
  w->addClassDeclImportDir("/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
  // w->addClassDeclImportDir("/afs/cern.ch/user/w/whitbeck/scratch0/HiggsStats/newHiggsStats/CMSSW_4_1_3/src/HiggsAnalysis/CombinedLimit/data/PDFs/");
  //w->addClassDeclImportDir("/afs/cern.ch/user/b/bonato/scratch0/PhysAnalysis/CMSSW_4_2_3_patch5/src/ZJetsAnalysis/ZJetsAnalysisV1/test/statistical_tools/PDFs/");
  //w->addClassDeclImportDir("/afs/cern.ch/user/s/sbologne/scratch0/CMSSW/CMSSW_4_2_4/src/HiggsAnalysis/CombinedLimit/test/rotatedEPSForLP/PDFs/");
  w->addClassDeclImportDir("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/HWWlvjj/HWWlvjjAnalyzer/test/analysis/FIT/PDFs");

  w->importClassCode(RooFermi::Class(),kTRUE);
  w->importClassCode("RooFermi",kTRUE);
  w->importClassCode(RooRelBW::Class(),kTRUE);
  w->importClassCode("RooRelBW",kTRUE);
  w->importClassCode(RooDoubleCB::Class(),kTRUE);
  w->importClassCode("RooDoubleCB",kTRUE);
  w->importClassCode(RooCB::Class(),kTRUE);
  w->importClassCode("RooCB",kTRUE);
  //fro roorodenbach!!!!!!!!!!!!!!!!!!!!!!
  //w->importClassCode(RooRodenbach::Class(),kTRUE);
  //w->importClassCode("RooRodenbach",kTRUE);
  w->import(background);
  w->import(signal);
  // w->import(signal_ggH);
  // w->import(signal_VBF);
  w->import(*dataset_obs);

  //string outFileName="datacards_20110803_epsrotatedRange/"+str_massH+"/hww2l2q_"+str_chan+str_btag+".input.root";
  string outFileName="datacards/"+str_massH+"/hwwlvqq_"+str_chan+".input.root";
  
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");
  w->Write();
  outFile->Close();

  //calculate expected bkg events
  //eps functions ////////////////////////////
  RooRealVar CMS_hwwlvqq_mWWfull("CMS_hwwlvqq_mWWfull", "ww inv mass",183.0 ,800.0);
  //expo functions
  //RooRealVar CMS_hwwlvqq_mWWfull("CMS_hwwlvqq_mWWfull", "ww inv mass",230.0 ,800.0);

  //expo////////////////////////////////////////////////////////////
  //RooExponential backgroundFull("backgroundFull","Exponential background over Full range",CMS_hwwlvqq_mWWfull,slope);
  
  //eps ////////////////////////////////////////////////////////////////////////////////////////
  RooGenericPdf fermiFull("fermiFull","fermi function","1/(1+exp((@1-@0)/@2))",RooArgList(CMS_hwwlvqq_mWWfull,cutOff_BKG,beta_BKG));
  RooCB CBbkgFull("CBbkgFull","Crystal ball for background",CMS_hwwlvqq_mWWfull,m,wdth,alpha,n, theta);
  RooProdPdf backgroundFull("backgroundFull","backgroundFull",RooArgSet(fermiFull,CBbkgFull));

  //new shape andrew////////////////////////////////////////////////////////////////////////////////
  /*RooRodenbach RodFull("RodFull","Rod",CMS_hwwlvqq_mWWfull,m,wdth,alpha);
  RooGenericPdf fermiFull("fermiFull","fermi function","1/(1+exp((@1-@0)/@2))",RooArgList(CMS_hwwlvqq_mWWfull,cutOff,beta));
  RooProdPdf backgroundFull("backgroundFull","backgroundFull",RooArgSet(fermiFull,RodFull));
  */
  //first muon then electrons ///for expo /////////////////////////////////
  
  /*EvtNorm.push_back(chan==1? 228.10 : 200.12 );  // 0btag 
  EvtNorm.push_back(chan==1? 230.80 : 195.80 );  // 1btag 
  EvtNorm.push_back(chan==1?  16.82 :  13.79 );  // 2btag
  */
  //for eps///////////////////////////
  /*EvtNorm.push_back(chan==1? 345.7 : 286.4 );  // 0btag 
  EvtNorm.push_back(chan==1? 376.4 : 334.7 );  // 1btag 
  EvtNorm.push_back(chan==1? 24.3 : 20.3 );  // 2btag*/

  //for LP
  double EvtNorm = (chan==1) ? 575.85 : 490.54;  // 0btag //changeeeeee

  string mwwcut2="CMS_hwwlvqq_mWWfull>"+ossmww1.str(); 
  mwwcut2+="&&CMS_hwwlvqq_mWWfull<"+ossmww2.str();
  RooDataHist *BkgHisto = backgroundFull.generateBinned(CMS_hwwlvqq_mWWfull,EvtNorm,kTRUE,kFALSE);  
  exp_yield=float( BkgHisto->sumEntries(mwwcut2.c_str() )  );
  cout<<"MH"<<massH<<"  With this cut: "<<mwwcut2.c_str()<<"  ===> "<<exp_yield<<"   TOT "<< BkgHisto->sumEntries( )<<"  mww<300 "<< BkgHisto->sumEntries("CMS_hwwlvqq_mWWfull<300.0" )<<endl;

  delete file;
  
}


