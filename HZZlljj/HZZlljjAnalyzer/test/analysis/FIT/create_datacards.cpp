#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"

#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "PdfDiagonalizer.h"

#include "SidebandFitter.h"

#include "DrawBase.h"


float mZZmin_ = 183.;


struct HiggsParameters {

  float mH;
  float CSgg;
  float CSgg_p;
  float CSgg_m;
  float CSpdfgg_p;
  float CSpdfgg_m;
  float CSvbf;
  float CSvbf_p;
  float CSvbf_m;
  float CSpdfvbf_p;
  float CSpdfvbf_m;
  float Gamma;
  float BRHZZ;
  float BRZZ2l2q;

};



double sign( double x ) {

  double returnSign = 0.;

  if( x>=0. ) returnSign =  1.;
  else returnSign =  -1.;

  return returnSign;

}


std::string get_datacardDir( const std::string& dataset, const std::string& init, const std::string& flags );
std::string leptType_datacards( const std::string& leptType_str );

void create_singleDatacard( const std::string& dataset, const std::string& PUType, const std::string& init, const std::string& flags, float mass, float lumi, const std::string& leptType_str, int nbtags, TF1* f1_eff_vs_mass );

HiggsParameters get_higgsParameters( float mass );

double linear_interp( double x, double x_old, double mass, double mH, double mH_old );
TF1* get_eff_vs_mass( const std::string& leptType_str, int nbtags, const std::string& PUType, float mZZmin );


double get_signalParameter(int btag, double massH, std::string varname);

std::string systString( std::pair<double,double> systPair, double maxDiff=0.01 );
std::pair<double,double> theorSyst( double errMinus, double errPlus, double addMinus=0., double addPlus=0. );
std::pair<double,double> theorSyst_HighmH( double mHVal );
std::pair<double,double> leptTriggerSyst( const std::string& leptType_str);
std::pair<double,double> leptEffSyst( const std::string& leptType_str);
std::pair<double,double> leptScaleSyst( const std::string& leptType_str);

std::pair<double,double> jetScaleSyst( double mass );
std::pair<double,double> bTagEffSyst( const std::string& leptType_str, int nbtags, double mass );

double backgroundNorm( const std::string& dataset, const std::string& leptType_str, int nbtags );








int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./create_datacards [dataset] [init=\"MC\"] [flags=\"\"]" << std::endl;
    exit(113);
  }

  
  std::string dataset(argv[1]);

  float lumi_ELE;
  float lumi_MU;
  if( dataset=="LP11") {
    lumi_ELE=1556.; //pb^-1
    lumi_MU =1615.; //pb^-1
  } else if( dataset=="Run2011A_FULL" ) {
    lumi_ELE=2100.; //pb^-1
    lumi_MU =2100.; //pb^-1
  } else if( dataset=="HR11" ) {
    lumi_ELE=4200.; //pb^-1
    lumi_MU =4200.; //pb^-1
  } else if( dataset=="HR11_v2" ) {
    lumi_ELE=4600.; //pb^-1
    lumi_MU =4600.; //pb^-1
  } else if( dataset=="Run2011B_v2" ) {
    lumi_ELE=2500.; //pb^-1
    lumi_MU =2500.; //pb^-1
  } else {
    std::cout << "Unknown dataset '" << dataset << "'. Exiting." << std::endl;
    exit(333);
  }


  std::string init = "MC";
  if( argc==3 ) {
    std::string init_str(argv[2]);
    init = init_str;
  }


  std::string flags = "";
  if( argc==4 ) {
    std::string flags_str(argv[3]);
    flags = "_" + flags_str;
  }



  TString dataset_tstr(dataset);
  std::string PUType = "Run2011A";
  if( dataset=="HR11" )
    PUType = "HR11";
  if( dataset=="HR11_v2" )
    PUType = "HR11_73pb";
  if( dataset_tstr.Contains("Run2011B") )
    PUType = "HR11_73pb"; 

  //first loop over available signal MC files to fit efficiency:
  TF1* f1_eff_vs_mass_MU_0btag = get_eff_vs_mass("MU", 0, PUType, mZZmin_);
  TF1* f1_eff_vs_mass_MU_1btag = get_eff_vs_mass("MU", 1, PUType, mZZmin_);
  TF1* f1_eff_vs_mass_MU_2btag = get_eff_vs_mass("MU", 2, PUType, mZZmin_);

  TF1* f1_eff_vs_mass_ELE_0btag = get_eff_vs_mass("ELE", 0, PUType, mZZmin_);
  TF1* f1_eff_vs_mass_ELE_1btag = get_eff_vs_mass("ELE", 1, PUType, mZZmin_);
  TF1* f1_eff_vs_mass_ELE_2btag = get_eff_vs_mass("ELE", 2, PUType, mZZmin_);




  std::ifstream ifs("masses.txt");
  
  while( ifs.good() ) {
    
    float mass;
    ifs >> mass;

    std::cout << std::endl << std::endl;;
    std::cout << "++++++++++++++++++++++" << std::endl;
    std::cout << "+++ MASS: " << mass << std::endl;
    std::cout << "++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::string datacardDir = get_datacardDir( dataset, init, flags );
    char mkdir_command[100];
    sprintf( mkdir_command, "mkdir -p %s/%.0f", datacardDir.c_str(), mass);
    system(mkdir_command);

    create_singleDatacard( dataset, PUType, init, flags, mass, lumi_ELE, "ELE", 0, f1_eff_vs_mass_ELE_0btag);
    create_singleDatacard( dataset, PUType, init, flags, mass, lumi_ELE, "ELE", 1, f1_eff_vs_mass_ELE_1btag);
    create_singleDatacard( dataset, PUType, init, flags, mass, lumi_ELE, "ELE", 2, f1_eff_vs_mass_ELE_2btag);
    create_singleDatacard( dataset, PUType, init, flags, mass, lumi_MU,   "MU", 0, f1_eff_vs_mass_MU_0btag);
    create_singleDatacard( dataset, PUType, init, flags, mass, lumi_MU,   "MU", 1, f1_eff_vs_mass_MU_1btag);
    create_singleDatacard( dataset, PUType, init, flags, mass, lumi_MU,   "MU", 2, f1_eff_vs_mass_MU_2btag);

  } //while masses

  return 0;

}



void create_singleDatacard( const std::string& dataset, const std::string& PUType, const std::string& init, const std::string& flags, float mass, float lumi, const std::string& leptType_str, int nbtags, TF1* f1_eff_vs_mass ) {

  if( leptType_str!="ELE" && leptType_str!="MU" ) {
    std::cout << "Unkown Lept Type '" << leptType_str << "'. Exiting." << std::endl;
    exit(12333);
  }


  SidebandFitter* sf = new SidebandFitter( dataset, PUType, init, flags );

  int leptType_int = SidebandFitter::convert_leptType( leptType_str );

  HiggsParameters hp = get_higgsParameters(mass);

  // open fitResults file (all lept types):
  std::string fitResultsFileName = sf->get_fitResultsRootFileName( nbtags, "ALL" );
  TFile* fitResultsFile = TFile::Open(fitResultsFileName.c_str());

  // get fit result:
  char fitResultName[200];
  sprintf( fitResultName, "fitResults_%dbtag_decorr", nbtags );
  RooFitResult* bgFitResult = (RooFitResult*)fitResultsFile->Get(fitResultName);

  // get workspace:
  char workspaceName[200];
  sprintf( workspaceName, "fitWorkspace_%dbtag", nbtags );
  RooWorkspace* bgws = (RooWorkspace*)fitResultsFile->Get(workspaceName);

  // get sidebands tree:
  TTree* treeSidebandsDATA_alphaCorr = (TTree*)fitResultsFile->Get("sidebandsDATA_alpha");




  //// get main variable from input workspace:
  RooRealVar* CMS_hzz2l2q_mZZ = bgws->var("CMS_hzz2l2q_mZZ");
  CMS_hzz2l2q_mZZ->setMin(mZZmin_); 


  //// redefine mZZ variable, so that we can change its range to 183-800:
  //RooRealVar* CMS_hzz2l2q_mZZ = new RooRealVar("CMS_hzz2l2q_mZZ", "CMS_hzz2l2q_mZZ", 183., 800.);


  
  char suffix[100];
  sprintf( suffix, "%s%s%db", (leptType_datacards(leptType_str)).c_str(), (leptType_datacards(leptType_str)).c_str(), nbtags);
  std::string suffix_str(suffix);

  std::string datacardDir = get_datacardDir( dataset, init, flags );

  char datacardName[400];
  sprintf( datacardName, "%s/%.0f/hzz2l2q_%s.%.0f.txt", datacardDir.c_str(), mass, suffix, mass);


  std::ofstream ofs(datacardName);
  ofs << "# Simple counting experiment, with one signal and one background process" << std::endl;
  ofs << "#imax 1  number of channels" << std::endl;
  ofs << "#jmax 1  number of backgrounds" << std::endl;
  ofs << "#kmax *  number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "shapes ggH CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root  w:signal" << std::endl;
  ofs << "shapes VBF CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root  w:signal" << std::endl;
  ofs << "shapes background CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root w:background_decorr" << std::endl;
  ofs << "shapes data_obs   CMS_hzz2l2q_" << suffix_str << " hzz2l2q_" << suffix_str << ".input.root w:dataset_obs" << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "bin         CMS_hzz2l2q_" << suffix << std::endl;


  RooDataSet* dataset_obs = sf->get_observedDataset( CMS_hzz2l2q_mZZ, leptType_str, nbtags );
  float observedYield = dataset_obs->sumEntries();

  ofs << "observation   " << observedYield << std::endl;
  ofs << "------------ " << std::endl;
  ofs << "bin                CMS_hzz2l2q_" << suffix << "\tCMS_hzz2l2q_" << suffix << "\tCMS_hzz2l2q_" << suffix << std::endl;
  ofs << "process            ggH\t\t\tVBF\t\t\tbackground" << std::endl;
  ofs << "process            -1\t\t\t0\t\t\t1" << std::endl;




  float eff = f1_eff_vs_mass->Eval(hp.mH);
  float rate_gg   = eff*hp.CSgg *hp.BRHZZ*hp.BRZZ2l2q*lumi*0.5; //xsect has both ee and mm
  float rate_vbf  = eff*hp.CSvbf*hp.BRHZZ*hp.BRZZ2l2q*lumi*0.5; //xsect has both ee and mm

  // compute expected BG yield from observed sideband events:
  Double_t rate_background = sf->get_backgroundNormalization( nbtags, leptType_str, "DATA", 183., 800. );



  ofs << "rate               " << rate_gg << "\t\t" << rate_vbf << "\t\t" << rate_background << std::endl;
  ofs << "------------ " << std::endl;


  // and now systematics:

  ofs << "lumi\t\t\tlnN\t1.045\t\t\t1.045\t\t\t1.0" << std::endl;

 
    std::pair<double,double> pdf_gg  = theorSyst( hp.CSpdfgg_m, hp.CSpdfgg_p, 0.04, 0.015 );
    ofs << "pdf_gg\t\tlnN\t" << systString(pdf_gg) << "\t1.0\t\t\t1.0" << std::endl;
    
    std::pair<double,double> pdf_vbf = theorSyst( hp.CSpdfvbf_m, hp.CSpdfvbf_p, 0.04, 0.015 );
    ofs << "pdf_qqbar\t\tlnN\t1.0\t\t\t" << systString(pdf_vbf) << "\t1.0" << std::endl;
    
    std::pair<double,double> QCDscale_ggH = theorSyst( hp.CSgg_m, hp.CSgg_p);
    ofs << "QCDscale_ggH\tlnN\t" << systString(QCDscale_ggH) << "\t1.0\t\t\t1.0" << std::endl;
    
    std::pair<double,double> QCDscale_qqH = theorSyst( hp.CSvbf_m, hp.CSvbf_p);
    ofs << "QCDscale_qqH\tlnN\t1.0\t\t\t" << systString(QCDscale_qqH) << "\t1.0" << std::endl;
    
    std::pair<double,double>theoryUncXS_HighmH  = theorSyst_HighmH( hp.mH);
    ofs << "theoryUncXS_HighMH\tlnN\t"<<systString(theoryUncXS_HighmH)<<"\t\t"<<systString(theoryUncXS_HighmH)<<" \t1.0"<<std::endl;


  ofs << "CMS_trigger_" << leptType_datacards(leptType_str) << "\tlnN\t" << systString(leptTriggerSyst(leptType_str)) << "\t" << systString(leptTriggerSyst(leptType_str)) << "\t1.0" << std::endl;

  ofs << "CMS_eff_" << leptType_datacards(leptType_str) << "\t\tlnN\t" << systString(leptEffSyst(leptType_str)) << "\t" << systString(leptEffSyst(leptType_str)) << "\t1.0" << std::endl;

  ofs << "CMS_scale_" << leptType_datacards(leptType_str) << "\t\tlnN\t" << systString(leptScaleSyst(leptType_str)) << "\t" << systString(leptScaleSyst(leptType_str)) << "\t1.0" << std::endl;

  ofs << "CMS_scale_j\t\tlnN\t" << systString(jetScaleSyst(hp.mH)) << "\t" << systString(jetScaleSyst(hp.mH)) << "\t1.0" << std::endl;

  ofs << "CMS_eff_b\t\tlnN\t" << systString(bTagEffSyst(leptType_str, nbtags, hp.mH)) << "\t" << systString(bTagEffSyst(leptType_str, nbtags, hp.mH)) << "\t1.0" << std::endl;

  ofs << "CMS_hzz2l2q_pu\t\tlnN\t1.04\t\t\t1.04\t\t\t1.0" << std::endl;

  ofs << "CMS_hzz2l2q_qgsep0b\t\tlnN\t1.046\t\t\t1.046\t\t\t1.0" << std::endl;

  

  // syst done. now finish with parameters:

  double bgNorm = backgroundNorm(dataset,leptType_str,nbtags);
  char bgNorm_char[100];
  sprintf( bgNorm_char, "%.0lf", bgNorm);
  std::string bgNorm_str(bgNorm_char);

  double alpha = rate_background/bgNorm;
  char alpha_char[100];
  sprintf( alpha_char, "%lf", alpha);
  std::string alpha_str(alpha_char);

  char bgNormName[200];
  sprintf( bgNormName, "CMS_hzz2l2q_bkg%db%s%sp0", nbtags, (leptType_datacards(leptType_str)).c_str(), (leptType_datacards(leptType_str)).c_str() );
  std::string bgNormName_str(bgNormName);
  ofs << bgNormName_str << "\tgmN " << bgNorm_str << "\t-----\t-----\t" << alpha_str << std::endl;


 
  RooArgList bgPars = bgFitResult->floatParsFinal();

  for( unsigned iVar=0; iVar<bgPars.getSize(); ++iVar ) {
    RooRealVar* thisVar = dynamic_cast<RooRealVar*>(bgPars.at(iVar));
//     std::string varn(thisVar->GetName());
//     if(varn=="CMS_hzz2l2q_bkg_1b_eig1" )
//       ofs << thisVar->GetName() << "\tparam\t\t" << thisVar->getVal() << "\t" << thisVar->getError()*2 << std::endl;
//     else
      ofs << thisVar->GetName() << "\tparam\t\t" << thisVar->getVal() << "\t" << thisVar->getError() << std::endl;
  }



  ofs.close();
  fitResultsFile->Close();


  
  std::cout << std::endl << std::endl;
  std::cout << "+++ DATACARD FOR MASS " << mass << " ( " << nbtags << " B-TAGS, " << leptType_str << " CHANNEL ) IS DONE." << std::endl;
  std::cout << std::endl;

  // datacard is done. now create output workspace and write it to rootfile

  char outfileName[900];
  sprintf( outfileName, "%s/%.0f/hzz2l2q_%s.input.root", datacardDir.c_str(), mass, suffix);
  TFile* outfile = TFile::Open( outfileName, "RECREATE");
  outfile->cd();


  RooWorkspace* w = new RooWorkspace("w","w");
  w->addClassDeclImportDir("/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
  //w->addClassDeclImportDir("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/HZZlljj/HZZlljjAnalyzer/test/analysis/FIT/PDFs");



  // import variable in output workspace:
  w->import(*CMS_hzz2l2q_mZZ);

  // import observed dataset:
  w->import(*dataset_obs);

  // get BG shape:
  RooAbsPdf* background_decorr = bgws->pdf("background_decorr");
 
  // and import it:
  w->import(*background_decorr, RooFit::RecycleConflictNodes());



//// now define signal shape:
//RooAbsPdf* signal = get_signalShape( CMS_hzz2l2q_mZZ, nbtags, hp.mH );


  // now define signal shape (didn manage to do use get_signalShape without a crash):

  // ------------------- Crystal Ball (matched) -------------------------------
  float massH = hp.mH;
  char sigp1name[200];
  char sigp2name[200];
  sprintf( sigp1name, "CMS_hzz2l2q_sig%dbp1", nbtags ); //m
  sprintf( sigp2name, "CMS_hzz2l2q_sig%dbp2", nbtags ); //width
  RooRealVar CB_mean(sigp1name,sigp1name, get_signalParameter(nbtags,massH,"matched_MassCBmean"));
  RooRealVar CB_sigma(sigp2name,sigp2name,get_signalParameter(nbtags,massH,"matched_MassCBsigma"));
  RooRealVar CB_alpha1("CB_alpha1","param 3 of CB",get_signalParameter(nbtags,massH,"matched_MassCBalpha1"));
  RooRealVar CB_n1("CB_n1","param 4 of CB",get_signalParameter(nbtags,massH,"matched_MassCBn1"));
  RooRealVar CB_alpha2("CB_alpha2","param 3 of CB",get_signalParameter(nbtags,massH,"matched_MassCBalpha2"));
  RooRealVar CB_n2("CB_n2","param 5 of CB",get_signalParameter(nbtags,massH,"matched_MassCBn2"));

  RooDoubleCB* CB_SIG = new RooDoubleCB("CB_SIG","Crystal Ball",*CMS_hzz2l2q_mZZ,CB_mean,CB_sigma,CB_alpha1,CB_n1,CB_alpha2,CB_n2);


  // ------------------- SmearedTriangle (un-matched) -------------------------------
  RooRealVar CB_UMmean( "CB_UMmean"," CB_UMmean", get_signalParameter(nbtags,massH,"unmatched_MassCBmean"));
  RooRealVar CB_UMsigma("CB_UMsigma","CB_UMsigma",get_signalParameter(nbtags,massH,"unmatched_MassCBsigma"));
  RooRealVar CB_UMalpha("CB_UMalpha","CB_UMalpha",get_signalParameter(nbtags,massH,"unmatched_MassCBalpha"));
  RooRealVar CB_UMn("CB_UMn","CB_UMn",get_signalParameter(nbtags,massH,"unmatched_MassCBn"));
  RooCBShape* CB_UM = new RooCBShape("CB_UM","Crystal Ball unmacthed",*CMS_hzz2l2q_mZZ,CB_UMmean,CB_UMsigma ,CB_UMalpha,CB_UMn);

  RooRealVar TRI_start("TRI_start","TRI_start", get_signalParameter(nbtags,massH,"unmatched_Mass_start"));
  RooRealVar TRI_turn("TRI_turn","TRI_turn", get_signalParameter(nbtags,massH,"unmatched_Mass_turn"));
  RooRealVar TRI_stop("TRI_stop","TRI_stop", get_signalParameter(nbtags,massH,"unmatched_Mass_stop"));
  Triangle* TRI = new Triangle("TRI","TRI",*CMS_hzz2l2q_mZZ,TRI_start,TRI_turn,TRI_stop);

  //------------------------ convolution -------------------------

  //CMS_hzz2l2q_mZZ->setBins(10000,"fft");

  RooFFTConvPdf* TRI_SMEAR = new RooFFTConvPdf("TRI_SMEAR","triangle (X) CB",*CMS_hzz2l2q_mZZ,*TRI,*CB_UM);
  TRI_SMEAR->setBufferFraction(1.0);


  //------------------------ add matched and unmatched -------------------------
  RooRealVar MATCH("MATCH","MATCH", get_signalParameter(nbtags,massH,"N_matched"));
  RooAddPdf* signal = new RooAddPdf("signal","signal",*CB_SIG,*TRI_SMEAR,MATCH);




  // and import it:
  w->import(*signal, RooFit::RecycleConflictNodes());



  // done. now save:
  w->Write();
  outfile->Close();

  delete sf;

}










TF1* get_eff_vs_mass( const std::string& leptType_str, int nbtags, const std::string& PUType, float mZZmin ) {


  TH1F::AddDirectory(kTRUE);

  int leptType_int = SidebandFitter::convert_leptType(leptType_str);

  ifstream ifsMC("massesMC.txt"); //the points at which we have MC samples

  TGraph* gr_eff_vs_mass = new TGraph(0);

  int iPoint=0;
  
  while( ifsMC.good() ) {

    double mass;
    ifsMC >> mass;

    HiggsParameters hp = get_higgsParameters(mass);


    char signalfileName[800];
    sprintf( signalfileName, "HZZlljjRM_GluGluToHToZZTo2L2Q_M-%.0f_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1_optLD_looseBTags_v2_PU%s_ALL.root", hp.mH, PUType.c_str() );

    TFile* signalFile = TFile::Open(signalfileName);
    TTree* signalTree = (TTree*)signalFile->Get("tree_passedEvents");

    char signalCut[500];
    sprintf( signalCut, "HLTSF*PUWeight*( mZjj>75. && mZjj<105. && CMS_hzz2l2q_mZZ>%f && CMS_hzz2l2q_mZZ<800. && nBTags==%d && leptType==%d)", mZZmin_, nbtags, leptType_int);
    TH1D* h1_mZZ_signal = new TH1D("mZZ_signal", "", 65, 150., 800.);
    h1_mZZ_signal->Sumw2();
    signalTree->Project( "mZZ_signal", "CMS_hzz2l2q_mZZ", signalCut );

    float signalYield = h1_mZZ_signal->Integral();

    TH1F* h1_generatedEventsPU = (TH1F*)signalFile->Get("nCounterPU");
    float generatedYield = h1_generatedEventsPU->Integral();

    float efficiency = signalYield/(generatedYield/3.); //three lept types in powheg

    gr_eff_vs_mass->SetPoint( iPoint++, mass, efficiency );

    signalFile->Close();

  } //while masses


  char functName[200];
  sprintf( functName, "eff_vs_mass_%s_%dbtag", leptType_str.c_str(), nbtags );
  TF1* f1_eff_vs_mass = new TF1(functName, "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 150., 700.);
  gr_eff_vs_mass->Fit(f1_eff_vs_mass, "RQN");
  f1_eff_vs_mass->SetLineStyle(2);
  f1_eff_vs_mass->SetLineColor(38);

  // define an instance of DrawBase to set the style
  DrawBase* db = new DrawBase("eff");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* axes = new TH2D("axes", "", 10, 160., 640., 10, 0., 0.115);
  axes->SetXTitle("m_{H} [GeV]");
  axes->SetYTitle("Efficiency");
  axes->Draw();

  gr_eff_vs_mass->SetMarkerStyle(20);
  gr_eff_vs_mass->SetMarkerSize(1.3);
  gr_eff_vs_mass->SetMarkerColor(46);
  f1_eff_vs_mass->Draw("same");
  gr_eff_vs_mass->Draw("Psame");

  TPaveText* labelCMS = db->get_labelCMS();
  TPaveText* labelSqrt = db->get_labelSqrt();

  labelCMS->Draw("same");
  labelSqrt->Draw("same");

  std::string leptType_forlabel = (leptType_str=="ELE") ? "electron" : "muon";
  char channelLabelText[300];
  sprintf(channelLabelText, "%d b-tag category (%s channel)", nbtags, leptType_forlabel.c_str() );

  TPaveText* labelChannel = new TPaveText( 0.4, 0.2, 0.85, 0.25, "brNDC");
  labelChannel->SetTextSize(0.035);
  labelChannel->SetFillColor(0);
  labelChannel->AddText(channelLabelText);
  labelChannel->Draw("same");

  char effDirName[200];
  sprintf( effDirName, "EfficiencyFits_PU%s", PUType.c_str() );
  char mkdirCommand[200];
  sprintf( mkdirCommand, "mkdir -p %s", effDirName );
  
  system( mkdirCommand ); 

  char canvasName[500];
  sprintf( canvasName, "%s/effFit_%s_%dbtag.eps", effDirName, leptType_str.c_str(), nbtags);
  c1->SaveAs(canvasName);

  delete c1;

  return f1_eff_vs_mass;

}




std::string leptType_datacards( const std::string& leptType_str ) {
 
  std::string returnString="";

  if( leptType_str=="MU" ) returnString = "m";
  if( leptType_str=="ELE" ) returnString = "e";

  return returnString;

}


HiggsParameters get_higgsParameters( float mass ) {

  std::string nameXsecFile = "xsect_higgs_173points_new.txt";
  std::ifstream xsect_file(nameXsecFile.c_str());

  if (! xsect_file.is_open()) { 
    std::cout << "Failed to open file with xsections"<<endl;
    exit(13111);
  }

  xsect_file.clear();
  xsect_file.seekg(0);
 
  float mH, CSgg, CSgg_p, CSgg_m, CSpdfgg_p,CSpdfgg_m,CSvbf, CSvbf_p, CSvbf_m,CSpdfvbf_p,CSpdfvbf_m, 
        Gamma, BRHZZ, BRZZ2l2q;

  float mH_old, CSgg_old, CSgg_p_old, CSgg_m_old, CSpdfgg_p_old,CSpdfgg_m_old,CSvbf_old, CSvbf_p_old, CSvbf_m_old,CSpdfvbf_p_old,CSpdfvbf_m_old, 
        Gamma_old, BRHZZ_old, BRZZ2l2q_old;

  HiggsParameters hp;

  while(xsect_file.good()) {

    mH_old = mH;
    CSgg_old = CSgg;
    CSgg_p_old = CSgg_p;
    CSgg_m_old = CSgg_m;
    CSpdfgg_p_old = CSpdfgg_p;
    CSpdfgg_m_old = CSpdfgg_m;
    CSvbf_old = CSvbf;
    CSvbf_p_old = CSvbf_p;
    CSvbf_m_old = CSvbf_m;
    CSpdfvbf_p_old = CSpdfvbf_p;
    CSpdfvbf_m_old = CSpdfvbf_m;
    Gamma_old = Gamma;
    BRHZZ_old = BRHZZ;
    BRZZ2l2q_old = BRZZ2l2q;

    xsect_file >> mH >> CSgg>> CSgg_p >> CSgg_m >> CSpdfgg_p >> CSpdfgg_m >> CSvbf >> CSvbf_p >> CSvbf_m >> CSpdfvbf_p >> CSpdfvbf_m >> Gamma >> BRHZZ >> BRZZ2l2q;

    if( mH == mass ) {

      hp.mH = mH;
      hp.CSgg = CSgg;
      hp.CSgg_p = CSgg_p;
      hp.CSgg_m = CSgg_m;
      hp.CSpdfgg_p = CSpdfgg_p;
      hp.CSpdfgg_m = CSpdfgg_m;
      hp.CSvbf = CSvbf;
      hp.CSvbf_p = CSvbf_p;
      hp.CSvbf_m = CSvbf_m;
      hp.CSpdfvbf_p = CSpdfvbf_p;
      hp.CSpdfvbf_m = CSpdfvbf_m;
      hp.Gamma = Gamma;
      hp.BRHZZ = BRHZZ;
      hp.BRZZ2l2q = BRZZ2l2q;
    
      break;

    } if( mass > mH_old && mass < mH ) {

      hp.mH = mass;
      hp.CSgg       = linear_interp( CSgg, CSgg_old, mass, mH, mH_old );
      hp.CSgg_p     = linear_interp( CSgg_p, CSgg_p_old, mass, mH, mH_old );
      hp.CSgg_m     = linear_interp( CSgg_m, CSgg_m_old, mass, mH, mH_old );
      hp.CSpdfgg_p  = linear_interp( CSpdfgg_p, CSpdfgg_p_old, mass, mH, mH_old );
      hp.CSpdfgg_m  = linear_interp( CSpdfgg_m, CSpdfgg_m_old, mass, mH, mH_old );
      hp.CSvbf      = linear_interp( CSvbf, CSvbf_old, mass, mH, mH_old );
      hp.CSvbf_p    = linear_interp( CSvbf_p, CSvbf_p_old, mass, mH, mH_old );
      hp.CSvbf_m    = linear_interp( CSvbf_m, CSvbf_m_old, mass, mH, mH_old );
      hp.CSpdfvbf_p = linear_interp( CSpdfvbf_p, CSpdfvbf_p_old, mass, mH, mH_old );
      hp.CSpdfvbf_m = linear_interp( CSpdfvbf_m, CSpdfvbf_m_old, mass, mH, mH_old );
      hp.Gamma      = linear_interp( Gamma, Gamma_old, mass, mH, mH_old );
      hp.BRHZZ      = linear_interp( BRHZZ, BRHZZ_old, mass, mH, mH_old );
      hp.BRZZ2l2q   = linear_interp( BRZZ2l2q, BRZZ2l2q_old, mass, mH, mH_old );
    
      break;

    } // if

  } //while ifs

  return hp;

}




double linear_interp( double x, double x_old, double mass, double mH, double mH_old ) {

  return (x_old + ( x-x_old ) * ( mass-mH_old ) / ( mH-mH_old ));

}



/*
RooAbsPdf* get_signalShape( RooRealVar* CMS_hzz2l2q_mZZ, int nbtags, float massH ) {


  // ------------------- Crystal Ball (matched) -------------------------------
  char sigp1name[200];
  char sigp2name[200];
  sprintf( sigp1name, "CMS_hzz2l2q_sig%dbp1", nbtags ); //m
  sprintf( sigp2name, "CMS_hzz2l2q_sig%dbp2", nbtags ); //width
  RooRealVar CB_mean(sigp1name,sigp1name, get_signalParameter(nbtags,massH,"matched_MassCBmean"));
  RooRealVar CB_sigma(sigp2name,sigp2name,get_signalParameter(nbtags,massH,"matched_MassCBsigma"));
  RooRealVar CB_alpha1("CB_alpha1","param 3 of CB",get_signalParameter(nbtags,massH,"matched_MassCBalpha1"));
  RooRealVar CB_n1("CB_n1","param 4 of CB",get_signalParameter(nbtags,massH,"matched_MassCBn1"));
  RooRealVar CB_alpha2("CB_alpha2","param 3 of CB",get_signalParameter(nbtags,massH,"matched_MassCBalpha2"));
  RooRealVar CB_n2("CB_n2","param 5 of CB",get_signalParameter(nbtags,massH,"matched_MassCBn2"));

  RooDoubleCB* CB_SIG = new RooDoubleCB("CB_SIG","Crystal Ball",*CMS_hzz2l2q_mZZ,CB_mean,CB_sigma,CB_alpha1,CB_n1,CB_alpha2,CB_n2);


  // ------------------- SmearedTriangle (un-matched) -------------------------------
  RooRealVar CB_UMmean( "CB_UMmean"," CB_UMmean", get_signalParameter(nbtags,massH,"unmatched_MassCBmean"));
  RooRealVar CB_UMsigma("CB_UMsigma","CB_UMsigma",get_signalParameter(nbtags,massH,"unmatched_MassCBsigma"));
  RooRealVar CB_UMalpha("CB_UMalpha","CB_UMalpha",get_signalParameter(nbtags,massH,"unmatched_MassCBalpha"));
  RooRealVar CB_UMn("CB_UMn","CB_UMn",get_signalParameter(nbtags,massH,"unmatched_MassCBn"));
  RooCBShape* CB_UM = new RooCBShape("CB_UM","Crystal Ball unmacthed",*CMS_hzz2l2q_mZZ,CB_UMmean,CB_UMsigma ,CB_UMalpha,CB_UMn);

  RooRealVar TRI_start("TRI_start","TRI_start", get_signalParameter(nbtags,massH,"unmatched_Mass_start"));
  RooRealVar TRI_turn("TRI_turn","TRI_turn", get_signalParameter(nbtags,massH,"unmatched_Mass_turn"));
  RooRealVar TRI_stop("TRI_stop","TRI_stop", get_signalParameter(nbtags,massH,"unmatched_Mass_stop"));
  Triangle* TRI = new Triangle("TRI","TRI",*CMS_hzz2l2q_mZZ,TRI_start,TRI_turn,TRI_stop);

  //------------------------ convolution -------------------------
  
  //CMS_hzz2l2q_mZZ.setBins(10000,"fft");

  RooFFTConvPdf* TRI_SMEAR = new RooFFTConvPdf("TRI_SMEAR","triangle (X) CB",*CMS_hzz2l2q_mZZ,*TRI,*CB_UM);
  TRI_SMEAR->setBufferFraction(1.0);
  

  //------------------------ add matched and unmatched -------------------------
  RooRealVar MATCH("MATCH","MATCH", get_signalParameter(nbtags,massH,"N_matched"));
  RooAddPdf* signal = new RooAddPdf("signal","signal",*CB_SIG,*TRI_SMEAR,MATCH);

  return signal;  

}
*/



double get_signalParameter(int btag, double massH, std::string varname) {

  int masses[18] = {190,200,210,230,250,275,300,325,350,375,400,425,475,500,525,550,575,600};

  RooRealVar var(varname.c_str(),varname.c_str(),0.);
  RooArgSet paramsup, paramslow;

  paramsup.add(var);
  paramslow.add(var);

  char filename[200];

  //which files to read?
  for(int i =0 ; i <18 ; i++){
    if(masses[i]==massH){//direct Match
      sprintf(filename,"signalFitResults/out-%d-%s-btag%d.config",masses[i],"EM",btag);
      //std::cout << filename << " : " << paramsup.readFromFile(filename, "READ", "Parameters") <<std::endl;
      paramsup.readFromFile(filename, "READ", "Parameters");
      return var.getVal();
    }
  }

  //no direct match => interpolation
  int indexlow = -1;
  int indexhigh= -1;
  for(int i =0 ; i <18 ; i++){
    if(masses[i]>massH){
      indexhigh=i;
      break;
    }
  }
  for(int i =17 ; i >-1 ; i--){
    if(masses[i]<massH){
      indexlow=i;
      break;
    }
  }
  if(indexlow==-1 || indexhigh== -1){
    std::cout << "Out of Range"<< std::endl;
    exit(1);
  }

  //std::cout << indexlow << " " << indexhigh <<std::endl;

  sprintf(filename,"signalFitResults/out-%d-%s-btag%d.config",masses[indexlow],"EM",btag);
  paramsup.readFromFile(filename, "READ", "Parameters");
  double low = var.getVal();
  sprintf(filename,"signalFitResults/out-%d-%s-btag%d.config",masses[indexhigh],"EM",btag);
  paramsup.readFromFile(filename, "READ", "Parameters");
  double high = var.getVal();
  
  double deltaM = masses[indexhigh] - masses[indexlow];
  
  return (massH-masses[indexlow])/deltaM*(high-low) + low;
  
}




std::string systString( std::pair<double,double> systPair, double maxDiff ) {

  double syst_ave = 1. + 0.5*(fabs(systPair.first-1.) + fabs(systPair.second-1.));
  
  char syst_char[100];
  if( fabs(syst_ave-systPair.second)/syst_ave < maxDiff )
    sprintf( syst_char, "%f    ", syst_ave );
  else
    sprintf( syst_char, "%f/%f", systPair.first, systPair.second );

  std::string syst_str(syst_char);

  return syst_str;

}
 


std::pair<double,double> theorSyst( double errMinus, double errPlus, double addMinus, double addPlus ) {

  float systPlus  = sign(errPlus) *sqrt(errPlus*errPlus   + addPlus*addPlus);
  float systMinus = sign(errMinus)*sqrt(errMinus*errMinus + addMinus*addMinus);

  systPlus  += 1.;
  systMinus += 1.;

  std::pair<double,double> returnPair;
  returnPair.first = systMinus;
  returnPair.second = systPlus;

  return returnPair;

}

std::pair<double,double> theorSyst_HighmH( double mHVal){

 //extra error due to Higgs width
  double theoryHighMass = 1.5*(mHVal/1000)*(mHVal/1000.0)*(mHVal/1000.0);//mHVal=Higgs mass in GeV
  double systPlus=1.0*theoryHighMass;
  double systMinus=-1.0*theoryHighMass;

  systPlus  += 1.;
  systMinus += 1.;
  
  std::pair<double,double> returnPair;
  returnPair.first = systMinus;
  returnPair.second = systPlus;

  return returnPair;

}

std::pair<double,double> leptTriggerSyst( const std::string& leptType_str) {

  double syst;

  if( leptType_str=="MU" )  syst = 1.02;
  if( leptType_str=="ELE" ) syst = 1.01;

  std::pair<double,double> returnPair;
  returnPair.first  = syst; //symmetrical for now
  returnPair.second = syst;

  return returnPair;

}

std::pair<double,double> leptEffSyst( const std::string& leptType_str) {

  double syst;

  if( leptType_str=="MU" )  syst = 1.008;
  if( leptType_str=="ELE" ) syst = 1.034;

  std::pair<double,double> returnPair;
  returnPair.first  = syst; //symmetrical for now
  returnPair.second = syst;

  return returnPair;

}

std::pair<double,double> leptScaleSyst( const std::string& leptType_str) {

  double syst;

  if( leptType_str=="MU" )  syst = 1.01;
  if( leptType_str=="ELE" ) syst = 1.03;

  std::pair<double,double> returnPair;
  returnPair.first  = syst; //symmetrical for now
  returnPair.second = syst;

  return returnPair;

}


std::pair<double,double> jetScaleSyst( double mass ) {

  float p0= 8.3  , p1=-0.0215 ;
  float m0=-8.6, m1=0.02 ;

  std::pair<double,double> returnPair;
  returnPair.first  = 1.0 + 0.01*(m0+m1*mass);
  returnPair.second = 1.0 + 0.01*(p0+p1*mass);

  return returnPair;

}


std::pair<double,double> bTagEffSyst( const std::string& leptType_str, int nbtags, double mass ) {

  float p0=0.0, p1=0.0;
  float m0=0.0, m1=0.0;

  if( leptType_str=="ELE" ) {
    if(nbtags==0){
      p0=0.98644 ;//0.983256647923;
      p1=-0.000100716 ;//-0.0000883532570978;
      m0= 1.0176 ;//1.02907356239;
      m1=8.12716e-05 ;//0.0000713061639147;
    }
    else if(nbtags==1){
      p0=1.04847 ;//1.04446110956;
      p1=-1.47685e-05 ;//-0.0000195508160829;
      m0=0.943369 ;//0.940063743731;
      m1=4.63272e-05 ;//0.0000737044467898;
    }
    else if(nbtags==2){
      p0=1.1035 ;//1.13365470372;
      p1=1.74066e-05 ;//0.00000584572717646;
      m0=0.885855 ;//0.82161771535;
      m1=-4.08034e-07 ;//-0.0000161054152592;
    }
  }
  else if(leptType_str=="MU" ) {
    if(nbtags==0){
      p0=0.990292 ;//0.984636818312;
      p1=-0.000108687;//-0.0000898705296203;
      m0= 1.02153;//1.02836905579;
      m1= 7.44788e-05 ;//0.0000726807344479;
    }
    else if(nbtags==1){
      p0=1.04438 ;//1.04385580002;
      p1=-9.26161e-06 ;//0.0000206096278947;
      m0=1.02153 ;//0.938029 ;//0.942713582987;
      m1=5.95741e-05 ;//0.0000719882385098;
    }
    else if(nbtags==2){
      p0=1.11229;//1.1333366687;
      p1=6.55678e-07 ;//0.00000462542786413;
      m0=0.897045 ;//0.813316607701;
      m1=-2.25839e-05;//-0.00000205840248842;
    }
  }

  std::pair<double,double> returnPair;
  returnPair.first  = m1*mass+m0;
  returnPair.second = p1*mass+p0;

  return returnPair;

}



double backgroundNorm( const std::string& dataset, const std::string& leptType_str, int nbtags ) {

  int leptType_int = SidebandFitter::convert_leptType( leptType_str );

  std::string fileName = "HZZlljjRM_DATA_" + dataset + "_optLD_looseBTags_v2_ALL.root";
  TFile* file_data = TFile::Open(fileName.c_str());
  TTree* tree = (TTree*)file_data->Get("tree_passedEvents");
  char selection[400];
  sprintf(selection, "isSidebands && leptType==%d && nBTags==%d", leptType_int, nbtags);
  float nEvents_sidebands = tree->GetEntries(selection);


  return nEvents_sidebands;

}


std::string get_datacardDir( const std::string& dataset, const std::string& init, const std::string& flags ) {

  char datacardDir[400];
  sprintf( datacardDir, "datacards_%s_fit%s%s", dataset.c_str(), init.c_str(), flags.c_str());

  std::string datacardDir_str(datacardDir);

  return datacardDir_str;

}
