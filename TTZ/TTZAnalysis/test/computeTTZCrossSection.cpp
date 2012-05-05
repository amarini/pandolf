#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "RooHistError.h"

#include "CommonTools/StatTools.h"



int main( int argc, char* argv[] ) {


  std::string sel = "sel1";
  if( argc>1 ) {
    std::string sel_tmp(argv[1]);
    sel = sel_tmp;
  }

  std::string suffix =  sel + "_TCHE_ALL.root";

  //std::string dataFileName = = "TTZTrilepton_DATA_Run2011_FULL_" + suffix;
  //TFile* dataFile = TFile::Open(dataFileName.c_str());

  std::string signalFileName = "TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi_" + suffix;
  TFile* signalFile = TFile::Open(signalFileName.c_str());
  TTree* signalTree = (TTree*)signalFile->Get("tree_passedEvents");

  //std::string mcZJetsFile = "TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11_" + suffix;
  //std::string mcTTbarFile = "TTZTrilepton_TTJ_Fall11_highstat_" + suffix;
  //std::string mcVVFile = "TTZTrilepton_VV_Summer11_" + suffix;
  //std::string mcTTWFile = "TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi_" + suffix;
  //bgTree->Add(dyTreeName.c_str());
  //bgTree->Add(ttTreeName.c_str());
  //bgTree->Add(dibosonTreeName.c_str());


  std::string yieldFileName = "TTZTrileptonPlots_DATA_Run2011_FULL_" + sel + "_TCHE_ALL/yields.txt";
  ifstream ifs(yieldFileName.c_str());

  std::cout << "-> Opened yield file '" << yieldFileName << "'." << std::endl;
  std::string channel;
  std::string obs_str;
  std::string s_str, b_pred_str, b_pred_err_str;
  bool found = false;
  while( ifs.good() ) {
    ifs >> channel >> obs_str >> s_str >> b_pred_str >> b_pred_err_str;
    if( channel == "Total" ) found = true;
  }

  int obs = atoi(obs_str.c_str());
  float s = atof(s_str.c_str());
  float b_pred = atof(b_pred_str.c_str());
  float b_pred_err = atof(b_pred_err_str.c_str());

  if( !found ) {
    std::cout << "There must be a problem. Didn't find total." << std::endl;
    exit(3431);
  } 
   

  // stat error on observed:
  double obs_plus, obs_minus;
  RooHistError::instance().getPoissonInterval(obs,obs_minus,obs_plus,1.);
  double obs_errPlus = obs_plus-obs;
  double obs_errMinus = obs-obs_minus;

  //double b_pred_err = 0.2*b_pred; //random value for now

  float ZBi = StatTools::computeZBi( s+b_pred, b_pred, b_pred_err );
  float ZBi_obs = StatTools::computeZBi( obs, b_pred, b_pred_err );

  float obs_ttz = obs - b_pred;

  TH1D* h1_nCounter_ttZ = (TH1D*)signalFile->Get("nCounter");
  float nTotal_ttz = h1_nCounter_ttZ->GetBinContent(1);

  float lumi_pb = 4980.;
  float crossSection_ttz = 0.139;

  float n_passed_ttz = nTotal_ttz*s/(crossSection_ttz*lumi_pb);
  float eff_ttz = n_passed_ttz/nTotal_ttz;

  float crossSection = obs_ttz / ( lumi_pb*eff_ttz );

  float xsecErr_stat_plus = obs_errPlus / ( lumi_pb*eff_ttz );
  float xsecErr_stat_minus = obs_errMinus / ( lumi_pb*eff_ttz );

  float lumi_err = 0.025*lumi_pb;
  float xsecErr_syst_lumi = lumi_err * crossSection/lumi_pb;

  float xsecErr_syst = xsecErr_syst_lumi; //just lumi for now

  float xsecErr_tot_plus = sqrt( xsecErr_stat_plus*xsecErr_stat_plus + xsecErr_syst*xsecErr_syst );
  float xsecErr_tot_minus = sqrt( xsecErr_stat_minus*xsecErr_stat_minus + xsecErr_syst*xsecErr_syst );


  std::cout << std::endl << std::endl << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
  std::cout << " Expected BG: " << b_pred << " +- " << b_pred_err << std::endl;
  std::cout << " Expected Signal: " << s << " (eff=" << eff_ttz*100. << "%)" << std::endl;
  std::cout << " Expected B+S: " << s+b_pred << std::endl;
  std::cout << " Expected ZBi: " << ZBi << std::endl;
  std::cout << " Observed Events: " << obs << std::endl;
  std::cout << " ZBi (observed): " << ZBi_obs << std::endl;
  std::cout << " Observed Signal (BG subtracted): " << obs_ttz << std::endl;
  std::cout << " Expected Cross Section: " << crossSection_ttz << " pb" << std::endl;
  std::cout << " Measured Cross Section: " << std::endl;
  std::cout << " " <<  crossSection << "  +" << xsecErr_stat_plus << "/-" << xsecErr_stat_minus << " (stat)   +-" << xsecErr_syst << " (syst)  pb" << std::endl;
  std::cout << " " <<  crossSection << "  +" << xsecErr_tot_plus << "/-" << xsecErr_tot_minus << " pb" << std::endl;
  std::cout << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;


  return 0;

}




