#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"


std::pair< float, float >  getSyst( const std::string& syst, const std::string& sel, const std::string& bg_signal );


int main(int argc, char* argv[]) {

  std::string selection = "optsel3";
  if( argc>1 ) {
    std::string sel_tmp(argv[1]);
    selection = sel_tmp;
  }

  std::string dir = "TTZTrileptonPlots_DATA_Run2011_FULL_" + selection + "_TCHE_ALL/";

  std::string yieldFileName = dir + "yields.txt";

  std::ifstream ifs(yieldFileName.c_str());

  std::cout << "-> Opened yield file '" << yieldFileName << "'." << std::endl;
  std::string channel;
  std::string obs_str;
  std::string s_str, b_pred_str, b_pred_err_str;
  bool found = false;
  while( ifs.good() ) {
    ifs >> channel >> obs_str >> s_str >> b_pred_str >> b_pred_err_str;
    if( channel == "Total" ) found = true;
  }

  if( !found ) {
    std::cout << "There must be a problem. Didn't find total." << std::endl;
    exit(3431);
  } 

  int obs = atoi(obs_str.c_str());
  float s = atof(s_str.c_str());
  float b_pred = atof(b_pred_str.c_str());
  float b_pred_err = atof(b_pred_err_str.c_str());
   
  
  // start creating datacard
  std::string datacardName = dir + "datacard.txt";
  std::ofstream datacard(datacardName.c_str());

  datacard << "#imax 1  number of channels" << std::endl;
  datacard << "#jmax 1  number of backgrounds" << std::endl;
  datacard << "#kmax *  number of nuisance parameters" << std::endl;
  datacard << "imax 1" << std::endl;
  datacard << "jmax 1" << std::endl;
  datacard << "kmax *" << std::endl;

  datacard << std::endl << std::endl;

  datacard << "bin         \t1" << std::endl;
  datacard << "observation \t" << obs << std::endl;

  datacard << std::endl << std::endl;

  datacard << "bin         \t1      \t\t1" << std::endl; 
  datacard << "process     \tttZ    \t\tbg" << std::endl; 
  datacard << "process     \t0      \t\t1" << std::endl; 
  datacard << "rate        \t" << s << "\t\t" << b_pred << std::endl; 

  datacard << std::endl << std::endl;

  datacard << "#syst" << std::endl;
  datacard << "lumi     lnN\t1.022  \t\t1.022" << std::endl; //taken from SMP-12-008
  datacard << "bgUncert lnN\t-      \t\t" << 1. + b_pred_err/b_pred << std::endl;

  std::pair< float, float >  leptSystBG = getSyst( "Lept", selection, "BG" );
  std::pair< float, float >  leptSystSignal = getSyst( "Lept", selection, "Signal" );
  datacard << "lept     lnN\t" <<  leptSystSignal.first << "/" << leptSystSignal.second << "\t" << leptSystBG.first << "/" << leptSystBG.second << std::endl;

  std::pair< float, float >  btagSystBG = getSyst( "BTag", selection, "BG" );
  std::pair< float, float >  btagSystSignal = getSyst( "BTag", selection, "Signal" );
  datacard << "btag     lnN\t" <<  btagSystSignal.first << "/" << btagSystSignal.second << "\t" << btagSystBG.first << "/" << btagSystBG.second << std::endl;

  std::pair< float, float >  jesSystBG = getSyst( "JES", selection, "BG" );
  std::pair< float, float >  jesSystSignal = getSyst( "JES", selection, "Signal" );
  datacard << "jes      lnN\t" <<  jesSystSignal.first << "/" << jesSystSignal.second << "\t" << jesSystBG.first << "/" << jesSystBG.second << std::endl;

  std::pair< float, float >  jerSystBG = getSyst( "JER", selection, "BG" );
  std::pair< float, float >  jerSystSignal = getSyst( "JER", selection, "Signal" );
  datacard << "jer      lnN\t" << jerSystSignal.second << "\t\t" << jerSystBG.second << std::endl;

  datacard.close();

  std::cout << "-> Created datacard: " << datacardName << std::endl;

  return 0;

}




std::pair<float, float>  getSyst( const std::string& syst, const std::string& sel, const std::string& bg_signal ) {


  if( bg_signal!="BG" && bg_signal!="Signal" ) {
    std::cout << "bg_signal must be either 'BG' or 'Signal'. Exiting." << std::endl;
    exit(73);
  }

  std::string dataset = (bg_signal=="BG") ? bg_signal : "TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  

  std::string systFile = "TTZTrilepton_" + dataset + "_" + sel + "_TCHE_ALL.root";
  std::string systFileUP = "TTZTrilepton_" + dataset + "_" + sel + "_TCHE_ALL_" + syst + "UP.root";
  std::string systFileDOWN = "TTZTrilepton_" + dataset + "_" + sel + "_TCHE_ALL_" + syst + "DOWN.root";

  TFile* file_systFile = TFile::Open( systFile.c_str() );
  TFile* file_systFileUP = TFile::Open( systFileUP.c_str() );
  TFile* file_systFileDOWN = TFile::Open( systFileDOWN.c_str() );

  TH1D* h1_mean = (TH1D*)file_systFile->Get("channelYields");
  TH1D* h1_systUP = (TH1D*)file_systFileUP->Get("channelYields");
  TH1D* h1_systDOWN = (file_systFileDOWN!=0) ? (TH1D*)file_systFileDOWN->Get("channelYields") : 0;

  float int_mean = h1_mean->Integral();
  float int_systUP = h1_systUP->Integral();
  float int_systDOWN = (h1_systDOWN!=0) ? h1_systDOWN->Integral() : -1.;

  float systUP = (int_systUP-int_mean)/int_mean;
  float systDOWN = (int_systDOWN>0.) ? (int_systDOWN-int_mean)/int_mean : 0.;

  std::pair<float, float> returnSyst;
  returnSyst.first = 1. + systDOWN;
  returnSyst.second = 1. + systUP;

  return returnSyst;

}
