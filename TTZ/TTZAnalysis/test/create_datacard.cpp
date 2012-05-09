#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"


std::pair< float, float >  getBGSyst( const std::string& syst, const std::string& sel );


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

  datacard << "bin         \t1      \t1" << std::endl; 
  datacard << "process     \tttZ    \tbg" << std::endl; 
  datacard << "process     \t0      \t1" << std::endl; 
  datacard << "rate        \t" << s << "\t" << b_pred << std::endl; 

  datacard << std::endl << std::endl;

  datacard << "lumi     lnN\t1.035  \t1.035" << std::endl;
  datacard << "bgUncert lnN\t-      \t" << 1. + b_pred_err/b_pred << std::endl;

  std::pair< float, float >  btagSyst = getBGSyst( "BTag", selection );
  datacard << "btag     lnN\t-      \t" << btagSyst.first << "/" << btagSyst.second << std::endl;

  std::pair< float, float >  jesSyst = getBGSyst( "JES", selection );
  datacard << "jes      lnN\t-      \t" << jesSyst.first << "/" << jesSyst.second << std::endl;

  datacard.close();

  std::cout << "-> Created datacard: " << datacardName << std::endl;

  return 0;

}




std::pair<float, float>  getBGSyst( const std::string& syst, const std::string& sel ) {

  std::string systFile = "TTZTrilepton_BG_" + sel + "_TCHE_ALL.root";
  std::string systFileUP = "TTZTrilepton_BG_" + sel + "_TCHE_ALL_" + syst + "UP.root";
  std::string systFileDOWN = "TTZTrilepton_BG_" + sel + "_TCHE_ALL_" + syst + "DOWN.root";

  TFile* file_systFile = TFile::Open( systFile.c_str() );
  TFile* file_systFileUP = TFile::Open( systFileUP.c_str() );
  TFile* file_systFileDOWN = TFile::Open( systFileDOWN.c_str() );

  TH1D* h1_mean = (TH1D*)file_systFile->Get("channelYields");
  TH1D* h1_systUP = (TH1D*)file_systFileUP->Get("channelYields");
  TH1D* h1_systDOWN = (TH1D*)file_systFileDOWN->Get("channelYields");

  float int_mean = h1_mean->Integral();
  float int_systUP = h1_systUP->Integral();
  float int_systDOWN = h1_systDOWN->Integral();

  float systUP = (int_systUP-int_mean)/int_mean;
  float systDOWN = (int_systDOWN-int_mean)/int_mean;

  std::cout << syst << " Syst: UP: " << systUP << " DOWN: " << systDOWN << std::endl;
//  float systValue = ( fabs(systUP)>fabs(systDOWN) ) ? systUP : systDOWN;
//  //float systValue = TMath::Max(systUP, systDOWN);
  

  std::pair<float, float> returnSyst;
  returnSyst.first = 1. + systDOWN;
  returnSyst.second = 1. + systUP;

  return returnSyst;

}
