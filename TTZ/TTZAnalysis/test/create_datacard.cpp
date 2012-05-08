#include <iostream>
#include <fstream>
#include <cstdlib>




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

  datacard << "bin         \t1\t1" << std::endl; 
  datacard << "process     \tttZ\tbg" << std::endl; 
  datacard << "process     \t0\t1" << std::endl; 
  datacard << "rate        \t" << s << "\t" << b_pred << std::endl; 

  datacard << std::endl << std::endl;

  datacard << "lumi lnN    \t1.035\t1.035" << std::endl;
  datacard << "bgUncert    \t-    \t" << 1. + b_pred_err/b_pred << std::endl;

  datacard.close();

  std::cout << "-> Created datacard: " << datacardName << std::endl;

  return 0;

}
