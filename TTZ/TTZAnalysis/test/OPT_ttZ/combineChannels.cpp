#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

#include "StatTools.h"




int main() {


  ifstream ifs_SSDL("yields_SSDL.txt");
  ifstream ifs_TTZTL("yields_TTZTL.txt");

  float obs_TTZTL, b_TTZTL, b_err_TTZTL;
  ifs_TTZTL >> obs_TTZTL >> b_TTZTL >> b_err_TTZTL;

  float obs_SSDL, b_SSDL, b_err_SSDL;
  ifs_SSDL >> obs_SSDL >> b_SSDL >> b_err_SSDL;


  float ZBi_SSDL  = StatTools::computeZBi( obs_SSDL, b_SSDL, b_err_SSDL );
  float ZBi_TTZTL = StatTools::computeZBi( obs_TTZTL, b_TTZTL, b_err_TTZTL );

  float ZPL_SSDL  = StatTools::computeZPL( obs_SSDL, b_SSDL, b_err_SSDL );
  float ZPL_TTZTL = StatTools::computeZPL( obs_TTZTL, b_TTZTL, b_err_TTZTL );


  std::cout << std::endl << std::endl;
  std::cout << "SSDL Channel: " << std::endl;
  std::cout << "Observed: " << obs_SSDL << "         Expected Background: " << b_SSDL << " +- " << b_err_SSDL << std::endl;
  std::cout << "Z_PL (SSDL): " << ZPL_SSDL << "   ( ZBi: " << ZBi_SSDL << " )" << std::endl;
  
  std::cout << std::endl << std::endl;
  std::cout << "TTZTL Channel: " << std::endl;
  std::cout << "Observed: " << obs_TTZTL << "         Expected Background: " << b_TTZTL << " +- " << b_err_TTZTL << std::endl;
  std::cout << "Z_PL (TTZTL): " << ZPL_TTZTL << "   ( ZBi: " << ZBi_TTZTL << " )" << std::endl;

  return 0;

}



