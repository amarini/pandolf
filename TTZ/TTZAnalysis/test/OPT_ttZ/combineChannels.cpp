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

  float obs_SSDL_ee, b_SSDL_ee, b_err_SSDL_ee;
  float obs_SSDL_em, b_SSDL_em, b_err_SSDL_em;
  float obs_SSDL_mm, b_SSDL_mm, b_err_SSDL_mm;
  float obs_SSDL_tot, b_SSDL_tot, b_err_SSDL_tot;
  std::string name;
  ifs_SSDL >> name >> obs_SSDL_mm >> b_SSDL_mm >> b_err_SSDL_mm;
  ifs_SSDL >> name >> obs_SSDL_em >> b_SSDL_em >> b_err_SSDL_em;
  ifs_SSDL >> name >> obs_SSDL_ee >> b_SSDL_ee >> b_err_SSDL_ee;
  ifs_SSDL >> name >> obs_SSDL_tot >> b_SSDL_tot >> b_err_SSDL_tot;

  StatChannel channel_SSDL_mm("mm", obs_SSDL_mm, b_SSDL_mm, b_err_SSDL_mm);
  StatChannel channel_SSDL_em("em", obs_SSDL_em, b_SSDL_em, b_err_SSDL_em);
  StatChannel channel_SSDL_ee("ee", obs_SSDL_ee, b_SSDL_ee, b_err_SSDL_ee);

  std::vector< StatChannel > channels_SSDL;
  channels_SSDL.push_back(channel_SSDL_mm);
  channels_SSDL.push_back(channel_SSDL_em);
  channels_SSDL.push_back(channel_SSDL_ee);
  

  float ZBi_SSDL  = StatTools::computeZBi( obs_SSDL_tot, b_SSDL_tot, b_err_SSDL_tot );
  float ZBi_TTZTL = StatTools::computeZBi( obs_TTZTL, b_TTZTL, b_err_TTZTL );

  float ZPL_SSDL_ee  = StatTools::computeZPL( obs_SSDL_ee, b_SSDL_ee, b_err_SSDL_ee );
  float ZPL_SSDL_em  = StatTools::computeZPL( obs_SSDL_em, b_SSDL_em, b_err_SSDL_em );
  float ZPL_SSDL_mm  = StatTools::computeZPL( obs_SSDL_mm, b_SSDL_mm, b_err_SSDL_mm );
  float ZPL_SSDL_tot  = StatTools::computeZPL( obs_SSDL_tot, b_SSDL_tot, b_err_SSDL_tot );
  float ZPL_SSDL_comb  = StatTools::computeZPL( channels_SSDL );

  std::cout <<  "obs_SSDL_ee  : " <<  obs_SSDL_ee <<  "  \tb_SSDL_ee  : " <<  b_SSDL_ee <<  " +- " <<  b_err_SSDL_ee <<  "   \tZPL_SSDL_ee  : " <<  ZPL_SSDL_ee << std::endl;
  std::cout <<  "obs_SSDL_em  : " <<  obs_SSDL_em <<  "  \tb_SSDL_em  : " <<  b_SSDL_em <<  " +- " <<  b_err_SSDL_em <<  "   \tZPL_SSDL_em  : " <<  ZPL_SSDL_em << std::endl;
  std::cout <<  "obs_SSDL_mm  : " <<  obs_SSDL_mm <<  "  \tb_SSDL_mm  : " <<  b_SSDL_mm <<  " +- " <<  b_err_SSDL_mm <<  "   \tZPL_SSDL_mm  : " <<  ZPL_SSDL_mm << std::endl;
  std::cout <<  "obs_SSDL_tot : " <<  obs_SSDL_tot<<  "  \tb_SSDL_tot : " <<  b_SSDL_tot<<  " +- " <<  b_err_SSDL_tot<<  "   \tZPL_SSDL_tot : " <<  ZPL_SSDL_tot<< std::endl;
  std::cout <<  "ZPL_SSDL_comb: " <<  ZPL_SSDL_comb<< std::endl;

  float ZPL_TTZTL = StatTools::computeZPL( obs_TTZTL, b_TTZTL, b_err_TTZTL );

  TF2* f2_lik_SSDL = StatTools::getLikelihoodFunction( "lik_SSDL", obs_SSDL_tot, b_SSDL_tot, b_err_SSDL_tot );
  TF2* f2_lik_TTZTL = StatTools::getLikelihoodFunction( "lik_TTZTL", obs_TTZTL, b_TTZTL, b_err_TTZTL );

  TF2* f2_lik_combined = new TF2("lik_combined", "lik_SSDL*lik_TTZTL", 0., obs_SSDL_tot+obs_TTZTL-b_SSDL_tot-b_TTZTL+5., 0., 2.*(b_TTZTL+b_SSDL_tot) );


  std::cout << std::endl << std::endl;
  std::cout << "SSDL Channel: " << std::endl;
  std::cout << "Observed: " << obs_SSDL_tot << "         Expected Background: " << b_SSDL_tot << " +- " << b_err_SSDL_tot << std::endl;
  std::cout << "Z_PL (SSDL): " << ZPL_SSDL_tot << "   ( ZBi: " << ZBi_SSDL << " )" << std::endl;
  
  std::cout << std::endl << std::endl;
  std::cout << "TTZTL Channel: " << std::endl;
  std::cout << "Observed: " << obs_TTZTL << "         Expected Background: " << b_TTZTL << " +- " << b_err_TTZTL << std::endl;
  std::cout << "Z_PL (TTZTL): " << ZPL_TTZTL << "   ( ZBi: " << ZBi_TTZTL << " )" << std::endl;

  std::cout << "Z_PL (combined): " << StatTools::computeZPL(f2_lik_combined) << std::endl;

  return 0;

}



