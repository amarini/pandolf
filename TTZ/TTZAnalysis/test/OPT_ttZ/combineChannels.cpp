#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"



TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err );
TF1* getLogLikelihoodRatio( const std::string& name, TF2* f2 );



int main() {


  ifstream ifs_SSDL("yields_SSDL.txt");
  ifstream ifs_TTZTL("yields_TTZTL.txt");

  float obs_TTZTL, b_TTZTL, b_err_TTZTL;
  ifs_TTZTL >> obs_TTZTL >> b_TTZTL >> b_err_TTZTL;

  float obs_SSDL, b_SSDL, b_err_SSDL;
  ifs_SSDL >> obs_SSDL >> b_SSDL >> b_err_SSDL;


  TF2* f2_SSDL  = getLikelihoodFunction( "likelihoodSSDL", obs_SSDL, b_SSDL, b_err_SSDL );
  TF2* f2_TTZTL = getLikelihoodFunction( "likelihoodTTZTL", obs_TTZTL, b_TTZTL, b_err_TTZTL );

  TFile* prova = TFile::Open("prova.root", "recreate");
  prova->cd();
  f2_SSDL->Write();
  f2_TTZTL->Write();
  prova->Close();

  TF1* f1_llr_SSDL  = getLogLikelihoodRatio( "llrSSDL", f2_SSDL );
  TF1* f1_llr_TTZTL = getLogLikelihoodRatio( "llrTTZTL", f2_SSDL );

  float ZPL_SSDL = sqrt( -2.*log(f1_llr_SSDL->Eval(0) ) );
  float ZPL_TTZTL = sqrt( -2.*log(f1_llr_TTZTL->Eval(0) ) );

  std::cout << "Z_PL (SSDL): " << ZPL_SSDL << std::endl;
  std::cout << "Z_PL (TTZTL): " << ZPL_TTZTL << std::endl;

  return 0;

}



TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err ) {

  //TF2* f2_new = new TF2( name.c_str(), "(pow(x+y, [0])*exp(-(x+y))/(Factorial([0]))*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2]))");
  TF2* f2_new = new TF2( name.c_str(), "TMath::Poisson( [0], x+y )*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2])");
  f2_new->SetParameter( 0, obs );
  f2_new->SetParameter( 1, b );
  f2_new->SetParameter( 2, b_err );

  return f2_new;

}



TF1* getLogLikelihoodRatio( const std::string& name, TF2* f2 ) {


  TF1* f1 = new TF1("oo", "1");

  return f1;

}
