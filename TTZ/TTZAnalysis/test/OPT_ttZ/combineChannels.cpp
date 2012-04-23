#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"



TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err );
TF1* getLogLikelihoodRatio( const std::string& name, TF2* f2 );
float findMaximum2D( TF2* f2, int nsteps=1000 );



int main() {


  ifstream ifs_SSDL("yields_SSDL.txt");
  ifstream ifs_TTZTL("yields_TTZTL.txt");

  float obs_TTZTL, b_TTZTL, b_err_TTZTL;
  ifs_TTZTL >> obs_TTZTL >> b_TTZTL >> b_err_TTZTL;

  float obs_SSDL, b_SSDL, b_err_SSDL;
  ifs_SSDL >> obs_SSDL >> b_SSDL >> b_err_SSDL;


  TF2* f2_SSDL  = getLikelihoodFunction( "likelihoodSSDL", obs_SSDL, b_SSDL, b_err_SSDL );
  TF2* f2_TTZTL = getLikelihoodFunction( "likelihoodTTZTL", obs_TTZTL, b_TTZTL, b_err_TTZTL );

  //TFile* prova = TFile::Open("prova.root", "recreate");
  //prova->cd();
  //f2_SSDL->Write();
  //f2_TTZTL->Write();
  //prova->Close();

  TF1* f1_llr_SSDL  = getLogLikelihoodRatio( "llrSSDL", f2_SSDL );
  TF1* f1_llr_TTZTL = getLogLikelihoodRatio( "llrTTZTL", f2_TTZTL );

  float ZPL_SSDL = sqrt( -2.*log(f1_llr_SSDL->Eval(0) ) );
  float ZPL_TTZTL = sqrt( -2.*log(f1_llr_TTZTL->Eval(0) ) );

  std::cout << "Z_PL (SSDL): " << ZPL_SSDL << std::endl;
  std::cout << "Z_PL (TTZTL): " << ZPL_TTZTL << std::endl;

  return 0;

}



TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err ) {

  TString name_tstr(name);

  // x is the epxected signal yield, y the expected bg yield
  double xmin = 0.;
  double xmax = (name_tstr.Contains("TTZTL")) ? 10. : 21.;
  double ymin = 0.;
  double ymax = (name_tstr.Contains("TTZTL")) ? 3. : 17.;

  //TF2* f2_new = new TF2( name.c_str(), "(pow(x+y, [0])*exp(-(x+y))/(Factorial([0]))*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2]))");
  TF2* f2_new = new TF2( name.c_str(), "TMath::Poisson( [0], x+y )*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2])", xmin, xmax, ymin, ymax);
  f2_new->SetParameter( 0, obs );
  f2_new->SetParameter( 1, b );
  f2_new->SetParameter( 2, b_err );

  return f2_new;

}



TF1* getLogLikelihoodRatio( const std::string& name, TF2* f2 ) {


  float L_min2d = findMaximum2D( f2 );

  TF1* f1 = new TF1("oo", "1");

  return f1;

}



float findMaximum2D( TF2* f2, int nsteps ) {

  std::cout << "-> Maximixing " << f2->GetName() << std::endl;

  float xmin = f2->GetXmin();
  float xmax = f2->GetXmax();
  float ymin = f2->GetYmin();
  float ymax = f2->GetYmax();

  float xstep = (xmax-xmin)/(float)nsteps;
  float ystep = (ymax-ymin)/(float)nsteps;

  float Lmax_found = 0.;
  float xmax_found = -1.;
  float ymax_found = -1.;

  for( unsigned istepx=0; istepx<nsteps; ++istepx ) {

    for( unsigned istepy=0; istepy<nsteps; ++istepy ) {

      float thisx = istepx*xstep;
      float thisy = istepy*ystep;

      float thisL = f2->Eval( thisx, thisy );

      if( thisL > Lmax_found ) {
        Lmax_found = thisL;
        xmax_found = thisx;
        ymax_found = thisy;
      }

    } // for y

  } // for x

  std::cout << f2->GetName() << " max: " << Lmax_found << " found in (" << xmax_found << "," << ymax_found << ")" << std::endl;

  return Lmax_found;

}
