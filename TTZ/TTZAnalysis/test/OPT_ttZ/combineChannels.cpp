#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"

#include "ZBiCalculator.h"



TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err );
float getLogLikelihoodRatio( const std::string& name, TF2* f2 );
float findMaximum2D( TF2* f2, int nsteps=1000, bool fix_x=false );



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

  float llr_SSDL  = getLogLikelihoodRatio( "llrSSDL", f2_SSDL );
  float llr_TTZTL = getLogLikelihoodRatio( "llrTTZTL", f2_TTZTL );

  float ZPL_SSDL = sqrt( -2.*log(llr_SSDL ) );
  float ZPL_TTZTL = sqrt( -2.*log(llr_TTZTL ) );

  ZBiCalculator ZbiCalc;

  std::cout << "Z_PL (SSDL): " << ZPL_SSDL << "   ( ZBi: " << ZbiCalc.computeZBi(obs_SSDL,b_SSDL,b_err_SSDL) << " )" << std::endl;
  std::cout << "Z_PL (TTZTL): " << ZPL_TTZTL << "   ( ZBi: " << ZbiCalc.computeZBi(obs_TTZTL,b_TTZTL,b_err_TTZTL) << " )" << std::endl;

  return 0;

}



TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err ) {

  TString name_tstr(name);

  // x is the epxected signal yield, y the expected bg yield
  double xmin = 0.;
  double xmax = fabs(obs-b) + 5.;
  double ymin = 0.;
  double ymax = 2.*b;

  //TF2* f2_new = new TF2( name.c_str(), "(pow(x+y, [0])*exp(-(x+y))/(Factorial([0]))*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2]))");
  TF2* f2_new = new TF2( name.c_str(), "TMath::Poisson( [0], x+y )*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2])", xmin, xmax, ymin, ymax);
  f2_new->SetParameter( 0, obs );
  f2_new->SetParameter( 1, b );
  f2_new->SetParameter( 2, b_err );

  return f2_new;

}



float getLogLikelihoodRatio( const std::string& name, TF2* f2 ) {


  int nsteps = 1000;
  float L_max2d = findMaximum2D( f2, nsteps );
  float L_max1d_x0 = findMaximum2D( f2, nsteps, true );

  return L_max1d_x0/L_max2d;

}



// finds max of likelihood function, scanning the full x-y phase space
// remember: x is the epxected signal yield, y the expected bg yield
// if fix_x is true, it will maximise only on y, with x=0
float findMaximum2D( TF2* f2, int nsteps, bool fix_x ) {

  if( fix_x ) 
    std::cout << "-> Maximixing " << f2->GetName() << " with x=0." << std::endl;
  else
    std::cout << "-> Maximixing (2D) " << f2->GetName() << std::endl;

  float xmin = f2->GetXmin();
  float xmax = f2->GetXmax();
  float ymin = f2->GetYmin();
  float ymax = f2->GetYmax();

  float xstep = (xmax-xmin)/(float)nsteps;
  float ystep = (ymax-ymin)/(float)nsteps;

  float Lmax_found = 0.;
  float xmax_found = -1.;
  float ymax_found = -1.;

  int nsteps_x = ( fix_x ) ? 1 : nsteps;

  for( unsigned istepx=0; istepx<nsteps_x; ++istepx ) {

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

  if( Lmax_found==0. || xmax_found < 0. || ymax_found < 0. ) {
    std::cout << "ERROR!!! Didn't find a max for function: " << f2->GetName() << std::endl;
    exit(11111);
  }

  std::cout << f2->GetName() << " max: " << Lmax_found << " found in (" << xmax_found << "," << ymax_found << ")" << std::endl;

  return Lmax_found;

}



/*
// finds max of likelihood function, with respect to y only
// (so result is a TF1, function of x)
// remember: x is the epxected signal yield, y the expected bg yield
TF1* findMaximum1D( TF2* f2, int nsteps ) {

  std::cout << "-> Maximixing (1D) " << f2->GetName() << std::endl;

  std::string newname(f2->GetName());
  newname += "_min1d";

  float ymin = f2->GetYmin();
  float ymax = f2->GetYmax();

  float ystep = (ymax-ymin)/(float)nsteps;

  float Lmax_found = 0.;
  float ymax_found = -1.;

  for( unsigned istepy=0; istepy<nsteps; ++istepy ) {

    float thisx = istepx*xstep;
    float thisy = istepy*ystep;

    float thisL = f2->Eval( thisx, thisy );

    if( thisL > Lmax_found ) {
      Lmax_found = thisL;
      ymax_found = thisy;
    }

  } // for y


  std::cout << f2->GetName() << " max: " << Lmax_found << " found in (" << xmax_found << "," << ymax_found << ")" << std::endl;

  return Lmax_found;

}
*/

