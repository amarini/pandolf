#include "StatTools.h"

#include "TMath.h"
#include "TF2.h"


float StatTools::computeZBi( float obs, float b_pred, float b_pred_err ) {

  float tau = b_pred / ( b_pred_err*b_pred_err );
  float n_off = tau*b_pred;
  float P_Bi = TMath::BetaIncomplete( 1./(1.+tau), obs, n_off+1. );
  float Z_Bi = sqrt(2)*TMath::ErfInverse( 1 - 2.*P_Bi );

  return Z_Bi;

}



float StatTools::computeZPL( float obs, float b_pred, float b_pred_err ) {

  TF2* f2_likelihood = getLikelihoodFunction( "likelihood", obs, b_pred, b_pred_err );

  float ZPL = computeZPL(f2_likelihood);

  delete f2_likelihood;

  return ZPL;

}



float StatTools::computeZPL( TF2* f2_likelihood ) {

  float llr  = getLogLikelihoodRatio( "llr", f2_likelihood );

  float ZPL = sqrt( -2.*log(llr ) );

  return ZPL;

}



float StatTools::computeZPL( const std::vector< StatChannel >& channels ) {

  std::vector<TF2*> functions;
  std::string combined_name;

  float xmax, ymax;

  for( unsigned ichan=0; ichan<channels.size(); ++ichan ) {

    std::string f_name = "likelihood_" + channels[ichan].name;
    TF2* f2_likelihood = getLikelihoodFunction( f_name, channels[ichan].obs, channels[ichan].b, channels[ichan].b_err );

    functions.push_back(f2_likelihood);

    if( ichan==0 )
      combined_name = f_name;
    else
      combined_name = combined_name + "*" + f_name;

    xmax += fabs(channels[ichan].obs-channels[ichan].b);
    ymax += 2.*channels[ichan].b;

  }

  
  TF2* combined_likelihood = new TF2( "combined_likelihood", combined_name.c_str(), 0., xmax, 0., ymax );
  
  return computeZPL( combined_likelihood );

}



TF2* StatTools::getLikelihoodFunction( const std::string& name, int obs, float b, float b_err ) {

  TString name_tstr(name);

  // x is the epxected signal yield, y the expected bg yield
  double xmin = 0.;
  double xmax = fabs(obs-b) + 5.;
  double ymin = 0.;
  double ymax = 2.*b;

  TF2* f2_new = new TF2( name.c_str(), "TMath::Poisson( [0], x+y )*exp(-0.5*((y-[1])/[2])**2)/(sqrt(2*pi)*[2])", xmin, xmax, ymin, ymax);
  f2_new->SetParameter( 0, obs );
  f2_new->SetParameter( 1, b );
  f2_new->SetParameter( 2, b_err );

  return f2_new;

}



float StatTools::getLogLikelihoodRatio( const std::string& name, TF2* f2, int nsteps ) {

  float L_max2d = findMaximum2D( f2, nsteps );
  float L_max1d_x0 = findMaximum2D( f2, nsteps, true );

  return L_max1d_x0/L_max2d;

}



// finds max of likelihood function, scanning the full x-y phase space
// remember: x is the epxected signal yield, y the expected bg yield
// if fix_x is true, it will maximise only on y, with x=0
float StatTools::findMaximum2D( TF2* f2, int nsteps, bool fix_x ) {

  if( fix_x ) 
    std::cout << "[StatTools]:: Maximixing " << f2->GetName() << " with x=0." << std::endl;
  else
    std::cout << "[StatTools]:: Maximixing (2D) " << f2->GetName() << std::endl;

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
    std::cout << "[StatTools]:: ERROR!!! Didn't find a max for function: " << f2->GetName() << std::endl;
    exit(11111);
  }

  //std::cout << f2->GetName() << " max: " << Lmax_found << " found in (" << xmax_found << "," << ymax_found << ")" << std::endl;

  return Lmax_found;

}




