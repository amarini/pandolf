#ifndef StatTools_h
#define StatTools_h

#include <string>

class TF2;

class StatTools {

 public: 

  StatTools(){};

  static float computeZBi( float obs, float b_pred, float b_pred_err );
  static float computeZPL( float obs, float b_pred, float b_pred_err );

  static TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err );
  static float getLogLikelihoodRatio( const std::string& name, TF2* f2 );
  static float findMaximum2D( TF2* f2, int nsteps, bool fix_x=false );

};


#endif
