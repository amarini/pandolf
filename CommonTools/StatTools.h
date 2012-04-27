#ifndef StatTools_h
#define StatTools_h

#include <string>

class TF2;


class StatChannel {

 public:

  StatChannel(const std::string& name_t, float obs_t, float b_t, float b_pred_t) {
    name = name_t;
    obs = obs_t;
    b = b_t;
    b_pred = b_pred_t;
  }


  std::string name;
  float obs;
  float b;
  float b_pred;

};


class StatTools {

 public: 

  StatTools(){};

  static float computeZBi( float obs, float b_pred, float b_pred_err );
  static float computeZPL( float obs, float b_pred, float b_pred_err );
  static float computeZPL( TF2* f2_likelihood );
  static float computeZPL( const std::vector<StatChannel>& channels );

  static TF2* getLikelihoodFunction( const std::string& name, int obs, float b, float b_err );
  static float getLogLikelihoodRatio( const std::string& name, TF2* f2 );
  static float findMaximum2D( TF2* f2, int nsteps, bool fix_x=false );

};


#endif
