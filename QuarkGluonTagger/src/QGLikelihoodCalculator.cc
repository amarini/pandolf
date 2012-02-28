#include "../interface/QGLikelihoodCalculator.h"




QGLikelihoodCalculator::QGLikelihoodCalculator( const std::string& filedirPath ) {

  std::string fileName_nCharged = filedirPath + "/nCharged.txt";
  std::string fileName_nNeutral = filedirPath + "/nNeutral.txt";
  std::string fileName_ptD = filedirPath + "/ptD.txt";

  jcp_nCharged_quark_ = new JetCorrectorParameters(fileName_nCharged, "quark");
  jcp_nCharged_gluon_ = new JetCorrectorParameters(fileName_nCharged, "gluon");

  jcp_nNeutral_quark_ = new JetCorrectorParameters(fileName_nNeutral, "quark");
  jcp_nNeutral_gluon_ = new JetCorrectorParameters(fileName_nNeutral, "gluon");

  jcp_ptD_quark_ = new JetCorrectorParameters(fileName_ptD, "quark");
  jcp_ptD_gluon_ = new JetCorrectorParameters(fileName_ptD, "gluon");


  sjc_nCharged_quark_ = new SimpleJetCorrector(*jcp_nCharged_quark_);
  sjc_nCharged_gluon_ = new SimpleJetCorrector(*jcp_nCharged_gluon_);

  sjc_nNeutral_quark_ = new SimpleJetCorrector(*jcp_nNeutral_quark_);
  sjc_nNeutral_gluon_ = new SimpleJetCorrector(*jcp_nNeutral_gluon_);

  sjc_ptD_quark_ = new SimpleJetCorrector(*jcp_ptD_quark_);
  sjc_ptD_gluon_ = new SimpleJetCorrector(*jcp_ptD_gluon_);

}




QGLikelihoodCalculator::~QGLikelihoodCalculator() {

  delete jcp_nCharged_quark_;
  delete jcp_nCharged_gluon_;
  delete jcp_nNeutral_quark_;
  delete jcp_nNeutral_gluon_;
  delete jcp_ptD_quark_;
  delete jcp_ptD_gluon_;

  delete sjc_nCharged_quark_;
  delete sjc_nCharged_gluon_;
  delete sjc_nNeutral_quark_;
  delete sjc_nNeutral_gluon_;
  delete sjc_ptD_quark_;
  delete sjc_ptD_gluon_;

}




float QGLikelihoodCalculator::computeQGLikelihood( float pt, float rhoPF, int nCharged, int nNeutral, float ptD ) {


  std::vector<float> v_pt_rho;
  v_pt_rho.push_back( pt );
  v_pt_rho.push_back( rhoPF );

  std::vector<float> v_nCharged;
  v_nCharged.push_back( (float)nCharged );

  std::vector<float> v_nNeutral;
  v_nNeutral.push_back( (float)nNeutral );

  std::vector<float> v_ptD;
  v_ptD.push_back( ptD );

  float quarkProb_nCharged = sjc_nCharged_quark_->correction(v_pt_rho, v_nCharged);
  float gluonProb_nCharged = sjc_nCharged_gluon_->correction(v_pt_rho, v_nCharged);

  float quarkProb_nNeutral = sjc_nNeutral_quark_->correction(v_pt_rho, v_nNeutral);
  float gluonProb_nNeutral = sjc_nNeutral_gluon_->correction(v_pt_rho, v_nNeutral);

  float quarkProb_ptD = sjc_ptD_quark_->correction(v_pt_rho, v_ptD);
  float gluonProb_ptD = sjc_ptD_gluon_->correction(v_pt_rho, v_ptD);


  float quarkProb = quarkProb_nCharged*quarkProb_nNeutral*quarkProb_ptD;
  float gluonProb = gluonProb_nCharged*gluonProb_nNeutral*gluonProb_ptD;

  float QGLikelihood = (gluonProb+quarkProb>0.) ? quarkProb / (gluonProb + quarkProb ) : -1.;


  return QGLikelihood;

}



