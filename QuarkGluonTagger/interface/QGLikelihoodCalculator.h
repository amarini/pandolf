// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//
// ------------------------------------------------------------

#ifndef QGLikelihoodCalculator_h
#define QGLikelihoodCalculator_h

#include <string>

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"




class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName_nCharged="QGTaggerConfig_nCharged_AK5PF.txt", const std::string& fileName_nNeutral="QGTaggerConfig_nNeutral_AK5PF.txt", const std::string& fileName_ptD="QGTaggerConfig_ptD_AK5PF.txt");
   ~QGLikelihoodCalculator();


  float computeQGLikelihood( float pt, float rhoPF, int nCharged, int nNeutral, float ptD );
  

 private:


  JetCorrectorParameters *jcp_nCharged_quark_;
  JetCorrectorParameters *jcp_nCharged_gluon_;
  JetCorrectorParameters *jcp_nNeutral_quark_;
  JetCorrectorParameters *jcp_nNeutral_gluon_;
  JetCorrectorParameters *jcp_ptD_quark_;
  JetCorrectorParameters *jcp_ptD_gluon_;

  SimpleJetCorrector *sjc_nCharged_quark_;
  SimpleJetCorrector *sjc_nCharged_gluon_;
  SimpleJetCorrector *sjc_nNeutral_quark_;
  SimpleJetCorrector *sjc_nNeutral_gluon_;
  SimpleJetCorrector *sjc_ptD_quark_;
  SimpleJetCorrector *sjc_ptD_gluon_;

};


#endif
