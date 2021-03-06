// ------------------------------------------------------------
//  
//    Ntp1Finalizer_HZZlljjRM - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/BTagSFUtil.h"



class Ntp1Finalizer_HZZlljjRM : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_HZZlljjRM( const std::string& dataset, const std::string& selectionType, const std::string& PUType="HR11", const std::string& leptType="ALL");
  virtual ~Ntp1Finalizer_HZZlljjRM() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );

  float get_helicityLD_thresh(float mass, int nBTags);


 private:

   std::string leptType_;
   std::string PUType_;
   std::string selectionType_;

   float  ptLept1_thresh_;
   float  ptLept2_thresh_;
   float  etaLept1_thresh_;
   float  etaLept2_thresh_;
   float  ptJet1_thresh_;
   float  ptJet2_thresh_;
   float  etaJet1_thresh_;
   float  etaJet2_thresh_;
   bool   invert_mZll_;
   float  mZll_threshLo_;
   float  mZll_threshHi_;
   float  mZjj_threshLo_;
   float  mZjj_threshHi_;
   float  helicityLD_slope_0btags_;
   float  helicityLD_slope_1btags_;
   float  helicityLD_slope_2btags_;
   float  helicityLD_intercept_0btags_;
   float  helicityLD_intercept_1btags_;
   float  helicityLD_intercept_2btags_;
   float  helicityLD_minThresh_0btags_;
   float  helicityLD_minThresh_1btags_;
   float  helicityLD_minThresh_2btags_;
   float  helicityLD_maxThresh_0btags_;
   float  helicityLD_maxThresh_1btags_;
   float  helicityLD_maxThresh_2btags_;
   float  QGLikelihoodProd_thresh_;
   float  metSignificance_thresh_;
   bool use_looseBTags_;

};

