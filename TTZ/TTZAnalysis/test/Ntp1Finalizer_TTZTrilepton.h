// ------------------------------------------------------------
//  
//    Ntp1Finalizer_TTZTriplepton - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/BTagSFUtil.h"



class Ntp1Finalizer_TTZTriplepton : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_TTZTriplepton( const std::string& dataset, const std::string& selectionType, const std::string& PUType="HR11", const std::string& leptType="ALL");
  virtual ~Ntp1Finalizer_TTZTriplepton() {};

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
   float  mZll_threshLo_;
   float  mZll_threshHi_;
   float  mZjj_threshLo_;
   float  mZjj_threshHi_;
   float  metSignificance_thresh_;
   bool use_looseBTags_;

};

