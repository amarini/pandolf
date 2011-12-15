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

   float  ptLeptZ1_thresh_;
   float  ptLeptZ2_thresh_;
   float  ptLept3_thresh_;
   float  etaLeptZ1_thresh_;
   float  etaLeptZ2_thresh_;
   float  etaLept3_thresh_;
   float  ptJetB1_thresh_;
   float  ptJetB2_thresh_;
   float  ptJet3_thresh_;
   float  ptJet4_thresh_;
   float  etaJetB1_thresh_;
   float  etaJetB2_thresh_;
   float  etaJet3_thresh_;
   float  etaJet4_thresh_;
   float  mZll_threshLo_;
   float  mZll_threshHi_;

};

