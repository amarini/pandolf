// ------------------------------------------------------------
//  
//    Ntp1Finalizer_TTZTrilepton - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/BTagSFUtil.h"



class Ntp1Finalizer_TTZTrilepton : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_TTZTrilepton( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType="SSVHE", const std::string& PUType="HR11", const std::string& leptType="ALL");
  virtual ~Ntp1Finalizer_TTZTrilepton() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );

  float get_btagThresh( const std::string& btag_OP_ );


 private:

   std::string leptType_;
   std::string PUType_;
   std::string selectionType_;
   std::string bTaggerType_;

   float  ptJet_thresh_;
   float  ptBJet_thresh_;
   float  etaJet_thresh_;

   std::string btagJetB1_OP_;
   std::string btagJetB2_OP_;

   float  ptLeptZ1_thresh_;
   float  ptLeptZ2_thresh_;
   float  ptLept3_thresh_;
   float  etaLeptZ1_thresh_;
   float  etaLeptZ2_thresh_;
   float  etaLept3_thresh_;

   float  combinedIsoRelLept3_thresh_;

   float  met_thresh_;
   float  ht_thresh_;

   int    njets_thresh_;

   float  mZll_threshLo_;
   float  mZll_threshHi_;

   float  ptZll_thresh_;

};

