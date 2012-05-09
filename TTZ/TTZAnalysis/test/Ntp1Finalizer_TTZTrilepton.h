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

  Ntp1Finalizer_TTZTrilepton( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType="TCHE", const std::string& leptType="ALL");
  virtual ~Ntp1Finalizer_TTZTrilepton() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );

  float get_btagThresh( const std::string& btag_OP_ );

  void set_jes( int jes ) { jes_ = jes; };
  void set_btagSyst( int btagSyst ) { btagSyst_ = btagSyst; };

  void set_jer( bool jer ) { jer_ = jer; };

 private:

   std::string leptType_;
   std::string selectionType_;
   std::string bTaggerType_;

   float  ptJet_thresh_;
   float  etaJet_thresh_;

   int njets_thresh_;
   int nBjets_loose_thresh_;
   int nBjets_medium_thresh_;

   float  ptLeptZ1_thresh_;
   float  ptLeptZ2_thresh_;
   float  ptLept3_thresh_;
   float  etaLeptZ1_thresh_;
   float  etaLeptZ2_thresh_;
   float  etaLept3_thresh_;

   float  combinedIsoRelLept3_thresh_;

   float  met_thresh_;
   float  ht_thresh_;
   float  mbjjZ_best_thresh_;

   float  mZll_threshLo_;
   float  mZll_threshHi_;

   float  ptZll_thresh_;

   // syst flags:
   int jes_; // number of sigmas
   int btagSyst_; // number of sigmas
   bool jer_;

};

