// ------------------------------------------------------------
//  
//    Ntp1Finalizer_TTW - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/interface/BTagSFUtil.h"



class Ntp1Finalizer_TTW : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_TTW( const std::string& dataset, const std::string& selectionType, const std::string& jetChoice, const std::string& bTaggerType="JP");
  virtual ~Ntp1Finalizer_TTW() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );

  bool passedBTag( float btag1, float btag2, const std::string& btagger );



 private:

   std::string selectionType_;
   std::string bTaggerType_;
   std::string leptType_;

   std::string jetChoice_;

   float mll_thresh_;

   float ht_thresh_;
   float pfMet_thresh_;

   int njets_thresh_;
   int nbjets_thresh_;
   int nbjetsmed_thresh_;

   float ptLept1_thresh_;
   float ptLept2_thresh_;
   float etaLept1_thresh_;
   float etaLept2_thresh_;

   float ptJet_thresh_;
   float etaJet_thresh_;

};

