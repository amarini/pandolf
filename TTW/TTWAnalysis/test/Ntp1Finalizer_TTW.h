// ------------------------------------------------------------
//  
//    Ntp1Finalizer_TTW - Derived class 
//    for the finalization of the H->ZZ->lljj analysis.
//
// ------------------------------------------------------------



#include "Ntp1Finalizer.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/BTagSFUtil.h"



class Ntp1Finalizer_TTW : public Ntp1Finalizer {

 public:

  Ntp1Finalizer_TTW( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType="TCHE", const std::string& leptType="ALL");
  virtual ~Ntp1Finalizer_TTW() {};

  virtual void finalize();
  void setSelectionType( const std::string& selectionType );

  float get_helicityLD_thresh(float mass, int nBTags);


 private:

   std::string selectionType_;
   std::string bTaggerType_;
   std::string leptType_;

   float  ptLept1_thresh_;
   float  ptLept2_thresh_;
   float  etaLept1_thresh_;
   float  etaLept2_thresh_;
   float  ptJet_thresh_;
   float  ptJet1_thresh_;
   float  ptJet2_thresh_;
   float  ptJet3_thresh_;
   float  ptJet4_thresh_;
   float  etaJet1_thresh_;
   float  etaJet2_thresh_;
   float  etaJet3_thresh_;
   float  etaJet4_thresh_;

};

