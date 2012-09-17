//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_TTZgen_h
#define Ntp1Analyzer_TTZgen_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_TTZgen : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_TTZgen( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_TTZgen();

   virtual void CreateOutputFile();
   virtual void Loop();




 private:


   Float_t eZllMC_;
   Float_t ptZllMC_;
   Float_t etaZllMC_;
   Float_t phiZllMC_;

   Int_t leptType_;

   Float_t eleptZ1MC_;
   Float_t ptleptZ1MC_;
   Float_t etaleptZ1MC_;
   Float_t phileptZ1MC_;

   Float_t eleptZ2MC_;
   Float_t ptleptZ2MC_;
   Float_t etaleptZ2MC_;
   Float_t phileptZ2MC_;

   Float_t et1_;
   Float_t ptt1_;
   Float_t etat1_;
   Float_t phit1_;

   Float_t et2_;
   Float_t ptt2_;
   Float_t etat2_;
   Float_t phit2_;




   Int_t nJets_;

   Float_t  ptJet_[50];
   Float_t   eJet_[50];
   Float_t phiJet_[50];
   Float_t etaJet_[50];
   Int_t pdgIdJet_[50];


   bool DEBUG_VERBOSE_;

};




#endif
