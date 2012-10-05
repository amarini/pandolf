//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_LeptonStudies_h
#define Ntp1Analyzer_LeptonStudies_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_LeptonStudies : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_LeptonStudies( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_LeptonStudies();

   virtual void CreateOutputFile();
   virtual void Loop();




 private:

   int   nEle_;
   float ptEle_[10];
   float etaEle_[10];
   float pfIsoEle_[10];
   bool  isLooseSUSYElectronEle_[10];
   bool  isTightSUSYElectronEle_[10];
   bool  isNotConversionEle_[10];
   bool  passedHLTEle_[10];
   float mvaidtrigEle_[10];
   float mvaidnontrigEle_[10];
   bool  matchedToGenEle_[10];



   bool DEBUG_VERBOSE_;

};




#endif
