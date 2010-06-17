//------------------------------------------------------------------
//
//    Derived TreeAnalyzer class. Inherits from TreeAnalyzer.
//    Reads output of GammaJetAnalyzer, and produces subtrees,
//    to be used in Higgs->ZZ->lljj analyses.
//
//------------------------------------------------------------------


#ifndef TreeAnalyzer_HZZlljj_h
#define TreeAnalyzer_HZZlljj_h

#include "TreeAnalyzer.h"
#include "TH1F.h"


class TreeAnalyzer_HZZlljj : public TreeAnalyzer {

 public:

   TreeAnalyzer_HZZlljj( std::string dataset, std::string recoType, std::string jetAlgo, TTree* tree=0);

   virtual void Loop();



 private:

   Float_t  ptZ_;
   Float_t   eZ_;
   Float_t phiZ_;
   Float_t etaZ_;

   Float_t  ptJetReco_;
   Float_t   eJetReco_;
   Float_t  ptCorrJetReco_;
   Float_t   eCorrJetReco_;
   Float_t phiJetReco_;
   Float_t etaJetReco_;

   Float_t   ptJetGen_;
   Float_t    eJetGen_;
   Float_t  phiJetGen_;
   Float_t  etaJetGen_;

   Float_t  eTracksReco_;
   Float_t  ePhotonsReco_;
   Float_t  eNeutralHadronsReco_;
   Float_t  eMuonsReco_;
   Float_t  eElectronsReco_;
   Float_t  eHFHadronsReco_;
   Float_t  eHFEMReco_;

   Int_t  nTracksReco_;
   Int_t  nPhotonsReco_;
   Int_t  nNeutralHadronsReco_;
   Int_t  nMuonsReco_;
   Int_t  nElectronsReco_;
   Int_t  nHFHadronsReco_;
   Int_t  nHFEMReco_;

   Float_t  eTracksGen_;
   Float_t  ePhotonsGen_;
   Float_t  eNeutralHadronsGen_;
   Float_t  eMuonsGen_;
   Float_t  eElectronsGen_;

   Int_t  nTracksGen_;
   Int_t  nPhotonsGen_;
   Int_t  nNeutralHadronsGen_;
   Int_t  nMuonsGen_;
   Int_t  nElectronsGen_;

   Float_t  pt2ndJetReco_;
   Float_t   e2ndJetReco_;
   Float_t  ptCorr2ndJetReco_;
   Float_t   eCorr2ndJetReco_;
   Float_t phi2ndJetReco_;
   Float_t eta2ndJetReco_;
   Float_t   pt2ndJetGen_;

   Float_t    e2ndJetGen_;
   Float_t  phi2ndJetGen_;
   Float_t  eta2ndJetGen_;

   Float_t  eTracks2ndReco_;
   Float_t  ePhotons2ndReco_;
   Float_t  eNeutralHadrons2ndReco_;
   Float_t  eMuons2ndReco_;
   Float_t  eElectrons2ndReco_;
   Float_t  eHFHadrons2ndReco_;
   Float_t  eHFEM2ndReco_;

   Int_t  nTracks2ndReco_;
   Int_t  nPhotons2ndReco_;
   Int_t  nNeutralHadrons2ndReco_;
   Int_t  nMuons2ndReco_;
   Int_t  nElectrons2ndReco_;
   Int_t  nHFHadrons2ndReco_;
   Int_t  nHFEM2ndReco_;

   Float_t  eTracks2ndGen_;
   Float_t  ePhotons2ndGen_;
   Float_t  eNeutralHadrons2ndGen_;
   Float_t  eMuons2ndGen_;
   Float_t  eElectrons2ndGen_;

   Int_t  nTracks2ndGen_;
   Int_t  nPhotons2ndGen_;
   Int_t  nNeutralHadrons2ndGen_;
   Int_t  nMuons2ndGen_;
   Int_t  nElectrons2ndGen_;

   Float_t epfMet_;
   Float_t phipfMet_;

   bool DEBUG_VERBOSE_;

  
};




#endif
