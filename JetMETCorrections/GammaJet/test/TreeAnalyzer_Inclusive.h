//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 12 18:22:37 2009 by ROOT version 5.22/00
// from TChain pippo/
//////////////////////////////////////////////////////////

#ifndef TreeAnalyzer_Inclusive_h
#define TreeAnalyzer_Inclusive_h

#include "TreeAnalyzer.h"

#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cmath>





class TreeAnalyzer_Inclusive : public TreeAnalyzer {
public :

   TreeAnalyzer_Inclusive( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", TTree* tree=0);

   virtual void Loop();

private :

   Int_t    nJetReco_;
   Float_t  ptJetReco_[40];
   Float_t   eJetReco_[40];
   Float_t  ptCorrJetReco_[40];
   Float_t   eCorrJetReco_[40];
   Float_t phiJetReco_[40];
   Float_t etaJetReco_[40];
   Float_t   ptJetGen_[40];

   Float_t    eJetGen_[40];
   Float_t  phiJetGen_[40];
   Float_t  etaJetGen_[40];

// Int_t   pdgIdPart_;
// Float_t    ptPart_;
// Float_t   phiPart_;
// Float_t   etaPart_;

   Float_t  emfJetReco_[40];

   Float_t  eTracksReco_[40];
   Float_t  ePhotonsReco_[40];
   Float_t  eNeutralHadronsReco_[40];
   Float_t  eMuonsReco_[40];
   Float_t  eElectronsReco_[40];
   Float_t  eHFHadronsReco_[40];
   Float_t  eHFEMReco_[40];

   Float_t  ptTracksReco_[40];
   Float_t  ptPhotonsReco_[40];
   Float_t  ptNeutralHadronsReco_[40];
   Float_t  ptMuonsReco_[40];
   Float_t  ptElectronsReco_[40];
   Float_t  ptHFHadronsReco_[40];
   Float_t  ptHFEMReco_[40];

   Int_t  nTracksReco_[40];
   Int_t  nPhotonsReco_[40];
   Int_t  nNeutralHadronsReco_[40];
   Int_t  nMuonsReco_[40];
   Int_t  nElectronsReco_[40];
   Int_t  nHFHadronsReco_[40];
   Int_t  nHFEMReco_[40];

   Float_t  eTracksGen_[40];
   Float_t  ePhotonsGen_[40];
   Float_t  eNeutralHadronsGen_[40];
   Float_t  eMuonsGen_[40];
   Float_t  eElectronsGen_[40];
   Float_t  eHFHadronsGen_[40];
   Float_t  eHFEMGen_[40];

   Int_t  nTracksGen_[40];
   Int_t  nPhotonsGen_[40];
   Int_t  nNeutralHadronsGen_[40];
   Int_t  nMuonsGen_[40];
   Int_t  nElectronsGen_[40];
   Int_t  nHFHadronsGen_[40];
   Int_t  nHFEMGen_[40];

   Float_t sumEt_pfMet_;
   Float_t pfMet_;

   bool DEBUG_VERBOSE_;


};

#endif // #ifdef TreeAnalyzer_Inclusive_cxx
