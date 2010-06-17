//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 12 18:22:37 2009 by ROOT version 5.22/00
// from TChain pippo/
//////////////////////////////////////////////////////////

#ifndef TreeAnalyzer_DiJet_h
#define TreeAnalyzer_DiJet_h

#include "TreeAnalyzer.h"

#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cmath>




class TreeAnalyzer_DiJet : public TreeAnalyzer {
public :

   TreeAnalyzer_DiJet( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, TTree* tree=0);

   virtual void Loop();

private:

   Float_t deltaR_phot_;

   Float_t  ptJetReco_;
   Float_t   eJetReco_;
   Float_t  ptCorrJetReco_;
   Float_t   eCorrJetReco_;
   Float_t phiJetReco_;
   Float_t etaJetReco_;

   Float_t    eJetGen_;
   Float_t   ptJetGen_;
   Float_t  phiJetGen_;
   Float_t  etaJetGen_;

   Int_t   pdgIdPart_;
   Float_t    ptPart_;
   Float_t   phiPart_;
   Float_t   etaPart_;

   Int_t   pdgIdPart2nd_;
   Float_t    ptPart2nd_;
   Float_t   phiPart2nd_;
   Float_t   etaPart2nd_;

   Float_t  eTracksReco_;
   Float_t  ePhotonsReco_;
   Float_t  eNeutralHadronsReco_;
   Float_t  eMuonsReco_;
   Float_t  eElectronsReco_;
   Float_t  eHFHadronsReco_;
   Float_t  eHFEMReco_;

   Float_t  ptTracksReco_;
   Float_t  ptPhotonsReco_;
   Float_t  ptNeutralHadronsReco_;
   Float_t  ptMuonsReco_;
   Float_t  ptElectronsReco_;
   Float_t  ptHFHadronsReco_;
   Float_t  ptHFEMReco_;

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
   Float_t  eHFHadronsGen_;
   Float_t  eHFEMGen_;

   Int_t  nTracksGen_;
   Int_t  nPhotonsGen_;
   Int_t  nNeutralHadronsGen_;
   Int_t  nMuonsGen_;
   Int_t  nElectronsGen_;
   Int_t  nHFHadronsGen_;
   Int_t  nHFEMGen_;

   Float_t e2ndJetReco_;
   Float_t pt2ndJetReco_;
   Float_t ptCorr2ndJetReco_;
   Float_t eta2ndJetReco_;
   Float_t phi2ndJetReco_;

   Float_t e2ndJetGen_;
   Float_t pt2ndJetGen_;
   Float_t eta2ndJetGen_;
   Float_t phi2ndJetGen_;

   Float_t  eTracks2ndReco_;
   Float_t  ePhotons2ndReco_;
   Float_t  eNeutralHadrons2ndReco_;
   Float_t  eMuons2ndReco_;
   Float_t  eElectrons2ndReco_;
   Float_t  eHFHadrons2ndReco_;
   Float_t  eHFEM2ndReco_;

   Float_t  ptTracks2ndReco_;
   Float_t  ptPhotons2ndReco_;
   Float_t  ptNeutralHadrons2ndReco_;
   Float_t  ptMuons2ndReco_;
   Float_t  ptElectrons2ndReco_;
   Float_t  ptHFHadrons2ndReco_;
   Float_t  ptHFEM2ndReco_;

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
   Float_t  eHFHadrons2ndGen_;
   Float_t  eHFEM2ndGen_;

   Int_t  nTracks2ndGen_;
   Int_t  nPhotons2ndGen_;
   Int_t  nNeutralHadrons2ndGen_;
   Int_t  nMuons2ndGen_;
   Int_t  nElectrons2ndGen_;
   Int_t  nHFHadrons2ndGen_;
   Int_t  nHFEM2ndGen_;

   Float_t sumEt_pfMet_;
   Float_t pfMet_;

   bool DEBUG_VERBOSE_;
   

};

#endif

