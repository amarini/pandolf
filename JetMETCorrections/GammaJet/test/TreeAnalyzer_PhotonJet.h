//------------------------------------------------------------------
//
//    Derived TreeAnalyzer class. Inherits from TreeAnalyzer.
//    Reads output of GammaJetAnalyzer, and produces subtrees,
//    to be used in Photon+Jet analyses.
//
//------------------------------------------------------------------


#ifndef TreeAnalyzer_PhotonJet_h
#define TreeAnalyzer_PhotonJet_h

#include "TreeAnalyzer.h"
#include "TH1F.h"


class TreeAnalyzer_PhotonJet : public TreeAnalyzer {

 public:

   TreeAnalyzer_PhotonJet( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", bool useGenJets=false, TTree* tree=0);
   virtual ~TreeAnalyzer_PhotonJet();

   virtual void CreateOutputFile();
   virtual void Loop();
   //float getCorrectionFactor( const std::string& fileName, float pt, float eta );
   //float getCorrectedPt( float pt, float eta, bool applyAlsoResidual );



 private:

   Float_t eventWeight_medium_;
   Float_t eventWeight_loose_;

   Bool_t isIsolated_veryloose_;

   Bool_t isIsolated_hcal_loose_;
   Bool_t isIsolated_hcal_medium_;
   Bool_t isIsolated_hcal_tight_;

   Bool_t isIsolated_ecal_loose_;
   Bool_t isIsolated_ecal_medium_;
   Bool_t isIsolated_ecal_tight_;

   Bool_t isIsolated_ptTracks_loose_;
   Bool_t isIsolated_ptTracks_medium_;
   Bool_t isIsolated_ptTracks_tight_;

   Bool_t isIsolated_nTracks_loose_;
   Bool_t isIsolated_nTracks_medium_;
   Bool_t isIsolated_nTracks_tight_;

   Bool_t clusterMajOK_loose_;
   Bool_t clusterMajOK_medium_;
   Bool_t clusterMajOK_tight_;

   Bool_t clusterMinOK_loose_;
   Bool_t clusterMinOK_medium_;
   Bool_t clusterMinOK_tight_;

   Bool_t passedPhotonID_loose_;
   Bool_t passedPhotonID_medium_;
   Bool_t passedPhotonID_tight_;

   Bool_t matchedToMC_;

   Float_t ePhotReco_;
   Float_t ptPhotReco_;
   Float_t etaPhotReco_;
   Float_t phiPhotReco_;
   Float_t hcalIsoPhotReco_;
   Float_t ecalIsoPhotReco_;
   Int_t   nTrkIsoPhotReco_;
   Float_t ptTrkIsoPhotReco_;
   Float_t clusterMajPhotReco_;
   Float_t clusterMinPhotReco_;
   Int_t   hasPixelSeedPhotReco_;

   Float_t ePhotGen_;
   Float_t ptPhotGen_;
   Float_t etaPhotGen_;
   Float_t phiPhotGen_;

   Float_t deltaR_phot_;

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

   Float_t pt2ndJetReco_;
   Float_t ptCorr2ndJetReco_;
   Float_t eta2ndJetReco_;
   Float_t phi2ndJetReco_;

   Float_t pt2ndJetGen_;
   Float_t eta2ndJetGen_;
   Float_t phi2ndJetGen_;

   Float_t ptSecondaryJetsReco_;
   Float_t ptSecondaryJetsGen_;

   // for mpf:
   Float_t epfMet_;
   Float_t epfMetCorr_;
   Float_t phipfMet_;
   Float_t eMet_;
   Float_t phiMet_;
   Float_t etcMet_;
   Float_t phitcMet_;

   bool DEBUG_VERBOSE_;
   bool useGenJets_;

   TH1F* h1_ptPhot;
   TH1F* h1_eff_denom_vs_pt;
   TH1F* h1_eff_num_medium_vs_pt;
   TH1F* h1_eff_num_loose_vs_pt;
  
};




#endif
