//------------------------------------------------------------------
//
//    Derived TreeAnalyzer class. Inherits from TreeAnalyzer.
//    Reads output of GammaJetAnalyzer, and produces subtrees,
//    to be used in Photon+Jet analyses.
//
//------------------------------------------------------------------


#ifndef TreeAnalyzer_PhotonJetInclusive_h
#define TreeAnalyzer_PhotonJetInclusive_h

#include "TreeAnalyzer.h"
#include "TH1F.h"


class TreeAnalyzer_PhotonJetInclusive : public TreeAnalyzer {

 public:

   TreeAnalyzer_PhotonJetInclusive( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", bool useGenJets=false, TTree* tree=0);
   virtual ~TreeAnalyzer_PhotonJetInclusive();

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

   Bool_t passed_Photon10_;
   Bool_t passed_Photon15_;
   Bool_t passed_Photon20_;
   Bool_t passed_Photon25_;
   Bool_t passed_Photon30_;
   Bool_t passed_Photon35_;
   Bool_t passed_Photon40_;
   Bool_t passed_Photon50_;
   Bool_t passed_Photon60_;
   Bool_t passed_Photon70_;

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
   Float_t pid_twrHCALPhotReco_;
   Float_t pid_HoverEPhotReco_;
   Float_t pid_jurECALPhotReco_;
   Float_t pid_sIEtaIEtaPhotReco_;
   Float_t pid_hlwTrackPhotReco_;
   Float_t pid_hlwTrackNoDzPhotReco_;

   Float_t ePhotGen_;
   Float_t ptPhotGen_;
   Float_t etaPhotGen_;
   Float_t phiPhotGen_;

   Float_t deltaR_phot_;

   Int_t  nJet_;
   Float_t  ptJet_[100];
   Float_t   eJet_[100];
   Float_t  ptCorrJet_[100];
   Float_t   eCorrJet_[100];
   Float_t phiJet_[100];
   Float_t etaJet_[100];
   Float_t  ptDJet_[100];
   Float_t  rmsCandJet_[100];
   Float_t  trackCountingHighEffBJetTagsJet_[100];

   Float_t   ptJetGen_[100];
   Float_t    eJetGen_[100];
   Float_t  phiJetGen_[100];
   Float_t  etaJetGen_[100];

   Int_t   pdgIdPart_[100];
   Float_t    ptPart_[100];
   Float_t   phiPart_[100];
   Float_t   etaPart_[100];

   Float_t eChargedHadronsJet_[100];
   Float_t ePhotonsJet_[100];
   Float_t eNeutralHadronsJet_[100];
   Float_t eMuonsJet_[100];
   Float_t eElectronsJet_[100];
   Float_t eHFHadronsJet_[100];
   Float_t eHFEMJet_[100];

   Int_t nChargedHadronsJet_[100];
   Int_t nPhotonsJet_[100];
   Int_t nNeutralHadronsJet_[100];
   Int_t nMuonsJet_[100];
   Int_t nElectronsJet_[100];
   Int_t nHFHadronsJet_[100];
   Int_t nHFEMJet_[100];



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
