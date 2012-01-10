//------------------------------------------------------------------
//
//    Derived TreeAnalyzer class. Inherits from TreeAnalyzer.
//    Reads output of GammaJetAnalyzer, and produces subtrees,
//    to be used in Photon+Jet analyses.
//
//------------------------------------------------------------------


#ifndef TreeAnalyzer_VHgg_h
#define TreeAnalyzer_VHgg_h

#include "TreeAnalyzer.h"
#include "TH1F.h"


class TreeAnalyzer_VHgg : public TreeAnalyzer {

 public:

   TreeAnalyzer_VHgg( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", bool useGenJets=false, TTree* tree=0);
   virtual ~TreeAnalyzer_VHgg();

   virtual void CreateOutputFile();
   virtual void Loop();
 //float getCorrectionFactor( const std::string& fileName, float pt, float eta );
 //float getCorrectedPt( float pt, float eta, bool applyAlsoResidual );



 private:


   Float_t eventWeight_medium_;
   Float_t eventWeight_loose_;

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
   Bool_t passed_Photon75_;

   Float_t ePhot1_;
   Float_t ptPhot1_;
   Float_t etaPhot1_;
   Float_t phiPhot1_;
   Float_t hcalIsoPhot1_;
   Float_t ecalIsoPhot1_;
   Int_t   nTrkIsoPhot1_;
   Float_t ptTrkIsoPhot1_;
   Float_t clusterMajPhot1_;
   Float_t clusterMinPhot1_;
   Int_t   hasPixelSeedPhot1_;
   Float_t pid_twrHCALPhot1_;
   Float_t pid_HoverEPhot1_;
   Float_t pid_jurECALPhot1_;
   Float_t pid_sIEtaIEtaPhot1_;
   Float_t pid_hlwTrackPhot1_;
   Float_t pid_hlwTrackNoDzPhot1_;

   Float_t ePhotGen1_;
   Float_t ptPhotGen1_;
   Float_t etaPhotGen1_;
   Float_t phiPhotGen1_;

   Float_t ePhot2_;
   Float_t ptPhot2_;
   Float_t etaPhot2_;
   Float_t phiPhot2_;
   Float_t hcalIsoPhot2_;
   Float_t ecalIsoPhot2_;
   Int_t   nTrkIsoPhot2_;
   Float_t ptTrkIsoPhot2_;
   Float_t clusterMajPhot2_;
   Float_t clusterMinPhot2_;
   Int_t   hasPixelSeedPhot2_;
   Float_t pid_twrHCALPhot2_;
   Float_t pid_HoverEPhot2_;
   Float_t pid_jurECALPhot2_;
   Float_t pid_sIEtaIEtaPhot2_;
   Float_t pid_hlwTrackPhot2_;
   Float_t pid_hlwTrackNoDzPhot2_;

   Float_t ePhotGen2_;
   Float_t ptPhotGen2_;
   Float_t etaPhotGen2_;
   Float_t phiPhotGen2_;


   Int_t nJet_;
   Float_t  ptJet_[20];
   Float_t   eJet_[20];
   Float_t  ptRawJet_[20];
   Float_t   eRawJet_[20];
   Float_t phiJet_[20];
   Float_t etaJet_[20];
   Float_t  ptDJet_[20];
   Float_t  rmsCandJet_[20];
   Float_t  trackCountingHighEffBJetTagsJet_[20];
   Float_t  simpleSecondaryVertexHighEffBJetTagsJet_[20];
   Float_t   ptJetGen_[20];

   Float_t    eJetGen_[20];
   Float_t  phiJetGen_[20];
   Float_t  etaJetGen_[20];

   Int_t   pdgIdPartJet_[20];
   Float_t    ptPartJet_[20];
   Float_t   phiPartJet_[20];
   Float_t   etaPartJet_[20];
   Float_t   ePartJet_[20];

   Float_t  eTracksJet_[20];
   Float_t  ePhotonsJet_[20];
   Float_t  eNeutralHadronsJet_[20];
   Float_t  eMuonsJet_[20];
   Float_t  eElectronsJet_[20];
   Float_t  eHFHadronsJet_[20];
   Float_t  eHFEMJet_[20];

   Int_t  nTracksJet_[20];
   Int_t  nPhotonsJet_[20];
   Int_t  nNeutralHadronsJet_[20];
   Int_t  nMuonsJet_[20];
   Int_t  nElectronsJet_[20];
   Int_t  nHFHadronsJet_[20];
   Int_t  nHFEMJet_[20];


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
  
};




#endif
