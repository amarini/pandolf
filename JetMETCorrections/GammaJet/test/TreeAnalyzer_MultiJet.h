//------------------------------------------------------------------
//
//    Derived TreeAnalyzer class. Inherits from TreeAnalyzer.
//    Reads output of GammaJetAnalyzer, and produces subtrees,
//    to be used in Photon+Jet analyses.
//
//------------------------------------------------------------------


#ifndef TreeAnalyzer_MultiJet_h
#define TreeAnalyzer_MultiJet_h

#include "TreeAnalyzer.h"
#include "TH1F.h"


class TreeAnalyzer_MultiJet : public TreeAnalyzer {

 public:

   TreeAnalyzer_MultiJet( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", TTree* tree=0);
   virtual ~TreeAnalyzer_MultiJet();

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

   Bool_t matchedToMC_;

   Int_t nJet_;
   Float_t  ptJet_[20];
   Float_t   eJet_[20];
   Float_t  ptRawJet_[20];
   Float_t   eRawJet_[20];
   Float_t phiJet_[20];
   Float_t etaJet_[20];
   Float_t  ptDJet_[20];
   Float_t  rmsCandJet_[20];
   Float_t  QGLikelihoodJet_[20];
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
   Int_t   pdgIdMomJet_[20];

   Int_t   pdgIdPartStatus3Jet_[20];
   Float_t    ptPartStatus3Jet_[20];
   Float_t   phiPartStatus3Jet_[20];
   Float_t   etaPartStatus3Jet_[20];
   Float_t     ePartStatus3Jet_[20];
   Int_t    pdgIdMomStatus3Jet_[20];

   Float_t  eChargedHadronsJet_[20];
   Float_t  ePhotonsJet_[20];
   Float_t  eNeutralHadronsJet_[20];
   Float_t  eMuonsJet_[20];
   Float_t  eElectronsJet_[20];
   Float_t  eHFHadronsJet_[20];
   Float_t  eHFEMJet_[20];

   Int_t  nChargedHadronsJet_[20];
   Int_t  nPhotonsJet_[20];
   Int_t  nNeutralHadronsJet_[20];
   Int_t  nMuonsJet_[20];
   Int_t  nElectronsJet_[20];
   Int_t  nHFHadronsJet_[20];
   Int_t  nHFEMJet_[20];


   Float_t epfMet_;
   Float_t epfMetCorr_;
   Float_t phipfMet_;
   Float_t eMet_;
   Float_t phiMet_;
   Float_t etcMet_;
   Float_t phitcMet_;

   Float_t ht_akt5_;

   Bool_t passed_HT250_;
   Bool_t passed_HT300_;
   Bool_t passed_HT350_;
   Bool_t passed_HT400_;
   Bool_t passed_HT450_;
   Bool_t passed_HT500_;
   Bool_t passed_HT550_;
   Bool_t passed_HT600_;
   Bool_t passed_HT650_;
   Bool_t passed_HT700_;
   Bool_t passed_HT750_;


   bool DEBUG_VERBOSE_;
  
};




#endif
