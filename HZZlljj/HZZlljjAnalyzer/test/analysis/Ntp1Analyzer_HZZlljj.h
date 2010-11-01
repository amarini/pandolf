//------------------------------------------------------------------
//
//    Derived Ntp1Analyzer class. Inherits from Ntp1Analyzer.
//    Reads output of e/c/p tree dumper, and produces subtrees,
//    to be used in the H->ZZ->lljj analysis.
//
//------------------------------------------------------------------


#ifndef Ntp1Analyzer_HZZlljj_h
#define Ntp1Analyzer_HZZlljj_h

#include "Ntp1Analyzer.h"
#include "TH1F.h"


class Ntp1Analyzer_HZZlljj : public Ntp1Analyzer {

 public:

   Ntp1Analyzer_HZZlljj( const std::string& dataset, const std::string& flags="", TTree* tree=0);
   virtual ~Ntp1Analyzer_HZZlljj();

   virtual void CreateOutputFile();
   virtual void Loop();



 private:

   int leptType_; //0: muon; 1: electron

   Float_t eZqqMC_;
   Float_t ptZqqMC_;
   Float_t etaZqqMC_;
   Float_t phiZqqMC_;

   Float_t eZllMC_;
   Float_t ptZllMC_;
   Float_t etaZllMC_;
   Float_t phiZllMC_;

   Float_t eHiggsMC_;
   Float_t ptHiggsMC_;
   Float_t etaHiggsMC_;
   Float_t phiHiggsMC_;

   Float_t eLept1_;
   Float_t ptLept1_;
   Float_t etaLept1_;
   Float_t phiLept1_;

   Float_t eLept1Gen_;
   Float_t ptLept1Gen_;
   Float_t etaLept1Gen_;
   Float_t phiLept1Gen_;

   Float_t eLept2_;
   Float_t ptLept2_;
   Float_t etaLept2_;
   Float_t phiLept2_;

   Float_t eLept2Gen_;
   Float_t ptLept2Gen_;
   Float_t etaLept2Gen_;
   Float_t phiLept2Gen_;

   Float_t  ptJetRecoil_;
   Float_t   eJetRecoil_;
   Float_t phiJetRecoil_;
   Float_t etaJetRecoil_;

   Float_t  ptJetLead_;
   Float_t   eJetLead_;
   Float_t phiJetLead_;
   Float_t etaJetLead_;

   Float_t  ptJetLead2_;
   Float_t   eJetLead2_;
   Float_t phiJetLead2_;
   Float_t etaJetLead2_;

   Float_t  ptJetLead3_;
   Float_t   eJetLead3_;
   Float_t phiJetLead3_;
   Float_t etaJetLead3_;

   Int_t nPairs_;

   Int_t  iJet1_[20];
   Float_t  ptJet1_[20];
   Float_t   eJet1_[20];
   Float_t phiJet1_[20];
   Float_t etaJet1_[20];

   Float_t   ptJetGen1_;
   Float_t    eJetGen1_;
   Float_t  phiJetGen1_;
   Float_t  etaJetGen1_;
   Int_t    partIdJetGen1_;

   Float_t   ptPart1_;
   Float_t    ePart1_;
   Float_t  phiPart1_;
   Float_t  etaPart1_;

   Float_t  eChargedHadronsJet1_[20];
   Float_t  ePhotonsJet1_[20];
   Float_t  eNeutralHadronsJet1_[20];
   Float_t  eMuonsJet1_[20];
   Float_t  eElectronsJet1_[20];
   Float_t  eHFHadronsJet1_[20];
   Float_t  eHFEMJet1_[20];

   Int_t  nChargedHadronsJet1_[20];
   Int_t  nPhotonsJet1_[20];
   Int_t  nNeutralHadronsJet1_[20];
   Int_t  nMuonsJet1_[20];
   Int_t  nElectronsJet1_[20];
   Int_t  nHFHadronsJet1_[20];
   Int_t  nHFEMJet1_[20];

   Int_t  nPFCand1_;
   Float_t  ePFCand1_[100];
   Float_t  ptPFCand1_[100];
   Float_t  etaPFCand1_[100];
   Float_t  phiPFCand1_[100];
   Int_t  particleTypePFCand1_[100];

   Int_t  iJet2_[20];
   Float_t  ptJet2_[20];
   Float_t   eJet2_[20];
   Float_t phiJet2_[20];
   Float_t etaJet2_[20];

   Float_t   ptJetGen2_;
   Float_t    eJetGen2_;
   Float_t  phiJetGen2_;
   Float_t  etaJetGen2_;
   Int_t    partIdJetGen2_;

   Float_t   ptPart2_;
   Float_t    ePart2_;
   Float_t  phiPart2_;
   Float_t  etaPart2_;

   Float_t  eChargedHadronsJet2_[20];
   Float_t  ePhotonsJet2_[20];
   Float_t  eNeutralHadronsJet2_[20];
   Float_t  eMuonsJet2_[20];
   Float_t  eElectronsJet2_[20];
   Float_t  eHFHadronsJet2_[20];
   Float_t  eHFEMJet2_[20];

   Int_t  nChargedHadronsJet2_[20];
   Int_t  nPhotonsJet2_[20];
   Int_t  nNeutralHadronsJet2_[20];
   Int_t  nMuonsJet2_[20];
   Int_t  nElectronsJet2_[20];
   Int_t  nHFHadronsJet2_[20];
   Int_t  nHFEMJet2_[20];

   Int_t  nPFCand2_;
   Float_t  ePFCand2_[100];
   Float_t  ptPFCand2_[100];
   Float_t  etaPFCand2_[100];
   Float_t  phiPFCand2_[100];
   Int_t  particleTypePFCand2_[100];


   Int_t   pdgIdPartJet1_;
   Int_t   pdgIdPartJet2_;

   Float_t epfMet_;
   Float_t phipfMet_;

   TH1F* h1_nEvents_vs_ptEle; 
   TH1F* h1_nEvents_vs_ptMuon; 
   TH1F* h1_passed_vs_ptEle; 
   TH1F* h1_passed_vs_ptMuon; 
   TH1F* h1_deltaRmatching_muons; 
   TH1F* h1_deltaRmatching_electrons; 
   TH1F* h1_deltaRmatching_jet_parton; 
   TH1F* h1_deltaRmatching_genjet_parton; 
   TH1F* h1_deltaRmatching_jet_genjet; 
   TH1F* h1_deltaRmatching_jet_leptonParton;
   TH1F* h1_indexMatchedJet;
   TH1F* h1_indexMatched05Jet;
   TH1F* h1_nMatched_per_event;
   TH1F* h1_nMatched05_per_event;
   TH1F* h1_pdgIdParton1;
   TH1F* h1_pdgIdParton2;
// TH1F* h1_ptHadronicZ; 
// TH1F* h1_deltaRqq; 

   bool DEBUG_VERBOSE_;

};




#endif
