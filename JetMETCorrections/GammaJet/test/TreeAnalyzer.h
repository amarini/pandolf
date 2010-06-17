//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 12 18:22:37 2009 by ROOT version 5.22/00
// from TChain pippo/
//////////////////////////////////////////////////////////


//------------------------------------------------------------------
//
//    Base TreeAnalyzer class. Abstract.
//    Reads output of GammaJetAnalyzer, and produces subtrees. 
//    Loop() and BookStuff() methods have to be implemented.
//    Additional data members should be assigned.
//
//------------------------------------------------------------------



#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>




class TreeAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         genpt;
   Bool_t          isMC;
   Int_t           store;
   Int_t           lbn;
   Int_t           bx;
   Int_t           orbit;
   Int_t           run;
   Int_t           event;
   Int_t           nMC;
   Int_t           pdgIdMC[150];   //[nMC]
   Int_t           motherIDMC[150];   //[nMC]
   Int_t           statusMC[150];   //[nMC]
   Float_t         ptMC [150];   //[nMC]
   Float_t         eMC  [150];   //[nMC]
   Float_t         etaMC[150];   //[nMC]
   Float_t         phiMC[150];   //[nMC]
   Int_t           nSIM;
   Int_t           pdgIdSIM[150];   //[nSIM]
   Int_t           statusSIM[150];   //[nSIM]
   Float_t         ptSIM [150];   //[nSIM]
   Float_t         eSIM  [150];   //[nSIM]
   Float_t         etaSIM[150];   //[nSIM]
   Float_t         phiSIM[150];   //[nSIM]
   Float_t         rSIM[150];   //[nSIM]
   Float_t         zSIM[150];   //[nSIM]
   Int_t           nPF;
   Int_t           pdgIdPF[150];   //[nPF]
   Float_t         ptPF [150];   //[nPF]
   Float_t         ePF  [150];   //[nPF]
   Float_t         etaPF[150];   //[nPF]
   Float_t         phiPF[150];   //[nPF]
   Int_t           nPhot;
   Float_t         ptPhot [40];   //[nPhot]
   Float_t         ePhot  [40];   //[nPhot]
   Float_t         escPhot  [40];   //[nPhot]
   Float_t         etaPhot[40];   //[nPhot]
   Float_t         phiPhot[40];   //[nPhot]
   Int_t           nconvPhot;
   Float_t         chi2convPhot[10];   //[nconvPhot]
   Float_t         ndofconvPhot[10];   //[nconvPhot]
   Float_t         rconvPhot[10];   //[nconvPhot]
   Float_t         phiconvPhot[10];   //[nconvPhot]
   Float_t         zconvPhot[10];   //[nconvPhot]
   Int_t           ntrkconvPhot[10];   //[nconvPhot]
   Float_t         eovpconvPhot[10];   //[nconvPhot]
   Float_t         etaecalconvPhot[10];   //[nconvPhot]
   Float_t         phiecalconvPhot[10];   //[nconvPhot]
   Float_t         eecalconvPhot[10];   //[nconvPhot]
   Int_t           algoconvPhot[10];   //[nconvPhot]
   Float_t         d0convPhot[10];   //[nconvPhot]
   Float_t         detaecalconvPhot[10];   //[nconvPhot]
   Float_t         dphiecalconvPhot[10];   //[nconvPhot]
   Float_t         dphivtxconvPhot[10];   //[nconvPhot]
   Float_t         pairsepconvPhot[10];   //[nconvPhot]
   Float_t         pairmassconvPhot[10];   //[nconvPhot]
   Float_t         trchi21convPhot[10];   //[nconvPhot]
   Float_t         trndof1convPhot[10];   //[nconvPhot]
   Int_t           trqual1convPhot[10];   //[nconvPhot]
   Float_t         trpt1convPhot[10];   //[nconvPhot]
   Float_t         trerr1convPhot[10];   //[nconvPhot]
   Float_t         trchi22convPhot[10];   //[nconvPhot]
   Float_t         trndof2convPhot[10];   //[nconvPhot]
   Int_t           trqual2convPhot[10];   //[nconvPhot]
   Float_t         trpt2convPhot[10];   //[nconvPhot]
   Float_t         trerr2convPhot[10];   //[nconvPhot]
   Bool_t          pid_isEM[40];   //[nPhot]
   Bool_t          pid_isLoose[40];   //[nPhot]
   Bool_t          pid_isTight[40];   //[nPhot]
   Float_t         pid_jurECAL[40];   //[nPhot]
   Float_t         pid_twrHCAL[40];   //[nPhot]
   Float_t         pid_HoverE[40];   //[nPhot]
   Float_t         pid_hlwTrack[40];   //[nPhot]
   Float_t         pid_etawid[40];   //[nPhot]
   Float_t         ptiso0015Phot[40];   //[nPhot]
   Int_t           ntrkiso0015Phot[40];   //[nPhot]
   Float_t         ptiso035Phot[40];   //[nPhot]
   Int_t           ntrkiso035Phot[40];   //[nPhot]
   Float_t         ptiso04Phot[40];   //[nPhot]
   Int_t           ntrkiso04Phot[40];   //[nPhot]
   Float_t         ptiso05Phot[40];   //[nPhot]
   Int_t           ntrkiso05Phot[40];   //[nPhot]
   Float_t         ptiso07Phot[40];   //[nPhot]
   Int_t           ntrkiso07Phot[40];   //[nPhot]
   Float_t         ptiso1Phot[40];   //[nPhot]
   Int_t           ntrkiso1Phot[40];   //[nPhot]
   Float_t         hcalovecal01Phot[40];   //[nPhot]
   Float_t         hcalovecal015Phot[40];   //[nPhot]
   Float_t         hcalovecal04Phot[40];   //[nPhot]
   Float_t         hcalovecal05Phot[40];   //[nPhot]
   Float_t         hcalovecal07Phot[40];   //[nPhot]
   Float_t         hcalovecal1Phot[40];   //[nPhot]
   Float_t         ecaliso01Phot[40];   //[nPhot]
   Float_t         ecaliso015Phot[40];   //[nPhot]
   Float_t         ecaliso04Phot[40];   //[nPhot]
   Float_t         ecaliso05Phot[40];   //[nPhot]
   Float_t         ecaliso07Phot[40];   //[nPhot]
   Float_t         ecaliso1Phot[40];   //[nPhot]
   Float_t         LATPhot[40];   //[nPhot]
   Float_t         sMajMajPhot[40];   //[nPhot]
   Float_t         sMinMinPhot[40];   //[nPhot]
   Float_t         alphaPhot[40];   //[nPhot]
   Float_t         sEtaEtaPhot[40];   //[nPhot]
   Float_t         sPhiPhiPhot[40];   //[nPhot]
   Float_t         E1Phot[40];   //[nPhot]
   Float_t         E9Phot[40];   //[nPhot]
   Float_t         E25Phot[40];   //[nPhot]
   Float_t         FisherPhot[40];   //[nPhot]
   Int_t           nJetGen;
   Float_t         ptJetGen[100];  
   Float_t         eJetGen[100];   
   Float_t         etaJetGen[100]; 
   Float_t         phiJetGen[100];
   Int_t           nMuonsGen[100];
   Int_t           nElectronsGen[100];
   Int_t           nPhotonsGen[100];
   Int_t           nTracksGen[100];
   Int_t           nNeutralHadronsGen[100];
   Int_t           nHFHadronsGen[100];
   Int_t           nHFEMGen[100];
   Float_t         eMuonsGen[100];
   Float_t         eElectronsGen[100];
   Float_t         ePhotonsGen[100];
   Float_t         eTracksGen[100];
   Float_t         eNeutralHadronsGen[100];
   Float_t         eHFHadronsGen[100];
   Float_t         eHFEMGen[100];
   Float_t         ptMuonsGen[100];
   Float_t         ptElectronsGen[100];
   Float_t         ptPhotonsGen[100];
   Float_t         ptTracksGen[100];
   Float_t         ptNeutralHadronsGen[100];
   Float_t         ptHFHadronsGen[100];
   Float_t         ptHFEMGen[100];
   Float_t         phiMuonsGen[100];
   Float_t         phiElectronsGen[100];
   Float_t         phiPhotonsGen[100];
   Float_t         phiTracksGen[100];
   Float_t         phiNeutralHadronsGen[100];
   Float_t         etaMuonsGen[100];
   Float_t         etaElectronsGen[100];
   Float_t         etaPhotonsGen[100];
   Float_t         etaTracksGen[100];
   Float_t         etaNeutralHadronsGen[100];
   Int_t           nJet;
   Float_t         ptJet[100];  
   Float_t         ptCorrJet[100];  
   Float_t         eJet[100];   
   Float_t         etaJet[100]; 
   Float_t         phiJet[100];
   Float_t         emfJet[100];
   Int_t           nChargedHadrons[100];
   Int_t           nPhotons[100];
   Int_t           nMuons[100];
   Int_t           nElectrons[100];
   Int_t           nNeutralHadrons[100];
   Int_t           nHFHadrons[100];
   Int_t           nHFEM[100];
   Float_t         eChargedHadrons[100];
   Float_t         ePhotons[100];
   Float_t         eMuons[100];
   Float_t         eElectrons[100];
   Float_t         eNeutralHadrons[100];
   Float_t         eHFHadrons[100];
   Float_t         eHFEM[100];
   Float_t         ptChargedHadrons[100];
   Float_t         ptPhotons[100];
   Float_t         ptMuons[100];
   Float_t         ptElectrons[100];
   Float_t         ptNeutralHadrons[100];
   Float_t         ptHFHadrons[100];
   Float_t         ptHFEM[100];
   Float_t         phiChargedHadrons[100];
   Float_t         phiPhotons[100];
   Float_t         phiMuons[100];
   Float_t         phiElectrons[100];
   Float_t         phiNeutralHadrons[100];
   Float_t         phiHFHadrons[100];
   Float_t         phiHFEM[100];
   Float_t         etaChargedHadrons[100];
   Float_t         etaPhotons[100];
   Float_t         etaMuons[100];
   Float_t         etaElectrons[100];
   Float_t         etaNeutralHadrons[100];
   Float_t         etaHFHadrons[100];
   Float_t         etaHFEM[100];
   Float_t         sMet  ;
   Float_t         eMet  ;
   Float_t         phiMet;
   Float_t         stcMet  ;
   Float_t         etcMet  ;
   Float_t         phitcMet;
   Float_t         spfMet  ;
   Float_t         epfMet  ;
   Float_t         phipfMet;
   Float_t         sMetGen  ;
   Float_t         eMetGen  ;
   Float_t         phiMetGen;
   Float_t         sMetGen2  ;
   Float_t         eMetGen2  ;
   Float_t         phiMetGen2;
   Bool_t          hltPass;
   Int_t           nHLT;
   Int_t           hltNamesLen;
   Char_t          HLTNames[6000];   //[hltNamesLen]
   Bool_t          HLTResults[200];   //[nHLT]
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         vntracks;
   Float_t         vchi2;


   std::string dataset_;

   std::string recoType_;
   std::string jetAlgo_;
   std::string algoType_;

   std::string analyzerType_;
   std::string flags_;

   typedef std::vector< std::pair<int,int> > GoodLSVector;
   typedef std::map< int, GoodLSVector  >    LSRange ;
   LSRange goodLS_;
   LSRange::const_iterator goodLSCache_; // ptr to list of good LS for last run
   bool filterGoodRuns_;

   typedef std::pair< int, int > RunLumiPair;
   typedef std::map< RunLumiPair, double > LSLumi;
   LSLumi LSLumimap_;
   double totalIntLumi_;

   TFile* outfile_;
   TTree* jetTree_;

   Int_t run_;
   Int_t oldrun_;
   Int_t currentLS_;
   Int_t event_;
   Int_t lbn_;
   Float_t ptHat_;
   Float_t ptHatMin_;
   Float_t ptHatMax_;
   Float_t eventWeight_;

   bool DEBUG_;


   // List of branches
   TBranch        *b_genpt;   //!
   TBranch        *b_isMC;   //!
   TBranch        *b_store;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_pdgIdMC;   //!
   TBranch        *b_motherIDMC;   //!
   TBranch        *b_statusMC;   //!
   TBranch        *b_ptMC ;   //!
   TBranch        *b_eMC  ;   //!
   TBranch        *b_etaMC;   //!
   TBranch        *b_phiMC;   //!
   TBranch        *b_nSIM;   //!
   TBranch        *b_pdgIdSIM;   //!
   TBranch        *b_statusSIM;   //!
   TBranch        *b_ptSIM ;   //!
   TBranch        *b_eSIM  ;   //!
   TBranch        *b_etaSIM;   //!
   TBranch        *b_phiSIM;   //!
   TBranch        *b_rSIM;   //!
   TBranch        *b_zSIM;   //!
   TBranch        *b_nPF;   //!
   TBranch        *b_pdgIdPF;   //!
   TBranch        *b_ptPF ;   //!
   TBranch        *b_ePF  ;   //!
   TBranch        *b_etaPF;   //!
   TBranch        *b_phiPF;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_ptPhot ;   //!
   TBranch        *b_ePhot  ;   //!
   TBranch        *b_escPhot  ;   //!
   TBranch        *b_etaPhot;   //!
   TBranch        *b_phiPhot;   //!
   TBranch        *b_nconvPhot;   //!
   TBranch        *b_chi2convPhot;   //!
   TBranch        *b_ndofconvPhot;   //!
   TBranch        *b_rconvPhot;   //!
   TBranch        *b_phiconvPhot;   //!
   TBranch        *b_zconvPhot;   //!
   TBranch        *b_ntrkconvPhot;   //!
   TBranch        *b_eovpconvPhot;   //!
   TBranch        *b_etaecalconvPhot;   //!
   TBranch        *b_phiecalconvPhot;   //!
   TBranch        *b_eecalconvPhot;   //!
   TBranch        *b_algoconvPhot;   //!
   TBranch        *b_d0convPhot;   //!
   TBranch        *b_detaecalconvPhot;   //!
   TBranch        *b_dphiecalconvPhot;   //!
   TBranch        *b_dphivtxconvPhot;   //!
   TBranch        *b_pairsepconvPhot;   //!
   TBranch        *b_pairmassconvPhot;   //!
   TBranch        *b_trchi21convPhot;   //!
   TBranch        *b_trndof1convPhot;   //!
   TBranch        *b_trqual1convPhot;   //!
   TBranch        *b_trpt1convPhot;   //!
   TBranch        *b_trerr1convPhot;   //!
   TBranch        *b_trchi22convPhot;   //!
   TBranch        *b_trndof2convPhot;   //!
   TBranch        *b_trqual2convPhot;   //!
   TBranch        *b_trpt2convPhot;   //!
   TBranch        *b_trerr2convPhot;   //!
   TBranch        *b_pid_isEM;   //!
   TBranch        *b_pid_isLoose;   //!
   TBranch        *b_pid_isTight;   //!
   TBranch        *b_pid_jurECAL;   //!
   TBranch        *b_pid_twrHCAL;   //!
   TBranch        *b_pid_HoverE;   //!
   TBranch        *b_pid_hlwTrack;   //!
   TBranch        *b_pid_etawid;   //!
   TBranch        *b_ptiso0015Phot;   //!
   TBranch        *b_ntrkiso0015Phot;   //!
   TBranch        *b_ptiso035Phot;   //!
   TBranch        *b_ntrkiso035Phot;   //!
   TBranch        *b_ptiso04Phot;   //!
   TBranch        *b_ntrkiso04Phot;   //!
   TBranch        *b_ptiso05Phot;   //!
   TBranch        *b_ntrkiso05Phot;   //!
   TBranch        *b_ptiso07Phot;   //!
   TBranch        *b_ntrkiso07Phot;   //!
   TBranch        *b_ptiso1Phot;   //!
   TBranch        *b_ntrkiso1Phot;   //!
   TBranch        *b_hcalovecal01Phot;   //!
   TBranch        *b_hcalovecal015Phot;   //!
   TBranch        *b_hcalovecal04Phot;   //!
   TBranch        *b_hcalovecal05Phot;   //!
   TBranch        *b_hcalovecal07Phot;   //!
   TBranch        *b_hcalovecal1Phot;   //!
   TBranch        *b_ecaliso01Phot;   //!
   TBranch        *b_ecaliso015Phot;   //!
   TBranch        *b_ecaliso04Phot;   //!
   TBranch        *b_ecaliso05Phot;   //!
   TBranch        *b_ecaliso07Phot;   //!
   TBranch        *b_ecaliso1Phot;   //!
   TBranch        *b_LATPhot;   //!
   TBranch        *b_sMajMajPhot;   //!
   TBranch        *b_sMinMinPhot;   //!
   TBranch        *b_alphaPhot;   //!
   TBranch        *b_sEtaEtaPhot;   //!
   TBranch        *b_sPhiPhiPhot;   //!
   TBranch        *b_E1Phot;   //!
   TBranch        *b_E9Phot;   //!
   TBranch        *b_E25Phot;   //!
   TBranch        *b_FisherPhot;   //!
   TBranch        *b_nJetGen;   //!
   TBranch        *b_ptJetGen;   //!
   TBranch        *b_eJetGen;   //!
   TBranch        *b_etaJetGen;   //!
   TBranch        *b_phiJetGen;   //!
   TBranch        *b_nMuonsGen;   //!
   TBranch        *b_nElectronsGen;   //!
   TBranch        *b_nPhotonsGen;   //!
   TBranch        *b_nTracksGen;   //!
   TBranch        *b_nNeutralHadronsGen;   //!
   TBranch        *b_nHFHadronsGen;   //!
   TBranch        *b_nHFEMGen;   //!
   TBranch        *b_eMuonsGen;   //!
   TBranch        *b_eElectronsGen;   //!
   TBranch        *b_ePhotonsGen;   //!
   TBranch        *b_eTracksGen;   //!
   TBranch        *b_eNeutralHadronsGen;   //!
   TBranch        *b_eHFHadronsGen;   //!
   TBranch        *b_eHFEMGen;   //!
   TBranch        *b_ptMuonsGen;   //!
   TBranch        *b_ptElectronsGen;   //!
   TBranch        *b_ptPhotonsGen;   //!
   TBranch        *b_ptTracksGen;   //!
   TBranch        *b_ptNeutralHadronsGen;   //!
   TBranch        *b_ptHFHadronsGen;   //!
   TBranch        *b_ptHFEMGen;   //!
   TBranch        *b_phiMuonsGen;   //!
   TBranch        *b_phiElectronsGen;   //!
   TBranch        *b_phiPhotonsGen;   //!
   TBranch        *b_phiTracksGen;   //!
   TBranch        *b_phiNeutralHadronsGen;   //!
   TBranch        *b_etaMuonsGen;   //!
   TBranch        *b_etaElectronsGen;   //!
   TBranch        *b_etaPhotonsGen;   //!
   TBranch        *b_etaTracksGen;   //!
   TBranch        *b_etaNeutralHadronsGen;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_ptJet;   //!
   TBranch        *b_ptCorrJet;   //!
   TBranch        *b_eJet;   //!
   TBranch        *b_etaJet;   //!
   TBranch        *b_phiJet;   //!
   TBranch        *b_emfJet;   //!
   TBranch        *b_nChargedHadrons;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_nNeutralHadrons;   //!
   TBranch        *b_nHFHadrons;   //!
   TBranch        *b_nHFEM;   //!
   TBranch        *b_eChargedHadrons;   //!
   TBranch        *b_ePhotons;   //!
   TBranch        *b_eMuons;   //!
   TBranch        *b_eElectrons;   //!
   TBranch        *b_eNeutralHadrons;   //!
   TBranch        *b_eHFHadrons;   //!
   TBranch        *b_eHFEM;   //!
   TBranch        *b_ptChargedHadrons;   //!
   TBranch        *b_ptPhotons;   //!
   TBranch        *b_ptMuons;   //!
   TBranch        *b_ptElectrons;   //!
   TBranch        *b_ptNeutralHadrons;   //!
   TBranch        *b_ptHFHadrons;   //!
   TBranch        *b_ptHFEM;   //!
   TBranch        *b_phiChargedHadrons;   //!
   TBranch        *b_phiPhotons;   //!
   TBranch        *b_phiMuons;   //!
   TBranch        *b_phiElectrons;   //!
   TBranch        *b_phiNeutralHadrons;   //!
   TBranch        *b_phiHFHadrons;   //!
   TBranch        *b_phiHFEM;   //!
   TBranch        *b_etaChargedHadrons;   //!
   TBranch        *b_etaPhotons;   //!
   TBranch        *b_etaMuons;   //!
   TBranch        *b_etaElectrons;   //!
   TBranch        *b_etaNeutralHadrons;   //!
   TBranch        *b_etaHFHadrons;   //!
   TBranch        *b_etaHFEM;   //!
   TBranch        *b_sMet;   //!
   TBranch        *b_eMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_stcMet;   //!
   TBranch        *b_etcMet;   //!
   TBranch        *b_phitcMet;   //!
   TBranch        *b_spfMet;   //!
   TBranch        *b_epfMet;   //!
   TBranch        *b_phipfMet;   //!
   TBranch        *b_sMetGen;   //!
   TBranch        *b_eMetGen;   //!
   TBranch        *b_phiMetGen;   //!
   TBranch        *b_sMetGen2;   //!
   TBranch        *b_eMetGen2;   //!
   TBranch        *b_phiMetGen2;   //!
   TBranch        *b_hltPass;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_hltNamesLen;   //!
   TBranch        *b_HLTNames;   //!
   TBranch        *b_HLTResults;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vntracks;   //!
   TBranch        *b_vchi2;   //!

   TreeAnalyzer( const std::string& analyzerType, const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", TTree *tree=0);
   virtual ~TreeAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop()=0;
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     ReadJSONFile(const std::string& json);
   virtual void     ReadCSVFile(const std::string& csv);
   virtual void     UpdateCache();
   virtual bool     isGoodLS ();

   //virtual void BookStuff()=0;

};

#endif
