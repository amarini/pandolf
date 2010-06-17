//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 12 18:22:37 2009 by ROOT version 5.22/00
// from TChain pippo/
//////////////////////////////////////////////////////////

#ifndef TreeAnalyzer_MinBias_h
#define TreeAnalyzer_MinBias_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <cmath>


bool DEBUG_ = false;
bool DEBUG_VERBOSE_ = false;



class TreeAnalyzer_MinBias {
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
   Float_t         eMuonsGen[100];
   Float_t         eElectronsGen[100];
   Float_t         ePhotonsGen[100];
   Float_t         eTracksGen[100];
   Float_t         eNeutralHadronsGen[100];
   Float_t         ptMuonsGen[100];
   Float_t         ptElectronsGen[100];
   Float_t         ptPhotonsGen[100];
   Float_t         ptTracksGen[100];
   Float_t         ptNeutralHadronsGen[100];
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
   Int_t           nEventTracks;
   Int_t           nEventChargedHadrons;


   std::string recoType_;
   std::string jetAlgo_;
   std::string algoType_;

   TFile* outfile_;
   TTree* jetTree_;

   Float_t trackRatio_;
   Float_t ptHatMax_;
   Float_t eventWeight_;

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

   Float_t e2ndJetReco_;
   Float_t pt2ndJetReco_;
   Float_t ptCorr2ndJetReco_;
   Float_t eta2ndJetReco_;
   Float_t phi2ndJetReco_;

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

   Int_t  nTracks2ndReco_;
   Int_t  nPhotons2ndReco_;
   Int_t  nNeutralHadrons2ndReco_;
   Int_t  nMuons2ndReco_;
   Int_t  nElectrons2ndReco_;
   Int_t  nHFHadrons2ndReco_;
   Int_t  nHFEM2ndReco_;

   Float_t sumEt_pfMet_;
   Float_t pfMet_;


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
   TBranch        *b_eMuonsGen;   //!
   TBranch        *b_eElectronsGen;   //!
   TBranch        *b_ePhotonsGen;   //!
   TBranch        *b_eTracksGen;   //!
   TBranch        *b_eNeutralHadronsGen;   //!
   TBranch        *b_ptMuonsGen;   //!
   TBranch        *b_ptElectronsGen;   //!
   TBranch        *b_ptPhotonsGen;   //!
   TBranch        *b_ptTracksGen;   //!
   TBranch        *b_ptNeutralHadronsGen;   //!
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
   TBranch        *b_nEventTracks;   //!
   TBranch        *b_nEventChargedHadrons;   //!

   TreeAnalyzer_MinBias(std::string dataset, const std::string& recoType, const std::string& jetAlgo, TTree *tree=0);
   virtual ~TreeAnalyzer_MinBias();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void BookStuff();

};

#endif

#ifdef TreeAnalyzer_MinBias_cxx
TreeAnalyzer_MinBias::TreeAnalyzer_MinBias(std::string dataset, const std::string& recoType, const std::string& jetAlgo, TTree* tree)
{

   recoType_ = recoType;
   jetAlgo_ = jetAlgo;
   algoType_ = (recoType=="pf") ? recoType+jetAlgo : jetAlgo;

   Int_t nFiles;
   Float_t xsection; //will cut on pthat, so will have to subtract x-section of next bin to have effective x-section

   if( dataset=="PhotonJet_Summer09_Pt15" ) {
     xsection = 288813. - 32203.8;
     nFiles = 222;
     ptHatMax_ = 30.;
   } else if( dataset=="PhotonJet_Summer09_Pt30" ) {
     xsection = 32203.8 - 1012.08;
     nFiles = 210;
     ptHatMax_ = 80.;
   } else if( dataset=="PhotonJet_Summer09_Pt80" ) {
     xsection = 1012.08 - 51.36;
     nFiles = 207;
     ptHatMax_ = 170.;
   } else if( dataset=="PhotonJet_Summer09_Pt170" ) {
     xsection = 51.36 - 4.193;
     nFiles = 310;
     ptHatMax_ = 300.;
   } else if( dataset=="PhotonJet_Summer09_Pt300" ) {
     xsection =  4.193 - 0.45125;
     nFiles = 220;
     ptHatMax_ = 470.;
   } else if( dataset=="PhotonJet_Summer09_Pt470" ) {
     xsection =  0.45125 - 0.02;
     nFiles = 206;
     ptHatMax_ = 800.;
   } else if( dataset=="PhotonJet_Summer09_Pt800" ) {
     xsection = 0.02 - 0.000268;
     nFiles = 256;
     ptHatMax_ = 1400.;
   } else if( dataset=="PhotonJet_Summer09_Pt1400" ) {
     xsection = 0.000268;
     nFiles = 214;
     ptHatMax_ = 14000.;
   } else if( dataset=="QCD_Summer09_Pt15" ) {
     nFiles = 1054;
     xsection = 1458126879.8;
     ptHatMax_ = 30.;
   } else if( dataset=="QCD_Summer09_Pt30" ) {
     nFiles = 1698;
     xsection = 109005537.31;
     ptHatMax_ = 80.;
   } else if( dataset=="QCD_Summer09_Pt80" ) {
     nFiles = 438;
     xsection = 1936120.4893;
     ptHatMax_ = 170.;
   } else if( dataset=="QCD_Summer09_Pt170" ) {
     nFiles = 416;
     xsection = 62508.776856;
     ptHatMax_ = 300.;
   } else if( dataset=="QCD_Summer09_Pt300" ) {
     nFiles = 420;
     xsection = 3669.4197667;
     ptHatMax_ = 470.;
   } else if( dataset=="QCD_Summer09_Pt470" ) {
     nFiles = 434;
     xsection = 315.32221016;
     ptHatMax_ = 800.;
   } else if( dataset=="QCD_Summer09_Pt800" ) {
     nFiles = 425;
     xsection = 11.94070485;
     ptHatMax_ = 1400.;
   } else if( dataset=="QCD_Summer09_Pt1400" ) {
     nFiles = 425;
     xsection = 0.17207350709;
     ptHatMax_ = 2200.;
   } else if( dataset=="MC_MinBias_900GeV" ) {
     nFiles = 335;
     xsection = 1.;
     ptHatMax_ = 10000.;
   } else if( dataset=="MC_900GeV_V8K" ) {
     nFiles = 368;
     xsection = 1.;
     ptHatMax_ = 10000.;
   } else if( dataset=="MC_900GeV_V8I"||dataset=="MC_900GeV_V8I_NOTRIG" ) {
     nFiles = 335;
     xsection = 1.;
     ptHatMax_ = 10000.;
   } else if( dataset=="MC_2360GeV_NOTRIG" ) {
     nFiles = 368;
     xsection = 1.;
     ptHatMax_ = 10000.;
   } else if( dataset=="DATA_900GeV" ) {
     nFiles = 253;
     xsection = 1.;
     ptHatMax_ = 10000.;
   } else if( dataset=="DATA_2360GeV" ) {
     nFiles = 3;
     xsection = 1.;
     ptHatMax_ = 10000.;
   } else {
     std::cout << "WARNING! Dataset '" << dataset << "' has no hardwired cross-section, setting it to -1." << std::endl;
     xsection = -1.;
     nFiles = 50;
     ptHatMax_ = 10000.;
   }


   if (tree == 0) {

#ifdef SINGLE_TREE

      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("pippo");

#else // SINGLE_TREE

      TChain * chain = new TChain("pippo","");

      for(int i=1; i<(nFiles+1); ++i) {
        char treePath[100];
        sprintf(treePath, "%s/output_%d.root/myanalysis/pippo", dataset.c_str(), i);
        chain->Add(treePath);
      }

  
      tree = chain;
      std::cout << "-> Added " << nFiles << " files. Tree has " << tree->GetEntries() << " entries." << std::endl;

#endif // SINGLE_TREE

   }


   //will cut on pt_hat, so have to divide only by correct number of events:
   char cutOnPtHat[60];
   sprintf( cutOnPtHat, "genpt<%f", ptHatMax_);
   eventWeight_ = xsection/(Float_t)tree->GetEntries(cutOnPtHat);

   Init(tree);

  
   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova2ndLevel_"+dataset;
   else {
    if(dataset!="") outfileName = "MinBias_2ndLevelTree_"+dataset;
    else outfileName = "MinBias_2ndLevelTree";
   }

   outfileName = outfileName + "_" + algoType_ + ".root";

   outfile_ = new TFile(outfileName.c_str(), "RECREATE");
   

}

TreeAnalyzer_MinBias::~TreeAnalyzer_MinBias()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   outfile_->cd();
   jetTree_->Write();
   outfile_->Write();
   outfile_->Close();

}

Int_t TreeAnalyzer_MinBias::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeAnalyzer_MinBias::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeAnalyzer_MinBias::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   std::string branchName;

   fChain->SetBranchAddress("genpt", &genpt, &b_genpt);
   fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
   fChain->SetBranchAddress("store", &store, &b_store);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", statusMC, &b_statusMC);
   fChain->SetBranchAddress("ptMC ", ptMC , &b_ptMC );
   fChain->SetBranchAddress("eMC  ", eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", phiMC, &b_phiMC);
   fChain->SetBranchAddress("nSIM", &nSIM, &b_nSIM);
   fChain->SetBranchAddress("pdgIdSIM", pdgIdSIM, &b_pdgIdSIM);
   fChain->SetBranchAddress("statusSIM", statusSIM, &b_statusSIM);
   fChain->SetBranchAddress("ptSIM ", ptSIM , &b_ptSIM );
   fChain->SetBranchAddress("eSIM  ", eSIM  , &b_eSIM  );
   fChain->SetBranchAddress("etaSIM", etaSIM, &b_etaSIM);
   fChain->SetBranchAddress("phiSIM", phiSIM, &b_phiSIM);
   fChain->SetBranchAddress("rSIM", rSIM, &b_rSIM);
   fChain->SetBranchAddress("zSIM", zSIM, &b_zSIM);
   fChain->SetBranchAddress("nPF", &nPF, &b_nPF);
   fChain->SetBranchAddress("pdgIdPF", pdgIdPF, &b_pdgIdPF);
   fChain->SetBranchAddress("ptPF ", ptPF , &b_ptPF );
   fChain->SetBranchAddress("ePF  ", ePF  , &b_ePF  );
   fChain->SetBranchAddress("etaPF", etaPF, &b_etaPF);
   fChain->SetBranchAddress("phiPF", phiPF, &b_phiPF);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("ptPhot ", ptPhot , &b_ptPhot );
   fChain->SetBranchAddress("ePhot  ", ePhot  , &b_ePhot  );
   fChain->SetBranchAddress("escPhot  ", escPhot  , &b_escPhot  );
   fChain->SetBranchAddress("etaPhot", etaPhot, &b_etaPhot);
   fChain->SetBranchAddress("phiPhot", phiPhot, &b_phiPhot);
   fChain->SetBranchAddress("nconvPhot", &nconvPhot, &b_nconvPhot);
   fChain->SetBranchAddress("chi2convPhot", chi2convPhot, &b_chi2convPhot);
   fChain->SetBranchAddress("ndofconvPhot", ndofconvPhot, &b_ndofconvPhot);
   fChain->SetBranchAddress("rconvPhot", rconvPhot, &b_rconvPhot);
   fChain->SetBranchAddress("phiconvPhot", phiconvPhot, &b_phiconvPhot);
   fChain->SetBranchAddress("zconvPhot", zconvPhot, &b_zconvPhot);
   fChain->SetBranchAddress("ntrkconvPhot", ntrkconvPhot, &b_ntrkconvPhot);
   fChain->SetBranchAddress("eovpconvPhot", eovpconvPhot, &b_eovpconvPhot);
   fChain->SetBranchAddress("etaecalconvPhot", etaecalconvPhot, &b_etaecalconvPhot);
   fChain->SetBranchAddress("phiecalconvPhot", phiecalconvPhot, &b_phiecalconvPhot);
   fChain->SetBranchAddress("eecalconvPhot", eecalconvPhot, &b_eecalconvPhot);
   fChain->SetBranchAddress("algoconvPhot", algoconvPhot, &b_algoconvPhot);
   fChain->SetBranchAddress("d0convPhot", d0convPhot, &b_d0convPhot);
   fChain->SetBranchAddress("detaecalconvPhot", detaecalconvPhot, &b_detaecalconvPhot);
   fChain->SetBranchAddress("dphiecalconvPhot", dphiecalconvPhot, &b_dphiecalconvPhot);
   fChain->SetBranchAddress("dphivtxconvPhot", dphivtxconvPhot, &b_dphivtxconvPhot);
   fChain->SetBranchAddress("pairsepconvPhot", pairsepconvPhot, &b_pairsepconvPhot);
   fChain->SetBranchAddress("pairmassconvPhot", pairmassconvPhot, &b_pairmassconvPhot);
   fChain->SetBranchAddress("trchi21convPhot", trchi21convPhot, &b_trchi21convPhot);
   fChain->SetBranchAddress("trndof1convPhot", trndof1convPhot, &b_trndof1convPhot);
   fChain->SetBranchAddress("trqual1convPhot", trqual1convPhot, &b_trqual1convPhot);
   fChain->SetBranchAddress("trpt1convPhot", trpt1convPhot, &b_trpt1convPhot);
   fChain->SetBranchAddress("trerr1convPhot", trerr1convPhot, &b_trerr1convPhot);
   fChain->SetBranchAddress("trchi22convPhot", trchi22convPhot, &b_trchi22convPhot);
   fChain->SetBranchAddress("trndof2convPhot", trndof2convPhot, &b_trndof2convPhot);
   fChain->SetBranchAddress("trqual2convPhot", trqual2convPhot, &b_trqual2convPhot);
   fChain->SetBranchAddress("trpt2convPhot", trpt2convPhot, &b_trpt2convPhot);
   fChain->SetBranchAddress("trerr2convPhot", trerr2convPhot, &b_trerr2convPhot);
   fChain->SetBranchAddress("pid_isEM", pid_isEM, &b_pid_isEM);
   fChain->SetBranchAddress("pid_isLoose", pid_isLoose, &b_pid_isLoose);
   fChain->SetBranchAddress("pid_isTight", pid_isTight, &b_pid_isTight);
   fChain->SetBranchAddress("pid_jurECAL", pid_jurECAL, &b_pid_jurECAL);
   fChain->SetBranchAddress("pid_twrHCAL", pid_twrHCAL, &b_pid_twrHCAL);
   fChain->SetBranchAddress("pid_HoverE", pid_HoverE, &b_pid_HoverE);
   fChain->SetBranchAddress("pid_hlwTrack", pid_hlwTrack, &b_pid_hlwTrack);
   fChain->SetBranchAddress("pid_etawid", pid_etawid, &b_pid_etawid);
   fChain->SetBranchAddress("ptiso0015Phot", ptiso0015Phot, &b_ptiso0015Phot);
   fChain->SetBranchAddress("ntrkiso0015Phot", ntrkiso0015Phot, &b_ntrkiso0015Phot);
   fChain->SetBranchAddress("ptiso035Phot", ptiso035Phot, &b_ptiso035Phot);
   fChain->SetBranchAddress("ntrkiso035Phot", ntrkiso035Phot, &b_ntrkiso035Phot);
   fChain->SetBranchAddress("ptiso04Phot", ptiso04Phot, &b_ptiso04Phot);
   fChain->SetBranchAddress("ntrkiso04Phot", ntrkiso04Phot, &b_ntrkiso04Phot);
   fChain->SetBranchAddress("ptiso05Phot", ptiso05Phot, &b_ptiso05Phot);
   fChain->SetBranchAddress("ntrkiso05Phot", ntrkiso05Phot, &b_ntrkiso05Phot);
   fChain->SetBranchAddress("ptiso07Phot", ptiso07Phot, &b_ptiso07Phot);
   fChain->SetBranchAddress("ntrkiso07Phot", ntrkiso07Phot, &b_ntrkiso07Phot);
   fChain->SetBranchAddress("ptiso1Phot", ptiso1Phot, &b_ptiso1Phot);
   fChain->SetBranchAddress("ntrkiso1Phot", ntrkiso1Phot, &b_ntrkiso1Phot);
   fChain->SetBranchAddress("hcalovecal01Phot", hcalovecal01Phot, &b_hcalovecal01Phot);
   fChain->SetBranchAddress("hcalovecal015Phot", hcalovecal015Phot, &b_hcalovecal015Phot);
   fChain->SetBranchAddress("hcalovecal04Phot", hcalovecal04Phot, &b_hcalovecal04Phot);
   fChain->SetBranchAddress("hcalovecal05Phot", hcalovecal05Phot, &b_hcalovecal05Phot);
   fChain->SetBranchAddress("hcalovecal07Phot", hcalovecal07Phot, &b_hcalovecal07Phot);
   fChain->SetBranchAddress("hcalovecal1Phot", hcalovecal1Phot, &b_hcalovecal1Phot);
   fChain->SetBranchAddress("ecaliso01Phot", ecaliso01Phot, &b_ecaliso01Phot);
   fChain->SetBranchAddress("ecaliso015Phot", ecaliso015Phot, &b_ecaliso015Phot);
   fChain->SetBranchAddress("ecaliso04Phot", ecaliso04Phot, &b_ecaliso04Phot);
   fChain->SetBranchAddress("ecaliso05Phot", ecaliso05Phot, &b_ecaliso05Phot);
   fChain->SetBranchAddress("ecaliso07Phot", ecaliso07Phot, &b_ecaliso07Phot);
   fChain->SetBranchAddress("ecaliso1Phot", ecaliso1Phot, &b_ecaliso1Phot);
   fChain->SetBranchAddress("LATPhot", LATPhot, &b_LATPhot);
   fChain->SetBranchAddress("sMajMajPhot", sMajMajPhot, &b_sMajMajPhot);
   fChain->SetBranchAddress("sMinMinPhot", sMinMinPhot, &b_sMinMinPhot);
   fChain->SetBranchAddress("alphaPhot", alphaPhot, &b_alphaPhot);
   fChain->SetBranchAddress("sEtaEtaPhot", sEtaEtaPhot, &b_sEtaEtaPhot);
   fChain->SetBranchAddress("sPhiPhiPhot", sPhiPhiPhot, &b_sPhiPhiPhot);
   fChain->SetBranchAddress("E1Phot", E1Phot, &b_E1Phot);
   fChain->SetBranchAddress("E9Phot", E9Phot, &b_E9Phot);
   fChain->SetBranchAddress("E25Phot", E25Phot, &b_E25Phot);
   fChain->SetBranchAddress("FisherPhot", FisherPhot, &b_FisherPhot);
   branchName = "nJetGen_"+jetAlgo_;
   fChain->SetBranchAddress(branchName.c_str(), &nJetGen, &b_nJetGen);
   branchName = "ptJetGen_"+jetAlgo_+" ";
   fChain->SetBranchAddress(branchName.c_str(), ptJetGen, &b_ptJetGen);
   branchName = "eJetGen_"+jetAlgo_+"  ";
   fChain->SetBranchAddress(branchName.c_str(), eJetGen, &b_eJetGen);
   branchName = "etaJetGen_"+jetAlgo_;
   fChain->SetBranchAddress(branchName.c_str(), etaJetGen, &b_etaJetGen);
   if( jetAlgo_=="akt5") {
     branchName = "nMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nMuonsGen, &b_nMuonsGen);
     branchName = "nElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nElectronsGen, &b_nElectronsGen);
     branchName = "nPhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nPhotonsGen, &b_nPhotonsGen);
     branchName = "nTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nTracksGen, &b_nTracksGen);
     branchName = "nNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nNeutralHadronsGen, &b_nNeutralHadronsGen);
     branchName = "eMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eMuonsGen, &b_eMuonsGen);
     branchName = "eElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eElectronsGen, &b_eElectronsGen);
     branchName = "ePhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ePhotonsGen, &b_ePhotonsGen);
     branchName = "eTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eTracksGen, &b_eTracksGen);
     branchName = "eNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eNeutralHadronsGen, &b_eNeutralHadronsGen);
     branchName = "ptMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ptMuonsGen, &b_ptMuonsGen);
     branchName = "ptElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ptElectronsGen, &b_ptElectronsGen);
     branchName = "ptPhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ptPhotonsGen, &b_ptPhotonsGen);
     branchName = "ptTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ptTracksGen, &b_ptTracksGen);
     branchName = "ptNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ptNeutralHadronsGen, &b_ptNeutralHadronsGen);
     branchName = "etaMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), etaMuonsGen, &b_etaMuonsGen);
     branchName = "etaElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), etaElectronsGen, &b_etaElectronsGen);
     branchName = "etaPhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), etaPhotonsGen, &b_etaPhotonsGen);
     branchName = "etaTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), etaTracksGen, &b_etaTracksGen);
     branchName = "etaNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), etaNeutralHadronsGen, &b_etaNeutralHadronsGen);
     branchName = "phiMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), phiMuonsGen, &b_phiMuonsGen);
     branchName = "phiElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), phiElectronsGen, &b_phiElectronsGen);
     branchName = "phiPhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), phiPhotonsGen, &b_phiPhotonsGen);
     branchName = "phiTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), phiTracksGen, &b_phiTracksGen);
     branchName = "phiNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), phiNeutralHadronsGen, &b_phiNeutralHadronsGen);
   }
   branchName = "nJet_"+algoType_;
   fChain->SetBranchAddress(branchName.c_str(), &nJet, &b_nJet);
   branchName = "ptJet_"+algoType_+" ";
   fChain->SetBranchAddress(branchName.c_str(), ptJet, &b_ptJet);
   branchName = "ptCorrJet_"+algoType_+" ";
   fChain->SetBranchAddress(branchName.c_str(), ptCorrJet, &b_ptCorrJet);
   branchName = "eJet_"+algoType_+"  ";
   fChain->SetBranchAddress(branchName.c_str(), eJet, &b_eJet);
   branchName = "etaJet_"+algoType_;
   fChain->SetBranchAddress(branchName.c_str(), etaJet, &b_etaJet);
   branchName = "phiJet_"+algoType_;
   fChain->SetBranchAddress(branchName.c_str(), phiJet, &b_phiJet);
   if( recoType_=="pf" ) {
     branchName = "nChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nChargedHadrons, &b_nChargedHadrons);
     branchName = "nPhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nPhotons, &b_nPhotons);
     branchName = "nMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nMuons, &b_nMuons);
     branchName = "nElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nElectrons, &b_nElectrons);
     branchName = "nNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nNeutralHadrons, &b_nNeutralHadrons);
     branchName = "nHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nHFHadrons, &b_nHFHadrons);
     branchName = "nHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nHFEM, &b_nHFEM);
     branchName = "eChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eChargedHadrons, &b_eChargedHadrons);
     branchName = "ePhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ePhotons, &b_ePhotons);
     branchName = "eMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eMuons, &b_eMuons);
     branchName = "eElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eElectrons, &b_eElectrons);
     branchName = "eNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eNeutralHadrons, &b_eNeutralHadrons);
     branchName = "eHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eHFHadrons, &b_eHFHadrons);
     branchName = "eHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eHFEM, &b_eHFEM);
     branchName = "ptChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptChargedHadrons, &b_ptChargedHadrons);
     branchName = "ptPhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptPhotons, &b_ptPhotons);
     branchName = "ptMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptMuons, &b_ptMuons);
     branchName = "ptElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptElectrons, &b_ptElectrons);
     branchName = "ptNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptNeutralHadrons, &b_ptNeutralHadrons);
     branchName = "ptHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptHFHadrons, &b_ptHFHadrons);
     branchName = "ptHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ptHFEM, &b_ptHFEM);
     branchName = "etaChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaChargedHadrons, &b_etaChargedHadrons);
     branchName = "etaPhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaPhotons, &b_etaPhotons);
     branchName = "etaMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaMuons, &b_etaMuons);
     branchName = "etaElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaElectrons, &b_etaElectrons);
     branchName = "etaNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaNeutralHadrons, &b_etaNeutralHadrons);
     branchName = "etaHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaHFHadrons, &b_etaHFHadrons);
     branchName = "etaHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), etaHFEM, &b_etaHFEM);
     branchName = "phiChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiChargedHadrons, &b_phiChargedHadrons);
     branchName = "phiPhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiPhotons, &b_phiPhotons);
     branchName = "phiMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiMuons, &b_phiMuons);
     branchName = "phiElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiElectrons, &b_phiElectrons);
     branchName = "phiNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiNeutralHadrons, &b_phiNeutralHadrons);
     branchName = "phiHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiHFHadrons, &b_phiHFHadrons);
     branchName = "phiHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), phiHFEM, &b_phiHFEM);
   }
   fChain->SetBranchAddress("sMet  ", &sMet  , &b_sMet);
   fChain->SetBranchAddress("eMet  ", &eMet  , &b_eMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("stcMet  ", &stcMet  , &b_stcMet);
   fChain->SetBranchAddress("etcMet  ", &etcMet  , &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("spfMet  ", &spfMet  , &b_spfMet);
   fChain->SetBranchAddress("epfMet  ", &epfMet  , &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("sMetGen  ", &sMetGen  , &b_sMetGen);
   fChain->SetBranchAddress("eMetGen  ", &eMetGen  , &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("sMetGen2  ", &sMetGen2  , &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2  ", &eMetGen2  , &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("hltPass", &hltPass, &b_hltPass);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("hltNamesLen", &hltNamesLen, &b_hltNamesLen);
   fChain->SetBranchAddress("HLTNames", HLTNames, &b_HLTNames);
   fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("vntracks", &vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", &vchi2, &b_vchi2);
   fChain->SetBranchAddress("nEventTracks", &nEventTracks, &b_nEventTracks);
   fChain->SetBranchAddress("nChargedHadrons", &nEventChargedHadrons, &b_nEventChargedHadrons);
   Notify();
}

Bool_t TreeAnalyzer_MinBias::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeAnalyzer_MinBias::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeAnalyzer_MinBias::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeAnalyzer_MinBias_cxx
