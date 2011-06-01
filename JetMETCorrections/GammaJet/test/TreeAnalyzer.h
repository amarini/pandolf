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
#include <TH1F.h>
#include <TFile.h>
#include <iostream>


struct GenEventParameters{

  Float_t crossSection;
  Float_t ptHatMin;
  Float_t ptHatMax;

};


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
   Float_t         rhoCalo;
   Float_t         rhoPF;
   Int_t           nMC;
   Int_t           pdgIdMC[150];   //[nMC]
   Int_t           statusMC[150];   //[nMC]
   Int_t           motherIDMC[150];   //[nMC]
   Float_t         ptMC [150];   //[nMC]
   Float_t         eMC  [150];   //[nMC]
   Float_t         etaMC[150];   //[nMC]
   Float_t         phiMC[150];   //[nMC]
   Int_t           nPhot;
   Float_t         ptPhot [40];   //[nPhot]
   Float_t         ePhot  [40];   //[nPhot]
   Float_t         escPhot  [40];   //[nPhot]
   Float_t         eseedPhot  [40];   //[nPhot]
   Float_t         etaPhot[40];   //[nPhot]
   Float_t         phiPhot[40];   //[nPhot]
   Float_t         timePhot[40];   //[nPhot]
   Float_t         e4SwissCrossPhot[40];   //[nPhot]
   Int_t           hasPixelSeedPhot[40];   //[nPhot]
   Bool_t          pid_isEM[40];   //[nPhot]
   Bool_t          pid_isLoose[40];   //[nPhot]
   Bool_t          pid_isTight[40];   //[nPhot]
   Float_t         pid_jurECAL[40];   //[nPhot]
   Float_t         pid_twrHCAL[40];   //[nPhot]
   Float_t         pid_HoverE[40];   //[nPhot]
   Float_t         pid_hlwTrack[40];   //[nPhot]
   Float_t         pid_hlwTrackNoDz[40];   //[nPhot]
   Float_t         pid_etawid[40];   //[nPhot]
   Float_t         ptiso035Phot[40];   //[nPhot]
   Int_t           ntrkiso035Phot[40];   //[nPhot]
   Float_t         ptiso04Phot[40];   //[nPhot]
   Int_t           ntrkiso04Phot[40];   //[nPhot]
   Float_t         hcalovecal04Phot[40];   //[nPhot]
   Float_t         ecaliso04Phot[40];   //[nPhot]
   Float_t         sMajMajPhot[40];   //[nPhot]
   Float_t         sMinMinPhot[40];   //[nPhot]
   Float_t         alphaPhot[40];   //[nPhot]
   Float_t         sEtaEtaPhot[40];   //[nPhot]
   Float_t         sEtaPhiPhot[40];   //[nPhot]
   Float_t         sPhiPhiPhot[40];   //[nPhot]
   Float_t         E1Phot[40];   //[nPhot]
   Float_t         E9Phot[40];   //[nPhot]
   Float_t         E25Phot[40];   //[nPhot]
   Int_t           nJet;
   Float_t         ptJet[100];   
   Float_t         ptCorrJet[100]; 
   Float_t         eJet[100];  
   Float_t         etaJet[100]; 
   Float_t         phiJet[100];
   Float_t         emfJet[100];
   Float_t         n90Jet[100];
   Float_t         n90HitsJet[100]; 
   Float_t         fHPDJet[100]; 
   Float_t         fRBXJet[100];
   Float_t         ptDJet[100];   
   Float_t         rmsCandJet[100];   
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
   Int_t           nJetGen;
   Float_t         ptJetGen[100];   //[nJetGen_akt5]
   Float_t         eJetGen[100];   //[nJetGen_akt5]
   Float_t         etaJetGen[100];   //[nJetGen_akt5]
   Float_t         phiJetGen[100];   //[nJetGen_akt5]
   Int_t           nMuonsGen[100];   //[nJetGen_akt5]
   Int_t           nElectronsGen[100];   //[nJetGen_akt5]
   Int_t           nPhotonsGen[100];   //[nJetGen_akt5]
   Int_t           nTracksGen[100];   //[nJetGen_akt5]
   Int_t           nNeutralHadronsGen[100];   //[nJetGen_akt5]
   Int_t           nHFHadronsGen[100];   //[nJetGen_akt5]
   Int_t           nHFEMGen[100];   //[nJetGen_akt5]
   Int_t           nNeutronsGen[100];   //[nJetGen_akt5]
   Int_t           nK0LGen[100];   //[nJetGen_akt5]
   Int_t           nK0SGen[100];   //[nJetGen_akt5]
   Int_t           nLambdasGen[100];   //[nJetGen_akt5]
   Int_t           nCsiGen[100];   //[nJetGen_akt5]
   Int_t           nOtherNeutralHadronsGen[100];   //[nJetGen_akt5]
   Float_t         eMuonsGen[100];   //[nJetGen_akt5]
   Float_t         eElectronsGen[100];   //[nJetGen_akt5]
   Float_t         ePhotonsGen[100];   //[nJetGen_akt5]
   Float_t         eTracksGen[100];   //[nJetGen_akt5]
   Float_t         eNeutralHadronsGen[100];   //[nJetGen_akt5]
   Float_t         eHFHadronsGen[100];   //[nJetGen_akt5]
   Float_t         eHFEMGen[100];   //[nJetGen_akt5]
   Float_t         eNeutronsGen[100];   //[nJetGen_akt5]
   Float_t         eK0LGen[100];   //[nJetGen_akt5]
   Float_t         eK0SGen[100];   //[nJetGen_akt5]
   Float_t         eLambdasGen[100];   //[nJetGen_akt5]
   Float_t         eCsiGen[100];   //[nJetGen_akt5]
   Float_t         eOtherNeutralHadronsGen[100];   //[nJetGen_akt5]
   Float_t         sMet  ;
   Float_t         eMet  ;
   Float_t         phiMet;
   Float_t         signifMet;
   Float_t         sCorrMet  ;
   Float_t         eCorrMet  ;
   Float_t         phiCorrMet;
   Float_t         signifCorrMet;
   Float_t         smuCorrMet  ;
   Float_t         emuCorrMet  ;
   Float_t         phimuCorrMet;
   Float_t         signifmuCorrMet;
   Float_t         sNoHFMet  ;
   Float_t         eNoHFMet  ;
   Float_t         phiNoHFMet;
   Float_t         signifNoHFMet;
   Float_t         stcMet  ;
   Float_t         etcMet  ;
   Float_t         phitcMet;
   Float_t         signiftcMet;
   Float_t         spfMet  ;
   Float_t         epfMet  ;
   Float_t         phipfMet;
   Float_t         signifpfMet;
   Float_t         sMetGen  ;
   Float_t         eMetGen  ;
   Float_t         phiMetGen;
   Float_t         signifMetGen;
   Float_t         sMetGen2  ;
   Float_t         eMetGen2  ;
   Float_t         phiMetGen2;
   Int_t           nvertex;
   Float_t         vxMC;
   Float_t         vyMC;
   Float_t         vzMC;
   Float_t         vx[100];   //[nvertex]
   Float_t         vy[100];   //[nvertex]
   Float_t         vz[100];   //[nvertex]
   Float_t         vntracks[100];   //[nvertex]
   Float_t         vchi2[100];   //[nvertex]
   Float_t         vndof[100];   //[nvertex]
   //Bool_t          hltPass;
   Int_t           nHLT;
   //Int_t           hltNamesLen;
   std::vector<std::string> *HLTNames;
   std::vector<bool>  *HLTResults;
 //Char_t          HLTNames[6000];   //[hltNamesLen]
 //Bool_t          HLTResults[50];   //[nHLT]
   Double_t        Xsec;


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
   TH1F* h1_lumi_;
   TH1F* h1_nCounter_;

   Int_t run_;
   Int_t LS_;
   Int_t nvertex_;
   Int_t oldrun_;
   Int_t currentLS_;
   Int_t event_;
   Float_t ptHat_;
   Float_t ptHatMin_;
   Float_t ptHatMax_;
   Float_t eventWeight_;
   Float_t rhoCalo_;
   Float_t rhoPF_;

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
   TBranch        *b_rhoCalo;   //!
   TBranch        *b_rhoPF;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_pdgIdMC;   //!
   TBranch        *b_statusMC;   //!
   TBranch        *b_motherIDMC;   //!
   TBranch        *b_ptMC ;   //!
   TBranch        *b_eMC  ;   //!
   TBranch        *b_etaMC;   //!
   TBranch        *b_phiMC;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_ptPhot ;   //!
   TBranch        *b_ePhot  ;   //!
   TBranch        *b_escPhot  ;   //!
   TBranch        *b_eseedPhot  ;   //!
   TBranch        *b_etaPhot;   //!
   TBranch        *b_phiPhot;   //!
   TBranch        *b_timePhot;   //!
   TBranch        *b_e4SwissCrossPhot;   //!
   TBranch        *b_hasPixelSeedPhot;   //!
   TBranch        *b_pid_isEM;   //!
   TBranch        *b_pid_isLoose;   //!
   TBranch        *b_pid_isTight;   //!
   TBranch        *b_pid_jurECAL;   //!
   TBranch        *b_pid_twrHCAL;   //!
   TBranch        *b_pid_HoverE;   //!
   TBranch        *b_pid_hlwTrack;   //!
   TBranch        *b_pid_hlwTrackNoDz;   //!
   TBranch        *b_pid_etawid;   //!
   TBranch        *b_ptiso035Phot;   //!
   TBranch        *b_ntrkiso035Phot;   //!
   TBranch        *b_ptiso04Phot;   //!
   TBranch        *b_ntrkiso04Phot;   //!
   TBranch        *b_hcalovecal04Phot;   //!
   TBranch        *b_ecaliso04Phot;   //!
   TBranch        *b_sMajMajPhot;   //!
   TBranch        *b_sMinMinPhot;   //!
   TBranch        *b_alphaPhot;   //!
   TBranch        *b_sEtaEtaPhot;   //!
   TBranch        *b_sEtaPhiPhot;   //!
   TBranch        *b_sPhiPhiPhot;   //!
   TBranch        *b_E1Phot;   //!
   TBranch        *b_E9Phot;   //!
   TBranch        *b_E25Phot;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_ptJet;   //!
   TBranch        *b_ptCorrJet;   //!
   TBranch        *b_eJet;   //!
   TBranch        *b_etaJet;   //!
   TBranch        *b_phiJet;   //!
   TBranch        *b_emfJet;   //!
   TBranch        *b_n90Jet;   //!
   TBranch        *b_n90HitsJet;   //!
   TBranch        *b_fHPDJet;   //!
   TBranch        *b_fRBXJet;   //!
   TBranch        *b_ptDJet;   //!
   TBranch        *b_rmsCandJet;   //!
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
   TBranch        *b_nNeutronsGen;   //!
   TBranch        *b_nK0LGen;   //!
   TBranch        *b_nK0SGen;   //!
   TBranch        *b_nLambdasGen;   //!
   TBranch        *b_nCsiGen;   //!
   TBranch        *b_nOtherNeutralHadronsGen;   //!
   TBranch        *b_eMuonsGen;   //!
   TBranch        *b_eElectronsGen;   //!
   TBranch        *b_ePhotonsGen;   //!
   TBranch        *b_eTracksGen;   //!
   TBranch        *b_eNeutralHadronsGen;   //!
   TBranch        *b_eHFHadronsGen;   //!
   TBranch        *b_eHFEMGen;   //!
   TBranch        *b_eNeutronsGen;   //!
   TBranch        *b_eK0LGen;   //!
   TBranch        *b_eK0SGen;   //!
   TBranch        *b_eLambdasGen;   //!
   TBranch        *b_eCsiGen;   //!
   TBranch        *b_eOtherNeutralHadronsGen;   //!
   TBranch        *b_sMet;   //!
   TBranch        *b_eMet;   //!
   TBranch        *b_phiMet;   //!
   TBranch        *b_signifMet;   //!
   TBranch        *b_sCorrMet;   //!
   TBranch        *b_eCorrMet;   //!
   TBranch        *b_phiCorrMet;   //!
   TBranch        *b_signifCorrMet;   //!
   TBranch        *b_smuCorrMet;   //!
   TBranch        *b_emuCorrMet;   //!
   TBranch        *b_phimuCorrMet;   //!
   TBranch        *b_signifmuCorrMet;   //!
   TBranch        *b_sNoHFMet;   //!
   TBranch        *b_eNoHFMet;   //!
   TBranch        *b_phiNoHFMet;   //!
   TBranch        *b_signifNoHFMet;   //!
   TBranch        *b_stcMet;   //!
   TBranch        *b_etcMet;   //!
   TBranch        *b_phitcMet;   //!
   TBranch        *b_signiftcMet;   //!
   TBranch        *b_spfMet;   //!
   TBranch        *b_epfMet;   //!
   TBranch        *b_phipfMet;   //!
   TBranch        *b_signifpfMet;   //!
   TBranch        *b_sMetGen;   //!
   TBranch        *b_eMetGen;   //!
   TBranch        *b_phiMetGen;   //!
   TBranch        *b_signifMetGen;   //!
   TBranch        *b_sMetGen2;   //!
   TBranch        *b_eMetGen2;   //!
   TBranch        *b_phiMetGen2;   //!
   TBranch        *b_nvertex;   //!
   TBranch        *b_vxMC;   //!
   TBranch        *b_vyMC;   //!
   TBranch        *b_vzMC;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vntracks;   //!
   TBranch        *b_vchi2;   //!
   TBranch        *b_vndof;   //!
   //TBranch        *b_hltPass;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_hltNamesLen;   //!
   TBranch        *b_HLTNames;   //!
   TBranch        *b_HLTResults;   //!
   TBranch        *b_Xsec;   //!

   TreeAnalyzer( const std::string& analyzerType, const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags="", TTree *tree=0);
   virtual ~TreeAnalyzer();

   virtual std::string GetFlags() { return flags_; };
   virtual void SetFlags( const std::string& flags ) { flags_ = flags; };

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     LoadInput();
   virtual void     LoadInputFromFile( const std::string& fileName );
   virtual void     CreateOutputFile();
   virtual void     Init(TTree *tree);
   virtual void     Loop()=0;
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     ReadJSONFile(const std::string& json);
   virtual void     ReadCSVFile(const std::string& csv);
   virtual void     UpdateCache();
   virtual bool     isGoodLS ();
   virtual bool     passedTrigger_regexp( const std::string& trigger );
   virtual GenEventParameters     getGenEventParameters ();

   //virtual void BookStuff()=0;

};

#endif
