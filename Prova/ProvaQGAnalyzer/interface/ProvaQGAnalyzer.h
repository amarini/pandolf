// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "RecoJets/JetProducers/interface/JetIDHelper.h"

#include "QuarkGluonAnalysis/QuarkGluonTagger2012/interface/QGLikelihoodCalculator.h"

//For HggVertexAnalysis


#include "TH1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

#include "TMVA/Reader.h"

#include <map>
#include <set>


using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//


class ProvaQGAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ProvaQGAnalyzer(const edm::ParameterSet&);
      ~ProvaQGAnalyzer();


   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      // Constants
      static const int kParton = 3;
      static const int kPhoton = 22;
      static const int kElectron = 11;


      // ----------member data ---------------------------
      bool _debug;

      edm::InputTag puSummaryInfo_;
      edm::InputTag JetPFsrcakt5_;
      string JetCorrector_pfakt5_; 
      double pfjetptthr_;
      double jptjetptthr_;
      int pfjetnmin_;
      double Xsec_;

//      edm::Service<TFileService> fs_;
      TFile* outfile;

      // Tree with multiple info
      TTree * m_tree ;

      reco::helper::JetIDHelper *jetID_;

      // Auxiliary event info will help to study correction stability
      // for different stores, as a function of instantaneous lumi (pile-up),
      // bunch-crossing (out-of-time pile-up), orbit number (beam heating) etc.
      Bool_t isMC;
      Int_t store;
      Int_t lbn;
      Int_t bx;
      Int_t orbit;
      Int_t run;
      Int_t event;

      Float_t rho;
      Int_t nJet_pfakt5;
      Float_t ptJet_pfakt5[100];
      Float_t ptCorrJet_pfakt5[100];
      Float_t eJet_pfakt5[100];
      Float_t etaJet_pfakt5[100];
      Float_t phiJet_pfakt5[100];
      Float_t ptDJet_pfakt5[100];
      Float_t rmsCandJet_pfakt5[100];
      Float_t rmsCandTrueJet_pfakt5[100];
      Float_t axis1Jet_pfakt5[100];
      Float_t axis2Jet_pfakt5[100];
      Float_t pullJet_pfakt5[100];
      Float_t tanaJet_pfakt5[100];
      Float_t ptD_QCJet_pfakt5[100];
      Float_t rmsCandTrue_QCJet_pfakt5[100];
      Float_t axis1_QCJet_pfakt5[100];
      Float_t axis2_QCJet_pfakt5[100];
      Float_t pull_QCJet_pfakt5[100];
      Float_t tana_QCJet_pfakt5[100];
      Float_t RchgJet_pfakt5[100];
      Float_t RneutralJet_pfakt5[100];
      Float_t RJet_pfakt5[100];
      Float_t Rchg_QCJet_pfakt5[100];
      Int_t   nChg_ptCutJet_pfakt5[100];
      Int_t   nChg_QCJet_pfakt5[100];
      Int_t   nChg_ptCut_QCJet_pfakt5[100];
      Int_t   nNeutral_ptCutJet_pfakt5[100];
      Int_t   nPFCand_QC_ptCutJet_pfakt5[100];
      Float_t pTMaxJet_pfakt5[100];
      Float_t pTMaxChgJet_pfakt5[100];
      Float_t pTMaxNeutralJet_pfakt5[100];
      Float_t pTMaxChg_QCJet_pfakt5[100];

      Float_t qglJet_pfakt5[100];

      Int_t   pdgIdPartJet_pfakt5[100];


      Float_t beta_pfakt5[100][100];
      Float_t betaStar_pfakt5[100][100];
      Float_t combinedSecondaryVertexBJetTags_pfakt5[100], 
              combinedSecondaryVertexMVABJetTags_pfakt5[100],
              jetBProbabilityBJetTags_pfakt5[100],
              jetProbabilityBJetTags_pfakt5[100],
              simpleSecondaryVertexHighEffBJetTags_pfakt5[100],
              simpleSecondaryVertexHighPurBJetTags_pfakt5[100],
              softMuonBJetTags_pfakt5[100],
              softMuonByIP3dBJetTags_pfakt5[100],
              softMuonByPtBJetTags_pfakt5[100],
              softElectronBJetTags_pfakt5[100],
              softElectronByIP3dBJetTags_pfakt5[100],
              softElectronByPtBJetTags_pfakt5[100],
              trackCountingHighPurBJetTags_pfakt5[100],
              trackCountingHighEffBJetTags_pfakt5[100];




      // Extra variables for PFlow studies
      Int_t nChargedHadrons_pfakt5[100];
      Int_t nPhotons_pfakt5[100];
      Int_t nElectrons_pfakt5[100];
      Int_t nMuons_pfakt5[100];
      Int_t nNeutralHadrons_pfakt5[100];
      Int_t nHFHadrons_pfakt5[100];
      Int_t nHFEM_pfakt5[100];

      Float_t eChargedHadrons_pfakt5[100];
      Float_t ePhotons_pfakt5[100];
      Float_t eElectrons_pfakt5[100];
      Float_t eMuons_pfakt5[100];
      Float_t eNeutralHadrons_pfakt5[100];
      Float_t eHFHadrons_pfakt5[100];
      Float_t eHFEM_pfakt5[100];



      QGLikelihoodCalculator* qglikeli;


      std::string outFileName;
};

