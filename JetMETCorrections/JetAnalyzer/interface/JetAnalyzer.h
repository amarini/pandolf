// system include files
#include <memory>
#include <iostream>
#include <string>

// user include files
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


#include "TH1.h"
#include "TFile.h"
#include "TTree.h"


//
// class decleration
//

class JetAnalyzer : public edm::EDAnalyzer {
   public:
      explicit JetAnalyzer(const edm::ParameterSet&);
      ~JetAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      edm::InputTag MCTruthCollection_; 
      edm::InputTag trackTags_; 
      std::string filename_; 
      edm::Service<TFileService> fs_;

      // Tree with multiple info
      TH1F* h1_mfpK0S_itCone5;
      TH1F* h1_mfpLambda_itCone5;

      TTree * m_tree ;

      Int_t event;

      Float_t ptHat;

      Int_t nMC;
      Int_t pdgIdMC[80];
      Int_t statusMC[80];
      Float_t massMC[80];
      Int_t motherIDMC[80];
      Float_t pxMC[80];
      Float_t pyMC[80];
      Float_t pzMC[80];
      Float_t eMC[80];
      Float_t etaMC[80];
      Float_t phiMC[80];

      Int_t nPhot;
      Float_t pxPhot[40];
      Float_t pyPhot[40];
      Float_t pzPhot[40];
      Float_t ePhot[40];
      Float_t escPhot[40];
      Float_t etaPhot[40];
      Float_t phiPhot[40];

      Float_t ptiso0015Phot[40];
      Int_t ntrkiso0015Phot[40];
      Float_t ptiso035Phot[40];
      Int_t ntrkiso035Phot[40];
      Float_t ptiso05Phot[40];
      Int_t ntrkiso05Phot[40];
      Float_t ptiso07Phot[40];
      Int_t ntrkiso07Phot[40];
      Float_t ptiso1Phot[40];
      Int_t ntrkiso1Phot[40];
      Float_t ptisoatecal02Phot[40]; 
      Int_t ntrkisoatecal02Phot[40]; 
      Float_t ptisoatecal025Phot[40];
      Int_t ntrkisoatecal025Phot[40];
      Float_t ptisoatecal03Phot[40]; 
      Int_t ntrkisoatecal03Phot[40]; 
      Float_t ptisoatecal035Phot[40];
      Int_t ntrkisoatecal035Phot[40];
      Float_t ptisoatecal04Phot[40]; 
      Int_t ntrkisoatecal04Phot[40]; 
      Float_t hcalovecal01Phot[40];
      Float_t hcalovecal015Phot[40];
      Float_t hcalovecal04Phot[40]; 
      Float_t hcalovecal05Phot[40]; 
      Float_t hcalovecal07Phot[40];
      Float_t hcalovecal1Phot[40]; 
      Float_t ecaliso01Phot[40]; 
      Float_t ecaliso015Phot[40]; 
      Float_t ecaliso04Phot[40];  
      Float_t ecaliso05Phot[40];  
      Float_t ecaliso07Phot[40];  
      Float_t ecaliso1Phot[40];  
      Float_t ecaliso04oldPhot[40];  
      Float_t LATPhot[40];
      Float_t sMinMinPhot[40];
      Float_t sMajMajPhot[40];
      Float_t sMinMinoldPhot[40];
      Float_t sMajMajoldPhot[40];
      Float_t FisherPhot[40];
      Float_t alphaPhot[40];
      Float_t sEtaEtaPhot[40];
      Float_t sPhiPhiPhot[40];
      Float_t E1Phot[40];
      Float_t E9Phot[40];
      Float_t E25Phot[40];

      Int_t nJetGen_itCone5;
      Float_t pxJetGen_itCone5[40];
      Float_t pyJetGen_itCone5[40];
      Float_t pzJetGen_itCone5[40];
      Float_t eJetGen_itCone5[40];
      Float_t etaJetGen_itCone5[40];
      Float_t phiJetGen_itCone5[40];

      Int_t nTracksGen_itCone5[40];
      Int_t nPhotonsGen_itCone5[40];
      Int_t nElectronsGen_itCone5[40];
      Int_t nMuonsGen_itCone5[40];
      Int_t nNeutralHadronsGen_itCone5[40];

      Float_t eTracksGen_itCone5[40];
      Float_t ePhotonsGen_itCone5[40];
      Float_t eElectronsGen_itCone5[40];
      Float_t eMuonsGen_itCone5[40];
      Float_t eNeutralHadronsGen_itCone5[40];

      Float_t pxTracksGen_itCone5[40];
      Float_t pxPhotonsGen_itCone5[40];
      Float_t pxElectronsGen_itCone5[40];
      Float_t pxMuonsGen_itCone5[40];
      Float_t pxNeutralHadronsGen_itCone5[40];

      Float_t pyTracksGen_itCone5[40];
      Float_t pyPhotonsGen_itCone5[40];
      Float_t pyElectronsGen_itCone5[40];
      Float_t pyMuonsGen_itCone5[40];
      Float_t pyNeutralHadronsGen_itCone5[40];

      Float_t pzTracksGen_itCone5[40];
      Float_t pzPhotonsGen_itCone5[40];
      Float_t pzElectronsGen_itCone5[40];
      Float_t pzMuonsGen_itCone5[40];
      Float_t pzNeutralHadronsGen_itCone5[40];

      Float_t averagePartPtGen_itCone5[40];

/*
      Int_t nJetGen_sisCone5;
      Float_t pxJetGen_sisCone5[40];
      Float_t pyJetGen_sisCone5[40];
      Float_t pzJetGen_sisCone5[40];
      Float_t eJetGen_sisCone5[40];
      Float_t etaJetGen_sisCone5[40];
      Float_t phiJetGen_sisCone5[40];
   
      Int_t nTracksGen_sisCone5[40];
      Int_t nPhotonsGen_sisCone5[40];
      Int_t nElectronsGen_sisCone5[40];
      Int_t nMuonsGen_sisCone5[40];
      Int_t nNeutralHadronsGen_sisCone5[40];

      Float_t eTracksGen_sisCone5[40];
      Float_t ePhotonsGen_sisCone5[40];
      Float_t eElectronsGen_sisCone5[40];
      Float_t eMuonsGen_sisCone5[40];
      Float_t eNeutralHadronsGen_sisCone5[40];

      Float_t pxTracksGen_sisCone5[40];
      Float_t pxPhotonsGen_sisCone5[40];
      Float_t pxElectronsGen_sisCone5[40];
      Float_t pxMuonsGen_sisCone5[40];
      Float_t pxNeutralHadronsGen_sisCone5[40];

      Float_t pyTracksGen_sisCone5[40];
      Float_t pyPhotonsGen_sisCone5[40];
      Float_t pyElectronsGen_sisCone5[40];
      Float_t pyMuonsGen_sisCone5[40];
      Float_t pyNeutralHadronsGen_sisCone5[40];

      Float_t pzTracksGen_sisCone5[40];
      Float_t pzPhotonsGen_sisCone5[40];
      Float_t pzElectronsGen_sisCone5[40];
      Float_t pzMuonsGen_sisCone5[40];
      Float_t pzNeutralHadronsGen_sisCone5[40];

      Int_t nJetGen_kt4;
      Float_t pxJetGen_kt4[40];
      Float_t pyJetGen_kt4[40];
      Float_t pzJetGen_kt4[40];
      Float_t eJetGen_kt4[40];
      Float_t etaJetGen_kt4[40];
      Float_t phiJetGen_kt4[40];
      
      Int_t nTracksGen_kt4[40];
      Int_t nPhotonsGen_kt4[40];
      Int_t nElectronsGen_kt4[40];
      Int_t nMuonsGen_kt4[40];
      Int_t nNeutralHadronsGen_kt4[40];

      Float_t eTracksGen_kt4[40];
      Float_t ePhotonsGen_kt4[40];
      Float_t eElectronsGen_kt4[40];
      Float_t eMuonsGen_kt4[40];
      Float_t eNeutralHadronsGen_kt4[40];

      Float_t pxTracksGen_kt4[40];
      Float_t pxPhotonsGen_kt4[40];
      Float_t pxElectronsGen_kt4[40];
      Float_t pxMuonsGen_kt4[40];
      Float_t pxNeutralHadronsGen_kt4[40];

      Float_t pyTracksGen_kt4[40];
      Float_t pyPhotonsGen_kt4[40];
      Float_t pyElectronsGen_kt4[40];
      Float_t pyMuonsGen_kt4[40];
      Float_t pyNeutralHadronsGen_kt4[40];

      Float_t pzTracksGen_kt4[40];
      Float_t pzPhotonsGen_kt4[40];
      Float_t pzElectronsGen_kt4[40];
      Float_t pzMuonsGen_kt4[40];
      Float_t pzNeutralHadronsGen_kt4[40];
*/
      Int_t nJetCalo_itCone5;
      Float_t pxJetCalo_itCone5[40];
      Float_t pyJetCalo_itCone5[40];
      Float_t pzJetCalo_itCone5[40];
      Float_t eJetCalo_itCone5[40];
      Float_t etaJetCalo_itCone5[40];
      Float_t phiJetCalo_itCone5[40];

 /*     Int_t nJetCalo_sisCone5;
      Float_t pxJetCalo_sisCone5[40];
      Float_t pyJetCalo_sisCone5[40];
      Float_t pzJetCalo_sisCone5[40];
      Float_t eJetCalo_sisCone5[40];
      Float_t etaJetCalo_sisCone5[40];
      Float_t phiJetCalo_sisCone5[40];

      Int_t nJetCalo_kt4;
      Float_t pxJetCalo_kt4[40];
      Float_t pyJetCalo_kt4[40];
      Float_t pzJetCalo_kt4[40];
      Float_t eJetCalo_kt4[40];
      Float_t etaJetCalo_kt4[40];
      Float_t phiJetCalo_kt4[40];
*/
      Int_t nJetPF_itCone5;
      Float_t pxJetPF_itCone5[40];
      Float_t pyJetPF_itCone5[40];
      Float_t pzJetPF_itCone5[40];
      Float_t eJetPF_itCone5[40];
      Float_t etaJetPF_itCone5[40];
      Float_t phiJetPF_itCone5[40];

      Int_t nChargedHadronsPF_itCone5[40];
      Int_t nPhotonsPF_itCone5[40];
      Int_t nElectronsPF_itCone5[40];
      Int_t nMuonsPF_itCone5[40];
      Int_t nNeutralHadronsPF_itCone5[40];
      Int_t nHFHadronsPF_itCone5[40];
      Int_t nHFEMPF_itCone5[40];

      Float_t eChargedHadronsPF_itCone5[40];
      Float_t ePhotonsPF_itCone5[40];
      Float_t eElectronsPF_itCone5[40];
      Float_t eMuonsPF_itCone5[40];
      Float_t eNeutralHadronsPF_itCone5[40];
      Float_t eHFHadronsPF_itCone5[40];
      Float_t eHFEMPF_itCone5[40];

      Float_t pxChargedHadronsPF_itCone5[40];
      Float_t pxPhotonsPF_itCone5[40];
      Float_t pxElectronsPF_itCone5[40];
      Float_t pxMuonsPF_itCone5[40];
      Float_t pxNeutralHadronsPF_itCone5[40];
      Float_t pxHFHadronsPF_itCone5[40];
      Float_t pxHFEMPF_itCone5[40];

      Float_t pyChargedHadronsPF_itCone5[40];
      Float_t pyPhotonsPF_itCone5[40];
      Float_t pyElectronsPF_itCone5[40];
      Float_t pyMuonsPF_itCone5[40];
      Float_t pyNeutralHadronsPF_itCone5[40];
      Float_t pyHFHadronsPF_itCone5[40];
      Float_t pyHFEMPF_itCone5[40];

      Float_t pzChargedHadronsPF_itCone5[40];
      Float_t pzPhotonsPF_itCone5[40];
      Float_t pzElectronsPF_itCone5[40];
      Float_t pzMuonsPF_itCone5[40];
      Float_t pzNeutralHadronsPF_itCone5[40];
      Float_t pzHFHadronsPF_itCone5[40];
      Float_t pzHFEMPF_itCone5[40];

      Float_t averagePartPtPF_itCone5[40];

/*
      Int_t nJetPF_sisCone5;
      Float_t pxJetPF_sisCone5[40];
      Float_t pyJetPF_sisCone5[40];
      Float_t pzJetPF_sisCone5[40];
      Float_t eJetPF_sisCone5[40];
      Float_t etaJetPF_sisCone5[40];
      Float_t phiJetPF_sisCone5[40];

      Int_t nJetPF_kt4;
      Float_t pxJetPF_kt4[40];
      Float_t pyJetPF_kt4[40];
      Float_t pzJetPF_kt4[40];
      Float_t eJetPF_kt4[40];
      Float_t etaJetPF_kt4[40];
      Float_t phiJetPF_kt4[40];
*/

  
};

