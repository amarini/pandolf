#define TreeAnalyzer_MinBias_cxx
#include "TreeAnalyzer_MinBias.h"
#include <iostream>
#include "TMath.h"


class AnalysisPhoton {
 public:
  AnalysisPhoton() {
    pt = 0.;
    eta = 20.;
    phi = 0.;
  }

  Float_t e;
  Float_t pt;
  Float_t eta;
  Float_t phi;

  Float_t eGen;
  Float_t ptGen;
  Float_t etaGen;
  Float_t phiGen;

  Float_t hcalIso;
  Float_t ecalIso;
  Float_t nTrkIso;
  Float_t ptTrkIso;
  Float_t clusterMaj;
  Float_t clusterMin;

};


class AnalysisJet {
 public:
  AnalysisJet() {
    eReco = 0.;
    ptReco = 0.;
    etaReco = 20.;
    phiReco = 0.;
    eGen = 0.;
    ptGen = 0.;
    etaGen = 20.;
    phiGen = 0.;
  }

  Float_t eReco;
  Float_t ptReco;
  Float_t etaReco;
  Float_t phiReco;
  //Float_t eCorrReco;
  Float_t ptCorrReco;

  Float_t thetaReco() { return 2.*atan(exp(-etaReco)); };
  Float_t pReco() { return ptReco/sin(thetaReco()); };
  Float_t pxReco() { return ptReco*cos(phiReco); };
  Float_t pyReco() { return ptReco*sin(phiReco); };
  Float_t pzReco() { return pReco()*cos(thetaReco()); };

  Float_t eGen;
  Float_t ptGen;
  Float_t etaGen;
  Float_t phiGen;

  Float_t thetaGen() { return 2.*atan(exp(-etaGen)); };
  Float_t pGen() { return ptGen/sin(thetaGen()); };
  Float_t pxGen() { return ptGen*cos(phiGen); };
  Float_t pyGen() { return ptGen*sin(phiGen); };
  Float_t pzGen() { return pGen()*cos(thetaGen()); };

  Float_t emfReco;

  Float_t eTracksReco;
  Float_t ePhotonsReco;
  Float_t eNeutralHadronsReco;
  Float_t eMuonsReco;
  Float_t eElectronsReco;
  Float_t eHFHadronsReco;
  Float_t eHFEMReco;

  Int_t nTracksReco;
  Int_t nPhotonsReco;
  Int_t nNeutralHadronsReco;
  Int_t nMuonsReco;
  Int_t nElectronsReco;
  Int_t nHFHadronsReco;
  Int_t nHFEMReco;

  Float_t eTracksGen;
  Float_t ePhotonsGen;
  Float_t eNeutralHadronsGen;
  Float_t eMuonsGen;
  Float_t eElectronsGen;

  Int_t nTracksGen;
  Int_t nPhotonsGen;
  Int_t nNeutralHadronsGen;
  Int_t nMuonsGen;
  Int_t nElectronsGen;

};




void TreeAnalyzer_MinBias::Loop()
{

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 10000;
   else nentries = fChain->GetEntries();

   BookStuff();


   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%10000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;


     //vertex requirements: select only good collision events:
     if( vntracks < 3 ) continue;
     //if( vchi2 > 5 ) continue;
     if( fabs(vz) > 15. ) continue;
     
     if( epfMet > 0.5*spfMet ) continue;

     //eliminate events with a lot of tracks:
     //trackRatio_ = ((Float_t)nEventTracks/(Float_t)nEventChargedHadrons);
     //if( trackRatio_ > 5. ) continue;

     AnalysisJet firstJet;
     AnalysisJet secondJet;

     //look for first jet:
     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
       thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];

       thisJet.emfReco = (recoType_=="pf") ? 0. : emfJet[iRecoJet];

       thisJet.eTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eChargedHadrons[iRecoJet] : 0.;
       thisJet.ePhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ePhotons[iRecoJet] : 0.;
       thisJet.eNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eNeutralHadrons[iRecoJet] : 0.;
       thisJet.eMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eMuons[iRecoJet] : 0.;
       thisJet.eElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eElectrons[iRecoJet] : 0.;
       thisJet.eHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eHFHadrons[iRecoJet] : 0.;
       thisJet.eHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eHFEM[iRecoJet] : 0.;

       thisJet.nTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nChargedHadrons[iRecoJet] : 0;
       thisJet.nPhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nPhotons[iRecoJet] : 0;
       thisJet.nNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nNeutralHadrons[iRecoJet] : 0;
       thisJet.nMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nMuons[iRecoJet] : 0;
       thisJet.nElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nElectrons[iRecoJet] : 0;
       thisJet.nHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFHadrons[iRecoJet] : 0;
       thisJet.nHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFEM[iRecoJet] : 0;

       if( thisJet.ptReco > firstJet.ptReco )
         firstJet = thisJet;

     } //for reco jets

     if( firstJet.eReco == 0. ) continue;


     //look for second jet:
     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
       thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];
       

       thisJet.eTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eChargedHadrons[iRecoJet] : 0.;
       thisJet.ePhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ePhotons[iRecoJet] : 0.;
       thisJet.eNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eNeutralHadrons[iRecoJet] : 0.;
       thisJet.eMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eMuons[iRecoJet] : 0.;
       thisJet.eElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eElectrons[iRecoJet] : 0.;
       thisJet.eHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eHFHadrons[iRecoJet] : 0.;
       thisJet.eHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eHFEM[iRecoJet] : 0.;

       thisJet.nTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nChargedHadrons[iRecoJet] : 0;
       thisJet.nPhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nPhotons[iRecoJet] : 0;
       thisJet.nNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nNeutralHadrons[iRecoJet] : 0;
       thisJet.nMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nMuons[iRecoJet] : 0;
       thisJet.nElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nElectrons[iRecoJet] : 0;
       thisJet.nHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFHadrons[iRecoJet] : 0;
       thisJet.nHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFEM[iRecoJet] : 0;

       if( (thisJet.etaReco == firstJet.etaReco)&&(thisJet.phiReco==firstJet.phiReco)&&
           (thisJet.ptReco==firstJet.ptReco) ) continue;


       if( (thisJet.ptReco < firstJet.ptReco) && (thisJet.ptReco > secondJet.ptReco) )
         secondJet = thisJet;


     } //for reco jets



       eJetReco_  =  firstJet.eReco;
      ptJetReco_  =  firstJet.ptReco;
   //eCorrJetReco_  =  firstJet.eCorrReco;
  ptCorrJetReco_  =  firstJet.ptCorrReco;
     phiJetReco_  =  firstJet.phiReco;
     etaJetReco_  =  firstJet.etaReco;

     eTracksReco_= firstJet.eTracksReco;
     ePhotonsReco_= firstJet.ePhotonsReco;
     eNeutralHadronsReco_= firstJet.eNeutralHadronsReco;
     eMuonsReco_= firstJet.eMuonsReco;
     eElectronsReco_= firstJet.eElectronsReco;
     eHFHadronsReco_= firstJet.eHFHadronsReco;
     eHFEMReco_= firstJet.eHFEMReco;

     nTracksReco_= firstJet.nTracksReco;
     nPhotonsReco_= firstJet.nPhotonsReco;
     nNeutralHadronsReco_= firstJet.nNeutralHadronsReco;
     nMuonsReco_= firstJet.nMuonsReco;
     nElectronsReco_= firstJet.nElectronsReco;
     nHFHadronsReco_= firstJet.nHFHadronsReco;
     nHFEMReco_= firstJet.nHFEMReco;


     e2ndJetReco_= secondJet.eReco;
     pt2ndJetReco_= secondJet.ptReco;
     ptCorr2ndJetReco_= secondJet.ptCorrReco;
     eta2ndJetReco_= secondJet.etaReco;
     phi2ndJetReco_= secondJet.phiReco;

     eTracks2ndReco_= secondJet.eTracksReco;
     ePhotons2ndReco_= secondJet.ePhotonsReco;
     eNeutralHadrons2ndReco_= secondJet.eNeutralHadronsReco;
     eMuons2ndReco_= secondJet.eMuonsReco;
     eElectrons2ndReco_= secondJet.eElectronsReco;
     eHFHadrons2ndReco_= secondJet.eHFHadronsReco;
     eHFEM2ndReco_= secondJet.eHFEMReco;

     nTracks2ndReco_= secondJet.nTracksReco;
     nPhotons2ndReco_= secondJet.nPhotonsReco;
     nNeutralHadrons2ndReco_= secondJet.nNeutralHadronsReco;
     nMuons2ndReco_= secondJet.nMuonsReco;
     nElectrons2ndReco_= secondJet.nElectronsReco;
     nHFHadrons2ndReco_= secondJet.nHFHadronsReco;
     nHFEM2ndReco_= secondJet.nHFEMReco;

     pfMet_ = epfMet;
     sumEt_pfMet_ = spfMet;

     jetTree_->Fill(); 


   } //for entries

} //loop



void TreeAnalyzer_MinBias::BookStuff() {

  outfile_->cd();

  jetTree_ = new TTree("jetTree", "Reduced Tree for Jet Studies");


  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)


  jetTree_->Branch("trackRatio",  &trackRatio_,  "trackRatio_/F");

  jetTree_->Branch("ePhotReco",  &ePhotReco_,  "ePhotReco_/F");
  jetTree_->Branch("ptPhotReco",  &ptPhotReco_,  "ptPhotReco_/F");
  jetTree_->Branch("etaPhotReco",  &etaPhotReco_,  "etaPhotReco_/F");
  jetTree_->Branch("phiPhotReco",  &phiPhotReco_,  "phiPhotReco_/F");
  jetTree_->Branch("hcalIsoPhotReco",  &hcalIsoPhotReco_,  "hcalIsoPhotReco_/F");
  jetTree_->Branch("ecalIsoPhotReco",  &ecalIsoPhotReco_,  "ecalIsoPhotReco_/F");
  jetTree_->Branch("nTrkIsoPhotReco",  &nTrkIsoPhotReco_,  "nTrkIsoPhotReco_/I");
  jetTree_->Branch("ptTrkIsoPhotReco",  &ptTrkIsoPhotReco_,  "ptTrkIsoPhotReco_/F");
  jetTree_->Branch("clusterMajPhotReco",  &clusterMajPhotReco_,  "clusterMajPhotReco_/F");
  jetTree_->Branch("clusterMinPhotReco",  &clusterMinPhotReco_,  "clusterMinPhotReco_/F");

  jetTree_->Branch("eJetReco",  &eJetReco_,  "eJetReco_/F");
  jetTree_->Branch( "ptCorrJetReco",  &ptCorrJetReco_,  "ptCorrJetReco_/F");
  jetTree_->Branch("eCorrJetReco",  &eCorrJetReco_,  "eCorrJetReco_/F");
  jetTree_->Branch( "ptJetReco",  &ptJetReco_,  "ptJetReco_/F");
  jetTree_->Branch("etaJetReco", &etaJetReco_, "etaJetReco_/F");
  jetTree_->Branch("phiJetReco", &phiJetReco_, "phiJetReco_/F");

  jetTree_->Branch("pt2ndJetReco", &pt2ndJetReco_, "pt2ndJetReco_/F");
  jetTree_->Branch("e2ndJetReco", &e2ndJetReco_, "e2ndJetReco_/F");
  jetTree_->Branch("ptCorr2ndJetReco", &ptCorr2ndJetReco_, "ptCorr2ndJetReco_/F");
  jetTree_->Branch("eta2ndJetReco", &eta2ndJetReco_, "eta2ndJetReco_/F");
  jetTree_->Branch("phi2ndJetReco", &phi2ndJetReco_, "phi2ndJetReco_/F");

  jetTree_->Branch("eTracksReco", &eTracksReco_, "eTracksReco_/F");
  jetTree_->Branch("ePhotonsReco", &ePhotonsReco_, "ePhotonsReco_/F");
  jetTree_->Branch("eNeutralHadronsReco", &eNeutralHadronsReco_, "eNeutralHadronsReco_/F");
  jetTree_->Branch("eMuonsReco", &eMuonsReco_, "eMuonsReco_/F");
  jetTree_->Branch("eElectronsReco", &eElectronsReco_, "eElectronsReco_/F");
  jetTree_->Branch("eHFHadronsReco", &eHFHadronsReco_, "eHFHadronsReco_/F");
  jetTree_->Branch("eHFEMReco", &eHFEMReco_, "eHFEMReco_/F");

  jetTree_->Branch("nTracksReco", &nTracksReco_, "nTracksReco_/I");
  jetTree_->Branch("nPhotonsReco", &nPhotonsReco_, "nPhotonsReco_/I");
  jetTree_->Branch("nNeutralHadronsReco", &nNeutralHadronsReco_, "nNeutralHadronsReco_/I");
  jetTree_->Branch("nMuonsReco", &nMuonsReco_, "nMuonsReco_/I");
  jetTree_->Branch("nElectronsReco", &nElectronsReco_, "nElectronsReco_/I");
  jetTree_->Branch("nHFHadronsReco", &nHFHadronsReco_, "nHFHadronsReco_/I");
  jetTree_->Branch("nHFEMReco", &nHFEMReco_, "nHFEMReco_/I");

  jetTree_->Branch("eTracks2ndReco", &eTracks2ndReco_, "eTracks2ndReco_/F");
  jetTree_->Branch("ePhotons2ndReco", &ePhotons2ndReco_, "ePhotons2ndReco_/F");
  jetTree_->Branch("eNeutralHadrons2ndReco", &eNeutralHadrons2ndReco_, "eNeutralHadrons2ndReco_/F");
  jetTree_->Branch("eMuons2ndReco", &eMuons2ndReco_, "eMuons2ndReco_/F");
  jetTree_->Branch("eElectrons2ndReco", &eElectrons2ndReco_, "eElectrons2ndReco_/F");
  jetTree_->Branch("eHFHadrons2ndReco", &eHFHadrons2ndReco_, "eHFHadrons2ndReco_/F");
  jetTree_->Branch("eHFEM2ndReco", &eHFEM2ndReco_, "eHFEM2ndReco_/F");

  jetTree_->Branch("nTracks2ndReco", &nTracks2ndReco_, "nTracks2ndReco_/I");
  jetTree_->Branch("nPhotons2ndReco", &nPhotons2ndReco_, "nPhotons2ndReco_/I");
  jetTree_->Branch("nNeutralHadrons2ndReco", &nNeutralHadrons2ndReco_, "nNeutralHadrons2ndReco_/I");
  jetTree_->Branch("nMuons2ndReco", &nMuons2ndReco_, "nMuons2ndReco_/I");
  jetTree_->Branch("nElectrons2ndReco", &nElectrons2ndReco_, "nElectrons2ndReco_/I");
  jetTree_->Branch("nHFHadrons2ndReco", &nHFHadrons2ndReco_, "nHFHadrons2ndReco_/I");
  jetTree_->Branch("nHFEM2ndReco", &nHFEM2ndReco_, "nHFEM2ndReco_/I");

} //bookstuff


