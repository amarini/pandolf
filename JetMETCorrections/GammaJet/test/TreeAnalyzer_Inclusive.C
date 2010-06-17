#ifndef TreeAnalyzer_Inclusive_cxx
#define TreeAnalyzer_Inclusive_cxx
#include "TreeAnalyzer_Inclusive.h"
#include <iostream>
#include "TMath.h"
#include "AnalysisJet.h"



TreeAnalyzer_Inclusive::TreeAnalyzer_Inclusive(const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, TTree* tree) :
     TreeAnalyzer( "Inclusive", dataset, recoType, jetAlgo, (std::string)"", tree ) {


  DEBUG_VERBOSE_ = false;

  outfile_->cd();

  jetTree_ = new TTree("jetTree", "Reduced Tree for Jet Studies");


  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)

  jetTree_->Branch("nJetReco",  &nJetReco_,  "nJetReco_/I");
  jetTree_->Branch("eJetReco",  eJetReco_,  "eJetReco_[nJetReco_]/F");
  jetTree_->Branch( "ptCorrJetReco",  ptCorrJetReco_,  "ptCorrJetReco_[nJetReco_]/F");
  jetTree_->Branch("eCorrJetReco",  eCorrJetReco_,  "eCorrJetReco_[nJetReco_]/F");
  jetTree_->Branch( "ptJetReco",  ptJetReco_,  "ptJetReco_[nJetReco_]/F");
  jetTree_->Branch("etaJetReco", etaJetReco_, "etaJetReco_[nJetReco_]/F");
  jetTree_->Branch("phiJetReco", phiJetReco_, "phiJetReco_[nJetReco_]/F");

  jetTree_->Branch("eTracksReco", eTracksReco_, "eTracksReco_[nJetReco_]/F");
  jetTree_->Branch("ePhotonsReco", ePhotonsReco_, "ePhotonsReco_[nJetReco_]/F");
  jetTree_->Branch("eNeutralHadronsReco", eNeutralHadronsReco_, "eNeutralHadronsReco_[nJetReco_]/F");
  jetTree_->Branch("eMuonsReco", eMuonsReco_, "eMuonsReco_[nJetReco_]/F");
  jetTree_->Branch("eElectronsReco", eElectronsReco_, "eElectronsReco_[nJetReco_]/F");
  jetTree_->Branch("eHFHadronsReco", eHFHadronsReco_, "eHFHadronsReco_[nJetReco_]/F");
  jetTree_->Branch("eHFEMReco", eHFEMReco_, "eHFEMReco_[nJetReco_]/F");

  jetTree_->Branch("ptTracksReco", ptTracksReco_, "ptTracksReco_[nJetReco_]/F");
  jetTree_->Branch("ptPhotonsReco", ptPhotonsReco_, "ptPhotonsReco_[nJetReco_]/F");
  jetTree_->Branch("ptNeutralHadronsReco", ptNeutralHadronsReco_, "ptNeutralHadronsReco_[nJetReco_]/F");
  jetTree_->Branch("ptMuonsReco", ptMuonsReco_, "ptMuonsReco_[nJetReco_]/F");
  jetTree_->Branch("ptElectronsReco", ptElectronsReco_, "ptElectronsReco_[nJetReco_]/F");
  jetTree_->Branch("ptHFHadronsReco", ptHFHadronsReco_, "ptHFHadronsReco_[nJetReco_]/F");
  jetTree_->Branch("ptHFEMReco", ptHFEMReco_, "ptHFEMReco_[nJetReco_]/F");

  jetTree_->Branch("nTracksReco", nTracksReco_, "nTracksReco_[nJetReco_]/I");
  jetTree_->Branch("nPhotonsReco", nPhotonsReco_, "nPhotonsReco_[nJetReco_]/I");
  jetTree_->Branch("nNeutralHadronsReco", nNeutralHadronsReco_, "nNeutralHadronsReco_[nJetReco_]/I");
  jetTree_->Branch("nMuonsReco", nMuonsReco_, "nMuonsReco_[nJetReco_]/I");
  jetTree_->Branch("nElectronsReco", nElectronsReco_, "nElectronsReco_[nJetReco_]/I");
  jetTree_->Branch("nHFHadronsReco", nHFHadronsReco_, "nHFHadronsReco_[nJetReco_]/I");
  jetTree_->Branch("nHFEMReco", nHFEMReco_, "nHFEMReco_[nJetReco_]/I");

  jetTree_->Branch("nJetReco",  &nJetReco_,  "nJetReco_/I");
  jetTree_->Branch("eJetGen",  eJetGen_,  "eJetGen_[nJetReco_]/F");
  jetTree_->Branch( "ptJetGen",  ptJetGen_,  "ptJetGen_[nJetReco_]/F");
  jetTree_->Branch("etaJetGen", etaJetGen_, "etaJetGen_[nJetReco_]/F");
  jetTree_->Branch("phiJetGen", phiJetGen_, "phiJetGen_[nJetReco_]/F");

  jetTree_->Branch("eTracksGen", eTracksGen_, "eTracksGen_[nJetReco_]/F");
  jetTree_->Branch("ePhotonsGen", ePhotonsGen_, "ePhotonsGen_[nJetReco_]/F");
  jetTree_->Branch("eNeutralHadronsGen", eNeutralHadronsGen_, "eNeutralHadronsGen_[nJetReco_]/F");
  jetTree_->Branch("eMuonsGen", eMuonsGen_, "eMuonsGen_[nJetReco_]/F");
  jetTree_->Branch("eElectronsGen", eElectronsGen_, "eElectronsGen_[nJetReco_]/F");
  jetTree_->Branch("eHFHadronsGen", eHFHadronsGen_, "eHFHadronsGen_[nJetReco_]/F");
  jetTree_->Branch("eHFEMGen", eHFEMGen_, "eHFEMGen_[nJetReco_]/F");

  jetTree_->Branch("nTracksGen", nTracksGen_, "nTracksGen_[nJetReco_]/I");
  jetTree_->Branch("nPhotonsGen", nPhotonsGen_, "nPhotonsGen_[nJetReco_]/I");
  jetTree_->Branch("nNeutralHadronsGen", nNeutralHadronsGen_, "nNeutralHadronsGen_[nJetReco_]/I");
  jetTree_->Branch("nMuonsGen", nMuonsGen_, "nMuonsGen_[nJetReco_]/I");
  jetTree_->Branch("nElectronsGen", nElectronsGen_, "nElectronsGen_[nJetReco_]/I");
  jetTree_->Branch("nHFHadronsGen", nHFHadronsGen_, "nHFHadronsGen_[nJetReco_]/I");
  jetTree_->Branch("nHFEMGen", nHFEMGen_, "nHFEMGen_[nJetReco_]/I");

}


void TreeAnalyzer_Inclusive::Loop()
{

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 10000;
   else nentries = fChain->GetEntries();


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
     

     if( !isGoodLS() ) continue; //this takes care also of integrated luminosity


     nJetReco_ = 0;

     //for reco jets:
     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       if( iRecoJet >= 40 ) continue;

       AnalysisJet thisJet;

       eJetReco_[iRecoJet]  =  eJet[iRecoJet];
       ptJetReco_[iRecoJet]  =  ptJet[iRecoJet];
       ptCorrJetReco_[iRecoJet]  =  ptCorrJet[iRecoJet];
       phiJetReco_[iRecoJet] = phiJet[iRecoJet];
       etaJetReco_[iRecoJet] = etaJet[iRecoJet];

       emfJetReco_[iRecoJet] = (recoType_=="pf") ? 0. : emfJet[iRecoJet];

       eTracksReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eChargedHadrons[iRecoJet] : 0.;
       ePhotonsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ePhotons[iRecoJet] : 0.;
       eNeutralHadronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eNeutralHadrons[iRecoJet] : 0.;
       eMuonsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eMuons[iRecoJet] : 0.;
       eElectronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eElectrons[iRecoJet] : 0.;
       eHFHadronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eHFHadrons[iRecoJet] : 0.;
       eHFEMReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? eHFEM[iRecoJet] : 0.;

       ptTracksReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptChargedHadrons[iRecoJet] : 0.;
       ptPhotonsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptPhotons[iRecoJet] : 0.;
       ptNeutralHadronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptNeutralHadrons[iRecoJet] : 0.;
       ptMuonsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptMuons[iRecoJet] : 0.;
       ptElectronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptElectrons[iRecoJet] : 0.;
       ptHFHadronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptHFHadrons[iRecoJet] : 0.;
       ptHFEMReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptHFEM[iRecoJet] : 0.;

       nTracksReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nChargedHadrons[iRecoJet] : 0;
       nPhotonsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nPhotons[iRecoJet] : 0;
       nNeutralHadronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nNeutralHadrons[iRecoJet] : 0;
       nMuonsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nMuons[iRecoJet] : 0;
       nElectronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nElectrons[iRecoJet] : 0;
       nHFHadronsReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFHadrons[iRecoJet] : 0;
       nHFEMReco_[iRecoJet] = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFEM[iRecoJet] : 0;


       Float_t deltaR_max = 999.;
       AnalysisJet foundJet;

       for( unsigned iGenJet=0; iGenJet<nJetGen; ++iGenJet ) {

         AnalysisJet thisJet;

         thisJet.eGen = eJetGen[iGenJet];
         thisJet.ptGen = ptJetGen[iGenJet];
         thisJet.etaGen = etaJetGen[iGenJet];
         thisJet.phiGen = phiJetGen[iGenJet];

         thisJet.eTracksGen = (jetAlgo_=="akt5") ? eTracksGen[iGenJet] : 0.;
         thisJet.ePhotonsGen = (jetAlgo_=="akt5") ? ePhotonsGen[iGenJet] : 0.;
         thisJet.eNeutralHadronsGen = (jetAlgo_=="akt5") ? eNeutralHadronsGen[iGenJet] : 0.;
         thisJet.eMuonsGen = (jetAlgo_=="akt5") ? eMuonsGen[iGenJet] : 0.;
         thisJet.eElectronsGen = (jetAlgo_=="akt5") ? eElectronsGen[iGenJet] : 0.;
         thisJet.eHFHadronsGen = (jetAlgo_=="akt5") ? eHFHadronsGen[iGenJet] : 0.;
         thisJet.eHFEMGen = (jetAlgo_=="akt5") ? eHFEMGen[iGenJet] : 0.;
       
         thisJet.nTracksGen = (jetAlgo_=="akt5") ? nTracksGen[iGenJet] : 0;
         thisJet.nPhotonsGen = (jetAlgo_=="akt5") ? nPhotonsGen[iGenJet] : 0;
         thisJet.nNeutralHadronsGen = (jetAlgo_=="akt5") ? nNeutralHadronsGen[iGenJet] : 0;
         thisJet.nMuonsGen = (jetAlgo_=="akt5") ? nMuonsGen[iGenJet] : 0;
         thisJet.nElectronsGen = (jetAlgo_=="akt5") ? nElectronsGen[iGenJet] : 0;
         thisJet.nHFHadronsGen = (jetAlgo_=="akt5") ? nHFHadronsGen[iGenJet] : 0;
         thisJet.nHFEMGen = (jetAlgo_=="akt5") ? nHFEMGen[iGenJet] : 0;

         Float_t deltaEta = etaJetReco_[iRecoJet] - thisJet.etaGen;
         Float_t deltaPhi = phiJetReco_[iRecoJet] - thisJet.phiGen;
         if( deltaPhi > TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
         if( deltaPhi < -TMath::Pi()) deltaPhi += 2.*TMath::Pi();
         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
 
         if( deltaR < deltaR_max ) {
           deltaR_max = deltaR;
           foundJet = thisJet;
         }

       } // for gen jets

       eJetGen_[iRecoJet]  =  foundJet.eGen;
       ptJetGen_[iRecoJet]  =  foundJet.ptGen;
       phiJetGen_[iRecoJet] = foundJet.phiGen;
       etaJetGen_[iRecoJet] = foundJet.etaGen;

       eTracksGen_[iRecoJet] =  foundJet.eTracksGen;
       ePhotonsGen_[iRecoJet] =  foundJet.ePhotonsGen;
       eNeutralHadronsGen_[iRecoJet] =  foundJet.eNeutralHadronsGen;
       eMuonsGen_[iRecoJet] = foundJet.eMuonsGen;
       eElectronsGen_[iRecoJet] =  foundJet.eElectronsGen;
       eHFHadronsGen_[iRecoJet] =  foundJet.eHFHadronsGen;
       eHFEMGen_[iRecoJet] =  foundJet.eHFEMGen;

       nTracksGen_[iRecoJet] =  foundJet.nTracksGen;
       nPhotonsGen_[iRecoJet] =  foundJet.nPhotonsGen;
       nNeutralHadronsGen_[iRecoJet] =  foundJet.nNeutralHadronsGen;
       nMuonsGen_[iRecoJet] =  foundJet.nMuonsGen;
       nElectronsGen_[iRecoJet] =  foundJet.nElectronsGen;
       nHFHadronsGen_[iRecoJet] =  foundJet.nHFHadronsGen;
       nHFEMGen_[iRecoJet] =  foundJet.nHFEMGen;

       nJetReco_++;
       
     } //for reco jets

     pfMet_ = epfMet;
     sumEt_pfMet_ = spfMet;

     jetTree_->Fill(); 


   } //for entries

} //loop


#endif
