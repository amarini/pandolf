#ifndef TreeAnalyzer_DiJet_cxx
#define TreeAnalyzer_DiJet_cxx
#include "TreeAnalyzer_DiJet.h"
#include <iostream>
#include "TMath.h"
#include "AnalysisJet.h"
#include "fitTools.C"



TreeAnalyzer_DiJet::TreeAnalyzer_DiJet(const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, TTree* tree) :
     TreeAnalyzer( "DiJet", dataset, recoType, jetAlgo, (std::string)"", tree ) {


  DEBUG_VERBOSE_ = false;

  outfile_->cd();

  jetTree_ = new TTree("jetTree", "Reduced Tree for Jet Studies");


  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)

  jetTree_->Branch("run",&run_,"run_/I");
  jetTree_->Branch("event",&event_,"event_/I");
  jetTree_->Branch("lbn",&lbn_,"lbn_/I");

  jetTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  jetTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  jetTree_->Branch("eJetReco",  &eJetReco_,  "eJetReco_/F");
  jetTree_->Branch( "ptCorrJetReco",  &ptCorrJetReco_,  "ptCorrJetReco_/F");
  jetTree_->Branch("eCorrJetReco",  &eCorrJetReco_,  "eCorrJetReco_/F");
  jetTree_->Branch( "ptJetReco",  &ptJetReco_,  "ptJetReco_/F");
  jetTree_->Branch("etaJetReco", &etaJetReco_, "etaJetReco_/F");
  jetTree_->Branch("phiJetReco", &phiJetReco_, "phiJetReco_/F");

  jetTree_->Branch(  "eJetGen",   &eJetGen_,   "eJetGen_/F");
  jetTree_->Branch( "ptJetGen",  &ptJetGen_,  "ptJetGen_/F");
  jetTree_->Branch("etaJetGen", &etaJetGen_, "etaJetGen_/F");
  jetTree_->Branch("phiJetGen", &phiJetGen_, "phiJetGen_/F");

  jetTree_->Branch("pt2ndJetReco", &pt2ndJetReco_, "pt2ndJetReco_/F");
  jetTree_->Branch("e2ndJetReco", &e2ndJetReco_, "e2ndJetReco_/F");
  jetTree_->Branch("ptCorr2ndJetReco", &ptCorr2ndJetReco_, "ptCorr2ndJetReco_/F");
  jetTree_->Branch("eta2ndJetReco", &eta2ndJetReco_, "eta2ndJetReco_/F");
  jetTree_->Branch("phi2ndJetReco", &phi2ndJetReco_, "phi2ndJetReco_/F");

  jetTree_->Branch(  "e2ndJetGen",   &e2ndJetGen_,   "e2ndJetGen_/F");
  jetTree_->Branch( "pt2ndJetGen",  &pt2ndJetGen_,  "pt2ndJetGen_/F");
  jetTree_->Branch("eta2ndJetGen", &eta2ndJetGen_, "eta2ndJetGen_/F");
  jetTree_->Branch("phi2ndJetGen", &phi2ndJetGen_, "phi2ndJetGen_/F");

  jetTree_->Branch("eTracksReco", &eTracksReco_, "eTracksReco_/F");
  jetTree_->Branch("ePhotonsReco", &ePhotonsReco_, "ePhotonsReco_/F");
  jetTree_->Branch("eNeutralHadronsReco", &eNeutralHadronsReco_, "eNeutralHadronsReco_/F");
  jetTree_->Branch("eMuonsReco", &eMuonsReco_, "eMuonsReco_/F");
  jetTree_->Branch("eElectronsReco", &eElectronsReco_, "eElectronsReco_/F");
  jetTree_->Branch("eHFHadronsReco", &eHFHadronsReco_, "eHFHadronsReco_/F");
  jetTree_->Branch("eHFEMReco", &eHFEMReco_, "eHFEMReco_/F");

  jetTree_->Branch("ptTracksReco", &ptTracksReco_, "ptTracksReco_/F");
  jetTree_->Branch("ptPhotonsReco", &ptPhotonsReco_, "ptPhotonsReco_/F");
  jetTree_->Branch("ptNeutralHadronsReco", &ptNeutralHadronsReco_, "ptNeutralHadronsReco_/F");
  jetTree_->Branch("ptMuonsReco", &ptMuonsReco_, "ptMuonsReco_/F");
  jetTree_->Branch("ptElectronsReco", &ptElectronsReco_, "ptElectronsReco_/F");
  jetTree_->Branch("ptHFHadronsReco", &ptHFHadronsReco_, "ptHFHadronsReco_/F");
  jetTree_->Branch("ptHFEMReco", &ptHFEMReco_, "ptHFEMReco_/F");

  jetTree_->Branch("nTracksReco", &nTracksReco_, "nTracksReco_/I");
  jetTree_->Branch("nPhotonsReco", &nPhotonsReco_, "nPhotonsReco_/I");
  jetTree_->Branch("nNeutralHadronsReco", &nNeutralHadronsReco_, "nNeutralHadronsReco_/I");
  jetTree_->Branch("nMuonsReco", &nMuonsReco_, "nMuonsReco_/I");
  jetTree_->Branch("nElectronsReco", &nElectronsReco_, "nElectronsReco_/I");
  jetTree_->Branch("nHFHadronsReco", &nHFHadronsReco_, "nHFHadronsReco_/I");
  jetTree_->Branch("nHFEMReco", &nHFEMReco_, "nHFEMReco_/I");

  jetTree_->Branch("eTracksGen", &eTracksGen_, "eTracksGen_/F");
  jetTree_->Branch("ePhotonsGen", &ePhotonsGen_, "ePhotonsGen_/F");
  jetTree_->Branch("eNeutralHadronsGen", &eNeutralHadronsGen_, "eNeutralHadronsGen_/F");
  jetTree_->Branch("eMuonsGen", &eMuonsGen_, "eMuonsGen_/F");
  jetTree_->Branch("eElectronsGen", &eElectronsGen_, "eElectronsGen_/F");
  jetTree_->Branch("eHFHadronsGen", &eHFHadronsGen_, "eHFHadronsGen_/F");
  jetTree_->Branch("eHFEMGen", &eHFEMGen_, "eHFEMGen_/F");

  jetTree_->Branch("nTracksGen", &nTracksGen_, "nTracksGen_/I");
  jetTree_->Branch("nPhotonsGen", &nPhotonsGen_, "nPhotonsGen_/I");
  jetTree_->Branch("nNeutralHadronsGen", &nNeutralHadronsGen_, "nNeutralHadronsGen_/I");
  jetTree_->Branch("nMuonsGen", &nMuonsGen_, "nMuonsGen_/I");
  jetTree_->Branch("nElectronsGen", &nElectronsGen_, "nElectronsGen_/I");
  jetTree_->Branch("nHFHadronsGen", &nHFHadronsGen_, "nHFHadronsGen_/I");
  jetTree_->Branch("nHFEMGen", &nHFEMGen_, "nHFEMGen_/I");

  jetTree_->Branch("eTracks2ndReco", &eTracks2ndReco_, "eTracks2ndReco_/F");
  jetTree_->Branch("ePhotons2ndReco", &ePhotons2ndReco_, "ePhotons2ndReco_/F");
  jetTree_->Branch("eNeutralHadrons2ndReco", &eNeutralHadrons2ndReco_, "eNeutralHadrons2ndReco_/F");
  jetTree_->Branch("eMuons2ndReco", &eMuons2ndReco_, "eMuons2ndReco_/F");
  jetTree_->Branch("eElectrons2ndReco", &eElectrons2ndReco_, "eElectrons2ndReco_/F");
  jetTree_->Branch("eHFHadrons2ndReco", &eHFHadrons2ndReco_, "eHFHadrons2ndReco_/F");
  jetTree_->Branch("eHFEM2ndReco", &eHFEM2ndReco_, "eHFEM2ndReco_/F");

  jetTree_->Branch("ptTracks2ndReco", &ptTracks2ndReco_, "eTracks2ndReco_/F");
  jetTree_->Branch("ptPhotons2ndReco", &ptPhotons2ndReco_, "ptPhotons2ndReco_/F");
  jetTree_->Branch("ptNeutralHadrons2ndReco", &ptNeutralHadrons2ndReco_, "ptNeutralHadrons2ndReco_/F");
  jetTree_->Branch("ptMuons2ndReco", &ptMuons2ndReco_, "ptMuons2ndReco_/F");
  jetTree_->Branch("ptElectrons2ndReco", &ptElectrons2ndReco_, "ptElectrons2ndReco_/F");
  jetTree_->Branch("ptHFHadrons2ndReco", &ptHFHadrons2ndReco_, "ptHFHadrons2ndReco_/F");
  jetTree_->Branch("ptHFEM2ndReco", &ptHFEM2ndReco_, "ptHFEM2ndReco_/F");

  jetTree_->Branch("nTracks2ndReco", &nTracks2ndReco_, "nTracks2ndReco_/I");
  jetTree_->Branch("nPhotons2ndReco", &nPhotons2ndReco_, "nPhotons2ndReco_/I");
  jetTree_->Branch("nNeutralHadrons2ndReco", &nNeutralHadrons2ndReco_, "nNeutralHadrons2ndReco_/I");
  jetTree_->Branch("nMuons2ndReco", &nMuons2ndReco_, "nMuons2ndReco_/I");
  jetTree_->Branch("nElectrons2ndReco", &nElectrons2ndReco_, "nElectrons2ndReco_/I");
  jetTree_->Branch("nHFHadrons2ndReco", &nHFHadrons2ndReco_, "nHFHadrons2ndReco_/I");
  jetTree_->Branch("nHFEM2ndReco", &nHFEM2ndReco_, "nHFEM2ndReco_/I");

  jetTree_->Branch("eTracks2ndGen", &eTracks2ndGen_, "eTracks2ndGen_/F");
  jetTree_->Branch("ePhotons2ndGen", &ePhotons2ndGen_, "ePhotons2ndGen_/F");
  jetTree_->Branch("eNeutralHadrons2ndGen", &eNeutralHadrons2ndGen_, "eNeutralHadrons2ndGen_/F");
  jetTree_->Branch("eMuons2ndGen", &eMuons2ndGen_, "eMuons2ndGen_/F");
  jetTree_->Branch("eElectrons2ndGen", &eElectrons2ndGen_, "eElectrons2ndGen_/F");
  jetTree_->Branch("eHFHadrons2ndGen", &eHFHadrons2ndGen_, "eHFHadrons2ndGen_/F");
  jetTree_->Branch("eHFEM2ndGen", &eHFEM2ndGen_, "eHFEM2ndGen_/F");

  jetTree_->Branch("nTracks2ndGen", &nTracks2ndGen_, "nTracks2ndGen_/I");
  jetTree_->Branch("nPhotons2ndGen", &nPhotons2ndGen_, "nPhotons2ndGen_/I");
  jetTree_->Branch("nNeutralHadrons2ndGen", &nNeutralHadrons2ndGen_, "nNeutralHadrons2ndGen_/I");
  jetTree_->Branch("nMuons2ndGen", &nMuons2ndGen_, "nMuons2ndGen_/I");
  jetTree_->Branch("nElectrons2ndGen", &nElectrons2ndGen_, "nElectrons2ndGen_/I");
  jetTree_->Branch("nHFHadrons2ndGen", &nHFHadrons2ndGen_, "nHFHadrons2ndGen_/I");
  jetTree_->Branch("nHFEM2ndGen", &nHFEM2ndGen_, "nHFEM2ndGen_/I");

}



void TreeAnalyzer_DiJet::Loop()
{

   if (fChain == 0) return;

//  DEBUG_ = true;

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


   ////vertex requirements: select only good collision events:
   //if( vntracks < 3 ) continue;
   ////if( vchi2 > 5 ) continue;
   //if( fabs(vz) > 15. ) continue;
     
     run_ = run;
     event_ = event;
     lbn_ = lbn;
     ptHat_ = genpt;

     if( !isGoodLS() ) continue; //this takes care also of integrated luminosity

     if( isMC )
       if( ptHat_ > ptHatMax_ ) continue;


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

       thisJet.ptTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptChargedHadrons[iRecoJet] : 0.;
       thisJet.ptPhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptPhotons[iRecoJet] : 0.;
       thisJet.ptNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptNeutralHadrons[iRecoJet] : 0.;
       thisJet.ptMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptMuons[iRecoJet] : 0.;
       thisJet.ptElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptElectrons[iRecoJet] : 0.;
       thisJet.ptHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptHFHadrons[iRecoJet] : 0.;
       thisJet.ptHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptHFEM[iRecoJet] : 0.;

       thisJet.nTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nChargedHadrons[iRecoJet] : 0;
       thisJet.nPhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nPhotons[iRecoJet] : 0;
       thisJet.nNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nNeutralHadrons[iRecoJet] : 0;
       thisJet.nMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nMuons[iRecoJet] : 0;
       thisJet.nElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nElectrons[iRecoJet] : 0;
       thisJet.nHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFHadrons[iRecoJet] : 0;
       thisJet.nHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? nHFEM[iRecoJet] : 0;

       if( thisJet.ptReco > firstJet.ptReco )
         firstJet = thisJet;

     }

     if( firstJet.eReco == 0. ) continue;



     Float_t deltaR_max = 999.;

     for( unsigned iGenJet=0; iGenJet<nJetGen; ++iGenJet ) {

       AnalysisJet otherJet;

       otherJet.eGen = eJetGen[iGenJet];
       otherJet.ptGen = ptJetGen[iGenJet];
       otherJet.etaGen = etaJetGen[iGenJet];
       otherJet.phiGen = phiJetGen[iGenJet];

       Float_t deltaEta = firstJet.etaReco - otherJet.etaGen;
       Float_t deltaPhi = fitTools::delta_phi(firstJet.phiReco, otherJet.phiGen);
       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR < deltaR_max ) {

         deltaR_max = deltaR;

         firstJet.eGen = eJetGen[iGenJet];
         firstJet.ptGen = ptJetGen[iGenJet];
         firstJet.etaGen = etaJetGen[iGenJet];
         firstJet.phiGen = phiJetGen[iGenJet];

         firstJet.eTracksGen = (jetAlgo_=="akt5") ? eTracksGen[iGenJet] : 0.;
         firstJet.ePhotonsGen = (jetAlgo_=="akt5") ? ePhotonsGen[iGenJet] : 0.;
         firstJet.eNeutralHadronsGen = (jetAlgo_=="akt5") ? eNeutralHadronsGen[iGenJet] : 0.;
         firstJet.eMuonsGen = (jetAlgo_=="akt5") ? eMuonsGen[iGenJet] : 0.;
         firstJet.eElectronsGen = (jetAlgo_=="akt5") ? eElectronsGen[iGenJet] : 0.;
         firstJet.eHFHadronsGen = (jetAlgo_=="akt5") ? eHFHadronsGen[iGenJet] : 0.;
         firstJet.eHFEMGen = (jetAlgo_=="akt5") ? eHFEMGen[iGenJet] : 0.;
       
         firstJet.nTracksGen = (jetAlgo_=="akt5") ? nTracksGen[iGenJet] : 0;
         firstJet.nPhotonsGen = (jetAlgo_=="akt5") ? nPhotonsGen[iGenJet] : 0;
         firstJet.nNeutralHadronsGen = (jetAlgo_=="akt5") ? nNeutralHadronsGen[iGenJet] : 0;
         firstJet.nMuonsGen = (jetAlgo_=="akt5") ? nMuonsGen[iGenJet] : 0;
         firstJet.nElectronsGen = (jetAlgo_=="akt5") ? nElectronsGen[iGenJet] : 0;
         firstJet.nHFHadronsGen = (jetAlgo_=="akt5") ? nHFHadronsGen[iGenJet] : 0;
         firstJet.nHFEMGen = (jetAlgo_=="akt5") ? nHFEMGen[iGenJet] : 0;

       } //if deltaR

     } // for gen jets


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

       thisJet.ptTracksReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptChargedHadrons[iRecoJet] : 0.;
       thisJet.ptPhotonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptPhotons[iRecoJet] : 0.;
       thisJet.ptNeutralHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptNeutralHadrons[iRecoJet] : 0.;
       thisJet.ptMuonsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptMuons[iRecoJet] : 0.;
       thisJet.ptElectronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptElectrons[iRecoJet] : 0.;
       thisJet.ptHFHadronsReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptHFHadrons[iRecoJet] : 0.;
       thisJet.ptHFEMReco = (recoType_=="pf" && (jetAlgo_=="akt5"||jetAlgo_=="akt7")) ? ptHFEM[iRecoJet] : 0.;

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


     for( unsigned iGenJet=0; iGenJet<nJetGen; ++iGenJet ) {

       AnalysisJet otherJet;

       otherJet.eGen = eJetGen[iGenJet];
       otherJet.ptGen = ptJetGen[iGenJet];
       otherJet.etaGen = etaJetGen[iGenJet];
       otherJet.phiGen = phiJetGen[iGenJet];

       Float_t deltaEta = secondJet.etaReco - otherJet.etaGen;
       Float_t deltaPhi = fitTools::delta_phi(secondJet.phiReco, otherJet.phiGen);
       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR < deltaR_max ) {

         deltaR_max = deltaR;

         secondJet.eGen = eJetGen[iGenJet];
         secondJet.ptGen = ptJetGen[iGenJet];
         secondJet.etaGen = etaJetGen[iGenJet];
         secondJet.phiGen = phiJetGen[iGenJet];

         secondJet.eTracksGen = (jetAlgo_=="akt5") ? eTracksGen[iGenJet] : 0.;
         secondJet.ePhotonsGen = (jetAlgo_=="akt5") ? ePhotonsGen[iGenJet] : 0.;
         secondJet.eNeutralHadronsGen = (jetAlgo_=="akt5") ? eNeutralHadronsGen[iGenJet] : 0.;
         secondJet.eMuonsGen = (jetAlgo_=="akt5") ? eMuonsGen[iGenJet] : 0.;
         secondJet.eElectronsGen = (jetAlgo_=="akt5") ? eElectronsGen[iGenJet] : 0.;
         secondJet.eHFHadronsGen = (jetAlgo_=="akt5") ? eHFHadronsGen[iGenJet] : 0.;
         secondJet.eHFEMGen = (jetAlgo_=="akt5") ? eHFEMGen[iGenJet] : 0.;
       
         secondJet.nTracksGen = (jetAlgo_=="akt5") ? nTracksGen[iGenJet] : 0;
         secondJet.nPhotonsGen = (jetAlgo_=="akt5") ? nPhotonsGen[iGenJet] : 0;
         secondJet.nNeutralHadronsGen = (jetAlgo_=="akt5") ? nNeutralHadronsGen[iGenJet] : 0;
         secondJet.nMuonsGen = (jetAlgo_=="akt5") ? nMuonsGen[iGenJet] : 0;
         secondJet.nElectronsGen = (jetAlgo_=="akt5") ? nElectronsGen[iGenJet] : 0;
         secondJet.nHFHadronsGen = (jetAlgo_=="akt5") ? nHFHadronsGen[iGenJet] : 0;
         secondJet.nHFEMGen = (jetAlgo_=="akt5") ? nHFEMGen[iGenJet] : 0;

       } //if deltaR

     } // for gen jets


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

     ptTracksReco_=         firstJet.ptTracksReco;
     ptPhotonsReco_=        firstJet.ptPhotonsReco;
     ptNeutralHadronsReco_= firstJet.ptNeutralHadronsReco;
     ptMuonsReco_=          firstJet.ptMuonsReco;
     ptElectronsReco_=      firstJet.ptElectronsReco;
     ptHFHadronsReco_=      firstJet.ptHFHadronsReco;
     ptHFEMReco_=           firstJet.ptHFEMReco;

     nTracksReco_= firstJet.nTracksReco;
     nPhotonsReco_= firstJet.nPhotonsReco;
     nNeutralHadronsReco_= firstJet.nNeutralHadronsReco;
     nMuonsReco_= firstJet.nMuonsReco;
     nElectronsReco_= firstJet.nElectronsReco;
     nHFHadronsReco_= firstJet.nHFHadronsReco;
     nHFEMReco_= firstJet.nHFEMReco;

     eJetGen_  =  firstJet.eGen;
     ptJetGen_  =  firstJet.ptGen;
     phiJetGen_ = firstJet.phiGen;
     etaJetGen_ = firstJet.etaGen;

     eTracksGen_ =  firstJet.eTracksGen;
     ePhotonsGen_ =  firstJet.ePhotonsGen;
     eNeutralHadronsGen_ =  firstJet.eNeutralHadronsGen;
     eMuonsGen_ = firstJet.eMuonsGen;
     eElectronsGen_ =  firstJet.eElectronsGen;
     eHFHadronsGen_ =  firstJet.eHFHadronsGen;
     eHFEMGen_ =  firstJet.eHFEMGen;

     nTracksGen_ =  firstJet.nTracksGen;
     nPhotonsGen_ =  firstJet.nPhotonsGen;
     nNeutralHadronsGen_ =  firstJet.nNeutralHadronsGen;
     nMuonsGen_ =  firstJet.nMuonsGen;
     nElectronsGen_ =  firstJet.nElectronsGen;
     nHFHadronsGen_ =  firstJet.nHFHadronsGen;
     nHFEMGen_ =  firstJet.nHFEMGen;


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

     ptTracks2ndReco_=         secondJet.ptTracksReco;
     ptPhotons2ndReco_=        secondJet.ptPhotonsReco;
     ptNeutralHadrons2ndReco_= secondJet.ptNeutralHadronsReco;
     ptMuons2ndReco_=          secondJet.ptMuonsReco;
     ptElectrons2ndReco_=      secondJet.ptElectronsReco;
     ptHFHadrons2ndReco_=      secondJet.ptHFHadronsReco;
     ptHFEM2ndReco_=           secondJet.ptHFEMReco;

     nTracks2ndReco_= secondJet.nTracksReco;
     nPhotons2ndReco_= secondJet.nPhotonsReco;
     nNeutralHadrons2ndReco_= secondJet.nNeutralHadronsReco;
     nMuons2ndReco_= secondJet.nMuonsReco;
     nElectrons2ndReco_= secondJet.nElectronsReco;
     nHFHadrons2ndReco_= secondJet.nHFHadronsReco;
     nHFEM2ndReco_= secondJet.nHFEMReco;

     e2ndJetGen_  =  secondJet.eGen;
     pt2ndJetGen_  =  secondJet.ptGen;
     phi2ndJetGen_ = secondJet.phiGen;
     eta2ndJetGen_ = secondJet.etaGen;

     eTracks2ndGen_ =  secondJet.eTracksGen;
     ePhotons2ndGen_ =  secondJet.ePhotonsGen;
     eNeutralHadrons2ndGen_ =  secondJet.eNeutralHadronsGen;
     eMuons2ndGen_ = secondJet.eMuonsGen;
     eElectrons2ndGen_ =  secondJet.eElectronsGen;
     eHFHadrons2ndGen_ =  secondJet.eHFHadronsGen;
     eHFEM2ndGen_ =  secondJet.eHFEMGen;

     nTracks2ndGen_ =  secondJet.nTracksGen;
     nPhotons2ndGen_ =  secondJet.nPhotonsGen;
     nNeutralHadrons2ndGen_ =  secondJet.nNeutralHadronsGen;
     nMuons2ndGen_ =  secondJet.nMuonsGen;
     nElectrons2ndGen_ =  secondJet.nElectronsGen;
     nHFHadrons2ndGen_ =  secondJet.nHFHadronsGen;
     nHFEM2ndGen_ =  secondJet.nHFEMGen;


     pfMet_ = epfMet;
     sumEt_pfMet_ = spfMet;

     jetTree_->Fill(); 


   } //for entries

} //loop



#endif
