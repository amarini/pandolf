#include "TreeAnalyzer_HZZlljj.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "fitTools.C"


TreeAnalyzer_HZZlljj::TreeAnalyzer_HZZlljj( std::string dataset, std::string recoType, std::string jetAlgo, TTree* tree ) :
     TreeAnalyzer( "HZZlljj", dataset, recoType, jetAlgo, (std::string)"", tree ) {



  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)

  jetTree_->Branch("run",&run_,"run_/I");
  jetTree_->Branch("event",&event_,"event_/I");

  jetTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  jetTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");

  jetTree_->Branch("eZ",  &eZ_,  "eZ_/F");
  jetTree_->Branch( "ptZ",  &ptZ_,  "ptZ_/F");
  jetTree_->Branch("etaZ", &etaZ_, "etaZ_/F");
  jetTree_->Branch("phiZ", &phiZ_, "phiZ_/F");

  jetTree_->Branch("eJetReco",  &eJetReco_,  "eJetReco_/F");
  jetTree_->Branch( "ptCorrJetReco",  &ptCorrJetReco_,  "ptCorrJetReco_/F");
  jetTree_->Branch( "ptJetReco",  &ptJetReco_,  "ptJetReco_/F");
  jetTree_->Branch("etaJetReco", &etaJetReco_, "etaJetReco_/F");
  jetTree_->Branch("phiJetReco", &phiJetReco_, "phiJetReco_/F");
  jetTree_->Branch(  "eJetGen",   &eJetGen_,   "eJetGen_/F");
  jetTree_->Branch(  "ptJetGen",   &ptJetGen_,   "ptJetGen_/F");
  jetTree_->Branch( "etaJetGen",  &etaJetGen_,  "etaJetGen_/F");
  jetTree_->Branch( "phiJetGen",  &phiJetGen_,  "phiJetGen_/F");

  jetTree_->Branch("e2ndJetReco",  &e2ndJetReco_,  "e2ndJetReco_/F");
  jetTree_->Branch( "ptCorr2ndJetReco",  &ptCorr2ndJetReco_,  "ptCorr2ndJetReco_/F");
  jetTree_->Branch( "pt2ndJetReco",  &pt2ndJetReco_,  "pt2ndJetReco_/F");
  jetTree_->Branch("eta2ndJetReco", &eta2ndJetReco_, "eta2ndJetReco_/F");
  jetTree_->Branch("phi2ndJetReco", &phi2ndJetReco_, "phi2ndJetReco_/F");
  jetTree_->Branch(  "e2ndJetGen",   &e2ndJetGen_,   "e2ndJetGen_/F");
  jetTree_->Branch(  "pt2ndJetGen",   &pt2ndJetGen_,   "pt2ndJetGen_/F");
  jetTree_->Branch( "eta2ndJetGen",  &eta2ndJetGen_,  "eta2ndJetGen_/F");
  jetTree_->Branch( "phi2ndJetGen",  &phi2ndJetGen_,  "phi2ndJetGen_/F");

  jetTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  jetTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");

} //constructor






void TreeAnalyzer_HZZlljj::Loop()
{

   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;

     run_ = run;
     event_ = event;
     ptHat_ = genpt;

     if( !isGoodLS() ) continue; //this takes care also of integrated luminosity

     epfMet_ = epfMet;
     phipfMet_ = phipfMet;

     Float_t ptZ = 0.;
     Float_t  eZ = 0.;
     Float_t etaZ, phiZ;

     //first step: look for Z->ll (for now using MC info)
     for( unsigned iMC=0; iMC<nMC; ++iMC ) {
       if( pdgIdMC[iMC]==23 && statusMC[iMC]==3 ) {
          ptZ =  ptMC[iMC];
           eZ =   eMC[iMC];
         phiZ = phiMC[iMC];
         etaZ = etaMC[iMC];
       }
     } //for Z

     if( eZ==0. ) continue;

     ptZ_  = ptZ;
     eZ_   = eZ;
     phiZ_ = phiZ;
     etaZ_ = etaZ;

     if( nJet<2 ) continue;

     AnalysisJet firstJet;
     AnalysisJet secondJet;

     //look for first jet:
     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];
       thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];

       thisJet.emfReco = (recoType_=="pf") ? 0. : emfJet[iRecoJet];

       thisJet.eTracksReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eChargedHadrons[iRecoJet] : 0.;
       thisJet.ePhotonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? ePhotons[iRecoJet] : 0.;
       thisJet.eNeutralHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eNeutralHadrons[iRecoJet] : 0.;
       thisJet.eMuonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eMuons[iRecoJet] : 0.;
       thisJet.eElectronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eElectrons[iRecoJet] : 0.;
       thisJet.eHFHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eHFHadrons[iRecoJet] : 0.;
       thisJet.eHFEMReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eHFEM[iRecoJet] : 0.;

       thisJet.nTracksReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nChargedHadrons[iRecoJet] : 0;
       thisJet.nPhotonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nPhotons[iRecoJet] : 0;
       thisJet.nNeutralHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nNeutralHadrons[iRecoJet] : 0;
       thisJet.nMuonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nMuons[iRecoJet] : 0;
       thisJet.nElectronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nElectrons[iRecoJet] : 0;
       thisJet.nHFHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nHFHadrons[iRecoJet] : 0;
       thisJet.nHFEMReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nHFEM[iRecoJet] : 0;

       if( thisJet.ptCorrReco > firstJet.ptCorrReco )
         firstJet = thisJet;

     } //for reco jets

     if( firstJet.eCorrReco() == 0. ) continue;

     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];
       thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];

       if( (thisJet.etaReco == firstJet.etaReco)&&(thisJet.phiReco==firstJet.phiReco)&&
           (thisJet.ptReco==firstJet.ptReco) ) continue;


       if( (thisJet.ptCorrReco < firstJet.ptCorrReco) && (thisJet.ptCorrReco > secondJet.ptCorrReco) )
         secondJet = thisJet;


     } //for reco jets


       eJetReco_  =  firstJet.eReco;
      ptJetReco_  =  firstJet.ptReco;
  ptCorrJetReco_  =  firstJet.ptCorrReco;
     phiJetReco_  =  firstJet.phiReco;
     etaJetReco_  =  firstJet.etaReco;
        eJetGen_  =  firstJet.eGen;
       ptJetGen_  =  firstJet.ptGen;
      phiJetGen_  =  firstJet.phiGen;
      etaJetGen_  =  firstJet.etaGen;

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

     eTracksGen_= firstJet.eTracksGen;
     ePhotonsGen_= firstJet.ePhotonsGen;
     eNeutralHadronsGen_= firstJet.eNeutralHadronsGen;
     eMuonsGen_= firstJet.eMuonsGen;

     nTracksGen_= firstJet.nTracksGen;
     nPhotonsGen_= firstJet.nPhotonsGen;
     nNeutralHadronsGen_= firstJet.nNeutralHadronsGen;
     nMuonsGen_= firstJet.nMuonsGen;

       e2ndJetReco_  =  secondJet.eReco;
      pt2ndJetReco_  =  secondJet.ptReco;
  ptCorr2ndJetReco_  =  secondJet.ptCorrReco;
     phi2ndJetReco_  =  secondJet.phiReco;
     eta2ndJetReco_  =  secondJet.etaReco;
        e2ndJetGen_  =  secondJet.eGen;
       pt2ndJetGen_  =  secondJet.ptGen;
      phi2ndJetGen_  =  secondJet.phiGen;
      eta2ndJetGen_  =  secondJet.etaGen;

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

     eTracks2ndGen_= secondJet.eTracksGen;
     ePhotons2ndGen_= secondJet.ePhotonsGen;
     eNeutralHadrons2ndGen_= secondJet.eNeutralHadronsGen;
     eMuons2ndGen_= secondJet.eMuonsGen;

     nTracks2ndGen_= secondJet.nTracksGen;
     nPhotons2ndGen_= secondJet.nPhotonsGen;
     nNeutralHadrons2ndGen_= secondJet.nNeutralHadronsGen;
     nMuons2ndGen_= secondJet.nMuonsGen;



     jetTree_->Fill(); 

     

   } //for entries


} //loop



