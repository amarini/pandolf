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
  jetTree_->Branch("eTracksReco", &eTracksReco_, "eTracksReco_/F");
  jetTree_->Branch(  "eJetGen",   &eJetGen_,   "eJetGen_/F");
  jetTree_->Branch(  "ptJetGen",   &ptJetGen_,   "ptJetGen_/F");
  jetTree_->Branch( "etaJetGen",  &etaJetGen_,  "etaJetGen_/F");
  jetTree_->Branch( "phiJetGen",  &phiJetGen_,  "phiJetGen_/F");

  jetTree_->Branch("e2ndJetReco",  &e2ndJetReco_,  "e2ndJetReco_/F");
  jetTree_->Branch( "ptCorr2ndJetReco",  &ptCorr2ndJetReco_,  "ptCorr2ndJetReco_/F");
  jetTree_->Branch( "pt2ndJetReco",  &pt2ndJetReco_,  "pt2ndJetReco_/F");
  jetTree_->Branch("eta2ndJetReco", &eta2ndJetReco_, "eta2ndJetReco_/F");
  jetTree_->Branch("phi2ndJetReco", &phi2ndJetReco_, "phi2ndJetReco_/F");
  jetTree_->Branch("eTracks2ndReco", &eTracks2ndReco_, "eTracks2ndReco_/F");
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
     if( nJet<2 ) continue;


     std::vector<AnalysisJet> foundJets;

     int Z_ID=0; //this will be needed later when looking for Z->ll

     //look for Z jets:
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


       //look for matched genjet:

       Float_t deltaR_max = 999.;
     
       for( unsigned iGenJet=0; iGenJet<nJetGen; ++iGenJet ) {
     
         AnalysisJet otherJet;
     
         otherJet.eGen = eJetGen[iGenJet];
         otherJet.ptGen = ptJetGen[iGenJet];
         otherJet.etaGen = etaJetGen[iGenJet];
         otherJet.phiGen = phiJetGen[iGenJet];
     
         Float_t deltaEta = thisJet.etaReco - otherJet.etaGen;
         Float_t deltaPhi = fitTools::delta_phi(thisJet.phiReco, otherJet.phiGen);
         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
     
         if( deltaR < deltaR_max ) {
     
           deltaR_max = deltaR;
     
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
     
         } //if deltaR
     
       } // for gen jets

       //look for matched parton:
       Float_t deltaRMCmin = 999.;
       Int_t   motherID_found = 0;
       Int_t   iPartMC_found = 0;
     
       for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {
     
         if( statusMC[iPartMC]!=3 ) continue;
     
         Float_t eta = etaMC[iPartMC];
         Float_t phi = phiMC[iPartMC];
         Int_t   pdgId = pdgIdMC[iPartMC];
       
         Float_t deltaEtaMC = eta-thisJet.etaGen;
         Float_t deltaPhiMC = phi-thisJet.phiGen;
         Float_t pi = 3.14159;
         if( deltaPhiMC >=  pi ) deltaPhiMC -= 2.*pi;
         if( deltaPhiMC <= -pi ) deltaPhiMC += 2.*pi;
       
         Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );
     
         bool goodPdgId = false;
         if( (fabs(pdgId)<=9) || (fabs(pdgId)==21) ) goodPdgId = true;
       
         if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
           iPartMC_found = iPartMC;
           deltaRMCmin = deltaRMC;
           motherID_found = motherIDMC[iPartMC];
         }
     
       } //for MC particles


       if( deltaRMCmin < 0.25 && pdgIdMC[motherID_found]==23 ) {
         Z_ID= motherID_found;
         foundJets.push_back(thisJet);
       }

     } //for reco jets


     if( foundJets.size()<2 ) continue;



     //for now simply fill with the first two found jets:

       eJetReco_  =  foundJets[0].eReco;
      ptJetReco_  =  foundJets[0].ptReco;
  ptCorrJetReco_  =  foundJets[0].ptCorrReco;
     phiJetReco_  =  foundJets[0].phiReco;
     etaJetReco_  =  foundJets[0].etaReco;
     eTracksReco_  =  foundJets[0].eTracksReco;
        eJetGen_  =  foundJets[0].eGen;
       ptJetGen_  =  foundJets[0].ptGen;
      phiJetGen_  =  foundJets[0].phiGen;
      etaJetGen_  =  foundJets[0].etaGen;


     eTracksReco_        = foundJets[0].eTracksReco;
     ePhotonsReco_       = foundJets[0].ePhotonsReco;
     eNeutralHadronsReco_= foundJets[0].eNeutralHadronsReco;
     eMuonsReco_         = foundJets[0].eMuonsReco;
     eElectronsReco_     = foundJets[0].eElectronsReco;
     eHFHadronsReco_     = foundJets[0].eHFHadronsReco;
     eHFEMReco_          = foundJets[0].eHFEMReco;

     nTracksReco_        = foundJets[0].nTracksReco;
     nPhotonsReco_       = foundJets[0].nPhotonsReco;
     nNeutralHadronsReco_= foundJets[0].nNeutralHadronsReco;
     nMuonsReco_         = foundJets[0].nMuonsReco;
     nElectronsReco_     = foundJets[0].nElectronsReco;
     nHFHadronsReco_     = foundJets[0].nHFHadronsReco;
     nHFEMReco_          = foundJets[0].nHFEMReco;

     eTracksGen_        = foundJets[0].eTracksGen;
     ePhotonsGen_       = foundJets[0].ePhotonsGen;
     eNeutralHadronsGen_= foundJets[0].eNeutralHadronsGen;
     eMuonsGen_         = foundJets[0].eMuonsGen;

     nTracksGen_        = foundJets[0].nTracksGen;
     nPhotonsGen_       = foundJets[0].nPhotonsGen;
     nNeutralHadronsGen_= foundJets[0].nNeutralHadronsGen;
     nMuonsGen_         = foundJets[0].nMuonsGen;

       e2ndJetReco_  =  foundJets[1].eReco;
      pt2ndJetReco_  =  foundJets[1].ptReco;
  ptCorr2ndJetReco_  =  foundJets[1].ptCorrReco;
     phi2ndJetReco_  =  foundJets[1].phiReco;
     eta2ndJetReco_  =  foundJets[1].etaReco;
    eTracks2ndReco_  =  foundJets[1].eTracksReco;
        e2ndJetGen_  =  foundJets[1].eGen;
       pt2ndJetGen_  =  foundJets[1].ptGen;
      phi2ndJetGen_  =  foundJets[1].phiGen;
      eta2ndJetGen_  =  foundJets[1].etaGen;

     eTracks2ndReco_        = foundJets[1].eTracksReco;
     ePhotons2ndReco_       = foundJets[1].ePhotonsReco;
     eNeutralHadrons2ndReco_= foundJets[1].eNeutralHadronsReco;
     eMuons2ndReco_         = foundJets[1].eMuonsReco;
     eElectrons2ndReco_     = foundJets[1].eElectronsReco;
     eHFHadrons2ndReco_     = foundJets[1].eHFHadronsReco;
     eHFEM2ndReco_          = foundJets[1].eHFEMReco;

     nTracks2ndReco_        = foundJets[1].nTracksReco;
     nPhotons2ndReco_       = foundJets[1].nPhotonsReco;
     nNeutralHadrons2ndReco_= foundJets[1].nNeutralHadronsReco;
     nMuons2ndReco_         = foundJets[1].nMuonsReco;
     nElectrons2ndReco_     = foundJets[1].nElectronsReco;
     nHFHadrons2ndReco_     = foundJets[1].nHFHadronsReco;
     nHFEM2ndReco_          = foundJets[1].nHFEMReco;

     eTracks2ndGen_         = foundJets[1].eTracksGen;
     ePhotons2ndGen_        = foundJets[1].ePhotonsGen;
     eNeutralHadrons2ndGen_ = foundJets[1].eNeutralHadronsGen;
     eMuons2ndGen_          = foundJets[1].eMuonsGen;

     nTracks2ndGen_         = foundJets[1].nTracksGen;
     nPhotons2ndGen_        = foundJets[1].nPhotonsGen;
     nNeutralHadrons2ndGen_ = foundJets[1].nNeutralHadronsGen;
     nMuons2ndGen_          = foundJets[1].nMuonsGen;


    
     //now look for Z->ll (for now using MC info)

     Float_t ptZ = 0.;
     Float_t  eZ = 0.;
     Float_t etaZ, phiZ;


     for( unsigned iMC=0; iMC<nMC; ++iMC ) {
       if( pdgIdMC[iMC]==23 && statusMC[iMC]==3 && iMC!=Z_ID ) {
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



     jetTree_->Fill(); 

     

   } //for entries


} //loop



