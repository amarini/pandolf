#include "TreeAnalyzer_VHgg.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisPhoton.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "TLorentzVector.h"





TreeAnalyzer_VHgg::TreeAnalyzer_VHgg( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags, bool useGenJets, TTree* tree ) :
     TreeAnalyzer( "VHgg", dataset, recoType, jetAlgo, flags, tree ) {


  useGenJets_=useGenJets;

} //constructor



void TreeAnalyzer_VHgg::CreateOutputFile() {

  if( useGenJets_ ) {
    std::string newflags = GetFlags() + "_GENJETS";
    SetFlags( newflags );
  }

  TreeAnalyzer::CreateOutputFile();

  

  jetTree_->Branch("run",&run_,"run_/I");
  jetTree_->Branch("event",&event_,"event_/I");
  jetTree_->Branch("LS",&LS_,"LS_/I");
  jetTree_->Branch("nvertex",&nvertex_,"nvertex_/I");

  jetTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  jetTree_->Branch("nPU",&nPU_,"nPU_/I");

  jetTree_->Branch("rhoCalo",&rhoCalo_,"rhoCalo_/F");
  jetTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");

  jetTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");


  jetTree_->Branch("passed_Photon10", &passed_Photon10_, "passed_Photon10_/O");
  jetTree_->Branch("passed_Photon15", &passed_Photon15_, "passed_Photon15_/O");
  jetTree_->Branch("passed_Photon20", &passed_Photon20_, "passed_Photon20_/O");
  jetTree_->Branch("passed_Photon25", &passed_Photon25_, "passed_Photon25_/O");
  jetTree_->Branch("passed_Photon30", &passed_Photon30_, "passed_Photon30_/O");
  jetTree_->Branch("passed_Photon35", &passed_Photon35_, "passed_Photon35_/O");
  jetTree_->Branch("passed_Photon40", &passed_Photon40_, "passed_Photon40_/O");
  jetTree_->Branch("passed_Photon50", &passed_Photon50_, "passed_Photon50_/O");
  jetTree_->Branch("passed_Photon60", &passed_Photon60_, "passed_Photon60_/O");
  jetTree_->Branch("passed_Photon70", &passed_Photon70_, "passed_Photon70_/O");
  jetTree_->Branch("passed_Photon75", &passed_Photon75_, "passed_Photon75_/O");

  jetTree_->Branch("ePhot1",  &ePhot1_,  "ePhot1_/F");
  jetTree_->Branch("ptPhot1",  &ptPhot1_,  "ptPhot1_/F");
  jetTree_->Branch("etaPhot1",  &etaPhot1_,  "etaPhot1_/F");
  jetTree_->Branch("phiPhot1",  &phiPhot1_,  "phiPhot1_/F");
  jetTree_->Branch("hcalIsoPhot1",  &hcalIsoPhot1_,  "hcalIsoPhot1_/F");
  jetTree_->Branch("ecalIsoPhot1",  &ecalIsoPhot1_,  "ecalIsoPhot1_/F");
  jetTree_->Branch("nTrkIsoPhot1",  &nTrkIsoPhot1_,  "nTrkIsoPhot1_/I");
  jetTree_->Branch("ptTrkIsoPhot1",  &ptTrkIsoPhot1_,  "ptTrkIsoPhot1_/F");
  jetTree_->Branch("clusterMajPhot1",  &clusterMajPhot1_,  "clusterMajPhot1_/F");
  jetTree_->Branch("clusterMinPhot1",  &clusterMinPhot1_,  "clusterMinPhot1_/F");
  jetTree_->Branch("hasPixelSeedPhot1",  &hasPixelSeedPhot1_,  "hasPixelSeedPhot1_/I");
  jetTree_->Branch("pid_twrHCALPhot1",  &pid_twrHCALPhot1_,  "pid_twrHCALPhot1_/F");
  jetTree_->Branch("pid_HoverEPhot1",  &pid_HoverEPhot1_,  "pid_HoverEPhot1_/F");
  jetTree_->Branch("pid_jurECALPhot1",  &pid_jurECALPhot1_,  "pid_jurECALPhot1_/F");
  jetTree_->Branch("pid_sIEtaIEtaPhot1",  &pid_sIEtaIEtaPhot1_,  "pid_sIEtaIEtaPhot1_/F");
  jetTree_->Branch("pid_hlwTrackPhot1",  &pid_hlwTrackPhot1_,  "pid_hlwTrackPhot1_/F");
  jetTree_->Branch("pid_hlwTrackNoDzPhot1",  &pid_hlwTrackNoDzPhot1_,  "pid_hlwTrackNoDzPhot1_/F");

  jetTree_->Branch("ePhotGen1",  &ePhotGen1_,  "ePhotGen1_/F");
  jetTree_->Branch("ptPhotGen1",  &ptPhotGen1_,  "ptPhotGen1_/F");
  jetTree_->Branch("etaPhotGen1",  &etaPhotGen1_,  "etaPhotGen1_/F");
  jetTree_->Branch("phiPhotGen1",  &phiPhotGen1_,  "phiPhotGen1_/F");

  jetTree_->Branch("ePhot2",  &ePhot2_,  "ePhot2_/F");
  jetTree_->Branch("ptPhot2",  &ptPhot2_,  "ptPhot2_/F");
  jetTree_->Branch("etaPhot2",  &etaPhot2_,  "etaPhot2_/F");
  jetTree_->Branch("phiPhot2",  &phiPhot2_,  "phiPhot2_/F");
  jetTree_->Branch("hcalIsoPhot2",  &hcalIsoPhot2_,  "hcalIsoPhot2_/F");
  jetTree_->Branch("ecalIsoPhot2",  &ecalIsoPhot2_,  "ecalIsoPhot2_/F");
  jetTree_->Branch("nTrkIsoPhot2",  &nTrkIsoPhot2_,  "nTrkIsoPhot2_/I");
  jetTree_->Branch("ptTrkIsoPhot2",  &ptTrkIsoPhot2_,  "ptTrkIsoPhot2_/F");
  jetTree_->Branch("clusterMajPhot2",  &clusterMajPhot2_,  "clusterMajPhot2_/F");
  jetTree_->Branch("clusterMinPhot2",  &clusterMinPhot2_,  "clusterMinPhot2_/F");
  jetTree_->Branch("hasPixelSeedPhot2",  &hasPixelSeedPhot2_,  "hasPixelSeedPhot2_/I");
  jetTree_->Branch("pid_twrHCALPhot2",  &pid_twrHCALPhot2_,  "pid_twrHCALPhot2_/F");
  jetTree_->Branch("pid_HoverEPhot2",  &pid_HoverEPhot2_,  "pid_HoverEPhot2_/F");
  jetTree_->Branch("pid_jurECALPhot2",  &pid_jurECALPhot2_,  "pid_jurECALPhot2_/F");
  jetTree_->Branch("pid_sIEtaIEtaPhot2",  &pid_sIEtaIEtaPhot2_,  "pid_sIEtaIEtaPhot2_/F");
  jetTree_->Branch("pid_hlwTrackPhot2",  &pid_hlwTrackPhot2_,  "pid_hlwTrackPhot2_/F");
  jetTree_->Branch("pid_hlwTrackNoDzPhot2",  &pid_hlwTrackNoDzPhot2_,  "pid_hlwTrackNoDzPhot2_/F");

  jetTree_->Branch("ePhotGen2",  &ePhotGen2_,  "ePhotGen2_/F");
  jetTree_->Branch("ptPhotGen2",  &ptPhotGen2_,  "ptPhotGen2_/F");
  jetTree_->Branch("etaPhotGen2",  &etaPhotGen2_,  "etaPhotGen2_/F");
  jetTree_->Branch("phiPhotGen2",  &phiPhotGen2_,  "phiPhotGen2_/F");

  

  jetTree_->Branch("nJet",  &nJet_,  "nJet_/I");
  jetTree_->Branch("eJet",  eJet_,  "eJet_[nJet_]/F");
  jetTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJet_]/F");
  jetTree_->Branch( "ptRawJet",  ptRawJet_,  "ptRawJet_[nJet_]/F");
  jetTree_->Branch("etaJet", etaJet_, "etaJet_[nJet_]/F");
  jetTree_->Branch("phiJet", phiJet_, "phiJet_[nJet_]/F");
  jetTree_->Branch( "ptDJet",  ptDJet_,  "ptDJet_[nJet_]/F");
  jetTree_->Branch( "rmsCandJet",  rmsCandJet_,  "rmsCandJet_[nJet_]/F");
  jetTree_->Branch("trackCountingHighEffBJetTagsJet",  trackCountingHighEffBJetTagsJet_,  "trackCountingHighEffBJetTagsJet_[nJet_]/F");
  jetTree_->Branch("simpleSecondaryVertexHighEffBJetTagsJet",  simpleSecondaryVertexHighEffBJetTagsJet_,  "simpleSecondaryVertexHighEffBJetTagsJet_[nJet_]/F");
  jetTree_->Branch(  "eJetGen",   eJetGen_,   "eJetGen_[nJet_]/F");
  jetTree_->Branch(  "ptJetGen",   ptJetGen_,   "ptJetGen_[nJet_]/F");
  jetTree_->Branch( "etaJetGen",  etaJetGen_,  "etaJetGen_[nJet_]/F");
  jetTree_->Branch( "phiJetGen",  phiJetGen_,  "phiJetGen_[nJet_]/F");
  jetTree_->Branch("pdgIdPartJet", pdgIdPartJet_, "pdgIdPartJet_[nJet_]/I");
  jetTree_->Branch("pdgIdMomJet", pdgIdMomJet_, "pdgIdMomJet_[nJet_]/I");
  jetTree_->Branch(   "ptPartJet",    ptPartJet_,    "ptPartJet_[nJet_]/F");
  jetTree_->Branch(  "etaPartJet",   etaPartJet_,   "etaPartJet_[nJet_]/F");
  jetTree_->Branch(  "phiPartJet",   phiPartJet_,   "phiPartJet_[nJet_]/F");

  jetTree_->Branch("eChargedHadronsJet", &eChargedHadronsJet_, "eChargedHadronsJet_[nJet_]/F");
  jetTree_->Branch("ePhotonsJet", &ePhotonsJet_, "ePhotonsJet_[nJet_]/F");
  jetTree_->Branch("eNeutralHadronsJet", &eNeutralHadronsJet_, "eNeutralHadronsJet_[nJet_]/F");
  jetTree_->Branch("eMuonsJet", &eMuonsJet_, "eMuonsJet_[nJet_]/F");
  jetTree_->Branch("eElectronsJet", &eElectronsJet_, "eElectronsJet_[nJet_]/F");
  jetTree_->Branch("eHFHadronsJet", &eHFHadronsJet_, "eHFHadronsJet_[nJet_]/F");
  jetTree_->Branch("eHFEMJet", &eHFEMJet_, "eHFEMJet_[nJet_]/F");

  jetTree_->Branch("nChargedHadronsJet", &nChargedHadronsJet_, "nChargedHadronsJet_[nJet_]/I");
  jetTree_->Branch("nPhotonsJet", &nPhotonsJet_, "nPhotonsJet_[nJet_]/I");
  jetTree_->Branch("nNeutralHadronsJet", &nNeutralHadronsJet_, "nNeutralHadronsJet_[nJet_]/I");
  jetTree_->Branch("nMuonsJet", &nMuonsJet_, "nMuonsJet_[nJet_]/I");
  jetTree_->Branch("nElectronsJet", &nElectronsJet_, "nElectronsJet_[nJet_]/I");
  jetTree_->Branch("nHFHadronsJet", &nHFHadronsJet_, "nHFHadronsJet_[nJet_]/I");
  jetTree_->Branch("nHFEMJet", &nHFEMJet_, "nHFEMJet_[nJet_]/I");

  jetTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  jetTree_->Branch("epfMetCorr",&epfMetCorr_,"epfMetCorr_/F");
  jetTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");
  jetTree_->Branch("eMet",&eMet_,"eMet_/F");
  jetTree_->Branch("phiMet",&phiMet_,"phiMet_/F");
  jetTree_->Branch("etcMet",&etcMet_,"etcMet_/F");
  jetTree_->Branch("phitcMet",&phitcMet_,"phitcMet_/F");

} 



TreeAnalyzer_VHgg::~TreeAnalyzer_VHgg() {

  outfile_->cd();

}



void TreeAnalyzer_VHgg::Loop()
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

     if( (jentry%20000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;

     run_ = run;
     LS_ = lbn;
     event_ = event;
     nvertex_ = nvertex;
     ptHat_ = genpt;

     nPU_ = nPU;

     rhoCalo_ = rhoCalo;
     rhoPF_ = rhoPF;

     if( !isGoodLS() ) continue; //this takes care also of integrated luminosity

     // good primary vertex requirement:
     if( nvertex==0 ) continue;
     bool goodVertex = (vndof[0] >= 4.0 && sqrt(vx[0]*vx[0]+vy[0]*vy[0]) < 2. && fabs(vz[0]) < 24. );
     if( !goodVertex ) continue;
 
     if( isMC )
       if( (ptHat_ > ptHatMax_) || (ptHat_ < ptHatMin_) ) continue;


     epfMet_ = epfMet;
     phipfMet_ = phipfMet;
     eMet_ = eMet;
     phiMet_ = phiMet;
     etcMet_ = etcMet;
     phitcMet_ = phitcMet;



     // this is a feasibility test: first look for MC higgs photons
     int iFoundPhotGen=0;
     for( unsigned int iMC=0; iMC<nMC; ++iMC ) {

       if( statusMC[iMC]!=3 ) continue;
       if( pdgIdMC[iMC]!=22 ) continue;
       //if( pdgIdMC[motherIDMC[iMC]]!=25 ) continue;

       if( iFoundPhotGen==0 ) {
         ptPhotGen1_ = ptMC[iMC];
         etaPhotGen1_ = etaMC[iMC];
         phiPhotGen1_ = phiMC[iMC];
         ePhotGen1_ = eMC[iMC];
         iFoundPhotGen++;
       } else if( iFoundPhotGen==1 ) {
         ptPhotGen2_ = ptMC[iMC];
         etaPhotGen2_ = etaMC[iMC];
         phiPhotGen2_ = phiMC[iMC];
         ePhotGen2_ = eMC[iMC];
         iFoundPhotGen++;
       } else {
         std::cout << "More than two gen photons found!!!" << std::endl;
       }

     } //for MC

     if( iFoundPhotGen!=2 ) {
       std::cout << "There is a problem. Didn't find the two Higgs photons." << std::endl;
     }

     TLorentzVector photGen1, photGen2;
     photGen1.SetPtEtaPhiE( ptPhotGen1_, etaPhotGen1_, phiPhotGen1_, ePhotGen1_ );
     photGen2.SetPtEtaPhiE( ptPhotGen2_, etaPhotGen2_, phiPhotGen2_, ePhotGen2_ );

 
     // match to reco:
     float bestDeltaR1=999.;
     float bestDeltaR2=999.;
     AnalysisPhoton phot1, phot2;

     if( nPhot<2 ) continue;

     for( unsigned iPhot=0; iPhot<nPhot; ++iPhot ) {

       AnalysisPhoton thisPhot;

       thisPhot.e = ePhot[iPhot];
       thisPhot.pt = ptPhot[iPhot];
       thisPhot.eta = etaPhot[iPhot];
       thisPhot.phi = phiPhot[iPhot];

       thisPhot.SetPtEtaPhiE( thisPhot.pt, thisPhot.eta, thisPhot.phi, thisPhot.e );

       thisPhot.hcalIso = hcalovecal04Phot[iPhot]; 
       thisPhot.ecalIso = ecaliso04Phot[iPhot]/thisPhot.e;
       thisPhot.nTrkIso = ntrkiso035Phot[iPhot];
       thisPhot.ptTrkIso = ptiso035Phot[iPhot]/thisPhot.pt;
       thisPhot.clusterMaj = sMajMajPhot[iPhot];
       thisPhot.clusterMin = sMinMinPhot[iPhot];
       thisPhot.hasPixelSeed = hasPixelSeedPhot[iPhot];
       thisPhot.pid_twrHCAL = pid_twrHCAL[iPhot];
       thisPhot.pid_jurECAL = pid_jurECAL[iPhot];
       thisPhot.pid_HoverE = pid_HoverE[iPhot];
       thisPhot.pid_hlwTrack = pid_hlwTrack[iPhot];
       thisPhot.pid_hlwTrackNoDz = pid_hlwTrackNoDz[iPhot];
       thisPhot.pid_etawid = pid_etawid[iPhot];

       float thisDeltaR1 = thisPhot.DeltaR(photGen1);
       float thisDeltaR2 = thisPhot.DeltaR(photGen2);

       if( thisDeltaR1 < bestDeltaR1 ) {
         phot1 = thisPhot;
         bestDeltaR1 = thisDeltaR1;
       }

       if( thisDeltaR2 < bestDeltaR2 ) {
         phot2 = thisPhot;
         bestDeltaR2 = thisDeltaR2;
       }

     } //for reco photons


     if( bestDeltaR1>0.1 || bestDeltaR2>0.1 || phot1.Pt()==phot2.Pt() ) continue;


     ePhot1_ = phot1.e;
     ptPhot1_ = phot1.pt;
     etaPhot1_ = phot1.eta;
     phiPhot1_ = phot1.phi;
     hcalIsoPhot1_ = phot1.hcalIso;
     ecalIsoPhot1_ = phot1.ecalIso;
     nTrkIsoPhot1_ = phot1.nTrkIso;
     ptTrkIsoPhot1_ = phot1.ptTrkIso;
     clusterMajPhot1_ = phot1.clusterMaj;
     clusterMinPhot1_ = phot1.clusterMin;
     hasPixelSeedPhot1_ = phot1.hasPixelSeed;
     pid_twrHCALPhot1_ = phot1.pid_twrHCAL;
     pid_HoverEPhot1_ = phot1.pid_HoverE;
     pid_jurECALPhot1_ = phot1.pid_jurECAL;
     pid_sIEtaIEtaPhot1_ = phot1.pid_etawid;
     pid_hlwTrackPhot1_ = phot1.pid_hlwTrack;
     pid_hlwTrackNoDzPhot1_ = phot1.pid_hlwTrackNoDz;

     ePhot2_ = phot2.e;
     ptPhot2_ = phot2.pt;
     etaPhot2_ = phot2.eta;
     phiPhot2_ = phot2.phi;
     hcalIsoPhot2_ = phot2.hcalIso;
     ecalIsoPhot2_ = phot2.ecalIso;
     nTrkIsoPhot2_ = phot2.nTrkIso;
     ptTrkIsoPhot2_ = phot2.ptTrkIso;
     clusterMajPhot2_ = phot2.clusterMaj;
     clusterMinPhot2_ = phot2.clusterMin;
     hasPixelSeedPhot2_ = phot2.hasPixelSeed;
     pid_twrHCALPhot2_ = phot2.pid_twrHCAL;
     pid_HoverEPhot2_ = phot2.pid_HoverE;
     pid_jurECALPhot2_ = phot2.pid_jurECAL;
     pid_sIEtaIEtaPhot2_ = phot2.pid_etawid;
     pid_hlwTrackPhot2_ = phot2.pid_hlwTrack;
     pid_hlwTrackNoDzPhot2_ = phot2.pid_hlwTrackNoDz;



     // and now jets:
     std::vector<AnalysisJet*> jets;
     nJet_ = 0;
     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       if( nJet_>=20 ) break;

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];

       thisJet.SetPtEtaPhiE( thisJet.ptReco, thisJet.etaReco, thisJet.phiReco, thisJet.eReco );

       if( thisJet.ptReco<20. ) continue;

       thisJet.eTracksReco = eChargedHadrons[iRecoJet];
       thisJet.ePhotonsReco = ePhotons[iRecoJet];
       thisJet.eNeutralHadronsReco = eNeutralHadrons[iRecoJet];
       thisJet.eMuonsReco = eMuons[iRecoJet];
       thisJet.eElectronsReco = eElectrons[iRecoJet];
       thisJet.eHFHadronsReco = eHFHadrons[iRecoJet];
       thisJet.eHFEMReco = eHFEM[iRecoJet];

       thisJet.ptD = ptDJet[iRecoJet];
       thisJet.rmsCand = rmsCandJet[iRecoJet];

       thisJet.trackCountingHighEffBJetTags = trackCountingHighEffBJetTags[iRecoJet];
       thisJet.simpleSecondaryVertexHighEffBJetTags = simpleSecondaryVertexHighEffBJetTags[iRecoJet];

       thisJet.nTracksReco = nChargedHadrons[iRecoJet];
       thisJet.nPhotonsReco = nPhotons[iRecoJet];
       thisJet.nNeutralHadronsReco = nNeutralHadrons[iRecoJet];
       thisJet.nMuonsReco = nMuons[iRecoJet];
       thisJet.nElectronsReco = nElectrons[iRecoJet];
       thisJet.nHFHadronsReco = nHFHadrons[iRecoJet];
       thisJet.nHFEMReco = nHFEM[iRecoJet];


       if( thisJet.DeltaR(phot1)>0.5 && thisJet.DeltaR(phot2)>0.5 ) {

         //match to parton:
         int i_foundPart=-1;
         float bestDeltaRPart=999.;
         for( unsigned iMC=0; iMC<nMC; ++iMC ) {
 
           if( statusMC[iMC]!=3 ) continue;
           if( !(fabs(pdgIdMC[iMC])<7 || pdgIdMC[iMC]==21) ) continue;

           TLorentzVector thisPart;
           thisPart.SetPtEtaPhiE( ptMC[iMC], etaMC[iMC], phiMC[iMC], eMC[iMC] );
               

           float thisDeltaR = thisPart.DeltaR(thisJet);
      
           if( thisDeltaR<bestDeltaRPart ) {
             bestDeltaRPart = thisDeltaR;
             i_foundPart = iMC;
           }

         } //for partons

         if( i_foundPart!=-1 ) {
           thisJet.ptPart = ptMC[i_foundPart];
           thisJet.etaPart = etaMC[i_foundPart];
           thisJet.phiPart = phiMC[i_foundPart];
           thisJet.ePart = eMC[i_foundPart];
           thisJet.pdgIdPart = pdgIdMC[i_foundPart];
           thisJet.pdgIdMom = pdgIdMC[motherIDMC[i_foundPart]];
         }

         AnalysisJet* newJet = new AnalysisJet(thisJet);
         jets.push_back(newJet);
         nJet_++;

       } //if not matched to photons
 

     } //for reco jets

     
     if( jets.size()<2 ) continue;


     for( unsigned iJet=0; iJet<jets.size(); iJet++ ) {

       ptJet_[iJet] = jets[iJet]->ptReco;

       eJet_[iJet]  =  jets[iJet]->eReco;
       ptJet_[iJet]  =  jets[iJet]->ptCorrReco;
       ptRawJet_[iJet]  =  jets[iJet]->ptReco;
       phiJet_[iJet] = jets[iJet]->phiReco;
       etaJet_[iJet] = jets[iJet]->etaReco;

       eChargedHadronsJet_[iJet] = jets[iJet]->eTracksReco;
       ePhotonsJet_[iJet] = jets[iJet]->ePhotonsReco;
       eNeutralHadronsJet_[iJet] = jets[iJet]->eNeutralHadronsReco;
       eMuonsJet_[iJet] = jets[iJet]->eMuonsReco;
       eElectronsJet_[iJet] = jets[iJet]->eElectronsReco;
       eHFHadronsJet_[iJet] = jets[iJet]->eHFHadronsReco;
       eHFEMJet_[iJet] = jets[iJet]->eHFEMReco;

       ptDJet_[iJet]= jets[iJet]->ptD;
       rmsCandJet_[iJet]= jets[iJet]->rmsCand;

       trackCountingHighEffBJetTagsJet_[iJet]= jets[iJet]->trackCountingHighEffBJetTags;
       simpleSecondaryVertexHighEffBJetTagsJet_[iJet]= jets[iJet]->simpleSecondaryVertexHighEffBJetTags;

       nChargedHadronsJet_[iJet] = jets[iJet]->nTracksReco;
       nPhotonsJet_[iJet] = jets[iJet]->nPhotonsReco;
       nNeutralHadronsJet_[iJet] = jets[iJet]->nNeutralHadronsReco;
       nMuonsJet_[iJet] = jets[iJet]->nMuonsReco;
       nElectronsJet_[iJet] = jets[iJet]->nElectronsReco;
       nHFHadronsJet_[iJet] = jets[iJet]->nHFHadronsReco;
       nHFEMJet_[iJet] = jets[iJet]->nHFEMReco;

       ePartJet_[iJet]  =  jets[iJet]->ePart;
       ptPartJet_[iJet]  =  jets[iJet]->ptPart;
       phiPartJet_[iJet] = jets[iJet]->phiPart;
       etaPartJet_[iJet] = jets[iJet]->etaPart;
       pdgIdPartJet_[iJet] = jets[iJet]->pdgIdPart;
       pdgIdMomJet_[iJet] = jets[iJet]->pdgIdMom;

     }



     passed_Photon10_ = passedTrigger_regexp("HLT_Photon10_");
     passed_Photon15_ = passedTrigger_regexp("HLT_Photon15_");
     passed_Photon20_ = passedTrigger_regexp("HLT_Photon20_");
     passed_Photon25_ = passedTrigger_regexp("HLT_Photon25_");
     passed_Photon30_ = passedTrigger_regexp("HLT_Photon30_");
     passed_Photon35_ = passedTrigger_regexp("HLT_Photon35_");
     passed_Photon40_ = passedTrigger_regexp("HLT_Photon40_");
     passed_Photon50_ = passedTrigger_regexp("HLT_Photon50_");
     passed_Photon60_ = passedTrigger_regexp("HLT_Photon60_");
     passed_Photon70_ = passedTrigger_regexp("HLT_Photon70_");
     passed_Photon75_ = passedTrigger_regexp("HLT_Photon75_");


     jetTree_->Fill(); 

   } //for entries


 }

