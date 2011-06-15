#include "TreeAnalyzer_PhotonJetInclusive.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisPhoton.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "fitTools.h"


//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"



TreeAnalyzer_PhotonJetInclusive::TreeAnalyzer_PhotonJetInclusive( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags, bool useGenJets, TTree* tree ) :
     TreeAnalyzer( "PhotonJetInclusive", dataset, recoType, jetAlgo, flags, tree ) {

  useGenJets_=useGenJets;

} //constructor



void TreeAnalyzer_PhotonJetInclusive::CreateOutputFile() {

  if( useGenJets_ ) {
    std::string newflags = GetFlags() + "_GENJETS";
    SetFlags( newflags );
  }

  TreeAnalyzer::CreateOutputFile();

//  int nBins_pt = 16;
// Double_t Lower[nBins_pt];
//  fitTools::getPtBins( nBins_pt, Lower ); 
  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Double_t ptPhotBinning_array[200]; //ugly! no more than 200 pt bins supported
  for( unsigned i=0; i<ptPhot_binning.size(); ++i )
    ptPhotBinning_array[i] = ptPhot_binning[i];


  h1_ptPhot = new TH1F("ptPhot", "", 100., 0., 150.);

  h1_eff_denom_vs_pt = new TH1F("eff_denom_vs_pt", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  h1_eff_denom_vs_pt->Sumw2();
  h1_eff_num_medium_vs_pt = new TH1F("eff_num_medium_vs_pt", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  h1_eff_num_medium_vs_pt->Sumw2();
  h1_eff_num_loose_vs_pt = new TH1F("eff_num_loose_vs_pt", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  h1_eff_num_loose_vs_pt->Sumw2();
  

  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)

  jetTree_->Branch("run",&run_,"run_/I");
  jetTree_->Branch("LS",&LS_,"LS_/I");
  jetTree_->Branch("event",&event_,"event_/I");
  jetTree_->Branch("nvertex",&nvertex_,"nvertex_/I");

  jetTree_->Branch("nPU",&nPU_,"nPU_/I");
  jetTree_->Branch("ptHat",&ptHat_,"ptHat_/F");
  jetTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");

  jetTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  jetTree_->Branch("eventWeight_medium",&eventWeight_medium_,"eventWeight_medium_/F");
  jetTree_->Branch("eventWeight_loose",&eventWeight_loose_,"eventWeight_loose_/F");
  jetTree_->Branch("isIsolated_veryloose",&isIsolated_veryloose_, "isIsolated_veryloose_/O");

  jetTree_->Branch("isIsolated_hcal_loose",&isIsolated_hcal_loose_, "isIsolated_hcal_loose_/O");
  jetTree_->Branch("isIsolated_hcal_medium",&isIsolated_hcal_medium_, "isIsolated_hcal_medium_/O");
  jetTree_->Branch("isIsolated_hcal_tight",&isIsolated_hcal_tight_, "isIsolated_hcal_tight_/O");

  jetTree_->Branch("isIsolated_ecal_loose", &isIsolated_ecal_loose_,  "isIsolated_ecal_loose_/O");
  jetTree_->Branch("isIsolated_ecal_medium",&isIsolated_ecal_medium_, "isIsolated_ecal_medium_/O");
  jetTree_->Branch("isIsolated_ecal_tight", &isIsolated_ecal_tight_,  "isIsolated_ecal_tight_/O");

  jetTree_->Branch("isIsolated_ptTracks_loose", &isIsolated_ptTracks_loose_,  "isIsolated_ptTracks_loose_/O");
  jetTree_->Branch("isIsolated_ptTracks_medium",&isIsolated_ptTracks_medium_, "isIsolated_ptTracks_medium_/O");
  jetTree_->Branch("isIsolated_ptTracks_tight", &isIsolated_ptTracks_tight_,  "isIsolated_ptTracks_tight_/O");

  jetTree_->Branch("isIsolated_nTracks_loose", &isIsolated_nTracks_loose_,  "isIsolated_nTracks_loose_/O");
  jetTree_->Branch("isIsolated_nTracks_medium",&isIsolated_nTracks_medium_, "isIsolated_nTracks_medium_/O");
  jetTree_->Branch("isIsolated_nTracks_tight", &isIsolated_nTracks_tight_,  "isIsolated_nTracks_tight_/O");

  jetTree_->Branch("clusterMajOK_loose", &clusterMajOK_loose_,  "clusterMajOK_loose_/O");
  jetTree_->Branch("clusterMajOK_medium",&clusterMajOK_medium_, "clusterMajOK_medium_/O");
  jetTree_->Branch("clusterMajOK_tight", &clusterMajOK_tight_,  "clusterMajOK_tight_/O");

  jetTree_->Branch("clusterMinOK_loose", &clusterMinOK_loose_,  "clusterMinOK_loose_/O");
  jetTree_->Branch("clusterMinOK_medium",&clusterMinOK_medium_, "clusterMinOK_medium_/O");
  jetTree_->Branch("clusterMinOK_tight", &clusterMinOK_tight_,  "clusterMinOK_tight_/O");

  jetTree_->Branch("passedPhotonID_loose",&passedPhotonID_loose_, "passedPhotonID_loose_/O");
  jetTree_->Branch("passedPhotonID_medium",&passedPhotonID_medium_,"passedPhotonID_medium_/O");
  jetTree_->Branch("passedPhotonID_tight",&passedPhotonID_tight_, "passedPhotonID_tight_/O");
  jetTree_->Branch("matchedToMC",&matchedToMC_, "matchedToMC_/O");

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
  jetTree_->Branch("hasPixelSeedPhotReco",  &hasPixelSeedPhotReco_,  "hasPixelSeedPhotReco_/I");
  jetTree_->Branch("pid_twrHCALPhotReco",  &pid_twrHCALPhotReco_,  "pid_twrHCALPhotReco_/F");
  jetTree_->Branch("pid_HoverEPhotReco",  &pid_HoverEPhotReco_,  "pid_HoverEPhotReco_/F");
  jetTree_->Branch("pid_jurECALPhotReco",  &pid_jurECALPhotReco_,  "pid_jurECALPhotReco_/F");
  jetTree_->Branch("pid_sIEtaIEtaPhotReco",  &pid_sIEtaIEtaPhotReco_,  "pid_sIEtaIEtaPhotReco_/F");
  jetTree_->Branch("pid_hlwTrackPhotReco",  &pid_hlwTrackPhotReco_,  "pid_hlwTrackPhotReco_/F");
  jetTree_->Branch("pid_hlwTrackNoDzPhotReco",  &pid_hlwTrackNoDzPhotReco_,  "pid_hlwTrackNoDzPhotReco_/F");

  jetTree_->Branch("ePhotGen",  &ePhotGen_,  "ePhotGen_/F");
  jetTree_->Branch("ptPhotGen",  &ptPhotGen_,  "ptPhotGen_/F");
  jetTree_->Branch("etaPhotGen",  &etaPhotGen_,  "etaPhotGen_/F");
  jetTree_->Branch("phiPhotGen",  &phiPhotGen_,  "phiPhotGen_/F");

  jetTree_->Branch("deltaR_phot",  &deltaR_phot_,  "deltaR_phot_/F");

  jetTree_->Branch("nJet",  &nJet_,  "nJet_/I");
  jetTree_->Branch("eJet",  eJet_,  "eJet_[nJet_]/F");
  jetTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJet_]/F");
  jetTree_->Branch( "ptCorrJet",  ptCorrJet_,  "ptCorrJet_[nJet_]/F");
  jetTree_->Branch("etaJet", etaJet_, "etaJet_[nJet_]/F");
  jetTree_->Branch("phiJet", phiJet_, "phiJet_[nJet_]/F");
  jetTree_->Branch( "ptDJet",  ptDJet_,  "ptDJet_[nJet_]/F");
  jetTree_->Branch( "ptDJet",  ptDJet_,  "ptDJet_[nJet_]/F");
  jetTree_->Branch( "rmsCandJet",  rmsCandJet_,  "rmsCandJet_[nJet_]/F");
  jetTree_->Branch("trackCountingHighEffBJetTagsJet",  trackCountingHighEffBJetTagsJet_,  "trackCountingHighEffBJetTagsJet_[nJet_]/F");
  jetTree_->Branch(  "eJetGen",   eJetGen_,   "eJetGen_[nJet_]/F");
  jetTree_->Branch(  "ptJetGen",   ptJetGen_,   "ptJetGen_[nJet_]/F");
  jetTree_->Branch( "etaJetGen",  etaJetGen_,  "etaJetGen_[nJet_]/F");
  jetTree_->Branch( "phiJetGen",  phiJetGen_,  "phiJetGen_[nJet_]/F");
  jetTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nJet_]/I");
  jetTree_->Branch(   "ptPart",    ptPart_,    "ptPart_[nJet_]/F");
  jetTree_->Branch(  "etaPart",   etaPart_,   "etaPart_[nJet_]/F");
  jetTree_->Branch(  "phiPart",   phiPart_,   "phiPart_[nJet_]/F");

  jetTree_->Branch("eChargedHadronsJet", eChargedHadronsJet_, "eChargedHadronsJet_[nJet_]/F");
  jetTree_->Branch("ePhotonsJet", ePhotonsJet_, "ePhotonsJet_[nJet_]/F");
  jetTree_->Branch("eNeutralHadronsJet", eNeutralHadronsJet_, "eNeutralHadronsJet_[nJet_]/F");
  jetTree_->Branch("eMuonsJet", eMuonsJet_, "eMuonsJet_[nJet_]/F");
  jetTree_->Branch("eElectronsJet", eElectronsJet_, "eElectronsJet_[nJet_]/F");
  jetTree_->Branch("eHFHadronsJet", eHFHadronsJet_, "eHFHadronsJet_[nJet_]/F");
  jetTree_->Branch("eHFEMJet", eHFEMJet_, "eHFEMJet_[nJet_]/F");

  jetTree_->Branch("nChargedHadronsJet", nChargedHadronsJet_, "nChargedHadronsJet_[nJet_]/I");
  jetTree_->Branch("nPhotonsJet", nPhotonsJet_, "nPhotonsJet_[nJet_]/I");
  jetTree_->Branch("nNeutralHadronsJet", nNeutralHadronsJet_, "nNeutralHadronsJet_[nJet_]/I");
  jetTree_->Branch("nMuonsJet", nMuonsJet_, "nMuonsJet_[nJet_]/I");
  jetTree_->Branch("nElectronsJet", nElectronsJet_, "nElectronsJet_[nJet_]/I");
  jetTree_->Branch("nHFHadronsJet", nHFHadronsJet_, "nHFHadronsJet_[nJet_]/I");
  jetTree_->Branch("nHFEMJet", nHFEMJet_, "nHFEMJet_[nJet_]/I");

  jetTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  jetTree_->Branch("epfMetCorr",&epfMetCorr_,"epfMetCorr_/F");
  jetTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");
  jetTree_->Branch("eMet",&eMet_,"eMet_/F");
  jetTree_->Branch("phiMet",&phiMet_,"phiMet_/F");
  jetTree_->Branch("etcMet",&etcMet_,"etcMet_/F");
  jetTree_->Branch("phitcMet",&phitcMet_,"phitcMet_/F");

} 



TreeAnalyzer_PhotonJetInclusive::~TreeAnalyzer_PhotonJetInclusive() {

  outfile_->cd();
  h1_eff_denom_vs_pt->Write();
  h1_eff_num_medium_vs_pt->Write();
  h1_eff_num_loose_vs_pt->Write();

}



void TreeAnalyzer_PhotonJetInclusive::Loop()
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
     LS_ = lbn;
     nvertex_ = nvertex;
     ptHat_ = genpt;
     nPU_ = nPU;
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


     //default values. if usegenjets, they will be updated afterwards:
     eventWeight_medium_ = 1.;
     eventWeight_loose_ = 1.;


     //foundPhot will be either a reco photon (signal) or a genjet (if usegenjets):
     AnalysisPhoton foundPhot;


     //foundRecoPhot is the highest-pt photon candidate of the event:
     AnalysisPhoton foundRecoPhot;

     for( unsigned int iPhot=0; iPhot<nPhot; ++iPhot ) {

       AnalysisPhoton thisPhot;

       thisPhot.e = ePhot[iPhot];
       thisPhot.pt = ptPhot[iPhot];
       thisPhot.eta = etaPhot[iPhot];
       thisPhot.phi = phiPhot[iPhot];
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

       //if( thisPhot.pt < 10. ) continue;

       if( thisPhot.pt > foundRecoPhot.pt )
         foundRecoPhot = thisPhot;

     } //for reco photons

     //match to MC:
     Float_t deltaRmin_Phot = 999.;
     for( unsigned int iMC=0; iMC<nMC; ++iMC ) {
   
       if( pdgIdMC[iMC]!=22 ) continue;
   
       Float_t deltaEta = foundRecoPhot.eta - etaMC[iMC];
       Float_t deltaPhi = foundRecoPhot.phi - phiMC[iMC];
       if( deltaPhi>=TMath::Pi() ) deltaPhi-=2.*TMath::Pi();
       if( deltaPhi<=-TMath::Pi() ) deltaPhi+=2.*TMath::Pi();
       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
   
       if( deltaR < deltaRmin_Phot ) {
         deltaRmin_Phot = deltaR;
         foundRecoPhot.eGen = eMC[iMC];
         foundRecoPhot.ptGen = ptMC[iMC];
         foundRecoPhot.etaGen = etaMC[iMC];
         foundRecoPhot.phiGen = phiMC[iMC];
       }
   
     } //for mc


     //when running on QCD, photon ID will cut away all but a handful of events
     //so will substitute photon with one of the two leading genjets:
     AnalysisJet firstGenJet, secondGenJet, foundGenJet;

     //look for first genjet:
     for(unsigned int iGenJet=0; iGenJet<nJet; ++iGenJet) {

       AnalysisJet thisJet;

       thisJet.eGen  =  eJetGen[iGenJet];
       thisJet.ptGen  =  ptJetGen[iGenJet];
       thisJet.phiGen = phiJetGen[iGenJet];
       thisJet.etaGen = etaJetGen[iGenJet];
     
       if( thisJet.ptGen > firstGenJet.ptGen )
         firstGenJet= thisJet;

     } //for gen jets

     //look for second genjet:
     for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

       AnalysisJet thisJet;

       thisJet.eGen  =  eJetGen[iGenJet];
       thisJet.ptGen  =  ptJetGen[iGenJet];
       thisJet.phiGen = phiJetGen[iGenJet];
       thisJet.etaGen = etaJetGen[iGenJet];
     
       if( (thisJet.ptGen < firstGenJet.ptGen)&&(thisJet.ptGen>secondGenJet.ptGen) )
         secondGenJet= thisJet;

     } //for gen jets

     Float_t coin = rand.Uniform(1.);
     if( coin<0.5 )
       foundGenJet = firstGenJet;
     else
       foundGenJet = secondGenJet;


     if( useGenJets_ ) {
       foundPhot.e = foundGenJet.eGen;
       foundPhot.pt = foundGenJet.ptGen;
       foundPhot.eta = foundGenJet.etaGen;
       foundPhot.phi = foundGenJet.phiGen;
       foundPhot.hcalIso = 0.;
       foundPhot.ecalIso = 0.;
       foundPhot.nTrkIso = 0;
       foundPhot.ptTrkIso = 0.;
       foundPhot.clusterMaj = 0.;
       foundPhot.clusterMin = 0.;
       foundPhot.hasPixelSeed = 0;
       foundPhot.pid_jurECAL = 0.;
       foundPhot.pid_HoverE = 0.;
       foundPhot.pid_twrHCAL = 0.;
       foundPhot.pid_etawid = 0.;
       foundPhot.pid_hlwTrack = 0.;
       foundPhot.pid_hlwTrackNoDz = 0.;
       foundPhot.eGen = foundGenJet.eGen;
       foundPhot.ptGen = foundGenJet.ptGen;
       foundPhot.etaGen = foundGenJet.etaGen;
       foundPhot.phiGen = foundGenJet.phiGen;
     } else {
       foundPhot = foundRecoPhot;
     }


     //if( foundPhot.pt < 10. ) continue;
     if( foundRecoPhot.pt < 1. ) continue;


     
     isIsolated_veryloose_ = true;
     if( foundPhot.hcalIso > foundPhot.e ) isIsolated_veryloose_ = false;
     if( foundPhot.ptTrkIso > foundPhot.pt ) isIsolated_veryloose_ = false;

     isIsolated_hcal_loose_ =  foundPhot.isIsolated_hcal("loose");
     isIsolated_hcal_medium_ = foundPhot.isIsolated_hcal("medium");
     isIsolated_hcal_tight_ =  foundPhot.isIsolated_hcal("tight");

     isIsolated_ecal_loose_ =  foundPhot.isIsolated_ecal("loose");
     isIsolated_ecal_medium_ = foundPhot.isIsolated_ecal("medium");
     isIsolated_ecal_tight_ =  foundPhot.isIsolated_ecal("tight");

     isIsolated_ptTracks_loose_ =  foundPhot.isIsolated_ptTracks("loose");
     isIsolated_ptTracks_medium_ = foundPhot.isIsolated_ptTracks("medium");
     isIsolated_ptTracks_tight_ =  foundPhot.isIsolated_ptTracks("tight");

     isIsolated_nTracks_loose_ =  foundPhot.isIsolated_nTracks("loose");
     isIsolated_nTracks_medium_ = foundPhot.isIsolated_nTracks("medium");
     isIsolated_nTracks_tight_ =  foundPhot.isIsolated_nTracks("tight");

     clusterMajOK_loose_ =  foundPhot.clusterMajOK("loose");
     clusterMajOK_medium_ = foundPhot.clusterMajOK("medium");
     clusterMajOK_tight_ =  foundPhot.clusterMajOK("tight");

     clusterMinOK_loose_ =  foundPhot.clusterMinOK("loose");
     clusterMinOK_medium_ = foundPhot.clusterMinOK("medium");
     clusterMinOK_tight_ =  foundPhot.clusterMinOK("tight");

     passedPhotonID_loose_ = foundPhot.passedPhotonID("loose");
     passedPhotonID_medium_ = foundPhot.passedPhotonID("medium");
     passedPhotonID_tight_ = foundPhot.passedPhotonID("tight");

if( DEBUG_VERBOSE_ && passedPhotonID_medium_==true) {
  std::cout << "HCAL iso: " << foundPhot.hcalIso << std::endl;
  std::cout << "ECAL iso: " << foundPhot.ecalIso << std::endl;
  std::cout << "ptTrk iso: " << foundPhot.ptTrkIso << std::endl;
  std::cout << "nTrk iso: " << foundPhot.nTrkIso << std::endl;
  std::cout << "sMajMaj: " << foundPhot.hcalIso << std::endl;
  std::cout << "sMinMin: " << foundPhot.hcalIso << std::endl;
}

     ePhotReco_ = foundPhot.e;
     ptPhotReco_ = foundPhot.pt;
     etaPhotReco_ = foundPhot.eta;
     phiPhotReco_ = foundPhot.phi;
     hcalIsoPhotReco_ = foundPhot.hcalIso;
     ecalIsoPhotReco_ = foundPhot.ecalIso;
     nTrkIsoPhotReco_ = foundPhot.nTrkIso;
     ptTrkIsoPhotReco_ = foundPhot.ptTrkIso;
     clusterMajPhotReco_ = foundPhot.clusterMaj;
     clusterMinPhotReco_ = foundPhot.clusterMin;
     hasPixelSeedPhotReco_ = foundPhot.hasPixelSeed;
     pid_twrHCALPhotReco_ = foundPhot.pid_twrHCAL;
     pid_HoverEPhotReco_ = foundPhot.pid_HoverE;
     pid_jurECALPhotReco_ = foundPhot.pid_jurECAL;
     pid_sIEtaIEtaPhotReco_ = foundPhot.pid_etawid;
     pid_hlwTrackPhotReco_ = foundPhot.pid_hlwTrack;
     pid_hlwTrackNoDzPhotReco_ = foundPhot.pid_hlwTrackNoDz;
     ePhotGen_ = foundPhot.eGen;
     ptPhotGen_ = foundPhot.ptGen;
     etaPhotGen_ = foundPhot.etaGen;
     phiPhotGen_ = foundPhot.phiGen;

     Float_t deltaEta = foundPhot.eta-foundPhot.etaGen;
     Float_t deltaPhi = foundPhot.phi-foundPhot.phiGen;
     Float_t pi = 3.14159;
     if( deltaPhi >  pi ) deltaPhi -= 2.*pi;
     if( deltaPhi < -pi ) deltaPhi += 2.*pi;
     deltaR_phot_ = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

     matchedToMC_ = (deltaR_phot_<0.1);


     nJet_=0;

     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       if( nJet>=10 ) {
         break;
       }

       if( ptCorrJet[iRecoJet] < 15. ) continue;


       ptJet_[nJet_] = ptJet[iRecoJet];
       eJet_[nJet_] = eJet[iRecoJet];
       ptCorrJet_[nJet_] = ptCorrJet[iRecoJet];
       phiJet_[nJet_] = phiJet[iRecoJet];
       etaJet_[nJet_] = etaJet[iRecoJet];
       ptDJet_[nJet_] = ptDJet[iRecoJet];
       rmsCandJet_[nJet_] = rmsCandJet[iRecoJet];

       trackCountingHighEffBJetTagsJet_[nJet_] = trackCountingHighEffBJetTags[iRecoJet];
   
       eChargedHadronsJet_[nJet_] = eChargedHadrons[iRecoJet];
       ePhotonsJet_[nJet_] = ePhotons[iRecoJet];
       eNeutralHadronsJet_[nJet_] = eNeutralHadrons[iRecoJet];
       eMuonsJet_[nJet_] = eMuons[iRecoJet];
       eElectronsJet_[nJet_] = eElectrons[iRecoJet];
       eHFHadronsJet_[nJet_] = eHFHadrons[iRecoJet];
       eHFEMJet_[nJet_] = eHFEM[iRecoJet];
   
       nChargedHadronsJet_[nJet_] = nChargedHadrons[iRecoJet];
       nPhotonsJet_[nJet_] = nPhotons[iRecoJet];
       nNeutralHadronsJet_[nJet_] = nNeutralHadrons[iRecoJet];
       nMuonsJet_[nJet_] = nMuons[iRecoJet];
       nElectronsJet_[nJet_] = nElectrons[iRecoJet];
       nHFHadronsJet_[nJet_] = nHFHadrons[iRecoJet];
       nHFEMJet_[nJet_] = nHFEM[iRecoJet];



       // match to closest genjet:

       TLorentzVector thisRawJet;
       thisRawJet.SetPtEtaPhiE( ptJet[iRecoJet], etaJet[iRecoJet], phiJet[iRecoJet], eJet[iRecoJet]);

       Float_t deltaRmin = 999.;

       for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

         Float_t  eJetGen_i = eJetGen[iGenJet];
         Float_t  ptJetGen_i  =  ptJetGen[iGenJet];
         Float_t  phiJetGen_i = phiJetGen[iGenJet];
         Float_t  etaJetGen_i = etaJetGen[iGenJet];

         TLorentzVector thisGenJet;
         thisGenJet.SetPtEtaPhiE( ptJetGen[iGenJet], etaJetGen[iGenJet], phiJetGen[iGenJet], eJetGen[iGenJet] );

         Float_t deltaR = thisRawJet.DeltaR(thisGenJet);

         if( deltaR < deltaRmin ) {
           deltaRmin = deltaR;
           ptJetGen_[nJet_] = ptJetGen_i;
           phiJetGen_[nJet_] = phiJetGen_i;
           etaJetGen_[nJet_] = etaJetGen_i;
           eJetGen_[nJet_] = eJetGen_i;
         }

       } // for gen jets


       nJet_++;

   
     } //for reco jets

     
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


     bool eventOK = ( matchedToMC_ || isIsolated_veryloose_);

     if( eventOK && (ptPhotReco_>15.) )
       jetTree_->Fill(); 

     h1_ptPhot->Fill( ptPhotReco_, eventWeight_ );

     Float_t ptPhotMin = h1_eff_denom_vs_pt->GetBinLowEdge(1);
     // to compute efficiencies:
     if( fabs(foundRecoPhot.eta)<1.3 )
       h1_eff_denom_vs_pt->Fill(ptPhotGen_);
     if( foundRecoPhot.pt>ptPhotMin && fabs(foundRecoPhot.eta)<1.3 && foundRecoPhot.passedPhotonID("medium") )
       h1_eff_num_medium_vs_pt->Fill(foundRecoPhot.ptGen);
     if( foundRecoPhot.pt>ptPhotMin && fabs(foundRecoPhot.eta)<1.3 && foundRecoPhot.passedPhotonID("loose") )
       h1_eff_num_loose_vs_pt->Fill(foundRecoPhot.ptGen);
     

   } //for entries


} //loop



