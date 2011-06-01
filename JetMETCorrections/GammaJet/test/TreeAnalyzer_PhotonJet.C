#include "TreeAnalyzer_PhotonJet.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisPhoton.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "fitTools.h"


//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
//#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"



TreeAnalyzer_PhotonJet::TreeAnalyzer_PhotonJet( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags, bool useGenJets, TTree* tree ) :
     TreeAnalyzer( "PhotonJet", dataset, recoType, jetAlgo, flags, tree ) {


  useGenJets_=useGenJets;

} //constructor



void TreeAnalyzer_PhotonJet::CreateOutputFile() {

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
  jetTree_->Branch("event",&event_,"event_/I");
  jetTree_->Branch("LS",&LS_,"LS_/I");
  jetTree_->Branch("nvertex",&nvertex_,"nvertex_/I");

  jetTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  jetTree_->Branch("rhoCalo",&rhoCalo_,"rhoCalo_/F");
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

  jetTree_->Branch("eJetReco",  &eJetReco_,  "eJetReco_/F");
  jetTree_->Branch( "ptJetReco",  &ptJetReco_,  "ptJetReco_/F");
  jetTree_->Branch( "ptCorrJetReco",  &ptCorrJetReco_,  "ptCorrJetReco_/F");
  jetTree_->Branch("etaJetReco", &etaJetReco_, "etaJetReco_/F");
  jetTree_->Branch("phiJetReco", &phiJetReco_, "phiJetReco_/F");
  jetTree_->Branch( "ptDJetReco",  &ptDJetReco_,  "ptDJetReco_/F");
  jetTree_->Branch( "rmsCandJetReco",  &rmsCandJetReco_,  "rmsCandJetReco_/F");
  jetTree_->Branch(  "eJetGen",   &eJetGen_,   "eJetGen_/F");
  jetTree_->Branch(  "ptJetGen",   &ptJetGen_,   "ptJetGen_/F");
  jetTree_->Branch( "etaJetGen",  &etaJetGen_,  "etaJetGen_/F");
  jetTree_->Branch( "phiJetGen",  &phiJetGen_,  "phiJetGen_/F");
  jetTree_->Branch("pdgIdPart", &pdgIdPart_, "pdgIdPart_/I");
  jetTree_->Branch(   "ptPart",    &ptPart_,    "ptPart_/F");
  jetTree_->Branch(  "etaPart",   &etaPart_,   "etaPart_/F");
  jetTree_->Branch(  "phiPart",   &phiPart_,   "phiPart_/F");
  jetTree_->Branch("pdgIdPart2nd", &pdgIdPart2nd_, "pdgIdPart2nd_/I");
  jetTree_->Branch(   "ptPart2nd",    &ptPart2nd_,    "ptPart2nd_/F");
  jetTree_->Branch(  "etaPart2nd",   &etaPart2nd_,   "etaPart2nd_/F");
  jetTree_->Branch(  "phiPart2nd",   &phiPart2nd_,   "phiPart2nd_/F");

  jetTree_->Branch("pt2ndJetReco", &pt2ndJetReco_, "pt2ndJetReco_/F");
  jetTree_->Branch("ptCorr2ndJetReco", &ptCorr2ndJetReco_, "ptCorr2ndJetReco_/F");
  jetTree_->Branch("eta2ndJetReco", &eta2ndJetReco_, "eta2ndJetReco_/F");
  jetTree_->Branch("phi2ndJetReco", &phi2ndJetReco_, "phi2ndJetReco_/F");

  jetTree_->Branch("pt2ndJetGen", &pt2ndJetGen_, "pt2ndJetGen_/F");
  jetTree_->Branch("eta2ndJetGen", &eta2ndJetGen_, "eta2ndJetGen_/F");
  jetTree_->Branch("phi2ndJetGen", &phi2ndJetGen_, "phi2ndJetGen_/F");

  jetTree_->Branch("ptSecondaryJetsReco", &ptSecondaryJetsReco_, "ptSecondaryJetsReco_/F");
  jetTree_->Branch("ptSecondaryJetsGen", &ptSecondaryJetsGen_, "ptSecondaryJetsGen_/F");

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

  jetTree_->Branch("eTracksGen", &eTracksGen_, "eTracksGen_/F");
  jetTree_->Branch("eMuonsGen", &eMuonsGen_, "eMuonsGen_/F");
  jetTree_->Branch("eElectronsGen", &eElectronsGen_, "eElectronsGen_/F");
  jetTree_->Branch("eNeutralHadronsGen", &eNeutralHadronsGen_, "eNeutralHadronsGen_/F");
  jetTree_->Branch("ePhotonsGen", &ePhotonsGen_, "ePhotonsGen_/F");
 
  jetTree_->Branch("nTracksGen", &nTracksGen_, "nTracksGen_/I");
  jetTree_->Branch("nMuonsGen", &nMuonsGen_, "nMuonsGen_/I");
  jetTree_->Branch("nElectronsGen", &nElectronsGen_, "nElectronsGen_/I");
  jetTree_->Branch("nNeutralHadronsGen", &nNeutralHadronsGen_, "nNeutralHadronsGen_/I");
  jetTree_->Branch("nPhotonsGen", &nPhotonsGen_, "nPhotonsGen_/I");

  jetTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  jetTree_->Branch("epfMetCorr",&epfMetCorr_,"epfMetCorr_/F");
  jetTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");
  jetTree_->Branch("eMet",&eMet_,"eMet_/F");
  jetTree_->Branch("phiMet",&phiMet_,"phiMet_/F");
  jetTree_->Branch("etcMet",&etcMet_,"etcMet_/F");
  jetTree_->Branch("phitcMet",&phitcMet_,"phitcMet_/F");

} 



TreeAnalyzer_PhotonJet::~TreeAnalyzer_PhotonJet() {

  outfile_->cd();
  h1_eff_denom_vs_pt->Write();
  h1_eff_num_medium_vs_pt->Write();
  h1_eff_num_loose_vs_pt->Write();

}



void TreeAnalyzer_PhotonJet::Loop()
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
   
       //if( statusMC[iMC]!=3 ) continue;
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

     AnalysisJet firstJet;
     AnalysisJet secondJet;

     TVector2 vpfmet( epfMet*cos(phipfMet), epfMet*sin(phipfMet) );
     vpfmet *= -1.; //now its just sum(pt vectors)
     TVector2 vpfmetCorr( 0., 0.);

     //look for first jet:
     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
     //if( recoType_=="jpt" ) {
     //  if( isMC ) thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)false );
     //  else       thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)true );
     //} else {
         thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
     //}
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];
     //thisJet.eCorrReco  =  eCorrJet[iRecoJet];

       thisJet.emfReco = (recoType_=="pf") ? 0. : emfJet[iRecoJet];

       thisJet.eTracksReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eChargedHadrons[iRecoJet] : 0.;
       thisJet.ePhotonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? ePhotons[iRecoJet] : 0.;
       thisJet.eNeutralHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eNeutralHadrons[iRecoJet] : 0.;
       thisJet.eMuonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eMuons[iRecoJet] : 0.;
       thisJet.eElectronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eElectrons[iRecoJet] : 0.;
       thisJet.eHFHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eHFHadrons[iRecoJet] : 0.;
       thisJet.eHFEMReco = (recoType_=="pf" && jetAlgo_=="akt5") ? eHFEM[iRecoJet] : 0.;

       thisJet.ptD = (recoType_=="pf" && jetAlgo_=="akt5") ? ptDJet[iRecoJet] : 0.;
       thisJet.rmsCand = (recoType_=="pf" && jetAlgo_=="akt5") ? rmsCandJet[iRecoJet] : 0.;

       thisJet.nTracksReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nChargedHadrons[iRecoJet] : 0;
       thisJet.nPhotonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nPhotons[iRecoJet] : 0;
       thisJet.nNeutralHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nNeutralHadrons[iRecoJet] : 0;
       thisJet.nMuonsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nMuons[iRecoJet] : 0;
       thisJet.nElectronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nElectrons[iRecoJet] : 0;
       thisJet.nHFHadronsReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nHFHadrons[iRecoJet] : 0;
       thisJet.nHFEMReco = (recoType_=="pf" && jetAlgo_=="akt5") ? nHFEM[iRecoJet] : 0;

       //jet has to be in "Mercedes area" in transverse plane wrt phot:
       Float_t deltaPhi = fitTools::delta_phi( foundPhot.phi, thisJet.phiReco);
       Float_t pi = 3.14159;
       //if( (fabs(deltaPhi) > 2.*pi/3.) && (thisJet.ptReco > firstJet.ptReco) )
       if( (fabs(deltaPhi) > pi/2.) && (thisJet.ptReco > firstJet.ptReco) ) {
         firstJet = thisJet;
       }


     //// correct pf met by hand and create corrected missing et:
     //if( thisJet.ptCorrReco>6. ) { // correct using only jets with pt corr > 6
     // TVector2 v( thisJet.ptReco*cos(thisJet.phiReco), thisJet.ptReco*sin(thisJet.phiReco) );
     // vpfmet -= v;
     // TVector2 v2( thisJet.ptCorrReco*cos(thisJet.phiReco), thisJet.ptCorrReco*sin(thisJet.phiReco) );
     // vpfmetCorr += v2;
     //}

     } //for reco jets

     //adding corrected met (with jets down to 6 gev) with the uncorrected soft leftovers:
     //epfMetCorr_ = vpfmetCorr.Mod() + vpfmet.Mod(); 
     epfMetCorr_ = 0.;

     if( firstJet.eReco == 0. ) continue;

     


     //look for second jet:
     Float_t pxSumReco = 0.;
     Float_t pySumReco = 0.;
     Float_t pzSumReco = 0.;

     for(unsigned int iRecoJet=0; iRecoJet<nJet; ++iRecoJet) {

       AnalysisJet thisJet;

       thisJet.eReco  =  eJet[iRecoJet];
       thisJet.ptReco  =  ptJet[iRecoJet];
     //if( recoType_=="jpt" ) {
     //  if( isMC ) thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)false );
     //  else       thisJet.ptCorrReco = getCorrectedPt( ptJet[iRecoJet], etaJet[iRecoJet], (bool)true );
     //} else {
         thisJet.ptCorrReco  =  ptCorrJet[iRecoJet];
     //}
       thisJet.phiReco = phiJet[iRecoJet];
       thisJet.etaReco = etaJet[iRecoJet];

       if( (thisJet.etaReco == firstJet.etaReco)&&(thisJet.phiReco==firstJet.phiReco)&&
           (thisJet.ptReco==firstJet.ptReco) ) continue;


       Float_t deltaEta = foundPhot.eta - thisJet.etaReco;
       Float_t deltaPhi = foundPhot.phi - thisJet.phiReco;
       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       float deltaR_thresh = (jetAlgo_=="akt7") ? 0.7 : 0.5;
       if( deltaR > deltaR_thresh ) { //far away from photon
    
         pxSumReco += thisJet.pxReco();
         pySumReco += thisJet.pyReco();
         pzSumReco += thisJet.pzReco();

         if( (thisJet.ptReco < firstJet.ptReco) && (thisJet.ptReco > secondJet.ptReco) ) {
           secondJet = thisJet;
         }

       }

     } //for reco jets

 
     Float_t pSumReco = sqrt( pxSumReco*pxSumReco + pySumReco*pySumReco + pzSumReco*pzSumReco );
     Float_t thetaSumReco = (pSumReco>0.) ? acos( pzSumReco/pSumReco ) : 0.;
     Float_t ptSumReco = pSumReco*sin(thetaSumReco);
     


     Float_t deltaRmin = 999.;

     for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

       Float_t  eJetGen_i = eJetGen[iGenJet];
       Float_t  ptJetGen_i  =  ptJetGen[iGenJet];
       Float_t  phiJetGen_i = phiJetGen[iGenJet];
       Float_t  etaJetGen_i = etaJetGen[iGenJet];

       Float_t  eTracksGen_i = (jetAlgo_=="akt5") ? eTracksGen[iGenJet] : 0.;
       Float_t  ePhotonsGen_i = (jetAlgo_=="akt5") ? ePhotonsGen[iGenJet] : 0.;
       Float_t  eNeutralHadronsGen_i = (jetAlgo_=="akt5") ? eNeutralHadronsGen[iGenJet] : 0.;
       Float_t  eMuonsGen_i = (jetAlgo_=="akt5") ? eMuonsGen[iGenJet] : 0.;
       Float_t  eElectronsGen_i = (jetAlgo_=="akt5") ? eElectronsGen[iGenJet] : 0.;

       Int_t  nTracksGen_i = (jetAlgo_=="akt5") ? nTracksGen[iGenJet] : 0;
       Int_t  nPhotonsGen_i = (jetAlgo_=="akt5") ? nPhotonsGen[iGenJet] : 0;
       Int_t  nNeutralHadronsGen_i = (jetAlgo_=="akt5") ? nNeutralHadronsGen[iGenJet] : 0;
       Int_t  nMuonsGen_i = (jetAlgo_=="akt5") ? nMuonsGen[iGenJet] : 0;
       Int_t  nElectronsGen_i = (jetAlgo_=="akt5") ? nElectronsGen[iGenJet] : 0;

       Float_t deltaEta = firstJet.etaReco - etaJetGen_i;
       Float_t deltaPhi = firstJet.phiReco - phiJetGen_i;

       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR < deltaRmin ) {
         deltaRmin = deltaR;
         firstJet.ptGen = ptJetGen_i;
         firstJet.eGen = eJetGen_i;
         firstJet.etaGen = etaJetGen_i;
         firstJet.phiGen = phiJetGen_i;

         firstJet.eTracksGen = eTracksGen_i;
         firstJet.ePhotonsGen = ePhotonsGen_i;
         firstJet.eNeutralHadronsGen = eNeutralHadronsGen_i;
         firstJet.eMuonsGen = eMuonsGen_i;
         firstJet.eElectronsGen = eElectronsGen_i;

         firstJet.nTracksGen = nTracksGen_i;
         firstJet.nPhotonsGen = nPhotonsGen_i;
         firstJet.nNeutralHadronsGen = nNeutralHadronsGen_i;
         firstJet.nMuonsGen = nMuonsGen_i;
         firstJet.nElectronsGen = nElectronsGen_i;
       }

     } // for gen jets


     Float_t pxSumGen = 0.;
     Float_t pySumGen = 0.;
     Float_t pzSumGen = 0.;


     //look for second gen jet (not necessarily matched to second reco jet)
     for(unsigned int iGenJet=0; iGenJet<nJetGen; ++iGenJet) {

       AnalysisJet thisJet;

       thisJet.eGen  =    eJetGen[iGenJet];
       thisJet.ptGen  =  ptJetGen[iGenJet];
       thisJet.phiGen = phiJetGen[iGenJet];
       thisJet.etaGen = etaJetGen[iGenJet];

       if( (thisJet.etaGen == firstJet.etaGen)&&(thisJet.phiGen==firstJet.phiGen)&&
           (thisJet.ptGen==firstJet.ptGen) ) continue;


       Float_t deltaEta = foundPhot.etaGen - thisJet.etaGen;
       Float_t deltaPhi = foundPhot.phiGen - thisJet.phiGen;
       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR > 0.3 ) {  //far away from photon

         pxSumGen += thisJet.pxGen();
         pySumGen += thisJet.pyGen();
         pzSumGen += thisJet.pzGen();

         if( (thisJet.ptGen < firstJet.ptGen) && (thisJet.ptGen>secondJet.ptGen) ) {
           secondJet.eGen   = thisJet.eGen;
           secondJet.ptGen  = thisJet.ptGen;
           secondJet.phiGen = thisJet.phiGen;
           secondJet.etaGen = thisJet.etaGen;
         }

       }

     } //for Gen jets



     Float_t pSumGen = sqrt( pxSumGen*pxSumGen + pySumGen*pySumGen + pzSumGen*pzSumGen );
     Float_t thetaSumGen = (pSumGen>0.) ? acos( pzSumGen/pSumGen ) : 0.;
     Float_t ptSumGen = pSumGen*sin(thetaSumGen);
   

       eJetReco_  =  firstJet.eReco;
      ptJetReco_  =  firstJet.ptReco;
  ptCorrJetReco_  =  firstJet.ptCorrReco;
     phiJetReco_  =  firstJet.phiReco;
     etaJetReco_  =  firstJet.etaReco;
     ptDJetReco_  =  firstJet.ptD;
 rmsCandJetReco_  =  firstJet.rmsCand;
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

     pt2ndJetReco_= secondJet.ptReco;
     ptCorr2ndJetReco_= secondJet.ptCorrReco;
     eta2ndJetReco_= secondJet.etaReco;
     phi2ndJetReco_= secondJet.phiReco;

      pt2ndJetGen_= secondJet.ptGen;
     eta2ndJetGen_= secondJet.etaGen;
     phi2ndJetGen_= secondJet.phiGen;

     ptSecondaryJetsReco_ = ptSumReco;
     ptSecondaryJetsGen_ = ptSumGen;

     //look for first jet parton:
     Float_t deltaRMCmin = 999.;
     Int_t pdgIdPart_found = 0;
     Float_t etaPart_found;
     Float_t phiPart_found;
     Float_t ptPart_found;


     for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

       if( statusMC[iPartMC]!=3 ) continue;

       Float_t pt = ptMC[iPartMC];
     
       Float_t eta = etaMC[iPartMC];
       Float_t phi = phiMC[iPartMC];
       Int_t   pdgId = pdgIdMC[iPartMC];
     
       Float_t deltaEtaMC = eta-firstJet.etaGen;
       Float_t deltaPhiMC = phi-firstJet.phiGen;
       if( deltaPhiMC >= TMath::Pi() ) deltaPhiMC -= 2.*TMath::Pi();
       if( deltaPhiMC <= -TMath::Pi() ) deltaPhiMC += 2.*TMath::Pi();
     
       Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );

       bool goodPdgId = false;
       if( (fabs(pdgId)<=9) || (fabs(pdgId)==21) ) goodPdgId = true;
     
       if( (deltaRMC < deltaRMCmin) && goodPdgId ) {
         deltaRMCmin = deltaRMC;
         pdgIdPart_found = pdgIdMC[iPartMC];
         etaPart_found = eta;
         phiPart_found = phi;
         ptPart_found = pt;
       }

     } //for MC particles


     //look for second jet parton:
     Float_t deltaRMCmin_2 = 999.;
     Int_t pdgIdPart2nd_found = 0;
     Float_t etaPart2nd_found = 0.;
     Float_t phiPart2nd_found = 0.;
     Float_t ptPart2nd_found = 0.;


     for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

       if( statusMC[iPartMC]!=3 ) continue;

       Float_t pt = ptMC[iPartMC];
     
       Float_t eta = etaMC[iPartMC];
       Float_t phi = phiMC[iPartMC];
       Int_t   pdgId = pdgIdMC[iPartMC];
     
       Float_t deltaEtaMC = eta-secondJet.etaGen;
       Float_t deltaPhiMC = phi-secondJet.phiGen;
       if( deltaPhiMC >= TMath::Pi() ) deltaPhiMC -= 2.*TMath::Pi();
       if( deltaPhiMC <= -TMath::Pi() ) deltaPhiMC += 2.*TMath::Pi();
     
       Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );

       bool goodPdgId = false;
       if( (fabs(pdgId)<=9) || (fabs(pdgId)==21) ) goodPdgId = true;
     
       if( (deltaRMC < deltaRMCmin_2) && goodPdgId ) {
         deltaRMCmin_2 = deltaRMC;
         pdgIdPart2nd_found = pdgIdMC[iPartMC];
         etaPart2nd_found = eta;
         phiPart2nd_found = phi;
         ptPart2nd_found = pt;
       }

     } //for MC particles



     pdgIdPart_=  pdgIdPart_found;
     ptPart_=  ptPart_found;
     phiPart_= phiPart_found;
     etaPart_= etaPart_found;

     pdgIdPart2nd_=  pdgIdPart2nd_found;
     ptPart2nd_=  ptPart2nd_found;
     phiPart2nd_= phiPart2nd_found;
     etaPart2nd_= etaPart2nd_found;


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


/*

   // now if using genjets
   // correct weights with photon ID efficiency

   if( useGenJets_ ) {

     std::cout << "-> Correcting eventweights with photon ID efficiencies." << std::endl;

     TH1F* h1_effloose = new TH1F(*h1_eff_num_loose_vs_pt);
     h1_effloose->SetName("effloose_vs_pt");
     h1_effloose->Divide( h1_eff_denom_vs_pt );

     TH1F* h1_effmedium = new TH1F(*h1_eff_num_medium_vs_pt);
     h1_effmedium->SetName("effmedium_vs_pt");
     h1_effmedium->Divide( h1_eff_denom_vs_pt );

     Float_t oldWeight=eventWeight_;
     //jetTree_->SetBranchStatus( "eventWeight", 0 );
     Float_t ptPhot_tmp;
     //jetTree_->SetBranchAddress( "ptPhotReco", &ptPhot_tmp );
     jetTree_->SetBranchAddress( "ptPhotGen", &ptPhot_tmp );

     TTree* newTree = jetTree_->CloneTree(0);


     int nentries = jetTree_->GetEntries();


     for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

       jetTree_->GetEntry(ientry);

       if( (ientry % 10000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

       int thebin = h1_eff_denom_vs_pt->FindBin( ptPhot_tmp );

       //std::cout << "ptPhot: " << ptPhot_tmp << "\tbin: " << thebin;

       if( thebin<1 || thebin>h1_eff_denom_vs_pt->GetNbinsX() ) {
         eventWeight_medium_= 0.; //scary
         eventWeight_loose_= 0.; //scary
       //std::cout << "\t-> continued." << std::endl;
         continue;
       }

       if( h1_eff_denom_vs_pt->GetBinContent(thebin)==0. ) {
         eventWeight_medium_= 0.; 
         eventWeight_loose_= 0.; 
       } else {
         //std::cout << "\th1_effloose(bin): " << h1_effloose->GetBinContent(thebin) << "\th1_effmedium(bin): " << h1_effmedium->GetBinContent( thebin ) << std::endl;
         eventWeight_loose_  = oldWeight*h1_effloose->GetBinContent( thebin );
         eventWeight_medium_ = oldWeight*h1_effmedium->GetBinContent( thebin );
       }

       newTree->Fill();

     } //for entries

     jetTree_ = newTree;


   } //if usegenjets
*/


} //loop



/*
float TreeAnalyzer_PhotonJet::getCorrectedPt( float pt, float eta, bool applyAlsoResidual ) {

  std::string acronym;
  if( recoType_=="calo" )
    acronym = "Calo";
  else if( recoType_=="pf" )
    acronym = "PF";
  else if( recoType_=="jpt" )
    acronym = "JPT";
  else {
    std::cout << "RECOTYPE: " << recoType_ << " not implemented. Exiting" << std::endl;
    exit(9191);
  }
  
  std::string L2_filename = "Spring10_L2Relative_AK5" + acronym + ".txt"; 
  std::string L3_filename = "Spring10_L3Absolute_AK5" + acronym + ".txt"; 
  std::string Residual_filename = "Spring10DataV1_L2L3Residual_AK5" + acronym + ".txt"; 

//std::cout << "eta: " << eta << " ptraw: " << pt;
  pt *= getCorrectionFactor(L2_filename, pt, eta);
//std::cout << " -> L2: " << pt;
  pt *= getCorrectionFactor(L3_filename, pt, eta);
//std::cout << " -> L3: " << pt;
  if( applyAlsoResidual )
    pt *= getCorrectionFactor(Residual_filename, pt, eta);
//std::cout << " -> Residual: " << pt << std::endl;;

  return pt;

}


float TreeAnalyzer_PhotonJet::getCorrectionFactor( const std::string& fileName, float pt, float eta ) {

  std::string filePath = "CondFormats/JetMETObjects/data/" + fileName;
  JetCorrectorParameters *JetCorPar = new JetCorrectorParameters(filePath);
  std::vector<JetCorrectorParameters> vParam;
  vParam.push_back(*JetCorPar);
  FactorizedJetCorrector* jec = new FactorizedJetCorrector(vParam);
  jec->setJetEta(eta);
  jec->setJetPt(pt);
  double corr = jec->getCorrection();
  delete JetCorPar;
  JetCorPar=0;
  delete jec;
  jec=0;

  return corr;

}
*/
