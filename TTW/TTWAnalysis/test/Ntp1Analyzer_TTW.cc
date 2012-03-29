#include "Ntp1Analyzer_TTW.h"


#include <iostream>
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRegexp.h"
#include "TMVA/Reader.h"

#include "AnalysisElectron.h"
#include "AnalysisMuon.h"
#include "AnalysisJet.h"

#include "PUWeight.h"


//#include "fitTools.h"


int DEBUG_EVENTNUMBER = 4213;


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
float getWeightPU(Int_t nPU);
//std::vector<AnalysisLepton> findLeptons( std::vector<AnalysisLepton> lPlus, std::vector<AnalysisLepton> lMinus );
bool findZ( std::vector<AnalysisLepton> leptPlus, std::vector<AnalysisLepton> leptMinus, float ZVetoWindow );






Ntp1Analyzer_TTW::Ntp1Analyzer_TTW( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "TTW", dataset, flags, tree ) {


  h1_nCounter_Zee_ = new TH1D("nCounter_Zee", "", 1, 0., 1.);
  h1_nCounter_Zmumu_ = new TH1D("nCounter_Zmumu", "", 1, 0., 1.);


} //constructor



void Ntp1Analyzer_TTW::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nPU",&nPU_,"nPU_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");
  reducedTree_->Branch("genWeight",&genWeight_,"genWeight_/F");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  reducedTree_->Branch("eventWeightPU",&eventWeightPU_,"eventWeightPU_/F");
  reducedTree_->Branch("eventWeightPU_ave",&eventWeightPU_ave_,"eventWeightPU_ave_/F");
  reducedTree_->Branch("leptTypeMC",&leptTypeMC_,"leptTypeMC_/I");

//// triggers:
//reducedTree_->Branch("HLT_Mu11", &HLT_Mu11_, "HLT_Mu11_/O");
//reducedTree_->Branch("HLT_Ele17_SW_EleId_L1R", &HLT_Ele17_SW_EleId_L1R_, "HLT_Ele17_SW_EleId_L1R_/O");
//reducedTree_->Branch("HLT_DoubleMu3", &HLT_DoubleMu3_, "HLT_DoubleMu3_/O");
  reducedTree_->Branch("passed_HLT_DoubleMu6", &passed_HLT_DoubleMu6_,"passed_HLT_DoubleMu6_/O)");
  reducedTree_->Branch("passed_HLT_DoubleMu7", &passed_HLT_DoubleMu7_,"passed_HLT_DoubleMu7_/O)");
  reducedTree_->Branch("passed_HLT_Mu13_Mu8",  &passed_HLT_Mu13_Mu8_, "passed_HLT_Mu13_Mu8_/O)");
  reducedTree_->Branch("passed_HLT_Mu17_Mu8",  &passed_HLT_Mu17_Mu8_, "passed_HLT_Mu17_Mu8_/O)");
  reducedTree_->Branch("passed_HLT_IsoMu17",   &passed_HLT_IsoMu17_,  "passed_HLT_IsoMu17_/O)");
  reducedTree_->Branch("passed_HLT_IsoMu24",   &passed_HLT_IsoMu24_,  "passed_HLT_IsoMu24_/O)");
  reducedTree_->Branch("passed_HLT_Mu8_Jet40",   &passed_HLT_Mu8_Jet40_,  "passed_HLT_Mu8_Jet40_/O)");
  reducedTree_->Branch("passed_HLT_L2DoubleMu23_NoVertex",   &passed_HLT_L2DoubleMu23_NoVertex_,  "passed_HLT_L2DoubleMu23_NoVertex_/O)");
  reducedTree_->Branch("passed_HLT_L2DoubleMu30_NoVertex",   &passed_HLT_L2DoubleMu30_NoVertex_,  "passed_HLT_L2DoubleMu30_NoVertex_/O)");
  reducedTree_->Branch("passed_HLT_TripleMu5",   &passed_HLT_TripleMu5_,  "passed_HLT_TripleMu5_/O)");
  reducedTree_->Branch("passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL", &passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_, "passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_/O");
  reducedTree_->Branch("passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_, "passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_/O");
  reducedTree_->Branch("passed_HLT_Mu17_Ele8_CaloIdL",   &passed_HLT_Mu17_Ele8_CaloIdL_,  "passed_HLT_Mu17_Ele8_CaloIdL_/O)");
  reducedTree_->Branch("passed_HLT_Mu8_Ele17_CaloIdL",   &passed_HLT_Mu8_Ele17_CaloIdL_,  "passed_HLT_Mu8_Ele17_CaloIdL_/O)");
  reducedTree_->Branch("passed_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",   &passed_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_,  "passed_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_/O)");


  reducedTree_->Branch("ptHat",&ptHat_,"ptHat_/F");

  
  reducedTree_->Branch("eLept1",  &eLept1_,  "eLept1_/F");
  reducedTree_->Branch("ptLept1",  &ptLept1_,  "ptLept1_/F");
  reducedTree_->Branch("etaLept1",  &etaLept1_,  "etaLept1_/F");
  reducedTree_->Branch("phiLept1",  &phiLept1_,  "phiLept1_/F");
  reducedTree_->Branch("chargeLept1",  &chargeLept1_,  "chargeLept1_/I");
  reducedTree_->Branch("leptTypeLept1",  &leptTypeLept1_,  "leptTypeLept1_/I");
  reducedTree_->Branch("combinedIsoRelLept1",  &combinedIsoRelLept1_,  "combinedIsoRelLept1_/F");

  reducedTree_->Branch("eLept1Gen",  &eLept1Gen_,  "eLept1Gen_/F");
  reducedTree_->Branch("ptLept1Gen",  &ptLept1Gen_,  "ptLept1Gen_/F");
  reducedTree_->Branch("etaLept1Gen",  &etaLept1Gen_,  "etaLept1Gen_/F");
  reducedTree_->Branch("phiLept1Gen",  &phiLept1Gen_,  "phiLept1Gen_/F");

  reducedTree_->Branch("eLept2",  &eLept2_,  "eLept2_/F");
  reducedTree_->Branch("ptLept2",  &ptLept2_,  "ptLept2_/F");
  reducedTree_->Branch("etaLept2",  &etaLept2_,  "etaLept2_/F");
  reducedTree_->Branch("phiLept2",  &phiLept2_,  "phiLept2_/F");
  reducedTree_->Branch("chargeLept2",  &chargeLept2_,  "chargeLept2_/I");
  reducedTree_->Branch("leptTypeLept2",  &leptTypeLept2_,  "leptTypeLept2_/I");
  reducedTree_->Branch("combinedIsoRelLept2",  &combinedIsoRelLept2_,  "combinedIsoRelLept2_/F");

  reducedTree_->Branch("eLept2Gen",  &eLept2Gen_,  "eLept2Gen_/F");
  reducedTree_->Branch("ptLept2Gen",  &ptLept2Gen_,  "ptLept2Gen_/F");
  reducedTree_->Branch("etaLept2Gen",  &etaLept2Gen_,  "etaLept2Gen_/F");
  reducedTree_->Branch("phiLept2Gen",  &phiLept2Gen_,  "phiLept2Gen_/F");

  //reducedTree_->Branch("nLept", &nLept_, "nLept_/I");

  //reducedTree_->Branch("leptTypeLept", leptTypeLept_, "leptTypeLept_[nLept_]/I");
  //reducedTree_->Branch("eLept",  eLept_,  "eLept_[nLept_]/F");
  //reducedTree_->Branch( "ptLept",  ptLept_,  "ptLept_[nLept_]/F");
  //reducedTree_->Branch("etaLept", etaLept_, "etaLept_[nLept_]/F");
  //reducedTree_->Branch("phiLept", phiLept_, "phiLept_[nLept_]/F");
  //reducedTree_->Branch("chargeLept", chargeLept_, "chargeLept_[nLept_]/I");


  reducedTree_->Branch("nJets", &nJets_, "nJets_/I");

  reducedTree_->Branch("eJet",  eJet_,  "eJet_[nJets_]/F");
  reducedTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJets_]/F");
  reducedTree_->Branch("etaJet", etaJet_, "etaJet_[nJets_]/F");
  reducedTree_->Branch("phiJet", phiJet_, "phiJet_[nJets_]/F");

  reducedTree_->Branch("ptDJet", ptDJet_, "ptDJet_[nJets_]/F");
  reducedTree_->Branch("rmsCandJet", rmsCandJet_, "rmsCandJet_[nJets_]/F");
  reducedTree_->Branch("nChargedJet", nChargedJet_, "nChargedJet_[nJets_]/F");
  reducedTree_->Branch("nNeutralJet", nNeutralJet_, "nNeutralJet_[nJets_]/F");
  reducedTree_->Branch("QGLikelihoodJet", QGLikelihoodJet_, "QGLikelihoodJet_[nJets_]/F");

  reducedTree_->Branch("eChargedHadronsJet", eChargedHadronsJet_, "eChargedHadronsJet_[nJets_]/F");
  reducedTree_->Branch("ePhotonsJet", ePhotonsJet_, "ePhotonsJet_[nJets_]/F");
  reducedTree_->Branch("eNeutralHadronsJet", eNeutralHadronsJet_, "eNeutralHadronsJet_[nJets_]/F");
  reducedTree_->Branch("eMuonsJet", eMuonsJet_, "eMuonsJet_[nJets_]/F");
  reducedTree_->Branch("eElectronsJet", eElectronsJet_, "eElectronsJet_[nJets_]/F");
  reducedTree_->Branch("eHFHadronsJet", eHFHadronsJet_, "eHFHadronsJet_[nJets_]/F");
  reducedTree_->Branch("eHFEMJet", eHFEMJet_, "eHFEMJet_[nJets_]/F");

  reducedTree_->Branch("nChargedHadronsJet", nChargedHadronsJet_, "nChargedHadronsJet_[nJets_]/I");
  reducedTree_->Branch("nPhotonsJet", nPhotonsJet_, "nPhotonsJet_[nJets_]/I");
  reducedTree_->Branch("nNeutralHadronsJet", nNeutralHadronsJet_, "nNeutralHadronsJet_[nJets_]/I");
  reducedTree_->Branch("nMuonsJet", nMuonsJet_, "nMuonsJet_[nJets_]/I");
  reducedTree_->Branch("nElectronsJet", nElectronsJet_, "nElectronsJet_[nJets_]/I");
  reducedTree_->Branch("nHFHadronsJet", nHFHadronsJet_, "nHFHadronsJet_[nJets_]/I");
  reducedTree_->Branch("nHFEMJet", nHFEMJet_, "nHFEMJet_[nJets_]/I");

  reducedTree_->Branch("trackCountingHighEffBJetTagJet", trackCountingHighEffBJetTagJet_, "trackCountingHighEffBJetTagJet_[nJets_]/F");
  reducedTree_->Branch("trackCountingHighPurBJetTagJet", trackCountingHighPurBJetTagJet_, "trackCountingHighPurBJetTagJet_[nJets_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighEffBJetTagJet", simpleSecondaryVertexHighEffBJetTagJet_, "simpleSecondaryVertexHighEffBJetTagJet_[nJets_]/F");
  reducedTree_->Branch("simpleSecondaryVertexHighPurBJetTagJet", simpleSecondaryVertexHighPurBJetTagJet_, "simpleSecondaryVertexHighPurBJetTagJet_[nJets_]/F");
  reducedTree_->Branch("jetBProbabilityBJetTagJet", jetBProbabilityBJetTagJet_, "jetBProbabilityBJetTagJet_[nJets_]/F");
  reducedTree_->Branch("jetProbabilityBJetTagJet", jetProbabilityBJetTagJet_, "jetProbabilityBJetTagJet_[nJets_]/F");

  reducedTree_->Branch("eGenJet",  eGenJet_,  "eGenJet_[nJets_]/F");
  reducedTree_->Branch( "ptGenJet",  ptGenJet_,  "ptGenJet_[nJets_]/F");
  reducedTree_->Branch("etaGenJet", etaGenJet_, "etaGenJet_[nJets_]/F");
  reducedTree_->Branch("phiGenJet", phiGenJet_, "phiGenJet_[nJets_]/F");

  reducedTree_->Branch("ePartJet",  ePartJet_,  "ePartJet_[nJets_]/F");
  reducedTree_->Branch( "ptPartJet",  ptPartJet_,  "ptPartJet_[nJets_]/F");
  reducedTree_->Branch("etaPartJet", etaPartJet_, "etaPartJet_[nJets_]/F");
  reducedTree_->Branch("phiPartJet", phiPartJet_, "phiPartJet_[nJets_]/F");
  reducedTree_->Branch("pdgIdPartJet", pdgIdPartJet_, "pdgIdPartJet_[nJets_]/I");


  reducedTree_->Branch("nPart", &nPart_, "nPart_/I");
  reducedTree_->Branch("ePart",  ePart_,  "ePart_[nPart_]/F");
  reducedTree_->Branch( "ptPart",  ptPart_,  "ptPart_[nPart_]/F");
  reducedTree_->Branch("etaPart", etaPart_, "etaPart_[nPart_]/F");
  reducedTree_->Branch("phiPart", phiPart_, "phiPart_[nPart_]/F");
  reducedTree_->Branch("pdgIdPart", pdgIdPart_, "pdgIdPart_[nPart_]/I");
  reducedTree_->Branch("motherPart", motherPart_, "motherPart_[nPart_]/I");


  reducedTree_->Branch("epfMet",&epfMet_,"epfMet_/F");
  reducedTree_->Branch("sumEtpfMet", &sumEtpfMet_,"sumEtpfMet_/F");
  reducedTree_->Branch("metSignificance", &metSignificance_,"metSignificance_/F");
  reducedTree_->Branch("mEtSig", &mEtSig_,"mEtSig_/F");
  reducedTree_->Branch("phipfMet",&phipfMet_,"phipfMet_/F");


  qglikeli_ = new QGLikelihoodCalculator("/afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_3_patch5/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");

} 



Ntp1Analyzer_TTW::~Ntp1Analyzer_TTW() {

  outfile_->cd();
  h1_nCounter_Zee_->Write();
  h1_nCounter_Zmumu_->Write();

  if( qglikeli_!=0 ) delete qglikeli_;

}



void Ntp1Analyzer_TTW::Loop()
{


   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   // to fix Z->ee BR bug, compute number of events having Z->ee and Z->mumu
   int nCounterZee = (isMC_) ? fChain->GetEntries("statusMc==3 && idMc==11") : 0;
   int nCounterZmumu = (isMC_) ? fChain->GetEntries("statusMc==3 && idMc==13") : 0;
   h1_nCounter_Zee_->SetBinContent( 1, nCounterZee );
   h1_nCounter_Zmumu_->SetBinContent( 1, nCounterZmumu );


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;

   // count number of events with PU reweighting:
   std::string puType = "Spring11_Flat10";
   std::string puType_ave = "Spring11_Flat10";
   TString dataset_tstr(dataset_);
   if( dataset_tstr.Contains("Summer11") && dataset_tstr.Contains("PU_S4") ) {
     puType = "Summer11_S4";
     puType_ave = "Summer11_S4_ave";
   }
  // PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
  // PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
  // //PUWeight* fPUWeight = new PUWeight(1089.2, "2011A", puType);
  // TFile* filePU = TFile::Open("Pileup_2011_to_173692_LPLumiScale_68mb.root");
  // TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
  // fPUWeight->SetDataHistogram(h1_nPU_data);
  // fPUWeight_ave->SetDataHistogram(h1_nPU_data);

   float nCounterPU=0.;
   float nCounterPU_ave=0.;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;



     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     genWeight_ = genWeight; //default
     eventWeight_ = -1.; //default
     leptTypeMC_ = -1;




     if( !isGoodEvent(jentry) ) continue; //this takes care also of trigger


     if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     if( !goodVertex ) continue;
  
     nPU_ = nPU[1]; //in time PU only

     nPU_ave_ = 0.;
     for( unsigned iBX=0; iBX<nBX; ++iBX ) {
       nPU_ave_ += nPU[iBX]; 
     }
     nPU_ave_ /= (float)nBX;

     // PU reweighting:
     eventWeightPU_=1.;
     eventWeightPU_ave_=1.;
    // if( isMC_ ) {
    //   eventWeightPU_ = fPUWeight->GetWeight(nPU_);
    //   eventWeightPU_ave_ = fPUWeight_ave->GetWeight(nPU_ave_);
    // }
     nCounterPU += eventWeightPU_;
     nCounterPU_ave += eventWeightPU_ave_;




     nvertex_ = nPV;
     rhoPF_ = rhoFastjet;


     // save trigger info:
     passed_HLT_DoubleMu6_ = this->PassedHLT( jentry, "HLT_DoubleMu6_v");
     passed_HLT_DoubleMu7_ = this->PassedHLT( jentry, "HLT_DoubleMu7_v");
     passed_HLT_Mu13_Mu8_ = this->PassedHLT( jentry, "HLT_Mu13_Mu8_v");
     passed_HLT_Mu17_Mu8_ = this->PassedHLT( jentry, "HLT_Mu17_Mu8_v");
     passed_HLT_IsoMu17_ = this->PassedHLT( jentry, "HLT_IsoMu17_v");
     passed_HLT_IsoMu24_ = this->PassedHLT( jentry, "HLT_IsoMu24_v");
     passed_HLT_IsoMu24_eta2p1_ = this->PassedHLT( jentry, "HLT_IsoMu24_eta2p1_v");
     passed_HLT_Mu8_Jet40_ = this->PassedHLT( jentry, "HLT_Mu8_Jet40_v");
     passed_HLT_L2DoubleMu23_NoVertex_ = this->PassedHLT( jentry, "HLT_L2DoubleMu23_NoVertex_v");
     passed_HLT_L2DoubleMu30_NoVertex_ = this->PassedHLT( jentry, "HLT_L2DoubleMu30_NoVertex_v");
     passed_HLT_TripleMu5_ = this->PassedHLT( jentry, "HLT_TripleMu5_v");
     passed_HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_ = this->PassedHLT( jentry, "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");
     passed_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_ = this->PassedHLT( jentry, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
     passed_HLT_Mu17_Ele8_CaloIdL_ = this->PassedHLT( jentry, "HLT_Mu17_Ele8_CaloIdL_v" );
     passed_HLT_Mu8_Ele17_CaloIdL_ = this->PassedHLT( jentry, "HLT_Mu8_Ele17_CaloIdL_v" );
     passed_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_ = this->PassedHLT( jentry, "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v" );
  



     ptHat_ = (isMC_) ? genPtHat : ptHat_;




     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------


     epfMet_ = energyPFMet[0];
     sumEtpfMet_ = sumEtPFMet[0];
     metSignificance_ = significancePFMet[0];
     mEtSig_ = mEtSigPFMet[0];
     phipfMet_ = phiPFMet[0];


     if( event_==DEBUG_EVENTNUMBER ) {
       std::cout << std::endl << std::endl;
       std::cout << "----- LOG for run: " << run_ << "    event: " << event_ << std::endl;
       std::cout << std::endl << "*** Muons:" << std::endl;
     }


     std::vector<AnalysisLepton> leptPlus;
     std::vector<AnalysisLepton> leptMinus;

     // ------------------
     // MUONS
     // ------------------

     std::vector<AnalysisLepton> muonsPlus;
     std::vector<AnalysisLepton> muonsMinus;
     std::vector<AnalysisLepton> muons;


     //for( unsigned int iMuon=0; iMuon<nMuon && (muons.size()<2); ++iMuon ) {
     for( unsigned int iMuon=0; iMuon<nMuon; ++iMuon ) {

       AnalysisMuon thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );
       thisMuon.charge = chargeMuon[iMuon];
       thisMuon.leptType = 0;

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisMuon.Pt: " << thisMuon.Pt() << std::endl;
         std::cout << "thisMuon.Eta: " << thisMuon.Eta() << std::endl;
         std::cout << "thisMuon.charge: " << thisMuon.charge << std::endl;
       }

       // --------------
       // kinematics:
       // --------------
       if( thisMuon.Pt() < 10. ) continue;
       if( fabs(thisMuon.Eta()) > 2.4 ) continue;

       thisMuon.isGlobalMuonPromptTight = (muonIdMuon[iMuon]>>8)&1;
       thisMuon.isAllTrackerMuon = (muonIdMuon[iMuon]>>11)&1;

       thisMuon.pixelHits = numberOfValidPixelBarrelHitsTrack[trackIndexMuon[iMuon]]+numberOfValidPixelEndcapHitsTrack[trackIndexMuon[iMuon]];
       thisMuon.trackerHits = trackValidHitsTrack[trackIndexMuon[iMuon]];

       thisMuon.nMatchedStations = numberOfMatchesMuon[iMuon];

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisMuon.isGlobalMuonPromptTight: " << thisMuon.isGlobalMuonPromptTight << std::endl;
         std::cout << "thisMuon.isAllTrackerMuon: " << thisMuon.isAllTrackerMuon << std::endl;
         std::cout << "thisMuon.pixelHits: " << thisMuon.pixelHits << std::endl;
         std::cout << "thisMuon.trackerHits: " << thisMuon.trackerHits << std::endl;
         std::cout << "thisMuon.nMatchedStations: " << thisMuon.nMatchedStations << std::endl;
       }


       // to compute dxy, look for primary vertex:
       int hardestPV = -1;
       float sumPtMax = 0.0;
       for(int v=0; v<nPV; v++) {
         if(SumPtPV[v] > sumPtMax) {
           sumPtMax = SumPtPV[v];
           hardestPV = v;
         }
       }  
   
       float dxy;
       if( hardestPV==-1 ) {
         dxy = 0.;
       } else {
         dxy = fabs(trackDxyPV(PVxPV[hardestPV], PVyPV[hardestPV], PVzPV[hardestPV],
                              trackVxTrack[trackIndexMuon[iMuon]], trackVyTrack[trackIndexMuon[iMuon]], trackVzTrack[trackIndexMuon[iMuon]],
                              pxTrack[trackIndexMuon[iMuon]], pyTrack[trackIndexMuon[iMuon]], pzTrack[trackIndexMuon[iMuon]]));
       }


       float dz = fabs(trackVzTrack[trackIndexMuon[iMuon]]-PVzPV[hardestPV]);

       thisMuon.dxy = dxy;
       thisMuon.dz = dz;

       // isolation:
       thisMuon.sumPt03 = sumPt03Muon[iMuon];
       thisMuon.emEt03  = emEt03Muon[iMuon];
       thisMuon.hadEt03 = hadEt03Muon[iMuon];
       thisMuon.isolation = thisMuon.combinedIsoRel();

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisMuon.dxy: " << thisMuon.dxy << std::endl;
         std::cout << "thisMuon.dz: " << thisMuon.dz << std::endl;
         std::cout << "thisMuon.sumPt03: " << thisMuon.sumPt03 << std::endl;
         std::cout << "thisMuon.emEt03: " << thisMuon.emEt03 << std::endl;
         std::cout << "thisMuon.hadEt03: " << thisMuon.hadEt03 << std::endl;
       }

       if( !thisMuon.passedVBTF() ) continue;

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "PASSED VBTF. ";
         if( thisMuon.charge > 0 ) std::cout << "Adding to collection of positive muons." << std::endl;
         else std::cout << "Adding to collection of negative muons." << std::endl;
       }

       if( thisMuon.charge > 0 ) {
         muonsPlus.push_back(thisMuon);
         leptPlus.push_back(thisMuon);
       } else {
         muonsMinus.push_back(thisMuon);
         leptMinus.push_back(thisMuon);
       }

       muons.push_back(thisMuon);

     } //for muons


     // Z veto:
     float ZVetoWindow=15.;
     if( findZ( muonsPlus, muonsMinus, ZVetoWindow ) ) continue;



     // ------------------
     // ELECTRONS
     // ------------------

     std::vector<AnalysisLepton> electronsPlus;
     std::vector<AnalysisLepton> electronsMinus;
     bool firstPassedVBTF80 = false;

     if( event_==DEBUG_EVENTNUMBER )
       std::cout << std::endl << "*** Electrons:" << std::endl;


     //for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {
     for( unsigned int iEle=0; (iEle<nEle); ++iEle ) {

       AnalysisElectron thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );
       thisEle.charge = chargeEle[iEle];
       thisEle.leptType = 1;

       // triple charge consistency:
       int scCharge = scPixChargeEle[iEle];
       int ctfCharge = chargeTrack[trackIndexEle[iEle]];
       int gsfCharge = chargeGsfTrack[gsfTrackIndexEle[iEle]];
       if( (scCharge!=ctfCharge) || (scCharge!=gsfCharge) || (gsfCharge!=ctfCharge) ) continue;

       float scEta = (superClusterIndexEle[iEle]>=0) ? etaSC[superClusterIndexEle[iEle]] : etaPFSC[PFsuperClusterIndexEle[iEle]];

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisEle.Pt: " << thisEle.Pt() << std::endl;
         std::cout << "thisEle.Eta: " << thisEle.Eta() << std::endl;
         std::cout << "thisEle.scEta: " << scEta << std::endl;
         std::cout << "thisEle.charge: " << thisEle.charge << std::endl;
       }

       // --------------
       // kinematics:
       // --------------
       if( thisEle.Pt() < 10. ) continue;
       if( fabs(scEta)>1.4442 && fabs(scEta)<1.566 ) continue; //crack region vetoed with SC eta
       if( fabs(thisEle.Eta()) > 2.5 ) continue; //acceptance cut with electron eta


       thisEle.etaSC = scEta;

       // isolation
       thisEle.dr03TkSumPt = dr03TkSumPtEle[iEle];
       thisEle.dr03EcalRecHitSumEt = dr03EcalRecHitSumEtEle[iEle];
       thisEle.dr03HcalTowerSumEt = dr03HcalTowerSumEtEle[iEle];
       thisEle.isolation = thisEle.combinedIsoRel();

       // electron ID
       thisEle.sigmaIetaIeta = (superClusterIndexEle[iEle]>=0) ? sqrt(covIEtaIEtaSC[superClusterIndexEle[iEle]]) : sqrt(covIEtaIEtaPFSC[PFsuperClusterIndexEle[iEle]]);
       thisEle.deltaPhiAtVtx = deltaPhiAtVtxEle[iEle];
       thisEle.deltaEtaAtVtx = deltaEtaAtVtxEle[iEle];
       thisEle.hOverE = hOverEEle[iEle];
       thisEle.fBrem = fbremEle[iEle];
       thisEle.eOverP = eSuperClusterOverPEle[iEle];

       // conversion rejection
       thisEle.expInnerLayersGsfTrack = expInnerLayersGsfTrack[gsfTrackIndexEle[iEle]];
       thisEle.convDist = convDistEle[iEle];
       thisEle.convDcot = convDcotEle[iEle];

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisEle.dr03TkSumPt: " << thisEle.dr03TkSumPt << std::endl;
         std::cout << "thisEle.dr03EcalRecHitSumEt: " << thisEle.dr03EcalRecHitSumEt << std::endl;
         std::cout << "thisEle.dr03HcalTowerSumEt: " << thisEle.dr03HcalTowerSumEt << std::endl;
         std::cout << "thisEle.sigmaIetaIeta: " << thisEle.sigmaIetaIeta << std::endl;
         std::cout << "thisEle.deltaPhiAtVtx: " << thisEle.deltaPhiAtVtx << std::endl;
         std::cout << "thisEle.deltaEtaAtVtx: " << thisEle.deltaEtaAtVtx << std::endl;
         std::cout << "thisEle.hOverE: " << thisEle.hOverE << std::endl;
         std::cout << "thisEle.expInnerLayersGsfTrack: " << thisEle.expInnerLayersGsfTrack << std::endl;
         std::cout << "thisEle.convDist: " << thisEle.convDist << std::endl;
         std::cout << "thisEle.convDcot: " << thisEle.convDcot << std::endl;
       }



       // electron ID used in AN 2011-466
       bool passed_electronID = thisEle.electronIDVBTF80() 
                             && thisEle.conversionRejectionVBTF80() 
                             && thisEle.combinedIsoRel()<0.15
                             && thisEle.dr03EcalRecHitSumEt/thisEle.Pt()<0.07
                             && thisEle.passedAdditionalCuts();


       if( !passed_electronID ) continue;

       if( event_==DEBUG_EVENTNUMBER ) std::cout << "Passed ElectronID." << std::endl;

      // // additional ID to be as tight as trigger (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL):
      // if( fabs(scEta)<1.4442 ) { //barrel
      //   if( fabs(thisEle.deltaPhiAtVtx) > 0.15 ) continue;
      // } else { //endcaps
      //   if( fabs(thisEle.deltaPhiAtVtx) > 0.1 ) continue;
      //   if( thisEle.hOverE > 0.1 ) continue;
      // }

      // if( event_==DEBUG_EVENTNUMBER ) std::cout << "Passed additional eleID cuts (HLT)." << std::endl;


       // check that not matched to muon (clean electrons faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<AnalysisLepton>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu )
         if( iMu->DeltaR(thisEle)<0.1 ) matchedtomuon=true;

       if( matchedtomuon ) continue;

       if( event_==DEBUG_EVENTNUMBER ) std::cout << "Not matched to any muon." << std::endl;

       if( event_==DEBUG_EVENTNUMBER ) {
         if( thisEle.charge > 0 ) std::cout << "Adding to collection of positive electrons." << std::endl;
         else std::cout << "Adding to collection of negative electrons." << std::endl;
       }

       if( thisEle.charge > 0 ) {
         electronsPlus.push_back(thisEle);
         leptPlus.push_back(thisEle);
       } else {
         electronsMinus.push_back(thisEle);
         leptMinus.push_back(thisEle);
       }

     } //for electrons

     
     // Z veto:
     if( findZ( electronsPlus, electronsMinus, ZVetoWindow ) ) continue;



     // now define lepton candidates
     std::vector<AnalysisLepton> leptons;


     if( leptPlus.size()>=2 && leptMinus.size()<2 ) {

       leptons.push_back( leptPlus[0] );
       leptons.push_back( leptPlus[1] );

     } else if( leptPlus.size()<2 && leptMinus.size()>=2 ) {

       leptons.push_back( leptMinus[0] );
       leptons.push_back( leptMinus[1] );

     } else if( leptPlus.size()>=2 && leptMinus.size()>=2 ) {

       float sumPtPlus = leptPlus[0].Pt()+leptPlus[1].Pt();
       float sumPtMinus = leptMinus[0].Pt()+leptMinus[1].Pt();

       if( sumPtPlus>sumPtMinus ) {

         leptons.push_back( leptPlus[0] );
         leptons.push_back( leptPlus[1] );

       } else {

         leptons.push_back( leptMinus[0] );
         leptons.push_back( leptMinus[1] );

       }

     } else {

       continue;

     }




     if( leptons.size()!=2 ) {
       std::cout << "leptons size is not 2. This shouldnt be possible." << std::endl;
       continue;
     }




     if( leptons[0].Pt()>leptons[1].Pt() ) {

       eLept1_ = leptons[0].Energy();
       ptLept1_ = leptons[0].Pt();
       etaLept1_ = leptons[0].Eta();
       phiLept1_ = leptons[0].Phi();
       chargeLept1_ = leptons[0].charge;
       leptTypeLept1_ = leptons[0].leptType;
       combinedIsoRelLept1_ = leptons[0].isolation;
       
       eLept2_ = leptons[1].Energy();
       ptLept2_ = leptons[1].Pt();
       etaLept2_ = leptons[1].Eta();
       phiLept2_ = leptons[1].Phi();
       chargeLept2_ = leptons[1].charge;
       leptTypeLept2_ = leptons[1].leptType;
       combinedIsoRelLept2_ = leptons[1].isolation;

     } else {

       eLept2_ = leptons[0].Energy();
       ptLept2_ = leptons[0].Pt();
       etaLept2_ = leptons[0].Eta();
       phiLept2_ = leptons[0].Phi();
       chargeLept2_ = leptons[0].charge;
       leptTypeLept2_ = leptons[0].leptType;
       combinedIsoRelLept2_ = leptons[0].isolation;
       
       eLept1_ = leptons[1].Energy();
       ptLept1_ = leptons[1].Pt();
       etaLept1_ = leptons[1].Eta();
       phiLept1_ = leptons[1].Phi();
       chargeLept1_ = leptons[1].charge;
       leptTypeLept1_ = leptons[1].leptType;
       combinedIsoRelLept1_ = leptons[1].isolation;

     }


/*
     // save other leptons:
     nLept_=0;
     for( unsigned iMuonPlus=0; iMuonPlus<muonsPlus.size() && nLept_<10; ++iMuonPlus ) {
       if( iMuonPlus!=iMuonPlus_found ) {
         leptTypeLept_[nLept_] = 0;
         eLept_[nLept_] = muonsPlus[iMuonPlus].Energy();
         ptLept_[nLept_] = muonsPlus[iMuonPlus].Pt();
         etaLept_[nLept_] = muonsPlus[iMuonPlus].Eta();
         phiLept_[nLept_] = muonsPlus[iMuonPlus].Phi();
         chargeLept_[nLept_] = muonsPlus[iMuonPlus].charge;
         nLept_++;
       }
     }
     for( unsigned iMuonMinus=0; iMuonMinus<muonsMinus.size() && nLept_<10; ++iMuonMinus ) {
       if( iMuonMinus!=iMuonMinus_found ) {
         leptTypeLept_[nLept_] = 0;
         eLept_[nLept_] = muonsMinus[iMuonMinus].Energy();
         ptLept_[nLept_] = muonsMinus[iMuonMinus].Pt();
         etaLept_[nLept_] = muonsMinus[iMuonMinus].Eta();
         phiLept_[nLept_] = muonsMinus[iMuonMinus].Phi();
         chargeLept_[nLept_] = muonsMinus[iMuonMinus].charge;
         nLept_++;
       }
     }
     for( unsigned iElectronPlus=0; iElectronPlus<electronsPlus.size() && nLept_<10; ++iElectronPlus ) {
       if( iElectronPlus!=iElePlus_found ) {
         leptTypeLept_[nLept_] = 1;
         eLept_[nLept_] = electronsPlus[iElectronPlus].Energy();
         ptLept_[nLept_] = electronsPlus[iElectronPlus].Pt();
         etaLept_[nLept_] = electronsPlus[iElectronPlus].Eta();
         phiLept_[nLept_] = electronsPlus[iElectronPlus].Phi();
         chargeLept_[nLept_] = electronsPlus[iElectronPlus].charge;
         nLept_++;
       }
     }
     for( unsigned iElectronMinus=0; iElectronMinus<electronsMinus.size() && nLept_<10; ++iElectronMinus ) {
       if( iElectronMinus!=iEleMinus_found ) {
         leptTypeLept_[nLept_] = 1;
         eLept_[nLept_] = electronsMinus[iElectronMinus].Energy();
         ptLept_[nLept_] = electronsMinus[iElectronMinus].Pt();
         etaLept_[nLept_] = electronsMinus[iElectronMinus].Eta();
         phiLept_[nLept_] = electronsMinus[iElectronMinus].Phi();
         chargeLept_[nLept_] = electronsMinus[iElectronMinus].charge;
         nLept_++;
       }
     }
*/


     // ------------------
     // JETS
     // ------------------

     float jetPt_thresh = 15.;

     // first save leading jets in event:
     std::vector<AnalysisJet> leadJets;
     std::vector<int> leadJetsIndex; //index in the event collection (needed afterwards for PFCandidates)

     for( unsigned int iJet=0; iJet<nAK5PFPUcorrJet; ++iJet ) {

       AnalysisJet thisJet( pxAK5PFPUcorrJet[iJet], pyAK5PFPUcorrJet[iJet], pzAK5PFPUcorrJet[iJet], energyAK5PFPUcorrJet[iJet] );

       thisJet.eChargedHadrons = chargedHadronEnergyAK5PFPUcorrJet[iJet];
       thisJet.ePhotons        = photonEnergyAK5PFPUcorrJet[iJet];
       thisJet.eNeutralHadrons = neutralHadronEnergyAK5PFPUcorrJet[iJet];
       thisJet.eElectrons      = electronEnergyAK5PFPUcorrJet[iJet];
       thisJet.eMuons          = muonEnergyAK5PFPUcorrJet[iJet];

       thisJet.nChargedHadrons = chargedHadronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nPhotons        = photonMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nNeutralHadrons = neutralHadronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nElectrons      = electronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nMuons          = muonMultiplicityAK5PFPUcorrJet[iJet];

       thisJet.nCharged = chargedHadronMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.nNeutral = neutralHadronMultiplicityAK5PFPUcorrJet[iJet]+photonMultiplicityAK5PFPUcorrJet[iJet];
       thisJet.rmsCand =  rmsCandAK5PFPUcorrJet[iJet];
       thisJet.ptD =  ptDAK5PFPUcorrJet[iJet];

       thisJet.trackCountingHighEffBJetTag = trackCountingHighEffBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.trackCountingHighPurBJetTag = trackCountingHighPurBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.simpleSecondaryVertexHighEffBJetTag = simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.simpleSecondaryVertexHighPurBJetTag = simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.jetBProbabilityBJetTag = jetBProbabilityBJetTagsAK5PFPUcorrJet[iJet];
       thisJet.jetProbabilityBJetTag = jetProbabilityBJetTagsAK5PFPUcorrJet[iJet];



       //// save at least 3 lead jets (if event has them) and all jets with pt>thresh:
       //if( leadJets.size()>=3 && thisJet.Pt()<jetPt_thresh ) break;

       // far away from leptons:
       bool matchedToLepton=false;
       for( unsigned iMuonPlus=0; iMuonPlus<muonsPlus.size() && !matchedToLepton; ++iMuonPlus )
         if( thisJet.DeltaR( muonsPlus[iMuonPlus] ) <= 0.5 ) matchedToLepton=true;
       for( unsigned iMuonMinus=0; iMuonMinus<muonsMinus.size() && !matchedToLepton; ++iMuonMinus )
         if( thisJet.DeltaR( muonsMinus[iMuonMinus] ) <= 0.5 ) matchedToLepton=true;
       for( unsigned iElectronPlus=0; iElectronPlus<electronsPlus.size() && !matchedToLepton; ++iElectronPlus )
         if( thisJet.DeltaR( electronsPlus[iElectronPlus] ) <= 0.5 ) matchedToLepton=true;
       for( unsigned iElectronMinus=0; iElectronMinus<electronsMinus.size() && !matchedToLepton; ++iElectronMinus )
         if( thisJet.DeltaR( electronsMinus[iElectronMinus] ) <= 0.5 ) matchedToLepton=true;
       if( matchedToLepton ) continue;
       

       // jet ID:
       int multiplicity = thisJet.nCharged +  thisJet.nNeutral + HFEMMultiplicityAK5PFPUcorrJet[iJet] + HFHadronMultiplicityAK5PFPUcorrJet[iJet];
       if( multiplicity < 2 ) continue;
       if( fabs(thisJet.Eta())<2.4 && thisJet.nChargedHadrons == 0 ) continue;
       if( thisJet.eNeutralHadrons >= 0.99*thisJet.Energy() ) continue;
       if( thisJet.ePhotons >= 0.99*thisJet.Energy() ) continue;

       // match to genjet:
       float bestDeltaR=999.;
       TLorentzVector matchedGenJet;
       for( unsigned iGenJet=0; iGenJet<nAK5GenJet; ++iGenJet ) {
         TLorentzVector thisGenJet(pxAK5GenJet[iGenJet], pyAK5GenJet[iGenJet], pzAK5GenJet[iGenJet], energyAK5GenJet[iGenJet]);
         if( thisGenJet.DeltaR(thisJet) < bestDeltaR ) {
           bestDeltaR=thisGenJet.DeltaR(thisJet);
           matchedGenJet=thisGenJet;
         }
       }

       thisJet.ptGen  = (isMC_ && matchedGenJet.Pt()>0.) ? matchedGenJet.Pt() : 0.;
       thisJet.etaGen = (isMC_ && matchedGenJet.Pt()>0.) ? matchedGenJet.Eta() : 20.;
       thisJet.phiGen = (isMC_ && matchedGenJet.Pt()>0.) ? matchedGenJet.Phi() : 0.;
       thisJet.eGen   = (isMC_ && matchedGenJet.Pt()>0.) ? matchedGenJet.Energy() : 0.;

       // match to parton:
       float bestDeltaR_part=999.;
       TLorentzVector matchedPart;
       int pdgIdPart=0;
       for( unsigned iPart=0; iPart<nMc; ++iPart ) {
         if( statusMc[iPart]!=3 ) continue; //partons
         if( idMc[iPart]!=21 && abs(idMc[iPart])>6 ) continue; //quarks or gluons
         if( pMc[iPart]*sin(thetaMc[iPart])<0.1 ) continue; 
         TLorentzVector thisPart;
         thisPart.SetPtEtaPhiE(pMc[iPart]*sin(thetaMc[iPart]), etaMc[iPart], phiMc[iPart], energyMc[iPart]);
         if( thisPart.Pt() < 0.1 ) continue;
         if( thisPart.DeltaR(thisJet) < bestDeltaR_part ) {
           bestDeltaR_part=thisPart.DeltaR(thisJet);
           matchedPart=thisPart;
           pdgIdPart=idMc[iPart];
         }
       }
   
       thisJet.ptPart  = (isMC_ && matchedPart.Pt()>0.) ? matchedPart.Pt() : 0.;
       thisJet.etaPart = (isMC_ && matchedPart.Pt()>0.) ? matchedPart.Eta() : 20.;
       thisJet.phiPart = (isMC_ && matchedPart.Pt()>0.) ? matchedPart.Phi() : 0.;
       thisJet.ePart   = (isMC_ && matchedPart.Pt()>0.) ? matchedPart.Energy() : 0.;
       thisJet.pdgIdPart   = pdgIdPart;

       
       leadJets.push_back(thisJet);
       leadJetsIndex.push_back(iJet);

     }


     if( leadJets.size()<2 ) continue;
     if( leadJets[1].Pt()<jetPt_thresh ) continue; //at least 2 jets over thresh



     nJets_ = 0;
     nPart_ = 0;


     for( unsigned iJet=0; iJet<leadJets.size() && nJets_<50; ++iJet ) {

       if( leadJets[iJet].Pt()<jetPt_thresh ) continue;

       eJet_[nJets_] = leadJets[iJet].Energy();
       ptJet_[nJets_] = leadJets[iJet].Pt();
       etaJet_[nJets_] = leadJets[iJet].Eta();
       phiJet_[nJets_] = leadJets[iJet].Phi();
       eChargedHadronsJet_[nJets_] = leadJets[iJet].eChargedHadrons;
       ePhotonsJet_[nJets_]        = leadJets[iJet].ePhotons;
       eNeutralHadronsJet_[nJets_] = leadJets[iJet].eNeutralHadrons;
       eElectronsJet_[nJets_]      = leadJets[iJet].eElectrons;
       eMuonsJet_[nJets_]          = leadJets[iJet].eMuons;
       nChargedHadronsJet_[nJets_] = leadJets[iJet].nChargedHadrons;
       nPhotonsJet_[nJets_]        = leadJets[iJet].nPhotons;
       nNeutralHadronsJet_[nJets_] = leadJets[iJet].nNeutralHadrons;
       nElectronsJet_[nJets_]      = leadJets[iJet].nElectrons;
       nMuonsJet_[nJets_]          = leadJets[iJet].nMuons;

       ptDJet_[nJets_] = leadJets[iJet].ptD;
       rmsCandJet_[nJets_] = leadJets[iJet].rmsCand;
       nChargedJet_[nJets_] = leadJets[iJet].nCharged;
       nNeutralJet_[nJets_] = leadJets[iJet].nNeutral;

       QGLikelihoodJet_[nJets_] = qglikeli_->computeQGLikelihoodPU( leadJets[iJet].Pt(), rhoPF_, leadJets[iJet].nCharged, leadJets[iJet].nNeutral, leadJets[iJet].ptD );

       trackCountingHighEffBJetTagJet_[nJets_] = leadJets[iJet].trackCountingHighEffBJetTag;
       trackCountingHighPurBJetTagJet_[nJets_] = leadJets[iJet].trackCountingHighPurBJetTag;
       simpleSecondaryVertexHighEffBJetTagJet_[nJets_] = leadJets[iJet].simpleSecondaryVertexHighEffBJetTag;
       simpleSecondaryVertexHighPurBJetTagJet_[nJets_] = leadJets[iJet].simpleSecondaryVertexHighPurBJetTag;
       jetBProbabilityBJetTagJet_[nJets_] = leadJets[iJet].jetBProbabilityBJetTag;
       jetProbabilityBJetTagJet_[nJets_] = leadJets[iJet].jetProbabilityBJetTag;


       eGenJet_[nJets_] = leadJets[iJet].eGen;
       ptGenJet_[nJets_] = leadJets[iJet].ptGen;
       etaGenJet_[nJets_] = leadJets[iJet].etaGen;
       phiGenJet_[nJets_] = leadJets[iJet].phiGen;
        
       ePartJet_[nJets_] = leadJets[iJet].ePart;
       ptPartJet_[nJets_] = leadJets[iJet].ptPart;
       etaPartJet_[nJets_] = leadJets[iJet].etaPart;
       phiPartJet_[nJets_] = leadJets[iJet].phiPart;
       pdgIdPartJet_[nJets_] = leadJets[iJet].pdgIdPart;

       nJets_++;

     }


     //// at least 4 jets in the event
     //if( nJets_<4 ) continue;



     if( isMC_ ) {

       // store event partons in tree:
       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         //if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {
         if( statusMc[iMc]==3 && pMc[iMc]*sin(thetaMc[iMc])>0.1 ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

           if( nPart_<20 ) {

             ptPart_[nPart_] = thisParticle->Pt();
             etaPart_[nPart_] = thisParticle->Eta();
             phiPart_[nPart_] = thisParticle->Phi();
             ePart_[nPart_] = thisParticle->Energy();
             pdgIdPart_[nPart_] = idMc[iMc];
             motherPart_[nPart_] = idMc[mothMc[iMc]];

             nPart_++;

           } else {
      
             std::cout << "Found more than 20 partons, skipping." << std::endl;

           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

     } // if is mc


     reducedTree_->Fill(); 


   } //for entries

   h1_nCounterPU_->SetBinContent( 1, nCounterPU );
   h1_nCounterPU_ave_->SetBinContent( 1, nCounterPU_ave );

} //loop



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}




/*
std::vector<AnalysisLepton> findLeptons( std::vector<AnalysisLepton> lPlus, std::vector<AnalysisLepton> lMinus ) {

  float maxPt = 0.;
  int iLept_found = -1;
  AnalysisLepton hardestLepton;

  for( unsigned iLept=0; iLept<lPlus.size(); ++iLept ) {
    if( lPlus[iLept].Pt()>maxPt ) {
      maxPt=lPlus[iLept].Pt();
      iLept_found = iLept;
      hardestLepton = lPlus[iLept];
    }
  } // for lPlus

  for( unsigned iLept=0; iLept<lMinus.size(); ++iLept ) {
    if( lMinus[iLept].Pt()>maxPt ) {
      maxPt=lMinus[iLept].Pt();
      iLept_found = iLept;
      hardestLepton = lMinus[iLept];
    }
  } // for lMinus


  if( iLept_found<0 ) {
    std::cout << "Didn't find leptons!! Exiting!" << std::endl;
    exit(33);
  }

  std::vector<AnalysisLepton> leptons;
  leptons.push_back( hardestLepton );

  if( hardestLepton.charge<0 ) {
  
    for( unsigned iLept=0; iLept<lMinus.size() && leptons.size()==1; ++iLept ) {
      if( iLept==iLept_found ) continue;
      leptons.push_back(lMinus[iLept]);
    } // for lMinus
    
    for( unsigned iLept=0; iLept<lPlus.size() && leptons.size()==1; ++iLept ) {
      leptons.push_back(lPlus[iLept]);
    } // for lPlus
    
  } else {
  
    for( unsigned iLept=0; iLept<lPlus.size() && leptons.size()==1; ++iLept ) {
      if( iLept==iLept_found ) continue;
      leptons.push_back(lPlus[iLept]);
    } // for lPlus
    
    for( unsigned iLept=0; iLept<lMinus.size() && leptons.size()==1; ++iLept ) {
      leptons.push_back(lMinus[iLept]);
    } // for lMinus

  }


  std::cout << "leptons size: " << leptons.size() << std::endl;
  std::cout << "pts: " << leptons[0].Pt() << " " << leptons[1].Pt() << std::endl;

  return leptons;

}
*/



bool findZ( std::vector<AnalysisLepton> leptPlus, std::vector<AnalysisLepton> leptMinus, float ZVetoWindow ) {

  bool foundZ = false;

  for( unsigned iLeptPlus=0; iLeptPlus<leptPlus.size() && !foundZ; ++iLeptPlus ) {

    for( unsigned iLeptMinus=0; iLeptMinus<leptMinus.size() && !foundZ; ++iLeptMinus ) {

      TLorentzVector diLepton = leptPlus[iLeptPlus] + leptMinus[iLeptMinus];
   
      if( fabs(diLepton.M()-91.19)<ZVetoWindow ) foundZ = true;

    } // for minus

  } // for plus

  return foundZ;

}
