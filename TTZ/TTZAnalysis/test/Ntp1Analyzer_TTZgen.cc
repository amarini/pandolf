#include "Ntp1Analyzer_TTZgen.h"


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


int DEBUG_EVENTNUMBER = 157480550;
float mZ = 91.1876;


double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);
float getWeightPU(Int_t nPU);
std::vector<AnalysisLepton> getBestZMassPair( const std::vector<AnalysisLepton>& leptPlus, const std::vector<AnalysisLepton>& leptMinus );






Ntp1Analyzer_TTZgen::Ntp1Analyzer_TTZgen( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "TTZgen", dataset, flags, tree ) {



} //constructor



void Ntp1Analyzer_TTZgen::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("eZllMC",  &eZllMC_,  "eZllMC_/F");
  reducedTree_->Branch("ptZllMC",  &ptZllMC_,  "ptZllMC_/F");
  reducedTree_->Branch("etaZllMC",  &etaZllMC_,  "etaZllMC_/F");
  reducedTree_->Branch("phiZllMC",  &phiZllMC_,  "phiZllMC_/F");
  
  reducedTree_->Branch("eleptZ1MC",  &eleptZ1MC_,  "eleptZ1MC_/F");
  reducedTree_->Branch("ptleptZ1MC",  &ptleptZ1MC_,  "ptleptZ1MC_/F");
  reducedTree_->Branch("etaleptZ1MC",  &etaleptZ1MC_,  "etaleptZ1MC_/F");
  reducedTree_->Branch("phileptZ1MC",  &phileptZ1MC_,  "phileptZ1MC_/F");
  
  reducedTree_->Branch("eleptZ2MC",  &eleptZ2MC_,  "eleptZ2MC_/F");
  reducedTree_->Branch("ptleptZ2MC",  &ptleptZ2MC_,  "ptleptZ2MC_/F");
  reducedTree_->Branch("etaleptZ2MC",  &etaleptZ2MC_,  "etaleptZ2MC_/F");
  reducedTree_->Branch("phileptZ2MC",  &phileptZ2MC_,  "phileptZ2MC_/F");
  
  reducedTree_->Branch("et1",  &et1_,  "et1_/F");
  reducedTree_->Branch("ptt1",  &ptt1_,  "ptt1_/F");
  reducedTree_->Branch("etat1",  &etat1_,  "etat1_/F");
  reducedTree_->Branch("phit1",  &phit1_,  "phit1_/F");

  reducedTree_->Branch("et2",  &et2_,  "et2_/F");
  reducedTree_->Branch("ptt2",  &ptt2_,  "ptt2_/F");
  reducedTree_->Branch("etat2",  &etat2_,  "etat2_/F");
  reducedTree_->Branch("phit2",  &phit2_,  "phit2_/F");


  reducedTree_->Branch("nJets", &nJets_, "nJets_/I");

  reducedTree_->Branch("eJet",  eJet_,  "eJet_[nJets_]/F");
  reducedTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJets_]/F");
  reducedTree_->Branch("etaJet", etaJet_, "etaJet_[nJets_]/F");
  reducedTree_->Branch("phiJet", phiJet_, "phiJet_[nJets_]/F");


  

} 



Ntp1Analyzer_TTZgen::~Ntp1Analyzer_TTZgen() {

  outfile_->cd();

}



void Ntp1Analyzer_TTZgen::Loop()
{


   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;



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
   PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
   PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
   //PUWeight* fPUWeight = new PUWeight(1089.2, "2011A", puType);
   TFile* filePU = TFile::Open("Pileup_2011_to_173692_LPLumiScale_68mb.root");
   TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
   fPUWeight->SetDataHistogram(h1_nPU_data);
   fPUWeight_ave->SetDataHistogram(h1_nPU_data);




   float nCounterPU=0.;
   float nCounterPU_ave=0.;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%100000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;



     leptType_ = -1;



     if( !isGoodEvent(jentry) ) continue; //this takes care also of trigger


//   if( nPV==0 ) continue;
//   bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
//   if( !goodVertex ) continue;
//
//   nPU_ave_ = 0.;
//   for( unsigned iBX=0; iBX<nBX; ++iBX ) {
//     nPU_ave_ += nPU[iBX]; 
//   }
//   nPU_ave_ /= (float)nBX;

//   // PU reweighting:
//   eventWeightPU_=1.;
//   eventWeightPU_ave_=1.;
//   if( isMC_ ) {
//     eventWeightPU_ = fPUWeight->GetWeight(nPU_);
//     eventWeightPU_ave_ = fPUWeight_ave->GetWeight(nPU_ave_);
//   }
//   nCounterPU += eventWeightPU_;
//   nCounterPU_ave += eventWeightPU_ave_;






     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     int zIndexqq=-1;
     int zIndexll=-1;


     if( isMC_ ) {


       // look for Z->ll

       std::vector<TLorentzVector> electronsMC;
       std::vector<TLorentzVector> muonsMC;
       std::vector<TLorentzVector> tops;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         // partons only
         if( statusMc[iMc] != 3 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       
         if( idMc[mothMc[iMc]]==23 ) {
           zIndexll = mothMc[iMc]; 
           if( fabs(idMc[iMc])==11 && idMc[mothMc[iMc]]==23 ) electronsMC.push_back( *thisParticle );
           if( fabs(idMc[iMc])==13 && idMc[mothMc[iMc]]==23 ) muonsMC.push_back( *thisParticle );
         }

         if( abs(idMc[mothMc[iMc]])==6 )  //top
           tops.push_back(*thisParticle);

         delete thisParticle;
         thisParticle = 0;

       }



       if( tops.size()!=2 ) {

         std::cout << "There must be a problem, didnt find 2 tops! (found " << tops.size() << ")" << std::endl;
         exit(77);

       } else {

         if( tops[0].Pt() > tops[1].Pt() ) {
      
           ptt1_  = tops[0].Pt();
           et1_   = tops[0].Energy();
           etat1_ = tops[0].Eta();
           phit1_ = tops[0].Phi();

           ptt2_  = tops[1].Pt();
           et2_   = tops[1].Energy();
           etat2_ = tops[1].Eta();
           phit2_ = tops[1].Phi();

         } else {
      
           ptt1_  = tops[1].Pt();
           et1_   = tops[1].Energy();
           etat1_ = tops[1].Eta();
           phit1_ = tops[1].Phi();

           ptt2_  = tops[0].Pt();
           et2_   = tops[0].Energy();
           etat2_ = tops[0].Eta();
           phit2_ = tops[0].Phi();

         }

       } //if found tops


       if( electronsMC.size()==2 ) {
         if( electronsMC[0].Pt() > electronsMC[1].Pt() ) {
           lept1MC = electronsMC[0];
           lept2MC = electronsMC[1];
         } else {
           lept1MC = electronsMC[1];
           lept2MC = electronsMC[0];
         }
       } else if( muonsMC.size()==2 ) {
         if( muonsMC[0].Pt() > muonsMC[1].Pt() ) {
           lept1MC = muonsMC[0];
           lept2MC = muonsMC[1];
         } else {
           lept1MC = muonsMC[1];
           lept2MC = muonsMC[0];
         }
       } else {
         //taus
         noLeptons = true;
       }



       if( !noLeptons ) {

         TLorentzVector ZllMC;
         ZllMC.SetPtEtaPhiE( pMc[zIndexll]*sin(thetaMc[zIndexll]), etaMc[zIndexll], phiMc[zIndexll], energyMc[zIndexll] );

         ptZllMC_  = ZllMC.Pt();
         eZllMC_   = ZllMC.Energy();
         etaZllMC_ = ZllMC.Eta();
         phiZllMC_ = ZllMC.Phi();

         if( muonsMC.size() > 0 ) {
      
           leptType_ = 0;

         } else if( electronsMC.size() > 0 ) {

           leptType_ = 1;

         }

         ptleptZ1MC_  = lept1MC.Pt();
         eleptZ1MC_   = lept1MC.Energy();
         etaleptZ1MC_ = lept1MC.Eta();
         phileptZ1MC_ = lept1MC.Phi();

         ptleptZ2MC_  = lept2MC.Pt();
         eleptZ2MC_   = lept2MC.Energy();
         etaleptZ2MC_ = lept2MC.Eta();
         phileptZ2MC_ = lept2MC.Phi();

       }


//     // now look for the higgs:
//     if( zIndexll!=-1 && zIndexqq!=-1 ) {

//       int higgsIndex = mothMc[zIndexll];

//       if( idMc[higgsIndex] == 25 ) {

//         TLorentzVector HiggsMC;
//         HiggsMC.SetPtEtaPhiE( pMc[higgsIndex]*sin(thetaMc[higgsIndex]), etaMc[higgsIndex], phiMc[higgsIndex], energyMc[higgsIndex] );

//         eHiggsMC_   = HiggsMC.Energy(); 
//         ptHiggsMC_  = HiggsMC.Pt(); 
//         etaHiggsMC_ = HiggsMC.Eta(); 
//         phiHiggsMC_ = HiggsMC.Phi(); 

//       } // if higgs

//     } //if found two Z's

//   } //if isMC


     nJets_ = 0;

     for( unsigned iGenJet=0; iGenJet<nAK5GenJet && nJets_<50; ++iGenJet ) {


       TLorentzVector thisGenJet(pxAK5GenJet[iGenJet], pyAK5GenJet[iGenJet], pzAK5GenJet[iGenJet], energyAK5GenJet[iGenJet]);

       if( thisGenJet.DeltaR(lept1MC) < 0.3 || thisGenJet.DeltaR(lept2MC) < 0.3 ) continue;

       // match to parton
       float bestDeltaR = 999.;
       int genJetPartonID_found = -999;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {
       
         //if( statusMc[iMc]==3 && (fabs(idMc[iMc])<=6 || idMc[iMc]==21) ) {
         if( statusMc[iMc]!=3 ) continue;
         if( pMc[iMc]*sin(thetaMc[iMc])<0.1 ) continue;

         int thisPdgId = idMc[iMc];
         if( thisPdgId>6 && thisPdgId!=21 ) continue; //quarks and gluons
       
         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

         float thisDeltaR = thisGenJet.DeltaR(*thisParticle);

         if( thisDeltaR < bestDeltaR ) {

           bestDeltaR=thisDeltaR;
           genJetPartonID_found=thisPdgId;

         } //if

       } //loop on partons


       eJet_[nJets_] = thisGenJet.Energy();
       ptJet_[nJets_] = thisGenJet.Pt();
       etaJet_[nJets_] = thisGenJet.Eta();
       phiJet_[nJets_] = thisGenJet.Phi();

       pdgIdJet_[nJets_] = genJetPartonID_found;

       nJets_++; 

     } //loop on genjets



/*
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

     // ------------------
     // MUONS
     // ------------------

     std::vector<AnalysisLepton> muonsPlus;
     std::vector<AnalysisLepton> muonsMinus;
     int chargeFirstMuon;


     //for( unsigned int iMuon=0; iMuon<nMuon && (muons.size()<2); ++iMuon ) {
     for( unsigned int iMuon=0; iMuon<nMuon; ++iMuon ) {

       AnalysisMuon thisMuon( pxMuon[iMuon], pyMuon[iMuon], pzMuon[iMuon], energyMuon[iMuon] );
       thisMuon.charge = chargeMuon[iMuon];

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

       thisMuon.dxy = 0.;
       thisMuon.dz = 0.;
       //thisMuon.dxy = dxy;
       //thisMuon.dz = dz;

h1_dxy->Fill( thisMuon.dxy );
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

       if( thisMuon.charge > 0 ) muonsPlus.push_back(thisMuon);
       else muonsMinus.push_back(thisMuon);

     } //for muons


     std::vector<AnalysisLepton> muons = getBestZMassPair( muonsPlus, muonsMinus );



     // ------------------
     // ELECTRONS
     // ------------------

     std::vector<AnalysisLepton> electronsPlus;
     std::vector<AnalysisLepton> electronsMinus;
     int chargeFirstEle = 0;
     bool firstPassedVBTF80 = false;

     if( event_==DEBUG_EVENTNUMBER )
       std::cout << std::endl << "*** Electrons:" << std::endl;


     //for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {
     for( unsigned int iEle=0; (iEle<nEle); ++iEle ) {

       AnalysisElectron thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );
       thisEle.charge = chargeEle[iEle];

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


       bool passed_VBTF95 = thisEle.passedVBTF95();
       bool passed_VBTF80 = thisEle.passedVBTF80();


       //if( !passed_VBTF95 ) continue;
       if( !passed_VBTF80 ) continue;

       if( event_==DEBUG_EVENTNUMBER ) std::cout << "Passed VBTF95." << std::endl;

       // additional ID to be as tight as trigger (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL):
       if( fabs(scEta)<1.4442 ) { //barrel
         if( fabs(thisEle.deltaPhiAtVtx) > 0.15 ) continue;
         if( thisEle.hOverE > 0.12 ) continue; //conforming to SSDL analysis
       } else { //endcaps
         if( fabs(thisEle.deltaPhiAtVtx) > 0.1 ) continue;
         if( thisEle.hOverE > 0.15 ) continue; // looks like (from ntuples) there's a cut at 0.15 at hlt, not 0.075 as CaloIdT says
       }

       if( event_==DEBUG_EVENTNUMBER ) std::cout << "Passed additional eleID cuts (HLT)." << std::endl;


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

       if( thisEle.charge > 0 ) electronsPlus.push_back(thisEle);
       else electronsMinus.push_back(thisEle);

     } //for electrons


     std::vector<AnalysisLepton> electrons = getBestZMassPair( electronsPlus, electronsMinus );


     if( event_==DEBUG_EVENTNUMBER ) 
       std::cout << "Found: " << muons.size() << " muons and " << electrons.size() << " electrons." << std::endl;



     std::vector< AnalysisLepton > leptons;

     // definition of leptType:
     // 0: mumu
     // 1: ee
     // 2: e+mu-
     // 3: e-mu+


     if( electrons.size() < 2 && muons.size() < 2 ) { //this is the ttbar opposite flavour control region:

       if( event_==DEBUG_EVENTNUMBER ) 
         std::cout << "Going in emu control region." << std::endl;

       // at least one opposite-sign pair:
       if( (electronsPlus.size()+muonsPlus.size())==0 || (electronsMinus.size()+muonsMinus.size())==0 ) continue;


       std::vector<AnalysisLepton> elePlus_muMinus = getBestZMassPair( electronsPlus, muonsMinus );
       std::vector<AnalysisLepton> eleMinus_muPlus = getBestZMassPair( electronsMinus, muonsPlus );

       if( elePlus_muMinus.size() == 2 && eleMinus_muPlus.size() == 2 ) { 

         TLorentzVector Zepmm = elePlus_muMinus[0] + elePlus_muMinus[1];
         TLorentzVector Zemmp = eleMinus_muPlus[0] + eleMinus_muPlus[1];

         if( fabs(Zepmm.M()-mZ) < fabs(Zemmp.M()-mZ) ) {

           if( event_==DEBUG_EVENTNUMBER ) 
             std::cout << "LeptType=2." << std::endl;

           leptType_=2;
           if( elePlus_muMinus[0].Pt() > elePlus_muMinus[1].Pt() ) {
             leptons.push_back( elePlus_muMinus[0] );
             leptons.push_back( elePlus_muMinus[1] );
           } else {
             leptons.push_back( elePlus_muMinus[1] );
             leptons.push_back( elePlus_muMinus[0] );
           }

         } else { 

           if( event_==DEBUG_EVENTNUMBER ) 
             std::cout << "LeptType=3." << std::endl;

           leptType_=3;
           if( eleMinus_muPlus[0].Pt() > eleMinus_muPlus[1].Pt() ) {
             leptons.push_back( eleMinus_muPlus[0] );
             leptons.push_back( eleMinus_muPlus[1] );
           } else {
             leptons.push_back( eleMinus_muPlus[1] );
             leptons.push_back( eleMinus_muPlus[0] );
           }

         }


       } else if( elePlus_muMinus.size() == 2 ) {

         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "LeptType=2." << std::endl;

         leptType_ = 2;

         if( elePlus_muMinus[0].Pt() > elePlus_muMinus[1].Pt() ) {

           leptons.push_back( elePlus_muMinus[0] );
           leptons.push_back( elePlus_muMinus[1] );

         } else {

           leptons.push_back( elePlus_muMinus[1] );
           leptons.push_back( elePlus_muMinus[0] );

         }

       } else if( eleMinus_muPlus.size() == 2 ) {

         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "LeptType=3." << std::endl;

         leptType_ = 3;

         if( eleMinus_muPlus[0].Pt() > eleMinus_muPlus[1].Pt() ) {

           leptons.push_back( eleMinus_muPlus[0] );
           leptons.push_back( eleMinus_muPlus[1] );

         } else {

           leptons.push_back( eleMinus_muPlus[1] );
           leptons.push_back( eleMinus_muPlus[0] );

         }

       } else {

         std::cout << "There must be an error this is not possible." << std::endl;
         exit(9101);

       }


     } else { // this is the (same flavor) signal region

       if( event_==DEBUG_EVENTNUMBER ) 
         std::cout << "This is the same flavor region." << std::endl;

       if( electrons.size() == 2 && muons.size() == 2 ) { 

         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "Going to choose pair with best Z mass." << std::endl;

         TLorentzVector Zee = electrons[0] + electrons[1];
         TLorentzVector Zmm = muons[0] + muons[1];

         if( fabs(Zee.M()-mZ) < fabs(Zmm.M()-mZ) ) {

           if( event_==DEBUG_EVENTNUMBER ) 
             std::cout << "Chose electron pair. LeptType=1." << std::endl;

           leptType_=1;
           if( electrons[0].Pt() > electrons[1].Pt() ) {
             leptons.push_back( electrons[0] );
             leptons.push_back( electrons[1] );
           } else {
             leptons.push_back( electrons[1] );
             leptons.push_back( electrons[0] );
           }

         } else { 

           if( event_==DEBUG_EVENTNUMBER ) 
             std::cout << "Chose muon pair. LeptType=0." << std::endl;

           leptType_=0;
           if( muons[0].Pt() > muons[1].Pt() ) {
             leptons.push_back( muons[0] );
             leptons.push_back( muons[1] );
           } else {
             leptons.push_back( muons[1] );
             leptons.push_back( muons[0] );
           }

         }


       } else if( electrons.size() == 2 ) {

         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "LeptType=1." << std::endl;

         leptType_ = 1;

         if( electrons[0].Pt() > electrons[1].Pt() ) {

           leptons.push_back( electrons[0] );
           leptons.push_back( electrons[1] );

         } else {

           leptons.push_back( electrons[1] );
           leptons.push_back( electrons[0] );

         }

       } else if( muons.size() == 2 ) {

         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "LeptType=0." << std::endl;

         leptType_ = 0;

         if( muons[0].Pt() > muons[1].Pt() ) {

           leptons.push_back( muons[0] );
           leptons.push_back( muons[1] );

         } else {

           leptons.push_back( muons[1] );
           leptons.push_back( muons[0] );

         }

       } else {

         std::cout << "There must be an error this is not possible." << std::endl;
         exit(9101);

       }

     } // if opposite flavour control region

 
     if( leptons.size()<2 ) continue;

  

     eLeptZ1_ = leptons[0].Energy();
     ptLeptZ1_ = leptons[0].Pt();
     etaLeptZ1_ = leptons[0].Eta();
     phiLeptZ1_ = leptons[0].Phi();
     chargeLeptZ1_ = leptons[0].charge;
     combinedIsoRelLeptZ1_ = leptons[0].isolation;
     matchedToHLTLeptZ1_ = isMatchedToHLT(leptons[0].Eta(),leptons[0].Phi(),0.3);
     
     eLeptZ2_ = leptons[1].Energy();
     ptLeptZ2_ = leptons[1].Pt();
     etaLeptZ2_ = leptons[1].Eta();
     phiLeptZ2_ = leptons[1].Phi();
     chargeLeptZ2_ = leptons[1].charge;
     combinedIsoRelLeptZ2_ = leptons[1].isolation;
     matchedToHLTLeptZ2_ = isMatchedToHLT(leptons[1].Eta(),leptons[1].Phi(),0.3);


     
     if( event_==DEBUG_EVENTNUMBER ) 
       std::cout << "Now adding in other leptons." << std::endl;

     // save other leptons:
     nLept_=0;
     for( unsigned iMuonPlus=0; iMuonPlus<muonsPlus.size() && nLept_<10; ++iMuonPlus ) {
       bool notSelected = (muonsPlus[iMuonPlus]!=leptons[0] && muonsPlus[iMuonPlus]!=leptons[1]); 
       if( notSelected ) {
         leptTypeLept_[nLept_] = 0;
         eLept_[nLept_] = muonsPlus[iMuonPlus].Energy();
         ptLept_[nLept_] = muonsPlus[iMuonPlus].Pt();
         etaLept_[nLept_] = muonsPlus[iMuonPlus].Eta();
         phiLept_[nLept_] = muonsPlus[iMuonPlus].Phi();
         chargeLept_[nLept_] = muonsPlus[iMuonPlus].charge;
         combinedIsoRelLept_[nLept_] = muonsPlus[iMuonPlus].isolation;
         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "Added a positive muon: pt=" << ptLept_[nLept_] << ", eta=" << etaLept_[nLept_] << std::endl;
         nLept_++;
       }
     }
     for( unsigned iMuonMinus=0; iMuonMinus<muonsMinus.size() && nLept_<10; ++iMuonMinus ) {
       bool notSelected = (muonsMinus[iMuonMinus]!=leptons[0] && muonsMinus[iMuonMinus]!=leptons[1]);
       if( notSelected ) {
         leptTypeLept_[nLept_] = 0;
         eLept_[nLept_] = muonsMinus[iMuonMinus].Energy();
         ptLept_[nLept_] = muonsMinus[iMuonMinus].Pt();
         etaLept_[nLept_] = muonsMinus[iMuonMinus].Eta();
         phiLept_[nLept_] = muonsMinus[iMuonMinus].Phi();
         chargeLept_[nLept_] = muonsMinus[iMuonMinus].charge;
         combinedIsoRelLept_[nLept_] = muonsMinus[iMuonMinus].isolation;
         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "Added a negative muon: pt=" << ptLept_[nLept_] << ", eta=" << etaLept_[nLept_] << std::endl;
         nLept_++;
       }
     }
     for( unsigned iElectronPlus=0; iElectronPlus<electronsPlus.size() && nLept_<10; ++iElectronPlus ) {
       bool notSelected = ( electronsPlus[iElectronPlus]!=leptons[0] && electronsPlus[iElectronPlus]!=leptons[1] );
       if( notSelected ) {
         leptTypeLept_[nLept_] = 1;
         eLept_[nLept_] = electronsPlus[iElectronPlus].Energy();
         ptLept_[nLept_] = electronsPlus[iElectronPlus].Pt();
         etaLept_[nLept_] = electronsPlus[iElectronPlus].Eta();
         phiLept_[nLept_] = electronsPlus[iElectronPlus].Phi();
         chargeLept_[nLept_] = electronsPlus[iElectronPlus].charge;
         combinedIsoRelLept_[nLept_] = electronsPlus[iElectronPlus].isolation;
         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "Added a positive electron: pt=" << ptLept_[nLept_] << ", eta=" << etaLept_[nLept_] << std::endl;
         nLept_++;
       }
     }
     for( unsigned iElectronMinus=0; iElectronMinus<electronsMinus.size() && nLept_<10; ++iElectronMinus ) {
       bool notSelected = ( electronsMinus[iElectronMinus]!=leptons[0] && electronsMinus[iElectronMinus]!=leptons[1] );
       if( notSelected ) {
         leptTypeLept_[nLept_] = 1;
         eLept_[nLept_] = electronsMinus[iElectronMinus].Energy();
         ptLept_[nLept_] = electronsMinus[iElectronMinus].Pt();
         etaLept_[nLept_] = electronsMinus[iElectronMinus].Eta();
         phiLept_[nLept_] = electronsMinus[iElectronMinus].Phi();
         chargeLept_[nLept_] = electronsMinus[iElectronMinus].charge;
         combinedIsoRelLept_[nLept_] = electronsMinus[iElectronMinus].isolation;
         if( event_==DEBUG_EVENTNUMBER ) 
           std::cout << "Added a negative electron: pt=" << ptLept_[nLept_] << ", eta=" << etaLept_[nLept_] << std::endl;
         nLept_++;
       }
     }

*/
/*
     // --------------------
     // match leptons to MC:
     // --------------------
     int correctIdMc = (leptType_==0 ) ? 13 : 11;

     for( unsigned iLept=0; iLept<leptons.size(); ++iLept ) {

       float deltaRmin = 100.;
       TLorentzVector matchedLeptonMC;

       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         if( statusMc[iMc]==1 && fabs(idMc[iMc])==correctIdMc && idMc[mothMc[mothMc[iMc]]]==23 ) {

           TLorentzVector* thisParticle = new TLorentzVector();
           thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
           float thisDeltaR = leptons[iLept].DeltaR( *thisParticle );
           if( thisDeltaR < deltaRmin ) {
             deltaRmin = thisDeltaR;
             matchedLeptonMC = *thisParticle;
           }

           delete thisParticle;
           thisParticle = 0;

         } //if correct id mc

       } // for i mc

       if( !noLeptons ) {
         if( leptType_==0 ) {
           h1_deltaRmatching_muons->Fill( deltaRmin );
           if( deltaRmin<0.1 ) {
             h1_passed_vs_ptMuon->Fill( matchedLeptonMC.Pt() );
           }
         } else if( leptType_==1 ) { 
           h1_deltaRmatching_electrons->Fill( deltaRmin );
           if( deltaRmin<0.1 ) {
             h1_passed_vs_ptEle->Fill( matchedLeptonMC.Pt() );
           }
         }  //if lept type
       } //if yes leptons


     } //for i leptons
*/


/*
     // ------------------
     // JETS
     // ------------------

     if( event_==DEBUG_EVENTNUMBER ) 
       std::cout << "Now adding in jets." << std::endl;

     float jetPt_thresh = 18.; //so that we can vary the JER uncert up one sigma (max value is ~10%)
     float jetEta_thresh = 2.5;

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

       thisJet.nCharged = chargedHadronMultiplicityAK5PFPUcorrJet[iJet]+electronMultiplicityAK5PFPUcorrJet[iJet]+muonMultiplicityAK5PFPUcorrJet[iJet];
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
         if( idMc[iPart]!=21 && abs(idMc[iPart])>5 ) continue; //quarks or gluons (excluding top)
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


     if( event_==DEBUG_EVENTNUMBER ) 
       std::cout << "leadJets.size(): " << leadJets.size() << std::endl;

     if( leadJets.size()<2 ) continue;
     if( leadJets[1].Pt()<jetPt_thresh ) continue; //at least 2 jets over thresh



     nJets_ = 0;
     nPart_ = 0;


     for( unsigned iJet=0; iJet<leadJets.size() && nJets_<50; ++iJet ) {

       if( leadJets[iJet].Pt()<jetPt_thresh ) continue;
       if( fabs(leadJets[iJet].Eta())>jetEta_thresh ) continue;

       eJet_[nJets_] = leadJets[iJet].Energy();
       ptJet_[nJets_] = leadJets[iJet].Pt();
       etaJet_[nJets_] = leadJets[iJet].Eta();
       phiJet_[nJets_] = leadJets[iJet].Phi();

       fJetCorrUnc->setJetPt(leadJets[iJet].Pt());   
       fJetCorrUnc->setJetEta(leadJets[iJet].Eta()); 
       ptUncertJet_[nJets_] = fJetCorrUnc->getUncertainty(true);

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
  
       if( event_==DEBUG_EVENTNUMBER ) 
         std::cout << "Adding jet: pt=" << ptJet_[nJets_] << " eta=" << etaJet_[nJets_] << std::endl;

       nJets_++;

     }

     if( event_==DEBUG_EVENTNUMBER ) 
       std::cout << "Found a total of " << nJets_ << " jets." << std::endl;

     // at least 3 jets in the event
     if( nJets_<3 ) continue;



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

*/

     } // if is mc



     if( event_==DEBUG_EVENTNUMBER ) 
       std::cout << "This is event passed all cuts." << std::endl;

     reducedTree_->Fill(); 


   } //for entries

} //finalize



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}




std::vector<AnalysisLepton> getBestZMassPair( const std::vector<AnalysisLepton>& leptPlus, const std::vector<AnalysisLepton>& leptMinus ) {

  std::vector<AnalysisLepton> returnLeptons;

  float bestMZ=999999999.;
  int iPlus_found = -1;
  int iMinus_found = -1;

  for( unsigned iPlus=0; iPlus<leptPlus.size(); ++iPlus ) {
    for( unsigned iMinus=0; iMinus<leptMinus.size(); ++iMinus ) {
      TLorentzVector l1( leptPlus[iPlus] );
      TLorentzVector l2( leptMinus[iMinus] );
      TLorentzVector dilepton = l1+l2;
      if( returnLeptons.size()==0 ) {
        returnLeptons.push_back(leptPlus [iPlus]);
        returnLeptons.push_back(leptMinus[iMinus]);
        iPlus_found = iPlus;
        iMinus_found = iMinus;
        bestMZ= dilepton.M();
      } else if( fabs(dilepton.M()-mZ) < fabs(bestMZ-mZ) ) { //already found a pair
        returnLeptons.clear();
        returnLeptons.push_back(leptPlus [iPlus]);
        returnLeptons.push_back(leptMinus[iMinus]);
        iPlus_found = iPlus;
        iMinus_found = iMinus;
        bestMZ= dilepton.M();
      }
    }  // for minus
  }  // for plus

  return returnLeptons;

}
