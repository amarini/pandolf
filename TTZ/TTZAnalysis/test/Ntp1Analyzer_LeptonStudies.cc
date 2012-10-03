#include "Ntp1Analyzer_LeptonStudies.h"


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

#include "JetCorrectionUncertainty.h"


//#include "fitTools.h"


int DEBUG_EVENTNUMBER = 157480550;
float mZ = 91.1876;



float getWeightPU(Int_t nPU);






Ntp1Analyzer_LeptonStudies::Ntp1Analyzer_LeptonStudies( const std::string& dataset, const std::string& flags, TTree* tree ) :
     Ntp1Analyzer( "LeptonStudies", dataset, flags, tree ) {



} //constructor



void Ntp1Analyzer_LeptonStudies::CreateOutputFile() {

  Ntp1Analyzer::CreateOutputFile();

  
  reducedTree_->Branch("run",&run_,"run_/I");
  reducedTree_->Branch("LS",&LS_,"LS_/I");
  reducedTree_->Branch("event",&event_,"event_/I");
  reducedTree_->Branch("nPU",&nPU_,"nPU_/I");
  reducedTree_->Branch("nPU_ave",&nPU_ave_,"nPU_ave_/I");
  reducedTree_->Branch("nvertex",&nvertex_,"nvertex_/I");
  reducedTree_->Branch("rhoPF",&rhoPF_,"rhoPF_/F");
  reducedTree_->Branch("genWeight",&genWeight_,"genWeight_/F");
  reducedTree_->Branch("eventWeight",&eventWeight_,"eventWeight_/F");
  reducedTree_->Branch("eventWeightPU",&eventWeightPU_,"eventWeightPU_/F");
  reducedTree_->Branch("eventWeightPU_ave",&eventWeightPU_ave_,"eventWeightPU_ave_/F");

  reducedTree_->Branch("nEle",&nEle_,"nEle_/I");
  reducedTree_->Branch("ptEle",ptEle_,"ptEle_[nEle_]/F");
  reducedTree_->Branch("etaEle",etaEle_,"etaEle_[nEle_]/F");
  reducedTree_->Branch("pfIsoEle",pfIsoEle_,"pfIsoEle_[nEle_]/F");
  reducedTree_->Branch("isLooseSUSYElectronEle",isLooseSUSYElectronEle_,"isLooseSUSYElectronEle_[nEle_]/O");
  reducedTree_->Branch("isTightSUSYElectronEle",isTightSUSYElectronEle_,"isTightSUSYElectronEle_[nEle_]/O");
  reducedTree_->Branch("isNotConversionEle",isNotConversionEle_,"isNotConversionEle_[nEle_]/O");
  reducedTree_->Branch("passedHLTEle",passedHLTEle_,"passedHLTEle_[nEle_]/O");
  reducedTree_->Branch("mvaidtrigEle",mvaidtrigEle_,"mvaidtrigEle_[nEle_]/F");
  reducedTree_->Branch("matchedToGenEle",matchedToGenEle_,"matchedToGenEle_[nEle_]/O");
  

} 



Ntp1Analyzer_LeptonStudies::~Ntp1Analyzer_LeptonStudies() {

  outfile_->cd();

}



void Ntp1Analyzer_LeptonStudies::Loop()
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
   //PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
   //PUWeight* fPUWeight_ave = new PUWeight(-1, "2011A", puType_ave);
   ////PUWeight* fPUWeight = new PUWeight(1089.2, "2011A", puType);
   //TFile* filePU = TFile::Open("Pileup_2011_to_173692_LPLumiScale_68mb.root");
   //TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
   //fPUWeight->SetDataHistogram(h1_nPU_data);
   //fPUWeight_ave->SetDataHistogram(h1_nPU_data);


   
   // this file is obtained with the instructions found in: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#GetTxtFiles
   //JetCorrectionUncertainty *fJetCorrUnc = new JetCorrectionUncertainty("AK5PF_Uncertainty_GR_R_42_V19.txt");
   JetCorrectionUncertainty *fJetCorrUnc = new JetCorrectionUncertainty("GR_R_52_V9_Uncertainty_AK5PF.txt");



   float nCounterPU=0.;
   float nCounterPU_ave=0.;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%10000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;



     run_ = runNumber;
     LS_ = lumiBlock;
     event_ = eventNumber;
     genWeight_ = genWeight; //default
     eventWeight_ = -1.; //default

     if( dataset_tstr.Contains("spadhi") )
       nPU_ = nPU[0]; //generated only with one nPU
     else
       nPU_ = nPU[1]; //in time PU only




     if( !isGoodEvent(jentry) ) continue; //this takes care also of trigger


     if( nPV==0 ) continue;
     bool goodVertex = (ndofPV[0] >= 4.0 && sqrt(PVxPV[0]*PVxPV[0]+PVyPV[0]*PVyPV[0]) < 2. && fabs(PVzPV[0]) < 24. );
     //if( !goodVertex ) continue;
  
     nPU_ave_ = 0.;
     for( unsigned iBX=0; iBX<nBX; ++iBX ) {
       nPU_ave_ += nPU[iBX]; 
     }
     nPU_ave_ /= (float)nBX;

     // PU reweighting:
     eventWeightPU_=1.;
     eventWeightPU_ave_=1.;
     //if( isMC_ ) {
     //  eventWeightPU_ = fPUWeight->GetWeight(nPU_);
     //  eventWeightPU_ave_ = fPUWeight_ave->GetWeight(nPU_ave_);
     //}
     nCounterPU += eventWeightPU_;
     nCounterPU_ave += eventWeightPU_ave_;




     nvertex_ = nPV;
     rhoPF_ = rhoFastjet;


     // save trigger info:


     //bool isMC = ( runNumber < 5 );


     ptHat_ = (isMC_) ? genPtHat : ptHat_;



     bool noLeptons = false;
     TLorentzVector lept1MC, lept2MC;
     int zIndexqq=-1;
     int zIndexll=-1;

/*
     if( isMC_ ) {


       // first look for Z->qq
       std::vector<TLorentzVector> quarksMC;

       for( unsigned iMc=0; iMc<nMc && quarksMC.size()<2; ++iMc ) {

         // quarks have status 3
         if( statusMc[iMc] != 3 ) continue;

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );
         if( thisParticle->Pt()<0.1 ) continue;

         if( fabs(idMc[iMc])<7 && idMc[mothMc[iMc]]==23 ) {
           zIndexqq = mothMc[iMc];
           quarksMC.push_back( *thisParticle );
         }

       }

       // (checked that always 2 quarks are found)
       if( quarksMC.size()==2 && zIndexqq!=-1 ) {

         TLorentzVector ZqqMC;
         ZqqMC.SetPtEtaPhiE( pMc[zIndexqq]*sin(thetaMc[zIndexqq]), etaMc[zIndexqq], phiMc[zIndexqq], energyMc[zIndexqq] );

         ptZqqMC_  = ZqqMC.Pt();
         eZqqMC_   = ZqqMC.Energy();
         etaZqqMC_ = ZqqMC.Eta();
         phiZqqMC_ = ZqqMC.Phi();

      // float ptZqq = pMc[zIndexqq]*sin(thetaMc[zIndexqq]);
      // h1_ptHadronicZ->Fill( ptZqq );

      // float deltaRqq = quarksMC[0].DeltaR(quarksMC[1]);
      // h1_deltaRqq->Fill(deltaRqq);

       }

       // now look for Z->ll

       std::vector<TLorentzVector> electronsMC;
       std::vector<TLorentzVector> muonsMC;

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

         delete thisParticle;
         thisParticle = 0;

       }

       if( electronsMC.size()==2 ) {
         if( electronsMC[0].Pt() > electronsMC[1].Pt() ) {
           lept1MC = electronsMC[0];
           lept2MC = electronsMC[1];
         } else {
           lept1MC = electronsMC[1];
           lept2MC = electronsMC[0];
         }
         if( (fabs(lept1MC.Eta()) < 2.5) && ( fabs(lept1MC.Eta())<1.4442 || fabs(lept1MC.Eta())>1.566) ) h1_nEvents_vs_ptEle->Fill( lept1MC.Pt() );
         if( (fabs(lept2MC.Eta()) < 2.5) && ( fabs(lept2MC.Eta())<1.4442 || fabs(lept2MC.Eta())>1.566) ) h1_nEvents_vs_ptEle->Fill( lept2MC.Pt() );
       } else if( muonsMC.size()==2 ) {
         if( muonsMC[0].Pt() > muonsMC[1].Pt() ) {
           lept1MC = muonsMC[0];
           lept2MC = muonsMC[1];
         } else {
           lept1MC = muonsMC[1];
           lept2MC = muonsMC[0];
         }
         if( fabs(lept1MC.Eta()) < 2.4 ) h1_nEvents_vs_ptMuon->Fill( lept1MC.Pt() );
         if( fabs(lept1MC.Eta()) < 2.1 || fabs(lept2MC.Eta()) < 2.1 ) h1_nEvents_vs_ptMuon->Fill( lept2MC.Pt() );
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

         if( muonsMC.size() > 0 ) leptTypeMC_ = 0;
         else if( electronsMC.size() > 0 ) leptTypeMC_ = 1;

       }


       // now look for the higgs:
       if( zIndexll!=-1 && zIndexqq!=-1 ) {

         int higgsIndex = mothMc[zIndexll];

         if( idMc[higgsIndex] == 25 ) {

           TLorentzVector HiggsMC;
           HiggsMC.SetPtEtaPhiE( pMc[higgsIndex]*sin(thetaMc[higgsIndex]), etaMc[higgsIndex], phiMc[higgsIndex], energyMc[higgsIndex] );

           eHiggsMC_   = HiggsMC.Energy(); 
           ptHiggsMC_  = HiggsMC.Pt(); 
           etaHiggsMC_ = HiggsMC.Eta(); 
           phiHiggsMC_ = HiggsMC.Phi(); 

         } // if higgs

       } //if found two Z's

     } //if isMC
*/





     // -----------------------------
     //      FROM NOW ON RECO
     // -----------------------------




     if( event_==DEBUG_EVENTNUMBER ) {
       std::cout << std::endl << std::endl;
       std::cout << "----- LOG for run: " << run_ << "    event: " << event_ << std::endl;
       std::cout << std::endl << "*** Muons:" << std::endl;
     }

     // ------------------
     // MUONS
     // ------------------

     std::vector<AnalysisLepton> muons;


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

       thisMuon.isGlobalMuon = (muonIdMuon[iMuon]>>13)&1;
       thisMuon.isGlobalMuonPromptTight = (muonIdMuon[iMuon]>>8)&1;
       thisMuon.isAllTrackerMuon = (muonIdMuon[iMuon]>>11)&1;
       thisMuon.isPFMuon = pfmuonIdMuon[iMuon];

       thisMuon.pixelHits = numberOfValidPixelBarrelHitsTrack[trackIndexMuon[iMuon]]+numberOfValidPixelEndcapHitsTrack[trackIndexMuon[iMuon]];
       thisMuon.trackerHits = trackValidHitsTrack[trackIndexMuon[iMuon]];

       thisMuon.nMatchedStations = numberOfMatchesMuon[iMuon];

       int globalMuonTrack = combinedTrackIndexMuon[iMuon];
       thisMuon.normChiSquare = (thisMuon.isGlobalMuon) ? trackNormalizedChi2GlobalMuonTrack[globalMuonTrack] : -1;
       thisMuon.nValidMuonHits = (thisMuon.isGlobalMuon) ? numberOfValidMuonHitsGlobalMuonTrack[globalMuonTrack] : -1;


       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisMuon.isGlobalMuonPromptTight: " << thisMuon.isGlobalMuonPromptTight << std::endl;
         std::cout << "thisMuon.isAllTrackerMuon: " << thisMuon.isAllTrackerMuon << std::endl;
         std::cout << "thisMuon.pixelHits: " << thisMuon.pixelHits << std::endl;
         std::cout << "thisMuon.trackerHits: " << thisMuon.trackerHits << std::endl;
         std::cout << "thisMuon.nMatchedStations: " << thisMuon.nMatchedStations << std::endl;
       }


       int ctfMuon = trackIndexMuon[iMuon]; 
       thisMuon.dxy = transvImpactParTrack[ctfMuon];
       thisMuon.dz = muonDzPV(iMuon,0);

       thisMuon.sumPt03 = sumPt03Muon[iMuon];
       thisMuon.emEt03  = emEt03Muon[iMuon];
       thisMuon.hadEt03 = hadEt03Muon[iMuon];
       thisMuon.isolation = thisMuon.combinedIsoRel();
       thisMuon.mvaisoMuon = mvaisoMuon[iMuon];

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "thisMuon.dxy: " << thisMuon.dxy << std::endl;
         std::cout << "thisMuon.dz: " << thisMuon.dz << std::endl;
         std::cout << "thisMuon.sumPt03: " << thisMuon.sumPt03 << std::endl;
         std::cout << "thisMuon.emEt03: " << thisMuon.emEt03 << std::endl;
         std::cout << "thisMuon.hadEt03: " << thisMuon.hadEt03 << std::endl;
       }

       //if( !thisMuon.passedVBTF() ) continue;
       if( !thisMuon.isTightMuon2012() ) continue;

       if( event_==DEBUG_EVENTNUMBER ) {
         std::cout << "PASSED VBTF. ";
         if( thisMuon.charge > 0 ) std::cout << "Adding to collection of positive muons." << std::endl;
         else std::cout << "Adding to collection of negative muons." << std::endl;
       }

       muons.push_back(thisMuon);

        
       // match to gen:
       float deltaRmin_muon = 999.;
       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         // partons only
         if( statusMc[iMc] != 3 ) continue;
         if( !(abs(idMc[iMc])==13) ) continue;
         if( !(idMc[mothMc[iMc]]==23 || abs(idMc[mothMc[iMc]]==24) ) ) continue; //muons from W/Z only

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       
         float thisDeltaR =  thisParticle->DeltaR(thisMuon);
         if( thisDeltaR < deltaRmin_muon ) {
           deltaRmin_muon = thisDeltaR;
         }

         delete thisParticle;
         thisParticle = 0;

       }

       thisMuon.matchedToGen = (deltaRmin_muon<0.1); 

     } //for muons





     // ------------------
     // ELECTRONS
     // ------------------

     std::vector<AnalysisLepton> electrons;

     if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ )
       std::cout << std::endl << "*** Electrons:" << std::endl;


     nEle_ = 0;

     //for( unsigned int iEle=0; (iEle<nEle) && (electrons.size()<2); ++iEle ) {
     for( unsigned int iEle=0; (iEle<nEle); ++iEle ) {

       if( nEle_ >= 10 ) continue;

       AnalysisElectron thisEle( pxEle[iEle], pyEle[iEle], pzEle[iEle], energyEle[iEle] );
       thisEle.charge = chargeEle[iEle];

       float scEta = (superClusterIndexEle[iEle]>=0) ? etaSC[superClusterIndexEle[iEle]] : etaPFSC[PFsuperClusterIndexEle[iEle]];

       if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) {
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

       
       int gsf = gsfTrackIndexEle[iEle];
       thisEle.dxy = transvImpactParGsfTrack[gsf];
       thisEle.dz = eleDzPV(iEle,0);


       // triple charge consistency:
       int scCharge = scPixChargeEle[iEle];
       int ctfCharge = chargeTrack[trackIndexEle[iEle]];
       int gsfCharge = chargeGsfTrack[gsfTrackIndexEle[iEle]];
       thisEle.tripleChargeConsist = !((scCharge!=ctfCharge) || (scCharge!=gsfCharge) || (gsfCharge!=ctfCharge));


       // isolation
       thisEle.dr03TkSumPt = dr03TkSumPtEle[iEle];
       thisEle.dr03EcalRecHitSumEt = dr03EcalRecHitSumEtEle[iEle];
       thisEle.dr03HcalTowerSumEt = dr03HcalTowerSumEtEle[iEle];
       thisEle.isolation = thisEle.combinedIsoRel();
       thisEle.pfCandChargedIso04 = pfCandChargedIso04Ele[iEle];
       thisEle.pfCandNeutralIso04 = pfCandNeutralIso04Ele[iEle];
       thisEle.pfCandPhotonIso04 = pfCandPhotonIso04Ele[iEle];
       thisEle.rhoJetsFastJet = rhoJetsFastjet;

       // electron ID
       thisEle.sigmaIetaIeta = (superClusterIndexEle[iEle]>=0) ? sqrt(covIEtaIEtaSC[superClusterIndexEle[iEle]]) : sqrt(covIEtaIEtaPFSC[PFsuperClusterIndexEle[iEle]]);
       thisEle.deltaPhiAtVtx = deltaPhiAtVtxEle[iEle];
       thisEle.deltaEtaAtVtx = deltaEtaAtVtxEle[iEle];
       thisEle.hOverE = hOverEEle[iEle];
       thisEle.eOverP = eSuperClusterOverPEle[iEle];
       thisEle.pAtVertex = sqrt( pxGsfTrack[gsf]*pxGsfTrack[gsf] + pyGsfTrack[gsf]*pyGsfTrack[gsf] + pzGsfTrack[gsf]*pzGsfTrack[gsf] );

       // conversion rejection
       thisEle.expInnerLayersGsfTrack = expInnerLayersGsfTrack[gsf];
       thisEle.convDist = convDistEle[iEle];
       thisEle.convDcot = convDcotEle[iEle];
       thisEle.hasMatchedConversion = hasMatchedConversionEle[iEle];

       // electron ID mva:
       thisEle.mvaidtrigEle = mvaidtrigEle[iEle];


       if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) {
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




       //bool passed = thisEle.isGoodElectron2012_CutsLoose();
       //bool passed = thisEle.passedHLT2012();
       //if( !passed ) continue;



       if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) std::cout << "Is good electron (2012)." << std::endl;

//     // additional ID to be as tight as trigger (HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL):
//     if( fabs(scEta)<1.4442 ) { //barrel
//       if( fabs(thisEle.deltaPhiAtVtx) > 0.15 ) continue;
//       if( thisEle.hOverE > 0.12 ) continue; //conforming to SSDL analysis
//     } else { //endcaps
//       if( fabs(thisEle.deltaPhiAtVtx) > 0.1 ) continue;
//       if( thisEle.hOverE > 0.15 ) continue; // looks like (from ntuples) there's a cut at 0.15 at hlt, not 0.075 as CaloIdT says
//     }

       if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) std::cout << "Passed additional eleID cuts (HLT)." << std::endl;


       // check that not matched to muon (clean electrons faked by muon MIP):
       bool matchedtomuon=false;
       for( std::vector<AnalysisLepton>::iterator iMu=muons.begin(); iMu!=muons.end(); ++iMu )
         if( iMu->DeltaR(thisEle)<0.1 ) matchedtomuon=true;

       if( matchedtomuon ) continue;

       if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) std::cout << "Not matched to any muon." << std::endl;

       if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) {
         if( thisEle.charge > 0 ) std::cout << "Adding to collection of positive electrons." << std::endl;
         else std::cout << "Adding to collection of negative electrons." << std::endl;
       }

       electrons.push_back(thisEle);
       
       // match to gen:
       float deltaRmin_Ele = 999.;
       for( unsigned iMc=0; iMc<nMc; ++iMc ) {

         // partons only
         if( statusMc[iMc] != 3 ) continue;
         if( !(abs(idMc[iMc])==11) ) continue;
         if( !(idMc[mothMc[iMc]]==23 || abs(idMc[mothMc[iMc]]==24) ) ) continue; //Eles from W/Z only

         TLorentzVector* thisParticle = new TLorentzVector();
         thisParticle->SetPtEtaPhiE( pMc[iMc]*sin(thetaMc[iMc]), etaMc[iMc], phiMc[iMc], energyMc[iMc] );

       
         float thisDeltaR =  thisParticle->DeltaR(thisEle);
         if( thisDeltaR < deltaRmin_Ele ) {
           deltaRmin_Ele = thisDeltaR;
         }

         delete thisParticle;
         thisParticle = 0;

       }

       thisEle.matchedToGen = (deltaRmin_Ele<0.1); 

       ptEle_[nEle_] = thisEle.Pt();
       etaEle_[nEle_] = thisEle.Eta();
       pfIsoEle_[nEle_] = thisEle.getPfIso_rhoCorrected();
       isLooseSUSYElectronEle_[nEle_] = thisEle.electronID2012_SUSYloose();
       isTightSUSYElectronEle_[nEle_] = thisEle.electronID2012_SUSYtight();
       isNotConversionEle_[nEle_] = thisEle.conversionRejection2012_CutsLoose();
       passedHLTEle_[nEle_] = thisEle.passedHLT2012();
       mvaidtrigEle_[nEle_] = thisEle.mvaidtrigEle;
       matchedToGenEle_[nEle_] = thisEle.matchedToGen;

       nEle_++;

     } //for electrons



     if( event_==DEBUG_EVENTNUMBER || DEBUG_VERBOSE_ ) 
       std::cout << "This is event passed all cuts." << std::endl;

     reducedTree_->Fill(); 


   } //for entries


} //loop



double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}




