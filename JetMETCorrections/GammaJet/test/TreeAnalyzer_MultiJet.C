#include "TreeAnalyzer_MultiJet.h"


#include <iostream>
#include "TMath.h"
#include "AnalysisJet.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "QGLikelihood/QGLikelihoodCalculator.h"




TreeAnalyzer_MultiJet::TreeAnalyzer_MultiJet( const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags, TTree* tree ) :
     TreeAnalyzer( "MultiJet", dataset, recoType, jetAlgo, flags, tree ) {


} //constructor



void TreeAnalyzer_MultiJet::CreateOutputFile() {


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



  jetTree_->Branch("nJet",  &nJet_,  "nJet_/I");
  jetTree_->Branch("eJet",  eJet_,  "eJet_[nJet_]/F");
  jetTree_->Branch( "ptJet",  ptJet_,  "ptJet_[nJet_]/F");
  jetTree_->Branch( "ptRawJet",  ptRawJet_,  "ptRawJet_[nJet_]/F");
  jetTree_->Branch("etaJet", etaJet_, "etaJet_[nJet_]/F");
  jetTree_->Branch("phiJet", phiJet_, "phiJet_[nJet_]/F");
  jetTree_->Branch( "ptDJet",  ptDJet_,  "ptDJet_[nJet_]/F");
  jetTree_->Branch( "rmsCandJet",  rmsCandJet_,  "rmsCandJet_[nJet_]/F");
  jetTree_->Branch( "QGLikelihoodJet",  QGLikelihoodJet_,  "QGLikelihoodJet_[nJet_]/F");
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

  jetTree_->Branch("ht_akt5",&ht_akt5_,"ht_akt5_/F");

  jetTree_->Branch("passed_HT250", &passed_HT250_, "passed_HT250_/O");
  jetTree_->Branch("passed_HT300", &passed_HT300_, "passed_HT300_/O");
  jetTree_->Branch("passed_HT350", &passed_HT350_, "passed_HT350_/O");
  jetTree_->Branch("passed_HT400", &passed_HT400_, "passed_HT400_/O");
  jetTree_->Branch("passed_HT450", &passed_HT450_, "passed_HT450_/O");
  jetTree_->Branch("passed_HT500", &passed_HT500_, "passed_HT500_/O");
  jetTree_->Branch("passed_HT550", &passed_HT550_, "passed_HT550_/O");
  jetTree_->Branch("passed_HT600", &passed_HT600_, "passed_HT600_/O");
  jetTree_->Branch("passed_HT650", &passed_HT650_, "passed_HT650_/O");
  jetTree_->Branch("passed_HT700", &passed_HT700_, "passed_HT700_/O");
  jetTree_->Branch("passed_HT750", &passed_HT750_, "passed_HT750_/O");


} 



TreeAnalyzer_MultiJet::~TreeAnalyzer_MultiJet() {

  outfile_->cd();

}



void TreeAnalyzer_MultiJet::Loop()
{

   DEBUG_VERBOSE_ = false;

   if (fChain == 0) return;


   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();


   Long64_t nbytes = 0, nb = 0;

   TRandom3 rand;


   Int_t nJet_akt5;
   fChain->SetBranchAddress( "nJet_akt5", &nJet_akt5 );
   Float_t ptCorrJet_akt5[100];
   fChain->SetBranchAddress( "ptCorrJet_akt5 ", ptCorrJet_akt5 );
   Float_t etaJet_akt5[100];
   fChain->SetBranchAddress( "etaJet_akt5", etaJet_akt5 );



   QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_3_patch5/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");


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

     
     ht_akt5_ = 0.;

     for(unsigned int iCaloJet=0; iCaloJet<nJet_akt5 && iCaloJet<100; ++iCaloJet) {

       if( ptCorrJet_akt5[iCaloJet] > 40. && fabs(etaJet_akt5[iCaloJet])<3. ) ht_akt5_ += ptCorrJet_akt5[iCaloJet];

     }  // for calojets


     // will rely on HT600 trigger on data, so preselect:
     if( ht_akt5_ < 500. ) continue;

     


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

       thisJet.SetPtEtaPhiE( thisJet.ptCorrReco, thisJet.etaReco, thisJet.phiReco, thisJet.eReco*thisJet.ptCorrReco/thisJet.ptReco );

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


     } //for reco jets

     
     if( jets.size()<4 ) continue; //at least 4 jets


//   // will be relying mostly on HLT_HT600 trigger, so the following requirement is ~100% efficient:
//   if( jets[0]->Pt() + jets[1]->Pt() + jets[2]->Pt() + jets[3]->Pt() < 375. ) continue; // preselection


     for( unsigned iJet=0; iJet<jets.size(); iJet++ ) {

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

       if( fabs(jets[iJet]->Eta())<2.4 ) 
         QGLikelihoodJet_[iJet] = qglikeli->computeQGLikelihoodPU( jets[iJet]->Pt(), rhoPF, jets[iJet]->nCharged(), jets[iJet]->nNeutral(), jets[iJet]->ptD );
       else
         QGLikelihoodJet_[iJet] = 0.;


     }



     passed_HT250_ = passedTrigger_regexp("HLT_HT250_");
     passed_HT300_ = passedTrigger_regexp("HLT_HT300_");
     passed_HT350_ = passedTrigger_regexp("HLT_HT350_");
     passed_HT400_ = passedTrigger_regexp("HLT_HT400_");
     passed_HT450_ = passedTrigger_regexp("HLT_HT450_");
     passed_HT500_ = passedTrigger_regexp("HLT_HT500_");
     passed_HT550_ = passedTrigger_regexp("HLT_HT550_");
     passed_HT600_ = passedTrigger_regexp("HLT_HT600_");
     passed_HT650_ = passedTrigger_regexp("HLT_HT650_");
     passed_HT700_ = passedTrigger_regexp("HLT_HT700_");
     passed_HT750_ = passedTrigger_regexp("HLT_HT750_");


     jetTree_->Fill(); 

   } //for entries


 }

