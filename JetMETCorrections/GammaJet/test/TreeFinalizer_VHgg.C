
#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "AnalysisJet.C"
#include "AnalysisPhoton.C"

#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/emanuele/CommonTools/src/PUWeight.C"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TFitParticleEtThetaPhi.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/TKinFitter.cc"
#include "/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/KinematicFit/DiJetKinFitter.cc"


bool isAOD_ = true;
bool DEBUG_ = false;
bool MCassoc_ = false;
std::string RECOTYPE_;
std::string ALGOTYPE_;
std::string PARTTYPE_;

TChain* tree;
Double_t totalLumi=0.;


void addInput(const std::string& dataset);
void finalize(const std::string& dataset);


void finalize_all() {

  finalize( "WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1" );
  finalize( "WH_ZH_HToGG_M-125_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1" );
  finalize( "DiPhotonBoxBorn" );

}




void finalize(const std::string& dataset) {


  std::string outfileName;
  if( DEBUG_ ) outfileName = "provaVHgg_"+dataset;
  else {
   if(dataset!="") outfileName = "VHgg_"+dataset;
   else outfileName = "VHgg";
  }
  if( MCassoc_ ) outfileName = outfileName + "_MCassoc";
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  tree = new TChain("jetTree");


  if( dataset=="DiPhotonBoxBorn" ) {
    addInput( "DiPhotonJets_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1" );
    addInput( "DiPhotonBox_Pt-25To250_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1");
  } else {
    addInput( dataset );
  }






  TH1F* h1_nvertex = new TH1F("nvertex", "", 26, -0.5, 25.5);
  h1_nvertex->Sumw2();
  TH1F* h1_nvertexPU = new TH1F("nvertexPU", "", 26, -0.5, 25.5);
  h1_nvertexPU->Sumw2();


  TH1D* h1_ptPhot1 = new TH1D("ptPhot1", "", 500, 0., 500.);
  h1_ptPhot1->Sumw2();
  TH1D* h1_etaPhot1 = new TH1D("etaPhot1", "", 200, -3., 3.);
  h1_etaPhot1->Sumw2();

  TH1D* h1_ptPhot2 = new TH1D("ptPhot2", "", 500, 0., 500.);
  h1_ptPhot2->Sumw2();
  TH1D* h1_etaPhot2 = new TH1D("etaPhot2", "", 200, -3., 3.);
  h1_etaPhot2->Sumw2();

  TH1D* h1_mgg = new TH1D("mgg", "", 1000, 0., 1000.);
  h1_mgg->Sumw2();

  TH1D* h1_kinfit_chiSquare = new TH1D("kinfit_chiSquare", "", 100, 0., 10.);
  h1_kinfit_chiSquare->Sumw2();
  TH1D* h1_kinfit_chiSquareProb = new TH1D("kinfit_chiSquareProb", "", 100, 0., 1.0001);
  h1_kinfit_chiSquareProb->Sumw2();


  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 500, 0., 500.);
  h1_ptJet1->Sumw2();
  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 500, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_pdgIdJet1 = new TH1D("pdgIdJet1", "", 38, -15.5, 22.5);
  h1_pdgIdJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();

  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 500, 0., 500.);
  h1_ptJet2->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 500, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_pdgIdJet2 = new TH1D("pdgIdJet2", "", 38, -15.5, 22.5);
  h1_pdgIdJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet2->Sumw2();

  TH1D* h1_mjj = new TH1D("mjj", "", 1000, 0., 1000.);
  h1_mjj->Sumw2();
  TH1D* h1_QGLikelihoodProd = new TH1D("QGLikelihoodProd", "", 100, 0., 1.0001);
  h1_QGLikelihoodProd->Sumw2();


  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree->SetBranchAddress("rhoPF", &rhoPF);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Float_t eMet;
  tree->SetBranchAddress("epfMet", &eMet);
  Float_t phiMet;
  tree->SetBranchAddress("phipfMet", &phiMet);


  Int_t nPU;
  tree->SetBranchAddress("nPU", &nPU);
  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t ePhot1;
  tree->SetBranchAddress("ePhot1", &ePhot1);
  Float_t ptPhot1;
  tree->SetBranchAddress("ptPhot1", &ptPhot1);
  Float_t etaPhot1;
  tree->SetBranchAddress("etaPhot1", &etaPhot1);
  Float_t phiPhot1;
  tree->SetBranchAddress("phiPhot1", &phiPhot1);

  Float_t ePhotGen1;
  tree->SetBranchAddress("ePhotGen1", &ePhotGen1);
  Float_t ptPhotGen1;
  tree->SetBranchAddress("ptPhotGen1", &ptPhotGen1);
  Float_t etaPhotGen1;
  tree->SetBranchAddress("etaPhotGen1", &etaPhotGen1);
  Float_t phiPhotGen1;
  tree->SetBranchAddress("phiPhotGen1", &phiPhotGen1);

  Float_t hcalIsoPhot1;
  tree->SetBranchAddress("hcalIsoPhot1", &hcalIsoPhot1);
  Float_t ecalIsoPhot1;
  tree->SetBranchAddress("ecalIsoPhot1", &ecalIsoPhot1);
  Int_t nTrkIsoPhot1;
  tree->SetBranchAddress("nTrkIsoPhot1", &nTrkIsoPhot1);
  Float_t ptTrkIsoPhot1;
  tree->SetBranchAddress("ptTrkIsoPhot1", &ptTrkIsoPhot1);
  Float_t clusterMajPhot1;
  tree->SetBranchAddress("clusterMajPhot1", &clusterMajPhot1);
  Float_t clusterMinPhot1;
  tree->SetBranchAddress("clusterMinPhot1", &clusterMinPhot1);
  Int_t hasPixelSeedPhot1;
  tree->SetBranchAddress("hasPixelSeedPhot1", &hasPixelSeedPhot1);
  Float_t pid_twrHCALPhot1;
  tree->SetBranchAddress("pid_twrHCALPhot1", &pid_twrHCALPhot1);
  Float_t pid_HoverEPhot1;
  tree->SetBranchAddress("pid_HoverEPhot1", &pid_HoverEPhot1);
  Float_t pid_jurECALPhot1;
  tree->SetBranchAddress("pid_jurECALPhot1", &pid_jurECALPhot1);

  Float_t ePhot2;
  tree->SetBranchAddress("ePhot2", &ePhot2);
  Float_t ptPhot2;
  tree->SetBranchAddress("ptPhot2", &ptPhot2);
  Float_t etaPhot2;
  tree->SetBranchAddress("etaPhot2", &etaPhot2);
  Float_t phiPhot2;
  tree->SetBranchAddress("phiPhot2", &phiPhot2);

  Float_t ePhotGen2;
  tree->SetBranchAddress("ePhotGen2", &ePhotGen2);
  Float_t ptPhotGen2;
  tree->SetBranchAddress("ptPhotGen2", &ptPhotGen2);
  Float_t etaPhotGen2;
  tree->SetBranchAddress("etaPhotGen2", &etaPhotGen2);
  Float_t phiPhotGen2;
  tree->SetBranchAddress("phiPhotGen2", &phiPhotGen2);

  Float_t hcalIsoPhot2;
  tree->SetBranchAddress("hcalIsoPhot2", &hcalIsoPhot2);
  Float_t ecalIsoPhot2;
  tree->SetBranchAddress("ecalIsoPhot2", &ecalIsoPhot2);
  Int_t nTrkIsoPhot2;
  tree->SetBranchAddress("nTrkIsoPhot2", &nTrkIsoPhot2);
  Float_t ptTrkIsoPhot2;
  tree->SetBranchAddress("ptTrkIsoPhot2", &ptTrkIsoPhot2);
  Float_t clusterMajPhot2;
  tree->SetBranchAddress("clusterMajPhot2", &clusterMajPhot2);
  Float_t clusterMinPhot2;
  tree->SetBranchAddress("clusterMinPhot2", &clusterMinPhot2);
  Int_t hasPixelSeedPhot2;
  tree->SetBranchAddress("hasPixelSeedPhot2", &hasPixelSeedPhot2);
  Float_t pid_twrHCALPhot2;
  tree->SetBranchAddress("pid_twrHCALPhot2", &pid_twrHCALPhot2);
  Float_t pid_HoverEPhot2;
  tree->SetBranchAddress("pid_HoverEPhot2", &pid_HoverEPhot2);
  Float_t pid_jurECALPhot2;
  tree->SetBranchAddress("pid_jurECALPhot2", &pid_jurECALPhot2);



  Int_t nJet;
  tree->SetBranchAddress("nJet", &nJet);
  Float_t eJet[20];
  tree->SetBranchAddress("eJet", eJet);
  Float_t ptJet[20];
  tree->SetBranchAddress("ptJet", ptJet);
  Float_t ptRawJet[20];
  tree->SetBranchAddress("ptRawJet", ptRawJet);
  Float_t etaJet[20];
  tree->SetBranchAddress("etaJet", etaJet);
  Float_t phiJet[20];
  tree->SetBranchAddress("phiJet", phiJet);
  Float_t eChargedHadronsJet[20];
  tree->SetBranchAddress("eChargedHadronsJet", eChargedHadronsJet);
  Float_t ePhotonsJet[20];
  tree->SetBranchAddress("ePhotonsJet", ePhotonsJet);
  Int_t eNeutralHadronsJet[20];
  tree->SetBranchAddress("eNeutralHadronsJet", eNeutralHadronsJet);
  Int_t nChargedHadronsJet[20];
  tree->SetBranchAddress("nChargedHadronsJet", nChargedHadronsJet);
  Int_t nPhotonsJet[20];
  tree->SetBranchAddress("nPhotonsJet", nPhotonsJet);
  Int_t nNeutralHadronsJet[20];
  tree->SetBranchAddress("nNeutralHadronsJet", nNeutralHadronsJet);
  Float_t ptDJet[20];
  tree->SetBranchAddress("ptDJet", ptDJet);

  Float_t ptPartJet[20];
  tree->SetBranchAddress("ptPartJet", ptPartJet);
  Float_t etaPartJet[20];
  tree->SetBranchAddress("etaPartJet", etaPartJet);
  Float_t phiPartJet[20];
  tree->SetBranchAddress("phiPartJet", phiPartJet);
  Int_t pdgIdPartJet[20];
  tree->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);
  Int_t pdgIdMomJet[20];
  tree->SetBranchAddress("pdgIdMomJet", pdgIdMomJet);


  Bool_t passed_Photon10;
  tree->SetBranchAddress("passed_Photon10", &passed_Photon10);
  Bool_t passed_Photon15;
  tree->SetBranchAddress("passed_Photon15", &passed_Photon15);
  Bool_t passed_Photon20;
  tree->SetBranchAddress("passed_Photon20", &passed_Photon20);
  Bool_t passed_Photon25;
  tree->SetBranchAddress("passed_Photon25", &passed_Photon25);
  Bool_t passed_Photon30;
  tree->SetBranchAddress("passed_Photon30", &passed_Photon30);
  Bool_t passed_Photon35;
  tree->SetBranchAddress("passed_Photon35", &passed_Photon35);
//Bool_t passed_Photon40;
//tree->SetBranchAddress("passed_Photon40", &passed_Photon40);
  Bool_t passed_Photon50;
  tree->SetBranchAddress("passed_Photon50", &passed_Photon50);
//Bool_t passed_Photon60;
//tree->SetBranchAddress("passed_Photon60", &passed_Photon60);
  Bool_t passed_Photon70;
  tree->SetBranchAddress("passed_Photon70", &passed_Photon70);
  Bool_t passed_Photon75;
  tree->SetBranchAddress("passed_Photon75", &passed_Photon75);

  float ptJet1_t, ptJet2_t;
  float etaJet1_t, etaJet2_t;
  float QGLikelihoodJet1_t, QGLikelihoodJet2_t;
  float mgg, mjj;
  float chiSquareProb;



  TTree* tree_passedEvents = new TTree("tree_passedEvents", "");
  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "ptPhot1", &ptPhot1, "ptPhot1/F" );
  tree_passedEvents->Branch( "ptPhot2", &ptPhot2, "ptPhot2/F" );
  tree_passedEvents->Branch( "etaPhot1", &etaPhot1, "etaPhot1/F" );
  tree_passedEvents->Branch( "etaPhot2", &etaPhot2, "etaPhot2/F" );
  tree_passedEvents->Branch( "ptJet1", &ptJet1_t, "ptJet1_t/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2_t, "ptJet2_t/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1_t, "etaJet1_t/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2_t, "etaJet2_t/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet1", &QGLikelihoodJet1_t, "QGLikelihoodJet1_t/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet2", &QGLikelihoodJet2_t, "QGLikelihoodJet2_t/F" );
  tree_passedEvents->Branch( "mjj", &mjj, "mjj/F" );
  tree_passedEvents->Branch( "mgg", &mgg, "mgg/F" );
  tree_passedEvents->Branch( "chiSquareProb", &chiSquareProb, "chiSquareProb/F" );


  gROOT->cd();

  PUWeight* fPUWeight = new PUWeight();

  QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator("/cmsrm/pc18/pandolf/CMSSW_4_2_3_patch1/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");

  DiJetKinFitter* fitter_jets = new DiJetKinFitter( "fitter_jets", "fitter_jets", 85. );




  float selectedCorrectPair = 0.;
  float allPairs = 0.;

  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;

    h1_nvertex->Fill( nvertex, eventWeight);

    bool isMC = run<5;
    if( isMC ) {

      // PU reweighting:
     eventWeight *= fPUWeight->GetWeight(nPU);
    }

    h1_nvertexPU->Fill( nvertex, eventWeight);


    // kinematic selection on GEN photons
    if( ptPhotGen1 < 30. ) continue;
    if( ptPhotGen2 < 30. ) continue;
    if( fabs(etaPhotGen1) > 2.5 ) continue;
    if( fabs(etaPhotGen2) > 2.5 ) continue;

    AnalysisPhoton phot1, phot2;
    if( ptPhot1>ptPhot2 ) {
      phot1.SetPtEtaPhiE( ptPhot1, etaPhot1, phiPhot1, ePhot1 );
      phot2.SetPtEtaPhiE( ptPhot2, etaPhot2, phiPhot2, ePhot2 );
    } else {
      phot2.SetPtEtaPhiE( ptPhot1, etaPhot1, phiPhot1, ePhot1 );
      phot1.SetPtEtaPhiE( ptPhot2, etaPhot2, phiPhot2, ePhot2 );
    }


    h1_ptPhot1->Fill( phot1.Pt(), eventWeight );
    h1_etaPhot1->Fill( phot1.Eta(), eventWeight );
    
    h1_ptPhot2->Fill( phot2.Pt(), eventWeight );
    h1_etaPhot2->Fill( phot2.Eta(), eventWeight );

    TLorentzVector diPhoton = phot1 + phot2;
    h1_mgg->Fill( diPhoton.M(), eventWeight );
    
    float bestMjj = -1.;
    AnalysisJet jet1, jet2;
    bool correctPair=false;

//  // first look for matched jets:
//  for( unsigned iJet=0; iJet<nJet; ++iJet ) {

//    for( unsigned jJet=iJet+1; jJet<nJet; ++jJet ) {

//      AnalysisJet firstJet, secondJet;
//    
//      firstJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], eJet[iJet] );
//      secondJet.SetPtEtaPhiE( ptJet[jJet], etaJet[jJet], phiJet[jJet], eJet[jJet] );

//      firstJet.ptD = ptDJet[iJet];
//      firstJet.pdgIdPart = pdgIdPartJet[iJet];
//      firstJet.pdgIdMom = pdgIdMomJet[iJet];

//      secondJet.ptD = ptDJet[jJet];
//      secondJet.pdgIdPart = pdgIdPartJet[jJet];
//      secondJet.pdgIdMom = pdgIdMomJet[jJet];

//      if( fabs(firstJet.pdgIdMom)==23 || fabs(firstJet.pdgIdMom)==24 ) {

//        correctPair=true;
//        jet1 = firstJet;
//        jet2 = secondJet;

//        break;

//      } //if found correct pair

//    } //for j jets

//  } //for i jets


//  if( !correctPair ) {

    for( unsigned iJet=0; iJet<nJet; ++iJet ) {

      for( unsigned jJet=iJet+1; jJet<nJet; ++jJet ) {

        AnalysisJet firstJet, secondJet;
      
        firstJet.SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], ptJet[iJet]/ptRawJet[iJet]*eJet[iJet] );
        secondJet.SetPtEtaPhiE( ptJet[jJet], etaJet[jJet], phiJet[jJet], ptJet[jJet]/ptRawJet[jJet]*eJet[jJet] );

        if( firstJet.Pt()<20. ) continue;
        if( secondJet.Pt()<20. ) continue;
        if( fabs(firstJet.Eta())>2.4 ) continue;
        if( fabs(secondJet.Eta())>2.4 ) continue;

        firstJet.ptD = ptDJet[iJet];
        firstJet.nTracksReco = nChargedHadronsJet[iJet];
        firstJet.nNeutralHadronsReco = nNeutralHadronsJet[iJet];
        firstJet.nPhotonsReco = nPhotonsJet[iJet];
        firstJet.pdgIdPart = pdgIdPartJet[iJet];
        firstJet.pdgIdMom = pdgIdMomJet[iJet];

        secondJet.ptD = ptDJet[jJet];
        secondJet.nTracksReco = nChargedHadronsJet[jJet];
        secondJet.nNeutralHadronsReco = nNeutralHadronsJet[jJet];
        secondJet.nPhotonsReco = nPhotonsJet[jJet];
        secondJet.pdgIdPart = pdgIdPartJet[jJet];
        secondJet.pdgIdMom = pdgIdMomJet[jJet];

        TLorentzVector thisDiJet = firstJet + secondJet;
        float thisMjj = thisDiJet.M();

//        if( (fabs(firstJet.pdgIdMom)==23 && fabs(secondJet.pdgIdMom)==23) || (fabs(firstJet.pdgIdMom)==24 && fabs(secondJet.pdgIdMom)==24) ) {

        if( bestMjj<0. || fabs(thisMjj-85.) < fabs(bestMjj-85.) ) {

          bestMjj = thisMjj;
          jet1 = firstJet;
          jet2 = secondJet;

          if( (fabs(firstJet.pdgIdMom)==23 && fabs(secondJet.pdgIdMom)==23) || (fabs(firstJet.pdgIdMom)==24 && fabs(secondJet.pdgIdMom)==24) ) correctPair=true;
          else correctPair=false;

        } //if best mass

      } //for j jets

    } //for i jets


    if( bestMjj<0. ) { //didnt find a suitable pair
      continue;
    }

    if( correctPair ) selectedCorrectPair += 1.;
    allPairs += 1.;


    std::pair<TLorentzVector,TLorentzVector> jets_kinfit = fitter_jets->fit(jet1, jet2);
    chiSquareProb = TMath::Prob(fitter_jets->getS(), fitter_jets->getNDF());
    h1_kinfit_chiSquare->Fill( fitter_jets->getS()/fitter_jets->getNDF(), eventWeight ); 
    h1_kinfit_chiSquareProb->Fill( chiSquareProb, eventWeight ); 


    float QGLikelihoodJet1 = qglikeli->computeQGLikelihoodPU( jet1.Pt(), rhoPF, jet1.nCharged(), jet1.nNeutral(), jet1.ptD );
    float QGLikelihoodJet2 = qglikeli->computeQGLikelihoodPU( jet2.Pt(), rhoPF, jet2.nCharged(), jet2.nNeutral(), jet2.ptD );

    h1_ptJet1->Fill( jet1.Pt(), eventWeight );
    h1_etaJet1->Fill( jet1.Eta(), eventWeight );
    h1_pdgIdJet1->Fill( jet1.pdgIdPart, eventWeight );
    h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );

    h1_ptJet2->Fill( jet2.Pt(), eventWeight );
    h1_etaJet2->Fill( jet2.Eta(), eventWeight );
    h1_pdgIdJet2->Fill( jet2.pdgIdPart, eventWeight );
    h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );



    TLorentzVector diJet = jet1 + jet2;
    h1_mjj->Fill( diJet.M(), eventWeight );

    h1_QGLikelihoodProd->Fill( QGLikelihoodJet1*QGLikelihoodJet2, eventWeight );


    // set tree branch values
    ptJet1_t = jet1.Pt();
    ptJet2_t = jet2.Pt();

    etaJet1_t = jet1.Eta();
    etaJet2_t = jet2.Eta();

    QGLikelihoodJet1_t = QGLikelihoodJet1;
    QGLikelihoodJet2_t = QGLikelihoodJet2;

    mgg = diPhoton.M();
    mjj = diJet.M();

    // fill tree:
    tree_passedEvents->Fill();

  } //for entries


  std::cout << std::endl << std::endl;
  std::cout << "Finished loop." << std::endl;
  std::cout << "Chose correct pair in " << 100.*selectedCorrectPair/allPairs << "% of events" << std::endl;
  

  outFile->cd();

  tree_passedEvents->Write();

  h1_nvertex->Write();
  h1_nvertexPU->Write();

  h1_ptPhot1->Write();
  h1_etaPhot1->Write();

  h1_ptPhot2->Write();
  h1_etaPhot2->Write();

  h1_mgg->Write();

  h1_kinfit_chiSquare->Write();
  h1_kinfit_chiSquareProb->Write();

  h1_ptJet1->Write();
  h1_etaJet1->Write();
  h1_pdgIdJet1->Write();
  h1_QGLikelihoodJet1->Write();

  h1_ptJet2->Write();
  h1_etaJet2->Write();
  h1_pdgIdJet2->Write();
  h1_QGLikelihoodJet2->Write();

  h1_mjj->Write();
  h1_QGLikelihoodProd->Write();


  outFile->Close();


  delete tree;
  tree = 0;

  delete h1_nvertex;
  delete h1_nvertexPU;
  delete h1_ptPhot1;
  delete h1_etaPhot1;
  delete h1_ptPhot2;
  delete h1_etaPhot2;
  delete h1_mgg;
  delete h1_kinfit_chiSquare;
  delete h1_kinfit_chiSquareProb;
  delete h1_ptJet1;
  delete h1_etaJet1;
  delete h1_pdgIdJet1;
  delete h1_QGLikelihoodJet1;
  delete h1_ptJet2;
  delete h1_etaJet2;
  delete h1_pdgIdJet2;
  delete h1_QGLikelihoodJet2;
  delete h1_mjj;
  delete h1_QGLikelihoodProd;


  totalLumi = 0.;

}


void addInput( const std::string& dataset ) {


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "VHgg_2ndLevelTreeW_" + dataset + "_pfakt5";
    infileName += ".root";
    TFile* infile = TFile::Open(infileName.c_str(), "read");
    if( infile==0 ) {
      std::cout << "Didn't find file '" << infileName << "'. Did you forget to finalize (i.e. the \"W\")?" << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(77);
    }
    std::string treeName = infileName +"/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infile->Close();

  /*} else {

    char singleLine[500];

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/jetTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
      }
      infile->Close();

    }
    fclose(iff);

  } */

} //addinput

