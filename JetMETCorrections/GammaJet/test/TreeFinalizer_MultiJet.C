
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
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

#include "AnalysisJet.C"
#include "AnalysisPhoton.C"

#include "/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/CommonTools/PUWeight.C"
#include "/cmsrm/pc25/pandolf/CMSSW_4_2_8_patch7/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C"


bool isAOD_ = true;
bool DEBUG_ = false;
bool MCassoc_ = false;
std::string RECOTYPE_;
std::string ALGOTYPE_;
std::string PARTTYPE_;

TChain* tree;
Double_t totalLumi=0.;


void addInput(const std::string& dataset);



void finalize(const std::string& dataset, bool dijet_selection=false) {

  TString dataset_tstr(dataset);


  tree = new TChain("jetTree");


  if( dataset=="QCD_HT_Summer11" ) {
    addInput( "QCD_TuneZ2_HT-100To250_7TeV-madgraph_Summer11-PU_S4_START42_V11-v2" );
    addInput( "QCD_TuneZ2_HT-250To500_7TeV-madgraph_Summer11-PU_S4_START42_V11-v3" );
    addInput( "QCD_TuneZ2_HT-500To1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1" );
    addInput( "QCD_TuneZ2_HT-1000_7TeV-madgraph_Summer11-PU_S4_START42_V11-v1" );
  } else if( dataset=="G_Summer11" ) {
    addInput( "G_Pt-170to300_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput( "G_Pt-300to470_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput( "G_Pt-800to1400_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput( "G_Pt-470to800_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput( "G_Pt-120to170_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput( "G_Pt-80to120_TuneZ2_7TeV_pythia6_Summer11-PU_S4_START42_V11-v1");
    addInput( "G_Pt-1400to1800_TuneZ2_7TeV_pythia_Summer11-PU_S4_START42_V11-v1");
  } else if( dataset=="HT_Run2011_FULL" ) {
    addInput( "HT_Run2011A-May10ReReco-v1_HLT" );
    addInput( "HT_Run2011A-PromptReco-v4_HLT" );
    addInput( "HT_Run2011A-PromptReco-v6_HLT" );
    addInput( "HT_Run2011B-PromptReco-v1_HLT" );
  } else {
    addInput( dataset );
  }






  TH1F* h1_nvertex = new TH1F("nvertex", "", 26, -0.5, 25.5);
  h1_nvertex->Sumw2();
  TH1F* h1_nvertexPU = new TH1F("nvertexPU", "", 26, -0.5, 25.5);
  h1_nvertexPU->Sumw2();

  TH1F* h1_rhoPF = new TH1F("rhoPF", "", 100, 0., 25. );
  h1_rhoPF->Sumw2();
  TH1F* h1_rhoPFPU = new TH1F("rhoPFPU", "", 100, 0., 25.);
  h1_rhoPFPU->Sumw2();

  TH1D* h1_ht_akt5 = new TH1D("ht_akt5", "", 300, 200., 3200.);
  h1_ht_akt5->Sumw2();
  TH1D* h1_htmet_akt5 = new TH1D("htmet_akt5", "", 300, 200., 3200.);
  h1_htmet_akt5->Sumw2();
  TH1D* h1_sumpt_pfakt5 = new TH1D("sumpt_pfakt5", "", 300, 200., 3200.);
  h1_sumpt_pfakt5->Sumw2();

  TH1D* h1_deltaR_part_jet0 = new TH1D("deltaR_part_jet0", "", 100, 0., 5.);
  h1_deltaR_part_jet0->Sumw2();
  TH1D* h1_deltaR_part_jet1 = new TH1D("deltaR_part_jet1", "", 100, 0., 5.);
  h1_deltaR_part_jet1->Sumw2();
  TH1D* h1_deltaR_part_jet2 = new TH1D("deltaR_part_jet2", "", 100, 0., 5.);
  h1_deltaR_part_jet2->Sumw2();
  TH1D* h1_deltaR_part_jet3 = new TH1D("deltaR_part_jet3", "", 100, 0., 5.);
  h1_deltaR_part_jet3->Sumw2();


  TH1D* h1_ptJet0 = new TH1D("ptJet0", "", 500, 0., 500.);
  h1_ptJet0->Sumw2();
  TH1D* h1_etaJet0 = new TH1D("etaJet0", "", 500, -5., 5.);
  h1_etaJet0->Sumw2();
  TH1D* h1_pdgIdJet0 = new TH1D("pdgIdJet0", "", 38, -15.5, 22.5);
  h1_pdgIdJet0->Sumw2();
  TH1D* h1_QGLikelihoodJet0 = new TH1D("QGLikelihoodJet0", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet0->Sumw2();
  TH1D* h1_ptDJet0 = new TH1D("ptDJet0", "", 100, 0., 1.0001);
  h1_ptDJet0->Sumw2();
  TH1D* h1_nChargedJet0 = new TH1D("nChargedJet0", "", 101, -0.5, 100.001);
  h1_nChargedJet0->Sumw2();
  TH1D* h1_nNeutralJet0 = new TH1D("nNeutralJet0", "", 101, -0.5, 100.001);
  h1_nNeutralJet0->Sumw2();

  TH1D* h1_ptJet1 = new TH1D("ptJet1", "", 500, 0., 500.);
  h1_ptJet1->Sumw2();
  TH1D* h1_etaJet1 = new TH1D("etaJet1", "", 500, -5., 5.);
  h1_etaJet1->Sumw2();
  TH1D* h1_pdgIdJet1 = new TH1D("pdgIdJet1", "", 38, -15.5, 22.5);
  h1_pdgIdJet1->Sumw2();
  TH1D* h1_QGLikelihoodJet1 = new TH1D("QGLikelihoodJet1", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet1->Sumw2();
  TH1D* h1_ptDJet1 = new TH1D("ptDJet1", "", 100, 0., 1.0001);
  h1_ptDJet1->Sumw2();
  TH1D* h1_nChargedJet1 = new TH1D("nChargedJet1", "", 101, -0.5, 100.001);
  h1_nChargedJet1->Sumw2();
  TH1D* h1_nNeutralJet1 = new TH1D("nNeutralJet1", "", 101, -0.5, 100.001);
  h1_nNeutralJet1->Sumw2();

  TH1D* h1_ptJet2 = new TH1D("ptJet2", "", 500, 0., 500.);
  h1_ptJet2->Sumw2();
  TH1D* h1_etaJet2 = new TH1D("etaJet2", "", 500, -5., 5.);
  h1_etaJet2->Sumw2();
  TH1D* h1_pdgIdJet2 = new TH1D("pdgIdJet2", "", 38, -15.5, 22.5);
  h1_pdgIdJet2->Sumw2();
  TH1D* h1_QGLikelihoodJet2 = new TH1D("QGLikelihoodJet2", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet2->Sumw2();
  TH1D* h1_ptDJet2 = new TH1D("ptDJet2", "", 100, 0., 1.0001);
  h1_ptDJet2->Sumw2();
  TH1D* h1_nChargedJet2 = new TH1D("nChargedJet2", "", 101, -0.5, 100.001);
  h1_nChargedJet2->Sumw2();
  TH1D* h1_nNeutralJet2 = new TH1D("nNeutralJet2", "", 101, -0.5, 100.001);
  h1_nNeutralJet2->Sumw2();

  TH1D* h1_ptJet3 = new TH1D("ptJet3", "", 500, 0., 500.);
  h1_ptJet3->Sumw2();
  TH1D* h1_etaJet3 = new TH1D("etaJet3", "", 500, -5., 5.);
  h1_etaJet3->Sumw2();
  TH1D* h1_pdgIdJet3 = new TH1D("pdgIdJet3", "", 38, -15.5, 22.5);
  h1_pdgIdJet3->Sumw2();
  TH1D* h1_QGLikelihoodJet3 = new TH1D("QGLikelihoodJet3", "", 100, 0., 1.0001);
  h1_QGLikelihoodJet3->Sumw2();
  TH1D* h1_ptDJet3 = new TH1D("ptDJet3", "", 100, 0., 1.0001);
  h1_ptDJet3->Sumw2();
  TH1D* h1_nChargedJet3 = new TH1D("nChargedJet3", "", 101, -0.5, 100.001);
  h1_nChargedJet3->Sumw2();
  TH1D* h1_nNeutralJet3 = new TH1D("nNeutralJet3", "", 101, -0.5, 100.001);
  h1_nNeutralJet3->Sumw2();




  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t LS;
  tree->SetBranchAddress("LS", &LS);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Float_t rhoPF;
  tree->SetBranchAddress("rhoPF", &rhoPF);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Float_t eMet;
  tree->SetBranchAddress("eMet", &eMet);
  Float_t ht_akt5;
  tree->SetBranchAddress("ht_akt5", &ht_akt5);
  Float_t epfMet;
  tree->SetBranchAddress("epfMet", &epfMet);
  Float_t phiMet;
  tree->SetBranchAddress("phipfMet", &phiMet);


  Int_t nPU;
  tree->SetBranchAddress("nPU", &nPU);
  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);


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
  Float_t eNeutralHadronsJet[20];
  tree->SetBranchAddress("eNeutralHadronsJet", eNeutralHadronsJet);
  Int_t nChargedHadronsJet[20];
  tree->SetBranchAddress("nChargedHadronsJet", nChargedHadronsJet);
  Int_t nPhotonsJet[20];
  tree->SetBranchAddress("nPhotonsJet", nPhotonsJet);
  Int_t nNeutralHadronsJet[20];
  tree->SetBranchAddress("nNeutralHadronsJet", nNeutralHadronsJet);
  Float_t ptDJet[20];
  tree->SetBranchAddress("ptDJet", ptDJet);
  Float_t QGLikelihoodJet[20];
  tree->SetBranchAddress("QGLikelihoodJet", QGLikelihoodJet);

  Float_t ePartJet[20];
  tree->SetBranchAddress("ePartJet", ePartJet);
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


  Bool_t passed_HT150;
  tree->SetBranchAddress("passed_HT150", &passed_HT150);
  Bool_t passed_HT200;
  tree->SetBranchAddress("passed_HT200", &passed_HT200);
  Bool_t passed_HT250;
  tree->SetBranchAddress("passed_HT250", &passed_HT250);
  Bool_t passed_HT300;
  tree->SetBranchAddress("passed_HT300", &passed_HT300);
  Bool_t passed_HT350;
  tree->SetBranchAddress("passed_HT350", &passed_HT350);
  Bool_t passed_HT400;
  tree->SetBranchAddress("passed_HT400", &passed_HT400);
  Bool_t passed_HT450;
  tree->SetBranchAddress("passed_HT450", &passed_HT450);
  Bool_t passed_HT500;
  tree->SetBranchAddress("passed_HT500", &passed_HT500);
  Bool_t passed_HT550;
  tree->SetBranchAddress("passed_HT550", &passed_HT550);
  Bool_t passed_HT600;
  tree->SetBranchAddress("passed_HT600", &passed_HT600);
  Bool_t passed_HT650;
  tree->SetBranchAddress("passed_HT650", &passed_HT650);
  Bool_t passed_HT700;
  tree->SetBranchAddress("passed_HT700", &passed_HT700);




  std::string puType = "Spring11_Flat10";
  if( dataset_tstr.Contains("Summer11") ) puType = "Summer11_S4";
  PUWeight* fPUWeight = new PUWeight(-1, "2011A", puType);
  std::string puFileName;
  //if( PUType_=="Run2011A_73pb" )
  //puFileName = "all2011A.pileup_v2.root";
  //puFileName = "PUTarget.Run2011B.175832-180252.root";
    puFileName = "PileUpTarget_runs_173236_178380.root";
  std::cout << std::endl << "-> Using data pileup file: " << puFileName << std::endl;
  TFile* filePU = TFile::Open(puFileName.c_str());
  TH1F* h1_nPU_data = (TH1F*)filePU->Get("pileup");
  fPUWeight->SetDataHistogram(h1_nPU_data);

  PUWeight* fPUWeightRunA = new PUWeight(-1, "2011A", puType);
  std::string puFileNameRunA = "all2011A.pileup_v2_73mb.root";
  TFile* filePURunA = TFile::Open(puFileNameRunA.c_str());
  TH1F* h1_nPU_dataRunA = (TH1F*)filePURunA->Get("pileup");
  fPUWeightRunA->SetDataHistogram(h1_nPU_dataRunA);






  std::string outfileName;
  if( DEBUG_ ) outfileName = "provaMultiJet_"+dataset;
  else {
   if(dataset!="") outfileName = "MultiJet_"+dataset;
   else outfileName = "MultiJet";
  }
  if( MCassoc_ ) outfileName = outfileName + "_MCassoc";
  if( dijet_selection ) outfileName = outfileName + "_DIJET";
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();




  float ptJet0, ptJet1, ptJet2, ptJet3;
  float etaJet0, etaJet1, etaJet2, etaJet3;
  int pdgIdPartJet0, pdgIdPartJet1, pdgIdPartJet2, pdgIdPartJet3;
  float QGLikelihoodJet0, QGLikelihoodJet1, QGLikelihoodJet2, QGLikelihoodJet3;
  float ptDJet0, ptDJet1, ptDJet2, ptDJet3;
  int nChargedJet0, nChargedJet1, nChargedJet2, nChargedJet3;
  int nNeutralJet0, nNeutralJet1, nNeutralJet2, nNeutralJet3;
  float eventWeight_noPU;
  float PUWeight(1.), PUWeightRunA(1.);


  TTree* tree_passedEvents = new TTree("tree_passedEvents", "");
  tree_passedEvents->Branch( "run", &run, "run/I" );
  tree_passedEvents->Branch( "event", &event, "event/I" );
  tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
  tree_passedEvents->Branch( "eventWeight_noPU", &eventWeight_noPU, "eventWeight_noPU/F" );
  tree_passedEvents->Branch( "PUWeight", &PUWeight, "PUWeight/F" );
  tree_passedEvents->Branch( "PUWeightRunA", &PUWeightRunA, "PUWeightRunA/F" );
  tree_passedEvents->Branch( "nvertex", &nvertex, "nvertex/I" );
  tree_passedEvents->Branch( "rhoPF", &rhoPF, "rhoPF/F" );
  tree_passedEvents->Branch( "ht_akt5", &ht_akt5, "ht_akt5/F" );
  
  tree_passedEvents->Branch("passed_HT150", &passed_HT150, "passed_HT150/O");
  tree_passedEvents->Branch("passed_HT200", &passed_HT200, "passed_HT200/O");
  tree_passedEvents->Branch("passed_HT250", &passed_HT250, "passed_HT250/O");
  tree_passedEvents->Branch("passed_HT300", &passed_HT300, "passed_HT300/O");
  tree_passedEvents->Branch("passed_HT350", &passed_HT350, "passed_HT350/O");
  tree_passedEvents->Branch("passed_HT400", &passed_HT400, "passed_HT400/O");
  tree_passedEvents->Branch("passed_HT450", &passed_HT450, "passed_HT450/O");
  tree_passedEvents->Branch("passed_HT500", &passed_HT500, "passed_HT500/O");
  tree_passedEvents->Branch("passed_HT550", &passed_HT550, "passed_HT550/O");
  tree_passedEvents->Branch("passed_HT600", &passed_HT600, "passed_HT600/O");
  tree_passedEvents->Branch("passed_HT650", &passed_HT650, "passed_HT650/O");
  tree_passedEvents->Branch("passed_HT700", &passed_HT700, "passed_HT700/O");

  tree_passedEvents->Branch( "ptJet0", &ptJet0, "ptJet0/F" );
  tree_passedEvents->Branch( "ptJet1", &ptJet1, "ptJet1/F" );
  tree_passedEvents->Branch( "ptJet2", &ptJet2, "ptJet2/F" );
  tree_passedEvents->Branch( "ptJet3", &ptJet3, "ptJet3/F" );

  tree_passedEvents->Branch( "etaJet0", &etaJet0, "etaJet0/F" );
  tree_passedEvents->Branch( "etaJet1", &etaJet1, "etaJet1/F" );
  tree_passedEvents->Branch( "etaJet2", &etaJet2, "etaJet2/F" );
  tree_passedEvents->Branch( "etaJet3", &etaJet3, "etaJet3/F" );

  tree_passedEvents->Branch( "pdgIdPartJet0", &pdgIdPartJet0, "pdgIdPartJet0/I" );
  tree_passedEvents->Branch( "pdgIdPartJet1", &pdgIdPartJet1, "pdgIdPartJet1/I" );
  tree_passedEvents->Branch( "pdgIdPartJet2", &pdgIdPartJet2, "pdgIdPartJet2/I" );
  tree_passedEvents->Branch( "pdgIdPartJet3", &pdgIdPartJet3, "pdgIdPartJet3/I" );

  tree_passedEvents->Branch( "QGLikelihoodJet0", &QGLikelihoodJet0, "QGLikelihoodJet0/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet1", &QGLikelihoodJet1, "QGLikelihoodJet1/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet2", &QGLikelihoodJet2, "QGLikelihoodJet2/F" );
  tree_passedEvents->Branch( "QGLikelihoodJet3", &QGLikelihoodJet3, "QGLikelihoodJet3/F" );

  tree_passedEvents->Branch( "ptDJet0", &ptDJet0, "ptDJet0/F" );
  tree_passedEvents->Branch( "ptDJet1", &ptDJet1, "ptDJet1/F" );
  tree_passedEvents->Branch( "ptDJet2", &ptDJet2, "ptDJet2/F" );
  tree_passedEvents->Branch( "ptDJet3", &ptDJet3, "ptDJet3/F" );

  tree_passedEvents->Branch( "nChargedJet0", &nChargedJet0, "nChargedJet0/I" );
  tree_passedEvents->Branch( "nChargedJet1", &nChargedJet1, "nChargedJet1/I" );
  tree_passedEvents->Branch( "nChargedJet2", &nChargedJet2, "nChargedJet2/I" );
  tree_passedEvents->Branch( "nChargedJet3", &nChargedJet3, "nChargedJet3/I" );

  tree_passedEvents->Branch( "nNeutralJet0", &nNeutralJet0, "nNeutralJet0/I" );
  tree_passedEvents->Branch( "nNeutralJet1", &nNeutralJet1, "nNeutralJet1/I" );
  tree_passedEvents->Branch( "nNeutralJet2", &nNeutralJet2, "nNeutralJet2/I" );
  tree_passedEvents->Branch( "nNeutralJet3", &nNeutralJet3, "nNeutralJet3/I" );



  gROOT->cd();



  bool debug=true;

  std::map< int, std::map<int, std::vector<int> > > run_lumi_ev_map;

  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    bool isMC = run<5;


    if( !isMC ) {

      std::map<int, std::map<int, std::vector<int> > >::iterator it;

      it = run_lumi_ev_map.find(run);


      if( it==run_lumi_ev_map.end() ) {

        std::vector<int> events;
        events.push_back(event);
        std::map<int, std::vector<int> > lumi_ev_map;
        lumi_ev_map.insert( std::pair<int,std::vector<int> >(LS, events));
        run_lumi_ev_map.insert( std::pair<int, std::map<int, std::vector<int> > > (run, lumi_ev_map) );

      } else { //run exists, look for LS


        std::map<int, std::vector<int> >::iterator it_LS;
        it_LS = it->second.find( LS );

        if( it_LS==(it->second.end())  ) {

          std::vector<int> events;
          events.push_back(event);
          it->second.insert( std::pair<int, std::vector<int> > (LS, events) );

        } else { //LS exists, look for event

          std::vector<int>::iterator ev;
          for( ev=it_LS->second.begin(); ev!=it_LS->second.end(); ++ev )
            if( *ev==event ) break;


          if( ev==it_LS->second.end() ) {

            it_LS->second.push_back(event);

          } else {

            std::cout << "DISCARDING DUPLICATE EVENT!! Run: " << run << " LS: " << LS << " event: " << event << std::endl;

            continue;

          }
        }
      }

    } //if is mc


    //initialize tree branches:
    ptJet0 = 0.;
    etaJet0 = 10.;
    pdgIdPartJet0 = -100;
    QGLikelihoodJet0 = -1.;
    ptDJet0 = -1.;
    nChargedJet0 = -1.;
    nNeutralJet0 = -1.;

    ptJet1 = 0.;
    etaJet1 = 10.;
    pdgIdPartJet1 = -100;
    QGLikelihoodJet1 = -1.;
    ptDJet1 = -1.;
    nChargedJet1 = -1.;
    nNeutralJet1 = -1.;

    ptJet2 = 0.;
    etaJet2 = 10.;
    pdgIdPartJet2 = -100;
    QGLikelihoodJet2 = -1.;
    ptDJet2 = -1.;
    nChargedJet2 = -1.;
    nNeutralJet2 = -1.;

    ptJet3 = 0.;
    etaJet3 = 10.;
    pdgIdPartJet3 = -100;
    QGLikelihoodJet3 = -1.;
    ptDJet3 = -1.;
    nChargedJet3 = -1.;
    nNeutralJet3 = -1.;



    if( ht_akt5 > 3500. ) continue;


    if( dijet_selection ) {

      if( ht_akt5 < 150. ) continue;

    } else {

      if( !isMC ) {
        if( !passed_HT600 ) continue; //trigger on data
        if( run<173236 || run>178380 ) continue; //run range in which HT600 was unprescaled: corresponds to 1894.3 pb-1
      }
  
      // avoid trigger turn-on:
      if( ht_akt5 < 650. ) continue;

    }
  

    if( dijet_selection ) {

      if( nJet<2 ) continue;
    
    } else {

      if( nJet<4 ) continue;
      if( ptJet[3]<20. ) continue; //speed it up a little

    }

 
    if( dijet_selection && nJet>2 ) {
      float ptAve = 0.5*(ptJet[0] + ptJet[1]);
      if( ptJet[2] > 0.5*ptAve ) continue;
    }


    if( eventWeight <= 0. ) eventWeight = 1.;

    h1_nvertex->Fill( nvertex, eventWeight);
    h1_rhoPF->Fill( rhoPF, eventWeight);

    if( isMC ) {
      // PU reweighting:
     eventWeight_noPU = eventWeight;
     PUWeight = fPUWeight->GetWeight(nPU);
     PUWeightRunA = fPUWeightRunA->GetWeight(nPU);
     eventWeight *= PUWeight;
    }

    h1_nvertexPU->Fill( nvertex, eventWeight);
    h1_rhoPFPU->Fill( rhoPF, eventWeight);

    h1_ht_akt5->Fill( ht_akt5, eventWeight );
    h1_htmet_akt5->Fill( ht_akt5 + eMet, eventWeight );


    std::vector<AnalysisJet*> jets;

    //for( unsigned iJet=0; iJet<1 && jets.size()<4; ++iJet ) {
    for( unsigned iJet=0; iJet<nJet && jets.size()<4; ++iJet ) {

      if( ptJet[iJet]<20. ) continue;
      if( fabs(etaJet[iJet])>2. ) continue;
      //int nCandidates = nChargedHadronsJet[iJet] + nNeutralHadronsJet[iJet] + nPhotonsJet[iJet];
      //if( nChargedHadronsJet[iJet]==0 || nCandidates<2 ) continue;

      AnalysisJet* thisJet = new AnalysisJet();
    
      thisJet->SetPtEtaPhiE( ptJet[iJet], etaJet[iJet], phiJet[iJet], ptJet[iJet]/ptRawJet[iJet]*eJet[iJet] );

      thisJet->QGLikelihood = QGLikelihoodJet[iJet];
      thisJet->ptD = ptDJet[iJet];
      thisJet->nTracksReco = nChargedHadronsJet[iJet];
      thisJet->nNeutralHadronsReco = nNeutralHadronsJet[iJet];
      thisJet->nPhotonsReco = nPhotonsJet[iJet];
      thisJet->eTracksReco = eChargedHadronsJet[iJet];
      thisJet->eNeutralHadronsReco = eNeutralHadronsJet[iJet];
      thisJet->ePhotonsReco = ePhotonsJet[iJet];

      if( !(thisJet->passedJetID()) ) continue;

      if( isMC ) {

        thisJet->pdgIdPart = pdgIdPartJet[iJet];
        thisJet->ptPart = ptPartJet[iJet];
        thisJet->etaPart = etaPartJet[iJet];
        thisJet->phiPart = phiPartJet[iJet];
        thisJet->ePart = ePartJet[iJet];
        thisJet->pdgIdMom = pdgIdMomJet[iJet];

        TLorentzVector* parton = new TLorentzVector();
        parton->SetPtEtaPhiE( thisJet->ptPart, thisJet->etaPart, thisJet->phiPart, thisJet->ePart );

        float deltaR_part_jet = thisJet->DeltaR( *parton );
        if( deltaR_part_jet>0.5 ) thisJet->pdgIdPart = -100;

        if( iJet==0 )
          h1_deltaR_part_jet0->Fill( deltaR_part_jet, eventWeight );
        else if( iJet==1 )
          h1_deltaR_part_jet1->Fill( deltaR_part_jet, eventWeight );
        else if( iJet==2 )
          h1_deltaR_part_jet2->Fill( deltaR_part_jet, eventWeight );
        else if( iJet==3 )
          h1_deltaR_part_jet3->Fill( deltaR_part_jet, eventWeight );

      } else {

        thisJet->pdgIdPart = 0;
        thisJet->ptPart = 0.;
        thisJet->etaPart = 0.;
        thisJet->phiPart = 0.;
        thisJet->ePart = 0.;
        thisJet->pdgIdMom = 0;

      }


      jets.push_back(thisJet);

    } //for i jets


    if( dijet_selection ) {
      if( jets.size()<2 ) continue;
    } else {
      if( jets.size()<4 ) continue;
    }


    float sumpt=0.;
    for( unsigned iJet=0; iJet<jets.size(); ++iJet ) sumpt += jets[iJet]->Pt();
    h1_sumpt_pfakt5->Fill( sumpt, eventWeight );


    if( jets.size()>0 ) {

      ptJet0 = jets[0]->Pt();
      etaJet0 = jets[0]->Eta();
      pdgIdPartJet0 = jets[0]->pdgIdPart;
      QGLikelihoodJet0 = jets[0]->QGLikelihood;
      ptDJet0 = jets[0]->ptD;
      nChargedJet0 = jets[0]->nCharged();
      nNeutralJet0 = jets[0]->nNeutral();

      h1_ptJet0->Fill( jets[0]->Pt(), eventWeight );
      h1_etaJet0->Fill( jets[0]->Eta(), eventWeight );
      h1_pdgIdJet0->Fill( jets[0]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet0->Fill( QGLikelihoodJet0, eventWeight );
      h1_ptDJet0->Fill( jets[0]->ptD, eventWeight );
      h1_nChargedJet0->Fill( jets[0]->nCharged(), eventWeight );
      h1_nNeutralJet0->Fill( jets[0]->nNeutral(), eventWeight );

    }

    if( jets.size()>1 ) {

      ptJet1 = jets[1]->Pt();
      etaJet1 = jets[1]->Eta();
      pdgIdPartJet1 = jets[1]->pdgIdPart;
      QGLikelihoodJet1 = jets[1]->QGLikelihood;
      ptDJet1 = jets[1]->ptD;
      nChargedJet1 = jets[1]->nCharged();
      nNeutralJet1 = jets[1]->nNeutral();

      h1_ptJet1->Fill( jets[1]->Pt(), eventWeight );
      h1_etaJet1->Fill( jets[1]->Eta(), eventWeight );
      h1_pdgIdJet1->Fill( jets[1]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet1->Fill( QGLikelihoodJet1, eventWeight );
      h1_ptDJet1->Fill( jets[1]->ptD, eventWeight );
      h1_nChargedJet1->Fill( jets[1]->nCharged(), eventWeight );
      h1_nNeutralJet1->Fill( jets[1]->nNeutral(), eventWeight );

    }

    if( jets.size()>2 ) {

      ptJet2 = jets[2]->Pt();
      etaJet2 = jets[2]->Eta();
      pdgIdPartJet2 = jets[2]->pdgIdPart;
      QGLikelihoodJet2 = jets[2]->QGLikelihood;
      ptDJet2 = jets[2]->ptD;
      nChargedJet2 = jets[2]->nCharged();
      nNeutralJet2 = jets[2]->nNeutral();

      h1_ptJet2->Fill( jets[2]->Pt(), eventWeight );
      h1_etaJet2->Fill( jets[2]->Eta(), eventWeight );
      h1_pdgIdJet2->Fill( jets[2]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet2->Fill( QGLikelihoodJet2, eventWeight );
      h1_ptDJet2->Fill( jets[2]->ptD, eventWeight );
      h1_nChargedJet2->Fill( jets[2]->nCharged(), eventWeight );
      h1_nNeutralJet2->Fill( jets[2]->nNeutral(), eventWeight );

    }

    if( jets.size()>3 ) {

      ptJet3 = jets[3]->Pt();
      etaJet3 = jets[3]->Eta();
      pdgIdPartJet3 = jets[3]->pdgIdPart;
      QGLikelihoodJet3 = jets[3]->QGLikelihood;
      ptDJet3 = jets[3]->ptD;
      nChargedJet3 = jets[3]->nCharged();
      nNeutralJet3 = jets[3]->nNeutral();

      h1_ptJet3->Fill( jets[3]->Pt(), eventWeight );
      h1_etaJet3->Fill( jets[3]->Eta(), eventWeight );
      h1_pdgIdJet3->Fill( jets[3]->pdgIdPart, eventWeight );
      h1_QGLikelihoodJet3->Fill( QGLikelihoodJet3, eventWeight );
      h1_ptDJet3->Fill( jets[3]->ptD, eventWeight );
      h1_nChargedJet3->Fill( jets[3]->nCharged(), eventWeight );
      h1_nNeutralJet3->Fill( jets[3]->nNeutral(), eventWeight );

    }


    // fill tree:
    tree_passedEvents->Fill();

  } //for entries


  std::cout << std::endl << std::endl;
  std::cout << "Finished loop." << std::endl;
  

  outFile->cd();

  tree_passedEvents->Write();

  h1_nvertex->Write();
  h1_nvertexPU->Write();
  h1_rhoPF->Write();
  h1_rhoPFPU->Write();

  h1_ht_akt5->Write();
  h1_htmet_akt5->Write();
  h1_sumpt_pfakt5->Write();

  h1_deltaR_part_jet0->Write();
  h1_deltaR_part_jet1->Write();
  h1_deltaR_part_jet2->Write();
  h1_deltaR_part_jet3->Write();

  h1_ptJet0->Write();
  h1_etaJet0->Write();
  h1_pdgIdJet0->Write();
  h1_QGLikelihoodJet0->Write();
  h1_ptDJet0->Write();
  h1_nChargedJet0->Write();
  h1_nNeutralJet0->Write();

  h1_ptJet1->Write();
  h1_etaJet1->Write();
  h1_pdgIdJet1->Write();
  h1_QGLikelihoodJet1->Write();
  h1_ptDJet1->Write();
  h1_nChargedJet1->Write();
  h1_nNeutralJet1->Write();

  h1_ptJet2->Write();
  h1_etaJet2->Write();
  h1_pdgIdJet2->Write();
  h1_QGLikelihoodJet2->Write();
  h1_ptDJet2->Write();
  h1_nChargedJet2->Write();
  h1_nNeutralJet2->Write();

  h1_ptJet3->Write();
  h1_etaJet3->Write();
  h1_pdgIdJet3->Write();
  h1_QGLikelihoodJet3->Write();
  h1_ptDJet3->Write();
  h1_nChargedJet3->Write();
  h1_nNeutralJet3->Write();


  outFile->Close();


  delete tree;
  tree = 0;

  delete h1_nvertex;
  delete h1_nvertexPU;
  delete h1_rhoPF;
  delete h1_rhoPFPU;

  delete h1_ptJet0;
  delete h1_etaJet0;
  delete h1_pdgIdJet0;
  delete h1_QGLikelihoodJet0;

  delete h1_ptJet1;
  delete h1_etaJet1;
  delete h1_pdgIdJet1;
  delete h1_QGLikelihoodJet1;

  delete h1_ptJet2;
  delete h1_etaJet2;
  delete h1_pdgIdJet2;
  delete h1_QGLikelihoodJet2;

  delete h1_ptJet3;
  delete h1_etaJet3;
  delete h1_pdgIdJet3;
  delete h1_QGLikelihoodJet3;


  totalLumi = 0.;

}


void addInput( const std::string& dataset ) {


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "MultiJet_2ndLevelTreeW_" + dataset + "_pfakt5";
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

