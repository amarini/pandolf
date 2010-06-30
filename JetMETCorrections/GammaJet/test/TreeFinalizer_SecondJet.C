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
#include <iostream>
#include <vector>
#include <cmath>

#include "fitTools.C"




bool DEBUG_ = false;
bool MCASSOC_ = false;
bool NOEVENTWEIGHTS_ = false;
bool NOHEMISPHERE_  = true;
std::string RECOTYPE_;
std::string PARTONTYPE = "ALL"; //this can be: "ALL", "GLUON", "QUARK"
bool USEGENJETS_ = true; //for QCD background


int NBINS_PT = 16;



std::vector< std::vector< TH1D* > > allocateHistoVector(const std::string& name, const std::string& abs_rel, int nPoints, bool is_pt=false);
void deleteHistoVector(std::vector< std::vector< TH1D* > > histoVector, int nPoints);


void secondJet(const std::string& dataset, std::string recoType, std::string jetAlgo, bool nohemisphere=true, std::string partonType="ALL") {

  RECOTYPE_ = recoType;

  std::string algoType = (recoType=="pf") ? recoType+jetAlgo : jetAlgo;

  NOHEMISPHERE_ = nohemisphere; 
  PARTONTYPE = partonType;

  TChain* tree = new TChain("jetTree");
  
  std::string infileName;

  std::string prefix = "PhotonJet_2ndLevelTree_";
  std::string suffix = "_"+algoType;
  std::string suffix_QCD = (USEGENJETS_) ? "_GENJETS_CORR" : "";


  if( dataset=="PhotonJet") {

    tree->Add("PhotonJet_2ndLevelTree_PhotonJet15.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet30.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet80.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet170.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet300.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet470.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet800.root/jetTree");
    tree->Add("PhotonJet_2ndLevelTree_PhotonJet1400.root/jetTree");

  } else if( dataset=="PhotonJet_Summer09") {

    infileName = prefix + "PhotonJet_Summer09_Pt15" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt30" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt80" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt170" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt300" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt470" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt800" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt1400" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="QCD_Summer09") {

    infileName = prefix + "QCD_Summer09_Pt15" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt30" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt80" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt170" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt300" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt470" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt800" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="PhotonJetQCD_Summer09") {

    infileName = prefix + "PhotonJet_Summer09_Pt15" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt30" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt80" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt170" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt300" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt470" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt800" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt1400" + suffix + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt15" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt30" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt80" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt170" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt300" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt470" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt800" + suffix + suffix_QCD + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else {

    if( dataset=="" )
      infileName = "PhotonJet_2ndLevelTree.root/jetTree";
    else
      infileName = prefix + dataset + suffix + ".root/jetTree";

    std::cout << "Opening file: " << infileName << std::endl;
    tree->Add(infileName.c_str());

  }

  if( NOEVENTWEIGHTS_ ) suffix = "_NOWEIGHTS" + suffix;


  bool draw=true;


//TH1D* h1_deltaR_phot_2ndJet = new TH1D("deltaR_phot_2ndJet", "", 100, 0., 1.);



  Double_t Lower[NBINS_PT];
  fitTools::getPtBins(NBINS_PT, Lower);



  TH1D* h1_ptPhot = new TH1D("ptPhot", "", 1390, 10., 1400.);
  TH1D* h1_ptPhot_cutOn2ndJet = new TH1D("ptPhot_cutOn2ndJet", "", 1390, 10., 1400.);
  TH2D* h2_pt2ndJetRel_vs_ptPhot = new TH2D("pt2ndJetRel_vs_ptPhot", "", NBINS_PT-1, Lower, 100, 0., 1.);
  h2_pt2ndJetRel_vs_ptPhot->Sumw2();

//std::vector<TH1D*> h1_pt2ndJet;
//for( int i_ptBin=0; i_ptBin<NBINS_PT; ++i_ptBin ) {
//    char histName[100];
//    sprintf(histName, "pt2ndJet_%d", i_ptBin);
//    Float_t maxPercent = ( Lower[i_ptBin+1] < 410. ) ? 0.3 : 0.15;
//    int histMax = (int)(Lower[i_ptBin+1]*maxPercent);
//    TH1D* tmp = new TH1D(histName, "", histMax, 0., (Float_t)histMax);
//    tmp->Sumw2();
//    h1_pt2ndJet.push_back(tmp);
//} //for ipt bins


  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);
  Bool_t matchedToMC;
  tree->SetBranchAddress("matchedToMC", &matchedToMC);
  Bool_t passedPhotonID_medium;
  tree->SetBranchAddress("passedPhotonID_medium", &passedPhotonID_medium);

  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t ePhotReco;
  tree->SetBranchAddress("ePhotReco", &ePhotReco);
  Float_t ptPhotReco;
  tree->SetBranchAddress("ptPhotReco", &ptPhotReco);
  Float_t etaPhotReco;
  tree->SetBranchAddress("etaPhotReco", &etaPhotReco);
  Float_t phiPhotReco;
  tree->SetBranchAddress("phiPhotReco", &phiPhotReco);

  Float_t ePhotGen;
  tree->SetBranchAddress("ePhotGen", &ePhotGen);
  Float_t ptPhotGen;
  tree->SetBranchAddress("ptPhotGen", &ptPhotGen);
  Float_t etaPhotGen;
  tree->SetBranchAddress("etaPhotGen", &etaPhotGen);
  Float_t phiPhotGen;
  tree->SetBranchAddress("phiPhotGen", &phiPhotGen);

  Float_t hcalIsoPhotReco;
  tree->SetBranchAddress("hcalIsoPhotReco", &hcalIsoPhotReco);
  Float_t ecalIsoPhotReco;
  tree->SetBranchAddress("ecalIsoPhotReco", &ecalIsoPhotReco);
  Int_t nTrkIsoPhotReco;
  tree->SetBranchAddress("nTrkIsoPhotReco", &nTrkIsoPhotReco);
  Float_t ptTrkIsoPhotReco;
  tree->SetBranchAddress("ptTrkIsoPhotReco", &ptTrkIsoPhotReco);
  Float_t clusterMajPhotReco;
  tree->SetBranchAddress("clusterMajPhotReco", &clusterMajPhotReco);
  Float_t clusterMinPhotReco;
  tree->SetBranchAddress("clusterMinPhotReco", &clusterMinPhotReco);

  Float_t eJetReco;
  tree->SetBranchAddress("eJetReco", &eJetReco);
  Float_t ptJetReco;
  tree->SetBranchAddress("ptJetReco", &ptJetReco);
  Float_t etaJetReco;
  tree->SetBranchAddress("etaJetReco", &etaJetReco);
  Float_t phiJetReco;
  tree->SetBranchAddress("phiJetReco", &phiJetReco);
  Float_t eTracksReco;
  tree->SetBranchAddress("eTracksReco", &eTracksReco);

  Float_t eJetGen;
  tree->SetBranchAddress("eJetGen", &eJetGen);
  Float_t ptJetGen;
  tree->SetBranchAddress("ptJetGen", &ptJetGen);
  Float_t etaJetGen;
  tree->SetBranchAddress("etaJetGen", &etaJetGen);
  Float_t phiJetGen;
  tree->SetBranchAddress("phiJetGen", &phiJetGen);
  Float_t eTracksGen;
  tree->SetBranchAddress("eTracksGen", &eTracksGen);

  Float_t ptPart;
  tree->SetBranchAddress("ptPart", &ptPart);
  Float_t etaPart;
  tree->SetBranchAddress("etaPart", &etaPart);
  Float_t phiPart;
  tree->SetBranchAddress("phiPart", &phiPart);
  Int_t pdgIdPart;
  tree->SetBranchAddress("pdgIdPart", &pdgIdPart);

  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);

  Float_t pt2ndJetGen;
  tree->SetBranchAddress("pt2ndJetGen", &pt2ndJetGen);
  Float_t eta2ndJetGen;
  tree->SetBranchAddress("eta2ndJetGen", &eta2ndJetGen);
  Float_t phi2ndJetGen;
  tree->SetBranchAddress("phi2ndJetGen", &phi2ndJetGen);

  Float_t ptSecondaryJetsReco;
  tree->SetBranchAddress("ptSecondaryJetsReco", &ptSecondaryJetsReco);
  Float_t ptSecondaryJetsGen;
  tree->SetBranchAddress("ptSecondaryJetsGen", &ptSecondaryJetsGen);



  Int_t nPoints = 5;

  //absolute 2nd jet pt binning histograms
  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_gen  = allocateHistoVector("RecoPhot_vs_Gen", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_gen   = allocateHistoVector("RecoGen_vs_Gen", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_gen   = allocateHistoVector("GenPart_vs_Gen", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_gen  = allocateHistoVector("PartPhot_vs_Gen", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_gen   = allocateHistoVector("GenPhot_vs_Gen", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetGenMean  = allocateHistoVector("pt2ndJetGenMean", "abs", nPoints, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_reco = allocateHistoVector("RecoPhot_vs_Reco", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_reco  = allocateHistoVector("RecoGen_vs_Reco", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_reco  = allocateHistoVector("GenPart_vs_Reco", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_reco = allocateHistoVector("PartPhot_vs_Reco", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_reco  = allocateHistoVector("GenPhot_vs_Reco", "abs", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoMean = allocateHistoVector("pt2ndJetRecoMean", "abs", nPoints, (bool)true);

  //relative 2nd jet pt binning histograms
  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_genRel  = allocateHistoVector("RecoPhot_vs_GenRel", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_genRel   = allocateHistoVector("RecoGen_vs_GenRel",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_genRel   = allocateHistoVector("GenPart_vs_GenRel",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_genRel  = allocateHistoVector("PartPhot_vs_GenRel", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_genRel   = allocateHistoVector("GenPhot_vs_GenRel",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetGenRelMean  = allocateHistoVector("pt2ndJetGenRelMean", "rel", nPoints, (bool)true);

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_recoRel = allocateHistoVector("RecoPhot_vs_RecoRel", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_recoRel  = allocateHistoVector("RecoGen_vs_RecoRel",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_recoRel  = allocateHistoVector("GenPart_vs_RecoRel",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_recoRel = allocateHistoVector("PartPhot_vs_RecoRel", "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_recoRel  = allocateHistoVector("GenPhot_vs_RecoRel",  "rel", nPoints);
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoRelMean = allocateHistoVector("pt2ndJetRecoRelMean", "rel", nPoints, (bool)true);





  TH2D* h2_response_vs_ptPhot = new TH2D("response_vs_ptPhot", "", NBINS_PT-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot->Sumw2();

  TH2D* h2_ptPhotMean_vs_pt_cutOn2ndJet = new TH2D("ptPhotMean_vs_pt_cutOn2ndJet", "", NBINS_PT-1, Lower, 2200, 0., 2200.);
  h2_ptPhotMean_vs_pt_cutOn2ndJet->Sumw2();


  //ptPhot TH2 for x-axis points on vs. pt spectra:
  TH2D* h2_ptPhotMean_vs_pt= new TH2D("ptPhotMean_vs_pt", "", NBINS_PT-1, Lower, 2200, 0., 2200.);
  h2_ptPhotMean_vs_pt->Sumw2();



  


  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;
    if( NOEVENTWEIGHTS_ )   eventWeight = 1.;

    //both photon and jet in barrel:
    if( fabs(etaPhotReco)>1.3 ) continue;
    if( fabs(etaJetReco)>1.3 ) continue;

    
    //leading jet and photon back2back in transverse plane
    bool back2back = true;
    Float_t deltaPhi_jet = phiPhotReco - phiJetReco;
    Float_t pi = 3.14159;
    if( fabs(deltaPhi_jet) < (pi - 0.2) ) back2back = false;
    //if( !back2back ) continue;  //09.11.09 deltaPhi cut is actually useless: proved it.


    if( MCASSOC_ ) {
      if( !matchedToMC ) continue;
    } else {
     if( !passedPhotonID_medium ) continue;
    }  
  
    if( PARTONTYPE=="GLUON" )
      if( fabs(pdgIdPart)!=21 ) continue;
    if( PARTONTYPE=="QUARK" )
      if( fabs(pdgIdPart)>9 ) continue;




    Double_t r_RecoPhot = (ptPhotReco!=0.) ? ptJetReco/ptPhotReco : 0.;
    Double_t r_RecoGen = (ptJetGen!=0.) ? ptJetReco/ptJetGen : 0.;
    Double_t r_GenPart = (ptPart!=0.) ? ptJetGen/ptPart : 0.;
    Double_t r_GenPhot = (ptPhotReco!=0.) ? ptJetGen/ptPhotReco : 0.;
    Double_t r_PartPhot = (ptPhotReco!=0.) ? ptPart/ptPhotReco : 0.;


    Int_t ptBin = -1;

    for( int i=0; i<NBINS_PT; ++i) {
      if( (ptPhotReco<Lower[i+1])&&(ptPhotReco>Lower[i]) )
        ptBin = i;
    }

    if( ptBin<0 ) continue;

    //fill pt spectrum histogram (needed for stat error):
    h1_ptPhot->Fill( ptPhotReco, eventWeight );
    h2_pt2ndJetRel_vs_ptPhot->Fill( ptPhotReco, pt2ndJetReco/ptPhotReco, eventWeight);
    //h1_pt2ndJet[ptBin]->Fill( pt2ndJetReco, eventWeight );


    //------------ STANDARD CUT ON SECOND JET ------------------------
    
    if( (pt2ndJetReco<0.1*ptPhotReco)||(pt2ndJetReco<10.) ) {

      if( back2back && fabs(etaJetReco-etaPhotReco)<1. ) {
        if(r_RecoPhot!=0.) h2_response_vs_ptPhot->Fill(ptPhotReco, r_RecoPhot, eventWeight); 
        h2_ptPhotMean_vs_pt_cutOn2ndJet->Fill(ptPhotReco, ptPhotReco, eventWeight);
        h1_ptPhot_cutOn2ndJet->Fill( ptPhotReco, eventWeight );
      }

    } 

   
    // --- REQUIRE THAT 2nd JET IS NOT IN PHOTON HEMISPHERE:
    Float_t deltaPhi_phot2nd = phi2ndJetReco - phiPhotReco;
    if( deltaPhi_phot2nd >= TMath::Pi() ) deltaPhi_phot2nd -= 2.*TMath::Pi();
    if( deltaPhi_phot2nd < -TMath::Pi() ) deltaPhi_phot2nd += 2.*TMath::Pi();
   
    if( !NOHEMISPHERE_ ) 
      if( fabs(deltaPhi_phot2nd) < (TMath::Pi()/2.) ) continue;


    h2_ptPhotMean_vs_pt->Fill( ptPhotReco, ptPhotReco, eventWeight );

    //------------ ABSOLUTE SECOND JET PT BINNING --------------------
    Float_t minPt = h1_pt2ndJetRecoMean[ptBin][0]->GetXaxis()->GetXmin();
    Float_t ptStep = h1_pt2ndJetRecoMean[ptBin][0]->GetXaxis()->GetXmax()  - minPt;

    int iGen = (int)floor((pt2ndJetGen-minPt)/ptStep);

    if( (iGen>=0)&&(iGen<5) ) {
      if(r_RecoPhot!=0.) h1_RecoPhot_vs_gen[ptBin][iGen]->Fill(r_RecoPhot, eventWeight);
      if(r_RecoGen!=0.)  h1_RecoGen_vs_gen[ptBin][iGen]->Fill(r_RecoGen, eventWeight);
      if(r_GenPart!=0.)  h1_GenPart_vs_gen[ptBin][iGen]->Fill(r_GenPart, eventWeight);
      if(r_PartPhot!=0.) h1_PartPhot_vs_gen[ptBin][iGen]->Fill(r_PartPhot, eventWeight);
      if(r_GenPhot!=0.)  h1_GenPhot_vs_gen[ptBin][iGen]->Fill(r_GenPhot, eventWeight);
      h1_pt2ndJetGenMean[ptBin][iGen]->Fill(pt2ndJetGen, eventWeight);
      //h1_pt2ndJetGenMean[ptBin][iGen]->Fill(ptSecondaryJetsGen, eventWeight);
    } 


    int iReco = (int)floor((pt2ndJetReco-minPt)/ptStep);

    if( (iReco>=0)&&(iReco<5) ) {
      if(r_RecoPhot!=0.) h1_RecoPhot_vs_reco[ptBin][iReco]->Fill(r_RecoPhot, eventWeight);
      if(r_RecoGen!=0.)  h1_RecoGen_vs_reco[ptBin][iReco]->Fill(r_RecoGen, eventWeight);
      if(r_GenPart!=0.)  h1_GenPart_vs_reco[ptBin][iReco]->Fill(r_GenPart, eventWeight);
      if(r_PartPhot!=0.) h1_PartPhot_vs_reco[ptBin][iReco]->Fill(r_PartPhot, eventWeight);
      if(r_GenPhot!=0.)  h1_GenPhot_vs_reco[ptBin][iReco]->Fill(r_GenPhot, eventWeight);
      h1_pt2ndJetRecoMean[ptBin][iReco]->Fill(pt2ndJetReco, eventWeight);
      //h1_pt2ndJetRecoMean[ptBin][iReco]->Fill(ptSecondaryJetsReco, eventWeight);
    } 



    //------------ RELATIVE SECOND JET PT BINNING --------------------

    Double_t pt2ndJetGenRel = 100.*pt2ndJetGen/ptPhotReco; //in percentage
    Float_t minPerc = h1_pt2ndJetRecoRelMean[ptBin][0]->GetXaxis()->GetXmin();
    Float_t percStep = h1_pt2ndJetRecoRelMean[ptBin][0]->GetXaxis()->GetXmax()  - minPerc;

    int iGenRel = (int)floor((pt2ndJetGenRel-minPerc)/percStep);

    if( DEBUG_ ) {
      std::cout << std::endl << "---------" << std::endl;
      std::cout << "ptSecondJet: " << pt2ndJetGen << "\tptPhot: " << ptPhotReco << "\tpt2ndJetRel: " << pt2ndJetGenRel << std::endl;
      std::cout << "pt2ndJetGenRel: " << pt2ndJetGenRel << " minPerc: " << minPerc << " percStep: " << percStep << std::endl;
      std::cout << "iGenRel: " << iGenRel << std::endl;
      std::cout << std::endl << "---------" << std::endl;
    }

    if( (iGenRel>=0)&&(iGenRel<5) ) {
      if(r_RecoPhot!=0.) h1_RecoPhot_vs_genRel[ptBin][iGenRel]->Fill(r_RecoPhot, eventWeight);
      if(r_RecoGen!=0.)  h1_RecoGen_vs_genRel[ptBin][iGenRel]->Fill(r_RecoGen, eventWeight);
      if(r_GenPart!=0.)  h1_GenPart_vs_genRel[ptBin][iGenRel]->Fill(r_GenPart, eventWeight);
      if(r_PartPhot!=0.) h1_PartPhot_vs_genRel[ptBin][iGenRel]->Fill(r_PartPhot, eventWeight);
      if(r_GenPhot!=0.)  h1_GenPhot_vs_genRel[ptBin][iGenRel]->Fill(r_GenPhot, eventWeight);
      h1_pt2ndJetGenRelMean[ptBin][iGenRel]->Fill(pt2ndJetGenRel, eventWeight);
      //std::cout << "Filling : iGenRel: " << iGenRel << "pt2ndJetGenRel: " << pt2ndJetGenRel << std::endl;
      //h1_pt2ndJetGenMean[ptBin][iGenRel]->Fill(ptSecondaryJetsGen, eventWeight);
    } 

    Double_t pt2ndJetRecoRel = 100.*pt2ndJetReco/ptPhotReco; //in percentage
    int iRecoRel = (int)floor((pt2ndJetRecoRel-minPerc)/percStep);

    if( (iRecoRel>=0)&&(iRecoRel<5) ) {
      if(r_RecoPhot!=0.) h1_RecoPhot_vs_recoRel[ptBin][iRecoRel]->Fill(r_RecoPhot, eventWeight);
      if(r_RecoGen!=0.)  h1_RecoGen_vs_recoRel[ptBin][iRecoRel]->Fill(r_RecoGen, eventWeight);
      if(r_GenPart!=0.)  h1_GenPart_vs_recoRel[ptBin][iRecoRel]->Fill(r_GenPart, eventWeight);
      if(r_PartPhot!=0.) h1_PartPhot_vs_recoRel[ptBin][iRecoRel]->Fill(r_PartPhot, eventWeight);
      if(r_GenPhot!=0.)  h1_GenPhot_vs_recoRel[ptBin][iRecoRel]->Fill(r_GenPhot, eventWeight);
      h1_pt2ndJetRecoRelMean[ptBin][iRecoRel]->Fill(pt2ndJetRecoRel, eventWeight);
      //h1_pt2ndJetGenMean[ptBin][iRecoRel]->Fill(ptSecondaryJetsGen, eventWeight);
    } 



  } //for entries



  //TFile* outFile = TFile::Open("secondJet_yesPtHat.root", "recreate");
  
  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaSecondJet_"+dataset;
  else {
   if(dataset!="") outfileName = "SecondJet_"+dataset;
   else outfileName = "SecondJet";
  }

  if( NOHEMISPHERE_ ) suffix = suffix + "_NOHEMISPHERE";
  if( MCASSOC_ ) suffix = suffix + "_MCassoc";
  if( PARTONTYPE!="ALL" ) suffix =  suffix + "_" + PARTONTYPE;
  outfileName = outfileName + suffix;
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h2_ptPhotMean_vs_pt->Write();
  h1_ptPhot->Write();
  h1_ptPhot_cutOn2ndJet->Write();
  h2_pt2ndJetRel_vs_ptPhot->Write();

  for( int i_pt=0; i_pt<(NBINS_PT-1); ++i_pt) {

    char dirName[40];
    sprintf(dirName, "ptBin_%d_%d", (int)Lower[i_pt], (int)Lower[i_pt+1]);
    outFile->mkdir(dirName);
    outFile->cd(dirName);
    //h1_pt2ndJet[i_pt]->Write();

    for( int i=0; i<nPoints; ++i) {
      
      TF1* gaussian = new TF1("gaussian", "gaus");

      Float_t nSigma = 2.5;
   
      fitTools::fitProjection(h1_RecoPhot_vs_gen[i_pt][i], gaussian, nSigma);
      h1_RecoPhot_vs_gen[i_pt][i]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_gen[i_pt][i], gaussian, nSigma);
      h1_RecoGen_vs_gen[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPart_vs_gen[i_pt][i], gaussian, nSigma);
      h1_GenPart_vs_gen[i_pt][i]->Write();
      fitTools::fitProjection(h1_PartPhot_vs_gen[i_pt][i], gaussian, nSigma);
      h1_PartPhot_vs_gen[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_gen[i_pt][i], gaussian, nSigma);
      h1_GenPhot_vs_gen[i_pt][i]->Write();
      h1_pt2ndJetGenMean[i_pt][i]->Write();

      fitTools::fitProjection(h1_RecoPhot_vs_reco[i_pt][i], gaussian, nSigma);
      h1_RecoPhot_vs_reco[i_pt][i]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_reco[i_pt][i], gaussian, nSigma);
      h1_RecoGen_vs_reco[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPart_vs_reco[i_pt][i], gaussian, nSigma);
      h1_GenPart_vs_reco[i_pt][i]->Write();
      fitTools::fitProjection(h1_PartPhot_vs_reco[i_pt][i], gaussian, nSigma);
      h1_PartPhot_vs_reco[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_reco[i_pt][i], gaussian, nSigma);
      h1_GenPhot_vs_reco[i_pt][i]->Write();
      h1_pt2ndJetRecoMean[i_pt][i]->Write();

      fitTools::fitProjection(h1_RecoPhot_vs_genRel[i_pt][i], gaussian, nSigma);
      h1_RecoPhot_vs_genRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_genRel[i_pt][i], gaussian, nSigma);
      h1_RecoGen_vs_genRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPart_vs_genRel[i_pt][i], gaussian, nSigma);
      h1_GenPart_vs_genRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_PartPhot_vs_genRel[i_pt][i], gaussian, nSigma);
      h1_PartPhot_vs_genRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_genRel[i_pt][i], gaussian, nSigma);
      h1_GenPhot_vs_genRel[i_pt][i]->Write();
      h1_pt2ndJetGenRelMean[i_pt][i]->Write();

      fitTools::fitProjection(h1_RecoPhot_vs_recoRel[i_pt][i], gaussian, nSigma);
      h1_RecoPhot_vs_recoRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_RecoGen_vs_recoRel[i_pt][i], gaussian, nSigma);
      h1_RecoGen_vs_recoRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPart_vs_recoRel[i_pt][i], gaussian, nSigma);
      h1_GenPart_vs_recoRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_PartPhot_vs_recoRel[i_pt][i], gaussian, nSigma);
      h1_PartPhot_vs_recoRel[i_pt][i]->Write();
      fitTools::fitProjection(h1_GenPhot_vs_recoRel[i_pt][i], gaussian, nSigma);
      h1_GenPhot_vs_recoRel[i_pt][i]->Write();
      h1_pt2ndJetRecoRelMean[i_pt][i]->Write();

    } //for n points

  } //for pt bins


  outFile->Close();

  std::string distribFileName = "SecondJet_barrel"+suffix+"_RECO";

  fitTools::fitDistribution_TGraph(h2_response_vs_ptPhot, h2_ptPhotMean_vs_pt_cutOn2ndJet, "response", "pt", "barrel", "", algoType, outfileName, distribFileName, 0.95);


  //clean up mess:

  deleteHistoVector(h1_RecoPhot_vs_gen, nPoints);
  deleteHistoVector(h1_RecoGen_vs_gen, nPoints);
  deleteHistoVector(h1_GenPart_vs_gen, nPoints);
  deleteHistoVector(h1_PartPhot_vs_gen, nPoints);
  deleteHistoVector(h1_GenPhot_vs_gen, nPoints);
  deleteHistoVector(h1_pt2ndJetGenMean, nPoints);

  deleteHistoVector(h1_RecoPhot_vs_reco, nPoints);
  deleteHistoVector(h1_RecoGen_vs_reco, nPoints);
  deleteHistoVector(h1_GenPart_vs_reco, nPoints);
  deleteHistoVector(h1_PartPhot_vs_reco, nPoints);
  deleteHistoVector(h1_GenPhot_vs_reco, nPoints);
  deleteHistoVector(h1_pt2ndJetRecoMean, nPoints);

  deleteHistoVector(h1_RecoPhot_vs_genRel, nPoints);
  deleteHistoVector(h1_RecoGen_vs_genRel, nPoints);
  deleteHistoVector(h1_GenPart_vs_genRel, nPoints);
  deleteHistoVector(h1_PartPhot_vs_genRel, nPoints);
  deleteHistoVector(h1_GenPhot_vs_genRel, nPoints);
  deleteHistoVector(h1_pt2ndJetGenRelMean, nPoints);

  deleteHistoVector(h1_RecoPhot_vs_recoRel, nPoints);
  deleteHistoVector(h1_RecoGen_vs_recoRel, nPoints);
  deleteHistoVector(h1_GenPart_vs_recoRel, nPoints);
  deleteHistoVector(h1_PartPhot_vs_recoRel, nPoints);
  deleteHistoVector(h1_GenPhot_vs_recoRel, nPoints);
  deleteHistoVector(h1_pt2ndJetRecoRelMean, nPoints);


//delete response_vs_ptPhot;
//response_vs_ptPhot = 0;
//delete ptPhotMean_vs_pt_cutOn2ndJet;
//ptPhotMean_vs_pt_cutOn2ndJet = 0;
//delete ptPhotMean_vs_pt;
//ptPhotMean_vs_pt = 0;
  delete h1_ptPhot;
  h1_ptPhot = 0;
  delete h1_ptPhot_cutOn2ndJet;
  h1_ptPhot_cutOn2ndJet = 0;

}
  


std::vector< std::vector< TH1D* > > allocateHistoVector(const std::string& name, const std::string& abs_rel, Int_t nPoints, bool is_pt) {


  std::vector< std::vector< TH1D* > > returnVector;

  Double_t Lower[NBINS_PT];
  fitTools::getPtBins(NBINS_PT, Lower);

  for( int i_ptBin=0; i_ptBin<NBINS_PT; ++i_ptBin ) {

    std::vector< TH1D* > histoVector;
    returnVector.push_back( histoVector );

    Float_t minPt;
    Float_t ptStep;
    if( abs_rel=="abs" ) {
      minPt = 5.;
      ptStep = 5.;
    } else { //minpt and step are in percent:

      if( RECOTYPE_=="pf" ) {
        if( Lower[i_ptBin]<=80.) {
          minPt = 8.;
          ptStep = 3.;
        } else if( Lower[i_ptBin]<=350.) {
          minPt = 5.;
          ptStep = 3.;
        } else {
          minPt = 2.;
          ptStep = 2.;
        }
      } else { //is calo
        if( Lower[i_ptBin]<=80.) {
          minPt = 8.;
          ptStep = 1.5;
        } else if( Lower[i_ptBin]<=350.) {
          minPt = 6.;
          ptStep = 2.;
        } else {
          minPt = 2.;
          ptStep = 2.;
        }
      } //if calo-pf
    }
 
    double iPt = minPt;
    Int_t nBins = 50;

    for( int i=0; i<nPoints; ++i) {

      Double_t xMin = (is_pt) ? iPt : 0.;
      Double_t xMax = (is_pt) ? (iPt+ptStep) : 2.;

      TH1D* tmp;

      std::string prefix = (is_pt) ? "" : "r_";
      char histName[100];
      sprintf(histName, "%s%s_%d_%d", prefix.c_str(), name.c_str(), i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, xMin, xMax); 
      tmp->Sumw2();
      returnVector[i_ptBin].push_back(tmp);

      iPt+=ptStep;

    } //for nPoints

  } //for ptBins


  return returnVector;

}
  

void deleteHistoVector(std::vector< std::vector< TH1D* > > histoVector, int nPoints) {


  for( int i_ptBin=0; i_ptBin<NBINS_PT; ++i_ptBin ) {

    for( int i=0; i<nPoints; ++i) {

      delete histoVector[i_ptBin][i];
      histoVector[i_ptBin][i] = 0;

    }
  }

}
