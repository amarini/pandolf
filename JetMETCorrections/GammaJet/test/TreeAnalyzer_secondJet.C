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
bool MCASSOC_ = true;
bool NOEVENTWEIGHTS_ = false;
bool NOHEMISPHERE_  = true;
std::string PARTONTYPE = "ALL"; //this can be: "ALL", "GLUON", "QUARK"


void secondJet(const std::string& dataset, std::string algoType, std::string partonType="ALL") {

  PARTONTYPE = partonType;

  std::string infile_dataset = dataset;

  TChain* tree = new TChain("jetTree");
  
  std::string infileName;

  std::string prefix = "PhotonJet_2ndLevelTree_";
  std::string suffix = "_"+algoType;

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

  } else {

    if( dataset=="" )
      infileName = "PhotonJet_2ndLevelTree.root/jetTree";
    else
      infileName = "PhotonJet_2ndLevelTree_"+infile_dataset+".root/jetTree";

    std::cout << "Opening file: " << infileName << std::endl;
    tree->Add(infileName.c_str());

  }

  if( NOEVENTWEIGHTS_ ) suffix = "_NOWEIGHTS" + suffix;


  bool draw=true;


  TH1D* h1_deltaR_phot_2ndJet = new TH1D("deltaR_phot_2ndJet", "", 100, 0., 1.);


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
  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_gen;
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_gen;
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_gen;
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_gen;
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_gen;
  std::vector< std::vector< TH1D* > > h1_pt2ndJetGenMean;

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_reco;
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_reco;
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_reco;
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_reco;
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_reco;
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoMean;

  //relative 2nd jet pt binning histograms
  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_genRel;
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_genRel;
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_genRel;
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_genRel;
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_genRel;
  std::vector< std::vector< TH1D* > > h1_pt2ndJetGenRelMean;

  std::vector< std::vector< TH1D* > > h1_RecoPhot_vs_recoRel;
  std::vector< std::vector< TH1D* > > h1_RecoGen_vs_recoRel;
  std::vector< std::vector< TH1D* > > h1_GenPart_vs_recoRel;
  std::vector< std::vector< TH1D* > > h1_PartPhot_vs_recoRel;
  std::vector< std::vector< TH1D* > > h1_GenPhot_vs_recoRel;
  std::vector< std::vector< TH1D* > > h1_pt2ndJetRecoRelMean;

  //pt bins:
  //int nBins_pt = 24;
  int nBins_pt = 16;
  int i=0;

  Double_t Lower[nBins_pt];
  Lower[i++] = 20.;
  Lower[i++] = 42.;
  Lower[i++] = 62.;
  Lower[i++] = 87.;
  Lower[i++] = 108.;
  Lower[i++] = 140.;
  Lower[i++] = 183.;
  Lower[i++] = 239.;
  Lower[i++] = 312.;
  Lower[i++] = 407.;
  Lower[i++] = 530.;
  Lower[i++] = 691.;
  Lower[i++] = 902.;
  Lower[i++] = 1176.;
  Lower[i++] = 1533.;
  Lower[i++] = 2000.;

  //Lower[i++] = 97.;
  //Lower[i++] = 108.;
  //Lower[i++] = 121.;
  //Lower[i++] = 135.;
  //Lower[i++] = 152.;
  //Lower[i++] = 172.;
  //Lower[i++] = 197.;
  //Lower[i++] = 231.;
  //Lower[i++] = 283.;
  //Lower[i++] = 397.;
  //Lower[i++] = 600.;
  //Lower[i++] = 1000.;
  //Lower[i++] = 2000.;

  //Lower[i++] = 20.;
  //Lower[i++] = 25.;
  //Lower[i++] = 30.;
  //Lower[i++] = 36.;
  //Lower[i++] = 42.;
  //Lower[i++] = 48.;
  //Lower[i++] = 55.;
  //Lower[i++] = 62.;
  //Lower[i++] = 70.;
  //Lower[i++] = 78.;
  //Lower[i++] = 87.;
  //Lower[i++] = 97.;
  //Lower[i++] = 108.;
  //Lower[i++] = 121.;
  //Lower[i++] = 135.;
  //Lower[i++] = 152.;
  //Lower[i++] = 172.;
  //Lower[i++] = 197.;
  //Lower[i++] = 231.;
  //Lower[i++] = 283.;
  //Lower[i++] = 397.;
  //Lower[i++] = 600.;
  //Lower[i++] = 1000.;
  //Lower[i++] = 2000.;
  
  TH2D* h2_response_vs_ptPhot = new TH2D("response_vs_ptPhot", "", nBins_pt-1, Lower, 120, 0., 2.);
  h2_response_vs_ptPhot->Sumw2();

  TH2D* h2_ptPhotMean_vs_pt_cutOn2ndJet = new TH2D("ptPhotMean_vs_pt_cutOn2ndJet", "", nBins_pt-1, Lower, 2200, 0., 2200.);
  h2_ptPhotMean_vs_pt_cutOn2ndJet->Sumw2();


  //ptPhot TH2 for x-axis points on vs. pt spectra:
  TH2D* h2_ptPhotMean_vs_pt= new TH2D("ptPhotMean_vs_pt", "", nBins_pt-1, Lower, 2200, 0., 2200.);
  h2_ptPhotMean_vs_pt->Sumw2();


  Float_t minPt = 5.;
  Float_t ptStep = 5.;


  for( int i_ptBin=0; i_ptBin<nBins_pt; ++i_ptBin ) {

    std::vector< TH1D* > RecoPhot_vs_gen;
    h1_RecoPhot_vs_gen.push_back( RecoPhot_vs_gen );
    std::vector< TH1D* > RecoGen_vs_gen;
    h1_RecoGen_vs_gen.push_back( RecoGen_vs_gen );
    std::vector< TH1D* > GenPart_vs_gen;
    h1_GenPart_vs_gen.push_back( GenPart_vs_gen );
    std::vector< TH1D* > PartPhot_vs_gen;
    h1_PartPhot_vs_gen.push_back( PartPhot_vs_gen );
    std::vector< TH1D* > GenPhot_vs_gen;
    h1_GenPhot_vs_gen.push_back( GenPhot_vs_gen );
    std::vector< TH1D* > pt2ndJetGenMean;
    h1_pt2ndJetGenMean.push_back( pt2ndJetGenMean );

    std::vector< TH1D* > RecoPhot_vs_reco;
    h1_RecoPhot_vs_reco.push_back( RecoPhot_vs_reco );
    std::vector< TH1D* > RecoGen_vs_reco;
    h1_RecoGen_vs_reco.push_back( RecoGen_vs_reco );
    std::vector< TH1D* > GenPart_vs_reco;
    h1_GenPart_vs_reco.push_back( GenPart_vs_reco );
    std::vector< TH1D* > PartPhot_vs_reco;
    h1_PartPhot_vs_reco.push_back( PartPhot_vs_reco );
    std::vector< TH1D* > GenPhot_vs_reco;
    h1_GenPhot_vs_reco.push_back( GenPhot_vs_reco );
    std::vector< TH1D* > pt2ndJetRecoMean;
    h1_pt2ndJetRecoMean.push_back( pt2ndJetRecoMean );

    std::vector< TH1D* > RecoPhot_vs_genRel;
    h1_RecoPhot_vs_genRel.push_back( RecoPhot_vs_genRel );
    std::vector< TH1D* > RecoGen_vs_genRel;
    h1_RecoGen_vs_genRel.push_back( RecoGen_vs_genRel );
    std::vector< TH1D* > GenPart_vs_genRel;
    h1_GenPart_vs_genRel.push_back( GenPart_vs_genRel );
    std::vector< TH1D* > PartPhot_vs_genRel;
    h1_PartPhot_vs_genRel.push_back( PartPhot_vs_genRel );
    std::vector< TH1D* > GenPhot_vs_genRel;
    h1_GenPhot_vs_genRel.push_back( GenPhot_vs_genRel );
    std::vector< TH1D* > pt2ndJetGenRelMean;
    h1_pt2ndJetGenRelMean.push_back( pt2ndJetGenRelMean );

    std::vector< TH1D* > RecoPhot_vs_recoRel;
    h1_RecoPhot_vs_recoRel.push_back( RecoPhot_vs_recoRel );
    std::vector< TH1D* > RecoGen_vs_recoRel;
    h1_RecoGen_vs_recoRel.push_back( RecoGen_vs_recoRel );
    std::vector< TH1D* > GenPart_vs_recoRel;
    h1_GenPart_vs_recoRel.push_back( GenPart_vs_recoRel );
    std::vector< TH1D* > PartPhot_vs_recoRel;
    h1_PartPhot_vs_recoRel.push_back( PartPhot_vs_recoRel );
    std::vector< TH1D* > GenPhot_vs_recoRel;
    h1_GenPhot_vs_recoRel.push_back( GenPhot_vs_recoRel );
    std::vector< TH1D* > pt2ndJetRecoRelMean;
    h1_pt2ndJetRecoRelMean.push_back( pt2ndJetRecoRelMean );

    int iPt = minPt;
    char histName[30];
    Int_t nBins = 50;
 
    //absolute second jet pt binning histograms:
    //absolute binning:  5-10GeV 10-15GeV 15-20GeV 20-25GeV 25-30GeV

    for( int i=0; i<nPoints; ++i) {
      TH1D* tmp;

      sprintf(histName, "r_RecoPhot_vs_Gen_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoPhot_vs_gen[i_ptBin].push_back(tmp);
      sprintf(histName, "r_RecoGen_vs_Gen_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoGen_vs_gen[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPart_vs_Gen_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPart_vs_gen[i_ptBin].push_back(tmp);
      sprintf(histName, "r_PartPhot_vs_Gen_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_PartPhot_vs_gen[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPhot_vs_Gen_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPhot_vs_gen[i_ptBin].push_back(tmp);
      sprintf(histName, "pt2ndJetGenMean_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, iPt, iPt+ptStep);
      tmp->Sumw2();
      h1_pt2ndJetGenMean[i_ptBin].push_back(tmp);

      sprintf(histName, "r_RecoPhot_vs_Reco_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoPhot_vs_reco[i_ptBin].push_back(tmp);
      sprintf(histName, "r_RecoGen_vs_Reco_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoGen_vs_reco[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPart_vs_Reco_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPart_vs_reco[i_ptBin].push_back(tmp);
      sprintf(histName, "r_PartPhot_vs_Reco_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_PartPhot_vs_reco[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPhot_vs_Reco_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPhot_vs_reco[i_ptBin].push_back(tmp);
      sprintf(histName, "pt2ndJetRecoMean_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, iPt, iPt+ptStep);
      tmp->Sumw2();
      h1_pt2ndJetRecoMean[i_ptBin].push_back(tmp);

      iPt+=ptStep;

    } //for nPoints

    //relative second jet pt binning histograms:
    //relative binning:  5-8%, 8-11%, 11-14%, 14-17%, 17-20% of pt_photon

    Float_t minPerc;
    Float_t percStep;

    //variable relative binning range in different pt regions
    //if( Lower[i_ptBin]<=60.) {
    //  minPerc = 15.;
    //  percStep = 3.;
    if( Lower[i_ptBin]<=80.) {
      minPerc = 8.;
      percStep = 3.;
    } else if( Lower[i_ptBin]<=350.) {
      minPerc = 5.;
      percStep = 3.;
    } else {
      minPerc = 2.;
      percStep = 2.;
    }


    int iPerc = minPerc;


    for( int i=0; i<nPoints; ++i) {
      TH1D* tmp;

      sprintf(histName, "r_RecoPhot_vs_GenRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoPhot_vs_genRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_RecoGen_vs_GenRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoGen_vs_genRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPart_vs_GenRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPart_vs_genRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_PartPhot_vs_GenRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_PartPhot_vs_genRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPhot_vs_GenRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPhot_vs_genRel[i_ptBin].push_back(tmp);
      sprintf(histName, "pt2ndJetGenRelMean_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, iPerc, iPerc+percStep);
      tmp->Sumw2();
      h1_pt2ndJetGenRelMean[i_ptBin].push_back(tmp);

      sprintf(histName, "r_RecoPhot_vs_RecoRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoPhot_vs_recoRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_RecoGen_vs_RecoRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_RecoGen_vs_recoRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPart_vs_RecoRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPart_vs_recoRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_PartPhot_vs_RecoRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_PartPhot_vs_recoRel[i_ptBin].push_back(tmp);
      sprintf(histName, "r_GenPhot_vs_RecoRel_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, 0., 2.); 
      tmp->Sumw2();
      h1_GenPhot_vs_recoRel[i_ptBin].push_back(tmp);
      sprintf(histName, "pt2ndJetRecoRelMean_%d_%d", i_ptBin, i);
      tmp = new TH1D(histName, "", nBins, iPerc, iPerc+percStep);
      tmp->Sumw2();
      h1_pt2ndJetRecoRelMean[i_ptBin].push_back(tmp);

      iPerc+=percStep;

    } //for nPoints

  } //for ptBins




  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;
    if( NOEVENTWEIGHTS_ )   eventWeight = 1.;

    //both photon and jet in barrel (maybe jet requirement could induce bias):
    if( fabs(etaPhotReco)>1.3 ) continue;
    if( fabs(etaJetReco)>1.4 ) continue;

    
    //leading jet and photon back2back in transverse plane
    bool back2back = true;
    Float_t deltaPhi_jet = phiPhotReco - phiJetReco;
    Float_t pi = 3.14159;
    if( fabs(deltaPhi_jet) < (pi - 0.2) ) back2back = false;
    //if( !back2back ) continue;  //09.11.09 deltaPhi cut is actually useless: proved it.


    if( MCASSOC_ )
      if( !matchedToMC ) continue;
    else
      if( !passedPhotonID_medium ) continue;
    
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

    for( int i=0; i<nBins_pt; ++i) {
      if( (ptPhotReco<Lower[i+1])&&(ptPhotReco>Lower[i]) )
        ptBin = i;
    }

    if( ptBin<0 ) continue;


    //------------ STANDARD CUT ON SECOND JET ------------------------
    
    if( (pt2ndJetReco<0.1*ptPhotReco)||(pt2ndJetReco<10.) ) {
      if(r_RecoPhot!=0.) h2_response_vs_ptPhot->Fill(ptPhotReco, r_RecoPhot, eventWeight); 
      h2_ptPhotMean_vs_pt_cutOn2ndJet->Fill(ptPhotReco, ptPhotReco, eventWeight);
    } 

   
    // --- REQUIRE THAT 2nd JET IS NOT IN PHOTON HEMISPHERE:
    Float_t deltaPhi_phot2nd = phi2ndJetReco - phiPhotReco;
    if( deltaPhi_phot2nd >= TMath::Pi() ) deltaPhi_phot2nd -= 2.*TMath::Pi();
    if( deltaPhi_phot2nd < -TMath::Pi() ) deltaPhi_phot2nd += 2.*TMath::Pi();
   
    if( !NOHEMISPHERE_ ) 
      if( fabs(deltaPhi_phot2nd) < (TMath::Pi()/2.) ) continue;


    h2_ptPhotMean_vs_pt->Fill( ptPhotReco, ptPhotReco, eventWeight );

    //------------ ABSOLUTE SECOND JET PT BINNING --------------------

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

  if( PARTONTYPE!="ALL" ) suffix = "_" + PARTONTYPE + suffix;
  if( MCASSOC_ ) suffix = "_MCassoc" + suffix;
  if( NOHEMISPHERE_ ) suffix = "_NOHEMISPHERE" + suffix;
  outfileName = outfileName + suffix;
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h2_ptPhotMean_vs_pt->Write();

  for( int i_pt=0; i_pt<(nBins_pt-1); ++i_pt) {

    char dirName[40];
    sprintf(dirName, "ptBin_%d_%d", (int)Lower[i_pt], (int)Lower[i_pt+1]);
    outFile->mkdir(dirName);
    outFile->cd(dirName);

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



}
  

