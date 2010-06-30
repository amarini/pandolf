#include <TH1F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>
#include <iostream>
#include <vector>
#include <cmath>




bool PFJetID( Float_t Rch, Float_t Rnh, Float_t Rgamma, Int_t N, Int_t Nch);



void finalize(const std::string& dataset, std::string algoType, Float_t pt_thresh=5., std::string raw_corr="raw", bool use_deltaPhi=true) {
  
  std::string infileName = "MinBias_2ndLevelTree_" + dataset + "_" + algoType + ".root";

  TFile* infile = TFile::Open(infileName.c_str());
  TTree* tree = (TTree*)infile->Get("jetTree");
  std::cout << "Added " + infileName + ". tree has " << tree->GetEntries() << " entries." << std::endl;
  
  std::cout << "Done adding." << std::endl;


  Float_t trackRatio;
  tree->SetBranchAddress("trackRatio", &trackRatio);

  Float_t ePhotReco;
  tree->SetBranchAddress("ePhotReco", &ePhotReco);
  Float_t ptPhotReco;
  tree->SetBranchAddress("ptPhotReco", &ptPhotReco);
  Float_t etaPhotReco;
  tree->SetBranchAddress("etaPhotReco", &etaPhotReco);
  Float_t phiPhotReco;
  tree->SetBranchAddress("phiPhotReco", &phiPhotReco);

  Float_t eJetReco;
  tree->SetBranchAddress("eJetReco", &eJetReco);
  Float_t ptJetReco;
  tree->SetBranchAddress("ptJetReco", &ptJetReco);
  Float_t ptCorrJetReco;
  tree->SetBranchAddress("ptCorrJetReco", &ptCorrJetReco);
  Float_t etaJetReco;
  tree->SetBranchAddress("etaJetReco", &etaJetReco);
  Float_t phiJetReco;
  tree->SetBranchAddress("phiJetReco", &phiJetReco);

  Float_t eTracksReco;
  tree->SetBranchAddress("eTracksReco", &eTracksReco);
  Float_t ePhotonsReco;
  tree->SetBranchAddress("ePhotonsReco", &ePhotonsReco);
  Float_t eNeutralHadronsReco;
  tree->SetBranchAddress("eNeutralHadronsReco", &eNeutralHadronsReco);
  Float_t eMuonsReco;
  tree->SetBranchAddress("eMuonsReco", &eMuonsReco);
  Float_t eElectronsReco;
  tree->SetBranchAddress("eElectronsReco", &eElectronsReco);
  Float_t eHFHadronsReco;
  tree->SetBranchAddress("eHFHadronsReco", &eHFHadronsReco);
  Float_t eHFEMReco;
  tree->SetBranchAddress("eHFEMReco", &eHFEMReco);

  Int_t nTracksReco;
  tree->SetBranchAddress("nTracksReco", &nTracksReco);
  Int_t nPhotonsReco;
  tree->SetBranchAddress("nPhotonsReco", &nPhotonsReco);
  Int_t nNeutralHadronsReco;
  tree->SetBranchAddress("nNeutralHadronsReco", &nNeutralHadronsReco);
  Int_t nElectronsReco;
  tree->SetBranchAddress("nElectronsReco", &nElectronsReco);
  Int_t nMuonsReco;
  tree->SetBranchAddress("nMuonsReco", &nMuonsReco);
  Int_t nHFHadronsReco;
  tree->SetBranchAddress("nHFHadronsReco", &nHFHadronsReco);
  Int_t nHFEMReco;
  tree->SetBranchAddress("nHFEMReco", &nHFEMReco);

  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t ptCorr2ndJetReco;
  tree->SetBranchAddress("ptCorr2ndJetReco", &ptCorr2ndJetReco);
  Float_t e2ndJetReco;
  tree->SetBranchAddress("e2ndJetReco", &e2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);

  Float_t eTracks2ndReco;
  tree->SetBranchAddress("eTracks2ndReco", &eTracks2ndReco);
  Float_t ePhotons2ndReco;
  tree->SetBranchAddress("ePhotons2ndReco", &ePhotons2ndReco);
  Float_t eNeutralHadrons2ndReco;
  tree->SetBranchAddress("eNeutralHadrons2ndReco", &eNeutralHadrons2ndReco);
  Float_t eMuons2ndReco;
  tree->SetBranchAddress("eMuons2ndReco", &eMuons2ndReco);
  Float_t eElectrons2ndReco;
  tree->SetBranchAddress("eElectrons2ndReco", &eElectrons2ndReco);
  Float_t eHFHadrons2ndReco;
  tree->SetBranchAddress("eHFHadrons2ndReco", &eHFHadrons2ndReco);
  Float_t eHFEM2ndReco;
  tree->SetBranchAddress("eHFEM2ndReco", &eHFEM2ndReco);

  Int_t nTracks2ndReco;
  tree->SetBranchAddress("nTracks2ndReco", &nTracks2ndReco);
  Int_t nPhotons2ndReco;
  tree->SetBranchAddress("nPhotons2ndReco", &nPhotons2ndReco);
  Int_t nNeutralHadrons2ndReco;
  tree->SetBranchAddress("nNeutralHadrons2ndReco", &nNeutralHadrons2ndReco);
  Int_t nElectrons2ndReco;
  tree->SetBranchAddress("nElectrons2ndReco", &nElectrons2ndReco);
  Int_t nMuons2ndReco;
  tree->SetBranchAddress("nMuons2ndReco", &nMuons2ndReco);
  Int_t nHFHadrons2ndReco;
  tree->SetBranchAddress("nHFHadrons2ndReco", &nHFHadrons2ndReco);
  Int_t nHFEM2ndReco;
  tree->SetBranchAddress("nHFEM2ndReco", &nHFEM2ndReco);

  Int_t nBins_eta = 40;
  Int_t nBins_pt = 30;
  Int_t nBins_R = 30;
  if( pt_thresh>8. ) {
    nBins_eta /= 2;
    nBins_pt /= 2;
    nBins_R /= 2;
  }
  
  TH1F* h1_deltaPhiJet = new TH1F("deltaPhiJet", "", 20, 1., 3.15);
  TH1F* h1_asymmJet = new TH1F("asymmJet", "", 20, -0.8, 0.8);

  TH1F* h1_massJet = new TH1F("massJet", "", 15., 0., 15.);
  TH1F* h1_diJetMass = new TH1F("diJetMass", "", 50., 0., 100.);

  TH1F* h1_etaJet = new TH1F("etaJet", "", nBins_eta, -2.4, 2.4);
  TH1F* h1_phiJet = new TH1F("phiJet", "", nBins_eta, -3.5, 3.5);
  TH1F* h1_ptJet = new TH1F("ptJet", "", 50, 0., 50.);
  TH1F* h1_ptCorrJet = new TH1F("ptCorrJet", "", 50, 0., 50.);
  TH1F* h1_nCandJet = new TH1F("nCandJet", "", 26, -0.5, 25.);

  TH1F* h1_RchJet = new TH1F("RchJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet = new TH1F("RnhJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet = new TH1F("RgammaJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet = new TH1F("RmuJet", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet = new TH1F("ReJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet = new TH1F("RhfhadJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet = new TH1F("RhfemJet", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet = new TH1F("EchJet", "", nBins_pt, 0., 40.);
  TH1F* h1_EnhJet = new TH1F("EnhJet", "", nBins_pt, 0., 20.);
  TH1F* h1_EnhJet_barrel = new TH1F("EnhJet_barrel", "", nBins_pt, 0., 20.);
  TH1F* h1_EnhJet_endcap = new TH1F("EnhJet_endcap", "", nBins_pt, 0., 20.);
  TH1F* h1_EgammaJet = new TH1F("EgammaJet", "", nBins_pt, 0., 20.);
  TH1F* h1_EeJet = new TH1F("EeJet", "", nBins_pt, 0., 40.);
  TH1F* h1_EmuJet = new TH1F("EmuJet", "", nBins_pt, 0., 40.);
  TH1F* h1_EhfhadJet = new TH1F("EhfhadJet", "", nBins_pt, 0., 40.);
  TH1F* h1_EhfemJet = new TH1F("EhfemJet", "", nBins_pt, 0., 40.);

  TH1F* h1_EphotAveJet = new TH1F("EphotAveJet", "", 40, 0., 10.);

  TH1F* h1_NchJet = new TH1F("NchJet", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet = new TH1F("NnhJet", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet = new TH1F("NgammaJet", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet = new TH1F("NeJet", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet = new TH1F("NmuJet", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet = new TH1F("NhfhadJet", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet = new TH1F("NhfemJet", "", 21, -0.5, 20.5);

  TProfile* hp_pt_vs_eta = new TProfile("pt_vs_eta", "pt vs. eta", nBins_eta, -2.4, 2.4, 0., 300.);
  TProfile* hp_ptCorr_vs_eta = new TProfile("ptCorr_vs_eta", "ptCorr vs. eta", nBins_eta, -2.4, 2.4, 0., 300.);

  TProfile* hp_Rch_vs_eta = new TProfile("Rch_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Rnh_vs_eta = new TProfile("Rnh_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Rgamma_vs_eta = new TProfile("Rgamma_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Rmu_vs_eta = new TProfile("Rmu_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Re_vs_eta = new TProfile("Re_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Rhfhad_vs_eta = new TProfile("Rhfhad_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Rhfem_vs_eta = new TProfile("Rhfem_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);

  TProfile* hp_Ech_vs_eta = new TProfile("Ech_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Enh_vs_eta = new TProfile("Enh_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Egamma_vs_eta = new TProfile("Egamma_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Ee_vs_eta = new TProfile("Ee_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Emu_vs_eta = new TProfile("Emu_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Ehfhad_vs_eta = new TProfile("Ehfhad_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Ehfem_vs_eta = new TProfile("Ehfem_vs_eta", "R_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);

  TProfile* hp_Nch_vs_eta = new TProfile("Nch_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Nnh_vs_eta = new TProfile("Nnh_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Ngamma_vs_eta = new TProfile("Ngamma_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Ne_vs_eta = new TProfile("Ne_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Nmu_vs_eta = new TProfile("Nmu_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Nhfhad_vs_eta = new TProfile("Nhfhad_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);
  TProfile* hp_Nhfem_vs_eta = new TProfile("Nhfem_vs_eta", "N_{ch} vs. eta", nBins_eta, -2.4, 2.4, 0., 100.);

  TProfile* hp_Rch_vs_pt = new TProfile("Rch_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Rnh_vs_pt = new TProfile("Rnh_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Rgamma_vs_pt = new TProfile("Rgamma_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Rmu_vs_pt = new TProfile("Rmu_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Re_vs_pt = new TProfile("Re_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Rhfhad_vs_pt = new TProfile("Rhfhad_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Rhfem_vs_pt = new TProfile("Rhfem_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);

  TProfile* hp_Ech_vs_pt = new TProfile("Ech_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Enh_vs_pt = new TProfile("Enh_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Egamma_vs_pt = new TProfile("Egamma_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Ee_vs_pt = new TProfile("Ee_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Emu_vs_pt = new TProfile("Emu_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Ehfhad_vs_pt = new TProfile("Ehfhad_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Ehfem_vs_pt = new TProfile("Ehfem_vs_pt", "R_{ch} vs. pt", 7, 7., 14., 0., 100.);

  TProfile* hp_Nch_vs_pt = new TProfile("Nch_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Nnh_vs_pt = new TProfile("Nnh_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Ngamma_vs_pt = new TProfile("Ngamma_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Ne_vs_pt = new TProfile("Ne_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Nmu_vs_pt = new TProfile("Nmu_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Nhfhad_vs_pt = new TProfile("Nhfhad_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);
  TProfile* hp_Nhfem_vs_pt = new TProfile("Nhfem_vs_pt", "N_{ch} vs. pt", 7, 7., 14., 0., 100.);

  TProfile* hp_Rch_vs_ptCorr = new TProfile("Rch_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Rnh_vs_ptCorr = new TProfile("Rnh_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Rgamma_vs_ptCorr = new TProfile("Rgamma_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Rmu_vs_ptCorr = new TProfile("Rmu_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Re_vs_ptCorr = new TProfile("Re_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Rhfhad_vs_ptCorr = new TProfile("Rhfhad_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Rhfem_vs_ptCorr = new TProfile("Rhfem_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);

  TProfile* hp_Ech_vs_ptCorr = new TProfile("Ech_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Enh_vs_ptCorr = new TProfile("Enh_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Egamma_vs_ptCorr = new TProfile("Egamma_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Ee_vs_ptCorr = new TProfile("Ee_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Emu_vs_ptCorr = new TProfile("Emu_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Ehfhad_vs_ptCorr = new TProfile("Ehfhad_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Ehfem_vs_ptCorr = new TProfile("Ehfem_vs_ptCorr", "R_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);

  TProfile* hp_Nch_vs_ptCorr = new TProfile("Nch_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Nnh_vs_ptCorr = new TProfile("Nnh_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Ngamma_vs_ptCorr = new TProfile("Ngamma_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Ne_vs_ptCorr = new TProfile("Ne_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Nmu_vs_ptCorr = new TProfile("Nmu_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Nhfhad_vs_ptCorr = new TProfile("Nhfhad_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);
  TProfile* hp_Nhfem_vs_ptCorr = new TProfile("Nhfem_vs_ptCorr", "N_{ch} vs. ptCorr", 7, 7., 14., 0., 100.);

  TProfile* hp_Egamma_vs_Ech = new TProfile("Egamma_vs_Ech", "", nBins_pt, 0., 15., 0., 100.);
  TProfile* hp_Enh_vs_Ech = new TProfile("Enh_vs_Ech", "", nBins_pt, 0., 15., 0., 100.);

  TProfile* hp_Rgamma_vs_Rch = new TProfile("Rgamma_vs_Rch", "", nBins_R, 0., 1.1, 0., 1.1);
  TProfile* hp_Rnh_vs_Rch = new TProfile("Rnh_vs_Rch", "", nBins_R, 0., 1.1, 0., 1.1);
  TProfile* hp_Rnh_vs_Rgamma = new TProfile("Rnh_vs_Rgamma", "", nBins_R, 0., 1.1, 0., 1.1);

  TProfile* hp_Ngamma_vs_Nch = new TProfile("Ngamma_vs_Nch", "", 21, -0.5, 20.5, -0.5, 20.5);
  TProfile* hp_Nnh_vs_Nch = new TProfile("Nnh_vs_Nch", "", 21, -0.5, 20.5, -0.5, 20.5);

  TRandom* rand = new TRandom();

  int nEntries = tree->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( raw_corr=="raw" )
      if( pt2ndJetReco<pt_thresh ) continue;
    else
      if( ptCorr2ndJetReco<pt_thresh ) continue;

    if( !(fabs(etaJetReco)<2.4&&fabs(eta2ndJetReco)<2.4) ) continue;

    //PFJetID:
    Int_t N = nTracksReco+nPhotonsReco+nNeutralHadronsReco+nElectronsReco+nMuonsReco+nHFHadronsReco+nHFEMReco;
    Int_t Nch = nTracksReco+nElectronsReco+nMuonsReco;
    if( !PFJetID( eTracksReco, eNeutralHadronsReco, ePhotonsReco, N, Nch) ) continue;
    
    Int_t N2nd = nTracks2ndReco+nPhotons2ndReco+nNeutralHadrons2ndReco+nElectrons2ndReco+nMuons2ndReco+nHFHadrons2ndReco+nHFEM2ndReco;
    Int_t Nch2nd = nTracks2ndReco+nElectrons2ndReco+nMuons2ndReco;
    if( !PFJetID( eTracks2ndReco, eNeutralHadrons2ndReco, ePhotons2ndReco, N2nd, Nch2nd) ) continue;

    //loose back-to-back
    Float_t deltaPhi =  phiJetReco-phi2ndJetReco;
    if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
    if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
    h1_deltaPhiJet->Fill( deltaPhi );
    if( use_deltaPhi )
      if( fabs(deltaPhi)<(TMath::Pi()-1.) ) continue;

    Float_t coin = rand->Uniform(1.);
    Float_t asymm = (ptJetReco-pt2ndJetReco)/(ptJetReco+pt2ndJetReco);
    if( coin <0.5 ) asymm*=-1.;
    h1_asymmJet->Fill( asymm );
  
    //Float_t ptAve = 0.5*(ptCorrJetReco+ptCorr2ndJetReco);

    h1_ptJet->Fill( ptJetReco );
    h1_ptJet->Fill( pt2ndJetReco );
    h1_ptCorrJet->Fill( ptCorrJetReco );
    h1_ptCorrJet->Fill( ptCorr2ndJetReco );

    hp_pt_vs_eta->Fill( etaJetReco, ptJetReco );
    hp_pt_vs_eta->Fill( eta2ndJetReco, pt2ndJetReco );
    hp_ptCorr_vs_eta->Fill( etaJetReco, ptCorrJetReco );
    hp_ptCorr_vs_eta->Fill( eta2ndJetReco, ptCorr2ndJetReco );
    

    Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
    Float_t pJetReco = ptJetReco/sin(thetaJetReco);
    Float_t pxJetReco = pJetReco*sin(thetaJetReco)*cos(phiJetReco);
    Float_t pyJetReco = pJetReco*sin(thetaJetReco)*sin(phiJetReco);
    Float_t pzJetReco = pJetReco*cos(thetaJetReco);
    Float_t massJet = sqrt( eJetReco*eJetReco - pJetReco*pJetReco );
    h1_massJet->Fill( massJet );

    Float_t theta2ndJetReco = 2.*atan(exp(-eta2ndJetReco));
    Float_t p2ndJetReco = pt2ndJetReco/sin(theta2ndJetReco);
    Float_t px2ndJetReco = p2ndJetReco*sin(theta2ndJetReco)*cos(phi2ndJetReco);
    Float_t py2ndJetReco = p2ndJetReco*sin(theta2ndJetReco)*sin(phi2ndJetReco);
    Float_t pz2ndJetReco = p2ndJetReco*cos(theta2ndJetReco);
    Float_t mass2ndJet = sqrt( e2ndJetReco*e2ndJetReco - p2ndJetReco*p2ndJetReco );
    h1_massJet->Fill( mass2ndJet );

    Float_t eSum = eJetReco+e2ndJetReco;
    Float_t pxSum = pxJetReco+px2ndJetReco;
    Float_t pySum = pyJetReco+py2ndJetReco;
    Float_t pzSum = pzJetReco+pz2ndJetReco;

    Float_t diJetMass = sqrt( eSum*eSum - pxSum*pxSum - pySum*pySum - pzSum*pzSum );
    h1_diJetMass->Fill( diJetMass ); 
  
    h1_etaJet->Fill( etaJetReco );
    h1_etaJet->Fill( eta2ndJetReco );
    h1_phiJet->Fill( phiJetReco );
    h1_phiJet->Fill( phi2ndJetReco );
    h1_nCandJet->Fill( N );
    h1_nCandJet->Fill( N2nd );
  
    h1_RchJet->Fill( eTracksReco );
    h1_RchJet->Fill( eTracks2ndReco );
    h1_RnhJet->Fill( eNeutralHadronsReco );
    h1_RnhJet->Fill( eNeutralHadrons2ndReco );
    h1_RgammaJet->Fill( ePhotonsReco );
    h1_RgammaJet->Fill( ePhotons2ndReco );
    h1_ReJet->Fill( eElectronsReco );
    h1_ReJet->Fill( eElectrons2ndReco );
    h1_RmuJet->Fill( eMuonsReco );
    h1_RmuJet->Fill( eMuons2ndReco );
    h1_RhfhadJet->Fill( eHFHadronsReco );
    h1_RhfhadJet->Fill( eHFHadrons2ndReco );
    h1_RhfemJet->Fill( eHFEMReco );
    h1_RhfemJet->Fill( eHFEM2ndReco );
  
    h1_EchJet->Fill( eTracksReco*eJetReco );
    h1_EchJet->Fill( eTracks2ndReco*e2ndJetReco );
    h1_EnhJet->Fill( eNeutralHadronsReco*eJetReco );
    h1_EnhJet->Fill( eNeutralHadrons2ndReco*e2ndJetReco );
    if( fabs(etaJetReco)<1.4 )
      h1_EnhJet_barrel->Fill( eNeutralHadronsReco*eJetReco );
    if( fabs(eta2ndJetReco)<1.4 )
      h1_EnhJet_barrel->Fill( eNeutralHadrons2ndReco*e2ndJetReco );
    if( fabs(etaJetReco)>1.4 )
      h1_EnhJet_endcap->Fill( eNeutralHadronsReco*eJetReco );
    if( fabs(eta2ndJetReco)>1.4 )
      h1_EnhJet_endcap->Fill( eNeutralHadrons2ndReco*e2ndJetReco );
    h1_EgammaJet->Fill( ePhotonsReco*eJetReco );
    h1_EgammaJet->Fill( ePhotons2ndReco*e2ndJetReco );
    h1_EeJet->Fill( eElectronsReco*eJetReco );
    h1_EeJet->Fill( eElectrons2ndReco*e2ndJetReco );
    h1_EmuJet->Fill( eMuonsReco*eJetReco );
    h1_EmuJet->Fill( eMuons2ndReco*e2ndJetReco );
    h1_EhfhadJet->Fill( eHFHadronsReco*eJetReco );
    h1_EhfhadJet->Fill( eHFHadrons2ndReco*e2ndJetReco );
    h1_EhfemJet->Fill( eHFEMReco*eJetReco );
    h1_EhfemJet->Fill( eHFEM2ndReco*e2ndJetReco );

    if( nPhotonsReco>0 )
      h1_EphotAveJet->Fill( ePhotonsReco*eJetReco/(Float_t)nPhotonsReco );
    if( nPhotons2ndReco>0 )
      h1_EphotAveJet->Fill( ePhotons2ndReco*e2ndJetReco/(Float_t)nPhotons2ndReco );
  
    h1_NchJet->Fill( nTracksReco );
    h1_NchJet->Fill( nTracks2ndReco );
    h1_NnhJet->Fill( nNeutralHadronsReco );
    h1_NnhJet->Fill( nNeutralHadrons2ndReco );
    h1_NgammaJet->Fill( nPhotonsReco );
    h1_NgammaJet->Fill( nPhotons2ndReco );
    h1_NeJet->Fill( nElectronsReco );
    h1_NeJet->Fill( nElectrons2ndReco );
    h1_NmuJet->Fill( nMuonsReco );
    h1_NmuJet->Fill( nMuons2ndReco );
    h1_NhfhadJet->Fill( nHFHadronsReco );
    h1_NhfhadJet->Fill( nHFHadrons2ndReco );
    h1_NhfemJet->Fill( nHFEMReco );
    h1_NhfemJet->Fill( nHFEM2ndReco );
  
    hp_Rch_vs_eta->Fill( etaJetReco, eTracksReco );
    hp_Rch_vs_eta->Fill( eta2ndJetReco, eTracks2ndReco );
    hp_Rnh_vs_eta->Fill( etaJetReco, eNeutralHadronsReco );
    hp_Rnh_vs_eta->Fill( eta2ndJetReco, eNeutralHadrons2ndReco );
    hp_Rgamma_vs_eta->Fill( etaJetReco, ePhotonsReco );
    hp_Rgamma_vs_eta->Fill( eta2ndJetReco, ePhotons2ndReco );
    hp_Rmu_vs_eta->Fill( etaJetReco, eMuonsReco );
    hp_Rmu_vs_eta->Fill( eta2ndJetReco, eMuons2ndReco );
    hp_Re_vs_eta->Fill( etaJetReco, eElectronsReco );
    hp_Re_vs_eta->Fill( eta2ndJetReco, eElectrons2ndReco );
    hp_Rhfhad_vs_eta->Fill( etaJetReco, eHFHadronsReco );
    hp_Rhfhad_vs_eta->Fill( eta2ndJetReco, eHFHadrons2ndReco );
    hp_Rhfem_vs_eta->Fill( etaJetReco, eHFEMReco );
    hp_Rhfem_vs_eta->Fill( eta2ndJetReco, eHFEM2ndReco );
  
    hp_Ech_vs_eta->Fill( etaJetReco, eTracksReco*eJetReco );
    hp_Ech_vs_eta->Fill( eta2ndJetReco, eTracks2ndReco*e2ndJetReco );
    hp_Enh_vs_eta->Fill( etaJetReco, eNeutralHadronsReco*eJetReco );
    hp_Enh_vs_eta->Fill( eta2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco );
    hp_Egamma_vs_eta->Fill( etaJetReco, ePhotonsReco*eJetReco );
    hp_Egamma_vs_eta->Fill( eta2ndJetReco, ePhotons2ndReco*e2ndJetReco );
    hp_Emu_vs_eta->Fill( etaJetReco, eMuonsReco*eJetReco );
    hp_Emu_vs_eta->Fill( eta2ndJetReco, eMuons2ndReco*e2ndJetReco );
    hp_Ee_vs_eta->Fill( etaJetReco, eElectronsReco*eJetReco );
    hp_Ee_vs_eta->Fill( eta2ndJetReco, eElectrons2ndReco*e2ndJetReco );
    hp_Ehfhad_vs_eta->Fill( etaJetReco, eHFHadronsReco*eJetReco );
    hp_Ehfhad_vs_eta->Fill( eta2ndJetReco, eHFHadrons2ndReco*e2ndJetReco );
    hp_Ehfem_vs_eta->Fill( etaJetReco, eHFEMReco*eJetReco );
    hp_Ehfem_vs_eta->Fill( eta2ndJetReco, eHFEM2ndReco*e2ndJetReco );

    hp_Nch_vs_eta->Fill( etaJetReco, nTracksReco );
    hp_Nch_vs_eta->Fill( eta2ndJetReco, nTracks2ndReco );
    hp_Nnh_vs_eta->Fill( etaJetReco, nNeutralHadronsReco );
    hp_Nnh_vs_eta->Fill( eta2ndJetReco, nNeutralHadrons2ndReco );
    hp_Ngamma_vs_eta->Fill( etaJetReco, nPhotonsReco );
    hp_Ngamma_vs_eta->Fill( eta2ndJetReco, nPhotons2ndReco );
    hp_Nmu_vs_eta->Fill( etaJetReco, nMuonsReco );
    hp_Nmu_vs_eta->Fill( eta2ndJetReco, nMuons2ndReco );
    hp_Ne_vs_eta->Fill( etaJetReco, nElectronsReco );
    hp_Ne_vs_eta->Fill( eta2ndJetReco, nElectrons2ndReco );
    hp_Nhfhad_vs_eta->Fill( etaJetReco, nHFHadronsReco );
    hp_Nhfhad_vs_eta->Fill( eta2ndJetReco, nHFHadrons2ndReco );
    hp_Nhfem_vs_eta->Fill( etaJetReco, nHFEMReco );
    hp_Nhfem_vs_eta->Fill( eta2ndJetReco, nHFEM2ndReco );

    hp_Rch_vs_pt->Fill( ptJetReco, eTracksReco );
    hp_Rch_vs_pt->Fill( pt2ndJetReco, eTracks2ndReco );
    hp_Rnh_vs_pt->Fill( ptJetReco, eNeutralHadronsReco );
    hp_Rnh_vs_pt->Fill( pt2ndJetReco, eNeutralHadrons2ndReco );
    hp_Rgamma_vs_pt->Fill( ptJetReco, ePhotonsReco );
    hp_Rgamma_vs_pt->Fill( pt2ndJetReco, ePhotons2ndReco );
    hp_Rmu_vs_pt->Fill( ptJetReco, eMuonsReco );
    hp_Rmu_vs_pt->Fill( pt2ndJetReco, eMuons2ndReco );
    hp_Re_vs_pt->Fill( ptJetReco, eElectronsReco );
    hp_Re_vs_pt->Fill( pt2ndJetReco, eElectrons2ndReco );
    hp_Rhfhad_vs_pt->Fill( ptJetReco, eHFHadronsReco );
    hp_Rhfhad_vs_pt->Fill( pt2ndJetReco, eHFHadrons2ndReco );
    hp_Rhfem_vs_pt->Fill( ptJetReco, eHFEMReco );
    hp_Rhfem_vs_pt->Fill( pt2ndJetReco, eHFEM2ndReco );
  
    hp_Ech_vs_pt->Fill( ptJetReco, eTracksReco*eJetReco );
    hp_Ech_vs_pt->Fill( pt2ndJetReco, eTracks2ndReco*e2ndJetReco );
    hp_Enh_vs_pt->Fill( ptJetReco, eNeutralHadronsReco*eJetReco );
    hp_Enh_vs_pt->Fill( pt2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco );
    hp_Egamma_vs_pt->Fill( ptJetReco, ePhotonsReco*eJetReco );
    hp_Egamma_vs_pt->Fill( pt2ndJetReco, ePhotons2ndReco*e2ndJetReco );
    hp_Emu_vs_pt->Fill( ptJetReco, eMuonsReco*eJetReco );
    hp_Emu_vs_pt->Fill( pt2ndJetReco, eMuons2ndReco*e2ndJetReco );
    hp_Ee_vs_pt->Fill( ptJetReco, eElectronsReco*eJetReco );
    hp_Ee_vs_pt->Fill( pt2ndJetReco, eElectrons2ndReco*e2ndJetReco );
    hp_Ehfhad_vs_pt->Fill( ptJetReco, eHFHadronsReco*eJetReco );
    hp_Ehfhad_vs_pt->Fill( pt2ndJetReco, eHFHadrons2ndReco*e2ndJetReco );
    hp_Ehfem_vs_pt->Fill( ptJetReco, eHFEMReco*eJetReco );
    hp_Ehfem_vs_pt->Fill( pt2ndJetReco, eHFEM2ndReco*e2ndJetReco );

    hp_Nch_vs_pt->Fill( ptJetReco, nTracksReco );
    hp_Nch_vs_pt->Fill( pt2ndJetReco, nTracks2ndReco );
    hp_Nnh_vs_pt->Fill( ptJetReco, nNeutralHadronsReco );
    hp_Nnh_vs_pt->Fill( pt2ndJetReco, nNeutralHadrons2ndReco );
    hp_Ngamma_vs_pt->Fill( ptJetReco, nPhotonsReco );
    hp_Ngamma_vs_pt->Fill( pt2ndJetReco, nPhotons2ndReco );
    hp_Nmu_vs_pt->Fill( ptJetReco, nMuonsReco );
    hp_Nmu_vs_pt->Fill( pt2ndJetReco, nMuons2ndReco );
    hp_Ne_vs_pt->Fill( ptJetReco, nElectronsReco );
    hp_Ne_vs_pt->Fill( pt2ndJetReco, nElectrons2ndReco );
    hp_Nhfhad_vs_pt->Fill( ptJetReco, nHFHadronsReco );
    hp_Nhfhad_vs_pt->Fill( pt2ndJetReco, nHFHadrons2ndReco );
    hp_Nhfem_vs_pt->Fill( ptJetReco, nHFEMReco );
    hp_Nhfem_vs_pt->Fill( pt2ndJetReco, nHFEM2ndReco );

    hp_Rch_vs_ptCorr->Fill( ptCorrJetReco, eTracksReco );
    hp_Rch_vs_ptCorr->Fill( ptCorr2ndJetReco, eTracks2ndReco );
    hp_Rnh_vs_ptCorr->Fill( ptCorrJetReco, eNeutralHadronsReco );
    hp_Rnh_vs_ptCorr->Fill( ptCorr2ndJetReco, eNeutralHadrons2ndReco );
    hp_Rgamma_vs_ptCorr->Fill( ptCorrJetReco, ePhotonsReco );
    hp_Rgamma_vs_ptCorr->Fill( ptCorr2ndJetReco, ePhotons2ndReco );
    hp_Rmu_vs_ptCorr->Fill( ptCorrJetReco, eMuonsReco );
    hp_Rmu_vs_ptCorr->Fill( ptCorr2ndJetReco, eMuons2ndReco );
    hp_Re_vs_ptCorr->Fill( ptCorrJetReco, eElectronsReco );
    hp_Re_vs_ptCorr->Fill( ptCorr2ndJetReco, eElectrons2ndReco );
    hp_Rhfhad_vs_ptCorr->Fill( ptCorrJetReco, eHFHadronsReco );
    hp_Rhfhad_vs_ptCorr->Fill( ptCorr2ndJetReco, eHFHadrons2ndReco );
    hp_Rhfem_vs_ptCorr->Fill( ptCorrJetReco, eHFEMReco );
    hp_Rhfem_vs_ptCorr->Fill( ptCorr2ndJetReco, eHFEM2ndReco );
  
    hp_Ech_vs_ptCorr->Fill( ptCorrJetReco, eTracksReco*eJetReco );
    hp_Ech_vs_ptCorr->Fill( ptCorr2ndJetReco, eTracks2ndReco*e2ndJetReco );
    hp_Enh_vs_ptCorr->Fill( ptCorrJetReco, eNeutralHadronsReco*eJetReco );
    hp_Enh_vs_ptCorr->Fill( ptCorr2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco );
    hp_Egamma_vs_ptCorr->Fill( ptCorrJetReco, ePhotonsReco*eJetReco );
    hp_Egamma_vs_ptCorr->Fill( ptCorr2ndJetReco, ePhotons2ndReco*e2ndJetReco );
    hp_Emu_vs_ptCorr->Fill( ptCorrJetReco, eMuonsReco*eJetReco );
    hp_Emu_vs_ptCorr->Fill( ptCorr2ndJetReco, eMuons2ndReco*e2ndJetReco );
    hp_Ee_vs_ptCorr->Fill( ptCorrJetReco, eElectronsReco*eJetReco );
    hp_Ee_vs_ptCorr->Fill( ptCorr2ndJetReco, eElectrons2ndReco*e2ndJetReco );
    hp_Ehfhad_vs_ptCorr->Fill( ptCorrJetReco, eHFHadronsReco*eJetReco );
    hp_Ehfhad_vs_ptCorr->Fill( ptCorr2ndJetReco, eHFHadrons2ndReco*e2ndJetReco );
    hp_Ehfem_vs_ptCorr->Fill( ptCorrJetReco, eHFEMReco*eJetReco );
    hp_Ehfem_vs_ptCorr->Fill( ptCorr2ndJetReco, eHFEM2ndReco*e2ndJetReco );

    hp_Nch_vs_ptCorr->Fill( ptCorrJetReco, nTracksReco );
    hp_Nch_vs_ptCorr->Fill( ptCorr2ndJetReco, nTracks2ndReco );
    hp_Nnh_vs_ptCorr->Fill( ptCorrJetReco, nNeutralHadronsReco );
    hp_Nnh_vs_ptCorr->Fill( ptCorr2ndJetReco, nNeutralHadrons2ndReco );
    hp_Ngamma_vs_ptCorr->Fill( ptCorrJetReco, nPhotonsReco );
    hp_Ngamma_vs_ptCorr->Fill( ptCorr2ndJetReco, nPhotons2ndReco );
    hp_Nmu_vs_ptCorr->Fill( ptCorrJetReco, nMuonsReco );
    hp_Nmu_vs_ptCorr->Fill( ptCorr2ndJetReco, nMuons2ndReco );
    hp_Ne_vs_ptCorr->Fill( ptCorrJetReco, nElectronsReco );
    hp_Ne_vs_ptCorr->Fill( ptCorr2ndJetReco, nElectrons2ndReco );
    hp_Nhfhad_vs_ptCorr->Fill( ptCorrJetReco, nHFHadronsReco );
    hp_Nhfhad_vs_ptCorr->Fill( ptCorr2ndJetReco, nHFHadrons2ndReco );
    hp_Nhfem_vs_ptCorr->Fill( ptCorrJetReco, nHFEMReco );
    hp_Nhfem_vs_ptCorr->Fill( ptCorr2ndJetReco, nHFEM2ndReco );


    //correlations between components:
    hp_Egamma_vs_Ech->Fill( eTracksReco*eJetReco, ePhotonsReco*eJetReco );
    hp_Egamma_vs_Ech->Fill( eTracks2ndReco*e2ndJetReco, ePhotons2ndReco*e2ndJetReco );
    hp_Enh_vs_Ech->Fill( eTracksReco*eJetReco, eNeutralHadronsReco*eJetReco );
    hp_Enh_vs_Ech->Fill( eTracks2ndReco*e2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco );

    hp_Rgamma_vs_Rch->Fill( eTracksReco, ePhotonsReco );
    hp_Rgamma_vs_Rch->Fill( eTracks2ndReco, ePhotons2ndReco );
    hp_Rnh_vs_Rch->Fill( eTracksReco, eNeutralHadronsReco );
    hp_Rnh_vs_Rch->Fill( eTracks2ndReco, eNeutralHadrons2ndReco );
    hp_Rnh_vs_Rgamma->Fill( ePhotonsReco, eNeutralHadronsReco );
    hp_Rnh_vs_Rgamma->Fill( ePhotons2ndReco, eNeutralHadrons2ndReco );

    hp_Ngamma_vs_Nch->Fill( nTracksReco, nPhotonsReco );
    hp_Ngamma_vs_Nch->Fill( nTracks2ndReco, nPhotons2ndReco );
    hp_Nnh_vs_Nch->Fill( nTracksReco, nNeutralHadronsReco );
    hp_Nnh_vs_Nch->Fill( nTracks2ndReco, nNeutralHadrons2ndReco );


  } //for entries


  std::string outfileName_str;

  if(dataset!="") outfileName_str = "MinBias_"+dataset;
  else outfileName_str = "MinBias";

  outfileName_str += "_" +algoType;

  char outfileName[150];
  if( use_deltaPhi )
    sprintf( outfileName, "%s_%dGeV%s.root", outfileName_str.c_str(), (Int_t)floor(pt_thresh), raw_corr.c_str());
  else 
    sprintf( outfileName, "%s_%dGeV%s_noDeltaPhi.root", outfileName_str.c_str(), (Int_t)floor(pt_thresh), raw_corr.c_str());

  TFile* outFile = new TFile(outfileName, "RECREATE");
  outFile->cd();

  h1_deltaPhiJet->Write();
  h1_asymmJet->Write();

  h1_massJet->Write();
  h1_diJetMass->Write();

  h1_etaJet->Write();
  h1_phiJet->Write();
  h1_ptJet->Write();
  h1_ptCorrJet->Write();
  h1_nCandJet->Write();

  h1_RchJet->Write();
  h1_RnhJet->Write();
  h1_RgammaJet->Write();
  h1_ReJet->Write();
  h1_RmuJet->Write();
  h1_RhfhadJet->Write();
  h1_RhfemJet->Write();

  h1_EchJet->Write();
  h1_EnhJet->Write();
  h1_EnhJet_barrel->Write();
  h1_EnhJet_endcap->Write();
  h1_EgammaJet->Write();
  h1_EeJet->Write();
  h1_EmuJet->Write();
  h1_EhfhadJet->Write();
  h1_EhfemJet->Write();

  h1_EphotAveJet->Write();

  h1_NchJet->Write();
  h1_NnhJet->Write();
  h1_NgammaJet->Write();
  h1_NeJet->Write();
  h1_NmuJet->Write();
  h1_NhfhadJet->Write();
  h1_NhfemJet->Write();

  hp_pt_vs_eta->Write();
  hp_ptCorr_vs_eta->Write();

  hp_Rch_vs_eta->Write();
  hp_Rnh_vs_eta->Write();
  hp_Rgamma_vs_eta->Write();
  hp_Rmu_vs_eta->Write();
  hp_Re_vs_eta->Write();
  hp_Rhfhad_vs_eta->Write();
  hp_Rhfem_vs_eta->Write();

  hp_Ech_vs_eta->Write();
  hp_Enh_vs_eta->Write();
  hp_Egamma_vs_eta->Write();
  hp_Emu_vs_eta->Write();
  hp_Ee_vs_eta->Write();
  hp_Ehfhad_vs_eta->Write();
  hp_Ehfem_vs_eta->Write();

  hp_Nch_vs_eta->Write();
  hp_Nnh_vs_eta->Write();
  hp_Ngamma_vs_eta->Write();
  hp_Nmu_vs_eta->Write();
  hp_Ne_vs_eta->Write();
  hp_Nhfhad_vs_eta->Write();
  hp_Nhfem_vs_eta->Write();

  hp_Rch_vs_pt->Write();
  hp_Rnh_vs_pt->Write();
  hp_Rgamma_vs_pt->Write();
  hp_Rmu_vs_pt->Write();
  hp_Re_vs_pt->Write();
  hp_Rhfhad_vs_pt->Write();
  hp_Rhfem_vs_pt->Write();

  hp_Ech_vs_pt->Write();
  hp_Enh_vs_pt->Write();
  hp_Egamma_vs_pt->Write();
  hp_Emu_vs_pt->Write();
  hp_Ee_vs_pt->Write();
  hp_Ehfhad_vs_pt->Write();
  hp_Ehfem_vs_pt->Write();

  hp_Nch_vs_pt->Write();
  hp_Nnh_vs_pt->Write();
  hp_Ngamma_vs_pt->Write();
  hp_Nmu_vs_pt->Write();
  hp_Ne_vs_pt->Write();
  hp_Nhfhad_vs_pt->Write();
  hp_Nhfem_vs_pt->Write();

  hp_Rch_vs_ptCorr->Write();
  hp_Rnh_vs_ptCorr->Write();
  hp_Rgamma_vs_ptCorr->Write();
  hp_Rmu_vs_ptCorr->Write();
  hp_Re_vs_ptCorr->Write();
  hp_Rhfhad_vs_ptCorr->Write();
  hp_Rhfem_vs_ptCorr->Write();

  hp_Ech_vs_ptCorr->Write();
  hp_Enh_vs_ptCorr->Write();
  hp_Egamma_vs_ptCorr->Write();
  hp_Emu_vs_ptCorr->Write();
  hp_Ee_vs_ptCorr->Write();
  hp_Ehfhad_vs_ptCorr->Write();
  hp_Ehfem_vs_ptCorr->Write();

  hp_Nch_vs_ptCorr->Write();
  hp_Nnh_vs_ptCorr->Write();
  hp_Ngamma_vs_ptCorr->Write();
  hp_Nmu_vs_ptCorr->Write();
  hp_Ne_vs_ptCorr->Write();
  hp_Nhfhad_vs_ptCorr->Write();
  hp_Nhfem_vs_ptCorr->Write();

  hp_Egamma_vs_Ech->Write();
  hp_Enh_vs_Ech->Write();

  hp_Rgamma_vs_Rch->Write();
  hp_Rnh_vs_Rch->Write();
  hp_Rnh_vs_Rgamma->Write();

  hp_Ngamma_vs_Nch->Write();
  hp_Nnh_vs_Nch->Write();

  outFile->Close();


}
  

bool PFJetID( Float_t Rch, Float_t Rnh, Float_t Rgamma, Int_t N, Int_t Nch) {

  bool pass = false;

  if( Rnh<1. && Rch>0. && Rgamma<1. ) pass = true;
  //if( (N>1) && (Rch>0.15) ) pass=true;
  //if( Nch>1 ) pass=true;

  return pass;

}
