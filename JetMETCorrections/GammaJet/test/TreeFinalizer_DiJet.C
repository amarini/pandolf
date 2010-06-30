#include <TH1F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "fitTools.C"




bool PFJetID( Float_t Rch, Float_t Rnh, Float_t Rgamma, Int_t N, Int_t Nch, Float_t eta);



void finalize(const std::string& dataset, std::string algoType, Float_t pt_thresh=5., std::string raw_corr="corr", bool useJetID=false) {

  Float_t etaMax = 3.;

  TChain* tree = new TChain("jetTree");

  std::string infileName, treeName;
  TFile* infile=0;
  TH1F* h1_lumi=0;
  Double_t totalLumi=0.;

  if( dataset=="DATA_7TeV_20100407") {

    //infileName = "DiJet_2ndLevelTree_DATA_7TeV_GOODCOLL_Apr1stSkim_" + algoType + ".root/jetTree";
    //tree->Add(infileName.c_str());
    //std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "DiJet_2ndLevelTree_DATA_MinimumBias_Commissioning10-GOODCOLL-v8_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "DiJet_2ndLevelTree_Commissioning10-GOODCOLL-v8_runs132605_606_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
  
  } else if( dataset=="QCD_Spring10") {

    infileName = "DiJet_2ndLevelTree_QCD_Spring10_Pt0to15_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "DiJet_2ndLevelTree_QCD_Spring10_Pt15_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "DiJet_2ndLevelTree_QCD_Spring10_Pt30_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "DiJet_2ndLevelTree_QCD_Spring10_Pt80_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "DiJet_2ndLevelTree_QCD_Spring10_Pt170_" + algoType + "_NOPTHATCUT.root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else {

    infileName = "files_DiJet_2ndLevel_" + dataset+".txt";

    //open from file.txt:
    FILE* iff = fopen(infileName.c_str(),"r");
    if(iff == 0) {
      std::cout << "cannot open input file " << infileName << " ... adding single file." << std::endl;
      infileName = "Inclusive_2ndLevelTree_" + dataset + "_" + algoType + ".root/jetTree";
      tree->Add(infileName.c_str());
      std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    } else {

      char singleLine[500];

      while( fscanf(iff, "%s", singleLine) !=EOF ) {

        std::string rootfilename(singleLine);
        std::string treename = rootfilename + "/jetTree";
        std::cout << "-> Adding " << treename;
        tree->Add(treename.c_str());
        TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
        h1_lumi = (TH1F*)infile->Get("lumi");
        if( h1_lumi!=0 ) {
          totalLumi += h1_lumi->GetBinContent(1);
          std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
        } else {
          std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
        }

      }
      fclose(iff);

    }
 
//  infileName = "DiJet_2ndLevelTree_" + dataset + "_" + algoType + ".root";
//  infile = TFile::Open(infileName.c_str(), "READ");
//  h1_lumi = (TH1F*)infile->Get("lumi");
//  if( h1_lumi!=0 ) {
//    totalLumi += h1_lumi->GetBinContent(1);
//  } else {
//    std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
//  }
//  treeName = infileName + "/jetTree";
//  tree->Add(treeName.c_str());
//  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  }
  
  std::cout << "Done adding." << std::endl;

  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);


  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Int_t lbn;
  tree->SetBranchAddress("lbn", &lbn);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);


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
  Float_t ePhotonsGen;
  tree->SetBranchAddress("ePhotonsGen", &ePhotonsGen);
  Float_t eNeutralHadronsGen;
  tree->SetBranchAddress("eNeutralHadronsGen", &eNeutralHadronsGen);
  Float_t eMuonsGen;
  tree->SetBranchAddress("eMuonsGen", &eMuonsGen);
  Float_t eElectronsGen;
  tree->SetBranchAddress("eElectronsGen", &eElectronsGen);
  Float_t eHFHadronsGen;
  tree->SetBranchAddress("eHFHadronsGen", &eHFHadronsGen);
  Float_t eHFEMGen;
  tree->SetBranchAddress("eHFEMGen", &eHFEMGen);

  Float_t ptTracksReco;
  tree->SetBranchAddress("ptTracksReco", &ptTracksReco);
  Float_t ptPhotonsReco;
  tree->SetBranchAddress("ptPhotonsReco", &ptPhotonsReco);
  Float_t ptNeutralHadronsReco;
  tree->SetBranchAddress("ptNeutralHadronsReco", &ptNeutralHadronsReco);
  Float_t ptMuonsReco;
  tree->SetBranchAddress("ptMuonsReco", &ptMuonsReco);
  Float_t ptElectronsReco;
  tree->SetBranchAddress("ptElectronsReco", &ptElectronsReco);
  Float_t ptHFHadronsReco;
  tree->SetBranchAddress("ptHFHadronsReco", &ptHFHadronsReco);
  Float_t ptHFEMReco;
  tree->SetBranchAddress("ptHFEMReco", &ptHFEMReco);

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

  Float_t e2ndJetGen;
  tree->SetBranchAddress("e2ndJetGen", &e2ndJetGen);
  Float_t pt2ndJetGen;
  tree->SetBranchAddress("pt2ndJetGen", &pt2ndJetGen);
  Float_t eta2ndJetGen;
  tree->SetBranchAddress("eta2ndJetGen", &eta2ndJetGen);
  Float_t phi2ndJetGen;
  tree->SetBranchAddress("phi2ndJetGen", &phi2ndJetGen);

  Float_t eTracks2ndGen;
  tree->SetBranchAddress("eTracks2ndGen", &eTracks2ndGen);
  Float_t ePhotons2ndGen;
  tree->SetBranchAddress("ePhotons2ndGen", &ePhotons2ndGen);
  Float_t eNeutralHadrons2ndGen;
  tree->SetBranchAddress("eNeutralHadrons2ndGen", &eNeutralHadrons2ndGen);
  Float_t eMuons2ndGen;
  tree->SetBranchAddress("eMuons2ndGen", &eMuons2ndGen);
  Float_t eElectrons2ndGen;
  tree->SetBranchAddress("eElectrons2ndGen", &eElectrons2ndGen);
  Float_t eHFHadrons2ndGen;
  tree->SetBranchAddress("eHFHadrons2ndGen", &eHFHadrons2ndGen);
  Float_t eHFEM2ndGen;
  tree->SetBranchAddress("eHFEM2ndGen", &eHFEM2ndGen);

  Float_t ptTracks2ndReco;
  tree->SetBranchAddress("ptTracks2ndReco", &ptTracks2ndReco);
  Float_t ptPhotons2ndReco;
  tree->SetBranchAddress("ptPhotons2ndReco", &ptPhotons2ndReco);
  Float_t ptNeutralHadrons2ndReco;
  tree->SetBranchAddress("ptNeutralHadrons2ndReco", &ptNeutralHadrons2ndReco);
  Float_t ptMuons2ndReco;
  tree->SetBranchAddress("ptMuons2ndReco", &ptMuons2ndReco);
  Float_t ptElectrons2ndReco;
  tree->SetBranchAddress("ptElectrons2ndReco", &ptElectrons2ndReco);
  Float_t ptHFHadrons2ndReco;
  tree->SetBranchAddress("ptHFHadrons2ndReco", &ptHFHadrons2ndReco);
  Float_t ptHFEM2ndReco;
  tree->SetBranchAddress("ptHFEM2ndReco", &ptHFEM2ndReco);


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
  Int_t nBins_E = 30;
  Int_t nBins_R = 30;

  Float_t eMax = 600.;
  Float_t ptMax = 300.;
  
  TH1F* h1_deltaPhiJet = new TH1F("deltaPhiJet", "", 40, 1., 3.152);
  TH1F* h1_asymmJet = new TH1F("asymmJet", "", 20, -0.8, 0.8);

  TH1F* h1_massJet = new TH1F("massJet", "", 80, 0., 40.);
  TH1F* h1_ptOverMJet = new TH1F("ptOverMJet", "", 80, 0., 60.);
  TH1F* h1_pOverMJet = new TH1F("pOverMJet", "", 80, 0., 200.);
  TH1F* h1_diJetMass = new TH1F("diJetMass", "", 50., 0., 300.);

  TH1F* h1_etaJet = new TH1F("etaJet", "", nBins_eta, -etaMax, etaMax);
  TH1F* h1_phiJet = new TH1F("phiJet", "", nBins_eta, -3.5, 3.5);
  TH1F* h1_ptJet = new TH1F("ptJet", "", 50, 0., 150.);
  TH1F* h1_ptCorrJet = new TH1F("ptCorrJet", "", 60, 0., ptMax);
  TH1F* h1_ptCorrJet_barrel = new TH1F("ptCorrJet_barrel", "", 60, 0., ptMax);
  TH1F* h1_ptCorrJet_endcap = new TH1F("ptCorrJet_endcap", "", 60, 0., ptMax);
  TH1F* h1_ptCorrJet_eta1430 = new TH1F("ptCorrJet_eta1430", "", 60, 0., ptMax);
  TH1F* h1_ptCorrJet_Rch050 = new TH1F("ptCorrJet_Rch050", "", 60, 0., ptMax);
  TH1F* h1_ptCorrJet_Rch5070 = new TH1F("ptCorrJet_Rch5070", "", 60, 0., ptMax);
  TH1F* h1_ptCorrJet_Rch70100 = new TH1F("ptCorrJet_Rch70100", "", 60, 0., ptMax);
  TH1F* h1_nCandJet = new TH1F("nCandJet", "", 41, -0.5, 40.);

  TH1F* h1_RchJet = new TH1F("RchJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet = new TH1F("RnhJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet = new TH1F("RgammaJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet = new TH1F("RmuJet", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet = new TH1F("ReJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet = new TH1F("RhfhadJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet = new TH1F("RhfemJet", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet = new TH1F("EchJet", "", nBins_E, 0., eMax);
  TH1F* h1_EnhJet = new TH1F("EnhJet", "", nBins_E, 0., eMax);
  TH1F* h1_EgammaJet = new TH1F("EgammaJet", "", nBins_E, 0., eMax);
  TH1F* h1_EeJet = new TH1F("EeJet", "", nBins_E, 0., eMax);
  TH1F* h1_EmuJet = new TH1F("EmuJet", "", nBins_E, 0., eMax);
  TH1F* h1_EhfhadJet = new TH1F("EhfhadJet", "", nBins_E, 0., eMax);
  TH1F* h1_EhfemJet = new TH1F("EhfemJet", "", nBins_E, 0., eMax);

  TH1F* h1_PTchJet = new TH1F("PTchJet", "", nBins_E, 0., ptMax);
  TH1F* h1_PTnhJet = new TH1F("PTnhJet", "", nBins_E, 0., ptMax);
  TH1F* h1_PTgammaJet = new TH1F("PTgammaJet", "", nBins_E, 0., ptMax);
  TH1F* h1_PTeJet = new TH1F("PTeJet", "", nBins_E, 0., ptMax);
  TH1F* h1_PTmuJet = new TH1F("PTmuJet", "", nBins_E, 0., ptMax);
  TH1F* h1_PThfhadJet = new TH1F("PThfhadJet", "", nBins_E, 0., ptMax);
  TH1F* h1_PThfemJet = new TH1F("PThfemJet", "", nBins_E, 0., ptMax);

  TH1F* h1_EphotAveJet = new TH1F("EphotAveJet", "", 40, 0., 10.);

  TH1F* h1_NchJet = new TH1F("NchJet", "", 31, -0.5, 30.5);
  TH1F* h1_NnhJet = new TH1F("NnhJet", "", 31, -0.5, 30.5);
  TH1F* h1_NgammaJet = new TH1F("NgammaJet", "", 31, -0.5, 30.5);
  TH1F* h1_NeJet = new TH1F("NeJet", "", 31, -0.5, 30.5);
  TH1F* h1_NmuJet = new TH1F("NmuJet", "", 31, -0.5, 30.5);
  TH1F* h1_NhfhadJet = new TH1F("NhfhadJet", "", 31, -0.5, 30.5);
  TH1F* h1_NhfemJet = new TH1F("NhfemJet", "", 31, -0.5, 30.5);

  TH1F* h1_RchJet_barrel = new TH1F("RchJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet_barrel = new TH1F("RnhJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet_barrel = new TH1F("RgammaJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet_barrel = new TH1F("RmuJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet_barrel = new TH1F("ReJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet_barrel = new TH1F("RhfhadJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet_barrel = new TH1F("RhfemJet_barrel", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet_barrel = new TH1F("EchJet_barrel", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EnhJet_barrel = new TH1F("EnhJet_barrel", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EgammaJet_barrel = new TH1F("EgammaJet_barrel", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EeJet_barrel = new TH1F("EeJet_barrel", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EmuJet_barrel = new TH1F("EmuJet_barrel", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EhfhadJet_barrel = new TH1F("EhfhadJet_barrel", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EhfemJet_barrel = new TH1F("EhfemJet_barrel", "", nBins_E, 0., eMax/2.);

  TH1F* h1_NchJet_barrel = new TH1F("NchJet_barrel", "", 31, -0.5, 30.5);
  TH1F* h1_NnhJet_barrel = new TH1F("NnhJet_barrel", "", 31, -0.5, 30.5);
  TH1F* h1_NgammaJet_barrel = new TH1F("NgammaJet_barrel", "", 31, -0.5, 30.5);
  TH1F* h1_NeJet_barrel = new TH1F("NeJet_barrel", "", 31, -0.5, 30.5);
  TH1F* h1_NmuJet_barrel = new TH1F("NmuJet_barrel", "", 31, -0.5, 30.5);
  TH1F* h1_NhfhadJet_barrel = new TH1F("NhfhadJet_barrel", "", 31, -0.5, 30.5);
  TH1F* h1_NhfemJet_barrel = new TH1F("NhfemJet_barrel", "", 31, -0.5, 30.5);

  TH1F* h1_RchJet_eta1425 = new TH1F("RchJet_eta1425", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet_eta1425 = new TH1F("RnhJet_eta1425", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet_eta1425 = new TH1F("RgammaJet_eta1425", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet_eta1425 = new TH1F("RmuJet_eta1425", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet_eta1425 = new TH1F("ReJet_eta1425", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet_eta1425 = new TH1F("RhfhadJet_eta1425", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet_eta1425 = new TH1F("RhfemJet_eta1425", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet_eta1425 = new TH1F("EchJet_eta1425", "", nBins_E, 0., eMax);
  TH1F* h1_EnhJet_eta1425 = new TH1F("EnhJet_eta1425", "", nBins_E, 0., eMax);
  TH1F* h1_EgammaJet_eta1425 = new TH1F("EgammaJet_eta1425", "", nBins_E, 0., eMax);
  TH1F* h1_EeJet_eta1425 = new TH1F("EeJet_eta1425", "", nBins_E, 0., eMax);
  TH1F* h1_EmuJet_eta1425 = new TH1F("EmuJet_eta1425", "", nBins_E, 0., eMax);
  TH1F* h1_EhfhadJet_eta1425 = new TH1F("EhfhadJet_eta1425", "", nBins_E, 0., eMax);
  TH1F* h1_EhfemJet_eta1425 = new TH1F("EhfemJet_eta1425", "", nBins_E, 0., eMax);

  TH1F* h1_NchJet_eta1425 = new TH1F("NchJet_eta1425", "", 31, -0.5, 30.5);
  TH1F* h1_NnhJet_eta1425 = new TH1F("NnhJet_eta1425", "", 31, -0.5, 30.5);
  TH1F* h1_NgammaJet_eta1425 = new TH1F("NgammaJet_eta1425", "", 31, -0.5, 30.5);
  TH1F* h1_NeJet_eta1425 = new TH1F("NeJet_eta1425", "", 31, -0.5, 30.5);
  TH1F* h1_NmuJet_eta1425 = new TH1F("NmuJet_eta1425", "", 31, -0.5, 30.5);
  TH1F* h1_NhfhadJet_eta1425 = new TH1F("NhfhadJet_eta1425", "", 31, -0.5, 30.5);
  TH1F* h1_NhfemJet_eta1425 = new TH1F("NhfemJet_eta1425", "", 31, -0.5, 30.5);

//Int_t iLower=0;
//Double_t etaMaxStack = 3.8;
//Double_t binWidth = 2.*etaMaxStack/45.;
//Double_t Lower_stack[2+45+1];
//Lower_stack[iLower++] = -5.;
//Lower_stack[iLower++] = -etaMaxStack;
//for( Double_t iEta=(-etaMaxStack+binWidth); iEta < etaMaxStack; iEta+=binWidth )
//  Lower_stack[iLower++] = iEta;
//Lower_stack[iLower++] = etaMaxStack;
//Lower_stack[iLower++] = 5.;

  int nBins_stack_eta = fitTools::getNbins_stack("eta");
  Double_t Lower_stack_eta[nBins_stack_eta];
  fitTools::getBins_stack( nBins_stack_eta, Lower_stack_eta, "eta");

  TProfile* hp_Rch_vs_eta_stack = new TProfile("Rch_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_Rnh_vs_eta_stack = new TProfile("Rnh_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_Rgamma_vs_eta_stack = new TProfile("Rgamma_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_Rmu_vs_eta_stack = new TProfile("Rmu_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_Re_vs_eta_stack = new TProfile("Re_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_Rhfhad_vs_eta_stack = new TProfile("Rhfhad_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_Rhfem_vs_eta_stack = new TProfile("Rhfem_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);

  TProfile* hp_RchGen_vs_eta_stack = new TProfile("RchGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_RnhGen_vs_eta_stack = new TProfile("RnhGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_RgammaGen_vs_eta_stack = new TProfile("RgammaGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_RmuGen_vs_eta_stack = new TProfile("RmuGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_ReGen_vs_eta_stack = new TProfile("ReGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_RhfhadGen_vs_eta_stack = new TProfile("RhfhadGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);
  TProfile* hp_RhfemGen_vs_eta_stack = new TProfile("RhfemGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 10.);

  Float_t etaMaxProfile = 5.;

  TProfile* hp_pt_vs_eta = new TProfile("pt_vs_eta", "pt vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 300.);
  TProfile* hp_ptCorr_vs_eta = new TProfile("ptCorr_vs_eta", "ptCorr vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 300.);

  TProfile* hp_Rch_vs_eta = new TProfile("Rch_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Rnh_vs_eta = new TProfile("Rnh_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Rgamma_vs_eta = new TProfile("Rgamma_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Rmu_vs_eta = new TProfile("Rmu_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Re_vs_eta = new TProfile("Re_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Rhfhad_vs_eta = new TProfile("Rhfhad_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Rhfem_vs_eta = new TProfile("Rhfem_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);

  TProfile* hp_Ech_vs_eta = new TProfile("Ech_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Enh_vs_eta = new TProfile("Enh_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Egamma_vs_eta = new TProfile("Egamma_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Ee_vs_eta = new TProfile("Ee_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Emu_vs_eta = new TProfile("Emu_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Ehfhad_vs_eta = new TProfile("Ehfhad_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Ehfem_vs_eta = new TProfile("Ehfem_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);

  TProfile* hp_PTch_vs_eta = new TProfile("PTch_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_PTnh_vs_eta = new TProfile("PTnh_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_PTgamma_vs_eta = new TProfile("PTgamma_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_PTe_vs_eta = new TProfile("PTe_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_PTmu_vs_eta = new TProfile("PTmu_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_PThfhad_vs_eta = new TProfile("PThfhad_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_PThfem_vs_eta = new TProfile("PThfem_vs_eta", "R_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);

  TProfile* hp_Nch_vs_eta = new TProfile("Nch_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Nnh_vs_eta = new TProfile("Nnh_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Ngamma_vs_eta = new TProfile("Ngamma_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Ne_vs_eta = new TProfile("Ne_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Nmu_vs_eta = new TProfile("Nmu_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Nhfhad_vs_eta = new TProfile("Nhfhad_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);
  TProfile* hp_Nhfem_vs_eta = new TProfile("Nhfem_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMaxProfile, etaMaxProfile, 0., 1000.);

  Int_t nBins_pt_lower = fitTools::getNbins_stack("pt");;
  Double_t Lower[nBins_pt_lower];
  fitTools::getBins_stack( nBins_pt_lower, Lower, "pt" );

  TProfile* hp_Rch_vs_pt_barrel = new TProfile("Rch_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rnh_vs_pt_barrel = new TProfile("Rnh_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rgamma_vs_pt_barrel = new TProfile("Rgamma_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rmu_vs_pt_barrel = new TProfile("Rmu_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Re_vs_pt_barrel = new TProfile("Re_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfhad_vs_pt_barrel = new TProfile("Rhfhad_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfem_vs_pt_barrel = new TProfile("Rhfem_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Ech_vs_pt_barrel = new TProfile("Ech_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Enh_vs_pt_barrel = new TProfile("Enh_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Egamma_vs_pt_barrel = new TProfile("Egamma_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ee_vs_pt_barrel = new TProfile("Ee_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Emu_vs_pt_barrel = new TProfile("Emu_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ehfhad_vs_pt_barrel = new TProfile("Ehfhad_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ehfem_vs_pt_barrel = new TProfile("Ehfem_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Nch_vs_pt_barrel = new TProfile("Nch_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nnh_vs_pt_barrel = new TProfile("Nnh_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ngamma_vs_pt_barrel = new TProfile("Ngamma_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ne_vs_pt_barrel = new TProfile("Ne_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nmu_vs_pt_barrel = new TProfile("Nmu_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nhfhad_vs_pt_barrel = new TProfile("Nhfhad_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nhfem_vs_pt_barrel = new TProfile("Nhfem_vs_pt_barrel", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_RchGen_vs_pt_barrel = new TProfile("RchGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_RnhGen_vs_pt_barrel = new TProfile("RnhGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_RgammaGen_vs_pt_barrel = new TProfile("RgammaGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_RmuGen_vs_pt_barrel = new TProfile("RmuGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_ReGen_vs_pt_barrel = new TProfile("ReGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_RhfhadGen_vs_pt_barrel = new TProfile("RhfhadGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_RhfemGen_vs_pt_barrel = new TProfile("RhfemGen_vs_pt_barrel", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Rch_vs_ptCorr_barrel = new TProfile("Rch_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rnh_vs_ptCorr_barrel = new TProfile("Rnh_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rgamma_vs_ptCorr_barrel = new TProfile("Rgamma_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rmu_vs_ptCorr_barrel = new TProfile("Rmu_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Re_vs_ptCorr_barrel = new TProfile("Re_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfhad_vs_ptCorr_barrel = new TProfile("Rhfhad_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfem_vs_ptCorr_barrel = new TProfile("Rhfem_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);

  int nBins_stack_phi = fitTools::getNbins_stack("phi");
  Double_t Lower_stack_phi[nBins_stack_phi];
  fitTools::getBins_stack( nBins_stack_phi, Lower_stack_phi, "phi");

  TProfile* hp_Rch_vs_phi_barrel = new TProfile("Rch_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Rnh_vs_phi_barrel = new TProfile("Rnh_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Rgamma_vs_phi_barrel = new TProfile("Rgamma_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Rmu_vs_phi_barrel = new TProfile("Rmu_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Re_vs_phi_barrel = new TProfile("Re_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Rhfhad_vs_phi_barrel = new TProfile("Rhfhad_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Rhfem_vs_phi_barrel = new TProfile("Rhfem_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);

  TProfile* hp_RchGen_vs_phi_barrel = new TProfile("RchGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_RnhGen_vs_phi_barrel = new TProfile("RnhGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_RgammaGen_vs_phi_barrel = new TProfile("RgammaGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_RmuGen_vs_phi_barrel = new TProfile("RmuGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_ReGen_vs_phi_barrel = new TProfile("ReGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_RhfhadGen_vs_phi_barrel = new TProfile("RhfhadGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_RhfemGen_vs_phi_barrel = new TProfile("RhfemGen_vs_phi_barrel", "R_{ch} vs. phi", nBins_stack_phi-1, Lower_stack_phi,  0., 10.);
  TProfile* hp_Egamma_vs_Ech = new TProfile("Egamma_vs_Ech", "", nBins_E, 0., 15., 0., 1000.);
  TProfile* hp_Enh_vs_Ech = new TProfile("Enh_vs_Ech", "", nBins_E, 0., 15., 0., 1000.);

  TProfile* hp_Rgamma_vs_Rch = new TProfile("Rgamma_vs_Rch", "", nBins_R, 0., 1.1, 0., 1.1);
  TProfile* hp_Rnh_vs_Rch = new TProfile("Rnh_vs_Rch", "", nBins_R, 0., 1.1, 0., 1.1);
  TProfile* hp_Rgamma_vs_Rnh = new TProfile("Rgamma_vs_Rnh", "", nBins_R, 0., 1.1, 0., 1.1);

  TProfile* hp_Ngamma_vs_Nch = new TProfile("Ngamma_vs_Nch", "", 21, -0.5, 20.5, -0.5, 20.5);
  TProfile* hp_Nnh_vs_Nch = new TProfile("Nnh_vs_Nch", "", 21, -0.5, 20.5, -0.5, 20.5);

  TRandom* rand = new TRandom();

  int nEntries = tree->GetEntries();

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;
    if( dataset=="MinimumBias_Commissioning10-GOODCOLL-v9" ) eventWeight=1.;

    //PFJetID:
    Int_t N = nTracksReco+nPhotonsReco+nNeutralHadronsReco+nElectronsReco+nMuonsReco+nHFHadronsReco+nHFEMReco;
    Int_t Nch = nTracksReco+nElectronsReco+nMuonsReco;
    bool pfJetIDOK = PFJetID( eTracksReco, eNeutralHadronsReco, ePhotonsReco, N, Nch, etaJetReco);
    
    Int_t N2nd = nTracks2ndReco+nPhotons2ndReco+nNeutralHadrons2ndReco+nElectrons2ndReco+nMuons2ndReco+nHFHadrons2ndReco+nHFEM2ndReco;
    Int_t Nch2nd = nTracks2ndReco+nElectrons2ndReco+nMuons2ndReco;
    bool pf2ndJetIDOK = PFJetID( eTracks2ndReco, eNeutralHadrons2ndReco, ePhotons2ndReco, N2nd, Nch2nd, eta2ndJetReco);

    //fill pt stacks before pt cut:
    if( (useJetID && pfJetIDOK) || (!useJetID) ) {

      if( ptJetGen > pt_thresh ) {

        hp_RchGen_vs_eta_stack->Fill( etaJetGen, eTracksGen, eventWeight );
        hp_RnhGen_vs_eta_stack->Fill( etaJetGen, eNeutralHadronsGen, eventWeight );
        hp_RgammaGen_vs_eta_stack->Fill( etaJetGen, ePhotonsGen, eventWeight );
        hp_RmuGen_vs_eta_stack->Fill( etaJetGen, eMuonsGen, eventWeight );
        hp_ReGen_vs_eta_stack->Fill( etaJetGen, eElectronsGen, eventWeight );
        hp_RhfhadGen_vs_eta_stack->Fill( etaJetGen, eHFHadronsGen, eventWeight );
        hp_RhfemGen_vs_eta_stack->Fill( etaJetGen, eHFEMGen, eventWeight );

      }

      if( fabs(etaJetReco)<1.4 ) {
    
        hp_Rch_vs_pt_barrel->Fill( ptJetReco, eTracksReco, eventWeight );
        hp_Rnh_vs_pt_barrel->Fill( ptJetReco, eNeutralHadronsReco, eventWeight );
        hp_Rgamma_vs_pt_barrel->Fill( ptJetReco, ePhotonsReco, eventWeight );
        hp_Rmu_vs_pt_barrel->Fill( ptJetReco, eMuonsReco, eventWeight );
        hp_Re_vs_pt_barrel->Fill( ptJetReco, eElectronsReco, eventWeight );
        hp_Rhfhad_vs_pt_barrel->Fill( ptJetReco, eHFHadronsReco, eventWeight );
        hp_Rhfem_vs_pt_barrel->Fill( ptJetReco, eHFEMReco, eventWeight );
      
        hp_Ech_vs_pt_barrel->Fill( ptJetReco, eTracksReco*eJetReco, eventWeight );
        hp_Enh_vs_pt_barrel->Fill( ptJetReco, eNeutralHadronsReco*eJetReco, eventWeight );
        hp_Egamma_vs_pt_barrel->Fill( ptJetReco, ePhotonsReco*eJetReco, eventWeight );
        hp_Emu_vs_pt_barrel->Fill( ptJetReco, eMuonsReco*eJetReco, eventWeight );
        hp_Ee_vs_pt_barrel->Fill( ptJetReco, eElectronsReco*eJetReco, eventWeight );
        hp_Ehfhad_vs_pt_barrel->Fill( ptJetReco, eHFHadronsReco*eJetReco, eventWeight );
        hp_Ehfem_vs_pt_barrel->Fill( ptJetReco, eHFEMReco*eJetReco, eventWeight );
    
        hp_Nch_vs_pt_barrel->Fill( ptJetReco, nTracksReco, eventWeight );
        hp_Nnh_vs_pt_barrel->Fill( ptJetReco, nNeutralHadronsReco, eventWeight );
        hp_Ngamma_vs_pt_barrel->Fill( ptJetReco, nPhotonsReco, eventWeight );
        hp_Nmu_vs_pt_barrel->Fill( ptJetReco, nMuonsReco, eventWeight );
        hp_Ne_vs_pt_barrel->Fill( ptJetReco, nElectronsReco, eventWeight );
        hp_Nhfhad_vs_pt_barrel->Fill( ptJetReco, nHFHadronsReco, eventWeight );
        hp_Nhfem_vs_pt_barrel->Fill( ptJetReco, nHFEMReco, eventWeight );
    
        hp_Rch_vs_ptCorr_barrel->Fill( ptCorrJetReco, eTracksReco, eventWeight );
        hp_Rnh_vs_ptCorr_barrel->Fill( ptCorrJetReco, eNeutralHadronsReco, eventWeight );
        hp_Rgamma_vs_ptCorr_barrel->Fill( ptCorrJetReco, ePhotonsReco, eventWeight );
        hp_Rmu_vs_ptCorr_barrel->Fill( ptCorrJetReco, eMuonsReco, eventWeight );
        hp_Re_vs_ptCorr_barrel->Fill( ptCorrJetReco, eElectronsReco, eventWeight );
        hp_Rhfhad_vs_ptCorr_barrel->Fill( ptCorrJetReco, eHFHadronsReco, eventWeight );
        hp_Rhfem_vs_ptCorr_barrel->Fill( ptCorrJetReco, eHFEMReco, eventWeight );
    

        if( fabs(etaJetGen)<1.4 ) {

          hp_RchGen_vs_pt_barrel->Fill( ptJetGen, eTracksGen, eventWeight );
          hp_RnhGen_vs_pt_barrel->Fill( ptJetGen, eNeutralHadronsGen, eventWeight );
          hp_RgammaGen_vs_pt_barrel->Fill( ptJetGen, ePhotonsGen, eventWeight );
          hp_RmuGen_vs_pt_barrel->Fill( ptJetGen, eMuonsGen, eventWeight );
          hp_ReGen_vs_pt_barrel->Fill( ptJetGen, eElectronsGen, eventWeight );
          hp_RhfhadGen_vs_pt_barrel->Fill( ptJetGen, eHFHadronsGen, eventWeight );
          hp_RhfemGen_vs_pt_barrel->Fill( ptJetGen, eHFEMGen, eventWeight );

          if( ptJetGen>pt_thresh ) {
            hp_RchGen_vs_phi_barrel->Fill( phiJetGen, eTracksGen, eventWeight );
            hp_RnhGen_vs_phi_barrel->Fill( phiJetGen, eNeutralHadronsGen, eventWeight );
            hp_RgammaGen_vs_phi_barrel->Fill( phiJetGen, ePhotonsGen, eventWeight );
            hp_RmuGen_vs_phi_barrel->Fill( phiJetGen, eMuonsGen, eventWeight );
            hp_ReGen_vs_phi_barrel->Fill( phiJetGen, eElectronsGen, eventWeight );
            hp_RhfhadGen_vs_phi_barrel->Fill( phiJetGen, eHFHadronsGen, eventWeight );
            hp_RhfemGen_vs_phi_barrel->Fill( phiJetGen, eHFEMGen, eventWeight );
          }
        }
      }

      //second jet:
      if( pt2ndJetGen > pt_thresh ) {

        hp_RchGen_vs_eta_stack->Fill( eta2ndJetGen, eTracks2ndGen, eventWeight );
        hp_RnhGen_vs_eta_stack->Fill( eta2ndJetGen, eNeutralHadrons2ndGen, eventWeight );
        hp_RgammaGen_vs_eta_stack->Fill( eta2ndJetGen, ePhotons2ndGen, eventWeight );
        hp_RmuGen_vs_eta_stack->Fill( eta2ndJetGen, eMuons2ndGen, eventWeight );
        hp_ReGen_vs_eta_stack->Fill( eta2ndJetGen, eElectrons2ndGen, eventWeight );
        hp_RhfhadGen_vs_eta_stack->Fill( eta2ndJetGen, eHFHadrons2ndGen, eventWeight );
        hp_RhfemGen_vs_eta_stack->Fill( eta2ndJetGen, eHFEM2ndGen, eventWeight );

      }
      if( fabs(eta2ndJetReco)<1.4 ) {
    
        hp_Rch_vs_pt_barrel->Fill( pt2ndJetReco, eTracks2ndReco, eventWeight );
        hp_Rnh_vs_pt_barrel->Fill( pt2ndJetReco, eNeutralHadrons2ndReco, eventWeight );
        hp_Rgamma_vs_pt_barrel->Fill( pt2ndJetReco, ePhotons2ndReco, eventWeight );
        hp_Rmu_vs_pt_barrel->Fill( pt2ndJetReco, eMuons2ndReco, eventWeight );
        hp_Re_vs_pt_barrel->Fill( pt2ndJetReco, eElectrons2ndReco, eventWeight );
        hp_Rhfhad_vs_pt_barrel->Fill( pt2ndJetReco, eHFHadrons2ndReco, eventWeight );
        hp_Rhfem_vs_pt_barrel->Fill( pt2ndJetReco, eHFEM2ndReco, eventWeight );
      
        hp_Ech_vs_pt_barrel->Fill( pt2ndJetReco, eTracks2ndReco*e2ndJetReco, eventWeight );
        hp_Enh_vs_pt_barrel->Fill( pt2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco, eventWeight );
        hp_Egamma_vs_pt_barrel->Fill( pt2ndJetReco, ePhotons2ndReco*e2ndJetReco, eventWeight );
        hp_Emu_vs_pt_barrel->Fill( pt2ndJetReco, eMuons2ndReco*e2ndJetReco, eventWeight );
        hp_Ee_vs_pt_barrel->Fill( pt2ndJetReco, eElectrons2ndReco*e2ndJetReco, eventWeight );
        hp_Ehfhad_vs_pt_barrel->Fill( pt2ndJetReco, eHFHadrons2ndReco*e2ndJetReco, eventWeight );
        hp_Ehfem_vs_pt_barrel->Fill( pt2ndJetReco, eHFEM2ndReco*e2ndJetReco, eventWeight );
    
        hp_Nch_vs_pt_barrel->Fill( pt2ndJetReco, nTracks2ndReco, eventWeight );
        hp_Nnh_vs_pt_barrel->Fill( pt2ndJetReco, nNeutralHadrons2ndReco, eventWeight );
        hp_Ngamma_vs_pt_barrel->Fill( pt2ndJetReco, nPhotons2ndReco, eventWeight );
        hp_Nmu_vs_pt_barrel->Fill( pt2ndJetReco, nMuons2ndReco, eventWeight );
        hp_Ne_vs_pt_barrel->Fill( pt2ndJetReco, nElectrons2ndReco, eventWeight );
        hp_Nhfhad_vs_pt_barrel->Fill( pt2ndJetReco, nHFHadrons2ndReco, eventWeight );
        hp_Nhfem_vs_pt_barrel->Fill( pt2ndJetReco, nHFEM2ndReco, eventWeight );
    
        hp_Rch_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, eTracks2ndReco, eventWeight );
        hp_Rnh_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, eNeutralHadrons2ndReco, eventWeight );
        hp_Rgamma_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, ePhotons2ndReco, eventWeight );
        hp_Rmu_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, eMuons2ndReco, eventWeight );
        hp_Re_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, eElectrons2ndReco, eventWeight );
        hp_Rhfhad_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, eHFHadrons2ndReco, eventWeight );
        hp_Rhfem_vs_ptCorr_barrel->Fill( ptCorr2ndJetReco, eHFEM2ndReco, eventWeight );
    

        if( fabs(eta2ndJetGen)<1.4 ) {

          hp_RchGen_vs_pt_barrel->Fill( pt2ndJetGen, eTracks2ndGen, eventWeight );
          hp_RnhGen_vs_pt_barrel->Fill( pt2ndJetGen, eNeutralHadrons2ndGen, eventWeight );
          hp_RgammaGen_vs_pt_barrel->Fill( pt2ndJetGen, ePhotons2ndGen, eventWeight );
          hp_RmuGen_vs_pt_barrel->Fill( pt2ndJetGen, eMuons2ndGen, eventWeight );
          hp_ReGen_vs_pt_barrel->Fill( pt2ndJetGen, eElectrons2ndGen, eventWeight );
          hp_RhfhadGen_vs_pt_barrel->Fill( pt2ndJetGen, eHFHadrons2ndGen, eventWeight );
          hp_RhfemGen_vs_pt_barrel->Fill( pt2ndJetGen, eHFEM2ndGen, eventWeight );

          if( pt2ndJetGen>pt_thresh ) {
            hp_RchGen_vs_phi_barrel->Fill( phi2ndJetGen, eTracks2ndGen, eventWeight );
            hp_RnhGen_vs_phi_barrel->Fill( phi2ndJetGen, eNeutralHadrons2ndGen, eventWeight );
            hp_RgammaGen_vs_phi_barrel->Fill( phi2ndJetGen, ePhotons2ndGen, eventWeight );
            hp_RmuGen_vs_phi_barrel->Fill( phi2ndJetGen, eMuons2ndGen, eventWeight );
            hp_ReGen_vs_phi_barrel->Fill( phi2ndJetGen, eElectrons2ndGen, eventWeight );
            hp_RhfhadGen_vs_phi_barrel->Fill( phi2ndJetGen, eHFHadrons2ndGen, eventWeight );
            hp_RhfemGen_vs_phi_barrel->Fill( phi2ndJetGen, eHFEM2ndGen, eventWeight );
          }
        } //if barrel gen
      } //if barrel reco
    } //if jetID

    if( raw_corr=="raw" ) {
      if( pt2ndJetReco<pt_thresh ) continue;
    } else {
      if( (ptCorrJetReco< pt_thresh) || (ptCorr2ndJetReco<pt_thresh) ) continue;
    }


    if( useJetID ) {
      if( !pfJetIDOK || !pf2ndJetIDOK ) continue;
    }

    //loose back-to-back
    Float_t deltaPhi =  phiJetReco-phi2ndJetReco;
    if( deltaPhi >= TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
    if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();
    bool back2back = (fabs(deltaPhi)>(TMath::Pi()-0.5));

    if( back2back ) { //fill these before etamax cut:
      hp_Rch_vs_eta_stack->Fill( etaJetReco, eTracksReco, eventWeight );
      hp_Rch_vs_eta_stack->Fill( eta2ndJetReco, eTracks2ndReco, eventWeight );
      hp_Rnh_vs_eta_stack->Fill( etaJetReco, eNeutralHadronsReco, eventWeight );
      hp_Rnh_vs_eta_stack->Fill( eta2ndJetReco, eNeutralHadrons2ndReco, eventWeight );
      hp_Rgamma_vs_eta_stack->Fill( etaJetReco, ePhotonsReco, eventWeight );
      hp_Rgamma_vs_eta_stack->Fill( eta2ndJetReco, ePhotons2ndReco, eventWeight );
      hp_Rmu_vs_eta_stack->Fill( etaJetReco, eMuonsReco, eventWeight );
      hp_Rmu_vs_eta_stack->Fill( eta2ndJetReco, eMuons2ndReco, eventWeight );
      hp_Re_vs_eta_stack->Fill( etaJetReco, eElectronsReco, eventWeight );
      hp_Re_vs_eta_stack->Fill( eta2ndJetReco, eElectrons2ndReco, eventWeight );
      hp_Rhfhad_vs_eta_stack->Fill( etaJetReco, eHFHadronsReco, eventWeight );
      hp_Rhfhad_vs_eta_stack->Fill( eta2ndJetReco, eHFHadrons2ndReco, eventWeight );
      hp_Rhfem_vs_eta_stack->Fill( etaJetReco, eHFEMReco, eventWeight );
      hp_Rhfem_vs_eta_stack->Fill( eta2ndJetReco, eHFEM2ndReco, eventWeight );

    } 

    if( !(fabs(etaJetReco)<etaMax&&fabs(eta2ndJetReco)<etaMax) ) continue;

    //fill deltaPhi plots before deltaPhi cut:
    h1_deltaPhiJet->Fill( deltaPhi, eventWeight );
    if( !back2back ) continue;

    if( ( eNeutralHadronsReco*eJetReco > 300. && fabs(etaJetReco)<2.5 && fabs(etaJetReco)>1.4 ) || (eNeutralHadrons2ndReco*e2ndJetReco > 300. && fabs(eta2ndJetReco)<2.5 && fabs(eta2ndJetReco)>1.4 ) ) 
      std::cout << "Run: " << run << "\tEvent: " << event << "\tLS: " << lbn << std::endl;
   
    Float_t coin = rand->Uniform(1.);
    Float_t asymm = (ptJetReco-pt2ndJetReco)/(ptJetReco+pt2ndJetReco);
    if( coin <0.5 ) asymm*=-1.;
    h1_asymmJet->Fill( asymm, eventWeight );
  
    //Float_t ptAve = 0.5*(ptCorrJetReco+ptCorr2ndJetReco);

    h1_ptJet->Fill( ptJetReco, eventWeight );
    h1_ptJet->Fill( pt2ndJetReco, eventWeight );
    h1_ptCorrJet->Fill( ptCorrJetReco, eventWeight );
    h1_ptCorrJet->Fill( ptCorr2ndJetReco, eventWeight );

    if( fabs( etaJetReco ) < 2.4 ) {
      if( eTracksReco < 0.5 ) {
        h1_ptCorrJet_Rch050->Fill( ptCorrJetReco, eventWeight );
      } else if( eTracksReco < 0.7 ) {
        h1_ptCorrJet_Rch5070->Fill( ptCorrJetReco, eventWeight );
      } else {
        h1_ptCorrJet_Rch70100->Fill( ptCorrJetReco, eventWeight );
      }
    }

    if( fabs( eta2ndJetReco ) < 2.4 ) {
      if( eTracks2ndReco < 0.5 ) {
        h1_ptCorrJet_Rch050->Fill( ptCorr2ndJetReco, eventWeight );
      } else if( eTracks2ndReco < 0.7 ) {
        h1_ptCorrJet_Rch5070->Fill( ptCorr2ndJetReco, eventWeight );
      } else {
        h1_ptCorrJet_Rch70100->Fill( ptCorr2ndJetReco, eventWeight );
      }
    }


    hp_pt_vs_eta->Fill( etaJetReco, ptJetReco, eventWeight );
    hp_pt_vs_eta->Fill( eta2ndJetReco, pt2ndJetReco, eventWeight );
    hp_ptCorr_vs_eta->Fill( etaJetReco, ptCorrJetReco, eventWeight );
    hp_ptCorr_vs_eta->Fill( eta2ndJetReco, ptCorr2ndJetReco, eventWeight );
    

    Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
    Float_t pJetReco = ptJetReco/sin(thetaJetReco);
    Float_t pxJetReco = pJetReco*sin(thetaJetReco)*cos(phiJetReco);
    Float_t pyJetReco = pJetReco*sin(thetaJetReco)*sin(phiJetReco);
    Float_t pzJetReco = pJetReco*cos(thetaJetReco);
    Float_t massJet = sqrt( eJetReco*eJetReco - pJetReco*pJetReco );
    h1_massJet->Fill( massJet, eventWeight );
    h1_ptOverMJet->Fill( ptJetReco/massJet, eventWeight );
    h1_pOverMJet->Fill( pJetReco/massJet, eventWeight );

    Float_t theta2ndJetReco = 2.*atan(exp(-eta2ndJetReco));
    Float_t p2ndJetReco = pt2ndJetReco/sin(theta2ndJetReco);
    Float_t px2ndJetReco = p2ndJetReco*sin(theta2ndJetReco)*cos(phi2ndJetReco);
    Float_t py2ndJetReco = p2ndJetReco*sin(theta2ndJetReco)*sin(phi2ndJetReco);
    Float_t pz2ndJetReco = p2ndJetReco*cos(theta2ndJetReco);
    Float_t mass2ndJet = sqrt( e2ndJetReco*e2ndJetReco - p2ndJetReco*p2ndJetReco );
    h1_massJet->Fill( mass2ndJet, eventWeight );
    h1_ptOverMJet->Fill( pt2ndJetReco/mass2ndJet, eventWeight );
    h1_pOverMJet->Fill( p2ndJetReco/mass2ndJet, eventWeight );

    Float_t eSum = eJetReco+e2ndJetReco;
    Float_t pxSum = pxJetReco+px2ndJetReco;
    Float_t pySum = pyJetReco+py2ndJetReco;
    Float_t pzSum = pzJetReco+pz2ndJetReco;

    Float_t diJetMass = sqrt( eSum*eSum - pxSum*pxSum - pySum*pySum - pzSum*pzSum );
    h1_diJetMass->Fill( diJetMass, eventWeight ); 
  
    h1_etaJet->Fill( etaJetReco, eventWeight );
    h1_etaJet->Fill( eta2ndJetReco, eventWeight );
    h1_phiJet->Fill( phiJetReco, eventWeight );
    h1_phiJet->Fill( phi2ndJetReco, eventWeight );
    h1_nCandJet->Fill( N, eventWeight );
    h1_nCandJet->Fill( N2nd, eventWeight );
  
    h1_RchJet->Fill( eTracksReco, eventWeight );
    h1_RchJet->Fill( eTracks2ndReco, eventWeight );
    h1_RnhJet->Fill( eNeutralHadronsReco, eventWeight );
    h1_RnhJet->Fill( eNeutralHadrons2ndReco, eventWeight );
    h1_RgammaJet->Fill( ePhotonsReco, eventWeight );
    h1_RgammaJet->Fill( ePhotons2ndReco, eventWeight );
    h1_ReJet->Fill( eElectronsReco, eventWeight );
    h1_ReJet->Fill( eElectrons2ndReco, eventWeight );
    h1_RmuJet->Fill( eMuonsReco, eventWeight );
    h1_RmuJet->Fill( eMuons2ndReco, eventWeight );
    h1_RhfhadJet->Fill( eHFHadronsReco, eventWeight );
    h1_RhfhadJet->Fill( eHFHadrons2ndReco, eventWeight );
    h1_RhfemJet->Fill( eHFEMReco, eventWeight );
    h1_RhfemJet->Fill( eHFEM2ndReco, eventWeight );
  
    h1_EchJet->Fill( eTracksReco*eJetReco, eventWeight );
    h1_EchJet->Fill( eTracks2ndReco*e2ndJetReco, eventWeight );
    h1_EnhJet->Fill( eNeutralHadronsReco*eJetReco, eventWeight );
    h1_EnhJet->Fill( eNeutralHadrons2ndReco*e2ndJetReco, eventWeight );
    h1_EgammaJet->Fill( ePhotonsReco*eJetReco, eventWeight );
    h1_EgammaJet->Fill( ePhotons2ndReco*e2ndJetReco, eventWeight );
    h1_EeJet->Fill( eElectronsReco*eJetReco, eventWeight );
    h1_EeJet->Fill( eElectrons2ndReco*e2ndJetReco, eventWeight );
    h1_EmuJet->Fill( eMuonsReco*eJetReco, eventWeight );
    h1_EmuJet->Fill( eMuons2ndReco*e2ndJetReco, eventWeight );
    h1_EhfhadJet->Fill( eHFHadronsReco*eJetReco, eventWeight );
    h1_EhfhadJet->Fill( eHFHadrons2ndReco*e2ndJetReco, eventWeight );
    h1_EhfemJet->Fill( eHFEMReco*eJetReco, eventWeight );
    h1_EhfemJet->Fill( eHFEM2ndReco*e2ndJetReco, eventWeight );

    h1_PTchJet->Fill( ptTracksReco*ptJetReco, eventWeight );
    h1_PTchJet->Fill( ptTracks2ndReco*pt2ndJetReco, eventWeight );
    h1_PTnhJet->Fill( ptNeutralHadronsReco*ptJetReco, eventWeight );
    h1_PTnhJet->Fill( ptNeutralHadrons2ndReco*pt2ndJetReco, eventWeight );
    h1_PTgammaJet->Fill( ptPhotonsReco*ptJetReco, eventWeight );
    h1_PTgammaJet->Fill( ptPhotons2ndReco*pt2ndJetReco, eventWeight );
    h1_PTeJet->Fill( ptElectronsReco*ptJetReco, eventWeight );
    h1_PTeJet->Fill( ptElectrons2ndReco*pt2ndJetReco, eventWeight );
    h1_PTmuJet->Fill( ptMuonsReco*ptJetReco, eventWeight );
    h1_PTmuJet->Fill( ptMuons2ndReco*pt2ndJetReco, eventWeight );
    h1_PThfhadJet->Fill( ptHFHadronsReco*ptJetReco, eventWeight );
    h1_PThfhadJet->Fill( ptHFHadrons2ndReco*pt2ndJetReco, eventWeight );
    h1_PThfemJet->Fill( ptHFEMReco*ptJetReco, eventWeight );
    h1_PThfemJet->Fill( ptHFEM2ndReco*pt2ndJetReco, eventWeight );

    h1_NchJet->Fill( nTracksReco, eventWeight );
    h1_NchJet->Fill( nTracks2ndReco, eventWeight );
    h1_NnhJet->Fill( nNeutralHadronsReco, eventWeight );
    h1_NnhJet->Fill( nNeutralHadrons2ndReco, eventWeight );
    h1_NgammaJet->Fill( nPhotonsReco, eventWeight );
    h1_NgammaJet->Fill( nPhotons2ndReco , eventWeight);
    h1_NeJet->Fill( nElectronsReco, eventWeight );
    h1_NeJet->Fill( nElectrons2ndReco, eventWeight );
    h1_NmuJet->Fill( nMuonsReco, eventWeight );
    h1_NmuJet->Fill( nMuons2ndReco, eventWeight );
    h1_NhfhadJet->Fill( nHFHadronsReco, eventWeight );
    h1_NhfhadJet->Fill( nHFHadrons2ndReco, eventWeight );
    h1_NhfemJet->Fill( nHFEMReco, eventWeight );
    h1_NhfemJet->Fill( nHFEM2ndReco, eventWeight );
  
    if( fabs(etaJetReco)<1.4 ) {

      h1_ptCorrJet_barrel->Fill( ptCorrJetReco, eventWeight );
  
      h1_RchJet_barrel->Fill( eTracksReco, eventWeight );
      h1_RnhJet_barrel->Fill( eNeutralHadronsReco, eventWeight );
      h1_RgammaJet_barrel->Fill( ePhotonsReco, eventWeight );
      h1_ReJet_barrel->Fill( eElectronsReco, eventWeight );
      h1_RmuJet_barrel->Fill( eMuonsReco, eventWeight );
      h1_RhfhadJet_barrel->Fill( eHFHadronsReco, eventWeight );
      h1_RhfemJet_barrel->Fill( eHFEMReco, eventWeight );

      h1_EchJet_barrel->Fill( eTracksReco*eJetReco, eventWeight );
      h1_EnhJet_barrel->Fill( eNeutralHadronsReco*eJetReco, eventWeight );
      h1_EgammaJet_barrel->Fill( ePhotonsReco*eJetReco, eventWeight );
      h1_EeJet_barrel->Fill( eElectronsReco*eJetReco, eventWeight );
      h1_EmuJet_barrel->Fill( eMuonsReco*eJetReco, eventWeight );
      h1_EhfhadJet_barrel->Fill( eHFHadronsReco*eJetReco, eventWeight );
      h1_EhfemJet_barrel->Fill( eHFEMReco*eJetReco, eventWeight );

      h1_NchJet_barrel->Fill( nTracksReco, eventWeight );
      h1_NnhJet_barrel->Fill( nNeutralHadronsReco, eventWeight );
      h1_NgammaJet_barrel->Fill( nPhotonsReco, eventWeight );
      h1_NeJet_barrel->Fill( nElectronsReco, eventWeight );
      h1_NmuJet_barrel->Fill( nMuonsReco, eventWeight );
      h1_NhfhadJet_barrel->Fill( nHFHadronsReco, eventWeight );
      h1_NhfemJet_barrel->Fill( nHFEMReco, eventWeight );

      hp_Rch_vs_phi_barrel->Fill( phiJetReco, eTracksReco, eventWeight );
      hp_Rnh_vs_phi_barrel->Fill( phiJetReco, eNeutralHadronsReco, eventWeight );
      hp_Rgamma_vs_phi_barrel->Fill( phiJetReco, ePhotonsReco, eventWeight );
      hp_Rmu_vs_phi_barrel->Fill( phiJetReco, eMuonsReco, eventWeight );
      hp_Re_vs_phi_barrel->Fill( phiJetReco, eElectronsReco, eventWeight );
      hp_Rhfhad_vs_phi_barrel->Fill( phiJetReco, eHFHadronsReco, eventWeight );
      hp_Rhfem_vs_phi_barrel->Fill( phiJetReco, eHFEMReco, eventWeight );

    } else if( fabs(etaJetReco)<2.5 ) {
  
      h1_ptCorrJet_endcap->Fill( ptCorrJetReco, eventWeight );

      h1_RchJet_eta1425->Fill( eTracksReco, eventWeight );
      h1_RnhJet_eta1425->Fill( eNeutralHadronsReco, eventWeight );
      h1_RgammaJet_eta1425->Fill( ePhotonsReco, eventWeight );
      h1_ReJet_eta1425->Fill( eElectronsReco, eventWeight );
      h1_RmuJet_eta1425->Fill( eMuonsReco, eventWeight );
      h1_RhfhadJet_eta1425->Fill( eHFHadronsReco, eventWeight );
      h1_RhfemJet_eta1425->Fill( eHFEMReco, eventWeight );

      h1_EchJet_eta1425->Fill( eTracksReco*eJetReco, eventWeight );
      h1_EnhJet_eta1425->Fill( eNeutralHadronsReco*eJetReco, eventWeight );
      h1_EgammaJet_eta1425->Fill( ePhotonsReco*eJetReco, eventWeight );
      h1_EeJet_eta1425->Fill( eElectronsReco*eJetReco, eventWeight );
      h1_EmuJet_eta1425->Fill( eMuonsReco*eJetReco, eventWeight );
      h1_EhfhadJet_eta1425->Fill( eHFHadronsReco*eJetReco, eventWeight );
      h1_EhfemJet_eta1425->Fill( eHFEMReco*eJetReco, eventWeight );

      h1_NchJet_eta1425->Fill( nTracksReco, eventWeight );
      h1_NnhJet_eta1425->Fill( nNeutralHadronsReco, eventWeight );
      h1_NgammaJet_eta1425->Fill( nPhotonsReco, eventWeight );
      h1_NeJet_eta1425->Fill( nElectronsReco, eventWeight );
      h1_NmuJet_eta1425->Fill( nMuonsReco, eventWeight );
      h1_NhfhadJet_eta1425->Fill( nHFHadronsReco, eventWeight );
      h1_NhfemJet_eta1425->Fill( nHFEMReco, eventWeight );

    } //if endcaps

    if( fabs(etaJetReco)>1.4 && fabs(etaJetReco)<3. ) {
      h1_ptCorrJet_eta1430->Fill( ptCorrJetReco, eventWeight );
    }


    if( fabs(eta2ndJetReco)<1.4 ) {
  
      h1_ptCorrJet_barrel->Fill( ptCorr2ndJetReco, eventWeight );

      h1_RchJet_barrel->Fill( eTracks2ndReco, eventWeight );
      h1_RnhJet_barrel->Fill( eNeutralHadrons2ndReco, eventWeight );
      h1_RgammaJet_barrel->Fill( ePhotons2ndReco, eventWeight );
      h1_ReJet_barrel->Fill( eElectrons2ndReco, eventWeight );
      h1_RmuJet_barrel->Fill( eMuons2ndReco, eventWeight );
      h1_RhfhadJet_barrel->Fill( eHFHadrons2ndReco, eventWeight );
      h1_RhfemJet_barrel->Fill( eHFEM2ndReco, eventWeight );

      h1_EchJet_barrel->Fill( eTracks2ndReco*e2ndJetReco, eventWeight );
      h1_EnhJet_barrel->Fill( eNeutralHadrons2ndReco*e2ndJetReco, eventWeight );
      h1_EgammaJet_barrel->Fill( ePhotons2ndReco*e2ndJetReco, eventWeight );
      h1_EeJet_barrel->Fill( eElectrons2ndReco*e2ndJetReco, eventWeight );
      h1_EmuJet_barrel->Fill( eMuons2ndReco*e2ndJetReco, eventWeight );
      h1_EhfhadJet_barrel->Fill( eHFHadrons2ndReco*e2ndJetReco, eventWeight );
      h1_EhfemJet_barrel->Fill( eHFEM2ndReco*e2ndJetReco, eventWeight );

      h1_NchJet_barrel->Fill( nTracks2ndReco, eventWeight );
      h1_NnhJet_barrel->Fill( nNeutralHadrons2ndReco, eventWeight );
      h1_NgammaJet_barrel->Fill( nPhotons2ndReco, eventWeight );
      h1_NeJet_barrel->Fill( nElectrons2ndReco, eventWeight );
      h1_NmuJet_barrel->Fill( nMuons2ndReco, eventWeight );
      h1_NhfhadJet_barrel->Fill( nHFHadrons2ndReco , eventWeight);
      h1_NhfemJet_barrel->Fill( nHFEM2ndReco, eventWeight );

      hp_Rch_vs_phi_barrel->Fill( phi2ndJetReco, eTracks2ndReco, eventWeight );
      hp_Rnh_vs_phi_barrel->Fill( phi2ndJetReco, eNeutralHadrons2ndReco, eventWeight );
      hp_Rgamma_vs_phi_barrel->Fill( phi2ndJetReco, ePhotons2ndReco, eventWeight );
      hp_Rmu_vs_phi_barrel->Fill( phi2ndJetReco, eMuons2ndReco, eventWeight );
      hp_Re_vs_phi_barrel->Fill( phi2ndJetReco, eElectrons2ndReco, eventWeight );
      hp_Rhfhad_vs_phi_barrel->Fill( phi2ndJetReco, eHFHadrons2ndReco, eventWeight );
      hp_Rhfem_vs_phi_barrel->Fill( phi2ndJetReco, eHFEM2ndReco, eventWeight );

    } else if( fabs(eta2ndJetReco)<2.5 ) {
  
      h1_ptCorrJet_endcap->Fill( ptCorr2ndJetReco, eventWeight );

      h1_RchJet_eta1425->Fill( eTracks2ndReco, eventWeight );
      h1_RnhJet_eta1425->Fill( eNeutralHadrons2ndReco, eventWeight );
      h1_RgammaJet_eta1425->Fill( ePhotons2ndReco, eventWeight );
      h1_ReJet_eta1425->Fill( eElectrons2ndReco , eventWeight);
      h1_RmuJet_eta1425->Fill( eMuons2ndReco, eventWeight );
      h1_RhfhadJet_eta1425->Fill( eHFHadrons2ndReco, eventWeight );
      h1_RhfemJet_eta1425->Fill( eHFEM2ndReco, eventWeight );

      h1_EchJet_eta1425->Fill( eTracks2ndReco*e2ndJetReco, eventWeight );
      h1_EnhJet_eta1425->Fill( eNeutralHadrons2ndReco*e2ndJetReco, eventWeight );
      h1_EgammaJet_eta1425->Fill( ePhotons2ndReco*e2ndJetReco, eventWeight );
      h1_EeJet_eta1425->Fill( eElectrons2ndReco*e2ndJetReco, eventWeight );
      h1_EmuJet_eta1425->Fill( eMuons2ndReco*e2ndJetReco, eventWeight );
      h1_EhfhadJet_eta1425->Fill( eHFHadrons2ndReco*e2ndJetReco, eventWeight );
      h1_EhfemJet_eta1425->Fill( eHFEM2ndReco*e2ndJetReco, eventWeight );

      h1_NchJet_eta1425->Fill( nTracks2ndReco, eventWeight );
      h1_NnhJet_eta1425->Fill( nNeutralHadrons2ndReco, eventWeight );
      h1_NgammaJet_eta1425->Fill( nPhotons2ndReco, eventWeight );
      h1_NeJet_eta1425->Fill( nElectrons2ndReco, eventWeight );
      h1_NmuJet_eta1425->Fill( nMuons2ndReco, eventWeight );
      h1_NhfhadJet_eta1425->Fill( nHFHadrons2ndReco, eventWeight );
      h1_NhfemJet_eta1425->Fill( nHFEM2ndReco, eventWeight );

    } //if endcaps

    if( fabs(eta2ndJetReco)>1.4 && fabs(eta2ndJetReco)<3. ) {
      h1_ptCorrJet_eta1430->Fill( ptCorr2ndJetReco, eventWeight );
    }


    hp_Rch_vs_eta->Fill( etaJetReco, eTracksReco, eventWeight );
    hp_Rch_vs_eta->Fill( eta2ndJetReco, eTracks2ndReco, eventWeight );
    hp_Rnh_vs_eta->Fill( etaJetReco, eNeutralHadronsReco, eventWeight );
    hp_Rnh_vs_eta->Fill( eta2ndJetReco, eNeutralHadrons2ndReco, eventWeight );
    hp_Rgamma_vs_eta->Fill( etaJetReco, ePhotonsReco, eventWeight );
    hp_Rgamma_vs_eta->Fill( eta2ndJetReco, ePhotons2ndReco, eventWeight );
    hp_Rmu_vs_eta->Fill( etaJetReco, eMuonsReco, eventWeight );
    hp_Rmu_vs_eta->Fill( eta2ndJetReco, eMuons2ndReco, eventWeight );
    hp_Re_vs_eta->Fill( etaJetReco, eElectronsReco, eventWeight );
    hp_Re_vs_eta->Fill( eta2ndJetReco, eElectrons2ndReco, eventWeight );
    hp_Rhfhad_vs_eta->Fill( etaJetReco, eHFHadronsReco, eventWeight );
    hp_Rhfhad_vs_eta->Fill( eta2ndJetReco, eHFHadrons2ndReco, eventWeight );
    hp_Rhfem_vs_eta->Fill( etaJetReco, eHFEMReco, eventWeight );
    hp_Rhfem_vs_eta->Fill( eta2ndJetReco, eHFEM2ndReco, eventWeight );
  
    hp_Ech_vs_eta->Fill( etaJetReco, eTracksReco*eJetReco, eventWeight );
    hp_Ech_vs_eta->Fill( eta2ndJetReco, eTracks2ndReco*e2ndJetReco, eventWeight );
    hp_Enh_vs_eta->Fill( etaJetReco, eNeutralHadronsReco*eJetReco, eventWeight );
    hp_Enh_vs_eta->Fill( eta2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco, eventWeight );
    hp_Egamma_vs_eta->Fill( etaJetReco, ePhotonsReco*eJetReco, eventWeight );
    hp_Egamma_vs_eta->Fill( eta2ndJetReco, ePhotons2ndReco*e2ndJetReco, eventWeight );
    hp_Emu_vs_eta->Fill( etaJetReco, eMuonsReco*eJetReco, eventWeight );
    hp_Emu_vs_eta->Fill( eta2ndJetReco, eMuons2ndReco*e2ndJetReco, eventWeight );
    hp_Ee_vs_eta->Fill( etaJetReco, eElectronsReco*eJetReco, eventWeight );
    hp_Ee_vs_eta->Fill( eta2ndJetReco, eElectrons2ndReco*e2ndJetReco , eventWeight);
    hp_Ehfhad_vs_eta->Fill( etaJetReco, eHFHadronsReco*eJetReco, eventWeight );
    hp_Ehfhad_vs_eta->Fill( eta2ndJetReco, eHFHadrons2ndReco*e2ndJetReco, eventWeight );
    hp_Ehfem_vs_eta->Fill( etaJetReco, eHFEMReco*eJetReco, eventWeight );
    hp_Ehfem_vs_eta->Fill( eta2ndJetReco, eHFEM2ndReco*e2ndJetReco, eventWeight );

    hp_PTch_vs_eta->Fill( etaJetReco, ptTracksReco*ptJetReco, eventWeight );
    hp_PTch_vs_eta->Fill( eta2ndJetReco, ptTracks2ndReco*pt2ndJetReco, eventWeight );
    hp_PTnh_vs_eta->Fill( etaJetReco, ptNeutralHadronsReco*ptJetReco, eventWeight );
    hp_PTnh_vs_eta->Fill( eta2ndJetReco, ptNeutralHadrons2ndReco*pt2ndJetReco, eventWeight );
    hp_PTgamma_vs_eta->Fill( etaJetReco, ptPhotonsReco*ptJetReco, eventWeight );
    hp_PTgamma_vs_eta->Fill( eta2ndJetReco, ptPhotons2ndReco*pt2ndJetReco, eventWeight );
    hp_PTmu_vs_eta->Fill( etaJetReco, ptMuonsReco*ptJetReco, eventWeight );
    hp_PTmu_vs_eta->Fill( eta2ndJetReco, ptMuons2ndReco*pt2ndJetReco, eventWeight );
    hp_PTe_vs_eta->Fill( etaJetReco, ptElectronsReco*ptJetReco, eventWeight );
    hp_PTe_vs_eta->Fill( eta2ndJetReco, ptElectrons2ndReco*pt2ndJetReco, eventWeight );
    hp_PThfhad_vs_eta->Fill( etaJetReco, ptHFHadronsReco*ptJetReco, eventWeight );
    hp_PThfhad_vs_eta->Fill( eta2ndJetReco, ptHFHadrons2ndReco*pt2ndJetReco, eventWeight );
    hp_PThfem_vs_eta->Fill( etaJetReco, ptHFEMReco*ptJetReco, eventWeight );
    hp_PThfem_vs_eta->Fill( eta2ndJetReco, ptHFEM2ndReco*pt2ndJetReco , eventWeight);

    hp_Nch_vs_eta->Fill( etaJetReco, nTracksReco , eventWeight);
    hp_Nch_vs_eta->Fill( eta2ndJetReco, nTracks2ndReco , eventWeight);
    hp_Nnh_vs_eta->Fill( etaJetReco, nNeutralHadronsReco , eventWeight);
    hp_Nnh_vs_eta->Fill( eta2ndJetReco, nNeutralHadrons2ndReco, eventWeight );
    hp_Ngamma_vs_eta->Fill( etaJetReco, nPhotonsReco, eventWeight );
    hp_Ngamma_vs_eta->Fill( eta2ndJetReco, nPhotons2ndReco, eventWeight );
    hp_Nmu_vs_eta->Fill( etaJetReco, nMuonsReco, eventWeight );
    hp_Nmu_vs_eta->Fill( eta2ndJetReco, nMuons2ndReco, eventWeight );
    hp_Ne_vs_eta->Fill( etaJetReco, nElectronsReco, eventWeight );
    hp_Ne_vs_eta->Fill( eta2ndJetReco, nElectrons2ndReco, eventWeight );
    hp_Nhfhad_vs_eta->Fill( etaJetReco, nHFHadronsReco, eventWeight );
    hp_Nhfhad_vs_eta->Fill( eta2ndJetReco, nHFHadrons2ndReco, eventWeight );
    hp_Nhfem_vs_eta->Fill( etaJetReco, nHFEMReco, eventWeight );
    hp_Nhfem_vs_eta->Fill( eta2ndJetReco, nHFEM2ndReco, eventWeight );


    //correlations between components:
    hp_Egamma_vs_Ech->Fill( eTracksReco*eJetReco, ePhotonsReco*eJetReco, eventWeight );
    hp_Egamma_vs_Ech->Fill( eTracks2ndReco*e2ndJetReco, ePhotons2ndReco*e2ndJetReco, eventWeight );
    hp_Enh_vs_Ech->Fill( eTracksReco*eJetReco, eNeutralHadronsReco*eJetReco, eventWeight );
    hp_Enh_vs_Ech->Fill( eTracks2ndReco*e2ndJetReco, eNeutralHadrons2ndReco*e2ndJetReco, eventWeight );

    hp_Rgamma_vs_Rch->Fill( eTracksReco, ePhotonsReco, eventWeight );
    hp_Rgamma_vs_Rch->Fill( eTracks2ndReco, ePhotons2ndReco, eventWeight );
    hp_Rnh_vs_Rch->Fill( eTracksReco, eNeutralHadronsReco, eventWeight );
    hp_Rnh_vs_Rch->Fill( eTracks2ndReco, eNeutralHadrons2ndReco, eventWeight );
    hp_Rgamma_vs_Rnh->Fill( eNeutralHadronsReco, ePhotonsReco, eventWeight );
    hp_Rgamma_vs_Rnh->Fill( eNeutralHadrons2ndReco, ePhotons2ndReco, eventWeight );

    hp_Ngamma_vs_Nch->Fill( nTracksReco, nPhotonsReco, eventWeight );
    hp_Ngamma_vs_Nch->Fill( nTracks2ndReco, nPhotons2ndReco, eventWeight );
    hp_Nnh_vs_Nch->Fill( nTracksReco, nNeutralHadronsReco, eventWeight );
    hp_Nnh_vs_Nch->Fill( nTracks2ndReco, nNeutralHadrons2ndReco, eventWeight );


  } //for entries


  std::string outfileName_str;

  if(dataset!="") outfileName_str = "DiJet_"+dataset;
  else outfileName_str = "DiJet";

  outfileName_str += "_" +algoType;

  char outfileName[150];
  if( useJetID )
    sprintf( outfileName, "%s_%dGeV%s_eta%d.root", outfileName_str.c_str(), (Int_t)floor(pt_thresh), raw_corr.c_str(), (int)(10.*etaMax));
  else
    sprintf( outfileName, "%s_%dGeV%s_eta%d_noID.root", outfileName_str.c_str(), (Int_t)floor(pt_thresh), raw_corr.c_str(), (int)(10.*etaMax));

  TFile* outFile = new TFile(outfileName, "RECREATE");
  outFile->cd();


  h1_totalLumi->Write();

  h1_deltaPhiJet->Write();
  h1_asymmJet->Write();

  h1_massJet->Write();
  h1_ptOverMJet->Write();
  h1_pOverMJet->Write();
  h1_diJetMass->Write();

  h1_etaJet->Write();
  h1_phiJet->Write();
  h1_ptJet->Write();
  h1_ptCorrJet->Write();
  h1_ptCorrJet_barrel->Write();
  h1_ptCorrJet_endcap->Write();
  h1_ptCorrJet_eta1430->Write();
  h1_ptCorrJet_Rch050->Write();
  h1_ptCorrJet_Rch5070->Write();
  h1_ptCorrJet_Rch70100->Write();
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
  h1_EgammaJet->Write();
  h1_EeJet->Write();
  h1_EmuJet->Write();
  h1_EhfhadJet->Write();
  h1_EhfemJet->Write();

  h1_PTchJet->Write();
  h1_PTnhJet->Write();
  h1_PTgammaJet->Write();
  h1_PTeJet->Write();
  h1_PTmuJet->Write();
  h1_PThfhadJet->Write();
  h1_PThfemJet->Write();

  h1_EphotAveJet->Write();

  h1_NchJet->Write();
  h1_NnhJet->Write();
  h1_NgammaJet->Write();
  h1_NeJet->Write();
  h1_NmuJet->Write();
  h1_NhfhadJet->Write();
  h1_NhfemJet->Write();

  h1_RchJet_barrel->Write();
  h1_RnhJet_barrel->Write();
  h1_RgammaJet_barrel->Write();
  h1_ReJet_barrel->Write();
  h1_RmuJet_barrel->Write();
  h1_RhfhadJet_barrel->Write();
  h1_RhfemJet_barrel->Write();

  h1_EchJet_barrel->Write();
  h1_EnhJet_barrel->Write();
  h1_EgammaJet_barrel->Write();
  h1_EeJet_barrel->Write();
  h1_EmuJet_barrel->Write();
  h1_EhfhadJet_barrel->Write();
  h1_EhfemJet_barrel->Write();

  h1_NchJet_barrel->Write();
  h1_NnhJet_barrel->Write();
  h1_NgammaJet_barrel->Write();
  h1_NeJet_barrel->Write();
  h1_NmuJet_barrel->Write();
  h1_NhfhadJet_barrel->Write();
  h1_NhfemJet_barrel->Write();

  h1_RchJet_eta1425->Write();
  h1_RnhJet_eta1425->Write();
  h1_RgammaJet_eta1425->Write();
  h1_ReJet_eta1425->Write();
  h1_RmuJet_eta1425->Write();
  h1_RhfhadJet_eta1425->Write();
  h1_RhfemJet_eta1425->Write();

  h1_EchJet_eta1425->Write();
  h1_EnhJet_eta1425->Write();
  h1_EgammaJet_eta1425->Write();
  h1_EeJet_eta1425->Write();
  h1_EmuJet_eta1425->Write();
  h1_EhfhadJet_eta1425->Write();
  h1_EhfemJet_eta1425->Write();

  h1_NchJet_eta1425->Write();
  h1_NnhJet_eta1425->Write();
  h1_NgammaJet_eta1425->Write();
  h1_NeJet_eta1425->Write();
  h1_NmuJet_eta1425->Write();
  h1_NhfhadJet_eta1425->Write();
  h1_NhfemJet_eta1425->Write();

  hp_pt_vs_eta->Write();
  hp_ptCorr_vs_eta->Write();

  hp_Rch_vs_eta_stack->Write();
  hp_Rnh_vs_eta_stack->Write();
  hp_Rgamma_vs_eta_stack->Write();
  hp_Rmu_vs_eta_stack->Write();
  hp_Re_vs_eta_stack->Write();
  hp_Rhfhad_vs_eta_stack->Write();
  hp_Rhfem_vs_eta_stack->Write();

  hp_RchGen_vs_eta_stack->Write();
  hp_RnhGen_vs_eta_stack->Write();
  hp_RgammaGen_vs_eta_stack->Write();
  hp_RmuGen_vs_eta_stack->Write();
  hp_ReGen_vs_eta_stack->Write();
  hp_RhfhadGen_vs_eta_stack->Write();
  hp_RhfemGen_vs_eta_stack->Write();

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

  hp_PTch_vs_eta->Write();
  hp_PTnh_vs_eta->Write();
  hp_PTgamma_vs_eta->Write();
  hp_PTmu_vs_eta->Write();
  hp_PTe_vs_eta->Write();
  hp_PThfhad_vs_eta->Write();
  hp_PThfem_vs_eta->Write();

  hp_Nch_vs_eta->Write();
  hp_Nnh_vs_eta->Write();
  hp_Ngamma_vs_eta->Write();
  hp_Nmu_vs_eta->Write();
  hp_Ne_vs_eta->Write();
  hp_Nhfhad_vs_eta->Write();
  hp_Nhfem_vs_eta->Write();

  hp_Rch_vs_pt_barrel->Write();
  hp_Rnh_vs_pt_barrel->Write();
  hp_Rgamma_vs_pt_barrel->Write();
  hp_Rmu_vs_pt_barrel->Write();
  hp_Re_vs_pt_barrel->Write();
  hp_Rhfhad_vs_pt_barrel->Write();
  hp_Rhfem_vs_pt_barrel->Write();

  hp_Ech_vs_pt_barrel->Write();
  hp_Enh_vs_pt_barrel->Write();
  hp_Egamma_vs_pt_barrel->Write();
  hp_Emu_vs_pt_barrel->Write();
  hp_Ee_vs_pt_barrel->Write();
  hp_Ehfhad_vs_pt_barrel->Write();
  hp_Ehfem_vs_pt_barrel->Write();

  hp_Nch_vs_pt_barrel->Write();
  hp_Nnh_vs_pt_barrel->Write();
  hp_Ngamma_vs_pt_barrel->Write();
  hp_Nmu_vs_pt_barrel->Write();
  hp_Ne_vs_pt_barrel->Write();
  hp_Nhfhad_vs_pt_barrel->Write();
  hp_Nhfem_vs_pt_barrel->Write();

  hp_Rch_vs_ptCorr_barrel->Write();
  hp_Rnh_vs_ptCorr_barrel->Write();
  hp_Rgamma_vs_ptCorr_barrel->Write();
  hp_Rmu_vs_ptCorr_barrel->Write();
  hp_Re_vs_ptCorr_barrel->Write();
  hp_Rhfhad_vs_ptCorr_barrel->Write();
  hp_Rhfem_vs_ptCorr_barrel->Write();

//hp_Ech_vs_ptCorr_barrel->Write();
//hp_Enh_vs_ptCorr_barrel->Write();
//hp_Egamma_vs_ptCorr_barrel->Write();
//hp_Emu_vs_ptCorr_barrel->Write();
//hp_Ee_vs_ptCorr_barrel->Write();
//hp_Ehfhad_vs_ptCorr_barrel->Write();
//hp_Ehfem_vs_ptCorr_barrel->Write();

//hp_Nch_vs_ptCorr_barrel->Write();
//hp_Nnh_vs_ptCorr_barrel->Write();
//hp_Ngamma_vs_ptCorr_barrel->Write();
//hp_Nmu_vs_ptCorr_barrel->Write();
//hp_Ne_vs_ptCorr_barrel->Write();
//hp_Nhfhad_vs_ptCorr_barrel->Write();
//hp_Nhfem_vs_ptCorr_barrel->Write();

  hp_Rch_vs_phi_barrel->Write();
  hp_Rnh_vs_phi_barrel->Write();
  hp_Rgamma_vs_phi_barrel->Write();
  hp_Rmu_vs_phi_barrel->Write();
  hp_Re_vs_phi_barrel->Write();
  hp_Rhfhad_vs_phi_barrel->Write();
  hp_Rhfem_vs_phi_barrel->Write();

  hp_RchGen_vs_phi_barrel->Write();
  hp_RnhGen_vs_phi_barrel->Write();
  hp_RgammaGen_vs_phi_barrel->Write();
  hp_RmuGen_vs_phi_barrel->Write();
  hp_ReGen_vs_phi_barrel->Write();
  hp_RhfhadGen_vs_phi_barrel->Write();
  hp_RhfemGen_vs_phi_barrel->Write();

  hp_Egamma_vs_Ech->Write();
  hp_Enh_vs_Ech->Write();

  hp_Rgamma_vs_Rch->Write();
  hp_Rnh_vs_Rch->Write();
  hp_Rgamma_vs_Rnh->Write();

  hp_Ngamma_vs_Nch->Write();
  hp_Nnh_vs_Nch->Write();

  outFile->Close();


}
  

bool PFJetID( Float_t Rch, Float_t Rnh, Float_t Rgamma, Int_t N, Int_t Nch, Float_t eta) {

  bool pass = false;

  if( fabs(eta)<2.4 ) {
    if( Rch>0. ) pass = true;
  } else {
    if( Rnh<1. && Rgamma<1. ) pass=true;
  }

  return pass;

}
