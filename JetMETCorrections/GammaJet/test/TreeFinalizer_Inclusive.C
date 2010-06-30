#include <TH1F.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "fitTools.C"




bool PFJetID( Float_t Rch, Float_t Rnh, Float_t Rgamma, Int_t N, Int_t Nch, Float_t eta);



void finalize(const std::string& dataset, std::string algoType, Float_t pt_thresh=5., std::string raw_corr="raw", bool useJetID=true ) {
  
  Float_t etaMax = 3.;

  std::string infileName, treeName;
  TFile* infile=0;
  TH1F* h1_lumi=0;
  Double_t totalLumi=0.;

  TChain* tree = new TChain("jetTree");

  if( dataset=="DATA_7TeV_20100407") {

    infileName = "Inclusive_2ndLevelTree_DATA_7TeV_GOODCOLL_Apr1stSkim_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "Inclusive_2ndLevelTree_DATA_MinimumBias_Commissioning10-GOODCOLL-v8_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "Inclusive_2ndLevelTree_Commissioning10-GOODCOLL-v8_runs132605_606_" + algoType + ".root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="MinBias_Spring10-START3X_V26A_356ReReco-v1" ) {

    infileName = "Inclusive_2ndLevelTree_" + dataset + "_" + algoType + "_0.root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = "Inclusive_2ndLevelTree_" + dataset + "_" + algoType + "_1.root/jetTree";
    tree->Add(infileName.c_str());
    std::cout << "-> Added " << infileName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else {

    infileName = "files_Inclusive_2ndLevel_" + dataset+".txt";

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

    } //if file

  }

  
  std::cout << "Done adding." << std::endl;

  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);

  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Int_t nJetReco;
  tree->SetBranchAddress("nJetReco", &nJetReco);
  Float_t eJetReco[40];
  tree->SetBranchAddress("eJetReco", eJetReco);
  Float_t ptJetReco[40];
  tree->SetBranchAddress("ptJetReco", ptJetReco);
  Float_t ptCorrJetReco[40];
  tree->SetBranchAddress("ptCorrJetReco", ptCorrJetReco);
  Float_t etaJetReco[40];
  tree->SetBranchAddress("etaJetReco", etaJetReco);
  Float_t phiJetReco[40];
  tree->SetBranchAddress("phiJetReco", phiJetReco);

  Float_t eTracksReco[40];
  tree->SetBranchAddress("eTracksReco", eTracksReco);
  Float_t ePhotonsReco[40];
  tree->SetBranchAddress("ePhotonsReco", ePhotonsReco);
  Float_t eNeutralHadronsReco[40];
  tree->SetBranchAddress("eNeutralHadronsReco", eNeutralHadronsReco);
  Float_t eMuonsReco[40];
  tree->SetBranchAddress("eMuonsReco", eMuonsReco);
  Float_t eElectronsReco[40];
  tree->SetBranchAddress("eElectronsReco", eElectronsReco);
  Float_t eHFHadronsReco[40];
  tree->SetBranchAddress("eHFHadronsReco", eHFHadronsReco);
  Float_t eHFEMReco[40];
  tree->SetBranchAddress("eHFEMReco", eHFEMReco);

  Float_t ptTracksReco[40];
  tree->SetBranchAddress("ptTracksReco", ptTracksReco);
  Float_t ptPhotonsReco[40];
  tree->SetBranchAddress("ptPhotonsReco", ptPhotonsReco);
  Float_t ptNeutralHadronsReco[40];
  tree->SetBranchAddress("ptNeutralHadronsReco", ptNeutralHadronsReco);
  Float_t ptMuonsReco[40];
  tree->SetBranchAddress("ptMuonsReco", ptMuonsReco);
  Float_t ptElectronsReco[40];
  tree->SetBranchAddress("ptElectronsReco", ptElectronsReco);
  Float_t ptHFHadronsReco[40];
  tree->SetBranchAddress("ptHFHadronsReco", ptHFHadronsReco);
  Float_t ptHFEMReco[40];
  tree->SetBranchAddress("ptHFEMReco", ptHFEMReco);

  Int_t nTracksReco[40];
  tree->SetBranchAddress("nTracksReco", nTracksReco);
  Int_t nPhotonsReco[40];
  tree->SetBranchAddress("nPhotonsReco", nPhotonsReco);
  Int_t nNeutralHadronsReco[40];
  tree->SetBranchAddress("nNeutralHadronsReco", nNeutralHadronsReco);
  Int_t nElectronsReco[40];
  tree->SetBranchAddress("nElectronsReco", nElectronsReco);
  Int_t nMuonsReco[40];
  tree->SetBranchAddress("nMuonsReco", nMuonsReco);
  Int_t nHFHadronsReco[40];
  tree->SetBranchAddress("nHFHadronsReco", nHFHadronsReco);
  Int_t nHFEMReco[40];
  tree->SetBranchAddress("nHFEMReco", nHFEMReco);

  Float_t eJetGen[40];
  tree->SetBranchAddress("eJetGen", eJetGen);
  Float_t ptJetGen[40];
  tree->SetBranchAddress("ptJetGen", ptJetGen);
  Float_t etaJetGen[40];
  tree->SetBranchAddress("etaJetGen", etaJetGen);
  Float_t phiJetGen[40];
  tree->SetBranchAddress("phiJetGen", phiJetGen);

  Float_t eTracksGen[40];
  tree->SetBranchAddress("eTracksGen", eTracksGen);
  Float_t ePhotonsGen[40];
  tree->SetBranchAddress("ePhotonsGen", ePhotonsGen);
  Float_t eNeutralHadronsGen[40];
  tree->SetBranchAddress("eNeutralHadronsGen", eNeutralHadronsGen);
  Float_t eMuonsGen[40];
  tree->SetBranchAddress("eMuonsGen", eMuonsGen);
  Float_t eElectronsGen[40];
  tree->SetBranchAddress("eElectronsGen", eElectronsGen);
  Float_t eHFHadronsGen[40];
  tree->SetBranchAddress("eHFHadronsGen", eHFHadronsGen);
  Float_t eHFEMGen[40];
  tree->SetBranchAddress("eHFEMGen", eHFEMGen);

  Int_t nBins_eta = 40;
  Int_t nBins_E = 50;
  Int_t nBins_R = 30;
  
  Float_t eMax = 700.;
  Float_t ptMax = 500.;

  TH1F* h1_nJets = new TH1F("nJets", "", 11, -0.5, 10.5);

  TH1F* h1_massJet = new TH1F("massJet", "", 80, 0., 40.);
  TH1F* h1_ptOverMJet = new TH1F("ptOverMJet", "", 80, 0., 40.);
  TH1F* h1_pOverMJet = new TH1F("pOverMJet", "", 80, 0., 40.);
  TH1F* h1_MoEJet = new TH1F("MoEJet", "", 15, 0., 0.6);
  TH1F* h1_nCandJet = new TH1F("nCandJet", "", 26, -0.5, 25.);

  TH1F* h1_etaJet = new TH1F("etaJet", "", nBins_eta, -etaMax, etaMax);
  TH1F* h1_phiJet = new TH1F("phiJet", "", nBins_eta, -3.5, 3.5);
  TH1F* h1_ptJet = new TH1F("ptJet", "", 50, 0., ptMax);
  TH1F* h1_ptCorrJet = new TH1F("ptCorrJet", "", 50, 0., ptMax);

  TH1F* h1_RchJet = new TH1F("RchJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet = new TH1F("RnhJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet = new TH1F("RgammaJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet = new TH1F("RmuJet", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet = new TH1F("ReJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet = new TH1F("RhfhadJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet = new TH1F("RhfemJet", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet = new TH1F("EchJet", "", nBins_E, 0., eMax);
  TH1F* h1_EnhJet = new TH1F("EnhJet", "", nBins_E, 0., eMax);
  TH1F* h1_EnhJet_endcap = new TH1F("EnhJet_endcap", "", nBins_E, 0., eMax);
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

  TH1F* h1_NchJet = new TH1F("NchJet", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet = new TH1F("NnhJet", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet = new TH1F("NgammaJet", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet = new TH1F("NeJet", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet = new TH1F("NmuJet", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet = new TH1F("NhfhadJet", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet = new TH1F("NhfemJet", "", 21, -0.5, 20.5);

  TH1F* h1_RchJet_eta02 = new TH1F("RchJet_eta02", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet_eta02 = new TH1F("RnhJet_eta02", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet_eta02 = new TH1F("RgammaJet_eta02", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet_eta02 = new TH1F("RmuJet_eta02", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet_eta02 = new TH1F("ReJet_eta02", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet_eta02 = new TH1F("RhfhadJet_eta02", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet_eta02 = new TH1F("RhfemJet_eta02", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet_eta02 = new TH1F("EchJet_eta02", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EnhJet_eta02 = new TH1F("EnhJet_eta02", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EgammaJet_eta02 = new TH1F("EgammaJet_eta02", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EeJet_eta02 = new TH1F("EeJet_eta02", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EmuJet_eta02 = new TH1F("EmuJet_eta02", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EhfhadJet_eta02 = new TH1F("EhfhadJet_eta02", "", nBins_E, 0., eMax/2.);
  TH1F* h1_EhfemJet_eta02 = new TH1F("EhfemJet_eta02", "", nBins_E, 0., eMax/2.);

  TH1F* h1_NchJet_eta02 = new TH1F("NchJet_eta02", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet_eta02 = new TH1F("NnhJet_eta02", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet_eta02 = new TH1F("NgammaJet_eta02", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet_eta02 = new TH1F("NeJet_eta02", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet_eta02 = new TH1F("NmuJet_eta02", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet_eta02 = new TH1F("NhfhadJet_eta02", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet_eta02 = new TH1F("NhfemJet_eta02", "", 21, -0.5, 20.5);

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

  TH1F* h1_NchJet_barrel = new TH1F("NchJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet_barrel = new TH1F("NnhJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet_barrel = new TH1F("NgammaJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet_barrel = new TH1F("NeJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet_barrel = new TH1F("NmuJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet_barrel = new TH1F("NhfhadJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet_barrel = new TH1F("NhfemJet_barrel", "", 21, -0.5, 20.5);

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

  TH1F* h1_NchJet_eta1425 = new TH1F("NchJet_eta1425", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet_eta1425 = new TH1F("NnhJet_eta1425", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet_eta1425 = new TH1F("NgammaJet_eta1425", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet_eta1425 = new TH1F("NeJet_eta1425", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet_eta1425 = new TH1F("NmuJet_eta1425", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet_eta1425 = new TH1F("NhfhadJet_eta1425", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet_eta1425 = new TH1F("NhfemJet_eta1425", "", 21, -0.5, 20.5);

  TH1F* h1_RchJet_eta163 = new TH1F("RchJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet_eta163 = new TH1F("RnhJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet_eta163 = new TH1F("RgammaJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet_eta163 = new TH1F("RmuJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet_eta163 = new TH1F("ReJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet_eta163 = new TH1F("RhfhadJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet_eta163 = new TH1F("RhfemJet_eta163", "", nBins_R, 0., 1.1);

  TH1F* h1_EgammanhJet        = new TH1F("EgammanhJet", "", nBins_E, 0., eMax);
  TH1F* h1_EgammanhJet_barrel = new TH1F("EgammanhJet_barrel", "", nBins_E, 0., eMax);
  TH1F* h1_EgammanhJet_eta163 = new TH1F("EgammanhJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EgammanhJet_eta1425 = new TH1F("EgammanhJet_eta1425", "", nBins_E, 0., eMax);

  TH1F* h1_RgammanhJet        = new TH1F("RgammanhJet", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammanhJet_barrel = new TH1F("RgammanhJet_barrel", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammanhJet_eta163 = new TH1F("RgammanhJet_eta163", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammanhJet_eta1425 = new TH1F("RgammanhJet_eta1425", "", nBins_R, 0., 1.1);

  TH1F* h1_NgammanhJet        = new TH1F("NgammanhJet", "", 21, -0.5, 20.5);
  TH1F* h1_NgammanhJet_barrel = new TH1F("NgammanhJet_barrel", "", 21, -0.5, 20.5);
  TH1F* h1_NgammanhJet_eta163 = new TH1F("NgammanhJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NgammanhJet_eta1425 = new TH1F("NgammanhJet_eta1425", "", 21, -0.5, 20.5);

  TH1F* h1_EchJet_eta163 = new TH1F("EchJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EnhJet_eta163 = new TH1F("EnhJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EgammaJet_eta163 = new TH1F("EgammaJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EeJet_eta163 = new TH1F("EeJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EmuJet_eta163 = new TH1F("EmuJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EhfhadJet_eta163 = new TH1F("EhfhadJet_eta163", "", nBins_E, 0., eMax);
  TH1F* h1_EhfemJet_eta163 = new TH1F("EhfemJet_eta163", "", nBins_E, 0., eMax);

  TH1F* h1_NchJet_eta163 = new TH1F("NchJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet_eta163 = new TH1F("NnhJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet_eta163 = new TH1F("NgammaJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet_eta163 = new TH1F("NeJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet_eta163 = new TH1F("NmuJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet_eta163 = new TH1F("NhfhadJet_eta163", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet_eta163 = new TH1F("NhfemJet_eta163", "", 21, -0.5, 20.5);

  TH1F* h1_RchJet_eta23 = new TH1F("RchJet_eta23", "", nBins_R, 0., 1.1);
  TH1F* h1_RnhJet_eta23 = new TH1F("RnhJet_eta23", "", nBins_R, 0., 1.1);
  TH1F* h1_RgammaJet_eta23 = new TH1F("RgammaJet_eta23", "", nBins_R, 0., 1.1);
  TH1F* h1_RmuJet_eta23 = new TH1F("RmuJet_eta23", "", nBins_R, 0., 1.1);
  TH1F* h1_ReJet_eta23 = new TH1F("ReJet_eta23", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfhadJet_eta23 = new TH1F("RhfhadJet_eta23", "", nBins_R, 0., 1.1);
  TH1F* h1_RhfemJet_eta23 = new TH1F("RhfemJet_eta23", "", nBins_R, 0., 1.1);

  TH1F* h1_EchJet_eta23 = new TH1F("EchJet_eta23", "", nBins_E, 0., eMax);
  TH1F* h1_EnhJet_eta23 = new TH1F("EnhJet_eta23", "", nBins_E, 0., eMax);
  TH1F* h1_EgammaJet_eta23 = new TH1F("EgammaJet_eta23", "", nBins_E, 0., eMax);
  TH1F* h1_EeJet_eta23 = new TH1F("EeJet_eta23", "", nBins_E, 0., eMax);
  TH1F* h1_EmuJet_eta23 = new TH1F("EmuJet_eta23", "", nBins_E, 0., eMax);
  TH1F* h1_EhfhadJet_eta23 = new TH1F("EhfhadJet_eta23", "", nBins_E, 0., eMax);
  TH1F* h1_EhfemJet_eta23 = new TH1F("EhfemJet_eta23", "", nBins_E, 0., eMax);

  TH1F* h1_NchJet_eta23 = new TH1F("NchJet_eta23", "", 21, -0.5, 20.5);
  TH1F* h1_NnhJet_eta23 = new TH1F("NnhJet_eta23", "", 21, -0.5, 20.5);
  TH1F* h1_NgammaJet_eta23 = new TH1F("NgammaJet_eta23", "", 21, -0.5, 20.5);
  TH1F* h1_NeJet_eta23 = new TH1F("NeJet_eta23", "", 21, -0.5, 20.5);
  TH1F* h1_NmuJet_eta23 = new TH1F("NmuJet_eta23", "", 21, -0.5, 20.5);
  TH1F* h1_NhfhadJet_eta23 = new TH1F("NhfhadJet_eta23", "", 21, -0.5, 20.5);
  TH1F* h1_NhfemJet_eta23 = new TH1F("NhfemJet_eta23", "", 21, -0.5, 20.5);

  Float_t etaMax_profiles = etaMax;

  TProfile* hp_pt_vs_eta = new TProfile("pt_vs_eta", "pt vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles,  0., 300.);
  TProfile* hp_ptCorr_vs_eta = new TProfile("ptCorr_vs_eta", "ptCorr vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles,  0., 300.);

  TProfile* hp_Rch_vs_eta = new TProfile("Rch_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles,  0., 1000.);
  TProfile* hp_Rnh_vs_eta = new TProfile("Rnh_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Rgamma_vs_eta = new TProfile("Rgamma_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Rmu_vs_eta = new TProfile("Rmu_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Re_vs_eta = new TProfile("Re_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Rhfhad_vs_eta = new TProfile("Rhfhad_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Rhfem_vs_eta = new TProfile("Rhfem_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);

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

  TProfile* hp_RchGen_vs_eta_stack = new TProfile("RchGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);
  TProfile* hp_RnhGen_vs_eta_stack = new TProfile("RnhGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);
  TProfile* hp_RgammaGen_vs_eta_stack = new TProfile("RgammaGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);
  TProfile* hp_RmuGen_vs_eta_stack = new TProfile("RmuGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);
  TProfile* hp_ReGen_vs_eta_stack = new TProfile("ReGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);
  TProfile* hp_RhfhadGen_vs_eta_stack = new TProfile("RhfhadGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);
  TProfile* hp_RhfemGen_vs_eta_stack = new TProfile("RhfemGen_vs_eta_stack", "R_{ch} vs. eta", nBins_stack_eta-1, Lower_stack_eta, 0., 1000.);

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

  TProfile* hp_Ech_vs_eta = new TProfile("Ech_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Enh_vs_eta = new TProfile("Enh_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Egamma_vs_eta = new TProfile("Egamma_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Ee_vs_eta = new TProfile("Ee_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Emu_vs_eta = new TProfile("Emu_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Ehfhad_vs_eta = new TProfile("Ehfhad_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Ehfem_vs_eta = new TProfile("Ehfem_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);

  TProfile* hp_PTch_vs_eta = new TProfile("PTch_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_PTnh_vs_eta = new TProfile("PTnh_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_PTgamma_vs_eta = new TProfile("PTgamma_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_PTe_vs_eta = new TProfile("PTe_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_PTmu_vs_eta = new TProfile("PTmu_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_PThfhad_vs_eta = new TProfile("PThfhad_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_PThfem_vs_eta = new TProfile("PThfem_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);

  TProfile* hp_Nch_vs_eta = new TProfile("Nch_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Nnh_vs_eta = new TProfile("Nnh_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Ngamma_vs_eta = new TProfile("Ngamma_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Ne_vs_eta = new TProfile("Ne_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Nmu_vs_eta = new TProfile("Nmu_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Nhfhad_vs_eta = new TProfile("Nhfhad_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Nhfem_vs_eta = new TProfile("Nhfem_vs_eta", "N_{ch} vs. eta", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);

  TProfile* hp_Egammanh_vs_eta = new TProfile("Egammanh_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Rgammanh_vs_eta = new TProfile("Rgammanh_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);
  TProfile* hp_Ngammanh_vs_eta = new TProfile("Ngammanh_vs_eta", "", nBins_eta, -etaMax_profiles, etaMax_profiles, 0., 1000.);

  Int_t nBins_pt_lower = fitTools::getNbins_stack("pt");;
  Double_t Lower[nBins_pt_lower];
  fitTools::getBins_stack( nBins_pt_lower, Lower, "pt" );

  TProfile* hp_Rch_vs_pt = new TProfile("Rch_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rnh_vs_pt = new TProfile("Rnh_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rgamma_vs_pt = new TProfile("Rgamma_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rmu_vs_pt = new TProfile("Rmu_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Re_vs_pt = new TProfile("Re_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfhad_vs_pt = new TProfile("Rhfhad_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfem_vs_pt = new TProfile("Rhfem_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Ech_vs_pt = new TProfile("Ech_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Enh_vs_pt = new TProfile("Enh_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Egamma_vs_pt = new TProfile("Egamma_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ee_vs_pt = new TProfile("Ee_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Emu_vs_pt = new TProfile("Emu_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ehfhad_vs_pt = new TProfile("Ehfhad_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ehfem_vs_pt = new TProfile("Ehfem_vs_pt", "R_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Nch_vs_pt = new TProfile("Nch_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nnh_vs_pt = new TProfile("Nnh_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ngamma_vs_pt = new TProfile("Ngamma_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Ne_vs_pt = new TProfile("Ne_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nmu_vs_pt = new TProfile("Nmu_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nhfhad_vs_pt = new TProfile("Nhfhad_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Nhfem_vs_pt = new TProfile("Nhfem_vs_pt", "N_{ch} vs. pt", nBins_pt_lower-1, Lower, 0., 1000.);

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

  TProfile* hp_Rch_vs_ptCorr = new TProfile("Rch_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rnh_vs_ptCorr = new TProfile("Rnh_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rgamma_vs_ptCorr = new TProfile("Rgamma_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rmu_vs_ptCorr = new TProfile("Rmu_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Re_vs_ptCorr = new TProfile("Re_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfhad_vs_ptCorr = new TProfile("Rhfhad_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfem_vs_ptCorr = new TProfile("Rhfem_vs_ptCorr", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Rch_vs_ptCorr_barrel = new TProfile("Rch_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rnh_vs_ptCorr_barrel = new TProfile("Rnh_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rgamma_vs_ptCorr_barrel = new TProfile("Rgamma_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rmu_vs_ptCorr_barrel = new TProfile("Rmu_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Re_vs_ptCorr_barrel = new TProfile("Re_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfhad_vs_ptCorr_barrel = new TProfile("Rhfhad_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);
  TProfile* hp_Rhfem_vs_ptCorr_barrel = new TProfile("Rhfem_vs_ptCorr_barrel", "R_{ch} vs. ptCorr", nBins_pt_lower-1, Lower, 0., 1000.);

  TProfile* hp_Egamma_vs_Ech = new TProfile("Egamma_vs_Ech", "", nBins_E, 0., 15., 0., 1000.);
  TProfile* hp_Enh_vs_Ech = new TProfile("Enh_vs_Ech", "", nBins_E, 0., 15., 0., 1000.);

  TProfile* hp_Rgamma_vs_Rch = new TProfile("Rgamma_vs_Rch", "", nBins_R, 0., 1.1, 0., 1.1);
  TProfile* hp_Rnh_vs_Rch = new TProfile("Rnh_vs_Rch", "", nBins_R, 0., 1.1, 0., 1.1);
  TProfile* hp_Rnh_vs_Rgamma = new TProfile("Rnh_vs_Rgamma", "", nBins_R, 0., 1.1, 0., 1.1);

//TH2D* h2_Rgamma_vs_Rch = new TH2D("Rgamma_vs_Rch", "", nBins_R, 0., 1.1, nBins_R, 0., 1.1);
//TH2D* h2_Rnh_vs_Rch    = new TH2D("Rnh_vs_Rch", "", nBins_R, 0., 1.1, nBins_R, 0., 1.1);
//TH2D* h2_Rnh_vs_Rgamma = new TH2D("Rnh_vs_Rgamma", "", nBins_R, 0., 1.1, nBins_R, 0., 1.1);

  TProfile* hp_Ngamma_vs_Nch = new TProfile("Ngamma_vs_Nch", "", 21, -0.5, 20.5, -0.5, 20.5);
  TProfile* hp_Nnh_vs_Nch = new TProfile("Nnh_vs_Nch", "", 21, -0.5, 20.5, -0.5, 20.5);

  TRandom* rand = new TRandom();

  int nEntries = tree->GetEntries();

  std::cout << "-> About to start loop. Total entries: " << nEntries << std::endl;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    eventWeight = 1.;

    h1_nJets->Fill(nJetReco);
 
    for( unsigned iJetReco=0; iJetReco<nJetReco; ++iJetReco ) {

      if( iJetReco >= 40 ) continue;

      //PFJetID:
      Int_t N = nTracksReco[iJetReco]+nPhotonsReco[iJetReco]+nNeutralHadronsReco[iJetReco]+nElectronsReco[iJetReco]+nMuonsReco[iJetReco]+nHFHadronsReco[iJetReco]+nHFEMReco[iJetReco];
      Int_t Nch = nTracksReco[iJetReco]+nElectronsReco[iJetReco]+nMuonsReco[iJetReco];
      bool pfJetIDOK = PFJetID( eTracksReco[iJetReco], eNeutralHadronsReco[iJetReco], ePhotonsReco[iJetReco], N, Nch, etaJetReco[iJetReco]);

      //fill pt stacks before pt cut:
      if( (useJetID && pfJetIDOK) || (!useJetID) ) {

        hp_Rch_vs_pt->Fill( ptJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
        hp_Rnh_vs_pt->Fill( ptJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
        hp_Rgamma_vs_pt->Fill( ptJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
        hp_Rmu_vs_pt->Fill( ptJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
        hp_Re_vs_pt->Fill( ptJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
        hp_Rhfhad_vs_pt->Fill( ptJetReco[iJetReco], eHFHadronsReco[iJetReco], eventWeight );
        hp_Rhfem_vs_pt->Fill( ptJetReco[iJetReco], eHFEMReco[iJetReco], eventWeight );
      
        hp_Ech_vs_pt->Fill( ptJetReco[iJetReco], eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Enh_vs_pt->Fill( ptJetReco[iJetReco], eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Egamma_vs_pt->Fill( ptJetReco[iJetReco], ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Emu_vs_pt->Fill( ptJetReco[iJetReco], eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Ee_vs_pt->Fill( ptJetReco[iJetReco], eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Ehfhad_vs_pt->Fill( ptJetReco[iJetReco], eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Ehfem_vs_pt->Fill( ptJetReco[iJetReco], eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        hp_Nch_vs_pt->Fill( ptJetReco[iJetReco], nTracksReco[iJetReco] , eventWeight);
        hp_Nnh_vs_pt->Fill( ptJetReco[iJetReco], nNeutralHadronsReco[iJetReco], eventWeight );
        hp_Ngamma_vs_pt->Fill( ptJetReco[iJetReco], nPhotonsReco[iJetReco] , eventWeight);
        hp_Nmu_vs_pt->Fill( ptJetReco[iJetReco], nMuonsReco[iJetReco], eventWeight );
        hp_Ne_vs_pt->Fill( ptJetReco[iJetReco], nElectronsReco[iJetReco], eventWeight );
        hp_Nhfhad_vs_pt->Fill( ptJetReco[iJetReco], nHFHadronsReco[iJetReco], eventWeight );
        hp_Nhfem_vs_pt->Fill( ptJetReco[iJetReco], nHFEMReco[iJetReco], eventWeight );

        hp_Rch_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
        hp_Rnh_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
        hp_Rgamma_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
        hp_Rmu_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
        hp_Re_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
        hp_Rhfhad_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], eHFHadronsReco[iJetReco], eventWeight );
        hp_Rhfem_vs_ptCorr->Fill( ptCorrJetReco[iJetReco], eHFEMReco[iJetReco] , eventWeight);
      
        if( fabs(etaJetReco[iJetReco])<1.4 ) {
      
          hp_Rch_vs_pt_barrel->Fill( ptJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
          hp_Rnh_vs_pt_barrel->Fill( ptJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
          hp_Rgamma_vs_pt_barrel->Fill( ptJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
          hp_Rmu_vs_pt_barrel->Fill( ptJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
          hp_Re_vs_pt_barrel->Fill( ptJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
          hp_Rhfhad_vs_pt_barrel->Fill( ptJetReco[iJetReco], eHFHadronsReco[iJetReco] , eventWeight);
          hp_Rhfem_vs_pt_barrel->Fill( ptJetReco[iJetReco], eHFEMReco[iJetReco], eventWeight );
        
          hp_Ech_vs_pt_barrel->Fill( ptJetReco[iJetReco], eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
          hp_Enh_vs_pt_barrel->Fill( ptJetReco[iJetReco], eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
          hp_Egamma_vs_pt_barrel->Fill( ptJetReco[iJetReco], ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
          hp_Emu_vs_pt_barrel->Fill( ptJetReco[iJetReco], eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
          hp_Ee_vs_pt_barrel->Fill( ptJetReco[iJetReco], eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
          hp_Ehfhad_vs_pt_barrel->Fill( ptJetReco[iJetReco], eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
          hp_Ehfem_vs_pt_barrel->Fill( ptJetReco[iJetReco], eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      
          hp_Nch_vs_pt_barrel->Fill( ptJetReco[iJetReco], nTracksReco[iJetReco], eventWeight );
          hp_Nnh_vs_pt_barrel->Fill( ptJetReco[iJetReco], nNeutralHadronsReco[iJetReco] , eventWeight);
          hp_Ngamma_vs_pt_barrel->Fill( ptJetReco[iJetReco], nPhotonsReco[iJetReco] , eventWeight);
          hp_Nmu_vs_pt_barrel->Fill( ptJetReco[iJetReco], nMuonsReco[iJetReco] , eventWeight);
          hp_Ne_vs_pt_barrel->Fill( ptJetReco[iJetReco], nElectronsReco[iJetReco], eventWeight );
          hp_Nhfhad_vs_pt_barrel->Fill( ptJetReco[iJetReco], nHFHadronsReco[iJetReco] , eventWeight);
          hp_Nhfem_vs_pt_barrel->Fill( ptJetReco[iJetReco], nHFEMReco[iJetReco] , eventWeight);
      
          hp_Rch_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
          hp_Rnh_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], eNeutralHadronsReco[iJetReco] , eventWeight);
          hp_Rgamma_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
          hp_Rmu_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
          hp_Re_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], eElectronsReco[iJetReco] , eventWeight);
          hp_Rhfhad_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], eHFHadronsReco[iJetReco] , eventWeight);
          hp_Rhfem_vs_ptCorr_barrel->Fill( ptCorrJetReco[iJetReco], eHFEMReco[iJetReco] , eventWeight);
      
      //}   else if( fabs(etaJetReco[iJetReco])<2.5 ) {
      
      //  hp_Rch_vs_pt_endcap->Fill( ptJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
      //  hp_Rnh_vs_pt_endcap->Fill( ptJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
      //  hp_Rgamma_vs_pt_endcap->Fill( ptJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
      //  hp_Rmu_vs_pt_endcap->Fill( ptJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
      //  hp_Re_vs_pt_endcap->Fill( ptJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
      //  hp_Rhfhad_vs_pt_endcap->Fill( ptJetReco[iJetReco], eHFHadronsReco[iJetReco], eventWeight );
      //  hp_Rhfem_vs_pt_endcap->Fill( ptJetReco[iJetReco], eHFEMReco[iJetReco], eventWeight );
      //
      //  hp_Ech_vs_pt_endcap->Fill( ptJetReco[iJetReco], eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      //  hp_Enh_vs_pt_endcap->Fill( ptJetReco[iJetReco], eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      //  hp_Egamma_vs_pt_endcap->Fill( ptJetReco[iJetReco], ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      //  hp_Emu_vs_pt_endcap->Fill( ptJetReco[iJetReco], eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      //  hp_Ee_vs_pt_endcap->Fill( ptJetReco[iJetReco], eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      //  hp_Ehfhad_vs_pt_endcap->Fill( ptJetReco[iJetReco], eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      //  hp_Ehfem_vs_pt_endcap->Fill( ptJetReco[iJetReco], eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );
      
      //  hp_Nch_vs_pt_endcap->Fill( ptJetReco[iJetReco], nTracksReco[iJetReco], eventWeight );
      //  hp_Nnh_vs_pt_endcap->Fill( ptJetReco[iJetReco], nNeutralHadronsReco[iJetReco], eventWeight );
      //  hp_Ngamma_vs_pt_endcap->Fill( ptJetReco[iJetReco], nPhotonsReco[iJetReco], eventWeight );
      //  hp_Nmu_vs_pt_endcap->Fill( ptJetReco[iJetReco], nMuonsReco[iJetReco], eventWeight );
      //  hp_Ne_vs_pt_endcap->Fill( ptJetReco[iJetReco], nElectronsReco[iJetReco], eventWeight );
      //  hp_Nhfhad_vs_pt_endcap->Fill( ptJetReco[iJetReco], nHFHadronsReco[iJetReco], eventWeight );
      //  hp_Nhfem_vs_pt_endcap->Fill( ptJetReco[iJetReco], nHFEMReco[iJetReco], eventWeight );
      
        }

        if( fabs(etaJetGen[iJetReco])<1.4 ) {

          hp_RchGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], eTracksGen[iJetReco], eventWeight );
          hp_RnhGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], eNeutralHadronsGen[iJetReco], eventWeight );
          hp_RgammaGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], ePhotonsGen[iJetReco], eventWeight );
          hp_RmuGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], eMuonsGen[iJetReco], eventWeight );
          hp_ReGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], eElectronsGen[iJetReco], eventWeight );
          hp_RhfhadGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], eHFHadronsGen[iJetReco] , eventWeight);
          hp_RhfemGen_vs_pt_barrel->Fill( ptJetGen[iJetReco], eHFEMGen[iJetReco], eventWeight );

          if( ptJetGen[iJetReco]>pt_thresh ) {
            hp_RchGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eTracksGen[iJetReco], eventWeight );
            hp_RnhGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eNeutralHadronsGen[iJetReco], eventWeight );
            hp_RgammaGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], ePhotonsGen[iJetReco], eventWeight );
            hp_RmuGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eMuonsGen[iJetReco], eventWeight );
            hp_ReGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eElectronsGen[iJetReco], eventWeight );
            hp_RhfhadGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eHFHadronsGen[iJetReco], eventWeight );
            hp_RhfemGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eHFEMGen[iJetReco] , eventWeight);
          }

        }
    

        if( ptJetGen[iJetReco] > pt_thresh ) {

          hp_RchGen_vs_eta_stack->Fill( etaJetGen[iJetReco], eTracksGen[iJetReco], eventWeight );
          hp_RnhGen_vs_eta_stack->Fill( etaJetGen[iJetReco], eNeutralHadronsGen[iJetReco], eventWeight );
          hp_RgammaGen_vs_eta_stack->Fill( etaJetGen[iJetReco], ePhotonsGen[iJetReco] , eventWeight);
          hp_RmuGen_vs_eta_stack->Fill( etaJetGen[iJetReco], eMuonsGen[iJetReco], eventWeight );
          hp_ReGen_vs_eta_stack->Fill( etaJetGen[iJetReco], eElectronsGen[iJetReco], eventWeight );
          hp_RhfhadGen_vs_eta_stack->Fill( etaJetGen[iJetReco], eHFHadronsGen[iJetReco], eventWeight );
          hp_RhfemGen_vs_eta_stack->Fill( etaJetGen[iJetReco], eHFEMGen[iJetReco], eventWeight );
          
        }
        

      } //if use pfjetid


      //pt cut:
      if( raw_corr=="raw" ) {
        if( ptJetReco[iJetReco]<pt_thresh ) continue;
      } else {
        if( ptCorrJetReco[iJetReco]<pt_thresh ) continue;
      }


      //jet id:
      if( (useJetID) && (!pfJetIDOK) ) continue;
      

      hp_Rch_vs_eta_stack->Fill( etaJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
      hp_Rnh_vs_eta_stack->Fill( etaJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
      hp_Rgamma_vs_eta_stack->Fill( etaJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
      hp_Rmu_vs_eta_stack->Fill( etaJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
      hp_Re_vs_eta_stack->Fill( etaJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
      hp_Rhfhad_vs_eta_stack->Fill( etaJetReco[iJetReco], eHFHadronsReco[iJetReco], eventWeight );
      hp_Rhfem_vs_eta_stack->Fill( etaJetReco[iJetReco], eHFEMReco[iJetReco], eventWeight );
    
      if( fabs( etaJetReco[iJetReco] ) < 1.4 ) {
      
        h1_NchJet_barrel->Fill( nTracksReco[iJetReco], eventWeight );
        h1_NnhJet_barrel->Fill( nNeutralHadronsReco[iJetReco], eventWeight );
        h1_NgammaJet_barrel->Fill( nPhotonsReco[iJetReco], eventWeight );
        h1_NeJet_barrel->Fill( nElectronsReco[iJetReco], eventWeight );
        h1_NmuJet_barrel->Fill( nMuonsReco[iJetReco], eventWeight );
        h1_NhfhadJet_barrel->Fill( nHFHadronsReco[iJetReco], eventWeight );
        h1_NhfemJet_barrel->Fill( nHFEMReco[iJetReco], eventWeight );
      
        h1_RchJet_barrel->Fill( eTracksReco[iJetReco], eventWeight );
        h1_RnhJet_barrel->Fill( eNeutralHadronsReco[iJetReco], eventWeight );
        h1_RgammaJet_barrel->Fill( ePhotonsReco[iJetReco], eventWeight );
        h1_ReJet_barrel->Fill( eElectronsReco[iJetReco], eventWeight );
        h1_RmuJet_barrel->Fill( eMuonsReco[iJetReco], eventWeight );
        h1_RhfhadJet_barrel->Fill( eHFHadronsReco[iJetReco], eventWeight );
        h1_RhfemJet_barrel->Fill( eHFEMReco[iJetReco], eventWeight );
      
        h1_EchJet_barrel->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EnhJet_barrel->Fill( eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EgammaJet_barrel->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EeJet_barrel->Fill( eElectronsReco[iJetReco]*eJetReco[iJetReco] , eventWeight);
        h1_EmuJet_barrel->Fill( eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfhadJet_barrel->Fill( eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfemJet_barrel->Fill( eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        h1_EgammanhJet_barrel->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco])*eJetReco[iJetReco], eventWeight );
        h1_RgammanhJet_barrel->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco]), eventWeight );
        h1_EgammanhJet_barrel->Fill( (nPhotonsReco[iJetReco]+nNeutralHadronsReco[iJetReco]), eventWeight );

        hp_Rch_vs_phi_barrel->Fill( phiJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
        hp_Rnh_vs_phi_barrel->Fill( phiJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
        hp_Rgamma_vs_phi_barrel->Fill( phiJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
        hp_Rmu_vs_phi_barrel->Fill( phiJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
        hp_Re_vs_phi_barrel->Fill( phiJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
        hp_Rhfhad_vs_phi_barrel->Fill( phiJetReco[iJetReco], eHFHadronsReco[iJetReco], eventWeight );
        hp_Rhfem_vs_phi_barrel->Fill( phiJetReco[iJetReco], eHFEMReco[iJetReco], eventWeight );
      
      } else if( fabs(etaJetReco[iJetReco])<2.5 ) {

        h1_NchJet_eta1425->Fill( nTracksReco[iJetReco] , eventWeight);
        h1_NnhJet_eta1425->Fill( nNeutralHadronsReco[iJetReco], eventWeight );
        h1_NgammaJet_eta1425->Fill( nPhotonsReco[iJetReco], eventWeight );
        h1_NeJet_eta1425->Fill( nElectronsReco[iJetReco], eventWeight );
        h1_NmuJet_eta1425->Fill( nMuonsReco[iJetReco], eventWeight );
        h1_NhfhadJet_eta1425->Fill( nHFHadronsReco[iJetReco], eventWeight );
        h1_NhfemJet_eta1425->Fill( nHFEMReco[iJetReco], eventWeight );
      
        h1_RchJet_eta1425->Fill( eTracksReco[iJetReco], eventWeight );
        h1_RnhJet_eta1425->Fill( eNeutralHadronsReco[iJetReco], eventWeight );
        h1_RgammaJet_eta1425->Fill( ePhotonsReco[iJetReco], eventWeight );
        h1_ReJet_eta1425->Fill( eElectronsReco[iJetReco], eventWeight );
        h1_RmuJet_eta1425->Fill( eMuonsReco[iJetReco], eventWeight );
        h1_RhfhadJet_eta1425->Fill( eHFHadronsReco[iJetReco], eventWeight );
        h1_RhfemJet_eta1425->Fill( eHFEMReco[iJetReco], eventWeight );
      
        h1_EchJet_eta1425->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EnhJet_eta1425->Fill( eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EgammaJet_eta1425->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EeJet_eta1425->Fill( eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EmuJet_eta1425->Fill( eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfhadJet_eta1425->Fill( eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfemJet_eta1425->Fill( eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        h1_EgammanhJet_eta1425->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco])*eJetReco[iJetReco], eventWeight );
        h1_RgammanhJet_eta1425->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco]), eventWeight );
        h1_EgammanhJet_eta1425->Fill( (nPhotonsReco[iJetReco]+nNeutralHadronsReco[iJetReco]), eventWeight );

      }

      if( fabs(etaJetReco[iJetReco]) < 3. && fabs(etaJetReco[iJetReco]) > 1.6 ) {

        h1_NchJet_eta163->Fill( nTracksReco[iJetReco], eventWeight );
        h1_NnhJet_eta163->Fill( nNeutralHadronsReco[iJetReco], eventWeight );
        h1_NgammaJet_eta163->Fill( nPhotonsReco[iJetReco], eventWeight );
        h1_NeJet_eta163->Fill( nElectronsReco[iJetReco], eventWeight );
        h1_NmuJet_eta163->Fill( nMuonsReco[iJetReco], eventWeight );
        h1_NhfhadJet_eta163->Fill( nHFHadronsReco[iJetReco], eventWeight );
        h1_NhfemJet_eta163->Fill( nHFEMReco[iJetReco], eventWeight );
      
        h1_RchJet_eta163->Fill( eTracksReco[iJetReco], eventWeight );
        h1_RnhJet_eta163->Fill( eNeutralHadronsReco[iJetReco], eventWeight );
        h1_RgammaJet_eta163->Fill( ePhotonsReco[iJetReco], eventWeight );
        h1_ReJet_eta163->Fill( eElectronsReco[iJetReco], eventWeight );
        h1_RmuJet_eta163->Fill( eMuonsReco[iJetReco], eventWeight );
        h1_RhfhadJet_eta163->Fill( eHFHadronsReco[iJetReco], eventWeight );
        h1_RhfemJet_eta163->Fill( eHFEMReco[iJetReco], eventWeight );
      
        h1_EchJet_eta163->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EnhJet_eta163->Fill( eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EgammaJet_eta163->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EeJet_eta163->Fill( eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EmuJet_eta163->Fill( eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfhadJet_eta163->Fill( eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfemJet_eta163->Fill( eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        h1_EgammanhJet_eta163->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco])*eJetReco[iJetReco] , eventWeight);
        h1_RgammanhJet_eta163->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco]), eventWeight );
        h1_EgammanhJet_eta163->Fill( (nPhotonsReco[iJetReco]+nNeutralHadronsReco[iJetReco]), eventWeight );

      } //if endcap

      if( fabs( etaJetReco[iJetReco] ) < 2. ) {

        h1_NchJet_eta02->Fill( nTracksReco[iJetReco], eventWeight );
        h1_NnhJet_eta02->Fill( nNeutralHadronsReco[iJetReco], eventWeight );
        h1_NgammaJet_eta02->Fill( nPhotonsReco[iJetReco], eventWeight );
        h1_NeJet_eta02->Fill( nElectronsReco[iJetReco], eventWeight );
        h1_NmuJet_eta02->Fill( nMuonsReco[iJetReco], eventWeight );
        h1_NhfhadJet_eta02->Fill( nHFHadronsReco[iJetReco], eventWeight );
        h1_NhfemJet_eta02->Fill( nHFEMReco[iJetReco], eventWeight );
      
        h1_RchJet_eta02->Fill( eTracksReco[iJetReco], eventWeight );
        h1_RnhJet_eta02->Fill( eNeutralHadronsReco[iJetReco], eventWeight );
        h1_RgammaJet_eta02->Fill( ePhotonsReco[iJetReco], eventWeight );
        h1_ReJet_eta02->Fill( eElectronsReco[iJetReco], eventWeight );
        h1_RmuJet_eta02->Fill( eMuonsReco[iJetReco], eventWeight );
        h1_RhfhadJet_eta02->Fill( eHFHadronsReco[iJetReco], eventWeight );
        h1_RhfemJet_eta02->Fill( eHFEMReco[iJetReco], eventWeight );
      
        h1_EchJet_eta02->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EnhJet_eta02->Fill( eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EgammaJet_eta02->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EeJet_eta02->Fill( eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EmuJet_eta02->Fill( eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfhadJet_eta02->Fill( eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfemJet_eta02->Fill( eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

      } else if( fabs(etaJetReco[iJetReco]) < 3. ) {

        h1_NchJet_eta23->Fill( nTracksReco[iJetReco], eventWeight );
        h1_NnhJet_eta23->Fill( nNeutralHadronsReco[iJetReco], eventWeight );
        h1_NgammaJet_eta23->Fill( nPhotonsReco[iJetReco], eventWeight );
        h1_NeJet_eta23->Fill( nElectronsReco[iJetReco], eventWeight );
        h1_NmuJet_eta23->Fill( nMuonsReco[iJetReco], eventWeight );
        h1_NhfhadJet_eta23->Fill( nHFHadronsReco[iJetReco], eventWeight );
        h1_NhfemJet_eta23->Fill( nHFEMReco[iJetReco], eventWeight );
      
        h1_RchJet_eta23->Fill( eTracksReco[iJetReco], eventWeight );
        h1_RnhJet_eta23->Fill( eNeutralHadronsReco[iJetReco], eventWeight );
        h1_RgammaJet_eta23->Fill( ePhotonsReco[iJetReco], eventWeight );
        h1_ReJet_eta23->Fill( eElectronsReco[iJetReco], eventWeight );
        h1_RmuJet_eta23->Fill( eMuonsReco[iJetReco], eventWeight );
        h1_RhfhadJet_eta23->Fill( eHFHadronsReco[iJetReco], eventWeight );
        h1_RhfemJet_eta23->Fill( eHFEMReco[iJetReco], eventWeight );
      
        h1_EchJet_eta23->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EnhJet_eta23->Fill( eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EgammaJet_eta23->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EeJet_eta23->Fill( eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EmuJet_eta23->Fill( eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfhadJet_eta23->Fill( eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfemJet_eta23->Fill( eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

      }

//    if( fabs(etaJetGen[iJetReco]) < 1.4 ) {

//      hp_RchGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eTracksGen[iJetReco] );
//      hp_RnhGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eNeutralHadronsGen[iJetReco] );
//      hp_RgammaGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], ePhotonsGen[iJetReco] );
//      hp_RmuGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eMuonsGen[iJetReco] );
//      hp_ReGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eElectronsGen[iJetReco] );
//      hp_RhfhadGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eHFHadronsGen[iJetReco] );
//      hp_RhfemGen_vs_phi_barrel->Fill( phiJetGen[iJetReco], eHFEMGen[iJetReco] );
//      
//    }


      if( fabs(etaJetReco[iJetReco])<etaMax ) {

        h1_ptJet->Fill( ptJetReco[iJetReco], eventWeight );
        h1_ptCorrJet->Fill( ptCorrJetReco[iJetReco], eventWeight );

        hp_pt_vs_eta->Fill( etaJetReco[iJetReco], ptJetReco[iJetReco], eventWeight );
        hp_ptCorr_vs_eta->Fill( etaJetReco[iJetReco], ptCorrJetReco[iJetReco], eventWeight );
        

        Float_t i_thetaJetReco = 2.*atan(exp(-etaJetReco[iJetReco]));
        Float_t i_pJetReco = ptJetReco[iJetReco]/sin(i_thetaJetReco);
        Float_t i_pxJetReco = i_pJetReco*sin(i_thetaJetReco)*cos(phiJetReco[iJetReco]);
        Float_t i_pyJetReco = i_pJetReco*sin(i_thetaJetReco)*sin(phiJetReco[iJetReco]);
        Float_t i_pzJetReco = i_pJetReco*cos(i_thetaJetReco);
        Float_t i_massJet = sqrt( eJetReco[iJetReco]*eJetReco[iJetReco] - i_pJetReco*i_pJetReco );
        h1_massJet->Fill( i_massJet , eventWeight);
        h1_MoEJet->Fill( i_massJet/eJetReco[iJetReco], eventWeight );
        h1_ptOverMJet->Fill( ptJetReco[iJetReco]/i_massJet, eventWeight );
        h1_pOverMJet->Fill( i_pJetReco/i_massJet, eventWeight );

      
        h1_etaJet->Fill( etaJetReco[iJetReco], eventWeight );
        h1_phiJet->Fill( phiJetReco[iJetReco], eventWeight );
        h1_nCandJet->Fill( N, eventWeight );
      
        h1_RchJet->Fill( eTracksReco[iJetReco], eventWeight );
        h1_RnhJet->Fill( eNeutralHadronsReco[iJetReco], eventWeight );
        h1_RgammaJet->Fill( ePhotonsReco[iJetReco], eventWeight );
        h1_ReJet->Fill( eElectronsReco[iJetReco], eventWeight );
        h1_RmuJet->Fill( eMuonsReco[iJetReco], eventWeight );
        h1_RhfhadJet->Fill( eHFHadronsReco[iJetReco], eventWeight );
        h1_RhfemJet->Fill( eHFEMReco[iJetReco], eventWeight );
      
        h1_EchJet->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EnhJet->Fill( eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EgammaJet->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EeJet->Fill( eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EmuJet->Fill( eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfhadJet->Fill( eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        h1_EhfemJet->Fill( eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        h1_PTchJet->Fill( ptTracksReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        h1_PTnhJet->Fill( ptNeutralHadronsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        h1_PTgammaJet->Fill( ptPhotonsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        h1_PTeJet->Fill( ptElectronsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        h1_PTmuJet->Fill( ptMuonsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        h1_PThfhadJet->Fill( ptHFHadronsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        h1_PThfemJet->Fill( ptHFEMReco[iJetReco]*ptJetReco[iJetReco], eventWeight );

        h1_EgammanhJet->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco])*eJetReco[iJetReco], eventWeight );
        h1_RgammanhJet->Fill( (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco]), eventWeight );
        h1_EgammanhJet->Fill( (nPhotonsReco[iJetReco]+nNeutralHadronsReco[iJetReco]), eventWeight );


        if( nPhotonsReco[iJetReco]>0 )
          h1_EphotAveJet->Fill( ePhotonsReco[iJetReco]*eJetReco[iJetReco]/(Float_t)nPhotonsReco[iJetReco], eventWeight );
      
        h1_NchJet->Fill( nTracksReco[iJetReco], eventWeight );
        h1_NnhJet->Fill( nNeutralHadronsReco[iJetReco], eventWeight );
        h1_NgammaJet->Fill( nPhotonsReco[iJetReco], eventWeight );
        h1_NeJet->Fill( nElectronsReco[iJetReco], eventWeight );
        h1_NmuJet->Fill( nMuonsReco[iJetReco], eventWeight );
        h1_NhfhadJet->Fill( nHFHadronsReco[iJetReco], eventWeight );
        h1_NhfemJet->Fill( nHFEMReco[iJetReco], eventWeight );
      
        hp_Rch_vs_eta->Fill( etaJetReco[iJetReco], eTracksReco[iJetReco], eventWeight );
        hp_Rnh_vs_eta->Fill( etaJetReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
        hp_Rgamma_vs_eta->Fill( etaJetReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
        hp_Rgammanh_vs_eta->Fill( etaJetReco[iJetReco], (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco]), eventWeight);
        hp_Rmu_vs_eta->Fill( etaJetReco[iJetReco], eMuonsReco[iJetReco], eventWeight );
        hp_Re_vs_eta->Fill( etaJetReco[iJetReco], eElectronsReco[iJetReco], eventWeight );
        hp_Rhfhad_vs_eta->Fill( etaJetReco[iJetReco], eHFHadronsReco[iJetReco], eventWeight );
        hp_Rhfem_vs_eta->Fill( etaJetReco[iJetReco], eHFEMReco[iJetReco], eventWeight );
      
        hp_Ech_vs_eta->Fill( etaJetReco[iJetReco], eTracksReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Enh_vs_eta->Fill( etaJetReco[iJetReco], eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Egamma_vs_eta->Fill( etaJetReco[iJetReco], ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Egammanh_vs_eta->Fill( etaJetReco[iJetReco], (ePhotonsReco[iJetReco]+eNeutralHadronsReco[iJetReco])*eJetReco[iJetReco], eventWeight );
        hp_Emu_vs_eta->Fill( etaJetReco[iJetReco], eMuonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Ee_vs_eta->Fill( etaJetReco[iJetReco], eElectronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Ehfhad_vs_eta->Fill( etaJetReco[iJetReco], eHFHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Ehfem_vs_eta->Fill( etaJetReco[iJetReco], eHFEMReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        hp_PTch_vs_eta->Fill( etaJetReco[iJetReco], ptTracksReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        hp_PTnh_vs_eta->Fill( etaJetReco[iJetReco], ptNeutralHadronsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        hp_PTgamma_vs_eta->Fill( etaJetReco[iJetReco], ptPhotonsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        hp_PTmu_vs_eta->Fill( etaJetReco[iJetReco], ptMuonsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        hp_PTe_vs_eta->Fill( etaJetReco[iJetReco], ptElectronsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        hp_PThfhad_vs_eta->Fill( etaJetReco[iJetReco], ptHFHadronsReco[iJetReco]*ptJetReco[iJetReco], eventWeight );
        hp_PThfem_vs_eta->Fill( etaJetReco[iJetReco], ptHFEMReco[iJetReco]*ptJetReco[iJetReco], eventWeight );

        hp_Nch_vs_eta->Fill( etaJetReco[iJetReco], nTracksReco[iJetReco], eventWeight );
        hp_Nnh_vs_eta->Fill( etaJetReco[iJetReco], nNeutralHadronsReco[iJetReco], eventWeight );
        hp_Ngamma_vs_eta->Fill( etaJetReco[iJetReco], nPhotonsReco[iJetReco], eventWeight );
        hp_Ngammanh_vs_eta->Fill( etaJetReco[iJetReco], (nPhotonsReco[iJetReco]+nNeutralHadronsReco[iJetReco]), eventWeight);
        hp_Nmu_vs_eta->Fill( etaJetReco[iJetReco], nMuonsReco[iJetReco], eventWeight );
        hp_Ne_vs_eta->Fill( etaJetReco[iJetReco], nElectronsReco[iJetReco], eventWeight );
        hp_Nhfhad_vs_eta->Fill( etaJetReco[iJetReco], nHFHadronsReco[iJetReco], eventWeight );
        hp_Nhfem_vs_eta->Fill( etaJetReco[iJetReco], nHFEMReco[iJetReco], eventWeight );

     

        //correlations between components:
        hp_Egamma_vs_Ech->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], ePhotonsReco[iJetReco]*eJetReco[iJetReco], eventWeight );
        hp_Enh_vs_Ech->Fill( eTracksReco[iJetReco]*eJetReco[iJetReco], eNeutralHadronsReco[iJetReco]*eJetReco[iJetReco], eventWeight );

        hp_Rgamma_vs_Rch->Fill( eTracksReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
        hp_Rnh_vs_Rch->Fill( eTracksReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
        hp_Rnh_vs_Rgamma->Fill( ePhotonsReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );

      //h2_Rgamma_vs_Rch->Fill( eTracksReco[iJetReco], ePhotonsReco[iJetReco], eventWeight );
      //h2_Rnh_vs_Rch->Fill( eTracksReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );
      //h2_Rnh_vs_Rgamma->Fill( ePhotonsReco[iJetReco], eNeutralHadronsReco[iJetReco], eventWeight );

        hp_Ngamma_vs_Nch->Fill( nTracksReco[iJetReco], nPhotonsReco[iJetReco], eventWeight );
        hp_Nnh_vs_Nch->Fill( nTracksReco[iJetReco], nNeutralHadronsReco[iJetReco] , eventWeight);

      } //if eta < etamax

    } //for reco jets

  } //for entries


  std::string outfileName_str;

  if(dataset!="") outfileName_str = "Inclusive_"+dataset;
  else outfileName_str = "Inclusive";

  outfileName_str += "_" +algoType;

  char outfileName[150];
  if( useJetID )
    sprintf( outfileName, "%s_%dGeV%s_eta%d.root", outfileName_str.c_str(), (Int_t)floor(pt_thresh), raw_corr.c_str(), (int)(10.*etaMax));
  else
    sprintf( outfileName, "%s_%dGeV%s_eta%d_noID.root", outfileName_str.c_str(), (Int_t)floor(pt_thresh), raw_corr.c_str(), (int)(10.*etaMax));

  TFile* outFile = new TFile(outfileName, "RECREATE");
  outFile->cd();

  h1_totalLumi->Write();

  h1_nJets->Write();

  h1_massJet->Write();
  h1_MoEJet->Write();
  h1_ptOverMJet->Write();
  h1_pOverMJet->Write();
  h1_nCandJet->Write();

  h1_etaJet->Write();
  h1_phiJet->Write();
  h1_ptJet->Write();
  h1_ptCorrJet->Write();

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

  h1_PTchJet->Write();
  h1_PTnhJet->Write();
  h1_PTgammaJet->Write();
  h1_PTeJet->Write();
  h1_PTmuJet->Write();
  h1_PThfhadJet->Write();
  h1_PThfemJet->Write();

  h1_EphotAveJet->Write();

  h1_EgammanhJet->Write();
  h1_RgammanhJet->Write();
  h1_NgammanhJet->Write();

  h1_NchJet->Write();
  h1_NnhJet->Write();
  h1_NgammaJet->Write();
  h1_NeJet->Write();
  h1_NmuJet->Write();
  h1_NhfhadJet->Write();
  h1_NhfemJet->Write();

  h1_RchJet_eta02->Write();
  h1_RnhJet_eta02->Write();
  h1_RgammaJet_eta02->Write();
  h1_ReJet_eta02->Write();
  h1_RmuJet_eta02->Write();
  h1_RhfhadJet_eta02->Write();
  h1_RhfemJet_eta02->Write();

  h1_EchJet_eta02->Write();
  h1_EnhJet_eta02->Write();
  h1_EgammaJet_eta02->Write();
  h1_EeJet_eta02->Write();
  h1_EmuJet_eta02->Write();
  h1_EhfhadJet_eta02->Write();
  h1_EhfemJet_eta02->Write();

  h1_NchJet_eta02->Write();
  h1_NnhJet_eta02->Write();
  h1_NgammaJet_eta02->Write();
  h1_NeJet_eta02->Write();
  h1_NmuJet_eta02->Write();
  h1_NhfhadJet_eta02->Write();
  h1_NhfemJet_eta02->Write();

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

  h1_EgammanhJet_barrel->Write();
  h1_RgammanhJet_barrel->Write();
  h1_NgammanhJet_barrel->Write();

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

  h1_EgammanhJet_eta1425->Write();
  h1_RgammanhJet_eta1425->Write();
  h1_NgammanhJet_eta1425->Write();

  h1_RchJet_eta163->Write();
  h1_RnhJet_eta163->Write();
  h1_RgammaJet_eta163->Write();
  h1_ReJet_eta163->Write();
  h1_RmuJet_eta163->Write();
  h1_RhfhadJet_eta163->Write();
  h1_RhfemJet_eta163->Write();

  h1_EchJet_eta163->Write();
  h1_EnhJet_eta163->Write();
  h1_EgammaJet_eta163->Write();
  h1_EeJet_eta163->Write();
  h1_EmuJet_eta163->Write();
  h1_EhfhadJet_eta163->Write();
  h1_EhfemJet_eta163->Write();

  h1_NchJet_eta163->Write();
  h1_NnhJet_eta163->Write();
  h1_NgammaJet_eta163->Write();
  h1_NeJet_eta163->Write();
  h1_NmuJet_eta163->Write();
  h1_NhfhadJet_eta163->Write();
  h1_NhfemJet_eta163->Write();

  h1_EgammanhJet_eta163->Write();
  h1_RgammanhJet_eta163->Write();
  h1_NgammanhJet_eta163->Write();

  h1_RchJet_eta23->Write();
  h1_RnhJet_eta23->Write();
  h1_RgammaJet_eta23->Write();
  h1_ReJet_eta23->Write();
  h1_RmuJet_eta23->Write();
  h1_RhfhadJet_eta23->Write();
  h1_RhfemJet_eta23->Write();

  h1_EchJet_eta23->Write();
  h1_EnhJet_eta23->Write();
  h1_EgammaJet_eta23->Write();
  h1_EeJet_eta23->Write();
  h1_EmuJet_eta23->Write();
  h1_EhfhadJet_eta23->Write();
  h1_EhfemJet_eta23->Write();

  h1_NchJet_eta23->Write();
  h1_NnhJet_eta23->Write();
  h1_NgammaJet_eta23->Write();
  h1_NeJet_eta23->Write();
  h1_NmuJet_eta23->Write();
  h1_NhfhadJet_eta23->Write();
  h1_NhfemJet_eta23->Write();

  hp_pt_vs_eta->Write();
  hp_ptCorr_vs_eta->Write();

  hp_Rch_vs_eta->Write();
  hp_Rnh_vs_eta->Write();
  hp_Rgamma_vs_eta->Write();
  hp_Rgammanh_vs_eta->Write();
  hp_Rmu_vs_eta->Write();
  hp_Re_vs_eta->Write();
  hp_Rhfhad_vs_eta->Write();
  hp_Rhfem_vs_eta->Write();

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

  hp_Ech_vs_eta->Write();
  hp_Enh_vs_eta->Write();
  hp_Egamma_vs_eta->Write();
  hp_Egammanh_vs_eta->Write();
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
  hp_Ngammanh_vs_eta->Write();
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

  hp_RchGen_vs_pt_barrel->Write();
  hp_RnhGen_vs_pt_barrel->Write();
  hp_RgammaGen_vs_pt_barrel->Write();
  hp_RmuGen_vs_pt_barrel->Write();
  hp_ReGen_vs_pt_barrel->Write();
  hp_RhfhadGen_vs_pt_barrel->Write();
  hp_RhfemGen_vs_pt_barrel->Write();

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

//hp_Rch_vs_ptCorr->Write();
//hp_Rnh_vs_ptCorr->Write();
//hp_Rgamma_vs_ptCorr->Write();
//hp_Rmu_vs_ptCorr->Write();
//hp_Re_vs_ptCorr->Write();
//hp_Rhfhad_vs_ptCorr->Write();
//hp_Rhfem_vs_ptCorr->Write();

//hp_Ech_vs_ptCorr->Write();
//hp_Enh_vs_ptCorr->Write();
//hp_Egamma_vs_ptCorr->Write();
//hp_Emu_vs_ptCorr->Write();
//hp_Ee_vs_ptCorr->Write();
//hp_Ehfhad_vs_ptCorr->Write();
//hp_Ehfem_vs_ptCorr->Write();

//hp_Nch_vs_ptCorr->Write();
//hp_Nnh_vs_ptCorr->Write();
//hp_Ngamma_vs_ptCorr->Write();
//hp_Nmu_vs_ptCorr->Write();
//hp_Ne_vs_ptCorr->Write();
//hp_Nhfhad_vs_ptCorr->Write();
//hp_Nhfem_vs_ptCorr->Write();

  hp_Egamma_vs_Ech->Write();
  hp_Enh_vs_Ech->Write();

  hp_Rgamma_vs_Rch->Write();
  hp_Rnh_vs_Rch->Write();
  hp_Rnh_vs_Rgamma->Write();

//h2_Rgamma_vs_Rch->Write();
//h2_Rnh_vs_Rch->Write();
//h2_Rnh_vs_Rgamma->Write();

  hp_Ngamma_vs_Nch->Write();
  hp_Nnh_vs_Nch->Write();

  outFile->Close();


}
  

bool PFJetID( Float_t Rch, Float_t Rnh, Float_t Rgamma, Int_t N, Int_t Nch, Float_t eta) {

  bool pass = false;

  if( Rch>0. ) pass=true;
  else if( Rnh<1. &&  Rgamma<1. ) pass = true;

//if( fabs(eta)<2.4 ) {
//  if( Rch>0. ) pass = true;
//} else {
//  if( Rnh<1. && Rgamma<1. ) pass=true;
//}


//  if( Rnh<1. && Rch>0. && Rgamma<1. ) pass = true;

  //if( (N>1) && (Rch>0.15) ) pass=true;
  //if( Nch>1 ) pass=true;

  return pass;

}
