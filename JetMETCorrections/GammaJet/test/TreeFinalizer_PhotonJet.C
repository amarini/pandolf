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
std::string RECOTYPE_;


std::vector<TH1F*> getResponseHistos(const std::string& name);

void finalize(const std::string& dataset, std::string recoType, std::string jetAlgo, bool useGenJets=false, bool MCassoc=false) {

  RECOTYPE_ = recoType;

  std::string algoType = (recoType=="calo") ? jetAlgo : recoType+jetAlgo;
  if( recoType=="jpt"&&jetAlgo=="akt5" ) algoType="jptak5"; 

  TChain* tree = new TChain("jetTree");
  
  std::string infileName, treeName;
  TFile* infile=0;
  TH1F* h1_lumi=0;
  Double_t totalLumi=0.;

  std::string prefix = "PhotonJet_2ndLevelTree_";
  std::string suffix = "_"+algoType;
  if( useGenJets ) suffix = suffix + "_GENJETS";


  if( dataset=="DATA_7TeV_20100507") {

    infileName = prefix + "DATA_MinimumBias_Commissioning10-GOODCOLL-v8" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "Commissioning10-GOODCOLL-v8_runs132605_606" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "MinimumBias_Commissioning10-GOODCOLL-v9" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;


  } else if( dataset=="PhotonJet_Summer09") {

    infileName = prefix + "PhotonJet_Summer09_Pt15" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt30" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt80" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt170" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt300" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt470" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt800" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt1400" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="QCD_Summer09") {

    infileName = prefix + "QCD_Summer09_Pt15" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt30" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt80" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt170" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt300" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt470" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt800" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="PhotonJetQCD_Summer09") {

    infileName = prefix + "PhotonJet_Summer09_Pt15" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt30" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt80" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt170" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt300" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt470" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt800" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Summer09_Pt1400" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt15" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt30" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt80" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt170" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt300" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt470" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "QCD_Summer09_Pt800" + suffix +  ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="PhotonJet_Spring10") {

  //infileName = prefix + "PhotonJet_Spring10_Pt0to15" + suffix + ".root";
  //infile = TFile::Open(infileName.c_str(), "READ");
  //h1_lumi = (TH1F*)infile->Get("lumi");
  //if( h1_lumi!=0 ) {
  //  totalLumi += h1_lumi->GetBinContent(1);
  //} else {
  //  std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
  //}
  //treeName = infileName + "/jetTree";
  //tree->Add(treeName.c_str());
  //std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "PhotonJet_Spring10_Pt5to15" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "PhotonJet_Spring10_Pt15" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infileName = prefix + "PhotonJet_Spring10_Pt30" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "PhotonJet_Spring10_Pt80" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="QCD_Spring10") {

//  infileName = prefix + "QCD_Spring10_Pt0to15" + suffix + ".root";
//  treeName = infileName + "/jetTree";
//  tree->Add(treeName.c_str());
//  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "MinBias_357ReReco_v3_Pt0to15" + suffix + ".root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

//  infileName = prefix + "QCD_Spring10_Pt15" + suffix + ".root";
//  treeName = infileName + "/jetTree";
//  tree->Add(treeName.c_str());
//  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "QCD_Spring10_Pt15to20" + suffix + ".root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "QCD_Spring10_Pt20to30" + suffix + ".root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

//  infileName = prefix + "QCD_Spring10_Pt30" + suffix + ".root";
//  treeName = infileName + "/jetTree";
//  tree->Add(treeName.c_str());
//  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "QCD_Spring10_Pt30to50" + suffix + ".root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "QCD_Spring10_Pt50to80" + suffix + ".root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "QCD_Spring10_Pt80" + suffix + ".root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "QCD_Spring10_Pt170" + suffix + "_NOPTHATCUT.root";
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="DATA_CMSWEEK_MAY10") {

    infileName = prefix + "MinimumBias_Commissioning10-Apr20Skim_GOODCOLL-v1" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "MinimumBias_Commissioning10-GOODCOLL-v9" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  } else if( dataset=="DATA_EG") {

    infileName = prefix + "MinimumBias_Commissioning10_May6thPDSkim2_SD_EG" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "MinimumBias_Commissioning10_SD_EG-v9" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "EG_Run2010A-PromptReco-v1" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

    infileName = prefix + "EG_Run2010A-PromptReco-v2" + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;


  } else {

    infileName = prefix + dataset + suffix + ".root";
    infile = TFile::Open(infileName.c_str(), "READ");
    h1_lumi = (TH1F*)infile->Get("lumi");
    if( h1_lumi!=0 ) {
      totalLumi += h1_lumi->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
    }
    treeName = infileName + "/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;

  }

  std::cout << "-> Total integrated luminosity: " << totalLumi << " ub-1." << std::endl;
  TH1F* h1_totalLumi = new TH1F("totalLumi", "", 1, 0., 1.);
  h1_totalLumi->SetBinContent(1, totalLumi);

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Float_t ptMax = 100.;

  TH1D* h1_ptPhot = new TH1D("ptPhot", "", 10, ptPhot_binning[0], ptMax);
  h1_ptPhot->Sumw2();
  TH1D* h1_ptJetReco = new TH1D("ptJetReco", "", 10, 0., ptMax);
  h1_ptJetReco->Sumw2();
  TH1D* h1_pt2ndJetReco = new TH1D("pt2ndJetReco", "", 10, 5., ptMax);
  h1_pt2ndJetReco->Sumw2();
  TH1D* h1_etaPhot = new TH1D("etaPhot", "", 15, -1.3, 1.3);
  h1_etaPhot->Sumw2();
  TH1D* h1_phiPhot = new TH1D("phiPhot", "", 15, -3.1416, 3.1416);
  h1_phiPhot->Sumw2();
  TH1D* h1_deltaPhi_Nm1 = new TH1D("deltaPhi_Nm1", "", 15, 0., 3.1416);
  h1_deltaPhi_Nm1->Sumw2();
  TH1D* h1_ptSecondJetRel_Nm1 = new TH1D("ptSecondJetRel_Nm1", "", 15, 0., 1.5);
  h1_ptSecondJetRel_Nm1->Sumw2();

  TH1D* h1_ptPhot_loose = new TH1D("ptPhot_loose", "", 10, ptPhot_binning[0], ptMax);
  h1_ptPhot_loose->Sumw2();
  TH1D* h1_etaPhot_loose = new TH1D("etaPhot_loose", "", 15, -1.3, 1.3);
  h1_etaPhot_loose->Sumw2();
  TH1D* h1_phiPhot_loose = new TH1D("phiPhot_loose", "", 15, -3.1416, 3.1416);
  h1_phiPhot_loose->Sumw2();

  TH1D* h1_ptPhot_clusterOK = new TH1D("ptPhot_clusterOK", "", 10, ptPhot_binning[0], ptMax);
  h1_ptPhot_clusterOK->Sumw2();
  TH1D* h1_etaPhot_clusterOK = new TH1D("etaPhot_clusterOK", "", 15, -1.3, 1.3);
  h1_etaPhot_clusterOK->Sumw2();
  TH1D* h1_phiPhot_clusterOK = new TH1D("phiPhot_clusterOK", "", 15, -3.1416, 3.1416);
  h1_phiPhot_clusterOK->Sumw2();
  TH1D* h1_deltaPhi_clusterOK = new TH1D("deltaPhi_clusterOK", "", 15, 0., 3.1416);
  h1_deltaPhi_clusterOK->Sumw2();

  TH1D* h1_ptPhot_clusterOK_isolated = new TH1D("ptPhot_clusterOK_isolated", "", 28, ptPhot_binning[0], ptMax);
  h1_ptPhot_clusterOK_isolated->Sumw2();
  TH1D* h1_etaPhot_clusterOK_isolated = new TH1D("etaPhot_clusterOK_isolated", "", 15, -1.3, 1.3);
  h1_etaPhot_clusterOK_isolated->Sumw2();
  TH1D* h1_phiPhot_clusterOK_isolated = new TH1D("phiPhot_clusterOK_isolated", "", 15, -3.1416, 3.1416);
  h1_phiPhot_clusterOK_isolated->Sumw2();

  TH1D* h1_hcalIsoPhotReco_Nm1 = new TH1D("hcalIsoPhotReco_Nm1", "", 10, 0., 0.5);
  h1_hcalIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_ecalIsoPhotReco_Nm1 = new TH1D("ecalIsoPhotReco_Nm1", "", 20, 0., 1.);
  h1_ecalIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_ptTrkIsoPhotReco_Nm1 = new TH1D("ptTrkIsoPhotReco_Nm1", "", 20, 0., 1.);
  h1_ptTrkIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_nTrkIsoPhotReco_Nm1 = new TH1D("nTrkIsoPhotReco_Nm1", "", 11, -0.5, 10.5);
  h1_nTrkIsoPhotReco_Nm1->Sumw2();
  TH1D* h1_clusterMajPhotReco_Nm1 = new TH1D("clusterMajPhotReco_Nm1", "", 50, 0., 1.5);
  h1_clusterMajPhotReco_Nm1->Sumw2();
  TH1D* h1_clusterMinPhotReco_Nm1 = new TH1D("clusterMinPhotReco_Nm1", "", 50, 0., 0.6);
  h1_clusterMinPhotReco_Nm1->Sumw2();

  TH1D* h1_hcalIsoPhotReco_clusterOK = new TH1D("hcalIsoPhotReco_clusterOK", "", 10, 0., 0.5);
  h1_hcalIsoPhotReco_clusterOK->Sumw2();
  TH1D* h1_ecalIsoPhotReco_clusterOK = new TH1D("ecalIsoPhotReco_clusterOK", "", 20, 0., 1.);
  h1_ecalIsoPhotReco_clusterOK->Sumw2();
  TH1D* h1_ptTrkIsoPhotReco_clusterOK = new TH1D("ptTrkIsoPhotReco_clusterOK", "", 20, 0., 1.);
  h1_ptTrkIsoPhotReco_clusterOK->Sumw2();
  TH1D* h1_nTrkIsoPhotReco_clusterOK = new TH1D("nTrkIsoPhotReco_clusterOK", "", 11, -0.5, 10.5);
  h1_nTrkIsoPhotReco_clusterOK->Sumw2();

  TH1D* h1_clusterMajPhotReco_isolated = new TH1D("clusterMajPhotReco_isolated", "", 30, 0., 1.5);
  h1_clusterMajPhotReco_isolated->Sumw2();
  TH1D* h1_clusterMinPhotReco_isolated = new TH1D("clusterMinPhotReco_isolated", "", 20, 0., 0.6);
  h1_clusterMinPhotReco_isolated->Sumw2();

  TH1D* h1_clusterMajPhotReco= new TH1D("clusterMajPhotReco", "", 30, 0., 1.5);
  h1_clusterMajPhotReco->Sumw2();
  TH1D* h1_clusterMinPhotReco= new TH1D("clusterMinPhotReco", "", 20, 0., 0.6);
  h1_clusterMinPhotReco->Sumw2();

  TH1D* h1_deltaPhi_clusterOK_isolated = new TH1D("deltaPhi_clusterOK_isolated", "", 15, 0., 3.1416);
  h1_deltaPhi_clusterOK_isolated->Sumw2();
  TH1D* h1_ptSecondJetRel_clusterOK_isolated = new TH1D("ptSecondJetRel_clusterOK_isolated", "", 15, 0., 1.5);
  h1_ptSecondJetRel_clusterOK_isolated->Sumw2();
//TH1D* h1_response = new TH1D("response", "", 15, 0., 2.);
//h1_response->Sumw2();

//std::vector<TH1F*> h1_ptPhotMean;
//for( unsigned i=0; i<(ptPhot_binning.size()-1); ++i ) {
//  Float_t ptMin = ptPhot_binning[i];
//  Float_t ptMax = ptPhot_binning[i+1];
//  char histoName[100];
//  sprintf( histoName, "ptPhot_%.0f_%.0f", ptMin, ptMax );
//  TH1F* tmp = new Th1F(histoName, "", 


  Double_t ptPhotBinning_array[200]; //ugly! no more than 200 pt bins supported
  for( unsigned i=0; i<ptPhot_binning.size(); ++i )
    ptPhotBinning_array[i] = ptPhot_binning[i];


  TProfile* hp_ptJetGenMean = new TProfile("ptJetGenMean", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  TProfile* hp_ptPhotMean = new TProfile("ptPhotMean", "", ptPhot_binning.size()-1, ptPhotBinning_array);
  TProfile* hp_ptPhotMean_loose = new TProfile("ptPhotMean_loose", "", ptPhot_binning.size()-1, ptPhotBinning_array);

  std::vector<TH1F*> h1_response              = getResponseHistos("response");
  std::vector<TH1F*> h1_responseGEN           = getResponseHistos("responseGEN");
  std::vector<TH1F*> h1_responseMPF           = getResponseHistos("responseMPF");

  std::vector<TH1F*> h1_response_clusterOK    = getResponseHistos("response_clusterOK");
  std::vector<TH1F*> h1_responseGEN_clusterOK = getResponseHistos("responseGEN_clusterOK");

  std::vector<TH1F*> h1_response_loose        = getResponseHistos("response_loose");
  std::vector<TH1F*> h1_responseGEN_loose     = getResponseHistos("responseGEN_loose");



  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);
  Float_t eventWeight_loose;
  if( useGenJets )
    tree->SetBranchAddress("eventWeight_loose", &eventWeight_loose);
  Float_t eventWeight_medium;
  if( useGenJets )
    tree->SetBranchAddress("eventWeight_medium", &eventWeight_medium);

  Float_t epfMet;
  tree->SetBranchAddress("epfMet", &epfMet);
  Float_t phipfMet;
  tree->SetBranchAddress("phipfMet", &phipfMet);

  Bool_t isIsolated_hcal_loose;
  tree->SetBranchAddress("isIsolated_hcal_loose", &isIsolated_hcal_loose);
  Bool_t isIsolated_ecal_loose;
  tree->SetBranchAddress("isIsolated_ecal_loose", &isIsolated_ecal_loose);
  Bool_t isIsolated_ptTracks_loose;
  tree->SetBranchAddress("isIsolated_ptTracks_loose", &isIsolated_ptTracks_loose);
  Bool_t isIsolated_nTracks_loose;
  tree->SetBranchAddress("isIsolated_nTracks_loose", &isIsolated_nTracks_loose);
  Bool_t clusterMajOK_loose;
  tree->SetBranchAddress("clusterMajOK_loose", &clusterMajOK_loose);
  Bool_t clusterMinOK_loose;
  tree->SetBranchAddress("clusterMinOK_loose", &clusterMinOK_loose);

  Bool_t isIsolated_hcal_medium;
  tree->SetBranchAddress("isIsolated_hcal_medium", &isIsolated_hcal_medium);
  Bool_t isIsolated_ecal_medium;
  tree->SetBranchAddress("isIsolated_ecal_medium", &isIsolated_ecal_medium);
  Bool_t isIsolated_ptTracks_medium;
  tree->SetBranchAddress("isIsolated_ptTracks_medium", &isIsolated_ptTracks_medium);
  Bool_t isIsolated_nTracks_medium;
  tree->SetBranchAddress("isIsolated_nTracks_medium", &isIsolated_nTracks_medium);
  Bool_t clusterMajOK_medium;
  tree->SetBranchAddress("clusterMajOK_medium", &clusterMajOK_medium);
  Bool_t clusterMinOK_medium;
  tree->SetBranchAddress("clusterMinOK_medium", &clusterMinOK_medium);

  Bool_t passedPhotonID_loose;
  tree->SetBranchAddress("passedPhotonID_loose", &passedPhotonID_loose);
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





  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;

    bool jetInBarrel = (fabs(etaJetReco)<1.3);

    if( ptPhotReco<ptPhot_binning[0] ) continue;
    if( fabs(etaPhotReco)>1.3 ) continue;
    if( clusterMinPhotReco<0.15 && !MCassoc && !useGenJets ) continue; //protection vs EB spikes


    //first find correct photon pt bin:
    int theBin = hp_ptPhotMean->FindBin( ptPhotReco );
    theBin -= 1; //because arrays start from 0


    //leading jet and photon back2back in transverse plane
    bool back2back = true;
    Float_t deltaPhi_jet = fabs(fitTools::delta_phi(phiPhotReco, phiJetReco));
    Float_t pi = TMath::Pi();
    if( fabs(deltaPhi_jet) < (pi - 1.) ) back2back = false; //loose back to back for now

    bool secondJetOK = ( pt2ndJetReco < 0.5*ptPhotReco );
  //Float_t secondJetPt_thresh = (recoType=="calo") ? 5. : 8.;
  //bool secondJetOK = ( pt2ndJetReco < 0.1*ptPhotReco || pt2ndJetReco<secondJetPt_thresh );
  

    //before selection fill N-1 isolation plots (no event topology for isolation variables):
    if(                           isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_hcalIsoPhotReco_Nm1->Fill( hcalIsoPhotReco, eventWeight);
    if( isIsolated_hcal_medium                           && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_ecalIsoPhotReco_Nm1->Fill( ecalIsoPhotReco, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium                               && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium  ) h1_ptTrkIsoPhotReco_Nm1->Fill( ptTrkIsoPhotReco, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium                              && clusterMajOK_medium && clusterMinOK_medium  ) h1_nTrkIsoPhotReco_Nm1->Fill( nTrkIsoPhotReco, eventWeight);
    //no cluster cuts on cluster N-1's:
    if( isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium   ) h1_clusterMajPhotReco_Nm1->Fill( clusterMajPhotReco, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium   ) h1_clusterMinPhotReco_Nm1->Fill( clusterMinPhotReco, eventWeight);
    // yes topology for topology variables:
    if( isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium              && secondJetOK && jetInBarrel) h1_deltaPhi_Nm1->Fill( deltaPhi_jet, eventWeight);
    if( isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium && clusterMajOK_medium && clusterMinOK_medium && back2back                && jetInBarrel) h1_ptSecondJetRel_Nm1->Fill( pt2ndJetReco/ptPhotReco, eventWeight);


    bool isIsolated_loose = (isIsolated_hcal_loose && isIsolated_ecal_loose && isIsolated_ptTracks_loose && isIsolated_nTracks_loose);
    bool isIsolated_medium = (isIsolated_hcal_medium && isIsolated_ecal_medium && isIsolated_ptTracks_medium && isIsolated_nTracks_medium);
    bool clusterShapeOK_medium = (clusterMajOK_medium && clusterMinOK_medium );



    if( clusterShapeOK_medium ) {
      h1_hcalIsoPhotReco_clusterOK->Fill( hcalIsoPhotReco, eventWeight);
      h1_ecalIsoPhotReco_clusterOK->Fill( ecalIsoPhotReco, eventWeight);
      h1_ptTrkIsoPhotReco_clusterOK->Fill( ptTrkIsoPhotReco, eventWeight);
      h1_nTrkIsoPhotReco_clusterOK->Fill( nTrkIsoPhotReco, eventWeight);
      h1_deltaPhi_clusterOK->Fill( deltaPhi_jet, eventWeight);

      if( isIsolated_medium ) {
        h1_deltaPhi_clusterOK_isolated->Fill( deltaPhi_jet, eventWeight);
        h1_ptSecondJetRel_clusterOK_isolated->Fill( pt2ndJetReco/ptPhotReco, eventWeight);
        h1_phiPhot_clusterOK_isolated->Fill( phiPhotReco, eventWeight );
        h1_etaPhot_clusterOK_isolated->Fill( etaPhotReco, eventWeight );
        h1_ptPhot_clusterOK_isolated->Fill( ptPhotReco, eventWeight );
      }

      if( back2back && secondJetOK && jetInBarrel) {
        h1_ptPhot_clusterOK->Fill( ptPhotReco, eventWeight );
        h1_phiPhot_clusterOK->Fill( phiPhotReco, eventWeight );
        h1_etaPhot_clusterOK->Fill( etaPhotReco, eventWeight );
        h1_response_clusterOK[theBin]->Fill( ptJetReco/ptPhotReco, eventWeight );
        h1_responseGEN_clusterOK[theBin]->Fill( ptJetReco/ptJetGen, eventWeight );
      } //if back2back
    }

    if( isIsolated_medium ) {
      h1_clusterMajPhotReco_isolated->Fill( clusterMajPhotReco, eventWeight );
      h1_clusterMinPhotReco_isolated->Fill( clusterMinPhotReco, eventWeight );
    }

    h1_clusterMajPhotReco->Fill( clusterMajPhotReco, eventWeight );
    h1_clusterMinPhotReco->Fill( clusterMinPhotReco, eventWeight );


    bool photonOK_medium = passedPhotonID_medium || MCassoc || useGenJets;
    bool photonOK_loose = (isIsolated_loose && clusterShapeOK_medium) || MCassoc || useGenJets;

    //event selection (medium):
    if( photonOK_medium && back2back && jetInBarrel) {

      Float_t correctWeight = (useGenJets) ? eventWeight_medium : eventWeight;

    //// no cut in 2nd jet to produce gen response plot to have more stat
    //h1_responseGEN[theBin]->Fill( ptJetReco/ptJetGen, correctWeight );
      
      if( secondJetOK ) {
        //fill responseGEN histos before any selection to avoid biases:
        int theBinGEN = hp_ptPhotMean->FindBin( ptJetGen );
        theBinGEN -= 1; //because arrays start from 0
      

        //h1_responseGEN[theBin]->Fill( ptJetReco/ptJetGen, correctWeight );
        
        if( ptJetGen>ptPhot_binning[0] ) {
          h1_responseGEN[theBinGEN]->Fill( ptJetReco/ptJetGen, eventWeight );
          hp_ptJetGenMean->Fill( ptJetGen, ptJetGen, eventWeight );
        }

        h1_phiPhot->Fill( phiPhotReco, correctWeight );
        h1_etaPhot->Fill( etaPhotReco, correctWeight );
        h1_ptPhot->Fill( ptPhotReco, correctWeight );
        h1_ptJetReco->Fill( ptJetReco, correctWeight );
        h1_pt2ndJetReco->Fill( pt2ndJetReco, correctWeight );
        hp_ptPhotMean->Fill( ptPhotReco, ptPhotReco, correctWeight );
        h1_response[theBin]->Fill( ptJetReco/ptPhotReco, correctWeight );

        if( recoType=="pf" ) {
          // compute mpf :
          Float_t phi_Phot_Met = fitTools::delta_phi( phiPhotReco, phipfMet );
          Float_t mpfResponse = 1. + epfMet*ptPhotReco*cos( phi_Phot_Met ) / (ptPhotReco*ptPhotReco);
          h1_responseMPF[theBin]->Fill( mpfResponse, correctWeight );
        } //if pf
          
      }

    } 

    //event selection (loose):
    if( photonOK_loose && back2back && jetInBarrel) {

      Float_t correctWeight = (useGenJets) ? eventWeight_loose : eventWeight;
      
      // no cut in 2nd jet to produce gen response plot to have more stat
      h1_responseGEN_loose[theBin]->Fill( ptJetReco/ptJetGen, correctWeight );
      
    
      if( secondJetOK ) {
        //if( ptPhotReco>28. && ptJetReco/ptPhotReco<0.1 )
        //  std::cout << " ptPhot: " << ptPhotReco << "\tptJet: " << ptJetReco << "\tetaJet: " << etaJetReco << "\tphiJet: " << phiJetReco << "\tclustMin: " << clusterMinPhotReco << "\tclustMaj: " << clusterMajPhotReco << "\tdeltaPhi: " << deltaPhi_jet << "\tnTrkIso:" << nTrkIsoPhotReco << std::endl;
        h1_ptPhot_loose->Fill( ptPhotReco, correctWeight );
        h1_phiPhot_loose->Fill( phiPhotReco, correctWeight );
        h1_etaPhot_loose->Fill( etaPhotReco, correctWeight );
        hp_ptPhotMean_loose->Fill( ptPhotReco, ptPhotReco, correctWeight );
        h1_response_loose[theBin]->Fill( ptJetReco/ptPhotReco, correctWeight );
      }

    } //if loose


  } //for entries



  //TFile* outFile = TFile::Open("secondJet_yesPtHat.root", "recreate");
  
  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaPhotonJet_"+dataset;
  else {
   if(dataset!="") outfileName = "PhotonJet_"+dataset;
   else outfileName = "PhotonJet";
  }

  outfileName = outfileName + suffix;
  if( MCassoc ) outfileName = outfileName + "_MCassoc";
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h1_totalLumi->Write();

  h1_ptPhot->Write();
  h1_phiPhot->Write();
  h1_etaPhot->Write();
  h1_ptJetReco->Write();
  h1_pt2ndJetReco->Write();

  h1_ptPhot_loose->Write();
  h1_phiPhot_loose->Write();
  h1_etaPhot_loose->Write();

  h1_ptPhot_clusterOK->Write();
  h1_phiPhot_clusterOK->Write();
  h1_etaPhot_clusterOK->Write();

  h1_ptPhot_clusterOK_isolated->Write();
  h1_phiPhot_clusterOK_isolated->Write();
  h1_etaPhot_clusterOK_isolated->Write();


  h1_deltaPhi_Nm1->Write();
  h1_ptSecondJetRel_Nm1->Write();

  h1_deltaPhi_clusterOK_isolated->Write();
  h1_ptSecondJetRel_clusterOK_isolated->Write();

  h1_hcalIsoPhotReco_Nm1->Write();
  h1_ecalIsoPhotReco_Nm1->Write();
  h1_ptTrkIsoPhotReco_Nm1->Write();
  h1_nTrkIsoPhotReco_Nm1->Write();
  h1_clusterMajPhotReco_Nm1->Write();
  h1_clusterMinPhotReco_Nm1->Write();

  h1_hcalIsoPhotReco_clusterOK->Write();
  h1_ecalIsoPhotReco_clusterOK->Write();
  h1_ptTrkIsoPhotReco_clusterOK->Write();
  h1_nTrkIsoPhotReco_clusterOK->Write();
  h1_deltaPhi_clusterOK->Write();

  h1_clusterMajPhotReco_isolated->Write();
  h1_clusterMinPhotReco_isolated->Write();

  h1_clusterMajPhotReco->Write();
  h1_clusterMinPhotReco->Write();

  hp_ptJetGenMean->Write();
  hp_ptPhotMean->Write();
  hp_ptPhotMean_loose->Write();

  for( unsigned i=0; i<h1_response.size(); ++i ) {
    h1_response[i]->Write();
    h1_responseGEN[i]->Write();
    h1_response_loose[i]->Write();
    h1_responseGEN_loose[i]->Write();
    h1_response_clusterOK[i]->Write();
    h1_responseGEN_clusterOK[i]->Write();

    if( recoType=="pf" )
      h1_responseMPF[i]->Write();
  }


  outFile->Close();

  delete h1_hcalIsoPhotReco_Nm1;
  h1_hcalIsoPhotReco_Nm1 = 0;
  delete h1_ecalIsoPhotReco_Nm1;
  h1_ecalIsoPhotReco_Nm1 = 0;
  delete h1_ptTrkIsoPhotReco_Nm1;
  h1_ptTrkIsoPhotReco_Nm1 = 0;
  delete h1_nTrkIsoPhotReco_Nm1;
  h1_nTrkIsoPhotReco_Nm1 = 0;
  delete h1_clusterMajPhotReco_Nm1;
  h1_clusterMajPhotReco_Nm1 = 0;
  delete h1_clusterMinPhotReco_Nm1;
  h1_clusterMinPhotReco_Nm1 = 0;

  delete h1_ptJetReco;
  h1_ptJetReco = 0;
  delete h1_pt2ndJetReco;
  h1_pt2ndJetReco = 0;
  delete h1_ptPhot;
  h1_ptPhot = 0;
  delete h1_ptPhot;
  h1_ptPhot = 0;
  delete h1_etaPhot;
  h1_etaPhot = 0;
  delete h1_phiPhot;
  h1_phiPhot = 0;
  delete h1_deltaPhi_Nm1;
  h1_deltaPhi_Nm1 = 0;
  delete h1_ptSecondJetRel_Nm1;
  h1_ptSecondJetRel_Nm1 = 0;

  delete h1_hcalIsoPhotReco_clusterOK;
  h1_hcalIsoPhotReco_clusterOK = 0;
  delete h1_ecalIsoPhotReco_clusterOK;
  h1_ecalIsoPhotReco_clusterOK = 0;
  delete h1_ptTrkIsoPhotReco_clusterOK;
  h1_ptTrkIsoPhotReco_clusterOK = 0;
  delete h1_nTrkIsoPhotReco_clusterOK;
  h1_nTrkIsoPhotReco_clusterOK = 0;

  delete h1_clusterMajPhotReco_isolated;
  h1_clusterMajPhotReco_isolated = 0;
  delete h1_clusterMinPhotReco_isolated;
  h1_clusterMinPhotReco_isolated = 0;
//delete h1_response;
//h1_response = 0;

}
  
std::vector<TH1F*> getResponseHistos(const std::string& name) {

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();
  std::vector<TH1F*> returnVector;

  for( unsigned i=0; i<(ptPhot_binning.size()-1); ++i ) {
    char histoName[100];
    sprintf( histoName, "%s_ptPhot_%.0f_%.0f", name.c_str(), ptPhot_binning[i], ptPhot_binning[i+1]);
    TH1F* newHisto = new TH1F(histoName, "", 15, 0., 2.);
    newHisto->Sumw2();
    returnVector.push_back(newHisto);
  }

  return returnVector;

}


