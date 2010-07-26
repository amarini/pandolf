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
bool MCASSOC_;


//Float_t computeInvMass( e, pt, eta, phi, e2nd, pt2nd, eta2nd, phi2nd);


void finalize(const std::string& dataset, std::string recoType, std::string jetAlgo, bool MCassoc=false) {

  RECOTYPE_ = recoType;
  MCASSOC_ = MCassoc;

  std::string algoType = (recoType=="calo") ? jetAlgo : recoType+jetAlgo;
  if( recoType=="jpt"&&jetAlgo=="akt5" ) algoType="jptak5"; 

  TChain* tree = new TChain("jetTree");
  
  std::string infileName, treeName;
  TFile* infile=0;
  TH1F* h1_lumi=0;
  Double_t totalLumi=0.;

  std::string prefix = "HZZlljj_2ndLevelTree_";
  std::string suffix = "_"+algoType;
  std::string suffix_QCD = "";


  if( dataset=="MC_PhotonJet_Summer09") {

    infileName = prefix + "MC_PhotonJet_Summer09_Pt0to15" + suffix + ".root";
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
    infileName = prefix + "MC_PhotonJet_Summer09_Pt15to20" + suffix + ".root";
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
    infileName = prefix + "MC_PhotonJet_Summer09_Pt20to30" + suffix + ".root";
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

  Float_t ptMax = 60.;


  float dijetmass_max = (dataset=="H500ZZlljj_gluonfusion") ? 200. : 150;
  TH1F* h1_diJetMassCorr = new TH1F("diJetMassCorr", "", 100, 40., dijetmass_max);

  float invMass_min, invMass_max;
  if( dataset=="H130ZZlljj_gluonfusion" ) {
    invMass_min = 80.;
    invMass_max = 220.;
  } else if( dataset=="H150ZZlljj_gluonfusion" ) {
    invMass_min = 80.;
    invMass_max = 220.;
  } else if( dataset=="H500ZZlljj_gluonfusion" ) {
    invMass_min = 300.;
    invMass_max = 800.;
  }
  TH1F* h1_higgsMassCorr = new TH1F("higgsMassCorr", "", 100, invMass_min, invMass_max);


  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);


  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);


  Float_t eZ;
  tree->SetBranchAddress("eZ", &eZ);
  Float_t ptZ;
  tree->SetBranchAddress("ptZ", &ptZ);
  Float_t etaZ;
  tree->SetBranchAddress("etaZ", &etaZ);
  Float_t phiZ;
  tree->SetBranchAddress("phiZ", &phiZ);

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


  Float_t e2ndJetReco;
  tree->SetBranchAddress("e2ndJetReco", &e2ndJetReco);
  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t ptCorr2ndJetReco;
  tree->SetBranchAddress("ptCorr2ndJetReco", &ptCorr2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);
  Float_t eTracks2ndReco;
  tree->SetBranchAddress("eTracks2ndReco", &eTracks2ndReco);

  Float_t e2ndJetGen;
  tree->SetBranchAddress("e2ndJetGen", &e2ndJetGen);
  Float_t pt2ndJetGen;
  tree->SetBranchAddress("pt2ndJetGen", &pt2ndJetGen);
  Float_t eta2ndJetGen;
  tree->SetBranchAddress("eta2ndJetGen", &eta2ndJetGen);
  Float_t phi2ndJetGen;
  tree->SetBranchAddress("phi2ndJetGen", &phi2ndJetGen);





  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    if( eventWeight <= 0. ) eventWeight = 1.;
    if( recoType=="pf" ) {
      if( eTracksReco < 0.1 ) continue;
      if( eTracks2ndReco < 0.1 ) continue;
    }

    Float_t eCorrJetReco = eJetReco*ptCorrJetReco/ptJetReco;

    Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
    Float_t pCorrJetReco  = ptCorrJetReco/sin(thetaJetReco);
    Float_t pxCorrJetReco = pCorrJetReco*cos(phiJetReco)*sin(thetaJetReco);
    Float_t pyCorrJetReco = pCorrJetReco*sin(phiJetReco)*sin(thetaJetReco);
    Float_t pzCorrJetReco = pCorrJetReco*cos(thetaJetReco);

    Float_t eCorr2ndJetReco = e2ndJetReco*ptCorr2ndJetReco/pt2ndJetReco;
    Float_t theta2ndJetReco = 2.*atan(exp(-eta2ndJetReco));
    Float_t pCorr2ndJetReco  = ptCorr2ndJetReco/sin(theta2ndJetReco);
    Float_t pxCorr2ndJetReco = pCorr2ndJetReco*cos(phi2ndJetReco)*sin(theta2ndJetReco);
    Float_t pyCorr2ndJetReco = pCorr2ndJetReco*sin(phi2ndJetReco)*sin(theta2ndJetReco);
    Float_t pzCorr2ndJetReco = pCorr2ndJetReco*cos(theta2ndJetReco);

    Float_t eCorrDiJet  =  eCorrJetReco+eCorr2ndJetReco;
    Float_t pxCorrDiJet = pxCorrJetReco+pxCorr2ndJetReco;
    Float_t pyCorrDiJet = pyCorrJetReco+pyCorr2ndJetReco;
    Float_t pzCorrDiJet = pzCorrJetReco+pzCorr2ndJetReco;

    Float_t diJetMassCorr = sqrt( eCorrDiJet*eCorrDiJet - pxCorrDiJet*pxCorrDiJet - pyCorrDiJet*pyCorrDiJet - pzCorrDiJet*pzCorrDiJet );

    h1_diJetMassCorr->Fill( diJetMassCorr );

    Float_t thetaZ = 2.*atan(exp(-etaZ));
    Float_t pZ  = ptZ/sin(thetaZ);
    Float_t pxZ = pZ*cos(phiZ)*sin(thetaZ);
    Float_t pyZ = pZ*sin(phiZ)*sin(thetaZ);
    Float_t pzZ = pZ*cos(thetaZ);

    Float_t  eHiggs =  eCorrDiJet+eZ;
    Float_t pxHiggs = pxCorrDiJet+pxZ;
    Float_t pyHiggs = pyCorrDiJet+pyZ;
    Float_t pzHiggs = pzCorrDiJet+pzZ;
    
    Float_t higgsMassCorr = sqrt( eHiggs*eHiggs - pxHiggs*pxHiggs - pyHiggs*pyHiggs - pzHiggs*pzHiggs );

    h1_higgsMassCorr->Fill( higgsMassCorr );

  } //for entries



  //TFile* outFile = TFile::Open("secondJet_yesPtHat.root", "recreate");
  
  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaHZZlljj_"+dataset;
  else {
   if(dataset!="") outfileName = "HZZlljj_"+dataset;
   else outfileName = "HZZlljj_";
  }

  outfileName = outfileName + suffix;
  if( MCASSOC_ ) outfileName = outfileName + "_MCassoc";
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();


  h1_diJetMassCorr->Write();
  h1_higgsMassCorr->Write();

  outFile->Close();

  delete h1_diJetMassCorr;
  h1_diJetMassCorr = 0;
  delete h1_higgsMassCorr;
  h1_higgsMassCorr = 0;

}
  


/*Float_t computeInvMass( e, pt, eta, phi, e2nd, pt2nd, eta2nd, phi2nd) {

    Float_t theta= 2.*atan(exp(-eta));
    Float_t p= pt/sin(theta);
    Float_t px= p*cos(phi)*sin(theta);
    Float_t py= p*sin(phi)*sin(theta);
    Float_t pz= p*cos(theta);

    Float_t e2nd= e2nd*pt2nd/pt2nd;
    Float_t theta2nd= 2.*atan(exp(-eta2nd));
    Float_t p2nd= pt2nd/sin(theta2nd);
    Float_t px2nd= p2nd*cos(phi2nd)*sin(theta2nd);
    Float_t py2nd= p2nd*sin(phi2nd)*sin(theta2nd);
    Float_t pz2nd= p2nd*cos(theta2nd);

    Float_t  ePair  =  e+ e2nd;
    Float_t pxPair = px+px2nd;
    Float_t pyPair = py+py2nd;
    Float_t pzPair = pz+pz2nd;

    Float_t mass = sqrt( ePair*ePair - pxPair*pxPair - pyPair*pyPair - pzPair*pzPair );

    return mass;

}*/
