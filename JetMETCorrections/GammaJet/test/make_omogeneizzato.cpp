#include <string>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"




int main() {

  std::string dataset = "HT_Run2011_FULL";

  std::string infileName = "DiJet_" + dataset + ".root";
  TFile* infile = TFile::Open(infileName.c_str());
  TTree* chain = (TTree*)infile->Get("tree_passedEvents");

  std::cout << "-> Opened file: " << infileName << std::endl;
  //chain->Add( "DiJet_HT_Run2011A-May10thRereco-v1_HLT.root/tree_passedEvents" );
  //chain->Add( "DiJet_HT_Run2011A-PromptReco-v4_HLT.root/tree_passedEvents" );
  //chain->Add( "DiJet_HT_Run2011A-PromptReco-v6_HLT.root/tree_passedEvents" );
  //chain->Add( "DiJet_HT_Run2011B-PromptReco-v1_HLT.root/tree_passedEvents" );

  


  Float_t rhoPF;
  chain->SetBranchAddress( "rhoPF", &rhoPF );

  Float_t ptJet0;
  chain->SetBranchAddress( "ptJet0", &ptJet0 );
  Float_t ptJet1;
  chain->SetBranchAddress( "ptJet1", &ptJet1 );
  Float_t ptJet2;
  chain->SetBranchAddress( "ptJet2", &ptJet2 );

  Float_t etaJet0;
  chain->SetBranchAddress( "etaJet0", &etaJet0 );
  Float_t etaJet1;
  chain->SetBranchAddress( "etaJet1", &etaJet1 );

  Int_t nChargedJet0;
  chain->SetBranchAddress( "nChargedJet0", &nChargedJet0 );
  Int_t nChargedJet1;
  chain->SetBranchAddress( "nChargedJet1", &nChargedJet1 );

  Int_t nNeutralJet0;
  chain->SetBranchAddress( "nNeutralJet0", &nNeutralJet0 );
  Int_t nNeutralJet1;
  chain->SetBranchAddress( "nNeutralJet1", &nNeutralJet1 );

  Float_t ptDJet0;
  chain->SetBranchAddress( "ptDJet0", &ptDJet0 );
  Float_t ptDJet1;
  chain->SetBranchAddress( "ptDJet1", &ptDJet1 );

  Float_t QGLikelihoodJet0;
  chain->SetBranchAddress( "QGLikelihoodJet0", &QGLikelihoodJet0 );
  Float_t QGLikelihoodJet1;
  chain->SetBranchAddress( "QGLikelihoodJet1", &QGLikelihoodJet1 );

  Float_t ht_akt5;
  chain->SetBranchAddress( "ht_akt5", &ht_akt5 );
  Bool_t passed_HT150;
  chain->SetBranchAddress( "passed_HT150", &passed_HT150);
  Bool_t passed_HT250;
  chain->SetBranchAddress( "passed_HT250", &passed_HT250);
  Bool_t passed_HT350;
  chain->SetBranchAddress( "passed_HT350", &passed_HT350);
  Bool_t passed_HT400;
  chain->SetBranchAddress( "passed_HT400", &passed_HT400);
  Bool_t passed_HT500;
  chain->SetBranchAddress( "passed_HT500", &passed_HT500);
  Bool_t passed_HT600;
  chain->SetBranchAddress( "passed_HT600", &passed_HT600);


  std::string outfileName = "Omog_" + dataset + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  outfile->cd();

  Float_t ptJet0_out;
  Int_t nChargedJet0_out;
  Int_t nNeutralJet0_out;
  Float_t ptDJet0_out;
  Float_t QGLikelihoodJet0_out;
  TTree* tree_omogeneizzato = new TTree("omogeneizzato", "");
  tree_omogeneizzato->Branch( "rhoPF", &rhoPF, "rhoPF/F" );
  tree_omogeneizzato->Branch( "ptJet0", &ptJet0_out, "ptJet0_out/F" );
  tree_omogeneizzato->Branch( "nChargedJet0", &nChargedJet0_out, "nChargedJet0_out/I" );
  tree_omogeneizzato->Branch( "nNeutralJet0", &nNeutralJet0_out, "nNeutralJet0_out/I" );
  tree_omogeneizzato->Branch( "ptDJet0", &ptDJet0_out, "ptDJet0_out/F" );
  tree_omogeneizzato->Branch( "QGLikelihoodJet0", &QGLikelihoodJet0_out, "QGLikelihoodJet0_out/F" );

  
  int nEntries = chain->GetEntries();

  for( unsigned iEntry=0; iEntry<nEntries; ++iEntry ) {

    if( iEntry % 500000 == 0 ) std::cout << "Entry: " << iEntry << " / " << nEntries << std::endl;

    chain->GetEntry(iEntry);


    if( ptJet2 > 0.2*(ptJet0+ptJet1)/2. ) continue; //dijet selection


    if( passed_HT150 ) {

      if( ht_akt5 < 160. ) continue;

      bool firstJetOK =  ( ptJet0>50. && ptJet0<100. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>50. && ptJet1<100. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT150


    else if( passed_HT250 ) {

      if( ht_akt5 < 265. ) continue;

      bool firstJetOK = ( ptJet0>100. && ptJet0<150. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>100. && ptJet1<150.&& fabs(etaJet1)<2.  );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT250


    else if( passed_HT350 ) {

      if( ht_akt5 < 365. ) continue;

      bool firstJetOK =  ( ptJet0>150. && ptJet0<200. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>150. && ptJet1<200. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT350


    else if( passed_HT400 ) {

      if( ht_akt5 < 420. ) continue;

      bool firstJetOK =  ( ptJet0>200. && ptJet0<250. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>200. && ptJet1<250. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT400


    else if( passed_HT500 ) {

      if( ht_akt5 < 525. ) continue;

      bool firstJetOK =  ( ptJet0>250. && ptJet0<300. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>250. && ptJet1<300. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT500


    else if( passed_HT600 ) {

      if( ht_akt5 < 640. ) continue;

      bool firstJetOK =  ( ptJet0>300. && fabs(etaJet0)<2. );
      bool secondJetOK = ( ptJet1>300. && fabs(etaJet1)<2. );

      if( firstJetOK ) {
  
        ptJet0_out = ptJet0;
        nChargedJet0_out = nChargedJet0;
        nNeutralJet0_out = nNeutralJet0;
        ptDJet0_out = ptDJet0;
        QGLikelihoodJet0_out = QGLikelihoodJet0;

        tree_omogeneizzato->Fill();

      }

      if( secondJetOK ) {
  
        ptJet0_out = ptJet1;
        nChargedJet0_out = nChargedJet1;
        nNeutralJet0_out = nNeutralJet1;
        ptDJet0_out = ptDJet1;
        QGLikelihoodJet0_out = QGLikelihoodJet1;

        tree_omogeneizzato->Fill();

      }

    } // if HT500

  } // for entries


  outfile->cd();
  
  tree_omogeneizzato->Write();

  outfile->Close();

  return 0;

}
