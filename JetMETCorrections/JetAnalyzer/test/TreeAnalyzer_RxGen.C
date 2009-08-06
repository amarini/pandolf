#define TreeAnalyzer_RxGen_cxx
#include "TreeAnalyzer_RxGen.h"
#include <iostream>





void TreeAnalyzer_RxGen::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L TreeAnalyzer_RxGen.C
//      Root > TreeAnalyzer_RxGen t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;



   Long64_t nentries;

   if( DEBUG_ ) nentries = 10000;
   else nentries = fChain->GetEntries();

   BookStuff();



   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%10000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;



     // ------------  ITERATIVE CONE 0.5 JETS

     for(unsigned int iGenJet=0; iGenJet<nJetGen_itCone5; ++iGenJet) {

       Float_t eJet = eJetGen_itCone5[iGenJet];
       Float_t etaJet = etaJetGen_itCone5[iGenJet];

       Float_t eTracks = eTracksGen_itCone5[iGenJet];
       Float_t ePhotons = ePhotonsGen_itCone5[iGenJet];
       Float_t eNeutralHadrons = eNeutralHadronsGen_itCone5[iGenJet];
       Float_t eMuons = eMuonsGen_itCone5[iGenJet];

       Float_t eNeutron= eNeutronGen_itCone5[iGenJet];
       Float_t eK0L= eK0LGen_itCone5[iGenJet];
       Float_t eK0S= eK0SGen_itCone5[iGenJet];
       Float_t eLambda= eLambdaGen_itCone5[iGenJet];
       Float_t eCsi= eCsiGen_itCone5[iGenJet];

       Float_t  Rch = eTracks/eJet;
       Float_t  Rgamma = ePhotons/eJet;
       Float_t  Rnh = eNeutralHadrons/eJet;
       Float_t  Rmu = eMuons/eJet;

       Float_t Rnh_neutron = (eNeutralHadrons>0.) ? eNeutron/eNeutralHadrons : -1.;
       Float_t Rnh_K0L = (eNeutralHadrons>0.) ? eK0L/eNeutralHadrons : -1.;
       Float_t Rnh_K0S = (eNeutralHadrons>0.) ? eK0S/eNeutralHadrons : -1.;
       Float_t Rnh_lambda = (eNeutralHadrons>0.) ? eLambda/eNeutralHadrons : -1.;
       Float_t Rnh_csi = (eNeutralHadrons>0.) ? eCsi/eNeutralHadrons : -1.;
  
       hp_RchGen_vs_eta_->Fill(etaJet, Rch);
       hp_RgammaGen_vs_eta_->Fill(etaJet, Rgamma);
       hp_RnhGen_vs_eta_->Fill(etaJet, Rnh);
       hp_RmuGen_vs_eta_->Fill(etaJet, Rmu);
       
       hp_Rnh_neutron_vs_eta_->Fill(etaJet, Rnh_neutron);
       hp_Rnh_K0L_vs_eta_->Fill(etaJet, Rnh_K0L);
       hp_Rnh_K0S_vs_eta_->Fill(etaJet, Rnh_K0S);
       hp_Rnh_lambda_vs_eta_->Fill(etaJet, Rnh_lambda);
       hp_Rnh_csi_vs_eta_->Fill(etaJet, Rnh_csi);

       h1_RchGen_->Fill(Rch);
       h1_RgammaGen_->Fill(Rgamma);
       h1_RnhGen_->Fill(Rnh);
       h1_RmuGen_->Fill(Rmu);

       if( fabs(etaJet)<1.4 ) {
         h1_RchGen_barrel_->Fill(Rch);
         h1_RgammaGen_barrel_->Fill(Rgamma);
         h1_RnhGen_barrel_->Fill(Rnh);
         h1_RmuGen_barrel_->Fill(Rmu);
       } else if( fabs(etaJet)<3. ) {
         h1_RchGen_endcap_->Fill(Rch);
         h1_RgammaGen_endcap_->Fill(Rgamma);
         h1_RnhGen_endcap_->Fill(Rnh);
         h1_RmuGen_endcap_->Fill(Rmu);
       }

     } // for gen jets

   } //for entries

} //loop



void TreeAnalyzer_RxGen::BookStuff() {

  outfile_->cd();

  hp_RchGen_vs_eta_ = new TProfile("RchGen_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_RgammaGen_vs_eta_ = new TProfile("RgammaGen_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_RnhGen_vs_eta_ = new TProfile("RnhGen_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_RmuGen_vs_eta_ = new TProfile("RmuGen_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);

  hp_Rnh_neutron_vs_eta_ = new TProfile("Rnh_neutron_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_K0L_vs_eta_ = new TProfile("Rnh_K0L_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_K0S_vs_eta_ = new TProfile("Rnh_K0S_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_lambda_vs_eta_ = new TProfile("Rnh_lambda_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);
  hp_Rnh_csi_vs_eta_ = new TProfile("Rnh_csi_vs_eta", "", 100, -5.5, 5.5, 0., 1000.);

  h1_RchGen_= new TH1F("RchGen", "", 100, 0., 1.1);
  h1_RgammaGen_= new TH1F("RgammaGen", "", 100, 0., 1.1);
  h1_RnhGen_= new TH1F("RnhGen", "", 100, 0., 1.1);
  h1_RmuGen_= new TH1F("RmuGen", "", 100, 0., 1.1);

  h1_RchGen_barrel_= new TH1F("RchGen_barrel", "", 100, 0., 1.1);
  h1_RgammaGen_barrel_= new TH1F("RgammaGen_barrel", "", 100, 0., 1.1);
  h1_RnhGen_barrel_= new TH1F("RnhGen_barrel", "", 100, 0., 1.1);
  h1_RmuGen_barrel_= new TH1F("RmuGen_barrel", "", 100, 0., 1.1);

  h1_RchGen_endcap_= new TH1F("RchGen_endcap", "", 100, 0., 1.1);
  h1_RgammaGen_endcap_= new TH1F("RgammaGen_endcap", "", 100, 0., 1.1);
  h1_RnhGen_endcap_= new TH1F("RnhGen_endcap", "", 100, 0., 1.1);
  h1_RmuGen_endcap_= new TH1F("RmuGen_endcap", "", 100, 0., 1.1);



} //bookstuff

