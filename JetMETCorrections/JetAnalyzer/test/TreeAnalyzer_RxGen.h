//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 28 18:20:31 2009 by ROOT version 5.18/00a
// from TChain pippo/
//////////////////////////////////////////////////////////

#ifndef TreeAnalyzer_RxGen_h
#define TreeAnalyzer_RxGen_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TProfile.h>
#include <iostream>
#include <cmath>



bool DEBUG_ = false;
bool DEBUG_VERBOSE_ = false;



class TreeAnalyzer_RxGen {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           nJetGen_itCone5;
   Float_t         pxJetGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         pyJetGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         pzJetGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eJetGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         etaJetGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         phiJetGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eTracksGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         ePhotonsGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eMuonsGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eNeutralHadronsGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eNeutronGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eK0LGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eK0SGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eLambdaGen_itCone5[40];   //[nJetGen_itCone5]
   Float_t         eCsiGen_itCone5[40];   //[nJetGen_itCone5]
   Int_t           nJetGen_sisCone5;
   Float_t         pxJetGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         pyJetGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         pzJetGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         eJetGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         etaJetGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         phiJetGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         eTracksGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         ePhotonsGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         eMuonsGen_sisCone5[40];   //[nJetGen_sisCone5]
   Float_t         eNeutralHadronsGen_sisCone5[40];   //[nJetGen_sisCone5]
   Int_t           nJetGen_kt4;
   Float_t         pxJetGen_kt4[40];   //[nJetGen_kt4]
   Float_t         pyJetGen_kt4[40];   //[nJetGen_kt4]
   Float_t         pzJetGen_kt4[40];   //[nJetGen_kt4]
   Float_t         eJetGen_kt4[40];   //[nJetGen_kt4]
   Float_t         etaJetGen_kt4[40];   //[nJetGen_kt4]
   Float_t         phiJetGen_kt4[40];   //[nJetGen_kt4]
   Float_t         eTracksGen_kt4[40];   //[nJetGen_kt4]
   Float_t         ePhotonsGen_kt4[40];   //[nJetGen_kt4]
   Float_t         eMuonsGen_kt4[40];   //[nJetGen_kt4]
   Float_t         eNeutralHadronsGen_kt4[40];   //[nJetGen_kt4]


   TFile* outfile_;

   TProfile* hp_RchGen_vs_eta_;
   TProfile* hp_RgammaGen_vs_eta_;
   TProfile* hp_RnhGen_vs_eta_;
   TProfile* hp_RmuGen_vs_eta_;

   TProfile* hp_Rnh_neutron_vs_eta_;
   TProfile* hp_Rnh_K0L_vs_eta_;
   TProfile* hp_Rnh_K0S_vs_eta_;
   TProfile* hp_Rnh_lambda_vs_eta_;
   TProfile* hp_Rnh_csi_vs_eta_;

   TH1F* h1_RchGen_;
   TH1F* h1_RgammaGen_;
   TH1F* h1_RnhGen_;
   TH1F* h1_RmuGen_;

   TH1F* h1_RchGen_barrel_;
   TH1F* h1_RgammaGen_barrel_;
   TH1F* h1_RnhGen_barrel_;
   TH1F* h1_RmuGen_barrel_;

   TH1F* h1_RchGen_endcap_;
   TH1F* h1_RgammaGen_endcap_;
   TH1F* h1_RnhGen_endcap_;
   TH1F* h1_RmuGen_endcap_;


   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_nJetGen_itCone5;   //!
   TBranch        *b_pxJetGen_itCone5;   //!
   TBranch        *b_pyJetGen_itCone5;   //!
   TBranch        *b_pzJetGen_itCone5;   //!
   TBranch        *b_eJetGen_itCone5;   //!
   TBranch        *b_etaJetGen_itCone5;   //!
   TBranch        *b_phiJetGen_itCone5;   //!
   TBranch        *b_eTracksGen_itCone5;   //!
   TBranch        *b_ePhotonsGen_itCone5;   //!
   TBranch        *b_eMuonsGen_itCone5;   //!
   TBranch        *b_eNeutralHadronsGen_itCone5;   //!
   TBranch        *b_eNeutronGen_itCone5;   //!
   TBranch        *b_eK0LGen_itCone5;   //!
   TBranch        *b_eK0SGen_itCone5;   //!
   TBranch        *b_eLambdaGen_itCone5;   //!
   TBranch        *b_eCsiGen_itCone5;   //!
   TBranch        *b_nJetGen_sisCone5;   //!
   TBranch        *b_pxJetGen_sisCone5;   //!
   TBranch        *b_pyJetGen_sisCone5;   //!
   TBranch        *b_pzJetGen_sisCone5;   //!
   TBranch        *b_eJetGen_sisCone5;   //!
   TBranch        *b_etaJetGen_sisCone5;   //!
   TBranch        *b_phiJetGen_sisCone5;   //!
   TBranch        *b_eTracksGen_sisCone5;   //!
   TBranch        *b_ePhotonsGen_sisCone5;   //!
   TBranch        *b_eMuonsGen_sisCone5;   //!
   TBranch        *b_eNeutralHadronsGen_sisCone5;   //!
   TBranch        *b_nJetGen_kt4;   //!
   TBranch        *b_pxJetGen_kt4;   //!
   TBranch        *b_pyJetGen_kt4;   //!
   TBranch        *b_pzJetGen_kt4;   //!
   TBranch        *b_eJetGen_kt4;   //!
   TBranch        *b_etaJetGen_kt4;   //!
   TBranch        *b_phiJetGen_kt4;   //!
   TBranch        *b_eTracksGen_kt4;   //!
   TBranch        *b_ePhotonsGen_kt4;   //!
   TBranch        *b_eMuonsGen_kt4;   //!
   TBranch        *b_eNeutralHadronsGen_kt4;   //!

   TreeAnalyzer_RxGen(TTree* tree=0);
   virtual ~TreeAnalyzer_RxGen();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void BookStuff();
   //void MatchReco2Gen(std::string algoType);

};

#endif

#ifdef TreeAnalyzer_RxGen_cxx
TreeAnalyzer_RxGen::TreeAnalyzer_RxGen(TTree* tree)
{


   if (tree == 0) {

#ifdef SINGLE_TREE

      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("pippo");

#else // SINGLE_TREE

      TChain * chain = new TChain("pippo","");

      
      for(int i=1; i<100; ++i) {
        char treePath[60];
        sprintf(treePath, "QCD_2_2_6_NH/output_%d.root/myanalysis/pippo", i); 
        chain->Add(treePath);
      }
/* 
      for(int i=1; i<2304; ++i) {
      //for(int i=1; i<2; ++i) {
        char treePath[60];
        sprintf(treePath, "QCD_2_2_3/output_%d.root/myanalysis/pippo", i); 
        chain->Add(treePath);
      }
*/

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);

  
   char outfileName[60];

   if( DEBUG_ ) sprintf(outfileName, "prova2ndLevel.root");
   else sprintf(outfileName, "RxGenOutputFile.root");
   outfile_ = new TFile(outfileName, "RECREATE");
   

}

TreeAnalyzer_RxGen::~TreeAnalyzer_RxGen()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   outfile_->cd();

   hp_RchGen_vs_eta_->Write();
   hp_RgammaGen_vs_eta_->Write();
   hp_RnhGen_vs_eta_->Write();
   hp_RmuGen_vs_eta_->Write();

   hp_Rnh_neutron_vs_eta_->Write();
   hp_Rnh_K0L_vs_eta_->Write();
   hp_Rnh_K0S_vs_eta_->Write();
   hp_Rnh_lambda_vs_eta_->Write();
   hp_Rnh_csi_vs_eta_->Write();

   h1_RchGen_->Write();
   h1_RgammaGen_->Write();
   h1_RnhGen_->Write();
   h1_RmuGen_->Write();

   h1_RchGen_barrel_->Write();
   h1_RgammaGen_barrel_->Write();
   h1_RnhGen_barrel_->Write();
   h1_RmuGen_barrel_->Write();

   h1_RchGen_endcap_->Write();
   h1_RgammaGen_endcap_->Write();
   h1_RnhGen_endcap_->Write();
   h1_RmuGen_endcap_->Write();


   outfile_->Write();
   outfile_->Close();

}

Int_t TreeAnalyzer_RxGen::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeAnalyzer_RxGen::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeAnalyzer_RxGen::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nJetGen_itCone5", &nJetGen_itCone5, &b_nJetGen_itCone5);
   fChain->SetBranchAddress("pxJetGen_itCone5", pxJetGen_itCone5, &b_pxJetGen_itCone5);
   fChain->SetBranchAddress("pyJetGen_itCone5", pyJetGen_itCone5, &b_pyJetGen_itCone5);
   fChain->SetBranchAddress("pzJetGen_itCone5", pzJetGen_itCone5, &b_pzJetGen_itCone5);
   fChain->SetBranchAddress("eJetGen_itCone5", eJetGen_itCone5, &b_eJetGen_itCone5);
   fChain->SetBranchAddress("etaJetGen_itCone5", etaJetGen_itCone5, &b_etaJetGen_itCone5);
   fChain->SetBranchAddress("phiJetGen_itCone5", phiJetGen_itCone5, &b_phiJetGen_itCone5);
   fChain->SetBranchAddress("eTracksGen_itCone5", eTracksGen_itCone5, &b_eTracksGen_itCone5);
   fChain->SetBranchAddress("ePhotonsGen_itCone5", ePhotonsGen_itCone5, &b_ePhotonsGen_itCone5);
   fChain->SetBranchAddress("eMuonsGen_itCone5", eMuonsGen_itCone5, &b_eMuonsGen_itCone5);
   fChain->SetBranchAddress("eNeutralHadronsGen_itCone5", eNeutralHadronsGen_itCone5, &b_eNeutralHadronsGen_itCone5);
   fChain->SetBranchAddress("eNeutronGen_itCone5", eNeutronGen_itCone5, &b_eNeutronGen_itCone5);
   fChain->SetBranchAddress("eK0LGen_itCone5", eK0LGen_itCone5, &b_eK0LGen_itCone5);
   fChain->SetBranchAddress("eK0SGen_itCone5", eK0SGen_itCone5, &b_eK0SGen_itCone5);
   fChain->SetBranchAddress("eLambdaGen_itCone5", eLambdaGen_itCone5, &b_eLambdaGen_itCone5);
   fChain->SetBranchAddress("eCsiGen_itCone5", eCsiGen_itCone5, &b_eCsiGen_itCone5);
   fChain->SetBranchAddress("nJetGen_sisCone5", &nJetGen_sisCone5, &b_nJetGen_sisCone5);
   fChain->SetBranchAddress("pxJetGen_sisCone5", pxJetGen_sisCone5, &b_pxJetGen_sisCone5);
   fChain->SetBranchAddress("pyJetGen_sisCone5", pyJetGen_sisCone5, &b_pyJetGen_sisCone5);
   fChain->SetBranchAddress("pzJetGen_sisCone5", pzJetGen_sisCone5, &b_pzJetGen_sisCone5);
   fChain->SetBranchAddress("eJetGen_sisCone5", eJetGen_sisCone5, &b_eJetGen_sisCone5);
   fChain->SetBranchAddress("etaJetGen_sisCone5", etaJetGen_sisCone5, &b_etaJetGen_sisCone5);
   fChain->SetBranchAddress("phiJetGen_sisCone5", phiJetGen_sisCone5, &b_phiJetGen_sisCone5);
   fChain->SetBranchAddress("eTracksGen_sisCone5", eTracksGen_sisCone5, &b_eTracksGen_sisCone5);
   fChain->SetBranchAddress("ePhotonsGen_sisCone5", ePhotonsGen_sisCone5, &b_ePhotonsGen_sisCone5);
   fChain->SetBranchAddress("eMuonsGen_sisCone5", eMuonsGen_sisCone5, &b_eMuonsGen_sisCone5);
   fChain->SetBranchAddress("eNeutralHadronsGen_sisCone5", eNeutralHadronsGen_sisCone5, &b_eNeutralHadronsGen_sisCone5);
   fChain->SetBranchAddress("nJetGen_kt4", &nJetGen_kt4, &b_nJetGen_kt4);
   fChain->SetBranchAddress("pxJetGen_kt4", pxJetGen_kt4, &b_pxJetGen_kt4);
   fChain->SetBranchAddress("pyJetGen_kt4", pyJetGen_kt4, &b_pyJetGen_kt4);
   fChain->SetBranchAddress("pzJetGen_kt4", pzJetGen_kt4, &b_pzJetGen_kt4);
   fChain->SetBranchAddress("eJetGen_kt4", eJetGen_kt4, &b_eJetGen_kt4);
   fChain->SetBranchAddress("etaJetGen_kt4", etaJetGen_kt4, &b_etaJetGen_kt4);
   fChain->SetBranchAddress("phiJetGen_kt4", phiJetGen_kt4, &b_phiJetGen_kt4);
   fChain->SetBranchAddress("eTracksGen_kt4", eTracksGen_kt4, &b_eTracksGen_kt4);
   fChain->SetBranchAddress("ePhotonsGen_kt4", ePhotonsGen_kt4, &b_ePhotonsGen_kt4);
   fChain->SetBranchAddress("eMuonsGen_kt4", eMuonsGen_kt4, &b_eMuonsGen_kt4);
   fChain->SetBranchAddress("eNeutralHadronsGen_kt4", eNeutralHadronsGen_kt4, &b_eNeutralHadronsGen_kt4);


   Notify();
}

Bool_t TreeAnalyzer_RxGen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeAnalyzer_RxGen::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeAnalyzer_RxGen::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeAnalyzer_RxGen_cxx
