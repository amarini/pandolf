//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 18 22:59:13 2009 by ROOT version 5.18/00a
// from TChain pippo/
//////////////////////////////////////////////////////////

#ifndef ProvaAnalyzer_h
#define ProvaAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ProvaAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           nMC;
   Int_t           pdgIdMC[40];   //[nMC]
   Int_t           statusMC[40];   //[nMC]
   Float_t         pxMC [40];   //[nMC]
   Float_t         pyMC [40];   //[nMC]
   Float_t         pzMC [40];   //[nMC]
   Float_t         eMC  [40];   //[nMC]
   Float_t         etaMC[40];   //[nMC]
   Float_t         phiMC[40];   //[nMC]
   Int_t           nJetGen_itCone5;
   Float_t         pxJetGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         pyJetGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         pzJetGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         eJetGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         etaJetGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         phiJetGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         eTracksGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         ePhotonsGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         eMuonsGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         eNeutralHadronsGen_itCone5[33];   //[nJetGen_itCone5]
   Float_t         eHFhadronsGen_itCone5[33];   //[nJetGen_itCone5]
   Int_t           nJetGen_sisCone5;
   Float_t         pxJetGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         pyJetGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         pzJetGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         eJetGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         etaJetGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         phiJetGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         eTracksGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         ePhotonsGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         eMuonsGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         eNeutralHadronsGen_sisCone5[32];   //[nJetGen_sisCone5]
   Float_t         eHFhadronsGen_sisCone5[32];   //[nJetGen_sisCone5]
   Int_t           nJetGen_kt4;
   Float_t         pxJetGen_kt4[31];   //[nJetGen_kt4]
   Float_t         pyJetGen_kt4[31];   //[nJetGen_kt4]
   Float_t         pzJetGen_kt4[31];   //[nJetGen_kt4]
   Float_t         eJetGen_kt4[31];   //[nJetGen_kt4]
   Float_t         etaJetGen_kt4[31];   //[nJetGen_kt4]
   Float_t         phiJetGen_kt4[31];   //[nJetGen_kt4]
   Float_t         eTracksGen_kt4[31];   //[nJetGen_kt4]
   Float_t         ePhotonsGen_kt4[31];   //[nJetGen_kt4]
   Float_t         eMuonsGen_kt4[31];   //[nJetGen_kt4]
   Float_t         eNeutralHadronsGen_kt4[31];   //[nJetGen_kt4]
   Float_t         eHFhadronsGen_kt4[31];   //[nJetGen_kt4]
   Int_t           nJetCalo_itCone5;
   Float_t         pxJetCalo_itCone5[14];   //[nJetCalo_itCone5]
   Float_t         pyJetCalo_itCone5[14];   //[nJetCalo_itCone5]
   Float_t         pzJetCalo_itCone5[14];   //[nJetCalo_itCone5]
   Float_t         eJetCalo_itCone5[14];   //[nJetCalo_itCone5]
   Float_t         etaJetCalo_itCone5[14];   //[nJetCalo_itCone5]
   Float_t         phiJetCalo_itCone5[14];   //[nJetCalo_itCone5]
   Int_t           nJetCalo_sisCone5;
   Float_t         pxJetCalo_sisCone5[12];   //[nJetCalo_sisCone5]
   Float_t         pyJetCalo_sisCone5[12];   //[nJetCalo_sisCone5]
   Float_t         pzJetCalo_sisCone5[12];   //[nJetCalo_sisCone5]
   Float_t         eJetCalo_sisCone5[12];   //[nJetCalo_sisCone5]
   Float_t         etaJetCalo_sisCone5[12];   //[nJetCalo_sisCone5]
   Float_t         phiJetCalo_sisCone5[12];   //[nJetCalo_sisCone5]
   Int_t           nJetCalo_kt4;
   Float_t         pxJetCalo_kt4[12];   //[nJetCalo_kt4]
   Float_t         pyJetCalo_kt4[12];   //[nJetCalo_kt4]
   Float_t         pzJetCalo_kt4[12];   //[nJetCalo_kt4]
   Float_t         eJetCalo_kt4[12];   //[nJetCalo_kt4]
   Float_t         etaJetCalo_kt4[12];   //[nJetCalo_kt4]
   Float_t         phiJetCalo_kt4[12];   //[nJetCalo_kt4]
   Int_t           nJetPF_itCone5;
   Float_t         pxJetPF_itCone5[22];   //[nJetPF_itCone5]
   Float_t         pyJetPF_itCone5[22];   //[nJetPF_itCone5]
   Float_t         pzJetPF_itCone5[22];   //[nJetPF_itCone5]
   Float_t         eJetPF_itCone5[22];   //[nJetPF_itCone5]
   Float_t         etaJetPF_itCone5[22];   //[nJetPF_itCone5]
   Float_t         phiJetPF_itCone5[22];   //[nJetPF_itCone5]
   Int_t           firstCandIndexJetPF_itCone5[22];   //[nJetPF_itCone5]
   Int_t           lastCandIndexJetPF_itCone5[22];   //[nJetPF_itCone5]
   Int_t           nCandPF_itCone5;
   Int_t           pdgIdCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         eCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         etaCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         phiCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         pxCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         pyCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         pzCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         pxTrackCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         pyTrackCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         pzTrackCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         ecalEnergyCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         hcalEnergyCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         ps1EnergyCandPF_itCone5[228];   //[nCandPF_itCone5]
   Float_t         ps2EnergyCandPF_itCone5[228];   //[nCandPF_itCone5]

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_pdgIdMC;   //!
   TBranch        *b_statusMC;   //!
   TBranch        *b_pxMC ;   //!
   TBranch        *b_pyMC ;   //!
   TBranch        *b_pzMC ;   //!
   TBranch        *b_eMC  ;   //!
   TBranch        *b_etaMC;   //!
   TBranch        *b_phiMC;   //!
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
   TBranch        *b_eHFhadronsGen_itCone5;   //!
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
   TBranch        *b_eHFhadronsGen_sisCone5;   //!
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
   TBranch        *b_eHFhadronsGen_kt4;   //!
   TBranch        *b_nJetCalo_itCone5;   //!
   TBranch        *b_pxJetCalo_itCone5;   //!
   TBranch        *b_pyJetCalo_itCone5;   //!
   TBranch        *b_pzJetCalo_itCone5;   //!
   TBranch        *b_eJetCalo_itCone5;   //!
   TBranch        *b_etaJetCalo_itCone5;   //!
   TBranch        *b_phiJetCalo_itCone5;   //!
   TBranch        *b_nJetCalo_sisCone5;   //!
   TBranch        *b_pxJetCalo_sisCone5;   //!
   TBranch        *b_pyJetCalo_sisCone5;   //!
   TBranch        *b_pzJetCalo_sisCone5;   //!
   TBranch        *b_eJetCalo_sisCone5;   //!
   TBranch        *b_etaJetCalo_sisCone5;   //!
   TBranch        *b_phiJetCalo_sisCone5;   //!
   TBranch        *b_nJetCalo_kt4;   //!
   TBranch        *b_pxJetCalo_kt4;   //!
   TBranch        *b_pyJetCalo_kt4;   //!
   TBranch        *b_pzJetCalo_kt4;   //!
   TBranch        *b_eJetCalo_kt4;   //!
   TBranch        *b_etaJetCalo_kt4;   //!
   TBranch        *b_phiJetCalo_kt4;   //!
   TBranch        *b_nJetPF_itCone5;   //!
   TBranch        *b_pxJetPF_itCone5;   //!
   TBranch        *b_pyJetPF_itCone5;   //!
   TBranch        *b_pzJetPF_itCone5;   //!
   TBranch        *b_eJetPF_itCone5;   //!
   TBranch        *b_etaJetPF_itCone5;   //!
   TBranch        *b_phiJetPF_itCone5;   //!
   TBranch        *b_firstCandIndexJetPF_itCone5;   //!
   TBranch        *b_lastCandIndexJetPF_itCone5;   //!
   TBranch        *b_nCandPF_itCone5;   //!
   TBranch        *b_pdgIdCandPF_itCone5;   //!
   TBranch        *b_eCandPF_itCone5;   //!
   TBranch        *b_etaCandPF_itCone5;   //!
   TBranch        *b_phiCandPF_itCone5;   //!
   TBranch        *b_pxCandPF_itCone5;   //!
   TBranch        *b_pyCandPF_itCone5;   //!
   TBranch        *b_pzCandPF_itCone5;   //!
   TBranch        *b_pxTrackCandPF_itCone5;   //!
   TBranch        *b_pyTrackCandPF_itCone5;   //!
   TBranch        *b_pzTrackCandPF_itCone5;   //!
   TBranch        *b_ecalEnergyCandPF_itCone5;   //!
   TBranch        *b_hcalEnergyCandPF_itCone5;   //!
   TBranch        *b_ps1EnergyCandPF_itCone5;   //!
   TBranch        *b_ps2EnergyCandPF_itCone5;   //!

   ProvaAnalyzer(TTree *tree=0);
   virtual ~ProvaAnalyzer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ProvaAnalyzer_cxx
ProvaAnalyzer::ProvaAnalyzer(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("pippo");

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("pippo","");
      chain->Add("QCD_2_2_3/output_147.root/myanalysis/pippo");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

ProvaAnalyzer::~ProvaAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ProvaAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ProvaAnalyzer::LoadTree(Long64_t entry)
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

void ProvaAnalyzer::Init(TTree *tree)
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
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", statusMC, &b_statusMC);
   fChain->SetBranchAddress("pxMC ", pxMC , &b_pxMC );
   fChain->SetBranchAddress("pyMC ", pyMC , &b_pyMC );
   fChain->SetBranchAddress("pzMC ", pzMC , &b_pzMC );
   fChain->SetBranchAddress("eMC  ", eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", phiMC, &b_phiMC);
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
   fChain->SetBranchAddress("eHFhadronsGen_itCone5", eHFhadronsGen_itCone5, &b_eHFhadronsGen_itCone5);
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
   fChain->SetBranchAddress("eHFhadronsGen_sisCone5", eHFhadronsGen_sisCone5, &b_eHFhadronsGen_sisCone5);
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
   fChain->SetBranchAddress("eHFhadronsGen_kt4", eHFhadronsGen_kt4, &b_eHFhadronsGen_kt4);
   fChain->SetBranchAddress("nJetCalo_itCone5", &nJetCalo_itCone5, &b_nJetCalo_itCone5);
   fChain->SetBranchAddress("pxJetCalo_itCone5", pxJetCalo_itCone5, &b_pxJetCalo_itCone5);
   fChain->SetBranchAddress("pyJetCalo_itCone5", pyJetCalo_itCone5, &b_pyJetCalo_itCone5);
   fChain->SetBranchAddress("pzJetCalo_itCone5", pzJetCalo_itCone5, &b_pzJetCalo_itCone5);
   fChain->SetBranchAddress("eJetCalo_itCone5", eJetCalo_itCone5, &b_eJetCalo_itCone5);
   fChain->SetBranchAddress("etaJetCalo_itCone5", etaJetCalo_itCone5, &b_etaJetCalo_itCone5);
   fChain->SetBranchAddress("phiJetCalo_itCone5", phiJetCalo_itCone5, &b_phiJetCalo_itCone5);
   fChain->SetBranchAddress("nJetCalo_sisCone5", &nJetCalo_sisCone5, &b_nJetCalo_sisCone5);
   fChain->SetBranchAddress("pxJetCalo_sisCone5", pxJetCalo_sisCone5, &b_pxJetCalo_sisCone5);
   fChain->SetBranchAddress("pyJetCalo_sisCone5", pyJetCalo_sisCone5, &b_pyJetCalo_sisCone5);
   fChain->SetBranchAddress("pzJetCalo_sisCone5", pzJetCalo_sisCone5, &b_pzJetCalo_sisCone5);
   fChain->SetBranchAddress("eJetCalo_sisCone5", eJetCalo_sisCone5, &b_eJetCalo_sisCone5);
   fChain->SetBranchAddress("etaJetCalo_sisCone5", etaJetCalo_sisCone5, &b_etaJetCalo_sisCone5);
   fChain->SetBranchAddress("phiJetCalo_sisCone5", phiJetCalo_sisCone5, &b_phiJetCalo_sisCone5);
   fChain->SetBranchAddress("nJetCalo_kt4", &nJetCalo_kt4, &b_nJetCalo_kt4);
   fChain->SetBranchAddress("pxJetCalo_kt4", pxJetCalo_kt4, &b_pxJetCalo_kt4);
   fChain->SetBranchAddress("pyJetCalo_kt4", pyJetCalo_kt4, &b_pyJetCalo_kt4);
   fChain->SetBranchAddress("pzJetCalo_kt4", pzJetCalo_kt4, &b_pzJetCalo_kt4);
   fChain->SetBranchAddress("eJetCalo_kt4", eJetCalo_kt4, &b_eJetCalo_kt4);
   fChain->SetBranchAddress("etaJetCalo_kt4", etaJetCalo_kt4, &b_etaJetCalo_kt4);
   fChain->SetBranchAddress("phiJetCalo_kt4", phiJetCalo_kt4, &b_phiJetCalo_kt4);
   fChain->SetBranchAddress("nJetPF_itCone5", &nJetPF_itCone5, &b_nJetPF_itCone5);
   fChain->SetBranchAddress("pxJetPF_itCone5", pxJetPF_itCone5, &b_pxJetPF_itCone5);
   fChain->SetBranchAddress("pyJetPF_itCone5", pyJetPF_itCone5, &b_pyJetPF_itCone5);
   fChain->SetBranchAddress("pzJetPF_itCone5", pzJetPF_itCone5, &b_pzJetPF_itCone5);
   fChain->SetBranchAddress("eJetPF_itCone5", eJetPF_itCone5, &b_eJetPF_itCone5);
   fChain->SetBranchAddress("etaJetPF_itCone5", etaJetPF_itCone5, &b_etaJetPF_itCone5);
   fChain->SetBranchAddress("phiJetPF_itCone5", phiJetPF_itCone5, &b_phiJetPF_itCone5);
   fChain->SetBranchAddress("firstCandIndexJetPF_itCone5", firstCandIndexJetPF_itCone5, &b_firstCandIndexJetPF_itCone5);
   fChain->SetBranchAddress("lastCandIndexJetPF_itCone5", lastCandIndexJetPF_itCone5, &b_lastCandIndexJetPF_itCone5);
   fChain->SetBranchAddress("nCandPF_itCone5", &nCandPF_itCone5, &b_nCandPF_itCone5);
   fChain->SetBranchAddress("pdgIdCandPF_itCone5", pdgIdCandPF_itCone5, &b_pdgIdCandPF_itCone5);
   fChain->SetBranchAddress("eCandPF_itCone5", eCandPF_itCone5, &b_eCandPF_itCone5);
   fChain->SetBranchAddress("etaCandPF_itCone5", etaCandPF_itCone5, &b_etaCandPF_itCone5);
   fChain->SetBranchAddress("phiCandPF_itCone5", phiCandPF_itCone5, &b_phiCandPF_itCone5);
   fChain->SetBranchAddress("pxCandPF_itCone5", pxCandPF_itCone5, &b_pxCandPF_itCone5);
   fChain->SetBranchAddress("pyCandPF_itCone5", pyCandPF_itCone5, &b_pyCandPF_itCone5);
   fChain->SetBranchAddress("pzCandPF_itCone5", pzCandPF_itCone5, &b_pzCandPF_itCone5);
   fChain->SetBranchAddress("pxTrackCandPF_itCone5", pxTrackCandPF_itCone5, &b_pxTrackCandPF_itCone5);
   fChain->SetBranchAddress("pyTrackCandPF_itCone5", pyTrackCandPF_itCone5, &b_pyTrackCandPF_itCone5);
   fChain->SetBranchAddress("pzTrackCandPF_itCone5", pzTrackCandPF_itCone5, &b_pzTrackCandPF_itCone5);
   fChain->SetBranchAddress("ecalEnergyCandPF_itCone5", ecalEnergyCandPF_itCone5, &b_ecalEnergyCandPF_itCone5);
   fChain->SetBranchAddress("hcalEnergyCandPF_itCone5", hcalEnergyCandPF_itCone5, &b_hcalEnergyCandPF_itCone5);
   fChain->SetBranchAddress("ps1EnergyCandPF_itCone5", ps1EnergyCandPF_itCone5, &b_ps1EnergyCandPF_itCone5);
   fChain->SetBranchAddress("ps2EnergyCandPF_itCone5", ps2EnergyCandPF_itCone5, &b_ps2EnergyCandPF_itCone5);
   Notify();
}

Bool_t ProvaAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ProvaAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ProvaAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ProvaAnalyzer_cxx
