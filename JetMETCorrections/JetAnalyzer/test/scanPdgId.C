#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include <cmath>


int main() {

  Int_t nMC;
  Int_t pdgId[2000];
  Int_t status[2000];

  TFile* file = TFile::Open("output.root");
  TTree* tree = (TTree*) file->Get("myanalysis/pippo");

  tree->SetBranchAddress("nMC", &nMC);
  tree->SetBranchAddress("pdgIdMC", pdgId);
  tree->SetBranchAddress("statusMC", status);

  TH1F* h1_maxIndex = new TH1F("maxIndex", "", 501, -0.5, 500.5);

  int nEntries = tree->GetEntries();

  for(int jentry=0; jentry<nEntries; ++jentry) {

    tree->GetEntry(jentry);

    int foundIndex=9999;

    for(int i=0; i<nMC; ++i)
      if( fabs(pdgId[i])<10 ) foundIndex=i;

    if( foundIndex!=9999 )
      h1_maxIndex->Fill(foundIndex);

  } //for

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  h1_maxIndex->Draw();
  c1->SaveAs("maxIndex.eps");

  return 0;

}
