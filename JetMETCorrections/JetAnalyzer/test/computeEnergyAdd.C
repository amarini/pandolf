#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"


void computeEnergyAdd() {

  TFile* file = TFile::Open("Rch0.root");
  
  TH1F* h1_PF = (TH1F*)file->Get("responseE_vs_E_MEAN_barrel_PF");
  TH1F* h1_calo = (TH1F*)file->Get("responseE_vs_E_MEAN_barrel_calo");

  int nBins = 21;
  Double_t Lower[nBins];

  for(int iLower=0; iLower<11; ++iLower)
    Lower[iLower] = iLower*10.;

  Lower[11] = 120.;
  Lower[12] = 140.;
  Lower[13] = 160.;
  Lower[14] = 180.;
  Lower[15] = 200.;
  Lower[16] = 250.;
  Lower[17] = 300.;
  Lower[18] = 400.;
  Lower[19] = 500.;
  Lower[20] = 600.;

  TH1F* h1_deltaEreco = new TH1F("deltaEreco", "", nBins-2, Lower);

  for(int i=1; i<nBins; ++i) {
    Float_t genEnergy = (Lower[i+1]-Lower[i])/2.;
    Float_t PFresponse = h1_PF->GetBinContent(i);
    Float_t caloResponse = h1_calo->GetBinContent(i);
    Float_t deltaEreco = (PFresponse-caloResponse)*genEnergy; 
    h1_deltaEreco->SetBinContent(i, deltaEreco);
  }



  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  h1_deltaEreco->Draw();
  c1->SaveAs("energyAdd.eps");

}
