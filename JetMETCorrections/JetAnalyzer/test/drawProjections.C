#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"


Int_t NBINS_ = 20;


void draw(const std::string& name);

void drawProjections() {

  TStyle *simpleStyle = new TStyle("simpleStyle","");
  simpleStyle->SetCanvasColor(0);
  simpleStyle->SetFrameFillColor(0);
  simpleStyle->SetStatColor(0);
  simpleStyle->SetOptStat(0);
  simpleStyle->SetTitleFillColor(0);
  simpleStyle->SetCanvasBorderMode(0);
  simpleStyle->SetPadBorderMode(0);
  simpleStyle->SetFrameBorderMode(0);
  simpleStyle->cd();


  draw("PFItCone5_barrel");
  draw("PFItCone5_endcap");
  draw("caloItCone5_barrel");
  draw("caloItCone5_endcap");
  draw("Rch90_100_endcap");


}



void draw(const std::string& name) {

  std::string fileName = "Projections/" + name + ".root";
  TFile* file = TFile::Open(fileName.c_str());

  for(int i = 1; i<NBINS_; ++i) {

    char histName[50];
    sprintf(histName, "projection_%d", i);
    TH1F* h1 = (TH1F*)file->Get(histName);

    h1->SetTitle("");
    h1->SetLineWidth(1.5);

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();
    std::string canvasName(histName);
    canvasName = "Projections/"+ name + "_" + canvasName + ".eps";
    h1->Draw();
    c1->SaveAs(canvasName.c_str());
    delete c1;

  } //for

} 
