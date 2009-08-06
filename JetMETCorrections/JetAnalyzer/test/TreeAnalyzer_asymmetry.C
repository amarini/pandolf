#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include <vector>
#include <string>





//void asymmetry(std::string algoType="caloItCone5") {
void asymmetry() {

  int nBins = 5;


  TChain* chain = new TChain("jetTree");
  int nFiles=49;

  for(int i=1; i<(nFiles+1); ++i) {
    char name[100];
    sprintf(name, "2ndLevelOutputFiles/2ndLevelOutputFile_%d.root/jetTree", i);
    chain->Add(name);
  }


 
  int nJetReco_caloItCone5;
  chain->SetBranchAddress("nJetReco_caloItCone5", &nJetReco_caloItCone5);
 
  Float_t ptReco_caloItCone5[40];
  chain->SetBranchAddress("ptReco_caloItCone5", ptReco_caloItCone5);
 
  Float_t phiReco_caloItCone5[40];
  chain->SetBranchAddress("phiReco_caloItCone5", phiReco_caloItCone5);
 

  std::vector< TH1F* > h1_asymmetry_less10;
  std::vector< TH1F* > h1_asymmetry_less12;
  std::vector< TH1F* > h1_asymmetry_less15;
  std::vector< TH1F* > h1_asymmetry_less17;
  std::vector< TH1F* > h1_asymmetry_less20;
  std::vector< TH1F* > h1_asymmetry_less30;


  for(int i=0; i<nBins; ++i) {

    char histName[40];
    sprintf(histName, "less10_bin%d", i);
    TH1F* less10 = new TH1F(histName, "", 100, -0.6, 0.6);
    h1_asymmetry_less10.push_back(less10);

    sprintf(histName, "less12_bin%d", i);
    TH1F* less12 = new TH1F(histName, "", 100, -0.6, 0.6);
    h1_asymmetry_less12.push_back(less12);

    sprintf(histName, "less15_bin%d", i);
    TH1F* less15 = new TH1F(histName, "", 100, -0.6, 0.6);
    h1_asymmetry_less15.push_back(less15);

    sprintf(histName, "less17_bin%d", i);
    TH1F* less17 = new TH1F(histName, "", 100, -0.6, 0.6);
    h1_asymmetry_less17.push_back(less17);

    sprintf(histName, "less20_bin%d", i);
    TH1F* less20 = new TH1F(histName, "", 100, -0.6, 0.6);
    h1_asymmetry_less20.push_back(less20);

    sprintf(histName, "less30_bin%d", i);
    TH1F* less30 = new TH1F(histName, "", 100, -0.6, 0.6);
    h1_asymmetry_less30.push_back(less30);

  }


  TRandom t;

  Int_t nEntries = chain->GetEntries();

  //for(int iEntry=0; iEntry<100000; ++iEntry) { 
  for(int iEntry=0; iEntry<nEntries; ++iEntry) { 

    if( (iEntry % 100000) == 0 ) std::cout << "Entry #" << iEntry << " /" << nEntries << std::endl;
    chain->GetEntry(iEntry);

    //------------ CALO IT CONE 0.5

    if( nJetReco_caloItCone5 > 2 ) {

      Float_t pt1(0.), pt2(0.), pt3(0.);
      Float_t phi1(0.), phi2(0.);


      for(int iJet = 0; iJet<nJetReco_caloItCone5; ++iJet ) {

        Float_t ptJet = ptReco_caloItCone5[iJet];

        if( ptJet > pt1 ) {
          pt3=pt2;
          pt2=pt1;
          pt1=ptJet;
          phi1=phiReco_caloItCone5[iJet];
        } else if( ptJet > pt2 ) {
          pt3=pt2;
          pt2=ptJet;
          phi2=phiReco_caloItCone5[iJet];
        } else if( ptJet > pt3 ) {
          pt3=ptJet;
        }
        
      } //for jets


      Float_t pi = 3.14159;
      //Float_t deltaPhi = acos(cos(phi1-phi2));
      Float_t deltaPhi = phi1-phi2;
      if( deltaPhi > pi ) deltaPhi-=2*pi;
      if( deltaPhi < -pi ) deltaPhi+=2*pi;

      if( (pt3 < 30.) && (deltaPhi>2.7) ) {

        Float_t asymmetry = (pt1-pt2)/(pt1+pt2);
        Float_t oneHalfProb = t.Uniform(1.);
        if(oneHalfProb<0.5) asymmetry *= -1.;

        Float_t averagePt = (pt1+pt2)/2.;
        
        int iBin=-1;
         
        if( (averagePt > 80.)&&(averagePt < 100.) ) iBin=0;
        else if( (averagePt > 100.)&&(averagePt < 125.) ) iBin=1;
        else if( (averagePt > 125.)&&(averagePt < 160.) ) iBin=2;
        else if( (averagePt > 160.)&&(averagePt < 190.) ) iBin=3;
        else if( (averagePt > 190.)&&(averagePt < 215.) ) iBin=4;

        if(iBin!=-1) {

          h1_asymmetry_less30[iBin]->Fill(asymmetry);

          if( pt3 < 20. ) h1_asymmetry_less20[iBin]->Fill(asymmetry);
          if( pt3 < 17. ) h1_asymmetry_less17[iBin]->Fill(asymmetry);
          if( pt3 < 15. ) h1_asymmetry_less15[iBin]->Fill(asymmetry);
          if( pt3 < 12. ) h1_asymmetry_less12[iBin]->Fill(asymmetry);
          if( pt3 < 10. ) h1_asymmetry_less10[iBin]->Fill(asymmetry);

        } // if ibin

      } //if pt3 < 30

    } //if nJet>2

  } //for entries

  TCanvas* c_prova = new TCanvas("c_prova", "c_prova", 800, 600);
  c_prova->cd();
  


  Double_t Lower[nBins];

  Lower[0] = 80.;  
  Lower[1] = 100.;  
  Lower[2] = 125.;  
  Lower[3] = 160.;  
  Lower[4] = 190.;  
  Lower[5] = 215.;  


  TH1F* h1_resolution_vs_pt = new TH1F("resolution_vs_pt", "Jet Energy Resolution vs. Average Jet p_{T}", nBins, Lower);
  h1_resolution_vs_pt->SetXTitle("Average p_{T}^{RECO} [GeV/c]");
  h1_resolution_vs_pt->SetYTitle("#sigma(p_{T})/p_{T}");


  for(int iBin=0; iBin<nBins; ++iBin) {

    int nPoints = 6;
    Float_t x[nPoints], y[nPoints];
    Float_t err_x[nPoints], err_y[nPoints];

    x[0] = 10.;
    x[1] = 12.;
    x[2] = 15.;
    x[3] = 17.;
    x[4] = 20.;
    x[5] = 30.;

    y[0] = sqrt(2)*h1_asymmetry_less10[iBin]->GetRMS();
    y[1] = sqrt(2)*h1_asymmetry_less12[iBin]->GetRMS();
    y[2] = sqrt(2)*h1_asymmetry_less15[iBin]->GetRMS();
    y[3] = sqrt(2)*h1_asymmetry_less17[iBin]->GetRMS();
    y[4] = sqrt(2)*h1_asymmetry_less20[iBin]->GetRMS();
    y[5] = sqrt(2)*h1_asymmetry_less30[iBin]->GetRMS();

    err_x[0] = 0.; 
    err_x[1] = 0.; 
    err_x[2] = 0.; 
    err_x[3] = 0.; 
    err_x[4] = 0.; 
    err_x[5] = 0.; 

    err_y[0] = y[0]/sqrt(h1_asymmetry_less10[iBin]->GetEntries()); 
    err_y[1] = y[1]/sqrt(h1_asymmetry_less12[iBin]->GetEntries()); 
    err_y[2] = y[2]/sqrt(h1_asymmetry_less15[iBin]->GetEntries()); 
    err_y[3] = y[3]/sqrt(h1_asymmetry_less17[iBin]->GetEntries()); 
    err_y[4] = y[4]/sqrt(h1_asymmetry_less20[iBin]->GetEntries()); 
    err_y[5] = y[5]/sqrt(h1_asymmetry_less30[iBin]->GetEntries()); 


    TGraphErrors* resolution_vs_pt3 = new TGraphErrors(nPoints, x, y, err_x, err_y);
    resolution_vs_pt3->SetMarkerStyle(20);

    TF1* line = new TF1("line", "[0] + [1]*x");
    line->SetParameter(0, y[0]);
    line->SetParameter(1, (y[1]-y[0])/(x[1]-x[0]));
    line->SetRange(0., 32.);

    resolution_vs_pt3->Fit(line, "R");

    char plotName[70];
    sprintf(plotName, "asymmetry/asymmetry_less10_%d.eps", iBin);
    h1_asymmetry_less10[iBin]->Draw();
    c_prova->SaveAs(plotName);

    sprintf(plotName, "asymmetry/asymmetry_less20_%d.eps", iBin);
    h1_asymmetry_less20[iBin]->Draw();
    c_prova->SaveAs(plotName);

    sprintf(plotName, "asymmetry/asymmetry_less30_%d.eps", iBin);
    h1_asymmetry_less30[iBin]->Draw();
    c_prova->SaveAs(plotName);


    sprintf(plotName, "asymmetry/resolution_vs_pt3_%d.eps", iBin);
    resolution_vs_pt3->Draw("AP");
    c_prova->SaveAs(plotName);


    h1_resolution_vs_pt->SetBinContent(iBin+1, line->GetParameter(0));
    h1_resolution_vs_pt->SetBinError(iBin+1, line->GetParError(0));

  } //for pt bins


  TFile* outFile = new TFile("asymmetryOutputFile.root", "RECREATE");

  outFile->cd();
  h1_resolution_vs_pt->Write();
  outFile->Close();
  outFile->Write();

} //asymmetry
     
  
  

