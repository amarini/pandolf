#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "DrawBase.h"



float computeChiSquare( TH1F* h1_DATA, TH1F* h1_MC );


int main() {

  // just to set the style:
  DrawBase* db = new DrawBase("TTSF");
  
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(4980.);

  TFile* file_DATA = TFile::Open("TTZTrilepton_DATA_Run2011_FULL_presel_TCHE_ALL.root");

  TFile* file_DYJets = TFile::Open("TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11_presel_TCHE_ALL.root");
  TFile* file_TTJets = TFile::Open("TTZTrilepton_TTJ_Fall11_highstat_presel_TCHE_ALL.root");
  TFile* file_VV     = TFile::Open("TTZTrilepton_VV_Summer11_presel_TCHE_ALL.root");


  TH1F* h1_DATA = (TH1F*)file_DATA->Get("mZll_OF_prepresel");

  TH1F* h1_DYJets = (TH1F*)file_DYJets->Get("mZll_OF_prepresel");
  TH1F* h1_TTJets = (TH1F*)file_TTJets->Get("mZll_OF_prepresel");
  TH1F* h1_VV     = (TH1F*)file_VV->Get("mZll_OF_prepresel");


  h1_DATA  ->Rebin(4);
  h1_DYJets->Rebin(4);
  h1_TTJets->Rebin(4);
  h1_VV    ->Rebin(4);

  h1_DYJets->Scale(4980.);
  h1_TTJets->Scale(4980.);
  h1_VV    ->Scale(4980.);

  std::vector< TH1F* > otherHistos;
  otherHistos.push_back( h1_DYJets );
  otherHistos.push_back( h1_VV );

  float minSF = 0.9;
  float maxSF = 1.5;

  int nSteps = 1000;
  float step = (maxSF-minSF)/(float)nSteps;

  TH1F* h1_chiSquare = new TH1F("chiSquare", "", nSteps, minSF, maxSF );


  float minChiSquare = 9999.;
  float foundSF = -1.;

  for( unsigned i=0; i<nSteps; ++i ) {

    float thisSF = minSF + i*step;

    TH1F* h1_TTJets_sf = new TH1F(*h1_TTJets);
    h1_TTJets_sf->Scale( thisSF );

    for( unsigned ihisto=0; ihisto<otherHistos.size(); ++ihisto )
      h1_TTJets_sf->Add( otherHistos[ihisto] );

    float thisChiSquare = computeChiSquare( h1_DATA, h1_TTJets_sf );

    h1_chiSquare->SetBinContent( i+1, thisChiSquare ); 

    if( thisChiSquare < minChiSquare ) {
      minChiSquare = thisChiSquare;
      foundSF = thisSF;
    }

  }

  h1_chiSquare->SetXTitle( "t#bar{t} Scale Factor" );
  h1_chiSquare->SetYTitle( "#chi^{2}" );
  h1_chiSquare->SetMarkerStyle( 20 );
  h1_chiSquare->SetMarkerSize( 1.6 );
  h1_chiSquare->SetMarkerColor( 46 );
  h1_chiSquare->GetYaxis()->SetRangeUser(0., h1_chiSquare->GetMaximum());

  std::cout << "Min Chi Square: " << minChiSquare << std::endl;
  std::cout << "Best SF: " << foundSF << std::endl;

  TPaveText* label_sqrt = db->get_labelSqrt();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  h1_chiSquare->Draw( "P" );
  label_sqrt->Draw( "P" );

  c1->SaveAs("chiSquare.eps");

  return 0;

}




float computeChiSquare( TH1F* h1_DATA, TH1F* h1_MC ) {

  float chiSquare=0.;
  int nbins = h1_DATA->GetNbinsX();

  for( unsigned ibin=1; ibin<nbins+1; ++ibin ) {

    float data = h1_DATA->GetBinContent(ibin);
    float mc = h1_MC->GetBinContent(ibin);

    float binDifference = fabs(data - mc);
    float err = sqrt(data); // fairly high stat in every bin anyways

    float addendum = binDifference / err;

    chiSquare += (addendum*addendum);
   
  } //for bins

  chiSquare /= (nbins-1);

  return chiSquare;

}

