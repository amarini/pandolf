#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>

#include "DrawBase.h"
#include "CommonTools/StatTools.h"





int main( int argc, char* argv[] ) {

  std::string selectionType = "May7_ht_noMET_2loosebjets";
  if( argc>1 ) {
    std::string selectionType_str(argv[1]);
    selectionType = selectionType_str;
  }


  // this sets the style:
  DrawBase* db = new DrawBase("OPT_ZBi");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(4980.);

  TPaveText* label_sqrt = db->get_labelSqrt();
    

  std::string optcutsdir = "optcuts_" + selectionType;
  std::string ZBiFileName = optcutsdir + "/ZBiScan.txt";

  ofstream ofs_ZBi(ZBiFileName.c_str());
  ofs_ZBi << "Expected for 5 fb-1:" << std::endl;
  ofs_ZBi << "Seff   \tS     \tB +- s(B)\tZBi" << std::endl;

  TGraphErrors* gr_ZBi = new TGraphErrors(0);
  float ZBi_max = 0.;
  float effS_ZBi_max = 0.;
  float effMax = 0.;

  TFile* signalFile = TFile::Open("../TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi_presel_TCHE_ALL.root");
  TTree* signalTree = (TTree*)signalFile->Get("tree_passedEvents");

  TH1D* h1_nCounter_signal = (TH1D*)signalFile->Get("nCounter");
  float nGen_signal = h1_nCounter_signal->GetBinContent(1);

  TChain* backgroundTree = new TChain("tree_passedEvents");
  backgroundTree->Add("../TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11_presel_TCHE_ALL.root/tree_passedEvents");
  backgroundTree->Add("../TTZTrilepton_TTJ_Fall11_highstat_presel_TCHE_ALL.root/tree_passedEvents");
  backgroundTree->Add("../TTZTrilepton_VV_Summer11_presel_TCHE_ALL.root/tree_passedEvents");
  backgroundTree->Add("../TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi_presel_TCHE_ALL.root/tree_passedEvents");




  for( unsigned iEff=1; iEff<10; ++iEff ) {

    char infileName[300];
    sprintf( infileName, "%s/cuts_Seff%d.txt", optcutsdir.c_str(), iEff*10);
    ifstream ifs(infileName);
    std::cout << "-> Opening Seff file: " << infileName << std::endl;
  
    std::vector<std::string> varNames;
    std::vector<float> cutsMin;
    std::vector<float> cutsMax;

    while( ifs.good() && !ifs.eof() ) {

      std::string varName;
      float cutMin, cutMax;

      ifs >> varName >> cutMin >> cutMax;

      varNames.push_back( varName );
      cutsMin.push_back( cutMin );
      cutsMax.push_back( cutMax );

    } //while file is good
  
    ifs.close();

    // eliminate last element (last line is read and is empty):
    varNames.pop_back();
    cutsMin.pop_back();
    cutsMax.pop_back();


    std::string selection = "eventWeight*( ";
    for( unsigned ivar=0; ivar<varNames.size(); ++ivar ) {
      if( ivar!=0 ) selection += " && ";
      char thisCut[200];
      sprintf( thisCut, "%s >= %f && %s < %f", varNames[ivar].c_str(), cutsMin[ivar], varNames[ivar].c_str(), cutsMax[ivar] );
      std::string thisCut_str(thisCut);
      selection += thisCut_str;
    }

    selection += " && leptType<2 )";

//std::cout << selection << std::endl;

    TH1F* h1_bg = new TH1F("bg", "", 2, 0, 2);
    h1_bg->Sumw2();
    TH1F* h1_signal = new TH1F("signal", "", 2, 0, 2);
    h1_signal->Sumw2();
   
    signalTree->Project( "signal", "leptType", selection.c_str() );
    backgroundTree->Project( "bg", "leptType", selection.c_str() );

    double lumi = 5000.;
    
    double signal = h1_signal->Integral();
    double background_error;
    double background = h1_bg->IntegralAndError( 0, 3, background_error );
   
    signal *= lumi;
    background *= lumi;
    background_error *= lumi;

    float ZBi = StatTools::computeZBi( signal+background, background, background_error );

std::cout << "signal: " << signal << " bg: " << background << " +- " << background_error << std::endl;


    float effS = (float)h1_signal->GetEntries()/nGen_signal/(0.22*0.22*0.67);
    //float effS = (float)h1_signal->GetEntries()/nGen_signal/(0.22*0.22*0.67);

    if( effS > effMax )
      effMax = effS;

    gr_ZBi->SetPoint( iEff-1, 100.*effS, ZBi );

    if( ZBi > ZBi_max ) {
      ZBi_max = ZBi;
      effS_ZBi_max = effS;
    }

    float yMax = h1_signal->GetMaximum() + h1_bg->GetMaximum();
    yMax*=1.5;

//  THStack* stack = new THStack();
//  stack->Add( h1_bg );
//  stack->Add( h1_signal );

//  TH2D* h2_axes = new TH2D("axes", "", 3, -0.5, 2.5, 10, 0., yMax);
//  h2_axes->GetXaxis()->SetLabelSize(0.085);
//  h2_axes->GetXaxis()->SetBinLabel(1, "#mu#mu");
//  h2_axes->GetXaxis()->SetBinLabel(2, "e#mu");
//  h2_axes->GetXaxis()->SetBinLabel(3, "ee");
//  h2_axes->SetYTitle("Events");


//  TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.88);
//  legend->SetFillColor(0);
//  legend->SetTextSize(0.035);
//  legend->AddEntry( h1_signal, "Signal", "F");
//  legend->AddEntry( h1_bg, "Background", "F");

//  char canvasName[250];
//  sprintf( canvasName, "%s/yieldPlot_Seff%d.eps", optcutsdir.c_str(), iEff*10);

    //TPaveText* label = new TPaveText( 0.15, 0.65, 0.45, 0.85, "brNDC");
    //label->SetFillColor(0);
    //label->SetTextSize(0.035);
    //label->AddText("L = 5 fb^{-1}");
    //char signalLabel[100];
    //sprintf( signalLabel, "s = %.2f (%d%%)", s, (int)(((float)h1_signal->GetEntries()/nTotal_s)*100) );
    //label->AddText( signalLabel );
    //char bgLabel[100];
    //sprintf( bgLabel, "b = %.2f", b_pred);
    //label->AddText( bgLabel );
    //char signifLabel[100];
    //sprintf( signifLabel, "ZBi = %.2f", ZBi);
    //label->AddText( signifLabel );

//  char b_text[500];
//  sprintf( b_text, "BG = %.2f #pm %.2f", b_pred, b_pred_err );
//  char obs_text[500];
//  sprintf( obs_text, "OBS = %.2f", obs );
//  char Zbi_text[100];
//  sprintf( Zbi_text, "ZBi = %.3f", ZBi );
//  TPaveText* label_ZBi = new TPaveText( 0.23, 0.73, 0.5, 0.88, "brNDC" );
//  label_ZBi->SetFillColor(0);
//  label_ZBi->SetTextSize(0.035);
//  label_ZBi->AddText(b_text);
//  label_ZBi->AddText(obs_text);
//  label_ZBi->AddText(Zbi_text);

//  TCanvas* c1 = new TCanvas("c1", "c1", 600., 600.);
//  c1->cd();
//  h2_axes->Draw();
//  //h1_signal->Draw("same");
//  stack->Draw("histo same");
//  h1_bg->Draw("0 E2 same");
//  legend->Draw("same");
//  //label->Draw("same");
//  label_ZBi->Draw("same");
//  label_sqrt->Draw("same");
//  gPad->RedrawAxis();
//  c1->SaveAs(canvasName);

//  delete c1;
//  delete legend;
//  delete h2_axes;
    //delete stack;
    

    ofs_ZBi << effS << "\t" << signal << "\t" << background << " +- " << background_error << "\t" << ZBi << std::endl;

    delete h1_signal;
    delete h1_bg;

std::cout << "### " << iEff << "   ZBi: " << ZBi << std::endl;
  } // for iEff

  std::cout << "> > >   BEST ZBi: " << ZBi_max << std::endl;
  std::cout << "> > >   signal eff: " << effS_ZBi_max << std::endl;

  ofs_ZBi.close();

  db->resetStyle();

  gr_ZBi->SetMarkerSize(2.);
  gr_ZBi->SetMarkerStyle(21);
  gr_ZBi->SetMarkerColor(kRed+3);


  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.3*effMax*100., 10, 0., 1.6*ZBi_max ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("ZBi (5 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->cd();

  
  h2_axes_gr->Draw();
  gr_ZBi->Draw("P same");
  label_sqrt->Draw("same");

  char ZBi_vs_Seff_name[250];
  sprintf(ZBi_vs_Seff_name, "%s/ZBi_vs_Seff.eps", optcutsdir.c_str() );
  c_gr->SaveAs(ZBi_vs_Seff_name);

}




