#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TGraphErrors.h"

#include "DrawBase.h"
#include "fitTools.h"


std::string FIT_RMS = "RMS";



void draw_vs_pt_plots( DrawBase* db, bool correctedPt=false );


int main( int argc, char* argv[] ) {

  if( argc != 6 && argc!=7 ) {
    std::cout << "USAGE: ./draw_all_methods_vs_pt [DATA_dataset] [mcSignal_dataset] [mcBG_dataset] [recoType] [jetAlgo] [flags=\"\"]" << std::endl;
    exit(53);
  }

  std::string data_dataset( argv[1] );
  std::string mc_photonjet( argv[2] );
  std::string mc_QCD( argv[3] );
  std::string recoType( argv[4] );
  std::string jetAlgo( argv[5] );
  std::string flags = "";
  if( argc==7 ) {
    std::string flags_str(argv[6]);
    flags = flags_str;
  }

  std::string algoType;
  if( recoType=="calo" )
    algoType = jetAlgo;
  else
    algoType = recoType+jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5") algoType = "jptak5";
  if( recoType=="jpt" && jetAlgo=="akt7") algoType = "jptak7";
  

  DrawBase* db = new DrawBase("PhotonJet", recoType, jetAlgo);
  db->set_flags( flags );

  char dataFileName[200];
  if( flags != "" ) {
    sprintf( dataFileName, "PhotonJet_%s_%s_%s.root", data_dataset.c_str(), algoType.c_str(), flags.c_str());
  } else {
    sprintf( dataFileName, "PhotonJet_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  }
  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->add_dataFile( dataFile, data_dataset );


  char mcPhotonJetFileName[200];
  std::string signalflags(flags);
  //if( useMCassoc_ ) signalflags = flags;
  if( signalflags!="" && flags!="GENJETS" ) {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s_%s.root", mc_photonjet.c_str(), algoType.c_str(), signalflags.c_str());
  } else {
    sprintf( mcPhotonJetFileName, "PhotonJet_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName);
  std::cout << "Opened mc file '" << mcPhotonJetFileName << "'." << std::endl;

  db->add_mcFile( mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);


  char mcQCDFile_char[200];
  std::string bgflags(flags);
  //if( useMCassoc_ && mc_QCD==mc_photonjet ) bgflags = signalflags;
  if( bgflags!="" ) {
    sprintf( mcQCDFile_char, "PhotonJet_%s_%s_%s.root", mc_QCD.c_str(), algoType.c_str(), bgflags.c_str());
  } else {
    sprintf( mcQCDFile_char, "PhotonJet_%s_%s.root", mc_QCD.c_str(), algoType.c_str());
  }
  TFile* mcQCDFile = TFile::Open(mcQCDFile_char);
  std::cout << "Opened mc file '" << mcQCDFile_char << "'." << std::endl;

  if( mc_QCD!=mc_photonjet )
    db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);




  db->set_shapeNormalization();

  db->set_outputdir();



  draw_vs_pt_plots(db);
  draw_vs_pt_plots(db, (bool)true);

  return 0;

}





void draw_vs_pt_plots( DrawBase* db, bool correctedPt ) {

  std::string file_noextrap_name = "PhotonJetGraphs_" + db->get_fullSuffix() + ".root";
  std::string file_extrap_name = "PhotonJetExtrapGraphs_" + db->get_fullSuffix() + "_" + FIT_RMS + ".root";

  //now open graph files and plot them on on top of the other:
  TFile* file_noextrap = TFile::Open( file_noextrap_name.c_str(), "read");
  if( file_noextrap!=0 ) std::cout << "-> Opened file: '" << file_noextrap_name << "'" << std::endl;
  TFile* file_extrap = TFile::Open( file_extrap_name.c_str(), "read");
  if( file_extrap!=0 ) std::cout << "-> Opened file: '" << file_extrap_name << "'" << std::endl;


  std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();
  float xMin = ptPhot_binning[0];
  float xMax = ptPhot_binning[ ptPhot_binning.size()-2 ];

  int markerSize = 2.;


  std::string responseGEN_name = (correctedPt) ? "GENresponse_L2L3_vs_pt" : "GENresponse_vs_pt";
  TGraphErrors* gr_responseGEN_vs_pt = (TGraphErrors*)file_noextrap->Get(responseGEN_name.c_str());
  gr_responseGEN_vs_pt->SetMarkerStyle(29);
  gr_responseGEN_vs_pt->SetMarkerSize(markerSize);
  gr_responseGEN_vs_pt->SetMarkerColor(kBlack);

  std::string funcType = (db->get_recoType()=="calo") ? "powerlaw" : "rpf";

  TF1* fit_responseGEN = fitTools::fitResponseGraph(gr_responseGEN_vs_pt, funcType, "responseGEN_vs_pt", "R");
  fit_responseGEN->SetLineWidth(2.);
  fit_responseGEN->SetRange(xMin, xMax);
  // now get fit error band:
  TH1D* band_responseGEN = fitTools::getBand(fit_responseGEN, "band_responseGEN");
  band_responseGEN->SetFillColor(kYellow-9);
  band_responseGEN->SetLineWidth(2.);


  std::string responseBALANCING_name = (correctedPt) ? "response_L2L3_vs_pt" : "response_vs_pt";
  TGraphErrors* gr_responseBALANCING_vs_pt = (TGraphErrors*)file_noextrap->Get(responseBALANCING_name.c_str());
  gr_responseBALANCING_vs_pt->SetMarkerStyle(21);
  gr_responseBALANCING_vs_pt->SetMarkerSize(markerSize);
  gr_responseBALANCING_vs_pt->SetMarkerColor(kGray+2);

  std::string responseBALANCINGMC_name = (correctedPt) ? "MCresponse_L2L3_vs_pt" : "MCresponse_vs_pt";
  TGraphErrors* gr_responseBALANCINGMC_vs_pt = (TGraphErrors*)file_noextrap->Get(responseBALANCINGMC_name.c_str());
  gr_responseBALANCINGMC_vs_pt->SetMarkerStyle(25);
  gr_responseBALANCINGMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseBALANCINGMC_vs_pt->SetMarkerColor(kGray+2);


  TGraphErrors* gr_responseMPF_vs_pt = (TGraphErrors*)file_noextrap->Get("responseMPF_vs_pt");
  gr_responseMPF_vs_pt->SetMarkerStyle(20);
  gr_responseMPF_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPF_vs_pt->SetMarkerColor(38);

  TGraphErrors* gr_responseMPFMC_vs_pt = (TGraphErrors*)file_noextrap->Get("MCresponseMPF_vs_pt");
  gr_responseMPFMC_vs_pt->SetMarkerStyle(24);
  gr_responseMPFMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPFMC_vs_pt->SetMarkerColor(38);


  std::string responseEXTRAP_name = (correctedPt) ? "gr_DATARespL2L3_vs_pt" : "gr_DATAResp_vs_pt";
  TGraphErrors* gr_responseEXTRAP_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAP_name.c_str());
  gr_responseEXTRAP_vs_pt->SetMarkerStyle(22);
  gr_responseEXTRAP_vs_pt->SetMarkerSize(markerSize);
  gr_responseEXTRAP_vs_pt->SetMarkerColor(46);

  std::string responseEXTRAPMC_name = (correctedPt) ? "gr_extrapRespL2L3_vs_pt" : "gr_extrapResp_vs_pt";
  TGraphErrors* gr_responseEXTRAPMC_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAPMC_name.c_str());
  gr_responseEXTRAPMC_vs_pt->SetMarkerStyle(26);
  gr_responseEXTRAPMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseEXTRAPMC_vs_pt->SetMarkerColor(46);

  
  //db->drawObjects( graphs, "all_vs_pt", "Photon p_{T} [GeV/c]", xMin, xMax, "Response", 0.65, 1.1, (bool)true );

  float ymin = (db->get_recoType()=="calo") ? 0.0 : 0.6;
  float ymax = 1.1;
  if( correctedPt ) {
    ymin = 0.7;
    ymax = 1.2;
  }

  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, ymin, ymax);
  axes->SetXTitle( "Photon p_{T} [GeV/c]" );
  if( correctedPt )
    axes->SetYTitle( "L2L3 Response" );
  else
    axes->SetYTitle( "Response" );
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();
  axes->GetXaxis()->SetTitleOffset(1.1);
  axes->GetYaxis()->SetTitleOffset(1.5);


  TLine* line_one = new TLine( xMin, 1., xMax, 1. );


  TLegend* legend = new TLegend( 0.48, 0.15, 0.87, 0.45, "  |#eta| < 1.3" );
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.033);
  legend->AddEntry( gr_responseBALANCING_vs_pt, "#gamma+Jet Balancing", "P");
  legend->AddEntry( gr_responseEXTRAP_vs_pt, "#gamma+Jet Extrapolation", "P");
  if( !correctedPt ) {
    legend->AddEntry( gr_responseMPF_vs_pt, "#gamma+Jet MPF", "P");
    legend->AddEntry( gr_responseMPFMC_vs_pt, "#gamma+Jet MPF (MC)", "P");
  }
  legend->AddEntry( band_responseGEN, "True Response", "FL");

  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();
  TPaveText* label_algo = db->get_labelAlgo();
  label_algo->SetTextSize(0.032);


  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetLogx();
  axes->Draw();
  if( correctedPt )
    line_one->Draw("same");
  band_responseGEN->Draw("C e3 same");
  fit_responseGEN->Draw("same");
  if( db->get_recoType()=="jpt" )
    gr_responseGEN_vs_pt->Draw("psame");
  gr_responseBALANCING_vs_pt->Draw("psame");
  if( !correctedPt ) {
    gr_responseMPF_vs_pt->Draw("psame");
    gr_responseMPFMC_vs_pt->Draw("psame");
  }
  gr_responseEXTRAP_vs_pt->Draw("psame");
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gPad->RedrawAxis();

  std::string name_base = db->get_outputdir() + "/all_";
  if( correctedPt ) name_base = name_base + "L2L3_";
  std::string name_eps = name_base + "vs_pt.eps";
  std::string name_png = name_base + "vs_pt.png";
  c1->SaveAs( name_eps.c_str() );
  c1->SaveAs( name_png.c_str() );


  // and now data/MC comparisons:

  c1->Clear();
  TH2D* axes2 = new TH2D( "axes2", "", 10, xMin, xMax, 10, 0.7, 1.2);
  axes2->SetXTitle( "Photon p_{T} [GeV/c]" );
  axes2->SetYTitle( "Data / MC" );
  axes2->GetXaxis()->SetMoreLogLabels();
  axes2->GetXaxis()->SetNoExponent();
  axes2->GetXaxis()->SetTitleOffset(1.1);
  axes2->GetYaxis()->SetTitleOffset(1.5);
  axes2->Draw();

  TGraphErrors* gr_dataMC_BALANCING = fitTools::get_graphRatio( gr_responseBALANCING_vs_pt, gr_responseBALANCINGMC_vs_pt);
  gr_dataMC_BALANCING->SetMarkerStyle(21); 
  gr_dataMC_BALANCING->SetMarkerSize(markerSize); 
  gr_dataMC_BALANCING->SetMarkerColor(kGray+2); 
  TGraphErrors* gr_dataMC_MPF = fitTools::get_graphRatio( gr_responseMPF_vs_pt, gr_responseMPFMC_vs_pt);
  gr_dataMC_MPF->SetMarkerStyle(20); 
  gr_dataMC_MPF->SetMarkerSize(markerSize); 
  gr_dataMC_MPF->SetMarkerColor(38);
  TGraphErrors* gr_dataMC_EXTRAP = fitTools::get_graphRatio( gr_responseEXTRAP_vs_pt, gr_responseEXTRAPMC_vs_pt);
  gr_dataMC_EXTRAP->SetMarkerStyle(22); 
  gr_dataMC_EXTRAP->SetMarkerSize(markerSize); 
  gr_dataMC_EXTRAP->SetMarkerColor(46);
  
  TF1* f_const_BALANCING = new TF1("const_BALANCING", "[0]", xMin, xMax);
  f_const_BALANCING->SetLineStyle(2);
  f_const_BALANCING->SetLineColor(kGray+2);
  gr_dataMC_BALANCING->Fit(f_const_BALANCING, "QR");

  TF1* f_const_MPF = new TF1("const_MPF", "[0]", xMin, xMax);
  f_const_MPF->SetLineStyle(2);
  f_const_MPF->SetLineColor(38);
  gr_dataMC_MPF->Fit(f_const_MPF, "QR");

  TF1* f_const_EXTRAP = new TF1("const_EXTRAP", "[0]", xMin, xMax);
  f_const_EXTRAP->SetLineStyle(2);
  f_const_EXTRAP->SetLineColor(46);
  gr_dataMC_EXTRAP->Fit(f_const_EXTRAP, "QR");

  char balancingText[200];
  sprintf( balancingText, "Balancing (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_BALANCING->GetParameter(0), f_const_BALANCING->GetParError(0), f_const_BALANCING->GetChisquare(), f_const_BALANCING->GetNDF() );
  char mpfText[200];
  sprintf( mpfText, "MPF (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_MPF->GetParameter(0), f_const_MPF->GetParError(0), f_const_MPF->GetChisquare(), f_const_MPF->GetNDF() );
  char extrapText[200];
  sprintf( extrapText, "Extrapolation (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_EXTRAP->GetParameter(0), f_const_EXTRAP->GetParError(0), f_const_EXTRAP->GetChisquare(), f_const_EXTRAP->GetNDF() );

  TLegend* legend2 = new TLegend( 0.14, 0.15, 0.7, 0.4, "   #gamma+Jet, |#eta| < 1.3" );
  legend2->SetFillColor(kWhite);
  legend2->SetTextSize(0.028);
  legend2->AddEntry( gr_dataMC_BALANCING, balancingText, "P");
  legend2->AddEntry( gr_dataMC_MPF, mpfText, "P");
  legend2->AddEntry( gr_dataMC_EXTRAP, extrapText, "P");


  TPaveText* label_algo2 = db->get_labelAlgo(1);
  label_algo2->SetTextSize(0.032);

  gr_dataMC_BALANCING->Draw("p same");
  gr_dataMC_MPF->Draw("p same");
  gr_dataMC_EXTRAP->Draw("p same");
  line_one->Draw("same");
  legend2->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo2->Draw("same");

  name_base = db->get_outputdir() + "/dataMC_";
  if( correctedPt ) name_base = name_base + "L2L3_";
  name_eps = name_base + "vs_pt.eps";
  name_png = name_base + "vs_pt.png";
  c1->SaveAs(name_eps.c_str());
  c1->SaveAs(name_png.c_str());

  file_noextrap->Close();
  file_extrap->Close();

  delete c1;
  c1 = 0;
  delete axes;
  axes = 0;
  delete legend;
  legend = 0;

}



