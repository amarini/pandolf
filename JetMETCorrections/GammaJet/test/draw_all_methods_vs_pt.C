#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TGraphErrors.h"

#include "DrawBase.h"
#include "fitTools.h"


bool ELIF_ = false;
bool FIXM_ = false;



void draw_vs_pt_plots( const std::string& resp_reso, const std::string& etaRegion, const std::string& FIT_RMS, DrawBase* db, bool correctedPt=false, const std::string& recoGen="RecoRel" );


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
  db->set_isCMSArticle( (bool)true );

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




  //db->set_lumi(33.8);
  db->set_lumiNormalization(36.);
  //db->set_shapeNormalization();

  db->set_outputdir();

  //std::string fit_rms = (db->get_recoType()=="calo") ? "FIT" : "RMS99";
  std::string fit_rms = "RMS99";

  draw_vs_pt_plots("response",   "eta011", fit_rms, db);
  draw_vs_pt_plots("resolution", "eta011", fit_rms, db);
  draw_vs_pt_plots("response",   "eta011", fit_rms, db, (bool)true);
  draw_vs_pt_plots("resolution", "eta011", fit_rms, db, (bool)true);
  draw_vs_pt_plots("response",   "eta011", fit_rms, db, (bool)true, "RecoRelRaw");
  draw_vs_pt_plots("resolution", "eta011", fit_rms, db, (bool)true, "RecoRelRaw");
  draw_vs_pt_plots("response",   "eta013", fit_rms, db);
  draw_vs_pt_plots("resolution", "eta013", fit_rms, db);
  draw_vs_pt_plots("response",   "eta013", fit_rms, db, (bool)true);
  draw_vs_pt_plots("resolution", "eta013", fit_rms, db, (bool)true);
  draw_vs_pt_plots("response",   "eta013", fit_rms, db, (bool)true, "RecoRelRaw");
  draw_vs_pt_plots("resolution", "eta013", fit_rms, db, (bool)true, "RecoRelRaw");
//draw_vs_pt_plots("response",   "eta243", "RMS99", db);
//draw_vs_pt_plots("resolution", "eta243", "RMS99", db);
//draw_vs_pt_plots("response",   "eta243", "RMS99", db, (bool)true);
//draw_vs_pt_plots("resolution", "eta243", "RMS99", db, (bool)true);
//draw_vs_pt_plots("response",   "eta23", "RMS99", db, (bool)true);
//draw_vs_pt_plots("resolution", "eta23", "RMS99", db, (bool)true);
//draw_vs_pt_plots("response",   "eta35", "RMS99", db);
//draw_vs_pt_plots("resolution", "eta35", "RMS99", db);
//draw_vs_pt_plots("response",   "eta35", "RMS99", db, (bool)true);
//draw_vs_pt_plots("resolution", "eta35", "RMS99", db, (bool)true);

  return 0;

}





void draw_vs_pt_plots( const std::string& resp_reso, const std::string& etaRegion, const std::string& FIT_RMS, DrawBase* db, bool correctedPt, const std::string& recoGen ) {

  std::string etaRegion_str;
  if( etaRegion=="eta013" ) etaRegion_str = "|#eta| < 1.3";
  else if( etaRegion=="eta011" ) etaRegion_str = "|#eta| < 1.1";
  else if( etaRegion=="eta009" ) etaRegion_str = "|#eta| < 0.9";
  else if( etaRegion=="eta132" ) etaRegion_str = "1.3 < |#eta| < 2";
  else if( etaRegion=="eta1524" ) etaRegion_str = "1.5 < |#eta| < 2.4";
  else if( etaRegion=="eta23" ) etaRegion_str = "2 < |#eta| < 3";
  else if( etaRegion=="eta243" ) etaRegion_str = "2.4 < |#eta| < 3";
  else if( etaRegion=="eta35" ) etaRegion_str = "3 < |#eta| < 5";
  else etaRegion_str = "[unknown eta region]";

  if( resp_reso!="response" && resp_reso!="resolution" ) {
    std::cout << "Only 'Response' and 'Resolution' supported. Exiting." << std::endl;
    exit(776);
  } 

  std::string file_noextrap_name = "PhotonJetGraphs_" + db->get_fullSuffix() + ".root";
  std::string file_extrap_name = "PhotonJetExtrapGraphs_" + db->get_fullSuffix();
  if( etaRegion!="" ) file_extrap_name += "_" + etaRegion;
  if( correctedPt ) {
    file_extrap_name += "L2L3";
    if( recoGen=="RecoRelRaw" ) file_extrap_name += "Raw";
  }
  file_extrap_name += "_" + FIT_RMS;
  if( FIXM_ ) file_extrap_name += "_FIXM";
  file_extrap_name += ".root";

  //now open graph files and plot them on on top of the other:
  TFile* file_noextrap = TFile::Open( file_noextrap_name.c_str(), "read");
  if( file_noextrap!=0 ) std::cout << "-> Opened file: '" << file_noextrap_name << "'" << std::endl;
  else {
    std::cout << "Didn't find file '" << file_noextrap_name << "'. Skipping." << std::endl;
    return;
  }
  TFile* file_extrap = TFile::Open( file_extrap_name.c_str(), "read");
  if( file_extrap!=0 ) std::cout << "-> Opened file: '" << file_extrap_name << "'" << std::endl;
  else {
    std::cout << "Didn't find file '" << file_extrap_name << "'. Skipping." << std::endl;
    return;
  }


  std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();
  float xMin = (resp_reso=="response") ? ptPhot_binning[1] : ptPhot_binning[2];
  //if( db->get_recoType()=="calo" ) xMin = ptPhot_binning[1];
  float xMax = ptPhot_binning[ ptPhot_binning.size()-3 ];

  int markerSize = 2.;


  std::string responseGEN_name = "GEN"+resp_reso;
  if( correctedPt ) responseGEN_name += "_L2L3";
  if( etaRegion!="" ) responseGEN_name += "_" + etaRegion;
  responseGEN_name += "_vs_pt";
  //std::cout << "Looking for '" << responseGEN_name << "' in file '" << file_noextrap->GetName() << "." << std::endl;
  //TGraphErrors* gr_responseGEN_vs_pt = (TGraphErrors*)file_noextrap->Get(responseGEN_name.c_str());
  std::string intrName = "gr_intr";
  if( resp_reso=="response" ) intrName += "Resp";
  else intrName += "Reso";
  intrName += "_vs_pt";
  std::cout << "Looking for '" << intrName << "' in file '" << file_extrap->GetName() << "." << std::endl;
  TGraphErrors* gr_responseGEN_vs_pt = (TGraphErrors*)file_extrap->Get(intrName.c_str());
  gr_responseGEN_vs_pt->SetMarkerStyle(29);
  std::cout << "Got it." << std::endl;
  gr_responseGEN_vs_pt->SetMarkerSize(markerSize);
  gr_responseGEN_vs_pt->SetMarkerColor(kBlack);
  if( /*resp_reso=="response" &&*/ (db->get_recoType()=="calo") ) {
    gr_responseGEN_vs_pt->RemovePoint(0); 
    gr_responseGEN_vs_pt->RemovePoint(0); 
  }

  std::string funcType;
  if(resp_reso=="response" ){
    funcType = (db->get_recoType()=="pf") ? "rpf" : "powerlaw";
  } else {
    funcType = (db->get_recoType()=="pf" /*|| db->get_recoType()=="jpt"*/) ? "NSCPF" : "NSC";
  }


  TF1* fit_responseGEN = ( resp_reso=="response" ) ? fitTools::fitResponseGraph(gr_responseGEN_vs_pt, funcType, "f1_responseGEN", "RN", 1000.)
                                                   : fitTools::fitResolutionGraph(gr_responseGEN_vs_pt, funcType, "f1_responseGEN", "RN", 1000.);
  fit_responseGEN->SetLineWidth(2.);
  fit_responseGEN->SetRange(xMin, xMax);
  // now get fit error band:
  TH1D* band_responseGEN = fitTools::getBand(fit_responseGEN, "band_responseGEN");
  band_responseGEN->SetFillColor(kYellow-9);
  band_responseGEN->SetLineWidth(2.);


  std::string responseBALANCING_name = resp_reso;
  if( correctedPt ) responseBALANCING_name += "_L2L3";
  if( etaRegion!="" ) responseBALANCING_name += "_" + etaRegion;
  responseBALANCING_name += "_vs_pt";
  std::cout << "Looking for '" << responseBALANCING_name << "'." << std::endl;
  TGraphErrors* gr_responseBALANCING_vs_pt = (TGraphErrors*)file_noextrap->Get(responseBALANCING_name.c_str());
  gr_responseBALANCING_vs_pt->SetMarkerStyle(21);
  gr_responseBALANCING_vs_pt->SetMarkerSize(markerSize);
  gr_responseBALANCING_vs_pt->SetMarkerColor(kGray+2);
  std::cout << "Got it." << std::endl;
std::cout << "responseBALANCING points: " << gr_responseBALANCING_vs_pt->GetN() << std::endl;

  std::string responseBALANCINGMC_name = "MC"+resp_reso;
  if( correctedPt ) responseBALANCINGMC_name += "_L2L3";
  if( etaRegion!="" ) responseBALANCINGMC_name += "_" + etaRegion;
  responseBALANCINGMC_name += "_vs_pt";
  std::cout << "Looking for '" << responseBALANCINGMC_name << "'." << std::endl;
  TGraphErrors* gr_responseBALANCINGMC_vs_pt = (TGraphErrors*)file_noextrap->Get(responseBALANCINGMC_name.c_str());
  gr_responseBALANCINGMC_vs_pt->SetMarkerStyle(25);
  gr_responseBALANCINGMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseBALANCINGMC_vs_pt->SetMarkerColor(kGray+2);
  std::cout << "Got it." << std::endl;
std::cout << "responseBALANCINGMC points: " << gr_responseBALANCINGMC_vs_pt->GetN() << std::endl;


  std::string responseMPF_name = resp_reso+"MPF";
  if( etaRegion!="" ) responseMPF_name += "_" + etaRegion;
  responseMPF_name += "_vs_pt";
  TGraphErrors* gr_responseMPF_vs_pt = (TGraphErrors*)file_noextrap->Get(responseMPF_name.c_str());
  gr_responseMPF_vs_pt->SetMarkerStyle(20);
  gr_responseMPF_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPF_vs_pt->SetMarkerColor(38);

  std::string responseMPFMC_name = "MC"+responseMPF_name;
  TGraphErrors* gr_responseMPFMC_vs_pt = (TGraphErrors*)file_noextrap->Get(responseMPFMC_name.c_str());
  gr_responseMPFMC_vs_pt->SetMarkerStyle(24);
  gr_responseMPFMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPFMC_vs_pt->SetMarkerColor(38);


  std::string resp_reso_short = (resp_reso=="response") ? "Resp" : "Reso";

  //std::string responseEXTRAP_name = (correctedPt) ? "gr_DATA"+resp_reso_short+"L2L3_vs_pt" : "gr_DATA"+resp_reso_short+"_vs_pt";
  std::string responseEXTRAP_name = "gr_DATA"+resp_reso_short+"_vs_pt";
  TGraphErrors* gr_responseEXTRAP_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAP_name.c_str());
//if( ELIF_ )
//  gr_responseEXTRAP_vs_pt->SetMarkerStyle(21);
//else
    gr_responseEXTRAP_vs_pt->SetMarkerStyle(22);
  gr_responseEXTRAP_vs_pt->SetMarkerSize(markerSize);
  gr_responseEXTRAP_vs_pt->SetMarkerColor(46);
  gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove first point (cant extrapolate at such low pt)
  if( resp_reso=="resolution" && !ELIF_ ) gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove second point also
  //if( (db->get_recoType()=="calo")||(db->get_recoType()=="jpt") )
  if( db->get_recoType()=="calo" )
    gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove also third point for calo

  //std::string responseEXTRAPMC_name = (correctedPt) ? "gr_extrap"+resp_reso_short+"L2L3_vs_pt" : "gr_extrap"+resp_reso_short+"_vs_pt";
  std::string responseEXTRAPMC_name = "gr_extrap"+resp_reso_short+"_vs_pt";
  TGraphErrors* gr_responseEXTRAPMC_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAPMC_name.c_str());
  gr_responseEXTRAPMC_vs_pt->SetMarkerStyle(26);
  gr_responseEXTRAPMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseEXTRAPMC_vs_pt->SetMarkerColor(46);
  gr_responseEXTRAPMC_vs_pt->RemovePoint(0); //remove first point (cant extrapolate at such low pt)
  if( resp_reso=="resolution" && !ELIF_ ) gr_responseEXTRAPMC_vs_pt->RemovePoint(0); //remove second point also
  //if( (db->get_recoType()=="calo")||(db->get_recoType()=="jpt") )
  if( db->get_recoType()=="calo" )
    gr_responseEXTRAPMC_vs_pt->RemovePoint(0); //remove also second point for calo

  std::string responseELIF_name = "gr_DATAReso_subtr_vs_pt";
  TGraphErrors* gr_responseELIF_vs_pt = (TGraphErrors*)file_extrap->Get(responseELIF_name.c_str());
  gr_responseELIF_vs_pt->SetMarkerStyle(25);
  gr_responseELIF_vs_pt->SetMarkerSize(markerSize);
  gr_responseELIF_vs_pt->SetMarkerColor(kMagenta+4);
  gr_responseELIF_vs_pt->RemovePoint(0); //remove first point (cant extrapolate at such low pt)
  

  float ymin, ymax;
  if( resp_reso=="response" ) {
    ymin = (db->get_recoType()=="calo") ? 0.0 : 0.7;
    ymax = (db->get_recoType()=="jpt") ? 1.15 : 1.05;
    if( correctedPt ) {
      ymin = 0.7;
      ymax = 1.25;
    //ymin = (db->get_recoType()=="calo") ? 0.3 : 0.7;
    //ymax = (db->get_recoType()=="calo") ? 1.3 : 1.2;
    }
  } else {
    ymin = 0.;
    ymax = 0.3;
    if( db->get_recoType()=="calo" && !correctedPt ) ymax = 0.6;
  }

  std::string plotVarName = (resp_reso=="response") ? "Response" : "Resolution";

  TH2D* axes = new TH2D( "axes", "", 10, xMin, xMax, 10, ymin, ymax);
  axes->SetXTitle( "Photon p_{T} [GeV/c]" );
  std::string yTitle = plotVarName;
  if( correctedPt ) { 
    if( resp_reso=="response" ) yTitle = "Corrected " + yTitle;
    else yTitle = "Jet p_{T} Resolution";
  }
  axes->SetYTitle( yTitle.c_str() );
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();


  TLine* line_one = new TLine( xMin, 1., xMax, 1. );

  float legend_xmin, legend_xmax, legend_ymin, legend_ymax;
  if( resp_reso=="response" ) {
    legend_xmin = 0.52;
    legend_ymin = 0.2;
    legend_xmax = 0.9;
    legend_ymax = (correctedPt) ? 0.45 : 0.5;
  } else {
    if( db->get_recoType()=="calo" ) {
      legend_xmin = 0.17;
      legend_ymin = 0.17;
      legend_xmax = 0.57;
      legend_ymax = 0.385;
    } else {
      legend_xmin = 0.40;
      legend_ymin = 0.615;
      legend_xmax = 0.80;
      legend_ymax = 0.83;
    }
  }


  //bool drawStars = !( resp_reso=="response" &&  db->get_recoType()!="jpt" ) || ( etaRegion!="eta011" && etaRegion!="eta013" && etaRegion!="eta1524" );
  bool drawStars = ( resp_reso=="response" &&  db->get_recoType()=="jpt" ) || ( etaRegion!="eta011" && etaRegion!="eta013" && etaRegion!="eta1524" );

  std::string legendTitle = "  " + etaRegion_str;
  TLegend* legend = new TLegend( legend_xmin, legend_ymin, legend_xmax, legend_ymax, legendTitle.c_str());
  legend->SetFillColor(kWhite);
  //legend->SetFillStyle(0);
  legend->SetTextSize(0.033);
  if( resp_reso=="response" ) {
    legend->AddEntry( gr_responseBALANCING_vs_pt, "#gamma+Jet Balancing", "P");
    if( !correctedPt ) {
      legend->AddEntry( gr_responseMPF_vs_pt, "#gamma+Jet MPF", "P");
      //legend->AddEntry( gr_responseMPFMC_vs_pt, "#gamma+Jet MPF (MC)", "P");
    }
  }
  legend->AddEntry( gr_responseEXTRAP_vs_pt, "#gamma+Jet Extrapolation", "P");
  legend->AddEntry( gr_responseEXTRAPMC_vs_pt, "#gamma+Jet Extrap. (MC)", "P");
  if( ELIF_ && resp_reso=="resolution" )
    legend->AddEntry( gr_responseELIF_vs_pt, "#gamma+Jet Subtraction", "P");
  std::string legendTrue = "True " + plotVarName;
  if( drawStars ) 
    legend->AddEntry( gr_responseGEN_vs_pt, legendTrue.c_str(), "P");
  else
    //legend->AddEntry( band_responseGEN, legendTrue.c_str(), "FL");
    legend->AddEntry( fit_responseGEN, legendTrue.c_str(), "L");

  TPaveText* label_cms = db->get_labelCMS(0);
  TPaveText* label_sqrt = db->get_labelSqrt(0);
  TPaveText* label_algo = (db->get_recoType()=="calo" && resp_reso=="resolution") ? db->get_labelAlgo(2) : db->get_labelAlgo();
  label_algo->SetTextSize(0.032);


  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  c1->SetLogx();
  axes->Draw();
  if( correctedPt )
    line_one->Draw("same");
  if( drawStars ) {
    gr_responseGEN_vs_pt->Draw("psame");
  } else {
    //band_responseGEN->Draw("C e3 same");
    fit_responseGEN->Draw("same");
  }
  //gr_responseGEN_vs_pt->Draw("psame");
//if( db->get_recoType()=="jpt" )
//  gr_responseGEN_vs_pt->Draw("psame");

  if( resp_reso=="response" ) {
    gr_responseBALANCING_vs_pt->Draw("psame");
    if( !correctedPt ) {
      gr_responseMPF_vs_pt->Draw("psame");
      //gr_responseMPFMC_vs_pt->Draw("psame");
    }
  }
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");
  gr_responseEXTRAPMC_vs_pt->Draw("psame");
  gr_responseEXTRAP_vs_pt->Draw("psame");
  if( ELIF_ && resp_reso=="resolution" )
    gr_responseELIF_vs_pt->Draw("psame");

  gPad->RedrawAxis();

  std::string name_base = db->get_outputdir();
  if( FIXM_ ) name_base += "_FIXM";
  name_base += "/" + resp_reso;
  name_base = name_base + FIT_RMS;
  if( correctedPt ) {
    name_base = name_base + "_L2L3";
    if( recoGen=="RecoRelRaw" ) name_base += "Raw";
  }
  if( etaRegion!="" ) name_base = name_base + "_" + etaRegion;
  std::string thisName = name_base + "_all_";
  std::string name_eps = thisName + "vs_pt.eps";
  std::string name_png = thisName + "vs_pt.png";
  c1->SaveAs( name_eps.c_str() );
  c1->SaveAs( name_png.c_str() );


  // and now data/MC comparisons:

  c1->Clear();
  float dataMC_ymin = (resp_reso=="response") ? 0.65 : 0.;
  float dataMC_ymax = (resp_reso=="response") ? 1.2 : 2.5;
  //if( resp_reso=="resolution" && db->get_recoType()=="calo") dataMC_ymax=3.;
  //float xMin_dataMC = (ELIF_) ? 22. : xMin;
  float xMin_dataMC = xMin;
  TH2D* axes2 = new TH2D( "axes2", "", 10, xMin_dataMC, xMax, 10, dataMC_ymin, dataMC_ymax);
  axes2->SetXTitle( "Photon p_{T} [GeV/c]" );
  axes2->SetYTitle( "Data / MC" );
  axes2->GetXaxis()->SetMoreLogLabels();
  axes2->GetXaxis()->SetNoExponent();
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

  TGraphErrors* gr_dataMC_ELIF=0;
  if( ELIF_ && resp_reso=="resolution" ) {
 // TFile* file_elif = TFile::Open("PhotonJet_elif.root");
 // gr_dataMC_ELIF = (TGraphErrors*) file_elif->Get("DataOverMCRatio");
    gr_dataMC_ELIF = (TGraphErrors*)file_extrap->Get("gr_reso_ratio_vs_pt");
    //gr_dataMC_ELIF = fitTools::get_graphRatio( gr_responseELIF_vs_pt, gr_MC_ELIF);
    gr_dataMC_ELIF->SetMarkerStyle(21); 
    gr_dataMC_ELIF->SetMarkerSize(markerSize); 
    gr_dataMC_ELIF->SetMarkerColor(kGray+2); 
  }

//// temporay fix!!!
//if( db->get_recoType()=="calo" ) 
//  gr_dataMC_EXTRAP->RemovePoint(gr_dataMC_EXTRAP->GetN()-1);
  
  float xMax_fit = xMax;

  TF1* f_const_BALANCING = new TF1("const_BALANCING", "[0]", xMin, xMax_fit);
  f_const_BALANCING->SetLineStyle(2);
  f_const_BALANCING->SetLineColor(kGray+2);
  gr_dataMC_BALANCING->Fit(f_const_BALANCING, "QR");

  TF1* f_const_MPF = new TF1("const_MPF", "[0]", xMin, xMax_fit);
  f_const_MPF->SetLineStyle(2);
  f_const_MPF->SetLineColor(38);
  gr_dataMC_MPF->Fit(f_const_MPF, "QR");

  if( etaRegion=="eta243"||etaRegion=="eta1524" ) xMax_fit = 200.;
  TF1* f_const_EXTRAP = new TF1("const_EXTRAP", "[0]", xMin_dataMC, xMax_fit);
  //TF1* f_const_EXTRAP = new TF1("const_EXTRAP", "[0]", 20., 150.);
  f_const_EXTRAP->SetParameter(0, 1.);
  f_const_EXTRAP->SetLineStyle(2);
  f_const_EXTRAP->SetLineColor(46);
  gr_dataMC_EXTRAP->Fit(f_const_EXTRAP, "QR");

  TF1* f_const_ELIF = new TF1("const_ELIF", "[0]", xMin_dataMC, xMax_fit);
  //TF1* f_const_EXTRAP = new TF1("const_EXTRAP", "[0]", 20., 150.);
  f_const_ELIF->SetParameter(0, 1.);
  f_const_ELIF->SetLineStyle(2);
  f_const_ELIF->SetLineColor(kGray+2);
  if( ELIF_ && resp_reso=="resolution" ) gr_dataMC_ELIF->Fit(f_const_ELIF, "QR");


  // get syst band from file:
  std::string systFile_name = "totalSyst_";
  if( resp_reso=="response") systFile_name += "resp";
  else                       systFile_name += "reso";
  systFile_name += "_"+db->get_algoType();
  systFile_name += ".root";
  TFile* file_syst = TFile::Open(systFile_name.c_str());
  TH1D* syst_band = (file_syst!=0) ? (TH1D*)file_syst->Get("syst_total") : 0;
  if( syst_band!=0 ) {
    syst_band->SetFillColor(kYellow-9);
    syst_band->SetFillStyle(1001);
    syst_band->SetLineColor(kBlack);
    syst_band->SetLineWidth(1);
    for( unsigned ibin_syst=0; ibin_syst<syst_band->GetNbinsX(); ++ibin_syst ) {
      syst_band->SetBinError( ibin_syst+1, syst_band->GetBinContent(ibin_syst+1)/100.); //not in percent
      syst_band->SetBinContent( ibin_syst+1, 1. ); //put it around one
    }
  }

  char balancingText[400];
  sprintf( balancingText, "Balancing (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_BALANCING->GetParameter(0), f_const_BALANCING->GetParError(0), f_const_BALANCING->GetChisquare(), f_const_BALANCING->GetNDF() );
  char mpfText[400];
  sprintf( mpfText, "MPF (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_MPF->GetParameter(0), f_const_MPF->GetParError(0), f_const_MPF->GetChisquare(), f_const_MPF->GetNDF() );
  char extrapText[400];
  if( ELIF_ && resp_reso=="resolution" )
    sprintf( extrapText, "Direct Extrap. (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_EXTRAP->GetParameter(0), f_const_EXTRAP->GetParError(0), f_const_EXTRAP->GetChisquare(), f_const_EXTRAP->GetNDF() );
  else
    sprintf( extrapText, "Extrapolation (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_EXTRAP->GetParameter(0), f_const_EXTRAP->GetParError(0), f_const_EXTRAP->GetChisquare(), f_const_EXTRAP->GetNDF() );
  char elifText[400];
  if( ELIF_ && resp_reso=="resolution" )
    sprintf( elifText, "Ratio Method (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_ELIF->GetParameter(0), f_const_ELIF->GetParError(0), f_const_ELIF->GetChisquare(), f_const_ELIF->GetNDF() );

  float legend_yMax2 = (ELIF_) ? 0.4 : 0.3;
  if( resp_reso=="response" ) legend_yMax2 = 0.4;

  legendTitle = "   #gamma+Jet, " + etaRegion_str;
  TLegend* legend2 = new TLegend( 0.15, 0.2, 0.73, legend_yMax2, legendTitle.c_str());
  legend2->SetFillColor(kWhite);
  //legend2->SetTextFont(42);
  legend2->SetTextSize(0.03);
  //legend2->AddEntry( gr_dataMC_BALANCING, balancingText, "P");
  if( resp_reso=="response" ) {
    legend2->AddEntry( gr_dataMC_BALANCING, balancingText, "P");
    if( !correctedPt )
      legend2->AddEntry( gr_dataMC_MPF, mpfText, "P");
  }
  legend2->AddEntry( gr_dataMC_EXTRAP, extrapText, "P");
  if( ELIF_ && resp_reso=="resolution" )
    legend2->AddEntry( gr_dataMC_ELIF, elifText, "P");

  //TLegend* legend_syst = new TLegend( 0.53, 0.7, 0.88, 0.88, "Anti-k_{T} 0.5 PFJets");
  float yMin_syst = (resp_reso=="response") ? 0.43 : 0.63;
  float yMax_syst = (resp_reso=="response") ? 0.50 : 0.7;
  TLegend* legend_syst = new TLegend( 0.35, yMin_syst, 0.65, yMax_syst, "");
  legend_syst->SetFillColor(kWhite);
  legend_syst->SetTextFont(42);
  legend_syst->SetTextSize(0.032);
  if( syst_band!=0 ) {
    if( resp_reso=="response" )
      legend_syst->AddEntry(syst_band, "Extrap. Syst. Uncertainty", "F");
    else
      legend_syst->AddEntry(syst_band, "Syst. Uncertainty", "F");
  }


  TPaveText* label_algo2 = db->get_labelAlgo(1);
  label_algo2->SetTextSize(0.032);

  legend2->Draw("same");
  if( ELIF_ && resp_reso=="resolution" ) {
    if( syst_band!=0 ) {
      syst_band->Draw("c l e3 same");
      legend_syst->Draw("same");
    }
  }
  if( ELIF_ ) {
    TLine* line_one2 = new TLine( xMin_dataMC, 1., xMax, 1. );
    line_one2->Draw("same");
  } else {
    line_one->Draw("same");
  }
  //gr_dataMC_BALANCING->Draw("p same");
  if( resp_reso=="response" ) {
    gr_dataMC_BALANCING->Draw("p same");
    if( !correctedPt )
      gr_dataMC_MPF->Draw("p same");
  }
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo2->Draw("same");
  gr_dataMC_EXTRAP->Draw("p same");
  if( ELIF_ && resp_reso=="resolution" ) gr_dataMC_ELIF->Draw("p same");
  gPad->RedrawAxis();

//name_base = db->get_outputdir() + "/" + resp_reso;
//name_base = name_base + FIT_RMS;
//if( correctedPt ) {
//  name_base = name_base + "_L2L3";
//  if( recoGen=="RecoRelRaw" ) name_base += "Raw";
//}
//if( etaRegion!="" ) name_base = name_base + "_" + etaRegion;
  thisName = name_base + "_dataMC_";
  if( ELIF_ && resp_reso=="resolution" ) thisName += "ELIF_";
  name_eps = thisName + "vs_pt.eps";
  name_png = thisName + "vs_pt.png";
  c1->SaveAs(name_eps.c_str());
  c1->SaveAs(name_png.c_str());



  if( resp_reso=="resolution" ) {

    // QUADRATURE DIFFERENCE BETWEEN DATA AND MC

    TGraphErrors* gr_squareDiff_BALANCING = new TGraphErrors(0);
    TGraphErrors* gr_squareDiff_MPF = new TGraphErrors(0);
    TGraphErrors* gr_squareDiff_EXTRAP = new TGraphErrors(0);

    for( unsigned iPoint=0; iPoint<gr_dataMC_BALANCING->GetN(); ++iPoint ) {

      Double_t r_BALANCING, x_BALANCING;
      gr_responseBALANCING_vs_pt->GetPoint(iPoint, x_BALANCING, r_BALANCING);
      Double_t r_BALANCINGMC, x_BALANCINGMC;
      gr_responseBALANCINGMC_vs_pt->GetPoint(iPoint, x_BALANCINGMC, r_BALANCINGMC);
      Double_t squareDiff_BALANCING = sqrt( fabs(r_BALANCING*r_BALANCING - r_BALANCINGMC*r_BALANCINGMC) );
      Double_t r_BALANCING_err = gr_responseBALANCING_vs_pt->GetErrorY(iPoint);
      Double_t r_BALANCINGMC_err = gr_responseBALANCINGMC_vs_pt->GetErrorY(iPoint);
      Double_t squareDiff_BALANCING_err = sqrt( r_BALANCING_err*r_BALANCING_err + r_BALANCINGMC_err*r_BALANCINGMC_err );
      gr_squareDiff_BALANCING->SetPoint( iPoint, x_BALANCINGMC, squareDiff_BALANCING );
      gr_squareDiff_BALANCING->SetPointError( iPoint, 0., squareDiff_BALANCING_err );

      Double_t r_MPF, x_MPF;
      gr_responseMPF_vs_pt->GetPoint(iPoint, x_MPF, r_MPF);
      Double_t r_MPFMC, x_MPFMC;
      gr_responseMPFMC_vs_pt->GetPoint(iPoint, x_MPFMC, r_MPFMC);
      Double_t squareDiff_MPF = sqrt( fabs(r_MPF*r_MPF - r_MPFMC*r_MPFMC) );
      Double_t r_MPF_err = gr_responseMPF_vs_pt->GetErrorY(iPoint);
      Double_t r_MPFMC_err = gr_responseMPFMC_vs_pt->GetErrorY(iPoint);
      Double_t squareDiff_MPF_err = sqrt( r_MPF_err*r_MPF_err + r_MPFMC_err*r_MPFMC_err );
      gr_squareDiff_MPF->SetPoint( iPoint, x_MPFMC, squareDiff_MPF );
      gr_squareDiff_MPF->SetPointError( iPoint, 0., squareDiff_MPF_err );

      Double_t r_EXTRAP, x_EXTRAP;
      gr_responseEXTRAP_vs_pt->GetPoint(iPoint, x_EXTRAP, r_EXTRAP);
      Double_t r_EXTRAPMC, x_EXTRAPMC;
      gr_responseEXTRAPMC_vs_pt->GetPoint(iPoint, x_EXTRAPMC, r_EXTRAPMC);
      Double_t squareDiff_EXTRAP = sqrt( fabs(r_EXTRAP*r_EXTRAP - r_EXTRAPMC*r_EXTRAPMC) );
      Double_t r_EXTRAP_err = gr_responseEXTRAP_vs_pt->GetErrorY(iPoint);
      Double_t r_EXTRAPMC_err = gr_responseEXTRAPMC_vs_pt->GetErrorY(iPoint);
      Double_t squareDiff_EXTRAP_err = sqrt( r_EXTRAP_err*r_EXTRAP_err + r_EXTRAPMC_err*r_EXTRAPMC_err );
      gr_squareDiff_EXTRAP->SetPoint( iPoint, x_EXTRAPMC, squareDiff_EXTRAP );
      gr_squareDiff_EXTRAP->SetPointError( iPoint, 0., squareDiff_EXTRAP_err );


    }

    gr_squareDiff_BALANCING->SetMarkerStyle(21); 
    gr_squareDiff_BALANCING->SetMarkerSize(markerSize); 
    gr_squareDiff_BALANCING->SetMarkerColor(kGray+2); 

    gr_squareDiff_MPF->SetMarkerStyle(20); 
    gr_squareDiff_MPF->SetMarkerSize(markerSize); 
    gr_squareDiff_MPF->SetMarkerColor(38);

    gr_squareDiff_EXTRAP->SetMarkerStyle(22); 
    gr_squareDiff_EXTRAP->SetMarkerSize(markerSize); 
    gr_squareDiff_EXTRAP->SetMarkerColor(46);

  //gr_squareDiff_BALANCING->Fit(f_const_BALANCING, "QR");
  //gr_squareDiff_MPF->Fit(f_const_MPF, "QR");
  //gr_squareDiff_EXTRAP->Fit(f_const_EXTRAP, "QR");

    TLegend* legend_squareDiff = new TLegend(0.5, 0.55, 0.85, 0.75);
    legend_squareDiff->SetTextSize(0.035);
    legend_squareDiff->SetFillColor(0);
    //legend_squareDiff->SetFillStyle(1);
    legend_squareDiff->AddEntry( gr_squareDiff_BALANCING, "Balancing", "P" );
    legend_squareDiff->AddEntry( gr_squareDiff_MPF, "MPF", "P" );
    legend_squareDiff->AddEntry( gr_squareDiff_EXTRAP, "Extrapolation", "P" );

    sprintf( balancingText, "Balancing (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_BALANCING->GetParameter(0), f_const_BALANCING->GetParError(0), f_const_BALANCING->GetChisquare(), f_const_BALANCING->GetNDF() );
    sprintf( mpfText, "MPF (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_MPF->GetParameter(0), f_const_MPF->GetParError(0), f_const_MPF->GetChisquare(), f_const_MPF->GetNDF() );
    sprintf( extrapText, "Extrapolation (FIT = %.3lf #pm %.3lf, #chi^{2}/NDF = %.2lf/%d)", f_const_EXTRAP->GetParameter(0), f_const_EXTRAP->GetParError(0), f_const_EXTRAP->GetChisquare(), f_const_EXTRAP->GetNDF() );

    TLegend* legend3 = new TLegend( 0.14, 0.15, 0.7, 0.4, legendTitle.c_str());
    legend3->SetFillColor(kWhite);
    legend3->SetTextSize(0.028);
    legend3->AddEntry( gr_squareDiff_BALANCING, balancingText, "P");
    legend3->AddEntry( gr_squareDiff_MPF, mpfText, "P");
    legend3->AddEntry( gr_squareDiff_EXTRAP, extrapText, "P");

    TH2D* axes3 = new TH2D( "axes3", "", 10, xMin, xMax, 10, 0., 0.3);
    axes3->SetXTitle( "Photon p_{T} [GeV/c]" );
    axes3->SetYTitle( "#sqrt{Data^{2} - MC^{2}}" );
    axes3->GetXaxis()->SetMoreLogLabels();
    axes3->GetXaxis()->SetNoExponent();
    axes3->Draw();

    gr_squareDiff_BALANCING->Draw("p same");
    gr_squareDiff_MPF->Draw("p same");
    gr_squareDiff_EXTRAP->Draw("p same");
    //legend3->Draw("same");
    legend_squareDiff->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo2->Draw("same");

  //name_base = db->get_outputdir() + "/" + resp_reso;
  //name_base = name_base + FIT_RMS;
  //if( correctedPt ) {
  //  name_base = name_base + "_L2L3";
  //  if( recoGen=="RecoRelRaw" ) name_base += "Raw";
  //}
  //if( etaRegion!="" ) name_base = name_base + "_" + etaRegion;
    thisName = name_base + "_squareDiff_vs_pt";
    name_eps = thisName + ".eps";
    name_png = thisName + ".png";
    c1->SaveAs(name_eps.c_str());
    c1->SaveAs(name_png.c_str());



    // compare NSC fit parameters:

    c1->Clear();

    axes->GetYaxis()->SetRangeUser(0., 0.35);
    axes->Draw();

    std::string reso_func_type = (db->get_recoType()=="calo") ? "NSC" : "NSCPF";

    TF1* fit_responseGEN_2   = fitTools::fitResolutionGraph(gr_responseGEN_vs_pt, reso_func_type, "f1_responseGEN_2", "RN");
    //fit_responseGEN_2->SetLineColor(46);
    fit_responseGEN_2->SetLineWidth(1);
    fit_responseGEN_2->SetLineStyle(2);

    TF1* fit_responseEXTRAP   = fitTools::fitResolutionGraph(gr_responseEXTRAP_vs_pt, reso_func_type, "f1_responseEXTRAP", "RN", 1400., 20.);
    fit_responseEXTRAP->SetLineColor(46);
    fit_responseEXTRAP->SetLineWidth(2);
    fit_responseEXTRAP->SetLineStyle(1);
    TH1D* band_responseEXTRAP = fitTools::getBand(fit_responseEXTRAP, "band_responseEXTRAP");
    band_responseEXTRAP->SetFillColor(kYellow-9);
    band_responseEXTRAP->SetLineWidth(2.);


    TF1* fit_responseEXTRAPMC = fitTools::fitResolutionGraph(gr_responseEXTRAPMC_vs_pt, reso_func_type, "f1_responseEXTRAPMC", "RN");
    fit_responseEXTRAPMC->SetLineColor(46);
    fit_responseEXTRAPMC->SetLineWidth(1);
    fit_responseEXTRAPMC->SetLineStyle(2);

    char labelText[200];

    TPaveText* labelGEN = new TPaveText(0.15, 0.67, 0.5, 0.88, "brNDC");
    labelGEN->SetTextSize(0.035);
    labelGEN->AddText("MC Truth");
    labelGEN->SetFillColor(0);
    labelGEN->SetFillStyle(0);
    sprintf( labelText, "N = %.4f #pm %.4f", fit_responseGEN_2->GetParameter(0), fit_responseGEN_2->GetParError(0) );
    labelGEN->AddText(labelText);
    sprintf( labelText, "S = %.4f #pm %.4f", fit_responseGEN_2->GetParameter(1), fit_responseGEN_2->GetParError(1) );
    labelGEN->AddText(labelText);
    sprintf( labelText, "C = %.4f #pm %.4f", fit_responseGEN_2->GetParameter(2), fit_responseGEN_2->GetParError(2) );
    labelGEN->AddText(labelText);

    TPaveText* labelEXTRAP = new TPaveText(0.51, 0.67, 0.88, 0.88, "brNDC");
    labelEXTRAP->SetTextSize(0.035);
    labelEXTRAP->SetTextColor(46);
    labelEXTRAP->SetFillColor(0);
    labelEXTRAP->SetFillStyle(0);
    labelEXTRAP->AddText("Extrapolation Data");
    sprintf( labelText, "N = %.3f #pm %.3f", fit_responseEXTRAP->GetParameter(0), fit_responseEXTRAP->GetParError(0) );
    labelEXTRAP->AddText(labelText);
    sprintf( labelText, "S = %.3f #pm %.3f", fit_responseEXTRAP->GetParameter(1), fit_responseEXTRAP->GetParError(1) );
    labelEXTRAP->AddText(labelText);
    sprintf( labelText, "C = %.3f #pm %.3f", fit_responseEXTRAP->GetParameter(2), fit_responseEXTRAP->GetParError(2) );
    labelEXTRAP->AddText(labelText);

    TPaveText* labelEXTRAPMC = new TPaveText(0.51, 0.46, 0.88, 0.66, "brNDC");
    labelEXTRAPMC->SetTextSize(0.035);
    labelEXTRAPMC->SetTextColor(46);
    labelEXTRAPMC->SetFillColor(0);
    labelEXTRAPMC->SetFillStyle(0);
    labelEXTRAPMC->AddText("Extrapolation MC");
    sprintf( labelText, "N = %.3f #pm %.3f", fit_responseEXTRAPMC->GetParameter(0), fit_responseEXTRAPMC->GetParError(0) );
    labelEXTRAPMC->AddText(labelText);
    sprintf( labelText, "S = %.3f #pm %.3f", fit_responseEXTRAPMC->GetParameter(1), fit_responseEXTRAPMC->GetParError(1) );
    labelEXTRAPMC->AddText(labelText);
    sprintf( labelText, "C = %.3f #pm %.3f", fit_responseEXTRAPMC->GetParameter(2), fit_responseEXTRAPMC->GetParError(2) );
    labelEXTRAPMC->AddText(labelText);

    gr_responseEXTRAP_vs_pt->Draw("Psame");
    gr_responseEXTRAPMC_vs_pt->Draw("Psame");
    gr_responseGEN_vs_pt->Draw("Psame");
    fit_responseEXTRAP->Draw("Lsame");
    fit_responseEXTRAPMC->Draw("Lsame");
    fit_responseGEN_2->Draw("Lsame");
    labelGEN->Draw("same");
    labelEXTRAP->Draw("same");
    labelEXTRAPMC->Draw("same");

  //name_base = db->get_outputdir() + "/" + resp_reso;
  //name_base = name_base + FIT_RMS;
  //if( correctedPt ) {
  //  name_base = name_base + "_L2L3";
  //  if( recoGen=="RecoRelRaw" ) name_base += "Raw";
  //}
  //if( etaRegion!="" ) name_base = name_base + "_" + etaRegion;
    thisName = name_base + "_NSC_vs_pt";
    if( FIXM_ ) thisName += "_FIXM";
    name_eps = thisName + ".eps";
    name_png = thisName + ".png";
    c1->SaveAs(name_eps.c_str());
    c1->SaveAs(name_png.c_str());


    // draw stat component vs. pt
    c1->Clear();

    TH2D* axes4 = new TH2D( "axes4", "", 10, 20., 1400., 10, 0., 30.);
    axes4->SetXTitle( "Photon p_{T} [GeV/c]" );
    axes4->SetYTitle( "Relative Uncertainty [%]");
    axes4->GetXaxis()->SetMoreLogLabels();
    axes4->GetXaxis()->SetNoExponent();
    axes4->Draw();

    TH1D* stat_errorEXTRAP = new TH1D(*band_responseEXTRAP);
    for( unsigned iBin=1; iBin<stat_errorEXTRAP->GetNbinsX(); ++iBin) {

      float pt = stat_errorEXTRAP->GetBinCenter(iBin);
      float reso_fit = fit_responseEXTRAP->Eval(pt);
      float reso_fit_err = band_responseEXTRAP->GetBinError(iBin);

      float syst_reso = 100.*reso_fit_err/reso_fit;

//std::cout << std::endl << "pt: " << pt << std::endl;
//std::cout << "reso_fit: " << reso_fit << std::endl;
//std::cout << "reso_fit_err: " << reso_fit_err << std::endl;
//std::cout << "syst_reso: " << syst_reso << std::endl;
      stat_errorEXTRAP->SetBinContent( iBin, syst_reso );

    }
    stat_errorEXTRAP->SetLineColor(46);
    stat_errorEXTRAP->SetLineWidth(2);
    stat_errorEXTRAP->SetFillColor(46);
    stat_errorEXTRAP->SetFillStyle(3004);
    
    axes4->Draw();
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    stat_errorEXTRAP->Draw("HC same");
    gPad->RedrawAxis();

  //name_base = db->get_outputdir() + "/" + resp_reso;
  //name_base = name_base + FIT_RMS;
  //if( correctedPt ) {
  //  name_base = name_base + "_L2L3";
  //  if( recoGen=="RecoRelRaw" ) name_base += "Raw";
  //}
  //if( etaRegion!="" ) name_base = name_base + "_" + etaRegion;
    thisName = name_base + "_STAT_vs_pt";
    name_eps = thisName + ".eps";
    name_png = thisName + ".png";
    c1->SaveAs(name_eps.c_str());
    c1->SaveAs(name_png.c_str());

  }

  file_noextrap->Close();
  file_extrap->Close();

  delete c1;
  c1 = 0;
  delete axes;
  axes = 0;
  delete legend;
  legend = 0;

}



