#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"


bool useMCassoc_ = false;
bool ONEVTX = false;



std::pair<TH1D,TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db );
//void drawSignalPtMix( std::pair<TH1D*,TH1D*> h1pair_3050, std::pair<TH1D*,TH1D*> h1pair_5080, std::pair<TH1D*,TH1D*> h1pair_80120, int mass, float frac_3050, float frac_5080, float frac_80120, DrawBase* db );


int main(int argc, char* argv[]) {

  if( argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawQGStudies [data_dataset] [photonID=\"medium\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset = "Photon_Run2011A-May10ReReco-v1";
  if( argc>1 ) {
    std::string dataset_tmp(argv[1]);
    data_dataset = dataset_tmp;
  }

  //std::string data_dataset = "DATA_Run2011A_1fb";
  std::string mc_photonjet = "G_Summer11";
  //std::string mc_photonjet = "G_Spring11";
  //std::string mc_photonjet = "G_Pt-80to120_Tune23_7TeV_herwigpp_Summer11-PU_S3_START42_V11-v2";
  //std::string mc_QCD = mc_photonjet;
  std::string mc_QCD = "QCD_EMEnriched_Spring11";
  //std::string mc_QCD = "QCD_EMEnriched_Summer11";
  std::string recoType = "pf";
  std::string jetAlgo = "akt5";
  std::string norm = "LUMI";
  if( norm!="LUMI" && norm!="SHAPE" ) {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  std::string photonID="medium";
  if( argc==3 ) {
    std::string photonID_tmp(argv[2]);
    photonID = photonID_tmp;
  }

  std::string algoType;
  if( recoType=="calo" )
    algoType = jetAlgo;
  else
    algoType = recoType+jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5") algoType = "jptak5";
  



  DrawBase* db = new DrawBase("QGStudies", recoType, jetAlgo);

  char dataFileName[150];
  if( photonID=="medium" )
    sprintf( dataFileName, "QGStudies_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  else
    sprintf( dataFileName, "QGStudies_%s_%s_%s.root", data_dataset.c_str(), algoType.c_str(), photonID.c_str());
  TFile* dataFile = TFile::Open(dataFileName);

  db->add_dataFile( dataFile, data_dataset );

  char mcPhotonJetFileName[150];
  if( photonID=="medium" )
    sprintf( mcPhotonJetFileName, "QGStudies_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
  else
    sprintf( mcPhotonJetFileName, "QGStudies_%s_%s_%s.root", mc_photonjet.c_str(), algoType.c_str(), photonID.c_str());
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName);

  db->add_mcFile( mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);


  char mcQCDFile_char[150];
  if( photonID=="medium" )
    sprintf( mcQCDFile_char, "QGStudies_%s_%s.root", mc_QCD.c_str(), algoType.c_str());
  else
    sprintf( mcQCDFile_char, "QGStudies_%s_%s_%s.root", mc_QCD.c_str(), algoType.c_str(), photonID.c_str());

  if( mc_QCD!=mc_photonjet ) {
    TFile* mcQCDFile = TFile::Open(mcQCDFile_char);
    db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);
  }



  if( norm=="LUMI" )
    db->set_lumiNormalization(191.);
  else
    db->set_shapeNormalization();

  if( data_dataset=="DATA_Run2011A_1fb" )
    db->set_lumi(1000.);
  if( data_dataset=="Photon_Run2011A-May10ReReco-v1" )
    db->set_lumi(191.);

  db->set_flags(photonID);

  db->set_outputdir();


  bool log = true;

  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "GeV/c", "Events", log);
  db->drawHisto( "nvertex_PUW", "Number of Reconstructed Vertexes", "GeV/c", "Events", log);

//db->drawHisto( "ptJetReco", "Jet Transverse Momentum", "GeV/c", "Events", log);
//db->set_rebin(5);
//db->drawHisto( "ptPhot", "Photon Transverse Momentum", "GeV/c", "Events", log);
//db->drawHisto( "ptPhot_3050", "Photon Transverse Momentum", "GeV/c", "Events", log);
//db->drawHisto( "ptPhot_5080", "Photon Transverse Momentum", "GeV/c", "Events", log);
//db->drawHisto( "ptPhot_80120", "Photon Transverse Momentum", "GeV/c", "Events", log);
//db->set_rebin(1);

//db->set_legendTitle("30 < p_{T} < 50 GeV");
//db->drawHisto( "nChargedJetReco_3050", "Jet Charged Multiplicity", "", "Events", log);
//db->drawHisto( "nNeutralJetReco_3050", "Jet Neutral Multiplicity", "", "Events", log);
//db->drawHisto( "ptDJetReco_3050", "Jet p_{T} Distribution", "", "Events", log);
//db->set_rebin(5);
//db->drawHisto( "QGLikelihoodJetReco_3050", "Jet Q-G Likelihood", "", "Events", log);
//db->set_rebin(1);

//db->drawHisto( "nChargedJetReco_antibtag_3050", "Jet Charged Multiplicity", "", "Events", log);
//db->drawHisto( "nNeutralJetReco_antibtag_3050", "Jet Neutral Multiplicity", "", "Events", log);
//db->drawHisto( "ptDJetReco_antibtag_3050", "Jet p_{T} Distribution", "", "Events", log);
//db->set_rebin(5);
//db->drawHisto( "QGLikelihoodJetReco_antibtag_3050", "Jet Q-G Likelihood", "", "Events", log);
//db->set_rebin(1);

//db->set_legendTitle("50 < p_{T} < 80 GeV");
//db->drawHisto( "nChargedJetReco_5080", "Jet Charged Multiplicity", "", "Events", log);
//db->drawHisto( "nNeutralJetReco_5080", "Jet Neutral Multiplicity", "", "Events", log);
//db->drawHisto( "ptDJetReco_5080", "Jet p_{T} Distribution", "", "Events", log);
//db->set_rebin(2);
//db->drawHisto( "QGLikelihoodJetReco_5080", "Jet Q-G Likelihood", "", "Events", log);
//db->set_rebin(1);

//db->drawHisto( "nChargedJetReco_antibtag_5080", "Jet Charged Multiplicity", "", "Events", log);
//db->drawHisto( "nNeutralJetReco_antibtag_5080", "Jet Neutral Multiplicity", "", "Events", log);
//db->drawHisto( "ptDJetReco_antibtag_5080", "Jet p_{T} Distribution", "", "Events", log);
//db->set_rebin(2);
//db->drawHisto( "QGLikelihoodJetReco_antibtag_5080", "Jet Q-G Likelihood", "", "Events", log);
//db->set_rebin(1);

//db->set_legendTitle("80 < p_{T} < 120 GeV");
//db->drawHisto( "nChargedJetReco_80120", "Jet Charged Multiplicity", "", "Events", log);
//db->drawHisto( "nNeutralJetReco_80120", "Jet Neutral Multiplicity", "", "Events", log);
//db->drawHisto( "ptDJetReco_80120", "Jet p_{T} Distribution", "", "Events", log);
//db->set_rebin(2);
//db->drawHisto( "QGLikelihoodJetReco_80120", "Jet Q-G Likelihood", "", "Events", log);
//db->set_rebin(1);

//db->drawHisto( "nChargedJetReco_antibtag_80120", "Jet Charged Multiplicity", "", "Events", log);
//db->drawHisto( "nNeutralJetReco_antibtag_80120", "Jet Neutral Multiplicity", "", "Events", log);
//db->drawHisto( "ptDJetReco_antibtag_80120", "Jet p_{T} Distribution", "", "Events", log);
//db->set_rebin(2);
//db->drawHisto( "QGLikelihoodJetReco_antibtag_80120", "Jet Q-G Likelihood", "", "Events", log);
//db->set_rebin(1);

  std::pair<TH1D,TH1D> h1pair_3050 = drawVariable_BGsubtr( "QGLikelihoodJetReco_antibtag", 30, 50, db );
  std::pair<TH1D,TH1D> h1pair_5080 = drawVariable_BGsubtr( "QGLikelihoodJetReco_antibtag", 50, 80, db );
  std::pair<TH1D,TH1D> h1pair_80120 = drawVariable_BGsubtr( "QGLikelihoodJetReco_antibtag", 80, 120, db );

//  drawSignalPtMix( h1pair_3050, h1pair_5080, h1pair_80120, 300, 0.03, 0.41, 0.56, db );


  delete db;
  db = 0;

  return 0;

}  




std::pair<TH1D, TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db ) {

  TFile* fileMC_photonjet = db->get_mcFile(0).file;
  TFile* fileMC_qcd = db->get_mcFile(1).file;

  if( fileMC_qcd==0 ) {
    std::cout << "Didn't find QCD file. Exiting." << std::endl;
    exit(199);
  }

  TFile* file_data = db->get_dataFile(0).file;

  char histoName[200];
  sprintf( histoName, "%s_%d%d", varName.c_str(), ptMin, ptMax);

  TH1D* h1_photonjet = (TH1D*)fileMC_photonjet->Get(histoName);
  TH1D* h1_qcd = (TH1D*)fileMC_qcd->Get(histoName);
  TH1D* h1_data = (TH1D*)file_data->Get(histoName);

  int nBins = h1_qcd->GetXaxis()->GetNbins();


  char histoName_quark[200];
  sprintf( histoName_quark, "%s_quark_noPhotID_%d%d", varName.c_str(), ptMin, ptMax);
  char histoName_gluon[200];
  sprintf( histoName_gluon, "%s_gluon_noPhotID_%d%d", varName.c_str(), ptMin, ptMax);

  //TH1D* h1_quark = (TH1D*)fileMC_qcd->Get(histoName_quark);
  //TH1D* h1_gluon = (TH1D*)fileMC_qcd->Get(histoName_gluon);
  TH1D* h1_quark = (TH1D*)fileMC_photonjet->Get(histoName_quark);
  TH1D* h1_gluon = (TH1D*)fileMC_photonjet->Get(histoName_gluon);

  // same area:
  float quark_integral = h1_quark->Integral(1, nBins);
  float gluon_integral = h1_gluon->Integral(1, nBins);
  h1_quark->Scale(1./quark_integral);
  h1_gluon->Scale(1./gluon_integral);


  char histoName_quarkFraction[200];
  sprintf( histoName_quarkFraction, "quarkFraction_%d%d", ptMin, ptMax );
  //sprintf( histoName_quarkFraction, "quarkFraction_antibtag_%d%d", ptMin, ptMax );

  TH1D* h1_quarkFraction_qcd = (TH1D*)fileMC_qcd->Get(histoName_quarkFraction);
  TH1D* h1_quarkFraction_photonjet = (TH1D*)fileMC_photonjet->Get(histoName_quarkFraction);

  std::cout <<  h1_quarkFraction_qcd << std::endl;
  std::cout <<  h1_quarkFraction_photonjet << std::endl;

  float quarkFraction_qcd = h1_quarkFraction_qcd->GetBinContent(1);
  float quarkFraction_photonjet = h1_quarkFraction_photonjet->GetBinContent(1);

  float qcd_integral = h1_qcd->Integral(1, nBins);
  float photonjet_integral = h1_photonjet->Integral(1, nBins);
  float data_integral = h1_data->Integral(1, nBins);

  TH1D* h1_qcd_fromFractions = new TH1D(*h1_qcd);

  for( unsigned iBin=1; iBin<nBins+1; ++iBin ) {

    float thisBinValue = quarkFraction_qcd*h1_quark->GetBinContent(iBin) + (1.-quarkFraction_qcd)*h1_gluon->GetBinContent(iBin);
    h1_qcd_fromFractions->SetBinContent( iBin, thisBinValue );

  } //for bins


  float qcd_fromFraction_integral = h1_qcd_fromFractions->Integral(1, nBins);
  h1_qcd_fromFractions->Scale( 1./qcd_fromFraction_integral );
  h1_photonjet->Scale( 1./photonjet_integral );

  float ratio = photonjet_integral/qcd_integral;
  h1_qcd_fromFractions->Scale( data_integral/(ratio + 1.) );
  h1_photonjet->Scale( data_integral*ratio/(ratio + 1.) );

  h1_photonjet->SetFillColor( 46 );
  h1_qcd_fromFractions->SetFillColor( 38 );
  h1_qcd->SetFillColor( kRed );
  h1_qcd->SetLineColor( kRed );
  h1_qcd->SetFillStyle( 3004 );

  h1_data->SetMarkerStyle(20);
  h1_data->SetMarkerSize(1.2);

  int rebin = 2;
  if( ptMin==30 ) rebin = 5;

  h1_quark->Rebin(rebin);
  h1_gluon->Rebin(rebin);
  h1_photonjet->Rebin(rebin);
  h1_qcd->Rebin(rebin);
  h1_qcd_fromFractions->Rebin(rebin);
  h1_data->Rebin(rebin);

  nBins = h1_quark->GetXaxis()->GetNbins();

  THStack* mcstack = new THStack();
  mcstack->Add( h1_qcd_fromFractions );
  mcstack->Add( h1_photonjet );

  float data_max = h1_data->GetMaximum();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();


  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., 1.5*data_max);
  h2_axes->SetXTitle("Q-G Likelihood");
  h2_axes->SetYTitle("Events");


  TLegend* legend_qcd = new TLegend( 0.2, 0.75, 0.5, 0.9 );
  legend_qcd->SetFillColor(0);
  legend_qcd->SetTextSize(0.04);
  legend_qcd->AddEntry( h1_qcd, "QCD MC", "F" );
  legend_qcd->AddEntry( h1_qcd_fromFractions, "QCD From Comp", "F" );

  
  //h2_axes->Draw(); 

  h1_qcd_fromFractions->DrawNormalized("histo");
  h1_qcd->DrawNormalized("histo same");

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  legend_qcd->Draw("same");

  gPad->RedrawAxis();

  char canvasName_qcd[500];
  sprintf( canvasName_qcd, "%s/compare_qcd_%s_%d%d.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );
  
  c1->SaveAs( canvasName_qcd );
  
  c1->Clear();

  char legendTitle[200];
  sprintf( legendTitle, "%d < p_{T} < %d GeV", ptMin , ptMax );
  TLegend* legend = new TLegend( 0.2, 0.65, 0.5, 0.9, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( h1_data, "Data", "P" );
  legend->AddEntry( h1_photonjet, "#gamma+Jet MC", "F" );
  legend->AddEntry( h1_qcd_fromFractions, "QCD MC", "F" );
  

  h2_axes->Draw(); 

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  legend->Draw("same");

  mcstack->Draw("histo same");

  h1_data->Draw("E same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/%s_%d%d_fromComponents.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );

  c1->SaveAs(canvasName);


  // and now subtract gluon component both from data and MC


  TH1D* h1_photonjet_quarkOnly = new TH1D(*h1_quark);
  h1_photonjet_quarkOnly->SetName("photonjet_quarkOnly");
  float photonjet_integral_quarkOnly = h1_photonjet_quarkOnly->Integral(1, nBins);
  h1_photonjet_quarkOnly->Scale( 1./photonjet_integral_quarkOnly );
  h1_photonjet_quarkOnly->Scale( quarkFraction_photonjet*h1_photonjet->Integral(1,nBins) );

  TH1D* h1_qcd_quarkOnly = new TH1D(*h1_quark);
  h1_qcd_quarkOnly->SetName("qcd_quarkOnly");
  float qcd_integral_quarkOnly = h1_qcd_quarkOnly->Integral(1, nBins);
  h1_qcd_quarkOnly->Scale( 1./qcd_integral_quarkOnly );
  h1_qcd_quarkOnly->Scale( quarkFraction_qcd*h1_qcd_fromFractions->Integral(1,nBins) );

  TH1D* h1_photonjet_gluonOnly = new TH1D(*h1_gluon);
  h1_photonjet_gluonOnly->SetName("photonjet_gluonOnly");
  float photonjet_integral_gluonOnly = h1_photonjet_gluonOnly->Integral(1, nBins);
  h1_photonjet_gluonOnly->Scale( 1./photonjet_integral_gluonOnly );
  h1_photonjet_gluonOnly->Scale( (1.-quarkFraction_photonjet)*h1_photonjet->Integral(1,nBins) );

  TH1D* h1_qcd_gluonOnly = new TH1D(*h1_gluon);
  h1_qcd_gluonOnly->SetName("qcd_gluonOnly");
  float qcd_integral_gluonOnly = h1_qcd_gluonOnly->Integral(1, nBins);
  h1_qcd_gluonOnly->Scale( 1./qcd_integral_gluonOnly );
  h1_qcd_gluonOnly->Scale( (1.-quarkFraction_qcd)*h1_qcd_fromFractions->Integral(1,nBins) );

//TH1D* h1_qcd_quark = new TH1D(*h1_qcd);
//h1_qcd_quark->SetName("qcd_quark");
//h1_qcd_quark->Add( h1_quark, -quarkFraction_qcd );
////h1_qcd_quark->Add( h1_quark, -quarkFraction_qcd*(h1_qcd_quark->Integral(1, nBins)) );
//TH1D* h1_qcd_gluon = new TH1D(*h1_qcd);
//h1_qcd_gluon->SetName("qcd_gluon");
//h1_qcd_gluon->Add( h1_gluon, -(1.-quarkFraction_qcd) );
////h1_qcd_gluon->Add( h1_gluon, -(1.-quarkFraction_qcd)*(h1_qcd_gluon->Integral(1, nBins)) );


  TH1D* h1_mc_gluon_only = new TH1D(*h1_photonjet_gluonOnly);
  h1_mc_gluon_only->SetName("mc_gluon_only");
  h1_mc_gluon_only->Add( h1_qcd_gluonOnly );

  h1_data->SetName("data");


  // subtract from mc:
//  TH1D* h1_mc_quark_only = new TH1D("mc_quark_only", "", nBins, 0., 1.0001);
  TH1D* h1_mc_quark_only = new TH1D(*h1_photonjet);
  h1_mc_quark_only->Add(h1_qcd_fromFractions);
//  h1_mc_quark_only = (mcstack->GetHistogram());
  h1_mc_quark_only->SetName("mc_quark_only");
  h1_mc_quark_only->Add( h1_mc_gluon_only, -1. );

  // subtract from data:
  TH1D* h1_data_quark_only = new TH1D(*h1_data);
  h1_data_quark_only->SetName("data_quark_only");
  h1_data_quark_only->Add( h1_mc_gluon_only, -1. );


  //h1_mc_quark_only->Scale( h1_data_quark_only->Integral(1,nBins)/h1_mc_quark_only->Integral(1,nBins) );

TFile* file_prova = TFile::Open("prova.root", "recreate");
file_prova->cd();
h1_qcd_quarkOnly->Write();
h1_qcd_gluonOnly->Write();
h1_photonjet_quarkOnly->Write();
h1_photonjet_gluonOnly->Write();
  h1_mc_gluon_only->Write();
  h1_mc_quark_only->Write();
h1_data->Write();
h1_data_quark_only->Write();
file_prova->Close();
  h1_mc_quark_only->SetFillColor(46);

  h1_data_quark_only->SetMarkerStyle(20);
  h1_data_quark_only->SetMarkerSize(1.1);


  c1->Clear();

  h2_axes->Draw();

  label_cms->Draw("same");
  label_sqrt->Draw("same");

  TLegend* legend_quark = new TLegend(0.2, 0.65, 0.5, 0.9, legendTitle );
  legend_quark->SetFillColor(0);
  legend_quark->SetTextSize(0.04);
  legend_quark->AddEntry( h1_data_quark_only, "Data (gluon subtracted)", "P" );
  legend_quark->AddEntry( h1_mc_quark_only, "MC (gluon subtracted)", "F" );
  legend_quark->Draw("same");

  h1_mc_quark_only->Draw("histo same"); 
  h1_data_quark_only->Draw("E same"); 

  gPad->RedrawAxis();

  char canvasName_quark[500];
  sprintf( canvasName_quark, "%s/%s_%d%d_quarkOnly.eps", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax );

  c1->SaveAs( canvasName_quark );

  
  //check efficiencies:

  int iBin_cut = h1_data_quark_only->FindBin(0.2);

  float eff_mc   = h1_mc_quark_only->Integral(iBin_cut+1, nBins)/h1_mc_quark_only->Integral(1, nBins);
  float eff_data = h1_data_quark_only->Integral(iBin_cut+1, nBins)/h1_data_quark_only->Integral(1, nBins);

  float effErr_mc   = sqrt( eff_mc*(1.-eff_mc)/h1_mc_quark_only->GetEntries());
  float effErr_data = sqrt( eff_data*(1.-eff_data)/h1_data_quark_only->GetEntries());

  std::cout << std::endl << "*** pt bin: " << ptMin << "-" << ptMax << std::endl;
  std::cout << "requiring QGLikelihood>0.2" << std::endl;
  std::cout << "eff(mc): " << eff_mc << " +/- " << effErr_mc << std::endl;
  std::cout << "eff(data): " << eff_data << " +/- " << effErr_data << std::endl;
  std::cout << std::endl;

  std::pair< TH1D, TH1D > returnPair;
  returnPair.first = *h1_mc_quark_only;
  returnPair.second = *h1_data_quark_only;

  return returnPair;

}



//void drawSignalPtMix( std::pair<TH1D*,TH1D*> h1pair_3050, std::pair<TH1D*,TH1D*> h1pair_5080, std::pair<TH1D*,TH1D*> h1pair_80120, int mass, float frac_3050, float frac_5080, float frac_80120, DrawBase* db ) {
//
//
//
//
//}
