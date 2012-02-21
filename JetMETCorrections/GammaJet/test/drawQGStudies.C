#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"


bool useMCassoc_ = false;
bool ONEVTX = false;



std::pair<TH1D,TH1D> drawVariable_BGsubtr( const std::string& varName, int ptMin, int ptMax, DrawBase* db );
std::pair<TH1D, TH1D> drawVariable_BGsubtr_fromTree( const std::string& varName, const std::string& triggerCondition, int ptMin, int ptMax, DrawBase* db );
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
  //std::string mc_QCD = "QCD_EMEnriched_Spring11";
  std::string mc_QCD = "QCD_EMEnriched_Summer11";
  //std::string mc_QCD = "QCD_Summer11";
  std::string recoType = "pf";
  std::string jetAlgo = "akt5";
  //std::string norm = "LUMI";
  std::string norm = "SHAPE";
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
  

  TString data_dataset_tstr(data_dataset);
  std::string puReweighing;
  if( data_dataset_tstr.Contains("Run2011A") )
    puReweighing = "Run2011A";
  else if( data_dataset_tstr.Contains("Run2011B") )
    puReweighing = "Run2011B";
  else if( data_dataset_tstr.Contains("Run2011_FULL") )
    //puReweighing = "Run2011A";
    puReweighing = "Run2011FULL";


  DrawBase* db = new DrawBase("QGStudies", recoType, jetAlgo);
  db->set_isCMSArticle(true);
  db->set_lumiOnRightSide(true);

  std::string dataFileName = "QGStudies_" + data_dataset + "_" + algoType;
  dataFileName += ".root";
  //if( photonID=="medium" )
  //  sprintf( dataFileName, "QGStudies_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  //else
  //  sprintf( dataFileName, "QGStudies_%s_%s_%s.root", data_dataset.c_str(), algoType.c_str(), photonID.c_str());
  TFile* dataFile = TFile::Open(dataFileName.c_str());

  db->add_dataFile( dataFile, data_dataset );

  //char mcPhotonJetFileName[150];
  //if( photonID=="medium" )
  //  sprintf( mcPhotonJetFileName, "QGStudies_%s_%s.root", mc_photonjet.c_str(), algoType.c_str());
  //else
  //  sprintf( mcPhotonJetFileName, "QGStudies_%s_%s_%s.root", mc_photonjet.c_str(), algoType.c_str(), photonID.c_str());
  std::string mcPhotonJetFileName = "QGStudies_" + mc_photonjet + "_" + algoType;
  mcPhotonJetFileName += "_PU" + puReweighing;
  mcPhotonJetFileName += ".root";
  TFile* mcPhotonJetFile = TFile::Open(mcPhotonJetFileName.c_str());

  db->add_mcFile( mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);


//char mcQCDFile_char[150];
//if( photonID=="medium" )
//  sprintf( mcQCDFile_char, "QGStudies_%s_%s.root", mc_QCD.c_str(), algoType.c_str());
//else
//  sprintf( mcQCDFile_char, "QGStudies_%s_%s_%s.root", mc_QCD.c_str(), algoType.c_str(), photonID.c_str());

  std::string mcQCDFileName = "QGStudies_" + mc_QCD + "_" + algoType;
  mcQCDFileName += "_PU" + puReweighing;
  mcQCDFileName += ".root";
  if( mc_QCD!=mc_photonjet ) {
    TFile* mcQCDFile = TFile::Open(mcQCDFileName.c_str());
    db->add_mcFile( mcQCDFile, mc_QCD, "QCD MC", 38);
  }



  if( data_dataset=="DATA_Run2011A_1fb" )
    db->set_lumi(1000.);
  if( data_dataset=="Photon_Run2011A-May10ReReco-v1" )
    db->set_lumi(191.);
  if( data_dataset=="Photon_Run2011B-PromptReco-v1" )
    db->set_lumi(2500.);
  if( data_dataset=="Photon_Run2011_FULL")
    db->set_lumi(4600.);
  if( data_dataset=="Photon_Run2011A_FULL")
    db->set_lumi(2100.);

  if( norm=="LUMI" )
    db->set_lumiNormalization(191.);
  else
    db->set_shapeNormalization();


  db->set_flags(photonID);

  db->set_outputdir();


  bool log = true;

  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "GeV/c", "Events", log);
  db->drawHisto( "nvertex_PUW", "Number of Reconstructed Vertexes", "GeV/c", "Events", log);

//db->drawHisto( "ptJetReco", "Jet Transverse Momentum", "GeV/c", "Events", log);
  db->set_rebin(5);
//db->drawHisto( "ptPhot", "Photon Transverse Momentum", "GeV/c", "Events", log);
  db->drawHisto( "ptPhot_3050", "Photon Transverse Momentum", "GeV", "Events", log);
  db->drawHisto( "ptPhot_5080", "Photon Transverse Momentum", "GeV", "Events", log);
  db->drawHisto( "ptPhot_80120", "Photon Transverse Momentum", "GeV", "Events", log);
  db->drawHisto( "ptJet_3050", "Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto( "ptJet_5080", "Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto( "ptJet_80120", "Jet Transverse Momentum", "GeV", "Events", log);

  db->set_rebin(2);
  db->drawHisto( "QGLikelihoodJetReco_antibtag_5080", "Quark-Gluon LD");
  db->drawHisto( "QGLikelihoodJetReco_antibtag_80120", "Quark-Gluon LD");

//db->drawHisto_fromTree("tree_passedEvents", "ptPhot", "eventWeight*((run<5 || passed_Photon90_CaloIdVL_IsoL) && ptPhot>100. && !btagged && (ptJet1<10. || ptJet1<0.2*ptPhot) && passedID_no2ndJet )", 100, 50., 550., "ptPhot_HLT", "Photon p_{T}", "GeV");
//db->drawHisto_fromTree("tree_passedEvents", "ptJet0", "eventWeight*((run<5 || passed_Photon75_CaloIdVL || passed_Photon75_CaloIdVL_IsoL) && !btagged && (ptJet1<10. || ptJet1<0.2*ptPhot) && passedID_no2ndJet && ptPhot>80. && ptPhot<100. )", 100, 0., 200., "ptJet_HLT_ptPhot80100", "Jet p_{T}", "GeV", "Events", true);
//db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet0", "eventWeight*(( run<5 ||  passed_Photon75_CaloIdVL || passed_Photon75_CaloIdVL_IsoL ) && ptPhot>80. && ptPhot < 100. && !btagged && passedPhotonID)", 100, 0., 1.001, "QGLikelihood_HLT_ptPhot80100", "Q-G Discriminant", "");
//db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet0", "eventWeight*(ptPhot>80. && ptPhot < 100. && !btagged && (ptJet1<10. || ptJet1<0.2*ptPhot) && passedID_no2ndJet )", 100, 0., 1.001, "QGLikelihood_ptPhot80100", "Q-G Discriminant", "");
//db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet0", "eventWeight*((run<5 || passed_Photon90_CaloIdVL_IsoL) && ptPhot>80. && ptPhot < 100. && !btagged && (ptJet1<10. || ptJet1<0.2*ptPhot) && passedID_no2ndJet )", 100, 0., 1.001, "QGLikelihood_HLT_ptPhot80100", "Q-G Discriminant", "");
//db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet0", "eventWeight*((run<5 || passed_Photon90_CaloIdVL_IsoL) && ptPhot>100. && !btagged && (ptJet1<10. || ptJet1<0.2*ptPhot) && passedID_no2ndJet )", 100, 0., 1.001, "QGLikelihood_HLT_ptPhot100", "Q-G Discriminant", "");


  db->set_rebin();
  db->drawHisto_fromTree("tree_passedEvents", "nvertex", "eventWeight*(run<5 || passed_Photon75_CaloIdVL || passed_Photon75_CaloIdVL_IsoL)", 25, 0.5, 25.5, "nvertex_PU_HLT", "Number of Reconstructed Vertixes", "", "Events", true);
  std::pair<TH1D,TH1D> h1pair_80120 = drawVariable_BGsubtr_fromTree( "QGLikelihoodJet0", "passed_Photon75_CaloIdVL || passed_Photon75_CaloIdVL_IsoL", 80, 100, db );

//std::pair<TH1D,TH1D> h1pair_3050 = drawVariable_BGsubtr( "QGLikelihoodJetReco_antibtag", 30, 50, db );
//std::pair<TH1D,TH1D> h1pair_5080 = drawVariable_BGsubtr( "QGLikelihoodJetReco_antibtag", 50, 80, db );
//std::pair<TH1D,TH1D> h1pair_80120 = drawVariable_BGsubtr( "QGLikelihoodJetReco_antibtag", 80, 120, db );

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


  // get template shapes for quark and gluons from events with no photonID (maximise stat)
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
  h2_axes->SetXTitle("Quark-Gluon LD");
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





std::pair<TH1D, TH1D> drawVariable_BGsubtr_fromTree( const std::string& varName, const std::string& triggerCondition, int ptMin, int ptMax, DrawBase* db ) {

  TFile* fileMC_photonjet = db->get_mcFile(0).file;
  TFile* fileMC_qcd = db->get_mcFile(1).file;

  if( fileMC_qcd==0 ) {
    std::cout << "Didn't find QCD file. Exiting." << std::endl;
    exit(199);
  }

  TTree* tree_photonjet = (TTree*)fileMC_photonjet->Get("tree_passedEvents");
  TTree* tree_qcd = (TTree*)fileMC_qcd->Get("tree_passedEvents");


  TFile* file_data = db->get_dataFile(0).file;
  TTree* tree_data = (TTree*)file_data->Get("tree_passedEvents");

  char histoName[200];
  sprintf( histoName, "%s_%d%d", varName.c_str(), ptMin, ptMax);

  int nBins = 50;

  TH1D* h1_photonjet = new TH1D("photonjet", "", nBins, 0., 1.0001);
  TH1D* h1_qcd = new TH1D("qcd", "", nBins, 0., 1.0001);
  TH1D* h1_data = new TH1D("data", "", nBins, 0., 1.0001);

  float ptJetMin = 80.;
  float ptJetMax = 120.;

  float secondJetThreshold = 0.2;

  char cuts[600];
  sprintf( cuts, "eventWeight*(      ptJet0>%f && ptJet0<%f && !btagged && ptPhot>%d. && ptPhot<%d. && (ptJet1<10. || ptJet1<%f*ptPhot) && passedID_no2ndJet )", ptJetMin, ptJetMax, ptMin, ptMax, secondJetThreshold );
  char cuts_data[800];
  sprintf( cuts_data, "eventWeight*( ptJet0>%f && ptJet0<%f && !btagged && ptPhot>%d. && ptPhot<%d. && (ptJet1<10. || ptJet1<%f*ptPhot) && passedID_no2ndJet  && (%s) )", ptJetMin, ptJetMax, ptMin, ptMax, secondJetThreshold, triggerCondition.c_str() );

  tree_photonjet->Project( "photonjet", varName.c_str(), cuts );
  tree_qcd->Project(       "qcd",       varName.c_str(), cuts );
  tree_data->Project(      "data",      varName.c_str(), cuts_data );


  //int nBins = h1_qcd->GetXaxis()->GetNbins();


  // get template shapes for quark and gluons from events with no photonID (maximise stat)
  char cuts_quark[600];
  sprintf( cuts_quark, "eventWeight*( ptJet0>%f&& ptJet0<%f && !btagged && ptPhot>%d. && ptPhot<%d. && abs(pdgIdPartJet0)<6)", ptJetMin, ptJetMax, ptMin, ptMax );
  char cuts_gluon[600];
  sprintf( cuts_gluon, "eventWeight*( ptJet0>%f&& ptJet0<%f && !btagged && ptPhot>%d. && ptPhot<%d. && pdgIdPartJet0==21)", ptJetMin, ptJetMax, ptMin, ptMax );

  TH1D* h1_quark = new TH1D("quark", "", nBins, 0., 1.0001 );
  TH1D* h1_gluon = new TH1D("gluon", "", nBins, 0., 1.0001 );
  tree_photonjet->Project( "quark", varName.c_str(), cuts_quark );
  tree_photonjet->Project( "gluon", varName.c_str(), cuts_gluon );


  // same area:
  float quark_integral = h1_quark->Integral(1, nBins);
  float gluon_integral = h1_gluon->Integral(1, nBins);
  h1_quark->Scale(1./quark_integral);
  h1_gluon->Scale(1./gluon_integral);


  char cuts_quark_full[600];
  sprintf( cuts_quark_full, "eventWeight*( ptJet0>%f && ptJet0<%f && !btagged && ptPhot>%d. && ptPhot<%d. && abs(pdgIdPartJet0)<6 && (ptJet1<10. || ptJet1<%f*ptPhot) && passedID_no2ndJet )", ptJetMin, ptJetMax, ptMin, ptMax, secondJetThreshold );
  char cuts_gluon_full[600];
  sprintf( cuts_gluon_full, "eventWeight*( ptJet0>%f && ptJet0<%f && !btagged && ptPhot>%d. && ptPhot<%d. && pdgIdPartJet0==21 && (ptJet1<10. || ptJet1<%f*ptPhot) && passedID_no2ndJet )", ptJetMin, ptJetMax, ptMin, ptMax, secondJetThreshold );

  TH1D* h1_photonJet_quark_full = new TH1D("photonJet_quark_full", "", nBins, 0., 1.0001 );
  TH1D* h1_photonJet_gluon_full = new TH1D("photonJet_gluon_full", "", nBins, 0., 1.0001 );
  tree_photonjet->Project( "photonJet_quark_full", varName.c_str(), cuts_quark_full );
  tree_photonjet->Project( "photonJet_gluon_full", varName.c_str(), cuts_gluon_full );

  TH1D* h1_qcd_quark_full = new TH1D("qcd_quark_full", "", nBins, 0., 1.0001 );
  TH1D* h1_qcd_gluon_full = new TH1D("qcd_gluon_full", "", nBins, 0., 1.0001 );
  tree_qcd->Project( "qcd_quark_full", varName.c_str(), cuts_quark_full );
  tree_qcd->Project( "qcd_gluon_full", varName.c_str(), cuts_gluon_full );


  float quarkFraction_qcd = h1_qcd_quark_full->Integral() / ( h1_qcd_quark_full->Integral() + h1_qcd_gluon_full->Integral() );
  float quarkFraction_photonjet = h1_photonJet_quark_full->Integral() / ( h1_photonJet_quark_full->Integral() + h1_photonJet_gluon_full->Integral() );

  float qcd_integral = h1_qcd->Integral(1, nBins);
  float photonjet_integral = h1_photonjet->Integral(1, nBins);
  float data_integral = h1_data->Integral(1, nBins);

  float quarkFraction_total = ( qcd_integral*quarkFraction_qcd + photonjet_integral*quarkFraction_photonjet) / ( qcd_integral + photonjet_integral );

std::cout << "quarkFraction_qcd: " << quarkFraction_qcd << std::endl;
std::cout << "quarkFraction_photonjet: " << quarkFraction_photonjet << std::endl;
std::cout << "quarkFraction_total: " << quarkFraction_total << std::endl;

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

  char yAxisTitle[200];
  if( (h1_data->GetBinWidth(1)) < 0.1 )
    sprintf( yAxisTitle, "Events / (%.2f)", h1_data->GetBinWidth(1) );
  else if( ((int)(10.*h1_data->GetBinWidth(1)) % 10) == 0 )
    sprintf( yAxisTitle, "Events / (%.0f)", h1_data->GetBinWidth(1) );
  else
    sprintf( yAxisTitle, "Events / (%.1f)", h1_data->GetBinWidth(1) );

  float yMax = 1.5*data_max;

  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.0001, 10, 0., yMax);
  h2_axes->SetXTitle("Quark-Gluon LD");
  h2_axes->SetYTitle(yAxisTitle);
  if( yMax>1000. ) {
    h2_axes->GetYaxis()->SetTitleOffset(1.55);
  }
  if( yMax>10000. ) {
    h2_axes->GetYaxis()->SetTitleOffset(1.55);
    h2_axes->GetYaxis()->SetLabelSize(0.04);
  }

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
  //sprintf( legendTitle, "%d < p_{T} < %d GeV", ptMin , ptMax );
  sprintf( legendTitle, "%.0f < p_{T} < %.0f GeV", ptJetMin , ptJetMax );
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


