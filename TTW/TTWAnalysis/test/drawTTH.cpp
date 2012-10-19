#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



void drawSingleVar( DrawBase* db, DrawBase* dbsig, const std::string& varname, int nbins, float xmin, float xmax, const std::string&  savename, const std::string& axisName, const std::string& units);
void drawSignalPlot( DrawBase* db, const std::string& varname, int nbins, float xmin, float xmax, const std::string& savename, const std::string& axisName, const std::string& units="", bool shapeNorm=false );


int main(int argc, char* argv[]) {

  if(  argc != 3 && argc != 4 ) {
    std::cout << "USAGE: ./drawTTW [(string)selType] [jetChoice] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);
  std::string jetChoice(argv[2]);

  std::string bTaggerType = "JP";
  if( argc>=4 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }



  DrawBase* db = new DrawBase("TTW");
  DrawBase* dbsig = new DrawBase("TTH");


  db->set_shapeNormalization();
  dbsig->set_lumiNormalization(30000.);


  std::string outputdir_str = "TTHPlots_MConly_" + selType + "_" + jetChoice + "_" + bTaggerType;
  db->set_outputdir(outputdir_str);
  dbsig->set_outputdir(outputdir_str);


  std::string mcTTHFileName = "TTW_TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1";
  mcTTHFileName += "_" + selType;
  mcTTHFileName += "_" + jetChoice;
  mcTTHFileName += "_" + bTaggerType;
  mcTTHFileName += ".root";
  TFile* mcTTHFile = TFile::Open(mcTTHFileName.c_str());
  db->add_mcFile( mcTTHFile, "TTH", "t#bar{t} + H", kBlack, 0);
  dbsig->add_mcFile( mcTTHFile, "TTH", "t#bar{t} + H", kBlack, 0);



  std::string mcTTWFileName = "TTW_TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + jetChoice;
  mcTTWFileName += "_" + bTaggerType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", 46, 3005);



  std::string mcTTZFileName = "TTW_TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1";
  mcTTZFileName += "_" + selType;
  mcTTZFileName += "_" + jetChoice;
  mcTTZFileName += "_" + bTaggerType;
  mcTTZFileName += ".root";
  TFile* mcTTZFile = TFile::Open(mcTTZFileName.c_str());
  db->add_mcFile( mcTTZFile, "TTZ", "t#bar{t} + Z", 38, 3004);



  drawSingleVar( db, dbsig, "nnbjets+nbjets", 9, -0.5, 8.5, "njets", "Jet Multiplicity (p_{T} > 20 GeV)", "");
  drawSingleVar( db, dbsig, "nnbjets", 9, -0.5, 8.5, "nnbjets", "Non b-Tagged Jet Multiplicity (p_{T} > 20 GeV)", "");
  drawSingleVar( db, dbsig, "nbjets", 9, -0.5, 8.5, "nbjets", "b-Tagged Jet Multiplicity (p_{T} > 20 GeV)", "");
  drawSingleVar( db, dbsig, "nnbjetsmed", 9, -0.5, 8.5, "nnbjetsmed", "Non b-Tagged Jet Multiplicity (p_{T} > 20 GeV)", "");
  drawSingleVar( db, dbsig, "nbjetsmed", 9, -0.5, 8.5, "nbjetsmed", "b-Tagged Jet Multiplicity (p_{T} > 20 GeV)", "");
  drawSingleVar( db, dbsig, "mjj1", 25, 0., 250., "mjj1", "Invariant Mass of First Jet Pair", "GeV");
  drawSingleVar( db, dbsig, "mjj2", 25, 0., 250., "mjj2", "Invariant Mass of Second Jet Pair", "GeV");
  drawSingleVar( db, dbsig, "pfMet", 25, 0., 250., "pfMet", "Particle Flow ME_{T}", "GeV");
  drawSingleVar( db, dbsig, "ptLept1", 25, 0., 350., "ptLept1", "Leading Lepton p_{T}", "GeV");
  drawSingleVar( db, dbsig, "ptLept2", 25, 0., 200., "ptLept2", "Subleading Lepton p_{T}", "GeV");
  drawSingleVar( db, dbsig, "ht", 25, 0., 1000., "ht", "H_{T}", "GeV");

  return 0;

}



void drawSingleVar( DrawBase* db, DrawBase* dbsig, const std::string& varname, int nbins, float xmin, float xmax, const std::string& savename, const std::string& axisName, const std::string& units ) {

  db->drawHisto_fromTree( "tree_passedEvents", varname, "", nbins, xmin, xmax, savename, axisName, units, "Events", true );
  drawSignalPlot( dbsig, varname, nbins, xmin, xmax, savename, axisName, units, true);

}  


void drawSignalPlot( DrawBase* db, const std::string& varname, int nbins, float xmin, float xmax, const std::string& savename, const std::string& axisName, const std::string& units, bool shapeNorm ) {

  if( shapeNorm )
    db->set_shapeNormalization();


  TH1F::AddDirectory(kTRUE);


  std::string axisName_units(axisName);
  if( units!="" ) axisName_units = axisName_units + " [" + units + "]";


  TFile* sigfile = db->get_mcFile(0).file;
  TTree* sigtree = (TTree*)sigfile->Get("tree_passedEvents");

 
  TH1D* h1_ttHtot = new TH1D("ttHtot", "", nbins, xmin, xmax); 
  h1_ttHtot->Sumw2();

  TH1D* h1_ttHWW = new TH1D("ttHWW", "", nbins, xmin, xmax); 
  h1_ttHWW->Sumw2();
  TH1D* h1_ttHZZ = new TH1D("ttHZZ", "", nbins, xmin, xmax); 
  h1_ttHZZ->Sumw2();
  TH1D* h1_ttHTauTau = new TH1D("ttHTauTau", "", nbins, xmin, xmax); 
  h1_ttHTauTau->Sumw2();
  TH1D* h1_ttHBB = new TH1D("ttHBB", "", nbins, xmin, xmax); 
  h1_ttHBB->Sumw2();


  h1_ttHtot->SetLineWidth(2);
  h1_ttHtot->SetLineColor(kBlack);
  h1_ttHtot->SetFillStyle(0);

  h1_ttHWW->SetLineWidth(2);
  h1_ttHWW->SetLineColor(46);
  h1_ttHWW->SetFillColor(46);
  h1_ttHWW->SetFillStyle(3004);

  h1_ttHTauTau->SetLineWidth(2);
  h1_ttHTauTau->SetLineColor(38);
  h1_ttHTauTau->SetFillColor(38);
  h1_ttHTauTau->SetFillStyle(3005);


  sigtree->Project( "ttHtot",    varname.c_str(), "30000.*eventWeight" );
  sigtree->Project( "ttHWW",     varname.c_str(), "30000.*eventWeight*is_ttHWWEvent" );
  sigtree->Project( "ttHTauTau", varname.c_str(), "30000.*eventWeight*is_ttHTauTauEvent" );
  sigtree->Project( "ttHZZ",     varname.c_str(), "30000.*eventWeight*is_ttHZZEvent" );


  float ymax = h1_ttHtot->GetMaximum()*1.3;
  if( shapeNorm )
    ymax /= h1_ttHtot->Integral();

  TH2D* h2_axes = new TH2D("axes", "", 10, xmin, xmax, 10, 0., ymax);
  h2_axes->SetXTitle(axisName_units.c_str());
  if( shapeNorm ) {
    h2_axes->SetYTitle("Normalized to Unity");
  } else {
    h2_axes->SetYTitle("Events");
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
 
  h2_axes->Draw();

  if( shapeNorm ) {
    h1_ttHTauTau->DrawNormalized("histo same");
    h1_ttHWW->DrawNormalized("histo same");
    h1_ttHtot->DrawNormalized("histo same");
  } else {
    h1_ttHTauTau->Draw("histo same");
    h1_ttHWW->Draw("histo same");
    h1_ttHtot->Draw("histo same");
  }


  TLegend* legend = new TLegend(0.61, 0.7, 0.91, 0.9);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->AddEntry( h1_ttHWW, "ttH, H #rightarrow WW", "F");
  legend->AddEntry( h1_ttHTauTau, "ttH, H #rightarrow #tau#tau", "F");
  legend->AddEntry( h1_ttHtot, "ttH, total", "F");
  legend->Draw("same");


  TPaveText* label_cms = db->get_labelTop();
  label_cms->Draw("same");

  gPad->RedrawAxis();


  std::string canvasName = db->get_outputdir() + "/ttHCompare_" + savename.c_str() + ".eps"; 

  c1->SaveAs(canvasName.c_str());
  
  delete c1;
  delete legend;
  delete h2_axes;

  delete h1_ttHtot;
  delete h1_ttHWW;
  delete h1_ttHZZ;
  delete h1_ttHTauTau;
  delete h1_ttHBB;
  

}
