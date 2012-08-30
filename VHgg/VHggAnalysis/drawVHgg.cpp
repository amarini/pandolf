#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"
#include "cl95cms.C"

bool separate_signals = false;


void printYields( DrawBase* db, const std::string& suffix, bool doUL=false );


int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawVHgg [(string)selType] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);


  std::string bTaggerType = "JP";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }




  DrawBase* db_nostack = new DrawBase("VHgg_nostack");
  DrawBase* db_stack = new DrawBase("VHgg_stack");

  db_nostack->set_lumiOnRightSide();
  db_stack->set_lumiOnRightSide();

  db_nostack->set_shapeNormalization();

  db_stack->set_lumiNormalization(30000.);
  db_stack->set_noStack(false);

  std::string outputdir_str = "VHggPlots_MConly_" + selType + "_" + bTaggerType;
  db_nostack->set_outputdir(outputdir_str);
  db_stack->set_outputdir(outputdir_str);

  int signalFillColor = 46;


  // VH only for shape comparisons
  std::string VHFileName = "VHgg_WH_ZH_HToGG_M-125_8TeV-pythia6";
  VHFileName += "_" + selType;
  VHFileName += "_" + bTaggerType;
  VHFileName += ".root";
  TFile* VHFile = TFile::Open(VHFileName.c_str());
  db_nostack->add_mcFile( VHFile, "VH", "VH (125)", signalFillColor, 0);

  if( separate_signals ) {

    db_stack->add_mcFile( VHFile, "VH", "VH (125)", signalFillColor, 0);

    std::string GluGluHFileName = "VHgg_GluGluToHToGG_M-125_8TeV-powheg-pythia6";
    GluGluHFileName += "_" + selType;
    GluGluHFileName += "_" + bTaggerType;
    GluGluHFileName += ".root";
    TFile* GluGluHFile = TFile::Open(GluGluHFileName.c_str());
    db_stack->add_mcFile( GluGluHFile, "GluGluH", "H (125)", signalFillColor+1, 0);

    std::string VBFHFileName = "VHgg_VBF_HToGG_M-125_8TeV-powheg-pythia6";
    VBFHFileName += "_" + selType;
    VBFHFileName += "_" + bTaggerType;
    VBFHFileName += ".root";
    TFile* VBFHFile = TFile::Open(VBFHFileName.c_str());
    db_stack->add_mcFile( VBFHFile, "VBFH", "H (125)", signalFillColor+2, 0);

  } else {

    // inclusive signal file for stack plot
    std::string HToGGFileName = "VHgg_HToGG_M-125_8TeV-pythia6";
    HToGGFileName += "_" + selType;
    HToGGFileName += "_" + bTaggerType;
    HToGGFileName += ".root";
    TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());
    db_stack->add_mcFile( HToGGFile, "HToGG", "H (125)", signalFillColor, 0);

  }


  std::string DiPhotonFileName = "VHgg_DiPhoton_8TeV-pythia6";
  DiPhotonFileName += "_" + selType;
  DiPhotonFileName += "_" + bTaggerType;
  DiPhotonFileName += ".root";
  TFile* DiPhotonFile = TFile::Open(DiPhotonFileName.c_str());
  db_nostack->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);
  db_stack->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);

  std::string GammaJetFileName = "VHgg_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  GammaJetFileName += "_" + selType;
  GammaJetFileName += "_" + bTaggerType;
  GammaJetFileName += ".root";
  TFile* GammaJetFile = TFile::Open(GammaJetFileName.c_str());
  db_nostack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 38);
  db_stack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 38);


  std::string DiBosonFileName = "VHgg_VV_8TeV";
  DiBosonFileName += "_" + selType;
  DiBosonFileName += "_" + bTaggerType;
  DiBosonFileName += ".root";
  TFile* DiBosonFile = TFile::Open(DiBosonFileName.c_str());
  db_stack->add_mcFile( DiBosonFile, "DiBoson", "Diboson", 39);

  std::string TriBosonFileName = "VHgg_VGG_8TeV";
  TriBosonFileName += "_" + selType;
  TriBosonFileName += "_" + bTaggerType;
  TriBosonFileName += ".root";
  TFile* TriBosonFile = TFile::Open(TriBosonFileName.c_str());
  db_stack->add_mcFile( TriBosonFile, "TriBoson", "V + #gamma#gamma", 40);

  std::string TTFileName = "VHgg_TT_8TeV";
  TTFileName += "_" + selType;
  TTFileName += "_" + bTaggerType;
  TTFileName += ".root";
  TFile* TTFile = TFile::Open(TTFileName.c_str());
  db_stack->add_mcFile( TTFile, "TT", "Top", 44);

  std::string QCDFileName = "VHgg_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  QCDFileName += "_" + selType;
  QCDFileName += "_" + bTaggerType;
  QCDFileName += ".root";
  TFile* QCDFile = TFile::Open(QCDFileName.c_str());
  db_stack->add_mcFile( QCDFile, "QCD", "QCD", 41);







  bool log = true;


  db_nostack->drawHisto("nvertex");
  db_nostack->drawHisto("nvertex_PUW");

  db_nostack->drawHisto("njets", "Number of Jets", "", "Events");
  db_nostack->drawHisto("nbjets_loose", "Number of b-Jets (Loose)", "", "Events");
  db_nostack->drawHisto("nbjets_medium", "Number of b-Jets (Medium)", "", "Events");

  db_nostack->drawHisto("mjj", "Dijet Mass", "GeV");
  db_nostack->drawHisto("mjj_0btag", "Dijet Mass", "GeV");
  db_nostack->drawHisto("mjj_1btag", "Dijet Mass", "GeV");
  db_nostack->drawHisto("mjj_2btag", "Dijet Mass", "GeV");
  db_nostack->drawHisto("qgljet0", "Lead Jet Q-G LD");
  db_nostack->drawHisto("qgljet1", "Sublead Jet Q-G LD");

  db_nostack->drawHisto("ptphot0", "Lead Photon p_{T}", "GeV");
  db_nostack->drawHisto("ptphot1", "Sublead Photon p_{T}", "GeV");
  db_nostack->drawHisto("ptjet0", "Lead Jet p_{T}", "GeV");
  db_nostack->drawHisto("ptjet1", "Sublead Jet p_{T}", "GeV");

  db_nostack->set_rebin(2);
  db_nostack->drawHisto("deltaPhi", "#Delta#Phi(dijet-diphoton)", "rad");
  db_nostack->drawHisto("ptDijet", "Dijet p_{T}", "GeV");
  db_nostack->drawHisto("ptDiphot" "Diphoton p_{T}", "GeV");
  db_nostack->drawHisto("ptRatio", "Dijet p_{T}) / Diphoton p_{T}");
  db_nostack->drawHisto("ptDifference", "Dijet p_{T}) - Diphoton p_{T}", "GeV");

  db_nostack->drawHisto("deltaEtaJets", "Jet-Jet #Delta#eta");
  db_nostack->drawHisto("deltaFabsEtaJets", "Jet-Jet #Delta|#eta|");
  db_nostack->drawHisto("zeppen", "Zeppenfeld Variable");



  db_stack->set_rebin(5);

  db_stack->drawHisto("mgg_prepresel", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "prepresel" );
  db_stack->drawHisto("mgg_presel", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "presel");

  db_stack->set_legendTitle( "0 b-tag Category" );
  db_stack->drawHisto("mgg_0btag", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "0tag", true );
  db_stack->set_legendTitle( "1 b-tag Category" );
  db_stack->drawHisto("mgg_1btag", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "1tag", true );
  db_stack->set_legendTitle( "2 b-tag Category" );
  db_stack->drawHisto("mgg_2btag", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "2tag", true );




  return 0;

}




void printYields( DrawBase* db, const std::string& suffix, bool doUL ) {

  std::string yieldsFileName = db->get_outputdir() + "/yields_" + suffix + ".txt";
  ofstream yieldsFile(yieldsFileName.c_str());


  float xMin = 120.;
  float xMax = 130.;

  std::vector<TH1D*> histos = db->get_lastHistos_mc();

  int binXmin = histos[0]->FindBin(xMin);
  int binXmax = histos[0]->FindBin(xMax) -1;

  bool foundSignal = false;
  float totalBG = 0.;
  float signal = 0.;

  yieldsFile << std::endl << "Yields (@ 30 fb-1): " << std::endl;
  for( unsigned int ii=0; ii<histos.size(); ++ii ) {
    yieldsFile << db->get_mcFile(ii).datasetName << " " << histos[ii]->Integral(binXmin, binXmax) << std::endl;
    if( db->get_mcFile(ii).datasetName != "HToGG" ) 
      totalBG += histos[ii]->Integral(binXmin, binXmax);
    else {
      foundSignal = true;
      signal = histos[ii]->Integral(binXmin, binXmax);
    }
  }


  yieldsFile << "Total BG: " << totalBG << std::endl;

  float signal_xsec = 2.28E-03*(19.37 + 1.573 + 0.6966 + 0.3943); // ttH missing for now
  float total_signal = signal_xsec*db->get_lumi();
  float effS = signal/total_signal;
  yieldsFile << "Signal efficiency: " << effS << std::endl;

  if( !foundSignal ) std::cout << "WARNING!!! DIDN'T FIND SIGNAL HToGG!" << std::endl; 

  
  if( doUL ) {

    float ul = CLA( db->get_lumi(), 0., effS, 0., totalBG, 0. );
    yieldsFile << std::endl << "No error on BG:" << std::endl;
    yieldsFile << "UL: " << ul << std::endl;
    yieldsFile << "UL/SM: " << ul/signal_xsec << std::endl;
    float ul_bgerr = CLA( db->get_lumi(), 0., effS, 0., totalBG, 0.05*totalBG );
    yieldsFile << std::endl << "5\% error on BG:" << std::endl;
    yieldsFile << "UL: " << ul_bgerr << std::endl;
    yieldsFile << "UL/SM: " << ul_bgerr/signal_xsec << std::endl;

  }



}
