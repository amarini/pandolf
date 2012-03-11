#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawTTW [(string)selType] [bTaggerType=\"TCHE\"]" << std::endl;
    exit(23);
  }

  std::string leptType = "ALL";

  std::string selType(argv[1]);

  std::string bTaggerType = "TCHE";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }




  DrawBase* db = new DrawBase("TTW");


  std::string outputdir_str = "TTWPlots_MConly_" + selType + "_" + bTaggerType + "_" + leptType;
  db->set_outputdir(outputdir_str);


  std::string mcTTWFileName = "TTW_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + leptType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", kRed+3, 3005);


  std::string mcTTbarFileName = "TTW_TTTo2L2Nu2B_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TT", "t#bar{t}", 39, 3003);

  std::string mcZJetsFileName = "TTW_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 30, 3001);

  std::string mcWZFileName = "TTW_VV_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1";
  mcWZFileName += "_" + selType;
  mcWZFileName += "_" + leptType;
  mcWZFileName += ".root";
  TFile* mcWZFile = TFile::Open(mcWZFileName.c_str());
  db->add_mcFile( mcWZFile, "WZtoAnything_TuneZ2", "WW / WZ / ZZ", 38, 3004);


  std::string TTZFileName = "TTW_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  TTZFileName += "_" + selType;
  TTZFileName += "_" + leptType;
  TTZFileName += ".root";
  TFile* TTZFile = TFile::Open(TTZFileName.c_str());
  db->add_mcFile( TTZFile, "ttZ", "t#bar{t} + Z", 41, 3002);






  db->set_shapeNormalization();




  bool log = true;


  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 20 GeV)", "", "Events", log);

  db->drawHisto("pfMet", "pfME_{T}", "GeV", "Events", log);
  db->drawHisto("metSignificance", "pfME_{T} Significance", "GeV", "Events", log);


  db->set_rebin(10);
  db->drawHisto("mbb", "b-b Invariant Mass", "GeV", "Jet Pairs", log);
  db->drawHisto("mqq", "Jet-Jet Invariant Mass", "GeV", "Jet Pairs", log);

  db->drawHisto("ptbb", "b-b Trasverse Momentum", "GeV", "Jet Pairs", log);
  db->drawHisto("ptqq", "Jet-Jet Trasverse Momentum", "GeV", "Jet Pairs", log);

  db->set_rebin(20);
  db->drawHisto("deltaRqq", "#DeltaR Between b-Jets", "" );
  db->drawHisto("deltaRbb", "#DeltaR Between Jets", "" );


  db->set_rebin(1);
  db->drawHisto("deltaRll", "#DeltaR Between Leptons", "", "Lepton Pairs");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaLept1", "Lead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaLept2", "Sublead Lepton Pseudorapidity", "", "Events");

  db->set_yAxisMaxScale( 1.1 );
  db->set_rebin(5);
  db->set_xAxisMax(250.);
  db->drawHisto("ptLept1", "Lead Lepton p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(150.);
  db->drawHisto("ptLept2", "Sublead Lepton p_{T}", "GeV", "Events", log);


  db->set_xAxisMax(250.);
  db->drawHisto("ptJet1", "Lead b-Jet p_{T}", "GeV", "Events", log);
  db->drawHisto("ptJet3", "Lead q-Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(150.);
  db->drawHisto("ptJet2", "Sublead b-Jet p_{T}", "GeV", "Events", log);
  db->drawHisto("ptJet4", "Sublead q-Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax();
  db->set_yAxisMaxScale( 1.6 );
  db->set_xAxisMax();
  db->drawHisto("etaJet1", "Lead b-Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("etaJet2", "Sublead b-Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("etaJet3", "Lead q-Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("etaJet4", "Sublead q-Jet Pseudorapidity", "", "Events", log);

  db->set_rebin();
  db->drawHisto("partFlavorJet1", "Lead b-Jet PDG ID", "", "Events");
  db->drawHisto("partFlavorJet2", "Sublead b-Jet PDG ID", "", "Events");
  db->drawHisto("partFlavorJet3", "Lead q-Jet PDG ID", "", "Events");
  db->drawHisto("partFlavorJet4", "Sublead q-Jet PDG ID", "", "Events");


  db->set_lumiNormalization(10000.);
  db->set_noStack(false);
  db->drawHisto("leptType", "leptType", "");




  delete db;
  db = 0;

  return 0;

}  


