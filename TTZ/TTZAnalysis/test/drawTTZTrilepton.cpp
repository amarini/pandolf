#include <stdlib.h>
#include <iostream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawTTZTrilepton [(string)selType] [bTaggerType=\"SSVHE\"]" << std::endl;
    exit(23);
  }

  std::string leptType = "ALL";

  std::string selType(argv[1]);

  std::string PUType = "PUHR11_73pb";

  std::string bTaggerType = "SSVHE";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }




  DrawBase* db = new DrawBase("TTZTrilepton");


  std::string outputdir_str = "TTZTrileptonPlots_MConly_" + selType + "_" + bTaggerType + "_" + leptType;
  db->set_outputdir(outputdir_str);

  int signalFillColor = 42;

  std::string TTZFileName = "TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  TTZFileName += "_" + selType;
  TTZFileName += "_" + bTaggerType;
  TTZFileName += "_" + PUType;
  TTZFileName += "_" + leptType;
  TTZFileName += ".root";
  TFile* TTZFile = TFile::Open(TTZFileName.c_str());
  db->add_mcFile( TTZFile, "ttZ", "t#bar{t} + Z", signalFillColor, 3005);


  std::string mcZJetsFileName = "TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + bTaggerType;
  mcZJetsFileName += "_" + PUType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 30, 3001);

  std::string mcTTbarFileName = "TTZTrilepton_TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + bTaggerType;
  mcTTbarFileName += "_" + PUType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TT", "t#bar{t}", 39, 3003);

  std::string mcWZFileName = "TTZTrilepton_WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcWZFileName += "_" + selType;
  mcWZFileName += "_" + bTaggerType;
  mcWZFileName += "_" + PUType;
  mcWZFileName += "_" + leptType;
  mcWZFileName += ".root";
  TFile* mcWZFile = TFile::Open(mcWZFileName.c_str());
  db->add_mcFile( mcWZFile, "WZtoAnything_TuneZ2", "WZ + jets", 38, 3004);


  std::string mcTTWFileName = "TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + bTaggerType;
  mcTTWFileName += "_" + PUType;
  mcTTWFileName += "_" + leptType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", 33, 3002);




  db->set_shapeNormalization();




  bool log = true;


  db->set_legendTitle("Trilepton channel");

  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 10 GeV)", "", "Events", log);
  db->drawHisto("nvertex", "Number of Reconstructed Vertexes", "", "Events", log);

  float lumifb = 20.;

  db->set_lumiNormalization(lumifb*1000.);
  db->set_noStack(false);

  std::vector<TH1D*> lastHistosMC;
  float signalYield;
  float channelYieldGen = lumifb * 139 * 0.06 * 0.2 * 0.67;

  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>70. && mZll<110.)", 30, 0., 300., "ptZll", "p_{T} (Z)", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>70. && mZll<110. && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20. )", 30, 0., 300., "ptZll_jetpt20", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>70. && mZll<110. && bTagJetB1>0. && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", 30, 0., 300., "ptZll_jetpt20_btag", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>70. && mZll<110. && pfMet>30. && bTagJetB1>0. && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", 30, 0., 300., "ptZll_jetpt20_btag_met", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;

  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(ptZll>100. && mZll>70. && mZll<110. && pfMet>30. && bTagJetB1>0. && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", 30, 0., 300., "ptZll_jetpt20_btag_met_ptZll", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;



  delete db;
  db = 0;

  return 0;

}  


