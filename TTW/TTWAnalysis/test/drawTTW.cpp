#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawTTW [(string)selType] [bTaggerType=\"TCHE\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);

  std::string bTaggerType = "TCHE";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }




  DrawBase* db = new DrawBase("TTW");


  std::string outputdir_str = "TTWPlots_MConly_" + selType + "_" + bTaggerType;
  db->set_outputdir(outputdir_str);


  std::string mcTTWFileName = "TTW_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", kRed+3, 3005);


  //std::string mcTTbarFileName = "TTW_TTTo2L2Nu2B_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1";
  std::string mcTTbarFileName = "TTW_TTJets_TuneZ2_7TeV-madgraph-tauola_Fall11";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TT", "t#bar{t}", 30, 3001);






  db->set_shapeNormalization();




  bool log = true;



  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 20 GeV)", "", "Events", log);

  db->set_rebin(10);
  db->drawHisto("pfMet", "pfME_{T}", "GeV", "Events", log);
  db->set_rebin();
  db->drawHisto("metSignificance", "pfME_{T} Significance", "GeV", "Events", log);


  db->set_rebin(10);
  db->drawHisto("mbb", "b-b Invariant Mass", "GeV", "Jet Pairs", log);
  db->drawHisto("mqq", "Jet-Jet Invariant Mass", "GeV", "Jet Pairs", log);

  db->drawHisto("ptbb", "b-b Trasverse Momentum", "GeV", "Jet Pairs", log);
  db->drawHisto("ptqq", "Jet-Jet Trasverse Momentum", "GeV", "Jet Pairs", log);

  db->set_rebin(20);
  db->drawHisto("deltaRqq", "#DeltaR Between b-Jets", "" );
  db->drawHisto("deltaRbb", "#DeltaR Between Jets", "" );

  db->drawHisto("deltaR_b_lept_min", "#DeltaR Between b and Closest Lepton", "" );
  db->drawHisto("deltaR_b_lept_max", "#DeltaR Between b and Farthest Lepton", "" );


  db->set_rebin(1);
  db->drawHisto("deltaRll", "#DeltaR Between Leptons", "", "Lepton Pairs");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaLept1", "Lead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaLept2", "Sublead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaMu", "Muon Pseudorapidity", "", "Events");
  db->drawHisto("etaEle", "Electron Pseudorapidity", "", "Events");

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

  db->drawHisto("QGLikelihoodJet1", "Lead b-Jet QG LD", "", "Events", log);
  db->drawHisto("QGLikelihoodJet2", "Sublead b-Jet QG LD", "", "Events", log);
  db->drawHisto("QGLikelihoodJet3", "Lead q-Jet QG LD", "", "Events", log);
  db->drawHisto("QGLikelihoodJet4", "Sublead q-Jet QG LD", "", "Events", log);

  db->set_rebin();
  db->drawHisto("partFlavorJet1", "Lead b-Jet PDG ID", "", "Events");
  db->drawHisto("partFlavorJet2", "Sublead b-Jet PDG ID", "", "Events");
  db->drawHisto("partFlavorJet3", "Lead q-Jet PDG ID", "", "Events");
  db->drawHisto("partFlavorJet4", "Sublead q-Jet PDG ID", "", "Events");




  std::string mcZJetsFileName = "TTW_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 39, 3001);


  std::string mcWZFileName = "TTW_VV_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1";
  mcWZFileName += "_" + selType;
  mcWZFileName += ".root";
  TFile* mcWZFile = TFile::Open(mcWZFileName.c_str());
  db->add_mcFile( mcWZFile, "WZtoAnything_TuneZ2", "WW / WZ / ZZ", 38, 3004);


  std::string TTZFileName = "TTW_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  TTZFileName += "_" + selType;
  TTZFileName += ".root";
  TFile* TTZFile = TFile::Open(TTZFileName.c_str());
  db->add_mcFile( TTZFile, "ttZ", "t#bar{t} + Z", 41, 3002);

  std::string singleTopFileName = "TTW_SingleTop_Summer11";
  singleTopFileName += "_" + selType;
  singleTopFileName += ".root";
  TFile* singleTopFile = TFile::Open(singleTopFileName.c_str());
  db->add_mcFile( singleTopFile, "Singletop", "Single Top", 42, 3003);

  db->set_lumiNormalization(5000.);
  db->set_noStack(false);
  db->set_yAxisMaxScale(2.3);
  db->set_getBinLabels(true);
  db->drawHisto("leptType", "", "", "Events", log);


  std::string yieldsFileName = "yields_"+selType+".txt";
  ofstream yieldsFile(yieldsFileName.c_str());

  yieldsFile << "------------------------" << std::endl;
  yieldsFile << "Expected Yields @ 5 fb-1" << std::endl;
  yieldsFile << "------------------------" << std::endl;
  yieldsFile << "Dataset \t\t& #mu#mu \t& ee \t\t& e#mu " << std::endl;

  float s_mumu = 0.;
  float s_ee = 0.;
  float s_emu = 0.;

  float b_mumu  = 0.;
  float b_ee  = 0.;
  float b_emu  = 0.;

  for( unsigned i=0; i<db->get_lastHistos_mc().size(); ++i )  {

    float mumu = db->get_lastHistos_mc()[i]->GetBinContent(1);
    float ee = db->get_lastHistos_mc()[i]->GetBinContent(2);
    float emu = db->get_lastHistos_mc()[i]->GetBinContent(3);

    if(  db->get_mcFiles()[i].legendName=="t#bar{t} + W" ) {

      s_mumu += mumu;
      s_ee += ee;
      s_emu += emu;

    } else {

      b_mumu += mumu;
      b_ee += ee;
      b_emu += emu;

    }

    yieldsFile << db->get_mcFiles()[i].legendName.c_str();
    if( db->get_mcFiles()[i].legendName.size()<12 ) yieldsFile << "\t";
    yieldsFile << Form("\t& %.5f \t& %.5f \t& %.5f \\\\", mumu, ee, emu)<< std::endl;

  }
    
  yieldsFile << "Total Background\t& " << b_mumu << "\t& " << b_ee << "\t& " << b_emu << "\\\\" << std::endl;
  yieldsFile << "s/sqrt(b)    \t& " << s_mumu/sqrt(b_mumu) << "\t& " << s_ee/sqrt(b_ee) << "\t& " << s_emu/sqrt(b_emu) << "\\\\" << std::endl;
  yieldsFile.close();


  delete db;
  db = 0;

  return 0;

}  


