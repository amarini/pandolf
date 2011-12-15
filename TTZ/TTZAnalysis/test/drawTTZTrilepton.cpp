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

  std::string bTaggerType = "SSVHE";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }




  DrawBase* db = new DrawBase("TTZTrilepton");


  std::string outputdir_str = "TTZTrileptonPlots_MConly_" + selType + "_" + bTaggerType + "_" + leptType;
  db->set_outputdir(outputdir_str);


  std::string TTZFileName = "TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  TTZFileName += "_" + selType;
  TTZFileName += "_" + bTaggerType;
  TTZFileName += "_" + leptType;
  TTZFileName += ".root";
  TFile* TTZFile = TFile::Open(TTZFileName.c_str());
  db->add_mcFile( TTZFile, "ttZ", "t#bar{t} + Z", kRed+3, 3005);


  std::string mcZJetsFileName = "TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + bTaggerType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 30, 3001);

  std::string mcWZFileName = "TTZTrilepton_WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcWZFileName += "_" + selType;
  mcWZFileName += "_" + bTaggerType;
  mcWZFileName += "_" + leptType;
  mcWZFileName += ".root";
  TFile* mcWZFile = TFile::Open(mcWZFileName.c_str());
  db->add_mcFile( mcWZFile, "WZtoAnything_TuneZ2", "WZ + jets", 38, 3004);

  std::string mcTTbarFileName = "TTZTrilepton_TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TTtW", "t#bar{t}", 39, 3003);


  std::string mcTTWFileName = "TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + leptType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", 41, 3002);




  db->set_shapeNormalization();




  bool log = true;


  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 20 GeV)", "", "Events", log);


  db->drawHisto("ptJet_all_presel", "Jet Transverse Momentum", "GeV", "Jets", log);

  db->set_rebin(10);
  db->drawHisto("mZjj_all_presel", "DiJet Invariant Mass", "GeV", "Jet Pairs", log);

  db->set_rebin(2);
  db->drawHisto("deltaRjj", "#DeltaR Between Jets (p_{T} > 30 GeV)", "", "Jet Pairs");
  db->set_rebin(1);
  db->drawHisto("deltaRjj_all_presel", "#DeltaR Between Jets (p_{T} > 30 GeV)", "", "Jet Pairs");
  db->drawHisto("deltaRll_presel", "#DeltaR Between Leptons", "", "Lepton Pairs");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaLept1_presel", "Lead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaLept2_presel", "Sublead Lepton Pseudorapidity", "", "Events");
  db->drawHisto("etaJet_all_presel", "Jet Pseudorapidity", "", "Jets");

  db->set_yAxisMaxScale( 1.1 );
  db->set_rebin(5);
  db->set_xAxisMax(250.);
  db->drawHisto("ptLept1_presel", "Lead Lepton p_{T}", "GeV", "Events", log);
  db->drawHisto("ptLept1", "Lead Lepton p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(150.);
  db->drawHisto("ptLept2_presel", "Sublead Lepton p_{T}", "GeV", "Events", log);
  db->drawHisto("ptLept2", "Sublead Lepton p_{T}", "GeV", "Events", log);

  db->set_xAxisMax(250.);
  db->drawHisto("ptJet1", "Lead Jet p_{T}", "GeV", "Events", log);
  db->drawHisto("ptJet1_prekin", "Lead Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(150.);
  db->drawHisto("ptJet2", "Sublead Jet p_{T}", "GeV", "Events", log);
  db->drawHisto("ptJet2_prekin", "Sublead Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax();
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("etaJet1", "Lead Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("etaJet2", "Sublead Jet Pseudorapidity", "", "Events", log);
  db->drawHisto("tcheJet1", "Lead Jet TCHE", "", "Events", log);
  db->drawHisto("tcheJet2", "Sublead Jet TCHE", "", "Events", log);

  db->set_rebin(10);
  db->drawHisto("ptZll_presel", "Dilepton Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("ptZjj_all_presel", "Dijet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("ptZll", "Dilepton Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("ptZjj", "Dijet Transverse Momentum", "GeV", "Events", log);

  db->set_rebin(5);
  db->drawHisto("mZll", "m_{ll}", "GeV", "Events", log);
  db->set_yAxisMaxScale( 1.3 );
  db->set_rebin(10);
  db->drawHisto("mZjj", "m_{jj}", "GeV", "Events", log);
  //db->drawHisto("mZjj_nogluetag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZjj_0btag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZjj_1btag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZjj_2btag", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("150 < m_{lljj} < 250 GeV");
  db->drawHisto("mZjj_loMass", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("250 < m_{lljj} < 400 GeV");
  db->drawHisto("mZjj_medMass", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("m_{lljj} > 400 GeV");
  db->drawHisto("mZjj_hiMass", "m_{jj}", "GeV", "Events", log);
  db->set_legendTitle("");

  db->set_rebin(1);
  db->drawHisto("mZll_presel", "m_{ll}", "GeV", "Events", log);

  db->set_rebin(20);
  db->drawHisto("mZZ_kinfit_hiMass_all", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_gluetag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag", "m_{lljj}", "GeV", "Events", log);

  db->set_legendTitle("0 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag", "m_{lljj}", "GeV", "Events", log);
  db->set_legendTitle("");

  db->set_rebin(1);
  db->drawHisto("pfMet", "Particle Flow Missing E_{T}", "GeV", "Events", log);
  db->drawHisto("mEtSig", "ME_{T} / Sum E_{T}", "", "Events", log);
  db->drawHisto("metSignificance", "ME_{T} Significance", "", "Events", log);
  db->drawHisto("metSignificance_2btag", "ME_{T} Significance", "", "Events", log);

  db->set_rebin(1);
  db->drawHisto("nChargedJet1", "Leading Jet Charged Multiplicity", "", "Events");
  db->drawHisto("nNeutralJet1", "Leading Jet Neutral Multiplicity", "", "Events");
  db->drawHisto("ptDJet1", "Leading Jet p_{T}D", "", "Events");
  db->drawHisto("nChargedJet2", "Subleading Jet Charged Multiplicity", "", "Events");
  db->drawHisto("nNeutralJet2", "Subleading Jet Neutral Multiplicity", "", "Events");
  db->drawHisto("ptDJet2", "Subleading Jet p_{T}D", "", "Events");

  db->set_rebin(4);
  db->set_yAxisMaxScale(1.6);
  db->drawHisto("QGLikelihoodJet1", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodJet2", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodProd", "Q-G Likelihood Product", "", "Events");
  db->set_legendTitle("282 < m_{lljj} < 330 GeV");
  db->drawHisto("QGLikelihoodJet1_MW300", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodJet2_MW300", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodProd_MW300", "Q-G Likelihood Product", "", "Events");
  db->set_legendTitle("376 < m_{lljj} < 440 GeV");
  db->drawHisto("QGLikelihoodJet1_MW400", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodJet2_MW400", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodProd_MW400", "Q-G Likelihood Product", "", "Events");
  db->set_legendTitle("");
  db->drawHisto("QGLikelihoodNoPUJet1", "Leading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodNoPUJet2", "Subleading Jet Q-G Likelihood", "", "Events", false, 2);
  db->drawHisto("QGLikelihoodNoPUProd", "Q-G Likelihood Product", "", "Events");
  db->set_yAxisMaxScale();

  db->set_yAxisMaxScale( 1.6 );
  db->set_rebin(3);
  db->drawHisto("cosThetaStar", "cos(#theta^{*})", "", "Events");
  db->drawHisto("cosTheta2", "cos(#theta_{2})", "", "Events");
  db->set_yAxisMaxScale( 1.8 );
  db->drawHisto("cosTheta1", "cos(#theta_{1})", "", "Events");
  db->drawHisto("phi", "#phi", "rad", "Events");
  db->drawHisto("phi1", "#phi_{1}", "rad", "Events");
  db->set_yAxisMaxScale( 1.6 );
  db->drawHisto("helicityLD", "Angular Likelihood Discriminant", "", "Events");
  //db->drawHisto("helicityLD_nogluetag", "Angular Likelihood Discriminant", "", "Events");
  db->set_yAxisMaxScale();

  db->set_rebin(1);
  db->set_legendTitle("");


  //------------
  // MUON PLOTS:
  //------------


  db->set_yAxisMaxScale( 1.6 );
  if( leptType=="ALL" || leptType=="MU" )
    db->drawHisto("mZmumu_presel", "m_{#mu#mu}", "GeV", "Events", log);
  db->set_yAxisMaxScale();

  db->set_rebin(10);
  db->set_legendTitle("Dimuon channel");
  db->drawHisto("mZjj_MU", "m_{jj}", "GeV", "Events", log);
  //db->drawHisto("mZjj_nogluetag_MU", "m_{jj}", "GeV", "Events", log);

  db->set_rebin(20);
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_gluetag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);

  db->set_legendTitle("0 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag_MU", "m_{#mu#mujj}", "GeV", "Events", log);

  db->set_rebin(1);
  db->set_legendTitle("");



  //----------------
  // ELECTRON PLOTS:
  //----------------

  
  db->set_yAxisMaxScale( 1.6 );
  if( leptType=="ALL" || leptType=="ELE" )
    db->drawHisto("mZee_presel", "m_{ee}", "GeV", "Events", log);
  db->set_yAxisMaxScale( );

  db->set_rebin(10);
  db->set_legendTitle("Dielectron channel");
  db->drawHisto("mZjj_ELE", "m_{jj}", "GeV", "Events", log);
  //db->drawHisto("mZjj_nogluetag_ELE", "m_{jj}", "GeV", "Events", log);

  db->set_rebin(20);
  db->set_legendTitle("Gluon-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_gluetag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("0 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_0btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_1btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Category");
  db->drawHisto("mZZ_kinfit_hiMass_2btag_ELE", "m_{eejj}", "GeV", "Events", log);

  db->set_legendTitle("0 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_0btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("1 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_1btag_ELE", "m_{eejj}", "GeV", "Events", log);
  db->set_legendTitle("2 b-tag Sidebands");
  db->drawHisto("mZZ_kinfit_hiMass_sidebands_2btag_ELE", "m_{eejj}", "GeV", "Events", log);




  delete db;
  db = 0;

  return 0;

}  


