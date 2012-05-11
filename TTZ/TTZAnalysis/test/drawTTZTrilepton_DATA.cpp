#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



struct ValueAndError {

  float val;
  float err;

};


ValueAndError get_ttbarSF( const DrawBase& db );
ValueAndError get_DYWZSF( const DrawBase& db );
float computeChiSquare( TH1D* h1_DATA, TH1D* h1_MC );

void drawChannelYieldPlot( DrawBase* db, const std::string& selName, char selection[], float lumi_fb, ValueAndError ttbarSF, ValueAndError DYWZSF );


int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawTTZTrilepton [(string)selType] [bTaggerType=\"TCHE\"]" << std::endl;
    exit(23);
  }

  std::string leptType = "ALL";

  std::string selType(argv[1]);

  std::string PUType = "PUHR11_73pb";

  std::string bTaggerType = "TCHE";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }

  std::string data_dataset = "DATA_Run2011_FULL";
  if( argc>=4 ) {
    std::string data_dataset_str(argv[3]);
    data_dataset = data_dataset_str;
  }


  float btag_thresh;
  if( bTaggerType == "TCHE" ) 
    btag_thresh = 3.3;
  else if( bTaggerType == "SSVHE" )
    btag_thresh = 1.74;
  else {
    std::cout << "Btagger: " << bTaggerType << " not supported. Exiting." << std::endl;
    exit(1111);
  }




  DrawBase* db = new DrawBase("TTZTrilepton");

  db->set_lumiOnRightSide(true);

  std::string outputdir_str = "TTZTrileptonPlots_" + data_dataset + "_" + selType + "_" + bTaggerType + "_" + leptType;
  db->set_outputdir(outputdir_str);

  std::string dataFileName = "TTZTrilepton_" + data_dataset;
  dataFileName += "_" + selType;
  dataFileName += "_" + bTaggerType;
  //dataFileName += "_" + PUType;
  dataFileName += "_" + leptType;
  dataFileName += ".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db->add_dataFile( dataFile, data_dataset );


  int signalFillColor = 42;

  std::string TTZFileName = "TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  TTZFileName += "_" + selType;
  TTZFileName += "_" + bTaggerType;
  //TTZFileName += "_" + PUType;
  TTZFileName += "_" + leptType;
  TTZFileName += ".root";
  TFile* TTZFile = TFile::Open(TTZFileName.c_str());
  db->add_mcFile( TTZFile, "ttZ", "t#bar{t} + Z", signalFillColor, 0);


  //std::string mcZJetsFileName = "TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  std::string mcZJetsFileName = "TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + bTaggerType;
  //mcZJetsFileName += "_" + PUType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 46, 0);

  //std::string mcTTbarFileName = "TTZTrilepton_TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  std::string mcTTbarFileName = "TTZTrilepton_TTJ_Fall11_highstat";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + bTaggerType;
  //mcTTWFileName += "_" + PUType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TTtW", "t#bar{t}", 39, 0);

  //std::string mcWZFileName = "TTZTrilepton_WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  std::string mcVVFileName = "TTZTrilepton_VV_Summer11";
  mcVVFileName += "_" + selType;
  mcVVFileName += "_" + bTaggerType;
  //mcVVFileName += "_" + PUType;
  mcVVFileName += "_" + leptType;
  mcVVFileName += ".root";
  TFile* mcVVFile = TFile::Open(mcVVFileName.c_str());
  db->add_mcFile( mcVVFile, "VV_Summer11", "Diboson", 38, 0);

//std::string mcZZFileName = "TTZTrilepton_ZZ_TuneZ2_7TeV_pythia6_tauola_Summer11-PU_S4_START42_V11-v1";
//mcZZFileName += "_" + selType;
//mcZZFileName += "_" + bTaggerType;
//mcZZFileName += "_" + leptType;
//mcZZFileName += ".root";
//TFile* mcZZFile = TFile::Open(mcZZFileName.c_str());
//db->add_mcFile( mcZZFile, "ZZtoAnything_TuneZ2", "ZZ + jets", kRed, 3004);


  std::string mcTTWFileName = "TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + bTaggerType;
  //mcTTWFileName += "_" + PUType;
  mcTTWFileName += "_" + leptType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", 33, 0);


  //db->set_shapeNormalization();
  //db->set_noStack(true);
  //db->drawHisto("nvertex_PUW", "Number of Reconstructed Vertexes", "", "Events", true);
  //db->set_noStack(false);
  //exit(1);


  float lumi_fb = 4.98;

  db->set_lumiNormalization(lumi_fb*1000.);
  db->set_noStack(false);

  db->set_drawZeros(false);


  //db->set_legendTitle("Trilepton channel");
  db->drawHisto("nvertex", "Number of Reconstructed Vertexes", "", "Events", true);
  db->drawHisto("nvertex_PUW", "Number of Reconstructed Vertexes", "", "Events", true);

  // prepresel (basically only dilepton + 3 jets):
  db->set_xAxisMin(3);
  db->drawHisto("nJets_prepresel", "Jet Multiplicity", "", "Events", true);
  db->set_xAxisMin(0);
  db->drawHisto("nBJets_loose_prepresel", "b-Jet Multiplicity (loose)", "", "Events", true);
  db->drawHisto("nBJets_medium_prepresel", "b-Jet Multiplicity (medium)", "", "Events", true);
  db->drawHisto("rhoPF_noPUW", "Particle Flow Energy Density", "GeV", "Events", true);
  db->drawHisto("rhoPF_prepresel", "Particle Flow Energy Density", "GeV", "Events", true);
  db->set_rebin(2);
  db->set_xAxisMax(130.);
  db->drawHisto("mZll_prepresel", "Dilepton Invariant Mass", "GeV", "Events", true, 2);
  db->drawHisto("mZll_prepresel_ELE", "Dielectron Invariant Mass", "GeV", "Events", true, 2);
  db->drawHisto("mZll_prepresel_MU", "Dimuon Invariant Mass", "GeV", "Events", true, 2);


  // opposite flavor leptons: control region for ttbar:
  db->set_rebin(5);
  db->drawHisto("mZll_OF_prepresel", "Opposite Flavor Dilepton Mass", "GeV", "Events");
  //db->drawHisto("mZll_OF2_prepresel", "Opposite Flavor Dilepton Mass", "GeV", "Events");
  //db->drawHisto("mZll_OF3_prepresel", "Opposite Flavor Dilepton Mass", "GeV", "Events");
  // scale ttbar MC to match data:
  ValueAndError ttbarSF = get_ttbarSF( *db );
  db->set_mcWeight( "TTtW", ttbarSF.val );
  db->drawHisto("mZll_OF_prepresel", "Opposite Flavor Dilepton Mass", "GeV", "Events", false, 1, "scaled");


  db->set_rebin();
  db->set_xAxisMax();
  db->drawHisto("nJets_presel", "Jet Multiplicity", "", "Events", true);

  db->set_rebin(4);
  db->set_xAxisMin(50.);
  db->set_xAxisMax(130.);
  db->drawHisto("mZll_presel", "Dilepton Invariant Mass", "GeV", "Events", true, 2, "noscaling" );

  // now add one lepton (prepresel -> presel)
  // anti-btag: control region for Z+Jets and WZ:
  db->drawHisto("mZll_presel_antibtag", "Dilepton Invariant Mass", "GeV", "Events", true, 2);
  ValueAndError DYWZSF = get_DYWZSF( *db );
  db->set_mcWeight( "ZJets", DYWZSF.val );
  db->set_mcWeight( "VV_Summer11", DYWZSF.val );

  db->drawHisto("mZll_presel_antibtag", "Dilepton Invariant Mass", "GeV", "Events", true, 2, "scaled");
  db->set_rebin();
  db->set_xAxisMax();
  db->set_xAxisMin();
  db->drawHisto("nJets_presel", "Jet Multiplicity", "", "Events", true, 1, "scaled");


  db->set_rebin(2);
  db->drawHisto("rhoPF_presel", "Particle Flow Energy Density", "GeV", "Events", true);


  db->set_rebin();
  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 20 GeV)", "", "Events");

  db->set_rebin(5);
  db->drawHisto("ptLept3_presel", "Third Lepton p_{T}", "GeV", "Events");
  db->drawHisto("etaLept3_presel", "Third Lepton #eta", "", "Events");
  db->drawHisto("leptTypeLept3_presel", "Third Lepton Flavor", "", "Events");
  db->drawHisto("combinedIsoRelLept3_presel", "Third Lepton Isolation", "", "Events");

  db->set_rebin(4);
  db->set_xAxisMax(130);
  db->drawHisto("mZll_presel", "Dilepton Invariant Mass", "GeV", "Events", true, 2, "scaled");
  db->drawHisto("mZll", "Dilepton Invariant Mass", "GeV", "Events", true, 2);
  //db->drawHisto_fromTree("tree_passedEvents", "mZll", "eventWeight*(nBjets_medium==0 && isMZllSignalRegion)", 100, 50., 130., "mZll_antibtag", "Dilepton Invariant Mass", "GeV");
  //float DYWZSF2 = get_DYWZSF( *db );
  db->set_xAxisMax();


  ValueAndError noCorr;
  noCorr.val = 1.;
  noCorr.err = 0.;

  drawChannelYieldPlot( db, "noScaling", "", lumi_fb, noCorr, noCorr );
  drawChannelYieldPlot( db, "", "", lumi_fb, ttbarSF, DYWZSF );
  drawChannelYieldPlot( db, "ptZll80", "eventWeight*(ptZll>80.)", lumi_fb, ttbarSF, DYWZSF );


/*
  std::vector<TH1D*> lastHistosMC;
  float signalYield;
  float channelYieldGen = lumi_fb * 139 * 0.06 * 0.2 * 0.67;

  db->set_rebin(2);
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>81. && mZll<101.)", 50, 0., 500., "ptZll", "p_{T} (Z)", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>81. && mZll<101. && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20. )", 50, 0., 500., "ptZll_jetpt20", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << 100.*signalYield/channelYieldGen << "%)" << std::endl;
  char selection[1000];
  sprintf( selection, "eventWeight*(mZll>81. && mZll<101. && bTagJetB1>%f && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", btag_thresh);
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", selection, 30, 0., 300., "ptZll_jetpt20_btag", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << 100.*signalYield/channelYieldGen << "%)" << std::endl;
  sprintf(selection, "eventWeight*(mZll>81. && mZll<101. && pfMet>30. && bTagJetB1>%f && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", btag_thresh);
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", selection, 30, 0., 300., "ptZll_jetpt20_btag_met", "p_{T} (Z)", "GeV", "Events");
  drawChannelYieldPlot( db, "jetpt20_btag_met", selection, lumi_fb, ttbarSF, DYWZSF );
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << 100.*signalYield/channelYieldGen << "%)" << std::endl;


  std::string yieldsFileName = "yieldsData_"+selType+"_"+bTaggerType+".txt";
  ofstream yieldsFile(yieldsFileName.c_str());

  yieldsFile << "------------------------" << std::endl;
  yieldsFile << "Yields @ " << lumi_fb << " fb-1" << std::endl;
  yieldsFile << "------------------------" << std::endl;

  float s = 0.;
  float b  = 0.;

  for( unsigned i=0; i<db->get_lastHistos_mc().size(); ++i )  {

    float yield = db->get_lastHistos_mc()[i]->Integral();

    if(  db->get_mcFiles()[i].legendName=="t#bar{t} + Z" ) {

      s += yield;

    } else {

      b += yield;

    }

    yieldsFile << db->get_mcFiles()[i].legendName.c_str();
    if( db->get_mcFiles()[i].legendName.size()<12 ) yieldsFile << "\t";
    yieldsFile << "\t" << yield << std::endl;

  }
    
  yieldsFile << "Total Background\t& " << b << std::endl;
  yieldsFile << "Total (S+B)\t\t& " << b+s << std::endl;

  //yieldsFile << "s/sqrt(b)    \t& " << s_mumu/sqrt(b_mumu) << "\t& " << s_ee/sqrt(b_ee) << "\t& " << s_emu/sqrt(b_emu) << "\\\\" << std::endl;
  yieldsFile << "Observed:    \t& " << db->get_lastHistos_data()[0]->Integral() << std::endl;



  sprintf( selection, "eventWeight*(ptZll>100. && mZll>81. && mZll<101. && pfMet>30. && bTagJetB1>%f && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", btag_thresh);
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", selection, 30, 0., 300., "ptZll_jetpt20_btag_met_ptZll", "p_{T} (Z)", "GeV", "Events");
  drawChannelYieldPlot( db, "jetpt20_btag_met_ptZll", selection, lumi_fb, ttbarSF, DYWZSF );
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << 100.*signalYield/channelYieldGen << "%)" << std::endl;



  s = 0.;
  b  = 0.;


  yieldsFile << std::endl << std::endl << "requiring pt(Zll)>100 GeV: " << std::endl;
  

  for( unsigned i=0; i<db->get_lastHistos_mc().size(); ++i )  {

    float yield = db->get_lastHistos_mc()[i]->Integral();

    if(  db->get_mcFiles()[i].legendName=="t#bar{t} + Z" ) {

      s += yield;

    } else {

      b += yield;

    }

    yieldsFile << db->get_mcFiles()[i].legendName.c_str();
    if( db->get_mcFiles()[i].legendName.size()<12 ) yieldsFile << "\t";
    yieldsFile << "\t" << yield << std::endl;

  }
    
  yieldsFile << "Total Background\t& " << b << std::endl;
  yieldsFile << "Total (S+B)\t\t& " << b+s << std::endl;

  //yieldsFile << "s/sqrt(b)    \t& " << s_mumu/sqrt(b_mumu) << "\t& " << s_ee/sqrt(b_ee) << "\t& " << s_emu/sqrt(b_emu) << "\\\\" << std::endl;
  yieldsFile << "Observed:    \t& " << db->get_lastHistos_data()[0]->Integral() << std::endl;
  yieldsFile.close();

*/





  delete db;
  db = 0;

  return 0;

}  



void drawChannelYieldPlot( DrawBase* db, const std::string& selName, char selection[], float lumi_fb, ValueAndError ttbarSF, ValueAndError DYWZSF ) {

  TH1F::AddDirectory(kTRUE);

  std::string selection_str(selection);

  if( selection_str=="" ) {

    selection_str = "eventWeight*(";

  } else {

    //get rid of last parenthesis:
    selection_str.erase(selection_str.end()-1);
    selection_str += " && ";

  }

  // and now lepton channel definitions:
  std::string sel_mmm = selection_str + " isMZllSignalRegion && passed_btag && leptType==0 && leptType3==0)";
  std::string sel_mme = selection_str + " isMZllSignalRegion && passed_btag && leptType==0 && leptType3==1)";
  std::string sel_eem = selection_str + " isMZllSignalRegion && passed_btag && leptType==1 && leptType3==0)";
  std::string sel_eee = selection_str + " isMZllSignalRegion && passed_btag && leptType==1 && leptType3==1)";


  TTree* tree_data = (TTree*)(db->get_dataFile(0).file->Get("tree_passedEvents"));

  TH1D* h1_data_mmm = new TH1D("data_mmm", "", 100, 0., 10000.);
  TH1D* h1_data_mme = new TH1D("data_mme", "", 100, 0., 10000.);
  TH1D* h1_data_eem = new TH1D("data_eem", "", 100, 0., 10000.);
  TH1D* h1_data_eee = new TH1D("data_eee", "", 100, 0., 10000.);

  tree_data->Project("data_mmm", "ptZll", sel_mmm.c_str());
  tree_data->Project("data_mme", "ptZll", sel_mme.c_str());
  tree_data->Project("data_eem", "ptZll", sel_eem.c_str());
  tree_data->Project("data_eee", "ptZll", sel_eee.c_str());

  float mmm_data = h1_data_mmm->Integral();
  float mme_data = h1_data_mme->Integral();
  float eem_data = h1_data_eem->Integral();
  float eee_data = h1_data_eee->Integral();
  float tot_data = mmm_data + mme_data + eem_data + eee_data;


  TH1D* h1_yields_data = new TH1D("yields_data", "", 5, 0., 5.);
  h1_yields_data->SetBinContent( 1, eee_data );
  h1_yields_data->SetBinContent( 2, eem_data );
  h1_yields_data->SetBinContent( 3, mme_data );
  h1_yields_data->SetBinContent( 4, mmm_data );
  h1_yields_data->SetBinContent( 5, tot_data );
  

  TGraphAsymmErrors* gr_data = fitTools::getGraphPoissonErrors(h1_yields_data);
  gr_data->SetMarkerStyle(20);
  gr_data->SetMarkerSize(1.4);


  //TLegend* legend = new TLegend( 0.6, 0.57, 0.92, 0.9 );
  TLegend* legend = new TLegend( 0.2, 0.57, 0.52, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.042);
  legend->AddEntry( gr_data, "Data", "P" );


  THStack* stackMC = new THStack();
  std::vector<TH1D*> vh1_yields_mc;
  TH1D* h1_yields_mc_totBG = new TH1D("yields_mc_totBG", "", 5, 0., 5.);
  h1_yields_mc_totBG->Sumw2();
  TH1D* h1_yields_mc_signal = new TH1D("yields_mc_signal", "", 5, 0., 5.);
  h1_yields_mc_signal->Sumw2();


  for( unsigned i=0; i<db->get_mcFiles().size(); ++i) {

    // reverse order:
    int iMC = db->get_mcFiles().size()-i-1;

    TTree* tree_mc = (TTree*)(db->get_mcFile(iMC).file->Get("tree_passedEvents"));

    TH1D* h1_mc_mmm = new TH1D("mc_mmm", "", 100, 0., 10000.);
    TH1D* h1_mc_mme = new TH1D("mc_mme", "", 100, 0., 10000.);
    TH1D* h1_mc_eem = new TH1D("mc_eem", "", 100, 0., 10000.);
    TH1D* h1_mc_eee = new TH1D("mc_eee", "", 100, 0., 10000.);

    h1_mc_mmm->Sumw2();
    h1_mc_mme->Sumw2();
    h1_mc_eem->Sumw2();
    h1_mc_eee->Sumw2();

    tree_mc->Project("mc_mmm", "ptZll", sel_mmm.c_str());
    tree_mc->Project("mc_mme", "ptZll", sel_mme.c_str());
    tree_mc->Project("mc_eem", "ptZll", sel_eem.c_str());
    tree_mc->Project("mc_eee", "ptZll", sel_eee.c_str());

    float mmm_mc = h1_mc_mmm->Integral();
    float mme_mc = h1_mc_mme->Integral();
    float eem_mc = h1_mc_eem->Integral();
    float eee_mc = h1_mc_eee->Integral();
    float tot_mc = mmm_mc + mme_mc + eem_mc + eee_mc;

    float scaling = lumi_fb*1000.;
    float scaling_err = 0.;
    if( db->get_mcFile(iMC).datasetName=="TTtW" ) {
      scaling *= ttbarSF.val;
      scaling_err = ttbarSF.err;
    }
    if( db->get_mcFile(iMC).datasetName=="ZJets" || db->get_mcFile(iMC).datasetName=="VV_Summer11" ) {
      scaling *= DYWZSF.val;
      scaling_err = DYWZSF.err;
    }

    eee_mc *= scaling;
    eem_mc *= scaling;
    mme_mc *= scaling;
    mmm_mc *= scaling;
    tot_mc *= scaling;

    char hname[100];
    sprintf( hname, "yields_mc_%d", iMC);
    TH1D* h1_yields_mc = new TH1D(hname, "", 5, 0., 5.);
    h1_yields_mc->SetBinContent( 1, eee_mc );
    h1_yields_mc->SetBinContent( 2, eem_mc );
    h1_yields_mc->SetBinContent( 3, mme_mc );
    h1_yields_mc->SetBinContent( 4, mmm_mc );
    h1_yields_mc->SetBinContent( 5, tot_mc );
  

    // add in quadrature scaling error:
    float oldErr_eee = h1_yields_mc->GetBinError(1);
    float oldErr_eem = h1_yields_mc->GetBinError(2);
    float oldErr_mme = h1_yields_mc->GetBinError(3);
    float oldErr_mmm = h1_yields_mc->GetBinError(4);
    float oldErr_tot = h1_yields_mc->GetBinError(5);
 
    float newErr_eee = sqrt( oldErr_eee*oldErr_eee + eee_mc*eee_mc*scaling_err*scaling_err );
    float newErr_eem = sqrt( oldErr_eem*oldErr_eem + eem_mc*eem_mc*scaling_err*scaling_err );
    float newErr_mme = sqrt( oldErr_mme*oldErr_mme + mme_mc*mme_mc*scaling_err*scaling_err );
    float newErr_mmm = sqrt( oldErr_mmm*oldErr_mmm + mmm_mc*mmm_mc*scaling_err*scaling_err );
    float newErr_tot = sqrt( oldErr_tot*oldErr_tot + tot_mc*tot_mc*scaling_err*scaling_err );

    h1_yields_mc->SetBinError( 1, newErr_eee );
    h1_yields_mc->SetBinError( 2, newErr_eem );
    h1_yields_mc->SetBinError( 3, newErr_mme );
    h1_yields_mc->SetBinError( 4, newErr_mmm );
    h1_yields_mc->SetBinError( 5, newErr_tot );
  

    h1_yields_mc->SetFillColor( db->get_mcFile(iMC).fillColor );

    vh1_yields_mc.push_back( (TH1D*)h1_yields_mc );

    stackMC->Add(h1_yields_mc, "HISTO");

    if( db->get_mcFile(iMC).datasetName!="ttZ" && db->get_mcFile(iMC).datasetName!="ttW" )
      h1_yields_mc_totBG->Add(h1_yields_mc);
    else 
      h1_yields_mc_signal->Add(h1_yields_mc);

    delete h1_mc_mmm;
    delete h1_mc_mme;
    delete h1_mc_eem;
    delete h1_mc_eee;
    
  }

 
  for( unsigned i=0; i<db->get_mcFiles().size(); ++i) {

    int inverseIndex = db->get_mcFiles().size()-i-1;
    if( vh1_yields_mc[inverseIndex]->Integral()>0 )
      legend->AddEntry( vh1_yields_mc[inverseIndex], db->get_mcFile(i).legendName.c_str(), "F" );

  } 



  float yMaxData = h1_yields_data->GetMaximum();
  float yMaxMC = stackMC->GetMaximum();
  float yMax = TMath::Max(yMaxData, yMaxMC);
  yMax *= 1.4;

  TH2D* h2_axes = new TH2D("axes", "", 5, 0., 5., 10, 0., yMax);
  h2_axes->GetXaxis()->SetLabelSize(0.085);
  h2_axes->GetXaxis()->SetBinLabel(1, "(ee)e");
  h2_axes->GetXaxis()->SetBinLabel(2, "(ee)#mu");
  h2_axes->GetXaxis()->SetBinLabel(3, "(#mu#mu)e");
  h2_axes->GetXaxis()->SetBinLabel(4, "(#mu#mu)#mu");
  h2_axes->GetXaxis()->SetBinLabel(5, "Total" );
  h2_axes->SetYTitle("Events");


  TPaveText* label_sqrt = db->get_labelSqrt();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  h2_axes->Draw();
  stackMC->Draw("histo same");
  legend->Draw("same");
  gr_data->Draw("P same");
  label_sqrt->Draw("same");

  gPad->RedrawAxis();
  
  char canvasName[500];
  if( selName!="" )
    sprintf( canvasName, "%s/channelYields_%s.eps", db->get_outputdir().c_str(), selName.c_str() );
  else
    sprintf( canvasName, "%s/channelYields.eps", db->get_outputdir().c_str() );
  c1->SaveAs(canvasName);

  
  std::string yieldfilename = db->get_outputdir() + "/yields";
  if( selName!="" ) yieldfilename = yieldfilename + "_" + selName;
  yieldfilename = yieldfilename + ".txt";

  ofstream ofs(yieldfilename.c_str());

  ofs << "channel\tobserved\tsignal\tb_pred\tb_pred_error" << std::endl;
  ofs << "(ee)e  \t" << h1_yields_data->GetBinContent(1) << "\t\t" << h1_yields_mc_signal->GetBinContent(1) << "\t" << h1_yields_mc_totBG->GetBinContent(1) << "\t" << h1_yields_mc_totBG->GetBinError(1) << std::endl;
  ofs << "(ee)m  \t" << h1_yields_data->GetBinContent(2) << "\t\t" << h1_yields_mc_signal->GetBinContent(2) << "\t" << h1_yields_mc_totBG->GetBinContent(2) << "\t" << h1_yields_mc_totBG->GetBinError(2) << std::endl;
  ofs << "(mm)e  \t" << h1_yields_data->GetBinContent(3) << "\t\t" << h1_yields_mc_signal->GetBinContent(3) << "\t" << h1_yields_mc_totBG->GetBinContent(3) << "\t" << h1_yields_mc_totBG->GetBinError(3) << std::endl;
  ofs << "(mm)m  \t" << h1_yields_data->GetBinContent(4) << "\t\t" << h1_yields_mc_signal->GetBinContent(4) << "\t" << h1_yields_mc_totBG->GetBinContent(4) << "\t" << h1_yields_mc_totBG->GetBinError(4) << std::endl;
  ofs << "Total  \t" << h1_yields_data->GetBinContent(5) << "\t\t" << h1_yields_mc_signal->GetBinContent(5) << "\t" << h1_yields_mc_totBG->GetBinContent(5) << "\t" << h1_yields_mc_totBG->GetBinError(5) << std::endl;


  ofs.close();


  delete c1;
  delete h2_axes;
  delete gr_data;
  delete h1_yields_data;
  delete h1_yields_mc_totBG;
  delete h1_yields_mc_signal;

  delete h1_data_mmm;
  delete h1_data_mme;
  delete h1_data_eem;
  delete h1_data_eee;

  for( unsigned i=0; i<vh1_yields_mc.size(); ++i ) 
    delete vh1_yields_mc[vh1_yields_mc.size()-i-1];

}




ValueAndError get_ttbarSF( const DrawBase& db ) {


  //TH1D* h1_DATA = new TH1D(*(db->get_lastHistos_data()[0]));
  TH1D* h1_DATA = new TH1D(*(db.get_lastHistos_data()[0]));

  std::vector< TH1D* > lastHistosMC = db.get_lastHistos_mc();
  TH1D* h1_TTJets;

  std::vector< TH1D* > otherHistosMC;
  // create short list of histos withous ttbar:
  for( unsigned iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {

    if( db.get_mcFiles()[iHisto].datasetName!="TTtW" ) {

      TH1D* h1_newHisto = new TH1D(*(lastHistosMC[iHisto]));
      otherHistosMC.push_back( h1_newHisto );

    } else {

      h1_TTJets = new TH1D(*(lastHistosMC[iHisto]));

    }

  }


  float minSF = 0.9;
  float maxSF = 1.5;

  int nSteps = 1000;
  float step = (maxSF-minSF)/(float)nSteps;

  TH1D* h1_chiSquare = new TH1D("chiSquare", "", nSteps, minSF, maxSF );


  float minChiSquare = 9999.;
  float minChiSquare_bin = 0;
  float foundSF = -1.;

  for( unsigned i=0; i<h1_chiSquare->GetNbinsX(); ++i ) {

    float thisSF = h1_chiSquare->GetBinCenter(i+1);

    TH1D* h1_TTJets_sf = new TH1D(*h1_TTJets);
    h1_TTJets_sf->Scale( thisSF );

    for( unsigned ihisto=0; ihisto<otherHistosMC.size(); ++ihisto )
      h1_TTJets_sf->Add( otherHistosMC[ihisto] );

    float thisChiSquare = computeChiSquare( h1_DATA, h1_TTJets_sf );

    h1_chiSquare->SetBinContent( i+1, thisChiSquare ); 

    if( thisChiSquare < minChiSquare ) {
      minChiSquare = thisChiSquare;
      minChiSquare_bin = i+1;
      foundSF = thisSF;
    }

    delete h1_TTJets_sf;

  }

  // now search for error on SF (minChiSquare+1)
  float minChiSquare_oneSigma = minChiSquare+1.;
  float bestDiff_right = 999.;
  float bestDiff_left = 999.;
  float foundSF_minus = -1.;
  float foundSF_plus = -1.;

  for( unsigned iBin=1; iBin<h1_chiSquare->GetNbinsX()+1; ++iBin ) {

    float thisSF = h1_chiSquare->GetBinCenter(iBin);
    float thisChiSquare = h1_chiSquare->GetBinContent(iBin);
    
    if( iBin<minChiSquare_bin ) { //left side of parabola

      if( fabs(thisChiSquare-minChiSquare_oneSigma) < bestDiff_left ) {
        bestDiff_left = fabs(thisChiSquare-minChiSquare_oneSigma);
        foundSF_minus = thisSF;
      }

    } else { //right side
  
      if( fabs(thisChiSquare-minChiSquare_oneSigma) < bestDiff_right ) {
        bestDiff_right = fabs(thisChiSquare-minChiSquare_oneSigma);
        foundSF_plus = thisSF;
      }

    } // left-right

  } // for bins


  float err_plus = foundSF_plus-foundSF;
  float err_minus = foundSF-foundSF_minus;

  float err = 0.5*(err_plus+err_minus);

  float yMax = 1.1*h1_chiSquare->GetMaximum();

  h1_chiSquare->SetXTitle( "t#bar{t} Scale Factor" );
  h1_chiSquare->SetYTitle( "Normalized #chi^{2} (Data-MC)" );
  h1_chiSquare->SetMarkerStyle( 20 );
  h1_chiSquare->SetMarkerSize( 1.6 );
  h1_chiSquare->SetMarkerColor( 46 );
  h1_chiSquare->GetYaxis()->SetRangeUser(0., yMax);

  std::cout << std::endl << "-> TTbar Scaling" << std::endl;
  std::cout << "Min Chi Square: " << minChiSquare << std::endl;
  std::cout << "Best SF: " << foundSF << " +/- " << err << std::endl;


  TLine* line_plus = new TLine( foundSF_plus, 0., foundSF_plus, yMax );
  TLine* line_minus = new TLine( foundSF_minus, 0., foundSF_minus, yMax );
  TLine* line_mean = new TLine( foundSF, 0., foundSF, yMax );

  line_plus->SetLineWidth( 2 );
  line_mean->SetLineWidth( 2 );
  line_minus->SetLineWidth( 2 );
  
  line_plus->SetLineStyle( 2 );
  line_minus->SetLineStyle( 2 );
  
  line_plus->SetLineColor( 38 );
  line_mean->SetLineColor( kBlack );
  line_minus->SetLineColor( 38 );

  TPaveText* label_sqrt = db.get_labelSqrt();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  
  
  h1_chiSquare->Draw( "P" );
  line_plus->Draw("same");
  line_mean->Draw("same");
  line_minus->Draw("same");
  label_sqrt->Draw( "Psame" );
  h1_chiSquare->Draw( "Psame" );
  
  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/ttbarChiSquareScan.eps", db.get_outputdir().c_str() );
  c1->SaveAs(canvasName);

  delete c1;
  delete h1_chiSquare;
  delete h1_TTJets;
  delete h1_DATA;

  ValueAndError ve;
  ve.val = foundSF;
  ve.err = err;

  return ve;

}




float computeChiSquare( TH1D* h1_DATA, TH1D* h1_MC ) {

  float chiSquare=0.;
  int nbins = h1_DATA->GetNbinsX();

  for( unsigned ibin=1; ibin<nbins+1; ++ibin ) {

    float data = h1_DATA->GetBinContent(ibin);
    float mc = h1_MC->GetBinContent(ibin);

    float binDifference = fabs(data - mc);
    float err = sqrt(data); // fairly high stat in every bin anyways

    float addendum = binDifference / err;

    chiSquare += (addendum*addendum);
   
  } //for bins

  chiSquare /= (nbins-1.);

  return chiSquare;

}




ValueAndError get_DYWZSF( const DrawBase& db ) {


  TH1D* h1_DATA = new TH1D(*(db.get_lastHistos_data()[0]));
  float dataIntegral = h1_DATA->Integral();

  std::vector< TH1D* > lastHistosMC = db.get_lastHistos_mc();

  float mcIntegral = 0.;
  for( unsigned iHisto=0; iHisto<lastHistosMC.size(); ++iHisto )
      mcIntegral += lastHistosMC[iHisto]->Integral();

  float sf = dataIntegral/mcIntegral;
  float sf_err = sqrt(dataIntegral)/mcIntegral;

  std::cout << std::endl << "-> DY/WZ Scaling" << std::endl;
  std::cout << "SF: " << sf << " +/- " << sf_err << std::endl;

  ValueAndError ve;
  ve.val = sf;
  ve.err = sf_err;

  return ve;

}
