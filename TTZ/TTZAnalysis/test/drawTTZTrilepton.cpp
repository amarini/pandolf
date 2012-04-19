#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



void drawChannelYieldPlot( DrawBase* db, const std::string& selType, const std::string& bTaggerType, float lumi_fb, const std::string& saveName, std::string additionalCuts );


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

  db->set_lumiOnRightSide();

  std::string outputdir_str = "TTZTrileptonPlots_MConly_" + selType + "_" + bTaggerType + "_" + leptType;
  db->set_outputdir(outputdir_str);

  int signalFillColor = 42;

  std::string TTZFileName = "TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  TTZFileName += "_" + selType;
  TTZFileName += "_" + bTaggerType;
  //TTZFileName += "_" + PUType;
  TTZFileName += "_" + leptType;
  TTZFileName += ".root";
  TFile* TTZFile = TFile::Open(TTZFileName.c_str());
  db->add_mcFile( TTZFile, "ttZ", "t#bar{t} + Z", signalFillColor, 3005);


  std::string mcZJetsFileName = "TTZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcZJetsFileName += "_" + selType;
  mcZJetsFileName += "_" + bTaggerType;
  //mcZJetsFileName += "_" + PUType;
  mcZJetsFileName += "_" + leptType;
  mcZJetsFileName += ".root";
  TFile* mcZJetsFile = TFile::Open(mcZJetsFileName.c_str());
  db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 46, 3004);
  //db->add_mcFile( mcZJetsFile, "ZJets", "Z + jets", 30, 3004);

  //std::string mcTTbarFileName = "TTZTrilepton_TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  std::string mcTTbarFileName = "TTZTrilepton_TTJ_Fall11_highstat";
  mcTTbarFileName += "_" + selType;
  mcTTbarFileName += "_" + bTaggerType;
  //mcTTWFileName += "_" + PUType;
  mcTTbarFileName += "_" + leptType;
  mcTTbarFileName += ".root";
  TFile* mcTTbarFile = TFile::Open(mcTTbarFileName.c_str());
  db->add_mcFile( mcTTbarFile, "TTtW", "t#bar{t}", 39, 0);




  db->set_shapeNormalization();




  bool log = true;


  //db->set_legendTitle("Trilepton channel");

  db->drawHisto("nJets", "Jet Multiplicity (p_{T} > 10 GeV)", "", "Events", log);
  db->drawHisto("mZll_prepresel", "Dilepton mass", "GeV", "Events");
  db->drawHisto("mZll", "Dilepton mass", "GeV", "Events");


  db->set_rebin(10); 
  db->set_xAxisMax(200.); 
  db->drawHisto("pfMet", "pfME_{T}", "GeV", "Events", log);
  db->drawHisto("pfMet_presel", "pfME_{T}", "GeV", "Events", log);
  db->set_rebin(4); 
  db->set_xAxisMax(250.); 
  db->set_xAxisMax();
  db->drawHisto("ptLept3_presel", "Third Lepton p_{T}", "GeV", "Events");
  db->set_xAxisMax(0.15);
  db->set_rebin(); 
  db->drawHisto("combinedIsoRelLept3_presel", "Third Lepton Isolation", "", "Events");


  db->set_xAxisMax();
  db->set_rebin(10);
  db->drawHisto_fromTree("tree_passedEvents", "TMath::Max( TMath::Max(ptJetB1, ptJetB2), TMath::Max(ptJet3, ptJet4) )", "eventWeight", 100, 20., 420., "ptJetMax", "Leading Jet p_{T}", "GeV");
  db->set_xAxisMax(250.);
  db->drawHisto("ptLeptZ1", "Lead Z Lepton p_{T}", "GeV", "Events");
  db->drawHisto("ptJetB1", "Leading b-Tagged Jet p_{T}", "GeV");
  db->set_xAxisMax(150.);
  db->drawHisto("ptLeptZ2", "Sublead Z Lepton p_{T}", "GeV", "Events");
  
  db->set_xAxisMax(5.);
  db->set_rebin(5);
  std::string axisName = "Leading b-Tag (" + bTaggerType + ")";
  db->drawHisto("bTagJetB1", axisName );
  axisName = "Subleading b-Tag (" + bTaggerType + ")";
  db->drawHisto("bTagJetB2", "Subleading b-Tag");

  db->set_xAxisMax();
  db->set_rebin(20);
  db->drawHisto("mTW", "", "GeV", "Events");
  db->drawHisto("mT_lZmet", "", "GeV", "Events");

  db->drawHisto("deltaRbb", "", "GeV", "Events");

  db->drawHisto("mb1jj", "", "GeV", "Events");
  db->drawHisto("mb2jj", "", "GeV", "Events");
  db->drawHisto("mbjj_best", "", "GeV", "Events");

  db->drawHisto("mb1jjZ", "", "GeV", "Events");
  db->drawHisto("mb2jjZ", "", "GeV", "Events");
  db->drawHisto("mbjjZ_best", "", "GeV", "Events");

  db->drawHisto("mTb1W", "", "GeV", "Events");
  db->drawHisto("mTb2W", "", "GeV", "Events");
  db->drawHisto("mTbW_best", "", "GeV", "Events");

  db->drawHisto("mTb1WZ", "", "GeV", "Events");
  db->drawHisto("mTb2WZ", "", "GeV", "Events");
  db->drawHisto("mTbWZ_best", "", "GeV", "Events");





  std::string mcWZFileName = "TTZTrilepton_WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1";
  mcWZFileName += "_" + selType;
  mcWZFileName += "_" + bTaggerType;
  //mcWZFileName += "_" + PUType;
  mcWZFileName += "_" + leptType;
  mcWZFileName += ".root";
  TFile* mcWZFile = TFile::Open(mcWZFileName.c_str());
  db->add_mcFile( mcWZFile, "WZtoAnything_TuneZ2", "WZ + jets", 38, 3004);


  std::string mcTTWFileName = "TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + bTaggerType;
  //mcTTWFileName += "_" + PUType;
  mcTTWFileName += "_" + leptType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", 33, 3002);


  float lumifb = 4.6;

  db->set_lumiNormalization(lumifb*1000.);
  db->set_noStack(false);

  db->set_rebin(2);
  db->drawHisto( "ptZll", "Z p_{T}", "GeV" );
  drawChannelYieldPlot( db, selType, bTaggerType, lumifb, "", "" );
  drawChannelYieldPlot( db, selType, bTaggerType, lumifb, "ptZll80", "ptZll > 80." );



  std::vector<TH1D*> lastHistosMC;
  float signalYield;
  float channelYieldGen = lumifb * 139 * 0.06 * 0.2 * 0.67;


  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>81. && mZll<101.)", 30, 0., 300., "ptZll", "p_{T} (Z)", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", "eventWeight*(mZll>81. && mZll<101. && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20. )", 30, 0., 300., "ptZll_jetpt20", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;
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
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;
  sprintf(selection, "eventWeight*(mZll>81. && mZll<101. && pfMet>30. && bTagJetB1>%f && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", btag_thresh);
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", selection, 30, 0., 300., "ptZll_jetpt20_btag_met", "p_{T} (Z)", "GeV", "Events");
  lastHistosMC = db->get_lastHistos_mc();
  for( unsigned int iHisto=0; iHisto<lastHistosMC.size(); ++iHisto ) {
    if( lastHistosMC[iHisto]->GetFillColor()==signalFillColor ) {
      signalYield = lastHistosMC[iHisto]->Integral(1, lastHistosMC[iHisto]->GetNbinsX()+1);
      break;
    }
  }
  std::cout << "Signal yield: " << signalYield << " (" << signalYield/channelYieldGen << "%)" << std::endl;

  sprintf( selection, "eventWeight*(ptZll>100. && mZll>81. && mZll<101. && pfMet>30. && bTagJetB1>%f && ptJetB1>20. && ptJetB2>20. && ptJet3>20. && ptJet4>20.)", btag_thresh);
  db->drawHisto_fromTree("tree_passedEvents", "ptZll", selection, 30, 0., 300., "ptZll_jetpt20_btag_met_ptZll", "p_{T} (Z)", "GeV", "Events");
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


void drawChannelYieldPlot( DrawBase* db, const std::string& selType, const std::string& bTaggerType, float lumi_fb, const std::string& saveName, std::string additionalCuts ) {


  if( additionalCuts!="" ) additionalCuts += " && ";

  // and now lepton channel definitions:
  std::string sel_mmm = "eventWeight*( " + additionalCuts + " isMZllSignalRegion && leptType==0 && leptType3==0)";
  std::string sel_mme = "eventWeight*( " + additionalCuts + " isMZllSignalRegion && leptType==0 && leptType3==1)";
  std::string sel_eem = "eventWeight*( " + additionalCuts + " isMZllSignalRegion && leptType==1 && leptType3==0)";
  std::string sel_eee = "eventWeight*( " + additionalCuts + " isMZllSignalRegion && leptType==1 && leptType3==1)";



  TLegend* legend = new TLegend( 0.6, 0.57, 0.92, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.042);


  THStack* stackMC = new THStack();
  std::vector<TH1D*> vh1_yields_mc;

  std::string yieldsFileName = "yields_"+selType;
  if( saveName!= "" ) yieldsFileName = yieldsFileName + "_" + saveName;
  yieldsFileName = yieldsFileName +"_"+bTaggerType+".txt";

  ofstream yieldsFile(yieldsFileName.c_str());

  yieldsFile << "------------------------" << std::endl;
  yieldsFile << "Yields @ " << lumi_fb << " fb-1" << std::endl;
  yieldsFile << "------------------------" << std::endl;

  yieldsFile << std::endl;
  yieldsFile << "                    (mm)m   \t(mm)e   \t(ee)m   \t(ee)e   \tTotal" << std::endl << std::endl;


  float s_mmm = 0.;
  float s_mme = 0.;
  float s_eem = 0.;
  float s_eee = 0.;
  float s = 0.;

  float b_mmm = 0.;
  float b_mme = 0.;
  float b_eem = 0.;
  float b_eee = 0.;
  float b  = 0.;


  for( unsigned i=0; i<db->get_mcFiles().size(); ++i) {

    // reverse order:
    int iMC = db->get_mcFiles().size()-i-1;

    TTree* tree_mc = (TTree*)(db->get_mcFile(iMC).file->Get("tree_passedEvents"));

    TH1D* h1_mc_mmm = new TH1D("mc_mmm", "", 100, 0., 10000.);
    TH1D* h1_mc_mme = new TH1D("mc_mme", "", 100, 0., 10000.);
    TH1D* h1_mc_eem = new TH1D("mc_eem", "", 100, 0., 10000.);
    TH1D* h1_mc_eee = new TH1D("mc_eee", "", 100, 0., 10000.);

    tree_mc->Project("mc_mmm", "ptZll", sel_mmm.c_str());
    tree_mc->Project("mc_mme", "ptZll", sel_mme.c_str());
    tree_mc->Project("mc_eem", "ptZll", sel_eem.c_str());
    tree_mc->Project("mc_eee", "ptZll", sel_eee.c_str());

    float mmm = lumi_fb*1000.*h1_mc_mmm->Integral();
    float mme = lumi_fb*1000.*h1_mc_mme->Integral();
    float eem = lumi_fb*1000.*h1_mc_eem->Integral();
    float eee = lumi_fb*1000.*h1_mc_eee->Integral();
    float total = mmm + mme + eem + eee;

    char hname[100];
    sprintf( hname, "yields_mc_%d", iMC);
    TH1D* h1_yields_mc = new TH1D(hname, "", 4, 0., 4.);
    h1_yields_mc->SetBinContent( 1, mmm );
    h1_yields_mc->SetBinContent( 2, mme );
    h1_yields_mc->SetBinContent( 3, eem );
    h1_yields_mc->SetBinContent( 4, eee );

    h1_yields_mc->SetFillColor( db->get_mcFile(iMC).fillColor );

    vh1_yields_mc.push_back( (TH1D*)h1_yields_mc );

    stackMC->Add(h1_yields_mc, "HISTO");



    if(  db->get_mcFiles()[iMC].legendName=="t#bar{t} + Z" ) {

      s_mmm += mmm;
      s_mme += mme;
      s_eem += eem;
      s_eee += eee;
      s += total;

    } else {

      yieldsFile << db->get_mcFiles()[iMC].legendName;
      for( unsigned ichar=0; ichar<20-db->get_mcFiles()[iMC].legendName.size(); ++ichar ) yieldsFile << " ";
      yieldsFile << Form("%.4f \t %.4f \t %.4f \t %.4f \t %.4f", mmm, mme, eem, eee, total) << std::endl;

      b_mmm += mmm;
      b_mme += mme;
      b_eem += eem;
      b_eee += eee;
      b += total;

    }

    delete h1_mc_mmm;
    delete h1_mc_mme;
    delete h1_mc_eem;
    delete h1_mc_eee;
    
  }
    

  yieldsFile << "Total BG            " << Form("%.4f \t %.4f \t %.4f \t %.4f \t %.4f", b_mmm, b_mme, b_eem, b_eee, b) << std::endl;
  yieldsFile << std::endl;
  yieldsFile << "Signal              " << Form("%.4f \t %.4f \t %.4f \t %.4f \t %.4f", s_mmm, s_mme, s_eem, s_eee, s) << std::endl;
  yieldsFile << "s / b               " << Form("%.4f \t %.4f \t %.4f \t %.4f \t %.4f", s_mmm/b_mmm, s_mme/b_mme, s_eem/b_eem, s_eee/b_eee, s/b) << std::endl;
  yieldsFile << "s / sqrt(b)         " << Form("%.4f \t %.4f \t %.4f \t %.4f \t %.4f", s_mmm/sqrt(b_mmm), s_mme/sqrt(b_mme), s_eem/sqrt(b_eem), s_eee/sqrt(b_eee), s/sqrt(b)) << std::endl;

  yieldsFile.close();

 
  for( unsigned i=0; i<db->get_mcFiles().size(); ++i) {

    int inverseIndex = db->get_mcFiles().size()-i-1;
    if( vh1_yields_mc[inverseIndex]->Integral()>0 )
      legend->AddEntry( vh1_yields_mc[inverseIndex], db->get_mcFile(i).legendName.c_str(), "F" );

  } 






  float yMax = stackMC->GetMaximum();
  yMax *= 2.2;

  TH2D* h2_axes = new TH2D("axes", "", 4, 0., 4., 10, 0., yMax);
  h2_axes->GetXaxis()->SetLabelSize(0.085);
  h2_axes->GetXaxis()->SetBinLabel(1, "(#mu#mu)#mu");
  h2_axes->GetXaxis()->SetBinLabel(2, "(#mu#mu)e");
  h2_axes->GetXaxis()->SetBinLabel(3, "(ee)#mu");
  h2_axes->GetXaxis()->SetBinLabel(4, "(ee)e");
  h2_axes->SetYTitle("Events");


  TPaveText* label_sqrt = db->get_labelSqrt();

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  h2_axes->Draw();
  stackMC->Draw("histo same");
  legend->Draw("same");
  label_sqrt->Draw("same");

  gPad->RedrawAxis();
  
  char canvasName[500];
  if( saveName!="" )
    sprintf( canvasName, "%s/channelYields_%s.eps", db->get_outputdir().c_str(), saveName.c_str() );
  else
    sprintf( canvasName, "%s/channelYields.eps", db->get_outputdir().c_str() );
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;

  for( unsigned i=0; i<vh1_yields_mc.size(); ++i ) 
    delete vh1_yields_mc[vh1_yields_mc.size()-i-1];

}

