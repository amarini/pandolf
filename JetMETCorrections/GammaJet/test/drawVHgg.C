#include <stdlib.h>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"


//void drawHistoWithShape( std::vector<TH1D*> refHistos, const std::string& plotName, DrawBase* db );
std::vector<TH1D*> get_normalizedHistos( std::vector<TH1D*> referenceHistos, std::vector<TH1D*> lastHistos, float luminorm );
void printYields( std::vector<TH1D*> histos );



int main(int argc, char* argv[]) {

//if( argc != 2 && argc != 3 ) {
//  std::cout << "USAGE: ./drawVHgg [LUMI_SHAPE] [photonID=\"medium\"]" << std::endl;
//  exit(23);
//}

  



  DrawBase* db = new DrawBase("VHgg");


  std::string mcSignalFileName = "VHgg_WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root";
  TFile* mcSignalFile = TFile::Open(mcSignalFileName.c_str());
  db->add_mcFile( mcSignalFile, "WH_ZH_HToGG_M-120_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1", "FP VH (120)", 46);

  //std::string mcBGFileName = "VHgg_DiPhotonBoxBorn.root";
  //TFile* mcBGFile = TFile::Open(mcBGFileName.c_str());
  //db->add_mcFile( mcBGFile, "DiPhotonBoxBorn", "SM DiPhoton", 38);

  std::string mcBGFileNameBorn = "VHgg_DiPhotonJets_7TeV-madgraph-Fall11-PU_S6_START42_V14B-v1.root";
  TFile* mcBGFileBorn = TFile::Open(mcBGFileNameBorn.c_str());
  db->add_mcFile( mcBGFileBorn, "DiPhotonBorn", "#gamma#gamma (Born)", 38);

  std::string mcBGFileNameBox = "VHgg_DiPhotonBox_Pt-25To250_7TeV-pythia6-Fall11-PU_S6_START42_V14B-v1.root";
  TFile* mcBGFileBox = TFile::Open(mcBGFileNameBox.c_str());
  db->add_mcFile( mcBGFileBox, "DiPhotonBox", "#gamma#gamma (Box)", 29);

  db->set_shapeNormalization();

  db->set_outputdir();


  bool log = true;

  db->drawHisto( "nvertex", "Number of Reconstructed Vertexes", "", "Events", log);

  db->drawHisto( "pdgIdJet1", "Lead Jet Parton PDG ID", "", "Events", log);
  db->drawHisto( "pdgIdJet2", "Sublead Jet Parton PDG ID", "", "Events", log);


  db->set_rebin(5);
  db->set_xAxisMin(30.);
  db->set_xAxisMax(200.);
  db->drawHisto( "ptPhot1", "Lead Photon p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(100.);
  db->drawHisto( "ptPhot2", "Sublead Photon p_{T}", "GeV", "Events", log);
  db->set_xAxisMin(20.);
  db->set_xAxisMax(200.);
  db->drawHisto( "ptJet1", "Lead Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMax(100.);
  db->drawHisto( "ptJet2", "Sublead Jet p_{T}", "GeV", "Events", log);
  db->set_xAxisMin(50.);
  db->drawHisto( "mjj", "DiJet Invariant Mass", "GeV", "Events", log);

  db->set_xAxisMin();
  db->set_xAxisMax();
  db->drawHisto( "etaPhot1", "Lead Photon #eta", "", "Events", log);
  db->drawHisto( "etaPhot2", "Sublead Photon #eta", "", "Events", log);

  db->drawHisto( "etaJet1", "Lead Jet #eta", "", "Events", log);
  db->drawHisto( "etaJet2", "Sublead Jet #eta", "", "Events", log);

  db->drawHisto( "QGLikelihoodJet1", "Lead Jet QG Likelihood", "", "Events", log);
  db->drawHisto( "QGLikelihoodJet2", "Sublead Jet QG Likelihood", "", "Events", log);
  db->drawHisto( "QGLikelihoodProd", "QG Likelihood Product", "", "Events", log);

  db->set_rebin();
  db->drawHisto( "kinfit_chiSquare", "Kinematic Fit #chi^{2}", "", "Events", log);
  //db->set_xAxisMax(0.01);
  db->drawHisto( "kinfit_chiSquareProb", "Kinematic Fit #chi^{2} Probability", "", "Events", log);
  db->drawHisto( "kinfit_chiSquareProbMax", "Kinematic Fit #chi^{2} Probability", "", "Events", log);
  db->set_xAxisMax();

  db->drawHisto_fromTree("tree_passedEvents", "etaJet1-etaJet2", "eventWeight", 100, -5., 5., "deltaEtajj", "#Delta#eta(jet-jet)", "");
  db->drawHisto_fromTree("tree_passedEvents", "chiSquareProbMax", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5)", 20, 0., 0.01, "chiSquareProbMax_fanelli", "Kinematic Fit #chi^{2} Probability", "GeV");

  float luminorm = 10000.;
  db->set_lumiNormalization(luminorm);
  db->set_noStack(false);


  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight", 100, 50., 550., "mgg", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>60.&&mjj<110.)", 100, 50., 550., "mgg_mjj", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>60.&&mjj<110.&&ptJet1>40.)", 100, 50., 550., "mgg_mjj_ptj", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>60.&&mjj<110.&&ptJet1>40.&&chiSquareProbMax>0.01)", 100, 50., 550., "mgg_mjj_ptj_chiSquareProb", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>60.&&mjj<110.&&ptJet1>40.&&chiSquareProbMax>0.01&&QGLikelihoodJet1>0.9)", 100, 50., 550., "mgg_mjj_ptj_chiSquareProb_QG1", "m_{#gamma#gamma}", "GeV");


  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5 && QGLikelihoodJet1>0.95)", 20, 90., 190., "mgg_fanelli_QG1", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5 && chiSquareProbMax>0.001)", 20, 90., 190., "mgg_fanelli_chiSquareProb", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5 && QGLikelihoodJet1>0.95 && chiSquareProbMax>0.001)", 20, 90., 190., "mgg_fanelli_QG1_chiSquareProb", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5 && QGLikelihoodJet1*QGLikelihoodJet2>0.8)", 20, 90., 190., "mgg_fanelli_QGP", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&chiSquareProbMax>0.1)", 20, 90., 190., "mgg_fanelli_chiSquareMaxInsteadOfMjj", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mjj", "eventWeight*(chiSquareProbMax>0.001)", 100, 30., 230., "mjj_chiSquareMax", "DiJet Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mjj", "eventWeight*(chiSquareProbMax<0.001)", 100, 30., 230., "mjj_chiSquareMaxInvert", "DiJet Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mjj", "eventWeight*(chiSquareProbMax<0.01)", 100, 30., 230., "mjj_chiSquareMax1Invert", "DiJet Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mjj", "eventWeight*(chiSquareProbMax<0.1)", 100, 30., 230., "mjj_chiSquareMax10Invert", "DiJet Mass", "GeV");


  db->set_legendTitle("");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&simpleSecondaryVertexHighEfficiencyJet1>0.&&simpleSecondaryVertexHighEfficiencyJet2>0.)", 20, 90., 190., "mgg_fanelli_nomjj_btag", "m_{#gamma#gamma}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5)", 20, 90., 190., "mgg_fanelli", "m_{#gamma#gamma}", "GeV");
  std::vector<TH1D*> lastHistos = db->get_lastHistos_mc();
  std::vector<TH1D*> referenceHistos;
  for(unsigned int iHisto=0; iHisto<lastHistos.size(); ++iHisto ) {
    TH1D* newHisto = new TH1D(*(lastHistos[iHisto]));
    referenceHistos.push_back(newHisto);
  }
  std::vector<TH1D*> normalizedHistos = get_normalizedHistos( referenceHistos, lastHistos, luminorm );
  printYields( normalizedHistos );
  db->drawHisto_fromHistos( db->get_lastHistos_data(), normalizedHistos, db->get_lastHistos_mc_superimp(), "mgg", "m_{#gamma#gamma}", "GeV", "Entries", false, 1, "fanelli_shape");
  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&((trackCountingHighEfficiencyJet1>3.3&&trackCountingHighEfficiencyJet2>1.7)||(trackCountingHighEfficiencyJet1>1.7&&trackCountingHighEfficiencyJet2>3.3)))", 20, 90., 190., "mgg_fanelli_btagTC_medloose", "m_{#gamma#gamma}", "GeV");
  printYields( normalizedHistos );
  lastHistos = db->get_lastHistos_mc();
  normalizedHistos = get_normalizedHistos( referenceHistos, lastHistos, luminorm );
  db->drawHisto_fromHistos( db->get_lastHistos_data(), normalizedHistos, db->get_lastHistos_mc_superimp(), "mgg", "m_{#gamma#gamma}", "GeV", "Entries", false, 1, "fanelli_btagTC_medloose_shape");
  printYields( normalizedHistos );

  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&trackCountingHighEfficiencyJet1>3.3&&trackCountingHighEfficiencyJet2>3.3)", 20, 90., 190., "mgg_fanelli_btagTC_medmed", "m_{#gamma#gamma}", "GeV");
  printYields( normalizedHistos );
  lastHistos = db->get_lastHistos_mc();
  normalizedHistos = get_normalizedHistos( referenceHistos, lastHistos, luminorm );
  db->drawHisto_fromHistos( db->get_lastHistos_data(), normalizedHistos, db->get_lastHistos_mc_superimp(), "mgg", "m_{#gamma#gamma}", "GeV", "Entries", false, 1, "fanelli_btagTC_medmed_shape");
  printYields( normalizedHistos );

  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&simpleSecondaryVertexHighEfficiencyJet1>1.74&&simpleSecondaryVertexHighEfficiencyJet2>1.74)", 20, 90., 190., "mgg_fanelli_btag", "m_{#gamma#gamma}", "GeV");
  printYields( normalizedHistos );
  lastHistos = db->get_lastHistos_mc();
  normalizedHistos = get_normalizedHistos( referenceHistos, lastHistos, luminorm );
  db->drawHisto_fromHistos( db->get_lastHistos_data(), normalizedHistos, db->get_lastHistos_mc_superimp(), "mgg", "m_{#gamma#gamma}", "GeV", "Entries", false, 1, "fanelli_btag_shape");
  printYields( normalizedHistos );

  db->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&trackCountingHighEfficiencyJet1>1.7&&trackCountingHighEfficiencyJet2>1.7)", 20, 90., 190., "mgg_fanelli_btagTC_looseloose", "m_{#gamma#gamma}", "GeV");
  lastHistos = db->get_lastHistos_mc();
  normalizedHistos = get_normalizedHistos( referenceHistos, lastHistos, luminorm );
  db->drawHisto_fromHistos( db->get_lastHistos_data(), normalizedHistos, db->get_lastHistos_mc_superimp(), "mgg", "m_{#gamma#gamma}", "GeV", "Entries", false, 1, "fanelli_btagTC_looseloose_shape");
  printYields( normalizedHistos );



  db->set_shapeNormalization();
  db->set_noStack(true);
  db->set_legendTitle("115 < m_{#gamma#gamma} < 125 GeV");
  db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet1", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&mgg>115.&&mgg<125.)", 20, 0., 1.0001, "QGLikelihoodJet1_fanelli", "Lead Jet QGLikelihood", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet2", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&mgg>115.&&mgg<125.)", 20, 0., 1.0001, "QGLikelihoodJet2_fanelli", "Sublead Jet QGLikelihood", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "QGLikelihoodJet1*QGLikelihoodJet2", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&mgg>115.&&mgg<125.)", 20, 0., 1.0001, "QGLikelihoodProd_fanelli", "QGLikelihoodProduct" "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "chiSquareProb", "eventWeight*(mjj>55.&&mjj<115.&&ptJet1>35.&&ptPhot1>60.&&abs(etaJet1-etaJet2)<2.5&&mgg>115.&&mgg<125.)", 20, 0., 1.0001, "chiSquareProb_fanelli", "Kinematic Fit #chi^{2} Probability", "GeV");




  delete db;
  db = 0;

  return 0;

}  





std::vector<TH1D*> get_normalizedHistos( std::vector<TH1D*> referenceHistos, std::vector<TH1D*> lastHistos, float luminorm ) {

  std::vector<TH1D*> returnHistos;

  for( unsigned iHisto=0; iHisto<referenceHistos.size(); ++iHisto ) {
  
    TH1D* newHisto = new TH1D(*(referenceHistos[iHisto]));

    if( newHisto->Integral()!=0. )
      newHisto->Scale( 1./newHisto->Integral() );
    newHisto->Scale( lastHistos[iHisto]->Integral()/luminorm );

    returnHistos.push_back(newHisto);

  }

  return returnHistos;

}


void printYields( std::vector<TH1D*> histos ) {

  float xMin = 115.;
  float xMax = 125.;

  int binXmin = histos[0]->FindBin(xMin);
  int binXmax = histos[0]->FindBin(xMax) -1;

  std::cout << std::endl << "Yields (@ 10 fb-1): " << std::endl;
  for( unsigned int ii=0; ii<histos.size(); ++ii )
    std::cout << histos[ii]->Integral(binXmin, binXmax) << std::endl;


}
