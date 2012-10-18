#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawTTW [(string)selType] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);

  std::string bTaggerType = "JP";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }



  float lumi_fb = 4.6;

  DrawBase* db = new DrawBase("TTW");


  std::string outputdir_str = "TTWPlots_MConly_" + selType + "_" + bTaggerType;
  db->set_outputdir(outputdir_str);


  std::string mcTTHFileName = "TTW_TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1";
  mcTTHFileName += "_" + selType;
  mcTTHFileName += "_" + bTaggerType;
  mcTTHFileName += ".root";
  TFile* mcTTHFile = TFile::Open(mcTTHFileName.c_str());
  db->add_mcFile( mcTTHFile, "TTH", "t#bar{t} + H", kBlack, 0);



  std::string mcTTWFileName = "TTW_TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1";
  mcTTWFileName += "_" + selType;
  mcTTWFileName += "_" + bTaggerType;
  mcTTWFileName += ".root";
  TFile* mcTTWFile = TFile::Open(mcTTWFileName.c_str());
  db->add_mcFile( mcTTWFile, "ttW", "t#bar{t} + W", 46, 3005);



  std::string mcTTZFileName = "TTW_TTZJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1";
  mcTTZFileName += "_" + selType;
  mcTTZFileName += "_" + bTaggerType;
  mcTTZFileName += ".root";
  TFile* mcTTZFile = TFile::Open(mcTTZFileName.c_str());
  db->add_mcFile( mcTTZFile, "TTZ", "t#bar{t} + Z", 38, 3004);





  db->set_shapeNormalization();


  bool log = true;

  db->drawHisto_fromTree( "tree_passedEvents", "nnbjets", "", 9, -0.5, 8.5, "nnbjets", "Non b-Tagged Jet Multiplicity (p_{T} > 20 GeV)", "", "Events", log);
  db->drawHisto_fromTree( "tree_passedEvents", "nbjets", "", 9, -0.5, 8.5, "nbjets", "b-Tagged Jet Multiplicity (p_{T} > 20 GeV)", "", "Events", log);
  db->drawHisto_fromTree( "tree_passedEvents", "nnbjets+nbjets", "", 9, -0.5, 8.5, "nnbjets+nbjets", "Jet Multiplicity (p_{T} > 20 GeV)", "", "Events", log);


  delete db;
  db = 0;

  return 0;

}  


