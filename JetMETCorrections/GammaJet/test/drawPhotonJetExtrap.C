#include <stdlib.h>
#include "DrawExtrap.h"
#include "fitTools.h"


bool useMCassoc_ = false;
bool NOQ=false;
bool FIXM=false;
bool EXCLUDE_FIRST_POINT=true;



int main(int argc, char* argv[]) {

  if( argc != 6 && argc!=7 && argc!=8 ) {
    std::cout << "USAGE: ./drawPhotonJetExtrap [data_dataset] [mc_dataset] [recoType] [jetAlgo] [FIT_RMS] [flags=\"\"] [GENJETS=\"\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset( argv[1] );
  std::string mc_dataset( argv[2] );
  std::string recoType( argv[3] );
  std::string jetAlgo( argv[4] );
  std::string FIT_RMS( argv[5] );
  std::string flags = "";
  if( argc==7 ) {
    std::string flags_str(argv[6]);
    flags = flags_str;
    std::cout << "flags set." << std::endl;
  }
  bool useGenJets=false;
  if( argc==8 ) {
    std::string flags_str(argv[7]);
    if( flags_str == "GENJETS" ) useGenJets=true;
  }


  std::string algoType;
  if( recoType=="calo" )
    algoType = jetAlgo;
  else
    algoType = recoType+jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5") algoType = "jptak5";
  

  std::string mcFlags = (useGenJets) ? "GENJETS_" + flags : flags;

  DrawExtrap* db = new DrawExtrap("PhotonJet", recoType, jetAlgo, mcFlags);
  db->set_pdf_aussi((bool)false);
  db->set_isCMSArticle((bool)false);

  db->set_FITRMS(FIT_RMS);

  std::string NOQtext = (NOQ) ? "_NOQ" : "";
  std::string FIXMtext = (FIXM) ? "_FIXM" : "";
  std::string NOFIRSTPtext = (EXCLUDE_FIRST_POINT) ? "_NOFIRSTP" : "";

  char outputdir_char[200];
  if( mcFlags!="" ) {
    sprintf( outputdir_char, "PhotonJetExtrapPlots_%s_vs_%s_%s_%s_%s%s%s%s", data_dataset.c_str(), mc_dataset.c_str(), algoType.c_str(), mcFlags.c_str(), FIT_RMS.c_str(), FIXMtext.c_str(), NOFIRSTPtext.c_str(), NOQtext.c_str());
  } else {
    sprintf( outputdir_char, "PhotonJetExtrapPlots_%s_vs_%s_%s_%s%s%s%s", data_dataset.c_str(), mc_dataset.c_str(), algoType.c_str(), FIT_RMS.c_str(), FIXMtext.c_str(), NOFIRSTPtext.c_str(), NOQtext.c_str());
  }
  std::string outputdir_str(outputdir_char);

  std::vector< float > ptPhot_binning = fitTools::getPtPhot_binning();

  db->set_outputdir(outputdir_str);

  char dataFileName[150];
  if( flags != "" ) {
    sprintf( dataFileName, "PhotonJet_%s_%s_%s.root", data_dataset.c_str(), algoType.c_str(), flags.c_str());
  } else {
    sprintf( dataFileName, "PhotonJet_%s_%s.root", data_dataset.c_str(), algoType.c_str());
  }
  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->add_dataFile( dataFile, data_dataset );


  char mcFileName[150];
  if( mcFlags!="" ) {
    sprintf( mcFileName, "PhotonJet_%s_%s_%s.root", mc_dataset.c_str(), algoType.c_str(), mcFlags.c_str());
  } else {
    sprintf( mcFileName, "PhotonJet_%s_%s.root", mc_dataset.c_str(), algoType.c_str());
  }
  TFile* mcFile = TFile::Open(mcFileName);
  std::cout << "Opened mc file '" << mcFileName << "'." << std::endl;

  db->add_mcFile( mcFile, mc_dataset.c_str(), "#gamma+jet + QCD MC", 46);


  TString data_dataset_tstr(data_dataset);
  float intLumi;
  if( data_dataset_tstr.Contains("2010") ) {
    intLumi = 36.;
  } else if( data_dataset_tstr.Contains("2011") ) {
    intLumi = 50.;
  } 
  db->set_lumiNormalization(intLumi);


  bool log = true;

  db->set_NOQ( NOQ );
  db->set_FIXM( FIXM );
  db->set_EXCLUDEFIRSTPOINT( EXCLUDE_FIRST_POINT );
 
  db->drawResponseExtrap("eta011");
  db->drawResponseExtrap("eta011", (bool)true);
  db->drawResponseExtrap("eta011", (bool)true, "RecoRelRaw");
std::cout << "prova" << std::endl;
  db->drawResponseExtrap("eta013");
  db->drawResponseExtrap("eta013", (bool)true);
  db->drawResponseExtrap("eta013", (bool)true, "RecoRelRaw");
//db->drawResponseExtrap("eta1524");
//db->drawResponseExtrap("eta1524", (bool)true);
//db->drawResponseExtrap("eta243");
//db->drawResponseExtrap("eta243", (bool)true);
//db->drawResponseExtrap("eta35");


  delete db;
  db = 0;

  return 0;

}  


