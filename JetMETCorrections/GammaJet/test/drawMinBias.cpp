#include <stdlib.h>
#include "DrawBase.h"




//void drawMinBias( std::string dijet_inclusive, std::string data_dataset, std::string mc_dataset, std::string algoType, Int_t pt_thresh=5, std::string raw_corr="raw", bool useID=true) {

int main(int argc, char* argv[]) {

  if( argc != 8 && argc != 9 ) {
    std::cout << "USAGE: ./drawMinBias [(string)DiJet/Inclusive] [(string)data_dataset] [(string)mc_dataset] [(string)recoType] [(string)jetAlgo] [(int)pt_thresh] [(string)raw/corr] [(bool)useID=true]" << std::endl;
    exit(23);
  }

  std::string dijet_inclusive( argv[1] );
  std::string data_dataset( argv[2] );
  std::string mc_dataset( argv[3] );
  std::string recoType( argv[4] );
  std::string jetAlgo( argv[5] );
  int pt_thresh = atoi( argv[6] );
  std::string raw_corr( argv[7] );
  Float_t etamax = 3.;

  std::string algoType = (recoType=="calo") ? jetAlgo : recoType+jetAlgo;

  bool useID = true;
  if( argc==9 ) {
    std::string bool_str( argv[8] );
    if( bool_str=="false" ) useID = false;
  }

  bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  DrawBase* db = new DrawBase("MinBias", "pf", "akt5");
  db->set_raw_corr( raw_corr );
  db->set_pt_thresh( pt_thresh );
  //db->set_etamax( etamax );
  db->set_etamax(3.);
  db->set_pdf_aussi((bool)false);

  char outputdir_char[200];
  sprintf( outputdir_char, "%sPlots_%s_vs_%s_%s_%dGeV%s_eta%d", dijet_inclusive.c_str(), data_dataset.c_str(), mc_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*db->get_etamax()));
  std::string outputdir_str(outputdir_char);
  if( !useID )
    outputdir_str = outputdir_str + "_noID";

  db->set_outputdir(outputdir_str);

  //std::string dataFileName = "MinBias_" + data_dataset + "_" + algoType;
  char dataFileName[150];
  if( useID )
    sprintf( dataFileName, "%s_%s_%s_%dGeV%s_eta%d.root", dijet_inclusive.c_str(), data_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*db->get_etamax()));
    //sprintf( dataFileName, "%s_DATA_%s_%s_%dGeV%s_eta%d.root", data_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*ETAMAX));
  else 
    sprintf( dataFileName, "%s_%s_%s_%dGeV%s_eta%d_noID.root", dijet_inclusive.c_str(), data_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*db->get_etamax()));
    //sprintf( dataFileName, "%s_DATA_%s_%s_%dGeV%s_eta%d_noDeltaPhi.root", data_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*ETAMAX));
  TFile* dataFile = TFile::Open(dataFileName);
  if( dataFile!=0 )
    std::cout << "-> Opened data file '" << dataFileName << "'." << std::endl;
  else
    std::cout << "-> ERROR opening data file '" << dataFileName << "'." << std::endl;

  db->set_dataFile( dataFile );

  char mcFileName[150];
//if( mc_dataset=="" ) {
//  if( useID )
//    sprintf( mcFileName, "%s_%s_%s_%dGeV%s_eta%d.root", dijet_inclusive.c_str(), data_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*ETAMAX));
//  else
//    sprintf( mcFileName, "%s_%s_%s_%dGeV%s_eta%d_noID.root", dijet_inclusive.c_str(), data_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*ETAMAX));
//} else {
    if( useID )
      sprintf( mcFileName, "%s_%s_%s_%dGeV%s_eta%d.root", dijet_inclusive.c_str(), mc_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*db->get_etamax()));
    else
      sprintf( mcFileName, "%s_%s_%s_%dGeV%s_eta%d_noID.root", dijet_inclusive.c_str(), mc_dataset.c_str(), algoType.c_str(), pt_thresh, raw_corr.c_str(), (int)(10.*db->get_etamax()));
//}
  TFile* mcFile = TFile::Open(mcFileName);
  if( mcFile!=0 )
    std::cout << "-> Opened mc file '" << mcFileName << "'." << std::endl;
  else
    std::cout << "-> ERROR opening mc file '" << mcFileName << "'." << std::endl;

  db->set_mcFile( mcFile );

//  db->set_crossSectionNormalization();

  db->set_shapeNormalization();

  //db->set_scaleFactor(0.000565489);
  std::cout << "-> Normalized." << std::endl;

//TH1F* h1_nJets_data = (TH1F*)dataFile->Get("nJets");
//TH1F* h1_nJets_mc = (TH1F*)mcFile->Get("nJets");

//TH1F* h1_phiJet_data = (TH1F*)dataFile->Get("phiJet");
//TH1F* h1_phiJet_mc = (TH1F*)mcFile->Get("phiJet");


//Float_t scaleFactor = ( sameEvents ) ? h1_nJets_data->Integral()/h1_nJets_mc->Integral() : h1_phiJet_data->Integral()/h1_phiJet_mc->Integral();

  bool log = true;

  if( dijet_inclusive == "DiJet" ) {
    db->drawHisto( "deltaPhiJet", "", "" ,0);
    db->drawHisto( "asymmJet", "", "" );
    db->drawHisto( "diJetMass", "", "", 1, log);
  }

  db->drawHisto( "massJet", "", "" ,1, log);
  db->drawHisto( "MoEJet", "", "", 1, log);
  db->drawHisto( "ptOverMJet", "", "", 1, log);
  db->drawHisto( "pOverMJet", "", "", 1, log);

  db->drawHisto( "ptJet", "", "", 1, log);
  db->drawHisto( "ptCorrJet", "", "" ,1, log);
  db->drawHisto( "ptCorrJet", "barrel", "", 1, log);
  db->drawHisto( "ptCorrJet", "endcap", "", 1, log);
  db->drawHisto( "ptCorrJet", "eta1430", "", 1, log);
  db->drawHisto( "ptCorrJet", "Rch050", "" ,1, log);
  db->drawHisto( "ptCorrJet", "Rch5070", "" ,1, log);
  db->drawHisto( "ptCorrJet", "Rch70100", "" ,1, log);
  db->drawHisto( "phiJet", "", "");
  db->drawHisto( "etaJet", "", "");

  db->drawHisto( "nCandJet", "", "", log);

  db->drawHisto( "RchJet", "", "", 1);
  db->drawHisto( "RchJet", "eta02", "" ,1);
  db->drawHisto( "RchJet", "eta23", "" ,1);
  db->drawHisto( "RchJet", "barrel", "" ,1);
  db->drawHisto( "RchJet", "eta163", "" ,1);
  db->drawHisto( "RchJet", "eta1425", "" ,1);
  db->drawHisto( "RnhJet", "", "", 1, log);
  db->drawHisto( "RnhJet", "eta02", "", 1, log);
  db->drawHisto( "RnhJet", "eta23", "", 1, log);
  db->drawHisto( "RnhJet", "barrel", "", 1, log);
  db->drawHisto( "RnhJet", "eta163", "", 1, log);
  db->drawHisto( "RnhJet", "eta1425", "", 1, log);
  db->drawHisto( "RgammaJet", "", "", 1, log);
  db->drawHisto( "RgammaJet", "eta02", "" ,1, log);
  db->drawHisto( "RgammaJet", "eta23", "" ,1, log);
  db->drawHisto( "RgammaJet", "barrel", "" ,1, log);
  db->drawHisto( "RgammaJet", "eta163", "" ,1, log);
  db->drawHisto( "RgammaJet", "eta1425", "" ,1, log);
  db->drawHisto( "RgammanhJet", "barrel", "" ,1, log);
  db->drawHisto( "RgammanhJet", "eta163", "" ,1, log);
  db->drawHisto( "RgammanhJet", "eta1425", "" ,1, log);

  db->drawHisto( "EchJet", "", "", 1, log);
  db->drawHisto( "EchJet", "eta02", "" ,1, log);
  db->drawHisto( "EchJet", "eta23", "" ,1, log);
  db->drawHisto( "EchJet", "barrel", "" ,1, log);
  db->drawHisto( "EchJet", "eta163", "" ,1, log);
  db->drawHisto( "EchJet", "eta1425", "" ,1, log);
  db->drawHisto( "EnhJet", "", "", 1, log);
  db->drawHisto( "EnhJet", "eta02", "" ,1, log);
  db->drawHisto( "EnhJet", "eta23", "" ,1, log);
  db->drawHisto( "EnhJet", "barrel", "" ,1, log);
  db->drawHisto( "EnhJet", "eta163", "" ,1, log);
  db->drawHisto( "EnhJet", "eta1425", "" ,1, log);
  db->drawHisto( "EgammaJet", "", "" ,1, log);
  db->drawHisto( "EgammaJet", "eta02", "" ,1, log);
  db->drawHisto( "EgammaJet", "eta23", "" ,1, log);
  db->drawHisto( "EgammaJet", "barrel", "" ,1, log);
  db->drawHisto( "EgammaJet", "eta163", "" ,1, log);
  db->drawHisto( "EgammaJet", "eta1425", "" ,1, log);
  db->drawHisto( "EgammanhJet", "barrel", "" ,1, log);
  db->drawHisto( "EgammanhJet", "eta163", "" ,1, log);
  db->drawHisto( "EgammanhJet", "eta1425", "" ,1, log);

  db->drawHisto( "PTchJet", "", "" ,1, log);
  db->drawHisto( "PTgammaJet", "", "" ,1, log);
  db->drawHisto( "PTnhJet", "", "", 1, log);

  db->drawHisto( "EphotAveJet", "", "" ,1, log);

  db->drawHisto( "NchJet", "", "", 1, log);
  db->drawHisto( "NchJet", "barrel", "", 1, log);
  db->drawHisto( "NchJet", "eta163", "", 1, log);
  db->drawHisto( "NchJet", "eta1425", "", 1, log);
  db->drawHisto( "NchJet", "eta02", "", 1, log);
  db->drawHisto( "NchJet", "eta23", "", 1, log);
  db->drawHisto( "NnhJet", "", "", 1, log);
  db->drawHisto( "NnhJet", "eta02", "" ,1, log);
  db->drawHisto( "NnhJet", "eta23", "" ,1, log);
  db->drawHisto( "NnhJet", "barrel", "" ,1, log);
  db->drawHisto( "NnhJet", "eta163", "" ,1, log);
  db->drawHisto( "NnhJet", "eta1425", "" ,1, log);
  db->drawHisto( "NgammaJet", "", "", 1, log);
  db->drawHisto( "NgammaJet", "eta02", "" ,1, log);
  db->drawHisto( "NgammaJet", "eta23", "" ,1, log);
  db->drawHisto( "NgammaJet", "barrel", "" ,1, log);
  db->drawHisto( "NgammaJet", "eta163", "" ,1, log);
  db->drawHisto( "NgammaJet", "eta1425", "" ,1, log);
  db->drawHisto( "NgammanhJet", "barrel", "" ,1, log);
  db->drawHisto( "NgammanhJet", "eta163", "" ,1, log);
  db->drawHisto( "NgammanhJet", "eta1425", "" ,1, log);

  std::string ptVar;
  if( raw_corr=="raw" )
    ptVar = "pt";
  else 
    ptVar = "ptCorr";

  db->drawProfile( "pt", "eta", 1);
  db->drawProfile( "ptCorr", "eta", 1);

  db->drawProfile( "Rch", "eta", 1);
  db->drawProfile( "Rgamma", "eta", 1);
  db->drawProfile( "Rnh", "eta", 1);
  db->drawProfile( "Rgammanh", "eta", 1);

  db->drawProfile( "Ech", "eta", 1);
  db->drawProfile( "Egamma", "eta", 1);
  db->drawProfile( "Enh", "eta", 1);
  db->drawProfile( "Egammanh", "eta", 1);

  db->drawProfile( "Ech", "pt", 1);
  db->drawProfile( "Egamma", "pt", 1);
  db->drawProfile( "Enh", "pt", 1);

  db->drawProfile( "PTch", "eta", 1);
  db->drawProfile( "PTgamma", "eta", 1);
  db->drawProfile( "PTnh", "eta", 1);

  db->drawProfile( "Nch", "eta", 1);
  db->drawProfile( "Ngamma", "eta", 1);
  db->drawProfile( "Nnh", "eta", 1);
  db->drawProfile( "Ngammanh", "eta", 1);

  db->drawProfile( "Rch", ptVar, 1);
  db->drawProfile( "Rgamma", ptVar, 1);
  db->drawProfile( "Rnh", ptVar, 1);

  db->drawProfile( "Ech", ptVar, 1);
  db->drawProfile( "Egamma", ptVar, 1);
  db->drawProfile( "Enh", ptVar, 1);

  db->drawProfile( "Nch", ptVar, 1);
  db->drawProfile( "Ngamma", ptVar, 1);
  db->drawProfile( "Nnh", ptVar, 1);

  db->drawProfile( "Egamma", "Ech", 1);
  db->drawProfile( "Enh",    "Ech", 1);

  db->drawProfile( "Rgamma", "Rch", 1);
  db->drawProfile( "Rnh",    "Rch", 1);
  db->drawProfile( "Rgamma",    "Rnh", 1);

//drawH2( dataFile, mcFile, "Rgamma", "Rch", 1);
//drawH2( dataFile, mcFile, "Rnh",    "Rch", 1);
//drawH2( dataFile, mcFile, "Rnh",    "Rgamma", 1);

  db->drawProfile( "Ngamma", "Nch", 1);
  db->drawProfile( "Nnh",    "Nch", 1);

  // stacks:
  db->drawStack( "R", "eta", "RECO", (bool)true);
  db->drawStack( "R", "eta", "RECO", (bool)false);

  db->drawStack( "R", ptVar, "barrel", "RECO", (bool)true);
  db->drawStack( "R", ptVar, "barrel", "RECO", (bool)false);

  db->drawStack( "R", "phi", "barrel", "RECO", (bool)true);
  db->drawStack( "R", "phi", "barrel", "RECO", (bool)false);

  db->drawStack( "R", "eta", "GEN", (bool)false);
  db->drawStack( "R", "pt", "barrel", "GEN", (bool)false);
  db->drawStack( "R", "phi", "barrel", "GEN", (bool)false);

  return 0;

}  


