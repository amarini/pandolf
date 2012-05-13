#include <iostream>
#include <string>
#include <cstdlib>

#include "DrawBase.h"


void drawSingleSyst( DrawBase* db, const std::string& syst, const std::string& sel, const std::string& bg_signal );


int main( int argc, char* argv[] ) {

  std::string selection="optsel3";
  std::string bTaggerType="TCHE";
  
  
  DrawBase* db = new DrawBase("TTZTrilepton");

  db->set_lumiOnRightSide(true);
  db->set_lumiNormalization(4980.);
  db->set_noStack();

  std::string outputdir_str = "TTZTrileptonPlots_DATA_Run2011_FULL_" + selection + "_" + bTaggerType + "_ALL";
  db->set_outputdir(outputdir_str);

  drawSingleSyst( db, "Lept", selection, "BG" );
  drawSingleSyst( db, "BTag", selection, "BG" );
  drawSingleSyst( db, "JES" , selection, "BG" );
  drawSingleSyst( db, "JER" , selection, "BG" );

  drawSingleSyst( db, "Lept", selection, "Signal" );
  drawSingleSyst( db, "BTag", selection, "Signal" );
  drawSingleSyst( db, "JES" , selection, "Signal" );
  drawSingleSyst( db, "JER" , selection, "Signal" );

  return 0;

}




void drawSingleSyst( DrawBase* db, const std::string& syst, const std::string& sel, const std::string& bg_signal ) {

  if( bg_signal!="BG" && bg_signal!="Signal" ) {
    std::cout << "bg_signal must be either 'BG' or 'Signal'. Exiting." << std::endl;
    exit(73);
  }

  std::string dataset = (bg_signal=="BG") ? bg_signal : "TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi";
  
  std::string systFile = "TTZTrilepton_" + dataset + "_" + sel + "_TCHE_ALL.root";
  std::string systFileUP = "TTZTrilepton_" + dataset + "_" + sel + "_TCHE_ALL_" + syst + "UP.root";
  std::string systFileDOWN = "TTZTrilepton_" + dataset + "_" + sel + "_TCHE_ALL_" + syst + "DOWN.root";

  TFile* file_systFile = TFile::Open( systFile.c_str() );
  TFile* file_systFileUP = TFile::Open( systFileUP.c_str() );
  TFile* file_systFileDOWN = TFile::Open( systFileDOWN.c_str() );

  std::string systUP_text = syst + " + 1 #sigma";
  std::string systDOWN_text = syst + " - 1 #sigma";


  DrawBase* newdb = new DrawBase( *db );

  std::string systFlagName = bg_signal + "Only" + syst + "Syst";
  newdb->set_flags( systFlagName );

  std::string legendTitle = syst + " Systematic";
  newdb->set_legendTitle( legendTitle );

  newdb->add_mcFile( file_systFile, "mean", "Mean", kBlack, 0);
  newdb->add_mcFile( file_systFileUP, "systUP", systUP_text, 30, 0 );
  if( file_systFileDOWN!=0 ) //JER has no down
    newdb->add_mcFile( file_systFileDOWN, "systDOWN", systDOWN_text, 50, 0 );

  std::string instanceName = (bg_signal=="BG") ? "Background" : bg_signal; 
  instanceName = instanceName + " Events";
  newdb->drawHisto( "nJets_presel", "Jet Multiplicity", "", instanceName );
  newdb->drawHisto( "nBJets_loose_presel", "b-Jet Multiplicity (Loose)", "", instanceName );
  newdb->drawHisto( "nBJets_medium_presel", "b-Jet Multiplicity (Medium)", "", instanceName );
  newdb->set_getBinLabels(true);
  newdb->drawHisto( "channelYields", "", "", instanceName, false, 2);
  newdb->set_getBinLabels(false);

}
