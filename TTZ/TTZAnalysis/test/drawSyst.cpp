#include <iostream>
#include <string>

#include "DrawBase.h"


void drawSingleSyst( DrawBase* db, const std::string& syst, const std::string& sel );


int main( int argc, char* argv[] ) {

  std::string selection="optsel3";
  std::string bTaggerType="TCHE";
  
  
  DrawBase* db = new DrawBase("TTZTrilepton");

  db->set_lumiOnRightSide(true);
  db->set_lumiNormalization(4980.);
  db->set_noStack();

  std::string outputdir_str = "TTZTrileptonPlots_DATA_Run2011_FULL_" + selection + "_" + bTaggerType + "_ALL";
  db->set_outputdir(outputdir_str);

  drawSingleSyst( db, "BTag", selection );
  drawSingleSyst( db, "JES", selection );
  drawSingleSyst( db, "JER", selection );

  return 0;

}




void drawSingleSyst( DrawBase* db, const std::string& syst, const std::string& sel ) {

  std::string systFile = "TTZTrilepton_BG_" + sel + "_TCHE_ALL.root";
  std::string systFileUP = "TTZTrilepton_BG_" + sel + "_TCHE_ALL_" + syst + "UP.root";
  std::string systFileDOWN = "TTZTrilepton_BG_" + sel + "_TCHE_ALL_" + syst + "DOWN.root";

  TFile* file_systFile = TFile::Open( systFile.c_str() );
  TFile* file_systFileUP = TFile::Open( systFileUP.c_str() );
  TFile* file_systFileDOWN = TFile::Open( systFileDOWN.c_str() );

  std::string systUP_text = syst + " + 1 #sigma";
  std::string systDOWN_text = syst + " - 1 #sigma";


  DrawBase* newdb = new DrawBase( *db );

  std::string systFlagName = syst + "Syst";
  newdb->set_flags( systFlagName );

  std::string legendTitle = syst + " Systematic";
  newdb->set_legendTitle( legendTitle );

  newdb->add_mcFile( file_systFile, "mean", "Mean", kBlack, 0);
  newdb->add_mcFile( file_systFileUP, "systUP", systUP_text, 30, 0 );
  if( file_systFileDOWN!=0 ) //JER has no down
    newdb->add_mcFile( file_systFileDOWN, "systDOWN", systDOWN_text, 50, 0 );

  newdb->drawHisto( "nJets_presel", "Jet Multiplicity", "", "Events");
  newdb->drawHisto( "nBJets_loose_presel", "b-Jet Multiplicity (Loose)", "", "Events");
  newdb->drawHisto( "nBJets_medium_presel", "b-Jet Multiplicity (Medium)", "", "Events");
  newdb->drawHisto( "channelYields", "", "", "Events", false, 2);


}
