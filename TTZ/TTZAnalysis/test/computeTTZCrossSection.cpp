#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "RooHistError.h"

#include "CommonTools/StatTools.h"




std::pair< float, float > getSystFromString( const std::string& systString );


int main( int argc, char* argv[] ) {


  std::string sel = "optsel3";
  if( argc>1 ) {
    std::string sel_tmp(argv[1]);
    sel = sel_tmp;
  }

  std::string suffix =  sel + "_TCHE_ALL.root";

  //std::string dataFileName = = "ttZTrilepton_DATA_Run2011_FULL_" + suffix;
  //TFile* dataFile = TFile::Open(dataFileName.c_str());

  std::string ttZFileName = "TTZTrilepton_TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi_" + suffix;
  TFile* ttZFile = TFile::Open(ttZFileName.c_str());
  //TTree* signalTree = (TTree*)ttZFile->Get("tree_passedEvents");

  std::string ttWFileName = "TTZTrilepton_TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi_" + suffix;
  TFile* ttWFile = TFile::Open(ttWFileName.c_str());

  //std::string mcZJetsFile = "ttZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11_" + suffix;
  //std::string mcTTbarFile = "ttZTrilepton_TTJ_Fall11_highstat_" + suffix;
  //std::string mcVVFile = "ttZTrilepton_VV_Summer11_" + suffix;
  //std::string mcttWFile = "ttZTrilepton_ttW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi_" + suffix;
  //bgTree->Add(dyTreeName.c_str());
  //bgTree->Add(ttTreeName.c_str());
  //bgTree->Add(dibosonTreeName.c_str());


  std::string dir = "TTZTrileptonPlots_DATA_Run2011_FULL_" + sel + "_TCHE_ALL/";
  std::string yieldFileName = dir + "yields.txt";
  std::string datacardName = dir + "datacard.txt";
  ifstream yieldsFile(yieldFileName.c_str());
  ifstream datacard(datacardName.c_str());

  std::cout << "-> Opened yield file '" << yieldFileName << "'." << std::endl;
  std::string channel;
  std::string obs_str;
  std::string s_str, ttZ_str, ttW_str, b_pred_str, b_pred_err_str;
  bool found = false;
  while( yieldsFile.good() ) {
    yieldsFile >> channel >> obs_str >> s_str >> ttZ_str >> ttW_str >> b_pred_str >> b_pred_err_str;
    if( channel == "Total" ) found = true;
  }

  if( !found ) {
    std::cout << "There must be a problem. Didn't find total." << std::endl;
    exit(3431);
  } 

  int obs = atoi(obs_str.c_str());
  float s = atof(s_str.c_str());
  float ttZ = atof(ttZ_str.c_str());
  float ttW = atof(ttW_str.c_str());
  float b_pred = atof(b_pred_str.c_str());
  float b_pred_err = atof(b_pred_err_str.c_str());
   

  // go get the syst on signal from datacard:
  std::cout << "-> Opened datacard file '" << datacardName << "'." << std::endl;
  bool go = false;
  //float lumiSystDOWN = 0.;
  //float lumiSystUP = 0.;
  float totalSignalSystDOWN = 0.;
  float totalSignalSystUP = 0.;
  float totalBGSystDOWN = 0.;
  float totalBGSystUP = 0.;
  while( datacard.good() ) {
    if( !go ) {
      char line[500];
      datacard.getline( line, 500 );
      TString line_tstr(line);
      if( (line_tstr.BeginsWith("#syst")) ) go=true;
    } else {
      std::string systName, shape, systSignal, systBG;
      datacard >> systName >> shape >> systSignal >> systBG;
      std::pair<float, float> pair_systSignal = getSystFromString(systSignal);
      std::pair<float, float> pair_systBG = getSystFromString(systBG);
      float systSignalDOWN = pair_systSignal.first;
      float systSignalUP = pair_systSignal.second;
      float systBGDOWN = pair_systBG.first;
      float systBGUP = pair_systBG.second;
   
      // last line reads zeroes dont know why
      if( systSignalDOWN==0. && systSignalUP==0. && systBGDOWN==0. && systBGUP==0. ) continue;

      systSignalDOWN = (systSignalDOWN>0.) ? fabs( 1.-systSignalDOWN ) : 0.;
      systSignalUP   = (systSignalUP>0.) ? fabs( 1.-systSignalUP ) : 0.;

      systBGDOWN = (systBGDOWN>0.) ? fabs( 1.-systBGDOWN ) : 0.;
      systBGUP   = (systBGUP>0.) ? fabs( 1.-systBGUP ) : 0.;

      totalBGSystDOWN += systBGDOWN*systBGDOWN;
      totalBGSystUP += systBGUP*systBGUP;

      totalSignalSystDOWN += systSignalDOWN*systSignalDOWN;
      totalSignalSystUP += systSignalUP*systSignalUP;
      std::cout << systName << " " << systSignalDOWN << " " << systSignalUP << " " << systBGDOWN << " " << systBGUP << std::endl;

    }
  }

  if( !go ) {
    std::cout << "Didn't find systematics in datacard!" << std::endl;
    exit(135);
  }

  totalBGSystDOWN = sqrt(totalBGSystDOWN);
  totalBGSystUP   = sqrt(totalBGSystUP);

  totalSignalSystDOWN = sqrt(totalSignalSystDOWN);
  totalSignalSystUP = sqrt(totalSignalSystUP);

  std::cout << "totalBGSystDOWN : " <<  totalBGSystDOWN << std::endl;
  std::cout << "totalBGSystUP   : " <<  totalBGSystUP   << std::endl;
  std::cout << "totalSignalSystDOWN : " << totalSignalSystDOWN << std::endl;
  std::cout << "totalSignalSystUP   : " << totalSignalSystUP << std::endl;



  // stat error on observed:
  double obs_plus, obs_minus;
  RooHistError::instance().getPoissonInterval(obs,obs_minus,obs_plus,1.);
  double obs_errPlus = obs_plus-obs;
  double obs_errMinus = obs-obs_minus;
  std::cout << "observed poissonian interval: " << obs_minus << "-" << obs_plus << std::endl;


  float ZBi = StatTools::computeZBi( s+b_pred, b_pred, b_pred_err );
  float ZBi_obs = StatTools::computeZBi( obs, b_pred, b_pred_err );

  float obs_ttV = obs - b_pred;

  TH1D* h1_nCounter_ttZ = (TH1D*)ttZFile->Get("nCounter");
  float nTotal_ttZ = h1_nCounter_ttZ->GetBinContent(1);

  TH1D* h1_nCounter_ttW = (TH1D*)ttWFile->Get("nCounter");
  float nTotal_ttW = h1_nCounter_ttW->GetBinContent(1);

  float lumi_pb = 4980.;
  float crossSection_ttZ = 0.139;
  float crossSection_ttW = 0.169;
  //float crossSection_ttW = 0.1633;
  float crossSection_ttV = crossSection_ttZ+crossSection_ttW;

  float n_passed_ttZ = nTotal_ttZ*ttZ/(crossSection_ttZ*lumi_pb);
  float eff_ttZ = n_passed_ttZ/nTotal_ttZ;

  float n_passed_ttW = nTotal_ttW*ttW/(crossSection_ttW*lumi_pb);
  float eff_ttW = n_passed_ttW/nTotal_ttW;

  float eff_ttV = ( crossSection_ttZ*eff_ttZ + crossSection_ttW*eff_ttW ) / ( crossSection_ttZ + crossSection_ttW );
std::cout << "n_passed_ttZ: " << n_passed_ttZ << " (eff: " << eff_ttZ*100. << "%   /BR: " << eff_ttZ*100./(2.*0.22*0.67*0.06) << "%)" << std::endl;
std::cout << "n_passed_ttW: " << n_passed_ttW << " (eff: " << eff_ttW*100. << "%   /BR: " << eff_ttW*100./(0.22*0.22*0.67) << "%)" << std::endl;
std::cout << "eff_ttV: " << eff_ttV*100. << " %" << std::endl;

  float crossSectionObs_ttZ = obs_ttV / ( lumi_pb*eff_ttZ );
  float crossSectionObs_ttV = obs_ttV / ( lumi_pb*eff_ttV );

  float xsecErr_ttZ_stat_plus  = obs_errPlus /  ( lumi_pb*eff_ttZ );
  float xsecErr_ttZ_stat_minus = obs_errMinus / ( lumi_pb*eff_ttZ );

  float xsecErr_ttV_stat_plus  = obs_errPlus /  ( lumi_pb*eff_ttV );
  float xsecErr_ttV_stat_minus = obs_errMinus / ( lumi_pb*eff_ttV );



  // xsec = obs_ttV / ( lumi*eff )
  // signal syst are on lumi and eff
  // so err(eff) -> err(xsec) = obs_ttV/( lumi*eff*eff ) * err(eff) = xsec*err(eff)/eff
  // so err(lumi) -> err(xsec) = obs_ttV/( lumi*lumi*eff ) * err(lumi) = xsec*err(lumi)/lumi
  // total error is xsec*(quadrature sum of all relative errors)
  float xsecErr_ttZ_systSignal_plus  = crossSection_ttZ*(totalSignalSystUP);
  float xsecErr_ttZ_systSignal_minus = crossSection_ttZ*(totalSignalSystDOWN);
  
  float xsecErr_ttV_systSignal_plus  = crossSection_ttV*(totalSignalSystUP);
  float xsecErr_ttV_systSignal_minus = crossSection_ttV*(totalSignalSystDOWN);

  
  // background systematics instead affect the numerator:
  // xsec = ( obs - BG ) / ( lumi*eff )
  // so err(BG) -> err(xsec) = err(BG) / ( lumi*eff )
  float xsecErr_ttZ_systBG_plus  = totalBGSystUP   * b_pred / ( lumi_pb*eff_ttZ );
  float xsecErr_ttZ_systBG_minus = totalBGSystDOWN * b_pred / ( lumi_pb*eff_ttZ );

  float xsecErr_ttV_systBG_plus  = totalBGSystUP   * b_pred / ( lumi_pb*eff_ttV );
  float xsecErr_ttV_systBG_minus = totalBGSystDOWN * b_pred / ( lumi_pb*eff_ttV );

  float xsecErr_ttZ_syst_plus  = sqrt( xsecErr_ttZ_systSignal_plus *xsecErr_ttZ_systSignal_plus  + xsecErr_ttZ_systBG_plus *xsecErr_ttZ_systBG_plus );
  float xsecErr_ttZ_syst_minus = sqrt( xsecErr_ttZ_systSignal_minus*xsecErr_ttZ_systSignal_minus + xsecErr_ttZ_systBG_minus*xsecErr_ttZ_systBG_minus );

  float xsecErr_ttV_syst_plus  = sqrt( xsecErr_ttV_systSignal_plus *xsecErr_ttV_systSignal_plus  + xsecErr_ttV_systBG_plus *xsecErr_ttV_systBG_plus );
  float xsecErr_ttV_syst_minus = sqrt( xsecErr_ttV_systSignal_minus*xsecErr_ttV_systSignal_minus + xsecErr_ttV_systBG_minus*xsecErr_ttV_systBG_minus );

  float xsecErr_ttZ_tot_plus  = sqrt( xsecErr_ttZ_stat_plus *xsecErr_ttZ_stat_plus  + xsecErr_ttZ_syst_plus *xsecErr_ttZ_syst_plus );
  float xsecErr_ttZ_tot_minus = sqrt( xsecErr_ttZ_stat_minus*xsecErr_ttZ_stat_minus + xsecErr_ttZ_syst_minus*xsecErr_ttZ_syst_minus );

  float xsecErr_ttV_tot_plus  = sqrt( xsecErr_ttV_stat_plus *xsecErr_ttV_stat_plus  + xsecErr_ttV_syst_plus *xsecErr_ttV_syst_plus );
  float xsecErr_ttV_tot_minus = sqrt( xsecErr_ttV_stat_minus*xsecErr_ttV_stat_minus + xsecErr_ttV_syst_minus*xsecErr_ttV_syst_minus );


  std::cout << std::endl << std::endl << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
  std::cout << " Expected BG: " << b_pred << " +- " << b_pred_err << std::endl;
  std::cout << " Expected Signal: " << s << " (eff=" << eff_ttZ*100./(2.*0.06*0.22*0.67) << "%)" << std::endl;
  std::cout << " Expected B+S: " << s+b_pred << std::endl;
  std::cout << " ZBi (expected): " << ZBi << std::endl;
  std::cout << " Observed Events: " << obs << std::endl;
  std::cout << " ZBi (observed): " << ZBi_obs << std::endl;
  std::cout << " Observed Signal (BG subtracted): " << obs_ttV << std::endl;
  std::cout << " Expected ttZ Cross Section: " << crossSection_ttZ << " pb" << std::endl;
  std::cout << " Measured ttZ Cross Section: " << std::endl;
  std::cout << " " <<  crossSectionObs_ttZ << "  +" << xsecErr_ttZ_stat_plus << "/-" << xsecErr_ttZ_stat_minus << " (stat)   +" << xsecErr_ttZ_syst_plus << "/-" << xsecErr_ttZ_syst_minus << " (syst)  pb" << std::endl;
  std::cout << " " <<  crossSectionObs_ttZ << "  +" << xsecErr_ttZ_tot_plus << "/-" << xsecErr_ttZ_tot_minus << " pb" << std::endl;
  std::cout << " Expected ttV Cross Section: " << crossSection_ttV << " pb" << std::endl;
  std::cout << " Measured ttV Cross Section: " << std::endl;
  std::cout << " " <<  crossSectionObs_ttV << "  +" << xsecErr_ttV_stat_plus << "/-" << xsecErr_ttV_stat_minus << " (stat)   +" << xsecErr_ttV_syst_plus << "/-" << xsecErr_ttV_syst_minus << " (syst)  pb" << std::endl;
  std::cout << " " <<  crossSectionObs_ttV << "  +" << xsecErr_ttV_tot_plus << "/-" << xsecErr_ttV_tot_minus << " pb" << std::endl;
  std::cout << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;


  return 0;

}



std::pair< float, float > getSystFromString( const std::string& systString ) {

  TString systString_tstr(systString);
  float systDOWN, systUP;
  if( systString_tstr.Contains("/") ) {
    sscanf( systString.c_str(), "%f/%f", &systDOWN, &systUP );
  } else if( systString=="-" ) {
    systUP = 0.;
    systDOWN = 0.;
  } else {
    systUP = atof(systString.c_str());
    systDOWN = systUP;
  }

  std::pair< float, float> returnPair;
  returnPair.first = systDOWN;
  returnPair.second = systUP;

  return returnPair;

}

