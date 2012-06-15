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



  //std::string mcZJetsFile = "ttZTrilepton_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11_" + suffix;
  //std::string mcTTbarFile = "ttZTrilepton_TTJ_Fall11_highstat_" + suffix;
  //std::string mcVVFile = "ttZTrilepton_VV_Summer11_" + suffix;
  //std::string mcttWFile = "ttZTrilepton_ttW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi_" + suffix;
  //bgTree->Add(dyTreeName.c_str());
  //bgTree->Add(ttTreeName.c_str());
  //bgTree->Add(dibosonTreeName.c_str());


  std::string datacardName = "datacard_TTWZ.txt";
  ifstream datacard(datacardName.c_str());


  int obs = 15;
  float ttZ = 2.245;
  float ttW = 5.164;
  float s = ttW+ttZ;
  float b_pred_err = 2.451;
   
  float b_fake=0.;
  float b_cmid=0.;
  float b_wz=0.;
  float b_rare=0.;

  // go get the syst on signal from datacard:
  std::cout << "-> Opened datacard file '" << datacardName << "'." << std::endl;
  bool go = false;
  //float lumiSystDOWN = 0.;
  //float lumiSystUP = 0.;
  float totalSignalSystDOWN = 0.;
  float totalSignalSystUP = 0.;
  float totalBG_fake_SystDOWN = 0.;
  float totalBG_fake_SystUP = 0.;
  float totalBG_cmid_SystDOWN = 0.;
  float totalBG_cmid_SystUP = 0.;
  float totalBG_wz_SystDOWN = 0.;
  float totalBG_wz_SystUP = 0.;
  float totalBG_rare_SystDOWN = 0.;
  float totalBG_rare_SystUP = 0.;
  while( datacard.good() ) {
    char line[500];
    datacard.getline( line, 500 );
    TString line_tstr(line);
    if( !go ) {
      char rate[100];
      if( (line_tstr.BeginsWith("rate")) ) 
        sscanf(line, "%s %f %f %f %f %f", rate, &s, &b_fake, &b_cmid, &b_wz, &b_rare);
      if( (line_tstr.BeginsWith("#syst")) ) go=true;
    } else {
      std::string systName, shape, systSignal, systBG_fake, systBG_cmid, systBG_wz, systBG_rare;
      datacard >> systName >> shape >> systSignal >> systBG_fake >> systBG_cmid >> systBG_wz >> systBG_rare;
      std::pair<float, float> pair_systSignal = getSystFromString(systSignal);

      float systSignalDOWN = pair_systSignal.first;
      float systSignalUP = pair_systSignal.second;

      std::pair<float, float> pair_systBG_fake = getSystFromString(systBG_fake);
      std::pair<float, float> pair_systBG_cmid = getSystFromString(systBG_cmid);
      std::pair<float, float> pair_systBG_wz = getSystFromString(systBG_wz);
      std::pair<float, float> pair_systBG_rare = getSystFromString(systBG_rare);

      float systBG_fake_DOWN = pair_systBG_fake.first;
      float systBG_fake_UP   = pair_systBG_fake.second;
      float systBG_cmid_DOWN = pair_systBG_cmid.first;
      float systBG_cmid_UP   = pair_systBG_cmid.second;
      float systBG_wz_DOWN   = pair_systBG_wz.first;
      float systBG_wz_UP     = pair_systBG_wz.second;
      float systBG_rare_DOWN = pair_systBG_rare.first;
      float systBG_rare_UP   = pair_systBG_rare.second;
   
      // last line reads zeroes dont know why
      if( systSignalDOWN==0. && systSignalUP==0. 
       && systBG_fake_DOWN==0. && systBG_fake_UP==0.
       && systBG_cmid_DOWN==0. && systBG_cmid_UP==0.
       && systBG_wz_DOWN==0. && systBG_wz_UP==0.
       && systBG_rare_DOWN==0. && systBG_rare_UP==0.
       ) continue;

      systSignalDOWN = (systSignalDOWN>0.) ? fabs( 1.-systSignalDOWN ) : 0.;
      systSignalUP   = (systSignalUP>0.) ? fabs( 1.-systSignalUP ) : 0.;

      totalSignalSystDOWN += systSignalDOWN*systSignalDOWN;
      totalSignalSystUP += systSignalUP*systSignalUP;


      systBG_fake_DOWN = (systBG_fake_DOWN>0.) ? fabs( 1.-systBG_fake_DOWN ) : 0.;
      systBG_fake_UP   = (systBG_fake_UP>0.) ? fabs( 1.-systBG_fake_UP ) : 0.;

      systBG_cmid_DOWN = (systBG_cmid_DOWN>0.) ? fabs( 1.-systBG_cmid_DOWN ) : 0.;
      systBG_cmid_UP   = (systBG_cmid_UP>0.) ? fabs( 1.-systBG_cmid_UP ) : 0.;

      systBG_wz_DOWN = (systBG_wz_DOWN>0.) ? fabs( 1.-systBG_wz_DOWN ) : 0.;
      systBG_wz_UP   = (systBG_wz_UP>0.) ? fabs( 1.-systBG_wz_UP ) : 0.;

      systBG_rare_DOWN = (systBG_rare_DOWN>0.) ? fabs( 1.-systBG_rare_DOWN ) : 0.;
      systBG_rare_UP   = (systBG_rare_UP>0.) ? fabs( 1.-systBG_rare_UP ) : 0.;

      //// do it again, as it's larger than 2:
      //systBG_rare_DOWN = (systBG_rare_DOWN>1.) ? fabs( 1.-systBG_rare_DOWN ) : systBG_rare_DOWN;
      //systBG_rare_UP   = (systBG_rare_UP>1.) ? fabs( 1.-systBG_rare_UP ) : systBG_rare_UP;

      std::cout << systName << std::endl;
      std::cout << "systBG_fake_DOWN: "  << systBG_fake_DOWN << std::endl;
      std::cout << "systBG_fake_UP  : "  << systBG_fake_UP   << std::endl;
                                        
      std::cout << "systBG_cmid_DOWN: "  << systBG_cmid_DOWN << std::endl;
      std::cout << "systBG_cmid_UP  : "  << systBG_cmid_UP   << std::endl;
                                        
      std::cout << "systBG_wz_DOWN  : "  << systBG_wz_DOWN << std::endl;
      std::cout << "systBG_wz_UP    : "  << systBG_wz_UP   << std::endl;
                                        
      std::cout << "systBG_rare_DOWN: "  << systBG_rare_DOWN << std::endl;
      std::cout << "systBG_rare_UP  : "  << systBG_rare_UP   << std::endl;

      totalBG_fake_SystDOWN += systBG_fake_DOWN*systBG_fake_DOWN;
      totalBG_fake_SystUP += systBG_fake_UP*systBG_fake_UP;

      totalBG_cmid_SystDOWN += systBG_cmid_DOWN*systBG_cmid_DOWN;
      totalBG_cmid_SystUP += systBG_cmid_UP*systBG_cmid_UP;

      totalBG_wz_SystDOWN += systBG_wz_DOWN*systBG_wz_DOWN;
      totalBG_wz_SystUP += systBG_wz_UP*systBG_wz_UP;

      totalBG_rare_SystDOWN += systBG_rare_DOWN*systBG_rare_DOWN;
      totalBG_rare_SystUP += systBG_rare_UP*systBG_rare_UP;

    }
  }

  if( !go ) {
    std::cout << "Didn't find systematics in datacard!" << std::endl;
    exit(135);
  }

   
  totalBG_fake_SystDOWN  = sqrt(totalBG_fake_SystDOWN);
  totalBG_cmid_SystDOWN  = sqrt(totalBG_cmid_SystDOWN);
  totalBG_wz_SystDOWN    = sqrt(totalBG_wz_SystDOWN  );
  totalBG_rare_SystDOWN  = sqrt(totalBG_rare_SystDOWN);
   
  totalBG_fake_SystUP  = sqrt(totalBG_fake_SystUP);
  totalBG_cmid_SystUP  = sqrt(totalBG_cmid_SystUP);
  totalBG_wz_SystUP    = sqrt(totalBG_wz_SystUP  );
  totalBG_rare_SystUP  = sqrt(totalBG_rare_SystUP);

  float totalBG_SystDOWN  = totalBG_fake_SystDOWN*totalBG_fake_SystDOWN*b_fake*b_fake +
                            totalBG_cmid_SystDOWN*totalBG_cmid_SystDOWN*b_cmid*b_cmid +
                            totalBG_wz_SystDOWN*totalBG_wz_SystDOWN*b_wz*b_wz +
                            totalBG_rare_SystDOWN*totalBG_rare_SystDOWN*b_rare*b_rare;
  float totalBG_SystUP  = totalBG_fake_SystUP*totalBG_fake_SystUP*b_fake*b_fake +
                          totalBG_cmid_SystUP*totalBG_cmid_SystUP*b_cmid*b_cmid +
                          totalBG_wz_SystUP*totalBG_wz_SystUP*b_wz*b_wz +
                          totalBG_rare_SystUP*totalBG_rare_SystUP*b_rare*b_rare;
   
  totalBG_SystDOWN = sqrt( totalBG_SystDOWN );
  totalBG_SystUP = sqrt( totalBG_SystUP );

  //totalBG_fake_SystDOWN  = b_fake*sqrt(totalBG_fake_SystDOWN);
  //totalBG_cmid_SystDOWN  = b_cmid*sqrt(totalBG_cmid_SystDOWN);
  //totalBG_wz_SystDOWN    = b_wz*sqrt(totalBG_wz_SystDOWN  );
  //totalBG_rare_SystDOWN  = b_rare*sqrt(totalBG_rare_SystDOWN);
  // 
  //totalBG_fake_SystUP  = b_fake*sqrt(totalBG_fake_SystUP);
  //totalBG_cmid_SystUP  = b_cmid*sqrt(totalBG_cmid_SystUP);
  //totalBG_wz_SystUP    = b_wz*sqrt(totalBG_wz_SystUP  );
  //totalBG_rare_SystUP  = b_rare*sqrt(totalBG_rare_SystUP);
   

  totalSignalSystDOWN = sqrt(totalSignalSystDOWN);
  totalSignalSystUP = sqrt(totalSignalSystUP);


  float b_pred = b_fake + b_cmid + b_wz + b_rare;


  // stat error on observed:
  double obs_plus, obs_minus;
  RooHistError::instance().getPoissonInterval(obs,obs_minus,obs_plus,1.);
  double obs_errPlus = obs_plus-obs;
  double obs_errMinus = obs-obs_minus;


  float ZBi = StatTools::computeZBi( s+b_pred, b_pred, b_pred_err );
  float ZBi_obs = StatTools::computeZBi( obs, b_pred, b_pred_err );

  float obs_ttV = obs - b_pred;
  float obs_ttW = obs_ttV - ttZ;

  //float nTotal_ttZ = 1467136.;
  //float nTotal_ttW = 1089608.;

  float lumi_pb = 4980.;
  float crossSection_ttZ = 0.139;
  float crossSection_ttW = 0.169;
  float crossSection_ttV = crossSection_ttZ+crossSection_ttW;

  //float n_passed_ttZ = nTotal_ttZ*ttZ/(crossSection_ttZ*lumi_pb);
  //float eff_ttZ = n_passed_ttZ/nTotal_ttZ;
  float eff_ttZ = ttZ/(crossSection_ttZ*lumi_pb);
std::cout << "eff_ttZ: " << 100.*eff_ttZ << "%" << std::endl;

  //float n_passed_ttW = nTotal_ttW*ttW/(crossSection_ttW*lumi_pb);
  //float eff_ttW = n_passed_ttW/nTotal_ttW;
  float eff_ttW = ttW/(crossSection_ttW*lumi_pb);
std::cout << "eff_ttW: " << 100.*eff_ttW << "%" << std::endl;

  float eff_ttV = ( crossSection_ttZ*eff_ttZ + crossSection_ttW*eff_ttW ) / ( crossSection_ttZ + crossSection_ttW );
std::cout << "eff_ttV: " << 100.*eff_ttV << "%" << std::endl;

  //float crossSectionObs_ttZ = obs_ttV / ( lumi_pb*eff_ttZ );
  float crossSectionObs_ttW = obs_ttW / ( lumi_pb*eff_ttW );
  float crossSectionObs_ttV = obs_ttV / ( lumi_pb*eff_ttV );

  //float xsecErr_ttZ_stat_plus  = obs_errPlus /  ( lumi_pb*eff_ttZ );
  //float xsecErr_ttZ_stat_minus = obs_errMinus / ( lumi_pb*eff_ttZ );

  float xsecErr_ttW_stat_plus  = obs_errPlus /  ( lumi_pb*eff_ttW );
  float xsecErr_ttW_stat_minus = obs_errMinus / ( lumi_pb*eff_ttW );

  float xsecErr_ttV_stat_plus  = obs_errPlus /  ( lumi_pb*eff_ttV );
  float xsecErr_ttV_stat_minus = obs_errMinus / ( lumi_pb*eff_ttV );



  // xsec = obs_ttV / ( lumi*eff )
  // signal syst are on lumi and eff
  // so err(eff) -> err(xsec) = obs_ttV/( lumi*eff*eff ) * err(eff) = xsec*err(eff)/eff
  // so err(lumi) -> err(xsec) = obs_ttV/( lumi*lumi*eff ) * err(lumi) = xsec*err(lumi)/lumi
  // total error is xsec*(quadrature sum of all relative errors)
  
  float xsecErr_ttW_systSignal_plus  = crossSection_ttW*(totalSignalSystUP);
  float xsecErr_ttW_systSignal_minus = crossSection_ttW*(totalSignalSystDOWN);

  float xsecErr_ttV_systSignal_plus  = crossSection_ttV*(totalSignalSystUP);
  float xsecErr_ttV_systSignal_minus = crossSection_ttV*(totalSignalSystDOWN);

  
  // background systematics instead affect the numerator:
  // xsec = ( obs - BG ) / ( lumi*eff )
  // so err(BG) -> err(xsec) = err(BG) / ( lumi*eff )

  //float xsecErr_ttV_systBG_plus  = (totalBG_fake_SystUP*b_fake + totalBG_cmid_SystUP*b_cmid + totalBG_wz_SystUP*b_wz + totalBG_rare_SystUP*b_rare) / ( lumi_pb*eff_ttV );
  //float xsecErr_ttV_systBG_minus  = (totalBG_fake_SystDOWN*b_fake + totalBG_cmid_SystDOWN*b_cmid + totalBG_wz_SystDOWN*b_wz + totalBG_rare_SystDOWN*b_rare) / ( lumi_pb*eff_ttV );

  float xsecErr_ttW_systBG_plus   = totalBG_SystUP   / ( lumi_pb*eff_ttW );
  float xsecErr_ttW_systBG_minus  = totalBG_SystDOWN / ( lumi_pb*eff_ttW );

  float xsecErr_ttV_systBG_plus   = totalBG_SystUP   / ( lumi_pb*eff_ttV );
  float xsecErr_ttV_systBG_minus  = totalBG_SystDOWN / ( lumi_pb*eff_ttV );


  float xsecErr_ttW_syst_plus  = sqrt( xsecErr_ttW_systSignal_plus *xsecErr_ttW_systSignal_plus  + xsecErr_ttW_systBG_plus *xsecErr_ttW_systBG_plus );
  float xsecErr_ttW_syst_minus = sqrt( xsecErr_ttW_systSignal_minus*xsecErr_ttW_systSignal_minus + xsecErr_ttW_systBG_minus*xsecErr_ttW_systBG_minus );

  float xsecErr_ttV_syst_plus  = sqrt( xsecErr_ttV_systSignal_plus *xsecErr_ttV_systSignal_plus  + xsecErr_ttV_systBG_plus *xsecErr_ttV_systBG_plus );
  float xsecErr_ttV_syst_minus = sqrt( xsecErr_ttV_systSignal_minus*xsecErr_ttV_systSignal_minus + xsecErr_ttV_systBG_minus*xsecErr_ttV_systBG_minus );


  float xsecErr_ttW_tot_plus  = sqrt( xsecErr_ttW_stat_plus *xsecErr_ttW_stat_plus  + xsecErr_ttW_syst_plus *xsecErr_ttW_syst_plus );
  float xsecErr_ttW_tot_minus = sqrt( xsecErr_ttW_stat_minus*xsecErr_ttW_stat_minus + xsecErr_ttW_syst_minus*xsecErr_ttW_syst_minus );

  float xsecErr_ttV_tot_plus  = sqrt( xsecErr_ttV_stat_plus *xsecErr_ttV_stat_plus  + xsecErr_ttV_syst_plus *xsecErr_ttV_syst_plus );
  float xsecErr_ttV_tot_minus = sqrt( xsecErr_ttV_stat_minus*xsecErr_ttV_stat_minus + xsecErr_ttV_syst_minus*xsecErr_ttV_syst_minus );


  std::cout << std::endl << std::endl << std::endl;
  std::cout << "=========================================" << std::endl;
  std::cout << std::endl;
  std::cout << " Expected BG: " << b_pred << " +- " << b_pred_err << std::endl;
  std::cout << " Expected Signal: " << s << std::endl;
  std::cout << " Expected B+S: " << s+b_pred << std::endl;
  std::cout << " Observed: " << obs << std::endl;
  std::cout << " Observed Poisson Interval: " << obs_minus << "-" << obs_plus << std::endl;
  std::cout << " Observed ttW Signal (BG subtracted): " << obs_ttW << std::endl;
  std::cout << " Expected ttW Cross Section: " << crossSection_ttW << " pb" << std::endl;
  std::cout << " Measured ttW Cross Section: " << std::endl;
  std::cout << " " <<  crossSectionObs_ttW << "  +" << xsecErr_ttW_stat_plus << "/-" << xsecErr_ttW_stat_minus << " (stat)   +" << xsecErr_ttW_syst_plus << "/-" << xsecErr_ttW_syst_minus << " (syst)  pb" << std::endl;
  std::cout << " " <<  crossSectionObs_ttW << "  +" << xsecErr_ttW_tot_plus << "/-" << xsecErr_ttW_tot_minus << " pb" << std::endl;
  std::cout << " Observed ttV Signal (BG subtracted): " << obs_ttV << std::endl;
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

