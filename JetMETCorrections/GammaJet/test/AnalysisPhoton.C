#include "AnalysisPhoton.h"
#include <iostream>
#include <cstdlib>

  
bool AnalysisPhoton::isIsolated_hcal( const std::string& tightness ) const {

  Float_t hcalIso_rel_thresh;
  Float_t hcalIso_abs_thresh;

  if( tightness == "loose" ) {
    hcalIso_rel_thresh = 0.1;
    hcalIso_abs_thresh = 4.;
  } else if( tightness == "medium" ) {
    hcalIso_rel_thresh = 0.05;
    hcalIso_abs_thresh = 2.4;
  } else if( tightness == "tight" ) {
    hcalIso_rel_thresh = 0.03;
    hcalIso_abs_thresh = 2.4;
  } else {
    std::cout << "WARNING: no isolation \"" << tightness << "\" exists. Exiting." << std::endl;
    exit(1);
  }

  bool pass = true;

  if( (this->hcalIso >= hcalIso_rel_thresh) && (this->hcalIso*this->e >= hcalIso_abs_thresh) ) pass=false;

  return pass;

}


bool AnalysisPhoton::isIsolated_ecal( const std::string& tightness ) const {

  Float_t ecalIso_rel_thresh;
  Float_t ecalIso_abs_thresh;

  if( tightness == "loose" ) {
    ecalIso_rel_thresh = 0.1;
    ecalIso_abs_thresh = 4.5;
  } else if( tightness == "medium" ) {
    ecalIso_rel_thresh = 0.05;
    ecalIso_abs_thresh = 3.;
  } else if( tightness == "tight" ) {
    ecalIso_rel_thresh = 0.03;
    ecalIso_abs_thresh = 3.;
  } else {
    std::cout << "WARNING: no isolation \"" << tightness << "\" exists. Exiting." << std::endl;
    exit(1);
  }

  bool pass = true;

  if( (this->ecalIso >= ecalIso_rel_thresh) && (this->ecalIso*this->e >= ecalIso_abs_thresh) ) pass=false;

  return pass;

}


bool AnalysisPhoton::isIsolated_ptTracks( const std::string& tightness ) const {

  Float_t ptTrkIso_thresh;

  if( tightness == "loose" ) {
    ptTrkIso_thresh = 0.2;
  } else if( tightness == "medium" ) {
    ptTrkIso_thresh = 0.1;
  } else if( tightness == "tight" ) {
    ptTrkIso_thresh = 0.05;
  } else {
    std::cout << "WARNING: no isolation \"" << tightness << "\" exists. Exiting." << std::endl;
    exit(1);
  }

  bool pass = true;

  if( this->ptTrkIso >= ptTrkIso_thresh ) pass=false;

  return pass;
}


bool AnalysisPhoton::isIsolated_nTracks( const std::string& tightness ) const {

  Int_t nTrkIso_thresh;

  if( tightness == "loose" ) {
    nTrkIso_thresh = 5;
  } else if( tightness == "medium" ) {
    nTrkIso_thresh = 3;
  } else if( tightness == "tight" ) {
    nTrkIso_thresh = 1;
  } else {
    std::cout << "WARNING: no isolation \"" << tightness << "\" exists. Exiting." << std::endl;
    exit(1);
  }

  bool pass=true;

  if( this->nTrkIso >= nTrkIso_thresh ) pass=false;

  return pass;

}


bool AnalysisPhoton::isIsolated( const std::string& tightness ) const {
  
  if( !this->isIsolated_hcal( tightness ) ) return false;
  if( !this->isIsolated_ecal( tightness ) ) return false;
  if( !this->isIsolated_nTracks( tightness ) ) return false;
  if( !this->isIsolated_ptTracks( tightness ) ) return false;

  return true;

}

bool AnalysisPhoton::clusterMajOK( const std::string& tightness ) const {

  Float_t clusterMaj_thresh;

  if( tightness == "loose" ) {
    clusterMaj_thresh = 0.6;
  } else if( tightness == "medium" ) {
    clusterMaj_thresh = 0.35;
  } else if( tightness == "tight" ) {
    clusterMaj_thresh = 0.35;
  } else {
    std::cout << "WARNING: no PhotonID \"" << tightness << "\" exists. Exiting." << std::endl;
    exit(1);
  }

  bool pass=true;

  if( this->clusterMaj >= clusterMaj_thresh || this->clusterMaj <= 0.15 ) pass=false;

  return pass;

}

bool AnalysisPhoton::clusterMinOK( const std::string& tightness ) const {

  Float_t clusterMin_thresh;

  if( tightness == "loose" ) {
    clusterMin_thresh = 0.5;
  } else if( tightness == "medium" ) {
    clusterMin_thresh = 0.3;
  } else if( tightness == "tight" ) {
    clusterMin_thresh = 0.3;
  } else {
    std::cout << "WARNING: no PhotonID \"" << tightness << "\" exists. Exiting." << std::endl;
    exit(1);
  }

  bool pass=true;

  if( this->clusterMin >= clusterMin_thresh || this->clusterMin <= 0.15 ) pass=false;

  return pass;

}

bool AnalysisPhoton::passedPhotonID( const std::string& tightness ) const {

  if( !this->isIsolated(tightness) ) return false;
  if( !this->clusterMajOK( tightness ) ) return false;
  if( !this->clusterMinOK( tightness ) ) return false;

  return true;

}
