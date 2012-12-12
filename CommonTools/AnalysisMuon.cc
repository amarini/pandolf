#include "AnalysisMuon.h"



bool AnalysisMuon::passedMuonID() {

  bool passed = ( isGlobalMuonPromptTight && isAllTrackerMuon 
                && pixelHits>0 && trackerHits>10 
                && nMatchedStations>=2
                && (dxy<0.02) && (dz<1.) );

  return passed;

}


bool AnalysisMuon::passedMuonID_tight2012() {

  bool passed = ( isGlobalMuon && isPFMuon 
                && normChiSquare < 10.
                && pixelHits>0 && trackerHits>5 
                && nValidMuonHits > 0
                && nMatchedStations>=2
                && (dxy<0.02) && (dz<0.5) ); 

  return passed;

}



bool AnalysisMuon::isIsolated() {

  bool isIsolated = this->combinedIsoRel() < 0.15;

  return isIsolated;

}


bool AnalysisMuon::isIsolated2012_CutsDeltaBeta() {

  bool isIsolated = (this->pfIsolationSumPUPtR04/this->Pt()) < 0.15;

  return isIsolated;

}





bool AnalysisMuon::isIsolated2012_MVA() {

  float abseta = fabs(this->Eta());
  float pt = this->Pt();
  float bdtiso = mvaisoMuon;

  bool isIsolated = (
                             ( pt <= 20 && abseta >= 0.000 && abseta < 0.479 && bdtiso > 0.86 ) ||
                             ( pt <= 20 && abseta >= 0.479 && abseta < 2.400 && bdtiso > 0.82 ) ||
                             ( pt >  20 && abseta >= 0.000 && abseta < 0.479 && bdtiso > 0.82 ) ||
                             ( pt >  20 && abseta >= 0.479 && abseta < 2.400 && bdtiso > 0.86 )
                             );


  return isIsolated;

}



float AnalysisMuon::combinedIsoRel() {

 return (sumPt03 + emEt03 + hadEt03)/this->Pt();

}


bool AnalysisMuon::passedVBTF() {

  bool isIsolated = this->isIsolated();
  bool passedMuonID = this->passedMuonID();

  return ( isIsolated && passedMuonID );

}


bool AnalysisMuon::isTightMuon2012() {

  bool isIsolated = this->isIsolated2012_CutsDeltaBeta();
  bool passedMuonID = this->passedMuonID_tight2012();

  return ( isIsolated && passedMuonID );

}

