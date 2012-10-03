// ---------------------------------------------------------------
//
//  AnalysisMuon - muon class used in the HZZlljj analysis
//
// ---------------------------------------------------------------


#ifndef AnalysisMuon_h
#define AnalysisMuon_h

#include "AnalysisLepton.h"




class AnalysisMuon : public AnalysisLepton {

 public:

  AnalysisMuon( float x=0., float y=0., float z=0., float t=0.) : AnalysisLepton( x, y, z, t ) {};
  AnalysisMuon( double x=0., double y=0., double z=0., double t=0.) : AnalysisLepton( x, y, z, t ) {};

  AnalysisMuon( const TLorentzVector &v) : AnalysisLepton( v ) {};

  bool passedMuonID();
  bool passedMuonID_tight2012();

  bool isIsolated();
  bool isIsolated2012();
  virtual float combinedIsoRel();

  bool passedVBTF();
  bool isTightMuon2012();


  // public data members:

  bool isGlobalMuon;
  bool isPFMuon;
  bool isGlobalMuonPromptTight;
  bool isAllTrackerMuon;

  float normChiSquare;

  int pixelHits;
  int trackerHits;

  int nMatchedStations;
  int nValidMuonHits;

  double dxy;
  double dz;

  // isolation:
  double sumPt03;
  double emEt03;
  double hadEt03;
  double mvaisoMuon;

};

#endif
