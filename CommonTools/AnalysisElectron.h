// ---------------------------------------------------------------
//
//  AnalysisElectron - electron class used in the HZZlljj analysis
//
// ---------------------------------------------------------------


#ifndef AnalysisElectron_h
#define AnalysisElectron_h

#include "AnalysisLepton.h"




class AnalysisElectron : public AnalysisLepton {

 public:

  AnalysisElectron( float x=0., float y=0., float z=0., float t=0.) : AnalysisLepton( x, y, z, t ) {};
  AnalysisElectron( double x=0., double y=0., double z=0., double t=0.) : AnalysisLepton( x, y, z, t ) {};

  AnalysisElectron( const TLorentzVector &v) : AnalysisLepton( v ) {};

  bool isIsolatedVBTF80();
  bool isIsolatedVBTF95();
  bool isIsolated2012_MVAWP95();
  bool isIsolated2012_CutsLoose();
  bool separatedIsoRel();//

  bool electronIDVBTF80();
  bool electronIDVBTF95();
  bool electronID2012_MVA();
  bool electronID2012_CutsLoose();
  bool separatedIDVBTF80();//
  bool passedHLT2012();

  bool conversionRejectionVBTF80();
  bool conversionRejectionVBTF95();
  bool conversionRejection2012_CutsLoose();

  bool passedVBTF80();
  bool passedVBTF95();
  bool passedTrigger80();//

  bool isGoodElectron2012_CutsLoose();
  bool isGoodElectron2012_MVA();

  bool passedAdditionalCuts();

  double combinedIsoRel();


  // public data members:

  double etaSC; //supercluster eta


  // impact parameter:
  double dxy;
  double dz;

  // isolation:
  double dr03TkSumPt;
  double dr03EcalRecHitSumEt;
  double dr03HcalTowerSumEt;

  double pfCandChargedIso04;
  double pfCandNeutralIso04;
  double pfCandPhotonIso04;
  double rhoJetsFastJet;

  // electron ID:
  double sigmaIetaIeta; 
  double deltaPhiAtVtx; 
  double deltaEtaAtVtx; 
  double hOverE; 
  double fBrem; 
  double eOverP; 
  double mvaidtrigEle;

  // conversion rejection:
  int expInnerLayersGsfTrack;
  double convDist;
  double convDcot;
  bool hasMatchedConversion;

};

#endif
