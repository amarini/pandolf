#include "AnalysisJet.h"
#include <cmath>

AnalysisJet::AnalysisJet() { 

  eReco = 0.;
  ptReco = 0.;
  ptCorrReco=0.;
  etaReco = 20.; 
  phiReco = 0.;
  eGen = 0.;
  ptGen = 0.;
  etaGen = 20.;
  phiGen = 0.;

}

Float_t AnalysisJet::thetaReco() const { return 2.*atan(exp(-etaReco)); };
Float_t AnalysisJet::pReco() const { return ptReco/sin(thetaReco()); };
Float_t AnalysisJet::pxReco() const { return ptReco*cos(phiReco); };
Float_t AnalysisJet::pyReco() const { return ptReco*sin(phiReco); };
Float_t AnalysisJet::pzReco() const { return pReco()*cos(thetaReco()); };
Float_t AnalysisJet::eCorrReco() const { 
  Float_t eCorrReco = ( ptReco>0. ) ? ptCorrReco*eReco/ptReco : 0.;
  return eCorrReco;
}

Float_t AnalysisJet::thetaGen() const { return 2.*atan(exp(-etaGen)); };
Float_t AnalysisJet::pGen() const { return ptGen/sin(thetaGen()); };
Float_t AnalysisJet::pxGen() const { return ptGen*cos(phiGen); };
Float_t AnalysisJet::pyGen() const { return ptGen*sin(phiGen); };
Float_t AnalysisJet::pzGen() const { return pGen()*cos(thetaGen()); };

