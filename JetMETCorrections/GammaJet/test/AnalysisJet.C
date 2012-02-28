#include "AnalysisJet.h"
#include <cmath>
#include <iostream>


Float_t AnalysisJet::thetaReco() const { return 2.*atan(exp(-etaReco)); };
Float_t AnalysisJet::pReco() const { return ptReco/sin(thetaReco()); };
Float_t AnalysisJet::pxReco() const { return ptReco*cos(phiReco); };
Float_t AnalysisJet::pyReco() const { return ptReco*sin(phiReco); };
Float_t AnalysisJet::pzReco() const { return pReco()*cos(thetaReco()); };
Float_t AnalysisJet::eCorrReco() const { 
  Float_t eCorr = ( ptReco>0. ) ? ptCorrReco*eReco/ptReco : 0.;
  return eCorr;
}

Float_t AnalysisJet::thetaGen() const { return 2.*atan(exp(-etaGen)); };
Float_t AnalysisJet::pGen() const { return ptGen/sin(thetaGen()); };
Float_t AnalysisJet::pxGen() const { return ptGen*cos(phiGen); };
Float_t AnalysisJet::pyGen() const { return ptGen*sin(phiGen); };
Float_t AnalysisJet::pzGen() const { return pGen()*cos(thetaGen()); };


Int_t AnalysisJet::nConstituents() const { return nTracksReco+nNeutralHadronsReco+nPhotonsReco+nElectronsReco+nMuonsReco+nHFEMReco+nHFHadronsReco; };
Int_t AnalysisJet::nCharged() const { return nTracksReco; };
Int_t AnalysisJet::nNeutral() const { return nNeutralHadronsReco + nPhotonsReco; };


Bool_t AnalysisJet::passedJetID( const std::string& strength ) const {

  Bool_t returnBool = false;

  if( strength=="minimal" ) {

    returnBool = (this->nConstituents()>1) &&
                 ((ePhotonsReco + eHFEMReco) <0.99) &&
                 ((eNeutralHadronsReco + eHFHadronsReco)<0.99);


    if( fabs(etaReco)<2.4 ) {
      returnBool = returnBool &&  (eTracksReco>0.);
    } 

  } else {

    std::cout << "Jet ID '" << strength << "' not implemented yet. Returning FALSE." << std::endl;

  }

  return returnBool;

}
