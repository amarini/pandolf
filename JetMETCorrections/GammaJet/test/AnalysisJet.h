#ifndef AnalysisJet_h
#define AnalysisJet_h

#include "TROOT.h"
#include "TLorentzVector.h"


class AnalysisJet : public TLorentzVector {

 public:
  AnalysisJet() : TLorentzVector(0.,0.,0.,0.) {
    ptReco=0.;
    ptCorrReco=0.;
    eReco=0.;
    etaReco=0.;
    phiReco=0.;

    nTracksReco=0;
    nNeutralHadronsReco=0;
    nPhotonsReco=0;
    nElectronsReco=0;
    nMuonsReco=0;
    nHFEMReco=0;
    nHFHadronsReco=0;

    eTracksReco=0.;
    eNeutralHadronsReco=0.;
    ePhotonsReco=0.;
    eElectronsReco=0.;
    eMuonsReco=0.;
    eHFEMReco=0.;
    eHFHadronsReco=0.;

  };
  ~AnalysisJet(){};

  Float_t eReco;
  Float_t ptReco;
  Float_t etaReco;
  Float_t phiReco;
  Float_t ptCorrReco;

  Float_t ptD;
  Float_t rmsCand;

  Float_t beta;
  Float_t betaStar;

  Float_t eGen;
  Float_t ptGen;
  Float_t etaGen;
  Float_t phiGen;

  Float_t ePart;
  Float_t ptPart;
  Float_t etaPart;
  Float_t phiPart;
  Int_t pdgIdPart;
  Int_t pdgIdMom;

  Float_t ePartStatus3;
  Float_t ptPartStatus3;
  Float_t etaPartStatus3;
  Float_t phiPartStatus3;
  Int_t pdgIdPartStatus3;
  Int_t pdgIdMomStatus3;

  Float_t thetaReco() const;
  Float_t pReco() const;
  Float_t pxReco() const;
  Float_t pyReco() const;
  Float_t pzReco() const;
  Float_t eCorrReco() const;
  
  Float_t thetaGen() const;
  Float_t pGen() const;
  Float_t pxGen() const;
  Float_t pyGen() const;
  Float_t pzGen() const;

  Float_t emfReco;

  Float_t eTracksReco;
  Float_t ePhotonsReco;
  Float_t eNeutralHadronsReco;
  Float_t eMuonsReco;
  Float_t eElectronsReco;
  Float_t eHFHadronsReco;
  Float_t eHFEMReco;

  Float_t ptTracksReco;
  Float_t ptPhotonsReco;
  Float_t ptNeutralHadronsReco;
  Float_t ptMuonsReco;
  Float_t ptElectronsReco;
  Float_t ptHFHadronsReco;
  Float_t ptHFEMReco;


  Int_t nTracksReco;
  Int_t nPhotonsReco;
  Int_t nNeutralHadronsReco;
  Int_t nMuonsReco;
  Int_t nElectronsReco;
  Int_t nHFHadronsReco;
  Int_t nHFEMReco;

  Int_t nCharged() const;
  Int_t nNeutral() const;
  Int_t nConstituents() const;

  Float_t eTracksGen;
  Float_t ePhotonsGen;
  Float_t eNeutralHadronsGen;
  Float_t eMuonsGen;
  Float_t eElectronsGen;
  Float_t eHFHadronsGen;
  Float_t eHFEMGen;

  Float_t ptTracksGen;
  Float_t ptPhotonsGen;
  Float_t ptNeutralHadronsGen;
  Float_t ptMuonsGen;
  Float_t ptElectronsGen;
  Float_t ptHFHadronsGen;
  Float_t ptHFEMGen;

  Int_t nTracksGen;
  Int_t nPhotonsGen;
  Int_t nNeutralHadronsGen;
  Int_t nMuonsGen;
  Int_t nElectronsGen;
  Int_t nHFHadronsGen;
  Int_t nHFEMGen;

  Float_t QGLikelihood;

  //btag:
  Float_t trackCountingHighEffBJetTags;
  Float_t simpleSecondaryVertexHighEffBJetTags;

  Bool_t passedJetID( const std::string& strength="minimal" ) const;

};



#endif
