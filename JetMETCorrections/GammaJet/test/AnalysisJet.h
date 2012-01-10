#ifndef AnalysisJet_h
#define AnalysisJet_h

#include "TROOT.h"
#include "TLorentzVector.h"


class AnalysisJet : public TLorentzVector {

 public:
  AnalysisJet() : TLorentzVector() {};
  ~AnalysisJet(){};

  Float_t eReco;
  Float_t ptReco;
  Float_t etaReco;
  Float_t phiReco;
//Float_t eCorrReco;
  Float_t ptCorrReco;

  Float_t ptD;
  Float_t rmsCand;

  Float_t eGen;
  Float_t ptGen;
  Float_t etaGen;
  Float_t phiGen;

  Float_t ePart;
  Float_t ptPart;
  Float_t etaPart;
  Float_t phiPart;
  Int_t pdgIdPart;

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

  //btag:
  Float_t trackCountingHighEffBJetTags;
  Float_t simpleSecondaryVertexHighEffBJetTags;

};



#endif
