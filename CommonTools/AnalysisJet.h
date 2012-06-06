#ifndef AnalysisJet_h
#define AnalysisJet_h

#include "TLorentzVector.h"

class AnalysisJet : public TLorentzVector {

 public:

  AnalysisJet( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {
    rmsCand=0.;
    ptD=0.;
    nCharged=0;
    nNeutral=0;
  }

  AnalysisJet( const TLorentzVector &v) : TLorentzVector( v ) {
    rmsCand=0.;
    ptD=0.;
    nCharged=0;
    nNeutral=0;
  }

  float rmsCand;
  float ptD;
  int nCharged;
  int nNeutral;
  float QGLikelihood;
  float QGLikelihoodNoPU;

  float eChargedHadrons;
  float ePhotons;
  float eNeutralHadrons;
  float eElectrons;
  float eMuons;
  float eHFHadrons;
  float eHFEM;

  int nChargedHadrons;
  int nPhotons;
  int nNeutralHadrons;
  int nElectrons;
  int nMuons;
  int nHFHadrons;
  int nHFEM;

  float pt_preKinFit;
  float eta_preKinFit;
  float phi_preKinFit;
  float e_preKinFit;

  float ptGen;
  float etaGen;
  float phiGen;
  float eGen;

  float ptPart;
  float etaPart;
  float phiPart;
  float ePart;
  int pdgIdPart;

  float ptPart_status2;
  float etaPart_status2;
  float phiPart_status2;
  float ePart_status2;
  int pdgIdPart_status2;

  bool btag_loose() const;
  bool btag_medium() const;

  //btags:
  float trackCountingHighEffBJetTag;
  float trackCountingHighPurBJetTag;
  float simpleSecondaryVertexHighEffBJetTag;
  float simpleSecondaryVertexHighPurBJetTag;
  float jetBProbabilityBJetTag;
  float jetProbabilityBJetTag;

};















#endif
