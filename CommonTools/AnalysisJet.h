#ifndef AnalysisJet_h
#define AnalysisJet_h

#include "TLorentzVector.h"
#include <string>

using namespace std;

class AnalysisJet : public TLorentzVector {

 public:

// ---- Constructor & Destructor ----------------------------------------------------------------------------
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
  
  virtual ~AnalysisJet(){}
// ---- Member functions ------------------------------------------------------------------------
   Bool_t passedJetID( const std::string& strength ) const;

// ---- Data Members ----------------------------------------------------------------------------
  float rmsCand;
  float ptD;
  int nCharged; //DON'T USE
  int nNeutral; //DON'T USE
  float QGLikelihood;
  float QGLikelihoodNoPU;
  
  //extra QGL variables
  float axis1;
  float axis2;
  float pull;
  float tana;

  float ptD_QC;
  float rmsCand_QC;
  float axis1_QC;
  float axis2_QC;
  float pull_QC;
  float tana_QC;

  float nChg_ptCut;
  float nChg_QC;
  float nChg_ptCut_QC;
  float nNeutral_ptCut;
  float Rchg;
  float Rneutral;
  float R;
  float Rchg_QC;
  //

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
  int pdgIdPartMom;
  int pdgIdPartMomMom;

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
  float combinedSecondaryVertexBJetTag;
  float simpleSecondaryVertexHighEffBJetTag;
  float simpleSecondaryVertexHighPurBJetTag;
  float jetBProbabilityBJetTag;
  float jetProbabilityBJetTag;

  float beta;
  float betaStar;


private:
  int nConstituents() const {return nChargedHadrons+nPhotons+nNeutralHadrons+nElectrons+nMuons+nHFHadrons+nHFEM ;}
};















#endif
