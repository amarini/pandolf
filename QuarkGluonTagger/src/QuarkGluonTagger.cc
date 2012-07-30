#include <memory>
#include <iostream>
#include <vector>

#include "../interface/QuarkGluonTagger.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"



QuarkGluonTagger::QuarkGluonTagger(const edm::ParameterSet& iConfig)
{
        src_        = iConfig.getParameter<edm::InputTag> ("jets");
        srcRho_     = iConfig.getParameter<edm::InputTag> ("rho");
        jecService_ = iConfig.getParameter<std::string>   ("jec");
	//        isPatJet_ = iConfig.getParameter<bool>   ("isPatJet");
	isPatJet_ = iConfig.existsAs<bool>("isPatJet") ? iConfig.getParameter<bool>("isPatJet") : false ; 
       
        produces<edm::ValueMap<float> >().setBranchAlias("qg");
        qglikeli_ = new QGLikelihoodCalculator();
}

QuarkGluonTagger::~QuarkGluonTagger()
{

}

// ------------ method called to produce the data  ------------
void
QuarkGluonTagger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;  

  edm::Handle<reco::PFJetCollection> pfjets;
  edm::Handle< vector<pat::Jet> > patjets;

  if(!isPatJet_)
    iEvent.getByLabel(src_,pfjets);
  else
    iEvent.getByLabel(src_,patjets);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  JEC_ = JetCorrector::getJetCorrector(jecService_,iSetup);

  std::vector<float> values;
  if(!isPatJet_)
    values.reserve(pfjets->size());
  else
    values.reserve(patjets->size());
  

  if(!isPatJet_) {
    for(reco::PFJetCollection::const_iterator ijet = pfjets->begin(); ijet != pfjets->end(); ++ijet) {
      
      //jet energy correction:
      double cor = JEC_->correction(*ijet,iEvent,iSetup);
      double corPt = cor*ijet->pt();
    
      // get the variables for the LD:
      double ptD = ijet->constituentPtDistribution();
      int nCharged = ijet->chargedHadronMultiplicity();
      int nNeutral = ijet->neutralHadronMultiplicity()+ijet->photonMultiplicity();
    
    
      // compute the LD:
      float qgl(-1.0);
      if (nCharged + nNeutral > 0 ) {
	if( fabs(ijet->eta())<2.4 )
	  qgl = qglikeli_->computeQGLikelihood(corPt,*rho,nCharged,nNeutral,ptD);
	else
	  qgl = -1.;
      }
      //cout<<corPt<<" "<<ijet->eta()<<" "<<nCharged<<" "<<nNeutral<<" "<<ptD<<" "<<qgl<<endl;

      // fill the value map:
      values.push_back(qgl);
    } // end loop PFjets
  } else {
    for(vector<pat::Jet>::const_iterator ijet = patjets->begin(); ijet != patjets->end(); ++ijet) {
      //NO jet energy correction applied (already done at PAT level)
      double corPt = ijet->pt();
    
      // get the variables for the LD:
      double ptD = ijet->constituentPtDistribution();
      int nCharged = ijet->chargedHadronMultiplicity();
      int nNeutral = ijet->neutralHadronMultiplicity()+ijet->photonMultiplicity();
    
    
      // compute the LD:
      float qgl(-1.0);
      if (nCharged + nNeutral > 0 ) {
	if( fabs(ijet->eta())<2.4 )
	  qgl = qglikeli_->computeQGLikelihood(corPt,*rho,nCharged,nNeutral,ptD);
	else
	  qgl = -1.;
      }
      //cout<<corPt<<" "<<ijet->eta()<<" "<<nCharged<<" "<<nNeutral<<" "<<ptD<<" "<<qgl<<endl;

      // fill the value map:
      values.push_back(qgl);
    } // end loop PATjets

  } // end if

  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  if(!isPatJet_)
    filler.insert(pfjets, values.begin(), values.end());
  else 
    filler.insert(patjets, values.begin(), values.end());

  filler.fill();
  
  // put value map into event
  iEvent.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void
QuarkGluonTagger::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
QuarkGluonTagger::endJob() {
}



//define this as a plug-in
DEFINE_FWK_MODULE(QuarkGluonTagger);

