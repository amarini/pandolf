#include <memory>
#include <iostream>
#include <vector>

#include "../interface/QuarkGluonTagger.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"





QuarkGluonTagger::QuarkGluonTagger(const edm::ParameterSet& iConfig)
{
        src_        = iConfig.getParameter<edm::InputTag> ("jets");
        srcRho_     = iConfig.getParameter<edm::InputTag> ("rho");
        jecService_ = iConfig.getParameter<std::string>   ("jec");
        
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
  iEvent.getByLabel(src_,pfjets);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  JEC_ = JetCorrector::getJetCorrector(jecService_,iSetup);

  std::vector<float> values;
  values.reserve(pfjets->size());
    
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
  }

  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(pfjets, values.begin(), values.end());
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

