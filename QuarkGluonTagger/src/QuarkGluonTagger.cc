// system include files
#include <memory>
#include <iostream>
#include <vector>

#include "pandolf/QuarkGluonTagger/interface/QuarkGluonTagger.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"





QuarkGluonTagger::QuarkGluonTagger(const edm::ParameterSet& iConfig)
{
        src_        = iConfig.getParameter<edm::InputTag> ("jets");
        srcRho_     = iConfig.getParameter<edm::InputTag> ("rho");
        jecService_ = iConfig.getParameter<std::string>   ("jec");
        fileName_   = iConfig.getParameter<std::string>   ("filename");
        
        produces<edm::ValueMap<float> >().setBranchAlias("qg");
        //edm::FileInPath fip(fileName);
        qglikeli_ = new QGLikelihoodCalculator(fileName_);
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
    //----- calculate the jec -------------------
    int index = ijet-pfjets->begin();
    edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfjets,index));
    double cor = JEC_->correction(*ijet,iEvent,iSetup);
    //----- calculate the ptD ------------
    vector<PFCandidatePtr> pfConst(ijet->getPFConstituents());
    double sumpt(0.0),sumpt2(0.0); 
    for(unsigned ipf=0;ipf<pfConst.size();ipf++) {
      sumpt  += pfConst[ipf]->pt();
      sumpt2 += pfConst[ipf]->pt() * pfConst[ipf]->pt();
    }
    double ptD = sqrt(sumpt2)/sumpt;
    //----- calculate the likelihood ------------
    int nCharged = ijet->chargedHadronMultiplicity();
    int nNeutral = ijet->neutralHadronMultiplicity()+ijet->photonMultiplicity();
    double corPt = cor*ijet->pt();
    float qgl(-1.0);
    if (nCharged + nNeutral > 0 ) {
      qgl = qglikeli_->computeQGLikelihoodPU(corPt,*rho,nCharged,nNeutral,ptD);
    }
    //cout<<corPt<<" "<<ijet->eta()<<" "<<nCharged<<" "<<nNeutral<<" "<<ptD<<" "<<qgl<<endl;
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

