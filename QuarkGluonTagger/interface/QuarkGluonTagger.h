#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "../interface/QGLikelihoodCalculator.h"




class QuarkGluonTagger : public edm::EDProducer {
   public:
      explicit QuarkGluonTagger(const edm::ParameterSet&);
      ~QuarkGluonTagger();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data --------------------------
      edm::InputTag src_,srcRho_;
      std::string jecService_;
      QGLikelihoodCalculator *qglikeli_;
      const JetCorrector *JEC_;           
};

