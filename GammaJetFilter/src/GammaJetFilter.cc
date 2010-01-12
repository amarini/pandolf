// -*- C++ -*-
//
// Package:    GammaJetFilter
// Class:      GammaJetFilter
// 
/**\class GammaJetFilter GammaJetFilter.cc JetMETCorrections/GammaJetFilter/src/GammaJetFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Francesco PANDOLFI
//         Created:  Tue Dec  1 19:53:14 CET 2009
// $Id$
//
//


// system include files
#include <memory>
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"


//
// class declaration
//

class GammaJetFilter : public edm::EDFilter {
   public:
      explicit GammaJetFilter(const edm::ParameterSet&);
      ~GammaJetFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      edm::InputTag Photonsrc_;
      edm::InputTag JetPFsrcakt5_;
      edm::InputTag JetPFsrcakt7_;
      double PtThreshold_;
      unsigned nFoundEvents_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GammaJetFilter::GammaJetFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   Photonsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("Photonsrc");
   JetPFsrcakt5_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfakt5");
   JetPFsrcakt7_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfakt7");

   nFoundEvents_=0;

}


GammaJetFilter::~GammaJetFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GammaJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   bool okForAntiKt5=false;
   bool okForAntiKt7=false;
   bool eventOK = false;

   Handle<PhotonCollection>  PhotonHandle;
   iEvent.getByLabel(Photonsrc_, PhotonHandle);

   Handle<PFJetCollection>  PFJetAntikt5Handle;
   iEvent.getByLabel(JetPFsrcakt5_, PFJetAntikt5Handle);

   Handle<PFJetCollection>  PFJetAntikt7Handle;
   iEvent.getByLabel(JetPFsrcakt7_, PFJetAntikt7Handle);

   bool foundPhot=false;
   Float_t ptPhot=0.;
   Float_t phiPhot=0.;
   Float_t etaPhot=0.;

   for (PhotonCollection::const_iterator it = PhotonHandle->begin(); (it!=PhotonHandle->end())&&(!foundPhot); ++it) {

     if( (it->pt()>5.)&&(fabs(it->eta())<=2.4) ) {

       ptPhot = it->pt();
       etaPhot = it->eta();
       phiPhot = it->phi();
       foundPhot=true;

     }

   }

   if( foundPhot ) {

     std::cout << "Found photon. pt: " << ptPhot << "\teta: " << etaPhot << "\tphi: " << phiPhot << std::endl;
     
     //look for matched jets:
     //anti-kt0.5 pfjets:

     Float_t deltaRmin5 = 999.;
     const PFJet* foundJet5=0;
 
     for (PFJetCollection::const_iterator iJet5 = PFJetAntikt5Handle->begin(); iJet5!=PFJetAntikt5Handle->end(); ++iJet5) {

       Float_t deltaEta = etaPhot - iJet5->eta();
       Float_t deltaPhi = phiPhot - iJet5->phi();
       if( deltaPhi > TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
       if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

       Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

       if( deltaR < deltaRmin5 ) {
         deltaRmin5 = deltaR;
         foundJet5 = &(*iJet5);
       }

     } //for antikt5 pfjets

     if( deltaRmin5<0.5 ) {

       std::cout << "Found anti-kt 0.5 jet. pt: " << foundJet5->pt() << "\teta: " << foundJet5->eta() << "\tphi: " << foundJet5->phi();
       Int_t nChargedHadrons=0;
       std::vector<const PFCandidate*> pfCandidates = foundJet5->getPFConstituents();

       for (std::vector<const PFCandidate*>::const_iterator jt = pfCandidates.begin(); jt != pfCandidates.end(); ++jt) {
         PFCandidate::ParticleType id = (*jt)->particleId();
         if (id==PFCandidate::h) nChargedHadrons++;
       }

       if( nChargedHadrons<3 ) okForAntiKt5=true;
       std::cout << "\tnchargedhadrons: " << nChargedHadrons << std::endl;

     }


     if( !okForAntiKt5 ) {

       //anti-kt0.7 pfjets:

       Float_t deltaRmin7 = 999.;
       const PFJet* foundJet7=0;
 
       for (PFJetCollection::const_iterator iJet7 = PFJetAntikt7Handle->begin(); iJet7!=PFJetAntikt7Handle->end(); ++iJet7) {
       
         Float_t deltaEta = etaPhot - iJet7->eta();
         Float_t deltaPhi = phiPhot - iJet7->phi();
         if( deltaPhi > TMath::Pi() ) deltaPhi -= 2.*TMath::Pi();
         if( deltaPhi < -TMath::Pi() ) deltaPhi += 2.*TMath::Pi();

         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

         if( deltaR < deltaRmin7 ) {
           deltaRmin7 = deltaR;
           foundJet7 = &(*iJet7);
         }

       } //for antikt7 pfjets

       if( deltaRmin7<0.7 ) {

         std::cout << "Found anti-kt 0.7 jet. pt: " << foundJet7->pt() << "\teta: " << foundJet7->eta() << "\tphi: " << foundJet7->phi();
         Int_t nChargedHadrons=0;
         std::vector<const PFCandidate*> pfCandidates = foundJet7->getPFConstituents();

         for (std::vector<const PFCandidate*>::const_iterator jt = pfCandidates.begin(); jt != pfCandidates.end(); ++jt) {
           PFCandidate::ParticleType id = (*jt)->particleId();
           if (id==PFCandidate::h) nChargedHadrons++;
         }

         if( nChargedHadrons<4 ) okForAntiKt7=true;
         std::cout << "\tnchargedhadrons: " << nChargedHadrons << std::endl;

       } //if deltaRmin

     } //if !okForAntiKt5


     if( okForAntiKt5 || okForAntiKt7 ) {
       eventOK=true;
       nFoundEvents_++;
     }

   } //if foundPhot

   return eventOK;

}

// ------------ method called once each job just before starting event loop  ------------
void 
GammaJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GammaJetFilter::endJob() {

  std::cout << "-----------------------> FOUND " << nFoundEvents_ << " GOOD EVENTS." << std::endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaJetFilter);
