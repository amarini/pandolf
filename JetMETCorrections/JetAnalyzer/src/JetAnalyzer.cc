// -*- C++ -*-
//
// Package:    JetAnalyzer
// Class:      JetAnalyzer
// 
/**\class JetAnalyzer JetAnalyzer.cc MyAnalysis/JetAnalyzer/src/JetAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Daniele del Re
// Modified by Francesco Pandolfi
//         Created:  Thu Sep 13 16:00:15 CEST 2007
// $Id: JetAnalyzer.cc,v 1.2 2008/05/28 18:22:13 delre Exp $
//
//

//
// constructors and destructor
//

#include "JetMETCorrections/JetAnalyzer/interface/JetAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"


#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

#include "RecoEcal/EgammaCoreTools/interface/ClusterShapeAlgo.h"

#include "MyAnalysis/IsolationTools/interface/SuperClusterHitsEcalIsolation.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"


using namespace edm;
using namespace std;
using namespace reco;

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig)
{
  MCTruthCollection_ = iConfig.getUntrackedParameter<edm::InputTag>("MCTruthCollection");
  trackTags_ = iConfig.getUntrackedParameter<edm::InputTag>("tracks");
  
}


JetAnalyzer::~JetAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   nMC = 0;
   nPhot=0;
   nJetGen_itCone5=0;
   //nJetGen_sisCone5=0;
   //nJetGen_kt4=0;
   nJetCalo_itCone5=0;
   //nJetCalo_sisCone5=0;
   //nJetCalo_kt4=0;
   nJetPF_itCone5=0;
   //nCandPF_itCone5=0;
   //nJetPF_sisCone5=0;
   //nJetPF_kt4=0;


   // get generated pt hat
   Handle<double> genEventScale; 
   iEvent.getByLabel( "genEventScale", genEventScale ); 

   Handle<GenEventInfoProduct> genEventInfo; 
   iEvent.getByLabel( "generator", genEventInfo ); 

   // get MC info from GenParticleCandidates 
   Handle<GenParticleCollection> genParticles;
   iEvent.getByLabel( "genParticles", genParticles );

   Handle<reco::PhotonCollection> photons_h;
   iEvent.getByLabel( "photons", photons_h );
   const reco::PhotonCollection* photons = photons_h.product();
   
   Handle<reco::TrackCollection> tracks_h;
   iEvent.getByLabel( "generalTracks", tracks_h);
   const reco::TrackCollection* tracks = tracks_h.product();

   // get HCAL info
   Handle<HBHERecHitCollection> hbhe;
   iEvent.getByLabel( "hbhereco", hbhe);
   const HBHERecHitMetaCollection mhbhe(*hbhe);
 
   // get ECAL reco hits
   Handle<EBRecHitCollection> ecalhits;
   const EBRecHitCollection* rhits=0;
   iEvent.getByLabel( "ecalRecHit", "EcalRecHitsEB", ecalhits);
   const EcalRecHitMetaCollection mecalhits(*ecalhits);    
   rhits = ecalhits.product(); // get a ptr to the product

   Handle<EERecHitCollection> ecalhitsee;
   const EERecHitCollection* rhitsee=0;
   iEvent.getByLabel( "ecalRecHit", "EcalRecHitsEE", ecalhitsee);
   rhitsee = ecalhits.product(); // get a ptr to the product

   // get geometry
   edm::ESHandle<CaloGeometry> geoHandle;
   //   iSetup.get<IdealGeometryRecord>().get(geoHandle);
   iSetup.get<CaloGeometryRecord>().get(geoHandle);
   const CaloGeometry* geometry = geoHandle.product();
   const CaloSubdetectorGeometry* geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

   // get topology
   const CaloSubdetectorTopology *topology_p;
   edm::ESHandle<CaloTopology> topoHandle;
   iSetup.get<CaloTopologyRecord>().get(topoHandle);
   topology_p = topoHandle->getSubdetectorTopology(DetId::Ecal, EcalBarrel);

   //   edm::ESHandle<TrackerGeometry> trackerHandle_;
   edm::ESHandle<MagneticField> theMagField;
   iSetup.get<IdealMagneticFieldRecord>().get(theMagField);
   //   iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle_);
  

   // get gen jet collections
   Handle<GenJetCollection> genJets_itCone5_h;
   iEvent.getByLabel( "iterativeCone5GenJets", genJets_itCone5_h);
   const GenJetCollection* genJets_itCone5 = genJets_itCone5_h.product(); 
   
   Handle<GenJetCollection> genJets_sisCone5_h;
   iEvent.getByLabel( "sisCone5GenJets", genJets_sisCone5_h);
   const GenJetCollection* genJets_sisCone5 = genJets_sisCone5_h.product(); 
   
   Handle<GenJetCollection> genJets_kt4_h;
   iEvent.getByLabel( "kt4GenJets", genJets_kt4_h);
   const GenJetCollection* genJets_kt4 = genJets_kt4_h.product(); 

  
 

   // get reco jet collections
   Handle<CaloJetCollection> caloJets_itCone5_h;
   iEvent.getByLabel( "iterativeCone5CaloJets", caloJets_itCone5_h);
   const CaloJetCollection* caloJets_itCone5 = caloJets_itCone5_h.product(); 
   
   Handle<CaloJetCollection> caloJets_sisCone5_h;
   iEvent.getByLabel( "sisCone5CaloJets", caloJets_sisCone5_h);
   const CaloJetCollection* caloJets_sisCone5 = caloJets_sisCone5_h.product(); 
   
   Handle<CaloJetCollection> caloJets_kt4_h;
   iEvent.getByLabel( "kt4CaloJets", caloJets_kt4_h);
   const CaloJetCollection* caloJets_kt4 = caloJets_kt4_h.product(); 

   Handle<PFJetCollection> pfJets_itCone5_h;
   iEvent.getByLabel( "iterativeCone5PFJets", pfJets_itCone5_h);
   const PFJetCollection* pfJets_itCone5 = pfJets_itCone5_h.product();
/*
   Handle<PFJetCollection> pfJets_sisCone5_h;
   iEvent.getByLabel( "sisCone5PFJets", pfJets_sisCone5_h);
   const PFJetCollection* pfJets_sisCone5 = pfJets_sisCone5_h.product();

   Handle<PFJetCollection> pfJets_kt4_h;
   iEvent.getByLabel( "kt4PFJets", pfJets_kt4_h);
   const PFJetCollection* pfJets_kt4 = pfJets_kt4_h.product();
  */ 

  ClusterShapeAlgo algo;

 
   // Loop over MC truth
   
   ptHat = genEventInfo->qScale();   

   for ( GenParticleCollection::const_iterator p = genParticles->begin();
	p != genParticles->end(); ++p ) {
    
 
     if( p->energy()<1. ) continue;

     //if(nMC>=2000) {cout << "number of MC particle is larger than 2000. Skipping" << endl; continue;}
     if(nMC>=80) {continue;}  // to reduce the root file size
  
     pdgIdMC[nMC] = p->pdgId();
     statusMC[nMC] = p->status();
     massMC[nMC] = p->mass();
     pxMC[nMC] = p->px();	 
     pyMC[nMC] = p->py();	 
     pzMC[nMC] = p->pz();	 
     eMC[nMC] = p->energy();	 
     etaMC[nMC] = p->eta();	 
     phiMC[nMC] = p->phi();	 

     nMC++;

   }
   
   

   // Loop over reco photons

   for ( reco::PhotonCollection::const_iterator it = photons->begin(); 
	 it != photons->end(); it++ ) {
     
     if(it->energy()<5.) continue;
     if(nPhot>=40) {cout << "number of photons is larger than 40. Skipping" << endl; continue;}
     
     pxPhot[nPhot] = it->px();	 
     pyPhot[nPhot] = it->py();	 
     pzPhot[nPhot] = it->pz();	 
     ePhot[nPhot] = it->energy();	 
     escPhot[nPhot] = it->superCluster()->energy();	 
     etaPhot[nPhot] = it->eta();	 
     phiPhot[nPhot] = it->phi();	      
     
     double ptiso0015(0.),  ptiso035(0.), ptiso05(0.), ptiso07(0.), ptiso1(0.);
     int ntrkiso0015(0), ntrkiso035(0), ntrkiso05(0), ntrkiso07(0), ntrkiso1(0);
//      double ptisoatecal02(0.), ptisoatecal025(0.), ptisoatecal03(0.), ptisoatecal035(0.), ptisoatecal04(0.);
//      int ntrkisoatecal02(0), ntrkisoatecal025(0), ntrkisoatecal03(0), ntrkisoatecal035(0), ntrkisoatecal04(0);
     
     // calculate track isolation for different cone values

     for(TrackCollection::const_iterator itTrack = tracks->begin();
	 itTrack != tracks->end();                      
	 ++itTrack) {
       double etaTrack = itTrack->innerMomentum().eta();  
       double phiTrack = itTrack->innerMomentum().phi();
       
       double deltaPhi = phiTrack-it->phi();
       double deltaEta = etaTrack-it->eta();
       if(deltaPhi > Geom::pi()) deltaPhi -= 2.*Geom::pi();
       if(deltaPhi < -Geom::pi()) deltaPhi += 2.*Geom::pi();
       double deltaR = std::sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
       
       if(deltaR < .015)  {ptiso0015  += sqrt(itTrack->innerMomentum().Mag2()); ntrkiso0015++; }
       if(deltaR < .35)   {ptiso035 += sqrt(itTrack->innerMomentum().Mag2()); ntrkiso035++;}
       if(deltaR < .5)    {ptiso05  += sqrt(itTrack->innerMomentum().Mag2()); ntrkiso05++; }
       if(deltaR < .7)    {ptiso07 += sqrt(itTrack->innerMomentum().Mag2()); ntrkiso07++;}
       if(deltaR < 1.)    {ptiso1  += sqrt(itTrack->innerMomentum().Mag2()); ntrkiso1++; }
       

       // Trk isolation at ECAL commented for now

//        const math::XYZVector tempvec = itTrack->innerMomentum();
//        const math::XYZPoint temppos = itTrack->innerPosition();
//        int tempcharge = itTrack->charge();
       
//        ECALPositionCalculator tempcalculator;      
       
//        double phitrk = tempcalculator.ecalPhi(tempvec,temppos,tempcharge);
//        double etatrk = tempcalculator.ecalEta(tempvec,temppos);
       
//        deltaPhi = phitrk-it->phi();
//        deltaEta = etatrk-it->eta();      
       
//        deltaR = std::sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);
       
//        if(deltaPhi > Geom::pi()) deltaPhi -= 2.*Geom::pi();
//        if(deltaPhi < -Geom::pi()) deltaPhi += 2.*Geom::pi();
       
//        if(deltaR < .2)  {ptisoatecal02  += sqrt(itTrack->innerMomentum().Mag2()); ntrkisoatecal02++; }
//        if(deltaR < .25) {ptisoatecal025 += sqrt(itTrack->innerMomentum().Mag2()); ntrkisoatecal025++;}
//        if(deltaR < .3)  {ptisoatecal03  += sqrt(itTrack->innerMomentum().Mag2()); ntrkisoatecal03++; }
//        if(deltaR < .35) {ptisoatecal035 += sqrt(itTrack->innerMomentum().Mag2()); ntrkisoatecal035++;}
//        if(deltaR < .4)  {ptisoatecal04  += sqrt(itTrack->innerMomentum().Mag2()); ntrkisoatecal04++; }
       
     }
     
     ptiso0015Phot[nPhot] = ptiso0015;
     ntrkiso0015Phot[nPhot] = ntrkiso0015;
     ptiso035Phot[nPhot] = ptiso035;
     ntrkiso035Phot[nPhot] = ntrkiso035;
     ptiso05Phot[nPhot] = ptiso05;
     ntrkiso05Phot[nPhot] = ntrkiso05;
     ptiso07Phot[nPhot] = ptiso07;
     ntrkiso07Phot[nPhot] = ntrkiso07;
     ptiso1Phot[nPhot] = ptiso1;
     ntrkiso1Phot[nPhot] = ntrkiso1;

//      ptisoatecal02Phot[nPhot] = ptisoatecal02;
//      ntrkisoatecal02Phot[nPhot] = ntrkisoatecal02;
//      ptisoatecal025Phot[nPhot] = ptisoatecal025;
//      ntrkisoatecal025Phot[nPhot] = ntrkisoatecal025;
//      ptisoatecal03Phot[nPhot] = ptisoatecal03;
//      ntrkisoatecal03Phot[nPhot] = ntrkisoatecal03;
//      ptisoatecal035Phot[nPhot] = ptisoatecal035;
//      ntrkisoatecal035Phot[nPhot] = ntrkisoatecal035;
//      ptisoatecal04Phot[nPhot] = ptisoatecal04;
//      ntrkisoatecal04Phot[nPhot] = ntrkisoatecal04;

     // calculate HCAL isolation
  
     double hcalEnergy = 0;
     reco::SuperClusterRef sc = it->get<reco::SuperClusterRef>();
     CaloConeSelector selector1(0.1, geometry, DetId::Hcal);
     std::auto_ptr<CaloRecHitMetaCollectionV> selected=selector1.select(sc->eta(),sc->phi(),mhbhe);
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) hcalEnergy += hit->energy();
     hcalovecal01Phot[nPhot] = hcalEnergy/it->energy();
     hcalEnergy = 0;
     CaloConeSelector selector15(0.15, geometry, DetId::Hcal);
     selected = selector15.select(sc->eta(),sc->phi(),mhbhe);
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) hcalEnergy += hit->energy();
     hcalovecal015Phot[nPhot] = hcalEnergy/it->energy();
     hcalEnergy = 0;
     CaloConeSelector selector4(0.4, geometry, DetId::Hcal); 
     selected = selector4.select(sc->eta(),sc->phi(),mhbhe); 
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) hcalEnergy += hit->energy(); 
     hcalovecal04Phot[nPhot] = hcalEnergy/it->energy(); 
     hcalEnergy = 0;  
     CaloConeSelector selector5(0.5, geometry, DetId::Hcal);  
     selected = selector5.select(sc->eta(),sc->phi(),mhbhe);  
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) hcalEnergy += hit->energy();  
     hcalovecal05Phot[nPhot] = hcalEnergy/it->energy();  
     CaloConeSelector selector7(0.7, geometry, DetId::Hcal);
     selected = selector7.select(sc->eta(),sc->phi(),mhbhe);
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) hcalEnergy += hit->energy();
     hcalovecal07Phot[nPhot] = hcalEnergy/it->energy();
     hcalEnergy = 0;  
     CaloConeSelector selector10(1., geometry, DetId::Hcal); 
     selected = selector10.select(sc->eta(),sc->phi(),mhbhe); 
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) hcalEnergy += hit->energy(); 
     hcalovecal1Phot[nPhot] = hcalEnergy/it->energy(); 
     hcalEnergy = 0; 
 
     // calculate ECAL isolation 

     // ecal isolation with SC rechits removal
     SuperClusterHitsEcalIsolation scBasedIsolation(rhits,rhitsee);
     scBasedIsolation.setExtRadius(0.1);
     scBasedIsolation.excludeHalo(false);
     ecaliso01Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
     scBasedIsolation.setExtRadius(0.15);
     ecaliso015Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
     scBasedIsolation.setExtRadius(0.4);
     ecaliso04Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
     scBasedIsolation.setExtRadius(0.5);
     ecaliso05Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
     scBasedIsolation.setExtRadius(0.7);
     ecaliso07Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
     scBasedIsolation.setExtRadius(1.);
     ecaliso1Phot[nPhot]  = scBasedIsolation.getSum(iEvent,iSetup,&(*sc));
     // leave also old implementation
     double ecalEnergy = 0;
     CaloConeSelector selectorecal2(0.4, geometry, DetId::Ecal);
     selected=selectorecal2.select(it->eta(),it->phi(),mecalhits);
     for (CaloRecHitMetaCollectionV::const_iterator hit=selected->begin(); hit != selected->end(); ++hit) ecalEnergy += hit->energy();
     ecaliso04oldPhot[nPhot] = ecalEnergy;
     

     // cluster shape variables
     
     CaloClusterPtr tempCluster = it->superCluster()->seed();
     
     if(TMath::Abs(tempCluster->eta())<1.47){
       reco::ClusterShape tempShape=algo.Calculate(*tempCluster, rhits, &(*geometry_p), &(*topology_p),4.7);
       sMajMajPhot[nPhot]=tempShape.sMajMaj();
       sMinMinPhot[nPhot]=tempShape.sMinMin();
       FisherPhot[nPhot]=tempShape.fisher();
       alphaPhot[nPhot]=tempShape.s_alpha();
       sEtaEtaPhot[nPhot]=tempShape.sEtaEta();
       sPhiPhiPhot[nPhot]=tempShape.sPhiPhi();
       E1Phot[nPhot]=tempShape.eMax();
       E9Phot[nPhot]=tempShape.e3x3();
       E25Phot[nPhot]=tempShape.e5x5();
       // leave also old implementation
       tempShape=algo.Calculate(*tempCluster, rhits, &(*geometry_p), &(*topology_p),4.2);
       sMajMajoldPhot[nPhot]=tempShape.sMajMaj();
       sMinMinoldPhot[nPhot]=tempShape.sMinMin();
     }else{
       sMajMajPhot[nPhot]=-100.;
       sMinMinPhot[nPhot]=-100.;
       FisherPhot[nPhot]=-100.;
       alphaPhot[nPhot]=-100.;
       sEtaEtaPhot[nPhot]=-100.;
       sPhiPhiPhot[nPhot]=-100.;
       E1Phot[nPhot]=-100.;
       E9Phot[nPhot]=-100.;
       E25Phot[nPhot]=-100.;
       sMajMajoldPhot[nPhot]=-100.;
       sMinMinoldPhot[nPhot]=-100.;
     }

     nPhot++;
    
   }



   // Loop over gen Jets


   for ( GenJetCollection::const_iterator it = genJets_itCone5->begin(); 
	it != genJets_itCone5->end(); it++ ) {

     if(nJetGen_itCone5>=40) {cout << "Number of It. Cone 5 GEN jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetGen_itCone5[nJetGen_itCone5] = it->px();	 
     pyJetGen_itCone5[nJetGen_itCone5] = it->py();	 
     pzJetGen_itCone5[nJetGen_itCone5] = it->pz();	 
     eJetGen_itCone5[nJetGen_itCone5] = it->energy();	 
     etaJetGen_itCone5[nJetGen_itCone5] = it->eta();	 
     phiJetGen_itCone5[nJetGen_itCone5] = it->phi();	      
    
      nMuonsGen_itCone5[nJetGen_itCone5] = 0;
      eMuonsGen_itCone5[nJetGen_itCone5] = 0.;
     pxMuonsGen_itCone5[nJetGen_itCone5] = 0.;
     pyMuonsGen_itCone5[nJetGen_itCone5] = 0.;
     pzMuonsGen_itCone5[nJetGen_itCone5] = 0.;

      nElectronsGen_itCone5[nJetGen_itCone5] = 0;
      eElectronsGen_itCone5[nJetGen_itCone5] = 0.;
     pxElectronsGen_itCone5[nJetGen_itCone5] = 0.;
     pyElectronsGen_itCone5[nJetGen_itCone5] = 0.;
     pzElectronsGen_itCone5[nJetGen_itCone5] = 0.;

      nPhotonsGen_itCone5[nJetGen_itCone5] = 0;
      ePhotonsGen_itCone5[nJetGen_itCone5] = 0.;
     pxPhotonsGen_itCone5[nJetGen_itCone5] = 0.;
     pyPhotonsGen_itCone5[nJetGen_itCone5] = 0.;
     pzPhotonsGen_itCone5[nJetGen_itCone5] = 0.;

      nTracksGen_itCone5[nJetGen_itCone5] = 0;
      eTracksGen_itCone5[nJetGen_itCone5] = 0.;
     pxTracksGen_itCone5[nJetGen_itCone5] = 0.;
     pyTracksGen_itCone5[nJetGen_itCone5] = 0.;
     pzTracksGen_itCone5[nJetGen_itCone5] = 0.;

      nNeutralHadronsGen_itCone5[nJetGen_itCone5] = 0;
      eNeutralHadronsGen_itCone5[nJetGen_itCone5] = 0.;
     pxNeutralHadronsGen_itCone5[nJetGen_itCone5] = 0.;
     pyNeutralHadronsGen_itCone5[nJetGen_itCone5] = 0.;
     pzNeutralHadronsGen_itCone5[nJetGen_itCone5] = 0.;

     Float_t sumPt=0.;
 
     std::vector< const GenParticle* > jetParticles = it->getGenConstituents();
     std::vector< const GenParticle* > K0S;
     std::vector< const GenParticle* > Lambda;

     Float_t energy = it->energy();
     for( std::vector< const GenParticle* >::iterator iPart = jetParticles.begin();
       iPart != jetParticles.end(); ++iPart) {

       sumPt += (*iPart)->pt();

       Int_t partPdgId = (*iPart)->pdgId();
       Float_t partE = (*iPart)->energy();
       Float_t partPx = (*iPart)->px();
       Float_t partPy = (*iPart)->py();
       Float_t partPz = (*iPart)->pz();


       if( fabs(partPdgId)==13 ) { //muons
         nMuonsGen_itCone5[nJetGen_itCone5] += 1;
         eMuonsGen_itCone5[nJetGen_itCone5] += partE;
         pxMuonsGen_itCone5[nJetGen_itCone5] += partPx;
         pyMuonsGen_itCone5[nJetGen_itCone5] += partPy;
         pzMuonsGen_itCone5[nJetGen_itCone5] += partPz;
       } else if( fabs(partPdgId)==11 ) { //electrons
          nElectronsGen_itCone5[nJetGen_itCone5] += 1;
          eElectronsGen_itCone5[nJetGen_itCone5] += partE;
         pxElectronsGen_itCone5[nJetGen_itCone5] += partPx;
         pyElectronsGen_itCone5[nJetGen_itCone5] += partPy;
         pzElectronsGen_itCone5[nJetGen_itCone5] += partPz;
       } else if( (*iPart)->charge() != 0 ) {
          nTracksGen_itCone5[nJetGen_itCone5] += 1;
          eTracksGen_itCone5[nJetGen_itCone5] += partE;
         pxTracksGen_itCone5[nJetGen_itCone5] += partPx;
         pyTracksGen_itCone5[nJetGen_itCone5] += partPy;
         pzTracksGen_itCone5[nJetGen_itCone5] += partPz;
       } else if( partPdgId==22 ) { //photons
          nPhotonsGen_itCone5[nJetGen_itCone5] += 1;
          ePhotonsGen_itCone5[nJetGen_itCone5] += partE;
         pxPhotonsGen_itCone5[nJetGen_itCone5] += partPx;
         pyPhotonsGen_itCone5[nJetGen_itCone5] += partPy;
         pzPhotonsGen_itCone5[nJetGen_itCone5] += partPz;
       } else if( (fabs(partPdgId) != 12)&&(fabs(partPdgId) != 14)&&(fabs(partPdgId) != 16) ) { //neutrinos
         if( fabs(partPdgId)==310 )  K0S.push_back( *iPart );
         else if( fabs(partPdgId)==3122 ) Lambda.push_back( *iPart );
         else {
            nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
            eNeutralHadronsGen_itCone5[nJetGen_itCone5] += partE;
           pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += partPx;
           pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += partPy;
           pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += partPz;
         }
       }

     } //for jetParticles

     averagePartPtGen_itCone5[nJetGen_itCone5] = sumPt/(Float_t)jetParticles.size();

     Handle<SimTrackContainer> simTracks_h;
     const SimTrackContainer* simTracks;
     iEvent.getByLabel("g4SimHits", simTracks_h);
     simTracks = simTracks_h.product();

     Handle<SimVertexContainer> simVert_h;
     const SimVertexContainer* simVertices;
     iEvent.getByLabel("g4SimHits", simVert_h);
     simVertices = simVert_h.product();


     // ------------------ BEGIN K0Short 
     
     std::vector<int> K0S_trackId;

     for( int i=0; i<K0S.size(); ++i ) {

       K0S_trackId.push_back(-1);

       for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

         if(iSim->genpartIndex() != -1) {

          const reco::Candidate* p = &(*genParticles)[iSim->genpartIndex()-1];

          if( (K0S[i]->pdgId() != p->pdgId()) ||
              (K0S[i]->eta() != p->eta()) ||
              (K0S[i]->phi() != p->phi()) ||
              (K0S[i]->energy() != p->energy()) ) continue;
              
            K0S_trackId[i] = iSim->trackId();

          } //if genpartindex

       } //for sim tracks

     } //for K0s


     for( int i=0; i<K0S_trackId.size(); ++i ) {

       bool first=true;

       int vertIndex = 0;

       Float_t remainingE = K0S[i]->energy();
       Float_t remainingPx = K0S[i]->px();
       Float_t remainingPy = K0S[i]->py();
       Float_t remainingPz = K0S[i]->pz();

       for(SimVertexContainer::const_iterator iVert = simVertices->begin(); iVert != simVertices->end(); ++iVert) {

         if( K0S_trackId[i] == iVert->parentIndex() ) {

           for(SimTrackContainer::const_iterator iTrack=simTracks->begin(); iTrack!=simTracks->end(); ++iTrack) {
      
             if( iTrack->vertIndex()==vertIndex ) {
            
               Float_t vertX = iVert->position().x();
               Float_t vertY = iVert->position().y();
               Float_t vertZ = iVert->position().z();
            
               Float_t eta = iTrack->momentum().eta();
               Float_t theta = 2.*atan( exp(-eta));

               Float_t R = sqrt(vertX*vertX+vertY*vertY);
               if( (first)&&(fabs(eta)<1.4) ) {
                 h1_mfpK0S_itCone5->Fill(R*sin(theta));
                 first=false;
               }
            
               bool decayedBeforeCalo=false;

               if( fabs(eta)<1.4 ) {
                 if( R<129. ) decayedBeforeCalo=true;
               } else { 
                 if( fabs(vertZ) < 304. ) decayedBeforeCalo=true;
               }

               double trackE =  iTrack->momentum().e();
               double trackPx = iTrack->momentum().px();
               double trackPy = iTrack->momentum().py();
               double trackPz = iTrack->momentum().pz();
               double trackP = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz );
               double trackPt = trackP*sin(theta);

               if( (decayedBeforeCalo) && (iTrack->type()==111) ) { //pizeros
                   nPhotonsGen_itCone5[nJetGen_itCone5] += 2;  //boh
                   ePhotonsGen_itCone5[nJetGen_itCone5]  += trackE;
                   pxPhotonsGen_itCone5[nJetGen_itCone5] += trackPx;
                   pyPhotonsGen_itCone5[nJetGen_itCone5] += trackPy;
                   pzPhotonsGen_itCone5[nJetGen_itCone5] += trackPz;
               } else if( (R<30.) && (fabs(eta)<2.5) && (iTrack->charge()!=0) && (trackPt>0.1) ) { //min it tracking (4 step) pt = 100 mev? it's 75MeV, but it's fine like this
                    nTracksGen_itCone5[nJetGen_itCone5] += 1;
                    eTracksGen_itCone5[nJetGen_itCone5] += trackE;
                   pxTracksGen_itCone5[nJetGen_itCone5] += trackPx;
                   pyTracksGen_itCone5[nJetGen_itCone5] += trackPy;
                   pzTracksGen_itCone5[nJetGen_itCone5] += trackPz;
               } else {
                 nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
                 eNeutralHadronsGen_itCone5[nJetGen_itCone5]  += trackE;
                 pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPx;
                 pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPy;
                 pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPz;
               }

               remainingE  -= trackE;
               remainingPx -= trackPx;
               remainingPy -= trackPy;
               remainingPz -= trackPz;

             } //if vertindex

           } //for simtracks
         
         } // if found track ID

         ++vertIndex;

       } //for sim vertex

       if( remainingE > 0.001 ) {
          nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
          eNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingE;
         pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPx;
         pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPy;
         pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPz;
       }
         
     } //for K0S track ID

     // ------------------ END K0Short 

     // ------------------ BEGIN Lambda 
     
     std::vector<int> Lambda_trackId;

     for( int i=0; i<Lambda.size(); ++i ) {

       Lambda_trackId.push_back(-1);

       for(SimTrackContainer::const_iterator iSim = simTracks->begin(); iSim != simTracks->end(); ++iSim) {

         if(iSim->genpartIndex() != -1) {

          const reco::Candidate* p = &(*genParticles)[iSim->genpartIndex()-1];

          if( (Lambda[i]->pdgId() != p->pdgId()) ||
              (Lambda[i]->eta() != p->eta()) ||
              (Lambda[i]->phi() != p->phi()) ||
              (Lambda[i]->energy() != p->energy()) ) continue;
              
            Lambda_trackId[i] = iSim->trackId();

          } //if genpartindex

       } //for sim tracks

     } //for Lambda

     for( int i=0; i<Lambda_trackId.size(); ++i ) {

       bool first=true;

       int vertIndex = 0;

       Float_t remainingE = Lambda[i]->energy();
       Float_t remainingPx = Lambda[i]->px();
       Float_t remainingPy = Lambda[i]->py();
       Float_t remainingPz = Lambda[i]->pz();

       for(SimVertexContainer::const_iterator iVert = simVertices->begin(); iVert != simVertices->end(); ++iVert) {

         if( Lambda_trackId[i] == iVert->parentIndex() ) {

           for(SimTrackContainer::const_iterator iTrack=simTracks->begin(); iTrack!=simTracks->end(); ++iTrack) {
      
             if( iTrack->vertIndex()==vertIndex ) {
            
               Float_t vertX = iVert->position().x();
               Float_t vertY = iVert->position().y();
               Float_t vertZ = iVert->position().z();
            
               Float_t eta = iTrack->momentum().eta();
               Float_t theta = 2.*atan( exp(-eta));

               Float_t R = sqrt(vertX*vertX+vertY*vertY);
               if( (first)&&(fabs(eta)<1.4) ) {
                 h1_mfpLambda_itCone5->Fill(R*sin(theta));
                 first=false;
               }
            
               bool decayedBeforeCalo=false;
               if( fabs(eta)<1.4 ) {
                 if( R<129. ) decayedBeforeCalo=true;
               } else if( fabs(eta)<3. ) {
                 if( fabs(vertZ) < 304. ) decayedBeforeCalo=true;
               }

               double trackE  = iTrack->momentum().e();
               double trackPx = iTrack->momentum().px();
               double trackPy = iTrack->momentum().py();
               double trackPz = iTrack->momentum().pz();
               double trackP = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz );
               double trackPt = trackP*sin(theta);

               if( (decayedBeforeCalo) && (iTrack->type()==111) ) { //pizeros
                   nPhotonsGen_itCone5[nJetGen_itCone5] += 2;  //boh
                   ePhotonsGen_itCone5[nJetGen_itCone5]  += trackE;
                   pxPhotonsGen_itCone5[nJetGen_itCone5] += trackPx;
                   pyPhotonsGen_itCone5[nJetGen_itCone5] += trackPy;
                   pzPhotonsGen_itCone5[nJetGen_itCone5] += trackPz;
               } else if( (R<30.) && (fabs(eta)<2.5) && (iTrack->charge()!=0) && (trackPt>0.1) ) { //min it tracking (4 step) pt = 100 mev? it's 75MeV, but it's fine like this
                    nTracksGen_itCone5[nJetGen_itCone5] += 1;
                    eTracksGen_itCone5[nJetGen_itCone5] += trackE;
                   pxTracksGen_itCone5[nJetGen_itCone5] += trackPx;
                   pyTracksGen_itCone5[nJetGen_itCone5] += trackPy;
                   pzTracksGen_itCone5[nJetGen_itCone5] += trackPz;
               } else {
                 nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
                 eNeutralHadronsGen_itCone5[nJetGen_itCone5]  += trackE;
                 pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPx;
                 pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPy;
                 pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += trackPz;
               }

               remainingE  -= trackE;
               remainingPx -= trackPx;
               remainingPy -= trackPy;
               remainingPz -= trackPz;

             } //if vertindex

           } //for simtracks
         
         } // if found track ID

         ++vertIndex;

       } //for sim vertex

       if( remainingE > 0.001 ) {
          nNeutralHadronsGen_itCone5[nJetGen_itCone5] += 1;
          eNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingE;
         pxNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPx;
         pyNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPy;
         pzNeutralHadronsGen_itCone5[nJetGen_itCone5] += remainingPz;
       }
         
     } //for Lambda track ID

     nJetGen_itCone5++;

   }


/*
   for ( GenJetCollection::const_iterator it = genJets_sisCone5->begin(); 
	it != genJets_sisCone5->end(); it++ ) {

     if(nJetGen_sisCone5>=40) {cout << "Number of SIS Cone 5 GEN jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetGen_sisCone5[nJetGen_sisCone5] = it->px();	 
     pyJetGen_sisCone5[nJetGen_sisCone5] = it->py();	 
     pzJetGen_sisCone5[nJetGen_sisCone5] = it->pz();	 
     eJetGen_sisCone5[nJetGen_sisCone5] = it->energy();	 
     etaJetGen_sisCone5[nJetGen_sisCone5] = it->eta();	 
     phiJetGen_sisCone5[nJetGen_sisCone5] = it->phi();	      
     
 
      nMuonsGen_sisCone5[nJetGen_sisCone5] = 0;
      eMuonsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pxMuonsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pyMuonsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pzMuonsGen_sisCone5[nJetGen_sisCone5] = 0.;

      nElectronsGen_sisCone5[nJetGen_sisCone5] = 0;
      eElectronsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pxElectronsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pyElectronsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pzElectronsGen_sisCone5[nJetGen_sisCone5] = 0.;

      nPhotonsGen_sisCone5[nJetGen_sisCone5] = 0;
      ePhotonsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pxPhotonsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pyPhotonsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pzPhotonsGen_sisCone5[nJetGen_sisCone5] = 0.;

      nTracksGen_sisCone5[nJetGen_sisCone5] = 0;
      eTracksGen_sisCone5[nJetGen_sisCone5] = 0.;
     pxTracksGen_sisCone5[nJetGen_sisCone5] = 0.;
     pyTracksGen_sisCone5[nJetGen_sisCone5] = 0.;
     pzTracksGen_sisCone5[nJetGen_sisCone5] = 0.;

      nNeutralHadronsGen_sisCone5[nJetGen_sisCone5] = 0;
      eNeutralHadronsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pxNeutralHadronsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pyNeutralHadronsGen_sisCone5[nJetGen_sisCone5] = 0.;
     pzNeutralHadronsGen_sisCone5[nJetGen_sisCone5] = 0.;



     std::vector< const GenParticle* > jetParticles = it->getGenConstituents();

     for( std::vector< const GenParticle* >::iterator iPart = jetParticles.begin();
       iPart != jetParticles.end(); ++iPart) {

       Int_t partPdgId = (*iPart)->pdgId();
       Float_t partE = (*iPart)->energy();
       Float_t partPx = (*iPart)->px();
       Float_t partPy = (*iPart)->py();
       Float_t partPz = (*iPart)->pz();


       if( fabs(partPdgId)==13 ) { //muons
         nMuonsGen_sisCone5[nJetGen_sisCone5] += 1;
         eMuonsGen_sisCone5[nJetGen_sisCone5] += partE;
         pxMuonsGen_sisCone5[nJetGen_sisCone5] += partPx;
         pyMuonsGen_sisCone5[nJetGen_sisCone5] += partPy;
         pzMuonsGen_sisCone5[nJetGen_sisCone5] += partPz;
       } else if( fabs(partPdgId)==11 ) {
          nElectronsGen_sisCone5[nJetGen_sisCone5] += 1;
          eElectronsGen_sisCone5[nJetGen_sisCone5] += partE;
         pxElectronsGen_sisCone5[nJetGen_sisCone5] += partPx;
         pyElectronsGen_sisCone5[nJetGen_sisCone5] += partPy;
         pzElectronsGen_sisCone5[nJetGen_sisCone5] += partPz;
       } else if( (*iPart)->charge() != 0 ) {
          nTracksGen_sisCone5[nJetGen_sisCone5] += 1;
          eTracksGen_sisCone5[nJetGen_sisCone5] += partE;
         pxTracksGen_sisCone5[nJetGen_sisCone5] += partPx;
         pyTracksGen_sisCone5[nJetGen_sisCone5] += partPy;
         pzTracksGen_sisCone5[nJetGen_sisCone5] += partPz;
       } else if( partPdgId==22 ) { //photons
          nPhotonsGen_sisCone5[nJetGen_sisCone5] += 1;
          ePhotonsGen_sisCone5[nJetGen_sisCone5] += partE;
         pxPhotonsGen_sisCone5[nJetGen_sisCone5] += partPx;
         pyPhotonsGen_sisCone5[nJetGen_sisCone5] += partPy;
         pzPhotonsGen_sisCone5[nJetGen_sisCone5] += partPz;
       } else if( (fabs(partPdgId) != 12)&&(fabs(partPdgId) != 14)&&(fabs(partPdgId) != 16) ) { //neutrinos
          nNeutralHadronsGen_sisCone5[nJetGen_sisCone5] += 1;
          eNeutralHadronsGen_sisCone5[nJetGen_sisCone5] += partE;
         pxNeutralHadronsGen_sisCone5[nJetGen_sisCone5] += partPx;
         pyNeutralHadronsGen_sisCone5[nJetGen_sisCone5] += partPy;
         pzNeutralHadronsGen_sisCone5[nJetGen_sisCone5] += partPz;
       }

     } //for jetParticles

     nJetGen_sisCone5++;
   }



   for ( GenJetCollection::const_iterator it = genJets_kt4->begin(); 
	it != genJets_kt4->end(); it++ ) {

     if(nJetGen_kt4>=40) {cout << "Number of kT 4  GEN jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetGen_kt4[nJetGen_kt4] = it->px();	 
     pyJetGen_kt4[nJetGen_kt4] = it->py();	 
     pzJetGen_kt4[nJetGen_kt4] = it->pz();	 
     eJetGen_kt4[nJetGen_kt4] = it->energy();	 
     etaJetGen_kt4[nJetGen_kt4] = it->eta();	 
     phiJetGen_kt4[nJetGen_kt4] = it->phi();	      
     
      nMuonsGen_kt4[nJetGen_kt4] = 0;
      eMuonsGen_kt4[nJetGen_kt4] = 0.;
     pxMuonsGen_kt4[nJetGen_kt4] = 0.;
     pyMuonsGen_kt4[nJetGen_kt4] = 0.;
     pzMuonsGen_kt4[nJetGen_kt4] = 0.;

      nElectronsGen_kt4[nJetGen_kt4] = 0;
      eElectronsGen_kt4[nJetGen_kt4] = 0.;
     pxElectronsGen_kt4[nJetGen_kt4] = 0.;
     pyElectronsGen_kt4[nJetGen_kt4] = 0.;
     pzElectronsGen_kt4[nJetGen_kt4] = 0.;

      nPhotonsGen_kt4[nJetGen_kt4] = 0;
      ePhotonsGen_kt4[nJetGen_kt4] = 0.;
     pxPhotonsGen_kt4[nJetGen_kt4] = 0.;
     pyPhotonsGen_kt4[nJetGen_kt4] = 0.;
     pzPhotonsGen_kt4[nJetGen_kt4] = 0.;

      nTracksGen_kt4[nJetGen_kt4] = 0;
      eTracksGen_kt4[nJetGen_kt4] = 0.;
     pxTracksGen_kt4[nJetGen_kt4] = 0.;
     pyTracksGen_kt4[nJetGen_kt4] = 0.;
     pzTracksGen_kt4[nJetGen_kt4] = 0.;

      nNeutralHadronsGen_kt4[nJetGen_kt4] = 0;
      eNeutralHadronsGen_kt4[nJetGen_kt4] = 0.;
     pxNeutralHadronsGen_kt4[nJetGen_kt4] = 0.;
     pyNeutralHadronsGen_kt4[nJetGen_kt4] = 0.;
     pzNeutralHadronsGen_kt4[nJetGen_kt4] = 0.;


     std::vector< const GenParticle* > jetParticles = it->getGenConstituents();

     for( std::vector< const GenParticle* >::iterator iPart = jetParticles.begin();
       iPart != jetParticles.end(); ++iPart) {

       Int_t partPdgId = (*iPart)->pdgId();
       Float_t partE = (*iPart)->energy();
       Float_t partPx = (*iPart)->px();
       Float_t partPy = (*iPart)->py();
       Float_t partPz = (*iPart)->pz();


       if( fabs(partPdgId)==13 ) { //muons
         nMuonsGen_kt4[nJetGen_kt4] += 1;
         eMuonsGen_kt4[nJetGen_kt4] += partE;
         pxMuonsGen_kt4[nJetGen_kt4] += partPx;
         pyMuonsGen_kt4[nJetGen_kt4] += partPy;
         pzMuonsGen_kt4[nJetGen_kt4] += partPz;
       } else if( fabs(partPdgId)==11 ) {
          nElectronsGen_kt4[nJetGen_kt4] += 1;
          eElectronsGen_kt4[nJetGen_kt4] += partE;
         pxElectronsGen_kt4[nJetGen_kt4] += partPx;
         pyElectronsGen_kt4[nJetGen_kt4] += partPy;
         pzElectronsGen_kt4[nJetGen_kt4] += partPz;
       } else if( (*iPart)->charge() != 0 ) {
          nTracksGen_kt4[nJetGen_kt4] += 1;
          eTracksGen_kt4[nJetGen_kt4] += partE;
         pxTracksGen_kt4[nJetGen_kt4] += partPx;
         pyTracksGen_kt4[nJetGen_kt4] += partPy;
         pzTracksGen_kt4[nJetGen_kt4] += partPz;
       } else if( partPdgId==22 ) { //photons
          nPhotonsGen_kt4[nJetGen_kt4] += 1;
          ePhotonsGen_kt4[nJetGen_kt4] += partE;
         pxPhotonsGen_kt4[nJetGen_kt4] += partPx;
         pyPhotonsGen_kt4[nJetGen_kt4] += partPy;
         pzPhotonsGen_kt4[nJetGen_kt4] += partPz;
       } else if( (fabs(partPdgId) != 12)&&(fabs(partPdgId) != 14)&&(fabs(partPdgId) != 16) ) { //neutrinos
          nNeutralHadronsGen_kt4[nJetGen_kt4] += 1;
          eNeutralHadronsGen_kt4[nJetGen_kt4] += partE;
         pxNeutralHadronsGen_kt4[nJetGen_kt4] += partPx;
         pyNeutralHadronsGen_kt4[nJetGen_kt4] += partPy;
         pzNeutralHadronsGen_kt4[nJetGen_kt4] += partPz;
       }

     } //for jetParticles

     nJetGen_kt4++;
   }
*/



   // Loop over reco Jets

   // calo jets:
   
   for ( CaloJetCollection::const_iterator it = caloJets_itCone5->begin(); 
         it != caloJets_itCone5->end(); it++ ) {

     if(nJetCalo_itCone5>=40) {cout << "Number of It. Cone 5 calo jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetCalo_itCone5[nJetCalo_itCone5] = it->px();	 
     pyJetCalo_itCone5[nJetCalo_itCone5] = it->py();	 
     pzJetCalo_itCone5[nJetCalo_itCone5] = it->pz();	 
     eJetCalo_itCone5[nJetCalo_itCone5] = it->energy();	 
     etaJetCalo_itCone5[nJetCalo_itCone5] = it->eta();	 
     phiJetCalo_itCone5[nJetCalo_itCone5] = it->phi();	      
     
     nJetCalo_itCone5++;
   }

/*
   for ( CaloJetCollection::const_iterator it = caloJets_sisCone5->begin(); 
         it != caloJets_sisCone5->end(); it++ ) {

     if(nJetCalo_sisCone5>=40) {cout << "Number of SIS Cone 5 calo jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetCalo_sisCone5[nJetCalo_sisCone5] = it->px();	 
     pyJetCalo_sisCone5[nJetCalo_sisCone5] = it->py();	 
     pzJetCalo_sisCone5[nJetCalo_sisCone5] = it->pz();	 
     eJetCalo_sisCone5[nJetCalo_sisCone5] = it->energy();	 
     etaJetCalo_sisCone5[nJetCalo_sisCone5] = it->eta();	 
     phiJetCalo_sisCone5[nJetCalo_sisCone5] = it->phi();	      
     
     nJetCalo_sisCone5++;
   }



   for ( CaloJetCollection::const_iterator it = caloJets_kt4->begin(); 
         it != caloJets_kt4->end(); it++ ) {

     if(nJetCalo_kt4>=40) {cout << "Number of kT 4 calo jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetCalo_kt4[nJetCalo_kt4] = it->px();	 
     pyJetCalo_kt4[nJetCalo_kt4] = it->py();	 
     pzJetCalo_kt4[nJetCalo_kt4] = it->pz();	 
     eJetCalo_kt4[nJetCalo_kt4] = it->energy();	 
     etaJetCalo_kt4[nJetCalo_kt4] = it->eta();	 
     phiJetCalo_kt4[nJetCalo_kt4] = it->phi();	      
     
     nJetCalo_kt4++;
   }
*/


   // PF Jets:
//std::cout << std::endl << "*** NEW EVENT" << std::endl;

   for ( PFJetCollection::const_iterator it = pfJets_itCone5->begin(); 
         it != pfJets_itCone5->end(); it++ ) {

     if(nJetPF_itCone5>=40) {cout << "Number of It. Cone 5 PF jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetPF_itCone5[nJetPF_itCone5] = it->px();	 
     pyJetPF_itCone5[nJetPF_itCone5] = it->py();	 
     pzJetPF_itCone5[nJetPF_itCone5] = it->pz();	 
     eJetPF_itCone5[nJetPF_itCone5] = it->energy();	 
     etaJetPF_itCone5[nJetPF_itCone5] = it->eta();	 
     phiJetPF_itCone5[nJetPF_itCone5] = it->phi();	      


     Int_t nChargedHadrons = 0;
     Int_t nPhotons = 0;
     Int_t nNeutralHadrons = 0;
     Int_t nElectrons = 0;
     Int_t nMuons = 0;
     Int_t nHFHadrons = 0;
     Int_t nHFEM = 0;

     Float_t eChargedHadrons = 0.;
     Float_t ePhotons = 0.;
     Float_t eNeutralHadrons = 0.;
     Float_t eElectrons = 0.;
     Float_t eMuons = 0.;
     Float_t eHFHadrons = 0.;
     Float_t eHFEM = 0.;

     Float_t pxChargedHadrons = 0.;
     Float_t pxPhotons = 0.;
     Float_t pxNeutralHadrons = 0.;
     Float_t pxElectrons = 0.;
     Float_t pxMuons = 0.;
     Float_t pxHFHadrons = 0.;
     Float_t pxHFEM = 0.;
    
     Float_t pyChargedHadrons = 0.;
     Float_t pyPhotons = 0.;
     Float_t pyNeutralHadrons = 0.;
     Float_t pyElectrons = 0.;
     Float_t pyMuons = 0.;
     Float_t pyHFHadrons = 0.;
     Float_t pyHFEM = 0.;
    
     Float_t pzChargedHadrons = 0.;
     Float_t pzPhotons = 0.;
     Float_t pzNeutralHadrons = 0.;
     Float_t pzElectrons = 0.;
     Float_t pzMuons = 0.;
     Float_t pzHFHadrons = 0.;
     Float_t pzHFEM = 0.;
    
     std::vector<const PFCandidate*> pfCandidates = it->getPFConstituents();

     
     Float_t sumPt=0.;

     for( std::vector<const PFCandidate*>::const_iterator iCandidate = pfCandidates.begin();
          iCandidate!=pfCandidates.end(); ++iCandidate ) {


       int particleType = (*iCandidate)->particleId();
       Float_t energy = (*iCandidate)->energy();
       Float_t px = (*iCandidate)->px();
       Float_t py = (*iCandidate)->py();
       Float_t pz = (*iCandidate)->pz();

       Float_t p = sqrt(px*px + py*py + pz*pz);
       Float_t theta = (*iCandidate)->theta();
       Float_t pt = p*sin(theta);

       sumPt += pt;
      
       if( particleType==1 ) {
         nChargedHadrons += 1;
         eChargedHadrons += energy;
         pxChargedHadrons += px;
         pyChargedHadrons += py;
         pzChargedHadrons += pz;
       }

       if( particleType==2 ) {
         nElectrons += 1;
          eElectrons += energy;
         pxElectrons += px;
         pyElectrons += py;
         pzElectrons += pz;
       }

       if( particleType==3 ) {
         nMuons += 1;
          eMuons += energy;
         pxMuons += px;
         pyMuons += py;
         pzMuons += pz;
       }

       if( particleType==4 ) {
         nPhotons += 1;
          ePhotons += energy;
         pxPhotons += px;
         pyPhotons += py;
         pzPhotons += pz;
       }

       if( particleType==5 ) {
         nNeutralHadrons += 1;
          eNeutralHadrons += energy;
         pxNeutralHadrons += px;
         pyNeutralHadrons += py;
         pzNeutralHadrons += pz;
       }

       if( particleType==6 ) {
         nHFHadrons += 1;
          eHFHadrons += energy;
         pxHFHadrons += px;
         pyHFHadrons += py;
         pzHFHadrons += pz;
       }

       if( particleType==7 ) {
         nHFEM += 1;
          eHFEM += energy;
         pxHFEM += px;
         pyHFEM += py;
         pzHFEM += pz;
       }

     } //for PFCandidates

     averagePartPtPF_itCone5[nJetPF_itCone5] = sumPt/pfCandidates.size();

      nChargedHadronsPF_itCone5[nJetPF_itCone5] =  nChargedHadrons;
      eChargedHadronsPF_itCone5[nJetPF_itCone5] =  eChargedHadrons;
     pxChargedHadronsPF_itCone5[nJetPF_itCone5] = pxChargedHadrons;
     pyChargedHadronsPF_itCone5[nJetPF_itCone5] = pyChargedHadrons;
     pzChargedHadronsPF_itCone5[nJetPF_itCone5] = pzChargedHadrons;
 
      nElectronsPF_itCone5[nJetPF_itCone5] =  nElectrons;
      eElectronsPF_itCone5[nJetPF_itCone5] =  eElectrons;
     pxElectronsPF_itCone5[nJetPF_itCone5] = pxElectrons;
     pyElectronsPF_itCone5[nJetPF_itCone5] = pyElectrons;
     pzElectronsPF_itCone5[nJetPF_itCone5] = pzElectrons;
 
      nMuonsPF_itCone5[nJetPF_itCone5] =  nMuons;
      eMuonsPF_itCone5[nJetPF_itCone5] =  eMuons;
     pxMuonsPF_itCone5[nJetPF_itCone5] = pxMuons;
     pyMuonsPF_itCone5[nJetPF_itCone5] = pyMuons;
     pzMuonsPF_itCone5[nJetPF_itCone5] = pzMuons;
 
      nPhotonsPF_itCone5[nJetPF_itCone5] =  nPhotons;
      ePhotonsPF_itCone5[nJetPF_itCone5] =  ePhotons;
     pxPhotonsPF_itCone5[nJetPF_itCone5] = pxPhotons;
     pyPhotonsPF_itCone5[nJetPF_itCone5] = pyPhotons;
     pzPhotonsPF_itCone5[nJetPF_itCone5] = pzPhotons;
 
      nNeutralHadronsPF_itCone5[nJetPF_itCone5] =  nNeutralHadrons;
      eNeutralHadronsPF_itCone5[nJetPF_itCone5] =  eNeutralHadrons;
     pxNeutralHadronsPF_itCone5[nJetPF_itCone5] = pxNeutralHadrons;
     pyNeutralHadronsPF_itCone5[nJetPF_itCone5] = pyNeutralHadrons;
     pzNeutralHadronsPF_itCone5[nJetPF_itCone5] = pzNeutralHadrons;
 
      nHFHadronsPF_itCone5[nJetPF_itCone5] =  nHFHadrons;
      eHFHadronsPF_itCone5[nJetPF_itCone5] =  eHFHadrons;
     pxHFHadronsPF_itCone5[nJetPF_itCone5] = pxHFHadrons;
     pyHFHadronsPF_itCone5[nJetPF_itCone5] = pyHFHadrons;
     pzHFHadronsPF_itCone5[nJetPF_itCone5] = pzHFHadrons;
 
      nHFEMPF_itCone5[nJetPF_itCone5] =  nHFEM;
      eHFEMPF_itCone5[nJetPF_itCone5] =  eHFEM;
     pxHFEMPF_itCone5[nJetPF_itCone5] = pxHFEM;
     pyHFEMPF_itCone5[nJetPF_itCone5] = pyHFEM;
     pzHFEMPF_itCone5[nJetPF_itCone5] = pzHFEM;
 
     nJetPF_itCone5++;

   } //for PFJets

/*
   for ( PFJetCollection::const_iterator it = pfJets_sisCone5->begin(); 
         it != pfJets_sisCone5->end(); it++ ) {

     if(nJetPF_sisCone5>=40) {cout << "Number of SIS Cone 5 PF jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetPF_sisCone5[nJetPF_sisCone5] = it->px();	 
     pyJetPF_sisCone5[nJetPF_sisCone5] = it->py();	 
     pzJetPF_sisCone5[nJetPF_sisCone5] = it->pz();	 
     eJetPF_sisCone5[nJetPF_sisCone5] = it->energy();	 
     etaJetPF_sisCone5[nJetPF_sisCone5] = it->eta();	 
     phiJetPF_sisCone5[nJetPF_sisCone5] = it->phi();	      
     chargedHadronEnergyJetPF_sisCone5[nJetPF_sisCone5] = it->chargedHadronEnergy();	      
     neutralHadronEnergyJetPF_sisCone5[nJetPF_sisCone5] = it->neutralHadronEnergy();	      
     chargedEmEnergyJetPF_sisCone5[nJetPF_sisCone5] = it->chargedEmEnergy();	      
     neutralEmEnergyJetPF_sisCone5[nJetPF_sisCone5] = it->neutralEmEnergy();	      
     chargedMuEnergyJetPF_sisCone5[nJetPF_sisCone5] = it->chargedMuEnergy();	      
     chargedMultiplicityJetPF_sisCone5[nJetPF_sisCone5] = it->chargedMultiplicity();	      
     neutralMultiplicityJetPF_sisCone5[nJetPF_sisCone5] = it->neutralMultiplicity();	      
     muonMultiplicityJetPF_sisCone5[nJetPF_sisCone5] = it->muonMultiplicity();	      
     
     nJetPF_sisCone5++;
   }


   for ( PFJetCollection::const_iterator it = pfJets_kt4->begin(); 
         it != pfJets_kt4->end(); it++ ) {

     if(nJetPF_kt4>=40) {cout << "Number of kT 4 PF jets is larger than 40. Skipping" << endl; continue;}
     if(it->pt()<5.) continue;
     pxJetPF_kt4[nJetPF_kt4] = it->px();	 
     pyJetPF_kt4[nJetPF_kt4] = it->py();	 
     pzJetPF_kt4[nJetPF_kt4] = it->pz();	 
     eJetPF_kt4[nJetPF_kt4] = it->energy();	 
     etaJetPF_kt4[nJetPF_kt4] = it->eta();	 
     phiJetPF_kt4[nJetPF_kt4] = it->phi();	      
     chargedHadronEnergyJetPF_kt4[nJetPF_kt4] = it->chargedHadronEnergy();	      
     neutralHadronEnergyJetPF_kt4[nJetPF_kt4] = it->neutralHadronEnergy();	      
     chargedEmEnergyJetPF_kt4[nJetPF_kt4] = it->chargedEmEnergy();	      
     neutralEmEnergyJetPF_kt4[nJetPF_kt4] = it->neutralEmEnergy();	      
     chargedMuEnergyJetPF_kt4[nJetPF_kt4] = it->chargedMuEnergy();	      
     chargedMultiplicityJetPF_kt4[nJetPF_kt4] = it->chargedMultiplicity();	      
     neutralMultiplicityJetPF_kt4[nJetPF_kt4] = it->neutralMultiplicity();	      
     muonMultiplicityJetPF_kt4[nJetPF_kt4] = it->muonMultiplicity();	      
     
     nJetPF_kt4++;
   }

*/


  
   event++;  
   m_tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
JetAnalyzer::beginJob(const edm::EventSetup&)
{


  h1_mfpK0S_itCone5 = fs_->make<TH1F>("mfpK0S_itCone5","", 100, 0., 170.);
  h1_mfpLambda_itCone5 = fs_->make<TH1F>("mfpLambda_itCone5","", 100, 0., 170.);

  m_tree = fs_->make<TTree>("pippo","Analysis tree");

  m_tree->Branch("event",&event,"event/I");
  m_tree->Branch("ptHat",&ptHat,"ptHat/F");
  m_tree->Branch("nMC",&nMC,"nMC/I");
  m_tree->Branch("pdgIdMC",pdgIdMC,"pdgIdMC[nMC]/I");
  m_tree->Branch("statusMC",statusMC,"statusMC[nMC]/I");
  m_tree->Branch("pxMC ",pxMC ,"pxMC[nMC]/F");
  m_tree->Branch("pyMC ",pyMC ,"pyMC[nMC]/F");
  m_tree->Branch("pzMC ",pzMC ,"pzMC[nMC]/F");
  m_tree->Branch("eMC  ",eMC  ,"eMC[nMC]/F");
  m_tree->Branch("etaMC",etaMC,"etaMC[nMC]/F");
  m_tree->Branch("phiMC",phiMC,"phiMC[nMC]/F");

  m_tree->Branch("nPhot",&nPhot,"nPhot/I");
  m_tree->Branch("pxPhot ",&pxPhot ,"pxPhot[nPhot]/F");
  m_tree->Branch("pyPhot ",&pyPhot ,"pyPhot[nPhot]/F");
  m_tree->Branch("pzPhot ",&pzPhot ,"pzPhot[nPhot]/F");
  m_tree->Branch("ePhot  ",&ePhot  ,"ePhot[nPhot]/F");
  m_tree->Branch("escPhot  ",&escPhot  ,"escPhot[nPhot]/F");
  m_tree->Branch("etaPhot",&etaPhot,"etaPhot[nPhot]/F");
  m_tree->Branch("phiPhot",&phiPhot,"phiPhot[nPhot]/F");

  m_tree->Branch("ptiso0015Phot",&ptiso0015Phot,"ptiso0015Phot[nPhot]/F");
  m_tree->Branch("ntrkiso0015Phot",&ntrkiso0015Phot,"ntrkiso0015Phot[nPhot]/I");
  m_tree->Branch("ptiso035Phot",&ptiso035Phot,"ptiso035Phot[nPhot]/F");
  m_tree->Branch("ntrkiso035Phot",&ntrkiso035Phot,"ntrkiso035Phot[nPhot]/I");
  m_tree->Branch("ptiso05Phot",&ptiso05Phot,"ptiso05Phot[nPhot]/F");
  m_tree->Branch("ntrkiso05Phot",&ntrkiso05Phot,"ntrkiso05Phot[nPhot]/I");
  m_tree->Branch("ptiso07Phot",&ptiso07Phot,"ptiso07Phot[nPhot]/F");
  m_tree->Branch("ntrkiso07Phot",&ntrkiso07Phot,"ntrkiso07Phot[nPhot]/I");
  m_tree->Branch("ptiso1Phot",&ptiso1Phot,"ptiso1Phot[nPhot]/F");
  m_tree->Branch("ntrkiso1Phot",&ntrkiso1Phot,"ntrkiso1Phot[nPhot]/I");
//   m_tree->Branch("ptisoatecal02Phot",&ptisoatecal02Phot,"ptisoatecal02Phot[nPhot]/F");
//   m_tree->Branch("ntrkisoatecal02Phot",&ntrkisoatecal02Phot,"ntrkisoatecal02Phot[nPhot]/I");
//   m_tree->Branch("ptisoatecal025Phot",&ptisoatecal02Phot,"ptisoatecal025Phot[nPhot]/F");
//   m_tree->Branch("ntrkisoatecal025Phot",&ntrkisoatecal02Phot,"ntrkisoatecal025Phot[nPhot]/I");
//   m_tree->Branch("ptisoatecal03Phot",&ptisoatecal02Phot,"ptisoatecal03Phot[nPhot]/F");
//   m_tree->Branch("ntrkisoatecal03Phot",&ntrkisoatecal03Phot,"ntrkisoatecal03Phot[nPhot]/I");
//   m_tree->Branch("ptisoatecal035Phot",&ptisoatecal035Phot,"ptisoatecal035Phot[nPhot]/F");
//   m_tree->Branch("ntrkisoatecal035Phot",&ntrkisoatecal035Phot,"ntrkisoatecal035Phot[nPhot]/I");
//   m_tree->Branch("ptisoatecal04Phot",&ptisoatecal04Phot,"ptisoatecal04Phot[nPhot]/F");
  m_tree->Branch("hcalovecal01Phot",&hcalovecal01Phot,"hcalovecal01Phot[nPhot]/F");
  m_tree->Branch("hcalovecal015Phot",&hcalovecal015Phot,"hcalovecal015Phot[nPhot]/F");
  m_tree->Branch("hcalovecal04Phot",&hcalovecal04Phot,"hcalovecal04Phot[nPhot]/F"); 
  m_tree->Branch("hcalovecal05Phot",&hcalovecal05Phot,"hcalovecal05Phot[nPhot]/F"); 
  m_tree->Branch("hcalovecal07Phot",&hcalovecal07Phot,"hcalovecal07Phot[nPhot]/F");
  m_tree->Branch("hcalovecal1Phot",&hcalovecal1Phot,"hcalovecal1Phot[nPhot]/F"); 
  m_tree->Branch("ecaliso01Phot",&ecaliso01Phot,"ecaliso01Phot[nPhot]/F"); 
  m_tree->Branch("ecaliso015Phot",&ecaliso015Phot,"ecaliso015Phot[nPhot]/F");  
  m_tree->Branch("ecaliso04Phot",&ecaliso04Phot,"ecaliso04Phot[nPhot]/F");  
  m_tree->Branch("ecaliso05Phot",&ecaliso05Phot,"ecaliso05Phot[nPhot]/F");  
  m_tree->Branch("ecaliso07Phot",&ecaliso07Phot,"ecaliso07Phot[nPhot]/F"); 
  m_tree->Branch("ecaliso1Phot",&ecaliso1Phot,"ecaliso1Phot[nPhot]/F");  
  m_tree->Branch("ecaliso04oldPhot",&ecaliso04oldPhot,"ecaliso04oldPhot[nPhot]/F");  
  m_tree->Branch("LATPhot",&LATPhot,"LATPhot[nPhot]/F");
  m_tree->Branch("sMajMajPhot",&sMajMajPhot,"sMajMaj2Phot[nPhot]/F");
  m_tree->Branch("sMinMinPhot",&sMinMinPhot,"sMinMin2Phot[nPhot]/F");
  m_tree->Branch("alphaPhot",&alphaPhot,"alphaPhot[nPhot]/F");
  m_tree->Branch("sEtaEtaPhot",&sEtaEtaPhot,"sEtaEtaPhot[nPhot]/F");
  m_tree->Branch("sPhiPhiPhot",&sPhiPhiPhot,"sPhiPhiPhot[nPhot]/F");
  m_tree->Branch("E1Phot",&E1Phot,"E1Phot[nPhot]/F");
  m_tree->Branch("E9Phot",&E9Phot,"E9Phot[nPhot]/F");
  m_tree->Branch("E25Phot",&E25Phot,"E25Phot[nPhot]/F");
  m_tree->Branch("FisherPhot",&FisherPhot,"FisherPhot[nPhot]/F");
  m_tree->Branch("sMajMajoldPhot",&sMajMajoldPhot,"sMajMaj2oldPhot[nPhot]/F");
  m_tree->Branch("sMinMinoldPhot",&sMinMinoldPhot,"sMinMin2oldPhot[nPhot]/F");

  m_tree->Branch("nJetGen_itCone5",&nJetGen_itCone5,"nJetGen_itCone5/I");
  m_tree->Branch("pxJetGen_itCone5",pxJetGen_itCone5,"pxJetGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pyJetGen_itCone5",pyJetGen_itCone5,"pyJetGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pzJetGen_itCone5",pzJetGen_itCone5,"pzJetGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("eJetGen_itCone5",eJetGen_itCone5,"eJetGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("etaJetGen_itCone5",etaJetGen_itCone5,"etaJetGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("phiJetGen_itCone5",phiJetGen_itCone5,"phiJetGen_itCone5[nJetGen_itCone5]/F");

  m_tree->Branch("nTracksGen_itCone5",        nTracksGen_itCone5,        "nTracksGen_itCone5[nJetGen_itCone5]/I");
  m_tree->Branch("nPhotonsGen_itCone5",       nPhotonsGen_itCone5,       "nPhotonsGen_itCone5[nJetGen_itCone5]/I");
  m_tree->Branch("nElectronsGen_itCone5",     nElectronsGen_itCone5,     "nElectronsGen_itCone5[nJetGen_itCone5]/I");
  m_tree->Branch("nMuonsGen_itCone5",         nMuonsGen_itCone5,         "nMuonsGen_itCone5[nJetGen_itCone5]/I");
  m_tree->Branch("nNeutralHadronsGen_itCone5",nNeutralHadronsGen_itCone5,"nNeutralHadronsGen_itCone5[nJetGen_itCone5]/I");

  m_tree->Branch("eTracksGen_itCone5",        eTracksGen_itCone5,        "eTracksGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("ePhotonsGen_itCone5",       ePhotonsGen_itCone5,       "ePhotonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("eElectronsGen_itCone5",     eElectronsGen_itCone5,     "eElectronsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("eMuonsGen_itCone5",         eMuonsGen_itCone5,         "eMuonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("eNeutralHadronsGen_itCone5",eNeutralHadronsGen_itCone5,"eNeutralHadronsGen_itCone5[nJetGen_itCone5]/F");

  m_tree->Branch("pxTracksGen_itCone5",        pxTracksGen_itCone5,        "pxTracksGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pxPhotonsGen_itCone5",       pxPhotonsGen_itCone5,       "pxPhotonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pxElectronsGen_itCone5",     pxElectronsGen_itCone5,     "pxElectronsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pxMuonsGen_itCone5",         pxMuonsGen_itCone5,         "pxMuonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pxNeutralHadronsGen_itCone5",pxNeutralHadronsGen_itCone5,"pxNeutralHadronsGen_itCone5[nJetGen_itCone5]/F");

  m_tree->Branch("pyTracksGen_itCone5",        pyTracksGen_itCone5,        "pyTracksGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pyPhotonsGen_itCone5",       pyPhotonsGen_itCone5,       "pyPhotonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pyElectronsGen_itCone5",     pyElectronsGen_itCone5,     "pyElectronsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pyMuonsGen_itCone5",         pyMuonsGen_itCone5,         "pyMuonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pyNeutralHadronsGen_itCone5",pyNeutralHadronsGen_itCone5,"pyNeutralHadronsGen_itCone5[nJetGen_itCone5]/F");

  m_tree->Branch("pzTracksGen_itCone5",        pzTracksGen_itCone5,        "pzTracksGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pzPhotonsGen_itCone5",       pzPhotonsGen_itCone5,       "pzPhotonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pzElectronsGen_itCone5",     pzElectronsGen_itCone5,     "pzElectronsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pzMuonsGen_itCone5",         pzMuonsGen_itCone5,         "pzMuonsGen_itCone5[nJetGen_itCone5]/F");
  m_tree->Branch("pzNeutralHadronsGen_itCone5",pzNeutralHadronsGen_itCone5,"pzNeutralHadronsGen_itCone5[nJetGen_itCone5]/F");

/*
  m_tree->Branch("nJetGen_sisCone5",&nJetGen_sisCone5,"nJetGen_sisCone5/I");
  m_tree->Branch("pxJetGen_sisCone5",pxJetGen_sisCone5,"pxJetGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pyJetGen_sisCone5",pyJetGen_sisCone5,"pyJetGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pzJetGen_sisCone5",pzJetGen_sisCone5,"pzJetGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("eJetGen_sisCone5",eJetGen_sisCone5,"eJetGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("etaJetGen_sisCone5",etaJetGen_sisCone5,"etaJetGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("phiJetGen_sisCone5",phiJetGen_sisCone5,"phiJetGen_sisCone5[nJetGen_sisCone5]/F");

  m_tree->Branch("nTracksGen_sisCone5",        nTracksGen_sisCone5,        "nTracksGen_sisCone5[nJetGen_sisCone5]/I");
  m_tree->Branch("nPhotonsGen_sisCone5",       nPhotonsGen_sisCone5,       "nPhotonsGen_sisCone5[nJetGen_sisCone5]/I");
  m_tree->Branch("nElectronsGen_sisCone5",     nElectronsGen_sisCone5,     "nElectronsGen_sisCone5[nJetGen_sisCone5]/I");
  m_tree->Branch("nMuonsGen_sisCone5",         nMuonsGen_sisCone5,         "nMuonsGen_sisCone5[nJetGen_sisCone5]/I");
  m_tree->Branch("nNeutralHadronsGen_sisCone5",nNeutralHadronsGen_sisCone5,"nNeutralHadronsGen_sisCone5[nJetGen_sisCone5]/I");

  m_tree->Branch("eTracksGen_sisCone5",        eTracksGen_sisCone5,        "eTracksGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("ePhotonsGen_sisCone5",       ePhotonsGen_sisCone5,       "ePhotonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("eElectronsGen_sisCone5",     eElectronsGen_sisCone5,     "eElectronsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("eMuonsGen_sisCone5",         eMuonsGen_sisCone5,         "eMuonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("eNeutralHadronsGen_sisCone5",eNeutralHadronsGen_sisCone5,"eNeutralHadronsGen_sisCone5[nJetGen_sisCone5]/F");

  m_tree->Branch("pxTracksGen_sisCone5",        pxTracksGen_sisCone5,        "pxTracksGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pxPhotonsGen_sisCone5",       pxPhotonsGen_sisCone5,       "pxPhotonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pxElectronsGen_sisCone5",     pxElectronsGen_sisCone5,     "pxElectronsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pxMuonsGen_sisCone5",         pxMuonsGen_sisCone5,         "pxMuonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pxNeutralHadronsGen_sisCone5",pxNeutralHadronsGen_sisCone5,"pxNeutralHadronsGen_sisCone5[nJetGen_sisCone5]/F");

  m_tree->Branch("pyTracksGen_sisCone5",        pyTracksGen_sisCone5,        "pyTracksGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pyPhotonsGen_sisCone5",       pyPhotonsGen_sisCone5,       "pyPhotonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pyElectronsGen_sisCone5",     pyElectronsGen_sisCone5,     "pyElectronsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pyMuonsGen_sisCone5",         pyMuonsGen_sisCone5,         "pyMuonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pyNeutralHadronsGen_sisCone5",pyNeutralHadronsGen_sisCone5,"pyNeutralHadronsGen_sisCone5[nJetGen_sisCone5]/F");

  m_tree->Branch("pzTracksGen_sisCone5",        pzTracksGen_sisCone5,        "pzTracksGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pzPhotonsGen_sisCone5",       pzPhotonsGen_sisCone5,       "pzPhotonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pzElectronsGen_sisCone5",     pzElectronsGen_sisCone5,     "pzElectronsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pzMuonsGen_sisCone5",         pzMuonsGen_sisCone5,         "pzMuonsGen_sisCone5[nJetGen_sisCone5]/F");
  m_tree->Branch("pzNeutralHadronsGen_sisCone5",pzNeutralHadronsGen_sisCone5,"pzNeutralHadronsGen_sisCone5[nJetGen_sisCone5]/F");


  m_tree->Branch("nJetGen_kt4",&nJetGen_kt4,"nJetGen_kt4/I");
  m_tree->Branch("pxJetGen_kt4",pxJetGen_kt4,"pxJetGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pyJetGen_kt4",pyJetGen_kt4,"pyJetGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pzJetGen_kt4",pzJetGen_kt4,"pzJetGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("eJetGen_kt4",eJetGen_kt4,"eJetGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("etaJetGen_kt4",etaJetGen_kt4,"etaJetGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("phiJetGen_kt4",phiJetGen_kt4,"phiJetGen_kt4[nJetGen_kt4]/F");

  m_tree->Branch("nTracksGen_kt4",        nTracksGen_kt4,        "nTracksGen_kt4[nJetGen_kt4]/I");
  m_tree->Branch("nPhotonsGen_kt4",       nPhotonsGen_kt4,       "nPhotonsGen_kt4[nJetGen_kt4]/I");
  m_tree->Branch("nElectronsGen_kt4",     nElectronsGen_kt4,     "nElectronsGen_kt4[nJetGen_kt4]/I");
  m_tree->Branch("nMuonsGen_kt4",         nMuonsGen_kt4,         "nMuonsGen_kt4[nJetGen_kt4]/I");
  m_tree->Branch("nNeutralHadronsGen_kt4",nNeutralHadronsGen_kt4,"nNeutralHadronsGen_kt4[nJetGen_kt4]/I");

  m_tree->Branch("eTracksGen_kt4",        eTracksGen_kt4,        "eTracksGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("ePhotonsGen_kt4",       ePhotonsGen_kt4,       "ePhotonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("eElectronsGen_kt4",     eElectronsGen_kt4,     "eElectronsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("eMuonsGen_kt4",         eMuonsGen_kt4,         "eMuonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("eNeutralHadronsGen_kt4",eNeutralHadronsGen_kt4,"eNeutralHadronsGen_kt4[nJetGen_kt4]/F");

  m_tree->Branch("pxTracksGen_kt4",        pxTracksGen_kt4,        "pxTracksGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pxPhotonsGen_kt4",       pxPhotonsGen_kt4,       "pxPhotonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pxElectronsGen_kt4",     pxElectronsGen_kt4,     "pxElectronsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pxMuonsGen_kt4",         pxMuonsGen_kt4,         "pxMuonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pxNeutralHadronsGen_kt4",pxNeutralHadronsGen_kt4,"pxNeutralHadronsGen_kt4[nJetGen_kt4]/F");

  m_tree->Branch("pyTracksGen_kt4",        pyTracksGen_kt4,        "pyTracksGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pyPhotonsGen_kt4",       pyPhotonsGen_kt4,       "pyPhotonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pyElectronsGen_kt4",     pyElectronsGen_kt4,     "pyElectronsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pyMuonsGen_kt4",         pyMuonsGen_kt4,         "pyMuonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pyNeutralHadronsGen_kt4",pyNeutralHadronsGen_kt4,"pyNeutralHadronsGen_kt4[nJetGen_kt4]/F");

  m_tree->Branch("pzTracksGen_kt4",        pzTracksGen_kt4,        "pzTracksGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pzPhotonsGen_kt4",       pzPhotonsGen_kt4,       "pzPhotonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pzElectronsGen_kt4",     pzElectronsGen_kt4,     "pzElectronsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pzMuonsGen_kt4",         pzMuonsGen_kt4,         "pzMuonsGen_kt4[nJetGen_kt4]/F");
  m_tree->Branch("pzNeutralHadronsGen_kt4",pzNeutralHadronsGen_kt4,"pzNeutralHadronsGen_kt4[nJetGen_kt4]/F");
*/

  m_tree->Branch("nJetCalo_itCone5",&nJetCalo_itCone5,"nJetCalo_itCone5/I");
  m_tree->Branch("pxJetCalo_itCone5",pxJetCalo_itCone5,"pxJetCalo_itCone5[nJetCalo_itCone5]/F");
  m_tree->Branch("pyJetCalo_itCone5",pyJetCalo_itCone5,"pyJetCalo_itCone5[nJetCalo_itCone5]/F");
  m_tree->Branch("pzJetCalo_itCone5",pzJetCalo_itCone5,"pzJetCalo_itCone5[nJetCalo_itCone5]/F");
  m_tree->Branch("eJetCalo_itCone5",eJetCalo_itCone5,"eJetCalo_itCone5[nJetCalo_itCone5]/F");
  m_tree->Branch("etaJetCalo_itCone5",etaJetCalo_itCone5,"etaJetCalo_itCone5[nJetCalo_itCone5]/F");
  m_tree->Branch("phiJetCalo_itCone5",phiJetCalo_itCone5,"phiJetCalo_itCone5[nJetCalo_itCone5]/F");

/*  m_tree->Branch("nJetCalo_sisCone5",&nJetCalo_sisCone5,"nJetCalo_sisCone5/I");
  m_tree->Branch("pxJetCalo_sisCone5",pxJetCalo_sisCone5,"pxJetCalo_sisCone5[nJetCalo_sisCone5]/F");
  m_tree->Branch("pyJetCalo_sisCone5",pyJetCalo_sisCone5,"pyJetCalo_sisCone5[nJetCalo_sisCone5]/F");
  m_tree->Branch("pzJetCalo_sisCone5",pzJetCalo_sisCone5,"pzJetCalo_sisCone5[nJetCalo_sisCone5]/F");
  m_tree->Branch("eJetCalo_sisCone5",eJetCalo_sisCone5,"eJetCalo_sisCone5[nJetCalo_sisCone5]/F");
  m_tree->Branch("etaJetCalo_sisCone5",etaJetCalo_sisCone5,"etaJetCalo_sisCone5[nJetCalo_sisCone5]/F");
  m_tree->Branch("phiJetCalo_sisCone5",phiJetCalo_sisCone5,"phiJetCalo_sisCone5[nJetCalo_sisCone5]/F");

  m_tree->Branch("nJetCalo_kt4",&nJetCalo_kt4,"nJetCalo_kt4/I");
  m_tree->Branch("pxJetCalo_kt4",pxJetCalo_kt4,"pxJetCalo_kt4[nJetCalo_kt4]/F");
  m_tree->Branch("pyJetCalo_kt4",pyJetCalo_kt4,"pyJetCalo_kt4[nJetCalo_kt4]/F");
  m_tree->Branch("pzJetCalo_kt4",pzJetCalo_kt4,"pzJetCalo_kt4[nJetCalo_kt4]/F");
  m_tree->Branch("eJetCalo_kt4",eJetCalo_kt4,"eJetCalo_kt4[nJetCalo_kt4]/F");
  m_tree->Branch("etaJetCalo_kt4",etaJetCalo_kt4,"etaJetCalo_kt4[nJetCalo_kt4]/F");
  m_tree->Branch("phiJetCalo_kt4",phiJetCalo_kt4,"phiJetCalo_kt4[nJetCalo_kt4]/F");
*/
  m_tree->Branch("nJetPF_itCone5",&nJetPF_itCone5,"nJetPF_itCone5/I");
  m_tree->Branch("pxJetPF_itCone5",pxJetPF_itCone5,"pxJetPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyJetPF_itCone5",pyJetPF_itCone5,"pyJetPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzJetPF_itCone5",pzJetPF_itCone5,"pzJetPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("eJetPF_itCone5",eJetPF_itCone5,"eJetPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("etaJetPF_itCone5",etaJetPF_itCone5,"etaJetPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("phiJetPF_itCone5",phiJetPF_itCone5,"phiJetPF_itCone5[nJetPF_itCone5]/F");

  m_tree->Branch("nChargedHadronsPF_itCone5",nChargedHadronsPF_itCone5,"nChargedHadronsPF_itCone5[nJetPF_itCone5]/I");
  m_tree->Branch("nPhotonsPF_itCone5",       nPhotonsPF_itCone5,       "nPhotonsPF_itCone5[nJetPF_itCone5]/I");
  m_tree->Branch("nMuonsPF_itCone5",         nMuonsPF_itCone5,         "nMuonsPF_itCone5[nJetPF_itCone5]/I");
  m_tree->Branch("nElectronsPF_itCone5",     nElectronsPF_itCone5,     "nElectronsPF_itCone5[nJetPF_itCone5]/I");
  m_tree->Branch("nNeutralHadronsPF_itCone5",nNeutralHadronsPF_itCone5,"nNeutralHadronsPF_itCone5[nJetPF_itCone5]/I");
  m_tree->Branch("nHFHadronsPF_itCone5",     nHFHadronsPF_itCone5,     "nHFHadronsPF_itCone5[nJetPF_itCone5]/I");
  m_tree->Branch("nHFEMPF_itCone5",     nHFEMPF_itCone5,     "nHFEMPF_itCone5[nJetPF_itCone5]/I");

  m_tree->Branch("eChargedHadronsPF_itCone5",eChargedHadronsPF_itCone5,"eChargedHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("ePhotonsPF_itCone5",ePhotonsPF_itCone5,"ePhotonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("eMuonsPF_itCone5",eMuonsPF_itCone5,"eMuonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("eElectronsPF_itCone5",eElectronsPF_itCone5,"eElectronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("eNeutralHadronsPF_itCone5",eNeutralHadronsPF_itCone5,"eNeutralHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("eHFHadronsPF_itCone5",eHFHadronsPF_itCone5,"eHFHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("eHFEMPF_itCone5",eHFEMPF_itCone5,"eHFEMPF_itCone5[nJetPF_itCone5]/F");

  m_tree->Branch("pxChargedHadronsPF_itCone5",pxChargedHadronsPF_itCone5,"pxChargedHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pxPhotonsPF_itCone5",       pxPhotonsPF_itCone5,       "pxPhotonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pxMuonsPF_itCone5",         pxMuonsPF_itCone5,         "pxMuonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pxElectronsPF_itCone5",     pxElectronsPF_itCone5,     "pxElectronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pxNeutralHadronsPF_itCone5",pxNeutralHadronsPF_itCone5,"pxNeutralHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pxHFHadronsPF_itCone5",     pxHFHadronsPF_itCone5,     "pxHFHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pxHFEMPF_itCone5",     pxHFEMPF_itCone5,     "pxHFEMPF_itCone5[nJetPF_itCone5]/F");

  m_tree->Branch("pyChargedHadronsPF_itCone5",pyChargedHadronsPF_itCone5,"pyChargedHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyPhotonsPF_itCone5",       pyPhotonsPF_itCone5,       "pyPhotonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyMuonsPF_itCone5",         pyMuonsPF_itCone5,         "pyMuonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyElectronsPF_itCone5",     pyElectronsPF_itCone5,     "pyElectronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyNeutralHadronsPF_itCone5",pyNeutralHadronsPF_itCone5,"pyNeutralHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyHFHadronsPF_itCone5",     pyHFHadronsPF_itCone5,     "pyHFHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pyHFEMPF_itCone5",     pyHFEMPF_itCone5,     "pyHFEMPF_itCone5[nJetPF_itCone5]/F");

  m_tree->Branch("pzChargedHadronsPF_itCone5",pzChargedHadronsPF_itCone5,"pzChargedHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzPhotonsPF_itCone5",       pzPhotonsPF_itCone5,       "pzPhotonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzMuonsPF_itCone5",         pzMuonsPF_itCone5,         "pzMuonsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzElectronsPF_itCone5",     pzElectronsPF_itCone5,     "pzElectronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzNeutralHadronsPF_itCone5",pzNeutralHadronsPF_itCone5,"pzNeutralHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzHFHadronsPF_itCone5",     pzHFHadronsPF_itCone5,     "pzHFHadronsPF_itCone5[nJetPF_itCone5]/F");
  m_tree->Branch("pzHFEMPF_itCone5",     pzHFEMPF_itCone5,     "pzHFEMPF_itCone5[nJetPF_itCone5]/F");

  m_tree->Branch("averagePartPtPF_itCone5",averagePartPtPF_itCone5,"averagePartPtPF_itCone5[nJetPF_itCone5]/F");


  event = 0;  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetAnalyzer::endJob() {

}


//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
