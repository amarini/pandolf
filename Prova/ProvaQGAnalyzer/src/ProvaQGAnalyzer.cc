// -*- C++ -*-
//
// Package:    ProvaQGAnalyzer
// Class:      ProvaQGAnalyzer
// 
/**\class ProvaQGAnalyzer ProvaQGAnalyzer.cc MyAnalysis/ProvaQGAnalyzer/src/ProvaQGAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Daniele del Re
//         Created:  Thu Sep 13 16:00:15 CEST 2007
// $Id: ProvaQGAnalyzer.cc,v 1.1 2013/01/17 17:25:42 pandolf Exp $
//
//

//
// constructors and destructor
//

#include "Prova/ProvaQGAnalyzer/interface/ProvaQGAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "DataFormats/VertexReco/interface/VertexCollection.h"
//#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"
#include "TRegexp.h"
#include "TString.h"
#include "TVector3.h"


#include <set>
#include <algorithm>

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/PhotonTkIsolation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"




using namespace edm;
using namespace reco;


//values for Jet Vertex PU ID
#define DEF_GOODVTX_NDOF 4
#define DEF_GOODVTX_Z 24




ProvaQGAnalyzer::ProvaQGAnalyzer(const edm::ParameterSet& iConfig)
{
  outFileName= iConfig.getUntrackedParameter<std::string>("outFileName","output.root");
  JetPFsrcakt5_ = iConfig.getUntrackedParameter<edm::InputTag>("jetspfakt5");
  JetCorrector_pfakt5_ = iConfig.getParameter<std::string>("JetCorrectionService_pfakt5");

  //qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_8_patch7/src/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root");
  //qglikeli = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/CMSSW_5_3_7_patch4_testQG/src/Prova/ProvaQGAnalyzer/test");
  qglikeli = new QGLikelihoodCalculator();

}


ProvaQGAnalyzer::~ProvaQGAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   //TFile* file_prova = TFile::Open("prova.root", "recreate");
   //file_prova->cd();
   //h1_hbherh_detid->Write();
   //h1_etaPhot->Write();
   //h2_n_vs_eta->Write();
   //file_prova->Close();
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ProvaQGAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   nJet_pfakt5 = 0;



   isMC = !iEvent.isRealData(); // separate MC processing
   store = iEvent.eventAuxiliary().storeNumber(); // study stability across a store
   lbn = iEvent.luminosityBlock(); // sum LBN lumi for normalization
   bx = iEvent.bunchCrossing(); // study effect of out-of-time pile-up
   orbit = iEvent.orbitNumber(); // study beam heating with time (longer bunches)
   run = iEvent.id().run(); // unique ID - part 1
   event = iEvent.id().event(); // unique ID - part 2


   edm::Handle<double> rhoH;
   if( iEvent.getByLabel(edm::InputTag("kt6PFJetsForIso","rho"),rhoH) )
     rho = *rhoH;
   else 
     rho = 0;
        

   // get primary vertices
   Handle<VertexCollection> VertexHandle;
   //Handle<vector<Vertex> > VertexHandle;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS", VertexHandle);
   //iEvent.getByLabel(Vertexsrc_, VertexHandle);

   Handle<GenParticleCollection> genParticles;
   if( isMC ) iEvent.getByLabel("genParticles", genParticles);

   Handle<ValueMap<float> > qglMap;
   iEvent.getByLabel("qglAK5PF",qglMap);

   Handle<PFJetCollection> pfjetsakt5;
   iEvent.getByLabel(JetPFsrcakt5_, pfjetsakt5);

   //get jet correctors
   const JetCorrector* corrector_pfakt5 = 0;

   corrector_pfakt5 = JetCorrector::getJetCorrector (JetCorrector_pfakt5_, iSetup);

      for (PFJetCollection::const_iterator it = pfjetsakt5->begin(); 
	   it != pfjetsakt5->end(); ++it) {
	
	if (nJet_pfakt5>=100) {cout << "number of reco jets pfakt5 is larger than 100. Skipping" << endl; continue;}
	  
	  ptJet_pfakt5[nJet_pfakt5] = it->pt();
	  eJet_pfakt5[nJet_pfakt5] = it->energy();	 
	  etaJet_pfakt5[nJet_pfakt5] = it->eta();	 
	  phiJet_pfakt5[nJet_pfakt5] = it->phi();	      

	  // Jet Energy Scale Corrections on-the-fly     
	  PFJet  correctedJet = *it;
	  edm::RefToBase<reco::Jet> jetRef(edm::Ref<PFJetCollection>(pfjetsakt5,nJet_pfakt5));
	  double scale = corrector_pfakt5->correction(*it,iEvent,iSetup);
	  correctedJet.scaleEnergy(scale);
	  ptCorrJet_pfakt5[nJet_pfakt5] = correctedJet.pt();

        if( ptCorrJet_pfakt5[nJet_pfakt5] < 20. ) continue;


	  
	  // Extra variables for PFlow studies
	  Int_t nChargedHadrons = 0;
	  Int_t nPhotons = 0;
	  Int_t nNeutralHadrons = 0;
	  Int_t nElectrons = 0;
	  Int_t nMuons = 0;
	  Int_t nHFHadrons = 0;
	  Int_t nHFEM = 0;
	  
	  TLorentzVector p4ChargedHadrons;
	  TLorentzVector p4Photons;
	  TLorentzVector p4NeutralHadrons;
	  TLorentzVector p4Electrons;
	  TLorentzVector p4Muons;
	  TLorentzVector p4HFHadrons;
	  TLorentzVector p4HFEM;
	  
	  vector<PFCandidatePtr> pfCandidates = it->getPFConstituents();
	  
	  float rms_cands_wrong=0.;

	  float SumW=0;
	  float SumW2=0;
	  float SumDeta=0;
	  float SumDeta2=0;
	  float SumDphi=0;
	  float SumDphi2=0;
	  float SumDetaDphi=0;
	  
	  float SumW_QC=0;
	  float SumW2_QC=0;
	  float SumDeta_QC=0;
	  float SumDeta2_QC=0;
	  float SumDphi_QC=0;
	  float SumDphi2_QC=0;
	  float SumDetaDphi_QC=0;
	  
	  float Eta0=it->eta();
	  float Phi0=it->phi();

        // initialize:
        rmsCandJet_pfakt5[nJet_pfakt5] =  -999.;
        ptDJet_pfakt5[nJet_pfakt5] =      -999.;
        axis1Jet_pfakt5[nJet_pfakt5] =    -999.;
        axis2Jet_pfakt5[nJet_pfakt5] =    -999.;
        pullJet_pfakt5[nJet_pfakt5] =     -999.;
        tanaJet_pfakt5[nJet_pfakt5]  =    -999.;

        rmsCandTrue_QCJet_pfakt5[nJet_pfakt5] =  -999.;
        ptD_QCJet_pfakt5[nJet_pfakt5] =      -999.;
        axis1_QCJet_pfakt5[nJet_pfakt5] =    -999.;
        axis2_QCJet_pfakt5[nJet_pfakt5] =    -999.;
        pull_QCJet_pfakt5[nJet_pfakt5] =     -999.;
        tana_QCJet_pfakt5[nJet_pfakt5]  =    -999.;

        RchgJet_pfakt5[nJet_pfakt5] = 0.;
        RneutralJet_pfakt5[nJet_pfakt5] = 0.;
        RJet_pfakt5[nJet_pfakt5] = 0.;
        Rchg_QCJet_pfakt5[nJet_pfakt5] = 0.;

	  pTMaxJet_pfakt5[nJet_pfakt5] = 0.;
	  pTMaxChgJet_pfakt5[nJet_pfakt5] = 0.;
	  pTMaxNeutralJet_pfakt5[nJet_pfakt5] = 0.;
        pTMaxChg_QCJet_pfakt5[nJet_pfakt5] = 0.;

        qglJet_pfakt5[nJet_pfakt5] = -1.;

        pdgIdPartJet_pfakt5[nJet_pfakt5] = 0;

        nChg_ptCutJet_pfakt5[nJet_pfakt5] = 0;
        nChg_QCJet_pfakt5[nJet_pfakt5] = 0;
        nChg_ptCut_QCJet_pfakt5[nJet_pfakt5] = 0;
        nNeutral_ptCutJet_pfakt5[nJet_pfakt5] = 0;

        nPFCand_QC_ptCutJet_pfakt5[nJet_pfakt5] = 0;

        std::vector<bool> jetPart_forMult,jetPart_forAxis;

	  
	  for(int i=0;i<it->nConstituents();++i)
		{

		reco::TrackRef itrk ;
		reco::PFCandidatePtr  part = it->getPFConstituent(i);

		double pt=part->pt();
		double eta=part->eta();
		double phi=part->phi();

		PFCandidate::ParticleType id = part->particleId();
		// Convert particle momentum to normal TLorentzVector, wrong type :(
		math::XYZTLorentzVectorD const& p4t = part->p4();
		TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
		TLorentzVector jetp4;
		jetp4.SetPtEtaPhiE(it->pt(), it->eta(), it->phi(), it->energy());

		if (id==PFCandidate::h) { // charged hadrons
		  nChargedHadrons += 1;
		  p4ChargedHadrons += p4;
		}
		if (id==PFCandidate::e) { // electrons
		  nElectrons += 1;
		  p4Electrons += p4;
		}
		if (id==PFCandidate::mu) { // muons
		  nMuons += 1;
		  p4Muons += p4;
		}
		if (id==PFCandidate::gamma) { // photons
		  nPhotons += 1;
		  p4Photons += p4;
		}
		if (id==PFCandidate::h0) { // neutral hadrons
		  nNeutralHadrons += 1;
		  p4NeutralHadrons += p4;
		}
		if (id==PFCandidate::h_HF) { // HF hadrons
		  nHFHadrons += 1;
		  p4HFHadrons += p4;
		}
		if (id==PFCandidate::egamma_HF) { // HF EM clusters
		  nHFEM += 1;
		  p4HFEM += p4;
		}



		if (part.isNonnull())
		  itrk = (*part).trackRef();
		if (pt > pTMaxJet_pfakt5[nJet_pfakt5]) 
		  pTMaxJet_pfakt5[nJet_pfakt5] = pt;
		if (itrk.isNonnull() && pt > pTMaxChgJet_pfakt5[nJet_pfakt5]) 
		  pTMaxChgJet_pfakt5[nJet_pfakt5] = pt;
		if (!itrk.isNonnull() && pt > pTMaxNeutralJet_pfakt5[nJet_pfakt5]) 
		  pTMaxNeutralJet_pfakt5[nJet_pfakt5] = pt;
		if (!itrk.isNonnull() && pt > 1.0) 
		  nNeutral_ptCutJet_pfakt5[nJet_pfakt5]++;
		
		bool trkForAxis = false;
		bool trkForMult = false;
		
		//-----matching with vertex tracks-------
		if (!itrk.isNonnull()) { 
		  trkForMult = true;
		  trkForAxis = true;
		}
		else {
		  if (pt > 1.0)
		    nChg_ptCutJet_pfakt5[nJet_pfakt5]++;
		  float dZmin = 999;
		  int index_min = 999;
		  reco::VertexCollection::const_iterator vtxClose;
		  for(unsigned ivtx = 0;ivtx < VertexHandle->size();ivtx++) {
		    float dZ_cut = fabs(itrk->dz((*VertexHandle)[ivtx].position()));
		    float sumpT = 0;
		    for(reco::Vertex::trackRef_iterator itk = (*VertexHandle)[ivtx].tracks_begin();itk!=(*VertexHandle)[ivtx].tracks_end(); ++itk) {
		      sumpT = sumpT + ((*itk)->pt())*((*itk)->pt());
		    }
		    if (dZ_cut < dZmin) {
		      dZmin = dZ_cut;
		      index_min = ivtx;
		        //  std::cout<<"dz=="<<dZ_cut<<std::endl;
		    }
		  }//Loop over vertices 
		  if (index_min == 0) {
		    float dz = itrk->dz((*VertexHandle)[0].position());
		    float d0 = itrk->dxy((*VertexHandle)[0].position());
		    float vtx_xError = (*VertexHandle)[0].xError();
		    float vtx_yError = (*VertexHandle)[0].yError();
		    float vtx_zError = (*VertexHandle)[0].zError();
		    float d0_sigma=sqrt(pow(itrk->d0Error(),2) + pow(vtx_xError,2) + pow(vtx_yError,2));
		    float dz_sigma=sqrt(pow(itrk->dzError(),2) + pow(vtx_zError,2));
		    if (itrk->quality(reco::TrackBase::qualityByName("highPurity")) && fabs(dz/dz_sigma) < 5.) {
		      trkForAxis = true;
		      if (fabs(d0/d0_sigma) < 5.)
		        trkForMult = true;
		    }//
		  }
		  if (trkForMult)
		    nChg_QCJet_pfakt5[nJet_pfakt5]++;
		  if (itrk.isNonnull() && trkForMult && pt > 1.0)
		    nChg_ptCut_QCJet_pfakt5[nJet_pfakt5]++;
		  if (pt > pTMaxChg_QCJet_pfakt5[nJet_pfakt5] && trkForAxis) 
		    pTMaxChg_QCJet_pfakt5[nJet_pfakt5] = pt;
		}// for charged particles only


            jetPart_forMult.push_back(trkForMult);
            jetPart_forAxis.push_back(trkForAxis);

		double dphi = 2*atan(tan((phi-Phi0)/2));      
		double deta = eta-Eta0;
		SumW+=pt;
		SumW2+=pt*pt;
		SumDeta+=pt*pt*deta;
		SumDeta2+=pt*pt*deta*deta;
		SumDphi+=pt*pt*dphi;
		SumDphi2+=pt*pt*dphi*dphi;
		SumDetaDphi+=pt*pt*deta*dphi;
		float deltaR = jetp4.DeltaR(p4);
		rms_cands_wrong += (p4.Pt()*p4.Pt()*deltaR*deltaR);
		if (trkForAxis) {
		  SumW_QC+=pt;
		  SumW2_QC+=pt*pt;
		  SumDeta_QC+=pt*pt*deta;
		  SumDeta2_QC+=pt*pt*deta*deta;
		  SumDphi_QC+=pt*pt*dphi;
		  SumDphi2_QC+=pt*pt*dphi*dphi;
		  SumDetaDphi_QC+=pt*pt*deta*dphi;
		}

	}
	float ave_deta = SumDeta/SumW2;
	float ave_dphi = SumDphi/SumW2;
	float  ave_deta2 = SumDeta2/SumW2;
	float  ave_dphi2 = SumDphi2/SumW2;
      float a = ave_deta2-ave_deta*ave_deta;
      float b = ave_dphi2-ave_dphi*ave_dphi;
      float c = -(SumDetaDphi/SumW2-ave_deta*ave_dphi);
      float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
      if (a+b+delta > 0) {
        axis1Jet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a+b+delta));
      }
      if (a+b-delta > 0) {  
        axis2Jet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a+b-delta));
      }
      if (c != 0) {
        tanaJet_pfakt5[nJet_pfakt5] = 0.5*(b-a+delta)/c;
      }	

	float ave_deta_QC = SumDeta_QC/SumW2_QC;
	float ave_dphi_QC = SumDphi_QC/SumW2_QC;
	float  ave_deta2_QC = SumDeta2_QC/SumW2_QC;
	float  ave_dphi2_QC = SumDphi2_QC/SumW2_QC;
      float a_QC = ave_deta2_QC-ave_deta_QC*ave_deta_QC;
      float b_QC = ave_dphi2_QC-ave_dphi_QC*ave_dphi_QC;
      float c_QC = -(SumDetaDphi_QC/SumW2_QC-ave_deta_QC*ave_dphi_QC);
      float delta_QC = sqrt(fabs((a_QC-b_QC)*(a_QC-b_QC)+4*c_QC*c_QC));
      if (a_QC+b_QC+delta_QC > 0) {
        axis1_QCJet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a_QC+b_QC+delta_QC));
      }
      if (a_QC+b_QC-delta_QC > 0) {  
        axis2_QCJet_pfakt5[nJet_pfakt5] = sqrt(0.5*(a_QC+b_QC-delta_QC));
      }
      if (c_QC != 0) {
        tana_QCJet_pfakt5[nJet_pfakt5] = 0.5*(b_QC-a_QC+delta_QC)/c_QC;
      }	

      ptDJet_pfakt5[nJet_pfakt5] =sqrt( SumW2/ (SumW*SumW));
      ptD_QCJet_pfakt5[nJet_pfakt5] =sqrt( SumW2_QC/ (SumW_QC*SumW_QC));

      // this is thw rong definition of rms, kept only for backwards compatibility:
	rmsCandJet_pfakt5[nJet_pfakt5] = rms_cands_wrong/SumW2;

	//-------calculate pull------
    	float ddetaR_sum(0.0), ddphiR_sum(0.0),ddetaR_sum_QC(0.0), ddphiR_sum_QC(0.0);
      float sum_ddR = 0.;
      float sum_ddR_QC = 0.;
    	for(int i=0; i<it->nConstituents(); ++i) {
			double pt=it->getJetConstituentsQuick()[i]->pt();
			double eta=it->getJetConstituentsQuick()[i]->eta();
			double phi=it->getJetConstituentsQuick()[i]->phi();
			double dphi = 2*atan(tan((phi-Phi0)/2));      
			double deta = eta-Eta0;
  		    float weight = pt*pt;
  		    float ddeta, ddphi,ddR;
  		    ddeta = deta - ave_deta ;//jetPart_deta[i] - ave_deta ; 
  		    ddphi = 2*atan(tan(( dphi - ave_dphi)/2.)) ;
  		    ddR = sqrt(ddeta*ddeta + ddphi*ddphi);
		    sum_ddR += ddR *ddR* weight;
  		    ddetaR_sum += ddR*ddeta*weight;
  		    ddphiR_sum += ddR*ddphi*weight;
                if (jetPart_forAxis[i]) { // this should be ave_deta_QC
  		      float ddeta_QC = deta - ave_deta_QC ;//jetPart_deta[i] - ave_deta ; 
  		      float ddphi_QC = 2*atan(tan(( dphi - ave_dphi_QC)/2.)) ;
  		      float ddR_QC = sqrt(ddeta_QC*ddeta_QC + ddphi_QC*ddphi_QC);
		      sum_ddR_QC += ddR_QC *ddR_QC* weight;
  		      ddetaR_sum_QC += ddR_QC*ddeta_QC*weight;
  		      ddphiR_sum_QC += ddR_QC*ddphi_QC*weight;
                }
  		  }//second loop over constituents  
  if (SumW2 > 0) {
    float ddetaR_ave = ddetaR_sum/SumW2;
    float ddphiR_ave = ddphiR_sum/SumW2;
    pullJet_pfakt5[nJet_pfakt5] = sqrt(ddetaR_ave*ddetaR_ave+ddphiR_ave*ddphiR_ave);
  }

  if (SumW2_QC > 0) {
    float ddetaR_ave_QC = ddetaR_sum_QC/SumW2_QC;
    float ddphiR_ave_QC = ddphiR_sum_QC/SumW2_QC;
    pull_QCJet_pfakt5[nJet_pfakt5] = sqrt(ddetaR_ave_QC*ddetaR_ave_QC+ddphiR_ave_QC*ddphiR_ave_QC);
  }

  nPFCand_QC_ptCutJet_pfakt5[nJet_pfakt5] = nChg_QCJet_pfakt5[nJet_pfakt5] + nNeutral_ptCutJet_pfakt5[nJet_pfakt5];

  rmsCandTrueJet_pfakt5[nJet_pfakt5] = sqrt( sum_ddR / SumW2);
  rmsCandTrue_QCJet_pfakt5[nJet_pfakt5] = sqrt( sum_ddR_QC / SumW2_QC);

  RchgJet_pfakt5[nJet_pfakt5] = pTMaxChgJet_pfakt5[nJet_pfakt5]/SumW;
  RneutralJet_pfakt5[nJet_pfakt5] = pTMaxNeutralJet_pfakt5[nJet_pfakt5]/SumW;
  RJet_pfakt5[nJet_pfakt5] = pTMaxJet_pfakt5[nJet_pfakt5]/SumW;
  Rchg_QCJet_pfakt5[nJet_pfakt5] = pTMaxChg_QCJet_pfakt5[nJet_pfakt5]/SumW_QC;


	    //float qglNEW = qglikeli->computeQGLikelihood( ptCorrJet_pfakt5[nJet_pfakt5], etaJet_pfakt5[nJet_pfakt5], rho, nNeutral_ptCutJet_pfakt5[nJet_pfakt5]+nChg_QCJet_pfakt5[nJet_pfakt5], ptD_QCJet_pfakt5[nJet_pfakt5], axis2_QCJet_pfakt5[nJet_pfakt5] ); 
          qglJet_pfakt5[nJet_pfakt5] = (*qglMap)[jetRef];
          //std::cout << ptCorrJet_pfakt5[nJet_pfakt5] << " " << etaJet_pfakt5[nJet_pfakt5] << " " << rho << " " << nNeutral_ptCutJet_pfakt5[nJet_pfakt5]+nChg_QCJet_pfakt5[nJet_pfakt5] << " "<< ptD_QCJet_pfakt5[nJet_pfakt5] << " " << axis2_QCJet_pfakt5[nJet_pfakt5] << std::endl;
          //std::cout << qgl << " " << qglNEW << std::endl;


//	  float sumPt_cands=0.;
//	  float sumPt2_cands=0.;
//	  
//	  for (vector<PFCandidatePtr>::const_iterator jt = pfCandidates.begin();
//	       jt != pfCandidates.end(); ++jt) {
//	    
//	    PFCandidate::ParticleType id = (*jt)->particleId();
//	    // Convert particle momentum to normal TLorentzVector, wrong type :(
//	    math::XYZTLorentzVectorD const& p4t = (*jt)->p4();
//	    TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
//	    TLorentzVector jetp4;
//	    jetp4.SetPtEtaPhiE(it->pt(), it->eta(), it->phi(), it->energy());
//	    if(p4.Pt()!=0){
//	      sumPt_cands += p4.Pt();
//	      sumPt2_cands += (p4.Pt()*p4.Pt());
//	      //float deltaR = it->p4().DeltaR(p4);
//	      float deltaR = jetp4.DeltaR(p4);
//	      rms_cands += (p4.Pt()*p4.Pt()*deltaR*deltaR);
//	    }

	    
	  //} //for PFCandidates

	  //ptDJet_pfakt5[nJet_pfakt5] = sqrt( sumPt2_cands )/sumPt_cands;
	  //rmsCandJet_pfakt5[nJet_pfakt5] = rms_cands/sumPt2_cands;


	 // const TLorentzVector *p = 0;
	  
	 // nChargedHadrons_pfakt5[nJet_pfakt5] =  nChargedHadrons;
	 // p = &p4ChargedHadrons;
	 // eChargedHadrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	 // 
	 // nElectrons_pfakt5[nJet_pfakt5] =  nElectrons;
	 // p = &p4Electrons;
	 // eElectrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	 // 
	 // nMuons_pfakt5[nJet_pfakt5] =  nMuons;
	 // p = &p4Muons;
	 // eMuons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	 // 
	 // nPhotons_pfakt5[nJet_pfakt5] =  nPhotons;
	 // p = &p4Photons;
	 // ePhotons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	 // 
	 // nNeutralHadrons_pfakt5[nJet_pfakt5] =  nNeutralHadrons;
	 // p = &p4NeutralHadrons;
	 // eNeutralHadrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	 // 
	 // nHFHadrons_pfakt5[nJet_pfakt5] =  nHFHadrons;
	 // p = &p4HFHadrons;
	 // eHFHadrons_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	 // 
	 // nHFEM_pfakt5[nJet_pfakt5] =  nHFEM;
	 // p = &p4HFEM;
	 // eHFEM_pfakt5[nJet_pfakt5] = p->E() / it->energy();
	  
	  
	  int index = nJet_pfakt5;
	  combinedSecondaryVertexBJetTags_pfakt5[index] = -999.;
	  combinedSecondaryVertexMVABJetTags_pfakt5[index] = -999.;
	  jetBProbabilityBJetTags_pfakt5[index] = -999.;
	  jetProbabilityBJetTags_pfakt5[index] =  -999.;
	  simpleSecondaryVertexHighEffBJetTags_pfakt5[index] =  -999.;
	  simpleSecondaryVertexHighPurBJetTags_pfakt5[index] =  -999.;
	  softMuonBJetTags_pfakt5[index] =  -999.;
	  softMuonByIP3dBJetTags_pfakt5[index] =  -999.;
	  softMuonByPtBJetTags_pfakt5[index] =  -999.;
	  //softElectronBJetTags_pfakt5[index] =  (*softElectronBJetTags)[index].second ;
	  softElectronByIP3dBJetTags_pfakt5[index] =   -999.;
	  softElectronByPtBJetTags_pfakt5[index]     = -999.;    
	  trackCountingHighPurBJetTags_pfakt5[index] = -999.;
	  trackCountingHighEffBJetTags_pfakt5[index] = -999.;
	  
	  //PU id
	  reco::TrackRefVector vTrks(it->getTrackRefs());
	  float sumTrkPt(0.0);
	  float sumTrkPtBeta[100],sumTrkPtBetaStar[100];
	  for (int ivtx=0;ivtx<100;++ivtx)
	    {
	      sumTrkPtBeta[ivtx]=0.;
	      sumTrkPtBetaStar[ivtx]=0.;
	    }

	  //---- loop over the tracks of the jet ----
	  for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
	    if ( VertexHandle->size() == 0) break;
	    sumTrkPt += (*i_trk)->pt();
	    int myVertex=-1;
	    //---- loop over all vertices ----------------------------
	    for(unsigned ivtx = 0;ivtx <  VertexHandle->size();ivtx++) {
	      //---- loop over the tracks associated with the vertex ---
	      if (!((* VertexHandle)[ivtx].isFake()) && (* VertexHandle)[ivtx].ndof() >= DEF_GOODVTX_NDOF && fabs((* VertexHandle)[ivtx].z()) <= DEF_GOODVTX_Z) {
		for(reco::Vertex::trackRef_iterator i_vtxTrk = (* VertexHandle)[ivtx].tracks_begin(); i_vtxTrk != (* VertexHandle)[ivtx].tracks_end(); ++i_vtxTrk) {
		  //---- match the jet track to the track from the vertex ----
		  reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
		  //---- check if the tracks match -------------------------
		  if (trkRef == (*i_trk)) 
		    {
		      myVertex=ivtx;
		      break;
		    }
		}
	      }
	    }
	    if (myVertex>-1)
	      {
		for(unsigned ivtx = 0;ivtx <  VertexHandle->size();ivtx++) {
		  if (ivtx== (unsigned int) myVertex)
		    sumTrkPtBeta[ivtx] += (*i_trk)->pt();
		  else
		    sumTrkPtBetaStar[ivtx] += (*i_trk)->pt();
		}
	      }
	  }
	  
	  
	  if (sumTrkPt > 0) {
	    for(unsigned ivtx = 0; ivtx <  VertexHandle->size();ivtx++) {
	      beta_pfakt5[index][ivtx]     = sumTrkPtBeta[ivtx]/sumTrkPt;
	      betaStar_pfakt5[index][ivtx] = sumTrkPtBetaStar[ivtx]/sumTrkPt;
	    }
	  }

	//if(it->pt() > 10. ){
	//  // B tagging
	//  combinedSecondaryVertexBJetTags_pfakt5[index] =  (*combinedSecondaryVertexBJetTags)[index].second ;
	//  combinedSecondaryVertexMVABJetTags_pfakt5[index] =  (*combinedSecondaryVertexMVABJetTags)[index].second ;
	//  jetBProbabilityBJetTags_pfakt5[index] =  (*jetBProbabilityBJetTags)[index].second ;
	//  jetProbabilityBJetTags_pfakt5[index] =  (*jetProbabilityBJetTags)[index].second ;
	//}
	

     TLorentzVector thisJet;
     thisJet.SetPtEtaPhiE( ptCorrJet_pfakt5[index], etaJet_pfakt5[index], phiJet_pfakt5[index], eJet_pfakt5[index] );

     float bestDeltaR = 999.;
     int pdgId_found = 0;

     int nMC=0;

     for (GenParticleCollection::const_iterator p = genParticles->begin(); p != genParticles->end(); ++p) {


       if (p->status() != 3) continue;
       if( !(abs(p->pdgId())<6 || p->pdgId()==21) ) continue;
       if( p->pt()<1. ) continue;

       nMC++;
       if( nMC > 100 ) break;

       TLorentzVector thisParton;
       thisParton.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );

       float thisDeltaR = thisJet.DeltaR(thisParton);

       if( thisDeltaR < bestDeltaR ) {
         bestDeltaR = thisDeltaR;
         pdgId_found = p->pdgId();
       }

     } // for partons


      pdgIdPartJet_pfakt5[index] = (bestDeltaR<0.3) ? pdgId_found : 0;

	  ++nJet_pfakt5;
	  
      } // pfakt5
 
  

   m_tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void 
ProvaQGAnalyzer::beginJob()
{

  //m_tree = fs_->make<TTree>("pippo","Analysis tree");
  outfile = TFile::Open(outFileName.c_str(), "RECREATE");
  outfile->mkdir("myanalysis");
  outfile->cd("myanalysis");

  m_tree = new TTree ("pippo","Analysis tree") ;
  //  m_tree->SetAutoSave (10000000) ;

  m_tree->Branch("isMC",&isMC,"isMC/O");
  m_tree->Branch("store",&store,"store/I");
  m_tree->Branch("lbn",&lbn,"lbn/I");
  m_tree->Branch("bx",&bx,"bx/I");
  m_tree->Branch("orbit",&orbit,"orbit/I");
  m_tree->Branch("run",&run,"run/I");
  m_tree->Branch("event",&event,"event/I");

  m_tree->Branch("rhoPF",&rho,"rhoPF/F");

      m_tree->Branch("nJet_pfakt5",&nJet_pfakt5,"nJet_pfakt5/I");
      m_tree->Branch("ptJet_pfakt5 ",ptJet_pfakt5 ,"ptJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptCorrJet_pfakt5 ",ptCorrJet_pfakt5 ,"ptCorrJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eJet_pfakt5  ",eJet_pfakt5  ,"eJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("etaJet_pfakt5",etaJet_pfakt5,"etaJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("phiJet_pfakt5",phiJet_pfakt5,"phiJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ptDJet_pfakt5",ptDJet_pfakt5,"ptDJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("rmsCandJet_pfakt5",rmsCandJet_pfakt5,"rmsCandJet_pfakt5[nJet_pfakt5]/F");

      m_tree->Branch("axis2_QCJet_pfakt5",axis2_QCJet_pfakt5,"axis2_QCJet_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("nPFCand_QC_ptCutJet_pfakt5",nPFCand_QC_ptCutJet_pfakt5,"nPFCand_QC_ptCutJet_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("ptD_QCJet_pfakt5",ptD_QCJet_pfakt5,"ptD_QCJet_pfakt5[nJet_pfakt5]/F");

      m_tree->Branch("qglJet_pfakt5",qglJet_pfakt5,"qglJet_pfakt5[nJet_pfakt5]/F");

      m_tree->Branch("pdgIdPartJet_pfakt5",pdgIdPartJet_pfakt5,"pdgIdPartJet_pfakt5[nJet_pfakt5]/I");

      

      //   // Extra variables for PFlow studies
      m_tree->Branch("nChargedHadrons_pfakt5",nChargedHadrons_pfakt5,"nChargedHadrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nPhotons_pfakt5",       nPhotons_pfakt5,       "nPhotons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nMuons_pfakt5",         nMuons_pfakt5,         "nMuons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nElectrons_pfakt5",     nElectrons_pfakt5,     "nElectrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nNeutralHadrons_pfakt5",nNeutralHadrons_pfakt5,"nNeutralHadrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nHFHadrons_pfakt5",     nHFHadrons_pfakt5,     "nHFHadrons_pfakt5[nJet_pfakt5]/I");
      m_tree->Branch("nHFEM_pfakt5",     nHFEM_pfakt5,     "nHFEM_pfakt5[nJet_pfakt5]/I");
      
      m_tree->Branch("eChargedHadrons_pfakt5",eChargedHadrons_pfakt5,"eChargedHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("ePhotons_pfakt5",ePhotons_pfakt5,"ePhotons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eMuons_pfakt5",eMuons_pfakt5,"eMuons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eElectrons_pfakt5",eElectrons_pfakt5,"eElectrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eNeutralHadrons_pfakt5",eNeutralHadrons_pfakt5,"eNeutralHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eHFHadrons_pfakt5",eHFHadrons_pfakt5,"eHFHadrons_pfakt5[nJet_pfakt5]/F");
      m_tree->Branch("eHFEM_pfakt5",eHFEM_pfakt5,"eHFEM_pfakt5[nJet_pfakt5]/F");

  m_tree->Branch("Xsec",  &Xsec_, "Xsec/D");


  event = 0;  
}

// ------------ method called once each job just after ending the event loop  ------------
void ProvaQGAnalyzer::endJob() {

  
std::cout << "in endjob" << std::endl;
  outfile->cd("myanalysis");
  m_tree->Write();
  outfile->Close();

std::cout << "finish endjob" << std::endl;
  //outfile->Delete();
  
//   //avoid writing the tree second time (automatically)
//  m_tree->Delete();

}



//define this as a plug-in
DEFINE_FWK_MODULE(ProvaQGAnalyzer);
