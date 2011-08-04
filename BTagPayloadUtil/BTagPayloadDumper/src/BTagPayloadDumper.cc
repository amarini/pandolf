// -*- C++ -*-
//
// Package:    BTagPayloadDumper
// Class:      BTagPayloadDumper
// 
/**\class BTagPayloadDumper BTagPayloadDumper.cc UserCode/BTagPayloadDumper/src/BTagPayloadDumper.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Pandolfi,32 4-C03,+41227672087,
//         Created:  Wed Aug  3 12:20:40 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoBTag/Records/interface/BTagPerformanceRecord.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TFile.h"
#include "TH2F.h"


//
// class declaration
//

class BTagPayloadDumper : public edm::EDAnalyzer {
   public:
      explicit BTagPayloadDumper(const edm::ParameterSet&);
      ~BTagPayloadDumper();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
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
BTagPayloadDumper::BTagPayloadDumper(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


BTagPayloadDumper::~BTagPayloadDumper()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BTagPayloadDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  //using namespace reco;


 ///////////////////////////////
 ///   Begin DB setup
 ///////////////////////////////

            //// This is needed for the DB
  std::map<std::string,PerformanceResult::ResultType> measureMap;
  measureMap["BTAGBEFF"]=PerformanceResult::BTAGBEFF;
  measureMap["BTAGBERR"]=PerformanceResult::BTAGBERR;
  measureMap["BTAGCEFF"]=PerformanceResult::BTAGCEFF;
  measureMap["BTAGCERR"]=PerformanceResult::BTAGCERR;
  measureMap["BTAGLEFF"]=PerformanceResult::BTAGLEFF;
  measureMap["BTAGLERR"]=PerformanceResult::BTAGLERR;
  measureMap["BTAGNBEFF"]=PerformanceResult::BTAGNBEFF;
  measureMap["BTAGNBERR"]=PerformanceResult::BTAGNBERR;
  measureMap["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR;
  measureMap["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  measureMap["BTAGCEFFCORR"]=PerformanceResult::BTAGCEFFCORR;
  measureMap["BTAGCERRCORR"]=PerformanceResult::BTAGCERRCORR;
  measureMap["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  measureMap["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;
  measureMap["BTAGNBEFFCORR"]=PerformanceResult::BTAGNBEFFCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["BTAGNBERRCORR"]=PerformanceResult::BTAGNBERRCORR;
  measureMap["MUEFF"]=PerformanceResult::MUEFF;
  measureMap["MUERR"]=PerformanceResult::MUERR;
  measureMap["MUFAKE"]=PerformanceResult::MUFAKE; 
  measureMap["MUEFAKE"]=PerformanceResult::MUEFAKE;

  edm::ESHandle<BtagPerformance> perfH_BTAG;
  edm::ESHandle<BtagPerformance> perfH_MISTAG;

  std::vector<std::string> measureName;
  std::vector<std::string> measureType;
  
  // Define which Btag and Mistag algorithm you want to use. These are not user defined and need to be exact
  measureName.push_back("MISTAGTCHEM");
  measureName.push_back("BTAGTCHEM");
  measureName.push_back("MISTAGTCHEL");
  measureName.push_back("BTAGTCHEL");
  measureName.push_back("MISTAGTCHEM");
  measureName.push_back("MISTAGTCHEL");

  // Tell DB you want the SF. These are not user defined and need to be exact
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFFCORR");
  measureType.push_back("BTAGBEFFCORR");
  measureType.push_back("BTAGLEFF");
  measureType.push_back("BTAGLEFF");

  // These are user defined maps that we will use to store the SF 
  std::map<std::string,float> ScaleFactors_j0;   //store the Btag and Mistag SF for jet0
  std::map<std::string,float> ScaleFactors_j1;   //store the Btag and Mistag SF for jet1
  std::map<std::string,float> ScaleFactorsEff_j0;   //store the Mistag eff for jet0
  std::map<std::string,float> ScaleFactorsEff_j1;   //store the Mistag eff for jet1


            ///////////////////////////////
            ///   End DB setup
            ///////////////////////////////


  float ptMin = 20.;
  float ptMax = 500.; 
  float ptBinWidth = 1.; 
  unsigned int nBinsPt = (int)((ptMax-ptMin)/ptBinWidth);

  float etaMax = 2.4; 
  float etaBinWidth = 0.1; 
  unsigned int nBinsEta = (int)(etaMax/etaBinWidth);


  std::vector<std::string> btagAlgos;
  btagAlgos.push_back("TCHEL");
  btagAlgos.push_back("TCHEM");



  //for( size_t iMeasure = 0; iMeasure < measureName.size(); iMeasure++ ) {
  for( unsigned int iBtagAlgo=0; iBtagAlgo<btagAlgos.size(); ++iBtagAlgo ) {

    std::string outFileName = "BTagPayloads_"+btagAlgos[iBtagAlgo]+".root";

    TFile* outfile = TFile::Open(outFileName.c_str(), "RECREATE");
    outfile->cd();


    std::string labelName_btag = "BTAG"+btagAlgos[iBtagAlgo];
    std::string labelName_mistag = "MISTAG"+btagAlgos[iBtagAlgo];
    
    TH2F* h2_BTAGLEFFCORR = new TH2F("BTAGLEFFCORR", "", nBinsPt, ptMin, ptMax, nBinsEta, 0., etaMax);
    TH2F* h2_BTAGBEFFCORR = new TH2F("BTAGBEFFCORR", "", nBinsPt, ptMin, ptMax, nBinsEta, 0., etaMax);
    TH2F* h2_BTAGLEFF = new TH2F("BTAGLEFF", "", nBinsPt, ptMin, ptMax, nBinsEta, 0., etaMax);


     for( unsigned int iEta=0; iEta<nBinsEta; ++iEta ) {

       float thisEta = etaBinWidth*iEta;

       for( unsigned int iPt=0; iPt<nBinsPt; ++iPt ) {
      
         float thisPt = ptBinWidth*iPt + ptMin;
      
         //Setup our measurement
         iSetup.get<BTagPerformanceRecord>().get( labelName_btag.c_str(), perfH_BTAG);
         const BtagPerformance & perf_BTAG = *(perfH_BTAG.product());
         iSetup.get<BTagPerformanceRecord>().get( labelName_mistag.c_str(), perfH_MISTAG);
         const BtagPerformance & perf_MISTAG = *(perfH_MISTAG.product());

         BinningPointByMap measurePoint;
         measurePoint.reset();
         measurePoint.insert(BinningVariables::JetEt, thisPt );                         ///// pass in the et of the jet
         measurePoint.insert(BinningVariables::JetAbsEta, abs(  thisEta ) );       ///// pass in the absolute eta of the jet

         h2_BTAGLEFFCORR->SetBinContent( iPt, iEta, perf_MISTAG.getResult( measureMap[ "BTAGLEFFCORR" ], measurePoint) );
         h2_BTAGBEFFCORR->SetBinContent( iPt, iEta, perf_BTAG.getResult( measureMap[ "BTAGBEFFCORR" ], measurePoint) );
         h2_BTAGLEFF->SetBinContent( iPt, iEta, perf_MISTAG.getResult( measureMap[ "BTAGLEFF" ], measurePoint) );
         
       
       } //for eta

     } //for pt

     h2_BTAGLEFFCORR->Write();
     h2_BTAGBEFFCORR->Write();
     h2_BTAGLEFF->Write();
     
     outfile->Close();

     delete h2_BTAGLEFFCORR;
     delete h2_BTAGBEFFCORR;
     delete h2_BTAGLEFF;

   } // for btag algos


}


// ------------ method called once each job just before starting event loop  ------------
void 
BTagPayloadDumper::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagPayloadDumper::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
BTagPayloadDumper::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
BTagPayloadDumper::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
BTagPayloadDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BTagPayloadDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BTagPayloadDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagPayloadDumper);
