#include <iostream>
#include "TString.h"
#include "RedNtpFinalizer.h"



RedNtpFinalizer::RedNtpFinalizer( const std::string& analyzerType, const std::string& dataset, const std::string& flags) {

  DEBUG_ = false;


  tree_ = new TChain("AnaTree");

  analyzerType_ = analyzerType;
  //inputAnalyzerType_ = analyzerType;
  redNtpDir_ = "";
  dataset_ = dataset;
  flags_ = flags;

  outFile_ = 0;
  outputDir_ = "";

  nGenEvents_ = 0.;
  xSection_ = get_crossSection(dataset);

  h1_nPU_gen_ = new TH1F("nPU_gen", "", 55, -0.5, 54.5);

} //constructor



RedNtpFinalizer::~RedNtpFinalizer() {

  std::cout << std::endl << "-> Histograms saved in: " << outFile_->GetName() << std::endl;

  if( tree_!=0 ) {
    delete tree_;
    tree_=0;
  }

  delete h1_nPU_gen_;

} //destructor



void RedNtpFinalizer::createOutputFile( const std::string& additionalFlags ) {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova_"+dataset_;
   else {
    if(dataset_!="") outfileName =  analyzerType_ + "_" + dataset_;
    else outfileName = analyzerType_;
   }


   if( flags_!="" )
     outfileName = outfileName + "_" + flags_;
   if( additionalFlags!="" )
     outfileName = outfileName + "_" + additionalFlags;

   outfileName = outfileName + ".root";

   if( outputDir_!="" ) {
     outfileName = outputDir_ + "/" + outfileName;
     std::string mkdir_command = "mkdir -p " + outputDir_;
     system( mkdir_command.c_str() );
   }

   outFile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outFile_->cd();


}




void RedNtpFinalizer::clear() {

  if( tree_!=0 ) delete tree_;
  tree_ = new TChain("AnaTree");

  if( h1_nPU_gen_!=0 ) delete h1_nPU_gen_;
  h1_nPU_gen_ = new TH1F("nPU_gen", "", 55, -0.5, 54.5);

}


void RedNtpFinalizer::addFile(const std::string& dataset, const std::string& selection) {

  std::string infileName = redNtpDir_ + "/redntp_" + dataset + ".root"; 
  TFile* infile = TFile::Open(infileName.c_str(), "READ");
  if( infile==0 ) {
    std::cout << "---> Didn't find file: " << infileName << std::endl;
    std::cout << "---> Exiting!!" << std::endl;
    exit(11);
  }
  infile->cd();
  std::string treeName = infileName +"/AnaTree";
  tree_->Add(treeName.c_str());
  std::cout << "-> Added " << treeName << ". Tree now has " << tree_->GetEntries() << " entries." << std::endl;
  TH1F* h1_nGenEvents = (TH1F*)infile->Get("ptphotgen1");
  //h1_nPU_gen_->Add((TH1F*)infile->Get("nPU_gen"));
  if( h1_nGenEvents!= 0 ) {
    nGenEvents_ += h1_nGenEvents->GetEntries();
  } else {
    std::cout << std::endl << std::endl << "WARNING!! Dataset '" << dataset << "' has no nGenEvents information!!!" << std::endl;
  }
  infile->Close();

}





int RedNtpFinalizer::get_nBTags( const AnalysisJet& jet1, const AnalysisJet& jet2, BTagSFUtil* btsfutil, bool loosebtags ) {

  int nBTags;

  bool jet1_tagged_medium = jet1.btag_medium();
  bool jet1_tagged_loose  = jet1.btag_loose();
  bool jet2_tagged_medium = jet2.btag_medium();
  bool jet2_tagged_loose  = jet2.btag_loose();

  btsfutil->modifyBTagsWithSF( jet1_tagged_loose, jet1_tagged_medium, jet1.Pt(), jet1.Eta(), jet1.pdgIdPart );
  btsfutil->modifyBTagsWithSF( jet2_tagged_loose, jet2_tagged_medium, jet2.Pt(), jet2.Eta(), jet2.pdgIdPart );

  if( loosebtags ) {

    bool twoBTags  = ( jet1_tagged_medium && jet2_tagged_loose  )
                  || ( jet1_tagged_loose  && jet2_tagged_medium );
    bool oneBTag   = (!twoBTags) && ( jet1_tagged_loose || jet2_tagged_loose );

    if( twoBTags ) nBTags=2;
    else if( oneBTag ) nBTags=1;
    else nBTags=0;

  } else {

    bool twoBTags  = ( jet1_tagged_medium && jet2_tagged_medium );
    bool oneBTag   = (!twoBTags) && ( jet1_tagged_medium || jet2_tagged_medium );

    if( twoBTags ) nBTags=2;
    else if( oneBTag ) nBTags=1;
    else nBTags=0;

  }

  return nBTags;

}


float RedNtpFinalizer::get_crossSection(const std::string& dataset) const {


  TString dataset_tstr(dataset);


  // k factors - same scale factors for 7 TeV and 8 TeV
  double kfactordiphot = 1.3;                // for prompt and box 
  double kfactordiphotmadgraph = 1.15;       // for madgraph di-jets
  double kfactorgamjet = 1.3;
  double kfactorqcd = 1;
  double kfactordy = 1.;
  //double kfactordy = 1.15;
  double kfactorwg = 1;       
  double kfactorzg = 1;       
  double kfactorwmgg = 1;     
  double kfactorwpgg = 1;     
  double kfactorzgg = 1;      
  double kfactorttgg = 1;     
  double kfactorttjets = 1;   
  double kfactorwjets = 1;    
  double kfactorww = 1;       
  double kfactorwz = 1;       
  double kfactorzz = 1;     



  float xSection = 1.; 


  if( dataset_tstr.Contains("DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph") ) {
    xSection = 3532.81; 
  } else if( dataset_tstr.Contains("DYJetsToLL_M-10To50filter_8TeV-madgraph") ) {
    xSection = 860.5013;
  } else if( dataset_tstr.Contains("T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola") ) {
    xSection = 11.1773;
  } else if( dataset_tstr.Contains("Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola") ) {
    xSection = 11.1773;
  } else if( dataset_tstr.Contains("TTJets_TuneZ2star_8TeV-madgraph-tauola") || dataset_tstr.BeginsWith("TT_CT10") ) {
    xSection = 225.1967;
  } else if( dataset_tstr.Contains("WJetsToLNu_TuneZ2Star_8TeV-madgraph") ) {
    xSection = 37509.;
  } else if( dataset_tstr.Contains("WW_TuneZ2star_8TeV_pythia6_tauola") ) {
    xSection = 5.8123/(2.*0.108*2.*0.108);
  } else if( dataset_tstr.Contains("WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola") || dataset_tstr.Contains("WWTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola") ) {
    xSection = 5.8123;
  } else if( dataset_tstr.Contains("WZ_TuneZ2star_8TeV_pythia6_tauola") ) {
    xSection = 22.4486;
  } else if( dataset_tstr.Contains("WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola") ) {
    xSection = 0.7346;
  } else if( dataset_tstr.Contains("ZZ_TuneZ2star_8TeV_pythia6_tauola") ) {
    xSection = 9.0314;
  } else if( dataset_tstr.Contains("ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola") ) {
    xSection = 0.3649;
  } else if( dataset_tstr.Contains("ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola") ) {
    xSection = 1.2752;
  } else if( dataset_tstr.Contains("ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola") ) {
    xSection = 0.0921;
  } else if( dataset_tstr.Contains("WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola") ) {
    xSection = 1.0575;
  } else if( dataset_tstr.Contains("TTTo2L2Nu2B_8TeV-powheg-pythia6") ) {
    xSection = 23.6402;
  } else if( dataset_tstr.Contains("WGToLNuG_TuneZ2star_8TeV-madgraph-tauola") ) {
    xSection = 553.92;
  } else if( dataset_tstr.Contains("ZGToLLG_8TeV-madgraph") ) {
    xSection = 181.338*0.06;
  } else if( dataset_tstr.Contains("ZG_Inclusive_8TeV-madgraph") ) {
    xSection = 181.338;
  } else if( dataset_tstr.Contains("TTWJets_8TeV-madgraph") ) {
    xSection = 0.232;
  } else if( dataset_tstr.Contains("TTZJets_8TeV-madgraph") ) {
    xSection = 0.208;
  } else if( dataset_tstr.Contains("TTGJets_8TeV-madgraph") ) {
    xSection = 1.444; //LO from PREP
  } else if( dataset_tstr.Contains("DiPhotonBox_Pt-10To25_8TeV-pythia6") ) {
    xSection = 424.8 * kfactordiphot;                
  } else if( dataset_tstr.Contains("DiPhotonBox_Pt-25To250_8TeV-pythia6") ) {
    xSection = 15.54 * kfactordiphot;                
  } else if( dataset_tstr.Contains("DiPhotonBox_Pt-250ToInf_8TeV-pythia6") ) {
    xSection = 0.029038 * kfactordiphot;  
  } else if( dataset_tstr.Contains("DiPhotonJets_8TeV-madgraph") ) {
    xSection = 81. * kfactordiphotmadgraph;  
  } else if( dataset_tstr.Contains("GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6") ) {
    xSection = 0.001835 * 81930.0 * kfactorgamjet; 
  } else if( dataset_tstr.Contains("GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6") ) {
    xSection = 0.05387 * 8884.0 * kfactorgamjet;
  } else if( dataset_tstr.Contains("QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6") ) {
    xSection = 0.000235 * 5.195e+07 * kfactorqcd;
  } else if( dataset_tstr.Contains("QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6") ) {
    xSection = 0.002175 * 2.365e+07 * kfactorqcd;
  } else if( dataset_tstr.Contains("WH_ZH_HToGG_M-125_8TeV-pythia6") ) {
    xSection = (0.6966 + 0.3943)*2.28E-03;
  } else if( dataset_tstr.Contains("TTH_HToGG_M-125_8TeV-pythia6") ) {
    xSection = (0.1302)*2.28E-03;
  } else if( dataset_tstr.Contains("VBF_HToGG_M-125_8TeV-powheg-pythia6") ) {
    xSection = (1.578)*2.28E-03;
  } else if( dataset_tstr.Contains("GluGluToHToGG_M-125_8TeV-powheg-pythia6") ) {
    xSection = (19.52)*2.28E-03;
  } else if( dataset_tstr.BeginsWith("ZGG") ) {
    xSection = 0.068 * kfactorzgg;        // Z+gg  - da AN FP
  } else if( dataset_tstr.BeginsWith("TTbarGG") ) {
    xSection = 0.001316 * kfactorttgg;    // tt+gg - da AN FP
  } else if( dataset_tstr.BeginsWith("WmGG") ) {
    xSection = 0.0504 * kfactorwmgg;      // W-gg  - da AN FP
  } else if( dataset_tstr.BeginsWith("WpGG") ) {
    xSection = 0.0667 * kfactorwpgg;      // W+gg  - da AN FP
  } else if( dataset_tstr.BeginsWith("tHqLeptonic") ) {
    xSection = 0.0152 * 0.108 * 3. * 2.28E-03;      // from http://arxiv.org/abs/1211.0499
    //xSection *= 1.4; // NLO scale factor taken from http://arxiv.org/pdf/1302.3856v2.pdf
  } else if( dataset_tstr.BeginsWith("tHqHadronic") ) {
    xSection = 0.0152 * (1. - 0.108 * 3.) * 2.28E-03;      // from http://arxiv.org/abs/1211.0499
  } else {
    std::cout << std::endl << std::endl;
    std::cout << "-> WARNING!! Dataset: '" << dataset << "' not present in database. Cross section unknown." << std::endl;
    std::cout << "-> Will set unitary weights." << std::endl;
  }


  return xSection;

}
