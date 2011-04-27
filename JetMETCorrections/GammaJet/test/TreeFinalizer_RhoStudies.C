#include <iostream>
#include <vector>
#include <cmath>

#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "TQObject.h"
//#include "RQ_OBJECT.h"

#include "fitTools.C"


struct photonidegcuts {

  float hovereiso;
  float hcaliso_rel;
  float hcaliso_abs;
  float ecaliso_rel;
  float ecaliso_abs;
  float trackiso_rel;
  float trackiso_abs;
  float setaetaEB;
  float setaetaEE;

};


class AnalysisPhoton : public TLorentzVector {

//RQ_OBJECT("AnalysisPhoton")
ClassDef (AnalysisPhoton, 0) 

public:

  AnalysisPhoton( float x=0., float y=0., float z=0., float t=0.) : TLorentzVector( x, y, z, t ) {};


  Float_t hcalIso;
  Float_t ecalIso;
  Int_t   nTrkIso;
  Float_t ptTrkIso;
  Float_t clusterMaj;
  Float_t clusterMin;
  Int_t hasPixelSeed;

  Float_t pid_hlwTrack;
  Float_t pid_twrHCAL;
  Float_t pid_HoverE;
  Float_t pid_jurECAL;
  Float_t pid_etawid;

};


bool DEBUG_ = false;

TChain* tree;


void addInput(const std::string& dataset);



bool passedPhotonID(const AnalysisPhoton & photon, photonidegcuts const& pid, vector<bool> *vpass);

void finalize(const std::string& dataset, std::string photonID="medium" ) {



  tree = new TChain("jetTree");



  std::string infileName, treeName;




  if( dataset=="G_TuneZ2_7TeV_pythia6" ) {

    addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
  //addInput( "G_Pt_800to1400_TuneZ2_7TeV_pythia6" );
  //addInput( "G_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="G_TuneZ2_7TeV_pythia6_CORR" ) {

    addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6_CORR" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6_CORR" );

  } else if( dataset=="QCD_TuneZ2_7TeV_pythia6" ) {

    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="G_QCD_TuneZ2_7TeV_pythia6" ) {

    //addInput( "G_Pt_0to15_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "G_Pt_470to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_15to30_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_30to50_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_50to80_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_80to120_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_120to170_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_170to300_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_300to470_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_470to600_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_600to800_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_800to1000_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6" );
    addInput( "QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6" );

  } else if( dataset=="DATA_EG_35X" ) {

    addInput( "MinimumBias_Commissioning10_May6thPDSkim2_SD_EG" );
    addInput( "MinimumBias_Commissioning10_SD_EG-v9" );
    addInput( "EG_Run2010A-PromptReco-v1" );
    addInput( "EG_Run2010A-PromptReco-v2" );

  } else if( dataset=="DATA_EG" ) {

    addInput( "MinimumBias_Commissioning10_SD_EG_Jun14thSkim_v1" );
    addInput( "EG_Run2010A-Jun14thReReco_v1" );
    addInput( "EG_Run2010A-PromptReco-v4" );
    addInput( "EG_Run2010A_PromptReco_v4_139347_139375" );
    addInput( "EG_Run2010A_PromptReco_v4_139376_139459" );

  } else if( dataset=="DATA_EG_37X" ) {

    addInput( "EG_Run2010A_Jul15thReReco_v1" );
    addInput( "EG_Run2010A_Jul26thReReco_v1" );

  } else if( dataset=="EG_3pb" ) {

    addInput( "EG_Run2010A-PromptReco-v4" );
    addInput( "EG_Run2010A-PromptReco-v4_runs143337-144114" );

  } else if( dataset=="DATA_Nov4ReReco" ) {

    addInput( "EG_Run2010A-Nov4ReReco_v1" );
    addInput( "Photon_Run2010B_Nov4ReReco" );

  } else {

    addInput( dataset );

  }




  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  Float_t ptMax = 100.;
  int nBins_photBinning = ptPhot_binning.size()-1;
  Double_t ptPhot_binning_array[100];
  for( unsigned i=0; i<nBins_photBinning; ++i)
    ptPhot_binning_array[i] = ptPhot_binning[i];



  Int_t run;
  tree->SetBranchAddress("run", &run);
  Int_t nvertex;
  tree->SetBranchAddress("nvertex", &nvertex);
  Int_t event;
  tree->SetBranchAddress("event", &event);
  Float_t eventWeight;
  tree->SetBranchAddress("eventWeight", &eventWeight);

  Float_t eMet;
  tree->SetBranchAddress("epfMet", &eMet);
  Float_t phiMet;
  tree->SetBranchAddress("phipfMet", &phiMet);


  Float_t rhoCalo;
  tree->SetBranchAddress("rhoCalo", &rhoCalo);
  Float_t rhoPF;
  tree->SetBranchAddress("rhoPF", &rhoPF);

  Float_t ptHat;
  tree->SetBranchAddress("ptHat", &ptHat);

  Float_t ePhotReco;
  tree->SetBranchAddress("ePhotReco", &ePhotReco);
  Float_t ptPhotReco;
  tree->SetBranchAddress("ptPhotReco", &ptPhotReco);
  Float_t etaPhotReco;
  tree->SetBranchAddress("etaPhotReco", &etaPhotReco);
  Float_t phiPhotReco;
  tree->SetBranchAddress("phiPhotReco", &phiPhotReco);

  Float_t ePhotGen;
  tree->SetBranchAddress("ePhotGen", &ePhotGen);
  Float_t ptPhotGen;
  tree->SetBranchAddress("ptPhotGen", &ptPhotGen);
  Float_t etaPhotGen;
  tree->SetBranchAddress("etaPhotGen", &etaPhotGen);
  Float_t phiPhotGen;
  tree->SetBranchAddress("phiPhotGen", &phiPhotGen);

  Float_t hcalIsoPhotReco;
  tree->SetBranchAddress("hcalIsoPhotReco", &hcalIsoPhotReco);
  Float_t ecalIsoPhotReco;
  tree->SetBranchAddress("ecalIsoPhotReco", &ecalIsoPhotReco);
  Int_t nTrkIsoPhotReco;
  tree->SetBranchAddress("nTrkIsoPhotReco", &nTrkIsoPhotReco);
  Float_t ptTrkIsoPhotReco;
  tree->SetBranchAddress("ptTrkIsoPhotReco", &ptTrkIsoPhotReco);
  Float_t clusterMajPhotReco;
  tree->SetBranchAddress("clusterMajPhotReco", &clusterMajPhotReco);
  Float_t clusterMinPhotReco;
  tree->SetBranchAddress("clusterMinPhotReco", &clusterMinPhotReco);
  Float_t hasPixelSeedPhotReco;
  tree->SetBranchAddress("hasPixelSeedPhotReco", &hasPixelSeedPhotReco);

  Float_t pid_towerHCALPhotReco;
  tree->SetBranchAddress("pid_towerHCALPhotReco", &pid_towerHCALPhotReco);
  Float_t pid_HoverEPhotReco;
  tree->SetBranchAddress("pid_HoverEPhotReco", &pid_HoverEPhotReco);
  Float_t pid_jurECALPhotReco;
  tree->SetBranchAddress("pid_jurECALPhotReco", &pid_jurECALPhotReco);
  Float_t pid_hlwTrackPhotReco;
  tree->SetBranchAddress("pid_hlwTrackPhotReco", &pid_hlwTrackPhotReco);
  Float_t pid_sIEtaIEtaPhotReco;
  tree->SetBranchAddress("pid_sIEtaIEtaPhotReco", &pid_sIEtaIEtaPhotReco);

  Bool_t matchedToMC;
  tree->SetBranchAddress("matchedToMC", &matchedToMC);

  Float_t eJetReco;
  tree->SetBranchAddress("eJetReco", &eJetReco);
  Float_t ptJetReco;
  tree->SetBranchAddress("ptJetReco", &ptJetReco);
  Float_t ptCorrJetReco;
  tree->SetBranchAddress("ptCorrJetReco", &ptCorrJetReco);
  Float_t etaJetReco;
  tree->SetBranchAddress("etaJetReco", &etaJetReco);
  Float_t phiJetReco;
  tree->SetBranchAddress("phiJetReco", &phiJetReco);
  Float_t eTracksReco;
  tree->SetBranchAddress("eTracksReco", &eTracksReco);

  Float_t eJetGen;
  tree->SetBranchAddress("eJetGen", &eJetGen);
  Float_t ptJetGen;
  tree->SetBranchAddress("ptJetGen", &ptJetGen);
  Float_t etaJetGen;
  tree->SetBranchAddress("etaJetGen", &etaJetGen);
  Float_t phiJetGen;
  tree->SetBranchAddress("phiJetGen", &phiJetGen);
  Float_t eTracksGen;
  tree->SetBranchAddress("eTracksGen", &eTracksGen);

  Float_t ptPart;
  tree->SetBranchAddress("ptPart", &ptPart);
  Float_t etaPart;
  tree->SetBranchAddress("etaPart", &etaPart);
  Float_t phiPart;
  tree->SetBranchAddress("phiPart", &phiPart);
  Int_t pdgIdPart;
  tree->SetBranchAddress("pdgIdPart", &pdgIdPart);

  Float_t pt2ndJetReco;
  tree->SetBranchAddress("pt2ndJetReco", &pt2ndJetReco);
  Float_t ptCorr2ndJetReco;
  tree->SetBranchAddress("ptCorr2ndJetReco", &ptCorr2ndJetReco);
  Float_t eta2ndJetReco;
  tree->SetBranchAddress("eta2ndJetReco", &eta2ndJetReco);
  Float_t phi2ndJetReco;
  tree->SetBranchAddress("phi2ndJetReco", &phi2ndJetReco);

  Float_t pt2ndJetGen;
  tree->SetBranchAddress("pt2ndJetGen", &pt2ndJetGen);
  Float_t eta2ndJetGen;
  tree->SetBranchAddress("eta2ndJetGen", &eta2ndJetGen);
  Float_t phi2ndJetGen;
  tree->SetBranchAddress("phi2ndJetGen", &phi2ndJetGen);

  Float_t ptSecondaryJetsReco;
  tree->SetBranchAddress("ptSecondaryJetsReco", &ptSecondaryJetsReco);
  Float_t ptSecondaryJetsGen;
  tree->SetBranchAddress("ptSecondaryJetsGen", &ptSecondaryJetsGen);

  Bool_t passed_Photon10;
  tree->SetBranchAddress("passed_Photon10", &passed_Photon10);
  Bool_t passed_Photon15;
  tree->SetBranchAddress("passed_Photon15", &passed_Photon15);
  Bool_t passed_Photon20;
  tree->SetBranchAddress("passed_Photon20", &passed_Photon20);
  Bool_t passed_Photon25;
  tree->SetBranchAddress("passed_Photon25", &passed_Photon25);
  Bool_t passed_Photon30;
  tree->SetBranchAddress("passed_Photon30", &passed_Photon30);
  Bool_t passed_Photon35;
  tree->SetBranchAddress("passed_Photon35", &passed_Photon35);


  int nBinsRho = 10;
  float rhoCalo_max = 6.;
  float rhoPF_max = 16.;


  TProfile* hp_hcalIsoMean_vs_rhoCalo = new TProfile("hcalIsoMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_hcalIsoMean_vs_rhoCalo->Sumw2();
  TProfile* hp_hcalIsoPIDMean_vs_rhoCalo = new TProfile("hcalIsoPIDMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_hcalIsoPIDMean_vs_rhoCalo->Sumw2();
  TProfile* hp_hcalIsoMean_vs_rhoPF = new TProfile("hcalIsoMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_hcalIsoMean_vs_rhoPF->Sumw2();
  TProfile* hp_hcalIsoPIDMean_vs_rhoPF = new TProfile("hcalIsoPIDMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_hcalIsoPIDMean_vs_rhoPF->Sumw2();

  TProfile* hp_ecalIsoMean_vs_rhoCalo = new TProfile("ecalIsoMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_ecalIsoMean_vs_rhoCalo->Sumw2();
  TProfile* hp_ecalIsoPIDMean_vs_rhoCalo = new TProfile("ecalIsoPIDMean_vs_rhoCalo", "", nBinsRho, 0., rhoCalo_max);
  hp_ecalIsoPIDMean_vs_rhoCalo->Sumw2();
  TProfile* hp_ecalIsoMean_vs_rhoPF = new TProfile("ecalIsoMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_ecalIsoMean_vs_rhoPF->Sumw2();
  TProfile* hp_ecalIsoPIDMean_vs_rhoPF = new TProfile("ecalIsoPIDMean_vs_rhoPF", "", nBinsRho, 0., rhoPF_max);
  hp_ecalIsoPIDMean_vs_rhoPF->Sumw2();



  int nEntries = tree->GetEntries();
//nEntries = 100000;

  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);


    if( eventWeight <= 0. ) eventWeight = 1.;


    AnalysisPhoton photon;
    photon.SetPtEtaPhiE( ptPhotReco, etaPhotReco, phiPhotReco, ePhotReco );
    photon.hcalIso = hcalIsoPhotReco;
    photon.ecalIso = ecalIsoPhotReco;
    photon.nTrkIso = nTrkIsoPhotReco;
    photon.ptTrkIso = ptTrkIsoPhotReco;
    photon.clusterMaj = clusterMajPhotReco;
    photon.clusterMin = clusterMinPhotReco;
    photon.hasPixelSeed = hasPixelSeedPhotReco;
   
    photon.pid_hlwTrack = pid_hlwTrackPhotReco;
    photon.pid_twrHCAL = pid_towerHCALPhotReco;
    photon.pid_HoverE = pid_HoverEPhotReco;
    photon.pid_jurECAL = pid_jurECALPhotReco;
    photon.pid_etawid = pid_sIEtaIEtaPhotReco;




    if( ptPhotReco<25. ) continue;
    //if( fabs(etaPhotReco)>1.3 ) continue;
    if( clusterMinPhotReco<0.15 ) continue; //protection vs EB spikes



    // match photon to MC:
    TLorentzVector photonGen;
    photonGen.SetPtEtaPhiE( ptPhotGen, etaPhotGen, phiPhotGen, ePhotGen );

    bool isMatched = photon.DeltaR(photonGen)<0.5;

    if( !isMatched ) continue;
    

    //------------
    // PHOTON ID:
    //------------

    photonidegcuts hggtightid;
    hggtightid.hovereiso=           0.02;
    hggtightid.hcaliso_rel=         0.0025;
    hggtightid.hcaliso_abs=         2.;
    hggtightid.ecaliso_rel=         0.006;
    hggtightid.ecaliso_abs=         2.;
    hggtightid.trackiso_rel=        0.001;
    hggtightid.trackiso_abs=        1.5;
    hggtightid.setaetaEB=           0.010;
    hggtightid.setaetaEE=           0.028;
    
    std::vector<bool> idpasseg(5);

    bool passedID = passedPhotonID(photon, hggtightid, &idpasseg);
    passedID = ( passedID && !hasPixelSeedPhotReco );

    //if( !passedID ) continue; 

    float pid_hcalIso = photon.pid_HoverE + photon.pid_twrHCAL / photon.Pt();
    float pid_ecalIso = photon.pid_jurECAL*cosh(photon.Eta()) / photon.Energy();
  
    hp_hcalIsoMean_vs_rhoCalo->Fill( rhoCalo, hcalIsoPhotReco, eventWeight );
    hp_hcalIsoPIDMean_vs_rhoCalo->Fill( rhoCalo, pid_hcalIso, eventWeight );
    hp_hcalIsoMean_vs_rhoPF->Fill( rhoPF, hcalIsoPhotReco, eventWeight );
    hp_hcalIsoPIDMean_vs_rhoPF->Fill( rhoPF, pid_hcalIso, eventWeight );

    hp_ecalIsoMean_vs_rhoCalo->Fill( rhoCalo, ecalIsoPhotReco, eventWeight );
    hp_ecalIsoPIDMean_vs_rhoCalo->Fill( rhoCalo, pid_ecalIso, eventWeight );
    hp_ecalIsoMean_vs_rhoPF->Fill( rhoPF, ecalIsoPhotReco, eventWeight );
    hp_ecalIsoPIDMean_vs_rhoPF->Fill( rhoPF, pid_ecalIso, eventWeight );



  

  } //for entries



  std::string outfileName;

  if( DEBUG_ ) outfileName = "provaRhoStudies_"+dataset;
  else {
   if(dataset!="") outfileName = "RhoStudies_"+dataset;
   else outfileName = "PhotonJet";
  }

  outfileName = outfileName;
  if( photonID!="medium" ) outfileName = outfileName + "_" + photonID;
  outfileName += ".root";

  TFile* outFile = new TFile(outfileName.c_str(), "RECREATE");
  outFile->cd();

  hp_hcalIsoMean_vs_rhoCalo->Write();
  hp_hcalIsoPIDMean_vs_rhoCalo->Write();
  hp_hcalIsoMean_vs_rhoPF->Write();
  hp_hcalIsoPIDMean_vs_rhoPF->Write();

  hp_ecalIsoMean_vs_rhoCalo->Write();
  hp_ecalIsoPIDMean_vs_rhoCalo->Write();
  hp_ecalIsoMean_vs_rhoPF->Write();
  hp_ecalIsoPIDMean_vs_rhoPF->Write();



  outFile->Close();

  delete tree;
  tree = 0;


}


void addInput( const std::string& dataset ) {


  // opening from filelist now deprecated (files have to be merged with merge_and_setWeights
  //std::string infileName = "files_PhotonJet_2ndLevel_" + dataset+"_" + ALGOTYPE_ +".txt";


  //FILE* iff = fopen(infileName.c_str(),"r");
  //if(iff == 0) {
    //std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::string infileName = "PhotonJet_2ndLevelTreeW_" + dataset + "_pfakt5";
    infileName += ".root";
    TFile* infile = TFile::Open(infileName.c_str(), "read");
    if( infile==0 ) {
      std::cout << "Didn't find file '" << infileName << "'. Did you forget to finalize (i.e. the \"W\")?" << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(77);
    }
    std::string treeName = infileName +"/jetTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    infile->Close();


} //addinput


  

bool passedPhotonID(const AnalysisPhoton& photon, photonidegcuts const& pid, vector<bool> *vpass) {

 // Use photon supercluster energy (would be e5x5 if r9>0.93 otherwise)
 bool ptiso =   (photon.pid_hlwTrack < photon.Pt() * pid.trackiso_rel + pid.trackiso_abs);
 bool ecaliso = (photon.pid_jurECAL  < photon.Pt() * pid.ecaliso_rel + pid.ecaliso_abs);
 bool hcaliso = (photon.pid_twrHCAL  < photon.Pt() * pid.hcaliso_rel + pid.hcaliso_abs);
 bool hoveiso = (photon.pid_HoverE < pid.hovereiso);
 bool setaeta =  photon.pid_etawid < pid.setaetaEB;
 if(TMath::Abs(photon.Eta()) > 1.44) {
   setaeta = photon.pid_etawid < pid.setaetaEE;
 }

 if (vpass) {
   //assert((*vpass).size()==7);
   if((*vpass).size()!=5) { cout << "major failure! (*vpass).size()!=7.. die!" << endl; exit(0) ; }
   (*vpass)[0] = ptiso;
   (*vpass)[1] = ecaliso;
   (*vpass)[2] = hcaliso;
   (*vpass)[3] = hoveiso;
   (*vpass)[4] = setaeta;
 }

 return (ptiso && hcaliso && ecaliso && hoveiso && setaeta);
}

