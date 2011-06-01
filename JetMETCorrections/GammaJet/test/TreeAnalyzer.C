#include "TreeAnalyzer.h"
#include "TH1F.h"
#include "TRegexp.h"
#include <cstdlib>
#include <fstream>



TreeAnalyzer::TreeAnalyzer(const std::string& analyzerType, const std::string& dataset, const std::string& recoType, const std::string& jetAlgo, const std::string& flags, TTree* tree)
{

   dataset_ = dataset;

   DEBUG_ = false;
   filterGoodRuns_ = false; //default: do not filter
   totalIntLumi_ = 0.;

   analyzerType_ = analyzerType;
   recoType_ = recoType;
   jetAlgo_ = jetAlgo;

   if( recoType!="calo" && recoType!="jpt" && recoType!="pf" ) {
     std::cout << "Unknown reco type (" << recoType_ << ")!! Exiting." << std::endl;
     exit(33);
   } 

   if( jetAlgo!="akt5" && jetAlgo!="akt7" && jetAlgo!="kt4" && jetAlgo!="kt6" && jetAlgo!="ic5" ) {
     std::cout << "Unknown jet algo (" << jetAlgo << ")!! Exiting." << std::endl;
     exit(35);
   } 

   if( recoType=="jpt" && jetAlgo=="akt5" ) {
     algoType_ = "jptak5"; //this is a bug in the ntuple producer
   } else {
     algoType_ = (recoType=="calo") ? jetAlgo : recoType+jetAlgo;
   }
   flags_ = flags;

}



// destructor
TreeAnalyzer::~TreeAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   std::cout << "Total Int Lumi: " << totalIntLumi_ << " ub-1." << std::endl;
   outfile_->cd();
   h1_lumi_->SetBinContent(1, totalIntLumi_);
   h1_lumi_->Write();
   h1_nCounter_->Write();
   jetTree_->Write();
   outfile_->Write();
   outfile_->Close();
   

}


void TreeAnalyzer::LoadInput() {

   std::cout << "-> Loading input... (could take a while)" << std::endl;

   char treePath[400];
   TChain * chain = new TChain("pippo","");
   if( dataset_=="MinimumBias_Commissioning10-GOODCOLL-v9" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/MinimumBias/Commissioning10-GOODCOLL-v9_2ndTry/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="MinimumBias_Commissioning10_May6thPDSkim2_SD_EG" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/MinimumBias/Commissioning10_May6thPDSkim2_SD_EG/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="MinimumBias_Commissioning10-Apr20Skim_GOODCOLL-v1" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/MinimumBias/Commissioning10-Apr20Skim_GOODCOLL-v1_3rdTry/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="MinimumBias_Commissioning10_SD_EG-v9" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/MinimumBias/Commissioning10_SD_EG-v9/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Spring10_Pt0to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Spring10/PhotonJetPt0to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Spring10_Pt5to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Spring10/PhotonJetPt5to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Spring10_Pt15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Spring10/PhotonJetPt15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Spring10_Pt30" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Spring10/PhotonJetPt30/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Spring10_Pt80" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Spring10/PhotonJetPt80/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Spring10_Pt170" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Spring10/PhotonJetPt170/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="MinBias_357ReReco_v3"||dataset_=="MinBias_357ReReco_v3_Pt0to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/MinBias_357ReReco_v3/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt0to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt0to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt5to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt0to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt5to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt15to20" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt15to20/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt20to30" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt20to30/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt30" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt30/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt30to50" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt30/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt30to50/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt50to80" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt30/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt50to80/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt80" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt80/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt170" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt170/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt300" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt300/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt470" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt470/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Spring10_Pt1400" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Spring10/QCDPt1400/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="QCD_Pythia8" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/QCD_Pythia8/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="EG_Run2010A-PromptReco-v1" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/EG/Run2010A-PromptReco-v1/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="EG_Run2010A-PromptReco-v2" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/EG/Run2010A-PromptReco-v2/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="EG_Run2010A-PromptReco-v4" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/DATA/EG/Run2010A-PromptReco-v4/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt0to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt0to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt5to15" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt5to15/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt15to20" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt15to20/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt20to30" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt20to30/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt30to50" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt30to50/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt50to80" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt50to80/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt80to120" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt80to120/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt120to170" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt120to170/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="PhotonJet_Summer1036X_Pt170to300" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/PhotonJet_Summer1036X/PhotonJetPt170to300/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="Wenu_Summer10_START37_V5_S09_v1" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/Wenu_Summer10_START37_V5_S09_v1/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="H130ZZlljj_gluonfusion" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/H130ZZlljj_gluonfusion/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="H150ZZlljj_gluonfusion" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/H150ZZlljj_gluonfusion/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else if( dataset_=="H500ZZlljj_gluonfusion" ) {
     sprintf(treePath, "/cmsrm/pc21_2/pandolf/MC/H500ZZlljj_gluonfusion/output_*.root/myanalysis/pippo");
     chain->Add(treePath);
   } else {
     sprintf(treePath, "%s/output_*.root/myanalysis/pippo", dataset_.c_str());
     chain->Add(treePath);
   }
     

   if( chain->GetEntries()==0 ) {
     std::cout << "-> Looking for a list file..." << std::endl;
     std::string fileListName = "files_PhotonJet_" + dataset_ + ".txt";
     this->LoadInputFromFile( fileListName );
     return;
   }


   TTree* tree = chain;
   std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
   //std::cout << "-> Added " << nFiles << " files. Tree has " << tree->GetEntries() << " entries." << std::endl;

   CreateOutputFile();
   Init(tree);

}


void TreeAnalyzer::LoadInputFromFile( const std::string& fileName ) {

   FILE* iff = fopen(fileName.c_str(),"r");
   if(iff == 0) {
     std::cout << "cannot open input file " << fileName << " ... now exiting." << std::endl;
     exit(-1);
   }

   TChain * chain = new TChain("pippo","");

   char singleLine[500];
 
   while( fscanf(iff, "%s", singleLine) !=EOF ) {
   
     std::string singleLine_str(singleLine);
     singleLine_str = singleLine_str + "/myanalysis/pippo";
     std::cout << "-> Adding " << singleLine_str << std::endl;
     chain->Add(singleLine_str.c_str());

   }
   fclose(iff);

   TTree* tree = chain;
   std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
   //std::cout << "-> Added " << nFiles << " files. Tree has " << tree->GetEntries() << " entries." << std::endl;

   CreateOutputFile();
   Init(tree);
  

}

void TreeAnalyzer::CreateOutputFile() {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova2ndLevel_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_2ndLevelTree_"+dataset_;
    else outfileName = analyzerType_ + "_2ndLevelTree";
   }


   outfileName = outfileName + "_" + algoType_ + flags_ + ".root";

   outfile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outfile_->cd();

   jetTree_ = new TTree("jetTree", "Reduced Tree for Jet Studies");
   jetTree_->SetMaxTreeSize(100000000000ULL); //setting max tree size to 100 GB

   h1_lumi_ = new TH1F("lumi", "", 1, 0., 1.);
   h1_nCounter_ = new TH1F("nCounter", "", 1, 0., 1.);

}


Int_t TreeAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   if (!tree) return;

   GenEventParameters genPars = this->getGenEventParameters();
 //Float_t xsection = genPars.crossSection;
   ptHatMin_ = genPars.ptHatMin;
   ptHatMax_ = genPars.ptHatMax;


 ////will cut on pt_hat, so have to divide only by correct number of events:
 //char cutOnPtHat[70];
 //sprintf( cutOnPtHat, "genpt>%f && genpt<%f", ptHatMin_, ptHatMax_);
 //Float_t nEntries_cut = (Float_t)tree->GetEntries(cutOnPtHat);
 //eventWeight_ = ( (nEntries_cut>0.) && (xsection>0.) ) ? xsection/nEntries_cut : 1.;

   //will cut on pt_hat, so nCounter needs to be set to number of actual events which will be analyzed:
   char cutOnPtHat[70];
   sprintf( cutOnPtHat, "genpt>%lf && genpt<%lf", (Double_t)ptHatMin_, (Double_t)ptHatMax_);
   Int_t nEntries_cut = tree->GetEntries(cutOnPtHat);
   h1_nCounter_->SetBinContent( 1, nEntries_cut );




   // Set branch addresses and branch pointers
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   std::string branchName;

   fChain->SetBranchAddress("genpt", &genpt, &b_genpt);
   fChain->SetBranchAddress("isMC", &isMC, &b_isMC);
   fChain->SetBranchAddress("store", &store, &b_store);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("rhoCalo", &rhoCalo, &b_rhoCalo);
   fChain->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("pdgIdMC", &pdgIdMC, &b_pdgIdMC);
   fChain->SetBranchAddress("statusMC", &statusMC, &b_statusMC);
   fChain->SetBranchAddress("motherIDMC", &motherIDMC, &b_motherIDMC);
   fChain->SetBranchAddress("ptMC ", &ptMC , &b_ptMC );
   fChain->SetBranchAddress("eMC  ", &eMC  , &b_eMC  );
   fChain->SetBranchAddress("etaMC", &etaMC, &b_etaMC);
   fChain->SetBranchAddress("phiMC", &phiMC, &b_phiMC);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("ptPhot ", ptPhot , &b_ptPhot );
   fChain->SetBranchAddress("ePhot  ", ePhot  , &b_ePhot  );
   fChain->SetBranchAddress("escPhot  ", escPhot  , &b_escPhot  );
   fChain->SetBranchAddress("eseedPhot  ", eseedPhot  , &b_eseedPhot  );
   fChain->SetBranchAddress("etaPhot", etaPhot, &b_etaPhot);
   fChain->SetBranchAddress("phiPhot", phiPhot, &b_phiPhot);
   fChain->SetBranchAddress("timePhot", timePhot, &b_timePhot);
   fChain->SetBranchAddress("e4SwissCrossPhot", e4SwissCrossPhot, &b_e4SwissCrossPhot);
   fChain->SetBranchAddress("hasPixelSeedPhot", hasPixelSeedPhot, &b_hasPixelSeedPhot);
   fChain->SetBranchAddress("pid_isEM", pid_isEM, &b_pid_isEM);
   fChain->SetBranchAddress("pid_isLoose", pid_isLoose, &b_pid_isLoose);
   fChain->SetBranchAddress("pid_isTight", pid_isTight, &b_pid_isTight);
   fChain->SetBranchAddress("pid_jurECAL", pid_jurECAL, &b_pid_jurECAL);
   fChain->SetBranchAddress("pid_twrHCAL", pid_twrHCAL, &b_pid_twrHCAL);
   fChain->SetBranchAddress("pid_HoverE", pid_HoverE, &b_pid_HoverE);
   fChain->SetBranchAddress("pid_hlwTrack", pid_hlwTrack, &b_pid_hlwTrack);
   fChain->SetBranchAddress("pid_hlwTrackNoDz", pid_hlwTrackNoDz, &b_pid_hlwTrackNoDz);
   fChain->SetBranchAddress("pid_etawid", pid_etawid, &b_pid_etawid);
   fChain->SetBranchAddress("ptiso035Phot", ptiso035Phot, &b_ptiso035Phot);
   fChain->SetBranchAddress("ntrkiso035Phot", ntrkiso035Phot, &b_ntrkiso035Phot);
   fChain->SetBranchAddress("ptiso04Phot", ptiso04Phot, &b_ptiso04Phot);
   fChain->SetBranchAddress("ntrkiso04Phot", ntrkiso04Phot, &b_ntrkiso04Phot);
   fChain->SetBranchAddress("hcalovecal04Phot", hcalovecal04Phot, &b_hcalovecal04Phot);
   fChain->SetBranchAddress("ecaliso04Phot", ecaliso04Phot, &b_ecaliso04Phot);
   fChain->SetBranchAddress("sMajMajPhot", sMajMajPhot, &b_sMajMajPhot);
   fChain->SetBranchAddress("sMinMinPhot", sMinMinPhot, &b_sMinMinPhot);
   fChain->SetBranchAddress("alphaPhot", alphaPhot, &b_alphaPhot);
   fChain->SetBranchAddress("sEtaEtaPhot", sEtaEtaPhot, &b_sEtaEtaPhot);
   fChain->SetBranchAddress("sEtaPhiPhot", sEtaPhiPhot, &b_sEtaPhiPhot);
   fChain->SetBranchAddress("sPhiPhiPhot", sPhiPhiPhot, &b_sPhiPhiPhot);
   fChain->SetBranchAddress("E1Phot", E1Phot, &b_E1Phot);
   fChain->SetBranchAddress("E9Phot", E9Phot, &b_E9Phot);
   fChain->SetBranchAddress("E25Phot", E25Phot, &b_E25Phot);
   branchName = "nJet_"+algoType_;
   fChain->SetBranchAddress(branchName.c_str(), &nJet, &b_nJet);
   branchName = "ptJet_"+algoType_+" ";
   fChain->SetBranchAddress(branchName.c_str(), ptJet, &b_ptJet);
   branchName = "ptCorrJet_"+algoType_+" ";
   fChain->SetBranchAddress(branchName.c_str(), ptCorrJet, &b_ptCorrJet);
   branchName = "eJet_"+algoType_+"  ";
   fChain->SetBranchAddress(branchName.c_str(), eJet, &b_eJet);
   branchName = "etaJet_"+algoType_;
   fChain->SetBranchAddress(branchName.c_str(), etaJet, &b_etaJet);
   branchName = "phiJet_"+algoType_;
   fChain->SetBranchAddress(branchName.c_str(), phiJet, &b_phiJet);
   if( recoType_=="calo" ) {
     branchName = "n90Jet_"+algoType_;
     fChain->SetBranchAddress("n90Jet", n90Jet, &b_n90Jet);
     branchName = "n90HitsJet_"+algoType_;
     fChain->SetBranchAddress("n90HitsJet", n90HitsJet, &b_n90HitsJet);
     branchName = "fHPDJet_"+algoType_;
     fChain->SetBranchAddress("fHPDJet", fHPDJet, &b_fHPDJet);
     branchName = "fRBXJet_"+algoType_;
     fChain->SetBranchAddress("fRBXJet", fRBXJet, &b_fRBXJet);
   }
   if( recoType_=="calo"||recoType_=="jpt" ) {
     branchName = "emfJet_"+algoType_;
     fChain->SetBranchAddress("emfJet", emfJet, &b_emfJet);
   }
   if( recoType_=="pf" && ( jetAlgo_=="akt5" || jetAlgo_=="akt7" ) ) {
     branchName = "nChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nChargedHadrons, &b_nChargedHadrons);
     branchName = "nPhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nPhotons, &b_nPhotons);
     branchName = "nMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nMuons, &b_nMuons);
     branchName = "nElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nElectrons, &b_nElectrons);
     branchName = "nNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nNeutralHadrons, &b_nNeutralHadrons);
     branchName = "nHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nHFHadrons, &b_nHFHadrons);
     branchName = "nHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), nHFEM, &b_nHFEM);
     branchName = "eChargedHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eChargedHadrons, &b_eChargedHadrons);
     branchName = "ePhotons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), ePhotons, &b_ePhotons);
     branchName = "eMuons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eMuons, &b_eMuons);
     branchName = "eElectrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eElectrons, &b_eElectrons);
     branchName = "eNeutralHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eNeutralHadrons, &b_eNeutralHadrons);
     branchName = "eHFHadrons_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eHFHadrons, &b_eHFHadrons);
     branchName = "eHFEM_"+algoType_;
     fChain->SetBranchAddress(branchName.c_str(), eHFEM, &b_eHFEM);
     if( recoType_=="pf" && jetAlgo_=="akt5" ) {
       branchName = "ptDJet_"+algoType_;
       fChain->SetBranchAddress(branchName.c_str(), ptDJet, &b_ptDJet);
       branchName = "rmsCandJet_"+algoType_;
       fChain->SetBranchAddress(branchName.c_str(), rmsCandJet, &b_rmsCandJet);
     }
   }
   branchName = "nJetGen_"+jetAlgo_;
   fChain->SetBranchAddress(branchName.c_str(), &nJetGen, &b_nJetGen);
   branchName = "ptJetGen_"+jetAlgo_+" ";
   fChain->SetBranchAddress(branchName.c_str(), ptJetGen, &b_ptJetGen);
   branchName = "eJetGen_"+jetAlgo_+"  ";
   fChain->SetBranchAddress(branchName.c_str(), eJetGen, &b_eJetGen);
   branchName = "etaJetGen_"+jetAlgo_;
   fChain->SetBranchAddress(branchName.c_str(), etaJetGen, &b_etaJetGen);
   branchName = "phiJetGen_"+jetAlgo_;
   fChain->SetBranchAddress(branchName.c_str(), phiJetGen, &b_phiJetGen);
   if( jetAlgo_=="akt5" ) {
     branchName = "nMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nMuonsGen, &b_nMuonsGen);
     branchName = "nElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nElectronsGen, &b_nElectronsGen);
     branchName = "nPhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nPhotonsGen, &b_nPhotonsGen);
     branchName = "nTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nTracksGen, &b_nTracksGen);
     branchName = "nNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nNeutralHadronsGen, &b_nNeutralHadronsGen);
     branchName = "nHFHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nHFHadronsGen, &b_nHFHadronsGen);
     branchName = "nHFEMGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), nHFEMGen, &b_nHFEMGen);
     branchName = "eMuonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eMuonsGen, &b_eMuonsGen);
     branchName = "eElectronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eElectronsGen, &b_eElectronsGen);
     branchName = "ePhotonsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), ePhotonsGen, &b_ePhotonsGen);
     branchName = "eTracksGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eTracksGen, &b_eTracksGen);
     branchName = "eNeutralHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eNeutralHadronsGen, &b_eNeutralHadronsGen);
     branchName = "eHFHadronsGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eHFHadronsGen, &b_eHFHadronsGen);
     branchName = "eHFEMGen_"+jetAlgo_;
     fChain->SetBranchAddress(branchName.c_str(), eHFEMGen, &b_eHFEMGen);
   }
   fChain->SetBranchAddress("sMet  ", &sMet  , &b_sMet);
   fChain->SetBranchAddress("eMet  ", &eMet  , &b_eMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
   fChain->SetBranchAddress("sCorrMet  ", &sCorrMet  , &b_sCorrMet);
   fChain->SetBranchAddress("eCorrMet  ", &eCorrMet  , &b_eCorrMet);
   fChain->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
   fChain->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
   fChain->SetBranchAddress("smuCorrMet  ", &smuCorrMet  , &b_smuCorrMet);
   fChain->SetBranchAddress("emuCorrMet  ", &emuCorrMet  , &b_emuCorrMet);
   fChain->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
   fChain->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
   fChain->SetBranchAddress("sNoHFMet  ", &sNoHFMet  , &b_sNoHFMet);
   fChain->SetBranchAddress("eNoHFMet  ", &eNoHFMet  , &b_eNoHFMet);
   fChain->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
   fChain->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
   fChain->SetBranchAddress("stcMet  ", &stcMet  , &b_stcMet);
   fChain->SetBranchAddress("etcMet  ", &etcMet  , &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   fChain->SetBranchAddress("spfMet  ", &spfMet  , &b_spfMet);
   fChain->SetBranchAddress("epfMet  ", &epfMet  , &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   fChain->SetBranchAddress("sMetGen  ", &sMetGen  , &b_sMetGen);
   fChain->SetBranchAddress("eMetGen  ", &eMetGen  , &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   fChain->SetBranchAddress("sMetGen2  ", &sMetGen2  , &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2  ", &eMetGen2  , &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("nvertex", &nvertex, &b_nvertex);
   fChain->SetBranchAddress("vxMC", &vxMC, &b_vxMC);
   fChain->SetBranchAddress("vyMC", &vyMC, &b_vyMC);
   fChain->SetBranchAddress("vzMC", &vzMC, &b_vzMC);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("vntracks", vntracks, &b_vntracks);
   fChain->SetBranchAddress("vchi2", vchi2, &b_vchi2);
   fChain->SetBranchAddress("vndof", vndof, &b_vndof);
   //fChain->SetBranchAddress("hltPass", &hltPass, &b_hltPass);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   //fChain->SetBranchAddress("hltNamesLen", &hltNamesLen, &b_hltNamesLen);
   //fChain->SetBranchAddress("HLTNames", HLTNames, &b_HLTNames);
   //fChain->SetBranchAddress("HLTResults", HLTResults, &b_HLTResults);
   fChain->SetBranchAddress("HLTNames", &HLTNames);
   fChain->SetBranchAddress("HLTResults", &HLTResults);
   fChain->SetBranchAddress("Xsec", &Xsec, &b_Xsec);
   Notify();
}

Bool_t TreeAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeAnalyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void TreeAnalyzer::UpdateCache() {
     // cache current run
     if( oldrun_ != run ) {
       oldrun_ = run;
       goodLSCache_ = goodLS_.find( run );
     }
}


void TreeAnalyzer::ReadCSVFile(const std::string& csv) {

  std::cout << "Reading CSV file of LS luminosities " << csv << std::endl;
  FILE* iff = fopen(csv.c_str(),"r");

  if(iff == 0) {
    std::cout << "cannot open CSV file " << csv << " ... now exiting." << std::endl;
    exit(-2);
  }

  int irun, iLS;
  double iDelivered, iRecorded;

  //while( fscanf(iff,"%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%lf,%*[ \t\n]%lf", &irun, &iLS, &iDelivered, &iRecorded) != EOF ) {
  while( fscanf(iff,"%d,%d,%lf,%lf", &irun, &iLS, &iDelivered, &iRecorded) != EOF ) {

    RunLumiPair rlpair( std::pair<int, int>(irun, iLS) );
    //double istLumi = iRecorded;
    //double intLumi = 23.*istLumi/1.0e30; //in microb-1 (1 LS = 23 seconds)
    LSLumimap_[rlpair] = iRecorded; //in microb-1 

  }
  fclose(iff);

}


void TreeAnalyzer::ReadJSONFile(const std::string& json) {


  std::cout << "Reading JSON file of good runs " << json << std::endl;
  FILE* iff = fopen(json.c_str(),"r");

  if(iff == 0) {
    std::cout << "cannot open JSON file " << json << " ... now exiting." << std::endl;
    exit(-1);
  }

  char c1, c2, c3;
  int run1, run2, LS1, LS2;

  std::cout << "Following LS will be used" << std::endl;
  std::cout << "-------------------------" << std::endl;
  while( fscanf(iff,"%*[ \t\n]%c%d:%d-%d:%d%c%c",&c1,&run1,&LS1,&run2,&LS2,&c2,&c3 ) != EOF ) {
      std::cout << "run: " << run1 << "  LS range: " << LS1
         << " --> " << LS2 << std::endl;
      goodLS_[run1].push_back(  std::pair<int,int>(LS1,LS2) );
  }
  fclose(iff);
  filterGoodRuns_ = true; // will run only on good runs/LS

}


bool TreeAnalyzer::isGoodLS() {

     bool returnBool = false;

   //  if( filterGoodRuns_ ) std::cout << "is true!!" << std::endl;
   //  else std::cout << "is FALSE" << std::endl;

     if(!filterGoodRuns_) returnBool = true; // if filtered not requested all events are good

     this->UpdateCache();

   //  std::cout << "run: " << run << "\tlbn: " << lbn << std::endl;

     // check whether this run is part of the good runs. else retrun false
     if( goodLSCache_ != goodLS_.end() ) {

   //  std::cout << "in the loop" << std::endl;
        // get list of LS intervals
        const GoodLSVector& lsvector =   goodLSCache_->second; 
        // loop over good LS intervals and return as soon as one interval contains this event
        for(GoodLSVector::const_iterator iLS = lsvector.begin(); (iLS != lsvector.end())&&(returnBool==false); iLS++) {
      //     std::cout << "good range: " << iLS->first << "-->" << iLS->second;
           if(lbn >= iLS->first && lbn <= iLS->second ) {
             //cout << "Accepting run: " << Run << " LS: " << LumiSection << endl;
             returnBool = true;
          } // check current LS being in the interval
        } // loop over good LS for this run
     }

     if( returnBool==true ) {
     //  std::cout << "currentLS_: " << currentLS_ << "\tlbn: " << lbn << std::endl;
       if( currentLS_ != lbn )  {
          currentLS_ = lbn;
          RunLumiPair rlpair = (std::pair<int, int>(run,lbn));
          totalIntLumi_ += LSLumimap_[rlpair];
       }
     }

     return returnBool;

}


GenEventParameters TreeAnalyzer::getGenEventParameters() {

   GenEventParameters returnGenPars;
   returnGenPars.ptHatMin = -1.;

   if( dataset_=="PhotonJet_Summer09_Pt15" ) {
     returnGenPars.crossSection = 288813. - 32203.8;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="PhotonJet_Summer09_Pt30" ) {
     returnGenPars.crossSection = 32203.8 - 1012.08;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="PhotonJet_Summer09_Pt80" ) {
     returnGenPars.crossSection = 1012.08 - 51.36;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="PhotonJet_Summer09_Pt170" ) {
     returnGenPars.crossSection = 51.36 - 4.193;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="PhotonJet_Summer09_Pt300" ) {
     returnGenPars.crossSection =  4.193 - 0.45125;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="PhotonJet_Summer09_Pt470" ) {
     returnGenPars.crossSection =  0.45125 - 0.02;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="PhotonJet_Summer09_Pt800" ) {
     returnGenPars.crossSection = 0.02 - 0.000268;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="PhotonJet_Summer09_Pt1400" ) {
     returnGenPars.crossSection = 0.000268;
     returnGenPars.ptHatMax = 14000.;
   } else if( dataset_=="QCD_Summer09_Pt15" ) {
     returnGenPars.crossSection = 1458126879.8;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="QCD_Summer09_Pt30" ) {
     returnGenPars.crossSection = 109005537.31;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="QCD_Summer09_Pt80" ) {
     returnGenPars.crossSection = 1936120.4893;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="QCD_Summer09_Pt170" ) {
     returnGenPars.crossSection = 62508.776856;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="QCD_Summer09_Pt300" ) {
     returnGenPars.crossSection = 3669.4197667;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="QCD_Summer09_Pt470" ) {
     returnGenPars.crossSection = 315.32221016;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="QCD_Summer09_Pt800" ) {
     returnGenPars.crossSection = 11.94070485;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="QCD_Summer09_Pt1400" ) {
     returnGenPars.crossSection = 0.17207350709;
     returnGenPars.ptHatMax = 2200.;
   } else if( dataset_=="MC_PhotonJet_Summer09_Pt0to15" ) {
     returnGenPars.crossSection = 84460000.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="MC_PhotonJet_Summer09_Pt15to20" ) {
     returnGenPars.crossSection = 114700.;
     returnGenPars.ptHatMax = 20.;
   } else if( dataset_=="MC_PhotonJet_Summer09_Pt20to30" ) {
     returnGenPars.crossSection = 57180.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="MinBias_Spring10-START3X_V26A_356ReReco-v1"||dataset_=="MinBias_Spring10-START3X_V26A_357ReReco-v3"||dataset_=="MinBias_357ReReco_v3"||dataset_=="MinBias_357ReReco_v3_Pt0to15" ) {
     returnGenPars.crossSection = 71260000000.;
     returnGenPars.ptHatMax = ( dataset_=="MinBias_357ReReco_v3_Pt0to15" ) ? 15. : 7000.;
   } else if( dataset_=="PhotonJet_Spring10_Pt0to15" || dataset_ == "PhotonJet_Summer1036X_Pt0to15" ) {
     returnGenPars.crossSection = 84460000.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="PhotonJet_Spring10_Pt5to15" || dataset_=="PhotonJet_Summer1036X_Pt5to15" ) {
     returnGenPars.crossSection = 4030000.;
     returnGenPars.ptHatMax = 5.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="PhotonJet_Spring10_Pt15" ) {
     returnGenPars.crossSection = 192200.-20070.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="PhotonJet_Spring10_Pt30" ) {
     returnGenPars.crossSection = 20070.-556.5;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="PhotonJet_Spring10_Pt80" ) {
     returnGenPars.crossSection = 556.5-24.37;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="PhotonJet_Spring10_Pt170" ) {
     returnGenPars.crossSection = 24.37-1.636;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="PhotonJet_Spring10_Pt300" ) {
     returnGenPars.crossSection = 1.636-0.136;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="PhotonJet_Spring10_Pt470" ) {
     returnGenPars.crossSection = 0.136-0.003477;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="PhotonJet_Spring10_Pt800" ) {
     returnGenPars.crossSection = 0.003477-0.00001286;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="PhotonJet_Spring10_Pt1400" ) {
     returnGenPars.crossSection = 0.00001286;
     returnGenPars.ptHatMax = 7000.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt15to20" ) {
     returnGenPars.crossSection = 114700.;
     returnGenPars.ptHatMin = 15.;
     returnGenPars.ptHatMax = 20.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt20to30" ) {
     returnGenPars.crossSection = 57180.;
     returnGenPars.ptHatMin = 20.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt30to50" ) {
     returnGenPars.crossSection = 16520.;
     returnGenPars.ptHatMin = 30.;
     returnGenPars.ptHatMax = 50.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt50to80" ) {
     returnGenPars.crossSection = 2723.;
     returnGenPars.ptHatMin = 50.;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt80to120" ) {
     returnGenPars.crossSection = 446.2;
     returnGenPars.ptHatMin = 80.;
     returnGenPars.ptHatMax = 120.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt120to170" ) {
     returnGenPars.crossSection = 84.43;
     returnGenPars.ptHatMin = 120.;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="PhotonJet_Summer1036X_Pt170to300" ) {
     returnGenPars.crossSection = 22.55;
     returnGenPars.ptHatMin = 170.;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="QCD_Spring10_Pt0to15" ) {
     returnGenPars.crossSection = 48445000000.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="QCD_Spring10_Pt5to15" ) {
     returnGenPars.crossSection = 36640000000.;
     returnGenPars.ptHatMin = 5.;
     returnGenPars.ptHatMax = 15.;
   } else if( dataset_=="QCD_Spring10_Pt15" ) {
     returnGenPars.crossSection = 876215000.-60411000.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="QCD_Spring10_Pt15to20" ) {
     returnGenPars.crossSection = 579411000.;
     returnGenPars.ptHatMin = 15.;
     returnGenPars.ptHatMax = 20.;
   } else if( dataset_=="QCD_Spring10_Pt20to30" ) {
     returnGenPars.crossSection = 236051000.;
     returnGenPars.ptHatMin = 20.;
     returnGenPars.ptHatMax = 30.;
   } else if( dataset_=="QCD_Spring10_Pt30" ) {
     returnGenPars.crossSection = 60411000.-923821.;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="QCD_Spring10_Pt30to50" ) {
     returnGenPars.crossSection = 53114800.;
     returnGenPars.ptHatMin = 30.;
     returnGenPars.ptHatMax = 50.;
   } else if( dataset_=="QCD_Spring10_Pt50to80" ) {
     returnGenPars.crossSection = 6358210.;
     returnGenPars.ptHatMin = 50.;
     returnGenPars.ptHatMax = 80.;
   } else if( dataset_=="QCD_Spring10_Pt80" ) {
     returnGenPars.crossSection = 923821.-25474.9;
     returnGenPars.ptHatMax = 170.;
   } else if( dataset_=="QCD_Spring10_Pt170" ) {
     returnGenPars.crossSection = 25474.9-1255.87;
     returnGenPars.ptHatMax = 300.;
   } else if( dataset_=="QCD_Spring10_Pt300" ) {
     returnGenPars.crossSection = 1255.87-87.9799;
     returnGenPars.ptHatMax = 470.;
   } else if( dataset_=="QCD_Spring10_Pt470" ) {
     returnGenPars.crossSection = 87.9799-2.18608;
     returnGenPars.ptHatMax = 800.;
   } else if( dataset_=="QCD_Spring10_Pt800" ) {
     returnGenPars.crossSection = 2.18608-0.0112233;
     returnGenPars.ptHatMax = 1400.;
   } else if( dataset_=="QCD_Spring10_Pt1400" ) {
     returnGenPars.crossSection = 0.0112233;
     returnGenPars.ptHatMax = 10000.;
   } else if( dataset_=="Wenu_Summer10_START37_V5_S09_v1" ) {
     returnGenPars.crossSection = 7899.;
     returnGenPars.ptHatMax = 10000.;
   } else {
     std::cout << "WARNING! Dataset '" << dataset_ << "' has no hardwired cross-section, setting it to -1." << std::endl;
     returnGenPars.crossSection = -1.;
     returnGenPars.ptHatMax = 10000.;
   }

   if( returnGenPars.crossSection != -1 ) 
     std::cout << "-> Dataset was in database. Cross-section correctly set." << std::endl;

   return returnGenPars;

} // getGenEventParameters



bool TreeAnalyzer::passedTrigger_regexp( const std::string& trigger ) {

  TRegexp regexp(trigger);
  bool returnBool = false;

  for(unsigned iTrig=0; iTrig<HLTNames->size(); ++iTrig ) {

    TString hltname_tstr( HLTNames->at(iTrig) );
    if( hltname_tstr.Contains(regexp) ) {
      if( HLTResults->at(iTrig)==true ) returnBool = true;
      break; // logical OR requires at least one OK
    }

  } //for itrig

  return returnBool;

}
