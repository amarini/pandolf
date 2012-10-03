#include "Ntp1Analyzer.h"
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "TH1F.h"
#include "TRandom.h"
#include "TLorentzVector.h"




Ntp1Analyzer::Ntp1Analyzer(const std::string& analyzerType, const std::string& dataset, const std::string& flags, TTree* tree)
{

   dataset_ = dataset;

   DEBUG_ = false;
   filterGoodRuns_ = false; //default: do not filter
   totalIntLumi_ = 0.;

   analyzerType_ = analyzerType;

   flags_ = flags;

   ptHatMin_ = 0.;
   ptHatMax_ = 10000.;

   cachedLS_ = 0;
   cachedRun_ = 0;

   rand_ = new TRandom();

   index_requiredTriggers_.clear();
   index_notTriggers_.clear();

   requiredTriggers_.clear();
   notTriggers_.clear();

}



// destructor
Ntp1Analyzer::~Ntp1Analyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   outfile_->cd();
   h1_nCounter_->Write();
   h1_nCounterPU_->Write();
   h1_nCounterPU_ave_->Write();
   reducedTree_->Write();
   outfile_->Write();
   outfile_->Close();
   

}


void Ntp1Analyzer::LoadInput() {

   std::cout << "-> Loading input... (could take a while)" << std::endl;

   std::string treeDir;
   char treePath[400];
   TChain * chain = new TChain("ntp1","");
// if( dataset_=="Wenu_Summer10_START37_V5_S09_v1" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/Wenu_Summer10_START37_V5_S09_v1";
// } else if( dataset_=="HZZ_qqll_gluonfusion_M200" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M200";
// } else if( dataset_=="HZZ_qqll_gluonfusion_M300" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M300";
// } else if( dataset_=="HZZ_qqll_gluonfusion_M400" ) {
//   treeDir = "/cmsrm/pc21_2/pandolf/MC/HZZ_qqll_gluonfusion_M400";
// } else {
// }


   sprintf(treePath, "%s/default_*.root/ntp1", treeDir.c_str());

   int addInt = chain->Add(treePath);

   if( addInt==0 ) {
     std::cout << "Didn't find files to add for dataset: '" << dataset_ << "'. Looking for a list..." << std::endl;
     std::string fileName = "files_" + dataset_ + ".txt";
     this->LoadInputFromFile(fileName);
   } 
   /*else {
     TTree* tree = chain;
     std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
     this->CreateOutputFile();
     Init(tree);
     //load trigger mask:
     std::string firstFileName = treeDir + "/default_1.root";
     if( dataset_=="HZZ_qqll_gluonfusion_M300_CANDS" || dataset_=="HZZ_qqll_gluonfusion_M400_CANDS" )
       firstFileName = treePath;
     TFile* firstFile = TFile::Open( firstFileName.c_str(), "read" );
     this->LoadTrigger(firstFile);
   }*/

}


void Ntp1Analyzer::LoadInputFromFile( const std::string& fileName ) {

   FILE* iff = fopen(fileName.c_str(),"r");
   if(iff == 0) {
     std::cout << "cannot open input file " << fileName << " ... now exiting." << std::endl;
     exit(-1);
   }

   TChain * chain = new TChain("ntp1","");

   char singleLine[3000];
   bool isFirstFile=true;

   TFile* firstFile = 0;

   while( fscanf(iff, "%s", singleLine) !=EOF ) {
  
     std::string singleLine_str(singleLine);
     TString singleLine_tstr(singleLine);
     std::string treeName_str = singleLine_str + "/ntp1";
     if( singleLine_tstr.BeginsWith("/eos/") )
       treeName_str = "root://eoscms//" + treeName_str;
     std::cout << "-> Adding " << treeName_str << std::endl;
     chain->Add(treeName_str.c_str());
     if( isFirstFile ) {
       firstFile = TFile::Open(treeName_str.c_str(), "read");
       isFirstFile=false;
     }

   }
   fclose(iff);

   TTree* tree = chain;
   std::cout << "-> Tree has " << tree->GetEntries() << " entries." << std::endl;
   this->CreateOutputFile();
   Init(tree);
   //this->LoadTrigger(firstFile);

}



void Ntp1Analyzer::AddRequiredTrigger( const std::string& trigger, int runMin, int runMax ) { 

  TriggerMask newTriggerMask;
  newTriggerMask.HLTName = trigger;
  newTriggerMask.runMin = runMin;
  newTriggerMask.runMax = runMax;
  
  requiredTriggers_.push_back(newTriggerMask); 

}



void Ntp1Analyzer::AddRequiredTriggerNOT( const std::string& trigger, int runMin, int runMax ) { 

  TriggerMask newTriggerMask;
  newTriggerMask.HLTName = trigger;
  newTriggerMask.runMin = runMin;
  newTriggerMask.runMax = runMax;
  
  notTriggers_.push_back(newTriggerMask); 

}




void Ntp1Analyzer::LoadTrigger( int iEntry, bool verbose, TFile* condFile ) {

  
  TTree* treeCond = (condFile==0) ? 0 : (TTree*)(condFile->Get("Conditions"));

  std::vector<std::string> foundTriggers;
  std::vector<std::string> foundTriggersNOT;

  //new version: trigger loaded directly from ntp1 tree:
  if( treeCond==0 ) { 

    fChain->GetEntry(iEntry);
    //fChain->GetEntry(0);

    // required triggers:
    std::vector<int> index_requiredTriggers_found;
    for (unsigned int iTriggerMask=0; iTriggerMask<requiredTriggers_.size(); ++iTriggerMask)
      {
        if( runNumber<requiredTriggers_[iTriggerMask].runMin ) continue;
        if( requiredTriggers_[iTriggerMask].runMax>0 && runNumber>requiredTriggers_[iTriggerMask].runMax ) continue; 
        std::string thisTrigger = requiredTriggers_[iTriggerMask].HLTName;
        bool foundThisTrigger = false;
//std::cout << std::cout << "looking for: " << thisTrigger << std::endl;
        for(unsigned int i=0; i<nameHLT->size(); i++) 
          {
//std::cout << std::endl << indexHLT[i] << " " << nameHLT->at(i);
            TString nameHLT_tstr(nameHLT->at(i));
            if( nameHLT_tstr.BeginsWith((thisTrigger)) )
              {
//std::cout << " <----- HERE IT IS!" << std::endl;
                foundThisTrigger = true;
                index_requiredTriggers_found.push_back( indexHLT[i] ) ;
                foundTriggers.push_back( nameHLT->at(i) ) ;
                break;
              }
          }
          if( !foundThisTrigger && verbose ) std::cout << "-> WARNING!! Didn't find HLT path: " << thisTrigger << ". Ignoring it." << std::endl;
      }
    index_requiredTriggers_ = index_requiredTriggers_found;

    // NOT triggers
    std::vector<int> index_requiredTriggersNOT_found;
    for (unsigned int iTriggerMask=0; iTriggerMask<notTriggers_.size(); ++iTriggerMask)
      {
        if( runNumber<notTriggers_[iTriggerMask].runMin ) continue;
        if( notTriggers_[iTriggerMask].runMax>0 && runNumber>notTriggers_[iTriggerMask].runMax ) continue; 
        std::string thisTrigger = notTriggers_[iTriggerMask].HLTName;
        bool foundThisTrigger = false;
        for(unsigned int i=0; i<nameHLT->size(); i++) 
          {
//std::cout << std::endl << nameHLT->at(i);
            TString nameHLT_tstr(nameHLT->at(i));
            if( nameHLT_tstr.BeginsWith(thisTrigger) )
            //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) ) 
              {
//std::cout << " <----- HERE IT IS!" << std::endl;
                foundThisTrigger = true;
                index_requiredTriggersNOT_found.push_back( indexHLT[i] ) ;
                foundTriggersNOT.push_back( nameHLT->at(i) ) ;
                break;
              }
          }
          if( !foundThisTrigger && verbose ) std::cout << "-> WARNING!! Didn't find HLT path: " << thisTrigger << ". Ignoring it." << std::endl;
      }
    index_notTriggers_ = index_requiredTriggersNOT_found;

  }

  // old version now deprecated:
//} else { //old version: Conditions

//  int           nHLT_;
//  std::vector<std::string>  *nameHLT_;
//  std::vector<unsigned int> *indexHLT_;

//  //To get the pointers for the vectors
//  nameHLT_=0;
//  indexHLT_=0;

//  treeCond->SetBranchAddress("nHLT", &nHLT_);
//  treeCond->SetBranchAddress("nameHLT", &nameHLT_);
//  treeCond->SetBranchAddress("indexHLT", &indexHLT_);
//  treeCond->GetEntry(0);

//  std::vector<int> triggerMask;
//  for (std::vector< std::string >::const_iterator fIter=requiredTriggers_.begin();fIter!=requiredTriggers_.end();++fIter)
//    {
//      bool foundThisTrigger = false;
//      for(unsigned int i=0; i<nameHLT_->size(); i++) 
//        {
//          TString nameHLT_tstr(nameHLT_->at(i));
//          //if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
//          if( nameHLT_tstr.BeginsWith((*fIter).c_str()) ) 
//            {
//              foundThisTrigger = true;
//              triggerMask.push_back( indexHLT_->at(i) ) ;
//              break;
//            }
//        }
//        if( !foundThisTrigger && verbose ) std::cout << "-> WARNING!! Didn't find HLT path: " << (*fIter).c_str() << ". Ignoring it." << std::endl;
//    }
//  index_requiredTriggers_ = triggerMask;

//}


  if( requiredTriggers_.size()==0 && notTriggers_.size()==0 && verbose )
    std::cout << "-> No trigger selection required." << std::endl;

  for (int i=0;i<index_requiredTriggers_.size();++i)
    if( verbose ) std::cout << "[ReloadTriggerMask]::Requiring bit " << index_requiredTriggers_[i] << " " << foundTriggers[i] << std::endl;

  for (int i=0;i<index_notTriggers_.size();++i)
    if( verbose ) std::cout << "[ReloadTriggerMask]::Vetoing bit " << index_notTriggers_[i] << " " << foundTriggersNOT[i] << std::endl;


} // LoadTrigger




bool Ntp1Analyzer::PassedHLT( int iEntry, const std::string& HLTName ) { //default is OR of all required triggers (HLTName=="")


  if ( index_requiredTriggers_.size()==0 && index_notTriggers_.size()==0 && HLTName=="" ) return true;


  bool rememberToReset = false;
  std::vector<int> index_requiredTriggers_tmp = index_requiredTriggers_;
  std::vector<TriggerMask> requiredTriggers_tmp = requiredTriggers_;
  if( HLTName!="" ) {
    index_requiredTriggers_.clear();
    requiredTriggers_.clear();
    TriggerMask newTriggerMask;
    newTriggerMask.HLTName = HLTName;
    newTriggerMask.runMin = -1;
    newTriggerMask.runMax = -1;
    requiredTriggers_.push_back(newTriggerMask);
    this->LoadTrigger(iEntry, false);
    rememberToReset = true;
  }



  // first NOT triggers:
  for( int i=0; i<index_notTriggers_.size(); i++ ) {

    int block_veto =  index_notTriggers_[i]/30;
    int pos_veto = index_notTriggers_[i]%30;
    int word_veto = firedTrg[block_veto];
    
    if( (word_veto >> pos_veto)%2 ) {
      if( rememberToReset ) {
        index_requiredTriggers_ = index_requiredTriggers_tmp;
        requiredTriggers_ = requiredTriggers_tmp;
      }
      return false;
    }

  } // for not triggers
  

  // now required triggers:
  for( int i=0; i<index_requiredTriggers_.size(); i++ ) {

    if( HLTName=="" || requiredTriggers_[i].HLTName==HLTName ) {

      int block_required =  index_requiredTriggers_[i]/30;
      int pos_required = index_requiredTriggers_[i]%30;
      int word_required = firedTrg[block_required];

      if( (word_required >> pos_required)%2 ) {
        if( rememberToReset ) {
          index_requiredTriggers_ = index_requiredTriggers_tmp;
          requiredTriggers_ = requiredTriggers_tmp;
        }
        return true;
      }

    } // if required

  } // required trigger loop


  if( rememberToReset ) {
    index_requiredTriggers_ = index_requiredTriggers_tmp;
    requiredTriggers_ = requiredTriggers_tmp;
  }


  return false;

}



bool Ntp1Analyzer::isMatchedToHLT( float eta, float phi, float deltaR_max ) {


  bool match=false;
  for( int i=0; i<index_requiredTriggers_.size(); i++ ) {  // loop over require trigger paths
    
    int pathIndex=index_requiredTriggers_[i];
    // std::cout << "testing trigger " << pathIndex << " with " << sizePassing[pathIndex] << " passing objects" << std::endl; 
    
    if( sizePassing[pathIndex]>  0 ) {  //some object has passed the required trigger 
      
      for(int np = 0; np < sizePassing[pathIndex]; np++ ){
        int iP = indexPassing[ indexPassingPerPath[pathIndex] +np];
        // std::cout << "passing object eta: " << triggerObsEta[iP] << " phi: " <<  triggerObsPhi[iP] << std::endl; 

        TLorentzVector recoObject;
        recoObject.SetPtEtaPhiE( 10., eta, phi, 10. );
        TLorentzVector hltObject;
        hltObject.SetPtEtaPhiE( 10., triggerObsEta[iP], triggerObsPhi[iP], 10. );
        if( recoObject.DeltaR( hltObject ) < deltaR_max ) {
          match=true;
          //std::cout << "MATCH!" <<std::endl;	
          break;
        }
      }
    }
    if(match)  //it's enough if one path matches	
      break;
  }

  return match;

}




void Ntp1Analyzer::CreateOutputFile() {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova2ndLevel_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_2ndLevelTree_"+dataset_;
    else outfileName = analyzerType_ + "_2ndLevelTree";
   }


   if( flags_!="" )
     outfileName = outfileName + "_" + flags_;
   outfileName = outfileName + ".root";

   outfile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outfile_->cd();

   reducedTree_ = new TTree("reducedTree", "Reduced Tree");
   reducedTree_->SetMaxTreeSize(100000000000ULL); //setting max tree size to 100 GB

   h1_nCounter_ = new TH1F("nCounter", "", 1, 0., 1.);
   h1_nCounterPU_ = new TH1F("nCounterPU", "", 1, 0., 1.);
   h1_nCounterPU_ave_ = new TH1F("nCounterPU_ave", "", 1, 0., 1.);

}


Int_t Ntp1Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ntp1Analyzer::LoadTree(Long64_t entry)
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

void Ntp1Analyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   if (!tree) return;

   // Set object pointer
   nameHLT = 0;
   // Set branch addresses and branch pointers
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->GetEntry(0);
   isMC_ = (runNumber < 160000);

   h1_nCounter_->SetBinContent( 1, fChain->GetEntries() );



   std::string branchName;

   fChain->SetBranchAddress("nl1Technical", &nl1Technical, &b_nl1Technical);
   fChain->SetBranchAddress("l1Technical", l1Technical, &b_l1Technical);
   fChain->SetBranchAddress("nl1Global", &nl1Global, &b_nl1Global);
   fChain->SetBranchAddress("l1Global", l1Global, &b_l1Global);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("orbitNumber", &orbitNumber, &b_orbitNumber);
   fChain->SetBranchAddress("rhoFastjet", &rhoFastjet, &b_rhoFastjet);
   fChain->SetBranchAddress("rhoJetsFastJet", &rhoJetsFastjet, &b_rhoJetsFastjet); //CAPITAL J in fastJet - different in line above
   if( isMC_ ) {
     fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
     fChain->SetBranchAddress("nPU", nPU, &b_nPU);
     fChain->SetBranchAddress("bxPU", bxPU, &b_bxPU);
     fChain->SetBranchAddress("nMc", &nMc, &b_nMc);
     fChain->SetBranchAddress("pMc", pMc, &b_pMc);
     fChain->SetBranchAddress("thetaMc", thetaMc, &b_thetaMc);
     fChain->SetBranchAddress("etaMc", etaMc, &b_etaMc);
     fChain->SetBranchAddress("phiMc", phiMc, &b_phiMc);
     fChain->SetBranchAddress("energyMc", energyMc, &b_energyMc);
     fChain->SetBranchAddress("idMc", idMc, &b_idMc);
     fChain->SetBranchAddress("mothMc", mothMc, &b_mothMc);
     fChain->SetBranchAddress("statusMc", statusMc, &b_statusMc);
   }
   fChain->SetBranchAddress("nTrg", &nTrg, &b_nTrg);
   fChain->SetBranchAddress("firedTrg", firedTrg, &b_firedTrg);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("nameHLT", &nameHLT, &b_nameHLT);
   fChain->SetBranchAddress("indexHLT", indexHLT, &b_indexHLT);
   fChain->SetBranchAddress("nTriggerPaths", &nTriggerPaths, &b_nTriggerPaths);
   fChain->SetBranchAddress("nTriggerObsPassing", &nTriggerObsPassing, &b_nTriggerObsPassing);
   fChain->SetBranchAddress("sizePassing", sizePassing, &b_sizePassing);
   fChain->SetBranchAddress("indexPassing", indexPassing, &b_indexPassing);
   fChain->SetBranchAddress("indexPassingPerPath", indexPassingPerPath, &b_indexPassingPerPath);
   fChain->SetBranchAddress("nTriggerObs", &nTriggerObs, &b_nTriggerObs);
   fChain->SetBranchAddress("triggerObsPt", triggerObsPt, &b_triggerObsPt);
   fChain->SetBranchAddress("triggerObsPhi", triggerObsPhi, &b_triggerObsPhi);
   fChain->SetBranchAddress("triggerObsEta", triggerObsEta, &b_triggerObsEta);
   fChain->SetBranchAddress("triggerObsMass", triggerObsMass, &b_triggerObsMass);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("chargeEle", chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("energyEle", energyEle, &b_energyEle);
   fChain->SetBranchAddress("thetaEle", thetaEle, &b_thetaEle);
   fChain->SetBranchAddress("etaEle", etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", phiEle, &b_phiEle);
   fChain->SetBranchAddress("pxEle", pxEle, &b_pxEle);
   fChain->SetBranchAddress("pyEle", pyEle, &b_pyEle);
   fChain->SetBranchAddress("pzEle", pzEle, &b_pzEle);
   fChain->SetBranchAddress("vertexXEle", vertexXEle, &b_vertexXEle);
   fChain->SetBranchAddress("vertexYEle", vertexYEle, &b_vertexYEle);
   fChain->SetBranchAddress("vertexZEle", vertexZEle, &b_vertexZEle);
   fChain->SetBranchAddress("fiducialFlagsEle", fiducialFlagsEle, &b_fiducialFlagsEle);
   fChain->SetBranchAddress("recoFlagsEle", recoFlagsEle, &b_recoFlagsEle);
   fChain->SetBranchAddress("energyCorrectionsEle", energyCorrectionsEle, &b_energyCorrectionsEle);
   //fChain->SetBranchAddress("esEnergyEle", esEnergyEle, &b_esEnergyEle);
   fChain->SetBranchAddress("superClusterIndexEle", superClusterIndexEle, &b_superClusterIndexEle);
   fChain->SetBranchAddress("PFsuperClusterIndexEle", PFsuperClusterIndexEle, &b_PFsuperClusterIndexEle);
   fChain->SetBranchAddress("trackIndexEle", trackIndexEle, &b_trackIndexEle);
   fChain->SetBranchAddress("gsfTrackIndexEle", gsfTrackIndexEle, &b_gsfTrackIndexEle);
   fChain->SetBranchAddress("convDistEle", convDistEle, &b_convDistEle);
   fChain->SetBranchAddress("convDcotEle", convDcotEle, &b_convDcotEle);
   fChain->SetBranchAddress("convRadiusEle", convRadiusEle, &b_convRadiusEle);
   fChain->SetBranchAddress("convTrackIndexEle", convTrackIndexEle, &b_convTrackIndexEle);
   fChain->SetBranchAddress("hasMatchedConversionEle", hasMatchedConversionEle, &b_hasMatchedConversionEle);
   //fChain->SetBranchAddress("convXEle", convXEle, &b_convXEle);
   //fChain->SetBranchAddress("convYEle", convYEle, &b_convYEle);
   //fChain->SetBranchAddress("convZEle", convZEle, &b_convZEle);
   //fChain->SetBranchAddress("convChi2ProbEle", convChi2ProbEle, &b_convChi2ProbEle);
   fChain->SetBranchAddress("scPixChargeEle", scPixChargeEle, &b_scPixChargeEle);
   fChain->SetBranchAddress("classificationEle", classificationEle, &b_classificationEle);
   fChain->SetBranchAddress("standardClassificationEle", standardClassificationEle, &b_standardClassificationEle);
   fChain->SetBranchAddress("fbremEle", fbremEle, &b_fbremEle);
   fChain->SetBranchAddress("nbremsEle", nbremsEle, &b_nbremsEle);
   fChain->SetBranchAddress("hOverEEle", hOverEEle, &b_hOverEEle);
   fChain->SetBranchAddress("eSuperClusterOverPEle", eSuperClusterOverPEle, &b_eSuperClusterOverPEle);
   fChain->SetBranchAddress("eSeedOverPoutEle", eSeedOverPoutEle, &b_eSeedOverPoutEle);
   fChain->SetBranchAddress("deltaEtaAtVtxEle", deltaEtaAtVtxEle, &b_deltaEtaAtVtxEle);
   fChain->SetBranchAddress("deltaPhiAtVtxEle", deltaPhiAtVtxEle, &b_deltaPhiAtVtxEle);
   fChain->SetBranchAddress("deltaEtaAtCaloEle", deltaEtaAtCaloEle, &b_deltaEtaAtCaloEle);
   fChain->SetBranchAddress("deltaPhiAtCaloEle", deltaPhiAtCaloEle, &b_deltaPhiAtCaloEle);
   //fChain->SetBranchAddress("tipEle", tipEle, &b_tipEle);
   fChain->SetBranchAddress("dr03TkSumPtEle", dr03TkSumPtEle, &b_dr03TkSumPtEle);
   fChain->SetBranchAddress("dr03EcalRecHitSumEtEle", dr03EcalRecHitSumEtEle, &b_dr03EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr03HcalTowerSumEtEle", dr03HcalTowerSumEtEle, &b_dr03HcalTowerSumEtEle);
   fChain->SetBranchAddress("dr04TkSumPtEle", dr04TkSumPtEle, &b_dr04TkSumPtEle);
   fChain->SetBranchAddress("dr04EcalRecHitSumEtEle", dr04EcalRecHitSumEtEle, &b_dr04EcalRecHitSumEtEle);
   fChain->SetBranchAddress("dr04HcalTowerSumEtEle", dr04HcalTowerSumEtEle, &b_dr04HcalTowerSumEtEle);
   fChain->SetBranchAddress("scBasedEcalSum03Ele", scBasedEcalSum03Ele, &b_scBasedEcalSum03Ele);
   fChain->SetBranchAddress("scBasedEcalSum04Ele", scBasedEcalSum04Ele, &b_scBasedEcalSum04Ele);
   fChain->SetBranchAddress("pfCandChargedIso04Ele", pfCandChargedIso04Ele, &b_pfCandChargedIso04Ele);
   fChain->SetBranchAddress("pfCandNeutralIso04Ele", pfCandNeutralIso04Ele, &b_pfCandNeutralIso04Ele);
   fChain->SetBranchAddress("pfCandPhotonIso04Ele", pfCandPhotonIso04Ele, &b_pfCandPhotonIso04Ele);
   fChain->SetBranchAddress("eleIdLikelihoodEle", eleIdLikelihoodEle, &b_eleIdLikelihoodEle);
   fChain->SetBranchAddress("pflowMVAEle", pflowMVAEle, &b_pflowMVAEle);
   fChain->SetBranchAddress("mvaidtrigEle", mvaidtrigEle, &b_mvaidtrigEle);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("nPerPVPho", &nPerPVPho, &b_nPerPVPho);
   fChain->SetBranchAddress("chargePho", chargePho, &b_chargePho);
   fChain->SetBranchAddress("energyPho", energyPho, &b_energyPho);
   fChain->SetBranchAddress("thetaPho", thetaPho, &b_thetaPho);
   fChain->SetBranchAddress("etaPho", etaPho, &b_etaPho);
   fChain->SetBranchAddress("phiPho", phiPho, &b_phiPho);
   fChain->SetBranchAddress("pxPho", pxPho, &b_pxPho);
   fChain->SetBranchAddress("pyPho", pyPho, &b_pyPho);
   fChain->SetBranchAddress("pzPho", pzPho, &b_pzPho);
   fChain->SetBranchAddress("vertexXPho", vertexXPho, &b_vertexXPho);
   fChain->SetBranchAddress("vertexYPho", vertexYPho, &b_vertexYPho);
   fChain->SetBranchAddress("vertexZPho", vertexZPho, &b_vertexZPho);
   fChain->SetBranchAddress("fiducialFlagsPho", fiducialFlagsPho, &b_fiducialFlagsPho);
   fChain->SetBranchAddress("recoFlagsPho", recoFlagsPho, &b_recoFlagsPho);
   fChain->SetBranchAddress("superClusterIndexPho", superClusterIndexPho, &b_superClusterIndexPho);
   fChain->SetBranchAddress("PFsuperClusterIndexPho", PFsuperClusterIndexPho, &b_PFsuperClusterIndexPho);
   fChain->SetBranchAddress("hOverEPho", hOverEPho, &b_hOverEPho);
   fChain->SetBranchAddress("hTowOverEPho", hTowOverEPho, &b_hTowOverEPho);
   fChain->SetBranchAddress("dr03TkSumPtPho", dr03TkSumPtPho, &b_dr03TkSumPtPho);
   fChain->SetBranchAddress("dr03HollowTkSumPtPho", dr03HollowTkSumPtPho, &b_dr03HollowTkSumPtPho);
   fChain->SetBranchAddress("dr03EcalRecHitSumEtPho", dr03EcalRecHitSumEtPho, &b_dr03EcalRecHitSumEtPho);
   fChain->SetBranchAddress("dr03HcalTowerSumEtPho", dr03HcalTowerSumEtPho, &b_dr03HcalTowerSumEtPho);
   fChain->SetBranchAddress("dr04TkSumPtPho", dr04TkSumPtPho, &b_dr04TkSumPtPho);
   fChain->SetBranchAddress("dr04HollowTkSumPtPho", dr04HollowTkSumPtPho, &b_dr04HollowTkSumPtPho);
   fChain->SetBranchAddress("dr04EcalRecHitSumEtPho", dr04EcalRecHitSumEtPho, &b_dr04EcalRecHitSumEtPho);
   fChain->SetBranchAddress("dr04HcalTowerSumEtPho", dr04HcalTowerSumEtPho, &b_dr04HcalTowerSumEtPho);
   fChain->SetBranchAddress("chargedHadronIsoPho", chargedHadronIsoPho, &b_chargedHadronIsoPho);
   fChain->SetBranchAddress("neutralHadronIsoPho", neutralHadronIsoPho, &b_neutralHadronIsoPho);
   fChain->SetBranchAddress("photonIsoPho", photonIsoPho, &b_photonIsoPho);
   fChain->SetBranchAddress("hasPixelSeedPho", hasPixelSeedPho, &b_hasPixelSeedPho);
   fChain->SetBranchAddress("hasMatchedConversionPho", hasMatchedConversionPho, &b_hasMatchedConversionPho);
   fChain->SetBranchAddress("dr01NeutralHadronPFIsoPho", dr01NeutralHadronPFIsoPho, &b_dr01NeutralHadronPFIsoPho);
   fChain->SetBranchAddress("dr01PhotonPFIsoPho", dr01PhotonPFIsoPho, &b_dr01PhotonPFIsoPho);
   fChain->SetBranchAddress("dr02NeutralHadronPFIsoPho", dr02NeutralHadronPFIsoPho, &b_dr02NeutralHadronPFIsoPho);
   fChain->SetBranchAddress("dr02PhotonPFIsoPho", dr02PhotonPFIsoPho, &b_dr02PhotonPFIsoPho);
   fChain->SetBranchAddress("dr03NeutralHadronPFIsoPho", dr03NeutralHadronPFIsoPho, &b_dr03NeutralHadronPFIsoPho);
   fChain->SetBranchAddress("dr03PhotonPFIsoPho", dr03PhotonPFIsoPho, &b_dr03PhotonPFIsoPho);
   fChain->SetBranchAddress("dr04NeutralHadronPFIsoPho", dr04NeutralHadronPFIsoPho, &b_dr04NeutralHadronPFIsoPho);
   fChain->SetBranchAddress("dr04PhotonPFIsoPho", dr04PhotonPFIsoPho, &b_dr04PhotonPFIsoPho);
   fChain->SetBranchAddress("dr05NeutralHadronPFIsoPho", dr05NeutralHadronPFIsoPho, &b_dr05NeutralHadronPFIsoPho);
   fChain->SetBranchAddress("dr05PhotonPFIsoPho", dr05PhotonPFIsoPho, &b_dr05PhotonPFIsoPho);
   fChain->SetBranchAddress("dr06NeutralHadronPFIsoPho", dr06NeutralHadronPFIsoPho, &b_dr06NeutralHadronPFIsoPho);
   fChain->SetBranchAddress("dr06PhotonPFIsoPho", dr06PhotonPFIsoPho, &b_dr06PhotonPFIsoPho);
   fChain->SetBranchAddress("dr01ChargedHadronPFIsoPho", dr01ChargedHadronPFIsoPho, &b_dr01ChargedHadronPFIsoPho);
   fChain->SetBranchAddress("dr02ChargedHadronPFIsoPho", dr02ChargedHadronPFIsoPho, &b_dr02ChargedHadronPFIsoPho);
   fChain->SetBranchAddress("dr03ChargedHadronPFIsoPho", dr03ChargedHadronPFIsoPho, &b_dr03ChargedHadronPFIsoPho);
   fChain->SetBranchAddress("dr04ChargedHadronPFIsoPho", dr04ChargedHadronPFIsoPho, &b_dr04ChargedHadronPFIsoPho);
   fChain->SetBranchAddress("dr05ChargedHadronPFIsoPho", dr05ChargedHadronPFIsoPho, &b_dr05ChargedHadronPFIsoPho);
   fChain->SetBranchAddress("dr06ChargedHadronPFIsoPho", dr06ChargedHadronPFIsoPho, &b_dr06ChargedHadronPFIsoPho);
   fChain->SetBranchAddress("nConv", &nConv, &b_nConv);
   fChain->SetBranchAddress("pxPairConv", pxPairConv, &b_pxPairConv);
   fChain->SetBranchAddress("pyPairConv", pyPairConv, &b_pyPairConv);
   fChain->SetBranchAddress("pzPairConv", pzPairConv, &b_pzPairConv);
   fChain->SetBranchAddress("pxRefittedPairConv", pxRefittedPairConv, &b_pxRefittedPairConv);
   fChain->SetBranchAddress("pyRefittedPairConv", pyRefittedPairConv, &b_pyRefittedPairConv);
   fChain->SetBranchAddress("pzRefittedPairConv", pzRefittedPairConv, &b_pzRefittedPairConv);
   fChain->SetBranchAddress("etaRefittedPairConv", etaRefittedPairConv, &b_etaRefittedPairConv);
   fChain->SetBranchAddress("phiRefittedPairConv", phiRefittedPairConv, &b_phiRefittedPairConv);
   fChain->SetBranchAddress("ptRefittedPairConv", ptRefittedPairConv, &b_ptRefittedPairConv);
   fChain->SetBranchAddress("energyRefittedPairConv", energyRefittedPairConv, &b_energyRefittedPairConv);
   fChain->SetBranchAddress("eOverPRefittedPairConv", eOverPRefittedPairConv, &b_eOverPRefittedPairConv);
   fChain->SetBranchAddress("zOfPVFromTracksConv", zOfPVFromTracksConv, &b_zOfPVFromTracksConv);
   fChain->SetBranchAddress("xVtxConv", xVtxConv, &b_xVtxConv);
   fChain->SetBranchAddress("yVtxConv", yVtxConv, &b_yVtxConv);
   fChain->SetBranchAddress("zVtxConv", zVtxConv, &b_zVtxConv);
   fChain->SetBranchAddress("chi2VtxConv", chi2VtxConv, &b_chi2VtxConv);
   fChain->SetBranchAddress("chi2ProbVtxConv", chi2ProbVtxConv, &b_chi2ProbVtxConv);
   fChain->SetBranchAddress("isValidVtxConv", isValidVtxConv, &b_isValidVtxConv);
   fChain->SetBranchAddress("nTracksVtxConv", nTracksVtxConv, &b_nTracksVtxConv);
   fChain->SetBranchAddress("mvaOutVtxConv", mvaOutVtxConv, &b_mvaOutVtxConv);
   fChain->SetBranchAddress("trk1DzConv", trk1DzConv, &b_trk1DzConv);
   fChain->SetBranchAddress("trk1DzErrorConv", trk1DzErrorConv, &b_trk1DzErrorConv);
   fChain->SetBranchAddress("trk1ChargeConv", trk1ChargeConv, &b_trk1ChargeConv);
   fChain->SetBranchAddress("trk1AlgoConv", trk1AlgoConv, &b_trk1AlgoConv);
   fChain->SetBranchAddress("trk1PtConv", trk1PtConv, &b_trk1PtConv);
   fChain->SetBranchAddress("trk1D0Conv", trk1D0Conv, &b_trk1D0Conv);
   fChain->SetBranchAddress("trk1PoutConv", trk1PoutConv, &b_trk1PoutConv);
   fChain->SetBranchAddress("trk1PinConv", trk1PinConv, &b_trk1PinConv);
   fChain->SetBranchAddress("trk2DzConv", trk2DzConv, &b_trk2DzConv);
   fChain->SetBranchAddress("trk2DzErrorConv", trk2DzErrorConv, &b_trk2DzErrorConv);
   fChain->SetBranchAddress("trk2ChargeConv", trk2ChargeConv, &b_trk2ChargeConv);
   fChain->SetBranchAddress("trk2AlgoConv", trk2AlgoConv, &b_trk2AlgoConv);
   fChain->SetBranchAddress("trk2PtConv", trk2PtConv, &b_trk2PtConv);
   fChain->SetBranchAddress("trk2D0Conv", trk2D0Conv, &b_trk2D0Conv);
   fChain->SetBranchAddress("trk2PoutConv", trk2PoutConv, &b_trk2PoutConv);
   fChain->SetBranchAddress("trk2PinConv", trk2PinConv, &b_trk2PinConv);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("nBCSC", nBCSC, &b_nBCSC);
   fChain->SetBranchAddress("nCrystalsSC", nCrystalsSC, &b_nCrystalsSC);
   fChain->SetBranchAddress("rawEnergySC", rawEnergySC, &b_rawEnergySC);
   fChain->SetBranchAddress("energySC", energySC, &b_energySC);
   fChain->SetBranchAddress("etaSC", etaSC, &b_etaSC);
   fChain->SetBranchAddress("thetaSC", thetaSC, &b_thetaSC);
   fChain->SetBranchAddress("phiSC", phiSC, &b_phiSC);
   fChain->SetBranchAddress("phiWidthSC", phiWidthSC, &b_phiWidthSC);
   fChain->SetBranchAddress("etaWidthSC", etaWidthSC, &b_etaWidthSC);
   fChain->SetBranchAddress("e3x3SC", e3x3SC, &b_e3x3SC);
   fChain->SetBranchAddress("e5x5SC", e5x5SC, &b_e5x5SC);
   fChain->SetBranchAddress("eMaxSC", eMaxSC, &b_eMaxSC);
   fChain->SetBranchAddress("e2x2SC", e2x2SC, &b_e2x2SC);
   fChain->SetBranchAddress("e2ndSC", e2ndSC, &b_e2ndSC);
   fChain->SetBranchAddress("e1x5SC", e1x5SC, &b_e1x5SC);
   fChain->SetBranchAddress("e2x5MaxSC", e2x5MaxSC, &b_e2x5MaxSC);
   fChain->SetBranchAddress("e4SwissCrossSC", e4SwissCrossSC, &b_e4SwissCrossSC);
   fChain->SetBranchAddress("covIEtaIEtaSC", covIEtaIEtaSC, &b_covIEtaIEtaSC);
   fChain->SetBranchAddress("covIEtaIPhiSC", covIEtaIPhiSC, &b_covIEtaIPhiSC);
   fChain->SetBranchAddress("covIPhiIPhiSC", covIPhiIPhiSC, &b_covIPhiIPhiSC);
   fChain->SetBranchAddress("hOverESC", hOverESC, &b_hOverESC);
   fChain->SetBranchAddress("recoFlagSC", recoFlagSC, &b_recoFlagSC);
   fChain->SetBranchAddress("timeSC", timeSC, &b_timeSC);
   fChain->SetBranchAddress("chi2SC", chi2SC, &b_chi2SC);
   fChain->SetBranchAddress("seedEnergySC", seedEnergySC, &b_seedEnergySC);
   fChain->SetBranchAddress("sMajSC", sMajSC, &b_sMajSC);
   fChain->SetBranchAddress("sMinSC", sMinSC, &b_sMinSC);
   fChain->SetBranchAddress("alphaSC", alphaSC, &b_alphaSC);
   fChain->SetBranchAddress("nPFSC", &nPFSC, &b_nPFSC);
   fChain->SetBranchAddress("nBCPFSC", nBCPFSC, &b_nBCPFSC);
   fChain->SetBranchAddress("nCrystalsPFSC", nCrystalsPFSC, &b_nCrystalsPFSC);
   fChain->SetBranchAddress("rawEnergyPFSC", rawEnergyPFSC, &b_rawEnergyPFSC);
   fChain->SetBranchAddress("energyPFSC", energyPFSC, &b_energyPFSC);
   fChain->SetBranchAddress("etaPFSC", etaPFSC, &b_etaPFSC);
   fChain->SetBranchAddress("thetaPFSC", thetaPFSC, &b_thetaPFSC);
   fChain->SetBranchAddress("phiPFSC", phiPFSC, &b_phiPFSC);
   fChain->SetBranchAddress("phiWidthPFSC", phiWidthPFSC, &b_phiWidthPFSC);
   fChain->SetBranchAddress("etaWidthPFSC", etaWidthPFSC, &b_etaWidthPFSC);
   fChain->SetBranchAddress("e3x3PFSC", e3x3PFSC, &b_e3x3PFSC);
   fChain->SetBranchAddress("e5x5PFSC", e5x5PFSC, &b_e5x5PFSC);
   fChain->SetBranchAddress("eMaxPFSC", eMaxPFSC, &b_eMaxPFSC);
   fChain->SetBranchAddress("e2x2PFSC", e2x2PFSC, &b_e2x2PFSC);
   fChain->SetBranchAddress("e2ndPFSC", e2ndPFSC, &b_e2ndPFSC);
   fChain->SetBranchAddress("e1x5PFSC", e1x5PFSC, &b_e1x5PFSC);
   fChain->SetBranchAddress("e2x5MaxPFSC", e2x5MaxPFSC, &b_e2x5MaxPFSC);
   fChain->SetBranchAddress("e4SwissCrossPFSC", e4SwissCrossPFSC, &b_e4SwissCrossPFSC);
   fChain->SetBranchAddress("covIEtaIEtaPFSC", covIEtaIEtaPFSC, &b_covIEtaIEtaPFSC);
   fChain->SetBranchAddress("covIEtaIPhiPFSC", covIEtaIPhiPFSC, &b_covIEtaIPhiPFSC);
   fChain->SetBranchAddress("covIPhiIPhiPFSC", covIPhiIPhiPFSC, &b_covIPhiIPhiPFSC);
   fChain->SetBranchAddress("hOverEPFSC", hOverEPFSC, &b_hOverEPFSC);
   fChain->SetBranchAddress("recoFlagPFSC", recoFlagPFSC, &b_recoFlagPFSC);
   fChain->SetBranchAddress("timePFSC", timePFSC, &b_timePFSC);
   fChain->SetBranchAddress("chi2PFSC", chi2PFSC, &b_chi2PFSC);
   fChain->SetBranchAddress("seedEnergyPFSC", seedEnergyPFSC, &b_seedEnergyPFSC);
   //fChain->SetBranchAddress("nBC", &nBC, &b_nBC);
   //fChain->SetBranchAddress("nCrystalsBC", nCrystalsBC, &b_nCrystalsBC);
   //fChain->SetBranchAddress("energyBC", energyBC, &b_energyBC);
   //fChain->SetBranchAddress("etaBC", etaBC, &b_etaBC);
   //fChain->SetBranchAddress("thetaBC", thetaBC, &b_thetaBC);
   //fChain->SetBranchAddress("phiBC", phiBC, &b_phiBC);
   //fChain->SetBranchAddress("e3x3BC", e3x3BC, &b_e3x3BC);
   //fChain->SetBranchAddress("e5x5BC", e5x5BC, &b_e5x5BC);
   //fChain->SetBranchAddress("eMaxBC", eMaxBC, &b_eMaxBC);
   //fChain->SetBranchAddress("e2x2BC", e2x2BC, &b_e2x2BC);
   //fChain->SetBranchAddress("e2ndBC", e2ndBC, &b_e2ndBC);
   //fChain->SetBranchAddress("covIEtaIEtaBC", covIEtaIEtaBC, &b_covIEtaIEtaBC);
   //fChain->SetBranchAddress("covIEtaIPhiBC", covIEtaIPhiBC, &b_covIEtaIPhiBC);
   //fChain->SetBranchAddress("covIPhiIPhiBC", covIPhiIPhiBC, &b_covIPhiIPhiBC);
   //fChain->SetBranchAddress("recoFlagBC", recoFlagBC, &b_recoFlagBC);
   //fChain->SetBranchAddress("timeBC", timeBC, &b_timeBC);
   //fChain->SetBranchAddress("chi2BC", chi2BC, &b_chi2BC);
   //fChain->SetBranchAddress("seedEnergyBC", seedEnergyBC, &b_seedEnergyBC);
   //fChain->SetBranchAddress("idClosProblBC", idClosProblBC, &b_idClosProblBC);
   //fChain->SetBranchAddress("sevClosProblBC", sevClosProblBC, &b_sevClosProblBC);
   //fChain->SetBranchAddress("fracClosProblBC", fracClosProblBC, &b_fracClosProblBC);
   //fChain->SetBranchAddress("indexSCBC", indexSCBC, &b_indexSCBC);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("pxTrack", pxTrack, &b_pxTrack);
   fChain->SetBranchAddress("pyTrack", pyTrack, &b_pyTrack);
   fChain->SetBranchAddress("pzTrack", pzTrack, &b_pzTrack);
   fChain->SetBranchAddress("vtxIndexTrack", vtxIndexTrack, &b_vtxIndexTrack);
   fChain->SetBranchAddress("vtxWeightTrack", vtxWeightTrack, &b_vtxWeightTrack);
   fChain->SetBranchAddress("chargeTrack", chargeTrack, &b_chargeTrack);
   fChain->SetBranchAddress("ptErrorTrack", ptErrorTrack, &b_ptErrorTrack);
   fChain->SetBranchAddress("trackValidHitsTrack", trackValidHitsTrack, &b_trackValidHitsTrack);
   fChain->SetBranchAddress("trackLostHitsTrack", trackLostHitsTrack, &b_trackLostHitsTrack);
   fChain->SetBranchAddress("trackNormalizedChi2Track", trackNormalizedChi2Track, &b_trackNormalizedChi2Track);
   fChain->SetBranchAddress("qualityMaskTrack", qualityMaskTrack, &b_qualityMaskTrack);
   fChain->SetBranchAddress("impactPar3DTrack", impactPar3DTrack, &b_impactPar3DTrack);
   fChain->SetBranchAddress("impactPar3DErrorTrack", impactPar3DErrorTrack, &b_impactPar3DErrorTrack);
   fChain->SetBranchAddress("transvImpactParTrack", transvImpactParTrack, &b_transvImpactParTrack);
   fChain->SetBranchAddress("transvImpactParErrorTrack", transvImpactParErrorTrack, &b_transvImpactParErrorTrack);
   fChain->SetBranchAddress("trackVxTrack", trackVxTrack, &b_trackVxTrack);
   fChain->SetBranchAddress("trackVyTrack", trackVyTrack, &b_trackVyTrack);
   fChain->SetBranchAddress("trackVzTrack", trackVzTrack, &b_trackVzTrack);
   //fChain->SetBranchAddress("pxAtOuterTrack", pxAtOuterTrack, &b_pxAtOuterTrack);
   //fChain->SetBranchAddress("pyAtOuterTrack", pyAtOuterTrack, &b_pyAtOuterTrack);
   //fChain->SetBranchAddress("pzAtOuterTrack", pzAtOuterTrack, &b_pzAtOuterTrack);
   //fChain->SetBranchAddress("xAtOuterTrack", xAtOuterTrack, &b_xAtOuterTrack);
   //fChain->SetBranchAddress("yAtOuterTrack", yAtOuterTrack, &b_yAtOuterTrack);
   //fChain->SetBranchAddress("zAtOuterTrack", zAtOuterTrack, &b_zAtOuterTrack);
   //fChain->SetBranchAddress("pxAtInnerTrack", pxAtInnerTrack, &b_pxAtInnerTrack);
   //fChain->SetBranchAddress("pyAtInnerTrack", pyAtInnerTrack, &b_pyAtInnerTrack);
   //fChain->SetBranchAddress("pzAtInnerTrack", pzAtInnerTrack, &b_pzAtInnerTrack);
   //fChain->SetBranchAddress("xAtInnerTrack", xAtInnerTrack, &b_xAtInnerTrack);
   //fChain->SetBranchAddress("yAtInnerTrack", yAtInnerTrack, &b_yAtInnerTrack);
   //fChain->SetBranchAddress("zAtInnerTrack", zAtInnerTrack, &b_zAtInnerTrack);
   //fChain->SetBranchAddress("recHitsSizeTrack", recHitsSizeTrack, &b_recHitsSizeTrack);
   fChain->SetBranchAddress("pixelHitsTrack", pixelHitsTrack, &b_pixelHitsTrack);
   fChain->SetBranchAddress("expInnerLayersTrack", expInnerLayersTrack, &b_expInnerLayersTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsTrack", numberOfValidPixelBarrelHitsTrack, &b_numberOfValidPixelBarrelHitsTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsTrack", numberOfValidPixelEndcapHitsTrack, &b_numberOfValidPixelEndcapHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsTrack", numberOfValidStripTIBHitsTrack, &b_numberOfValidStripTIBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsTrack", numberOfValidStripTIDHitsTrack, &b_numberOfValidStripTIDHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsTrack", numberOfValidStripTOBHitsTrack, &b_numberOfValidStripTOBHitsTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsTrack", numberOfValidStripTECHitsTrack, &b_numberOfValidStripTECHitsTrack);
   fChain->SetBranchAddress("nGsfTrack", &nGsfTrack, &b_nGsfTrack);
   fChain->SetBranchAddress("pxGsfTrack", pxGsfTrack, &b_pxGsfTrack);
   fChain->SetBranchAddress("pyGsfTrack", pyGsfTrack, &b_pyGsfTrack);
   fChain->SetBranchAddress("pzGsfTrack", pzGsfTrack, &b_pzGsfTrack);
   fChain->SetBranchAddress("vtxIndexGsfTrack", vtxIndexGsfTrack, &b_vtxIndexGsfTrack);
   fChain->SetBranchAddress("vtxWeightGsfTrack", vtxWeightGsfTrack, &b_vtxWeightGsfTrack);
   fChain->SetBranchAddress("chargeGsfTrack", chargeGsfTrack, &b_chargeGsfTrack);
   fChain->SetBranchAddress("ptErrorGsfTrack", ptErrorGsfTrack, &b_ptErrorGsfTrack);
   fChain->SetBranchAddress("trackValidHitsGsfTrack", trackValidHitsGsfTrack, &b_trackValidHitsGsfTrack);
   fChain->SetBranchAddress("trackLostHitsGsfTrack", trackLostHitsGsfTrack, &b_trackLostHitsGsfTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GsfTrack", trackNormalizedChi2GsfTrack, &b_trackNormalizedChi2GsfTrack);
   fChain->SetBranchAddress("qualityMaskGsfTrack", qualityMaskGsfTrack, &b_qualityMaskGsfTrack);
   fChain->SetBranchAddress("impactPar3DGsfTrack", impactPar3DGsfTrack, &b_impactPar3DGsfTrack);
   fChain->SetBranchAddress("impactPar3DErrorGsfTrack", impactPar3DErrorGsfTrack, &b_impactPar3DErrorGsfTrack);
   fChain->SetBranchAddress("transvImpactParGsfTrack", transvImpactParGsfTrack, &b_transvImpactParGsfTrack);
   fChain->SetBranchAddress("transvImpactParErrorGsfTrack", transvImpactParErrorGsfTrack, &b_transvImpactParErrorGsfTrack);
   fChain->SetBranchAddress("trackVxGsfTrack", trackVxGsfTrack, &b_trackVxGsfTrack);
   fChain->SetBranchAddress("trackVyGsfTrack", trackVyGsfTrack, &b_trackVyGsfTrack);
   fChain->SetBranchAddress("trackVzGsfTrack", trackVzGsfTrack, &b_trackVzGsfTrack);
   //fChain->SetBranchAddress("pxAtOuterGsfTrack", pxAtOuterGsfTrack, &b_pxAtOuterGsfTrack);
   //fChain->SetBranchAddress("pyAtOuterGsfTrack", pyAtOuterGsfTrack, &b_pyAtOuterGsfTrack);
   //fChain->SetBranchAddress("pzAtOuterGsfTrack", pzAtOuterGsfTrack, &b_pzAtOuterGsfTrack);
   //fChain->SetBranchAddress("xAtOuterGsfTrack", xAtOuterGsfTrack, &b_xAtOuterGsfTrack);
   //fChain->SetBranchAddress("yAtOuterGsfTrack", yAtOuterGsfTrack, &b_yAtOuterGsfTrack);
   //fChain->SetBranchAddress("zAtOuterGsfTrack", zAtOuterGsfTrack, &b_zAtOuterGsfTrack);
   //fChain->SetBranchAddress("pxAtInnerGsfTrack", pxAtInnerGsfTrack, &b_pxAtInnerGsfTrack);
   //fChain->SetBranchAddress("pyAtInnerGsfTrack", pyAtInnerGsfTrack, &b_pyAtInnerGsfTrack);
   //fChain->SetBranchAddress("pzAtInnerGsfTrack", pzAtInnerGsfTrack, &b_pzAtInnerGsfTrack);
   //fChain->SetBranchAddress("xAtInnerGsfTrack", xAtInnerGsfTrack, &b_xAtInnerGsfTrack);
   //fChain->SetBranchAddress("yAtInnerGsfTrack", yAtInnerGsfTrack, &b_yAtInnerGsfTrack);
   //fChain->SetBranchAddress("zAtInnerGsfTrack", zAtInnerGsfTrack, &b_zAtInnerGsfTrack);
   //fChain->SetBranchAddress("recHitsSizeGsfTrack", recHitsSizeGsfTrack, &b_recHitsSizeGsfTrack);
   fChain->SetBranchAddress("pixelHitsGsfTrack", pixelHitsGsfTrack, &b_pixelHitsGsfTrack);
   fChain->SetBranchAddress("expInnerLayersGsfTrack", expInnerLayersGsfTrack, &b_expInnerLayersGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGsfTrack", numberOfValidPixelBarrelHitsGsfTrack, &b_numberOfValidPixelBarrelHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGsfTrack", numberOfValidPixelEndcapHitsGsfTrack, &b_numberOfValidPixelEndcapHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGsfTrack", numberOfValidStripTIBHitsGsfTrack, &b_numberOfValidStripTIBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGsfTrack", numberOfValidStripTIDHitsGsfTrack, &b_numberOfValidStripTIDHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGsfTrack", numberOfValidStripTOBHitsGsfTrack, &b_numberOfValidStripTOBHitsGsfTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGsfTrack", numberOfValidStripTECHitsGsfTrack, &b_numberOfValidStripTECHitsGsfTrack);
   fChain->SetBranchAddress("chargeModeGsfTrack", chargeModeGsfTrack, &b_chargeModeGsfTrack);
   fChain->SetBranchAddress("pxModeGsfTrack", pxModeGsfTrack, &b_pxModeGsfTrack);
   fChain->SetBranchAddress("pyModeGsfTrack", pyModeGsfTrack, &b_pyModeGsfTrack);
   fChain->SetBranchAddress("pzModeGsfTrack", pzModeGsfTrack, &b_pzModeGsfTrack);
   //fChain->SetBranchAddress("recoFlagsGsfTrack", recoFlagsGsfTrack, &b_recoFlagsGsfTrack);
   fChain->SetBranchAddress("nGlobalMuonTrack", &nGlobalMuonTrack, &b_nGlobalMuonTrack);
   fChain->SetBranchAddress("pxGlobalMuonTrack", pxGlobalMuonTrack, &b_pxGlobalMuonTrack);
   fChain->SetBranchAddress("pyGlobalMuonTrack", pyGlobalMuonTrack, &b_pyGlobalMuonTrack);
   fChain->SetBranchAddress("pzGlobalMuonTrack", pzGlobalMuonTrack, &b_pzGlobalMuonTrack);
   fChain->SetBranchAddress("vtxIndexGlobalMuonTrack", vtxIndexGlobalMuonTrack, &b_vtxIndexGlobalMuonTrack);
   fChain->SetBranchAddress("vtxWeightGlobalMuonTrack", vtxWeightGlobalMuonTrack, &b_vtxWeightGlobalMuonTrack);
   fChain->SetBranchAddress("chargeGlobalMuonTrack", chargeGlobalMuonTrack, &b_chargeGlobalMuonTrack);
   fChain->SetBranchAddress("ptErrorGlobalMuonTrack", ptErrorGlobalMuonTrack, &b_ptErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackValidHitsGlobalMuonTrack", trackValidHitsGlobalMuonTrack, &b_trackValidHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackLostHitsGlobalMuonTrack", trackLostHitsGlobalMuonTrack, &b_trackLostHitsGlobalMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2GlobalMuonTrack", trackNormalizedChi2GlobalMuonTrack, &b_trackNormalizedChi2GlobalMuonTrack);
   fChain->SetBranchAddress("qualityMaskGlobalMuonTrack", qualityMaskGlobalMuonTrack, &b_qualityMaskGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DGlobalMuonTrack", impactPar3DGlobalMuonTrack, &b_impactPar3DGlobalMuonTrack);
   fChain->SetBranchAddress("impactPar3DErrorGlobalMuonTrack", impactPar3DErrorGlobalMuonTrack, &b_impactPar3DErrorGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParGlobalMuonTrack", transvImpactParGlobalMuonTrack, &b_transvImpactParGlobalMuonTrack);
   fChain->SetBranchAddress("transvImpactParErrorGlobalMuonTrack", transvImpactParErrorGlobalMuonTrack, &b_transvImpactParErrorGlobalMuonTrack);
   fChain->SetBranchAddress("trackVxGlobalMuonTrack", trackVxGlobalMuonTrack, &b_trackVxGlobalMuonTrack);
   fChain->SetBranchAddress("trackVyGlobalMuonTrack", trackVyGlobalMuonTrack, &b_trackVyGlobalMuonTrack);
   fChain->SetBranchAddress("trackVzGlobalMuonTrack", trackVzGlobalMuonTrack, &b_trackVzGlobalMuonTrack);
   //fChain->SetBranchAddress("pxAtOuterGlobalMuonTrack", pxAtOuterGlobalMuonTrack, &b_pxAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("pyAtOuterGlobalMuonTrack", pyAtOuterGlobalMuonTrack, &b_pyAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("pzAtOuterGlobalMuonTrack", pzAtOuterGlobalMuonTrack, &b_pzAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("xAtOuterGlobalMuonTrack", xAtOuterGlobalMuonTrack, &b_xAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("yAtOuterGlobalMuonTrack", yAtOuterGlobalMuonTrack, &b_yAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("zAtOuterGlobalMuonTrack", zAtOuterGlobalMuonTrack, &b_zAtOuterGlobalMuonTrack);
   //fChain->SetBranchAddress("pxAtInnerGlobalMuonTrack", pxAtInnerGlobalMuonTrack, &b_pxAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("pyAtInnerGlobalMuonTrack", pyAtInnerGlobalMuonTrack, &b_pyAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("pzAtInnerGlobalMuonTrack", pzAtInnerGlobalMuonTrack, &b_pzAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("xAtInnerGlobalMuonTrack", xAtInnerGlobalMuonTrack, &b_xAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("yAtInnerGlobalMuonTrack", yAtInnerGlobalMuonTrack, &b_yAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("zAtInnerGlobalMuonTrack", zAtInnerGlobalMuonTrack, &b_zAtInnerGlobalMuonTrack);
   //fChain->SetBranchAddress("recHitsSizeGlobalMuonTrack", recHitsSizeGlobalMuonTrack, &b_recHitsSizeGlobalMuonTrack);
   fChain->SetBranchAddress("pixelHitsGlobalMuonTrack", pixelHitsGlobalMuonTrack, &b_pixelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("expInnerLayersGlobalMuonTrack", expInnerLayersGlobalMuonTrack, &b_expInnerLayersGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsGlobalMuonTrack", numberOfValidPixelBarrelHitsGlobalMuonTrack, &b_numberOfValidPixelBarrelHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsGlobalMuonTrack", numberOfValidPixelEndcapHitsGlobalMuonTrack, &b_numberOfValidPixelEndcapHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsGlobalMuonTrack", numberOfValidStripTIBHitsGlobalMuonTrack, &b_numberOfValidStripTIBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsGlobalMuonTrack", numberOfValidStripTIDHitsGlobalMuonTrack, &b_numberOfValidStripTIDHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsGlobalMuonTrack", numberOfValidStripTOBHitsGlobalMuonTrack, &b_numberOfValidStripTOBHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsGlobalMuonTrack", numberOfValidStripTECHitsGlobalMuonTrack, &b_numberOfValidStripTECHitsGlobalMuonTrack);
   fChain->SetBranchAddress("numberOfValidMuonHitsGlobalMuonTrack", numberOfValidMuonHitsGlobalMuonTrack, &b_numberOfValidMuonHitsGlobalMuonTrack);
   fChain->SetBranchAddress("nSTAMuonTrack", &nSTAMuonTrack, &b_nSTAMuonTrack);
   fChain->SetBranchAddress("pxSTAMuonTrack", pxSTAMuonTrack, &b_pxSTAMuonTrack);
   fChain->SetBranchAddress("pySTAMuonTrack", pySTAMuonTrack, &b_pySTAMuonTrack);
   fChain->SetBranchAddress("pzSTAMuonTrack", pzSTAMuonTrack, &b_pzSTAMuonTrack);
   //fChain->SetBranchAddress("vtxIndexSTAMuonTrack", vtxIndexSTAMuonTrack, &b_vtxIndexSTAMuonTrack);
   //fChain->SetBranchAddress("vtxWeightSTAMuonTrack", vtxWeightSTAMuonTrack, &b_vtxWeightSTAMuonTrack);
   fChain->SetBranchAddress("chargeSTAMuonTrack", chargeSTAMuonTrack, &b_chargeSTAMuonTrack);
   fChain->SetBranchAddress("ptErrorSTAMuonTrack", ptErrorSTAMuonTrack, &b_ptErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackValidHitsSTAMuonTrack", trackValidHitsSTAMuonTrack, &b_trackValidHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackLostHitsSTAMuonTrack", trackLostHitsSTAMuonTrack, &b_trackLostHitsSTAMuonTrack);
   fChain->SetBranchAddress("trackNormalizedChi2STAMuonTrack", trackNormalizedChi2STAMuonTrack, &b_trackNormalizedChi2STAMuonTrack);
   fChain->SetBranchAddress("qualityMaskSTAMuonTrack", qualityMaskSTAMuonTrack, &b_qualityMaskSTAMuonTrack);
   //fChain->SetBranchAddress("impactPar3DSTAMuonTrack", impactPar3DSTAMuonTrack, &b_impactPar3DSTAMuonTrack);
   //fChain->SetBranchAddress("impactPar3DErrorSTAMuonTrack", impactPar3DErrorSTAMuonTrack, &b_impactPar3DErrorSTAMuonTrack);
   //fChain->SetBranchAddress("transvImpactParSTAMuonTrack", transvImpactParSTAMuonTrack, &b_transvImpactParSTAMuonTrack);
   //fChain->SetBranchAddress("transvImpactParErrorSTAMuonTrack", transvImpactParErrorSTAMuonTrack, &b_transvImpactParErrorSTAMuonTrack);
   fChain->SetBranchAddress("trackVxSTAMuonTrack", trackVxSTAMuonTrack, &b_trackVxSTAMuonTrack);
   fChain->SetBranchAddress("trackVySTAMuonTrack", trackVySTAMuonTrack, &b_trackVySTAMuonTrack);
   fChain->SetBranchAddress("trackVzSTAMuonTrack", trackVzSTAMuonTrack, &b_trackVzSTAMuonTrack);
   //fChain->SetBranchAddress("pxAtOuterSTAMuonTrack", pxAtOuterSTAMuonTrack, &b_pxAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("pyAtOuterSTAMuonTrack", pyAtOuterSTAMuonTrack, &b_pyAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("pzAtOuterSTAMuonTrack", pzAtOuterSTAMuonTrack, &b_pzAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("xAtOuterSTAMuonTrack", xAtOuterSTAMuonTrack, &b_xAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("yAtOuterSTAMuonTrack", yAtOuterSTAMuonTrack, &b_yAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("zAtOuterSTAMuonTrack", zAtOuterSTAMuonTrack, &b_zAtOuterSTAMuonTrack);
   //fChain->SetBranchAddress("pxAtInnerSTAMuonTrack", pxAtInnerSTAMuonTrack, &b_pxAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("pyAtInnerSTAMuonTrack", pyAtInnerSTAMuonTrack, &b_pyAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("pzAtInnerSTAMuonTrack", pzAtInnerSTAMuonTrack, &b_pzAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("xAtInnerSTAMuonTrack", xAtInnerSTAMuonTrack, &b_xAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("yAtInnerSTAMuonTrack", yAtInnerSTAMuonTrack, &b_yAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("zAtInnerSTAMuonTrack", zAtInnerSTAMuonTrack, &b_zAtInnerSTAMuonTrack);
   //fChain->SetBranchAddress("recHitsSizeSTAMuonTrack", recHitsSizeSTAMuonTrack, &b_recHitsSizeSTAMuonTrack);
   fChain->SetBranchAddress("pixelHitsSTAMuonTrack", pixelHitsSTAMuonTrack, &b_pixelHitsSTAMuonTrack);
   fChain->SetBranchAddress("expInnerLayersSTAMuonTrack", expInnerLayersSTAMuonTrack, &b_expInnerLayersSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelBarrelHitsSTAMuonTrack", numberOfValidPixelBarrelHitsSTAMuonTrack, &b_numberOfValidPixelBarrelHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidPixelEndcapHitsSTAMuonTrack", numberOfValidPixelEndcapHitsSTAMuonTrack, &b_numberOfValidPixelEndcapHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIBHitsSTAMuonTrack", numberOfValidStripTIBHitsSTAMuonTrack, &b_numberOfValidStripTIBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTIDHitsSTAMuonTrack", numberOfValidStripTIDHitsSTAMuonTrack, &b_numberOfValidStripTIDHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTOBHitsSTAMuonTrack", numberOfValidStripTOBHitsSTAMuonTrack, &b_numberOfValidStripTOBHitsSTAMuonTrack);
   fChain->SetBranchAddress("numberOfValidStripTECHitsSTAMuonTrack", numberOfValidStripTECHitsSTAMuonTrack, &b_numberOfValidStripTECHitsSTAMuonTrack);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVxPV", PVxPV, &b_PVxPV);
   fChain->SetBranchAddress("PVyPV", PVyPV, &b_PVyPV);
   fChain->SetBranchAddress("PVzPV", PVzPV, &b_PVzPV);
   fChain->SetBranchAddress("PVErrxPV", PVErrxPV, &b_PVErrxPV);
   fChain->SetBranchAddress("PVErryPV", PVErryPV, &b_PVErryPV);
   fChain->SetBranchAddress("PVErrzPV", PVErrzPV, &b_PVErrzPV);
   fChain->SetBranchAddress("SumPtPV", SumPtPV, &b_SumPtPV);
   fChain->SetBranchAddress("ndofPV", ndofPV, &b_ndofPV);
   fChain->SetBranchAddress("chi2PV", chi2PV, &b_chi2PV);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("chargeMuon", chargeMuon, &b_chargeMuon);
   fChain->SetBranchAddress("energyMuon", energyMuon, &b_energyMuon);
   fChain->SetBranchAddress("thetaMuon", thetaMuon, &b_thetaMuon);
   fChain->SetBranchAddress("etaMuon", etaMuon, &b_etaMuon);
   fChain->SetBranchAddress("phiMuon", phiMuon, &b_phiMuon);
   fChain->SetBranchAddress("pxMuon", pxMuon, &b_pxMuon);
   fChain->SetBranchAddress("pyMuon", pyMuon, &b_pyMuon);
   fChain->SetBranchAddress("pzMuon", pzMuon, &b_pzMuon);
   fChain->SetBranchAddress("vertexXMuon", vertexXMuon, &b_vertexXMuon);
   fChain->SetBranchAddress("vertexYMuon", vertexYMuon, &b_vertexYMuon);
   fChain->SetBranchAddress("vertexZMuon", vertexZMuon, &b_vertexZMuon);
   fChain->SetBranchAddress("trackIndexMuon", trackIndexMuon, &b_trackIndexMuon);
   fChain->SetBranchAddress("standAloneTrackIndexMuon", standAloneTrackIndexMuon, &b_standAloneTrackIndexMuon);
   fChain->SetBranchAddress("combinedTrackIndexMuon", combinedTrackIndexMuon, &b_combinedTrackIndexMuon);
   fChain->SetBranchAddress("muonIdMuon", muonIdMuon, &b_muonIdMuon);
   fChain->SetBranchAddress("sumPt03Muon", sumPt03Muon, &b_sumPt03Muon);
   fChain->SetBranchAddress("emEt03Muon", emEt03Muon, &b_emEt03Muon);
   fChain->SetBranchAddress("hadEt03Muon", hadEt03Muon, &b_hadEt03Muon);
   fChain->SetBranchAddress("hoEt03Muon", hoEt03Muon, &b_hoEt03Muon);
   fChain->SetBranchAddress("nTrk03Muon", nTrk03Muon, &b_nTrk03Muon);
   fChain->SetBranchAddress("nJets03Muon", nJets03Muon, &b_nJets03Muon);
   fChain->SetBranchAddress("sumPt05Muon", sumPt05Muon, &b_sumPt05Muon);
   fChain->SetBranchAddress("emEt05Muon", emEt05Muon, &b_emEt05Muon);
   fChain->SetBranchAddress("hadEt05Muon", hadEt05Muon, &b_hadEt05Muon);
   fChain->SetBranchAddress("hoEt05Muon", hoEt05Muon, &b_hoEt05Muon);
   fChain->SetBranchAddress("nTrk05Muon", nTrk05Muon, &b_nTrk05Muon);
   fChain->SetBranchAddress("nJets05Muon", nJets05Muon, &b_nJets05Muon);
   fChain->SetBranchAddress("EcalExpDepoMuon", EcalExpDepoMuon, &b_EcalExpDepoMuon);
   fChain->SetBranchAddress("HcalExpDepoMuon", HcalExpDepoMuon, &b_HcalExpDepoMuon);
   fChain->SetBranchAddress("HoExpDepoMuon", HoExpDepoMuon, &b_HoExpDepoMuon);
   fChain->SetBranchAddress("emS9Muon", emS9Muon, &b_emS9Muon);
   fChain->SetBranchAddress("hadS9Muon", hadS9Muon, &b_hadS9Muon);
   fChain->SetBranchAddress("hoS9Muon", hoS9Muon, &b_hoS9Muon);
   fChain->SetBranchAddress("CaloCompMuon", CaloCompMuon, &b_CaloCompMuon);
   fChain->SetBranchAddress("numberOfMatchesMuon", numberOfMatchesMuon, &b_numberOfMatchesMuon);
   fChain->SetBranchAddress("mvaisoMuon", mvaisoMuon, &b_mvaisoMuon);
   fChain->SetBranchAddress("pfmuonIdMuon", pfmuonIdMuon, &b_pfmuonIdMuon);
   fChain->SetBranchAddress("nMet", &nMet, &b_nMet);
   fChain->SetBranchAddress("chargeMet", chargeMet, &b_chargeMet);
   fChain->SetBranchAddress("energyMet", energyMet, &b_energyMet);
   fChain->SetBranchAddress("thetaMet", thetaMet, &b_thetaMet);
   fChain->SetBranchAddress("etaMet", etaMet, &b_etaMet);
   fChain->SetBranchAddress("phiMet", phiMet, &b_phiMet);
   fChain->SetBranchAddress("pxMet", pxMet, &b_pxMet);
   fChain->SetBranchAddress("pyMet", pyMet, &b_pyMet);
   fChain->SetBranchAddress("pzMet", pzMet, &b_pzMet);
   fChain->SetBranchAddress("vertexXMet", vertexXMet, &b_vertexXMet);
   fChain->SetBranchAddress("vertexYMet", vertexYMet, &b_vertexYMet);
   fChain->SetBranchAddress("vertexZMet", vertexZMet, &b_vertexZMet);
   fChain->SetBranchAddress("nTCMet", &nTCMet, &b_nTCMet);
   fChain->SetBranchAddress("chargeTCMet", chargeTCMet, &b_chargeTCMet);
   fChain->SetBranchAddress("energyTCMet", energyTCMet, &b_energyTCMet);
   fChain->SetBranchAddress("thetaTCMet", thetaTCMet, &b_thetaTCMet);
   fChain->SetBranchAddress("etaTCMet", etaTCMet, &b_etaTCMet);
   fChain->SetBranchAddress("phiTCMet", phiTCMet, &b_phiTCMet);
   fChain->SetBranchAddress("pxTCMet", pxTCMet, &b_pxTCMet);
   fChain->SetBranchAddress("pyTCMet", pyTCMet, &b_pyTCMet);
   fChain->SetBranchAddress("pzTCMet", pzTCMet, &b_pzTCMet);
   fChain->SetBranchAddress("vertexXTCMet", vertexXTCMet, &b_vertexXTCMet);
   fChain->SetBranchAddress("vertexYTCMet", vertexYTCMet, &b_vertexYTCMet);
   fChain->SetBranchAddress("vertexZTCMet", vertexZTCMet, &b_vertexZTCMet);
   fChain->SetBranchAddress("nPFMet", &nPFMet, &b_nPFMet);
   fChain->SetBranchAddress("chargePFMet", chargePFMet, &b_chargePFMet);
   fChain->SetBranchAddress("energyPFMet", energyPFMet, &b_energyPFMet);
   fChain->SetBranchAddress("thetaPFMet", thetaPFMet, &b_thetaPFMet);
   fChain->SetBranchAddress("etaPFMet", etaPFMet, &b_etaPFMet);
   fChain->SetBranchAddress("phiPFMet", phiPFMet, &b_phiPFMet);
   fChain->SetBranchAddress("pxPFMet", pxPFMet, &b_pxPFMet);
   fChain->SetBranchAddress("pyPFMet", pyPFMet, &b_pyPFMet);
   fChain->SetBranchAddress("pzPFMet", pzPFMet, &b_pzPFMet);
   fChain->SetBranchAddress("vertexXPFMet", vertexXPFMet, &b_vertexXPFMet);
   fChain->SetBranchAddress("vertexYPFMet", vertexYPFMet, &b_vertexYPFMet);
   fChain->SetBranchAddress("vertexZPFMet", vertexZPFMet, &b_vertexZPFMet);
   fChain->SetBranchAddress("significancePFMet", significancePFMet, &b_significancePFMet);
   fChain->SetBranchAddress("mEtSigPFMet", mEtSigPFMet, &b_mEtSigPFMet);
   fChain->SetBranchAddress("sumEtPFMet", sumEtPFMet, &b_sumEtPFMet);
   if( isMC_ ) {
     fChain->SetBranchAddress("nGenMet", &nGenMet, &b_nGenMet);
     fChain->SetBranchAddress("chargeGenMet", chargeGenMet, &b_chargeGenMet);
     fChain->SetBranchAddress("energyGenMet", energyGenMet, &b_energyGenMet);
     fChain->SetBranchAddress("thetaGenMet", thetaGenMet, &b_thetaGenMet);
     fChain->SetBranchAddress("etaGenMet", etaGenMet, &b_etaGenMet);
     fChain->SetBranchAddress("phiGenMet", phiGenMet, &b_phiGenMet);
     fChain->SetBranchAddress("pxGenMet", pxGenMet, &b_pxGenMet);
     fChain->SetBranchAddress("pyGenMet", pyGenMet, &b_pyGenMet);
     fChain->SetBranchAddress("pzGenMet", pzGenMet, &b_pzGenMet);
     fChain->SetBranchAddress("vertexXGenMet", vertexXGenMet, &b_vertexXGenMet);
     fChain->SetBranchAddress("vertexYGenMet", vertexYGenMet, &b_vertexYGenMet);
     fChain->SetBranchAddress("vertexZGenMet", vertexZGenMet, &b_vertexZGenMet);
   }
 //fChain->SetBranchAddress("nPFCand", &nPFCand, &b_nPFCand);
 //fChain->SetBranchAddress("chargePFCand", chargePFCand, &b_chargePFCand);
 //fChain->SetBranchAddress("energyPFCand", energyPFCand, &b_energyPFCand);
 //fChain->SetBranchAddress("thetaPFCand", thetaPFCand, &b_thetaPFCand);
 //fChain->SetBranchAddress("etaPFCand", etaPFCand, &b_etaPFCand);
 //fChain->SetBranchAddress("phiPFCand", phiPFCand, &b_phiPFCand);
 //fChain->SetBranchAddress("pxPFCand", pxPFCand, &b_pxPFCand);
 //fChain->SetBranchAddress("pyPFCand", pyPFCand, &b_pyPFCand);
 //fChain->SetBranchAddress("pzPFCand", pzPFCand, &b_pzPFCand);
 //fChain->SetBranchAddress("vertexXPFCand", vertexXPFCand, &b_vertexXPFCand);
 //fChain->SetBranchAddress("vertexYPFCand", vertexYPFCand, &b_vertexYPFCand);
 //fChain->SetBranchAddress("vertexZPFCand", vertexZPFCand, &b_vertexZPFCand);
 //fChain->SetBranchAddress("particleTypePFCand", particleTypePFCand, &b_particleTypePFCand);
 //fChain->SetBranchAddress("iPFJetPFCand", iPFJetPFCand, &b_iPFJetPFCand);
   fChain->SetBranchAddress("nAK5Jet", &nAK5Jet, &b_nAK5Jet);
   fChain->SetBranchAddress("chargeAK5Jet", chargeAK5Jet, &b_chargeAK5Jet);
   fChain->SetBranchAddress("energyAK5Jet", energyAK5Jet, &b_energyAK5Jet);
   fChain->SetBranchAddress("thetaAK5Jet", thetaAK5Jet, &b_thetaAK5Jet);
   fChain->SetBranchAddress("etaAK5Jet", etaAK5Jet, &b_etaAK5Jet);
   fChain->SetBranchAddress("phiAK5Jet", phiAK5Jet, &b_phiAK5Jet);
   fChain->SetBranchAddress("pxAK5Jet", pxAK5Jet, &b_pxAK5Jet);
   fChain->SetBranchAddress("pyAK5Jet", pyAK5Jet, &b_pyAK5Jet);
   fChain->SetBranchAddress("pzAK5Jet", pzAK5Jet, &b_pzAK5Jet);
   fChain->SetBranchAddress("nAK5PFPUcorrJet", &nAK5PFPUcorrJet, &b_nAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargeAK5PFPUcorrJet", chargeAK5PFPUcorrJet, &b_chargeAK5PFPUcorrJet);
   fChain->SetBranchAddress("energyAK5PFPUcorrJet", energyAK5PFPUcorrJet, &b_energyAK5PFPUcorrJet);
   fChain->SetBranchAddress("thetaAK5PFPUcorrJet", thetaAK5PFPUcorrJet, &b_thetaAK5PFPUcorrJet);
   fChain->SetBranchAddress("etaAK5PFPUcorrJet", etaAK5PFPUcorrJet, &b_etaAK5PFPUcorrJet);
   fChain->SetBranchAddress("phiAK5PFPUcorrJet", phiAK5PFPUcorrJet, &b_phiAK5PFPUcorrJet);
   fChain->SetBranchAddress("pxAK5PFPUcorrJet", pxAK5PFPUcorrJet, &b_pxAK5PFPUcorrJet);
   fChain->SetBranchAddress("pyAK5PFPUcorrJet", pyAK5PFPUcorrJet, &b_pyAK5PFPUcorrJet);
   fChain->SetBranchAddress("pzAK5PFPUcorrJet", pzAK5PFPUcorrJet, &b_pzAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexXAK5PFPUcorrJet", vertexXAK5PFPUcorrJet, &b_vertexXAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexYAK5PFPUcorrJet", vertexYAK5PFPUcorrJet, &b_vertexYAK5PFPUcorrJet);
   fChain->SetBranchAddress("vertexZAK5PFPUcorrJet", vertexZAK5PFPUcorrJet, &b_vertexZAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFPUcorrJet", chargedHadronEnergyAK5PFPUcorrJet, &b_chargedHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFPUcorrJet", neutralHadronEnergyAK5PFPUcorrJet, &b_neutralHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("photonEnergyAK5PFPUcorrJet", photonEnergyAK5PFPUcorrJet, &b_photonEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("electronEnergyAK5PFPUcorrJet", electronEnergyAK5PFPUcorrJet, &b_electronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("muonEnergyAK5PFPUcorrJet", muonEnergyAK5PFPUcorrJet, &b_muonEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFHadronEnergyAK5PFPUcorrJet", HFHadronEnergyAK5PFPUcorrJet, &b_HFHadronEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFEMEnergyAK5PFPUcorrJet", HFEMEnergyAK5PFPUcorrJet, &b_HFEMEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedHadronMultiplicityAK5PFPUcorrJet", chargedHadronMultiplicityAK5PFPUcorrJet, &b_chargedHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralHadronMultiplicityAK5PFPUcorrJet", neutralHadronMultiplicityAK5PFPUcorrJet, &b_neutralHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("photonMultiplicityAK5PFPUcorrJet", photonMultiplicityAK5PFPUcorrJet, &b_photonMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("electronMultiplicityAK5PFPUcorrJet", electronMultiplicityAK5PFPUcorrJet, &b_electronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFPUcorrJet", muonMultiplicityAK5PFPUcorrJet, &b_muonMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFHadronMultiplicityAK5PFPUcorrJet", HFHadronMultiplicityAK5PFPUcorrJet, &b_HFHadronMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("HFEMMultiplicityAK5PFPUcorrJet", HFEMMultiplicityAK5PFPUcorrJet, &b_HFEMMultiplicityAK5PFPUcorrJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFPUcorrJet", chargedEmEnergyAK5PFPUcorrJet, &b_chargedEmEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFPUcorrJet", neutralEmEnergyAK5PFPUcorrJet, &b_neutralEmEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5PFPUcorrJet", combinedSecondaryVertexBJetTagsAK5PFPUcorrJet, &b_combinedSecondaryVertexBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet", combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet, &b_combinedSecondaryVertexMVABJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5PFPUcorrJet", jetBProbabilityBJetTagsAK5PFPUcorrJet, &b_jetBProbabilityBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5PFPUcorrJet", jetProbabilityBJetTagsAK5PFPUcorrJet, &b_jetProbabilityBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet", simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet, &b_simpleSecondaryVertexHighEffBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet", simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet, &b_simpleSecondaryVertexHighPurBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5PFPUcorrJet", softMuonBJetTagsAK5PFPUcorrJet, &b_softMuonBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonByIP3dBJetTagsAK5PFPUcorrJet", softMuonByIP3dBJetTagsAK5PFPUcorrJet, &b_softMuonByIP3dBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softMuonByPtBJetTagsAK5PFPUcorrJet", softMuonByPtBJetTagsAK5PFPUcorrJet, &b_softMuonByPtBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronBJetTagsAK5PFPUcorrJet", softElectronBJetTagsAK5PFPUcorrJet, &b_softElectronBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronByIP3dBJetTagsAK5PFPUcorrJet", softElectronByIP3dBJetTagsAK5PFPUcorrJet, &b_softElectronByIP3dBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("softElectronByPtBJetTagsAK5PFPUcorrJet", softElectronByPtBJetTagsAK5PFPUcorrJet, &b_softElectronByPtBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5PFPUcorrJet", trackCountingHighPurBJetTagsAK5PFPUcorrJet, &b_trackCountingHighPurBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5PFPUcorrJet", trackCountingHighEffBJetTagsAK5PFPUcorrJet, &b_trackCountingHighEffBJetTagsAK5PFPUcorrJet);
   fChain->SetBranchAddress("uncorrEnergyAK5PFPUcorrJet", uncorrEnergyAK5PFPUcorrJet, &b_uncorrEnergyAK5PFPUcorrJet);
   fChain->SetBranchAddress("ptDAK5PFPUcorrJet", ptDAK5PFPUcorrJet, &b_ptDAK5PFPUcorrJet);
   fChain->SetBranchAddress("rmsCandAK5PFPUcorrJet", rmsCandAK5PFPUcorrJet, &b_rmsCandAK5PFPUcorrJet);
   fChain->SetBranchAddress("axis1AK5PFPUcorrJet", axis1AK5PFPUcorrJet, &b_axis1AK5PFPUcorrJet);
   fChain->SetBranchAddress("axis2AK5PFPUcorrJet", axis2AK5PFPUcorrJet, &b_axis2AK5PFPUcorrJet);
   fChain->SetBranchAddress("pullAK5PFPUcorrJet", pullAK5PFPUcorrJet, &b_pullAK5PFPUcorrJet);
   fChain->SetBranchAddress("tanaAK5PFPUcorrJet", tanaAK5PFPUcorrJet, &b_tanaAK5PFPUcorrJet);
   fChain->SetBranchAddress("ptD_QCAK5PFPUcorrJet", ptD_QCAK5PFPUcorrJet, &b_ptD_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("rmsCand_QCAK5PFPUcorrJet", rmsCand_QCAK5PFPUcorrJet, &b_rmsCand_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("axis1_QCAK5PFPUcorrJet", axis1_QCAK5PFPUcorrJet, &b_axis1_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("axis2_QCAK5PFPUcorrJet", axis2_QCAK5PFPUcorrJet, &b_axis2_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("pull_QCAK5PFPUcorrJet", pull_QCAK5PFPUcorrJet, &b_pull_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("tana_QCAK5PFPUcorrJet", tana_QCAK5PFPUcorrJet, &b_tana_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("nChg_ptCutAK5PFPUcorrJet", nChg_ptCutAK5PFPUcorrJet, &b_nChg_ptCutAK5PFPUcorrJet);
   fChain->SetBranchAddress("nChg_QCAK5PFPUcorrJet", nChg_QCAK5PFPUcorrJet, &b_nChg_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("nChg_ptCut_QCAK5PFPUcorrJet", nChg_ptCut_QCAK5PFPUcorrJet, &b_nChg_ptCut_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("nNeutral_ptCutAK5PFPUcorrJet", nNeutral_ptCutAK5PFPUcorrJet, &b_nNeutral_ptCutAK5PFPUcorrJet);
   fChain->SetBranchAddress("RchgAK5PFPUcorrJet", RchgAK5PFPUcorrJet, &b_RchgAK5PFPUcorrJet);
   fChain->SetBranchAddress("Rchg_QCAK5PFPUcorrJet", Rchg_QCAK5PFPUcorrJet, &b_Rchg_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("RneutralAK5PFPUcorrJet", RneutralAK5PFPUcorrJet, &b_RneutralAK5PFPUcorrJet);
   fChain->SetBranchAddress("RAK5PFPUcorrJet", RAK5PFPUcorrJet, &b_RAK5PFPUcorrJet);
   fChain->SetBranchAddress("RchgAK5PFPUcorrJet", RchgAK5PFPUcorrJet, &b_RchgAK5PFPUcorrJet);
   fChain->SetBranchAddress("Rchg_QCAK5PFPUcorrJet", Rchg_QCAK5PFPUcorrJet, &b_Rchg_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("pTMaxAK5PFPUcorrJet", pTMaxAK5PFPUcorrJet, &b_pTMaxAK5PFPUcorrJet);
   fChain->SetBranchAddress("pTMaxChgAK5PFPUcorrJet", pTMaxChgAK5PFPUcorrJet, &b_pTMaxChgAK5PFPUcorrJet);
   fChain->SetBranchAddress("pTMaxNeutralAK5PFPUcorrJet", pTMaxNeutralAK5PFPUcorrJet, &b_pTMaxNeutralAK5PFPUcorrJet);
   fChain->SetBranchAddress("pTMaxChg_QCAK5PFPUcorrJet", pTMaxChg_QCAK5PFPUcorrJet, &b_pTMaxChg_QCAK5PFPUcorrJet);
   fChain->SetBranchAddress("betaAK5PFPUcorrJet", betaAK5PFPUcorrJet, &b_betaAK5PFPUcorrJet);
   fChain->SetBranchAddress("betastarAK5PFPUcorrJet", betastarAK5PFPUcorrJet, &b_betastarAK5PFPUcorrJet);
   fChain->SetBranchAddress("rmsCandsHandAK5PFPUcorrJet", rmsCandsHandAK5PFPUcorrJet, &b_rmsCandsHandAK5PFPUcorrJet);
   fChain->SetBranchAddress("nAK5PFNoPUJet", &nAK5PFNoPUJet, &b_nAK5PFNoPUJet);
   fChain->SetBranchAddress("chargeAK5PFNoPUJet", chargeAK5PFNoPUJet, &b_chargeAK5PFNoPUJet);
   fChain->SetBranchAddress("energyAK5PFNoPUJet", energyAK5PFNoPUJet, &b_energyAK5PFNoPUJet);
   fChain->SetBranchAddress("thetaAK5PFNoPUJet", thetaAK5PFNoPUJet, &b_thetaAK5PFNoPUJet);
   fChain->SetBranchAddress("etaAK5PFNoPUJet", etaAK5PFNoPUJet, &b_etaAK5PFNoPUJet);
   fChain->SetBranchAddress("phiAK5PFNoPUJet", phiAK5PFNoPUJet, &b_phiAK5PFNoPUJet);
   fChain->SetBranchAddress("pxAK5PFNoPUJet", pxAK5PFNoPUJet, &b_pxAK5PFNoPUJet);
   fChain->SetBranchAddress("pyAK5PFNoPUJet", pyAK5PFNoPUJet, &b_pyAK5PFNoPUJet);
   fChain->SetBranchAddress("pzAK5PFNoPUJet", pzAK5PFNoPUJet, &b_pzAK5PFNoPUJet);
   fChain->SetBranchAddress("vertexXAK5PFNoPUJet", vertexXAK5PFNoPUJet, &b_vertexXAK5PFNoPUJet);
   fChain->SetBranchAddress("vertexYAK5PFNoPUJet", vertexYAK5PFNoPUJet, &b_vertexYAK5PFNoPUJet);
   fChain->SetBranchAddress("vertexZAK5PFNoPUJet", vertexZAK5PFNoPUJet, &b_vertexZAK5PFNoPUJet);
   fChain->SetBranchAddress("chargedHadronEnergyAK5PFNoPUJet", chargedHadronEnergyAK5PFNoPUJet, &b_chargedHadronEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("neutralHadronEnergyAK5PFNoPUJet", neutralHadronEnergyAK5PFNoPUJet, &b_neutralHadronEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("photonEnergyAK5PFNoPUJet", photonEnergyAK5PFNoPUJet, &b_photonEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("electronEnergyAK5PFNoPUJet", electronEnergyAK5PFNoPUJet, &b_electronEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("muonEnergyAK5PFNoPUJet", muonEnergyAK5PFNoPUJet, &b_muonEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("HFHadronEnergyAK5PFNoPUJet", HFHadronEnergyAK5PFNoPUJet, &b_HFHadronEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("HFEMEnergyAK5PFNoPUJet", HFEMEnergyAK5PFNoPUJet, &b_HFEMEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("chargedHadronMultiplicityAK5PFNoPUJet", chargedHadronMultiplicityAK5PFNoPUJet, &b_chargedHadronMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("neutralHadronMultiplicityAK5PFNoPUJet", neutralHadronMultiplicityAK5PFNoPUJet, &b_neutralHadronMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("photonMultiplicityAK5PFNoPUJet", photonMultiplicityAK5PFNoPUJet, &b_photonMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("electronMultiplicityAK5PFNoPUJet", electronMultiplicityAK5PFNoPUJet, &b_electronMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("muonMultiplicityAK5PFNoPUJet", muonMultiplicityAK5PFNoPUJet, &b_muonMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("HFHadronMultiplicityAK5PFNoPUJet", HFHadronMultiplicityAK5PFNoPUJet, &b_HFHadronMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("HFEMMultiplicityAK5PFNoPUJet", HFEMMultiplicityAK5PFNoPUJet, &b_HFEMMultiplicityAK5PFNoPUJet);
   fChain->SetBranchAddress("chargedEmEnergyAK5PFNoPUJet", chargedEmEnergyAK5PFNoPUJet, &b_chargedEmEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("neutralEmEnergyAK5PFNoPUJet", neutralEmEnergyAK5PFNoPUJet, &b_neutralEmEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("combinedSecondaryVertexBJetTagsAK5PFNoPUJet", combinedSecondaryVertexBJetTagsAK5PFNoPUJet, &b_combinedSecondaryVertexBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("combinedSecondaryVertexMVABJetTagsAK5PFNoPUJet", combinedSecondaryVertexMVABJetTagsAK5PFNoPUJet, &b_combinedSecondaryVertexMVABJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("jetBProbabilityBJetTagsAK5PFNoPUJet", jetBProbabilityBJetTagsAK5PFNoPUJet, &b_jetBProbabilityBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("jetProbabilityBJetTagsAK5PFNoPUJet", jetProbabilityBJetTagsAK5PFNoPUJet, &b_jetProbabilityBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet", simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet, &b_simpleSecondaryVertexHighEffBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("simpleSecondaryVertexHighPurBJetTagsAK5PFNoPUJet", simpleSecondaryVertexHighPurBJetTagsAK5PFNoPUJet, &b_simpleSecondaryVertexHighPurBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("softMuonBJetTagsAK5PFNoPUJet", softMuonBJetTagsAK5PFNoPUJet, &b_softMuonBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("softMuonByIP3dBJetTagsAK5PFNoPUJet", softMuonByIP3dBJetTagsAK5PFNoPUJet, &b_softMuonByIP3dBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("softMuonByPtBJetTagsAK5PFNoPUJet", softMuonByPtBJetTagsAK5PFNoPUJet, &b_softMuonByPtBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("softElectronBJetTagsAK5PFNoPUJet", softElectronBJetTagsAK5PFNoPUJet, &b_softElectronBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("softElectronByIP3dBJetTagsAK5PFNoPUJet", softElectronByIP3dBJetTagsAK5PFNoPUJet, &b_softElectronByIP3dBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("softElectronByPtBJetTagsAK5PFNoPUJet", softElectronByPtBJetTagsAK5PFNoPUJet, &b_softElectronByPtBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("trackCountingHighPurBJetTagsAK5PFNoPUJet", trackCountingHighPurBJetTagsAK5PFNoPUJet, &b_trackCountingHighPurBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("trackCountingHighEffBJetTagsAK5PFNoPUJet", trackCountingHighEffBJetTagsAK5PFNoPUJet, &b_trackCountingHighEffBJetTagsAK5PFNoPUJet);
   fChain->SetBranchAddress("uncorrEnergyAK5PFNoPUJet", uncorrEnergyAK5PFNoPUJet, &b_uncorrEnergyAK5PFNoPUJet);
   fChain->SetBranchAddress("ptDAK5PFNoPUJet", ptDAK5PFNoPUJet, &b_ptDAK5PFNoPUJet);
   fChain->SetBranchAddress("rmsCandAK5PFNoPUJet", rmsCandAK5PFNoPUJet, &b_rmsCandAK5PFNoPUJet);
   fChain->SetBranchAddress("axis1AK5PFNoPUJet", axis1AK5PFNoPUJet, &b_axis1AK5PFNoPUJet);
   fChain->SetBranchAddress("axis2AK5PFNoPUJet", axis2AK5PFNoPUJet, &b_axis2AK5PFNoPUJet);
   fChain->SetBranchAddress("pullAK5PFNoPUJet", pullAK5PFNoPUJet, &b_pullAK5PFNoPUJet);
   fChain->SetBranchAddress("tanaAK5PFNoPUJet", tanaAK5PFNoPUJet, &b_tanaAK5PFNoPUJet);
   fChain->SetBranchAddress("ptD_QCAK5PFNoPUJet", ptD_QCAK5PFNoPUJet, &b_ptD_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("rmsCand_QCAK5PFNoPUJet", rmsCand_QCAK5PFNoPUJet, &b_rmsCand_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("axis1_QCAK5PFNoPUJet", axis1_QCAK5PFNoPUJet, &b_axis1_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("axis2_QCAK5PFNoPUJet", axis2_QCAK5PFNoPUJet, &b_axis2_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("pull_QCAK5PFNoPUJet", pull_QCAK5PFNoPUJet, &b_pull_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("tana_QCAK5PFNoPUJet", tana_QCAK5PFNoPUJet, &b_tana_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("nChg_ptCutAK5PFNoPUJet", nChg_ptCutAK5PFNoPUJet, &b_nChg_ptCutAK5PFNoPUJet);
   fChain->SetBranchAddress("nChg_QCAK5PFNoPUJet", nChg_QCAK5PFNoPUJet, &b_nChg_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("nChg_ptCut_QCAK5PFNoPUJet", nChg_ptCut_QCAK5PFNoPUJet, &b_nChg_ptCut_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("nNeutral_ptCutAK5PFNoPUJet", nNeutral_ptCutAK5PFNoPUJet, &b_nNeutral_ptCutAK5PFNoPUJet);
   fChain->SetBranchAddress("RchgAK5PFNoPUJet", RchgAK5PFNoPUJet, &b_RchgAK5PFNoPUJet);
   fChain->SetBranchAddress("RneutralAK5PFNoPUJet", RneutralAK5PFNoPUJet, &b_RneutralAK5PFNoPUJet);
   fChain->SetBranchAddress("RAK5PFNoPUJet", RAK5PFNoPUJet, &b_RAK5PFNoPUJet);
   fChain->SetBranchAddress("RchgAK5PFNoPUJet", RchgAK5PFNoPUJet, &b_RchgAK5PFNoPUJet);
   fChain->SetBranchAddress("pTMaxAK5PFNoPUJet", pTMaxAK5PFNoPUJet, &b_pTMaxAK5PFNoPUJet);
   fChain->SetBranchAddress("pTMaxChgAK5PFNoPUJet", pTMaxChgAK5PFNoPUJet, &b_pTMaxChgAK5PFNoPUJet);
   fChain->SetBranchAddress("pTMaxNeutralAK5PFNoPUJet", pTMaxNeutralAK5PFNoPUJet, &b_pTMaxNeutralAK5PFNoPUJet);
   fChain->SetBranchAddress("pTMaxChg_QCAK5PFNoPUJet", pTMaxChg_QCAK5PFNoPUJet, &b_pTMaxChg_QCAK5PFNoPUJet);
   fChain->SetBranchAddress("betaAK5PFNoPUJet", betaAK5PFNoPUJet, &b_betaAK5PFNoPUJet);
   fChain->SetBranchAddress("betastarAK5PFNoPUJet", betastarAK5PFNoPUJet, &b_betastarAK5PFNoPUJet);
   fChain->SetBranchAddress("rmsCandsHandAK5PFNoPUJet", rmsCandsHandAK5PFNoPUJet, &b_rmsCandsHandAK5PFNoPUJet);
   if( isMC_ ) {
     fChain->SetBranchAddress("nAK5GenJet", &nAK5GenJet, &b_nAK5GenJet);
     fChain->SetBranchAddress("chargeAK5GenJet", chargeAK5GenJet, &b_chargeAK5GenJet);
     fChain->SetBranchAddress("energyAK5GenJet", energyAK5GenJet, &b_energyAK5GenJet);
     fChain->SetBranchAddress("thetaAK5GenJet", thetaAK5GenJet, &b_thetaAK5GenJet);
     fChain->SetBranchAddress("etaAK5GenJet", etaAK5GenJet, &b_etaAK5GenJet);
     fChain->SetBranchAddress("phiAK5GenJet", phiAK5GenJet, &b_phiAK5GenJet);
     fChain->SetBranchAddress("pxAK5GenJet", pxAK5GenJet, &b_pxAK5GenJet);
     fChain->SetBranchAddress("pyAK5GenJet", pyAK5GenJet, &b_pyAK5GenJet);
     fChain->SetBranchAddress("pzAK5GenJet", pzAK5GenJet, &b_pzAK5GenJet);
     fChain->SetBranchAddress("vertexXAK5GenJet", vertexXAK5GenJet, &b_vertexXAK5GenJet);
     fChain->SetBranchAddress("vertexYAK5GenJet", vertexYAK5GenJet, &b_vertexYAK5GenJet);
     fChain->SetBranchAddress("vertexZAK5GenJet", vertexZAK5GenJet, &b_vertexZAK5GenJet);
     fChain->SetBranchAddress("genPtHat", &genPtHat, &b_genPtHat);
     fChain->SetBranchAddress("genProcessId", &genProcessId, &b_genProcessId);
     fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
     fChain->SetBranchAddress("genAlphaQCD", &genAlphaQCD, &b_genAlphaQCD);
     fChain->SetBranchAddress("genAlphaQED", &genAlphaQED, &b_genAlphaQED);
   }

   Notify();
}

Bool_t Ntp1Analyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ntp1Analyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ntp1Analyzer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


void Ntp1Analyzer::UpdateCache() {
   // cache current run for lumi measurement:
   if( oldrun_ != runNumber ) {
     oldrun_ = runNumber;
     goodLSCache_ = goodLS_.find( runNumber );
   }
}


void Ntp1Analyzer::ReadCSVFile(const std::string& csv) {

  std::cout << "Reading CSV file of LS luminosities " << csv << std::endl;
  FILE* iff = fopen(csv.c_str(),"r");

  if(iff == 0) {
    std::cout << "cannot open CSV file " << csv << " ... now exiting." << std::endl;
    exit(-2);
  }

  int irun, iLS, iHFcounts, iVtxcounts;
  double iHFlumi, iVtxlumi;

  while( fscanf(iff,"%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%d,%*[ \t\n]%lf,%*[ \t\n]%lf", &irun, &iLS, &iHFcounts, &iVtxcounts, &iHFlumi, &iVtxlumi) != EOF ) {

    RunLumiPair rlpair( std::pair<int, int>(irun, iLS) );
    double istLumi = (iVtxlumi>0.) ? iVtxlumi : iHFlumi; //using Vtx lumi as default
    double intLumi = 23.*istLumi/1.0e30; //in microb-1 (1 LS = 23 seconds)
    LSLumimap_[rlpair] = intLumi; 

  }
  fclose(iff);

}


void Ntp1Analyzer::ReadJSONFile(const std::string& json) {

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


bool Ntp1Analyzer::isGoodEvent( int iEntry ) {

     bool okForJSON = false;
     bool okForHLT = false;

     // FIRST STEP: check if LS is a good one in the JSON file:

     if(!filterGoodRuns_) okForJSON = true; // if filtered not requested all events are good

     this->UpdateCache();


     // check whether this run is part of the good runs. else retrun false
     if( goodLSCache_ != goodLS_.end() ) {

        // get list of LS intervals
        const GoodLSVector& lsvector =   goodLSCache_->second; 
        // loop over good LS intervals and return as soon as one interval contains this event
        for(GoodLSVector::const_iterator iLS = lsvector.begin(); (iLS != lsvector.end())&&(okForJSON==false); iLS++) {
     
           if(lumiBlock >= iLS->first && lumiBlock <= iLS->second ) {
    
             okForJSON = true;
          } // check current LS being in the interval
        } // loop over good LS for this run
     }

     // SECOND STEP: if it's ok in the JSON, check if it has passed the trigger
     // once per run, reload trigger mask:
     if( cachedRun_!=runNumber ) {
       if( cachedRun_==0 )
         std::cout << "-> Loading Trigger Mask for run " << runNumber << "." << std::endl;
       else 
         std::cout << "-> Passing from run " << cachedRun_ << " to run " << runNumber << ". Reloading Trigger Mask." << std::endl;
       cachedRun_ = runNumber;
       this->LoadTrigger( iEntry );
     }

     okForHLT = this->PassedHLT(iEntry);

  // if( okForJSON && okForHLT ) { //will take lumi, so (once per LS) increment luminosity
  
  //   if( cachedLS_ != lumiBlock )  {
  //      cachedLS_ = lumiBlock;
  //      RunLumiPair rlpair = (std::pair<int, int>(runNumber,lumiBlock));
  //      totalIntLumi_ += LSLumimap_[rlpair];
  //   }
  // }

     bool returnBool = ( okForJSON && okForHLT );
     return returnBool;

}







double Ntp1Analyzer::trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz) {
  float elePt = sqrt(elePx*elePx + elePy*elePy);
  return ( - (eleVx-PVx)*elePy + (eleVy-PVy)*elePx ) / elePt;
}



double Ntp1Analyzer::eleDzPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDzPV(PVPos,lepVPos,lepMom);
}


double Ntp1Analyzer::muonDzPV(int imu, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int ctfMuon   = trackIndexMuon[imu];
  TVector3 lepVPos(trackVxTrack[ctfMuon],trackVyTrack[ctfMuon],trackVzTrack[ctfMuon]);
  TVector3 lepMom(pxMuon[imu],pyMuon[imu],pzMuon[imu]);
  return trackDzPV(PVPos,lepVPos,lepMom);
}


/// dz parameter with respect to PV for tracks
double Ntp1Analyzer::trackDzPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  float trackPt = trackMom.Pt();
  return (trackVPos.Z()-PVPos.Z()) - ((trackVPos.X()-PVPos.X())*trackMom.X()+(trackVPos.Y()-PVPos.Y())*trackMom.Y())/trackPt *trackMom.Pz()/trackPt;
}
