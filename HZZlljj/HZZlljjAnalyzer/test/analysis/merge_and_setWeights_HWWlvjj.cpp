#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TRegexp.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <cstdlib>


TChain* tree = 0;
std::string analysisType_;
std::string flags_;




float addInput( const std::string& dataset );
float getWeight( const std::string& dataset, int nEvents );


int main( int argc, char* argv[] ) {

  if( argc!=2 && argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./merge_and_setWeights [dataset] [analysisType=\"HWWlvjj\"] [flags=\"\"]" << std::endl;
    exit(917);
  }

  std::string dataset = argv[1];

  analysisType_ = "HWWlvjj";
  if( argc>=3 ) {
    std::string analysisType_str(argv[2]);
    analysisType_ = analysisType_str;
  }

  flags_ = "";
  if( argc==4 ) {
    std::string flags_str(argv[3]);
    flags_ = flags_str;
  }


  tree = new TChain("reducedTree");

  float nTotalEvents = addInput( dataset );

  std::cout << std::endl << "-> Finished adding. Total entries: " << tree->GetEntries() << std::endl;

  float weight = getWeight( dataset, nTotalEvents );

  float nTotalEventsW = (float)nTotalEvents*weight;

  // and now set the weights
  tree->SetBranchStatus( "eventWeight", 0 );
  
  std::string outfilename = analysisType_ + "_2ndLevelTreeW_"+dataset;
  if( flags_!="" ) outfilename += "_" + flags_;
  outfilename += +".root";
  TFile* outfile = new TFile(outfilename.c_str(), "recreate");
  outfile->cd();

  TH1F* h1_nCounter = new TH1F("nCounter", "", 1, 0., 1.);
  h1_nCounter->SetBinContent(1, nTotalEvents);
  TH1F* h1_nCounterW = new TH1F("nCounterW", "", 1, 0., 1.);
  h1_nCounterW->SetBinContent(1, nTotalEventsW);

  TTree* newTree = tree->CloneTree(0);
  Float_t newWeight;
  newTree->Branch( "eventWeight", &newWeight, "newWeight/F" );

  int nentries = tree->GetEntries();
  for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( (ientry % 100000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

    newWeight = weight;

    if( dataset=="MU_Run2010B_PromptReco_v2_runs146240_146733" ) newWeight = 0.5;

    newTree->Fill();

  } //for entries

  h1_nCounter->Write();
  h1_nCounterW->Write();
  newTree->Write();
  outfile->Write();
  outfile->Close();

  return 0;

}


float addInput( const std::string& dataset ) {

  std::string infileName = "files_"+analysisType_+"_2ndLevel_" + dataset;
  if( flags_!="" ) infileName += "_" + flags_;
  infileName += ".txt";
  TH1F* h1_nCounter;
  TH1F* h1_nCounter_Zee;
  TH1F* h1_nCounter_Zmumu;

  float totalEvents = 0;
  // these are needed to fix the Spring11 ZJets Zee BR bug:
  int totalEvents_Zee = 0;
  int totalEvents_Zmumu = 0;

  //open from file.txt:
  FILE* iff = fopen(infileName.c_str(),"r");
  if(iff == 0) {
    std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    infileName = analysisType_+"_2ndLevelTree_" + dataset;
    if( flags_!="" ) infileName += "_" + flags_;
    infileName += ".root";
    std::string treeName = infileName +"/reducedTree";
    tree->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
    TFile* infile = TFile::Open(infileName.c_str(), "READ");
    h1_nCounter = (TH1F*)infile->Get("nCounter");
    if( h1_nCounter!=0 ) {
      totalEvents += h1_nCounter->GetBinContent(1);
    } else {
      std::cout << " WARNING! File '" << infileName << "' has no nCounter information. Skipping." << std::endl;
    }
    h1_nCounter_Zee = (TH1F*)infile->Get("nCounter_Zee");
    if( h1_nCounter_Zee!=0 ) {
      totalEvents_Zee += h1_nCounter_Zee->GetBinContent(1);
    }
    h1_nCounter_Zmumu = (TH1F*)infile->Get("nCounter_Zmumu");
    if( h1_nCounter_Zmumu!=0 ) {
      totalEvents_Zmumu += h1_nCounter_Zmumu->GetBinContent(1);
    }
    std::cout << std::endl;
    infile->Close();

  } else { //if file is good:

    char singleLine[500];
    std::cout << "-> Correctly opened file: '" << infileName << "'." << std::endl;

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/reducedTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_nCounter = (TH1F*)infile->Get("nCounter");
      if( h1_nCounter!=0 ) {
        totalEvents += h1_nCounter->GetBinContent(1);
      } else {
        std::cout << std::endl << " WARNING! File '" << rootfilename << "' has no nCounter information. Skipping." << std::endl;
      }
      h1_nCounter_Zee = (TH1F*)infile->Get("nCounter_Zee");
      if( h1_nCounter_Zee!=0 ) {
        totalEvents_Zee += h1_nCounter_Zee->GetBinContent(1);
      }
      h1_nCounter_Zmumu = (TH1F*)infile->Get("nCounter_Zmumu");
      if( h1_nCounter_Zmumu!=0 ) {
        totalEvents_Zmumu += h1_nCounter_Zmumu->GetBinContent(1);
      }
      std::cout << std::endl;
      infile->Close();

    }
    fclose(iff);

  }

  //correct bug in Zee BR (Spring11 alpgen ZJets) only:
  TString dataset_tstr(dataset);
  if( dataset_tstr.Contains("Spring11") && ( dataset_tstr.BeginsWith("Z0Jet") 
                                          || dataset_tstr.BeginsWith("Z1Jet") 
                                          || dataset_tstr.BeginsWith("Z2Jet")
                                          || dataset_tstr.BeginsWith("Z3Jet")
                                          || dataset_tstr.BeginsWith("Z4Jet")
                                          || dataset_tstr.BeginsWith("Z5Jet") ) ) {

    std::cout << std::endl << "!!! Correcting bug in Z->ee BR (Spring11 Alpgen Z+Jets samples) " << std::endl;
    // reduces to = totalEvents when totalEvents_Zmumu+totalEvents_Zee = 2/3 totalEvents:
    //totalEvents = 2.*totalEvents*totalEvents / ( 3.* (float)(totalEvents_Zmumu+totalEvents_Zee) );
    totalEvents = 3.*(totalEvents_Zmumu+totalEvents_Zee)/2.;

  }

  return totalEvents;

} //addinput


float getWeight( const std::string& dataset, int nEvents ) {

  TString dataset_tstr(dataset);
  float xSection = -1.;

  bool isAlpgenZJets = false;

  // all cross sections in pb-1:
  if( dataset=="ZJets_madgraph" || dataset=="DYJetsToLL_TuneZ2_M-50_madgraph_Fall10" ) {
    xSection = 3048.; //NNLO see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset=="Z0Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 2350.*0.853 ; // sigma x filter efficiency taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionReProcessingSpring10#ALPGEN
  } else if( dataset=="Z1Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 870.*0.447 ; // sigma x filter efficiency
  } else if( dataset=="Z1Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 19.3*0.504; // sigma x filter efficiency
  } else if( dataset=="Z1Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.226*0.352; // sigma x filter efficiency
  } else if( dataset=="Z1Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.000528*0.266; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 372.*0.26; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 27.0*0.315; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.457*0.244; // sigma x filter efficiency
  } else if( dataset=="Z2Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.00132*0.203; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 140.*0.157; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 20.3*0.189; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.465*0.162; // sigma x filter efficiency
  } else if( dataset=="Z3Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.00152*0.149; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 46.1*0.0939; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 10.7*0.115; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.319*0.104; // sigma x filter efficiency
  } else if( dataset=="Z4Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.0011*0.106; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt0to100-alpgen_Spring10" ) {
    xSection = 13.9*0.0727; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt100to300-alpgen_Spring10" ) {
    xSection = 4.42*0.0956; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt300to800-alpgen_Spring10" ) {
    xSection = 0.164*0.103; // sigma x filter efficiency
  } else if( dataset=="Z5Jets_Pt800to1600-alpgen_Spring10" ) {
    xSection = 0.000588*0.109; // sigma x filter efficiency
  } else if( dataset_tstr.BeginsWith("Z0Jets") ) {
    xSection = 1929. ; // cross sections taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionFall2010#ALPGEN
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z1Jets_ptZ-0to100") ) {
    xSection = 380.8; 
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z1Jets_ptZ-100to300") ) {
    xSection = 8.721;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z1Jets_ptZ-300to800") ) {
    xSection = 0.07386;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z1Jets_ptZ-800to1600") ) {
    xSection = 0.0001374;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z2Jets_ptZ-0to100") ) {
    xSection = 103.5;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z2Jets_ptZ-100to300") ) {
    xSection = 8.534;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z2Jets_ptZ-300to800") ) {
    xSection = 0.1151;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z2Jets_ptZ-800to1600") ) {
    xSection = 0.0003023;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z3Jets_ptZ-0to100") ) {
    xSection = 22.89; 
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z3Jets_ptZ-100to300") ) {
    xSection = 3.951;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z3Jets_ptZ-300to800") ) {
    xSection = 0.08344;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z3Jets_ptZ-800to1600") ) {
    xSection = 0.0002480;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z4Jets_ptZ-0to100") ) {
    xSection = 4.619;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z4Jets_ptZ-100to300") ) {
    xSection = 1.298;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z4Jets_ptZ-300to800") ) {
    xSection = 3.935e-02;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z4Jets_ptZ-800to1600") ) {
    xSection = 1.394e-04;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z5Jets_ptZ-0to100") ) {
    xSection = 1.135;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z5Jets_ptZ-100to300") ) {
    xSection = 4.758e-01;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z5Jets_ptZ-300to800") ) {
    xSection = 1.946e-02;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("Z5Jets_ptZ-800to1600") ) {
    xSection = 1.374e-04;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZBB0Jets") ) {
    xSection = 1.703;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZBB1Jets") ) {
    xSection = 9.620e-01;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZBB2Jets") ) {
    xSection = 3.639e-01;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZBB3Jets") ) {
    xSection = 1.598e-01;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZCC0Jets") ) {
    xSection = 1.707;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZCC1Jets") ) {
    xSection = 9.526e-01;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZCC2Jets") ) {
    xSection = 3.659e-01;
    isAlpgenZJets = true;
  } else if( dataset_tstr.BeginsWith("ZCC3Jets") ) {
    xSection = 1.643e-01;
    isAlpgenZJets = true;
  } else if( dataset=="HZZ_qqll_gluonfusion_M130" ) {
    xSection = 25.560*0.03913*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M150" ) {
    xSection = 19.568*0.08234*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M200" ) {
    //xSection = 10.361*0.2537*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
    xSection = 10.361*0.2537*0.10097*0.7*2.*40.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2 ENHANCE SIGMA BY 40!!!!!
  } else if( dataset=="HZZ_qqll_gluonfusion_M300" ) {
    xSection = 5.2728*0.3053*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M400" ) {
    xSection = 4.8236*0.2664*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="HZZ_qqll_gluonfusion_M500" ) {
    xSection = 2.1914*0.2602*0.10097*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM200_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-200_7TeV-jhu-pythia6") ) {
    xSection = (5.249+0.6371+0.1032+0.06096+0.01849)*0.2537*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM250_2l2j_FASTSIM"  || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-250_7TeV-jhu-pythia6")) {
    xSection = (3.312+0.4304+0.04308+0.02540+0.008593)*0.2951*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM300_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-300_7TeV-jhu-pythia6") ) {
    xSection = (2.418+0.3010+0.02018+0.01169+0.004719)*0.3053*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM350_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-350_7TeV-jhu-pythia6") ) {
    xSection = (2.299599+0.21635+0.0021911+0.0055876+0.010794)*0.3023*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM400_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-400_7TeV-jhu-pythia6") ) {
    xSection = (2.035+0.1619+0.0054837+0.00030308+0.0059045)*0.2664*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM450_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-450_7TeV-jhu-pythia6") ) {
    xSection = (1.356+0.1235+0.00087374+0.0017291+0.00033991)*0.2582*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM500_2l2j"|| dataset=="JHUgen_HiggsSM500_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-500_7TeV-jhu-pythia6") ) {
    xSection = (0.8497+0.09491+0.00058401+0.0010272+0.0020377)*0.2602*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM550_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-550_7TeV-jhu-pythia6") ) {
    xSection = (0.5259+0.07356+0.00039969+0.00063066+0.0012623)*0.2657*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="JHUgen_HiggsSM600_2l2j_FASTSIM" || dataset_tstr.Contains("SMHiggsToZZTo2L2Q_M-600_7TeV-jhu-pythia6") ) {
    xSection = (0.3275+0.05763+0.00027833+0.00039795+0.00080323)*0.2724*0.067316*0.7*2.; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->jj) x 2
  } else if( dataset=="GluGluToHToZZTo4L_M-400_7TeV-powheg-pythia6_Fall10" ) {
    xSection = (2.0608)*0.2724*0.100974*0.100974; //sigma x BR(H->ZZ) x BR(Z->ll) x BR(Z->ll) (l=e,m,t)
  } else if( dataset_tstr.BeginsWith("TTJets") || dataset_tstr.BeginsWith("TT_") ) {
    xSection = 157.4; //NLO see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset_tstr.BeginsWith("ZZtoAnything") ) {
    xSection = 5.9*1.3; //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt plus factor 1.3 to account for glu-glu
  } else if( dataset=="Zmumu_Pythia" ) {
    xSection = 3048./3.; //NNLO see https://twiki.cern.ch/twiki/pub/CMS/GeneratorMain/ShortXsec.pdf
  } else if( dataset=="PhotonJet_Summer1036X_Pt5to15_pfakt5" ) {
    xSection = 4030000.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt15to20_pfakt5" ) {
    xSection = 114700.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt20to30_pfakt5" ) {
    xSection = 57180.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt30to50_pfakt5" ) { 
    xSection = 16520.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt50to80_pfakt5" ) {
    xSection = 2723.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt80to120_pfakt5" ) {
    xSection = 446.2;
  } else if( dataset=="PhotonJet_Summer1036X_Pt120to170_pfakt5" ) {
    xSection = 84.43;
  } else if( dataset=="PhotonJet_Summer1036X_Pt170to300_pfakt5" ) {
    xSection = 22.55;
  } else if( dataset=="PhotonJet_Summer1036X_Pt300to500_pfakt5" ) {
    xSection = 1.545;
  } else if( dataset=="PhotonJet_Summer1036X_Pt500toInf_pfakt5" ) {
    xSection = 0.0923;
  } else if( dataset=="QCD_Pt_120to170_TuneZ2_7TeV_pythia6" ) {
    xSection = 1.151e+05;
  } else if( dataset=="QCD_Pt_170to300_TuneZ2_7TeV_pythia6" ) {
    xSection = 2.426e+04;
  } else if( dataset=="WW200" ) {//##
    xSection = 1.1202;//##
  } else if( dataset=="WW300" ) {//##
    xSection = 0.4823;//##
  } else if( dataset=="WW400" ) {//##
    xSection = 0.341;//##
  } else if( dataset=="WW500" ) {//##
    xSection = 0.13372;//##
  } else if( dataset=="GluGluToHToWWToLNuQQ_M-250_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1" ) {//##
    xSection = 0.668653;//##
  } else if( dataset=="GluGluToHToWWToLNuQQ_M-350_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1" ) {//##
    xSection = 0.44895;//##
  } else if( dataset=="GluGluToHToWWToLNuQQ_M-450_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1" ) {//##
    xSection = 0.288;//##
  } else if( dataset=="GluGluToHToWWToLNuQQ_M-550_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1" ) {//##
    xSection = 0.07673;//##
  } else if( dataset=="GluGlu170" ) {//##
    xSection = 2.148;//## Other
  } else if( (dataset=="TT_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_3") || (dataset=="TTJets_TuneZ2_7TeV-madgraph-tauola") ) {//## t tbar
    xSection = 157.5;//## //TTbar To All
  } else if( (dataset=="TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="TToBLNu_TuneZ2_s-channel_7TeV-madgraph") ) {//## single top, s-channel
    xSection = 4.6;//##
  } else if( (dataset=="TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="TToBLNu_TuneZ2_t-channel_7TeV-madgraph") ) {//## single top, t-channel
    xSection = 62.8;//##
  } else if( (dataset=="TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="TToBLNu_TuneZ2_tW-channel_7TeV-madgraph") ) {//## single top, tW-channel
    xSection = 10.56;//##
  } else if( (dataset=="WWtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="WWtoAnything_TuneZ2_7TeV-pythia6-tauola") ) {//##  WW
    xSection = 42.9;//## //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt
  } else if( (dataset=="WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="WZtoAnything_TuneZ2_7TeV-pythia6-tauola") ) {//##  WZ
    xSection = 18.3;//## //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt
  } else if( (dataset=="ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1" ) || (dataset=="ZZtoAnything_TuneZ2_7TeV-pythia6-tauola") ) {//##  ZZ
    xSection = 5.9;//## //MCFM NLO see http://ceballos.web.cern.ch/ceballos/hwwlnln/cross_sections_backgrounds.txt
  } else if( (dataset=="WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1_2") || ( dataset=="WJetsToLNu_TuneZ2_7TeV-madgraph-tauola" ) ) {//## W+Jets
    xSection = 31314;//##
  } else if( (dataset=="DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="DYToEE_M-10To20_TuneZ2_7TeV-pythia6") ) {//## DY. EE 10to20
    xSection = 3457./3.;//##
  } else if( (dataset=="DYToEE_M-20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="DYToEE_M-20_TuneZ2_7TeV-pythia6") ) {//## DY. EE >20
    xSection = 4819.6/3.;//##
  } else if( (dataset=="DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6") ) {//## DY. MUMU 10to20
    xSection = 3457./3.;//##
  } else if( (dataset=="DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_2") || (dataset=="DYToMuMu_M-20_TuneZ2_7TeV-pythia6") ) {//## DY. MUMU >20
    xSection = 4819.6/3.;//##
  } else if( dataset=="QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6" ) {//## QCD_MUEnriched Pt20
    xSection =  296600000.*0.0002855 ;//##
  } else if( dataset=="QCD_Pt-20to30_EMEnriched_TuneZ2_7TeV-pythia6_3" ) {//## QCD_EMEnriched Pt20to30
    xSection = 236000000. * 0.0104;//##
  } else if( dataset=="QCD_Pt-30to80_EMEnriched_TuneZ2_7TeV-pythia6_3" ) {//## QCD_EMEnriched Pt30to80
    xSection = 59480000. * 0.065;//##
  } else if( dataset=="QCD_Pt-80to170_EMEnriched_TuneZ2_7TeV-pythia6_3" ) {//## QCD_EMEnriched Pt80to170
    xSection = 900000. * 0.155;//##
  } else if( dataset=="QCD_Pt-20to30_BCtoE_TuneZ2_7TeV-pythia6" ) {//## QCD_ Pt20to30
    xSection = 236000000. * 0.00056;//##
  } else if( dataset=="QCD_Pt-30to80_BCtoE_TuneZ2_7TeV-pythia6" ) {//## QCD_ Pt30to80
    xSection = 59480000. * 0.00230;//##
  } else if( dataset=="QCD_Pt-80to170_BCtoE_TuneZ2_7TeV-pythia6" ) {//## QCD_ Pt80to170
    xSection = 900000 * 0.0104;//##
  } else if( dataset=="G1Jet_Pt-120to180_TuneZ2_7TeV-alpgen" ) {//##  G1Jet_Pt-120to180
    xSection = 40.11;//##
  } else if( dataset=="G1Jet_Pt-180to240_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 5.612;//##
  } else if( dataset=="G1Jet_Pt-20to60_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 22700;//##
  } else if( dataset=="G1Jet_Pt-300to5000_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 0.5450;//##
  } else if( dataset=="G1Jet_Pt-60to120_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 799.3;//##
  } else if( dataset=="G2Jets_Pt-120to180_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 35.56;//##
  } else if( dataset=="G2Jets_Pt-180to240_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 5.941;//##
  } else if( dataset=="G2Jets_Pt-20to60_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 4080.;//##
  } else if( dataset=="G2Jets_Pt-240to300_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 1.448;//##
  } else if( dataset=="G2Jets_Pt-300to5000_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 0.7233;//##
  } else if( dataset=="G2Jets_Pt-60to120_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 425.8;//##
  } else if( dataset=="G3Jets_Pt-120to180_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 15.27;//##
  } else if( dataset=="G3Jets_Pt-180to240_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 3.09;//##
  } else if( dataset=="G3Jets_Pt-20to60_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 794.7;//##
  } else if( dataset=="G3Jets_Pt-240to300_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 0.7587;//##
  } else if( dataset=="G3Jets_Pt-300to5000_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 0.4663;//##
  } else if( dataset=="G3Jets_Pt-60to120_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 118.7;//##
  } else if( dataset=="G4Jets_Pt-120to180_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 5.732;//##
  } else if( dataset=="G4Jets_Pt-180to240_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 1.398;//##
  } else if( dataset=="G4Jets_Pt-20to60_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 148.3;//##
  } else if( dataset=="G4Jets_Pt-240to300_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 0.4199;//##
  } else if( dataset=="G4Jets_Pt-300to5000_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 0.2488;//##
  } else if( dataset=="G4Jets_Pt-60to120_TuneZ2_7TeV-alpgen" ) {//##
    xSection = 32.08;//##
  } else {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Dataset: '" << dataset << "' not present in database. Cross section unknown." << std::endl;
    std::cout << "-> Will set unitary weights." << std::endl;
    return 1.;
  }


  if( isAlpgenZJets ) {
    std::cout << "-> Scaling LO alpgen cross-section to NNLO." << std::endl;
    xSection*=1.27; // K factor
  }


  std::cout << "-> Total Events Analyzed: " << nEvents << ". " << std::endl;
  std::cout << "-> Dataset cross-section: " << xSection << " pb" << std::endl;
  float weight = xSection/((float)nEvents);
  std::cout << "=> Event weight: " << weight << std::endl;

  return weight;

}
