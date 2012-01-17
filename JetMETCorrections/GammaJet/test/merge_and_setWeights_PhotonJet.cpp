#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <cstdlib>


TChain* tree = 0;
TH1F* h1_eff_vs_pt = 0;

struct EventsAndLumi {
  int nTotalEvents;
  float totalLumi;
};



EventsAndLumi addInput( const std::string& analyzerType, const std::string& dataset, const std::string& algoName, const std::string& flags );
float getWeight( const std::string& dataset, int nEvents );


int main( int argc, char* argv[] ) {

  if( argc!=4 && argc!=5 && argc!=6 ) {
    std::cout << "USAGE: ./merge_and_setWeights_PhotonJet [dataset] [recoType] [jetAlgo] [analyzerType=\"PhotonJet\"] [flags=\"\"]" << std::endl;
    exit(917);
  }

  std::string dataset = argv[1];
  std::string recoType = argv[2];
  std::string jetAlgo = argv[3];

  std::string analyzerType="PhotonJet";
  if( argc==5 ) {
    std::string analyzerType_str(argv[4]);
    analyzerType = analyzerType_str;
  }

  std::string flags="";
  if( argc==6 ) {
    std::string flags_str(argv[5]);
    flags = flags_str;
  }

  std::string algoName = recoType+jetAlgo;
  if( recoType=="calo") algoName=jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5" ) algoName="jptak5";
  if( recoType=="jpt" && jetAlgo=="akt7" ) algoName="jptak7";
  
  tree = new TChain("jetTree");

  EventsAndLumi evlu;
  evlu = addInput( analyzerType, dataset, algoName, flags );

  float weight = getWeight( dataset, evlu.nTotalEvents );

  // and now set the weights
  tree->SetBranchStatus( "eventWeight", 0 );

  Float_t ptPhotReco;
  tree->SetBranchAddress( "ptPhotReco", &ptPhotReco );

  std::string outfilename;
  if( flags!="" )
    outfilename = analyzerType + "_2ndLevelTreeW_"+dataset+"_"+algoName+"_"+flags+".root";
  else
    outfilename = analyzerType + "_2ndLevelTreeW_"+dataset+"_"+algoName+".root";

  TFile* outfile = new TFile(outfilename.c_str(), "recreate");
  outfile->cd();

  TH1F* h1_lumi = new TH1F("lumi", "", 1, 0., 1.);
  h1_lumi->SetBinContent(1, evlu.totalLumi);

  TTree* newTree = tree->CloneTree(0);
  newTree->SetMaxTreeSize(100000000000ULL); //setting max tree size to 100 GB
  Float_t newWeight;
  newTree->Branch( "eventWeight", &newWeight, "newWeight/F" );

  int nentries = tree->GetEntries();
  for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( (ientry % 100000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

    newWeight = weight;

    if( flags=="GENJETS" && h1_eff_vs_pt!=0 ) {
      Int_t iBin = h1_eff_vs_pt->FindBin( ptPhotReco );
      newWeight *= h1_eff_vs_pt->GetBinContent(iBin);
    }

    newTree->Fill();

  } //for entries

  h1_lumi->Write();
  newTree->Write();
  outfile->Write();
  if( h1_eff_vs_pt!= 0 )
    h1_eff_vs_pt->Write();

  outfile->Close();

  return 0;

}


EventsAndLumi addInput( const std::string& analyzerType, const std::string& dataset, const std::string& algoName, const std::string& flags ) {

  TH1::AddDirectory(kFALSE);

  std::string infileName;
  if( flags!="" )
    infileName = "files_" + analyzerType + "_2ndLevel_" + dataset + "_" + algoName + "_" + flags + ".txt";
  else
    infileName = "files_" + analyzerType + "_2ndLevel_" + dataset + "_" + algoName + ".txt";
  
  TH1F* h1_lumi;
  TH1F* h1_nCounter;

  int totalEvents = 0;
  float totalLumi = 0.;

  //open from file.txt:
  FILE* iff = fopen(infileName.c_str(),"r");
  if(iff == 0) {
    std::cout << "cannot open input file '" << infileName << "' ... adding single file." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9172);
//  infileName = "HZZlljj_2ndLevelTree_" + dataset + ".root";
//  std::string treeName = infileName +"/jetTree";
//  tree->Add(treeName.c_str());
//  std::cout << "-> Added " << treeName << ". Tree has " << tree->GetEntries() << " entries." << std::endl;
//  TFile* infile = TFile::Open(infileName.c_str(), "READ");
//  h1_nCounter = (TH1F*)infile->Get("nCounter");
//  if( h1_nCounter!=0 ) {
//    totalEvents += h1_nCounter->GetBinContent(1);
//  } else {
//    std::cout << " WARNING! File '" << infileName << "' has no nCounter information. Skipping." << std::endl;
//  }
//  h1_lumi = (TH1F*)infile->Get("lumi");
//  if( h1_lumi!=0 ) {
//    totalLumi += h1_lumi->GetBinContent(1);
//    std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
//  } else {
//    std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
//  }
//  infile->Close();

  } else {

    char singleLine[500];

    std::cout << "-> Opened file: '" << infileName << std::endl;
    bool firstFile = true;
    TH1F* h1_eff_num_vs_pt=0;
    TH1F* h1_eff_denom_vs_pt=0;

    while( fscanf(iff, "%s", singleLine) !=EOF ) {

      std::string rootfilename(singleLine);
      std::string treename = rootfilename + "/jetTree";
      std::cout << "-> Added " << treename;
      tree->Add(treename.c_str());
      TFile* infile = TFile::Open(rootfilename.c_str(), "READ");
      h1_nCounter = (TH1F*)infile->Get("nCounter");
      if( h1_nCounter!=0 ) {
        totalEvents += h1_nCounter->GetBinContent(1);
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no nCounter information. Skipping." << std::endl;
      }
      h1_lumi = (TH1F*)infile->Get("lumi");
      if( h1_lumi!=0 ) {
        totalLumi += h1_lumi->GetBinContent(1);
        std::cout << "\tTotal lumi: " << totalLumi << " ub-1" << std::endl;
      } else {
        std::cout << " WARNING! File '" << infileName << "' has no lumi information. Skipping." << std::endl;
      }

      TH1F* h1_eff_num_vs_pt_tmp = (TH1F*)infile->Get("eff_num_medium_vs_pt"); //medium is default
      TH1F* h1_eff_denom_vs_pt_tmp = (TH1F*)infile->Get("eff_denom_vs_pt"); //medium is default
      if( h1_eff_num_vs_pt_tmp!=0 && h1_eff_denom_vs_pt_tmp!=0 ) { //genjets case
        if( firstFile ) {
          h1_eff_num_vs_pt = new TH1F(*h1_eff_num_vs_pt_tmp);
          h1_eff_denom_vs_pt = new TH1F(*h1_eff_denom_vs_pt_tmp);
          firstFile = false;
        } else {
          h1_eff_num_vs_pt->Add(h1_eff_num_vs_pt_tmp);
          h1_eff_denom_vs_pt->Add(h1_eff_denom_vs_pt_tmp);
        }
      } else {
        if( flags=="GENJETS" ) std::cout << "GENJETS but didnt find eff histos!!" << std::endl;
      }
      infile->Close();

    }

    fclose(iff);

    if( h1_eff_num_vs_pt!=0 && h1_eff_denom_vs_pt!=0 ) {
      h1_eff_vs_pt = new TH1F(*h1_eff_num_vs_pt);
      h1_eff_vs_pt->Divide(h1_eff_denom_vs_pt);
    }

  }

  EventsAndLumi evlu;
  evlu.nTotalEvents = totalEvents;
  evlu.totalLumi = totalLumi;

  return evlu;

} //addinput


float getWeight( const std::string& dataset, int nEvents ) {

  float xSection = -1.;
  TString dataset_tstr(dataset);

  if( dataset=="Zjets-madgraph" ) {
    xSection = 2800.;
  } else if( dataset=="HZZ_qqll_gluonfusion_M130" ) {
    xSection = 25.560*0.03913; //sigma x BR
  } else if( dataset=="HZZ_qqll_gluonfusion_M150" ) {
    xSection = 19.568*0.08234; //sigma x BR
  } else if( dataset=="HZZ_qqll_gluonfusion_M500" ) {
    xSection = 2.1914*0.2602; //sigma x BR
  } else if( dataset=="QCD_Spring10_Pt0to15" ) {
    xSection = 48445000000.;
  } else if( dataset=="QCD_Spring10_Pt5to15" ) {
    xSection = 36640000000.;
  } else if( dataset=="QCD_Spring10_Pt15" ) {
    xSection = 876215000.-60411000.;
  } else if( dataset=="QCD_Spring10_Pt15to20" ) {
    xSection = 579411000.;
  } else if( dataset=="QCD_Spring10_Pt20to30" ) {
    xSection = 236051000.;
  } else if( dataset=="QCD_Spring10_Pt30" ) {
    xSection = 60411000.-923821.;
  } else if( dataset=="QCD_Spring10_Pt30to50" ) {
    xSection = 53114800.;
  } else if( dataset=="QCD_Spring10_Pt50to80" ) {
    xSection = 6358210.;
  } else if( dataset=="QCD_Spring10_Pt80" ) {
    xSection = 923821.-25474.9;
  } else if( dataset=="QCD_Spring10_Pt170" ) {
    xSection = 25474.9-1255.87;
  } else if( dataset=="QCD_Spring10_Pt300" ) {
    xSection = 1255.87-87.9799;
  } else if( dataset=="QCD_Spring10_Pt470" ) {
    xSection = 87.9799-2.18608;
  } else if( dataset=="QCD_Spring10_Pt800" ) {
    xSection = 2.18608-0.0112233;
  } else if( dataset=="QCD_Spring10_Pt1400" ) {
    xSection = 0.0112233;
  } else if( dataset_tstr.BeginsWith("G_Pt_0to15") ) {
    xSection = 8.420e+07;
  //} else if( dataset=="G_Pt_15to30_TuneZ2_7TeV_pythia6" || dataset=="G_Pt_15to30_TuneZ2_7TeV_pythia6_CORR" ) {
  } else if( dataset_tstr.BeginsWith("G_Pt_15to30") || dataset_tstr.BeginsWith("G_Pt-15to30") ) {
    xSection = 1.717e+05;
  } else if( dataset_tstr.BeginsWith("G_Pt_30to50") || dataset_tstr.BeginsWith("G_Pt-30to50") ) {
    xSection = 1.669e+04;
  } else if( dataset_tstr.BeginsWith("G_Pt_50to80") || dataset_tstr.BeginsWith("G_Pt-50to80") ) {
    xSection = 2.722e+03;
  } else if( dataset_tstr.BeginsWith("G_Pt_80to120") || dataset_tstr.BeginsWith("G_Pt-80to120") ) {
    xSection = 4.472e+02;
  } else if( dataset_tstr.BeginsWith("G_Pt_120to170") || dataset_tstr.BeginsWith("G_Pt-120to170") ) {
    xSection = 8.417e+01;
  } else if( dataset_tstr.BeginsWith("G_Pt_170to300") || dataset_tstr.BeginsWith("G_Pt-170to300") ) {
    xSection = 2.264e+01;
  } else if( dataset_tstr.BeginsWith("G_Pt_300to470") || dataset_tstr.BeginsWith("G_Pt-300to470") ) {
    xSection = 1.493e+00;
  } else if( dataset_tstr.BeginsWith("G_Pt_470to800") || dataset_tstr.BeginsWith("G_Pt-470to800") ) {
    xSection = 1.323e-01;
  } else if( dataset_tstr.BeginsWith("G_Pt_800to1400") || dataset_tstr.BeginsWith("G_Pt-800to1400") ) {
    xSection = 3.481e-03;
  } else if( dataset_tstr.BeginsWith("G_Pt_1400to1800") || dataset_tstr.BeginsWith("G_Pt-1400to1800") ) {
    xSection = 1.270e-05;
  } else if( dataset_tstr.BeginsWith("G_Pt_1800") || dataset_tstr.BeginsWith("G_Pt-1800") ) {
    xSection = 2.936e-07;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_15to30") || dataset_tstr.BeginsWith("QCD_Pt-15to30") ) {
    xSection = 8.159e+08;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_30to50") || dataset_tstr.BeginsWith("QCD_Pt-30to50")) {
    xSection = 5.312e+07;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_50to80") || dataset_tstr.BeginsWith("QCD_Pt-50to80")) {
    xSection = 6.359e+06;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_80to120") || dataset_tstr.BeginsWith("QCD_Pt-80to120")) {
    xSection = 7.843e+05;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_120to170") || dataset_tstr.BeginsWith("QCD_Pt-120to170")) {
    xSection = 1.151e+05;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_170to300") || dataset_tstr.BeginsWith("QCD_Pt-170to300")) {
    xSection = 2.426e+04;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_300to470") || dataset_tstr.BeginsWith("QCD_Pt-300to470")) {
    xSection = 1.168e+03;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_470to600") || dataset_tstr.BeginsWith("QCD_Pt-470to600")) {
    xSection = 7.022e+01;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_600to800") || dataset_tstr.BeginsWith("QCD_Pt-600to800")) {
    xSection = 1.555e+01;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_800to1000") || dataset_tstr.BeginsWith("QCD_Pt-800to1000")) {
    xSection = 1.844e+00;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_1000to1400") || dataset_tstr.BeginsWith("QCD_Pt-1000to1400")) {
    xSection = 3.321e-01;
  } else if( dataset_tstr.BeginsWith("QCD_Pt_1400to1800") || dataset_tstr.BeginsWith("QCD_Pt-1400to1800")) {
    xSection = 1.087e-02;
  } else if( dataset=="PhotonJet_Summer1036X_Pt5to15" ) {
    xSection = 4030000.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt15to20" ) {
    xSection = 114700.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt20to30" ) {
    xSection = 57180.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt30to50" ) {
    xSection = 16520.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt50to80" ) {
    xSection = 2723.;
  } else if( dataset=="PhotonJet_Summer1036X_Pt80to120" ) {
    xSection = 446.2;
  } else if( dataset=="PhotonJet_Summer1036X_Pt120to170" ) {
    xSection = 84.43;
  } else if( dataset=="PhotonJet_Summer1036X_Pt170to300" ) {
    xSection = 22.55;
  } else if( dataset=="PhotonJet_Summer1036X_Pt300to500" ) {
    xSection = 1.545;
  } else if( dataset=="PhotonJet_Summer1036X_Pt500toInf" ) {
    xSection = 0.0923;
  } else if( dataset_tstr.BeginsWith("QCD_Pt-20to30_BCtoE") ) {
    xSection = 236000000.*0.00056;
  } else if( dataset_tstr.BeginsWith("QCD_Pt-20to30_EMEnriched") ) {
    xSection = 236000000.*0.0104;
  } else if( dataset_tstr.BeginsWith("QCD_Pt-30to80_BCtoE") ) {
    xSection = 59480000.*0.00230;
  } else if( dataset_tstr.BeginsWith("QCD_Pt-30to80_EMEnriched") ) {
    xSection = 59480000.*0.065;
  } else if( dataset_tstr.BeginsWith("QCD_Pt-80to170_BCtoE") ) {
    xSection = 900000.*0.0104;
  } else if( dataset_tstr.BeginsWith("QCD_Pt-80to170_EMEnriched") ) {
    xSection = 900000.*0.155;
  } else if( dataset_tstr.BeginsWith("WH_ZH_HToGG_M-120") ) {
    xSection = (0.6561 + 0.3551)*2.25E-03;
  } else if( dataset_tstr.BeginsWith("WH_ZH_HToGG_M-125") ) {
    xSection = (0.5729 + 0.3158)*2.25E-03;
  } else if( dataset_tstr.BeginsWith("DiPhotonJets") ) {
    xSection = 134 * 1.15;

  } else if( dataset_tstr.BeginsWith("G1Jet_Pt-20to60") ) {
    xSection = 2.270e+04;
  } else if( dataset_tstr.BeginsWith("G1Jet_Pt-60to120") ) {
    xSection = 7.993e+02;
  } else if( dataset_tstr.BeginsWith("G1Jet_Pt-120to180") ) {
    xSection = 4.011e+01 ;
  } else if( dataset_tstr.BeginsWith("G1Jet_Pt-180to240") ) {
    xSection = 5.612e+00 ;
  } else if( dataset_tstr.BeginsWith("G1Jet_Pt-240to300") ) {
    xSection = 1.224e+00 ;
  } else if( dataset_tstr.BeginsWith("G1Jet_Pt-300to5000") ) {
    xSection = 5.450e-01 ;


  } else if( dataset_tstr.BeginsWith("G2Jet_Pt-20to60") ) {
    xSection = 4.080e+03;
  } else if( dataset_tstr.BeginsWith("G2Jet_Pt-60to120") ) {
    xSection = 4.258e+02 ;
  } else if( dataset_tstr.BeginsWith("G2Jet_Pt-120to180") ) {
    xSection = 3.556e+01 ;
  } else if( dataset_tstr.BeginsWith("G2Jet_Pt-180to240") ) {
    xSection = 5.941e+00 ;
  } else if( dataset_tstr.BeginsWith("G2Jet_Pt-240to300") ) {
    xSection = 1.448e+00 ;
  } else if( dataset_tstr.BeginsWith("G2Jet_Pt-300to5000") ) {
    xSection = 7.233e-01 ;

  } else {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Dataset: '" << dataset << "' not present in database. Cross section unknown." << std::endl;
    std::cout << "-> Will set unitary weights." << std::endl;
    return 1.;
  }

  float weight = xSection/((float)nEvents);

  return weight;

}
