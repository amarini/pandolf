#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>
#include <cstdlib>


TChain* tree = 0;

struct EventsAndLumi {
  int nTotalEvents;
  float totalLumi;
};



EventsAndLumi addInput( const std::string& dataset, const std::string& algoName, const std::string& flags );
float getWeight( const std::string& dataset, int nEvents );


int main( int argc, char* argv[] ) {

  if( argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./merge_and_setWeights [dataset] [recoType] [jetAlgo] [flags=\"\"]" << std::endl;
    exit(917);
  }

  std::string dataset = argv[1];
  std::string recoType = argv[2];
  std::string jetAlgo = argv[3];

  std::string flags="";
  if( argc==5 ) {
    std::string flags_str(argv[4]);
    flags = flags_str;
  }

  std::string algoName = recoType+jetAlgo;
  if( recoType=="calo") algoName=jetAlgo;
  if( recoType=="jpt" && jetAlgo=="akt5" ) algoName="jptak5";
  if( recoType=="jpt" && jetAlgo=="akt7" ) algoName="jptak7";
  
  tree = new TChain("jetTree");

  EventsAndLumi evlu;
  evlu = addInput( dataset, algoName, jetAlgo );

  float weight = getWeight( dataset, evlu.nTotalEvents );

  // and now set the weights
  tree->SetBranchStatus( "eventWeight", 0 );

  std::string outfilename;
  if( flags!="" )
    outfilename = "PhotonJet_2ndLevelTreeW_"+dataset+"_"+algoName+"_"+flags+".root";
  else
    outfilename = "PhotonJet_2ndLevelTreeW_"+dataset+"_"+algoName+".root";

  TFile* outfile = new TFile(outfilename.c_str(), "recreate");
  outfile->cd();

  TH1F* h1_lumi = new TH1F("lumi", "", 1, 0., 1.);
  h1_lumi->SetBinContent(1, evlu.totalLumi);

  TTree* newTree = tree->CloneTree(0);
  Float_t newWeight;
  newTree->Branch( "eventWeight", &newWeight, "newWeight/F" );

  int nentries = tree->GetEntries();
  for( unsigned ientry = 0; ientry<nentries; ++ientry ) {

    tree->GetEntry(ientry);

    if( (ientry % 10000) ==0 ) std::cout << "Entry: " << ientry << " /" << nentries << std::endl;

    newWeight = weight;

    newTree->Fill();

  } //for entries

  h1_lumi->Write();
  newTree->Write();
  outfile->Write();
  outfile->Close();

  return 0;

}


EventsAndLumi addInput( const std::string& dataset, const std::string& algoName, const std::string& flags ) {

  std::string infileName;
  if( flags!="" )
    infileName = "files_PhotonJet_2ndLevel_" + dataset + "_" + algoName + ".txt";
  else
    infileName = "files_PhotonJet_2ndLevel_" + dataset + "_" + algoName + "_" + flags + ".txt";
  
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
      infile->Close();

    }
    fclose(iff);

  }

  EventsAndLumi evlu;
  evlu.nTotalEvents = totalEvents;
  evlu.totalLumi = totalLumi;

  return evlu;

} //addinput


float getWeight( const std::string& dataset, int nEvents ) {

  float xSection = -1.;

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
  } else {
    std::cout << std::endl;
    std::cout << "-> WARNING!! Dataset: '" << dataset << "' not present in database. Cross section unknown." << std::endl;
    std::cout << "-> Will set unitary weights." << std::endl;
    return 1.;
  }

  float weight = xSection/((float)nEvents);

  return weight;

}
