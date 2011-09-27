#include "TFile.h"
#include <cstdlib>
#include <string>
#include "TTree.h"
#include "TH1F.h"

// the following script converts the bool, isSidebands, into an int, isSB, and 
//rewrites the tree.  It also, in this version, adds a wght according to alpha(mZZ)
// NOTE: this was meant to use on DATA sidebands only!!!

void readWriteTree(string fileName, string treeName){

  TFile* f = new TFile(fileName.c_str());
  TTree* t = (TTree*) f->Get(treeName.c_str());

  TFile* f_alpha = new TFile("alpha.root");
  TH1F *alpha = (TH1F*) f_alpha->Get("alpha_MADGRAPH");

  float mWW;
  float mJJ;
  float wght;
  bool isSidebands;
  int isSB;

  t->SetBranchAddress("mWW",&mWW);
  t->SetBranchAddress("mJJ",&mJJ);
  t->SetBranchAddress("eventWeight",&wght);
  
  string tempFile = "NEW_"+fileName;
  TFile *outFile  = new TFile(tempFile.c_str(),"RECREATE");
  TTree* selectedEvents = new TTree("selectedEvents","selectedEvents");

  selectedEvents->Branch("mWW",&mWW);
  selectedEvents->Branch("mJJ",&mJJ);
  selectedEvents->Branch("isSB",&isSB);
  selectedEvents->Branch("wght",&wght);
  
  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);

    if( (mJJ>40. && mJJ<60.)||(mJJ>100.&&mJJ<160.) ) isSB=1;
    else isSB=0;

    if(!(mJJ>75 && mJJ<105)){

      // 150 and 10 should not be hard coded... get from histograms!

	wght=alpha->GetBinContent(alpha->FindBin(mWW));
    }

    selectedEvents->Fill();

  }
  outFile->cd();
  selectedEvents->Write();

}
