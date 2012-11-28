#include "TreeFinalizer.h"
#include <iostream>
#include "TROOT.h"
#include <cstdlib>



TreeFinalizer::TreeFinalizer( const std::string& analyzerType, const std::string& dataset, const std::string& flags ) {

  DEBUG_ = false;


  tree_ = new TChain("reducedTree");

  analyzerType_ = analyzerType;
  inputAnalyzerType_ = analyzerType;
  inputFileDir_ = "";
  dataset_ = dataset;
  flags_ = flags;

  outFile_ = 0;

  nCounter_ = 0.;
  nCounterW_ = 0.;
  nCounterPU_ = 0.;

} //constructor



TreeFinalizer::~TreeFinalizer() {

  std::cout << std::endl << "-> Histograms saved in: " << outFile_->GetName() << std::endl;

  if( tree_!=0 ) {
    delete tree_;
    tree_=0;
  }

} //destructor



void TreeFinalizer::createOutputFile( const std::string& additionalFlags ) {

   std::string outfileName;

   if( DEBUG_ ) outfileName = "prova_"+dataset_;
   else {
    if(dataset_!="") outfileName = analyzerType_ + "_" + dataset_;
    else outfileName = analyzerType_;
   }


   if( flags_!="" )
     outfileName = outfileName + "_" + flags_;
   if( additionalFlags!="" )
     outfileName = outfileName + "_" + additionalFlags;

   outfileName = outfileName ;
   outfileName = outfileName + ".root";

   outFile_ = TFile::Open(outfileName.c_str(), "RECREATE");
   
   outFile_->cd();


}





void TreeFinalizer::addInput( const std::string& dataset ) {

    std::string infileName = inputAnalyzerType_ + "_2ndLevelTreeW_" + dataset;
    if( inputFileDir_ != "" )
      infileName = inputFileDir_ + "/" + infileName;
    infileName += ".root";
    TFile* infile = TFile::Open(infileName.c_str(), "read");
    if( infile==0 ) {
      std::cout << "Didn't find file '" << infileName << "'. Did you forget to finalize (i.e. the \"W\")?" << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(77);
    }
    std::string treeName = infileName +"/reducedTree";
    tree_->Add(treeName.c_str());
    std::cout << "-> Added " << treeName << ". Tree has " << tree_->GetEntries() << " entries." << std::endl;
    infile->Close();


} //addinput





std::vector<TH1F*> TreeFinalizer::getResponseHistos(const std::string& name, unsigned binArraySize, Double_t* ptBins) {

  std::vector<TH1F*> returnVector;

  for( unsigned i=0; i<(binArraySize-1); ++i ) {
    char histoName[100];
    sprintf( histoName, "%s_ptBin_%.0f_%.0f", name.c_str(), ptBins[i], ptBins[i+1]);
    int nbins = 50;
    float xmin = 0.5;
    float xmax = 2.;
    TH1F* newHisto = new TH1F(histoName, "", nbins, xmin, xmax);
    newHisto->Sumw2();
    returnVector.push_back(newHisto);
  }

  return returnVector;

}


void TreeFinalizer::writeResponseHistos( TFile* file, std::vector<TH1F*> h1_response, std::string dirName ) {

  file->mkdir( dirName.c_str() );
  file->cd( dirName.c_str() );

  for( unsigned iHisto=0; iHisto<h1_response.size(); ++iHisto ) h1_response[iHisto]->Write();

  file->cd();

}



