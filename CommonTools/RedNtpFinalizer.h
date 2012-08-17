// ------------------------------------------------------------
//  
//    RedNtpFinalizer - Base class for finalizing the analysis
//    Used in the Ntp1 tree workflow, after a Ntp1Analyzer
//    class has produced 2nd level trees 
//    (and merge_and_setWeights has done its job)
//
//    It is an abstract class: the method 'finalize' has to be
//    implemented.
//
// ------------------------------------------------------------


#include <vector>
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "AnalysisJet.h"
#include "BTagSFUtil/interface/BTagSFUtil.h"



class RedNtpFinalizer {


 public:

  RedNtpFinalizer( const std::string& analyzerType, const std::string& dataset, const std::string& flags="" );
  virtual ~RedNtpFinalizer();

  void createOutputFile( const std::string& additionalFlags="" );
  virtual void addFile(const std::string& dataseti, const std::string& selection="");

  TChain* get_tree() { return tree_; };
  TFile* get_outFile() { return outFile_; };
  bool get_DEBUG() { return DEBUG_; };
  int get_nBTags( const AnalysisJet& jet1, const AnalysisJet& jet2, BTagSFUtil* btsfutil, bool losebtags=true );
  float get_crossSection(const std::string& dataset) const;

  void clear();

  void set_outFile( const std::string& fileName="", const std::string& suffix="" );
  void set_dataset( const std::string& dataset ) { dataset_ = dataset; };
  void set_inputAnalyzerType( const std::string& analyzerType ) { inputAnalyzerType_ = analyzerType; };
  void set_flags( const std::string& flags ) { flags_ = flags; };
  void set_DEBUG( bool DEBUG ) { DEBUG_ = DEBUG; };

  virtual void finalize() = 0;


  TChain* tree_;

  float nGenEvents_;
  float xSection_;

  TH1F* h1_nPU_gen_;

  std::string analyzerType_;
  //std::string inputAnalyzerType_;
  std::string dataset_;
  std::string flags_;

  TFile* outFile_;

  bool DEBUG_;

 private:

};
