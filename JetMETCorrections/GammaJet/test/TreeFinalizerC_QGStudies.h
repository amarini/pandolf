
#include "TreeFinalizer.h"
#include "AnalysisJet.h"



class TreeFinalizerC_QGStudies : public TreeFinalizer {

 public:

  TreeFinalizerC_QGStudies( const std::string& dataset, float secondJetThreshold=0.2, const std::string& recoType="pf", const std::string& jetAlgo="akt5", int iBlock=1, int nBlocks=1 ) : TreeFinalizer("QGStudies", dataset) {
    dataset_ = dataset;
    secondJetThreshold_ = secondJetThreshold;
    recoType_ = recoType;
    jetAlgo_ = jetAlgo;
    iBlock_ = iBlock;
    nBlocks_ = nBlocks;
  }
  virtual ~TreeFinalizerC_QGStudies() {};

  virtual void finalize();


 private:

   float secondJetThreshold_;
   std::string recoType_;
   std::string jetAlgo_;

   int iBlock_;
   int nBlocks_;

};


