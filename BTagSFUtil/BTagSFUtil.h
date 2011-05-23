#ifndef BTagSFUtil_h
#define BTagSFUtil_h

#include <Riostream.h>
#include "TRandom3.h"


struct BTagScaleFactor{

  float SF;
  float SF_err;
  float eff;
  float eff_err;

};



class BTagSFUtil {

 public:

  BTagSFUtil( int seed=0 );

  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, const std::string& tagger="TCHE" );
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt, float jeteta );
  void setSFFileName(const std::string fileName){ sfFileName_= fileName;std::cout<<"New SF Filename: "<<sfFileName_.c_str()<<std::endl;};

 private:

  TRandom3* rand_;
  std::string sfFileName_;

};

#endif
