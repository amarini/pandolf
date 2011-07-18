#ifndef BTagSFUtil_h
#define BTagSFUtil_h

#include <Riostream.h>
#include "TRandom3.h"

#include <memory>
#include <string>
#include <map>

struct BTagScaleFactor{

  float etamin;
  float etamax;
  float ptmin;
  float ptmax;

  float SF;
  float SF_err;
  float eff;
  float eff_err;

};



class BTagSFUtil {

 public:

  BTagSFUtil( int seed=0 );

  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, double Btageff_SF = 0.95, double Btagmistag_SF = 1.0, const std::string& tagger="TCHE", bool verbose = true);
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt, float jeteta , bool verbose);
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt );
  void setSFFileName(const std::string fileName){ sfFileName_= fileName;};

 private:
  
  TRandom3* rand_;
  std::string sfFileName_;

  std::vector<BTagScaleFactor> cache_[5];
  std::map<std::string,int> cacheAssoc_;
  int cachesUsed_;

};

#endif
