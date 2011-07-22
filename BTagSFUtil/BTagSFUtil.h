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

  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, double Btageff_SF = 0.95, double Btagmistag_SF = 1.0, const std::string& tagger="TCHE");
  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, int pdgIdPart, float Btageff_SF_l = 0.95,float Btageff_SF_m = 0.94, float Btagmistag_SF_l = 1.0,float Btagmistag_SF_m = 1.0, float Btagmistag_eff_l = 1.0, float  Btagmistag_eff_m = 1.0);
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt, float jeteta );
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt );
  void setSFFileName(const std::string fileName){ sfFileName_= fileName;};

 private:

  TRandom3* rand_;
  std::string sfFileName_;

};

#endif
