#ifndef BTagSFUtil_h
#define BTagSFUtil_h

//#include <Riostream.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TH2D.h"


struct BTagScaleFactor{

  float SF;
  float SF_err;
  float eff;
  float eff_err;

};



class BTagSFUtil {

 public:

  BTagSFUtil( int seed=0 );

  void set_fileMedium( TFile* file );
  void set_fileLoose( TFile* file );

  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, float sysSF=1.0);
  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, int pdgIdPart, float Btageff_SF_l = 0.95,float Btageff_SF_m = 0.94, float Btagmistag_SF_l = 1.0,float Btagmistag_SF_m = 1.0, float Btagmistag_eff_l = 1.0, float  Btagmistag_eff_m = 1.0);

  // deprecated:
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt, float jeteta );
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt );
  void setSFFileName(const std::string fileName){ sfFileName_= fileName;};

 private:

  TRandom3* rand_;
  std::string sfFileName_;

  TFile* fileMedium_;
  TFile* fileLoose_;

  TH2D* h2_Loose_BTAGBEFFCORR_;
  TH2D* h2_Loose_BTAGLEFFCORR_;
  TH2D* h2_Loose_BTAGLEFF_;

  TH2D* h2_Medium_BTAGBEFFCORR_;
  TH2D* h2_Medium_BTAGLEFFCORR_;
  TH2D* h2_Medium_BTAGLEFF_;

};

#endif
