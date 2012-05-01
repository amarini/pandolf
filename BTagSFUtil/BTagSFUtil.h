#ifndef BTagSFUtil_h
#define BTagSFUtil_h

//#include <Riostream.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TH2D.h"
#include <iostream>
#include "SFlightFuncs.h"
#include "MistagFuncs.h"




// ............................................................................................................
/*
Usage:  

In the analysis macro:

BTagSFUtil* btsfutil = new BTagSFUtil(13);
// 13 is just the seed

btsfutil->modifyBTagsWithSF(taggerName, jet1_tagged_loose, jet1_tagged_medium, Jet1pt_, Jet1eta_, jet1flav_ );

// taggerName can be: "TCHE", "CSV", "JP"

*/
// .............................................................................................................


struct BTagScaleFactor{

  float SF;
  float SF_err;
  float eff;
  float eff_err;

};



class BTagSFUtil {

 public:

  BTagSFUtil( int seed=0 );
  BTagSFUtil( const std::string& btagAlgo, int seed=0 );

  void init( const std::string& wp );
  void InitSFLight( const std::string& wp );
  void InitMistag( const std::string& wp );
  void InitSFb( const std::string& wp );
  void checkInit( const std::string& wp );

  float GetSFLight( float pt, float eta, const std::string& wp, const std::string& meanminmax="mean" );
  float GetMistag( float pt, float eta, const std::string& wp, const std::string& meanminmax="mean" );
  float GetSFb( float pt, float eta, const std::string& wp );

  TF1* GetFunctionEtaBins( float eta, const std::vector<float>& etaBins, const std::vector<TF1*>& functions ) const;


/***********************************/
/* UP TO DATE FUNCTIONS*/
/***********************************/
  void SF(const std::string& btagAlgo, const std::string& wp, float pt = 60., float eta = 1.6);
  void modifyBTagsWithSF_fast( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, const std::string& meanMinMax="mean" );
  void modifyBTagsWithSF(const std::string& btagAlgo , bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, float sysSF=1.0);


/***********************************/
/* FUNCTIONS NOT UP TO DATE */
/***********************************/
  void set_fileMedium( TFile* file );
  void set_fileLoose( TFile* file );
  void modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, int pdgIdPart, float Btageff_SF_l = 0.95,float Btageff_SF_m = 0.94, float Btagmistag_SF_l = 1.0,float Btagmistag_SF_m = 1.0, float Btagmistag_eff_l = 1.0, float  Btagmistag_eff_m = 1.0);

  // deprecated:
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt, float jeteta );
  BTagScaleFactor  getSF( const std::string& fileName, float jetpt );
  void setSFFileName(const std::string fileName){ sfFileName_= fileName;};

 private:

  std::string btagAlgo_;

  bool init_loose_;
  bool init_medium_;
  bool init_tight_;

  float SFb_;
  float SFlight_, Mistag_;

  std::vector<float> etaBins_SFLight_loose_;
  std::vector<TF1*> sfLightFunct_loose_;
  std::vector<TF1*> sfLightFunct_loose_min_;
  std::vector<TF1*> sfLightFunct_loose_max_;

  std::vector<float> etaBins_SFLight_medium_;
  std::vector<TF1*> sfLightFunct_medium_;
  std::vector<TF1*> sfLightFunct_medium_min_;
  std::vector<TF1*> sfLightFunct_medium_max_;

  std::vector<float> etaBins_SFLight_tight_;
  std::vector<TF1*> sfLightFunct_tight_;
  std::vector<TF1*> sfLightFunct_tight_min_;
  std::vector<TF1*> sfLightFunct_tight_max_;

  std::vector<float> etaBins_mistag_loose_;
  std::vector<TF1*> mistagFunct_loose_;

  std::vector<float> etaBins_mistag_medium_;
  std::vector<TF1*> mistagFunct_medium_;

  std::vector<float> etaBins_mistag_tight_;
  std::vector<TF1*> mistagFunct_tight_;

  TF1* SFbFunct_loose_;
  TF1* SFbFunct_medium_;


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
