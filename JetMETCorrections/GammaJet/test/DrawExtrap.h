// derived class from DrawBase
// draws the response and resolution extrapolations
// to pt (2nd Jet) -> 0

#ifndef DrawExtrap_h
#define DrawExtrap_h

#include "DrawBase.h"
#include "TFile.h"




class DrawExtrap : public DrawBase {

 public:

  DrawExtrap( const std::string& analysisType, const std::string& recoType="", const std::string& jetAlgo="", const std::string& flags="" );
  virtual ~DrawExtrap() {};
  
  void set_FITRMS( const std::string& fit_rms ) { FIT_RMS_ = fit_rms; };
  void set_NOQ( bool noq ) { NOQ_ = noq; };
  void set_INTPERC( float intperc ) { INTPERC_ = intperc; };
  void set_FIXM( bool fixm ) { FIXM_ = fixm; };

  void drawResponseExtrap(const std::string& etaRegion, bool corrected=false, const std::string& recoGen="RecoRel");


 private:

  void drawExtrapSinglePtBin(int iPtBin, Float_t& DATAReso, Float_t& DATAReso_err, Float_t& extrapReso, Float_t& extrapReso_err, Float_t& trueReso, Float_t& trueReso_err, Float_t& intrReso, Float_t& intrReso_err,
                                         Float_t& DATAResp, Float_t& DATAResp_err, Float_t& extrapResp, Float_t& extrapResp_err, Float_t& trueResp, Float_t& trueResp_err, Float_t& intrResp, Float_t& intrResp_err, 
                                         Float_t& imbalanceResp, Float_t& imbalanceResp_err, const std::string& recoGen) const;
  void getXPoints( TFile * file, const char* xHistoName, Int_t nPoints, Float_t* x, Float_t* x_err) const;
  void getYPoints( TFile * file, const char* yHistoName, Int_t nPoints, Float_t* y_resp, Float_t* y_resp_err,  Float_t* y_reso, Float_t* y_reso_err) const;

  bool NOQ_;
  std::string FIT_RMS_;
  float INTPERC_;
  bool FIXM_;

};


#endif
