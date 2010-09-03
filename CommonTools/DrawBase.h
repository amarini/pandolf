//************************************************************************
//
//      DrawBase is the Base Class from which all DrawX.C inherit
//
//      Draws Data-MC comparisons
//
//************************************************************************

#ifndef DrawBase_h
#define DrawBase_h

#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"


struct MCFile {
  TFile* file;
  std::string name;
  int fillColor;
  int fillStyle;
};


class DrawBase {

 public:

  DrawBase( const std::string& analysisType, const std::string& recoType="", const std::string& jetAlgo="" );
  ~DrawBase();

  void set_shapeNormalization();
  void set_lumiNormalization( float givenLumi=-1.);
  void set_sameEventNormalization();
  void set_sameInstanceNormalization();

  void drawHisto( std::string name, std::string etaRegion, std::string flags, std::string axisName="", int legendQuadrant=1, bool log_aussi=false);
  void drawProfile( std::string yVar, std::string xVar, int legendQuadrant=1);
  void drawStack(const std::string& varY, const std::string& varX, const std::string& RECO_GEN, bool isData) const { this->drawStack( varY, varX, "", RECO_GEN, isData); };
  void drawStack(const std::string& varY, const std::string& varX, const std::string& etaRegion, const std::string& RECO_GEN, bool isData) const;

  void set_analysisType( const std::string analysisType ) { analysisType_ = analysisType; };
  void set_dataFile( TFile* dataFile );
  void add_mcFile( TFile* mcFile, const std::string& name, int fillColor, int fillStyle=-1 );
  void set_outputdir( const std::string& outputdir ) { outputdir_ = outputdir; };
  void set_pt_thresh( Int_t pt_thresh ) { pt_thresh_ = pt_thresh; };
  void set_etamax( Float_t etamax ) { etamax_ = etamax; };
  void set_raw_corr( const std::string& raw_corr ) { raw_corr_ = raw_corr; };
  void set_pdf_aussi( bool pdf_aussi ) { pdf_aussi_ = pdf_aussi; };
  void set_scaleFactor( float scaleFactor ) { scaleFactor_ = scaleFactor;};
  void set_noStack( bool set ) { noStack_ = set; };

  std::string get_CMSText() const;
  std::string get_analysisType() const { return analysisType_; };
  TFile* get_dataFile() const { return dataFile_; };
  TFile* get_mcFile( int i ) const { return mcFiles_[i].file; };
  std::string get_outputdir() const { return outputdir_; };
  Int_t get_pt_thresh() const { return pt_thresh_; };
  Float_t get_etamax() const { return etamax_; };
  std::string get_raw_corr() const { return raw_corr_; };
  bool get_pdf_aussi() const { return pdf_aussi_; };

 private:

  std::string getEtaRangeText( const std::string& etaRegion ) const;
  std::string getSqrtText() const;
  std::string getAlgoName() const;
  void shrinkPad(double b=0.1, double l=0.1, double r=0.1, double t=0.1);
  std::string getAxisName(std::string name);
  TGraphErrors* getGraphRatio( TGraphErrors* gr_data, TGraphErrors* gr_MC);


  std::string analysisType_;
  std::string recoType_;
  std::string jetAlgo_;

  TFile* dataFile_;
  std::vector< MCFile > mcFiles_;
  
  Float_t scaleFactor_;
  Float_t lumi_;

  std::string outputdir_;
  Int_t pt_thresh_;
  Float_t etamax_;
  std::string raw_corr_;
  bool pdf_aussi_;

  bool noStack_;

};

#endif
