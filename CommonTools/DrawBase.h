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


struct InputFile {
  TFile* file;
  std::string datasetName;
  std::string legendName;
  int fillColor;
  int fillStyle;
};


struct HistoAndName {
  std::string histoName;
  std::string legendName;
};


struct LegendBox {

  float xMin;
  float xMax;
  float yMin;
  float yMax;

};


class DrawBase {

 public:

  DrawBase( const std::string& analysisType, const std::string& recoType="", const std::string& jetAlgo="", const std::string& flags="" );
  virtual ~DrawBase();

  void set_shapeNormalization();
  void set_lumiNormalization( float givenLumi=-1.);
  void set_sameEventNormalization();
  void set_sameInstanceNormalization();

  void drawHisto( const std::string& name, const std::string& etaRegion, const std::string& flags, const std::string& axisName="", int legendQuadrant=1, bool log_aussi=false);
  void drawProfile( const std::string& yVar, const std::string& xVar, int legendQuadrant=1);
  void drawStack(const std::string& varY, const std::string& varX, const std::string& RECO_GEN, bool isData) const { this->drawStack( varY, varX, "", RECO_GEN, isData); };
  void drawStack(const std::string& varY, const std::string& varX, const std::string& etaRegion, const std::string& RECO_GEN, bool isData) const;
  void compareDifferentHistos( const std::vector< HistoAndName > histosandnames, const std::string xAxisName, const std::string saveVarName, bool normalized=true, int legendQuadrant=1 );
  void compareDifferentHistos_singleFile( InputFile file, const std::vector< HistoAndName > histosandnames, const std::string xAxisName, const std::string saveVarName, bool normalized=true, int legendQuadrant=1 );
  void drawObjects( const std::vector< TObject* > objects, const std::string& name, 
                   const std::string& xAxisName, float xMin, float xMax, 
                   const std::string& yAxisName, float yMin, float yMax, 
                   bool logx=false, bool logy=false);

  void set_analysisType( const std::string analysisType ) { analysisType_ = analysisType; };
  void add_dataFile( TFile* dataFile, const std::string& datasetName );
  void add_mcFile( TFile* mcFile, const std::string& datasetName, const std::string& legendName, int fillColor, int fillStyle=-1 );
  void set_outputdir( const std::string& outputdir="" ); //if "" is passed, default outputdir is set (recommended)
  void set_flags( const std::string& flags ) { flags_ = flags; };
  void set_pt_thresh( Int_t pt_thresh ) { pt_thresh_ = pt_thresh; };
  void set_etamax( Float_t etamax ) { etamax_ = etamax; };
  void set_raw_corr( const std::string& raw_corr ) { raw_corr_ = raw_corr; };
  void set_pdf_aussi( bool pdf_aussi ) { pdf_aussi_ = pdf_aussi; };
  void set_scaleFactor( float scaleFactor ) { scaleFactor_ = scaleFactor;};
  void set_noStack( bool set ) { noStack_ = set; };
  void set_rebin( int rebin ) { rebin_ = rebin; };

  LegendBox get_legendBox( int legendQuadrant=1, const std::vector<std::string>* legendNames=0 ) const;
  TPaveText* get_labelCMS( int legendQuadrant=2 ) const;
  TPaveText* get_labelSqrt( int legendQuadrant=2 ) const;
  TPaveText* get_labelAlgo( int legendQuadrant=3 ) const;
  std::string get_CMSText() const;
  std::string get_analysisType() const { return analysisType_; };
  std::string get_recoType() const { return recoType_; };
  std::string get_flags() const { return flags_; };
  TFile* get_dataFile() const { return dataFile_.file; };
  TFile* get_mcFile( int i ) const { return mcFiles_[i].file; };
  std::string get_outputdir() const { return outputdir_; };
  Int_t get_pt_thresh() const { return pt_thresh_; };
  Float_t get_etamax() const { return etamax_; };
  std::string get_raw_corr() const { return raw_corr_; };
  bool get_pdf_aussi() const { return pdf_aussi_; };

  std::string get_etaRangeText( const std::string& etaRegion ) const;
  std::string get_sqrtText() const;
  std::string get_algoName() const;
  std::string get_axisName(std::string name);
  std::string get_outputSuffix() const;
  std::string get_fullSuffix() const;

 private:

  TGraphErrors* get_graphRatio( TGraphErrors* gr_data, TGraphErrors* gr_MC);


  std::string analysisType_;
  std::string recoType_;
  std::string jetAlgo_;

  std::string flags_;

  InputFile dataFile_;
  std::vector< InputFile > mcFiles_;
  
  Float_t scaleFactor_;
  Float_t lumi_;

  Int_t rebin_;

  std::string outputdir_;
  Int_t pt_thresh_;
  Float_t etamax_;
  std::string raw_corr_;
  bool pdf_aussi_;

  bool noStack_;

};

#endif
