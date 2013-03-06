//************************************************************************
//
//      DrawBase is the Base Class from which all DrawX.C inherit
//
//      Draws Data-MC comparisons
//
//************************************************************************

#ifndef DrawBase_h
#define DrawBase_h

#include <vector>

#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"


struct InputFile {
  TFile* file;
  std::string datasetName;
  float weight;
  std::string legendName;
  int fillColor;
  int fillStyle;
  int markerStyle;
  int lineColor;
  int lineWidth;
};


class HistoAndName {
 public:
  HistoAndName(int markerstyle=-1) { markerStyle=markerstyle; };
  std::string histoName;
  std::string legendName;
  int markerStyle;
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

  //void drawHisto( const std::string& name, const std::string& etaRegion, const std::string& flags, const std::string& axisName="", const std::string& units="", int legendQuadrant=1, bool log_aussi=false);
  void drawHisto_vs_pt( int nBinsPt, float* ptBins, const std::string& name, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log_aussi=false, int legendQuadrant=1, std::string flags="", const std::string& labelText="" );
  void drawHisto_vs_pt( std::vector<float> ptBins, const std::string& name, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log_aussi=false, int legendQuadrant=1, std::string flags="", const std::string& labelText="" );
  TCanvas* drawHisto( const std::string& name, const std::string& axisName="", const std::string& units="", const std::string& instanceName="Events", bool log_aussi=false, int legendQuadrant=1, const std::string& flags="", const std::string& labelText="", bool add_jetAlgoText=false );
  TCanvas* drawHisto_fromHistos( std::vector<TH1D*> dataHistos, std::vector<TH1D*> mcHistos, std::vector<TH1D*> mcHistos_superimp, const std::string& name, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log_aussi=false, int legendQuadrant=1, const std::string& flags="", const std::string& labelText="", bool add_jetAlgoText=false );
  TCanvas* drawHisto_fromTree( const std::string& treeName, const std::string& varName, const std::string& selection, int nBins, float xMin, float xMax, const std::string& name, const std::string& axisName, const std::string& units="", const std::string& instanceName="Entries", bool log_aussi=false, int legendQuadrant=1, const std::string& flags="", const std::string& labelText="", bool add_jetAlgoText=false );
  void draw2DHisto_fromTree( const std::string& treeName, const std::string& varName1, const std::string& varName2, const std::string& selection, int nBins1, float xMin1, float xMax1, int nBins2, float xMin2, float xMax2, const std::string& name, const std::string& axisName1, const std::string& axisName2, const std::string& units1="", const std::string& units2="", int legendQuadrant=1, const std::string& flags="", const std::string& labelText="" );
  void drawProfile( const std::string& yVar, const std::string& xVar, int legendQuadrant=1);
  void drawStack(const std::string& varY, const std::string& varX, const std::string& RECO_GEN, bool isData) const { this->drawStack( varY, varX, "", RECO_GEN, isData); };
  void drawStack(const std::string& varY, const std::string& varX, const std::string& etaRegion, const std::string& RECO_GEN, bool isData) const;
  void compareDifferentHistos( const std::vector< HistoAndName > histosandnames, const std::string saveVarName, const std::string xAxisName="", const std::string& units="", const std::string& instanceName="Entries", bool stacked=false, int legendQuadrant=1 );
  void compareDifferentHistos_singleFile( InputFile file, const std::vector< HistoAndName > histosandnames, const std::string saveVarName, const std::string xAxisName="", const std::string& units="", const std::string& instanceName="Entries", bool stacked=false, int legendQuadrant=1 );
  void drawObjects( const std::vector< TObject* > objects, const std::string& name, 
                   const std::string& xAxisName, float xMin, float xMax, 
                   const std::string& yAxisName, float yMin, float yMax, 
                   bool logx=false, bool logy=false);



  void clearFiles(); 
  void reset(); 

  void add_dataFile( TFile* dataFile, const std::string& datasetName, const std::string& legendName="Data", int markerColor=-1, int markerStyle=-1, int fillStyle=-1 );
  void add_dataFile( TFile* dataFile, float weight, const std::string& datasetName, const std::string& legendName="Data", int markerColor=-1, int markerStyle=-1, int fillStyle=-1 );
  void add_mcFile( TFile* mcFile, const std::string& datasetName, const std::string& legendName, int fillColor=-1, int fillStyle=-1, int markerStyle=-1, int lineColor=-1, int lineWidth=-1 );
  void add_mcFile_superimp( TFile* mcFile, const std::string& datasetName, const std::string& legendName, float multFactor=1., int lineColor=-1, int lineWidth=-1 );
  // in the following function weight must be cross_section(in pb) / Nevents:
  void add_mcFile( TFile* mcFile, float weight, const std::string& datasetName, const std::string& legendName, int fillColor=-1, int fillStyle=-1, int markerStyle=-1, int lineColor=-1, int lineWidth=-1 );

  void set_analysisType( const std::string analysisType ) { analysisType_ = analysisType; };
  void set_lumi( float lumi ) { lumi_ = lumi; };
  void set_outputdir( const std::string& outputdir="" ); //if "" is passed, default outputdir is set
  void set_flags( const std::string& flags ) { flags_ = flags; };
  void set_pt_thresh( Int_t pt_thresh ) { pt_thresh_ = pt_thresh; };
  void set_etamax( Float_t etamax ) { etamax_ = etamax; };
  void set_raw_corr( const std::string& raw_corr ) { raw_corr_ = raw_corr; };
  void set_pdf_aussi( bool pdf_aussi ) { pdf_aussi_ = pdf_aussi; };
  void set_root_aussi( bool root_aussi ) { root_aussi_ = root_aussi; };
  void set_logx( bool logx=true ) { logx_ = logx; };
  void set_scaleFactor( float scaleFactor ) { scaleFactor_ = scaleFactor;};
  void set_xAxisRange( float xAxisMin=9999., float xAxisMax=9999. ) { 
    xAxisMin_ = xAxisMin;
    xAxisMax_ = xAxisMax;};
  void set_xAxisMin( float xAxisMin=9999. ) { xAxisMin_ = xAxisMin;};
  void set_xAxisMax( float xAxisMax=9999. ) { xAxisMax_ = xAxisMax;};
  void set_yAxisMax( float yAxisMax=9999. ) { yAxisMax_ = yAxisMax;};
  void set_yAxisMaxScale( float yAxisMaxScale=1.4 ) { yAxisMaxScale_ = yAxisMaxScale;};
  void set_yAxisMaxScaleLog( float yAxisMaxScale=5. ) { yAxisMaxScaleLog_ = yAxisMaxScale;};
  void set_noStack( bool set=true ) { noStack_ = set; };
  void set_isCMSArticle( bool set=true );
  void set_drawZeros( bool set=true ) { drawZeros_=set; };
  void set_lumiOnRightSide( bool set=true );
  void set_rebin( int rebin=1 ) { rebin_ = rebin; };
  void set_mcMarkers( bool set=true );
  void set_markerSize( float markerSize ) { markerSize_ = markerSize; };
  void set_getBinLabels( bool getBinL=true ) { getBinLabels_ = getBinL; };
  void set_legendTitle( const std::string& title ) { legendTitle_ = title; };
  void set_mcWeight( int iFile, float weight );
  void set_mcWeight( const std::string& datasetName, float weight );
  void set_is7TeV( bool is7TeV ) { is7TeV_ = is7TeV; };
  void set_displayEmptyDatasets( bool displayThem ) { displayEmptyDatasets_ = displayThem; };
  void set_legendBox_xMax( float legendBox_xMax ) { legendBox_xMax_ = legendBox_xMax; };
  void set_legendBox_xMin( float legendBox_xMin ) { legendBox_xMin_ = legendBox_xMin; };
  void set_legendBox_yMax( float legendBox_yMax ) { legendBox_yMax_ = legendBox_yMax; };
  void set_legendBox_yMin( float legendBox_yMin ) { legendBox_yMin_ = legendBox_yMin; };
  void set_noMarkerBarsX( bool noMarkerBarsX=true ) { noMarkerBarsX_ = noMarkerBarsX; };
  void set_graphLineWidth( float graphLineWidth=1. ) { graphLineWidth_ = graphLineWidth; };


  void add_label( const std::string& text, float xmin=0.23, float ymin=0.87, float xmax=0.36, float ymax=0.9, float textSize=0.035 );
  void delete_label();

  void resetStyle() const { style_->cd(); };

  float get_lumi() const { return lumi_; };
  LegendBox get_legendBox( int legendQuadrant=1, const std::vector<std::string>* legendNames=0 ) const;
  TPaveText* get_labelTop() { return this->get_labelSqrt(); };
  TPaveText* get_labelCMS( int legendQuadrant=0 ) const;
  TPaveText* get_labelCMStop( bool wide=false ) const;
  TPaveText* get_labelSqrt( int legendQuadrant=0 ) const;
  TPaveText* get_labelAlgo( int legendQuadrant=3 ) const;
  std::string get_CMSText() const;
  std::string get_analysisType() const { return analysisType_; };
  std::string get_recoType() const { return recoType_; };
  std::string get_flags() const { return flags_; };
  std::string get_legendTitle() const { return legendTitle_; };
  InputFile get_dataFile( int i ) const { return dataFiles_[i]; };
  InputFile get_mcFile( int i ) const { return mcFiles_[i]; };
  std::vector< InputFile >  get_dataFiles() const { return dataFiles_; };
  std::vector< InputFile >  get_mcFiles() const { return mcFiles_; };
  std::string get_outputdir() const { return outputdir_; };
  Int_t get_pt_thresh() const { return pt_thresh_; };
  Float_t get_etamax() const { return etamax_; };
  std::string get_raw_corr() const { return raw_corr_; };
  bool get_pdf_aussi() const { return pdf_aussi_; };
  int get_rebin() const { return rebin_; };

  std::string get_etaRangeText( const std::string& etaRegion ) const;
  std::string get_sqrtText() const;
  std::string get_lumiText() const;
  std::string get_algoName() const;
  std::string get_algoType() const;
  std::string get_axisName(std::string name);
  std::string get_outputSuffix() const;
  std::string get_fullSuffix() const;

  TH2D* get_lastHistos_axes() const { return lastHistos_axes_; };
  std::vector< TH1D* > get_lastHistos_data() const { return lastHistos_data_; };
  TGraphAsymmErrors* get_lastHistos_dataGraph() const { return lastHistos_dataGraph_; };
  std::vector< TH1D* > get_lastHistos_mc() const { return lastHistos_mc_; };
  std::vector< TH1D* > get_lastHistos_mc_superimp() const { return lastHistos_mc_superimp_; };
  THStack* get_lastHistos_mcStack() const { return lastHistos_mcStack_; };
  //TList* get_lastHistos_primitives() const { return lastHistos_primitives_; };

  float get_xAxisMin() const { return xAxisMin_; };
  float get_xAxisMax() const { return xAxisMax_; };
  float get_yAxisMaxScale() const { return yAxisMaxScale_; };



 private:

  TGraphErrors* get_graphRatio( TGraphErrors* gr_data, TGraphErrors* gr_MC);

  
  TStyle* style_;

  std::string analysisType_;
  std::string recoType_;
  std::string jetAlgo_;

  std::string flags_;

  std::vector< InputFile > dataFiles_;
  std::vector< InputFile > mcFiles_;
  std::vector< InputFile > mcFiles_superimp_;
  
  TH1D* lastHistos_mcHistoSum_;
  std::vector< TH1D* > lastHistos_data_;
  TGraphAsymmErrors* lastHistos_dataGraph_;
  std::vector< TH1D* > lastHistos_mc_;
  std::vector< TH1D* > lastHistos_mc_superimp_;
  TH2D* lastHistos_axes_;
  THStack* lastHistos_mcStack_;
  //TList* lastHistos_primitives_;
  
  Float_t scaleFactor_;

  Float_t xAxisMin_;
  Float_t xAxisMax_;
  Float_t yAxisMax_;
  Float_t yAxisMaxScale_;
  Float_t yAxisMaxScaleLog_;

  Float_t markerSize_;
  Float_t graphLineWidth_;

  Float_t lumi_;
  bool is7TeV_;

  Int_t rebin_;

  std::string outputdir_;
  Int_t pt_thresh_;
  Float_t etamax_;
  std::string raw_corr_;
  bool pdf_aussi_;
  bool root_aussi_;
  bool logx_;
  bool getBinLabels_;
  bool noMarkerBarsX_;

  std::string legendTitle_;
  float legendTextSize_;

  bool poissonAsymmErrors_;
  bool drawZeros_;

  TPaveText* additionalLabel_;

  bool noStack_;

  bool isCMSArticle_;
  bool lumiOnRightSide_;

  bool displayEmptyDatasets_;

  float legendBox_xMin_;
  float legendBox_xMax_;
  float legendBox_yMin_;
  float legendBox_yMax_;


};

#endif
