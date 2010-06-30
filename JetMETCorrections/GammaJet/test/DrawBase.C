#include "DrawBase.h"
#include "fitTools.h"


DrawBase::DrawBase( const std::string& analysisType, const std::string& recoType, const std::string& jetAlgo ) {

  TStyle *simpleStyle = new TStyle("simpleStyle","");
  simpleStyle->SetCanvasColor(0);
  simpleStyle->SetPadColor(0);
  simpleStyle->SetFrameFillColor(0);
  simpleStyle->SetStatColor(0);
  simpleStyle->SetOptStat(0);
  simpleStyle->SetTitleFillColor(0);
  simpleStyle->SetCanvasBorderMode(0);
  simpleStyle->SetPadBorderMode(0);
  simpleStyle->SetFrameBorderMode(0);
  simpleStyle->cd();

  analysisType_ = analysisType;
  recoType_ = recoType;
  jetAlgo_ = jetAlgo;

  dataFile_ = 0;
  mcFile_ = 0;
  mcFile2_ = 0;

}


/*DrawBase::DrawBase(  const std::string analysisType, TFile* dataFile, TFile* mcFile, std::string outputdir, Int_t pt_thresh, Float_t etamax, std::string raw_corr, bool pdf_aussi ) {

  TStyle *simpleStyle = new TStyle("simpleStyle","");
  simpleStyle->SetCanvasColor(0);
  simpleStyle->SetPadColor(0);
  simpleStyle->SetFrameFillColor(0);
  simpleStyle->SetStatColor(0);
  simpleStyle->SetOptStat(0);
  simpleStyle->SetTitleFillColor(0);
  simpleStyle->SetCanvasBorderMode(0);
  simpleStyle->SetPadBorderMode(0);
  simpleStyle->SetFrameBorderMode(0);
  simpleStyle->cd();

  dataFile_ = dataFile;
  mcFile_ = mcFile;

  outputdir_ = outputdir;
  pt_thresh_ = pt_thresh;
  etamax_ = etamax;
  raw_corr_ = raw_corr;
  pdf_aussi_ = pdf_aussi;

}
*/

DrawBase::~DrawBase() {

  if( dataFile_!=0 ) {
    delete dataFile_;
    dataFile_=0;
  }

  if( mcFile_!=0 ) {
    delete mcFile_;
    mcFile_=0;
  }

  if( mcFile2_!=0 ) {
    delete mcFile2_;
    mcFile2_=0;
  }

}


void DrawBase::set_crossSectionNormalization() {

  if( dataFile_==0 || mcFile_==0 ) {
    std::cout << "Data/MC files not properly initialized. Cannot compute scale factor. Exiting." << std::endl;
    exit(132);
  }

  if( lumi_==0. ) {
    std::cout << "lumi_ is equal to zero. Cannot compute cross section normalization. Exiting." << std::endl;
    exit(97);
  }

  scaleFactor_ = lumi_/1000000.; //eventweights are set so that histos are number of expected evets @ 1 pb-1. lumi is in mub-1


}

void DrawBase::set_shapeNormalization() {

  scaleFactor_ = -1.;

}

void DrawBase::set_sameEventNormalization() {

  if( dataFile_==0 || mcFile_==0 ) {
    std::cout << "Data/MC files not properly initialized. Cannot compute scale factor. Exiting." << std::endl;
    exit(133);
  }

  TH1F* h1_nJets_data = (TH1F*)dataFile_->Get("nJets");

  TH1F* h1_nJets_mc = (TH1F*)mcFile_->Get("nJets");
  if( mcFile2_!=0 ) {
    TH1F* h1_nJets_mc2 = (TH1F*)mcFile2_->Get("nJets"); 
    h1_nJets_mc->Add( h1_nJets_mc2 );
  }

  scaleFactor_ =  h1_nJets_data->Integral()/h1_nJets_mc->Integral();

}

void DrawBase::set_sameInstanceNormalization() {

  if( dataFile_==0 || mcFile_==0 ) {
    std::cout << "Data/MC files not properly initialized. Cannot compute scale factor. Exiting." << std::endl;
    exit(133);
  }

  TH1F* h1_phi_data;
  TH1F* h1_phi_mc;
  std::string hname;
  if( analysisType_=="MinBias" ) {
    hname = "phiJet";
  } else if( analysisType_=="PhotonJet" ) {
    hname = "phiPhot";
  } else {
    std::cout << "Analysis type '" << analysisType_ << "' not implemented yet. Exiting." << std::endl;
    exit(765);
  }

  h1_phi_data = (TH1F*)dataFile_->Get(hname.c_str());
  h1_phi_mc = (TH1F*)mcFile_->Get(hname.c_str());
  if( mcFile2_!=0 ) {
    TH1F* h1_phi_mc2 = (TH1F*)mcFile2_->Get(hname.c_str());
    h1_phi_mc->Add( h1_phi_mc2 );
  }
  scaleFactor_ =  h1_phi_data->Integral()/h1_phi_mc->Integral();
  std::cout << "ScaleFactor: " << scaleFactor_ << std::endl;

}


void DrawBase::drawHisto( std::string name, std::string etaRegion, std::string flags, int legendQuadrant, bool log_aussi) {

  std::string histoName = name;
  if( etaRegion!="" ) histoName = name + "_" + etaRegion;

  TH1F* dataHisto = (TH1F*)dataFile_->Get(histoName.c_str());
  TH1F* mcHisto = (TH1F*)mcFile_->Get(histoName.c_str());

  if( dataHisto==0 || mcHisto==0 ) {
    std::cout << "Didn't find histo '" << histoName << "'. Continuing." << std::endl;
    return;
  }

  mcHisto->SetFillColor(38);
//mcHisto->SetFillColor(kRed-7);
//mcHisto->SetLineColor(kRed);
//mcHisto->SetLineWidth(2);

  if( scaleFactor_ > 0. ) 
    mcHisto->Scale(scaleFactor_);
  else {
    Float_t dataIntegral = dataHisto->Integral(0, dataHisto->GetNbinsX()+1);
    Float_t mcIntegral = mcHisto->Integral(0, mcHisto->GetNbinsX()+1);
    mcHisto->Scale( dataIntegral/mcIntegral );
  }

  dataHisto->SetMarkerStyle(20);

  Float_t yAxisMaxScale = (name=="phiJet" || name=="etaJet") ? 1.8 : 1.6;
  Float_t xMin = dataHisto->GetXaxis()->GetXmin();
  Float_t xMax = dataHisto->GetXaxis()->GetXmax();
  Float_t yMax_data = dataHisto->GetMaximum();
  Float_t yMax_mc = mcHisto->GetMaximum();
  Float_t yMax = (yMax_data>yMax_mc) ? yAxisMaxScale*yMax_data : yAxisMaxScale*yMax_mc;
  Float_t yMin = 0.;


  std::string xAxis = getAxisName(name);

  std::string instanceName = (analysisType_=="PhotonJet") ? "Events" : "Jets";
  std::string yAxis = instanceName;
  if( name=="ptJet" || name=="ptCorrJet" ) {
    char yAxis_char[60];
    sprintf(yAxis_char, "%s/(%d GeV/c)", instanceName.c_str(), (Int_t)dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  

  if( name=="EchJet" || name=="EgammaJet" || name=="EnhJet" ) {
    char yAxis_char[60];
    if( dataHisto->GetBinWidth(1) < 1. ) sprintf(yAxis_char, "%s/(%.1f GeV)", instanceName.c_str(), dataHisto->GetBinWidth(1));
    else                                 sprintf(yAxis_char, "%s/(%.0f GeV)", instanceName.c_str(), dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  


  if( name=="diJetMass" ) {
    char yAxis_char[60];
    sprintf(yAxis_char, "Events/(%.2f GeV/c^{2})", dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  
  if( name=="massJet" ) {
    char yAxis_char[60];
    sprintf(yAxis_char, "Jets/(%.2f GeV/c^{2})", dataHisto->GetBinWidth(1));
    //sprintf(yAxis_char, "Jets/(%d GeV/c^{2})", (Int_t)dataHisto->GetBinWidth(1));
    std::string yAxis_tmp(yAxis_char);
    yAxis=yAxis_tmp;
  }  
  if( name=="deltaPhiJet" || name=="asymmJet" ) yAxis = "Events";

  std::string etaRange = getEtaRangeText(etaRegion);

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax);
  h2_axes->SetXTitle(xAxis.c_str());
  h2_axes->SetYTitle(yAxis.c_str());
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.5);

  if( legendQuadrant<0 || legendQuadrant>5 ) {
    std::cout << "Invalid legendQuadrant! Exiting!" << std::endl;
    exit(1);
  }

  Float_t legend_xMin;
  Float_t legend_yMin;
  Float_t legend_xMax;
  Float_t legend_yMax;

  if( legendQuadrant==1 ) {
    legend_xMin = 0.6;
    legend_yMin = 0.73;
    legend_xMax = 0.88;
    legend_yMax = 0.88;
  } else if( legendQuadrant==0 ) {
    legend_xMin = 0.5;
    legend_yMin = 0.73;
    legend_xMax = 0.73;
    legend_yMax = 0.88;
  } else if( legendQuadrant==2 ) {
    legend_xMin = 0.18;
    legend_yMin = 0.73;
    legend_xMax = 0.41;
    legend_yMax = 0.88;
  } else if( legendQuadrant==5 ) {
    legend_xMin = 0.4;
    legend_yMin = 0.15;
    legend_xMax = 0.6;
    legend_yMax = 0.25;
  }


  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.035);
  legend->AddEntry(dataHisto, "Data", "P");
  legend->AddEntry(mcHisto, "Simulation", "F");

  TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.038);
  label_cms->SetTextFont(62);
  label_cms->AddText("CMS Preliminary 2010");

  TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  std::string label_sqrt_text = this->getSqrtText();
  label_sqrt->AddText(label_sqrt_text.c_str());
  

  TPaveText* label_cuts;

  if( analysisType_ == "MinBias" ) {

    Float_t label_cuts_xMin = 0.63;
    Float_t label_cuts_yMin = 0.5;
    Float_t label_cuts_xMax = 0.84;
    Float_t label_cuts_yMax = 0.65;

    if( name=="asymmJet" || name=="deltaPhiJet" || name=="RchJet" || name=="etaJet" || name=="phiJet" ) {
     label_cuts_xMin = 0.25;
     label_cuts_yMin = 0.60;
     label_cuts_xMax = 0.36;
     label_cuts_yMax = 0.73;
    }
    if( name=="RnhJet" ) {
     label_cuts_xMin = 0.4;
     label_cuts_yMin = 0.6;
     label_cuts_xMax = 0.6;
     label_cuts_yMax = 0.7;
    }

    label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax,  "brNDC");
    label_cuts->SetFillColor(kWhite);
    label_cuts->SetTextSize(0.035);
    label_cuts->SetTextFont(42);
    std::string jetAlgoName = this->getAlgoName();
    label_cuts->AddText(jetAlgoName.c_str());
    if( name != "etaJet" )
      label_cuts->AddText(etaRange.c_str());
    if( name != "ptJet" && name != "ptCorrJet" ) {
      char labelText[70];
      sprintf( labelText, "p_{T}^{%s} > %d GeV/c", raw_corr_.c_str(), pt_thresh_);
      label_cuts->AddText(labelText);
    }

  } else if( analysisType_=="PhotonJet" ) {

    Float_t label_cuts_xMin = 0.63;
    Float_t label_cuts_yMin = 0.5;
    Float_t label_cuts_xMax = 0.84;
    Float_t label_cuts_yMax = 0.65;

    if( name=="deltaPhi" ) {
     label_cuts_xMin = 0.25;
     label_cuts_yMin = 0.60;
     label_cuts_xMax = 0.36;
     label_cuts_yMax = 0.73;
    }
    label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax,  "brNDC");
    label_cuts->SetFillColor(kWhite);
    label_cuts->SetTextSize(0.035);
    label_cuts->SetTextFont(42);
    label_cuts->AddText("|#eta_{#gamma}| < 1.3");
    label_cuts->AddText("Anti-k_{T} 0.5 PFJets");

  } else {

    std::cout << "AnalysisType: '" << analysisType_ << "' not yet implemented. Exiting." << std::endl;
    exit(991);

  }

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();


  h2_axes->Draw("");
  mcHisto->Draw("histo same");
  dataHisto->Draw("E same");
  gPad->RedrawAxis();
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_cuts->Draw("same");



  std::string canvasName = outputdir_ + "/" + name;
  if( etaRegion!="" )
    canvasName = canvasName + "_" + etaRegion;
  if( flags!="" )
    canvasName = canvasName + "_" + flags;
  std::string canvasName_eps = canvasName + ".eps";
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SaveAs(canvasName_eps.c_str());
  std::string canvasName_png = canvasName + ".png";
  c1->SaveAs(canvasName_png.c_str());
  std::string canvasName_pdf = canvasName + ".pdf";
  if( pdf_aussi_ )
    c1->SaveAs(canvasName_pdf.c_str());

  if( log_aussi ) {
    TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.5, 5.*yMax);
    h2_axes_log->SetXTitle(xAxis.c_str());
    h2_axes_log->SetYTitle(yAxis.c_str());
    h2_axes_log->GetXaxis()->SetTitleOffset(1.1);
    h2_axes_log->GetYaxis()->SetTitleOffset(1.5);
    c1->SetLogy();
    h2_axes_log->Draw("");
    mcHisto->Draw("histo same");
    dataHisto->Draw("E same");
    gPad->RedrawAxis();
    legend->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cuts->Draw("same");
    std::string canvasName_log = canvasName + "_log";
    canvasName_eps = canvasName_log + ".eps";
    c1->SaveAs(canvasName_eps.c_str());
    canvasName_png = canvasName_log + ".png";
    c1->SaveAs(canvasName_png.c_str());
    canvasName_pdf = canvasName_log + ".pdf";
    if( pdf_aussi_ )
      c1->SaveAs(canvasName_pdf.c_str());
    delete h2_axes_log;
  }
    
  delete c1;
  delete legend;
  delete h2_axes;

} //drawHisto


void DrawBase::drawHisto_onlyData( std::string name, std::string etaRegion, std::string flags, int legendQuadrant, bool log_aussi) {

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  // if response will have to do a plot for each pt bin:
  int number_of_plots = (name!="response") ? 1 : (ptPhot_binning.size()-1);

  std::string ptPhotMean_name = "ptPhotMean";
  if( flags!= "" ) ptPhotMean_name = ptPhotMean_name + "_" + flags;

  TProfile* hp_ptPhot_data = (TProfile*)dataFile_->Get(ptPhotMean_name.c_str());

  if( name=="response" && (hp_ptPhot_data==0) ) {
    std::cout << "Didn't find '" << ptPhotMean_name << "' profile. Skipping" << std::endl;
    return;
  }

  TGraphErrors* gr_response_vs_pt = new TGraphErrors(0);
  gr_response_vs_pt->SetName("response_vs_pt");

  TGraphErrors* gr_resolution_vs_pt = new TGraphErrors(0);
  gr_resolution_vs_pt->SetName("resolution_vs_pt");



  for( unsigned iplot=0; iplot<number_of_plots; ++iplot ) {

    std::string histoName = name;
    if( etaRegion!="" ) histoName = histoName + "_" + etaRegion;
    if( flags!="" ) histoName = histoName + "_" + flags;

    std::string ptRange_str; 
    if( name=="response" ) {
      char ptRange[60];
      sprintf( ptRange, "ptPhot_%.0f_%.0f", ptPhot_binning[iplot], ptPhot_binning[iplot+1]);
      std::string ptRange_str_tmp(ptRange);
      ptRange_str = ptRange_str_tmp;
      histoName = histoName + "_" + ptRange_str;
    }


    TH1F* dataHisto = (TH1F*)dataFile_->Get(histoName.c_str());

    if( dataHisto==0 ) {
      std::cout << "Didn't find histo '" << histoName << "'. Continuing." << std::endl;
      return;
    }

    dataHisto->SetMarkerStyle(20);



    Float_t yAxisMaxScale = (name=="phiJet" || name=="etaJet" || name=="ptSecondJetRel" || name=="phiPhot" || name=="etaPhot" ) ? 1.8 : 1.6;
    if( name=="phiPhot" || name=="etaPhot" ) yAxisMaxScale=2.;
    if(name=="clusterMajPhotReco" || name=="clusterMinPhotReco") yAxisMaxScale = 2.;
    Float_t xMin = dataHisto->GetXaxis()->GetXmin();
    Float_t xMax = dataHisto->GetXaxis()->GetXmax();
    Float_t yMax = dataHisto->GetMaximum();
    Float_t yMin = 0.;


    std::string xAxis = getAxisName(name);

    std::string instanceName = (analysisType_=="PhotonJet") ? "Events" : "Jets";
    std::string yAxis = instanceName;
    if( name=="ptJet" || name=="ptCorrJet"||name=="ptPhot" ) {
      char yAxis_char[50];
      sprintf(yAxis_char, "%s/(%d GeV/c)", instanceName.c_str(), (Int_t)dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  

    if( name=="EchJet" || name=="EgammaJet" || name=="EnhJet" ) {
      char yAxis_char[50];
      if( dataHisto->GetBinWidth(1) < 1. ) sprintf(yAxis_char, "%s/(%.1f GeV)", instanceName.c_str(), dataHisto->GetBinWidth(1));
      else                                 sprintf(yAxis_char, "%s/(%.0f GeV)", instanceName.c_str(), dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  


    if( name=="diJetMass" ) {
      char yAxis_char[50];
      sprintf(yAxis_char, "Events/(%.2f GeV/c^{2})", dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  
    if( name=="massJet" ) {
      char yAxis_char[50];
      sprintf(yAxis_char, "Jets/(%.2f GeV/c ^{2})", dataHisto->GetBinWidth(1));
      //sprintf(yAxis_char, "Jets/(%d GeV/c^{2})", (Int_t)dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  
    if( name=="deltaPhiJet" || name=="asymmJet" ) yAxis = "Events";

    std::string etaRange = getEtaRangeText( etaRegion );

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax*1.6);
    h2_axes->SetXTitle(xAxis.c_str());
    h2_axes->SetYTitle(yAxis.c_str());
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.5);

    if( legendQuadrant<0 || legendQuadrant>5 ) {
      std::cout << "Invalid legendQuadrant! Exiting!" << std::endl;
      exit(1);
    }

    Float_t legend_xMin;
    Float_t legend_yMin;
    Float_t legend_xMax;
    Float_t legend_yMax;

    if( legendQuadrant==1 ) {
      legend_xMin = 0.6;
      legend_yMin = 0.73;
      legend_xMax = 0.88;
      legend_yMax = 0.88;
    } else if( legendQuadrant==0 ) {
      legend_xMin = 0.5;
      legend_yMin = 0.73;
      legend_xMax = 0.73;
      legend_yMax = 0.88;
    } else if( legendQuadrant==2 ) {
      legend_xMin = 0.18;
      legend_yMin = 0.73;
      legend_xMax = 0.41;
      legend_yMax = 0.88;
    } else if( legendQuadrant==5 ) {
      legend_xMin = 0.4;
      legend_yMin = 0.15;
      legend_xMax = 0.6;
      legend_yMax = 0.25;
    }


    TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax);
    legend->SetFillColor(kWhite);
    legend->SetTextSize(0.035);
    legend->AddEntry(dataHisto, "Data", "P");

    TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
    label_cms->SetFillColor(kWhite);
    label_cms->SetTextSize(0.038);
    label_cms->SetTextFont(62);
    label_cms->AddText("CMS Preliminary 2010");

    TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
    label_sqrt->SetFillColor(kWhite);
    label_sqrt->SetTextSize(0.038);
    label_sqrt->SetTextFont(42);
    std::string label_sqrt_text = this->getSqrtText();
    label_sqrt->AddText(label_sqrt_text.c_str());

    TPaveText* label_cuts;
    if( analysisType_ == "MinBias" ) {

      Float_t label_cuts_xMin = 0.63;
      Float_t label_cuts_yMin = 0.5;
      Float_t label_cuts_xMax = 0.84;
      Float_t label_cuts_yMax = 0.65;

      if( name=="asymmJet" || name=="deltaPhiJet" || name=="RchJet" || name=="etaJet" || name=="phiJet" ) {
       label_cuts_xMin = 0.25;
       label_cuts_yMin = 0.60;
       label_cuts_xMax = 0.36;
       label_cuts_yMax = 0.73;
      }
      if( name=="RnhJet" ) {
       label_cuts_xMin = 0.4;
       label_cuts_yMin = 0.6;
       label_cuts_xMax = 0.6;
       label_cuts_yMax = 0.7;
      }

      label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax,  "brNDC");
      label_cuts->SetFillColor(kWhite);
      label_cuts->SetTextSize(0.035);
      label_cuts->SetTextFont(42);
      std::string jetAlgoName = getAlgoName();
      label_cuts->AddText(jetAlgoName.c_str());
      if( name != "etaJet" )
        label_cuts->AddText(etaRange.c_str());
      if( name != "ptJet" && name != "ptCorrJet" ) {
        char labelText[70];
        sprintf( labelText, "p_{T}^{%s} > %d GeV/c", raw_corr_.c_str(), pt_thresh_);
        label_cuts->AddText(labelText);
      }

    } else if( analysisType_=="PhotonJet" ) {

      Float_t label_cuts_xMin = 0.63;
      Float_t label_cuts_yMin = 0.53;
      Float_t label_cuts_xMax = 0.84;
      Float_t label_cuts_yMax = 0.68;

      if( name=="deltaPhi" ) {
       label_cuts_xMin = 0.25;
       label_cuts_yMin = 0.60;
       label_cuts_xMax = 0.36;
       label_cuts_yMax = 0.73;
      }
      label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax,  "brNDC");
      label_cuts->SetFillColor(kWhite);
      label_cuts->SetTextSize(0.035);
      label_cuts->SetTextFont(42);
      
      if( name != "ptPhot" ) {
        char ptCutName[60];
        if( name != "response" )
          sprintf( ptCutName, "p_{T}^{#gamma} > %.0f GeV/c", ptPhot_binning[0]);
        else
          sprintf( ptCutName, "%.0f < p_{T}^{#gamma} < %.0f GeV/c", ptPhot_binning[iplot], ptPhot_binning[iplot+1] );

        label_cuts->AddText(ptCutName);
      }

      if( name != "etaPhot" )
        label_cuts->AddText("|#eta_{#gamma}| < 1.3");

      if( name=="response"||name=="deltaPhi"||name=="ptSecondJetRel" ) {
        std::string jetAlgoName = getAlgoName();
        label_cuts->AddText(jetAlgoName.c_str());
      }

    } else {

      std::cout << "AnalysisType: '" << analysisType_ << "' not yet implemented. Exiting." << std::endl;
      exit(991);

    }


    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();

    h2_axes->Draw("");
    dataHisto->Draw("E same");
    gPad->RedrawAxis();
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cuts->Draw("same");

    if( name=="response" ) {

      //label_mean->Draw("same");

      Float_t dataResponse = dataHisto->GetMean();
      Float_t dataResponseErr = dataHisto->GetMeanError();
      Float_t dataRMS = dataHisto->GetRMS();
      Float_t dataRMSErr = dataHisto->GetMeanError();
      Float_t dataResolution = dataRMS / dataResponse;
      Float_t dataResolutionErr = sqrt( dataRMSErr*dataRMSErr/(dataResponse*dataResponse) + dataResolution*dataResolution*dataResponseErr*dataResponseErr/(dataResponse*dataResponse*dataResponse*dataResponse) );

      gr_response_vs_pt->SetPoint( iplot, hp_ptPhot_data->GetBinContent( iplot+1 ), dataResponse );
      gr_response_vs_pt->SetPointError( iplot, hp_ptPhot_data->GetBinError( iplot+1 ), dataResponseErr );

      gr_resolution_vs_pt->SetPoint( iplot, hp_ptPhot_data->GetBinContent( iplot+1 ), dataResolution );
      gr_resolution_vs_pt->SetPointError( iplot, hp_ptPhot_data->GetBinError( iplot+1 ), dataResolutionErr );

    } //if reponse


    legend->Draw("same");
      

    std::string canvasName = outputdir_ + "/" + name;
    if( etaRegion!="" )
      canvasName = canvasName + "_" + etaRegion;
    if( flags!="" )
      canvasName = canvasName + "_" + flags;

    if( name=="response" )
      canvasName = canvasName + "_" + ptRange_str;

    canvasName = canvasName + "_ONLYDATA";

    std::string canvasName_eps = canvasName + ".eps";
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SaveAs(canvasName_eps.c_str());
    std::string canvasName_png = canvasName + ".png";
    c1->SaveAs(canvasName_png.c_str());
    std::string canvasName_pdf = canvasName + ".pdf";
    if( pdf_aussi_ )
      c1->SaveAs(canvasName_pdf.c_str());

    if( log_aussi ) {
      TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.5, 5.*yMax);
      h2_axes_log->SetXTitle(xAxis.c_str());
      h2_axes_log->SetYTitle(yAxis.c_str());
      h2_axes_log->GetXaxis()->SetTitleOffset(1.1);
      h2_axes_log->GetYaxis()->SetTitleOffset(1.5);
      c1->SetLogy();
      h2_axes_log->Draw("");
      dataHisto->Draw("E same");
      gPad->RedrawAxis();
      legend->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      label_cuts->Draw("same");
      std::string canvasName_log = canvasName + "_log";
      canvasName_eps = canvasName_log + ".eps";
      c1->SaveAs(canvasName_eps.c_str());
      canvasName_png = canvasName_log + ".png";
      c1->SaveAs(canvasName_png.c_str());
      canvasName_pdf = canvasName_log + ".pdf";
      if( pdf_aussi_ )
        c1->SaveAs(canvasName_pdf.c_str());
      delete h2_axes_log;
    }

      
    delete c1;
    delete legend;
    delete h2_axes;

  } //for iplots


  if( name=="response" ) {

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
    c1->cd();

    TPad* pad_hi = new TPad("pad_hi", "", 0., 0.33, 0.99, 0.99);
    pad_hi->Draw();
    pad_hi->SetLogx();
    pad_hi->SetLeftMargin(0.12);
    pad_hi->SetBottomMargin(0.015);
    TPad* pad_lo = new TPad("pad_lo", "", 0., 0., 0.99, 0.33);
    pad_lo->Draw();
    pad_lo->SetLogx();
    pad_lo->SetLeftMargin(0.12);
    pad_lo->SetTopMargin(1.);
    pad_lo->SetBottomMargin(0.3);


    pad_lo->cd();
    //shrinkPad(0.001, 0.2);

    Float_t scale_uncert = (recoType_=="calo") ? 0.1 : 0.05;

    Float_t ptPhotMax = 70.;

    //TH2* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 10., ptPhotMax, 10, (1.-3.*scale_uncert), (1.+3.*scale_uncert));
    TH2* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, ptPhot_binning[0], ptPhotMax, 10, 0.65, 1.35);
    h2_axes_lo_resp->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes_lo_resp->SetYTitle("Data / MC");
    h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.2);
    h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.55);
    h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);
    h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();
    h2_axes_lo_resp->GetXaxis()->SetNoExponent();
    //h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.);
    h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.085);
    h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);
    h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.09);
    h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.08);
    h2_axes_lo_resp->GetYaxis()->SetNdivisions(5, kTRUE);
    h2_axes_lo_resp->Draw("");
    
    TLine* line_one = new TLine( ptPhot_binning[0], 1., ptPhotMax, 1. );
    line_one->Draw("same");
    
    //TLine* line_plus_resp = new TLine( 10., 1.+scale_uncert, ptPhotMax, 1.+scale_uncert );
    TLine* line_plus_resp = new TLine( ptPhot_binning[0], 1.05, ptPhotMax, 1.05 );
    line_plus_resp->SetLineColor(46);
    line_plus_resp->SetLineWidth(2);
    //line_plus_resp->SetLineStyle(2);
    line_plus_resp->Draw("same");
    
    //TLine* line_minus_resp = new TLine( 10., 1.-scale_uncert, ptPhotMax, 1.-scale_uncert );
    TLine* line_minus_resp = new TLine( ptPhot_binning[0], 0.95, ptPhotMax, 0.95 );
    line_minus_resp->SetLineColor(46);
    line_minus_resp->SetLineWidth(2);
    //line_minus_resp->SetLineStyle(2);
    line_minus_resp->Draw("same");
    
    TLine* line_plus2_resp = new TLine( ptPhot_binning[0], 1.1, ptPhotMax, 1.1 );
    line_plus2_resp->SetLineColor(46);
    line_plus2_resp->SetLineWidth(2);
    line_plus2_resp->SetLineStyle(2);
    line_plus2_resp->Draw("same");
    
    TLine* line_minus2_resp = new TLine( ptPhot_binning[0], 0.9, ptPhotMax, 0.9 );
    line_minus2_resp->SetLineColor(46);
    line_minus2_resp->SetLineWidth(2);
    line_minus2_resp->SetLineStyle(2);
    line_minus2_resp->Draw("same");
    


    gPad->RedrawAxis();

    
    pad_hi->cd();
    //shrinkPad(0.001, 0.2);

    TH2D* h2_axes = new TH2D("axes_again", "", 10, ptPhot_binning[0], ptPhotMax, 10, 0., 1.25);
    //h2_axes->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes->SetYTitle("Response");
    //h2_axes->SetYTitle("< p_{T}^{jet} / p_{T}^{#gamma} >");
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.2);
    h2_axes->GetYaxis()->SetTitleSize(0.045);
    //h2_axes->GetXaxis()->SetMoreLogLabels();
    //h2_axes->GetXaxis()->SetNoExponent();
    h2_axes->GetXaxis()->SetLabelSize(0.);

    h2_axes->Draw();

    Float_t labelTextSize = 0.035;
    std::string jetAlgoName = getAlgoName();
    TPaveText* label_algo = new TPaveText(0.27, 0.03, 0.32, 0.10, "brNDC");
    label_algo->SetFillColor(kWhite);
    label_algo->SetTextSize(labelTextSize);
    label_algo->AddText(jetAlgoName.c_str());

    TLegend* legend = new TLegend(0.5, 0.06, 0.85, 0.3, "|#eta^{Jet}| < 1.3");
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(0);
    legend->SetTextSize(labelTextSize);
    legend->AddEntry( gr_response_vs_pt, "Data (#gamma+Jet)", "P");
    legend->Draw("same");

    Float_t cmsTextSize = 0.043;
    TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
    label_cms->SetFillColor(kWhite);
    label_cms->SetTextSize(cmsTextSize);
    label_cms->SetTextFont(62);
    label_cms->AddText("CMS Preliminary 2010");

    Float_t sqrtTextSize = 0.041;
    TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
    label_sqrt->SetFillColor(kWhite);
    label_sqrt->SetTextSize(sqrtTextSize);
    label_sqrt->SetTextFont(42);
    std::string label_sqrt_text = this->getSqrtText();
    label_sqrt->AddText(label_sqrt_text.c_str());
  
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");


    gr_response_vs_pt->SetMarkerStyle(20);
    gr_response_vs_pt->SetMarkerSize(1.8);
    gr_response_vs_pt->SetMarkerColor(kBlack);
    gr_response_vs_pt->Draw("Psame");

    gPad->RedrawAxis();

    std::string canvName;
    if( flags!="" )
      canvName = outputdir_ + "/response_vs_pt_"+flags;
    else
      canvName = outputdir_ + "/response_vs_pt";

    std::string canvName_eps = canvName + ".eps";
    c1->SaveAs(canvName_eps.c_str());
    std::string canvName_png = canvName + ".png";
    c1->SaveAs(canvName_png.c_str());
     

    // ----------------------------------------------------
    //             and now resolutions:
    // ----------------------------------------------------

    pad_hi->cd();

    TH2D* h2_axes2 = new TH2D("axes_again2", "", 10, ptPhot_binning[0], ptPhotMax, 10, 0., 1.);
    //h2_axes2->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes2->SetYTitle("Resolution");
    //h2_axes2->GetXaxis()->SetTitleOffset(1.1);
    h2_axes2->GetYaxis()->SetTitleOffset(1.5);
    //h2_axes2->GetXaxis()->SetMoreLogLabels();
    //h2_axes2->GetXaxis()->SetNoExponent();
    h2_axes2->GetXaxis()->SetLabelSize(0.);

    h2_axes2->Draw();

    TPaveText* label_cms2 = new TPaveText(0.58, 0.83, 0.75, 0.87, "brNDC");
    label_cms2->SetFillColor(kWhite);
    label_cms2->SetTextSize(cmsTextSize);
    label_cms2->SetTextFont(62);
    label_cms2->AddText("CMS Preliminary 2010");

    TPaveText* label_sqrt2 = new TPaveText(0.58, 0.78, 0.75, 0.82, "brNDC");
    label_sqrt2->SetFillColor(kWhite);
    label_sqrt2->SetTextSize(sqrtTextSize);
    label_sqrt2->SetTextFont(42);
    label_sqrt2->AddText(label_sqrt_text.c_str());
  
    TPaveText* label_algo2 = new TPaveText(0.27, 0.82, 0.32, 0.86, "brNDC");
    label_algo2->SetFillColor(kWhite);
    label_algo2->SetTextSize(labelTextSize);
    label_algo2->AddText(jetAlgoName.c_str());

    TLegend* legend2 = new TLegend(0.5, 0.5, 0.85, 0.73, "|#eta^{Jet}| < 1.3");
    legend2->SetFillColor(kWhite);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(labelTextSize);
    legend2->AddEntry( gr_resolution_vs_pt, "Data (#gamma+Jet)", "P");
    legend2->Draw("same");


    gr_resolution_vs_pt->SetMarkerStyle(20);
    gr_resolution_vs_pt->SetMarkerSize(1.8);
    gr_resolution_vs_pt->SetMarkerColor(kBlack);
    gr_resolution_vs_pt->Draw("Psame");

    label_cms2->Draw("same");
    label_sqrt2->Draw("same");
    label_algo2->Draw("same");

    gPad->RedrawAxis();

    pad_lo->cd();

    TH2* h2_axes_lo_reso = new TH2D("axes_lo_reso", "", 10, ptPhot_binning[0], ptPhotMax, 10, (1.-6.*scale_uncert), (1.+6.*scale_uncert));
    h2_axes_lo_reso->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes_lo_reso->SetYTitle("Data / MC");
    h2_axes_lo_reso->GetXaxis()->SetTitleOffset(1.2);
    h2_axes_lo_reso->GetYaxis()->SetTitleOffset(0.55);
    h2_axes_lo_reso->GetXaxis()->SetTickLength(0.06);
    h2_axes_lo_reso->GetXaxis()->SetMoreLogLabels();
    h2_axes_lo_reso->GetXaxis()->SetNoExponent();
    h2_axes_lo_reso->GetXaxis()->SetLabelSize(0.085);
    h2_axes_lo_reso->GetYaxis()->SetLabelSize(0.07);
    h2_axes_lo_reso->GetXaxis()->SetTitleSize(0.09);
    h2_axes_lo_reso->GetYaxis()->SetTitleSize(0.08);
    h2_axes_lo_reso->GetYaxis()->SetNdivisions(5, kTRUE);
    h2_axes_lo_reso->Draw("");
    
    line_one->Draw("same");

    TLine* line_plus_reso = new TLine( ptPhot_binning[0], 1.+2.*scale_uncert, ptPhotMax, 1.+2.*scale_uncert );
    line_plus_reso->SetLineColor(46);
    line_plus_reso->SetLineWidth(2);
    line_plus_reso->SetLineStyle(2);
    line_plus_reso->Draw("same");
    
    TLine* line_minus_reso = new TLine( ptPhot_binning[0], 1.-2.*scale_uncert, ptPhotMax, 1.-2.*scale_uncert );
    line_minus_reso->SetLineColor(46);
    line_minus_reso->SetLineWidth(2);
    line_minus_reso->SetLineStyle(2);
    line_minus_reso->Draw("same");
    

    gPad->RedrawAxis();

    if( flags!="" )
      canvName = outputdir_ + "/resolution_vs_pt_"+flags;
    else
      canvName = outputdir_ + "/resolution_vs_pt";

    canvName_eps = canvName + ".eps";
    c1->SaveAs(canvName_eps.c_str());
    canvName_png = canvName + ".png";
    c1->SaveAs(canvName_png.c_str());

    delete h2_axes;
    h2_axes = 0;
    delete h2_axes2;
    h2_axes2 = 0;
    delete h2_axes_lo_resp;
    h2_axes_lo_resp = 0;
    delete h2_axes_lo_reso;
    h2_axes_lo_reso = 0;
    delete c1;
    c1 = 0;

    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);


  } //if response

}



void DrawBase::drawHisto_2bkg( std::string name, std::string name_bkg1, std::string name_bkg2, std::string etaRegion, std::string flags, int legendQuadrant, bool log_aussi) {

  if( mcFile2_==0 ) {
    std::cout << "Set mcFile2 first!" << std::endl;
    exit(918);
  }


  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

  bool drawResponseGraphs = ( name=="response" || name=="responseMPF" );

  // if response will have to do a plot for each pt bin:
  int number_of_plots = (drawResponseGraphs) ? (ptPhot_binning.size()-1) : 1;

  std::string ptPhotMean_name = "ptPhotMean";
  if( flags!= "" ) ptPhotMean_name = ptPhotMean_name + "_" + flags;

  TProfile* hp_ptPhot_data = (TProfile*)dataFile_->Get(ptPhotMean_name.c_str());
  TProfile* hp_ptPhot_mc = (TProfile*)mcFile_->Get(ptPhotMean_name.c_str());
  TProfile* hp_ptPhot_mc2 = (TProfile*)mcFile2_->Get(ptPhotMean_name.c_str());
  TProfile* hp_ptJetGen_mc = (TProfile*)mcFile_->Get("ptJetGenMean");
  TProfile* hp_ptJetGen_mc2 = (TProfile*)mcFile2_->Get("ptJetGenMean");

  if( drawResponseGraphs && (hp_ptPhot_data==0 || hp_ptPhot_mc==0 || hp_ptPhot_mc2==0) ) {
    std::cout << "Didn't find '" << ptPhotMean_name << "' profile. Skipping" << std::endl;
    return;
  }

  TGraphErrors* gr_response_vs_pt = new TGraphErrors(0);
  gr_response_vs_pt->SetName("response_vs_pt");

  TGraphErrors* gr_responseMC_vs_pt = new TGraphErrors(0);
  gr_responseMC_vs_pt->SetName("responseMC_vs_pt");

  TGraphErrors* gr_responseGEN_vs_pt = new TGraphErrors(0);
  gr_responseGEN_vs_pt->SetName("responseGEN_vs_pt");

  TGraphErrors* gr_resolution_vs_pt = new TGraphErrors(0);
  gr_resolution_vs_pt->SetName("resolution_vs_pt");

  TGraphErrors* gr_resolutionMC_vs_pt = new TGraphErrors(0);
  gr_resolutionMC_vs_pt->SetName("resolutionMC_vs_pt");

  TGraphErrors* gr_resolutionGEN_vs_pt = new TGraphErrors(0);
  gr_resolutionGEN_vs_pt->SetName("resolutionGEN_vs_pt");


  for( unsigned iplot=0; iplot<number_of_plots; ++iplot ) {

    std::string histoName = name;
    if( etaRegion!="" ) histoName = histoName + "_" + etaRegion;
    if( flags!="" ) histoName = histoName + "_" + flags;

    std::string ptRange_str; 
    if( drawResponseGraphs ) {
      char ptRange[60];
      sprintf( ptRange, "ptPhot_%.0f_%.0f", ptPhot_binning[iplot], ptPhot_binning[iplot+1]);
      std::string ptRange_str_tmp(ptRange);
      ptRange_str = ptRange_str_tmp;
      histoName = histoName + "_" + ptRange_str;
    }


    TH1F* dataHisto = (TH1F*)dataFile_->Get(histoName.c_str());
    TH1F* mcHisto = (TH1F*)mcFile_->Get(histoName.c_str());
    TH1F* mcHisto2 = (TH1F*)mcFile2_->Get(histoName.c_str());

    if( dataHisto==0 || mcHisto==0 || mcHisto2==0 ) {
      std::cout << "Didn't find histo '" << histoName << "'. Continuing." << std::endl;
      return;
    }

    dataHisto->SetMarkerStyle(20);

    mcHisto->SetFillColor(38);

    mcHisto2->SetFillColor(46);

    TH1F* mcHisto_sum = new TH1F(*mcHisto);
    mcHisto_sum->Add(mcHisto2);
    mcHisto_sum->SetFillColor(46);
    if( scaleFactor_ > 0. ) {
      mcHisto_sum->Scale(scaleFactor_);
      mcHisto->Scale(scaleFactor_);
    } else {
      Float_t dataIntegral = dataHisto->Integral(0, dataHisto->GetNbinsX()+1);
      Float_t mcIntegral = mcHisto_sum->Integral(0, mcHisto_sum->GetNbinsX()+1);
      mcHisto_sum->Scale( dataIntegral/mcIntegral );
      mcHisto->Scale( dataIntegral/mcIntegral );
    }


    Float_t yAxisMaxScale = (name=="phiJet" || name=="etaJet" || name=="ptSecondJetRel" || name=="phiPhot" || name=="etaPhot" ) ? 1.8 : 1.6;
    if( name=="phiPhot" || name=="etaPhot" ) yAxisMaxScale=2.;
    if(name=="clusterMajPhotReco" || name=="clusterMinPhotReco") yAxisMaxScale = 2.;
    Float_t xMin = dataHisto->GetXaxis()->GetXmin();
    Float_t xMax = dataHisto->GetXaxis()->GetXmax();
    Float_t yMax_data = dataHisto->GetMaximum();
////Float_t yMax_mc = mcStack->GetMaximum();
    Float_t yMax_mc = mcHisto_sum->GetMaximum();
    Float_t yMax = (yMax_data>yMax_mc) ? yAxisMaxScale*yMax_data : yAxisMaxScale*yMax_mc;
    Float_t yMin = 0.;


    std::string xAxis = getAxisName(name);

    std::string instanceName = (analysisType_=="PhotonJet") ? "Events" : "Jets";
    std::string yAxis = instanceName;
    if( name=="ptJet" || name=="ptCorrJet"||name=="ptPhot" ) {
      char yAxis_char[50];
      sprintf(yAxis_char, "%s/(%d GeV/c)", instanceName.c_str(), (Int_t)dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  

    if( name=="EchJet" || name=="EgammaJet" || name=="EnhJet" ) {
      char yAxis_char[50];
      if( dataHisto->GetBinWidth(1) < 1. ) sprintf(yAxis_char, "%s/(%.1f GeV)", instanceName.c_str(), dataHisto->GetBinWidth(1));
      else                                 sprintf(yAxis_char, "%s/(%.0f GeV)", instanceName.c_str(), dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  


    if( name=="diJetMass" ) {
      char yAxis_char[50];
      sprintf(yAxis_char, "Events/(%.2f GeV/c^{2})", dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  
    if( name=="massJet" ) {
      char yAxis_char[50];
      sprintf(yAxis_char, "Jets/(%.2f GeV/c ^{2})", dataHisto->GetBinWidth(1));
      //sprintf(yAxis_char, "Jets/(%d GeV/c^{2})", (Int_t)dataHisto->GetBinWidth(1));
      std::string yAxis_tmp(yAxis_char);
      yAxis=yAxis_tmp;
    }  
    if( name=="deltaPhiJet" || name=="asymmJet" ) yAxis = "Events";

    std::string etaRange = getEtaRangeText( etaRegion );

    TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax);
    h2_axes->SetXTitle(xAxis.c_str());
    h2_axes->SetYTitle(yAxis.c_str());
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.5);

    if( legendQuadrant<0 || legendQuadrant>5 ) {
      std::cout << "Invalid legendQuadrant! Exiting!" << std::endl;
      exit(1);
    }

    Float_t legend_xMin;
    Float_t legend_yMin;
    Float_t legend_xMax;
    Float_t legend_yMax;

    if( legendQuadrant==1 ) {
      legend_xMin = 0.6;
      legend_yMin = 0.73;
      legend_xMax = 0.88;
      legend_yMax = 0.88;
    } else if( legendQuadrant==0 ) {
      legend_xMin = 0.5;
      legend_yMin = 0.73;
      legend_xMax = 0.73;
      legend_yMax = 0.88;
    } else if( legendQuadrant==2 ) {
      legend_xMin = 0.18;
      legend_yMin = 0.73;
      legend_xMax = 0.41;
      legend_yMax = 0.88;
    } else if( legendQuadrant==5 ) {
      legend_xMin = 0.4;
      legend_yMin = 0.15;
      legend_xMax = 0.6;
      legend_yMax = 0.25;
    }


    TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax);
    legend->SetFillColor(kWhite);
    legend->SetTextSize(0.035);
    legend->AddEntry(dataHisto, "Data", "P");
    legend->AddEntry(mcHisto, name_bkg1.c_str(), "F");
    legend->AddEntry(mcHisto2, name_bkg2.c_str(), "F");

    TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
    label_cms->SetFillColor(kWhite);
    label_cms->SetTextSize(0.038);
    label_cms->SetTextFont(62);
    label_cms->AddText("CMS Preliminary 2010");

    TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
    label_sqrt->SetFillColor(kWhite);
    label_sqrt->SetTextSize(0.038);
    label_sqrt->SetTextFont(42);
    std::string label_sqrt_text = this->getSqrtText();
    label_sqrt->AddText(label_sqrt_text.c_str());

    TPaveText* label_cuts;
    if( analysisType_ == "MinBias" ) {

      Float_t label_cuts_xMin = 0.63;
      Float_t label_cuts_yMin = 0.5;
      Float_t label_cuts_xMax = 0.84;
      Float_t label_cuts_yMax = 0.65;

      if( name=="asymmJet" || name=="deltaPhiJet" || name=="RchJet" || name=="etaJet" || name=="phiJet" ) {
       label_cuts_xMin = 0.25;
       label_cuts_yMin = 0.60;
       label_cuts_xMax = 0.36;
       label_cuts_yMax = 0.73;
      }
      if( name=="RnhJet" ) {
       label_cuts_xMin = 0.4;
       label_cuts_yMin = 0.6;
       label_cuts_xMax = 0.6;
       label_cuts_yMax = 0.7;
      }

      label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax,  "brNDC");
      label_cuts->SetFillColor(kWhite);
      label_cuts->SetTextSize(0.035);
      label_cuts->SetTextFont(42);
      std::string jetAlgoName = getAlgoName();
      label_cuts->AddText(jetAlgoName.c_str());
      if( name != "etaJet" )
        label_cuts->AddText(etaRange.c_str());
      if( name != "ptJet" && name != "ptCorrJet" ) {
        char labelText[70];
        sprintf( labelText, "p_{T}^{%s} > %d GeV/c", raw_corr_.c_str(), pt_thresh_);
        label_cuts->AddText(labelText);
      }

    } else if( analysisType_=="PhotonJet" ) {

      Float_t label_cuts_xMin = 0.63;
      Float_t label_cuts_yMin = 0.53;
      Float_t label_cuts_xMax = 0.84;
      Float_t label_cuts_yMax = 0.68;

      if( name=="deltaPhi" ) {
       label_cuts_xMin = 0.25;
       label_cuts_yMin = 0.60;
       label_cuts_xMax = 0.36;
       label_cuts_yMax = 0.73;
      }
      label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax,  "brNDC");
      label_cuts->SetFillColor(kWhite);
      label_cuts->SetTextSize(0.035);
      label_cuts->SetTextFont(42);
      
      if( name != "ptPhot" ) {
        char ptCutName[60];
        if( drawResponseGraphs )
          sprintf( ptCutName, "%.0f < p_{T}^{#gamma} < %.0f GeV/c", ptPhot_binning[iplot], ptPhot_binning[iplot+1] );
        else
          sprintf( ptCutName, "p_{T}^{#gamma} > %.0f GeV/c", ptPhot_binning[0]);

        label_cuts->AddText(ptCutName);
      }

      if( name != "etaPhot" )
        label_cuts->AddText("|#eta_{#gamma}| < 1.3");

      if( drawResponseGraphs||name=="deltaPhi"||name=="ptSecondJetRel" ) {
        std::string jetAlgoName = getAlgoName();
        label_cuts->AddText(jetAlgoName.c_str());
      }

    } else {

      std::cout << "AnalysisType: '" << analysisType_ << "' not yet implemented. Exiting." << std::endl;
      exit(991);

    }

    TPaveText* label_mean = new TPaveText(0.63, 0.37, 0.84, 0.5, "brNDC");
    label_mean->SetFillColor(kWhite);
    label_mean->SetTextSize(0.035);
    label_mean->SetTextFont(42);
    label_mean->AddText("Mean:");
    char dataMean[100];
//    sprintf( dataMean, "Data Mean = %.2f #pm %.2f", dataHisto->GetMean(), dataHisto->GetMeanError());
    sprintf( dataMean, "%.2f #pm %.2f (DATA)", dataHisto->GetMean(), dataHisto->GetMeanError());
    label_mean->AddText(dataMean);
    char mcMean[100];
//    sprintf( mcMean, "MC Mean = %.2f #pm %.2f", mcHisto_sum->GetMean(), mcHisto_sum->GetMeanError());
    sprintf( mcMean, "%.2f #pm %.2f (MC)", mcHisto_sum->GetMean(), mcHisto_sum->GetMeanError());
    label_mean->AddText(mcMean);


    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    c1->cd();

    h2_axes->Draw("");
    mcHisto_sum->Draw("histo same");
    mcHisto->Draw("histo same");
    dataHisto->Draw("E same");
    gPad->RedrawAxis();
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cuts->Draw("same");

    if( drawResponseGraphs ) {

      //label_mean->Draw("same");

      Float_t dataResponse = dataHisto->GetMean();
      Float_t dataResponseErr = dataHisto->GetMeanError();
      Float_t dataRMS = dataHisto->GetRMS();
      Float_t dataRMSErr = dataHisto->GetMeanError();
      Float_t dataResolution = dataRMS / dataResponse;
      Float_t dataResolutionErr = sqrt( dataRMSErr*dataRMSErr/(dataResponse*dataResponse) + dataResolution*dataResolution*dataResponseErr*dataResponseErr/(dataResponse*dataResponse*dataResponse*dataResponse) );

      gr_response_vs_pt->SetPoint( iplot, hp_ptPhot_data->GetBinContent( iplot+1 ), dataResponse );
      gr_response_vs_pt->SetPointError( iplot, hp_ptPhot_data->GetBinError( iplot+1 ), dataResponseErr );

      gr_resolution_vs_pt->SetPoint( iplot, hp_ptPhot_data->GetBinContent( iplot+1 ), dataResolution );
      gr_resolution_vs_pt->SetPointError( iplot, hp_ptPhot_data->GetBinError( iplot+1 ), dataResolutionErr );

      //weighted average of two MC pt's
      Float_t xmc1 = hp_ptPhot_mc->GetBinContent(iplot+1);
      Float_t xmc2 = hp_ptPhot_mc2->GetBinContent(iplot+1);
      Float_t errmc1 = hp_ptPhot_mc->GetBinError(iplot+1);
      Float_t errmc2 = hp_ptPhot_mc2->GetBinError(iplot+1);
      Float_t denommc = ( (1./errmc1) + (1./errmc2) );
      Float_t ptMeanMC = ( (xmc1/errmc1) + (xmc2/errmc2) ) / denommc;
      Float_t ptMeanErrMC = 1./denommc;
     
      Float_t mcResponse = mcHisto_sum->GetMean();
      Float_t mcResponseErr = mcHisto_sum->GetMeanError();
      Float_t mcRMS = mcHisto_sum->GetRMS();
      Float_t mcRMSErr = mcHisto_sum->GetMeanError();
      Float_t mcResolution = mcRMS / mcResponse;
      Float_t mcResolutionErr = sqrt( mcRMSErr*mcRMSErr/(mcResponse*mcResponse) + mcResolution*mcResolution*mcResponseErr*mcResponseErr/(mcResponse*mcResponse*mcResponse*mcResponse) );

      gr_responseMC_vs_pt->SetPoint( iplot, ptMeanMC, mcResponse );
      gr_responseMC_vs_pt->SetPointError( iplot, ptMeanErrMC, mcResponseErr );

      gr_resolutionMC_vs_pt->SetPoint( iplot, ptMeanMC, mcResolution );
      gr_resolutionMC_vs_pt->SetPointError( iplot, ptMeanErrMC, mcResolutionErr );

   
      char responseGEN_name[100];
      //if( flags!="" )
      //  sprintf( responseGEN_name, "%sGEN_%s_%s", name.c_str(), flags.c_str(), ptRange_str.c_str());
      //else
        sprintf( responseGEN_name, "responseGEN_%s", ptRange_str.c_str());

      TH1F* responseGEN = (TH1F*)mcFile_->Get(responseGEN_name);
      TH1F* responseGEN2 = (TH1F*)mcFile2_->Get(responseGEN_name);

      TH1F* responseGEN_sum = new TH1F(*responseGEN);
      responseGEN_sum->Add(responseGEN2);
      responseGEN_sum->Scale(scaleFactor_);
      //responseGEN_sum->SetLineStyle(3);
      responseGEN_sum->SetLineWidth(3);
      //responseGEN_sum->Draw("histo C same");
      //legend->AddEntry(responseGEN_sum, "Intrinsic", "L");

      Float_t genResponse = responseGEN_sum->GetMean();
      Float_t genResponseErr = responseGEN_sum->GetMeanError();
      Float_t genRMS = responseGEN_sum->GetRMS();
      Float_t genRMSErr = responseGEN_sum->GetMeanError();

    ////compute responseGEN only from QCD BG:
    //Float_t genResponse = responseGEN->GetMean();
    //Float_t genResponseErr = responseGEN->GetMeanError();
    //Float_t genRMS = responseGEN->GetRMS();
    //Float_t genRMSErr = responseGEN->GetMeanError();

    ////compute responseGEN only from Photon+Jet
    //Float_t genResponse = responseGEN2->GetMean();
    //Float_t genResponseErr = responseGEN2->GetMeanError();
    //Float_t genRMS = responseGEN2->GetRMS();
    //Float_t genRMSErr = responseGEN2->GetMeanError();

      Float_t genResolution = genRMS / genResponse;
      Float_t genResolutionErr = sqrt( genRMSErr*genRMSErr/(genResponse*genResponse) + genResolution*genResolution*genResponseErr*genResponseErr/(genResponse*genResponse*genResponse*genResponse) );

      Float_t xgen1 = hp_ptJetGen_mc->GetBinContent(iplot+1);
      Float_t xgen2 = hp_ptJetGen_mc2->GetBinContent(iplot+1);
      Float_t errgen1 = hp_ptJetGen_mc->GetBinError(iplot+1);
      Float_t errgen2 = hp_ptJetGen_mc2->GetBinError(iplot+1);
    //Float_t xgen1 = hp_ptPhot_mc->GetBinContent(iplot+1);
    //Float_t xgen2 = hp_ptPhot_mc2->GetBinContent(iplot+1);
    //Float_t errgen1 = hp_ptPhot_mc->GetBinError(iplot+1);
    //Float_t errgen2 = hp_ptPhot_mc2->GetBinError(iplot+1);
      Float_t denomgen = ( (1./errgen1) + (1./errgen2) );
      Float_t ptMeanGEN = ( (xgen1/errgen1) + (xgen2/errgen2) ) / denomgen;
      Float_t ptMeanErrGEN = 1./denomgen;
      
      gr_responseGEN_vs_pt->SetPoint( iplot, ptMeanGEN, genResponse );
      gr_responseGEN_vs_pt->SetPointError( iplot, ptMeanErrGEN, genResponseErr );

      gr_resolutionGEN_vs_pt->SetPoint( iplot, ptMeanGEN, genResolution );
      gr_resolutionGEN_vs_pt->SetPointError( iplot, ptMeanErrGEN, genResolutionErr );


    } //if reponse

std::cout << "fin qui OK. " << std::endl;
    legend->Draw("same");
      


    std::string canvasName = outputdir_ + "/" + name;
    if( etaRegion!="" )
      canvasName = canvasName + "_" + etaRegion;
    if( flags!="" )
      canvasName = canvasName + "_" + flags;

    if( drawResponseGraphs )
      canvasName = canvasName + "_" + ptRange_str;

    std::string canvasName_eps = canvasName + ".eps";
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);
    c1->SaveAs(canvasName_eps.c_str());
    std::string canvasName_png = canvasName + ".png";
    c1->SaveAs(canvasName_png.c_str());
    std::string canvasName_pdf = canvasName + ".pdf";
    if( pdf_aussi_ )
      c1->SaveAs(canvasName_pdf.c_str());

    if( log_aussi ) {
      TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, 0.5, 5.*yMax);
      h2_axes_log->SetXTitle(xAxis.c_str());
      h2_axes_log->SetYTitle(yAxis.c_str());
      h2_axes_log->GetXaxis()->SetTitleOffset(1.1);
      h2_axes_log->GetYaxis()->SetTitleOffset(1.5);
      c1->SetLogy();
      h2_axes_log->Draw("");
      mcHisto_sum->Draw("histo same");
      mcHisto->Draw("histo same");
      dataHisto->Draw("E same");
      gPad->RedrawAxis();
      legend->Draw("same");
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      label_cuts->Draw("same");
      std::string canvasName_log = canvasName + "_log";
      canvasName_eps = canvasName_log + ".eps";
      c1->SaveAs(canvasName_eps.c_str());
      canvasName_png = canvasName_log + ".png";
      c1->SaveAs(canvasName_png.c_str());
      canvasName_pdf = canvasName_log + ".pdf";
      if( pdf_aussi_ )
        c1->SaveAs(canvasName_pdf.c_str());
      delete h2_axes_log;
    }

      
    delete c1;
    delete legend;
    delete h2_axes;

  } //for n plots


  if( drawResponseGraphs ) {

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
    c1->cd();

    TPad* pad_hi = new TPad("pad_hi", "", 0., 0.33, 0.99, 0.99);
    pad_hi->Draw();
    pad_hi->SetLogx();
    pad_hi->SetLeftMargin(0.12);
    pad_hi->SetBottomMargin(0.015);
    TPad* pad_lo = new TPad("pad_lo", "", 0., 0., 0.99, 0.33);
    pad_lo->Draw();
    pad_lo->SetLogx();
    pad_lo->SetLeftMargin(0.12);
    pad_lo->SetTopMargin(1.);
    pad_lo->SetBottomMargin(0.3);


    pad_lo->cd();
    //shrinkPad(0.001, 0.2);

    Float_t scale_uncert = (recoType_=="calo") ? 0.1 : 0.05;

    Float_t ptPhotMax = 70.;

    //TH2* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 10., ptPhotMax, 10, (1.-3.*scale_uncert), (1.+3.*scale_uncert));
    TH2* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, ptPhot_binning[0], ptPhotMax, 10, 0.65, 1.35);
    h2_axes_lo_resp->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes_lo_resp->SetYTitle("Data / MC");
    h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.2);
    h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.55);
    h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);
    h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();
    h2_axes_lo_resp->GetXaxis()->SetNoExponent();
    //h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.);
    h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.085);
    h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);
    h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.09);
    h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.08);
    h2_axes_lo_resp->GetYaxis()->SetNdivisions(5, kTRUE);
    h2_axes_lo_resp->Draw("");
    
    TLine* line_one = new TLine( ptPhot_binning[0], 1., ptPhotMax, 1. );
    line_one->Draw("same");
    
    //TLine* line_plus_resp = new TLine( 10., 1.+scale_uncert, ptPhotMax, 1.+scale_uncert );
    TLine* line_plus_resp = new TLine( ptPhot_binning[0], 1.05, ptPhotMax, 1.05 );
    line_plus_resp->SetLineColor(46);
    line_plus_resp->SetLineWidth(2);
    //line_plus_resp->SetLineStyle(2);
    line_plus_resp->Draw("same");
    
    //TLine* line_minus_resp = new TLine( 10., 1.-scale_uncert, ptPhotMax, 1.-scale_uncert );
    TLine* line_minus_resp = new TLine( ptPhot_binning[0], 0.95, ptPhotMax, 0.95 );
    line_minus_resp->SetLineColor(46);
    line_minus_resp->SetLineWidth(2);
    //line_minus_resp->SetLineStyle(2);
    line_minus_resp->Draw("same");
    
    TLine* line_plus2_resp = new TLine( ptPhot_binning[0], 1.1, ptPhotMax, 1.1 );
    line_plus2_resp->SetLineColor(46);
    line_plus2_resp->SetLineWidth(2);
    line_plus2_resp->SetLineStyle(2);
    line_plus2_resp->Draw("same");
    
    TLine* line_minus2_resp = new TLine( ptPhot_binning[0], 0.9, ptPhotMax, 0.9 );
    line_minus2_resp->SetLineColor(46);
    line_minus2_resp->SetLineWidth(2);
    line_minus2_resp->SetLineStyle(2);
    line_minus2_resp->Draw("same");
    

    TGraphErrors* gr_resp_ratio = this->getGraphRatio( gr_response_vs_pt, gr_responseMC_vs_pt );
    gr_resp_ratio->SetName("response_ratio");
    gr_resp_ratio->SetMarkerStyle(20);
    gr_resp_ratio->SetMarkerSize(1.8);
    gr_resp_ratio->Draw("P");

    gPad->RedrawAxis();

    
    pad_hi->cd();
    //shrinkPad(0.001, 0.2);

    TH2D* h2_axes = new TH2D("axes_again", "", 10, ptPhot_binning[0], ptPhotMax, 10, 0., 1.25);
    //h2_axes->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes->SetYTitle("Response");
    //h2_axes->SetYTitle("< p_{T}^{jet} / p_{T}^{#gamma} >");
    h2_axes->GetXaxis()->SetTitleOffset(1.1);
    h2_axes->GetYaxis()->SetTitleOffset(1.2);
    h2_axes->GetYaxis()->SetTitleSize(0.045);
    //h2_axes->GetXaxis()->SetMoreLogLabels();
    //h2_axes->GetXaxis()->SetNoExponent();
    h2_axes->GetXaxis()->SetLabelSize(0.);

    h2_axes->Draw();

    Float_t labelTextSize = 0.035;
    std::string jetAlgoName = getAlgoName();
    TPaveText* label_algo = new TPaveText(0.27, 0.03, 0.32, 0.10, "brNDC");
    label_algo->SetFillColor(kWhite);
    label_algo->SetTextSize(labelTextSize);
    label_algo->AddText(jetAlgoName.c_str());

    TLegend* legend = new TLegend(0.5, 0.06, 0.85, 0.3, "|#eta^{Jet}| < 1.3");
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(0);
    legend->SetTextSize(labelTextSize);
    if( name=="response" ) {
      legend->AddEntry( gr_response_vs_pt, "Data (#gamma+Jet)", "P");
      legend->AddEntry( gr_responseMC_vs_pt, "Simulation (#gamma+Jet)", "P");
    } else if( name=="responseMPF" ) {
      legend->AddEntry( gr_response_vs_pt, "Data (MPF)", "P");
      legend->AddEntry( gr_responseMC_vs_pt, "Simulation (MPF)", "P");
    }
    legend->AddEntry( gr_responseGEN_vs_pt, "True Response", "P");
    legend->Draw("same");

    Float_t cmsTextSize = 0.043;
    TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
    label_cms->SetFillColor(kWhite);
    label_cms->SetTextSize(cmsTextSize);
    label_cms->SetTextFont(62);
    label_cms->AddText("CMS Preliminary 2010");

    Float_t sqrtTextSize = 0.041;
    TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
    label_sqrt->SetFillColor(kWhite);
    label_sqrt->SetTextSize(sqrtTextSize);
    label_sqrt->SetTextFont(42);
    std::string label_sqrt_text = this->getSqrtText();
    label_sqrt->AddText(label_sqrt_text.c_str());
  
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");


    gr_responseGEN_vs_pt->SetMarkerStyle(29);
    gr_responseGEN_vs_pt->SetMarkerColor(46);
    gr_responseGEN_vs_pt->SetMarkerSize(2.);
    gr_responseGEN_vs_pt->Draw("Psame");

    gr_responseMC_vs_pt->SetMarkerStyle(24);
    gr_responseMC_vs_pt->SetMarkerSize(1.8);
    gr_responseMC_vs_pt->SetMarkerColor(kBlack);
    gr_responseMC_vs_pt->SetLineColor(kBlack);
    gr_responseMC_vs_pt->Draw("Psame");

    gr_response_vs_pt->SetMarkerStyle(20);
    gr_response_vs_pt->SetMarkerSize(1.8);
    gr_response_vs_pt->SetMarkerColor(kBlack);
    gr_response_vs_pt->Draw("Psame");

    gPad->RedrawAxis();

    std::string mpf = (name=="responseMPF") ? "MPF" : "";

    std::string canvName;
    if( flags!="" )
      canvName = outputdir_ + "/response" + mpf + "_vs_pt_"+flags;
    else
      canvName = outputdir_ + "/response" + mpf + "_vs_pt";

    std::string canvName_eps = canvName + ".eps";
    c1->SaveAs(canvName_eps.c_str());
    std::string canvName_png = canvName + ".png";
    c1->SaveAs(canvName_png.c_str());
     

    // ----------------------------------------------------
    //             and now resolutions:
    // ----------------------------------------------------

    pad_hi->cd();

    TH2D* h2_axes2 = new TH2D("axes_again2", "", 10, ptPhot_binning[0], ptPhotMax, 10, 0., 1.);
    //h2_axes2->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes2->SetYTitle("Resolution");
    //h2_axes2->GetXaxis()->SetTitleOffset(1.1);
    h2_axes2->GetYaxis()->SetTitleOffset(1.5);
    //h2_axes2->GetXaxis()->SetMoreLogLabels();
    //h2_axes2->GetXaxis()->SetNoExponent();
    h2_axes2->GetXaxis()->SetLabelSize(0.);

    h2_axes2->Draw();

    TPaveText* label_cms2 = new TPaveText(0.58, 0.83, 0.75, 0.87, "brNDC");
    label_cms2->SetFillColor(kWhite);
    label_cms2->SetTextSize(cmsTextSize);
    label_cms2->SetTextFont(62);
    label_cms2->AddText("CMS Preliminary 2010");

    TPaveText* label_sqrt2 = new TPaveText(0.58, 0.78, 0.75, 0.82, "brNDC");
    label_sqrt2->SetFillColor(kWhite);
    label_sqrt2->SetTextSize(sqrtTextSize);
    label_sqrt2->SetTextFont(42);
    label_sqrt2->AddText(label_sqrt_text.c_str());
  
    TPaveText* label_algo2 = new TPaveText(0.27, 0.82, 0.32, 0.86, "brNDC");
    label_algo2->SetFillColor(kWhite);
    label_algo2->SetTextSize(labelTextSize);
    label_algo2->AddText(jetAlgoName.c_str());

    TLegend* legend2 = new TLegend(0.5, 0.5, 0.85, 0.73, "|#eta^{Jet}| < 1.3");
    legend2->SetFillColor(kWhite);
    legend2->SetFillStyle(0);
    legend2->SetTextSize(labelTextSize);
    if( name=="response" ) {
      legend2->AddEntry( gr_resolution_vs_pt, "Data (#gamma+Jet)", "P");
      legend2->AddEntry( gr_resolutionMC_vs_pt, "Simulation (#gamma+Jet)", "P");
    } else if( name == "responseMPF" ) {
      legend2->AddEntry( gr_resolution_vs_pt, "Data (#gamma+Jet, MPF)", "P");
      legend2->AddEntry( gr_resolutionMC_vs_pt, "Simulation (#gamma+Jet, MPF)", "P");
    }
    legend2->AddEntry( gr_resolutionGEN_vs_pt, "True Resolution", "P");
    legend2->Draw("same");

    gr_resolutionGEN_vs_pt->SetMarkerStyle(29);
    gr_resolutionGEN_vs_pt->SetMarkerColor(46);
    gr_resolutionGEN_vs_pt->SetMarkerSize(2.);
    gr_resolutionGEN_vs_pt->Draw("Psame");

    gr_resolutionMC_vs_pt->SetMarkerStyle(24);
    gr_resolutionMC_vs_pt->SetMarkerSize(1.8);
    gr_resolutionMC_vs_pt->SetMarkerColor(kBlack);
    gr_resolutionMC_vs_pt->SetLineColor(kBlack);
    gr_resolutionMC_vs_pt->Draw("Psame");

    gr_resolution_vs_pt->SetMarkerStyle(20);
    gr_resolution_vs_pt->SetMarkerSize(1.8);
    gr_resolution_vs_pt->SetMarkerColor(kBlack);
    gr_resolution_vs_pt->Draw("Psame");

    label_cms2->Draw("same");
    label_sqrt2->Draw("same");
    label_algo2->Draw("same");

    gPad->RedrawAxis();

    pad_lo->cd();

    TH2* h2_axes_lo_reso = new TH2D("axes_lo_reso", "", 10, ptPhot_binning[0], ptPhotMax, 10, (1.-6.*scale_uncert), (1.+6.*scale_uncert));
    h2_axes_lo_reso->SetXTitle("Photon p_{T} [GeV/c]");
    h2_axes_lo_reso->SetYTitle("Data / MC");
    h2_axes_lo_reso->GetXaxis()->SetTitleOffset(1.2);
    h2_axes_lo_reso->GetYaxis()->SetTitleOffset(0.55);
    h2_axes_lo_reso->GetXaxis()->SetTickLength(0.06);
    h2_axes_lo_reso->GetXaxis()->SetMoreLogLabels();
    h2_axes_lo_reso->GetXaxis()->SetNoExponent();
    h2_axes_lo_reso->GetXaxis()->SetLabelSize(0.085);
    h2_axes_lo_reso->GetYaxis()->SetLabelSize(0.07);
    h2_axes_lo_reso->GetXaxis()->SetTitleSize(0.09);
    h2_axes_lo_reso->GetYaxis()->SetTitleSize(0.08);
    h2_axes_lo_reso->GetYaxis()->SetNdivisions(5, kTRUE);
    h2_axes_lo_reso->Draw("");
    
    line_one->Draw("same");

    TLine* line_plus_reso = new TLine( ptPhot_binning[0], 1.+2.*scale_uncert, ptPhotMax, 1.+2.*scale_uncert );
    line_plus_reso->SetLineColor(46);
    line_plus_reso->SetLineWidth(2);
    line_plus_reso->SetLineStyle(2);
    line_plus_reso->Draw("same");
    
    TLine* line_minus_reso = new TLine( ptPhot_binning[0], 1.-2.*scale_uncert, ptPhotMax, 1.-2.*scale_uncert );
    line_minus_reso->SetLineColor(46);
    line_minus_reso->SetLineWidth(2);
    line_minus_reso->SetLineStyle(2);
    line_minus_reso->Draw("same");
    
    TGraphErrors* gr_reso_ratio = this->getGraphRatio( gr_resolution_vs_pt, gr_resolutionMC_vs_pt );
    gr_reso_ratio->SetName("reso_ratio");
    gr_reso_ratio->SetMarkerStyle(20);
    gr_reso_ratio->SetMarkerSize(1.8);
    gr_reso_ratio->Draw("P");

    gPad->RedrawAxis();

    if( flags!="" )
      canvName = outputdir_ + "/resolution" + mpf + "_vs_pt_"+flags;
    else
      canvName = outputdir_ + "/resolution" + mpf + "_vs_pt";

    canvName_eps = canvName + ".eps";
    c1->SaveAs(canvName_eps.c_str());
    canvName_png = canvName + ".png";
    c1->SaveAs(canvName_png.c_str());

    delete h2_axes;
    h2_axes = 0;
    delete h2_axes2;
    h2_axes2 = 0;
    delete h2_axes_lo_resp;
    h2_axes_lo_resp = 0;
    delete h2_axes_lo_reso;
    h2_axes_lo_reso = 0;
    delete c1;
    c1 = 0;

    gStyle->SetPadTickX(0);
    gStyle->SetPadTickY(0);


  } //if response


} //drawHisto_2bkg





void DrawBase::drawProfile( std::string yVar, std::string xVar, int legendQuadrant) {

  std::string name = yVar + "_vs_" + xVar;
  if( xVar == "pt" || xVar == "ptCorr" ) name = name + "_barrel"; //ugly fix for now

  TProfile* dataProfile = (TProfile*)dataFile_->Get(name.c_str());
  TProfile* mcProfile = (TProfile*)mcFile_->Get(name.c_str());

  if( dataProfile==0 || mcProfile==0 ) {
    std::cout << "Didn't find profile '" << name << "'. Continuing." << std::endl;
    return;
  }

  Float_t profile_xMin = dataProfile->GetXaxis()->GetXmin();
  Float_t profile_xMax = dataProfile->GetXaxis()->GetXmax();

  mcProfile->SetFillColor(38);
//mcProfile->SetFillColor(kRed-7);
//mcProfile->SetLineWidth(2);
//mcProfile->SetLineColor(kRed);

  dataProfile->SetMarkerStyle(20);

  //Float_t etamax__rounded = (etamax_>2.5) ? 3. : 2.5;
  Float_t etamax_rounded = 3.;
  Float_t xMin = (xVar=="eta") ? -etamax_rounded : profile_xMin;
  Float_t xMax = (xVar=="eta") ?  etamax_rounded : profile_xMax;
//Float_t xMin = profile_xMin;
//Float_t xMax = profile_xMax;

  if( xVar=="Nch" ) xMax = 11.5;
  
  Float_t dataMax = dataProfile->GetMaximum();
  Float_t mcMax = mcProfile->GetMaximum();
  Float_t plotMax = (dataMax>mcMax) ? dataMax : mcMax;

  std::string xAxisName = getAxisName(xVar);
  std::string yAxisName = getAxisName(yVar);

  Float_t yAxisMaxScale = 1.5;
  if( yVar=="pt" || yVar=="ptCorr" || yVar=="Rch" || yVar=="PTch" || yVar=="Nch" || yVar=="Ngamma" || yVar=="Nnh" || yVar=="Rnh" || yVar=="Rgamma" || yVar=="Rgammanh" ) yAxisMaxScale=1.8;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yAxisMaxScale*plotMax);
  h2_axes->SetXTitle(xAxisName.c_str());
  h2_axes->SetYTitle(yAxisName.c_str());
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);


  if( legendQuadrant<0 || legendQuadrant>5 ) {
    std::cout << "Invalid legendQuadrant! Exiting!" << std::endl;
    exit(1);
  }

  Float_t legend_xMin;
  Float_t legend_yMin;
  Float_t legend_xMax;
  Float_t legend_yMax;

  if( legendQuadrant==1 ) {
    legend_xMin = 0.6;
    legend_yMin = 0.73;
    legend_xMax = 0.88;
    legend_yMax = 0.88;
  } else if( legendQuadrant==2 ) {
    legend_xMin = 0.20;
    legend_yMin = 0.7;
    legend_xMax = 0.41;
    legend_yMax = 0.88;
  } else if( legendQuadrant==5 ) {
    legend_xMin = 0.4;
    legend_yMin = 0.15;
    legend_xMax = 0.6;
    legend_yMax = 0.25;
  }


  TLegend* legend = new TLegend(legend_xMin, legend_yMin, legend_xMax, legend_yMax);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.035);
  legend->AddEntry(dataProfile, "Data", "P");
  legend->AddEntry(mcProfile, "Simulation", "F");

  TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.038);
  label_cms->SetTextFont(62);
  label_cms->AddText("CMS Preliminary 2010");

  TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  std::string label_sqrt_text = this->getSqrtText();
  label_sqrt->AddText(label_sqrt_text.c_str());

  Float_t label_cuts_xMin = 0.4;
  Float_t label_cuts_yMin = 0.55;
  Float_t label_cuts_xMax = 0.6;
  Float_t label_cuts_yMax = 0.7;

//if( yVar=="pt" || yVar=="ptCorr" || yVar=="Rch" || yVar=="Rgamma" || yVar=="ETch" ) {
//  label_cuts_xMin = 0.4;
//  label_cuts_yMin = 0.65;
//  label_cuts_xMax = 0.6;
//  label_cuts_yMax = 0.8;
//}

  if( xVar=="pt" || xVar=="ptCorr" ) {
    label_cuts_xMin = 0.20;
    label_cuts_xMax = 0.35;
  }



  TPaveText* label_cuts = new TPaveText(label_cuts_xMin, label_cuts_yMin, label_cuts_xMax, label_cuts_yMax, "brNDC");
  label_cuts->SetFillColor(kWhite);
  label_cuts->SetTextSize(0.035);
  label_cuts->SetTextFont(42);
  label_cuts->AddText("Anti-k_{T} 0.5 PFJets");
  if( xVar != "eta" ) {
    char etaRange_ch[100];
    sprintf(etaRange_ch, "|#eta| < %.1f", etamax_);
    label_cuts->AddText(etaRange_ch);
  }
  if( yVar != "pt" && yVar != "ptCorr" ) {
    char labelText[70];
    sprintf( labelText, "p_{T}^{%s} > %d GeV/c", raw_corr_.c_str(), pt_thresh_);
    label_cuts->AddText(labelText);
  }


  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  //if(name=="ptJet") c1->SetLogy();
  if(xVar=="pt"||xVar=="ptCorr") {
    c1->SetLogx();
    h2_axes->GetXaxis()->SetMoreLogLabels();
    h2_axes->GetXaxis()->SetNoExponent();
  }
  h2_axes->Draw("");
  mcProfile->Draw("histo same");
  dataProfile->Draw("E same");
  gPad->RedrawAxis();
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_cuts->Draw("same");
 
  std::string canvasName = outputdir_ + "/" + name;
  std::string canvasName_eps = canvasName + ".eps";
  c1->SaveAs(canvasName_eps.c_str());
  std::string canvasName_png = canvasName + ".png";
  c1->SaveAs(canvasName_png.c_str());
  std::string canvasName_pdf = canvasName + ".pdf";
    if( pdf_aussi_ )
  c1->SaveAs(canvasName_pdf.c_str());
 

  delete c1;
  delete legend;
  delete h2_axes;

} //drawProfile



std::string DrawBase::getAxisName(std::string name) { 
  
  std::string axisName = "";

  if( name=="ptJet"|| name=="pt" ) axisName = "Jet p_{T}^{raw} [GeV/c]";
  if( name=="ptCorrJet"|| name=="ptCorr" ) axisName = "Corrected p_{T} [GeV/c]";
  if( name=="etaJet"|| name=="eta" ) axisName = "Jet #eta";
  if( name=="phiJet"|| name=="phi" ) axisName = "Jet #phi";
  if( name=="nCandJet" ) axisName = "N Candidates";
  if( name=="RchJet"|| name=="Rch" ) axisName = "R_{ch}";
  if( name=="RnhJet"|| name=="Rnh" ) axisName = "R_{nh}";
  if( name=="RgammaJet"|| name=="Rgamma" ) axisName = "R_{#gamma}";
  if( name=="RgammanhJet"|| name=="Rgammanh" ) axisName = "R_{#gamma+nh}";
  if( name=="EchJet"|| name=="Ech" ) axisName = "E_{ch} [GeV]";
  if( name=="EnhJet"|| name=="Enh" ) axisName = "E_{nh} [GeV]";
  if( name=="EgammaJet"|| name=="Egamma" ) axisName = "E_{#gamma} [GeV]";
  if( name=="EgammanhJet" || name=="Egammanh" ) axisName = "Energy of Photons + Neutral Hadrons in Jet";
  if( name=="PTchJet"|| name=="PTch" ) axisName = "p_{T}^{ch} [GeV/c]";
  if( name=="PTnhJet"|| name=="PTnh" ) axisName = "p_{T}^{nh} [GeV/c]";
  if( name=="PTgammaJet"|| name=="PTgamma" ) axisName = "p_{T}^{#gamma} [GeV/c]";
  if( name=="asymmJet" ) axisName = "Asymmetry";
  if( name=="deltaPhiJet"||name=="deltaPhi" ) axisName = "#Delta#Phi";
  if( name=="massJet" ) axisName = "Jet Invariant Mass [GeV/c^{2}]";
  if( name=="MoEJet" ) axisName = "Jet Invariant Mass / Energy";
  if( name=="ptOverMJet" ) axisName = "Jet p_{T} / Invariant Mass";
  if( name=="pOverMJet" ) axisName = "Jet Momentum / Invariant Mass";
  if( name=="diJetMass" ) axisName = "DiJet Invariant Mass [GeV/c ^{2}]";
  if( name=="EphotAveJet" ) axisName = "Average Photon Energy in Jet [GeV]";
  if( name=="NchJet" || name=="Nch" ) axisName = "Number of Charged Hadrons in Jet";
  if( name=="NgammaJet" || name=="Ngamma" ) axisName = "Number of Photons in Jet";
  if( name=="NnhJet" || name=="Nnh" ) axisName = "Number of Neutral Hadrons in Jet";
  if( name=="NgammanhJet" || name=="Ngammanh" ) axisName = "Number of Photons + Neutral Hadrons in Jet";
  if( name=="hcalIsoPhotReco" ) axisName = "H / E (#DeltaR < 0.4)";
  if( name=="ecalIsoPhotReco" ) axisName = "ECAL Isolation (GT,  #DeltaR < 0.4)";
  if( name=="nTrkIsoPhotReco" ) axisName = "Number of Tracks in #DeltaR < 0.35";
  if( name=="ptTrkIsoPhotReco" ) axisName = "#Sigma p_{T}^{tracks in #DeltaR < 0.35} / p_{T}^{#gamma}";
  if( name=="clusterMajPhotReco" ) axisName = "Photon Cluster Major Axis";
  if( name=="clusterMinPhotReco" ) axisName = "Photon Cluster Minor Axis";
  if( name=="ptPhot" ) axisName = "Photon p_{T} [GeV/c]";
  if( name=="phiPhot" ) axisName = "Photon #phi";
  if( name=="etaPhot" ) axisName = "Photon #eta";
  if( name=="pt2ndJet" ) axisName = "Second Jet p_{T}^{raw} [GeV/c]";
  if( name=="ptSecondJetRel" ) axisName = "p_{T}^{2ndJet} / p_{T}^{#gamma}";
  if( name=="response" || name=="response_loose" || name=="response_clusterOK" ) axisName = "p_{T}^{jet} / p_{T}^{#gamma}";

  return axisName;

}


void DrawBase::shrinkPad(double b, double l, double r, double t) {
 gPad->SetBottomMargin(b); 
 gPad->SetLeftMargin(l);
 gPad->SetRightMargin(r);
 gPad->SetTopMargin(t);
}



void DrawBase::drawStack(const std::string& varY, const std::string& varX, const std::string& etaRegion, const std::string& RECO_GEN, bool isData) const {

  TFile* file;
  if( isData ) file = dataFile_;
  else file = mcFile_;

  std::string histoName;

  std::string suffix;

  if( RECO_GEN=="RECO" )
    suffix = "_vs_" + varX;
  else if( RECO_GEN=="GEN" )
    suffix = "Gen_vs_" + varX;

  if( varX=="eta" )
    suffix = suffix + "_stack";

  if( etaRegion!="" )
    suffix = suffix + "_" + etaRegion;
  
  histoName = "Rch" + suffix;
//  std::cout << "-> Trying to get " << histoName << std::endl << std::flush;
  TH1D* h1_Rch= (TH1D*)file->Get(histoName.c_str());
  histoName = "Rgamma" + suffix;
//  std::cout << "-> Trying to get " << histoName << std::endl << std::flush;
  TH1D* h1_Rgamma= (TH1D*)file->Get(histoName.c_str());
  histoName = "Rnh" + suffix;
//  std::cout << "-> Trying to get " << histoName << std::endl << std::flush;
  TH1D* h1_Rnh= (TH1D*)file->Get(histoName.c_str());
  histoName = "Rmu" + suffix;
//  std::cout << "-> Trying to get " << histoName << std::endl << std::flush;
  TH1D* h1_Rmu= (TH1D*)file->Get(histoName.c_str());
  histoName = "Re" + suffix;
//  std::cout << "-> Trying to get " << histoName << std::endl << std::flush;
  TH1D* h1_Re= (TH1D*)file->Get(histoName.c_str());
  histoName = "Rhfhad" + suffix;
//  std::cout << "-> Trying to get " << histoName << std::endl;
  TH1D* h1_Rhfhad= (TH1D*)file->Get(histoName.c_str());
//  std::cout << "-> Trying to get " << histoName << std::endl;
  histoName = "Rhfem" + suffix;
  TH1D* h1_Rhfem= (TH1D*)file->Get(histoName.c_str());


  int nBins_stack=fitTools::getNbins_stack(varX);
  Double_t Lower_stack[nBins_stack];
  fitTools::getBins_stack(nBins_stack, Lower_stack, varX);


  histoName = varY+"ch_vs_" + varX + "_bis";
  TH1D* h1_Rch_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Rch_stack->SetFillColor(kRed);

  histoName = varY+"gamma_vs_" + varX + "_bis";
  TH1D* h1_Rgamma_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Rgamma_stack->SetFillColor(kBlue);

  histoName = varY+"nh_vs_" + varX + "_bis";
  TH1D* h1_Rnh_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Rnh_stack->SetFillColor(kGreen);

  histoName = varY+"mu_vs_" + varX + "_bis";
  TH1D* h1_Rmu_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Rmu_stack->SetFillColor(kOrange);

  histoName = varY+"e_vs_" + varX + "_bis";
  TH1D* h1_Re_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Re_stack->SetFillColor(kCyan);

  histoName = varY+"hfhad_vs_" + varX + "_bis";
  TH1D* h1_Rhfhad_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Rhfhad_stack->SetFillColor(kMagenta-9);

  histoName = varY+"hfhem_vs_" + varX + "_bis";
  TH1D* h1_Rhfem_stack = new TH1D(histoName.c_str(), "", nBins_stack-1, Lower_stack);
  h1_Rhfem_stack->SetFillColor(kBlue-5);


  for(int i=1; i<(nBins_stack); ++i) {
    h1_Rch_stack->SetBinContent( i, h1_Rch->GetBinContent(i) );
    h1_Rgamma_stack->SetBinContent( i, h1_Rgamma->GetBinContent(i) );
    h1_Rnh_stack->SetBinContent( i, h1_Rnh->GetBinContent(i) );
    h1_Rmu_stack->SetBinContent( i, h1_Rmu->GetBinContent(i) );
    h1_Re_stack->SetBinContent( i, h1_Re->GetBinContent(i) );
    h1_Rhfhad_stack->SetBinContent( i, h1_Rhfhad->GetBinContent(i) );
    h1_Rhfem_stack->SetBinContent( i, h1_Rhfem->GetBinContent(i) );
  }


  THStack* stack = new THStack("stack", "");
  stack->Add(h1_Rch_stack);
  stack->Add(h1_Rgamma_stack);
  stack->Add(h1_Re_stack);
  stack->Add(h1_Rnh_stack);
  stack->Add(h1_Rhfhad_stack);
  stack->Add(h1_Rhfem_stack);
  stack->Add(h1_Rmu_stack);
  stack->Draw(); //this draw is magically needed

//THStack* stack = new THStack("stack", "");
//stack->Add(h1_Rch_vs_eta, "HISTO");
//stack->Add(h1_Rgamma_vs_eta, "HISTO");
//stack->Add(h1_Re_vs_eta, "HISTO");
//stack->Add(h1_Rnh_vs_eta, "HISTO");
//if( RECO_GEN=="RECO" ) {
//  stack->Add(h1_Rhfhad_vs_eta, "HISTO");
//  stack->Add(h1_Rhfem_vs_eta, "HISTO");
//}
//stack->Add(h1_Rmu_vs_eta, "HISTO");
//stack->Draw(); //this draw is magically needed

  std::string varX_name;
  if( varX=="eta" ) {
    varX_name = "#eta";
  } else if( varX=="pt" ) {
    if( RECO_GEN=="RECO" ) {
      varX_name = "Raw p_{T}";
    } else {
      varX_name = "p_{T}";
    }
  } else if( varX=="ptCorr" ) {
    varX_name = "Corrected p_{T}";
  } else if( varX=="phi" ) {
    varX_name = "#phi";
  }

  std::string xTitle = (RECO_GEN=="RECO") ? "PFJet " : "GenJet ";
  xTitle = xTitle + varX_name;

  std::string yTitle;
  if( varY == "R" )
    yTitle = "Mean Fraction of Jet Energy";
  else if( varY == "N" )
    yTitle = "Particle Multiplicity";
  else if( varY == "Pt" ) {
    if( RECO_GEN=="RECO" )
      yTitle = "<p_{T}^{RECO}> [GeV/c]";
    else if( RECO_GEN=="GEN" )
      yTitle = "<p_{T}^{GEN}> [GeV/c]";
  }
  stack->GetXaxis()->SetTitle(xTitle.c_str());
  stack->GetYaxis()->SetTitle(yTitle.c_str());

  TPaveText* label_cms = new TPaveText(0.25, 0.83, 0.42, 0.87, "brNDC");
  label_cms->SetFillColor(kWhite);
  label_cms->SetTextSize(0.038);
  label_cms->SetTextFont(62);
  label_cms->AddText("CMS Preliminary 2010");


  TPaveText* label_sqrt = new TPaveText(0.25, 0.78, 0.42, 0.82, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  std::string label_sqrt_text;
  if( isData )
    label_sqrt_text = "#sqrt{s} = 7 TeV, DATA";
  else
    label_sqrt_text = "#sqrt{s} = 7 TeV, MC";
  label_sqrt->AddText(label_sqrt_text.c_str());

  Float_t yMin_cuts = ( varX!="eta" && varX!="pt" && varX!="ptCorr" ) ? 0.77 : 0.78;
  Float_t yMax_cuts = ( varX!="eta" && varX!="pt" && varX!="ptCorr" ) ? 0.89 : 0.88;
  Float_t textSize_cuts = ( varX!="eta" && varX!="pt" && varX!="ptCorr" ) ? 0.03 : 0.035;

  TPaveText* label_cuts = new TPaveText(0.6, yMin_cuts, 0.8, yMax_cuts,  "brNDC");
  label_cuts->SetFillColor(kWhite);
  label_cuts->SetTextSize(textSize_cuts);
  label_cuts->SetTextFont(42);
  label_cuts->AddText("Anti-k_{T} R=0.5");
  std::string apexText = ( RECO_GEN=="RECO" ) ? raw_corr_ : "GEN";
  //std::string apexText = raw_corr_;
  if( varX != "eta" ) {
    std::string etaRange = getEtaRangeText( etaRegion );
    label_cuts->AddText(etaRange.c_str());
  }
  if( (varX != "pt") && (varX != "ptCorr") ) {
    char label_pt_text[100];
    sprintf( label_pt_text, "p_{T}^{%s} > %d GeV/c", apexText.c_str(), pt_thresh_);
    label_cuts->AddText(label_pt_text);
  }

  Float_t xMin, xMax;
  if( varX=="eta" ) {
    xMin=-5.5;
    xMax= 5.5;
  } else {
    xMin = Lower_stack[0];
    xMax = Lower_stack[nBins_stack-1];
  }


  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 1.2);
  h2_axes->SetXTitle(xTitle.c_str());
  h2_axes->SetYTitle(yTitle.c_str());
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);
  if( varX=="pt" || varX=="ptCorr" ) {
    h2_axes->GetXaxis()->SetMoreLogLabels();
    h2_axes->GetXaxis()->SetNoExponent();
  }

  if( varX != "eta" )
    gStyle->SetPadTickY(1);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->SetLeftMargin(1.2);
  c1->SetBottomMargin(1.2);
  c1->cd();
  if( varX=="pt" || varX=="ptCorr" )
    c1->SetLogx();
  h2_axes->Draw();
  stack->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_cuts->Draw("same");
  gPad->RedrawAxis();


  std::string canvasName;
  if( RECO_GEN=="RECO" ) {
    if( isData ) {
      canvasName = outputdir_ + "/stack_" + varY +"_vs_" + varX;
    } else {
      canvasName = outputdir_ + "/stackMC_" + varY +"_vs_" + varX;
    }
  } else if( RECO_GEN=="GEN" ) {
    canvasName = outputdir_ +"/stackGEN_" + varY +"_vs_" + varX;
  }
  std::string canvasName_eps = canvasName + ".eps";
  std::string canvasName_png = canvasName + ".png";
  std::string canvasName_pdf = canvasName + ".pdf";
  c1->SaveAs(canvasName_eps.c_str());
  c1->SaveAs(canvasName_png.c_str());
  if( pdf_aussi_ )
    c1->SaveAs(canvasName_pdf.c_str());

  delete c1;
  delete h2_axes;

  delete h1_Rch_stack;
  delete h1_Rgamma_stack;
  delete h1_Rnh_stack;
  delete h1_Re_stack;
  delete h1_Rmu_stack;
  delete h1_Rhfhad_stack;
  delete h1_Rhfem_stack;

  gStyle->SetPadTickY(0);


} //drawStack


void DrawBase::set_dataFile( TFile* dataFile ) {

  dataFile_ = dataFile;

  TH1F* h1_lumi = (TH1F*)dataFile_->Get("totalLumi");
  if( h1_lumi==0 ) {
    std::cout << "WARNING! Lumi histogram not found!" << std::endl;
    lumi_ = 0;
  } else {
    lumi_ = h1_lumi->GetBinContent(1);
  }

}


std::string DrawBase::getEtaRangeText( const std::string& etaRegion ) const {

  char etaRange_ch[100];
  sprintf(etaRange_ch, "|#eta| < %.1f", etamax_);
  std::string etaRange(etaRange_ch);
  if( etaRegion=="barrel") etaRange = "|#eta| < 1.4";
  if( etaRegion=="endcap") etaRange = "1.4 < |#eta| < 2.4";
  if( etaRegion=="eta02") etaRange = "|#eta| < 2";
  if( etaRegion=="eta23") etaRange = "2 < |#eta| < 3";
  if( etaRegion=="eta163") etaRange = "1.6 < |#eta| < 3";
  if( etaRegion=="eta1425") etaRange = "1.4 < |#eta| < 2.5";
  if( etaRegion=="eta1430") etaRange = "1.4 < |#eta| < 3.0";
  if( etaRegion=="Rch050") etaRange = "|#eta| < 2.5, R_{ch} < 50%";
  if( etaRegion=="Rch5070") etaRange = "|#eta| < 2.5, 50 < R_{ch} < 70%";
  if( etaRegion=="Rch70100") etaRange = "|#eta| < 2.5, R_{ch} > 70%";
  if( etaRegion=="ptPhot_10_15") etaRange = "10 < p_{T}^{#gamma} < 15 GeV/c";
  if( etaRegion=="ptPhot_15_3500") etaRange = "p_{T}^{#gamma} > 15 GeV/c";

  return etaRange;

}


std::string DrawBase::getSqrtText() const {

  if( lumi_==0. ) {
    std::cout << "WARNING! Lumi histogram not found!" << std::endl;
    return std::string("#sqrt{s} = 7 TeV");
  }
  float lumi4Text(lumi_);
  bool onlyOneDecimal=false;
  std::string units = "#mub ^{-1}";
  if( lumi4Text > 10. ) {
    lumi4Text /= 1000.;
    units = "nb ^{-1}";
  }
  if( lumi4Text > 50. ) {
    lumi4Text /= 1000.;
    units = "pb ^{-1}";
  } else if(  lumi4Text > 10. ) {
    onlyOneDecimal=true;
  }

  char label_sqrt_text[150];
  if( onlyOneDecimal )
    sprintf( label_sqrt_text, "#sqrt{s} = 7 TeV, L = %.0f %s", lumi4Text, units.c_str());
  else
    sprintf( label_sqrt_text, "#sqrt{s} = 7 TeV, L = %.1f %s", lumi4Text, units.c_str());

  std::string returnString(label_sqrt_text);

  return returnString;

}

std::string DrawBase::getAlgoName() const {

  std::string algoName;

  if( jetAlgo_=="akt5" ) {
    algoName = "Anti-k_{T} 0.5 ";
  } else {
    std::cout << "Jet algo '" << jetAlgo_ << "' currently not supported. Exiting." << std::endl;
    exit( 918 );
  }

  if( recoType_=="pf" ) {
    algoName += "PFJets";
  } else if( recoType_=="calo" ) {
    algoName += "CaloJets"; 
  } else if( recoType_=="jpt" ) {
    algoName += "JPTJets";
  } else {
    std::cout << "Reco Type '" << recoType_ << "' currently not supported. Exiting." << std::endl;
    exit( 919 );
  }

  return algoName;

}


TGraphErrors* DrawBase::getGraphRatio( TGraphErrors* gr_data, TGraphErrors* gr_MC ) {

  TGraphErrors* gr_ratio = new TGraphErrors(0);

  for( unsigned i= 0; i<gr_MC->GetN(); ++i ) {

    Double_t datax, datay;
    gr_data->GetPoint( i, datax, datay );
    Double_t dataxerr = gr_data->GetErrorX(i);
    Double_t datayerr = gr_data->GetErrorY(i);

    Double_t mcx, mcy;
    gr_MC->GetPoint( i, mcx, mcy );
    Double_t mcxerr = gr_MC->GetErrorX(i);
    Double_t mcyerr = gr_MC->GetErrorY(i);

    Double_t ratiox = mcx;
    Double_t ratioxerr = mcxerr;

    Double_t ratioy = datay / mcy;
    Double_t ratioyerr = sqrt( datayerr*datayerr/(mcy*mcy) + datay*datay*mcyerr*mcyerr/(mcy*mcy*mcy*mcy) );


    gr_ratio->SetPoint( i, ratiox, ratioy );
    gr_ratio->SetPointError( i, ratioxerr, ratioyerr );

  } //for points

  return gr_ratio;

}



