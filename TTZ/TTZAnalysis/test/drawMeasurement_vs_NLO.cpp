#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "DrawBase.h"



int main() {

  // set the style:
  DrawBase* db = new DrawBase("meas_vs_NLO");
  db->set_lumi(4980.);
  db->set_lumiOnRightSide();
  TFile* f = new TFile();
  db->add_dataFile(f, "p");



  TStyle* style_ = new TStyle("DrawBaseStyle", "");
  style_->SetCanvasColor(0);
  style_->SetPadColor(0);
  style_->SetFrameFillColor(0);
  style_->SetStatColor(0);
  style_->SetOptStat(0);
  style_->SetTitleFillColor(0);
  style_->SetCanvasBorderMode(0);
  style_->SetPadBorderMode(0);
  style_->SetFrameBorderMode(0);
  style_->cd();

  // For the canvas:
  style_->SetCanvasBorderMode(0);
  style_->SetCanvasColor(kWhite);
  style_->SetCanvasDefH(600); //Height of canvas
  style_->SetCanvasDefW(600); //Width of canvas
  style_->SetCanvasDefX(0);   //POsition on screen
  style_->SetCanvasDefY(0);

  // For the Pad:
  style_->SetPadBorderMode(0);
  // style_->SetPadBorderSize(Width_t size = 1);
  style_->SetPadColor(kWhite);
  style_->SetPadGridX(false);
  style_->SetPadGridY(false);
  style_->SetGridColor(0);
  style_->SetGridStyle(3);
  style_->SetGridWidth(1);

  // For the frame:
  style_->SetFrameBorderMode(0);
  style_->SetFrameBorderSize(1);
  style_->SetFrameFillColor(0);
  style_->SetFrameFillStyle(0);
  style_->SetFrameLineColor(1);
  style_->SetFrameLineStyle(1);
  style_->SetFrameLineWidth(1);

//// For the histo:
//  // style_->SetHistFillColor(1);
//  // style_->SetHistFillStyle(0);
//  style_->SetHistLineColor(1);
//  style_->SetHistLineStyle(0);
//  style_->SetHistLineWidth(1);
//  // style_->SetLegoInnerR(Float_t rad = 0.5);
//  // style_->SetNumberContours(Int_t number = 20);

//  style_->SetEndErrorSize(2);
////  style_->SetErrorMarker(20);
//  style_->SetErrorX(0.);
//  
//  style_->SetMarkerStyle(20);

////For the fit/function:
//  style_->SetOptFit(1);
//  style_->SetFitFormat("5.4g");
//  style_->SetFuncColor(2);
//  style_->SetFuncStyle(1);
//  style_->SetFuncWidth(1);

////For the date:
//  style_->SetOptDate(0);
//  // style_->SetDateX(Float_t x = 0.01);
//  // style_->SetDateY(Float_t y = 0.01);

//// For the statistics box:
//  style_->SetOptFile(0);
//  style_->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
//  style_->SetStatColor(kWhite);
//  style_->SetStatFont(42);
//  style_->SetStatFontSize(0.025);
//  style_->SetStatTextColor(1);
//  style_->SetStatFormat("6.4g");
//  style_->SetStatBorderSize(1);
//  style_->SetStatH(0.1);
//  style_->SetStatW(0.15);
//  // style_->SetStatStyle(Style_t style = 1001);
//  // style_->SetStatX(Float_t x = 0);
//  // style_->SetStatY(Float_t y = 0);

  // Margins:
  style_->SetPadTopMargin(0.05);
  style_->SetPadBottomMargin(0.15);//0.13);
//style_->SetPadLeftMargin(0.);//0.16);
//style_->SetPadRightMargin(0.0);//0.02);

  // For the Global title:

  style_->SetOptTitle(0);
  style_->SetTitleFont(42);
  style_->SetTitleColor(1);
  style_->SetTitleTextColor(1);
  style_->SetTitleFillColor(10);
  style_->SetTitleFontSize(0.05);
  // style_->SetTitleH(0); // Set the height of the title box
  // style_->SetTitleW(0); // Set the width of the title box
  // style_->SetTitleX(0); // Set the position of the title box
  // style_->SetTitleY(0.985); // Set the position of the title box
  // style_->SetTitleStyle(Style_t style = 1001);
  // style_->SetTitleBorderSize(2);

  // For the axis titles:

  style_->SetTitleColor(1, "XYZ");
  style_->SetTitleFont(42, "XYZ");
  style_->SetTitleSize(0.05, "XYZ");
  // style_->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // style_->SetTitleYSize(Float_t size = 0.02);
  style_->SetTitleXOffset(1.15);//0.9);
  style_->SetTitleYOffset(1.4); // => 1.15 if exponents
  // style_->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  style_->SetLabelColor(1, "XYZ");
  style_->SetLabelFont(42, "XYZ");
  style_->SetLabelOffset(0.007, "XYZ");
  style_->SetLabelSize(0.045, "XYZ");

  // For the axis:

  style_->SetAxisColor(1, "XYZ");
  style_->SetStripDecimals(kTRUE);
  style_->SetTickLength(0.03, "XYZ");
  style_->SetNdivisions(510, "XYZ");
  style_->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  style_->SetPadTickY(1);



  TGraphAsymmErrors* gr_trilept = new TGraphAsymmErrors();
  TGraphAsymmErrors* gr_ssdl = new TGraphAsymmErrors();
  TGraphAsymmErrors* gr_combined = new TGraphAsymmErrors();

  gr_trilept->SetPoint(0, 0.565937, 2.5);
  gr_trilept->SetPointError(0, 0.285454, 0.390798, 0., 0.);
  gr_trilept->SetMarkerStyle(22);
  gr_trilept->SetMarkerSize(2.);
  gr_trilept->SetMarkerColor(38);

  gr_ssdl->SetPoint(0, 0.270711, 1.5);
  gr_ssdl->SetPointError(0, 0.189347, 0.2302, 0., 0.);
  gr_ssdl->SetMarkerStyle(21);
  gr_ssdl->SetMarkerSize(2.);
  gr_ssdl->SetMarkerColor(39);

  gr_combined->SetPoint(0, 0.386399224, 0.5);
  gr_combined->SetPointError(0, 0.180300611, 0.180300611, 0., 0.);
  gr_combined->SetMarkerStyle(20);
  gr_combined->SetMarkerSize(2.);
  gr_combined->SetMarkerColor(1);


  TLine* line_nlo = new TLine(0.308, 0., 0.308, 3.);
  line_nlo->SetLineWidth(2);
  line_nlo->SetLineColor(46);

  TLine* line_nlo_plus = new TLine(0.308+0.062, 0., 0.308+0.062, 3.);
  line_nlo_plus->SetLineStyle(2);
  line_nlo_plus->SetLineWidth(2);
  line_nlo_plus->SetLineColor(14);

  TLine* line_nlo_minus = new TLine(0.308-0.062, 0., 0.308-0.062, 3.);
  line_nlo_minus->SetLineStyle(2);
  line_nlo_minus->SetLineWidth(2);
  line_nlo_minus->SetLineColor(14);

  TLegend* legend = new TLegend(0.45, 0.2, 0.82, 0.5);
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( gr_trilept, "Trilepton Channel", "P" );
  legend->AddEntry( gr_ssdl, "Dilepton Channel", "P" );
  legend->AddEntry( gr_combined, "Combination", "P" );
  legend->AddEntry( line_nlo, "NLO Calculation", "L" );

  //TPaveText* label_sqrt = db->get_labelSqrt();
  TPaveText* label_sqrt = new TPaveText( 0.32, 0.953, 0.93, 0.975, "brNDC");
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetFillColor(0);
  label_sqrt->SetTextFont(62);
  label_sqrt->SetTextAlign(31);
  label_sqrt->AddText("CMS Preliminary, L = 4.98 fb ^{-1} at  #sqrt{s} = 7 TeV");
  
  TH2D* h2_axes = new TH2D("axes", "", 10, 0., 1.7, 1, 0., 3.);
  h2_axes->GetYaxis()->SetBinLabel(1, "");
  h2_axes->SetXTitle("Cross Section [pb]");


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);

  h2_axes->Draw();
  line_nlo->Draw("same");
  line_nlo_plus->Draw("same");
  line_nlo_minus->Draw("same");
  gr_trilept->Draw("P same");
  gr_ssdl->Draw("P same");
  gr_combined->Draw("P same");
  legend->Draw("same");
  label_sqrt->Draw("same");

  c1->SaveAs("measurement_vs_NLO.eps");

  return 0;

}
