#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"



void drawSingleFit( DrawBase* db, const std::string& varName, int ptMin, int ptMax, int rhoBin );


int main() {

  std::string fitFileName = "Fit_newBin.root";
  TFile* fitFile = TFile::Open(fitFileName.c_str());

  DrawBase* db = new DrawBase("QGFits");
  //db->set_lumiOnRightSide();

  db->add_mcFile( fitFile, "QCD_Flat_Summer11_newBins", "Gluons", 38 );
  db->add_mcFile( fitFile, "QCD_Flat_Summer11_newBins", "Quarks", 46 );

  db->set_shapeNormalization();
  db->set_outputdir();
  

  const int nPtBins = 18;
  Double_t ptBins[nPtBins+1];
  fitTools::getBins_int( nPtBins, ptBins, 20., 1000. );
  ptBins[nPtBins] = 3500.;

  for( unsigned ipt=0; ipt<nPtBins; ++ipt ) {
    for( unsigned irho=0; irho<20; ++irho ) {

      drawSingleFit( db, "nCharged", ptBins[ipt], ptBins[ipt+1], irho );

    }
  }

  return 0;

} 


void drawSingleFit( DrawBase* db, const std::string& varName, int ptMin, int ptMax, int rhoBin ) {

  TFile* fitFile = db->get_mcFile(0).file;

  char canvasName[400];
  sprintf( canvasName, "%s_pt_%d_%d_rho%d", varName.c_str(), ptMin, ptMax, rhoBin );

  char gluonHistoName[400];
  sprintf( gluonHistoName, "%s_gluon_pt%d_%d_rho%d", varName.c_str(), ptMin, ptMax, rhoBin );
  char quarkHistoName[400];
  sprintf( quarkHistoName, "%s_quark_pt%d_%d_rho%d", varName.c_str(), ptMin, ptMax, rhoBin );

  TCanvas* thisCanvas = (TCanvas*)fitFile->Get(canvasName);
  thisCanvas->cd();


  TList* list = thisCanvas->GetListOfPrimitives();

  TH1D* h1_gluon = (TH1D*)list->FindObject(gluonHistoName);
  TH1D* h1_quark = (TH1D*)list->FindObject(quarkHistoName);

  int quarkColor = 46;
  int gluonColor = 38;

  h1_gluon->SetMarkerStyle(21);
  h1_gluon->SetMarkerColor(gluonColor);
  h1_gluon->SetLineColor(gluonColor);
  h1_gluon->SetLineWidth(2);
  h1_gluon->SetFillColor(gluonColor);
  h1_gluon->SetFillStyle(3004);

  h1_quark->SetMarkerStyle(20);
  h1_quark->SetLineColor(quarkColor);
  h1_quark->SetMarkerColor(quarkColor);
  h1_quark->SetLineWidth(2);
  h1_quark->SetFillColor(quarkColor);
  h1_quark->SetFillStyle(3005);


  std::string functionName = (varName=="ptD") ? "functionPtD" : "gamma";
  std::string functionName_gluon = functionName + "_gluon";
  std::string functionName_quark = functionName + "_quark";

  TF1* f1_gluon = (TF1*)list->FindObject(functionName_gluon.c_str());
  TF1* f1_quark = (TF1*)list->FindObject(functionName_quark.c_str());

  f1_gluon->SetLineColor(gluonColor);
  f1_quark->SetLineColor(quarkColor);
  

  std::string axisName;
  if( varName=="nCharged" )
    axisName = "Jet Charged Multiplicity";
  else if( varName=="nNeutral" )
    axisName = "Jet Neutral Multiplicity";
  else if( varName=="ptD" )
    axisName = "Jet p_{T}D";


  float xMin = h1_quark->GetXaxis()->GetXmin();
  float xMax = (varName=="ptD") ? h1_quark->GetXaxis()->GetXmax() : 50.;
  float yMax_gluon = h1_gluon->GetMaximum();
  float yMax_quark = h1_quark->GetMaximum();
  float yMax = (yMax_gluon>yMax_quark) ? yMax_gluon : yMax_quark;
  yMax *= 1.3;

  std::string yTitle = "Normalized to Unity";

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., yMax);
  h2_axes->SetXTitle(axisName.c_str());
  h2_axes->SetYTitle(yTitle.c_str());

  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();
  
  int rhoMin = rhoBin;
  int rhoMax = rhoBin+1;

  char ptText[400];
  sprintf( ptText, "%d < p_{T} < %d GeV", ptMin, ptMax );
  char rhoText[400];
  sprintf( rhoText, "%d < #rho_{PF} < %d GeV", rhoMin, rhoMax );

  TPaveText* pt_rho_label = new TPaveText( 0.25, 0.79, 0.45, 0.91, "brNDC" );
  pt_rho_label->SetTextSize(0.035);
  pt_rho_label->SetFillColor(0);
  pt_rho_label->AddText( ptText );
  pt_rho_label->AddText( rhoText );

  char legendTitle[500];
  sprintf( legendTitle, "%d < p_{T} < %d, %d < #rho < %d", ptMin, ptMax, rhoMin, rhoMax );

  //TLegend* legend = new TLegend( 0.55, 0.65, 0.9, 0.9, legendTitle );
  TLegend* legend = new TLegend( 0.65, 0.75, 0.9, 0.9 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_gluon, "Gluons", "PL" );
  legend->AddEntry( h1_quark, "Quarks", "PL" );

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();
  h2_axes->Draw();
  legend->Draw("same");
  h1_quark->Draw("SAME");
  h1_gluon->Draw("SAME");
  f1_quark->Draw("SAME");
  f1_gluon->Draw("SAME");
  label_sqrt->Draw("same");
  label_cms->Draw("same");
  pt_rho_label->Draw("same");

  gPad->RedrawAxis();
  
  char saveCanvasName_char[900];
  sprintf( saveCanvasName_char, "%s/%s_pt%d_%d_rho%d", db->get_outputdir().c_str(), varName.c_str(), ptMin, ptMax, rhoBin );
  std::string saveCanvasName(saveCanvasName_char);
  std::string saveCanvasName_eps = saveCanvasName + ".eps";

  c1->SaveAs(saveCanvasName_eps.c_str());



  TCanvas* c1_log = new TCanvas("c1_log", "", 600, 600);
  c1_log->cd();
  c1_log->SetLogy();
 
  float yMin_log = 0.000001;  
  float yMax_log = 50.*yMax;

  TH2D* h2_axes_log = new TH2D("axes_log", "", 10, xMin, xMax, 10, yMin_log, yMax_log);
  h2_axes_log->SetXTitle(axisName.c_str());
  h2_axes_log->SetYTitle(yTitle.c_str());

  h2_axes_log->Draw();
  legend->Draw("same");
  h1_quark->Draw("SAME");
  h1_gluon->Draw("SAME");
  f1_quark->Draw("SAME");
  f1_gluon->Draw("SAME");
  label_sqrt->Draw("same");
  label_cms->Draw("same");
  pt_rho_label->Draw("same");

  std::string saveCanvasName_log = saveCanvasName + "_log";
  std::string saveCanvasName_log_eps = saveCanvasName_log + ".eps";
  c1_log->SaveAs(saveCanvasName_log_eps.c_str());

  gPad->RedrawAxis();

  delete h2_axes;
  delete h2_axes_log;
  delete legend;
  delete c1;
  delete c1_log;

}
