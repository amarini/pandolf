#include "DrawBase.h"



void draw_vs_pt( DrawBase* db, const std::string& varname, float etaMin, float etaMax, TFile* file_PF, TFile* file_calo );


int main() {

  DrawBase* db = new DrawBase("PF_vs_Calo");

  TFile* file_PF = TFile::Open("PhotonJetGraphs_DATA_Run2011A_extend1_vs_G_Summer11_pfakt5_LUMI_2ndJet10.root");
  TFile* file_calo = TFile::Open("PhotonJetGraphs_DATA_Run2011A_extend1_vs_G_Summer11_akt5_LUMI_2ndJet10.root");

  draw_vs_pt(db, "response", 0., 1.3, file_PF, file_calo);
  draw_vs_pt(db, "resolution", 0., 1.3, file_PF, file_calo);

  return 0;

}



void draw_vs_pt( DrawBase* db, const std::string& varname, float etaMin, float etaMax, TFile* file_PF, TFile* file_calo ) {

  TString varname_tstr(varname);

  char graphName[200];
  sprintf( graphName, "MC%s_eta%.0f%.0f_vs_pt", varname.c_str(), 10.*etaMin, 10.*etaMax );

  TGraphErrors* gr_PF   = (TGraphErrors*)file_PF->Get(graphName);
  TGraphErrors* gr_calo = (TGraphErrors*)file_calo->Get(graphName);

  gr_PF->SetMarkerSize(1.6);
  gr_calo->SetMarkerSize(1.6);

  gr_PF->SetMarkerColor(46);
  gr_calo->SetMarkerColor(38);

  gr_PF->SetMarkerStyle(20);
  gr_calo->SetMarkerStyle(21);


  double xMin = 20.;
  double xMax = 600.;
  double yMin = (varname_tstr.Contains("response")) ? 0.3 : 0.;
  double yMax = (varname_tstr.Contains("response")) ? 1.1 : 0.7;

  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, yMin, yMax);
  h2_axes->GetXaxis()->SetMoreLogLabels();
  h2_axes->GetXaxis()->SetNoExponent();
  h2_axes->SetXTitle("p_{T} [GeV]");
  if( varname_tstr.Contains("response") )
    h2_axes->SetYTitle("Response");
  else if( varname_tstr.Contains("resolution"))
    h2_axes->SetYTitle("Resolution [%]");

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  c1->SetLogx();

  TPaveText* label_cms = db->get_labelCMS();
  TPaveText* label_sqrt = db->get_labelSqrt();

  char legendTitle[200];
  if( etaMin!=0. )
    sprintf( legendTitle, "%.1f < |#eta| < %.1f", etaMin, etaMax );
  else
    sprintf( legendTitle, "|#eta| < %.1f", etaMax );

  double legend_xMin = 0.46;
  double legend_xMax = 0.91;
  double legend_yMin = (varname_tstr.Contains("response")) ? 0.15 : 0.55;
  double legend_yMax = (varname_tstr.Contains("response")) ? 0.5 : 0.9;


  TLegend* legend = new TLegend( legend_xMin, legend_yMin, legend_xMax, legend_yMax, legendTitle );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  legend->AddEntry( gr_calo, "Calorimeter Jets", "P" );
  legend->AddEntry( gr_PF, "Particle Flow Jets", "P" );

  h2_axes->Draw();
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  legend->Draw("same");
  
  gr_calo->Draw("P same");
  gr_PF->Draw("P same");

  char canvasName[300];
  sprintf( canvasName, "PF_vs_Calo/%s.eps", graphName );
  c1->SaveAs(canvasName);
  
}
