#include <cstdlib>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "DrawBase.h"



std::vector<float> drawSinglePlot( DrawBase* db, TFile* file, const std::string& varName, const std::string& units, const std::string& yAxisTitle, const std::string& rhoType, const std::string& EB_EE );
void drawSinglePlot_Mean( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisName, const std::string& units, const std::string& rhoType, const std::string& EB_EE );
void drawCompare( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisName, const std::string& rhoType, const std::string& EB_EE );
void drawCompare_fullSelection( DrawBase* db, TFile* file, const std::string& rhoType, const std::string& EB_EE, bool vs_nvertex=false );


int main( int argc, char* argv[] )  {


  if( argc != 2 ) {
    std::cout << "USAGE: ./drawRhoStudies [dataset]" << std::endl;
    exit(13);
  }


  std::string dataset(argv[1]);

  std::string fileName = "RhoStudies_"+dataset+".root";
  TFile* file = TFile::Open(fileName.c_str());

  DrawBase* db = new DrawBase("RhoStudies");
  std::string outputdirName = "RhoStudiesPlots_"+dataset;
  db->set_outputdir(outputdirName);

//drawSinglePlot( db, file, "hcalIsoPID", "HCAL Isolation", "Calo", "EB");
//drawSinglePlot( db, file, "ecalIsoPID", "ECAL Isolation", "Calo", "EB");

  char ofsName[400];
  sprintf( ofsName, "parameters_%s.txt", dataset.c_str());
  ofstream ofs(ofsName);

  ofs <<  "if(TMath::Abs(etaPhot[i]) < 1.44) {" << std::endl;

  std::vector<float> parameters;
  parameters = drawSinglePlot( db, file, "twrHcalIso", "HCAL Isolation", "GeV", "PF", "EB");
  ofs << "  hcaliso = (pid_twrHCAL[i] < ptPhot[i] * pid.hcaliso_rel + " << parameters[0] << " + " << parameters[1] << "*rhoPF - 2.0 + pid.hcaliso_abs );" << std::endl;
  parameters = drawSinglePlot( db, file, "HoverE", "H/E", "", "PF", "EB");
  ofs << "  hovereiso = (pid_HoverE[i] < " << parameters[0] << " + " << parameters[1] << "*rhoPF - 0.02 + pid.hovereiso );" << std::endl;
  parameters = drawSinglePlot( db, file, "jurEcalIso", "ECAL Isolation", "GeV", "PF", "EB");
  ofs << "  ecaliso = (pid_jurECAL[i] < ptPhot[i] * pid.ecaliso_rel + " << parameters[0] << " + " << parameters[1] << "*rhoPF - 2.0 + pid.ecaliso_abs );" << std::endl;
  parameters = drawSinglePlot( db, file, "jurEcalIso", "ECAL Isolation", "GeV", "PF", "EB1");
  parameters = drawSinglePlot( db, file, "jurEcalIso", "ECAL Isolation", "GeV", "PF", "EB2");
  parameters = drawSinglePlot( db, file, "hlwTrackIso", "Track Isolation", "GeV", "PF", "EB");
  ofs << "  ptiso = (pid_hlwTrackNoDz[i] < ptPhot[i] * pid.trackiso_rel + " << parameters[0] << " + " << parameters[1] << "*rhoPF - 1.5 + pid.trackiso_abs);" << std::endl;
  parameters = drawSinglePlot( db, file, "etawid", "#sigma i#eta i#eta", "", "PF", "EB");

  ofs << " }else{" << std::endl;
  parameters = drawSinglePlot( db, file, "twrHcalIso", "HCAL Isolation", "GeV", "PF", "EE");
  ofs << "  hcaliso = (pid_twrHCAL[i] < ptPhot[i] * pid.hcaliso_rel + " << parameters[0] << " + " << parameters[1] << "*rhoPF - 2.0 + pid.hcaliso_abs );" << std::endl;
  parameters = drawSinglePlot( db, file, "HoverE", "H/E", "", "PF", "EE");
  ofs << "  hovereiso = (pid_HoverE[i] < " << parameters[0] << " + " << parameters[1] << "*rhoPF - 0.02 + pid.hovereiso );" << std::endl;
  parameters = drawSinglePlot( db, file, "jurEcalIso", "ECAL Isolation", "GeV", "PF", "EE");
  ofs << "  ecaliso = (pid_jurECAL[i] < ptPhot[i] * pid.ecaliso_rel + " << parameters[0] << " + " << parameters[1] << "*rhoPF - 2.0 + pid.ecaliso_abs );" << std::endl;
  parameters = drawSinglePlot( db, file, "hlwTrackIso", "Track Isolation", "GeV", "PF", "EE");
  ofs << "  ptiso = (pid_hlwTrackNoDz[i] < ptPhot[i] * pid.trackiso_rel + " << parameters[0] << " + " << parameters[1] << "*rhoPF - 1.5 + pid.trackiso_abs);" << std::endl;
  parameters = drawSinglePlot( db, file, "etawid", "#sigma i#eta i#eta", "", "PF", "EE");

  ofs << "}" << std::endl;
  ofs.close();

  drawSinglePlot_Mean( db, file, "twrHcalIso", "HCAL Isolation", "GeV", "PF", "EB");
//  drawSinglePlot_Mean( db, file, "HoverE", "H/E", "", "PF", "EB");
  drawSinglePlot_Mean( db, file, "jurEcalIso", "ECAL Isolation", "GeV", "PF", "EB");
  drawSinglePlot_Mean( db, file, "hlwTrackIso", "Track Isolation", "GeV", "PF", "EB");

  drawSinglePlot_Mean( db, file, "twrHcalIso", "HCAL Isolation", "GeV", "PF", "EE");
//  drawSinglePlot_Mean( db, file, "HoverE", "H/E", "", "PF", "EE");
  drawSinglePlot_Mean( db, file, "jurEcalIso", "ECAL Isolation", "GeV", "PF", "EE");
  drawSinglePlot_Mean( db, file, "hlwTrackIso", "Track Isolation", "GeV", "PF", "EE");

  drawCompare( db, file, "twrHcalIso", "HCAL Isolation", "PF", "EB");
  drawCompare( db, file, "jurEcalIso", "ECAL Isolation", "PF", "EB");
  drawCompare( db, file, "hlwTrackIso", "Track Isolation", "PF", "EB");
  drawCompare( db, file, "HoverE", "H/E", "PF", "EB");

  drawCompare( db, file, "twrHcalIso", "HCAL Isolation", "PF", "EE");
  drawCompare( db, file, "jurEcalIso", "ECAL Isolation", "PF", "EE");
  drawCompare( db, file, "hlwTrackIso", "Track Isolation", "PF", "EE");
  drawCompare( db, file, "HoverE", "H/E", "PF", "EE");

  drawCompare_fullSelection( db, file, "PF", "EB" );
  drawCompare_fullSelection( db, file, "PF", "EE" );

  drawCompare_fullSelection( db, file, "PF", "EB", (bool)true );
  drawCompare_fullSelection( db, file, "PF", "EE", (bool)true );

  return 0;

}



std::vector<float> drawSinglePlot( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisName, const std::string& units, const std::string& rhoType, const std::string& EB_EE ) {

  std::vector<float> return_vector;

  std::string histoThreshName = varName+"Thresh" + EB_EE + "_vs_rho"+rhoType;
  std::string histoEffName = varName+"Eff" + EB_EE + "_vs_rho"+rhoType;


  TH1D* histoThresh = (TH1D*)file->Get(histoThreshName.c_str());
  TH1D* histoEff = (TH1D*)file->Get(histoEffName.c_str());


  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();

  TF1* constLine = new TF1("constLine", "[0]");
  constLine->SetLineStyle(2);
  constLine->SetLineWidth(1);
  if( varName=="etawid" && EB_EE=="EE" ) {
    constLine->SetRange(0., 14.);
  }

  histoEff->Fit( constLine, "Q" );


  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();


  // first draw efficiency:

  TH2D* h2_axes_eff = new TH2D("axes_eff", "", 10, histoThresh->GetXaxis()->GetXmin(), histoThresh->GetXaxis()->GetXmax(), 10, 0.7, 1.1 );
  if( rhoType=="PF" )
    h2_axes_eff->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes_eff->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  h2_axes_eff->SetYTitle("Efficiency");
  h2_axes_eff->Draw();

  histoEff->SetMarkerStyle(29);
  histoEff->SetMarkerSize(1.8);
  histoEff->SetMarkerColor(38);
  
  histoEff->Draw("p same");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");

  TPaveText* label_explain = new TPaveText(0.25, 0.31, 0.85, 0.45, "brNDC");
  label_explain->SetFillColor(0);
  label_explain->SetTextSize(0.035);
  std::string explainText =  "Efficency of dynamic cut on " + yAxisName;
  label_explain->AddText(explainText.c_str());
  label_explain->Draw("same");
  
  TPaveText* label_EBEE = new TPaveText(0.7, 0.85, 0.83, 0.88, "brNDC");
  label_EBEE->SetFillColor(0);
  label_EBEE->SetTextSize(0.04);
  if( EB_EE == "EB" ) 
    label_EBEE->AddText( "ECAL Barrel" );
  else if( EB_EE == "EE" )
    label_EBEE->AddText( "ECAL Endcaps" );
  else if( EB_EE == "EB1" )
    label_EBEE->AddText( "|#eta| < 1" );
  else if( EB_EE == "EB2" )
    label_EBEE->AddText( "1 < |#eta| < 1.4" );
  label_EBEE->Draw("same");

  TPaveText* label_fitResults_eff = new TPaveText(0.3, 0.22, 0.7, 0.3, "brNDC");
  label_fitResults_eff->SetFillColor(0);
  label_fitResults_eff->SetTextSize(0.035);
  char label_fitResultsText_eff[300];
  sprintf( label_fitResultsText_eff, "FIT = %.3f #pm %.3f", constLine->GetParameter(0), constLine->GetParError(0) );
  label_fitResults_eff->AddText(label_fitResultsText_eff);
  label_fitResults_eff->Draw("same");

  std::string canvasName_eff = db->get_outputdir() + "/" + varName + "Eff" + EB_EE + "_vs_" + rhoType;
  std::string canvasName_eff_png = canvasName_eff + ".png";
  std::string canvasName_eff_eps = canvasName_eff + ".eps";

  c1->SaveAs(canvasName_eff_png.c_str());
  c1->SaveAs(canvasName_eff_eps.c_str());

  
  c1->Clear();
  

  // then the thresholds:

  std::string lineName = "line"+rhoType;
  TF1* fitLine = new TF1("fitLine", "[0] + [1]*x", 0., histoThresh->GetXaxis()->GetXmax());
  fitLine->SetLineStyle(2);
  fitLine->SetLineWidth(1);
  // exclude last point, fit more stable:
  fitLine->SetRange(0., 14.);

  //if( varName=="twrHcalIso" && EB_EE=="EB" ) {
//}
//if( varName=="etawid" && EB_EE=="EE" ) {
//  fitLine->SetRange(0., 14.);
//}
  std::cout << std::endl << "#### " << varName << " " << EB_EE << " rho" << rhoType << std::endl;
  histoThresh->Fit(fitLine, "R");


  TH2D* h2_axes = new TH2D("axes", "", 10, histoThresh->GetXaxis()->GetXmin(), histoThresh->GetXaxis()->GetXmax(), 10, 0., 1.5*histoThresh->GetMaximum() );
  if( rhoType=="PF" )
    h2_axes->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  std::string yAxisTitle = yAxisName + " Threshold";
  if( units!="" ) yAxisTitle += " [" + units + "]";
  h2_axes->SetYTitle(yAxisTitle.c_str());
  h2_axes->Draw();

  histoThresh->SetMarkerStyle(29);
  histoThresh->SetMarkerSize(1.8);
  histoThresh->SetMarkerColor(kRed+1);
  
  histoThresh->Draw("Psame");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");

  label_EBEE->Draw("same");

  //TPaveText* label_fitResults = new TPaveText(0.3, 0.65, 0.6, 0.75, "brNDC");
  TPaveText* label_fitResults = new TPaveText(0.55, 0.25, 0.85, 0.35, "brNDC");
  label_fitResults->SetFillColor(0);
  label_fitResults->SetTextSize(0.04);
  char label_fitResultsText[300];
  sprintf( label_fitResultsText, "y = %.3f x + %.3f %s", fitLine->GetParameter(1), fitLine->GetParameter(0), units.c_str() );
  label_fitResults->AddText(label_fitResultsText);
  label_fitResults->Draw("same");


  std::string canvasName = db->get_outputdir() + "/" + varName + "Thresh" + EB_EE + "_vs_" + rhoType;
  std::string canvasName_png = canvasName + ".png";
  std::string canvasName_eps = canvasName + ".eps";

  c1->SaveAs(canvasName_png.c_str());
  c1->SaveAs(canvasName_eps.c_str());

  return_vector.push_back(fitLine->GetParameter(0));
  return_vector.push_back(fitLine->GetParameter(1));

  return return_vector;

}


void drawSinglePlot_Mean( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisName, const std::string& units, const std::string& rhoType, const std::string& EB_EE ) {

  std::string histoNameMean = varName+"Mean" + EB_EE + "_vs_rho"+rhoType;
  std::string histoNameEff = "eff_"+varName+EB_EE + rhoType+"_stoeckli";

  TH1D* histoMean = (TH1D*)file->Get(histoNameMean.c_str());
  TH1D* histoEff = (TH1D*)file->Get(histoNameEff.c_str());

  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();



  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();


  TPaveText* label_EBEE = new TPaveText(0.72, 0.86, 0.85, 0.89, "brNDC");
  label_EBEE->SetFillColor(0);
  label_EBEE->SetTextSize(0.04);
  if( EB_EE == "EB" ) 
    label_EBEE->AddText( "ECAL Barrel" );
  else
    label_EBEE->AddText( "ECAL Endcaps" );
  label_EBEE->Draw("same");


  std::string lineName = "line"+rhoType;
  TF1* fitLine = new TF1("fitLine", "[0] + [1]*x", 0., histoMean->GetXaxis()->GetXmax());
  fitLine->SetLineStyle(2);
  fitLine->SetLineWidth(1);
  if( varName=="etawid" && EB_EE=="EE" ) {
    fitLine->SetRange(0., 14.);
  }
  std::cout << std::endl << "#### " << varName << " " << EB_EE << " rho" << rhoType << std::endl;
  histoMean->Fit(fitLine, "R");


  TH2D* h2_axes = new TH2D("axes", "", 10, histoMean->GetXaxis()->GetXmin(), histoMean->GetXaxis()->GetXmax(), 10, 0., 1.5*histoMean->GetMaximum() );
  if( rhoType=="PF" )
    h2_axes->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  std::string yAxisTitle = yAxisName + " Mean";
  if( units!="" ) yAxisTitle += " [" + units + "]";
  h2_axes->SetYTitle(yAxisTitle.c_str());
  h2_axes->Draw();

  histoMean->SetMarkerStyle(29);
  histoMean->SetMarkerSize(1.8);
  histoMean->SetMarkerColor(38);
  
  histoMean->Draw("Psame");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");

  label_EBEE->Draw("same");

  //TPaveText* label_fitResults = new TPaveText(0.3, 0.65, 0.6, 0.75, "brNDC");
  TPaveText* label_fitResults = new TPaveText(0.55, 0.25, 0.85, 0.35, "brNDC");
  label_fitResults->SetFillColor(0);
  label_fitResults->SetTextSize(0.04);
  char label_fitResultsText[300];
  sprintf( label_fitResultsText, "y = %.3f x + %.3f %s", fitLine->GetParameter(1), fitLine->GetParameter(0), units.c_str() );
  label_fitResults->AddText(label_fitResultsText);
  label_fitResults->Draw("same");


  std::string canvasName = db->get_outputdir() + "/" + varName + "Mean" + EB_EE + "_vs_" + rhoType;
  std::string canvasName_png = canvasName + ".png";
  std::string canvasName_eps = canvasName + ".eps";

  c1->SaveAs(canvasName_png.c_str());
  c1->SaveAs(canvasName_eps.c_str());


  c1->Clear();


  // now efficiency:

  TH2D* h2_axes_eff = new TH2D("axes_eff", "", 10, histoEff->GetXaxis()->GetXmin(), histoEff->GetXaxis()->GetXmax(), 10, 0.7, 1.1 );
  if( rhoType=="PF" )
    h2_axes_eff->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes_eff->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  h2_axes_eff->SetYTitle("Efficiency");
  h2_axes_eff->Draw();

  histoEff->SetMarkerStyle(29);
  histoEff->SetMarkerSize(1.8);
  histoEff->SetMarkerColor(38);
  
  histoEff->Draw("p same");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");
  
  label_EBEE->Draw("same");

  std::string canvasName_eff = db->get_outputdir() + "/" + varName + "StoeckliEff" + EB_EE + "_vs_" + rhoType;
  std::string canvasName_eff_png = canvasName_eff + ".png";
  std::string canvasName_eff_eps = canvasName_eff + ".eps";

  c1->SaveAs(canvasName_eff_png.c_str());
  c1->SaveAs(canvasName_eff_eps.c_str());

}



void drawCompare( DrawBase* db, TFile* file, const std::string& varName, const std::string& yAxisName, const std::string& rhoType, const std::string& EB_EE ) {

  std::string histoNameEff_noCorr = "eff_"+varName+EB_EE + rhoType+"_noCorr";
  std::string histoNameEff_stoeckli = "eff_"+varName+EB_EE + rhoType+"_stoeckli";
  std::string histoNameEff_isoEff = "eff_"+varName+EB_EE + rhoType+"_isoEff";

  TH1D* histoEff_noCorr = (TH1D*)file->Get(histoNameEff_noCorr.c_str());
  TH1D* histoEff_stoeckli = (TH1D*)file->Get(histoNameEff_stoeckli.c_str());
  TH1D* histoEff_isoEff = (TH1D*)file->Get(histoNameEff_isoEff.c_str());

  float markerSize = 1.6;

  histoEff_noCorr->SetMarkerStyle(21);
  histoEff_noCorr->SetMarkerSize(markerSize);
  histoEff_noCorr->SetMarkerColor(kGray+2);

  histoEff_stoeckli->SetMarkerStyle(20);
  histoEff_stoeckli->SetMarkerSize(markerSize);
  histoEff_stoeckli->SetMarkerColor(38);

  histoEff_isoEff->SetMarkerStyle(22);
  histoEff_isoEff->SetMarkerSize(markerSize);
  histoEff_isoEff->SetMarkerColor(46);


  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();



  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();


  TPaveText* label_EBEE = new TPaveText(0.72, 0.87, 0.85, 0.9, "brNDC");
  label_EBEE->SetFillColor(0);
  label_EBEE->SetTextSize(0.04);
  if( EB_EE == "EB" ) 
    label_EBEE->AddText( "ECAL Barrel" );
  else
    label_EBEE->AddText( "ECAL Endcaps" );
  label_EBEE->Draw("same");



  float effMin = 0.75;
  if( EB_EE=="EB" && varName=="jurEcalIso" ) effMin=0.5;
  if( varName=="hlwTrackIso" ) effMin=0.;
  TH2D* h2_axes_eff = new TH2D("axes_eff", "", 10, histoEff_isoEff->GetXaxis()->GetXmin(), histoEff_isoEff->GetXaxis()->GetXmax(), 10, effMin, 1.1 );
  if( rhoType=="PF" )
    h2_axes_eff->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
  else
    h2_axes_eff->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  h2_axes_eff->SetYTitle("Efficiency");
  h2_axes_eff->Draw();


  TLegend* legend = new TLegend(0.2, 0.2, 0.55, 0.45, yAxisName.c_str());
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( histoEff_noCorr, "No Correction", "P" );
  legend->AddEntry( histoEff_stoeckli, "Mean Correction", "P" );
  legend->AddEntry( histoEff_isoEff, "Iso-Efficiency", "P" );
  legend->Draw("same");

  histoEff_noCorr->Draw("Psame");
  histoEff_stoeckli->Draw("Psame");
  histoEff_isoEff->Draw("Psame");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");

  label_EBEE->Draw("same");


  std::string canvasName = db->get_outputdir() + "/compare_" + varName + EB_EE + "_vs_" + rhoType;
  std::string canvasName_png = canvasName + ".png";
  std::string canvasName_eps = canvasName + ".eps";

  c1->SaveAs(canvasName_png.c_str());
  c1->SaveAs(canvasName_eps.c_str());


}



void drawCompare_fullSelection( DrawBase* db, TFile* file, const std::string& rhoType, const std::string& EB_EE, bool vs_nvertex ) {


  TH1D* h1_fullEff_noCorr;
  TH1D* h1_fullEff_stoeckli;
  TH1D* h1_fullEff_isoEff;

  std::vector<std::string> varNames;
  varNames.push_back( "twrHcalIso" );
  varNames.push_back( "jurEcalIso" );
  varNames.push_back( "hlwTrackIso" );
  varNames.push_back( "HoverE" );

  for( unsigned iVar=0; iVar<varNames.size(); ++iVar ) {

    std::string histoName_noCorr = "eff_" + varNames[iVar] + EB_EE + rhoType + "_noCorr";
    std::string histoName_stoeckli = "eff_" + varNames[iVar] + EB_EE + rhoType + "_stoeckli";
    std::string histoName_isoEff = "eff_" + varNames[iVar] + EB_EE + rhoType + "_isoEff";

    if( vs_nvertex ) {
      histoName_noCorr   += "_vs_nvertex";
      histoName_stoeckli += "_vs_nvertex";
      histoName_isoEff   += "_vs_nvertex";
    }

    TH1D* h1_noCorr = (TH1D*)file->Get(histoName_noCorr.c_str());
    TH1D* h1_stoeckli = (TH1D*)file->Get(histoName_stoeckli.c_str());
    TH1D* h1_isoEff = (TH1D*)file->Get(histoName_isoEff.c_str());
    
    if( iVar==0 ) {

      h1_fullEff_noCorr = new TH1D("fullEff_noCorr", "", h1_noCorr->GetXaxis()->GetNbins(), h1_noCorr->GetXaxis()->GetXmin(), h1_noCorr->GetXaxis()->GetXmax() );
      h1_fullEff_stoeckli = new TH1D("fullEff_stoeckli", "", h1_noCorr->GetXaxis()->GetNbins(), h1_noCorr->GetXaxis()->GetXmin(), h1_noCorr->GetXaxis()->GetXmax() );
      h1_fullEff_isoEff = new TH1D("fullEff_isoEff", "", h1_noCorr->GetXaxis()->GetNbins(), h1_noCorr->GetXaxis()->GetXmin(), h1_noCorr->GetXaxis()->GetXmax() );
 
      for( unsigned iBin=0; iBin<h1_noCorr->GetXaxis()->GetNbins(); ++iBin ) {
        h1_fullEff_noCorr->SetBinContent(iBin+1, 1.);
        h1_fullEff_stoeckli->SetBinContent(iBin+1, 1.);
        h1_fullEff_isoEff->SetBinContent(iBin+1, 1.);
      }

    }

    for( unsigned iBin=0; iBin<h1_noCorr->GetXaxis()->GetNbins(); ++iBin ) {
      h1_fullEff_noCorr->SetBinContent(iBin+1, h1_fullEff_noCorr->GetBinContent(iBin+1)*h1_noCorr->GetBinContent(iBin+1) );
      h1_fullEff_stoeckli->SetBinContent(iBin+1, h1_fullEff_stoeckli->GetBinContent(iBin+1)*h1_stoeckli->GetBinContent(iBin+1) );
      h1_fullEff_isoEff->SetBinContent(iBin+1, h1_fullEff_isoEff->GetBinContent(iBin+1)*h1_isoEff->GetBinContent(iBin+1) );
    }

  }
    


  float markerSize = 1.6;

  h1_fullEff_noCorr->SetMarkerStyle(21);
  h1_fullEff_noCorr->SetMarkerSize(markerSize);
  h1_fullEff_noCorr->SetMarkerColor(kGray+2);

  h1_fullEff_stoeckli->SetMarkerStyle(20);
  h1_fullEff_stoeckli->SetMarkerSize(markerSize);
  h1_fullEff_stoeckli->SetMarkerColor(38);

  h1_fullEff_isoEff->SetMarkerStyle(22);
  h1_fullEff_isoEff->SetMarkerSize(markerSize);
  h1_fullEff_isoEff->SetMarkerColor(46);


  TPaveText* cmsLabel = db->get_labelCMS();
  TPaveText* sqrtLabel = db->get_labelSqrt();



  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();


  TPaveText* label_EBEE = new TPaveText(0.72, 0.87, 0.85, 0.9, "brNDC");
  label_EBEE->SetFillColor(0);
  label_EBEE->SetTextSize(0.04);
  if( EB_EE == "EB" ) 
    label_EBEE->AddText( "ECAL Barrel" );
  else
    label_EBEE->AddText( "ECAL Endcaps" );
  label_EBEE->Draw("same");



  TH2D* h2_axes_eff = new TH2D("axes_eff", "", 10, h1_fullEff_isoEff->GetXaxis()->GetXmin(), h1_fullEff_isoEff->GetXaxis()->GetXmax(), 10, 0., 1.1 );
  if( vs_nvertex ) {
    h2_axes_eff->SetXTitle("Number of Reconstructed Vertexes");
  } else {
    if( rhoType=="PF" )
      h2_axes_eff->SetXTitle("Particle Flow Energy Density (#rho) [GeV]");
    else
      h2_axes_eff->SetXTitle("Calorimeter Energy Density (#rho) [GeV]");
  }
  h2_axes_eff->SetYTitle("Efficiency");
  h2_axes_eff->Draw();


  TLegend* legend = new TLegend(0.2, 0.2, 0.55, 0.45, "Full Selection" );
  legend->SetFillColor(0);
  legend->SetTextSize(0.035);
  legend->AddEntry( h1_fullEff_noCorr, "No Correction", "P" );
  legend->AddEntry( h1_fullEff_stoeckli, "Mean Correction", "P" );
  legend->AddEntry( h1_fullEff_isoEff, "Iso-Efficiency", "P" );
  legend->Draw("same");

  h1_fullEff_noCorr->Draw("Psame");
  h1_fullEff_stoeckli->Draw("Psame");
  h1_fullEff_isoEff->Draw("Psame");

  cmsLabel->Draw("same");
  sqrtLabel->Draw("same");

  label_EBEE->Draw("same");


  std::string canvasName = db->get_outputdir() + "/compareFullEff_" + EB_EE + "_vs_";
  if( vs_nvertex ) canvasName += "nvertex";
  else             canvasName += rhoType;
  std::string canvasName_png = canvasName + ".png";
  std::string canvasName_eps = canvasName + ".eps";

  c1->SaveAs(canvasName_png.c_str());
  c1->SaveAs(canvasName_eps.c_str());


}

