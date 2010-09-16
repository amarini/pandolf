#include "DrawExtrap.h"
#include "fitTools.h"

#include "TGraphErrors.h"
#include "TH1F.h"




DrawExtrap::DrawExtrap( const std::string& analysisType, const std::string& recoType, const std::string& jetAlgo ) : DrawBase( analysisType, recoType, jetAlgo ) {

  FIT_RMS_ = "FIT";
  NOQ_ = false;
  INTPERC_ = 95.;

}



void DrawExtrap::drawResponseExtrap( const std::string& recoGen ) const {

  std::vector<float> ptPhot_binning = fitTools::getPtPhot_binning();

//Float_t ptPhotRecoDATA[NBINS_PT-1];
//Float_t ptPhotReco_errDATA[NBINS_PT-1];

//Float_t ptPhotReco[NBINS_PT-1];
//Float_t ptPhotReco_err[NBINS_PT-1];

//Float_t extrapReso_RecoRelDATA[NBINS_PT-1];
//Float_t extrapReso_err_RecoRelDATA[NBINS_PT-1];
//Float_t extrapReso_RecoRel[NBINS_PT-1];
//Float_t extrapReso_err_RecoRel[NBINS_PT-1];
//Float_t trueReso_RecoRel[NBINS_PT-1];
//Float_t trueReso_err_RecoRel[NBINS_PT-1];
//Float_t intrReso_RecoRel[NBINS_PT-1];
//Float_t intrReso_err_RecoRel[NBINS_PT-1];
//Float_t extrapResp_RecoRelDATA[NBINS_PT-1];
//Float_t extrapResp_err_RecoRelDATA[NBINS_PT-1];
//Float_t extrapResp_RecoRel[NBINS_PT-1];
//Float_t extrapResp_err_RecoRel[NBINS_PT-1];
//Float_t trueResp_RecoRel[NBINS_PT-1];
//Float_t trueResp_err_RecoRel[NBINS_PT-1];
//Float_t intrResp_RecoRel[NBINS_PT-1];
//Float_t intrResp_err_RecoRel[NBINS_PT-1];
//Float_t imbalanceResp_RecoRel[NBINS_PT-1];
//Float_t imbalanceResp_err_RecoRel[NBINS_PT-1];
//Float_t pullResp_RecoRel[NBINS_PT-1];
//Float_t pullResp_err_RecoRel[NBINS_PT-1];
//Float_t pullReso_RecoRel[NBINS_PT-1];
//Float_t pullReso_err_RecoRel[NBINS_PT-1];

  TGraphErrors* gr_DATAResp_vs_pt = new TGraphErrors(0);
  gr_DATAResp_vs_pt->SetName("gr_DATAResp_vs_pt");
  gr_DATAResp_vs_pt->SetMarkerStyle(25);
  gr_DATAResp_vs_pt->SetMarkerColor(kBlack);
  gr_DATAResp_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_extrapResp_vs_pt = new TGraphErrors(0);
  gr_extrapResp_vs_pt->SetName("gr_extrapResp_vs_pt");
  gr_extrapResp_vs_pt->SetMarkerStyle(25);
  gr_extrapResp_vs_pt->SetMarkerColor(kBlack);
  gr_extrapResp_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_intrResp_vs_pt = new TGraphErrors(0);
  gr_intrResp_vs_pt->SetName("gr_intrResp_vs_pt");
  gr_intrResp_vs_pt->SetMarkerStyle(29);
  gr_intrResp_vs_pt->SetMarkerColor(kBlue);
  gr_intrResp_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_trueResp_vs_pt = new TGraphErrors(0);
  gr_trueResp_vs_pt->SetName("gr_trueResp_vs_pt");
  gr_trueResp_vs_pt->SetMarkerStyle(23);
  gr_trueResp_vs_pt->SetMarkerColor(kGray+2);
  gr_trueResp_vs_pt->SetMarkerSize(1.5);



  TGraphErrors* gr_DATAReso_vs_pt = new TGraphErrors(0);
  gr_DATAReso_vs_pt->SetName("gr_DATAReso_vs_pt");
  gr_DATAReso_vs_pt->SetMarkerStyle(25);
  gr_DATAReso_vs_pt->SetMarkerColor(kBlack);
  gr_DATAReso_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_extrapReso_vs_pt = new TGraphErrors(0);
  gr_extrapReso_vs_pt->SetName("gr_extrapReso_vs_pt");
  gr_extrapReso_vs_pt->SetMarkerStyle(25);
  gr_extrapReso_vs_pt->SetMarkerColor(kBlack);
  gr_extrapReso_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_intrReso_vs_pt = new TGraphErrors(0);
  gr_intrReso_vs_pt->SetName("gr_intrReso_vs_pt");
  gr_intrReso_vs_pt->SetMarkerStyle(29);
  gr_intrReso_vs_pt->SetMarkerColor(kBlue);
  gr_intrReso_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_trueReso_vs_pt = new TGraphErrors(0);
  gr_trueReso_vs_pt->SetName("gr_trueReso_vs_pt");
  gr_trueReso_vs_pt->SetMarkerStyle(23);
  gr_trueReso_vs_pt->SetMarkerColor(kGray+2);
  gr_trueReso_vs_pt->SetMarkerSize(1.5);

  TGraphErrors* gr_imbalance_vs_pt = new TGraphErrors(0);
  gr_imbalance_vs_pt->SetName("gr_imbalance_vs_pt");
  gr_imbalance_vs_pt->SetMarkerStyle(20);
  gr_imbalance_vs_pt->SetMarkerColor(kBlack);
  gr_imbalance_vs_pt->SetMarkerSize(1.5);


  TH2D* h2_ptPhotReco_vs_pt = (TH2D*)(get_mcFile(0))->Get("ptPhotMean_no2ndJet");
  TH2D* h2_ptPhotReco_vs_ptDATA = (get_dataFile()==0) ? 0 : (TH2D*)(get_dataFile())->Get("ptPhotMean_no2ndJet");


  for( int iPtBin=0; iPtBin<(ptPhot_binning.size()-3); //-3 instead of -1 (extrap reaches up to ~2 less bins in pt wrt balancing)
       ++iPtBin) {

    char projName[50];
    sprintf(projName, "projection_%d",iPtBin);

    TH1D* h1_proj = h2_ptPhotReco_vs_pt->ProjectionY(projName, iPtBin+1, iPtBin+1);
    float ptPhotReco_thisBin = h1_proj->GetMean();
    float ptPhotReco_err_thisBin = (h1_proj->GetEntries()>1.) ? h1_proj->GetRMS()/sqrt(h1_proj->GetEntries()) : h1_proj->GetRMS();

    TH1D* h1_projDATA = h2_ptPhotReco_vs_ptDATA->ProjectionY(projName, iPtBin+1, iPtBin+1);
    float ptPhotReco_thisBinDATA = h1_projDATA->GetMean();
    float ptPhotReco_err_thisBinDATA = (h1_projDATA->GetEntries()>1.) ? h1_projDATA->GetRMS()/sqrt(h1_projDATA->GetEntries()) : h1_projDATA->GetRMS();



    Double_t ptMin = ptPhot_binning[iPtBin];
    Double_t ptMax = ptPhot_binning[iPtBin+1];

   
    Int_t nPoints = 5;

    char xHistoName[90];
    sprintf( xHistoName, "extrap_ptBin_%d_%d/pt2ndJet%sMean_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);

    Float_t x[nPoints];
    Float_t x_err[nPoints];
    getXPoints( get_mcFile(0), xHistoName, nPoints, x, x_err);

    Float_t xDATA[nPoints];
    Float_t x_errDATA[nPoints];
    getXPoints( get_dataFile(), xHistoName, nPoints, xDATA, x_errDATA);

    Float_t y_resp_DATA[nPoints];
    Float_t y_resp_err_DATA[nPoints];

    Float_t y_resp_recoPhot[nPoints];
    Float_t y_resp_recoPhot_err[nPoints];

    Float_t y_resp_genPhot[nPoints];
    Float_t y_resp_genPhot_err[nPoints];

    Float_t y_resp_recoGen[nPoints];
    Float_t y_resp_recoGen_err[nPoints];

    Float_t y_reso_DATA[nPoints];
    Float_t y_reso_err_DATA[nPoints];

    Float_t y_reso_recoPhot[nPoints];
    Float_t y_reso_recoPhot_err[nPoints];

    Float_t y_reso_genPhot[nPoints];
    Float_t y_reso_genPhot_err[nPoints];

    Float_t y_reso_recoGen[nPoints];
    Float_t y_reso_recoGen_err[nPoints];

    char yHistoName[100];
    sprintf( yHistoName, "extrap_ptBin_%d_%d/r_RecoPhot_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
    getYPoints( get_dataFile(), yHistoName, nPoints, y_resp_DATA, y_resp_err_DATA,  y_reso_DATA, y_reso_err_DATA);
    sprintf( yHistoName, "extrap_ptBin_%d_%d/r_RecoPhot_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
    getYPoints( get_mcFile(0), yHistoName, nPoints, y_resp_recoPhot, y_resp_recoPhot_err,  y_reso_recoPhot, y_reso_recoPhot_err);
    sprintf( yHistoName, "extrap_ptBin_%d_%d/r_GenPhot_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
    getYPoints( get_mcFile(0), yHistoName, nPoints, y_resp_genPhot, y_resp_genPhot_err,  y_reso_genPhot, y_reso_genPhot_err);
    sprintf( yHistoName, "extrap_ptBin_%d_%d/r_RecoGen_vs_%s_%d", (int)ptMin, (int)ptMax, recoGen.c_str(), iPtBin);
    getYPoints( get_mcFile(0), yHistoName, nPoints, y_resp_recoGen, y_resp_recoGen_err,  y_reso_recoGen, y_reso_recoGen_err);

   
    //draw response histograms:

    TGraphErrors* gr_resp_DATA = new TGraphErrors(nPoints, xDATA, y_resp_DATA, x_errDATA, y_resp_err_DATA);
    gr_resp_DATA->SetMarkerStyle(20);
    gr_resp_DATA->SetMarkerColor(kRed);

    TGraphErrors* gr_resp_recoPhot = new TGraphErrors(nPoints, x, y_resp_recoPhot, x_err, y_resp_recoPhot_err);
    gr_resp_recoPhot->SetMarkerStyle(24);
    gr_resp_recoPhot->SetMarkerColor(kRed);

    TGraphErrors* gr_resp_genPhot = new TGraphErrors(nPoints, x, y_resp_genPhot, x_err, y_resp_genPhot_err);
    gr_resp_genPhot->SetMarkerStyle(22);
    gr_resp_genPhot->SetMarkerColor(kBlack);

    TGraphErrors* gr_resp_recoGen = new TGraphErrors(nPoints, x, y_resp_recoGen, x_err, y_resp_recoGen_err);
    gr_resp_recoGen->SetMarkerStyle(21);
    gr_resp_recoGen->SetMarkerColor(kBlue);


    Float_t lastX = x[nPoints-1];
    Float_t xMax_fit = lastX+2.;
    Float_t xMax_fit_green = x[2]+2.;

    Float_t xMax_axis;
    if( lastX <=12. ) 
      xMax_axis = 15.;
    else if( lastX <= 20. )
      xMax_axis = 25.;
    else if( lastX <= 25. )
      xMax_axis = 30.;
    else if( lastX <= 35. )
      xMax_axis = 40.;
    

    std::string xTitle;
    if( recoGen=="Reco" ) {
      xTitle = "p_{T}^{2ndJet reco} [GeV/c]";
    } else if( recoGen=="Gen" ) {
      xTitle = "p_{T}^{2ndJet gen} [GeV/c]";
    } else if( recoGen=="RecoRel" ) {
      xTitle = "p_{T}^{2ndJet reco}/p_{T}^{#gamma} [%]";
    } else if( recoGen=="GenRel" ) {
      xTitle = "p_{T}^{2ndJet gen}/p_{T}^{#gamma} [%]";
    }


    std::string fitFunct_name;
////  if( ptMin <=150. ) 
      fitFunct_name = "[0] - x*x*[1]";
////  else 
////    fitFunct_name = "[0] - x*[1]";


    TF1* fit_resp_genPhot = new TF1("fit_resp_genPhot", fitFunct_name.c_str());
    fit_resp_genPhot->SetRange(0., xMax_fit);
    fit_resp_genPhot->SetLineWidth(0.5);
    fit_resp_genPhot->SetLineColor(kBlack);
    gr_resp_genPhot->Fit( fit_resp_genPhot, "RQ");

    TF1* fit_resp_recoGen = new TF1("fit_resp_recoGen", "[0]");
    fit_resp_recoGen->SetRange(0., xMax_fit);
    fit_resp_recoGen->SetLineWidth(0.5);
    fit_resp_recoGen->SetLineColor(kBlue);
    gr_resp_recoGen->Fit( fit_resp_recoGen, "RQ");


    std::string total_resp_str = "fit_resp_recoGen*fit_resp_genPhot";
    TF1* total_resp = new TF1("total_resp", total_resp_str.c_str());
    total_resp->SetRange(0., xMax_fit);
    total_resp->SetLineColor(kGray+2);

   
    fitFunct_name = "[0]*[1] - x*x*[2]";
    TF1* fit_respParabola = new TF1("fit_respParabola", fitFunct_name.c_str());
    fit_respParabola->SetRange(0., xMax_fit);
    if( NOQ_ ) //to evaluate syst
      fit_respParabola->FixParameter(1, 0.5*fit_resp_genPhot->GetParameter(0));
    else
      fit_respParabola->FixParameter(1, fit_resp_genPhot->GetParameter(0));
    fit_respParabola->SetLineColor(2);
    fit_respParabola->SetLineColor(kRed);
    fit_respParabola->SetLineWidth(0.5);
    gr_resp_recoPhot->Fit( fit_respParabola, "RQ" );

    TF1* fit_respParabola_DATA = new TF1("fit_respParabola_DATA", fitFunct_name.c_str());
    fit_respParabola_DATA->SetRange(0., xMax_fit);
    if( NOQ_ ) //to evaluate syst
      fit_respParabola_DATA->FixParameter(1, 0.5*fit_resp_genPhot->GetParameter(0));
    else
      fit_respParabola_DATA->FixParameter(1, fit_resp_genPhot->GetParameter(0));
    fit_respParabola_DATA->SetLineColor(kRed);
    fit_respParabola_DATA->SetLineWidth(3);
    fit_respParabola_DATA->SetLineWidth(0.5);
    gr_resp_DATA->Fit( fit_respParabola_DATA, "RQ" );



    // set response graph points:

    gr_DATAResp_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBinDATA, fit_respParabola_DATA->GetParameter(0) );
    gr_DATAResp_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBinDATA, fit_respParabola_DATA->GetParError(0) );

    gr_extrapResp_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBin, fit_respParabola->GetParameter(0) );
    gr_extrapResp_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, fit_respParabola->GetParError(0) );

    gr_intrResp_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBin, fit_resp_recoGen->GetParameter(0) );
    gr_intrResp_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, fit_resp_recoGen->GetParError(0) );

    gr_imbalance_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBin, 100.*(fit_resp_genPhot->GetParameter(0) - 1.) ); //in percent
    gr_imbalance_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, 100.*(fit_resp_genPhot->GetParError(0)) );

    Float_t y0_recoGen = fit_resp_recoGen->GetParameter(0);
    Float_t y0_genPhot = fit_resp_genPhot->GetParameter(0);
    Float_t y0_err_recoGen = fit_resp_recoGen->GetParError(0);
    Float_t y0_err_genPhot = fit_resp_genPhot->GetParError(0);
    
    gr_trueResp_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBin, y0_genPhot*y0_recoGen );
    float trueResp_err = y0_err_genPhot*y0_err_genPhot*y0_recoGen*y0_recoGen;
    trueResp_err += y0_genPhot*y0_genPhot*y0_err_recoGen*y0_err_recoGen;
    trueResp_err = sqrt( trueResp_err );
    gr_trueResp_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, trueResp_err );


    Float_t yMin_axis = (ptPhot_binning[iPtBin]<50.) ? 0.7 : 0.8;
    if( ptPhot_binning[iPtBin]<20. ) yMin_axis = 0.6;
    TH2D* h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, yMin_axis, 1.1);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->GetXaxis()->SetTitleOffset(1.1);
    h2_axes_resp->SetYTitle("Response");
    h2_axes_resp->GetYaxis()->SetTitleOffset(1.2);

    TLegend* legend_resp = new TLegend(0.6, 0.65, 0.88, 0.88);
    legend_resp->SetTextSize(0.04);
    legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_resp_recoGen, "MC Intrinsic", "P");
    legend_resp->AddEntry(gr_resp_genPhot, "MC Imbalance", "P");
    legend_resp->AddEntry(total_resp, "Intr #oplus Imb", "L");
    legend_resp->AddEntry(gr_resp_recoPhot, "MC (#gamma + jet)", "P");
    legend_resp->AddEntry(gr_resp_DATA, "DATA (#gamma + jet)", "P");

    char labeltext[50];
    sprintf(labeltext, "%d < p_{T}^{#gamma} < %d GeV/c", (int)ptMin, (int)ptMax);  
    TPaveText* label_resp = new TPaveText(0.15, 0.15, 0.4, 0.2, "brNDC");
    label_resp->SetFillColor(kWhite);
    label_resp->SetTextSize(0.035);
    label_resp->AddText(labeltext);

    TCanvas* c1_resp = new TCanvas("c1_resp", "c1_resp", 800, 600);
    c1_resp->cd();
    h2_axes_resp->Draw();
    total_resp->Draw("same");
    gr_resp_recoGen->Draw("Psame");
    gr_resp_genPhot->Draw("Psame");
    gr_resp_recoPhot->Draw("Psame");
    gr_resp_DATA->Draw("Psame");
    legend_resp->Draw("same");
    label_resp->Draw("same");

    char canvasName_resp_eps[150];
    char canvasName_resp_png[150];
    sprintf(canvasName_resp_eps, "%s/response_ptBin_%d_%d_%s.eps", get_outputdir().c_str(), (int)ptMin, (int)ptMax, recoGen.c_str());
    sprintf(canvasName_resp_png, "%s/response_ptBin_%d_%d_%s.png", get_outputdir().c_str(), (int)ptMin, (int)ptMax, recoGen.c_str());
    c1_resp->SaveAs(canvasName_resp_eps);
    c1_resp->SaveAs(canvasName_resp_png);

    delete c1_resp;
    delete h2_axes_resp;
    delete legend_resp;
    delete label_resp;



    //draw resolution histograms:

    TGraphErrors* gr_reso_DATA = new TGraphErrors(nPoints, xDATA, y_reso_DATA, x_errDATA, y_reso_err_DATA);
    gr_reso_DATA->SetMarkerStyle(20);
    gr_reso_DATA->SetMarkerColor(kRed);

    TGraphErrors* gr_reso_recoPhot = new TGraphErrors(nPoints, x, y_reso_recoPhot, x_err, y_reso_recoPhot_err);
    gr_reso_recoPhot->SetMarkerStyle(24);
    gr_reso_recoPhot->SetMarkerColor(kRed);

    TGraphErrors* gr_reso_genPhot = new TGraphErrors(nPoints, x, y_reso_genPhot, x_err, y_reso_genPhot_err);
    gr_reso_genPhot->SetMarkerStyle(22);
    gr_reso_genPhot->SetMarkerColor(kBlack);

    TGraphErrors* gr_reso_recoGen = new TGraphErrors(nPoints, x, y_reso_recoGen, x_err, y_reso_recoGen_err);
    gr_reso_recoGen->SetMarkerStyle(21);
    gr_reso_recoGen->SetMarkerColor(kBlue);


    Double_t x1, x2, y1, y2;

    TF1* fit_reso_genPhot = new TF1("fit_reso_genPhot", "[0] + x*[1]");
    fit_reso_genPhot->SetRange(0., xMax_fit);
    fit_reso_genPhot->SetLineWidth(0.5);
    fit_reso_genPhot->SetLineColor(kBlack);
    gr_reso_genPhot->Fit(fit_reso_genPhot, "RQ");


    TF1* fit_reso_recoGen = new TF1("fit_reso_recoGen", "[0]");
    fit_reso_recoGen->SetRange(0., xMax_fit);
    gr_reso_recoGen->GetPoint(1, x1, y1);
    gr_reso_recoGen->GetPoint(2, x2, y2);
    fit_reso_recoGen->SetParameter(0, y1);
    fit_reso_recoGen->SetLineWidth(0.5);
    fit_reso_recoGen->SetLineColor(kBlue);
    gr_reso_recoGen->Fit(fit_reso_recoGen, "RQ");




    std::string sum_str;
    sum_str = "sqrt( fit_reso_recoGen*fit_reso_recoGen + fit_reso_genPhot*fit_reso_genPhot )";

    TF1* sum= new TF1("sum", sum_str.c_str());
    sum->SetRange(0., xMax_fit);
    sum->SetLineColor(kGray+2);

    Double_t x1_reco, y1_reco;
    gr_reso_recoPhot->GetPoint(1, x1_reco, y1_reco);
    Double_t x2_reco, y2_reco;
    gr_reso_recoPhot->GetPoint(2, x2_reco, y2_reco);

    // genPhot: y = mx + q
    // recoGen: y = c
    Float_t c = fit_reso_recoGen->GetParameter(0);
    Float_t q = fit_reso_genPhot->GetParameter(0);
    Float_t m = fit_reso_genPhot->GetParameter(1);
    // [0] = c; [1] = q; [2] = m
    TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "sqrt([0]*[0] + [1]*[1] + 2.*[1]*[2]*x + [2]*[2]*x*x)");
    fit_extrapToZero_sqrt->SetRange(0., xMax_fit);
    fit_extrapToZero_sqrt->SetParameter(0, c);
    fit_extrapToZero_sqrt->SetParLimits(0, 0.01, 0.3);
    if( NOQ_ )
      fit_extrapToZero_sqrt->FixParameter(1, 0.5*q); //to evaluate syst
    else
      fit_extrapToZero_sqrt->FixParameter(1, q); //fixed
    //fit_extrapToZero_sqrt->SetParLimits(1, 0., 0.1);
    fit_extrapToZero_sqrt->SetParameter(2, m);
    fit_extrapToZero_sqrt->SetParLimits(2, 0., 0.05);
    fit_extrapToZero_sqrt->SetLineStyle(2);
    fit_extrapToZero_sqrt->SetLineColor(kRed);
    fit_extrapToZero_sqrt->SetLineWidth(0.5);

    TF1* fit_extrapToZero_sqrt_DATA = new TF1(*fit_extrapToZero_sqrt);
    fit_extrapToZero_sqrt_DATA->SetName("fit_extrapToZero_sqrt_DATA");
    fit_extrapToZero_sqrt_DATA->SetLineStyle(1);
    fit_extrapToZero_sqrt_DATA->SetLineWidth(3);

    TF1* fit_extrapToZero_line = new TF1("fit_extrapToZero_line", "[0]+[1]*x");
    fit_extrapToZero_line->SetRange(0., xMax_fit);
    fit_extrapToZero_line->SetParameter(0, y1_reco);
    fit_extrapToZero_line->SetParameter(1, (y2_reco-y1_reco)/(x2_reco-x1_reco));
    fit_extrapToZero_line->SetParLimits(0, 0., 0.4);
    fit_extrapToZero_line->SetLineColor(kRed);
    fit_extrapToZero_line->SetLineWidth(0.5);

    gr_reso_recoPhot->Fit(fit_extrapToZero_sqrt, "RQ");

    gr_reso_DATA->Fit(fit_extrapToZero_sqrt_DATA, "RQ");

    gr_DATAReso_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBinDATA,  fit_extrapToZero_sqrt_DATA->GetParameter(0) );
    gr_DATAReso_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBinDATA, fit_extrapToZero_sqrt_DATA->GetParError(0) );

    gr_extrapReso_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt->GetParameter(0) );
    gr_extrapReso_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt->GetParError(0) );
    //extrapReso = sqrt(fit_extrapToZero_sqrt->GetParameter(0));
    //extrapReso_err = ( 1. / (2.*sqrt(extrapReso)) * fit_extrapToZero_sqrt->GetParError(0)); //error propag

    gr_intrReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_reso_recoGen->GetParameter(0) );
    gr_intrReso_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, fit_reso_recoGen->GetParError(0) );

    gr_trueReso_vs_pt->SetPoint( iPtBin, ptPhotReco_thisBin, sum->Eval(0.) );
    gr_trueReso_vs_pt->SetPointError( iPtBin, ptPhotReco_err_thisBin, sqrt( fit_reso_genPhot->GetParError(0)*fit_reso_genPhot->GetParError(0) + fit_reso_recoGen->GetParError(0)*fit_reso_recoGen->GetParError(0)) );



    Float_t ymax;
    if( ptPhot_binning[iPtBin] < 40. ) ymax = (get_recoType()=="pf") ? 0.5 : 0.8;
    else if( ptPhot_binning[iPtBin] < 80. ) ymax = (get_recoType()=="pf") ? 0.4 : 0.6;
    else ymax = 0.3;

    TH2D* h2_axes_reso = new TH2D("axes_reso", "", 10, 0., xMax_axis, 10, 0., ymax);
    h2_axes_reso->SetXTitle(xTitle.c_str());
    h2_axes_reso->GetXaxis()->SetTitleOffset(1.1);
    h2_axes_reso->SetYTitle("Resolution");
    h2_axes_reso->GetYaxis()->SetTitleOffset(1.2);

    TPaveText* label_reso = new TPaveText(0.67, 0.83, 0.8, 0.86, "brNDC");
    label_reso->SetFillColor(kWhite);
    label_reso->SetTextSize(0.035);
    label_reso->AddText(labeltext);

    Float_t minLegend = 0.2;
    TLegend* legend_reso = new TLegend(minLegend, 0.6, 0.45, 0.85);
    legend_reso->SetTextSize(0.04);
    legend_reso->SetFillStyle(0);
    legend_reso->SetFillColor(kWhite);
    legend_reso->AddEntry(gr_reso_recoGen, "MC Intrinsic", "P");
    legend_reso->AddEntry(gr_reso_genPhot, "MC Imbalance", "P");
    legend_reso->AddEntry(sum, "Intr #oplus Imb", "L");
    legend_reso->AddEntry(gr_reso_recoPhot, "MC (#gamma + jet)", "P");
    legend_reso->AddEntry(gr_reso_DATA, "DATA (#gamma + jet)", "P");

    TCanvas* c1_reso = new TCanvas("c1_reso", "c1_reso", 800, 600);
    c1_reso->cd();
    h2_axes_reso->Draw();
    sum->Draw("same");
    gr_reso_recoPhot->Draw("Psame");
    gr_reso_recoGen->Draw("Psame");
    gr_reso_genPhot->Draw("Psame");
    gr_reso_DATA->Draw("Psame");
    legend_reso->Draw("same");
    label_reso->Draw("same");


    char canvasName_reso_eps[150];
    char canvasName_reso_png[150];
    sprintf(canvasName_reso_eps, "%s/resolution_ptBin_%d_%d_%s.eps", get_outputdir().c_str(), (int)ptMin, (int)ptMax, recoGen.c_str());
    sprintf(canvasName_reso_png, "%s/resolution_ptBin_%d_%d_%s.png", get_outputdir().c_str(), (int)ptMin, (int)ptMax, recoGen.c_str());
    c1_reso->SaveAs(canvasName_reso_eps);
    c1_reso->SaveAs(canvasName_reso_png);

    delete c1_reso;
    delete h2_axes_reso;
    delete legend_reso;
    delete label_reso;


    delete gr_resp_recoPhot;
    delete gr_resp_recoGen;

    delete gr_reso_recoPhot;
    delete gr_reso_recoGen;

  //gr_pullResp_vs_pt->SetPoint( iPtBin, 100.*(extrapResp_RecoRel[iPtBin] - intrResp_RecoRel[iPtBin])/intrResp_RecoRel[iPtBin];
  //pullResp_err_RecoRel[iPtBin] = 100.*sqrt( extrapResp_err_RecoRel[iPtBin]*extrapResp_err_RecoRel[iPtBin] + intrResp_RecoRel[iPtBin]*intrResp_RecoRel[iPtBin] );

  //pullReso_RecoRel[iPtBin] = 100.*(extrapReso_RecoRel[iPtBin] - intrReso_RecoRel[iPtBin])/intrReso_RecoRel[iPtBin];
  //pullReso_err_RecoRel[iPtBin] = 100.*sqrt( extrapReso_err_RecoRel[iPtBin]*extrapReso_err_RecoRel[iPtBin] + intrReso_RecoRel[iPtBin]*intrReso_RecoRel[iPtBin] );

  //imbalance_RecoRel[iPtBin] *= 100.; //in percent
  //imbalance_err_RecoRel[iPtBin] *= 100.;

  } //for iPtBin


  std::string graphFileName = get_outputSuffix();

  TFile* graphFile = new TFile( graphFileName.c_str(), "recreate" );

} //drawExtrap






void DrawExtrap::getXPoints( TFile * file, const char* xHistoName, Int_t nPoints, Float_t* x, Float_t* x_err) const {


  for( int i=0; i<nPoints; ++i ) {
    char fullName[100];
    sprintf( fullName, "%s_%d", xHistoName,  i);
    TH1D* h1_pt2ndJetMean = (TH1D*)file->Get( fullName );
    x[i] = h1_pt2ndJetMean->GetMean();
    x_err[i] =  h1_pt2ndJetMean->GetRMS()/sqrt((Float_t)h1_pt2ndJetMean->GetEntries());
  }


} //getxpoints




void DrawExtrap::getYPoints( TFile * file, const char* yHistoName, Int_t nPoints, Float_t* y_resp, Float_t* y_resp_err,  Float_t* y_reso, Float_t* y_reso_err) const {

  for( int i=0; i<nPoints; ++i ) {

    char fullName[100];
    sprintf( fullName, "%s_%d", yHistoName, i);
//std::cout << fullName << " i: " << i << std::endl;

    TH1D* h1_r = (TH1D*)file->Get(fullName);

    //this ugly fix saves from empty relative pt binning plots
    if( h1_r->GetEntries()==0 ) {
      y_resp[i]=-1.;
      y_resp_err[i]=0.;
      y_reso[i]=-1.;
      y_reso_err[i]=0.;
      continue;
    }


    Float_t mean = h1_r->GetMean(); 
    Float_t rms = h1_r->GetRMS();
    Float_t mean_err = h1_r->GetMeanError();
    //Float_t rms_err = h1_r->GetRMSError();
    //Float_t rms_err = (LUMI>0) ? rms/sqrt((Float_t)LUMI*h1_r->Integral()*(Float_t)INT_PERC/100.) : h1_r->GetRMSError();
    Float_t rms_err = h1_r->GetRMSError();

    //Float_t mean, mean_err, rms, rms_err;
    //fitTools::getProjectionMeanAndRMS(h1_r, mean, mean_err, rms, rms_err, 1., 0.9);

    TF1* gaussian = new TF1("gaussian", "gaus");
    gaussian->SetLineColor(kRed);
    TH1D* newhisto = new TH1D(*h1_r);
    fitTools::fitProjection_sameArea(h1_r, gaussian, newhisto, (Float_t)(INTPERC_/100.), "RQ");
    //fitTools::fitProjection(h1_r, gaussian, 2., "RQ");

    Float_t mu = gaussian->GetParameter(1);
    Float_t sigma = gaussian->GetParameter(2);
    //TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
    //fitTools::fitProjection(h1_r, gaussian_chi, 1.5, "RQN");
    Float_t mu_err = gaussian->GetParError(1);
    //Float_t sigma_err = gaussian->GetParError(2);
    //Float_t sigma_err = (LUMI>0) ? sigma/sqrt((Float_t)LUMI*h1_r->Integral()*(Float_t)INT_PERC_/100.) : h1_r->GetRMSError();
    Float_t sigma_err = h1_r->GetRMSError();

    if( FIT_RMS_ == "FIT" ) {

      y_resp[i] = mu;
      y_resp_err[i] = mu_err;
    
      y_reso[i] = sigma/mu;
      y_reso_err[i] = sqrt( sigma_err*sigma_err/(mu*mu) + sigma*sigma*mu_err*mu_err/(mu*mu*mu*mu) ); 
//std::cout << "resp: " << y_resp[i] << " err: " << y_resp_err[i] << std::endl;

    } else if( FIT_RMS_ == "RMS" ) {

      y_resp[i] = mean;
      y_resp_err[i] = mean_err;
      
      y_reso[i] = rms/mean;
      y_reso_err[i] = sqrt( rms_err*rms_err/(mean*mean) + rms*rms*mean_err*mean_err/(mean*mean*mean*mean) ); 

    } else if( FIT_RMS_ == "MIX" ) {

      y_resp[i] = mean;
      y_resp_err[i] = mean_err;
      
      y_reso[i] = sigma/mean;
      y_reso_err[i] = sqrt( sigma_err*sigma_err/(mean*mean) + sigma*sigma*mean_err*mean_err/(mean*mean*mean*mean) ); 

    } else {

      std::cout << "WARNING!! FIT_RMS type '" << FIT_RMS_ << "' currently not supported. Exiting." << std::endl;
      exit(66);
      
    }

  } //for

} //getYPoints

