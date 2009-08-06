
void fitProjection(TH1D* h1_projection, TF1* gaussian) {

     
   Float_t histMean = h1_projection->GetMean();
   Float_t histRMS = h1_projection->GetRMS();
 
   gaussian->SetParameter(0, h1_projection->GetMaximum());
   gaussian->SetParameter(1, histMean);
   gaussian->SetParameter(2, histRMS);
   if( histRMS==0. ) return;

   gaussian->SetParLimits(1, 0., 2.*histMean);

   Float_t nSigma = 5.;

   Float_t lowerBound = histMean - nSigma*histRMS;
   Float_t upperBound = histMean + nSigma*histRMS;
   gaussian->SetRange( lowerBound, upperBound );
 
   h1_projection->Fit(gaussian, "R");

   int n_iter = 3;

   for(int i=0; i<n_iter; ++i) {

     Float_t lowerBound = gaussian->GetParameter(1) - nSigma*gaussian->GetParameter(2);
     Float_t upperBound = gaussian->GetParameter(1) + nSigma*gaussian->GetParameter(2);
     gaussian->SetRange( lowerBound, upperBound );

     h1_projection->Fit(gaussian, "R");

   }

}



void fillProfile(TH1F* h1_response_FIT, TH1F* h1_resolution_FIT, TH1F* h1_response_MEAN, TH1F* h1_resolution_RMS, TH2D* h2, std::string name="") {

  std::string fileName = "Projections/"+name+".root";
  TFile* projectionFile;
  if( name!= "" ) {
    projectionFile = TFile::Open(fileName.c_str(), "RECREATE");
    projectionFile->cd();
  }

  for(int iBin=1; iBin<(h2->GetNbinsX()+1); ++iBin) {

    char histName[50];
    sprintf(histName, "projection_%d",iBin);
    TH1D* h1_projection = h2->ProjectionY(histName, iBin, iBin);

    TF1* gaussian = new TF1("gaussian", "gaus");

    fitProjection(h1_projection, gaussian);

    if( name!="" ) {
        h1_projection->Write();
    }

    Float_t mu = gaussian->GetParameter(1);
    Float_t mu_err = gaussian->GetParError(1);
    h1_response_FIT->SetBinContent(iBin, mu);
    h1_response_FIT->SetBinError(iBin, mu_err);

    Float_t sigma = gaussian->GetParameter(2);
    Float_t resolution = (mu!=0.) ? sigma/mu : -1.;
    h1_resolution_FIT->SetBinContent(iBin, resolution);

    Float_t sigma_err = gaussian->GetParError(2);
    Float_t res_err = (mu!=0.) ? sqrt( sigma_err*sigma_err/(mu*mu) + mu_err*mu_err*sigma*sigma/(mu*mu*mu*mu) ) : 0.;
    h1_resolution_FIT->SetBinError(iBin, res_err);

    delete gaussian;
    gaussian = 0;

    Float_t n = h1_projection->GetEntries();
    Float_t mean = h1_projection->GetMean();
    Float_t mean_err = (n!=0) ? h1_projection->GetRMS()/sqrt(n) : 0.;
    h1_response_MEAN->SetBinContent(iBin, mean);
    h1_response_MEAN->SetBinError(iBin, mean_err);

    Float_t rms = h1_projection->GetRMS();
    Float_t rms_err = (n!=0) ? h1_projection->GetRMS()/sqrt(n) : 0.;
    resolution = (mean!=0.) ? rms/mean : -1.;
    res_err = (mean!=0.) ? sqrt( rms_err*rms_err/(mean*mean) + mean_err*mean_err*rms*rms/(mean*mean*mean*mean) ) : 0.;
    if( resolution != 0. ) {
      h1_resolution_RMS->SetBinContent(iBin, resolution);
      h1_resolution_RMS->SetBinError(iBin, res_err);
    }

    h1_projection = 0;

  } //for bins

  if(name!="") {
    projectionFile->Write();
    projectionFile->Close();
    delete projectionFile;
  }
  projectionFile = 0;

} //fill profile


void fillPositionResolution(TH1F* h1_sigmaEta, TH1F* h1_sigmaPhi, TH2D* h2_deltaEta, TH2D* h2_deltaPhi) {

  for(int iBin=1; iBin<(h2_deltaEta->GetNbinsX()+1); ++iBin) {

    TH1D* h1_projection = h2_deltaEta->ProjectionY("projectiony", iBin, iBin);

    TF1* gaussian = new TF1("gaussian", "gaus");
  
    fitProjection(h1_projection, gaussian);

    h1_sigmaEta->SetBinContent(iBin, gaussian->GetParameter(2));
    h1_sigmaEta->SetBinError(iBin, gaussian->GetParError(2));
  
    delete gaussian;
    gaussian = 0;

    h1_projection = 0;

  } //for bins eta

  for(int iBin=1; iBin<(h2_deltaPhi->GetNbinsX()+1); ++iBin) {

    TH1D* h1_projection = h2_deltaPhi->ProjectionY("projectiony", iBin, iBin);
  
    TF1* gaussian = new TF1("gaussian", "gaus");
  
    fitProjection(h1_projection, gaussian);

    h1_sigmaPhi->SetBinContent(iBin, gaussian->GetParameter(2));
    h1_sigmaPhi->SetBinError(iBin, gaussian->GetParError(2));
  
    delete gaussian;
    gaussian = 0;

    h1_projection = 0;

  } //for bins

} //fill position resolution


