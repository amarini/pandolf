#include <TH2F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>


void fillProfile(TH1F* h1_resolution_FIT, TH1F* h1_resolution_RMS, TH2D* h2, char name[40]="");
void fillPositionResolution(TH1F* h1_sigmaEta, TH1F* h1_sigmaPhi, TH2D* h2_deltaEta, TH2D* h2_deltaPhi);
void fitProjection(TH1D* h1_projection, TF1* gaussian);



void Rx() {

  bool DEBUG_ = false;

  TChain* tree = new TChain("jetTree");

  int n_files = 49;
  if( DEBUG_ ) n_files = 1;

  for(int i=1; i<(n_files+1); ++i) {

    char treeName[70];
    sprintf(treeName, "2ndLevelOutputFiles/2ndLevelOutputFile_%d.root/jetTree", i);
    tree->Add(treeName);

  } //for



  TFile* outfile = new TFile("RxOutputFile.root", "RECREATE");
  outfile->cd();


  std::string algoType = "PFItCone5";

  int nReco;
  std::string nReco_name = "nJetReco_"+algoType;
  tree->SetBranchAddress(nReco_name.c_str(), &nReco);
  Float_t ptReco[40];
  std::string ptReco_name = "ptJetReco_"+algoType;
  tree->SetBranchAddress(ptReco_name.c_str(), ptReco);
  Float_t etaReco[40];
  std::string etaReco_name = "etaJetReco_"+algoType;
  tree->SetBranchAddress(etaReco_name.c_str(), etaReco);
  Float_t phiReco[40];
  std::string phiReco_name = "phiJetReco_"+algoType;
  tree->SetBranchAddress(phiReco_name.c_str(), phiReco);

  //PF only variables:
  Float_t eReco[40];
  std::string eReco_name = "eJetReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eReco_name.c_str(), eReco);
  Float_t ptTracksReco[40];
  std::string ptTracksReco_name = "ptTracksReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ptTracksReco_name.c_str(), ptTracksReco);
  Float_t ptPhotonsReco[40];
  std::string ptPhotonsReco_name = "ptPhotonsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ptPhotonsReco_name.c_str(), ptPhotonsReco);
  Float_t eTracksReco[40];
  std::string eTracksReco_name = "eTracksReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eTracksReco_name.c_str(), eTracksReco);
  Float_t ePhotonsReco[40];
  std::string ePhotonsReco_name = "ePhotonsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ePhotonsReco_name.c_str(), ePhotonsReco);
  Float_t eNeutralHadronsReco[40];
  std::string eNeutralHadronsReco_name = "eNeutralHadronsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eNeutralHadronsReco_name.c_str(), eNeutralHadronsReco);
  Float_t eMuonsReco[40];
  std::string eMuonsReco_name = "eMuonsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eMuonsReco_name.c_str(), eMuonsReco);
  Float_t eHFhadronsReco[40];
  std::string eHFhadronsReco_name = "eHFhadronsReco_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eHFhadronsReco_name.c_str(), eHFhadronsReco);

  Float_t eGen[40];
  std::string eGen_name = "eJetGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eGen_name.c_str(), eGen);
  Float_t eTracksGen[40];
  std::string eTracksGen_name = "eTracksGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eTracksGen_name.c_str(), eTracksGen);
  Float_t ePhotonsGen[40];
  std::string ePhotonsGen_name = "ePhotonsGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(ePhotonsGen_name.c_str(), ePhotonsGen);
  Float_t eNeutralHadronsGen[40];
  std::string eNeutralHadronsGen_name = "eNeutralHadronsGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eNeutralHadronsGen_name.c_str(), eNeutralHadronsGen);
  Float_t eMuonsGen[40];
  std::string eMuonsGen_name = "eMuonsGen_"+algoType;
  if( algoType=="PFItCone5")
    tree->SetBranchAddress(eMuonsGen_name.c_str(), eMuonsGen);

  Float_t ptGen[40];
  std::string ptGen_name = "ptJetGen_"+algoType;
  tree->SetBranchAddress(ptGen_name.c_str(), ptGen);
  Float_t etaGen[40];
  std::string etaGen_name = "etaJetGen_"+algoType;
  tree->SetBranchAddress(etaGen_name.c_str(), etaGen);
  Float_t phiGen[40];
  std::string phiGen_name = "phiJetGen_"+algoType;
  tree->SetBranchAddress(phiGen_name.c_str(), phiGen);


  std::string histoName;

  //for pT spectrums, variable binning to avoid huge error bars at high pT:
  int nBins = 21;
  Double_t Lower[nBins];

  for(int iLower=0; iLower<11; ++iLower)
    Lower[iLower] = iLower*10.;
  
  Lower[11] = 120.;
  Lower[12] = 140.;
  Lower[13] = 160.;
  Lower[14] = 180.;
  Lower[15] = 200.;
  Lower[16] = 250.;
  Lower[17] = 300.;
  Lower[18] = 400.;
  Lower[19] = 500.;
  Lower[20] = 600.;


  histoName = "E_resolution_"+algoType;
  TH1F* h1_E_resolution = new TH1F("E_resolution", "R_{ch} Resolution", 100, -1.1, 1.1);
  h1_E_resolution->SetXTitle("(R_{ch}^{RECO}-R_{ch}^{GEN})/R_{ch}^{GEN} [%]");
  h1_E_resolution->SetYTitle("entries");

  histoName = "Ech_resolution_"+algoType;
  TH1F* h1_Ech_resolution = new TH1F("Ech_resolution", "R_{ch} Resolution", 100, -1.1, 1.1);
  h1_Ech_resolution->SetXTitle("(R_{ch}^{RECO}-R_{ch}^{GEN})/R_{ch}^{GEN} [%]");
  h1_Ech_resolution->SetYTitle("entries");

  histoName = "Egamma_resolution_"+algoType;
  TH1F* h1_Egamma_resolution = new TH1F("Egamma_resolution", "R_{#gamma} Resolution", 100, -1.1, 1.1);
  h1_Egamma_resolution->SetXTitle("(R_{#gamma}^{RECO}-R_{#gamma}^{GEN})/R_{#gamma}^{GEN} [%]");
  h1_Egamma_resolution->SetYTitle("entries");

  histoName = "Enh_resolution_"+algoType;
  TH1F* h1_Enh_resolution = new TH1F("Enh_resolution", "R_{nh} Resolution", 100, -1.1, 1.1);
  h1_Enh_resolution->SetXTitle("(R_{nh}^{RECO}-R_{nh}^{GEN})/R_{nh}^{GEN} [%]");
  h1_Enh_resolution->SetYTitle("entries");

  histoName = "E_resolution_50_80_"+algoType;
  TH1F* h1_E_resolution_50_80 = new TH1F("E_resolution_50_80", "R_{ch} Resolution", 100, -1.1, 1.1);
  h1_E_resolution_50_80->SetXTitle("(R_{ch}^{RECO}-R_{ch}^{GEN})/R_{ch}^{GEN} [%]");
  h1_E_resolution_50_80->SetYTitle("entries");

  histoName = "Ech_resolution_50_80_"+algoType;
  TH1F* h1_Ech_resolution_50_80 = new TH1F("Ech_resolution_50_80", "R_{ch} Resolution", 100, -1.1, 1.1);
  h1_Ech_resolution_50_80->SetXTitle("(R_{ch}^{RECO}-R_{ch}^{GEN})/R_{ch}^{GEN} [%]");
  h1_Ech_resolution_50_80->SetYTitle("entries");

  histoName = "Egamma_resolution_50_80_"+algoType;
  TH1F* h1_Egamma_resolution_50_80 = new TH1F("Egamma_resolution_50_80", "R_{#gamma} Resolution", 100, -1.1, 1.1);
  h1_Egamma_resolution_50_80->SetXTitle("(R_{#gamma}^{RECO}-R_{#gamma}^{GEN})/R_{#gamma}^{GEN} [%]");
  h1_Egamma_resolution_50_80->SetYTitle("entries");

  histoName = "Enh_resolution_50_80_"+algoType;
  TH1F* h1_Enh_resolution_50_80 = new TH1F("Enh_resolution_50_80", "R_{nh} Resolution", 100, -1.1, 1.1);
  h1_Enh_resolution_50_80->SetXTitle("(R_{nh}^{RECO}-R_{nh}^{GEN})/R_{nh}^{GEN} [%]");
  h1_Enh_resolution_50_80->SetYTitle("entries");

  histoName = "deltaEch_vs_pt_"+algoType;
  TH2D* h2_deltaEch_vs_pt = new TH2D(histoName.c_str(), "R_{ch}^{RECO}-R_{ch}^{GEN} (|#eta_{GEN}|<2.6)", nBins-2, Lower, 100, -1.1, 1.1);
  h2_deltaEch_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEch_vs_pt->SetYTitle("#Delta R_{ch}");

  histoName = "deltaEgamma_vs_pt_"+algoType;
  TH2D* h2_deltaEgamma_vs_pt = new TH2D(histoName.c_str(), "R_{#gamma}^{RECO}-R_{#gamma}^{GEN} (|#eta_{GEN}|<2.6)", nBins-2, Lower, 100, -1.1, 1.1);
  h2_deltaEgamma_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEgamma_vs_pt->SetYTitle("#Delta R_{#gamma}");

  histoName = "deltaEnh_vs_pt_"+algoType;
  TH2D* h2_deltaEnh_vs_pt = new TH2D(histoName.c_str(), "R_{nh}^{RECO}-R_{nh}^{GEN} (|#eta_{GEN}|<2.6)", nBins-2, Lower, 100, -1.1, 1.1);
  h2_deltaEnh_vs_pt->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h2_deltaEnh_vs_pt->SetYTitle("#Delta R_{nh}");

  histoName = "deltaEch_vs_eta_"+algoType;
  TH2D* h2_deltaEch_vs_eta = new TH2D(histoName.c_str(), "R_{ch}^{RECO}-R_{ch}^{GEN}", 100, -3., 3., 100, -1.1, 1.1);
  h2_deltaEch_vs_eta->SetXTitle("#eta^{GEN}");
  h2_deltaEch_vs_eta->SetYTitle("#Delta R_{ch}");

  histoName = "deltaEgamma_vs_eta_"+algoType;
  TH2D* h2_deltaEgamma_vs_eta = new TH2D(histoName.c_str(), "R_{#gamma}^{RECO}-R_{#gamma}^{GEN}", 100, -3., 3., 100, -1.1, 1.1);
  h2_deltaEgamma_vs_eta->SetXTitle("#eta^{GEN}");
  h2_deltaEgamma_vs_eta->SetYTitle("#Delta R_{#gamma}");

  histoName = "deltaEnh_vs_eta_"+algoType;
  TH2D* h2_deltaEnh_vs_eta = new TH2D(histoName.c_str(), "R_{nh}^{RECO}-R_{nh}^{GEN}", 100, -3., 3., 100, -1.1, 1.1);
  h2_deltaEnh_vs_eta->SetXTitle("#eta^{GEN}");
  h2_deltaEnh_vs_eta->SetYTitle("#Delta R_{nh}");



  int nEntries = tree->GetEntries();

  //for(int iEntry=0; iEntry<1000000; ++iEntry) {
  for(int iEntry=0; iEntry<nEntries; ++iEntry) {

    if( (iEntry % 100000)==0 ) std::cout << "Entry: " << iEntry << " /" << nEntries << std::endl;

    tree->GetEntry(iEntry);

    for(int iRecoJet = 0; iRecoJet<nReco; ++iRecoJet) {

      if( iRecoJet >= 40 ) continue;
      //if( ptReco[iRecoJet]<15. ) continue;
      if( ptGen[iRecoJet]<15. ) continue;
      //if( ptGen[iRecoJet]<20. ) continue;

      Float_t iEtaReco = etaReco[iRecoJet];
      Float_t iPhiReco = phiReco[iRecoJet];
      Float_t iEtaGen = etaGen[iRecoJet];
      Float_t iPhiGen = phiGen[iRecoJet];
      Float_t iThetaReco = 2*atan(exp(-iEtaReco));

      Float_t iPtReco = ptReco[iRecoJet];
      Float_t iEReco = eReco[iRecoJet];
      Float_t iEtReco = iEReco*sin(iThetaReco);
      Float_t iPtGen = ptGen[iRecoJet];
      Float_t iEGen = eGen[iRecoJet];
      Float_t response = iPtReco/iPtGen;

      //physically important is resolution on energy, not R_x:
      Float_t Ech = eTracksReco[iRecoJet];
      Float_t Egamma = ePhotonsReco[iRecoJet];
      Float_t Enh = eNeutralHadronsReco[iRecoJet];
      Float_t Emu = eMuonsReco[iRecoJet];
      Float_t Rhf = eHFhadronsReco[iRecoJet];

      Float_t EchGen = eTracksGen[iRecoJet];
      Float_t EgammaGen = ePhotonsGen[iRecoJet];
      Float_t EnhGen = eNeutralHadronsGen[iRecoJet];
      Float_t EmuGen = eMuonsGen[iRecoJet];

      Float_t deltaE= iEReco-iEGen;
      Float_t deltaEch = Ech-EchGen;
      Float_t deltaEgamma = Egamma-EgammaGen;
      Float_t deltaEnh = Enh-EnhGen;

      Float_t deltaEta = iEtaReco-iEtaGen;
      Float_t deltaPhi = iPhiReco-iPhiGen;
      Float_t pi = 3.14159;
      if(deltaPhi < -pi) deltaPhi += 2.*pi;
      if(deltaPhi > pi) deltaPhi -= 2.*pi;

      Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );


      if( deltaR <0.3 ) { //matched
    
        deltaEta = fabs(iEtaReco) - fabs(iEtaGen);
    
        if( fabs(iEtaGen)<2.6 ) { //tracker region

          h1_E_resolution->Fill( deltaE/iEGen );
          h1_Ech_resolution->Fill( deltaEch/EchGen );
          h1_Egamma_resolution->Fill( deltaEgamma/EgammaGen );
          h1_Enh_resolution->Fill( deltaEnh/EnhGen );

          if( (iPtGen>50.) && (iPtGen<80.) ) {
            h1_E_resolution_50_80->Fill( deltaE/iEGen );
            h1_Ech_resolution_50_80->Fill( deltaEch/EchGen );
            h1_Egamma_resolution_50_80->Fill( deltaEgamma/EgammaGen );
            h1_Enh_resolution_50_80->Fill( deltaEnh/EnhGen );
          }

          h2_deltaEch_vs_pt->Fill(iPtGen, deltaEch);
          h2_deltaEgamma_vs_pt->Fill(iPtGen, deltaEgamma);
          h2_deltaEnh_vs_pt->Fill(iPtGen, deltaEnh);

          h2_deltaEch_vs_eta->Fill(iEtaGen, deltaEch);
          h2_deltaEgamma_vs_eta->Fill(iEtaGen, deltaEgamma);
          h2_deltaEnh_vs_eta->Fill(iEtaGen, deltaEnh);

        }

      } //for recoJets

    } //if matched

  } //for entries

 
  histoName = "Ech_resolution_vs_pt_FIT_"+algoType;
  TH1F* h1_Ech_resolution_vs_pt_FIT= new TH1F(histoName.c_str(), "R_{ch} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", nBins-2, Lower);
  h1_Ech_resolution_vs_pt_FIT->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Ech_resolution_vs_pt_FIT->SetYTitle("#sigma(R_{ch}^{RECO}-R_{ch}^{GEN}) [%]");
  h1_Ech_resolution_vs_pt_FIT->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Egamma_resolution_vs_pt_FIT_"+algoType;
  TH1F* h1_Egamma_resolution_vs_pt_FIT= new TH1F(histoName.c_str(), "R_{#gamma} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", nBins-2, Lower);
  h1_Egamma_resolution_vs_pt_FIT->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Egamma_resolution_vs_pt_FIT->SetYTitle("#sigma(R_{#gamma}^{RECO}-R_{#gamma}^{GEN}) [%]");
  h1_Egamma_resolution_vs_pt_FIT->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Enh_resolution_vs_pt_FIT_"+algoType;
  TH1F* h1_Enh_resolution_vs_pt_FIT= new TH1F(histoName.c_str(), "R_{nh} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", nBins-2, Lower);
  h1_Enh_resolution_vs_pt_FIT->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Enh_resolution_vs_pt_FIT->SetYTitle("#sigma(R_{nh}^{RECO}-R_{nh}^{GEN}) [%]");
  h1_Enh_resolution_vs_pt_FIT->GetYaxis()->SetRangeUser(0., 1.);

  histoName = "Ech_resolution_vs_pt_RMS_"+algoType;
  TH1F* h1_Ech_resolution_vs_pt_RMS= new TH1F(histoName.c_str(), "R_{ch} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", nBins-2, Lower);
  h1_Ech_resolution_vs_pt_RMS->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Ech_resolution_vs_pt_RMS->SetYTitle("#sigma(R_{ch}^{RECO}-R_{ch}^{GEN}) [%]");
  h1_Ech_resolution_vs_pt_RMS->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Egamma_resolution_vs_pt_RMS_"+algoType;
  TH1F* h1_Egamma_resolution_vs_pt_RMS= new TH1F(histoName.c_str(), "R_{#gamma} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", nBins-2, Lower);
  h1_Egamma_resolution_vs_pt_RMS->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Egamma_resolution_vs_pt_RMS->SetYTitle("#sigma(R_{#gamma}^{RECO}-R_{#gamma}^{GEN}) [%]");
  h1_Egamma_resolution_vs_pt_RMS->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Enh_resolution_vs_pt_RMS_"+algoType;
  TH1F* h1_Enh_resolution_vs_pt_RMS= new TH1F(histoName.c_str(), "R_{nh} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", nBins-2, Lower);
  h1_Enh_resolution_vs_pt_RMS->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Enh_resolution_vs_pt_RMS->SetYTitle("#sigma(R_{nh}^{RECO}-R_{nh}^{GEN}) [%]");
  h1_Enh_resolution_vs_pt_RMS->GetYaxis()->SetRangeUser(0., 1.);

  histoName = "Ech_resolution_vs_eta_FIT_"+algoType;
  TH1F* h1_Ech_resolution_vs_eta_FIT= new TH1F(histoName.c_str(), "R_{ch} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", 100, -3., 3.);
  h1_Ech_resolution_vs_eta_FIT->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Ech_resolution_vs_eta_FIT->SetYTitle("#sigma(R_{ch}^{RECO}-R_{ch}^{GEN}) [%]");
  h1_Ech_resolution_vs_eta_FIT->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Egamma_resolution_vs_eta_FIT_"+algoType;
  TH1F* h1_Egamma_resolution_vs_eta_FIT= new TH1F(histoName.c_str(), "R_{#gamma} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", 100, -3., 3.);
  h1_Egamma_resolution_vs_eta_FIT->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Egamma_resolution_vs_eta_FIT->SetYTitle("#sigma(R_{#gamma}^{RECO}-R_{#gamma}^{GEN}) [%]");
  h1_Egamma_resolution_vs_eta_FIT->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Enh_resolution_vs_eta_FIT_"+algoType;
  TH1F* h1_Enh_resolution_vs_eta_FIT= new TH1F(histoName.c_str(), "R_{nh} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", 100, -3., 3.);
  h1_Enh_resolution_vs_eta_FIT->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Enh_resolution_vs_eta_FIT->SetYTitle("#sigma(R_{nh}^{RECO}-R_{nh}^{GEN}) [%]");
  h1_Enh_resolution_vs_eta_FIT->GetYaxis()->SetRangeUser(0., 1.);

  histoName = "Ech_resolution_vs_eta_RMS_"+algoType;
  TH1F* h1_Ech_resolution_vs_eta_RMS= new TH1F(histoName.c_str(), "R_{ch} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", 100, -3., 3.);
  h1_Ech_resolution_vs_eta_RMS->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Ech_resolution_vs_eta_RMS->SetYTitle("#sigma(R_{ch}^{RECO}-R_{ch}^{GEN}) [%]");
  h1_Ech_resolution_vs_eta_RMS->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Egamma_resolution_vs_eta_RMS_"+algoType;
  TH1F* h1_Egamma_resolution_vs_eta_RMS= new TH1F(histoName.c_str(), "R_{#gamma} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", 100, -3., 3.);
  h1_Egamma_resolution_vs_eta_RMS->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Egamma_resolution_vs_eta_RMS->SetYTitle("#sigma(R_{#gamma}^{RECO}-R_{#gamma}^{GEN}) [%]");
  h1_Egamma_resolution_vs_eta_RMS->GetYaxis()->SetRangeUser(0., 1.);
  
  histoName = "Enh_resolution_vs_eta_RMS_"+algoType;
  TH1F* h1_Enh_resolution_vs_eta_RMS= new TH1F(histoName.c_str(), "R_{nh} Resolution vs. True Jet p_{T} (0.0<|#eta_{GEN}|<2.6)", 100, -3., 3.);
  h1_Enh_resolution_vs_eta_RMS->SetXTitle("p_{T}^{GEN} [GeV/c]");
  h1_Enh_resolution_vs_eta_RMS->SetYTitle("#sigma(R_{nh}^{RECO}-R_{nh}^{GEN}) [%]");
  h1_Enh_resolution_vs_eta_RMS->GetYaxis()->SetRangeUser(0., 1.);

  fillProfile(h1_Ech_resolution_vs_pt_FIT, h1_Ech_resolution_vs_pt_RMS, h2_deltaEch_vs_pt);
  fillProfile(h1_Egamma_resolution_vs_pt_FIT, h1_Egamma_resolution_vs_pt_RMS, h2_deltaEgamma_vs_pt);
  fillProfile(h1_Enh_resolution_vs_pt_FIT, h1_Enh_resolution_vs_pt_RMS, h2_deltaEnh_vs_pt);

  fillProfile(h1_Ech_resolution_vs_eta_FIT, h1_Ech_resolution_vs_eta_RMS, h2_deltaEch_vs_eta);
  fillProfile(h1_Egamma_resolution_vs_eta_FIT, h1_Egamma_resolution_vs_eta_RMS, h2_deltaEgamma_vs_eta);
  fillProfile(h1_Enh_resolution_vs_eta_FIT, h1_Enh_resolution_vs_eta_RMS, h2_deltaEnh_vs_eta);
 

  delete h2_deltaEch_vs_pt;
  h2_deltaEch_vs_pt = 0;
  delete h2_deltaEgamma_vs_pt;
  h2_deltaEgamma_vs_pt= 0;
  delete h2_deltaEnh_vs_pt;
  h2_deltaEnh_vs_pt= 0;

  delete h2_deltaEch_vs_eta;
  h2_deltaEch_vs_eta = 0;
  delete h2_deltaEgamma_vs_eta;
  h2_deltaEgamma_vs_eta= 0;
  delete h2_deltaEnh_vs_eta;
  h2_deltaEnh_vs_eta= 0;

  outfile->cd();

  h1_E_resolution->Draw();
  h1_E_resolution_50_80->Draw();

  h1_Ech_resolution->Draw();
  h1_Egamma_resolution->Draw();
  h1_Enh_resolution->Draw();

  h1_Ech_resolution_50_80->Draw();
  h1_Egamma_resolution_50_80->Draw();
  h1_Enh_resolution_50_80->Draw();

  h1_Ech_resolution_vs_pt_FIT->Draw();
  h1_Egamma_resolution_vs_pt_FIT->Draw();
  h1_Enh_resolution_vs_pt_FIT->Draw();

  h1_Ech_resolution_vs_pt_RMS->Draw();
  h1_Egamma_resolution_vs_pt_RMS->Draw();
  h1_Enh_resolution_vs_pt_RMS->Draw();

  h1_Ech_resolution_vs_eta_FIT->Draw();
  h1_Egamma_resolution_vs_eta_FIT->Draw();
  h1_Enh_resolution_vs_eta_FIT->Draw();

  h1_Ech_resolution_vs_eta_RMS->Draw();
  h1_Egamma_resolution_vs_eta_RMS->Draw();
  h1_Enh_resolution_vs_eta_RMS->Draw();

  outfile->Write();
  outfile->Close();
  delete outfile;
  outfile = 0;


} //fillFile


void fillProfile(TH1F* h1_resolution_FIT, TH1F* h1_resolution_RMS, TH2D* h2, char name[40]) {

  for(int iBin=1; iBin<(h2->GetNbinsX()+1); ++iBin) {

    TH1D* h1_projection = h2->ProjectionY("projectiony", iBin, iBin);

    TF1* gaussian = new TF1("gaussian", "gaus");

    fitProjection(h1_projection, gaussian);

    if( name!="" ) {
        TCanvas* c = new TCanvas("c1", "c1", 800, 600);
        h1_projection->Draw();
        char canvasName[50];
        sprintf(canvasName, "Projections/%s_fit%d.eps", name, iBin);
        c->SaveAs(canvasName);
    }

    Float_t mu = gaussian->GetParameter(1);
    Float_t mu_err = gaussian->GetParError(1);

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


void fitProjection(TH1D* h1_projection, TF1* gaussian) {

   Float_t histMean = h1_projection->GetMean();
   Float_t histRMS = h1_projection->GetRMS();
 
   gaussian->SetParameter(0, h1_projection->GetMaximum());
   gaussian->SetParameter(1, histMean);
   gaussian->SetParameter(2, histRMS);
   Float_t lowerBound = histMean - 1.5*histRMS;
   Float_t upperBound = histMean + 1.5*histRMS;
   gaussian->SetRange( lowerBound, upperBound );
 
   h1_projection->Fit(gaussian, "R");

}
