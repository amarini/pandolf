#include <cmath>

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"

using namespace std;



int main() {

  TFile* file = TFile::Open("output_SAVE.root");
  TTree* tree = (TTree*)file->Get("myanalysis/pippo");

  float rhoPF;
  tree->SetBranchAddress("rhoPF", &rhoPF);
  int nJet;
  tree->SetBranchAddress("nJet", &nJet);
  
  float ptJet[10];
  tree->SetBranchAddress("ptJet", ptJet);
  float etaJet[10];
  tree->SetBranchAddress("etaJet", etaJet);
  float phiJet[10];
  tree->SetBranchAddress("phiJet", phiJet);
  float axis2_QCJet[10];
  tree->SetBranchAddress("axis2_QCJet", axis2_QCJet);
  float ptD_QCJet[10];
  tree->SetBranchAddress("ptD_QCJet", ptD_QCJet);
  int nPFCand_QC_ptCutJet[10];
  tree->SetBranchAddress("nPFCand_QC_ptCutJet", nPFCand_QC_ptCutJet);

  float qglJet[10];
  tree->SetBranchAddress("qglJet", qglJet);

  int pdgIdPartJet[10];
  tree->SetBranchAddress("pdgIdPartJet", pdgIdPartJet);


  QGLikelihoodCalculator* qglc = new QGLikelihoodCalculator("/afs/cern.ch/work/p/pandolf/public/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root");

  TH1D* h1_qglOLD_quark_pt3050 = new TH1D("qglOLD_quark_pt3050", "", 100, 0., 1.0001);
  TH1D* h1_qglNEW_quark_pt3050 = new TH1D("qglNEW_quark_pt3050", "", 100, 0., 1.0001);
  TH1D* h1_qglOLD_gluon_pt3050 = new TH1D("qglOLD_gluon_pt3050", "", 100, 0., 1.0001);
  TH1D* h1_qglNEW_gluon_pt3050 = new TH1D("qglNEW_gluon_pt3050", "", 100, 0., 1.0001);

  TH1D* h1_qglOLD_quark_pt80120 = new TH1D("qglOLD_quark_pt80120", "", 100, 0., 1.0001);
  TH1D* h1_qglNEW_quark_pt80120 = new TH1D("qglNEW_quark_pt80120", "", 100, 0., 1.0001);
  TH1D* h1_qglOLD_gluon_pt80120 = new TH1D("qglOLD_gluon_pt80120", "", 100, 0., 1.0001);
  TH1D* h1_qglNEW_gluon_pt80120 = new TH1D("qglNEW_gluon_pt80120", "", 100, 0., 1.0001);

  TH1D* h1_qglOLD_quark_pt200300 = new TH1D("qglOLD_quark_pt200300", "", 100, 0., 1.0001);
  TH1D* h1_qglNEW_quark_pt200300 = new TH1D("qglNEW_quark_pt200300", "", 100, 0., 1.0001);
  TH1D* h1_qglOLD_gluon_pt200300 = new TH1D("qglOLD_gluon_pt200300", "", 100, 0., 1.0001);
  TH1D* h1_qglNEW_gluon_pt200300 = new TH1D("qglNEW_gluon_pt200300", "", 100, 0., 1.0001);



  int nentries = tree->GetEntries();

  for( unsigned int iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);


    if( nJet==0 ) continue;
    if( ptJet[0]<20. ) continue;

    float pt = ptJet[0];
    float eta = etaJet[0];
    int pdgid = pdgIdPartJet[0];
    float qglNEW = qglJet[0];

    float qglOLD = qglc->computeQGLikelihood2012( pt, eta, rhoPF, nPFCand_QC_ptCutJet[0], ptD_QCJet[0], axis2_QCJet[0]);


    if( fabs(eta)<2. ) {

      if( pt>30. && pt<50. ) {

        if( abs(pdgid)<4 && abs(pdgid)>0 ) {
          h1_qglOLD_quark_pt3050->Fill(qglOLD);
          h1_qglNEW_quark_pt3050->Fill(qglNEW);
        } else if( pdgid==21 ) {
          h1_qglOLD_gluon_pt3050->Fill(qglOLD);
          h1_qglNEW_gluon_pt3050->Fill(qglNEW);
        }

      } else if( pt>80 && pt<120. ) {

        if( abs(pdgid)<4 && abs(pdgid)>0 ) {
          h1_qglOLD_quark_pt80120->Fill(qglOLD);
          h1_qglNEW_quark_pt80120->Fill(qglNEW);
        } else if( pdgid==21 ) {
          h1_qglOLD_gluon_pt80120->Fill(qglOLD);
          h1_qglNEW_gluon_pt80120->Fill(qglNEW);
        }

      } else if( pt>200. && pt<300.) {

        if( abs(pdgid)<4 && abs(pdgid)>0 ) {
          h1_qglOLD_quark_pt200300->Fill(qglOLD);
          h1_qglNEW_quark_pt200300->Fill(qglNEW);
        } else if( pdgid==21 ) {
          h1_qglOLD_gluon_pt200300->Fill(qglOLD);
          h1_qglNEW_gluon_pt200300->Fill(qglNEW);
        }

      }


    } else if( fabs(eta)>3. ) {



    }


  } //entries



  TFile* outfile = TFile::Open("provaQG.root", "RECREATE");
  outfile->cd();

  h1_qglOLD_quark_pt3050->Write();
  h1_qglNEW_quark_pt3050->Write();
  h1_qglOLD_gluon_pt3050->Write();
  h1_qglNEW_gluon_pt3050->Write();

  h1_qglOLD_quark_pt80120->Write();
  h1_qglNEW_quark_pt80120->Write();
  h1_qglOLD_gluon_pt80120->Write();
  h1_qglNEW_gluon_pt80120->Write();

  h1_qglOLD_quark_pt200300->Write();
  h1_qglNEW_quark_pt200300->Write();
  h1_qglOLD_gluon_pt200300->Write();
  h1_qglNEW_gluon_pt200300->Write();

  outfile->Close();

  return 0;

}
