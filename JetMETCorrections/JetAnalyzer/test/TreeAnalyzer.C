#define TreeAnalyzer_cxx
#include "TreeAnalyzer.h"
#include <iostream>





void TreeAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L TreeAnalyzer.C
//      Root > TreeAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;



   Long64_t nentries;

   if( DEBUG_ ) nentries = 100000;
   else nentries = fChain->GetEntries();

   BookStuff();



   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;

if( DEBUG_VERBOSE_ ) std::cout << "entry n." << jentry << std::endl;

     if( (jentry%10000) == 0 ) std::cout << "Event #" << jentry  << " of " << nentries << std::endl;

     nJetReco_caloItCone5_ = 0;
     nJetReco_caloSISCone5_ = 0;
     nJetReco_calokt4_ = 0;
     nJetReco_PFItCone5_ = 0;


     // ------------  CALO ITERATIVE CONE 0.5 JETS

     for(unsigned int iRecoJet=0; iRecoJet<nJetCalo_itCone5; ++iRecoJet) {

       if( nJetReco_caloItCone5_>=40 ) continue;

       Float_t  pxJetReco  =  pxJetCalo_itCone5[iRecoJet];
       Float_t  pyJetReco  =  pyJetCalo_itCone5[iRecoJet];
       Float_t  pzJetReco  =  pzJetCalo_itCone5[iRecoJet];
       Float_t  phiJetReco = phiJetCalo_itCone5[iRecoJet];
       Float_t  etaJetReco = etaJetCalo_itCone5[iRecoJet];
       Float_t  eJetReco = eJetCalo_itCone5[iRecoJet];

       Float_t pJetReco = sqrt( pxJetReco*pxJetReco + pyJetReco*pyJetReco + pzJetReco*pzJetReco );
       Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
       Float_t ptJetReco = pJetReco*sin(thetaJetReco);

       if ( ptJetReco < 4. ) continue;


       Float_t deltaRmin = 999.;
       Float_t etaJetGen_found;
       Float_t phiJetGen_found;
       Float_t ptJetGen_found;
       Float_t eJetGen_found;

       for(unsigned int iGenJet=0; iGenJet<nJetGen_itCone5; ++iGenJet) {

         Float_t  pxJetGen  =  pxJetGen_itCone5[iGenJet];
         Float_t  pyJetGen  =  pyJetGen_itCone5[iGenJet];
         Float_t  pzJetGen  =  pzJetGen_itCone5[iGenJet];
         Float_t  phiJetGen = phiJetGen_itCone5[iGenJet];
         Float_t  etaJetGen = etaJetGen_itCone5[iGenJet];
         Float_t  eJetGen = eJetGen_itCone5[iGenJet];

         Float_t deltaEta = etaJetReco - etaJetGen;
         Float_t deltaPhi = phiJetReco - phiJetGen;

         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

         if( deltaR < deltaRmin ) {
           deltaRmin = deltaR;
           etaJetGen_found = etaJetGen;
           phiJetGen_found = phiJetGen;
           Float_t pJetGen = sqrt( pxJetGen*pxJetGen + pyJetGen*pyJetGen + pzJetGen*pzJetGen );
           Float_t thetaJetGen = 2.*atan(exp(-etaJetGen));
           ptJetGen_found = pJetGen*sin(thetaJetGen);
           eJetGen_found = eJetGen;
         }

       } // for gen jets

       if( deltaRmin<999. ) {

           eJetReco_caloItCone5_[nJetReco_caloItCone5_]  =   eJetReco;
          ptJetReco_caloItCone5_[nJetReco_caloItCone5_]  =  ptJetReco;
         phiJetReco_caloItCone5_[nJetReco_caloItCone5_]  = phiJetReco;
         etaJetReco_caloItCone5_[nJetReco_caloItCone5_]  = etaJetReco;
           eJetGen_caloItCone5_[nJetReco_caloItCone5_]  =   eJetGen_found;
           ptJetGen_caloItCone5_[nJetReco_caloItCone5_]  =  ptJetGen_found;
          phiJetGen_caloItCone5_[nJetReco_caloItCone5_]  = phiJetGen_found;
          etaJetGen_caloItCone5_[nJetReco_caloItCone5_]  = etaJetGen_found;

         
         Float_t deltaRMCmin = 999.;
         Int_t pdgIdPart_found=0;
         Float_t etaPart_found;
         Float_t phiPart_found;
         Float_t ptPart_found;    


         for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

           if( statusMC[iPartMC]!=3 ) continue;

           Float_t px = pxMC[iPartMC];
           Float_t py = pyMC[iPartMC];
           Float_t pz = pzMC[iPartMC];

           Float_t eta = etaMC[iPartMC];
           Float_t phi = phiMC[iPartMC];

           Float_t deltaEtaMC = eta-etaJetGen_found;
           Float_t deltaPhiMC = phi-phiJetGen_found;

           Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );

           if( deltaRMC < deltaRMCmin ) {
             deltaRMCmin = deltaRMC;
             pdgIdPart_found = pdgIdMC[iPartMC];
             etaPart_found = eta;
             phiPart_found = phi;
             Float_t p = sqrt( px*px+ py*py+ pz*pz);
             Float_t theta = 2.*atan(exp(-eta));
             ptPart_found = p*sin(theta);
           }

         } //for MC particles

         pdgIdPart_caloItCone5_[nJetReco_caloItCone5_]  =  pdgIdPart_found;
         ptPart_caloItCone5_[nJetReco_caloItCone5_]  =  ptPart_found;
         phiPart_caloItCone5_[nJetReco_caloItCone5_]  = phiPart_found;
         etaPart_caloItCone5_[nJetReco_caloItCone5_]  = etaPart_found; 
           
         nJetReco_caloItCone5_++;

       } //if deltaRmin<999

     } //for recoJets

if( DEBUG_VERBOSE_ ) std::cout << "Finito calo It cone 5" << std::endl;
     // ------------  CALO SIS CONE 0.5 JETS

     for(unsigned int iRecoJet=0; iRecoJet<nJetCalo_sisCone5; ++iRecoJet) {

       if( nJetReco_caloSISCone5_>=40 ) continue;

       Float_t  pxJetReco  =  pxJetCalo_sisCone5[iRecoJet];
       Float_t  pyJetReco  =  pyJetCalo_sisCone5[iRecoJet];
       Float_t  pzJetReco  =  pzJetCalo_sisCone5[iRecoJet];
       Float_t  phiJetReco = phiJetCalo_sisCone5[iRecoJet];
       Float_t  etaJetReco = etaJetCalo_sisCone5[iRecoJet];
       Float_t  eJetReco = eJetCalo_sisCone5[iRecoJet];

       Float_t pJetReco = sqrt( pxJetReco*pxJetReco + pyJetReco*pyJetReco + pzJetReco*pzJetReco );
       Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
       Float_t ptJetReco = pJetReco*sin(thetaJetReco);

       if ( ptJetReco < 4. ) continue;

       Float_t deltaRmin = 999.;
       Float_t etaJetGen_found;
       Float_t phiJetGen_found;
       Float_t ptJetGen_found;
       Float_t eJetGen_found;

       for(unsigned int iGenJet=0; iGenJet<nJetGen_sisCone5; ++iGenJet) {

         Float_t  pxJetGen  =  pxJetGen_sisCone5[iGenJet];
         Float_t  pyJetGen  =  pyJetGen_sisCone5[iGenJet];
         Float_t  pzJetGen  =  pzJetGen_sisCone5[iGenJet];
         Float_t  phiJetGen = phiJetGen_sisCone5[iGenJet];
         Float_t  etaJetGen = etaJetGen_sisCone5[iGenJet];
         Float_t  eJetGen = eJetGen_sisCone5[iGenJet];

         Float_t deltaEta = etaJetReco - etaJetGen;
         Float_t deltaPhi = phiJetReco - phiJetGen;

         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

         if( deltaR < deltaRmin ) {
           deltaRmin = deltaR;
           etaJetGen_found = etaJetGen;
           phiJetGen_found = phiJetGen;
           Float_t pJetGen = sqrt( pxJetGen*pxJetGen + pyJetGen*pyJetGen + pzJetGen*pzJetGen );
           Float_t thetaJetGen = 2.*atan(exp(-etaJetGen));
           ptJetGen_found = pJetGen*sin(thetaJetGen);
           eJetGen_found = eJetGen;
         }

       } // for gen jets

       if( deltaRmin<999. ) {

           eJetReco_caloSISCone5_[nJetReco_caloSISCone5_]  =   eJetReco;
          ptJetReco_caloSISCone5_[nJetReco_caloSISCone5_]  =  ptJetReco;
         phiJetReco_caloSISCone5_[nJetReco_caloSISCone5_]  = phiJetReco;
         etaJetReco_caloSISCone5_[nJetReco_caloSISCone5_]  = etaJetReco;
           eJetGen_caloSISCone5_[nJetReco_caloSISCone5_]  =   eJetGen_found;
           ptJetGen_caloSISCone5_[nJetReco_caloSISCone5_]  =  ptJetGen_found;
          phiJetGen_caloSISCone5_[nJetReco_caloSISCone5_]  = phiJetGen_found;
          etaJetGen_caloSISCone5_[nJetReco_caloSISCone5_]  = etaJetGen_found;
           
         Float_t deltaRMCmin = 999.;
         Int_t pdgIdPart_found=0;
         Float_t etaPart_found;
         Float_t phiPart_found;
         Float_t ptPart_found;


         for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

           if( statusMC[iPartMC]!=3 ) continue;

           Float_t px = pxMC[iPartMC];
           Float_t py = pyMC[iPartMC];
           Float_t pz = pzMC[iPartMC];
         
           Float_t eta = etaMC[iPartMC];
           Float_t phi = phiMC[iPartMC];
         
           Float_t deltaEtaMC = eta-etaJetGen_found;
           Float_t deltaPhiMC = phi-phiJetGen_found;
         
           Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );
         
           if( deltaRMC < deltaRMCmin ) {
             deltaRMCmin = deltaRMC;
             pdgIdPart_found = pdgIdMC[iPartMC];
             etaPart_found = eta;
             phiPart_found = phi;
             Float_t p = sqrt( px*px+ py*py+ pz*pz);
             Float_t theta = 2.*atan(exp(-eta));
             ptPart_found = p*sin(theta);
           }

         } //for MC particles

         pdgIdPart_caloSISCone5_[nJetReco_caloSISCone5_]  =  pdgIdPart_found;
         ptPart_caloSISCone5_[nJetReco_caloSISCone5_]  =  ptPart_found;
         phiPart_caloSISCone5_[nJetReco_caloSISCone5_]  = phiPart_found;
         etaPart_caloSISCone5_[nJetReco_caloSISCone5_]  = etaPart_found;

         nJetReco_caloSISCone5_++;

         nJetReco_caloSISCone5_++;

       } //if deltaRmin<999

     } //for recoJets

if( DEBUG_VERBOSE_ ) std::cout << "Finito calo sis cone 5" << std::endl;

     // ------------  CALO kT 0.4 JETS

     for(unsigned int iRecoJet=0; iRecoJet<nJetCalo_kt4; ++iRecoJet) {

       if( nJetReco_calokt4_>=40 ) continue;

       Float_t  pxJetReco  =  pxJetCalo_kt4[iRecoJet];
       Float_t  pyJetReco  =  pyJetCalo_kt4[iRecoJet];
       Float_t  pzJetReco  =  pzJetCalo_kt4[iRecoJet];
       Float_t  phiJetReco = phiJetCalo_kt4[iRecoJet];
       Float_t  etaJetReco = etaJetCalo_kt4[iRecoJet];
       Float_t  eJetReco = eJetCalo_kt4[iRecoJet];

       Float_t pJetReco = sqrt( pxJetReco*pxJetReco + pyJetReco*pyJetReco + pzJetReco*pzJetReco );
       Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
       Float_t ptJetReco = pJetReco*sin(thetaJetReco);

       if ( ptJetReco < 4. ) continue;

       Float_t deltaRmin = 999.;
       Float_t etaJetGen_found;
       Float_t phiJetGen_found;
       Float_t ptJetGen_found;
       Float_t eJetGen_found;

       for(unsigned int iGenJet=0; iGenJet<nJetGen_kt4; ++iGenJet) {

         Float_t  pxJetGen  =  pxJetGen_kt4[iGenJet];
         Float_t  pyJetGen  =  pyJetGen_kt4[iGenJet];
         Float_t  pzJetGen  =  pzJetGen_kt4[iGenJet];
         Float_t  phiJetGen = phiJetGen_kt4[iGenJet];
         Float_t  etaJetGen = etaJetGen_kt4[iGenJet];
         Float_t  eJetGen = eJetGen_kt4[iGenJet];

         Float_t deltaEta = etaJetReco - etaJetGen;
         Float_t deltaPhi = phiJetReco - phiJetGen;

         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

         if( deltaR < deltaRmin ) {
           deltaRmin = deltaR;
           etaJetGen_found = etaJetGen;
           phiJetGen_found = phiJetGen;
           Float_t pJetGen = sqrt( pxJetGen*pxJetGen + pyJetGen*pyJetGen + pzJetGen*pzJetGen );
           Float_t thetaJetGen = 2.*atan(exp(-etaJetGen));
           ptJetGen_found = pJetGen*sin(thetaJetGen);
           eJetGen_found = eJetGen;
         }

       } // for gen jets

       if( deltaRmin<999. ) {

           eJetReco_calokt4_[nJetReco_calokt4_]  =   eJetReco;
          ptJetReco_calokt4_[nJetReco_calokt4_]  =  ptJetReco;
         phiJetReco_calokt4_[nJetReco_calokt4_]  = phiJetReco;
         etaJetReco_calokt4_[nJetReco_calokt4_]  = etaJetReco;
            eJetGen_calokt4_[nJetReco_calokt4_]  =   eJetGen_found;
           ptJetGen_calokt4_[nJetReco_calokt4_]  =  ptJetGen_found;
          phiJetGen_calokt4_[nJetReco_calokt4_]  = phiJetGen_found;
          etaJetGen_calokt4_[nJetReco_calokt4_]  = etaJetGen_found;
           
         Float_t deltaRMCmin = 999.;
         Int_t pdgIdPart_found=0;
         Float_t etaPart_found;
         Float_t phiPart_found;
         Float_t ptPart_found;


         for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

           if( statusMC[iPartMC]!=3 ) continue;

           Float_t px = pxMC[iPartMC];
           Float_t py = pyMC[iPartMC];
           Float_t pz = pzMC[iPartMC];
         
           Float_t eta = etaMC[iPartMC];
           Float_t phi = phiMC[iPartMC];
         
           Float_t deltaEtaMC = eta-etaJetGen_found;
           Float_t deltaPhiMC = phi-phiJetGen_found;
         
           Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );
         
           if( deltaRMC < deltaRMCmin ) {
             deltaRMCmin = deltaRMC;
             pdgIdPart_found = pdgIdMC[iPartMC];
             etaPart_found = eta;
             phiPart_found = phi;
             Float_t p = sqrt( px*px+ py*py+ pz*pz);
             Float_t theta = 2.*atan(exp(-eta));
             ptPart_found = p*sin(theta);
           }

         } //for MC particles

         pdgIdPart_calokt4_[nJetReco_calokt4_]  =  pdgIdPart_found;
         ptPart_calokt4_[nJetReco_calokt4_]  =  ptPart_found;
         phiPart_calokt4_[nJetReco_calokt4_]  = phiPart_found;
         etaPart_calokt4_[nJetReco_calokt4_]  = etaPart_found;

         nJetReco_calokt4_++;

         nJetReco_calokt4_++;

       } //if deltaRmin<999

     } //for recoJets


if(DEBUG_VERBOSE_  ) std::cout << "Finito calo kt 4" << std::endl;


     // ------------  PF ITERATIVE CONE 0.5 JETS

     for(unsigned int iRecoJet=0; iRecoJet<nJetPF_itCone5; ++iRecoJet) {

       if( nJetReco_PFItCone5_>=40 ) continue;

       Float_t  eJetReco  =  eJetPF_itCone5[iRecoJet];
       Float_t  pxJetReco  =  pxJetPF_itCone5[iRecoJet];
       Float_t  pyJetReco  =  pyJetPF_itCone5[iRecoJet];
       Float_t  pzJetReco  =  pzJetPF_itCone5[iRecoJet];
       Float_t  phiJetReco = phiJetPF_itCone5[iRecoJet];
       Float_t  etaJetReco = etaJetPF_itCone5[iRecoJet];

       Float_t pJetReco = sqrt( pxJetReco*pxJetReco + pyJetReco*pyJetReco + pzJetReco*pzJetReco );
       Float_t thetaJetReco = 2.*atan(exp(-etaJetReco));
       Float_t ptJetReco = pJetReco*sin(thetaJetReco);

       if ( ptJetReco < 4. ) continue;

       Float_t deltaRmin = 999.;
       Float_t etaJetGen_found;
       Float_t phiJetGen_found;
       Float_t ptJetGen_found;
       Float_t eJetGen_found;
       Float_t eTracksGen_found;
       Float_t ePhotonsGen_found;
       Float_t eNeutralHadronsGen_found;
       Float_t eMuonsGen_found;

       for(unsigned int iGenJet=0; iGenJet<nJetGen_itCone5; ++iGenJet) {

         Float_t  pxJetGen  =  pxJetGen_itCone5[iGenJet];
         Float_t  pyJetGen  =  pyJetGen_itCone5[iGenJet];
         Float_t  pzJetGen  =  pzJetGen_itCone5[iGenJet];
         Float_t  phiJetGen = phiJetGen_itCone5[iGenJet];
         Float_t  etaJetGen = etaJetGen_itCone5[iGenJet];
         Float_t  eJetGen = eJetGen_itCone5[iGenJet];
         Float_t  eTracksGen = eTracksGen_itCone5[iGenJet];
         Float_t  ePhotonsGen = ePhotonsGen_itCone5[iGenJet];
         Float_t  eNeutralHadronsGen = eNeutralHadronsGen_itCone5[iGenJet];
         Float_t  eMuonsGen = eMuonsGen_itCone5[iGenJet];

         Float_t deltaEta = etaJetReco - etaJetGen;
         Float_t deltaPhi = phiJetReco - phiJetGen;

         Float_t deltaR = sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );

         if( deltaR < deltaRmin ) {
           deltaRmin = deltaR;
           etaJetGen_found = etaJetGen;
           phiJetGen_found = phiJetGen;
           Float_t pJetGen = sqrt( pxJetGen*pxJetGen + pyJetGen*pyJetGen + pzJetGen*pzJetGen );
           Float_t thetaJetGen = 2.*atan(exp(-etaJetGen));
           ptJetGen_found = pJetGen*sin(thetaJetGen);
           eJetGen_found = eJetGen;
           eTracksGen_found = eTracksGen;
           ePhotonsGen_found = ePhotonsGen;
           eNeutralHadronsGen_found = eNeutralHadronsGen;
           eMuonsGen_found = eMuonsGen;
         }

       } // for gen jets

       if( deltaRmin<999. ) {

           eJetReco_PFItCone5_[nJetReco_PFItCone5_]  =  eJetReco;
          ptJetReco_PFItCone5_[nJetReco_PFItCone5_]  =  ptJetReco;
         phiJetReco_PFItCone5_[nJetReco_PFItCone5_]  = phiJetReco;
         etaJetReco_PFItCone5_[nJetReco_PFItCone5_]  = etaJetReco;
           ptJetGen_PFItCone5_[nJetReco_PFItCone5_]  =  ptJetGen_found;
          phiJetGen_PFItCone5_[nJetReco_PFItCone5_]  = phiJetGen_found;
          etaJetGen_PFItCone5_[nJetReco_PFItCone5_]  = etaJetGen_found;
            eJetGen_PFItCone5_[nJetReco_PFItCone5_]  =   eJetGen_found;

         eTracksGen_PFItCone5_[nJetReco_PFItCone5_] = eTracksGen_found;
         ePhotonsGen_PFItCone5_[nJetReco_PFItCone5_] = ePhotonsGen_found;
         eNeutralHadronsGen_PFItCone5_[nJetReco_PFItCone5_] = eNeutralHadronsGen_found;
         eMuonsGen_PFItCone5_[nJetReco_PFItCone5_] = eMuonsGen_found;

         Float_t deltaRMCmin = 999.;
         Int_t pdgIdPart_found=0;
         Float_t etaPart_found;
         Float_t phiPart_found;
         Float_t ptPart_found;


         for(Int_t iPartMC=0; iPartMC<nMC; ++iPartMC) {

           if( statusMC[iPartMC]!=3 ) continue;

           Float_t px = pxMC[iPartMC];
           Float_t py = pyMC[iPartMC];
           Float_t pz = pzMC[iPartMC];
         
           Float_t eta = etaMC[iPartMC];
           Float_t phi = phiMC[iPartMC];
         
           Float_t deltaEtaMC = eta-etaJetGen_found;
           Float_t deltaPhiMC = phi-phiJetGen_found;
         
           Float_t deltaRMC = sqrt( deltaEtaMC*deltaEtaMC + deltaPhiMC*deltaPhiMC );
         
           if( deltaRMC < deltaRMCmin ) {
             deltaRMCmin = deltaRMC;
             pdgIdPart_found = pdgIdMC[iPartMC];
             etaPart_found = eta;
             phiPart_found = phi;
             Float_t p = sqrt( px*px+ py*py+ pz*pz);
             Float_t theta = 2.*atan(exp(-eta));
             ptPart_found = p*sin(theta);
           }

         } //for MC particles

         pdgIdPart_PFItCone5_[nJetReco_PFItCone5_]  =  pdgIdPart_found;
         ptPart_PFItCone5_[nJetReco_PFItCone5_]  =  ptPart_found;
         phiPart_PFItCone5_[nJetReco_PFItCone5_]  = phiPart_found;
         etaPart_PFItCone5_[nJetReco_PFItCone5_]  = etaPart_found;


         int  firstCand    = firstCandIndexJetPF_itCone5[iRecoJet];
         int  lastCand    = lastCandIndexJetPF_itCone5[iRecoJet];

         //found a bug. ugly but efficient bug fix:
         if( lastCand==249 ) continue;

         int nTracks = 0;
         Float_t trackEnergy = 0.;
         Float_t trackEt = 0.;
         int nPhotons = 0;
         Float_t photonEnergy = 0.;
         Float_t photonEt = 0.;
         int nNeutralHadrons = 0;
         Float_t neutralHadronEnergy = 0.;
         Float_t neutralHadronEt = 0.;
         int nMuons = 0;
         Float_t muonEnergy = 0.;
         Float_t muonEt = 0.;
         int nHFhadrons = 0;
         Float_t HFhadronEnergy = 0.;
         Float_t HFhadronEt = 0.;
    
         Float_t trackPx=0.;
         Float_t trackPy=0.;
         Float_t trackPz=0.;
         Float_t photonPx=0.;
         Float_t photonPy=0.;
         Float_t photonPz=0.;
         Float_t nhPx=0.;
         Float_t nhPy=0.;
         Float_t nhPz=0.;

         for(int iCand = firstCand; iCand<=lastCand; ++iCand) {

           Float_t etaCand = etaCandPF_itCone5[iCand];
           Float_t thetaCand = 2.*atan(exp(-etaCand));

           Float_t eCand = eCandPF_itCone5[iCand];
           Float_t etCand = eCand*sin(thetaCand);

           Float_t pxCandTrk = pxTrackCandPF_itCone5[iCand];
           Float_t pyCandTrk = pyTrackCandPF_itCone5[iCand];
           Float_t pzCandTrk = pzTrackCandPF_itCone5[iCand];

           Float_t pCandTrk = sqrt(pxCandTrk*pxCandTrk + pyCandTrk*pyCandTrk + pzCandTrk*pzCandTrk);
           Float_t ptCandTrk = pCandTrk*sin(thetaCand);

           Float_t pxCand = pxCandPF_itCone5[iCand];
           Float_t pyCand = pyCandPF_itCone5[iCand];
           Float_t pzCand = pzCandPF_itCone5[iCand];

           Float_t pCand = sqrt(pxCand*pxCand + pyCand*pyCand + pzCand*pzCand);
           Float_t ptCand = pCand*sin(thetaCand);

                     
           if( DEBUG_VERBOSE_ ) {
             std::cout << std::endl << "pxTrackCand: " << pxCand << std::endl;
	       std::cout << std::endl << "pyTrackCand: " << pyCand << std::endl;
             std::cout << std::endl << "pzTrackCand: " << pzCand << std::endl;
             std::cout << std::endl << "pCand: " << pCand << std::endl;
             std::cout << std::endl << "ptCand: " << ptCand << std::endl;
           }


           trackPx += pxCandTrk;
           trackPy += pyCandTrk;
           trackPz += pzCandTrk;
           if( pCandTrk>0. ) {
             trackEnergy+=eCand;
             trackEt+=etCand;
             nTracks+=1;
           }

           int particleID = pdgIdCandPF_itCone5[iCand];

           if(particleID==3) { //muon
             nMuons+=1;
             muonEnergy += eCand;
             muonEt += etCand;
           }

           if(particleID==4) { //photon
             nPhotons+=1;
             photonEnergy += eCand;
             photonEt += etCand;
             photonPx += pxCand;
             photonPy += pyCand;
             photonPz += pzCand;
           }

           if(particleID==5) { //neutral hadron
             nNeutralHadrons+=1;
             neutralHadronEnergy += eCand;
             neutralHadronEt += etCand;
             nhPx += pxCand;
             nhPy += pyCand;
             nhPz += pzCand;
           }

           if(particleID==6) { //HF hadron
             nHFhadrons+=1;
             HFhadronEnergy += eCand;
             HFhadronEt += etCand;
           }

         } //for candidates
         
         Float_t p_tracks = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz);
         Float_t theta_tracks = acos( trackPz/p_tracks );
         Float_t pt_tracks = p_tracks*sin(theta_tracks);
 
         Float_t p_photons = sqrt( photonPx*photonPx + photonPy*photonPy + photonPz*photonPz);
         Float_t theta_photons = acos( photonPz/p_photons );
         Float_t pt_photons = p_photons*sin(theta_photons);
 
         Float_t p_nh = sqrt( nhPx*nhPx + nhPy*nhPy + nhPz*nhPz);
         Float_t theta_nh = acos( nhPz/p_nh );
         Float_t pt_nh = p_nh*sin(theta_nh);
 
         Float_t pxJetReco_noNH = pxJetReco - nhPx;
         Float_t pyJetReco_noNH = pyJetReco - nhPy;
         Float_t pzJetReco_noNH = pzJetReco - nhPz;
         Float_t pJetReco_noNH = sqrt( pxJetReco_noNH*pxJetReco_noNH + pyJetReco_noNH*pyJetReco_noNH + pzJetReco_noNH*pzJetReco_noNH); 
         Float_t thetaJetReco_noNH = acos( pzJetReco_noNH / pJetReco_noNH );

         etaJetReco_noNH_PFItCone5_[nJetReco_PFItCone5_] = -log( tan( thetaJetReco_noNH/2. ) );
         phiJetReco_noNH_PFItCone5_[nJetReco_PFItCone5_] = atan( pyJetReco_noNH / pxJetReco_noNH );
         ptJetReco_noNH_PFItCone5_[nJetReco_PFItCone5_] = pJetReco_noNH*sin(thetaJetReco_noNH);

         nTracksReco_PFItCone5_[nJetReco_PFItCone5_] = nTracks;
         eTracksReco_PFItCone5_[nJetReco_PFItCone5_] = trackEnergy;
         etTracksReco_PFItCone5_[nJetReco_PFItCone5_] = trackEt;
         ptTracksReco_PFItCone5_[nJetReco_PFItCone5_] = pt_tracks;

         nPhotonsReco_PFItCone5_[nJetReco_PFItCone5_] = nPhotons;
         ePhotonsReco_PFItCone5_[nJetReco_PFItCone5_] = photonEnergy;
         etPhotonsReco_PFItCone5_[nJetReco_PFItCone5_] = photonEt;
         ptPhotonsReco_PFItCone5_[nJetReco_PFItCone5_] = pt_photons;

         nNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_] = nNeutralHadrons;
         eNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_] = neutralHadronEnergy;
         etNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_] = neutralHadronEt;
         ptNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_] = pt_nh;

         nMuonsReco_PFItCone5_[nJetReco_PFItCone5_] = nMuons;
         eMuonsReco_PFItCone5_[nJetReco_PFItCone5_] = muonEnergy;
         etMuonsReco_PFItCone5_[nJetReco_PFItCone5_] = muonEt;

         nHFhadronsReco_PFItCone5_[nJetReco_PFItCone5_] = nHFhadrons;
         eHFhadronsReco_PFItCone5_[nJetReco_PFItCone5_] = HFhadronEnergy;
         etHFhadronsReco_PFItCone5_[nJetReco_PFItCone5_] = HFhadronEt;
   
         nJetReco_PFItCone5_++;

         //std::cout << "-------------------> jentry: " << jentry << "\tjet n.: " << iRecoJet << std::endl;
         /*if( (pt_tracks/ptJetReco)>3. ) {
           std::cout << std::endl <<  "ptTracks/ptJet: " << pt_tracks/ptJetReco << "\tnCand: " << (lastCand-firstCand) << std::endl;
           std::cout << "jentry: " << jentry << "\tjet n.: " << iRecoJet << "/" << nJetPF_itCone5 << std::endl;
         }*/

       } //if deltaRmin<999

     } //for recoJets

if( DEBUG_VERBOSE_ ) std::cout << "finito pf it cone 5" << std::endl;

     jetTree_->Fill(); 

if( DEBUG_VERBOSE_ ) std::cout << "fatto fill" << std::endl;

   } //for entries

} //loop



void TreeAnalyzer::BookStuff() {

  outfile_->cd();

  jetTree_ = new TTree("jetTree", "Reduced Tree for Jet Studies");


  //each reco jet is matched to closest gen jet
  //two vectors are saved so that genJet[i] is the genJet matched to recoJet[i]
  //(this is repeated for every algorithm)

  jetTree_->Branch("nJetReco_caloItCone5",&nJetReco_caloItCone5_,"nJetReco_caloItCone5_/I");
  jetTree_->Branch( "eJetReco_caloItCone5",  eJetReco_caloItCone5_,  "eJetReco_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch( "ptJetReco_caloItCone5",  ptJetReco_caloItCone5_,  "ptJetReco_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch("etaJetReco_caloItCone5", etaJetReco_caloItCone5_, "etaJetReco_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch("phiJetReco_caloItCone5", phiJetReco_caloItCone5_, "phiJetReco_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch( "eJetGen_caloItCone5",  eJetGen_caloItCone5_,  "eJetGen_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch(  "ptJetGen_caloItCone5",   ptJetGen_caloItCone5_,   "ptJetGen_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch( "etaJetGen_caloItCone5",  etaJetGen_caloItCone5_,  "etaJetGen_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch( "phiJetGen_caloItCone5",  phiJetGen_caloItCone5_,  "phiJetGen_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch("pdgIdPart_caloItCone5", pdgIdPart_caloItCone5_, "pdgIdPart_caloItCone5_[nJetReco_caloItCone5_]/I");
  jetTree_->Branch(   "ptPart_caloItCone5",    ptPart_caloItCone5_,    "ptPart_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch(  "etaPart_caloItCone5",   etaPart_caloItCone5_,   "etaPart_caloItCone5_[nJetReco_caloItCone5_]/F");
  jetTree_->Branch(  "phiPart_caloItCone5",   phiPart_caloItCone5_,   "phiPart_caloItCone5_[nJetReco_caloItCone5_]/F");

  jetTree_->Branch("nJetReco_caloSISCone5",&nJetReco_caloSISCone5_,"nJetReco_caloSISCone5_/I");
  jetTree_->Branch( "eJetReco_caloSISCone5",  eJetReco_caloSISCone5_,  "eJetReco_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch( "ptJetReco_caloSISCone5",  ptJetReco_caloSISCone5_,  "ptJetReco_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch("etaJetReco_caloSISCone5", etaJetReco_caloSISCone5_, "etaJetReco_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch("phiJetReco_caloSISCone5", phiJetReco_caloSISCone5_, "phiJetReco_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch( "eJetGen_caloSISCone5",  eJetGen_caloSISCone5_,  "eJetGen_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch(  "ptJetGen_caloSISCone5",   ptJetGen_caloSISCone5_,   "ptJetGen_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch( "etaJetGen_caloSISCone5",  etaJetGen_caloSISCone5_,  "etaJetGen_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch( "phiJetGen_caloSISCone5",  phiJetGen_caloSISCone5_,  "phiJetGen_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch("pdgIdPart_caloSISCone5", pdgIdPart_caloSISCone5_, "pdgIdPart_caloSISCone5_[nJetReco_caloSISCone5_]/I");
  jetTree_->Branch(   "ptPart_caloSISCone5",    ptPart_caloSISCone5_,    "ptPart_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch(  "etaPart_caloSISCone5",   etaPart_caloSISCone5_,   "etaPart_caloSISCone5_[nJetReco_caloSISCone5_]/F");
  jetTree_->Branch(  "phiPart_caloSISCone5",   phiPart_caloSISCone5_,   "phiPart_caloSISCone5_[nJetReco_caloSISCone5_]/F");

  jetTree_->Branch("nJetReco_calokt4",&nJetReco_calokt4_,"nJetReco_calokt4_/I");
  jetTree_->Branch( "eJetReco_calokt4",  eJetReco_calokt4_,  "eJetReco_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch( "ptJetReco_calokt4",  ptJetReco_calokt4_,  "ptJetReco_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch("etaJetReco_calokt4", etaJetReco_calokt4_, "etaJetReco_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch("phiJetReco_calokt4", phiJetReco_calokt4_, "phiJetReco_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch( "eJetGen_calokt4",  eJetGen_calokt4_,  "eJetGen_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch(  "ptJetGen_calokt4",   ptJetGen_calokt4_,   "ptJetGen_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch( "etaJetGen_calokt4",  etaJetGen_calokt4_,  "etaJetGen_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch( "phiJetGen_calokt4",  phiJetGen_calokt4_,  "phiJetGen_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch("pdgIdPart_calokt4", pdgIdPart_calokt4_, "pdgIdPart_calokt4_[nJetReco_calokt4_]/I");
  jetTree_->Branch(   "ptPart_calokt4",    ptPart_calokt4_,    "ptPart_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch(  "etaPart_calokt4",   etaPart_calokt4_,   "etaPart_calokt4_[nJetReco_calokt4_]/F");
  jetTree_->Branch(  "phiPart_calokt4",   phiPart_calokt4_,   "phiPart_calokt4_[nJetReco_calokt4_]/F");

  jetTree_->Branch("nJetReco_PFItCone5",&nJetReco_PFItCone5_,"nJetReco_PFItCone5_/I");
  jetTree_->Branch("eJetReco_PFItCone5",  eJetReco_PFItCone5_,  "eJetReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch( "ptJetReco_PFItCone5",  ptJetReco_PFItCone5_,  "ptJetReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("etaJetReco_PFItCone5", etaJetReco_PFItCone5_, "etaJetReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("phiJetReco_PFItCone5", phiJetReco_PFItCone5_, "phiJetReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch(  "eJetGen_PFItCone5",   eJetGen_PFItCone5_,   "eJetGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch(  "ptJetGen_PFItCone5",   ptJetGen_PFItCone5_,   "ptJetGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch( "etaJetGen_PFItCone5",  etaJetGen_PFItCone5_,  "etaJetGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch( "phiJetGen_PFItCone5",  phiJetGen_PFItCone5_,  "phiJetGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("pdgIdPart_PFItCone5", pdgIdPart_PFItCone5_, "pdgIdPart_PFItCone5_[nJetReco_PFItCone5_]/I");
  jetTree_->Branch(   "ptPart_PFItCone5",    ptPart_PFItCone5_,    "ptPart_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch(  "etaPart_PFItCone5",   etaPart_PFItCone5_,   "etaPart_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch(  "phiPart_PFItCone5",   phiPart_PFItCone5_,   "phiPart_PFItCone5_[nJetReco_PFItCone5_]/F");

  jetTree_->Branch("etaJetReco_noNH_PFItCone5", etaJetReco_noNH_PFItCone5_, "etaJetReco_noNH_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("phiJetReco_noNH_PFItCone5", phiJetReco_noNH_PFItCone5_, "phiJetReco_noNH_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("ptJetReco_noNH_PFItCone5", ptJetReco_noNH_PFItCone5_, "ptJetReco_noNH_PFItCone5_[nJetReco_PFItCone5_]/F");

  jetTree_->Branch("nTracksReco_PFItCone5", nTracksReco_PFItCone5_, "nTracksReco_PFItCone5_[nJetReco_PFItCone5_]/I");
  jetTree_->Branch("eTracksReco_PFItCone5", eTracksReco_PFItCone5_, "eTracksReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("etTracksReco_PFItCone5", etTracksReco_PFItCone5_, "etTracksReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("ptTracksReco_PFItCone5", ptTracksReco_PFItCone5_, "ptTracksReco_PFItCone5_[nJetReco_PFItCone5_]/F");

  jetTree_->Branch("nPhotonsReco_PFItCone5", nPhotonsReco_PFItCone5_, "nPhotonsReco_PFItCone5_[nJetReco_PFItCone5_]/I");
  jetTree_->Branch("ePhotonsReco_PFItCone5", ePhotonsReco_PFItCone5_, "ePhotonsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("etPhotonsReco_PFItCone5", etPhotonsReco_PFItCone5_, "etPhotonsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("ptPhotonsReco_PFItCone5", ptPhotonsReco_PFItCone5_, "ptPhotonsReco_PFItCone5_[nJetReco_PFItCone5_]/F");

  jetTree_->Branch("nNeutralHadronsReco_PFItCone5", nNeutralHadronsReco_PFItCone5_, "nNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_]/I");
  jetTree_->Branch("eNeutralHadronsReco_PFItCone5", eNeutralHadronsReco_PFItCone5_, "eNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("etNeutralHadronsReco_PFItCone5", etNeutralHadronsReco_PFItCone5_, "etNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("ptNeutralHadronsReco_PFItCone5", ptNeutralHadronsReco_PFItCone5_, "ptNeutralHadronsReco_PFItCone5_[nJetReco_PFItCone5_]/F");

  jetTree_->Branch("nMuonsReco_PFItCone5", nMuonsReco_PFItCone5_, "nMuonsReco_PFItCone5_[nJetReco_PFItCone5_]/I");
  jetTree_->Branch("eMuonsReco_PFItCone5", eMuonsReco_PFItCone5_, "eMuonsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("etMuonsReco_PFItCone5", etMuonsReco_PFItCone5_, "etMuonsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
 
  jetTree_->Branch("nHFhadronsReco_PFItCone5", nHFhadronsReco_PFItCone5_, "nHFhadronsReco_PFItCone5_[nJetReco_PFItCone5_]/I");
  jetTree_->Branch("eHFhadronsReco_PFItCone5", eHFhadronsReco_PFItCone5_, "eHFhadronsReco_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("etHFhadronsReco_PFItCone5", etHFhadronsReco_PFItCone5_, "etHFhadronsReco_PFItCone5_[nJetReco_PFItCone5_]/F");

  jetTree_->Branch("eTracksGen_PFItCone5", eTracksGen_PFItCone5_, "eTracksGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("ePhotonsGen_PFItCone5", ePhotonsGen_PFItCone5_, "ePhotonsGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("eNeutralHadronsGen_PFItCone5", eNeutralHadronsGen_PFItCone5_, "eNeutralHadronsGen_PFItCone5_[nJetReco_PFItCone5_]/F");
  jetTree_->Branch("eMuonsGen_PFItCone5", eMuonsGen_PFItCone5_, "eMuonsGen_PFItCone5_[nJetReco_PFItCone5_]/F");
 
} //bookstuff

