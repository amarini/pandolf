#include "RedNtpFinalizer_VHgg.h"





RedNtpFinalizer_VHgg::RedNtpFinalizer_VHgg( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer( "VHgg", dataset ) 
{

  bTaggerType_ = bTaggerType;

  setSelectionType(selectionType);

}

RedNtpFinalizer_VHgg::~RedNtpFinalizer_VHgg()
{
   if (!tree_) return;
   delete tree_->GetCurrentFile();
}





void RedNtpFinalizer_VHgg::finalize()
{


   this->Init();


   std::string fullFlags = selectionType_ + "_" + bTaggerType_;
   this->set_flags(fullFlags); //this is for the outfile name
   this->createOutputFile();

   outFile_->cd();



   TH1D*  h1_njets = new TH1D("njets", "", 11, -0.5, 10.5);
   h1_njets->Sumw2();
   TH1D*  h1_nbjets_loose = new TH1D("nbjets_loose", "", 11, -0.5, 10.5);
   h1_nbjets_loose->Sumw2();
   TH1D*  h1_nbjets_medium = new TH1D("nbjets_medium", "", 11, -0.5, 10.5);
   h1_nbjets_medium->Sumw2();

   TH1D*  h1_posMatchedJet = new TH1D("posMatchedJet", "", 11, -0.5, 10.5);
   h1_posMatchedJet->Sumw2();
   TH1D*  h1_ptMatchedJet = new TH1D("ptMatchedJet", "", 1000, 0., 1000.);
   h1_ptMatchedJet->Sumw2();
   TH1D*  h1_etaMatchedJet = new TH1D("etaMatchedJet", "", 100, -5., 5.);
   h1_etaMatchedJet->Sumw2();
   TH1D*  h1_phiMatchedJet = new TH1D("phiMatchedJet", "", 100, -3.1416, 3.1416);
   h1_phiMatchedJet->Sumw2();



   if (tree_ == 0) return;

   Long64_t nentries = tree_->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if( ientry % 10000 == 0 ) std::cout << "-> Entry " << ientry << " / " << nentries << std::endl;

      nb = tree_->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      bool isMC = ( run<5 );

      
      float eventWeight = 1.;


      if( isMC ) {

        if( nGenEvents_<=0 ) {
          std::cout << std::endl << "-> There must be a problem: nGenEvents is not positive!! Exiting" << std::endl;
          exit(99);
        }

        eventWeight = xSection_ / nGenEvents_ ;
 
        // pu/pt reeventWeighting
        if(dopureeventWeight_) eventWeight *= pu_weight; 
        if(doptreeventWeight_) eventWeight *= pt_weight; 

      }



      // analysis cuts:

      if(npu>=30) continue;

      // if(massggnewvtx<90 || massggnewvtx>190) continue;
      if(massggnewvtx<100 || massggnewvtx>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
         || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      //     if(ptphot1<ptphot1cut) continue; //pt first photon
      if(ptphot2<ptphot2cut_) continue; //pt second photon


      if(ptphot1<ptphot1cut_* massggnewvtx/120.) continue; //pt first photon
//       if(ptphot2<ptphot2cut* massggnewvtx/120.) continue; //pt second photon

      if(pthiggsmincut_>0 && ptgg< pthiggsmincut_) continue; //pt higgs min
      if(pthiggsmaxcut_>0 && ptgg>=pthiggsmaxcut_) continue; //pt higgs max


//    if(ptjet1cut>0 && (ptcorrjet1<ptjet1cut_ || TMath::Abs(etajet1)>4.7)) continue; //pt first jet
//    if(ptjet2cut>0 && (ptcorrjet2<ptjet2cut_ || TMath::Abs(etajet2)>4.7)) continue; //pt second jet


       // photon identification

       // r9 categories
       bool isr9phot1(0), isr9phot2(0);
 
       if(TMath::Abs(etascphot1)<1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)<1.4442 && r9phot2>.94) isr9phot2 = 1;
       if(TMath::Abs(etascphot1)>1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)>1.4442 && r9phot2>.94) isr9phot2 = 1;
 
       if(r9cat_ == 1) {
         if(!isr9phot1 || !isr9phot2) continue;
       } else if (r9cat_ == 0){
         if(isr9phot1 && isr9phot2) continue;
       } 
 
       // photon id
       bool idphot1(0), idphot2(0), looseidphot1(0), looseidphot2(0), pxlphot1(1), pxlphot2(1);
 
 //       if(pixelseedcut) { 
 // 	pxlphot1 = !pid_haspixelseedphot1;
 // 	pxlphot2 = !pid_haspixelseedphot2;
 //       }
       
       idphot1 = (idcicphot1 >= cicselection_);
       idphot2 = (idcicphot2 >= cicselection_);
 
       if(!cs_){ // photon id no control sample
 
       if(cicselection_>0) {
         if(!(idphot1)) continue;
         if(!(idphot2)) continue;
       }else{
         if(!(idphot1 && pxlphot1)) continue;
         if(!(idphot2 && pxlphot2)) continue;
       }
       
       }else{ // photon id for control sample
        
         looseidphot1 = (idcicphot1 > 0 );
         looseidphot2 = (idcicphot2 > 0 );
         //	  if( !( (idphot1 && looseidphot2 && !idphot2) || (idphot2 && looseidphot1 && !idphot1) ) ) continue; 
         // Not perfect should be using the same electronVeto wrt CiC selection (now using matchedPromptEle veto)
         if( !( (idphot1 && !idphot2 && !pid_hasMatchedPromptElephot2) || (idphot2 && !idphot1 && !pid_hasMatchedPromptElephot1) ) ) continue; 
 
       }



       // jets

       int njets_selected = 0;
       int njets_selected_btagloose = 0;
       int njets_selected_btagmedium = 0;

       for( unsigned ijet=0; ijet<njets; ++ijet ) {

         if( ptcorrjet[ijet] < ptjetthresh_ ) continue;
         if( fabs(etajet[ijet]) > etajetthresh_ ) continue;

         if( isMC ) {
           if( partMomPdgIDjet[ijet] == 23 || abs( partMomPdgIDjet[ijet] ) == 24 ) {
             h1_posMatchedJet->Fill( ijet, eventWeight );
             h1_ptMatchedJet->Fill( ptcorrjet[ijet], eventWeight );
             h1_etaMatchedJet->Fill( etajet[ijet], eventWeight );
             h1_phiMatchedJet->Fill( phijet[ijet], eventWeight );
             //index_MatchedJet->push_back(ijet);
           }
         }

         if( btagjprobjet[ijet]>0.275 ) njets_selected_btagloose++;
         if( btagjprobjet[ijet]>0.545 ) njets_selected_btagmedium++;

         njets_selected++;

         //AnalysisJet thisJet;
         //thisJet.set
         //selectedJets.push_back( thisJet );

       } //for jets

       h1_njets->Fill( njets_selected, eventWeight );
       h1_nbjets_loose->Fill( njets_selected_btagloose, eventWeight );
       h1_nbjets_medium->Fill( njets_selected_btagmedium, eventWeight );


   } //for entries



   outFile_->cd();

   h1_njets->Write();
   h1_nbjets_loose->Write();
   h1_nbjets_medium->Write();
   h1_posMatchedJet->Write();
   h1_ptMatchedJet->Write();
   h1_etaMatchedJet->Write();
   h1_phiMatchedJet->Write();

} //finalize


















Int_t RedNtpFinalizer_VHgg::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!tree_) return 0;
   return tree_->GetEntry(entry);
}
Long64_t RedNtpFinalizer_VHgg::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!tree_) return -5;
   Long64_t centry = tree_->LoadTree(entry);
   if (centry < 0) return centry;
   if (tree_->GetTreeNumber() != fCurrent) {
      fCurrent = tree_->GetTreeNumber();
   }
   return centry;
}

void RedNtpFinalizer_VHgg::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   fCurrent = -1;
   tree_->SetMakeClass(1);

   tree_->SetBranchAddress("run", &run, &b_run);
   tree_->SetBranchAddress("event", &event, &b_event);
   tree_->SetBranchAddress("lumi", &lumi, &b_lumi);
   tree_->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   tree_->SetBranchAddress("massgg", &massgg, &b_massgg);
   tree_->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   tree_->SetBranchAddress("ptggnewvtx", &ptggnewvtx, &b_ptggnewvtx);
   tree_->SetBranchAddress("phigg", &phigg, &b_phigg);
   tree_->SetBranchAddress("etagg", &etagg, &b_etagg);
   tree_->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
   tree_->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
   tree_->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
   tree_->SetBranchAddress("deltaRToTrackphot1", &deltaRToTrackphot1, &b_deltaRToTrackphot1);
   tree_->SetBranchAddress("deltaRToTrackphot2", &deltaRToTrackphot2, &b_deltaRToTrackphot2);
   tree_->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
   tree_->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
   tree_->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
   tree_->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
   tree_->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
   tree_->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
   tree_->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
   tree_->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
   tree_->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
   tree_->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
   tree_->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
   tree_->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
   tree_->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
   tree_->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
   tree_->SetBranchAddress("energyErrphot1", &energyErrphot1, &b_energyErrphot1);
   tree_->SetBranchAddress("energyErrphot2", &energyErrphot2, &b_energyErrphot2);
   tree_->SetBranchAddress("energySmearingphot1", &energySmearingphot1, &b_energySmearingphot1);
   tree_->SetBranchAddress("energySmearingphot2", &energySmearingphot2, &b_energySmearingphot2);
   tree_->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
   tree_->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
   tree_->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
   tree_->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
   tree_->SetBranchAddress("promptGamma", &promptGamma, &b_promptGamma);
   tree_->SetBranchAddress("LOGamma", &LOGamma, &b_LOGamma);
   tree_->SetBranchAddress("ISRGamma", &ISRGamma, &b_ISRGamma);
   tree_->SetBranchAddress("FSRGamma", &FSRGamma, &b_FSRGamma);
   tree_->SetBranchAddress("idmvaphot1", &idmvaphot1, &b_idmvaphot1);
   tree_->SetBranchAddress("idmvaphot2", &idmvaphot2, &b_idmvaphot2);
   tree_->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
   tree_->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
   tree_->SetBranchAddress("idcicnoelvetophot1", &idcicnoelvetophot1, &b_idcicnoelvetophot1);
   tree_->SetBranchAddress("idcicnoelvetophot2", &idcicnoelvetophot2, &b_idcicnoelvetophot2);
   tree_->SetBranchAddress("idcicpfphot1", &idcicpfphot1, &b_idcicpfphot1);
   tree_->SetBranchAddress("idcicpfphot2", &idcicpfphot2, &b_idcicpfphot2);
   tree_->SetBranchAddress("idcicpfnoelvetophot1", &idcicpfnoelvetophot1, &b_idcicpfnoelvetophot1);
   tree_->SetBranchAddress("idcicpfnoelvetophot2", &idcicpfnoelvetophot2, &b_idcicpfnoelvetophot2);
   tree_->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
   tree_->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
   tree_->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
   tree_->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
   tree_->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
   tree_->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
   tree_->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
   tree_->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
   tree_->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
   tree_->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
   tree_->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
   tree_->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
   tree_->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
   tree_->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
   tree_->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
   tree_->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
   tree_->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
   tree_->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
   tree_->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
   tree_->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
   tree_->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
   tree_->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
   tree_->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
   tree_->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
   tree_->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
   tree_->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
   tree_->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
   tree_->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
   tree_->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
   tree_->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
   tree_->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
   tree_->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
   tree_->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
   tree_->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
   tree_->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
   tree_->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
   tree_->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
   tree_->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
   tree_->SetBranchAddress("njets", &njets, &b_njets);
   tree_->SetBranchAddress("ptjet", ptjet, &b_ptjet);
   tree_->SetBranchAddress("ptcorrjet", ptcorrjet, &b_ptcorrjet);
   tree_->SetBranchAddress("ecorrjet", ecorrjet, &b_ecorrjet);
   tree_->SetBranchAddress("etajet", etajet, &b_etajet);
   tree_->SetBranchAddress("phijet", phijet, &b_phijet);
   tree_->SetBranchAddress("betajet", betajet, &b_betajet);
   tree_->SetBranchAddress("betastarjet", betastarjet, &b_betastarjet);
   tree_->SetBranchAddress("btagvtxjet", btagvtxjet, &b_btagvtxjet);
   tree_->SetBranchAddress("btagtrkjet", btagtrkjet, &b_btagtrkjet);
   tree_->SetBranchAddress("btagjprobjet", btagjprobjet, &b_btagjprobjet);
   tree_->SetBranchAddress("ptDjet", ptDjet, &b_ptDjet);
   tree_->SetBranchAddress("rmsjet", rmsjet, &b_rmsjet);
   tree_->SetBranchAddress("ntrkjet", ntrkjet, &b_ntrkjet);
   tree_->SetBranchAddress("nneutjet", nneutjet, &b_nneutjet);
   tree_->SetBranchAddress("jetIdSimple_mvajet", jetIdSimple_mvajet, &b_jetIdSimple_mvajet);
   tree_->SetBranchAddress("jetIdFull_mvajet", jetIdFull_mvajet, &b_jetIdFull_mvajet);
   tree_->SetBranchAddress("jetId_dR2Meanjet", jetId_dR2Meanjet, &b_jetId_dR2Meanjet);
   tree_->SetBranchAddress("jetId_betaStarClassicjet", jetId_betaStarClassicjet, &b_jetId_betaStarClassicjet);
   tree_->SetBranchAddress("jetId_frac01jet", jetId_frac01jet, &b_jetId_frac01jet);
   tree_->SetBranchAddress("jetId_frac02jet", jetId_frac02jet, &b_jetId_frac02jet);
   tree_->SetBranchAddress("jetId_frac03jet", jetId_frac03jet, &b_jetId_frac03jet);
   tree_->SetBranchAddress("jetId_frac04jet", jetId_frac04jet, &b_jetId_frac04jet);
   tree_->SetBranchAddress("jetId_frac05jet", jetId_frac05jet, &b_jetId_frac05jet);
   tree_->SetBranchAddress("jetId_betajet", jetId_betajet, &b_jetId_betajet);
   tree_->SetBranchAddress("jetId_betaStarjet", jetId_betaStarjet, &b_jetId_betaStarjet);
   tree_->SetBranchAddress("jetIdCutBased_wpjet", jetIdCutBased_wpjet, &b_jetIdCutBased_wpjet);
   tree_->SetBranchAddress("jetIdSimple_wpjet", jetIdSimple_wpjet, &b_jetIdSimple_wpjet);
   tree_->SetBranchAddress("jetIdFull_wpjet", jetIdFull_wpjet, &b_jetIdFull_wpjet);
   tree_->SetBranchAddress("assjet", assjet, &b_assjet);
   tree_->SetBranchAddress("partPdgIDjet", partPdgIDjet, &b_partPdgIDjet);
   tree_->SetBranchAddress("partMomPdgIDjet", partMomPdgIDjet, &b_partMomPdgIDjet);
   tree_->SetBranchAddress("deltaeta", &deltaeta, &b_deltaeta);
   tree_->SetBranchAddress("zeppenjet", &zeppenjet, &b_zeppenjet);
   tree_->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   tree_->SetBranchAddress("deltaphinewvtx", &deltaphinewvtx, &b_deltaphinewvtx);
   tree_->SetBranchAddress("deltaphigg", &deltaphigg, &b_deltaphigg);
   tree_->SetBranchAddress("invmassjet", &invmassjet, &b_invmassjet);
   tree_->SetBranchAddress("invmass2g1j", &invmass2g1j, &b_invmass2g1j);
   tree_->SetBranchAddress("invmass2g2j", &invmass2g2j, &b_invmass2g2j);
   tree_->SetBranchAddress("pt2g2j", &pt2g2j, &b_pt2g2j);
   tree_->SetBranchAddress("eta2j", &eta2j, &b_eta2j);
   tree_->SetBranchAddress("phi2j", &phi2j, &b_phi2j);
   tree_->SetBranchAddress("pt2j", &pt2j, &b_pt2j);
   tree_->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   tree_->SetBranchAddress("vtxId", &vtxId, &b_vtxId);
   tree_->SetBranchAddress("vtxPos_x", &vtxPos_x, &b_vtxPos_x);
   tree_->SetBranchAddress("vtxPos_y", &vtxPos_y, &b_vtxPos_y);
   tree_->SetBranchAddress("vtxPos_z", &vtxPos_z, &b_vtxPos_z);
   tree_->SetBranchAddress("vtxIdMVA", &vtxIdMVA, &b_vtxIdMVA);
   tree_->SetBranchAddress("vtxIdEvtProb", &vtxIdEvtProb, &b_vtxIdEvtProb);
   tree_->SetBranchAddress("sMet", &sMet, &b_sMet);
   tree_->SetBranchAddress("eMet", &eMet, &b_eMet);
   tree_->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   tree_->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
   tree_->SetBranchAddress("eSmearedMet", &eSmearedMet, &b_eSmearedMet);
   tree_->SetBranchAddress("phiSmearedMet", &phiSmearedMet, &b_phiSmearedMet);
   tree_->SetBranchAddress("eShiftedMet", &eShiftedMet, &b_eShiftedMet);
   tree_->SetBranchAddress("phiShiftedMet", &phiShiftedMet, &b_phiShiftedMet);
   tree_->SetBranchAddress("eShiftedScaledMet", &eShiftedScaledMet, &b_eShiftedScaledMet);
   tree_->SetBranchAddress("phiShiftedScaledMet", &phiShiftedScaledMet, &b_phiShiftedScaledMet);
   tree_->SetBranchAddress("eSmearedShiftedMet", &eSmearedShiftedMet, &b_eSmearedShiftedMet);
   tree_->SetBranchAddress("phiSmearedShiftedMet", &phiSmearedShiftedMet, &b_phiSmearedShiftedMet);
   tree_->SetBranchAddress("eShiftedScaledMetPUcorr", &eShiftedScaledMetPUcorr, &b_eShiftedScaledMetPUcorr);
   tree_->SetBranchAddress("phiShiftedScaledMetPUcorr", &phiShiftedScaledMetPUcorr, &b_phiShiftedScaledMetPUcorr);
   tree_->SetBranchAddress("eSmearedShiftedMePUcorrt", &eSmearedShiftedMePUcorrt, &b_eSmearedShiftedMetPUcorr);
   tree_->SetBranchAddress("phiSmearedShiftedMetPUcorr", &phiSmearedShiftedMetPUcorr, &b_phiSmearedShiftedMetPUcorr);
   tree_->SetBranchAddress("sCorrMet", &sCorrMet, &b_sCorrMet);
   tree_->SetBranchAddress("eCorrMet", &eCorrMet, &b_eCorrMet);
   tree_->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
   tree_->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
   tree_->SetBranchAddress("smuCorrMet", &smuCorrMet, &b_smuCorrMet);
   tree_->SetBranchAddress("emuCorrMet", &emuCorrMet, &b_emuCorrMet);
   tree_->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
   tree_->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
   tree_->SetBranchAddress("sNoHFMet", &sNoHFMet, &b_sNoHFMet);
   tree_->SetBranchAddress("eNoHFMet", &eNoHFMet, &b_eNoHFMet);
   tree_->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
   tree_->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
   tree_->SetBranchAddress("stcMet", &stcMet, &b_stcMet);
   tree_->SetBranchAddress("etcMet", &etcMet, &b_etcMet);
   tree_->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   tree_->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   tree_->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
   tree_->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
   tree_->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
   tree_->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
   tree_->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
   tree_->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
   tree_->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
   tree_->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
   tree_->SetBranchAddress("eassocPfMet", &eassocPfMet, &b_eassocPfMet);
   tree_->SetBranchAddress("phiassocPfMet", &phiassocPfMet, &b_phiassocPfMet);
   tree_->SetBranchAddress("signifassocPfMet", &signifassocPfMet, &b_signifassocPfMet);
   tree_->SetBranchAddress("eassocOtherVtxPfMet", &eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
   tree_->SetBranchAddress("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
   tree_->SetBranchAddress("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
   tree_->SetBranchAddress("etrkPfMet", &etrkPfMet, &b_etrkPfMet);
   tree_->SetBranchAddress("phitrkPfMet", &phitrkPfMet, &b_phitrkPfMet);
   tree_->SetBranchAddress("signiftrkPfMet", &signiftrkPfMet, &b_signiftrkPfMet);
   tree_->SetBranchAddress("ecleanPfMet", &ecleanPfMet, &b_ecleanPfMet);
   tree_->SetBranchAddress("phicleanPfMet", &phicleanPfMet, &b_phicleanPfMet);
   tree_->SetBranchAddress("signifcleanPfMet", &signifcleanPfMet, &b_signifcleanPfMet);
   tree_->SetBranchAddress("ecleanedSaclayPfMet", &ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
   tree_->SetBranchAddress("phicleanedSaclayPfMet", &phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
   tree_->SetBranchAddress("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
   tree_->SetBranchAddress("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
   tree_->SetBranchAddress("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
   tree_->SetBranchAddress("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
   tree_->SetBranchAddress("globalPfSums", &globalPfSums, &b_globalPfSums);
   tree_->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
   tree_->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
   tree_->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   tree_->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   tree_->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
   tree_->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
   tree_->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
   tree_->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
   tree_->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
   tree_->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
   tree_->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   tree_->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   tree_->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
   tree_->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
   tree_->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   tree_->SetBranchAddress("npu", &npu, &b_npu);
   tree_->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
   tree_->SetBranchAddress("xsection", &xsection, &b_xsection);
   tree_->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
   tree_->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
   tree_->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   tree_->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
   tree_->SetBranchAddress("gen_custom_processId", &gen_custom_processId, &b_gen_custom_processId);
   tree_->SetBranchAddress("gen_pt_gamma1", &gen_pt_gamma1, &b_gen_pt_gamma1);
   tree_->SetBranchAddress("gen_pt_gamma2", &gen_pt_gamma2, &b_gen_pt_gamma2);
   tree_->SetBranchAddress("gen_eta_gamma1", &gen_eta_gamma1, &b_gen_eta_gamma1);
   tree_->SetBranchAddress("gen_eta_gamma2", &gen_eta_gamma2, &b_gen_eta_gamma2);
   tree_->SetBranchAddress("gen_phi_gamma1", &gen_phi_gamma1, &b_gen_phi_gamma1);
   tree_->SetBranchAddress("gen_phi_gamma2", &gen_phi_gamma2, &b_gen_phi_gamma2);
   tree_->SetBranchAddress("gen_pt_genjet1", &gen_pt_genjet1, &b_gen_pt_genjet1);
   tree_->SetBranchAddress("gen_pt_genjet2", &gen_pt_genjet2, &b_gen_pt_genjet2);
   tree_->SetBranchAddress("gen_eta_genjet1", &gen_eta_genjet1, &b_gen_eta_genjet1);
   tree_->SetBranchAddress("gen_eta_genjet2", &gen_eta_genjet2, &b_gen_eta_genjet2);
   tree_->SetBranchAddress("gen_phi_genjet1", &gen_phi_genjet1, &b_gen_phi_genjet1);
   tree_->SetBranchAddress("gen_phi_genjet2", &gen_phi_genjet2, &b_gen_phi_genjet2);
   tree_->SetBranchAddress("gen_mass_diphoton", &gen_mass_diphoton, &b_gen_mass_diphoton);
   tree_->SetBranchAddress("gen_pt_diphoton", &gen_pt_diphoton, &b_gen_pt_diphoton);
   tree_->SetBranchAddress("gen_eta_diphoton", &gen_eta_diphoton, &b_gen_eta_diphoton);
   tree_->SetBranchAddress("gen_phi_diphoton", &gen_phi_diphoton, &b_gen_phi_diphoton);
   tree_->SetBranchAddress("gen_mass_dijet", &gen_mass_dijet, &b_gen_mass_dijet);
   tree_->SetBranchAddress("gen_pt_dijet", &gen_pt_dijet, &b_gen_pt_dijet);
   tree_->SetBranchAddress("gen_eta_dijet", &gen_eta_dijet, &b_gen_eta_dijet);
   tree_->SetBranchAddress("gen_phi_dijet", &gen_phi_dijet, &b_gen_phi_dijet);
   tree_->SetBranchAddress("gen_zeppenfeld", &gen_zeppenfeld, &b_gen_zeppenfeld);
   tree_->SetBranchAddress("gen_pt_lep1", &gen_pt_lep1, &b_gen_pt_lep1);
   tree_->SetBranchAddress("gen_pt_lep2", &gen_pt_lep2, &b_gen_pt_lep2);
   tree_->SetBranchAddress("gen_eta_lep1", &gen_eta_lep1, &b_gen_eta_lep1);
   tree_->SetBranchAddress("gen_eta_lep2", &gen_eta_lep2, &b_gen_eta_lep2);
   tree_->SetBranchAddress("gen_phi_lep1", &gen_phi_lep1, &b_gen_phi_lep1);
   tree_->SetBranchAddress("gen_phi_lep2", &gen_phi_lep2, &b_gen_phi_lep2);
   tree_->SetBranchAddress("gen_pid_lep1", &gen_pid_lep1, &b_gen_pid_lep1);
   tree_->SetBranchAddress("gen_pid_lep2", &gen_pid_lep2, &b_gen_pid_lep2);
   tree_->SetBranchAddress("ptele1", &ptele1, &b_ptele1);
   tree_->SetBranchAddress("ptele2", &ptele2, &b_ptele2);
   tree_->SetBranchAddress("etaele1", &etaele1, &b_etaele1);
   tree_->SetBranchAddress("etaele2", &etaele2, &b_etaele2);
   tree_->SetBranchAddress("phiele1", &phiele1, &b_phiele1);
   tree_->SetBranchAddress("phiele2", &phiele2, &b_phiele2);
   tree_->SetBranchAddress("eneele1", &eneele1, &b_eneele1);
   tree_->SetBranchAddress("eneele2", &eneele2, &b_eneele2);
   tree_->SetBranchAddress("sIeIeele1", &sIeIeele1, &b_sIeIeele1);
   tree_->SetBranchAddress("sIeIeele2", &sIeIeele2, &b_sIeIeele2);
   tree_->SetBranchAddress("dphiele1", &dphiele1, &b_dphiele1);
   tree_->SetBranchAddress("dphiele2", &dphiele2, &b_dphiele2);
   tree_->SetBranchAddress("detaele1", &detaele1, &b_detaele1);
   tree_->SetBranchAddress("detaele2", &detaele2, &b_detaele2);
   tree_->SetBranchAddress("mhitsele1", &mhitsele1, &b_mhitsele1);
   tree_->SetBranchAddress("mhitsele2", &mhitsele2, &b_mhitsele2);
   tree_->SetBranchAddress("d0ele1", &d0ele1, &b_d0ele1);
   tree_->SetBranchAddress("d0ele2", &d0ele2, &b_d0ele2);
   tree_->SetBranchAddress("dzele1", &dzele1, &b_dzele1);
   tree_->SetBranchAddress("dzele2", &dzele2, &b_dzele2);
   tree_->SetBranchAddress("invMassele1g1", &invMassele1g1, &b_invMassele1g1);
   tree_->SetBranchAddress("invMassele1g2", &invMassele1g2, &b_invMassele1g2);
   tree_->SetBranchAddress("invMassele2g1", &invMassele2g1, &b_invMassele2g1);
   tree_->SetBranchAddress("invMassele2g2", &invMassele2g2, &b_invMassele2g2);
   tree_->SetBranchAddress("oEmoPele1", &oEmoPele1, &b_oEmoPele1);
   tree_->SetBranchAddress("oEmoPele2", &oEmoPele2, &b_oEmoPele2);
   tree_->SetBranchAddress("mvanotrigele1", &mvanotrigele1, &b_mvanotrigele1);
   tree_->SetBranchAddress("mvanotrigele2", &mvanotrigele2, &b_mvanotrigele2);
   tree_->SetBranchAddress("mvatrigele1", &mvatrigele1, &b_mvatrigele1);
   tree_->SetBranchAddress("mvatrigele2", &mvatrigele2, &b_mvatrigele2);
   tree_->SetBranchAddress("matchconvele1", &matchconvele1, &b_matchconvele1);
   tree_->SetBranchAddress("matchconvele2", &matchconvele2, &b_matchconvele2);
   tree_->SetBranchAddress("chHadIso03ele1", &chHadIso03ele1, &b_chHadIso03ele1);
   tree_->SetBranchAddress("chHadIso03ele2", &chHadIso03ele2, &b_chHadIso03ele2);
   tree_->SetBranchAddress("nHadIso03ele1", &nHadIso03ele1, &b_nHadIso03ele1);
   tree_->SetBranchAddress("nHadIso03ele2", &nHadIso03ele2, &b_nHadIso03ele2);
   tree_->SetBranchAddress("photIso03ele1", &photIso03ele1, &b_photIso03ele1);
   tree_->SetBranchAddress("photIso03ele2", &photIso03ele2, &b_photIso03ele2);
   tree_->SetBranchAddress("pteleloose1", &pteleloose1, &b_pteleloose1);
   tree_->SetBranchAddress("pteleloose2", &pteleloose2, &b_pteleloose2);
   tree_->SetBranchAddress("etaeleloose1", &etaeleloose1, &b_etaeleloose1);
   tree_->SetBranchAddress("etaeleloose2", &etaeleloose2, &b_etaeleloose2);
   tree_->SetBranchAddress("phieleloose1", &phieleloose1, &b_phieleloose1);
   tree_->SetBranchAddress("phieleloose2", &phieleloose2, &b_phieleloose2);
   tree_->SetBranchAddress("eneeleloose1", &eneeleloose1, &b_eneeleloose1);
   tree_->SetBranchAddress("eneeleloose2", &eneeleloose2, &b_eneeleloose2);
   tree_->SetBranchAddress("sIeIeeleloose1", &sIeIeeleloose1, &b_sIeIeeleloose1);
   tree_->SetBranchAddress("sIeIeeleloose2", &sIeIeeleloose2, &b_sIeIeeleloose2);
   tree_->SetBranchAddress("dphieleloose1", &dphieleloose1, &b_dphieleloose1);
   tree_->SetBranchAddress("dphieleloose2", &dphieleloose2, &b_dphieleloose2);
   tree_->SetBranchAddress("detaeleloose1", &detaeleloose1, &b_detaeleloose1);
   tree_->SetBranchAddress("detaeleloose2", &detaeleloose2, &b_detaeleloose2);
   tree_->SetBranchAddress("mhitseleloose1", &mhitseleloose1, &b_mhitseleloose1);
   tree_->SetBranchAddress("mhitseleloose2", &mhitseleloose2, &b_mhitseleloose2);
   tree_->SetBranchAddress("d0eleloose1", &d0eleloose1, &b_d0eleloose1);
   tree_->SetBranchAddress("d0eleloose2", &d0eleloose2, &b_d0eleloose2);
   tree_->SetBranchAddress("dzeleloose1", &dzeleloose1, &b_dzeleloose1);
   tree_->SetBranchAddress("dzeleloose2", &dzeleloose2, &b_dzeleloose2);
   tree_->SetBranchAddress("invMasseleloose1g1", &invMasseleloose1g1, &b_invMasseleloose1g1);
   tree_->SetBranchAddress("invMasseleloose1g2", &invMasseleloose1g2, &b_invMasseleloose1g2);
   tree_->SetBranchAddress("invMasseleloose2g1", &invMasseleloose2g1, &b_invMasseleloose2g1);
   tree_->SetBranchAddress("invMasseleloose2g2", &invMasseleloose2g2, &b_invMasseleloose2g2);
   tree_->SetBranchAddress("oEmoPeleloose1", &oEmoPeleloose1, &b_oEmoPeleloose1);
   tree_->SetBranchAddress("oEmoPeleloose2", &oEmoPeleloose2, &b_oEmoPeleloose2);
   tree_->SetBranchAddress("mvanotrigeleloose1", &mvanotrigeleloose1, &b_mvanotrigeleloose1);
   tree_->SetBranchAddress("mvanotrigeleloose2", &mvanotrigeleloose2, &b_mvanotrigeleloose2);
   tree_->SetBranchAddress("mvatrigeleloose1", &mvatrigeleloose1, &b_mvatrigeleloose1);
   tree_->SetBranchAddress("mvatrigeleloose2", &mvatrigeleloose2, &b_mvatrigeleloose2);
   tree_->SetBranchAddress("matchconveleloose1", &matchconveleloose1, &b_matchconveleloose1);
   tree_->SetBranchAddress("matchconveleloose2", &matchconveleloose2, &b_matchconveleloose2);
   tree_->SetBranchAddress("chHadIso03eleloose1", &chHadIso03eleloose1, &b_chHadIso03eleloose1);
   tree_->SetBranchAddress("chHadIso03eleloose2", &chHadIso03eleloose2, &b_chHadIso03eleloose2);
   tree_->SetBranchAddress("nHadIso03eleloose1", &nHadIso03eleloose1, &b_nHadIso03eleloose1);
   tree_->SetBranchAddress("nHadIso03eleloose2", &nHadIso03eleloose2, &b_nHadIso03eleloose2);
   tree_->SetBranchAddress("photIso03eleloose1", &photIso03eleloose1, &b_photIso03eleloose1);
   tree_->SetBranchAddress("photIso03eleloose2", &photIso03eleloose2, &b_photIso03eleloose2);
   tree_->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   tree_->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   tree_->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   tree_->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   tree_->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   tree_->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   tree_->SetBranchAddress("enemu1", &enemu1, &b_enemu1);
   tree_->SetBranchAddress("enemu2", &enemu2, &b_enemu2);
   tree_->SetBranchAddress("pixhitsmu1", &pixhitsmu1, &b_pixhitsmu1);
   tree_->SetBranchAddress("pixhitsmu2", &pixhitsmu2, &b_pixhitsmu2);
   tree_->SetBranchAddress("trkhitsmu1", &trkhitsmu1, &b_trkhitsmu1);
   tree_->SetBranchAddress("trkhitsmu2", &trkhitsmu2, &b_trkhitsmu2);
   tree_->SetBranchAddress("hitsmu1", &hitsmu1, &b_hitsmu1);
   tree_->SetBranchAddress("hitsmu2", &hitsmu2, &b_hitsmu2);
   tree_->SetBranchAddress("chi2mu1", &chi2mu1, &b_chi2mu1);
   tree_->SetBranchAddress("chi2mu2", &chi2mu2, &b_chi2mu2);
   tree_->SetBranchAddress("matchmu1", &matchmu1, &b_matchmu1);
   tree_->SetBranchAddress("matchmu2", &matchmu2, &b_matchmu2);
   tree_->SetBranchAddress("d0mu1", &d0mu1, &b_d0mu1);
   tree_->SetBranchAddress("d0mu2", &d0mu2, &b_d0mu2);
   tree_->SetBranchAddress("dzmu1", &dzmu1, &b_dzmu1);
   tree_->SetBranchAddress("dzmu2", &dzmu2, &b_dzmu2);
   tree_->SetBranchAddress("chHadmu1", &chHadmu1, &b_chHadmu1);
   tree_->SetBranchAddress("chHadmu2", &chHadmu2, &b_chHadmu2);
   tree_->SetBranchAddress("nHadmu1", &nHadmu1, &b_nHadmu1);
   tree_->SetBranchAddress("nHadmu2", &nHadmu2, &b_nHadmu2);
   tree_->SetBranchAddress("photmu1", &photmu1, &b_photmu1);
   tree_->SetBranchAddress("photmu2", &photmu2, &b_photmu2);
   tree_->SetBranchAddress("puptmu1", &puptmu1, &b_puptmu1);
   tree_->SetBranchAddress("puptmu2", &puptmu2, &b_puptmu2);
   tree_->SetBranchAddress("ptmuloose1", &ptmuloose1, &b_ptmuloose1);
   tree_->SetBranchAddress("ptmuloose2", &ptmuloose2, &b_ptmuloose2);
   tree_->SetBranchAddress("etamuloose1", &etamuloose1, &b_etamuloose1);
   tree_->SetBranchAddress("etamuloose2", &etamuloose2, &b_etamuloose2);
   tree_->SetBranchAddress("phimuloose1", &phimuloose1, &b_phimuloose1);
   tree_->SetBranchAddress("phimuloose2", &phimuloose2, &b_phimuloose2);
   tree_->SetBranchAddress("enemuloose1", &enemuloose1, &b_enemuloose1);
   tree_->SetBranchAddress("enemuloose2", &enemuloose2, &b_enemuloose2);
   tree_->SetBranchAddress("pixhitsmuloose1", &pixhitsmuloose1, &b_pixhitsmuloose1);
   tree_->SetBranchAddress("pixhitsmuloose2", &pixhitsmuloose2, &b_pixhitsmuloose2);
   tree_->SetBranchAddress("trkhitsmuloose1", &trkhitsmuloose1, &b_trkhitsmuloose1);
   tree_->SetBranchAddress("trkhitsmuloose2", &trkhitsmuloose2, &b_trkhitsmuloose2);
   tree_->SetBranchAddress("hitsmuloose1", &hitsmuloose1, &b_hitsmuloose1);
   tree_->SetBranchAddress("hitsmuloose2", &hitsmuloose2, &b_hitsmuloose2);
   tree_->SetBranchAddress("chi2muloose1", &chi2muloose1, &b_chi2muloose1);
   tree_->SetBranchAddress("chi2muloose2", &chi2muloose2, &b_chi2muloose2);
   tree_->SetBranchAddress("matchmuloose1", &matchmuloose1, &b_matchmuloose1);
   tree_->SetBranchAddress("matchmuloose2", &matchmuloose2, &b_matchmuloose2);
   tree_->SetBranchAddress("d0muloose1", &d0muloose1, &b_d0muloose1);
   tree_->SetBranchAddress("d0muloose2", &d0muloose2, &b_d0muloose2);
   tree_->SetBranchAddress("dzmuloose1", &dzmuloose1, &b_dzmuloose1);
   tree_->SetBranchAddress("dzmuloose2", &dzmuloose2, &b_dzmuloose2);
   tree_->SetBranchAddress("chHadmuloose1", &chHadmuloose1, &b_chHadmuloose1);
   tree_->SetBranchAddress("chHadmuloose2", &chHadmuloose2, &b_chHadmuloose2);
   tree_->SetBranchAddress("nHadmuloose1", &nHadmuloose1, &b_nHadmuloose1);
   tree_->SetBranchAddress("nHadmuloose2", &nHadmuloose2, &b_nHadmuloose2);
   tree_->SetBranchAddress("photmuloose1", &photmuloose1, &b_photmuloose1);
   tree_->SetBranchAddress("photmuloose2", &photmuloose2, &b_photmuloose2);
   tree_->SetBranchAddress("puptmuloose1", &puptmuloose1, &b_puptmuloose1);
   tree_->SetBranchAddress("puptmuloose2", &puptmuloose2, &b_puptmuloose2);
}




void RedNtpFinalizer_VHgg::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;


  if( selectionType=="sel0" ) {

   dopureeventWeight_ = true;
   doptreeventWeight_ = true;
   r9cat_ = 1;
   cicselection_ = 4;
   cs_ = false;
   ptphot1cut_ = 50.;
   ptphot2cut_ = 30.;
   pthiggsmincut_ = 0.;
   pthiggsmaxcut_ = 10000.;

   ptjetthresh_ = 20.;
   etajetthresh_ = 20.;

  } else {

    std::cout << std::endl << std::endl << "Selection '" << selectionType << "' currently not implemented. Exiting." << std::endl;
    exit(12345);

  }


}

