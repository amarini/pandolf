#include "AnalysisElectron.h"

#include "UserCode/emanuele/CommonTools/include/ElectronEffectiveArea.h"


#include <iostream>

// isolation criteria:


bool AnalysisElectron::isIsolatedVBTF80() {

  double combinedIsoRel = this->combinedIsoRel();

  double combinedIsoRelThresh;
  if( fabs(this->Eta())<1.4442 )
    combinedIsoRelThresh = 0.07;
  else
    combinedIsoRelThresh = 0.06;


  bool isIsolated = (combinedIsoRel < combinedIsoRelThresh);

  return isIsolated;

}


bool AnalysisElectron::isIsolatedVBTF95() {

  double combinedIsoRel = this->combinedIsoRel();

  double combinedIsoRelThresh;
  if( fabs(this->Eta())<1.4442 )
    combinedIsoRelThresh = 0.15;
  else
    combinedIsoRelThresh = 0.1;

  bool isIsolated = (combinedIsoRel < combinedIsoRelThresh);

  return isIsolated;

}


// Cut based ID taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
bool AnalysisElectron::isIsolated2012_CutsLoose() {

  float abseta = fabs(this->Eta());

  // calculate the PU subtracted isolation
  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2012;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGammaAndNeutralHad_ = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;

  float eff_area_ga  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
  float eff_area_ganh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);

  float iso = pfCandChargedIso04;
  //  iso += max<float>(0.,pfCandNeutralIso04Ele[eleIndex]-eff_area_nh*rhoJetsFastJet + pfCandPhotonIso04Ele[eleIndex]-eff_area_ga*rhoJetsFastJet);
  iso += TMath::Max(0.,pfCandNeutralIso04+pfCandPhotonIso04-eff_area_ganh*rhoJetsFastJet);


  bool isIsolated = (  iso/this->Pt() < 0.15 );

  return isIsolated; 

}



bool AnalysisElectron::isIsolated2012_MVAWP95() {

  float abseta = fabs(this->Eta());

  // calculate the PU subtracted isolation
  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2012;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGammaAndNeutralHad_ = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;

  float eff_area_ga  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
  float eff_area_ganh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);

  float iso = pfCandChargedIso04;
  //  iso += max<float>(0.,pfCandNeutralIso04Ele[eleIndex]-eff_area_nh*rhoJetsFastJet + pfCandPhotonIso04Ele[eleIndex]-eff_area_ga*rhoJetsFastJet);
  iso += TMath::Max(0.,pfCandNeutralIso04+pfCandPhotonIso04-eff_area_ganh*rhoJetsFastJet);


  // WP95 taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012#Triggering_electrons_MVA
  bool isIsolated = (  ( this->Pt()<20.  && abseta>=0.0   && abseta<0.8   && iso/this->Pt() < 0.408 ) || 
                       ( this->Pt()<20.  && abseta>=0.8   && abseta<1.479 && iso/this->Pt() < 0.500 ) || 
                       ( this->Pt()<20.  && abseta>=1.479 && abseta<2.5   && iso/this->Pt() < 0.463 ) ||
                       ( this->Pt()>=20. && abseta>=0.0   && abseta<0.8   && iso/this->Pt() < 0.253 ) || 
                       ( this->Pt()>=20. && abseta>=0.8   && abseta<1.479 && iso/this->Pt() < 0.225 ) || 
                       ( this->Pt()>=20. && abseta>=1.479 && abseta<2.5   && iso/this->Pt() < 0.308 ) 
                    );

  //bool isIsolated = (iso/this->Pt() < 0.15);

  return isIsolated; 

}




// electron ID criteria:


bool AnalysisElectron::electronIDVBTF80() {

  double sigmaIetaIeta_thresh80;
  double deltaPhiAtVtx_thresh80;
  double deltaEtaAtVtx_thresh80;
  double hOverE_thresh80;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh80 = 0.01;
    deltaPhiAtVtx_thresh80 = 0.06;
    deltaEtaAtVtx_thresh80 = 0.004;
    hOverE_thresh80 = 0.04;
  } else {
    sigmaIetaIeta_thresh80 = 0.03;
    deltaPhiAtVtx_thresh80 = 0.03; 
    deltaEtaAtVtx_thresh80 = 0.007;
    hOverE_thresh80 = 0.15;
  }


  bool eleID_VBTF80 = (sigmaIetaIeta < sigmaIetaIeta_thresh80) &&
                      (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh80) &&
                      (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh80) &&
                      (hOverE < hOverE_thresh80);

  return eleID_VBTF80;

}



bool AnalysisElectron::passedAdditionalCuts() {

  bool passedAC = ( fBrem>0.15 || (fabs(etaSC)<1. && eOverP>0.95) );

  return passedAC;

}


bool AnalysisElectron::separatedIDVBTF80() {

  double sigmaIetaIeta_thresh80;
  double deltaPhiAtVtx_thresh80;
  double deltaEtaAtVtx_thresh80;
  double hOverE_thresh80;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh80 = 0.01;
    deltaPhiAtVtx_thresh80 = 0.06;
    deltaEtaAtVtx_thresh80 = 0.004;
    hOverE_thresh80 = 0.04;
  } else {
    sigmaIetaIeta_thresh80 = 0.03;
    deltaPhiAtVtx_thresh80 = 0.03;//
    deltaEtaAtVtx_thresh80 = 0.007;
    hOverE_thresh80 = 0.05;//
  }

  bool eleID_VBTF80 = (sigmaIetaIeta < sigmaIetaIeta_thresh80) &&
                      (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh80) &&
                      (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh80) &&
                      (hOverE < hOverE_thresh80);

  return eleID_VBTF80;

}




bool AnalysisElectron::electronIDVBTF95() {

  double sigmaIetaIeta_thresh95;
  double deltaPhiAtVtx_thresh95;
  double deltaEtaAtVtx_thresh95;
  double hOverE_thresh95;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh95 = 0.01;
    deltaPhiAtVtx_thresh95 = 0.8; 
    deltaEtaAtVtx_thresh95 = 0.007;
    hOverE_thresh95 = 0.15; 
  } else {
    sigmaIetaIeta_thresh95 = 0.03;
    deltaPhiAtVtx_thresh95 = 0.7; 
    deltaEtaAtVtx_thresh95 = 0.01; 
    hOverE_thresh95 = 0.15;
  }


  bool eleID_VBTF95 = (sigmaIetaIeta < sigmaIetaIeta_thresh95) &&
                      (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh95) &&
                      (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh95) &&
                      (hOverE < hOverE_thresh95);

  return eleID_VBTF95;

}




bool AnalysisElectron::passedHLT2012() {

  double sigmaIetaIeta_thresh;
  double deltaPhiAtVtx_thresh;
  double deltaEtaAtVtx_thresh;
  double hOverE_thresh;
  double ecalIso_thresh;
  double hcalIso_thresh;
  double trkIso_thresh;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh = 0.01;
    deltaPhiAtVtx_thresh = 0.15; 
    deltaEtaAtVtx_thresh = 0.007;
    hOverE_thresh = 0.12; 
    ecalIso_thresh = 0.2;
    hcalIso_thresh = 0.2;
    trkIso_thresh = 0.2;
  } else {
    sigmaIetaIeta_thresh = 0.03;
    deltaPhiAtVtx_thresh = 0.1; 
    deltaEtaAtVtx_thresh = 0.01; 
    hOverE_thresh = 0.10;
    ecalIso_thresh = 0.2;
    hcalIso_thresh = 0.2;
    trkIso_thresh = 0.2;
  }


  // isolation
  float ecalIsolAbs = 0.0;
  if ( fabs(this->Eta())<1.4442 ) ecalIsolAbs = TMath::Max(0.0, dr03EcalRecHitSumEt-1.0);
  else ecalIsolAbs = dr03EcalRecHitSumEt;
  float ecalIsol = ecalIsolAbs/this->Pt();
  float hcalIsol = dr03HcalTowerSumEt/this->Pt();
  float trkIsol = dr03TkSumPt/this->Pt();


  bool eleID_HLT = (sigmaIetaIeta < sigmaIetaIeta_thresh) &&
                      (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh) &&
                      (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh) &&
                      ( ecalIsol < ecalIso_thresh ) &&
                      ( hcalIsol < hcalIso_thresh ) &&
                      ( trkIsol < trkIso_thresh ) &&
                      (hOverE < hOverE_thresh);

  return eleID_HLT;


}




/*
  bool passedHLT = true;
  bool passedHLTID, passedHLTIso;
  passedHLTID = passedHLTIso = true;
  
  // acceptance
  if( fabs(this->Eta()) > 2.5 ) { passedHLT = false; passedHLTID = false; passedHLTIso = false; }
  if( this->Pt() < 10.  )       { passedHLT = false; passedHLTID = false; passedHLTIso = false; }

  // taking shower shape
  float thisSigmaIeIe = sigmaIetaIeta;
  float scEta = etaSC;

  // barrel or endcap
  bool isEleEB = false;
  if (fabs(scEta)<1.479) isEleEB = true;

  // sigmaIetaIeta
  if ( isEleEB && thisSigmaIeIe>0.01) { passedHLT = false; passedHLTID = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { passedHLT = false; passedHLTID = false; }

  // isolation
  float ecalIsolAbs = 0.0;
  if ( isEleEB ) ecalIsolAbs = TMath::Max(0.0, dr03EcalRecHitSumEt-1.0);
  else ecalIsolAbs = dr03EcalRecHitSumEt;
  float ecalIsol = ecalIsolAbs/this->Pt();
  float hcalIsol = dr03HcalTowerSumEt/this->Pt();
  float trackerIsol = dr03TkSumPt/this->Pt();

  if (ecalIsol>0.2) { passedHLT = false; passedHLTIso = false; }
  if (hcalIsol>0.2) { passedHLT = false; passedHLTIso = false; }
  if (trackerIsol>0.2) { passedHLT = false; passedHLTIso = false; }

  // H/E
  if ( isEleEB && hOverE>0.12) { passedHLT = false; passedHLTID = false; }
  if (!isEleEB && hOverE>0.10) { passedHLT = false; passedHLTID = false; }

  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtx)>0.007) ) { passedHLT = false; passedHLTID = false; }
  if (!isEleEB && (fabs(deltaEtaAtVtx)>0.009) ) { passedHLT = false; passedHLTID = false; }

  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtx)>0.15) ) { passedHLT = false; passedHLTID = false; }
  if (!isEleEB && (fabs(deltaPhiAtVtx)>0.10) ) { passedHLT = false; passedHLTID = false; }

  // full conversion rejection
  if (expInnerLayersGsfTrack>0 || hasMatchedConversion) passedHLT = false;

  // impact parameter cuts
  if (fabs(dxy)>0.02) passedHLT = false;
  if (fabs(dz)>0.10)  passedHLT = false;

  return passedHLT;

}
*/



// Cut based ID taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification

bool AnalysisElectron::electronID2012_CutsLoose() {


  if( !(this->passedHLT2012()) ) return false;

  double sigmaIetaIeta_thresh;
  double deltaPhiAtVtx_thresh;
  double deltaEtaAtVtx_thresh;
  double hOverE_thresh;
  double dxy_thresh;
  double dz_thresh;
  double EpInverseDiff_thresh;

  if( fabs(this->Eta())<1.4442 ) {
    sigmaIetaIeta_thresh = 0.01;
    deltaPhiAtVtx_thresh = 0.15; 
    deltaEtaAtVtx_thresh = 0.007;
    hOverE_thresh = 0.12; 
    dxy_thresh = 0.02; 
    dz_thresh = 0.2; 
    EpInverseDiff_thresh = 0.05;
  } else {
    sigmaIetaIeta_thresh = 0.03;
    deltaPhiAtVtx_thresh = 0.7; 
    deltaEtaAtVtx_thresh = 0.01; 
    hOverE_thresh = 100.; //no cut
    dxy_thresh = 0.04; 
    dz_thresh = 0.2; 
    EpInverseDiff_thresh = 0.05;
  }


  double p_inv = 1./pAtVertex;
  double ecalEnergy = pAtVertex*eOverP;
  double E_inv = 1./ecalEnergy;
  double EpInverseDiff = fabs( E_inv-p_inv );

  bool eleID = (sigmaIetaIeta < sigmaIetaIeta_thresh) &&
                   (fabs(deltaPhiAtVtx) < deltaPhiAtVtx_thresh) &&
                   (fabs(deltaEtaAtVtx) < deltaEtaAtVtx_thresh) &&
                   (dxy < dxy_thresh) &&
                   (dz < dz_thresh) &&
                   (EpInverseDiff < EpInverseDiff_thresh) &&
                   (hOverE < hOverE_thresh);

  return eleID;


}




bool AnalysisElectron::electronID2012_MVA() {

  if( !this->passedHLT2012() ) return false;

  float pt = this->Pt();
  float abseta = fabs(this->Eta());
  float bdt = mvaidtrigEle;


  // WP95 taken from https://twiki.cern.ch/twiki/bin/view/Main/HVVElectronId2012#Triggering_electrons_MVA
  bool passed =  (
                 ( pt <= 20. && abseta >= 0.000 && abseta < 0.800 && bdt > -0.288 ) ||
                 ( pt <= 20. && abseta >= 0.800 && abseta < 1.479 && bdt > 0.039 ) ||
                 ( pt <= 20. && abseta >= 1.479 && abseta < 2.500 && bdt > -0.017 ) ||
                 ( pt >  20. && abseta >= 0.000 && abseta < 0.800 && bdt > 0.858 ) ||
                 ( pt >  20. && abseta >= 0.800 && abseta < 1.479 && bdt > 0.425 ) ||
                 ( pt >  20. && abseta >= 1.479 && abseta < 2.500 && bdt > 0.759 )
                 );

// this is WP HWW:
//bool passed =  (
//               ( pt <= 20. && abseta >= 0.000 && abseta < 0.800 && bdt > 0.00 ) ||
//               ( pt <= 20. && abseta >= 0.800 && abseta < 1.479 && bdt > 0.10 ) ||
//               ( pt <= 20. && abseta >= 1.479 && abseta < 2.500 && bdt > 0.62 ) ||
//               ( pt >  20. && abseta >= 0.000 && abseta < 0.800 && bdt > 0.94 ) ||
//               ( pt >  20. && abseta >= 0.800 && abseta < 1.479 && bdt > 0.85 ) ||
//               ( pt >  20. && abseta >= 1.479 && abseta < 2.500 && bdt > 0.92 )
//               );

  return passed;

}







// conversion rejection criteria: 


bool AnalysisElectron::conversionRejectionVBTF80() {

   bool convRej_VBTF80 = (expInnerLayersGsfTrack<=0) && (fabs(convDist)>0.02 || fabs(convDcot)>0.02);
   
   return convRej_VBTF80;

}

bool AnalysisElectron::conversionRejectionVBTF95() {

   bool convRej_VBTF95 = (expInnerLayersGsfTrack<=1);
   
   return convRej_VBTF95;

}

bool AnalysisElectron::conversionRejection2012_CutsLoose() {

   return (!hasMatchedConversion && expInnerLayersGsfTrack<=1);

}





// electron selection criteria:


bool AnalysisElectron::passedVBTF80() {

  bool isIsolated = this->isIsolatedVBTF80();
  bool electronID = this->electronIDVBTF80();
  bool conversionRejection = this->conversionRejectionVBTF80();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}

bool AnalysisElectron::passedTrigger80() {

  bool isIsolated = this->separatedIsoRel();
  bool electronID = this->separatedIDVBTF80();
  bool conversionRejection = this->conversionRejectionVBTF80();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}

bool AnalysisElectron::passedVBTF95() {

  bool isIsolated = this->isIsolatedVBTF95();
  bool electronID = this->electronIDVBTF95();
  bool conversionRejection = this->conversionRejectionVBTF95();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}


bool AnalysisElectron::isGoodElectron2012_CutsLoose() {

  bool isIsolated = this->isIsolated2012_CutsLoose();
  bool electronID = this->electronID2012_CutsLoose();
  bool conversionRejection = this->conversionRejection2012_CutsLoose();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}

bool AnalysisElectron::isGoodElectron2012_MVA() {

  bool isIsolated = this->isIsolated2012_MVAWP95();
  bool electronID = this->electronID2012_MVA();
  bool conversionRejection = this->conversionRejection2012_CutsLoose();

  bool passed = ( isIsolated && electronID && conversionRejection );

  return passed;

}



double AnalysisElectron::combinedIsoRel() {

  double combinedIsoRel;

  if( fabs(this->Eta())<1.4442 )
    combinedIsoRel = ( dr03TkSumPt + TMath::Max(0., dr03EcalRecHitSumEt - 1.) + dr03HcalTowerSumEt ) / this->Pt();
  else
    combinedIsoRel = ( dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt ) / this->Pt();

  return combinedIsoRel;

}

bool AnalysisElectron::separatedIsoRel() {

  bool separedIsoRel=false;

  if( fabs(this->Eta())<1.4442 ){
   if ( dr03TkSumPt/this->Pt()<0.09 && dr03EcalRecHitSumEt/this->Pt()<0.07 && dr03HcalTowerSumEt/this->Pt()<0.10 ) separedIsoRel=true; }
  else{
   if ( dr03TkSumPt/this->Pt()<0.04 && dr03EcalRecHitSumEt/this->Pt()<0.05 && dr03HcalTowerSumEt/this->Pt()<0.025 ) separedIsoRel=true; }

  return separedIsoRel;

}

