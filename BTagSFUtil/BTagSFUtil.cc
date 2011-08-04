#include <iostream>
#include <fstream>
#include "BTagSFUtil.h"
#include "TMath.h"




BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom3(seed);
  setSFFileName("");

}



void BTagSFUtil::set_fileMedium( TFile* file ) {

  fileMedium_ = file;

  if( fileMedium_==0 ) {
    std::cout << "WARNING!!! File '" << file->GetName() << "' does not exist!! Will not set histograms!!" << std::endl;
    return;
  }

  h2_Medium_BTAGBEFFCORR_ = (TH2D*)fileMedium_->Get("BTAGBEFFCORR");
  h2_Medium_BTAGLEFFCORR_ = (TH2D*)fileMedium_->Get("BTAGLEFFCORR");
  h2_Medium_BTAGLEFF_     = (TH2D*)fileMedium_->Get("BTAGLEFF");

  if( h2_Medium_BTAGBEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGBEFFCORR' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Medium_BTAGLEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFFCORR' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Medium_BTAGLEFF_     == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFF' in file '" << fileMedium_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }

}




void BTagSFUtil::set_fileLoose( TFile* file ) {

  fileLoose_ = file;

  if( fileLoose_==0 ) {
    std::cout << "WARNING!!! File '" << file->GetName() << "' does not exist!! Will not set histograms!!" << std::endl;
    return;
  }

  h2_Loose_BTAGBEFFCORR_ = (TH2D*)fileLoose_->Get("BTAGBEFFCORR");
  h2_Loose_BTAGLEFFCORR_ = (TH2D*)fileLoose_->Get("BTAGLEFFCORR");
  h2_Loose_BTAGLEFF_     = (TH2D*)fileLoose_->Get("BTAGLEFF");

  if( h2_Loose_BTAGBEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGBEFFCORR' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Loose_BTAGLEFFCORR_ == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFFCORR' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }
  if( h2_Loose_BTAGLEFF_     == 0 ) {
    std::cout << "WARNING!!! Didn't find histogram 'BTAGLEFF' in file '" << fileLoose_->GetName() << "'!!! Program will crash soon!" << std::endl;
  }

}

  

void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, int pdgIdPart, float Btageff_SF_l, float Btageff_SF_m, float Btagmistag_SF_l,float Btagmistag_SF_m, float Btagmistag_eff_l, float Btagmistag_eff_m) {
  

  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > Btageff_SF_m ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > Btageff_SF_l ) {isBTagged_loose=false; }//
    }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    float  Btagmistag_SF = 1.0;
    float  Btagmistag_eff = 1.0;

    if( isBTagged_loose ) {
      Btagmistag_SF = Btagmistag_SF_m;
      Btagmistag_eff = Btagmistag_eff_m;
    } else  {
      Btagmistag_SF = Btagmistag_SF_l;
      Btagmistag_eff = Btagmistag_eff_l;
    }

    float mistagPercent = ( Btagmistag_SF*Btagmistag_eff - Btagmistag_eff ) / ( 1. - Btagmistag_eff );

    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:

    if( coin < mistagPercent ) {
      if( !isBTagged_loose ) {isBTagged_loose = true;}
      else if( !isBTagged_medium ) {isBTagged_medium = true; }

    }
    
  } //if light quark

} //modifyBTagsWithSF












void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, float sysSF) {
  

  int iBin_pt = (jetpt<240.) ? h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(jetpt) : h2_Medium_BTAGBEFFCORR_->GetXaxis()->FindBin(239.);
  int iBin_eta = h2_Medium_BTAGBEFFCORR_->GetYaxis()->FindBin(jeteta);



  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 ||  abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    // SF for b's
    float b_SF_Medium = h2_Medium_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
    float b_SF_Loose = h2_Loose_BTAGBEFFCORR_->GetBinContent( iBin_pt, iBin_eta );

    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > b_SF_Medium ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > b_SF_Loose ) isBTagged_loose=false; //
    }


  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    // SF for light quarks
    float light_SF_Medium = h2_Medium_BTAGLEFFCORR_->GetBinContent( iBin_pt, iBin_eta );
    float light_SF_Loose = h2_Loose_BTAGLEFFCORR_->GetBinContent( iBin_pt, iBin_eta );

    float light_eff_Medium = h2_Medium_BTAGLEFF_->GetBinContent( iBin_pt, iBin_eta );
    float light_eff_Loose = h2_Loose_BTAGLEFF_->GetBinContent( iBin_pt, iBin_eta );


    float mistagPercent_Loose = ( light_SF_Loose*light_eff_Loose - light_eff_Loose ) / ( 1. - light_eff_Loose );
    float mistagPercent_Medium = ( light_SF_Medium*light_eff_Medium - light_eff_Medium ) / ( 1. - light_eff_Medium );


    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:
    if( !isBTagged_loose ) {
      if( coin < mistagPercent_Loose ) isBTagged_loose = true;
    } else if( !isBTagged_medium ) {
      if( coin < mistagPercent_Medium ) isBTagged_medium = true; 
    }

    
  } //if light quark

} //modifyBTagsWithSF




//This method of getSF uses the functional form of the mistag SF based on 2011 data. We do not need an eta information becasuse the dependence is flat in eta. 
BTagScaleFactor BTagSFUtil::getSF( const std::string& type, float jetpt ) {

  BTagScaleFactor btsf;

  float SFlight = -1;

  if ( type == "medium"){
  
    SFlight = 1.23344 - 0.00033012*jetpt - 0.000000227661*jetpt*jetpt + 0.00000000154075*jetpt*jetpt*jetpt;
    btsf.SF = SFlight;

  }
  else{
    
    SFlight = 1.14062 - 0.000456396*jetpt + 0.0000000838135*jetpt*jetpt;
    btsf.SF = SFlight;

  }

  btsf.SF_err = 0.;
  btsf.eff = 1.;
  btsf.eff_err = 0.;

  return btsf;


}




BTagScaleFactor BTagSFUtil::getSF( const std::string& fileName, float jetpt, float jeteta ) {

   BTagScaleFactor btsf;

   bool foundSF = false;

   ifstream ifs(fileName.c_str());

   if(!ifs.good()){
     std::cout<<"WARNING ! Didn't find file "<<fileName.c_str()<<std::endl;
   }

   while( ifs.good() && !foundSF ) {

     float etaMin, etaMax, ptMin, ptMax, eff, eff_err, SF, SF_err;
     ifs >> etaMin >> etaMax >> ptMin >> ptMax >> eff >> eff_err >> SF >> SF_err;

     if( fabs(jeteta)>=etaMin && fabs(jeteta)<etaMax && jetpt>=ptMin && (jetpt<ptMax||ptMax==999.) ) {
       btsf.SF = SF;
       btsf.SF_err = SF_err;
       btsf.eff = eff;
       btsf.eff_err = eff_err;
       foundSF = true;
     }

   } // while

   ifs.close();

   if( !foundSF ) {
     std::cout << "WARNING! Didn't find SF in file '" << fileName << "'. Setting it to 1." << std::endl;
     btsf.SF = 1.;
     btsf.SF_err = 0.;
     btsf.eff = 1.;
     btsf.eff_err = 0.;
   }

   return btsf;

}
