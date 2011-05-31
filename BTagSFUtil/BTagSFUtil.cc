#include <fstream>
#include "BTagSFUtil.h"
#include "TMath.h"




BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom3(seed);
  setSFFileName("");
}



void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, double Btageff_SF, double Btagmistag_SF, const std::string& tagger) {


  // b quarks:
  if( abs( pdgIdPart ) == 5 ) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    // SF for b's
    float b_SF = Btageff_SF;

    float coin = rand_->Uniform(1.);
    
    float eff_m = 0.7;
    float eff_l = 0.8;
    float b_SF_l = (b_SF*eff_l-eff_m)/(eff_l-eff_m);
    

    if( isBTagged_medium ){ 
      if( coin > b_SF ) isBTagged_medium=false; //turn medium off, loose is still on
    }
    else if( isBTagged_loose ){
      if( coin > b_SF_l ) isBTagged_loose=false; //
    }

  // light quarks:
  } else if( abs( pdgIdPart)>0 ) { //in data it is 0 (save computing time)

    // no need to upgrade if is light and medium tagged
    if( isBTagged_medium ) return;

    BTagScaleFactor btsf;

    // SF for light quarks:
    // remember: use medium efficiencies if it is loose tagged,
    // use loose efficiencies if it not tagged
    // (don't do anything if it is medium tagged)
    // this is because the jet has a probability of being upgraded
    if(sfFileName_==""){
      //std::cout<<"Empty string with SF filename. Setting automatically to "<<std::flush;
      sfFileName_="SF_light_"+tagger;
      //std::cout<<sfFileName_.c_str()<<std::endl;
    }    

    if( isBTagged_loose ) {
      btsf = getSF(sfFileName_+"M.txt", jetpt, jeteta);
    } else  {
      btsf = getSF(sfFileName_+"L.txt", jetpt, jeteta);
    }

    float sysSF = Btagmistag_SF;

    float mistagPercent = ( (sysSF*btsf.SF)*btsf.eff - btsf.eff ) / ( 1. - btsf.eff );

    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:
    if( coin < mistagPercent ) {

      if( !isBTagged_loose ) {isBTagged_loose = true;}
      else if( !isBTagged_medium ) {isBTagged_medium = true; }

    }
    

  } //if b-light quark

} //modifyBTagsWithSF




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
