#include <iostream>
#include <fstream>
#include "BTagSFUtil.h"
#include "TMath.h"




BTagSFUtil::BTagSFUtil( int seed ) {

  rand_ = new TRandom(seed);

}



void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, const std::string& tagger ) {


  // b quarks:
  if( abs( pdgIdPart ) == 5 ) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    // SF for b's
    float b_SF = 0.9;

    float coin = rand_->Uniform(1.);

    if( coin > b_SF ) {

      // scale it down one btag category:
      if( isBTagged_medium ) isBTagged_medium=false;
      else if( isBTagged_loose ) isBTagged_loose=false;
 
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
    if( isBTagged_loose ) {
      std::string fileName = "SF_light_"+tagger+"M.txt";
      btsf = getSF(fileName, jetpt, jeteta);
    } else  {
      std::string fileName = "SF_light_"+tagger+"L.txt";
      btsf = getSF(fileName, jetpt, jeteta);
    }

    float mistagPercent = ( btsf.SF*btsf.eff - btsf.eff ) / ( 1. - btsf.eff );

    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:
    if( coin < mistagPercent ) {

      if( !isBTagged_loose ) isBTagged_loose = true;
      else if( !isBTagged_medium ) isBTagged_medium = true;

    }

  } //if b-light quark

} //modifyBTagsWithSF




BTagScaleFactor BTagSFUtil::getSF( const std::string& fileName, float jetpt, float jeteta ) {

   BTagScaleFactor btsf;

   bool foundSF = false;

   ifstream ifs(fileName.c_str());

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


