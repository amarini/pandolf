#include <fstream>
#include "BTagSFUtil.h"
#include "TMath.h"



BTagSFUtil::BTagSFUtil( int seed ) {
  rand_ = new TRandom3(seed);
  setSFFileName("");
  cachesUsed_=0;
}



void BTagSFUtil::modifyBTagsWithSF( bool& isBTagged_loose, bool& isBTagged_medium, float jetpt, float jeteta, int pdgIdPart, double Btageff_SF, double Btagmistag_SF, const std::string& tagger, bool verbose) {


  // b quarks and c quarks:
  if( abs( pdgIdPart ) == 5 || abs( pdgIdPart ) == 4) { 

    // no need to downgrade if is a b and not tagged
    if( !isBTagged_loose ) return;

    // SF for b's
    float b_SF = Btageff_SF;

    float coin = rand_->Uniform(1.);
    
    if( isBTagged_medium ){ 
      if( coin > b_SF ) {isBTagged_medium=false; isBTagged_loose=false;} //turn medium and loose off, 
    }
    else if( isBTagged_loose && !isBTagged_medium ){
      if( coin > b_SF ) isBTagged_loose=false; //
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
      //btsf = getSF(sfFileName_+"M.txt", jetpt, jeteta, verbose);
      btsf = getSF("medium", jetpt);
    } else  {
      //btsf = getSF(sfFileName_+"L.txt", jetpt, jeteta, verbose);
      btsf = getSF("loose", jetpt);
    }

    float sysSF = Btagmistag_SF;
    float mistagPercent = sysSF*btsf.SF - 1.0;
    float coin = rand_->Uniform(1.);

    // for light quarks, the jet has to be upgraded:
    if( coin < mistagPercent ) {

      if( !isBTagged_loose ) {isBTagged_loose = true;}
      else if( !isBTagged_medium ) {isBTagged_medium = true; }

    }
    

  } //if light quark

} //modifyBTagsWithSF




//This method of getSF uses the functional form of the mistag SF based on 2011 data.
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
  btsf.etamax = -999.;
  btsf.etamin = -999.;
  btsf.ptmax  = -999.;
  btsf.ptmin  = -999.;

  return btsf;


}




BTagScaleFactor BTagSFUtil::getSF( const std::string& fileName, float jetpt, float jeteta ,bool verbose) {
  // set default values in case we abort or don't find proper values
  BTagScaleFactor btsf_default;
  btsf_default.SF = 1.;
  btsf_default.SF_err = 0.;
  btsf_default.eff = 1.;
  btsf_default.eff_err = 0.;
  btsf_default.etamax = -999.;
  btsf_default.etamin = -999.;
  btsf_default.ptmax  = -999.;
  btsf_default.ptmin  = -999.;
  
  
  std::map<std::string,int>::iterator cacheLine = cacheAssoc_.find(fileName);
  
  if(cacheLine == cacheAssoc_.end() ){// we haven't read this file yet.
    if (cachesUsed_ == 5){// but we are out of further caches to use
      std::cout<<"WARNING !Out of file chaches for Btag scaling "<<fileName.c_str()<<std::endl;
      return btsf_default;
    }

    ifstream ifs(fileName.c_str());
    if(!ifs.good()){ // Warn and return default if file is not found.
      std::cout<<"WARNING ! Didn't find file "<<fileName.c_str()<<std::endl;
       return btsf_default;
    }
    BTagScaleFactor btsf;
    // read the file and fill the cache     
    while( ifs.good() ) {
      float etaMin, etaMax, ptMin, ptMax, eff, eff_err, SF, SF_err;
      ifs >> etaMin >> etaMax >> ptMin >> ptMax >> eff >> eff_err >> SF >> SF_err;
      btsf.etamax =  etaMax;
      btsf.etamin =  etaMin;
      btsf.ptmax =  ptMax;
      btsf.ptmin =  ptMin;
      
      btsf.SF = SF;
      btsf.SF_err = SF_err;
      btsf.eff = eff;
      btsf.eff_err = eff_err;
      
      cache_[cachesUsed_].push_back(btsf);
    }
    ifs.close();
    
    // new cache to the table
    cacheAssoc_.insert(std::pair<std::string,int>(fileName,cachesUsed_));
    cachesUsed_++;
    cacheLine = cacheAssoc_.find(fileName); // refresh
    
  }
  
  int cn = (*cacheLine).second;
    
  for(unsigned int i = 0 ; i < cache_[cn].size() ; ++i){
    BTagScaleFactor btsf = cache_[cn].at(i);
    if( fabs(jeteta)>=btsf.etamin && fabs(jeteta)<btsf.etamax && jetpt>=btsf.ptmin && (jetpt<btsf.ptmax||btsf.ptmax==999.) ) 
      return btsf;    
  }
    
  // This point can only be reached if we have the proper table but didn' find a proper entry
  if(verbose) std::cout << "WARNING! Didn't find SF in file '" << fileName << "'. Setting it to 1." << std::endl; 
  return btsf_default;  
}
