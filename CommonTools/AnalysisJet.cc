#include "AnalysisJet.h"

#include <map>
#include <iostream>

bool AnalysisJet::btag_loose() const {

  return trackCountingHighEffBJetTag>1.7;

}


bool AnalysisJet::btag_medium() const {

  return trackCountingHighEffBJetTag>3.3;

}
Bool_t AnalysisJet::passedJetID( const std::string& strength ) const {

 Bool_t returnBool = false;
  
 std::map<string,long> numbers;
   	numbers["minimal"]=1;

switch ( numbers[strength] ) {
	case 1: //minimal
	 	returnBool = ( this->nConstituents() > 1) && 
			( (ePhotons + eHFEM)/Energy() < 0.99 ) &&
			( (eNeutralHadrons)/Energy() < 0.99 );
	
		if( fabs(Eta())<2.4 ) {returnBool = returnBool && ( eChargedHadrons > 0. );}
		break;
	default:
	std::cout << "Jet ID '" << strength << "' not implemented yet. Returning FALSE." << std::endl;
	}

 return returnBool;

}
