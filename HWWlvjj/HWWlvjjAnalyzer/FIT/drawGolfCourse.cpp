#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"


#include "CommonTools/DrawBase.h"



int main() {

  std::vector<float> masses;
  masses.push_back(300.);
  masses.push_back(325.);
  masses.push_back(350.);
  masses.push_back(375.);
  masses.push_back(400.);
  masses.push_back(420.);
  masses.push_back(440.);
  masses.push_back(460.);
  masses.push_back(480.);
  masses.push_back(500.);
  masses.push_back(520.);
  masses.push_back(540.);
  masses.push_back(560.);
  masses.push_back(580.);
  masses.push_back(600.);


  for( unsigned int iMass=0; iMass<masses.size(); ++iMass ) {

    char expectedFileName[500];
    sprintf( expectedFileName, "datacards/PL_%.0f.log", masses[iMass] );

    ifstream ifs(expectedFileName);
    float limit;

    while( ifs.good() ) {

      //std::string thisLine;
      char thisLine[1000];
      ifs.getline(thisLine, 1000);
      //ifs >> thisLine;
      std::string thisLine_str(thisLine_tstr);

      TString thisLine_tstr(thisLine_str);
      if( thisLine_tstr.BeginsWith("median expected limit") ) {
        std::stringstream sstream(thisLine_str);
        std::string uselessStuff;
        //median expected limit: r < 5.7839 @ 95%CL (97 toyMC)
        sstream >> uselessStuff >> uselessStuff >> uselessStuff >> uselessStuff >> limit;
      } else if( thisLine_tstr.BeginsWith("   68% expected band : ") ) {
        
      }


      std::cout << thisLine << std::endl;

    } //while expected ifs good

    std::cout << "+++ THIS LIMIT: " << limit << std::endl;
  }

}
