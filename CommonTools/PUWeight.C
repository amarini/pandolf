///////////////////////////////////////////////////////////////////////
//
//    FILE: PUWeight.h
//   CLASS: PUWeight
// AUTHORS: I. Gonzalez Caballero
//    DATE: 09/03/2011
//
///////////////////////////////////////////////////////////////////////
#include "PUWeight.h"

// ROOT Includes
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"

// C++ includes
#include <iostream>
using namespace std;


//Set DEBUGPUWEIGHT to 1 to get some debug information. Set it to 2 for more
//detail debug information.
#define DEBUGPUWEIGHT 2

#ifdef DEBUG
#define DEBUGPUWEIGHT 1
#endif


PUWeight::PUWeight(float luminosity, const char* year, const std::string& idealMCPU):
  fData(0),
  fMC(0),
  fWeight(0) {

  //Load Data histogram
  if (!LoadDataHistogram(luminosity, year))
    return;

  //No MC given. Take ideal MC
  IdealMCHistogram(idealMCPU);


  //Calculate Weight
  //CalculateWeight();
}


PUWeight::PUWeight(const char* mcfolder, const char* mcproccess, 
		   float luminosity, const char* year):
  fData(0),
  fMC(0),
  fWeight(0) {

  //Load Data histogram
  if (!LoadDataHistogram(luminosity, year))
    return;
  
  //Load MC Histogram
  if (!LoadMCHistogram(mcfolder, mcproccess))
    return;



  //Calculate Weight
  //CalculateWeight();
}




TH1F* PUWeight::LoadMCHistogram(const char* mcfolder, const char* mcproccess) {
#ifdef DEBUGPUWEIGHT
  cout << ">> Getting pileup for the MC " << mcproccess 
       << " inside " << mcfolder << "..." << endl;
#endif
  
  TString dsfile;
  dsfile.Form("http://www.hep.uniovi.es/jfernan/PUhistos/%s/%s.root", 
	      mcfolder, mcproccess);
#if (DEBUGPUWEIGHT > 1)
  cout << "   + Opening " << dsfile << endl;
#endif
  
  TFile* fds = TFile::Open(dsfile);
  if (!fds) {
    cerr << "ERROR [PUWeight]: Could not open file " << dsfile << "!"  << endl
	 << "                  Revise dataset name (" << mcproccess 
	 << ") or internet connection" << endl;
    return 0;
  }
  
  //Read dataset histogram...
#if (DEBUGPUWEIGHT > 1)
  cout << "   + Looking for histogram..." << endl;
#endif
  
  fMC = (TH1F*) fds->Get("htemp")->Clone("PU_MC");
  if (!fMC) {
    cerr << "ERROR [PUWeight]: Could not find histogram for dataset " << mcproccess << "!"
	 << endl;
    return 0;
  }
  fMC->SetDirectory(0);

  if (fMC->Integral() != 1) {
    cout << "NOTE [PUWeight]: MC histogram is not normalized to 1! Normalizing..."
	 << endl;
    fMC->Scale(1./fMC->Integral());
  }

  fds->Close();
  return fMC;
  
}

void PUWeight::SetMCHistogram(const TH1F* mcHisto) {

  //Read MC histogram...
  TH1F* fMC_tmp = (TH1F*)mcHisto->Clone("PU_MC_tmp");

  if( fMC!=0 ) delete fMC;
  int nBins = fMC_tmp->GetNbinsX();
  float xMax = (float)nBins-0.5;
  fMC = new TH1F("PU_MC", "", nBins, -0.5, xMax);
  for( unsigned iBin=1; iBin<nBins+1; ++iBin )
    fMC->SetBinContent( iBin, fMC_tmp->GetBinContent(iBin) );


  if (fMC->Integral() != 1) {
    cout << "NOTE [PUWeight]: MC histogram is not normalized to 1! Normalizing..."
	 << endl;
    fMC->Scale(1./fMC->Integral());
  }

  //CalculateWeight();
  
}


void PUWeight::SetDataHistogram(const TH1F* dataHisto) {

#ifdef DEBUGPUWEIGHT
  cout << ">> Switching to new data histogram: " << dataHisto->GetName()
       << endl;
#endif

  //Read data histogram...
  TH1F* fData_tmp = (TH1F*)dataHisto->Clone("PU_Data_tmp");

  if( fData!=0 ) delete fData;
  int nBins = fData_tmp->GetNbinsX();
  float xMax = (float)nBins-0.5;
  fData = new TH1F("PU_Data", "", nBins, -0.5, xMax);
  for( unsigned iBin=1; iBin<nBins+1; ++iBin )
    fData->SetBinContent( iBin, fData_tmp->GetBinContent(iBin) );

  if (fData->Integral() != 1) {
    cout << "NOTE [PUWeight]: Data histogram is not normalized to 1! Normalizing..."  << endl;
    fData->Scale(1./fData->Integral());
  }




  fData->SetDirectory(0);
  //CalculateWeight();

}


TH1F* PUWeight::LoadDataHistogram(float luminosity, const char* year) {

#ifdef DEBUGPUWEIGHT
  cout << ">> Getting pileup for the " << luminosity << " pb-1 of data..." 
       << endl;
#endif
  
  TString dtfile;
  TFile* fdt = 0;
  if (luminosity > 0) {
    dtfile.Form("http://www.hep.uniovi.es/jfernan/PUhistos/Data%s/PUdata_%.1f.root", 
		year, luminosity);

  
#if (DEBUGPUWEIGHT > 1)
    cout << "   + Opening " << dtfile << endl;
#endif

    fdt = TFile::Open(dtfile);
    if (!fdt) {
      cerr << "NOTE [PUWeight]: Could not find file " << dtfile << "!"  << endl;
      cerr << "                 Trying default PU profile for data..." << endl;
    }
  }

  if (!fdt) {
    dtfile="http://www.hep.uniovi.es/jfernan/PUhistos/Data2011A/PUdata.root";

#if (DEBUGPUWEIGHT > 1)
    cout << "   + Opening " << dtfile << endl;
#endif

    fdt = TFile::Open(dtfile);
    if (!fdt) {
      cerr << "ERROR [PUWeight]: Could not find default profile in \"" 
	   << dtfile << "\"!"  << endl
	   << "                  Is your internet connection working?" << endl;
      return 0;
    }
  }
  
  //Read data histogram...
  TH1F* fData_tmp = (TH1F*) fdt->Get("pileup")->Clone("PU_Data_tmp");
  if (!fData_tmp) {
    cerr << "ERROR [PUWeight]: Could not find histogram for data!" << endl;
    return 0;
  }

  if( fData_tmp->GetXaxis()->GetNbins()==36 ) {
    fData = new TH1F("PU_Data", "", 35, -0.5, 34.5);
    for( unsigned iBin=1; iBin<26; ++iBin ) 
      fData->SetBinContent( iBin, fData_tmp->GetBinContent(iBin) );
  } else {
    fData = (TH1F*)fData_tmp->Clone("PU_Data");
  }
  
  fData->SetDirectory(0);
  
  if (fData->Integral() != 1) {
    cout << "NOTE [PUWeight]: Data histogram is not normalized to 1! Normalizing..."
	 << endl;
    fData->Scale(1./fData->Integral());
  }

  fdt->Close();

  return fData;
}


/*
TH1F* PUWeight::CalculateWeight() {
  if (fData && fMC) {
    unsigned int nbins = fData->GetXaxis()->GetNbins();
    float xmin  = fData->GetXaxis()->GetXmin();
    float xmax  = fData->GetXaxis()->GetXmax();
    fWeight = new TH1F("PUWeight", "PU Weight", nbins, xmin, xmax);
    fWeight->SetDirectory(0);
    fWeight->Divide(fData, fMC);
  }
  else {
    cerr << "ERROR [PUWeight]: Something weird happened when trying to calculate the weights."
	 << endl 
	 << "                  I could not find the data and/or mc histograms!"
	 << endl;
  }


  return fWeight;
}
*/


float PUWeight::GetWeight(unsigned int pu) const {


  int bin = pu+1;

  int nBinsData = fData->GetNbinsX();
  int nBinsMC = fMC->GetNbinsX();

  float num = 0.;
  float denom = 0.;

  if( bin<nBinsData+1 && bin<nBinsMC+1 ) {
    num = fData->GetBinContent(bin);
    denom = fMC->GetBinContent(bin);
  } else if( bin<nBinsData+1 ) {
    num = fData->GetBinContent(bin);
    denom = fMC->GetBinContent(nBinsMC);
  } else if( bin<nBinsMC+1 ) {
    num = fData->GetBinContent(nBinsMC);
    denom = fMC->GetBinContent(bin);
  } else {
    num = 0.;
    denom = 0.;
  }


  float returnWeight = (denom>0.) ? num/denom : 0.;

  return returnWeight;

}



TH1F* PUWeight::IdealMCHistogram( const std::string& puType) {
  unsigned int nbins = 36;
  float xmin = -0.5;
  float xmax = 35.5;

//if (fData) {
//  nbins = fData->GetXaxis()->GetNbins();
//  xmin  = fData->GetXaxis()->GetXmin();
//  xmax  = fData->GetXaxis()->GetXmax();
//}


  fMC = new TH1F("PU_MC", "PU^{MC} Weight", nbins, xmin, xmax);
 
/*
  float   idealpu[]  = {0.0698146584, 0.0698146584, 0.0698146584, 
		     0.0698146584, 0.0698146584, 0.0698146584,
		     0.0698146584, 0.0698146584, 0.0698146584,
		     0.0698146584, 0.0698146584, 0.0630151648,
		     0.0526654164, 0.0402754482, 0.0292988928,
		     0.0194384503, 0.0122016783, 0.007207042,
		     0.004003637,  0.0020278322, 0.0010739954,
		     0.0004595759, 0.0002229748, 0.0001028162,
		     4.58337152809607E-05};
*/

  float   idealpu[55];

  if( puType=="Spring11_Flat10" ) {

      idealpu[0]   =   0.0698146584;
      idealpu[1]   =   0.0698146584;
      idealpu[2]   =   0.0698146584;
      idealpu[3]   =   0.0698146584;
      idealpu[4]   =   0.0698146584;
      idealpu[5]   =   0.0698146584;
      idealpu[6]   =   0.0698146584;
      idealpu[7]   =   0.0698146584;
      idealpu[8]   =   0.0698146584;
      idealpu[9]   =   0.0698146584;
      idealpu[10]   =   0.0698146584;
      idealpu[11]   =   0.0630151648;
      idealpu[12]   =   0.0526654164;
      idealpu[13]   =   0.0402754482;
      idealpu[14]   =   0.0292988928;
      idealpu[15]   =   0.0194384503;
      idealpu[16]   =   0.0122016783;
      idealpu[17]   =   0.007207042;
      idealpu[18]   =   0.004003637;
      idealpu[19]   =   0.0020278322;
      idealpu[20]   =   0.0010739954;
      idealpu[21]   =   0.0004595759;
      idealpu[22]   =   0.0002229748;
      idealpu[23]   =   0.0001028162;
      idealpu[24]   =   4.58337152809607E-05;
      idealpu[25]   =   0.;
      idealpu[26]   =   0.;
      idealpu[27]   =   0.;
      idealpu[28]   =   0.;
      idealpu[29]   =   0.;
      idealpu[30]   =   0.;
      idealpu[31]   =   0.;
      idealpu[32]   =   0.;
      idealpu[33]   =   0.;
      idealpu[34]   =   0.;
      idealpu[35]   =   0.;
      idealpu[36]   =   0.;
      idealpu[37]   =   0.;
      idealpu[38]   =   0.;
      idealpu[39]   =   0.;
      idealpu[40]   =   0.;
      idealpu[41]   =   0.;
      idealpu[42]   =   0.;
      idealpu[43]   =   0.;
      idealpu[44]   =   0.;
      idealpu[45]   =   0.;
      idealpu[46]   =   0.;
      idealpu[47]   =   0.;
      idealpu[48]   =   0.;
      idealpu[49]   =   0.;
      idealpu[50]   =   0.;
      idealpu[51]   =   0.;
      idealpu[52]   =   0.;
      idealpu[53]   =   0.;
      idealpu[54]   =   0.;

  } else if( puType=="Summer11_S4" ) {

      idealpu[0] = 1.45346E-01;
      idealpu[1] = 6.42802E-02;
      idealpu[2] = 6.95255E-02;
      idealpu[3] = 6.96747E-02;
      idealpu[4] = 6.92955E-02;
      idealpu[5] = 6.84997E-02;
      idealpu[6] = 6.69528E-02;
      idealpu[7] = 6.45515E-02;
      idealpu[8] = 6.09865E-02;
      idealpu[9] = 5.63323E-02;
      idealpu[10] = 5.07322E-02;
      idealpu[11] = 4.44681E-02;
      idealpu[12] = 3.79205E-02;
      idealpu[13] = 3.15131E-02;
      idealpu[14] = 2.54220E-02;
      idealpu[15] = 2.00184E-02;
      idealpu[16] = 1.53776E-02;
      idealpu[17] = 1.15387E-02;
      idealpu[18] = 8.47608E-03;
      idealpu[19] = 6.08715E-03;
      idealpu[20] = 4.28255E-03;
      idealpu[21] = 2.97185E-03;
      idealpu[22] = 2.01918E-03;
      idealpu[23] = 1.34490E-03;
      idealpu[24] = 8.81587E-04;
      idealpu[25] = 5.69954E-04;
      idealpu[26] = 3.61493E-04;
      idealpu[27] = 2.28692E-04;
      idealpu[28] = 1.40791E-04;
      idealpu[29] = 8.44606E-05;
      idealpu[30] = 5.10204E-05;
      idealpu[31] = 3.07802E-05;
      idealpu[32] = 1.81401E-05;
      idealpu[33] = 1.00201E-05;
      idealpu[34] = 5.80004E-06;
      idealpu[35] = 0.;
      idealpu[36]   =   0.;
      idealpu[37]   =   0.;
      idealpu[38]   =   0.;
      idealpu[39]   =   0.;
      idealpu[40]   =   0.;
      idealpu[41]   =   0.;
      idealpu[42]   =   0.;
      idealpu[43]   =   0.;
      idealpu[44]   =   0.;
      idealpu[45]   =   0.;
      idealpu[46]   =   0.;
      idealpu[47]   =   0.;
      idealpu[48]   =   0.;
      idealpu[49]   =   0.;
      idealpu[50]   =   0.;
      idealpu[51]   =   0.;
      idealpu[52]   =   0.;
      idealpu[53]   =   0.;
      idealpu[54]   =   0.;

      //idealpu[0] =  0.14551;
      //idealpu[1] =  0.0644453;
      //idealpu[2] =  0.0696412;
      //idealpu[3] =  0.0700311;
      //idealpu[4] =  0.0694257;
      //idealpu[5] =  0.0685655;
      //idealpu[6] =  0.0670929;
      //idealpu[7] =  0.0646049;
      //idealpu[8] =  0.0609383;
      //idealpu[9] =  0.0564597;
      //idealpu[10] =  0.0508014;
      //idealpu[11] =  0.0445226;
      //idealpu[12] =  0.0378796;
      //idealpu[13] =  0.0314746;
      //idealpu[14] =  0.0254139;
      //idealpu[15] =  0.0200091;
      //idealpu[16] =  0.0154191;
      //idealpu[17] =  0.0116242;
      //idealpu[18] =  0.00846857;
      //idealpu[19] =  0.00614328;
      //idealpu[20] =  0.00426355;
      //idealpu[21] =  0.00300632;
      //idealpu[22] =  0.00203485;
      //idealpu[23] =  0.00133045;
      //idealpu[24] =  0.000893794;

  } else if( puType=="Summer11_S4_ave" ) {

      idealpu[0] =  0.104109;
      idealpu[1] =  0.0703573;
      idealpu[2] =  0.0698445;
      idealpu[3] =  0.0698254;
      idealpu[4] =  0.0697054;
      idealpu[5] =  0.0697907;
      idealpu[6] =  0.0696751;
      idealpu[7] =  0.0694486;
      idealpu[8] =  0.0680332;
      idealpu[9] =  0.0651044;
      idealpu[10] =  0.0598036;
      idealpu[11] =  0.0527395;
      idealpu[12] =  0.0439513;
      idealpu[13] =  0.0352202;
      idealpu[14] =  0.0266714;
      idealpu[15] =  0.019411;
      idealpu[16] =  0.0133974;
      idealpu[17] =  0.00898536;
      idealpu[18] =  0.0057516;
      idealpu[19] =  0.00351493;
      idealpu[20] =  0.00212087;
      idealpu[21] =  0.00122891;
      idealpu[22] =  0.00070592;
      idealpu[23] =  0.000384744;
      idealpu[24] =  0.000219377;
      idealpu[25]   =   0.;
      idealpu[26]   =   0.;
      idealpu[27]   =   0.;
      idealpu[28]   =   0.;
      idealpu[29]   =   0.;
      idealpu[30]   =   0.;
      idealpu[31]   =   0.;
      idealpu[32]   =   0.;
      idealpu[33]   =   0.;
      idealpu[34]   =   0.;
      idealpu[35]   =   0.;
      idealpu[36]   =   0.;
      idealpu[37]   =   0.;
      idealpu[38]   =   0.;
      idealpu[39]   =   0.;
      idealpu[40]   =   0.;
      idealpu[41]   =   0.;
      idealpu[42]   =   0.;
      idealpu[43]   =   0.;
      idealpu[44]   =   0.;
      idealpu[45]   =   0.;
      idealpu[46]   =   0.;
      idealpu[47]   =   0.;
      idealpu[48]   =   0.;
      idealpu[49]   =   0.;
      idealpu[50]   =   0.;
      idealpu[51]   =   0.;
      idealpu[52]   =   0.;
      idealpu[53]   =   0.;
      idealpu[54]   =   0.;


  } else if( puType=="Fall11" ) {

     idealpu[0]   =   0.014583699;
     idealpu[1]   =   0.025682975;
     idealpu[2]   =   0.038460562;
     idealpu[3]   =   0.049414536;
     idealpu[4]   =   0.056931087;
     idealpu[5]   =   0.061182816;
     idealpu[6]   =   0.062534625;
     idealpu[7]   =   0.061476918;
     idealpu[8]   =   0.058677499;
     idealpu[9]   =   0.055449877;
     idealpu[10]   =   0.051549051;
     idealpu[11]   =   0.047621024;
     idealpu[12]   =   0.043923799;
     idealpu[13]   =   0.040569076;
     idealpu[14]   =   0.037414654;
     idealpu[15]   =   0.034227033;
     idealpu[16]   =   0.031437714;
     idealpu[17]   =   0.028825596;
     idealpu[18]   =   0.026218978;
     idealpu[19]   =   0.023727061;
     idealpu[20]   =   0.021365645;
     idealpu[21]   =   0.01918743 ;
     idealpu[22]   =   0.016972815;
     idealpu[23]   =   0.014920601;
     idealpu[24]   =   0.013038989;
     idealpu[25]   =   0.011293777;
     idealpu[26]   =   0.009612465;
     idealpu[27]   =   0.008193556;
     idealpu[28]   =   0.006888047;
     idealpu[29]   =   0.005715239;
     idealpu[30]   =   0.004711232;
     idealpu[31]   =   0.003869926;
     idealpu[32]   =   0.003154521;
     idealpu[33]   =   0.002547417;
     idealpu[34]   =   0.002024714;
     idealpu[35]   =   0.001574411;
     idealpu[36]   =   0.001245808;
     idealpu[37]   =   0.000955206;
     idealpu[38]   =   0.000735305;
     idealpu[39]   =   0.000557304;
     idealpu[40]   =   0.000412503;
     idealpu[41]   =   0.000305502;
     idealpu[42]   =   0.000231002;
     idealpu[43]   =   0.000165701;
     idealpu[44]   =   0.000121201;
     idealpu[45]   =   9.30006E-05;
     idealpu[46]   =   6.40004E-05;
     idealpu[47]   =   4.22003E-05;
     idealpu[48]   =   2.85002E-05;
     idealpu[49]   =   1.96001E-05;
     idealpu[50]   =   1.59001E-05;
     idealpu[51]   =   1.01001E-05;
     idealpu[52]   =   8.50006E-06;
     idealpu[53]   =   6.60004E-06;
     idealpu[54]   =   2.70002E-06;

  } else if( puType=="Fall11Truth" ) {

    idealpu[0]   =   0.003388501;
    idealpu[1]   =   0.010357558;
    idealpu[2]   =   0.024724258;
    idealpu[3]   =   0.042348605;
    idealpu[4]   =   0.058279812;
    idealpu[5]   =   0.068851751;
    idealpu[6]   =   0.072914824;
    idealpu[7]   =   0.071579609;
    idealpu[8]   =   0.066811668;
    idealpu[9]   =   0.060672356;
    idealpu[10]   =   0.054528356;
    idealpu[11]   =   0.04919354;
    idealpu[12]   =   0.044886042;
    idealpu[13]   =   0.041341896;
    idealpu[14]   =   0.0384679;
    idealpu[15]   =   0.035871463;
    idealpu[16]   =   0.03341952;
    idealpu[17]   =   0.030915649;
    idealpu[18]   =   0.028395374;
    idealpu[19]   =   0.025798107;
    idealpu[20]   =   0.023237445;
    idealpu[21]   =   0.020602754;
    idealpu[22]   =   0.0180688;
    idealpu[23]   =   0.015559693;
    idealpu[24]   =   0.013211063;
    idealpu[25]   =   0.010964293;
    idealpu[26]   =   0.008920993;
    idealpu[27]   =   0.007080504;
    idealpu[28]   =   0.005499239;
    idealpu[29]   =   0.004187022;
    idealpu[30]   =   0.003096474;
    idealpu[31]   =   0.002237361;
    idealpu[32]   =   0.001566428;
    idealpu[33]   =   0.001074149;
    idealpu[34]   =   0.000721755;
    idealpu[35]   =   0.000470838;
    idealpu[36]   =   0.00030268;
    idealpu[37]   =   0.000184665;
    idealpu[38]   =   0.000112883;
    idealpu[39]   =   6.74043E-05;
    idealpu[40]   =   3.82178E-05;
    idealpu[41]   =   2.22847E-05;
    idealpu[42]   =   1.20933E-05;
    idealpu[43]   =   6.96173E-06;
    idealpu[44]   =   3.4689E-06;
    idealpu[45]   =   1.96172E-06;
    idealpu[46]   =   8.49283E-07;
    idealpu[47]   =   5.02393E-07;
    idealpu[48]   =   2.15311E-07;
    idealpu[49]   =   9.56938E-08;
      idealpu[50]   =   0.;
      idealpu[51]   =   0.;
      idealpu[52]   =   0.;
      idealpu[53]   =   0.;
      idealpu[54]   =   0.;

  } else if( puType=="Summer12" ) {

    idealpu[0] =  2.344E-05,
    idealpu[1] =  2.344E-05,
    idealpu[2] =  2.344E-05,
    idealpu[3] =  2.344E-05,
    idealpu[4] =  4.687E-04,
    idealpu[5] =  4.687E-04,
    idealpu[6] =  7.032E-04,
    idealpu[7] =  9.414E-04,
    idealpu[8] =  1.234E-03,
    idealpu[9] =  1.603E-03,
    idealpu[10] =  2.464E-03,
    idealpu[11] =  3.250E-03,
    idealpu[12] =  5.021E-03,
    idealpu[13] =  6.644E-03,
    idealpu[14] =  8.502E-03,
    idealpu[15] =  1.121E-02,
    idealpu[16] =  1.518E-02,
    idealpu[17] =  2.033E-02,
    idealpu[18] =  2.608E-02,
    idealpu[19] =  3.171E-02,
    idealpu[20] =  3.667E-02,
    idealpu[21] =  4.060E-02,
    idealpu[22] =  4.338E-02,
    idealpu[23] =  4.520E-02,
    idealpu[24] =  4.641E-02,
    idealpu[25] =  4.735E-02,
    idealpu[26] =  4.816E-02,
    idealpu[27] =  4.881E-02,
    idealpu[28] =  4.917E-02,
    idealpu[29] =  4.909E-02,
    idealpu[30] =  4.842E-02,
    idealpu[31] =  4.707E-02,
    idealpu[32] =  4.501E-02,
    idealpu[33] =  4.228E-02,
    idealpu[34] =  3.896E-02,
    idealpu[35] =  3.521E-02,
    idealpu[36] =  3.118E-02,
    idealpu[37] =  2.702E-02,
    idealpu[38] =  2.287E-02,
    idealpu[39] =  1.885E-02,
    idealpu[40] =  1.508E-02,
    idealpu[41] =  1.166E-02,
    idealpu[42] =  8.673E-03,
    idealpu[43] =  6.190E-03,
    idealpu[44] =  4.222E-03,
    idealpu[45] =  2.746E-03,
    idealpu[46] =  1.698E-03,
    idealpu[47] =  9.971E-04,
    idealpu[48] =  5.549E-04,
    idealpu[49] =  2.924E-04,
    idealpu[50] =  1.457E-04,
    idealpu[51] =  6.864E-05,
    idealpu[52] =  3.054E-05,
    idealpu[53] =  1.282E-05,
    idealpu[54] =  5.081E-06
    //idealpu[55] =  1.898E-06,
    //idealpu[56] =  6.688E-07,
    //idealpu[57] =  2.221E-07,
    //idealpu[58] =  6.947E-08,
    //idealpu[59] =  2.047E-08

  } else {

    std::cout << std::endl << "-----> WARNING!!!!!!!!!!" << std::endl;
    std::cout << "Unknown puType: '" << puType << "'. Using Summer11 default." << std::endl;
      idealpu[0] = 1.45346E-01;
      idealpu[1] = 6.42802E-02;
      idealpu[2] = 6.95255E-02;
      idealpu[3] = 6.96747E-02;
      idealpu[4] = 6.92955E-02;
      idealpu[5] = 6.84997E-02;
      idealpu[6] = 6.69528E-02;
      idealpu[7] = 6.45515E-02;
      idealpu[8] = 6.09865E-02;
      idealpu[9] = 5.63323E-02;
      idealpu[10] = 5.07322E-02;
      idealpu[11] = 4.44681E-02;
      idealpu[12] = 3.79205E-02;
      idealpu[13] = 3.15131E-02;
      idealpu[14] = 2.54220E-02;
      idealpu[15] = 2.00184E-02;
      idealpu[16] = 1.53776E-02;
      idealpu[17] = 1.15387E-02;
      idealpu[18] = 8.47608E-03;
      idealpu[19] = 6.08715E-03;
      idealpu[20] = 4.28255E-03;
      idealpu[21] = 2.97185E-03;
      idealpu[22] = 2.01918E-03;
      idealpu[23] = 1.34490E-03;
      idealpu[24] = 8.81587E-04;
      idealpu[25] = 5.69954E-04;
      idealpu[26] = 3.61493E-04;
      idealpu[27] = 2.28692E-04;
      idealpu[28] = 1.40791E-04;
      idealpu[29] = 8.44606E-05;
      idealpu[30] = 5.10204E-05;
      idealpu[31] = 3.07802E-05;
      idealpu[32] = 1.81401E-05;
      idealpu[33] = 1.00201E-05;
      idealpu[34] = 5.80004E-06;
      idealpu[35] =   0.;
      idealpu[36]   =   0.;
      idealpu[37]   =   0.;
      idealpu[38]   =   0.;
      idealpu[39]   =   0.;
      idealpu[40]   =   0.;
      idealpu[41]   =   0.;
      idealpu[42]   =   0.;
      idealpu[43]   =   0.;
      idealpu[44]   =   0.;
      idealpu[45]   =   0.;
      idealpu[46]   =   0.;
      idealpu[47]   =   0.;
      idealpu[48]   =   0.;
      idealpu[49]   =   0.;
      idealpu[50]   =   0.;
      idealpu[51]   =   0.;
      idealpu[52]   =   0.;
      idealpu[53]   =   0.;
      idealpu[54]   =   0.;

  }

  for (unsigned int i = 0; i < nbins; i++) {
    if (i < nbins)
      fMC->Fill(i, idealpu[i]);
    else
      fMC->Fill(i, idealpu[nbins-1]);
  }
  return fMC;
}
