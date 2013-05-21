
#include <map>
#include "TFile.h"
#include <string>
#include "TChain.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include <iostream>
#include <stdlib.h>

#define MAXARRAY 1023

using namespace std;

int DrawComparisonFromNtp1(int job=0, int nJobs=1)
{
map<string,TChain*> trees;
map<string,float > xSec;

bool flat=true;

if(flat)
{
trees["flat"] = new TChain("ntp1");
cout << "Added " << trees["flat"]->Add("root://eoscms///store/group/phys_higgs/pandolf/MC/Summer12/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_2/default_MC*.root") << " files to the flat tree"<<endl;
}
else
{
trees["15to30"]=new TChain("ntp1");
trees["30to50"]=new TChain("ntp1");
trees["50to80"]=new TChain("ntp1");
trees["80to120"]=new TChain("ntp1");
trees["120to170"]=new TChain("ntp1");
trees["170to300"]=new TChain("ntp1");
trees["300to470"]=new TChain("ntp1");
trees["470to600"]=new TChain("ntp1");
trees["600to800"]=new TChain("ntp1");
trees["800to1000"]=new TChain("ntp1");

for(map<string,TChain*>::iterator it=trees.begin();it != trees.end();it++)
	{
	cout<<"Added "<<it->second->Add(Form("root://eoscms///store/group/phys_higgs/cmshww/emanuele/ProdSummer12/V00-5_3_X/MCSummer12_DR53X/QCD_Pt-%s_TuneZ2star_8TeV_pythia6/default_MC*.root",it->first.c_str())) <<" files to "<<it->first<<endl;
	}
}


xSec["flat"]=1.0;

xSec["15to30"]=9.8828742E8;
xSec["30to50"]=6.6285328E7;
xSec["50to80"]=8148778.0;
xSec["80to120"]=1033680.0;
xSec["120to170"]=156293.3;
xSec["170to300"]=34138.15;
xSec["300to470"]=1759.549;
xSec["470to600"]=113.8791;
xSec["600to800"]=26.9921;
xSec["800to1000"]=3.550036;

Double_t Bins[]={30,40,50,60,70,80,100,120,140,160,180,200,250,300,350,400,500,600,700,800,1000,1200,1400,1600,1800,2000};
TFile *out=TFile::Open(Form("Plot/output_Comparison_%d_%d.root",job,nJobs),"RECREATE");
system(Form("[ -f Plot/done_%d_%d.txt ] && rm Plot/done_%d_%d.txt",job,nJobs,job,nJobs));
TH1F* PtJet_Q=new TH1F("PtJet_Q","PtJet_Q",25,Bins);
TH1F* PtJet_G=new TH1F("PtJet_G","PtJet_G",25,Bins);
TH1F* PtJet_O=new TH1F("PtJet_O","PtJet_O",25,Bins);

PtJet_Q->Sumw2();PtJet_G->Sumw2();PtJet_O->Sumw2();

TH1F* PtJetFwd_Q=new TH1F("PtJetFwd_Q","PtJetFwd_Q",25,Bins);
TH1F* PtJetFwd_G=new TH1F("PtJetFwd_G","PtJetFwd_G",25,Bins);
TH1F* PtJetFwd_O=new TH1F("PtJetFwd_O","PtJetFwd_O",25,Bins);

PtJetFwd_Q->Sumw2();PtJetFwd_G->Sumw2();PtJetFwd_O->Sumw2();

for(map<string,TChain*>::iterator it=trees.begin();it!=trees.end();it++)
	{
	Int_t nMc;it->second->SetBranchAddress("nMc",&nMc);
	Int_t idMc[MAXARRAY];it->second->SetBranchAddress("idMc",idMc);
	Int_t statusMc[MAXARRAY];it->second->SetBranchAddress("statusMc",statusMc);
	Float_t etaMc[MAXARRAY];it->second->SetBranchAddress("etaMc",etaMc);
	Float_t thetaMc[MAXARRAY];it->second->SetBranchAddress("thetaMc",thetaMc);
	Float_t phiMc[MAXARRAY];it->second->SetBranchAddress("phiMc",phiMc);
	Float_t energyMc[MAXARRAY];it->second->SetBranchAddress("energyMc",energyMc);
	
	Int_t           nAK5PFPUcorrJet;                                  
	//Int_t           chargeAK5PFPUcorrJet[300];
	Float_t         energyAK5PFPUcorrJet[300];
	//Float_t         thetaAK5PFPUcorrJet[300];
	Float_t         etaAK5PFPUcorrJet[300];
	Float_t         phiAK5PFPUcorrJet[300];
	Float_t         pxAK5PFPUcorrJet[300];
	Float_t         pyAK5PFPUcorrJet[300];
	Float_t         pzAK5PFPUcorrJet[300];
	Float_t         chargedHadronEnergyAK5PFPUcorrJet[300];
	Float_t         neutralHadronEnergyAK5PFPUcorrJet[300];
	Float_t         photonEnergyAK5PFPUcorrJet[300];
	Float_t         electronEnergyAK5PFPUcorrJet[300];
	Float_t         muonEnergyAK5PFPUcorrJet[300];
	Float_t         HFHadronEnergyAK5PFPUcorrJet[300];
	Float_t         HFEMEnergyAK5PFPUcorrJet[300];
	Int_t           chargedHadronMultiplicityAK5PFPUcorrJet[300];
	Int_t           neutralHadronMultiplicityAK5PFPUcorrJet[300];
	Int_t           photonMultiplicityAK5PFPUcorrJet[300];
	Int_t           electronMultiplicityAK5PFPUcorrJet[300];
	Int_t           muonMultiplicityAK5PFPUcorrJet[300];
	Int_t           HFHadronMultiplicityAK5PFPUcorrJet[300];
	Int_t           HFEMMultiplicityAK5PFPUcorrJet[300];
	Float_t         chargedEmEnergyAK5PFPUcorrJet[300];
	Float_t         neutralEmEnergyAK5PFPUcorrJet[300];
	
it->second->SetBranchAddress("nAK5PFPUcorrJet"				,&nAK5PFPUcorrJet);
//it->second->SetBranchAddress("chargeAK5PFPUcorrJet"			,chargeAK5PFPUcorrJet);
it->second->SetBranchAddress("energyAK5PFPUcorrJet"			,energyAK5PFPUcorrJet);
//it->second->SetBranchAddress("thetaAK5PFPUcorrJet"			,thetaAK5PFPUcorrJet);
it->second->SetBranchAddress("etaAK5PFPUcorrJet"			,etaAK5PFPUcorrJet);
it->second->SetBranchAddress("phiAK5PFPUcorrJet"			,phiAK5PFPUcorrJet);
it->second->SetBranchAddress("pxAK5PFPUcorrJet"				,pxAK5PFPUcorrJet);
it->second->SetBranchAddress("pyAK5PFPUcorrJet"				,pyAK5PFPUcorrJet);
it->second->SetBranchAddress("pzAK5PFPUcorrJet"				,pzAK5PFPUcorrJet);
it->second->SetBranchAddress("chargedHadronEnergyAK5PFPUcorrJet"	,chargedHadronEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("neutralHadronEnergyAK5PFPUcorrJet"	,neutralHadronEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("photonEnergyAK5PFPUcorrJet"		,photonEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("electronEnergyAK5PFPUcorrJet"		,electronEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("muonEnergyAK5PFPUcorrJet"			,muonEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("HFHadronEnergyAK5PFPUcorrJet"		,HFHadronEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("HFEMEnergyAK5PFPUcorrJet"			,HFEMEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("chargedHadronMultiplicityAK5PFPUcorrJet"	,chargedHadronMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("neutralHadronMultiplicityAK5PFPUcorrJet"	,neutralHadronMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("photonMultiplicityAK5PFPUcorrJet"		,photonMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("electronMultiplicityAK5PFPUcorrJet"	,electronMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("muonMultiplicityAK5PFPUcorrJet"		,muonMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("HFHadronMultiplicityAK5PFPUcorrJet"	,HFHadronMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("HFEMMultiplicityAK5PFPUcorrJet"		,HFEMMultiplicityAK5PFPUcorrJet);
it->second->SetBranchAddress("chargedEmEnergyAK5PFPUcorrJet"		,chargedEmEnergyAK5PFPUcorrJet);
it->second->SetBranchAddress("neutralEmEnergyAK5PFPUcorrJet"		,neutralEmEnergyAK5PFPUcorrJet);


long long nEntries=it->second->GetEntries();
double weight=xSec[it->first]/nEntries;

long long entryProcess=nEntries/nJobs + 1;
long long entryMin=job*entryProcess;
long long entryMax=TMath::Min(nEntries,(job+1)*entryProcess);

cout<<"nEntries="<<nEntries<<" xSec="<<xSec[it->first]<<" weight="<<weight<<endl;
cout<<"entryMin="<<entryMin<<" entryMax="<<entryMax<< "job="<<job<<" nJobs="<<nJobs<<endl;

for(int iEntry=entryMin;iEntry<entryMax;iEntry++)
{ 	
	it->second->GetEntry(iEntry);
	//find the first jet ided
	int Jet[3];
	int JetMc[3];
	int nJet=0; 

	int iJet=0; 
	for(iJet=0;iJet<nAK5PFPUcorrJet;iJet++)
	{
	//jet Id
	float phf= photonEnergyAK5PFPUcorrJet[iJet]/energyAK5PFPUcorrJet[iJet];
	float chf= chargedHadronEnergyAK5PFPUcorrJet[iJet]/energyAK5PFPUcorrJet[iJet];
	float nhf= neutralHadronEnergyAK5PFPUcorrJet[iJet]/energyAK5PFPUcorrJet[iJet];
	float elf= electronEnergyAK5PFPUcorrJet[iJet]/energyAK5PFPUcorrJet[iJet];
	int   chm= chargedHadronMultiplicityAK5PFPUcorrJet[iJet];
	int   npr= chargedHadronMultiplicityAK5PFPUcorrJet[iJet] + neutralHadronMultiplicityAK5PFPUcorrJet[iJet] + photonMultiplicityAK5PFPUcorrJet[iJet] + electronMultiplicityAK5PFPUcorrJet[iJet] + muonMultiplicityAK5PFPUcorrJet[iJet] + HFHadronMultiplicityAK5PFPUcorrJet[iJet] + HFEMMultiplicityAK5PFPUcorrJet[iJet] ;
	float eta= etaAK5PFPUcorrJet[iJet];
	bool id=(npr>1 && phf<0.99 && nhf<0.99 && ((fabs(eta)<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0) || fabs(eta)>2.4));
		
	if (!id ) continue; //require id for the first jet
	//matching parton	
	TLorentzVector j;
	j.SetPxPyPzE(pxAK5PFPUcorrJet[iJet],pyAK5PFPUcorrJet[iJet],pzAK5PFPUcorrJet[iJet],energyAK5PFPUcorrJet[iJet]);
	int i_foundPart=-1;
       	float bestDeltaRPart=.3;
       	for( int iMC=0; iMC<nMc; ++iMC ) {

         if( statusMc[iMC]!=3 ) continue;
         if( !(fabs(idMc[iMC])<7 || idMc[iMC]==21) ) continue;

         TLorentzVector thisPart;
         thisPart.SetPtEtaPhiE( energyMc[iMC]*sin(thetaMc[iMC]), etaMc[iMC], phiMc[iMC], energyMc[iMC] );
         if( thisPart.Pt()<0.1 ) continue;


         float thisDeltaR = thisPart.DeltaR(j);

         if( thisDeltaR<bestDeltaRPart ) {
           bestDeltaRPart = thisDeltaR;
           i_foundPart = iMC;
         }

       } //for partons
	
	Jet[nJet]=iJet;	
	JetMc[nJet]=i_foundPart;
	nJet++;
	if(nJet>=4) break; //nJet starts from 0
	} //Jet Loop
	
	

	bool dijet=false;
	if(nJet<2) continue;
	
	TLorentzVector j1,j2,j3;
	j1.SetPxPyPzE(pxAK5PFPUcorrJet[Jet[0]],pyAK5PFPUcorrJet[Jet[0]],pzAK5PFPUcorrJet[Jet[0]],energyAK5PFPUcorrJet[Jet[0]]);
	j2.SetPxPyPzE(pxAK5PFPUcorrJet[Jet[1]],pyAK5PFPUcorrJet[Jet[1]],pzAK5PFPUcorrJet[Jet[1]],energyAK5PFPUcorrJet[Jet[1]]);
	if(nJet>2)j3.SetPxPyPzE(pxAK5PFPUcorrJet[Jet[2]],pyAK5PFPUcorrJet[Jet[2]],pzAK5PFPUcorrJet[Jet[2]],energyAK5PFPUcorrJet[Jet[2]]);
	dijet= (fabs(j1.DeltaPhi(j2)>2.5)) ;
	
	dijet= dijet && ((nJet==2) || ((  j3.Pt()/(j1.Pt()+j2.Pt())  ) < .10 ) || (j3.Pt()<10)  ); //20% of the mean

	if(!dijet) continue;	

	bool isQ=false, isG=false;	
	if((JetMc[0]>0) && ( fabs(idMc[JetMc[0]])<5   ) && (idMc[JetMc[0]]!=0) ) isQ=true;
	if((JetMc[0]>0) && ( fabs(idMc[JetMc[0]])==21 ) && (idMc[JetMc[0]]!=0) ) isG=true;
	//exclude the ones with low weight
	//Double_t Bins[]={30,40,50,60,70,80,100,120,140,160,180,200,250,300,350,400,500,600,700,800,1000,1200,1400,1600,1800,2000};
	bool fill=true;
	for(int bin=0;bin<20;bin++)
		{
		if( Bins[bin] > j1.Pt() || j1.Pt()> Bins[bin+1])  continue; //not in the right bin
		if( energyMc[JetMc[0]]*sin(thetaMc[JetMc[0]]) < Bins[bin] ) fill=false;
		}
	if (!fill) continue;

	if( fabs(j1.Eta())<2.0){	
	if(isQ) PtJet_Q->Fill(j1.Pt(),weight);
	else if(isG) PtJet_G->Fill(j1.Pt(),weight);
	else PtJet_O->Fill(j1.Pt(),weight);
	}
	
	if( 3.0<fabs(j1.Eta()) && fabs(j1.Eta() <4.7) )
	{
	if(isQ) PtJetFwd_Q->Fill(j1.Pt(),weight);
	else if(isG) PtJetFwd_G->Fill(j1.Pt(),weight);
	else PtJetFwd_O->Fill(j1.Pt(),weight);
	}
	
	//cout<<"Status leading Jet="<<idMc[JetMc[0]]<<" SubLeading="<<idMc[JetMc[1]]<<endl;
} //loop on all entries
}//loop on all mc

PtJet_Q->SetMarkerStyle(20);PtJet_Q->SetMarkerColor(kBlue-4);
PtJet_G->SetMarkerStyle(20);PtJet_G->SetMarkerColor(kRed-4);
PtJet_O->SetMarkerStyle(20);PtJet_O->SetMarkerColor(kGreen-4);
PtJet_G->GetXaxis()->SetTitle("P_{T} Jet");
PtJet_G->GetYaxis()->SetTitle("fraction");

PtJetFwd_Q->SetMarkerStyle(20);PtJetFwd_Q->SetMarkerColor(kBlue-4);
PtJetFwd_G->SetMarkerStyle(20);PtJetFwd_G->SetMarkerColor(kRed-4);
PtJetFwd_O->SetMarkerStyle(20);PtJetFwd_O->SetMarkerColor(kGreen-4);
PtJetFwd_G->GetXaxis()->SetTitle("P_{T} Jet");
PtJetFwd_G->GetYaxis()->SetTitle("fraction");


PtJet_Q->Write();
PtJet_G->Write();
PtJet_O->Write();

PtJetFwd_Q->Write();
PtJetFwd_G->Write();
PtJetFwd_O->Write();

TH1F* tot=(TH1F*)PtJet_Q->Clone("tot"); tot->Add(PtJet_G); tot->Add(PtJet_O);
TH1F* r_q =(TH1F*)PtJet_Q->Clone("Qfrac"); r_q->Divide(tot);
TH1F* r_g =(TH1F*)PtJet_G->Clone("Gfrac"); r_g->Divide(tot);
TH1F* r_o =(TH1F*)PtJet_O->Clone("Ofrac"); r_o->Divide(tot);

TCanvas *c=new TCanvas("fraction","fraction",800,800);
r_q->SetMarkerStyle(20);r_q->SetMarkerColor(kBlue-4);
r_g->SetMarkerStyle(20);r_g->SetMarkerColor(kRed-4);
r_o->SetMarkerStyle(20);r_o->SetMarkerColor(kGreen-4);
r_g->GetXaxis()->SetTitle("P_{T} Jet");
r_g->GetYaxis()->SetTitle("fraction");
r_g->Draw("P");
r_q->Draw("P SAME");
r_o->Draw("P SAME");

if(nJobs==1)c->Write();
system(Form("touch Plot/done_%d_%d.txt",job,nJobs));
} 


