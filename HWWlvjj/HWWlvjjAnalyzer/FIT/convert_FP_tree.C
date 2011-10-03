void convert_FP_tree(){

  // TFile *dfile = new TFile("/afs/cern.ch/cms/CAF/CMSALCA/ALCA_TRACKERALIGN/HIP/bonato/HtoZZto2L2J-crab/trees_FP_Summer11/HZZlljjRM_DATA_800pb_optLD_looseBTags_v2_ALL_FP.root");
  // TFile *dfile = new TFile("/afs/cern.ch/user/w/whitbeck/scratch0/HZZlljjRM_DATA_1fb_optLD_looseBTags_v2_ALL.root");//859
 TFile *dfile = new TFile("HWWlvjj_DATA_6july_helicity_ALL.root","READ");//1000
 TTree *t=(TTree*)dfile->Get("Tree_FITUL");
 cout<<"Loaded tree "<<t->GetName()<<endl;
 int leptTypeIN;
 float mWWIN, mJJIN,eventWeightIN;
 cout<<"addressing"<<endl;
 t->SetBranchAddress("mWW", &mWWIN);
 t->SetBranchAddress("mJJ", &mJJIN);
 t->SetBranchAddress("leptType", &leptTypeIN);
 t->SetBranchAddress("eventWeight", &eventWeightIN);


 double leptTypeOUT;
 double mWWOUT, mJJOUT,eventWeightOUT;
 TFile *fout=new TFile("./convertedTree_PROVA.root","RECREATE");
 TTree *tout=new TTree("Tree_FITUL","Converted from FP");
 tout->Branch("CMS_hwwlvqq_mWW",&mWWOUT,"CMS_hwwlvqq_mWW/D");
 tout->Branch("mJJ",&mJJOUT,"mJJ/D");
 tout->Branch("leptType",&leptTypeOUT,"leptType/D");
 tout->Branch("eventWeight",&eventWeightOUT,"eventWeight/D");

 cout<<"start"<<endl;
 for(int i=0;i<t->GetEntries();i++){
   t->GetEntry(i);
   leptTypeOUT=double(leptTypeIN);
   mWWOUT=double(mWWIN);
   mJJOUT=double(mJJIN);
   eventWeightOUT=double(eventWeightIN);
   tout->Fill();
 }
 cout<<"finished"<<endl;
 tout->Write();
 delete tout;
 delete fout;

}
