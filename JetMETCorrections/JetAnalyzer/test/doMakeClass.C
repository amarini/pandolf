void doMakeClass() {

  TChain* chain = new TChain("pippo");
 /* 
  for(int i=1; i<462; ++i) {

    char fileName[70];

    sprintf(fileName, "JetsEt20_allPF/output_%d.root", i);

    chain->Add(fileName);

  }
*/

  chain->Add("QCD_2_2_3/output_147.root/myanalysis/pippo");
  chain->MakeClass("ProvaAnalyzer");

}
