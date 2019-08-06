#include <TGrid.h>
void DoMerge(Int_t nrun=265344){

  //  TString path=Form("/alice/cern.ch/user/e/epiga/KLambdapPb_MC44/output"); //High intensity run v2
     TString path=Form("/alice/cern.ch/user/c/chdemart/K0sjetpp_LHC16l_10runs_3rdtry/output/"); //High intensity run v2

  TString pattern="AnalysisResults.root";

  TString outfile="AnalysisResults2016l_10runs_3rdtry.";
  //outfile+=nrun;
  outfile+="root";

  cout << "I'm going to merge the following files:" << endl;
  cout << path << " \n into this file:\n" << outfile << endl;
 
  gSystem->Load("libVMC.so");
  printf("*** Connect to GRID ***\n");
  

  if (!gGrid) {
    TGrid::Connect("alien://");
  }

  TGridResult* result = gGrid->Query(path,pattern);
  result->Print();
  
  TFileMerger m;
  
  if (outfile) m.OutputFile(outfile);
  Int_t i=0;
  while (result->GetKey(i,"turl")) {
    m.AddFile(result->GetKey(i,"turl"));
    i++;
  }
  if (i)
    m.Merge();
}

