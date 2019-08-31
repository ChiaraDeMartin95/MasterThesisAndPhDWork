#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TFile.h>

void ACDiffPeriodCheck( Bool_t isMC=0, Bool_t isEfficiency=0,Int_t sysTrigger=0, Int_t sysV0=5, Int_t type=0){

  TString PathInBis;
 TString PathOut;
  TFile *fileinbis;
 

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;

  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};
  TString SPeriod[numPeriod]={"2016", "ALL"};
  TString SRun[numPeriod]={"2016k", "All"};

  Int_t Marker[numSysV0]={7,4,20,22,29, 35};
  Int_t Color[numSysV0]={2,3,4,6,7,9, 27};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};

 
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassME[nummolt+1][numzeta][numPtV0][numPtTrigger];

  TH1D *hInvMassK0Short_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_SEbins_master[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_SEbins_ratio[nummolt+1][numzeta][numPtV0][numPtTrigger][numPeriod];
  TH1D *hInvMassK0Short_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_MEbins_master[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_MEbins_ratio[nummolt+1][numzeta][numPtV0][numPtTrigger][numPeriod];

 TH2D *hMassvsPt_SEbins[nummolt+1][numzeta][numPtTrigger];
 TH2D *hMassvsPt_SEbins_master[nummolt+1][numzeta][numPtTrigger];
 TH2D *hMassvsPt_SEbins_ratio[nummolt+1][numzeta][numPtTrigger][numPeriod];
 TH2D *hMassvsPt_MEbins[nummolt+1][numzeta][numPtTrigger];
 TH2D *hMassvsPt_MEbins_master[nummolt+1][numzeta][numPtTrigger];
 TH2D *hMassvsPt_MEbins_ratio[nummolt+1][numzeta][numPtTrigger][numPeriod];


 TH2D *hDeltaEtaDeltaPhi_SEbins_before[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_SEbins_master_before[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_SEbins_ratio_before[nummolt+1][numzeta][numPtV0][numPtTrigger][numPeriod];
 TH2D *hDeltaEtaDeltaPhi_MEbins_before[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_MEbins_master_before[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_MEbins_ratio_before[nummolt+1][numzeta][numPtV0][numPtTrigger][numPeriod];

 
 TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_SEbins_master[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_SEbins_ratio[nummolt+1][numzeta][numPtV0][numPtTrigger][numPeriod];
 TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_MEbins_master[nummolt+1][numzeta][numPtV0][numPtTrigger];
 TH2D *hDeltaEtaDeltaPhi_MEbins_ratio[nummolt+1][numzeta][numPtV0][numPtTrigger][numPeriod];

   PathOut="FinalOutput/histo/ACDiffPeriodCheck";
	 
	  if(isMC && isEfficiency){
	 
	    PathOut+="_MCEff";
	  }
	  if(isMC && !isEfficiency){
	 
	    PathOut+="_MCTruth";
	  }

	  TString PathOut2=PathOut;
	  TString PathOut3=PathOut;	 
	  PathOut +=Form("MassDistr_SysT%i_SysV0%i.root",sysTrigger, sysV0); 

	  PathOut2+=Form("ACBeforeCuts_SysT%i_SysV0%i.root",sysTrigger, sysV0); 
	  PathOut3+=Form("ACAfterCuts_SysT%i_SysV0%i.root",sysTrigger, sysV0); 

	  TFile *fileOut = new TFile(PathOut, "RECREATE"); 
	  TFile *fileOut2 = new TFile(PathOut2, "RECREATE"); 
	  TFile *fileOut3 = new TFile(PathOut3, "RECREATE"); 

  //Form("hMassvsPt_"+tipo[type]+"_%i",molt
  for(Int_t m=0; m<nummolt; m++){
    cout << "\n\n molt " << m << endl;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t p=0; p<numPeriod; p++){
	  TString PathIn="FinalOutput/DATA";
	  TString PathIn2;
	  PathIn+=SPeriod[p];
	  PathIn +="/histo/AngularCorrelation";
	  PathIn+=SRun[p];  
	  if(isMC && isEfficiency){
	    PathIn+="_MCEff";
	 
	  }
	  if(isMC && !isEfficiency){
	    PathIn+="_MCTruth";
	
	  }
	  PathIn2=PathIn;
	  PathIn+=Form("_MassDistr_SysT%i_SysV0%i.root",sysTrigger, sysV0); 
	  PathIn2+=Form("_SysT%i_SysV0%i.root",sysTrigger, sysV0); 
	

	  TFile *filein=new TFile(PathIn, "");
	  TFile *filein2=new TFile(PathIn2, "");
	  cout << "\n\n" << PathIn << endl;
	  cout << PathIn2 << endl;

	  if (p==0){
	    hMassvsPt_SEbins_master[m][z][tr]=(TH2D*)filein->Get(Form("SE_hMassvsPt_"+tipo[type]+"_%i", m));
	    hMassvsPt_MEbins_master[m][z][tr]=(TH2D*)filein->Get(Form("ME_hMassvsPt_"+tipo[type]+"_%i", m));
	    
	   	  
	  }
	  else {
	    hMassvsPt_SEbins[m][z][tr]=(TH2D*)filein->Get(Form("SE_hMassvsPt_"+tipo[type]+"_%i", m));
	    hMassvsPt_MEbins[m][z][tr]=(TH2D*)filein->Get(Form("ME_hMassvsPt_"+tipo[type]+"_%i", m));	   
	 
	    hMassvsPt_SEbins_ratio[m][z][tr][p]=(TH2D*)hMassvsPt_SEbins[m][z][tr]->Clone(Form("SE_hMassvsPt_ratio"+tipo[type]+"_%i", m));
	    hMassvsPt_SEbins_ratio[m][z][tr][p]->Divide(  hMassvsPt_SEbins[m][z][tr],   hMassvsPt_SEbins_master[m][z][tr]);
	   
	    hMassvsPt_MEbins_ratio[m][z][tr][p]=(TH2D*)hMassvsPt_MEbins[m][z][tr]->Clone(Form("ME_hMassvsPt_"+tipo[type]+"_%i", m));
	    hMassvsPt_MEbins_ratio[m][z][tr][p]->Divide(  hMassvsPt_MEbins[m][z][tr],   hMassvsPt_MEbins_master[m][z][tr]);

	   
	    fileOut->cd();
	    hMassvsPt_SEbins_ratio[m][z][tr][p]->Write();
	    hMassvsPt_MEbins_ratio[m][z][tr][p]->Write();
	  }

	  cout << "ora distrb. 1D" << endl;

	  for(Int_t v=0; v<numPtV0; v++){
	    nameSE[m][z][v][tr]="SE_";
	    namemassSE[m][z][v][tr]="InvMassSE_";
	    nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	    namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	    nameME[m][z][v][tr]="ME_";
	    namemassME[m][z][v][tr]="InvMassME_";
	    nameME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	    namemassME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	   	   
	    if(p==0){
	    hInvMassK0Short_SEbins_master[m][z][v][tr]= (TH1D*)filein->Get(namemassSE[m][z][v][tr]);
	    hInvMassK0Short_SEbins_master[m][z][v][tr]->Sumw2();
	    hInvMassK0Short_MEbins_master[m][z][v][tr]= (TH1D*)filein->Get(namemassME[m][z][v][tr]);
	    hInvMassK0Short_MEbins_master[m][z][v][tr]->Sumw2();
	    hDeltaEtaDeltaPhi_SEbins_master_before[m][z][v][tr]=(TH2D*)filein->Get("SE_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_MEbins_master_before[m][z][v][tr]=(TH2D*)filein->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_SEbins_master[m][z][v][tr]=(TH2D*)filein2->Get("SE_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_MEbins_master[m][z][v][tr]=(TH2D*)filein2->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[v]);
	  
	    }
	    else {
	    hInvMassK0Short_SEbins[m][z][v][tr]= (TH1D*)filein->Get(namemassSE[m][z][v][tr]);
	    hInvMassK0Short_SEbins[m][z][v][tr]->Sumw2();
	    hInvMassK0Short_SEbins_ratio[m][z][v][tr][p]= (TH1D*) hInvMassK0Short_SEbins[m][z][v][tr]->Clone(namemassSE[m][z][v][tr] + "_ratio");
	    hInvMassK0Short_SEbins_ratio[m][z][v][tr][p]-> Divide(hInvMassK0Short_SEbins[m][z][v][tr], hInvMassK0Short_SEbins_master[m][z][v][tr]);
	    hInvMassK0Short_MEbins[m][z][v][tr]= (TH1D*)filein->Get(namemassME[m][z][v][tr]);
	    hInvMassK0Short_MEbins[m][z][v][tr]->Sumw2();
	    hInvMassK0Short_MEbins_ratio[m][z][v][tr][p]= (TH1D*) hInvMassK0Short_MEbins[m][z][v][tr]->Clone(namemassME[m][z][v][tr] + "_ratio");
	    hInvMassK0Short_MEbins_ratio[m][z][v][tr][p]-> Divide(hInvMassK0Short_MEbins[m][z][v][tr], hInvMassK0Short_MEbins_master[m][z][v][tr]);
	    cout << "ok " << endl;


	    hDeltaEtaDeltaPhi_SEbins_before[m][z][v][tr]=(TH2D*)filein->Get("SE_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_MEbins_before[m][z][v][tr]=(TH2D*)filein->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]=(TH2D*)filein2->Get("SE_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]=(TH2D*)filein2->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    cout << "ok " << endl;
	    hDeltaEtaDeltaPhi_SEbins_ratio_before[m][z][v][tr][p]=(TH2D*) hDeltaEtaDeltaPhi_SEbins_before[m][z][v][tr]->Clone("SE_ratio_before_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    cout << "ok " << endl;
	    hDeltaEtaDeltaPhi_MEbins_ratio_before[m][z][v][tr][p]=(TH2D*) hDeltaEtaDeltaPhi_MEbins_before[m][z][v][tr]->Clone("ME_ratio_before_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    cout << "ok " << endl;
	    hDeltaEtaDeltaPhi_SEbins_ratio[m][z][v][tr][p]=(TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone("SE_ratio_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    cout << "ok " << endl;
	    hDeltaEtaDeltaPhi_MEbins_ratio[m][z][v][tr][p]=(TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Clone("ME_ratio_m"+ Smolt[m]+"_v"+SPtV0[v]);
	    cout << "dividing " << endl;

	    hDeltaEtaDeltaPhi_SEbins_ratio_before[m][z][v][tr][p]->Divide(hDeltaEtaDeltaPhi_SEbins_before[m][z][v][tr], 	    hDeltaEtaDeltaPhi_SEbins_master_before[m][z][v][tr]);
	    hDeltaEtaDeltaPhi_MEbins_ratio_before[m][z][v][tr][p]->Divide(	    hDeltaEtaDeltaPhi_MEbins_before[m][z][v][tr], 	    hDeltaEtaDeltaPhi_MEbins_master_before[m][z][v][tr]);
	    hDeltaEtaDeltaPhi_SEbins_ratio[m][z][v][tr][p]->Divide(	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr], 	    hDeltaEtaDeltaPhi_SEbins_master[m][z][v][tr]);
	    hDeltaEtaDeltaPhi_MEbins_ratio[m][z][v][tr][p]->Divide(	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr], 	    hDeltaEtaDeltaPhi_MEbins_master[m][z][v][tr]);

	    cout << "writing in file " << endl;
	    fileOut->cd();
	    hInvMassK0Short_SEbins_ratio[m][z][v][tr][p]->Write();
	    hInvMassK0Short_MEbins_ratio[m][z][v][tr][p]->Write();
	    hDeltaEtaDeltaPhi_SEbins_ratio_before[m][z][v][tr][p]->Write();
	    hDeltaEtaDeltaPhi_MEbins_ratio_before[m][z][v][tr][p]->Write();
	    hDeltaEtaDeltaPhi_SEbins_ratio[m][z][v][tr][p]->Write();
	    hDeltaEtaDeltaPhi_MEbins_ratio[m][z][v][tr][p]->Write();
	  }
	    
	  }
	}
      }
    }
  }

 cout << " ho creato il file " << PathOut<< endl;

}

