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
#include <TSpline.h>

void CfrDataMCRatio(Float_t PtTrigMin=3.0, Int_t isDataOrMC=2, TString year0="2016", Int_t ishhCorr=0, Bool_t MasterThesis=0, Bool_t isSE=0, Int_t NumberOfRegions=6,   Float_t ScaleFactorJet = 0.84/2.08*2.13628/2.63894, Float_t ScaleFactorPions=1./0.8767, Bool_t EtaPhiEnlargedStudy=1 ){
  //if isSE is true the ratio of K0s to hadrons is computed, in order to study the strangeness enhancement effect

  if (!EtaPhiEnlargedStudy) {
    cout << "you are using the K0s results obtained from a reduced DeltaEta and DeltaPhi regions with respect to the regions used in hhCorrelation. Are you sure you want to continue?" << endl;
    //    return;
  }

  const Int_t NumberTypeAnalysis=6; 
  Int_t WhatKindOf=0;
  //  Float_t ScaleFactorJet = 0.84/2.08; //DeltaEta(hK0s) =~ 0.8; DeltaEta(hh) =~ 2; The TOTAL in-jet production for hh and hK0s should be compared: this is why we need this scale factor!
  //ScaleFactorPions mi dice la percentuale di pioni sul totale di particelle associate: è l'84% se pT, Trig, min > 3 GeV/c (questo per hhCorr, per hhCorr_New è 87.7%)

  if (NumberOfRegions>6) {cout << "There are only 6 different region " << endl; return;}
  if (!ishhCorr && !isSE) WhatKindOf=0;
  if (ishhCorr) WhatKindOf=1;
  if (!ishhCorr && isSE) WhatKindOf=2;

  Bool_t isSEBis=0;
  if (ishhCorr==0) isSEBis=0;
  if (ishhCorr==1) isSEBis=0;
  if (isSE==1) {ishhCorr=0; isSEBis=1;}
  cout << isSE << endl;
  Int_t MinLoop=0;
  Int_t MaxLoop=0;
  if (ishhCorr) {MinLoop=1; MaxLoop=2;}
  else if (isSE) {MinLoop=0; MaxLoop=2;}
  else {MinLoop=0; MaxLoop=1;}

  if (isDataOrMC > 2){
    cout << "is DataOrMC 0 for data, 1 for MC, 2 for both " << endl;
    return;
  }
  Int_t DataOrMCMin=0;
  Int_t DataOrMCMax=0;
  if (isDataOrMC==0) {DataOrMCMin=0; DataOrMCMax=1;}
  else   if (isDataOrMC==1) {DataOrMCMin=1; DataOrMCMax=2;} 
  else {DataOrMCMin=0; DataOrMCMax=2;}

  TString PathIn[2];
  TString PathOut;
  TFile *filein[2];
  TFile *fileout;

  TString ishhCorrOrNot[3]= {"", "_hhCorr_New", "SE_New"}; 
  TString nomefileoutput = "CfrDataMCRatioOutput/CfrDataMCRatioOutput" + ishhCorrOrNot[ishhCorr]+ Form("_PtTrigMin%.1f_NumberOfRegions%i", PtTrigMin, NumberOfRegions);
  if (isSEBis==1) nomefileoutput ="CfrDataMCRatioOutput/CfrDataMCRatioOutput_" + ishhCorrOrNot[2]+Form("_PtTrigMin%.1f_NumberOfRegions%i", PtTrigMin, NumberOfRegions);
  if (EtaPhiEnlargedStudy){
    nomefileoutput = "CfrDataMCRatioOutput/CfrDataMCRatioOutput" + ishhCorrOrNot[ishhCorr]+ Form("_PtTrigMin%.1f_DeltaEtaPhiEnlarged_NumberOfRegions%i", PtTrigMin, NumberOfRegions);
    if (isSEBis==1) nomefileoutput ="CfrDataMCRatioOutput/CfrDataMCRatioOutput_" + ishhCorrOrNot[2]+Form("_PtTrigMin%.1f_DeltaEtaPhiEnlarged_NumberOfRegions%i", PtTrigMin, NumberOfRegions);
  }

  TString  nomefileoutputPDF=nomefileoutput;
  nomefileoutput += ".root";
  cout << "nome file output " << nomefileoutput << endl;
  fileout = new TFile (nomefileoutput, "RECREATE");
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;
  //colori e marker diversi per jet, OJ, J+OJ *********************************
  Float_t Color[NumberTypeAnalysis]={1,801,2, 909,881,  868};
  if (NumberOfRegions==3) {
    Color[0]=2;
    Color[1]=860;
    Color[2]=418;
  }
  Int_t MarkerOrigin[NumberTypeAnalysis]={20,21,33, 20, 21, 33}; 
  Int_t MarkerMC[NumberTypeAnalysis]={24,25,27, 24, 25, 27}; 
  //   Int_t MarkerMC[NumberTypeAnalysis]={20,21,33,20,21,33}; 
  TString Sys[2*NumberTypeAnalysis]={"Jet DATA", "Jet MC", "OJ DATA", "OJ MC", "J+OJ DATA", "J+OJ MC",  "J+OJ NotScaled DATA", "J+OJ NotScaled MC", "AwaySide DATA", "AwaySide MC", "AllButJet DATA", "AllButJet MC"};
  TString StatErrBis[NumberTypeAnalysis]={"stat. in-jet", "stat. out-of-jet", "stat. inclusive", "stat. inclusive not scaled", "stat. away side","stat. all but jet"};  
  TString SysErrBis[NumberTypeAnalysis]={"sist. in-jet", "sist. out-of-jet", "sist. inclusive",  "stat. inclusive not scaled", "stat. away si", "stat. all but jet"};
  TString Region[NumberTypeAnalysis]={"Jet", "Bulk", "All","AllNS", "AwaySide", "AllButJet"};
  TString JetOrNot[NumberTypeAnalysis]={"Jet ", "OJ " , "J+OJ ", "J+OJ NS", "AwaySide", "AllButJet"};
  TString SysErr[2]={"syst. DATA", "syst. MC"};
  TString StatErr[2]={"stat. DATA", "stat. MC"};
  TString SysErrMB[2]={"syst. DATA 0-100 %", "syst. MC 0-100 %"};
  TString StatErrMB[2]={"stat. DATA 0-100 %", "stat. MC 0-100 %"};
  TString DATAorMC[2]={"DATA", "MC"};
  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  TString SmoltBis[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  TString SPtV0[numPtV0]={"0-1", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString SPeriod[numPeriod]={"2016", "2017"};//, "2018"};
  TString SRun[numPeriod]={"2016k", "2017h"};//, "2018m"};
  //20,21,33
  Int_t Marker[6]={20,21,20,21, 20, 21};

  Int_t MarkerBis[2]={23,24};
  Int_t MarkerTris[2]={20, 33};
  Int_t MarkerMBStat[6]={20,4,20,4, 20, 4};
  Int_t MarkerMBSist[6]={21,25,21,25, 21, 25};
  TF1* pol0[2][NumberOfRegions];
  TF1* pol1[2][NumberOfRegions];

  Double_t   mult[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};

  TH1D*  fHistYieldvsErrSoloStat[6][6];
  TH1D*  fHistYieldvsErrSoloSist[6][6];
  TH1D*  fHistYieldvsErrSoloStatMB[6][6];
  TH1D*  fHistYieldvsErrSoloSistMB[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatio[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatio[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatioMB[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatioMB[6][6];

  TH1D*  fHistYieldvsErrErrSoloStat[6][6];
  TH1D*  fHistYieldvsErrErrSoloSist[6][6];
  TH1D*  fHistYieldvsErrErrSoloStatMB[6][6];
  TH1D*  fHistYieldvsErrErrSoloSistMB[6][6];

  TH1D*  fHistYieldvsErrSoloStat_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSist_SE[6][6];
  TH1D*  fHistYieldvsErrSoloStatMB_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSistMB_SE[6][6];

  TH1D*  fHistYieldvsErrSoloStatRatio_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatio_SE[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatioMB_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatioMB_SE[6][6];

  TH1D*  fHistYieldvsErrErrSoloStat_SE[6][6];
  TH1D*  fHistYieldvsErrErrSoloSist_SE[6][6];
  TH1D*  fHistYieldvsErrErrSoloStatMB_SE[6][6];
  TH1D*  fHistYieldvsErrErrSoloSistMB_SE[6][6];

  TGraph*  fGraphYieldDatiPubblicati;
  TH1D*    fHistYieldDatiPubblicati;
  TH1D*    fHistYieldErroriDatiPubblicati;
  TH1D*    fHistYieldDatiPubblicatiSE;
  TH1D*    fHistYieldDatiPubblicatiSE13TeV;
  TH1D*    fHistYieldErroriDatiPubblicatiSE;
  TH1D*    fHistYieldErroriDatiPubblicatiSE13TeV;

  TFile *filedatipubbl;
  TFile *filedatipubblSE;
  TFile *filedatipubblSE13TeV;
  TDirectoryFile *dir;
  TDirectoryFile *dirSE;
  TDirectoryFile *dirSE13TeV;
  filedatipubbl= new TFile("HEPData-1574358449-v1-Table_8c.root", "");
  filedatipubblSE= new TFile("HEPData-ins1471838-v1-Table_36.root", "");
  filedatipubblSE13TeV= new TFile("HEPData-1584713694-v1-Table_18.root", "");
  if (!filedatipubbl || !filedatipubblSE || !filedatipubblSE13TeV) return;
  cout << " I got the files where published data are stored " << endl;
  dir  = (TDirectoryFile*)filedatipubbl->Get("Table 8c");
  dirSE  = (TDirectoryFile*)filedatipubblSE->Get("Table 36");
  dirSE13TeV  = (TDirectoryFile*)filedatipubblSE13TeV->Get("Table 18");
  if (!dir || !dirSE || !dirSE13TeV) return;
  cout << " I got the directory where published data are stored " << endl;
  fGraphYieldDatiPubblicati=(TGraph*)dir->Get("Graph1D_y1");
  fHistYieldDatiPubblicati=(TH1D*)dir->Get("Hist1D_y1");
  fHistYieldErroriDatiPubblicati=(TH1D*)dir->Get("Hist1D_y1_e2"); //this is only the total systematic error (the statistic is 1/10 and the uncorrelated systematic is approx half)
  fHistYieldDatiPubblicatiSE=(TH1D*)dirSE->Get("Hist1D_y1");
  fHistYieldErroriDatiPubblicatiSE=(TH1D*)dirSE->Get("Hist1D_y1_e3"); //uncorrelated systematic error (the statistical error is ~10 times smaller
  fHistYieldDatiPubblicatiSE13TeV=(TH1D*)dirSE13TeV->Get("Hist1D_y1");
  fHistYieldErroriDatiPubblicatiSE13TeV=(TH1D*)dirSE13TeV->Get("Hist1D_y1_e3"); //uncorrelated systematic error (the statistical error is ~10 times smaller
  cout << " I got the histos  where published data are stored " << endl;
  TLegend* legend[NumberTypeAnalysis];
  TLegend* legendtype[NumberTypeAnalysis];
  TLegend* legendbis[NumberTypeAnalysis];
  Int_t msize=2;

  auto legenderr=new TLegend(0.6,0.6, 0.9, 0.9);
  auto legenderrbis=new TLegend(0.6,0.6, 0.9, 0.9);

  TString  PathInJet[2];
  TString  PathInJetMC[2];
  TString  PathInBulk[2];
  TString  PathInBulkMC[2];
  TString  PathInAll[2];
  TString  PathInAllMC[2];
  TString  PathInAllNS[2];
  TString  PathInAllNSMC[2];
  TString  PathInAwaySide[2];
  TString  PathInAwaySideMC[2];
  TString  PathInAllButJet[2];
  TString  PathInAllButJetMC[2];

  TString  hhCorr[3]={"", "_New_hhCorr", "_SE_New"}; //"_New_hhCorr"
  cout << "I'm getting the files " << endl;
  for(Int_t ishhCorrLoop=MinLoop; ishhCorrLoop<MaxLoop; ishhCorrLoop++){
    if (!MasterThesis){

      PathInJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_JetData_PtMin%.1f.root", PtTrigMin);	 
      PathInJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_JetMC_PtMin%.1f.root", PtTrigMin);	 
      PathInBulk[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_BulkData_PtMin%.1f.root", PtTrigMin); 	 
      PathInBulkMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_BulkMC_PtMin%.1f.root", PtTrigMin);
      PathInAll[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllData_PtMin%.1f.root", PtTrigMin); 	 
      PathInAllMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllMC_PtMin%.1f.root", PtTrigMin);  

      PathInAllNS[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllNSData_PtMin%.1f.root", PtTrigMin); 	 
      PathInAllNSMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllNSMC_PtMin%.1f.root", PtTrigMin);  
      PathInAwaySide[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AwaySideData_PtMin%.1f.root", PtTrigMin); 	 
      PathInAwaySideMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AwaySideMC_PtMin%.1f.root", PtTrigMin);  
      PathInAllButJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllButJetData_PtMin%.1f.root", PtTrigMin); 	
      PathInAllButJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllButJetMC_PtMin%.1f.root", PtTrigMin);

      if (ishhCorrLoop==0 && EtaPhiEnlargedStudy){
      PathInJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_JetData_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);	 
      PathInJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_JetMC_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);	 
      PathInBulk[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_BulkData_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin); 	 
      PathInBulkMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_BulkMC_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);
      PathInAll[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllData_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin); 	 
      PathInAllMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllMC_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);  

      PathInAllNS[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllNSData_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin); 	 
      PathInAllNSMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllNSMC_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);  
      PathInAwaySide[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AwaySideData_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin); 	 
      PathInAwaySideMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AwaySideMC_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);  
      PathInAllButJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllButJetData_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin); 	
      PathInAllButJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllButJetMC_PtMin%.1f_DeltaEtaPhiEnlarged.root", PtTrigMin);

      }
    }
    //per analsisi presentata nella tesi
    else{
      if (ishhCorrLoop=1) continue;
      PathInJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_Jet.root";	 
      PathInJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_JetMC.root";	 
      PathInBulk[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_Bulk.root"; 	 
      PathInBulkMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_BulkMC.root";
      PathInAll[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_All.root"; 	 
      PathInAllMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_AllMC.root";  
    }
  }
  cout << "...Done\n\n" << endl;
  cout << "These are the input files for K0s: " << PathInJet[0]<< endl;
  cout << "I'm getting the histos " << endl;
  //**************** prendo gli istogrammi **********************************************
  for(Int_t l=0; l<NumberOfRegions; l++){ //loop sulla region (jet, OJ, inclusive)
    for(Int_t j=DataOrMCMin; j<DataOrMCMax; j++){ //loop data or MC
      cout << " region " << l << "data/MC " << j << endl;
      for (Int_t ishhCorrLoop=MinLoop; ishhCorrLoop<MaxLoop; ishhCorrLoop++){
	if(l==0){
	  if (j==0)   PathIn[ishhCorrLoop]=PathInJet[ishhCorrLoop];
	  if (j==1)   PathIn[ishhCorrLoop]=PathInJetMC[ishhCorrLoop];
	}
	if(l==1){
	  if (j==0)  PathIn[ishhCorrLoop]=PathInBulk[ishhCorrLoop];
	  if (j==1)  PathIn[ishhCorrLoop]=PathInBulkMC[ishhCorrLoop];
	}
	if (l==2){
	  if (j==0)  PathIn[ishhCorrLoop]=PathInAll[ishhCorrLoop];
	  if (j==1)  PathIn[ishhCorrLoop]=PathInAllMC[ishhCorrLoop];
	}
	if(l==3){
	  if (j==0)   PathIn[ishhCorrLoop]=PathInAllNS[ishhCorrLoop];
	  if (j==1)   PathIn[ishhCorrLoop]=PathInAllNSMC[ishhCorrLoop];
	}
	if(l==4){
	  if (j==0)   PathIn[ishhCorrLoop]=PathInAwaySide[ishhCorrLoop];
	  if (j==1)   PathIn[ishhCorrLoop]=PathInAwaySideMC[ishhCorrLoop];
	}
	if(l==5){
	  if (j==0)   PathIn[ishhCorrLoop]=PathInAllButJet[ishhCorrLoop];
	  if (j==1)   PathIn[ishhCorrLoop]=PathInAllButJetMC[ishhCorrLoop];
	}

	filein[ishhCorrLoop]=new TFile(PathIn[ishhCorrLoop], "");
      }
      fHistYieldvsErrSoloSistMB[j][l]=   (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloSistMB");   
      fHistYieldvsErrSoloStatMB[j][l]=   (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloStatMB");   
      fHistYieldvsErrSoloStat[j][l]=     (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloStat");	    
      fHistYieldvsErrSoloSist[j][l]=     (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloSist");	    
      fHistYieldvsErrSoloStatRatio[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloStatRatio");
      fHistYieldvsErrSoloSistRatio[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloSistRatio");
      fHistYieldvsErrErrSoloStat[j][l]=  (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloStat");  
      fHistYieldvsErrErrSoloSist[j][l]=  (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloSist");  
      fHistYieldvsErrErrSoloStatMB[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloStatMB");
      fHistYieldvsErrErrSoloSistMB[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloSistMB");

      if (l==0 && EtaPhiEnlargedStudy && ishhCorr==0 && !isSEBis){ //in order to compare results from EtaPhiEnlarged Study for hK0s (DeltaEtaJet < 1, DeltaPhi < 1.1) to default results (DeltaEtaJet < 0.4, DeltaPhi < 1), the production in the jet should be rescaled  
	cout << "I'm scaling " << endl;
	/*
	  fHistYieldvsErrSoloSistMB[j][l]   ->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrSoloStatMB[j][l]   ->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrSoloStat[j][l]     ->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrSoloSist[j][l]     ->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrSoloStatRatio[j][l]->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrSoloSistRatio[j][l]->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrErrSoloStat[j][l]  ->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrErrSoloSist[j][l]  ->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrErrSoloStatMB[j][l]->Scale(1./ScaleFactorJet);
	  fHistYieldvsErrErrSoloSistMB[j][l]->Scale(1./ScaleFactorJet);
	*/
	  }

      cout << "isSE " << isSEBis << endl;
      if (isSEBis){
	cout << "ishhCorr vale " << ishhCorr << endl;
	cout << "l (region ) " << l << " data or MC j " << j << " path hh" << PathIn[1]<< endl;   
	fHistYieldvsErrSoloSistMB_SE[j][l]=   (TH1D*)filein[1]->Get("fHistYieldvsErrSoloSistMB");   
	fHistYieldvsErrSoloStatMB_SE[j][l]=   (TH1D*)filein[1]->Get("fHistYieldvsErrSoloStatMB");   
	fHistYieldvsErrSoloStat_SE[j][l]=     (TH1D*)filein[1]->Get("fHistYieldvsErrSoloStat");	    
	fHistYieldvsErrSoloSist_SE[j][l]=     (TH1D*)filein[1]->Get("fHistYieldvsErrSoloSist");	    
	fHistYieldvsErrSoloStatRatio_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrSoloStatRatio");
	fHistYieldvsErrSoloSistRatio_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrSoloSistRatio");
	fHistYieldvsErrErrSoloStat_SE[j][l]=  (TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloStat");  
	fHistYieldvsErrErrSoloSist_SE[j][l]=  (TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloSist");  
	fHistYieldvsErrErrSoloStatMB_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloStatMB");
	fHistYieldvsErrErrSoloSistMB_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloSistMB");

	// fHistYieldvsErrSoloSistMB_SE[j][l]     ->Sumw2();
	// fHistYieldvsErrSoloStatMB_SE[j][l]   ->Sumw2();
	// fHistYieldvsErrSoloStat_SE[j][l]     ->Sumw2();
	// fHistYieldvsErrSoloSist_SE[j][l]  	  ->Sumw2();
	// fHistYieldvsErrSoloStatRatio_SE[j][l]->Sumw2();
	// fHistYieldvsErrSoloSistRatio_SE[j][l]->Sumw2();
	// fHistYieldvsErrErrSoloStat_SE[j][l]  ->Sumw2();
	// fHistYieldvsErrErrSoloSist_SE[j][l]  ->Sumw2();
	// fHistYieldvsErrErrSoloStatMB_SE[j][l]->Sumw2();
	//fHistYieldvsErrErrSoloSistMB_SE[j][l]->Sumw2();
      
	fHistYieldvsErrSoloSistMB[j][l]   ->Divide(fHistYieldvsErrSoloSistMB_SE[j][l]    );
	fHistYieldvsErrSoloStatMB[j][l]   ->Divide(fHistYieldvsErrSoloStatMB_SE[j][l]    );
	fHistYieldvsErrSoloStat[j][l]     ->Divide(fHistYieldvsErrSoloStat_SE[j][l]      );
	fHistYieldvsErrSoloSist[j][l]     ->Divide(fHistYieldvsErrSoloSist_SE[j][l]      );
	fHistYieldvsErrSoloStatRatio[j][l]->Divide(fHistYieldvsErrSoloStatRatio_SE[j][l] );
	fHistYieldvsErrSoloSistRatio[j][l]->Divide(fHistYieldvsErrSoloSistRatio_SE[j][l] );

	// fHistYieldvsErrErrSoloStat[j][l]  ->Divide(fHistYieldvsErrErrSoloStat_SE[j][l]   );
	// fHistYieldvsErrErrSoloSist[j][l]  ->Divide(fHistYieldvsErrErrSoloSist_SE[j][l]   );
	// fHistYieldvsErrErrSoloStatMB[j][l]->Divide(fHistYieldvsErrErrSoloStatMB_SE[j][l] );
	// fHistYieldvsErrErrSoloSistMB[j][l]->Divide(fHistYieldvsErrErrSoloSistMB_SE[j][l] );
        if (EtaPhiEnlargedStudy) ScaleFactorJet=1;
	if (l==0){
	  fHistYieldvsErrSoloSistMB[j][l]   ->Scale(ScaleFactorJet);
	  fHistYieldvsErrSoloStatMB[j][l]   ->Scale(ScaleFactorJet);
	  fHistYieldvsErrSoloStat[j][l]     ->Scale(ScaleFactorJet);
	  fHistYieldvsErrSoloSist[j][l]     ->Scale(ScaleFactorJet);
	  fHistYieldvsErrSoloStatRatio[j][l]->Scale(ScaleFactorJet);
	  fHistYieldvsErrSoloSistRatio[j][l]->Scale(ScaleFactorJet);
	  fHistYieldvsErrErrSoloStat[j][l]  ->Scale(ScaleFactorJet);
	  fHistYieldvsErrErrSoloSist[j][l]  ->Scale(ScaleFactorJet);
	  fHistYieldvsErrErrSoloStatMB[j][l]->Scale(ScaleFactorJet);
	  fHistYieldvsErrErrSoloSistMB[j][l]->Scale(ScaleFactorJet);
	}

	for (Int_t m =0; m < nummolt; m++){
	  cout << "region under study " << l << "Data or MC " << DATAorMC[j]  << "bin content of mult : " << mult[m] << "  " <<  	  fHistYieldvsErrSoloStat[j][l]->GetBinContent(fHistYieldvsErrSoloStat[j][l]->GetXaxis()->FindBin(mult[m]))<< endl;
	}

      }

    }//Data or MC
  } //end loop number of regions




  
  cout << "...Done\n\n" << endl;

  //********cfr data MC separatamente per J, OJ, J+OJ*********************************************************
  cout << "I'm comparing data to MC " << endl;
  Int_t col=0;
  TCanvas *canvas[6];
  for(Int_t l=0; l<NumberOfRegions; l++){ //loop sulla region (jet, OJ, inclusive)
    Int_t count=0;
    canvas[l] = new TCanvas(Form("canvas%i", l), Region[l], 1300, 1000);
    legend[l]=new TLegend(0.7,0.1, 0.9, 0.3);
    legendtype[l]=new TLegend(0.7,0.35, 0.9, 0.45);
    legendbis[l]=new TLegend(0.1,0.7, 0.3, 0.9);
    cout << "region " << l << endl;
    for(Int_t j=DataOrMCMin; j<DataOrMCMax; j++){ //loop data or MC
      cout << "DataorMC" << DataOrMCMin << endl;
      if (!ishhCorr)      fHistYieldvsErrSoloStat[j][l]->SetTitle ("Yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloStat[j][l]->SetTitle ("Yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloStat[j][l]->SetTitle ("Ratio of yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (ishhCorr)  fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetTitle ("N_{charged}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (isSE)  fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{charged}");      

      fHistYieldvsErrSoloStat[j][l]->SetMarkerColor(Color[l]);
      fHistYieldvsErrSoloStat[j][l]->SetMarkerStyle(MarkerOrigin[l]);
      if (j==1)      fHistYieldvsErrSoloStat[j][l]->SetMarkerStyle(MarkerMC[l]);
      fHistYieldvsErrSoloStat[j][l]->SetLineColor(Color[l]);
      if (ishhCorr)      fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0.0001,7);
      else if (!ishhCorr && !isSEBis){
	fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0.0001,0.35);
	if (l==3) 	fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0.0001,5);
      }
      else fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0.0001,0.14);
      fHistYieldvsErrSoloStat[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloStat[j][l]->SetMarkerSize(msize);
      if (l==2)    fHistYieldvsErrSoloStat[j][l]->SetMarkerSize(3);
      fHistYieldvsErrSoloStat[j][l]->Draw("samee");
    
      if (!ishhCorr)      fHistYieldvsErrSoloSist[j][l]->SetTitle ("Yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloSist[j][l]->SetTitle ("Yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloSist[j][l]->SetTitle ("Ratio of yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (ishhCorr)  fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetTitle ("N_{charged}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (isSE)  fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{charged}");      
      fHistYieldvsErrSoloSist[j][l]->SetMarkerColor(Color[l]);
      fHistYieldvsErrSoloSist[j][l]->SetMarkerStyle(MarkerOrigin[l]);
      if (j==1)       fHistYieldvsErrSoloSist[j][l]->SetMarkerStyle(MarkerMC[l]);
      fHistYieldvsErrSoloSist[j][l]->SetLineColor(Color[l]);
      fHistYieldvsErrSoloSist[j][l]->SetFillStyle(0);
      if (ishhCorr)      fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0.0001,7);
      else if (!ishhCorr && !isSEBis){
	fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0.0001,0.35);
	if (l==3) 	fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0.0001,5);
      }
      else fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0.0001,0.14);

      fHistYieldvsErrSoloSist[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloSist[j][l]->SetMarkerSize(msize);
      if (l==2)    fHistYieldvsErrSoloSist[j][l]->SetMarkerSize(3);
      fHistYieldvsErrSoloSist[j][l]->Draw("samee2");
      
      legend[l]->AddEntry(fHistYieldvsErrSoloStat[j][l], JetOrNot[l]+DATAorMC[j], "pl");
      if(j==1)    legendtype[l]->AddEntry(fHistYieldvsErrSoloStatMB[j][l], "stat. ", "el");
      if(j==1)    legendtype[l]->AddEntry(fHistYieldvsErrSoloSistMB[j][l], "syst. ", "f");
      if(j==1)    legend[l]->Draw();
      //      if(j==1)    legendtype[l]->Draw();


    }


    // for(Int_t l=0; l<NumberOfRegions; l++){
    fileout->WriteTObject(canvas[l]);
    canvas[l]->SaveAs(nomefileoutputPDF+".pdf(");
    //   fileout->WriteTObject(canvaserr[l]);
    // }

    if (isDataOrMC==2 && NumberOfRegions<=3)   canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/CfrDataMCRatioOutput/Yield"+Region[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
    if (isDataOrMC==2 && NumberOfRegions>3 && l>2)   canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/CfrDataMCRatioOutput/Yield"+Region[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));

    //  if (MasterThesis)    canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/Yield"+Region[l]+".pdf");
  }
  cout << "...Done \n\n "<< endl;

  //cfr J, OJ, J+OJ separatamente per data e MC****************************************************
  TCanvas *canvasratiobis[2];

  for(Int_t l=DataOrMCMin; l<DataOrMCMax; l++){ //loop data or MC
    Int_t count=0;
    canvasratiobis[l] = new TCanvas(Form("canvasratiobis%i", l), "Yieldratio_"+DATAorMC[l], 1300, 1000);

    for(Int_t j=0; j<NumberOfRegions; j++){ //loop J, OJ, J+OJ

      if (!ishhCorr)      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle ("Relative yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle ("Relative yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle ("Ratio of relative yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetTitle ("(N_{K^{0}_{S}}/N_{Trigg}) / (N_{K^{0}_{S}}/N_{Trigg})_{0-100 %}");
      if (ishhCorr)  fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetTitle ("(N_{charged}/N_{Trigg}) / (N_{charged}/N_{Trigg})_{0-100 %}");
      if (isSE)  fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetTitle ("");      
      //      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle("");
      fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetRangeUser(0, 2.2); //..
      fHistYieldvsErrSoloStatRatio[l][j]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerSize(msize);
      if (j==2)    fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerSize(3);
      fHistYieldvsErrSoloStatRatio[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrSoloStatRatio[l][j]->Draw("samee");
 
      //      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle("");
      if (!ishhCorr)      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle ("Relative yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle ("Relative yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle ("Ratio of relative yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetTitle ("(N_{K^{0}_{S}}/N_{Trigg}) / (N_{K^{0}_{S}}/N_{Trigg})_{0-100 %}");
      if (ishhCorr)  fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetTitle ("(N_{charged}/N_{Trigg}) / (N_{charged}/N_{Trigg})_{0-100 %}");
      if (isSE)  fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetTitle ("");      
      fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetRangeUser(0, 2.2);//..
      fHistYieldvsErrSoloSistRatio[l][j]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      fHistYieldvsErrSoloSistRatio[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrSoloSistRatio[l][j]->SetFillStyle(0);
      fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerSize(msize);
      if (j==2)    fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerSize(3);
      fHistYieldvsErrSoloSistRatio[l][j]->Draw("samee2");
        
      legendbis[l]->AddEntry(   fHistYieldvsErrSoloStatRatio[l][j], JetOrNot[j]+DATAorMC[l], "pl");
      if(j==0)legendbis[l]->Draw();
      //      if(j==0)legendtype[l]->Draw();
      // fHistYieldvsErrSoloSist[l][j]->SetMarkerColor(Color[j]);
      // fHistYieldvsErrSoloSist[l][j]->Draw("samee");
    }
    canvasratiobis[l]->Close();
    //   canvasratiobis[l]->SaveAs("FinalOutput/DATA2016/YieldRatio"+DATAorMC[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
    //  if (MasterThesis)        canvasratiobis[l]->SaveAs("FinalOutput/DATA2016/YieldRatio"+DATAorMC[l]+".pdf");
  }


    if (!ishhCorr){
      for(Int_t k=1; k <= fHistYieldDatiPubblicati->GetNbinsX(); k++){
	fHistYieldDatiPubblicati->SetBinError(k,     fHistYieldErroriDatiPubblicati->GetBinContent(k));
      }
      fHistYieldDatiPubblicati->Scale(1.6);
      for (int i=0;i< fGraphYieldDatiPubblicati->GetN();i++)  fGraphYieldDatiPubblicati->GetY()[i] *= 1.6; //I'm scaling the TGraph

      for(Int_t k=1; k <= fHistYieldDatiPubblicatiSE->GetNbinsX(); k++){
	fHistYieldDatiPubblicatiSE->SetBinError(k,     fHistYieldErroriDatiPubblicatiSE->GetBinContent(k));
      }
      for(Int_t k=1; k <= fHistYieldDatiPubblicatiSE13TeV->GetNbinsX(); k++){
	fHistYieldDatiPubblicatiSE13TeV->SetBinError(k,     fHistYieldErroriDatiPubblicatiSE13TeV->GetBinContent(k));
      }

    }
    //    fHistYieldDatiPubblicati->Scale(1.6/2*2.2);
    fHistYieldDatiPubblicatiSE->Sumw2();
    fHistYieldDatiPubblicatiSE->SetLineColor(801);
    fHistYieldDatiPubblicatiSE->Scale(0.5); // i dati pubblicati sono relativi a 2*Nk0s/(Npi+ + Npi-)
    fHistYieldDatiPubblicatiSE13TeV->Sumw2();
    fHistYieldDatiPubblicatiSE13TeV->SetLineColor(881);
    fHistYieldDatiPubblicatiSE13TeV->Scale(0.5); // i dati pubblicati sono relativi a 2*Nk0s/(Npi+ + Npi-)

  TSpline3 *splinePub;
  if (isSEBis)  splinePub = new TSpline3( fHistYieldDatiPubblicatiSE13TeV, "splineK0sPionRatio");
  else splinePub = new TSpline3("splineK0s", fGraphYieldDatiPubblicati);
  splinePub->SetLineColor(909);
  TF1 * PublishedDataFit= new TF1( "PublishedDataFit","pol1", 2, 25);
  cout << "\n fit of K0s yield at 13 TeV " << endl;
  fHistYieldDatiPubblicati->Fit(PublishedDataFit,"R");

  TCanvas *canvasnotratio[2];
  for(Int_t l=DataOrMCMin; l<DataOrMCMax; l++){ //loop data or MC
    Int_t count=0;
    canvasnotratio[l] = new TCanvas(Form("canvasnotratio%i", l), "Yield_"+DATAorMC[l], 1300, 1000);
    cout << l << endl;
    for(Int_t j=0; j<NumberOfRegions; j++){
      if (!ishhCorr && !isSEBis)  {
	fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 0.35);
	if (j==3) fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 5);
      }
      else if (ishhCorr)       fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 7);
      else {
	fHistYieldvsErrSoloStat[l][j]->Scale(ScaleFactorPions);
	fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0,0.1);
      }
      fHistYieldvsErrSoloStat[l][j]->GetXaxis()->SetRangeUser(0,25);
      // fHistYieldvsErrSoloStat[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrSoloStat[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      // fHistYieldvsErrSoloStat[l][j]->SetMarkerSize(msize);
      // if (j==2)    fHistYieldvsErrSoloStat[l][j]->SetMarkerSize(3);
      // fHistYieldvsErrSoloStat[l][j]->SetLineColor(Color[j]);
      pol0[l][j]= new TF1(Form("pol0_%i_%i", l, j), "pol0", 0,25);
      pol1[l][j]= new TF1(Form("pol1_%i_%i", l, j), "pol1", 0,25);
      pol0[l][j]->SetParName(0, "q");      pol1[l][j]->SetParName(0, "q");
      pol1[l][j]->SetParName(1, "m");
      pol1[l][j]->SetLineColor(Color[j]);
      pol1[l][j]->SetLineWidth(1);
      pol0[l][j]->SetLineColor(kMagenta);
      fHistYieldvsErrSoloStat[l][j]->Draw("samee");
      cout << "pol0 fit for " << DATAorMC[l] << " region " << j << endl;
      if (isSEBis)       fHistYieldvsErrSoloStat[l][j]->Fit(pol0[l][j], "R0");
      cout << "pol1 fit for " << DATAorMC[l] << " region " << j << endl;
      if (isSEBis)       fHistYieldvsErrSoloStat[l][j]->Fit(pol1[l][j], "R0");
              
      if (!ishhCorr && !isSEBis)  {
	fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 0.35);
	if (j==3) fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 5);
      }

      else if (ishhCorr)       fHistYieldvsErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0, 7);
      else {
	fHistYieldvsErrSoloSist[l][j]->Scale(ScaleFactorPions);
	fHistYieldvsErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0,0.1);
      }
      fHistYieldvsErrSoloSist[l][j]->GetXaxis()->SetRangeUser(0,25);
      // fHistYieldvsErrSoloSist[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrSoloSist[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      // fHistYieldvsErrSoloSist[l][j]->SetLineColor(Color[j]);
      // fHistYieldvsErrSoloSist[l][j]->SetFillStyle(0);
      // fHistYieldvsErrSoloSist[l][j]->SetMarkerSize(msize);
      // if (j==2)    fHistYieldvsErrSoloSist[l][j]->SetMarkerSize(3);
      fHistYieldvsErrSoloSist[l][j]->Draw("samee2");

      if (!ishhCorr && !isSEBis && l==0 &&NumberOfRegions!=3 && j==0) {
	fHistYieldDatiPubblicati->Draw("same");
	//	splinePub->Draw("same");
	PublishedDataFit->Draw("same");
	legendbis[l]->AddEntry(  fHistYieldDatiPubblicati , "13 TeV |y| < 0.5 (err. syst.total)", "pl");
      }
      if (!ishhCorr && isSEBis && l==0 && j==0){
	fHistYieldDatiPubblicatiSE->Draw("same");
	fHistYieldDatiPubblicatiSE13TeV->Draw("same");
	legendbis[l]->AddEntry(  fHistYieldDatiPubblicatiSE , "7 TeV |y| < 0.5 (err. syst.uncorr)", "pl");
	legendbis[l]->AddEntry(  fHistYieldDatiPubblicatiSE13TeV , "13 TeV |y| < 0.5 (err. syst.uncorr)", "pl");
	//	splinePub->Draw("same");
      }
      if(j==0)legendbis[l]->Draw();
      //      if(j==0)legendtype[l]->Draw();
    }
    canvasnotratio[l]->SaveAs("FinalOutput/DATA2016/CfrDataMCRatioOutput/Yield"+DATAorMC[l]+hhCorr[WhatKindOf]+Form("NumberRegions%i_PtMin%.1f.pdf",NumberOfRegions, PtTrigMin));
    canvasnotratio[l]->SaveAs(nomefileoutputPDF+".pdf");
    //  if (MasterThesis)         canvasnotratio[l]->SaveAs("FinalOutput/DATA2016/Yield"+DATAorMC[l]+".pdf");
  }

  //    TCanvas *canvasCfrDatiPubblicati;
  // canvasCfrDatiPubblicati = new TCanvas("canvasCfrDatiPubblicati", "YieldNotScaled", 1300, 1000);
  //fHistYieldvsErrSoloSist[0][5]->Draw("samee2");
  // fHistYieldvsErrSoloStat[0][5]->Draw("samee");

  cout << "\n\n\ndisegno errori relativi in tre canvas (1. JET, 2. OJ 3.J+OJ)" << endl;
  TCanvas *canvaserr[NumberTypeAnalysis];

  
  for(Int_t l=0; l<NumberOfRegions; l++){
    Int_t count=0;
    canvaserr[l] = new TCanvas(Form("canvaserr%i", l), "ErrorOnRegion" +Region[l], 1300, 1000);
    cout << l << endl;
   
    for(Int_t j=DataOrMCMin; j<DataOrMCMax; j++){ //loop data or MC   
      //    cout << l << "  " << j << endl;
      count++;    
      //      fHistYieldvsErrErrSoloStat[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloStat[j][l]->SetMarkerColor(2);
      fHistYieldvsErrErrSoloStat[j][l]->SetMarkerStyle(Marker[j]);
      fHistYieldvsErrErrSoloStat[j][l]->SetLineColor(2);
      fHistYieldvsErrErrSoloStat[j][l]->Draw("sameep");

      fHistYieldvsErrErrSoloStat[j][l]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[j][l]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[j][l]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[j][l]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrErrSoloStat[j][l]->SetMarkerSize(msize);
      
      //      fHistYieldvsErrErrSoloStatMB[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloStatMB[j][l]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloStatMB[j][l]->SetMarkerStyle(MarkerMBStat[j]);
      fHistYieldvsErrErrSoloStatMB[j][l]->SetLineColor(1);
      fHistYieldvsErrErrSoloStatMB[j][l]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloStatMB[j][l]->GetXaxis()->SetRangeUser(0,25);
      //      fHistYieldvsErrErrSoloStatMB[j][l]->Draw("sameep");

      //      fHistYieldvsErrErrSoloSist[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloSist[j][l]->SetMarkerColor(4);
      fHistYieldvsErrErrSoloSist[j][l]->SetMarkerStyle(Marker[j]);
      fHistYieldvsErrErrSoloSist[j][l]->SetLineColor(4);
      fHistYieldvsErrErrSoloSist[j][l]->SetFillStyle(0);
      fHistYieldvsErrErrSoloSist[j][l]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSist[j][l]->Draw("sameep");

      fHistYieldvsErrErrSoloSist[j][l]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[j][l]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[j][l]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[j][l]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0,0.14);
      fHistYieldvsErrErrSoloSist[j][l]->GetXaxis()->SetRangeUser(0,25);

      //      fHistYieldvsErrErrSoloSistMB[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloSistMB[j][l]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloSistMB[j][l]->SetMarkerStyle(MarkerMBSist[j]);
      fHistYieldvsErrErrSoloSistMB[j][l]->SetLineColor(1);
      fHistYieldvsErrErrSoloSistMB[j][l]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSistMB[j][l]->GetXaxis()->SetRangeUser(0,25);
      //      fHistYieldvsErrErrSoloSistMB[j][l]->Draw("sameep");

        
      if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloStat[j][l], StatErr[j], "pl");
      if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloSist[j][l], SysErr[j], "pl");
      // if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloStatMB[j][l], StatErrMB[j], "pl");
      // if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloSistMB[j][l], SysErrMB[j], "pl");
      if(j ==DataOrMCMin)legenderr->Draw();
    
    }
    //    canvaserr[l]->Close();
    //    canvaserr[l] ->SaveAs("FinalOutput/DATA2016/RelErr"+Region[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
    if (MasterThesis)     canvaserr[l] ->SaveAs("FinalOutput/DATA2016/RelErr"+Region[l]+".pdf");
  }


  if (isSEBis){
    TH1F*    fHistRatioBulkJet[2];
    TH1F*    fHistRatioInclusiveJet[2];
    TCanvas *canvasDoubleRatio;
    //I perform ratio (Bulk ratio)/(Jet ratio)
    for(Int_t j=DataOrMCMin; j<DataOrMCMax; j++){ //loop data or MC
      //canvasDoubleRatio[j]= new TCanvas ("SEDoubleRatio"+DATAorMC[j], "SEDoubleRatio"+DATAorMC[j], 1300, 1000);
      canvasDoubleRatio= new TCanvas ("SEDoubleRatio", "SEDoubleRatio", 1300, 1000);
      fHistRatioBulkJet[j] = (TH1F*)	fHistYieldvsErrSoloStat[j][1] -> Clone("fHistRatioBulkJet"+DATAorMC[j]);
      fHistRatioBulkJet[j]->Sumw2();
      fHistRatioBulkJet[j]->Divide(	fHistYieldvsErrSoloStat[j][0]);

      fHistRatioInclusiveJet[j] = (TH1F*)	fHistYieldvsErrSoloStat[j][2] -> Clone("fHistRatioInclusiveJet"+DATAorMC[j]);
      fHistRatioInclusiveJet[j]->Sumw2();
      fHistRatioInclusiveJet[j]->Divide(	fHistYieldvsErrSoloStat[j][0]);

      fHistRatioBulkJet[j]->Draw("same");
      fHistRatioInclusiveJet[j]->Draw("same");
      fHistRatioInclusiveJet[j]->GetYaxis()->SetRangeUser(1.3,1.7);
      fHistRatioBulkJet[j]->GetYaxis()->SetRangeUser(0.9,1.1);

    }
  }
 
  cout << "\n\n\ndisegno errori relativi in due canvas (1. DATA, 2. MC)" << endl;
  TCanvas *canvaserrall[2];
  for(Int_t l=DataOrMCMin; l<DataOrMCMax; l++){ //loop data or MC   
    Int_t count=0;
    canvaserrall[l] = new TCanvas(Form("canvaserrall%i", l), Form("canvaserrall%i", l), 1300, 1000);
    cout << l << endl;
   
    for(Int_t j=0; j<NumberOfRegions; j++){
      //    cout << l << "  " << j << endl;
      count++;
      //      fHistYieldvsErrErrSoloStat[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloStat[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrErrSoloStat[l][j]->SetMarkerStyle(MarkerTris[0]);
      fHistYieldvsErrErrSoloStat[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrErrSoloStat[l][j]->Draw("sameep");

      fHistYieldvsErrErrSoloStat[l][j]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[l][j]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[l][j]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[l][j]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[l][j]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrErrSoloStat[l][j]->SetMarkerSize(msize);
            
      //      fHistYieldvsErrErrSoloStatMB[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloStatMB[l][j]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloStatMB[l][j]->SetMarkerStyle(MarkerTris[0]);
      fHistYieldvsErrErrSoloStatMB[l][j]->SetLineColor(1);
      fHistYieldvsErrErrSoloStatMB[l][j]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloStatMB[l][j]->GetXaxis()->SetRangeUser(0,5);
      //fHistYieldvsErrErrSoloStatMB[l][j]->Draw("sameep");

      //      fHistYieldvsErrErrSoloSist[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloSist[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrErrSoloSist[l][j]->SetMarkerStyle(MarkerTris[1]);
      fHistYieldvsErrErrSoloSist[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrErrSoloSist[l][j]->SetFillStyle(0);
      fHistYieldvsErrErrSoloSist[l][j]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSist[l][j]->Draw("sameep");

      fHistYieldvsErrErrSoloSist[l][j]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[l][j]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[l][j]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[l][j]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0,0.14);
      fHistYieldvsErrErrSoloSist[l][j]->GetXaxis()->SetRangeUser(0,25);

      //      fHistYieldvsErrErrSoloSistMB[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloSistMB[l][j]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloSistMB[l][j]->SetMarkerStyle(MarkerTris[1]);
      fHistYieldvsErrErrSoloSistMB[l][j]->SetLineColor(1);
      fHistYieldvsErrErrSoloSistMB[l][j]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSistMB[l][j]->GetXaxis()->SetRangeUser(0,25);
      //fHistYieldvsErrErrSoloSistMB[l][j]->Draw("sameep");
       
      if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloStat[l][j], StatErrBis[j], "pl");
      if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloSist[l][j], SysErrBis[j], "pl");
      // if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloStatMB[l][j], StatErrMB[j], "pl");
      // if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloSistMB[l][j], SysErrMB[j], "pl");
      if(count==3)legenderrbis->Draw();

    }

    //   canvaserrall[l]->Close();
  }


  //Canvas where I compare my results to the published ones (ratio my/published)
  TCanvas *canvasComparison = new TCanvas("canvasComparison","canvasComparison",800,500);
  canvasComparison->Divide(2,1);
  TH1F*  HistRatioToPublishedl[2];
  if (NumberOfRegions==6){
    for(Int_t l=DataOrMCMin; l<DataOrMCMax; l++){ //loop data or MC
      cout << l << endl;
      if (isSEBis)	HistRatioToPublishedl[l]=(TH1F*)	fHistYieldvsErrSoloStat[l][2]->Clone(Form("fHistK0sPionRatioToPublished_%i", l));
      else if (!ishhCorr )HistRatioToPublishedl[l]=(TH1F*)	fHistYieldvsErrSoloStat[l][3]->Clone(Form("fHistK0sYieldRatioToPublished_%i", l));
      else HistRatioToPublishedl[l]=(TH1F*)	fHistYieldvsErrSoloStat[l][3]->Clone(Form("fHistHadronYieldRatioToPublished_%i", l));
      for(Int_t b=1; b <=    HistRatioToPublishedl[l]->GetNbinsX(); b++){
	if ((   HistRatioToPublishedl[l]->GetBinContent(b))!=0){
	  if (isSEBis)	  HistRatioToPublishedl[l]->SetBinContent(b,fHistYieldvsErrSoloStat[l][2]->GetBinContent(b)/splinePub->Eval(   HistRatioToPublishedl[l]->GetBinCenter(b)));
	  else 	if (!ishhCorr)  HistRatioToPublishedl[l]->SetBinContent(b,fHistYieldvsErrSoloStat[l][3]->GetBinContent(b)/PublishedDataFit->Eval(   HistRatioToPublishedl[l]->GetBinCenter(b)));
	  cout <<"spline eval " << splinePub->Eval(   HistRatioToPublishedl[l]->GetBinCenter(b))<<   " fit eval " <<  PublishedDataFit->Eval(   HistRatioToPublishedl[l]->GetBinCenter(b))<< " bin content of hist ratio " << HistRatioToPublishedl[l]->GetBinContent(b)<< endl;
	}

      }
      canvasComparison->cd(1);
      if (isSEBis)  {
	fHistYieldvsErrSoloStat[l][2]->SetMarkerSize(1);
	fHistYieldvsErrSoloStat[l][2]->Draw("");
	fHistYieldDatiPubblicatiSE->Draw("same");
	//	splinePub->Draw("same");
      }
      else if (!ishhCorr){
	fHistYieldvsErrSoloStat[l][3]->GetYaxis()->SetRangeUser(0,7);
	fHistYieldvsErrSoloStat[l][3]->Draw("");
	//	fHistYieldDatiPubblicati->Draw("same");
	fGraphYieldDatiPubblicati->Draw("same");
	//	splinePub->Draw("same");
      }
      canvasComparison->cd(2);
      HistRatioToPublishedl[l]->SetMarkerSize(1);
      HistRatioToPublishedl[l]->GetYaxis()->SetRangeUser(0.9,3);
      if (l==0)      HistRatioToPublishedl[l]->Draw(""); //I only compare data
    }//end loop on DataMC
  }//end if on numberofregions
  //end of canvas where I compare my results to the published ones (ratio my/published)


  //  canvaserrall[0] ->SaveAs("FinalOutput/DATA2016/RelErrAllDATA0" +hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin)); 
  //  if (MasterThesis)  canvaserrall[0] ->SaveAs("FinalOutput/DATA2016/RelErrAllDATA0.pdf"); 
  //canvaserrall[1] ->SaveAs(Form("FinalOutput/DATA2016/RelErrAllDATA1_PtMin%.1f.pdf", PtTrigMin)); 
  canvasComparison->SaveAs(nomefileoutputPDF+".pdf)");
  fileout->WriteTObject(canvasComparison);
  for(Int_t l=DataOrMCMin; l<DataOrMCMax; l++){
    //  fileout->WriteTObject(canvasratiobis[l]);
    fileout->WriteTObject(canvasnotratio[l]);
    //  fileout->WriteTObject(canvaserrall[l]);
  }
    fileout->WriteTObject(fHistYieldvsErrSoloStat[0][0]);
  fileout->Close();
  cout << "\n\n ho creato il file " << nomefileoutput << endl;
  cout << "e il file pdf " << nomefileoutputPDF+".pdf" << endl;
  cout << "name of saved histo " <<     fHistYieldvsErrSoloStat[0][0]->GetName()<< endl;
  if (isSEBis){
    cout << "\n\n***************************************************************"<< endl;
    cout<< "Be careful! K0s/h ratios have been scaled to properly compare with 7TeV pulbished data (both my data and published ones have been scaled)" << endl;
  }
}

