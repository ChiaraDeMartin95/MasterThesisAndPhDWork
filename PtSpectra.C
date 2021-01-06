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
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>


void StyleHisto(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.8);
  histo->SetTitle(title);
}

void PtSpectra( Int_t type=0,  Int_t TypeAnalysis=1,Int_t numsysV0index=400, Int_t numsysTriggerindex=99, Int_t sysPhi =0, Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0,TString year="1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/,  TString Path1 ="_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, Bool_t isEfficiency=1,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=1, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,   Int_t sys=0, Bool_t SysTuvaWay=0){

  if (TypeAnalysis>2) {cout << "sys errors not yet implemented for these regions " << endl; return;}

  if (TypeAnalysis==0 && sysPhi>2) return;
  if (TypeAnalysis==1 && sysPhi>2) return;
  if (TypeAnalysis==2 && sysPhi>0) return;

  if (type==0){
    PtBinning=1;
    if (!isMC) year= "1617_hK0s";
    else  year= "1617MC_hK0s";
  }
  else if (type==8){
    PtBinning=0;
    if (!isMC) year="Run2DataRed_MECorr_hXi";
    else  year=  "AllMC_hXi";
  }

  TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
  TString RegionTypeOld[3] = {"Jet", "Bulk", "All"};

  if (ishhCorr && type!=0){
    type=0;
  }// {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  const   Int_t NumberTypeAnalysis=11;
  cout << "here's the meaning of different values of TypeAnalysis:" << endl;
  cout << 0 << "in-jet production " << endl;
  cout << 1 << "out-of-jet production  (Delta Phi between 1 and 2)" << endl;
  cout << 2 << "inclusive production (from JetBulkEffCorr)" << endl;
  cout << 3 << "inclusive production not scaled by DeltaEta and DeltaPhi (for comparison with published data) (from JetBulkEffCorr)" << endl;

  if (TypeAnalysis>10) return;

  if (year != "2018f1_extra" && year != "2016k" && year != "2018f1_extra_onlyTriggerWithHighestPt" && year != "2016k_onlyTriggerWithHighestPt") {
    //    cout << "output file should be changed: it must include the year name to avoid overwriting output files " << endl;
    //    return;
  } 

  gStyle->SetOptStat(0);

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *filein;

  TFile *fileSignalSys;
  TFile *fileSysTuva;
  TFile *fileDCAzTrigger;
  TFile *fileinDeltaEtaSys;
  TString PathInDeltaEtaSys;
  TString PathInSignalSys;
  TString PathInDCAzSys;

  TString PathIn0;
  TString PathIn1;
  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);
  file+=Path1;
  TString PathInBis =  "FinalOutput/AnalysisResults" + year/* + Path1*/  + ".root";
  if (isMC && isEfficiency)  PathInBis =  "FinalOutput/AnalysisResults" + year  + "_MCEff" /*+ Path1*/ +".root";
  if (ishhCorr && !isMC)  PathInBis =  "FinalOutput/AnalysisResults" + year  +"_hhCorr" +Path1 + ".root";
  if (ishhCorr && isMC)  PathInBis =  "FinalOutput/AnalysisResults" + year  + "_hhCorr_MCEff" + Path1 + ".root";
  cout << "path in (task output) "<< PathInBis << endl;
  fileinbis=new TFile(PathInBis,"");
  if (!fileinbis){cout << PathInBis << " does not exist " << endl; return;}


  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;

  const Int_t numSysV0 = 7;//num eff sistematici associati alla selezione delle particelle associate
  const Int_t numSysV0hh = 3; //num eff sistematici associati alla selezione delle particelle associate (nel caso ishhCorr)
  const Int_t numSyst  = 13; //all systematics together(numsysv0 + numsystphi + systematic related to bin counting region (1))
  const Int_t numSysthh  = 9;
  const Int_t numsystPhi  = 5; //includo tutti i sistematici anche se i primi 3 non vengono usati in analisi hh (si tratta sia dei sistematici legati alla scelta delle regioni di massa invariante ('centrale' e sidebands) (non per hh) che dei sistematici legati alla scelta delle regioni in DeltaEta di jet e out of jet) 
  cout << "ok " << endl;
  Int_t PtV0Min = 0; //0 el
  if (type>0 || ishhCorr )   PtV0Min =1;
  if (type==0 && PtBinning!=0) PtV0Min=1;
  Int_t numSysV0Global = 0; 
  Int_t numSystGlobal  = 0;
  Int_t numsystPhiGlobal  = 0;
  if (!ishhCorr){
    /*    if (type==0){
	  numSysV0Global = numSysV0; 
	  numSystGlobal  = numSyst;
	  numsystPhiGlobal  = numsystPhi;
	  }*/
    //    if (type!=0){    
    numSysV0Global = 1; 
    numSystGlobal  = 2;
    numsystPhiGlobal  = 0;
    //    }
  }
  if (ishhCorr){
    numSysV0Global = numSysV0hh;
    numSystGlobal  = numSyst;
    numsystPhiGlobal  = numsystPhi;
  }

  cout << "ok " << endl;
  TString JetOrBulk[NumberTypeAnalysis]={"Jet", "Bulk", "All", "AllNS", "AwaySide", "AllButJet", "JetBFit", "JetZYAM", "ASBFit", "ASZYAM", "JetNS"};
  TString DataOrMC[2]={"Data", "MC"};
  Int_t jet=TypeAnalysis;

  cout << "ok " << endl;
  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  cout << "ok " << endl;
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  Float_t SpectrumSup[numtipo][NumberTypeAnalysis]={{0.015,0.2, 0.2,0.8,0.2,0.2,0.015, 0.015, 0.03, 0.03, 0.03},{0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001, 0.002, 0.002, 0.002},{0.015,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001}};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};

  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};

  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};

  for (Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Nmolt[m] = Nmolt0[m];
      Smolt[m] = Smolt0[m];
      SmoltLegend[m] = SmoltLegend0[m];
    }
    else     if (MultBinning==1){
      Nmolt[m] = Nmolt1[m];
      Smolt[m] = Smolt1[m];
      SmoltLegend[m] = SmoltLegend1[m];
    }
    else     if (MultBinning==2){
      Nmolt[m] = Nmolt2[m];
      Smolt[m] = Smolt2[m];
      SmoltLegend[m] = SmoltLegend2[m];
    }

  }



  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", "", ""};  
  //TString SPtV0[numPtV0]={"", "0.5-1", "0.5-1",  "1-1.5","1.5-2", "2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0-0.5"};
    SPtV0[1]={"0.5-1"};
  }

  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", "", ""};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }


  TString SPtV01[numPtV0]={"0-0.1","0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0, 0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0-0.1", "0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0,0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};
  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1 || PtBinning==2) numPtV0Max = numPtV0;
  else numPtV0Max = numPtV0-2;

  if (PtBinning==1){
    for(Int_t v=0; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=0; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }
  cout << "ok " << endl;
  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};
  TString SFit[3]={"exponential", "power law", "retta"};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};
  Int_t inutile[numSyst]=     {0,1,    1, 0,  1,  0,  1,  1,  1,  0,  1,  1,  0};
  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  Int_t MarkerBetter[numSyst]={1,22, 32, 1, 29, 1,   3,  34, 33, 1, 20, 21, 22};
  //  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  //  Int_t Color[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 628,868,1, 909, 1,801,418,860,  1, 881, 7, 1};
  Int_t ColorBetter[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[3] ={628, 418, 600};
  TCanvas *canvasPhiSys[4*(nummolt+1)];
  TCanvas *canvasSpectrumSys[2*(nummolt+1)];
  TCanvas *canvasSpectrumSysBis[2*(nummolt+1)];
  TLegend *legend_Bpassed[nummolt+1][numPtV0];
  TLegend *legend_Bpassed_Spectrum[nummolt+1];      
  TLegend *legend_corr[nummolt+1][numPtV0];      
  TLegend *legend_Corr_Spectrum[nummolt+1];      
  TLegend *legend_UnCorr[nummolt+1][numPtV0];      
  TLegend *legendCorrBis[nummolt+1];      
  TLegend *legendErrorSpectrum;
  TLegend *legendBulk;
  auto legend3= new TLegend(0.55, 0.55, 0.9, 0.9);
  legend3->SetHeader("Selezioni applicate");     
  // auto legend4= new TLegend(0.55, 0.55, 0.9, 0.9); //ter
  auto legend4= new TLegend(0.55, 0.55, 0.9, 0.9); //noter
  //  legend4->SetHeader("Selezioni applicate");     
  
  auto legendYield= new TLegend(0.7, 0.1, 0.9, 0.3);
  legendYield->SetHeader("Fit functions");     
  auto legendYieldError= new TLegend(0.7, 0.1, 0.9, 0.3);
  legendYieldError->SetHeader("Type of error");     

  TLegend *legendError = new TLegend(0.6, 0.6, 0.9, 0.9);

  Float_t YSup=0.012;
  Float_t YInf=-0.001;
  if (type==0 && TypeAnalysis==0) YSup=0.008;
  if (type==8){
    YSup=0.0006;
    YInf = -0.0001;
  }

  TH1D* fHistSpectrum_master[nummolt+1];
  TH1D* fHistSpectrum_masterOnlyStat[nummolt+1];
  TH1D* fHistSpectrum_masterOnlyStatRel[nummolt+1];
  TH1D* fHistSpectrum_masterSystCorr[nummolt+1];
  TH1D* fHistSpectrum_masterSystUnCorr[nummolt+1];
  TH1D* fHistSpectrum[nummolt+1][numSyst];  

  TH1D* fHistPhiDistr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosist[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistRelError[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistRelErrSE[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistSE[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistSETuva[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistDCAz[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistPurity[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistDeltaEta[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solostat[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_master[nummolt+1][numPtV0];

  TF1*  gauss[nummolt+1][numPtV0][numSyst];
  TF1*  gaussint[nummolt+1][numPtV0][numSyst];

  Double_t YieldPerErrore[nummolt+1][numSyst]={0};
  Double_t YieldDefault[nummolt+1][numSyst]= {0};
  
  TH1F*   fHistV0EfficiencyPtBins[nummolt];
  TH1F*   HistContV0PtBins[nummolt];

  TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Lambda", "Xi","Xi", "Omega", "Omega", "Xi", "Omega"};
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask"+dirinputtype[type]);
  if (!dir) {cout << " input dir not found " << endl; return;}

  TList *list = (TList*)dir->Get("MyOutputContainer");
  if (!list) {cout << " input list not found " << endl; return;}


  TString stringout;
  stringout = Dir+"/DATA"+year0+"/PtSpectra" +hhCorr[ishhCorr];
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_SysPhi%i_PtMin%.1f_", sysPhi, PtTrigMin);
  stringout+= RegionType[TypeAnalysis];
  if (SysTuvaWay) stringout+= "_SysTuvaWay";
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  cout << " definition integral regions " << endl;
  //********************************************************************* 
  //definition of integral regions in deltaPhi distribution
  //********************************************************************* 
  Float_t ALowBin[3]={-1}; 
  Float_t AUpBin[3]={1};
  Float_t ALowBinFit[3]={-1};
  Float_t AUpBinFit[3]={1};
  Float_t DeltaPhiWidth[3]={0};
  Float_t DeltaPhiWidthApprox[3]={0};

  if (TypeAnalysis==0 || TypeAnalysis==10){
    if (!ishhCorr && (type==0 || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-0.8};
      AUpBinFit[0]=	AUpBin[0]={0.8};
      
      ALowBinFit[1]=	ALowBin[1]={-1.05};
      AUpBinFit[1]=	AUpBin[1]={1.05};

      ALowBinFit[2]=	ALowBin[2]={-1.2};
      AUpBinFit[2]=	AUpBin[2]={1.2};

    }
    else if (!ishhCorr && (type==4 || type==5 || type==8)){
      ALowBinFit[0]=	ALowBin[0]={-0.7};
      AUpBinFit[0]=	AUpBin[0]={0.7};
      
      ALowBinFit[1]=	ALowBin[1]={-1.1};
      AUpBinFit[1]=	AUpBin[1]={1.1};
    }
    if (ishhCorr){
      ALowBinFit[0]=	ALowBin[0]={-1.2};
      AUpBinFit[0]=	AUpBin[0]={1.2};
      
      ALowBinFit[1]=	ALowBin[1]={-1.4};
      AUpBinFit[1]=	AUpBin[1]={1.4};
    }
    
  }
  else if (TypeAnalysis==1){
    ALowBinFit[0]=	ALowBin[0]={1.1};
    AUpBinFit[0]=	AUpBin[0]={2};

    ALowBinFit[1]=	ALowBin[1]={1};
    AUpBinFit[1]=	AUpBin[1]={2};

    ALowBinFit[2]=	ALowBin[2]={1.1};
    AUpBinFit[2]=	AUpBin[2]={2.1};

    /*
    ALowBinFit[1]=	ALowBin[1]={2};
    AUpBinFit[1]=	AUpBin[1]={4.28};

    ALowBinFit[2]=	ALowBin[2]={1};
    AUpBinFit[2]=	AUpBin[2]={4.28};
    */
  }
  else if (TypeAnalysis==2 || TypeAnalysis==3){
    ALowBinFit[0]=	ALowBin[0]=-0.5*TMath::Pi()+0.01;
    AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
  }
  else if (TypeAnalysis==4){
    ALowBinFit[0]=	ALowBin[0]=TMath::Pi()-1;
    AUpBinFit[0]=	AUpBin[0]=TMath::Pi()+1;
  }
  else if (TypeAnalysis==5){
    ALowBinFit[0]=	ALowBin[0]=1;
    AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
  }
  else if (TypeAnalysis==6 || TypeAnalysis==7){
    if (type==0){
      ALowBinFit[0]=	ALowBin[0]=-0.8;
      AUpBinFit[0]=	AUpBin[0]=0.8;
    }
    else if (type==8){
      ALowBinFit[0]=	ALowBin[0]=-0.7;
      AUpBinFit[0]=	AUpBin[0]=0.7;
    }
  }
  else if (TypeAnalysis==8 || TypeAnalysis==9){
    if (type==0){
      ALowBinFit[0]=	ALowBin[0]=1.6;
      AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
    }
    else if (type==8){
      ALowBinFit[0]=	ALowBin[0]=1.6;
      AUpBinFit[0]=	AUpBin[0]=1.5*TMath::Pi()-0.01;
    }
  }

  //low and upper ranges for the fit****************************
  //************************************************************

  Double_t LowRange[nummolt+1]= {2, 2, 1, 1.5, 1.5, 1.5}; 
  Double_t UpRange[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeSpectrumPart[nummolt+1]= {0};
  Double_t UpRangeSpectrumPart[nummolt+1]= {0};

  Double_t LowRangeJet[nummolt+1]= {1, 1, 1, 1, 1, 1}; 
  Double_t UpRangeJet[nummolt+1]= {8,8,8,8,8,8}; 

  //  Double_t LowRangeBulk[nummolt+1]= {1, 1, 1, 1, 1, 1}; 
  Double_t LowRangeBulk[nummolt+1]= {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; 
  Double_t UpRangeBulk[nummolt+1]= {4,4,4,4,4,4};

  //  Double_t LowRangeAll[nummolt+1]=  {1, 1, 1, 1, 1, 1}; 
  Double_t LowRangeAll[nummolt+1]=  {0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; 
  Double_t UpRangeAll[nummolt+1]= {4,4,4,4,4,4};

  Double_t LowRangeAS[nummolt+1]= {1, 1, 1, 1, 1, 1}; 
  Double_t UpRangeAS[nummolt+1]=  {4,4,4,4,4,4};

  Double_t LowRangeAllButJet[nummolt+1]=  {1, 1, 1, 1, 1, 1}; 
  Double_t UpRangeAllButJet[nummolt+1]={4,4,4,4,4,4};

  Double_t LowRangeJetBFit[nummolt+1]= {1.5, 1.5, 1,1, 1.5,1}; 
  Double_t UpRangeJetBFit[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeJetZYAM[nummolt+1]= {1.5 , 1.5, 0.5,1, 1.5,1}; 
  Double_t UpRangeJetZYAM[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeASBFit[nummolt+1]= {1.5 , 1.5, 0.5,1, 1.5,1}; 
  Double_t UpRangeASBFit[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeASZYAM[nummolt+1]= {1.5 , 1.5, 0.5,1, 1.5,1}; 
  Double_t UpRangeASZYAM[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeJet1[nummolt+1]= {2, 2, 1.5, 1, 1, 1}; 
  Double_t UpRangeJet1[nummolt+1]= {8,8,8,8,8,8}; 
  Double_t LowRangeJet1MC[nummolt+1]= {2, 1.5, 1.5, 1.5, 1.5, 1.5}; 
  Double_t UpRangeJet1MC[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeJet2[nummolt+1]= {2, 1.5, 1, 1, 1, 1}; 
  Double_t UpRangeJet2[nummolt+1]= {8,8,8,8,8,8}; 
  Double_t LowRangeJet2MC[nummolt+1]= {2, 1.5, 1.5, 1.5, 1.5, 1.5}; 
  Double_t UpRangeJet2MC[nummolt+1]= {8,8,8,8,8,8}; 

  if (!ishhCorr && type==0){
    for (Int_t m =0; m<nummolt+1; m++){
      if (m==4 && !isMC)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1.0;
      else   LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.5;
      LowRangeBulk[m]= 0; //0
      LowRangeAll[m]= 0; //0 
      UpRangeAll[m]= 4; 
    }
    if (PtBinning==1 && year=="1617_hK0s") {
      for (Int_t m =0; m<nummolt+1; m++){
	if (PtTrigMin==3){
	  if (!isMC){
	    //	    if (m==5|| m==2 || m==4 )     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    //	    else  {LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;}
	    LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;	  
	    if (m==4) UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;
	    LowRangeBulk[m]= 0.1; 
	    LowRangeAll[m]= 0.1; 
	    UpRangeAll[m]= 2; 
	    UpRangeBulk[m]= 2; 
	    LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 0.1; 
	    UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 2; 
	  }
	} //end 3 
	else if (PtTrigMin==4){
	  if (!isMC){
	    if (m==5)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;
	    else  {LowRangeJet[m] = 0.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 0.1;}
	    UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;	  
	    LowRangeBulk[m]= 0.1; 
	    LowRangeAll[m]= 0.1; 
	    UpRangeAll[m]= 2; 
	    UpRangeBulk[m]= 2; 
	    LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 0.1; 
	    UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 2; 
	  }
	}
      }
    }
  }

  if (MultBinning==1 && type==8){
    for (Int_t m=0; m<nummolt+1; m++){
      if (!isMC)      LowRangeJet[m] = LowRangeJet1[m];
      else if (isMC)      LowRangeJet[m] = LowRangeJet1MC[m];
    }
  }

  if (MultBinning==2 && type==8){
    for (Int_t m=0; m<nummolt+1; m++){
      if (!isMC)      LowRangeJet[m] = LowRangeJet2[m];
      else if (isMC)      LowRangeJet[m] = LowRangeJet2MC[m];
    }
  }

  if (MultBinning==0 && type==8){
    for (Int_t m=0; m<nummolt+1; m++){
      if (isMC)  {
	if (m==3)     LowRangeJet[m] =     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1.5;
	else  if (m==0 || m==4 || m==5)  {LowRangeJet[m] = 1.5;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1;}
	else  if (m==1) {LowRangeJet[m] = 2.0;     LowRangeJetBFit[m] =     LowRangeJetZYAM[m]= 1;}
	else  if (m==2) {LowRangeJet[m] = LowRangeJetBFit[m]=1.5;       LowRangeJetZYAM[m]= 1;}
	UpRangeJet[m] = UpRangeJetBFit[m] = UpRangeJetZYAM[m]=4;

	LowRangeBulk[m]= 1; //not enought statistics below 
	LowRangeAll[m]= 1; 
	UpRangeAll[m]= 4; 
	UpRangeBulk[m]= 4; 
	if (m==5) {
	  LowRangeBulk[m]= 0.5;
	  LowRangeAll[m]= 0.5; 
	  UpRangeAll[m]= 8; 
	  UpRangeBulk[m]= 8; 
	}
      }
      else {
	LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 1.5; 
	if (m==5 || m==3) 	LowRangeASBFit[m]= 	  LowRangeASZYAM[m]= 1; 
	UpRangeASBFit[m]= 	  UpRangeASZYAM[m]= 4; 
      }
    }
  }


  //end of low and upper ranges for the fit****************************
  //************************************************************

  Int_t PtBinMin[nummolt+1]={0};
  for (Int_t m=0; m<nummolt+1; m++){
    if (TypeAnalysis==0 || TypeAnalysis==10){
      LowRange[m] =  LowRangeJet[m];
      UpRange[m] =  UpRangeJet[m];
    }
    if (TypeAnalysis==1){
      LowRange[m] =  LowRangeBulk[m];
      UpRange[m] =  UpRangeBulk[m];
    }
    if (TypeAnalysis==2 || TypeAnalysis==3){
      LowRange[m] =  LowRangeAll[m];
      UpRange[m] =  UpRangeAll[m];
    }
    if (TypeAnalysis==4){
      LowRange[m] =  LowRangeAS[m];
      UpRange[m] =  UpRangeAS[m];
    }
    if (TypeAnalysis==5){
      LowRange[m] =  LowRangeAllButJet[m];
      UpRange[m] =  UpRangeAllButJet[m];
    }
    if (TypeAnalysis==6){
      LowRange[m] =  LowRangeJetBFit[m];
      UpRange[m] =  UpRangeJetBFit[m];
    }
    if (TypeAnalysis==7){
      LowRange[m] =  LowRangeJetZYAM[m];
      UpRange[m] =  UpRangeJetZYAM[m];
    }
    if (TypeAnalysis==8){
      LowRange[m] =  LowRangeASBFit[m];
      UpRange[m] =  UpRangeASBFit[m];
    }
    if (TypeAnalysis==9){
      LowRange[m] =  LowRangeASZYAM[m];
      UpRange[m] =  UpRangeASZYAM[m];
    }

    //bin counting is done at pt> LowRangeSpectrumPart, integral of fit function at pt <  LowRangeSpectrumPart
    //low range spectrum part might be different from the lower limit of fit range
    if  (TypeAnalysis==1 || TypeAnalysis==2){
      if (type>0 && !isMC)    LowRangeSpectrumPart[m] = 0.5;
      else if (type>0 && isMC && m!=5)    LowRangeSpectrumPart[m] = 1.;
      else     LowRangeSpectrumPart[m] = 0.1;
    }
    else  LowRangeSpectrumPart[m] = LowRange[m];
    if (type==8 && TypeAnalysis==1 && m==4) UpRangeSpectrumPart[m] = 4.;
    else if (type==8 && TypeAnalysis==0 && m==4) UpRangeSpectrumPart[m] = 4.;
    else UpRangeSpectrumPart[m] = 8.;
    //    else  LowRangeSpectrumPart[m] = 8;
    cout << " ok " << endl;
    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
    }
  }

  //syst1_limit indica gli effetti sistematici sullo spettro in pT
  Int_t syst1_limit=0;
  if(!TypeAnalysis==1) syst1_limit= 1; //2 for BC systematic
  if (TypeAnalysis==1) syst1_limit= 1;
  if (TypeAnalysis==2) syst1_limit= 1; //non considero sistematico associato a bin counting per total
  if (TypeAnalysis>2) syst1_limit=1;


  Int_t sysV0=0;    
  if(isMC && isEfficiency) file = year + "_MCEff" +Path1;
  //  file+=Form("_PtBinning%i", PtBinning);
  //  file+= Path1;

  //here I start the newly  created part!!

  TString titleX=  "p_{T} (GeV/c)";
  TString titleY=  "1/#Delta #eta #Delta #phi 1/N_{trigg} dN/dp_{T}";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 

  Float_t  YieldPt[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrStat[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSist[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistSE[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistSETuva[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistDCAz[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistPurity[nummolt+1][numPtV0] ={0};
  Float_t  YieldPtErrSistDeltaEta[nummolt+1][numPtV0] ={0};
  TH1D* fHistSpectrumStat[nummolt+1];
  TH1D* fHistSpectrumSist[nummolt+1];
  TH1D* fHistSpectrumSistSE[nummolt+1];
  TH1D* fHistSpectrumSistSETuva[nummolt+1];
  TH1D* fHistSpectrumSistDCAz[nummolt+1];
  TH1D* fHistSpectrumSistPurity[nummolt+1];
  TH1D* fHistSpectrumSistDeltaEta[nummolt+1];
  TH1D* fHistSpectrumStatRelError[nummolt+1];
  TH1D* fHistSpectrumSistRelError[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorSE[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorSETuva[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorDeltaEta[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorDCAz[nummolt+1];
  TH1D* fHistSpectrumSistRelErrorPurity[nummolt+1];

  TF1 * lineat1 = new TF1("pol0","pol0", -1./2*TMath::Pi(), 3./2*TMath::Pi());
  lineat1->SetLineColor(kBlack);
  lineat1->FixParameter(0,0);

  TCanvas* canvasPlotProj[nummolt+1];
  TCanvas* canvasPtSpectra = new TCanvas ("canvasPtSpectra", "canvasPtSpectra", 1300, 800);
  canvasPtSpectra->Divide(3,2);
  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  canvasPtSpectraRelErrorAll->Divide(3,2);

  Float_t LimSup=0.2;
  if (type==0){
    if (TypeAnalysis==0) LimSup =0.02;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSup =0.001;
    else   if (TypeAnalysis==1) LimSup =0.01;
    else   if (TypeAnalysis==2) LimSup =0.01;
  }

  Float_t LimSupError=0.01;
  if (type==0) {
    if (TypeAnalysis==0) LimSupError =0.2;
    else   if (TypeAnalysis==1) LimSupError =0.08;
    else   if (TypeAnalysis==2) LimSupError =0.02;
  }
  else  if (type==8) {
    if (TypeAnalysis==0) LimSupError =0.8;
    else   if (TypeAnalysis==1) LimSupError =0.1;
    else   if (TypeAnalysis==2) LimSupError =0.05;
  }

  for(Int_t m=0; m<nummolt+1; m++){
    cout << " m " << m << endl;
    canvasPlotProj[m] = new TCanvas (Form("canvasPlot_m%i", m), Form("canvasPlot_m%i", m), 1300, 800);
    canvasPlotProj[m]-> Divide((float)numPtV0/2+1, 2);
    fHistSpectrumStat[m]=new TH1D ("fHistSpectrum_"+Smolt[m],"fHistSpectrumStat_"+Smolt[m], numPtV0Max, NPtV0) ;
    fHistSpectrumStat[m]->Sumw2();
    fHistSpectrumSist[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSist_"+Smolt[m]);
    fHistSpectrumSistSE[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistSE_"+Smolt[m]);
    fHistSpectrumSistSETuva[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistSETuva_"+Smolt[m]);
    fHistSpectrumSistDCAz[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDCAz_"+Smolt[m]);
    fHistSpectrumSistPurity[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistPurity_"+Smolt[m]);
    if (TypeAnalysis!=2)    fHistSpectrumSistDeltaEta[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistDeltaEta_"+Smolt[m]);

    //PathIn1 definition
    if (isMC && !isEfficiency) PathIn1=Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file+  Form("_Sys%i", sys)+"_Output.root";
    else {
      PathIn1 = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file;
      if(type>=0){
	PathIn1 +="_"+tipo[type];
	PathIn1 +=Srap[israp];
	PathIn1 +=SSkipAssoc[SkipAssoc];
      }
      PathIn0=PathIn1;
      PathIn1+= hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", 0, sysV0, sys, PtTrigMin)+"_Output.root";
    }
    //    cout << "\n\n" << PathIn1 << endl;

    //PathInSignalSys definition
    TString PathInSys=Dir+"/DATA2016/histo/SignalExtractionStudy";
    TString PathInSysTuva=Dir+"/DATA2016/LoosestTightestTopoSel";
    PathInSys+=year ;
    PathInSysTuva+=year ;
    TString PathInSysDeltaEta=Dir+"/DATA2016/SysDeltaEta";
    if (PtBinning>0) {
      PathInSys +=Form("_PtBinning%i",PtBinning);
      PathInSysTuva +=Form("_PtBinning%i",PtBinning);
      PathInSysDeltaEta +=Form("_PtBinning%i",PtBinning);
    }
    if(type>=0){
      if (!ishhCorr) {
	PathInSys +="_"+tipo[type];
	PathInSysTuva +="_"+tipo[type];
	PathInSysDeltaEta +="_"+tipo[type];
      }
      PathInSys +=Srap[israp];
      PathInSys +=SSkipAssoc[SkipAssoc];
      PathInSysTuva +=Srap[israp];
      PathInSysTuva +=SSkipAssoc[SkipAssoc];
      PathInSysDeltaEta +=Srap[israp];
      PathInSysDeltaEta +=SSkipAssoc[SkipAssoc];
    }

    PathInSignalSys= PathInSys+  Form("_SysT%i_SysV0Num%i_Sys%i_PtMin%.1f.root", 0, numsysV0index, sys, PtTrigMin);
    PathInSysTuva= PathInSysTuva+  Form("_PtMin%.1f_", PtTrigMin)+RegionTypeOld[TypeAnalysis]+".root" ;
    PathInDCAzSys= PathInSys +  Form("_SysTNum%i_SysV0%i_Sys%i_PtMin%.1f.root", numsysTriggerindex, 0, sys, PtTrigMin);

    //    filein = new TFile(PathIn1, "");
    cout <<"systematics assoc to signal extraction from file " <<  PathInSignalSys << endl;
    fileSignalSys = new TFile(PathInSignalSys, "");
    if (SysTuvaWay)    fileSysTuva = new TFile(PathInSysTuva, "");
    fileDCAzTrigger = new TFile(PathInDCAzSys, "");
    
    //    if (!filein) {cout << filein << "does not exist " << endl ; return; }
    if (!fileSignalSys) {cout << PathInSignalSys << "does not exist " << endl ; return; }
    if (!fileSysTuva && SysTuvaWay) {cout << PathInSysTuva << "does not exist " << endl ; return; }
    if (!fileDCAzTrigger) {cout << PathInDCAzSys << "does not exist " << endl ; return; }
    
    if(TypeAnalysis!=2){
      PathInDeltaEtaSys = PathInSysDeltaEta + hhCorr[ishhCorr]+Form("_PtMin%.1f_", PtTrigMin)+RegionType[TypeAnalysis]+".root";
      fileinDeltaEtaSys = new TFile(PathInDeltaEtaSys, "");
      if (!fileinDeltaEtaSys) {cout << PathInDeltaEtaSys << "does not exist " << endl ; return; }
      cout << "syst associated to DeltaEta region from file " <<       PathInDeltaEtaSys<<endl;
    }

    //    cout <<     "I start the loop on pT " << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      //      cout << v << endl;
      if (TypeAnalysis==0){
	if (type==8 && NPtV0[v]<2.5-0.001)  fHistPhiDistr_master[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmoothBis");
	else fHistPhiDistr_master[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");

      }
      if (TypeAnalysis==1) fHistPhiDistr_master[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");
      if (TypeAnalysis==2 || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_master[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr");
      if (TypeAnalysis==3) fHistPhiDistr_master[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled");
      if(!fHistPhiDistr_master[m][v]) {cout << "histo master in file signalsys  not present " << endl; return;}
      fHistPhiDistr_master[m][v]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);

      //getting phi distributions with systematic errors from signal extraction
      //T      if (!SysTuvaWay){
      if (TypeAnalysis==0) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefSys");
      if (TypeAnalysis==1) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_DefSys");
      if (TypeAnalysis==2 || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr_DefSys");
      if (TypeAnalysis==3) fHistPhiDistr_solosistSE[m][v]=(TH1D*)fileSignalSys->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled_DefSys");
      if(!fHistPhiDistr_solosistSE[m][v]){ cout << " histo with signal extraction error not presetn " << endl; return;}
      //T      }
      //T else{
      if (SysTuvaWay){ //T
	fHistPhiDistr_solosistSETuva[m][v]=(TH1D*)fileSysTuva->Get("PhiDistr_solosistSETuvaWay_m"+Smolt[m]+"_v"+SPtV0[v]);
	if(!fHistPhiDistr_solosistSETuva[m][v])  {cout << "histo SETuva not present " << endl; return;}
      } //T
	//T      }
	if (TypeAnalysis==0) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefSys");
	if (TypeAnalysis==1) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_DefSys");
	if (TypeAnalysis==2 || TypeAnalysis==4 || TypeAnalysis==5) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr_DefSys");
	if (TypeAnalysis==3) fHistPhiDistr_solosistDCAz[m][v]=(TH1D*)fileDCAzTrigger->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled_DefSys");
	if(!fHistPhiDistr_solosistDCAz[m][v]){ cout << " histo with DCAzTrigger error not presetn " << endl; return;}
      
      if (TypeAnalysis!=2){
	fHistPhiDistr_solosistDeltaEta[m][v]=(TH1D*)fileinDeltaEtaSys->Get("PhiDistr_solosistDeltaEta_m"+Smolt[m]+"_v"+SPtV0[v]);
	if(!fHistPhiDistr_solosistDeltaEta[m][v])  {cout << "histo DeltaEta not present " << endl; return;}
	fHistPhiDistr_solosistDeltaEta[m][v]->SetName("PhiDistrSysDEta_m"+Smolt[m]+"_v"+SPtV0[v]);
      }

      cout << " I might end if histograms not present...." << endl;


      //      cout << " Histogram is present!" << endl;

      fHistPhiDistr[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solostat[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solostat_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosist[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosist_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistPurity[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistPurity_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistRelError[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistRelError_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistRelErrSE[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistRelErrSE_m"+Smolt[m]+"_v"+SPtV0[v]);

      //setting 2% error on sist purity
      for (Int_t b=1; b<=fHistPhiDistr_master[m][v]->GetNbinsX(); b++){
	//	if (TypeAnalysis!=0) 	fHistPhiDistr_solosistPurity[m][v]->SetBinError(b, fHistPhiDistr_solosistPurity[m][v]->GetBinContent(b)*0.02);
	fHistPhiDistr_solosistPurity[m][v]->SetBinError(b, fHistPhiDistr_solosistPurity[m][v]->GetBinContent(b)*0.02);
	//	else 	fHistPhiDistr_solosistPurity[m][v]->SetBinError(b,0);
      }
      //

      for(Int_t i =0; i <3; i++){
	DeltaPhiWidth[i]=TMath::Abs(fHistPhiDistr_master[m][v]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_master[m][v]->FindBin(ALowBinFit[sysPhi])) - fHistPhiDistr_master[m][v]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_master[m][v]->FindBin(AUpBinFit[sysPhi])));
	DeltaPhiWidthApprox[i]=TMath::Abs(AUpBinFit[i]-ALowBinFit[i]);
	if (TypeAnalysis==0)            DeltaPhiWidth[i]=TMath::Abs(fHistPhiDistr_master[m][v]->GetXaxis()->GetBinLowEdge(fHistPhiDistr_master[m][v]->FindBin(ALowBinFit[0])) - fHistPhiDistr_master[m][v]->GetXaxis()->GetBinUpEdge(fHistPhiDistr_master[m][v]->FindBin(AUpBinFit[0])));
	if (TypeAnalysis==8 || TypeAnalysis==9 || TypeAnalysis==10) DeltaPhiWidth[i]=1; //not scaled for AS and JetNS                                                                            
	cout << " ciao " << endl;
      }
      cout << " v " << v << endl;
    } //end loop v
    cout << " end loop " << endl;
  } //end loop m

  cout << " new loop " << endl;
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      //      cout << "I've set errors " << endl;
      for(Int_t dphi=1; dphi<=fHistPhiDistr_solostat[m][v]->GetNbinsX(); dphi++){
	//cout << "dphi " << dphi << endl;
	//	cout << fHistPhiDistr_solosistDeltaEta[m][v]->GetBinContent(dphi)<< endl;
	/* old
	if (TypeAnalysis==2)	fHistPhiDistr_solosist[m][v]->SetBinError(dphi,       fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi));
	else 	if (TypeAnalysis==1 || TypeAnalysis==0) fHistPhiDistr_solosist[m][v]->SetBinError(dphi,    sqrt(pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)+ pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2)));
	//	cout << "I've set errors " << endl;
	*/ 
	if (TypeAnalysis==2) fHistPhiDistr_solosist[m][v]->SetBinError(dphi,    sqrt(pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)+ pow(fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi),2) + pow(fHistPhiDistr_solosistPurity[m][v]->GetBinError(dphi),2)));
	else 	if (TypeAnalysis==1 || TypeAnalysis==0) fHistPhiDistr_solosist[m][v]->SetBinError(dphi,    sqrt(pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2)+ pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2) + pow(fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi),2) +  pow(fHistPhiDistr_solosistPurity[m][v]->GetBinError(dphi),2)));

	  fHistPhiDistr_solosistRelError[m][v]->SetBinContent( dphi,   TMath::Abs(  fHistPhiDistr_solosist[m][v]->GetBinError(dphi)/fHistPhiDistr_solosist[m][v]->GetBinContent(dphi)));
	  fHistPhiDistr_solosistRelErrSE[m][v]->SetBinContent( dphi,   TMath::Abs(  fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi)/fHistPhiDistr_solosist[m][v]->GetBinContent(dphi)));
      }
      fHistPhiDistr_solosistRelError[m][v]->Smooth();
      fHistPhiDistr_solosistRelErrSE[m][v]->Smooth();
      for(Int_t dphi=1; dphi<=fHistPhiDistr_solostat[m][v]->GetNbinsX(); dphi++){
	  if (type==8 && TypeAnalysis==0 && m==4 && v==PtV0Min+1){
	    cout << "\nnon smoothed " << 	  fHistPhiDistr_solosist[m][v]->GetBinError(dphi);
	    cout << " rel error " << 	  fHistPhiDistr_solosistRelError[m][v]->GetBinContent(dphi) << " content " << fHistPhiDistr_solosist[m][v]->GetBinContent(dphi);
	    fHistPhiDistr_solosist[m][v]->SetBinError(dphi,      fHistPhiDistr_solosistRelError[m][v]->GetBinContent(dphi)*      fHistPhiDistr_solosist[m][v]->GetBinContent(dphi));
	    fHistPhiDistr_solosistSE[m][v]->SetBinError(dphi,      fHistPhiDistr_solosistRelErrSE[m][v]->GetBinContent(dphi)*      fHistPhiDistr_solosist[m][v]->GetBinContent(dphi));
	    cout << " smoothed " << 	  fHistPhiDistr_solosist[m][v]->GetBinError(dphi)<< endl;
	  }

	  //	  cout << "DCAz " << fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi) << " " <<  fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi) << " DeltaEta " <<  fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi) << endl;
	  //	  cout << " total " << fHistPhiDistr_solosist[m][v]->GetBinError(dphi) << endl;
	
      }
      
      YieldPt[m][v]=0;
      YieldPtErrStat[m][v]=0;
      YieldPtErrSist[m][v]=0;
      YieldPtErrSistSE[m][v]=0;
      YieldPtErrSistSETuva[m][v]=0;
      YieldPtErrSistDCAz[m][v]=0;
      YieldPtErrSistPurity[m][v]=0;
      YieldPtErrSistDeltaEta[m][v]=0;
      for(Int_t dphi=fHistPhiDistr_solostat[m][v]->FindBin(ALowBinFit[sysPhi]); dphi<=fHistPhiDistr_solostat[m][v]->FindBin(AUpBinFit[sysPhi]); dphi++){
	YieldPt[m][v]+= fHistPhiDistr_solostat[m][v]->GetBinContent(dphi);
	YieldPtErrStat[m][v]+= pow(fHistPhiDistr_solostat[m][v]->GetBinError(dphi),2);
	YieldPtErrSist[m][v]+= pow(fHistPhiDistr_solosist[m][v]->GetBinError(dphi),2);
	YieldPtErrSistSE[m][v]+= pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2);
	if (SysTuvaWay)	YieldPtErrSistSETuva[m][v]+= pow(fHistPhiDistr_solosistSETuva[m][v]->GetBinError(dphi),2);
	else	YieldPtErrSistSETuva[m][v]+= pow(fHistPhiDistr_solosistSE[m][v]->GetBinError(dphi),2);
	YieldPtErrSistDCAz[m][v]+= pow(fHistPhiDistr_solosistDCAz[m][v]->GetBinError(dphi),2);
	YieldPtErrSistPurity[m][v]+= pow(fHistPhiDistr_solosistPurity[m][v]->GetBinError(dphi),2);
	if (TypeAnalysis!=2)	YieldPtErrSistDeltaEta[m][v]+= pow(fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi),2);
      }
      YieldPtErrStat[m][v]= sqrt(YieldPtErrStat[m][v]);
      YieldPtErrSist[m][v]= sqrt(YieldPtErrSist[m][v]);
      YieldPtErrSistSE[m][v]= sqrt(YieldPtErrSistSE[m][v]);
      YieldPtErrSistSETuva[m][v]= sqrt(YieldPtErrSistSETuva[m][v]);
      YieldPtErrSistDCAz[m][v]= sqrt(YieldPtErrSistDCAz[m][v]);
      YieldPtErrSistPurity[m][v]= sqrt(YieldPtErrSistPurity[m][v]);
      YieldPtErrSistDeltaEta[m][v]= sqrt(YieldPtErrSistDeltaEta[m][v]);

      cout << " m " << m << " v " << NPtV0[v] << "Purity " <<     YieldPtErrSistPurity[m][v]/YieldPt[m][v] << " DCAz " <<       YieldPtErrSistDCAz[m][v]/YieldPt[m][v] << " DeltaEta " <<       YieldPtErrSistDeltaEta[m][v]/YieldPt[m][v] << " SE " <<       YieldPtErrSistSE[m][v]/YieldPt[m][v] << " tot: " <<       YieldPtErrSist[m][v]/YieldPt[m][v]<<endl;
      //      for(Int_t v=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]-0.001); v<fHistSpectrumStat[m]->GetNbinsX() ; v++ ){
      for(Int_t v=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]-0.001); v<fHistSpectrumStat[m]->FindBin(UpRangeSpectrumPart[m]-0.001); v++){
      fHistSpectrumStat[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumStat[m]->SetBinError(v+1, YieldPtErrStat[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumSist[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumSist[m]->SetBinError(v+1, YieldPtErrSist[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

      fHistSpectrumSistSE[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumSistSE[m]->SetBinError(v+1, YieldPtErrSistSE[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

      fHistSpectrumSistSETuva[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumSistSETuva[m]->SetBinError(v+1, YieldPtErrSistSETuva[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

      fHistSpectrumSistDCAz[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumSistDCAz[m]->SetBinError(v+1, YieldPtErrSistDCAz[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

      fHistSpectrumSistPurity[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      fHistSpectrumSistPurity[m]->SetBinError(v+1, YieldPtErrSistPurity[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));

      if (TypeAnalysis!=2){
	fHistSpectrumSistDeltaEta[m]->SetBinContent(v+1, YieldPt[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
	fHistSpectrumSistDeltaEta[m]->SetBinError(v+1, YieldPtErrSistDeltaEta[m][v]/fHistSpectrumStat[m]->GetBinWidth(v+1));
      }

      }
      /*
      for(Int_t b=1;  b<fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); b++){
	fHistSpectrumStat[m]->SetBinContent(b,0);
	fHistSpectrumStat[m]->SetBinError(b,0);
	fHistSpectrumSist[m]->SetBinContent(b,0);
	fHistSpectrumSist[m]->SetBinError(b,0);
	fHistSpectrumSistSE[m]->SetBinContent(b,0);
	fHistSpectrumSistSE[m]->SetBinError(b,0);
	if (TypeAnalysis!=2){
	  fHistSpectrumSistDeltaEta[m]->SetBinContent(b,0);
	  fHistSpectrumSistDeltaEta[m]->SetBinError(b,0);
	}
      }
      */
      if (TypeAnalysis!=0) YInf=0;
      canvasPlotProj[m]->cd(v+1);
      fHistPhiDistr_solostat[m][v]->GetYaxis()->SetRangeUser(YInf, YSup);
      fHistPhiDistr_solosist[m][v]->GetYaxis()->SetRangeUser(YInf, YSup);
      fHistPhiDistr_solostat[m][v]->Draw("same");
      fHistPhiDistr_solosist[m][v]->SetFillStyle(0);
      fHistPhiDistr_solosist[m][v]->Draw("same e2");
      if (TypeAnalysis==0) lineat1->Draw("same");
      
    } //end loop on pt v0
    //    fHistSpectrumStat[m]->Scale(1./NTrigger[m]);
    //    fHistSpectrumSist[m]->Scale(1./NTrigger[m]);
    fHistSpectrumStat[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSist[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistSE[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistSETuva[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistDCAz[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    fHistSpectrumSistPurity[m]->Scale(1./DeltaPhiWidth[sysPhi]);
    if (TypeAnalysis!=2)    fHistSpectrumSistDeltaEta[m]->Scale(1./DeltaPhiWidth[sysPhi]);

    fHistSpectrumStatRelError[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumStatRelError_"+Smolt[m]);
    fHistSpectrumSistRelError[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelError_"+Smolt[m]);
    fHistSpectrumSistRelErrorSE[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
    fHistSpectrumSistRelErrorSETuva[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorSETuva_"+Smolt[m]);
    fHistSpectrumSistRelErrorDeltaEta[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorDeltaEta_"+Smolt[m]);
    fHistSpectrumSistRelErrorDCAz[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorDCAz_"+Smolt[m]);
    fHistSpectrumSistRelErrorPurity[m]=(TH1D*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorPurity_"+Smolt[m]);

    //    for(Int_t  b=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
    for(Int_t  b=fHistSpectrumStat[m]->FindBin(LowRangeSpectrumPart[m]+0.001); b<=fHistSpectrumStat[m]->FindBin(UpRangeSpectrumPart[m]-0.001) ; b++){
      fHistSpectrumStatRelError[m]->SetBinContent(b,    fHistSpectrumStat[m]->GetBinError(b)/    fHistSpectrumStat[m] ->GetBinContent(b));
      fHistSpectrumSistRelError[m]->SetBinContent(b,    fHistSpectrumSist[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
      fHistSpectrumSistRelErrorSE[m]->SetBinContent(b,    fHistSpectrumSistSE[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
      fHistSpectrumSistRelErrorSETuva[m]->SetBinContent(b,    fHistSpectrumSistSETuva[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
      fHistSpectrumSistRelErrorDCAz[m]->SetBinContent(b,    fHistSpectrumSistDCAz[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
      fHistSpectrumSistRelErrorPurity[m]->SetBinContent(b,    fHistSpectrumSistPurity[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
      if (TypeAnalysis!=2)    fHistSpectrumSistRelErrorDeltaEta[m]->SetBinContent(b,    fHistSpectrumSistDeltaEta[m]->GetBinError(b)/    fHistSpectrumSist[m] ->GetBinContent(b));
      fHistSpectrumStatRelError[m]->SetBinError(b,0);
      fHistSpectrumSistRelError[m]->SetBinError(b,0);
      fHistSpectrumSistRelErrorSE[m]->SetBinError(b,0);
      fHistSpectrumSistRelErrorSETuva[m]->SetBinError(b,0);
      fHistSpectrumSistRelErrorDCAz[m]->SetBinError(b,0);
      fHistSpectrumSistRelErrorPurity[m]->SetBinError(b,0);
      fHistSpectrumSistRelErrorDeltaEta[m]->SetBinError(b,0);
    }

    canvasPtSpectra->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSist[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY, title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumStat[m], 0, LimSup, Color[TypeAnalysis], 6, titleX, titleY,  title+SmoltLegend[m]);
    fHistSpectrumSist[m]    ->SetFillStyle(0);
    fHistSpectrumSist[m]->Draw("same e2");
    fHistSpectrumStat[m]->Draw("same e");

    canvasPtSpectraRelErrorAll->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumStatRelError[m], 0, LimSupError, Color[TypeAnalysis], 33, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelError[m], 0, LimSupError, Color[TypeAnalysis], 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorSE[m], 0, LimSupError, 807, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorSETuva[m], 0, LimSupError, 630, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorDCAz[m], 0, LimSupError, 867, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorPurity[m], 0, LimSupError, 909, 27, titleX, titleYRel,  title+SmoltLegend[m]);
    StyleHisto(fHistSpectrumSistRelErrorDeltaEta[m], 0, LimSupError, 881, 27, titleX, titleYRel,  title+SmoltLegend[m]);

    if(m==0)    legendError->AddEntry(fHistSpectrumStatRelError[m], "stat.", "pl");
    if(m==0 && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelError[m], "syst.", "pl");
    if(m==0)    legendError->AddEntry(fHistSpectrumSistRelErrorSE[m], "syst. topol. sel. M1", "l");
    if(m==0 && SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorSETuva[m], "syst topol. sel. M2", "l");
    if(m==0 && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorDCAz[m], "syst. DCAzTrigger", "l");
    if(m==0 && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorPurity[m], "syst. purity", "l");
    if(m==0 && TypeAnalysis!=2 && !SysTuvaWay)    legendError->AddEntry(fHistSpectrumSistRelErrorDeltaEta[m], "syst. #Delta #eta region", "l");
    if (!SysTuvaWay)    fHistSpectrumSistRelError[m]->Draw("same p");
    fHistSpectrumSistRelErrorSE[m]->Draw("same");
    if (SysTuvaWay)    fHistSpectrumSistRelErrorSETuva[m]->Draw("same");
    if (!SysTuvaWay){
    fHistSpectrumSistRelErrorDCAz[m]->Draw("same");
    fHistSpectrumSistRelErrorPurity[m]->Draw("same");
    if (TypeAnalysis!=2)  fHistSpectrumSistRelErrorDeltaEta[m]->Draw("same");
    }
    fHistSpectrumStatRelError[m]->Draw("same p");
    legendError->Draw("");

    canvasPtSpectraRelError->cd(m+1);
    gPad->SetLeftMargin(0.15);
    fHistSpectrumSistRelError[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");
  }//end loop on m
  cout << " going to write on file " << endl;   
  fileout->WriteTObject(canvasPtSpectra);
  fileout->WriteTObject(canvasPtSpectraRelError);
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject(canvasPlotProj[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorSE[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorDeltaEta[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorDCAz[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorPurity[m]);
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fileout->WriteTObject(fHistPhiDistr_solostat[m][v]);
      fileout->WriteTObject(fHistPhiDistr_solosist[m][v]);
    }
    fileout->WriteTObject(      fHistSpectrumSist[m]);
    fileout->WriteTObject(      fHistSpectrumStat[m]);

  }

  TH1F * DeltaPhiLimit = new TH1F ("DeltaPhiLimit", "DeltaPhiLimit", 2,0,2);
  DeltaPhiLimit->SetBinContent(1, fHistPhiDistr_master[1][1]->GetXaxis()->GetBinLowEdge( fHistPhiDistr_master[1][1]->FindBin(ALowBin[sysPhi] )));
  DeltaPhiLimit->SetBinContent(2, fHistPhiDistr_master[1][1]->GetXaxis()->GetBinUpEdge( fHistPhiDistr_master[1][1]->FindBin(AUpBin[sysPhi] )));
  DeltaPhiLimit->GetXaxis()->SetBinLabel(1,"Low DeltaPhi");
  DeltaPhiLimit->GetXaxis()->SetBinLabel(2,"Up DeltaPhi");
  fileout->WriteTObject(  DeltaPhiLimit);

  fileout->Close();
  cout << "pt spectra value for m==5 " << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      cout << "v: " << v << " " << fHistSpectrumStat[5]->GetBinContent(v+1);
    }
  cout << "DeltaPhi width of one histo bin " << fHistPhiDistr_master[1][1]->GetBinWidth(1)<< endl;
  cout << "\nDeltaPhiWidth " << DeltaPhiWidth[sysPhi] << "("<< ALowBin[sysPhi]<<" ; "<<AUpBin[sysPhi] << ")"<<endl;
  cout <<  "effective range ("<<  fHistPhiDistr_master[1][1]->GetXaxis()->GetBinLowEdge( fHistPhiDistr_master[1][1]->FindBin(ALowBin[sysPhi] ))<<" ; "<< fHistPhiDistr_master[1][1]->GetXaxis()->GetBinUpEdge( fHistPhiDistr_master[1][1]->FindBin(AUpBin[sysPhi] ))<< ")"<<endl;
  cout << "starting from the files "  << PathInSignalSys << ", " << PathInDeltaEtaSys << " and " << PathInDCAzSys << endl;
  cout << " I have created the file " << stringout << endl;
}

