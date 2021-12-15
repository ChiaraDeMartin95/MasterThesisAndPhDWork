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

void StyleHisto(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString Title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetTitle(Title);
  histo->SetTitleSize(1.5);
}

void SysDeltaEta(Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Int_t TypeAnalysis=0, Bool_t isMC=0,   Int_t israp=0,TString year="161718Full_AOD234_hXi"/*"1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/,  TString Path1 =""/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, Int_t type=8,  Bool_t isEfficiency=1,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=0, Int_t sysTrigger=0, Int_t sys=0,/*leave it like this*/ Float_t PtTrigMin1=0.15, Bool_t SysTuvaWay=0, Bool_t isPreliminary=0, TString yearLowPtTrig=""){

  //  if (SESysTuvaWay) I calculate the sys uncertainty associated to topo selection taking |max yield - min yield|/2 among the values obtained with the tightest and the loosest selections

  if (TypeAnalysis>1 && !SysTuvaWay) {cout << "DeltaEta Sys errors not implemented for these regions " << endl; return;}


  if (ishhCorr && type!=0){
    type=0;
  }// {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  if (type==0){
    if (isPreliminary)    year = "1617_hK0s"; //for preliminary results
    else year = "1617_AOD234_hK0s";
    yearLowPtTrig= "";
  }
  else if (type==8){
    if (isPreliminary)  {
      year = "Run2DataRed_MECorr_hXi";
      yearLowPtTrig = "";
    }
    else   {
      year = "161718Full_AOD234_hXi";
      yearLowPtTrig = "_161718Full_AOD234_hXi";
    }
    PtBinning=0;
  }

  const   Int_t NumberTypeAnalysis=11;
  gStyle->SetOptStat(0);
  Int_t numsysDEta=6;
  if (SysTuvaWay) numsysDEta=5;
  Int_t minDEta=4;
  TString RegionType[3] ={"Jet", "Bulk", "All"} ;

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *fileinbisPart1;
  TFile *fileinbisPart2;
  TFile *filein;
  TFile *fileinOOJ;

  TFile *fileSignalSys;
  TFile *fileinDeltaEtaSys[numsysDEta+1];
  TFile *fileinDeltaEtaSysOOJ[numsysDEta+1];
  TString PathInSignalSys;
  TString PathInDeltaEtaSys[numsysDEta+1];
  TString PathInDeltaEtaSysOOJ[numsysDEta+1];
  Int_t NSign=5;

  TString PathIn0;
  TString PathIn1;
  TString PathInOOJ0;
  TString PathInOOJ1;
  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);
  file+=Path1;

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  TString JetOrBulk[NumberTypeAnalysis]={"Jet", "Bulk", "All", "AllNS", "AwaySide", "AllButJet", "JetBFit", "JetZYAM", "ASBFit", "ASZYAM", "JetNS"};
  TString DataOrMC[2]={"Data", "MC"};
  Int_t jet=TypeAnalysis;

  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

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
  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};
  TString SFit[3]={"exponential", "power law", "retta"};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[3] ={628, 418, 600};
  Int_t ColorsysDEta[12] = { 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};

  TH1D* fHistPhiDistr[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosist[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistSE[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_DeltaEta[nummolt+1][numPtV0][numsysDEta+1];
  TH1D* fHistPhiDistr_DeltaEtaRatio[nummolt+1][numPtV0][numsysDEta+1];
  TH1D* fHistPhiDistr_solosistDeltaEta[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_solosistDeltaEtaRelError[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_masterRelError[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_master[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_masterRebin[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_min[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_max[nummolt+1][numPtV0];

  TString stringout;
  stringout = Dir+"/DATA"+year0;
  if (!SysTuvaWay)  stringout +="/SysDeltaEta";
  else  stringout +="/LoosestTightestTopoSel";
  stringout +=  year;
  stringout += hhCorr[ishhCorr];
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  //  stringout+=   Form("_SysPhi%i_PtMin%.1f_", sysPhi, PtTrigMin);
  stringout+=   Form("_PtMin%.1f_",  PtTrigMin);
  stringout+= RegionType[TypeAnalysis];
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  //canvas definition
  TCanvas*  canvasDeltaEta[nummolt+1];
  TCanvas*  canvasBarlow[nummolt+1];
  TCanvas*  canvasProjSys[nummolt+1];
  TCanvas*  canvasDeltaEtaRatio[nummolt+1];
  TCanvas*  canvasRelError[nummolt+1];

  //histo definition
  TH1D*  hBarlowVar[nummolt+1][numPtV0][numsysDEta+1];

  //variables definition
  ///////////////////
  Float_t BarlowVar[nummolt+1][numPtV0][numsysDEta+1][100]={0};
  Int_t BarlowSign[nummolt+1][numPtV0][numsysDEta+1]={0};
  Bool_t BarlowPassed[nummolt+1][numPtV0][numsysDEta+1]={0};
  Float_t SysError[nummolt+1][numPtV0][100]={0};

  //old  TLegend *legend = new TLegend (0.5, 0.6, 0.9, 0.9);
  TLegend *legend = new TLegend (0.1, 0.2, 0.9, 0.9);
  legend->SetHeader("#Delta#eta selections");
  //  legend->SetTextSize(0.03);
  TString SDeltaEta[numsysDEta+1] ={""};
  TH1F *  fHistEtaLimitsOfRegion;
  Float_t LowEtaLimit[numsysDEta+1]={0};
  Float_t UpEtaLimit[numsysDEta+1]={0};
  Float_t LimSupRatio =1.1;
  Float_t LimInfRatio =0.9;
  if (TypeAnalysis==0) {
    LimSupRatio = 1.2;
    LimInfRatio = 0.8;
  }

  Float_t NTrigger[nummolt+1]={0};
  Float_t NTriggerAll=0;
  TString   PathInNTrig = Dir+"/DATA2016/SystematicAnalysis" +year;
  if (!isPreliminary && PtBinning>0) PathInNTrig += "_PtBinning1";
  PathInNTrig +=Path1; //it was without year                                                          
  if (isPreliminary)  PathInNTrig +="_Jet0.75";
  PathInNTrig +="_"+tipo[type];
  PathInNTrig +=Srap[israp];
  PathInNTrig +=SSkipAssoc[SkipAssoc];
  PathInNTrig +=  Form("_JetData_PtMin%.1f.root", PtTrigMin);
  TFile *fileNTrig = new TFile(PathInNTrig, "");
  if (!fileNTrig) return;
  TH1F* NTriggerHisto = (TH1F*)  fileNTrig->Get("NTriggerHisto");

  Int_t PtV0Min=1;
  //  if (!ishhCorr && type==0) PtV0Min=0;
  Float_t YSup=0.01;
  Float_t YInf=-0.001;
  if (type==8) {
    YSup = 0.0006;
    YInf = -0.0002;
  }

  Int_t sysV0=0;    
  if(isMC && isEfficiency) file = year + "_MCEff" +Path1;

  TF1* lineat1= new TF1("pol0", "pol0", -TMath::Pi()/2, TMath::Pi()*3./2);
  lineat1->SetLineColor(kBlue);
  lineat1->FixParameter(0,1);
  TF1* lineat0= new TF1("pol0", "pol0At0", -TMath::Pi()/2, TMath::Pi()*3./2);
  lineat0->SetLineColor(kBlack);
  lineat0->FixParameter(0,0);

  for(Int_t m=0; m<nummolt+1; m++){
    //    if (m!=nummolt) continue;
    cout << "\n\n\e[35mAnalysing multiplicity... " << Smolt[m] << "\e[39m" <<  endl;
    if (m!=nummolt) {
      NTrigger[m] =    NTriggerHisto->GetBinContent(NTriggerHisto->FindBin(Nmolt[m]+0.0001));
      NTriggerAll +=    NTriggerHisto->GetBinContent(NTriggerHisto->FindBin(Nmolt[m]+0.0001));
    }
    else NTrigger[m] = NTriggerAll;
    //    cout << "NTrigger= " << NTrigger[m] << endl;

    canvasDeltaEta[m]= new TCanvas(Form("canvasProj_m%i",m), Form("canvasProj_m%i",m), 1200, 800);
    canvasDeltaEta[m]->Divide((float)numPtV0/2+1, 2);
    
    canvasBarlow[m]= new TCanvas(Form("canvasBarlowVar_m%i",m), Form("canvasBarlowVar_m%i",m), 1200, 800);
    canvasBarlow[m]->Divide((float)numPtV0/2+1, 2);

    canvasProjSys[m]= new TCanvas(Form("canvasProjSys_m%i",m), Form("canvasProjSys_m%i",m), 1200, 800);
    canvasProjSys[m]->Divide((float)numPtV0/2+1, 2);

    canvasDeltaEtaRatio[m]= new TCanvas(Form("canvasProjRatio_m%i",m), Form("canvasProjRatio_m%i",m), 1200, 800);
    canvasDeltaEtaRatio[m]->Divide((float)numPtV0/2+1, 2);

    canvasRelError[m]= new TCanvas(Form("canvasRelError_m%i",m), Form("canvasRelError_m%i",m), 1200, 800);
    canvasRelError[m]->Divide((float)numPtV0/2+1, 2);

    //PathIn1 definition
    if (isMC && !isEfficiency) {
      PathIn1=Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file+  Form("_Sys%i", sys)+"_Output";
      if (!isPreliminary)       PathIn1 += "_IsEtaEff"; 
      PathIn1 += ".root";
      PathInOOJ1="OOJComparison"+file+Form("_Sys%i",sys)+"_Output.root";
    }
    else {
      PathIn1 = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file;
      PathInOOJ1= "OOJComparison" + file;
      PathInOOJ1 += yearLowPtTrig;
      if(type>=0){
	PathIn1 +="_"+tipo[type];
	PathIn1 +=Srap[israp];
	PathIn1 +=SSkipAssoc[SkipAssoc];
	PathInOOJ1 +="_"+tipo[type];
	PathInOOJ1 +=Srap[israp];
	PathInOOJ1 +=SSkipAssoc[SkipAssoc];
      }
      PathIn0=PathIn1;
      PathInOOJ0=PathInOOJ1;
      PathIn1+= hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_Output";
      PathInOOJ1+= hhCorr[ishhCorr]+Form("_sys%i_PtTrigMin%.1f_PtTrigMin%.1f",  sys, PtTrigMin, PtTrigMin1)+"_Output";
      if (!isPreliminary)       PathIn1 += "_IsEtaEff"; 
      if (!isPreliminary)       PathInOOJ1 += "_IsEtaEff"; 
      PathIn1+=".root";
      PathInOOJ1+=".root";
    }
    cout << "\e[36mInput file: \e[39m" << PathIn1 << endl;
    if (TypeAnalysis==0 && type==8) cout << "\e[36mInput file for jet distribution: \e[39m" << PathInOOJ1 << endl;

    filein = new TFile(PathIn1, "");
    if (!filein) {cout << filein << "does not exist " << endl ; return; }

    if (type==8 && TypeAnalysis==0){
      fileinOOJ = new TFile(PathInOOJ1, "");
      if (!fileinOOJ) {cout << fileinOOJ << "does not exist " << endl ; return; }
    }    

    Int_t sysDEtaCorr=0;
    for (Int_t sysDEta=0; sysDEta<=numsysDEta; sysDEta++){
      if (isPreliminary && sysDEta==5)      sysDEtaCorr=0;
      else  sysDEtaCorr=sysDEta;
      if (sysDEta>0 && sysDEta<4) continue;
      if (type==8 && sysDEta==5 && TypeAnalysis==1 &&!SysTuvaWay) continue;
      cout << "sysDEta " <<       sysDEta << endl;
      //      if (sysDEta>numsysDEta) return;
      if (!SysTuvaWay || sysDEta==0)   {
	PathInDeltaEtaSys[sysDEta] = PathIn0 + hhCorr[ishhCorr]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysDEtaCorr, PtTrigMin)+"_Output";
	PathInDeltaEtaSysOOJ[sysDEta] = PathInOOJ0 + hhCorr[ishhCorr]+Form("_sys%i_PtTrigMin%.1f_PtTrigMin%.1f", sysDEta, PtTrigMin, PtTrigMin1)+"_Output";
	if (!isPreliminary)	PathInDeltaEtaSys[sysDEta] += "_IsEtaEff";
	if (!isPreliminary)	PathInDeltaEtaSysOOJ[sysDEta] += "_IsEtaEff";
	PathInDeltaEtaSys[sysDEta] += ".root";
	PathInDeltaEtaSysOOJ[sysDEta] += ".root";
      }
      else{
	if (sysDEta==4){
	  PathInDeltaEtaSys[sysDEta] = PathIn0 + hhCorr[ishhCorr]+Form("_SysT%i_SysV0Loosest_Sys%i_PtMin%.1f", sysTrigger,0, PtTrigMin)+"_Output.root";
	  PathInDeltaEtaSysOOJ[sysDEta] = PathInOOJ0 + hhCorr[ishhCorr]+Form("_sysV0Loosest_PtTrigMin%.1f_PtTrigMin%.1f", PtTrigMin, PtTrigMin1)+"_Output.root";
	}
	else	if (sysDEta==5){
	  PathInDeltaEtaSys[sysDEta] = PathIn0 + hhCorr[ishhCorr]+Form("_SysT%i_SysV0Tightest_Sys%i_PtMin%.1f", sysTrigger,0, PtTrigMin)+"_Output.root";
	  PathInDeltaEtaSysOOJ[sysDEta] = PathInOOJ0 + hhCorr[ishhCorr]+Form("_sysV0Tightest_PtTrigMin%.1f_PtTrigMin%.1f", PtTrigMin, PtTrigMin1)+"_Output.root";
	}
      }

      fileinDeltaEtaSys[sysDEta] = new TFile(PathInDeltaEtaSys[sysDEta], "");
      if (type==8){
	fileinDeltaEtaSysOOJ[sysDEta] = new TFile(PathInDeltaEtaSysOOJ[sysDEta], "");
      }

      if (!fileinDeltaEtaSys[sysDEta]) {cout << PathInDeltaEtaSys[sysDEta] << "does not exist " << endl ; return;}
      if (type==8 && TypeAnalysis ==0 && !fileinDeltaEtaSysOOJ[sysDEta]) {cout << PathInDeltaEtaSysOOJ[sysDEta] << "does not exist " << endl ; return;}

      fHistEtaLimitsOfRegion=(TH1F*)       fileinDeltaEtaSys[sysDEta]->Get("fHistEtaLimitsOfRegion");
      if (TypeAnalysis==0){
	LowEtaLimit[sysDEta]=      fHistEtaLimitsOfRegion->GetBinContent(1);
	UpEtaLimit[sysDEta]=      fHistEtaLimitsOfRegion->GetBinContent(2);
      }
      else {
	LowEtaLimit[sysDEta]=      fHistEtaLimitsOfRegion->GetBinContent(3);
	UpEtaLimit[sysDEta]=      fHistEtaLimitsOfRegion->GetBinContent(4);
      }
      cout << "I get DeltaEta ranges from files: " << LowEtaLimit[sysDEta] <<"-"<<UpEtaLimit[sysDEta]<< endl;
    }

    cout << "\n\n\e[35mI got the files... Starting the loop on pT\e[39m" << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      cout << "\n\n----> Pt interval: " << SPtV0[v] << " GeV/c " << endl;
      if (NPtV0[v] < 2.5-0.001 && TypeAnalysis==0 && type==8){
	if (isPreliminary)	fHistPhiDistr_master[m][v]=(TH1D*)fileinOOJ->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmoothBis");
	else 	fHistPhiDistr_master[m][v]=(TH1D*)fileinOOJ->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmooth");
      }
      else if (TypeAnalysis==0)      fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
      if (TypeAnalysis==1) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");
      if (TypeAnalysis==2) fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr");
      if(!fHistPhiDistr_master[m][v]) { cout << "histo " << "ME_m"<<Smolt[m]<<"_v"<<SPtV0[v]<<"_AC_phi_V0Sub_Bulk_EffCorr" <<" not found " << endl; return;}
      fHistPhiDistr_master[m][v]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);
      cout << "I got default histo " << endl;

      for (Int_t sysDEta=minDEta; sysDEta<=numsysDEta; sysDEta++){
	if (type==8 && sysDEta==5 &&TypeAnalysis==1 && !SysTuvaWay) continue;
	if (NPtV0[v] < 2.5-0.001 && TypeAnalysis==0 && type==8){
	  if (isPreliminary)	  fHistPhiDistr_DeltaEta[m][v][sysDEta]=(TH1D*)fileinDeltaEtaSysOOJ[sysDEta]->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmoothBis");
	  else 	  fHistPhiDistr_DeltaEta[m][v][sysDEta]=(TH1D*)fileinDeltaEtaSysOOJ[sysDEta]->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmooth");
	}
	else if (TypeAnalysis==0){
	  if (type==0) fHistPhiDistr_DeltaEta[m][v][sysDEta]=(TH1D*)fileinDeltaEtaSys[sysDEta]->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	  else fHistPhiDistr_DeltaEta[m][v][sysDEta]=(TH1D*)fileinDeltaEtaSysOOJ[sysDEta]->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr_DefaultMethod");
	}
	if (TypeAnalysis==1) fHistPhiDistr_DeltaEta[m][v][sysDEta]=(TH1D*)fileinDeltaEtaSys[sysDEta]->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");
	if (TypeAnalysis==2) fHistPhiDistr_DeltaEta[m][v][sysDEta]=(TH1D*)fileinDeltaEtaSys[sysDEta]->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_JetBulkEffCorr");
	if(!fHistPhiDistr_DeltaEta[m][v][sysDEta]) {cout << "histo fHistPhiDistr_DeltaEta for sysDEta " << sysDEta << " not present. It sould have been in  " <<fileinDeltaEtaSysOOJ[sysDEta]->GetName() << endl; return;}
	fHistPhiDistr_DeltaEta[m][v][sysDEta]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]+Form("_sysDEta%i",sysDEta));
	cout << "I got histo for sysDEta " << sysDEta << endl;
      }

      if (!(type==8 && (NPtV0[v] <2.5-0.001) && TypeAnalysis==0))       fHistPhiDistr_master[m][v]->Scale(1./NTrigger[m]);
      fHistPhiDistr_masterRelError[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_masterRelError_m"+Smolt[m]+"_v"+SPtV0[v]);
      if (!SysTuvaWay)      fHistPhiDistr_solosistDeltaEta[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistDeltaEta_m"+Smolt[m]+"_v"+SPtV0[v]);
      else       fHistPhiDistr_solosistDeltaEta[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistSETuvaWay_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosistDeltaEtaRelError[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("PhiDistr_solosistDeltaEtaRelError_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_min[m][v] =(TH1D*)      fHistPhiDistr_master[m][v]->Clone("PhiDistr_min_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_max[m][v] =(TH1D*)      fHistPhiDistr_master[m][v]->Clone("PhiDistr_max_m"+Smolt[m]+"_v"+SPtV0[v]);
      
      //1. I insert different deltaphi proj in a canvas (jey + OJ together)
      canvasDeltaEta[m]->cd(v+1);
      StyleHisto(fHistPhiDistr_master[m][v], YInf, YSup, Color[TypeAnalysis] ,1, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
      if (!SysTuvaWay)       SDeltaEta[0] =Form("%.2f-%.2f", LowEtaLimit[0], UpEtaLimit[0]);
      else    SDeltaEta[0] ="default sel";
      if (m==nummolt && v==PtV0Min)	legend->AddEntry(	fHistPhiDistr_master[m][v], SDeltaEta[0], "pl");
      //      cout << "ciao chiara.." << SDeltaEta[0] << endl;
      fHistPhiDistr_master[m][v]->DrawClone("");
      
      for (Int_t sysDEta=minDEta; sysDEta<=numsysDEta; sysDEta++){
	if (type==8 && sysDEta==5 &&TypeAnalysis==1 && !SysTuvaWay) continue;
	if (!SysTuvaWay) SDeltaEta[sysDEta] =Form("%.2f-%.2f", LowEtaLimit[sysDEta], UpEtaLimit[sysDEta]);
	else{
	  if (sysDEta==4)	  SDeltaEta[sysDEta] ="Loosest sel";
	  else if (sysDEta==5) 	  SDeltaEta[sysDEta] ="Tightest sel";
	}
	//cout << 	SDeltaEta[sysDEta] << endl;
	fHistPhiDistr_DeltaEta[m][v][sysDEta]->Sumw2();
	if (!(type==8 && TypeAnalysis==0))  	fHistPhiDistr_DeltaEta[m][v][sysDEta]->Scale(1./NTrigger[m]);
	StyleHisto(	fHistPhiDistr_DeltaEta[m][v][sysDEta], YInf, YSup, ColorsysDEta[sysDEta] ,1, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	if (m==nummolt && v==PtV0Min)	legend->AddEntry(	fHistPhiDistr_DeltaEta[m][v][sysDEta], SDeltaEta[sysDEta], "pl");
	fHistPhiDistr_DeltaEta[m][v][sysDEta]->DrawClone("same");
	fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta] = (TH1D*) fHistPhiDistr_DeltaEta[m][v][sysDEta]->Clone("PhiDistr_DeltaEtaRatio_m"+Smolt[m]+"_v"+SPtV0[v]);
      }
      
      //      legend->Draw("");

      canvasDeltaEta[m]->cd(12);
      legend->Draw("");
      
      cout << "\n\e[35m2. I evaluate if they are Barlow significant (and plot sigmaBarlow), are they correlated? \e[39m " << endl;
      for (Int_t sysDEta=minDEta; sysDEta<=numsysDEta; sysDEta++){
	if (type==8 && sysDEta==5 &&TypeAnalysis==1 && !SysTuvaWay) continue;
	if (fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetNbinsX()==26) NSign = 5;
	else if (fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetNbinsX()==52) NSign =8;
	else {cout << "NSign need to be redefined " << endl; return;}
	if (TypeAnalysis==0) NSign =4;
	BarlowSign[m][v][sysDEta]=0;
	BarlowPassed[m][v][sysDEta]=0;
	hBarlowVar[m][v][sysDEta]=(TH1D*) 	  fHistPhiDistr_DeltaEta[m][v][sysDEta]->Clone("BarlowVar_m"+Smolt[m]+"_v"+SPtV0[v]);
	cout << "nbins " << fHistPhiDistr_master[m][v]->GetNbinsX()<< endl;
	for (Int_t dphi=0; dphi<fHistPhiDistr_master[m][v]->GetNbinsX(); dphi++){
	  BarlowVar[m][v][sysDEta][dphi]=( fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetBinContent(dphi+1) -  fHistPhiDistr_master[m][v]->GetBinContent(dphi+1))/sqrt( TMath::Abs(pow(fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetBinError(dphi+1),2) -  pow(fHistPhiDistr_master[m][v]->GetBinError(dphi+1),2)));
	  hBarlowVar[m][v][sysDEta] ->SetBinContent(dphi+1, BarlowVar[m][v][sysDEta][dphi]);
 	  hBarlowVar[m][v][sysDEta] ->SetBinError(dphi+1, 0);
	  if (TypeAnalysis==0 && TMath::Abs(fHistPhiDistr_master[m][v]->GetBinCenter(dphi+1)) > 1.3) continue;
	  if (TMath::Abs(BarlowVar[m][v][sysDEta][dphi])>2) BarlowSign[m][v][sysDEta]++;
	  //	  cout << " BarlowContatore " << BarlowSign[m][v][sysDEta]<< endl;
	}

	if (BarlowSign[m][v][sysDEta]>NSign) BarlowPassed[m][v][sysDEta]=1;
	cout << "NPtV0 " << NPtV0[v] <<" SysDEta " << sysDEta << " " <<  BarlowSign[m][v][sysDEta]<< " ->Barlow Passed (yes =1) : " << BarlowPassed[m][v][sysDEta]<< endl;
	StyleHisto(hBarlowVar[m][v][sysDEta], -5, 5, ColorsysDEta[sysDEta], 33, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	if ( BarlowPassed[m][v][sysDEta]) 	StyleHisto(hBarlowVar[m][v][sysDEta], -5, 5, ColorsysDEta[sysDEta], 27, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	canvasBarlow[m]->cd(v+1);
	hBarlowVar[m][v][sysDEta] ->Draw("same p");
	lineat0->Draw("same");

	canvasBarlow[m]->cd(12);
	legend->Draw("");

	canvasDeltaEtaRatio[m]->cd(v+1);
	if (sysDEta==minDEta){
	  fHistPhiDistr_masterRebin[m][v] = (TH1D*) fHistPhiDistr_master[m][v]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_Rebin");
	  fHistPhiDistr_masterRebin[m][v]->Rebin(2);
	}
	fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta]->Rebin(2);
	fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta] ->Divide(fHistPhiDistr_masterRebin[m][v]);
	//	cout << fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta]->GetNbinsX() << endl;
	//	cout << fHistPhiDistr_masterRebin[m][v]->GetNbinsX() << endl;
	if (fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta]->GetNbinsX() != fHistPhiDistr_masterRebin[m][v]->GetNbinsX())	return;
	for (Int_t dphi=0; dphi<fHistPhiDistr_master[m][v]->GetNbinsX(); dphi++){
	  fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta] ->SetBinError(dphi+1, 0);
	}
	if (BarlowPassed[m][v][sysDEta]) StyleHisto(fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta], LimInfRatio, LimSupRatio, ColorsysDEta[sysDEta], 33, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	else StyleHisto(fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta], LimInfRatio, LimSupRatio, ColorsysDEta[sysDEta], 27, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	if (TypeAnalysis==0) fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta]->GetXaxis()->SetRangeUser(-1.3, 1.3);
	//	if (BarlowPassed[m][v][sysDEta]) {
	fHistPhiDistr_DeltaEtaRatio[m][v][sysDEta]->Draw("same p");
	lineat1->DrawClone("same");
	//	}
      }//sysDEta

      canvasDeltaEtaRatio[m]->cd(12);
      legend->Draw("");
      

      cout << "\n\e[35m3. If Barlow significant, I calculate the systematic error associated \e[39m" << endl;
      Int_t IsSystematic=0;
      for (Int_t sysDEta=minDEta; sysDEta<=numsysDEta; sysDEta++){
	if (type==8 && sysDEta==5 &&TypeAnalysis==1 && !SysTuvaWay) continue;
	if (!BarlowPassed[m][v][sysDEta]) continue;
	IsSystematic ++;
	for (Int_t dphi=0; dphi<fHistPhiDistr_master[m][v]->GetNbinsX(); dphi++){
	  if( fHistPhiDistr_min[m][v]->GetBinContent(dphi+1) >	  fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetBinContent(dphi+1)){
	    fHistPhiDistr_min[m][v]->SetBinContent(dphi+1, fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetBinContent(dphi+1));
	  }
	  if( fHistPhiDistr_max[m][v]->GetBinContent(dphi+1) <	  fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetBinContent(dphi+1)){
	    fHistPhiDistr_max[m][v]->SetBinContent(dphi+1, fHistPhiDistr_DeltaEta[m][v][sysDEta]->GetBinContent(dphi+1));
	  }
	}//phi
      }//sysDEta
      canvasProjSys[m]->cd(v+1);
      StyleHisto(      fHistPhiDistr_min[m][v], YInf, YSup, kRed, 1, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
      //      fHistPhiDistr_min[m][v]->Draw("same");
      StyleHisto(      fHistPhiDistr_max[m][v], YInf, YSup, kBlue, 1, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
      //      fHistPhiDistr_max[m][v]->Draw("same");

      for (Int_t dphi=0; dphi<fHistPhiDistr_master[m][v]->GetNbinsX(); dphi++){
	SysError[m][v][dphi] = (fHistPhiDistr_max[m][v]->GetBinContent(dphi+1) - fHistPhiDistr_min[m][v]->GetBinContent(dphi+1))/2;
	/*	
		cout << "\n" << endl;
		cout << fHistPhiDistr_max[m][v]->GetBinContent(dphi+1) << " >> " << fHistPhiDistr_master[m][v]->GetBinContent(dphi+1) << " >> " << fHistPhiDistr_min[m][v]->GetBinContent(dphi+1) << endl;
		cout << fHistPhiDistr_master[m][v]->GetBinContent(dphi+1) << " var " << fHistPhiDistr_DeltaEta[m][v][5]->GetBinContent(dphi+1) << endl;
		cout << " sys error " << 	SysError[m][v][dphi]<< endl;
		cout << "relative sys serror " << 	SysError[m][v][dphi]/fHistPhiDistr_master[m][v]->GetBinContent(dphi+1) << endl;
	*/
	if (!IsSystematic) 	fHistPhiDistr_solosistDeltaEta[m][v]->SetBinError(dphi+1, 0);
	else	fHistPhiDistr_solosistDeltaEta[m][v]->SetBinError(dphi+1, SysError[m][v][dphi]);
      }
      canvasProjSys[m]->cd(v+1);
      fHistPhiDistr_master[m][v]->Draw("same");
      StyleHisto(      fHistPhiDistr_solosistDeltaEta[m][v], YInf, YSup, Color[TypeAnalysis], 1, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
      fHistPhiDistr_solosistDeltaEta[m][v]->SetFillStyle(0);
      fHistPhiDistr_solosistDeltaEta[m][v]->Draw("same e2");

      canvasRelError[m]->cd(v+1);
      for (Int_t dphi=0; dphi<fHistPhiDistr_master[m][v]->GetNbinsX(); dphi++){
	if (TypeAnalysis==1){
	  fHistPhiDistr_solosistDeltaEtaRelError[m][v]->SetBinContent(dphi+1,       fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi+1)/     fHistPhiDistr_solosistDeltaEta[m][v]->GetBinContent(dphi+1));
	  fHistPhiDistr_solosistDeltaEtaRelError[m][v]->SetBinError(dphi+1, 0);
	  fHistPhiDistr_masterRelError[m][v]->SetBinContent(dphi+1,       fHistPhiDistr_master[m][v]->GetBinError(dphi+1)/     fHistPhiDistr_master[m][v]->GetBinContent(dphi+1));
	  fHistPhiDistr_masterRelError[m][v]->SetBinError(dphi+1,0);
	}
	else {
	  fHistPhiDistr_solosistDeltaEtaRelError[m][v]->SetBinContent(dphi+1,       fHistPhiDistr_solosistDeltaEta[m][v]->GetBinError(dphi+1));
	  fHistPhiDistr_solosistDeltaEtaRelError[m][v]->SetBinError(dphi+1, 0);
	  fHistPhiDistr_masterRelError[m][v]->SetBinContent(dphi+1,       fHistPhiDistr_master[m][v]->GetBinError(dphi+1));
	  fHistPhiDistr_masterRelError[m][v]->SetBinError(dphi+1,0);
	}
	//	cout << "stat rel error" << fHistPhiDistr_masterRelError[m][v]->GetBinContent(dphi+1)<<endl;
	//	cout << "syst rel error" << fHistPhiDistr_solosistDeltaEtaRelError[m][v]->GetBinContent(dphi+1)<<endl;
      }
      if (TypeAnalysis==1){
	StyleHisto(      fHistPhiDistr_solosistDeltaEtaRelError[m][v], 0, 0.3,Color[TypeAnalysis], 27,Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	StyleHisto(      fHistPhiDistr_masterRelError[m][v]  , 0, 0.3,Color[TypeAnalysis], 33, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
      }
      else {
	StyleHisto(      fHistPhiDistr_solosistDeltaEtaRelError[m][v], 0, YSup/10,Color[TypeAnalysis], 27, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	StyleHisto(      fHistPhiDistr_masterRelError[m][v]  , 0, YSup/10,Color[TypeAnalysis], 33, Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
      }
      fHistPhiDistr_solosistDeltaEtaRelError[m][v]->Draw("same p");
      fHistPhiDistr_masterRelError[m][v]  ->Draw("same p");
    } //end loop v
  } //end loop m

  cout << "\n\e[35mGoing to write on file...\e[39m"<< endl;  
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject(canvasDeltaEta[m]);
    fileout->WriteTObject(canvasBarlow[m]);
    fileout->WriteTObject(canvasProjSys[m]);
    fileout->WriteTObject(canvasDeltaEtaRatio[m]);
    fileout->WriteTObject(canvasRelError[m]);
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      fileout->WriteTObject(fHistPhiDistr_solosistDeltaEta[m][v]);
      fileout->WriteTObject(fHistPhiDistr_master[m][v]);
      for (Int_t sysDEta=minDEta; sysDEta<=numsysDEta; sysDEta++){
	if (BarlowPassed[m][v][sysDEta]) {
	  cout << " m " << m << " v " << v << " sysdEta " << sysDEta << " barlow passed? " << BarlowPassed[m][v][sysDEta]<< endl;
	}
      }
    } //end loop on pt v0
  }//end loop on m
  
  fileout->Close();

  cout << "\n\nDefault limit: " << LowEtaLimit[0] << "- " <<     UpEtaLimit[0]<< "\n" <<endl;
  for (Int_t sysDEta=minDEta; sysDEta<=numsysDEta; sysDEta++){
    if (type==8 && sysDEta==5 &&TypeAnalysis==1 && !SysTuvaWay) continue;
    cout << LowEtaLimit[sysDEta] << "- " <<     UpEtaLimit[sysDEta]<< endl;
    if (TypeAnalysis==0 && type==8) cout << "Input file: " <<fileinDeltaEtaSysOOJ[sysDEta]-> GetName()<< "\n"<<endl;
    else cout << "Input file: " <<fileinDeltaEtaSys[sysDEta]-> GetName()<< "\n"<<endl;
  }

  cout << "\n\e[35mI have created the file:\e[39m " << stringout << endl;
  cout << "The effect is considered significant if |Barlow variable| > 2 for more than " << NSign << " out of ";
  if (TypeAnalysis==0) {
    Int_t counter=0;
    for (Int_t n=1; n<= fHistPhiDistr_DeltaEta[0][1][minDEta]->GetNbinsX(); n++){
      if (TMath::Abs(fHistPhiDistr_master[0][1]->GetBinCenter(n)) <= 1.3) counter++;
    }
    cout << counter << " bins" << endl;
  }
  else cout << fHistPhiDistr_DeltaEta[0][1][minDEta]->GetNbinsX() << " bins" << endl;
}

