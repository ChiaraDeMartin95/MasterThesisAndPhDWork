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
#include <TFile.h>
#include <TLegend.h>

void SignalExtractionSys(Int_t numsysV0index =100, Int_t sysTrigger=0, Bool_t ishhCorr=0, Bool_t SkipAssoc=1,Float_t PtTrigMin=3,  Float_t PtTrigMinFit=3, Int_t sysV0=0, Int_t sys=0,Int_t type=0,  Int_t israp=0,Bool_t isMC=0, Bool_t isEfficiency=1,   TString year0 = "2016",TString year="1617_AOD234_hK0s"/*"1617_hK0s"/*"2016kehjl_hK0s"/"AllMC_hXi"/*"2018f1_extra_hK0s"/*"2016k_hK0s"/*"2016k_MECorr"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"*/, TString yearMC="1617_GP_AOD235_With18c12b"/*"1617MC_hK0s"/*"AllMC_hXi"/*"1617GP_hXi"/*"2016kl_hXi"/*"2018f1_extra_MECorr"/"2018f1_extra_hK0s"/*"2018f1_extra_hK0s_30runs_150MeV"*/,  TString Path1 =""/*"_PtTrigMax2.5"/*"_NewMultClassBis"*/, TString Dir ="FinalOutput/",  Float_t ptjmax=15, Int_t rebin=2,  Int_t rebinx=2,  Int_t rebiny=2, Bool_t MasterThesisAnalysis=0, Bool_t isEnlargedDeltaEtaPhi=0,Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Int_t MultBinning=0, Int_t PtBinning=1, Float_t PtTrigMin2=0.15, Bool_t isForPreliminary=0){ 

  //  gStyle->SetOptStat(0);

  if (type==8){
    year = "Run2DataRed_MECorr_hXi";
    yearMC = "AllMC_hXi";
    PtBinning=0;
  }
  else if (type==0){
    if (isForPreliminary){  //for Preliminary results
      year = "1617_hK0s";
      yearMC = "1617MC_hK0s";
    }
    year = "1617_AOD234_hK0s";
    yearMC = "1617_GP_AOD235";
    PtBinning=1;
  }
  const Int_t numtipo=10;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  Dir+="DATA"+year0;
  TString file = year;
  if (isMC) file += "_MCEff";
  //  file +=Path1; 
  //  if(isMC && isEfficiency) file = year + "_MCEff";
  if(isMC && !isEfficiency) file = yearMC + "_MCTruth";

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  const Int_t numPtTrigger=1;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
  Int_t PtV0Min=1; //el 1
  if (!ishhCorr && type==0) PtV0Min=0;

  Int_t Marker[nummolt+1]={7,4,20,22,29,28};
  Int_t Color[nummolt+1]={2,3,4,6,7, 9};
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 

  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Smolt[m] = Smolt0[m];
      Nmolt[m] = Nmolt0[m];
    }
    if (MultBinning==1){
      Smolt[m] = Smolt1[m];
      Nmolt[m] = Nmolt1[m];
    }
    if (MultBinning==2){
      Smolt[m] = Smolt2[m];
      Nmolt[m] = Nmolt2[m];
    }

  }

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  //  TString SPtV0[numPtV0]={"", "","0.5-1", "1-1.5","1.5-2", "2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};
  //TString SNPtV0[numPtV0+1]={"0.0","0.0", "0.5","1.0","1.5","2.0","3.0","4.0","8.0"};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }

  TString SPtV01[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};
  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1) numPtV0Max = numPtV0;
  else numPtV0Max = numPtV0-1;
  if (PtBinning==1){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }

  TString PathInFirstBis;
  TString PathInFirst="./FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";
  PathInFirst+=year;
  if(isMC && isEfficiency){
    PathInFirst+="_MCEff";
  }
  if(isMC && !isEfficiency){
    PathInFirst+="_MCTruth";
  }
  PathInFirst+=Path1;
  if (PtBinning>0)  PathInFirst+=Form("_PtBinning%i",PtBinning);
  PathInFirst+="_";
  if (!ishhCorr){
    PathInFirst +=tipo[type];
  }
  PathInFirst +=Srap[israp];
  if (!SkipAssoc)  PathInFirst +="_AllAssoc";
  if (ishhCorr) PathInFirst +="_hhCorr";

  TString PathIn;
  TString PathInDef;
  TFile *filein;
  TFile *fileinOOJNew;
  TFile *fileinFirst;
  TFile *fileinDef;

  TString   PathInNTrig = Dir+"/SystematicAnalysis" +year;
  if (!isForPreliminary) PathInNTrig +=Form("_PtBinning%i", PtBinning);
  PathInNTrig +=Path1; //it was without year                                                       
  if (isForPreliminary)  PathInNTrig +="_Jet0.75";
  PathInNTrig +="_"+tipo[type];
  PathInNTrig +=Srap[israp];
  PathInNTrig +=SSkipAssoc[SkipAssoc];
  PathInNTrig +=  Form("_JetData_PtMin%.1f.root", PtTrigMin);
  TFile *fileNTrig = new TFile(PathInNTrig, "");
  if (!fileNTrig) return;
  TH1F* NTriggerHisto = (TH1F*)  fileNTrig->Get("NTriggerHisto");
  if (!NTriggerHisto) {cout << "ntrigger histo not present " << endl; return;}

  TString PathOut=Dir+"/histo/SignalExtractionStudy"  +file+Path1 ;
  if (PtBinning>0)    PathOut +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      PathOut +="_"+tipo[type];
    PathOut +=Srap[israp];
    PathOut +=SSkipAssoc[SkipAssoc];
  }
  if (sysTrigger==0)  PathOut+=   Form("_SysT%i_SysV0Num%i_Sys%i_PtMin%.1f.root", 0, numsysV0index, sys, PtTrigMin);
  else   PathOut+=   Form("_SysTNum%i_SysV0%i_Sys%i_PtMin%.1f.root", numsysV0index, 0, sys, PtTrigMin);

  TFile *fileout = new TFile (PathOut , "RECREATE");

  const Int_t numRegions=3;
  TString SRegion[numRegions] = {"Jet", "Bulk", "Inclusive"};
  TString SRegionInput[numRegions] = {"BulkSub_", "Bulk_", "JetBulk"};
  Int_t ColorRegion[numRegions] = {628, 413, 600};
  Int_t ColorRegionBis[numRegions] = {634, 418, 602};
  TH1F * hSpread[nummolt+1][numPtV0][60][numRegions]; 
  TH1F * hDeltaPhi[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiSist[nummolt+1][numPtV0][numRegions][numsysV0index]; 
  TH1F * hDeltaPhiDef[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiDefClone[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiMean[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiDefSys[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiRatio[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiRelStat[nummolt+1][numPtV0][numRegions]; 
  TH1F * hDeltaPhiRelSist[nummolt+1][numPtV0][numRegions]; 
  Float_t hYieldDef[nummolt+1][numPtV0][60][numRegions]; 
  Float_t hYieldEStatDef[nummolt+1][numPtV0][60][numRegions]; 
  Float_t hYield[nummolt+1][numPtV0][60][numRegions]; 
  Float_t hYieldMean[nummolt+1][numPtV0][60][numRegions]; 
  Float_t hYieldSpread[nummolt+1][numPtV0][60][numRegions]; 
  TF1 * lineat1 = new TF1("pol0", "pol0", -TMath::Pi()/2, 3./2*TMath::Pi());
  lineat1->FixParameter(0,1);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineWidth(0.1);
  TF1 * lineat0 = new TF1("pol0", "pol0", -TMath::Pi()/2, 3./2*TMath::Pi());
  lineat0->FixParameter(0,0);
  lineat0->SetLineColor(kBlack);
  lineat0->SetLineWidth(0.1);


  TCanvas *canvasPlotProj[nummolt+1];
  TCanvas *canvasPlotProjAllSist[nummolt+1];
  TCanvas *canvasRelError[nummolt+1];
  TCanvas *canvasPlotRatio[nummolt+1];
  TCanvas *canvasSys[nummolt+1][numPtV0][3];
  Float_t NTrigger[nummolt+1]={0};
  Float_t NTriggerAll=0;

  Float_t YSup=0.01;
  Float_t YInf=-0.001;
  if (type==8){
    YSup=0.0006;
    YInf = -0.0001;
  }
  //I define same variables (common to all files) and get some default info
  TString nameME[nummolt+1][numPtV0];
  TString PathIn0;
  TString PathIn0Def;
  PathIn0 = Dir+"/histo/AngularCorrelation"  +file+Path1 ;
  if (PtBinning>0)    PathIn0 +=Form("_PtBinning%i",PtBinning);
  //  if (sysTrigger==0)
  if (isForPreliminary)  PathIn0 +=Form("_Jet%.2f",0.75); //just because of a difference in input file name
  PathIn0Def=PathIn0;//+Form("_Jet%.2f",0.75);
  if(type>=0){
    if (!ishhCorr)    {  PathIn0 +="_"+tipo[type]; PathIn0Def +="_"+tipo[type];}
    PathIn0 +=Srap[israp];
    PathIn0 +=SSkipAssoc[SkipAssoc];
    PathIn0Def +=Srap[israp];
    PathIn0Def +=SSkipAssoc[SkipAssoc];
  }
  PathInDef= PathIn0Def + Form("_SysT%i_SysV0Default_Sys%i_PtMin%.1f_Output.root", 0,  sys, PtTrigMin);

  TString PathInOOJNew;
  TString PathInOOJ= "OOJComparison"+year;
  if (PtBinning>0) PathInOOJ+=Form("_PtBinning%i",PtBinning);
  PathInOOJ+= Path1;
  //  PathInOOJ+="_Jet0.75";
  if(type>=0){
    PathInOOJ +="_"+tipo[type];
    PathInOOJ +=Srap[israp];
    PathInOOJ +=SSkipAssoc[SkipAssoc];
  }

  TString PathInDefOOJNew=PathInOOJ;
  PathInDefOOJNew+=Form("_SysV0Default_PtTrigMin%.1f_PtTrigMin%.1f_Output.root", PtTrigMin, PtTrigMin2); 

  TFile* finDef = new TFile(PathInDef, "");
  if (!finDef) {cout << " default file not present " << endl; return;}
  TFile* finDefOOJNew;
  if (type==8){
    finDefOOJNew = new TFile(PathInDefOOJNew, "");
    if (!finDefOOJNew) return;
  }

  for(Int_t m=0; m<nummolt+1; m++){
    cout << "m " << m << endl;
    if (m!=nummolt) {
      NTrigger[m] =    NTriggerHisto->GetBinContent(NTriggerHisto->FindBin(Nmolt[m]+0.0001));
      NTriggerAll +=    NTriggerHisto->GetBinContent(NTriggerHisto->FindBin(Nmolt[m]+0.0001));
    }
    else NTrigger[m] = NTriggerAll;

    canvasPlotProj[m] = new TCanvas (Form("canvasPlot_m%i", m), Form("canvasPlot_m%i", m), 1300, 800);
    canvasPlotProj[m]-> Divide((float)numPtV0/2+1, 2);
    
    canvasPlotProjAllSist[m] = new TCanvas (Form("canvasPlotProjAllSist_m%i", m), Form("canvasPlotProjAllSist_m%i", m), 1300, 800);
    canvasPlotProjAllSist[m]-> Divide((float)numPtV0/2+1, 2);

    canvasRelError[m] = new TCanvas (Form("canvasRelError_m%i", m), Form("canvasRelError_m%i", m), 1300, 800);
    canvasRelError[m]-> Divide((float)numPtV0/2+1, 2);

    canvasPlotRatio[m] = new TCanvas (Form("canvasPlotRatio_m%i", m), Form("canvasPlotRatio_m%i", m), 1300, 800);
    canvasPlotRatio[m]-> Divide((float)numPtV0/2+1, 2);

    for(Int_t v=PtV0Min; v<numPtV0Max; v++){
      //      cout << " v " << v << endl;     
      for (Int_t i=0; i<3; i++){	
	canvasSys[m][v][i] = new TCanvas (Form("canvasSys_m%i_v%i_i%i", m, v,i), Form("canvasSys_m%i_v%i_i%i", m,v,i), 1300, 800);
	canvasSys[m][v][i]-> Divide(5,2);
      }
      nameME[m][v]="ME_";
      nameME[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v];
      for (Int_t iregion=0; iregion<numRegions; iregion++){
	//	cout << "iregion " << iregion << endl;
	if (type==8 && (NPtV0[v] <2.5-0.001) && iregion==0){
	  hDeltaPhiDef[m][v][iregion]=(TH1F*)finDefOOJNew->Get(nameME[m][v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmoothBis");
	  cout << " hey there! I'm getting jet distr obtained with new OOJ subtraction method" << endl;
	}
	else	hDeltaPhiDef[m][v][iregion]=(TH1F*)finDef->Get(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr");
	if (iregion>0) hDeltaPhiDef[m][v][iregion]->Rebin(2);
	if (!(type==8 && (NPtV0[v] <2.5-0.001) && iregion==0)) hDeltaPhiDef[m][v][iregion]->Scale(1./NTrigger[m]);
	hDeltaPhiDefClone[m][v][iregion]=(TH1F*)	hDeltaPhiDef[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr_Clone");
	hDeltaPhiMean[m][v][iregion]=(TH1F*)	hDeltaPhiDef[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr_Mean");
	hDeltaPhiDefSys[m][v][iregion]=(TH1F*)	hDeltaPhiDef[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr_DefSys");
	hDeltaPhiRatio[m][v][iregion]=(TH1F*)	hDeltaPhiDef[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr_Ratio");
	hDeltaPhiRelStat[m][v][iregion]=(TH1F*)	hDeltaPhiDef[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr_RelStat");
	hDeltaPhiRelSist[m][v][iregion]=(TH1F*)	hDeltaPhiDef[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr_RelSist");
	if (m==0 && v==PtV0Min)	cout << "number dphi bins " << hDeltaPhiDef[m][v][iregion]->GetNbinsX()<< endl;
	for (Int_t dphi =0; dphi<hDeltaPhiDef[m][v][iregion]->GetNbinsX(); dphi++){
	  //	  cout << dphi << endl;
	  hSpread[m][v][dphi][iregion]= new TH1F ("hSpread"+ SRegion[iregion] + "_m" + Smolt[m] + "_v" + SPtV0[v]+ Form("DPhi%i", dphi),  "hSpread"+ SRegion[iregion] + "_m" + Smolt[m] + "_v" + SPtV0[v]+ Form("DPhi%.2f%-.2f", hDeltaPhiDef[m][v][iregion]->GetXaxis()->GetBinLowEdge(dphi+1), hDeltaPhiDef[m][v][iregion]->GetXaxis()->GetBinUpEdge(dphi+1)), 400, -10, 10);
	  hYieldDef[m][v][dphi][iregion] = 	hDeltaPhiDef[m][v][iregion]->GetBinContent(dphi+1);
	  hYieldEStatDef[m][v][dphi][iregion] = 	hDeltaPhiDef[m][v][iregion]->GetBinError(dphi+1);
	}
	//	cout << "end of dphi loop " << endl;
      }//iregion
    }//v
  }
  
  cout << " I got default histogram from file " << PathInDef << endl;

  TH1F *     histoTopoSel;
  TH2F *  histoCosCtau = new TH2F("histoCosCtau", "histoCosCtau", 100, 0.9, 1,100, 0, 30);
  for (Int_t sysV0index =0; sysV0index < numsysV0index; sysV0index++){
    if (sysV0index==57 && type==0 && sysTrigger==1) continue;
    if (sysTrigger==0)     PathInFirstBis =PathInFirst + Form("_MassDistr_SysT%i_SysV0index%i_PtMin%.1f.root",sysTrigger, sysV0index,  PtTrigMin);
    else   PathInFirstBis =PathInFirst + Form("_MassDistr_SysTindex%i_SysV0%i_PtMin%.1f.root",sysV0index, 0, PtTrigMin);
    fileinFirst = new TFile(PathInFirstBis, "");
    //    cout << sysV0index << endl;
    if (!fileinFirst) continue;
    histoTopoSel = (TH1F*) fileinFirst->Get("histoTopoSel");
    histoCosCtau ->Fill(histoTopoSel->GetBinContent(1), histoTopoSel->GetBinContent(5));
  }
  
  fileout->WriteTObject(histoCosCtau );

  cout << "\n\n second part of the macro " << endl;
  Int_t	NumberFilesUsed=0;
  for(Int_t m=nummolt; m>=0; m--){
    cout << "\n\n\n" << m << endl;
  //  for(Int_t m=0; m>=0; m--){
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){
      NumberFilesUsed=0;
      for (Int_t sysV0index =0; sysV0index < numsysV0index; sysV0index++){
	if (sysV0index==57 && type==0 && sysTrigger==1) continue;
	//	cout << "sysV0index" << sysV0index << endl;
	PathInOOJNew=PathInOOJ;
	if (sysTrigger==0) PathInOOJNew+=Form("_SysV0index%i", sysV0index);
	else  PathInOOJNew+=Form("_SysTindex%i", sysV0index);
	PathInOOJNew+=Form("_PtTrigMin%.1f_PtTrigMin%.1f_Output.root", PtTrigMin, PtTrigMin2); 	

	PathIn = PathIn0;
	if (sysTrigger==0)	PathIn+=   Form("_SysT%i_SysV0index%i_Sys%i_PtMin%.1f_Output.root", sysTrigger, sysV0index, sys, PtTrigMin);
	else 	PathIn+=   Form("_SysTindex%i_SysV0%i_Sys%i_PtMin%.1f_Output.root", sysV0index, 0, sys, PtTrigMin);
	filein = new TFile(PathIn, "");
	if (!filein) continue;
	//	cout << "I've found the file" << endl;

	if (type==8){
	  fileinOOJNew = new TFile(PathInOOJNew, "");
	  if (!fileinOOJNew) continue;
	}

	for (Int_t iregion=0; iregion<numRegions; iregion++){
	  if (type==8 && (NPtV0[v] <2.5-0.001) && iregion==0){
	    hDeltaPhi[m][v][iregion]=(TH1F*)fileinOOJNew->Get(nameME[m][v]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmoothBis");
	    //	    cout << " hey there! I'm using new OOJ distr!" << endl;
	  }
	  else{
	  hDeltaPhi[m][v][iregion]=(TH1F*)filein->Get(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr");
	  }

	  if (!hDeltaPhi[m][v][iregion]) {
	    cout << PathIn << endl;
	    cout << "m " << m << " v " << v << " index " << sysV0index << endl;
	    cout << "histogram not present for m " << m << " v " << v << " and " << SRegion[iregion] <<  endl; 
	    cout << "the name of the histo you were looking for is " << nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+"EffCorr" << endl;
	    continue;
	  }
	  NumberFilesUsed++;
	  if (iregion>0) hDeltaPhi[m][v][iregion]->Rebin(2);
	  if (!(type==8 && (NPtV0[v] <2.5-0.001) && iregion==0)) hDeltaPhi[m][v][iregion]->Scale(1./NTrigger[m]);

	  if (iregion==2 && sysV0index<10){
	    canvasPlotProjAllSist[m]->cd(v+1);
	 
	    hDeltaPhiSist[m][v][iregion][sysV0index]=(TH1F*)	  hDeltaPhi[m][v][iregion]->Clone(nameME[m][v]+"_AC_phi_V0Sub_"+SRegionInput[iregion]+Form("EffCorr_%i",sysV0index));
	    hDeltaPhiSist[m][v][iregion][sysV0index]->SetLineColor(kBlack);
	    if (!(type==8 && (NPtV0[v] <2.5-0.001) && iregion==0))	    hDeltaPhiSist[m][v][iregion][sysV0index]->Scale(1./NTrigger[m]);
	    hDeltaPhiSist[m][v][iregion][sysV0index]->SetTitle("v "+SPtV0[v]);
	    if (sysV0index==0)	    hDeltaPhiSist[m][v][iregion][sysV0index]->DrawCopy("");
	    hDeltaPhiSist[m][v][iregion][sysV0index]->DrawCopy("same");

	    if (sysV0index==0) {	    
	      hDeltaPhiDefClone[m][v][iregion]->SetLineColor(kRed);
	      hDeltaPhiDefClone[m][v][iregion]->SetMarkerColor(kRed);
	      hDeltaPhiDefClone[m][v][iregion]->SetMarkerStyle(1);
	      hDeltaPhiDefClone[m][v][iregion]->SetTitle("v "+SPtV0[v]);
	      hDeltaPhiDefClone[m][v][iregion]->Draw("same");
	    }
	    //	    cout << " I've drawn the histo" << 	    hDeltaPhiSist[m][v][iregion][sysV0index]->GetBinContent(1) << " vs def: " << 	    hDeltaPhiDefClone[m][v][iregion]->GetBinContent(1) << endl;
	  }
	  //	  cout << "number dphi bins " << hDeltaPhi[m][v][iregion]->GetNbinsX()<< endl;

	  for (Int_t dphi =0; dphi<hDeltaPhi[m][v][iregion]->GetNbinsX(); dphi++){
	    hYield[m][v][dphi][iregion] = 	hDeltaPhi[m][v][iregion]->GetBinContent(dphi+1);

	    /*
	    if (v==7 && m==4 && iregion==1 && (dphi==1 || dphi==12 || dphi==2)){
	      cout << " dphi " << dphi << " yield def " << hYieldDef[m][v][dphi][iregion] << " yield " << hYield[m][v][dphi][iregion]<< endl;
	    }
	    if ( hYield[m][v][dphi][iregion]==0) {
	      cout << " dphi " << dphi << " v " << v << " iregion " << iregion << " yield def " << hYieldDef[m][v][dphi][iregion] << " yield " << hYield[m][v][dphi][iregion]<< endl;
	    }
	   
	    if (TMath::Abs((hYield[m][v][dphi][iregion] - hYieldDef[m][v][dphi][iregion])/hYieldEStatDef[m][v][dphi][iregion])>3){
	      cout << " v " << v << "dphi " << dphi << " " <<	    hYield[m][v][dphi][iregion]<< " " << 	    hYieldDef[m][v][dphi][iregion]<< " " << (hYield[m][v][dphi][iregion] - hYieldDef[m][v][dphi][iregion])/hYieldEStatDef[m][v][dphi][iregion] << endl;
	    }
	    */
	    hSpread[m][v][dphi][iregion]->Fill((hYield[m][v][dphi][iregion] - hYieldDef[m][v][dphi][iregion])/hYieldEStatDef[m][v][dphi][iregion]);
	    //	    if (iregion==2 && dphi==10) cout   << hYield[m][v][dphi][iregion]<<endl;
	    //	    hSpread[m][v][dphi][iregion]->Fill(hYield[m][v][dphi][iregion]);
	  } //phi
	}//region
	filein ->Close();
	if (type==8) 	  fileinOOJNew ->Close();
      } //sysV0index
      //      cout << "\n\n I got all info from files " << endl;
      
     
      for (Int_t iregion=0; iregion<numRegions; iregion++){
	if (hDeltaPhi[m][v][iregion]){
	  //	  cout << " m " << m << " v " << v << " iregion " << iregion << endl;
	  for (Int_t dphi =0; dphi<hDeltaPhiDef[m][v][iregion]->GetNbinsX(); dphi++){
	    hYieldMean[m][v][dphi][iregion] = hSpread[m][v][dphi][iregion]->GetMean()*hYieldEStatDef[m][v][dphi][iregion] + hYieldDef[m][v][dphi][iregion];
	    hYieldSpread[m][v][dphi][iregion] = hSpread[m][v][dphi][iregion]->GetRMS()*hYieldEStatDef[m][v][dphi][iregion];
	    //hYieldMean[m][v][dphi][iregion] = hSpread[m][v][dphi][iregion]->GetMean();
	    hSpread[m][v][dphi][iregion]->SetLineColor(ColorRegion[iregion]);
	    if (iregion==0)   hSpread[m][v][dphi][iregion]->Rebin(2);
	    else    hSpread[m][v][dphi][iregion]->Rebin(8);
	    hSpread[m][v][dphi][iregion]->GetXaxis()->SetRangeUser(-10 ,10);//was -3, 3
	    hSpread[m][v][dphi][iregion]->GetXaxis()->SetTitle("(Yield_{i} - Yield_{Def})/#sigma^{stat}_{YieldDef}");
	    hSpread[m][v][dphi][iregion]   ->SetTitle(Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	    /*
	    if (hSpread[m][v][dphi][iregion]->GetRMS()>1){
	      cout << "large spread " << " v " << v << " dphi " << dphi << " " << hSpread[m][v][dphi][iregion]->GetRMS()<< endl;
	    }
	    if (TMath::Abs(hSpread[m][v][dphi][iregion]->GetMean())>1){
	      cout <<"large mean " <<  " v " << v << " dphi " << dphi << " " << hSpread[m][v][dphi][iregion]->GetMean()<< endl;
	    }
	    */
	    if (dphi<10){
	      canvasSys[m][v][0]->cd(dphi+1);
	      hSpread[m][v][dphi][iregion]->Draw("same");
	    }
	    else   if (dphi<20){
	      canvasSys[m][v][1]->cd(dphi+1-10);
	      hSpread[m][v][dphi][iregion]->Draw("same");
	    }
	    else   if (dphi<30){
	      canvasSys[m][v][2]->cd(dphi+1-20);
	      hSpread[m][v][dphi][iregion]->Draw("same");
	    }
	    
	    if (iregion==2){
	      //	    cout << "dphi " << endl;
	      //	    cout <<hSpread[m][v][dphi][iregion]->GetMean()<<" " << hYieldMean[m][v][dphi][iregion] << " def " << hYieldDef[m][v][dphi][iregion]<< " +- " << hYieldEStatDef[m][v][dphi][iregion]<< endl;
	    }

	    //hYieldSpread[m][v][dphi][iregion] = hSpread[m][v][dphi][iregion]->GetRMS();
	    hDeltaPhiMean[m][v][iregion]->SetBinContent(dphi+1, 	  hYieldMean[m][v][dphi][iregion]);
	    hDeltaPhiMean[m][v][iregion]->SetBinError(dphi+1,	hDeltaPhiDef[m][v][iregion]->GetBinError(dphi+1) );

	    if (hDeltaPhiDef[m][v][iregion]->GetBinContent(dphi+1) != 0) {
	    hDeltaPhiDefSys[m][v][iregion]->SetBinContent(dphi+1, 	  hYieldDef[m][v][dphi][iregion]);
	    hDeltaPhiDefSys[m][v][iregion]->SetBinError(dphi+1, 	  hYieldSpread[m][v][dphi][iregion]/hDeltaPhiMean[m][v][iregion]->GetBinContent(dphi+1)*hDeltaPhiDef[m][v][iregion]->GetBinContent(dphi+1));
	    }
	    else{
	      hDeltaPhiDefSys[m][v][iregion]->SetBinError(dphi+1, 0);
	    }
	    //	    hDeltaPhiRatio[m][v][iregion]->SetBinContent(dphi+1,hSpread[m][v][dphi][iregion]->GetMean()*hYieldEStatDef[m][v][dphi][iregion]/hSpread[m][v][dphi][iregion]->GetRMS());
	    hDeltaPhiRatio[m][v][iregion]->SetBinContent(dphi+1, (hYieldMean[m][v][dphi][iregion]-hYieldDef[m][v][dphi][iregion])/hYieldSpread[m][v][dphi][iregion]);
	    hDeltaPhiRatio[m][v][iregion]->SetBinError(dphi+1,0);

	    hDeltaPhiRelStat[m][v][iregion]->SetBinContent(dphi+1,hDeltaPhiDef[m][v][iregion]->GetBinError(dphi+1) / hDeltaPhiDef[m][v][iregion]->GetBinContent(dphi+1) );
	    hDeltaPhiRelStat[m][v][iregion]->SetBinError(dphi+1,0);

	    hDeltaPhiRelSist[m][v][iregion]->SetBinContent(dphi+1,hDeltaPhiDefSys[m][v][iregion]->GetBinError(dphi+1) / hDeltaPhiDef[m][v][iregion]->GetBinContent(dphi+1) );
	    hDeltaPhiRelSist[m][v][iregion]->SetBinError(dphi+1,0);

	    /*
	    if (v==7 && m==4 && iregion==1){
	      cout << " dphi " << dphi << " mean " << hSpread[m][v][dphi][iregion]->GetMean() << " rms " << hSpread[m][v][dphi][iregion]->GetRMS()<< endl;
	      cout <<   hDeltaPhiDefSys[m][v][iregion]->GetBinContent(dphi+1) << " " <<   hDeltaPhiDefSys[m][v][iregion]->GetBinError(dphi+1)<< endl;
	    }

	    if (	    hDeltaPhiRelSist[m][v][iregion]->GetBinContent(dphi+1) > 0.5){
	      cout <<"m " << m << " v " << v << " iregion " << iregion <<  " rel sist " << 	    hDeltaPhiRelSist[m][v][iregion]->GetBinContent(dphi+1)<< endl;
	    }
	    if (	    hDeltaPhiRelSist[m][v][iregion]->GetBinContent(dphi+1) ==0){
	      cout <<" zero m " << m << " v " << v << " iregion " << iregion <<  " rel sist " << 	    hDeltaPhiRelSist[m][v][iregion]->GetBinContent(dphi+1)<< endl;
	    }
	    */
	  }//phi

	  /*
	hDeltaPhiDefSys[m][v][iregion]->Scale(1./NTrigger[m]);
	hDeltaPhiDef[m][v][iregion]   ->Scale(1./NTrigger[m]);
	hDeltaPhiMean[m][v][iregion]  ->Scale(1./NTrigger[m]);
	  */
	if (m==0 && type==8){
	  YSup=0.001;
	  YInf = -0.0001;
	}

	hDeltaPhiDefSys[m][v][iregion]->GetYaxis()->SetRangeUser(YInf,YSup);
	hDeltaPhiDef[m][v][iregion]   ->GetYaxis()->SetRangeUser(YInf,YSup);
	hDeltaPhiMean[m][v][iregion]  ->GetYaxis()->SetRangeUser(YInf,YSup);
	hDeltaPhiRelStat[m][v][iregion]  ->GetYaxis()->SetRangeUser(0, 0.1);
	hDeltaPhiRelSist[m][v][iregion]  ->GetYaxis()->SetRangeUser(0, 0.1);
	if (type==8){
	hDeltaPhiRelStat[m][v][iregion]  ->GetYaxis()->SetRangeUser(0, 0.4);
	hDeltaPhiRelSist[m][v][iregion]  ->GetYaxis()->SetRangeUser(0, 0.4);
	}
	hDeltaPhiRatio[m][v][iregion]  ->GetYaxis()->SetRangeUser(-3,3);

	hDeltaPhiDefSys[m][v][iregion]->SetLineColor(ColorRegion[iregion]);
	hDeltaPhiDef[m][v][iregion]   ->SetLineColor(ColorRegion[iregion]);
	hDeltaPhiMean[m][v][iregion]  ->SetLineColor(ColorRegionBis[iregion]);
	hDeltaPhiRatio[m][v][iregion]  ->SetLineColor(ColorRegionBis[iregion]);
	hDeltaPhiRelStat[m][v][iregion]   ->SetLineColor(ColorRegion[iregion]);
	hDeltaPhiRelSist[m][v][iregion]   ->SetLineColor(ColorRegion[iregion]);

	hDeltaPhiDefSys[m][v][iregion]->SetMarkerColor(ColorRegion[iregion]);
	hDeltaPhiDef[m][v][iregion]   ->SetMarkerColor(ColorRegion[iregion]);
	hDeltaPhiMean[m][v][iregion]  ->SetMarkerColor(ColorRegionBis[iregion]);
	hDeltaPhiRatio[m][v][iregion]  ->SetMarkerColor(ColorRegionBis[iregion]);
	hDeltaPhiRelStat[m][v][iregion]   ->SetMarkerColor(ColorRegion[iregion]);
	hDeltaPhiRelSist[m][v][iregion]   ->SetMarkerColor(ColorRegion[iregion]);
	hDeltaPhiRelStat[m][v][iregion]   ->SetMarkerStyle(33);
	hDeltaPhiRelSist[m][v][iregion]   ->SetMarkerStyle(27);
	hDeltaPhiRelStat[m][v][iregion]   ->GetYaxis()->SetTitle("Relative uncertainty");
	hDeltaPhiRelSist[m][v][iregion]   ->GetYaxis()->SetTitle("Relative uncertainty");
	hDeltaPhiRelStat[m][v][iregion]   ->SetTitle(Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	hDeltaPhiRelSist[m][v][iregion]   ->SetTitle(Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));

	hDeltaPhiMean[m][v][iregion]  ->SetTitle("v "+SPtV0[v]);
	hDeltaPhiRatio[m][v][iregion]  ->SetTitle("v "+SPtV0[v]);

	if (iregion>0){
	  //	  hDeltaPhiDefSys[m][v][iregion]->Rebin(2);
	  //	  hDeltaPhiDef[m][v][iregion]   ->Rebin(2);
	  //	  hDeltaPhiMean[m][v][iregion]  ->Rebin(2);
	}

	canvasPlotProj[m]->cd(v+1);
	gPad->SetLeftMargin(0.15);
	hDeltaPhiDef[m][v][iregion]   ->SetTitle(Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	hDeltaPhiDefSys[m][v][iregion]   ->SetTitle(Form("%.1f < p_{T} < %.1f", NPtV0[v], NPtV0[v+1]));
	hDeltaPhiDef[m][v][iregion]->SetMarkerStyle(1);
	hDeltaPhiDef[m][v][iregion]->DrawClone("same e");
	//	hDeltaPhiMean[m][v][iregion]->DrawClone("same e");
	hDeltaPhiDefSys[m][v][iregion]->SetMarkerStyle(1);
	hDeltaPhiDefSys[m][v][iregion]->SetFillStyle(0);
	hDeltaPhiDefSys[m][v][iregion]->DrawClone("same e2");
	lineat1->Draw("same");

	//	hDeltaPhiDefSys[m][v][iregion]->Scale(NTrigger[m]);
	//	hDeltaPhiDef[m][v][iregion]   ->Scale(NTrigger[m]);
	//	hDeltaPhiMean[m][v][iregion]  ->Scale(NTrigger[m]);	

	canvasPlotRatio[m]->cd(v+1);
	gPad->SetLeftMargin(0.15);
	hDeltaPhiRatio[m][v][iregion]->GetYaxis()->SetTitle("(Yield_{Def}-Yield_{Mean})/RMS");
	hDeltaPhiRatio[m][v][iregion]->GetYaxis()->SetTitleSize(0.05);
	hDeltaPhiRatio[m][v][iregion]->DrawClone("same");
	lineat0->Draw("same");

	canvasRelError[m]->cd(v+1);
	if (iregion!=0){
	  hDeltaPhiRelStat[m][v][iregion]->Draw("same p");
	  hDeltaPhiRelSist[m][v][iregion]->Draw("same p");
	}
	}
      }//region 
      for (Int_t i=0; i<3; i++){      
	fileout->WriteTObject(canvasSys[m][v][i]);
      }
    }//v

    TString OutputFile = "PictureForNote";
    if (!isForPreliminary) OutputFile = "NewResultsForNote";
    canvasPlotProj[m]->SaveAs(OutputFile + "DeltaPhiProjSE"+tipo[type]+"_m"+Smolt[m]);
    fileout->WriteTObject(canvasPlotProj[m]);
    fileout->WriteTObject(canvasPlotProjAllSist[m]);
    fileout->WriteTObject(canvasPlotRatio[m]);
    fileout->WriteTObject(canvasRelError[m]);

  }//m 

  for(Int_t m=nummolt; m>=0; m--){
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){
      for (Int_t iregion=0; iregion<numRegions; iregion++){
	    fileout->WriteTObject(hDeltaPhiDef[m][v][iregion]);
	    fileout->WriteTObject(hDeltaPhiMean[m][v][iregion]);
	    fileout->WriteTObject(hDeltaPhiDefSys[m][v][iregion]);
	for (Int_t dphi =0; dphi<hDeltaPhiDef[m][v][iregion]->GetNbinsX(); dphi++){
	  if (m==nummolt){
	    fileout->WriteTObject(hSpread[m][v][dphi][iregion]);
	  }
	}
      }
    } //v
  }//m


  fileout->Close();

  cout << "where is the bin content =0?" << endl;
  for(Int_t m=nummolt; m>=0; m--){
    for(Int_t v=PtV0Min; v<numPtV0Max; v++){
      for (Int_t iregion=0; iregion<numRegions; iregion++){
	for (Int_t dphi =0; dphi<hDeltaPhiDef[m][v][iregion]->GetNbinsX(); dphi++){

	  if (hDeltaPhiDef[m][v][iregion]->GetBinContent(dphi+1) == 0){
	    cout << " m "<< m << " v " << v << " iregion " << iregion << " dphi " << dphi << endl;
	  }
	}
      }
    }
  }

  cout << " number of files used " << NumberFilesUsed<< endl;
  cout << "\npartendo dai file (esempio) " << PathIn <<  " ho creato: "<< endl;
  if (type==8)   cout << "\npartendo dai file (esempio) per la distribuzione in jet " << PathInOOJNew <<  " ho creato: "<< endl;
  cout << "\nil file " << PathOut << endl;

}
