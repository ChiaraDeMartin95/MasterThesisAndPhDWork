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
void readTreePLChiarahK0s_first( Int_t sysV0=0,Int_t sysV0Index=0, Int_t type=0 /*type = 0 for K0s */,Bool_t SkipAssoc=0 ,Int_t israp=1, Bool_t ishhCorr=0, Float_t PtTrigMin=0.15, Float_t ptjmax=30, bool isMC = 1,Bool_t isEfficiency=1, Int_t sysTrigger=0, TString year= "LHC19h11_HM_INELgt0Bis286380"/*"LHC17o_HM_INELgt0281961"/*"hK0s_pp_18p_AOD264"/*"17pq_hK0s"/*"17d20bEPOS_hK0s_EtaEff"/*"2019h11b+extra_hK0s"/"2019h11_HM_hK0s"/"15g3c3_hK0s"/"2018f1_extra_15runs_Trig0Pt"/*"17pq_hK0s_pttrig0.15"/"1617_AOD234_hK0s"/*"1617_GP_AOD235_With18c12b"/*"LHC16_GP_AOD235_17d16"/*"2017e5_extra_AOD235_hK0s"/*"2019h11a_extra_HM_hK0s"/"17pq_hK0s_FAST"/"AllhK0sHM_RedNo16k"/*"2018f1_extra_30runs_Mylabel_hK0s"/*"2018f1_extra_Nature_15runs"/"2019h11c_extra_HM_hK0s"/*"2016k_HM_hK0s"/*"2018g4_extra_EtaEff_hK0s"/*"2016k_HM_hK0s"/*"2018g4_extra_EtaEff_hK0s"/*"LHC17_hK0s"/*"17d20bEPOS_hK0s"/"1617MC_hK0s"/*"2018f1_extra_Reco_hK0s"/*"2016k_hK0s"/*"2016kehjl_hK0s"/*"2016k_hK0s"/"17anch17_hK0s"/"1617_hK0s"/*"2018f1_extra_hK0s_CP"*/, TString year0="2016", TString Path1 ="", Bool_t CommonParton=0, Int_t PtBinning=1, Bool_t isSysStudy=0, Int_t numsysV0index=20, Bool_t IsPtTrigMultDep = 0, Int_t Region=0, Bool_t isEfficiencyMassSel=0, Bool_t isHM=1, Bool_t isNewInputPath=1, Int_t MultBinning=0, Bool_t ispp5TeV=0, Bool_t isINEL=0)
{
  
  //Region=0; default, no selections on  dEta and dPhi
  //Region=1; Jet selection //-0.85 < dPhi < 0.85 ; -0.75 < dEta<0.75
  //Region=2; OOJ selection //0.85 < dPhi < 2.05 ; 0.75 < dEta<1.18 

  if (year == "LHC17o_HM_INELgt0281961" || year == "LHC19h11_HM_INELgt0Bis286380") isINEL=1;
  if (isINEL){
    SkipAssoc =0;
    PtTrigMin=0.;
    ptjmax = 30;
    israp =1;
  }
  if (sysV0Index>20) return; 
  //SysStudy: if 1, values different from default ones will be used; the topological variable changed will sysV0Index > 20 ) return;be indicated in the output file, together with its value  

  //sysV0 indicated the type of topological variable which is varied. (trial iintervals are indicated below)
  //***  sysV0=1: cosinePointingAngle : 0.970 - 0.999
  //***  sysV0=2: DCAPosToPV          : 0.05 - 0.14
  //***  sysV0=3: DCANegToPV          : 0.05 - 0.14
  //***  sysV0=4: DCAV0ToPV           : 0.2 - 0.6
  //***  sysV0=5: ctau                : 10 - 40
  //***  sysV0=6: LambdaMass          : 3 - 9  

  //sysV0Index is the index of the array of the set of selections

  //  Float_t CosinePAngleExtr[2] = {0.970, 0.998};
  Float_t CosinePAngleExtr[2] = {0.995, 1};
  Float_t CosinePAngle = 0;
  if (sysV0==1) CosinePAngle = (CosinePAngleExtr[1] -  CosinePAngleExtr[0])/numsysV0index * sysV0Index +  CosinePAngleExtr[0];
  else  CosinePAngle=0.995; //0.995;
  
  Float_t DCAPosToPVExtr[2] = {0.05, 0.14};
  Float_t DCAPosToPV = 0;
  if (sysV0==2) DCAPosToPV = (DCAPosToPVExtr[1] -  DCAPosToPVExtr[0])/numsysV0index * sysV0Index +  DCAPosToPVExtr[0];
  else DCAPosToPV = 0.06;

  Float_t DCANegToPVExtr[2] = {0.05, 0.14};
  Float_t DCANegToPV = 0;
  if (sysV0==3) DCANegToPV = (DCANegToPVExtr[1] -  DCANegToPVExtr[0])/numsysV0index * sysV0Index +  DCANegToPVExtr[0];
  else DCANegToPV = 0.06;

  Float_t DCAV0ToPVExtr[2] = {0.2, 0.6};
  Float_t DCAV0ToPV = 0;
  if (sysV0==4) DCAV0ToPV = (DCAV0ToPVExtr[1] -  DCAV0ToPVExtr[0])/numsysV0index * sysV0Index +  DCAV0ToPVExtr[0];
  else DCAV0ToPV = 0.5;

  //  Float_t ctauExtr[2] = {10, 40};
  Float_t ctauExtr[2] = {2, 20};
  Float_t ctau = 0;
  if (sysV0==5) ctau = (ctauExtr[1] -  ctauExtr[0])/numsysV0index * sysV0Index +  ctauExtr[0];
  else ctau = 20;

  //  Float_t LambdaMassWindowExtr[2] = {0.003, 0.009};
  Float_t LambdaMassWindowExtr[2] = {0.000, 0.03};
  Float_t LambdaMassWindow =0; 
  if (sysV0==6)  LambdaMassWindow= (LambdaMassWindowExtr[1] -  LambdaMassWindowExtr[0])/numsysV0index * sysV0Index +  LambdaMassWindowExtr[0];
  else LambdaMassWindow = 0.005;

  if (sysV0==0 && sysV0Index==1 && isSysStudy){ //I set the tighest selections to check effect on signal yield (the tighest selections are those which, when applied together with default selections for all other variables, give a -5% in signal yield)
    LambdaMassWindow = 0.017;
    ctau = 7.5;
    DCAV0ToPV = 0.3;
    DCAPosToPV = 0.14;
    DCANegToPV = 0.14;
    CosinePAngle = 0.99875;
  }

  if (sysV0==0 && sysV0Index==2 && isSysStudy){ //I set the loosest selections to check effect on signal yield (the loosest selections are those applied in task, for most of te variables)
    LambdaMassWindow = 0;
    ctau = 30;
    DCAV0ToPV = 0.6;
    DCAPosToPV = 0.03;
    DCANegToPV = 0.03;
    CosinePAngle = 0.97;
  }

  if (sysV0==0 && sysV0Index==3 && isSysStudy){ //I set tighter selections (less tight than those for sysV0Index ==1) to check effect on signal yield (the tighest selections are those which, when applied together with default selections for all other variables, give a -2% in signal yield)
    LambdaMassWindow = 0.010;
    ctau = 10.5;
    DCAV0ToPV = 0.4;
    DCAPosToPV = 0.09;
    DCANegToPV = 0.09;
    CosinePAngle = 0.998;
  }

  //rap=0 : no rapidity window chsen for cascade, |Eta| < 0.8; rap=1 |y| < 0.5

  if (ishhCorr && !isEfficiency) {
    cout << "This macro should not be run is hh correlation is studied; go directly to readTreePLChiara_second " << endl;
    cout << "for hh correlation it is necessary only to compute efficiency of associated particles " << endl;
    return;
  }

  if (israp>1) return;
  if (sysV0>6) return;
  if (sysV0>2 && ishhCorr) return;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=4;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda", "LambdaAntiLambda"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};

 TString PathInMass= "FinalOutput/DATA" + year0 + "/invmass_distribution_thesis/invmass_dist\
ribution";
 TString PathInMassDef;

  TString PathIn="./FinalOutput/AnalysisResults";
  TString PathOut="./FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";

  PathIn+=year;
  PathOut+=year;  
  if (ishhCorr)  PathIn+="_hhCorr";

  if(isMC && isEfficiency){
    PathIn+="_MCEff";
    PathOut+="_MCEff";
    PathInMass+="_MCEff";
  }
  if(isMC && !isEfficiency){
    PathIn+="_MCTruth";
    PathOut+="_MCTruth";
  }
 
  //  PathIn+=Path1;
  PathIn+=".root";
  PathOut+=Path1;

  if(PtBinning)  PathInMass+=Form("_PtBinning%i",PtBinning);
  PathInMass+=Path1;
  if (PtBinning>0)  PathOut+=Form("_PtBinning%i",PtBinning);
  PathOut+="_";
  if (!ishhCorr){
    PathOut +=tipo[type];
  }
  PathOut +=Srap[israp];
  if (!SkipAssoc)  PathOut +="_AllAssoc";
  if (ishhCorr) PathOut +="_hhCorr";  
  if (isSysStudy)  PathOut +=Form("_MassDistr_SysT%i_SysV0%i_index%i_PtMin%.1f",sysTrigger, sysV0, sysV0Index, PtTrigMin); 
  else   PathOut +=Form("_MassDistr_SysT%i_SysV0%i_PtMin%.1f",sysTrigger, sysV0, PtTrigMin); 
  if (IsPtTrigMultDep) PathOut += "_IsPtTrigMultDep";
  if(Region==2)  PathOut+="_dEtadPhiOOJRegion"; 
  else  if (Region==1)  PathOut+="_dEtadPhiJetRegion";
  if (isEfficiencyMassSel) PathOut+= "_isEfficiencyMassSel";
  if (MultBinning!=0) PathOut+= Form("_MultBinning%i", MultBinning);
  if (isINEL) PathOut+= "_INEL";
  //  PathOut +="_HM";
  //  PathOut+="_LooseCosinePAngle";
  //  PathOut += "_NoZeroMult";
  PathOut+= ".root";


  TH1F* histoSigma;
  TH1F* histoMean;
  Double_t sigma[numtipo][nummolt+1][numPtV0]={0};
  Double_t mass[numtipo][nummolt+1][numPtV0]={0};
  TFile *fileMassSigma;
   Int_t PtBinMin=0; 
  if (!ishhCorr && type!=0) PtBinMin=1; //for associated particles different from hadrons I do not start from 0
  cout << "PtBinMin " << PtBinMin << endl;

  if (isEfficiencyMassSel){
    //  PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab];
    PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Srap[israp];
    if (!SkipAssoc)  PathInMassDef +="_AllAssoc";
    PathInMassDef+="_isMeanFixedPDG_BkgRetta";
    if(IsPtTrigMultDep) PathInMassDef  +="_IsPtTrigMultDep" ;
    PathInMassDef+= Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", 5, sysTrigger, sysV0, 0,PtTrigMin);
    fileMassSigma= new TFile(PathInMassDef);
    cout << PathInMassDef << endl;
    histoSigma=(TH1F*)fileMassSigma->Get("histo_sigma");
    histoMean=(TH1F*)fileMassSigma->Get("histo_mean");
    for (Int_t m=0; m< nummolt+1; m++){
      for(Int_t v=PtBinMin; v<numPtV0; v++){
	mass[type][m][v]=histoMean->GetBinContent(v+1);
	sigma[type][m][v]=histoSigma->GetBinContent(v+1);
	//	cout << 	sigma[type][m][v]<< endl;
      }
    }
    fileMassSigma->Close();
  }
 

  Bool_t MassInPeakK0s =0; //variable for a rough cut on K0s candidates' invariant mass (only done to fill histograms needed for deltaphi projections of K0s with or without an ancestor parton in common with trigger particle)
  TString dirinputtype[4] = {"", "Lambda", "Lambda", "Lambda"};
  TFile *fin = new TFile(PathIn);
  if (!fin) {cout << "file input not available " ; return;}
  TDirectoryFile *d;
  if (isNewInputPath){
    if (year.Index("Hybrid")!=-1)   {
      d = (TDirectoryFile*)fin->Get("MyTask_MCHybrid_PtTrigMin3.0_PtTrigMax30.0");
    }
    else if (isMC) {
      d = (TDirectoryFile*)fin->Get("MyTask_MCTruth_PtTrigMin3.0_PtTrigMax15.0");
    }
    else if (year == "2016k_HM_hK0s" ||  year.Index("2019h11")!=-1 ){
      d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin3.0_PtTrigMax30.0");
    }
    else if (year.Index("17pq_hK0s")!=-1 && isMC) {
      cout << " here I am " << endl;
      d = (TDirectoryFile*)fin->Get("MyTask_MCTruth_PtTrigMin3.0_PtTrigMax15.0");
    }
    else if (year.Index("17pq_hK0s_pttrig0.15")!=-1 && !isMC) {
      cout << " here I am " << endl;
      d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin0.2_PtTrigMax15.0");
    }
    else if (PtTrigMin == 0.15){
      d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin0.2_PtTrigMax15.0");
    }
    else if (PtTrigMin == 0){
      d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin0.0_PtTrigMax15.0");
    }
    else if (year.Index("AOD235")!=-1 && isMC) {
      d = (TDirectoryFile*)fin->Get("MyTask_MCTruth_PtTrigMin3.0_PtTrigMax15.0");
    }
    else  d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin3.0_PtTrigMax15.0");
  }
  else {
    d = (TDirectoryFile*)fin->Get("MyTask"+dirinputtype[type]);
  }
  if (year == "LHC17o_HM_INELgt0281961" || year == "LHC19h11_HM_INELgt0Bis286380") {
    d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin0.2_PtTrigMax30.0");
  }
  if (!d)  {cout << "dir input not available " ; return;}

  cout << "hey there ! "<< endl;

  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Double_t ctauK0s = 2.6844;

  Double_t     fSignTreeVariablePtTrigger;		       
  Int_t        fSignTreeVariableChargeTrigger;		       
  Double_t     fSignTreeVariableEtaTrigger;		       
  Double_t     fSignTreeVariablePhiTrigger;		       
  Double_t     fSignTreeVariableDCAz;			       
  Double_t     fSignTreeVariableDCAxy;			  
  Int_t        fSignTreeVariableisPrimaryTrigger;			   
  Int_t        fSignTreeVariablePDGCodeTrigger;                        
  Int_t        fSignTreeVariableChargeAssoc;		       
  Int_t        fSignTreeVariableisPrimaryV0;			  
  Int_t        fSignTreeVariablePDGCodeAssoc;                        
  Double_t     fSignTreeVariableRapK0Short;		       
  Bool_t        fSignTreeVariableSkipAssoc;
  Bool_t       fSignTreeVariableIsCommonParton;
  Int_t        fSignTreeVariablePdgCommonParton;                        
  Double_t     fSignTreeVariablePtV0;			       
  Double_t     fSignTreeVariableEtaV0;			       
  Double_t     fSignTreeVariablePhiV0;			       
  Double_t     fSignTreeVariableAssocDCAz;
  Double_t     fSignTreeVariableAssocDCAxy;
  Double_t     fSignTreeVariableDcaV0ToPrimVertex;
  Double_t     fSignTreeVariableDcaPosToPrimVertex;
  Double_t     fSignTreeVariableDcaNegToPrimVertex;
  Double_t     fSignTreeVariableV0CosineOfPointingAngle;
  Double_t     fSignTreeVariablectau;			       
  Double_t     fSignTreeVariableInvMassK0s;
  Double_t     fSignTreeVariableInvMassLambda;
  Double_t     fSignTreeVariableInvMassAntiLambda;
  Double_t     fSignTreeVariablePtArmenteros;
  Double_t     fSignTreeVariableAlpha;
  Double_t     fSignTreeVariableDeltaEta;		       
  Double_t     fSignTreeVariableDeltaPhi;		       
  Double_t     fSignTreeVariableDeltaTheta;		       
  Double_t     fSignTreeVariableMultiplicity;		       
  Double_t     fSignTreeVariableZvertex;                        


  Double_t     fBkgTreeVariablePtTrigger;		       
  Int_t        fBkgTreeVariableChargeTrigger;		       
  Double_t     fBkgTreeVariableEtaTrigger;		       
  Double_t     fBkgTreeVariablePhiTrigger;		       
  Double_t     fBkgTreeVariableDCAz;			       
  Double_t     fBkgTreeVariableDCAxy;			  
  Int_t        fBkgTreeVariableisPrimaryTrigger;			   
  Int_t        fBkgTreeVariablePDGCodeTrigger;                        
  Int_t        fBkgTreeVariableChargeAssoc;		       
  Int_t        fBkgTreeVariableisPrimaryV0;			  
  Int_t        fBkgTreeVariablePDGCodeAssoc;                        
  Double_t     fBkgTreeVariableRapK0Short;		       
  Bool_t        fBkgTreeVariableSkipAssoc;                        
  Bool_t       fBkgTreeVariableIsCommonParton;
  Int_t        fBkgTreeVariablePdgCommonParton;
  Double_t     fBkgTreeVariablePtV0;			       
  Double_t     fBkgTreeVariableEtaV0;			       
  Double_t     fBkgTreeVariablePhiV0;			       
  Double_t     fBkgTreeVariableAssocDCAz;
  Double_t     fBkgTreeVariableAssocDCAxy;
  Double_t     fBkgTreeVariableDcaV0ToPrimVertex;
  Double_t     fBkgTreeVariableDcaPosToPrimVertex;
  Double_t     fBkgTreeVariableDcaNegToPrimVertex;
  Double_t     fBkgTreeVariableV0CosineOfPointingAngle;
  Double_t     fBkgTreeVariablectau;			       
  Double_t     fBkgTreeVariableInvMassK0s;
  Double_t     fBkgTreeVariableInvMassLambda;
  Double_t     fBkgTreeVariableInvMassAntiLambda;
  Double_t     fBkgTreeVariablePtArmenteros;
  Double_t     fBkgTreeVariableAlpha;
  Double_t     fBkgTreeVariableDeltaEta;		       
  Double_t     fBkgTreeVariableDeltaPhi;		       
  Double_t     fBkgTreeVariableDeltaTheta;		       
  Double_t     fBkgTreeVariableMultiplicity;		       
  Double_t     fBkgTreeVariableZvertex;                        

  Int_t CounterSignPairsAfterPtMinCut=0; 
  Int_t CounterBkgPairsAfterPtMinCut=0; 
  Int_t TrueCounterSignPairsAfterPtMinCut=0; 
  Int_t TrueCounterBkgPairsAfterPtMinCut=0; 
  Int_t CounterSignPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t CounterBkgPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t TrueCounterSignPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t TrueCounterBkgPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t CounterSignPairsAllAssocMult[nummolt+1]={0}; 
  Int_t TrueCounterSignPairsAllAssocMult[nummolt+1]={0}; 


  Float_t ctauCasc[numtipo] = {2.6844,7.89, 7.89, 7.89};
  Float_t PDGCode[numtipo-1] = {310, 3122, -3122};
  Float_t LimInfMass[numtipo]= {0.4, 1, 1, 1};
  Float_t LimSupMass[numtipo]= {0.6, 1.2, 1.2, 1.2};
  Float_t LimInfMassTight[numtipo]= {0.491, 1, 1, 1};
  Float_t LimSupMassTight[numtipo]= {0.505, 1.2, 1.2, 1.2};

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  TString Smolt5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "all"};
  Double_t Nmolt5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};

  if (isHM){
    Nmolt[1] = 0.001;
    Nmolt[2] = 0.005;
    Nmolt[3] = 0.01;
    Nmolt[4] = 0.05;
    Nmolt[5] = 0.1;
    Smolt[0] = "0-0.001";
    Smolt[1] = "0.001-0.005";
    Smolt[2] = "0.005-0.01";
    Smolt[3] = "0.01-0.05";
    Smolt[4] = "0.05-0.1";
    Smolt[5] = "0-0.1";
    if (MultBinning==1){
      Nmolt[1] = 0;
      Nmolt[2] = 0;
      Nmolt[3] = 0.01;
      Nmolt[4] = 0.05;
      Nmolt[5] = 0.1;
      Smolt[0] = "0-0a";
      Smolt[1] = "0-0b";
      Smolt[2] = "0-0.01";
      Smolt[3] = "0.01-0.05";
      Smolt[4] = "0.05-0.1";
      Smolt[5] = "0-0.1";
    }
  }
  if (ispp5TeV && MultBinning==3) {
    for (Int_t m=0; m<nummolt+1; m++){
    Nmolt[m] = Nmolt5TeV[m];
    Smolt[m] = Smolt5TeV[m];
    }
  }

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  TString SPtV00[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};

  if (type>0)SPtV0[1]={"0.5-1"};
  else {
    SPtV0[0]={"0-0.5"}; 
    SPtV0[1]={"0.5-1"};
  }
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8,100};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }


  TString SPtV01[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};

  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1) numPtV0Max = numPtV0;
  else numPtV0Max = numPtV0-1;
  if (PtBinning==1){
    for(Int_t v=PtBinMin; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=PtBinMin; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }

  Float_t PtTrigMultDep[nummolt+1] = {3.0, 3.05, 3.05, 3.1, 3.1, 3};
  Double_t NPtTrigger[numPtTrigger+1]={PtTrigMin,ptjmax};
  for (Int_t m=0; m<nummolt+1; m++){
    if (IsPtTrigMultDep) PtTrigMultDep[m] = PtTrigMin;
  }
 

  TString ContainerName ="";
  if (isNewInputPath) {
    if (isMC && isEfficiency) ContainerName= "_hK0s_Task_RecoAndEfficiency";
    else     if (year == "2016k_HM_hK0s"|| year.Index("2019h11")!=-1 )  ContainerName= "_hK0s_Task_suffix";
    else if (year.Index("17pq_hK0s")!=-1 && isMC) ContainerName= "_hK0s_Task_RecoAndEfficiency";
    //    else if (year.Index("AOD235")!=-1 && isMC) ContainerName= "_hK0s_Task_RecoAndEfficiency";
    //else if (isMC && isEfficiency) ContainerName= "_hK0s_Task_RecoAndEfficiency";
    else  ContainerName= "_hK0s_Task_suffix";
    ContainerName= "_hK0s_Task_suffix";
    //else  ContainerName= "_hK0s_Task_";
  }

  cout << ContainerName<< endl;
  TTree *tSign;
  TTree *tBkg;
  if (year == "17pq_hK0s_pttrig0.15" || year.Index("NoTrackLength")!=-1 || year.Index("NoOOBPileUp")!=-1 || year.Index("15g3c3")!=-1 || year.Index("NSigma5")!= -1 ||  year.Index("V0RadiusSelection")!= -1 || year.Index("Trig0Pt")!=-1 || year =="17d20bEPOS_hK0s_EtaEff" || year.Index("AOD234")!=-1 || year.Index("AOD264")!=-1 || year == "LHC17o_HM_INELgt0281961" || year == "LHC19h11_HM_INELgt0Bis286380"){
    tSign = (TTree *)d->Get("MyOutputContainer1"+ContainerName);
    tBkg  = (TTree *)d->Get("MyOutputContainer2"+ContainerName);
  }
  else{
    tSign = (TTree *)d->Get("fSignalTree");
    tBkg  = (TTree *)d->Get("fBkgTree");
  }
  if (!tSign || !tBkg) {cout << "Input tree not found " << endl; return;}
  //what is the fraction of AC events in each multiplicity class?
  TList *d1 = (TList*)d->Get("MyOutputContainer" + ContainerName);
  if (!d1) {cout << " container named MyOutputContainer" << ContainerName << " does not exist " << endl; return;}
  TH1F* fHistEventMult=(TH1F*)  d1->FindObject("fHistEventMult");
  if (!fHistEventMult) {cout << "no info about total number of INT7 events analyzed" << endl; return;}

  cout <<" bin label: " << fHistEventMult->GetXaxis()->GetBinLabel(7) << endl;
  Double_t TotEvtINT7 = fHistEventMult->GetBinContent(7);
  if (isHM) TotEvtINT7 = fHistEventMult->GetBinContent(9);

  TH2F* hMultiplicity2D=(TH2F*)  d1->FindObject("fHistPtMaxvsMult");
  TH2F* hMultiplicity2DBefAll=(TH2F*)  d1->FindObject("fHistPtMaxvsMultBefAll");
  TH1F* hMultiplicityBefAll;
  TH1F* hMultiplicity;
  TH2F* hMultvsNumberAssoc=(TH2F*)  d1->FindObject("fHistMultvsV0All");
  TH1F* hMultvsNumberAssoc_Proj[nummolt];
  if (!hMultiplicity2D) cout << " no info about multiplicity distribution of AC events available fHistPtMaxvsMult" << endl;
  if (!hMultiplicity2DBefAll) cout << " no info about multiplicity distribution of AC events available fHistPtMaxvsMultBefAll" << endl;
  if (!hMultvsNumberAssoc) cout << " no info about multiplicity distribution of AC events available fHistMultvsV0" << endl;
  Float_t ACcounter[nummolt+1];

  if (hMultiplicity2D && hMultvsNumberAssoc && hMultiplicity2DBefAll){
    hMultiplicity=(TH1F*)  hMultiplicity2D->ProjectionY("fHistPtMaxvsMult1D",     hMultiplicity2D->GetXaxis()->FindBin(PtTrigMin+0.0001),  hMultiplicity2D->GetXaxis()->FindBin(ptjmax-0.0001));
    hMultiplicityBefAll=(TH1F*)  hMultiplicity2DBefAll->ProjectionY("fHistPtMaxvsMult1DBefAll",     hMultiplicity2DBefAll->GetXaxis()->FindBin(PtTrigMin+0.0001),  hMultiplicity2DBefAll->GetXaxis()->FindBin(ptjmax-0.0001));
    if (year == "LHC17o_HM_INELgt0281961")     hMultiplicityBefAll=(TH1F*)d1->FindObject("fHist_multiplicity_EvwTrigger");
    ACcounter[5]= hMultiplicity->GetEntries();

    cout <<"total number of events with NT>0  " << hMultiplicityBefAll->GetEntries() << " " <<   (Float_t)hMultiplicityBefAll->GetEntries()/TotEvtINT7 << endl;
    cout <<"\ntotal number of events used in the AC (no selections on topo var and on skipAssoc!) " << hMultiplicity->GetEntries() << endl;

    for (Int_t m=0; m< nummolt; m++){ 
      hMultvsNumberAssoc_Proj[m] = (TH1F*)       hMultvsNumberAssoc->ProjectionX(Form("hMultvsNumberAssoc_%i", m), hMultvsNumberAssoc->GetYaxis()->FindBin(Nmolt[m]+0.0001), hMultvsNumberAssoc->GetYaxis()->FindBin(Nmolt[m+1]-0.0001));
      hMultvsNumberAssoc_Proj[m]->GetXaxis()->SetRangeUser(1,      hMultvsNumberAssoc_Proj[m]->GetXaxis()->GetXmax());
      cout << " m " << m << endl;
      ACcounter[m] =0;
      for (Int_t b= hMultiplicity->GetXaxis()->FindBin(Nmolt[m]+0.0001); b <=  hMultiplicity->GetXaxis()->FindBin(Nmolt[m+1]-0.0001); b++){
	ACcounter[m] +=  hMultiplicity->GetBinContent(b);
      }
      cout << "fraction of events in mult bin (for ptTrig> PtTrigMin) " << Smolt[m] << ": " << ACcounter[m]/ACcounter[5] <<  " ~average V0 number (for pTTrig> 0.150 GeV usually) " <<     hMultvsNumberAssoc_Proj[m]->GetMean() <<endl;
    }
    
  }

  
  //Signal

  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);		       
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);		       
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);		       
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);		       
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);			       
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);			       
  tSign->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger);               
  tSign->SetBranchAddress("fTreeVariablePDGCodeTrigger"            ,&fSignTreeVariablePDGCodeTrigger);          

  tSign->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0);                   
  tSign->SetBranchAddress("fTreeVariableChargeAssoc"               ,&fSignTreeVariableChargeAssoc);		       
  tSign->SetBranchAddress("fTreeVariablePDGCodeAssoc"              ,&fSignTreeVariablePDGCodeAssoc);             
  tSign->SetBranchAddress("fTreeVariableSkipAssoc"                 ,&fSignTreeVariableSkipAssoc);                     
  tSign->SetBranchAddress("fTreeVariableIsCommonParton"                   ,&fSignTreeVariableIsCommonParton);
  tSign->SetBranchAddress("fTreeVariablePdgCommonParton"                   ,&fSignTreeVariablePdgCommonParton);
  tSign->SetBranchAddress("fTreeVariablePtV0"                      ,&fSignTreeVariablePtV0);
  tSign->SetBranchAddress("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0);
  tSign->SetBranchAddress("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0);

  tSign->SetBranchAddress("fTreeVariableAssocDCAz"                 ,&fSignTreeVariableAssocDCAz);         
  tSign->SetBranchAddress("fTreeVariableAssocDCAxy"                ,&fSignTreeVariableAssocDCAxy);              

  tSign->SetBranchAddress("fTreeVariableRapK0Short"                ,&fSignTreeVariableRapK0Short);
  tSign->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex"         ,&fSignTreeVariableDcaV0ToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableDcaPosToPrimVertex"        ,&fSignTreeVariableDcaPosToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableDcaNegToPrimVertex"        ,&fSignTreeVariableDcaNegToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fSignTreeVariableV0CosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariablectau"                      ,&fSignTreeVariablectau);
  tSign->SetBranchAddress("fTreeVariableInvMassK0s"                ,&fSignTreeVariableInvMassK0s);
  tSign->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fSignTreeVariableInvMassLambda);
  tSign->SetBranchAddress("fTreeVariableInvMassAntiLambda"         ,&fSignTreeVariableInvMassAntiLambda);      
  tSign->SetBranchAddress("fTreeVariablePtArmenteros"              ,&fSignTreeVariablePtArmenteros);
  tSign->SetBranchAddress("fTreeVariableAlpha"                     ,&fSignTreeVariableAlpha);

  tSign->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta);
  tSign->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi);
  tSign->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta);

  tSign->SetBranchAddress("fTreeVariableMultiplicity"              ,&fSignTreeVariableMultiplicity);
  tSign->SetBranchAddress("fTreeVariableZvertex"                   ,&fSignTreeVariableZvertex);
   


  //BackGround

  tBkg->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger);		       
  tBkg->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz);			       
  tBkg->SetBranchAddress("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy);			       
  tBkg->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger);               
  tBkg->SetBranchAddress("fTreeVariablePDGCodeTrigger"            ,&fBkgTreeVariablePDGCodeTrigger);          

  tBkg->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0);                   
  tBkg->SetBranchAddress("fTreeVariableChargeAssoc"               ,&fBkgTreeVariableChargeAssoc);		       
  tBkg->SetBranchAddress("fTreeVariablePDGCodeAssoc"              ,&fBkgTreeVariablePDGCodeAssoc);             
  tBkg->SetBranchAddress("fTreeVariableSkipAssoc"                 ,&fBkgTreeVariableSkipAssoc);                     
  tBkg->SetBranchAddress("fTreeVariableIsCommonParton"                   ,&fBkgTreeVariableIsCommonParton);
  tBkg->SetBranchAddress("fTreeVariablePdgCommonParton"                   ,&fBkgTreeVariablePdgCommonParton);
  tBkg->SetBranchAddress("fTreeVariablePtV0"                      ,&fBkgTreeVariablePtV0);
  tBkg->SetBranchAddress("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0);
  tBkg->SetBranchAddress("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0);

  tBkg->SetBranchAddress("fTreeVariableAssocDCAz"                 ,&fBkgTreeVariableAssocDCAz);         
  tBkg->SetBranchAddress("fTreeVariableAssocDCAxy"                ,&fBkgTreeVariableAssocDCAxy);              

  tBkg->SetBranchAddress("fTreeVariableRapK0Short"                ,&fBkgTreeVariableRapK0Short);
  tBkg->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex"         ,&fBkgTreeVariableDcaV0ToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableDcaPosToPrimVertex"        ,&fBkgTreeVariableDcaPosToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableDcaNegToPrimVertex"        ,&fBkgTreeVariableDcaNegToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fBkgTreeVariableV0CosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariablectau"                      ,&fBkgTreeVariablectau);
  tBkg->SetBranchAddress("fTreeVariableInvMassK0s"                ,&fBkgTreeVariableInvMassK0s);
  tBkg->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fBkgTreeVariableInvMassLambda);
  tBkg->SetBranchAddress("fTreeVariableInvMassAntiLambda"         ,&fBkgTreeVariableInvMassAntiLambda);      
  tBkg->SetBranchAddress("fTreeVariablePtArmenteros"              ,&fBkgTreeVariablePtArmenteros);
  tBkg->SetBranchAddress("fTreeVariableAlpha"                     ,&fBkgTreeVariableAlpha);

  tBkg->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta);
  tBkg->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi);
  tBkg->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta);

  tBkg->SetBranchAddress("fTreeVariableMultiplicity"              ,&fBkgTreeVariableMultiplicity);
  tBkg->SetBranchAddress("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex);
   
  

  Int_t EntriesSign = 0; 
  Int_t EntriesBkg  = 0; 
  
  TFile *fout = new TFile(PathOut,"RECREATE");

  Int_t numMultBins =100;
  Float_t UpperLimitMult=100;
  if (isHM) UpperLimitMult = 0.1;

  TH1F* hMultiplicityhK0sPairs = new TH1F("hMultiplicityhK0sPairs", "hMultiplicityhK0sPairs", 100, 0, UpperLimitMult);
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[nummolt+1][numzeta][numPtTrigger];

  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];


  TH1F * HistoInfo = new TH1F("HistoInfo", "HistoInfo",25,0.5, 25.5);
  HistoInfo->GetXaxis()->SetBinLabel(1, "% Ev. NT>0");
  HistoInfo->GetXaxis()->SetBinLabel(2, "NV0/NInt7");
  HistoInfo->GetXaxis()->SetBinLabel(3, "<p_{T,Xi}>");
  HistoInfo->GetXaxis()->SetBinLabel(4,"<p_{T,Trig}>");
  HistoInfo->GetXaxis()->SetBinLabel(5,"SE pairs");
  HistoInfo->GetXaxis()->SetBinLabel(6,"ME pairs");
  for (Int_t m=0; m< nummolt; m++){
    HistoInfo->GetXaxis()->SetBinLabel(7+m,Form("SEpairs Mult distr m%i",m));
    HistoInfo->GetXaxis()->SetBinLabel(12+m, Form("NV0/ACev m%i",m));
    HistoInfo->GetXaxis()->SetBinLabel(17+m,Form("fraction of non skipped V0 m%i",m));
  }


  TH2F*  fHistPtAssocvsPtTrig[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPtAssocvsPtTrig[j] = new TH2F (Form("fHistPtAssocvsPtTrig%i",j),Form("fHistPtAssocvsPtTrig%i",j),100, 0, 30, 100, 0, 10);
  }

  //--------histograms to get some info from K0s candidate having a common origin to the trigger particle----------------     
  TH1D *hSign_PtAssoc_CPTrue       = new TH1D("hSign_PtAssoc_CPTrue",       "hSign_PtAssoc_CPTrue",       300, 0, 30);
  hSign_PtAssoc_CPTrue->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssoc_CPTrue->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssoc_CPTrue->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssoc_CPTrue->GetYaxis()->SetLabelSize(0.05);

  TH1D *hSign_PtAssoc_NOCPTrue       = new TH1D("hSign_PtAssoc_NOCPTrue",       "hSign_PtAssoc_NOCPTrue",       300, 0, 30);
  hSign_PtAssoc_NOCPTrue->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssoc_NOCPTrue->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssoc_NOCPTrue->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssoc_NOCPTrue->GetYaxis()->SetLabelSize(0.05);

  //------------------Histograms os selected particles (V0) for future efficiency calculation ----------------
  TH3F*    fHistSelectedV0PtTMaxPhi=new TH3F(Form("fHistSelectedV0PtTMaxPhi_%i", sysV0), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi() ,  numMultBins, 0, UpperLimitMult);
  fHistSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F*    fHistCPSelectedV0PtTMaxPhi=new TH3F(Form("fHistCPSelectedV0PtTMaxPhi_%i", sysV0), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi() ,  numMultBins, 0, UpperLimitMult);
  fHistCPSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistCPSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F*    fHistNOCPSelectedV0PtTMaxPhi=new TH3F(Form("fHistNOCPSelectedV0PtTMaxPhi_%i", sysV0), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi() ,  numMultBins, 0, UpperLimitMult);
  fHistNOCPSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistNOCPSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F *    fHistSelectedV0PtTMaxEta=new TH3F(Form("fHistSelectedV0PtTMaxEta_%i", sysV0), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2,  numMultBins, 0, UpperLimitMult );
  fHistSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistCPSelectedV0PtTMaxEta=new TH3F(Form("fHistCPSelectedV0PtTMaxEta_%i", sysV0), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles with common parton (Casc, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2,  numMultBins, 0,UpperLimitMult );
  fHistCPSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistCPSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistNOCPSelectedV0PtTMaxEta=new TH3F(Form("fHistNOCPSelectedV0PtTMaxEta_%i", sysV0), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles with no common parton(Casc, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2, numMultBins, 0, UpperLimitMult );
  fHistNOCPSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistNOCPSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistSelectedV0PtPtTMax=new TH3F(Form("fHistSelectedV0PtPtTMax_%i",sysV0), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles (Casc, primary, events w T>0)", 300, 0, 30, 120, -30,30,  numMultBins, 0, UpperLimitMult );
  fHistSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F *    fHistCPSelectedV0PtPtTMax=new TH3F(Form("fHistCPSelectedV0PtPtTMax_%i",sysV0), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles with common parton (Casc, primary, events w T>0)", 300, 0, 30, 120, -30,30,  numMultBins, 0, UpperLimitMult  );
  fHistCPSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistCPSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F *    fHistNOCPSelectedV0PtPtTMax=new TH3F(Form("fHistNOCPSelectedV0PtPtTMax_%i",sysV0), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles with no common parton (Casc, primary, events w T>0)", 300, 0, 30, 120, -30,30,  numMultBins, 0, UpperLimitMult );
  fHistNOCPSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistNOCPSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F*    fHistSelectedV0PtEta=new TH3F(Form("fHistSelectedV0PtEta_%i", sysV0), "p_{T} and #eta distribution of selected V0 particles (Casc, primary, events w T>0)", 300, 0, 30, /*400*/480, -1.2, 1.2,  numMultBins, 0, UpperLimitMult);
  fHistSelectedV0PtEta->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtEta->GetYaxis()->SetTitle("#eta");

  TH3F*    fHistPrimaryV0[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPrimaryV0[j]=new TH3F(Form("fHistPrimaryV0_%i_cut%i",j, sysV0), "V0 MC (Casc, selected)", 4, 0.5, 4.5, 160, 0, 16, 60, 0, 30);
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s");
    fHistPrimaryV0[j]->GetYaxis()->SetTitle("p_{T}");
    fHistPrimaryV0[j]->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }

  //------------------Histograms os selected particles (trigger) for future efficiency calculation ----------------
  TH3F*      fHistSelectedTriggerPtPhi=new TH3F(Form("fHistSelectedTriggerPtPhi_%i",sysTrigger), "p_{T} and #phi distribution of selected trigger particles (primary)", 600, 0, 30, 400,0   , 2*TMath::Pi() ,  numMultBins, 0, UpperLimitMult);
  fHistSelectedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  TH3F*     fHistSelectedTriggerPtEta=new TH3F(Form("fHistSelectedTriggerPtEta_%i",sysTrigger), "p_{T} and #eta distribution of selected trigger particles (primary)", 600, 0, 30, 400,   -1.2, 1.2,  numMultBins, 0, UpperLimitMult);
  fHistSelectedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtEta->GetYaxis()->SetTitle("#eta");


  TH2F*    fHistPrimaryTrigger[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPrimaryTrigger[j]=new TH2F(Form("fHistPrimaryTrigger_%i_cut%i", j, sysTrigger), "Trigger MC (selected)", 4, 0.5, 4.5, 100, 0,30 );
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(1,"Primary selected triggers");
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected triggers");
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected triggers");
    fHistPrimaryTrigger[j]->GetYaxis()->SetTitle("p_{T}");
  }
 
  /*-----------------------Pt Trigger  --------------------------- */
  TH1D *hSign_PtTrigger       = new TH1D("hSign_PtTrigger",       "hSign_PtTrigger",       300, 0, 30);
  hSign_PtTrigger->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtTrigger->GetXaxis()->SetTitleSize(0.05);
  hSign_PtTrigger->GetXaxis()->SetLabelSize(0.05);
  hSign_PtTrigger->GetYaxis()->SetLabelSize(0.05);
  TH1D *hBkg_PtTrigger         = new TH1D("hBkg_PtTrigger",        "hBkg_PtTrigger",        300, 0, 30);
  hBkg_PtTrigger->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hBkg_PtTrigger->GetXaxis()->SetTitleSize(0.05);
  hBkg_PtTrigger->GetXaxis()->SetLabelSize(0.05);
  hBkg_PtTrigger->GetYaxis()->SetLabelSize(0.05);
  /*-----------------------Pt Assoc  --------------------------- */
  TH1D *hSign_PtAssoc       = new TH1D("hSign_PtAssoc",       "hSign_PtAssoc",       300, 0, 30);
  hSign_PtAssoc->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssoc->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssoc->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssoc->GetYaxis()->SetLabelSize(0.05);

  TH1D *hSign_PtAssocTrue       = new TH1D("hSign_PtAssocTrue",       "hSign_PtAssocTrue",       300, 0, 30);
  hSign_PtAssocTrue->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssocTrue->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssocTrue->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssocTrue->GetYaxis()->SetLabelSize(0.05);

  TH2D *hSign_PtTriggerPtAssoc = new TH2D ("hSign_PtTriggerPtAssoc", "hSign_PtTriggerPtAssoc", 100, 0, 10, 300, 0 ,30);
  hSign_PtTriggerPtAssoc->GetXaxis()->SetTitle("p_{T, Assoc} (Gev/c)");
  hSign_PtTriggerPtAssoc->GetYaxis()->SetTitle("p_{T, Trig} (Gev/c)");

  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hMassvsPt_SEbins[nummolt+1][numzeta][numPtTrigger];
  TH2D *hMassvsPt_SEbins_true[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[nummolt+1][numzeta][numPtTrigger];

  //Form("hMassvsPt_"+tipo[type]+"_%i",molt
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	hMassvsPt_SEbins[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i", m), 100, LimInfMass[type], LimSupMass[type], 100, 0, 10);
	hMassvsPt_SEbins_true[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m), 100, LimInfMass[type], LimSupMass[type], 100, 0, 10);
	hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_CPTrue_PtInt", "SE_m"+Smolt[m]+"_CPTrue_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_NOCPTrue_PtInt", "SE_m"+Smolt[m]+"_NOCPTrue_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hDeltaEtaDeltaPhi_SEbins_CPPtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_CP_PtInt", "SE_m"+Smolt[m]+"_CP_PtInt", 56, -1.5, 1.5,104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
        hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_NOCP_PtInt", "SE_m"+Smolt[m]+"_NOCP_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  nameSE[m][z][v][tr]="SE_";
	  namemassSE[m][z][v][tr]="InvMassSE_";
	  //nameSE[m][z][v][tr]+=Form("m%i_z%i_v%i_tr%i",m,z,v,tr);
	  nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_SEbins[m][z][v][tr]= new TH1D(namemassSE[m][z][v][tr], namemassSE[m][z][v][tr],100, LimInfMass[type], LimSupMass[type]);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr], nameSE[m][z][v][tr],   56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi()); 
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]= (TH2D*)   hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_CPTrue");
          hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]->SetTitle("AC for K0s (true) cand. with a common parton" +  SPtV0[v]);
          hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]= (TH2D*)       hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_CP");
          hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]->SetTitle("AC for K0s cand. with a common parton"+ SPtV0[v]);
          hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]= (TH2D*)         hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_NOCPTrue");
          hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]->SetTitle("AC for K0s (true) cand. without a common parton" + SPtV0[v]);
          hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]= (TH2D*)     hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_NOCP");
          hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->SetTitle("AC for K0s cand. without a common parton" +  SPtV0[v]);

	}
      }
    }
  }

  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hMassvsPt_MEbins[nummolt+1][numzeta][numPtTrigger];
  TH2D *hMassvsPt_MEbins_true[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	hMassvsPt_MEbins[m][z][tr]= new TH2D(Form("ME_hMassvsPt_"+tipo[type]+"_%i", m),Form("ME_hMassvsPt_"+tipo[type]+"_%i", m),100, LimInfMass[type], LimSupMass[type],100, 0, 10);
	hMassvsPt_MEbins_true[m][z][tr]= new TH2D(Form("ME_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("ME_hMassvsPt_"+tipo[type]+"_%i_true", m),100, LimInfMass[type], LimSupMass[type],100, 0, 10);
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  nameME[m][z][v][tr]="ME_";
	  namemassME[m][z][v][tr]="InvMassME_";
	  nameME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_MEbins[m][z][v][tr]= new TH1D(namemassME[m][z][v][tr], namemassME[m][z][v][tr],   100, LimInfMass[type], LimSupMass[type]);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]= new TH2D(nameME[m][z][v][tr], nameME[m][z][v][tr],  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);
	
	}
      }
    }
  }

  EntriesSign =  tSign->GetEntries();
  EntriesBkg  =  tBkg ->GetEntries();
    
  Int_t InJet[nummolt+1]={0};
  Int_t OutJet[nummolt+1]={0};

  Bool_t MoltSel=kFALSE; 
  Float_t     fSignTreeVariableInvMass= 0;
  Bool_t isParticleTrue=kFALSE;
  Int_t fSignTreeVariablePAPAssoc=0;
  // cout << "  entries Sign: " << EntriesSign<<endl;
  Int_t l=0;
  for(Int_t k = 0; k<EntriesSign; k++){
    //    cout << k << endl;
    //    if (k>100000) continue;
    tSign->GetEntry(k);  
    if (k==100000*l){
      l++;     
      cout << "SE ----" << (Float_t)k/EntriesSign<< endl;
    }
    //    if (fSignTreeVariableMultiplicity == 0) cout << "Mult ==0 " << endl;
    if (fSignTreeVariableMultiplicity == 0) continue;
    fSignTreeVariableDeltaEta=fSignTreeVariableEtaV0-fSignTreeVariableEtaTrigger;

    //    cout << " mult: " <<fSignTreeVariableMultiplicity << endl;
    //charge selection: not done, since both K0s and Lambdas have charge =0
    /*
      if ((type==0 || type ==2) && fSignTreeVariableChargeAssoc==1) continue;
      else if ((type==1 || type ==3) && fSignTreeVariableChargeAssoc==-1) continue;
    */

    //particle - antiparticle selection
    if (type==0 || type ==1) fSignTreeVariablePAPAssoc=1;
    else if (type==2)  fSignTreeVariablePAPAssoc=-1;
    //else ?

    //inv mass definition
    fSignTreeVariableInvMass= 0;
    if (type==0)     fSignTreeVariableInvMass= fSignTreeVariableInvMassK0s;
    else if (type==2)      fSignTreeVariableInvMass= fSignTreeVariableInvMassLambda;
    else if (type==3)      fSignTreeVariableInvMass= fSignTreeVariableInvMassAntiLambda;
    //    else ?

    //    cout << " 1 " << endl;
    MassInPeakK0s=0;
    if (CommonParton){
      if ( fSignTreeVariableInvMass> LimInfMassTight[type] &&  fSignTreeVariableInvMass<LimSupMassTight[type]) MassInPeakK0s=1;
    }
    //    cout << "MassInPeakK0s" << MassInPeakK0s << " " << fSignTreeVariableInvMass <<  " " << LimSupMassTight[type] << endl;
    //    cout << " 1 " << endl;
    //rapidity selection
    if (israp==0 && TMath::Abs(fSignTreeVariableEtaV0)>0.8)continue; 
    else if (israp==1 && TMath::Abs(fSignTreeVariableRapK0Short)>0.5)continue;

    //definition of true particle
    if (ishhCorr)  isParticleTrue=kTRUE; //the request of being primary will be made afterwards
    else if (type<=2) isParticleTrue= (fSignTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==3) isParticleTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[1]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[1])  );

    //    if (!isParticleTrue) continue;
    //************cuts on pT trigger min*********************************

    if (!isINEL){
      if (!IsPtTrigMultDep){
	if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMin) continue;
      }
      else {
	Bool_t IsSelectionPtMinNotPassed=0;
	for(Int_t m=0; m<nummolt; m++){
	  if ( fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]){
	    if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMultDep[m]) IsSelectionPtMinNotPassed=1;
	  }
	}
	if (IsSelectionPtMinNotPassed) continue;
      }

      if(TMath::Abs(fSignTreeVariablePtTrigger)>ptjmax) continue;
    }
    if(isMC==0 || (isMC==1 && isEfficiency==1)){
      if (!isINEL){
	//cuts on DCAz trigger*******************
	if(sysTrigger==0){
	  if(TMath::Abs(fSignTreeVariableDCAz)>1) continue;
	}
	if(sysTrigger==1){
	  if(TMath::Abs(fSignTreeVariableDCAz)>2) continue;
	}
	if(sysTrigger==2){
	  if(TMath::Abs(fSignTreeVariableDCAz)>0.5) continue;
	}
      }
      //******************* some other cuts for sys studies**************************
      if (!ishhCorr){
	  if (type==0){
	    if(TMath::Abs((fSignTreeVariableInvMassLambda - massLambda))<= LambdaMassWindow) continue;
	    if(TMath::Abs((fSignTreeVariableInvMassAntiLambda - massLambda))<= LambdaMassWindow) continue;
	  }
	  if(fSignTreeVariableV0CosineOfPointingAngle<CosinePAngle)            continue;
	  if(fSignTreeVariableDcaNegToPrimVertex < DCANegToPV)      continue;
	  if(fSignTreeVariableDcaPosToPrimVertex < DCAPosToPV)      continue;
	  if(fSignTreeVariableDcaV0ToPrimVertex > DCAV0ToPV)      continue;
	  if(fSignTreeVariablectau> ctau)                   continue;
      }
      else if (ishhCorr){
	if(sysV0==0){
	  if(TMath::Abs(fSignTreeVariableAssocDCAz)>1) continue;
	}
	if(sysV0==1){
	  if(TMath::Abs(fSignTreeVariableAssocDCAz)>2) continue;
	}
	if(sysV0==2){
	  if(TMath::Abs(fSignTreeVariableAssocDCAz)>0.5) continue;
	}
      }
    }

    //pT trigger vs pT assoc histogram (before pt selection on assoc)
    hMultiplicityhK0sPairs->Fill(fSignTreeVariableMultiplicity);
    for(Int_t m=0; m<nummolt+1; m++){
      if (m < nummolt){
	MoltSel = (fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]);
      }
      else if(m==nummolt) MoltSel=kTRUE;
      // cout << m << "  " << nummolt << "  " << MoltSel<< endl;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  if (MoltSel) {
	    fHistPtAssocvsPtTrig[m] ->Fill(fSignTreeVariablePtTrigger, fSignTreeVariablePtV0);
	    CounterSignPairsAllAssocMult[m]++;  
	    if (isParticleTrue)      TrueCounterSignPairsAllAssocMult[m]++;
	  }
	}
      }
    }


    if (SkipAssoc){    if (fSignTreeVariableSkipAssoc==1) continue;}

    if (isParticleTrue)      TrueCounterSignPairsAfterPtMinCut++;  
    CounterSignPairsAfterPtMinCut++;  
    //**********************************************************************************************

    fSignTreeVariableDeltaPhi= -fSignTreeVariableDeltaPhi;
    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();


    if (Region>0){
      if (Region==1){
	if (TMath::Abs(fSignTreeVariableDeltaPhi)> 0.85) continue; 
	if (TMath::Abs(fSignTreeVariableDeltaEta)> 0.75) continue; 
      }
      else       if (Region==2){
	if (TMath::Abs(fSignTreeVariableDeltaPhi)< 0.85 || TMath::Abs(fSignTreeVariableDeltaPhi)> 2.05) continue; 
	if (TMath::Abs(fSignTreeVariableDeltaEta)< 0.75 || TMath::Abs(fSignTreeVariableDeltaEta)> 1.18) continue;       }
    }

    if (isMC && isEfficiency) {
      if(fSignTreeVariableisPrimaryTrigger==1){
	fHistSelectedTriggerPtPhi->Fill(fSignTreeVariablePtTrigger,fSignTreeVariablePhiTrigger, fSignTreeVariableMultiplicity);
	fHistSelectedTriggerPtEta->Fill(fSignTreeVariablePtTrigger,fSignTreeVariableEtaTrigger, fSignTreeVariableMultiplicity);
      }

      for (Int_t m =0; m<5;m++){
	if(fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]){
	  for(Int_t p=1; p<=4; p++){
	    if (fSignTreeVariableisPrimaryTrigger==p){
	      fHistPrimaryTrigger[m]->Fill(p,fSignTreeVariablePtTrigger );
	    }
	  }
	}
      }
      for(Int_t p=1; p<=4; p++){
	if (fSignTreeVariableisPrimaryTrigger==p){
	  fHistPrimaryTrigger[5]->Fill(p,fSignTreeVariablePtTrigger );
	}
      }

      Bool_t OutMass=kFALSE;
      if(isParticleTrue &&  fSignTreeVariableisPrimaryV0==1){

        if (CommonParton){
          if ( fSignTreeVariableIsCommonParton)      {
            fHistCPSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
            fHistCPSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc ,  fSignTreeVariableMultiplicity);
            fHistCPSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
          }
          else     {
            fHistNOCPSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
            fHistNOCPSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc  , fSignTreeVariableMultiplicity);
            fHistNOCPSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
          }
        }

	//mass selections
	if (isEfficiencyMassSel){
	  for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	    if(fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      OutMass= (TMath::Abs(fSignTreeVariableInvMass - mass[type][5][v])>=4*sigma[type][5][v]);
	      //	      cout << mass[type][5][v] << " " << sigma[type][5][v] << " " << fSignTreeVariableInvMass<<endl;
	    }
	  }
	}

	if (!OutMass){
	  fHistSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
	  fHistSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
	  fHistSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger*fSignTreeVariablePAPAssoc , fSignTreeVariableMultiplicity);
	  fHistSelectedV0PtEta->Fill(fSignTreeVariablePtV0,fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
	}
      }
    

      for (Int_t m =0; m<5;m++){
	if(fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]){
	  for(Int_t p=1; p<=4; p++){
	    if (fSignTreeVariableisPrimaryV0==p){
	      if (!OutMass){
	      fHistPrimaryV0[m]->Fill(p,fSignTreeVariablePtV0, fSignTreeVariablePtTrigger );
	      }
	    }
	  }
	}
      }
      for(Int_t p=1; p<=4; p++){
	if (fSignTreeVariableisPrimaryV0==p){
	  if (!OutMass){
	    fHistPrimaryV0[5]->Fill(p,fSignTreeVariablePtV0, fSignTreeVariablePtTrigger);
	  }
	}
      }
    }

    hSign_PtAssoc->Fill(fSignTreeVariablePtV0);
    hSign_PtTrigger->Fill(fSignTreeVariablePtTrigger);
    hSign_PtTriggerPtAssoc->Fill(fSignTreeVariablePtV0, fSignTreeVariablePtTrigger);
    if (isParticleTrue)hSign_PtAssocTrue->Fill(fSignTreeVariablePtV0);
    if ((isMC && isParticleTrue) || !isMC){
      if (CommonParton){
        if (fSignTreeVariableIsCommonParton)      hSign_PtAssoc_CPTrue->Fill(fSignTreeVariablePtV0);
        else      hSign_PtAssoc_NOCPTrue->Fill(fSignTreeVariablePtV0);
      }
    }


    for(Int_t m=0; m<nummolt+1; m++){
      if (m < nummolt){
	MoltSel = (fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]);
      }
      else if(m==nummolt) MoltSel=kTRUE;
      // cout << m << "  " << nummolt << "  " << MoltSel<< endl;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  if(MoltSel){
	    CounterSignPairsAfterPtMinCutMult[m]++;  
	    if (isParticleTrue)      TrueCounterSignPairsAfterPtMinCutMult[m]++;  
	    if (!ishhCorr){
	      hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMass, fSignTreeVariablePtV0); 
	      if(isParticleTrue) hMassvsPt_SEbins_true[m][z][tr]->Fill(fSignTreeVariableInvMass, fSignTreeVariablePtV0); 
	    }
	    if (CommonParton){
              if ((isMC && isParticleTrue)){

                if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
                else   hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
              }
	      if (MassInPeakK0s){
		//	      if (!isParticleTrue){
		//	      if (kTRUE){
              if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CPPtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta,fSignTreeVariableDeltaPhi);
              else   hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }
              if (fSignTreeVariableIsCommonParton ){
                if (TMath::Abs(fSignTreeVariableDeltaPhi)<1) InJet[m]++;
                else OutJet[m]++;
              }
            }
	  }

	  for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	    if(MoltSel && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      if (!ishhCorr)	      hInvMassK0Short_SEbins[m][z][v][tr]->Fill(fSignTreeVariableInvMass);

	      if (CommonParton){
                if ((isMC && isParticleTrue) || !isMC){
                  if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
                  else   hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
                }
		if (MassInPeakK0s){
		  //if (!isParticleTrue){
		  //if (kTRUE){
                if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta,fSignTreeVariableDeltaPhi);
                else   hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		}
              }

	    }
	  }
	}
      }
    }
  }


  Float_t     fBkgTreeVariableInvMass= 0;
  l=0;
  if (!isSysStudy & !isINEL){
  for(Int_t k = 0; k<EntriesBkg; k++){
    // for(Int_t k = 0; k<1; k++){
    //if (k>100000) continue;
    tBkg->GetEntry(k);
    if (k==100000*l){
      l++;     
      cout << "ME ----" << (Float_t)k/EntriesBkg<< endl;
    }
    fBkgTreeVariableDeltaEta=fBkgTreeVariableEtaV0-fBkgTreeVariableEtaTrigger;

    //charge selection
    /*
      if ((type==0 || type ==2) && fBkgTreeVariableChargeAssoc==1) continue;
      else if ((type==1 || type ==3) && fBkgTreeVariableChargeAssoc==-1) continue;
    */
    //inv mass definition
    fBkgTreeVariableInvMass= 0;
    if (type==0)     fBkgTreeVariableInvMass= fBkgTreeVariableInvMassK0s;
    else if (type==2)      fBkgTreeVariableInvMass= fBkgTreeVariableInvMassLambda;
    else if (type==3)      fBkgTreeVariableInvMass= fBkgTreeVariableInvMassAntiLambda;

    //rapidity selection
    if (israp==0 && TMath::Abs(fBkgTreeVariableEtaV0)>0.8)continue;
    else if (israp==1 && TMath::Abs(fBkgTreeVariableRapK0Short)>0.5)continue;

    //definition of true particle
    if (ishhCorr) isParticleTrue =kTRUE;
    else     if (type<=2) isParticleTrue= (fBkgTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==3) isParticleTrue= ((fBkgTreeVariablePDGCodeAssoc==PDGCode[1]) ||(fBkgTreeVariablePDGCodeAssoc==PDGCode[1])  );
 
    //    if (!isParticleTrue) continue;
    if (SkipAssoc){    if (fBkgTreeVariableSkipAssoc==1) continue;}

    //************cuts on pT trigger min*********************************

    if (!IsPtTrigMultDep){
      if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMin) continue;
    }
    else {
      Bool_t IsSelectionPtMinNotPassedBkg=0;
      for(Int_t m=0; m<nummolt; m++){
	if ( fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1]){
	  if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMultDep[m]) IsSelectionPtMinNotPassedBkg=1;
	}
      }
      if (IsSelectionPtMinNotPassedBkg) continue;
    }


    if(TMath::Abs(fBkgTreeVariablePtTrigger)>ptjmax) continue;

    if(isMC==0 || (isMC==1 && isEfficiency==1)){

      //cuts on DCAz trigger*******************
      if(sysTrigger==0){
	if(TMath::Abs(fBkgTreeVariableDCAz)>1) continue;
      }
      if(sysTrigger==1){
	if(TMath::Abs(fBkgTreeVariableDCAz)>2) continue;
      }
      if(sysTrigger==2){
	if(TMath::Abs(fBkgTreeVariableDCAz)>0.5) continue;
      }

      //******************* some other cuts for sys studies **************************
      if (!ishhCorr){
	if (sysV0==0){
	  //the values are valid for K0s
	  if (type==0){
	    if(TMath::Abs((fBkgTreeVariableInvMassLambda - massLambda))<= LambdaMassWindow) continue;
	    if(TMath::Abs((fBkgTreeVariableInvMassAntiLambda - massLambda))<= LambdaMassWindow) continue;
	  }
	  if(fBkgTreeVariableV0CosineOfPointingAngle<CosinePAngle)            continue;
	  if(fBkgTreeVariableDcaNegToPrimVertex < DCANegToPV)      continue;
	  if(fBkgTreeVariableDcaPosToPrimVertex < DCAPosToPV)      continue;
	  if(fBkgTreeVariableDcaV0ToPrimVertex > DCAV0ToPV)      continue;
	  if(fBkgTreeVariablectau> ctau)                   continue;
	}
      }
      else if (ishhCorr){
	if(sysV0==0){
	  if(TMath::Abs(fBkgTreeVariableAssocDCAz)>1) continue;
	}
	if(sysV0==1){
	  if(TMath::Abs(fBkgTreeVariableAssocDCAz)>2) continue;
	}
	if(sysV0==2){
	  if(TMath::Abs(fBkgTreeVariableAssocDCAz)>0.5) continue;
	}
      }

    }
    if (isParticleTrue)      TrueCounterBkgPairsAfterPtMinCut++;  
    CounterBkgPairsAfterPtMinCut++;  
    //**********************************************************************************************

    fBkgTreeVariableDeltaPhi= -fBkgTreeVariableDeltaPhi;
    if (fBkgTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fBkgTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fBkgTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fBkgTreeVariableDeltaPhi += 2.0*TMath::Pi();
    /*
    if (EffRegion>0){
      if (EffRegion==1){
	if (TMath::Abs(fBkgTreeVariableDeltaPhi)> 0.85) continue; 
	if (TMath::Abs(fBkgTreeVariableDeltaEta)> 0.75) continue; 
	}
      else       if (EffRegion==2){
	if (TMath::Abs(fBkgTreeVariableDeltaPhi)< 0.85 || TMath::Abs(fBkgTreeVariableDeltaPhi)> 2.05) continue; 
	if (TMath::Abs(fBkgTreeVariableDeltaEta)< 0.75 || TMath::Abs(fBkgTreeVariableDeltaEta)> 1.18) continue; 
      }
    }
*/ 
    for(Int_t m=0; m<nummolt+1; m++){
      if (m < nummolt){
	MoltSel = (fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1]);
      }
      else  if(m==nummolt) MoltSel=kTRUE;
  
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  if(MoltSel){
	    CounterBkgPairsAfterPtMinCutMult[m]++;  
	    if (isParticleTrue)      TrueCounterBkgPairsAfterPtMinCutMult[m]++;  
	    if (!ishhCorr){
	      hMassvsPt_MEbins[m][z][tr]->Fill(fBkgTreeVariableInvMass, fBkgTreeVariablePtV0);
	      if(isParticleTrue) 	    hMassvsPt_MEbins_true[m][z][tr]->Fill(fBkgTreeVariableInvMass, fBkgTreeVariablePtV0);
	    }
	  }
	  for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	    if(MoltSel && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      if (!ishhCorr)	      hInvMassK0Short_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableInvMass);
	      hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	    }
	  }
	}
      }
    }
  }
  }


  cout << " ratio between K0s with CommonParton (CP) and K0s without in the different pt bins " << endl;
  Float_t      CPTrue[numPtV0] ={0};
  Float_t      NOCPTrue[numPtV0]={0};

  for (Int_t v=PtBinMin; v< numPtV0Max; v++){
    for(Int_t b= hSign_PtAssoc_CPTrue->GetXaxis()->FindBin(NPtV0[v]+0.001); b<= hSign_PtAssoc_CPTrue->GetXaxis()->FindBin(NPtV0[v+1]-0.001); b++ ){
      CPTrue[v] +=  hSign_PtAssoc_CPTrue->GetBinContent(b);
      NOCPTrue[v] +=  hSign_PtAssoc_NOCPTrue->GetBinContent(b);
    }
    cout << "v " << v << " " <<  CPTrue[v]/NOCPTrue[v]<< endl;
  }

  cout << "Pt Min delle particelle trigger " << PtTrigMin<< endl;

  cout << "\nsignal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterSignPairsAfterPtMinCut <<" true: " <<   TrueCounterSignPairsAfterPtMinCut << " -> all entries in sign Tree were : " <<   EntriesSign << " " <<(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign << endl;
  cout << "bkg pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterBkgPairsAfterPtMinCut << " true: " <<   TrueCounterBkgPairsAfterPtMinCut << " -> all entries in sign Tree were : " <<   EntriesBkg << " " <<(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg <<  endl;

  for (Int_t m=0; m< nummolt; m++){
    cout << m << endl;  
    cout << "signal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterSignPairsAfterPtMinCutMult[m] << ", " << (Float_t)CounterSignPairsAfterPtMinCutMult[m]/CounterSignPairsAfterPtMinCut<<endl;
    cout << "bkgal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterBkgPairsAfterPtMinCutMult[m] << ", " << (Float_t)CounterBkgPairsAfterPtMinCutMult[m]/CounterBkgPairsAfterPtMinCut<<endl;
  }

  cout << "\nsignal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterSignPairsAfterPtMinCut <<" percentage of INT7 events with at least one selected V0 (calculated assuming 1V0 per event in which there is a V0) "<<(Float_t)CounterSignPairsAfterPtMinCut/TotEvtINT7  <<  endl;

  cout << "entries of V0 selected histogram (all Casc primary true ) "<<  fHistSelectedV0PtTMaxPhi->GetEntries()<< endl;
  cout << "entries of V0 primary histogram (all Casc true ) "<<  fHistPrimaryV0[5]->GetEntries()<< endl;

  cout << "\n Other useful information; " << endl;
  cout << "average pT cascade " <<   hSign_PtAssoc->GetMean()<< " entries: " <<  endl;
  cout << "average pT trigger particles in AC events selected " <<     hSign_PtTrigger->GetMean()<< endl;
  //  cout << "average pT trigger particles in events with NT>0  " << << endl;
  cout << "\n\npartendo dal file " << PathIn << " ho creato il file " << PathOut<< endl;

  cout <<"\n\nINT7     " << " Ev. NT>0    " << "NV0/INT7  " << "<pT,Xi>  " << "<pT,Trig>   " << "SE pairs   " << "ME pairs   " << "Mult. distr                   " << "NV0/ev mult " << endl;  
  cout << std::setprecision(2);
  cout << "    " << TotEvtINT7;
  cout << "    " <<   (Float_t)hMultiplicityBefAll->GetEntries()/TotEvtINT7;
  cout << "    " <<(Float_t)CounterSignPairsAfterPtMinCut/TotEvtINT7;
  cout << std::setprecision(3);
  cout << "    " << hSign_PtAssoc->GetMean();
  cout << "    " <<hSign_PtTrigger->GetMean();
  cout << std::setprecision(2);
  cout << "      " <<(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign;
  cout << "      " <<(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg;
  cout << "      " ;
  for (Int_t m=0; m<nummolt; m++){
    if (m<nummolt-1)    cout << ACcounter[m]/ACcounter[5]<< "-";
    else     cout << ACcounter[m]/ACcounter[5]<< "      ";
  }
  cout << std::setprecision(3);
  for (Int_t m=0; m<nummolt; m++){
    if (m<nummolt-1)    cout <<  hMultvsNumberAssoc_Proj[m]->GetMean() <<"-";
    else    cout <<  hMultvsNumberAssoc_Proj[m]->GetMean() <<endl;
  }

  HistoInfo->SetBinContent(1,(Float_t)hMultiplicityBefAll->GetEntries()/TotEvtINT7);
  HistoInfo->SetBinContent(2,(Float_t)CounterSignPairsAfterPtMinCut/TotEvtINT7);
  HistoInfo->SetBinContent(3,hSign_PtAssoc->GetMean());
  HistoInfo->SetBinContent(4,hSign_PtTrigger->GetMean());
  HistoInfo->SetBinContent(5,(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign);
  HistoInfo->SetBinContent(6,(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg);
  for (Int_t m=0; m< nummolt; m++){
    HistoInfo->SetBinContent(7+m,(Float_t)CounterSignPairsAfterPtMinCutMult[m]/CounterSignPairsAfterPtMinCut);
    HistoInfo->SetBinContent(12+m, hMultvsNumberAssoc_Proj[m]->GetMean());
    HistoInfo->SetBinContent(17+m,	   (float)CounterSignPairsAfterPtMinCutMult[m]/CounterSignPairsAllAssocMult[m]);
  }


  for (Int_t m=0; m< nummolt+1; m++){
    cout << "Ratio between |#Delta#phi|<1 and all the rest for K0s cand. with origin in common with trigger  " << (Float_t)InJet[m]/OutJet[m] << endl;
  }


  cout <<"check if you have to change histogram binning for selected trigger particles (double the number of bins) " << endl;

  //  gStyle->SetOptStat(1011);
  TCanvas *DeltaPhiProj;
  TCanvas *DeltaPhiProjMult[nummolt+1];
  TLegend * legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  TLegend * legendPtBins = new TLegend(0.6, 0.7, 0.9, 0.9);
  
  if (CommonParton){
    DeltaPhiProj= new TCanvas("DeltaPhiProj", "DeltaPhiProj", 1300, 800);
    DeltaPhiProj->Divide(3,2);
    for(Int_t m=0; m<nummolt+1; m++){
      DeltaPhiProj->cd(m+1);
      for(Int_t z=0; z<numzeta; z++){
        for(Int_t tr=0; tr<numPtTrigger; tr++){
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_CPPtInt[m][z][tr]->ProjectionY( "SE_m"+Smolt[m]+"_CP_PtInt_py",0, -1, "E");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->GetXaxis()->SetTitle("#Delta #phi (rad)");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->GetXaxis()->SetLabelSize(0.05);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->GetXaxis()->SetLabelOffset(1);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->SetLineColor(kRed);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->SetMarkerColor(kRed);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->Rebin(2);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr]->Draw("");
	  if (m==0)	  legend->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[m][z][tr], "CP cand.", "pel");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[m][z][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[m][z][tr]->ProjectionY("SE_m"+Smolt[m]+"_CPTrue_PtInt_py",0, -1, "E");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[m][z][tr]->SetLineColor(kMagenta);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[m][z][tr]->SetMarkerColor(kMagenta);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[m][z][tr]->Rebin(2);
          //hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[m][z][tr]->Draw("same");               
	  //	  if (m==0)	  legend->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[m][z][tr], "CPTrue", "pel"); 

	  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[m][z][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->ProjectionY( "SE_m"+Smolt[m]+"_NOCP_PtInt_BulkReg_py",hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->GetXaxis()->FindBin(0.8),hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->GetXaxis()->FindBin(1.1) , "E");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[m][z][tr]->SetLineColor(kAzure+1);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[m][z][tr]->SetMarkerColor(kAzure+1);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[m][z][tr]->Rebin(2);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[m][z][tr]->Draw("same");
	  if (m==0)	  legend->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[m][z][tr], "NOCP cand. bulk", "pel");

          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[m][z][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->ProjectionY( "SE_m"+Smolt[m]+"_NOCP_PtInt_JetReg_py",hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->GetXaxis()->FindBin(-0.7),hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->GetXaxis()->FindBin(0.7) , "E");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[m][z][tr]->SetLineColor(kViolet+1);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[m][z][tr]->SetMarkerColor(kViolet+1);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[m][z][tr]->Rebin(2);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[m][z][tr]->Draw("same");
	  if (m==0)	  legend->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[m][z][tr], "NOCP cand. jet", "pel");

          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[m][z][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[m][z][tr]->ProjectionY( "SE_m"+Smolt[m]+"_NOCPTrue_PtInt_py",0, -1, "E");
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[m][z][tr]->SetLineColor(kBlue+3);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[m][z][tr]->SetMarkerColor(kBlue+3);
          hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[m][z][tr]->Rebin(2);
          // hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[m][z][tr]->Draw("same");  
	  //	  if (m==0)	  legend->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[m][z][tr], "NOCPTrue", "pel");

	  legend->Draw("");


	  DeltaPhiProjMult[m] = new TCanvas(Form("DeltaPhiProj_%i",m), Form("DeltaPhiProj_%i",m), 1300, 800);
          DeltaPhiProjMult[m] ->Divide(4,2);

          for (Int_t v=PtBinMin; v< numPtV0Max; v++){
	    DeltaPhiProjMult[m] ->cd(v+1);
	    hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]->ProjectionY( "SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_CP_py",0, -1, "E");
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->SetLineColor(kRed);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->SetMarkerColor(kRed);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->Rebin(2);
            if (v<=5)  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->Rebin(2);

            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]->ProjectionY("SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_CPTrue_py",0, -1, "E");
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->SetLineColor(kMagenta);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->SetMarkerColor(kMagenta);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->Rebin(2);
            if (v<=5)  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->Rebin(2);

            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->ProjectionY( "SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_NOCP_py",0, -1, "E");
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->SetLineColor(kBlue);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->SetMarkerColor(kBlue);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->Rebin(2);
            if (v<=5)  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->Rebin(2);

            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->ProjectionY( "SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_NOCP_BulkReg_py",hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->GetXaxis()->FindBin(0.8), hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->GetXaxis()->FindBin(1.3), "E");
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->SetLineColor(kAzure+1); //861                   
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->SetMarkerColor(kAzure+1);
            //      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->Scale(10./0.5);                         
 
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->Rebin(2);
            if (v<=5)  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->Rebin(2);

            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->ProjectionY( "SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_NOCP_JetReg_py",hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->GetXaxis()->FindBin(-0.7), hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->GetXaxis()->FindBin(0.7), "E");
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->SetLineColor(kViolet+1); //881                    
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->SetMarkerColor(kViolet+1);
            //      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->Scale(1./1.4);                             
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->Rebin(2);
            if (v<=5)  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->Rebin(2);
	    hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr] = (TH1D*)hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]->ProjectionY( "SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_NOCPTrue_py",0, -1, "E");
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->SetLineColor(kBlue+3);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->SetMarkerColor(kBlue+3);
            hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->Rebin(2);
            if (v<=5)  hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->Rebin(2);

            if (hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->GetMaximum() > hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->GetMaximum()){
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1.3*hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->GetMaximum());
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->Draw("");
              //              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->Draw("same");                          
	      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->Draw("same");
            }
            else {
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1.3*hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->GetMaximum());
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1.3*hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->GetMaximum());
              //              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr]->Draw("");                               
	      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->Draw("");
	      //              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->Draw("same");
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr]->Draw("same");
            }


            if (hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->GetMaximum() > hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->GetMaximum()){
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1.3*hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->GetMaximum());
	      //	      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->Draw("same");                            
	      //	      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->Draw("same");                       
            }
            else {
              hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1.3*hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->GetMaximum());

	      //	      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr]->Draw("same");                     //
	      //	      hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr]->Draw("same");                       

            }

	    hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->Draw("same");
	    hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->Draw("same");

	    if (m==0 && v==0)	  legendPtBins->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[m][z][v][tr], "CPTrue", "pel");
	    if (m==0 && v==0)	  legendPtBins->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[m][z][v][tr], "NOCPTrue", "pel");

	    if (m==0 && v==0)	  legendPtBins->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[m][z][v][tr], "CP cand.", "pel");
	    if (m==0 && v==0)	  legendPtBins->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[m][z][v][tr], "NOCP cand.", "pel");
	    if (m==0 && v==0)	  legendPtBins->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr], "NOCP cand. jet", "pel");
	    if (m==0 && v==0)	  legendPtBins->AddEntry(hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr], "NOCP cand. bulk", "pel");

	    legendPtBins->Draw("");

          }
        }
      }
    }
  }

  if (CommonParton){
    fout->WriteTObject(DeltaPhiProj);
    for (Int_t m=0; m< nummolt+1; m++){
      fout->WriteTObject(DeltaPhiProjMult[m]);
    }

    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	cout << " some infos about deltaphi projections: " << endl;
	cout << " ratio between max value and min value" << endl;
	for (Int_t m=0; m<nummolt+1; m++){
	  cout << "\nn m " << m << endl;
	  for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	    cout << " pt bin " << NPtV0[v] << " - " << NPtV0[v+1] << " NOCP_JetReg " <<     hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->GetMaximum()<< " " <<   hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->GetMinimum()<< endl;
	    cout << " pt bin " << NPtV0[v] << " - " << NPtV0[v+1] << " NOCP_JetReg " <<     hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->GetMaximum()/    hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[m][z][v][tr]->GetMinimum()<< " NOCP_BulkReg " <<     hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->GetMaximum()/    hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[m][z][v][tr]->GetMinimum()<< endl;
	  }
	}
      }
    }
  }

  TH1D *  hSign_PtTriggerPtAssoc_Proj[numPtV0];
  cout << "average pT of trigger particles associated to associated particles in pt bins: " << endl;
  for(Int_t v=PtBinMin; v<numPtV0Max; v++){
    hSign_PtTriggerPtAssoc_Proj[v]= (TH1D*)	  hSign_PtTriggerPtAssoc->ProjectionY("hSign_PtTriggerPtAssoc_Proj_v"+SPtV0[v], hSign_PtTriggerPtAssoc->GetXaxis()->FindBin(NPtV0[v]+0.0001) ,  hSign_PtTriggerPtAssoc->GetXaxis()->FindBin(NPtV0[v+1]-0.0001),"E");
    cout << SPtV0[v] << " " <<  hSign_PtTriggerPtAssoc_Proj[v]->GetMean() << endl;
  }
  cout << " writing " << endl;
  fout->Write();
  cout << "cut values " << endl;
  cout << "CosinePAngle  :" << CosinePAngle    <<endl;
  cout << "DCANegToPV    :"  <<DCANegToPV      <<endl;
  cout << "DCAPosToPV    :"  <<DCAPosToPV      <<endl;
  cout << "DCAV0ToPV     :"  <<DCAV0ToPV       <<endl;
  cout << "ctau          :"  <<ctau            <<endl;
  cout << " number of entries " << CounterSignPairsAfterPtMinCut<< endl;

}

