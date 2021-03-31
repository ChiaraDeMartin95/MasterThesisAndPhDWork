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

Double_t SetEfficiencyError(Int_t k, Int_t n){
  return sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
}

void readTreePLChiarahK0s_second(Bool_t ishhCorr=0, Int_t type=0, Bool_t SkipAssoc=1, Int_t israp=0, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Float_t PtTrigMin=3,Float_t PtTrigMinFit=3, Int_t sysTrigger=0, Int_t sysV0=0,Int_t syst=0,bool isMC = 0, Bool_t isEfficiency=0,TString year0="2016", TString year=/*"AllhK0sHM_RedNo16k"/*"1617GP_hK0s"/*"2018f1_extra_15runs"/*"2018f1_extra_RlabelBis_15runs_hK0s_Hybrid"/*"2018f1_extra_MylabelBis_15runs_hK0s_Hybrid"*/"2016k_HM_hK0s"/*"2018f1_extra_15runs_NohDaughtersofK0s_hK0s_Hybrid"/*"2018g4_extra_EtaEff_hK0s"/*"2018f1_extra_5runs_label_Hybrid_hK0s"/"1617GP_hK0s_Hybrid"/*"2018f1_Reco_hK0s"/"2018f1_extra_hK0s"/*"2018f1_extra_hK0s"/*"LHC17_hK0s"/"1617_hK0s"/*"2016kehjl_hK0s"/"2018f1_extra_hK0s_30runs_Hybrid"/"2016k_hK0s"/"2018f1_extra_hK0s"/"2018g4_extra_EtaEff_hK0s"*/,  TString Path1 ="",  Double_t ptjmax =15, Double_t nsigmamax=10, Bool_t isSigma=kFALSE, Int_t PtBinning=1, Bool_t IsTrueParticle=0, Bool_t IsPtTrigMultDep =0, Int_t TriggerPDGChoice=0, Bool_t isEtaEff=1, TString yearMC/*path of file where efficiency is stored*/ = "2019h11c_extra_HM_hK0s"/*"2018g4_extra_EtaEff_hK0s"*/, Bool_t isEta05=0, Bool_t isPrimaryTrigger=0, Bool_t isNewInputPath=1, Bool_t isHM=1, Int_t MultBinning=1){

  if (!isMC) isPrimaryTrigger=0;
  if (isMC && isEfficiency) isPrimaryTrigger=0;
  //isPrimaryTrigger is used when dealing with hybrid!

  if (isMC && !isEfficiency) {
    isEtaEff=0; //no efficiency correction needed in this case
    IsTrueParticle=0; //if I deal with MCTruth/MCHybrid, all associated particles are true by default
  }
  if (!isMC) {
    IsTrueParticle=0;
  }

  cout << year.Index("EtaEff")<< endl;
  //  if (year.find("Hybrid")!=string::npos && isEfficiency) {cout << " when analyzing hybrid isEfficiency should be zero " << endl; return;}
  if (!(year.Index("Hybrid")==-1) && isEfficiency) {cout << " when analyzing hybrid isEfficiency should be zero " << endl; return;}
  //***************************************************************************
  //******** What particles to choose as trigger whan analysing MC truth? *****

  // Since in data we are able to identify basically only pi, K, p, mu, e, also when analyzing the MC truth we want to take only these particles as trigger particles. And we want to exclude all the others already WHEN CHOOSING the trigger particle, so at the level of the TASK! 
  //Here some offline selections are implemented as well, but they are useless if the selection is done at the level of the task 
  //Check if fHistTriggerPDGCode and fHistTriggerPDGCodeAfterSel are the same. If it is so, the selections are correctly implemtned in the task

  //***************************************************************************

  TString sTriggerPDGChoice[3] = {"", "_IsOnlypiKpemu", "_IsNotSigmaOnly"};

  //isEtaEff ==1 : efficiency correction done here, taking into account the eta dependence of the efficiency

  //isMeanFixedPDG and isBkgParab are characteristics of the fit to the inv mass distributions 
  cout << isMC << endl;
  cout << " Pt Trigg Min è = " << PtTrigMin << endl;

  if (israp>1) return;
  if (israp==1 && ishhCorr) {cout << "in hh correlation associated hadrons are not identified" << endl; return;}
  if (type>3) {cout << "type value not allowed" << endl; return;}
  //lista degli effetti  sistematici studiati in questa macro
  //sys=1 nsigmamin=5 (def:4)
  //sys=2 sigmacentral =4 (def:3)
  if (sysV0>6) return;

  if (syst!=0 && syst!=1 && syst!=2) {
    cout << "syst should be changed " << endl;
    return;
  }
  if((sysV0!=0||sysTrigger!=0) && syst!=0){
    cout << "syst conflicting " << endl;
    return;
  }

  TF1 * lineat1 = new TF1("pol0", "pol0", 0, 30);
  lineat1->FixParameter(0,1);
  lineat1->SetLineColor(kBlack);

  Int_t PtBinMin=0;
  if (!ishhCorr && type!=0) PtBinMin=1; //for associated particles different from hadrons I do not start from 0            

  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  const Float_t ctauK0s=2.6844;

  TString file;

  const Int_t numtipo=4;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda", "LambdaAntiLambda"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  if (isEta05) Srap[0] = "_Eta0.5";
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
  Float_t ctauCasc[numtipo] = {2.6844,7.89, 7.89, 7.89};
  Float_t PDGCode[numtipo-1] = {310, 3122, -3122};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  const Int_t numPtTrigger=1;

  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString SBismolt[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  TString Szeta[numzeta]={""};

  Int_t numMultBins=100;
  Float_t UpperLimitMult = 100;
  if (isHM) UpperLimitMult = 0.1;

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
    SBismolt[0] = "0-0.001 %";
    SBismolt[1] = "0.001-0.005 %";
    SBismolt[2] = "0.005-0.01 %";
    SBismolt[3] = "0.01-0.05 %";
    SBismolt[4] = "0.05-0.1 %";
    SBismolt[5] = "0-0.1 %";
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
      SBismolt[0] = "0-0a %";
      SBismolt[1] = "0-0b %";
      SBismolt[2] = "0-0.01 %";
      SBismolt[3] = "0.01-0.05 %";
      SBismolt[4] = "0.05-0.1 %";
      SBismolt[5] = "0-0.1 %";
    }
  }

  TString PathIn="FinalOutput/AnalysisResults";
  TString PathOut="FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";
  TString PathInMass= "FinalOutput/DATA" + year0 + "/invmass_distribution_thesis/invmass_distribution";
  TString PathInMassDef;
  PathIn+=year;
  PathOut+=year;  
  // PathInMass+=year;
  
  if (ishhCorr){
    PathIn+="_hhCorr";
  }
 
  if(isMC && isEfficiency){ 
    PathIn+="_MCEff";
    PathOut+="_MCEff";
    PathInMass+="_MCEff";
  }
  if(isMC && !isEfficiency){
    PathIn+="_MCTruth";
    PathOut+="_MCTruth";
  }
 
  //  PathIn+=Path1; //change
  if(PtBinning)  PathInMass+=Form("_PtBinning%i",PtBinning);
  PathInMass+=Path1;
  PathIn+=".root";


  if(PtBinning)  PathOut+=Form("_PtBinning%i",PtBinning);
  PathOut+=Path1; //+"Bis";
  PathOut+="_"; 
  if (!ishhCorr) PathOut +=tipo[type];
  PathOut +=Srap[israp];
  PathOut +=SSkipAssoc[SkipAssoc];
  if (!ishhCorr) PathOut +=Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f",sysTrigger, sysV0, syst, PtTrigMin); 
  if (ishhCorr&& (!isMC ||(isMC && isEfficiency))) PathOut +=Form("_hhCorr_SysT%i_SysV0%i_Sys%i_PtMin%.1f",sysTrigger, sysV0, syst, PtTrigMin);
  if (IsTrueParticle) PathOut+="_IsParticleTrue";
  if (IsPtTrigMultDep) PathOut += "_IsPtTrigMultDep";
  PathOut+= sTriggerPDGChoice[TriggerPDGChoice];
  //  PathOut += "_EtaStudy";
  //  PathOut+="_TryisEta05";
  if (isEtaEff) PathOut+="_isEtaEff";
  //  PathOut+="_isEtaEff";
  //  PathOut+= "_DCAz0.5";
  if (isPrimaryTrigger)   PathOut+="_isPrimaryTrigger";
  if (MultBinning!=0) PathOut += Form("_MultBinning%i", MultBinning);
  //  PathOut+="_Try";
  //  PathOut += "_AllEta";
  TString pathoutpdf = PathOut;
  PathOut+= ".root";

  cout << "file di input " << PathIn << endl;
  TFile *fin = new TFile(PathIn);
  if (!fin) {cout << PathIn << " not available " << endl; return;} 

  //input file where efficiency is found
  TString PathInEfficiency = "FinalOutput/DATA2016/Efficiency/Efficiency" +yearMC;
  //  if (PtBinning>0)    PathInEfficiency+=Form("_PtBinning%i",PtBinning);
  PathInEfficiency+=Form("_PtBinning%i",1);
  PathInEfficiency+=Path1;// change                                                                   
  if(type>=0){
    PathInEfficiency +="_";
    if (!ishhCorr)             PathInEfficiency +=tipo[type];
    PathInEfficiency +=Srap[israp];
    PathInEfficiency +=SSkipAssoc[SkipAssoc];
    if (ishhCorr)     PathInEfficiency +="_hhCorr";
  }

  PathInEfficiency+= Form("_SysT%i_SysV0%i_PtMin%.1f", sysTrigger, sysV0, PtTrigMin);
  /*
  if (IsPtTrigMultDep)    PathInEfficiency+= "_IsPtTrigMultDep";
  if (IsEfficiencyMassSel)    PathInEfficiency+= "_isEfficiencyMassSel";
  */
  if (MultBinning!=0) PathInEfficiency += Form("_MultBinning%i", MultBinning);
  PathInEfficiency+= ".root";
  
  TFile *fileinEfficiency = new TFile (PathInEfficiency, "");

  TH2F * fHistEfficiencyV0PtEta[nummolt+1];
  TH1F * fHistEfficiencyV0PtPtBins[nummolt+1];
  if (isEtaEff){
    if (!fileinEfficiency) {cout << " input file with efficiency not found " << endl; return;}
    for(Int_t molt=0; molt<nummolt+1; molt++){
      if (MultBinning==1 && isHM){
	if (molt==0 || molt==1) continue;
      }
      fHistEfficiencyV0PtEta[molt] = (TH2F*) fileinEfficiency->Get("fHistV0EfficiencyPtV0EtaV0PtBins_"+ Smolt[molt]);
      cout << Smolt[molt] << endl;
      if (!fHistEfficiencyV0PtEta[molt]) {cout << "histogram 2D V0 efficiency pt vs eta " << endl; return;}
      fHistEfficiencyV0PtPtBins[molt] = (TH1F*) fileinEfficiency->Get("fHistV0EfficiencyPtBins_"+ Smolt[molt]);
      if (!fHistEfficiencyV0PtPtBins[molt]) {cout << "histogram 1D V0 efficiency pt " << endl; return;}
    }
  }

  Float_t   EffRelErrSign[nummolt+1][numPtV0]={0};
  Float_t   EffRelErrBkg[nummolt+1][numPtV0]={0};
  TFile *fileMassSigma;
  TString dirinputtype[4] = {"", "Lambda", "Lambda", "Lambda"};
  TDirectoryFile *d;
  if (isNewInputPath){
    if (!(year.Index("Hybrid")==-1))   {
      //      d = (TDirectoryFile*)fin->Get("MyTask_MCHybrid_PtTrigMin3.0_PtTrigMax15.0");
      d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin3.0_PtTrigMax15.0");
    }
    else if (isMC && !isEfficiency){
      d = (TDirectoryFile*)fin->Get("MyTask_MCTruth_PtTrigMin3.0_PtTrigMax30.0");
    }
    else {
      if (year=="2016k_HM_hK0s") d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin3.0_PtTrigMax30.0");
      else d = (TDirectoryFile*)fin->Get("MyTask_PtTrigMin3.0_PtTrigMax15.0");
    }
  } else {
    d = (TDirectoryFile*)fin->Get("MyTask"+dirinputtype[type]);
  }  
  if (!d) {cout << "directory not available " << endl; return;}

  TString NameContainer= "";
  if (isNewInputPath){
    if (!(year.Index("Hybrid")==-1))  NameContainer = "_hK0s_Task_Hybrid";
    else if (year=="2016k_HM_hK0s") NameContainer = "_hK0s_Task_suffix";
    else NameContainer = "_hK0s_Task_";
  }

  TDirectoryFile *d1 = (TDirectoryFile*)d->Get("MyOutputContainer" + NameContainer);
  if (!d1) {cout << "container not available " << endl; return;}


  TTree *tSign = (TTree *)d->Get("fSignalTree");
  TTree *tBkg  = (TTree *)d->Get("fBkgTree");
  if (!tSign) {cout << "Sign Tree is not there! " << endl; return; }
  if (!tBkg) {cout << "Bkg Tree is not there! " << endl; return; }

  
  Int_t Color[nummolt+1]={1,2,8,4,6,868};
  Float_t LimInfMass[numtipo][nummolt+1][numPtV0]={0};
  Float_t LimSupMass[numtipo][nummolt+1][numPtV0]={0}; 

  //identità particelle di trigger nel MC                                                                 
  TH1F* fHistTriggerPDGCode = new TH1F("fHistTriggerPDGCode", "fHistTriggerPDGCode", 11,0, 11);
  fHistTriggerPDGCode->SetTitle("PDGcode trigger particle");
  Int_t PDGCodeParticles[10] = {211,2212 , 321, 11, 13, 15, 3312, 3334, 3112, 3222}; //sigma- and sigma + respectively                                                                                                
  TString ParticlePDG[11] = {"#pi", "p", "K", "e", "#mu", "#tau", "Xi", "#Omega", "#Sigma^{+}", "#Sigma^{-}", ""};
  for (Int_t i =0; i<11; i++){
    fHistTriggerPDGCode->GetXaxis()->SetBinLabel(i+1,ParticlePDG[i]);
    if (i==10)     fHistTriggerPDGCode->GetXaxis()->SetBinLabel(i+1,"other");
  }

  TH1F* fHistTriggerPDGCodeAfterSel = (TH1F*) fHistTriggerPDGCode->Clone("fHistTriggerPDGCodeAfterSel");

  const Int_t  numPtTrInt=12;
  Float_t NPtTrInt[numPtTrInt+1] ={3,3.5,4,4.5,5,5.5,6,7,8,9,10,12,15};

  const Int_t numQAhisto=8;
  TH1F * fHistQA[numQAhisto];
  TH1F * fHistPtTriggervsPtAssoc[nummolt+1];
  TH1F * fHistNonPrimaryTrigger[nummolt+1];
  TH1F * fHistEfficiencyReduction[nummolt+1];
  TH2F * fHistDCATrigger = new TH2F("fHistDCATrigger", "fHistDCATrigger", 100, -0.02, 0.02, 400, -0.5, 0.5);

  TCanvas* canvasQA[numQAhisto];
  TString TitleQAhisto[numQAhisto] = { "Average pT trigger in AC events", "average pT trigger in AC events (each trigger counted only once", "mult distribution of AC events (each trigger counted only once)", "NV0/AC event", "fraction of V0 with pt< pt,trig", "average pT trigger events NT>0", "fraction of non primary trigger vs p_{T}^{K^{0}_{S}}", "effect of DCAPos/Neg>0.1 on K0s efficiency" };
  for (Int_t i=0; i<numQAhisto; i++){
    canvasQA[i]= new TCanvas(Form("canvasQA%i",i),TitleQAhisto[i], 800, 500);
    if (i==0 || i==1 || i==5) canvasQA[i]->Divide(2,2);
    if (i!=2){
      fHistQA[i]= new TH1F (Form("fHistQA%i", i), TitleQAhisto[i], nummolt,Nmolt ); //multiplicity on the x axis               
      fHistQA[i]->GetXaxis()->SetTitle("Multiplicity class");
    }
    else {
      fHistQA[i]= new TH1F (Form("fHistQA%i", i), TitleQAhisto[i], numMultBins,0,UpperLimitMult);
      fHistQA[i]->GetXaxis()->SetTitle("Multiplicity class");
    }
    //    if (i==numQAhisto-1) fHistQA[i]->SetTitle("fraction of V0 with pT< pT,Trig");
  }

  cout << d->GetName() << endl;
  cout << d1->GetName() << endl;
  TH2F * fHistPtMaxvsMultBefAll = (TH2F*) d1->FindObject("fHistPtMaxvsMultBefAll");
  if (!d1->FindObject("fHistPtMaxvsMultBefAll")) return;
  TH1F * fHistPtMaxvsMultBefAllProj[nummolt+1];
  TH1F * fHistPtMaxvsMultBefAllProjRatio[nummolt+1];
  fHistPtMaxvsMultBefAll->GetXaxis()->SetRangeUser(PtTrigMin, ptjmax);
  cout << "hola " << endl;

  TH1F* fHistEventMult=(TH1F*)  d1->FindObject("fHistEventMult");
  if (!fHistEventMult) cout << "no info about total number of INT7 events analyzed" << endl;
  cout <<" bin label: " << fHistEventMult->GetXaxis()->GetBinLabel(7) << endl;
  Double_t TotEvtINT7 = fHistEventMult->GetBinContent(7);

  for (Int_t m=0; m<nummolt+1; m++){
    cout << "hola " << m <<  endl;
    if (m!=nummolt)    fHistPtMaxvsMultBefAllProj[m] = (TH1F*) fHistPtMaxvsMultBefAll->ProjectionX(Form("fHistPtMaxvsMultBefAllProj_m%i",m),   fHistPtMaxvsMultBefAll->GetYaxis()->FindBin(Nmolt[m]+0.0001),  fHistPtMaxvsMultBefAll->GetYaxis()->FindBin(Nmolt[m+1]-0.0001));
    else     fHistPtMaxvsMultBefAllProj[m] = (TH1F*) fHistPtMaxvsMultBefAll->ProjectionX(Form("fHistPtMaxvsMultBefAllProj_m%i",m), 0, UpperLimitMult);
    fHistQA[5]->SetBinContent(m+1, fHistPtMaxvsMultBefAllProj[m]->GetMean());
    fHistQA[5]->SetBinError(m+1, fHistPtMaxvsMultBefAllProj[m]->GetMeanError());
  }

  cout << "hola " << endl;
  canvasQA[5] ->cd(1);
  fHistQA[5]->Draw("");

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};
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

  //  TString SPtTrigger[numPtTrigger]={"2-10"};
  Float_t PtTrigMultDep[nummolt+1] = {3.0, 3.05, 3.05, 3.1, 3.1, 3};
  Double_t NPtTrigger[numPtTrigger+1]={PtTrigMin,ptjmax};

  for (Int_t m=0; m<nummolt+1; m++){
    if (!IsPtTrigMultDep) PtTrigMultDep[m] = PtTrigMin;
    fHistPtTriggervsPtAssoc[m]=new TH1F (Form("fHistPtTriggervsPtAssoc%i",m), Form("fHistPtTriggervsPtAssoc%i",m), numPtV0,NPtV0);
    fHistNonPrimaryTrigger[m]=new TH1F (Form("fHistNonPrimaryTrigger%i",m), Form("fHistNonPrimaryTrigger%i",m), numPtV0,NPtV0);
    fHistEfficiencyReduction[m]=new TH1F (Form("fHistEfficiencyReduction%i",m), Form("fHistEfficiencyReduction%i",m), numPtV0,NPtV0);
  }

  Double_t sigma[numtipo][nummolt+1][numPtV0];
  Double_t mass[numtipo][nummolt+1][numPtV0];
  Double_t nsigmamin[numtipo][nummolt+1][numPtV0];
  Double_t sigmacentral[numtipo][nummolt+1][numPtV0];
  Double_t meansigma;
  Double_t meanmass;
  TH1F* histoSigma;
  TH1F* histoMean;
  TH1F* histo_ULsideB;
  TH1F* histo_LLsideB;
  TH1F* histo_NSigmaPeak;
  TH1F* histo_NSigmasideB;

  if (!ishhCorr){ 
    if(isMC==0 || (isMC==1 && isEfficiency==1)){
      for(Int_t m=0; m<nummolt+1; m++){
	if (MultBinning==1 && isHM){
	  if (m==0 || m==1) continue;
	}
	//      if (m!= nummolt) continue;
	//      PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, syst, PtTrigMin);
	PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab];
	if(IsPtTrigMultDep) PathInMassDef  +="_IsPtTrigMultDep" ;
	if (MultBinning!=0) PathIn += Form("_MultBinning%i", MultBinning);
	PathInMassDef+= Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, syst,PtTrigMinFit);
	fileMassSigma= new TFile(PathInMassDef);
	histoSigma=(TH1F*)fileMassSigma->Get("histo_sigma");
	histoMean=(TH1F*)fileMassSigma->Get("histo_mean");
	histo_ULsideB=(TH1F*)fileMassSigma->Get("histo_ULsideB");
	histo_LLsideB=(TH1F*)fileMassSigma->Get("histo_LLsideB");
	histo_NSigmasideB=(TH1F*)fileMassSigma->Get("histo_NSigmasideB");
	histo_NSigmaPeak=(TH1F*)fileMassSigma->Get("histo_NSigmaPeak");
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  sigmacentral[type][m][v]=      histo_NSigmaPeak->GetBinContent(v+1);
	  nsigmamin[type][m][v]=      histo_NSigmasideB->GetBinContent(v+1);
	  mass[type][m][v]=histoMean->GetBinContent(v+1);
	  sigma[type][m][v]=histoSigma->GetBinContent(v+1);
	  LimSupMass[type][m][v]=histo_ULsideB->GetBinContent(v+1);
	  LimInfMass[type][m][v]=histo_LLsideB->GetBinContent(v+1);
	  cout <<"mult interval " <<  m << " PtV0 interval " << v << " mean " << mass[type][m][v] << " sigma "<< sigma[type][m][v] << endl;
	}
      }
    }
  }
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

  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsEffw[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsEffwRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsEffwErrors[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hSign_PtTrigger[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerRatio[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerCountedOnce[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerCountedOnceBis[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerCountedOnceRatio[nummolt+1][numzeta];
  TH1D *hSign_PtV0[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerPtV0bins[nummolt+1][numzeta][numPtV0];
  Float_t CounterTriggerCountedOnce[nummolt+1][numzeta]={0};
  Float_t CounterACPairs[nummolt+1][numzeta]={0};
  Float_t CounterACPairsPtTr3[nummolt+1][numzeta]={0};
  Float_t CounterV0NotSkipped[nummolt+1]={0};
  Float_t CounterV0NotSkippedPtTr3[nummolt+1]={0};

  Float_t CounterNotPrimaryTrigger[nummolt+1][numPtV0] ={0};
  Float_t CounterAllTrigger[nummolt+1][numPtV0] ={0};

  Float_t CounterPairsBeforeSel[nummolt+1][numPtV0] ={0};
  Float_t CounterPairsAfterSel1[nummolt+1][numPtV0] ={0};

  Float_t CounterPairsPtK0sGr3[nummolt+1] ={0};
  Float_t CounterAllPairs[nummolt+1] ={0};

  for(Int_t m=0; m<nummolt+1; m++){
    CounterV0NotSkipped[m]=0;
    CounterV0NotSkippedPtTr3[m]=0;
    for(Int_t z=0; z<numzeta; z++){
      CounterACPairs[m][z]=0;
      CounterACPairsPtTr3[m][z]=0;
      CounterTriggerCountedOnce[m][z]=0;
    }
  }

  const   Int_t numDeltaEta=4;
  TString SDeltaEta[numDeltaEta]={"_|DeltaEta|<1.6","_|DeltaEta|<1.2","_|DeltaEta|<0.75", "_0.75<|DeltaEta|<1.2" };
  Float_t DeltaEtaLimitLow[numDeltaEta]={-1.6, -1.2, -0.75, 0.75};
  Float_t DeltaEtaLimitUp[numDeltaEta]={1.6, 1.2, 0.75, 1.2};
  TH1D *hSign_EtaV0[nummolt+1][numDeltaEta];
  TH1D *hSign_EtaV0MidRap[nummolt+1][numDeltaEta];
  TH1D *hSign_RapV0[nummolt+1][numDeltaEta];
  TH2D *hSign_DeltaEtaEtaV0[nummolt+1];
  TH2D *hSign_DeltaEtaEtaTrigger[nummolt+1];

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==1 && isHM){
      if (m==0 || m==1) continue;
    }
    hSign_DeltaEtaEtaV0[m]=new TH2D("hSign_DeltaEtaEtaV0_"+Smolt[m], "hSign_DeltaEtaEtaV0_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    hSign_DeltaEtaEtaTrigger[m]=new TH2D("hSign_DeltaEtaEtaTrigger_"+Smolt[m], "hSign_DeltaEtaEtaTrigger_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    for (Int_t DeltaEta=0; DeltaEta< numDeltaEta; DeltaEta++){
      hSign_EtaV0[m][DeltaEta]=new TH1D("hSign_EtaV0_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_EtaV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hSign_RapV0[m][DeltaEta]=new TH1D("hSign_RapV0_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_RapV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hSign_EtaV0MidRap[m][DeltaEta]=new TH1D("hSign_EtaV0_MidRap_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_EtaV0_MidRap_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
    }
  }

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==1 && isHM){
      if (m==0 || m==1) continue;
    }

    for(Int_t z=0; z<numzeta; z++){

      hSign_PtTrigger[m][z]=new TH1D("hSign_PtTrigger"+Smolt[m], "hSign_PtTrigger"+Smolt[m], 300,0,30);
      hSign_PtTrigger[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTrigger[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTrigger[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTrigger[m][z]->GetYaxis()->SetLabelSize(0.05);
      hSign_PtTrigger[m][z]->SetLineColor(Color[m]);
      hSign_PtTrigger[m][z]->SetMarkerColor(Color[m]);
      hSign_PtTrigger[m][z]->SetMarkerStyle(33);

      hSign_PtTriggerCountedOnceBis[m][z]=new TH1D("hSign_PtTriggerCountedOnceBis"+Smolt[m], "hSign_PtTriggerCountedOnce"+Smolt[m], numPtTrInt, NPtTrInt);
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnceBis[m][z]->GetYaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnceBis[m][z]->SetLineColor(Color[m]);
      hSign_PtTriggerCountedOnceBis[m][z]->SetMarkerColor(Color[m]);
      hSign_PtTriggerCountedOnceBis[m][z]->SetMarkerStyle(33);
      hSign_PtTriggerCountedOnce[m][z]=new TH1D("hSign_PtTriggerCountedOnce"+Smolt[m], "hSign_PtTriggerCountedOnce"+Smolt[m], 600, 0, 30);
      hSign_PtTriggerCountedOnce[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTriggerCountedOnce[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTriggerCountedOnce[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnce[m][z]->GetYaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnce[m][z]->SetLineColor(Color[m]);
      hSign_PtTriggerCountedOnce[m][z]->SetMarkerColor(Color[m]);
      hSign_PtTriggerCountedOnce[m][z]->SetMarkerStyle(33);

      hSign_PtV0[m][z]=new TH1D("hSign_PtV0"+Smolt[m], "hSign_PtV0"+Smolt[m], 300,0,30);
      hSign_PtV0[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtV0[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtV0[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtV0[m][z]->GetYaxis()->SetLabelSize(0.05);

      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  hSign_PtTriggerPtV0bins[m][z][v]=new TH1D("hSign_PtTriggerPtV0bins"+Smolt[m]+"_v"+SPtV0[v], "hSign_PtTriggerPtV0bins"+Smolt[m]+"_v"+SPtV0[v], 300, 0, 30);
	  nameSE[m][z][v][tr]="SE_";
	  namemassSE[m][z][v][tr]="InvMassSE_";
	  nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr], nameSE[m][z][v][tr],   56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_SEbinsRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_RelErr");

	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr]+ "_Effw", nameSE[m][z][v][tr]+ " eff corr",  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_SEbinsEffwRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_EffwRelErr");

	  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_EffwErrors");
	  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->SetTitle(nameSE[m][z][v][tr]+ " relative error of efficiency");

	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr]+"_SB", nameSE[m][z][v][tr]+"_SB",   56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	}
      }
    }
  }

  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbinsRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbinsEffw[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbinsEffwRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbinsEffwErrors[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hBkg_PtTrigger[nummolt+1][numzeta];
  TH1D *hBkg_PtV0[nummolt+1][numzeta];

  TH1D *hBkg_EtaV0[nummolt+1][numDeltaEta];
  TH1D *hBkg_EtaV0MidRap[nummolt+1][numDeltaEta];
  TH1D *hBkg_RapV0[nummolt+1][numDeltaEta];
  TH2D *hBkg_DeltaEtaEtaV0[nummolt+1];
  TH2D *hBkg_DeltaEtaEtaTrigger[nummolt+1];


  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==1 && isHM){
      if (m==0 || m==1) continue;
    }
    hBkg_DeltaEtaEtaV0[m]=new TH2D("hBkg_DeltaEtaEtaV0_"+Smolt[m], "hBkg_DeltaEtaEtaV0_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    hBkg_DeltaEtaEtaTrigger[m]=new TH2D("hBkg_DeltaEtaEtaTrigger_"+Smolt[m], "hBkg_DeltaEtaEtaTrigger_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    for (Int_t DeltaEta=0; DeltaEta< numDeltaEta; DeltaEta++){
      hBkg_EtaV0[m][DeltaEta]=new TH1D("hBkg_EtaV0_"+Smolt[m]+Form("_%i",DeltaEta), "hBkg_EtaV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hBkg_RapV0[m][DeltaEta]=new TH1D("hBkg_RapV0_"+Smolt[m]+Form("_%i",DeltaEta), "hBkg_RapV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hBkg_EtaV0MidRap[m][DeltaEta]=new TH1D("hBkg_EtaV0_MidRap_"+Smolt[m]+Form("_%i",DeltaEta), "hBkg_EtaV0_MidRap_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
    }
  }

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==1 && isHM){
      if (m==0 || m==1) continue;
    }
    for(Int_t z=0; z<numzeta; z++){
      hBkg_PtTrigger[m][z]=new TH1D("hBkg_PtTrigger"+Smolt[m], "hBkg_PtTrigger"+Smolt[m], 300,0,30);
      hBkg_PtV0[m][z]=new TH1D("hBkg_PtV0"+Smolt[m], "hBkg_PtV0"+Smolt[m], 300,0,30);
      hBkg_PtTrigger[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hBkg_PtTrigger[m][z]->GetXaxis()->SetTitleSize(0.05);
      hBkg_PtTrigger[m][z]->GetXaxis()->SetLabelSize(0.05);
      hBkg_PtTrigger[m][z]->GetYaxis()->SetLabelSize(0.05);
      hBkg_PtV0[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hBkg_PtV0[m][z]->GetXaxis()->SetTitleSize(0.05);
      hBkg_PtV0[m][z]->GetXaxis()->SetLabelSize(0.05);
      hBkg_PtV0[m][z]->GetYaxis()->SetLabelSize(0.05);

      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  nameME[m][z][v][tr]="ME_";
	  namemassME[m][z][v][tr]="InvMassME_";
	  nameME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]= new TH2D(nameME[m][z][v][tr], nameME[m][z][v][tr],  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_MEbinsRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_RelErr");

	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]= new TH2D(nameME[m][z][v][tr]+ "_Effw", nameME[m][z][v][tr]+ " eff corr",  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_MEbinsEffwRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_EffwRelErr");

	  hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_EffwErrors");
	  hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->SetTitle(nameME[m][z][v][tr]+ " relative error of efficiency");

	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]= new TH2D(nameME[m][z][v][tr]+"_SB", nameME[m][z][v][tr]+"_SB",  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	}
      }
    }
  }

  EntriesSign =  tSign->GetEntries();
  EntriesBkg  =  tBkg ->GetEntries();
     
  Bool_t BoolVar=kFALSE;
  Bool_t BoolMC=kFALSE;
  Bool_t MassLimit=kFALSE;
  Float_t     fSignTreeVariableInvMass= 0;
  Bool_t isParticleTrue=kFALSE;
  Int_t fSignTreeVariablePAPAssoc=0;
  Float_t  fSignTreeVariablePtTriggerTemp=0;
  Double_t effSign=0;
  Double_t sigmaEffSign=0;

  cout << "\n\n I will process " << EntriesSign << " entries for theSE correlation " << endl;
  for(Int_t k = 0; k<EntriesSign; k++){
    if (k>10000000) continue;
    tSign->GetEntry(k);
    //  for(Int_t k = 0; k<1000000; k++){
    for (Int_t l=0; l<10000; l++){
      //      if (k ==100000*l) cout << "k = " << k << " over a total of " << EntriesSign << endl;
      if (k ==100000*l) {
	Float_t kFloat= k;
	//	cout << " etav0:" << fSignTreeVariableEtaV0 <<  " etatrigger: " << fSignTreeVariableEtaTrigger << " deltaeta: " << fSignTreeVariableDeltaEta << "  new: " << fSignTreeVariableEtaV0-fSignTreeVariableEtaTrigger <<  endl;
	//	cout << " phiv0:" << fSignTreeVariablePhiV0 <<  " phitrigger: " << fSignTreeVariablePhiTrigger << " deltaphi: " << fSignTreeVariableDeltaPhi << "  new: " << fSignTreeVariablePhiV0-fSignTreeVariablePhiTrigger <<  endl;
	cout << "PtTrigMin = " << PtTrigMin << " sysV0 " << sysV0 << " SE, processing..." << kFloat/EntriesSign << endl;
      }
    }


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

    //rapidity selection
    if (israp==0){    
      if (isEta05) {
	if (TMath::Abs(fSignTreeVariableEtaV0)>0.5)continue;
      }
      else{
	if (TMath::Abs(fSignTreeVariableEtaV0)>0.8)continue;
      }
    }
    else if (israp==1){
      if (TMath::Abs(fSignTreeVariableRapK0Short)>0.5)continue;
    }

    //definition of true particle
    if (type<=2) isParticleTrue= (fSignTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==3) isParticleTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[1]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[1])  );
  

    if (IsTrueParticle && isMC && isEfficiency==1){
      if (!isParticleTrue) continue;
      if (fSignTreeVariableisPrimaryV0!=1) continue;
    } 


    //************cuts on pT trigger min*********************************
    if (!IsPtTrigMultDep){
      if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMin) continue;
    }
    else {
      Bool_t IsSelectionPtMinNotPassed=0;
      for(Int_t m=0; m<nummolt; m++){
	if (MultBinning==1 && isHM){
	  if (m==0 || m==1) continue;
	}
	if ( fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]){
	  if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMultDep[m]) IsSelectionPtMinNotPassed=1;
	}
      }
      if (IsSelectionPtMinNotPassed) continue;
    }

    if(TMath::Abs(fSignTreeVariablePtTrigger)>ptjmax) continue;

    BoolVar = kFALSE;

    if(TMath::Abs(fSignTreeVariableDCAz)>0.04) continue;

    if(isMC==0 || (isMC==1 && isEfficiency==1)){
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

      //******************* some other cuts for sys studies**************************
      if (!ishhCorr){
	if (sysV0==0){
	
	  //the values are valid for K0s
	  if (type==0){
	    if(TMath::Abs((fSignTreeVariableInvMassLambda - massLambda))<= 0.005) continue;
	    if(TMath::Abs((fSignTreeVariableInvMassAntiLambda - massLambda))<= 0.005) continue;
	  }
	  if(fSignTreeVariableV0CosineOfPointingAngle<0.995)            continue;
	  if(fSignTreeVariableDcaNegToPrimVertex < 0.06)      continue;
	  if(fSignTreeVariableDcaPosToPrimVertex < 0.06)      continue;
	  if(fSignTreeVariableDcaV0ToPrimVertex > 0.5)                               continue;
	  if(fSignTreeVariablectau> 20)                   continue;

	}
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
  
    for(Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      if(m< nummolt) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	BoolMC =TMath::Abs((fSignTreeVariableInvMass - mass[type][m][v]))<sigmacentral[type][m][v]*sigma[type][m][v]; 
	if((isMC && !isEfficiency) || ishhCorr) {
	  BoolMC = kTRUE;
	}
	if(BoolMC && BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	  CounterV0NotSkipped[m] ++;
	  if (fSignTreeVariablePtV0> PtTrigMultDep[m]) CounterV0NotSkippedPtTr3[m] ++;

	  //some checks to see effect of DCA selections on purity of trigger particles:
	  CounterPairsBeforeSel[m][v] ++;
	  if(fSignTreeVariableDcaPosToPrimVertex < 0.1 ||  fSignTreeVariableDcaNegToPrimVertex <0.1) {
	    //	    continue;
	  }
	  CounterPairsAfterSel1[m][v] ++;
	  //***************************************************************************
	}
      }
      if (BoolVar){
	CounterAllPairs[m]++;
	if (fSignTreeVariablePtV0 > 3)      CounterPairsPtK0sGr3[m]++;
      }
    }
    if (SkipAssoc){    if (fSignTreeVariableSkipAssoc==1) continue;}


    Int_t counterpdg=0;
    for (Int_t i=0; i<10; i++){
      if (TMath::Abs(fSignTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])   {
        fHistTriggerPDGCode->Fill(i);
        counterpdg++;
      }
    }
if (counterpdg==0)   fHistTriggerPDGCode->Fill(10);

    //choose what particles use as trigger particles in MCtruth  (a previous selection might have been applied in the task!)
    Bool_t ispiKpemu=0;
    if (isMC && !isEfficiency){
      if (TriggerPDGChoice==1){
	for (Int_t i=0; i<5; i++){
	  if (TMath::Abs(fSignTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])  ispiKpemu = kTRUE;
	}
	if(!ispiKpemu)	  continue;
      }
      else if (TriggerPDGChoice==2){
	if(TMath::Abs(fSignTreeVariablePDGCodeTrigger) == 3222 || TMath::Abs(fSignTreeVariablePDGCodeTrigger)== 3112)	  continue;
      }
    }

    counterpdg=0;
    for (Int_t i=0; i<10; i++){
      if (TMath::Abs(fSignTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])   {
        fHistTriggerPDGCodeAfterSel->Fill(i);
        counterpdg++;
      }
    }
    if (counterpdg==0)   fHistTriggerPDGCodeAfterSel->Fill(10);
    
    //**********************************************************************************************

    fSignTreeVariableDeltaPhi = fSignTreeVariablePhiV0-fSignTreeVariablePhiTrigger; 
    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();


    for(Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      if(m< nummolt) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	if(BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	  CounterAllTrigger[m][v]++;
	}
      }
    }
    
    fHistDCATrigger->Fill(fSignTreeVariableDCAxy,  fSignTreeVariableDCAz);

    if (isMC){
      if (fSignTreeVariableisPrimaryTrigger !=1) {
	for(Int_t m=0; m<nummolt+1; m++){
	  if (MultBinning==1 && isHM){
	    if (m==0 || m==1) continue;
	  }
	  if(m< nummolt) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
	  else BoolVar=kTRUE;
	  for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	    if(BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      CounterNotPrimaryTrigger[m][v]++;
	    }
	  }
	}
	if (isPrimaryTrigger && !isEfficiency) continue;
      }
    }


    for(Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      //      if (m!= nummolt) continue;
      //      cout << " m " << m << endl;
      if(m< nummolt) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      hSign_DeltaEtaEtaV0[m]->Fill(fSignTreeVariableEtaV0, fSignTreeVariableDeltaEta);
      hSign_DeltaEtaEtaTrigger[m]->Fill(fSignTreeVariableEtaTrigger, fSignTreeVariableDeltaEta);
      for (Int_t DeltaEta=0; DeltaEta<numDeltaEta; DeltaEta++){
	if (BoolVar && fSignTreeVariableDeltaEta < DeltaEtaLimitUp[DeltaEta] && fSignTreeVariableDeltaEta > DeltaEtaLimitLow[DeltaEta]){
	  hSign_EtaV0[m][DeltaEta]->Fill(fSignTreeVariableEtaV0);
	  hSign_RapV0[m][DeltaEta]->Fill(fSignTreeVariableRapK0Short);
	  if (TMath::Abs(fSignTreeVariableRapK0Short) < 0.5)	    hSign_EtaV0MidRap[m][DeltaEta]->Fill(fSignTreeVariableEtaV0);
	}
      
      }
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  for(Int_t z=0; z<numzeta; z++){
	    /* defined above in a less error-prone way
	       if (type==4 || type==5 || type==8){
	       LimInfMass[type]=1.30;
	       LimSupMass[type]=1.342;

	       if (v >4)  {
	       LimInfMass[type]=1.29;
	       LimSupMass[type]=1.349;
	       }
	       }
	    */
	    if((isMC && !isEfficiency) || ishhCorr) {
	      BoolMC = kTRUE;
	      MassLimit=kTRUE;
	    }
	    else {
	      BoolMC =TMath::Abs((fSignTreeVariableInvMass - mass[type][m][v]))<sigmacentral[type][m][v]*sigma[type][m][v]; 
	      if(isSigma) MassLimit=(TMath::Abs((fSignTreeVariableInvMass - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && TMath::Abs((fSignTreeVariableInvMass - mass[type][m][v]))<nsigmamax*sigma[type][m][v]);
	      if(!isSigma)MassLimit=(TMath::Abs((fSignTreeVariableInvMass - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && fSignTreeVariableInvMass>LimInfMass[type][m][v] && fSignTreeVariableInvMass<LimSupMass[type][m][v]);
	    }	    
	    if(BoolMC && BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      hSign_PtTrigger[m][z]->Fill(fSignTreeVariablePtTrigger);
	      CounterACPairs[m][z]++;	
	      if (fSignTreeVariablePtV0>PtTrigMultDep[m]) 	      CounterACPairsPtTr3[m][z]++;	
	      if (fSignTreeVariablePtTrigger!= fSignTreeVariablePtTriggerTemp)	     {
		//		cout << "\nm " << m << " counted! " << endl;
		CounterTriggerCountedOnce[m][z]++;	
		hSign_PtTriggerCountedOnce[m][z]->Fill(fSignTreeVariablePtTrigger);
		hSign_PtTriggerCountedOnceBis[m][z]->Fill(fSignTreeVariablePtTrigger);
		fHistQA[2]->Fill(fSignTreeVariableMultiplicity);
	      }
	      hSign_PtV0[m][z]->Fill(fSignTreeVariablePtV0);
	    }	  
	    if(BoolVar && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      if(BoolMC){
		hSign_PtTriggerPtV0bins[m][z][v]->Fill(fSignTreeVariablePtTrigger);
		hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);

		effSign =0;	       
		if (isEtaEff){
		  for (Int_t pt=1; pt<= fHistEfficiencyV0PtEta[m]->GetNbinsX(); pt++){
		    if (fSignTreeVariablePtV0< fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinLowEdge(pt) || fSignTreeVariablePtV0>= fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinUpEdge(pt)) continue;
		    for (Int_t eta=1; eta<=fHistEfficiencyV0PtEta[m]->GetNbinsY(); eta++){
		      if (fSignTreeVariableEtaV0< fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinLowEdge(eta) || fSignTreeVariableEtaV0>= fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinUpEdge(eta)) continue;

		      effSign = fHistEfficiencyV0PtEta[m]->GetBinContent(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta)); 
		      sigmaEffSign = fHistEfficiencyV0PtEta[m]->GetBinError(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta)); 
		      /*
		      cout <<  "\neff sign "<<   fHistEfficiencyV0PtEta[m]->GetBinContent(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta)) << endl;
		      cout << " pt: " << fSignTreeVariablePtV0 << " eta: " << fSignTreeVariableEtaV0 << endl;
		      cout <<  fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinLowEdge(pt) << " - " << fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinUpEdge(pt) << endl;
		      cout <<  fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinLowEdge(eta) << " - " << fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinUpEdge(eta) << endl;
		      cout << "effw Sign: " << effSign << endl;
		      */

		    }
		  }
		  if (effSign!=0){
		    hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi, 1./effSign);
		    hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi, sigmaEffSign/effSign);
		  }
		}
	      }
	      if((!isMC || (isMC &&isEfficiency)) && MassLimit){
		hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }

	    }
	  }
	}
      }
    }
    //    cout << fSignTreeVariablePtTrigger << " temp " << fSignTreeVariablePtTriggerTemp << endl;
    fSignTreeVariablePtTriggerTemp = fSignTreeVariablePtTrigger;
  }


  BoolVar=kFALSE;
  BoolMC=kFALSE;
  MassLimit=kFALSE;

  Float_t     fBkgTreeVariableInvMass= 0;
  cout << "ciao " << endl;
  cout << "\n\n I will process " << EntriesBkg << " entries for theSE correlation " << endl;
  Double_t effBkg=0;
  Double_t sigmaEffBkg=0;
  for(Int_t k = 0; k<EntriesBkg; k++){
    if (k>1000000) continue;
    tBkg->GetEntry(k);     
    //  for(Int_t k = 0; k<1000000; k++){
    for (Int_t l=0; l<10000; l++){
      //if (k ==100000*l) cout << "k = " << k << " over a total of " << EntriesBkg << endl;
      if (k ==100000*l) {
	Float_t kFloat= k;
	//	cout << " etav0:" << fBkgTreeVariableEtaV0 <<  " etatrigger: " << fBkgTreeVariableEtaTrigger << " deltaeta: " << fBkgTreeVariableDeltaEta << "  new: " << fBkgTreeVariableEtaV0-fBkgTreeVariableEtaTrigger <<  endl;
	//	cout << " phiv0:" << fBkgTreeVariablePhiV0 <<  " phitrigger: " << fBkgTreeVariablePhiTrigger << " deltaphi: " << fBkgTreeVariableDeltaPhi << "  new: " << fBkgTreeVariablePhiV0-fBkgTreeVariablePhiTrigger <<  endl;
	cout << "PtTrigMin "<< PtTrigMin << " sysV0 "<< sysV0 << " ME processing..." << kFloat/EntriesBkg << endl;
      }
    }


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
    if (israp==0){    
      if (isEta05) {
	if (TMath::Abs(fBkgTreeVariableEtaV0)>0.5)continue;
      }
      else{
	if (TMath::Abs(fBkgTreeVariableEtaV0)>0.8)continue;
      }
    }
    else if (israp==1){
      if (TMath::Abs(fBkgTreeVariableRapK0Short)>0.5)continue;
    }

    //definition of true particle
    if (type<=2) isParticleTrue= (fBkgTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==3) isParticleTrue= ((fBkgTreeVariablePDGCodeAssoc==PDGCode[1]) ||(fBkgTreeVariablePDGCodeAssoc==PDGCode[1])  );
 
    if (SkipAssoc){    if (fBkgTreeVariableSkipAssoc==1) continue;}

    if (IsTrueParticle && isMC && isEfficiency==1){
      if (!isParticleTrue) continue;
      if (fBkgTreeVariableisPrimaryV0!=1) continue;
    } 

    //************cuts on pT trigger min*********************************
    if (!IsPtTrigMultDep){
      if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMin) continue;
    }
    else {
      Bool_t IsSelectionPtMinNotPassedBkg=0;
      for(Int_t m=0; m<nummolt; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
	if ( fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1]){
	  if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMultDep[m]) IsSelectionPtMinNotPassedBkg=1;
	}
      }
      if (IsSelectionPtMinNotPassedBkg) continue;
    }

    if(TMath::Abs(fBkgTreeVariablePtTrigger)>ptjmax) continue;

    if (isMC && !isEfficiency && isPrimaryTrigger){
      if (fBkgTreeVariableisPrimaryTrigger !=1) continue;
    }

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
	    if(TMath::Abs((fBkgTreeVariableInvMassLambda - massLambda))<= 0.005) continue;
	    if(TMath::Abs((fBkgTreeVariableInvMassAntiLambda - massLambda))<= 0.005) continue;
	  }
	  if(fBkgTreeVariableV0CosineOfPointingAngle<0.995)            continue;
	  if(fBkgTreeVariableDcaNegToPrimVertex < 0.06)      continue;
	  if(fBkgTreeVariableDcaPosToPrimVertex < 0.06)      continue;
	  if(fBkgTreeVariableDcaV0ToPrimVertex > 0.5)                               continue;
	  if(fBkgTreeVariablectau> 20)                   continue;

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

    //do not take sigma as trigger particles in MCtruth 
    Bool_t ispiKpemu=0;
    if (isMC && !isEfficiency){
      if (TriggerPDGChoice==1){
	for (Int_t i=0; i<5; i++){
	  if (TMath::Abs(fBkgTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])  ispiKpemu = kTRUE;
	}
	if(!ispiKpemu)	  continue;
      }
      else if (TriggerPDGChoice==2){
	if(TMath::Abs(fBkgTreeVariablePDGCodeTrigger) == 3222 || TMath::Abs(fBkgTreeVariablePDGCodeTrigger)== 3112)	  continue;
       
      }
    }

    //**********************************************************************************************

    Double_t fBkgTreeVariableDeltaPhiMinus = -fBkgTreeVariableDeltaPhi ;
    fBkgTreeVariableDeltaPhi = fBkgTreeVariablePhiV0-fBkgTreeVariablePhiTrigger; 

    if (fBkgTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fBkgTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fBkgTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fBkgTreeVariableDeltaPhi += 2.0*TMath::Pi();

    /*
      cout << " deltaphi: " << fBkgTreeVariableDeltaPhi << "  new: " << fBkgTreeVariablePhiV0-fBkgTreeVariablePhiTrigger <<  endl;

      if (fBkgTreeVariableDeltaPhiMinus >  (1.5*TMath::Pi())) fBkgTreeVariableDeltaPhiMinus -= 2.0*TMath::Pi();
      if (fBkgTreeVariableDeltaPhiMinus < (-0.5*TMath::Pi())) fBkgTreeVariableDeltaPhiMinus += 2.0*TMath::Pi();

      cout << " New: deltaphi: " << fBkgTreeVariableDeltaPhiMinus << "  new: " << fBkgTreeVariablePhiV0-fBkgTreeVariablePhiTrigger <<  endl;
    */
    for(Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      //      if (m!= nummolt) continue;
      if(m< nummolt) BoolVar =  fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      hBkg_DeltaEtaEtaV0[m]->Fill(fBkgTreeVariableEtaV0, fBkgTreeVariableDeltaEta);
      hBkg_DeltaEtaEtaTrigger[m]->Fill(fBkgTreeVariableEtaTrigger, fBkgTreeVariableDeltaEta);
      for (Int_t DeltaEta=0; DeltaEta<numDeltaEta; DeltaEta++){
	if (BoolVar && fBkgTreeVariableDeltaEta < DeltaEtaLimitUp[DeltaEta] && fBkgTreeVariableDeltaEta > DeltaEtaLimitLow[DeltaEta]){
	  hBkg_EtaV0[m][DeltaEta]->Fill(fBkgTreeVariableEtaV0);
	  hBkg_RapV0[m][DeltaEta]->Fill(fBkgTreeVariableRapK0Short);
	  if (TMath::Abs(fBkgTreeVariableRapK0Short) < 0.5)	    hBkg_EtaV0MidRap[m][DeltaEta]->Fill(fBkgTreeVariableEtaV0);
	}
      }

      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=PtBinMin; v<numPtV0Max; v++){

	    /*	    if (type==4 || type==5 || type==8){
		    LimInfMass[type]=1.30;
		    LimSupMass[type]=1.342;

		    if (v >4)  {
		    LimInfMass[type]=1.29;
		    LimSupMass[type]=1.352;
		    }
		    }*/
	    if((isMC && !isEfficiency) || ishhCorr) {
	      BoolMC = kTRUE;
	      MassLimit=kTRUE;
	    }
	    else {
	      BoolMC =TMath::Abs((fBkgTreeVariableInvMass - mass[type][m][v]))<sigmacentral[type][m][v]*sigma[type][m][v]; 
	      if(isSigma) MassLimit=(TMath::Abs((fBkgTreeVariableInvMass - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && TMath::Abs((fBkgTreeVariableInvMass - mass[type][m][v]))<nsigmamax*sigma[type][m][v]);
	      if(!isSigma)MassLimit=(TMath::Abs((fBkgTreeVariableInvMass - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && fBkgTreeVariableInvMass>LimInfMass[type][m][v] && fBkgTreeVariableInvMass<LimSupMass[type][m][v]);
	    }	    
 
	    if(BoolMC && BoolVar &&  fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      hBkg_PtTrigger[m][z]->Fill(fBkgTreeVariablePtTrigger);
	      hBkg_PtV0[m][z]->Fill(fBkgTreeVariablePtV0);
	    }	  
	    if(BoolVar && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      if(BoolMC){

		hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);

		effBkg =0;	       
		if (isEtaEff){
		  for (Int_t pt=1; pt<= fHistEfficiencyV0PtEta[m]->GetNbinsX(); pt++){
		    if (fBkgTreeVariablePtV0< fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinLowEdge(pt) || fBkgTreeVariablePtV0>= fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinUpEdge(pt)) continue;
		    for (Int_t eta=1; eta<=fHistEfficiencyV0PtEta[m]->GetNbinsY(); eta++){
		      if (fBkgTreeVariableEtaV0< fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinLowEdge(eta) || fBkgTreeVariableEtaV0>= fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinUpEdge(eta)) continue;

		      effBkg = fHistEfficiencyV0PtEta[m]->GetBinContent(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta)); 
		      sigmaEffBkg = fHistEfficiencyV0PtEta[m]->GetBinError(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta)); 
		      /*
		      cout <<    fHistEfficiencyV0PtEta[m]->GetBinContent(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta)) << endl;
		      cout << " pt: " << fBkgTreeVariablePtV0 << " eta: " << fBkgTreeVariableEtaV0 << endl;
		      cout <<  fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinLowEdge(pt) << " - " << fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinUpEdge(pt) << endl;
		      cout <<  fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinLowEdge(eta) << " - " << fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinUpEdge(eta) << endl;
		      cout << "effw Bkg: " << effBkg << endl;
		      */

		    }
		  }

		  if (effBkg!=0){
		    hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, 1./effBkg);
		    hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, sigmaEffBkg/effBkg);
		  }
		}
	      }
	      if((!isMC || (isMC &&isEfficiency)) && MassLimit){
		hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	      }
	   
	    }
	  }
	}
      }
    }
  }


  //put correct errors on 2D histograms
  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==1 && isHM){
      if (m==0 || m==1) continue;
    }
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	  Int_t bin=0;
	  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]);
          hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]);
	  for (Int_t dphi=1; dphi<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetNbinsY(); dphi++){ 
	    for (Int_t deta=1; deta<= hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetNbinsX(); deta++){
	      bin = hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBin(deta, dphi);
	      //	      if (hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin) !=0) {
	      /*
		cout << " bin " << bin << endl;
		cout << " sqrt of content : " << hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin) << endl;
		cout << " error set by ROOT: " << hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinError(bin) << endl;
	      */
	      if (isEtaEff){
		EffRelErrSign[m][v] = hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->GetBinContent(bin);
		//		EffRelErr[m][v] = fHistEfficiencyV0PtPtBins[m]->GetBinError(v+1)/fHistEfficiencyV0PtPtBins[m]->GetBinContent(v+1); 
		//not correct! 		hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->SetBinError(bin,sqrt(hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->GetBinContent(bin)));
		if (hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin)!=0) {
		  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->SetBinError(bin,sqrt(1./hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(bin) + pow(EffRelErrSign[m][v],2)) * hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		  //		  cout << "err assoc to efficiency: " << EffRelErrSign[m][v]<< " err2 associated to counts " << 1./hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(bin) << endl;
		}
		else{
		  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->SetBinError(bin,0);
		}
		//not correct!		hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->SetBinError(bin,sqrt(hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->GetBinContent(bin)));
		EffRelErrBkg[m][v] = hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->GetBinContent(bin);
		if (hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin)!=0) {
		  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->SetBinError(bin,sqrt(1./hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(bin) + pow(EffRelErrBkg[m][v],2)) * hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		  //		  cout << "err2 assoc to efficiency: " << pow(EffRelErrBkg[m][v],2)<< " error associated to counts " << 1./hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(bin) << endl;
		}
		else{
		  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->SetBinError(bin,0);
		}

		if (hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin)!=0) {
		  hDeltaEtaDeltaPhi_SEbinsEffwRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		if (hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin)!=0) {
		  hDeltaEtaDeltaPhi_MEbinsEffwRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		}
	      }
	      if (hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(bin) !=0){
		hDeltaEtaDeltaPhi_SEbinsRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(bin));
	      }
	      if (hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(bin) !=0){
		hDeltaEtaDeltaPhi_MEbinsRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(bin));
	      }
	    }
	  }
	}
      }
    }
  }


  fHistTriggerPDGCode->Scale(1./fHistTriggerPDGCode->GetEntries());
  fHistTriggerPDGCodeAfterSel->Scale(1./fHistTriggerPDGCodeAfterSel->GetEntries());
  fout->WriteTObject(fHistTriggerPDGCode);
  fout->WriteTObject(fHistTriggerPDGCodeAfterSel);

  TLegend * legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  legend->SetHeader("Multiplicity classes");
  TLegend * legendLow = new TLegend(0.6, 0.1, 0.9, 0.4);

  for(Int_t z=0; z<numzeta; z++){
    for (Int_t m =0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      hSign_PtTriggerCountedOnceBis[m][z]->Sumw2();
      hSign_PtTrigger[m][z]->Sumw2();

      hSign_PtTriggerCountedOnceBis[5][z]->Sumw2();
      hSign_PtTrigger[5][z]->Sumw2();


      if (m==0){
	//      hSign_PtTriggerCountedOnceBis[5][z]->Rebin(4);
	hSign_PtTrigger[5][z]->Rebin(4);
	hSign_PtTrigger[5][z]->Scale(1./ CounterACPairs[5][z]/4);
	hSign_PtTriggerCountedOnceBis[5][z]->Scale(1./ CounterTriggerCountedOnce[5][z]);
	for (Int_t i = 1; i<= hSign_PtTriggerCountedOnceBis[5][z]->GetNbinsX(); i++){
	  hSign_PtTriggerCountedOnceBis[5][z]->SetBinContent(i, hSign_PtTriggerCountedOnceBis[5][z]->GetBinContent(i)/hSign_PtTriggerCountedOnceBis[5][z]->GetBinWidth(i));
	}

      }

      if (m!=nummolt){
	//      hSign_PtTriggerCountedOnceBis[m][z]->Rebin(4);
	hSign_PtTrigger[m][z]->Rebin(4);
	hSign_PtTrigger[m][z]->Scale(1./ CounterACPairs[m][z]/4);
	hSign_PtTriggerCountedOnceBis[m][z]->Scale(1./ CounterTriggerCountedOnce[m][z]);
	for (Int_t i = 1; i<= hSign_PtTriggerCountedOnceBis[5][z]->GetNbinsX(); i++){
	  hSign_PtTriggerCountedOnceBis[m][z]->SetBinContent(i, hSign_PtTriggerCountedOnceBis[m][z]->GetBinContent(i)/hSign_PtTriggerCountedOnceBis[m][z]->GetBinWidth(i));
	}

      }


      hSign_PtTriggerCountedOnceRatio[m][z]=(TH1D*)      hSign_PtTriggerCountedOnceBis[m][z]->Clone("hSign_PtTriggerCountedOnceRatio"+Smolt[m]);
      hSign_PtTriggerCountedOnceRatio[m][z]->Divide(      hSign_PtTriggerCountedOnceBis[5][z]);
      hSign_PtTriggerRatio[m][z]=(TH1D*)      hSign_PtTrigger[m][z]->Clone("hSign_PtTriggerRatio"+Smolt[m]);


      cout <<       hSign_PtTrigger[m][z]->GetBinContent(50) << endl;
      cout <<       hSign_PtTriggerRatio[m][z]->GetBinContent(50) << endl;
      hSign_PtTriggerRatio[m][z]->Divide(      hSign_PtTrigger[5][z]);
      cout <<       hSign_PtTriggerRatio[m][z]->GetBinContent(50) << endl;
      cout << m << endl;
      fHistQA[0]->SetBinContent(m+1, hSign_PtTrigger[m][z]->GetMean());
      fHistQA[0]->SetBinError(m+1, hSign_PtTrigger[m][z]->GetMeanError());
      fHistQA[1]->SetBinContent(m+1, hSign_PtTriggerCountedOnce[m][z]->GetMean());
      fHistQA[1]->SetBinError(m+1, hSign_PtTriggerCountedOnce[m][z]->GetMeanError());
      fHistQA[3]->SetBinContent(m+1,  CounterACPairs[m][z]/CounterTriggerCountedOnce[m][z]);
      //fHistQA[4]->SetBinContent(m+1, CounterACPairs[m][z]/ CounterV0NotSkipped[m]);
      fHistQA[4]->SetBinContent(m+1, CounterACPairsPtTr3[m][z]/ CounterV0NotSkippedPtTr3[m]);

      canvasQA[0]->cd(1);
      fHistQA[0]->GetYaxis()->SetRangeUser(PtTrigMin+1, PtTrigMin+1.5);
      fHistQA[0]->Draw("");

      canvasQA[0]->cd(2);
      gPad->SetLogy();
      hSign_PtTrigger[m][z]->GetXaxis()->SetRangeUser(PtTrigMin, 6);
      hSign_PtTrigger[m][z]->GetYaxis()->SetRangeUser(0.01, 0.15);
      hSign_PtTrigger[m][z]->Draw("same p");
      if (m==nummolt) legend->Draw("");

      canvasQA[0]->cd(4);
      hSign_PtTriggerRatio[m][z]->GetYaxis()->SetRangeUser(0.98, 1.12);
      if (m!=nummolt)      hSign_PtTriggerRatio[m][z]->Draw("same l");
      lineat1->Draw("same");

      canvasQA[1]->cd(2);
      gPad->SetLogy();
      legend->AddEntry( hSign_PtTriggerCountedOnce[m][z], SBismolt[m], "pl");
      legendLow->AddEntry( hSign_PtTriggerCountedOnce[m][z], SBismolt[m], "pl");
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetRangeUser(PtTrigMin, 15);
      hSign_PtTriggerCountedOnceBis[m][z]->GetYaxis()->SetRangeUser(0.001, 1);
      hSign_PtTriggerCountedOnceBis[m][z]->Draw("same p");
      if (m==nummolt) legend->Draw("");

      canvasQA[1]->cd(4);
      hSign_PtTriggerCountedOnceRatio[m][z]->GetYaxis()->SetRangeUser(0.9, 1.1);
      if (m!=nummolt)      hSign_PtTriggerCountedOnceRatio[m][z]->Draw("same l");
      lineat1->Draw("same");

      canvasQA[1]->cd(1);
      fHistQA[1]->GetYaxis()->SetRangeUser(PtTrigMin+1, PtTrigMin+1.5);
      fHistQA[1]->Draw("");

      for(Int_t v=PtBinMin; v<numPtV0Max; v++){
	fHistPtTriggervsPtAssoc[m]->SetBinContent(v+1,hSign_PtTriggerPtV0bins[m][z][v]->GetMean());
	fHistPtTriggervsPtAssoc[m]->SetBinError(v+1,hSign_PtTriggerPtV0bins[m][z][v]->GetMeanError());
	fHistNonPrimaryTrigger[m]->SetBinContent(v+1, CounterNotPrimaryTrigger[m][v]/CounterAllTrigger[m][v]);
	fHistNonPrimaryTrigger[m]->SetBinError(v+1, SetEfficiencyError(CounterNotPrimaryTrigger[m][v],CounterAllTrigger[m][v]));
	fHistEfficiencyReduction[m]->SetBinContent(v+1, CounterPairsAfterSel1[m][v]/CounterPairsBeforeSel[m][v]);
	fHistEfficiencyReduction[m]->SetBinError(v+1, SetEfficiencyError(CounterPairsAfterSel1[m][v],CounterPairsBeforeSel[m][v]));

      }
      fHistPtTriggervsPtAssoc[m]->GetYaxis()->SetRangeUser(3.5, 7);
      fHistPtTriggervsPtAssoc[m]->SetMarkerColor(Color[m]);
      fHistPtTriggervsPtAssoc[m]->SetLineColor(Color[m]);
      fHistPtTriggervsPtAssoc[m]->SetMarkerStyle(33);
      fHistNonPrimaryTrigger[m]->GetYaxis()->SetRangeUser(0, 0.15);
      fHistNonPrimaryTrigger[m]->SetMarkerColor(Color[m]);
      fHistNonPrimaryTrigger[m]->SetLineColor(Color[m]);
      fHistNonPrimaryTrigger[m]->SetMarkerStyle(33);
      fHistEfficiencyReduction[m]->GetYaxis()->SetRangeUser(0.7, 1);
      fHistEfficiencyReduction[m]->SetMarkerColor(Color[m]);
      fHistEfficiencyReduction[m]->SetLineColor(Color[m]);
      fHistEfficiencyReduction[m]->SetMarkerStyle(33);

      canvasQA[4]->cd();
      fHistPtTriggervsPtAssoc[m]->Draw("same p");
      if (m==nummolt) legendLow->Draw("");

      canvasQA[6]->cd();
      fHistNonPrimaryTrigger[m]->Draw("same p");
      if (m==nummolt) legendLow->Draw("");

      canvasQA[7]->cd();
      fHistEfficiencyReduction[m]->Draw("same p");
      if (m==nummolt) legendLow->Draw("");

      canvasQA[6] -> SaveAs(pathoutpdf + "_NonPrimaryTriggerFraction.pdf");
      fout->WriteTObject(fHistPtTriggervsPtAssoc[m]);
      fout->WriteTObject(fHistNonPrimaryTrigger[m]);
      fout->WriteTObject(fHistEfficiencyReduction[m]);
    }
    
    fout->WriteTObject(fHistDCATrigger);
    canvasQA[5] ->cd(2);
    gPad->SetLogy();
    for (Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      fHistPtMaxvsMultBefAllProj[m]->Sumw2();
      fHistPtMaxvsMultBefAllProj[m]->SetLineColor(Color[m]);
      fHistPtMaxvsMultBefAllProj[m]->SetMarkerColor(Color[m]);
      fHistPtMaxvsMultBefAllProj[m]->SetMarkerStyle(33);
      fHistPtMaxvsMultBefAllProj[m]->GetXaxis()->SetRangeUser(PtTrigMin, 15);
      fHistPtMaxvsMultBefAllProj[m]->GetYaxis()->SetRangeUser(0.001, 1);
      fHistPtMaxvsMultBefAllProj[m]->Scale(fHistPtMaxvsMultBefAllProj[m]->GetBinWidth(1));
      cout << fHistPtMaxvsMultBefAllProj[m]->Integral()<< endl;
      fHistPtMaxvsMultBefAllProj[m]->Scale(1./fHistPtMaxvsMultBefAllProj[m]->Integral());
      fHistPtMaxvsMultBefAllProj[m]->Draw("same e");
      legend->Draw("");
    }
    canvasQA[5] ->cd(4);
    for (Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      if(m!=nummolt)      fHistPtMaxvsMultBefAllProj[m]->Rebin(5);
      if (m==0)      fHistPtMaxvsMultBefAllProj[5]->Rebin(5);
      fHistPtMaxvsMultBefAllProjRatio[m]= (TH1F*)     fHistPtMaxvsMultBefAllProj[m]->Clone(Form("fHistPtMaxvsMultBefAllProjRatio_m%i",m));
      fHistPtMaxvsMultBefAllProjRatio[m]->Divide(    fHistPtMaxvsMultBefAllProj[5]);
      fHistPtMaxvsMultBefAllProjRatio[m]->GetYaxis()->SetRangeUser(0.8, 1.2);
      fHistPtMaxvsMultBefAllProjRatio[m]->GetYaxis()->SetTitle("");
      if (m!=nummolt)      fHistPtMaxvsMultBefAllProjRatio[m]->Draw("same ep");
      //      legend->Draw("");
    }


    canvasQA[2]->cd();
    fHistQA[2]->Scale(1./fHistQA[2]->GetEntries());
    fHistQA[2]->Draw();
    canvasQA[3]->cd();
    fHistQA[3]->Draw();

    for (Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==1 && isHM){
	if (m==0 || m==1) continue;
      }
      for(Int_t z=0; z<numzeta; z++){
	fout->WriteTObject(hSign_PtTriggerCountedOnceBis[m][z]);
	fout->WriteTObject(hSign_PtTriggerCountedOnce[m][z]);
      }
    }
    for (Int_t i = 0; i<numQAhisto; i++){
      fout->WriteTObject(fHistQA[i]);
      fout->WriteTObject(canvasQA[i]);
    }
  }

  fout->Write();

  for (Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==1 && isHM){
      if (m==0 || m==1) continue;
    }
    cout << "m " << m << " fracion of pairs with pt,K0s > 3 GeV/c" <<  CounterPairsPtK0sGr3[m]/CounterAllPairs[m]<< endl;
  }

  cout << "partendo dal file " << PathIn << "\n e " << PathInMassDef << " (per le diverse molteplicità) \ne dal file delle efficienze "<<   PathInEfficiency << ",  \nho creato il file " << PathOut<< endl;
}

