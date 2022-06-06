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
Double_t SetEfficiencyError(Int_t k, Int_t n){
  return sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
}

void readTreePLChiaraCasc_second(Int_t type=4, Bool_t SkipAssoc=0, Int_t israp=0, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Float_t PtTrigMin=3,Float_t PtTrigMinFit=3, Int_t sysTrigger=0, Int_t sysV0=0,Int_t syst=0,bool isMC = 1, Bool_t isEfficiency=0,TString year0="2016", TString year="PythiaRopes"/*"17pq_hXi"/*"161718_HM_hXi_WithFlat16k_No18p"/*"161718_HM_hXi"/*"161718Full_AOD234_hXi"/*"17pq_pp5TeV_hXi_pttrig0.15"/*"161718_hXi"/*"LHC17_AOD234_Red"/*"2018f1_extra_10runs_hXi_Rlabel_Hybrid"/*"2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2"/*"161718_MD_EtaEff_hXi"/*"2018f1g4_extra_EtaEff_hXi"/"161718_MD_New_hXi"/*"LHC16_17_GP_Hybrid_hXi"/"2018g4_extra_hXi_SelTrigger"/*"AllMC_hXi"/"Run2DataRed_MECorr_hXi"/*"2016k_pass2_TOFOOBPileUp"*/,  TString Path1 =""/*"_PtTrigMax2.5"/*"NewMultClassBis"*/,  Double_t ptjmax =15, Double_t nsigmamax=10, Bool_t isSigma=kFALSE, Int_t MultBinning=3, Bool_t IsTrueParticle=0,Int_t TriggerPDGChoice=0, Bool_t isEtaEff=1, TString yearMC /*for efficiency*/="17pq_hXi"/*"161718_HM_hXi"/*"161718_HM_hXi_LowPtTrig"/*"161718_MD_EtaEff_LowPtTrig_hXi"/"161718Full_AOD235_hXi"/*"17pq_pp5TeV_hXi_pttrig0.15"/*"17pq_hXi"/*"161718_MD_EtaEff_LowPtTrig_hXi"/"161718_MD_EtaEff_hXi"/*"2018f1g4_extra_EtaEff_hXi"*/, Bool_t isEstimateRun3=0, Bool_t isNewInputPath=1, Bool_t isHM=0, Bool_t isMCForNorm=0, Int_t SidebandsSide =0, Bool_t isGenOnTheFly=1){

  //isGenOnTheFly --> events were generated on the fly and only the kinematic part is saved; the multiplicity distribution in percentile classes is not abvailable, instead classes based on the number of particles in the V0 acceptance are used
  if (!isMC && isGenOnTheFly) return;
  if (isGenOnTheFly) { //these variabes have no meaning for the MCtruth analysis -- they are set to zero in order not to appear in output file name
    isBkgParab = 0;
    MultBinning = 0;
    isHM =0;
  }

  //SidebandsSide =0-> left + right SB
  //SidebandsSide =1-> left SB
  //SidebandsSide =2-> right SB

  //isMCForNorm =1 ! for computation of normalisation factor (same binning used for data)
  if (year=="17pq_hXi") MultBinning=3;
  if (year=="17pq_pp5TeV_hXi_pttrig0.15") MultBinning=3;
  //isEstimateRun3=1: I use smaller mult classes
  //isMeanFixedPDG and isBkgParab are characteristics of the fit to the inv mass distributions 
  cout << isMC << endl;
  cout << " Pt Trigg Min è = " << PtTrigMin << endl;

  if (!isMC){
    IsTrueParticle=0;
  }
  if (isMC && !isEfficiency) isEtaEff =0;
  //***************************************************************************                               
  //******** What particles to choose as trigger whan analysing MC truth? *****   

  // Since in data we are able to identify basically only pi, K, p, mu, e, also when analyzing the MC truth we want to take only these particles as trigger particles. And we want to exclude all the others already WHEN CHOOSING the trigger particle, so at the level of the TASK!                                                     
  //Here some offline selections are implemented as well, but they are useless if the selection is done at the level of the task                                                                                              
  //Check if fHistTriggerPDGCode and fHistTriggerPDGCodeAfterSel are the same. If it is so, the selections are correctly implemtned in the task 

  //***************************************************************************                               

  TString sTriggerPDGChoice[3] = {"", "_IsOnlypiKpemu", "_IsNotSigmaOnly"};
  if (!isMC) TriggerPDGChoice=0;

  if (israp>1) return;
  if (type>5) {cout << "type value not allowed" << endl; return;}
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

  if (isEstimateRun3 && PtTrigMinFit!=3) { cout << " PtTrigMinFit should be 3 because no fit was performed at larger values, and they are not even necessary for the estimate I'm going to do " << endl; return; }

  TF1 * lineat1 = new TF1("pol0", "pol0", 0, 30);
  lineat1->FixParameter(0,1);
  lineat1->SetLineColor(kBlack);

  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  const Float_t ctauK0s=2.6844;

  TString file;

  const Int_t numtipo=6;
  Int_t nummolt=12;
  Int_t nummoltMax =12;
  if (!isGenOnTheFly) nummoltMax = 5; //11 for estimate run3 (isEstimateRun3)
  const Int_t numzeta=1;
  const Int_t numPtV0=8; //14; //was 8 in my analysis
  const Int_t numPtTrigger=1;
  Int_t numPtV0Min =1;
  if (isMC && !isEfficiency && !isMCForNorm) numPtV0Min=0;

  if (!isEstimateRun3 && nummoltMax!=5 && !isGenOnTheFly) {cout << " maybe you forgot the multiplciity classes you used for Run3 estimates" << endl; return;}
  if (isEstimateRun3) PtTrigMinFit=3; 
  if (isEstimateRun3 && nummoltMax!=11) {cout << " maybe you forgot the multiplciity classes you use for your analysis" << endl; return;}

  TString SmoltBin0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t NmoltBin0[nummolt+1]={0, 5, 10, 30, 50, 100};
  TString SmoltBin1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  Double_t NmoltBin1[nummolt+1]={0, 1, 5, 15, 30, 100};
  TString SmoltBin2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  Double_t NmoltBin2[nummolt+1]={0, 2, 7, 15, 30, 100};
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};
  TString SmoltGenOnTheFly[nummolt+1]={"0-15", "15-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-160", "160-300", "0-300"};
  TString SBismoltGenOnTheFly[nummolt+1]={"0-15", "15-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-160", "160-300", "0-300"};
  Double_t NmoltGenOnTheFly[nummolt+1]={0, 15, 30, 45, 54, 63, 72, 81, 90, 105, 120, 160, 300};

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  
  for (Int_t m=0; m<nummoltMax+1; m++){
    if (MultBinning==0){
      Smolt[m] = SmoltBin0[m];
      Nmolt[m] = NmoltBin0[m];
    }
    else if (MultBinning==1){
      Smolt[m] = SmoltBin1[m];
      Nmolt[m] = NmoltBin1[m];
    }
    else if (MultBinning==2){
      Smolt[m] = SmoltBin2[m];
      Nmolt[m] = NmoltBin2[m]; 
    }
    else if (MultBinning==3){
      Smolt[m] = Smoltpp5TeV[m];
      Nmolt[m] = Nmoltpp5TeV[m]; 
    }
  }
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

  TString SmoltRun3[12] = {"0-1", "1-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-40", "40-50", "50-70", "70-100", "all"} ;
  Double_t NmoltRun3[12] = {0,1,5,10,15,20,25,30,40,50,70,100};

  if (isEstimateRun3){
    for(Int_t molt=0; molt<nummoltMax+1; molt++){
      Smolt[molt] = SmoltRun3[molt];
      Nmolt[molt] = NmoltRun3[molt];
    }
  }
  if (isGenOnTheFly){
    for (Int_t m=0; m<nummoltMax+1; m++){
      Smolt[m] = SmoltGenOnTheFly[m];
      Nmolt[m] = NmoltGenOnTheFly[m];
    }
  }

  TString tipo[numtipo]={"XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
  Float_t ctauCasc[numtipo] = {4.91,4.91,  2.461, 2.461, 4.91, 2.461}; //cm , average ctau of Xi and Omega                                
  Float_t PDGCode[numtipo-2] = {3312, -3312, 3334, -3334};

  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};

  TString PathIn="FinalOutput/AnalysisResults";
  TString PathOut="FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";
  TString PathInMass= "FinalOutput/DATA" + year0 + "/invmass_distribution_thesis/invmass_distribution";
  TString PathInMassDef;
  PathIn+=year;
  PathOut+=year;  
  // PathInMass+=year;
  
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
  //  PathInMass+=Path1;
  PathIn+=".root";

  PathOut+=Path1;
  //  PathOut+="_NewSB";
  PathOut+="_"; 
  PathOut +=tipo[type];
  PathOut +=Srap[israp];
  PathOut +=SSkipAssoc[SkipAssoc];
  PathOut +=Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f",sysTrigger, sysV0, syst, PtTrigMin); 
  if (IsTrueParticle) PathOut+="_IsParticleTrue";
  //  PathOut+="_ciao";
  PathOut+= sTriggerPDGChoice[TriggerPDGChoice];
  if (isEstimateRun3) PathOut += "_IsEstimateRun3";
  //  PathOut+="_thinptbins"; 
  //  PathOut+= "_Try1";
  if (SidebandsSide ==1) PathOut+= "_LeftSB";
  else if (SidebandsSide ==2) PathOut+= "_RightSB";
  if (isEtaEff) PathOut+="_isEtaEff";
  //  PathOut+="_isEtaEff";
  if (MultBinning!=0) PathOut+=Form("_MultBinning%i",MultBinning);
  if (isMC && !isEfficiency && !isMCForNorm) PathOut += "_MCPrediction";
  TString PathOutPdf = PathOut;
  PathOut+= ".root";

  cout << "file di input " << PathIn << endl;
  TFile *fin = new TFile(PathIn);

  //input file where efficiency is found                                                               
  TString PathInEfficiency = "FinalOutput/DATA2016/Efficiency/Efficiency" +yearMC;
  PathInEfficiency+=Path1;// change                                                                   
 
  if(type>=0){
    PathInEfficiency +="_";
    PathInEfficiency +=tipo[type];
    PathInEfficiency +=Srap[israp];
    PathInEfficiency +=SSkipAssoc[SkipAssoc];
  }

  PathInEfficiency+= Form("_SysT%i_SysV0%i_PtMin%.1f", sysTrigger, sysV0, PtTrigMin);
  /*    												       if (IsPtTrigMultDep)    PathInEfficiency+= "_IsPtTrigMultDep"; 
 if (IsEfficiencyMassSel)    PathInEfficiency+= "_isEfficiencyMassSel";                              
  */
  //  PathInEfficiency+= "_thinptbins";
  if (MultBinning!=0) PathInEfficiency+=Form("_MultBinning%i",MultBinning);
  PathInEfficiency+= ".root";

  TFile *fileinEfficiency = new TFile (PathInEfficiency, "");
  TH2F * fHistEfficiencyV0PtEta[nummolt+1];
  TH1F * fHistEfficiencyV0PtPtBins[nummolt+1];
  if (isEtaEff){
    if (!fileinEfficiency) {cout << " input file with efficiency not found " << endl; return;}
  }
  if (isEtaEff){
    for(Int_t molt=0; molt<nummoltMax+1; molt++){
      if (MultBinning==3 && (molt==2 || molt==3 || molt==4)) continue; 
      if (MultBinning==1 && isHM && molt <= 1) continue;
      cout << molt << endl;
      cout << Smolt[molt] << endl;
      fHistEfficiencyV0PtEta[molt] = (TH2F*) fileinEfficiency->Get("fHistV0EfficiencyPtV0EtaV0PtBins_"+ Smolt[molt]);
      if (!fHistEfficiencyV0PtEta[molt]) {cout << "histogram 2D V0 efficiency pt vs eta " << endl; return;}
      fHistEfficiencyV0PtPtBins[molt] = (TH1F*) fileinEfficiency->Get("fHistV0EfficiencyPtBins_"+ Smolt[molt]);
      if (!fHistEfficiencyV0PtPtBins[molt]) {cout << "histogram 1D V0 efficiency pt " << endl; return;}
    }
  }

  Float_t   EffRelErrSign[nummolt+1][numPtV0]={0};
  Float_t   EffRelErrBkg[nummolt+1][numPtV0]={0};
  Float_t   EffRelErrSignSB[nummolt+1][numPtV0]={0};
  Float_t   EffRelErrBkgSB[nummolt+1][numPtV0]={0};
  TFile *fileMassSigma;
  TString dirinputtype[6] = {"Xi", "Xi", "Omega", "Omega", "Xi", "Omega"};
  TDirectoryFile *d;
  TString InputDir = "";
  if (isNewInputPath){
    if (!(year.Index("Hybrid")==-1))   {
      InputDir = "MyTask" +dirinputtype[type]+ Form("_MCTruth_PtTrigMin%.1f_PtTrigMax%.1f", PtTrigMin, ptjmax);
    }
    else if (isMC) {
      InputDir = "MyTask" +dirinputtype[type]+ Form("_MCTruth_PtTrigMin%.1f_PtTrigMax%.1f", PtTrigMin, ptjmax);
      if (isGenOnTheFly)  InputDir = Form("MyTask_PtTrigMin%.1f_PtTrigMax%.1f", PtTrigMin, ptjmax);
    }
    else {
      InputDir = "MyTask" +dirinputtype[type]+ Form("_PtTrigMin%.1f_PtTrigMax%.1f", PtTrigMin, ptjmax);
      if (TMath::Abs(PtTrigMin - 0.15) < 0.001 )     InputDir = "MyTask" +dirinputtype[type]+ Form("_PtTrigMin%.1f_PtTrigMax%.1f", 3., 15.);
    }
  } else {
    InputDir = "MyTask"+dirinputtype[type];
  }
  d = (TDirectoryFile*)fin->Get(InputDir);
  cout << "InputDir " << InputDir << endl;
  if (!d)  {cout << "dir input not available " ; return;}

  TString NameContainer ="";
  if (isNewInputPath) {
    //    NameContainer = "_hK0s_Task_RecoAndEfficiency";
    //    NameContainer = "_hXi_Task_name";
    if (year.Index("Hybrid")!=-1)    NameContainer = "_hXi_Task_Hybrid";
    else if (isMC && !isEfficiency) {
      NameContainer = "_hXi_Task_MCTruth";
      if (isGenOnTheFly) NameContainer = "_hK0s_Task_Xi";
    }
    else  {
      if (year=="17pq_hXi")      NameContainer = "_hXi_Task_";
      else    NameContainer = "_hXi_Task_Default";
      if (TMath::Abs(PtTrigMin - 0.15) < 0.001)      NameContainer = "_hXi_Task_LowPtTrig"; 
    }
  }
  TList *list = (TList*)d->Get("MyOutputContainer" + NameContainer);
  cout << "NameContainer " << NameContainer << endl;
  if (!list) {cout << " list does not exist " << endl; return;}

  //identità particelle di trigger nel MC
  TH1F* fHistTriggerPDGCode = new TH1F("fHistTriggerPDGCode", "fHistTriggerPDGCode", 11,0, 11);
  fHistTriggerPDGCode->SetTitle("PDGcode trigger particle");
  Int_t PDGCodeParticles[10] = {211,2212 , 321, 11, 13, 15, 3312, 3334, 3112, 3222}; //sigma- and sigma + respectively
  TString ParticlePDG[11] = {"#pi", "p", "K", "e", "#mu", "#tau", "Xi", "#Omega", "#Sigma^{+}", "#Sigma^{-}", ""};
  for (Int_t i =0; i<11; i++){
    fHistTriggerPDGCode->GetXaxis()->SetBinLabel(i+1,ParticlePDG[i]);
    if (i==10)     fHistTriggerPDGCode->GetXaxis()->SetBinLabel(i+1,"other");
  }

  //calcolo ntrigger **************
  Float_t NTrigger[nummolt+1]={0};
  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
  TH1D *fHistTriggervsMult_MultProj=(TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMin+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(ptjmax-0.00001));

  TH1F*  fHistTriggerXi = new TH1F ("fHistTriggerXi", "fHistTriggerXi", nummoltMax, Nmolt);
  TH1F*  fHistTriggerAll = new TH1F ("fHistTriggerAll", "fHistTriggerAll", nummoltMax, Nmolt);

  for(Int_t m=0; m<nummoltMax+1; m++){
    if(m<nummoltMax){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.0001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.0001); j++ ){
        NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
      fHistTriggerAll->SetBinContent(m+1, NTrigger[m]);
    }
    else {
      for(Int_t j=1; j<=fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
        NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
  }

  TH1F*  fHistTriggerAllButXi = (TH1F*)fHistTriggerAll->Clone("fHistTriggerAllButXi");
  //******************************

  //calcolo numero eventi INT7 in classi di molteplicità
  Float_t CounterINT7Events[nummolt+1]={0};
  TH3F * fHistNumberChargedAllEvents = (TH3F*) list->FindObject("fHistNumberChargedAllEvents");
  if (!fHistNumberChargedAllEvents){cout << " ho histo fHistNumberChargedAllEvents" << endl;  return;}
  TH1F * fHistNumberChargedAllEventsMult = (TH1F*) fHistNumberChargedAllEvents->ProjectionX("fHistNumberChargedAllEventsMolt",0,-1, 0, -1 );

  for (Int_t m=0; m<=nummoltMax; m++){
    if (m<nummoltMax){
      CounterINT7Events[m] =  fHistNumberChargedAllEventsMult->Integral(fHistNumberChargedAllEventsMult->FindBin(Nmolt[m]), fHistNumberChargedAllEventsMult->FindBin(Nmolt[m+1]));
    }
    else {
      CounterINT7Events[m] =  fHistNumberChargedAllEventsMult->Integral(0,100);
    }
  }

  TTree *tSign; 
  TTree *tBkg;
  if (year == "17pq_pp5TeV_hXi_pttrig0.15" || year == "161718Full_AOD234_hXi" || year == "161718_HM_hXi" || year == "161718_HM_hXi_WithFlat16k_No18p" || year == "Train2387_001" || year == "PythiaRopes"){
    tSign = (TTree *)d->Get("MyOutputContainer1"+ NameContainer);
    tBkg  = (TTree *)d->Get("MyOutputContainer2"+ NameContainer);
  }
  else {
    tSign = (TTree *)d->Get("fSignalTree");
    tBkg  = (TTree *)d->Get("fBkgTree");
  }
  
  if (!tSign) {cout << "Signal tree not present" << endl; return;}
  if (!tBkg) {cout << "Bkg tree not present" << endl; return;}

  Float_t LimInfMass[numtipo][nummolt+1][numPtV0]={0};
  Float_t LimSupMass[numtipo][nummolt+1][numPtV0]={0}; 
  Int_t Color[nummolt+1]={1,2,8,4,6,868};


  const Int_t  numPtTrInt=12;
  Float_t NPtTrInt[numPtTrInt+1] ={3,3.5,4,4.5,5,5.5,6,7,8,9,10,12,15};

  const Int_t numQAhisto=12;
  TH1F * fHistQA[numQAhisto];
  TH1F * fHistPtTriggervsPtAssoc[nummolt+1];
  TH1F * fHistNonPrimaryTrigger[nummolt+1];
  TH1F*  fHistSelectedMCTruth[nummolt+1];
  TH1F*  fHistGeneratedMCTruth[nummolt+1];
  TH1F*  fHistEfficiencyMCTruth[nummolt+1];
  TH1F * fHistV0vsMult = new TH1F ("fHistV0vsMult", "NV0 after topological selections and 3sigma mass selection", 100,0,100);
  fHistV0vsMult->GetXaxis()->SetTitle("Multiplicity class");

  TCanvas* canvasQA[numQAhisto];
  TString TitleQAhisto[numQAhisto] = {"Average pT trigger in AC events", "average pT trigger in AC events (each trigger counted only once","mult distribution of AC events (each trigger counted only once)", "NV0/AC event", "fraction of V0 with pT< pT,Trig"/*"Average Pt Trigg vs Pt assoc "*/, "NV0/INT7 event", "Trigger+Xi events / INT7 events"};
  for (Int_t i=0; i<numQAhisto; i++){
    canvasQA[i]= new TCanvas(Form("canvasQA%i",i),TitleQAhisto[i], 800, 500);
    if (i==6) canvasQA[i]->SetTitle("NV0/INT7 events for different mass selections");
    if (i==0 || i==1 || i==5 || i==6) canvasQA[i]->Divide(2,2);
    if (i!=2 && i<numQAhisto){
      fHistQA[i]= new TH1F (Form("fHistQA%i", i), TitleQAhisto[i], nummoltMax,Nmolt ); //multiplicity on the x axis                            
      fHistQA[i]->GetXaxis()->SetTitle("Multiplicity class");
    }
    else {
      fHistQA[i]= new TH1F (Form("fHistQA%i", i), TitleQAhisto[i], 100,0,100);
      fHistQA[i]->GetXaxis()->SetTitle("Multiplicity class");
    }
    //    if (i==4) fHistQA[i]->SetTitle("fraction of V0 with pT< pT,Trig");
  }  

  TString SPtV0[numPtV0]={"0-0.5","0.5-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,0.5, 1,1.5, 2,2.5,3,4,8};

  if (!isEstimateRun3 && NPtV0[numPtV0]!=8){cout << " maximum pt value not set to 8 " << endl;  return;}
  //  if (isEstimateRun3 && NPtV0[numPtV0]!=12){cout << " maximum pt value not set to 12 " << endl;  return;}
  //  TString SPtV0[numPtV0]={"0-0", "0-0.5","0.5-1", "1-1.5", "1.5-2", "2-3", "3-4", "4-8"};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5, 1,1.5,2,3,4,8};
 
  //  TString SPtV0[numPtV0]={"0-0.6", "0.6-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0", "2.0-2.2","2.2-2.5", "2.5-2.9", "2.9-3.4", "3.4-4", "4-5", "5-6.5"};
  //  Double_t NPtV0[numPtV0+1]={0,0.6, 1,1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4,5,6.5};

  Double_t NPtTrigger[numPtTrigger+1]={PtTrigMin,ptjmax};

  for (Int_t m=0; m<nummoltMax+1; m++){
    fHistPtTriggervsPtAssoc[m]=new TH1F (Form("fHistPtTriggervsPtAssoc%i",m), Form("fHistPtTriggervsPtAssoc%i",m), numPtV0, NPtV0);
    fHistNonPrimaryTrigger[m]=new TH1F (Form("fHistNonPrimaryTrigger%i",m), Form("fHistNonPrimaryTrigger%i",m), numPtV0,NPtV0);
    fHistSelectedMCTruth[m]=new TH1F (Form("fHistSelectedMCTruth%i",m), Form("fHistSelectedMCTruth%i",m), numPtV0, NPtV0);
    fHistGeneratedMCTruth[m]=new TH1F (Form("fHistGeneratedMCTruth%i",m), Form("fHistGeneratedMCTruth%i",m), numPtV0, NPtV0);
    fHistEfficiencyMCTruth[m]=new TH1F (Form("fHistEfficiencyMCTruth%i",m), Form("fHistEfficiencyMCTruth%i",m), numPtV0, NPtV0);
  }

  Double_t sigma[numtipo][nummolt+1][numPtV0]={0};
  Double_t mass[numtipo][nummolt+1][numPtV0]={0};
  Double_t nsigmamin[numtipo][nummolt+1][numPtV0]={0};
  Double_t sigmacentral[numtipo][nummolt+1][numPtV0]={0};
  Double_t meansigma=0;
  Double_t meanmass=0;
  TH1F* histoSigma;
  TH1F* histoMean;
  TH1F* histo_ULsideB;
  TH1F* histo_LLsideB;
  TH1F* histo_NSigmaPeak;
  TH1F* histo_NSigmasideB;

  if(isMC==0 || (isMC==1 && isEfficiency==1)){
    for(Int_t m=0; m<nummoltMax+1; m++){
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (MultBinning==1 && isHM && m <= 1) continue;
      //      if (m==0) continue;
      //      PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, syst, PtTrigMin);
      PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab] +Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f", m, sysTrigger, sysV0, syst,PtTrigMinFit);
      if (MultBinning!=0) PathInMassDef+=Form("_MultBinning%i",MultBinning);
      //      PathInMassDef+="_Try1";
      //      PathInMassDef+="_VarRange1";
      PathInMassDef+=".root";
      if (isEstimateRun3) {
	PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab] +Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", 5, sysTrigger, sysV0, syst,PtTrigMinFit);
      }
      if (!IsTrueParticle) {
	fileMassSigma= new TFile(PathInMassDef);
	histoSigma=(TH1F*)fileMassSigma->Get("histo_sigma");
	histoMean=(TH1F*)fileMassSigma->Get("histo_mean");
	histo_ULsideB=(TH1F*)fileMassSigma->Get("histo_ULsideB");
	histo_LLsideB=(TH1F*)fileMassSigma->Get("histo_LLsideB");
	histo_NSigmasideB=(TH1F*)fileMassSigma->Get("histo_NSigmasideB");
	histo_NSigmaPeak=(TH1F*)fileMassSigma->Get("histo_NSigmaPeak");
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
	  sigmacentral[type][m][v]=      histo_NSigmaPeak->GetBinContent(v+1);
	  nsigmamin[type][m][v]=      histo_NSigmasideB->GetBinContent(v+1);
	  mass[type][m][v]=histoMean->GetBinContent(v+1);
	  sigma[type][m][v]=histoSigma->GetBinContent(v+1);
	  LimSupMass[type][m][v]=histo_ULsideB->GetBinContent(v+1);
	  LimInfMass[type][m][v]=histo_LLsideB->GetBinContent(v+1);
	  cout <<"mult interval " <<  m << " PtV0 interval " << SPtV0[v] << " mean " << mass[type][m][v] << " sigma "<< sigma[type][m][v] << " LimSupMass " << LimSupMass[type][m][v] << " liminfmass " << LimInfMass[type][m][v] <<  " nsigmamin " << nsigmamin[type][m][v] << endl;
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
  Int_t        fSignTreeVariableChargeAssoc;
  Int_t        fSignTreeVariableisPrimaryTrigger;
  Int_t        fSignTreeVariableisPrimaryV0;
  Double_t     fSignTreeVariableRapAssoc;
  Double_t     fSignTreeVariableDcaXiDaughters;
  Double_t     fSignTreeVariableV0CosineOfPointingAngle;
  Double_t     fSignTreeVariableXiCosineOfPointingAngle;
  Double_t     fSignTreeVariablectau;
  Double_t     fSignTreeVariablePtV0;
  Double_t     fSignTreeVariableInvMassLambda;
  Double_t     fSignTreeVariableInvMassXi;
  Double_t     fSignTreeVariableInvMassOmega;
  Double_t     fSignTreeVariableEtaV0;
  Double_t     fSignTreeVariablePhiV0;
  Double_t     fSignTreeVariableDeltaEta;
  Double_t     fSignTreeVariableDeltaPhi;
  Double_t     fSignTreeVariableDeltaTheta;
  Double_t     fSignTreeVariableMultiplicity;
  Double_t     fSignTreeVariableZvertex;
  Int_t        fSignTreeVariablePDGCodeTrigger;
  Int_t        fSignTreeVariablePDGCodeAssoc;
  Bool_t        fSignTreeVariableSkipAssoc;


  Double_t     fBkgTreeVariablePtTrigger;
  Int_t        fBkgTreeVariableChargeTrigger;
  Double_t     fBkgTreeVariableEtaTrigger;
  Double_t     fBkgTreeVariablePhiTrigger;
  Double_t     fBkgTreeVariableDCAz;
  Double_t     fBkgTreeVariableDCAxy;
  Int_t        fBkgTreeVariableChargeAssoc;
  Int_t        fBkgTreeVariableisPrimaryTrigger;
  Int_t        fBkgTreeVariableisPrimaryV0;
  Double_t     fBkgTreeVariableRapAssoc;
  Double_t     fBkgTreeVariableDcaXiDaughters;
  Double_t     fBkgTreeVariableV0CosineOfPointingAngle;
  Double_t     fBkgTreeVariableXiCosineOfPointingAngle;
  Double_t     fBkgTreeVariablectau;
  Double_t     fBkgTreeVariablePtV0;
  Double_t     fBkgTreeVariableInvMassLambda;
  Double_t     fBkgTreeVariableInvMassXi;
  Double_t     fBkgTreeVariableInvMassOmega;
  Double_t     fBkgTreeVariableEtaV0;
  Double_t     fBkgTreeVariablePhiV0;
  Double_t     fBkgTreeVariableDeltaEta;
  Double_t     fBkgTreeVariableDeltaPhi;
  Double_t     fBkgTreeVariableDeltaTheta;
  Double_t     fBkgTreeVariableMultiplicity;
  Double_t     fBkgTreeVariableZvertex;
  Int_t        fBkgTreeVariablePDGCodeTrigger;
  Int_t        fBkgTreeVariablePDGCodeAssoc;
  Bool_t        fBkgTreeVariableSkipAssoc;
 
  Float_t CounterNotPrimaryTrigger[nummolt+1][numPtV0] ={0};
  Float_t CounterAllTrigger[nummolt+1][numPtV0] ={0};

  Int_t CounterSignPairsAfterPtMinCut=0;
  Int_t CounterBkgPairsAfterPtMinCut=0;
  Int_t TrueCounterSignPairsAfterPtMinCut=0;
  Int_t TrueCounterBkgPairsAfterPtMinCut=0;
  Int_t CounterSignPairsAfterPtMinCutMult[nummolt+1]={0};
  Int_t CounterSignPairsAfterPtMinCutMult3Sigma[nummolt+1]={0};
  Int_t CounterSignPairsAfterPtMinCutMult4Sigma[nummolt+1]={0};
  Int_t CounterSignPairsAfterPtMinCutMult5Sigma[nummolt+1]={0};
  Int_t CounterSignPairsAfterPtMinCutMultPtTr3[nummolt+1]={0};
  Int_t CounterBkgPairsAfterPtMinCutMult[nummolt+1]={0};
  Int_t TrueCounterSignPairsAfterPtMinCutMult[nummolt+1]={0};
  Int_t TrueCounterBkgPairsAfterPtMinCutMult[nummolt+1]={0};

  //Signal
  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);
  tSign->SetBranchAddress("fTreeVariableChargeAssoc"               ,&fSignTreeVariableChargeAssoc);
  tSign->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger);
  tSign->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0);
  tSign->SetBranchAddress("fTreeVariableRapAssoc"                  ,&fSignTreeVariableRapAssoc);
  tSign->SetBranchAddress("fTreeVariableDcaXiDaughters"            ,&fSignTreeVariableDcaXiDaughters);
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fSignTreeVariableV0CosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariableXiCosineOfPointingAngle"   ,&fSignTreeVariableXiCosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariablectau"                      ,&fSignTreeVariablectau);
  tSign->SetBranchAddress("fTreeVariablePtV0"                      ,&fSignTreeVariablePtV0);
  tSign->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fSignTreeVariableInvMassLambda);
  tSign->SetBranchAddress("fTreeVariableInvMassXi"                 ,&fSignTreeVariableInvMassXi);
  tSign->SetBranchAddress("fTreeVariableInvMassOmega"              ,&fSignTreeVariableInvMassOmega);
  tSign->SetBranchAddress("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0);
  tSign->SetBranchAddress("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0);
  tSign->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta);
  tSign->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi);
  tSign->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta);
  tSign->SetBranchAddress("fTreeVariableMultiplicity"              ,&fSignTreeVariableMultiplicity);
  tSign->SetBranchAddress("fTreeVariableZvertex"                   ,&fSignTreeVariableZvertex);
  tSign->SetBranchAddress("fTreeVariablePDGCodeTrigger"            ,&fSignTreeVariablePDGCodeTrigger);
  tSign->SetBranchAddress("fTreeVariablePDGCodeAssoc"              ,&fSignTreeVariablePDGCodeAssoc);
  tSign->SetBranchAddress("fTreeVariableSkipAssoc"                 ,&fSignTreeVariableSkipAssoc);

  //BackGround                                                                                                                                                                       
  tBkg->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger);
  tBkg->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger);
  tBkg->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger);
  tBkg->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger);
  tBkg->SetBranchAddress("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz);
  tBkg->SetBranchAddress("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy);    
  tBkg->SetBranchAddress("fTreeVariableChargeAssoc"               ,&fBkgTreeVariableChargeAssoc);
  tBkg->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger);
  tBkg->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0);
  tBkg->SetBranchAddress("fTreeVariableRapAssoc"                  ,&fBkgTreeVariableRapAssoc);
  tBkg->SetBranchAddress("fTreeVariableDcaXiDaughters"            ,&fBkgTreeVariableDcaXiDaughters);
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fBkgTreeVariableV0CosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariableXiCosineOfPointingAngle"   ,&fBkgTreeVariableXiCosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariablectau"                      ,&fBkgTreeVariablectau);
  tBkg->SetBranchAddress("fTreeVariablePtV0"                      ,&fBkgTreeVariablePtV0);
  tBkg->SetBranchAddress("fTreeVariableInvMassXi"                 ,&fBkgTreeVariableInvMassXi);
  tBkg->SetBranchAddress("fTreeVariableInvMassOmega"              ,&fBkgTreeVariableInvMassOmega);
  tBkg->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fBkgTreeVariableInvMassLambda);
  tBkg->SetBranchAddress("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0);
  tBkg->SetBranchAddress("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0);
  tBkg->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta);
  tBkg->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi);
  tBkg->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta);
  tBkg->SetBranchAddress("fTreeVariableMultiplicity"              ,&fBkgTreeVariableMultiplicity);
  tBkg->SetBranchAddress("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex);
  tBkg->SetBranchAddress("fTreeVariablePDGCodeTrigger"            ,&fBkgTreeVariablePDGCodeTrigger);
  tBkg->SetBranchAddress("fTreeVariablePDGCodeAssoc"              ,&fBkgTreeVariablePDGCodeAssoc);
  tBkg->SetBranchAddress("fTreeVariableSkipAssoc"                 ,&fBkgTreeVariableSkipAssoc);

  Int_t EntriesSign = 0; 
  Int_t EntriesBkg  = 0; 

  
  TFile *fout = new TFile(PathOut,"RECREATE");
  TDirectory  *dirSign= fout->mkdir("SE");
  TDirectory  *dirBkg= fout->mkdir("ME");

  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsEffw[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsEffwRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbinsEffwErrors[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMass_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebandsEffwRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebandsEffwErrors[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hSign_PtTrigger[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerRatio[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerCountedOnce[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerCountedOnceBis[nummolt+1][numzeta];
  TH1D *hSign_PtTriggerCountedOnceRatio[nummolt+1][numzeta];
  TH1D*  hSign_PtTriggerPtV0bins[nummolt+1][numzeta][numPtV0];
  TH1D *hSign_PtV0[nummolt+1][numzeta];
  Float_t CounterTriggerCountedOnce[nummolt+1][numzeta]={0};
  Float_t CounterTriggerCountedOnce3Sigma[nummolt+1]={0};
  Float_t CounterTriggerCountedOnce4Sigma[nummolt+1]={0};
  Float_t CounterTriggerCountedOnce5Sigma[nummolt+1]={0};
  Float_t CounterV0NotSkipped[nummolt+1]={0};
  Float_t CounterV0NotSkippedPtTr3[nummolt+1]={0};
  Float_t     CounterTriggerAlsoXi=0;
  Float_t     CounterTriggerNoXi=0;

  for(Int_t m=0; m<nummoltMax+1; m++){
    CounterV0NotSkipped[m]=0;
    CounterTriggerCountedOnce3Sigma[m]=0;
    CounterTriggerCountedOnce4Sigma[m]=0;
    CounterTriggerCountedOnce5Sigma[m]=0;
    for(Int_t z=0; z<numzeta; z++){
      CounterTriggerCountedOnce[m][z]=0;
    }
  }


  const   Int_t numDeltaEta=4;
  TString SDeltaEta[numDeltaEta]={"_|DeltaEta|<1.6","_|DeltaEta|<1.4", "_|DeltaEta|<1.1","_|DeltaEta|<0.8" };
  Float_t DeltaEtaLimit[numDeltaEta]={1.6, 1.4, 1.1, 0.8};
  TH1D *hSign_EtaV0[nummolt+1][numDeltaEta];
  TH1D *hSign_EtaV0MidRap[nummolt+1][numDeltaEta];
  TH1D *hSign_RapV0[nummolt+1][numDeltaEta];
  TH2D *hSign_DeltaEtaEtaV0[nummolt+1];
  TH2D *hSign_DeltaEtaEtaTrigger[nummolt+1];

  for(Int_t m=0; m<nummoltMax+1; m++){
    hSign_DeltaEtaEtaV0[m]=new TH2D("hSign_DeltaEtaEtaV0_"+Smolt[m], "hSign_DeltaEtaEtaV0_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    hSign_DeltaEtaEtaTrigger[m]=new TH2D("hSign_DeltaEtaEtaTrigger_"+Smolt[m], "hSign_DeltaEtaEtaTrigger_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    for (Int_t DeltaEta=0; DeltaEta< numDeltaEta; DeltaEta++){
      hSign_EtaV0[m][DeltaEta]=new TH1D("hSign_EtaV0_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_EtaV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hSign_RapV0[m][DeltaEta]=new TH1D("hSign_RapV0_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_RapV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hSign_EtaV0MidRap[m][DeltaEta]=new TH1D("hSign_EtaV0_MidRap_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_EtaV0_MidRap_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
    }
  }

  for(Int_t m=0; m<nummoltMax+1; m++){
    //      if (m==0) continue;
    for(Int_t z=0; z<numzeta; z++){
      hSign_PtTrigger[m][z]=new TH1D("hSign_PtTrigger"+Smolt[m], "hSign_PtTrigger"+Smolt[m], 300,0,30);
      hSign_PtV0[m][z]=new TH1D("hSign_PtV0"+Smolt[m], "hSign_PtV0"+Smolt[m], 300,0,30);
      hSign_PtTrigger[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTrigger[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTrigger[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTrigger[m][z]->GetYaxis()->SetLabelSize(0.05);

      hSign_PtTriggerCountedOnce[m][z]=new TH1D("hSign_PtTriggerCountedOnce"+Smolt[m], "hSign_PtTriggerCountedOnce"+Smolt[m], 300,0,30);
      hSign_PtTriggerCountedOnce[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTriggerCountedOnce[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTriggerCountedOnce[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnce[m][z]->GetYaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnce[m][z]->SetLineColor(Color[m]);
      hSign_PtTriggerCountedOnce[m][z]->SetMarkerColor(Color[m]);
      hSign_PtTriggerCountedOnce[m][z]->SetMarkerStyle(33);

      hSign_PtTriggerCountedOnceBis[m][z]=new TH1D("hSign_PtTriggerCountedOnceBis"+Smolt[m], "hSign_PtTriggerCountedOnceBis"+Smolt[m],numPtTrInt, NPtTrInt);
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnceBis[m][z]->GetYaxis()->SetLabelSize(0.05);
      hSign_PtTriggerCountedOnceBis[m][z]->SetLineColor(Color[m]);
      hSign_PtTriggerCountedOnceBis[m][z]->SetMarkerColor(Color[m]);
      hSign_PtTriggerCountedOnceBis[m][z]->SetMarkerStyle(33);

      hSign_PtV0[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtV0[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtV0[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtV0[m][z]->GetYaxis()->SetLabelSize(0.05);

      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
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
	  hDeltaEtaDeltaPhi_SEbinsRelErr[m][z][v][tr]->SetTitle(nameSE[m][z][v][tr]+ " rel errors");

	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_Effw");
	  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->SetTitle(nameSE[m][z][v][tr]+ " eff corr");

	  hDeltaEtaDeltaPhi_SEbinsEffwRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_EffwRelErr");
	  hDeltaEtaDeltaPhi_SEbinsEffwRelErr[m][z][v][tr]->SetTitle(nameSE[m][z][v][tr]+ " eff corr rel error");

	  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_EffwErrors");
	  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->SetTitle(nameSE[m][z][v][tr]+ " relative error of efficiency");

	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr]+"_SB", nameSE[m][z][v][tr]+"_SB",   56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_SB_Effw");
	  hDeltaEtaDeltaPhi_SEbins_sidebandsEffwRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_SB_EffwRelErr");
	  hDeltaEtaDeltaPhi_SEbins_sidebandsEffwErrors[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]-> Clone(nameSE[m][z][v][tr]+ "_SB_EffwErrors");

	  hInvMass_sidebands[m][z][v][tr]= new TH1D(nameSE[m][z][v][tr]+"_InvMassSidebands", nameSE[m][z][v][tr]+"_InvMassSidebands",  100, 1.29, 1.35);
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
  TH2D *hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins_sidebandsEffwRelErr[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins_sidebandsEffwErrors[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hBkg_PtTrigger[nummolt+1][numzeta];
  TH1D *hBkg_PtV0[nummolt+1][numzeta];

  TH1D *hBkg_EtaV0[nummolt+1][numDeltaEta];
  TH1D *hBkg_EtaV0MidRap[nummolt+1][numDeltaEta];
  TH1D *hBkg_RapV0[nummolt+1][numDeltaEta];
  TH2D *hBkg_DeltaEtaEtaV0[nummolt+1];
  TH2D *hBkg_DeltaEtaEtaTrigger[nummolt+1];


  for(Int_t m=0; m<nummoltMax+1; m++){
    hBkg_DeltaEtaEtaV0[m]=new TH2D("hBkg_DeltaEtaEtaV0_"+Smolt[m], "hBkg_DeltaEtaEtaV0_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    hBkg_DeltaEtaEtaTrigger[m]=new TH2D("hBkg_DeltaEtaEtaTrigger_"+Smolt[m], "hBkg_DeltaEtaEtaTrigger_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    for (Int_t DeltaEta=0; DeltaEta< numDeltaEta; DeltaEta++){
      hBkg_EtaV0[m][DeltaEta]=new TH1D("hBkg_EtaV0_"+Smolt[m]+Form("_%i",DeltaEta), "hBkg_EtaV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hBkg_RapV0[m][DeltaEta]=new TH1D("hBkg_RapV0_"+Smolt[m]+Form("_%i",DeltaEta), "hBkg_RapV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hBkg_EtaV0MidRap[m][DeltaEta]=new TH1D("hBkg_EtaV0_MidRap_"+Smolt[m]+Form("_%i",DeltaEta), "hBkg_EtaV0_MidRap_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
    }
  }

  for(Int_t m=0; m<nummoltMax+1; m++){
    //      if (m==0) continue;
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
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
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
	  hDeltaEtaDeltaPhi_MEbinsRelErr[m][z][v][tr]->SetTitle(nameME[m][z][v][tr]+ " rel error");

	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_Effw");
	  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->SetTitle(nameME[m][z][v][tr]+ " eff corr");

	  hDeltaEtaDeltaPhi_MEbinsEffwRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_EffwRelErr");
	  hDeltaEtaDeltaPhi_MEbinsEffwRelErr[m][z][v][tr]->SetTitle(nameME[m][z][v][tr]+ " eff corr rel error");

	  hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_EffwErrors");
	  hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->SetTitle(nameME[m][z][v][tr]+ " relative error on efficiency");


	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]= new TH2D(nameME[m][z][v][tr]+"_SB", nameME[m][z][v][tr]+"_SB",  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_SB_Effw");
	  hDeltaEtaDeltaPhi_MEbins_sidebandsEffwRelErr[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_SB_EffwRelErr");
	  hDeltaEtaDeltaPhi_MEbins_sidebandsEffwErrors[m][z][v][tr]= (TH2D*) hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]-> Clone(nameME[m][z][v][tr]+ "_SB_EffwErrors");
	}
      }
    }
  }

  EntriesSign =  tSign->GetEntries();
  EntriesBkg  =  tBkg ->GetEntries();
     
  Bool_t BoolVar=kFALSE;
  Bool_t BoolMC=kFALSE;
  Bool_t MassLimit=kFALSE;
  Float_t     fSignTreeVariableInvMassCasc= 0;
  Bool_t isCascTrue=kFALSE;
  Float_t  fSignTreeVariablePtTriggerTemp=0;
  Double_t effSign=0;
  Double_t sigmaEffSign=0;

  dirSign->cd();
  cout << "\n\n I will process " << EntriesSign << " entries for the SE correlation " << endl;
  for(Int_t k = 0; k<EntriesSign; k++){
    //    if (k>100000) continue;
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

    //charge selection                                                                                                                                    
    if ((type==0 || type ==2) && fSignTreeVariableChargeAssoc==1) continue;
    else if ((type==1 || type ==3) && fSignTreeVariableChargeAssoc==-1) continue;

    //inv mass definition                                                                                                                                 
    fSignTreeVariableInvMassCasc= 0;
    if (type==0 || type==1 || type==4)      fSignTreeVariableInvMassCasc= fSignTreeVariableInvMassXi;
    else if (type==2 || type==3 || type==5)      fSignTreeVariableInvMassCasc= fSignTreeVariableInvMassOmega;

    //rapidity selection                                                                                                                                  
    if (israp==0 && TMath::Abs(fSignTreeVariableEtaV0)>0.8)continue;
    else if (israp==1 && TMath::Abs(fSignTreeVariableRapAssoc)>0.5)continue;

    //definition of true cascade                                                                                                                          
    if (type<=3) isCascTrue= (fSignTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==4) isCascTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[0]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[1])  );
    else if (type==5) isCascTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[2]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[3])  );


    if (IsTrueParticle && isMC && isEfficiency==1){
      if (!isCascTrue) continue;
      if (fSignTreeVariableisPrimaryV0!=1) continue;
    }

    //************cuts on pT trigger min*********************************
    if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMin) continue;
    if(TMath::Abs(fSignTreeVariablePtTrigger)>ptjmax) continue;

    //I don't want Xi particles to be trigger particles
    if (isMC && !isEfficiency){

      for(Int_t m=0; m<nummoltMax+1; m++){
	if(m< nummoltMax) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
	else BoolVar=kTRUE;
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
	  if(BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	    fHistGeneratedMCTruth[m]->Fill(fSignTreeVariablePtV0);
	  }
	}
      }

      if (fSignTreeVariablePtV0>3)    CounterTriggerAlsoXi++;
      if(fSignTreeVariablePDGCodeTrigger == PDGCode[0] || fSignTreeVariablePDGCodeTrigger== PDGCode[1]){
	fHistTriggerXi->Fill(fSignTreeVariableMultiplicity);
	continue;
      }
      if (fSignTreeVariablePtV0>3)      CounterTriggerNoXi++;

      for(Int_t m=0; m<nummoltMax+1; m++){
	if(m< nummoltMax) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
	else BoolVar=kTRUE;
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
	  if(BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	    fHistSelectedMCTruth[m]->Fill(fSignTreeVariablePtV0);
	  }
	}
      }
    }

    if(isMC==0 || (isMC==1 && isEfficiency==1)){
      
      //******************* some other cuts for sys studies**************************                  
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
                                                                            
      if (sysV0==0){

	if (fSignTreeVariableDcaXiDaughters> 0.8) continue;
	if (fSignTreeVariableV0CosineOfPointingAngle<0.97)continue;
	if (fSignTreeVariableXiCosineOfPointingAngle < 0.995)continue;
	if (TMath::Abs(fSignTreeVariableInvMassLambda -massLambda) > 0.006)continue;
	if (fSignTreeVariablectau  > 3* ctauCasc[type])continue;

      } 

    }

    for(Int_t m=0; m<nummoltMax+1; m++){
      if(m< nummoltMax) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for(Int_t v=numPtV0Min; v<numPtV0; v++){
	BoolMC =TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))<sigmacentral[type][m][v]*sigma[type][m][v];
	//old choice 	if (isEstimateRun3)   BoolMC = (fSignTreeVariableInvMassCasc > 1.315 && fSignTreeVariableInvMassCasc < 1.330);
	if (isEstimateRun3)   BoolMC = kTRUE; 
	if((isMC && !isEfficiency) || (isMC && isEfficiency && IsTrueParticle)) {
	  BoolMC = kTRUE;
	}
	if(BoolMC && BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	  CounterV0NotSkipped[m] ++;
	  if (fSignTreeVariablePtV0>PtTrigMin)          CounterV0NotSkippedPtTr3[m] ++;
	}
      }
    }
    if(SkipAssoc){    if (fSignTreeVariableSkipAssoc==1) continue;}

    if (isCascTrue)      TrueCounterSignPairsAfterPtMinCut++;

    CounterSignPairsAfterPtMinCut++;

    Int_t counterpdg=0;
    for (Int_t i=0; i<10; i++){
      if (TMath::Abs(fSignTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])	{
	fHistTriggerPDGCode->Fill(i);
	counterpdg++;
      }
    }
    if (counterpdg==0)   fHistTriggerPDGCode->Fill(10);

    Bool_t ispiKpemu=0; //choose species type for trigger particles (!be careful! a selection might have been applied at the level of the task!)
    if (isMC && !isEfficiency){
      if (TriggerPDGChoice==1){
        for (Int_t i=0; i<5; i++){
          if (TMath::Abs(fSignTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])  ispiKpemu = kTRUE;
        }
        if(!ispiKpemu)    continue;
      }
      else if (TriggerPDGChoice==2){
        if(TMath::Abs(fSignTreeVariablePDGCodeTrigger) == 3222 || TMath::Abs(fSignTreeVariablePDGCodeTrigger)== 3112)     continue;
      }
    }

    //**********************************************************************************************
    fSignTreeVariableDeltaPhi = fSignTreeVariablePhiV0-fSignTreeVariablePhiTrigger; 
    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();


    for(Int_t m=0; m<nummoltMax+1; m++){
      if(m< nummoltMax) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
        if(BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
          CounterAllTrigger[m][v]++;
        }
      }
    }

    if (isMC){
      if (fSignTreeVariableisPrimaryTrigger !=1) {
        for(Int_t m=0; m<nummoltMax+1; m++){
          if(m< nummoltMax) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
          else BoolVar=kTRUE;
	  for(Int_t v=numPtV0Min; v<numPtV0; v++){
            if(BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
              CounterNotPrimaryTrigger[m][v]++;
            }
          }
        }
      }
    }


    if (fSignTreeVariableInvMassCasc > 1.315 && fSignTreeVariableInvMassCasc < 1.330){
      fHistV0vsMult->Fill(fSignTreeVariableMultiplicity);
    }

    for(Int_t m=0; m<nummoltMax+1; m++){
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      if (MultBinning==1 && isHM && m <= 1) continue;
      //      if (m==0) continue;
      if(m< nummoltMax) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      hSign_DeltaEtaEtaV0[m]->Fill(fSignTreeVariableEtaV0, fSignTreeVariableDeltaEta);
      hSign_DeltaEtaEtaTrigger[m]->Fill(fSignTreeVariableEtaTrigger, fSignTreeVariableDeltaEta);
      for (Int_t DeltaEta=0; DeltaEta<numDeltaEta; DeltaEta++){
	if (BoolVar && TMath::Abs(fSignTreeVariableDeltaEta) < DeltaEtaLimit[DeltaEta]){
	  hSign_EtaV0[m][DeltaEta]->Fill(fSignTreeVariableEtaV0);
	  hSign_RapV0[m][DeltaEta]->Fill(fSignTreeVariableRapAssoc);
	  if (TMath::Abs(fSignTreeVariableRapAssoc) < 0.5)	    hSign_EtaV0MidRap[m][DeltaEta]->Fill(fSignTreeVariableEtaV0);
	}
      }
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=numPtV0Min; v<numPtV0; v++){
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
	    if((isMC && !isEfficiency) || (isMC && isEfficiency && IsTrueParticle)) {
	      BoolMC = kTRUE;
	      MassLimit=kTRUE;
	    }
	    else {
	      BoolMC =TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))<sigmacentral[type][m][v]*sigma[type][m][v]; 
	      if(isSigma) MassLimit=(TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))<nsigmamax*sigma[type][m][v]);
	      if(!isSigma) {
		MassLimit=(TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && fSignTreeVariableInvMassCasc>LimInfMass[type][m][v] && fSignTreeVariableInvMassCasc<LimSupMass[type][m][v]);
		if (SidebandsSide==1) MassLimit=((mass[type][m][v] - fSignTreeVariableInvMassCasc)>nsigmamin[type][m][v]*sigma[type][m][v] && fSignTreeVariableInvMassCasc>LimInfMass[type][m][v]);
		else if (SidebandsSide==2) MassLimit=((fSignTreeVariableInvMassCasc- mass[type][m][v])>nsigmamin[type][m][v]*sigma[type][m][v] && fSignTreeVariableInvMassCasc<LimSupMass[type][m][v]);
	      }
	      //old choice	      if (isEstimateRun3) BoolMC = (fSignTreeVariableInvMassCasc > 1.315 && fSignTreeVariableInvMassCasc < 1.330);
	      if (isEstimateRun3) BoolMC =kTRUE;
	    }
	    if(BoolMC && BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      if (isEstimateRun3){
		if (TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][5][v]))<3*sigma[type][5][v]) CounterSignPairsAfterPtMinCutMult[m]++;
		if (TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][5][v]))<4*sigma[type][5][v]) CounterSignPairsAfterPtMinCutMult4Sigma[m]++;
		if (TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][5][v]))<5*sigma[type][5][v]) CounterSignPairsAfterPtMinCutMult5Sigma[m]++;
	      }
	      else 	      CounterSignPairsAfterPtMinCutMult[m]++; //number of signal pairs == number of V0 
	      if (fSignTreeVariablePtV0>PtTrigMin) CounterSignPairsAfterPtMinCutMultPtTr3[m]++; 
	      if (isCascTrue)      TrueCounterSignPairsAfterPtMinCutMult[m]++;
	      hSign_PtTrigger[m][z]->Fill(fSignTreeVariablePtTrigger);
	      hSign_PtV0[m][z]->Fill(fSignTreeVariablePtV0);

	      if (fSignTreeVariablePtTrigger!= fSignTreeVariablePtTriggerTemp)       {
		if (isEstimateRun3){
		  if (TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][5][v]))<3*sigma[type][5][v]) CounterTriggerCountedOnce[m][z]++; 
		  if (TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][5][v]))<4*sigma[type][5][v]) CounterTriggerCountedOnce4Sigma[m]++;
		  if (TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][5][v]))<5*sigma[type][5][v]) CounterTriggerCountedOnce5Sigma[m]++;
		}
		else                 CounterTriggerCountedOnce[m][z]++;
		hSign_PtTriggerCountedOnce[m][z]->Fill(fSignTreeVariablePtTrigger);
		hSign_PtTriggerCountedOnceBis[m][z]->Fill(fSignTreeVariablePtTrigger);
		fHistQA[2]->Fill( fSignTreeVariableMultiplicity);
	      }
	    }
	    if(BoolVar && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      effSign =0;
	      if (isEtaEff){
		for (Int_t pt=1; pt<= fHistEfficiencyV0PtEta[m]->GetNbinsX(); pt++){
		  if (fSignTreeVariablePtV0< fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinLowEdge(pt)|| fSignTreeVariablePtV0>= fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinUpEdge(pt)) continue;
		  for (Int_t eta=1; eta<=fHistEfficiencyV0PtEta[m]->GetNbinsY(); eta++){
		    if (fSignTreeVariableEtaV0< fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinLowEdge(eta) || fSignTreeVariableEtaV0>= fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinUpEdge(eta)) continue;
		    effSign = fHistEfficiencyV0PtEta[m]->GetBinContent(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta));
		    sigmaEffSign = fHistEfficiencyV0PtEta[m]->GetBinError(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta));

		  }
		}
	      }

	      if(BoolMC){
		hSign_PtTriggerPtV0bins[m][z][v]->Fill(fSignTreeVariablePtTrigger);
		hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		if (isEtaEff && effSign!=0){
		  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi, 1./effSign);
		  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi, sigmaEffSign/effSign); //rel efficiency
		}
	      }
	
	      if((!isMC || (isMC &&isEfficiency)) && MassLimit && !isEstimateRun3){
		hInvMass_sidebands[m][z][v][tr]->Fill(fSignTreeVariableInvMassCasc);
		hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		if (isEtaEff && effSign!=0){
                  hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi, 1./effSign);
		  hDeltaEtaDeltaPhi_SEbins_sidebandsEffwErrors[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi, sigmaEffSign/effSign); //rel efficiency
                }
	      }
 
	    }
	  }
	}
      }
    }
    fSignTreeVariablePtTriggerTemp = fSignTreeVariablePtTrigger;
  }


  BoolVar=kFALSE;
  BoolMC=kFALSE;
  MassLimit=kFALSE;

  Float_t     fBkgTreeVariableInvMassCasc= 0;
  Double_t effBkg=0;
  Double_t sigmaEffBkg=0;
  cout << "ciao " << endl;
  dirBkg->cd();
  cout << "\n\n I will process " << EntriesBkg << " entries for the ME correlation " << endl;
  if (!isEstimateRun3){
    for(Int_t k = 0; k<EntriesBkg; k++){
      //      if (k>100000) continue;

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
      fBkgTreeVariableDeltaEta=fBkgTreeVariableEtaV0-fBkgTreeVariableEtaTrigger;
      if ((type==0 || type ==2) && fBkgTreeVariableChargeAssoc==1) continue;
      else if ((type==1 || type ==3) && fBkgTreeVariableChargeAssoc==-1) continue;

      //inv mass definition                                                                                                                                 
      fBkgTreeVariableInvMassCasc= 0;
      if (type==0 || type==1 || type==4)      fBkgTreeVariableInvMassCasc= fBkgTreeVariableInvMassXi;
      else if (type==2 || type==3 || type==5)      fBkgTreeVariableInvMassCasc= fBkgTreeVariableInvMassOmega;

      //rapidity selection                                                                                                                                  
      if (israp==0 && TMath::Abs(fBkgTreeVariableEtaV0)>0.8)continue;
      else if (israp==1 && TMath::Abs(fBkgTreeVariableRapAssoc)>0.5)continue;

      //definition of true cascade                                                                                                                         
      if (type<=3) isCascTrue= (fBkgTreeVariablePDGCodeAssoc==PDGCode[type] );
      else if (type==4) isCascTrue= ((fBkgTreeVariablePDGCodeAssoc==PDGCode[0]) ||(fBkgTreeVariablePDGCodeAssoc==PDGCode[1])  );
      else if (type==5) isCascTrue= ((fBkgTreeVariablePDGCodeAssoc==PDGCode[2]) ||(fBkgTreeVariablePDGCodeAssoc==PDGCode[3])  );

      if(SkipAssoc){    if (fBkgTreeVariableSkipAssoc==1) continue;}

      if (IsTrueParticle && isMC && isEfficiency==1){
	if (!isCascTrue) continue;
	if (fBkgTreeVariableisPrimaryV0!=1) continue;
      }

      //I don't want Xi particles to be trigger particles
      if (isMC && !isEfficiency){
	if (fBkgTreeVariablePtV0>3)    CounterTriggerAlsoXi++;
	if(fBkgTreeVariablePDGCodeTrigger == PDGCode[0] || fBkgTreeVariablePDGCodeTrigger== PDGCode[1]) continue;
	if (fBkgTreeVariablePtV0>3)      CounterTriggerNoXi++;

      }

      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMin) continue;
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
	if (sysV0==0){
	  if (fBkgTreeVariableDcaXiDaughters> 0.8) continue;
	  if (fBkgTreeVariableV0CosineOfPointingAngle<0.97)continue;
	  if (fBkgTreeVariableXiCosineOfPointingAngle < 0.995)continue;
	  if (TMath::Abs(fBkgTreeVariableInvMassLambda -massLambda) > 0.006)continue;
	  if (fBkgTreeVariablectau  > 3* ctauCasc[type])continue;
	}
      }

      Bool_t ispiKpemu=0; //choose species type for trigger particles (!be careful! a selection might have been applied at the level of the task!)
      if (isMC && !isEfficiency){
	if (TriggerPDGChoice==1){
	  for (Int_t i=0; i<5; i++){
	    if (TMath::Abs(fBkgTreeVariablePDGCodeTrigger) == PDGCodeParticles[i])  ispiKpemu = kTRUE;
	  }
	  if(!ispiKpemu)    continue;
	}
	else if (TriggerPDGChoice==2){
	  if(TMath::Abs(fBkgTreeVariablePDGCodeTrigger) == 3222 || TMath::Abs(fBkgTreeVariablePDGCodeTrigger)==3112)       continue;

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
      for(Int_t m=0; m<nummoltMax+1; m++){
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue; 
      if (MultBinning==1 && isHM && m <= 1) continue;
	//      if (m==0) continue;
	if(m< nummoltMax) BoolVar =  fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1];
	else BoolVar=kTRUE;
	hBkg_DeltaEtaEtaV0[m]->Fill(fBkgTreeVariableEtaV0, fBkgTreeVariableDeltaEta);
	hBkg_DeltaEtaEtaTrigger[m]->Fill(fBkgTreeVariableEtaTrigger, fBkgTreeVariableDeltaEta);
	for (Int_t DeltaEta=0; DeltaEta<numDeltaEta; DeltaEta++){
	  if (BoolVar && TMath::Abs(fBkgTreeVariableDeltaEta) < DeltaEtaLimit[DeltaEta]){
	    hBkg_EtaV0[m][DeltaEta]->Fill(fBkgTreeVariableEtaV0);
	    hBkg_RapV0[m][DeltaEta]->Fill(fBkgTreeVariableRapAssoc);
	    if (TMath::Abs(fBkgTreeVariableRapAssoc) < 0.5)	    hBkg_EtaV0MidRap[m][DeltaEta]->Fill(fBkgTreeVariableEtaV0);
	  }
	}

	for(Int_t z=0; z<numzeta; z++){
	  for(Int_t tr=0; tr<numPtTrigger; tr++){
	    for(Int_t v=numPtV0Min; v<numPtV0; v++){

	      /*	    if (type==4 || type==5 || type==8){
			    LimInfMass[type]=1.30;
			    LimSupMass[type]=1.342;

			    if (v >4)  {
			    LimInfMass[type]=1.29;
			    LimSupMass[type]=1.352;
			    }
			    }*/
	      if((isMC && !isEfficiency) || (isMC && isEfficiency && IsTrueParticle)) {
		BoolMC = kTRUE;
		MassLimit=kTRUE;
	      }
	      else {
		BoolMC =TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))<sigmacentral[type][m][v]*sigma[type][m][v]; 
		if(isSigma) MassLimit=(TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))<nsigmamax*sigma[type][m][v]);
		if(!isSigma)
		  {
		    MassLimit=(TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin[type][m][v]*sigma[type][m][v] && fBkgTreeVariableInvMassCasc>LimInfMass[type][m][v] && fBkgTreeVariableInvMassCasc<LimSupMass[type][m][v]);
		    if (SidebandsSide==1) MassLimit=(( mass[type][m][v] - fBkgTreeVariableInvMassCasc)>nsigmamin[type][m][v]*sigma[type][m][v] && fBkgTreeVariableInvMassCasc>LimInfMass[type][m][v]);
		    else if (SidebandsSide==2) MassLimit=((fBkgTreeVariableInvMassCasc - mass[type][m][v])>nsigmamin[type][m][v]*sigma[type][m][v] && fBkgTreeVariableInvMassCasc<LimSupMass[type][m][v]);
		  }
	      }	    
 
	      if(BoolMC && BoolVar &&  fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
		hBkg_PtTrigger[m][z]->Fill(fBkgTreeVariablePtTrigger);
		hBkg_PtV0[m][z]->Fill(fBkgTreeVariablePtV0);
	      }	  
	      if(BoolVar && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
		effBkg =0;
		if (isEtaEff){
		  for (Int_t pt=1; pt<= fHistEfficiencyV0PtEta[m]->GetNbinsX(); pt++){
		    if (fBkgTreeVariablePtV0< fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinLowEdge(pt) 	|| fBkgTreeVariablePtV0>= fHistEfficiencyV0PtEta[m]->GetXaxis()->GetBinUpEdge(pt)) continue;
		    for (Int_t eta=1; eta<=fHistEfficiencyV0PtEta[m]->GetNbinsY(); eta++){
		      if (fBkgTreeVariableEtaV0< fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinLowEdge(eta) || fBkgTreeVariableEtaV0>= fHistEfficiencyV0PtEta[m]->GetYaxis()->GetBinUpEdge(eta)) continue;
		      effBkg = fHistEfficiencyV0PtEta[m]->GetBinContent(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta));
		      sigmaEffBkg = fHistEfficiencyV0PtEta[m]->GetBinError(fHistEfficiencyV0PtEta[m]->GetBin(pt, eta));
		    }
		  }
		}

		if(BoolMC){
		  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
		  if (effBkg!=0 && isEtaEff){
		    hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, 1./effBkg);
		    //old		      hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, pow(sigmaEffBkg/pow(effBkg,2),2));
		    hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, sigmaEffBkg/effBkg);
		  }
		}
		if((!isMC || (isMC &&isEfficiency)) && MassLimit){
		  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
		  if(isEtaEff && effBkg!=0){
		    hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, 1./effBkg);
		    hDeltaEtaDeltaPhi_MEbins_sidebandsEffwErrors[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi, sigmaEffBkg/effBkg);
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //put correct errors on 2D histograms                                                           
  cout << "\nput errors on SE/ME histograms " << endl;
  for(Int_t m=0; m<nummoltMax+1; m++){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (MultBinning==1 && isHM && m <= 1) continue;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=numPtV0Min; v<numPtV0; v++){
	  hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_SEbins_sidebandsEffwErrors[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_MEbins_sidebandsEffwErrors[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]);
	  Int_t bin=0;
	  for (Int_t dphi=1; dphi<=hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetNbinsY(); dphi++){
	    for (Int_t deta=1; deta<= hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetNbinsX(); deta++){
	      bin = hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBin(deta, dphi);		
	      if (isEtaEff){
		//		EffRelErr[m][v] = fHistEfficiencyV0PtPtBins[m]->GetBinError(v+1)/fHistEfficiencyV0PtPtBins[m]->GetBinContent(v+1); old
		EffRelErrSign[m][v] = hDeltaEtaDeltaPhi_SEbinsEffwErrors[m][z][v][tr]->GetBinContent(bin);
		EffRelErrSignSB[m][v] = hDeltaEtaDeltaPhi_SEbins_sidebandsEffwErrors[m][z][v][tr]->GetBinContent(bin);

		if (hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin)!=0) {
		  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->SetBinError(bin,sqrt(1./hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(bin) + pow(EffRelErrSign[m][v],2)) * hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		else{
		  hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->SetBinError(bin,0);
		}
		if (hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin)!=0) {
		  hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->SetBinError(bin,sqrt(1./hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetBinContent(bin) + pow(EffRelErrSignSB[m][v],2)) * hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		else{
		  hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->SetBinError(bin,0);
		}

		EffRelErrBkg[m][v] = hDeltaEtaDeltaPhi_MEbinsEffwErrors[m][z][v][tr]->GetBinContent(bin);
		EffRelErrBkgSB[m][v] = hDeltaEtaDeltaPhi_MEbins_sidebandsEffwErrors[m][z][v][tr]->GetBinContent(bin);

		if (hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin) !=0){
		  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->SetBinError(bin,sqrt(1./hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(bin) + pow(EffRelErrBkg[m][v],2)) * hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		//cout << "* m:" << m << " v:" << v << " bin:" << bin <<endl;
		//cout << "rel error due to counts " << 1./hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(bin)<< endl;
		//cout << "rel error due to eff " << EffRelErrBkg[m][v] << endl;
		//cout << " error assigned to corrected ME distribution " << hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinError(bin)/ hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin)<< endl;
		}
		else {
		  hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->SetBinError(bin,0);
		}
		if (hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin) !=0){
		  hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->SetBinError(bin,sqrt(1./hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetBinContent(bin) + pow(EffRelErrBkgSB[m][v],2)) * hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		else {
		  hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->SetBinError(bin,0);
		}

		if (hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin) !=0){
		  hDeltaEtaDeltaPhi_SEbinsEffwRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_SEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		if (hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin) !=0){
		  hDeltaEtaDeltaPhi_MEbinsEffwRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_MEbinsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		if (hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin) !=0){
		  hDeltaEtaDeltaPhi_SEbins_sidebandsEffwRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_SEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin));
		}
		if (hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin) !=0){
		  hDeltaEtaDeltaPhi_MEbins_sidebandsEffwRelErr[m][z][v][tr]->SetBinContent(bin, hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->GetBinError(bin)/hDeltaEtaDeltaPhi_MEbins_sidebandsEffw[m][z][v][tr]->GetBinContent(bin));
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
 
  fHistTriggerAllButXi ->Add(fHistTriggerXi,-1);
  fHistTriggerPDGCode->Scale(1./fHistTriggerPDGCode->GetEntries());
  fout->WriteTObject(fHistTriggerPDGCode);
  fout->WriteTObject(fHistTriggerAllButXi);
  fout->WriteTObject(fHistTriggerXi);
  fout->WriteTObject(fHistTriggerAll);

  //new
  TLegend * legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend * legendLow = new TLegend(0.6, 0.1, 0.9, 0.4);
  for(Int_t z=0; z<numzeta; z++){
    for (Int_t m =0; m<nummoltMax+1; m++){

      hSign_PtTriggerCountedOnceBis[m][z]->Sumw2();
      hSign_PtTrigger[m][z]->Sumw2();

      hSign_PtTriggerCountedOnceBis[5][z]->Sumw2();
      hSign_PtTrigger[5][z]->Sumw2();

      if (m==0){
	hSign_PtTrigger[nummoltMax][z]->Rebin(5);
	hSign_PtTrigger[nummoltMax][z]->Scale(1./ CounterSignPairsAfterPtMinCutMult[nummoltMax]/4);
	hSign_PtTriggerCountedOnceBis[nummoltMax][z]->Scale(1./ CounterTriggerCountedOnce[nummoltMax][z]);
	for (Int_t i = 1; i<= hSign_PtTriggerCountedOnceBis[nummoltMax][z]->GetNbinsX(); i++){
	  hSign_PtTriggerCountedOnceBis[nummoltMax][z]->SetBinContent(i, hSign_PtTriggerCountedOnceBis[nummoltMax][z]->GetBinContent(i)/hSign_PtTriggerCountedOnceBis[nummoltMax][z]->GetBinWidth(i));
	}

      }

      if (m!=nummoltMax){
	hSign_PtTrigger[m][z]->Rebin(5);
	hSign_PtTrigger[m][z]->Scale(1./  CounterSignPairsAfterPtMinCutMult[m]/4);
	hSign_PtTriggerCountedOnceBis[m][z]->Scale(1./ CounterTriggerCountedOnce[m][z]);
	for (Int_t i = 1; i<= hSign_PtTriggerCountedOnceBis[nummoltMax][z]->GetNbinsX(); i++){
	  hSign_PtTriggerCountedOnceBis[m][z]->SetBinContent(i, hSign_PtTriggerCountedOnceBis[m][z]->GetBinContent(i)/hSign_PtTriggerCountedOnceBis[m][z]->GetBinWidth(i));
	}

      }


      hSign_PtTriggerCountedOnceRatio[m][z]=(TH1D*)      hSign_PtTriggerCountedOnceBis[m][z]->Clone("hSign_PtTriggerCountedOnceRatio"+Smolt[m]);

      hSign_PtTriggerCountedOnceRatio[m][z]->Divide(      hSign_PtTriggerCountedOnceBis[nummoltMax][z]);

      hSign_PtTriggerRatio[m][z]=(TH1D*)      hSign_PtTrigger[m][z]->Clone("hSign_PtTriggerRatio"+Smolt[m]);

      hSign_PtTriggerRatio[m][z]->Divide(      hSign_PtTrigger[5][z]);
      //      cout <<  hSign_PtTriggerRatio[m][z] -> GetNbinsX() << " " <<  hSign_PtTrigger[nummoltMax][z]->GetNbinsX() << endl;

      fHistQA[0]->SetBinContent(m+1, hSign_PtTrigger[m][z]->GetMean());
      fHistQA[0]->SetBinError(m+1, hSign_PtTrigger[m][z]->GetMeanError());
      fHistQA[1]->SetBinContent(m+1, hSign_PtTriggerCountedOnce[m][z]->GetMean());
      fHistQA[1]->SetBinError(m+1, hSign_PtTriggerCountedOnce[m][z]->GetMeanError());
      fHistQA[3]->SetBinContent(m+1,  CounterSignPairsAfterPtMinCutMult[m]/CounterTriggerCountedOnce[m][z]);
      //      fHistQA[4]->SetBinContent(m+1, CounterSignPairsAfterPtMinCutMult[m]/ CounterV0NotSkipped[m]);
      fHistQA[4]->SetBinContent(m+1, CounterSignPairsAfterPtMinCutMultPtTr3[m]/ CounterV0NotSkippedPtTr3[m]);
      fHistQA[5]->SetBinContent(m+1,  CounterSignPairsAfterPtMinCutMult[m]/CounterINT7Events[m]);
      fHistQA[5]->SetBinError(m+1, sqrt(1./ CounterSignPairsAfterPtMinCutMult[m] + 1./CounterINT7Events[m])* fHistQA[5]->GetBinContent(m+1));

      fHistQA[6]->SetLineColor(kBlack);
      fHistQA[6]->SetMarkerColor(kBlack);
      fHistQA[6]->SetMarkerStyle(33);
      fHistQA[5]->SetLineColor(kBlack);
      fHistQA[5]->SetMarkerColor(kBlack);
      fHistQA[5]->SetMarkerStyle(33);

      fHistQA[6]->SetBinContent(m+1,  CounterTriggerCountedOnce[m][0]/CounterINT7Events[m]);
      fHistQA[6]->SetBinError(m+1, sqrt(1./ CounterTriggerCountedOnce[m][0] + 1./CounterINT7Events[m])* fHistQA[6]->GetBinContent(m+1));


      if (isEstimateRun3){
	fHistQA[5] ->SetTitle("NV0/INT7 event (3sigma)");
	fHistQA[7] ->SetTitle("NV0/INT7 event (4sigma)");
	fHistQA[9] ->SetTitle("NV0/INT7 event (5sigma)");
	fHistQA[6] ->SetTitle("Trigger+Xi events / INT7 events (3sigma)");
	fHistQA[8] ->SetTitle("Trigger+Xi events / INT7 events (4 sigma)");
	fHistQA[10] ->SetTitle("Trigger+Xi events / INT7 events (5 sigma)");

	fHistQA[7] ->SetLineColor(kBlue);
	fHistQA[8] ->SetLineColor(kBlue);

	fHistQA[9] ->SetLineColor(kRed);
	fHistQA[10] ->SetLineColor(kRed);

	fHistQA[7]->SetBinContent(m+1,  CounterSignPairsAfterPtMinCutMult4Sigma[m]/CounterINT7Events[m]);
	fHistQA[7]->SetBinError(m+1, sqrt(1./ CounterSignPairsAfterPtMinCutMult4Sigma[m] + 1./CounterINT7Events[m])* fHistQA[5]->GetBinContent(m+1));

	fHistQA[8]->SetBinContent(m+1,  CounterTriggerCountedOnce4Sigma[m]/CounterINT7Events[m]);
	fHistQA[8]->SetBinError(m+1, sqrt(1./ CounterTriggerCountedOnce4Sigma[m] + 1./CounterINT7Events[m])* fHistQA[6]->GetBinContent(m+1));

	fHistQA[9]->SetBinContent(m+1,  CounterSignPairsAfterPtMinCutMult5Sigma[m]/CounterINT7Events[m]);
	fHistQA[9]->SetBinError(m+1, sqrt(1./ CounterSignPairsAfterPtMinCutMult5Sigma[m] + 1./CounterINT7Events[m])* fHistQA[5]->GetBinContent(m+1));

	fHistQA[10]->SetBinContent(m+1,  CounterTriggerCountedOnce5Sigma[m]/CounterINT7Events[m]);
	fHistQA[10]->SetBinError(m+1, sqrt(1./ CounterTriggerCountedOnce5Sigma[m] + 1./CounterINT7Events[m])* fHistQA[6]->GetBinContent(m+1));

      }

      canvasQA[0]->cd(1);
      fHistQA[0]->GetYaxis()->SetRangeUser(PtTrigMin+1, PtTrigMin+1.5);
      fHistQA[0]->Draw("");

      canvasQA[0]->cd(2);
      gPad->SetLogy();
      hSign_PtTrigger[m][z]->GetXaxis()->SetRangeUser(PtTrigMin, 6);
      hSign_PtTrigger[m][z]->GetYaxis()->SetRangeUser(0.01, 0.15);
      hSign_PtTrigger[m][z]->Draw("same p");
      if (m==nummoltMax) legend->Draw("");

      canvasQA[0]->cd(4);
      hSign_PtTriggerRatio[m][z]->GetYaxis()->SetRangeUser(0.98, 1.12);
      if (m!=nummoltMax)      hSign_PtTriggerRatio[m][z]->Draw("same l");
      lineat1->Draw("same");

      canvasQA[1]->cd(1);
      fHistQA[1]->GetYaxis()->SetRangeUser(PtTrigMin+1, PtTrigMin+1.5);
      fHistQA[1]->Draw("");

      canvasQA[1]->cd(2);
      gPad->SetLogy();
      legend->AddEntry( hSign_PtTriggerCountedOnceBis[m][z], Smolt[m], "pl");
      legendLow->AddEntry( hSign_PtTriggerCountedOnceBis[m][z], Smolt[m], "pl");
      hSign_PtTriggerCountedOnceBis[m][z]->GetXaxis()->SetRangeUser(PtTrigMin, 15);
      hSign_PtTriggerCountedOnceBis[m][z]->GetYaxis()->SetRangeUser(0.001, 1);
      hSign_PtTriggerCountedOnceBis[m][z]->Draw("same p");
      if (m==nummoltMax) legend->Draw("");

      canvasQA[1]->cd(4);
      hSign_PtTriggerCountedOnceRatio[m][z]->GetYaxis()->SetRangeUser(0.9, 1.1);
      if (m!=nummoltMax)      hSign_PtTriggerCountedOnceRatio[m][z]->Draw("same l");
      lineat1->Draw("same");

      for(Int_t v=numPtV0Min; v<numPtV0; v++){
	fHistPtTriggervsPtAssoc[m]->SetBinContent(v+1,hSign_PtTriggerPtV0bins[m][z][v]->GetMean());
	fHistPtTriggervsPtAssoc[m]->SetBinError(v+1,hSign_PtTriggerPtV0bins[m][z][v]->GetMeanError());
	fHistNonPrimaryTrigger[m]->SetBinContent(v+1, CounterNotPrimaryTrigger[m][v]/CounterAllTrigger[m][v]);
        fHistNonPrimaryTrigger[m]->SetBinError(v+1, SetEfficiencyError(CounterNotPrimaryTrigger[m][v],CounterAllTrigger[m][v]));
      }
      fHistPtTriggervsPtAssoc[m]->GetYaxis()->SetRangeUser(3.5, 7);
      fHistPtTriggervsPtAssoc[m]->SetMarkerColor(Color[m]);
      fHistPtTriggervsPtAssoc[m]->SetLineColor(Color[m]);
      fHistPtTriggervsPtAssoc[m]->SetMarkerStyle(33);
      fHistNonPrimaryTrigger[m]->GetYaxis()->SetRangeUser(0, 0.15);
      fHistNonPrimaryTrigger[m]->SetMarkerColor(Color[m]);
      fHistNonPrimaryTrigger[m]->SetLineColor(Color[m]);
      fHistNonPrimaryTrigger[m]->SetMarkerStyle(33);

      canvasQA[4]->cd();
      fHistPtTriggervsPtAssoc[m]->Draw("same p");
      if (m==nummoltMax) legendLow->Draw("");

      canvasQA[11]->cd();
      fHistNonPrimaryTrigger[m]->Draw("same p");
      if (m==nummoltMax) legendLow->Draw("");
      canvasQA[6] -> SaveAs(PathOutPdf + "_NonPrimaryTriggerFraction.pdf");

      fout->WriteTObject(fHistPtTriggervsPtAssoc[m]); 
      fout->WriteTObject(fHistNonPrimaryTrigger[m]);
      if (isMC && !isEfficiency){
	fout->WriteTObject(fHistGeneratedMCTruth[m]);
	fout->WriteTObject(fHistSelectedMCTruth[m]);
	for (Int_t b=1; b< fHistSelectedMCTruth[m]->GetNbinsX(); b++){
	  fHistEfficiencyMCTruth[m]->SetBinContent(b+1, (float)fHistSelectedMCTruth[m]->GetBinContent(b+1)/fHistGeneratedMCTruth[m]->GetBinContent(b+1));
	  fHistEfficiencyMCTruth[m]->SetBinError(b+1, SetEfficiencyError(fHistSelectedMCTruth[m]->GetBinContent(b+1),fHistGeneratedMCTruth[m]->GetBinContent(b+1)));
	  cout << fHistEfficiencyMCTruth[m]->GetBinContent(b+1) << endl;
	}
	fout->WriteTObject(fHistEfficiencyMCTruth[m]);
      }
    }
    canvasQA[2]->cd();
    fHistQA[2]->Scale(1./fHistQA[2]->GetEntries());
    fHistQA[2]->Draw();
    canvasQA[3]->cd();
    fHistQA[3]->Draw();
    canvasQA[5]->cd(1);
    fHistQA[5]->Draw();
    canvasQA[5]->cd(2);
    fHistQA[6]->Draw();
    canvasQA[5]->cd(3);
    fHistV0vsMult->GetYaxis()->SetRangeUser(0, 1.2 *  fHistV0vsMult->GetMaximum());
    fHistV0vsMult->Draw();
    canvasQA[5]->cd(4);
    fHistNumberChargedAllEventsMult->GetYaxis()->SetRangeUser(0, 1.2 *  fHistNumberChargedAllEventsMult->GetMaximum());
    fHistNumberChargedAllEventsMult->Draw();

    if (isEstimateRun3){
      TLegend * legendsigmaInterval = new TLegend (0.7, 0.7, 0.9, 0.9);
      legendsigmaInterval->AddEntry(fHistQA[5], "|m-#mu|<3#sigma","pl"); 
      legendsigmaInterval->AddEntry(fHistQA[7], "|m-#mu|<4#sigma","pl"); 
      legendsigmaInterval->AddEntry(fHistQA[9], "|m-#mu|<5#sigma","pl"); 

      canvasQA[6]->cd(1);
      fHistQA[5]->Draw("");
      fHistQA[7]->Draw("same");
      fHistQA[9]->Draw("same");
      legendsigmaInterval->Draw("");
      canvasQA[6]->cd(2);
      fHistQA[6]->Draw("");
      fHistQA[8]->Draw("same");
      fHistQA[10]->Draw("same");
      legendsigmaInterval->Draw("");
      canvasQA[6]->cd(4);
      TH1F * fHistQARatio5To3= (TH1F*) fHistQA[10]->Clone("fHistQARatio5to3");
      TH1F * fHistQARatio5To4= (TH1F*) fHistQA[8]->Clone("fHistQARatio5to4");
      fHistQARatio5To3 ->Divide(fHistQA[6]);
      fHistQARatio5To4 ->Divide(fHistQA[6]);
      fHistQARatio5To3->GetYaxis()->SetRangeUser(1,1.2);      
      fHistQARatio5To4->GetYaxis()->SetRangeUser(1,1.2);      
      fHistQARatio5To3 ->Draw("");
      fHistQARatio5To4 ->Draw("same");
      legendsigmaInterval->Draw("");
    }

    for (Int_t i = 0; i<numQAhisto; i++){
      fout->WriteTObject(fHistQA[i]);
      fout->WriteTObject(canvasQA[i]);
    }
  }

  //new 

  if (isMC && !isEfficiency){
    cout << "\nTrigger particles with associated Xi with pT>3: " << endl;
    cout << " are or not Xi  " << CounterTriggerAlsoXi <<endl;
    cout << " are not Xi " << CounterTriggerNoXi << endl;
  }
  fout->Write();
  fout ->Close();

  cout << "partendo dal file " << PathIn << " e " << PathInMassDef << " (per le diverse molteplicità) ho creato il file " << PathOut<< endl;

  if (isEstimateRun3){
    for (Int_t m=0; m< nummoltMax+1; m++){
      cout << "\n m " << Smolt[m] << " NV0/INT7: " << CounterSignPairsAfterPtMinCutMult[m]/CounterINT7Events[m] << " +- " <<  sqrt(1./ CounterSignPairsAfterPtMinCutMult[m] + 1./CounterINT7Events[m])* fHistQA[5]->GetBinContent(m+1) << endl;
      cout << " m " << Smolt[m] << " Trigger-Xi events /INT7 events: " << CounterTriggerCountedOnce[m][0]/CounterINT7Events[m] << " +- " <<  sqrt(1./CounterTriggerCountedOnce[m][0] + 1./CounterINT7Events[m])* fHistQA[5]->GetBinContent(m+1) << " rel.uncert: " << fHistQA[6]->GetBinError(m+1) / fHistQA[6]->GetBinContent(m+1)<< endl;


    }
  }
}

