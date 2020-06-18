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

void AngularCorrelation_firstCasc(Bool_t ishhCorr=1, Bool_t SkipAssoc=1,Float_t ptjmin=3,  Float_t PtTrigMinFit=3, Int_t sysV0=0, Int_t sysTrigger=0,Int_t sys=0,Int_t type=0,  Int_t israp=0,Bool_t isMC=0, Bool_t isEfficiency=1,   TString year0 = "2016",TString year="2016k_MECorr"/*"Run2DataRed_MECorr_hXi""2016k_hK0s_30runs_150MeV"*/, TString yearMC="2018f1_extra_MECorr"/*"2016kl_hXi"/*"2018f1_extra_hK0s_30runs_150MeV"*/,  TString Path1 ="", TString Dir ="FinalOutput/",  Float_t ptjmax=15, Int_t rebin=2,  Int_t rebinx=2,  Int_t rebiny=2, Bool_t MasterThesisAnalysis=0, Bool_t isEnlargedDeltaEtaPhi=0,Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1){ 

  //masterthesis analysis è usato quando devo fare analisi presentata per la tesi
  if (ishhCorr && type!=0) {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  const  Float_t PtTrigMin=ptjmin;
  //  gStyle->SetOptStat(0);

  //lista degli effetti  sistematici studiati in questa macro
  if (sys==3 || sys>5) return;
  if (sysV0>6)return;
  if (ishhCorr && sys<4 && sys!=0) return; //se faccio correlazione hh non posso studiare sistematico associato a scelta regione Sidebands e peak nella distribuzione di massa invariante
  if (sysV0>2 && ishhCorr) return;
  if (sysTrigger!=0) return;
  if (sysV0!=0 && sys!=0) return;

  Int_t sysang=0;
  if (sys==1) sysang=1;
  if (sys==2) sysang=2;

  const Int_t numtipo=10;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
 
  Float_t JetValue=0.74;
  Float_t JetValueDefault=0.74; //questo valore è indipendente dalla scelta del valore di sys ed è il valore per cui viene diviso lo yield nel jet, in modo da essere confrontatato con lo yield OJ e JOJ

  if (ishhCorr || (!ishhCorr && isEnlargedDeltaEtaPhi)) {//rendo ampiezza Jet uguale per hh e hK0s (l'output di hK0s è nel file DeltaEtaPhiEnlarged)
    JetValueDefault =1.;
    if ((sys==0 || sys==5)) JetValue=1.;
    else  if (sys==4) JetValue=0.9;
  }
  else {
    if (sys==4) JetValue=0.9;
  }


  Float_t BulkLowValue=0.74;
  Float_t BulkUpValue=1.1;
  Float_t InclusiveUpValue=1.1;

  if (sys==5 && !ishhCorr){
    BulkLowValue=0.9;
    BulkUpValue=1.1;
  }
  if ((sys==0 || sys==4) && ishhCorr){
    BulkLowValue=1.05;
    BulkUpValue=1.3;
  }
  if (sys==5 && ishhCorr){
    BulkLowValue=1.15;
    BulkUpValue=1.4;
  }

  Dir+="DATA"+year0;
  TString file = year;
  TString file2 = year;
  if (isMC) file += "_MCEff";
  //  file +=Path1; 
  if (ishhCorr) {
    file2+="_hhCorr";
    //    file+="_hhCorr";
    if (isMC)     file2+="_MCEff";
  }
  file2+=Path1;
  //  if(isMC && isEfficiency) file = year + "_MCEff";
  if(isMC && !isEfficiency) file = yearMC + "_MCTruth";
  TString fileMC = yearMC;
  Int_t MC=0;
  if (isMC) MC=1;

  TString PathIn;
  if (isMC && !isEfficiency) PathIn= Dir+"/histo/AngularCorrelation" + file + ".root";
  else{
    PathIn= Dir+"/histo/AngularCorrelation" + file+Path1;
    if(type>=0){
      PathIn +="_";
      if (!ishhCorr)	  PathIn+=tipo[type];
      PathIn +=Srap[israp];
      PathIn +=SSkipAssoc[SkipAssoc];
      if (ishhCorr)      PathIn +="_hhCorr";
    }
    PathIn+= Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f.root", sysTrigger, sysV0, sysang, PtTrigMin);
  }
  if (MasterThesisAnalysis){
    PathIn= Dir+"/histo/AngularCorrelation" + file + Form("_SysT%i_SysV0%i_Sys%i", sysTrigger, sysV0, sysang)+".root";
  }
 
  TString PathInBis =  "FinalOutput/AnalysisResults" + file  + ".root";
  if (ishhCorr) PathInBis="FinalOutput/AnalysisResults" + file2  + ".root"; 


  TString PathInEfficiency = Dir+ "/Efficiency/Efficiency" +fileMC+Path1;// change
  if(type>=0){
    PathInEfficiency +="_";
    if (!ishhCorr)             PathInEfficiency +=tipo[type];
    PathInEfficiency +=Srap[israp];
    PathInEfficiency +=SSkipAssoc[SkipAssoc];
    if (ishhCorr)     PathInEfficiency +="_hhCorr";
  }

  PathInEfficiency+= Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
  //  if (ishhCorr) PathInEfficiency = Dir+ "/Efficiency/Efficiency" +fileMC  +"_hhCorr" + Path1+ Form("_SysT%i_SysV0%i_PtMin%.1f", sysTrigger, sysV0, PtTrigMin)+".root";
  if (MasterThesisAnalysis) PathInEfficiency = Dir+ "/Efficiency/Efficiency" +fileMC  +Form("_MCEff_Efficiency_SysT%i_SysV0%i", sysTrigger, sysV0)+".root";

  TString PathOut1;
  if (isMC && !isEfficiency) PathOut1=Dir+"/histo/AngularCorrelation" + file  +"_Output.root";
  else {
    PathOut1 = Dir+"/histo/AngularCorrelation"  +file+Path1; //+ "_Jet0.75";
    if(type>=0){
      if (!ishhCorr)      PathOut1 +="_"+tipo[type];
      PathOut1 +=Srap[israp];
      PathOut1 +=SSkipAssoc[SkipAssoc];
    }
    PathOut1+=   Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f_Output.root", sysTrigger, sysV0, sys, PtTrigMin);
    if (isEnlargedDeltaEtaPhi)    PathOut1 = Dir+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_DeltaEtaPhiEnlarged_Output.root";
    //    PathOut1 = Dir+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sys, PtTrigMin)+"_Output.root";
    if (MasterThesisAnalysis)    PathOut1 = Dir+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i", sysTrigger, sysV0, sys)+"_Output.root";
  }
 
  TString PathInMass= Dir+"/invmass_distribution_thesis/invmass_distribution";
  if (isMC && isEfficiency) PathInMass+="MCEff_";
  TString Title;
  PathInMass+= Path1; //+"_" ;
  //  PathInMass+= "Eta05_" ; //change
  PathInMass+="_";
  PathInMass+= year;

  cout << "path in: " << PathInBis << endl;
  cout << "path in angular correlation: " << PathIn << endl;
  cout << "path efficiency: " <<PathInEfficiency << endl;
  cout << "pathout: " << PathOut1 << endl;
  cout << "path in mass first part: " << PathInMass << endl;
  TString PathInMassDef;

  TFile *filepurezza;
  TFile *filein = new TFile(PathIn);
  if (!filein) return;
  TFile *fileinbis = new TFile(PathInBis);
  if (!fileinbis) return;
  TFile *fileinEfficiency = new TFile(PathInEfficiency);
  if (!fileinEfficiency) return;
  TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Xi", "Lambda","Xi", "Omega", "Omega", "Xi", "Omega"};
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask"+dirinputtype[type]);
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TList *list2 = (TList*)dir->Get("MyOutputContainer3");

  cout << "\n********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "REMEMBER TO RUN READTREEPLCHIARA_first.C and READTREEPLCHIARA_second.C FIRST! "  << endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=8;
  const Int_t numPtTrigger=1;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
  Int_t Marker[nummolt+1]={7,4,20,22,29,28};
  Int_t Color[nummolt+1]={2,3,4,6,7, 9};
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};

  Float_t ScaleFactor[nummolt+1][numzeta][numPtV0][numPtTrigger]={0};
  Float_t ScaleFactorBulk[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 
  Float_t ScaleFactorJet[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 
  Float_t ScaleFactorJetBulk[nummolt+1][numzeta][numPtV0][numPtTrigger]={0}; 

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  //  TString Smolt[nummolt+1]={"0-0", "0-10", "10-30", "30-50", "50-100", "_all"};
  //  Double_t Nmolt[nummolt+1]={0,0,10,30,50,100}; 

  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  //  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  //  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  //  if (type>0)SPtV0[1]={"0.5-1"};
  if (type>=0)SPtV0[1]={"0.5-1"};
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8};
  //  if (type>0) NPtV0[1]=0.5;
 if (type>=0) NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0"};
  //  if (type>0) SNPtV0[1]={"0.5"};
  if (type>=0) SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }

  TString SPtTrigger[numPtTrigger]={"3-30"};
  Double_t NPtTrigger[numPtTrigger+1]={ptjmin,ptjmax};
  TString Seta[numetabis]={"_eta<0.5", "_allEta"};
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Int_t numbin=2; //questo numero deve essere pari

  Float_t binwx;
  Float_t binwy;
  Bool_t Bool;
 
  TH1F *fHistEtaLimitsOfRegion = new TH1F("fHistEtaLimitsOfRegion", "fHistEtaLimitsOfRegion",6, 0, 6);

  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
  TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj");
  TH1D *HistoTriggerEfficiency     = (TH1D*)fileinEfficiency->Get("HistoTriggerEfficiency"); //trigger efficiency vs mult
  if (!HistoTriggerEfficiency) {cout << " no histotriggerefficiency" << endl; return;}
  TH1D *HistContTriggerMolt        = (TH1D*)fileinEfficiency->Get("HistContTriggerMolt"); //trigger contamination factor vs mult
  if (!HistContTriggerMolt) {cout << " no histoConttriggermolt" << endl; return;}

  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_SEbins_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins_rapMolt[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MErap[nummolt+1][numzeta][numPtV0][numPtTrigger];
  // Float_t ssbFactor[nummolt+1][numPtV0];    
  TH1F * histoSSB[nummolt+1];
  TH1F * histo_Bcentral[nummolt+1];
  TH1F * histo_Bside[nummolt+1];
  TH1F*  fHistPtvsMult_mult[nummolt+1];
  TH1D* HistBI[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D* HistBII[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];

  TH2D *hDeltaEtaDeltaPhi_ME_normbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB][numeta]; //SE/ME norm proiettato in delta Phi
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB][numeta]; //relative errors of histo projected above
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; //relative errors of histos above
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[nummolt+1][numzeta][numPtV0][numPtTrigger]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[nummolt+1][numzeta][numPtV0][numPtTrigger]; 
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata  + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[nummolt+1][numzeta][numPtV0][numPtTrigger]; //relative error of histo above
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata  + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EtaAll[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis]; //istogramma finale su cui effettuo fit

  Float_t norm_MEbins[nummolt+1][numPtV0][numSB]={0}; //valore ME isto in (0,0)
  Float_t norm_MEbins_norm[nummolt+1][numPtV0][numSB]={0}; //valore ME isto normalizzato in (0,0)
  Int_t   NTrigger[nummolt+1]={0}; //total number of trigger particles 
  Int_t   NTriggerV0[nummolt+1]={0}; //total number of trigger particles only in events with V > 0
  Float_t NSpectrum[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles
  Float_t NSpectrumError[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles
  Float_t NSpectrumV0[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles in events with V> 0
  TH1D*  fHistYvsMult[numetabis];

  TH1D*  fHistV0EfficiencyPtBins[nummolt+1]; //efficienza selezione V0 in PtV0 bins (vedi SPtV0)
  TH1D*  HistContV0PtBins[nummolt+1]; //contamination factor V0 in PtV0 bins (vedi SPtV0)

  TF1* 	  gauss[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis];
  TF1* 	  gaussint[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis];
  TF1*    pol0ZYAM[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TF1*    pol0Bulk[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TF1*    pol0BulkBis[nummolt+1][numzeta][numPtV0][numPtTrigger];

  Bool_t isProjectionPhi[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString Proj[3]={Form("|#Delta#eta| < %.1f, ",JetValue), Form("%.1f < |#Delta#eta| < %.1f, ",BulkLowValue, BulkUpValue),Form("|#Delta#eta| < %.1f, ",InclusiveUpValue)};

  TString TitleString[nummolt+1][numPtV0][2];
  TString TitleStringBis[nummolt+1][numPtV0];
  TString TitleStringTris[numPtV0];


  //********************************************************************* 
  //**************calcolo numero particelle di trigger*******************
  //********************************************************************* 
  for(Int_t m=0; m<nummolt+1; m++){
    //    if (m==0) continue;
    fHistV0EfficiencyPtBins[m]= (TH1D*)fileinEfficiency->Get("fHistV0EfficiencyPtBins_" + Smolt[m]);
    if (!fHistV0EfficiencyPtBins[m]) {cout << "no V0EfficiencyPtBins in file " << PathInEfficiency<<endl; return;}; 
    HistContV0PtBins[m]= (TH1D*)fileinEfficiency->Get("HistContV0PtBins_" + Smolt[m]);
    if (!HistContV0PtBins[m]) {cout << "no HistContV0PtBins " << PathInEfficiency<<endl; return;}; 
    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<=fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    cout << "\nn trigger in mult range (all triggers) [the one which should be used]    " << m << "  " <<  NTrigger[m] <<   endl;
  }
    
  for(Int_t sb=0; sb< numSB; sb++){
    if(sb==1 && isMC && !isEfficiency) continue;
    if(sb==1 && ishhCorr) continue;
    //    for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t m=nummolt; m>=0; m--){
    //if (m==0) continue;
      if (!ishhCorr){
	if(!isMC || (isMC && isEfficiency)){

	  PathInMassDef=PathInMass+ "_"+tipo[type]+Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab] +Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, sysang,PtTrigMinFit);

	  //	  if(type==0)	  PathInMassDef=PathInMass+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, sysang, PtTrigMinFit);
	  if (MasterThesisAnalysis) 	  PathInMassDef=PathInMass+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i.root", m, sysTrigger, sysV0, sysang);
	  cout << "******************************" << endl;
	  cout << "******************************" << endl;
	  cout << "\npath in mass def completo " << PathInMassDef << endl;
	  filepurezza= new TFile(PathInMassDef);
	  histoSSB[m]=(TH1F*)filepurezza->Get("histo_SSB");
	  histo_Bcentral[m]=(TH1F*)filepurezza->Get("histo_Bcentral");
	  histo_Bside[m]=(TH1F*)filepurezza->Get("histo_Bside");
	
	}
      }
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=1; v<numPtV0; v++){
	    for (Int_t sb=0; sb<2; sb++){
	      TitleString[m][v][sb]= "in mult. " + Smolt[m] +"% and " + SNPtV0[v] + " GeV/c < p_{T}^{Assoc} < " + SNPtV0[v+1]+" GeV/c"+ Ssidebandtitle[sb];
	    }
	    TitleStringBis[m][v]= "in mult. " + Smolt[m] +"% and " + SNPtV0[v] + " GeV/c < p_{T}^{Assoc} < " + SNPtV0[v+1]+" GeV/c";
	    TitleStringTris[v]= "in " + SNPtV0[v] + " GeV/c < p_{T}^{Assoc} < " + SNPtV0[v+1]+" GeV/c";

	    cout << "\n\n *********************************************" << endl;
	    cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v]<< " per sideband(1) or not(0) = " << sb << endl;
	  
	    nameME[m][z][v][tr][sb]="ME_";
	    nameME[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];
	    nameSE[m][z][v][tr][sb]="SE_";
	    nameSE[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];

	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameSE[m][z][v][tr][sb]);
	    if (!hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]){cout << "missing histo: " << nameSE[m][z][v][tr][sb] << endl;  return;}
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameME[m][z][v][tr][sb]);
	    if (!hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]) {cout << "missing histo: " << nameME[m][z][v][tr][sb] << endl;  return;}
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Rebin2D(rebinx,rebiny);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Rebin2D(rebinx,rebiny);

	    //	    cout << "comparison of errors " << endl; //I obtain same results
	    hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]= (TH2D*)	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Clone(nameSE[m][z][v][tr][sb]+"_RelError");
	    hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]= (TH2D*)	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb]+"_RelError");

	    for (Int_t j=1; j<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsX(); j++){
	      for (Int_t l=1; l<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsY(); l++){
	      Int_t Binuniv=  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBin(j,l);
	      hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]->SetBinContent(Binuniv,         hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinError(Binuniv)/hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinContent(Binuniv));
	      hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]->SetBinContent(Binuniv,         hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinError(Binuniv)/  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(Binuniv));
	      hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]->SetBinError(Binuniv,0);
	      hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]->SetBinError(Binuniv,0);

	      //	      cout << j << "  " << l << "  " << Binuniv << endl;
	      //	      cout <<    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinError( Binuniv) << " vs " << sqrt(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinContent(Binuniv))<<  endl;
	      }
	      }

	    //	    cout << "\nho preso istogrammi Me e Se" << endl;

	    /////////// primo modo per normalizzare ME 
	    binwx= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->GetBinWidth(1);
	    binwy= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->GetBinWidth(1);

	    cout << "\nafter the rebin: bin x (delta phi) " << binwx <<  " bin y (delta eta) " << binwy << endl;
	    Int_t cont=0;
	    for(Int_t i=0; i < numbin ; i++ ){
	      for(Int_t j=0; j < numbin ; j++ ){
		norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001+(i-numbin/2)*binwx), hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001+(j-numbin/2)*binwy));
		cont++;
	      }
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/cont;  	    
  
	    cout << "\nvalore histo in DeltaEta=0 e DeltaPhi=0 (quadrato) " <<  norm_MEbins[m][v][sb] << endl;

	    /////////// secondo modo per normalizzare ME 
	    norm_MEbins[m][v][sb]=0;
	    cont=0;
	    for(Int_t i=0; i <hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetNbinsY()  ; i++ ){
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001),i+1);
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-0.000001),i+1);
	      cont++;
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/2/cont;
	    
	    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 (media in delta phi, 2 righe) " <<  norm_MEbins[m][v][sb] << endl;

	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_norm");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC");
	    if (norm_MEbins[m][v][sb]==0){
	      cout << "sb " << sb << endl;
	      cout << "normalizzazione ME non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	      //norm continue;  
	    }

	    //	    cout << "\n divido SE distribution per ME disatribution normalized" << endl;
	    if (norm_MEbins[m][v][sb]!=0) hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]-> Scale(1./norm_MEbins[m][v][sb]);
	    //norm	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]);
	    //****************provo a utilizzare pair acceptance della classe di molteplicità 0-100% per tutte le molt************
	    cout << "I divide " << endl;
	    if (!ishhCorr)	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]);
	    if (ishhCorr)	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]);
	    cout << "I have divided " << endl; //error are propagated assuming non correlated histograms
	    hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC_RelError");

	    for (Int_t j=1; j<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsX(); j++){
	      for (Int_t l=1; l<=hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetNbinsY(); l++){
	      Int_t Binuniv=  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBin(j,l);
	      //	      cout <<	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinContent(Binuniv) << " +-  " << 	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetBinError(Binuniv) << endl;
	      //	      cout <<	      hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]->GetBinContent(Binuniv) << " +-  " << 	      hDeltaEtaDeltaPhi_ME_normbins[5][z][v][tr][0]->GetBinError(Binuniv) << endl;
	      //	      cout <<	      hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinContent(Binuniv) << " +-  " << 	      hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinError(Binuniv) << endl;
	      hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]->SetBinContent(Binuniv,         hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinError(Binuniv)/  hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetBinContent(Binuniv));
	      hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]->SetBinError(Binuniv,0);
	      }
	      }


	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->SetTitle("Raw SE distribution "+ TitleString[m][v][sb]);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->SetTitle("ME distribution "+ TitleString[m][v][sb]);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->SetTitle("Normalized ME distribution "+ TitleString[m][v][sb]);
    	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->SetTitle("SE corrected for pair acceptance "+TitleString[m][v][sb]);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->SetTitleOffset(1.3);
	    norm_MEbins_norm[m][v][sb]= hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001), hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001));
	    cout << "valore histo  in DeltaEta=0 e DeltaPhi=0, histo normalizzato  " <<  norm_MEbins_norm[m][v][sb] <<"\n"<< endl;

	    //histograms of relative errors of SE, ME and AC (not yet here)


	    //*****************************************************************************************************************
	    //proietto in Delta Phi
	    //****************************************************************************************************************

	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][0]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaJet", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-JetValue)+1, hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(JetValue)-1, "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][2]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaAll",  hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-InclusiveUpValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(InclusiveUpValue), "E"));

	    HistBI[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaBulkI", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-BulkUpValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-BulkLowValue), "E"));
	    HistBII[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaBulkII", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(BulkLowValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(BulkUpValue), "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]= (TH1D*)	    HistBI[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb]+"_AC_phi_etaBI");
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->Add(	    HistBII[m][z][v][tr][sb]);

	    /*//the following gives consistent result to the Add() method used just above; the errors given by root after Add(A+B) method is the sqrt(errA^2 + errB^2)
	      for(Int_t j=1; j < HistBI[m][z][v][tr][sb]->GetNbinsX(); j++){
		cout << "error before " <<  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->GetBinError(j) << endl;
		//	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->SetBinContent(j, HistBI[m][z][v][tr][sb]->GetBinContent(j)+ HistBII[m][z][v][tr][sb]->GetBinContent(j));
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->SetBinError(j, sqrt(pow(HistBI[m][z][v][tr][sb]->GetBinError(j),2)+ pow(HistBII[m][z][v][tr][sb]->GetBinError(j),2)));	    
	      cout << "error after " <<  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->GetBinError(j) << endl;
	      }
	    */


	    for(Int_t eta=0; eta<numeta; eta++){
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetTitle("#Delta#phi projection in " +Proj[eta] + TitleString[m][v][sb]);
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetXaxis()->SetTitle("#Delta#phi (radians)");
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetYaxis()->SetTitle("Counts");
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetXaxis()->SetTitleOffset(1);
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetYaxis()->SetTitleOffset(1);
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetYaxis()->SetTitleSize(0.05);
	      if (sb==1)	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetLineColor(kRed);
	      Bool_t NotScale=kFALSE;

	      if (!ishhCorr){
		if (sb==1){
		  if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) NotScale=kTRUE;
		  if(!NotScale){
		    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->Scale((Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001))/histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001)));
		  }
		  else {
		    cout << " I could not scale the SB distribution due to empty sidebands " << endl;
		  }
		  //		  cout << "Bcentral " << (Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001)) << " Bside " << histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001)) << " B central/Bside " << (Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001))/histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))<< endl;
		}
	      }

	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]= (TH1D*)	 	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->Clone(nameME[m][z][v][tr][sb]+Form("_AC_phi_eta%i_RelError",eta));
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetTitle("Relative errors of #Delta#phi projection in " +Proj[eta] + TitleString[m][v][sb]);
	      for(Int_t j=1; j < hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->GetNbinsX(); j++){
		hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->SetBinContent(j,TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetBinError(j)/hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetBinContent(j)) );
		hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->SetBinError(j,0);
	      }
	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->GetYaxis()->SetRangeUser(0,1);
	      hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->SetLineColor(kBlack);
	    }//chiusura ciclo eta
	  } //chiusura ciclo su pt della V0
	}//chiusura ciclo su Pt Trigger
      }//chiusura ciclo su z vertice
    }//chiusura ciclo su molteplcità
  }//chiusura ciclo su sideband or not

  cout << "********************************************" << endl;
  cout << "**                                       **" << endl;
  cout << "  **                                   **" << endl;
  cout << "    **                               **" << endl;
  cout << "      **                           **" << endl;
  cout << "        **                       **" << endl;
  cout << "          **                    **" << endl;
  cout << "            **                **" << endl;
  cout << "              **            **" << endl;
  cout << " \n\n ********** second part of the program *******************" << endl;
  for(Int_t m=0; m<nummolt+1; m++){
    //if (m==0) continue;
    // for(Int_t m=nummolt; m>=0; m--){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=1; v<numPtV0; v++){
	  if (!ishhCorr){
	    cout << "\n\n *********************************************" << endl;
	    cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v] << endl;
	    if(!isMC || (isMC && isEfficiency)){
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][0]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_rap");
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][1], hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][0]);
	      Title="ME sidebands/ME peak "+TitleStringBis[m][v];
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->SetTitle(Title);
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->GetXaxis()->SetTitle("#Delta#eta");
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->GetXaxis()->SetTitleOffset(1.3);
	      hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->GetYaxis()->SetTitle("#Delta#phi (radians)");

	    }
	  }

	  if((!isMC || (isMC && isEfficiency)) && !ishhCorr) Bool= (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0);
	  if((!isMC || (isMC && isEfficiency)) && ishhCorr) 	   Bool= (norm_MEbins[m][v][0]==0);
	  else  Bool= (norm_MEbins[m][v][0]==0);
	  if (Bool){
	    cout << "normalizzazione ME (SB o peak region) non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	    //norm	    continue;  
	  }

	  //***************************************************************
	  //sottraggo fondo data da "finte" K0s se analizzo hV0
	  //***************************************************************
	  for(Int_t eta=0; eta< numeta; eta++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta] =  (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub", eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta] =  (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub_RelError", eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->SetName(nameME[m][z][v][tr][0]+Form("_eta%i_AC_phi_V0Sub",eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->SetName(nameME[m][z][v][tr][0]+Form("_eta%i_AC_phi_V0Sub_RelError",eta));

	    if (isMC && !isEfficiency) continue;
	    if (ishhCorr) continue;

	    //	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Add(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta],-1); 	    //this for the correct fake removal

	    //error check: both methods below give the same results 
	    /*
	      for(Int_t i=1; i<hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetNbinsX(); i++ ){
	      cout << " bin " << i << ": " << 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinContent(i)<< " +- " <<     hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinError(i)<<endl;
	      cout << "bin in old way " << 	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->GetBinContent(i) << " +- "<<  sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->GetBinError(i),2) +pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->GetBinError(i),2)) << endl;
	      }
	    */
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Scale(histoSSB[m]->GetBinContent(v+1));//this for a fast fake removal
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetXaxis()->SetTitle("#Delta#phi (radians)");

	    for (Int_t j=1; j<= hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetNbinsX() ; j++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->SetBinContent(j, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinError(j)/  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetBinContent(j));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->SetBinError(j,0);
	    }
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->GetYaxis()->SetRangeUser(0,1);

	  }

	  //***************************************************************
	  //sottraggo distribuzione del bulk 
	  //***************************************************************
	  ScaleFactor[m][z][v][tr]=(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue)))/2./(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)));
	  //	  cout << (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue))) << "   " <<(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)))<<"  " << (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue))) << "  " << hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue));
	  //	  cout << "\n" <<ScaleFactor << endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactor[m][z][v][tr]);

	  //**********here I rebin the two histograms before subtraction*************** 
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub_Rebin", 1));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub_Rebin", 0));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->Rebin(rebin);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Rebin(rebin);
	  //***************************************************************************

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->SetLineColor(kGreen+2);

	  //norebin	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSubFit_EffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr");

	  //no rebin 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Add(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1],-1);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Add(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1],-1);

	  //fit of the bulk in the delta phi jet region to subtract bulk
	  pol0Bulk[m][z][v][tr]= new TF1("pol0",Form("pol0Bulk_m%i_v%i",m,v), -1, 1);
	  pol0Bulk[m][z][v][tr]->SetLineColor(kGreen+2);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->Fit(  pol0Bulk[m][z][v][tr], "R+");
	  for (Int_t b=1; b< hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetNbinsX(); b++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetBinContent(b, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinContent(b)-  pol0Bulk[m][z][v][tr]->GetParameter(0));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetBinError(b, sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinError(b),2)+ pow(pol0Bulk[m][z][v][tr]->GetParError(0),2)));
	    }
	  //	  error check //!!!!!ROOT computes errors assuming uncorrelated histograms 
/* 
	    for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetNbinsX(); i++){
	    cout << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinContent(i)<< endl;
	    cout << "norm factor value (should be zero) " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetNormFactor()<< endl;
	    //	    cout << "error squared unchecked " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinErrorSqUnchecked(i) << " error squared " <<  pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2)<< endl;
	    cout << "error of the first histogram of the sum " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinContent(i)<< " +- " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i)<< endl;
	    cout << "error of the second histogram of the sum " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinContent(i)<< " +- " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i)<< endl;
	    cout << "errors calculated by root: " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< " +- " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinError(i)<< " rel: " << hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinError(i)/hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	    cout << "errors assuming zero correlation: " << sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2) +  pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2)) << " rel error = " << sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2) +  pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2))/ 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	    cout << "errors assuming total correlation (ro =1): " << TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i)- hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i)) << " rel error = " << TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i)- hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i)) / hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	    //the following formula is not correct	    cout << "errors assuming 'partial' correlation  " <<  sqrt(TMath::Abs(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2))) << " rel error = " <<  sqrt(TMath::Abs(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinError(i),2))) / 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)<< endl;
	    cout << "sqrt(bin content) of BulkSub histo " << sqrt(TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)))<< " rel error " <<1./ sqrt(TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(i)))<<endl;
	    }
*/

	    for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetNbinsX(); i++){
	      //	    cout << " I set errors " << endl;
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->SetBinError(i, TMath::Abs(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->GetBinError(i)- hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->GetBinError(i)));
	    }
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr] = (TH1D*) 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorrNotScaled");

	  //riscalo per deltaEta width tutte e tre le proiezioni in delta phi (jet, bulk, inclusive) e clono histo. L'histo clonato sarà diviso per efficienza

	  //*****bulk*********
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(1./ScaleFactor[m][z][v][tr]);
	  ScaleFactorBulk[m][z][v][tr]=(1./2/(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)))); 
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactorBulk[m][z][v][tr]);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr");	 
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetLineColor(kGreen+2);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetTitle("UE distribution in " +Proj[1] + TitleStringBis[m][v]);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_Bulk_EffCorr_RelError");	 

	  //********jet*****
	  ScaleFactorJet[m][z][v][tr]=(1./(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValueDefault)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValueDefault))));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->Scale(ScaleFactorJet[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Scale(ScaleFactorJet[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Scale(ScaleFactorJet[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Scale(ScaleFactorJet[m][z][v][tr]);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]= (TH1D*)	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_RelError");
	  for (Int_t j=1; j<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->GetNbinsX() ; j++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->SetBinContent(j,hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinError(j)/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->GetBinContent(j));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->SetBinError(j,0);
	  }
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1);

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetLineColor(kBlue);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetLineColor(kGreen+3);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetLineColor(kRed);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->SetLineColor(kRed);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetTitle("Jet distribution in " +Proj[0] + TitleStringBis[m][v]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetTitle("Jet distribution in " +Proj[0] + TitleStringBis[m][v] + "Not scaled by #Delta#eta");

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorr_RelError");

	  //********inclusive*******
	  ScaleFactorJetBulk[m][z][v][tr]=1./2/(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(InclusiveUpValue)));
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Scale(ScaleFactorJetBulk[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulkEffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulkEffCorr_RelError");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->SetTitle("Jet + UE distribution in " +Proj[0] + TitleStringBis[m][v]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_JetBulkEffCorrNotScaled");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Scale(1./ScaleFactorJetBulk[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->SetTitle("Jet + UE distribution in " +Proj[0] + TitleStringBis[m][v]);


	  //****************************************************************************************************************
	  //divido angular correlation proietatta in deltaphi per efficienza di selezione Trigger e V0 e per contamination factors e per SSB
	  //****************************************************************************************************************
	  cout << PathInEfficiency << endl;
	  cout << "\n \n Trigger selection efficiency: " <<HistoTriggerEfficiency->GetBinContent(m+1) << endl;
	  cout << "\n V0 selection efficiency:  " << fHistV0EfficiencyPtBins[m]->GetBinContent(v+1) << endl;
	  cout << "\n Trigger contamination factor:  " << (HistContTriggerMolt->GetBinContent(m+1)) << endl;
	  cout << "\n V0 contamination factor:  " << (HistContV0PtBins[m]->GetBinContent(v+1))<< endl;

	  if (!isMC || (isMC && isEfficiency)){
	    if((fHistV0EfficiencyPtBins[m]->GetBinContent(v+1))!=0 ){

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));

	      //setto errore all'istogramma scalato 
	      Float_t Error_BulkSub_EffCorr=0;
	      Float_t Error_BulkSub_EffCorrNotScaled=0;
	      Float_t Error_BulkSubFit_EffCorr=0;
	      Float_t Error_EffCorr=0;
	      Float_t Error_Bulk_EffCorr=0;
	      Float_t Error_JetBulk_EffCorr=0;
	      Float_t Error_JetBulk_EffCorrNotScaled=0;

	      for(Int_t i=1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++ ){
	      Error_BulkSub_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

	      Error_BulkSubFit_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

	      Error_BulkSub_EffCorrNotScaled= 	      Error_BulkSub_EffCorr*ScaleFactorJet[m][z][v][tr];

	      Error_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

	      Error_Bulk_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

	      Error_JetBulk_EffCorr= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));

	      Error_JetBulk_EffCorrNotScaled= sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->GetBinError(i),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->GetBinContent(i)/(1-HistContV0PtBins[m]->GetBinContent(v+1))*HistContV0PtBins[m]->GetBinError(v+1),2) + pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->GetBinContent(i)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*fHistV0EfficiencyPtBins[m]->GetBinError(v+1),2));


	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetBinError(i, Error_BulkSub_EffCorr);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetBinError(i, Error_BulkSubFit_EffCorr);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetBinError(i, Error_BulkSub_EffCorrNotScaled);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->SetBinError(i, Error_EffCorr);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetBinError(i, Error_Bulk_EffCorr);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->SetBinError(i, Error_JetBulk_EffCorr);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->SetBinError(i, Error_JetBulk_EffCorrNotScaled);

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->SetBinContent(i, Error_BulkSub_EffCorr/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->GetBinContent(i));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->SetBinContent(i, Error_Bulk_EffCorr/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->GetBinContent(i));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->SetBinContent(i, Error_JetBulk_EffCorr/ hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->GetBinContent(i));
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->SetBinError(i,0);

	      }

	      //fit for ZYAM*******
	      pol0ZYAM[m][z][v][tr]= new TF1("pol0",Form("pol0ZYAM_m%i_v%i",m,v), 1,2);
	      pol0ZYAM[m][z][v][tr]->SetLineColor(kRed);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Fit(pol0ZYAM[m][z][v][tr], "R+");

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]= (TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr_BulkSubZYAM");
	      for(Int_t i=1; i<	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++){
	      //	    cout << " I set errors " << endl;
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->SetBinContent(i, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)-pol0ZYAM[m][z][v][tr]->GetParameter(0));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->SetBinError(i,sqrt( pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinError(i),2)+pow(pol0ZYAM[m][z][v][tr]->GetParError(0),2)));
	
	    }

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->Rebin(rebin);
	      //*******************

	      //fit for BulkSubFitBis*******

	      pol0BulkBis[m][z][v][tr]= new TF1("pol0",Form("pol0BulkBis_m%i_v%i",m,v), -1,1);
	      pol0BulkBis[m][z][v][tr]->SetLineColor(kGreen+3);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Fit(pol0BulkBis[m][z][v][tr], "R+");

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]= (TH1D*)	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EffCorr_BulkSubBulkFit");
	      for(Int_t i=1; i<	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetNbinsX(); i++){
	      //	    cout << " I set errors " << endl;
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetBinContent(i, hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinContent(i)-pol0BulkBis[m][z][v][tr]->GetParameter(0));
		hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetBinError(i,sqrt( pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->GetBinError(i),2)+pow(pol0BulkBis[m][z][v][tr]->GetParError(0),2)));
	
	    }

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->Rebin(rebin);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetLineColor(kGreen+3);
	      //*******************

	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]   ->GetYaxis()->SetRangeUser(0,1);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1);
	      isProjectionPhi[m][z][v][tr]=kTRUE;
	    }
	    else {
	      cout << "\n ****** efficienza pari a zero********* " << endl;
	      isProjectionPhi[m][z][v][tr]=kFALSE;
	      continue;
	    }
	    /*
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->Scale(1./NTrigger[m]);
	      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Scale(1./NTrigger[m]);
	    */
	  }

	  //la seguente parte era utilizzata per sottrarre il contributo dell'UE nella regione di jet, utilizzando il valore medio di alcuni bin simmetricamente disposti intorno a zero. Si è scelto di cambiare metodo. Per riutilizzarlo, alcuni nomi di histo vanno cambiati. 
	  //****************************************************************************************************************
	  //bkg subtraction
	  //****************************************************************************************************************
	  /*
	    for(Int_t etabis=0; etabis<2; etabis++){
	    Float_t DeltaPhiBkg=0;
	    Float_t counter=0;
	    Float_t DeltaPhiBkgError=0;

	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->FindBin(-1.1); i++){
	    counter++;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i);
	    }
	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.84); i++){
	    counter++;
	    //	    cout <<" content " <<  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i) << endl;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i);
	    }
	    DeltaPhiBkg=DeltaPhiBkg/counter;

	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-1.1); i++){
	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i)-DeltaPhiBkg,2);
	    }
	    for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.84); i++){//simmetrico attorno a pi/mezzi

	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i)-DeltaPhiBkg,2);
	    }
	    DeltaPhiBkgError=sqrt(DeltaPhiBkgError/counter/(counter+1));


	    cout << "numero di bin su cui si media per sottrarre fondo: " << counter << "\n valore del fondo di coppie non correlate: " << DeltaPhiBkg << " +- " << DeltaPhiBkgError<<  endl;
	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]= (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_BkgSub"+Seta[etabis]);
	    //hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Sumw2();

	    if(etabis==1){
	    for(Int_t i =1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->SetBinContent(i, (Float_t)(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i))-DeltaPhiBkg);
	    //the errors below are exactly the same
	    //	    cout << "\nbin n: " << i << endl;
	    //	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinError(i)<< endl;
	    //	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->GetBinError(i)<< endl;
	    }//
	    }


	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Rebin(rebin);
	    } //etabis
	  */
	}//v
      }//tr
    }//z
  }//m


  cout << "\n\n\n Studio la dipendenza della pair correction dalla classe di molteplicità"<< endl;
  for(Int_t sb=0; sb< numSB; sb++){
    if(sb==1 && isMC && !isEfficiency && !ishhCorr) continue;
    if(sb==1 && ishhCorr) continue;
    for(Int_t m=0; m<nummolt+1; m++){
    //if (m==0) continue;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=1; v<numPtV0; v++){

	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb] + "_rapMolt");
	    if (!	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]) return;
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->SetName(nameME[m][z][v][tr][sb] + "_rapMolt");
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->SetTitle("ME " +Smolt[m]+"% / ME 0-100%"+TitleStringTris[v] );
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta#eta");
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta#phi (radians)");
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->Divide( hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_MEbins[5][z][v][tr][sb]);
	  }
	}
      }
    }
  }

 
  auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->SetHeader("Selezioni applicate");     


  TCanvas *canvasJetFit[nummolt+1];
  for(Int_t m=0; m< nummolt+1; m++){
    //if (m==0) continue;
    canvasJetFit[m]= new TCanvas (Form("canvasJet_%i",m),Form("canvasJet_%i",m) , 1200, 1000);
    canvasJetFit[m]->Divide(4,2);
    for(Int_t v=1; v< numPtV0; v++){
      //      cout << "v" << v << endl;
      canvasJetFit[m]->cd(v+1);

      if((!isMC || (isMC && isEfficiency)) && !ishhCorr) Bool= (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0);
      if((!isMC || (isMC && isEfficiency)) && ishhCorr) Bool= (norm_MEbins[m][v][0]==0);
      else Bool= (norm_MEbins[m][v][0]==0);
    
      if (Bool) {
	cout << " normalization not performed " << endl;	   
	//norm 	continue;
      }
      if(!isProjectionPhi[m][0][v][0]){
	cout << "projection not performed " << endl;
	continue;
      }	   
      hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][0][v][0]->DrawCopy("");
    
    }
  }

  cout << "\n sto per scrivere su file " << endl;
  // salvo solo istogrammi con norm != 0 sia per sideband che per regione centrale, con efficienze diverse da zero e con integrale siiideband non =0 (ossia solo histo con significato fisico) 
  TFile *fileout = new TFile(PathOut1, "RECREATE");
  for(Int_t m=0; m<nummolt+1; m++){
    //if (m==0) continue;
    cout << "\n" <<m << endl;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=1; v<numPtV0; v++){
	  cout << "v " << v << endl;
	  if(!isMC || (isMC && isEfficiency)){

	    if (!ishhCorr){
	      if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) {
		cout << " B side =0 " << endl;
		continue;
	      }
	    }
	    
	    if (!ishhCorr)	     hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Write();
	  }
	  
	  if((!isMC || (isMC && isEfficiency)) && !ishhCorr) Bool= (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0);
	  if((!isMC || (isMC && isEfficiency)) && ishhCorr) Bool= (norm_MEbins[m][v][0]==0);

	  if (Bool) {
	    cout << " pt interval" << v << endl;
	    cout << " normalization not performed " << endl;	   
	    //norm 	    continue;
	   
	  }
	  if(!isProjectionPhi[m][z][v][tr]){
	    cout << "projection not performed " << endl;
	    continue;
	  }	   
	  for(Int_t sb=0; sb<numSB; sb++){
	    if(sb==1 && isMC && !isEfficiency && !ishhCorr) continue;
	    if(sb==1 && ishhCorr) continue;
		
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_MEbins_Error[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_SEbins_Error[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins_Error[m][z][v][tr][sb]->Write();
	    for (Int_t eta=0; eta<=2; eta++){
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->Write();
	    hDeltaEtaDeltaPhi_ACbins_phi_Error[m][z][v][tr][sb][eta]->Write();
	    }
	  }
	  //decommentare per avere proiezione in out-of-jet region sovrapposto a proiezione in jet-region
	  // hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(1./ScaleFactorBulk[m][z][v][tr]);
	  // hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactor[m][z][v][tr]);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->SetTitle("");
	  for (Int_t eta=0; eta<=2; eta++){
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Error[m][z][v][tr][eta]->Write();
	  }
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->SetTitle("Jet region projection corrected for fake "+tipoTitle[type] + ", efficiency and DeltaEta width");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->SetTitle("Jet region projection corrected for fake "+ tipoTitle[type]+", bulk subtracted and DeltaEta width");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->SetTitle("Raltive errors of Jet region projection corrected for fake "+ tipoTitle[type]+", bulk subtracted and DeltaEta width");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_Error[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->SetTitle("Jet region projection bulk-sub, (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr_Error[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->SetTitle("Jet region projection bulk-sub (bulk fit), (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][0]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Rebin[m][z][v][tr][1]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSubFit_EffCorr[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->SetTitle("Jet region projection bulk-sub (ZYAM), (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSub[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->SetTitle("Jet region projection bulk-sub (bulk fit), (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EffCorr_BulkSubBulkFit[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->SetTitle("Jet region projection bulk-sub, (fake-"+tipoTitle[type]+" + Eff) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorrNotScaled[m][z][v][tr]->Write();


	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->SetTitle("OJ region projection, (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_Bulk_EffCorr_Error[m][z][v][tr]->Write();

	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->SetTitle("JOJ region projection, (fake-"+tipoTitle[type]+" + Eff + DeltaEta width) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->SetTitle("JOJ region projection, (fake-"+tipoTitle[type]+" + Eff) corrected");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorr_Error[m][z][v][tr]->Write();
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_JetBulk_EffCorrNotScaled[m][z][v][tr]->Write();
	}
      }
    }
  }

  TCanvas *canvasDraw[nummolt+1][numzeta][numPtTrigger][numPtV0][5];
  for(Int_t m=0; m<nummolt+1; m++){
    //if (m==0) continue;
    cout << "\n" <<m << endl;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=1; v<numPtV0; v++){
	  cout << "la normalizzazione del ME non è stato effettuata per: " << endl;
	  if (norm_MEbins[m][v][0]==0) cout << " m " << m << " PtV0 " << v << " sideband 0"<< endl;
	  if (norm_MEbins[m][v][1]==0) cout << " m " << m << " PtV0 " << v << " sideband 1"<< endl;
 
	  if (!ishhCorr){
	    if(!isMC || (isMC && isEfficiency)){
	      if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) {
		cout << " B side =0 " << endl;
		continue;
	      
	      }
	    }
	    canvasDraw[m][z][v][tr][0]=new TCanvas (Form("CanvasDraw%i_%i_%i_%i_0", m,z,v,tr),Form("CanvasDraw%i_%i_%i_%i_0", m,z,v,tr), 800, 600 );
	    canvasDraw[m][z][v][tr][0]->cd();
	    if (ishhCorr) 	    hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Draw();
	      
	    if (m==0 && v ==1) canvasDraw[m][z][v][tr][0]->SaveAs("MESBRatio05.pdf");
	    canvasDraw[m][z][v][tr][0]->Close();	      
	    //fileout->WriteTObject(canvasDraw[m][z][v][tr][0]);
	  }
	}

      }
    }
  }


  fHistEtaLimitsOfRegion->SetBinContent(1,hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(-JetValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(1, "JetLowValue");
  fHistEtaLimitsOfRegion->SetBinContent(2,hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(JetValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(2, "JetUpValue");
  fHistEtaLimitsOfRegion->SetBinContent(3,hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(BulkLowValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(3, "BulkLowValue");
  fHistEtaLimitsOfRegion->SetBinContent(4,hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(BulkUpValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(4, "BulkUpValue");
  fHistEtaLimitsOfRegion->SetBinContent(5,hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(-InclusiveUpValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(5, "InclusiveLowValue");
  fHistEtaLimitsOfRegion->SetBinContent(6,hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(InclusiveUpValue)));
  fHistEtaLimitsOfRegion->GetXaxis()->SetBinLabel(6, "InclusiveUpValue");

  fHistEtaLimitsOfRegion->Write();

  fileout->Close(); 
  
  cout << "******************************************************************"<< endl;
  cout << "****** sto utilizzando pair acceptance della classe molt 0-100%*****, per utilizzare pair acceptance delle diverse classi, cambiare righe commentate con //norm" << endl;

  cout << "\npartendo dai file " << PathIn << " e "<< PathInBis << " e "<< PathInEfficiency <<  " ho creato: "<< endl;
  cout << "\nil file " << PathOut1 << endl;
  if (ishhCorr) cout << " use sysV0 = 0,1,2 and sys==0 " << endl;
  cout << "\nbin x (delta phi) " << binwx <<  " bin y (delta eta) " << binwy << endl;
  if (MasterThesisAnalysis) cout << " *****************************Be aware master thesis analysis has been selected!!" << endl;

  cout << "DeltaEta Jet interval " << hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(-JetValue)) << " - " << hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(JetValue)) << endl;

  cout << "DeltaEta Bulk interval " << hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(BulkLowValue)) << " - " << hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(BulkUpValue)) << endl;

  cout << "DeltaEta Inclusive interval " << hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(-InclusiveUpValue)) << " - " << hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[1][0][1][0][0]->GetXaxis()->FindBin(InclusiveUpValue)) << endl;

  cout << "!!!!!!!!! se valori si sovrappongono, modificare intervallo deltaEta" << endl;

}
