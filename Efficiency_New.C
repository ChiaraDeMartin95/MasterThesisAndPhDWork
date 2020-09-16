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

Double_t SetEfficiencyError(Int_t k, Int_t n){
  return sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
}

void Efficiency_New(Bool_t CommonParton=0, Bool_t isCPEff=0, Bool_t ishhCorr =0, Int_t type=6, Int_t israp=0, TString ResoHisto="2D",Bool_t SkipAssoc=1 ,   Float_t ptjmin=3,  Float_t ptjmax=15, Int_t sysTrigger=0, Int_t sysV0=0, TString data="1617GP_hXi" /*"2018f1_extra_MECorr"/"2018f1_extra_hK0s_CP"/* "2018f1_extra_hK0s_CP"/*"2018f1_extra_hK0s_CP_10runs_Bis"/"AllMC_hXi"*/, TString year0="2016", TString path1=""/*"_New"*/ /*"_10runs_FB128_TrackLengthCut"/"_NewMultClassBis"/"_PtTrigMax2.5"*/, Int_t MultBinning=0, Int_t PtBinning=1){

  //  if (type==6 && israp!=0) return; //so far I haven't appended any info on the rapidity to efficiency output files for K0s; if I decide to do as ofr the cascade, just remove this line

  if (sysTrigger> 3) return;
  if (!ishhCorr){
    if (sysV0> 6) return;
  }
  if (ishhCorr){
    if (sysV0> 2) return;
  }
  if (type>7){
    cout << "type can only be 0=XiNeg, 1=XiPos, 2=OmegaNeg, 3=OmegaPos, 4=Xi, 5=Omega , 6=K0s, 7=h " << endl;
    return;
  }
  if (ishhCorr && type!=7) {cout << " if you want to study hh correlation, you have to choose type=7 " << endl; return;}

  // ResHisto è pari a 2D se gli istogrammi della risoluzione sono 2D (Delta vs PtTrig, min)

  //2018f1_extra_onlyTriggerWithHighestPt"
  const   Float_t PtTrigMin=ptjmin;
  const   Float_t PtTrigMax=ptjmax;
  const Int_t numtipo = 8; 

  TString tipo[numtipo]={"XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega", "K0s", "h"};
  TString tipoPart[numtipo]={"Xi", "Xi", "Omega", "Omega", "Xi", "Omega", "", ""};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};

  TString file = "Efficiency";
  TString PathInBis =  "FinalOutput/AnalysisResults";
  file+=data;
  //  file+="_MCEff";
  PathInBis+=data;
  if (!ishhCorr) {
    PathInBis+="_MCEff";
    //    PathInBis+=path1;
    PathInBis+=".root";
  }
  if (ishhCorr) {

    PathInBis+="_hhCorr_MCEff";
    //    PathInBis+=path1;
    //    PathInBis+= "_MCEff.root";
    PathInBis+= ".root";
  }
  cout << PathInBis << endl;
  //  TString PathInBis =  "AnalysisResults.root";
  TString PathOut2="FinalOutput/DATA" + year0 + "/Efficiency/" + file + path1;
  if (PtBinning>0) PathOut2 += Form("_PtBinning%i",PtBinning);
  if (CommonParton && isCPEff)  PathOut2 += "_CPEff";
  if (CommonParton && !isCPEff)  PathOut2 += "_NOCPEff";
  PathOut2 += "_"+tipo[type]+Srap[israp]+Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
  if (!ishhCorr && !SkipAssoc)PathOut2="FinalOutput/DATA" + year0 + "/Efficiency/" + file + path1+  Form("_PtBinning%i",PtBinning)+ "_"+tipo[type]+Srap[israp]+Form("_AllAssoc_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
  if (ishhCorr) PathOut2="FinalOutput/DATA" + year0 + "/Efficiency/" + file  + path1 + Form("_PtBinning%i",PtBinning)+"_"+Srap[israp]+ "_hhCorr"+Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
  if (ishhCorr &&!SkipAssoc) PathOut2="FinalOutput/DATA" + year0 + "/Efficiency/" + file + "_hhCorr" + path1 + Form("_PtBinning%i",PtBinning)+Srap[israp]+ "_hhCorr"+Form("_AllAssoc_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
  cout << "new input file " << endl;

  cout << PathInBis << endl;
  TFile *fileinbis = new TFile(PathInBis, "");

  if (!fileinbis){
    cout << "input file does not exist " << endl; 
    return;
  }

  TString PathInSel="./FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";
  PathInSel+=data;
  PathInSel+="_MCEff";
  PathInSel+=path1;
  if (PtBinning>0)  PathInSel+=  Form("_PtBinning%i",PtBinning);
  PathInSel+="_";
  if (!ishhCorr) {
    PathInSel +=tipo[type];
  }
  PathInSel +=Srap[israp];
  if (!SkipAssoc)  PathInSel +="_AllAssoc";
  if (ishhCorr) PathInSel+="_hhCorr";
  PathInSel +=Form("_MassDistr_SysT%i_SysV0%i_PtMin%.1f",sysTrigger, sysV0, PtTrigMin);
  PathInSel+= ".root";
  cout << "fileinputSel :" << PathInSel << endl;

  TFile *fileinputSel = new TFile(PathInSel, "");
  if (!fileinputSel){
    cout << "input sel file does not exist; run readTreePLChiaraCasc_first.C first " << endl; 
    return;
  }
  cout << "fileinputSel :" << PathInSel << endl;

  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask"+ tipoPart[type]);
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TList *list2 = (TList*)dir->Get("MyOutputContainer3"); //contiene info V0
  TList *list3 = (TList*)dir->Get("MyOutputContainer4"); //contiene info trigger particle
  TList *listRisoluzione = (TList*)dir->Get("Risoluzione"); //contiene info risoluzione nei task più recenti
  if (!listRisoluzione) {
    cout << "non c'è la lista degli istogrammi che si dovrebbe chiamare Risoluzione, contolla in che lista sono gli istogrammi della risoluzione " << endl;
    cout << "Ed elimina lo studio del TProfile dalla macro: se non c'è lista non c'è neanche TProfile "<< endl;
    //    return; 
  }
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9; //14;//was 8
  const Int_t numPtTrigger=19;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  Int_t sysV0Gen=0;
  if (israp==0) sysV0Gen=0;
  else if (israp==1) sysV0Gen=1;



  Float_t ptjminGen=0;
  Float_t ptjmaxGen=0;
    Float_t ptjminBisFixed[nummolt+1]={0.150,1,2, 3,4,5};
    Float_t ptjminBis[nummolt+1];
    Float_t ptjmaxRes[nummolt+1]={ptjmax};

    cout << " I will print ... " << endl;
  if (type==0 || type ==2) {
    ptjminGen = -ptjmax; ptjmaxGen=-ptjmin;
    for (Int_t i =0; i< nummolt+1; i++){    ptjminBis[i]=-ptjmax;}
    for (Int_t i =0; i< nummolt+1; i++){    ptjmaxRes[i]=-ptjminBisFixed[i];}
    cout << ptjminBis[0] << "  " << ptjminBis[3] << endl;    
    cout << ptjmaxRes[0] << "  " << ptjmaxRes[3] << endl;    
}
  else  if (type==1 || type ==3 || type >5) {
    ptjminGen = ptjmin;  ptjmaxGen=ptjmax;
    for (Int_t i =0; i< nummolt+1; i++){    ptjminBis[i]=ptjminBisFixed[i];}
    for (Int_t i =0; i< nummolt+1; i++){    ptjmaxRes[i]=ptjmax;}
}
  else if (type==4 || type==5) {
    ptjminGen = -ptjmax;  ptjmaxGen=ptjmax;
    for (Int_t i =0; i< nummolt+1; i++){    ptjminBis[i]=-ptjmax;}
    for (Int_t i =0; i< nummolt+1; i++){    ptjmaxRes[i]=ptjmax;}

}

  
  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 

  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 

  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 

  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  
  for (Int_t m =0; m< nummolt+1; m++){  
    if (MultBinning==0){
    Smolt[m] = Smolt0[m];
    SmoltLegend[m] = SmoltLegend0[m];
    Nmolt[m] = Nmolt0[m];
    }
    else     if (MultBinning==1){
    Smolt[m] = Smolt1[m];
    SmoltLegend[m] = SmoltLegend1[m];
    Nmolt[m] = Nmolt1[m];
    }
    else    if (MultBinning==2){
    Smolt[m] = Smolt2[m];
    SmoltLegend[m] = SmoltLegend2[m];
    Nmolt[m] = Nmolt2[m];
    }
  }

  TString Szeta[numzeta]={""};

  /*
  TString SPtV0[numPtV0]={"0-0.6", "0.6-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0", "2.0-2.2", "2.2-2.5", "2.5-2.9", "2.9-3.4", "3.4-4", "4-5", "5-6.5"};
  Double_t NPtV0[numPtV0+1]={0,0.6, 1,1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4,5,6.5};
  */
  
  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  //TString SPtV0[numPtV0]={"", "", "0.5-1", "1-1.5","1.5-2","2-3", "3-4", "4-8"};
  if (type>=0)SPtV0[1]={"0.5-1"};
  if (type==6) {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8,100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  //el 
  if (type>=0) NPtV0[1]=0.5;
  if (type==6) SPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};
  if (type>=0) SNPtV0[1]={"0.5"};
  if (type==6)  SNPtV0[1]="0.5";
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }
  
  TString SPtV01[numPtV0]={ "0.1-0.5", "0.5-0.8", "0.8-1.2","1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8,1.2,1.6,2,2.5,3,4,8};
  TString SNPtV01[numPtV0+1]={"0.1","0.5","0.8","1.2","1.6","2.0","2.5","3.0","4.0","8.0"};
  if (PtBinning==1){
    for(Int_t j=0; j<numPtV0+1; j++){
      cout << " j " << j << endl;
      if (j<numPtV0)      SPtV0[j] = SPtV01[j];
      NPtV0[j] = NPtV01[j];
      SNPtV0[j] = SNPtV01[j];
    }
  }
  Int_t numPtV0Max =numPtV0;
  if (PtBinning==0) numPtV0Max =numPtV0-1;

  Double_t NPtTrigger[numPtTrigger+1]={0,0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 13.0, 20.0, 30.0};
  //  Int_t Marker[nummolt+1]={7,4,20,22,29, 35};
  Int_t Marker[nummolt+1]={7,20,20,22,29,25};
  //  Int_t Color[nummolt+1]={1,2,418,4,909,868};
  Int_t Color[nummolt+1]={1,2,8,4,6,868};
  Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  Int_t ColorSysV0[numSysV0]={2,419,4, kAzure+8,8,881};
 
  TCanvas *canvasEff=new TCanvas ("canvasEff", "canvasEff", 1500, 800);
  TCanvas *canvasEffBis[6];
  TString nomecanvaseff[6]={"TPt", "Tphi", "Teta", "V0Pt", "V0phi", "V0eta"};
  /*  
      for(Int_t j=0; j<6; j++){
      canvasEffBis[j]=new TCanvas (Form("canvasEffBis%i",j), "canvasEff"+nomecanvaseff[j], 800, 600);
      }
  */
  TCanvas *canvasRes=new TCanvas ("canvasRes", "canvasRes", 1500, 800);
  TCanvas *canvasCont=new TCanvas ("canvasCont", "canvasCont", 1000, 700);
  TCanvas *canvasUsed=new TCanvas ("canvasUsed", "canvasUsed", 1000, 700);
  // TCanvas *canvascontv0=new TCanvas ("canvascontv0", "canvascontv0", 800, 600);
 
  canvasEff->Divide(3,2);
  canvasRes->Divide(3,2);
  canvasCont->Divide(2,1);
  canvasUsed->Divide(2,2);
 
 
  //I exchange ordering to set option DCAz < 1 as default one ******************: not necessary anymore
  Int_t sysV0hh= sysV0;
  /*
  if (ishhCorr){
    if (sysV0 == 0) sysV0hh =0;
    else if (sysV0 == 1) sysV0hh =0;
    else  sysV0hh =2;
  }
  */
  //***************istogrammi delle efficienze**********************************************************
  TH3D*   fHistGeneratedTriggerPtPhi= (TH3D*)list3->FindObject("fHistGeneratedTriggerPtPhi");
  if (!fHistGeneratedTriggerPtPhi) {cout << "histo gen not found " << endl; return;}
  TH3D*   fHistGeneratedTriggerPtEta=     (TH3D*)list3->FindObject("fHistGeneratedTriggerPtEta");
  if (!fHistGeneratedTriggerPtEta) {cout << "histo gen not found " << endl; return;}

  //  TH3D*   fHistSelectedTriggerPtPhi=  (TH3D*)fileinputSel->Get(Form("fHistSelectedTriggerPtPhi_%i",sysTrigger));
  TH3D*   fHistSelectedTriggerPtPhi=  (TH3D*)list3->FindObject(Form("fHistSelectedTriggerPtPhi_%i",sysTrigger));
  if (!fHistSelectedTriggerPtPhi) {cout << "histo sel not found " << endl; return;}
  TH3D*   fHistSelectedTriggerPtEta=  (TH3D*)list3->FindObject(Form("fHistSelectedTriggerPtEta_%i",sysTrigger));
  if (!fHistSelectedTriggerPtEta) {cout << "histo sel not found " << endl; return;}

  TH3D*   fHistSelectedGenTriggerPtPhi;
  if ((TH3D*)list3->FindObject(Form("fHistSelectedGenTriggerPtPhi_%i", sysTrigger))) {
      fHistSelectedGenTriggerPtPhi=  (TH3D*)list3->FindObject(Form("fHistSelectedGenTriggerPtPhi_%i", sysTrigger));
    }
    else {
      cout << " fHistSelectedGenTriggerPtPhi non esiste " << endl;
      fHistSelectedGenTriggerPtPhi=   (TH3D*) fHistSelectedTriggerPtPhi->Clone(Form("fHistSelectedGenTriggerPtPhi", sysTrigger));
    }

  TH3D*   fHistGeneratedV0PtTMaxPhi;
  TH3D*   fHistGeneratedV0PtTMaxPhiInt1;
  TH3D*   fHistGeneratedV0PtTMaxPhiInt2;

  if (CommonParton && isCPEff){
    fHistGeneratedV0PtTMaxPhi=      (TH3D*)list2->FindObject(Form("fHistCPGeneratedV0PtTMaxPhi_%i", sysV0Gen));
  }
  else if (CommonParton && !isCPEff){
    fHistGeneratedV0PtTMaxPhiInt1=       (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtTMaxPhi_%i", sysV0Gen));
    fHistGeneratedV0PtTMaxPhiInt2=       (TH3D*)list2->FindObject(Form("fHistCPGeneratedV0PtTMaxPhi_%i", sysV0Gen));
    if (fHistGeneratedV0PtTMaxPhiInt1)    fHistGeneratedV0PtTMaxPhi= (TH3D*)     fHistGeneratedV0PtTMaxPhiInt1->Clone(Form("fHistNOCPGeneratedV0PtTMaxPhi_%i", sysV0Gen));
  }
  else {
    fHistGeneratedV0PtTMaxPhi=       (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtTMaxPhi_%i", sysV0Gen));
  }
  if (!fHistGeneratedV0PtTMaxPhi) {cout << "histo gen not found " << endl; return;}

  if (CommonParton && !isCPEff){
    if (!fHistGeneratedV0PtTMaxPhiInt2)  {cout << "histo gen not found " << endl; return;}
    for (Int_t b=1; b<= fHistGeneratedV0PtTMaxPhi->GetBin(fHistGeneratedV0PtTMaxPhi->GetNbinsX(), fHistGeneratedV0PtTMaxPhi->GetNbinsY(), fHistGeneratedV0PtTMaxPhi->GetNbinsZ()); b++){
      //      cout << " b " << b<< endl;
      fHistGeneratedV0PtTMaxPhi->SetBinContent(b, fHistGeneratedV0PtTMaxPhiInt1->GetBinContent(b)-fHistGeneratedV0PtTMaxPhiInt2->GetBinContent(b));
    }
  }


  TH3D*   fHistGeneratedV0PtTMaxEta;
  TH3D*   fHistGeneratedV0PtTMaxEtaInt1;
  TH3D*   fHistGeneratedV0PtTMaxEtaInt2;
  cout << " hey there 1 " << endl;
  if (CommonParton && isCPEff){
    fHistGeneratedV0PtTMaxEta=      (TH3D*)list2->FindObject(Form("fHistCPGeneratedV0PtTMaxEta_%i", sysV0Gen));
  }
  else if (CommonParton && !isCPEff){
    fHistGeneratedV0PtTMaxEtaInt1=       (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtTMaxEta_%i", sysV0Gen));
    fHistGeneratedV0PtTMaxEtaInt2=       (TH3D*)list2->FindObject(Form("fHistCPGeneratedV0PtTMaxEta_%i", sysV0Gen));
    if (fHistGeneratedV0PtTMaxEtaInt1)     fHistGeneratedV0PtTMaxEta= (TH3D*)     fHistGeneratedV0PtTMaxEtaInt1->Clone(Form("fHistNOCPGeneratedV0PtTMaxEta_%i", sysV0Gen));  
  }
  else {
    fHistGeneratedV0PtTMaxEta=       (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtTMaxEta_%i", sysV0Gen));
  }
  if (!fHistGeneratedV0PtTMaxEta) {cout << "histo gen not found " << endl; return;}

  if (CommonParton && !isCPEff){
    if (!fHistGeneratedV0PtTMaxEtaInt2) {cout << "histo gen not found " << endl; return;}
    for (Int_t b=1; b<= fHistGeneratedV0PtTMaxEta->GetBin(fHistGeneratedV0PtTMaxEta->GetNbinsX(), fHistGeneratedV0PtTMaxEta->GetNbinsY(), fHistGeneratedV0PtTMaxEta->GetNbinsZ()); b++){
      //      cout << " b " << b<< endl;
      fHistGeneratedV0PtTMaxEta->SetBinContent(b, fHistGeneratedV0PtTMaxEtaInt1->GetBinContent(b)-fHistGeneratedV0PtTMaxEtaInt2->GetBinContent(b));
    }
  }

  TH3D*   fHistGeneratedV0PtPtTMax;
  TH3D*   fHistGeneratedV0PtPtTMaxInt1;
  TH3D*   fHistGeneratedV0PtPtTMaxInt2;
  cout << " hey there 1 " << endl;
  if (CommonParton && isCPEff){
    fHistGeneratedV0PtPtTMax=       (TH3D*)list2->FindObject(Form("fHistCPGeneratedV0PtPtTMax_%i", sysV0Gen));
  }
  else if (CommonParton && !isCPEff){
    fHistGeneratedV0PtPtTMaxInt1=       (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtPtTMax_%i", sysV0Gen));
    fHistGeneratedV0PtPtTMaxInt2=       (TH3D*)list2->FindObject(Form("fHistCPGeneratedV0PtPtTMax_%i", sysV0Gen));
    if (fHistGeneratedV0PtPtTMaxInt1)    fHistGeneratedV0PtPtTMax= (TH3D*)     fHistGeneratedV0PtPtTMaxInt1->Clone(Form("fHistNOCPGeneratedV0PtPtTMax_%i", sysV0Gen));
  }
  else {
    fHistGeneratedV0PtPtTMax=       (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtPtTMax_%i", sysV0Gen));
  }
  if (!fHistGeneratedV0PtPtTMax) {cout << "histo gen not found " << endl; return;}

  if (CommonParton && !isCPEff){
    if (!fHistGeneratedV0PtPtTMaxInt2) {cout << "histo gen not found " << endl; return;}
    for (Int_t b=1; b<= fHistGeneratedV0PtPtTMax->GetBin(fHistGeneratedV0PtPtTMax->GetNbinsX(), fHistGeneratedV0PtPtTMax->GetNbinsY(), fHistGeneratedV0PtPtTMax->GetNbinsZ()); b++){
      //      cout << " b " << b<< endl;
      fHistGeneratedV0PtPtTMax->SetBinContent(b, fHistGeneratedV0PtPtTMaxInt1->GetBinContent(b)-fHistGeneratedV0PtPtTMaxInt2->GetBinContent(b));
    }
  }

  TH3D*   fHistSelectedV0PtTMaxPhi=       (TH3D*)fileinputSel->Get(Form("fHistSelectedV0PtTMaxPhi_%i", sysV0hh));
  if (CommonParton && isCPEff){
    fHistSelectedV0PtTMaxPhi=   (TH3D*)fileinputSel->Get(Form("fHistCPSelectedV0PtTMaxPhi_%i" , sysV0hh)); 
  }
  else if (CommonParton && !isCPEff){
    fHistSelectedV0PtTMaxPhi=   (TH3D*)fileinputSel->Get(Form("fHistNOCPSelectedV0PtTMaxPhi_%i" , sysV0hh)); 
  }
  else   fHistSelectedV0PtTMaxPhi  =        (TH3D*)fileinputSel->Get(Form("fHistSelectedV0PtTMaxPhi_%i" , sysV0hh)); 
  if (!fHistSelectedV0PtTMaxPhi) {cout << "histo phi assoc sel not found " << endl; return;}

  TH3D*   fHistSelectedV0PtTMaxEta=       (TH3D*)fileinputSel->Get(Form("fHistSelectedV0PtTMaxEta_%i", sysV0hh));
  if (CommonParton && isCPEff){
    fHistSelectedV0PtTMaxEta=   (TH3D*)fileinputSel->Get(Form("fHistCPSelectedV0PtTMaxEta_%i" , sysV0hh)); 
  }
  else if (CommonParton && !isCPEff){
    fHistSelectedV0PtTMaxEta=   (TH3D*)fileinputSel->Get(Form("fHistNOCPSelectedV0PtTMaxEta_%i" , sysV0hh)); 
  }
  else   fHistSelectedV0PtTMaxEta  =        (TH3D*)fileinputSel->Get(Form("fHistSelectedV0PtTMaxEta_%i" , sysV0hh)); 
  if (!fHistSelectedV0PtTMaxEta) {cout << "histo eta assoc sel not found " << endl; return;}

  TH3D*   fHistSelectedV0PtPtTMax;
  if (CommonParton && isCPEff){
    fHistSelectedV0PtPtTMax=   (TH3D*)fileinputSel->Get(Form("fHistCPSelectedV0PtPtTMax_%i" , sysV0hh)); 
  }
  else if (CommonParton && !isCPEff){
    fHistSelectedV0PtPtTMax=   (TH3D*)fileinputSel->Get(Form("fHistNOCPSelectedV0PtPtTMax_%i" , sysV0hh)); 
  }
  else   fHistSelectedV0PtPtTMax  =        (TH3D*)fileinputSel->Get(Form("fHistSelectedV0PtPtTMax_%i" , sysV0hh)); 
  if (!fHistSelectedV0PtPtTMax) {cout << "histo pt assoc sel not found " << endl; return;}

    TH3D*   fHistSelectedGenV0PtPtTMax;

 
    if ( (TH3D*)list2->FindObject("fHistSelectedGenV0PtPtTMax")) { //I don't want it to take this histogram, it's not the correct one since it's taken from the output of the task and not all selections were applied
	fHistSelectedGenV0PtPtTMax=  (TH3D*)list2->FindObject(Form("fHistSelectedGenV0PtPtTMax", sysV0));
    } else {
      cout << " fHistSelectedGenV0PtPtTMax non esiste " << endl;
      fHistSelectedGenV0PtPtTMax=  (TH3D*)fHistSelectedV0PtPtTMax->Clone(Form("fHistSelectedGenV0PtPtTMax", sysV0hh));

    }
    TH3D*   fHistReconstructedV0PtMass= (TH3D*)list->FindObject("fHistReconstructedV0PtMass");
    TH3D*   fHistSelectedV0PtMass=      (TH3D*)list->FindObject("fHistSelectedV0PtMass");
  
    //***************istogrammi delle risoluzioni**********************************************************
    TH2D*   fHistResolution2D[3][2];
    TH3D*   fHistResolution3D[3][2];
    TString nameRes[3][2];
    TH1D*   fHistResolution_1D[nummolt+1][3][2];
    TH1D*   fHistResolution_1D_Int1[nummolt+1][3][2];
    TH1D*   fHistResolution_1D_Int2[nummolt+1][3][2];
    TString nameRes_1D[nummolt+1][3][2]; 
    TString nameRes_2D[nummolt+1][3][2]; 
    TString TorV[2]={"Trigger", "V0"};
    TString Var[3]={"Pt", "Phi", "Eta"};
    cout << "\nI'm taking resolution histos ";
    for(Int_t m=0; m< 3; m++){
      for(Int_t t=0; t< 2; t++){
	nameRes[m][t]="fHistResolution" + TorV[t] + Var[m];
	/* version to be used if resolution hiostograms are 2D, comment if not */
	if (ResoHisto=="2D"){
	  if (listRisoluzione)	fHistResolution2D[m][t]=   (TH2D*)listRisoluzione->FindObject(nameRes[m][t]); //for studies of efficiency with different FB

	  else  fHistResolution2D[m][t]=   (TH2D*)list->FindObject(nameRes[m][t]);       //for thesis + onlyTriggerWithHighestPt
	}
	/* version to be used if resolution hiostograms are 2D, comment if not */

	/* version to be used if resolution hiostograms are 3D, comment if not */
	if (ResoHisto=="3D"){
	  fHistResolution3D[m][t]=   (TH3D*)list->FindObject(nameRes[m][t]); //for latest ones (hhCorr for example)
	}
	/* version to be used if resolution hiostograms are 3D, comment if not */
      }
    }
    cout << "...I've taken resolution histos " << endl;
  auto legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  legend->SetHeader("Multiplicity classes");     

  auto legendPtMin = new TLegend(0.6, 0.6, 0.9, 0.9);
  legendPtMin->SetHeader("Pt Trig Min");     

  auto legenddown = new TLegend(0.6, 0.1, 0.9, 0.4);
  legenddown->SetHeader("Multiplicity classes");     

  //distribuzioni  Pt, Phi, Eta in 2D delle selezionate e generate
  TH2D*  fHistSelected_2D[nummolt+1];
  TH2D*  fHistSelectedGen_2D[nummolt+1];
  TH2D*  fHistGenerated_2D[nummolt+1];
  
  TH2D*  fHistSelected_2D_TriggerPtPhi[nummolt+1];
  TH2D*  fHistSelectedGen_2D_TriggerPtPhi[nummolt+1];
  TH2D*  fHistGenerated_2D_TriggerPtPhi[nummolt+1];
  TH2D*  fHistSelected_2D_TriggerPtEta[nummolt+1];
  TH2D*  fHistGenerated_2D_TriggerPtEta[nummolt+1];

  TH2D*  fHistSelected_2D_V0PtTMaxPhi[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtTMaxPhi[nummolt+1];
  TH2D*  fHistSelected_2D_V0PtTMaxEta[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtTMaxEta[nummolt+1];
  TH2D*   fHistSelected_2D_V0PtPtTMax[nummolt+1];
  TH2D*   fHistSelectedGen_2D_V0PtPtTMax[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtPtTMax[nummolt+1];
  
  TH2D*  fHistSelectedMass_2D[nummolt+1];
  TH2D*  fHistRecoMass_2D[nummolt+1];

  //distribuzioni Pt, Phi, Eta in 1D delle selezionate e generate
  TH1D*  fHistSelected_1D_TriggerPt[nummolt+1];
  TH1D*  fHistSelectedGen_1D_TriggerPt[nummolt+1];
  TH1D*  fHistGenerated_1D_TriggerPt[nummolt+1];
  TH1D*  fHistSelected_1D_V0Pt[nummolt+1];
  TH1D*  fHistSelectedGen_1D_V0Pt[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Pt[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Pt_Int1[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Pt_Int2[nummolt+1];
  TH1D*  fHistSelected_1D_TriggerPhi[nummolt+1];
  TH1D*  fHistGenerated_1D_TriggerPhi[nummolt+1];
  TH1D*  fHistSelected_1D_V0Phi[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Phi[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Phi_Int1[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Phi_Int2[nummolt+1];
  TH1D*  fHistSelected_1D_TriggerEta[nummolt+1];
  TH1D*  fHistGenerated_1D_TriggerEta[nummolt+1];
  TH1D*  fHistSelected_1D_V0Eta[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Eta[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Eta_Int1[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Eta_Int2[nummolt+1];
  TH1D*  fHistSelectedMass_1D[nummolt+1];
  TH1D*  fHistRecoMass_1D[nummolt+1];
  
  TH2D*  fHistTriggerEfficiencyPtPhi[nummolt+1];                                 
  TH2D*  fHistTriggerEfficiencyPtEta[nummolt+1]; 
  TH2D*  fHistV0EfficiencyPtPhi[nummolt+1];
  TH2D*  fHistV0EfficiencyPtEta[nummolt+1];
  TH2D*  fHistV0EfficiencyPtPtTMax[nummolt+1];
  TH2D*  fHistEfficiencyV0Selection[nummolt+1];

  TH1D*  fHistTriggerEfficiencyPt[nummolt+1];
  TH1D*  fHistTriggerEfficiencyGenPt[nummolt+1];
  TH1D*  fHistTriggerEfficiencyPhi[nummolt+1];
  TH1D*  fHistTriggerEfficiencyEta[nummolt+1];
  TH1D*  fHistTriggerEfficiencyPtBins[nummolt+1];
  TH1D*  fHistTriggerEfficiencyGenPtBins[nummolt+1];
  TH1D*  fHistTriggerGeneratedPtBins[nummolt+1];
  TH1D*  fHistTriggerSelectedPtBins[nummolt+1];
  TH1D*  fHistTriggerSelectedGenPtBins[nummolt+1];

  TH1D*  fHistV0EfficiencyPt[nummolt+1];
  TH1D*  fHistV0EfficiencyPhi[nummolt+1];
  TH1D*  fHistV0EfficiencyEta[nummolt+1];
  TH1D*  fHistV0EfficiencyPtBins[nummolt+1];
  TH1D*  fHistV0EfficiencyGenPtBins[nummolt+1];
  TH1D*  fHistV0GeneratedPtBins[nummolt+1];
  TH1D*  fHistV0SelectedPtBins[nummolt+1];
  TH1D*  fHistV0SelectedGenPtBins[nummolt+1];
  
  TH1D*  HistoTriggerEfficiency= new TH1D("HistoTriggerEfficiency", "Trigger selection efficiency vs centrality", nummolt,Nmolt );

  //  TH1D*  fHistV0MeanEfficiencyPt[nummolt+1];

  TH2D* fHistV0EfficiencyReco[nummolt+1];
  TH1D* fHistV0EfficiencyRecoPt[nummolt+1];
  // TH2D* fHistV0EfficiencyRecoMolt;
  // TH1D* fHistV0EfficiencyRecoPtMolt;
  
  //TH1D* fHistV0MeanEfficiencyRecoPt[nummolt+1];

  //istogrammi per contaminazioni
  TH2D * HistContTrigger[nummolt+1];
  TH3D * HistContV0PtTMax[nummolt+1];
  TH2D * HistContV0[nummolt+1];
  TH1D * HistContTriggerPt[nummolt+1];
  TH1D * HistContV0Pt[nummolt+1];
  TH1D * HistContV0PtBins[nummolt+1];
  TH1D * HistInt[nummolt+1];
  TH1D * HistContTriggerMolt = new TH1D("HistContTriggerMolt", "Trigger contamination factor vs multiplicity", nummolt, Nmolt);

  Float_t ContTrigger[nummolt+1];
  Float_t ContV0[nummolt+1];
  Float_t ContTriggerInt;
  Float_t ContV0Int;
  Float_t ContTriggerIntError;
  Float_t ContV0IntError;

  Float_t TriggerEfficiency[nummolt+1];

  Double_t SelEntriespT3[nummolt+1];
  Double_t SelEntries[nummolt+1];
  Double_t GenEntries[nummolt+1];
  
  cout << "*************************************************************************************" << endl;
  cout << "Sto effettuando istogrammi efficienze in range di molteplicita'" << endl;
    
  for(Int_t molt=0; molt < nummolt+1; molt++){
    //    if (molt==0) continue;
    cout << "\n\n I'm analyzing multiplcity interval n. " << molt<< endl;
  
    cout << "Trigger 2D projection in Phi and Pt " << endl;
    if(molt < nummolt){
      fHistSelectedTriggerPtPhi->GetZaxis()->SetRange(fHistSelectedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistSelectedGenTriggerPtPhi->GetZaxis()->SetRange(fHistSelectedGenTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedGenTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistGeneratedTriggerPtPhi->GetZaxis()->SetRange(fHistGeneratedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    else{
      fHistSelectedTriggerPtPhi->GetZaxis()->SetRange(0,100);
      fHistSelectedGenTriggerPtPhi->GetZaxis()->SetRange(0,100); 
      fHistGeneratedTriggerPtPhi->GetZaxis()->SetRange(0,100);
    }
    fHistSelected_2D_TriggerPtPhi[molt] = (TH2D*)fHistSelectedTriggerPtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistSelectedGen_2D_TriggerPtPhi[molt] = (TH2D*)fHistSelectedGenTriggerPtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis //!! 

    fHistGenerated_2D_TriggerPtPhi[molt] = (TH2D*)fHistGeneratedTriggerPtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistSelected_2D_TriggerPtPhi[molt] ->SetName("fHistSelected_2D_TriggerPtPhi_"+ Smolt[molt] );
    fHistSelectedGen_2D_TriggerPtPhi[molt] ->SetName("fHistSelectedGen_2D_TriggerPtPhi_"+ Smolt[molt] );
    fHistGenerated_2D_TriggerPtPhi[molt]->SetName("fHistGenerated_2D_TriggerPtPhi_"+ Smolt[molt]);
    
    cout << "\nsel2D " <<     fHistSelected_2D_TriggerPtPhi[molt]->GetNbinsX() <<   "  " <<   fHistSelected_2D_TriggerPtPhi[molt]->GetNbinsY() <<" " <<    fHistSelected_2D_TriggerPtPhi[molt]->GetName() <<  endl;
    cout << "selgen2D " <<     fHistSelectedGen_2D_TriggerPtPhi[molt]->GetNbinsX() <<   "  " <<   fHistSelectedGen_2D_TriggerPtPhi[molt]->GetNbinsY() << "  " <<fHistSelectedGen_2D_TriggerPtPhi[molt]->GetName() << endl;
    cout << "gen2D " <<     fHistGenerated_2D_TriggerPtPhi[molt]->GetNbinsX() <<   "  " <<   fHistGenerated_2D_TriggerPtPhi[molt]->GetNbinsY() << "\n" << endl;

    fHistTriggerEfficiencyPtPhi[molt]= new TH2D("fHistTriggerEfficiencyPtPhi_"+ Smolt[molt],"fHistTriggerEfficiencyPtPhi_"+ Smolt[molt],fHistSelectedTriggerPtPhi->GetNbinsX(),fHistSelectedTriggerPtPhi->GetXaxis()->GetXmin(), fHistSelectedTriggerPtPhi->GetXaxis()->GetXmax(),fHistSelectedTriggerPtPhi->GetNbinsY(),fHistSelectedTriggerPtPhi->GetYaxis()->GetBinLowEdge(1), fHistSelectedTriggerPtPhi->GetYaxis()->GetBinUpEdge(fHistSelectedTriggerPtPhi->GetNbinsY()) );
    fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->SetTitle("p_{T}^{Trigg} (GeV/c)");      
    fHistTriggerEfficiencyPtPhi[molt]->GetYaxis()->SetTitle("#phi_{Trigg}");   

    fHistSelected_2D_TriggerPtPhi[molt] ->RebinX(2);   
    fHistSelected_2D_TriggerPtPhi[molt] ->RebinY(5);   

    fHistSelectedGen_2D_TriggerPtPhi[molt] ->RebinX(2); 
    fHistSelectedGen_2D_TriggerPtPhi[molt] ->RebinY(5); 

    if (type==6)    fHistGenerated_2D_TriggerPtPhi[molt]->RebinX(2);    
    else     fHistGenerated_2D_TriggerPtPhi[molt]->RebinX(4);    
    fHistGenerated_2D_TriggerPtPhi[molt]->RebinY(5);    
    
    fHistTriggerEfficiencyPtPhi[molt]->RebinX(2);    
    fHistTriggerEfficiencyPtPhi[molt]->RebinY(5);    
   
    
    cout << "after reebinning " << endl;
    cout << "\nsel2D " <<     fHistSelected_2D_TriggerPtPhi[molt]->GetNbinsX() <<   "  " <<   fHistSelected_2D_TriggerPtPhi[molt]->GetNbinsY() << endl;
    cout << "nselgen2D " <<     fHistSelectedGen_2D_TriggerPtPhi[molt]->GetNbinsX() <<   "  " <<   fHistSelectedGen_2D_TriggerPtPhi[molt]->GetNbinsY() << endl;
    cout << "gen2D " <<     fHistGenerated_2D_TriggerPtPhi[molt]->GetNbinsX() <<   "  " <<   fHistGenerated_2D_TriggerPtPhi[molt]->GetNbinsY() << endl;
    cout << "eff2D " <<        fHistTriggerEfficiencyPtPhi[molt]->GetNbinsX() <<   "  " <<      fHistTriggerEfficiencyPtPhi[molt]->GetNbinsY() << endl;
    
    fHistTriggerEfficiencyPtPhi[molt]->Divide(fHistSelected_2D_TriggerPtPhi[molt], fHistGenerated_2D_TriggerPtPhi[molt]); 
  
    cout << "Trigger 1D projection in Phi and Pt " << endl;
    fHistSelected_1D_TriggerPt[molt]=(TH1D*)fHistSelected_2D_TriggerPtPhi[molt]->ProjectionX("fHistSelected_1D_TriggerPt_"+ Smolt[molt]) ;
    fHistSelectedGen_1D_TriggerPt[molt]=(TH1D*)fHistSelectedGen_2D_TriggerPtPhi[molt]->ProjectionX("fHistSelectedGen_1D_TriggerPt_"+ Smolt[molt]) ; 
    // fHistSelected_1D_TriggerPt[molt]->Rebin(1);
    fHistGenerated_1D_TriggerPt[molt]=(TH1D*)fHistGenerated_2D_TriggerPtPhi[molt]->ProjectionX("fHistGenerated_1D_TriggerPt_"+ Smolt[molt]) ;
    //    fHistGenerated_1D_TriggerPt[molt]->Rebin(1);

    //  fHistTriggerEfficiencyPt[molt]= new TH1D("fHistTriggerEfficiencyPt_"+ Smolt[molt] , "fHistTriggerEfficiencyPt_"+ Smolt[molt] ,  fHistTriggerEfficiencyPtPhi[molt]->GetNbinsX(), fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->GetXmin(), fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->GetXmax() );
    fHistTriggerEfficiencyPt[molt]= (TH1D*)fHistSelected_1D_TriggerPt[molt]->Clone("fHistTriggerEfficiencyPt_"+ Smolt[molt]);
    fHistTriggerEfficiencyGenPt[molt]= (TH1D*)fHistSelectedGen_1D_TriggerPt[molt]->Clone("fHistTriggerEfficiencyGenPt_"+ Smolt[molt]); 
   
    if ( !   fHistTriggerEfficiencyGenPt[molt]) cout << "non ho histo efficiency genpt" << endl;
    fHistTriggerEfficiencyPt[molt]->GetXaxis()->SetTitle("p_{T}^{Trigg} (GeV/c)");      
    fHistTriggerEfficiencyPt[molt]->GetXaxis()->SetTitleSize(0.039);
    fHistTriggerEfficiencyPt[molt]->GetXaxis()->SetTitleOffset(1.2);
    fHistTriggerEfficiencyPt[molt]->SetStats(0);
    //    fHistTriggerEfficiencyPt[molt]->SetTitle("fHistTriggerEfficiencyPt_"+ Smolt[molt]);

    fHistTriggerEfficiencyGenPt[molt]->GetXaxis()->SetTitle("p_{T}^{Trigg} (GeV/c)");      
    fHistTriggerEfficiencyGenPt[molt]->GetXaxis()->SetTitleSize(0.039);
    fHistTriggerEfficiencyGenPt[molt]->GetXaxis()->SetTitleOffset(1.2);
    fHistTriggerEfficiencyGenPt[molt]->SetStats(0);

    fHistTriggerEfficiencyPt[molt] ->Divide(   fHistGenerated_1D_TriggerPt[molt] );
    cout << "I'm performing a division of histograms " << endl;    
    cout << "numerator number of bins " << fHistTriggerEfficiencyPt[molt]->GetNbinsX()<< endl;
    cout << "denominator number of bins " << fHistGenerated_1D_TriggerPt[molt]->GetNbinsX()<< endl;
    fHistTriggerEfficiencyGenPt[molt] ->Divide(  fHistGenerated_1D_TriggerPt[molt] );
    cout << "numerator number of bins " << fHistTriggerEfficiencyGenPt[molt]->GetNbinsX()<< endl;
    cout << "denominator number of bins " << fHistGenerated_1D_TriggerPt[molt]->GetNbinsX()<< endl;
    
    canvasEff->cd(1);
    fHistTriggerEfficiencyPt[molt]->GetYaxis()->SetRangeUser(0,1.6);
    fHistTriggerEfficiencyPt[molt]->GetYaxis()->SetTitle("#epsilon_{Trigg}");
    fHistTriggerEfficiencyPt[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyPt[molt]->GetYaxis()->SetTitleOffset(1);
    fHistTriggerEfficiencyPt[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyPt[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyPt[molt]->SetMarkerColor(Color[molt]);
    legend->AddEntry(fHistTriggerEfficiencyPt[molt],SmoltLegend[molt],"pel");   
    legenddown->AddEntry(fHistTriggerEfficiencyPt[molt],SmoltLegend[molt],"pel");   
    if (ResoHisto=="2D")    legendPtMin->AddEntry(fHistTriggerEfficiencyPt[molt],Form("p_{T}^{Trig, min} > %.0f", ptjminBisFixed[molt]),"pel");   

    fHistTriggerEfficiencyPt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();

    /*
      canvasEffBis[0]->cd();
      fHistTriggerEfficiencyPt[molt]->GetYaxis()->SetRangeUser(0,1.2);
      fHistTriggerEfficiencyPt[molt]->Draw("same");
      if(molt ==nummolt)     legend->Draw();
    */


    fHistTriggerEfficiencyPhi[molt]= new TH1D("fHistTriggerEfficiencyPhi_"+ Smolt[molt] , "fHistTriggerEfficiencyPhi_"+ Smolt[molt] ,  fHistTriggerEfficiencyPtPhi[molt]->GetNbinsY(), fHistTriggerEfficiencyPtPhi[molt]->GetYaxis()->GetBinLowEdge(1), fHistTriggerEfficiencyPtPhi[molt]->GetYaxis()->GetBinUpEdge(fHistTriggerEfficiencyPtPhi[molt]->GetNbinsY()) );
    fHistTriggerEfficiencyPhi[molt]->GetXaxis()->SetTitle("#phi_{Trigg}");      
    fHistTriggerEfficiencyPhi[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyPhi[molt]->GetXaxis()->SetTitleOffset(0.8);
    //    fHistTriggerEfficiencyPhi[molt]->SetTitle("fHistTriggerEfficiencyPhi_"+ Smolt[molt] );  
    fHistSelected_1D_TriggerPhi[molt]=(TH1D*)fHistSelected_2D_TriggerPtPhi[molt]->ProjectionY("fHistSelected_1D_TriggerPhi_"+ Smolt[molt],fHistSelected_2D_TriggerPtPhi[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistSelected_2D_TriggerPtPhi[molt]->GetXaxis()->FindBin(ptjmax -0.0001) );
    fHistGenerated_1D_TriggerPhi[molt]=(TH1D*)fHistGenerated_2D_TriggerPtPhi[molt]->ProjectionY("fHistGenerated_1D_TriggerPhi_"+ Smolt[molt],fHistSelected_2D_TriggerPtPhi[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistSelected_2D_TriggerPtPhi[molt]->GetXaxis()->FindBin(ptjmax -0.0001));
    fHistTriggerEfficiencyPhi[molt] ->Divide (  fHistSelected_1D_TriggerPhi[molt], fHistGenerated_1D_TriggerPhi[molt]);
    /*
      for(Int_t j=0; j<fHistTriggerEfficiencyPhi[molt]->GetNbinsX(); j++ ){

      fHistTriggerEfficiencyPhi[molt]->SetBinError(j+1, SetEfficiencyError(    fHistSelected_1D_TriggerPhi[molt]->GetBinContent(j+1),     fHistGenerated_1D_TriggerPhi[molt]->GetBinContent(j+1)));
      }
    */

    canvasEff->cd(2);
    fHistTriggerEfficiencyPhi[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyPhi[molt]->GetYaxis()->SetTitle("#epsilon_{Trigg}");
    fHistTriggerEfficiencyPhi[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyPhi[molt]->GetYaxis()->SetTitleOffset(1);
    fHistTriggerEfficiencyPhi[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyPhi[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyPhi[molt]->SetMarkerColor(Color[molt]);
    fHistTriggerEfficiencyPhi[molt]->SetStats(0);
    fHistTriggerEfficiencyPhi[molt]->Draw("same");
    if(molt ==nummolt)     legenddown->Draw();
  
    /*
      canvasEffBis[1]->cd();
      fHistTriggerEfficiencyPhi[molt]->Draw("same");
      if(molt ==nummolt)     legenddown->Draw();
    */


    //*************Trigger efficiency in pT bins of different width************************
    cout << "Trigger efficiency in Pt bins used in the analysis " << endl;
    fHistTriggerSelectedPtBins[molt]= new TH1D("fHistTriggerSelectedPtBins_" + Smolt[molt], "fHistTriggerSelectedPtBins_" + Smolt[molt],numPtTrigger, NPtTrigger );
    fHistTriggerSelectedGenPtBins[molt]= new TH1D("fHistTriggerSelectedGenPtBins_" + Smolt[molt], "fHistTriggerSelectedGenPtBins_" + Smolt[molt],numPtTrigger, NPtTrigger );
    fHistTriggerGeneratedPtBins[molt]= new TH1D("fHistTriggerGeneratedPtBins_" + Smolt[molt], "fHistTriggerGeneratedPtBins_" + Smolt[molt],numPtTrigger, NPtTrigger );
    fHistTriggerEfficiencyPtBins[molt]= new TH1D("fHistTriggerEfficiencyPtBins_" + Smolt[molt], "fHistTriggerEfficiencyPtBins_" + Smolt[molt],numPtTrigger, NPtTrigger );
    fHistTriggerEfficiencyGenPtBins[molt]= new TH1D("fHistTriggerEfficiencyGenPtBins_" + Smolt[molt], "fHistTriggerEfficiencyGenPtBins_" + Smolt[molt],numPtTrigger, NPtTrigger );
    fHistTriggerEfficiencyPtBins[molt]->GetXaxis()->SetTitle("p_{T}^{Trigger} (GeV/c)");      
    fHistTriggerEfficiencyPtBins[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyPtBins[molt]->GetXaxis()->SetTitleOffset(1);
    Float_t NumberOfSelectedTrigger;
    Float_t NumberOfSelectedGenTrigger;
    Float_t NumberOfGeneratedTrigger;
    cout << "\n Trigger efficiency in Pt bins " << endl;
    for(Int_t j=0; j<numPtTrigger; j++){
      NumberOfSelectedTrigger=0;
      NumberOfGeneratedTrigger=0;
      NumberOfSelectedGenTrigger=0;
      for(Int_t i=fHistSelected_1D_TriggerPt[molt]->GetXaxis()->FindBin(NPtTrigger[j]+0.00001); i<=fHistSelected_1D_TriggerPt[molt]->GetXaxis()->FindBin(NPtTrigger[j+1]-0.00001); i++){
	NumberOfGeneratedTrigger+=fHistGenerated_1D_TriggerPt[molt]->GetBinContent(i);
	NumberOfSelectedTrigger+=fHistSelected_1D_TriggerPt[molt]->GetBinContent(i);
	NumberOfSelectedGenTrigger+=fHistSelectedGen_1D_TriggerPt[molt]->GetBinContent(i);
      }
      fHistTriggerGeneratedPtBins[molt]->SetBinContent(j+1, NumberOfGeneratedTrigger);
      fHistTriggerSelectedPtBins[molt]->SetBinContent(j+1, NumberOfSelectedTrigger);
      fHistTriggerSelectedGenPtBins[molt]->SetBinContent(j+1, NumberOfSelectedGenTrigger);
      fHistTriggerEfficiencyPtBins[molt]->SetBinContent(j+1, NumberOfSelectedTrigger/NumberOfGeneratedTrigger);
      fHistTriggerEfficiencyPtBins[molt]->SetBinError(j+1, SetEfficiencyError(NumberOfSelectedTrigger,NumberOfGeneratedTrigger));

      fHistTriggerEfficiencyGenPtBins[molt]->SetBinContent(j+1, NumberOfSelectedGenTrigger/NumberOfGeneratedTrigger);
      fHistTriggerEfficiencyGenPtBins[molt]->SetBinError(j+1, SetEfficiencyError(NumberOfSelectedGenTrigger,NumberOfGeneratedTrigger));

      cout << 	"bin at pT,Trigger: " <<NPtTrigger[j] << "  " <<fHistTriggerEfficiencyPtBins[molt]->GetBinContent(j+1) << "+-" << 	fHistTriggerEfficiencyPtBins[molt]->GetBinError(j+1)<< endl;

    }

    fHistTriggerEfficiencyPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyGenPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    if (ishhCorr)     fHistTriggerEfficiencyPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    if (ishhCorr)     fHistTriggerEfficiencyGenPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyPtBins[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyPtBins[molt]->GetYaxis()->SetTitle("#epsilon_{Trigger}");
    fHistTriggerEfficiencyPtBins[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyPtBins[molt]->GetYaxis()->SetTitleOffset(1);
    fHistTriggerEfficiencyPtBins[molt]->SetStats(0);
    fHistTriggerEfficiencyPtBins[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyPtBins[molt]->SetMarkerColor(Color[molt]);
    fHistTriggerEfficiencyGenPtBins[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyGenPtBins[molt]->SetMarkerColor(Color[molt]);

    //*************End of rigger efficiency in pT bins of different width************************


    SelEntriespT3[molt]=0;
    SelEntries[molt]=0;
    GenEntries[molt]=0;

    //Trigger efficiency integrated in all variables except multiplicity (ptmin < pt < ptmax)
    for(Int_t j=fHistSelected_1D_TriggerPt[molt]->FindBin(3); j<fHistSelected_1D_TriggerPt[molt]->FindBin(ptjmax); j++ ){
      SelEntriespT3[molt]+=fHistSelected_1D_TriggerPt[molt]->GetBinContent(j);
    }
    for(Int_t j=fHistSelected_1D_TriggerPt[molt]->FindBin(ptjmin); j<fHistSelected_1D_TriggerPt[molt]->FindBin(ptjmax); j++ ){
      SelEntries[molt]+=fHistSelected_1D_TriggerPt[molt]->GetBinContent(j);
    }
    for(Int_t j=fHistGenerated_1D_TriggerPt[molt]->FindBin(ptjmin); j<fHistGenerated_1D_TriggerPt[molt]->FindBin(ptjmax); j++ ){
      GenEntries[molt]+=fHistGenerated_1D_TriggerPt[molt]->GetBinContent(j);
    }
    TriggerEfficiency[molt]=(SelEntries[molt]/GenEntries[molt]);
    HistoTriggerEfficiency->SetBinContent(molt+1, TriggerEfficiency[molt]);
    HistoTriggerEfficiency->SetBinError(molt+1, SetEfficiencyError(SelEntries[molt], GenEntries[molt]));

    HistoTriggerEfficiency->GetYaxis()->SetRangeUser(0,1);
    HistoTriggerEfficiency->SetMarkerStyle(ColorSysTrigger[sysTrigger]);
    HistoTriggerEfficiency->SetLineColor(ColorSysTrigger[sysTrigger]);
    HistoTriggerEfficiency->SetMarkerColor(ColorSysTrigger[sysTrigger]);
    canvasUsed->cd(1);
    HistoTriggerEfficiency->Draw("e");
      
    //per la particella di trigger - PtEta
    if(molt < nummolt){
      fHistSelectedTriggerPtEta->GetZaxis()->SetRange(fHistSelectedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistGeneratedTriggerPtEta->GetZaxis()->SetRange(fHistGeneratedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    else{
      fHistSelectedTriggerPtEta->GetZaxis()->SetRange(0,100);
      fHistGeneratedTriggerPtEta->GetZaxis()->SetRange(0,100);
    }
    fHistSelected_2D_TriggerPtEta[molt] = (TH2D*)fHistSelectedTriggerPtEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_TriggerPtEta[molt] = (TH2D*)fHistGeneratedTriggerPtEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axi
    fHistSelected_2D_TriggerPtEta[molt] ->SetName("fHistSelected_2D_TriggerPtEta_"+ Smolt[molt] );
    fHistGenerated_2D_TriggerPtEta[molt]->SetName("fHistGenerated_2D_TriggerPtEta_"+ Smolt[molt]);
    fHistTriggerEfficiencyPtEta[molt]= new TH2D("fHistTriggerEfficiencyPtEta_"+ Smolt[molt],"fHistTriggerEfficiencyPtEta_"+ Smolt[molt],fHistSelectedTriggerPtEta->GetNbinsX(),fHistSelectedTriggerPtEta->GetXaxis()->GetXmin(), fHistSelectedTriggerPtEta->GetXaxis()->GetXmax(),fHistSelectedTriggerPtEta->GetNbinsY(),fHistSelectedTriggerPtEta->GetYaxis()->GetBinLowEdge(1), fHistSelectedTriggerPtEta->GetYaxis()->GetBinUpEdge(fHistSelectedTriggerPtEta->GetNbinsY()) );
    fHistTriggerEfficiencyPtEta[molt]->GetXaxis()->SetTitle("p_{T}^{Trigg} (GeV/c)");      
    fHistTriggerEfficiencyPtEta[molt]->GetYaxis()->SetTitle("#phi_{Trigg}");      
    
    fHistSelected_2D_TriggerPtEta[molt] ->RebinX(2);   
    fHistSelected_2D_TriggerPtEta[molt] ->RebinY(10);   
    
    if (type==6)    fHistGenerated_2D_TriggerPtEta[molt]->RebinX(2);    
    else fHistGenerated_2D_TriggerPtEta[molt]->RebinX(4);    
    fHistGenerated_2D_TriggerPtEta[molt]->RebinY(10);    
    
    fHistTriggerEfficiencyPtEta[molt]->RebinX(2);    
    fHistTriggerEfficiencyPtEta[molt]->RebinY(10);    
    
    fHistTriggerEfficiencyPtEta[molt]->Divide(fHistSelected_2D_TriggerPtEta[molt], fHistGenerated_2D_TriggerPtEta[molt]); 
    
    fHistTriggerEfficiencyEta[molt]= new TH1D("fHistTriggerEfficiencyEta_"+ Smolt[molt] , "fHistTriggerEfficiencyEta_"+ Smolt[molt] ,  fHistTriggerEfficiencyPtEta[molt]->GetNbinsY(), fHistTriggerEfficiencyPtEta[molt]->GetYaxis()->GetBinLowEdge(1), fHistTriggerEfficiencyPtEta[molt]->GetYaxis()->GetBinUpEdge(fHistTriggerEfficiencyPtEta[molt]->GetNbinsY()) );
    fHistTriggerEfficiencyEta[molt]->GetXaxis()->SetTitle("#eta_{Trigg}");
    fHistTriggerEfficiencyEta[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyEta[molt]->GetXaxis()->SetTitleOffset(0.9);      
    fHistTriggerEfficiencyEta[molt]->SetStats(0);
    fHistTriggerEfficiencyEta[molt]->SetTitle("fHistTriggerEfficiencyEta_"+ Smolt[molt] );  
    fHistSelected_1D_TriggerEta[molt]=(TH1D*)fHistSelected_2D_TriggerPtEta[molt]->ProjectionY("fHistSelected_1D_TriggerEta_"+ Smolt[molt],fHistSelected_2D_TriggerPtEta[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistSelected_2D_TriggerPtEta[molt]->GetXaxis()->FindBin(ptjmax -0.0001));
    fHistGenerated_1D_TriggerEta[molt]=(TH1D*)fHistGenerated_2D_TriggerPtEta[molt]->ProjectionY("fHistGenerated_1D_TriggerEta_"+ Smolt[molt], fHistGenerated_2D_TriggerPtEta[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistGenerated_2D_TriggerPtEta[molt]->GetXaxis()->FindBin(ptjmax -0.0001));
    
    fHistTriggerEfficiencyEta[molt] ->Divide (fHistSelected_1D_TriggerEta[molt], fHistGenerated_1D_TriggerEta[molt]);
    /*
      for(Int_t j=0; j<fHistTriggerEfficiencyEta[molt]->GetNbinsX(); j++ ){

      fHistTriggerEfficiencyEta[molt]->SetBinError(j+1, SetEfficiencyError(    fHistSelected_1D_TriggerEta[molt]->GetBinContent(j+1),     fHistGenerated_1D_TriggerEta[molt]->GetBinContent(j+1)));
      }
    */
    canvasEff->cd(3);
    fHistTriggerEfficiencyEta[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyEta[molt]->GetYaxis()->SetTitle("#epsilon_{Trigg}");
    fHistTriggerEfficiencyEta[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistTriggerEfficiencyEta[molt]->GetYaxis()->SetTitleOffset(1);
    fHistTriggerEfficiencyEta[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyEta[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyEta[molt]->SetMarkerColor(Color[molt]);
    fHistTriggerEfficiencyEta[molt]->Draw("same");
    if(molt ==nummolt)     legenddown->Draw();
  
    /*
      canvasEffBis[2]->cd();
      fHistTriggerEfficiencyEta[molt]->GetXaxis()->SetRangeUser(-1,1);
      fHistTriggerEfficiencyEta[molt]->Draw("same");
      if(molt ==nummolt)     legenddown->Draw();
    */

    /* per calcolo errore efficienza
       for(Int_t i=1 ; i< fHistTriggerEfficiencyEta[molt]->GetNbinsX(); i++){
       //cfHistTriggerEfficiencyEta[molt]->SetBinContent(i,(Float_t)((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i)/((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i));
       if( ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i) != 0){
       //ccout <<  ((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i) << " " <<  ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i) << " " <<SetEfficiencyError(((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i), ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i))<< endl;
       fHistTriggerEfficiencyEta[molt]->SetBinError(i, SetEfficiencyError(((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i), ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i)));
       }
       }
    */

    cout << "V0 2D projection in Phi and PtTMax (Pt stands for PtTMax) " << endl;
    if(molt < nummolt){
      fHistSelectedV0PtTMaxPhi->GetZaxis()->SetRange(fHistSelectedV0PtTMaxPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtTMaxPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistGeneratedV0PtTMaxPhi->GetZaxis()->SetRange(fHistGeneratedV0PtTMaxPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedV0PtTMaxPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    else{
      fHistSelectedV0PtTMaxPhi->GetZaxis()->SetRange(0,100);
      fHistGeneratedV0PtTMaxPhi->GetZaxis()->SetRange(0,100);
    }
    fHistSelected_2D_V0PtTMaxPhi[molt] = (TH2D*)fHistSelectedV0PtTMaxPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_V0PtTMaxPhi[molt] = (TH2D*)fHistGeneratedV0PtTMaxPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistSelected_2D_V0PtTMaxPhi[molt] ->SetName("fHistSelected_2D_V0PtTMaxPhi_"+ Smolt[molt]);
    fHistGenerated_2D_V0PtTMaxPhi[molt]->SetName("fHistGenerated_2D_V0PtTMaxPhi_"+ Smolt[molt]);
    fHistV0EfficiencyPtPhi[molt]= new TH2D("fHistV0EfficiencyPtPhi_"+ Smolt[molt],"fHistV0EfficiencyPtPhi_"+ Smolt[molt],fHistSelectedV0PtTMaxPhi->GetNbinsX(),fHistSelectedV0PtTMaxPhi->GetXaxis()->GetXmin(), fHistSelectedV0PtTMaxPhi->GetXaxis()->GetXmax(),fHistSelectedV0PtTMaxPhi->GetNbinsY(),fHistSelectedV0PtTMaxPhi->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtTMaxPhi->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtTMaxPhi->GetNbinsY()) );
    fHistV0EfficiencyPtPhi[molt]->GetXaxis()->SetTitle("p_{T}^{Trigger, max} (GeV/c)");      
    fHistV0EfficiencyPtPhi[molt]->GetYaxis()->SetTitle("#phi_{Assoc}");
    fHistSelected_2D_V0PtTMaxPhi[molt] ->RebinX(1);   
    fHistSelected_2D_V0PtTMaxPhi[molt] ->RebinY(5);   

    fHistGenerated_2D_V0PtTMaxPhi[molt]->RebinX(1);    
    fHistGenerated_2D_V0PtTMaxPhi[molt]->RebinY(5);    

    fHistV0EfficiencyPtPhi[molt]->RebinX(1);    
    fHistV0EfficiencyPtPhi[molt]->RebinY(5);    

    cout << "sel2D " <<     fHistSelected_2D_V0PtTMaxPhi[molt]->GetNbinsX() <<   "  " <<   fHistSelected_2D_V0PtTMaxPhi[molt]->GetNbinsY() << endl;
    cout << "gen2D " <<     fHistGenerated_2D_V0PtTMaxPhi[molt]->GetNbinsX() <<   "  " <<   fHistGenerated_2D_V0PtTMaxPhi[molt]->GetNbinsY() << endl;
    cout << "eff2D " <<     fHistV0EfficiencyPtPhi[molt]->GetNbinsX() <<   "  " <<      fHistV0EfficiencyPtPhi[molt]->GetNbinsY() << endl;
    
    fHistV0EfficiencyPtPhi[molt]->Divide(fHistSelected_2D_V0PtTMaxPhi[molt], fHistGenerated_2D_V0PtTMaxPhi[molt]); 

    cout << "V0 1D projection in Phi " << endl;
    fHistV0EfficiencyPhi[molt]= new TH1D("fHistV0EfficiencyPhi_"+ Smolt[molt] , "fHistV0EfficiencyPhi_"+ Smolt[molt] , fHistSelected_2D_V0PtTMaxPhi[molt]->GetNbinsY(), fHistSelected_2D_V0PtTMaxPhi[molt]->GetYaxis()->GetBinLowEdge(1), fHistSelected_2D_V0PtTMaxPhi[molt]->GetYaxis()->GetBinUpEdge(fHistSelected_2D_V0PtTMaxPhi[molt]->GetNbinsY()) );
    fHistV0EfficiencyPhi[molt]->GetXaxis()->SetTitle("#phi_{Assoc}");      
    fHistV0EfficiencyPhi[molt]->GetXaxis()->SetTitleOffset(0.8);
    fHistV0EfficiencyPhi[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyPhi[molt]->SetTitle("fHistV0EfficiencyPhi_"+ Smolt[molt] );

    fHistSelected_1D_V0Phi[molt]=(TH1D*)fHistSelected_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistSelected_1D_V0Phi_"+ Smolt[molt],fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001) ); //***


    //!!!!!!!!!this way of adding the two projections return a wrong histogram (too many entries, it's not the sum of the two projections)!!!!!!!!!!!! I should have changes the name of the second projection to make it work!!
/*
  fHistGenerated_1D_V0Phi[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_"+ Smolt[molt],fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjmin+0.0001) ,fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001));
  cout << "generated histogram entries " << endl;
  cout <<     fHistGenerated_1D_V0Phi[molt]->GetEntries()  << endl;
  cout << ((TH1D*)fHistGenerated_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_"+ Smolt[molt],fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(-ptjmin-0.0001)))->GetEntries() << endl;
  fHistGenerated_1D_V0Phi[molt]->Add((TH1D*)fHistGenerated_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_"+ Smolt[molt],fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(-ptjmin-0.0001)));
  cout <<     fHistGenerated_1D_V0Phi[molt]->GetEntries() << endl;
*/
    if (type>=4 && type!=6){
    fHistGenerated_1D_V0Phi_Int1[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_Int1"+ Smolt[molt],fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjmin+0.0001) ,fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001));
    fHistGenerated_1D_V0Phi_Int2[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_Int2"+ Smolt[molt],fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistSelected_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(-ptjmin-0.0001));
    fHistGenerated_1D_V0Phi[molt]=(TH1D*)     fHistGenerated_1D_V0Phi_Int2[molt]->Clone("fHistGenerated_1D_V0Phi_"+ Smolt[molt]);
    fHistGenerated_1D_V0Phi[molt]->Add(fHistGenerated_1D_V0Phi_Int1[molt]);
    }
    else {
    fHistGenerated_1D_V0Phi[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_"+ Smolt[molt],fHistGenerated_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistGenerated_2D_V0PtTMaxPhi[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001));
    }
    fHistV0EfficiencyPhi[molt] ->Divide(  fHistSelected_1D_V0Phi[molt],  fHistGenerated_1D_V0Phi[molt]);
  
    canvasEff->cd(5);
    fHistV0EfficiencyPhi[molt]->GetYaxis()->SetRangeUser(0.00001,0.2);
    if (type==6)     fHistV0EfficiencyPhi[molt]->GetYaxis()->SetRangeUser(0.00001,1);
    if (ishhCorr)   fHistV0EfficiencyPhi[molt]->GetYaxis()->SetRangeUser(0.00001,1);
    fHistV0EfficiencyPhi[molt]->SetMarkerStyle(Marker[molt]);
    fHistV0EfficiencyPhi[molt]->GetYaxis()->SetTitle("#epsilon_{Assoc}");
    fHistV0EfficiencyPhi[molt]->GetYaxis()->SetTitleOffset(1);
    fHistV0EfficiencyPhi[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyPhi[molt]->SetStats(0);
    fHistV0EfficiencyPhi[molt]->SetLineColor(Color[molt]);
    fHistV0EfficiencyPhi[molt]->SetMarkerColor(Color[molt]);
    fHistV0EfficiencyPhi[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();
  
    /*
      canvasEffBis[4]->cd();
      fHistV0EfficiencyPhi[molt]->Draw("same");
      if(molt ==nummolt)     legend->Draw();
    */

    cout << "V0 2D projection in Pt and PtTMax " << endl;
    if(molt < nummolt){
      fHistSelectedV0PtPtTMax->GetZaxis()->SetRange(fHistSelectedV0PtPtTMax->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtPtTMax->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistSelectedGenV0PtPtTMax->GetZaxis()->SetRange(fHistSelectedGenV0PtPtTMax->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedGenV0PtPtTMax->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistGeneratedV0PtPtTMax->GetZaxis()->SetRange(fHistGeneratedV0PtPtTMax->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedV0PtPtTMax->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    else{
      fHistSelectedV0PtPtTMax->GetZaxis()->SetRange(0,100);
      fHistSelectedGenV0PtPtTMax->GetZaxis()->SetRange(0,100);
      fHistGeneratedV0PtPtTMax->GetZaxis()->SetRange(0,100);
    }
    fHistSelected_2D_V0PtPtTMax[molt] = (TH2D*)fHistSelectedV0PtPtTMax->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistSelectedGen_2D_V0PtPtTMax[molt] = (TH2D*)fHistSelectedGenV0PtPtTMax->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_V0PtPtTMax[molt] = (TH2D*)fHistGeneratedV0PtPtTMax->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistSelected_2D_V0PtPtTMax[molt] ->SetName("fHistSelected_2D_V0PtPtTMax_"+ Smolt[molt]);
    fHistSelectedGen_2D_V0PtPtTMax[molt] ->SetName("fHistSelectedGen_2D_V0PtPtTMax_"+ Smolt[molt]);
    fHistGenerated_2D_V0PtPtTMax[molt]->SetName("fHistGenerated_2D_V0PtPtTMax_"+ Smolt[molt]);

    fHistV0EfficiencyPtPtTMax[molt]= new TH2D("fHistV0EfficiencyPtPtTMax_"+ Smolt[molt],"fHistV0EfficiencyPtPtTMax_"+ Smolt[molt],fHistSelectedV0PtPtTMax->GetNbinsX(),fHistSelectedV0PtPtTMax->GetXaxis()->GetXmin(), fHistSelectedV0PtPtTMax->GetXaxis()->GetXmax(),fHistSelectedV0PtPtTMax->GetNbinsY(),fHistSelectedV0PtPtTMax->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtPtTMax->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtPtTMax->GetNbinsY()) );
    fHistV0EfficiencyPtPtTMax[molt]->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c)");      
    fHistV0EfficiencyPtPtTMax[molt]->GetYaxis()->SetTitle("#phi_{Assoc}");

    fHistSelected_2D_V0PtPtTMax[molt] ->RebinX(1);   
    fHistSelected_2D_V0PtPtTMax[molt] ->RebinX(1);   
    fHistSelectedGen_2D_V0PtPtTMax[molt] ->RebinY(1);   

    fHistGenerated_2D_V0PtPtTMax[molt]->RebinX(1);    
    fHistGenerated_2D_V0PtPtTMax[molt]->RebinY(1);    

    fHistV0EfficiencyPtPtTMax[molt]->RebinX(1);    
    fHistV0EfficiencyPtPtTMax[molt]->RebinY(1);    

    fHistV0EfficiencyPtPtTMax[molt]->Divide(fHistSelected_2D_V0PtPtTMax[molt], fHistGenerated_2D_V0PtPtTMax[molt]); 
 
    cout << "V0 1D projection in Pt with correct cut on pT(Trigger, soglia)" << endl;
    
    fHistSelected_1D_V0Pt[molt]=(TH1D*)fHistSelected_2D_V0PtPtTMax[molt]->ProjectionX("fHistSelected_1D_V0Pt_"+ Smolt[molt],fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjminGen+0.0001) ,fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjmaxGen-0.0001));
    fHistSelectedGen_1D_V0Pt[molt]=(TH1D*)fHistSelectedGen_2D_V0PtPtTMax[molt]->ProjectionX("fHistSelectedGen_1D_V0Pt_"+ Smolt[molt],fHistSelectedGen_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjminGen+0.0001) ,fHistSelectedGen_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjmaxGen-0.0001));

    if (type>=4 && type!=6){
    fHistGenerated_1D_V0Pt_Int1[molt]=(TH1D*)fHistGenerated_2D_V0PtPtTMax[molt]->ProjectionX("fHistGenerated_1D_V0Pt_Int1"+ Smolt[molt],fHistGenerated_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjmin+0.0001) ,fHistGenerated_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjmaxGen-0.0001));
    fHistGenerated_1D_V0Pt_Int2[molt]=(TH1D*)fHistGenerated_2D_V0PtPtTMax[molt]->ProjectionX("fHistGenerated_1D_V0Pt_Int2"+ Smolt[molt],fHistGenerated_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjminGen+0.0001) ,fHistGenerated_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(-ptjmin-0.0001));
    fHistGenerated_1D_V0Pt[molt]=(TH1D*)    fHistGenerated_1D_V0Pt_Int2[molt]->Clone("fHistGenerated_1D_V0Pt_"+ Smolt[molt]);
    fHistGenerated_1D_V0Pt[molt]->Add( fHistGenerated_1D_V0Pt_Int1[molt]);
    }
    else {
    fHistGenerated_1D_V0Pt[molt]=(TH1D*)fHistGenerated_2D_V0PtPtTMax[molt]->ProjectionX("fHistGenerated_1D_V0Pt_"+ Smolt[molt],fHistGenerated_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjminGen+0.0001) ,fHistGenerated_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(ptjmaxGen-0.0001));
    }

    fHistV0EfficiencyPt[molt]= (TH1D*) fHistSelected_1D_V0Pt[molt]->Clone("fHistV0EfficiencyPt_"+ Smolt[molt]);
    fHistV0EfficiencyPt[molt]->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c)");      
    fHistV0EfficiencyPt[molt]->GetXaxis()->SetTitleSize(0.039);      
    fHistV0EfficiencyPt[molt]->GetXaxis()->SetTitleOffset(1);
    fHistV0EfficiencyPt[molt]->GetYaxis()->SetTitleSize(0.05);      
    fHistV0EfficiencyPt[molt]->GetYaxis()->SetTitleOffset(1);
    fHistV0EfficiencyPt[molt]->SetTitle("fHistV0EfficiencyPt_"+ Smolt[molt] );

cout << "sel2D " <<     fHistSelected_1D_V0Pt[molt]->GetNbinsX() << endl;
    cout << "gen2D " <<     fHistGenerated_1D_V0Pt[molt]->GetNbinsX() <<endl;
    cout << "eff2D " <<     fHistV0EfficiencyPt[molt]->GetNbinsX() <<   endl;

    fHistV0EfficiencyPt[molt]->Divide(fHistSelected_1D_V0Pt[molt],  fHistGenerated_1D_V0Pt[molt]);

    /*
    cout << "\n\n\n*************************" << endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(1+0.0001)<< endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(2+0.0001)<< endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(3+0.0001)<< endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(4+0.0001)<< endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(5+0.0001)<< endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(6+0.0001)<< endl;
    cout <<   fHistSelected_2D_V0PtPtTMax[molt]->GetYaxis()->FindBin(7+0.0001)<< endl;
    cout << "*************************\n\n\n" << endl;
    */

    canvasEff->cd(4);
    fHistV0EfficiencyPt[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistV0EfficiencyPt[molt]->SetMarkerStyle(Marker[molt]);
    fHistV0EfficiencyPt[molt]->SetLineColor(Color[molt]);
    fHistV0EfficiencyPt[molt]->SetMarkerColor(Color[molt]);
    fHistV0EfficiencyPt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();

    cout << "V0 efficiency in Pt bins used in the analysis " << endl;
    fHistV0EfficiencyPtBins[molt]= new TH1D("fHistV0EfficiencyPtBins_" + Smolt[molt], "fHistV0EfficiencyPtBins_" + Smolt[molt],numPtV0, NPtV0 );
    fHistV0EfficiencyPtBins[molt]->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c)");      
    fHistV0EfficiencyPtBins[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyPtBins[molt]->GetXaxis()->SetTitleOffset(1);
    fHistV0EfficiencyGenPtBins[molt]= new TH1D("fHistV0EfficiencyGenPtBins_" + Smolt[molt], "fHistV0EfficiencyGenPtBins_" + Smolt[molt],numPtV0, NPtV0 );
    fHistV0EfficiencyGenPtBins[molt]->GetXaxis()->SetTitle("p_{T}^{Assoc} (GeV/c)");      
    fHistV0EfficiencyGenPtBins[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyGenPtBins[molt]->GetXaxis()->SetTitleOffset(1);

    Float_t NumberOfSelected;
    Float_t NumberOfSelectedGen;
    Float_t NumberOfGenerated;
    cout << "\n V0 efficiency in Pt bins " << endl;
    for(Int_t j=0; j<numPtV0Max; j++){
      NumberOfSelected=0;
      NumberOfSelectedGen=0;
      NumberOfGenerated=0;
      for(Int_t i=fHistSelected_1D_V0Pt[molt]->GetXaxis()->FindBin(NPtV0[j]+0.00001); i<=fHistSelected_1D_V0Pt[molt]->GetXaxis()->FindBin(NPtV0[j+1]-0.00001); i++){
	NumberOfGenerated+=fHistGenerated_1D_V0Pt[molt]->GetBinContent(i);
	NumberOfSelected+=fHistSelected_1D_V0Pt[molt]->GetBinContent(i);
	NumberOfSelectedGen+=fHistSelectedGen_1D_V0Pt[molt]->GetBinContent(i);
      }
      //      cout << "molt " << molt << "PtMin " << j << "NumberSelected " << NumberOfSelected << "NumberGenerated " << NumberOfGenerated << "error " << SetEfficiencyError(NumberOfSelected,NumberOfGenerated) << " first part " << ((Float_t)NumberOfSelected+1)*((Float_t)NumberOfSelected+2)/(NumberOfGenerated+2)/(NumberOfGenerated+3) << " second part " << pow((Float_t)NumberOfSelected+1,2)/pow(NumberOfGenerated+2,2) << endl;
      fHistV0EfficiencyPtBins[molt]->SetBinContent(j+1, NumberOfSelected/NumberOfGenerated);
      fHistV0EfficiencyPtBins[molt]->SetBinError(j+1, SetEfficiencyError(NumberOfSelected,NumberOfGenerated));
      fHistV0EfficiencyGenPtBins[molt]->SetBinContent(j+1, NumberOfSelectedGen/NumberOfGenerated);
      fHistV0EfficiencyGenPtBins[molt]->SetBinError(j+1, SetEfficiencyError(NumberOfSelectedGen,NumberOfGenerated));
      cout <<" Ptbinassoc: " << NPtV0[j] << "  " << 	fHistV0EfficiencyPtBins[molt]->GetBinContent(j+1) << "+-" << 	fHistV0EfficiencyPtBins[molt]->GetBinError(j+1) << " rel: " << fHistV0EfficiencyPtBins[molt]->GetBinError(j+1)/fHistV0EfficiencyPtBins[molt]->GetBinContent(j+1)<< endl;
    }
    canvasUsed->cd(2);
    fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetRangeUser(0,0.25);
    if (type==6)    fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetRangeUser(0,0.5);
    if (ishhCorr)     fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistV0EfficiencyPtBins[molt]->SetMarkerStyle(Marker[molt]);
    fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetTitle("#epsilon_{Assoc}");
    fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetTitleOffset(1);
    fHistV0EfficiencyPtBins[molt]->SetStats(0);
    fHistV0EfficiencyPtBins[molt]->SetLineColor(Color[molt]);
    fHistV0EfficiencyPtBins[molt]->SetMarkerColor(Color[molt]);
    fHistV0EfficiencyPtBins[molt]->Draw("samee");
    if(molt ==nummolt)     legend->Draw();

    fHistV0EfficiencyGenPtBins[molt]->GetYaxis()->SetRangeUser(0,0.25);
    if (type==6)    fHistV0EfficiencyGenPtBins[molt]->GetYaxis()->SetRangeUser(0,0.5);
    if (ishhCorr)     fHistV0EfficiencyGenPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistV0EfficiencyGenPtBins[molt]->SetMarkerStyle(Marker[molt]);
    fHistV0EfficiencyGenPtBins[molt]->GetYaxis()->SetTitle("#epsilon_{Assoc}");
    fHistV0EfficiencyGenPtBins[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyGenPtBins[molt]->GetYaxis()->SetTitleOffset(1);
    fHistV0EfficiencyGenPtBins[molt]->SetStats(0);
    fHistV0EfficiencyGenPtBins[molt]->SetLineColor(Color[molt]);
    fHistV0EfficiencyGenPtBins[molt]->SetMarkerColor(Color[molt]);

    /*
      canvasEffBis[3]->cd();
      fHistV0EfficiencyPtBins[molt]->Draw("same");
      if(molt ==nummolt)     legend->Draw();
    */

    cout << "V0 2D projection in Eta and PtTMax (Pt stands for PtTMax) " << endl;
    if(molt < nummolt){
      fHistSelectedV0PtTMaxEta->GetZaxis()->SetRange(fHistSelectedV0PtTMaxEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtTMaxEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistGeneratedV0PtTMaxEta->GetZaxis()->SetRange(fHistGeneratedV0PtTMaxEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedV0PtTMaxEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    else{
      fHistSelectedV0PtTMaxEta->GetZaxis()->SetRange(0,100);
      fHistGeneratedV0PtTMaxEta->GetZaxis()->SetRange(0,100);
    }

    fHistSelected_2D_V0PtTMaxEta[molt]  = (TH2D*)fHistSelectedV0PtTMaxEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_V0PtTMaxEta[molt] = (TH2D*)fHistGeneratedV0PtTMaxEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistSelected_2D_V0PtTMaxEta[molt] ->SetName("fHistSelected_2D_V0PtTMaxEta_"+ Smolt[molt]);
    fHistGenerated_2D_V0PtTMaxEta[molt]->SetName("fHistGenerated_2D_V0PtTMaxEta_"+ Smolt[molt]);

    fHistV0EfficiencyPtEta[molt]= new TH2D("fHistV0EfficiencyPtEta_"+ Smolt[molt],"fHistV0EfficiencyPtEta_"+ Smolt[molt],fHistSelectedV0PtTMaxEta->GetNbinsX(),fHistSelectedV0PtTMaxEta->GetXaxis()->GetXmin(), fHistSelectedV0PtTMaxEta->GetXaxis()->GetXmax(),fHistSelectedV0PtTMaxEta->GetNbinsY(),fHistSelectedV0PtTMaxEta->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtTMaxEta->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtTMaxEta->GetNbinsY()) );
    fHistV0EfficiencyPtEta[molt]->GetXaxis()->SetTitle("p_{T}^{Trigger, max} (GeV/c)");      
    fHistV0EfficiencyPtEta[molt]->GetYaxis()->SetTitle("#eta_{Assoc}");      

    fHistSelected_2D_V0PtTMaxEta[molt] ->RebinX(1);   
    fHistSelected_2D_V0PtTMaxEta[molt] ->RebinY(10);   

    fHistGenerated_2D_V0PtTMaxEta[molt]->RebinX(1);    
    fHistGenerated_2D_V0PtTMaxEta[molt]->RebinY(10);    

    fHistV0EfficiencyPtEta[molt]->RebinX(1);    
    fHistV0EfficiencyPtEta[molt]->RebinY(10);    

    fHistV0EfficiencyPtEta[molt]->Divide(fHistSelected_2D_V0PtTMaxEta[molt], fHistGenerated_2D_V0PtTMaxEta[molt]); 

    cout << "V0 1D projection in Eta " << endl;
    fHistV0EfficiencyEta[molt]= new TH1D("fHistV0EfficiencyEta_"+ Smolt[molt] , "fHistV0EfficiencyEta_"+ Smolt[molt] ,    fHistGenerated_2D_V0PtTMaxEta[molt]->GetNbinsY(),   fHistGenerated_2D_V0PtTMaxEta[molt]->GetYaxis()->GetBinLowEdge(1),   fHistGenerated_2D_V0PtTMaxEta[molt]->GetYaxis()->GetBinUpEdge(  fHistGenerated_2D_V0PtTMaxEta[molt]->GetNbinsY()) );
    fHistV0EfficiencyEta[molt]->GetXaxis()->SetTitle("#eta_{Assoc}");      
    fHistV0EfficiencyEta[molt]->SetTitle( "fHistV0EfficiencyEta_"+ Smolt[molt] );

    fHistSelected_1D_V0Eta[molt]=(TH1D*)fHistSelected_2D_V0PtTMaxEta[molt]->ProjectionY("fHistSelected_1D_V0Eta_"+ Smolt[molt], fHistSelected_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistSelected_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001));//***

    if (type>=4 && type!=6){
    fHistGenerated_1D_V0Eta_Int1[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxEta[molt]->ProjectionY("fHistGenerated_1D_V0Eta_Int1"+ Smolt[molt],fHistGenerated_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjmin+0.0001) ,fHistGenerated_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001));
    fHistGenerated_1D_V0Eta_Int2[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxEta[molt]->ProjectionY("fHistGenerated_1D_V0Eta_Int2"+ Smolt[molt],fHistGenerated_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistGenerated_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(-ptjmin-0.0001));
    fHistGenerated_1D_V0Eta[molt]=(TH1D*)fHistGenerated_1D_V0Eta_Int2[molt]->Clone("fHistGenerated_1D_V0Eta_"+ Smolt[molt]);
    fHistGenerated_1D_V0Eta[molt]->Add(fHistGenerated_1D_V0Eta_Int1[molt]);
    }
    else {
    fHistGenerated_1D_V0Eta[molt]=(TH1D*)fHistGenerated_2D_V0PtTMaxEta[molt]->ProjectionY("fHistGenerated_1D_V0Eta_"+ Smolt[molt],fHistGenerated_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjminGen+0.0001) ,fHistGenerated_2D_V0PtTMaxEta[molt]->GetXaxis()->FindBin(ptjmaxGen-0.0001));
    }

    fHistV0EfficiencyEta[molt] ->Divide(  fHistSelected_1D_V0Eta[molt],  fHistGenerated_1D_V0Eta[molt]);
    canvasEff->cd(6);
    fHistV0EfficiencyEta[molt]->GetYaxis()->SetRangeUser(0.001,0.2);
    if (type==6)     fHistV0EfficiencyEta[molt]->GetYaxis()->SetRangeUser(0.001,1);
    if (ishhCorr)   fHistV0EfficiencyEta[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistV0EfficiencyEta[molt]->GetYaxis()->SetTitle("#epsilon_{Assoc}");
    fHistV0EfficiencyEta[molt]->GetXaxis()->SetTitleOffset(0.85);
    fHistV0EfficiencyEta[molt]->GetXaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyEta[molt]->GetYaxis()->SetTitleOffset(1);
    fHistV0EfficiencyEta[molt]->GetYaxis()->SetTitleSize(0.05);
    fHistV0EfficiencyEta[molt]->SetStats(0);
    fHistV0EfficiencyEta[molt]->GetYaxis()->SetRangeUser(0,0.24);
    if (ishhCorr)   fHistV0EfficiencyEta[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistV0EfficiencyEta[molt]->SetMarkerStyle(Marker[molt]);
    fHistV0EfficiencyEta[molt]->SetLineColor(Color[molt]);
    fHistV0EfficiencyEta[molt]->SetMarkerColor(Color[molt]);
    fHistV0EfficiencyEta[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();

    /*  
	canvasEffBis[5]->cd();
	fHistV0EfficiencyEta[molt]->Draw("same");
	if(molt ==nummolt)     legend->Draw();
    */

    if (fHistSelectedV0PtMass && fHistReconstructedV0PtMass) {
    //per efficienza selezioni V0 
    if(molt < nummolt){
      fHistSelectedV0PtMass->GetZaxis()->SetRange(fHistSelectedV0PtMass->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtMass->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
      fHistReconstructedV0PtMass->GetZaxis()->SetRange(fHistReconstructedV0PtMass->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistReconstructedV0PtMass->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    else{
      fHistSelectedV0PtMass->GetZaxis()->SetRange(0,100);
      fHistReconstructedV0PtMass->GetZaxis()->SetRange(0,100);
    }
 

    fHistSelectedMass_2D[molt] = (TH2D*)fHistSelectedV0PtMass->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistRecoMass_2D[molt] = (TH2D*)fHistReconstructedV0PtMass->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistV0EfficiencyReco[molt]= new TH2D("fHistV0EfficiencyReco_"+ Smolt[molt],"fHistV0EfficiencyReco_"+ Smolt[molt],fHistSelectedV0PtMass->GetNbinsX(),fHistSelectedV0PtMass->GetXaxis()->GetXmin(), fHistSelectedV0PtMass->GetXaxis()->GetXmax(),fHistSelectedV0PtMass->GetNbinsY(),fHistSelectedV0PtMass->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtMass->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtMass->GetNbinsY()) );
    fHistV0EfficiencyReco[molt]->GetXaxis()->SetTitle("Mass");      
    fHistV0EfficiencyReco[molt]->GetYaxis()->SetTitle("p_{T} (GeV/c)");      
    fHistV0EfficiencyReco[molt]->Divide(fHistSelectedMass_2D[molt], fHistRecoMass_2D[molt]);

 
    fHistV0EfficiencyRecoPt[molt]= new TH1D("fHistV0EfficiencyRecoPt_"+ Smolt[molt] , "fHistV0EfficiencyRecoPt_"+ Smolt[molt],  fHistV0EfficiencyReco[molt]->GetNbinsY(), fHistV0EfficiencyReco[molt]->GetYaxis()->GetBinLowEdge(1), fHistV0EfficiencyReco[molt]->GetYaxis()->GetBinUpEdge(fHistV0EfficiencyReco[molt]->GetNbinsY()));
    fHistSelectedMass_1D[molt]= (TH1D*)fHistSelectedMass_2D[molt]->ProjectionY("fHistSelectedMass_1D_"+ Smolt[molt]);
    fHistRecoMass_1D[molt]= (TH1D*)fHistRecoMass_2D[molt]->ProjectionY("fHistRecoMass_1D_"+ Smolt[molt]);  
    fHistV0EfficiencyRecoPt[molt] ->Divide(fHistSelectedMass_1D[molt], fHistRecoMass_1D[molt]);
    }

    cout << "resoltiuon histos" << endl;
    //**********************resolution histograms*****************************
    for(Int_t m=0; m< 3; m++){
      for(Int_t t=0; t< 2; t++){
	/* version to be used if resolution hiostograms are 3D, comment if not */
	   if (ResoHisto=="3D"){
	fHistResolution3D[m][t]->GetZaxis()->SetRange(fHistResolution3D[m][t]->GetZaxis()->FindBin(PtTrigMin+0.001),fHistResolution3D[m][t]->GetZaxis()->FindBin(PtTrigMax-0.001) );
	fHistResolution2D[m][t]=(TH2D*)fHistResolution3D[m][t]->Project3D("yxo");
	nameRes_2D[molt][m][t]="fHistResolution2D" + TorV[t] + Var[m] + Smolt[molt];
	fHistResolution2D[m][t]->SetName(nameRes_2D[molt][m][t]);
	   }
	// version to be used if resolution hiostograms are 3D, comment if not */
      
	nameRes_1D[molt][m][t]="fHistResolution" + TorV[t] + Var[m] + Smolt[molt];
	if (ResoHisto=="3D"){
	  if(molt< nummolt){
	    fHistResolution_1D[molt][m][t]=(TH1D*)fHistResolution2D[m][t]->ProjectionX(nameRes_1D[molt][m][t], fHistResolution2D[m][t]->GetYaxis()->FindBin(Nmolt[molt]+0.0001), fHistResolution2D[m][t]->GetYaxis()->FindBin(Nmolt[molt+1]-0.0001));
	  }
	  else {
	    fHistResolution_1D[molt][m][t]=(TH1D*)fHistResolution2D[m][t]->ProjectionX(nameRes_1D[molt][m][t]);
	  }
	}
	else {
	  if (t==0)  fHistResolution_1D[molt][m][t]=(TH1D*)fHistResolution2D[m][t]->ProjectionX(nameRes_1D[molt][m][t], fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjminBisFixed[molt]+0.0001), fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjmax-0.0001));
	  else if (t==1){
	    if (type>=4){
	      fHistResolution_1D_Int1[molt][m][t]=(TH1D*)fHistResolution2D[m][t]->ProjectionX(nameRes_1D[molt][m][t]+"_Int1", fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjminBisFixed[molt]+0.0001), fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjmaxRes[molt]-0.0001));
	      fHistResolution_1D_Int2[molt][m][t]=(TH1D*)fHistResolution2D[m][t]->ProjectionX(nameRes_1D[molt][m][t]+"_Int2", fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjminBis[molt]+0.0001), fHistResolution2D[m][t]->GetYaxis()->FindBin(-ptjminBisFixed[molt]-0.0001));
	      fHistResolution_1D[molt][m][t]=(TH1D*)  fHistResolution_1D_Int2[molt][m][t]->Clone(nameRes_1D[molt][m][t]);
	      fHistResolution_1D[molt][m][t]->Add(  fHistResolution_1D_Int1[molt][m][t]);
	    }
	    else {
	      fHistResolution_1D[molt][m][t]=(TH1D*)fHistResolution2D[m][t]->ProjectionX(nameRes_1D[molt][m][t], fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjminBis[molt]+0.0001), fHistResolution2D[m][t]->GetYaxis()->FindBin(ptjmaxRes[molt]-0.0001));
	    }
	  }
	}

	Int_t normal = fHistResolution_1D[molt][m][t]->GetEntries();
	fHistResolution_1D[molt][m][t]->Scale(1./normal);
	if (t==0)      canvasRes->cd(m+1);
	else  canvasRes->cd(m+4);
	fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.4);
	fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.2, 0.2); //was 0.05
	if (type==0 || type==1 || type==4) {
	  if (t==1){
	    if (m==0) fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.2, 0.2);
	    if (m==1 || m==2) {
	      fHistResolution_1D[molt][m][t]->Rebin(2);
	      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.03, 0.03);
	      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.5);
	    }
	  }
	  else if (t==0){
	    if (m==0) fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.2, 0.2);
	    if (m==1 || m==2) {
	      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.03, 0.03);
	      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.5);
	    }
	  }
	}
	if (type==6) {
	  if (t==1){
	    if (m==0){
	      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.2, 0.2);
	      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.5);
	    }
	    else if (m==1 || m==2) {
	      fHistResolution_1D[molt][m][t]->Rebin(2);
	      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.03, 0.03);
	      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.5);
	    }
	  }
	  else if (t==0){
	    if (m==0){
	      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.2, 0.2);
	      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.5);
	    }
	    else   if (m==1 || m==2 ) {
	      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.03, 0.03);
	      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.5);
	    }
	  }
	}

	if (m==0 && t==0){
	  //	  fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.5, 0.5);
	  //	  fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 0.1);
	}
	fHistResolution_1D[molt][m][t]->SetMarkerStyle(Marker[molt]);
	fHistResolution_1D[molt][m][t]->SetLineColor(Color[molt]);
	fHistResolution_1D[molt][m][t]->SetMarkerColor(Color[molt]);
	fHistResolution_1D[molt][m][t]->Draw("same");
	if(ResoHisto=="3D" && molt ==nummolt)     legend->Draw();
	if(ResoHisto=="2D" && molt ==nummolt)     legendPtMin->Draw();
	//canvasRes->Close();
	//	cout << "molt: " << molt << " sto riempiendo histo " <<  nameRes_1D[molt][m][t]<< endl;
      }
    }

    //***************contamination factor for trigger and V0 in multiplicity intervals***************

    HistContTrigger[molt]=(TH2D*)list3->FindObject(Form("fHistPrimaryTrigger_%i_cut%i", molt, sysTrigger)); //histo tipo di trigger (primario, secondario) vs pT trigger
    //    cout << "deficit " << endl;
    HistContV0PtTMax[molt]=(TH3D*)fileinputSel->Get(Form("fHistPrimaryV0_%i_cut%i", molt, sysV0hh));//histo tipo di V0 (primario, secondario) vs pT V0

    HistInt[molt]=HistContTrigger[molt]->ProjectionX("", HistContTrigger[molt]->GetYaxis()->FindBin(ptjmin+0.00001),HistContTrigger[molt]->GetYaxis()->FindBin(ptjmax-0.0001));
    if((HistInt[molt]->GetBinContent(1) + HistInt[molt]->GetBinContent(2) + HistInt[molt]->GetBinContent(3)+ HistInt[molt]->GetBinContent(4)) != 0){
      ContTrigger[molt]=1. - (Double_t)(HistInt[molt]->GetBinContent(1))/(HistInt[molt]->GetBinContent(1) + HistInt[molt]->GetBinContent(2) + HistInt[molt]->GetBinContent(3)+ HistInt[molt]->GetBinContent(4));
      HistContTriggerMolt->SetBinContent(molt+1, ContTrigger[molt]);
    }
    else     HistContTriggerMolt->SetBinContent(molt+1, 0);

    canvasUsed->cd(3);
    HistContTriggerMolt->GetYaxis()->SetRangeUser(0,0.03);
    HistContTriggerMolt->GetYaxis()->SetTitle("C_{Trigg}");
    HistContTriggerMolt->GetYaxis()->SetTitleSize(0.5);
    HistContTriggerMolt->GetYaxis()->SetTitleOffset(1.2);
    HistContTriggerMolt->GetXaxis()->SetTitle("Multiplicity class");
    HistContTriggerMolt->GetXaxis()->SetTitleSize(0.5);
    HistContTriggerMolt->GetXaxis()->SetTitleOffset(1.2);
    HistContTriggerMolt->SetMarkerStyle(ColorSysTrigger[sysTrigger]);
    HistContTriggerMolt->SetLineColor(ColorSysTrigger[sysTrigger]);
    HistContTriggerMolt->SetMarkerColor(ColorSysTrigger[sysTrigger]);
    // HistContTriggerMolt->Draw("e");
    HistContTriggerMolt->Draw();
    
    HistContV0PtTMax[molt]->GetZaxis()->SetRange(    HistContV0PtTMax[molt]->GetZaxis()->FindBin(PtTrigMin+0.001),    HistContV0PtTMax[molt]->GetZaxis()->FindBin(PtTrigMax-0.0001));
    HistContV0[molt]=(TH2D*)HistContV0PtTMax[molt]->Project3D("yxo");
    HistContV0[molt]->SetName("HistContV0_"+Smolt[molt]);
    HistInt[molt]=      HistContV0[molt]->ProjectionX();
    ContV0[molt]=1. - (Double_t)(HistInt[molt]->GetBinContent(1))/(HistInt[molt]->GetBinContent(1) + HistInt[molt]->GetBinContent(2) + HistInt[molt]->GetBinContent(3)+ HistInt[molt]->GetBinContent(4));

    cout << " cont trigger " << ContTrigger[molt] << " cont v0 " << ContV0[molt] << endl;
      
    HistContTriggerPt[molt]=new TH1D("HistContTriggerPt_"+Smolt[molt], Form("HistContTriggerPt_%i", molt),HistContTrigger[molt]->GetNbinsY(),HistContTrigger[molt]->GetYaxis()->GetBinLowEdge(1), HistContTrigger[molt]->GetYaxis()->GetBinUpEdge(HistContTrigger[molt]->GetNbinsY()) );
    HistContTriggerPt[molt]->SetTitle("HistContTriggerPt_"+Smolt[molt]);
    HistContTriggerPt[molt]->GetXaxis()->SetTitle("p_T^{Trigg}");
    HistContTriggerPt[molt]->GetXaxis()->SetTitleSize(1.5);
    HistContTriggerPt[molt]->GetYaxis()->SetTitle("C factor");
    HistContTriggerPt[molt]->SetStats(0);

    HistContV0Pt[molt]=new TH1D("HistContV0Pt_"+Smolt[molt], Form("HistContV0Pt_%i", molt),HistContV0[molt]->GetNbinsY(),HistContV0[molt]->GetYaxis()->GetBinLowEdge(1), HistContV0[molt]->GetYaxis()->GetBinUpEdge(HistContV0[molt]->GetNbinsY()) );
    HistContV0Pt[molt]->SetTitle("HistContV0Pt_"+Smolt[molt]);
    HistContV0Pt[molt]->GetXaxis()->SetTitle("p_T^{Assoc}");
    HistContV0Pt[molt]->GetXaxis()->SetTitleSize(1.5);
    HistContV0Pt[molt]->GetYaxis()->SetTitle("C factor");

    HistContV0PtBins[molt]=new TH1D("HistContV0PtBins_"+Smolt[molt], Form("HistContV0PtBins_%i", molt),numPtV0, NPtV0 );
    HistContV0PtBins[molt]->SetTitle("");
    HistContV0PtBins[molt]->GetXaxis()->SetTitle("p_{T}^{Assoc}");
    HistContV0PtBins[molt]->GetYaxis()->SetTitle("C_{Assoc}");
    HistContV0PtBins[molt]->GetXaxis()->SetTitleSize(0.045);
    HistContV0PtBins[molt]->GetYaxis()->SetTitleSize(0.045);
    HistContV0PtBins[molt]->GetYaxis()->SetTitleOffset(1);
    HistContV0PtBins[molt]->SetStats(0);
    /*
      canvascontv0->cd();
      HistContV0PtBins[molt]->Draw("same");
      if (molt==nummolt)    legend->Draw("same");
    */
    Double_t denom=0;
    Double_t num=0;

    //histo contamination factor for trigger particles in pt bins of trigger particle
    for(Int_t pt=0; pt<HistContTrigger[molt]->GetNbinsY(); pt++ ){
      num= (Double_t)(HistContTrigger[molt]->GetBinContent(1, pt+1)); //number of primary trigger particles
      denom= (HistContTrigger[molt]->GetBinContent(1, pt+1) + HistContTrigger[molt]->GetBinContent(2, pt+1) + HistContTrigger[molt]->GetBinContent(3, pt+1)+ HistContTrigger[molt]->GetBinContent(4, pt+1)); //total number of trigger particles
      if (denom!=0){
	ContTriggerInt=1. - num/denom;
	//	ContTriggerIntError=sqrt(pow(ContTriggerInt,2)*(denom-num)+pow((num-denom)/pow(denom,2),2)*num);
	//ContTriggerIntError=sqrt(pow((num-denom)/pow(denom,2),2)*num);
	HistContTriggerPt[molt]->SetBinContent(pt+1, ContTriggerInt);
	//	HistContTriggerPt[molt]->SetBinError(pt+1, ContTriggerIntError);
	//	cout << "pt bin " << pt << " cont factor trigger " << ContTriggerInt <<" +- " << ContTriggerIntError <<endl;
      }
    }


    //histo contamination factor for V0 particles in pt bins of V0 particle
    for(Int_t pt=0; pt<HistContV0[molt]->GetNbinsY(); pt++ ){
      num= (Double_t)(HistContV0[molt]->GetBinContent(1, pt+1)); //number of primary V0 
      denom= (HistContV0[molt]->GetBinContent(1, pt+1) + HistContV0[molt]->GetBinContent(2, pt+1) + HistContV0[molt]->GetBinContent(3, pt+1)+ HistContV0[molt]->GetBinContent(4, pt+1)); //total number of V0
      if(denom!=0){
	ContV0Int =1. - num/denom;
	//ContV0IntError=sqrt(pow(ContV0Int,2)*(denom-num)+pow((num-denom)/pow(denom,2),2)*num);
	HistContV0Pt[molt]->SetBinContent(pt+1, ContV0Int);
	//HistContV0Pt[molt]->SetBinContent(pt+1, ContV0IntError);
	//	cout << "pt bin " << pt << " cont factor V0 " << ContV0Int <<" +- " << ContV0IntError << endl;
      }
      else{
	HistContV0Pt[molt]->SetBinContent(pt+1, 0);
	HistContV0Pt[molt]->SetBinContent(pt+1, 0);
      } 
    }

    cout << "\n\n contamination factor for V0 particles in Pt bins " << endl;
    //histo contamination factor for V0 particles in pt bins of V0 particle (bins used in analysis)
    for(Int_t binpt=0; binpt<numPtV0Max; binpt++){
      num=0;
      denom=0;
      for(Int_t pt=HistContV0[molt]->GetYaxis()->FindBin(NPtV0[binpt]+0.0001); pt<HistContV0[molt]->GetYaxis()->FindBin(NPtV0[binpt+1]-0.0001); pt++){
	num+= (Double_t)(HistContV0[molt]->GetBinContent(1, pt+1)); //number of primary V0 
	denom+= (HistContV0[molt]->GetBinContent(1, pt+1) + HistContV0[molt]->GetBinContent(2, pt+1) + HistContV0[molt]->GetBinContent(3, pt+1)+ HistContV0[molt]->GetBinContent(4, pt+1)); //total number of V0
      }
      if(denom!=0){
	ContV0Int =1. - num/denom;
	//	ContV0IntError=sqrt(pow(ContV0Int,2)*(denom-num)+pow((num-denom)/pow(denom,2),2)*num);
	HistContV0PtBins[molt]->SetBinContent(binpt+1, ContV0Int);
	HistContV0PtBins[molt]->SetBinError(binpt+1, 0);
	//HistContV0PtBins[molt]->SetBinContent(binpt+1, ContV0IntError);
	//	cout << "pt bin " << binpt << " cont factor V0 " << ContV0Int <<" +- " << ContV0IntError << endl;
      }
      else{
	HistContV0PtBins[molt]->SetBinContent(binpt+1, 0);
	HistContV0PtBins[molt]->SetBinError(binpt+1, 0);

	
      }
    }


    canvasUsed->cd(4);
    HistContV0PtBins[molt]->GetYaxis()->SetRangeUser(0,0.005);
    if (ishhCorr)     HistContV0PtBins[molt]->GetYaxis()->SetRangeUser(0,0.05);
    HistContV0PtBins[molt]->SetMarkerStyle(Marker[molt]);
    HistContV0PtBins[molt]->SetLineColor(Color[molt]);
    HistContV0PtBins[molt]->SetMarkerColor(Color[molt]);
    HistContV0PtBins[molt]->Draw("samee");
    if(molt ==nummolt)     legend->Draw();

    canvasCont->cd(1);
    HistContTriggerPt[molt]->Rebin(2);
    HistContTriggerPt[molt]->GetYaxis()->SetRangeUser(0, 0.15);
    HistContTriggerPt[molt]->SetMarkerStyle(Marker[molt]);
    HistContTriggerPt[molt]->SetLineColor(Color[molt]);
    HistContTriggerPt[molt]->SetMarkerColor(Color[molt]);
    HistContTriggerPt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();

    //    canvasUsed->Close();
    canvasCont->cd(2);
    HistContV0Pt[molt]->Rebin(2);
    HistContV0Pt[molt]->GetYaxis()->SetRangeUser(0, 0.15);
    HistContV0Pt[molt]->SetMarkerStyle(Marker[molt]);
    HistContV0Pt[molt]->SetLineColor(Color[molt]);
    HistContV0Pt[molt]->SetMarkerColor(Color[molt]);
    HistContV0Pt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();
    //    canvasCont->Close();
  
  }


  TProfile *fHistTriggerPtRecovsPtGen_PtBinspfx;
  TH2F* fHistTriggerPtRecovsPtGen;  
  TProfile *fHistTriggerPtRecovsPtGen_pfx;
  TH2F* fHistAssocPtRecovsPtGen;
  TH2F* fHistAssocPtRecovsPtGenInt1;
  TH2F* fHistAssocPtRecovsPtGenInt2;
  TProfile *fHistAssocPtRecovsPtGen_pfx;
  TH2F* fHistTriggerPtRecovsPtGen_PtBins;
  TH2F* fHistAssocPtRecovsPtGen_PtBins;
  TProfile *fHistAssocPtRecovsPtGen_PtBinspfx;
if (listRisoluzione){
  //********TProfile (used only for tasks investigating the  origin of strange behaviour of resolution)**************
  fHistTriggerPtRecovsPtGen = (TH2F*) listRisoluzione->FindObject("fHistTriggerPtRecovsPtGen");
  if (!fHistTriggerPtRecovsPtGen) {cout << "manca istogramma pt reco vs pt gen , non posso proseguire " << endl; return;}
  fHistTriggerPtRecovsPtGen_pfx = (TProfile*) fHistTriggerPtRecovsPtGen->ProfileX("fHistTriggerPtRecovsPtGen_pfx");
  if (type==0 || type==2)  fHistAssocPtRecovsPtGen = (TH2F*) listRisoluzione->FindObject("fHistAssocPtRecovsPtGenNeg");
  else if (type==1 || type==3)  fHistAssocPtRecovsPtGen = (TH2F*) listRisoluzione->FindObject("fHistAssocPtRecovsPtGenPos");
  else  if (type==4 || type==5) {
    fHistAssocPtRecovsPtGenInt1 = (TH2F*) listRisoluzione->FindObject("fHistAssocPtRecovsPtGenPos");
    fHistAssocPtRecovsPtGenInt2 = (TH2F*) listRisoluzione->FindObject("fHistAssocPtRecovsPtGenNeg");
    fHistAssocPtRecovsPtGen = (TH2F*)     fHistAssocPtRecovsPtGenInt2->Clone("fHistAssocPtRecovsPtGen");
    fHistAssocPtRecovsPtGen->Add(fHistAssocPtRecovsPtGenInt1);
  }
  else  if (type==6 || type==7) {
    fHistAssocPtRecovsPtGen= (TH2F*) listRisoluzione->FindObject("fHistAssocPtRecovsPtGen");
  }
  if (!fHistAssocPtRecovsPtGen) {cout << "manca istogramma pt reco vs pt gen per associati, non posso proseguire " << endl; return;}
  fHistAssocPtRecovsPtGen_pfx = (TProfile*) fHistAssocPtRecovsPtGen->ProfileX("fHistAssocPtRecovsPtGen_pfx");

  // TF1* rettaprofile = new TF1("pol1", "pol1", 0,30);
  // rettaprofile->FixParameter(0,0);
  // rettaprofile->FixParameter(0,1);
  // fHistTriggerPtRecovsPtGen_PtBinspfx->Fit(rettaprofile, "RB+");

  //Large bins TProfile
  fHistTriggerPtRecovsPtGen_PtBins = new TH2F("fHistTriggerPtRecovsPtGen_PtBins", "fHistTriggerPtRecovsPtGen_PtBins", numPtTrigger, NPtTrigger, numPtTrigger, NPtTrigger);
  Int_t CountForTProfile=0;
  for (Int_t pTTrigger=0; pTTrigger<numPtTrigger; pTTrigger++){
  for (Int_t pTTrigger1=0; pTTrigger1<numPtTrigger; pTTrigger1++){
  CountForTProfile=0;   
  //    cout << "\n\npTTrigger " << pTTrigger << " " << pTTrigger1<< " " << CountForTProfile<<  endl;
  for(Int_t i=fHistTriggerPtRecovsPtGen->GetXaxis()->FindBin(NPtTrigger[pTTrigger]+0.001); i<=fHistTriggerPtRecovsPtGen->GetXaxis()->FindBin(NPtTrigger[pTTrigger+1]-0.001); i++ ){
  for(Int_t j=fHistTriggerPtRecovsPtGen->GetXaxis()->FindBin(NPtTrigger[pTTrigger1]+0.001); j<=fHistTriggerPtRecovsPtGen->GetXaxis()->FindBin(NPtTrigger[pTTrigger1+1]-0.001); j++ ){
  // cout << " i j" << i << "  " << j <<"  " << CountForTProfile << endl;
  CountForTProfile +=fHistTriggerPtRecovsPtGen->GetBinContent(fHistTriggerPtRecovsPtGen->GetBin(i,j));
	  
  }
  }
  //  cout << "\n\n\npTTrigger " << pTTrigger << " " << pTTrigger1<< endl;
  //  cout << "low and up edge of global bin (x direction) "<< fHistTriggerPtRecovsPtGen_PtBins->GetXaxis()->GetBinLowEdge(pTTrigger+1)<< "  " << fHistTriggerPtRecovsPtGen_PtBins->GetXaxis()->GetBinUpEdge(pTTrigger+1)<<endl;
  //  cout << "global bin number " << fHistTriggerPtRecovsPtGen_PtBins->GetBin(pTTrigger+1,pTTrigger1+1) << " and its content " << CountForTProfile << endl;
  //  fHistTriggerPtRecovsPtGen_PtBins->AddBinContent(fHistTriggerPtRecovsPtGen_PtBins->GetBin(pTTrigger+1,pTTrigger1+1),   CountForTProfile );
  fHistTriggerPtRecovsPtGen_PtBins->Fill(NPtTrigger[pTTrigger]+0.001, NPtTrigger[pTTrigger1]+0.001, CountForTProfile );
  }
  }

  fHistTriggerPtRecovsPtGen_PtBinspfx = (TProfile*) fHistTriggerPtRecovsPtGen_PtBins->ProfileX("fHistTriggerPtRecovsPtGen_PtBinspfx");  


  //Large bins TProfile for assoc particles
  fHistAssocPtRecovsPtGen_PtBins = new TH2F("fHistAssocPtRecovsPtGen_PtBins", "fHistAssocPtRecovsPtGen_PtBins", numPtV0, NPtV0, numPtV0, NPtV0);
  for (Int_t pTAssoc=0; pTAssoc<numPtV0Max; pTAssoc++){
  for (Int_t pTAssoc1=0; pTAssoc1<numPtV0Max; pTAssoc1++){
  CountForTProfile=0;   
  //    cout << "\n\npTAssoc " << pTAssoc << " " << pTAssoc1<< " " << CountForTProfile<<  endl;
  for(Int_t i=fHistAssocPtRecovsPtGen->GetXaxis()->FindBin(NPtV0[pTAssoc]+0.001); i<=fHistAssocPtRecovsPtGen->GetXaxis()->FindBin(NPtV0[pTAssoc+1]-0.001); i++ ){
  for(Int_t j=fHistAssocPtRecovsPtGen->GetXaxis()->FindBin(NPtV0[pTAssoc1]+0.001); j<=fHistAssocPtRecovsPtGen->GetXaxis()->FindBin(NPtV0[pTAssoc1+1]-0.001); j++ ){
  // cout << " i j" << i << "  " << j <<"  " << CountForTProfile << endl;
  CountForTProfile +=fHistAssocPtRecovsPtGen->GetBinContent(fHistAssocPtRecovsPtGen->GetBin(i,j));
	  
  }
  }
  //  cout << "\n\n\npTAssoc " << pTAssoc << " " << pTAssoc1<< endl;
  //  cout << "low and up edge of global bin (x direction) "<< fHistAssocPtRecovsPtGen_PtBins->GetXaxis()->GetBinLowEdge(pTAssoc+1)<< "  " << fHistAssocPtRecovsPtGen_PtBins->GetXaxis()->GetBinUpEdge(pTAssoc+1)<<endl;
  cout << "global bin number " << fHistAssocPtRecovsPtGen_PtBins->GetBin(pTAssoc+1,pTAssoc1+1) << " and its content " << CountForTProfile << endl;
  //  fHistAssocPtRecovsPtGen_PtBins->AddBinContent(fHistAssocPtRecovsPtGen_PtBins->GetBin(pTAssoc+1,pTAssoc1+1),   CountForTProfile );
  fHistAssocPtRecovsPtGen_PtBins->Fill(NPtV0[pTAssoc]+0.001, NPtV0[pTAssoc1]+0.001, CountForTProfile );
  }
  }

  fHistAssocPtRecovsPtGen_PtBinspfx = (TProfile*) fHistAssocPtRecovsPtGen_PtBins->ProfileX("fHistAssocPtRecovsPtGen_PtBinspfx");  

  //********End of TProfile (used only for tasks investigating the  origin of strange behaviour of resolution)**************
  }
 
  cout << "\n \n sto per scriver sul file efficienza " << endl;
  
  TFile *fileoutbis;  
  fileoutbis = new TFile(PathOut2, "RECREATE");
  // fHistSelected_2D_TriggerPtPhi->Write();
  // fHistGenerated_2D_TriggerPtPhi->Write();

  /*
    fHistTriggerEfficiencyPtPhiMolt->Write();
    fHistTriggerEfficiencyPtEtaMolt->Write();
    fHistSelected_1D_TriggerPt->Write();
    fHistGenerated_1D_TriggerPt->Write();
    fHistTriggerEfficiencyPtMolt->Write();
    fHistSelected_1D_TriggerPhi->Write();
    fHistGenerated_1D_TriggerPhi->Write();
    fHistTriggerEfficiencyPhiMolt->Write();
    fHistSelected_1D_TriggerEta->Write();
    fHistGenerated_1D_TriggerEta->Write();
    fHistTriggerEfficiencyEtaMolt->Write();
    
    fHistV0EfficiencyPtPhiMolt->Write();
    fHistV0EfficiencyPtEtaMolt->Write();
    fHistSelected_1D_V0Pt->Write();
    fHistGenerated_1D_V0Pt->Write();
    fHistV0EfficiencyPtMolt ->Write();
    fHistSelected_1D_V0Phi->Write();
    fHistGenerated_1D_V0Phi->Write();
    fHistV0EfficiencyPhiMolt->Write();
    fHistSelected_1D_V0Eta->Write();
    fHistGenerated_1D_V0Eta->Write();
    fHistV0EfficiencyEtaMolt->Write();

    fHistSelectedMass_1D->Write();
    fHistRecoMass_1D->Write();
    fHistV0EfficiencyRecoMolt->Write();
    fHistV0EfficiencyRecoPtMolt->Write();
  */

  /*Resolution studies*/
  if (listRisoluzione){
    fHistTriggerPtRecovsPtGen_pfx->Write();
    fHistTriggerPtRecovsPtGen_PtBinspfx->Write();
    fHistTriggerPtRecovsPtGen->Write();
    fHistTriggerPtRecovsPtGen_PtBins->Write();
    fHistAssocPtRecovsPtGen_pfx->Write();
    fHistAssocPtRecovsPtGen_PtBinspfx->Write();
    fHistAssocPtRecovsPtGen->Write();
    fHistAssocPtRecovsPtGen_PtBins->Write();
  }  //end of resolution studies

  HistoTriggerEfficiency->Write();
  HistContTriggerMolt->Write();   
  fileoutbis->WriteTObject( canvasEff);
  fileoutbis->WriteTObject( canvasRes);
  fileoutbis->WriteTObject( canvasCont);
  fileoutbis->WriteTObject( canvasUsed);

  for(Int_t m=0; m<nummolt+1; m++){
    //    if (m==0) continue;


    fHistTriggerEfficiencyPtPhi[m]->Write();
    fHistTriggerEfficiencyPtEta[m]->Write();
    fHistSelected_1D_TriggerPt[m]->Write();
    fHistSelectedGen_1D_TriggerPt[m]->Write();
    fHistGenerated_1D_TriggerPt[m]->Write();
    fHistTriggerEfficiencyPt[m]->Write();
    fHistTriggerEfficiencyGenPt[m]->Write();
    fHistTriggerSelectedPtBins[m] ->Write();
    fHistTriggerSelectedGenPtBins[m] ->Write();
    fHistTriggerGeneratedPtBins[m] ->Write();
    fHistTriggerEfficiencyPtBins[m] ->Write();
    fHistTriggerEfficiencyGenPtBins[m]->Write();
    fHistSelected_1D_TriggerPhi[m]->Write();
    fHistGenerated_1D_TriggerPhi[m]->Write();
    fHistTriggerEfficiencyPhi[m]->Write();
    fHistSelected_1D_TriggerEta[m]->Write();
    fHistGenerated_1D_TriggerEta[m]->Write();
    fHistTriggerEfficiencyEta[m]->Write();
      
    fHistSelected_2D_V0PtPtTMax[m]->Write();
    fHistGenerated_2D_V0PtPtTMax[m]->Write();
    fHistSelected_2D_V0PtTMaxPhi[m]->Write();
    fHistGenerated_2D_V0PtTMaxPhi[m]->Write();
    fHistV0EfficiencyPtPhi[m]->Write();
    fHistV0EfficiencyPtEta[m]->Write();
    fHistSelected_1D_V0Pt[m]->Write();
    fHistSelectedGen_1D_V0Pt[m]->Write();
    fHistGenerated_1D_V0Pt[m]->Write();
    fHistV0EfficiencyPt[m]->Write();
    fHistSelected_1D_V0Phi[m]->Write();
    fHistGenerated_1D_V0Phi[m]->Write();
    fHistV0EfficiencyPhi[m]->Write();
    fHistSelected_1D_V0Eta[m]->Write();
    fHistGenerated_1D_V0Eta[m]->Write();
    fHistV0EfficiencyEta[m]->Write();
    
    fHistV0EfficiencyGenPtBins[m]->Write();
    fHistV0EfficiencyPtBins[m] ->Write();
      
    //    fHistV0EfficiencyReco[m]->Write();
    //    fHistV0EfficiencyRecoPt[m]->Write();
     
  
    HistContV0[m]->Write();
    HistContTrigger[m]->Write();
    HistContV0Pt[m]->Write();
    HistContV0PtBins[m]->Write();
    HistContTriggerPt[m]->Write();
  

    cout << m << endl;
    for(Int_t y=0; y< 3; y++){
      for(Int_t t=0; t< 2; t++){
	fHistResolution_1D[m][y][t]->Write();
      }
    }

  }
  fileoutbis->Close();
    
  TCanvas *TriggerEfficiencyC = new TCanvas("triggeredff", "triggereff", 800, 600);
  TriggerEfficiencyC->cd();
  HistoTriggerEfficiency->GetXaxis()->SetTitle("Multiplicity class");
  HistoTriggerEfficiency->GetYaxis()->SetTitle("#epsilon_{Trigg}");
  HistoTriggerEfficiency->GetYaxis()->SetTitleSize(0.05);
  HistoTriggerEfficiency->GetYaxis()->SetTitleOffset(0.8);
  HistoTriggerEfficiency->GetXaxis()->SetTitleSize(0.04);
  HistoTriggerEfficiency->GetXaxis()->SetTitleOffset(1.2);
  HistoTriggerEfficiency->SetLineColor(1);
  HistoTriggerEfficiency->SetMarkerColor(1);
  HistoTriggerEfficiency->Draw();
  TriggerEfficiencyC->Close();

  TCanvas *TriggerEfficiencyD = new TCanvas("triggeredffd", "trigger cont factor", 800, 600);
  TriggerEfficiencyD->cd();
  HistContTriggerMolt->GetYaxis()->SetRangeUser(0,0.015);
  HistContTriggerMolt->GetYaxis()->SetTitle("C_{Trigg}");
  HistContTriggerMolt->GetYaxis()->SetTitleSize(0.03);
  HistContTriggerMolt->GetYaxis()->SetTitleOffset(0.9);
  HistContTriggerMolt->GetXaxis()->SetTitle("Multiplicity class");
  HistContTriggerMolt->GetXaxis()->SetTitleSize(0.04);
  HistContTriggerMolt->GetXaxis()->SetTitleOffset(1.2);

  HistContTriggerMolt->SetLineColor(1);
  HistContTriggerMolt->SetMarkerColor(1);

  HistContTriggerMolt->Draw();
  TriggerEfficiencyD->Close();

  cout << "*************Resolution information********************" << endl;

  for(Int_t t=0; t< 2; t++){
    if (t==0) cout << "Trigger: " << endl;
    if (t==1) cout << "V0: " << endl;
    for(Int_t m=0; m< 3; m++){
      if (m==0) cout << "Pt resolution" << endl;
      if (m==1) cout << "Eta resolution" << endl;
      if (m==2) cout << "Phi resolution" << endl;
      for (Int_t molt=0; molt < nummolt+1; molt++){
    //if (molt==0) continue;
	if (ResoHisto=="3D")	cout << " Molt class: " << SmoltLegend[molt] << ", Mean: "<<      fHistResolution_1D[molt][m][t]->GetMean() << " +- " << fHistResolution_1D[molt][m][t]->GetMeanError()/sqrt(fHistResolution_1D[molt][m][t]->GetEntries()) << " RMS: "<< fHistResolution_1D[molt][m][t]->GetRMS()<< endl ;
	else cout << " p_{T}^{Trig, min}: " << ptjminBisFixed[molt] << ", Mean: "<<      fHistResolution_1D[molt][m][t]->GetMean()<< "+- " << fHistResolution_1D[molt][m][t]->GetMeanError()/sqrt(fHistResolution_1D[molt][m][t]->GetEntries()) << " RMS: "<< fHistResolution_1D[molt][m][t]->GetRMS()<< endl ;
      }
    }
  }
  cout << "\n\n" << endl;
  cout << "bin width of EffTriggerPt: " <<     fHistTriggerEfficiencyPt[1]->GetXaxis()->GetBinWidth(1)<< " GeV/c" << endl;
  cout << "bin width  of EffTriggerEta: " <<     fHistTriggerEfficiencyEta[1]->GetXaxis()->GetBinWidth(1) << endl;
  cout << "bin width  of EffTriggerPhi: " <<     fHistTriggerEfficiencyPhi[1]->GetXaxis()->GetBinWidth(1)<< endl;
  cout << "bin width  of EffV0Pt: " <<     fHistV0EfficiencyPt[1]->GetXaxis()->GetBinWidth(1)<< " GeV/c " <<endl;
  cout << "bin width  of EffV0Eta: " <<     fHistV0EfficiencyEta[1]->GetXaxis()->GetBinWidth(1)<< endl;
  cout << "bin width  of EffV0Phi: " <<     fHistV0EfficiencyPhi[1]->GetXaxis()->GetBinWidth(1)<< endl;

  cout << "\n\n" << endl;
  cout << "bin width in #sigmas resolution of EffTriggerPt: " <<     fHistTriggerEfficiencyPt[1]->GetXaxis()->GetBinWidth(1)/fHistResolution_1D[1][0][0]->GetRMS()<< endl;
  cout << "bin width in #sigmas resolution of EffTriggerEta: " <<     fHistTriggerEfficiencyEta[1]->GetXaxis()->GetBinWidth(1)/fHistResolution_1D[1][1][0]->GetRMS()<< endl;
  cout << "bin width in #sigmas resolution of EffTriggerPhi: " <<     fHistTriggerEfficiencyPhi[1]->GetXaxis()->GetBinWidth(1)/fHistResolution_1D[1][2][0]->GetRMS()<< endl;
  cout << "bin width in #sigmas resolution of EffV0Pt: " <<     fHistV0EfficiencyPt[1]->GetXaxis()->GetBinWidth(1)/fHistResolution_1D[1][0][1]->GetRMS()<< endl;
  cout << "bin width in #sigmas resolution of EffV0Eta: " <<     fHistV0EfficiencyEta[1]->GetXaxis()->GetBinWidth(1)/fHistResolution_1D[1][1][1]->GetRMS()<< endl;
  cout << "bin width in #sigmas resolution of EffV0Phi: " <<     fHistV0EfficiencyPhi[1]->GetXaxis()->GetBinWidth(1)/fHistResolution_1D[1][2][1]->GetRMS()<< endl;

  cout << endl;  
  for (Int_t molt=0; molt < nummolt+1; molt++){
    //if (molt==0) continue;
    cout << "molt " << molt << " numero di particelle di trigger selezionate: " << SelEntries[molt]  <<  endl;
    cout << "molt " << molt << " numero di particelle di trigger selezionate con pT > 3 GeV/c: " << SelEntriespT3[molt]  <<  endl;
    cout << "molt " << molt << " numero di particelle di trigger selezionate sul totale di particelle di trigger con pT>3 GeV/c: " << SelEntries[molt]/SelEntriespT3[molt] << endl;
  }

  cout << "******************************************************************"<< endl;
  cout << "partendo dai file "  << PathInBis << "e " << PathInSel << "  ho creato: "<< endl;
  cout << "il file " << PathOut2 << endl;

  cout << "\n\n Run it for sysV0=0,..,6 if !ishhCorr, else for sysV0=0,1,2" << endl;
  cout << "Run it for sysTrigger=0,1,2 (although results won't be used)" << endl;
}
