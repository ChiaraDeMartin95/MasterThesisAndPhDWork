#include <TLine.h>
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

void AngularCorrelationPlot(Bool_t ishhCorr=0, Float_t PtTrigChosen=0.15,Bool_t SkipAssoc=1,Int_t israp=1, Int_t sysV0=0, Int_t sysTrigger=0,Int_t sys=0,Int_t type=0, Int_t PtIntervalShown=1,   TString year0 = "2016",TString year="2016k_hK0s_30runs_150MeV", TString yearMC="2018f1_extra_hK0s_30runs_150MeV",  TString Path1 ="", TString Path2 ="", TString Dir ="FinalOutput/", Bool_t isEnlargedDeltaEta=0){ 

  if (PtIntervalShown>6){
    cout << " Thera are not so many Pt assoc intervals " << endl;
    return;
  }
  if (ishhCorr==1){ year = "2016k"; Path1="_New"; Path2="_New"; yearMC = "2018f1_extra";}

  if (!ishhCorr && type==0 && (((year!="2016k_onlyTriggerWithHighestPt"&& year!="2018f1_extra_onlyTriggerWithHighestPt") || yearMC!="2018f1_extra_onlyTriggerWithHighestPt"))) {
    //    cout << "for hV0 correlation you should use year = 2016k_onlyTriggerWithHighestPt together with the MC yearMC = 2018f1_extra_onlyTriggerWithHighestPt" << endl;
    // return; 
  }

  Float_t PtTrigMin=0;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  //lista degli effetti  sistematici studiati in questa macro
  if (sys==3 || sys>5) return;
  if (sysV0>6)return;
  if (ishhCorr && sys<4 && sys!=0) return; //se faccio correlazione hh non posso studiare sistematico associato a scelta regione Sidebands e peak nella distribuzione di massa invariante
  if (sysV0>2 && ishhCorr) return;
  if (sysTrigger!=0) return;
  if (sysV0!=0 && sys!=0) return;

  Double_t LimSupY[2] = {0.0012, 0.05}; //0.02 for K0s 
  TString hhCorr[2]= {"", "_hhCorr"};
  Int_t sysang=0;
  if (sys==1) sysang=1;
  if (sys==2) sysang=2;

  TLegend *legendPt = new TLegend(0.7, 0.7, 0.9, 0.9);
  legendPt->SetHeader("p_T^{Assoc} intervals");

  Dir+="DATA"+year0;
  TString file[2];
  file[0] = year+Path1;
  file[1] = yearMC + "_MCEff" +Path2;
 
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=8;
  const Int_t numPtTrig=10;
  const Int_t numtipo=10;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
 
  Int_t ColorPt[numPtV0]= {401,801,628,909,881,860,868,842};

  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  if (type>=0)SPtV0[1]={"0.5-1"};
  if (ishhCorr){
    SPtV0[0]={"0-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8};
  if (type>=0) NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0"};
  if (type>=0) SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }

  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  TLegend * legendDataMC = new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t ColorWidth[2]={868, 418};
  TString MCOrNot[2]={" Data" , " Pythia8" };
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};

  TString nameSE[nummolt+1][numPtV0];
  TString nameME[nummolt+1][numPtV0];
  TString nameMEEtaProj[nummolt+1][numPtV0];
  TString nameMEPhiProj[nummolt+1][numPtV0];
  TString namePhiProjJet[nummolt+1][numPtV0];
  TString namePhiProjBulk[nummolt+1][numPtV0];
  TString namePhiProjJetBulk[nummolt+1][numPtV0];
  TString nameEtaProj[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[nummolt+1][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[nummolt+1][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulkRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProjPHalf[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProjNHalf[nummolt+1][numPtV0][numPtTrig];
  TString PathIn[numPtTrig];

  TF1* Baseline[nummolt+1][numPtV0][numPtTrig];
  TF1* Gaussian[nummolt+1][numPtV0][numPtTrig];
  TF1* GaussianAS[nummolt+1][numPtV0][numPtTrig];
  TF1* GaussianEta[nummolt+1][numPtV0][numPtTrig];

  //  TLine *tlineEtaSx=new TLine(-1.04, -TMath::Pi()/2, -1.04,  3*TMath::Pi()/2);
  //  TLine *tlineEtaDx=new TLine(1.04, -TMath::Pi()/2, 1.04,  3*TMath::Pi()/2);
  TLine *tlineEtaSx=new TLine(-0.78, -TMath::Pi()/2, -0.78,  3*TMath::Pi()/2);
  TLine *tlineEtaDx=new TLine(0.78, -TMath::Pi()/2, 0.78,  3*TMath::Pi()/2);
  TLine *tlineEtaInclSx=new TLine(-1.14, -TMath::Pi()/2, -1.14,  3*TMath::Pi()/2);
  TLine *tlineEtaInclDx=new TLine(1.14, -TMath::Pi()/2, 1.14,  3*TMath::Pi()/2);
  TLine *tlinePhiSx=new TLine(-1.5, -1.06, 1.5,  -1.06); //was 1.32
  TLine *tlinePhiDx=new TLine(-1.5, +1.06, 1.5,  +1.06);
  TLine *tlinePhiBase=new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2,0);

  TCanvas *  canvasWidthGaussian[2];
  TCanvas *  canvasWidthGaussianEta[2];
  TCanvas *  canvasWidthGaussianAS[2];
  for (Int_t i =0; i<2; i++){  
    canvasWidthGaussian[i] = new TCanvas (Form("canvasWidthGaussian%i",i), Form("canvasWidthGaussian%i",i), 800, 500);
    canvasWidthGaussian[i]->Divide(3,2);
    canvasWidthGaussianEta[i] = new TCanvas (Form("canvasWidthGaussianEta%i",i), Form("canvasWidthGaussianEta%i",i), 800, 500);
    canvasWidthGaussianEta[i]->Divide(3,2);
    canvasWidthGaussianAS[i] = new TCanvas (Form("canvasWidthGaussianAS%i",i), Form("canvasWidthGaussianAS%i",i), 800, 500);
    canvasWidthGaussianAS[i]->Divide(3,2);
  }
  TCanvas *canvasSPectrum;
  TCanvas *canvasPlot[nummolt+1];
  TCanvas *canvasPlotME[nummolt+1];
  TCanvas *canvasPlotMEEtaProj[2];
  TCanvas *canvasPlotMEPhiProj[2];
  TCanvas *canvasPlotProj[nummolt+1];
  TCanvas *canvasPlotProjEta[nummolt+1][3];
  TH1F *HistoWidthGaussian[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianAS[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianEta[nummolt+1][2][numPtTrig];
  TH1F *RatioHistoWidthGaussian[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianAS[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianEta[nummolt+1][numPtTrig];
  TString PhiRegion[3]={"All", "DeltaPhi>Pi/2", "DeltaPhi < Pi/2"};

  //*********************************************************************                                                                                                  
  //**************calcolo numero particelle di trigger*******************                                                                                                
  Int_t NTrigger[nummolt+1][2];
  TFile *fileinbis[2]; //Data and MC
  TString PathInBis[2];
  PathInBis[0] =  "FinalOutput/AnalysisResults" + year +".root"; //+ Path1  + ".root"; change
  PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC  + "_MCEff" + Path1 +".root";
  if (ishhCorr){  
    PathInBis[0] =  "FinalOutput/AnalysisResults" + year +Path1+ "_hhCorr.root";
    PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC  + Path1+"_hhCorr_MCEff.root";
  }

  for (Int_t IntisMC=0; IntisMC<=0; IntisMC++){ //change
    fileinbis[IntisMC]=new TFile(PathInBis[IntisMC],"");
    cout <<"file from task: " << fileinbis[IntisMC] << endl;
    if (!fileinbis[IntisMC]) return;
    TString dirinputtype[numtipo] = {"", "Lambda", "Lambda", "Xi", "Lambda","Xi", "Omega", "Omega", "Xi", "Omega"};
    TDirectoryFile *dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTask"+dirinputtype[type]);
    TList *list = (TList*)dir->Get("MyOutputContainer");
    TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
    if (!fHistTriggervsMult ) return;
    cout << "ok up to here " << endl;
    TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigChosen+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(30-0.00001) );

    //*********************************************************************                                                                                                
 
    for(Int_t m=0; m<nummolt+1; m++){
      NTrigger[m][IntisMC]=0;
      if(m<nummolt){
	for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	  NTrigger[m][IntisMC]+= fHistTriggervsMult_MultProj->GetBinContent(j);
	}
      }
      else {
	for(Int_t j=1; j<=fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	  NTrigger[m][IntisMC]+= fHistTriggervsMult_MultProj->GetBinContent(j);
	}
      }
      cout << "n trigger in mult range (all triggers)    " << m << "  " <<  " is MC " << IntisMC << "  " << NTrigger[m][IntisMC] <<   endl;
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////              

  gStyle->SetOptStat(0);
 
  TFile *filein[numPtTrig];
  TString nomefileoutput= "AngularCorrelationPlot"+year;
  if(type>=0){
    nomefileoutput +="_"+tipo[type];
    nomefileoutput +=Srap[israp];
    nomefileoutput +=SSkipAssoc[SkipAssoc];
  }

  nomefileoutput += hhCorr[ishhCorr] + Form("_PtTrigMin%.1f_Output.root", PtTrigChosen);
  if (!ishhCorr && isEnlargedDeltaEta)  nomefileoutput= "AngularCorrelationPlot" + hhCorr[ishhCorr] + Form("_PtTrigMin%i_DeltaEtaPhiEnlarged_Output.root", PtTrigChosen);
   
  TFile *fileout = new TFile(nomefileoutput, "RECREATE");

  //disegno spettri in pTassoc per le varie molteplicità
  for(Int_t m=nummolt; m>=0; m--){


  }
  //*****************************************************
  Int_t counter=0;
  TCanvas*   canvasSummedPt=new TCanvas ("canvasSummedPt", "canvasSummedPt", 800, 500);
  canvasSummedPt->Divide(nummolt+1, 2);
	for (Int_t IntisMC=0; IntisMC<=0; IntisMC++){
	  canvasPlotMEEtaProj[IntisMC]=new TCanvas(Form("canvasPlotMEEtaProj_%i", IntisMC), Form("canvasPlotMEEtaProj_%i", IntisMC), 1300, 800);
    canvasPlotMEEtaProj[IntisMC]->Divide(6,2);
    canvasPlotMEPhiProj[IntisMC]=new TCanvas(Form("canvasPlotMEPhiProj_%i", IntisMC), Form("canvasPlotMEPhiProj_%i", IntisMC), 1300, 800);
    canvasPlotMEPhiProj[IntisMC]->Divide(6,2);
	}
  for(Int_t m=nummolt; m>=0; m--){
    cout << "\n\n m " << m << endl;
    canvasPlot[m]=new TCanvas(Form("canvasPlot_m%i",m), "canvasPlot"+Smolt[m], 1300, 800);
    canvasPlot[m]->Divide(4,2);
    canvasPlotME[m]=new TCanvas(Form("canvasPlotME_m%i",m), "canvasPlotME"+Smolt[m], 1300, 800);
    canvasPlotME[m]->Divide(4,2);
    canvasPlotProj[m]=new TCanvas(Form("canvasPlotProj_m%i",m), "canvasPlotProj"+Smolt[m], 1300, 800);
    canvasPlotProj[m]->Divide(4,2);
    for (Int_t t=0; t<3; t++){
      canvasPlotProjEta[m][t]=new TCanvas(Form("canvasPlotProjEta_m%i_%i",m,t), "canvasPlotProjEta"+Smolt[m]+"_"+PhiRegion[t], 1300, 800);
    canvasPlotProjEta[m][t]->Divide(4,2);
    }

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){

      PtTrigMin=PtTrig+PtTrigChosen;
      if (PtTrigMin!=PtTrigChosen) continue;
      counter++;
      cout << "PtTrig " << PtTrig << endl;

      HistoWidthGaussian[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussian[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussianAS[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianASm%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussianAS[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianASm%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussianEta[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianEtam%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussianEta[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianEtam%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);

    }


    for(Int_t v=1; v<numPtV0; v++){  
      //    if (v!=PtIntervalShown) continue;
      //	if (v>4) continue;
      nameSE[m][v]="ME_";
      nameSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_AC";
      nameME[m][v]="ME_m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_norm";
      nameMEEtaProj[m][v] =nameME[m][v]+"_EtaProj";
      nameMEPhiProj[m][v] =nameME[m][v]+"_PhiProj"; 
      namePhiProjJet[m][v]= nameSE[m][v] + "_phi_V0Sub_BulkSub_EffCorr";
      namePhiProjBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_Bulk_EffCorr";
      namePhiProjJetBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_JetBulkEffCorr";
      //      cout << nameSE[m][v]<< endl;

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	PtTrigMin=PtTrig+PtTrigChosen;
	if (PtTrigMin!=PtTrigChosen) continue;
	cout << "PtTrig " << PtTrig << endl;
	cout << "PtTrigMin " << PtTrigMin << endl;
	cout << "PtTrigChosen " << PtTrigChosen << endl;
	nameEtaProj[m][v][PtTrig]= nameSE[m][v] + Form("_EtaProj_PtTrig%i", PtTrig);
	//	if (PtTrigMin==4 || PtTrigMin==5 || PtTrigMin>10)continue;
	if (PtTrigMin>7) continue;
	for (Int_t IntisMC=0; IntisMC<=0; IntisMC++){
	
	  PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC];
	  if(type>=0){
	      PathIn[PtTrig] +="_"+tipo[type];
	      PathIn[PtTrig] +=Srap[israp];
	      PathIn[PtTrig] +=SSkipAssoc[SkipAssoc];
	    }
	   PathIn[PtTrig]+= hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_Output.root";
	  
	  if (!ishhCorr && isEnlargedDeltaEta)	  PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC] + hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_DeltaEtaPhiEnlarged_Output.root";
	  cout <<"path in : " <<  PathIn[PtTrig] << endl;
	  filein[PtTrig]= new TFile(PathIn[PtTrig], "");
	
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSE[m][v]);
	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameME[m][v]);
	  if (v==1)  hDeltaEtaDeltaPhi_MEbins[m][5][PtTrig]= (TH2F*)filein[PtTrig]->Get("ME_m"+ Smolt[m]+"_v"+SPtV0[5]+Ssideband[0]+"_norm");
	  if (!hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]){cout << "no SE 2D histo " << endl;  return;}
	  if (!hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]){cout << "no ME 2D histo " << endl;  return;}
	  if (!hDeltaEtaDeltaPhi_MEbins[m][5][PtTrig]){cout << "no ME 2D histo for v==5" << endl;  return;}

	  //projection along eta of the Mixed Event distribution
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionX(nameMEEtaProj[m][v],0,-1, "E");
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());

	  hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Clone(nameMEEtaProj[m][v]+"_Ratio");
	  //	  if (v==1)	{
	    hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][5][PtTrig]->ProjectionX(nameMEEtaProj[m][5]+ "_Master",0,-1, "E");
	    hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->Scale(1./ hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->GetNbinsY());
	    //	  }
	  hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]);

	  for (Int_t i=1; i<=   hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->SetBinError(i, sqrt(TMath::Abs(pow( hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetBinError(i),2)-pow( hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->GetBinError(i),2)))/ hDeltaEtaDeltaPhi_MEbins_EtaProjMaster[m][PtTrig]->GetBinContent(i));
	  }


	  cout << "ho preso isto" << endl;

	  //projection along phi of the Mixed Event distribution
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->ProjectionY(nameMEPhiProj[m][v],0,-1, "E");

	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->SetLineColor(ColorPt[v]);
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->SetMarkerColor(ColorPt[v]);

	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Clone(nameMEPhiProj[m][v]+"_Ratio");
	  if (v==1)	  hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]= 	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]= (TH1F*)	  hDeltaEtaDeltaPhi_MEbins[m][5][PtTrig]->ProjectionY(nameMEPhiProj[m][5]+"_Master",0,-1, "E");
	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->Divide( hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]);

	  for (Int_t i=1; i<=   hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->SetBinError(i,sqrt( TMath::Abs( pow(hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->GetBinError(i),2)- pow(hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->GetBinError(i),2)))/hDeltaEtaDeltaPhi_MEbins_PhiProjMaster[m][PtTrig]->GetBinContent(i));
	  }
	  cout << "ho preso isto" << endl;

	  //other projections
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]+"_RelError");
	  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]+"_RelError");
	  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]+"_RelError");
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig], 0,-1, "E");
	  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig],hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(+TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(+3./2*TMath::Pi()- 0.001)  , "E");
	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig],hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(-TMath::Pi()/2+ 0.001) ,hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->GetXaxis()->FindBin(+1./2*TMath::Pi()- 0.001)  , "E");

	  cout << "ho preso isto" << endl;

	  cout << "this rebin is performed for visualization purposes only " << endl;
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Rebin(2);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Rebin(2);

	  //canvasPlot[m]->cd(PtTrig+1);
	  if (IntisMC==0)	canvasPlot[m]->cd(v);
	  else	canvasPlot[m]->cd(v+7);
	  cout << "scelto cd " << endl;
	  //	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetLogz();
	  //	gPad->SetLogz();
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->SetTitle("Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin) +MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->SetTitle("Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin)+ MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetTitle("Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin)+ MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->Draw("colz");
	  tlinePhiDx->Draw();
	  tlinePhiSx->Draw();
	  tlineEtaDx->Draw();
	  tlineEtaSx->Draw();
	  tlineEtaInclDx->Draw();
	  tlineEtaInclSx->Draw();

	  if (IntisMC==0)	canvasPlotME[m]->cd(v);
	  else	canvasPlotME[m]->cd(v+7);
	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_MEbins[m][v][PtTrig]->Draw("colz");

	  cout << "m " << m << " v " << v << " "  << counter << endl;
	  if (counter==1)	  legendPt->AddEntry(	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig], SPtV0[v], "pl");
	  canvasPlotMEEtaProj[IntisMC]->cd(m+1);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);
	  hDeltaEtaDeltaPhi_MEbins_EtaProj[m][v][PtTrig]->Draw("same");
	  if (v==numPtV0-1) legendPt->Draw("same");
	  canvasPlotMEEtaProj[IntisMC]->cd(m+1+6);
	  hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	  hDeltaEtaDeltaPhi_MEbins_EtaProjRatio[m][v][PtTrig]->Draw("same");
	  if (v==numPtV0-1) legendPt->Draw("same");

	  canvasPlotMEPhiProj[IntisMC]->cd(m+1);
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0,20);
	  hDeltaEtaDeltaPhi_MEbins_PhiProj[m][v][PtTrig]->Draw("same");
	  if (v==numPtV0-1) legendPt->Draw("same");
	  canvasPlotMEPhiProj[IntisMC]->cd(m+1+6);
	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->GetYaxis()->SetRangeUser(0.8, 1.2);
	  hDeltaEtaDeltaPhi_MEbins_PhiProjRatio[m][v][PtTrig]->Draw("same");
	  if (v==numPtV0-1) legendPt->Draw("same");

	  //	canvasPlot[m]->cd(PtTrig+1+7);
	  cout << "setting to zero last bin " << endl;
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetLineColor(868);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetMarkerColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetMarkerColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetMarkerColor(868);

	  if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
	  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->SetLineColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->SetLineColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->SetLineColor(868);
	  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->SetMarkerColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->SetMarkerColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->SetMarkerColor(868);
	  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.07);
	  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->GetYaxis()->SetLabelSize(0.07);
	  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->GetYaxis()->SetLabelSize(0.07);
	  }
	  //division by number of trigger particles
	  cout << "division by trigger particles " << endl;
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Scale(1./NTrigger[m][IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Scale(1./NTrigger[m][IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Scale(1./NTrigger[m][IntisMC]);

	  //sum of different ptV0 bins
	  if (v==1){
	    cout << "Adding histos to make ptv0 sum" << endl;
	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    =	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Clone(Form("PhiProjJet_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   =	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Clone(Form("PhiProjBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]=	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Clone(Form("PhiProjJetBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
    
	  }

	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->SetTitle(Form("PhiProjJet_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->SetTitle(Form("PhiProjBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));	  
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->SetTitle(Form("PhiProjJetBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));

	  if (v!=1){
	    cout << "Adding histos to make ptv0 sum" << endl;
	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]);
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]);
	  }


	  if (IntisMC==0)	canvasPlotProj[m]->cd(v);
	  else canvasPlotProj[m]->cd(v+7);
	  cout << "cd chosen " << endl;
	  gPad->SetLeftMargin(0.2);

	  Gaussian[m][v][PtTrig]= new TF1 ( Form("GaussianFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  -0.8, 0.8);
	  GaussianAS[m][v][PtTrig]= new TF1 ( Form("GaussianASFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  TMath::Pi()-1.2,   TMath::Pi()+1.2);
	  GaussianEta[m][v][PtTrig]= new TF1 ( Form("GaussianEtaFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  -0.8,   0.8);
	  Baseline[m][v][PtTrig]= new TF1 ( Form("Baseline_m%i_v%i_PtTrig%i", m,v,PtTrig),"pol0",  1.5,  4);
	  Baseline[m][v][PtTrig]->SetLineColor(881);
	  GaussianAS[m][v][PtTrig]->SetLineColor(419);
	  GaussianEta[m][v][PtTrig]->SetLineColor(860);

	  cout << "m " << m << " v " << v << " PtTrig " << PtTrig << "  " << MCOrNot[IntisMC] << endl;
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Fit(Gaussian[m][v][PtTrig], "R");
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Fit(Baseline[m][v][PtTrig], "R0");
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Fit(GaussianAS[m][v][PtTrig], "R0");
	  if (v>1 && v <4) GaussianEta[m][v][PtTrig]->SetRange(-0.5, 0.5);
	  else if (v>=4)	  GaussianEta[m][v][PtTrig]->SetRange(-0.3, 0.3);
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Fit(GaussianEta[m][v][PtTrig], "R");

	  cout << "\nbaseline value" <<  Baseline[m][v][PtTrig]->GetParameter(0)<< " +- " <<  Baseline[m][v][PtTrig]->GetParError(0)<< endl;
	  cout << "\njet fit: chi square " << Gaussian[m][v][PtTrig]->GetChisquare() << " NDF " << Gaussian[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << Gaussian[m][v][PtTrig]->GetChisquare()/Gaussian[m][v][PtTrig]->GetNDF() << endl;
	  cout << "\nAwaySide fit: chi square " << GaussianAS[m][v][PtTrig]->GetChisquare() << " NDF " << GaussianAS[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << GaussianAS[m][v][PtTrig]->GetChisquare()/GaussianAS[m][v][PtTrig]->GetNDF() << endl;

	  //	  if (	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetMaximum() > 	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetMaximum()) 	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Draw("");


	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetRangeUser(-0.0001, LimSupY[ishhCorr]); //-0.004 for hh and hK0s
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(-0.0001, LimSupY[ishhCorr]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetRangeUser(-0.0001, LimSupY[ishhCorr]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetLabelSize(0.04);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.04);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetLabelSize(0.04);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetTitle(" "); //"N/N_{Trigg} per #Delta#eta"); 
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetTitle(" "); //"N/N_{Trigg} per #Delta#eta"); 
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetTitle(" "); //"N/N_{Trigg} per #Delta#eta"); 
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetTitleOffset(1.8); 
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetTitleOffset(1.8); 
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetTitleOffset(1.8); 

	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Draw("");
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Draw("same");

	  tlinePhiBase->Draw("");

	  if (IntisMC==0){
	    //	  canvasPlotProj[m]->cd(v+7);
	  if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
	    /*
	  hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]->Draw("");
	  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]   ->Draw("same");
	    */
	    cout << "here an average value of the relative error of the jet deltaPhi projection in the range ~[-1; 1]" << endl;
	    Int_t Counter=0;
	    Float_t Average=0;
	    for (Int_t b=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()    ->FindBin(-1); b<=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()->FindBin(1); b++){
	      Counter++;
	      Average+=  hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->GetBinError(b);
	    }
	    Average = Average/Counter;
	    cout << "average " << Average << endl;
	  }
	  }

	  for (Int_t t=0; t<3; t++){
	  if (IntisMC==0)	canvasPlotProjEta[m][t]->cd(v);
	  else canvasPlotProjEta[m][t]->cd(v+7);
	  gPad->SetLeftMargin(0.15);
	  if (t==0){
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Draw("");
	  }
	  else if (t==1){
	  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
	  hDeltaEtaDeltaPhi_EtaProjPHalf[m][v][PtTrig]->Draw("");
	  }
	  else if (t==2){
	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.5* hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
	  hDeltaEtaDeltaPhi_EtaProjNHalf[m][v][PtTrig]->Draw("");
	  }

	  }
	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinContent(v+1, Gaussian[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinError(v+1, Gaussian[m][v][PtTrig]->GetParError(2));
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianAS[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianAS[m][v][PtTrig]->GetParError(2));
	  HistoWidthGaussianEta[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianEta[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussianEta[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianEta[m][v][PtTrig]->GetParError(2));

	} //end of IntisMC loop
      } //end of PtTrigMin loop
    } //end of v loop


    //I draw delta-phi projection done summing over pT assoc bins*****
    for (Int_t IntisMC=0; IntisMC<=0; IntisMC++){
      if (!IntisMC)      canvasSummedPt->cd(1+m);
      else       canvasSummedPt->cd(1+m+nummolt+1);
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
	for(Int_t v=1; v<numPtV0; v++){  
	  cout <<" m " << m <<  " v " << v << " "  <<  	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetMaximum() << endl;
	}
	cout <<" m " << m <<  " " <<	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetMaximum() << endl;
	hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Scale(2.28);
	hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->Scale(2.28);
	hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->Scale(2.28);
	cout <<" m " << m <<  " " <<	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetMaximum() << endl;
	hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetYaxis()->SetRangeUser(-0.004, 0.007);//0.5 for hh
	hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->GetYaxis()->SetRangeUser(-0.004, 0.007);
	hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(-0.004, 0.007);

	hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetYaxis()->SetLabelSize(0.05);
	hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->GetYaxis()->SetLabelSize(0.05);
	hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->GetYaxis()->SetLabelSize(0.05);

	hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Draw("");
	hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->Draw("same");
	hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->Draw("same");
      }
    }

    //**************************************************************
    for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;

	HistoWidthGaussian[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0,0.7);
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);
	if (m==0)	legendDataMC->AddEntry(	HistoWidthGaussian[m][IntisMC][PtTrig],MCOrNot[IntisMC], "pl"); 

	HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0,5);
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);

	canvasWidthGaussian[0]->cd(m+1);
	HistoWidthGaussian[m][IntisMC][PtTrig]->Draw("same");
	if (IntisMC==1)	legendDataMC->Draw("same");

	canvasWidthGaussianAS[0]->cd(m+1);
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->Draw("same");
	if (IntisMC==1)	legendDataMC->Draw("same");

	canvasWidthGaussianEta[0]->cd(m+1);
	HistoWidthGaussianEta[m][IntisMC][PtTrig]->Draw("same");
	if (IntisMC==1)	legendDataMC->Draw("same");
      } 
    } //end second intisMC loop

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      PtTrigMin=PtTrig+PtTrigChosen;
      if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
      RatioHistoWidthGaussian[m][PtTrig]=(TH1F*)	HistoWidthGaussian[m][0][PtTrig]->Clone(Form("RatioDataPythiaWidthm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussian[m][PtTrig]->SetTitle(Form("RatioDataPythiaWidthm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussian[m][PtTrig]->Sumw2();
      RatioHistoWidthGaussian[m][PtTrig]->Divide(HistoWidthGaussian[m][1][PtTrig]);
      RatioHistoWidthGaussian[m][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);

      RatioHistoWidthGaussianAS[m][PtTrig]=(TH1F*)	HistoWidthGaussianAS[m][0][PtTrig]->Clone(Form("RatioDataPythiaWidthASm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussianAS[m][PtTrig]->SetTitle(Form("RatioDataPythiaWidthASm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussianAS[m][PtTrig]->Sumw2();
      RatioHistoWidthGaussianAS[m][PtTrig]->Divide(HistoWidthGaussianAS[m][1][PtTrig]);
      RatioHistoWidthGaussianAS[m][PtTrig]->GetYaxis()->SetRangeUser(0,2);

      canvasWidthGaussian[1]->cd(m+1);
      RatioHistoWidthGaussian[m][PtTrig]->Draw("same");

      canvasWidthGaussianAS[1]->cd(m+1);
      RatioHistoWidthGaussianAS[m][PtTrig]->Draw("same");	  
    }

    //    if (m==0 || m>3){
    fileout->WriteTObject(canvasPlot[m]); 
    fileout->WriteTObject(canvasPlotME[m]); 
    fileout->WriteTObject(canvasPlotProj[m]); 
    for (Int_t t=0; t<3; t++){
    fileout->WriteTObject(canvasPlotProjEta[m][t]); 
    }
    //    }

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
      
	fileout->WriteTObject(RatioHistoWidthGaussian[m][PtTrig]);
	fileout->WriteTObject(HistoWidthGaussian[m][0][PtTrig]);
	fileout->WriteTObject(HistoWidthGaussian[m][1][PtTrig]);
	fileout->WriteTObject(RatioHistoWidthGaussianAS[m][PtTrig]);
	fileout->WriteTObject(HistoWidthGaussianAS[m][0][PtTrig]);
	fileout->WriteTObject(HistoWidthGaussianAS[m][1][PtTrig]);
	fileout->WriteTObject(HistoWidthGaussianEta[m][0][PtTrig]);
	fileout->WriteTObject(HistoWidthGaussianEta[m][1][PtTrig]);
      
    }
  } //end of mult loop
    fileout->WriteTObject(canvasPlotMEEtaProj[0]); 
    fileout->WriteTObject(canvasPlotMEPhiProj[0]); 
  fileout->WriteTObject(canvasSummedPt); 
  fileout->WriteTObject(canvasWidthGaussian[0]);
  fileout->WriteTObject(canvasWidthGaussianAS[0]);
  fileout->WriteTObject(canvasWidthGaussian[1]);
  fileout->WriteTObject(canvasWidthGaussianAS[1]);
  fileout->WriteTObject(canvasWidthGaussianEta[0]);
  fileout->WriteTObject(canvasWidthGaussianEta[1]);
  fileout->Close();

  cout << "baseline fits " << endl;
  for(Int_t m=nummolt; m>=0; m--){
    cout << "\n\n" << endl;
    for(Int_t v=1; v<numPtV0; v++){  
    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
      cout << " m " << m << " v " << v << " PtTrig " << PtTrig  << endl;
      cout << "baseline value" <<  Baseline[m][v][PtTrig]->GetParameter(0)<< " +- " <<  Baseline[m][v][PtTrig]->GetParError(0)<< " number of sigmas from zero " <<  TMath::Abs(Baseline[m][v][PtTrig]->GetParameter(0)/ Baseline[m][v][PtTrig]->GetParError(0)) << endl;
    }
    }
  }

  for(Int_t m=nummolt; m>=0; m--){
    cout << "\n\n" << endl;
    for(Int_t v=1; v<numPtV0; v++){  
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+PtTrigChosen)!=PtTrigChosen) continue;
 
	if (hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig] && hDeltaEtaDeltaPhi_PhiProjJetBulkRelError[m][v][PtTrig]&& hDeltaEtaDeltaPhi_PhiProjBulkRelError[m][v][PtTrig]){
	  cout << " m " << m << " v " << v << " PtTrig " << PtTrig  << endl;
	  //	  cout << "here an average value of the relative error of the jet deltaPhi projection in the range ~[-1; 1]" << endl;
	  Int_t Counter=0;
	  Float_t Average=0;
	  for (Int_t b=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()    ->FindBin(-1); b<=hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]->GetXaxis()->FindBin(1); b++){
	    Counter++;
	    Average+= TMath::Abs( hDeltaEtaDeltaPhi_PhiProjJetRelError[m][v][PtTrig]    ->GetBinContent(b));
	  }
	  Average = Average/Counter;
	  cout << "average " << Average << endl;
	}
      }
    }
  }
  cout << "\n\n ho creato il file " << nomefileoutput << endl;
}

