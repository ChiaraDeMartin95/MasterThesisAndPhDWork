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

void BarlowSys( Int_t avoidthissyst=3,Int_t type=0, Bool_t isMC=0, Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year="2016k",TString year0="2016", TString yearMC="2018d8",  TString Path1 ="", TString Path2 ="",TString Path3="_15runs_6thtry", Int_t sysTrigger=0, Float_t numSigmaCorr=2){

  if (sysTrigger!=0){
    cout << "sysTrigger must be zero" << endl;
    return;
  }
 
  TFile *fileinbis;
  TFile *filein;
  TString PathIn1;
  TString file = year+Path1;
  TString PathInBis =  "FinalOutput/AnalysisResults" + file  + ".root";
  TString PathInTris =  "FinalOutput/DATA" + year0 + "/histo/AngularCorrelation" +file  + "_SysT0_SysV00_Sys0_Output.root";
  fileinbis=new TFile(PathInBis,"");

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numSyst  = 13;
  const Int_t numsystPhi  = 5;

  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
   TString SSysT[numSysTrigger]={"DCAz < 1","DCAz < 2","DCAz < 0.5"};
  TString SSysV0[numSysV0]={"default", "cosTP> 0.997", "ctau <3 ", "YK0s < 0.5", "Lrejection 10 MeV","no Lrejection 5 MeV", "V0dca< 0.3"}; //all except no Lrejecxtion 5 MeV are done with the default cut on (i.e. Lrejection 5 MeV)
  TString SSyst[numSyst]={"default", "cosTP> 0.997", "ctau <3 ", "YK0s < 0.5", "Lrejection 10 MeV","no Lrejection 5 MeV", "V0dca< 0.3", "Sideband 5sigma", "central 4sigma", "", "jet0.5", "bulk",  "BC 1.2"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-1.5","1.5-2", "2-2.5", "2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,1.5,2,2.5,3,4,8};
  TString SErrorSpectrum[3]={"stat","sist_uncorr","sist_corr"};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};

  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};

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
  auto legend3= new TLegend(0.6, 0.6, 1, 1);
  legend3->SetHeader("Selezioni applicate 3");     

  Int_t sys=0;
  TH1D* fHistSpectrum_master[nummolt+1];
  TH1D* fHistSpectrum_master_bis[nummolt+1];
  TH1D* fHistSpectrum_masterSystCorr[nummolt+1];
  TH1D* fHistSpectrum_masterSystUnCorr[nummolt+1];
  TH1D* fHistSpectrum_masterTotalUnCorr[nummolt+1];
  TH1D* fHistSpectrum[nummolt+1][numSyst];  
  TH1D* fHistSpectrum_Corr[nummolt+1][numSyst];
  TH1D* fHistSpectrum_ratio[nummolt+1][numSyst];
  TH1D* fHistSpectrum_Barlow[nummolt+1][numSyst];
  TH1D* fHistSigmaBarlowSpectrum[nummolt+1][numSyst];

  TH1D* fHistPhiDistr[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_solostat[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_master[nummolt+1][numPtV0];
  TH1D* fHistPhiDistr_ratio[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_Barlow[nummolt+1][numPtV0][numsystPhi + numSysV0];
  TH1D* fHistPhiDistr_Corr[nummolt+1][numPtV0][numsystPhi + numSysV0];  
  TH1D* fHistPhiDistr_CorrClone[nummolt+1][numPtV0][numsystPhi + numSysV0];

  TF1*  gauss[nummolt+1][numPtV0][numSyst];
  TF1*  gaussint[nummolt+1][numPtV0][numSyst];

  Double_t MiaNonna[nummolt+1][numPtV0]={0};
  Double_t Sigma[nummolt+1][numPtV0][numSyst]={0};
  Double_t SigmaSpectrum[nummolt+1][numSyst]={0};
  Double_t SigmaSystSpectrumUnCorr[nummolt+1][numSyst]={0};
  Double_t SigmaSystSpectrumCorr[nummolt+1][numSyst]={0};
  Double_t SigmaTotalSpectrumUnCorr[nummolt+1][numSyst]={0};
  Double_t Mean[nummolt+1][numPtV0][numSyst]={0}; 
  Double_t MeanSpectrum[nummolt+1][numSyst]={0};
  Bool_t   Correlated[nummolt+1][numPtV0][numSyst]={0};  
  Bool_t   CorrelatedSpectrum[nummolt+1][numSyst]={0};  
  Bool_t   CorrelatedBis[nummolt+1][numSyst]={0};
  Int_t    CorrelatedBisSys[nummolt+1]={0};
  Bool_t   BarlowPassed[nummolt+1][numPtV0][numSyst]={0};
  Bool_t   BarlowPassedSpectrum[nummolt+1][numSyst]={0};
  Int_t    NumberBarlowPassed[nummolt+1][numPtV0]={0};
  Int_t    NumberBarlowPassedSpectrum[nummolt+1]={0};
  Int_t    NumberCorr[nummolt+1][numPtV0]={0};  Int_t    NumberCorrSpectrum[nummolt+1]={0};
  Double_t SigmaBarlow[nummolt+1][numPtV0][numSyst][50]={0};

  Int_t   NTrigger[nummolt+1]={0}; //total number of trigger particles 
  Int_t   NTriggerV0[nummolt+1]={0}; //total number of trigger particles only in events with V > 0
  Float_t NSpectrum[nummolt+1][numPtV0][numSyst]={0}; //total number of trigger particles
  Float_t NSpectrumFinal[nummolt+1][numPtV0][100]={0}; //total number of trigger particles
  Float_t NSpectrumError[nummolt+1][numPtV0][numSyst]={0}; //total number of trigger particle
  Float_t NSpectrumErrorFinal[nummolt+1][numPtV0][100]={0}; //total number of trigger particle
  Float_t NSpectrumErrorSistUnCorrPhi[nummolt+1][numPtV0][numSyst]={0}; //total number of trigger particle
  Float_t NSpectrumErrorSistUnCorrPhiFinal[nummolt+1][numPtV0][100]={0}; //total number of trigger particles[numSyst]={0}; //total number of trigger particle
  Float_t NSpectrumErrorSoloStat[nummolt+1][numPtV0][numSyst]={0}; //total number of trigger particles
  Float_t NSpectrumErrorSoloStatFinal[nummolt+1][numPtV0][100]={0}; //total number of trigger particles

  TH1D* fHistSigmaSyst[nummolt+1][numPtV0][numSyst];
  TH1D* fHistSigmaSystNSmoothed[nummolt+1][numPtV0][numSyst];
  Double_t SigmaSyst[nummolt+1][numPtV0][50]={0};
  
  TH1F*   fHistV0EfficiencyPtBins[nummolt];
  TH1F*   HistContV0PtBins[nummolt];

  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask");
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtvsMultBefAll");
  TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj");

  TString stringout = Dir+"/DATA"+year0+"/SystematicAnalysis.root";



  // TCanvas *canvasJetFit[nummolt+1][2];
  // for(Int_t m=0; m< nummolt+1; m++){
  //   for(Int_t i=0; i< 2; i++){
  //   canvasJetFit[m][i]= new TCanvas (Form("canvasJetFit_%i_sys%i",m,12 ),Form("canvasJetFit_%i_sys%i",m,12) , 1200, 1000);
  //   canvasJetFit[m][i]->Divide(4,2);
  //   cout << m << endl;
  //   }
  // }

  //********************************************************************* 
  //definitio of integral regions in deltaPhi distribution
  //********************************************************************* 
  Float_t ALowBin[numSyst-numSysV0-numsystPhi+1]={-1}; 
  Float_t AUpBin[numSyst-numSysV0-numsystPhi+1]={1};
  Float_t ALowBinFit[numSyst-numSysV0-numsystPhi+1]={-1};
  Float_t AUpBinFit[numSyst-numSysV0-numsystPhi+1]={1};

  ALowBinFit[0]=	ALowBin[0]={-1};
  ALowBinFit[1]=	ALowBin[1]={-1.4};

  AUpBinFit[0]=	AUpBin[0]={1};
  AUpBinFit[1]=	AUpBin[1]={1.4};

  //********************************************************************* 
  //**************calcolo numero particelle di trigger*******************
  //********************************************************************* 
  TFile *fileAngularCorrMaster = new TFile(PathInTris, "");
  cout << PathInTris << endl;
  for(Int_t m=0; m<nummolt+1; m++){

    for(Int_t syst=0; syst<numSyst; syst++){
      fHistSpectrum[m][syst]=new TH1D ("fHistSpectrum_"+Smolt[m]+Form("_Sys%i",syst),"fHistSpectrum_"+Smolt[m]+Form("_Sys%i",syst), numPtV0, NPtV0) ;      

    }

    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    cout << "n trigger in mult range (all triggers)    " << m << "  " <<  NTrigger[m] <<   endl;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////

  Int_t sysV0=0;    
  if(isMC && isEfficiency) file = year + "_MCEff" + Path2;
  if(isMC && !isEfficiency) file = yearMC + "_MCTruth" + Path3;
 
  for(Int_t m=0; m<nummolt+1; m++){
    //cout << "\n\n****************************************************\n" <<m << endl;
    //    cout << "\n****************************************************\n" << endl;
   
    for(Int_t syst=0; syst<numSysV0+numsystPhi; syst++){
      //      cout << "********************************* syst = " << syst << endl;
      if (syst==9 || syst ==avoidthissyst) continue;
      if (syst< numSysV0){
	sysV0=syst;
	sys=0;
      }
      else {
	sys=syst-numSysV0+1;
	sysV0=0;
      }
      if (isMC && !isEfficiency) PathIn1=Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file+  Form("_Sys%i", sys)+"_Output.root";
      else PathIn1 = Dir+"/DATA"+year0+"/histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i_Sys%i", sysTrigger, sysV0, sys)+"_Output.root";
      cout << "\n\n" << PathIn1 << endl;
      filein = new TFile(PathIn1, "");
     
      for(Int_t v=0; v < numPtV0; v++){
	if(syst==0){
	  fHistPhiDistr_master[m][v]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	  fHistPhiDistr_master[m][v]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v]);
	}
	fHistPhiDistr[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	fHistPhiDistr_solostat[m][v][syst]=(TH1D*)filein->Get("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	fHistPhiDistr[m][v][syst]->SetName("PhiDistr_m"+Smolt[m]+"_v"+SPtV0[v] +Form("_syst%i", syst));
	fHistPhiDistr_solostat[m][v][syst]->SetName("PhiDistr_solostat_m"+Smolt[m]+"_v"+SPtV0[v] +Form("_syst%i", syst));
	fHistPhiDistr_ratio[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistPhiDistr_ratio_m%i_v%i_syst%i",m,v,syst));
	fHistPhiDistr_Barlow[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistPhiDistr_Barlow_m%i_v%i_syst%i",m,v,syst));

	if(syst==0){	
	  fHistPhiDistr_ratio[m][v][syst]={0};
	  fHistPhiDistr_Barlow[m][v][syst]={0};
	}

	if(syst!=0){	
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX(); j++){

	    if(fHistPhiDistr_master[m][v]->GetBinContent(j)!=0){
	      fHistPhiDistr_ratio[m][v][syst]->SetBinContent(j,fHistPhiDistr[m][v][syst]->GetBinContent(j)/fHistPhiDistr_master[m][v]->GetBinContent(j));
	      fHistPhiDistr_ratio[m][v][syst]->SetBinError(j,0);
	    }
	    if(sqrt(TMath::Abs(pow(fHistPhiDistr[m][v][syst]->GetBinError(j),2)-pow(fHistPhiDistr_master[m][v]->GetBinError(j),2)))!=0){
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinContent(j,(fHistPhiDistr[m][v][syst]->GetBinContent(j)-fHistPhiDistr_master[m][v]->GetBinContent(j))/sqrt(TMath::Abs(pow(fHistPhiDistr[m][v][syst]->GetBinError(j),2)-pow(fHistPhiDistr_master[m][v]->GetBinError(j),2))));
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinError(j,0);
	    }
	    else{
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinContent(j,0);
	      fHistPhiDistr_Barlow[m][v][syst]->SetBinError(j,0);
	      //	      cout << "mult " << m << " v " << v << " sys" << syst<< "j " << j << " denominatore Barlow = 0" << endl;
	    }
	  }
	}
	
	Int_t count=0;
	Correlated[m][v][0]=kTRUE; 
	BarlowPassed[m][v][0]=kTRUE;
	fHistPhiDistr_Corr[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistPhiDistr_Corr_m%i_v%i_syst%i",m,v,syst));
	//	fHistPhiDistr_Corr[m][v][syst]->SetBinContent(j,0);
	if (syst!=0) {
	  for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){ //escludo bin ai bordi
	    if (TMath::Abs(fHistPhiDistr_Barlow[m][v][syst]->GetBinContent(j)) >= 2) count++;
	  }
	  if (count > 0.05* (fHistPhiDistr_Barlow[m][v][syst]->GetNbinsX()-5)) BarlowPassed[m][v][syst]=kTRUE;
	  if ( BarlowPassed[m][v][syst]==kTRUE) {
	    //	    cout << "systematic effect n. " << SSyst[syst]<< " for m " << Smolt[m] << " and v " << SPtV0[v] << " has failed the Barlow check within 2 sigmas and is not significant" << endl;
	    NumberBarlowPassed[m][v]+=1;
	    for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){
	      Mean[m][v][syst] += fHistPhiDistr_ratio[m][v][syst]->GetBinContent(j);
	    }
	    Mean[m][v][syst]=Mean[m][v][syst]/fHistPhiDistr[m][v][syst]->GetNbinsX();
	    for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){
	      Sigma[m][v][syst] += pow(fHistPhiDistr_ratio[m][v][syst]->GetBinContent(j)-Mean[m][v][syst],2);
	    }
	    Sigma[m][v][syst]=sqrt(Sigma[m][v][syst]/(fHistPhiDistr[m][v][syst]->GetNbinsX()-5)/(fHistPhiDistr[m][v][syst]->GetNbinsX()-5+1));
	    if(TMath::Abs((Mean[m][v][syst]-1))>numSigmaCorr*Sigma[m][v][syst]) Correlated[m][v][syst]=kTRUE;
	    if (Correlated[m][v][syst]==kTRUE){
	      NumberCorr[m][v]+=1;
	    }
	  
	    Int_t numPhiBins = fHistPhiDistr[m][v][syst]->GetNbinsX();


	    for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX()-5; j++){
	      SigmaBarlow[m][v][syst][j]=(fHistPhiDistr[m][v][syst]->GetBinContent(j)-fHistPhiDistr_master[m][v]->GetBinContent(j))/sqrt(12.);
	      //	if (TMath::Abs(SigmaBarlow[m][v][syst][j]) > 1000) cout << m << " " << v << " " << syst <<  " " << j<< " " <<fHistPhiDistr[m][v][syst]->GetBinContent(j)<< "  " << fHistPhiDistr_master[m][v]->GetBinContent(j) << "  " << SigmaBarlow[m][v][syst][j] << endl;
	    }
	    if (Correlated[m][v][syst]==kFALSE){
	      for(Int_t j=1; j<= fHistPhiDistr[m][v][syst]->GetNbinsX(); j++){//errore pari a zero in tutti i bin
		fHistPhiDistr_Corr[m][v][syst]->SetBinContent(j,0);
	      }
	    }
	      
	  }
	}
      } //end loop on pt v0
    }//end loop on syst

    Int_t numPhiBins = fHistPhiDistr[0][0][0]->GetNbinsX();
 
    cout << "\n\ncalcolo l'errore sistematico (non correlato in deltaPhi) e lo sommo all'errore statistico in quadratura..." << endl;
    for(Int_t v=0; v < numPtV0; v++){
      for(Int_t syst=0; syst<numSysV0+numsystPhi; syst++){
	//	cout << " v " << v << " syst " << syst << endl; 
	if (syst==9 || syst == avoidthissyst) continue;
	fHistSigmaSyst[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistSigmaSyst_m%i_v%i_syst%i",m,v,syst));
	fHistSigmaSystNSmoothed[m][v][syst]=(TH1D*)fHistPhiDistr[m][v][syst]->Clone(Form("fHistSigmaSyst_notsmoothed_m%i_v%i_syst%i",m,v,syst));
	for(Int_t j=1; j <= numPhiBins; j++){
	  fHistSigmaSyst[m][v][syst]->SetBinContent(j, 0);
	  fHistSigmaSyst[m][v][syst]->SetBinError(j, 0);
	  fHistSigmaSystNSmoothed[m][v][syst]->SetBinContent(j, 0);
	  fHistSigmaSystNSmoothed[m][v][syst]->SetBinError(j, 0);
	}
	if (Correlated[m][v][syst]==kTRUE)continue; //riempio istogrammi solo se errore è non correlato in delta phi
	//	cout << " v " << v << " syst " << syst << endl; 
	for(Int_t j=1; j <= numPhiBins; j++){
	  //	  cout << " m" << m << " v " << v << " syst " << syst << "j " << j << " stat error " << fHistPhiDistr[m][v][syst]->GetBinError(j)/fHistPhiDistr[m][v][syst]->GetBinContent(j) << " sist error " << TMath::Abs(SigmaBarlow[m][v][syst][j])/fHistPhiDistr[m][v][syst]->GetBinContent(j) << endl; 
	  fHistSigmaSyst[m][v][syst]->SetBinContent(j, 0);
	  fHistSigmaSyst[m][v][syst]->SetBinError(j, 0);
	  //cout << " j" << j << " v " << v << " syst " << syst << endl; 
	  if (j < numPhiBins-5){
	    fHistSigmaSystNSmoothed[m][v][syst]->SetBinContent(j, TMath::Abs(SigmaBarlow[m][v][syst][j]/fHistPhiDistr_master[m][v]->GetBinContent(j)) );
	    fHistSigmaSyst[m][v][syst]->SetBinContent(j, TMath::Abs(SigmaBarlow[m][v][syst][j]));
	    fHistSigmaSyst[m][v][syst]->SetBinError(j, 0);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetBinError(j, 0);
	  }
	}
	fHistSigmaSyst[m][v][syst]->Smooth();
      }//end loop on syst

      for(Int_t j=1; j < numPhiBins-5; j++){
	for(Int_t syst=0; syst<numSysV0+numsystPhi; syst++){
	  if (syst==9 || syst == avoidthissyst) continue;
	  if (Correlated[m][v][syst]==kTRUE)continue;
	  SigmaSyst[m][v][j]+= pow(fHistSigmaSyst[m][v][syst]->GetBinContent(j),2);
	}
      
	SigmaSyst[m][v][j]= sqrt(SigmaSyst[m][v][j]);
	fHistPhiDistr[m][v][0]->SetBinError(j,sqrt(pow(fHistPhiDistr_master[m][v]->GetBinError(j),2) +pow(SigmaSyst[m][v][j],2)));  
	fHistSigmaSystNSmoothed[m][v][0]->SetBinContent(j,sqrt(pow(fHistPhiDistr_master[m][v]->GetBinError(j),2) +pow(SigmaSyst[m][v][j],2))/TMath::Abs(fHistPhiDistr_master[m][v]->GetBinContent(j) ));
	//	fHistSigmaSystNSmoothed[m][v][0]->SetBinContent(j,sqrt(pow(fHistPhiDistr_master[m][v]->GetBinError(j),2) )/TMath::Abs(fHistPhiDistr_master[m][v]->GetBinContent(j) ));
	fHistSigmaSystNSmoothed[m][v][0]->SetBinError(j,0);
      }
    }//end loop on v

    for(Int_t syst=0; syst< numSysV0 + numsystPhi; syst++){    
      Int_t count_corr=0;
      for(Int_t v=0; v< numPtV0; v++){    
	if(Correlated[m][v][syst]==kFALSE) continue; 
	count_corr++;
	//	cout << count_corr << endl;
      }
      if (count_corr==numPtV0) {
	CorrelatedBis[m][syst]=kTRUE;
	CorrelatedBisSys[m]++;
      }
    }
    
    fHistSpectrum_master_bis[m]=(TH1D*)fileAngularCorrMaster->Get("fHistSpectrum_"+Smolt[m]+"_eta<0.5");
    cout << "\n\ncalcolo l'errore sistematico (correlato in deltaPhi)" << endl;
    //    cout << numSyst- numsystPhi-numSysV0+1 << endl;
    
    for (Int_t syst1=0; syst1< numSyst- numsystPhi-numSysV0+1; syst1++){ //to consider systematics which are only related to bin counting
  
      //      cout << " syst1 " << syst1 << endl;
      
      for(Int_t v=0; v< numPtV0; v++){ 
	//	cout << "syst " << syst1 << "  " << 	NSpectrumFinal[m][v][syst1]<<"  " <<	NSpectrumErrorSoloStatFinal[m][v][syst1]<<  endl;	
	// NSpectrumFinal[m][v][syst1]=0;      
	// NSpectrumErrorFinal[m][v][syst1]=0;	
	// NSpectrumErrorSoloStatFinal[m][v][syst1]=0;
	
	for(Int_t syst=0; syst< numsystPhi + numSysV0; syst++){       
	  //	  cout << "syst " << syst << "  " << 	NSpectrum[m][v][0]<<"  " <<	NSpectrumErrorSoloStat[m][v][0]<<  endl;	
	  NSpectrum[m][v][syst]=0;      
	   NSpectrumError[m][v][syst]=0;	
	  NSpectrumErrorSoloStat[m][v][syst]=0;
	//  cout << "syst " << syst << "  " << 	NSpectrum[m][v][0]<<endl;//"  " <<	NSpectrumErrorSoloStat[m][v][syst]<<  endl;	
	}
	
	//	 cout <<  	NSpectrum[m][v][0]<<"  " <<	NSpectrumErrorSoloStat[m][v][0]<<  endl;	
		MiaNonna[m][v]=0;
	
	for(Int_t syst=0; syst< numsystPhi + numSysV0; syst++){    
	  cout << "\n\n\n ************************************syst "<<syst << endl;
	  if (syst==9 || syst == avoidthissyst) continue;
	  //inserisco errori sist correlati in delta phi nello spettro relativo a syst =0
	  
	  if ((Correlated[m][v][syst]==kTRUE && syst< 12) ){
	    // cout << " sono un sist correlato " << syst << endl;
	    fHistPhiDistr_CorrClone[m][v][syst]=(TH1D*)fHistPhiDistr_Corr[m][v][syst]->Clone(Form("fHistPhiDistr_Corr_Clone_m%i_v%i_syst%i",m,v,syst));
	    for(Int_t j=fHistPhiDistr_Corr[m][v][syst]->FindBin(ALowBinFit[syst1]); j<= fHistPhiDistr_Corr[m][v][syst]->FindBin(AUpBinFit[syst1]); j++){
	   fHistPhiDistr_CorrClone[m][v][syst]->SetBinContent(j, fHistPhiDistr_Corr[m][v][0]->GetBinContent(j)+ SigmaBarlow[m][v][syst][j]); //per syst = 0 è zero
	   NSpectrum[m][v][syst]+=fHistPhiDistr_CorrClone[m][v][syst]->GetBinContent(j);
	   //	   if (syst==0)	cout << "  v " << v << "syst " << syst << "  " << 	NSpectrum[m][v][0]<< endl;	
	       if (syst ==0) NSpectrumError[m][v][syst]+=pow(fHistPhiDistr[m][v][syst]->GetBinError(j),2);
	       //	        if (syst ==0) NSpectrumErrorSoloStat[m][v][syst]+=pow(fHistPhiDistr_solostat[m][v][0]->GetBinError(j),2);
	       if (syst ==0) NSpectrumErrorSoloStat[m][v][syst]+=pow(fHistPhiDistr_master[m][v]->GetBinError(j),2);
	    }
	     MiaNonna[m][v]+=pow(TMath::Abs(NSpectrum[m][v][0] -NSpectrum[m][v][syst]),2);
	  }
	 
	} //end loop syst
	
	NSpectrumErrorSistUnCorrPhi[m][v][syst1]=sqrt(NSpectrumError[m][v][0]);
	// Int_t syst=0;	
	NSpectrumErrorSoloStatFinal[m][v][syst1]=sqrt(NSpectrumErrorSoloStat[m][v][0]);
	
       	NSpectrumFinal[m][v][syst1]= NSpectrum[m][v][0];
	//	cout <<	NSpectrumFinal[m][v][syst1] << endl;
	
	NSpectrumErrorFinal[m][v][syst1]=sqrt(NSpectrumError[m][v][0]);
		
	NSpectrumErrorFinal[m][v][syst1]= sqrt(pow(NSpectrumErrorFinal[m][v][syst1],2) +MiaNonna[m][v]);
	MiaNonna[m][v]=MiaNonna[m][v]/NTrigger[m];
	 NSpectrumFinal[m][v][syst1]= NSpectrumFinal[m][v][syst1]/NTrigger[m];
	 NSpectrumErrorFinal[m][v][syst1]= NSpectrumErrorFinal[m][v][syst1]/NTrigger[m];
	 NSpectrumErrorSistUnCorrPhi[m][v][syst1]= NSpectrumErrorSistUnCorrPhi[m][v][syst1]/NTrigger[m];
	 NSpectrumErrorSoloStatFinal[m][v][syst1]=	NSpectrumErrorSoloStatFinal[m][v][syst1]/NTrigger[m];
	
	if(v==numPtV0-1){//l'ultimo bin è 4 volte più largo degli altri
	  NSpectrumFinal[m][v][syst1]= NSpectrumFinal[m][v][syst1]/4;
	  NSpectrumErrorFinal[m][v][syst1]= NSpectrumErrorFinal[m][v][syst1]/4;
	  NSpectrumErrorSistUnCorrPhi[m][v][syst1]= 	NSpectrumErrorSistUnCorrPhi[m][v][syst1]/4;
	  NSpectrumErrorSoloStatFinal[m][v][syst1]= 	NSpectrumErrorSoloStatFinal[m][v][syst1]/4;
	  MiaNonna[m][v]=MiaNonna[m][v]/4;
	}
	
	if(v==1 || v==2 || v==3 || v==4){//bin larghi 0.5
	  NSpectrumFinal[m][v][syst1]= NSpectrumFinal[m][v][syst1]*2;
	  NSpectrumErrorFinal[m][v][syst1]= NSpectrumErrorFinal[m][v][syst1]*2;
	  NSpectrumErrorSistUnCorrPhi[m][v][syst1]= 	NSpectrumErrorSistUnCorrPhi[m][v][syst1]*2;
	  NSpectrumErrorSoloStatFinal[m][v][syst1]= 	NSpectrumErrorSoloStatFinal[m][v][syst1]*2;
	  MiaNonna[m][v]=MiaNonna[m][v]*2;
	}
	fHistSpectrum[m][syst1]->SetBinContent(v+1, NSpectrumFinal[m][v][syst1]);
	fHistSpectrum[m][syst1]->SetBinError(v+1, NSpectrumErrorFinal[m][v][syst1]);
	
	cout << "errore statistico " <<    fHistSpectrum_master_bis[m]->GetBinError(v+1)<< endl;
	cout << "errore statistico qui calcolato" <<     NSpectrumErrorSoloStatFinal[m][v][syst1]/ fHistSpectrum_master_bis[m]->GetBinContent(v+1)<< endl;
	cout << "errore sist scorrelato in delta phi (+ stat)" << NSpectrumErrorSistUnCorrPhi[m][v][syst1]/ fHistSpectrum_master_bis[m]->GetBinContent(v+1)<< endl;
	cout << "errore sist correlato in delta phi " << sqrt(MiaNonna[m][v])/ fHistSpectrum_master_bis[m]->GetBinContent(v+1)<< endl; 
	 cout << "errore sist totale(da phi) + stat" << NSpectrumErrorFinal[m][v][syst1]/ fHistSpectrum_master_bis[m]->GetBinContent(v+1)<< endl; 
	
	 cout << "m " << m << " v " << v << " syst1 " << syst1<<endl;
	 //	cout <<  "  " << 	NSpectrumFinal[m][v][syst1] << endl;
	 cout <<"  " <<  	NSpectrumErrorFinal[m][v][syst1] << endl;
	
      }//end loop on v
    }//end loop syst 1

  
    fHistSpectrum_master[m]=(TH1D*)    fHistSpectrum[m][0]->Clone(Form("fHistSpectrum_master_%i",m));
    fHistSpectrum_masterSystCorr[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterSystCorr"+Smolt[m]);
    fHistSpectrum_masterSystUnCorr[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterSystUnCorr"+Smolt[m]);
    fHistSpectrum_masterTotalUnCorr[m]=      (TH1D*)fHistSpectrum[m][0]->Clone("fHistSpectrum_masterTotalUnCorr"+Smolt[m]);
    for(Int_t syst1=0; syst1<numSyst-numSysV0-numsystPhi+1; syst1++){
      fHistSpectrum_ratio[m][syst1]=      (TH1D*)fHistSpectrum[m][syst1]->Clone("fHistSpectrum_ratio"+Smolt[m]+Form("_SysBC%i",syst1));
      fHistSpectrum_Barlow[m][syst1]=      (TH1D*)fHistSpectrum[m][syst1]->Clone("fHistSpectrum_Barlow"+Smolt[m]+Form("_SysBC%i",syst1));
    }

   
    //Barlow check for spectrum
    for(Int_t syst=0; syst< numSyst-numsystPhi - numSysV0+1; syst++){    
      //      cout << "\n\n\n ************************************syst "<<syst << endl;
      if(syst==0){	
	fHistSpectrum_ratio[m][syst]={0};
	fHistSpectrum_Barlow[m][syst]={0};
      }
      if(syst!=0){	
	for(Int_t j=1; j<= numPtV0; j++){

	  if(fHistSpectrum_master[m]->GetBinContent(j)!=0){
	    fHistSpectrum_ratio[m][syst]->SetBinContent(j,fHistSpectrum[m][syst]->GetBinContent(j)/fHistSpectrum_master[m]->GetBinContent(j));
	    fHistSpectrum_ratio[m][syst]->SetBinError(j,0);
	  }
	  if(sqrt(TMath::Abs(pow(fHistSpectrum[m][syst]->GetBinError(j),2)-pow(fHistSpectrum_master[m]->GetBinError(j),2)))!=0){
	    fHistSpectrum_Barlow[m][syst]->SetBinContent(j,(fHistSpectrum[m][syst]->GetBinContent(j)-fHistSpectrum_master[m]->GetBinContent(j))/sqrt(TMath::Abs(pow(fHistSpectrum[m][syst]->GetBinError(j),2)-pow(fHistSpectrum_master[m]->GetBinError(j),2))));
	    fHistSpectrum_Barlow[m][syst]->SetBinError(j,0);
	  }
	  else{
	    fHistSpectrum_Barlow[m][syst]->SetBinContent(j,0);
	    fHistSpectrum_Barlow[m][syst]->SetBinError(j,0);
	  }

	}
      }
	
      Int_t count=0;
      CorrelatedSpectrum[m][0]=kTRUE;
      fHistSpectrum_Corr[m][syst]=(TH1D*)fHistSpectrum[m][syst]->Clone(Form("fHistSpectrum_Corr_m%i_syst%i",m,syst));
      //	fHistSpectrum_Corr[m][syst]->SetBinContent(j,0);
      if (syst!=0) {
	for(Int_t j=1; j<= numPtV0; j++){
	  if (TMath::Abs(fHistSpectrum_Barlow[m][syst]->GetBinContent(j)) >= 2) count++;
	}
	if (count >=1) BarlowPassedSpectrum[m][syst]=kTRUE;

	if ( BarlowPassedSpectrum[m][syst]==kTRUE) {
	  NumberBarlowPassedSpectrum[syst]+=1;
	  for(Int_t j=1; j<= fHistSpectrum[m][syst]->GetNbinsX(); j++){
	    MeanSpectrum[m][syst] += fHistSpectrum_ratio[m][syst]->GetBinContent(j);
	  }
	  MeanSpectrum[m][syst]=MeanSpectrum[m][syst]/fHistSpectrum[m][syst]->GetNbinsX();
	  for(Int_t j=1; j<= fHistSpectrum[m][syst]->GetNbinsX(); j++){
	    SigmaSpectrum[m][syst] += pow(fHistSpectrum_ratio[m][syst]->GetBinContent(j)-MeanSpectrum[m][syst],2);
	  }
	  SigmaSpectrum[m][syst]=sqrt(SigmaSpectrum[m][syst]/(fHistSpectrum[m][syst]->GetNbinsX())/(fHistSpectrum[m][syst]->GetNbinsX()+1));
	  if(TMath::Abs((MeanSpectrum[m][syst]-1))>numSigmaCorr*SigmaSpectrum[m][syst]) CorrelatedSpectrum[m][syst]=kTRUE;
	  
	  //	    if (CorrelatedSpectrum[m][syst]==kFALSE){
	  //	  cout << "  m " << m<< "  syst " << syst << endl;
	  fHistSigmaBarlowSpectrum[m][syst]=(TH1D*)fHistSpectrum[m][syst]->Clone(Form("fHistSpectrum_SigmaBarlow_m%i_syst%i",m,syst));
	  for(Int_t j=1; j<= fHistSpectrum[m][syst]->GetNbinsX(); j++){
	    fHistSigmaBarlowSpectrum[m][syst]->SetBinContent(j,(fHistSpectrum[m][syst]->GetBinContent(j)-fHistSpectrum_master[m]->GetBinContent(j))/sqrt(12.));
	    fHistSpectrum_Corr[m][syst]->SetBinContent(j,0);
	  }
	    
	  //}
	  if (CorrelatedSpectrum[m][syst]==kTRUE){
	    NumberCorrSpectrum[m]+=1;
	  }
	}
      }

    }//end loop syst


    
    cout << "setting sist errors to spectrum " << endl;
    cout << "\n\nmult "<< m << endl;
    for(Int_t j=0; j< numPtV0; j++){
	SigmaSystSpectrumUnCorr[m][j]={0};
	SigmaTotalSpectrumUnCorr[m][j]={0};
	SigmaSystSpectrumCorr[m][j]={0};
      //cout << "************\n\n" << endl;
      for(Int_t syst=1; syst< numSyst-numsystPhi - numSysV0+1; syst++){    

	  cout << " j " <<j << " syst  " <<syst << endl;
	//if (syst==9 || syst == avoidthissyst) continue;
	//	if (CorrelatedBis[m][syst]==kFALSE && syst< 12) continue; //remember: for syst=0 was put kTRUE
	if (CorrelatedSpectrum[m][syst]==kTRUE) continue;
	if ( BarlowPassedSpectrum[m][syst]==kFALSE) continue;
	cout << " componente non correlata m" << m << " v " << j << " syst " << syst << endl;
	cout << fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1)        	<< endl;
	SigmaSystSpectrumUnCorr[m][j]+= pow(fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1),2);
      }

      SigmaSystSpectrumUnCorr[m][j]+= pow(NSpectrumErrorFinal[m][j][0],2); //total uncorr = stat + sist uncorr

      SigmaTotalSpectrumUnCorr[m][j]=sqrt(SigmaSystSpectrumUnCorr[m][j]); //total uncorr

      SigmaSystSpectrumUnCorr[m][j]=       SigmaSystSpectrumUnCorr[m][j] -pow(  NSpectrumErrorSoloStatFinal[m][j][0],2); //- stat

      SigmaSystSpectrumUnCorr[m][j]= sqrt(SigmaSystSpectrumUnCorr[m][j]);    
      cout << " m " << m << " j " <<       SigmaSystSpectrumUnCorr[m][j] << "  " <<       SigmaTotalSpectrumUnCorr[m][j]<< endl;
       cout << "\n\n" << endl;
      for(Int_t syst=1; syst< numSyst-numsystPhi - numSysV0+1; syst++){    
	cout <<m << "  " << j << "  " <<syst << endl;
	if (CorrelatedSpectrum[m][syst]==kFALSE) continue;
	if ( BarlowPassedSpectrum[m][syst]==kFALSE) continue;
	cout << " componente correlata " << m << " v " << j << " syst " << syst << endl;
	//cout << fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1)        	<< endl;
	SigmaSystSpectrumCorr[m][j]+= pow(fHistSigmaBarlowSpectrum[m][syst]->GetBinContent(j+1),2);
      }
      SigmaSystSpectrumCorr[m][j]=sqrt(SigmaSystSpectrumCorr[m][j]);
      fHistSpectrum_masterSystCorr[m]->SetBinError(j+1,SigmaSystSpectrumCorr[m][j]);
      fHistSpectrum_masterSystUnCorr[m]->SetBinError(j+1,SigmaSystSpectrumUnCorr[m][j]);
      fHistSpectrum_masterTotalUnCorr[m]->SetBinError(j+1,SigmaTotalSpectrumUnCorr[m][j]);
      cout << " m " << m << " j " <<       SigmaSystSpectrumCorr[m][j] << endl;
    }
    
  }//end loop m


  cout << " i draw canvas " << endl;
  TCanvas *canvasSys[4];
  for (Int_t i=0; i < 3; i++){
    canvasSys[i]=new TCanvas(Form("canvas%i",i),Form("canvas%i",i), 1300, 1000);
    canvasSys[i]->Divide(3,2);
    if (i==2)     canvasSys[i]->Divide(1,1);
  }


  //yield vs multiplicity
  Double_t Yield[nummolt]={0};
  Double_t YieldErrTotal[nummolt]={0};
  Double_t YieldErrStat[nummolt]={0};
  Double_t YieldErrSistUnCorr[nummolt]={0};
  Double_t YieldErrSistCorr[nummolt]={0};
  TH1D* fHistYieldvsErrTotal=new TH1D ("fHistYieldvsErrTot","fHistYieldvsErrTot",300,0,30);
  TH1D* fHistYieldvsErrSoloStat=new TH1D ("fHistYieldvsErrSoloStat","fHistYieldvsErrSoloStat",300,0,30);
  TH1D* fHistYieldvsErrSoloSist=new TH1D ("fHistYieldvsErrSoloSist","fHistYieldvsErrSoloSist",300, 0, 30);
  Double_t mult[nummolt]={21.2, 16.17, 11.4625, 7.135, 3.33};
  for(Int_t m=0; m < nummolt; m++){
    for(Int_t v =0; v < numPtV0; v++){
      Yield[m]+=   fHistSpectrum_master[m]->GetBinContent(v+1);
      YieldErrStat[m]+= pow(   fHistSpectrum_master_bis[m]->GetBinError(v+1),2);
      YieldErrSistUnCorr[m]+=pow(    fHistSpectrum_masterSystUnCorr[m]->GetBinError(v+1),2);
      YieldErrSistCorr[m]+=pow(    fHistSpectrum_masterSystCorr[m]->GetBinError(v+1),2);
    
    }
    YieldErrTotal[m]= YieldErrStat[m]+ YieldErrSistUnCorr[m]+ YieldErrSistCorr[m];
    YieldErrStat[m]=sqrt( YieldErrStat[m]);
    YieldErrSistUnCorr[m]=sqrt( YieldErrSistUnCorr[m]);
    YieldErrSistCorr[m]=sqrt( YieldErrSistCorr[m]);
    YieldErrTotal[m]=sqrt( YieldErrTotal[m]);
    fHistYieldvsErrSoloSist->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);
    fHistYieldvsErrSoloStat->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);  
    fHistYieldvsErrTotal->SetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]), Yield[m]);
    fHistYieldvsErrSoloSist->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), sqrt(pow(YieldErrSistUnCorr[m],2) + pow(YieldErrSistCorr[m],2)));
    fHistYieldvsErrSoloStat->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrStat[m]);  
    fHistYieldvsErrTotal->SetBinError(fHistYieldvsErrTotal->FindBin(mult[m]), YieldErrTotal[m]);
    cout << "************m  " << m << endl;
    cout << "yield"<<     fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
    cout << "error stat rel"<<     fHistYieldvsErrSoloStat->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
    cout << "error sist rel"<<     fHistYieldvsErrSoloSist->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
 cout << "error total rel"<<     fHistYieldvsErrTotal->GetBinError(fHistYieldvsErrTotal->FindBin(mult[m]))/fHistYieldvsErrSoloSist->GetBinContent(fHistYieldvsErrTotal->FindBin(mult[m]))<< endl;
  }
  TF1 *expo = new TF1("expo", "expo", 0,30);
  fHistYieldvsErrTotal->Fit(expo,"R+");

  canvasSys[2]->cd();
    fHistYieldvsErrSoloStat->GetYaxis()->SetRangeUser(0,0.2);
    fHistYieldvsErrSoloStat->SetMarkerStyle(20);
    fHistYieldvsErrSoloStat->SetLineColor(1);
    fHistYieldvsErrSoloStat->SetMarkerColor(1);
    fHistYieldvsErrSoloStat->Draw("samee");
    fHistYieldvsErrSoloSist->GetYaxis()->SetRangeUser(0,0.2);
    fHistYieldvsErrSoloSist->SetMarkerStyle(20);
    fHistYieldvsErrSoloSist->SetLineColor(1);
    fHistYieldvsErrSoloSist->SetMarkerColor(2);
    fHistYieldvsErrSoloSist->SetFillStyle(0);
    fHistYieldvsErrSoloSist->Draw("samee2");
   // fHistYieldvsErrTotal->GetYaxis()->SetRangeUser(0,0.2);
   // fHistYieldvsErrTotal->SetMarkerStyle(20);
   // fHistYieldvsErrTotal->SetLineColor(1);
   // fHistYieldvsErrTotal->SetMarkerColor(3);
   // fHistYieldvsErrTotal->SetFillStyle(0);
   // fHistYieldvsErrTotal->Draw("samee2");
    expo->Draw("same");
  
  cout << "\n\n\ndrawing canvas " << endl;
  ///////////////////////////////////////////
  //********************************drawing canvas***************************************************************************


  TF1* retta = new TF1("retta", "pol0", -10, 10);
  retta->FixParameter(0,1);
  retta->SetLineColor(1);
  retta->SetLineWidth(0.1);
  TF1* retta1 = new TF1("retta1", "pol0", -10, 10);
  retta1->FixParameter(0,0);
  retta1->SetLineColor(1);
  retta1->SetLineWidth(0.1);



  for(Int_t m=0; m< nummolt +1; m++){
    legendCorrBis[m]= new TLegend(0.6, 0.6, 1, 1);
    legendCorrBis[m]->SetHeader("Selezioni applicate");
    canvasSys[0]->cd(m+1);
    //    cout << "m " << m << endl;
    Int_t count4=0;
    for(Int_t syst=0; syst< numSyst-numsystPhi - numSysV0 +1; syst++){
      //      cout << "syst " << syst << "  " <<CorrelatedBis[m][syst] << endl;
      // if (CorrelatedBis[m][syst]==kFALSE) continue;
      // if (syst==9 || syst == avoidthissyst) continue; 
      //      if (syst>0 && syst < 12) continue; 
      fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.04);
      fHistSpectrum[m][syst]->SetMarkerStyle(Marker[syst]);
      fHistSpectrum[m][syst]->SetLineColor(Color[syst]);
      fHistSpectrum[m][syst]->SetMarkerColor(Color[syst]);
      fHistSpectrum[m][syst]->Draw("samee");
      //if (CorrelatedBis[m][syst]==kTRUE || syst>=12){
	legendCorrBis[m]->AddEntry(fHistSpectrum[m][syst],SSyst[syst],"pel");   
	count4++;
	//}
      //cout << count4 << "  " << CorrelatedBisSys[m]+numSyst-numsystPhi-numSysV0 << endl;
      if(syst==1) legendCorrBis[m]->Draw();
    }
    
  }

  legendErrorSpectrum= new TLegend(0.6, 0.6, 1, 1);
  legendErrorSpectrum->SetHeader("Type of error");
  
  for(Int_t m=0; m< nummolt +1; m++){
    canvasSys[1]->cd(m+1);
    fHistSpectrum_master_bis[m]->GetYaxis()->SetRangeUser(0,0.04);
    fHistSpectrum_master_bis[m]->SetMarkerStyle(20);
    fHistSpectrum_master_bis[m]->SetLineColor(1);
    fHistSpectrum_master_bis[m]->SetMarkerColor(1);
    fHistSpectrum_master_bis[m]->Draw("samee");
    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,0.04);
    fHistSpectrum_masterSystCorr[m]->SetMarkerStyle(20);
    fHistSpectrum_masterSystCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystCorr[m]->SetMarkerColor(2);
    fHistSpectrum_masterSystCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterSystCorr[m]->Draw("samee2");
    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetRangeUser(0,0.04);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerStyle(20);
    fHistSpectrum_masterSystUnCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerColor(3);
    fHistSpectrum_masterSystUnCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterSystUnCorr[m]->Draw("samee2");
    if(m==0){
      legendErrorSpectrum->AddEntry(fHistSpectrum_master_bis[m],SErrorSpectrum[0], "pel");  
      legendErrorSpectrum->AddEntry(fHistSpectrum_masterSystCorr[m],SErrorSpectrum[1], "pel");   
      legendErrorSpectrum->AddEntry(fHistSpectrum_masterSystUnCorr[m],SErrorSpectrum[2], "pel");   
    }
    legendErrorSpectrum->Draw();
  }

  // cout << "\n\n\n ******************drawing canvases******************** " << endl;
  TFile *fileout = new TFile(stringout ,"RECREATE");
  fileout->WriteTObject( canvasSys[0]);
  fileout->WriteTObject( canvasSys[1]);
  fileout->WriteTObject( canvasSys[2]);
  
 
  canvasSpectrumSysBis[0]=new TCanvas("canvasSpectrumSysBis0","canvasSpectrumSysBis0", 1300, 1000);
  canvasSpectrumSysBis[1]=new TCanvas("canvasSpectrumSysBis1","canvasSpectrumSysBis1", 1300, 1000);
  canvasSpectrumSysBis[0]->Divide(3,3);
  canvasSpectrumSysBis[1]->Divide(3,3);

  for(Int_t m=0; m< nummolt +1; m++){
    canvasSpectrumSys[2*m]=new TCanvas(Form("canvasSpectrumSys0_%i",m),Form("canvasSpectrumSys0_%i",m), 1300, 1000);
    canvasSpectrumSys[2*m+1]=new TCanvas(Form("canvasSpectrumSys1_%i",m),Form("canvasSpectrumSys1_%i",m), 1300, 1000);
    canvasSpectrumSys[2*m]->Divide(4,3);
    canvasSpectrumSys[2*m+1]->Divide(4,3);

    canvasPhiSys[4*m]=new TCanvas(Form("canvasSys0_%i",m),Form("canvasSys0_%i",m), 1300, 1000);
    canvasPhiSys[4*m+1]=new TCanvas(Form("canvasSys1_%i",m),Form("canvasSys1_%i",m), 1300, 1000);
    canvasPhiSys[4*m+2]=new TCanvas(Form("canvasSys2_%i",m),Form("canvasSys2_%i",m), 1300, 1000);
    canvasPhiSys[4*m+3]=new TCanvas(Form("canvasSys3_%i",m),Form("canvasSys3_%i",m), 1300, 1000);

    canvasPhiSys[4*m]->Divide(2,3);
    canvasPhiSys[4*m+1]->Divide(2,3);
    canvasPhiSys[4*m+2]->Divide(2,3);
    canvasPhiSys[4*m+3]->Divide(2,3);


  }
  for(Int_t m=0; m< nummolt +1; m++){
    legend_Bpassed_Spectrum[m]= new TLegend(0.6, 0.6, 1, 1);
    legend_Corr_Spectrum[m]= new TLegend(0.6, 0.6, 1, 1);

    //
    Int_t counter=0;
    Int_t counter2=0;
    for(Int_t syst=0; syst<numSyst-numsystPhi - numSysV0+1; syst++){
	  
      // if (syst==9 || syst == avoidthissyst) continue;
      // if (CorrelatedBis[m][syst]==kFALSE) continue;
      // if (syst>0 && syst < 12) continue; 
      if (m<3)	  canvasSpectrumSysBis[0]->cd(m+1);
      else 	  canvasSpectrumSysBis[1]->cd(m+1-3);
      fHistSpectrum[m][syst]->GetYaxis()->SetRangeUser(0,0.04);
      fHistSpectrum[m][syst]->SetMarkerStyle(Marker[syst]);
      fHistSpectrum[m][syst]->SetLineColor(Color[syst]);
      fHistSpectrum[m][syst]->SetMarkerColor(Color[syst]);
      fHistSpectrum[m][syst]->Draw("samee");
      if(syst==numSyst-1) legendCorrBis[m]->Draw();

      if(BarlowPassedSpectrum[m][syst]==kTRUE) {
	counter++;
	if(syst!=0){
	  if (m<3) canvasSpectrumSysBis[0]->cd(6+m+1);
	  else canvasSpectrumSysBis[1]->cd(6+m+1-3);
	  fHistSpectrum_Barlow[m][syst]->GetYaxis()->SetRangeUser(-15,15);
	  fHistSpectrum_Barlow[m][syst]->SetMarkerStyle(Marker[syst]);
	  fHistSpectrum_Barlow[m][syst]->SetLineColor(Color[syst]);
	  fHistSpectrum_Barlow[m][syst]->SetMarkerColor(Color[syst]);
	  fHistSpectrum_Barlow[m][syst]->Draw("samep");
	  retta1->SetLineColor(1);
	  retta1->SetLineWidth(0.1);
	  retta1->Draw("same");
	  legend_Bpassed_Spectrum[m]->AddEntry(fHistSpectrum_Barlow[m][syst],SSyst[syst],"pel");   
	  if(counter==NumberBarlowPassedSpectrum[m]) legend_Bpassed_Spectrum[m]->Draw();

	  if (CorrelatedSpectrum[m][syst]==kTRUE){	  
	    counter2++;
	    if (m<3) canvasSpectrumSysBis[0]->cd(3+m+1);
	    else canvasSpectrumSysBis[1]->cd(3+m+1-3);
	    fHistSpectrum_ratio[m][syst]->GetYaxis()->SetRangeUser(-2,4);
	    fHistSpectrum_ratio[m][syst]->SetMarkerStyle(Marker[syst]);
	    fHistSpectrum_ratio[m][syst]->SetLineColor(Color[syst]);
	    fHistSpectrum_ratio[m][syst]->SetMarkerColor(Color[syst]);
	    fHistSpectrum_ratio[m][syst]->Draw("samep");
	    retta->SetLineColor(1);
	    retta->SetLineWidth(0.1);
	    retta->Draw("same");
	    legend_Corr_Spectrum[m]->AddEntry(fHistSpectrum_ratio[m][syst],SSyst[syst],"pel");   
	    //  cout << counter2 << " "<<NumberCorr[m] << endl;
	    if(counter2==NumberCorrSpectrum[m]) legend_Corr_Spectrum[m]->Draw();
	  }
	  
	}
      }
    }
  }


  for(Int_t m=0; m< nummolt +1; m++){
    for(Int_t v=0; v < numPtV0; v++){
      legend_corr[m][v]= new TLegend(0.6, 0.6, 1, 1);
      legend_corr[m][v]->SetHeader("Selezioni applicate");
      legend_Bpassed[m][v]= new TLegend(0.6, 0.6, 1, 1);
      legend_Bpassed[m][v]->SetHeader("Selezioni applicate");
      legend_UnCorr[m][v]= new TLegend(0.6, 0.6, 1, 1);
      legend_UnCorr[m][v]->SetHeader("Selezioni applicate");
	
      Int_t counter=0;
      Int_t counter2=0;
      Int_t counter_unc=0;
      //      cout << "******************** " << m << "  " << v << endl;
      for(Int_t syst=0; syst<numSysV0+numsystPhi; syst++){
	  
	if (syst==9 || syst == avoidthissyst) continue;
	if (v<2)	  canvasPhiSys[4*m]->cd(v+1);
	else 	  if (v<4)	  canvasPhiSys[4*m+1]->cd(v+1-2);
	else 	  if (v<6)	  canvasPhiSys[4*m+2]->cd(v+1-4);
	else	  canvasPhiSys[4*m+3]->cd(v+1-6);
	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,8000);
	if (m==5)	fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,8000);
	fHistPhiDistr[m][v][syst]->SetMarkerStyle(Marker[syst]);
	fHistPhiDistr[m][v][syst]->SetLineColor(Color[syst]);
	fHistPhiDistr[m][v][syst]->SetMarkerColor(Color[syst]);
	if (m==0 && v==0) legend3->AddEntry(fHistPhiDistr[m][v][syst],SSyst[syst],"pel");   
	fHistPhiDistr[m][v][syst]->Draw("samee");
	if(syst == numSysV0+numsystPhi-1) legend3->Draw();

	if(BarlowPassed[m][v][syst]==kTRUE) {
	  counter++;
	  if(syst!=0){
	    if (v<2)	  canvasPhiSys[4*m]->cd(4+v+1);
	    else 	  if (v<4)	  canvasPhiSys[4*m+1]->cd(4+v+1-2);
	    else 	  if (v<6)	  canvasPhiSys[4*m+2]->cd(4+v+1-4);
	    else	  canvasPhiSys[4*m+3]->cd(4+v+1-6);
	    fHistPhiDistr_Barlow[m][v][syst]->GetYaxis()->SetRangeUser(-15,15);
	    fHistPhiDistr_Barlow[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistPhiDistr_Barlow[m][v][syst]->SetLineColor(Color[syst]);
	    fHistPhiDistr_Barlow[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistPhiDistr_Barlow[m][v][syst]->Draw("samep");
	    retta1->Draw("same"); 
	    legend_Bpassed[m][v]->AddEntry(fHistPhiDistr_ratio[m][v][syst],SSyst[syst],"pel");   
	    if(counter==NumberBarlowPassed[m][v]) legend_Bpassed[m][v]->Draw();
	    //cout << counter << " "<<NumberBarlowPassed[m][v] << endl;

	    if (Correlated[m][v][syst]==kTRUE){	  
	      counter2++;
	      if (v<2)	  canvasPhiSys[4*m]->cd(2+v+1);
	      else 	  if (v<4)	  canvasPhiSys[4*m+1]->cd(2+v+1-2);
	      else 	  if (v<6)	  canvasPhiSys[4*m+2]->cd(2+v+1-4);
	      else	  canvasPhiSys[4*m+3]->cd(2+v+1-6);
	      fHistPhiDistr_ratio[m][v][syst]->GetYaxis()->SetRangeUser(-2,4);
	      fHistPhiDistr_ratio[m][v][syst]->SetMarkerStyle(Marker[syst]);
	      fHistPhiDistr_ratio[m][v][syst]->SetLineColor(Color[syst]);
	      fHistPhiDistr_ratio[m][v][syst]->SetMarkerColor(Color[syst]);
	      fHistPhiDistr_ratio[m][v][syst]->Draw("samep");
	      retta->SetLineColor(1);
	      retta->SetLineWidth(0.1);
	      retta->Draw("same");
	      legend_corr[m][v]->AddEntry(fHistPhiDistr_ratio[m][v][syst],SSyst[syst],"pel");   
	      //  cout << counter2 << " "<<NumberCorr[m][v] << endl;
	      if(counter2==NumberCorr[m][v]) legend_corr[m][v]->Draw();
	    }
	  
	  }
	

	  if (Correlated[m][v][syst]==kTRUE){ //remember: for syst=0 was put kTRUE
	    // cout << "\n\ndrawing canvas " << endl;
	    // cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	    if (v<4)  canvasSpectrumSys[2*m]->cd(v+1);
	    else  canvasSpectrumSys[2*m+1]->cd(v+1-4);
	    fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,8000);
	    if (m==5)	  fHistPhiDistr[m][v][syst]->GetYaxis()->SetRangeUser(-1000,20000);
	    fHistPhiDistr[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistPhiDistr[m][v][syst]->SetLineColor(Color[syst]);
	    fHistPhiDistr[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistPhiDistr[m][v][syst]->Draw("samee");
	    //	  if(counter2==NumberCorr[m][v]) legend_corr[m][v]->Draw();
	    if(syst==0) legend_corr[m][v]->Draw();
	 
	    /* 
	       if(syst!=0){
	       if (v<4)	  canvasSpectrumSys[2*m]->cd(4+v+1);
	       else	  canvasSpectrumSys[2*m+1]->cd(4+v+1-4);
	       fHistPhiDistr_ratio[m][v][syst]->GetYaxis()->SetRangeUser(-5,5);
	       fHistPhiDistr_ratio[m][v][syst]->SetMarkerStyle(Marker[syst]);
	       fHistPhiDistr_ratio[m][v][syst]->SetLineColor(Color[syst]);
	       fHistPhiDistr_ratio[m][v][syst]->SetMarkerColor(Color[syst]);
	       fHistPhiDistr_ratio[m][v][syst]->Draw("samee");

	       }
	    */
	    /*
	      if(syst==0){
	      //	    cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	      if (v<4)	  canvasSpectrumSys[2*m]->cd(4+v+1);
	      else	  canvasSpectrumSys[2*m+1]->cd(4+v+1-4);
	      fHistPhiDistr_master[m][v]->GetYaxis()->SetRangeUser(-1000,8000);
	      if (m==5) fHistPhiDistr_master[m][v]->GetYaxis()->SetRangeUser(-1000,20000);
	      fHistPhiDistr_master[m][v]->SetMarkerStyle(Marker[syst]);
	      fHistPhiDistr_master[m][v]->SetLineColor(Color[syst]);
	      fHistPhiDistr_master[m][v]->SetMarkerColor(Color[syst]);
	      fHistPhiDistr_master[m][v]->Draw("samee");

	      }
	    */


  
	  }

	   if (syst==0 || Correlated[m][v][syst]==kFALSE){ //remember: for syst=0 was put kTRUE
	    //  cout << "syst " << syst << endl;
	    counter_unc++;
	    if (v<4)	  canvasSpectrumSys[2*m]->cd(8+v+1);
	    else	  canvasSpectrumSys[2*m+1]->cd(8+v+1-4);
	    fHistSigmaSystNSmoothed[m][v][syst]->GetYaxis()->SetRangeUser(0,2);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetLineColor(Color[syst]);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistSigmaSystNSmoothed[m][v][syst]->SetTitle("Relative uncorrelated systematic uncertainty");
	    fHistSigmaSystNSmoothed[m][v][syst]->GetXaxis()->SetTitle("#Delta #phi");
	    fHistSigmaSystNSmoothed[m][v][syst]->GetYaxis()->SetTitle("Relative uncertainty");
	    fHistSigmaSystNSmoothed[m][v][syst]->Draw("same");
	    if (syst !=0 &&  Correlated[m][v][syst]==kFALSE)	    legend_UnCorr[m][v]->AddEntry(fHistSigmaSystNSmoothed[m][v][syst],SSyst[syst],"pel");   
	    else  legend_UnCorr[m][v]->AddEntry(fHistSigmaSystNSmoothed[m][v][syst],"total uncorrelated uncertainty","pel");   
	    if(counter_unc==(NumberBarlowPassed[m][v]-NumberCorr[m][v])) legend_UnCorr[m][v]->Draw();

	    //	cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	    if (v<4)	  canvasSpectrumSys[2*m]->cd(4+v+1);
	    else	  canvasSpectrumSys[2*m+1]->cd(4+v+1-4);
	    fHistSigmaSyst[m][v][syst]->GetYaxis()->SetRangeUser(0,400);
	    fHistSigmaSyst[m][v][syst]->SetMarkerStyle(Marker[syst]);
	    fHistSigmaSyst[m][v][syst]->SetLineColor(Color[syst]);
	    fHistSigmaSyst[m][v][syst]->SetMarkerColor(Color[syst]);
	    fHistSigmaSyst[m][v][syst]->Draw("sameep");
	    //	cout << " m " <<  m << " v " << v<< " syst " << syst << endl;
	   }
	}

      }
    }
  }//end loop on mult



  cout << "master val, master_bis val, master_bis error(stat), errore stat , errore stat + sist uncorr(deltaphi), errore stat +sist(delta phi) , errore uncorr (sist), error uncorr sist + stat, error corr (sist) "<< endl;
  for(Int_t m=0; m< nummolt +1; m++){
    cout << "*********************************mult " <<m << endl;
    for(Int_t j=0; j < numPtV0; j++){
      for(Int_t syst=0; syst<2; syst++){
	if(syst==9) continue;
	//	cout << "syst " << syst << "  " << fHistSpectrum_master[m]->GetBinContent(j+1) << "  " << fHistSpectrum[m][syst]->GetBinContent(j+1) <<  "  " <<  ALowBinFit[syst] << "  " << AUpBinFit[syst] << "  " << endl;//<<fHistPhiDistr_Corr[m][j][syst]->FindBin(ALowBinFit[syst])<<  endl;

      }
      cout << "^^^^" << j<< "        " << fHistSpectrum_master[m]->GetBinContent(j+1)<< "  " <<  "     " << fHistSpectrum_master_bis[m]->GetBinContent(j+1)<< "     " <<fHistSpectrum_master_bis[m]->GetBinError(j+1) << "  " << NSpectrumErrorSoloStatFinal[m][j][0]<< "  <<   " <<NSpectrumErrorSistUnCorrPhi[m][j][0] << "  " <<NSpectrumErrorFinal[m][j][0]<< " <= " << SigmaSystSpectrumUnCorr[m][j] << " < " << SigmaTotalSpectrumUnCorr[m][j] << " !=0 but not consid.     " << SigmaSystSpectrumCorr[m][j] << "          "   << endl;
    
      cout << "^^^^" << j<< "        " << fHistSpectrum_master[m]->GetBinContent(j+1)<< " stat rel:  " <<fHistSpectrum_master_bis[m]->GetBinError(j+1)/fHistSpectrum_master[m]->GetBinContent(j+1) <<" sist uncorr phi (rel)  " <<sqrt(pow(NSpectrumErrorSistUnCorrPhi[m][j][0],2) - pow(fHistSpectrum_master_bis[m]->GetBinError(j+1),2))/fHistSpectrum_master[m]->GetBinContent(j+1) << " sist phi (rel) " <<sqrt(pow(NSpectrumErrorFinal[m][j][0],2) - pow(fHistSpectrum_master_bis[m]->GetBinError(j+1),2))/fHistSpectrum_master[m]->GetBinContent(j+1)<< " sist uncorr pt (rel)" << SigmaSystSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1) << " sist tot (rel)" << SigmaTotalSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1) << "sist corr pt (rel)" << SigmaSystSpectrumCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1)  << endl;
  cout << "*************" << endl;
      //  break;
    }
  }
  cout << "end loop on m " << endl;


    canvasSys[3]=new TCanvas("canvas3","canvas3", 1300, 1000);
    canvasSys[3]->Divide(3,2);
    TLegend *legenderrorspectrum = new TLegend(0.6, 0.6,0.9,0.9);

  for(Int_t m=0; m < nummolt+1; m++){
    for(Int_t j =0; j < numPtV0; j++){
      fHistSpectrum_master_bis[m]->SetBinContent(j+1,NSpectrumErrorSoloStatFinal[m][j][0]/fHistSpectrum_master[m]->GetBinContent(j+1));
      //      fHistSpectrum_master_bis[m]->SetBinError(j+1,sqrt(pow(NSpectrumErrorFinal[m][j][0],2) - pow(fHistSpectrum_master_bis[m]->GetBinError(j+1),2))/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterSystCorr[m]->SetBinContent(j+1,SigmaSystSpectrumCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterSystUnCorr[m]->SetBinContent(j+1,SigmaSystSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1));
      fHistSpectrum_masterTotalUnCorr[m]->SetBinContent(j+1,SigmaTotalSpectrumUnCorr[m][j]/fHistSpectrum_master[m]->GetBinContent(j+1));
    }

    canvasSys[3]->cd(m+1);
    fHistSpectrum_master_bis[m]->GetYaxis()->SetRangeUser(0,0.35);
    fHistSpectrum_master_bis[m]->SetMarkerStyle(20);
    fHistSpectrum_master_bis[m]->SetLineColor(1);
    fHistSpectrum_master_bis[m]->SetMarkerColor(2);
    fHistSpectrum_master_bis[m]->Draw("same");
    if (m==0)legenderrorspectrum->AddEntry(    fHistSpectrum_master_bis[m], "stat " , "pel");
    fHistSpectrum_masterSystCorr[m]->GetYaxis()->SetRangeUser(0,0.35);
    fHistSpectrum_masterSystCorr[m]->SetMarkerStyle(20);
    fHistSpectrum_masterSystCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystCorr[m]->SetMarkerColor(2);
    fHistSpectrum_masterSystCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterSystCorr[m]->Draw("same");
    //    legenderrorspectrum->AddEntry(    fHistSpectrum_masterSystCorr[m], "sist corr " , "pel");
    fHistSpectrum_masterSystUnCorr[m]->GetYaxis()->SetRangeUser(0,0.35);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerStyle(20);
    fHistSpectrum_masterSystUnCorr[m]->SetLineColor(1);
    fHistSpectrum_masterSystUnCorr[m]->SetMarkerColor(3);
    fHistSpectrum_masterSystUnCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterSystUnCorr[m]->Draw("same");
    if (m==0)legenderrorspectrum->AddEntry(    fHistSpectrum_masterSystUnCorr[m], "sist uncorr " , "pel");
    fHistSpectrum_masterTotalUnCorr[m]->GetYaxis()->SetRangeUser(0,0.35);
    fHistSpectrum_masterTotalUnCorr[m]->SetMarkerStyle(20);
    fHistSpectrum_masterTotalUnCorr[m]->SetLineColor(1);
    fHistSpectrum_masterTotalUnCorr[m]->SetMarkerColor(4);
    fHistSpectrum_masterTotalUnCorr[m]->SetFillStyle(0);
    fHistSpectrum_masterTotalUnCorr[m]->Draw("same");
    if (m==0) legenderrorspectrum->AddEntry(    fHistSpectrum_masterTotalUnCorr[m], "total sist uncorr+ stat " , "pel");
    legenderrorspectrum->Draw();
  }

  //***************************+saving canvas to file *********************************
  
  fileout->WriteTObject( canvasSys[3]);
  fileout->WriteTObject( canvasSpectrumSysBis[0]);
  fileout->WriteTObject( canvasSpectrumSysBis[1]);

  for(Int_t m=0; m< nummolt +1; m++){
    fileout->WriteTObject( canvasSpectrumSys[2*m]);
    fileout->WriteTObject( canvasSpectrumSys[2*m+1]);
  }
  for(Int_t m=0; m< nummolt +1; m++){
    fileout->WriteTObject( canvasPhiSys[4*m]);
    fileout->WriteTObject( canvasPhiSys[4*m+1]);
    fileout->WriteTObject( canvasPhiSys[4*m+2]);
    fileout->WriteTObject( canvasPhiSys[4*m+3]);
  }
  fileout->Close();
  cout << "\n\n I've produced the file " << stringout << endl; 

}


