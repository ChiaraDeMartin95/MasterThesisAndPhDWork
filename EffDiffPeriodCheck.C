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
#include <TF1.h>

void ErrRatioCorr(TH1F* hNum, TH1F* hDenom, TH1F* hRatio){
  Float_t Err1=0;
  Float_t Err2=0;
  Float_t ErrC=0;
  Float_t Err=0;
  for (Int_t b=1; b<=hNum->GetNbinsX();b++){
    if (hNum->GetBinContent(b)==0 ||hDenom->GetBinContent(b)==0){
      hRatio->SetBinError(b,0);
      continue;
    }    
    Err1=pow(hNum->GetBinError(b)/hNum->GetBinContent(b),2);
    Err2=pow(hDenom->GetBinError(b)/hDenom->GetBinContent(b),2);
    ErrC=pow(hDenom->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    Err=sqrt(Err1+Err2-ErrC);
    hRatio->SetBinError(b,Err*hRatio->GetBinContent(b));  
}
  //  return hRatio;
}
void EffDiffPeriodCheck(Int_t type=1, Bool_t isEta=0, Bool_t isPhi=0, Bool_t isPt=1){

  if (isPt!=0 && (isPhi!=0 || isEta!=0)) return;
  if (isPhi && isEta) return;
 
  TString inputDir="FinalOutput/DATA";
  TString PathIn;
  TString file;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;
  const Int_t numCfr=4;

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  TString SPeriod[numPeriod]={"2016", "ALL"};
  TString SRun[numPeriod]={"2016k", "All"};
  TString Speriod[numCfr]={"2018f1_extra_DEtaEff_50runs_", "2018c12_extra", "2018j4_extra"};
 TString Syear0[numCfr]={"2016", "2016", "2016", "2016"};
 //  TString Speriod[numCfr]={"2018c12_extra", "2018j4_extra", "2018f1_extra"};
 // TString Syear0[numCfr]={"2016", "2017", "2018"};

  Int_t Marker[numCfr]={7,4,20};
  Int_t Color[numCfr]={861,601 ,628,417};

  TH1F* histoEff;
  TH1F* histoEffRatio[numCfr];
  TH1F* histoEffMaster;

  TCanvas *canvasSysV0[numSysV0];
  TLegend * legendV0 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legendV0->SetHeader("MC period");

  TString nomehisto;
  if(isPt) nomehisto="fHistV0EfficiencyPtBins_";
  if(isPhi) nomehisto="fHistV0EfficiencyPhi_";
  if(isEta) nomehisto="fHistV0EfficiencyEta_";
  /*
  PathOut="FinalOutput/histo/EffDiffPeriodCheck";
	 
  if(isMC && isEfficiency){
	 
    PathOut+="_MCEff";
  }
  if(isMC && !isEfficiency){
	 
    PathOut+="_MCTruth";
  }

  TString PathOut2=PathOut;
  TString PathOut3=PathOut;	 
  PathOut +=Form("SysT%i_SysV0%i.root",sysTrigger, sysV0); 

  PathOut2+=Form("ACBeforeCuts_SysT%i_SysV0%i.root",sysTrigger, sysV0); 
  PathOut3+=Form("ACAfterCuts_SysT%i_SysV0%i.root",sysTrigger, sysV0); 


  TFile *fileOut = new TFile(PathOut, "RECREATE"); 
  // TFile *fileOut2 = new TFile(PathOut2, "RECREATE"); 
  // TFile *fileOut3 = new TFile(PathOut3, "RECREATE"); 
  */
  TString tipo[2] = {"K0s", "Xi"};
  TString nameFinal[numCfr] = {"", "_Incl", "_Jet", "_OOJ"};
  for(Int_t i=0; i < numSysV0; i++){ 
    if (i>0) continue;
    cout << "\n\n*******sysV0 " << i << endl;
    canvasSysV0[i] = new TCanvas(Form("canvasSysV0_%i",i),Form("canvasSysV0_%i",i), 1300, 1000);
    canvasSysV0[i]->Divide(3,2);
    for(Int_t molt=0; molt< nummolt+1; molt++){
      cout << "   molt " << molt << endl;
      for(Int_t cfr=0; cfr<numCfr; cfr++){
	file=Syear0[cfr]+"/Efficiency/Efficiency"+Speriod[0];
	PathIn=inputDir + file + tipo[type]+Form("_Eta0.8_SysT%i_SysV0%i_PtMin3.0",0, i)+nameFinal[cfr]+ ".root";
	TFile* 	filein = new TFile(PathIn);
	cout <<"       cfr " << Speriod[cfr] << "  " << PathIn << endl; 
	if(cfr==0){
	  histoEffMaster=(TH1F*)filein->Get(nomehisto + Smolt[molt]);
	  if (!histoEffMaster) return;
	  histoEffMaster->Sumw2();
	}
	histoEff=(TH1F*)filein->Get(nomehisto + Smolt[molt]);
	if (!histoEff) return;
	histoEff->Sumw2();
	histoEffRatio[cfr]=(TH1F*)histoEff->Clone(Form("HistoV0EfficiencyPtBinsRatio_cfr%i", cfr));
	histoEffRatio[cfr]->Divide(histoEffMaster);
	for(Int_t j=0; j<histoEffRatio[cfr]->GetNbinsX(); j++){
	  cout << " bin j-esimo: " << j << "   errore: " <<histoEffRatio[cfr]->GetBinError(j+1)<< endl;
	}

	ErrRatioCorr(histoEff, histoEffMaster, histoEffRatio[cfr]);
	for(Int_t j=0; j<histoEffRatio[cfr]->GetNbinsX(); j++){
	  cout << " bin j-esimo: " << j << "   errore: " <<histoEffRatio[cfr]->GetBinError(j+1)<< endl;
	}
	canvasSysV0[i]->cd(molt+1);    
	histoEffRatio[cfr]->GetYaxis()->SetRangeUser(0.8,1.2);
	histoEffRatio[cfr]->SetMarkerStyle(Marker[cfr]);
	histoEffRatio[cfr]->SetLineColor(Color[cfr]);
	histoEffRatio[cfr]->SetMarkerColor(Color[cfr]);
	if (molt==0 && i==0) legendV0->AddEntry(histoEffRatio[cfr],nameFinal[cfr],"pel");
	cout << " m " << molt << " " << nameFinal[cfr] << endl;
	TF1 * pol0 = new TF1 ("pol0", "pol0", 0.5,8);
	histoEffRatio[cfr]->Fit("pol0", "R+");
	if (cfr!=0)	histoEffRatio[cfr]->Draw("samee");
	if(cfr==numCfr-1) legendV0->Draw();
      }
    }
    canvasSysV0[i]->SaveAs("EffDiffPeriodCheck.pdf");
  }

}

