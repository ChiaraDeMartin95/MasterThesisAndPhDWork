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
void YieldCutVariation( Int_t type=8 /*type = 0 for K0s */,Bool_t SkipAssoc=1,  Bool_t ishhCorr=0, bool isMC = 0,Bool_t isEfficiency=1, Int_t sysTrigger=0, TString year="Run2DataRed_MECorr_hXi"/*"2016kehjl_hK0s"/*"2016k_hK0s"/"17anch17_hK0s"/"1617_hK0s"/*"2018f1_extra_hK0s_CP"*/, TString year0="2016", TString Path1 ="", Int_t PtBinning=0, Int_t molt=5, Float_t PtTrigMin=3,Float_t PtTrigMax=15, Int_t rap=0, Int_t syst=0, Double_t nsigmamax=9, Bool_t isSignalFromIntegral=0, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1,   Double_t sigmacentral=4, Double_t nsigmamin=5, Int_t MultBinning=0, Int_t isSysStudy=1)
{

  //SysStudy: if 1, values different from default ones will be used; the topological variable changed will sysV0Index > 20 ) return;be ind  icated in the output file, together with its value                                                                                    
  //sysV0 indicated the type of topological variable which is varied. (trial iintervals are indicated below)                              
  //sysV0Index is the index of the array of the set of selections                                                                         

  cout << "ciao ciao " << endl;
  Int_t numSysVar=0;
  const  Int_t numSysVarK0s=6;  
  const   Int_t numSysVarXi=6;
  if (type==0) numSysVar = numSysVarK0s;
  else if (type==8) numSysVar = numSysVarXi;

  cout << "ciao ciao " << endl;
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};

  Int_t numsysV0index=20;
  if (type==8) numsysV0index=40;

  cout << "ciao ciao " << endl;
  TString SSysV0[10] = {"", "CosinePAngle", "DCAPosToPV", "DCANegToPV", "DCAV0ToPV", "ctau (cm)", "LambdaRejection (MeV/c^{2})"};
  TString SSysV0Xi[10] = {"", "DCAXiDaughters", "CosinePAngleXiToPV", "CosinePAngleV0ToXi", "InvMassLambdaWindow", "ctau (cm)", "DCAzTrigger"};

  //old  Float_t MinSysV0[10] = {0, 0.970, 0.05, 0.05, 0.2, 10, 0.003} ;
  //old  Float_t MaxSysV0[10] = {0, 0.998, 0.14, 0.14, 0.6, 40, 0.009} ;
  Float_t MinSysV0[10] = {0, 0.995, 0.05, 0.05, 0.2, 2, 0.001} ;
  Float_t MaxSysV0[10] = {0, 1, 0.14, 0.14, 0.6, 20, 0.03} ;
  Float_t MinSysV0Xi[10] = {0, 0.3, 0.95, 0.95, 0.003, 2*4.91, 0} ;
  Float_t MaxSysV0Xi[10] = {0, 1.8,    1,    1, 0.009, 5*4.91, 2} ;

  if (type==8){
    for (Int_t i=0; i<10; i++){
      MinSysV0[i] = MinSysV0Xi[i];
      MaxSysV0[i] = MaxSysV0Xi[i];
      SSysV0[i] = SSysV0Xi[i];
    }
  }
  cout << "ciao ciao " << endl;
  TCanvas *canvasYield = new TCanvas("canvasYield", "canvasYield", 1300, 800);
  canvasYield->Divide(3,2);
  TCanvas *canvasYieldRatio = new TCanvas("canvasYieldRatio", "canvasYieldRatio", 1300, 800);
  canvasYieldRatio->Divide(3,2);

  cout << "ciao ciao " << endl;
  TString nome_file_input;
  TH1F * histoYield[numSysVar];
  TH1F * histoYieldRatio[numSysVar];
  Float_t RawYield=0;
  Float_t RawYieldError=0;
  Float_t RawYieldDefault=0;
  Float_t RawYieldErrorDefault=0;
  TH1F*   histo_RawYieldDef; 
  TF1* lineat1[numSysVar];
  cout << "ciao ciao " << endl;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=9;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=10;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg","XiPos", "OmegaNeg", "OmegaPos>", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString    SSkipAssoc[2]={"_AllAssoc", ""};

  TString nome_file_output = "YieldCutVariation";
  if(isMC && isEfficiency){
    nome_file_output+="_MCEff";
  }
  nome_file_output+=Path1;
  if (PtBinning>0)    nome_file_output+=Form("_PtBinning%i",PtBinning);
  nome_file_output +=Form("_"+year+"_"+tipo[type]+Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+BkgType[isBkgParab]+"_molt%i_PtMin%.1f.root", molt, PtTrigMin);

  TFile * fileout = new TFile(nome_file_output, "RECREATE");

  cout << "ciao ciao " << endl;
  for (Int_t sysV0=0; sysV0<=numSysVar; sysV0++){
    //    if (sysV0>1) continue;
    cout << "\n\nsysV0 " << sysV0 << endl;
    lineat1[sysV0] = new TF1("pol0", "pol0", MinSysV0[sysV0], MaxSysV0[sysV0]);
    lineat1[sysV0] -> SetParameter(0,1);
    lineat1[sysV0] -> SetLineColor(kBlack);
    histoYield[sysV0] = new TH1F (Form("histoYield_%i",sysV0), SSysV0[sysV0], numsysV0index, MinSysV0[sysV0], MaxSysV0[sysV0]);
    histoYieldRatio[sysV0]= (TH1F*)       histoYield[sysV0]->Clone(Form("histoYieldRatio_%i",sysV0));
    for (Int_t sysV0index = 0; sysV0index < numsysV0index; sysV0index++){
      if (sysV0==0 && sysV0index!=0) continue;
      if (type==8 && sysV0==6 && sysV0index==0) continue;
      //if (sysV0==0 && sysV0index>3) continue;
      //      if (sysV0==6 && sysV0index>14) continue;
      //      if (sysV0index>20) continue;
      cout << " sysV0index " << sysV0index << endl;
      nome_file_input ="FinalOutput/DATA"+year0+"/invmass_distribution_thesis/invmass_distribution";
      if(isMC && isEfficiency){
	nome_file_input+="_MCEff";
      }
      nome_file_input+=Path1;
      if (PtBinning>0)    nome_file_input+=Form("_PtBinning%i",PtBinning);
      if (sysV0>0 || (sysV0==0 && sysV0index>0)){
	if (type==0){
	if (sysV0 ==1  || sysV0==5 || sysV0==6)       nome_file_input +=Form("_"+year+"_"+tipo[type]+Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_indexBis%i_Sys%i_PtMin%.1f.root", molt, sysTrigger, sysV0, sysV0index, syst,PtTrigMin);
	else       nome_file_input +=Form("_"+year+"_"+tipo[type]+Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_index%i_Sys%i_PtMin%.1f.root", molt, sysTrigger, sysV0, sysV0index, syst,PtTrigMin);
	}
	else       nome_file_input +=Form("_"+year+"_"+tipo[type]+Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_index%i_Sys%i_PtMin%.1f.root", molt, sysTrigger, sysV0, sysV0index, syst,PtTrigMin);
      }
      else      nome_file_input +=Form("_"+year+"_"+tipo[type]+Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", molt, sysTrigger, sysV0, syst,PtTrigMin);
      TFile *fin = new TFile(nome_file_input,"");
      if (!fin) return;
      cout << "il nome del file in input è " << nome_file_input << endl;
      TH1F *histo_RawYield = (TH1F*) fin->Get("histo_S");
      if (!histo_RawYield) {cout << "histogram histo_S not present " << endl; return;}
      //      if (sysV0==0) histo_RawYieldDef=(TH1F*) filein->Get("histo_S");
      //compute pT-integrated yield:
      RawYield=0;
      RawYieldError=0;
      for (Int_t b=1; b<=histo_RawYield->GetNbinsX(); b++){
	RawYield+= histo_RawYield->GetBinContent(b);
	RawYieldError+= pow(histo_RawYield->GetBinError(b),2);
      }
      RawYieldError= sqrt(RawYieldError);
      if (sysV0==0 && sysV0index==0){
	RawYieldErrorDefault = RawYieldError;
	RawYieldDefault = RawYield;
      }
      cout << "raw yield default " << RawYieldDefault << " +- " << RawYieldErrorDefault << endl;
      cout << "raw yield " << RawYield << " +- " << RawYieldError << endl;
      if (sysV0!= 0){
      histoYield[sysV0]->SetBinContent(sysV0index+1, RawYield);
      histoYieldRatio[sysV0]->SetBinContent(sysV0index+1, RawYield/RawYieldDefault);
      histoYieldRatio[sysV0]->GetYaxis()->SetRangeUser(0.8,1.2);
      }
    } //end of loop on sysV0index

    cout << "\nraw yield default " << RawYieldDefault << " +- " << RawYieldErrorDefault << endl;
    if (sysV0!=0){
    canvasYield->cd(sysV0);
    histoYield[sysV0]->GetXaxis()-> SetTitle(SSysV0[sysV0]);
    histoYield[sysV0]->Draw("");
    canvasYieldRatio->cd(sysV0);
    histoYieldRatio[sysV0]->GetXaxis()-> SetTitle(SSysV0[sysV0]);
    //    histoYieldRatio[sysV0]->GetYaxis()-> SetTitle("Yield ratio to default");
    histoYieldRatio[sysV0]->Draw("");
    lineat1[sysV0]->Draw("same");
    }
  }

  fileout->WriteTObject(canvasYield);
  fileout->WriteTObject(canvasYieldRatio); 
  fileout->Close();
  cout << " ho creato il file " << nome_file_output << endl;
}
