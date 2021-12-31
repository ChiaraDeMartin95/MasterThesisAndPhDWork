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
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.8);
  histo->SetTitle(title);
}


void SysPurity( Int_t VariableOfInt=0, Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15,  Bool_t isMC=0,   Int_t israp=0,TString year="161718Full_AOD234_hXi"/* "1617_AOD234_hK0s"/*"1617_hK0s"/*"AllMC_hXi"/*"Run2DataRed_hXi"/,  TString Path1 ="_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, Int_t type=8,  Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=1, Int_t PtBinning=0, Bool_t isHM =1, Bool_t ispp5TeV =0){

  const Int_t numInputFiles =6;
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;                                                                        
  const Int_t numtipo=10;

  if (type==0){
    year="1617_AOD234_hK0s";
    if (isHM) year = "AllhK0sHM_RedNo16k";
    else if (ispp5TeV) year = "17pq_hK0s";
  }
  else {
    year="161718Full_AOD234_hXi";
    PtBinning=0;
    if (isHM) year="161718_HM_hXi_WithFlat16k_No18p";
    else if (ispp5TeV) year = "17pq_hXi";
  }
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100};
  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};

  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==3){
      Smolt[m] = Smoltpp5TeV[m];
      Nmolt[m] = Nmoltpp5TeV[m];
    }
  }
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
      Smolt[0] = "0-0a";
      Smolt[1] = "0-0b";
      Smolt[2] = "0-0.01";
    }
  }

  TString TypeBkg[numInputFiles];
  TString TypeRange[numInputFiles];
  Int_t FileColor[numInputFiles] = {401,801,628, 909, 881,860};
  TString SLegend[numInputFiles] = {"pol1", "pol1_Smaller", "pol1_Wider", "pol2", "pol2_Smaller", "pol2_Wider" };
  Int_t ColorMult[21]={1, 401, 801, 628, 909, 881, 860, 868, 841, 418, 628, 909, 881, 867, 921, 401, 841, 862, 866 , 865, 864};

  //input files
  TString CommonInputPath= "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution";
  if (PtBinning>0) CommonInputPath += "_PtBinning1";
  CommonInputPath += "_";
  CommonInputPath += year;
  CommonInputPath +="_"+tipo[type];
  CommonInputPath +=Srap[israp];
  CommonInputPath +=SSkipAssoc[SkipAssoc];
  CommonInputPath +="_isMeanFixedPDG";

  //histos
  TH1F * hPurity[nummolt+1][numInputFiles];
  TH1F * hPurityDenom[nummolt+1];
  TH1F * hRatioPurity[nummolt+1][numInputFiles];
  TH1F * hRelError[nummolt+1];
  TH1F * hPurityMax[nummolt+1];
  TH1F * hPurityMin[nummolt+1];

  //canvases
  TCanvas *canvasPurity[nummolt+1];
  TCanvas *canvasRelError = new TCanvas("canvasRelError","canvasRelError", 1000, 800);
  for (Int_t m=0; m<=nummolt; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    canvasPurity[m] = new TCanvas( Form("canvasPurity_m%i",m), Form("canvasPurity_m%i",m), 1600, 800);
    canvasPurity[m]->Divide(2,1);
  }

  TString Variable ="";
  TString NameHisto ="";
  if (VariableOfInt==0) {
    NameHisto = "histo_SSB";
    Variable = "Purity";
  }
  else if (VariableOfInt==1) {
    NameHisto = "histo_mean";
    Variable = "InvMass";
  }
  else if (VariableOfInt==2) {
    NameHisto = "histo_sigma";
    Variable = "Sigma";
  }
  Float_t LowRange[3] = {0.98, 0.8, 0.8};
  Float_t UpRange[3] = {1.02, 1.2, 1.2};

  //output file 
  TString stringout;
  stringout = "FinalOutput/DATA2016/SysPurity";
  stringout +=  year;
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  stringout +="_"+tipo[type];
  stringout +=Srap[israp];
  stringout +=SSkipAssoc[SkipAssoc];
  stringout+=   Form("_PtMin%.1f",  PtTrigMin);
  if (VariableOfInt>0)  stringout+=   "_"+ Variable;
  if (MultBinning!=0) stringout += Form("_MultBinning%i", MultBinning);
  TString stringoutpdf = stringout;
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  TLegend * legend = new TLegend(0.6, 0.1, 0.9, 0.2);
  TF1 * pol0 = new TF1("pol0", "pol0", 0,8);
  pol0->SetLineColor(1);
  pol0->SetParameter(0,1);

  //Comparison of purity obtained with different selections and determination of syst uncertainty
  for (Int_t m=0; m<=nummolt; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    cout << "m" << m << endl;
    for (Int_t file = 0; file < numInputFiles; file++){
      if (!isHM && !ispp5TeV && type==0 && (file == 1 || file ==2 || file ==4)) continue; 
      if (!isHM && !ispp5TeV && type==8 && (file == 3 || file == 4)) continue; 
      cout << "file " << file << endl;
      if (file<numInputFiles/2) TypeBkg[file] = "_BkgRetta";
      else  TypeBkg[file] = "_BkgParab";
      if (file==0 || file == numInputFiles/2) TypeRange[file] ="";
      else if (file<numInputFiles/2) TypeRange[file] = Form("_VarRange%i", file);
      else TypeRange[file] = Form("_VarRange%i", file-numInputFiles/2);

      TString InputPath = CommonInputPath + TypeBkg[file] + Form("_molt%i_sysT0_sysV00_Sys0_PtMin%.1f", m, PtTrigMin);
      if (MultBinning!=0) InputPath += Form("_MultBinning%i", MultBinning);
      InputPath += TypeRange[file];
      InputPath+=".root";
      cout <<"Input file: " << InputPath<< endl;
      TFile * inputFile = new TFile(InputPath, "");
      if (!inputFile) return;

      hPurity[m][file] = (TH1F*)inputFile->Get(NameHisto);
      if (!hPurity[m][file]) return;
      hPurity[m][file]->SetName(Form("histoSSB_%i", file));
      if (m==nummolt)      legend->AddEntry(hPurity[m][file], SLegend[file], "pl");
      StyleHisto(hPurity[m][file], 0.8, 1, FileColor[file], 33, "p_{T} (GeV/c)", Variable, Variable + " " + Smolt[m]+"%", 0, 0,0);

      if (file==0) {
	hPurityDenom[m] = (TH1F*) hPurity[m][file]->Clone(Form("histoSSB_Default_m%i",m));
	hPurityMax[m] = (TH1F*) hPurity[m][file]->Clone(Form("histoSSB_Max_m%i",m));
	hPurityMin[m] = (TH1F*) hPurity[m][file]->Clone(Form("histoSSB_Min_m%i",m));
      }
      else {
	hRatioPurity[m][file] = (TH1F*) hPurity[m][file]->Clone(Form("histoSSBRatio_%i", file));
	hRatioPurity[m][file]->Divide(hPurityDenom[m]);
	StyleHisto(hRatioPurity[m][file], LowRange[VariableOfInt], UpRange[VariableOfInt], FileColor[file], 33, "p_{T} (GeV/c)", "", "Ratio to default purity", 0, 0,0);
	for (Int_t b =1; b<= hPurity[m][file]->GetNbinsX(); b++){
	  cout <<  hPurity[m][file]->GetBinContent(b) << " ratio: " <<hRatioPurity[m][file]->GetBinContent(b)<< endl;
	  if (hPurity[m][file]->GetBinContent(b)>hPurityMax[m]->GetBinContent(b)) hPurityMax[m]->SetBinContent(b, hPurity[m][file]->GetBinContent(b));
	  if (hPurity[m][file]->GetBinContent(b)<hPurityMin[m]->GetBinContent(b)) hPurityMin[m]->SetBinContent(b, hPurity[m][file]->GetBinContent(b));
	}
      }

      canvasPurity[m]->cd(1);
      gStyle->SetOptStat(0);
      hPurity[m][file]->Draw("same");
      if (file==numInputFiles-1)      legend->Draw("same");

      canvasPurity[m]->cd(2);
      gStyle->SetOptStat(0);
      if (file!=0) {
	hRatioPurity[m][file]->Draw("same");
	if (file==numInputFiles-1)      legend->Draw("same");
	pol0->Draw("same");
      }
    }

    hRelError[m] = (TH1F*) hPurityMax[m]->Clone(Form("hRelError_m%i",m));
    hRelError[m]->Add(hPurityMin[m],-1);
    hRelError[m]->Scale(0.5);
    StyleHisto(hRelError[m], 0, 0.01, ColorMult[m], 33, "p_{T} (GeV/c)", "", "Relative syst. uncertainty", 0, 0,0);
    canvasRelError->cd();
    gStyle->SetOptStat(0);
    if (m==nummolt)    hRelError[m]->Draw("same");

    fileout->WriteTObject(hRelError[m]);
    fileout->WriteTObject(canvasPurity[m]);
  }
  fileout->WriteTObject(canvasRelError);

  fileout->Close();

  Int_t counterm =0;
  for (Int_t m=0; m<=nummolt; m++){
    if (isHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    counterm ++;
    cout << "m" << m << endl;

    if (counterm==1)    canvasPurity[m]->SaveAs(stringoutpdf + "_" + Variable + ".pdf(");
    else if (m==nummolt) canvasPurity[m]->SaveAs(stringoutpdf + "_" + Variable + ".pdf)");
    else    canvasPurity[m]->SaveAs(stringoutpdf + "_" + Variable + ".pdf");

    if (m==nummolt) canvasPurity[m]->SaveAs(stringoutpdf + "_" + Variable + +"_" + Smolt[m] +".pdf");
  }

  canvasRelError->SaveAs(stringoutpdf + "_RelError.pdf");
  cout << "\nHo prodotto il file: " << stringout << endl;
}
