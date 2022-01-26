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
#include "TGraphAsymmErrors.h"
#include <TRandom.h>
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include </data/dataalice/cdemart/AliPhysicsChiara/AliPhysics/PWGLF/SPECTRA/UTILS/YieldMean.C>
//#include </data/dataalice/cdemart/ALICE_analysis_tutorial/YieldMean.C>

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

void PtSpectraInclusive(Int_t type=0, Bool_t isMC=0, TString year="", Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year0="2016",  Int_t MultBinning=1, Int_t PtBinning=1, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,  Bool_t TwoFitFunctions=0, Bool_t isNormCorrFullyComputed=0, Bool_t isMeanMacro=0, Bool_t isErrorAssumedPtCorr=1, Bool_t isSyst = 0, Bool_t isBkgParab =0){

  if (type!=0 && type!=8) {cout << "Macro not working for particles different from K0s (type=0) and Xi (type=8)" << endl; return;}


  //******* PUBLISHED YIELD VS MULT (MB) *********
  TFile *filedatipubbl;
  TDirectoryFile *dir;
  filedatipubbl= new TFile("HEPData-1574358449-v1-Table_8c.root", "");
  if (!filedatipubbl) return;
  //  cout << " I got the files where published data are stored " << endl;
  dir  = (TDirectoryFile*)filedatipubbl->Get("Table 8c");
  if (!dir) return;

  TH1F * fHistYieldDatiPubblicati;
  TH1F * fHistYieldDatiPubblicatiStat;
  TH1F * fHistYieldDatiPubblicatiSist;
  TH1F * fHistYieldDatiPubblicatiSistUncorr;
  TH1F * fHistYieldErroriDatiPubblicatiStat;
  TH1F * fHistYieldErroriDatiPubblicatiSist;
  TH1F * fHistYieldErroriDatiPubblicatiSistUncorr;
  TString yParticle = "y1"; //K0s
  if (type==8) yParticle = "y3"; //Xi
  fHistYieldDatiPubblicati=(TH1F*)dir->Get("Hist1D_"+yParticle);
  fHistYieldErroriDatiPubblicatiStat=(TH1F*)dir->Get("Hist1D_"+yParticle+"_e1");  //this is the statistical uncertainty                                
  fHistYieldErroriDatiPubblicatiSistUncorr=(TH1F*)dir->Get("Hist1D_"+yParticle+"_e3"); //this is only the uncorr systematic error                       
  fHistYieldErroriDatiPubblicatiSist=(TH1F*)dir->Get("Hist1D_"+yParticle+"_e2"); //this is the total systematic error (the statistic is 1/10 and the uncorrelated systematic is approx half)
  fHistYieldDatiPubblicatiStat=(TH1F*)  fHistYieldDatiPubblicati->Clone("Hist1D_Stat_"+yParticle);
  fHistYieldDatiPubblicatiSistUncorr=(TH1F*)  fHistYieldDatiPubblicati->Clone("Hist1D_SistUncorr_"+yParticle);
  fHistYieldDatiPubblicatiSist=(TH1F*)  fHistYieldDatiPubblicati->Clone("Hist1D_Sist_"+yParticle);
  for(Int_t k=1; k <= fHistYieldDatiPubblicati->GetNbinsX(); k++){
    fHistYieldDatiPubblicatiStat->SetBinError(k,     fHistYieldErroriDatiPubblicatiStat->GetBinContent(k));
    fHistYieldDatiPubblicatiSist->SetBinError(k,     fHistYieldErroriDatiPubblicatiSist->GetBinContent(k));
    fHistYieldDatiPubblicatiSistUncorr->SetBinError(k,     fHistYieldErroriDatiPubblicatiSistUncorr->GetBinContent(k));
  }

  //******* PUBLISHED <pT> VS MULT (MB) *********
  TFile *filedatipubblPt;
  TDirectoryFile *dirPt;
  filedatipubblPt= new TFile("HEPData-ins1748157-v1-Table_5c.root", "");
  if (!filedatipubblPt) return;
  //  cout << " I got the files where published data are stored " << endl;
  dirPt  = (TDirectoryFile*)filedatipubblPt->Get("Table 5c");
  if (!dirPt) return;

  TH1F * fHistAvgPtDatiPubblicati;
  TH1F * fHistAvgPtDatiPubblicatiStat;
  TH1F * fHistAvgPtDatiPubblicatiSist;
  TH1F * fHistAvgPtDatiPubblicatiSistUncorr;
  TH1F * fHistAvgPtErroriDatiPubblicatiStat;
  TH1F * fHistAvgPtErroriDatiPubblicatiSist;
  TH1F * fHistAvgPtErroriDatiPubblicatiSistUncorr;
  fHistAvgPtDatiPubblicati=(TH1F*)dirPt->Get("Hist1D_"+yParticle);
  fHistAvgPtErroriDatiPubblicatiStat=(TH1F*)dirPt->Get("Hist1D_"+yParticle+"_e1");  //this is the statistical uncertainty                                
  fHistAvgPtErroriDatiPubblicatiSistUncorr=(TH1F*)dirPt->Get("Hist1D_"+yParticle+"_e3"); //this is only the uncorr systematic error                       
  fHistAvgPtErroriDatiPubblicatiSist=(TH1F*)dirPt->Get("Hist1D_"+yParticle+"_e2"); //this is the total systematic error (the statistic is 1/10 and the uncorrelated systematic is approx half)
  fHistAvgPtDatiPubblicatiStat=(TH1F*)  fHistAvgPtDatiPubblicati->Clone("Hist1D_Stat_"+yParticle);
  fHistAvgPtDatiPubblicatiSistUncorr=(TH1F*)  fHistAvgPtDatiPubblicati->Clone("Hist1D_SistUncorr_"+yParticle);
  fHistAvgPtDatiPubblicatiSist=(TH1F*)  fHistAvgPtDatiPubblicati->Clone("Hist1D_Sist_"+yParticle);
  for(Int_t k=1; k <= fHistAvgPtDatiPubblicati->GetNbinsX(); k++){
    fHistAvgPtDatiPubblicatiStat->SetBinError(k,     fHistAvgPtErroriDatiPubblicatiStat->GetBinContent(k));
    fHistAvgPtDatiPubblicatiSist->SetBinError(k,     fHistAvgPtErroriDatiPubblicatiSist->GetBinContent(k));
    fHistAvgPtDatiPubblicatiSistUncorr->SetBinError(k,     fHistAvgPtErroriDatiPubblicatiSistUncorr->GetBinContent(k));
  }

  //files from where norm factor is taken
  TString SfileNormCorrFC ="";

  Bool_t  SkipAssoc =0;
  Float_t  PtTrigMin=0.;
  Float_t  PtTrigMax = 30;
  Int_t   israp =1;
  Bool_t  isppHM =1;
  MultBinning=1;
  if (type==0) {
    year = "161718_HM_hK0s_INELgt0";
    //    yearMC = "LHC19h11aPlus_hK0s_INELgt0";
  }
  PtBinning = 0;

  TString RegionType[4] = {"Jet", "Bulk", "Inclusive", "Bulk"};
  TString RegionTypeOld[4] = {"Jet", "Bulk", "All", "Bulk"};
  TString RegionTypeNew[4] = {"Jet", "Bulk", "Full", "Bulk"};

  gStyle->SetOptStat(0);

  TF1 * lineat2= new TF1 ("lineat2", "pol0", 0,8);
  lineat2->FixParameter(0,2);
  lineat2->SetLineColor(kBlack);
  lineat2->SetLineWidth(1);
  TF1 * lineatm2= new TF1 ("lineatm2", "pol0", 0,8);
  lineatm2->FixParameter(0,-2);
  lineatm2->SetLineColor(1);
  lineatm2->SetLineWidth(1);
  TF1 * lineat0= new TF1 ("lineat0", "pol0", 0,8);
  lineat0->FixParameter(0,0);
  lineat0->SetLineColor(kBlack);
  lineat0->SetLineWidth(1);

  TString hhCorr[2]={"", "_hhCorr"};
  TFile *fileinbis;
  TFile *filein;

  TString PathIn1;
  TString file = year;
  if (PtBinning!=0)  file+=Form("_PtBinning%i", PtBinning);

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=20;
  const Int_t numPtTrigger=1;
  const Int_t numsysDPhi=4;
  const Int_t numtipo=10;

  Float_t LowPtLimitForAvgPtFS[nummolt+1] = {0};

  Int_t PtV0Min = 0; //0 el
  if (type>0)   PtV0Min =1;
  if (type==0 && PtBinning!=0) PtV0Min=1;

  TString DataOrMC[2]={"Data", "MC"};
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  Float_t SpectrumSup[numtipo] = {0.2, 0.01};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 

  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  Double_t Nmoltpp5TeV[nummolt+1]={0, 10, 100, 100, 100, 100};
  Double_t NmoltHM[nummolt+1]={0, 0, 0, 0.01, 0.05, 0.1}; 
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  TString SmoltHM[nummolt+1]={"0-0a", "0-0b", "0-0.01", "0.01-0.05", "0.05-0.1", "0-0.1"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegendHM[nummolt+1]={"0-0a %", "0-0b %", "0-0.01 %", "0.01-0.05 %", "0.05-0.1 %", "0-0.1 %"};
  TString SmoltLegendpp5TeV[nummolt+1]={"0-10 %", "10-100 %", "100-100 %", "100-100 %", "100-100 %", "0-100 %"};

  for (Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Nmolt[m] = Nmolt0[m];
      Smolt[m] = Smolt0[m];
      SmoltLegend[m] = SmoltLegend0[m];
    }
    else     if (MultBinning==1){
      Nmolt[m] = Nmolt1[m];
      Smolt[m] = Smolt1[m];
      SmoltLegend[m] = SmoltLegend1[m];
    }
    else     if (MultBinning==2){
      Nmolt[m] = Nmolt2[m];
      Smolt[m] = Smolt2[m];
      SmoltLegend[m] = SmoltLegend2[m];
    }
    else if (MultBinning==3){
      Nmolt[m] = Nmoltpp5TeV[m];
      Smolt[m] = Smoltpp5TeV[m];
      SmoltLegend[m] = SmoltLegendpp5TeV[m];
    }
    if (isppHM){
      Nmolt[m] = NmoltHM[m];
      Smolt[m] = SmoltHM[m];
      SmoltLegend[m] = SmoltLegendHM[m];
    }
  }

  TString SPtV0[numPtV0]={"0.1-0.3", "0.3-0.5", "0.5-0.7", "0.7-0.9", "0.9-1.1", "1.1-1.3", "1.3-1.5", "1.5-1.7", "1.7-1.9", "1.9-2.1", "2.1-2.3", "2.3-2.5", "2.5-2.7", "2.7-3.0", "3.0-3.4", "3.2-3.8", "3.8-4.2", "4.2-4.8", "4.8-6.0", "6.0-8.0"};
  Double_t NPtV0[numPtV0+1]={0.1,0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.4, 3.8, 4.2, 4.8, 6.0, 8};

  Int_t numPtV0Max=numPtV0;

  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[4] ={628, 418, 600, 418};
  TString BkgType[2]={"BkgRetta", "BkgParab"};

  TString stringout;
  TString PathIn0;
  TString PathIn;
  PathIn0 = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_" + year;
  if(type>=0){
    PathIn0 +="_"+tipo[type];
    PathIn0 +=Srap[israp];
    PathIn0 +=SSkipAssoc[SkipAssoc];
  }
  PathIn0 += "_isMeanFixedPDG_";
  PathIn0 += BkgType[isBkgParab];

  stringout = Dir+"/DATA"+year0+"/";
  stringout += "PtSpectraInclusive" +hhCorr[0];
  if (isppHM)   stringout += "_pp13TeVHM";
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  //  if (year!="1617_hK0s" && year!="Run2DataRed_MECorr_hXi") stringout+= "_"+year;
  stringout+= "_"+year;
  if(type>=0){
    stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_PtMin%.1f", PtTrigMin);
  if (TwoFitFunctions) stringout += "_TwoFitFunctions";
  if (isNormCorrFullyComputed) stringout+="_isNormCorrFullyComputed";
  //stringout += "_Fixed";
  if (isMeanMacro)  stringout += "_YieldMeanMacro"; 
  if (isErrorAssumedPtCorr)  stringout += "_isErrorAssumedPtCorr"; 
  TString PathOutPictures = stringout;
  if (MultBinning!=0) stringout += Form("_MultBinning%i", MultBinning);
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  cout << "\nDefinition integral regions... " << endl;
  //low and upper ranges for the fit****************************
  //************************************************************

  Double_t LowRange[nummolt+1]= {0.1, 0.1, 0.1, 0.1, 0.1, 0.1}; 
  Double_t UpRange[nummolt+1]= {8,8,8,8,8,8}; 

  Double_t LowRangeSpectrumPart[nummolt+1]= {0};
  Double_t UpRangeSpectrumPart[nummolt+1]= {0};

  Float_t LowPtLimitForAvgPtFSAllMult = LowPtLimitForAvgPtFS[nummolt];
  //end of low and upper ranges for the fit****************************
  //************************************************************

  Int_t PtBinMin[nummolt+1]={0};
  for (Int_t m=0; m<nummolt+1; m++){
    //defining region where bin counting is performed (LowRangeSpectrumPart < pt < UpRangeSpectrumPart)
    //LowRangeSpectrumPart:
    LowRangeSpectrumPart[m] = LowRange[m];

    for(Int_t v=0; v < numPtV0Max; v++){
      if (LowRangeSpectrumPart[m]<=NPtV0[v]) break;
      PtBinMin[m]++;
      //      cout << "LowRange " << LowRangeSpectrumPart[m] << " v " << v << " PtBin " << NPtV0[v]<< endl;
      //      cout << "PtBinMin " << PtBinMin[m] << endl;
    }
    UpRangeSpectrumPart[m] = 8;
  }

  for (Int_t m=0; m<nummolt+1; m++){
    cout << "\n\e[32m*******Multiplicity: " << SmoltLegend[m] << " *************\e[39m" << endl;
    cout << "fit range " << LowRange[m] << "- " << UpRange[m] << endl;
    cout << "!extrapolation range " << LowRangeSpectrumPart[m] << "- " << UpRangeSpectrumPart[m] << endl;
  }

  if(isMC && isEfficiency) file = year + "_MCEff";
  //  file+=Form("_PtBinning%i", PtBinning);

  TString titleX=  "p_{T} (GeV/c)";
  TString titleY=  "1/#Delta#eta #Delta#phi 1/N_{trigg} dN/dp_{T}";
  titleY = "";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 

  TLegend *legendPhi = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legendError = new TLegend(0.6, 0.7, 0.9, 0.9);
  TLegend *legendError2 = new TLegend(0.6, 0.7, 0.9, 0.9);

  TH1F* fHistSpectrumStat[nummolt+1];
  TH1F* fHistSpectrumStatNotNorm[nummolt+1];
  TH1F* fHistSpectrumStatUp[nummolt+1];
  TH1F* fHistSpectrumStatDown[nummolt+1];
  TH1F* fHistSpectrumStatHard[nummolt+1];
  TH1F* fHistSpectrumStatSoft[nummolt+1];
  TH1F* fHistSpectrumTemp[nummolt+1];
  TH1F* fHistSpectrumStatpol0[nummolt+1];
  TH1F* fHistSpectrumStatOOJSubDef[nummolt+1];
  TH1F* fHistSpectrumStatLeadTrackDef[nummolt+1];
  TH1F* fHistSpectrumStatMCChoiceDef[nummolt+1];
  TH1F* fHistSpectrumStatFakeSBDef[nummolt+1];
  TH1F* fHistSpectrumStatEtaEff[nummolt+1];
  TH1F* fHistSpectrumStatEtaEffRatio[nummolt+1];
  TH1F* fHistSpectrumStatEtaEffRatioRef[nummolt+1];
  TH1F* fHistSpectrumStatNormCorr[nummolt+1];
  TH1F* fHistSpectrumStatNormCorrRatio[nummolt+1];
  TH1F* fHistSpectrumStatNormCorrRatioRef[nummolt+1];
  TH1F* fHistSpectrumSist[nummolt+1];
  TH1F* fHistSpectrum_max[nummolt+1];
  TH1F* fHistSpectrum_min[nummolt+1];

  TH1F* fHistSpectrumStatDPhi[nummolt+1][numsysDPhi];
  TH1F* fHistSpectrumSistDPhi[nummolt+1];
  TH1F* fHistSpectrumSistOOJ[nummolt+1];
  TH1F* fHistSpectrumSistAll[nummolt+1];
  TH1F* fHistSpectrumSistLeadTrack[nummolt+1];
  TH1F* fHistSpectrumSistMCChoice[nummolt+1];
  TH1F* fHistSpectrumSistFakeSB[nummolt+1];
  TH1F* fHistSpectrumSistOOJSubDef[nummolt+1];

  //histos for relative uncertainty
  TH1F* fHistSpectrumStatRelError[nummolt+1]; //stat
  TH1F* fHistSpectrumSistRelError[nummolt+1]; //sist on the DeltaPhi projections
  TH1F* fHistSpectrumSistRelErrorDPhi[nummolt+1]; //sist assoc to choice of DeltaPhi interval
  TH1F* fHistSpectrumSistRelErrorAll[nummolt+1]; //total sist
  TH1F* fHistSpectrumSistRelErrorDCAz[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorPurity[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorSE[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorDeltaEta[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorOOBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorIBPU[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorLeadTrack[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorMCChoice[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorFakeSB[nummolt+1]; 
  TH1F* fHistSpectrumSistRelErrorOOJSubDef[nummolt+1]; 

  TCanvas* canvasYield = new TCanvas ("canvasYield", "canvasYield", 1300, 800);
  TCanvas* canvasYieldNormCorr = new TCanvas ("canvasYieldNormCorr", "canvasYieldNormCorr", 1300, 800);
  TCanvas* canvasYieldNormCorrRatio = new TCanvas ("canvasYieldNormCorrRatio", "canvasYieldNormCorrRatio", 1300, 800);
  TCanvas* canvasPtvsMult = new TCanvas ("canvasPtvsMult", "canvasPtvsMult", 1300, 800);
  TCanvas* canvasYieldErr = new TCanvas ("canvasYieldErr", "canvasYieldErr", 1300, 800);
  TCanvas* canvasPtErr = new TCanvas ("canvasPtErr", "canvasPtErr", 1300, 800);
  TCanvas* canvasNormFactor = new TCanvas ("canvasNormFactor", "canvasNormFactor", 1300, 800);
  if (isppHM || MultBinning==3)   canvasNormFactor->Divide(2,2);
  else   canvasNormFactor->Divide(3,2);
  TCanvas* canvasPtSpectra = new TCanvas ("canvasPtSpectra", "canvasPtSpectra", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectra->Divide(2,2);
  else   canvasPtSpectra->Divide(3,2);
  TCanvas* canvaspol0 = new TCanvas ("canvaspol0", "canvaspol0", 1300, 800);
  if (isppHM || MultBinning==3)   canvaspol0->Divide(2,2);
  else canvaspol0->Divide(3,2);
  TCanvas* canvasOOJSubDef = new TCanvas ("canvasOOJSubDef", "canvasOOJSubDef", 1300, 800);
  canvasOOJSubDef->Divide(3,2);
  TCanvas* canvasLeadTrackDef = new TCanvas ("canvasLeadTrackDef", "canvasLeadTrackDef", 1300, 800);
  canvasLeadTrackDef->Divide(3,2);
  TCanvas* canvasMCChoiceDef = new TCanvas ("canvasMCChoiceDef", "canvasMCChoiceDef", 1300, 800);
  canvasMCChoiceDef->Divide(3,2);
  TCanvas* canvasFakeSBDef = new TCanvas ("canvasFakeSBDef", "canvasFakeSBDef", 1300, 800);
  canvasFakeSBDef->Divide(3,2);

  TCanvas* canvasNormCorr = new TCanvas ("canvasNormCorr", "canvasNormCorr", 1300, 800);
  canvasNormCorr->Divide(3,2);
  TCanvas* canvasNormCorrRatio = new TCanvas ("canvasNormCorrRatio", "canvasNormCorrRatio", 1300, 800);
  canvasNormCorrRatio->Divide(3,2);

  TCanvas* canvasPtSpectraFit = new TCanvas ("canvasPtSpectraFit", "canvasPtSpectraFit", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraFit->Divide(2,2);
  else   canvasPtSpectraFit->Divide(3,2);
  TCanvas* canvasPtSpectraFitUp = new TCanvas ("canvasPtSpectraFitUp", "canvasPtSpectraFitUp", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraFitUp->Divide(2,2);
  else  canvasPtSpectraFitUp->Divide(3,2);
  TCanvas* canvasPtSpectraFitDown = new TCanvas ("canvasPtSpectraFitDown", "canvasPtSpectraFitDown", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraFitDown->Divide(2,2);
  else  canvasPtSpectraFitDown->Divide(3,2);
  TCanvas* canvasPtSpectraFitHard = new TCanvas ("canvasPtSpectraFitHard", "canvasPtSpectraFitHard", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraFitHard->Divide(2,2);
  else  canvasPtSpectraFitHard->Divide(3,2);
  TCanvas* canvasPtSpectraFitSoft = new TCanvas ("canvasPtSpectraFitSoft", "canvasPtSpectraFitSoft", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraFitSoft->Divide(2,2);
  else  canvasPtSpectraFitSoft->Divide(3,2);
  TCanvas* canvasDummy = new TCanvas ("canvasDummy", "canvasDummy", 1300, 800);
  TCanvas* canvasPtSpectraFitRatio = new TCanvas ("canvasPtSpectraFitRatio", "canvasPtSpectraFitRatio", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraFitRatio->Divide(2,2);
  else  canvasPtSpectraFitRatio->Divide(3,2);

  TCanvas* canvasPtSpectraFitBis = new TCanvas ("canvasPtSpectraFitBis", "canvasPtSpectraFitBis", 1300, 800);
  canvasPtSpectraFitBis->Divide(3,2);

  TCanvas* canvasPtSpectraAll = new TCanvas ("canvasPtSpectraAll", "canvasPtSpectraAll", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraAll->Divide(2,2);
  else   canvasPtSpectraAll->Divide(3,2);
  TCanvas* canvasBarlow = new TCanvas ("canvasBarlow", "canvasBarlow", 1300, 800);
  if (isppHM || MultBinning==3)   canvasBarlow->Divide(2,2);
  else  canvasBarlow->Divide(3,2);
  TCanvas* canvasBarlowpol0 = new TCanvas ("canvasBarlowpol0", "canvasBarlowpol0", 1300, 800);
  canvasBarlowpol0->Divide(3,2);
  TCanvas* canvasBarlowOOJSubDef = new TCanvas ("canvasBarlowOOJSubDef", "canvasBarlowOOJSubDef", 1300, 800);
  canvasBarlowOOJSubDef->Divide(3,2);
  TCanvas* canvasBarlowLeadTrackDef = new TCanvas ("canvasBarlowLeadTrackDef", "canvasBarlowLeadTrackDef", 1300, 800);
  canvasBarlowLeadTrackDef->Divide(3,2);
  TCanvas* canvasBarlowMCChoiceDef = new TCanvas ("canvasBarlowMCChoiceDef", "canvasBarlowMCChoiceDef", 1300, 800);
  canvasBarlowMCChoiceDef->Divide(3,2);
  TCanvas* canvasBarlowFakeSBDef = new TCanvas ("canvasBarlowFakeSBDef", "canvasBarlowFakeSBDef", 1300, 800);
  canvasBarlowFakeSBDef->Divide(3,2);


  TCanvas* canvasPtSpectraRelError = new TCanvas ("canvasPtSpectraRelError", "canvasPtSpectraRelError", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraRelError->Divide(2,2);  
  else  canvasPtSpectraRelError->Divide(3,2);
  TCanvas* canvasPtSpectraRelErrorAll = new TCanvas ("canvasPtSpectraRelErrorAll", "canvasPtSpectraRelErrorAll", 1300, 800);
  if (isppHM || MultBinning==3)   canvasPtSpectraRelErrorAll->Divide(2,2);  
  else  canvasPtSpectraRelErrorAll->Divide(3,2);

  Float_t HighPtRatio[numtipo] = {1.15, 0, 0, 0, 0, 0, 0, 0, 1.2, 0};
  Float_t LowPtRatio[numtipo] = {0.85, 0, 0, 0, 0, 0, 0, 0, 0.8, 0};

  cout << "\nprendo histo per confronto con dati pubblicati " << endl;                                        
  TString PathDatiPubblicati ="";                                                                               
  if (type==0) PathDatiPubblicati = "HEPData-ins1748157-v1-Table_1.root";                                       
  else if (type==8) PathDatiPubblicati = "HEPData-1583750454-v1-Table_3.root";                                  
  TFile *filedatipubblicati = new TFile(PathDatiPubblicati, "");                                                
  if (!filedatipubblicati) {cout << "file dati pubblicati not there " << endl; return;}                            TString STable = "Table 3";                                                                                      if (type==0) STable = "Table 1";                                                                              
  TDirectoryFile *dirspectra = (TDirectoryFile*)filedatipubblicati->Get(STable);                                 
  if (!dirspectra)  {cout << "directory dati pubblicati not there " << endl; return;}              
  TH1F *   hspectrum[11];
  TH1F *   hspectrum1[11];
  TH1F *   hspectrum2[11];
  TH1F *   hspectrum3[11];
  for (Int_t i=0; i<11; i++){
    hspectrum[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i", i+1));                                            
    hspectrum1[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e1", i+1)); //stat
    hspectrum2[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e2", i+1)); //syst total                           
    hspectrum3[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e3", i+1)); //sist Uncorr                           
    if (!hspectrum[i] ||     !hspectrum1[i] || !hspectrum2[i]|| !hspectrum3[i] ) { cout << "histo is missing " << endl; return;}                                                                                              
  }


  TH1F*  fHistSpectrumSistRelErrorPublished = (TH1F*) hspectrum2[10]->Clone("fHistSpectrumSistRelErrorPublished");
  fHistSpectrumSistRelErrorPublished->Divide(hspectrum[10]);
  for (Int_t b=0; b<= fHistSpectrumSistRelErrorPublished->GetNbinsX(); b++){
    fHistSpectrumSistRelErrorPublished->SetBinError(b,0);
  }

  Bool_t BarlowPassed[nummolt+1][numsysDPhi]={0};
  Int_t BarlowSign[nummolt+1][numsysDPhi]={0};
  TH1F* hBarlowVar[nummolt+1][numsysDPhi];
  Float_t BarlowVar[nummolt+1][numPtV0][numsysDPhi]={0};
  Bool_t BarlowPassedpol0[nummolt+1]={0};
  Int_t BarlowSignpol0[nummolt+1]={0};
  TH1F* hBarlowVarpol0[nummolt+1];
  Float_t BarlowVarpol0[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedOOJSubDef[nummolt+1]={0};
  Int_t BarlowSignOOJSubDef[nummolt+1]={0};
  TH1F* hBarlowVarOOJSubDef[nummolt+1];
  Float_t BarlowVarOOJSubDef[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedLeadTrackDef[nummolt+1]={0};
  Int_t BarlowSignLeadTrackDef[nummolt+1]={0};
  TH1F* hBarlowVarLeadTrackDef[nummolt+1];
  Float_t BarlowVarLeadTrackDef[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedMCChoiceDef[nummolt+1]={0};
  Int_t BarlowSignMCChoiceDef[nummolt+1]={0};
  TH1F* hBarlowVarMCChoiceDef[nummolt+1];
  Float_t BarlowVarMCChoiceDef[nummolt+1][numPtV0]={0};
  Bool_t BarlowPassedFakeSBDef[nummolt+1]={0};
  Int_t BarlowSignFakeSBDef[nummolt+1]={0};
  TH1F* hBarlowVarFakeSBDef[nummolt+1];
  Float_t BarlowVarFakeSBDef[nummolt+1][numPtV0]={0};

  Float_t ErrDPhi[nummolt+1][numPtV0]={0};
  Int_t ErrDPhiCounter[nummolt+1]={0};

  Int_t ColorsysDPhi[12] = {807,  881, 867,909, 881,860,868,841,  418, 881, 7, 1};

  Float_t LimSupYield=3;
  Float_t LimInfYield=0.1;

  Float_t LimInfPtvsMult=10e-5;
  Float_t LimSupPtvsMult=1.8;
  Float_t LimInfPtvsMultTight=10e-5;
  Float_t LimSupPtvsMultTight=2;
  Int_t nrebin[nummolt+1] ={8, 8, 8, 8, 8, 8};
  Float_t LimSupYieldErr=0.04;
  Float_t LimSupPtErr=0.04;

  Float_t LimSup=2;
  Float_t LimInf=0;

  Float_t LimSupNormComp=1.2;
  Float_t LimInfNormComp=0.8;

  Float_t LimSupError=0.01;
  Float_t LimInfError=0;
  Float_t LimSupErrorLog=0.01;

  TH1F*  hDeltaPhiLimit;
  Float_t  LowBinDPhi[numsysDPhi] ={0};
  Float_t  UpBinDPhi[numsysDPhi]={0};

  //first part: I get default spectra with stat uncertainty
  cout << "\n**************************"<<endl;
  cout << "First part: I get default spectra with stat uncertainty";
  TString  PathInDef="";
  for(Int_t m=0; m<nummolt+1; m++){
    if (isppHM && m<2) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    PathInDef=PathIn0;
    PathInDef += Form("_molt%i_sysT0_sysV00_Sys0_PtMin%.1f", m, PtTrigMin);
    if (MultBinning!=0) PathInDef += Form("_MultBinning%i", MultBinning);
    PathInDef += "_INEL.root";
    cout << " from the file: " << PathInDef << endl;
    TFile *  fileinDef = new TFile(PathInDef, "");
    if (!fileinDef) {cout << PathInDef << " does not exist" << endl; return;}

    //    cout << Smolt[m] << endl;
    fHistSpectrumStat[m]=(TH1F*)fileinDef->Get("histo_SEffCorr");
    if (!fHistSpectrumStat[m]) {cout << "fHistSpectrumStat does not exist" << endl; return;}
    fHistSpectrumStat[m]->SetName("histo_SEffCorr" + Smolt[m]);
    fHistSpectrumStatUp[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumUp_"+Smolt[m]);
    fHistSpectrumStatDown[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumDown_"+Smolt[m]);
    fHistSpectrumStatHard[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumHard_"+Smolt[m]);
    fHistSpectrumStatSoft[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSoft_"+Smolt[m]);
    fHistSpectrumTemp[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumTemp_"+Smolt[m]);

    if (isSyst){
      fHistSpectrumSistRelErrorDCAz[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorDCAz_"+Smolt[m]); 
      fHistSpectrumSistRelErrorPurity[m]    =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorPurity_"+Smolt[m]); 
      fHistSpectrumSistRelErrorSE[m]      =(TH1F*)fileinDef->Get("fHistSpectrumSistRelErrorSE_"+Smolt[m]);
      fHistSpectrumSist[m]=(TH1F*)fileinDef->Get("fHistSpectrumSist_"+Smolt[m]);
      if (!fHistSpectrumSistRelErrorDCAz[m]) {cout << "Rel syst uncertainty DCAz " << endl; return;}
      if (!fHistSpectrumSistRelErrorPurity[m]) {cout << "Rel syst uncertainty purity " << endl; return;}
      if (!fHistSpectrumSistRelErrorSE[m]) {cout << "Rel syst uncertainty topological selections " << endl; return;}
      if (!fHistSpectrumSist[m]) {cout << "fHistSpectrumSist does not exist" << endl; return;}
    }
    else       fHistSpectrumSist[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSist_"+Smolt[m]);

    if (isSyst)    fHistSpectrumSistRelErrorAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistRelErrorAll_"+Smolt[m]);
    else     fHistSpectrumSistRelErrorAll[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorAll_"+Smolt[m]);
    fHistSpectrumStatRelError[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumStatRelError_"+Smolt[m]);
    fHistSpectrumSistRelErrorMB[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMB_"+Smolt[m]);
    fHistSpectrumSistRelErrorOOBPU[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorOOBPU_"+Smolt[m]);
    fHistSpectrumSistRelErrorIBPU[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorIBPU_"+Smolt[m]);
    fHistSpectrumSistRelErrorMCChoice[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistRelErrorMCChoice_"+Smolt[m]);

    fHistSpectrum_max[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMax_"+Smolt[m]);
    fHistSpectrum_min[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumMin_"+Smolt[m]);
    if (isSyst){
      fHistSpectrumSistAll[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistAll_"+Smolt[m]);
      fHistSpectrumSistMCChoice[m]=(TH1F*)fHistSpectrumSist[m]->Clone("fHistSpectrumSistMCChoice_"+Smolt[m]);
    }
    else {
      fHistSpectrumSistAll[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistAll_"+Smolt[m]);
      fHistSpectrumSistMCChoice[m]=(TH1F*)fHistSpectrumStat[m]->Clone("fHistSpectrumSistMCChoice_"+Smolt[m]);
    }

    if (isSyst){
      for(Int_t v=PtV0Min; v <    fHistSpectrumSistAll[m]->GetNbinsX() ; v++){
	//      cout << "v " << v << endl;
	if (fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSist[m]->GetBinError(v+1)/ fHistSpectrumSist[m]->GetBinContent(v+1)); //this works for inclusive, for jet and OOJ I will change it
	  fHistSpectrumStatRelError[m]->SetBinContent(v+1, fHistSpectrumStat[m]->GetBinError(v+1)/ fHistSpectrumStat[m]->GetBinContent(v+1));
	}
	else       if (type==0 && NPtV0[v] == 0.1){ //?
	  fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
	  //      cout << "sist rel error " <<    fHistSpectrumSistRelErrorAll[m]->GetBinContent(v+1) << endl;
	  fHistSpectrumStatRelError[m]->SetBinError(v+1,0);
	}
	else{
	  fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
	  fHistSpectrumStatRelError[m]->SetBinError(v+1,0);
	}
      }
    }
  } //end loop m

  //Rel uncertainties taken from Fiorella's analysis
  Float_t  Sigma2OOBPU[nummolt+1] ={0};
  Float_t  Sigma2IBPU[nummolt+1] ={0};
  Float_t  Sigma2MB[nummolt+1]={0};   
  Float_t OOBPU = 0.012; //relative syst. uncertainty on pT spectrum associated to OOB PileUp taken from Fiorella's paper and "averaged" over pT
  Float_t IBPU=0.02; //same but for Inbunch pileup
  Float_t MB = 0.01; //same but for material budget
  if (type==8){
    OOBPU = IBPU = MB = 0.02;
  }


  cout << "\n*******************************************"<<endl;
  cout<< "Systematic effect associated to MC used to calculate efficiency  " << endl;
  Float_t   YieldSpectrumErrMCChoice[nummolt+1]={0};
  TString    PathInMCChoiceDef;
  TLegend * legendMCChoice= new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t IsSignMCChoice=0;
  Int_t Varm = 0;
  TString PathMCChoiceSyst="";

  for(Int_t m=0; m<nummolt+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      //      cout << "v " << SPtV0[v] << endl;
      fHistSpectrumSistMCChoice[m]->SetBinError(v+1,0);
    }
  }

  if (isSyst){
    //Summing errors altogether********************
    for(Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	Sigma2OOBPU[m] =pow(OOBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
	Sigma2IBPU[m] =pow(IBPU *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;
	Sigma2MB[m]= pow(MB *fHistSpectrumStat[m]->GetBinContent(v+1),2) ;

	fHistSpectrumSistAll[m]->SetBinError(v+1,sqrt(pow(fHistSpectrumSist[m]->GetBinError(v+1),2) + Sigma2OOBPU[m]+ Sigma2IBPU[m] + Sigma2MB[m] + pow(fHistSpectrumSistMCChoice[m]->GetBinError(v+1),2))); 

	if (fHistSpectrumSistAll[m]->GetBinContent(v+1) == 0) fHistSpectrumSistAll[m]->SetBinError(v+1, 0);
	if ( fHistSpectrumStat[m]->GetBinContent(v+1)!=0){
	  fHistSpectrumSistRelErrorAll[m]->SetBinContent(v+1, fHistSpectrumSistAll[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	  fHistSpectrumSistRelErrorMB[m]->SetBinContent(v+1,MB);
	  fHistSpectrumSistRelErrorIBPU[m]->SetBinContent(v+1,IBPU);
	  fHistSpectrumSistRelErrorOOBPU[m]->SetBinContent(v+1,OOBPU);
	  fHistSpectrumSistRelErrorMCChoice[m]->SetBinContent(v+1,fHistSpectrumSistMCChoice[m]->GetBinError(v+1)/fHistSpectrumStat[m]->GetBinContent(v+1));
	}

	fHistSpectrumSistRelErrorAll[m]->SetBinError(v+1,0);
	fHistSpectrumSistRelErrorMB[m]->SetBinError(v+1,0);
	fHistSpectrumSistRelErrorIBPU[m]->SetBinError(v+1,0);
	fHistSpectrumSistRelErrorOOBPU[m]->SetBinError(v+1,0);
	fHistSpectrumSistRelErrorMCChoice[m]->SetBinError(v+1,0);
      }
  
      //  cout << "Drawing sist rel error " << endl;
      if (isppHM)     canvasPtSpectraRelErrorAll->cd(m+1-2);
      else     canvasPtSpectraRelErrorAll->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError, Color[type], 33, titleX, titleYRel ,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError, Color[type], 27, titleX, titleYRel ,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorSE[m], LimInfError, LimSupError, 807, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorDCAz[m], LimInfError, LimSupError, 867, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorPurity[m], LimInfError, LimSupError, 834, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorMB[m], LimInfError, LimSupError, 401, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorOOBPU[m], LimInfError, LimSupError, 825, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorIBPU[m], LimInfError, LimSupError, 631, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      StyleHisto(fHistSpectrumSistRelErrorMCChoice[m], LimInfError, LimSupError, 907, 27, titleX, titleYRel,  title+SmoltLegend[m],0,0,0);
      gPad->SetLogy();

      TLegend *legendErrorAll= new TLegend(0.6, 0.1, 0.9, 0.4);
      legendErrorAll->AddEntry(fHistSpectrumStatRelError[m], "stat.", "pl");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorSE[m], "syst. signal extr.", "l");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorDCAz[m], "syst. DCAzTrigger", "l");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorPurity[m], "syst. purity", "l");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMB[m], "Material Budget", "l");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorOOBPU[m], "Out-of-bunch PU", "l");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorIBPU[m], "In-bunch PU", "l");
      legendErrorAll->AddEntry(fHistSpectrumSistRelErrorMCChoice[m], "MCChoice", "l");

      fHistSpectrumSistRelErrorAll[m]->Draw("same p");
      fHistSpectrumSistRelErrorMB[m]->Draw("same");
      fHistSpectrumSistRelErrorIBPU[m]->Draw("same");
      fHistSpectrumSistRelErrorOOBPU[m]->Draw("same");
      fHistSpectrumSistRelErrorSE[m]->Draw("same");
      fHistSpectrumSistRelErrorDCAz[m]->Draw("same");
      fHistSpectrumSistRelErrorPurity[m]->Draw("same");
      fHistSpectrumSistRelErrorMCChoice[m]->Draw("same");
      fHistSpectrumStatRelError[m]->Draw("same p");
      legendErrorAll->Draw("");
    }
  }

  cout << "\n //*************spectra normalization ******************" << endl;

  /*
  TFile* fileNormCorrFC;
  TH1F * fHistNormCorrFC[nummolt+1];
  TH1F * fHistNormCorrAllMultFC[nummolt+1];
  TF1 * pol0NormFactor[nummolt+1];
  Float_t pol0RelError=0;
  Float_t SpectrumStatValueBeforeScaling[nummolt+1][numPtV0] = {0};
  Float_t SpectrumStatErrorBeforeScaling[nummolt+1][numPtV0] = {0};

  if (isNormCorrFullyComputed==1){ // In this case I perform a correction with the *comnplete* normalisation factor.
    cout << "Normalization factor for comparison with event generators ---- fully computed!" << endl;
    cout << "From file: " << SfileNormCorrFC << endl;
    fileNormCorrFC=new TFile(SfileNormCorrFC,"");
    if (!fileNormCorrFC) {cout << "file norm not there " << endl; return;}
    for(Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;

      fHistSpectrumStatNotNorm[m] = (TH1F*) fHistSpectrumStat[m]->Clone(""+ Smolt[m] + "_NotNorm");
      if (isppHM)     canvasNormFactor->cd(m+1-2);
      else if (MultBinning==3 && ispp5TeV && m==nummolt)  canvasNormFactor->cd(3);
      else     canvasNormFactor->cd(m+1);

      if (TypeAnalysis==3) RegionTypeOld[TypeAnalysis] = "BulkBlue";

      if (isppHM)      fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",0));
      else if (ispp5TeV) {
	if (type==8 || (type==0 && TypeAnalysis!=2)) fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",5));
	else 	fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));
      }
      else fHistNormCorrFC[m] = (TH1F*)fileNormCorrFC->Get("SpectrumRatio"+RegionTypeOld[TypeAnalysis]+Form("_m%i",m));

      if (TypeAnalysis==3) RegionTypeOld[TypeAnalysis] = "Bulk";
      pol0NormFactor[m] = new TF1(Form("pol0NormFactor_m%i", m), "pol0", 0, 8);
      if (TypeAnalysis == 0 && type==8) {
	cout << "pol0 fit to normalization factor (only for Xi in jet)" << endl;
	fHistNormCorrFC[m]-> Fit(pol0NormFactor[m], "R+");
      }

      fHistNormCorrFC[m]-> Draw("");
      if (TypeAnalysis == 0 && type==8)  {
	for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ 
	  SpectrumStatValueBeforeScaling[m][b-1] = 	fHistSpectrumStat[m]->GetBinContent(b);
	  SpectrumStatErrorBeforeScaling[m][b-1] = 	fHistSpectrumStat[m]->GetBinError(b);
	}
 	fHistSpectrumStat[m]->Scale(1./pol0NormFactor[m]->GetParameter(0));
	for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ 
	  pol0RelError = pol0NormFactor[m]->GetParError(0)/pol0NormFactor[m]->GetParameter(0);
	  fHistSpectrumStat[m]->SetBinError(b, fHistSpectrumStat[m]->GetBinContent(b) * sqrt(pow(SpectrumStatErrorBeforeScaling[m][b-1]/SpectrumStatValueBeforeScaling[m][b-1], 2) + pow(pol0RelError, 2)));
	  cout << "pol0Rel error " << pol0RelError << endl;
	  cout << fHistSpectrumStat[m]->GetBinError(b) / fHistSpectrumStat[m]->GetBinContent(b) << endl;
	}
      }
      else 	fHistSpectrumStat[m]->Divide(fHistNormCorrFC[m]);
      
      for (Int_t b=PtBinMin[m]+1; b<= fHistSpectrumSistAll[m]->GetNbinsX(); b++){ //I only have to *scale* sist error
	//	cout <<       fHistNormCorrFC[m]->GetBinContent(b) << endl;
	if (TypeAnalysis == 0 && type==8)  {
	  fHistSpectrumSistAll[m]->SetBinContent(b,fHistSpectrumSistAll[m]->GetBinContent(b)/pol0NormFactor[m]->GetParameter(0));
	  fHistSpectrumSistAll[m]->SetBinError(b,fHistSpectrumSistAll[m]->GetBinError(b)/pol0NormFactor[m]->GetParameter(0));
	}
	else {
	  fHistSpectrumSistAll[m]->SetBinContent(b,fHistSpectrumSistAll[m]->GetBinContent(b)/fHistNormCorrFC[m]->GetBinContent(b));
	  fHistSpectrumSistAll[m]->SetBinError(b,fHistSpectrumSistAll[m]->GetBinError(b)/fHistNormCorrFC[m]->GetBinContent(b));
	  if (fHistNormCorrFC[m]->GetBinContent(b) ==0 ) {
	    fHistSpectrumSistAll[m]->SetBinContent(b,0);
	    fHistSpectrumSistAll[m]->SetBinError(b,0);
	  }
	}
      }
    }
  }

  cout << "...end of spectra normalization *********************** \n" << endl;
  */
  for(Int_t m=0; m<nummolt+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;

    if (isppHM)     canvasPtSpectra->cd(m+1-2);
    else     canvasPtSpectra->cd(m+1);
    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistAll[m], LimInf, LimSup, Color[type], 1, titleX, titleY, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumStat[m], LimInf, LimSup, Color[type], 1, titleX, titleY,  title+SmoltLegend[m], 0, 0, 0);
    fHistSpectrumStat[m]->DrawClone("same");
    fHistSpectrumSistAll[m]->SetFillStyle(0);
    fHistSpectrumSistAll[m]->DrawClone("same e2");
    if (m==nummolt){
      legendError2->AddEntry(fHistSpectrumStat[m], "stat.", "ple");
      legendError2->AddEntry(fHistSpectrumSist[m], "syst.", "fe");
    }
    legendError2->Draw("");
    cout << "\n\e[32mMultiplicity: " << SmoltLegend[m]  << "\e[39m" << endl;
    for (Int_t b= PtBinMin[m] ; b< fHistSpectrumStat[m]->GetNbinsX(); b++){
      cout << " " << NPtV0[b] << " < pt < " << NPtV0[b+1] << " GeV/c , Yield: " << fHistSpectrumStat[m]->GetBinContent(b+1) << ", rel. stat. error: " << fHistSpectrumStat[m]->GetBinError(b+1)/ fHistSpectrumStat[m]->GetBinContent(b+1) << endl;
      cout << " " << NPtV0[b] << "< pt < " << NPtV0[b+1] << " GeV/c, Yield: " << fHistSpectrumSistAll[m]->GetBinContent(b+1) << ", rel. sist. error: " << fHistSpectrumSistAll[m]->GetBinError(b+1)/ fHistSpectrumSistAll[m]->GetBinContent(b+1) << endl;
    }

    if (isppHM)     canvasPtSpectraRelErrorAll->cd(m+1-2);
    else     canvasPtSpectraRelErrorAll->cd(m+1);

    gPad->SetLeftMargin(0.15);
    StyleHisto(fHistSpectrumSistRelErrorAll[m], LimInfError, LimSupError,Color[type] , 27, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumStatRelError[m], LimInfError, LimSupError,Color[type] , 33, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);
    StyleHisto(fHistSpectrumSistRelErrorPublished, LimInfError, LimSupError, 922, 33, titleX, titleYRel, title+SmoltLegend[m], 0, 0, 0);

    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");

    if (isppHM)     canvasPtSpectraRelError->cd(m+1-2);
    else     canvasPtSpectraRelError->cd(m+1);
    gPad->SetLeftMargin(0.15);
    if (m==nummolt) fHistSpectrumSistRelErrorPublished->Draw("same");
    fHistSpectrumSistRelErrorAll[m]->Draw("same p");
    fHistSpectrumStatRelError[m]->Draw("same p");
    if (m==nummolt){
      legendError->AddEntry(    fHistSpectrumStatRelError[m], "stat.", "pl");
      legendError->AddEntry(    fHistSpectrumSistRelErrorAll[m], "syst.", "pl");
    }
    if (m==nummolt)      legendError->AddEntry(fHistSpectrumSistRelErrorPublished, "syst. published", "l");
    legendError->Draw("");
  } //end loop m

  //fourth part: fit to obtain pt-integrated yield vs mult
  AliPWGFunc pwgfunc;
  const Int_t numfittipo=4;
  Int_t numfittipoEff =0;
  TLegend *legendfit=new TLegend(0.6, 0.6, 0.9, 0.9);
  TString   nameMTscaling[nummolt+1][numfittipo];

  TF1 * rettaUno= new TF1("rettaUno", "pol0",0,8);
  rettaUno->FixParameter(0,1);

  TH1F* fHistSpectrumRatioFit[nummolt+1][numfittipo];
  TH1F* hhoutYield[nummolt+1];
  TH1F* hhoutYieldMy[nummolt+1];
  TH1F* hhoutYieldRatioToMine[nummolt+1];
  TH1F* hhoutRelStat[nummolt+1];
  TH1F* hhoutRelStatMy[nummolt+1];
  TH1F* hhoutRelSystHigh[nummolt+1];
  TH1F* hhoutRelSystHighMy[nummolt+1];
  TH1F* hhoutRelSystLow[nummolt+1];
  TH1F* hhoutRelSystLowMy[nummolt+1];

  TH1F* hhoutAvgPt[nummolt+1];
  TH1F* hhoutAvgPtMy[nummolt+1];
  TH1F* hhoutAvgPtRatioToMine[nummolt+1];
  TH1F* hhoutAvgPtRelStat[nummolt+1];
  TH1F* hhoutAvgPtRelStatMy[nummolt+1];
  TH1F* hhoutAvgPtRelSystHigh[nummolt+1];
  TH1F* hhoutAvgPtRelSystHighMy[nummolt+1];
  TH1F* hhoutAvgPtRelSystLow[nummolt+1];
  TH1F* hhoutAvgPtRelSystLowMy[nummolt+1];

  Float_t hhoutYieldAvg[nummolt+1]={0};
  Float_t hhoutRelStatAvg[nummolt+1]={0};
  Float_t hhoutRelSystHighAvg[nummolt+1]={0};
  Float_t hhoutRelSystLowAvg[nummolt+1]={0};

  Float_t hhoutAvgPtAvg[nummolt+1]={0};
  Float_t hhoutAvgPtRelStatAvg[nummolt+1]={0};
  Float_t hhoutAvgPtRelSystHighAvg[nummolt+1]={0};
  Float_t hhoutAvgPtRelSystLowAvg[nummolt+1]={0};

  TH1 *hhout[numfittipo][nummolt+1];
  TString  Titlehhout[9] = {"kYield",
    "kYieldStat",
    "kYieldSysHi",
    "kYieldSysLo",
    "kMean",
    "kMeanStat",
    "kMeanSysHi",
    "kMeanSysLo",
    "kExtra"};

  Int_t ColorFit[numfittipo+1]={860, 881, 868, 628, 419};
  TFitResultPtr fFitResultPtr0[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrUp[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrDown[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrHard[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrSoft[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtr1[nummolt+1][numfittipo];
  TFitResultPtr fFitResultPtrTemp[nummolt+1][numfittipo];
  TF1* fit_MTscaling[nummolt+1][numfittipo];
  TF1* fit_MTscalingUp[nummolt+1][numfittipo];
  TF1* fit_MTscalingDown[nummolt+1][numfittipo];
  TF1* fit_MTscalingHard[nummolt+1][numfittipo];
  TF1* fit_MTscalingSoft[nummolt+1][numfittipo];
  TF1* fit_MTscalingBis[nummolt+1][numfittipo];
  TF1* fit_MTscalingTemp[nummolt+1][numfittipo];
  TString       nameFit[numfittipo]={"mT-scaling", "Boltzmann", "Fermi-Dirac", "Levi"};
  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
  Int_t factor=1;

  Float_t   AvgPt[nummolt+1]={0};
  Float_t   AvgPtFS[nummolt+1]={0};
  Float_t   AvgPtMeanMacro[nummolt+1]={0};

  Float_t   AvgPtDiscrErr[nummolt+1] = {0};
  Float_t   AvgPtFSNum[nummolt+1]={0};
  Float_t   AvgPtFSDenom[nummolt+1]={0};
  Float_t   AvgPtFSTemp[nummolt+1]={0};
  Float_t   AvgPtFSNumTemp[nummolt+1]={0};
  Float_t   AvgPtFSDenomTemp[nummolt+1]={0};
  Float_t   AvgPtFSHard[nummolt+1]={0};
  Float_t   AvgPtFSNumHard[nummolt+1]={0};
  Float_t   AvgPtFSDenomHard[nummolt+1]={0};
  Float_t   AvgPtFSSoft[nummolt+1]={0};
  Float_t   AvgPtFSNumSoft[nummolt+1]={0};
  Float_t   AvgPtFSDenomSoft[nummolt+1]={0};
  Float_t   AvgPtHard[nummolt+1]={0};
  Float_t   AvgPtSoft[nummolt+1]={0};
  Float_t   AvgPtMax[nummolt+1]={0};
  Float_t   AvgPtMin[nummolt+1]={0};
  Float_t   AvgTemp[3] = {0};
  Float_t   AvgPtFit[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitMinFit[nummolt+1]={0};
  Float_t   AvgPtFitMaxFit[nummolt+1]={0};
  Float_t   AvgPtMeanMacroMinFit[nummolt+1]={0};
  Float_t   AvgPtMeanMacroMaxFit[nummolt+1]={0};
  Float_t   AvgPtFitUp[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitDown[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitHard[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitSoft[nummolt+1][numfittipo]={0};
  Float_t   AvgPtFitTemp[nummolt+1][numfittipo]={0};

  Float_t   AvgPtSistErr[nummolt+1]={0};
  Float_t   AvgPtStatErr[nummolt+1]={0};
  Float_t   AvgPtSistErrFS[nummolt+1]={0};
  Float_t   AvgPtStatErrFS[nummolt+1]={0};
  Float_t   AvgPtSistErrMeanMacro[nummolt+1]={0};
  Float_t   AvgPtStatErrMeanMacro[nummolt+1]={0};
  Float_t   AvgPtSistErrFit[nummolt+1]={0};
  Float_t   AvgPtSistErrFitMeanMacro[nummolt+1]={0};
  Float_t   AvgPtStatErrFit[nummolt+1][numfittipo]={0};

  Float_t   YieldExtrHighPt[nummolt+1][numfittipo]={0};
  Float_t   YieldExtrLowPt[nummolt+1][numfittipo]={0};

  Float_t    YieldExtrHighPtAvg[nummolt+1]={0};
  Float_t    YieldExtrHighPtAvgUp[nummolt+1]={0};
  Float_t    YieldExtrHighPtAvgDown[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvg[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvgUp[nummolt+1]={0};
  Float_t    YieldExtrLowPtAvgDown[nummolt+1]={0};
  Float_t    YieldErrStatHighPtAvg[nummolt+1]={0};
  Float_t    YieldErrStatLowPtAvg[nummolt+1]={0};
  Float_t    YieldExtr[nummolt+1]={0};
  Float_t    YieldExtrMaxLowPt[nummolt+1]={0};
  Float_t    YieldExtrMinLowPt[nummolt+1]={0};
  Float_t    YieldExtrMaxHighPt[nummolt+1]={0};
  Float_t    YieldExtrMinHighPt[nummolt+1]={0};
  Float_t    YieldExtrErrStat[nummolt+1]={0};
  Float_t    YieldExtrErrSist[nummolt+1]={0};
  Float_t    YieldExtrErrSistFourFit[nummolt+1]={0};
  Float_t    Yield4FitErrSistLowPt[nummolt+1]={0};
  Float_t    Yield4FitErrSistHighPt[nummolt+1]={0};
  Float_t    YieldExtrErrSistLow[nummolt+1]={0};
  Float_t    YieldErrSystExtrLowPt[nummolt+1]={0};
  Float_t    YieldErrSystExtrHighPt[nummolt+1]={0};
  Float_t    YieldSpectrum[nummolt+1]={0};
  Float_t    YieldSpectrumErrStat[nummolt+1]={0};
  Float_t    YieldSpectrumErrSist[nummolt+1]={0};
  Float_t    Yield[nummolt+1]={0};
  Float_t    YieldErrStat[nummolt+1]={0};
  Float_t    YieldErrSist[nummolt+1]={0};
  Float_t    YieldErrSistMy[nummolt+1]={0};
  Float_t    YieldErrSistLow[nummolt+1]={0};
  Float_t    YieldErrRelNormFactor[nummolt+1]={0, 0, 0, 0.015, 0.138, 0.08};

  TCanvas* canvasFitResult = new TCanvas ("canvasFitResult", "canvasFitResult", 1300, 800);
  TH1F * hFitResult[numfittipo];
  TH1F*  fHistAvgPtDistr[nummolt+1][numfittipo];
  TH1F*  fHistAvgPtFSDistr[nummolt+1];

  cout << "\n\e[35m*******************************************************\n" << endl;
  cout << "Fitting the pt spectra... " <<endl; 
  cout << "\n*******************************************************\n\e[39m" << endl;
  for (Int_t typefit =0; typefit<numfittipo; typefit++){
    hFitResult[typefit] = new TH1F ("hFitResult"+nameFit[typefit], "hFitResult",nummolt, Nmolt);
  }
  for(Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && MultBinning==1 && m<=1) continue;
    cout << "\n\e[35m*** Multiplicity: " << SmoltLegend[m] << " ***\e[39m" << endl;
    hhoutYield[m] = new TH1F(Form("hhoutYield_m%i", m),Form("hhoutYield_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutYieldMy[m] = new TH1F(Form("hhoutYieldMy_m%i", m),Form("hhoutYieldMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelStat[m] = new TH1F(Form("hhoutRelStat_m%i", m),Form("hhoutRelStat_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelStatMy[m] = new TH1F(Form("hhoutRelStatMy_m%i", m),Form("hhoutRelStatMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystHigh[m] = new TH1F(Form("hhoutRelSystHigh_m%i", m),Form("hhoutRelSystHigh_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystHighMy[m] = new TH1F(Form("hhoutRelSystHighMy_m%i", m),Form("hhoutRelSystHighMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystLow[m] = new TH1F(Form("hhoutRelSystLow_m%i", m),Form("hhoutRelSystLow_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutRelSystLowMy[m] = new TH1F(Form("hhoutRelSystLowMy_m%i", m),Form("hhoutRelSystLowMy_m%i", m) , numfittipo+1, 0, numfittipo+1);

    hhoutAvgPt[m] = new TH1F(Form("hhoutAvgPt_m%i", m),Form("hhoutAvgPt_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtMy[m] = new TH1F(Form("hhoutAvgPtMy_m%i", m),Form("hhoutAvgPtMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelStat[m] = new TH1F(Form("hhoutAvgPtRelStat_m%i", m),Form("hhoutAvgPtRelStat_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelStatMy[m] = new TH1F(Form("hhoutAvgPtRelStatMy_m%i", m),Form("hhoutAvgPtRelStatMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystHigh[m] = new TH1F(Form("hhoutAvgPtRelSystHigh_m%i", m),Form("hhoutAvgPtRelSystHigh_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystHighMy[m] = new TH1F(Form("hhoutAvgPtRelSystHighMy_m%i", m),Form("hhoutAvgPtRelSystHighMy_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystLow[m] = new TH1F(Form("hhoutAvgPtRelSystLow_m%i", m),Form("hhoutAvgPtRelSystLow_m%i", m) , numfittipo+1, 0, numfittipo+1);
    hhoutAvgPtRelSystLowMy[m] = new TH1F(Form("hhoutAvgPtRelSystLowMy_m%i", m),Form("hhoutAvgPtRelSystLowMy_m%i", m) , numfittipo+1, 0, numfittipo+1);

    fHistAvgPtFSDistr[m] = new TH1F (Form("fHistAvgPtFSDistr_m%i", m), Form("fHistAvgPtFSDistr_m%i", m), 10000, LimInfPtvsMult, LimSupPtvsMult);

    AvgPtMeanMacroMaxFit[m] = 0;
    AvgPtMeanMacroMinFit[m] = 1000;
    numfittipoEff=0;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (nameFit[typefit]!="Levi") continue;
      numfittipoEff ++;
      pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
      nameMTscaling[m][typefit] = Form("fitMTscaling_m%i_fit%i",m, typefit);
      cout << "\n***********Fitting pt spectra with: " << nameFit[typefit]<< " (Quiet mode selected) "<< endl;
      if (typefit==0)      fit_MTscaling[m][typefit]=    pwgfunc.GetMTExp(massParticle[type], 0.1, 0.04*factor, nameMTscaling[m][typefit]); //mass, T, norm, name                                
      if (typefit==1) {
	fit_MTscaling[m][typefit]=    pwgfunc.GetBoltzmann(massParticle[type],0.1, 0.04*factor, nameMTscaling[m][typefit]);
      }
      if (typefit==2)      fit_MTscaling[m][typefit]=    pwgfunc.GetFermiDirac(massParticle[type],0.1, 0.04*factor, nameMTscaling[m][typefit]);
      if (typefit==3)  {
	fit_MTscaling[m][typefit]=    pwgfunc.GetLevi(massParticle[type],0.1, 0.03, 0.04*factor, nameMTscaling[m][typefit]);  //norm, n, T, mass (but the function must be called with these parameters in inverse order) 
	fit_MTscaling[m][typefit]->SetParLimits(0, 0,fHistSpectrumStat[m]->GetBinContent(fHistSpectrumStat[m]->GetMaximumBin())*0.5*10); //norm
	fit_MTscaling[m][typefit]->SetParLimits(2, 0.1, 10); //T
	fit_MTscaling[m][typefit]->SetParLimits(1, 2, 30); //n
	if (type==8)         fit_MTscaling[m][typefit]->SetParLimits(1, 15,30);
	fit_MTscaling[m][typefit]->SetParameter(2, 0.7);
      }
      if (m==nummolt)    legendfit->AddEntry( fit_MTscaling[m][typefit], nameFit[typefit],"l");

      fit_MTscaling[m][typefit]->SetLineColor(ColorFit[typefit]);

      fit_MTscalingBis[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Bis");

      fit_MTscaling[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fit_MTscalingBis[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      if (TwoFitFunctions) {
	fit_MTscaling[m][typefit]->SetRange(LowRange[m], 1.5);
	fit_MTscalingBis[m][typefit]->SetRange(2, UpRange[m]);
      }

      if (isMeanMacro){
	hhout[typefit][m] = YieldMean(fHistSpectrumStat[m], fHistSpectrumSistAll[m],  fit_MTscalingBis[m][typefit], 0, 20, 0.01, 0.1, "0qI", "log.root", LowRange[m], UpRange[m]);
	hhout[typefit][m] -> SetLineColor(ColorFit[typefit]);
	hhout[typefit][m] -> SetName("hhout_"+ nameFit[typefit]+Form("_m%i", m));
	hhout[typefit][m]->GetYaxis()->SetRangeUser(0,2);
	for (Int_t b=1; b<= hhout[typefit][m]->GetNbinsX(); b++){
	  hhout[typefit][m] ->GetXaxis() ->SetBinLabel(b,Titlehhout[b-1] );
	}
	cout << "m " << m << " typefit " << typefit << endl;
	hhoutYield[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutRelStat[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutRelSystHigh[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutRelSystLow[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

	hhoutAvgPt[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutAvgPtRelStat[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutAvgPtRelSystHigh[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
	hhoutAvgPtRelSystLow[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

	cout << "m " << m << " typefit " << typefit << endl;
	hhoutYield[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(9)); //extrapolatedfraction
	hhoutYieldAvg[m] +=       hhoutYield[m] ->GetBinContent(typefit+1); 

	hhoutRelStat[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(2)); // stat error
	hhoutRelStatAvg[m] +=       hhoutRelStat[m] ->GetBinContent(typefit+1); 

	hhoutRelSystHigh[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(3)); // syst error
	hhoutRelSystHighAvg[m] +=       hhoutRelSystHigh[m] ->GetBinContent(typefit+1); 

	hhoutRelSystLow[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(4)); // syst error
	hhoutRelSystLowAvg[m] +=       hhoutRelSystLow[m] ->GetBinContent(typefit+1); 

	hhoutAvgPt[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(5));  //avg pt
	hhoutAvgPtAvg[m] +=       hhoutAvgPt[m] ->GetBinContent(typefit+1); 

	hhoutAvgPtRelStat[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(6)); // stat error of avg pt
	hhoutAvgPtRelStatAvg[m] +=       hhoutAvgPtRelStat[m] ->GetBinContent(typefit+1); 

	hhoutAvgPtRelSystHigh[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(7)); // syst error of avg pt
	hhoutAvgPtRelSystHighAvg[m] +=       hhoutAvgPtRelSystHigh[m] ->GetBinContent(typefit+1); 

	hhoutAvgPtRelSystLow[m] ->SetBinContent(typefit+1, hhout[typefit][m]->GetBinContent(8)); // syst error of avg pt
	hhoutAvgPtRelSystLowAvg[m] +=       hhoutAvgPtRelSystLow[m] ->GetBinContent(typefit+1); 

	//2.C SYSTEMATIC RELATED TO CHOICE OF FIT FUNCTION
	if(AvgPtMeanMacroMaxFit[m] < hhout[typefit][m]->GetBinContent(5)){
	  AvgPtMeanMacroMaxFit[m] = hhout[typefit][m]->GetBinContent(5);
	}
	if(AvgPtMeanMacroMinFit[m] > hhout[typefit][m]->GetBinContent(5)){
	  AvgPtMeanMacroMinFit[m] = hhout[typefit][m]->GetBinContent(5);
	}
	cout <<  AvgPtMeanMacroMaxFit[m] << " " <<   AvgPtMeanMacroMinFit[m]<< endl;

      }

      fFitResultPtr0[m][typefit] = fHistSpectrumStat[m]->Fit(fit_MTscaling[m][typefit],"SR0IQ");
      fit_MTscaling[m][typefit]->SetRange(0,20);

      if (isppHM)     canvasPtSpectraFit->cd(m+1-2);
      else     canvasPtSpectraFit->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscaling[m][typefit]->Draw("same");
      //      legendfit->Draw("");

      fFitResultPtr1[m][typefit]=       fHistSpectrumStat[m]->Fit(    fit_MTscalingBis[m][typefit],"SR0IQ");
      fit_MTscalingBis[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitBis->cd(m+1-2);
      else     canvasPtSpectraFitBis->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStat[m]->Draw("same");
      fit_MTscalingBis[m][typefit]->Draw("same");
      legendfit->Draw("");

      fHistSpectrumRatioFit[m][typefit] = (TH1F*) fHistSpectrumStat[m]->Clone(Form("HistRatioFit_m%i_typefit%i", m, typefit));
      fHistSpectrumRatioFit[m][typefit]->Divide(fit_MTscaling[m][typefit]);
      if (isppHM)     canvasPtSpectraFitRatio->cd(m+1-2);
      else     canvasPtSpectraFitRatio->cd(m+1);
      gPad->SetLeftMargin(0.15);
      StyleHisto(fHistSpectrumRatioFit[m][typefit], 0.4, 2, ColorFit[typefit], 33, "p_{T} (GeV/c)", "", "Spectra/Fit ratio", 0, 0, 0);
      rettaUno->SetLineColor(kBlack);
      fHistSpectrumRatioFit[m][typefit]->Draw("samee");
      if (typefit==0)      rettaUno->Draw("same");
      legendfit->Draw("");

      hFitResult[typefit] ->SetBinContent(m+1,fit_MTscaling[m][typefit]->GetChisquare()/fit_MTscaling[m][typefit]->GetNDF());

      cout << "DEFINE + FIT SPECTRA TO COMPUTE SYST UNCERTAINTY" << endl;
      //DEFINE + FIT SPECTRA TO COMPUTE SYST UNCERTAINTY
      for(Int_t v=PtV0Min; v < numPtV0Max; v++){
	fHistSpectrumStatUp[m]->SetBinContent(v+1,fHistSpectrumSistAll[m]->GetBinError(v+1)+fHistSpectrumStat[m]->GetBinContent(v+1));
	fHistSpectrumStatDown[m]->SetBinContent(v+1,-fHistSpectrumSistAll[m]->GetBinError(v+1)+fHistSpectrumStat[m]->GetBinContent(v+1));
      }

      fHistSpectrumStatHard[m] = (TH1F*)YieldMean_ReturnExtremeHardHisto(fHistSpectrumSist[m]);
      fHistSpectrumStatHard[m]->SetName("fHistSpectrumHard_"+Smolt[m]);
      fHistSpectrumStatSoft[m] = (TH1F*)YieldMean_ReturnExtremeSoftHisto(fHistSpectrumSist[m]);
      fHistSpectrumStatSoft[m]->SetName("fHistSpectrumSoft_"+Smolt[m]);

      //fit +1sigma sistematica
      fit_MTscalingUp[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Up");    
      fit_MTscalingUp[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrUp[m][typefit]=       fHistSpectrumStatUp[m]->Fit(    fit_MTscalingUp[m][typefit],"SR0IQ");
      fit_MTscalingUp[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitUp->cd(m+1-2);
      else     canvasPtSpectraFitUp->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatUp[m]->Draw("same");
      fit_MTscalingUp[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit -1sigma sistematica
      fit_MTscalingDown[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Down");
      fit_MTscalingDown[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrDown[m][typefit]=       fHistSpectrumStatDown[m]->Fit(    fit_MTscalingDown[m][typefit],"SR0IQ");
      fit_MTscalingDown[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitDown->cd(m+1-2);
      else     canvasPtSpectraFitDown->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatDown[m]->Draw("same");
      fit_MTscalingDown[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit of softened spectrum
      fit_MTscalingSoft[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Soft");
      fit_MTscalingSoft[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrSoft[m][typefit]=       fHistSpectrumStatSoft[m]->Fit(    fit_MTscalingSoft[m][typefit],"SR0IQ");
      fit_MTscalingSoft[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitSoft->cd(m+1-2);
      else     canvasPtSpectraFitSoft->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatSoft[m]->Draw("same");
      fit_MTscalingSoft[m][typefit]->Draw("same");
      legendfit->Draw("");

      //fit of hardened spectrum
      fit_MTscalingHard[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Hard");
      fit_MTscalingHard[m][typefit]->SetRange(LowRange[m], UpRange[m]);
      fFitResultPtrHard[m][typefit]=       fHistSpectrumStatHard[m]->Fit(    fit_MTscalingHard[m][typefit],"SR0IQ");
      fit_MTscalingHard[m][typefit]->SetRange(0,20);
      if (isppHM)     canvasPtSpectraFitHard->cd(m+1-2);
      else     canvasPtSpectraFitHard->cd(m+1);
      gPad->SetLeftMargin(0.15);
      fHistSpectrumStatHard[m]->Draw("same");
      fit_MTscalingHard[m][typefit]->Draw("same");
      //      legendfit->Draw("");

      canvasDummy->cd();

      cout << "CALCULATE YIELDS"<< endl;
      //1A. CALCULATE YIELDS
      YieldExtrLowPt[m][typefit]= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPt[m][typefit] = fit_MTscalingBis[m][typefit]->Integral(UpRangeSpectrumPart[m],20);
      YieldExtrLowPtAvg[m]+= fit_MTscaling[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvg[m]+= fit_MTscalingBis[m][typefit]->Integral(UpRangeSpectrumPart[m],20);     

      hhoutYieldMy[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutRelStatMy[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutRelSystHighMy[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutRelSystLowMy[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

      hhoutAvgPtMy[m]       ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutAvgPtRelStatMy[m]     ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutAvgPtRelSystHighMy[m] ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);
      hhoutAvgPtRelSystLowMy[m]  ->GetXaxis()->SetBinLabel(typefit+1, nameFit[typefit]);

      hhoutYieldMy[m] ->SetBinContent(typefit+1, YieldExtrLowPt[m][typefit] + YieldExtrHighPt[m][typefit]);
      //      hhoutRelStatMy[m] ->SetBinContent(typefit+1, sqrt(pow(    fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ), 2) + pow ( fit_MTscalingBis[m][typefit]->IntegralError(UpRangeSpectrumPart[m],20,fFitResultPtr1[m][typefit] ->GetParams(),(fFitResultPtr1[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() ),2)) ); //this is only the stat uncertainty related to the extrapolated part!
      hhoutRelStatMy[m] ->SetBinContent(typefit+1, 0);
      hhoutRelSystHighMy[m] ->SetBinContent(typefit+1, 0);
      hhoutRelSystLowMy[m] ->SetBinContent(typefit+1, 0);

      hhoutAvgPtRelStatMy[m] ->SetBinContent(typefit+1, 0);
      hhoutAvgPtRelSystHighMy[m] ->SetBinContent(typefit+1, 0);
      hhoutAvgPtRelSystLowMy[m] ->SetBinContent(typefit+1, 0);

      //1.B CALCULATE STATISTICAL ERROR OF EXTRAPOLATED YIELD
      /* WRONG: I am doing 4 measurements of the same quantity, but 4 fits with 4 different functions. The statistical error on the extraplated part should not be sigma / sqrt(4) but ~sigma (I can take an average sigma of the 4 functions) */
      YieldErrStatLowPtAvg[m]+= fit_MTscaling[m][typefit]->IntegralError(0,LowRangeSpectrumPart[m],fFitResultPtr0[m][typefit] ->GetParams(),(fFitResultPtr0[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() );
      YieldErrStatHighPtAvg[m]+= fit_MTscalingBis[m][typefit]->IntegralError(UpRangeSpectrumPart[m],20,fFitResultPtr1[m][typefit] ->GetParams(),(fFitResultPtr1[m][typefit]->GetCovarianceMatrix()).GetMatrixArray() );

      //1.C CALCULATE SYST UNCERTAIANTY OF EXTRAPOLATED YIELD
      YieldExtrLowPtAvgUp[m]+= fit_MTscalingUp[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvgUp[m]+= fit_MTscalingUp[m][typefit]->Integral(UpRangeSpectrumPart[m],20);     
      YieldExtrLowPtAvgDown[m]+= fit_MTscalingDown[m][typefit]->Integral(0,LowRangeSpectrumPart[m]);
      YieldExtrHighPtAvgDown[m]+= fit_MTscalingDown[m][typefit]->Integral(UpRangeSpectrumPart[m],20);

      //1.D SYSTEMATIC RELATED TO CHOICE OF FIT FUNCTION
      if (typefit==0){
	YieldExtrMaxLowPt[m] =  YieldExtrLowPt[m][typefit];
	YieldExtrMaxHighPt[m] =  YieldExtrHighPt[m][typefit];
	YieldExtrMinLowPt[m] =  YieldExtrLowPt[m][typefit];
	YieldExtrMinHighPt[m] =  YieldExtrHighPt[m][typefit];
      }
      if ((YieldExtrLowPt[m][typefit]) > YieldExtrMaxLowPt[m]) {
	YieldExtrMaxLowPt[m] = YieldExtrLowPt[m][typefit];
      }
      if ((YieldExtrLowPt[m][typefit]) < YieldExtrMinLowPt[m]) {
	YieldExtrMinLowPt[m] = YieldExtrLowPt[m][typefit];
      }

      if ((YieldExtrHighPt[m][typefit]) > YieldExtrMaxHighPt[m]) {
	YieldExtrMaxHighPt[m] = YieldExtrHighPt[m][typefit];
      }
      if ((YieldExtrHighPt[m][typefit]) < YieldExtrMinHighPt[m]) {
	YieldExtrMinHighPt[m] = YieldExtrHighPt[m][typefit];
      }

      //2.A CALCULATE AVERAGE PT VS MULT FROM FIT
      Int_t numInt = 600; //number of intervals in [0, UpRangePtInterval]
      Float_t Pti=0;
      Float_t UpRangePtInterval = 20;
      Float_t DeltaPt=	UpRangePtInterval/numInt; //interval width

      for (Int_t i=0; i<= numInt; i++){
	Pti = UpRangePtInterval/numInt*i;
	AvgPtFit[m][typefit] += fit_MTscalingBis[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingBis[m][typefit]->Integral(0,UpRangePtInterval);
	//2.C CALCULATE SYST UNCERTAIANTY OF AVG PT FROM FIT
	AvgPtFitHard[m][typefit] += fit_MTscalingHard[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingHard[m][typefit]->Integral(0,UpRangePtInterval);
        AvgPtFitSoft[m][typefit] += fit_MTscalingSoft[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingSoft[m][typefit]->Integral(0,UpRangePtInterval);
      } 

      AvgPt[m] += AvgPtFit[m][typefit];
      AvgPtHard[m] += AvgPtFitHard[m][typefit];
      AvgPtSoft[m] += AvgPtFitSoft[m][typefit];
      hhoutAvgPtMy[m] ->SetBinContent(typefit+1,  AvgPtFit[m][typefit]);

      //2.B CALCULATE STATISTICAL ERROR OF AVG PT FROM FIT
      fit_MTscalingTemp[m][typefit]= (TF1*)       fit_MTscaling[m][typefit]->Clone(      nameMTscaling[m][typefit]+ "_Temp");
      Int_t IterNum =200; //500
      fHistAvgPtDistr[m][typefit] = new TH1F (Form("fHistAvgPtDistr_m%i_typefit%i", m, typefit), Form("fHistAvgPtDistr_m%i_typefit%i", m, typefit), 10000, LimInfPtvsMult, LimSupPtvsMult);

      for (Int_t i= 0; i<IterNum;i++ ){
	Float_t Temp=0;
	for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
	  //extract a random number
	  gRandom->SetSeed(i*fHistSpectrumStat[m]->GetNbinsX()+b);
	  Temp=	  gRandom->Gaus(0, fHistSpectrumStat[m]->GetBinError(b));
	  fHistSpectrumTemp[m]->SetBinContent(b,fHistSpectrumStat[m]->GetBinContent(b) + Temp);
	  fHistSpectrumTemp[m]->SetBinError(b,fHistSpectrumStat[m]->GetBinError(b)/fHistSpectrumStat[m]->GetBinContent(b)*fHistSpectrumTemp[m]->GetBinContent(b));
	  //	  cout <<" stat " <<  fHistSpectrumStat[m]->GetBinContent(b) << " temp " <<  fHistSpectrumTemp[m]->GetBinContent(b) << endl;
	}
	fit_MTscalingTemp[m][typefit]->SetRange(LowRange[m], UpRange[m]);
	fFitResultPtrTemp[m][typefit]=fHistSpectrumTemp[m]->Fit(fit_MTscalingTemp[m][typefit],"SR0QI");
	fit_MTscalingTemp[m][typefit]->SetRange(0,20);

	if (typefit==0){ //--------------------------------
	  //to compute avg pt from spectrum and stat error of avg pt obtained from spectrum!
	  for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
	    if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowPtLimitForAvgPtFS[m]) continue;
	    AvgPtFSNum[m] +=  fHistSpectrumStat[m]->GetBinContent(b) *  fHistSpectrumStat[m]->GetBinCenter(b)* fHistSpectrumStat[m]->GetBinWidth(b);
	    AvgPtFSDenom[m] +=  fHistSpectrumStat[m]->GetBinContent(b) *  fHistSpectrumStat[m]->GetBinWidth(b);
	    AvgPtFSNumTemp[m] +=  fHistSpectrumTemp[m]->GetBinContent(b) *  fHistSpectrumTemp[m]->GetBinCenter(b)* fHistSpectrumTemp[m]->GetBinWidth(b);
	    AvgPtFSDenomTemp[m] +=  fHistSpectrumTemp[m]->GetBinContent(b) *  fHistSpectrumTemp[m]->GetBinWidth(b);
	  }
	  AvgPtFS[m] = AvgPtFSNum[m]/AvgPtFSDenom[m];
	  AvgPtFSTemp[m] = AvgPtFSNumTemp[m]/AvgPtFSDenomTemp[m];
	}//------------------------------------------------

	AvgPtFitTemp[m][typefit]=0;
	for (Int_t l=0; l<= numInt; l++){
	  Pti = UpRangePtInterval/numInt*l;
	  AvgPtFitTemp[m][typefit] += fit_MTscalingTemp[m][typefit]->Eval(Pti)*Pti*DeltaPt/fit_MTscalingTemp[m][typefit]->Integral(0,UpRangePtInterval);
	  //	cout << " AVERAGE PT at iteration " << i << " of "<< numInt << " (pt = "  << Pti << "): " << 	AvgPtFit[m][typefit] << endl;
	}

	fHistAvgPtDistr[m][typefit] ->Fill(AvgPtFitTemp[m][typefit]);
	if (typefit ==0) fHistAvgPtFSDistr[m] ->Fill(AvgPtFSTemp[m]);
      } //end number of iterations

      TF1 *gaus = (TF1 *)gROOT->GetFunction("gaus");
      fHistAvgPtDistr[m][typefit]->Rebin(nrebin[m]);
      gaus->SetParameter(2, fHistAvgPtDistr[m][typefit]->GetRMS());
      fHistAvgPtDistr[m][typefit]->Fit(gaus, "q");
      AvgPtStatErrFit[m][typefit] = 	fHistAvgPtDistr[m][typefit]->GetRMS()/	fHistAvgPtDistr[m][typefit]->GetMean()* AvgPtFit[m][typefit];
      //AvgPtStatErrFit[m][typefit] = 	gaus->GetParameter(2)/gaus->GetParameter(1)* AvgPtFit[m][typefit];
      AvgPtStatErr[m]+= AvgPtStatErrFit[m][typefit];
    } //end loop typefit
    
    //2.C SYSTEMATIC RELATED TO CHOICE OF FIT FUNCTION
    AvgPtFitMaxFit[m] = 0;
    AvgPtFitMinFit[m] = 1000;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (nameFit[typefit]!="Levi") continue;
      if(AvgPtFitMaxFit[m] < AvgPtFit[m][typefit]){
	AvgPtFitMaxFit[m] = AvgPtFit[m][typefit];
      }
      if(AvgPtFitMinFit[m] > AvgPtFit[m][typefit]){
	AvgPtFitMinFit[m] = AvgPtFit[m][typefit];
      }
    }

    YieldExtrHighPtAvg[m]=YieldExtrHighPtAvg[m]/numfittipoEff;
    YieldExtrLowPtAvg[m]=YieldExtrLowPtAvg[m]/numfittipoEff;
    YieldExtrHighPtAvgUp[m]=YieldExtrHighPtAvgUp[m]/numfittipoEff;
    YieldExtrHighPtAvgDown[m]=YieldExtrHighPtAvgDown[m]/numfittipoEff;
    YieldExtrLowPtAvgUp[m]=YieldExtrLowPtAvgUp[m]/numfittipoEff;
    YieldExtrLowPtAvgDown[m]=YieldExtrLowPtAvgDown[m]/numfittipoEff;
    YieldErrStatHighPtAvg[m]=YieldErrStatHighPtAvg[m]/numfittipoEff;
    YieldErrStatLowPtAvg[m]=YieldErrStatLowPtAvg[m]/numfittipoEff;
    /*WRONG
    YieldErrStatHighPtAvg[m]=sqrt(YieldErrStatHighPtAvg[m])/numfittipoEff;
    YieldErrStatLowPtAvg[m]=sqrt(YieldErrStatLowPtAvg[m])/numfittipoEff;
    */

    //1.A EXTRAPOLATED YIELD
    YieldExtr[m] =     YieldExtrHighPtAvg[m]+    YieldExtrLowPtAvg[m]; 

    //1.B STATISTICAL ERROR OF EXTRAPOLATED YIELD
    YieldExtrErrStat[m] = sqrt(    pow(YieldErrStatHighPtAvg[m],2)+    pow(YieldErrStatLowPtAvg[m],2)); 

    //1.C SYSTEMATIC ERROR OF EXTRAPOLATED YIELD
    YieldErrSystExtrLowPt[m] = (YieldExtrLowPtAvgUp[m]- YieldExtrLowPtAvgDown[m])/2;
    YieldErrSystExtrHighPt[m] = (YieldExtrHighPtAvgUp[m]- YieldExtrHighPtAvgDown[m])/2;
    if (!isErrorAssumedPtCorr){
      YieldExtrErrSist[m]=sqrt(pow(  YieldErrSystExtrLowPt[m],2) + pow(  YieldErrSystExtrHighPt[m] ,2)) ;
    }
    else {
      YieldExtrErrSist[m]= YieldErrSystExtrLowPt[m]  + YieldErrSystExtrHighPt[m] ;
    }

    //1.D SYSTEMATIC ERROR OF YIELD RELATED TO CHOICE OF FIT FUNCTION
    Yield4FitErrSistHighPt[m] = (YieldExtrMaxHighPt[m]-YieldExtrMinHighPt[m])/2; 
    Yield4FitErrSistLowPt[m] = (YieldExtrMaxLowPt[m]-YieldExtrMinLowPt[m])/2; 
    YieldExtrErrSistFourFit[m]=sqrt(pow(    Yield4FitErrSistHighPt[m],2) + pow(    Yield4FitErrSistLowPt[m],2));

    //STATISTICAL AND SYSTEMATIC ERROR OF MEASURED YIELD
    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      if (NPtV0[v] >= UpRangeSpectrumPart[m]) continue; 
      YieldSpectrum[m] +=  ( fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1));
      YieldSpectrumErrStat[m] += pow ( fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1),2);
      if (!isErrorAssumedPtCorr){
      YieldSpectrumErrSist[m] += pow ( fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1),2);
      }
      else {
	YieldSpectrumErrSist[m] += fHistSpectrumSistAll[m]->GetBinError(v+1)*fHistSpectrumSistAll[m]->GetBinWidth(v+1);
      }
    }

    YieldSpectrumErrStat[m]=sqrt( YieldSpectrumErrStat[m]);
    if (!isErrorAssumedPtCorr){
      YieldSpectrumErrSist[m]=sqrt( YieldSpectrumErrSist[m]);
    }

    //FINAL RESULTS FOR YIELD
    Yield[m] =  YieldSpectrum[m]+YieldExtr[m];
    YieldErrStat[m] = sqrt(pow(YieldSpectrumErrStat[m],2) + pow(YieldExtrErrStat[m],2));

    if (!isErrorAssumedPtCorr){
      YieldErrSist[m] = sqrt(pow(YieldSpectrumErrSist[m],2) +  pow(YieldErrSystExtrLowPt[m],2) + pow( YieldErrSystExtrHighPt[m],2));
    }
    else {
      YieldErrSist[m] = YieldSpectrumErrSist[m] +  YieldErrSystExtrLowPt[m] + YieldErrSystExtrHighPt[m];
    }
    YieldErrSistMy[m] = YieldErrSist[m];
    YieldErrSist[m]= sqrt(pow(YieldErrSist[m],2) + pow(Yield4FitErrSistLowPt[m],2) + pow(Yield4FitErrSistHighPt[m],2) );

    //******************************************************************
    //FINAL VALUES TO AVERAGE PT VS MULT OBTAINED FROM FIT
    AvgPt[m] = AvgPt[m]/numfittipoEff; //AVERAGE PT
    AvgPtStatErr[m] = AvgPtStatErr[m]/numfittipoEff; //STATISTICAL ERROR OF AVG PT
    AvgPtSistErrFit[m] =  (AvgPtFitMaxFit[m] - AvgPtFitMinFit[m])/2; //SYST ERROR OF AVG PT RELATED TO CHOICE OF FIT FUNCTION
    if (isMeanMacro)     AvgPtSistErrFitMeanMacro[m] =  (AvgPtMeanMacroMaxFit[m] - AvgPtMeanMacroMinFit[m])/2; //SYST ERROR OF AVG PT CALCULATED BY YIELDMEAN MACRO RELATED TO CHOICE OF FIT FUNCTION

    //SYST ERROR OF AVG PT:
    AvgPtSoft[m] = AvgPtSoft[m]/numfittipoEff; 
    AvgPtHard[m] = AvgPtHard[m]/numfittipoEff;
    AvgPtSistErr[m] = TMath::Abs((AvgPtHard[m]- AvgPtSoft[m]))/2;
    //******************************************************************

    //******************************************************************
    //CALCULATE AVERAGE PT VS MULT FROM SPECTRUM
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowPtLimitForAvgPtFS[m]) continue;
      AvgPtFSNumHard[m] +=  fHistSpectrumStatHard[m]->GetBinContent(b) *  fHistSpectrumStatHard[m]->GetBinCenter(b)* fHistSpectrumStatHard[m]->GetBinWidth(b);
      AvgPtFSDenomHard[m] +=  fHistSpectrumStatHard[m]->GetBinContent(b) *  fHistSpectrumStatHard[m]->GetBinWidth(b);
      AvgPtFSNumSoft[m] +=  fHistSpectrumStatSoft[m]->GetBinContent(b) *  fHistSpectrumStatSoft[m]->GetBinCenter(b)* fHistSpectrumStatSoft[m]->GetBinWidth(b);
      AvgPtFSDenomSoft[m] +=  fHistSpectrumStatSoft[m]->GetBinContent(b) *  fHistSpectrumStatSoft[m]->GetBinWidth(b);
    }

    Float_t DeltaP =0;
    Float_t AvP =0;
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowPtLimitForAvgPtFS[m]) continue;
      DeltaP = fHistSpectrumStat[m]->GetBinWidth(b);
      AvP = fHistSpectrumStat[m]->GetBinCenter(b);
      //...not sure how I computed the errors below....:)
      /*
      AvgPtStatErrFS[m] += pow((DeltaP*AvP*AvgPtFSDenom[m]-pow(DeltaP,2)*AvP*fHistSpectrumStat[m]->GetBinContent(b))/pow(AvgPtFSDenom[m],2) *fHistSpectrumStat[m]->GetBinError(b) ,2);
      AvgPtSistErrFS[m] += pow((DeltaP*AvP*AvgPtFSDenom[m]-pow(DeltaP,2)*AvP*fHistSpectrumSist[m]->GetBinContent(b))/pow(AvgPtFSDenom[m],2) *fHistSpectrumSist[m]->GetBinError(b) ,2);
      */
    }

    AvgPtFSHard[m] = AvgPtFSNumHard[m]/AvgPtFSDenomHard[m];
    AvgPtFSSoft[m] = AvgPtFSNumSoft[m]/AvgPtFSDenomSoft[m];
    AvgPtSistErrFS[m] = TMath::Abs(AvgPtFSHard[m]-AvgPtFSSoft[m])/2;
    AvgPtStatErrFS[m] = fHistAvgPtFSDistr[m]->GetRMS()/	fHistAvgPtFSDistr[m]->GetMean()* AvgPtFS[m];
    //***************************************************

    //CALCULATE DISCRETIZATION ERROR
    for (Int_t b=1; b<=  fHistSpectrumStat[m]->GetNbinsX(); b++){
      if (fHistSpectrumStat[m]->GetXaxis()->GetBinLowEdge(b) < LowRange[m]) continue;
      //      AvgPtDiscrErr[m] += pow(fHistSpectrumStat[m]->GetBinContent(b),2) * pow(fHistSpectrumStat[m]->GetBinWidth(b),4)/12; 
      AvgPtDiscrErr[m] += fHistSpectrumStat[m]->GetBinContent(b) * pow(fHistSpectrumStat[m]->GetBinWidth(b),2); 
    }
    //    AvgPtDiscrErr[m] = sqrt(AvgPtDiscrErr[m])/YieldSpectrum[m];
    AvgPtDiscrErr[m] = AvgPtDiscrErr[m]/YieldSpectrum[m]/sqrt(12);

    //FILLING HISTOGRAMS
    hhoutYield[m] ->SetBinContent(numfittipo+1, hhoutYieldAvg[m]/numfittipoEff);
    hhoutYieldMy[m] ->SetBinContent(numfittipo+1, YieldExtr[m]);
    hhoutAvgPt[m] ->SetBinContent(numfittipo+1, hhoutAvgPtAvg[m]/numfittipoEff);
    hhoutAvgPtMy[m] ->SetBinContent(numfittipo+1, AvgPt[m]);

    hhoutYieldRatioToMine[m] = (TH1F*)    hhoutYield[m]->Clone(Form("hhoutYieldRatioToMine_m%i", m));
    hhoutYieldRatioToMine[m] ->Divide(hhoutYieldMy[m]);
    hhoutYieldRatioToMine[m] ->GetYaxis()->SetRangeUser(0.9,1.1);
    hhoutYieldRatioToMine[m]->SetLineColor(ColorMult[m]);

    hhoutAvgPtRatioToMine[m] = (TH1F*)    hhoutAvgPt[m]->Clone(Form("hhoutAvgPtRatioToMine_m%i", m));
    hhoutAvgPtRatioToMine[m] ->Divide(hhoutAvgPtMy[m]);
    hhoutAvgPtRatioToMine[m] ->GetYaxis()->SetRangeUser(0.9,1.1);
    hhoutAvgPtRatioToMine[m]->SetLineColor(ColorMult[m]);

    hhoutRelStat[m] ->SetBinContent(numfittipo+1, hhoutRelStatAvg[m]/numfittipoEff);
    hhoutRelSystHigh[m] ->SetBinContent(numfittipo+1, hhoutRelSystHighAvg[m]/numfittipoEff);
    hhoutRelSystLow[m] ->SetBinContent(numfittipo+1, hhoutRelSystLowAvg[m]/numfittipoEff);
    hhoutAvgPtRelStat[m] ->SetBinContent(numfittipo+1, hhoutAvgPtRelStatAvg[m]/numfittipoEff);
    hhoutAvgPtRelSystHigh[m] ->SetBinContent(numfittipo+1, hhoutAvgPtRelSystHighAvg[m]/numfittipoEff);
    hhoutAvgPtRelSystLow[m] ->SetBinContent(numfittipo+1, hhoutAvgPtRelSystLowAvg[m]/numfittipoEff);

    hhoutRelStat[m] ->Scale(1./Yield[m]);
    hhoutRelSystHigh[m] ->Scale(1./Yield[m]);
    hhoutRelSystLow[m] ->Scale(1./Yield[m]);
    hhoutAvgPtRelStat[m] ->Scale(1./AvgPt[m]);
    hhoutAvgPtRelSystHigh[m] ->Scale(1./AvgPt[m]);
    hhoutAvgPtRelSystLow[m] ->Scale(1./AvgPt[m]);

    hhoutRelStat[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutRelSystHigh[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutRelSystLow[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelStat[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutAvgPtRelSystHigh[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelSystLow[m]->GetYaxis()->SetRangeUser(0,0.1);

    hhoutRelStat[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystHigh[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystLow[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelStat[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystHigh[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystLow[m]->SetLineColor(ColorMult[m]);

    hhoutRelStatMy[m] ->SetBinContent(numfittipo+1, YieldErrStat[m]);
    hhoutRelStatMy[m] ->Scale(1./Yield[m]);
    hhoutRelSystHighMy[m] ->SetBinContent(numfittipo+1, YieldErrSistMy[m]/Yield[m]);
    hhoutRelSystLowMy[m] ->SetBinContent(numfittipo+1, YieldErrSistMy[m]/ Yield[m]);

    hhoutRelStatMy[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutRelSystHighMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutRelSystLowMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutRelStatMy[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystHighMy[m]->SetLineColor(ColorMult[m]);
    hhoutRelSystLowMy[m]->SetLineColor(ColorMult[m]);

    hhoutAvgPtRelStatMy[m] ->SetBinContent(numfittipo+1, AvgPtStatErr[m]/AvgPt[m]);
    hhoutAvgPtRelSystHighMy[m] ->SetBinContent(numfittipo+1, AvgPtSistErr[m]/AvgPt[m]);
    hhoutAvgPtRelSystLowMy[m] ->SetBinContent(numfittipo+1, AvgPtSistErr[m]/ AvgPt[m]);

    hhoutAvgPtRelStatMy[m]->GetYaxis()->SetRangeUser(0,0.05);
    hhoutAvgPtRelSystHighMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelSystLowMy[m]->GetYaxis()->SetRangeUser(0,0.1);
    hhoutAvgPtRelStatMy[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystHighMy[m]->SetLineColor(ColorMult[m]);
    hhoutAvgPtRelSystLowMy[m]->SetLineColor(ColorMult[m]);

    cout << "\n\n\e[35m*** Yields *** " << SmoltLegend[m] << "\e[39m " << endl;
    cout << "\nSpectra yield: " << endl;
    for(Int_t v = PtBinMin[m]; v < numPtV0Max; v++){
      cout << NPtV0[v] << " < pt < " << NPtV0[v+1] << ": " << fHistSpectrumStat[m]->GetBinContent(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) <<  " +- " << fHistSpectrumStat[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) << " (stat.) +- " <<  fHistSpectrumSist[m]->GetBinError(v+1)*fHistSpectrumStat[m]->GetBinWidth(v+1) << " (sist.)" << endl;
    }

    cout << "\npT integrated yield: " << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) +- "<< YieldErrSist[m] << " (sist) "<< endl;
    cout << " (stat rel: " <<  YieldErrStat[m]/ Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldErrSist[m]/ Yield[m] << ") " << endl;

    cout << "\nYield from spectrum: " <<   YieldSpectrum[m] << " +- " << YieldSpectrumErrStat[m]<<" (stat) +- " <<YieldSpectrumErrSist[m]<< " (sist)"  << endl;
    cout << " (stat rel: " <<  YieldSpectrumErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldSpectrumErrSist[m]/Yield[m] << ") "<<endl;
    cout << "Errors relative to yield from spectrum:" << endl;
    cout << " (stat rel: " <<  YieldSpectrumErrStat[m]/YieldSpectrum[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldSpectrumErrSist[m]/YieldSpectrum[m] << ") "<<endl;

    cout << "\nYield from extrapolation: " <<     YieldExtr[m] << " +- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSistFourFit[m]<< " (sist 4 fit) " << YieldExtrErrSist[m] << " (sist extr) "<< endl;
    cout << " (stat rel: " <<  YieldExtrErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist 4 fit rel: " <<  YieldExtrErrSistFourFit[m]/Yield[m] << ") "<<endl;
    cout << " (sist rel: " <<  YieldExtrErrSist[m]/Yield[m] << ") "<<endl;
    cout << "Errors relative to extrapolated yield:" << endl;
    cout << " (stat rel: " <<  YieldExtrErrStat[m]/YieldExtr[m] << ") " << endl;
    cout << " (sist 4 fit rel: " <<  YieldExtrErrSistFourFit[m]/YieldExtr[m] << ") "<<endl;
    cout << " (sist rel: " <<  YieldExtrErrSist[m]/YieldExtr[m] << ") "<<endl;

    cout << "\nFraction of yield from low pt extrapolation " <<  YieldExtrLowPtAvg[m]/Yield[m] << endl;
    cout << "Fraction of yield from high pt extrapolation " <<  YieldExtrHighPtAvg[m]/Yield[m] << endl;
    cout << "\nFit range: " << LowRange[m] << " - " << UpRange[m] << endl;
    cout << "Bin content range: " << LowRangeSpectrumPart[m] << " - " << UpRangeSpectrumPart[m] << endl;

    cout << "\e[35m**** Avg pt obtained with different fit functions **** " << SmoltLegend[m] << "\e[39m" << endl;
    for (Int_t typefit =0; typefit<numfittipo; typefit++){
      if (nameFit[typefit]!="Levi") continue;
      cout << "\n ***" << endl;
      cout << "Avg pt (from fit) " << nameFit[typefit]<< ": " << AvgPtFit[m][typefit]<< " GeV/c" <<endl;
      cout << "Avg pt (from distribution of <pt> obtained with random variation of spectrum within statistical uncertainty): " << endl;
      //      cout << "Avg pt (from last random variation of spectrum) " << AvgPtFitTemp[m][typefit]<< " GeV/c " << endl;
      cout << " MEAN: " << fHistAvgPtDistr[m][typefit]->GetMean()<< ", RMS: " << fHistAvgPtDistr[m][typefit]->GetRMS() << endl;
      cout << " bin width of <pt> distribution (check it is much smaller than RMS): " << fHistAvgPtDistr[m][typefit]->GetXaxis()->GetBinWidth(1)<< endl;
    }
    cout << "\n\n\e[35m**** Avg pt from avg of different fit functions **** " << SmoltLegend[m] <<"\e[39m" << endl;
    cout << "Average pt: " << AvgPt[m]<< " +- " << AvgPtStatErr[m] << " (stat.) +- " << AvgPtSistErr[m] << " (syst. no function choice) +- "<<AvgPtSistErrFit[m] <<" (syst. function choice)"<<  endl;
    cout << "Rel errors:  " <<  AvgPtStatErr[m]/AvgPt[m] << " (stat.) +- " << AvgPtSistErr[m]/AvgPt[m] << " (syst. no function choice) +- "<< AvgPtSistErrFit[m]/AvgPt[m]<< " (syst. function choice)"<<  endl;
    
    cout << "\nAvg pt from spectrum (no fit, calculated for pt > " << 	LowPtLimitForAvgPtFS[m]  << "): " << AvgPtFS[m]<< " +- " << AvgPtStatErrFS[m] << " (stat.) +- " << AvgPtSistErrFS[m] << " (syst. no function choice) "<< endl; 
    cout << " MEAN: " << fHistAvgPtFSDistr[m]->GetMean()<< ", RMS: " << fHistAvgPtFSDistr[m]->GetRMS() << endl;  
    cout << "Rel errors:  " <<  AvgPtStatErrFS[m]/AvgPtFS[m] << " (stat.) +- " << AvgPtSistErrFS[m]/AvgPtFS[m] << " (syst. no function choice) +- "<< endl;
    cout << "Mean of fHistSpectrumStat: " << fHistSpectrumStat[m]->GetMean() << endl;

    if (isMeanMacro){
      AvgPtMeanMacro[m] =  hhoutAvgPt[m]->GetBinContent(numfittipo+1);
      AvgPtStatErrMeanMacro[m] = hhoutAvgPtRelStat[m]->GetBinContent(numfittipo+1);
      AvgPtSistErrMeanMacro[m] = (hhoutAvgPtRelSystHigh[m]->GetBinContent(numfittipo+1) + hhoutAvgPtRelSystLow[m]->GetBinContent(numfittipo+1))/2;
    cout << "\nAverage pt from Mean Macro: " << AvgPtMeanMacro[m] << " +- " <<  AvgPtStatErrMeanMacro[m] << " (stat.) +- " << AvgPtSistErrMeanMacro[m]<< " (syst. no function choice) +- "<< "-" <<" (syst. function choice)"<<  endl;
    cout << "Rel errors:  " <<  AvgPtStatErrMeanMacro[m]/AvgPtMeanMacro[m] << " (stat.) +- " << AvgPtSistErrMeanMacro[m]/AvgPtMeanMacro[m] << " (syst. no function choice) " <<  AvgPtSistErrFitMeanMacro[m]/AvgPtMeanMacro[m]<< " (syst. function choice)"<<  endl;
    }

    cout << "Discretization error: " << AvgPtDiscrErr[m] << " (Rel to AvgPt from mean macro: " << AvgPtDiscrErr[m]/AvgPtMeanMacro[m] << ") " << endl;
  }//end loop m

  cout << "\n\nYields : " << endl;
  for (Int_t m=0; m<nummolt +1; m++){
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    if (isppHM && m<2) continue;
    cout <<"\n\e[32m************ Multiplicity: " << SmoltLegend[m] << " *****************\e[39m"  << endl;

    cout << "pT integrated yield: " << Yield[m] << "+-" << YieldErrStat[m]<< " (stat) +- "<< YieldErrSist[m] << " (sist) "<< endl;
    cout << " (stat rel: " <<  YieldErrStat[m]/ Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldErrSist[m]/ Yield[m] << ") " << endl;

    cout << "\nYield from spectrum: " <<   YieldSpectrum[m] << " +- " << YieldSpectrumErrStat[m]<<" (stat) +- " <<YieldSpectrumErrSist[m]<< " (sist)"  << endl;
    cout << " (stat rel: " <<  YieldSpectrumErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist rel: " <<  YieldSpectrumErrSist[m]/Yield[m] << ") "<<endl;

    cout << "\nYield from extrapolation: " <<     YieldExtr[m] << " +- " <<     YieldExtrErrStat[m]<< " (stat)  +- " << YieldExtrErrSistFourFit[m]<< " (sist 4 fit) " << YieldExtrErrSist[m] << " (sist extr) "<< endl;
    cout << " (stat rel: " <<  YieldExtrErrStat[m]/Yield[m] << ") " << endl;
    cout << " (sist 4 fit rel: " <<  YieldExtrErrSistFourFit[m]/Yield[m] << ") "<<endl;
    cout << " (sist rel: " <<  YieldExtrErrSist[m]/Yield[m] << ") "<<endl;

    cout << "\nFraction of yield from low pt extrapolation " <<  YieldExtrLowPtAvg[m]/Yield[m] << endl;
    cout << "Fraction of yield from high pt extrapolation " <<  YieldExtrHighPtAvg[m]/Yield[m] << endl;
    cout << "\nFit range: " << LowRange[m] << " - " << UpRange[m] << endl;
    cout << "Bin content range: " << LowRangeSpectrumPart[m] << " - " << UpRangeSpectrumPart[m] << endl;
  }

  TLegend *legendYield=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendPtvsMult=new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend *legendYieldErr=new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legendPtErr=new TLegend(0.6, 0.6, 0.9, 0.9);

  TF1 * pol0 = new TF1 ("pol0", "pol0", 0, 30);
  TString titleYieldX="dN_{ch}/d#eta";
  TString titleYieldYType[2]={"N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi", "N_{#Xi}/N_{Trigg} 1/#Delta#eta #Delta#phi"};
  TString titlePtvsMultYType[2]={"<p^{K^{0}_{S}}_{T}> (GeV/c)", "<p^{#Xi}_{T} (GeV/c)>"};
  TString titleYieldY;
  if(type==0) titleYieldY=titleYieldYType[0];
  else if(type==8) titleYieldY=titleYieldYType[1];
  TString titlePtvsMultY;
  if(type==0) titlePtvsMultY=titlePtvsMultYType[0];
  else if(type==8) titlePtvsMultY=titlePtvsMultYType[1];
  TString titleYield[3]={"In-jet", "Out-of-jet", "Inclusive"};
  TH1F* fHistYieldStat=new TH1F ("fHistYieldStatD","fHistYieldStatD",200,0,40);
  TH1F* fHistYieldSist=new TH1F ("fHistYieldSist","fHistYieldSist",200,0,40);

  TH1F* fHistPtvsMultStat = new TH1F ("fHistPtvsMultStat","fHistPtvsMultStat",200,0,40);
  TH1F* fHistPtvsMultSist = new TH1F ("fHistPtvsMultSist","fHistPtvsMultSist",200,0,40);
  TH1F* fHistPtvsMultStatFromSpectrum = new TH1F ("fHistPtvsMultStatFromSpectrum","fHistPtvsMultStatFromSpectrum",200,0,40);
  TH1F* fHistPtvsMultSistFromSpectrum = new TH1F ("fHistPtvsMultSistFromSpectrum","fHistPtvsMultSistFromSpectrum",200,0,40);
  TH1F* fHistPtvsMultStatMeanMacro = new TH1F ("fHistPtvsMultStatMeanMacro","fHistPtvsMultStatMeanMacro",200,0,40);
  TH1F* fHistPtvsMultSistMeanMacro = new TH1F ("fHistPtvsMultSistMeanMacro","fHistPtvsMultSistMeanMacro",200,0,40);

  TH1F* fHistPtStatRelErr = new TH1F ("fHistPtStatRelErr","fHistPtStatRelErr",200,0,40);
  TH1F* fHistPtSistRelErr = new TH1F ("fHistPtSistRelErr","fHistPtSistRelErr",200,0,40);
  TH1F* fHistPtStatRelErrFromSpectrum = new TH1F ("fHistPtStatRelErrFromSpectrum","fHistPtStatRelErrFromSpectrum",200,0,40);
  TH1F* fHistPtSistRelErrFromSpectrum = new TH1F ("fHistPtSistRelErrFromSpectrum","fHistPtSistRelErrFromSpectrum",200,0,40);
  TH1F* fHistPtStatRelErrMeanMacro = new TH1F ("fHistPtStatRelErrMeanMacro","fHistPtStatRelErrMeanMacro",200,0,40);
  TH1F* fHistPtSistRelErrMeanMacro = new TH1F ("fHistPtSistRelErrMeanMacro","fHistPtSistRelErrMeanMacro",200,0,40);
  TH1F* fHistPtSistRelErrFit = new TH1F ("fHistPtSistRelErrFit","fHistPtSistRelErrFit",200,0,40);

  TH1F* fHistYieldSistNoExtr=new TH1F ("fHistYieldSistNoExtr","fHistYieldSistNoExtr",200,0,40);
  TH1F* fHistYieldStatRelErr=new TH1F ("fHistYieldStatRelErr","fHistYieldStatRelErr",200,0,40);
  TH1F* fHistYieldSistRelErr=new TH1F ("fHistYieldSistRelErr","fHistYieldSistRelErr",200,0,40);
  TH1F* fHistYieldSistNoExtrRelErr=new TH1F ("fHistYieldSistNoExtrRelErr","fHistYieldSistNoExtrRelErr",200,0,40);
  TH1F* fHistYieldStatEtaEff=new TH1F ("fHistYieldStatEtaEff","fHistYieldStatEtaEff",200,0,40);
  TH1F* fHistYieldStatEtaEffRatio;
  TH1F* fHistYieldStatEtaEffRatioRef;
  TH1F* fHistYieldStatNormCorr=new TH1F ("fHistYieldStatNormCorr","fHistYieldStatNormCorr",200,0,40);
  TH1F* fHistYieldStatNormCorrRatio;
  TH1F* fHistYieldStatNormCorrRatioRef;

  Float_t   dNdEta[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  if (isppHM) {
    dNdEta[0] = 39.40;
    dNdEta[1] = 36.89;
    dNdEta[2] = 35.16;
    dNdEta[3] = 32.57;
    dNdEta[4] = 30.43;
    dNdEta[5] = 31.5;
    if (MultBinning==1){
      dNdEta[0] =0;
      dNdEta[1] =0;
      dNdEta[2] =36.29; //values from 16l with 18d8 MC
      dNdEta[3] =32.57;
      dNdEta[4] =30.43;
      dNdEta[5] = 31.5;
    }
  }

  Int_t LimSupChi=120;
  if (type==1) LimSupChi=100;
  canvasFitResult->cd();
  for (Int_t typefit=0; typefit<numfittipo; typefit++){
      if (nameFit[typefit]!="Levi") continue;
    StyleHisto( hFitResult[typefit], 0,LimSupChi, ColorFit[typefit], 33, "Multiplicity class", "#chi^{2}/NDF", "#chi^{2} vs multiplicity", 0, 0, 0);
    hFitResult[typefit]->GetYaxis()->SetTitleOffset(1.2);
    hFitResult[typefit]->Draw("same");
    if (typefit==0)    legendfit->Draw("");
  }

  for(Int_t m=0; m<nummolt; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fHistPtvsMultStat->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPt[m]);
    fHistPtvsMultSist->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPt[m]);
    fHistPtvsMultStat->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErr[m]);
    fHistPtvsMultSist->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),sqrt(pow(AvgPtSistErr[m],2) + pow(AvgPtSistErrFit[m],2)));
    fHistPtvsMultStatFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtFS[m]);
    fHistPtvsMultSistFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtFS[m]);
    fHistPtvsMultStatFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErrFS[m]);
    fHistPtvsMultSistFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtSistErrFS[m]);
    if (isMeanMacro){
    fHistPtvsMultStatMeanMacro->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtMeanMacro[m]);
    fHistPtvsMultSistMeanMacro->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtMeanMacro[m]);
    fHistPtvsMultStatMeanMacro->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]), AvgPtStatErrMeanMacro[m]);
    fHistPtvsMultSistMeanMacro->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]), sqrt(pow(AvgPtSistErrMeanMacro[m],2) + pow(AvgPtSistErrFitMeanMacro[m]/AvgPt[m]*AvgPtMeanMacro[m],2)));
    }

    fHistPtStatRelErr->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErr[m]/AvgPt[m]);
    fHistPtStatRelErr->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),0);
    fHistPtSistRelErr->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErr[m]/AvgPt[m]);
    fHistPtSistRelErr->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);

    fHistPtStatRelErrFromSpectrum->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErrFS[m]/AvgPtFS[m]);
    fHistPtStatRelErrFromSpectrum->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),0);
    fHistPtSistRelErrFromSpectrum->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErrFS[m]/AvgPtFS[m]);
    fHistPtSistRelErrFromSpectrum->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);

    if (isMeanMacro){
      fHistPtStatRelErrMeanMacro->SetBinContent(fHistPtvsMultStat->FindBin(dNdEta[m]),AvgPtStatErrMeanMacro[m]/AvgPtMeanMacro[m]);
      fHistPtStatRelErrMeanMacro->SetBinError(fHistPtvsMultStat->FindBin(dNdEta[m]),0);
      fHistPtSistRelErrMeanMacro->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErrMeanMacro[m]/AvgPtMeanMacro[m]);
      fHistPtSistRelErrMeanMacro->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);
    }

    fHistPtSistRelErrFit->SetBinContent(fHistPtvsMultSist->FindBin(dNdEta[m]),AvgPtSistErrFitMeanMacro[m]/AvgPtMeanMacro[m]);
    fHistPtSistRelErrFit->SetBinError(fHistPtvsMultSist->FindBin(dNdEta[m]),0);
  }

  for(Int_t m=0; m<nummolt; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fHistYieldStat->SetBinContent(fHistYieldStat->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldStat->SetBinError(fHistYieldStat->FindBin(dNdEta[m]),YieldErrStat[m]);
    fHistYieldSist->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldSist->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),YieldErrSist[m]);
    fHistYieldSistNoExtr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),Yield[m]);
    fHistYieldSistNoExtr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]), YieldSpectrumErrSist[m]);

    fHistYieldStatRelErr->SetBinContent(fHistYieldStat->FindBin(dNdEta[m]),YieldErrStat[m]/Yield[m]);
    fHistYieldStatRelErr->SetBinError(fHistYieldStat->FindBin(dNdEta[m]),0);
    fHistYieldSistRelErr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]),YieldErrSist[m]/Yield[m]);
    fHistYieldSistRelErr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);
    fHistYieldSistNoExtrRelErr->SetBinContent(fHistYieldSist->FindBin(dNdEta[m]), YieldSpectrumErrSist[m]/Yield[m]);
    fHistYieldSistNoExtrRelErr->SetBinError(fHistYieldSist->FindBin(dNdEta[m]),0);

    fHistYieldStat->SetMarkerSize(2);
    fHistYieldSist->SetMarkerSize(2);
    fHistYieldSistNoExtr->SetMarkerSize(2);
    fHistYieldStatRelErr->SetMarkerSize(2);
    fHistYieldSistRelErr->SetMarkerSize(2);
    fHistYieldSistNoExtrRelErr->SetMarkerSize(2);

  }


  Float_t Xl[nummolt+1]= {0};
  Float_t Xh[nummolt+1]= {0};
  Float_t multctrbin[nummolt+1] ={0} ;

  Float_t LimSupdNdEtaAxis = 25;
  if (isppHM) LimSupdNdEtaAxis = 40;

  for (Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    multctrbin[m] =   fHistYieldStat->GetXaxis()->GetBinCenter(  fHistYieldStat->FindBin(dNdEta[m]));
  }

  canvasYield->cd();
  StyleHisto(fHistYieldStat, LimInfYield, LimSupYield, Color[type], 1, titleYieldX, titleYieldY, "Yield vs multiplicity",1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSist, LimInfYield, LimSupYield, Color[type], 1, titleYieldX, titleYieldY, "Yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSistNoExtr, LimInfYield, LimSupYield, Color[type], 1, titleYieldX, titleYieldY, "Yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
     
  fHistYieldStat->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSist->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldSistNoExtr->GetYaxis()->SetTitleOffset(1.2);

  fHistYieldStat->DrawClone("e");
  //  fHistYieldSistNoExtr->SetFillStyle(9);
  fHistYieldSistNoExtr->SetFillStyle(3001);
  fHistYieldSistNoExtr->SetFillColorAlpha(Color[type], 1);
  //  fHistYieldSistNoExtr->Draw("same e2");
  fHistYieldSist->SetFillStyle(0);
  //  fHistYieldSist->Draw("same e2");
  fHistYieldDatiPubblicatiStat->SetMarkerColor(881);
  fHistYieldDatiPubblicatiStat->SetLineColor(881);
  fHistYieldDatiPubblicatiSist->SetMarkerColor(881);
  fHistYieldDatiPubblicatiSist->SetLineColor(881);
  fHistYieldDatiPubblicatiStat->Draw("same");
  fHistYieldDatiPubblicatiSist->SetFillStyle(0);
  fHistYieldDatiPubblicatiSist->Draw("same e2");
  legendYield->AddEntry(fHistYieldSist, "syst.", "ef");
  legendYield->AddEntry(fHistYieldSistNoExtr, "syst. (no extr)", "ef");
  legendYield->AddEntry(fHistYieldStat, "stat.", "pel");
  legendYield->Draw("");
  canvasYield->SaveAs("provay.pdf");

  fHistYieldStat->SetName("fHistYieldStat");
  fHistYieldStat->SetTitle("fHistYieldStat");

  canvasPtvsMult->cd();
  StyleHisto(fHistPtvsMultStat, LimInfPtvsMult, LimSupPtvsMult, Color[type], 33, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity ",1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultSist, LimInfPtvsMult, LimSupPtvsMult, Color[type], 33, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" , 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultStatFromSpectrum, LimInfPtvsMult, LimSupPtvsMult, Color[type], 4, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity ",1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultSistFromSpectrum, LimInfPtvsMult, LimSupPtvsMult, Color[type], 4, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" , 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultStatMeanMacro, LimInfPtvsMult, LimSupPtvsMult, Color[type]+2, 22, titleYieldX, titlePtvsMultY, "<p_{T}> vs multiplicity ",1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtvsMultSistMeanMacro, LimInfPtvsMult, LimSupPtvsMult, Color[type]+2, 22, titleYieldX, titlePtvsMultY, "<p_{T} vs multiplicity" , 1, 0, LimSupdNdEtaAxis);
     
  fHistPtvsMultStat->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSist->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultStatFromSpectrum->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSistFromSpectrum->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultStatMeanMacro->GetYaxis()->SetTitleOffset(1.2);
  fHistPtvsMultSistMeanMacro->GetYaxis()->SetTitleOffset(1.2);

  fHistPtvsMultStat->Draw("e");
  fHistPtvsMultSist->SetFillStyle(0);
  //  fHistPtvsMultSist->Draw("same e2");
  //  fHistPtvsMultStatFromSpectrum->Draw("same e");
  fHistPtvsMultSistFromSpectrum->SetFillStyle(0);
  //  fHistPtvsMultSistFromSpectrum->Draw("same e2");
  fHistPtvsMultStatMeanMacro->Draw("same e");
  fHistPtvsMultSistMeanMacro->SetFillStyle(0);
  //  fHistPtvsMultSistMeanMacro->Draw("same e2");
  fHistAvgPtDatiPubblicatiStat->SetMarkerColor(881);
  fHistAvgPtDatiPubblicatiStat->SetLineColor(881);
  fHistAvgPtDatiPubblicatiSist->SetMarkerColor(881);
  fHistAvgPtDatiPubblicatiSist->SetLineColor(881);
  fHistAvgPtDatiPubblicatiStat->Draw("same");
  fHistAvgPtDatiPubblicatiSist->SetFillStyle(0);
  fHistAvgPtDatiPubblicatiSist->Draw("same e2");
  legendPtvsMult->AddEntry(fHistPtvsMultSist, "syst. from fit", "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStat, "stat. from fit", "pel");
  legendPtvsMult->AddEntry(fHistPtvsMultSistFromSpectrum, Form("syst. from spectrum for pt>%.1f", LowPtLimitForAvgPtFSAllMult), "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStatFromSpectrum, Form("stat. from spectrum for pt>%.1f", LowPtLimitForAvgPtFSAllMult), "pel");
  legendPtvsMult->AddEntry(fHistPtvsMultSistMeanMacro, "syst. from MeanMacro", "ef");
  legendPtvsMult->AddEntry(fHistPtvsMultStatMeanMacro, "stat. from MeanMacro", "pel");
  legendPtvsMult->Draw("");

  canvasYieldErr->cd();
  StyleHisto(fHistYieldStatRelErr, 10e-5, LimSupYieldErr, Color[type], 33, titleYieldX, titleYRel,titleYRel+ " yield vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  fHistYieldStatRelErr->GetYaxis()->SetTitleOffset(1.2);
  StyleHisto(fHistYieldSistRelErr, 10e-5, LimSupYieldErr, Color[type], 27, titleYieldX, titleYRel, titleYRel+" yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistYieldSistNoExtrRelErr, 10e-5, LimSupYieldErr, 881, 27,  titleYieldX, titleYRel, titleYRel+ " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  legendYieldErr->AddEntry(fHistYieldStatRelErr, "stat.", "pl");
  //  if (!(TypeAnalysis==0 && type==8))   legendYieldErr->AddEntry(fHistYieldSistRelErr, "syst.", "pl");
  legendYieldErr->AddEntry(fHistYieldSistRelErr, "syst.", "pl");
  legendYieldErr->AddEntry(fHistYieldSistNoExtrRelErr, "syst. (no extr)", "pl");
  fHistYieldSistRelErr->GetYaxis()->SetTitleOffset(1.2);
  fHistYieldStatRelErr->Draw("e p");
  fHistYieldSistNoExtrRelErr->Draw("same p");
  fHistYieldSistRelErr->Draw("same p");
  legendYieldErr->Draw("");

  canvasPtErr->cd();
  StyleHisto(fHistPtStatRelErr, 10e-5, LimSupPtErr, 628, 33, titleYieldX, titleYRel,titleYRel+" of " + " <pt> vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtStatRelErrFromSpectrum, 10e-5, LimSupPtErr, 630, 27, titleYieldX, titleYRel, titleYRel+" of " + " <pt> vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtStatRelErrMeanMacro, 10e-5, LimSupPtErr, 881, 21,  titleYieldX, titleYRel, titleYRel+" of " + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);

  StyleHisto(fHistPtSistRelErr, 10e-5, LimSupPtErr, 418, 33, titleYieldX, titleYRel,titleYRel+" of " + " <pt> vs multiplicity", 1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtSistRelErrFromSpectrum, 10e-5, LimSupPtErr, 420, 27, titleYieldX, titleYRel, titleYRel+" of " + " <pt> vs multiplicity",  1, 0, LimSupdNdEtaAxis);
  StyleHisto(fHistPtSistRelErrMeanMacro, 10e-5, LimSupPtErr, 424, 21,  titleYieldX, titleYRel, titleYRel+" of " + " yield vs multiplicity",  1, 0, LimSupdNdEtaAxis);

  StyleHisto(fHistPtSistRelErrFit, 10e-5, LimSupPtErr, 1, 22, titleYieldX, titleYRel,titleYRel+" of " + " <pt> vs multiplicity", 1, 0, LimSupdNdEtaAxis);

  legendPtErr->AddEntry(fHistPtStatRelErr, "stat. from fit", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErr, "syst. from fit", "pl");
  legendPtErr->AddEntry(fHistPtStatRelErrFromSpectrum, "stat. from spectrum", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErrFromSpectrum, "syst. from spectrum", "pl");
  legendPtErr->AddEntry(fHistPtStatRelErrMeanMacro, "stat. from mean macro", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErrMeanMacro, "syst. from mean macro", "pl");
  legendPtErr->AddEntry(fHistPtSistRelErrFit, "syst. from fit function choice (Mean Macro)", "pl");

  fHistPtStatRelErr->Draw("e p");
  fHistPtSistRelErr->Draw("same p");
  fHistPtStatRelErrFromSpectrum->Draw("same e p");
  fHistPtSistRelErrFromSpectrum->Draw("same p");
  fHistPtStatRelErrMeanMacro->Draw("same e p");
  fHistPtSistRelErrMeanMacro->Draw("same p");
  fHistPtSistRelErrFit->Draw("same p");
  legendPtErr->Draw("");

  cout << "\n\n going to write on file " << endl;   
  if (isMeanMacro){
    for(Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      fileout->WriteTObject(hhoutYield[m]);
      fileout->WriteTObject(hhoutYieldMy[m]);
      fileout->WriteTObject(hhoutYieldRatioToMine[m]);
      fileout->WriteTObject(hhoutRelStat[m]);
      fileout->WriteTObject(hhoutRelStatMy[m]);
      fileout->WriteTObject(hhoutRelSystHigh[m]);
      fileout->WriteTObject(hhoutRelSystHighMy[m]);
      fileout->WriteTObject(hhoutRelSystLow[m]);
      fileout->WriteTObject(hhoutRelSystLowMy[m]);

      fileout->WriteTObject(hhoutAvgPt[m]);
      fileout->WriteTObject(hhoutAvgPtMy[m]);
      fileout->WriteTObject(hhoutAvgPtRatioToMine[m]);
      fileout->WriteTObject(hhoutAvgPtRelStat[m]);
      fileout->WriteTObject(hhoutAvgPtRelStatMy[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystHigh[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystHighMy[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystLow[m]);
      fileout->WriteTObject(hhoutAvgPtRelSystLowMy[m]);

      for (Int_t typefit=0; typefit<numfittipo; typefit++){
	if (nameFit[typefit]!="Levi") continue;
	fileout->WriteTObject(hhout[typefit][m]);
	//	fileout->WriteTObject(fHistAvgPtDistr[m][typefit]);
      }
    }
  }
  cout << " ok " << endl;
  fileout->WriteTObject(canvasYield);
  fileout->WriteTObject(canvasYieldErr);
  fileout->WriteTObject(canvasPtvsMult);
  fileout->WriteTObject(canvasPtErr);
  fileout->WriteTObject(fHistYieldStat);
  fileout->WriteTObject(fHistYieldSist);
  fileout->WriteTObject(fHistYieldSistNoExtr);
  fileout->WriteTObject(fHistPtvsMultStat);
  fileout->WriteTObject(fHistPtvsMultSist);
  fileout->WriteTObject(fHistPtvsMultStatFromSpectrum);
  fileout->WriteTObject(fHistPtvsMultSistFromSpectrum);
  cout << " ok " << endl;
  if (isMeanMacro){
    fileout->WriteTObject(fHistPtvsMultStatMeanMacro);
    fileout->WriteTObject(fHistPtvsMultSistMeanMacro);
  }
  fileout->WriteTObject(fHistYieldStatRelErr);
  fileout->WriteTObject(fHistYieldSistRelErr);
  fileout->WriteTObject(fHistYieldSistNoExtrRelErr);
  cout << " ok " << endl;
  for(Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
      fHistAvgPtFSDistr[m] ->GetXaxis()->SetRangeUser(LimInfPtvsMultTight, LimSupPtvsMultTight);
      //      fHistAvgPtDistr[m][0] ->GetXaxis()->SetRangeUser(LimInfPtvsMultTight, LimSupPtvsMultTight);
      //      fileout->WriteTObject( fHistAvgPtDistr[m][0]);
      fileout->WriteTObject( fHistAvgPtFSDistr[m]);
  }
  cout << " ok " << endl;
  fileout->WriteTObject(canvasMCChoiceDef);
  fileout->WriteTObject(canvasBarlowMCChoiceDef);
  fileout->WriteTObject(canvasPtSpectra);
  fileout->WriteTObject(canvasPtSpectraFit);
  fileout->WriteTObject(canvasPtSpectraFitRatio);
  fileout->WriteTObject(canvasPtSpectraFitUp);
  fileout->WriteTObject(canvasPtSpectraFitDown);
  fileout->WriteTObject(canvasPtSpectraFitHard);
  fileout->WriteTObject(canvasPtSpectraFitSoft);
  fileout->WriteTObject(canvasPtSpectraFitBis);
  fileout->WriteTObject(canvasPtSpectraRelError);
  fileout->WriteTObject(canvasPtSpectraRelErrorAll);
  fileout->WriteTObject(canvasFitResult);
  if (isNormCorrFullyComputed==1) fileout->WriteTObject(canvasNormFactor);
  cout << " ok " << endl;
  for(Int_t m=0; m<nummolt+1; m++){
    if (isppHM && MultBinning==1 && m<=1) continue;
    if (MultBinning==3 && (m==2 || m==3 || m==4)) continue;
    fileout->WriteTObject(fHistSpectrumSistAll[m]);
    //    if (isNormCorrFullyComputed)    fileout->WriteTObject(fHistSpectrumStatNotNorm[m]);
    fileout->WriteTObject(fHistSpectrumStat[m]);
    fileout->WriteTObject(fHistSpectrumStatUp[m]);
    fileout->WriteTObject(fHistSpectrumStatDown[m]);
    fileout->WriteTObject(fHistSpectrumStatHard[m]);
    fileout->WriteTObject(fHistSpectrumStatSoft[m]);
    fileout->WriteTObject(fHistSpectrumSistRelErrorMCChoice[m]);
  }
  cout << " ok " << endl;
  //  TString DirPicture = "PictureForNote/";
  TString DirPicture = PathOutPictures;
  canvasYield->SaveAs(DirPicture+"YieldvsMult"+tipo[type]+".pdf");
  canvasYieldErr->SaveAs(DirPicture+"YieldvsMultErr"+tipo[type]+".pdf");
  canvasPtErr->SaveAs(DirPicture+"PtvsMultErr"+tipo[type]+".pdf");
  canvasPtvsMult->SaveAs(DirPicture+"AvgPtvsMult"+tipo[type]+".pdf");
  canvasPtSpectra->SaveAs(DirPicture+"PtSpectra"+tipo[type]+".pdf");
  canvasPtSpectraFit->SaveAs(DirPicture+"PtSpectraFit"+tipo[type]+".pdf");
  canvasPtSpectraFitBis->SaveAs(DirPicture+"PtSpectraFitBis"+tipo[type]+".pdf");
  canvasPtSpectraRelError->SaveAs(DirPicture+"PtSpectraRelErr"+tipo[type]+".pdf");
  canvasPtSpectraRelErrorAll->SaveAs(DirPicture+"PtSpectraRelErrAll"+tipo[type]+".pdf");
  canvasFitResult->SaveAs(DirPicture+"ChiSquare"+tipo[type]+".pdf");

  fileout->Close();

  cout << "\n\e[35mStarting from the file(s):"<< endl;
  cout  << "->for default spectra: \e[39m";
  cout << PathInDef << endl;

  if (isNormCorrFullyComputed) cout << "\n\e[35mFile from where normalisation factor is taken:\e[39m " << SfileNormCorrFC << endl;

  cout << "\nWith respect to preliminaries: " << endl;
  cout << "1) The final spectra are obtained using PYTHIA8 efficiency, and are not an average of the spectra obtaiend with EPOS and PYTHIA8, becasue I couldn't compute 2D efficiency using EPOS. The relative uncertainty associated to the choice of MC is however computed from 0-100% class, starting from the file:\n" <<PathInMCChoiceDef << endl;

  cout << "numfittipoEff " << numfittipoEff << endl;
  cout << "\nI have created the file:\n" << stringout << "\n" <<endl;

}

