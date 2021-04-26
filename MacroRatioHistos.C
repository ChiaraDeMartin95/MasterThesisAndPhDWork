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


void ErrRatioCorr(TH1F* hNum, TH1F* hDenom, TH1F* hRatio, Bool_t FullCorr){
  //FullCorr == 1 means ro = 1;                                                                                                           
  //FullCorr == 0 means ro =                                                                                                               
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
    if (FullCorr==0){
      ErrC=pow(hDenom->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    else {
      ErrC=hDenom->GetBinError(b) * hNum->GetBinError(b)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    Err=sqrt(Err1+Err2-ErrC);
    hRatio->SetBinError(b,Err*hRatio->GetBinContent(b));
  }
  //  return hRatio;                                                                                                                      
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); //1.2                                                                                       
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void MacroRatioHistos(Int_t RunVar=4, Int_t AnalysisType =0, TString  CommonFileName = ""/*name common to all files */,  TString OutputName =""){

  //RunVar should be increased when you want to do a new comparison; also, OutputName And CommonFileName should be updated below!!
  TString TypeAnalysis[3] = {"Jet", "Bulk", "All"};
  TString TypeAnalysisBis[3] = {"Jet", "Bulk", "Inclusive"};
  Int_t numFiles=0;
  Int_t TypeOfComparison = 1; 

  //TypeOfComparison:
  //1) Comparison with respect to one plot (e.g., with respect to plot in file number 0, which might be a specific mult class or a specific pt,trig min). The number of the file with respect to which the ratio is computed is defined below
  //2) Comparison between a variable calculated in one way and the same variable calculated in another way -- done in CLASSES (e.g. multiplicity classes, pt,trig min classes)
  //3) Like 1, but for two different selections (e.g. I compare spectra calculated in one way (way A) to 0-100% calculated in way A, and spectra calculated in another way (way B) to 0-100% (way B) )

  if (RunVar==0){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelationRun2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3";
    numFiles=8;
  }
  else if (RunVar==1){
    CommonFileName= "FinalOutput/DATA2016/histo/AngularCorrelation2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName="hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOB";
    numFiles=4;
  }
  else if (RunVar==2){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOB_Comparison";
    numFiles=8;
    TypeOfComparison = 2; 
  }

  else if (RunVar==3){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius_Comparison";
    numFiles=12;
    TypeOfComparison = 2; 
  }
  else if (RunVar==4){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius";
    numFiles=6;
  }
  else if (RunVar==5){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius_SkipAllAssocComp";
    numFiles=8;
    TypeOfComparison = 2; 
  }
  else if (RunVar==6){
    //comparison between Xi spectra obtained with default selections, only-TOF OOB pileup rejection and +RV0 < 34 cm
    CommonFileName = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_SpectraComparison";
    numFiles=3;
  }
  else if (RunVar==7){ 
    //*****TYPE 1: comparison bertween efficiency calculated in different ways (same multiplicity class)*****
    //comparison between Xi efficiency obtained with old AODs and with new (refiltered) AODs (AOD235)
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = "FinalOutput/DATA2016/Efficiency/EfficiencyComparisonOldNewAODsXi";
    numFiles=12; //2 x 6 mult classes
    TypeOfComparison = 2; 
  }
  else if (RunVar==8){
    //comparison norm factor for K0s obtained with 200M INT7 events to the one obtained with 300M INT7 events
    CommonFileName = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid";
    OutputName = "FinalOutput/DATA2016/MCClosureCheckComparisonNormFactor";
    numFiles=12; //2 x 6 mult classes
    TypeOfComparison = 2; 
  }
  else if (RunVar==9){
    //multiplicity comparison norm factor for K0s obtained with 200M INT7 events to the one obtained with 300M INT7 events
    CommonFileName = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid";
    OutputName = "FinalOutput/DATA2016/MCClosureCheckMultComparisonNormFactor";
    numFiles=6; //6 mult classes
  }
  else if (RunVar==10){
    //comparison norm factor for K0s obtained with 200M INT7 events to the one obtained with 300M INT7 events
    CommonFileName = "CompareYieldDifferentCollisions";
    OutputName = "CompareYieldHMK0s_All_16k_vs_AllBut16k";
    numFiles=2; //
  }
  else if (RunVar==11){
    //comparison between K0s efficiency from 18g4_extra to K0s efficiency from 1617MC
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = "CompareK0sEfficiency_18g4_extra_vs_1617MC";
    numFiles=12; //2+6 mult
    TypeOfComparison = 2; 
  }
  else if (RunVar==12){
    //comparison MB yield obtained with 18g4_extra eff (pt) to those obtained with 18g4_extra eff (pt, eta)
    CommonFileName = "CompareYieldDifferentCollisions";
    OutputName = "CompareYieldMB_EtaEff_vs_NoEtaEff_18g4extra_"+ TypeAnalysis[AnalysisType];
    numFiles=2; //
  }
  else if (RunVar==13){
    //comparison HM yields vs mult obtained from all but 16k data; comparison is between results corrected by eff 2019h11c_extra and those corrected by efficiency from 2019h11a,b,c_extra
    CommonFileName = "FinalOutput/DATA2016/SystematicAnalysisAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_"+  TypeAnalysis[AnalysisType]+"Data_PtMin3.0_IsEtaEff";
    OutputName = "CompareYieldHMK0s_AllBut16k_Eff2019h11c_vs_2019h11abc_"+ TypeAnalysis[AnalysisType];
    numFiles=2; //
  }
  else if (RunVar==14){
    //comparison between 5 TeV and 13 TeV K0s efficiency
    CommonFileName ="FinalOutput/DATA2016/Efficiency/Efficiency"; 
    OutputName = "CompareK0sEfficiencyAcrossEnergies";
    numFiles=12; //2*6 mult
    TypeOfComparison = 2; 
  }
  else if (RunVar==15){
    //comparison K0s at 13 TeV vs K0s at 5 TeV (yields)
    CommonFileName = "FinalOutput/DATA2016/";
    OutputName = "CompareYieldK0s_5TeV_vs_13TeV_"+ TypeAnalysis[AnalysisType];
    numFiles=2; //
  }
  else if (RunVar==16){
    //comparison K0s at 13 TeV vs K0s at 5 TeV (spectra)
    CommonFileName = "FinalOutput/DATA2016/";
    OutputName = "CompareSpectraK0s_5TeV_vs_13TeV_"+ TypeAnalysis[AnalysisType];
    numFiles=12; //
    TypeOfComparison=3;
  }
  else if (RunVar==17){
    //comparison K0s purity at 13 TeV vs K0s at 5 TeV
    CommonFileName = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_";
    OutputName = "ComparePurityK0s_5TeV_vs_13TeV";
    numFiles=12; //
    TypeOfComparison=3;
  }
  else if (RunVar==18){
    //comparison K0s efficiency at 5 TeV CENTwoSDD vs FAST
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hK0s_";
    OutputName = "CompareFASTCENTwoSDDeff_K0s_5TeV_vs_13TeV";
    numFiles=12; //
    TypeOfComparison = 2; 
  }
  else if (RunVar==19){
    //COMPARE EFFICIENCIES ACROSS MULTIPLICITIES
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_hK0s_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    OutputName = "CompareEfficiencyAcrossMult_2018f1_extra_hK0s_y0.5_AllAssoc_PtMin0.15";
    numFiles=6; //
  }
  else if (RunVar==20){
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_";
    OutputName = "CompareEfficiencyAcrossMult_And_OldNewAODs";
    numFiles = 12;
    TypeOfComparison=3;
  }
  else if (RunVar==21){
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_";
    OutputName = "CompareEfficiency_OldNewAODs";
    numFiles = 12;
    TypeOfComparison=2;
  }

  //this macros superimpose in a same canvas pad same histos found in different files. A loop over the files is done. The ratio of the histos to the histo found in the first file is performed. You can choose if histos have to be considered fully correlated, uncorrelated, or if the histos are obtaiend from a subsample of the data used to obtain the histo found in the first file (i.e. partial correlation)

  TString VarName[numFiles] = {""}; 
  TString NameHisto[numFiles] = {""}; 
  TString histoName= "fHistQA6";
  if (RunVar==6) histoName = "histo_S";
  else   if (RunVar==7 || RunVar==11 || RunVar==14 || RunVar==18 || RunVar==19 || RunVar==20 || RunVar==21) histoName = "fHistV0EfficiencyPtBins";
  else   if (RunVar==8 || RunVar==9) histoName = "SpectrumRatio" +  TypeAnalysis[AnalysisType]+"_m"; //All -> Jet, Bulk
  else if (RunVar==10 || RunVar==12) histoName = "histoYieldComparison";
  else if (RunVar==13) histoName = "fHistYieldvsErrSoloStat";
  else if (RunVar==14) histoName = "fHistV0EfficiencyPtBins";
  else if (RunVar==15 || RunVar==16) histoName = "";
  else if (RunVar==17) histoName = "histo_SSB";

  Int_t numDef=3;
  if (RunVar>=7) numDef=0;
  Float_t num=0;
  Int_t numEff=0;

  Float_t LimSupPol0=8;
  Float_t LimInfPol0 = 0;
  if (RunVar==10 || RunVar==13) LimSupPol0 = 45;
  else if (RunVar==15) LimSupPol0 = 30;
  if (RunVar==13) LimInfPol0 = 25;
  TF1 * pol0At1 = new TF1("pol0", "pol0", LimInfPol0, LimSupPol0);
  pol0At1->SetParameter(0,1);
  pol0At1->SetLineColor(1);
  pol0At1->SetLineWidth(0.5);
  TF1 * pol1[numFiles]; 

  //histo style selections
  Float_t Low=10e-8;
  Float_t Up=0.003;
  Float_t LowRatio=10e-8;
  Float_t UpRatio=0.6;
  Int_t color[10]={1,402 , 628, /*905, 881,*/601, 867, 418, 905, 881 };
  Int_t style =33;
  TString titleX = "Multiplicity class";
  TString titleY="#hXi events / #INT7 events";
  TString title="Fraction of INT7 events containing trigger particle and Xi";
  TString titleRatio="Ratio to #it{p}_{T}^{trigg} > 3 GeV/#it{c}";
  TString Smolt[6] = {"_0-5", "_5-10", "_10-30", "_30-50", "_50-100", "__all"};  
  TString SmoltBis[6] = {" 0-5%", " 5-10%", " 10-30%", " 30-50%", " 50-100%", " 0-100%"};  

  if (RunVar==2){
    titleRatio = "Ratio of only TOF to (SPD+TOF) pileup rej";
    LowRatio = 0.5;
    UpRatio = 1;
  }

  else  if (RunVar==3){
    titleRatio = "Ratio of Run3 to Run 2 selections";
    LowRatio = 0.5;
    UpRatio = 1;
  }

  else  if (RunVar==5){
    titleRatio = "Ratio of SkipAssoc to AllAssoc with Run3Sel";
    LowRatio = 0.8;
    UpRatio = 1.3;
  }
  else  if (RunVar==6){
    titleRatio = "Ratio to default Run2Sel";
    Low = 10e-6;
    Up = 5*10e-3;
    LowRatio = 0;
    UpRatio = 1.5;
    titleY = "1/N_{trigg} dN^{#Xi}/dp_{T}";
    titleX = "p_{T} (GeV/c)";
    title = "Raw Xi spectrum 0-100%, p_{T}^{trigg} > 3 GeV/c";
  }
  else  if (RunVar==7){
    titleRatio = "Ratio to new AODs";
    Low = 0;
    Up = 0.1;//0.3
    LowRatio = 0.8;
    UpRatio = 1;
    titleY = "Xi efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "Xi efficiency";
  }
  else  if (RunVar==8){
    titleRatio = "Ratio to 300M events (new) norm factor";
    Low = 0.9;
    Up = 1.2;
    LowRatio = 0.95;
    UpRatio = 1.1;
    titleY = "K0s norm factor";
    titleX = "p_{T} (GeV/c)";
    title = "K0s norm factor";
  }
  else  if (RunVar==9){
    titleRatio = "Ratio to 0-100% norm factor";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s norm factor";
    titleX = "p_{T} (GeV/c)";
    title = "K0s norm factor";
  }
  else  if (RunVar==10){
    titleRatio = "Ratio to 16k only";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==11){
    titleRatio = "Ratio to 1617MC";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s efficiency";
    titleX = "p_{T}";
    title = "K0s efficiency";
  }
  else  if (RunVar==12){
    titleRatio = "Ratio to non-eta-dependent eff";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==13){
    titleRatio = "Ratio to 2019h11c_extra eff";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==14){
    titleRatio = "Ratio to 13 TeV efficiency";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s efficiency";
    titleX = "p_{T}";
    title = "K0s efficiency";
  }
  else  if (RunVar==15){
    titleRatio = "Ratio to 13 TeV";
    if (AnalysisType==0){
      Low = 0.02;
      Up = 0.04;
    }
    else {
      Low = 0;
      Up = 0.3;
    }
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==16){
    titleRatio = "Ratio to 0-100%";
    if (AnalysisType==0){
      Low = 0.0005;
      Up = 0.02;
    }
    else {
      Low = 0.00005;
      Up = 0.3;
    }
    LowRatio = 0;
    UpRatio = 2;
    titleY = "1/N_{trigg} dN/dp_{T} (1/(GeV/c))";
    titleX = "p_{T} (GeV/c)";
    title = "K0s yield";
  }
  else  if (RunVar==17){
    titleRatio = "Ratio to 0-100%";
    Low = 0.95;
    Up = 1;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "Purity";
    titleX = "p_{T} (GeV/c)";
    title = "K0s purity";
  }
  else  if (RunVar==18){
    titleRatio = "Ratio to CENTwoSDD";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "Efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "K0s efficiency";
  }
  else  if (RunVar==19 || RunVar==20 || RunVar==21){
    titleRatio = "Ratio to 0-100%";
    if (RunVar==21)     titleRatio = "Ratio to old AODs";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.9;
    UpRatio = 1.1;
    titleY = "K0s efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "K0s efficiency";
  }

  TLegend * legend = new TLegend (0.6, 0.7, 0.9, 0.9);
  if (RunVar<6) legend->SetHeader("#it{p}_{T}^{trigg} > ");
  TString LegendName[numFiles]={""};

  TH1F * histo[numFiles];
  TH1F * hdummy[numFiles];
  TH1F * histoRatio[numFiles];

  TString InputName="";
  TFile * InputFile;
  TString OutputNameRoot= OutputName +".root";
  TString OutputNamepdf= OutputName +".pdf";
  TFile * OutputFile= new TFile (OutputNameRoot, "RECREATE");

  gStyle->SetOptStat(0);
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvas->Divide(2,1);

  cout << "TypeOfComparison: " <<     TypeOfComparison<< endl;

  for (Int_t i=0; i<numFiles; i++){
    cout << histoName << endl;
    NameHisto[i] = "";
    num = numDef+i;
    numEff= num;

    if (RunVar==9 || RunVar==19) numEff = numFiles-num-1; 

    if (TypeOfComparison==2){
      if (i<numFiles/2)      numEff = num;
      else if (i>=numFiles/2) numEff=num-numFiles/2;
    }
    else if (TypeOfComparison==3){
      if (i<numFiles/2)      numEff = numFiles/2-num-1;
      else if (i>=numFiles/2) numEff=numFiles-num-1;
    }

    cout << "\nmulteplicity " << numEff << endl;
    if (RunVar==0 || RunVar==1 || RunVar==4)    VarName[i] = Form("PtMin%.1f", num);
    else if (RunVar==2){
      if (i>=0 && i<numFiles/2) VarName[i] = Form("Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
      else if (i>=numFiles/2) VarName[i] =  Form("2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==3){
      if (i>=0 && i<numFiles/2) VarName[i] = Form("Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
      else if (i>=numFiles/2) VarName[i] =  Form("2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==5){
      if (i<=numFiles/2) VarName[i]+= "_AllAssoc";
      VarName[i] += Form("_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==6){
      if (i==0) VarName[i] = "Run2DataRed_MECorr_hXi_Xi_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      if (i==1) VarName[i] = "2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      if (i==2) VarName[i] = "2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
    }
    else if (RunVar==7){
      if (i>=0 && i<numFiles/2) VarName[i] = "161718_MD_EtaEff_LowPtTrig_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      //"161718_MD_EtaEff_PtTrig3_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      //"161718_MD_EtaEff_LowPtTrig_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      else if (i>=numFiles/2) VarName[i] = "AllMC_hXi_PtTrigMax2.5_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      //"AllMC_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      //"AllMC_hXi_PtTrigMax2.5_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==8){
      if (i>=0 && i<numFiles/2) VarName[i] = "_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
      else if (i>=numFiles/2) VarName[i] = "_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
      NameHisto[i] = Form("%i",numEff);
    }
    else if (RunVar==9){
      VarName[i] = "_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
      NameHisto[i] = Form("%i",numEff);
    }
    else if (RunVar==10){
      if (i==0)   VarName[i] = "_HM16kAll.root";
      else    VarName[i] = "All.root";
      NameHisto[i] = "";
    }
    else if (RunVar==11){
      if (i>=0 && i<numFiles/2) VarName[i] = "1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      else if (i>=numFiles/2) VarName[i] = "2018g4_extra_EtaEff_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==12){
      if (i==0)   VarName[i] = "_MBEff18g4extra_NoEtaEff" +  TypeAnalysis[AnalysisType]+".root";
      else    VarName[i] =  TypeAnalysis[AnalysisType]+".root";
      NameHisto[i] = "";
    }
    else if (RunVar==13){
      if (i==0)   VarName[i] = ".root";
      else    VarName[i] = "_AllEff2019.root";
      NameHisto[i] = "";
    }
    else if (RunVar==14){
      if (i>=0 && i<numFiles/2)   VarName[i] = "1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      else if (i>=numFiles/2)   VarName[i] = "17pq_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==15){
      if (i==0)  {
	VarName[i] = "PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_" +  TypeAnalysisBis[AnalysisType]+ ".root";
	NameHisto[i] = "fHistYieldStat";
      }
      else {
	VarName[i] = "SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_" +  TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	NameHisto[i] = "fHistYieldvsErrSoloStat";
      }
    }
    else if (RunVar==16){
      if (i<numFiles/2)  {
	VarName[i] = "PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_" +  TypeAnalysisBis[AnalysisType]+ ".root";
	NameHisto[i] = "fHistSpectrum"+Smolt[numEff];
      }
      else {
	VarName[i] = "SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_" +  TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	NameHisto[i] = Form("fHistSpectrumPart_m%i_syst0",numEff );
      }
      cout << " molt " << Smolt[numEff] << endl;
    }
    else if (RunVar==17){
      if (i<numFiles/2)  {
	VarName[i] = Form("PtBinning1_1617_hK0s_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt%i_sysT0_sysV00_Sys0_PtMin3.0.root", numEff);
	NameHisto[i] = "";
      }
      else {
	VarName[i] = Form("PtBinning1_17pq_hK0s_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt%i_sysT0_sysV00_Sys0_PtMin3.0.root", numEff);
	NameHisto[i] = "";
      }
    }
    else if (RunVar==18){
      if (i<numFiles/2)  {
	VarName[i] = "CENTwoSDD_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
      else {
	VarName[i] = "FAST_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==19){
      VarName[i] = "";
      NameHisto[i] = Smolt[numEff];
      cout << histoName << endl;
      cout << NameHisto[i] << endl;
    }
    else if (RunVar==20){
      if (i<numFiles/2)  {
	VarName[i] = "hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	NameHisto[i] = Smolt[numEff];
      }
      else {
	VarName[i] = "AOD235_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	NameHisto[i] = Smolt[numEff];
      }
    }
    else if (RunVar==21){
      if (i<numFiles/2)  {
	VarName[i] = "hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
      else {
	VarName[i] = "AOD235_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
      NameHisto[i] = Smolt[numEff];
    }

    InputName = CommonFileName + VarName[i];
    //    if (RunVar==4) InputName += "_IsOnlypiKpemu";
    if (RunVar<6) InputName+="_IsEstimateRun3.root";
    cout << " loop n. " << i << " file name: " << InputName << endl;
    InputFile = new TFile (InputName, "");
    if (!InputFile) return;
    histo[i] = (TH1F*)InputFile->Get(histoName + NameHisto[i]);
    cout << " histo name " << histoName << NameHisto[i]<< endl;
    if (!histo[i]) {cout << "histogram is not there: " << histoName << NameHisto[i] << endl; return;}
    histo[i] ->SetName(Form("histoName%i", i));
    histo[i]->Sumw2();

    canvas->cd(1);
    if (RunVar<7 || RunVar==16)    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    if (TypeOfComparison==2 && i>=numFiles/2) style = 27;
    if (TypeOfComparison==3 && i>=numFiles/2) style = 27;
    if (RunVar==19) style = 33;
    if (RunVar==10 || RunVar==12 || RunVar==13 || RunVar==15) color[numEff-numDef]=1;
    StyleHistoYield(histo[i], Low, Up, color[numEff-numDef], style, titleX, titleY, title, 1, 1.2, 1.4);
    if (RunVar<=6) LegendName[i] = Form("%.1f GeV/#it{c}", num);
    if (RunVar==2){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}, only TOF PU rej", (float)numEff);
    }
    else if (RunVar==3){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}, Run3 sel", (float)numEff);
    }
    else if (RunVar==5){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}, AllAsso", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
    }
    else if (RunVar==6){
      if (i==0)  LegendName[i] = " Run2 sel";
      if (i==1)  LegendName[i] = " only-TOF OOBpu rej";
      if (i==2)  LegendName[i] = " Run3 sel";
    }
    else if (RunVar==7){
      if (i<numFiles/2)    LegendName[i] = "New AODs " + SmoltBis[numEff];
      else    LegendName[i] = "Old AODs " + SmoltBis[numEff];
    }
    else if (RunVar==8){
      if (i<numFiles/2)    LegendName[i] = "New " + SmoltBis[numEff];
      else    LegendName[i] = "Old " + SmoltBis[numEff];
    }
    else if (RunVar==9){
      LegendName[i] = "New " + SmoltBis[numEff];
    }
    else if (RunVar==10){
      if (i==0)   LegendName[i] = "16k ";
      else if (i==1)    LegendName[i] = "AllBut16k ";
    }
    else if (RunVar==11){
      if (i<numFiles/2)    LegendName[i] = "1617MC " + SmoltBis[numEff];
      else    LegendName[i] = "18g4_extra " + SmoltBis[numEff];
    }
    else if (RunVar==12){
      if (i==0)   LegendName[i] = "18g4_extra "; //denominator
      else if (i==1)    LegendName[i] = "18g4_extra -- eta-dep eff"; //numerator
    }
    else if (RunVar==13){
      if (i==0)   LegendName[i] = "Eff: 2019h11c_extra "; //denominator
      else if (i==1)    LegendName[i] = "Eff: 2019h11abc_extra"; //numerator
    }
    else if (RunVar==14){
      if (i<numFiles/2)   LegendName[i] = "Eff: 13 TeV "+ SmoltBis[numEff]; //denominator
      else   LegendName[i] = "Eff: 5 TeV" + SmoltBis[numEff]; //numerator
    }
    else if (RunVar==15){
      if (i==0)   LegendName[i] = "13 TeV "; //denominator
      else if (i==1)    LegendName[i] = "5 TeV "; //numerator
    }
    else if (RunVar==16){
      if (i<numFiles/2)    LegendName[i] = "13 TeV" + SmoltBis[numEff];
      else    LegendName[i] = "5 TeV " + SmoltBis[numEff];
    }
    else if (RunVar==17){
      if (i<numFiles/2)    LegendName[i] = "13 TeV" + SmoltBis[numEff];
      else    LegendName[i] = "5 TeV " + SmoltBis[numEff];
    }
    else if (RunVar==18){
      if (i<numFiles/2)    LegendName[i] = "CENTwoSDD " + SmoltBis[numEff];
      else    LegendName[i] = "FAST " + SmoltBis[numEff];
    }
    else if (RunVar==19){
      LegendName[i] = SmoltBis[numEff];
    }
    else if (RunVar==20 || RunVar==21){
      if (i<numFiles/2)    LegendName[i] = "Old AODs " + SmoltBis[numEff];
      else    LegendName[i] = "AOD 235 " + SmoltBis[numEff];
    }

    legend->AddEntry(histo[i], LegendName[i], "pl");
    cout << " n bins " <<     histo[i]->GetNbinsX() << endl;
    pol1[i] = new TF1(Form("pol1_%i", i), "pol1", 0,30);
    if (TypeOfComparison==3){
      if (i<numFiles/2){
	//	if (Smolt[numEff] == "_10-30"|| Smolt[numEff] == "__all" || Smolt[numEff] == "_30-50" )     histo[i]->Draw("same");
	//	if (Smolt[numEff] == "_10-30"|| Smolt[numEff] == "__all" )     histo[i]->Draw("same");
	//	if (Smolt[numEff] == "_50-100" || Smolt[numEff] == "_30-50"  )     histo[i]->Draw("same");
      }
      else {
	//	if (Smolt[numEff] == "_5-10"|| Smolt[numEff] == "__all" || Smolt[numEff] == "_10-30" )     histo[i]->Draw("same");
	//	if (Smolt[numEff] == "_50-100" || Smolt[numEff] == "_30-50" || Smolt[numEff] == "_10-30")     histo[i]->Draw("same");
      }
      histo[i]->Draw("same");
    }
    else histo[i]->Draw("same");

    if (RunVar==15 && AnalysisType!=0) {
      if (i==1) pol1[i]->SetLineStyle(8);
      pol1[i]->SetLineColor(1);
      pol1[i]->SetLineWidth(0.3);
      histo[i]->Fit(pol1[i], "R+");
    }
    if (i==numFiles-1) legend->Draw("");

    cout << " first canvas ok " << endl;
    canvas->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    histoRatio[i] = (TH1F*) histo[i]->Clone(histoName + "_Ratio");

    if (TypeOfComparison==2){
      if (i>=numFiles/2){
	histoRatio[i]->Divide(histo[i-numFiles/2]);
	ErrRatioCorr(histo[i], histo[i-numFiles/2], histoRatio[i], 0);
      }
    }
    else if (TypeOfComparison==3){
      if (i<numFiles/2){
	cout << " I divide " <<   histo[i]->GetName() << " by " << histo[0]->GetName() << endl;
	histoRatio[i]->Divide(histo[0]);
	ErrRatioCorr(histo[i], histo[0], histoRatio[i], 0);
      }
      else {
	cout << " I divide " <<   histo[i]->GetName() << " by " << histo[numFiles/2]->GetName() << endl;
	histoRatio[i]->Divide(histo[numFiles/2]);
	//	histoRatio[i]->Divide(histo[0]);
	ErrRatioCorr(histo[i], histo[numFiles/2], histoRatio[i], 0);
      }
      for (Int_t b=1; b<= histoRatio[i]->GetNbinsX() ; b++){
	cout << histoRatio[i]->GetBinContent(b) << " ";
      }
    }
    else {
      if (i!=0){
	histoRatio[i]->Divide(histo[0]); //no corr
	if (RunVar==12 || RunVar==13)	ErrRatioCorr(histo[i], histo[0], histoRatio[i], 1); //full corr
	else if (RunVar!=15) ErrRatioCorr(histo[i], histo[0], histoRatio[i], 0); //partial corr
      }
    }

    cout << " ok up to here " << endl;
    StyleHistoYield(histoRatio[i], LowRatio, UpRatio, color[numEff-numDef], style, titleX, "Ratio", titleRatio, 1, 1.2, 1.4);

    //old    if (RunVar==2 || RunVar==3 || RunVar==5 || RunVar==7 || RunVar==8 || RunVar==9 || RunVar==11 || RunVar==14 || RunVar==18 || RunVar==19) {
    if (TypeOfComparison==2) {
      if (i>=numFiles/2) {
	histoRatio[i]->Draw("same");
	if (RunVar==8 || RunVar==14 || RunVar==18) pol0At1->Draw("same");
      }
    }
    else if (TypeOfComparison==3){
      if (numEff!=5 && numEff!=11) histoRatio[i]->Draw("same");
      if (i==numFiles-1) {
	legend->Draw("");
	pol0At1->Draw("same");
      }
    }
    else {
      if (i!=0) histoRatio[i]->Draw("same");
      if (i==numFiles-1) legend->Draw("");
      if (i!=0 && (RunVar==10 || RunVar==13 || RunVar==15 || RunVar==19)) pol0At1->Draw("same");
    }

    hdummy[i]= new TH1F (Form("hdummy%i", i), Form("hdummy%i", i), 1000, 0, 30);
    if (RunVar==15 && AnalysisType!=0) {
      hdummy[i]->Add(pol1[i]);
      if (i==1)  {
	hdummy[1]->Divide(hdummy[0]);
	hdummy[1]->Draw("same");
      }
    }
  }
  canvas->SaveAs(OutputNamepdf);
  OutputFile->WriteTObject(canvas);
  OutputFile->Close();

  cout << "I produced the output file " << OutputNameRoot << " and " << OutputNamepdf << endl;

}
