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
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>
#include </data/dataalice/cdemart/ALICE_analysis_tutorial/Macros/ErrRatioCorr.C>
#include "Macros/constants.h"

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
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

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin){
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}
//take spectra in input (MB + HM 13 TeV)
//produces ratio of spectra wrt 0-100% multiplciity class

//takes Xi/K0s ratios in input (MB + HM 13 TeV) 
//produces doubles ratios wrt 0-100% multiplciity class

Int_t nummoltMax = 5;
const Int_t numtipo=10;
const Int_t numregions=3;

TString TitleYPtRatio = "N_{#Xi}/N_{K^{0}_{S}}";
TString TitleYPtDRatio = "N_{#Xi}/N_{K^{0}_{S}}_{mult} / N_{#Xi}/N_{K^{0}_{S}}_{0-100%}_{0-100%}";
TString titleX=  "#it{p}_{T} (GeV/#it{c})";
TString titleY=  "1/(#Delta#it{#eta} #Delta#it{#varphi}) 1/#it{N}_{trigg} d#it{N}/d#it{p}_{T} [(GeV/#it{c})^{-1}]";
TString title = "Multiplicity class "; 
TString titleYSpectraRatio = "(d#it{N}^{#Xi}/d#it{p}_{T}) / (d#it{N}^{K^{0}_{S}}/d#it{p}_{T})";
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};
TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString SRegionTypeBis[3] = {"Jet", "OOj", "Full"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};

TString sRegion1K0s[3]={"|#Delta#it{#eta}| < 0.86, |#Delta#it{#varphi}| < 1.1", "0.86 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1Xi[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.1", "0.75 < |#Delta#it{#eta}| < 1.2, 0.96 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1K0sGen[3]={"|#Delta#it{#eta}| < 0.96, |#Delta#it{#varphi}| < 1.09", "0.86 < |#Delta#it{#eta}| < 1.2, 0.97 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1XiGen[3]={"|#Delta#it{#eta}| < 0.96, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.97 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"}; //like the ones used for Preliminaries in 202
TString sRegion1[3][2] ={""};

Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
TString tipo[numtipo]={"K0s", "Xi"};
TString Stipo[numtipo]={"K^{0}_{S}", "#Xi"};
TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
TString Srap[2] = {"_Eta0.8", "_y0.5"};
TString SSkipAssoc[2] = {"_AllAssoc", ""};

void PlotSpectraRatios(Bool_t ispp5TeV = 1, Float_t PtTrigMin =3, Float_t PtTrigMax=15, TString Dir="FinalOutput",TString year0="2016",Int_t ChosenMult=8 /*0-100%*/){

  Bool_t  MultOK=0;
  if (ispp5TeV) ChosenMult = 2;
  gStyle->SetOptStat(0);

  //latexnames
  TString SmoltBis[nummoltMax+1+3]={"0#minus0.01", "0.01#minus0.05", "0.05#minus0.1", "0#minus5", "5#minus10", "10#minus30", "30#minus50", "50#minus100", "0#minus100"};
  TString Smolt[nummoltMax+1+3]={ "0-0.01", "0.01-0.05", "0.05-0.1", "0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString SmoltLegend[nummoltMax+1+3]={"0-0.01%", "0.01-0.05%", "0.05-0.1%", "0-5%", "5-10%", "10-30%", "30-50%", "50-100%", "0-100%"};

  TString SmoltBis5TeV[nummoltMax+1+3]={"0#minus10", "10#minus100","0#minus100"};
  TString Smolt5TeV[nummoltMax+1+3]={ "0-10", "10-100", "_all"};
  TString SmoltLegend5TeV[nummoltMax+1+3]={"0-10%", "10-100%", "0-100%"};

  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  Float_t YoffsetSpectra[2] = {1.6, 1.65};
  Float_t YoffsetSpectraRatio[2] = {0.8, 0.8};

  Int_t ColorMult[nummoltMax+1+3] ={634, 628, 797,815,418, 429, 867,601,1};
  Int_t ColorMult5TeV[nummoltMax+1+3] ={628, 797, 1,815,418, 429, 867,601,1};
  Int_t ColorMultGenOnTheFly[nummolt+1] ={634, 628, 807, 797, 815, 418, 429, 867, 856, 601, 1};
  Int_t Color[3] ={628, 418, 600};

  Float_t size[nummoltMax+1+3] ={2, 2, 2.8, 2.5, 2.8, 2, 2.8, 2.5, 2};
  Float_t size5TeV[nummoltMax+1+3] ={2, 2, 2, 2.5, 2.8, 2, 2.8, 2.5, 2};
  Int_t MarkerMult[nummoltMax+1+3] ={20, 21, 33, 34, 29, 24, 27, 28, 25};
  Int_t MarkerMult5TeV[nummoltMax+1+3] ={20, 21, 25, 34, 29, 24, 27, 28, 25};

  Float_t ScaleFactor[nummoltMax+1+3]={};
  TString sScaleFactor[nummoltMax+1+3]={};

  Float_t ScaleFactorK0s[nummoltMax+1+3]={1024, 512, 256, 128,64, 32,16,8,1};
  TString sScaleFactorK0s[nummoltMax+1+3]={" (x2^{10})", " (x2^{9})", " (x2^{8})"," (x2^{7})"," (x2^{6})"," (x2^{5})"," (x2^{4})"," (x2^{3})",""};

  Float_t ScaleFactorXi[nummoltMax+1+3]={2048, 1024, 512, 256, 128,64, 32,16,1};
  TString sScaleFactorXi[nummoltMax+1+3]={" (x2^{11})", " (x2^{10})", " (x2^{9})", " (x2^{8})"," (x2^{7})"," (x2^{6})"," (x2^{5})"," (x2^{4})",""};

  Float_t ScaleFactorK0s5TeV[nummoltMax+1+3]={8,4,1};
  TString sScaleFactorK0s5TeV[nummoltMax+1+3]={" (x2^{3})"," (x2^{3})",""};

  Float_t ScaleFactorXi5TeV[nummoltMax+1+3]={8,4,1};
  TString sScaleFactorXi5TeV[nummoltMax+1+3]={" (x2^{3})"," (x2^{3})",""};
  
  //filein
  TString  PathInYieldMB[2][3];
  TString  PathInYieldHM[2][3]; 

  TFile * fileInYield[2];
  TFile * fileInYieldMB[2][3];
  TFile * fileInYieldHM[2][3];
  TFile * fileInYieldFit[2];
  //  TFile * fileInDPhiProj[2];
  for (Int_t type=0; type<2; type++){
    for (Int_t iregion=0; iregion<3; iregion++){

      PathInYieldMB[type][iregion] += "FinalOutput/DATA2016/PtSpectraBisNew";
      if (ispp5TeV){
	if (type==0)  PathInYieldMB[type][iregion] += "_pp5TeV_PtBinning1_17pq_hK0s_K0s";
	else  PathInYieldMB[type][iregion] += "_pp5TeV_17pq_hXi_Xi";
      }
      else {
	if (type==0)  PathInYieldMB[type][iregion] += "_PtBinning1_1617_AOD234_hK0s_K0s";
	else  PathInYieldMB[type][iregion] += "_161718Full_AOD234_hXi_Xi";
      }
      PathInYieldMB[type][iregion] += "_Eta0.8_AllAssoc_PtMin3.0_";
      PathInYieldMB[type][iregion] += RegionType[iregion];
      PathInYieldMB[type][iregion] += "_isNormCorrFullyComputed_isErrorAssumedPtCorr_ChangesIncluded_isdNdEtaTriggered";
      if (ispp5TeV)       PathInYieldMB[type][iregion] += "_MultBinning3";
      if (type==0) PathInYieldMB[type][iregion] += "_EffCorr_isWingsCorrectionAppliedNew";
      if (type==0 && !ispp5TeV) PathInYieldMB[type][iregion] +="_MatBudgetCorr.root";
      else PathInYieldMB[type][iregion] += "_MatBudgetCorrFAST.root";

      PathInYieldHM[type][iregion] += "FinalOutput/DATA2016/PtSpectraBisNew";
      if (type==0)  PathInYieldHM[type][iregion] += "_pp13TeVHM_PtBinning1_AllhK0sHM_RedNo16k_K0s";
      else  PathInYieldHM[type][iregion] += "_pp13TeVHM_161718_HM_hXi_WithFlat16k_No18p_Xi";
      PathInYieldHM[type][iregion] += "_Eta0.8_AllAssoc_PtMin3.0_";
      PathInYieldHM[type][iregion] += RegionType[iregion];
      PathInYieldHM[type][iregion] += "_isNormCorrFullyComputed_isErrorAssumedPtCorr_ChangesIncluded_isdNdEtaTriggered_MultBinning1";
      if (type==0) PathInYieldHM[type][iregion] += "_EffCorr_isWingsCorrectionAppliedNew_MatBudgetCorrFAST.root";
      else PathInYieldHM[type][iregion] += "_MatBudgetCorrFAST.root";

      fileInYieldMB[type][iregion] = new TFile(PathInYieldMB[type][iregion], "");
      fileInYieldHM[type][iregion] = new TFile(PathInYieldHM[type][iregion], "");
    }
  }

  //fileout name
  TString stringout;
  TString stringoutpdf;
  stringout = Dir+"/DATA"+year0+"/PlotSpectraRatios";
  stringout +=Srap[0];
  stringout +=SSkipAssoc[0];
  stringout+=   Form("_PtMin%.1f", PtTrigMin);
  stringout+="_isErrorAssumedPtCorr";;
  if (ispp5TeV) stringout += "_pp5TeV";
  stringoutpdf = stringout;
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  //canvases
  TCanvas *canvasPtSpectra[3][2];
  TCanvas *canvasPtSpectraMultRatio[3][2];
  TCanvas *canvasPtSpectraRatio= new TCanvas("canvasPtSpectraRatio", "canvasPtSpectraRatio", 1300, 800);
  canvasPtSpectraRatio->Divide(3,2);
  TCanvas *canvasPtSpectraRatioAllMult= new TCanvas("canvasPtSpectraRatioAllMult", "canvasPtSpectraRatioAllMult", 1300, 400);
  canvasPtSpectraRatioAllMult->Divide(3,1);
  TCanvas *canvasPtSpectraK0s= new TCanvas("canvasPtSpectraK0s", "canvasPtSpectraK0s", 1300, 800);
  canvasPtSpectraK0s->Divide(3,2);
  for (Int_t type=0; type<2; type++){
    for (Int_t iregion=0; iregion<3; iregion++){
      canvasPtSpectra[iregion][type] = new TCanvas("canvasPtSpectra"+RegionType[iregion]+tipo[type], "canvasPtSpectra"+RegionType[iregion]+tipo[type], 900, 1000);
      canvasPtSpectraMultRatio[iregion][type] = new TCanvas("canvasPtSpectraMultRatio"+RegionType[iregion]+tipo[type], "canvasPtSpectraMultRatio"+RegionType[iregion]+tipo[type], 1300, 800);
    }
  }

  //limits spectra and spectra ratios
  Float_t LimSupMultRatio[2]= {3.1,3.9};
  LimSupMultRatio[0] = 8.-10e-4;
  LimSupMultRatio[1] = 8.-10e-4;
  if (ispp5TeV) {
    LimSupMultRatio[0] = 4.-10e-4;
    LimSupMultRatio[1] = 4.-10e-4;
  }

  Float_t LimSupSpectra[2][3] ={{19999.999, 19999.999, 19999.999},{19999.999, 19999.999, 19999.999}};
  Float_t LimInfSpectra[2][3] ={{5*10e-8, 5*10e-8, 5*10e-8},{1.01*10e-8, 1.01*10e-8, 1.01*10e-8}};
  if (ispp5TeV){
    for (Int_t i=0; i<3; i++){
      LimInfSpectra[0][i] = 5*10e-8;
      LimSupSpectra[0][i] = 99.999;
      LimInfSpectra[1][i] = 1.01*10e-8;
      LimSupSpectra[1][i] = 2.9999;
    }
  }

  Int_t bcorr=0;
  TH1F *fHistSpectrumStat[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumSist[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumStatScaled[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumSistScaled[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumSistScaledB[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumSistScaledForLegend[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumStatRatio[nummoltMax+1+3][numregions];
  TH1F *fHistSpectrumSistRatio[nummoltMax+1+3][numregions];
  TH1F *fHistSpectrumStatDoubleRatio[nummoltMax+1+3][numregions];
  TH1F *fHistSpectrumSistDoubleRatio[nummoltMax+1+3][numregions];
  TH1F *fHistSpectrumStatMultRatio[nummoltMax+1+3][2][numregions];
  TH1F *fHistSpectrumSistMultRatio[nummoltMax+1+3][2][numregions];

  TSpline3*  splineK0s[nummoltMax+1+3][numregions];
  TF1 * lineat1 = new TF1 ("pol0", "pol0",0,8);
  lineat1->SetParameter(0,1);
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineWidth(0.5);
  Float_t AvgStat[nummoltMax+1+3][numregions]={0};
  Float_t AvgSist[nummoltMax+1+3][numregions]={0};

  //  gStyle->SetLegendBorderSize(1);
  TLegend *legendMult = new TLegend(0.22, 0.55, 0.77, 0.7);
  legendMult->SetHeader("V0M Multiplicity Percentile");
  legendMult->SetNColumns(3);
  legendMult->SetFillStyle(0);
  TLegendEntry *lheaderMult = (TLegendEntry*)legendMult->GetListOfPrimitives()->First();
  lheaderMult-> SetTextSize(0.04);

  TLegend *legendAllMult = new TLegend(0.22, 0.45, 0.77, 0.7);
  legendAllMult->SetHeader("V0M Multiplicity Percentile");
  legendAllMult->SetNColumns(2);
  legendAllMult->SetFillStyle(0);
  TLegendEntry *lheaderAllMult = (TLegendEntry*)legendAllMult->GetListOfPrimitives()->First();
  lheaderAllMult-> SetTextSize(0.04);

  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendBorderSize(0);

  TLegend *Legend1BisRatio=new TLegend(0.13,0.77,0.53,0.94);
  Legend1BisRatio->SetFillStyle(0);
  Legend1BisRatio->SetTextAlign(13);
  Legend1BisRatio->SetTextSize(0.04);
  Legend1BisRatio->AddEntry("", "#bf{This work}", "");
  if (ispp5TeV)   Legend1BisRatio->AddEntry("", "pp, #sqrt{#it{s}} = 5.02 TeV", "");
  else   Legend1BisRatio->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  Legend1BisRatio->AddEntry("", "#it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");

  for (Int_t i=0; i<3; i++) {
    sRegion1[i][0] = sRegion1K0s[i];
    sRegion1[i][1] = sRegion1Xi[i];
  }
  cout << "\n**************************"<<endl;
  Int_t numMoltLoop = nummoltMax+1+2;
  if (ispp5TeV) numMoltLoop = 2;

  for (Int_t type=1; type>=0; type--){

    TLegend *Legend1Bis=new TLegend(0.3,0.8,0.92,0.92);
    Legend1Bis->SetFillStyle(0);
    Legend1Bis->SetTextAlign(32);
    Legend1Bis->SetTextSize(0.034);
    //      Legend1Bis->AddEntry("", "#bf{ALICE Preliminary}", "");
    Legend1Bis->AddEntry("", "#bf{This work}", "");
    //      Legend1Bis->AddEntry("", "#color[0]{#bf{ALICE Preliminary}}", "");
    if (ispp5TeV)    Legend1Bis->AddEntry("", "pp, #sqrt{#it{s}} = 5.02 TeV", "");
    else     Legend1Bis->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
    Legend1Bis->AddEntry("", NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
    //      Legend1Bis->AddEntry("", sRegionBlack[iregion], "");

    for (Int_t iregion=0; iregion<numregions; iregion++){
      //      for (Int_t m=0; m<nummolt+1; m++){
      for (Int_t m=numMoltLoop; m>=0; m--){
	if (ispp5TeV) {
	  Smolt[m] = Smolt5TeV[m];
	  SmoltLegend[m] = SmoltLegend5TeV[m];
	  SmoltBis[m] = SmoltBis5TeV[m];
	  ColorMult[m] = ColorMult5TeV[m];
	  MarkerMult[m] = MarkerMult5TeV[m];
	  size[m] = size5TeV[m];
	}
	cout << "Multiplicity : " << Smolt[m] << endl;
	if (ispp5TeV){
	  fHistSpectrumSist[m][type][iregion]=(TH1F*)fileInYieldMB[type][iregion]->Get("fHistSpectrumSistAll_"+Smolt[m]);
	  fHistSpectrumStat[m][type][iregion]=(TH1F*)fileInYieldMB[type][iregion]->Get("fHistSpectrum_"+Smolt[m]);
	}
	else {
	  if (m>=3){
	    fHistSpectrumSist[m][type][iregion]=(TH1F*)fileInYieldMB[type][iregion]->Get("fHistSpectrumSistAll_"+Smolt[m]);
	    fHistSpectrumStat[m][type][iregion]=(TH1F*)fileInYieldMB[type][iregion]->Get("fHistSpectrum_"+Smolt[m]);
	  }
	  else {
	    fHistSpectrumSist[m][type][iregion]=(TH1F*)fileInYieldHM[type][iregion]->Get("fHistSpectrumSistAll_"+Smolt[m]);
	    fHistSpectrumStat[m][type][iregion]=(TH1F*)fileInYieldHM[type][iregion]->Get("fHistSpectrum_"+Smolt[m]);
	  }
	}
	if (!fHistSpectrumStat[m][type][iregion]) { cout << " no hist spectrum stat" << endl;  return;}
	if (!fHistSpectrumSist[m][type][iregion]) { cout << " no hist spectrum sist" << endl;  return;}

	fHistSpectrumSist[m][type][iregion]->SetName("fHistSpectrumSist_"+Smolt[m]+Form("_type%i_region%i", type, iregion));
	fHistSpectrumStat[m][type][iregion]->SetName("fHistSpectrumStat_"+Smolt[m]+Form("_type%i_region%i", type, iregion));
	fHistSpectrumSistScaled[m][type][iregion]=(TH1F*)	fHistSpectrumSist[m][type][iregion]->Clone("fHistSpectrumSistScaled_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
	fHistSpectrumStatScaled[m][type][iregion]=(TH1F*)	fHistSpectrumStat[m][type][iregion]->Clone("fHistSpectrumStatScaled_"+Smolt[m]+Form("_type%i_region%i", type,iregion));

	//Xi/K0s ratio of spectra:
	if (type==1){
	  fHistSpectrumSistRatio[m][iregion]=(TH1F*)	fHistSpectrumSist[m][type][iregion]->Clone("fHistSpectrumSistRatio_"+Smolt[m]+Form("_region%i", iregion));
	  fHistSpectrumStatRatio[m][iregion]=(TH1F*)	fHistSpectrumStat[m][type][iregion]->Clone("fHistSpectrumStatRatio_"+Smolt[m]+Form("_region%i", iregion));
	}
	else {
	  splineK0s[m][iregion] = new TSpline3(fHistSpectrumStat[m][type][iregion] , Form("splineK0s_m%i_region%i",m, iregion),0.1, 3 );
	     
	  //I calculate the average stat and sist error in the pt region where I do the spline
	  Int_t Counter=0;
	  AvgStat[m][iregion] = 0;
	  AvgSist[m][iregion] = 0;
	  for (Int_t b=1; b<=  fHistSpectrumStat[m][type][iregion]->GetNbinsX(); b++){
	    if (fHistSpectrumStat[m][type][iregion]->GetBinContent(b)==0) continue;
	    if (fHistSpectrumStat[m][type][iregion]->GetXaxis()->GetBinLowEdge(b) <2.){
	      Counter++;
	      AvgStat[m][iregion] += fHistSpectrumStat[m][0][iregion]->GetBinError(b)/ fHistSpectrumStat[m][0][iregion]->GetBinContent(b);
	      AvgSist[m][iregion] += fHistSpectrumSist[m][0][iregion]->GetBinError(b)/ fHistSpectrumStat[m][0][iregion]->GetBinContent(b);
	    }
	  }
	  AvgStat[m][iregion]=  AvgStat[m][iregion]/Counter;
	  AvgSist[m][iregion] =  AvgSist[m][iregion]/Counter;
  
	  for (Int_t b=2; b<=  fHistSpectrumStatRatio[m][iregion]->GetNbinsX(); b++){
	    if (fHistSpectrumStat[m][1][iregion]->GetBinContent(b)==0){
	      fHistSpectrumStatRatio[m][iregion]-> SetBinContent(b,0);
	      fHistSpectrumSistRatio[m][iregion]-> SetBinContent(b,0);
	      fHistSpectrumStatRatio[m][iregion]-> SetBinError(b,0);
	      fHistSpectrumSistRatio[m][iregion]-> SetBinError(b,0);
	    }
	    else if ( fHistSpectrumStatRatio[m][iregion]->GetXaxis()->GetBinLowEdge(b) <2.){
	      fHistSpectrumStatRatio[m][iregion]-> SetBinContent(b,fHistSpectrumStat[m][1][iregion]->GetBinContent(b)/splineK0s[m][iregion]->Eval(fHistSpectrumStat[m][1][iregion]->GetBinCenter(b)));
	      fHistSpectrumStatRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumStat[m][1][iregion]->GetBinError(b)/fHistSpectrumStat[m][1][iregion]->GetBinContent(b),2) + pow(AvgStat[m][iregion],2)) *  fHistSpectrumStatRatio[m][iregion]->GetBinContent(b) ); 
	      fHistSpectrumSistRatio[m][iregion]-> SetBinContent(b,fHistSpectrumSist[m][1][iregion]->GetBinContent(b)/splineK0s[m][iregion]->Eval(fHistSpectrumSist[m][1][iregion]->GetBinCenter(b)));
	      fHistSpectrumSistRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumSist[m][1][iregion]->GetBinError(b)/fHistSpectrumSist[m][1][iregion]->GetBinContent(b),2) + pow(AvgSist[m][iregion],2)) *   fHistSpectrumSistRatio[m][iregion]->GetBinContent(b) ); 

	    }
	    else {
	      bcorr = fHistSpectrumStat[m][0][iregion]->FindBin(fHistSpectrumStat[m][1][iregion]->GetBinCenter(b));
	      fHistSpectrumStatRatio[m][iregion]-> SetBinContent(b,fHistSpectrumStat[m][1][iregion]->GetBinContent(b)/fHistSpectrumStat[m][0][iregion]->GetBinContent(bcorr));
	      fHistSpectrumStatRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumStat[m][1][iregion]->GetBinError(b)/fHistSpectrumStat[m][1][iregion]->GetBinContent(b),2) + pow(fHistSpectrumStat[m][0][iregion]->GetBinError(bcorr)/fHistSpectrumStat[m][0][iregion]->GetBinContent(bcorr),2))* fHistSpectrumStatRatio[m][iregion]->GetBinContent(b));
	      fHistSpectrumSistRatio[m][iregion]-> SetBinContent(b,fHistSpectrumSist[m][1][iregion]->GetBinContent(b)/fHistSpectrumSist[m][0][iregion]->GetBinContent(bcorr));
	      fHistSpectrumSistRatio[m][iregion]-> SetBinError(b,sqrt(pow(fHistSpectrumSist[m][1][iregion]->GetBinError(b)/fHistSpectrumSist[m][1][iregion]->GetBinContent(b),2) + pow(fHistSpectrumSist[m][0][iregion]->GetBinError(bcorr)/fHistSpectrumSist[m][0][iregion]->GetBinContent(bcorr),2))* fHistSpectrumSistRatio[m][iregion]->GetBinContent(b));

	    }
	  } // end loop pt bin
	 
	  /*
	  canvasPtSpectraK0s->cd(m+1);
	  fHistSpectrumStat[m][0][iregion]->Draw("same ");
	  fHistSpectrumSist[m][0][iregion]->SetFillStyle(0);
	  fHistSpectrumSist[m][0][iregion]->Draw("same e2");
	  splineK0s[m][iregion]->SetLineColor(Color[iregion]);
	  //splineK0s[m][iregion]->Draw("same");
	  */
	  if (m==0 || m ==1 || m==2 || m==6 || m==7){
	    fHistSpectrumStatRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(0.6),0);
	    fHistSpectrumStatRatio[m][iregion]->SetBinError(fHistSpectrumStatRatio[m][iregion]->FindBin(0.6),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(0.6),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinError(fHistSpectrumSistRatio[m][iregion]->FindBin(0.6),0);
	  }
	  if (m==7){
	    fHistSpectrumStatRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(5),0);
	    fHistSpectrumStatRatio[m][iregion]->SetBinError(fHistSpectrumStatRatio[m][iregion]->FindBin(5),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinContent(fHistSpectrumStatRatio[m][iregion]->FindBin(5),0);
	    fHistSpectrumSistRatio[m][iregion]->SetBinError(fHistSpectrumSistRatio[m][iregion]->FindBin(5),0);
	  }

	  StyleHistoYield(fHistSpectrumStatRatio[m][iregion], 10e-4, 0.4-10e-4,  ColorMult[m], MarkerMult[m], titleX, titleYSpectraRatio,  "", 0.8,1.15, YoffsetSpectra[type] );
	  StyleHistoYield(fHistSpectrumSistRatio[m][iregion], 10e-4, 0.4-10e-4,  ColorMult[m], MarkerMult[m], titleX, titleYSpectraRatio, "", 0.8 ,1.15, YoffsetSpectra[type] );
	  //	  fHistSpectrumStatRatio[m][iregion]->SetMarkerSize(0.8);
	  //	  fHistSpectrumSistRatio[m][iregion]->SetMarkerSize(0.8);

	  /*
	  for (Int_t b=1; b<=  fHistSpectrumStatRatio[m][iregion]->GetNbinsX(); b++){
	    bcorr = fHistSpectrumStat[m][0][iregion]->FindBin(fHistSpectrumStat[m][1][iregion]->GetBinCenter(b));	
	    cout << "\n iregion " << iregion << " m: " << m << " low edge bin: " <<	    fHistSpectrumStatRatio[m][iregion]->GetXaxis()->GetBinLowEdge(b)<< endl;
	    cout << " num info " << fHistSpectrumStat[m][1][iregion]->GetBinContent(b) << "+-" << fHistSpectrumStat[m][1][iregion]->GetBinError(b) << endl;
	    cout << " denum info " << fHistSpectrumStat[m][0][iregion]->GetBinContent(bcorr) << " or spline " << splineK0s[m][iregion]->Eval(fHistSpectrumSist[m][1][iregion]->GetBinCenter(b))<< "+-" << AvgStat[m][iregion] << endl;
	    cout << fHistSpectrumStatRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumStatRatio[m][iregion]->GetBinError(b)<< endl;
	    cout << fHistSpectrumSistRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumSistRatio[m][iregion]->GetBinError(b)<< endl;
	  }
	  */

	  fHistSpectrumStatDoubleRatio[m][iregion]=(TH1F*) fHistSpectrumStatRatio[m][iregion]->Clone("fHistSpectrumStatDoubleRatio_"+Smolt[m]+Form("_region%i", iregion));
	  fHistSpectrumSistDoubleRatio[m][iregion]=(TH1F*) fHistSpectrumSistRatio[m][iregion]->Clone("fHistSpectrumSistDoubleRatio_"+Smolt[m]+Form("_region%i", iregion));
	  
	  fHistSpectrumSistDoubleRatio[m][iregion]->Sumw2();
	  fHistSpectrumSistDoubleRatio[m][iregion]->Divide(fHistSpectrumSistRatio[ChosenMult][iregion]);
	  fHistSpectrumStatDoubleRatio[m][iregion]->Sumw2();
	  fHistSpectrumStatDoubleRatio[m][iregion]->Divide(fHistSpectrumStatRatio[ChosenMult][iregion]);

	  /*
	  if (iregion==0){
	    for (Int_t b=1; b<=  fHistSpectrumStatDoubleRatio[m][iregion]->GetNbinsX(); b++){
	      cout << "\n iregion " << iregion << " m: " << m << " low edge bin: " <<	    fHistSpectrumStatDoubleRatio[m][iregion]->GetXaxis()->GetBinLowEdge(b)<< endl;
	      cout << fHistSpectrumStatDoubleRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumStatDoubleRatio[m][iregion]->GetBinError(b)<< endl;
	      cout << fHistSpectrumSistDoubleRatio[m][iregion]->GetBinContent(b) << "+-" << fHistSpectrumSistDoubleRatio[m][iregion]->GetBinError(b)<< endl;
	    }
	  }
	  */
	}//end type 0 else
      }//end loop on mult

      for (Int_t m=0; m<=numMoltLoop; m++){
	if (type==0){
	  cout <<"\nDrawing pt spectra ratio Xi/K0s\e[39m" << endl;
	  canvasPtSpectraRatio->cd(iregion+1);
	  StyleCanvas(canvasPtSpectraRatio, 0.01, 0.14, 0.18, 0.01);
	  MultOK = (m==0 || m == numMoltLoop-1 || m == ChosenMult);
	  if (MultOK){
	    if (iregion==0) legendMult->AddEntry(fHistSpectrumStatRatio[m][iregion], SmoltLegend[m], "pl");
	    fHistSpectrumStatRatio[m][iregion]->DrawClone("same e");
	    fHistSpectrumSistRatio[m][iregion]->SetFillStyle(0);
	    fHistSpectrumSistRatio[m][iregion]->DrawClone("same e2");
	    if (m==numMoltLoop && iregion==0) legendMult->Draw("");
	  }
	  TLegend *legendRegionBlackRatio=new TLegend(0.67, 0.82, 0.94, 0.94);
	  TLegendEntry* lReAll1bRatio;
	  TLegendEntry* lReAll2bRatio;
	  lReAll1bRatio = legendRegionBlackRatio->AddEntry("", sRegionBlack[iregion], "");
	  lReAll1bRatio->SetTextSize(0.044);
	  lReAll1bRatio->SetTextAlign(32);
	  lReAll2bRatio = legendRegionBlackRatio->AddEntry("", sRegion1[iregion][type], "");
	  lReAll2bRatio->SetTextSize(0.039);
	  lReAll2bRatio->SetTextAlign(32);

	  legendRegionBlackRatio->Draw("");
	  if (iregion==0) Legend1BisRatio->Draw("");

	  canvasPtSpectraRatioAllMult->cd(iregion+1);
	  StyleCanvas(canvasPtSpectraRatioAllMult, 0.01, 0.14, 0.18, 0.01);
	  //	  StyleCanvas(canvasPtSpectraRatioAllMult, 0.05, 0.1, 0.1, 0.05);
	  if (iregion==0) legendAllMult->AddEntry(fHistSpectrumStatRatio[m][iregion], SmoltLegend[m], "pl");
	  fHistSpectrumStatRatio[m][iregion]->DrawClone("same e");
	  fHistSpectrumSistRatio[m][iregion]->SetFillStyle(0);
	  fHistSpectrumSistRatio[m][iregion]->DrawClone("same e2");
	  if (m==numMoltLoop && iregion==0) legendAllMult->Draw("");
	  legendRegionBlackRatio->Draw("");
	  if (iregion==0) Legend1BisRatio->Draw("");

	  canvasPtSpectraRatio->cd(iregion+4);
	  StyleCanvas(canvasPtSpectraRatio, 0.01, 0.14, 0.18, 0.01);
	  StyleHistoYield(fHistSpectrumStatDoubleRatio[m][iregion], 10e-4, 2.2-10e-4,  ColorMult[m], MarkerMult[m], titleX, "Ratio to 0-100%", "",0.8 ,1.15, YoffsetSpectra[type] );
	  StyleHistoYield(fHistSpectrumSistDoubleRatio[m][iregion], 10e-4, 2.2-10e-4,  ColorMult[m], MarkerMult[m], titleX, "Ratio to 0-100%", "", 0.8 ,1.15, YoffsetSpectra[type] );
	  fHistSpectrumStatDoubleRatio[m][iregion]->SetMarkerSize(1.2);
	  fHistSpectrumSistDoubleRatio[m][iregion]->SetMarkerSize(1.2);
	  MultOK = (m==0 || m==numMoltLoop-1);
	  if (MultOK){
	    fHistSpectrumStatDoubleRatio[m][iregion]->Draw("same p");
	    fHistSpectrumSistDoubleRatio[m][iregion]->SetFillStyle(0);
	    fHistSpectrumSistDoubleRatio[m][iregion]->Draw("same e2");
	    lineat1->Draw("same");
	  }
	}
      }

      cout << "\n\e[35mDrawing spectra " << endl;
      //drawing spectra
      canvasPtSpectra[iregion][type]->cd();
      gPad->SetLogy();
      canvasPtSpectra[iregion][type]->SetFillColor(0);
      canvasPtSpectra[iregion][type]->SetTickx(1);
      canvasPtSpectra[iregion][type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.2);
      gPad->SetBottomMargin(0.13);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      TLegend *legendRegionBlack=new TLegend(0.65, 0.7, 0.92, 0.78);
      TLegendEntry* lReAll1b;
      TLegendEntry* lReAll2b;
      lReAll1b = legendRegionBlack->AddEntry("", sRegionBlack[iregion], "");
      lReAll1b->SetTextSize(0.036);
      lReAll1b->SetTextAlign(32);
      lReAll2b = legendRegionBlack->AddEntry("", sRegion1[iregion][type], "");
      lReAll2b->SetTextSize(0.029);
      lReAll2b->SetTextAlign(32);

      TLegend* legendFitFunction = new TLegend(0.65, 0.6, 0.92, 0.68);
      TLegendEntry* lFitFunction1;

      fHistSpectrumSistScaledB[nummoltMax][type][iregion]= (TH1F*) 	fHistSpectrumSistScaled[0][type][iregion]->Clone("sistblack");
      fHistSpectrumSistScaledB[nummoltMax][type][iregion]->SetLineColor(1);
      fHistSpectrumSistScaledB[nummoltMax][type][iregion]->SetMarkerColor(1);
      fHistSpectrumSistScaledB[nummoltMax][type][iregion]->SetMarkerStyle(20);
      TLegend *legendStatBoxK0s=new TLegend(0.25, 0.33, 0.47, 0.42);
      legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[nummoltMax][type][iregion], "stat. error", "pe");
      legendStatBoxK0s->AddEntry(fHistSpectrumSistScaledB[nummoltMax][type][iregion], "syst. error", "ef");

      TLegend *legendStatBoxXi=new TLegend(0.23, 0.83, 0.42, 0.92);
      legendStatBoxXi->AddEntry(fHistSpectrumSistScaledB[nummoltMax][type][iregion], "stat. error", "pe");
      legendStatBoxXi->AddEntry(fHistSpectrumSistScaledB[nummoltMax][type][iregion], "syst. error", "ef");

      TLegend *legendStatBoxRatio=new TLegend(0.12, 0.58, 0.32, 0.67);
      legendStatBoxRatio->AddEntry(fHistSpectrumSistScaledB[nummoltMax][type][iregion], "stat. error", "pe");
      legendStatBoxRatio->AddEntry(fHistSpectrumSistScaledB[nummoltMax][type][iregion], "syst. error", "ef");

      TLegend *LegendMolt=new TLegend(0.25,0.16,0.9,0.31);
      LegendMolt->SetNColumns(3);
      LegendMolt->SetFillStyle(0);
      LegendMolt->SetHeader("V0M Multiplicity Percentile");
      TLegendEntry *lheader = (TLegendEntry*)LegendMolt->GetListOfPrimitives()->First();
      lheader-> SetTextSize(0.033);
      //LegendMolt->SetTextAlign(32);

      TLegend *LegendMoltRatio;
      if (ispp5TeV) LegendMoltRatio =new TLegend(0.12,0.72,0.52,0.9);
      else LegendMoltRatio =new TLegend(0.12,0.72,0.62,0.9);
      LegendMoltRatio->SetNColumns(2);
      LegendMoltRatio->SetFillStyle(0);
      LegendMoltRatio->SetHeader("V0M Multiplicity Percentile");
      TLegendEntry *lheaderRatio = (TLegendEntry*)LegendMoltRatio->GetListOfPrimitives()->First();
      lheaderRatio-> SetTextSize(0.033);

      TLine * lineat1Mult = new TLine(0, 1,8,1);
      lineat1Mult->SetLineColor(1);
      lineat1Mult->SetLineStyle(2);

      for (Int_t m=0; m<=numMoltLoop; m++){
	if(type==0){
	  if (ispp5TeV){
	    ScaleFactor[m] = ScaleFactorK0s5TeV[m];
	    sScaleFactor[m] = sScaleFactorK0s5TeV[m];
	  }
	  else {
	    ScaleFactor[m] = ScaleFactorK0s[m];
	    sScaleFactor[m] = sScaleFactorK0s[m];
	  }
	}
	else {
	  if (ispp5TeV){
	    ScaleFactor[m] = ScaleFactorXi5TeV[m];
	    sScaleFactor[m] = sScaleFactorXi5TeV[m];
	  }
	  else {
	    ScaleFactor[m] = ScaleFactorXi[m];
	    sScaleFactor[m] = sScaleFactorXi[m];
	  }
	}
	fHistSpectrumStatScaled[m][type][iregion]->Scale(ScaleFactor[m]);
	fHistSpectrumSistScaled[m][type][iregion]->Scale(ScaleFactor[m]);
	StyleHistoYield(fHistSpectrumStatScaled[m][type][iregion], LimInfSpectra[type][iregion], LimSupSpectra[type][iregion], ColorMult[m], MarkerMult[m], titleX, titleY,"", size[m],1.15, YoffsetSpectra[type] );
	StyleHistoYield(fHistSpectrumSistScaled[m][type][iregion], LimInfSpectra[type][iregion], LimSupSpectra[type][iregion], ColorMult[m], MarkerMult[m], titleX, titleY,"", size[m], 1.15,  YoffsetSpectra[type]);

	fHistSpectrumStatScaled[m][type][iregion]->DrawClone("same e0x0");  
	fHistSpectrumSistScaled[m][type][iregion]->SetFillStyle(0);
	fHistSpectrumSistScaled[m][type][iregion]->DrawClone("same e2");
	fHistSpectrumStatScaled[m][type][iregion]->Scale(1./ScaleFactor[m]);
	fHistSpectrumSistScaled[m][type][iregion]->Scale(1./ScaleFactor[m]);
	fHistSpectrumSistScaledForLegend[m][type][iregion] = (TH1F*) 	fHistSpectrumSistScaled[m][type][iregion]->Clone("fHistSpectrumSistScaledForLegend_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
 	LegendMolt->AddEntry(fHistSpectrumSistScaledForLegend[m][type][iregion],SmoltBis[m]+"%"+sScaleFactor[m]+" ","pef");
	if (m!=ChosenMult) 	LegendMoltRatio->AddEntry(fHistSpectrumSistScaledForLegend[m][type][iregion],SmoltBis[m]+"% ","pef");
      }

      LegendMolt->Draw("");
      Legend1Bis->Draw("");
      legendRegionBlack->Draw("");
      if(type==0)legendStatBoxK0s->Draw("");
      else if(type==1)legendStatBoxXi->Draw("");

      canvasPtSpectraMultRatio[iregion][type]->cd();
      canvasPtSpectraMultRatio[iregion][type]->SetFillColor(0);
      canvasPtSpectraMultRatio[iregion][type]->SetTickx(1);
      canvasPtSpectraMultRatio[iregion][type]->SetTicky(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.08);
      gPad->SetBottomMargin(0.13);
      gPad->SetRightMargin(0.05);
      gStyle->SetLegendBorderSize(0);
      gStyle->SetLegendFillColor(0);
      gStyle->SetLegendFont(42);

      for (Int_t m=0; m<=numMoltLoop; m++){
	fHistSpectrumStatMultRatio[m][type][iregion]= (TH1F*) fHistSpectrumStat[m][type][iregion]->Clone("fHistSpectrumStatMultRatio_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
	fHistSpectrumSistMultRatio[m][type][iregion]= (TH1F*) fHistSpectrumSist[m][type][iregion]->Clone("fHistSpectrumSistMultRatio_"+Smolt[m]+Form("_type%i_region%i", type,iregion));
	fHistSpectrumStatMultRatio[m][type][iregion]->Divide(fHistSpectrumStat[ChosenMult][type][iregion]);
	fHistSpectrumSistMultRatio[m][type][iregion]->Divide(fHistSpectrumSist[ChosenMult][type][iregion]);
	StyleHistoYield(fHistSpectrumStatMultRatio[m][type][iregion], 10e-4, LimSupMultRatio[type], ColorMult[m], MarkerMult[m], titleX, "Ratio to 0-100%","", size[m],1.15, YoffsetSpectraRatio[type]);
	StyleHistoYield(fHistSpectrumSistMultRatio[m][type][iregion], 10e-4, LimSupMultRatio[type], ColorMult[m], MarkerMult[m], titleX,  "Ratio to 0-100%","", size[m], 1.15,  YoffsetSpectraRatio[type]);

	if (m!=ChosenMult){
	  fHistSpectrumStatMultRatio[m][type][iregion]->Draw("same e0x0");  
	  fHistSpectrumSistMultRatio[m][type][iregion]->SetFillStyle(0);
	  fHistSpectrumSistMultRatio[m][type][iregion]->Draw("same e2");
	}
      }
      LegendMoltRatio->Draw("");
      Legend1Bis->Draw("");
      legendRegionBlack->Draw("");
      legendStatBoxRatio->Draw("");
      lineat1Mult->Draw("same");
    }//end loop region
  }

  TString DirPicture = "PictureForThesis/";
  if (ispp5TeV) DirPicture += "pp5TeV/5TeV";
  else DirPicture += "pp13TeV/";

  for (Int_t type=0; type<2; type++){
    cout << "\e[32m\n\n************* " << Stipo[type] << " *********************" << endl;
    for (Int_t iregion=0; iregion<3; iregion++){
      cout <<"\n\e[32mRegion: " << SRegionType[iregion]<< "\e[39m" <<  endl;
      cout << "\e[36m\n" << Stipo[type] << " pt spectra\e[39m" << endl;	
      canvasPtSpectra[iregion][type]->SaveAs(DirPicture+"FinalPtSpectraAllMult"+tipo[type]+RegionType[iregion]+".eps");
      canvasPtSpectra[iregion][type]->SaveAs(DirPicture+"FinalPtSpectraAlMult"+tipo[type]+RegionType[iregion]+".pdf");
      if (type==0 && iregion==0)      canvasPtSpectra[iregion][type]->SaveAs(stringoutpdf+"_Plots.pdf(");
      else canvasPtSpectra[iregion][type]->SaveAs(stringoutpdf+"_Plots.pdf");
      //      canvasPtSpectraK0s->SaveAs(stringoutpdf+"_Plots.pdf");
      cout << "\e[36m\n" << Stipo[type] << " pt spectra to 0-100% multiplicity class\e[39m" << endl;	
      canvasPtSpectraMultRatio[iregion][type]->SaveAs(DirPicture+"FinalPtSpectraMultRatio"+tipo[type]+RegionType[iregion]+".pdf");
      canvasPtSpectraMultRatio[iregion][type]->SaveAs(stringoutpdf+"_Plots.pdf");
      fileout->WriteTObject(canvasPtSpectra[iregion][type]);
      fileout->WriteTObject(canvasPtSpectraMultRatio[iregion][type]);
  
    }
  }

  cout << "\n\e[32mXi/K0s pt spectra ratio\e[39m" << endl;	
  canvasPtSpectraRatio->SaveAs(DirPicture+"PtSpectraRatioK0sXi.pdf");  	
  canvasPtSpectraRatio->SaveAs(stringoutpdf+"_Plots.pdf");
  fileout->WriteTObject(canvasPtSpectraRatio);
  cout << "\n\e[32mXi/K0s pt spectra ratio in all multiplicity classes\e[39m" << endl;	
  canvasPtSpectraRatioAllMult->SaveAs(DirPicture+"PtSpectraRatioK0sXiAllMult.pdf");  	
  canvasPtSpectraRatioAllMult->SaveAs(stringoutpdf+"_Plots.pdf)");
  fileout->WriteTObject(canvasPtSpectraRatioAllMult);

  fileout->Close();

  for (Int_t type=0; type<2; type++){
    cout << "\nStarting from the file(s) for MB (for the 3 diff. regions): "  <<  PathInYieldMB[type][0] << endl;
    cout << "\nStarting from the file(s) for HM (for the 3 diff. regions): "  <<  PathInYieldHM[type][0] << endl;
  }
  cout << "\nI have created the file:\n " << stringout << "\nand the file:\n" << stringoutpdf << "_Plots.pdf\n\n" << endl;
}

