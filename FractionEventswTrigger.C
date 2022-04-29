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
#include "TFitResult.h"
#include "Macros/SetEfficiencyError.C"

void canvasStyle(TCanvas *canvas){
  canvas->cd();
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.13);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.05);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp, Float_t xOffset, Float_t yOffset, Float_t mSize){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset);
  histo->SetTitle(title);

}

TString TitleYieldRatio="#Xi/K^{0}_{S} yield ratio vs multiplicity";
TString titleRelUnc = "Relative uncertainty";
TString titleMult = "Multiplicity class ";
TString titledNdeta="#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString titleYield[2]={"K^{0}_{S}", "#Xi"};
TString TitleYYieldRatio="#it{N}_{#Xi}/#it{N}_{K^{0}_{S}}";

//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta} #Delta#it{#varphi}"};
TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})", "#it{N}_{#Xi}/#it{N}_{trigg} 1/(#Delta#it{#eta} #Delta#it{#varphi})"};
//TString titleYieldYType[2]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)", "#it{N}_{#Xi}/#it{N}_{trigg} 1/#left(#Delta#it{#eta} #Delta#it{#varphi}#right)"}; //to much space around the parenthesis
TString titlePtvsMultYType[2]={"#LT#it{p}^{K^{0}_{S}}_{T}#GT (GeV/c)", "#LT#it{p}^{#Xi}_{T}#GT (GeV/c)"};

TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};
TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Full"};
TString NameP[2]={"h#minusK_{S}^{0}", "h#minus#Xi"};
//TString sRegion[3]={"#color[628]{Near#minusside jet}","#color[418]{Out#minusof#minusjet}","#color[600]{Full}"};
TString sRegion[3]={"#color[628]{Toward leading}","#color[418]{Transverse to leading}","#color[600]{Full}"};
//TString sRegionBlack[3]={"#color[1]{Near#minusside jet}","#color[1]{Out#minusof#minusjet}","#color[1]{Full}"};
TString sRegionBlack[3]={"#color[1]{Toward leading}","#color[1]{Transverse to leading}","#color[1]{Full}"};
//TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 1.09", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 1.8", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 0.85", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 2.0", "|#Delta#it{#eta}| < 1.2, #minus#pi/2 < #Delta#it{#varphi} < 3#pi/2"};


void FractionEventswTrigger(Int_t MinRange=0){

  Int_t nummolt = 5;

  TString SInputFile[3] = {"FinalOutput/AnalysisResults1617_AOD234_hK0s.root", "FinalOutput/AnalysisResults17pq_hK0s.root", "FinalOutput/AnalysisResultsAllhK0sHM_RedNo16k.root"};
  TString DirName[3] = {"MyTask_PtTrigMin3.0_PtTrigMax15.0", "MyTask_PtTrigMin3.0_PtTrigMax15.0", "MyTask_PtTrigMin3.0_PtTrigMax15.0"};
  TString ListName[3] = {"MyOutputContainer_hK0s_Task_", "MyOutputContainer_hK0s_Task_", "MyOutputContainer_hK0s_Task_"};

  TH1F * hEventsvsMult[3];
  TH1F * hEventsvsMult_EvwTrigger[3];
  TH1F * hEventsvsMult_EvwTriggerOne[3];
  TH2F * hEventsvsMult_EvwTriggerBis[3];

  TH2F * hTriggervsMult[3];
  TH1F * hTriggervsMult1D[3][nummolt];
  TH1F * hTriggervsMult1DAllMult[3];
  TH1F * hTriggervsMult1DAllMultNorm[3];

  TF1 * Diagonal = new TF1 ( "Diagonal", "pol1", 0, 40);
  Diagonal->SetParameter(1, 0.6/40);
  Diagonal->SetParameter(0, 0);
  Diagonal->SetLineStyle(8);
  Diagonal->SetLineColor(881);

  TF1 * DiagonalNorm = new TF1 ( "DiagonalNorm", "pol1", 0, 40);
  DiagonalNorm->SetParameter(1, 1./40);
  DiagonalNorm->SetParameter(0, 0);
  DiagonalNorm->SetLineStyle(8);
  DiagonalNorm->SetLineColor(881);

  Float_t Nmolt[3][nummolt+1]={{0,5,10,30,50,100}, {0, 10, 100, 100, 100, 100}, {0, 0.01, 0.05, 0.1, 0.1, 0.1}};
  Float_t dNdEta[3][nummolt+1]={{21.2, 16.17, 11.4625, 7.135, 3.33, 6.94}, {13.89, 6.95, 0, 0, 0, 0}, {36.29, 32.57, 30.43, 0, 0, 0}};
  TH1F*  hRatioPtMin[5];
  TH1F*  hRatioPtMinGrey[5];
  TH1F*  hRatio5TeVPtMin[5];
  TH1F*  hRatioPtMinNorm[5];
  TH1F*  hRatio5TeVPtMinNorm[5];
  Float_t LimInf =0; 
  Float_t LimSup =0; 
  Float_t FracEvWTrigger =0; 
  TString System[3] = {"ppMB13TeV",  "pp5TeV", "ppHM13TeV"};
  TString NameEventsvsMult[3] = {"AllEvents_ppMB13TeV",  "AllEvents_pp5TeV", "AllEvents_ppHM13TeV"};
  TString NameEventsvsMult_EvwTrigger[3] = {"EventsWTrigger_ppMB13TeV", "EventsWTrigger_pp5TeV",  "EventsWTrigger_ppHM13TeV"};

  TFile * InputFile;


  //Set titles
  //  TString titleY = "Fraction of events with trigger particle";
  TString titleY = "#it{N}_{trigg}/#it{N}_{ev}";
  TString titleYNorm = "Fraction of events w trigger particle wrt 50-100% class";
  TString titleYBis = "#LT#trigger particles#GT / #LT#trigger particles#GT (0-100%)";
  TString titleYTer = "#LT#trigger particles#GT / event";
  TString  titleX = titledNdeta;
  Float_t xOffset =1.2;
  Float_t yOffset =0.9; //1.25
  Float_t MarkerSize[3] ={2, 2, 3};
  Int_t MarkerType[3] = {20, 21, 33};

  gStyle->SetOptStat(0);
  Int_t Color = 628; 
  Int_t Color5TeV[5] = {867, 628, 829, 878, 800};
  Int_t Color13TeV[5] = {601, 634, 418,  881, 807};
  Int_t Marker[5] = {33, 21, 20, 27, 28};
  Float_t PtMin[5] = {3, 3.5, 4, 5, 6};
  //  Float_t Size[numTypes] = {2, 1.4};

  TLegend *LegendTitle=new TLegend(0.18,0.80,0.5,0.92);
  LegendTitle->SetMargin(0);
  LegendTitle->AddEntry("", "#bf{ALICE Preliminary}", "");
  LegendTitle->SetTextSize(0.06);

  TLegend *LegendEnergy=new TLegend(0.14,0.65,0.48,0.8);
  LegendEnergy->SetTextSize(0.05);

  TLegend *legendPtMin=new TLegend(0.18, 0.48, 0.51, 0.58);
  legendPtMin->SetMargin(0.15);
  legendPtMin->SetTextSize(0.045);

  for (Int_t pt=0; pt<5; pt++){
    cout <<"\n\nPtMin: " << PtMin[pt] << endl;
    hRatioPtMin[pt] = new TH1F (Form("hRatioPtMin%i", pt), Form("hRatioPtMin%i", pt),450, 0, 45);
    hRatioPtMinGrey[pt] = new TH1F (Form("hRatioPtMinGrey%i", pt), Form("hRatioPtMinGrey%i", pt),450, 0, 45);
    hRatio5TeVPtMin[pt] = new TH1F (Form("hRatio5TeVPtMin%i", pt), Form("hRatio5TeVPtMin%i", pt),450, 0, 45);
    hRatioPtMinNorm[pt] = new TH1F (Form("hRatioPtMinNorm%i", pt), Form("hRatioPtMinNorm%i", pt),450, 0, 45);
    hRatio5TeVPtMinNorm[pt] = new TH1F (Form("hRatio5TeVPtMinNorm%i", pt), Form("hRatio5TeVPtMinNorm%i", pt),450, 0, 45);
    for (Int_t bin=1; bin <=  hRatioPtMin[pt]->GetNbinsX(); bin++){
      hRatioPtMin[pt]->SetBinContent(hRatioPtMin[pt]->FindBin(bin), 0);
      hRatioPtMin[pt]->SetBinError(hRatioPtMin[pt]->FindBin(bin), 0);
      hRatio5TeVPtMin[pt]->SetBinContent(hRatioPtMin[pt]->FindBin(bin), 0);
      hRatio5TeVPtMin[pt]->SetBinError(hRatioPtMin[pt]->FindBin(bin), 0);  
    }
  }
  TH1F *hNumberTriggervsMult = new TH1F ("hNumberTriggervsMult", "hNumberTriggervsMult" ,450, 0, 45);
  TH1F *hNumberTriggervsMult5TeV = new TH1F ("hNumberTriggervsMult5TeV", "hNumberTriggervsMult5TeV" ,450, 0, 45);
  TH1F *hNumberTriggervsMultNorm = new TH1F ("hNumberTriggervsMultNorm", "hNumberTriggervsMultNorm" ,450, 0, 45);
  TH1F *hNumberTriggervsMult5TeVNorm = new TH1F ("hNumberTriggervsMult5TeVNorm", "hNumberTriggervsMult5TeVNorm" ,450, 0, 45);

  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvasStyle(canvas);
  //  canvas->SetGrid();

  TCanvas * canvasNorm = new TCanvas("canvasNorm", "canvasNorm", 1300, 800);
  canvasStyle(canvasNorm);

  //  TLegend *legendPtMin=new TLegend(0.6, 0.2, 0.9, 0.4);

  Float_t MinFracEvwTrigger[5]={0};
  for (Int_t i =0; i<3; i++){
    cout <<"Analysing : " << System[i] << endl;
    InputFile = new TFile(SInputFile[i], "");
    if (!InputFile) return;
    TDirectoryFile * dir = (TDirectoryFile*)InputFile->Get(DirName[i]); 
    if (!dir) return;
    TList * list = (TList*) dir->Get(ListName[i]);
    if (!list) return;
    hEventsvsMult[i] = (TH1F*)list->FindObject("fHist_multiplicityAllSelEvents");
    if (!hEventsvsMult[i]) return;
    hEventsvsMult[i]->SetName(NameEventsvsMult[i]);
    hEventsvsMult_EvwTriggerOne[i] = (TH1F*)list->FindObject("fHist_multiplicity_EvwTrigger");
    if (!hEventsvsMult_EvwTriggerOne[i]) return;
    hEventsvsMult_EvwTriggerBis[i] = (TH2F*)list->FindObject("fHistPtMaxvsMultBefAll");
    if (!hEventsvsMult_EvwTriggerBis[i]) return;
    hTriggervsMult[i] = (TH2F*)list->FindObject("fHistMultvsTriggerBefAll");
    if (!hTriggervsMult[i]) return;
    //    hEventsvsMult_EvwTrigger[i] = (TH1F*)     hEventsvsMult_EvwTriggerOne[i]->Clone(NameEventsvsMult_EvwTrigger[i]);

    //number of particles with pt > 3 GeV in events with at least one such particle
    hTriggervsMult1DAllMult[i] =  (TH1F*)     hTriggervsMult[i]->ProjectionX(Form("hTriggervsMult1DAllMult%i", i), 0, 100);

    for (Int_t m = 0; m< nummolt; m++){
      if (i==1 && m>1) continue;
      if (i==2 && m>2) continue;
      LimInf = hTriggervsMult[i]->GetYaxis()->FindBin(Nmolt[i][m]);
      LimSup = hTriggervsMult[i]->GetYaxis()->FindBin(Nmolt[i][m]);
      hTriggervsMult1D[i][m]= (TH1F*)     hTriggervsMult[i]->ProjectionX(Form("hTriggervsMult1D_%i_%m", i, m), LimInf, LimSup);
      hTriggervsMult1D[i][m]->GetXaxis()->SetRangeUser(MinRange, 30);
      if (i==0 || i ==2) hNumberTriggervsMult ->SetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m]), hTriggervsMult1D[i][m]->GetMean()/hTriggervsMult1DAllMult[i]->GetMean());
      else if (i==1)  hNumberTriggervsMult5TeV ->SetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m]), hTriggervsMult1D[i][m]->GetMean()/hTriggervsMult1DAllMult[i]->GetMean());
    }


    for (Int_t pt=0; pt<5; pt++){
      if (pt>0) continue;
      cout <<"\n\nPtMin: " << PtMin[pt] << endl;
      hEventsvsMult_EvwTrigger[i] = (TH1F*)     hEventsvsMult_EvwTriggerBis[i]->ProjectionY(NameEventsvsMult_EvwTrigger[i], hEventsvsMult_EvwTriggerBis[i]->GetXaxis()->FindBin(PtMin[pt]+0.0001), hEventsvsMult_EvwTriggerBis[i]->GetXaxis()->FindBin(30));

      hEventsvsMult_EvwTrigger[i]->SetName(NameEventsvsMult_EvwTrigger[i]);

      if (i==0)  MinFracEvwTrigger[pt] = hEventsvsMult_EvwTrigger[0]->Integral( hEventsvsMult_EvwTrigger[0]->FindBin(Nmolt[0][4]),  hEventsvsMult_EvwTrigger[0]->FindBin(Nmolt[0][5]))/hEventsvsMult[0]->Integral(hEventsvsMult_EvwTrigger[0]->FindBin(Nmolt[0][4]), hEventsvsMult_EvwTrigger[0]->FindBin(Nmolt[0][5]));

      for (Int_t m = 0; m< nummolt; m++){
	if (i==1 && m>1) continue;
	if (i==2 && m>2) continue;
	LimInf = hEventsvsMult_EvwTrigger[i]->FindBin(Nmolt[i][m]);
	LimSup = hEventsvsMult_EvwTrigger[i]->FindBin(Nmolt[i][m+1]);
	FracEvWTrigger = hEventsvsMult_EvwTrigger[i]->Integral(LimInf, LimSup)/hEventsvsMult[i]->Integral(LimInf, LimSup);

	if (i==0 || i ==2){
	  hRatioPtMin[pt]->SetBinContent(hRatioPtMin[pt]->FindBin(dNdEta[i][m]), FracEvWTrigger);
	  hRatioPtMin[pt]->SetBinError(hRatioPtMin[pt]->FindBin(dNdEta[i][m]), SetEfficiencyError(hEventsvsMult_EvwTrigger[i]->Integral(LimInf, LimSup), hEventsvsMult[i]->Integral(LimInf, LimSup)));
	  hRatioPtMinNorm[pt]->SetBinContent(hRatioPtMin[pt]->FindBin(dNdEta[i][m]), FracEvWTrigger/MinFracEvwTrigger[pt]);
	  hRatioPtMinNorm[pt]->SetBinError(hRatio5TeVPtMin[pt]->FindBin(dNdEta[i][m]), 0);
	  cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ")  Fraction of events w trigger particle wrt lowest mult class: " << hRatioPtMinNorm[pt]->GetBinContent(hRatioPtMin[pt]->FindBin(dNdEta[i][m])) << endl;
	}
	else {
	  hRatio5TeVPtMin[pt]->SetBinContent(hRatio5TeVPtMin[pt]->FindBin(dNdEta[i][m]), FracEvWTrigger);
	  hRatio5TeVPtMin[pt]->SetBinError(hRatio5TeVPtMin[pt]->FindBin(dNdEta[i][m]), SetEfficiencyError(hEventsvsMult_EvwTrigger[i]->Integral(LimInf, LimSup), hEventsvsMult[i]->Integral(LimInf, LimSup)));
	  hRatio5TeVPtMinNorm[pt]->SetBinContent(hRatio5TeVPtMin[pt]->FindBin(dNdEta[i][m]), FracEvWTrigger/MinFracEvwTrigger[pt]);
	  hRatio5TeVPtMinNorm[pt]->SetBinError(hRatio5TeVPtMin[pt]->FindBin(dNdEta[i][m]), 0);
	  cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ")  Fraction of events w trigger particle wrt lowest mult class: " << hRatio5TeVPtMinNorm[pt]->GetBinContent(hRatioPtMin[pt]->FindBin(dNdEta[i][m])) << endl;
	}
	cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ")  Fraction of events w trigger particle: " << FracEvWTrigger << endl;
	cout << "Min fraction of events w trigger particle: " << MinFracEvwTrigger[pt] << endl;
	// cout << hRatioPtMin[pt]->GetBinContent(hRatioPtMin[pt]->FindBin(dNdEta[i][m])) << " +- " << hRatioPtMin[pt]->GetBinError(hRatioPtMin[pt]->FindBin(dNdEta[i][m])) << endl;
      }
      StyleHisto(hRatioPtMin[pt], 10e-5, 0.6-10e-5, Color13TeV[pt], Marker[pt], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, 3);
      StyleHisto(hRatio5TeVPtMin[pt], 10e-5, 0.6-10e-5, Color5TeV[pt], Marker[pt], titleX, titleY, "" , 1,0, 40, xOffset, yOffset, 3);
      StyleHisto(hRatioPtMinNorm[pt], 0.9, 75, Color13TeV[pt], Marker[pt], titleX, titleYNorm, "" , 1,0, 40, xOffset, yOffset, 3);
      StyleHisto(hRatio5TeVPtMinNorm[pt], 0.9, 75, Color5TeV[pt], Marker[pt], titleX, titleYNorm, "" , 1,0, 40, xOffset, yOffset, 3);
      hRatioPtMinNorm[pt]->GetYaxis()->SetTitleSize(0.03);
      hRatio5TeVPtMinNorm[pt]->GetYaxis()->SetTitleSize(0.03);

      hRatioPtMin[pt]->GetYaxis()->SetTitleSize(0.06);
      hRatio5TeVPtMin[pt]->GetYaxis()->SetTitleSize(0.06);

      if (i==0)  {
	if (pt==0){
	  hRatioPtMinGrey[pt] = (TH1F*) hRatioPtMin[pt]->Clone("hRatioPtMinGrey");
	  hRatioPtMinGrey[pt]->SetMarkerColor(922);
	  hRatioPtMinGrey[pt]->SetLineColor(922);
	  legendPtMin->AddEntry(hRatioPtMinGrey[pt], Form("#it{p}_{T}^{trigg} > %.0f GeV/#it{c}", PtMin[pt]), "pe");
	  LegendEnergy->AddEntry(hRatioPtMin[pt], "pp, #sqrt{#it{s}} = 13 TeV", "pe");
	}
	else  legendPtMin->AddEntry(hRatioPtMin[pt], Form("#it{p}_{T}^{trigg} > %.0f GeV/#it{c}", PtMin[pt]), "pe");
      }
      else if (i==1){
	if (pt==0) LegendEnergy->AddEntry(hRatio5TeVPtMin[pt], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
      }

      canvas->cd();
      hRatioPtMin[pt]->DrawClone("same p");
      hRatio5TeVPtMin[pt]->DrawClone("same p");
      //      Diagonal->Draw("same");

      canvasNorm->cd();
      hRatioPtMinNorm[pt]->DrawClone("same p");
      hRatio5TeVPtMinNorm[pt]->DrawClone("same p");
    } 
  }

  canvas->cd();
  LegendTitle->Draw("");
  legendPtMin->Draw("");
  LegendEnergy->Draw("");
  canvasNorm->cd();
  legendPtMin->Draw("");

  TLegend *legendEnergyBoxColor=new TLegend(0.16, 0.62, 0.39, 0.74);
  legendEnergyBoxColor->AddEntry(hRatio5TeVPtMin[0], "pp, #sqrt{#it{s}} = 5.02 TeV", "pe");
  legendEnergyBoxColor->AddEntry(hRatioPtMin[0], "pp, #sqrt{#it{s}} = 13 TeV", "pe");
  //  legendEnergyBoxColor->Draw("");

  cout <<"\n\n\n"<< endl;
  TF1 * pol1Fit = new TF1("pol1", "pol1", 0, 25);
  //number of particles with pt > 3 GeV in events with at least one such particle
  for (Int_t i =0; i<3; i++){
    cout <<"\nAnalysing : " << System[i] << endl;
    InputFile = new TFile(SInputFile[i], "");
    if (!InputFile) return;
    TDirectoryFile * dir = (TDirectoryFile*)InputFile->Get(DirName[i]); 
    if (!dir) return;
    TList * list = (TList*) dir->Get(ListName[i]);
    if (!list) return;
    hTriggervsMult[i] = (TH2F*)list->FindObject("fHistMultvsTriggerBefAll");
    if (!hTriggervsMult[i]) return;
 
    hTriggervsMult1DAllMult[i] =  (TH1F*)     hTriggervsMult[i]->ProjectionX(Form("hTriggervsMult1DAllMult%i", i), 0, 100);
    hTriggervsMult1DAllMult[i]->GetXaxis()->SetRangeUser(MinRange, 30);

    for (Int_t m = 0; m< nummolt; m++){
      if (i==1 && m>1) continue;
      if (i==2 && m>2) continue;
      LimInf = hTriggervsMult[i]->GetYaxis()->FindBin(Nmolt[i][m]);
      LimSup = hTriggervsMult[i]->GetYaxis()->FindBin(Nmolt[i][m]);
      hTriggervsMult1D[i][m]= (TH1F*)     hTriggervsMult[i]->ProjectionX(Form("hTriggervsMult1D_%i_%m", i, m), LimInf, LimSup);
      hTriggervsMult1D[i][m]->GetXaxis()->SetRangeUser(MinRange, 30);
      if (i==0 || i ==2) {
	hNumberTriggervsMultNorm ->SetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m]), hTriggervsMult1D[i][m]->GetMean()/hTriggervsMult1DAllMult[0]->GetMean());
	hNumberTriggervsMult ->SetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m]), hTriggervsMult1D[i][m]->GetMean());
	hNumberTriggervsMultNorm ->SetBinError(hNumberTriggervsMult->FindBin(dNdEta[i][m]), 0);
	hNumberTriggervsMult ->SetBinError(hNumberTriggervsMult->FindBin(dNdEta[i][m]),hTriggervsMult1D[i][m]->GetMeanError());
	cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ") Average number of charged tracks with pt> 3 GeV/c (norm to MB) " << hNumberTriggervsMultNorm->GetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m])) << endl;
	cout << "error  " << hNumberTriggervsMult->GetBinError(hNumberTriggervsMult->FindBin(dNdEta[i][m]))<< endl;
      }
      else if (i==1){
	hNumberTriggervsMult5TeVNorm ->SetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m]), hTriggervsMult1D[i][m]->GetMean()/hTriggervsMult1DAllMult[1]->GetMean());
	hNumberTriggervsMult5TeVNorm ->SetBinError(hNumberTriggervsMult->FindBin(dNdEta[i][m]),0);
	hNumberTriggervsMult5TeV ->SetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m]), hTriggervsMult1D[i][m]->GetMean());
	hNumberTriggervsMult5TeV ->SetBinError(hNumberTriggervsMult->FindBin(dNdEta[i][m]),hTriggervsMult1D[i][m]->GetMeanError());
	cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ") Average number of charged tracks with pt> 3 GeV/c (norm to MB)" << hNumberTriggervsMult5TeVNorm->GetBinContent(hNumberTriggervsMult->FindBin(dNdEta[i][m])) << endl;
      }
      cout << "m " << Nmolt[i][m] << " (dNdEta = " << dNdEta[i][m]<< ") Average number of charged tracks with pt> 3 GeV/c " << hTriggervsMult1D[i][m]->GetMean() << endl;
    }
  }

  TCanvas * canvasBis = new TCanvas("canvasBis", "canvasBis", 1300, 800);
  canvasStyle(canvasBis);
  canvasBis->SetGrid();

  Float_t LowLimitNorm = 0.7;
  Float_t UpLimitNorm = 1.5;
  if (MinRange==0){
    LowLimitNorm =0;
    UpLimitNorm =20;
  }
  StyleHisto(hNumberTriggervsMultNorm, LowLimitNorm, UpLimitNorm, Color13TeV[0], Marker[0], titleX, titleYBis, "" , 1,0, 40, xOffset, yOffset, 2);
  StyleHisto(hNumberTriggervsMult5TeVNorm, LowLimitNorm, UpLimitNorm, Color5TeV[0], Marker[0], titleX, titleYBis, "" , 1,0, 40, xOffset, yOffset, 2);
  hNumberTriggervsMultNorm->GetYaxis()->SetTitleSize(0.035);
  hNumberTriggervsMult5TeVNorm->GetYaxis()->SetTitleSize(0.035);
  hNumberTriggervsMultNorm->DrawClone("same p");
  hNumberTriggervsMult5TeVNorm->DrawClone("same p");
  legendEnergyBoxColor->Draw("");

  TCanvas * canvasTer = new TCanvas("canvasTer", "canvasTer", 1300, 800);
  canvasStyle(canvasTer);
  canvasTer->SetGrid();

  Float_t LowLimit = 1;
  Float_t UpLimit = 2;
  if (MinRange==0){
    LowLimit =0;
    UpLimit =1.2;
  }

  StyleHisto(hNumberTriggervsMult, LowLimit, UpLimit, Color13TeV[0], Marker[0], titleX, titleYTer, "" , 1,0, 40, xOffset, yOffset, 2);
  StyleHisto(hNumberTriggervsMult5TeV, LowLimit, UpLimit, Color5TeV[0], Marker[0], titleX, titleYTer, "" , 1,0, 40, xOffset, yOffset, 2);


  hNumberTriggervsMult->DrawClone("same p");
  hNumberTriggervsMult5TeV->DrawClone("same p");

  TLegend* legendParameters = new TLegend(0.5, 0.7, 0.8, 0.85);

  pol1Fit->SetLineWidth(0.3);
  pol1Fit->SetRange(0,25);
  hNumberTriggervsMult->Fit(pol1Fit, "R0");
  pol1Fit->SetRange(0,45);
  pol1Fit->SetLineStyle(9);
  pol1Fit->SetLineColor(Color13TeV[0]);
  pol1Fit->DrawClone("same");
  TF1 * pol1FitClone1 = (TF1*) pol1Fit->Clone("pol1FitClone1");
  legendParameters->AddEntry(pol1FitClone1, Form("y = %.3f x + %.3f (fit for dNdeta < 25)", pol1Fit->GetParameter(1), pol1Fit->GetParameter(0)) , "l");

  pol1Fit->SetRange(0,45);
  hNumberTriggervsMult->Fit(pol1Fit, "R0");
  pol1Fit->SetLineStyle(1);
  pol1Fit->SetLineColor(Color13TeV[0]+1);
  pol1Fit->DrawClone("same");
  TF1 * pol1FitClone2 = (TF1*) pol1Fit->Clone("pol1FitClone2");
  legendParameters->AddEntry(pol1FitClone2, Form("y = %.3f x + %.3f (fit for dNdeta < 45)", pol1Fit->GetParameter(1), pol1Fit->GetParameter(0)), "l");

  pol1Fit->SetRange(25, 45);
  hNumberTriggervsMult->Fit(pol1Fit, "R0");
  pol1Fit->SetRange(0,45);
  pol1Fit->SetLineStyle(9);
  pol1Fit->SetLineColor(Color5TeV[0]);
  pol1Fit->DrawClone("same");
  TF1 * pol1FitClone3 = (TF1*) pol1Fit->Clone("pol1FitClone3");
  legendParameters->AddEntry(pol1FitClone3, Form("y = %.3f x + %.3f (fit for dNdeta > 25)", pol1Fit->GetParameter(1), pol1Fit->GetParameter(0)), "l");


  //  legendParameters->AddEntry("", "Fit for "+ titledNdeta +" < 25", "");
  legendParameters->Draw("");

  //  DiagonalNorm->Draw("same");
  legendEnergyBoxColor->Draw("");
  
  cout << "Saving canvas" << endl;
  canvas->SaveAs("FractionOfEventsWithTrigger.pdf");
  canvasNorm->SaveAs("FractionOfEventsWithTriggerNorm.pdf");
  canvasBis->SaveAs("NumberofTriggerParticlesPerEventNormalised.pdf");
  canvasTer->SaveAs("NumberofTriggerParticlesPerEvent.pdf");

  TString nameCanvas = Form("FractionEventsWithTrigger_%i.pdf", MinRange);
  canvas->SaveAs(nameCanvas + "(");
  //  canvasNorm->SaveAs("FractionEventsWithTrigger.pdf");
  canvasTer->SaveAs(nameCanvas);
  canvasBis->SaveAs(nameCanvas + ")");
  
  TString fileName = "FractionOfEventsWithTrigger.root";
  cout << "Writing canvas on file " << fileName <<  endl;
  TFile* fileout = new TFile(fileName, "RECREATE");
  canvas->Write();
  canvasNorm->Write();
  canvasBis->Write();
  canvasTer->Write();
  for (Int_t pt=0; pt<5; pt++){
    hRatioPtMin[pt]->Write();
    hRatio5TeVPtMin[pt]->Write();
    hRatioPtMinNorm[pt]->Write();
    hRatio5TeVPtMinNorm[pt]->Write();
  }
  fileout->Close();
  
}
