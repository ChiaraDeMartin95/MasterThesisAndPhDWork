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
#include "Macros/constants.h"

TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
Double_t Nmolt[nummolt+1]={0,5,10,30,50,100};

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

void PlotProjectionThesis(Int_t type=1, Int_t multChosen = 5, Int_t ptchosen = 2, Int_t isProj=0, Int_t DataSample=0, Bool_t isMERatio=1){ 

  //type = 0 for K0s, 1 for Xi

  //isProj == 0  //final dPhi projections (i == 2)
  //isProj == 1  //dPhi projections without SB subtraction (i == 0)
  //isProj == 2  //dPhi projections of SB (i == 1)

  //DataSample = 0 --> MB pp 13 TeV
  //DataSample = 1 --> HM pp 13 TeV
  //DataSample = 2 --> MB pp 5.02 TeV

  const Int_t numPtV0 = 9;
  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  if (type>0){
    SPtV0[0]="0-0.5";
    SPtV0[1]="0.5-1";
  }
  else {SPtV0[0]="0-0.5"; SPtV0[1]="0.5-1";}
  Double_t NPtV0[numPtV0+1]={0.1,0.5,1,1.5,2,2.5,3,4,8, 100};
  TString SNPtV0[numPtV0+1]={"0.1","0.5","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};

  TString SPtV01[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString SNPtV01[numPtV0]={"0.1", "0.5", "0.8", "1.2","1.6", "2.0","2.5", "3.0", "4.0"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  if (type==0){
    for (Int_t b=0; b<numPtV0+1; b++){
      NPtV0[b] = NPtV01[b];
      if (b<numPtV0){
	SNPtV0[b] = SNPtV01[b];
	SPtV0[b] = SPtV01[b];
	cout << SPtV0[b] << endl;
      }
    }
  }
  Int_t color[3] = {628,418,601};

  cout << "Chosen pt interval " << SPtV0[ptchosen] << endl;
  TString stringInput = "";
  TString year="";
  TFile* inputFile;
  TString PType[2] = {"_K0s", "_Xi"}; 
  stringInput = "FinalOutput/DATA2016/histo/AngularCorrelation";
  if (type==0){
    if (DataSample == 0) year = "1617_AOD234_hK0s";
    else if (DataSample == 1) year = "AllhK0sHM_RedNo16k";
    else if (DataSample == 2) year = "17pq_hK0s";
  }
  else {
    if (DataSample == 0) year = "161718Full_AOD234_hXi";
    else if (DataSample == 1) year = "161718_HM_hXi_WithFlat16k_No18p";
    else if (DataSample == 2) year = "17pq_hXi";
  }
  stringInput += year;
  if (type==0) stringInput += "_PtBinning1";
  stringInput += PType[type];
  stringInput += "_Eta0.8_AllAssoc";
  if (type==0 && DataSample==1) stringInput += "_isBkgParab";
  stringInput += "_SysT0_SysV00_Sys0_PtMin3.0_Output";
  if (type==0)  {
    if (DataSample==0) stringInput += "_Sidebands_IsEtaEff_NewdEtaChoice_isWingsCorrectionAppliedNew_MatBudgetCorr";
    else if (DataSample==1) stringInput += "_Sidebands_IsEtaEff_MultBinning1_NewdEtaChoice_EffCorr_isWingsCorrectionAppliedNew";
    else stringInput += "";
  } 
  else {
    //    if (DataSample==0) stringInput += "_IsEtaEff_NewdEtaChoice_isWingsCorrectionAppliedNew";
    if (DataSample==0) stringInput += "_IsEtaEff";
    else stringInput += "_Sidebands_IsEtaEff_MultBinning1";
  }
  //  if (type==0 && DataSample!=2)  stringInput += "_PlotForThesis";
  stringInput += "_PlotForThesis";
  stringInput += ".root";
  inputFile = new TFile(stringInput, "");
  cout << stringInput << endl;
  if (!inputFile) {cout << "No input file " << endl; return;}

  TString SinputFileJetXi = "";
  if (DataSample==0) {
    SinputFileJetXi = "OOJComparison161718Full_AOD234_hXi_161718Full_AOD234_hXi_Xi_Eta0.8_AllAssoc_sys0_PtTrigMin3.0_PtTrigMin0.2_Output_IsEtaEff_isOOJFromAllMult_PeakyErrorProp_PlotForThesis.root";
  }
  else if (DataSample==1) {
    SinputFileJetXi = "OOJComparison161718_HM_hXi_WithFlat16k_No18p_161718_HM_hXi_WithFlat16k_No18p_Xi_Eta0.8_AllAssoc_sys0_PtTrigMin3.0_PtTrigMin3.0_Output_IsEtaEff_MultBinning1_isOOJFromAllMult_Sidebands_PeakyErrorProp_PlotForThesis.root";
  }

  TFile *inputFileJetXi = new TFile(SinputFileJetXi, "");
  if (!inputFileJetXi) {cout << "Input file for Xi OOJ not present. File name: " << SinputFileJetXi << endl; return;} 

  TString SinputFileSistFull = "";
  TString SinputFileSistBulk = "";
  TString SinputFileSistJet = "";
  TString  SInitial="";
  TString  SFinal="";
  if (DataSample==0) {
    if (type==0){
      SInitial = "FinalOutput/DATA2016/PtSpectraNew_1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_AllAssoc_SysPhi0_PtMin3.0";
      SFinal = "_EffCorr_isWingsCorrectionAppliedNew.root";
    }
    else if (type==1){
      SInitial = "FinalOutput/DATA2016/PtSpectraNew_161718Full_AOD234_hXi_Xi_Eta0.8_AllAssoc_SysPhi0_PtMin3.0";
      SFinal = ".root";
    }
  }
  else if (DataSample==1) {
    if (type==0){
      SInitial = "";
      SFinal = "";
    }
    else if (type==1){
      SInitial = "FinalOutput/DATA2016/PtSpectraNew_161718_HM_hXi_WithFlat16k_No18p_Xi_Eta0.8_AllAssoc_SysPhi0_PtMin3.0";
      SFinal = "_MultBinning1.root";
    }
  }

  SinputFileSistFull = SInitial + "_Inclusive" + SFinal;
  SinputFileSistBulk = SInitial + "_Bulk" + SFinal;
  SinputFileSistJet = SInitial + "_Jet" + SFinal;
  TFile *inputFileSystFull = new TFile(SinputFileSistFull, "");
  TFile *inputFileSystBulk = new TFile(SinputFileSistBulk, "");
  TFile *inputFileSystJet = new TFile(SinputFileSistJet, "");
  if (!inputFileSystFull) {cout << "Input file for syst errors FULL. File name: " << SinputFileSistFull << endl; return;}
  if (!inputFileSystBulk) {cout << "Input file for syst errors FULL. File name: " << SinputFileSistBulk << endl; return;}
  if (!inputFileSystJet) {cout << "Input file for syst errors FULL. File name: " << SinputFileSistJet << endl; return;}
  if (DataSample==1) Smolt[multChosen] = "0-0.1";

  TString nameSE = "SE_m"+Smolt[multChosen]+"_v"+SPtV0[ptchosen]+"_Effw";
  TString nameME = "ME_m"+Smolt[multChosen]+"_v"+SPtV0[ptchosen]+"_norm";
  TString nameAC = "ME_m"+Smolt[multChosen]+"_v"+SPtV0[ptchosen]+"_AC";

  TH2F* hSE = (TH2F*)inputFile->Get(nameSE);
  TH2F* hME = (TH2F*)inputFile->Get(nameME);
  TH2F* hAC = (TH2F*)inputFile->Get(nameAC);
  if (!hSE) {cout << "No SE histo " << nameSE << endl; return;}
  if (!hME) {cout << "No ME histo " << nameME << endl; return;}
  if (!hAC) {cout << "No AC histo " << nameAC << endl; return;}

  //  StyleHisto(hSE, Low,  Up,  color, style,  titleX,  titleY,  title, 0,  XLow,  XUp,  xOffset,  yOffset,  mSize);

  TString  SLegendMolt[6] = {"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  if (DataSample==1) SLegendMolt[multChosen] = "0-0.1%";
  TCanvas* canvasPairAcceptance = new TCanvas("canvasPairAcceptance", "canvasPairAcceptance", 1000, 700);

  TLegend *LegendSist = new TLegend(0.65, 0.6, 0.95, 0.75);
  LegendSist->SetFillStyle(0);
  LegendSist->SetTextAlign(12);
  LegendSist->SetTextSize(0.04);
  LegendSist->SetBorderSize(0);

  TLegend *LegendALICE=new TLegend(0.06,0.77,0.28,0.97);
  LegendALICE->SetFillStyle(0);
  LegendALICE->SetTextAlign(12);
  LegendALICE->SetTextSize(0.06);
  LegendALICE->SetBorderSize(0);
  LegendALICE->AddEntry("", "#bf{This work}", "");
  //  if (ispp5TeV)   LegendALICE->AddEntry("", "pp, #sqrt{#it{s}} = 5.02 TeV", "");
  //  else LegendALICE->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  LegendALICE->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  LegendALICE->AddEntry("", "V0M Multiplicity Percentile "+SLegendMolt[multChosen], "");
  if (type==0)  LegendALICE->AddEntry("", "h-K_{S}^{0} correlations", "");
  else  if (type==1)  LegendALICE->AddEntry("", "h-#Xi^{#pm} correlations", "");

  TLegend *LegendALICEProj=new TLegend(0.27,0.14,0.47,0.34);
  LegendALICEProj->SetFillStyle(0);
  LegendALICEProj->SetTextAlign(12); //left align
  if (DataSample!=1)  LegendALICEProj->SetTextSize(0.042); 
  else   LegendALICEProj->SetTextSize(0.045); 
  LegendALICEProj->SetBorderSize(0);
  LegendALICEProj->AddEntry("", "#bf{This work}", "");
  //  if (ispp5TeV)   LegendALICEProj->AddEntry("", "pp, #sqrt{#it{s}} = 5.02 TeV", "");
  //  else LegendALICEProj->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  LegendALICEProj->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  LegendALICEProj->AddEntry("", "V0M Multiplicity Percentile "+SLegendMolt[multChosen], "");
  if (type==0) {
    if (isProj==2)  LegendALICEProj->AddEntry("", "h-fake K_{S}^{0} correlations", "");
    else  LegendALICEProj->AddEntry("", "h-K_{S}^{0} correlations", "");
  }
  else  if (type==1)  LegendALICEProj->AddEntry("", "h-#Xi^{#pm} correlations", "");
  //  if (isProj==2)  LegendALICEProj->AddEntry("", "|| > 4#sigma", "");

  TLegend *LegendALICEMERatio=new TLegend(0.27,0.14,0.47,0.34);
  LegendALICEMERatio->SetFillStyle(0);
  LegendALICEMERatio->SetTextAlign(12); //left align
  if (DataSample!=1)  LegendALICEMERatio->SetTextSize(0.042); 
  else   LegendALICEMERatio->SetTextSize(0.045); 
  LegendALICEMERatio->SetBorderSize(0);
  LegendALICEMERatio->AddEntry("", "#bf{This work}", "");
  LegendALICEMERatio->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  if (type==0)  LegendALICEMERatio->AddEntry("", "h-K_{S}^{0} correlations", "");
  else  if (type==1)  LegendALICEMERatio->AddEntry("", "h-#Xi^{#pm} correlations", "");

  TLegend *LegendALICEJet=new TLegend(0.52,0.21,0.72,0.41);
  LegendALICEJet->SetFillStyle(0);
  LegendALICEJet->SetTextAlign(12); //left align
  LegendALICEJet->SetTextSize(0.045);
  LegendALICEJet->SetBorderSize(0);
  LegendALICEJet->AddEntry("", "#bf{This work}", "");
  LegendALICEJet->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  //  LegendALICEJet->AddEntry("", "V0M Multiplicity Percentile "+SLegendMolt[multChosen], "");
  LegendALICEJet->AddEntry("", "V0M "+SLegendMolt[multChosen], "");
  if (type==0)  LegendALICEJet->AddEntry("", "h-K_{S}^{0} correlations", "");
  else  if (type==1)  LegendALICEJet->AddEntry("", "h-#Xi^{#pm} correlations", "");

  TLegend *LegendALICEJetXi;
  if (DataSample==0) LegendALICEJetXi=new TLegend(0.55,0.31,0.75,0.51);
  else LegendALICEJetXi=new TLegend(0.25,0.31,0.45,0.51);
  LegendALICEJetXi->SetFillStyle(0);
  LegendALICEJetXi->SetTextAlign(12); //left align
  LegendALICEJetXi->SetTextSize(0.04);
  LegendALICEJetXi->SetBorderSize(0);
  LegendALICEJetXi->AddEntry("", "#bf{This work}", "");
  LegendALICEJetXi->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  if (DataSample==1)  LegendALICEJetXi->AddEntry("", "V0M Multiplicity Percentile 0.01-0.05%", "");
  else if (DataSample==0)   LegendALICEJetXi->AddEntry("", "V0M 10-30%", "");
  if (type==0)  LegendALICEJetXi->AddEntry("", "h-K_{S}^{0} correlations", "");
  else  if (type==1)  LegendALICEJetXi->AddEntry("", "h-#Xi^{#pm} correlations", "");

  TLegend *LegendALICEwSB=new TLegend(0.27,0.14,0.47,0.34);
  LegendALICEwSB->SetFillStyle(0);
  LegendALICEwSB->SetTextAlign(12); //left align
  if (DataSample!=1)  LegendALICEwSB->SetTextSize(0.045);
  else   LegendALICEwSB->SetTextSize(0.042);
  LegendALICEwSB->SetBorderSize(0);
  LegendALICEwSB->AddEntry("", "#bf{This work}", "");
  LegendALICEwSB->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
  LegendALICEwSB->AddEntry("", "V0M Multiplicity Percentile "+SLegendMolt[multChosen], "");
  //  LegendALICEwSB->AddEntry("", "#color[628]{|#Delta#eta| < 0.86}", "");
  //  LegendALICEwSB->AddEntry("", "#color[1]{|#Delta#eta| < 0.86}", "");
  LegendALICEwSB->AddEntry("", "#color[1]{|#Delta#eta| < 1.2}", "");

  TLegend *LegendRegions=new TLegend(0.33,0.13,0.53,0.3); //y2 = 0.33
  LegendRegions->SetFillStyle(0);
  LegendRegions->SetTextAlign(12); //left align
  LegendRegions->SetTextSize(0.06);
  LegendRegions->SetBorderSize(0);

  TLegend *LegendMolt=new TLegend(0.3,0.13,0.7,0.3); //y2 = 0.33
  LegendMolt->SetFillStyle(0);
  LegendMolt->SetTextAlign(12); //left align
  LegendMolt->SetTextSize(0.043);
  LegendMolt->SetBorderSize(0);
  LegendMolt->SetHeader("V0M Multiplicity Percentile");
  LegendMolt->SetNColumns(2);

  TLegend *LegendJetRegions=new TLegend(0.54,0.5,0.91,0.71);
  LegendJetRegions->SetFillStyle(0);
  LegendJetRegions->SetTextAlign(12); //left align
  LegendJetRegions->SetTextSize(0.045);
  LegendJetRegions->SetBorderSize(0);

  TLegend *LegendJetRegionsXi;
  if (DataSample==0) LegendJetRegionsXi=new TLegend(0.57,0.27,0.94,0.48);
  else LegendJetRegionsXi=new TLegend(0.32,0.32,0.69,0.53);
  LegendJetRegionsXi->SetFillStyle(0);
  LegendJetRegionsXi->SetTextAlign(12); //left align
  LegendJetRegionsXi->SetTextSize(0.04);
  LegendJetRegionsXi->SetBorderSize(0);

  TLegend *LegendRegionsSist;
  if (type==0)LegendRegionsSist=new TLegend(0.54,0.5,0.91,0.71);
  else LegendRegionsSist=new TLegend(0.54,0.55,0.91,0.76);
  LegendRegionsSist->SetFillStyle(0);
  LegendRegionsSist->SetTextAlign(12); //left align
  LegendRegionsSist->SetTextSize(0.045);
  LegendRegionsSist->SetBorderSize(0);

  TLegend *LegendSB=new TLegend(0.33,0.15,0.53,0.35);
  LegendSB->SetFillStyle(0);
  LegendSB->SetTextAlign(12); //left align
  LegendSB->SetTextSize(0.045);
  LegendSB->SetBorderSize(0);

  //  TLegend *LegendCorr=new TLegend(0.45,0.78,0.8,0.97);
  TLegend *LegendCorr=new TLegend(0.1,0.78,0.35,0.97);
  LegendCorr->SetFillStyle(0);
  LegendCorr->SetTextAlign(12);
  LegendCorr->SetTextSize(0.06);
  LegendCorr->SetBorderSize(0);
  LegendCorr->AddEntry("", "3 < #it{p}_{T}^{trigg} < 15 GeV/#it{c}", "");
  LegendCorr->AddEntry("", Form("%.1f < #it{p}_{T}^{assoc} < %.1f GeV/#it{c}", NPtV0[ptchosen], NPtV0[ptchosen+1]), "");

  TString Directory = "CanvasAC_"+year + "_m"+Smolt[multChosen] + "_v"+ SPtV0[ptchosen];

  canvasPairAcceptance->SetFillColor(0);
  canvasPairAcceptance->SetTickx(1);
  canvasPairAcceptance->SetTicky(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  //  gStyle->SetPalette(ncolors,colors);
  hME->GetXaxis()->SetTitle("#Delta#eta");
  hME->GetYaxis()->SetTitle("#Delta#varphi");
  hME->Draw("lego2");
  LegendALICE->Draw("");
  LegendCorr->Draw("");
  cout <<"\n\e[35m 2D Pair acceptance\e[39m" << endl;
  canvasPairAcceptance->SaveAs(Directory+"_PairAcceptance.pdf");

  TCanvas* canvasThreePads = new TCanvas("canvasThreePads", "canvasThreePads", 1500, 600);
  canvasThreePads->Divide(3, 1);
  canvasThreePads->SetFillColor(0);
  canvasThreePads->SetTickx(1);
  canvasThreePads->SetTicky(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetPalette(55, 0); //55

  canvasThreePads->cd(1);
  gPad->SetTopMargin(0.22);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.04);
  hSE->Scale(1./hSE->GetXaxis()->GetBinWidth(1)/hSE->GetYaxis()->GetBinWidth(1)/1000);
  hSE->SetTitle("");
  hSE->GetZaxis()->SetTitle("#frac{d^{2}#it{N}_{assoc}^{SE}}{d#Delta#eta d#Delta#varphi} #frac{1}{1000}");
  hSE->GetZaxis()->SetTitleSize(0.041);  
  hSE->GetZaxis()->SetTitleOffset(2.2); 
  hSE->GetZaxis()->SetNdivisions(6);
  hSE->GetXaxis()->SetTitle("#Delta#eta");
  hSE->GetYaxis()->SetTitle("#Delta#varphi");
  //  hSE->GetXaxis()->SetRangeUser(-2, 2);
  if (type==1)  hSE->Rebin2D(2);
  hSE->DrawClone("SURF2");
  LegendALICE->Draw("");

  canvasThreePads->cd(2);
  gPad->SetTopMargin(0.22);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.04);
  hME->GetZaxis()->SetTitle("#varepsilon_{pair}");
  hME->GetZaxis()->SetTitleSize(0.06); 
  hME->GetZaxis()->SetTitleOffset(1.2); 
  //  hME->GetXaxis()->SetRangeUser(-2, 2);
  hME->SetTitle("");
  hME->GetXaxis()->SetTitle("#Delta#eta");
  hME->GetYaxis()->SetTitle("#Delta#varphi");
  if (type==1) {
    hME->Rebin2D(2);
    hME->Scale(1./4);
  }
  hME->DrawClone("SURF2");

  canvasThreePads->cd(3);
  gPad->SetTopMargin(0.22);
  gPad->SetLeftMargin(0.2);
  gPad->SetRightMargin(0.04);
  hAC->Scale(1./hAC->GetXaxis()->GetBinWidth(1)/hAC->GetYaxis()->GetBinWidth(1)/1000);
  hAC->GetZaxis()->SetTitle("#frac{d^{2}#it{N}_{assoc}^{SE}}{d#Delta#eta d#Delta#varphi} #frac{1}{#varepsilon_{pair}} #frac{1}{1000}");
  hAC->GetZaxis()->SetTitleOffset(2.3);
  hAC->GetZaxis()->SetTitleSize(0.041);  
  hAC->GetZaxis()->SetNdivisions(6);
  hAC->GetXaxis()->SetTitle("#Delta#eta");
  hAC->GetYaxis()->SetTitle("#Delta#varphi");
  //  hAC->GetXaxis()->SetRangeUser(-2, 2);
  hAC->SetTitle("");
  if (type==1)  hAC->Rebin2D(2);
  hAC->DrawClone("SURF2");
  LegendCorr->Draw("");
  cout <<"\n\e[35m 2D angular correlation distributions\e[39m" << endl;
  canvasThreePads->SaveAs(Directory+"_AC.pdf");

  TH1F *hProjJetwBulk[3];
  TH1F *hProjBulk[3];
  TH1F *hProjFull[3];
  TH1F *hProjJet[3];

  TString NameJetwBulk = "";
  TString NameBulk = "";
  TString NameFull = "";
  TString NameJet = "";
  TString NameStat = "";
  TString NameSist = "";

  TH1F *hProjBulkStat[3];
  TH1F *hProjFullStat[3];
  TH1F *hProjJetStat[3];
  TH1F *hProjJetStatGrey[3];
  TH1F *hProjJetSistGrey[3];
  TH1F *hProjBulkSist[3];
  TH1F *hProjFullSist[3];
  TH1F *hProjJetSist[3];

  Int_t veff[3]={0};

  TCanvas* canvasProjections = new TCanvas("canvasProjections", "canvasProjections", 1500, 600);
  canvasProjections->Divide(3, 1);
  canvasProjections->SetFillColor(0);
  canvasProjections->SetTickx(1);
  canvasProjections->SetTicky(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetPalette(55, 0); //55

  Int_t icounter=0;

  for (Int_t i=0; i<3; i++){
    if (isProj==0 && i!=2) continue;
    if (isProj==1 && i!=0) continue;
    if (isProj==2 && i!=1) continue;
    //    if (isProj==3 && i==0) continue;
    icounter++;
    for (Int_t v=0; v<3; v++){
      /*
      if (v==0) veff[v] = 1;
      else if (v==1) veff[v] = 3;
      else if (v==2) veff[v] = 5;
      */
      
      if (v==0) veff[v] = 5;
      else if (v==1) {
	if (type==0)	veff[v] = 4;
	else 	veff[v] = 6;
      }
      else if (v==2) veff[v] = 7;
      
      if (i==0){ //dPhi projections without SB subtraction
	NameJetwBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_eta0_AC_phi_BeforeSBSub";
	NameBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_eta1_AC_phi_BeforeSBSub";
	NameFull = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_eta2_AC_phi_BeforeSBSub";
      }
      else if (i==1){ //dPhi projections of SB
	NameJetwBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_eta0_AC_phi_SB_NTrigg";
	NameBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_eta1_AC_phi_SB_NTrigg";
	NameFull = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_eta2_AC_phi_SB_NTrigg";
      }
      else if (i==2){ //final dPhi projections
	NameJetwBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_EffCorr_TrCorr";
	NameBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_Bulk_EffCorr_TrCorr";
	NameFull = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_JetBulk_EffCorr_TrCorr";
      }

      hProjJetwBulk[v] = (TH1F*)inputFile->Get(NameJetwBulk);
      hProjBulk[v] = (TH1F*)inputFile->Get(NameBulk);
      hProjFull[v] = (TH1F*)inputFile->Get(NameFull);

      if (!hProjJetwBulk[v]) {cout << "no hProjJetwBulk" << endl; return; }
      if (!hProjBulk[v]) {cout << "no hProjBulk" << endl; return; }
      if (!hProjFull[v]) {cout << "no hProjFull" << endl; return; }

      hProjJetwBulk[v]->SetName(NameJetwBulk + Form("_%i", i)); 
      hProjBulk[v]->SetName(NameBulk + Form("_%i", i)); 
      hProjFull[v]->SetName(NameFull + Form("_%i", i)); 

      canvasProjections->cd(v+1);
      gPad->SetLeftMargin(0.27);
      gPad->SetRightMargin(0.02);
      /*
	gPad->SetLeftMargin(0.27);
	if (v!=0)     gPad->SetLeftMargin(0);
	gPad->SetRightMargin(0);
					       */
      hProjJetwBulk[v]->Scale(1./hProjJetwBulk[v]->GetXaxis()->GetBinWidth(1));
      StyleHisto(hProjJetwBulk[v], 0, 1.25* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[0], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      if (isProj==3){
	if (i==2) {
	  hProjJetwBulk[v]->SetMarkerColor(634);
	  hProjJetwBulk[v]->SetLineColor(634);
	}
	else if (i==1) {
	  hProjJetwBulk[v]->SetMarkerColor(1);
	  hProjJetwBulk[v]->SetLineColor(1);
	}
      }     
      hProjJetwBulk[v]->GetYaxis()->SetTitleSize(0.05);  
      hProjJetwBulk[v]->GetYaxis()->SetNdivisions(6);
      if (isProj!=3) hProjJetwBulk[v]->DrawClone("same e");

      hProjBulk[v]->Scale(1./hProjBulk[v]->GetXaxis()->GetBinWidth(1));
      //StyleHisto(hProjBulk[v], 0, 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[1], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      StyleHisto(hProjBulk[v], 0, 1.3* hProjFull[0]->GetMaximum(hProjFull[0]->GetMaximumBin()),  color[1], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      if (isProj==3){
	if (i==2) {
	  hProjBulk[v]->SetMarkerColor(829);
	  hProjBulk[v]->SetLineColor(829);
	}
	else if (i==1) {
	  hProjBulk[v]->SetMarkerColor(1);
	  hProjBulk[v]->SetLineColor(1);
	}
      }     
      hProjBulk[v]->GetYaxis()->SetTitleSize(0.05);  
      hProjBulk[v]->GetYaxis()->SetNdivisions(6);
      if (isProj!=3) hProjBulk[v]->DrawClone("same e");

      hProjFull[v]->Scale(1./hProjFull[v]->GetXaxis()->GetBinWidth(1));
      //StyleHisto(hProjFull[v], 0, 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[2], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      StyleHisto(hProjFull[v], 0, 1.3* hProjFull[v]->GetMaximum(hProjFull[0]->GetMaximumBin()),  color[2], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      if (isProj==3){
	if (i==2) {
	  hProjFull[v]->SetMarkerColor(867);
	  hProjFull[v]->SetLineColor(867);
	}
	else if (i==1) {
	  hProjFull[v]->SetMarkerColor(1);
	  hProjFull[v]->SetLineColor(1);
	}
      }     
      hProjFull[v]->GetYaxis()->SetTitleSize(0.05);  
      hProjFull[v]->GetYaxis()->SetNdivisions(6);
      hProjFull[v]->DrawClone("same e");

      TLegend *LegendPt = new TLegend(0.3, 0.8, 0.5, 0.84);
      LegendPt->AddEntry("", Form("%.1f < #it{p}_{T}^{assoc} < %.1f GeV/#it{c}", NPtV0[veff[v]],  NPtV0[veff[v]+1]), "");
      LegendPt->SetTextSize(0.05);
      LegendPt->Draw("same");
      if (v==0 && icounter==1){
	LegendRegions->AddEntry(hProjJetwBulk[v], "#color[628]{|#Delta#eta| < 0.86}", "pl");
	LegendRegions->AddEntry(hProjBulk[v], "#color[418]{0.86 < |#Delta#eta| < 1.2}", "pl");
	LegendRegions->AddEntry(hProjFull[v], "#color[601]{|#Delta#eta| < 1.2}", "pl");
      }
      if (isProj==3) {
	if (v==0 && i==0) LegendSB->AddEntry(hProjFull[v], "h-K_{S}^{0} correlations", "pl");
	if (v==0 && i==1) LegendSB->AddEntry(hProjFull[v], "h-fake K_{S}^{0} correlations (SB)", "pl");
	if (v==0 && i==2) LegendSB->AddEntry(hProjFull[v], "h-K_{S}^{0} correlations - SB", "pl");
	if (v==0)   LegendALICEwSB->Draw("same");
	else if (v==1) LegendSB->Draw("same");
      }
      else {
	if (v==0)   LegendALICEProj->Draw("same");
	else if (v==1) LegendRegions->Draw("same");
      }
    }
  }

  TString ScanvasProjections = "";
  if (isProj==0) ScanvasProjections = "Canvas"+year+"_dPhi.pdf";
  else if (isProj==1) ScanvasProjections = "Canvas"+year+"_dPhi_BeforeSBSub.pdf";
  else if (isProj==2) ScanvasProjections = "Canvas"+year+"_dPhi_SB_dPhi.pdf";
  else if (isProj==3) ScanvasProjections = "Canvas"+year+"_dPhi_SBAndNot_dPhi.pdf";
  cout <<"\n\e[35m 1D dPhi projections (from Peak of inv. mass., or from SB, or difference between the 2 \e[39m" << endl;
  canvasProjections->SaveAs(ScanvasProjections);

  TCanvas* canvasJetProjections = new TCanvas("canvasJetProjections", "canvasJetProjections", 1500, 600);
  canvasJetProjections->Divide(3, 1);
  canvasJetProjections->SetFillColor(0);
  canvasJetProjections->SetTickx(1);
  canvasJetProjections->SetTicky(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetPalette(55, 0); //55

  Float_t LowLimitJet = -0.001;
  for (Int_t v=0; v<3; v++){
    /*
      if (v==0) veff[v] = 1;
      else if (v==1) veff[v] = 3;
      else if (v==2) veff[v] = 5;
    */
      
    if (v==0) veff[v] = 3;
    else if (v==1) {
      if (type==0)	veff[v] = 4;
      else 	veff[v] = 5;
    }
    else if (v==2) veff[v] = 7;

    NameJetwBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_EffCorr_TrCorr";
    NameBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_Bulk_EffCorr_TrCorr";
    NameJet = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_BulkSub_EffCorr_TrCorr";

    hProjJetwBulk[v] = (TH1F*)inputFile->Get(NameJetwBulk);
    hProjBulk[v] = (TH1F*)inputFile->Get(NameBulk);
    hProjJet[v] = (TH1F*)inputFile->Get(NameJet);

    if (!hProjJetwBulk[v]) {cout << "no hProjJetwBulk" << endl; return; }
    if (!hProjBulk[v]) {cout << "no hProjBulk" << endl; return; }
    if (!hProjJet[v]) {cout << "no hProjJet" << endl; return; }

    canvasJetProjections->cd(v+1);
    gPad->SetLeftMargin(0.27);
    gPad->SetRightMargin(0.02);

    //    cout << "BC: " << hProjJetwBulk[v]->GetBinContent(2) << endl;
    hProjJetwBulk[v]->Rebin(2);
    //    cout << "BC reb: " << hProjJetwBulk[v]->GetBinContent(2) << endl;
    //    cout << hProjJetwBulk[v]->GetXaxis()->GetBinWidth(1) << endl;
    hProjJetwBulk[v]->Scale(1./hProjJetwBulk[v]->GetXaxis()->GetBinWidth(1));
    //    cout << "BC scaled: " << hProjJetwBulk[v]->GetBinContent(2) << endl;
    StyleHisto(hProjJetwBulk[v], LowLimitJet , 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  634 , 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjJetwBulk[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjJetwBulk[v]->GetYaxis()->SetNdivisions(6);
    hProjJetwBulk[v]->DrawClone("same e");

    //    cout << "BC: " << hProjBulk[v]->GetBinContent(2) << endl;
    hProjBulk[v]->Rebin(2);
    //    cout << "BC reb: " << hProjBulk[v]->GetBinContent(2) << endl;
    //    cout << hProjBulk[v]->GetXaxis()->GetBinWidth(1) << endl;
    hProjBulk[v]->Scale(1./hProjBulk[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjBulk[v], LowLimitJet, 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[1], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjBulk[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjBulk[v]->GetYaxis()->SetNdivisions(6);
    hProjBulk[v]->DrawClone("same e");

    //    cout << "BC: " << hProjJet[v]->GetBinContent(2) << endl;
    //    cout << hProjJet[v]->GetXaxis()->GetBinWidth(1) << endl;
    hProjJet[v]->Scale(1./hProjJet[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjJet[v], LowLimitJet, 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[0], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjJet[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjJet[v]->GetYaxis()->SetNdivisions(6);
    hProjJet[v]->DrawClone("same e");

    TLegend *LegendPt = new TLegend(0.3, 0.75, 0.5, 0.87);
    LegendPt->AddEntry("", Form("%.1f < #it{p}_{T}^{assoc} < %.1f GeV/#it{c}", NPtV0[veff[v]],  NPtV0[veff[v]+1]), "");
    LegendPt->SetTextSize(0.05);
    LegendPt -> Draw("same");
    if (v==0){
      LegendJetRegions->AddEntry(hProjJetwBulk[v], "#color[634]{|#Delta#eta| < 0.86}", "pl");
      LegendJetRegions->AddEntry(hProjBulk[v], "#color[418]{0.86 < |#Delta#eta| < 1.2}", "pl");
      //      LegendJetRegions->AddEntry(hProjJet[v], "#color[628]{|#Delta#eta| < 0.86 - 0.86 < |#Delta#eta| < 1.2}", "pl");
      //      LegendJetRegions->AddEntry(hProjJet[v], "#color[628]{Near-side jet}", "pl");
      LegendJetRegions->AddEntry(hProjJet[v], "#color[628]{Toward-leading}", "pl");
    }
    if (v==0)   LegendALICEJet->Draw("same");
    else if (v==1) LegendJetRegions->Draw("same");
  }
 
  TString ScanvasJetProjections = "";
  ScanvasJetProjections = "Canvas"+year+"_dPhi_JetSub.pdf";
  cout <<"\n\e[35mJet dPhi projections (before and after ooj subtraction) \e[39m" << endl;
  canvasJetProjections->SaveAs(ScanvasJetProjections);

  TCanvas* canvasXiJetProjections;
  canvasXiJetProjections = new TCanvas("canvasXiJetProjections", "canvasXiJetProjections", 1500, 600);
  if (type==1){
    if (DataSample==0)    Smolt[multChosen]= "10-30";
    else     Smolt[multChosen]= "0.01-0.05";
    canvasXiJetProjections->Divide(3, 1);
    canvasXiJetProjections->SetFillColor(0);
    canvasXiJetProjections->SetTickx(1);
    canvasXiJetProjections->SetTicky(1);
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetPalette(55, 0); //55

    Float_t LowLimitJet = -0.0002;
    Int_t vcounter=-1;
    TF1 * LineAtZero = new TF1("pol0", "pol0", -TMath::Pi()/2, 3./2*TMath::Pi());
    LineAtZero->SetParameter(0, 0);
    LineAtZero->SetLineColor(kGray+3);
    LineAtZero->SetLineStyle(3);
    for (Int_t v=0; v<3; v++){
      //      if (DataSample==1 && v==0) continue;
      vcounter++;
      if (v==0) veff[v] = 3;
      else if (v==1) veff[v] = 4;
      else if (v==2) veff[v] = 6;

      cout << Smolt[multChosen] << endl;
      cout << "SPtV0[veff[v]] " << SPtV0[veff[v]]<< endl;

      if (NPtV0[veff[v]] < 2.5){
	cout << "Hola ! " << endl;
	NameJetwBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_EffCorr_Rebin";
	NameBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_Bulk_EffCorrScaledAllMultSmoothed";
	NameJet = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_BulkSub_EffCorr_RebSmooth";
	hProjJetwBulk[v] = (TH1F*)inputFileJetXi->Get(NameJetwBulk);
	hProjBulk[v] = (TH1F*)inputFileJetXi->Get(NameBulk);
	hProjJet[v] = (TH1F*)inputFileJetXi->Get(NameJet);
      }
      else {
	NameJetwBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_EffCorr_TrCorr";
	NameBulk = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_Bulk_EffCorr_TrCorr";
	NameJet = "ME_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]]+"_AC_phi_V0Sub_BulkSub_EffCorr_TrCorr";
	hProjJetwBulk[v] = (TH1F*)inputFile->Get(NameJetwBulk);
	hProjBulk[v] = (TH1F*)inputFile->Get(NameBulk);
	hProjJet[v] = (TH1F*)inputFile->Get(NameJet);
      }

      if (!hProjJetwBulk[v]) {cout << "no hProjJetwBulk Jet Xi" << endl; return; }
      if (!hProjBulk[v]) {cout << "no hProjBulk Jet Xi" << endl; return; }
      if (!hProjJet[v]) {cout << "no hProjJet Jet xi" << endl; return; }

      /*
      if (DataSample==0) canvasXiJetProjections->cd(v+1);
      else if (DataSample==1) canvasXiJetProjections->cd(v);
      */
      canvasXiJetProjections->cd(v+1);
      gPad->SetLeftMargin(0.27);
      gPad->SetRightMargin(0.02);

      if (NPtV0[veff[v]] >= 2.5) hProjJetwBulk[v]->Rebin(2);
      hProjJetwBulk[v]->Scale(1./hProjJetwBulk[v]->GetXaxis()->GetBinWidth(1));
      StyleHisto(hProjJetwBulk[v], LowLimitJet , 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  634 , 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      hProjJetwBulk[v]->GetYaxis()->SetTitleSize(0.05);  
      hProjJetwBulk[v]->GetYaxis()->SetNdivisions(6);
      hProjJetwBulk[v]->DrawClone("same e");

      if (NPtV0[veff[v]] >= 2.5) hProjBulk[v]->Rebin(2);
      hProjBulk[v]->Scale(1./hProjBulk[v]->GetXaxis()->GetBinWidth(1));
      StyleHisto(hProjBulk[v], LowLimitJet, 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[1], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      hProjBulk[v]->GetYaxis()->SetTitleSize(0.05);  
      hProjBulk[v]->GetYaxis()->SetNdivisions(6);
      hProjBulk[v]->DrawClone("same e");

      hProjJet[v]->Scale(1./hProjJet[v]->GetXaxis()->GetBinWidth(1));
      StyleHisto(hProjJet[v], LowLimitJet, 1.3* hProjJetwBulk[v]->GetMaximum(hProjJetwBulk[v]->GetMaximumBin()),  color[0], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
      hProjJet[v]->GetYaxis()->SetTitleSize(0.05);  
      hProjJet[v]->GetYaxis()->SetNdivisions(6);
      hProjJet[v]->DrawClone("same e");

      LineAtZero->DrawClone("same e");

      TLegend *LegendPt = new TLegend(0.3, 0.75, 0.5, 0.87);
      LegendPt->AddEntry("", Form("%.1f < #it{p}_{T}^{assoc} < %.1f GeV/#it{c}", NPtV0[veff[v]],  NPtV0[veff[v]+1]), "");
      LegendPt->SetTextSize(0.05);
      LegendPt -> Draw("same");
      if (vcounter==0){
	LegendJetRegionsXi->AddEntry(hProjJetwBulk[v], "#color[634]{|#Delta#eta| < 0.86}", "pl");
	//	LegendJetRegionsXi->AddEntry(hProjBulk[v], "#color[418]{0.86 < |#Delta#eta| < 1.2 in V0M 0-0.01%}", "pl");
	LegendJetRegionsXi->AddEntry(hProjBulk[v], "#color[418]{0.86 < |#Delta#eta| < 1.2}", "pl");
	//	LegendJetRegionsXi->AddEntry(hProjJet[v], "#color[628]{Near-side jet}", "pl");
	LegendJetRegionsXi->AddEntry(hProjJet[v], "#color[628]{Toward-leading}", "pl");
      }
      if (vcounter==0)   LegendALICEJetXi->Draw("same");
      else if (vcounter==1) LegendJetRegionsXi->Draw("same");
    }
 
    TString ScanvasXiJetProjections = "";
    ScanvasXiJetProjections = "Canvas"+year+"_dPhi_JetSubXi.pdf";
    cout <<"\n\e[35mJet dPhi projections (before and after ooj subtraction) for Xi \e[39m" << endl;
    canvasXiJetProjections->SaveAs(ScanvasXiJetProjections);
  }

  TCanvas* canvasProjectionsSyst;
  canvasProjectionsSyst = new TCanvas("canvasProjectionsSyst", "canvasProjectionsSyst", 1500, 600);
  if (type==0){
    if (DataSample==0)    Smolt[multChosen]= "_all";
    else     Smolt[multChosen]= "0.01-0.05";
  }
  else {
    if (DataSample==0)    Smolt[multChosen]= "_all";
    else     Smolt[multChosen]= "0-0.1";
  }
  canvasProjectionsSyst->Divide(3, 1);
  canvasProjectionsSyst->SetFillColor(0);
  canvasProjectionsSyst->SetTickx(1);
  canvasProjectionsSyst->SetTicky(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetPalette(55, 0); //55

  Int_t  vcounter=-1;
  if (type==1) LowLimitJet = -0.0002;
  TF1 * LineAtZero = new TF1("pol0", "pol0", -TMath::Pi()/2, 3./2*TMath::Pi());
  LineAtZero->SetParameter(0, 0);
  LineAtZero->SetLineColor(kGray+3);
  LineAtZero->SetLineStyle(3);
  for (Int_t v=0; v<3; v++){
    vcounter++;
    if (type==0){
      if (v==0) veff[v] = 3;
      else if (v==1) veff[v] = 4;
      else if (v==2) veff[v] = 7;
    }
    else {
      if (v==0) veff[v] = 5;
      else if (v==1) veff[v] = 6;
      else if (v==2) veff[v] = 7;
    }

    NameStat = "PhiDistr_solostat_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]];
    NameSist = "PhiDistr_solosist_m"+Smolt[multChosen]+"_v"+ SPtV0[veff[v]];

    cout << NameStat << " " << NameSist << endl;
    hProjFullStat[v] = (TH1F*)inputFileSystFull->Get(NameStat);
    hProjFullSist[v] = (TH1F*)inputFileSystFull->Get(NameSist);
    if (!hProjFullStat[v]) {cout << "no hProjFull stat" << endl; return; }
    if (!hProjFullSist[v]) {cout << "no hProjFull syst" << endl; return; }
    hProjFullStat[v] ->SetName(NameStat+"_Full");
    hProjFullSist[v] ->SetName(NameSist+"_Full");

    hProjBulkStat[v] = (TH1F*)inputFileSystBulk->Get(NameStat);
    hProjBulkSist[v] = (TH1F*)inputFileSystBulk->Get(NameSist);
    if (!hProjBulkStat[v]) {cout << "no hProjBulk stat" << endl; return; }
    if (!hProjBulkSist[v]) {cout << "no hProjBulk syst" << endl; return; }
    hProjBulkStat[v] ->SetName(NameStat+"_Bulk");
    hProjBulkSist[v] ->SetName(NameSist+"_Bulk");

    hProjJetStat[v] = (TH1F*)inputFileSystJet->Get(NameStat);
    hProjJetSist[v] = (TH1F*)inputFileSystJet->Get(NameSist);
    if (!hProjJetStat[v]) {cout << "no hProjJet stat" << endl; return; }
    if (!hProjJetSist[v]) {cout << "no hProjJet syst" << endl; return; }
    hProjJetStat[v] ->SetName(NameStat+"_Jet");
    hProjJetSist[v] ->SetName(NameSist+"_Jet");

    canvasProjectionsSyst->cd(v+1);
    gPad->SetLeftMargin(0.27);
    gPad->SetRightMargin(0.02);

    hProjFullStat[v]->Rebin(2);
    hProjFullStat[v]->Scale(1./hProjFullStat[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjFullStat[v], LowLimitJet , 1.45* hProjFullStat[v]->GetMaximum(hProjFullStat[v]->GetMaximumBin()),  601, 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjFullStat[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjFullStat[v]->GetYaxis()->SetNdivisions(6);
    hProjFullStat[v]->DrawClone("same e");
    hProjFullSist[v]->Rebin(2);
    hProjFullSist[v]->Scale(1./hProjFullSist[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjFullSist[v], LowLimitJet , 1.45* hProjFullSist[v]->GetMaximum(hProjFullSist[v]->GetMaximumBin()), 601 , 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjFullSist[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjFullSist[v]->GetYaxis()->SetNdivisions(6);
    hProjFullSist[v]->SetFillStyle(0);
    hProjFullSist[v]->SetFillColorAlpha(color[2], 1);
    hProjFullSist[v]->DrawClone("same e2");

    hProjBulkStat[v]->Rebin(2);
    hProjBulkStat[v]->Scale(1./hProjBulkStat[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjBulkStat[v], LowLimitJet, 1.45* hProjFullStat[v]->GetMaximum(hProjFullStat[v]->GetMaximumBin()),  color[1], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjBulkStat[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjBulkStat[v]->GetYaxis()->SetNdivisions(6);
    hProjBulkStat[v]->DrawClone("same e");
    hProjBulkSist[v]->Rebin(2);
    hProjBulkSist[v]->Scale(1./hProjBulkSist[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjBulkSist[v], LowLimitJet, 1.45* hProjFullSist[v]->GetMaximum(hProjFullSist[v]->GetMaximumBin()),  color[1], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjBulkSist[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjBulkSist[v]->GetYaxis()->SetNdivisions(6);
    hProjBulkSist[v]->SetFillStyle(0);
    hProjBulkSist[v]->SetFillColorAlpha(color[1], 1);
    hProjBulkSist[v]->DrawClone("same e2");

    hProjJetStat[v]->Scale(1./hProjJetStat[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjJetStat[v], LowLimitJet, 1.45* hProjFullStat[v]->GetMaximum(hProjFullStat[v]->GetMaximumBin()),  color[0], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjJetStat[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjJetStat[v]->GetYaxis()->SetNdivisions(6);
    hProjJetStat[v]->DrawClone("same e");
    hProjJetSist[v]->Scale(1./hProjJetSist[v]->GetXaxis()->GetBinWidth(1));
    StyleHisto(hProjJetSist[v], LowLimitJet, 1.45* hProjFullSist[v]->GetMaximum(hProjFullSist[v]->GetMaximumBin()),  color[0], 33,  "#Delta#varphi",  "#frac{1}{#it{N}_{trigg}} #frac{1}{#Delta#eta} #frac{d#it{N}_{assoc}}{d#Delta#varphi}", "", 0,  0,  0, 1,  2.5, 1.2);
    hProjJetSist[v]->GetYaxis()->SetTitleSize(0.05);  
    hProjJetSist[v]->GetYaxis()->SetNdivisions(6);
    hProjJetSist[v]->SetFillStyle(0);
    hProjJetSist[v]->SetFillColorAlpha(color[0], 1);
    hProjJetSist[v]->DrawClone("same e2");

    LineAtZero->DrawClone("same e");

    if (v==0 && icounter==1){
      LegendRegionsSist->AddEntry(hProjJetStat[v], "#color[628]{Toward-leading}", "pl");
      LegendRegionsSist->AddEntry(hProjBulkStat[v], "#color[418]{0.86 < |#Delta#eta| < 1.2}", "pl");
      LegendRegionsSist->AddEntry(hProjFullStat[v], "#color[601]{|#Delta#eta| < 1.2}", "pl");

      hProjJetStatGrey[v] = (TH1F*)hProjJetStat[v]->Clone("GreyProjStat");
      hProjJetStatGrey[v]->SetLineColor(1);
      hProjJetStatGrey[v]->SetMarkerColor(1);
      hProjJetSistGrey[v] = (TH1F*)hProjJetSist[v]->Clone("GreyProjSist");
      hProjJetSistGrey[v]->SetLineColor(1);
      hProjJetSistGrey[v]->SetMarkerColor(1);
      LegendSist->AddEntry(hProjJetStatGrey[v], "stat. error", "pe");
      LegendSist->AddEntry(hProjJetSistGrey[v], "syst. error", "ef");
    }

    TLegend *LegendPt = new TLegend(0.3, 0.8, 0.5, 0.87);
    LegendPt->AddEntry("", Form("%.1f < #it{p}_{T}^{assoc} < %.1f GeV/#it{c}", NPtV0[veff[v]],  NPtV0[veff[v]+1]), "");
    LegendPt->SetTextSize(0.05);
    LegendPt -> Draw("same");

    if (v==0)   LegendALICEJet->Draw("same");
    else if (v==2) LegendRegionsSist->Draw("same");
    else if (v==1)    LegendSist->Draw("");
  }
 
  TString ScanvasProjectionsSyst = "";
  ScanvasProjectionsSyst = "Canvas"+year+"_dPhi_WithErrors.pdf";
  cout <<"\n\e[35mdPhi projections with errors \e[39m" << endl;
  canvasProjectionsSyst->SaveAs(ScanvasProjectionsSyst);

  TCanvas* canvasRegions = new TCanvas("canvasRegions", "canvasRegions", 1000, 700);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.15);
  gStyle->SetPalette(1, 0); //55
  hAC->GetXaxis()->SetTitleOffset(1.3);
  hAC->GetYaxis()->SetTitleOffset(1);
  hAC->DrawClone("colz");
  cout <<"\n\e[35m 2D angular correlation distribution FLAT\e[39m" << endl;
  canvasRegions->SaveAs(Directory+"_Regions.pdf");

  TCanvas* canvasMERatios = new TCanvas("canvasMERatios", "canvasMERatios", 1500, 600);
  canvasMERatios->Divide(3, 1);
  canvasMERatios->SetFillColor(0);
  canvasMERatios->SetTickx(1);
  canvasMERatios->SetTicky(1);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetPalette(55, 0); //55

  TString NameMENorm = ""; 
  TH2F *  hMENorm2D[nummolt+1][numPtV0];
  TH1F *  hMENorm[nummolt+1][numPtV0];
  TH1F *  hMENormRatio[nummolt+1][numPtV0];
  TH2F *  hMENormAllMult2D[numPtV0];
  TH1F *  hMENormAllMult[numPtV0];

  Int_t nummoltMax = 5;
  Int_t ColorMult[nummoltMax+1] = {1, 801, 628, 867, 909,  881};
  TString SmoltLegend[nummoltMax+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  if (DataSample==1){
    SmoltLegend[0] = "0-0.01 %";
    SmoltLegend[1] = "0.01-0.05 %";
    SmoltLegend[2] = "0.05-0.1 %";
    SmoltLegend[5] = "0-0.1 %";
    Smolt[0] = "0-0.01";
    Smolt[1] = "0.01-0.05";
    Smolt[2] = "0.05-0.1";
    Smolt[5] = "0-0.1";
  }
  for (Int_t m=0; m<=nummoltMax; m++){
    if (isMERatio && m==nummoltMax) continue;
    if (type==1 && (m==4 || m==3)) continue;
    if (DataSample==1 && m>2 && m <5) continue;
    for (Int_t v=0; v<3; v++){      
      if (v==0) {
	if (type==0)	veff[v] = 2;
	else 	veff[v] = 3;
      }
      else if (v==1) {
	if (type==0)	veff[v] = 4;
	else 	veff[v] = 6;
      }
      else if (v==2) veff[v] = 7;

      NameMENorm = "ME_m"+Smolt[5]+"_v"+ SPtV0[veff[v]]+"_norm";
      hMENormAllMult2D[v] = (TH2F*)inputFile->Get(NameMENorm);
      if (!hMENormAllMult2D[v]) {cout << "no ME all mult" << endl; return; }
      hMENormAllMult2D[v]->SetName("ME_m"+Smolt[5]+"_v"+ SPtV0[veff[v]]+"_norm_Denominator");
      NameMENorm = "ME_m"+Smolt[m]+"_v"+ SPtV0[veff[v]]+"_norm";
      hMENorm2D[m][v] = (TH2F*)inputFile->Get(NameMENorm);
      if (!hMENorm2D[m][v]) {cout << "no " <<  NameMENorm << endl; return; }
      hMENorm[m][v] = (TH1F*) hMENorm2D[m][v]->ProjectionX();
      hMENorm[m][v]->Scale(1./hMENorm2D[m][v]->GetYaxis()->GetNbins());
      hMENormAllMult[v] = (TH1F*) hMENormAllMult2D[v]->ProjectionX();
      hMENormAllMult[v]->Scale(1./hMENormAllMult2D[v]->GetYaxis()->GetNbins());

      hMENormRatio[m][v] = (TH1F*) hMENorm[m][v]->Clone(NameMENorm+ "MultRatio");
      if (isMERatio) hMENormRatio[m][v]->Divide(hMENormAllMult[v]);

      canvasMERatios->cd(v+1);
      gPad->SetLeftMargin(0.18);
      gPad->SetRightMargin(0.04);
      if (isMERatio) {
	if (type==0)	StyleHisto(hMENormRatio[m][v], 0.5+10e-3, 1.5-10e-3, ColorMult[m], 33,  "#Delta#eta",  "Pair acceptance ratio to 0-100%", "", 0,  0,  0, 1, 1.8, 1.2);
	else 	StyleHisto(hMENormRatio[m][v], 0.5+10e-3, 1.5-10e-3, ColorMult[m], 33,  "#Delta#eta",  "Pair acceptance ratio to 0-100%", "", 0,  0,  0, 1, 1.8, 1.2);
      }
      else       StyleHisto(hMENormRatio[m][v], 0, 1.2, ColorMult[m], 33,  "#Delta#eta",  "Pair acceptance", "", 0,  0,  0, 1,1.5, 1.2);
      hMENormRatio[m][v]->GetXaxis()->SetRangeUser(-1.5 + 10e-2, 1.5-10e-2);
      hMENormRatio[m][v]->GetYaxis()->SetTitleSize(0.05);  

      hMENormRatio[m][v]->DrawClone("same e");

      TLegend *LegendPt = new TLegend(0.3, 0.75, 0.5, 0.87);
      LegendPt->AddEntry("", Form("%.1f < #it{p}_{T}^{assoc} < %.1f GeV/#it{c}", NPtV0[veff[v]],  NPtV0[veff[v]+1]), "");
      LegendPt->SetTextSize(0.05);
      LegendPt -> Draw("same");
      if (v==0)  LegendMolt->AddEntry(hMENormRatio[m][v], SmoltLegend[m], "pl");
      if (v==0)   LegendALICEMERatio->Draw("same");
      else if (v==1) LegendMolt->Draw("same");
    }
  }

  if (isMERatio)  canvasMERatios->SaveAs("Canvas"+year+"_MEMultRatios.pdf");
  else  canvasMERatios->SaveAs("Canvas"+year+"_MEMult.pdf");

  TFile * outputf = new TFile("outputf.root", "RECREATE");
  canvasThreePads->Write();
  canvasProjections->Write();
  canvasJetProjections->Write();
  if (type==1)  canvasXiJetProjections->Write();
  canvasProjectionsSyst->Write();
  canvasMERatios->Write();
  outputf->Close();

  cout << "\nPartendo dai file \n" << stringInput << endl;
}
