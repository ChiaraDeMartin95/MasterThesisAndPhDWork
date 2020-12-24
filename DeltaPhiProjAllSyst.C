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
#include "TFitResult.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

void StyleHistoYield(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset){
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


void StyleHisto(TH1D *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title){
  histo->GetYaxis()->SetRangeUser(Low, Up);
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

void DeltaPhiProjAllSyst( Int_t type=8,Int_t numsysV0index=400, Int_t numsysTriggerindex=99, Int_t sysPhi =0, Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Bool_t isMC=0,   Int_t israp=0,TString year=/*"1617_hK0s"/*"AllMC_hXi"/*"2018f1_extra_hK0s"/"2016k_hK0s"*/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"2016k_New"*//*"Run2DataRed_hXi"/"2016kehjl_hK0s"*/,  TString Path1 =""/*"_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/, Bool_t isEfficiency=0,   TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=1, Int_t MultBinning=0, Int_t PtBinning=0, TString FitFixed="Boltzmann"/*"Fermi-Dirac"*/,   Int_t sys=0, Bool_t ANote=1){


  if (type==0){
    PtBinning=1;
    if (!isMC) year= "1617_hK0s";
    else  year= "1617MC_hK0s";
  }
  else if (type==8){
    PtBinning=0;
    if (!isMC) year="Run2DataRed_MECorr_hXi";
    else  year=  "AllMC_hXi";
  }

  TString RegionType[3] = {"Jet", "Bulk", "Inclusive"};

  if (ishhCorr && type!=0){
    type=0;
  }// {cout << " for hh correlation type should be set to zero in order to get the correct files " << endl; return;}

  const Int_t numregions=3;
  const   Int_t NumberTypeAnalysis=11;
  cout << "here's the meaning of different values of TypeAnalysis:" << endl;
  cout << 0 << "in-jet production " << endl;
  cout << 1 << "out-of-jet production  (Delta Phi between 1 and 2)" << endl;
  cout << 2 << "inclusive production (from JetBulkEffCorr)" << endl;
  cout << 3 << "inclusive production not scaled by DeltaEta and DeltaPhi (for comparison with published data) (from JetBulkEffCorr)" << endl;

  if (year != "2018f1_extra" && year != "2016k" && year != "2018f1_extra_onlyTriggerWithHighestPt" && year != "2016k_onlyTriggerWithHighestPt") {
    //    cout << "output file should be changed: it must include the year name to avoid overwriting output files " << endl;
    //    return;
  } 

  gStyle->SetOptStat(0);

  TString isMCOrData[3] = {"_Data", "_MC", "_DataMC"};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;

  const Int_t numSysV0 = 7;//num eff sistematici associati alla selezione delle particelle associate
  const Int_t numSysV0hh = 3; //num eff sistematici associati alla selezione delle particelle associate (nel caso ishhCorr)
  const Int_t numSyst  = 13; //all systematics together(numsysv0 + numsystphi + systematic related to bin counting region (1))
  const Int_t numSysthh  = 9;
  const Int_t numsystPhi  = 5; //includo tutti i sistematici anche se i primi 3 non vengono usati in analisi hh (si tratta sia dei sistematici legati alla scelta delle regioni di massa invariante ('centrale' e sidebands) (non per hh) che dei sistematici legati alla scelta delle regioni in DeltaEta di jet e out of jet) 
  cout << "ok " << endl;
  Int_t PtV0Min = 0; //0 el
  if (type>0 || ishhCorr )   PtV0Min =1;
  if (type==0 && PtBinning!=0) PtV0Min=1;
  Int_t numSysV0Global = 0; 
  Int_t numSystGlobal  = 0;
  Int_t numsystPhiGlobal  = 0;
  if (!ishhCorr){
    /*    if (type==0){
	  numSysV0Global = numSysV0; 
	  numSystGlobal  = numSyst;
	  numsystPhiGlobal  = numsystPhi;
	  }*/
    //    if (type!=0){    
    numSysV0Global = 1; 
    numSystGlobal  = 2;
    numsystPhiGlobal  = 0;
    //    }
  }
  if (ishhCorr){
    numSysV0Global = numSysV0hh;
    numSystGlobal  = numSyst;
    numsystPhiGlobal  = numsystPhi;
  }

  cout << "ok " << endl;
  TString JetOrBulk[NumberTypeAnalysis]={"Jet", "Bulk", "All", "AllNS", "AwaySide", "AllButJet", "JetBFit", "JetZYAM", "ASBFit", "ASZYAM", "JetNS"};
  TString DataOrMC[2]={"Data", "MC"};

  cout << "ok " << endl;
  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi", "Omega"};
  cout << "ok " << endl;
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
  TString titleYieldYType[9]={"#it{N}_{K^{0}_{S}}/#it{N}_{trigg} 1/#Delta#it{#eta}", "", "", "", "", "", "", "","#it{N}_{#Xi}/#it{N}_{trigg} 1/#Delta#it{#eta}"};
  TString SRegionType[3] = {"Near-side Jet", "Out-of-jet", "Inclusive"};
  TString NameP[9]={"h#minusK_{S}^{0}","", "", "", "", "", "", "",  "h#minus#Xi"};
  TString sRegion[3]={"#color[628]{Near-side jet}","#color[418]{Out-of-jet}","#color[600]{Inclusive}"};
  TString sRegionDEta[3]={"|#Delta#it{#eta}| < 0.75", "0.75 < |#Delta#it{#eta}| < 1.2", "|#Delta#it{#eta}| < 1.2"};
  TString sRegionBlack[3]={"#color[1]{Near-side jet}","#color[1]{Out-of-jet}","#color[1]{Inclusive}"};
  TString sRegion1[3]={"|#Delta#it{#eta}| < 0.75, |#Delta#it{#varphi}| < 0.85", "0.75 < |#Delta#it{#eta}| < 1.2, 0.85 < #Delta#it{#varphi} < 2", "|#Delta#it{#eta}| < 1.2, -#pi/2 < #Delta#it{#varphi} < 3#pi/2"};
  TString xTitle="#Delta#it{#varphi}";
Int_t  MarkerRegion[3]={20,21,33};
 Float_t  MarkerSize[3]={1,1,1.5};
 Float_t yOffset[10]={1.9,1,1,1,1,1,1,1,1.4,1};
  Float_t SpectrumSup[numtipo][NumberTypeAnalysis]={{0.015,0.2, 0.2,0.8,0.2,0.2,0.015, 0.015, 0.03, 0.03, 0.03},{0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05}, {0.05,0.2, 0.2,0.8,0.2,0.2,0.05, 0.05},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001},{0.001,0.01, 0.01,0.08,0.01,0.01,0.001, 0.001, 0.002, 0.002, 0.002},{0.015,0.01, 0.01,0.08,0.01,0.01,0.05, 0.001}};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};

  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};

  Double_t Nmolt0[nummolt+1]={0,5,10,30,50,100}; 
  Double_t Nmolt1[nummolt+1]={0,1,5,15,30,100}; 
  Double_t Nmolt2[nummolt+1]={0,2,7,15,30,100}; 
  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "_all"};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "_all"};
  TString SmoltLegend0[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString SmoltLegend1[nummolt+1]={"0-1 %", "1-5 %", "5-15 %", "15-30 %", "30-100 %", "0-100 %"};
  TString SmoltLegend2[nummolt+1]={"0-2 %", "2-7 %", "7-15 %", "15-30 %", "30-100 %", "0-100 %"};

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

  }



  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", "", ""};  
  //TString SPtV0[numPtV0]={"", "0.5-1", "0.5-1",  "1-1.5","1.5-2", "2-3", "3-4", "4-8"};
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0-0.5"};
    SPtV0[1]={"0.5-1"};
  }

  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }
  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", "", ""};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }


  TString SPtV01[numPtV0]={"0-0.1","0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0, 0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0-0.1", "0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0,0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};
  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1 || PtBinning==2) numPtV0Max = numPtV0;
  else numPtV0Max = numPtV0-2;

  if (PtBinning==1){
    for(Int_t v=0; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=0; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }
  cout << "ok " << endl;
  TString SErrorSpectrum[3]={"stat.","syst. uncorr.","syst. corr."};
  TString SSystJet[2]={"BC [-1.0, 1.0]", "BC [-1.2, 1.2]"};
  TString SSystBulk[3]={"BC [1.0, 2.0]", "BC [2.0, 4.28]", "BC [1.0, 4.28]"};
  TString SFit[3]={"exponential", "power law", "retta"};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};
  Int_t inutile[numSyst]=     {0,1,    1, 0,  1,  0,  1,  1,  1,  0,  1,  1,  0};
  Int_t Marker[numSyst]={20,21, 22, 23, 20, 25, 26, 25, 25, 29, 29, 31, 32};
  Int_t MarkerBetter[numSyst]={1,22, 32, 1, 29, 1,   3,  34, 33, 1, 20, 21, 22};
  //  Int_t Color[numSyst]= {1,  2,  3,  4,  5,  6,  7,  7,  4, 10,  6,  1,  2};
  //  Int_t Color[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  //  Int_t ColorBetter[numSyst]= {1, 628,868,1, 909, 1,801,418,860,  1, 881, 7, 1};
  Int_t ColorBetter[numSyst]= {1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  // Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  // Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};
  Int_t ColorMult[nummolt+1] ={1,2,8,4,6,868};
  Int_t Color[3] ={628, 418, 600};

  TLegend *legendRegionAllF=new TLegend(0.66, 0.6, 0.78, 0.73);
  //TLegend *legendRegionAllF=new TLegend(0.66, 0.3, 0.78, 0.43);
  //TLegend *legendRegionAllF=new TLegend(0.5, 0.43, 0.9, 0.71);                                                 
  legendRegionAllF->SetFillStyle(0);
  TLegendEntry * lReAll1[3];
  TLegendEntry *lReAll2[3];

  TLegend *legendError = new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegendEntry * lReAll1Bis[3];
  TLegendEntry *lReAll2Bis[3];

  TLegend *legendStatBoxBis=new TLegend(0.76, 0.82, 0.92, 0.92);
  //TLegend *legendStatBoxBis=new TLegend(0.68, 0.62, 0.84, 0.72);

  Float_t YSup=0.01;
  Float_t YInf=-0.001;

  TH1D* fHistPhiDistr_solosist[nummolt+1][numPtV0][numregions];
  TH1D* fHistPhiDistr_solostat[nummolt+1][numPtV0][numregions];
  TH1D* fHistPhiDistr_solosistBlack[nummolt+1][numPtV0][numregions];
  TH1D* fHistPhiDistr_solostatBlack[nummolt+1][numPtV0][numregions];
  TH1D* fHistPhiDistr_solostatMarker[nummolt+1][numPtV0][numregions];

  TLine * lineat0 = new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2,0);
  lineat0->SetLineColor(1);
  lineat0->SetLineStyle(2);

  TString stringout;
  stringout = Dir+"/DATA"+year0+"/DeltaPhiProjAllSyst";
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_SysPhi%i_PtMin%.1f_", sysPhi, PtTrigMin);
  stringout += ".root";
  TFile * fileout = new TFile(stringout, "RECREATE");

  //here I start the newly  created part!!

  TString titleX=  "p_{T} (GeV/c)";
  TString titleYRel = "Relative uncertainty";
  TString title = "Multiplicity class "; 


  TF1 * lineat1 = new TF1("pol0","pol0", -1./2*TMath::Pi(), 3./2*TMath::Pi());
  lineat1->SetLineColor(kBlack);
  lineat1->SetLineWidth(0.15);
  lineat1->FixParameter(0,0);

  TCanvas* canvasPlotProj[nummolt+1];
  TCanvas* canvasApp1 = new TCanvas("canvasApp1", "canvasApp1", 800,800);

  TString PathIn0;
  TString PathIn;
  PathIn0 = Dir+"/DATA"+year0+"/PtSpectra";
  if (PtBinning>0) PathIn0 +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      PathIn0 +="_"+tipo[type];
    PathIn0 +=Srap[israp];
    PathIn0 +=SSkipAssoc[SkipAssoc];
  }
  PathIn0+=   Form("_SysPhi%i_PtMin%.1f_", sysPhi, PtTrigMin);
  TFile * filein[3];

  TString TitleNew;
  TString titleY=  "N/N_{trigg} per #Delta#eta";
  for(Int_t m=0; m<nummolt+1; m++){
    cout << " m " << m << endl;
    if (ANote){
    canvasPlotProj[m] = new TCanvas (Form("canvasPlot_m%i", m), Form("canvasPlot_m%i", m), 800, 1300);
    canvasPlotProj[m]-> Divide(3,3);
    }
    else {
      canvasPlotProj[m] = new TCanvas (Form("canvasPlot_m%i", m), Form("canvasPlot_m%i", m), 1300, 800);
    canvasPlotProj[m]-> Divide(5,2);
    }
  for (Int_t iregion=0; iregion<3; iregion++){
    PathIn = PathIn0 + RegionType[iregion];
    PathIn += ".root";
    cout<< iregion << " iregion "  << PathIn << endl;
    filein[iregion] = new TFile(PathIn, "");
    //    cout <<     "I start the loop on pT " << endl;
    for(Int_t v=PtV0Min; v < numPtV0Max; v++){
      if (type==0) TitleNew ="p^{K^{0}_{S}}_{T}  [" + SPtV0[v]+") GeV/c" ;
      else TitleNew ="p^{#Xi}_{T}  [" + SPtV0[v]+") GeV/c" ;
      //      cout << v << endl;
      fHistPhiDistr_solostat[m][v][iregion] = (TH1D*)filein[iregion]->Get("PhiDistr_solostat_m"+Smolt[m]+"_v"+SPtV0[v]);
      fHistPhiDistr_solosist[m][v][iregion] = (TH1D*)filein[iregion]->Get("PhiDistr_solosist_m"+Smolt[m]+"_v"+SPtV0[v]);

      if (!fHistPhiDistr_solosist[m][v][iregion]) return;
      if (!fHistPhiDistr_solostat[m][v][iregion]) return;
      fHistPhiDistr_solosist[m][v][iregion]->SetName("PhiDistr_solosist_m"+Smolt[m]+"_v"+SPtV0[v]+Form("_region%i",iregion));
      fHistPhiDistr_solostat[m][v][iregion]->SetName("PhiDistr_solostat_m"+Smolt[m]+"_v"+SPtV0[v]+Form("_region%i",iregion));

      Float_t AvgErr=0;
      for (Int_t dphi =1; dphi <=       fHistPhiDistr_solostat[m][v][iregion]->GetNbinsX(); dphi++){
	AvgErr +=       fHistPhiDistr_solosist[m][v][iregion]->GetBinError(dphi);
      }
      AvgErr=AvgErr/       fHistPhiDistr_solostat[m][v][iregion]->GetNbinsX();
      cout << "AvgErr " << AvgErr << endl;
      for (Int_t dphi =1; dphi <=       fHistPhiDistr_solostat[m][v][iregion]->GetNbinsX(); dphi++){
	if(fHistPhiDistr_solosist[m][v][iregion]->GetBinError(dphi) > 7*AvgErr){
	fHistPhiDistr_solosist[m][v][iregion]->SetBinError(dphi, AvgErr);
	cout << " done ! " << endl;
	}
      }

      fHistPhiDistr_solostat[m][v][iregion]->SetTitle(TitleNew);
      fHistPhiDistr_solostat[m][v][iregion]->GetYaxis()->SetTitle(titleY);
      if (!ANote)       fHistPhiDistr_solostat[m][v][iregion]->GetYaxis()->SetTitle("");
      fHistPhiDistr_solostat[m][v][iregion]->GetYaxis()->SetTitleOffset(1.6);
      fHistPhiDistr_solosist[m][v][iregion]->SetTitle(TitleNew);
      fHistPhiDistr_solosist[m][v][iregion]->GetYaxis()->SetTitle(titleY);
      if (!ANote)      fHistPhiDistr_solosist[m][v][iregion]->GetYaxis()->SetTitle("");
      fHistPhiDistr_solosist[m][v][iregion]->GetYaxis()->SetTitleOffset(1.6);

      //      cout << iregion << " " <<       fHistPhiDistr_solostat[m][v][iregion]->GetBinContent(1)<< endl;

      canvasPlotProj[m]->cd(v);
      gPad->SetLeftMargin(0.2);
      if (m==0 && type==0){
	if (NPtV0[v] < 1.2) YSup = 0.014;
	else if (NPtV0[v] < 2.5) YSup = 0.01;
	else  YSup = 0.006;
      }
      if (m==4 && type==0){
	if (NPtV0[v] < 1.2) YSup = 0.006;
	else if (NPtV0[v] < 2.5) YSup = 0.004;
	else  YSup = 0.004;
      }
      else {
	if (NPtV0[v] < 1.2) YSup = 0.014;
	else if (NPtV0[v] < 2.5) YSup = 0.0099;
	else  YSup = 0.005;
      }
      if (type==8){
	YSup=0.0005999;
	YInf = -0.000069999;

	if (m==0){
	  if (NPtV0[v]< 2) 	  YSup = 0.0012;
	  else 	  if (NPtV0[v]<4) YSup = 0.0006;
	  else 	  if (NPtV0[v]<5) YSup = 0.0006;
	}

	if (m==2 ||m==3){
	  if (NPtV0[v]< 2) 	  YSup = 0.0006;
	  else 	  if (NPtV0[v]<4) YSup = 0.0003;
	  else 	  if (NPtV0[v]<5) YSup = 0.0003;
	}
	if (m==4){
	  if (NPtV0[v]< 2) 	  YSup = 0.0003;
	  else 	  if (NPtV0[v]<4) YSup = 0.0002;
	  else 	  if (NPtV0[v]<5) YSup = 0.0002;
	}

      }


      //      TF1 * pol0 = (TF1*)       fHistPhiDistr_solostat[m][v][iregion]->GetFunction(Form("pol0ZYAM_m%i_v%i",m,v));
      //      pol0->SetLineWidth(0);

      fHistPhiDistr_solostat[m][v][iregion]->GetYaxis()->SetRangeUser(YInf, YSup);
      fHistPhiDistr_solosist[m][v][iregion]->GetYaxis()->SetRangeUser(YInf, YSup);
      //      fHistPhiDistr_solostat[m][v][iregion]->DrawClone("same hist e p");
      fHistPhiDistr_solostat[m][v][iregion]->DrawClone("same e p");
      //      fHistPhiDistr_solostat[m][v][iregion]->DrawClone("same e");
      fHistPhiDistr_solosist[m][v][iregion]->SetFillStyle(0);
      fHistPhiDistr_solosist[m][v][iregion]->DrawClone("same hist e2 p");
      //      fHistPhiDistr_solosist[m][v][iregion]->DrawClone("same  e2");
      lineat1->Draw("same");
      
      fileout->WriteTObject(fHistPhiDistr_solostat[m][v][iregion]);

      if( (type==0 && NPtV0[v] == 1.2 && m==5) || (type==8 && NPtV0[v]==2.0 && m==5)){
	canvasApp1->cd();
	canvasApp1->SetFillColor(0);
	canvasApp1->SetTickx(1);
	canvasApp1->SetTicky(1);
	//      if (type==0)      gPad->SetTopMargin(0.03);
	gPad->SetTopMargin(0.06);
	gPad->SetLeftMargin(0.2);
	if(type==8)	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.1);
	gPad->SetRightMargin(0.03);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(0);
	gStyle->SetLegendFont(42);

	fHistPhiDistr_solostatMarker[m][v][iregion]=(TH1D*)      fHistPhiDistr_solostat[m][v][iregion]->Clone("fHistPhi"+sRegion[iregion]);
	//      fHistPhiDistr_solosistMarker[m][v][iregion]=(TH1D*)      fHistPhiDistr_solosist[m][v][iregion]->Clone("fHistPhiSist"+sRegion[iregion]);
	fHistPhiDistr_solostatMarker[m][v][iregion]->SetMarkerStyle(20);
	fHistPhiDistr_solostatMarker[m][v][iregion]->SetMarkerSize(1.5);

	//      fHistPhiDistr_solosistMarker[m][v][iregion]->SetMarkerStyle(20);
	//      fHistPhiDistr_solosistMarker[m][v][iregion]->SetMarkerColor(1);

	//      lReAll2[iregion]=      legendRegionAllF->AddEntry("", sRegion1[iregion], "");
	//      lReAll2[iregion]->SetTextSize(0.04);

	//      TLegend *LegendXi=new TLegend(0.06,0.72,0.76,0.92);

	/*
	  TLegend *legendStatBoxBis=new TLegend(0.7, 0.7, 0.8, 0.8);
	  legendStatBoxBis->AddEntry(fHistPhiDistr_solosistMarker[m][v][iregion], "syst. error", "ef");
	  legendStatBoxBis->AddEntry(fHistPhiDistr_solosistMarker[m][v][iregion], "stat. error", "pe");
	*/
	//	TLegend *LegendXi=new TLegend(0.08,0.72,0.7,0.92);
	TLegend *LegendXi=new TLegend(0.03,0.72,0.65,0.92);
	//Legend1->SetTextAlign(13);                                                                                
	LegendXi->SetFillStyle(0);
	LegendXi->AddEntry("", "#bf{ALICE Preliminary}", "");
	LegendXi->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
	LegendXi->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
	LegendXi->AddEntry("", "2.0 < #it{p}_{T}^{#Xi} < 2.5 GeV/#it{c}", "");


	//      TLegend *LegendKs=new TLegend(0.06,0.75,0.76,0.95);
	TLegend *LegendKs=new TLegend(0.08,0.72,0.7,0.92);
	//Legend1->SetTextAlign(13);                                                                                
	LegendKs->SetFillStyle(0);
	LegendKs->AddEntry("", "#bf{ALICE Preliminary}", "");
	LegendKs->AddEntry("", "pp, #sqrt{#it{s}} = 13 TeV", "");
	LegendKs->AddEntry(""/*(TObject*)0*/, NameP[type]+ " correlation, #it{p}_{T}^{trigg} > 3 GeV/#it{c}", "");
	LegendKs->AddEntry("", "1.2 < #it{p}_{T}^{K^{0}_{S}} < 1.6 GeV/#it{c}", "");

	if (iregion==0){
	  fHistPhiDistr_solostatBlack[m][v][iregion]= (TH1D*) fHistPhiDistr_solostat[m][v][iregion]->Clone("fHistStatBlack");
	  fHistPhiDistr_solosistBlack[m][v][iregion]= (TH1D*) fHistPhiDistr_solosist[m][v][iregion]->Clone("fHistSistBlack");
	  fHistPhiDistr_solosistBlack[m][v][iregion]->SetMarkerColor(1);
	  fHistPhiDistr_solosistBlack[m][v][iregion]->SetLineColor(1);
	  fHistPhiDistr_solosistBlack[m][v][iregion]->SetMarkerStyle(20);

	  legendStatBoxBis->AddEntry(fHistPhiDistr_solosistBlack[m][v][iregion], "stat. error", "ep");
	  legendStatBoxBis->AddEntry(fHistPhiDistr_solosistBlack[m][v][iregion], "syst. error", "ef");
	}

	StyleHistoYield(fHistPhiDistr_solostat[m][v][iregion], YInf, YSup, Color[iregion], MarkerRegion[iregion], xTitle, titleYieldYType[type],"", MarkerSize[iregion], 0.8, yOffset[type]);
	StyleHistoYield(fHistPhiDistr_solosist[m][v][iregion], YInf, YSup, Color[iregion], MarkerRegion[iregion], xTitle, titleYieldYType[type],"", MarkerSize[iregion], 0.8, yOffset[type]);


	fHistPhiDistr_solostat[m][v][iregion]->DrawClone("same hist e0x0 p");
	fHistPhiDistr_solosist[m][v][iregion]->SetFillStyle(0);
	fHistPhiDistr_solosist[m][v][iregion]->DrawClone("same hist e2 p");
	lineat0->Draw("same");

	fHistPhiDistr_solostat[m][v][iregion]->SetMarkerSize(1.9);
	if(iregion==2)	fHistPhiDistr_solostat[m][v][iregion]->SetMarkerSize(2.5);
	lReAll1[iregion] = legendRegionAllF->AddEntry(fHistPhiDistr_solostat[m][v][iregion], " "+sRegionDEta[iregion], "p");
	lReAll1[iregion]->SetTextSize(0.03);

	if (iregion==2){
	  if (type==0)      LegendKs->Draw("");
	  else if (type==8)      LegendXi->Draw("");
	  legendStatBoxBis->Draw("");
	  legendRegionAllF->Draw("");
	}
      }
    } //end loop on pt v0
  }
  }
  cout << " going to write on file " << endl;   
  TString SANote[2] = {"_Hor", ""};
    fileout->WriteTObject(canvasApp1);
  for(Int_t m=0; m<nummolt+1; m++){
    fileout->WriteTObject(canvasPlotProj[m]);
    canvasPlotProj[m]->SaveAs("PictureForNote/DeltaPhiProjAllSyst_"+Smolt[m] + "_"+tipo[type]+SANote[ANote]+".pdf");

  }

  if (type==0)  canvasApp1->SaveAs("PictureForNote/canvasApp1_K0s.pdf");
  if (type==8)  canvasApp1->SaveAs("PictureForNote/canvasApp1_Xi.pdf");
  if (type==0)  canvasApp1->SaveAs("PictureForNote/canvasApp1_K0s.eps");
  if (type==8)  canvasApp1->SaveAs("PictureForNote/canvasApp1_Xi.eps");
  fileout->Close();
  cout << " I've produced the file " << stringout << endl;
}

