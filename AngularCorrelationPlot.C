#include <TLine.h>
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

void AngularCorrelationPlot(Bool_t ishhCorr=0, Int_t PtTrigChosen=3,  Int_t isProj=2, Int_t sysV0=0, Int_t sysTrigger=0,Int_t sys=0,Int_t type=0, Int_t PtIntervalShown=1,   TString year0 = "2016",TString year="2016k", TString yearMC="2018f1_extra",  TString Path1 ="", TString Path2 ="", TString Dir ="FinalOutput/"){ 

  //isJetProj==1: the fully corrected DeltaPhi projection in the in-jet region is shown; if ==0, the same for the Out-of_jet region
  if (PtIntervalShown>6){
    cout << " Thera are not so many Pt assoc intervals " << endl;
    return;
  }
  if (ishhCorr==1){ year = "2016k"; Path1="_New"; Path2="_New"; yearMC = "2018f1_extra";}

  if (!ishhCorr && (((year!="2016k_onlyTriggerWithHighestPt"&& year!="2018f1_extra_onlyTriggerWithHighestPt") || yearMC!="2018f1_extra_onlyTriggerWithHighestPt"))) {
    cout << "for hV0 correlation you should use year = 2016k_onlyTriggerWithHighestPt together with the MC yearMC = 2018f1_extra_onlyTriggerWithHighestPt" << endl;
    return; 
  }

  Float_t PtTrigMin;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  //lista degli effetti  sistematici studiati in questa macro
  if (sys==3 || sys>5) return;
  if (sysV0>6)return;
  if (ishhCorr && sys<4 && sys!=0) return; //se faccio correlazione hh non posso studiare sistematico associato a scelta regione Sidebands e peak nella distribuzione di massa invariante
  if (sysV0>2 && ishhCorr) return;
  if (sysTrigger!=0) return;
  if (sysV0!=0 && sys!=0) return;

  TString hhCorr[2]= {"", "_hhCorr"};
  Int_t sysang=0;
  if (sys==1) sysang=1;
  if (sys==2) sysang=2;

  Dir+="DATA"+year0;
  TString file[2];
  file[0] = year+Path1;
  file[1] = yearMC + "_MCEff" +Path2;
 
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrig=10;
  const Int_t numtipo=2;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
 
  TLegend * legendDataMC = new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t ColorWidth[2]={868, 418};
  TString MCOrNot[2]={" Data" , " Pythia8" };
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};
  TString SPtV0[numPtV0]={"0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,1.5,2,2.5,3,4,8};
  TString SNPtV0[numPtV0+1]={"0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0"};

  TString nameSE[nummolt+1][numPtV0];
  TString namePhiProjJet[nummolt+1][numPtV0];
  TString namePhiProjBulk[nummolt+1][numPtV0];
  TString namePhiProjJetBulk[nummolt+1][numPtV0];
  TString nameEtaProj[nummolt+1][numPtV0][numPtTrig];
  TH2F *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk[nummolt+1][numPtV0][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[nummolt+1][2][numPtTrig];
  TH1F *hDeltaEtaDeltaPhi_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TString PathIn[numPtTrig];

  TF1* Gaussian[nummolt+1][numPtV0][numPtTrig];
  TF1* GaussianAS[nummolt+1][numPtV0][numPtTrig];

  TLine *tlineEtaSx=new TLine(-0.5, -TMath::Pi()/2, -0.5,  3*TMath::Pi()/2);
  TLine *tlineEtaDx=new TLine(0.5, -TMath::Pi()/2, 0.5,  3*TMath::Pi()/2);
  TLine *tlinePhiSx=new TLine(-1.5, -1, 1.5,  -1);
  TLine *tlinePhiDx=new TLine(-1.5, +1, 1.5,  +1);
  TLine *tlinePhiBase=new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2,0);

  TCanvas *  canvasWidthGaussian[2];
  TCanvas *  canvasWidthGaussianAS[2];
  for (Int_t i =0; i<2; i++){  
    canvasWidthGaussian[i] = new TCanvas (Form("canvasWidthGaussian%i",i), Form("canvasWidthGaussian%i",i), 800, 500);
    canvasWidthGaussian[i]->Divide(3,2);
    canvasWidthGaussianAS[i] = new TCanvas (Form("canvasWidthGaussianAS%i",i), Form("canvasWidthGaussianAS%i",i), 800, 500);
    canvasWidthGaussianAS[i]->Divide(3,2);
  }
  TCanvas *canvasSPectrum;
  TCanvas *canvasPlot[nummolt+1];
  TCanvas *canvasPlotProj[nummolt+1];
  TH1F *HistoWidthGaussian[nummolt+1][2][numPtTrig];
  TH1F *HistoWidthGaussianAS[nummolt+1][2][numPtTrig];
  TH1F *RatioHistoWidthGaussian[nummolt+1][numPtTrig];
  TH1F *RatioHistoWidthGaussianAS[nummolt+1][numPtTrig];


  //*********************************************************************                                                                                                  
  //**************calcolo numero particelle di trigger*******************                                                                                                
  Int_t NTrigger[nummolt+1][2];
  TFile *fileinbis[2]; //Data and MC
  TString PathInBis[2];
  PathInBis[0] =  "FinalOutput/AnalysisResults" + year + Path1  + ".root";
  PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC  + "_MCEff" + Path1 +".root";
  if (ishhCorr){  
  PathInBis[0] =  "FinalOutput/AnalysisResults" + year +Path1+ "_hhCorr.root";
  PathInBis[1] =  "FinalOutput/AnalysisResults" + yearMC  + Path1+"_hhCorr_MCEff.root";
  }

  for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
  fileinbis[IntisMC]=new TFile(PathInBis[IntisMC],"");
  if (!fileinbis[IntisMC]) return;
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis[IntisMC]->Get("MyTask");
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
  if (!fHistTriggervsMult ) return;
  cout << "ok up to here " << endl;
  TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigChosen+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(30-0.00001) );

  //*********************************************************************                                                                                                
 
  for(Int_t m=0; m<nummolt+1; m++){
    NTrigger[m][IntisMC]=0;
    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
        NTrigger[m][IntisMC]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
        NTrigger[m][IntisMC]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    cout << "n trigger in mult range (all triggers)    " << m << "  " <<  " is MC " << IntisMC << "  " << NTrigger[m][IntisMC] <<   endl;
  }
  }
  ///////////////////////////////////////////////////////////////////////////////////////////              
 
 TFile *filein[numPtTrig];
  TString nomefileoutput= "AngularCorrelationPlot" + hhCorr[ishhCorr] + Form("_PtTrigMin%i_Output.root", PtTrigChosen);
  TFile *fileout = new TFile(nomefileoutput, "RECREATE");

  //disegno spettri in pTassoc per le varie molteplicità
  for(Int_t m=nummolt; m>=0; m--){


  }
  //*****************************************************
   TCanvas*   canvasSummedPt=new TCanvas ("canvasSummedPt", "canvasSummedPt", 800, 500);
   canvasSummedPt->Divide(nummolt+1, 2);

  for(Int_t m=nummolt; m>=0; m--){
    cout << "\n\n m " << m << endl;
    canvasPlot[m]=new TCanvas(Form("canvasPlot_m%i",m), "canvasPlot"+Smolt[m], 1300, 800);
    canvasPlot[m]->Divide(7,2);
    canvasPlotProj[m]=new TCanvas(Form("canvasPlotProj_m%i",m), "canvasPlotProj"+Smolt[m], 1300, 800);
    canvasPlotProj[m]->Divide(7,2);

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      PtTrigMin=PtTrig+3;
      if (PtTrigMin!=PtTrigChosen) continue;
      cout << "PtTrig " << PtTrig << endl;

      HistoWidthGaussian[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussian[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussianAS[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianASm%i_PtTrig%.0f_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      HistoWidthGaussianAS[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianASm%i_PtTrig%.0f_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);

    }


    for(Int_t v=0; v<numPtV0; v++){  
      //    if (v!=PtIntervalShown) continue;
      //	if (v>4) continue;
      nameSE[m][v]="ME_";
      nameSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_AC";
      namePhiProjJet[m][v]= nameSE[m][v] + "_phi_V0Sub_BulkSub_EffCorr";
      namePhiProjBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_Bulk_EffCorr";
      namePhiProjJetBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_JetBulkEffCorr";
      //      cout << nameSE[m][v]<< endl;

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	PtTrigMin=PtTrig+3;
	if (PtTrigMin!=PtTrigChosen) continue;
	cout << "PtTrig " << PtTrig << endl;
	nameEtaProj[m][v][PtTrig]= nameSE[m][v] + Form("_EtaProj_PtTrig%i", PtTrig);
	//	if (PtTrigMin==4 || PtTrigMin==5 || PtTrigMin>10)continue;
	if (PtTrigMin>7) continue;
	for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
	  PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC] + hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_Output.root";
	  cout << PathIn[PtTrig] << endl;
	  filein[PtTrig]= new TFile(PathIn[PtTrig], "");
	
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSE[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]);
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig]);
	  //	  cout << "ho preso isto" << endl;
	  //	canvasPlot[m]->cd(PtTrig+1);
	  if (IntisMC==0)	canvasPlot[m]->cd(v+1);
	  else	canvasPlot[m]->cd(v+1+7);
	  //	  cout << "scelto cd " << endl;
	  //	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetLogz();
	  //	gPad->SetLogz();
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->SetTitle("AC v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f ", PtTrigMin) + MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->SetTitle("Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin) +MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->SetTitle("Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin)+ MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetTitle("Phi Proj v"+SPtV0[v]+ Form(" for p_{T}^{Trigg, min} %.1f", PtTrigMin)+ MCOrNot[IntisMC]);
	  hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->Draw("colz");
	  tlinePhiDx->Draw();
	  tlinePhiSx->Draw();
	  tlineEtaDx->Draw();
	  tlineEtaSx->Draw();

	  //	canvasPlot[m]->cd(PtTrig+1+7);
	  if (IntisMC==0)	canvasPlotProj[m]->cd(v+1);
	  else canvasPlotProj[m]->cd(v+7+1);

	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(628);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(418);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetLineColor(868);

	  //division by number of trigger particles
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Scale(1./NTrigger[m][IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Scale(1./NTrigger[m][IntisMC]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Scale(1./NTrigger[m][IntisMC]);

	  //sum of different ptV0 bins
	  if (v==0){
	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    =	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjJet[m][0][PtTrig]    ->Clone(Form("PhiProjJet_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   =	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjBulk[m][0][PtTrig]   ->Clone(Form("PhiProjBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]=	  (TH1F*) hDeltaEtaDeltaPhi_PhiProjJetBulk[m][0][PtTrig]->Clone(Form("PhiProjJetBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
    
	  }

	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->SetTitle(Form("PhiProjJet_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->SetTitle(Form("PhiProjBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));	  
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->SetTitle(Form("PhiProjJetBulk_PtTrigMin%.1f_m%i_isMC%i", PtTrigMin, m,IntisMC));

	  if (v!=0){
	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]);
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]    ->Add(	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]);
	  }

	  Gaussian[m][v][PtTrig]= new TF1 ( Form("GaussianFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  -0.8, 0.8);
	  GaussianAS[m][v][PtTrig]= new TF1 ( Form("GaussianASFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  TMath::Pi()-1.2,   TMath::Pi()+1.2);
	  GaussianAS[m][v][PtTrig]->SetLineColor(419);

	  cout << "m " << m << " v " << v << " PtTrig " << PtTrig << "  " << MCOrNot[IntisMC] << endl;
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Fit(Gaussian[m][v][PtTrig], "R");
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Fit(GaussianAS[m][v][PtTrig], "R");

	  cout << "\njet fit: chi square " << Gaussian[m][v][PtTrig]->GetChisquare() << " NDF " << Gaussian[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << Gaussian[m][v][PtTrig]->GetChisquare()/Gaussian[m][v][PtTrig]->GetNDF() << endl;
	  cout << "\nAwaySide fit: chi square " << GaussianAS[m][v][PtTrig]->GetChisquare() << " NDF " << GaussianAS[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << GaussianAS[m][v][PtTrig]->GetChisquare()/GaussianAS[m][v][PtTrig]->GetNDF() << endl;

	  //	  if (	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetMaximum() > 	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetMaximum()) 	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Draw("");
//	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetRangeUser(-0.004, 0.05);
//	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetRangeUser(-0.004, 0.05);
//	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetRangeUser(-0.004, 0.05);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->GetYaxis()->SetLabelSize(0.06);
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetYaxis()->SetLabelSize(0.06);
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->GetYaxis()->SetLabelSize(0.06);
	  gPad->SetLeftMargin(0.15);
	  hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->Draw("same");



	  tlinePhiBase->Draw("");
	  hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2*	hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
	  //	hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Draw("");

	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinContent(v+1, Gaussian[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinError(v+1, Gaussian[m][v][PtTrig]->GetParError(2));
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinContent(v+1, GaussianAS[m][v][PtTrig]->GetParameter(2));
	  HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetBinError(v+1, GaussianAS[m][v][PtTrig]->GetParError(2));

	} //end of IntisMC loop
      } //end of PtTrigMin loop
    } //end of v loop


    //I draw delta-phi projection done summing over pT assoc bins*****
    for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
      if (!IntisMC)      canvasSummedPt->cd(1+m);
      else       canvasSummedPt->cd(1+m+nummolt+1);
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+3)!=PtTrigChosen) continue;
    for(Int_t v=0; v<numPtV0; v++){  
      cout <<" m " << m <<  " v " << v << " "  <<  	  hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->GetMaximum() << endl;
    }
    cout <<" m " << m <<  " " <<	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetMaximum() << endl;
            hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Scale(2.28);
    	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->Scale(2.28);
    	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->Scale(2.28);
    cout <<" m " << m <<  " " <<	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetMaximum() << endl;
	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetYaxis()->SetRangeUser(-0.004, 0.5);
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->GetYaxis()->SetRangeUser(-0.004, 0.5);
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(-0.004, 0.5);

	    hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->GetYaxis()->SetLabelSize(0.05);
	    hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->GetYaxis()->SetLabelSize(0.05);
	    hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_PhiProjJet_PtV0Summed[m][IntisMC][PtTrig]    ->Draw("");
	  hDeltaEtaDeltaPhi_PhiProjBulk_PtV0Summed[m][IntisMC][PtTrig]   ->Draw("same");
	  hDeltaEtaDeltaPhi_PhiProjJetBulk_PtV0Summed[m][IntisMC][PtTrig]->Draw("same");
      }
    }

    //**************************************************************
    for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	if ((PtTrig+3)!=PtTrigChosen) continue;

	HistoWidthGaussian[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0,0.7);
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);
	if (m==0)	legendDataMC->AddEntry(	HistoWidthGaussian[m][IntisMC][PtTrig],MCOrNot[IntisMC], "pl"); 

	HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->GetYaxis()->SetRangeUser(0,5);
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);

	canvasWidthGaussian[0]->cd(m+1);
	HistoWidthGaussian[m][IntisMC][PtTrig]->Draw("same");
	if (IntisMC==1)	legendDataMC->Draw("same");

	canvasWidthGaussianAS[0]->cd(m+1);
	HistoWidthGaussianAS[m][IntisMC][PtTrig]->Draw("same");
	if (IntisMC==1)	legendDataMC->Draw("same");

      } 
    } //end second intisMC loop

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      PtTrigMin=PtTrig+3;
      if ((PtTrig+3)!=PtTrigChosen) continue;
      RatioHistoWidthGaussian[m][PtTrig]=(TH1F*)	HistoWidthGaussian[m][0][PtTrig]->Clone(Form("RatioDataPythiaWidthm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussian[m][PtTrig]->SetTitle(Form("RatioDataPythiaWidthm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussian[m][PtTrig]->Sumw2();
      RatioHistoWidthGaussian[m][PtTrig]->Divide(HistoWidthGaussian[m][1][PtTrig]);
      RatioHistoWidthGaussian[m][PtTrig]->GetYaxis()->SetRangeUser(0,1.2);

      RatioHistoWidthGaussianAS[m][PtTrig]=(TH1F*)	HistoWidthGaussianAS[m][0][PtTrig]->Clone(Form("RatioDataPythiaWidthASm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussianAS[m][PtTrig]->SetTitle(Form("RatioDataPythiaWidthASm%i_PtTrig%.0f", m, PtTrigMin));
      RatioHistoWidthGaussianAS[m][PtTrig]->Sumw2();
      RatioHistoWidthGaussianAS[m][PtTrig]->Divide(HistoWidthGaussianAS[m][1][PtTrig]);
      RatioHistoWidthGaussianAS[m][PtTrig]->GetYaxis()->SetRangeUser(0,2);

      canvasWidthGaussian[1]->cd(m+1);
      RatioHistoWidthGaussian[m][PtTrig]->Draw("same");

      canvasWidthGaussianAS[1]->cd(m+1);
      RatioHistoWidthGaussianAS[m][PtTrig]->Draw("same");	  
    }

    //    if (m==0 || m>3){
      //	fileout->WriteTObject(canvasPlot[m]); 
      fileout->WriteTObject(canvasPlotProj[m]); 

      //    }

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      if ((PtTrig+3)!=PtTrigChosen) continue;
      /*
      fileout->WriteTObject(RatioHistoWidthGaussian[m][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussian[m][0][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussian[m][1][PtTrig]);
      fileout->WriteTObject(RatioHistoWidthGaussianAS[m][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussianAS[m][0][PtTrig]);
      fileout->WriteTObject(HistoWidthGaussianAS[m][1][PtTrig]);
      */
    }
  } //end of mult loop
  fileout->WriteTObject(canvasSummedPt); 
  fileout->WriteTObject(canvasWidthGaussian[0]);
  fileout->WriteTObject(canvasWidthGaussianAS[0]);
  fileout->WriteTObject(canvasWidthGaussian[1]);
  fileout->WriteTObject(canvasWidthGaussianAS[1]);
  fileout->Close();
  cout << "\n\n ho cretao il file " << nomefileoutput << endl;
}

