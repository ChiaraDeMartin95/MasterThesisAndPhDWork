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

void AngularCorrelationPlot(Bool_t ishhCorr=0, Int_t isProj=2, Int_t sysV0=0, Int_t sysTrigger=0,Int_t sys=0,Int_t type=0, Int_t PtIntervalShown=1,   TString year0 = "2016",TString year="2016k_onlyTriggerWithHighestPt", TString yearMC="2018f1_extra_onlyTriggerWithHighestPt",  TString Path1 ="", TString Path2 ="",TString Path3="_15runs_6thtry", TString Dir ="FinalOutput/", Int_t PtTrigChosen=3){ 

  //isJetProj==1: the fully corrected DeltaPhi projection in the in-jet region is shown; if ==0, the same for the Out-of_jet region
  if (PtIntervalShown>6){
    cout << " Thera are not so many Pt intervals " << endl;
    return;
  }
  if (ishhCorr==1){ year = "2016k"; Path1="_New"; Path2="_New"; yearMC = "2018f1_extra";}

  //masterthesis analysis è usato quando devo fare analisi presentata per la tesi
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

  TLegend * legendDataMC = new TLegend(0.6, 0.1, 0.9, 0.4);
  Int_t ColorWidth[2]={868, 418};
  TString MCOrNot[2]={" Data" , " Pythia8" };
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Ssideband[numSB]={"", "_SB"};
  TString Ssidebandtitle[numSB]={"", " (sidebands)"};
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  //  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  //  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
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
  TH1F *hDeltaEtaDeltaPhi_EtaProj[nummolt+1][numPtV0][numPtTrig];
  TString PathIn[numPtTrig];

  TF1* Gaussian[nummolt+1][numPtV0][numPtTrig];

  TLine *tlineEtaSx=new TLine(-0.5, -TMath::Pi()/2, -0.5,  3*TMath::Pi()/2);
  TLine *tlineEtaDx=new TLine(0.5, -TMath::Pi()/2, 0.5,  3*TMath::Pi()/2);
  TLine *tlinePhiSx=new TLine(-1.5, -1, 1.5,  -1);
  TLine *tlinePhiDx=new TLine(-1.5, +1, 1.5,  +1);
  TLine *tlinePhiBase=new TLine(-TMath::Pi()/2, 0, 3*TMath::Pi()/2,0);

  TCanvas *  canvasWidthGaussian = new TCanvas ("canvasWidthGaussian", "canvasWidthGaussian", 800, 500);
  canvasWidthGaussian->Divide(3,2);
  TCanvas *canvasSPectrum;
  TCanvas *canvasPlot[nummolt+1];
  TCanvas *canvasPlotProj[nummolt+1];
  TH1F *HistoWidthGaussian[nummolt+1][2][numPtTrig];
  TFile *filein[numPtTrig];
  TString nomefileoutput= "AngularCorrelationPlot" + hhCorr[ishhCorr] + "_Output.root";
  TFile *fileout = new TFile(nomefileoutput, "RECREATE");

  //disegno spettri in pTassoc per le varie molteplicità
    for(Int_t m=nummolt; m>=0; m--){


    }
  //*****************************************************

    for(Int_t m=nummolt; m>=0; m--){
      cout << "\n\n m " << m << endl;
      canvasPlot[m]=new TCanvas(Form("canvasPlot_m%i",m), "canvasPlot"+Smolt[m], 1300, 800);
      canvasPlot[m]->Divide(7,2);
      canvasPlotProj[m]=new TCanvas(Form("canvasPlotProj_m%i",m), "canvasPlotProj"+Smolt[m], 1300, 800);
      canvasPlotProj[m]->Divide(7,2);

      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){

	cout << "PtTrig " << PtTrig << endl;
	PtTrigMin=PtTrig+3;
	if (PtTrigMin!=PtTrigChosen) continue;
	HistoWidthGaussian[m][0][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%i_Data", m, PtTrigMin), "JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
	HistoWidthGaussian[m][1][PtTrig]= new TH1F(Form("HistoWidthGaussianm%i_PtTrig%i_Pythia8", m, PtTrigMin),"JetWidth vs p_{T, assoc} for mult " +Smolt[m]+ Form(" PtTrig %i", PtTrig+3), numPtV0, NPtV0);
      }

      for(Int_t v=0; v<numPtV0; v++){  
	//    if (v!=PtIntervalShown) continue;
	//	if (v>4) continue;
	nameSE[m][v]="ME_";
	nameSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[0]+"_AC";
	namePhiProjJet[m][v]= nameSE[m][v] + "_phi_V0Sub_BulkSub_EffCorr";
	namePhiProjBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_Bulk_EffCorr";
	namePhiProjJetBulk[m][v]= nameSE[m][v] + "_phi_V0Sub_JetBulkEffCorr";


      cout << nameSE[m][v]<< endl;
      for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){

	cout << "PtTrig " << PtTrig << endl;
	PtTrigMin=PtTrig+3;
	if (PtTrigMin!=PtTrigChosen) continue;
	nameEtaProj[m][v][PtTrig]= nameSE[m][v] + Form("_EtaProj_PtTrig%i", PtTrig);


	//	if (PtTrigMin==4 || PtTrigMin==5 || PtTrigMin>10)continue;
	if (PtTrigMin>7) continue;
	for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
	PathIn[PtTrig]= Dir+"/histo/AngularCorrelation" + file[IntisMC] + hhCorr[ishhCorr]+ Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f", sysTrigger, sysV0, sysang, PtTrigMin)+"_Output.root";
	cout << PathIn[PtTrig] << endl;
	filein[PtTrig]= new TFile(PathIn[PtTrig]);
	
	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]= (TH2F*)filein[PtTrig]->Get(nameSE[m][v]);
	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJet[m][v]);
	hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjBulk[m][v]);
	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]= (TH1F*)filein[PtTrig]->Get(namePhiProjJetBulk[m][v]);
	hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]= (TH1F*)	hDeltaEtaDeltaPhi_SEbins[m][v][PtTrig]->ProjectionX(nameEtaProj[m][v][PtTrig]);
	cout << "ho preso isto" << endl;
	//	canvasPlot[m]->cd(PtTrig+1);
	if (IntisMC==0)	canvasPlot[m]->cd(v+1);
	else	canvasPlot[m]->cd(v+1+7);
	cout << "scelto cd " << endl;
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
	cout << "scelto cd " << endl;

	//	canvasPlot[m]->cd(PtTrig+1+7);
	if (IntisMC==0)	canvasPlotProj[m]->cd(v+1);
	else canvasPlotProj[m]->cd(v+7+1);

	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetBinContent(	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->GetNbinsX(), 0);
	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]    ->SetLineColor(628);
	hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]   ->SetLineColor(418);
	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->SetLineColor(868);

	Gaussian[m][v][PtTrig]= new TF1 ( Form("GaussianFit_m%i_v%i_PtTrig%i", m,v,PtTrig),"gaus",  -1.2, 1.2);
	cout << "m " << m << " v " << v << " PtTrig " << PtTrig << "  " << MCOrNot[IntisMC] << endl;
	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Fit(Gaussian[m][v][PtTrig], "R");
	cout << "chi square " << Gaussian[m][v][PtTrig]->GetChisquare() << " NDF " << Gaussian[m][v][PtTrig]->GetNDF() << " chisquare/NDF " << Gaussian[m][v][PtTrig]->GetChisquare()/Gaussian[m][v][PtTrig]->GetNDF() << endl;
	hDeltaEtaDeltaPhi_PhiProjJetBulk[m][v][PtTrig]->Draw("");
	hDeltaEtaDeltaPhi_PhiProjJet[m][v][PtTrig]->Draw("same");
	hDeltaEtaDeltaPhi_PhiProjBulk[m][v][PtTrig]->Draw("same");


	tlinePhiBase->Draw("");
	hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetYaxis()->SetRangeUser(0, 1.2*	hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->GetMaximum());
	//	hDeltaEtaDeltaPhi_EtaProj[m][v][PtTrig]->Draw("");


	HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinContent(v+1, Gaussian[m][v][PtTrig]->GetParameter(2));
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetBinError(v+1, Gaussian[m][v][PtTrig]->GetParError(2));
	} //end of IntisMC loop
      } //end of PtTrigMin loop
      } //end of v loop
      canvasWidthGaussian->cd(m+1);
      for (Int_t IntisMC=0; IntisMC<=1; IntisMC++){
	for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
	  if ((PtTrig+3)!=PtTrigChosen) continue;
	HistoWidthGaussian[m][IntisMC][PtTrig]->GetXaxis()->SetTitle("p_{T, assoc}");
	HistoWidthGaussian[m][IntisMC][PtTrig]->GetYaxis()->SetTitle("#sigma (rad)");
	if (m==0)	legendDataMC->AddEntry(	HistoWidthGaussian[m][IntisMC][PtTrig],MCOrNot[IntisMC], "pl"); 
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetLineColor(ColorWidth[IntisMC]);
	HistoWidthGaussian[m][IntisMC][PtTrig]->SetMarkerColor(ColorWidth[IntisMC]);
	HistoWidthGaussian[m][IntisMC][PtTrig]->Draw("same");
	if (IntisMC==1)	legendDataMC->Draw("same");
	}
      }

      if (m==0 || m>3){
	//	fileout->WriteTObject(canvasPlot[m]); 
	fileout->WriteTObject(canvasPlotProj[m]); 
      }
    } //end of mult loop
    fileout->WriteTObject(canvasWidthGaussian);
  fileout->Close();
  cout << "\n\n ho cretao il file " << nomefileoutput << endl;
}

