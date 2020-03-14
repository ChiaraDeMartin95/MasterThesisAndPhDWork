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

void AngularCorrelationPlotBis(Int_t isDataOrMC=0, Bool_t isAS=0){ 
  cout << "isDataOrMC =0 for data, =1 for MC = 2 for Data/MC " << endl;

  Int_t PtTrigMin;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TString hhCorr[2]= {"", "_hhCorr"};
  TString SDataOrMC[3]= {"Data", "Pythia8", "Ratio"};
  TString SisAS[2]={"", "AS"};
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrig=5;
  const Int_t numtipo=2;
  const Int_t numeta=3;
  const Int_t numetabis=3;
  const Int_t numSB=2;
 
  TLegend * legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  Int_t ColorPt[5]={1, 401, 801, 628, 909};
  Int_t Marker[2]={33, 27};
  TString MCOrNot[2]={" Data" , " Pythia8" };
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString SPtV0[numPtV0]={"0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,1.5,2,2.5,3,4,8};
  TString SNPtV0[numPtV0+1]={"0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0"};

  TH1F* HistoGaussianWidth[nummolt+1][numPtTrig][2];
  TH1F* RatioHistoGaussianWidth[nummolt+1][numPtTrig][2];
  TFile *filein[numPtTrig][2];

  TCanvas *canvasWidthCfr[numPtTrig];
  TCanvas *  canvasCfrPtTrigMinhK0s=new TCanvas("canvasCfrPtTrigMinhK0s", "canvasCfrPtTrigMinhK0s", 800, 500);
  TCanvas *  canvasCfrPtTrigMinhh=new TCanvas("canvasCfrPtTrigMinhh", "canvasCfrPtTrigMinhh", 800, 500);
  canvasCfrPtTrigMinhh->Divide(3,2);
  canvasCfrPtTrigMinhK0s->Divide(3,2);

  TString nomefileoutput = "AngularCorrelationPlotBisOutput" + SDataOrMC[isDataOrMC]+SisAS[isAS]+ ".root";
  TFile *fileout = new TFile(nomefileoutput, "RECREATE");

    for(Int_t PtTrig=0; PtTrig<numPtTrig; PtTrig++){
      PtTrigMin=PtTrig+3;
      cout << " PtTrigMin " << PtTrigMin << endl;
      canvasWidthCfr[PtTrig]	=new TCanvas (Form("canvasWidthCfrPtTrigMin%i", PtTrigMin) ,Form("canvasWidthCfrPtTrigMin%i", PtTrigMin), 800, 500);
      canvasWidthCfr[PtTrig]->Divide(3,2);
  for (Int_t ishhCorr=0; ishhCorr<=1; ishhCorr++){
    if (ishhCorr==0)    cout << "\n\nI'm analyzing hK0s" <<endl; 
    else    cout << "\n\nI'm analyzing hh" <<endl; 

      if (ishhCorr==0 && PtTrigMin>5) continue;
      if (ishhCorr==1 && PtTrigMin>7) continue;
     TString nomefilein= "AngularCorrelationPlot" + hhCorr[ishhCorr] + Form("_PtTrigMin%i_Output.root", PtTrigMin);
      filein[PtTrig][ishhCorr]= new TFile(nomefilein, "");
      if (!filein[PtTrig][ishhCorr]) return;
      for(Int_t m=nummolt; m>=0; m--){
	cout << "   m " << m << endl;
	TString nameHisto;
	if (isDataOrMC<=1) nameHisto = 	"HistoWidthGaussian" + SisAS[isAS]+Form("m%i_PtTrig%i_", m, PtTrigMin)+SDataOrMC[isDataOrMC];
	else nameHisto = "RatioDataPythiaWidth"+ SisAS[isAS]+ Form("m%i_PtTrig%i",  m, PtTrigMin);
	HistoGaussianWidth[m][PtTrig][ishhCorr]= (TH1F*)filein[PtTrig][ishhCorr]->Get(nameHisto);
	if (!HistoGaussianWidth[m][PtTrig][ishhCorr]) {cout << "histogram " <<nameHisto<< "not found " << endl; return;}

	HistoGaussianWidth[m][PtTrig][ishhCorr]->SetLineColor(ColorPt[PtTrig]);
	HistoGaussianWidth[m][PtTrig][ishhCorr]->SetMarkerColor(ColorPt[PtTrig]);
	HistoGaussianWidth[m][PtTrig][ishhCorr]->SetMarkerStyle(Marker[ishhCorr]);

	if (m==nummolt)	legend->AddEntry(HistoGaussianWidth[m][PtTrig][ishhCorr], Form("PtTrigMin%i_", PtTrigMin)+hhCorr[ishhCorr], "pl"); 

	canvasWidthCfr[PtTrig]->cd(m+1);
	HistoGaussianWidth[m][PtTrig][ishhCorr]->Draw("same");
	if (PtTrig==0 && ishhCorr==0)	legend->Draw("same");

	if (ishhCorr==0){
	canvasCfrPtTrigMinhK0s->cd(m+1);
	HistoGaussianWidth[m][PtTrig][ishhCorr]->Draw("same");
	if (PtTrig==0 && ishhCorr==0)	legend->Draw("same");
	}
	else{
	canvasCfrPtTrigMinhh->cd(m+1);
	HistoGaussianWidth[m][PtTrig][ishhCorr]->Draw("same");
	if (PtTrig==0 && ishhCorr==0)	legend->Draw("same");

	}
      } //end of mult loop
  }//end of hh hK0s loop
  fileout->WriteTObject(canvasWidthCfr[PtTrig]);
    }//end of PtTrigLoop

  fileout->WriteTObject(canvasCfrPtTrigMinhK0s);
  fileout->WriteTObject(canvasCfrPtTrigMinhh);
  fileout->Close();
  cout << "\n\n ho cretao il file " << nomefileoutput << endl;
}

