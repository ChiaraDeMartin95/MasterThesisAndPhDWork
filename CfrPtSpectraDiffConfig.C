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
#include <TSpline.h>
#include <TLine.h>
#include <Macros/ErrRatioCorr.C>

TSpline3 *sp3;//= (TSpline3*) sp ->Clone("splineFioClone");
Double_t spline(Double_t *x, Double_t* p) {
  Double_t xx = x[0];
  return sp3->Eval(xx);
}


void CfrPtSpectraDiffConfig(Bool_t CfrDiffPeriods=0, TString year0="2016", Int_t ishhCorr=0, Int_t NumberOfRegions=3,  Int_t type=0, Bool_t AllRun2=0, Bool_t CfrSkipAssoc=0, TString Path1=""/*"_FioPtBins"*/, Bool_t ispp5TeV=1, Bool_t isNormFactor=0){// Int_t HistoType=0){

  Bool_t Comp13or5TeV;
  if (ispp5TeV){
    cout << "Do you want to compare with 13 TeV results (with appropriate mult classes) (option = 0) or with 5 TeV results (with corresponding mult classes) (option = 1) ? " << endl;
    cin >> Comp13or5TeV; 
  }

  if (CfrSkipAssoc) {AllRun2=1; CfrDiffPeriods=0;}
  if (CfrDiffPeriods) AllRun2=1;
  TString yearMC;
  TString yeardati;
  if (AllRun2==0){
    yearMC="2015g3b1";
    yeardati="2015f";
    Path1="";
  }
  else {
    yearMC="2016kl";
    yeardati="Run2DataRed_hXi";
  }

  if (type==0){
    if (ispp5TeV){    
      yearMC = "17pq_hK0s_pttrig0.15";
      yeardati ="17pq_hK0s_pttrig0.15";
    }
    else {
      yearMC="2018f1_extra_hK0s";
      yeardati="2016k_hK0s";
    }
  }

  Int_t NumberOfRegionsFinal = NumberOfRegions;
  if (AllRun2==0)  NumberOfRegionsFinal = NumberOfRegions-2;
  if (Path1 == "_FioPtBins")   NumberOfRegionsFinal = 1; //NumberOfRegions-3;

  TString yeardati15f = "2015f";
  TString yearMC15g3b1 = "2015g3b1";

  gStyle->SetOptStat(0);

  TString PathIn;
  TFile *filein;
  TFile *fileout;
  const Int_t numhistoType=3;

  TString NameHisto[numhistoType]={"S", "Efficiency", "SEffCorr"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
  const Int_t numtipo=10;
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString nomefileoutput = "CfrPtSpectraDiffConfig_"+yeardati;
  if (CfrDiffPeriods) nomefileoutput += "_CfrDiffPeriods";
  if (CfrSkipAssoc) nomefileoutput += "_CfrSkipAssoc";
  nomefileoutput += Path1 + "_Plot";
  if (ispp5TeV)  nomefileoutput += "_pp5TeV";
  if (isNormFactor) nomefileoutput += "_IsNormFactor";
  TString  nomefileoutputPDF=nomefileoutput;
  nomefileoutput += "_Try";
  nomefileoutput += ".root";
  fileout = new TFile (nomefileoutput, "RECREATE");

  TString innamedati="FinalOutput/DATA" +year0+"/invmass_distribution_thesis/invmass_distribution";
  TString innamecompleto[20];
  TString innameMC = innamedati;
  innameMC+="_MCEff";
  if (type==0) innamedati += "_PtBinning1";
  if (type==0) innameMC += "_PtBinning1";
  innameMC+=Path1;
  innamedati+=Path1;

  TString innamedataorMC[2]={innamedati, innameMC};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;
  //colori e marker diversi per jet, OJ, J+OJ *********************************
  TString Smolt[nummolt+1] = {"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  Int_t Marker[6]={20,21,20,21, 20, 21};
  Int_t Color[6]={628,797, 433,881,991,1};
  Int_t MarkerBis[2]={23,24};
  Int_t MarkerTris[2]={20, 33};
  Int_t MarkerMBStat[6]={20,4,20,4, 20, 4};
  Int_t MarkerMBSist[6]={21,25,21,25, 21, 25};

  //cfr spettri pubblicati per le Xi
  cout << "\nPrendo histo per confronto con dati pubblicati\n" << endl;
  TString PathDatiPubblicati ="";
  if (type==0) {
    PathDatiPubblicati = "HEPData-ins1748157-v1-Table_1.root";
  }
  else if (type==8) PathDatiPubblicati = "HEPData-1583750454-v1-Table_3.root";
  TFile *filedatipubblicati = new TFile(PathDatiPubblicati, "");
  if (!filedatipubblicati) {cout << "file dati pubblicati not there " << endl; return;}
  TString STable = "Table 3";
  if (type==0) STable = "Table 1";
  TDirectoryFile *dirspectra = (TDirectoryFile*)filedatipubblicati->Get(STable);
  if (!dirspectra)  {cout << "directory dati pubblicati not there " << endl; return;}

  TString  PathDatiPubblicati5TeV = "CorrectedSpectra_K0S.root";
  TString  PathDatiPubblicati5TeVMB = "CorrectedSpectra_K0S_IntMult.root";
  TFile *filedatipubblicati5TeV = new TFile(PathDatiPubblicati5TeV, "");
  if (!filedatipubblicati5TeV) {cout << "file dati pubblicati 5 TeV not there " << endl; return;}
  TFile *filedatipubblicati5TeVMB = new TFile(PathDatiPubblicati5TeVMB, "");
  if (!filedatipubblicati5TeVMB) {cout << "file dati pubblicati 5 TeV MB not there " << endl; return;}

  TH1F* hspectrum[12];
  TH1F* hspectrum1[12];
  TH1F* hspectrum2[12];
  TH1F* hspectrum3[12];
  TH1F* hspectrumetot[12];
  TH1F* hspectrumCfr[6];
  TH1F* hspectrumCfrSys[6];
  TH1F* hspectrumCfrSum;
  TH1F * hspectrumCfrFioComparisonMB;
  TSpline3 *splineFio[6];
  TF1*   fsplineFio[6];
  TF1*   fsplineFioSum;

  for (Int_t i=0; i<12; i++){
    if (!ispp5TeV){
      if (i==11) continue;
      hspectrum[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i", i+1));
      hspectrum1[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e1", i+1)); //stat error
      hspectrum2[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e2", i+1)); //total sist error
      hspectrum3[i] = (TH1F*)dirspectra->Get(Form("Hist1D_y%i_e3", i+1)); //sist uncorr error
    }
    else {
      if (i==11) {
	hspectrum[i] = (TH1F*)filedatipubblicati5TeVMB->Get("CorrectedSpectra_K0S");
	if (!hspectrum[i]) {cout << "histogram not there " << endl; return;}
	hspectrum1[i] = (TH1F*)hspectrum[i]->Clone("CorrectedSpectra_K0S[%i]_StatErr");
	hspectrum2[i] = (TH1F*)filedatipubblicati5TeVMB->Get("CorrectedSpectra_Syst_K0S");
	hspectrum3[i] =  (TH1F*)filedatipubblicati5TeVMB->Get("CorrectedSpectra_Syst_Uncor_K0S");
      }
      else {
	hspectrum[i] = (TH1F*)filedatipubblicati5TeV->Get(Form("Corrected_RawSpectra_K0S[%i]", i));
	if (!hspectrum[i]) {cout << "histogram not there " << endl; return;}
	hspectrum1[i] = (TH1F*)hspectrum[i]->Clone(Form("Corrected_RawSpectra_K0S[%i]_StatErr", i));
	hspectrum2[i] = (TH1F*)filedatipubblicati5TeV->Get(Form("Corrected_RawSpectra_Syst_K0S[%i]", i));
	hspectrum3[i] =  (TH1F*)filedatipubblicati5TeV->Get(Form("Corrected_RawSpectra_Syst_UncK0S[%i]", i));
      }
    }

    if (!hspectrum[i] || !hspectrum1[i] || !hspectrum2[i]|| !hspectrum3[i] ) { cout << "histo is missing " << endl; return;}
    if(ispp5TeV){
      for (Int_t l=1; l <= hspectrum[i]->GetNbinsX() ; l++){
	hspectrum1[i]->SetBinContent(l, hspectrum[i]->GetBinError(l));
	hspectrum2[i]->SetBinContent(l, hspectrum2[i]->GetBinError(l));
	hspectrum3[i]->SetBinContent(l, hspectrum3[i]->GetBinError(l));
      }
    }
    hspectrumetot[i]= (TH1F*)    hspectrum3[i]->Clone(Form("Hist1D_y%i_etot", i+1));

    if (i<=nummolt) {
      hspectrumCfr[i] = (TH1F*)hspectrum[i]->Clone(Form("Spectrum_multStat%i",i+1));
      hspectrumCfrSys[i] = (TH1F*)hspectrum[i]->Clone(Form("Spectrum_multSys%i",i+1));
    }
    //    cout << " I clone a histo " << endl;                                                                                                     
    for (Int_t b=1; b<= hspectrumCfr[0]->GetNbinsX();b++){
      //      hspectrumetot[i]->SetBinContent(b,   sqrt(pow(hspectrum1[i]->GetBinContent(b),2) + pow(hspectrum2[i]->GetBinContent(b),2) + pow(hspectrum3[i]->GetBinContent(b),2)) );
      hspectrumetot[i]->SetBinContent(b,   sqrt(pow(hspectrum1[i]->GetBinContent(b),2) + pow(hspectrum2[i]->GetBinContent(b),2)) );
    }
  }

  hspectrumCfrSum= (TH1F*)hspectrum[0]->Clone("Spectrum_multSum");
  for (Int_t b=1; b<= hspectrumCfr[0]->GetNbinsX();b++){

    if (ispp5TeV){
      //Silvia's mult classes: [0-0.01],[0.01-0.1],[0.1-1],[1-5],[5-10],[10-20],[20-30],[30-40],[40-50],[50-70],[70-100]
      hspectrumCfr[0]->SetBinContent(b,1./5 * (0.01*hspectrum[0]->GetBinContent(b)+ 0.09*hspectrum[1]->GetBinContent(b) + 0.9*hspectrum[2]->GetBinContent(b) + 4* hspectrum[3]->GetBinContent(b)));
      hspectrumCfrSys[0]->SetBinContent(b,1./5 * (0.01*hspectrum[0]->GetBinContent(b)+ 0.09*hspectrum[1]->GetBinContent(b) + 0.9*hspectrum[2]->GetBinContent(b) + 4* hspectrum[3]->GetBinContent(b)));

      hspectrumCfr[0]->SetBinError(b, sqrt(pow(0.01/5*hspectrum1[0]->GetBinContent(b), 2)+ pow(0.09/5*hspectrum1[1]->GetBinContent(b),2) + pow(0.9/5*hspectrum1[2]->GetBinContent(b),2) + pow(4/5* hspectrum1[3]->GetBinContent(b),2)));
      hspectrumCfrSys[0]->SetBinError(b, sqrt(pow(0.01/5*hspectrum2[0]->GetBinContent(b), 2)+ pow(0.09/5*hspectrum2[1]->GetBinContent(b),2) + pow(0.9/5*hspectrum2[2]->GetBinContent(b),2) + pow(4/5* hspectrum2[3]->GetBinContent(b),2)));

      hspectrumCfr[1]->SetBinContent(b,hspectrum[4]->GetBinContent(b));
      hspectrumCfr[1]->SetBinError(b, hspectrum1[4]->GetBinContent(b));
      hspectrumCfrSys[1]->SetBinContent(b,hspectrum[4]->GetBinContent(b));
      hspectrumCfrSys[1]->SetBinError(b, hspectrum2[4]->GetBinContent(b));

      hspectrumCfr[2]->SetBinContent(b,1./20 * (hspectrum[5]->GetBinContent(b)*10+ hspectrum[6]->GetBinContent(b)*10));
      hspectrumCfr[2]->SetBinError(b, sqrt(pow(hspectrum1[5]->GetBinContent(b)*10./20,2) + pow(hspectrum1[6]->GetBinContent(b)*10./20,2)));
      hspectrumCfrSys[2]->SetBinContent(b,1./20 * (hspectrum[5]->GetBinContent(b)*10+ hspectrum[6]->GetBinContent(b)*10));
      hspectrumCfrSys[2]->SetBinError(b, sqrt(pow(hspectrum2[5]->GetBinContent(b)*10./20,2) + pow(hspectrum2[6]->GetBinContent(b)*10./20,2)));

      hspectrumCfr[3]->SetBinContent(b,1./20 * (hspectrum[7]->GetBinContent(b)*10+ hspectrum[8]->GetBinContent(b)*10));
      hspectrumCfr[3]->SetBinError(b, sqrt(pow(hspectrum1[7]->GetBinContent(b) * 10./20,2) + pow(hspectrum1[8]->GetBinContent(b)*10./20,2)));
      hspectrumCfrSys[3]->SetBinContent(b,1./20 * (hspectrum[7]->GetBinContent(b)*10+ hspectrum[8]->GetBinContent(b)*10));
      hspectrumCfrSys[3]->SetBinError(b, sqrt(pow(hspectrum2[7]->GetBinContent(b) * 10./20,2) + pow(hspectrum2[8]->GetBinContent(b)*10./20,2)));

      hspectrumCfr[4]->SetBinContent(b,1./50 * (hspectrum[9]->GetBinContent(b)*20+ hspectrum[10]->GetBinContent(b)*30));
      hspectrumCfr[4]->SetBinError(b, sqrt(pow(hspectrum1[9]->GetBinContent(b) * 20./50,2) + pow(hspectrum1[10]->GetBinContent(b)*30./50,2)));
      hspectrumCfrSys[4]->SetBinContent(b,1./50 * (hspectrum[9]->GetBinContent(b)*20+ hspectrum[10]->GetBinContent(b)*30));
      hspectrumCfrSys[4]->SetBinError(b, sqrt(pow(hspectrum2[9]->GetBinContent(b) * 20./50,2) + pow(hspectrum2[10]->GetBinContent(b)*30./50,2)));

      hspectrumCfr[5]->SetBinContent(b,hspectrum[11]->GetBinContent(b));
      hspectrumCfr[5]->SetBinError(b, hspectrum1[11]->GetBinContent(b));
      hspectrumCfrSys[5]->SetBinContent(b,hspectrum[11]->GetBinContent(b));
      hspectrumCfrSys[5]->SetBinError(b, hspectrum2[11]->GetBinContent(b));

      //OPTION A: 0-100% spectrum as average of published mult classes
      hspectrumCfrSum->SetBinContent(b, 0.0001*hspectrum[0]->GetBinContent(b)+ 0.0009*hspectrum[1]->GetBinContent(b)+0.009*hspectrum[2]->GetBinContent(b)+ 0.04*hspectrum[3]->GetBinContent(b)+ 0.05*hspectrum[4]->GetBinContent(b) +0.1 *hspectrum[5]->GetBinContent(b) +0.1*hspectrum[6]->GetBinContent(b) +0.1*hspectrum[7]->GetBinContent(b) +0.1*hspectrum[8]->GetBinContent(b) +0.2*hspectrum[9]->GetBinContent(b)+0.3*hspectrum[10]->GetBinContent(b));
      hspectrumCfrSum->SetBinError(b, sqrt(pow(0.0001*hspectrum[0]->GetBinError(b),2)+ pow(0.0009*hspectrum[1]->GetBinError(b),2)+pow(0.009*hspectrum[2]->GetBinError(b),2)+ pow(0.04*hspectrum[3]->GetBinError(b),2)+ pow(0.05*hspectrum[4]->GetBinError(b),2) +pow(0.1 *hspectrum[5]->GetBinError(b),2) +pow(0.1*hspectrum[6]->GetBinError(b),2) +pow(0.1*hspectrum[7]->GetBinError(b),2) +pow(0.1*hspectrum[8]->GetBinError(b),2) +pow(0.2*hspectrum[9]->GetBinError(b),2) + pow(0.3*hspectrum[10]->GetBinError(b),2)));

      /* OPTION B: 0-100% spectrum as average of average of published mult classes
      //can be taken from 13 TeV and modified
      */

    }
    else {
      hspectrumCfr[0]->SetBinContent(b,1./5 * (hspectrum[0]->GetBinContent(b)+ hspectrum[1]->GetBinContent(b)*4));
      hspectrumCfrSys[0]->SetBinContent(b,1./5 * (hspectrum[0]->GetBinContent(b)+ hspectrum[1]->GetBinContent(b)*4));

      hspectrumCfr[0]->SetBinError(b, sqrt(pow(hspectrum1[0]->GetBinContent(b) * 1./5,2) + pow(hspectrum1[1]->GetBinContent(b)*4./5,2)));
      hspectrumCfrSys[0]->SetBinError(b, sqrt(pow(hspectrum2[0]->GetBinContent(b) * 1./5,2) + pow(hspectrum2[1]->GetBinContent(b)*4./5,2)));

      hspectrumCfr[1]->SetBinContent(b,hspectrum[2]->GetBinContent(b));
      hspectrumCfr[1]->SetBinError(b, hspectrum1[2]->GetBinContent(b));
      hspectrumCfrSys[1]->SetBinContent(b,hspectrum[2]->GetBinContent(b));
      hspectrumCfrSys[1]->SetBinError(b, hspectrum2[2]->GetBinContent(b));

      hspectrumCfr[2]->SetBinContent(b,1./20 * (hspectrum[3]->GetBinContent(b)*5+ hspectrum[4]->GetBinContent(b)*5 + hspectrum[5]->GetBinContent(b)*10));
      hspectrumCfr[2]->SetBinError(b, sqrt(pow(hspectrum1[3]->GetBinContent(b)*5./20,2) + pow(hspectrum1[4]->GetBinContent(b)*5./20,2) + pow(hspectrum1[5]->GetBinContent(b)*10./20,2)));
      hspectrumCfrSys[2]->SetBinContent(b,1./20 * (hspectrum[3]->GetBinContent(b)*5+ hspectrum[4]->GetBinContent(b)*5 + hspectrum[5]->GetBinContent(b)*10));
      hspectrumCfrSys[2]->SetBinError(b, sqrt(pow(hspectrum2[3]->GetBinContent(b)*5./20,2) + pow(hspectrum2[4]->GetBinContent(b)*5./20,2) + pow(hspectrum2[5]->GetBinContent(b)*10./20,2)));

      hspectrumCfr[3]->SetBinContent(b,1./20 * (hspectrum[6]->GetBinContent(b)*10+ hspectrum[7]->GetBinContent(b)*10));
      hspectrumCfr[3]->SetBinError(b, sqrt(pow(hspectrum1[6]->GetBinContent(b) * 10./20,2) + pow(hspectrum1[7]->GetBinContent(b)*10./20,2)));
      hspectrumCfrSys[3]->SetBinContent(b,1./20 * (hspectrum[6]->GetBinContent(b)*10+ hspectrum[7]->GetBinContent(b)*10));
      hspectrumCfrSys[3]->SetBinError(b, sqrt(pow(hspectrum2[6]->GetBinContent(b) * 10./20,2) + pow(hspectrum2[7]->GetBinContent(b)*10./20,2)));

      hspectrumCfr[4]->SetBinContent(b,1./50 * (hspectrum[8]->GetBinContent(b)*20+ hspectrum[9]->GetBinContent(b)*30));
      hspectrumCfr[4]->SetBinError(b, sqrt(pow(hspectrum1[8]->GetBinContent(b) * 20./50,2) + pow(hspectrum1[9]->GetBinContent(b)*30./50,2)));
      hspectrumCfrSys[4]->SetBinContent(b,1./50 * (hspectrum[8]->GetBinContent(b)*20+ hspectrum[9]->GetBinContent(b)*30));
      hspectrumCfrSys[4]->SetBinError(b, sqrt(pow(hspectrum2[8]->GetBinContent(b) * 20./50,2) + pow(hspectrum2[9]->GetBinContent(b)*30./50,2)));

      hspectrumCfr[5]->SetBinContent(b,hspectrum[10]->GetBinContent(b));
      hspectrumCfr[5]->SetBinError(b, hspectrum1[10]->GetBinContent(b));
      hspectrumCfrSys[5]->SetBinContent(b,hspectrum[10]->GetBinContent(b));
      hspectrumCfrSys[5]->SetBinError(b, hspectrum2[10]->GetBinContent(b));

      //OPTION A: 0-100% spectrum as average of published mult classes
      hspectrumCfrSum->SetBinContent(b, 0.009*hspectrum[0]->GetBinContent(b)+ 0.036*hspectrum[1]->GetBinContent(b)+0.044*hspectrum[2]->GetBinContent(b)+ 0.046*hspectrum[3]->GetBinContent(b)+ 0.045*hspectrum[4]->GetBinContent(b) +0.09 *hspectrum[5]->GetBinContent(b) +0.091*hspectrum[6]->GetBinContent(b) +0.092*hspectrum[7]->GetBinContent(b) +0.192*hspectrum[8]->GetBinContent(b) +0.355*hspectrum[9]->GetBinContent(b));
      hspectrumCfrSum->SetBinError(b, sqrt(pow(0.009*hspectrum[0]->GetBinError(b),2)+ pow(0.036*hspectrum[1]->GetBinError(b),2)+pow(0.044*hspectrum[2]->GetBinError(b),2)+ pow(0.046*hspectrum[3]->GetBinError(b),2)+ pow(0.045*hspectrum[4]->GetBinError(b),2) +pow(0.09 *hspectrum[5]->GetBinError(b),2) +pow(0.091*hspectrum[6]->GetBinError(b),2) +pow(0.092*hspectrum[7]->GetBinError(b),2) +pow(0.192*hspectrum[8]->GetBinError(b),2) +pow(0.355*hspectrum[9]->GetBinError(b),2)));

      /* OPTION B: 0-100% spectrum as average of average of published mult classes
      //    hspectrumCfrSum->SetBinContent(b, 0.05*hspectrumCfr[0]->GetBinContent(b)+ 0.05*hspectrumCfr[1]->GetBinContent(b)+0.2*hspectrumCfr[2]->GetBinContent(b)+ 0.2*hspectrumCfr[3]->GetBinContent(b)+ 0.5*hspectrumCfr[4]->GetBinContent(b));
      //hspectrumCfrSum->SetBinError(b, sqrt(pow(0.05*hspectrumCfr[0]->GetBinError(b),2)+ pow(0.05*hspectrumCfr[1]->GetBinError(b),2)+pow(0.2*hspectrumCfr[2]->GetBinError(b),2)+ pow(0.2*hspectrumCfr[3]->GetBinError(b),2)+ pow(0.5*hspectrumCfr[4]->GetBinError(b),2)));
      */
    }
  }

  if (ispp5TeV){ //I need to compare differetn multiplicity class to have similar dNdeta
    if (Comp13or5TeV ==0){
      //5-10% 13 TeV ~ 0-5% 5 TeV
      hspectrumCfr[0] = hspectrumCfr[1];
      hspectrumCfrSys[0] = hspectrumCfrSys[1];
      //10-30% 13 TeV ~ 5-10% 5 TeV
      hspectrumCfr[1] = hspectrumCfr[2];
      hspectrumCfrSys[1] = hspectrumCfrSys[2];
    }
  }

  for (Int_t i=0; i<6; i++){
    if (type==4 || type==5)    hspectrumCfr[i]->Scale(1./2.);// this was to consider only Xi+ or Xi-
    splineFio[i] = new TSpline3(hspectrumCfr[i],Form("splineFio_%i",i));
  }

  cout << " \nSpline of spectra obtained from Fiorella's ones and evaluation of their mean relative error " << endl;
  Float_t  MeanRelErr[6]={0};
  Float_t  SigmaRelErr[6]={0};
  for (Int_t i=0; i<6; i++){
    //    splineFio[i] = new TSpline3(hspectrumCfr[i],Form("splineFio_%i",i));
    MeanRelErr[i]=0;
    SigmaRelErr[i]=0;

    for(Int_t b=1; b <= hspectrumCfr[i]->GetNbinsX(); b++){
      MeanRelErr[i]+=sqrt(pow(hspectrumCfr[i]->GetBinError(b),2) + pow(hspectrumCfrSys[i]->GetBinError(b),2))/hspectrumCfr[i]->GetBinContent(b);
    }
    MeanRelErr[i]= MeanRelErr[i]/hspectrumCfr[i]->GetNbinsX();
    for(Int_t b=1; b <= hspectrumCfr[i]->GetNbinsX(); b++){
      SigmaRelErr[i]+=pow(MeanRelErr[i]-(hspectrumCfr[i]->GetBinError(b)/hspectrumCfr[i]->GetBinContent(b)),2);
    }
    SigmaRelErr[i]=sqrt( SigmaRelErr[i]/(hspectrumCfr[i]->GetNbinsX()-1));
  }

  TString Region[4] = {"AllAssoc0.15", "SkipAssoc0.15", "AllAssoc3", "SkipAssoc3"};
  TString RegionCfrDiffPeriods[4] = {"AllAssocRun2", "SkipAssocRun2", "AllAssocLHC15", "SkipAssocLHC15"};
  if (CfrDiffPeriods){
    for (Int_t i=0; i<4; i++){
      Region[i] = RegionCfrDiffPeriods[i] ;
    }
  }
  TString isDataorMC[2]={"Data", "MC"};
  TString StringSystV0Analysis[2]={"", "_SystV0Analysis"};
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};
  TH1F * fHistSpectrum_master[NumberOfRegions][nummolt+1];
  TH1F * fHistSpectrum_master_scaled[NumberOfRegions][nummolt+1];
  TH1F * fHistSpectrum_SumMult[NumberOfRegions];
  TH1F * fHistSpectrum_SumMultApprox[NumberOfRegions];
  TH1F * fHistSpectrum_ratio[NumberOfRegions][nummolt+1];
  TH1F * histo_NTrigger[nummolt+1];
  Int_t NTrigger[nummolt+1];
  Float_t EffTrigger[nummolt+1];
  Float_t NScale[nummolt+1]={0.045, 0.044, 0.181, 0.183, 0.547, 1};

  TCanvas * canvasSpectrum[2][numhistoType];
  for (Int_t i=0; i<2; i++){
    for (Int_t HistoType=0; HistoType<numhistoType; HistoType++){
      canvasSpectrum[i][HistoType] = new TCanvas(Form("canvasSpectrum_%i_", i)+NameHisto[HistoType], Form("canvasSpectrum_%i_", i)+NameHisto[HistoType], 1300, 800);
      canvasSpectrum[i][HistoType]->Divide(3,2);
      //      gPad->SetLeftMargin(0.0015);
    }
  }

  //******Getting the Normalisation factor if present******
  TH1F*      fHistNormFactor[nummolt+1];
  Float_t Denom[nummolt+1] = {0};
  Float_t ErrDenom[nummolt+1] = {0};
  Float_t NewContent[nummolt+1] = {0};
  TH1F *     fHistEventLossMultAll;
  TH1F *     fHistEventLoss;
  TString  PathNormFactor = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrueLHC16kl_pass2_GP_Fio_Hybrid_vs_LHC16kl_pass2_GP_Fio_PtBinning1_K0s_y0.5_PtMin0.2.root";
  if (isNormFactor){
    TFile * fileNormFactor = new TFile(PathNormFactor, "");
    if (!fileNormFactor) {cout << PathNormFactor << " does not exist " << endl; return;}
    fHistEventLoss = (TH1F*) fileNormFactor->Get("fHistEventLoss");
    fHistEventLossMultAll = (TH1F*) fileNormFactor->Get("fHistEventLossMultAll");
    if (!fHistEventLoss || !fHistEventLossMultAll) {cout << "no histos trigger efficiency " << endl;return;}
    for (Int_t m=0; m< nummolt+1; m++){
      if (m<nummolt)      EffTrigger[m] =  fHistEventLoss->GetBinContent(m+1);
      fHistNormFactor[m] = (TH1F*) fileNormFactor->Get(Form("SpectrumRatioAll_m%i", m));
    }
    EffTrigger[nummolt] =  fHistEventLossMultAll->GetBinContent(1);
    for (Int_t m=0; m< nummolt+1; m++){
      //      cout <<"Trigger efficiency: " <<   EffTrigger[m]<< " (mult = " << m << endl;
    }
  }


  //**************** prendo gli istogrammi **********************************************
  Int_t israp=0;
  Bool_t SkipAssoc=0;
  Float_t PtTrigMin=0;
  TLegend *basiclegend = new TLegend(0.6,0.7, 0.9, 0.9);
  TLegend *legendError = new TLegend(0.6,0.55, 0.9, 0.7);
  TLegend *  legendRegion=new TLegend(0.6,0.7, 0.9, 0.9);
  TLegend *  legendRatio=new TLegend(0.6,0.7, 0.9, 0.9);
  TLegend *  legendRatioCompFio=new TLegend(0.6,0.7, 0.9, 0.9);
  TLine * lineOne;
  //  if (type==8) lineOne = new TLine (0,1,6.5,1);
  //  else  if (type==0) lineOne = new TLine (0,1,8,1);
  lineOne = new TLine (0,1,8,1);
  if (Path1=="_FioPtBins")   lineOne = new TLine (0,1,6.5,1);
  TF1 * rettaOne = new TF1 ("pl0", "pol0", 0,8);
  rettaOne->FixParameter(0,1);
  Float_t Max[numhistoType] = {0.003, 0.30, 0.03};
  //  cout << " Starting the loop " << endl;
  for (Int_t HistoType=0; HistoType<numhistoType; HistoType++){  
    //    cout << "HistoType " << HistoType<< endl;
    for(Int_t l=0; l<NumberOfRegionsFinal; l++){
      if (ispp5TeV) {israp =1; SkipAssoc=0; PtTrigMin =0.15;}
      if (l==0) {israp=1; SkipAssoc=0; PtTrigMin=0.15;}
      else     if (l==1) {israp=1; SkipAssoc=1; PtTrigMin=0.15;}
      else     if (l==2) {israp=1; SkipAssoc=0; if (CfrDiffPeriods)PtTrigMin=0.15; else PtTrigMin=3; }
      else     if (l==3) {israp=1; SkipAssoc=1; if (CfrDiffPeriods)PtTrigMin=0.15; else PtTrigMin=3; }
      if (CfrDiffPeriods){
	if (l==0 || l==1){    innamedataorMC[0]=innamedati+"_"+yeardati;    innamedataorMC[1]=innameMC+"_"+yearMC;}
	if (l==2 || l==3){     innamedataorMC[0]=innamedati+"_"+yeardati15f;    innamedataorMC[1]= innameMC+"_"+yearMC15g3b1;}
      }
      else{
	if (l==0 && HistoType==0)	innamedati+="_"+yeardati;    innameMC+="_"+yearMC;
      }

      //************************** Drawing the pt spectra
      for (Int_t m=0; m< nummolt+1; m++){
	//	cout << " m " << m << endl;
	Int_t sysTrigger=0;
	Int_t sysV0=0;
	Int_t sys=0;
	Int_t dataorMC =0;
	Int_t isMeanFixedPDG=1;
	Int_t isBkgParab=0;
	if (!CfrDiffPeriods){
	  innamedataorMC[0]=innamedati;
	  innamedataorMC[1]= innameMC;
	}
	PathIn= Form(innamedataorMC[dataorMC]+"_" +tipo[type] +Srap[israp]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, sys, PtTrigMin);
	if (HistoType>0)	cout << "\nNome file input (m = " << m << "):\n " <<PathIn << endl;
	filein=new TFile(PathIn, "");   
	fHistSpectrum_master[l][m] = (TH1F*) filein->Get("histo_"+NameHisto[HistoType]);
	if (!fHistSpectrum_master[l][m]) return;

	if (m==nummolt){
	  fHistSpectrum_SumMultApprox[l] = (TH1F*)fHistSpectrum_master[l][nummolt-1]->Clone("histo_"+NameHisto[HistoType]+ "_SumMultApprox");
	  for (Int_t b=1; b<= fHistSpectrum_SumMultApprox[l]->GetNbinsX(); b++){
	    fHistSpectrum_SumMultApprox[l]->SetBinContent(b, 0.05*fHistSpectrum_master[l][0]->GetBinContent(b)+ 0.05*fHistSpectrum_master[l][1]->GetBinContent(b)+0.2*fHistSpectrum_master[l][2]->GetBinContent(b)+ 0.2*fHistSpectrum_master[l][3]->GetBinContent(b)+ 0.5*fHistSpectrum_master[l][4]->GetBinContent(b));
	    fHistSpectrum_SumMultApprox[l]->SetBinError(b, sqrt(pow(0.05*fHistSpectrum_master[l][0]->GetBinError(b),2)+ pow(0.05*fHistSpectrum_master[l][1]->GetBinError(b),2)+pow(0.2*fHistSpectrum_master[l][2]->GetBinError(b),2)+ pow(0.2*fHistSpectrum_master[l][3]->GetBinError(b),2)+ pow(0.5*fHistSpectrum_master[l][4]->GetBinError(b),2)));
	  }
	  fHistSpectrum_SumMultApprox[l]->SetLineColor(881);
	  fHistSpectrum_SumMultApprox[l]->SetMarkerColor(881);
	}
	histo_NTrigger[m] = (TH1F*)filein->Get("histo_NTrigger");
	if ( histo_NTrigger[m]){
	  NTrigger[m] = histo_NTrigger[m] ->GetBinContent(1);
	  if (m==nummolt && HistoType!=1) {
	    fHistSpectrum_SumMult[l] = (TH1F*)fHistSpectrum_master[l][nummolt-1]->Clone("histo_"+NameHisto[HistoType]+ "_SumMult");
	    if (isNormFactor) fHistSpectrum_SumMult[l]->Scale(1./EffTrigger[nummolt-1]);
	    fHistSpectrum_SumMult[l]->Scale(NTrigger[nummolt-1]);
	    //fHistSpectrum_SumMult[l]->Scale(NScale[nummolt-1]);//alternative to scale by EffTrigger and NTrigger
	    for (Int_t m=0; m<nummolt-1; m++){
	      fHistSpectrum_master_scaled[l][m] = (TH1F*)fHistSpectrum_master[l][m]->Clone("histo_"+NameHisto[HistoType]+ "_scale");
	      if (isNormFactor) fHistSpectrum_master_scaled[l][m]->Scale(1./EffTrigger[m]);
	      fHistSpectrum_master_scaled[l][m]->Scale(NTrigger[m]);
	      //  fHistSpectrum_master_scaled[l][m]->Scale(NScale[m]); //alternative to scale by EffTrigger and NTrigger
	      fHistSpectrum_SumMult[l]->Add(fHistSpectrum_master_scaled[l][m]);
	    }
	    fHistSpectrum_SumMult[l]->Scale(1./NTrigger[nummolt] );
	    if (isNormFactor){
	      fHistSpectrum_SumMult[l]->Scale(EffTrigger[nummolt] );
	      for (Int_t m=0; m<nummolt; m++){
		//cout << "Sigma/SigmaTOT m:" << m << " " << (float)NTrigger[m]/NTrigger[nummolt] * EffTrigger[nummolt]/EffTrigger[m] << endl;
	      }
	    }
	    fHistSpectrum_SumMult[l]->SetLineColor(867);
	    fHistSpectrum_SumMult[l]->SetMarkerColor(867);
	    fHistSpectrum_SumMult[l]->SetMarkerStyle(33);
	  }
	  //	  cout << "Total trigg "<< NTrigger[nummolt] << endl;
	  Int_t Tmp=0;
	  for (Int_t m=0; m<nummolt; m++){
	    //cout << "Ntrigger m " << m << " " << NTrigger[m]<< endl;
	    Tmp+=NTrigger[m];
	  }
	  //	  cout << "Tmp " << Tmp << endl;
	}
	fHistSpectrum_master[l][m]->SetLineColor(Color[l]);
	fHistSpectrum_master[l][m]->SetMarkerColor(Color[l]);
	fHistSpectrum_master[l][m]->SetMarkerStyle(20);
	if (m==0 && HistoType==0) legendRegion->AddEntry(  fHistSpectrum_master[l][m], Region[l], "pl");
	if (type==4 || type==5 || type==8)   {
	  fHistSpectrum_master[l][m]->GetYaxis()->SetRangeUser(0,Max[HistoType]);
	}
	if (type==8 && NameHisto[HistoType]!="Efficiency"){
	  fHistSpectrum_master[l][m]->GetYaxis()->SetRangeUser(0,2*Max[HistoType]);
	}
	if ( NameHisto[HistoType]=="Efficiency") {
	  fHistSpectrum_master[l][m]->GetYaxis()->SetTitleOffset(1);
	  if (type==8)	  fHistSpectrum_master[l][m]->GetYaxis()->SetRangeUser(0.00001, 0.3);
	  if (type==0)	  fHistSpectrum_master[l][m]->GetYaxis()->SetRangeUser(0.00001, 0.5);
	}
	fHistSpectrum_master[l][m]->GetXaxis()->SetTitleSize(0.05);
	if ( NameHisto[HistoType]=="SEffCorr") {
	  //	  fHistSpectrum_master[l][m]->GetYaxis()->SetTitleOffset(1);
	  fHistSpectrum_master[l][m]->GetYaxis()->SetTitle("1/N_{Evt} dN/dp_{T}"); 
	  fHistSpectrum_master[l][m]->GetYaxis()->SetTitleOffset(1.2);
	  fHistSpectrum_master[l][m]->GetYaxis()->SetTitleSize(0.05);
	  if (m==0 && l==0) basiclegend->AddEntry(	  fHistSpectrum_master[l][m] ,"This analysis", "pl");
	  if (isNormFactor) {
	    for (Int_t b=1; b<=fHistSpectrum_master[l][m]->GetNbinsX() ; b++){
	      Denom[m] = fHistNormFactor[m]->GetBinContent(fHistNormFactor[m]->FindBin(fHistSpectrum_master[l][m]->GetBinCenter(b)));
	      ErrDenom[m] = fHistNormFactor[m]->GetBinError(fHistNormFactor[m]->FindBin(fHistSpectrum_master[l][m]->GetBinCenter(b)));
	      NewContent[m] =  fHistSpectrum_master[l][m]->GetBinContent(b)/Denom[m];
	      fHistSpectrum_master[l][m]->SetBinError(b, sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b)/fHistSpectrum_master[l][m]->GetBinContent(b), 2) + pow (ErrDenom[m]/Denom[m], 2))* NewContent[m]);
	      fHistSpectrum_master[l][m]->SetBinContent(b, NewContent[m]);
	    }
	    //	    fHistSpectrum_master[l][m]->Divide(fHistNormFactor[m]);
	    //	    cout << fHistSpectrum_master[l][m]->GetNbinsX() << " " << fHistNormFactor[m]->GetNbinsX() << endl;
	  }
	}
	canvasSpectrum[0][HistoType]->cd(m+1);
	if (HistoType==2) gPad->SetLogy();
	gPad->SetLeftMargin(0.15);
	//	gStyle->SetOptStat(0);
	fHistSpectrum_master[l][m]->SetTitle("Multiplicity class " + Smolt[m]);
	if (type==0 && HistoType==2)	fHistSpectrum_master[l][m]->GetYaxis()->SetRangeUser(10e-5,4);
	if (type==8 && HistoType==2)	fHistSpectrum_master[l][m]->GetYaxis()->SetRangeUser(10e-6,0.5);
	//	if (!((m==nummolt && HistoType==2)||(m==nummolt && HistoType==2)) )	fHistSpectrum_master[l][m]->Draw("ep");
	//if (!(m==nummolt && HistoType!=0))	fHistSpectrum_master[l][m]->Draw("ep");
	fHistSpectrum_master[l][m]->Draw("ep");
	if (m==nummolt && HistoType!=1 &&  histo_NTrigger[m]) {
	  //	  cout << "hey, I'm drawing your histo!"<< endl;
	  for (Int_t b=1; b< fHistSpectrum_SumMult[l]->GetNbinsX(); b++){
	    //	    cout <<	    fHistSpectrum_SumMult[l]->GetBinContent(b)<< " " << fHistSpectrum_master[l][m]->GetBinContent(b)<<endl;
	  }
	  //	  fHistSpectrum_SumMult[l]->Draw("same ep");
	  TSpline3* splineFioSum = new TSpline3(hspectrumCfrSum,"splineFio0100Weghted");
	  sp3= (TSpline3*) splineFioSum ->Clone("splineFioSum");
	  fsplineFioSum= new TF1("fSplineFioSum", spline, 0, 8); 
	  hspectrumCfrSum->SetLineColor(kPink+6);
	  hspectrumCfrSum->SetMarkerColor(kPink+6);
	  hspectrumCfrSum->Draw("same ep");
	}
	hspectrumCfr[m]->Sumw2();
	if (l==0)      hspectrumCfr[m]->Scale(1);
	hspectrumCfr[m]->GetYaxis()->SetRangeUser(10e-6,0.5);
	if (type==0) 	hspectrumCfr[m]->GetYaxis()->SetRangeUser(10e-5,4); //1.4 if not log
	hspectrumCfr[m]->GetXaxis()->SetRangeUser(0,8);
	hspectrumCfr[m]->SetMarkerStyle(33);
	hspectrumCfrSys[m]->GetYaxis()->SetRangeUser(10e-6,0.5);
	if (type==0) 	hspectrumCfrSys[m]->GetYaxis()->SetRangeUser(10e-5,4);
	hspectrumCfrSys[m]->GetXaxis()->SetRangeUser(0,8);
	hspectrumCfrSys[m]->SetMarkerStyle(33);

	sp3= (TSpline3*) splineFio[m] ->Clone(Form("splineFioClone_%i",m));
	fsplineFio[m]= new TF1(Form("fSplineFio_%i",m), spline, 0, 8); 
	fsplineFio[m]->SetLineColor(kBlue);

	if (!CfrDiffPeriods && NameHisto[HistoType]=="SEffCorr"){
	  if (m==0 && l==0)	  {
	    basiclegend->AddEntry(	  hspectrumCfr[m],"Published", "pl");
	    legendError->AddEntry(	  hspectrumCfrSys[m],"syst.", "ef");
	    legendError->AddEntry(	  hspectrumCfr[m],"stat.", "pel");
	  }
	  //	  if (m!=nummolt ){//|| ( m==nummolt && type==0)){
	  if (kTRUE){
	    hspectrumCfr[m]->Draw("same ep");
	    hspectrumCfrSys[m]->SetFillStyle(0);
	    hspectrumCfrSys[m]->Draw("same p e2");
	    basiclegend->Draw("");
	    legendError->Draw("");
	  }
	  //	  if (type==0 && m!=nummolt) {
	  //	  if (type==0) {
	  if (Path1!="_FioPtBins"){
	    splineFio[m]->Draw("same");
	    fsplineFio[m]->DrawClone("same");
	  }
	  //	  }
	}

	canvasSpectrum[1][HistoType]->cd(m+1);
	gPad->SetLeftMargin(0.15);
	//	gStyle->SetOptStat(0);
	fHistSpectrum_ratio[l][m] = (TH1F*)       fHistSpectrum_master[l][m] ->Clone(NameHisto[HistoType]+"_Ratio");

	if ( NameHisto[HistoType]=="SEffCorr") {
	  //	  fHistSpectrum_master[l][m]->GetYaxis()->SetTitleOffset(1);
	  fHistSpectrum_ratio[l][m]->GetYaxis()->SetTitle("Ratio to published spectrum"); 
	  fHistSpectrum_ratio[l][m]->GetYaxis()->SetTitleOffset(1.2);
	}

	if (!CfrDiffPeriods && NameHisto[HistoType]=="SEffCorr" && !CfrSkipAssoc){
	  if (m==0) legendRatioCompFio->AddEntry(  fHistSpectrum_ratio[l][m],"Ratio " + Region[l] + "/Fiorella's" , "pl");
	  if (Path1 == "_FioPtBins" ) {
	    for(Int_t b=2; b <=  fHistSpectrum_ratio[l][m]->GetNbinsX(); b++){
	      // if (m<nummolt){
	      
	      if (kTRUE){
		//cmparison between my spectra and average of Fiorella's ones
		fHistSpectrum_ratio[l][m]->SetBinContent(b,  fHistSpectrum_master[l][m]->GetBinContent(b)/hspectrumCfr[m]->GetBinContent(b-1));
		fHistSpectrum_ratio[l][m]->SetBinError(b,  sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2)+pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*hspectrumCfr[m]->GetBinError(b-1),2))/hspectrumCfr[m]->GetBinContent(b-1));
	      }
	      
	      /*
		else { //comparison between Fiorella's spectrum in 0-100% and the weighted average of Fio's spectra in the different mult classes 
		fHistSpectrum_ratio[l][m]->SetBinContent(b,  hspectrumCfr[m]->GetBinContent(b-1)/hspectrumCfrSum->GetBinContent(b-1));
		fHistSpectrum_ratio[l][m]->SetBinError(b,  sqrt(pow(hspectrumCfr[m]->GetBinError(b-1),2)+pow(fHistSpectrum_ratio[l][m]->GetBinContent(b-1)*hspectrumCfrSum->GetBinError(b-1),2))/hspectrumCfrSum->GetBinContent(b-1));
		}
	      */
	      /*
		else { //comparison between my spectrum in 0-100% and the weighted average of Fio's spectra in the different mult classes  
		fHistSpectrum_ratio[l][m]->SetBinContent(b,  fHistSpectrum_master[l][m]->GetBinContent(b)/hspectrumCfrSum->GetBinContent(b-1));
		fHistSpectrum_ratio[l][m]->SetBinError(b,  sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2)+pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*hspectrumCfrSum->GetBinError(b-1),2))/hspectrumCfrSum->GetBinContent(b-1));
		}
	      */
	      /*	      
			      else { //comparison between my spectrum in 0-100% and the weighted average of my spectra in the different mult classes
			      fHistSpectrum_ratio[l][m]->SetBinContent(b,  fHistSpectrum_master[l][m]->GetBinContent(b)/fHistSpectrum_SumMult[l]->GetBinContent(b));
			      fHistSpectrum_ratio[l][m]->SetBinError(b,  sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2)+pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*fHistSpectrum_SumMult[l]->GetBinError(b),2))/fHistSpectrum_SumMult[l]->GetBinContent(b));
			      }
	      */
	    }
	  }
	  else {
	    for(Int_t b=1; b <=  fHistSpectrum_ratio[l][m]->GetNbinsX(); b++){
	      //	      fHistSpectrum_ratio[l][m]->SetBinContent(b, fHistSpectrum_master[l][m]->GetBinContent(b)/splineFio[m]->Eval(fHistSpectrum_master[l][m]->GetBinCenter(b)));

	      fHistSpectrum_ratio[l][m]->SetBinContent(b, fHistSpectrum_master[l][m]->GetBinContent(b)/fsplineFio[m]->Integral(fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinUpEdge(b))*fHistSpectrum_ratio[l][m]->GetBinWidth(b));

	      //	      fHistSpectrum_ratio[l][m]->SetBinError(b, sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2) +pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*MeanRelErr[m]* fHistSpectrum_master[l][m]->GetBinContent(b),2)) /splineFio[m]->Eval(fHistSpectrum_master[l][m]->GetBinCenter(b))); //error calculated assuming zero correlation between my spectra and Fiorella's ; maybe wrong the right one should be:

	      //	      fHistSpectrum_ratio[l][m]->SetBinError(b, sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2) +pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*MeanRelErr[m]*splineFio[m]->Eval(fHistSpectrum_master[l][m]->GetBinCenter(b))  ,2)) /splineFio[m]->Eval(fHistSpectrum_master[l][m]->GetBinCenter(b))); //error calculated assuming zero correlation between my spectra and Fiorella's ; 
	      fHistSpectrum_ratio[l][m]->SetBinError(b, sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2) +pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*MeanRelErr[m]*splineFio[m]->Eval(fHistSpectrum_master[l][m]->GetBinCenter(b))  ,2)) / (fsplineFio[m]->Integral(fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinUpEdge(b))/fHistSpectrum_ratio[l][m]->GetBinWidth(b)  )); //error calculated assuming zero correlation between my spectra and Fiorella's ; 
	      //	fHistSpectrum_ratio[l][m]->SetBinError(b, TMath::Abs(fHistSpectrum_master[l][m]->GetBinError(b) -fHistSpectrum_ratio[l][m]->GetBinContent(b)*MeanRelErr[m]* fHistSpectrum_master[l][m]->GetBinContent(b)) /splineFio[m]->Eval(fHistSpectrum_master[l][m]->GetBinCenter(b))); //error calculated assuming full correlation between my spectra and Fiorella's ;  

	      if (m==nummolt){
		//TRY		fHistSpectrum_ratio[l][m]->SetBinContent(b, fHistSpectrum_SumMult[l]->GetBinContent(b)/fsplineFio[m]->Integral(fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinUpEdge(b))*fHistSpectrum_ratio[l][m]->GetBinWidth(b));
		//TRY2		fHistSpectrum_ratio[l][m]->SetBinContent(b, fHistSpectrum_SumMult[l]->GetBinContent(b)/fsplineFioSum->Integral(fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinLowEdge(b), fHistSpectrum_ratio[l][m]->GetXaxis()->GetBinUpEdge(b))*fHistSpectrum_ratio[l][m]->GetBinWidth(b));
	      }
	      //comparison between my spectrum in 0-100% and the weighted average of my spectra in the different mult classes
	      if (m==nummolt){
		//		fHistSpectrum_ratio[l][m]->SetBinContent(b,  fHistSpectrum_master[l][m]->GetBinContent(b)/fHistSpectrum_SumMult[l]->GetBinContent(b));                                                                      
		//		fHistSpectrum_ratio[l][m]->SetBinError(b,  sqrt(pow(fHistSpectrum_master[l][m]->GetBinError(b),2)+pow(fHistSpectrum_ratio[l][m]->GetBinContent(b)*fHistSpectrum_SumMult[l]->GetBinError(b),2))/fHistSpectrum_SumMult[l]->GetBinContent(b));     
	      }
	     
	    }
	  }
	}
	else if (CfrSkipAssoc){
	  if (l==1 || l==3 ){
	    if (m==0 && HistoType==0) legendRatio->AddEntry(  fHistSpectrum_ratio[l][m],"Ratio " + Region[l]+"/"+Region[l-1] , "pl");
	    fHistSpectrum_ratio[l][m]->Divide( fHistSpectrum_master[l-1][m]);
	    for(Int_t b=1; b <=  fHistSpectrum_ratio[l][m]->GetNbinsX(); b++){
	      fHistSpectrum_ratio[l][m]->SetBinError(b, sqrt(TMath::Abs(pow(fHistSpectrum_master[l-1][m]->GetBinError(b),2) -pow(fHistSpectrum_master[l][m]->GetBinError(b),2)))/fHistSpectrum_master[l-1][m]->GetBinContent(b)); //error calculated according to Barlow's prescription (SkipAssoc sample i a subsample of AllAssoc)
	      if (b==1) {fHistSpectrum_ratio[l][m]->SetBinContent(b,0); fHistSpectrum_ratio[l][m]->SetBinError(b,0);} //first bin (p< 0.5) is meaningless
	      if (fHistSpectrum_ratio[l][m]->GetBinContent(b) ==1)  fHistSpectrum_ratio[l][m]->SetBinError(b,0);
	    }
	  }
	}
	else {
	  //	  fHistSpectrum_ratio[l][m]->Sumw2(); //maybe already done?
	  if (l==2 || l==3){
	    if (m==0 && HistoType==0) legendRatio->AddEntry(  fHistSpectrum_ratio[l][m],"Ratio " + Region[l]+"/"+Region[l-2] , "pl");
	    fHistSpectrum_ratio[l][m]->Divide( fHistSpectrum_master[l-2][m]);
	  }
	}
	//	if (m!=nummolt || (m==nummolt && type==0)){
	//if (m!=nummolt){
	if (kTRUE){
	  if (!CfrSkipAssoc)     {
	    fHistSpectrum_ratio[l][m] ->GetYaxis()->SetRangeUser(0,10);
	    if (NameHisto[HistoType]=="SEffCorr") 	fHistSpectrum_ratio[l][m] ->GetYaxis()->SetRangeUser(0.8, 1.2);
	  }
	  else{
	    fHistSpectrum_ratio[l][m] ->GetYaxis()->SetRangeUser(0,2);
	    if (NameHisto[HistoType]=="SEffCorr")  	fHistSpectrum_ratio[l][m] ->GetYaxis()->SetRangeUser(0.4, 1.5);
	  }

	  if (CfrDiffPeriods){
	    fHistSpectrum_ratio[l][m] ->GetYaxis()->SetRangeUser(0,1.5);
	    if (l==2 || l==3)	fHistSpectrum_ratio[l][m] ->Draw("same ep");
	  }
	  else if (CfrSkipAssoc){
	    if (l==1 || l==3 ){
	      fHistSpectrum_ratio[l][m] ->Draw("same ep");
	      lineOne ->Draw("");
	    }
	  }
	  else{
	    if (l==2 || l==3 ||  NameHisto[HistoType]=="SEffCorr"){
	      fHistSpectrum_ratio[l][m] ->Draw("same ep");
	      lineOne ->Draw("");
	    }
	  }
	}
	if (m==nummolt){
	  hspectrumCfrFioComparisonMB = (TH1F*) hspectrumCfr[5]->Clone("FiorellaMBSpectraRatio");
	  hspectrumCfrFioComparisonMB->Divide(hspectrumCfrSum);
	  ErrRatioCorr(hspectrumCfr[5], hspectrumCfrSum, hspectrumCfrFioComparisonMB, 1);
	  //	  cout <<"Hola!" <<hspectrumCfr[5]->GetBinContent(1) << "/" << hspectrumCfrSum->GetBinContent(1)<< " = " <<   hspectrumCfrFioComparisonMB->GetBinContent(1) << endl;
	  hspectrumCfrFioComparisonMB->SetLineColor(kBlue);
	  hspectrumCfrFioComparisonMB->SetMarkerColor(kBlue);
	  //	  hspectrumCfrFioComparisonMB->Draw("same ep");
	}
	//	lineOne ->Draw("");
	//	if (l==3 && NameHisto[HistoType]=="SEffCorr" && !CfrSkipAssoc)        legendRatioCompFio->Draw("");
	//	if (l==3 && NameHisto[HistoType]!="SEffCorr")        legendRatio->Draw("");
	//	if (l==3 && NameHisto[HistoType]=="SEffCorr" && CfrSkipAssoc)        legendRatio->Draw("");

      } //end loop mult
    }//endl loop histo type
  }

  //drawin the pt spectra

  for (Int_t i=0; i<2; i++){
    for (Int_t HistoType=0; HistoType<numhistoType; HistoType++){
      if (i==0 && HistoType==2) canvasSpectrum[i][HistoType]->SaveAs("PictureForNote/CfrPtSpectraDiffConfig_"+ yeardati+"_SEffCorr.pdf");
      if (i==1 && HistoType==2) canvasSpectrum[i][HistoType]->SaveAs("PictureForNote/CfrPtSpectraDiffConfig_"+ yeardati+"_SEffCorrRatio.pdf");
      if (i==0 && HistoType==1) canvasSpectrum[i][HistoType]->SaveAs("PictureForNote/CfrPtSpectraDiffConfig_"+ yeardati+"_Efficiency.pdf");
      fileout->WriteTObject(canvasSpectrum[i][HistoType]);
      //      if (i==0 && HistoType==0) canvasSpectrum[i][HistoType]->SaveAs(nomefileoutputPDF+".pdf(");
      //      else   if (i!=2 && HistoType!=numhistoType-1)  canvasSpectrum[i][HistoType]->SaveAs(nomefileoutputPDF+".pdf");
      //      else  canvasSpectrum[i][HistoType]->SaveAs(nomefileoutputPDF+".pdf)");
    }
  }
  fileout->Close();
  if (isNormFactor) cout << "\n\e[35mNormalisation factor taken from file:\n" << PathNormFactor << endl;
  cout << "\n\e[35mI've produced the file:\n" << nomefileoutput << "\e[39m"<<endl;

}


