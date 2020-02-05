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

void CfrDiffPt(TString year0="2016", Int_t ishhCorr=1,  TString year="2016k", Bool_t isMultCalculatedByMe=0, Int_t numDiffPtToShow=5){

  //isMultCalculatedByMe=1 is Yield histograms are to be provided as a function of Nch (by me calculated)
  TString PathIn;
  TString PathOut;
  TFile *filein;
  TFile *fileinbis;
  TFile *fileout;
 
  const Int_t nummolt=5;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;
  const Int_t numDiffPt=8;

  if (numDiffPtToShow>numDiffPt) {
    cout << " You chose too many PtTrig values, numDiffPtToShow should be smaller" << endl;
    return;
  }
  //  TString DiffPt[numDiffPt]={"MasterThesis", "3 GeV/c"};
  Float_t PtTrigMin[numDiffPt]={3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  Float_t ColorPt[numDiffPt]={1,801,2, 909,881,  868, 3,4};
  Int_t Color[6]={2,2,418,418, 4, 4};
  Int_t ColorBis[6]={2,418,4,2,418,4};
  TString Sys[6]={"Jet DATA", "Jet MC", "OJ DATA", "OJ MC", "J+OJ DATA", "J+OJ MC"};
  TString SysErr[2]={"syst. DATA", "syst. MC"};
  TString StatErr[2]={"stat. DATA", "stat. MC"};
  TString StatErrBis[3]={"stat. in-jet", "stat. out-of-jet", "stat. inclusive"};  
  TString SysErrBis[3]={"sist. in-jet", "sist. out-of-jet", "sist. inclusive"};
  TString SysErrMB[2]={"syst. DATA 0-100 %", "syst. MC 0-100 %"};
  TString StatErrMB[2]={"stat. DATA 0-100 %", "stat. MC 0-100 %"};
  TString Region[3]={"Jet", "Out of jet", "All"};
  TString RegionFile[3]={"Jet", "OutOfJet", "All"};
  TString DATAorMC[2]={"DATA", "MC"};


  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  TString SmoltBis[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  Double_t MultTab[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  TString SPtV0[numPtV0]={"0-1", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString SPeriod[numPeriod]={"2016", "2017"};
  TString SRun[numPeriod]={"2016k", "2017h"};
  Int_t Marker[6]={33,4,33,4, 33, 4};
  Int_t MarkerOrigin[6]={20,21,33,20,21,33};
  Int_t MarkerBis[2]={23,24};
  Int_t MarkerTris[2]={31, 33};
  Int_t MarkerMBStat[6]={20,4,20,4, 20, 4};
  Int_t MarkerMBSist[6]={21,25,21,25, 21, 25};
  Double_t   mult[nummolt+1]={21.2, 16.17, 11.4625, 7.135, 3.33, 6.94};
  Double_t   MeanJet[nummolt+1][numDiffPt+1]={0};
  Double_t   Mean[nummolt+1][numDiffPt+1]={0};
  Double_t   ErrorMeanJet[nummolt+1][numDiffPt+1]={0};
  Double_t   MeanNoJet[nummolt+1]={0};
  Double_t   ErrorMeanNoJet[nummolt+1]={0};
  Int_t   NTrigger[nummolt+1][numDiffPt+1]={0};
  TF1 *pol1[numDiffPt+1];
  TH1D*  fHistYieldDistance[6];
  TH1D*  fHistYield[6];

  TH1D*  fHistYieldvsErrSoloStat[6];
  TH1F*  fHistYieldvsErrSoloStatNewMult[6];
  TH1D*  fHistYieldvsErrSoloStatMaster[6];
  TH1D*  fHistYieldvsErrSoloSist[6];
  TH1D*  fHistYieldvsErrSoloStatMB[6];
  TH1D*  fHistYieldvsErrSoloSistMB[6];
  TH1D*  fHistYieldvsErrTotal[6];

  TH1D*  fHistYieldvsErrSoloStatRatio[6];
  TH1D*  fHistYieldvsErrSoloSistRatio[6];
  TH1D*  fHistYieldvsErrSoloStatRatioMB[6];
  TH1D*  fHistYieldvsErrSoloSistRatioMB[6];
  TH1D*  fHistYieldvsErrTotalRatio[6];

  TH1D*  fHistYieldvsErrErrSoloStat[6];
  TH1D*  fHistYieldvsErrErrSoloSist[6];
  TH1D*  fHistYieldvsErrErrSoloStatMB[6];
  TH1D*  fHistYieldvsErrErrSoloSistMB[6];
  TH1D*  fHistYieldvsErrErrTotal[6];

  TH1D*    fHistYieldDatiPubblicati;
  TH1D*    fHistYieldErroriDatiPubblicati;

  TH1D* fHistChargedParticlesInJetEvents[numDiffPt];
  TH1D* fHistNTriggervsMult[numDiffPt];
  TH1D* fHistChargedParticlesInNoJetEvents;
  TH1D * fHistChargedvsMultJetEv1D;
  TH1D * fHistChargedvsMultNoJetEv1D;
  TH2D * fHistChargedvsMultJetEv;
  TH2D * fHistChargedvsMultNoJetEv;

  TFile *filedatipubbl= new TFile("HEPData-1574358449-v1-Table_8c.root", "");
  TDirectoryFile *dir3 = (TDirectoryFile*)filedatipubbl->Get("Table 8c");

  TFile *filecharged= new TFile("FinalOutput/AnalysisResults2016k_hhCorr.root", "");
  //  TFile *filecharged= new TFile("FinalOutput/AnalysisResults2018f1_extra_hhCorr_MCEff.root", "");
  TDirectoryFile *dir2 = (TDirectoryFile*)filecharged->Get("MyTask");
  TList *list2 = (TList*)dir2->Get("MyOutputContainer");
  TH3D *fHist3DJetEv= (TH3D*)list2->FindObject("fHistNumberChargedTrigger");
  TH3D *fHist3DNoJetEv= (TH3D*)list2->FindObject("fHistNumberChargedNoTrigger");
  TCanvas *canvasprova = new TCanvas ("canvasprova", "prova " , 1300, 1000);
  canvasprova->Divide(4,2);
  //fHist3DNoJetEv->Draw("");
  TLegend* legend[3];
  TLegend* legendtype[3];
  TLegend* legendbis[3];
  Int_t msize=3;

  TString TypeParticle[2]={"K^{0}_{S}","h"};

  auto legenderr=new TLegend(0.7,0.8, 0.9, 0.9);
  auto legendPtMin=new TLegend(0.1,0.7, 0.3, 0.9);
  auto legendDataMC=new TLegend(0.8,0.1, 0.9, 0.2);

  TString  PathInJet;
  TString  PathInJetMC;
  TString  PathInBulk;
  TString  PathInBulkMC;
  TString  PathInAll;
  TString  PathInAllMC;
  TString hhCorr[2]={"", "_hhCorr"};

  fHistYieldDatiPubblicati=(TH1D*)dir3->Get("Hist1D_y1");
  fHistYieldErroriDatiPubblicati=(TH1D*)dir3->Get("Hist1D_y1_e3");
  for(Int_t k=1; k < fHistYieldDatiPubblicati->GetNbinsX(); k++){
    fHistYieldDatiPubblicati->SetBinError(k,     fHistYieldErroriDatiPubblicati->GetBinContent(k));
  }

  //*********** calcolo numero di particelle di trigger per dati O MC
  TCanvas *canvastrigger = new TCanvas ("canvastrigger", "trigger " , 1300, 1000);
  TString file = year;
  TString PathInBis =  "FinalOutput/AnalysisResults" + file  + ".root";
  if (ishhCorr) PathInBis =  "FinalOutput/AnalysisResults" + file  + "_hhCorr.root";
  //  if (isMC && isEfficiency) PathInBis =  "FinalOutput/AnalysisResults" + file  + "_MCEff.root";
  fileinbis=new TFile(PathInBis,"");
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask");
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtMaxvsMultBefAll");
  TH1D *fHistTriggervsMult_MultProj;

 for (Int_t IPt = 0; IPt< numDiffPtToShow; IPt++){
   cout << "\nPt Trig Min " <<   PtTrigMin[IPt] << endl;
    fHistNTriggervsMult[IPt]= new TH1D(Form("fHistNTriggervsMult_PtMin%.1f", PtTrigMin[IPt]), "fHistNTriggervsMult", 300,0,30);
fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMin[IPt]+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(30-0.00001) );
  for(Int_t m=0; m<nummolt+1; m++){
    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	NTrigger[m][IPt]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	NTrigger[m][IPt]+= fHistTriggervsMult_MultProj->GetBinContent(j);
      }
    }
    cout << "n trigger in mult range (all triggers)    " << m << "  " <<  NTrigger[m][IPt] <<   endl;
  fHistNTriggervsMult[IPt]->SetBinContent(fHistNTriggervsMult[IPt]->GetXaxis()->FindBin(mult[m]),  NTrigger[m][IPt]);
  }

    fHistNTriggervsMult[IPt]->GetYaxis()->SetRangeUser(0, 1500000);
    fHistNTriggervsMult[IPt]->SetMarkerColor(ColorPt[IPt]);
    fHistNTriggervsMult[IPt]->SetLineColor(ColorPt[IPt]);
    fHistNTriggervsMult[IPt]->SetMarkerStyle(Marker[0]);
    fHistNTriggervsMult[IPt]->SetMarkerSize(3); 
    legendPtMin->AddEntry(fHistNTriggervsMult[IPt], Form("%.1f GeV",PtTrigMin[IPt]), "pl");
    canvastrigger->cd();
    //   fHistNTriggervsMult[IPt]->Fit(pol1[IPt], "R+");
    //    if (IPt==numDiffPt-1) legendPtMin->Draw("same");
    if (IPt==2) legendPtMin->Draw("same");
    fHistNTriggervsMult[IPt]    -> Draw("samep");
 }

  //*******calcolo numero particelle cariche in ogni intervallo di molteplicità **********
  TCanvas *canvasChargedvsMultJet = new TCanvas ("ChvsMultJet", "ChvsMultJet", 1300, 1000);
  TCanvas *canvasChargedvsMultNoJet = new TCanvas ("ChvsMultNoJet", "ChvsMultNoJet", 1300, 1000);
  fHist3DNoJetEv->GetZaxis()->SetRange(0, 30);
  fHistChargedvsMultNoJetEv=(TH2D*)fHist3DNoJetEv->Project3D("yxoe");
  fHistChargedvsMultNoJetEv-> SetName("fHistChargedvsMultNoJetEv");
  fHistChargedParticlesInNoJetEvents= new TH1D("fHistChargedParticlesInNoJetEvents", "fHistChargedParticlesInNoJetEvents", 300,0,30);

   for (Int_t IPt=0; IPt<numDiffPtToShow; IPt++){
     //    cout << "\n\nPtTrigMin " << PtTrigMin[IPt] << endl;
    fHist3DJetEv->GetZaxis()->SetRangeUser(PtTrigMin[IPt]+0.001, 30);
    fHistChargedvsMultJetEv=(TH2D*)fHist3DJetEv->Project3D( "yxoe");
    fHistChargedvsMultJetEv-> SetName(Form("fHistChargedvsMultJetEv_MinPtTrig%.1f", PtTrigMin[IPt]));
    canvasprova->cd(IPt+1);
    fHistChargedvsMultJetEv->Draw("colz");
    
    fHistChargedParticlesInJetEvents[IPt]= new TH1D(Form("fHistChargedParticlesInJetEvents_PtMin%.1f", PtTrigMin[IPt]), "fHistChargedParticlesInJetEvents", 300,0,30);

    for (Int_t m=0; m< nummolt+1; m++){
      if (m< nummolt){
	fHistChargedvsMultJetEv1D=(TH1D*)fHistChargedvsMultJetEv->ProjectionY("fHistChargedvsMultJetEv1D", fHistChargedvsMultJetEv->GetXaxis()->FindBin(Nmolt[m]+0.001), fHistChargedvsMultJetEv->GetXaxis()->FindBin(Nmolt[m+1]-0.001), "E");
	if (IPt==0)	fHistChargedvsMultNoJetEv1D =(TH1D*)fHistChargedvsMultNoJetEv->ProjectionY("fHistChargedvsMultNoJetEv1D", fHistChargedvsMultNoJetEv->GetXaxis()->FindBin(Nmolt[m]+0.001), fHistChargedvsMultNoJetEv->GetXaxis()->FindBin(Nmolt[m+1]-0.001), "E");
      }
      else {
	fHistChargedvsMultJetEv1D=(TH1D*)fHistChargedvsMultJetEv->ProjectionY("fHistChargedvsMultJetEv1D", 0, 100, "E");
	if (IPt==0)	fHistChargedvsMultNoJetEv1D =(TH1D*)fHistChargedvsMultNoJetEv->ProjectionY("fHistChargedvsMultNoJetEv1D", 0,100, "E");
      }
	// canvaschargedTrig->cd(m+1);      
	// fHistChargedvsMultJetEv1D->Draw("");
	// canvaschargedTrig->cd(m+1);      
	// fHistChargedvsMultJetEv1D->Draw("");
	
	//cout << "content of 10th bin " << fHistChargedvsMultJetEv1D->GetBinContent(10) << endl;
	MeanJet[m][IPt]=	fHistChargedvsMultJetEv1D->GetMean(1);
	ErrorMeanJet[m][IPt]=	fHistChargedvsMultJetEv1D->GetMeanError(1);
	cout << "mean number of Nch (by me calculated) for mult " << Nmolt[m] << ":  " << MeanJet[m][IPt] << endl;
	if (IPt==0){
	MeanNoJet[m]= fHistChargedvsMultNoJetEv1D->GetMean(1);
	ErrorMeanNoJet[m]=	fHistChargedvsMultNoJetEv1D->GetMeanError(1);
	cout << "mean number of Nch (by me calculated, in no-jet events, i.e. with pT,Trig < 3) for mult " << Nmolt[m] << ":  " << MeanNoJet[m] << endl;
	}

	fHistChargedParticlesInJetEvents[IPt]->SetBinContent(fHistChargedParticlesInJetEvents[IPt]->GetXaxis()->FindBin(mult[m]),  MeanJet[m][IPt]);
	fHistChargedParticlesInJetEvents[IPt]->SetBinError(fHistChargedParticlesInJetEvents[IPt]->GetXaxis()->FindBin(mult[m]),  ErrorMeanJet[m][IPt]);
	if (IPt==0)	fHistChargedParticlesInNoJetEvents->SetBinContent(fHistChargedParticlesInNoJetEvents->GetXaxis()->FindBin(mult[m]),  MeanNoJet[m]);
	if (IPt==0)	fHistChargedParticlesInNoJetEvents->SetBinError(fHistChargedParticlesInNoJetEvents->GetXaxis()->FindBin(mult[m]),  ErrorMeanNoJet[m]);
      }

    fHistChargedParticlesInJetEvents[IPt]->GetXaxis()->SetTitle("dN_{ch}/d#eta_{|#eta|<0.5}");
    fHistChargedParticlesInJetEvents[IPt]->GetYaxis()->SetTitle("N_{ch} (da me calcolato)");
    fHistChargedParticlesInJetEvents[IPt]->GetYaxis()->SetRangeUser(0, 50);
    fHistChargedParticlesInJetEvents[IPt]->SetMarkerColor(ColorPt[IPt]);
    fHistChargedParticlesInJetEvents[IPt]->SetLineColor(ColorPt[IPt]);
    fHistChargedParticlesInJetEvents[IPt]->SetMarkerStyle(Marker[0]);
    fHistChargedParticlesInJetEvents[IPt]->SetMarkerSize(3); 
    pol1[IPt]=new TF1(Form("pol1_PtTrigMin%.1f", PtTrigMin[IPt]), "pol1", 0, 30);
    pol1[IPt]->SetLineColor(ColorPt[IPt]);

    //    legendPtMin->AddEntry(fHistChargedParticlesInJetEvents[IPt], Form("%.1f GeV",PtTrigMin[IPt]), "pl");
    canvasChargedvsMultJet->cd();
    //   fHistChargedParticlesInJetEvents[IPt]->Fit(pol1[IPt], "R+");
    fHistChargedParticlesInJetEvents[IPt]    -> Draw("samep");

  }

    fHistChargedParticlesInNoJetEvents->GetYaxis()->SetRangeUser(0, 50);
    fHistChargedParticlesInNoJetEvents->SetMarkerColor(kGray);
    fHistChargedParticlesInNoJetEvents->SetLineColor(kGray);
    fHistChargedParticlesInNoJetEvents->SetMarkerStyle(Marker[0]);
    fHistChargedParticlesInNoJetEvents->SetMarkerSize(3);
    pol1[numDiffPt]=new TF1("pol1_NoJet", "pol1", 0, 30);
    pol1[numDiffPt]->SetLineColor(kGray);
    //    fHistChargedParticlesInNoJetEvents->Fit(pol1[numDiffPt], "R+");
    // cout << "ok " << endl;
    legendPtMin->AddEntry(fHistChargedParticlesInNoJetEvents, "NoJetEvent", "pl");
    canvasChargedvsMultJet->cd();
    legendPtMin->Draw("same");
    //    canvasChargedvsMultJet->cd();
    fHistChargedParticlesInNoJetEvents    -> Draw("samep"); 

  
  //********cfr data MC separatamente per J, OJ, J+OJ*********************************************************
  TCanvas *canvas[3];
  TCanvas *canvaserr[3];
  TCanvas *canvasdistance[3];
  Int_t col=0;
  for(Int_t l=0; l<3; l++){ //loop sulla region (jet, OJ, inclusive)
    Int_t count=0;
    canvas[l] = new TCanvas(Form("canvas%i", l), Region[l], 1300, 1000);
    canvaserr[l] = new TCanvas(Form("canvaserr%i", l), "Errori rtelativi_" + Region[l], 1300, 1000);
  canvasdistance[l] = new TCanvas(Form("canvasdistance%i", l), "Differenze percentuali_" + Region[l], 1300, 1000);
    for (Int_t IPt=0; IPt<numDiffPtToShow; IPt++){
      PathInJet="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +Form("_JetData_PtMin%.1f.root", PtTrigMin[IPt]);	 
      PathInJetMC="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +Form("_JetMC_PtMin%.1f.root", PtTrigMin[IPt]);	 
      PathInBulk="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +Form("_BulkData_PtMin%.1f.root", PtTrigMin[IPt]); 	 
      PathInBulkMC="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +Form("_BulkMC_PtMin%.1f.root", PtTrigMin[IPt]);
      PathInAll="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +Form("_AllData_PtMin%.1f.root", PtTrigMin[IPt]); 	 
      PathInAllMC="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +Form("_AllMC_PtMin%.1f.root", PtTrigMin[IPt]);  

      /*
      if (IPt==0){
	PathInJet="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_Jet.root";	 
	PathInJetMC="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_JetMC.root";	 
	PathInBulk="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_Bulk.root"; 	 
	PathInBulkMC="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_BulkMC.root";
	PathInAll="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_All.root"; 	 
	PathInAllMC="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_AllMC.root";  
      }
      */
      for(Int_t j=0; j<1; j++){ //loop data or MC
	if(l==0){
	  if (j==0)   PathIn=PathInJet;
	  if (j==1)   PathIn=PathInJetMC;
	  col=j;
	}
	if(l==1){
	  if (j==0)  PathIn=PathInBulk;
	  if (j==1)  PathIn=PathInBulkMC;
	  col=j+2;
	}
	if (l==2){
	  if (j==0)  PathIn=PathInAll;
	  if (j==1)  PathIn=PathInAllMC;
	  col=j+4;
	}
	filein=new TFile(PathIn, "");
	// fHistYieldvsErrSoloSistMB[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSistMB");
	// fHistYieldvsErrSoloStatMB[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStatMB");
	fHistYieldvsErrSoloStat[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStat");
	fHistYieldvsErrSoloStatNewMult[j]=new TH1F("fHistYieldvsErrSoloStatNewMult", "fHistYieldvsErrSoloStatNewMult", 1000, 0,100);

	cout << "\n\nPt index: " << IPt << endl; 
	for(Int_t m=0; m< nummolt+1; m++){
	  Mean[m][IPt]= MeanJet[m][IPt];
	  cout << "mean " << Mean[m][IPt]<< endl;	
  fHistYieldvsErrSoloStatNewMult[j]->SetBinContent(fHistYieldvsErrSoloStatNewMult[j]->GetXaxis()->FindBin(Mean[m][IPt]), 	  fHistYieldvsErrSoloStat[j]->GetBinContent(fHistYieldvsErrSoloStat[j]->GetXaxis()->FindBin(MultTab[m])));
  fHistYieldvsErrSoloStatNewMult[j]->SetBinError(fHistYieldvsErrSoloStatNewMult[j]->GetXaxis()->FindBin(Mean[m][IPt]), 	  fHistYieldvsErrSoloStat[j]->GetBinError(fHistYieldvsErrSoloStat[j]->GetXaxis()->FindBin(MultTab[m])));
  cout << "region under study " << l << "bin content of mult : " <<m << "  " <<  	  fHistYieldvsErrSoloStatNewMult[j]->GetBinContent(fHistYieldvsErrSoloStatNewMult[j]->GetXaxis()->FindBin(Mean[m][IPt]))<< endl;
	}

	fHistYieldDistance[j]=(TH1D*)fHistYieldvsErrSoloStat[j]->Clone("Relative difference");
	if (IPt==0) fHistYieldvsErrSoloStatMaster[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloStat");
	fHistYieldvsErrErrSoloStat[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloStat");
	fHistYieldvsErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrSoloSist");
	fHistYieldvsErrErrSoloSist[j]=(TH1D*)filein->Get("fHistYieldvsErrErrSoloSist");
	
	for (Int_t bin=1; bin<=fHistYieldDistance[j]->GetNbinsX(); bin++){
	  if (fHistYieldvsErrSoloStatMaster[j]->GetBinContent(bin) ==0) continue;
	  //old def	  if (IPt>0)    fHistYieldDistance[j]->SetBinContent(bin,(fHistYieldvsErrSoloStat[j]->GetBinContent(bin)-fHistYieldvsErrSoloStatMaster[j]->GetBinContent(bin))/fHistYieldvsErrSoloStatMaster[j]->GetBinContent(bin) );
	  if (IPt>0)    fHistYieldDistance[j]->SetBinContent(bin,fHistYieldvsErrSoloStat[j]->GetBinContent(bin)/fHistYieldvsErrSoloStatMaster[j]->GetBinContent(bin));
	  else  fHistYieldDistance[j]->SetBinContent(bin,0);
	  cout << "Pt " << IPt << " bin " << bin << "  "<< fHistYieldDistance[j]->GetBinContent(bin) << endl;;

	}
	fHistYieldvsErrSoloStat[j]->SetTitle ("Relative yield of " + TypeParticle[ishhCorr]+ " in " + Region[l] + " per trigger particle vs V0M multiplicity");
	fHistYieldvsErrSoloStat[j]->GetXaxis()->SetTitle ("<dN_{ch}/d#eta>_{|#eta|<0.5}");
	fHistYieldvsErrSoloStat[j]->GetYaxis()->SetTitle ("N_{"+TypeParticle[ishhCorr]+"}/N_{Trigg} 1/#Delta#eta #Delta#phi");
	fHistYieldvsErrSoloStat[j]->SetMarkerColor(ColorPt[IPt]);
	fHistYieldvsErrSoloStat[j]->SetMarkerStyle(Marker[j]);
	fHistYieldvsErrSoloStat[j]->SetLineColor(ColorPt[IPt]);
	fHistYieldvsErrSoloStat[j]->GetYaxis()->SetRangeUser(0,0.4);
	if (l==0) fHistYieldvsErrSoloStat[j]->GetYaxis()->SetRangeUser(0,0.2);
	fHistYieldvsErrSoloStat[j]->GetXaxis()->SetRangeUser(0,25);
	fHistYieldvsErrSoloStat[j]->SetMarkerSize(msize);
	if (j==0)    fHistYieldvsErrSoloStat[j]->SetMarkerSize(3);

	fHistYieldvsErrSoloStatNewMult[j]->SetTitle ("Yield of "+ TypeParticle[ishhCorr]+" in " + Region[l] + " per trigger particle vs N_{charged}");
	fHistYieldvsErrSoloStatNewMult[j]->GetXaxis()->SetTitle ("N_{charged}");
	fHistYieldvsErrSoloStatNewMult[j]->GetYaxis()->SetTitle ("N_{"+ TypeParticle[ishhCorr]+"}/N_{Trigg} 1/#Delta#eta #Delta#phi");
	fHistYieldvsErrSoloStatNewMult[j]->SetMarkerColor(ColorPt[IPt]);
	fHistYieldvsErrSoloStatNewMult[j]->SetMarkerStyle(Marker[j]);
	fHistYieldvsErrSoloStatNewMult[j]->SetLineColor(ColorPt[IPt]);
	fHistYieldvsErrSoloStatNewMult[j]->GetYaxis()->SetRangeUser(0,0.4);
	if (l==0) fHistYieldvsErrSoloStatNewMult[j]->GetYaxis()->SetRangeUser(0,0.2);
	fHistYieldvsErrSoloStatNewMult[j]->GetXaxis()->SetRangeUser(0,50);
	fHistYieldvsErrSoloStatNewMult[j]->SetMarkerSize(msize);
	if (j==0)    fHistYieldvsErrSoloStatNewMult[j]->SetMarkerSize(3);

	canvas[l]-> cd();
	if (isMultCalculatedByMe)	fHistYieldvsErrSoloStatNewMult[j]->Draw("samee");
	else		fHistYieldvsErrSoloStat[j]->Draw("samee");

	fHistYieldvsErrErrSoloStat[j]->SetMarkerColor(ColorPt[IPt]);
	fHistYieldvsErrErrSoloStat[j]->SetMarkerStyle(Marker[j]);
	fHistYieldvsErrErrSoloStat[j]->SetLineColor(ColorPt[IPt]);
	fHistYieldvsErrErrSoloStat[j]->SetMarkerSize(msize);
	canvaserr[l]-> cd();
	fHistYieldvsErrErrSoloStat[j]->Draw("samep");

    	fHistYieldvsErrSoloSist[j]->SetTitle("");
	fHistYieldvsErrSoloSist[j]->SetMarkerColor(ColorPt[IPt]);
	fHistYieldvsErrSoloSist[j]->SetMarkerStyle(Marker[j]);
	fHistYieldvsErrSoloSist[j]->SetLineColor(ColorPt[IPt]);
	fHistYieldvsErrSoloSist[j]->SetFillStyle(0);
	fHistYieldvsErrSoloSist[j]->GetYaxis()->SetRangeUser(0,0.4);
	if (l==0) fHistYieldvsErrSoloSist[j]->GetYaxis()->SetRangeUser(0,0.2);
	fHistYieldvsErrSoloSist[j]->GetXaxis()->SetRangeUser(0,25);
	fHistYieldvsErrSoloSist[j]->SetMarkerSize(msize);
	if (j==0)    fHistYieldvsErrSoloSist[j]->SetMarkerSize(3);
	canvas[l]-> cd();
	//	fHistYieldvsErrSoloSist[j]->Draw("samee2");

	fHistYieldvsErrErrSoloSist[j]->SetMarkerColor(ColorPt[IPt]);
	fHistYieldvsErrErrSoloSist[j]->SetMarkerStyle(Marker[j+1]);
	fHistYieldvsErrErrSoloSist[j]->SetLineColor(ColorPt[IPt]);
	fHistYieldvsErrErrSoloSist[j]->SetMarkerSize(msize);
	canvaserr[l]-> cd();
	//	fHistYieldvsErrErrSoloSist[j]->Draw("samep");

	fHistYieldDistance[j]->SetMarkerColor(ColorPt[IPt]);
	fHistYieldDistance[j]->SetMarkerStyle(Marker[j]);
	fHistYieldDistance[j]->SetLineColor(ColorPt[IPt]);
	fHistYieldDistance[j]->SetMarkerSize(msize);
	if (l==0)	fHistYieldDistance[j]->GetYaxis()->SetRangeUser(1,1.7);
	else fHistYieldDistance[j]->GetYaxis()->SetRangeUser(1,1.2);
	cout << "Pt " << IPt << " bin " << 34 << "  "<< fHistYieldDistance[j]->GetBinContent(34) << endl;;
	canvasdistance[l]-> cd();
	//	if (IPt>0) fHistYieldDistance[j]->Draw("samep");      
	fHistYieldDistance[j]->Draw("samep");      
	//	if (j==0 && l==0)      legendPtMin->AddEntry(fHistYieldvsErrSoloStat[j], Form("%.1f GeV",PtTrigMin[IPt]), "pl");
	if (IPt==0 && l==0)      legendDataMC->AddEntry(fHistYieldvsErrSoloStat[j], DATAorMC[j], "pl");
	if (IPt==0 && j==0 && l==0)      legenderr->AddEntry(fHistYieldvsErrSoloStat[j],"stat" , "el");	
	if (IPt==0 && j==0 && l==0)      legenderr->AddEntry(fHistYieldvsErrSoloSist[j],"sist" , "f");

	//	if (IPt==numDiffPt-1){
	if (IPt==0){
	canvas[l]-> cd();
	legendPtMin->Draw("same");
	legendDataMC->Draw();
	legenderr->Draw();

	canvaserr[l]-> cd();
	legendPtMin->Draw();
	legendDataMC->Draw();
	legenderr->Draw();

	canvasdistance[l]-> cd();
	legendPtMin->Draw();
	legendDataMC->Draw();
	legenderr->Draw();

	}
      }
      canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/Yield"+RegionFile[l]+hhCorr[ishhCorr]+Form("_PtMin%.1f.pdf", PtTrigMin[IPt]));
    }
  }

  TString stringout = "FinalOutput/DATA"+year0+"/OutputCfrDiffPt_" + hhCorr[ishhCorr]+".root";
  fileout= new TFile(stringout, "RECREATE");
  for(Int_t l=0; l<3; l++){ 
    fileout->WriteTObject(canvas[l]);
    fileout->WriteTObject(canvaserr[l]);
    fileout->WriteTObject(canvasdistance[l]);
  }

  cout << "I've produced the file " << stringout << endl;
  // fileout->WriteTObject(canvasChargedvsMultJet);
  // fileout->WriteTObject(canvasChargedvsMultNoJet);
}

