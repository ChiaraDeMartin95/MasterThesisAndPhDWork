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
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

void CfrCommonParton( Bool_t ishhCorr =0, Int_t type=6, Int_t israp=0,Bool_t SkipAssoc=1 ,   Float_t ptjmin=3,  Float_t ptjmax=15, Int_t sysTrigger=0, Int_t sysV0=0, TString data="2018f1_extra_hK0s_CP"/*"1617GP_hXi"/*"2018f1_extra_MECorr"*/ /*"2018f1_extra_hK0s_30runs_150MeV"*/, TString year0="2016", TString path1="", Bool_t FitGen=1, Int_t rebin=5, Bool_t IsAllMult=0 /*0 if we wnat 0-100% only*/, Int_t PtBinning=1){

  TString PathOut;
  TFile *filein;
  TFile *fileout;

  const Int_t numhistoType=3;
  TString CPString[numhistoType]={"", "_CPEff", "_NOCPEff"};
  TString CPStringLegend[numhistoType]={"All ", "CP", "NOCP"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};
  const Int_t numtipo=10;
  Float_t massParticle[numtipo]= {1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[numtipo]={"XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega", "K0s", "h"};
  TString tipoPart[numtipo]={"Xi", "Xi", "Omega", "Omega", "Xi", "Omega", "", ""};

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;

  TString Smult[nummolt+1] = {"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  //  Int_t Marker[6]={20,21,20,21, 20, 21};
  Int_t Marker[3]={20,28,33};

  Int_t Color[6]={859, 887, 633, 807, 814, 923};
  Int_t ColorAllMult[6]={859, 887, 633, 807, 814, 923};
  Int_t ColorOneMult[6]={628, 797, 867, 1,1,1};

  for (Int_t i=0; i<nummolt+1; i++){
    if (IsAllMult)    Color[i] = ColorAllMult[i];
    else    Color[i] = ColorOneMult[i];
  }
  Int_t MarkerBis[2]={23,24};
  Int_t MarkerTris[2]={20, 33};
  Int_t MarkerMBStat[6]={20,4,20,4, 20, 4};
  Int_t MarkerMBSist[6]={21,25,21,25, 21, 25};
  Double_t PtTrigMin = ptjmin;

  TString nomefileoutput = "CfrCommonParton_"+data;
  nomefileoutput += path1;
  nomefileoutput += Form("_PtTrigMin%.1f", PtTrigMin);
  TString  nomefileoutputPDF=nomefileoutput;
  nomefileoutput += ".root";
  cout << "nome file output " << nomefileoutput << endl;

  TString nomehistoSel[3] ={"fHistSelected_1D_V0Pt_", "fHistSelected_1D_V0Eta_", "fHistSelected_1D_V0Phi_"};
  TString nomehistoGen[3] ={"fHistGenerated_1D_V0Pt_", "fHistGenerated_1D_V0Eta_", "fHistGenerated_1D_V0Phi_"};
  TString nomehistoEff[3] ={"fHistV0EfficiencyPtBins_", "fHistV0EfficiencyEta_", "fHistV0EfficiencyPhi_"};

  TH1F* HistoSel[3][nummolt+1];
  TH1F* HistoGen[3][nummolt+1];
  TH1F* HistoEff[3][nummolt+1];
  TH1F* HistoEffMultAll[3];
  TH1F* HistoEffRatioMult[3][nummolt+1];
  TH1F* HistoSelRatio[3][nummolt+1];
  TH1F* HistoGenRatio[3][nummolt+1];
  TH1F* HistoEffRatio[3][nummolt+1];
  TH1F* HistoSelDenum[3][nummolt+1];
  TH1F* HistoGenDenum[3][nummolt+1];
  TH1F* HistoEffDenum[3][nummolt+1];

  TCanvas *canvas[6];
  TCanvas *canvasPtEff[6];

  //variables for the fit 
  AliPWGFunc pwgfunc;
  Int_t numfittipo=5;
  Int_t ColorFit[numfittipo]={860, 881, 868, 628, 419};
  TLegend *legendfit=new TLegend(0.6, 0.35, 0.9, 0.6);
  TString       nameFit[numfittipo]={"mT-scaling", "pT-scaling", "Boltzmann", "Fermi-Dirac", "Levi"};//"Bose-Einstein"}; 
  TH1F* HistoGenRatioFit[3][numfittipo];
  TF1* fit_scaling[2][numfittipo];
  TString   namescaling[2][numfittipo];


  TString file = "Efficiency";
  file+=data;

  TString PathIn;
  for (Int_t var=0; var<3; var++){ //loop over pt, phi, eta
    canvas[var] = new TCanvas (Form("canvas%i",var), Form("canvas%i",var), 1300, 800);
    canvas[var]->Divide(3,2);
    canvasPtEff[var] = new TCanvas (Form("canvasPtEff%i",var), Form("canvasPtEff%i", var), 1300, 800);
    canvasPtEff[var]->Divide(3,1);
  }

  TCanvas *canvasptSpectrum = new TCanvas ("canvasptSpectrum", "canvasptSpectrum", 1300, 800);
  canvasptSpectrum->Divide(3,2);


  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  TLegend* legendmult = new TLegend(0.7, 0.5, 0.9, 0.7);
  //  legend->SetHeader("");
  TF1 * lineat1= new TF1("pol0", "pol0",0,8);
  lineat1->FixParameter(0,1);
  lineat1->SetLineColor(kBlack);


  fileout = new TFile (nomefileoutput, "RECREATE");

  for (Int_t isCP=0; isCP<numhistoType; isCP++){
  //  for (Int_t isCP=0; isCP<1; isCP++){
    PathIn="FinalOutput/DATA" + year0 + "/Efficiency/" + file + path1+ CPString[isCP]+"_"+tipo[type]+Srap[israp]+Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
    if (PtBinning==1)     PathIn="FinalOutput/DATA" + year0 + "/Efficiency/" + file + path1+ Form("_PtBinning%i", PtBinning)+CPString[isCP]+"_"+tipo[type]+Srap[israp]+Form("_SysT%i_SysV0%i_PtMin%.1f.root", sysTrigger, sysV0, PtTrigMin);
    cout << "\n " << PathIn << endl;
    filein =new TFile(PathIn, "");
    if (!filein) {cout << "input file does not exist " << endl; return;}
    for (Int_t var=0; var<3; var++){ //loop over pt, phi, eta
      cout << var << endl;
      for (Int_t m=0; m< nummolt+1; m++){
	if (!IsAllMult)	if (m!=nummolt) continue;
	//	if (m>=3 && m <=4) continue;
      //      for (Int_t m=0; m< ; m++){
      HistoSel[var][m]=(TH1F*)filein->Get(nomehistoSel[var]+Smult[m]);
      HistoGen[var][m]=(TH1F*)filein->Get(nomehistoGen[var]+Smult[m]);
      HistoEff[var][m]=(TH1F*)filein->Get(nomehistoEff[var]+Smult[m]);
      if (m==0)      HistoEffMultAll[var]=(TH1F*)filein->Get(nomehistoEff[var]+Smult[nummolt]);
      //if (m==0)      HistoEffMultAll[var]=(TH1F*)filein->Get(nomehistoEff[var]+Smult[0]);
      if (!HistoSel[var][m]) {cout << " selected histo does not exist " << nomehistoSel[var]   +Smult[m]<< endl; return;}
      if (!HistoGen[var][m]) {cout << " genrated histo does not exist " <<  nomehistoGen[var]  +Smult[m]<<endl; return;}
      if (!HistoEff[var][m]) {cout << " eff histo does not exist " <<  nomehistoEff[var]       +Smult[m]<<endl; return;}


      if (IsAllMult){
	HistoEffMultAll[var]->Sumw2();
	HistoEffRatioMult[var][m]= (TH1F*)HistoEff[var][m]->Clone(nomehistoEff[var]+Smult[m]+"_Ratio");
	HistoEffRatioMult[var][m]->Sumw2();
	HistoEffRatioMult[var][m]->Divide(HistoEffMultAll[var]);
	HistoEffRatioMult[var][m]->SetLineColor(Color[m]);
	HistoEffRatioMult[var][m]->SetMarkerColor(Color[m]);
	HistoEffRatioMult[var][m]->SetMarkerStyle(Marker[isCP]);
	HistoEffRatioMult[var][m]->GetYaxis()->SetRangeUser(0.85,1.35);
	HistoEffRatioMult[var][m]->GetYaxis()->SetTitle("");
	canvasPtEff[var]->cd(isCP+1);
	if (m==0)	HistoEffRatioMult[var][m]->Draw("p");
	else if (m!=nummolt)	HistoEffRatioMult[var][m]->Draw("p same");
	lineat1->Draw("same");
	//	if (m==nummolt) legendmult->Draw("same");
      }

      if (var==0)    HistoSel[var][m]->GetXaxis()->SetRangeUser(0,10);
      if (var==0)    HistoGen[var][m]->GetXaxis()->SetRangeUser(0,10);
      if (IsAllMult){
      HistoSel[var][m]->SetLineColor(Color[m]);
      HistoGen[var][m]->SetLineColor(Color[m]);
      HistoEff[var][m]->SetLineColor(Color[m]);
      HistoSel[var][m]->SetMarkerColor(Color[m]);
      HistoGen[var][m]->SetMarkerColor(Color[m]);
      HistoEff[var][m]->SetMarkerColor(Color[m]);
      HistoSel[var][m]->SetMarkerStyle(Marker[isCP]);
      HistoGen[var][m]->SetMarkerStyle(Marker[isCP]);
      HistoEff[var][m]->SetMarkerStyle(Marker[isCP]);
      }
      else {
	HistoSel[var][m]->SetLineColor(Color[isCP]);
      HistoGen[var][m]->SetLineColor(Color[isCP]);
      HistoEff[var][m]->SetLineColor(Color[isCP]);
      HistoSel[var][m]->SetMarkerColor(Color[isCP]);
      HistoGen[var][m]->SetMarkerColor(Color[isCP]);
      HistoEff[var][m]->SetMarkerColor(Color[isCP]);
      HistoSel[var][m]->SetMarkerStyle(33);
      HistoGen[var][m]->SetMarkerStyle(33);
      HistoEff[var][m]->SetMarkerStyle(33);
      
      }

      if(var==0){
      HistoSel[var][m]->Rebin(rebin);
      HistoGen[var][m]->Rebin(rebin);
      }
      HistoSel[var][m]->GetYaxis()->SetRangeUser(0, 1.2*HistoSel[var][m]->GetMaximum());
      HistoGen[var][m]->GetYaxis()->SetRangeUser(0, 1.2*HistoGen[var][m]->GetMaximum());
      HistoEff[var][m]->GetYaxis()->SetRangeUser(0, 1.2*HistoEff[var][m]->GetMaximum());

      //fit to the pT distribution of selected particles
      if (var==0 ) {
	for (Int_t typefit =0; typefit<numfittipo; typefit++){
	 
	  pwgfunc.SetVarType(AliPWGFunc::kdNdpt);
	  namescaling[var][typefit] = Form("fitscaling_var%i_fit%i",var, typefit);
	  cout <<"\n\nisCP " << isCP << " fit: " <<  nameFit[typefit]<< endl;
	  if (typefit==0)   fit_scaling[var][typefit]=    pwgfunc.GetMTExp(massParticle[type], 0.1, 0.04*rebin, namescaling[var][typefit]); //mass, T, norm, name                                                                                                   
	  if (typefit==1)   fit_scaling[var][typefit]=    pwgfunc.GetPTExp(0.1, 0.04*rebin, namescaling[var][typefit]); //mass, T, norm, name                                                                                                                       
	  if (typefit==2)   fit_scaling[var][typefit]=    pwgfunc.GetBoltzmann(massParticle[type],0.1, 0.04*rebin, namescaling[var][typefit]);
	  if (typefit==3)   fit_scaling[var][typefit]=    pwgfunc.GetFermiDirac(massParticle[type],0.1, 0.04*rebin, namescaling[var][typefit]);  
	  if (typefit==4)  {
	    fit_scaling[var][typefit]=    pwgfunc.GetLevi(massParticle[type],0.1, 0.03, 0.04*rebin, namescaling[var][typefit]);
	    fit_scaling[var][typefit]->SetParLimits(0, 0, 1000000*rebin);
	    fit_scaling[var][typefit]->SetParLimits(2, 0.1, 10);
	    fit_scaling[var][typefit]->SetParLimits(1, 2, 100000);
	    fit_scaling[var][typefit]->SetParameter(2, 0.7);
	  }

	  if (nameFit[typefit]== "pT-scaling") continue;

	  if (m==5 && var==0 && isCP==0)    legendfit->AddEntry( fit_scaling[var][typefit], nameFit[typefit], "l");
	  fit_scaling[var][typefit]->SetLineColor(ColorFit[typefit]);
	  if (type==4){
	  if (isCP==1)     fit_scaling[var][typefit]->SetRange(1,8);//1
	  else      fit_scaling[var][typefit]->SetRange(0.5,8);//0.5
	  }
	  else if (type==6){
	  if (isCP==1)     fit_scaling[var][typefit]->SetRange(0.5,8);//1
	  else      fit_scaling[var][typefit]->SetRange(0,8);//0.5
	  }
	  if (!FitGen)	  HistoSel[var][m]->Fit(    fit_scaling[var][typefit],"R0");
	  else HistoGen[var][m]->Fit(    fit_scaling[var][typefit],"R0");
	  fit_scaling[var][typefit]->SetRange(0,10);

	  //print some info of the fit 
	  cout <<  "chisquare/NDF " << fit_scaling[var][typefit]->GetChisquare() << "/" << fit_scaling[var][typefit]->GetNDF() <<"\n" <<  endl;                                          
	  cout << " T " << fit_scaling[var][typefit]->GetParameter(1)<< endl;                                                        
	  cout << "integral percentage for pT < 0.5 GeV/c " << fit_scaling[var][typefit]->Integral(0, 0.5)/fit_scaling[var][typefit]->Integral(0, 30)<< endl;                                                                                               
	  cout << "integral percentage for pT < 1.0 GeV/c " << fit_scaling[var][typefit]->Integral(0, 1)/fit_scaling[var][typefit]->Integral(0, 30)<< endl;                                                                                                   

	  cout << "integral percentage for pT < 1.5 GeV/c " << fit_scaling[var][typefit]->Integral(0, 1.5)/fit_scaling[var][typefit]->Integral(0, 30)<< endl;                                                                                               
     
	  cout << "Integral up to 30 GeV/c " << fit_scaling[var][typefit]->Integral(0, 30)<< endl;                          
	}
      }

      if (isCP==0){
	HistoSelDenum[var][m]=(TH1F*)       HistoSel[var][m]->Clone(nomehistoSel[var]+Smult[m]+"_Denum");
	HistoGenDenum[var][m]=(TH1F*)       HistoGen[var][m]->Clone(nomehistoGen[var]+Smult[m]+"_Denum");
	HistoEffDenum[var][m]=(TH1F*)       HistoEff[var][m]->Clone(nomehistoEff[var]+Smult[m]+"_Denum");
      }

      else{
	HistoSelRatio[var][m]=(TH1F*)       HistoSel[var][m]->Clone(nomehistoSel[var]+Smult[m]+"_Ratio");
	HistoGenRatio[var][m]=(TH1F*)       HistoGen[var][m]->Clone(nomehistoGen[var]+Smult[m]+"_Ratio");
	HistoEffRatio[var][m]=(TH1F*)       HistoEff[var][m]->Clone(nomehistoEff[var]+Smult[m]+"_Ratio");
	HistoSelRatio[var][m]->Divide(   HistoSelDenum[var][m]);
	HistoGenRatio[var][m]->Divide(   HistoGenDenum[var][m]);
	HistoEffRatio[var][m]->Divide(   HistoEffDenum[var][m]);
			  
	HistoSelRatio[var][m]->GetYaxis()->SetRangeUser(0,2);
	HistoGenRatio[var][m]->GetYaxis()->SetRangeUser(0,2);
	HistoEffRatio[var][m]->GetYaxis()->SetRangeUser(0,2);

      }

      if (m==5 && var==0) legend->AddEntry(HistoSel[var][m], CPStringLegend[isCP], "pl");
      if (IsAllMult && isCP==0 && var==0) legendmult->AddEntry(HistoSel[var][m], Smult[m], "pl");

      //      gStyle->SetOptStat("n");

      cout << " ok up to here " << endl; 

      canvas[var]->cd(1);
      HistoSel[var][m]->Draw("same p");
      if (var==0){
      for (Int_t typefit =0; typefit<numfittipo; typefit++){
	if (    nameFit[typefit]== "pT-scaling") continue;
	if (isCP==0){
	if(!FitGen) 	fit_scaling[var][typefit]->Draw("same");
	}
      }
      if (isCP==0 && m==nummolt)      legendfit->Draw("same");
      }
      if (isCP==2 && m==nummolt )     {
	legend->Draw("same");
	if ( IsAllMult)	legendmult->Draw("same");
      }

      canvas[var]->cd(2);
      HistoGen[var][m]->Draw("same p");
      if (var==0){
      for (Int_t typefit =0; typefit<numfittipo; typefit++){
		if (    nameFit[typefit]== "pT-scaling") continue;

	if (FitGen)	fit_scaling[var][typefit]->Draw("same");
      }
      if (isCP==0 && m==nummolt)      legendfit->Draw("same");
      }
            
      if (isCP==2 && m==nummolt)    {
	legend->Draw("same");
	if(	IsAllMult)	legendmult->Draw("same");
      }

      canvas[var]->cd(3);
      HistoEff[var][m]->Draw("same p ");
      if (isCP==2 && m==nummolt)  {
    legend->Draw("same");
    if ( IsAllMult)   legendmult->Draw("same");
      }

      if (isCP!=0){

	canvas[var]->cd(4);
	HistoSelRatio[var][m]->Draw("same p");
	if (m==nummolt)	legend->Draw("same");

	canvas[var]->cd(5);
	HistoGenRatio[var][m]->Draw("same p");
	if (m==nummolt)	legend->Draw("same");

	canvas[var]->cd(6);
	HistoEffRatio[var][m]->Draw("same p");
	if (m==nummolt)	legend->Draw("same");
      }

      if (var==0){
	canvasptSpectrum->cd(isCP+1);
	if (var==0)    HistoGen[var][m]->GetXaxis()->SetRangeUser(0,10);
	HistoGen[var][m]->Draw("pe");
	for (Int_t typefit =0; typefit<numfittipo; typefit++){
	  if (nameFit[typefit]== "pT-scaling") continue;
	  if (FitGen) fit_scaling[var][typefit]->Draw("same");
	}
	if (m==nummolt)	legendfit->Draw("same");

	canvasptSpectrum->cd(isCP+1+3);
	for (Int_t typefit =0; typefit<numfittipo; typefit++){
	  if (nameFit[typefit]== "pT-scaling") continue;
	HistoGenRatioFit[var][typefit]=(TH1F*)       HistoGen[var][m]->Clone(nomehistoGen[var]+Form("_RatioFit_%i", typefit));
	HistoGenRatioFit[var][typefit]->Divide(fit_scaling[var][typefit]);
	HistoGenRatioFit[var][typefit]->GetYaxis()->SetRangeUser(0,3);
	HistoGenRatioFit[var][typefit]->SetLineColor(ColorFit[typefit]);
	HistoGenRatioFit[var][typefit]->SetMarkerColor(ColorFit[typefit]);
	HistoGenRatioFit[var][typefit]->Draw("p same");
	}
      }

      fileout->WriteTObject(      HistoSel[var][m]);
      fileout->WriteTObject(      HistoGen[var][m]);

      fileout->WriteTObject(      HistoEff[var][m]);
      
      fileout->WriteTObject(      HistoSelRatio[var][m]);
     
      fileout->WriteTObject(      HistoGenRatio[var][m]);
      fileout->WriteTObject(      HistoEffRatio[var][m]);
     
    } //end of var
    } //end of mult
  } //end of isCP


  for (Int_t var=0; var<3; var++){ //loop over pt, phi, eta            
    fileout->WriteTObject(canvas[var]);    

  }
    fileout->WriteTObject(canvasptSpectrum);    
  fileout->Close();
  cout << " I've produced the file " << nomefileoutput << endl;

}


