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
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"
#include </data/dataalice/AliceSoftware/aliphysics/master_chiara/src/PWG/Tools/AliPWGFunc.h>

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Bool_t XRange, Float_t XLow, Float_t XUp){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  if (XRange)   histo->GetXaxis()->SetRangeUser(XLow, XUp);
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

void MCClosureCheck( Int_t isEventLoss=2, Int_t type=0, Int_t sys=0, Int_t sysPhi=0, Int_t ishhCorr=0, Float_t PtTrigMin =3, Float_t PtTrigMax=15, Int_t israp=0,TString yearMC="2018g4_extra_EtaEff_hK0s"/*"1617_hK0s"/*"AllMC_hXi"/"2018f1_extra_Reco_hK0s"/*"2016k_hK0s"/"Run2DataRed_MECorr_hXi"/*"2016k_hK0s_30runs_150MeV"/*"  2016k_New"/"Run2DataRed_hXi"/*"2016kehjl_hK0s"*/, TString yearMCTruth="2018f1_extra_hK0s", TString yearMCHyb="2018g4_extra_EtaEff_Hybrid_hK0s"/*"2018f1_extra_Hybrid_hK0s"*/, TString Path1 =""/*"_Jet0.75"/*"_Jet0.75"/*"_NewMultClassBis_Jet0.75"*/,  TString Dir="FinalOutput",TString year0="2016", Bool_t SkipAssoc=0, Int_t MultBinning=0, Int_t PtBinning=1, Bool_t IsParticleTrue=0, Bool_t isEffMassSel=0, Bool_t IsNotSigmaTrigger=0, Bool_t isMEFromHybrid=1, Bool_t isMEFromCorrectCentrality=0, Bool_t isEtaEff=1, Bool_t isMEFromK0s=0, Bool_t isFioComp=0, Bool_t ispp5TeV=0, Bool_t isppHM =0, Bool_t isBulkBlue=0, Bool_t isTrigEff=1, Bool_t isSidebands =0, Bool_t isBkgParab=0){

  if (isEventLoss==1) IsParticleTrue=0;
  if (isMEFromK0s && type==0) {cout << "the option isMEFromK0s is meant to be used when Xi is being analyzed" << endl; return;}
  if (sys!=0 && sys!=4 && sys!=5 && sys!=6) {cout << "Not implemented for that sys value"<< endl; return;}
  if (sysPhi>2) {cout << "Not implemented for that sysPhi value"<< endl; return;}
  // if isEventLoss==0 I compare (MCreco/MCTruth) ---> for closure test as Lucia wants
  // if isEventLoss==1 I compare (MCHybrid/MCTruth) ---> for norm factor
  // if isEventLoss==2 I compare (MCReco/MCHybrid) ---> for closure test

  if (isEventLoss==1) {
    isEtaEff = 0;
    isMEFromHybrid=0;
  }

  Bool_t isRecoNormCorr=0;
  if (isEventLoss==0){
    cout << "Do you want to use the reco spectra corrected by the normalisation factor?" << endl;
    cin >> isRecoNormCorr;
  }

  //********** K0s *******************
  if (type==0){ //used for MC closure for preliminary
    //    yearMC = "2018f1_extra_Reco_hK0s";
    //    yearMCHyb = "2018f1_extra_Hybrid_hK0s";
    /*
    Used for 2020 Preliminaries
    yearMC = "2018g4_extra_EtaEff_hK0s";
    yearMCHyb="2018g4_extra_EtaEff_Hybrid_hK0s";
    */
    if (isEventLoss!=2){
      yearMCHyb = "1617GP_hK0s_Hybrid_New"; //norm factor
      yearMCTruth = "1617GP_hK0s"; //norm factor
    }
    //the following are used as tests:
    //    yearMCHyb = "2018f1_extra_MylabelBis_15runs_hK0s_Hybrid"; 
    //    yearMCHyb = "2018f1_extra_RlabelBis_15runs_hK0s_Hybrid"; 
    //    yearMCHyb = "2018f1_extra_RlabelCorrected_15runs_hK0s_Hybrid"; 
    //    yearMCTruth = "2018f1_extra_15runs"; 
    //    yearMCHyb = "2018f1_extra_15runs_NohDaughtersofK0s_hK0s_Hybrid";
   
    if (isFioComp){
      yearMCHyb = "LHC16kl_pass2_GP_Fio_Hybrid"; //norm factor
      yearMCTruth = "LHC16kl_pass2_GP_Fio"; //norm factor
      PtTrigMin =0.15;
    }
    if (ispp5TeV){
      /*
      yearMCHyb = "17pq_pp5TeV_Hybrid"; //norm factor
      yearMCTruth = "17pq_pp5TeV"; //norm factor
      */
      yearMC = "17l3b_hK0s";
      yearMCHyb = "17l3b_hK0s_Hybrid";
      yearMCTruth = "17l3b_hK0s";
    }
    else if (isppHM){
      yearMC = "161718HM_hK0s"; 
      yearMCHyb="161718HM_hK0s_Hybrid";
      yearMCTruth="161718HM_hK0s";
      isBkgParab =1;
    }
    else {
    //Used for 2022 Preliminaries (for MC closure only)
    yearMC = "16kl_hK0s"; //43 M events
    yearMCHyb="16kl_hK0s_Hybrid";
    yearMCTruth="16kl_hK0s";
    }
  }
  //**********************************

  //********** Xi *******************
  if (type==8){
    if (PtTrigMin==3){    
      //    yearMC = "AllMC_hXi";
      yearMC = "2018f1g4_extra_EtaEff_hXi";
      if (isEventLoss==2)      yearMCHyb="2018f1g4_extra_hXi_Hybrid";
      //OLD - not enough stat      if (isEventLoss==1)      yearMCHyb="2018f1g4_extra_hXi_Hybrid";
      //OLD - not enough stat      if (isEventLoss==1)      yearMCHyb="161718_MD_New_hXi_Hybrid";
      if (isEventLoss==1)      yearMCHyb="161718_hXi_Hybrid";
    }
    //    yearMCHyb="LHC16_17_GP_Hybrid_hXi";
    else {    
      yearMCHyb="161718_MD_hXi_Hybrid";
      yearMC="161718_MD_EtaEff_hXi";
    }
    yearMCTruth = "161718_hXi";
    //OLD not enough stat yearMCTruth = "161718_MD_New_hXi";
    //OLD not enough stat yearMCTruth="2018g4_extra_hXi_SelTrigger";
    PtBinning=0;
    if (isEventLoss==2 || isEventLoss==0)    IsParticleTrue=1;
    isEffMassSel=0;
  }
  //**********************************

  TString PathIn0;
  TString PathInAC;

  TString AllPathIn0[2][3];
  TString AllPathInAC[2][3];

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=10;//9         
  const Int_t numtipo=10;
  TString RegionTypeOld[3] = {"Jet", "Bulk", "All"};
  if (isBulkBlue) RegionTypeOld[1] = "BulkBlue";
  TString tipo[numtipo]={"K0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString tipoTitle[numtipo]={"K0s", "#Lambda", "#AntiLambda","#Lambda + #AntiLambda", "Xi^{-}", "Xi^{+}", "Omega^{-}", "Omega^{+}", "Xi"  , "Omega"};
  cout << "ok " << endl;
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString SSkipAssoc[2] = {"_AllAssoc", ""};

  TString titleX=  "p_{T} (GeV/c)";
  TString titleXPhi=  "#Delta #phi";
  TString titleYEff=  "#epsilon_{part}";
  TString titleY=  "1/#Delta#eta #Delta#phi 1/N_{trigg} dN/dp_{T}";
  TString titleYPhi=  "N/N_{trigg} per #Delta#eta";
  //  TString titleYRatio = "Ratio to MC truth";
  //  TString titleYRatio = "N/N_{Trigg} (MC reco) / N/N_{Trigg} (MC gen)";
  //  if (isEventLoss) titleYRatio = "#epsilon_{part}/#epsilon_{Trigg event}"; 
  TString titleYRatio = "#epsilon_{part} / #epsilon_{Trigg event}";
  if (isEventLoss==0) {
    if (isRecoNormCorr)  titleYRatio = "MC reco NFcorr / MC truth";
    else  titleYRatio = "MC reco / MC truth";
  }
  else if (isEventLoss==1) titleYRatio = "MC hybrid / MC truth";
  if (isEventLoss==2) titleYRatio = "MC reco / MC hybrid";
  TString title = "Multiplicity class ";

  Float_t LimInfProjRatio=0.8;
  Float_t LimSupProjRatio=1.2;
  if (type==0 && isEventLoss==1){
    LimInfProjRatio=0.9;
    LimSupProjRatio=1.1;
  }
  Int_t Color[3] ={628, 418, 600};
  Int_t ColorC[2] ={1,907};
  Int_t Colormult[21]={1, 401, 801, 628, 909, 881, 860, 868, 841, 418, 628, 909, 881, 867, 921, 401, 841, 862, 866, 865, 864};
  TLegend *legend= new TLegend(0.6, 0.6, 0.9, 0.9);
  TLegend *legendMult= new TLegend(0.65, 0.65, 0.9, 0.9);
  legendMult->SetHeader("Multiplicity class");
  TLegend *legendFit;
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  TString SmoltLegend[nummolt+1]={"0-5 %", "5-10 %", "10-30 %", "30-50 %", "50-100 %", "0-100 %"};
  TString Smoltpp5TeV[nummolt+1]={"0-10", "10-100", "100-100", "100-100", "100-100", "_all"};
  TString SmoltLegendpp5TeV[nummolt+1]={"0-10 %", "10-100 %", "100-100 %", "100-100 %", "100-100 %", "0-100 %"};
  TString SmoltppHM[nummolt+1]={"0-0.001", "0.001-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.1", "0-0.1"};
  TString SmoltLegendppHM[nummolt+1]={"0-0.001 %", "0.001-0.005 %", "0.005-0.01 %", "0.01-0.05 %", "0.05-0.1 %", "0-0.1 %"};

  if (ispp5TeV && MultBinning==3){
    for (Int_t m=0; m<nummolt+1; m++){
      Smolt[m] = Smoltpp5TeV[m];
      SmoltLegend[m] = SmoltLegendpp5TeV[m];
    }
  }
  if (isppHM){
    for (Int_t m=0; m<nummolt+1; m++){
      Smolt[m] = SmoltppHM[m];
      SmoltLegend[m] = SmoltLegendppHM[m];
    }
    if (MultBinning==1){
      Smolt[2] = "0-0.01";
      SmoltLegend[2] = "0-0.01 %";
    }
  }

  Double_t Nmolt[nummolt+1]={0,1,5,15,30,100};
  TString hhCorr[2]={"", "_hhCorr"};

  Float_t LimSup=0.1;
  Float_t LimInf=0;
  Float_t LimInfEff=0.8;
  Float_t LimSupEff=1;
  Float_t LimInfRatioSpectra[3] = {0.9, 0.9, 0.9};
  Float_t LimSupRatioSpectra[3] = {1.1, 1.1, 1.1};
  if (type==8){
    LimInfRatioSpectra[0] = 0.75;
    LimInfRatioSpectra[1] = 0.75;
    LimInfRatioSpectra[2] = 0.85;
    LimSupRatioSpectra[0] = 1.25;
    LimSupRatioSpectra[1] = 1.25;
    LimSupRatioSpectra[2] = 1.15;
  }
  else if (type==0){
    LimInfRatioSpectra[0] = 0.75;
    LimInfRatioSpectra[1] = 0.75;
    LimSupRatioSpectra[0] = 1.25;
    LimSupRatioSpectra[1] = 1.25;
  }

  TString  stringout = Dir+"/DATA"+year0+"/MCClosureCheck";
  if (isEventLoss==0) stringout += "_RecotoTrue";
  else if (isEventLoss==1) stringout += "_HybridtoTrue";
  else if (isEventLoss==2) stringout += "_RecoToHybrid";
  stringout += hhCorr[ishhCorr];
  if (isEventLoss==0)  stringout += yearMC + "_vs_"+ yearMCTruth;
  else if (isEventLoss==1) stringout += yearMCHyb + "_vs_"+ yearMCTruth;
  else stringout += yearMC + "_vs_"+ yearMCHyb;
  if (PtBinning>0) stringout +=Form("_PtBinning%i",PtBinning);
  if(type>=0){
    if (!ishhCorr)      stringout +="_"+tipo[type];
    stringout +=Srap[israp];
    //    stringout += "_SkipAssoc";
    //stringout +=SSkipAssoc[SkipAssoc];
  }
  stringout+=   Form("_PtMin%.1f", PtTrigMin);
  if (IsParticleTrue) stringout+= "_IsParticleTrue";
  if (isEffMassSel) stringout+= "_isEffMassSel";
  if (isMEFromHybrid) stringout+= "_IsMEFromHybrid";
  if (isMEFromK0s) stringout+= "_IsMEFromK0s";
  if (isMEFromCorrectCentrality) stringout+= "_IsMEFromCorrectCentrality";
  if (IsNotSigmaTrigger) stringout+= "_IsNotSigmaTrigger";
  if (isSidebands && isEventLoss!=1) stringout += "_Sidebands";
  if (isEtaEff) stringout+= "_IsEtaEff";
  if (isTrigEff && isEventLoss!=1) stringout+= "_IsTrigEff";
  if (MultBinning!=0) stringout += Form("_MultBinning%i", MultBinning);
  //  if (yearMCTruth.Index("Fio")!=-1)  stringout += "_isAllDeltaEta";
  //  stringout+="_isPrimaryTrigger";
  if (sys!=0) stringout+=Form("_sys%i", sys);
  if (sysPhi!=0) stringout+=Form("_sysPhi%i", sysPhi);
  if (isBulkBlue) stringout+="_isBulkBlue";
  if (isRecoNormCorr) stringout+="_isRecoNormCorr";
  //  stringout += "_Test";
  TString stringoutpdf=stringout; 
  //  stringout+="_thinptbins";
  stringout+= ".root";
  TString stringoutshortpdf = "ClosureTest_";
  if (ispp5TeV) stringoutshortpdf += "_pp5TeV";

  TCanvas *  canvasTriggerEfficiency = new TCanvas ("canvasTriggerEfficiency", "canvasTriggerEfficiency", 1300, 800);
  TCanvas *  canvasSpectrum[3];
  TCanvas *  canvasSpectrumRatio[3];
  TCanvas *  canvasPtEfficiency[3];
  TCanvas *  canvasSpectrumRelErrorClosure[3];
  TCanvas *  canvasFitResult[3];
  TCanvas *  canvasProj[nummolt+1][3];
  TCanvas *  canvasProjRatio[nummolt+1][3];
  Float_t   Sigma[3][nummolt+1][numPtV0]={0};
  TH1F *hSpectrum[3][2][nummolt+1];
  TH1F *hSpectrumPtEff[3][nummolt+1];
  TH1F *hSpectrumratio[3][nummolt+1];
  TH1F *hSpectrumRelErrorClosure[3][nummolt+1];
  TF1* pol0NormFactor[nummolt+1];
  TH1F*  hPhiFit[3][nummolt+1];
  TH1F*   hDeltaEtaDeltaPhi_PhiProj[3][2][nummolt+1][numPtV0];
  TH1F*   hDeltaEtaDeltaPhi_PhiProjRatio[3][nummolt+1][numPtV0];
  TF1* lineat1 = new TF1 ("pol0", "pol0", 0,8);
  lineat1->SetLineWidth(1);
  lineat1->SetLineColor(1);
  lineat1->FixParameter(0,1);
  TF1* lineat1Phi = new TF1 ("lineat1Phi", "pol0", -TMath::Pi()/2, 3./2*TMath::Pi());
  lineat1Phi->SetLineWidth(1);
  lineat1Phi->SetLineColor(1);
  lineat1Phi->FixParameter(0,1);
  TF1* PhiFit[3][nummolt+1][numPtV0];
  

  Int_t PtV0Min=1; //el 1                                                                 
  if (!ishhCorr && type==0) PtV0Min=0;
  Int_t numPtV0Max=numPtV0;
  if (PtBinning==1) numPtV0Max = numPtV0-1;
  //  if (PtBinning==1) numPtV0Max = numPtV0;
  else {
    if (type==0) numPtV0Max = numPtV0-2;
    else numPtV0Max = numPtV0-2;
    //numPtV0Max = numPtV0-2;
  }


  TString SPtV0[numPtV0]={"", "0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8", ""};
  //  TString SPtV0[numPtV0]={"", "","0.5-1", "1-1.5","1.5-2", "2-3", "3-4", "4-8"};                              
  if (type>0)SPtV0[1]={"0.5-1"};
  else {SPtV0[0]={"0-0.5"}; SPtV0[1]={"0.5-1"};}
  if (ishhCorr){
    SPtV0[0]={"0.1-0.5"};
    SPtV0[1]={"0.5-1"};
  }
  Double_t NPtV0[numPtV0+1]={0,0,1,1.5,2,2.5,3,4,8, 100};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5,1,1.5,2,3,4,8};                                                            
  NPtV0[1]=0.5;
  if (ishhCorr) {
    NPtV0[0]=0.1;
    NPtV0[1]=0.5;
  }

  TString SNPtV0[numPtV0+1]={"0.0","0.0","1.0","1.5","2.0","2.5","3.0","4.0","8.0", ""};
  SNPtV0[1]={"0.5"};
  if (ishhCorr){
    SNPtV0[0]={"0.1"};
    SNPtV0[1]={"0.5"};
  }

  TString SPtV01[numPtV0]={"0.1-0.5", "0.5-0.8", "0.8-1.2", "1.2-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV01[numPtV0+1]={0.1,0.5,0.8, 1.2,1.6,2,2.5,3,4,8};
  TString SPtV02[numPtV0]={"0.1-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.6","1.6-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV02[numPtV0+1]={0.1,0.4,0.6, 0.8,1.6,2,2.5,3,4,8};

  if (PtBinning==1){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV01[v];
      NPtV0[v] = NPtV01[v];
    }
  }
  if (PtBinning==2){
    for(Int_t v=PtV0Min; v<numPtV0Max+1; v++){
      if (v<numPtV0Max)      SPtV0[v] = SPtV02[v];
      NPtV0[v] = NPtV02[v];
    }
  }


  TFile *filein[2];
  TFile *fileinAC[2];
  TString nameSE[nummolt+1][numPtV0];
  TString namePhiProj[3][nummolt+1][numPtV0];
  //  TString namePhiProjReg[3]= {"_phi_V0Sub_BulkSub_EffCorr", "_phi_V0Sub_Bulk_EffCorr", "_phi_V0Sub_JetBulk_EffCorr"};
  TString namePhiProjReg[3]= {"_phi_V0Sub_EffCorr", "_phi_V0Sub_Bulk_EffCorr", "_phi_V0Sub_JetBulk_EffCorr"};

  TH1F*  fHistEventLoss;
  TH1F*  fHistEventLossMultAll;
  Float_t TriggerEff[nummolt+1]={0};

  Float_t YSup=0.01;
  Float_t YInf=-0.001;

  gStyle->SetOptStat(0);
  //  gStyle->SetOptFit(0);

  for (Int_t r =0; r<3; r++){
    cout << "\n\n\e[35m**************************************************"<< endl;
    cout << "Region: "<<RegionTypeOld[r] << "\e[39m"<< endl; 
    if ((sys!=0 || sysPhi!=0) && r==2) continue;
    if (sys==5 && r==0) continue;
    canvasSpectrumRatio[r] = new TCanvas( Form("canvasSpectrumRatio%i",r), 	"canvasSpectrumRatio_"+ RegionTypeOld[r], 1300, 800);
    canvasPtEfficiency[r] = new TCanvas( Form("canvasPtEfficiency%i",r), 	"canvasPtEfficiency_"+ RegionTypeOld[r], 1300, 800);
    canvasSpectrumRelErrorClosure[r] = new TCanvas( Form("canvasSpectrumRelErrorClosure%i",r), 	"canvasSpectrumRelErrorClosure_"+ RegionTypeOld[r], 1300, 800);
    canvasSpectrum[r] = new TCanvas( Form("canvasSpectrum%i",r), 	"canvasSpectrum_"+ RegionTypeOld[r], 1300, 800);
    if (ispp5TeV && MultBinning==3){
      canvasSpectrum[r]->Divide(2,2);
      canvasSpectrumRatio[r]->Divide(2,2);
      canvasSpectrumRelErrorClosure[r]->Divide(2,2);
    }
    else if (isppHM && MultBinning==1){
      canvasSpectrum[r]->Divide(2,2);
      canvasSpectrumRatio[r]->Divide(2,2);
      canvasSpectrumRelErrorClosure[r]->Divide(2,2);
    } else { 
      canvasSpectrum[r]->Divide(3,2);
      canvasSpectrumRatio[r]->Divide(3,2);
      canvasSpectrumRelErrorClosure[r]->Divide(3,2);
    }

    canvasFitResult[r] = new TCanvas( Form("canvasFitResult%i",r), 	"canvasFitResult_"+ RegionTypeOld[r], 1300, 800);
    for(Int_t m=0; m<nummolt+1; m++){
      if (MultBinning==3 && (m==2 || m==3 || m==4) ) continue;
      if (isppHM && MultBinning==1 && m<=1) continue;
      canvasProj[m][r]=new TCanvas(Form("canvasProj_m%i_MC%i",m,  r), "canvasProj"+Smolt[m], 1300, 800);
      canvasProj[m][r]->Divide(3,3);
      canvasProjRatio[m][r]=new TCanvas(Form("canvasProjRatio_m%i_MC%i",m,  r), "canvasProjRatio"+Smolt[m], 1300, 800);
      canvasProjRatio[m][r]->Divide(3,3);

    }

    for (Int_t i =0; i<2; i++){
      if (type==0){
	if (r==0) LimSup =0.025;
	else if (r==1 || r==2) LimSup = 0.1;
	LimInf =10e-5;
      }
      else  if (type==8) {
	if (r==0) LimSup =0.001;
	else   if (r==1) LimSup =0.01;
	else   if (r==2) LimSup =0.01;
	LimInf =10e-10;
      }
 
      PathIn0 = Dir+"/DATA"+year0+"/SystematicAnalysis";
      if (i==0){
	if(isEventLoss==1) PathIn0+= yearMCHyb;
	else PathIn0+= yearMC;
      }
      else if (i==1){
	if (isEventLoss!=2) PathIn0+= yearMCTruth;
	else PathIn0+=yearMCHyb;
      }
      if (PtBinning>0) PathIn0 +=Form("_PtBinning%i",PtBinning);
      if(type>=0){
	if (!ishhCorr)      PathIn0 +="_"+tipo[type];
	PathIn0 +=Srap[israp];
	if (!(i==0 && (isEventLoss==0 || isEventLoss==2)))	PathIn0 +="_AllAssoc";
	else if (SkipAssoc==0) PathIn0 +="_AllAssoc";
      }
      PathIn0+="_";
      PathIn0+= RegionTypeOld[r];
      if (i==0){
	if (isEventLoss==1) PathIn0+= "_MCTruth";
	else  PathIn0+= "MC";
      }
      else PathIn0+= "_MCTruth";
    
      PathIn0+=   Form("_PtMin%.1f", PtTrigMin);
      if (i==0 && (isEventLoss==0 || isEventLoss==2)){
	if (IsParticleTrue) PathIn0+= "_IsParticleTrue";
	if (isEffMassSel) PathIn0+= "_IsEfficiencyMassSel";
	if (isMEFromHybrid) PathIn0+= "_IsMEFromHybrid";
	if (isMEFromK0s) PathIn0+= "_IsMEFromK0s";
	if (isMEFromCorrectCentrality) PathIn0+= "_IsMEFromCorrectCentrality";
	if (isSidebands) PathIn0+= "_Sidebands";
	if (isEtaEff) PathIn0+= "_IsEtaEff";
	if (isTrigEff) PathIn0+= "_IsTrigEff";
      }
      if (MultBinning!=0) PathIn0 += Form("_MultBinning%i", MultBinning);
      //      if (yearMCTruth.Index("Fio")!=-1)  PathIn0 += "_isAllDeltaEta";
      //      if (PtBinning==0)      PathIn0+="_EL.root";
      //      else      PathIn0+=".root";
      if (i==1 && (isEventLoss==0 || isEventLoss==1)){
      if (IsNotSigmaTrigger) PathIn0+= "_IsNotSigmaTrigger";
      }
      if (i==1 && (isEventLoss==0 || isEventLoss==2)){
      if (isMEFromCorrectCentrality) PathIn0+= "_IsMEFromCorrectCentrality";
      }
      //      if (i==0)      PathIn0+="_thinptbins";
      //      if (i==0)      PathIn0+="_DCAz0.5";
      //      if (i==0)      PathIn0+="_isPrimaryTrigger";
      if (sys!=0) PathIn0+=Form("_sys%i", sys);
      if (sysPhi!=0) PathIn0+=Form("_sysPhi%i", sysPhi);
      if (isRecoNormCorr && i==0) PathIn0 += "_isNormCorr";
      if (isBkgParab) PathIn0+="_isBkgParab";
      PathIn0+= "_NewdEtaChoice";
      PathIn0+=".root";

      PathInAC = Dir+"/DATA"+year0+"/histo/AngularCorrelation";
      if (i==0){
	if(isEventLoss==1) PathInAC+= yearMCHyb;
	else PathInAC+= yearMC;
      }
      else if (i==1){
	if (isEventLoss!=2) PathInAC+= yearMCTruth;
	else PathInAC+=yearMCHyb;
      }

      if (i==0){
	if(isEventLoss==1) PathInAC+= "_MCTruth";
	else PathInAC+= "_MCEff";
      }
      else {
	PathInAC+= "_MCTruth";
      }
      if (PtBinning>0) PathInAC +=Form("_PtBinning%i",PtBinning);
      if(type>=0){
	if (!ishhCorr)      PathInAC +="_"+tipo[type];
	PathInAC +=Srap[israp];
       	if (!(i==0 && (isEventLoss==0 || isEventLoss==2)))	PathInAC +="_AllAssoc";
	else if (SkipAssoc==0) PathInAC +="_AllAssoc";
      }
      if (isBkgParab) PathInAC+="_isBkgParab";
      PathInAC+=   "_SysT0_SysV00";
      PathInAC+=   Form("_Sys%i_PtMin%.1f_Output", sys, PtTrigMin);
      if (i==0 && (isEventLoss==0 || isEventLoss==2)){
      if (IsParticleTrue) PathInAC+= "_IsParticleTrue";
      if (isEffMassSel) PathInAC+= "_IsEfficiencyMassSel";
      if (isMEFromHybrid) PathInAC+= "_IsMEFromHybrid";
      if (isMEFromK0s) PathInAC+= "_IsMEFromK0s";
      if (isMEFromCorrectCentrality) PathInAC+= "_IsMEFromCorrectCentrality";
      if (isSidebands) PathInAC+= "_Sidebands";
      if (isEtaEff) PathInAC+= "_IsEtaEff";
      if (isTrigEff) PathInAC+= "_IsTrigEff";
      }
      if (MultBinning!=0) PathInAC += Form("_MultBinning%i", MultBinning);
      //      if (yearMCTruth.Index("Fio")!=-1)  PathInAC += "_isAllDeltaEta";
      if (i==1 && (isEventLoss==0 || isEventLoss==1)){
      if (IsNotSigmaTrigger) PathInAC+= "_IsNotSigmaTrigger";
      }
      if (i==1 && (isEventLoss==0 || isEventLoss==2)){
      if (isMEFromCorrectCentrality) PathInAC+= "_IsMEFromCorrectCentrality";
      }
      //      if (i==0)      PathInAC+="_thinptbins";
      //      if (i==0)      PathInAC+="_DCAz0.5";
      //      if (i==0)      PathInAC+="_isPrimaryTrigger";
      PathInAC+= "_NewdEtaChoice";
      PathInAC+=".root";

      //      cout << "region " << r << " file " << i << " " << PathIn0 <<  " and " << PathInAC << endl;
      AllPathIn0[i][r] = PathIn0;
      AllPathInAC[i][r] = PathInAC;

      filein[i] = new TFile(PathIn0, "");
      cout << "Pathin " << PathIn0 << endl;
      if (!filein[i]) return;

      fileinAC[i] = new TFile(PathInAC, "");
      cout << "Pathin AC " << PathInAC << endl;
      if (!fileinAC[i]) return;

      //************* here I take the trigger efficiency calculated in BarlowSys.C when running Hybrid ***********
      if (isEventLoss==1 && i==0){
	fHistEventLoss = (TH1F*) filein[i]->Get("fHistEventLoss");
	fHistEventLossMultAll = (TH1F*) filein[i]->Get("fHistEventLossMultAll");
	if (!fHistEventLoss || !fHistEventLossMultAll) {cout << "event loss histogram not there!" << endl; return;}
	for(Int_t m=0; m<nummolt+1; m++){
	  if (MultBinning==3 && (m==2 || m==3 || m==4) ) continue;
	  if (isppHM && MultBinning==1 && m<=1) continue;
	  if (m!=nummolt) TriggerEff[m] = fHistEventLoss->GetBinContent(m+1);
	  else TriggerEff[m] = fHistEventLossMultAll->GetBinContent(1);
	}
	canvasTriggerEfficiency->cd();
	fHistEventLoss->Draw("");
      }
      //************* I have taken the trigger efficiency***********

      for(Int_t m=0; m<nummolt+1; m++){
	cout << "\n\e[32mMultiplicity "<< SmoltLegend[m] << "\e[39m" << endl;
	if (MultBinning==3 && (m==2 || m==3 || m==4) ) continue;
	if (isppHM && MultBinning==1 && m<=1) continue;
	if (m!=nummolt && isEventLoss==2 && type==8) continue;
	for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	  //	  cout << SPtV0[v] << " < pt < "<< SPtV0[v+1] << " GeV/c"<< endl;
	  nameSE[m][v]="ME_";
	  nameSE[m][v]+="m"+ Smolt[m]+"_v"+SPtV0[v]+"_AC";
	  namePhiProj[r][m][v]= nameSE[m][v] + namePhiProjReg[r];
	  //	  cout << 	  namePhiProj[r][m][v] << endl;

	  if(r==0)	  PhiFit[r][m][v] = new TF1 (Form("PhiFit_m%i_v%i_r%i",m,v,r), "pol0", -1.08,1.08);
	  //	  else if (r==1) 	  PhiFit[r][m][v] = new TF1 (Form("PhiFit_m%i_v%i_r%i",m,v,r), "pol0", -1.08, 1.08);
	  else  	  PhiFit[r][m][v] = new TF1 (Form("PhiFit_m%i_v%i_r%i",m,v,r), "pol0", -TMath::Pi()/2, 3./2*TMath::Pi());
	  PhiFit[r][m][v]->SetLineWidth(1);
	  PhiFit[r][m][v]->SetLineColor(kBlue);
	  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]= (TH1F*)fileinAC[i]->Get(namePhiProj[r][m][v]+ "_TrCorr");
	  if (!hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]){cout << "histo not there " << endl; return;}
	  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]->SetName(namePhiProj[r][m][v]+ Form("_TrCorr_%i", i));
	  //	  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]->Sumw2();
	  if (r!=0)	  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]->Rebin(2);
	  if (i==0)	  hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v] = (TH1F*) 	  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]->Clone(namePhiProj[r][m][v] + "_Ratio");
	  else   hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v]->Divide(  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]);
	}
	//	cout << "Getting the spectra... " << endl;
	hSpectrum[r][i][m] =(TH1F*) filein[i] ->Get(Form("fHistSpectrumPart_m%i_syst0", m));
	if (!hSpectrum[r][i][m]) {cout << " spectrum not there " << endl; return;}
	if(i==1)	hPhiFit[r][m]= (TH1F*) 	hSpectrum[r][i][m]->Clone(Form("FirResult_m%i_r%i", m, r));
 	hSpectrum[r][i][m]->Sumw2();
	if (i==0)      hSpectrumratio[r][m] = (TH1F*)      hSpectrum[r][i][m]->Clone("SpectrumRatio"+RegionTypeOld[r]+ Form("_m%i",m));
	else    hSpectrumratio[r][m]->Divide(   hSpectrum[r][i][m]);
	if (i==0)  {
	  hSpectrumRelErrorClosure[r][m] = (TH1F*)      hSpectrum[r][i][m]->Clone("SpectrumRelErrorClosure"+RegionTypeOld[r]+ Form("_m%i",m));
	  //hSpectrumRelErrorClosure[r][m]->Sumw2();
	}
	else {
	  hSpectrumRelErrorClosure[r][m]-> Add(hSpectrum[r][1][m],-1);
	  hSpectrumRelErrorClosure[r][m]-> Divide(hSpectrum[r][0][m]);
	}

	if (isEventLoss==1){
	  if (i==1){
	    for (Int_t b=1; b<=hSpectrumratio[r][m]->GetNbinsX(); b++){
	      //cout << "\n err uncorr" <<	    hSpectrumratio[r][m]->GetBinError(b)<<endl;
	      Sigma[r][m][b-1] =  hSpectrumratio[r][m]->GetBinContent(b)* sqrt(pow(hSpectrum[r][0][m]->GetBinError(b)/ hSpectrum[r][0][m]->GetBinContent(b),2)+ pow(hSpectrum[r][1][m]->GetBinError(b)/ hSpectrum[r][1][m]->GetBinContent(b),2)-1./(hSpectrum[r][0][m]->GetBinContent(b)*hSpectrum[r][1][m]->GetBinContent(b))*pow(hSpectrum[r][1][m]->GetBinError(b),2));
	      if (hSpectrumratio[r][m]->GetBinError(b)!=0){
		hSpectrumratio[r][m]->SetBinError(b, Sigma[r][m][b-1]);
	      }
	      //cout << " err corr" <<	    hSpectrumratio[r][m]->GetBinError(b)<<endl;
	    }
	    if (r==0 && type==8) {
	      pol0NormFactor[m] = new TF1(Form("pol0NormFactor_m%i", m), "pol0", 0, 8);
	      pol0NormFactor[m]->SetLineColor(634);
	      pol0NormFactor[m]->SetLineWidth(0.4);
	      cout << "pol0 fit to normalization factor (only for Xi in jet)" << endl;
	      hSpectrumratio[r][m]-> Fit(pol0NormFactor[m], "R+");
	      cout << "Fit to normalisation factor:\nConst = "<< pol0NormFactor[m]->GetParameter(0) << " +- " << pol0NormFactor[m]->GetParError(0) << endl;
	    }
	  }
       	}

	for(Int_t v=PtV0Min; v<numPtV0Max; v++){
	  canvasProj[m][r]->cd(v+1);
	  //	  if (m==0 && type==0) YSup = 0.014;
	  //	  else YSup=0.006;
	  YInf = -0.001;
	  YSup=0.006;
	  if (r==1 || r==2) {
	    YInf=0;
	    YSup =0.01;
	    if (NPtV0[v]>=2.5) YSup = 0.003;
	  }
	  if (type==8){
	    //	    YSup=0.0006;
	    YSup=0.0001*2;
	    //	    YInf = -0.0001;
	    YInf =0;
	    if (m==0) YSup = 0.0012;
	  }
	  StyleHisto(hDeltaEtaDeltaPhi_PhiProj[r][i][m][v], YInf, YSup, ColorC[i], 1, titleXPhi, titleYPhi,  Form("%.1f < p_{T} < %.1f",NPtV0[v], NPtV0[v+1]),  0, 0, 0);
	  hDeltaEtaDeltaPhi_PhiProj[r][i][m][v]->Draw("same");
	  //	  if(i==0)legend->Draw("");

	  canvasProjRatio[m][r]->cd(v+1);
	  gPad->SetLeftMargin(0.15);
	  StyleHisto(hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v], LimInfProjRatio, LimSupProjRatio, Color[r], 1, titleXPhi, titleYRatio, SmoltLegend[m] + Form(",  %.1f < p_{T} < %.1f GeV/c",NPtV0[v], NPtV0[v+1]) ,  0, 0, 0);

	  hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v]->GetYaxis()->SetTitleSize(0.04);
	  hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v]->Fit(PhiFit[r][m][v], "RQ+");
	  legendFit	  = new TLegend(0.4, 0.7, 0.9, 0.9);
	  //	  legendFit->AddEntry(PhiFit[r][m][v], Form("m=%.3f +- %.3f", PhiFit[r][m][v]->GetParameter(1), PhiFit[r][m][v]->GetParError(1)), "pl"); 
	  legendFit->AddEntry(PhiFit[r][m][v], Form("q=%.3f +- %.3f", PhiFit[r][m][v]->GetParameter(0), PhiFit[r][m][v]->GetParError(0)), "pl"); 
	  if (r==0) hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v]->GetXaxis()->SetRangeUser(-1.2, 1.2);
	  if (i==1)	{
	    hDeltaEtaDeltaPhi_PhiProjRatio[r][m][v]->Draw("same");
	    lineat1Phi->Draw("same");
	    legendFit->Draw("");
	  }
	  
	  if (i==1){
	  hPhiFit[r][m]->SetBinContent(	  hPhiFit[r][m]->FindBin(NPtV0[v]+0.0001), PhiFit[r][m][v]->GetParameter(0));
	  hPhiFit[r][m]->SetBinError(	  hPhiFit[r][m]->FindBin(NPtV0[v]+0.0001), PhiFit[r][m][v]->GetParError(0));
	  if(r==0 && NPtV0[v]==0)	  hPhiFit[r][m]->SetBinError(	  hPhiFit[r][m]->FindBin(NPtV0[v]+0.0001),0);
	  if(r==0 && NPtV0[v]==0)	  hPhiFit[r][m]->SetBinContent(	  hPhiFit[r][m]->FindBin(NPtV0[v]+0.0001),0);
	  }

	}
	if (i==1){
	canvasFitResult[r]->cd();
	gStyle->SetOptStat(0);
	//	gStyle->SetOptFit(0);
	hPhiFit[r][m]->SetMarkerColor(Colormult[m]);
	if(isEventLoss==2)	hPhiFit[r][m]->GetYaxis()->SetTitle(" pol0 parameter MC reco/MC gen");
	hPhiFit[r][m]->SetTitle("");
	hPhiFit[r][m]->GetYaxis()->SetTitleOffset(1.2);
	hPhiFit[r][m]->SetLineColor(Colormult[m]);
	hPhiFit[r][m]->GetYaxis()->SetRangeUser(0.85, 1.15);
	if (r==0) legendMult->AddEntry(hPhiFit[r][m], SmoltLegend[m], "pl");
	//	if (m==nummolt)	hPhiFit[r][m]->Draw("same");
	hPhiFit[r][m]->Draw("same");
	if (m==nummolt)	legendMult->Draw("");
	if (m==nummolt) lineat1->Draw("same");
	}
	if (isppHM && MultBinning==1) 	canvasSpectrum[r]->cd(m-1);
	else if (ispp5TeV && MultBinning==3 && m==nummolt) 	canvasSpectrum[r]->cd(3);
	else 	canvasSpectrum[r]->cd(m+1);
	gPad->SetLeftMargin(0.15);
	StyleHisto(hSpectrum[r][i][m], 0, LimSup, ColorC[i], 1, titleX, titleY,  title+SmoltLegend[m],  0, 0, 0);
	if(m==0 && r==0){
	  if(i==0){
	    if (isEventLoss!=1)	legend->AddEntry(hSpectrum[r][i][m],"MC reco","pel");
	    else legend->AddEntry(hSpectrum[r][i][m],"MC Hybrid","pel");
	  }
	  else  {
	    if (isEventLoss!=2)  legend->AddEntry(hSpectrum[r][i][m],"MC Truth","pel");
	    else   legend->AddEntry(hSpectrum[r][i][m],"MC Hybrid","pel");
	  }
	}
	hSpectrum[r][i][m]->Draw("same");
	if(i==0)legend->Draw("");

	if (ispp5TeV && MultBinning==3 && m==nummolt) 	canvasSpectrumRatio[r]->cd(3);
	else if (isppHM && MultBinning==1) 	canvasSpectrumRatio[r]->cd(m-1);
	else 	canvasSpectrumRatio[r]->cd(m+1);
	gPad->SetLeftMargin(0.15);
	gStyle->SetOptStat(0);
	StyleHisto(hSpectrumratio[r][m], LimInfRatioSpectra[r], LimSupRatioSpectra[r], Color[r], 1, titleX, titleYRatio,  title+SmoltLegend[m],  0, 0, 0);
	hSpectrumratio[r][m]->GetYaxis()->SetTitleSize(0.06);
	hSpectrumratio[r][m]->GetYaxis()->SetTitleOffset(1.2);
	if (i==1)	{
	  hSpectrumratio[r][m]->Draw("same");
	  lineat1->Draw("same");
	}

	if (isEventLoss==1 && i==1){
	  canvasPtEfficiency[r]->cd();
	  gPad->SetLeftMargin(0.15);
	  gStyle->SetOptStat(0);
	  hSpectrumPtEff[r][m] = (TH1F*) hSpectrumratio[r][m]->Clone("SpectrumEff"+RegionTypeOld[r]+ Form("_m%i",m));
	  cout << "\nTrigger eff " <<	    TriggerEff[m] << endl;
	  for (Int_t b=1; b<=  hSpectrumPtEff[r][m]->GetNbinsX(); b++){
	    cout <<hSpectrumPtEff[r][m]->GetXaxis()->GetBinLowEdge(b) << " < pt < "<< hSpectrumPtEff[r][m]->GetXaxis()->GetBinUpEdge(b) << " Norm factor: "<<   hSpectrumratio[r][m]->GetBinContent(b)<< endl;
	  }
	  hSpectrumPtEff[r][m]->Scale(TriggerEff[m]);
	  StyleHisto(hSpectrumPtEff[r][m], LimInfEff, LimSupEff, Colormult[m], 33, titleX, titleYEff, ""/* title+SmoltLegend[m]*/,  0, 0, 0);
	  hSpectrumPtEff[r][m]->GetYaxis()->SetTitleSize(0.04);
	  hSpectrumPtEff[r][m]->GetYaxis()->SetTitleOffset(1.2);
	  cout << "\n"<< endl;
	  for (Int_t b=1; b<=  hSpectrumPtEff[r][m]->GetNbinsX(); b++){
	    cout <<hSpectrumPtEff[r][m]->GetXaxis()->GetBinLowEdge(b) << " < pt < "<< hSpectrumPtEff[r][m]->GetXaxis()->GetBinUpEdge(b) << " Epsilon_part: "<<   hSpectrumPtEff[r][m]->GetBinContent(b)<< endl;
	  }
	  hSpectrumPtEff[r][m]->Draw("same ep");
	  lineat1->Draw("same");
	}

	if (ispp5TeV && MultBinning==3 && m==nummolt) 	canvasSpectrumRelErrorClosure[r]->cd(3);
	else if (isppHM && MultBinning==1) 	canvasSpectrumRelErrorClosure[r]->cd(m-1);
	else 	canvasSpectrumRelErrorClosure[r]->cd(m+1);

	gPad->SetLeftMargin(0.15);
	gStyle->SetOptStat(0);
	StyleHisto(hSpectrumRelErrorClosure[r][m], -0.3, 0.3, Color[r], 1, titleX, "uncertainty closure",  title+SmoltLegend[m],  0, 0, 0);
	hSpectrumRelErrorClosure[r][m]->GetYaxis()->SetTitleSize(0.08);
	hSpectrumRelErrorClosure[r][m]->GetYaxis()->SetTitleOffset(1.2);
	if (i==1)	{
	  hSpectrumRelErrorClosure[r][m]->Draw("same");
	  lineat1->Draw("same");
	}
      }
    }//end loop on file
  }//end loop on regiontype

  cout << "\n\nSaving canvases"<< endl;
  TFile * fileout = new TFile(stringout, "RECREATE");
  fileout->WriteTObject(fHistEventLossMultAll);
  fileout->WriteTObject(fHistEventLoss);
  if (isTrigEff) {
    canvasTriggerEfficiency->SaveAs(stringoutpdf+"_TriggerEfficiency.pdf");
    fileout->WriteTObject(canvasTriggerEfficiency);
  }
  cout << "ciao " << endl;
  for (Int_t r =0; r<3; r++){ 
    if ((sys!=0 || sysPhi!=0) && r==2) continue;
    if (sys==5 && r==0) continue;
    canvasFitResult[r]->SaveAs(stringoutpdf+"_fitresult"+RegionTypeOld[r]+".pdf");
    fileout->WriteTObject(canvasSpectrum[r]);
    fileout->WriteTObject(canvasSpectrumRatio[r]);
    if (isEventLoss==1) {
      fileout->WriteTObject(canvasPtEfficiency[r]);
      canvasPtEfficiency[r]->SaveAs(stringoutpdf+"_PtEff_"+RegionTypeOld[r]+".pdf");
    }
    fileout->WriteTObject(canvasSpectrumRelErrorClosure[r]);
    fileout->WriteTObject(canvasFitResult[r]);
    canvasSpectrumRatio[r]->SaveAs(stringoutpdf+"_PtRatio_"+RegionTypeOld[r]+".pdf");
    for(Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4) ) continue;
      if (m!=nummolt && isEventLoss==2 && type==8) continue;
      fileout->WriteTObject(      hSpectrumratio[r][m]);
      fileout->WriteTObject(      hSpectrumRelErrorClosure[r][m]);
      canvasProjRatio[m][r]->SaveAs(stringoutpdf+"_PhiRatio_"+Smolt[m] + "_"+RegionTypeOld[r]+".pdf");
      fileout->WriteTObject(canvasProj[m][r]);
      fileout->WriteTObject(canvasProjRatio[m][r]);
      fileout->WriteTObject(hPhiFit[r][m]);
    }
  }
  
  cout << "Saving in pdf file: " << endl;
  Int_t counterm=-1;
  Int_t counterr=-1;
  for (Int_t r =0; r<3; r++){ 
    if ((sys!=0 || sysPhi!=0) && r==2) continue;
    if (sys==5 && r==0) continue;
    counterr++;
    for(Int_t m=0; m<nummolt+1; m++){
      if (isppHM && MultBinning==1 && m<=1) continue;
      if (MultBinning==3 && (m==2 || m==3 || m==4) ) continue;
      if (m!=nummolt && isEventLoss==2 && type==8) continue;
      counterm++;
      canvasProjRatio[m][r]->SaveAs(stringoutshortpdf + Form("PhiRatio_m%i_r%i.pdf", m, r));
      if (counterm==0 && counterr==0) canvasProjRatio[m][r]->SaveAs(stringoutpdf+".pdf(");
      else canvasProjRatio[m][r]->SaveAs(stringoutpdf+".pdf");
      canvasProj[m][r]->SaveAs(stringoutpdf+".pdf");
    }
    canvasFitResult[r]->SaveAs(stringoutpdf+".pdf");
    canvasFitResult[r]->SaveAs(stringoutshortpdf + Form("FitResult_r%i.pdf", r));
    if (r==2)    canvasSpectrumRatio[r]->SaveAs(stringoutpdf+".pdf)");
    else     canvasSpectrumRatio[r]->SaveAs(stringoutpdf+".pdf");
    canvasSpectrumRatio[r]->SaveAs(stringoutshortpdf + Form("SpectrumRatio_r%i.pdf", r));
  }

  for (Int_t i =0; i<2; i++){ 
    cout << "\n" << endl;
    for (Int_t r =0; r<3; r++){ 
      if ((sys!=0 || sysPhi!=0) && r==2) continue;
      if (sys==5 && r==0) continue;
      cout << " region " << r << endl;
      cout << " partendo da " <<       AllPathInAC[i][r] << " e " <<       AllPathIn0[i][r] << endl;
    }
  }
  cout << "\n\e[36mHo creato il file:\e[39m " << stringout << "\n (anche il pdf con parte delle canvas)\n"<< endl;
  fileout->Close();
}

