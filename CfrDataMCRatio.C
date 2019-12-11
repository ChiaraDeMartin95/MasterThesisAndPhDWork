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

void CfrDataMCRatio(Float_t PtTrigMin=3.0, TString year0="2016", Int_t ishhCorr=1, Bool_t MasterThesis=0, Bool_t isSE=1){
  //if isSE is true the ratio of K0s to hadrons is computed, in order to study the strangeness enhancement effect


  Int_t WhatKindOf=0;

  if (!ishhCorr && !isSE) WhatKindOf=0;
  if (ishhCorr) WhatKindOf=1;
  if (!ishhCorr && isSE) WhatKindOf=2;

  if (isSE) ishhCorr=0;
  TString PathIn[2];
  TString PathOut;
  TFile *filein[2];
 
  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  const Int_t numPeriod=2;

  //colori e marker diversi per jet, OJ, J+OJ
  Int_t Color[3]={2,418, 4}; 
  Int_t MarkerOrigin[3]={20,21,33}; 
  Int_t MarkerMC[3]={24,25,27}; 

  TString Sys[6]={"Jet DATA", "Jet MC", "OJ DATA", "OJ MC", "J+OJ DATA", "J+OJ MC"};
  TString SysErr[2]={"syst. DATA", "syst. MC"};
  TString StatErr[2]={"stat. DATA", "stat. MC"};
  TString StatErrBis[3]={"stat. in-jet", "stat. out-of-jet", "stat. inclusive"};  
  TString SysErrBis[3]={"sist. in-jet", "sist. out-of-jet", "sist. inclusive"};
  TString SysErrMB[2]={"syst. DATA 0-100 %", "syst. MC 0-100 %"};
  TString StatErrMB[2]={"stat. DATA 0-100 %", "stat. MC 0-100 %"};

  TString Region[3]={"Jet", "Bulk", "All"};
  TString DATAorMC[2]={"DATA", "MC"};
  TString JetOrNot[3]={"Jet ", "OJ " , "J+OJ "};

  TString tipo[numtipo]={"kK0s", "bo"};
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  TString SmoltBis[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  TString SPtV0[numPtV0]={"0-1", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString SPeriod[numPeriod]={"2016", "2017"};//, "2018"};
  TString SRun[numPeriod]={"2016k", "2017h"};//, "2018m"};

  //20,21,33
  Int_t Marker[6]={33,4,33,4, 33, 4};

  Int_t MarkerBis[2]={23,24};
  Int_t MarkerTris[2]={31, 33};
  Int_t MarkerMBStat[6]={20,4,20,4, 20, 4};
  Int_t MarkerMBSist[6]={21,25,21,25, 21, 25};
  TF1* pol0[nummolt+1][numPtV0][numPeriod];
  TF1* pol0bis[nummolt+1][numPtV0][numPeriod];

  TH1D*  fHistYieldvsErrSoloStat[6][6];
  TH1D*  fHistYieldvsErrSoloSist[6][6];
  TH1D*  fHistYieldvsErrSoloStatMB[6][6];
  TH1D*  fHistYieldvsErrSoloSistMB[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatio[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatio[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatioMB[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatioMB[6][6];

  TH1D*  fHistYieldvsErrErrSoloStat[6][6];
  TH1D*  fHistYieldvsErrErrSoloSist[6][6];
  TH1D*  fHistYieldvsErrErrSoloStatMB[6][6];
  TH1D*  fHistYieldvsErrErrSoloSistMB[6][6];

  TH1D*  fHistYieldvsErrSoloStat_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSist_SE[6][6];
  TH1D*  fHistYieldvsErrSoloStatMB_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSistMB_SE[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatio_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatio_SE[6][6];
  TH1D*  fHistYieldvsErrSoloStatRatioMB_SE[6][6];
  TH1D*  fHistYieldvsErrSoloSistRatioMB_SE[6][6];

  TH1D*  fHistYieldvsErrErrSoloStat_SE[6][6];
  TH1D*  fHistYieldvsErrErrSoloSist_SE[6][6];
  TH1D*  fHistYieldvsErrErrSoloStatMB_SE[6][6];
  TH1D*  fHistYieldvsErrErrSoloSistMB_SE[6][6];

  TH1D*    fHistYieldDatiPubblicati;
  TH1D*    fHistYieldErroriDatiPubblicati;

  TFile *filedatipubbl;
  TDirectoryFile *dir;
  if (!ishhCorr){
  filedatipubbl= new TFile("HEPData-1574358449-v1-Table_8c.root", "");
  dir  = (TDirectoryFile*)filedatipubbl->Get("Table 8c");
  fHistYieldDatiPubblicati=(TH1D*)dir->Get("Hist1D_y1");
  fHistYieldErroriDatiPubblicati=(TH1D*)dir->Get("Hist1D_y1_e3");
  }

  TLegend* legend[3];
  TLegend* legendtype[3];
  TLegend* legendbis[3];
  Int_t msize=2;

  auto legenderr=new TLegend(0.6,0.6, 0.9, 0.9);
  auto legenderrbis=new TLegend(0.6,0.6, 0.9, 0.9);

  TString  PathInJet[2];
  TString  PathInJetMC[2];
  TString  PathInBulk[2];
  TString  PathInBulkMC[2];
  TString  PathInAll[2];
  TString  PathInAllMC[2];
  TString  hhCorr[3]={"", "_hhCorr", "_SE"};

  for(Int_t ishhCorrLoop=0; ishhCorrLoop<2; ishhCorrLoop++){
    if (!MasterThesis){
      PathInJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_Jet_PtMin%.1f.root", PtTrigMin);	 
      PathInJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_JetMC_PtMin%.1f.root", PtTrigMin);	 
      PathInBulk[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_Bulk_PtMin%.1f.root", PtTrigMin); 	 
      PathInBulkMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_BulkMC_PtMin%.1f.root", PtTrigMin);
      PathInAll[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_All_PtMin%.1f.root", PtTrigMin); 	 
      PathInAllMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorrLoop] +Form("_AllMC_PtMin%.1f.root", PtTrigMin);  
    }
    //per analsisi presentata nella tesi
    else{
      if (ishhCorrLoop=1) continue;
      PathInJet[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_Jet.root";	 
      PathInJetMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_JetMC.root";	 
      PathInBulk[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_Bulk.root"; 	 
      PathInBulkMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_BulkMC.root";
      PathInAll[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_All.root"; 	 
      PathInAllMC[ishhCorrLoop]="FinalOutput/DATA" + year0 + "/SystematicAnalysis" + hhCorr[ishhCorr] +"_AllMC.root";  
    }
  }

  //**************** prendo gli istogrammi **********************************************
  for(Int_t l=0; l<3; l++){ //loop sulla region (jet, OJ, inclusive)
    for(Int_t j=0; j<2; j++){ //loop data or MC
      for (Int_t ishhCorrLoop=0; ishhCorrLoop<2; ishhCorrLoop++){
	if(l==0){
	  if (j==0)   PathIn[ishhCorrLoop]=PathInJet[ishhCorrLoop];
	  if (j==1)   PathIn[ishhCorrLoop]=PathInJetMC[ishhCorrLoop];
	}
	if(l==1){
	  if (j==0)  PathIn[ishhCorrLoop]=PathInBulk[ishhCorrLoop];
	  if (j==1)  PathIn[ishhCorrLoop]=PathInBulkMC[ishhCorrLoop];
	}
	if (l==2){
	  if (j==0)  PathIn[ishhCorrLoop]=PathInAll[ishhCorrLoop];
	  if (j==1)  PathIn[ishhCorrLoop]=PathInAllMC[ishhCorrLoop];
	}
	filein[ishhCorrLoop]=new TFile(PathIn[ishhCorrLoop], "");
      }
      fHistYieldvsErrSoloSistMB[j][l]=   (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloSistMB");   
      fHistYieldvsErrSoloStatMB[j][l]=   (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloStatMB");   
      fHistYieldvsErrSoloStat[j][l]=     (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloStat");	    
      fHistYieldvsErrSoloSist[j][l]=     (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloSist");	    
      fHistYieldvsErrSoloStatRatio[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloStatRatio");
      fHistYieldvsErrSoloSistRatio[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrSoloSistRatio");
      fHistYieldvsErrErrSoloStat[j][l]=  (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloStat");  
      fHistYieldvsErrErrSoloSist[j][l]=  (TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloSist");  
      fHistYieldvsErrErrSoloStatMB[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloStatMB");
      fHistYieldvsErrErrSoloSistMB[j][l]=(TH1D*)filein[ishhCorr]->Get("fHistYieldvsErrErrSoloSistMB");

      if (isSE){
	cout << "ishhCorr vale " << ishhCorr << endl;
	cout << "l (region ) " << l << " data or MC j " << j << " path hh" << PathIn[1]<< endl;   
      fHistYieldvsErrSoloSistMB_SE[j][l]=   (TH1D*)filein[1]->Get("fHistYieldvsErrSoloSistMB");   
      fHistYieldvsErrSoloStatMB_SE[j][l]=   (TH1D*)filein[1]->Get("fHistYieldvsErrSoloStatMB");   
      fHistYieldvsErrSoloStat_SE[j][l]=     (TH1D*)filein[1]->Get("fHistYieldvsErrSoloStat");	    
      fHistYieldvsErrSoloSist_SE[j][l]=     (TH1D*)filein[1]->Get("fHistYieldvsErrSoloSist");	    
      fHistYieldvsErrSoloStatRatio_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrSoloStatRatio");
      fHistYieldvsErrSoloSistRatio_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrSoloSistRatio");
      fHistYieldvsErrErrSoloStat_SE[j][l]=  (TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloStat");  
      fHistYieldvsErrErrSoloSist_SE[j][l]=  (TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloSist");  
      fHistYieldvsErrErrSoloStatMB_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloStatMB");
      fHistYieldvsErrErrSoloSistMB_SE[j][l]=(TH1D*)filein[1]->Get("fHistYieldvsErrErrSoloSistMB");

      // fHistYieldvsErrSoloSistMB_SE[j][l]     ->Sumw2();
      // fHistYieldvsErrSoloStatMB_SE[j][l]   ->Sumw2();
      // fHistYieldvsErrSoloStat_SE[j][l]     ->Sumw2();
      // fHistYieldvsErrSoloSist_SE[j][l]  	  ->Sumw2();
      // fHistYieldvsErrSoloStatRatio_SE[j][l]->Sumw2();
      // fHistYieldvsErrSoloSistRatio_SE[j][l]->Sumw2();
      // fHistYieldvsErrErrSoloStat_SE[j][l]  ->Sumw2();
      // fHistYieldvsErrErrSoloSist_SE[j][l]  ->Sumw2();
      // fHistYieldvsErrErrSoloStatMB_SE[j][l]->Sumw2();
      //fHistYieldvsErrErrSoloSistMB_SE[j][l]->Sumw2();

      fHistYieldvsErrSoloSistMB[j][l]   ->Divide(fHistYieldvsErrSoloSistMB_SE[j][l]    );
      fHistYieldvsErrSoloStatMB[j][l]   ->Divide(fHistYieldvsErrSoloStatMB_SE[j][l]    );
      fHistYieldvsErrSoloStat[j][l]     ->Divide(fHistYieldvsErrSoloStat_SE[j][l]      );
      fHistYieldvsErrSoloSist[j][l]     ->Divide(fHistYieldvsErrSoloSist_SE[j][l]      );
      fHistYieldvsErrSoloStatRatio[j][l]->Divide(fHistYieldvsErrSoloStatRatio_SE[j][l] );
      fHistYieldvsErrSoloSistRatio[j][l]->Divide(fHistYieldvsErrSoloSistRatio_SE[j][l] );
      fHistYieldvsErrErrSoloStat[j][l]  ->Divide(fHistYieldvsErrErrSoloStat_SE[j][l]   );
      fHistYieldvsErrErrSoloSist[j][l]  ->Divide(fHistYieldvsErrErrSoloSist_SE[j][l]   );
      fHistYieldvsErrErrSoloStatMB[j][l]->Divide(fHistYieldvsErrErrSoloStatMB_SE[j][l] );
      fHistYieldvsErrErrSoloSistMB[j][l]->Divide(fHistYieldvsErrErrSoloSistMB_SE[j][l] );

      }
    }
  }


  //********cfr data MC separatamente per J, OJ, J+OJ*********************************************************

  Int_t col=0;
  TCanvas *canvas[3];
  for(Int_t l=0; l<3; l++){ //loop sulla region (jet, OJ, inclusive)
    Int_t count=0;
    canvas[l] = new TCanvas(Form("canvas%i", l), Region[l], 1300, 1000);
    legend[l]=new TLegend(0.7,0.1, 0.9, 0.3);
    legendtype[l]=new TLegend(0.7,0.35, 0.9, 0.45);
    legendbis[l]=new TLegend(0.1,0.7, 0.3, 0.9);

    for(Int_t j=0; j<2; j++){ //loop data or MC
      if (!ishhCorr)      fHistYieldvsErrSoloStat[j][l]->SetTitle ("Yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloStat[j][l]->SetTitle ("Yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloStat[j][l]->SetTitle ("Ratio of yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (ishhCorr)  fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetTitle ("N_{charged}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (isSE)  fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{charged}");      

      fHistYieldvsErrSoloStat[j][l]->SetMarkerColor(Color[l]);
      fHistYieldvsErrSoloStat[j][l]->SetMarkerStyle(MarkerOrigin[l]);
      if (j==1)      fHistYieldvsErrSoloStat[j][l]->SetMarkerStyle(MarkerMC[l]);
      fHistYieldvsErrSoloStat[j][l]->SetLineColor(Color[l]);
      fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0,0.4);
      if (l==0 && !ishhCorr) fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0,0.12);
      if (l==0 && ishhCorr) fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0,0.8);
      if (isSE) fHistYieldvsErrSoloStat[j][l]->GetYaxis()->SetRangeUser(0,0.1);
      fHistYieldvsErrSoloStat[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloStat[j][l]->SetMarkerSize(msize);
      if (l==2)    fHistYieldvsErrSoloStat[j][l]->SetMarkerSize(3);
      fHistYieldvsErrSoloStat[j][l]->Draw("samee");
    
      if (!ishhCorr)      fHistYieldvsErrSoloSist[j][l]->SetTitle ("Yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloSist[j][l]->SetTitle ("Yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloSist[j][l]->SetTitle ("Ratio of yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (ishhCorr)  fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetTitle ("N_{charged}/N_{Trigg} 1/#Delta#eta #Delta#phi");      
      if (isSE)  fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetTitle ("N_{K^{0}_{S}}/N_{charged}");      
      fHistYieldvsErrSoloSist[j][l]->SetMarkerColor(Color[l]);
      fHistYieldvsErrSoloSist[j][l]->SetMarkerStyle(MarkerOrigin[l]);
      if (j==1)       fHistYieldvsErrSoloSist[j][l]->SetMarkerStyle(MarkerMC[l]);
      fHistYieldvsErrSoloSist[j][l]->SetLineColor(Color[l]);
      fHistYieldvsErrSoloSist[j][l]->SetFillStyle(0);
      fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0,0.4);
      if (l==0 && !ishhCorr) fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0,0.12);
      if (l==0 && ishhCorr) fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0,0.8);
      if (isSE) fHistYieldvsErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0,0.1);
      fHistYieldvsErrSoloSist[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloSist[j][l]->SetMarkerSize(msize);
      if (l==2)    fHistYieldvsErrSoloSist[j][l]->SetMarkerSize(3);
      fHistYieldvsErrSoloSist[j][l]->Draw("samee2");
      
      legend[l]->AddEntry(fHistYieldvsErrSoloStat[j][l], JetOrNot[l]+DATAorMC[j], "pl");
      if(j==1)    legendtype[l]->AddEntry(fHistYieldvsErrSoloStatMB[j][l], "stat. ", "el");
      if(j==1)    legendtype[l]->AddEntry(fHistYieldvsErrSoloSistMB[j][l], "syst. ", "f");
      if(j==1)    legend[l]->Draw();
      if(j==1)    legendtype[l]->Draw();
    }
    //    canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/Yield"+Region[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
    //  if (MasterThesis)    canvas[l]->SaveAs("FinalOutput/DATA"+year0+"/Yield"+Region[l]+".pdf");
  }


  //cfr J, OJ, J+OJ separatamente per data e MC****************************************************
  TCanvas *canvasratiobis[2];
  for(Int_t l=0; l<2; l++){ //loop su data e MC
    Int_t count=0;
    canvasratiobis[l] = new TCanvas(Form("canvasratiobis%i", l), "Yieldratio_"+DATAorMC[l], 1300, 1000);
      for(Int_t j=0; j<3; j++){ //loop J, OJ, J+OJ

      if (!ishhCorr)      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle ("Relative yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle ("Relative yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle ("Ratio of relative yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetTitle ("(N_{K^{0}_{S}}/N_{Trigg}) / (N_{K^{0}_{S}}/N_{Trigg})_{0-100 %}");
      if (ishhCorr)  fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetTitle ("(N_{charged}/N_{Trigg}) / (N_{charged}/N_{Trigg})_{0-100 %}");
      if (isSE)  fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetTitle ("");      
      //      fHistYieldvsErrSoloStatRatio[l][j]->SetTitle("");
      fHistYieldvsErrSoloStatRatio[l][j]->GetYaxis()->SetRangeUser(0, 2.2); //..
      fHistYieldvsErrSoloStatRatio[l][j]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerSize(msize);
      if (j==2)    fHistYieldvsErrSoloStatRatio[l][j]->SetMarkerSize(3);
      fHistYieldvsErrSoloStatRatio[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrSoloStatRatio[l][j]->Draw("samee");
 
      //      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle("");
      if (!ishhCorr)      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle ("Relative yield of K^{0}_{S} in jet per trigger particle vs V0M multiplicity");
      if (ishhCorr)      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle ("Relative yield of charged unidentified in jet per trigger particle vs V0M multiplicity");
      if (isSE)      fHistYieldvsErrSoloSistRatio[l][j]->SetTitle ("Ratio of relative yields in jet per trigger particle vs V0M multiplicity");
      if (!ishhCorr)  fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetTitle ("(N_{K^{0}_{S}}/N_{Trigg}) / (N_{K^{0}_{S}}/N_{Trigg})_{0-100 %}");
      if (ishhCorr)  fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetTitle ("(N_{charged}/N_{Trigg}) / (N_{charged}/N_{Trigg})_{0-100 %}");
      if (isSE)  fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetTitle ("");      
      fHistYieldvsErrSoloSistRatio[l][j]->GetYaxis()->SetRangeUser(0, 2.2);//..
      fHistYieldvsErrSoloSistRatio[l][j]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      fHistYieldvsErrSoloSistRatio[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrSoloSistRatio[l][j]->SetFillStyle(0);
      fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerSize(msize);
      if (j==2)    fHistYieldvsErrSoloSistRatio[l][j]->SetMarkerSize(3);
      fHistYieldvsErrSoloSistRatio[l][j]->Draw("samee2");
        
      legendbis[l]->AddEntry(   fHistYieldvsErrSoloStatRatio[l][j], JetOrNot[j]+DATAorMC[l], "pl");
      if(j==0)legendbis[l]->Draw();
      if(j==0)legendtype[l]->Draw();
      // fHistYieldvsErrSoloSist[l][j]->SetMarkerColor(Color[j]);
      // fHistYieldvsErrSoloSist[l][j]->Draw("samee");
    }
    //    canvasratiobis[l]->SaveAs("FinalOutput/DATA2016/YieldRatio"+DATAorMC[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
      //  if (MasterThesis)        canvasratiobis[l]->SaveAs("FinalOutput/DATA2016/YieldRatio"+DATAorMC[l]+".pdf");
  }

  TCanvas *canvasnotratio[2];
  for(Int_t l=0; l<2; l++){
    Int_t count=0;
    canvasnotratio[l] = new TCanvas(Form("canvasnotratio%i", l), "Yield_"+DATAorMC[l], 1300, 1000);
    cout << l << endl;
    if (!ishhCorr){
    for(Int_t k=1; k < fHistYieldDatiPubblicati->GetNbinsX(); k++){
    fHistYieldDatiPubblicati->SetBinError(k,     fHistYieldErroriDatiPubblicati->GetBinContent(k));
    }
    }
    //    fHistYieldDatiPubblicati->Scale(1.6);
    for(Int_t j=0; j<3; j++){
      fHistYieldvsErrSoloStat[l][j]->Scale(1/1.6);
      if (!isSE)     fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 1.5);
      if (ishhCorr)       fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0, 4);
       if (isSE) fHistYieldvsErrSoloStat[l][j]->GetYaxis()->SetRangeUser(0,0.08);
      fHistYieldvsErrSoloStat[l][j]->GetXaxis()->SetRangeUser(0,25);
      // fHistYieldvsErrSoloStat[l][j]->SetMarkerColor(Color[j]);
      //fHistYieldvsErrSoloStat[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      // fHistYieldvsErrSoloStat[l][j]->SetMarkerSize(msize);
      // if (j==2)    fHistYieldvsErrSoloStat[l][j]->SetMarkerSize(3);
      // fHistYieldvsErrSoloStat[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrSoloStat[l][j]->Draw("samee");
        
      fHistYieldvsErrSoloSist[l][j]->Scale(1/1.6);
      if (!isSE)      fHistYieldvsErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0, 1.5);
      if (ishhCorr)      fHistYieldvsErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0, 4);
      if (isSE) fHistYieldvsErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0,0.08);
      fHistYieldvsErrSoloSist[l][j]->GetXaxis()->SetRangeUser(0,25);
      // fHistYieldvsErrSoloSist[l][j]->SetMarkerColor(Color[j]);
      //       fHistYieldvsErrSoloSist[l][j]->SetMarkerStyle(MarkerOrigin[j]);
      // fHistYieldvsErrSoloSist[l][j]->SetLineColor(Color[j]);
      // fHistYieldvsErrSoloSist[l][j]->SetFillStyle(0);
      // fHistYieldvsErrSoloSist[l][j]->SetMarkerSize(msize);
      // if (j==2)    fHistYieldvsErrSoloSist[l][j]->SetMarkerSize(3);
      fHistYieldvsErrSoloSist[l][j]->Draw("samee2");

      if (!ishhCorr && !isSE && l==0) fHistYieldDatiPubblicati->Draw("same");
        
      if(j==0)legendbis[l]->Draw();
      if(j==0)legendtype[l]->Draw();
    }
    //    canvasnotratio[l]->SaveAs("FinalOutput/DATA2016/Yield"+DATAorMC[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
    //  if (MasterThesis)         canvasnotratio[l]->SaveAs("FinalOutput/DATA2016/Yield"+DATAorMC[l]+".pdf");
  }

  cout << "\n\n\ndisegno errori relativi in tre canvas (1. JET, 2. OJ 3.J+OJ)" << endl;
  TCanvas *canvaserr[3];


  for(Int_t l=0; l<3; l++){
    Int_t count=0;
    canvaserr[l] = new TCanvas(Form("canvaserr%i", l), Region[l], 1300, 1000);
    cout << l << endl;
   
    for(Int_t j=0; j<2; j++){
      //    cout << l << "  " << j << endl;
      count++;    
      //      fHistYieldvsErrErrSoloStat[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloStat[j][l]->SetMarkerColor(2);
      fHistYieldvsErrErrSoloStat[j][l]->SetMarkerStyle(Marker[j]);
      fHistYieldvsErrErrSoloStat[j][l]->SetLineColor(2);
      fHistYieldvsErrErrSoloStat[j][l]->Draw("sameep");

      fHistYieldvsErrErrSoloStat[j][l]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[j][l]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[j][l]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[j][l]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrErrSoloStat[j][l]->SetMarkerSize(msize);
      
      //      fHistYieldvsErrErrSoloStatMB[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloStatMB[j][l]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloStatMB[j][l]->SetMarkerStyle(MarkerMBStat[j]);
      fHistYieldvsErrErrSoloStatMB[j][l]->SetLineColor(1);
      fHistYieldvsErrErrSoloStatMB[j][l]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloStatMB[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrErrSoloStatMB[j][l]->Draw("sameep");

      //      fHistYieldvsErrErrSoloSist[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloSist[j][l]->SetMarkerColor(4);
      fHistYieldvsErrErrSoloSist[j][l]->SetMarkerStyle(Marker[j]);
      fHistYieldvsErrErrSoloSist[j][l]->SetLineColor(4);
      fHistYieldvsErrErrSoloSist[j][l]->SetFillStyle(0);
      fHistYieldvsErrErrSoloSist[j][l]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSist[j][l]->Draw("sameep");

      fHistYieldvsErrErrSoloSist[j][l]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[j][l]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[j][l]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[j][l]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[j][l]->GetYaxis()->SetRangeUser(0,0.14);
      fHistYieldvsErrErrSoloSist[j][l]->GetXaxis()->SetRangeUser(0,25);

      //      fHistYieldvsErrErrSoloSistMB[j][l]->SetTitle("");
      fHistYieldvsErrErrSoloSistMB[j][l]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloSistMB[j][l]->SetMarkerStyle(MarkerMBSist[j]);
      fHistYieldvsErrErrSoloSistMB[j][l]->SetLineColor(1);
      fHistYieldvsErrErrSoloSistMB[j][l]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSistMB[j][l]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrErrSoloSistMB[j][l]->Draw("sameep");

        
      if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloStat[j][l], StatErr[j], "pl");
      if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloSist[j][l], SysErr[j], "pl");
      if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloStatMB[j][l], StatErrMB[j], "pl");
      if(l==0)legenderr->AddEntry(    fHistYieldvsErrErrSoloSistMB[j][l], SysErrMB[j], "pl");
      if(count==2)legenderr->Draw();
    
    }
    //    canvaserr[l] ->SaveAs("FinalOutput/DATA2016/RelErr"+Region[l]+hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin));
  if (MasterThesis)     canvaserr[l] ->SaveAs("FinalOutput/DATA2016/RelErr"+Region[l]+".pdf");
  }


  cout << "\n\n\ndisegno errori relativi in due canvas (1. DATA, 2. MC)" << endl;
  TCanvas *canvaserrall[2];
  
  for(Int_t l=0; l<2; l++){
    Int_t count=0;
    canvaserrall[l] = new TCanvas(Form("canvaserrall%i", l), Form("canvaserrall%i", l), 1300, 1000);
    cout << l << endl;
   
    for(Int_t j=0; j<3; j++){
      //    cout << l << "  " << j << endl;
      count++;
      //      fHistYieldvsErrErrSoloStat[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloStat[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrErrSoloStat[l][j]->SetMarkerStyle(MarkerTris[0]);
      fHistYieldvsErrErrSoloStat[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrErrSoloStat[l][j]->Draw("sameep");

      fHistYieldvsErrErrSoloStat[l][j]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[l][j]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[l][j]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloStat[l][j]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloStat[l][j]->GetXaxis()->SetRangeUser(0,25);
      fHistYieldvsErrErrSoloStat[l][j]->SetMarkerSize(msize);
            
      //      fHistYieldvsErrErrSoloStatMB[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloStatMB[l][j]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloStatMB[l][j]->SetMarkerStyle(MarkerTris[0]);
      fHistYieldvsErrErrSoloStatMB[l][j]->SetLineColor(1);
      fHistYieldvsErrErrSoloStatMB[l][j]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloStatMB[l][j]->GetXaxis()->SetRangeUser(0,5);
      //fHistYieldvsErrErrSoloStatMB[l][j]->Draw("sameep");

      //      fHistYieldvsErrErrSoloSist[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloSist[l][j]->SetMarkerColor(Color[j]);
      fHistYieldvsErrErrSoloSist[l][j]->SetMarkerStyle(MarkerTris[1]);
      fHistYieldvsErrErrSoloSist[l][j]->SetLineColor(Color[j]);
      fHistYieldvsErrErrSoloSist[l][j]->SetFillStyle(0);
      fHistYieldvsErrErrSoloSist[l][j]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSist[l][j]->Draw("sameep");

      fHistYieldvsErrErrSoloSist[l][j]->GetXaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[l][j]->GetXaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[l][j]->GetYaxis()->SetTitleOffset(1.1);
      fHistYieldvsErrErrSoloSist[l][j]->GetYaxis()->SetTitleSize(0.042);
      fHistYieldvsErrErrSoloSist[l][j]->GetYaxis()->SetRangeUser(0,0.14);
      fHistYieldvsErrErrSoloSist[l][j]->GetXaxis()->SetRangeUser(0,25);

      //      fHistYieldvsErrErrSoloSistMB[l][j]->SetTitle("");
      fHistYieldvsErrErrSoloSistMB[l][j]->SetMarkerColor(1);
      fHistYieldvsErrErrSoloSistMB[l][j]->SetMarkerStyle(MarkerTris[1]);
      fHistYieldvsErrErrSoloSistMB[l][j]->SetLineColor(1);
      fHistYieldvsErrErrSoloSistMB[l][j]->SetMarkerSize(msize);
      fHistYieldvsErrErrSoloSistMB[l][j]->GetXaxis()->SetRangeUser(0,25);
      //fHistYieldvsErrErrSoloSistMB[l][j]->Draw("sameep");
       
      if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloStat[l][j], StatErrBis[j], "pl");
      if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloSist[l][j], SysErrBis[j], "pl");
      // if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloStatMB[l][j], StatErrMB[j], "pl");
      // if(l==0)legenderrbis->AddEntry(    fHistYieldvsErrErrSoloSistMB[l][j], SysErrMB[j], "pl");
      if(count==3)legenderrbis->Draw();

    }


  }
  
  //  canvaserrall[0] ->SaveAs("FinalOutput/DATA2016/RelErrAllDATA0" +hhCorr[WhatKindOf]+Form("_PtMin%.1f.pdf", PtTrigMin)); 
  //  if (MasterThesis)  canvaserrall[0] ->SaveAs("FinalOutput/DATA2016/RelErrAllDATA0.pdf"); 
  //canvaserrall[1] ->SaveAs(Form("FinalOutput/DATA2016/RelErrAllDATA1_PtMin%.1f.pdf", PtTrigMin)); 
 
}

