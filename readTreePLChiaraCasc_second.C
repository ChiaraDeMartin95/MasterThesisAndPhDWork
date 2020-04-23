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
void readTreePLChiaraCasc_second(Int_t type=0, Int_t israp=1, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Float_t PtTrigMin=0.15,Float_t PtTrigMinFit=0.15, Int_t sysTrigger=0, Int_t sysV0=0,Int_t syst=0,bool isMC = 1, Bool_t isEfficiency=1,TString year0="2016", TString year="2018f1_extra_hXi_25runs_Bis",  TString Path1 ="",  Double_t ptjmax =15, Double_t nsigmamax=10, Bool_t isSigma=kFALSE){

  //isMeanFixedPDG and isBkgParab are characteristics of the fit to the inv mass distributions 
  cout << isMC << endl;
  cout << " Pt Trigg Min è = " << PtTrigMin << endl;

  if (israp>1) return;
  if (type>5) {cout << "type value not allowed" << endl; return;}
  //lista degli effetti  sistematici studiati in questa macro
  //sys=1 nsigmamin=5 (def:4)
  //sys=2 sigmacentral =4 (def:3)
  if (sysV0>6) return;

  if (syst!=0 && syst!=1 && syst!=2) {
    cout << "syst should be changed " << endl;
    return;
  }
  if((sysV0!=0||sysTrigger!=0) && syst!=0){
    cout << "syst conflicting " << endl;
    return;
  }

  Double_t sigmacentral=4;
  Double_t nsigmamin=4;
  if (syst==1) {
    nsigmamin=5;
  }
  if (syst==2) {
    sigmacentral=3;
  }

  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  const Float_t ctauK0s=2.6844;

  TString file;

  const Int_t numtipo=6;
  TString tipo[numtipo]={"XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  Float_t ctauCasc[numtipo] = {4.91,4.91,  2.461, 2.461, 4.91, 2.461}; //cm , average ctau of Xi and Omega                                                
  Float_t PDGCode[numtipo-2] = {3312, -3312, 3334, -3334};
  Float_t LimInfMass[numtipo]= {1.29, 1.29, 1.5, 1.5, 1.29, 1.5};
  Float_t LimSupMass[numtipo]= {1.35, 1.35, 1.8, 1.8, 1.35, 1.5};
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};

  TString PathIn="FinalOutput/AnalysisResults";
  TString PathOut="FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";
  TString PathInMass= "FinalOutput/DATA" + year0 + "/invmass_distribution_thesis/invmass_distribution";
  TString PathInMassDef;
  PathIn+=year;
  PathOut+=year;  
  // PathInMass+=year;
  
  if(isMC && isEfficiency){ 
    PathIn+="_MCEff";
    PathOut+="_MCEff_";
    PathInMass+="_MCEff";
  }
  if(isMC && !isEfficiency){
    PathIn+="_MCTruth";
    PathOut+="_MCTruth";
  }
 

  PathIn+=Path1;
  PathInMass+=Path1;
  PathIn+=".root";


  PathOut+=Path1;
  PathOut +=tipo[type];
  PathOut +=Srap[israp];
  if ((!isMC ||(isMC && isEfficiency))) PathOut +=Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f",sysTrigger, sysV0, syst, PtTrigMin); 
  PathOut+= ".root";

  cout << "file di input " << PathIn << endl;
  TFile *fin = new TFile(PathIn);
  TFile *fileMassSigma;
  TString dirinputtype[6] = {"Xi", "Xi", "Omega", "Omega", "Xi", "Omega"};
  TDirectoryFile *d = (TDirectoryFile*)fin->Get("MyTask"+dirinputtype[type]);
  TTree *tSign = (TTree *)d->Get("fSignalTree");
  TTree *tBkg  = (TTree *)d->Get("fBkgTree");

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=8;
  const Int_t numPtTrigger=1;
  
 
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-0.5","0.5-1.0", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,0.5, 1,1.5, 2,2.5,3,4,8};
  //  TString SPtTrigger[numPtTrigger]={"2-10"};
  Double_t NPtTrigger[numPtTrigger+1]={PtTrigMin,ptjmax};


  Double_t sigma[numtipo][nummolt+1][numPtV0];
  Double_t mass[numtipo][nummolt+1][numPtV0];
  Double_t meansigma;
  Double_t meanmass;
  TH1F* histoSigma;
  TH1F* histoMean;

  if(isMC==0 || (isMC==1 && isEfficiency==1)){
    for(Int_t m=0; m<nummolt+1; m++){
      //      PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, syst, PtTrigMin);
      PathInMassDef=PathInMass+ "_"+year+"_"+tipo[type]+Srap[israp]+"_"+MassFixedPDG[isMeanFixedPDG] + BkgType[isBkgParab] +Form("_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", m, sysTrigger, sysV0, syst,PtTrigMinFit);
      fileMassSigma= new TFile(PathInMassDef);
      histoSigma=(TH1F*)fileMassSigma->Get("histo_sigma");
      histoMean=(TH1F*)fileMassSigma->Get("histo_mean");
      for(Int_t v=1; v<numPtV0; v++){
	mass[type][m][v]=histoMean->GetBinContent(v+1);
	sigma[type][m][v]=histoSigma->GetBinContent(v+1);
	cout <<"mult interval " <<  m << " PtV0 interval " << v << " mean " << mass[type][m][v] << " sigma "<< sigma[type][m][v] << endl;
      }
    }
  }

  Double_t     fSignTreeVariablePtTrigger;
  Int_t        fSignTreeVariableChargeTrigger;
  Double_t     fSignTreeVariableEtaTrigger;
  Double_t     fSignTreeVariablePhiTrigger;
  Double_t     fSignTreeVariableDCAz;
  Double_t     fSignTreeVariableDCAxy;
  Int_t        fSignTreeVariableChargeAssoc;
  Int_t        fSignTreeVariableisPrimaryTrigger;
  Int_t        fSignTreeVariableisPrimaryV0;
  Double_t     fSignTreeVariableRapAssoc;
  Double_t     fSignTreeVariableDcaXiDaughters;
  Double_t     fSignTreeVariableV0CosineOfPointingAngle;
  Double_t     fSignTreeVariableXiCosineOfPointingAngle;
  Double_t     fSignTreeVariablectau;
  Double_t     fSignTreeVariablePtV0;
  Double_t     fSignTreeVariableInvMassLambda;
  Double_t     fSignTreeVariableInvMassXi;
  Double_t     fSignTreeVariableInvMassOmega;
  Double_t     fSignTreeVariableEtaV0;
  Double_t     fSignTreeVariablePhiV0;
  Double_t     fSignTreeVariableDeltaEta;
  Double_t     fSignTreeVariableDeltaPhi;
  Double_t     fSignTreeVariableDeltaTheta;
  Double_t     fSignTreeVariableMultiplicity;
  Double_t     fSignTreeVariableZvertex;
  Int_t        fSignTreeVariablePDGCodeTrigger;
  Int_t        fSignTreeVariablePDGCodeAssoc;
  Bool_t        fSignTreeVariableSkipAssoc;


  Double_t     fBkgTreeVariablePtTrigger;
  Int_t        fBkgTreeVariableChargeTrigger;
  Double_t     fBkgTreeVariableEtaTrigger;
  Double_t     fBkgTreeVariablePhiTrigger;
  Double_t     fBkgTreeVariableDCAz;
  Double_t     fBkgTreeVariableDCAxy;
  Int_t        fBkgTreeVariableChargeAssoc;
  Int_t        fBkgTreeVariableisPrimaryTrigger;
  Int_t        fBkgTreeVariableisPrimaryV0;
  Double_t     fBkgTreeVariableRapAssoc;
  Double_t     fBkgTreeVariableDcaXiDaughters;
  Double_t     fBkgTreeVariableV0CosineOfPointingAngle;
  Double_t     fBkgTreeVariableXiCosineOfPointingAngle;
  Double_t     fBkgTreeVariablectau;
  Double_t     fBkgTreeVariablePtV0;
  Double_t     fBkgTreeVariableInvMassLambda;
  Double_t     fBkgTreeVariableInvMassXi;
  Double_t     fBkgTreeVariableInvMassOmega;
  Double_t     fBkgTreeVariableEtaV0;
  Double_t     fBkgTreeVariablePhiV0;
  Double_t     fBkgTreeVariableDeltaEta;
  Double_t     fBkgTreeVariableDeltaPhi;
  Double_t     fBkgTreeVariableDeltaTheta;
  Double_t     fBkgTreeVariableMultiplicity;
  Double_t     fBkgTreeVariableZvertex;
  Int_t        fBkgTreeVariablePDGCodeTrigger;
  Int_t        fBkgTreeVariablePDGCodeAssoc;
  Bool_t        fBkgTreeVariableSkipAssoc;

 
  //Signal
  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);
  tSign->SetBranchAddress("fTreeVariableChargeAssoc"               ,&fSignTreeVariableChargeAssoc);
  tSign->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger);
  tSign->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0);
  tSign->SetBranchAddress("fTreeVariableRapAssoc"                  ,&fSignTreeVariableRapAssoc);
  tSign->SetBranchAddress("fTreeVariableDcaXiDaughters"            ,&fSignTreeVariableDcaXiDaughters);
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fSignTreeVariableV0CosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariableXiCosineOfPointingAngle"   ,&fSignTreeVariableXiCosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariablectau"                      ,&fSignTreeVariablectau);
  tSign->SetBranchAddress("fTreeVariablePtV0"                      ,&fSignTreeVariablePtV0);
  tSign->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fSignTreeVariableInvMassLambda);
  tSign->SetBranchAddress("fTreeVariableInvMassXi"                 ,&fSignTreeVariableInvMassXi);
  tSign->SetBranchAddress("fTreeVariableInvMassOmega"              ,&fSignTreeVariableInvMassOmega);
  tSign->SetBranchAddress("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0);
  tSign->SetBranchAddress("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0);
  tSign->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta);
  tSign->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi);
  tSign->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta);
  tSign->SetBranchAddress("fTreeVariableMultiplicity"              ,&fSignTreeVariableMultiplicity);
  tSign->SetBranchAddress("fTreeVariableZvertex"                   ,&fSignTreeVariableZvertex);
  tSign->SetBranchAddress("fTreeVariablePDGCodeTrigger"            ,&fSignTreeVariablePDGCodeTrigger);
  tSign->SetBranchAddress("fTreeVariablePDGCodeAssoc"              ,&fSignTreeVariablePDGCodeAssoc);
  tSign->SetBranchAddress("fTreeVariableSkipAssoc"                 ,&fSignTreeVariableSkipAssoc);

  //BackGround                                                                                                                                                                       
  tBkg->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger);
  tBkg->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger);
  tBkg->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger);
  tBkg->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger);
  tBkg->SetBranchAddress("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz);
  tBkg->SetBranchAddress("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy);    
  tBkg->SetBranchAddress("fTreeVariableChargeAssoc"               ,&fBkgTreeVariableChargeAssoc);
  tBkg->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger);
  tBkg->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0);
  tBkg->SetBranchAddress("fTreeVariableRapAssoc"                  ,&fBkgTreeVariableRapAssoc);
  tBkg->SetBranchAddress("fTreeVariableDcaXiDaughters"            ,&fBkgTreeVariableDcaXiDaughters);
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fBkgTreeVariableV0CosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariableXiCosineOfPointingAngle"   ,&fBkgTreeVariableXiCosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariablectau"                      ,&fBkgTreeVariablectau);
  tBkg->SetBranchAddress("fTreeVariablePtV0"                      ,&fBkgTreeVariablePtV0);
  tBkg->SetBranchAddress("fTreeVariableInvMassXi"                 ,&fBkgTreeVariableInvMassXi);
  tBkg->SetBranchAddress("fTreeVariableInvMassOmega"              ,&fBkgTreeVariableInvMassOmega);
  tBkg->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fBkgTreeVariableInvMassLambda);
  tBkg->SetBranchAddress("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0);
  tBkg->SetBranchAddress("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0);
  tBkg->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta);
  tBkg->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi);
  tBkg->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta);
  tBkg->SetBranchAddress("fTreeVariableMultiplicity"              ,&fBkgTreeVariableMultiplicity);
  tBkg->SetBranchAddress("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex);
  tBkg->SetBranchAddress("fTreeVariablePDGCodeTrigger"            ,&fBkgTreeVariablePDGCodeTrigger);
  tBkg->SetBranchAddress("fTreeVariablePDGCodeAssoc"              ,&fBkgTreeVariablePDGCodeAssoc);
  tBkg->SetBranchAddress("fTreeVariableSkipAssoc"                 ,&fBkgTreeVariableSkipAssoc);

  Int_t EntriesSign = 0; 
  Int_t EntriesBkg  = 0; 

  
  TFile *fout = new TFile(PathOut,"RECREATE");
  TDirectory  *dirSign= fout->mkdir("SE");
  TDirectory  *dirBkg= fout->mkdir("ME");

  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hSign_PtTrigger[nummolt+1][numzeta];
  TH1D *hSign_PtV0[nummolt+1][numzeta];

  const   Int_t numDeltaEta=4;
  TString SDeltaEta[numDeltaEta]={"_|DeltaEta|<1.6","_|DeltaEta|<1.4", "_|DeltaEta|<1.1","_|DeltaEta|<0.8" };
  Float_t DeltaEtaLimit[numDeltaEta]={1.6, 1.4, 1.1, 0.8};
  TH1D *hSign_EtaV0[nummolt+1][numDeltaEta];
  TH1D *hSign_EtaV0MidRap[nummolt+1][numDeltaEta];
  TH1D *hSign_RapV0[nummolt+1][numDeltaEta];
  TH2D *hSign_DeltaEtaEtaV0[nummolt+1];

  for(Int_t m=0; m<nummolt+1; m++){
    hSign_DeltaEtaEtaV0[m]=new TH2D("hSign_DeltaEtaEtaV0_"+Smolt[m], "hSign_DeltaEtaEtaV0_"+Smolt[m], 100, -2, 2, 100, -2, 2);
    for (Int_t DeltaEta=0; DeltaEta< numDeltaEta; DeltaEta++){
      hSign_EtaV0[m][DeltaEta]=new TH1D("hSign_EtaV0_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_EtaV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hSign_RapV0[m][DeltaEta]=new TH1D("hSign_RapV0_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_RapV0_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
      hSign_EtaV0MidRap[m][DeltaEta]=new TH1D("hSign_EtaV0_MidRap_"+Smolt[m]+Form("_%i",DeltaEta), "hSign_EtaV0_MidRap_"+Smolt[m]+SDeltaEta[DeltaEta], 100, -2, 2);
    }
  }

  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      hSign_PtTrigger[m][z]=new TH1D("hSign_PtTrigger"+Smolt[m], "hSign_PtTrigger"+Smolt[m], 300,0,30);
      hSign_PtV0[m][z]=new TH1D("hSign_PtV0"+Smolt[m], "hSign_PtV0"+Smolt[m], 300,0,30);
      hSign_PtTrigger[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtTrigger[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtTrigger[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtTrigger[m][z]->GetYaxis()->SetLabelSize(0.05);
      hSign_PtV0[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hSign_PtV0[m][z]->GetXaxis()->SetTitleSize(0.05);
      hSign_PtV0[m][z]->GetXaxis()->SetLabelSize(0.05);
      hSign_PtV0[m][z]->GetYaxis()->SetLabelSize(0.05);

      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=1; v<numPtV0; v++){
	  nameSE[m][z][v][tr]="SE_";
	  namemassSE[m][z][v][tr]="InvMassSE_";
	  nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr], nameSE[m][z][v][tr],   50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr]+"_SB", nameSE[m][z][v][tr]+"_SB",   50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	}
      }
    }
  }

  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hBkg_PtTrigger[nummolt+1][numzeta];
  TH1D *hBkg_PtV0[nummolt+1][numzeta];

  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      hBkg_PtTrigger[m][z]=new TH1D("hBkg_PtTrigger"+Smolt[m], "hBkg_PtTrigger"+Smolt[m], 300,0,30);
      hBkg_PtV0[m][z]=new TH1D("hBkg_PtV0"+Smolt[m], "hBkg_PtV0"+Smolt[m], 300,0,30);
      hBkg_PtTrigger[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hBkg_PtTrigger[m][z]->GetXaxis()->SetTitleSize(0.05);
      hBkg_PtTrigger[m][z]->GetXaxis()->SetLabelSize(0.05);
      hBkg_PtTrigger[m][z]->GetYaxis()->SetLabelSize(0.05);
      hBkg_PtV0[m][z]->GetXaxis()->SetTitle("p_{T} (Gev/c)");
      hBkg_PtV0[m][z]->GetXaxis()->SetTitleSize(0.05);
      hBkg_PtV0[m][z]->GetXaxis()->SetLabelSize(0.05);
      hBkg_PtV0[m][z]->GetYaxis()->SetLabelSize(0.05);

      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=1; v<numPtV0; v++){
	  nameME[m][z][v][tr]="ME_";
	  namemassME[m][z][v][tr]="InvMassME_";
	  nameME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]= new TH2D(nameME[m][z][v][tr], nameME[m][z][v][tr],  50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]= new TH2D(nameME[m][z][v][tr]+"_SB", nameME[m][z][v][tr]+"_SB",  50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	}
      }
    }
  }

  EntriesSign =  tSign->GetEntries();
  EntriesBkg  =  tBkg ->GetEntries();
     
  Bool_t BoolVar=kFALSE;
  Bool_t BoolMC=kFALSE;
  Bool_t MassLimit=kFALSE;
  Float_t     fSignTreeVariableInvMassCasc= 0;
  Bool_t isCascTrue=kFALSE;

  dirSign->cd();
  cout << "\n\n I will process " << EntriesSign << " entries for theSE correlation " << endl;
  for(Int_t k = 0; k<EntriesSign; k++){
    //  for(Int_t k = 0; k<1000000; k++){
    for (Int_t l=0; l<10000; l++){
      //      if (k ==100000*l) cout << "k = " << k << " over a total of " << EntriesSign << endl;
      if (k ==100000*l) {
	Float_t kFloat= k;
	cout << "PtTrigMin = " << PtTrigMin << " sysV0 " << sysV0 << " SE, processing..." << kFloat/EntriesSign << endl;
      }
    }

    tSign->GetEntry(k);

    //charge selection                                                                                                                                    
    if ((type==0 || type ==2) && fSignTreeVariableChargeAssoc==1) continue;
    else if ((type==1 || type ==3) && fSignTreeVariableChargeAssoc==-1) continue;

    //inv mass definition                                                                                                                                 
    fSignTreeVariableInvMassCasc= 0;
    if (type==0 || type==1 || type==4)      fSignTreeVariableInvMassCasc= fSignTreeVariableInvMassXi;
    else if (type==2 || type==3 || type==5)      fSignTreeVariableInvMassCasc= fSignTreeVariableInvMassOmega;

    //rapidity selection                                                                                                                                  
    if (israp==0 && TMath::Abs(fSignTreeVariableEtaV0)>0.8)continue;
    else if (israp==1 && TMath::Abs(fSignTreeVariableRapAssoc)>0.5)continue;

    //definition of true cascade                                                                                                                          
    if (type<=3) isCascTrue= (fSignTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==4) isCascTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[0]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[1])  );
    else if (type==5) isCascTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[2]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[3])  );

    if (fSignTreeVariableSkipAssoc==1) continue;

    if(isMC==0 || (isMC==1 && isEfficiency==1)){

      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMin) continue;
      if(TMath::Abs(fSignTreeVariablePtTrigger)>ptjmax) continue;
      
      //cuts on DCAz trigger*******************
      if(sysTrigger==0){
	if(TMath::Abs(fSignTreeVariableDCAz)>1) continue;
      }
      if(sysTrigger==1){
	if(TMath::Abs(fSignTreeVariableDCAz)>2) continue;
      }
      if(sysTrigger==2){
	if(TMath::Abs(fSignTreeVariableDCAz)>0.5) continue;
      }
      
      //******************* some other cuts for sys studies**************************                                                                                                  
      if (sysV0==0){

	if (fSignTreeVariableDcaXiDaughters> 0.8) continue;
	if (fSignTreeVariableV0CosineOfPointingAngle<0.97)continue;
	if (fSignTreeVariableXiCosineOfPointingAngle < 0.995)continue;
	if (TMath::Abs(fSignTreeVariableInvMassLambda -massLambda) > 0.006)continue;
	if (fSignTreeVariablectau  > 3* ctauCasc[type])continue;

      } 

    }
    //**********************************************************************************************

    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();

    for(Int_t m=0; m<nummolt+1; m++){
      if(m< nummolt) BoolVar = fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for (Int_t DeltaEta=0; DeltaEta<numDeltaEta; DeltaEta++){
	hSign_DeltaEtaEtaV0[m]->Fill(fSignTreeVariableEtaV0, fSignTreeVariableDeltaEta);
	if (BoolVar && TMath::Abs(fSignTreeVariableDeltaEta) < DeltaEtaLimit[DeltaEta]){
	  hSign_EtaV0[m][DeltaEta]->Fill(fSignTreeVariableEtaV0);
	  hSign_RapV0[m][DeltaEta]->Fill(fSignTreeVariableRapAssoc);
	  if (TMath::Abs(fSignTreeVariableRapAssoc) < 0.5)	    hSign_EtaV0MidRap[m][DeltaEta]->Fill(fSignTreeVariableEtaV0);
	}
      }
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=1; v<numPtV0; v++){
	    if((isMC && !isEfficiency)) {
	      BoolMC = kTRUE;
	      MassLimit=kTRUE;
	    }
	    else {
	      BoolMC =TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))<sigmacentral*sigma[type][m][v]; 
	      if(isSigma) MassLimit=(TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin*sigma[type][m][v] && TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))<nsigmamax*sigma[type][m][v]);
	      if(!isSigma)MassLimit=(TMath::Abs((fSignTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin*sigma[type][m][v] && fSignTreeVariableInvMassCasc>LimInfMass[type] && fSignTreeVariableInvMassCasc<LimSupMass[type]);
	    }	    
	    if(BoolMC && BoolVar && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      hSign_PtTrigger[m][z]->Fill(fSignTreeVariablePtTrigger);
	      hSign_PtV0[m][z]->Fill(fSignTreeVariablePtV0);
	    }	  
	    if(BoolVar && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      if(BoolMC){
		hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }
	      if((!isMC || (isMC &&isEfficiency))&& MassLimit){
		hDeltaEtaDeltaPhi_SEbins_sidebands[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }

	    }
	  }
	}
      }
    }
  }


  BoolVar=kFALSE;
  BoolMC=kFALSE;
  MassLimit=kFALSE;

  Float_t     fBkgTreeVariableInvMassCasc= 0;
  cout << "ciao " << endl;
  dirBkg->cd();
  cout << "\n\n I will process " << EntriesBkg << " entries for theSE correlation " << endl;
  for(Int_t k = 0; k<EntriesBkg; k++){
    //  for(Int_t k = 0; k<1000000; k++){
    for (Int_t l=0; l<10000; l++){
      //if (k ==100000*l) cout << "k = " << k << " over a total of " << EntriesBkg << endl;
      if (k ==100000*l) {
	Float_t kFloat= k;
	cout << "PtTrigMin "<< PtTrigMin << " sysV0 "<< sysV0 << " ME processing..." << kFloat/EntriesBkg << endl;
      }
    }

    tBkg->GetEntry(k);

    //charge selection                                                                                                                                    
    if ((type==0 || type ==2) && fBkgTreeVariableChargeAssoc==1) continue;
    else if ((type==1 || type ==3) && fBkgTreeVariableChargeAssoc==-1) continue;

    //inv mass definition                                                                                                                                 
    fBkgTreeVariableInvMassCasc= 0;
    if (type==0 || type==1 || type==4)      fBkgTreeVariableInvMassCasc= fBkgTreeVariableInvMassXi;
    else if (type==2 || type==3 || type==5)      fBkgTreeVariableInvMassCasc= fBkgTreeVariableInvMassOmega;

    //rapidity selection                                                                                                                                  
    if (israp==0 && TMath::Abs(fBkgTreeVariableEtaV0)>0.8)continue;
    else if (israp==1 && TMath::Abs(fBkgTreeVariableRapAssoc)>0.5)continue;

    //definition of true cascade                                                                                                                         
    if (type<=3) isCascTrue= (fBkgTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==4) isCascTrue= ((fBkgTreeVariablePDGCodeAssoc==PDGCode[0]) ||(fBkgTreeVariablePDGCodeAssoc==PDGCode[1])  );
    else if (type==5) isCascTrue= ((fBkgTreeVariablePDGCodeAssoc==PDGCode[2]) ||(fBkgTreeVariablePDGCodeAssoc==PDGCode[3])  );

    if (fBkgTreeVariableSkipAssoc==1) continue;

    if(isMC==0 || (isMC==1 && isEfficiency==1)){
      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMin) continue;
      if(TMath::Abs(fBkgTreeVariablePtTrigger)>ptjmax) continue;

      //cuts on DCAz trigger*******************
      if(sysTrigger==0){
	if(TMath::Abs(fBkgTreeVariableDCAz)>1) continue;
      }
      if(sysTrigger==1){
	if(TMath::Abs(fBkgTreeVariableDCAz)>2) continue;
      }
      if(sysTrigger==2){
	if(TMath::Abs(fBkgTreeVariableDCAz)>0.5) continue;
      }

      //******************* some other cuts for sys studies **************************
      if (sysV0==0){
        if (fBkgTreeVariableDcaXiDaughters> 0.8) continue;
        if (fBkgTreeVariableV0CosineOfPointingAngle<0.97)continue;
        if (fBkgTreeVariableXiCosineOfPointingAngle < 0.995)continue;
        if (TMath::Abs(fBkgTreeVariableInvMassLambda -massLambda) > 0.006)continue;
        if (fBkgTreeVariablectau  > 3* ctauCasc[type])continue;
      }
    }
     
    //**********************************************************************************************

    if (fBkgTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fBkgTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fBkgTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fBkgTreeVariableDeltaPhi += 2.0*TMath::Pi();
    for(Int_t m=0; m<nummolt+1; m++){
      if(m< nummolt) BoolVar =  fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=1; v<numPtV0; v++){

	    if((isMC && !isEfficiency) ) {
	      BoolMC = kTRUE;
	      MassLimit=kTRUE;
	    }
	    else {
	      BoolMC =TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))<sigmacentral*sigma[type][m][v]; 
	      if(isSigma) MassLimit=(TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin*sigma[type][m][v] && TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))<nsigmamax*sigma[type][m][v]);
	      if(!isSigma)MassLimit=(TMath::Abs((fBkgTreeVariableInvMassCasc - mass[type][m][v]))>nsigmamin*sigma[type][m][v] && fBkgTreeVariableInvMassCasc>LimInfMass[type] && fBkgTreeVariableInvMassCasc<LimSupMass[type]);
	    }	    
 
	    if(BoolMC && BoolVar &&  fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      hBkg_PtTrigger[m][z]->Fill(fBkgTreeVariablePtTrigger);
	      hBkg_PtV0[m][z]->Fill(fBkgTreeVariablePtV0);
	    }	  
	    if(BoolVar && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      if(BoolMC){
		hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	      }
	      if((!isMC || (isMC &&isEfficiency)) && MassLimit){
		hDeltaEtaDeltaPhi_MEbins_sidebands[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "ciao " << endl;
  fout->Write();

  cout << "partendo dal file " << PathIn << " e " << PathInMassDef << " (per le diverse molteplicità) ho creato il file " << PathOut<< endl;
}

