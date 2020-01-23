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
void readTreePLChiara_Cascades(Int_t type=0, Bool_t ishhCorr=0, Float_t PtTrigMin=3.0, Int_t sysV0=0,Int_t rap=1, bool isMC = 1,Bool_t isEfficiency=1,Int_t sysTrigger=0,
			    //TString PathIn ="./AnalysisResults2018d8_allrunswMult_MC.root", Int_t type=0 //type = 0 per K0short
TString year="2018f1_extra", TString year0="2016", TString Path1 ="_25runs"	    
			    )
{
  if (type>1) {
    cout << "type=0 or 1 for Xi " << endl;
    return; 
  }

  if (sysV0>9) {
    cout << "sysV0 must be smaller than 10" << endl;
return;
  }
  //2016k_onlyTriggerWithHighestPt
  //  if (sysV0>6) return;
  //  if (sysV0>2 && ishhCorr) return;

  //systematics for Xi
  //SysV0:
  // 0: default cuts Fiorella Fionda 
  // 1: fSignTreeVariableV0CosineOfPointingAngleSpecial; 0.99
  // 2: fSignTreeVariableCascCosineOfPointingAngle;      0.99
  // 3: fSignTreeVariableDcaV0ToPrimVertex;              0.15
  // 4: fSignTreeVariableDcaMesonToPrimVertex;           0.30
  // 5: fSignTreeVariableDcaBaryonToPrimVertex           0.11
  // 6: fSignTreeVariableDcaBachToPrimVertex;            0.17
  // 7: fSignTreeVariableDcaV0Daughters;                 1
  // 8: fSignTreeVariableDcaCascDaughters;               0.8
  // 9: fSignTreeVariableInvMassLambda;                  0.006

  Float_t RapCasc[2]={100, 0.5};
  Float_t V0CosineOfPointingAngleSpecial[2]={0.97, 0.97}; 
  Float_t CascCosineOfPointingAngle[2]={0.97, 0.97};      
  Float_t DcaV0ToPrimVertex[2]={0.06, 0.06};              
  Float_t DcaMesonToPrimVertex[2]={0.04, 0.04};           
  Float_t DcaBaryonToPrimVertex[2]={0.03, 0.03};           
  Float_t DcaBachToPrimVertex[2]={0.04, 0.04};                       
  Float_t DcaV0Daughters[2]={1.5, 1.5};                      
  Float_t DcaCascDaughters[2]={1.3, 1.3};                 
  Float_t InvMassLambda[2]={0.008, 0.008};       

  if (sysV0==1) {
    V0CosineOfPointingAngleSpecial[0]=0.99; 
    V0CosineOfPointingAngleSpecial[1]=0.99; 
  }

  if (sysV0==2) {
  CascCosineOfPointingAngle[0]=0.99;
  CascCosineOfPointingAngle[1]=0.99;
  }
  if (sysV0==3) {
  DcaV0ToPrimVertex[0]=0.15;
  DcaV0ToPrimVertex[1]=0.10;
  }
  if (sysV0==4) {
  DcaMesonToPrimVertex[0]=0.30;
  DcaMesonToPrimVertex[1]=0.30;
  }
  if (sysV0==5) {
  DcaBaryonToPrimVertex[0]=0.11;
  DcaBaryonToPrimVertex[1]=0.10;
  }
  if (sysV0==6) {
  DcaBachToPrimVertex[0]=0.17;                       
  DcaBachToPrimVertex[1]=0.10;                       
  }
  if (sysV0==7) {
  DcaV0Daughters[0]=1.;
  DcaV0Daughters[1]=1.;
  }
  if (sysV0==8) {
  DcaCascDaughters[0]=0.8;
  DcaCascDaughters[1]=0.6;
  }
  if (sysV0==9) {
  InvMassLambda[0]=0.006;
  InvMassLambda[1]=0.006;
  }

  const Int_t numtipo=4;
  TString tipo[numtipo]={"XiNeg", "XiPos", "OmegaPos", "OmegaNeg"};
  Int_t CascadesPDGCode[numtipo]={3312, -3312, 3334, -3334};
  Float_t LimInfHisto[numtipo]={1.28, 1.28, 1.62, 1.62};
  Float_t LimSupHisto[numtipo]={1.36, 1.36, 1.72, 1.72};
  TString PathIn="./FinalOutput/AnalysisResults";
  TString PathOut="./FinalOutput/DATA" + year0 + "/histo/AngularCorrelation";
  PathIn+=year;
  PathOut+=year;  
  
  if(isMC && isEfficiency){
    PathIn+="_MCEff";
    PathOut+="_MCEff";
  }
  if(isMC && !isEfficiency){
    PathIn+="_MCTruth";
    PathOut+="_MCTruth";
  }

  PathIn+="_Cascades";
  PathOut+="_Cascades_";
  PathOut+=tipo[type];
  PathIn+=Path1;
  PathIn+=".root";
  PathOut+=Path1;
  if (ishhCorr) PathOut +="_hhCorr";
  PathOut +=Form("_MassDistr_SysT%i_SysV0%i_PtMin%.1f_Rap%i",sysTrigger, sysV0, PtTrigMin, rap); 
  PathOut+= ".root";

  TFile *fin = new TFile(PathIn);
  TDirectoryFile *d = (TDirectoryFile*)fin->Get("MyTask");
  TTree *tSign = (TTree *)d->Get("fSignalTree");
  TTree *tBkg  = (TTree *)d->Get("fBkgTree");
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Double_t ctauK0s = 2.6844;
  Double_t ctauCasc[numtipo] = {4.91, 4.91, 2.461, 2.461};
  Double_t massXi = 1.32171;
  Double_t massOmega = 1.67245;

  Double_t fSignTreeVariableMultiplicity;
  Double_t fSignTreeVariableZvertex;
  Double_t fSignTreeVariablePDGCode;
  Double_t fSignTreeVariablePDGCodeBach;
  Double_t fSignTreeVariablePDGCodePos;
  Double_t fSignTreeVariablePDGCodeNeg;
  Double_t fSignTreeVariablePDGCodeLambda;
  Double_t fSignTreeVariablePDGCodeMotherLambda;
  Double_t fSignTreeVariableRunNumber;
  Double_t fSignTreeVariableBunchCrossNumber;
  Double_t fSignTreeVariableNegNSigmaPion;
  Double_t fSignTreeVariableNegNSigmaProton;
  Double_t fSignTreeVariablePosNSigmaPion;
  Double_t fSignTreeVariablePosNSigmaProton;
  Double_t fSignTreeVariableBachNSigmaPion;
  Double_t fSignTreeVariableBachNSigmaKaon;
  Double_t fSignTreeVariableDcaXiToPrimVertex;
  Double_t fSignTreeVariableXYDcaXiToPrimVertex;
  Double_t fSignTreeVariableZDcaXiToPrimVertex;
  Double_t fSignTreeVariableDcaV0ToPrimVertex;
  Double_t fSignTreeVariableDcaPosToPrimVertex;
  Double_t fSignTreeVariableDcaNegToPrimVertex;
  Double_t fSignTreeVariableDcaV0Daughters;
  Double_t fSignTreeVariableDcaCascDaughters;
  Double_t fSignTreeVariableDcaBachToPrimVertex;
  Double_t fSignTreeVariableV0CosineOfPointingAngle;
  Double_t fSignTreeVariableV0CosineOfPointingAngleSpecial;
  Double_t fSignTreeVariableCascCosineOfPointingAngle;
  Double_t fSignTreeVariablePtCasc;
  Double_t fSignTreeVariableChargeCasc;
  Double_t fSignTreeVariableEtaCasc;
  Double_t fSignTreeVariablectau;
  Double_t fSignTreeVariableInvMassXi;
  Double_t fSignTreeVariableInvMassOmega;
  Double_t fSignTreeVariableInvMassLambda;
  Double_t fSignTreeVariableRapXi;
  Double_t fSignTreeVariableRapOmega;
  Double_t fSignTreeVariableCascRadius;
  Double_t fSignTreeVariableV0Radius;
  Double_t fSignTreeVariableLeastNbrClusters;
  Double_t fSignTreeVariableV0Lifetime;
  Double_t  fSignTreeVariableIsPrimaryXi;
  Double_t  fSignTreeVariableIsPrimaryOmega;

  Double_t fBkgTreeVariableMultiplicity;
  Double_t fBkgTreeVariableZvertex;
  Double_t fBkgTreeVariablePDGCode;
  Double_t fBkgTreeVariableRunNumber;
  Double_t fBkgTreeVariableBunchCrossNumber;
  Double_t fBkgTreeVariableNegNSigmaPion;
  Double_t fBkgTreeVariableNegNSigmaProton;
  Double_t fBkgTreeVariablePosNSigmaPion;
  Double_t fBkgTreeVariablePosNSigmaProton;
  Double_t fBkgTreeVariableBachNSigmaPion;
  Double_t fBkgTreeVariableBachNSigmaKaon;
  Double_t fBkgTreeVariableDcaXiToPrimVertex;
  Double_t fBkgTreeVariableXYDcaXiToPrimVertex;
  Double_t fBkgTreeVariableZDcaXiToPrimVertex;
  Double_t fBkgTreeVariableDcaV0ToPrimVertex;
  Double_t fBkgTreeVariableDcaPosToPrimVertex;
  Double_t fBkgTreeVariableDcaNegToPrimVertex;
  Double_t fBkgTreeVariableDcaV0Daughters;
  Double_t fBkgTreeVariableDcaCascDaughters;
  Double_t fBkgTreeVariableDcaBachToPrimVertex;
  Double_t fBkgTreeVariableV0CosineOfPointingAngle;
  Double_t fBkgTreeVariableV0CosineOfPointingAngleSpecial;
  Double_t fBkgTreeVariableCascCosineOfPointingAngle;
  Double_t fBkgTreeVariablePtCasc;
  Double_t fBkgTreeVariableChargeCasc;
  Double_t fBkgTreeVariableEtaCasc;
  Double_t fBkgTreeVariablectau;
  Double_t fBkgTreeVariableInvMassXi;
  Double_t fBkgTreeVariableInvMassOmega;
  Double_t fBkgTreeVariableInvMassLambda;
  Double_t fBkgTreeVariableRapXi;
  Double_t fBkgTreeVariableRapOmega;
  Double_t fBkgTreeVariableCascRadius;
  Double_t fBkgTreeVariableV0Radius;
  Double_t fBkgTreeVariableLeastNbrClusters;
  Double_t fBkgTreeVariableV0Lifetime;
  Double_t  fBkgTreeVariableIsPrimaryXi;
  Double_t  fBkgTreeVariableIsPrimaryOmega;


  Int_t CounterSignPairsAfterPtMinCut=0; 
  Int_t CounterBkgPairsAfterPtMinCut=0; 

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;

  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  //TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,1.5,2,2.5,3,4,8};
  //   TString SPtV0[numPtV0]={"0-1", "1-2", "2,3", "3-4", "4-8"};
  //   Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  TString SPtTrigger[numPtTrigger]={"2-10"};
  Double_t NPtTrigger[numPtTrigger+1]={2,10};

  //Signal
  tSign->SetBranchAddress("fTreeVariableMultiplicity",                   &fSignTreeVariableMultiplicity);
  tSign->SetBranchAddress("fTreeVariableZvertex",                        &fSignTreeVariableZvertex);
  tSign->SetBranchAddress("fTreeVariablePDGCode",                        &fSignTreeVariablePDGCode);
  tSign->SetBranchAddress("fTreeVariablePDGCodeBach",                    &fSignTreeVariablePDGCodeBach);
  tSign->SetBranchAddress("fTreeVariablePDGCodePos",                     &fSignTreeVariablePDGCodePos);
  tSign->SetBranchAddress("fTreeVariablePDGCodeNeg",                     &fSignTreeVariablePDGCodeNeg);
  tSign->SetBranchAddress("fTreeVariablePDGCodeLambda",                  &fSignTreeVariablePDGCodeLambda);
  tSign->SetBranchAddress("fTreeVariablePDGCodeMotherLambda",            &fSignTreeVariablePDGCodeMotherLambda);
  tSign->SetBranchAddress("fTreeVariableRunNumber",                      &fSignTreeVariableRunNumber);
  tSign->SetBranchAddress("fTreeVariableBunchCrossNumber",               &fSignTreeVariableBunchCrossNumber);
  tSign->SetBranchAddress("fTreeVariableNegNSigmaPion",                  &fSignTreeVariableNegNSigmaPion);
  tSign->SetBranchAddress("fTreeVariableNegNSigmaProton",                &fSignTreeVariableNegNSigmaProton);
  tSign->SetBranchAddress("fTreeVariablePosNSigmaPion",                  &fSignTreeVariablePosNSigmaPion);
  tSign->SetBranchAddress("fTreeVariablePosNSigmaProton",                &fSignTreeVariablePosNSigmaProton);
  tSign->SetBranchAddress("fTreeVariableBachNSigmaPion",                 &fSignTreeVariableBachNSigmaPion);
  tSign->SetBranchAddress("fTreeVariableBachNSigmaKaon",                 &fSignTreeVariableBachNSigmaKaon);
  tSign->SetBranchAddress("fTreeVariableDcaXiToPrimVertex",              &fSignTreeVariableDcaXiToPrimVertex); 
  tSign->SetBranchAddress("fTreeVariableXYDcaXiToPrimVertex",            &fSignTreeVariableXYDcaXiToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableZDcaXiToPrimVertex",             &fSignTreeVariableZDcaXiToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex",              &fSignTreeVariableDcaV0ToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",             &fSignTreeVariableDcaPosToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",             &fSignTreeVariableDcaNegToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableDcaV0Daughters",                 &fSignTreeVariableDcaV0Daughters);
  tSign->SetBranchAddress("fTreeVariableDcaCascDaughters",               &fSignTreeVariableDcaCascDaughters);
  tSign->SetBranchAddress("fTreeVariableDcaBachToPrimVertex",            &fSignTreeVariableDcaBachToPrimVertex);
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",        &fSignTreeVariableV0CosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngleSpecial", &fSignTreeVariableV0CosineOfPointingAngleSpecial);
  tSign->SetBranchAddress("fTreeVariableCascCosineOfPointingAngle",      &fSignTreeVariableCascCosineOfPointingAngle);
  tSign->SetBranchAddress("fTreeVariablePtCasc",                         &fSignTreeVariablePtCasc);
  tSign->SetBranchAddress("fTreeVariableChargeCasc",                     &fSignTreeVariableChargeCasc);
  tSign->SetBranchAddress("fTreeVariableEtaCasc",                        &fSignTreeVariableEtaCasc);
  tSign->SetBranchAddress("fTreeVariablectau",                           &fSignTreeVariablectau);
  tSign->SetBranchAddress("fTreeVariableInvMassXi",                      &fSignTreeVariableInvMassXi);
  tSign->SetBranchAddress("fTreeVariableInvMassOmega",                   &fSignTreeVariableInvMassOmega);
  tSign->SetBranchAddress("fTreeVariableInvMassLambda",                  &fSignTreeVariableInvMassLambda);
  tSign->SetBranchAddress("fTreeVariableRapXi",                          &fSignTreeVariableRapXi);
  tSign->SetBranchAddress("fTreeVariableRapOmega",                       &fSignTreeVariableRapOmega);
  tSign->SetBranchAddress("fTreeVariableCascRadius",                     &fSignTreeVariableCascRadius);
  tSign->SetBranchAddress("fTreeVariableV0Radius",                       &fSignTreeVariableV0Radius);
  tSign->SetBranchAddress("fTreeVariableLeastNbrClusters",               &fSignTreeVariableLeastNbrClusters);
  tSign->SetBranchAddress("fTreeVariableV0Lifetime",                     &fSignTreeVariableV0Lifetime);
  tSign->SetBranchAddress("fTreeVariableIsPrimaryXi",                    &fSignTreeVariableIsPrimaryXi);
  tSign->SetBranchAddress("fTreeVariableIsPrimaryOmega",                 &fSignTreeVariableIsPrimaryOmega);

  //other variables to be used when AC is performed
  /*
  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);		       
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);		       
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);		       
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);		       
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);			       
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);			       
  //hh tSign->SetBranchAddress("fTreeVariableChargeAssoc"             ,&fSignTreeVariableChargeAssoc);		       	   
  //hh tSign->SetBranchAddress("fTreeVariableAssocDCAz"                      ,&fSignTreeVariableAssocDCAz);			       
  //hh tSign->SetBranchAddress("fTreeVariableAssocDCAxy"                     ,&fSignTreeVariableAssocDCAxy);	
  tSign->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger);			       
  tSign->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0);			       
  tSign->SetBranchAddress("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0);			       
  tSign->SetBranchAddress("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0);			       
  tSign->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta);		       
  tSign->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi);		       
  tSign->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta);		       
  */

  //BackGround
  /*
  tBkg->SetBranchAddress("fTreeVariableMultiplicity",                   &fBkgTreeVariableMultiplicity);
  tBkg->SetBranchAddress("fTreeVariableZvertex",                        &fBkgTreeVariableZvertex);
  tBkg->SetBranchAddress("fTreeVariablePDGCode",                        &fBkgTreeVariablePDGCode);
  tBkg->SetBranchAddress("fTreeVariableRunNumber",                      &fBkgTreeVariableRunNumber);
  tBkg->SetBranchAddress("fTreeVariableBunchCrossNumber",               &fBkgTreeVariableBunchCrossNumber);
  tBkg->SetBranchAddress("fTreeVariableNegNSigmaPion",                  &fBkgTreeVariableNegNSigmaPion);
  tBkg->SetBranchAddress("fTreeVariableNegNSigmaProton",                &fBkgTreeVariableNegNSigmaProton);
  tBkg->SetBranchAddress("fTreeVariablePosNSigmaPion",                  &fBkgTreeVariablePosNSigmaPion);
  tBkg->SetBranchAddress("fTreeVariablePosNSigmaProton",                &fBkgTreeVariablePosNSigmaProton);
  tBkg->SetBranchAddress("fTreeVariableBachNSigmaPion",                 &fBkgTreeVariableBachNSigmaPion);
  tBkg->SetBranchAddress("fTreeVariableBachNSigmaKaon",                 &fBkgTreeVariableBachNSigmaKaon);
  tBkg->SetBranchAddress("fTreeVariableDcaXiToPrimVertex",              &fBkgTreeVariableDcaXiToPrimVertex); 
  tBkg->SetBranchAddress("fTreeVariableXYDcaXiToPrimVertex",            &fBkgTreeVariableXYDcaXiToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableZDcaXiToPrimVertex",             &fBkgTreeVariableZDcaXiToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex",              &fBkgTreeVariableDcaV0ToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableDcaPosToPrimVertex",             &fBkgTreeVariableDcaPosToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableDcaNegToPrimVertex",             &fBkgTreeVariableDcaNegToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableDcaV0Daughters",                 &fBkgTreeVariableDcaV0Daughters);
  tBkg->SetBranchAddress("fTreeVariableDcaCascDaughters",               &fBkgTreeVariableDcaCascDaughters);
  tBkg->SetBranchAddress("fTreeVariableDcaBachToPrimVertex",            &fBkgTreeVariableDcaBachToPrimVertex);
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle",        &fBkgTreeVariableV0CosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngleSpecial", &fBkgTreeVariableV0CosineOfPointingAngleSpecial);
  tBkg->SetBranchAddress("fTreeVariableCascCosineOfPointingAngle",      &fBkgTreeVariableCascCosineOfPointingAngle);
  tBkg->SetBranchAddress("fTreeVariablePtCasc",                         &fBkgTreeVariablePtCasc);
  tBkg->SetBranchAddress("fTreeVariableChargeCasc",                     &fBkgTreeVariableChargeCasc);
  tBkg->SetBranchAddress("fTreeVariablectau",                           &fBkgTreeVariablectau);
  tBkg->SetBranchAddress("fTreeVariableInvMassXi",                      &fBkgTreeVariableInvMassXi);
  tBkg->SetBranchAddress("fTreeVariableInvMassOmega",                   &fBkgTreeVariableInvMassOmega);
  tBkg->SetBranchAddress("fTreeVariableInvMassLambda",                  &fBkgTreeVariableInvMassLambda);
  tBkg->SetBranchAddress("fTreeVariableRapXi",                          &fBkgTreeVariableRapXi);
  tBkg->SetBranchAddress("fTreeVariableRapOmega",                       &fBkgTreeVariableRapOmega);
  tBkg->SetBranchAddress("fTreeVariableCascRadius",                     &fBkgTreeVariableCascRadius);
  tBkg->SetBranchAddress("fTreeVariableV0Radius",                       &fBkgTreeVariableV0Radius);
  tBkg->SetBranchAddress("fTreeVariableLeastNbrClusters",               &fBkgTreeVariableLeastNbrClusters);
  tBkg->SetBranchAddress("fTreeVariableV0Lifetime",                     &fBkgTreeVariableV0Lifetime);
*/

  //other variables to be used when AC is performed
  /*
  tBkg->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger);		       
  tBkg->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz);			       
  tBkg->SetBranchAddress("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy);	
  //hh tBkg->SetBranchAddress("fTreeVariableChargeAssoc"             ,&fBkgTreeVariableChargeAssoc);		       
  //hh tBkg->SetBranchAddress("fTreeVariableAssocDCAz"                      ,&fBkgTreeVariableAssocDCAz);			       
  //hh tBkg->SetBranchAddress("fTreeVariableAssocDCAxy"                     ,&fBkgTreeVariableAssocDCAxy);			      
  tBkg->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger);			       
  tBkg->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0);
  tBkg->SetBranchAddress("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0);			       
  tBkg->SetBranchAddress("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0);			       
  tBkg->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta);		       
  tBkg->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi);		       
  tBkg->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta);		       
  */
  
  Int_t EntriesSign = 0; 
  Int_t EntriesBkg  = 0; 
  
  TFile *fout = new TFile(PathOut,"RECREATE");
  TDirectory  *dirSign= fout->mkdir("SE");
  //corr  TDirectory  *dirBkg= fout->mkdir("ME");

  /*-----------------------Pt Trigger  --------------------------- */
  /*
  TH1D *hSign_PtTrigger       = new TH1D("hSign_PtTrigger",       "hSign_PtTrigger",       300, 0, 30);
  hSign_PtTrigger->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtTrigger->GetXaxis()->SetTitleSize(0.05);
  hSign_PtTrigger->GetXaxis()->SetLabelSize(0.05);
  hSign_PtTrigger->GetYaxis()->SetLabelSize(0.05);
  TH1D *hBkg_PtTrigger         = new TH1D("hBkg_PtTrigger",        "hBkg_PtTrigger",        300, 0, 30);
  hBkg_PtTrigger->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hBkg_PtTrigger->GetXaxis()->SetTitleSize(0.05);
  hBkg_PtTrigger->GetXaxis()->SetLabelSize(0.05);
  hBkg_PtTrigger->GetYaxis()->SetLabelSize(0.05);
  */

  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hMassvsPt_SEbins[nummolt+1][numzeta][numPtTrigger];
  TH2D *hMassvsPt_SEbins_true[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1F *hXiAfterCuts= new TH1F ("hXiAfterCuts", "hXiAfterCuts", 18,0.5, 18.5);

  hXiAfterCuts->GetXaxis()->SetBinLabel(17,"All Xi");
  hXiAfterCuts->GetXaxis()->SetBinLabel(1,"NumberTPCClusters");
  hXiAfterCuts->GetXaxis()->SetBinLabel(2,"PID pos/neg");
  hXiAfterCuts->GetXaxis()->SetBinLabel(3,"PID Bach");
  hXiAfterCuts->GetXaxis()->SetBinLabel(4,"DCA pos/neg to PV");
  hXiAfterCuts->GetXaxis()->SetBinLabel(5,"DCA V0 to PV");
  hXiAfterCuts->GetXaxis()->SetBinLabel(6,"DCA V0 Daught.");
  hXiAfterCuts->GetXaxis()->SetBinLabel(7,"DCA Casc Daught.");
  hXiAfterCuts->GetXaxis()->SetBinLabel(8,"DCA Bach to PV");
  hXiAfterCuts->GetXaxis()->SetBinLabel(9,"Radius casc.");
  hXiAfterCuts->GetXaxis()->SetBinLabel(10,"Radius V0");
  hXiAfterCuts->GetXaxis()->SetBinLabel(11,"CosinePA V0");
  hXiAfterCuts->GetXaxis()->SetBinLabel(12,"CosinePA Casc");
  hXiAfterCuts->GetXaxis()->SetBinLabel(13,"Invmass Lambda");
  hXiAfterCuts->GetXaxis()->SetBinLabel(14,"V0 lifetime");
  hXiAfterCuts->GetXaxis()->SetBinLabel(15,"Casc lifetime");
  hXiAfterCuts->GetXaxis()->SetBinLabel(16,"Trues Casc");
  hXiAfterCuts->GetXaxis()->SetBinLabel(18,"Mass cuts");

  TH3F *hXiTrueSelected= new TH3F ("hXiTrueSelected", "hXiTrueSelected", 100,0, 100, 300, 0, 30, 80, LimInfHisto[type], LimSupHisto[type]);
  hXiTrueSelected->SetTitle("p_{T} and mult.class distribution of selected true Xi");
  hXiTrueSelected->GetXaxis()->SetTitle("Multiplicity class");
  hXiTrueSelected->GetYaxis()->SetTitle("p_{T, #Xi}");
  hXiTrueSelected->GetZaxis()->SetTitle("Inv mass");

  TH1F *hXiAfterCutsBkg= new TH1F ("hXiAfterCutsBkg", "PDG code della mamma del bachelor per le non-Xi", 7200,-3600, 3600);


  TH3F *hXiFalseSelected = new TH3F ("hXiFalseSelected", "hXiFalseSelected", 100,0, 100, 300, 0, 30, 80, LimInfHisto[type], LimSupHisto[type]);
  hXiFalseSelected->SetTitle("p_{T} and mult.class distribution of selected non-Xi");
  hXiFalseSelected->GetXaxis()->SetTitle("Multiplicity class");
  hXiFalseSelected->GetYaxis()->SetTitle("p_{T, #Xi}");
  hXiFalseSelected->GetZaxis()->SetTitle("Inv mass");

  TH1F *hXiFalseSelected1D = new TH1F ("hXiFalseSelected1D", "hXiFalseSelected1D",80, LimInfHisto[type], LimSupHisto[type]);

  TH3F *hXiFalseSelectedV0Particle = new TH3F ("hXiFalseSelectedV0Particle", "Distribuzione per candidate non-Xi con pos e neg da stessa mamma", 100,0, 100, 300, 0, 30, 80, LimInfHisto[type], LimSupHisto[type]);
  hXiFalseSelectedV0Particle->SetTitle("p_{T} and mult.class distribution of selected non-Xi");
  hXiFalseSelectedV0Particle->GetXaxis()->SetTitle("Multiplicity class");
  hXiFalseSelectedV0Particle->GetYaxis()->SetTitle("p_{T, #Xi}");
  hXiFalseSelectedV0Particle->GetZaxis()->SetTitle("Inv mass");

  TH1F *hXiFalseSelectedV0Particle1D = new TH1F ("hXiFalseSelectedV0Particle1D", "hXiFalseSelectedV0Particle1D",80, LimInfHisto[type], LimSupHisto[type]);

  TH3F *hXiSelected= new TH3F ("hXiSelected", "hXiSelected", 100,0, 100, 300, 0, 30, 80, LimInfHisto[type], LimSupHisto[type]);
  hXiSelected->SetTitle("p_{T} and mult.class distribution of selected true Xi");
  hXiSelected->GetXaxis()->SetTitle("Multiplicity class");
  hXiSelected->GetYaxis()->SetTitle("p_{T, #Xi}");
  hXiSelected->GetZaxis()->SetTitle("Inv mass");


  //Form("hMassvsPt_"+tipo[type]+"_%i",molt
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	hMassvsPt_SEbins[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i", m), 100,  LimInfHisto[type], LimSupHisto[type],  100, 0, 10);
	hMassvsPt_SEbins_true[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m), 100,  LimInfHisto[type], LimSupHisto[type], 100, 0, 10);
	for(Int_t v=0; v<numPtV0; v++){
	  nameSE[m][z][v][tr]="SE_";
	  namemassSE[m][z][v][tr]="InvMassSE_";
	  //nameSE[m][z][v][tr]+=Form("m%i_z%i_v%i_tr%i",m,z,v,tr);
	  nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_SEbins[m][z][v][tr]= new TH1D(namemassSE[m][z][v][tr], namemassSE[m][z][v][tr],100,  LimInfHisto[type], LimSupHisto[type]);
	  //	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr], nameSE[m][z][v][tr],   50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  // hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	}
      }
    }
  }

  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassME[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hMassvsPt_MEbins[nummolt+1][numzeta][numPtTrigger];
  TH2D *hMassvsPt_MEbins_true[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	hMassvsPt_MEbins[m][z][tr]= new TH2D(Form("ME_hMassvsPt_"+tipo[type]+"_%i", m),Form("ME_hMassvsPt_"+tipo[type]+"_%i", m),100, 0.4, 0.6,100, 0, 10);
	hMassvsPt_MEbins_true[m][z][tr]= new TH2D(Form("ME_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("ME_hMassvsPt_"+tipo[type]+"_%i_true", m),100, 0.4, 0.6,100, 0, 10);
	for(Int_t v=0; v<numPtV0; v++){
	  nameME[m][z][v][tr]="ME_";
	  namemassME[m][z][v][tr]="InvMassME_";
	  nameME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_MEbins[m][z][v][tr]= new TH1D(namemassME[m][z][v][tr], namemassME[m][z][v][tr],   100, 0.4, 0.6);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]= new TH2D(nameME[m][z][v][tr], nameME[m][z][v][tr],  50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);
	
	}
      }
    }
  }

  EntriesSign =  tSign->GetEntries();
  //  EntriesBkg  =  tBkg ->GetEntries();
    
  Bool_t MoltSel=kFALSE; 

  dirSign->cd();
  // cout << "  entries Sign: " << EntriesSign<<endl;
  for(Int_t k = 0; k<EntriesSign; k++){
    tSign->GetEntry(k);  
    if(isMC==0 || (isMC==1 && isEfficiency==1)){

      /*
      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMin) continue;
      CounterSignPairsAfterPtMinCut++;
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
      */
      //******************* some other cuts for sys studies**************************
      if (!ishhCorr){
	if ((type==0 || type==2) && fSignTreeVariableChargeCasc==1) continue; //type zero seleziona solo cascades negative
	if ((type==1 || type==3) && fSignTreeVariableChargeCasc==-1) continue; //type uno seleziona solo cascades positive
	hXiAfterCuts->Fill(17);
	if (fSignTreeVariableLeastNbrClusters<70) continue;
	hXiAfterCuts->Fill(1);
	if (fSignTreeVariableChargeCasc==1){
	  if (fSignTreeVariablePosNSigmaPion>4) continue;
	  if (fSignTreeVariableNegNSigmaProton>4) continue;
	}
	if (fSignTreeVariableChargeCasc==-1){
	  if (fSignTreeVariablePosNSigmaProton>4) continue;
	  if (fSignTreeVariableNegNSigmaPion>4) continue;
	}
	hXiAfterCuts->Fill(2);
	
	if ((type==0 || type==1) && fSignTreeVariableBachNSigmaPion>4) continue;
	if ((type==2 || type==3) && fSignTreeVariableBachNSigmaKaon>4) continue;
	
	hXiAfterCuts->Fill(3);
	if (fSignTreeVariableChargeCasc==1){
	  if (fSignTreeVariableDcaPosToPrimVertex<DcaMesonToPrimVertex[type]) continue; //mesone  def 0.04, 0.30 tight
	  if (fSignTreeVariableDcaNegToPrimVertex<DcaBaryonToPrimVertex[type]) continue; //barione def 0.03, 0.11 tight
	}
	
	if (fSignTreeVariableChargeCasc==-1){
	  if (fSignTreeVariableDcaPosToPrimVertex<DcaBaryonToPrimVertex[type]) continue; //barione
	  if (fSignTreeVariableDcaNegToPrimVertex<DcaMesonToPrimVertex[type]) continue; //mesone
	}
	
	hXiAfterCuts->Fill(4);
	// if (fSignTreeVariableDcaXiToPrimVertex) continue;
	// if (fSignTreeVariableXYDcaXiToPrimVertex) continue;
	//  if (fSignTreeVariableZDcaXiToPrimVertex) continue;
	if (fSignTreeVariableDcaV0ToPrimVertex<DcaV0ToPrimVertex[type]) continue; 
	hXiAfterCuts->Fill(5);
	if (fSignTreeVariableDcaV0Daughters>DcaV0Daughters[type]) continue; 
	hXiAfterCuts->Fill(6);
	if (fSignTreeVariableDcaCascDaughters>DcaCascDaughters[type]) continue; 
	hXiAfterCuts->Fill(7);
	if (fSignTreeVariableDcaBachToPrimVertex<DcaBachToPrimVertex[type]) continue; 
	hXiAfterCuts->Fill(8);
	
	if(fSignTreeVariableCascRadius<0.6 && (type ==0 || type==1)) continue;
	if(fSignTreeVariableCascRadius<0.5 && (type ==2 || type==3)) continue;
	hXiAfterCuts->Fill(9);
	if(fSignTreeVariableV0Radius<1.2  && (type ==0 || type==1)) continue;
	if(fSignTreeVariableV0Radius<1.1  && (type ==2 || type==3)) continue;
	hXiAfterCuts->Fill(10);

	//  if (fSignTreeVariableV0CosineOfPointingAngle) continue;
	if (fSignTreeVariableV0CosineOfPointingAngleSpecial<V0CosineOfPointingAngleSpecial[type]) continue; 
	hXiAfterCuts->Fill(11);
	if (fSignTreeVariableCascCosineOfPointingAngle<CascCosineOfPointingAngle[type]) continue; 
	hXiAfterCuts->Fill(12);
  
	if (TMath::Abs(fSignTreeVariableInvMassLambda - massLambda) > InvMassLambda[type]) continue;
	if (TMath::Abs(fSignTreeVariableInvMassXi - massXi) < 0.008   && (type ==2 || type==3)) continue;
	if (TMath::Abs(fSignTreeVariableInvMassOmega - massOmega) < 0.005   && (type ==0 || type==1)) continue;
	hXiAfterCuts->Fill(13);
	//  if (fSignTreeVariableRapXi) continue;
	//  if (fSignTreeVariableRapOmega) continue;
	if (fSignTreeVariableV0Lifetime>30) continue;
	hXiAfterCuts->Fill(14);
	if (fSignTreeVariablectau> 3*ctauCasc[type]) continue;
	hXiAfterCuts->Fill(15);
	if (TMath::Abs(fSignTreeVariableRapXi)>= RapCasc[rap]) continue;
	if ((fSignTreeVariableInvMassXi< 1.28 || fSignTreeVariableInvMassXi> 1.36) && (type==0 || type==1)) continue;
	if ((fSignTreeVariableInvMassOmega< 1.62 || fSignTreeVariableInvMassOmega> 1.72) && (type==2 || type==3)) continue;

	hXiAfterCuts->Fill(18);

	if (type==0 || type==1)  hXiSelected->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassXi);
	if (type==2 || type==3)  hXiSelected->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassOmega);

	if (fSignTreeVariablePDGCode==CascadesPDGCode[type]){
	  if (type==0 || type==1)  hXiTrueSelected->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassXi);
	  if (type==2 || type==3)  hXiTrueSelected->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassOmega);
	}
	if (fSignTreeVariablePDGCode!=CascadesPDGCode[type]){
	  if (type==0 || type==1)  {
	    hXiFalseSelected->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassXi);
	    hXiFalseSelected1D->Fill(fSignTreeVariableInvMassXi);
	  }
	  if (type==2 || type==3)  {
	    hXiFalseSelected->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassOmega);
	    hXiFalseSelected1D->Fill(fSignTreeVariableInvMassOmega);
	  }
	}

	if (fSignTreeVariablePDGCode!=CascadesPDGCode[type] && fSignTreeVariablePDGCodeLambda!=0){
	  if (type==0 || type==1) {
	    hXiFalseSelectedV0Particle->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassXi);
	    hXiFalseSelectedV0Particle1D->Fill(fSignTreeVariableInvMassXi);
	  }
	  if (type==2 || type==3)  {
	    hXiFalseSelectedV0Particle->Fill(fSignTreeVariableMultiplicity,fSignTreeVariablePtCasc,  fSignTreeVariableInvMassOmega);
	    hXiFalseSelectedV0Particle1D->Fill(fSignTreeVariableInvMassOmega);
	  }
	}

	if (fSignTreeVariablePDGCode==CascadesPDGCode[type] ){
	  hXiAfterCuts->Fill(16);
	}
	if (fSignTreeVariablePDGCode!=CascadesPDGCode[type] && 1.315< fSignTreeVariableInvMassXi &&    fSignTreeVariableInvMassXi< 1.327 ){
	  hXiAfterCutsBkg->Fill(fSignTreeVariablePDGCode);
	}
      }
      

      /*
	if (ishhCorr){
	if(sysV0==0){
	if(TMath::Abs(fSignTreeVariableAssocDCAz)>1) continue;
	}
	if(sysV0==1){
	if(TMath::Abs(fSignTreeVariableAssocDCAz)>2) continue;
	}
	if(sysV0==2){
	if(TMath::Abs(fSignTreeVariableAssocDCAz)>0.5) continue;
	}
	}
      */

    }
  
    //**********************************************************************************************

    // if(TMath::Abs((fSignTreeVariableInvMassK0s - massK0s))> 0.01){
    //   hSign_PtTrigger->Fill(fSignTreeVariablePtTrigger);
    // }
   
    //    cout << "*********************"<< endl;
//corr    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
//corr    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();

    for(Int_t m=0; m<nummolt+1; m++){
      if (m < nummolt){
	MoltSel = (fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]);
      }
      else if(m==nummolt) MoltSel=kTRUE;
      //      cout << m << "  " << nummolt << "  " << MoltSel<< endl;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  if(MoltSel){
	    if (type==0 || type==1){
	      if (type==0 && fSignTreeVariableChargeCasc==-1)	    hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMassXi, fSignTreeVariablePtCasc); 
	      if (type==1 && fSignTreeVariableChargeCasc==1) hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMassXi, fSignTreeVariablePtCasc); 
	    if(fSignTreeVariablePDGCode==CascadesPDGCode[type]) hMassvsPt_SEbins_true[m][z][tr]->Fill(fSignTreeVariableInvMassXi, fSignTreeVariablePtCasc); 
	  }
	    if (type==2 || type==3){
	      if (type==2 && fSignTreeVariableChargeCasc==-1)	      hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMassOmega, fSignTreeVariablePtCasc); 
	      if (type==3 && fSignTreeVariableChargeCasc==1)	      hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMassOmega, fSignTreeVariablePtCasc); 
	      if(fSignTreeVariablePDGCode==CascadesPDGCode[type]) hMassvsPt_SEbins_true[m][z][tr]->Fill(fSignTreeVariableInvMassOmega, fSignTreeVariablePtCasc); 
	    }
	  }
	  for(Int_t v=0; v<numPtV0; v++){
	    //corr if(MoltSel && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	    if(MoltSel  && fSignTreeVariablePtCasc>=NPtV0[v]&& fSignTreeVariablePtCasc<NPtV0[v+1]){
	      //corr	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	    if (type==0 || type==1)	      hInvMassK0Short_SEbins[m][z][v][tr]->Fill(fSignTreeVariableInvMassXi);
	    if (type==2 || type==3)	      hInvMassK0Short_SEbins[m][z][v][tr]->Fill(fSignTreeVariableInvMassOmega);
	    }
	  }
	}
      }
    }
  } //end loop on entries

  /* bakg if AC is performed
  cout << MoltSel<< endl;
  dirBkg->cd();
  for(Int_t k = 0; k<EntriesBkg; k++){
    // for(Int_t k = 0; k<1; k++){
    tBkg->GetEntry(k);
    if(isMC==0 || (isMC==1 && isEfficiency==1)){
      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMin) continue;
      CounterBkgPairsAfterPtMinCut++;
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
      if (!ishhCorr){
	if(sysV0!=5){
	  if(TMath::Abs((fBkgTreeVariableInvMassLambda - massLambda))<= 0.005) continue;
	  if(TMath::Abs((fBkgTreeVariableInvMassAntiLambda - massLambda))<= 0.005) continue;
	  if(sysV0==1){   
	    if(fBkgTreeVariableV0CosineOfPointingAngle <= 0.997) continue;
	  }
	  if(sysV0==2){   
	    if(fBkgTreeVariablectau >=8 )continue;
	  }
	  if(sysV0==3){   
	    if(fBkgTreeVariableRapK0Short >= 0.5)continue;
	  }
	  if(sysV0==4){   
	    if(TMath::Abs((fBkgTreeVariableInvMassLambda - massLambda))<= 0.010) continue;
	    if(TMath::Abs((fBkgTreeVariableInvMassAntiLambda - massLambda))<= 0.010) continue;
	  }
	  if(sysV0==6){   
	    if(fBkgTreeVariableDcaV0ToPrimVertex>=0.3 )continue;
	  }
   
	}
      }

      if (ishhCorr){
	if (ishhCorr){
	  if(sysV0==0){
	    if(TMath::Abs(fBkgTreeVariableAssocDCAz)>1) continue;
	  }
	  if(sysV0==1){
	    if(TMath::Abs(fBkgTreeVariableAssocDCAz)>2) continue;
	  }
	  if(sysV0==2){
	    if(TMath::Abs(fBkgTreeVariableAssocDCAz)>0.5) continue;
	  }
	}

      }
    }
  
    //**********************************************************************************************
    // if(TMath::Abs((fBkgTreeVariableInvMassK0s - massK0s))> 0.01){
    // hBkg_PtTrigger->Fill(fBkgTreeVariablePtTrigger);
    // }
    if (fBkgTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fBkgTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fBkgTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fBkgTreeVariableDeltaPhi += 2.0*TMath::Pi();
 
    for(Int_t m=0; m<nummolt+1; m++){
      if (m < nummolt){
	MoltSel = (fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1]);
      }
      else  if(m==nummolt) MoltSel=kTRUE;
  
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  if(MoltSel){
	    hMassvsPt_MEbins[m][z][tr]->Fill(fBkgTreeVariableInvMassK0s, fBkgTreeVariablePtV0);
	    if(fBkgTreeVariablePDGCode==310) 	    hMassvsPt_MEbins_true[m][z][tr]->Fill(fBkgTreeVariableInvMassK0s, fBkgTreeVariablePtV0);
	  }
	  for(Int_t v=0; v<numPtV0; v++){
	    if(MoltSel && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      hInvMassK0Short_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableInvMassK0s);
	      hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	    }
	  }
XS	}
      }
    }
  }
  */
   
  fout->Write();
//corr  cout << "Pt Min delle particelle trigger " << PtTrigMin<< endl;
//corr  cout << "signal pairs trigger-associated after Pt min cut " << 	CounterSignPairsAfterPtMinCut << endl;
//corr  cout << "bkg pairs trigger-associated after Pt min cut " << 	CounterBkgPairsAfterPtMinCut << endl;
  cout << "partendo dal file " << PathIn << " ho creato il file " << PathOut<< endl;
  cout << CascadesPDGCode[type]<< endl;
}

