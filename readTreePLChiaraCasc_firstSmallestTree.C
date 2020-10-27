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
void readTreePLChiaraCasc_firstSmallestTree(Int_t type=4 /*type = 0 for XiMinus, =1, for XiPlus, =2 for OmegaMinus, =3 for OmegaPlus */,Bool_t SkipAssoc=1 ,Int_t israp=0, Bool_t ishhCorr=0, Float_t PtTrigMin=3, Float_t ptjmax=15, Int_t sysV0=0, bool isMC = 0,Bool_t isEfficiency=1,Int_t sysTrigger=0,			    TString year=/*"AllMC_hXi"/*"2015f"/*"2015g3b1"*/"Run2DataRed_MECorr_hXi"/*"2016kl_hXi"/*"2016k_hXi_MECorr_25runs"/"1617GP_hXi"/*"2018f1_extra_hXi_CommonParton"*/, TString year0="2016", TString Path1 ="", Bool_t CommonParton=0, Int_t MultBinning=0)
{

  //rap=0 : no rapidity window chsen for cascade, |Eta| < 0.8; rap=1 |y| < 0.5
  if (ishhCorr) {
    cout << "This macro should not be run is hh correlation is studied; go directly to readTreePLChiara_second " << endl;
    return;
  }

  if (israp>1) return;
  if (sysV0>6) return;
  if (sysV0>2 && ishhCorr) return;

  //I set the loostest topological selections (an also the loosest DCAz for trigger)
  Float_t DCAXiDaughters=1.8;
  Float_t CosinePAngleXiToPV=0.95;
  Float_t CosinePAngleV0ToXi=0.95;
  Float_t InvMassLambdaWindow = 0.009;
  Float_t ctau = 5* 4.91;
  Float_t DCAzTrigger = 2;
  //////////////////////////////////////////////////////////////////////////////////

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=8; //14;//was 8 for my analysis
  const Int_t numPtTrigger=1;
  const Int_t numtipo=6;
  TString tipo[numtipo]={"XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega"};
  TString Srap[2] = {"_Eta0.8", "_y0.5"};
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
 
  //  PathIn+=Path1;
  PathIn+=".root";
  PathOut+=Path1;
  PathOut+="_";
  PathOut +=tipo[type];
  PathOut +=Srap[israp];
  if (!SkipAssoc)  PathOut +="_AllAssoc";
  PathOut +=Form("_MassDistr_SysT%i_SysV0%i_index%i_PtMin%.1f",sysTrigger, sysV0,0, PtTrigMin); 
  PathOut+= ".root";

  TString dirinputtype[6] = {"Xi", "Xi", "Omega", "Omega", "Xi", "Omega"};
  TFile *fin = new TFile(PathIn);
  if (!fin) {cout << "file input not available " ; return;}
  TDirectoryFile *d = (TDirectoryFile*)fin->Get("MyTask"+dirinputtype[type]);
  if (!d)  {cout << "dir input not available " ; return;}

  TTree *tSign = (TTree *)d->Get("fSignalTree");
  TTree *tBkg  = (TTree *)d->Get("fBkgTree");
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Bool_t MassInPeakCasc=0;
  Double_t ctauK0s = 2.6844;

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
  Bool_t       fSignTreeVariableSkipAssoc;                        
  Bool_t       fSignTreeVariableIsCommonParton;                        
  Int_t        fSignTreeVariablePdgCommonParton;                        


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
  Bool_t       fBkgTreeVariableSkipAssoc;                        
  Bool_t       fBkgTreeVariableIsCommonParton;                        
  Int_t        fBkgTreeVariablePdgCommonParton;                        

  Int_t CounterSignPairsAfterPtMinCut=0; 
  Int_t CounterSignPairsIntermediate=0; 
  Int_t CounterBkgPairsAfterPtMinCut=0; 
  Int_t CounterBkgPairsIntermediate=0; 
  Int_t TrueCounterSignPairsAfterPtMinCut=0; 
  Int_t TrueCounterBkgPairsAfterPtMinCut=0; 
  Int_t CounterSignPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t CounterBkgPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t TrueCounterSignPairsAfterPtMinCutMult[nummolt+1]={0}; 
  Int_t TrueCounterBkgPairsAfterPtMinCutMult[nummolt+1]={0}; 

  Float_t ctauCasc[numtipo] = {4.91,4.91,  2.461, 2.461, 4.91, 2.461}; //cm , average ctau of Xi and Omega     
  Float_t PDGCode[numtipo-2] = {3312, -3312, 3334, -3334};
  Float_t LimInfMass[numtipo]= {1.29, 1.29, 1.5, 1.5, 1.29, 1.5};
  Float_t LimSupMass[numtipo]= {1.35, 1.35, 1.8, 1.8, 1.35, 1.5};

  TString Smolt0[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  Double_t Nmolt0[nummolt+1]={0, 5, 10, 30, 50, 100};
  TString Smolt1[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "all"};
  Double_t Nmolt1[nummolt+1]={0, 1, 5, 15, 30, 100};
  TString Smolt2[nummolt+1]={"0-2", "2-7", "7-15", "15-30", "30-100", "all"};
  Double_t Nmolt2[nummolt+1]={0, 2, 7, 15, 30, 100};

  TString Smolt[nummolt+1]={"0-1", "1-5", "5-15", "15-30", "30-100", "all"};
  Double_t Nmolt[nummolt+1]={0, 1, 5, 15, 30, 100};
  
  for (Int_t m=0; m<nummolt+1; m++){
    if (MultBinning==0){
      Smolt[m] = Smolt0[m];
      Nmolt[m] = Nmolt0[m];
    }
    else if (MultBinning==1){
      Smolt[m] = Smolt1[m];
      Nmolt[m] = Nmolt1[m];
    }
    else if (MultBinning==2){
      Smolt[m] = Smolt2[m];
      Nmolt[m] = Nmolt2[m];
    }

  }

  //  TString Smolt[nummolt+1]={"0-0", "0-10", "10-30", "30-50", "50-100", "all"};
  //  Double_t Nmolt[nummolt+1]={0, 0, 10, 30, 50, 100};

  TString SPtV0[numPtV0]={"0-0.5","0.5-1.0", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  //  TString SPtV0[numPtV0]={"0-0", "0-0.5","0.5-1.0", "1-1.5", "1.5-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,0.5, 1,1.5,2,2.5,3,4,8};
  //Double_t NPtV0[numPtV0+1]={0,0,0.5, 1,1.5,2,3,4,8};


  //  TString SPtV0[numPtV0]={"0-0.6", "0.6-1.0", "1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0", "2.0-2.2", "2.2-2.5", "2.5-2.9", "2.9-3.4", "3.4-4", "4-5", "5-6.5"};
  //  Double_t NPtV0[numPtV0+1]={0,0.6, 1,1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4,5,6.5};
  Double_t NPtTrigger[numPtTrigger+1]={PtTrigMin,ptjmax};

  
  //what is the fraction of AC events in each multiplicity class?
  TList *d1 = (TList*)d->Get("MyOutputContainer");
  if (!d1) return;
  TH1F* fHistEventMult=(TH1F*)  d1->FindObject("fHistEventMult");
  if (!fHistEventMult) cout << "no info about total number of INT7 events analyzed" << endl;
  cout <<" bin label: " << fHistEventMult->GetXaxis()->GetBinLabel(7) << endl;
  Double_t TotEvtINT7 = fHistEventMult->GetBinContent(7);

  TH2F* hMultiplicity2D=(TH2F*)  d1->FindObject("fHistPtMaxvsMult");
  TH2F* hMultiplicity2DBefAll=(TH2F*)  d1->FindObject("fHistPtMaxvsMultBefAll");
  TH1F* hMultiplicity;
  TH1F* hMultiplicityBefAll;
  TH2F* hMultvsNumberAssoc=(TH2F*)  d1->FindObject("fHistMultvsV0");
  TH1F*  hMultvsNumberAssoc_Proj[nummolt];
  if (!hMultiplicity2D) cout << " no info about multiplicity distribution of AC events available " << endl;
  if (!hMultvsNumberAssoc) cout << " no info about multiplicity distribution of AC events available " << endl;
  Float_t ACcounter[nummolt+1];

  if (hMultiplicity2D && hMultvsNumberAssoc && hMultiplicity2DBefAll){
    hMultiplicity=(TH1F*)  hMultiplicity2D->ProjectionY("fHistPtMaxvsMult1D",     hMultiplicity2D->GetXaxis()->FindBin(PtTrigMin+0.0001),  hMultiplicity2D->GetXaxis()->FindBin(ptjmax-0.0001));
    hMultiplicityBefAll=(TH1F*)  hMultiplicity2DBefAll->ProjectionY("fHistPtMaxvsMult1DBefAll",     hMultiplicity2DBefAll->GetXaxis()->FindBin(PtTrigMin+0.0001),  hMultiplicity2DBefAll->GetXaxis()->FindBin(ptjmax-0.0001));
    ACcounter[5]= hMultiplicity->GetEntries();

    cout <<"total number of events with NT>0  " << hMultiplicityBefAll->GetEntries() << " " <<   (Float_t)hMultiplicityBefAll->GetEntries()/TotEvtINT7 << endl;
    cout <<"\ntotal number of events used in the AC (no selections on topo var and on skipAssoc!) " << hMultiplicity->GetEntries() << endl;

    for (Int_t m=0; m< nummolt; m++){ 
      hMultvsNumberAssoc_Proj[m] = (TH1F*)       hMultvsNumberAssoc->ProjectionX(Form("hMultvsNumberAssoc_%i", m), hMultvsNumberAssoc->GetYaxis()->FindBin(Nmolt[m]+0.001), hMultvsNumberAssoc->GetYaxis()->FindBin(Nmolt[m+1]-0.001));

      cout << " m " << m << endl;
      ACcounter[m] =0;
      for (Int_t b= hMultiplicity->GetXaxis()->FindBin(Nmolt[m]+0.001); b <=  hMultiplicity->GetXaxis()->FindBin(Nmolt[m+1]-0.001); b++){
	ACcounter[m] +=  hMultiplicity->GetBinContent(b);
      }
      cout << "fraction of events in mult bin (for ptTrig> PtTrigMin) " << Smolt[m] << ": " << ACcounter[m]/ACcounter[5] <<  " ~average V0 number (for pTTrig> 0.150 GeV usually) " <<     hMultvsNumberAssoc_Proj[m]->GetMean() <<endl;
    }
    
  }

  
  //Signal
  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);		       
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);		       
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);		       
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);		       
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);			       
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);			       
  tSign->SetBranchAddress("fTreeVariableChargeAssoc"             ,&fSignTreeVariableChargeAssoc);		       	   
  tSign->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger);			       
  tSign->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0);			       
  tSign->SetBranchAddress("fTreeVariableRapAssoc"                ,&fSignTreeVariableRapAssoc);		       
  tSign->SetBranchAddress("fTreeVariableDcaXiDaughters"         ,&fSignTreeVariableDcaXiDaughters);	       
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fSignTreeVariableV0CosineOfPointingAngle);	       
  tSign->SetBranchAddress("fTreeVariableXiCosineOfPointingAngle"   ,&fSignTreeVariableXiCosineOfPointingAngle);	       
  tSign->SetBranchAddress("fTreeVariablectau"                      ,&fSignTreeVariablectau);			       
  tSign->SetBranchAddress("fTreeVariablePtV0"                      ,&fSignTreeVariablePtV0);			       
  tSign->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fSignTreeVariableInvMassLambda);		       
  tSign->SetBranchAddress("fTreeVariableInvMassXi"             ,&fSignTreeVariableInvMassXi);		       
  tSign->SetBranchAddress("fTreeVariableInvMassOmega"             ,&fSignTreeVariableInvMassOmega);		       
  tSign->SetBranchAddress("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0);			       
  tSign->SetBranchAddress("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0);			       
  tSign->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta);		       
  tSign->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi);		       
  tSign->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta);		       
  tSign->SetBranchAddress("fTreeVariableMultiplicity"              ,&fSignTreeVariableMultiplicity);		       
  tSign->SetBranchAddress("fTreeVariableZvertex"                   ,&fSignTreeVariableZvertex);                        
  tSign->SetBranchAddress("fTreeVariablePDGCodeTrigger"                   ,&fSignTreeVariablePDGCodeTrigger);                         
  tSign->SetBranchAddress("fTreeVariablePDGCodeAssoc"                   ,&fSignTreeVariablePDGCodeAssoc);                         
  tSign->SetBranchAddress("fTreeVariableSkipAssoc"                   ,&fSignTreeVariableSkipAssoc);                         
  tSign->SetBranchAddress("fTreeVariableIsCommonParton"                   ,&fSignTreeVariableIsCommonParton);                        
  tSign->SetBranchAddress("fTreeVariablePdgCommonParton"                   ,&fSignTreeVariablePdgCommonParton);                        



  //BackGround
  tBkg->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger);		       
  tBkg->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz);			       
  tBkg->SetBranchAddress("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy);	 tBkg->SetBranchAddress("fTreeVariableChargeAssoc"             ,&fBkgTreeVariableChargeAssoc);		       
  tBkg->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger);			       
  tBkg->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0);
  tBkg->SetBranchAddress("fTreeVariableRapAssoc"                ,&fBkgTreeVariableRapAssoc);		       
  tBkg->SetBranchAddress("fTreeVariableDcaXiDaughters"         ,&fBkgTreeVariableDcaXiDaughters);	       
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fBkgTreeVariableV0CosineOfPointingAngle);	        
  tBkg->SetBranchAddress("fTreeVariableXiCosineOfPointingAngle"   ,&fBkgTreeVariableXiCosineOfPointingAngle);	        
  tBkg->SetBranchAddress("fTreeVariablectau"                      ,&fBkgTreeVariablectau);			       
  tBkg->SetBranchAddress("fTreeVariablePtV0"                      ,&fBkgTreeVariablePtV0);			       
  tBkg->SetBranchAddress("fTreeVariableInvMassXi"                ,&fBkgTreeVariableInvMassXi);		       
  tBkg->SetBranchAddress("fTreeVariableInvMassOmega"                ,&fBkgTreeVariableInvMassOmega);		       
  tBkg->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fBkgTreeVariableInvMassLambda);		       
  tBkg->SetBranchAddress("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0);			       
  tBkg->SetBranchAddress("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0);			       
  tBkg->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta);		       
  tBkg->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi);		       
  tBkg->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta);		       
  tBkg->SetBranchAddress("fTreeVariableMultiplicity"              ,&fBkgTreeVariableMultiplicity);		       
  tBkg->SetBranchAddress("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex);
  tBkg->SetBranchAddress("fTreeVariablePDGCodeTrigger"                   ,&fBkgTreeVariablePDGCodeTrigger);                        
  tBkg->SetBranchAddress("fTreeVariablePDGCodeAssoc"                   ,&fBkgTreeVariablePDGCodeAssoc);                        
  tBkg->SetBranchAddress("fTreeVariableSkipAssoc"                   ,&fBkgTreeVariableSkipAssoc);                        
  tBkg->SetBranchAddress("fTreeVariableIsCommonParton"                   ,&fBkgTreeVariableIsCommonParton);                        
  tBkg->SetBranchAddress("fTreeVariablePdgCommonParton"                   ,&fBkgTreeVariablePdgCommonParton);                        

  Int_t EntriesSign = 0; 
  Int_t EntriesBkg  = 0; 


  
  TFile *fout = new TFile(PathOut,"RECREATE");
  TTree * tSignO = new TTree("tSignO", "tSignO");
  //Signal
  tSignO->Branch("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger,		       	"fSignTreeVariablePtTrigger/D");		       	
  tSignO->Branch("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger,		"fSignTreeVariableChargeTrigger/I");		       
  tSignO->Branch("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger,		       	"fSignTreeVariableEtaTrigger/D");		       	
  tSignO->Branch("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger,		       	"fSignTreeVariablePhiTrigger/D");		       	
  tSignO->Branch("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz,			"fSignTreeVariableDCAz/D");			       
  tSignO->Branch("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy,			"fSignTreeVariableDCAxy/D");			
  tSignO->Branch("fTreeVariableChargeAssoc"               ,&fSignTreeVariableChargeAssoc,		"fSignTreeVariableChargeAssoc/I");		       
  tSignO->Branch("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger,		"fSignTreeVariableisPrimaryTrigger/I");		    
  tSignO->Branch("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0,		"fSignTreeVariableisPrimaryV0/I");		
  tSignO->Branch("fTreeVariableRapAssoc"                  ,&fSignTreeVariableRapAssoc,		       	"fSignTreeVariableRapAssoc/D");		       	
  tSignO->Branch("fTreeVariableDcaXiDaughters"            ,&fSignTreeVariableDcaXiDaughters,	       	"fSignTreeVariableDcaXiDaughters/D");	       	
  tSignO->Branch("fTreeVariableV0CosineOfPointingAngle"   ,&fSignTreeVariableV0CosineOfPointingAngle,	"fSignTreeVariableV0CosineOfPointingAngle/D");	        
  tSignO->Branch("fTreeVariableXiCosineOfPointingAngle"   ,&fSignTreeVariableXiCosineOfPointingAngle,	"fSignTreeVariableXiCosineOfPointingAngle/D");	        
  tSignO->Branch("fTreeVariablectau"                      ,&fSignTreeVariablectau,			"fSignTreeVariablectau/D");			       
  tSignO->Branch("fTreeVariablePtV0"                      ,&fSignTreeVariablePtV0,			"fSignTreeVariablePtV0/D");			       
  tSignO->Branch("fTreeVariableInvMassXi"                 ,&fSignTreeVariableInvMassXi,		       	"fSignTreeVariableInvMassXi/D");		       	
  tSignO->Branch("fTreeVariableInvMassOmega"              ,&fSignTreeVariableInvMassOmega,		"fSignTreeVariableInvMassOmega/D");		       
  tSignO->Branch("fTreeVariableInvMassLambda"             ,&fSignTreeVariableInvMassLambda,		"fSignTreeVariableInvMassLambda/D");		       
  tSignO->Branch("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0,			"fSignTreeVariableEtaV0/D");			       
  tSignO->Branch("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0,			"fSignTreeVariablePhiV0/D");			       
  tSignO->Branch("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta,		       	"fSignTreeVariableDeltaEta/D");		       	
  tSignO->Branch("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi,		       	"fSignTreeVariableDeltaPhi/D");		       	
  tSignO->Branch("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta,		       	"fSignTreeVariableDeltaTheta/D");		       	
  tSignO->Branch("fTreeVariableMultiplicity"              ,&fSignTreeVariableMultiplicity,		"fSignTreeVariableMultiplicity/D");		       
  tSignO->Branch("fTreeVariableZvertex"                   ,&fSignTreeVariableZvertex,			"fSignTreeVariableZvertex/D");			
  tSignO->Branch("fTreeVariablePDGCodeTrigger"            ,&fSignTreeVariablePDGCodeTrigger,            	"fSignTreeVariablePDGCodeTrigger/I");            	
  tSignO->Branch("fTreeVariablePDGCodeAssoc"              ,&fSignTreeVariablePDGCodeAssoc,               "fSignTreeVariablePDGCodeAssoc/I");                 
  tSignO->Branch("fTreeVariableSkipAssoc"                 ,&fSignTreeVariableSkipAssoc,                  "fSignTreeVariableSkipAssoc/O");                        
  tSignO->Branch("fTreeVariableIsCommonParton"            ,&fSignTreeVariableIsCommonParton,             "fSignTreeVariableIsCommonParton/O");             
  tSignO->Branch("fTreeVariablePdgCommonParton"           ,&fSignTreeVariablePdgCommonParton,            "fSignTreeVariablePdgCommonParton/I");             


  TTree * tBkgO = new TTree("tBkgO", "tBkgO");
  //BackGround
  tBkgO->Branch("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger,		       	"fBkgTreeVariablePtTrigger/D");		       	
  tBkgO->Branch("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger,		"fBkgTreeVariableChargeTrigger/I");		       
  tBkgO->Branch("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger,		       	"fBkgTreeVariableEtaTrigger/D");		       	
  tBkgO->Branch("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger,		       	"fBkgTreeVariablePhiTrigger/D");		       	
  tBkgO->Branch("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz,			"fBkgTreeVariableDCAz/D");			       
  tBkgO->Branch("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy,			"fBkgTreeVariableDCAxy/D");			
  tBkgO->Branch("fTreeVariableChargeAssoc"               ,&fBkgTreeVariableChargeAssoc,		"fBkgTreeVariableChargeAssoc/I");		       
  tBkgO->Branch("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger,		"fBkgTreeVariableisPrimaryTrigger/I");		    
  tBkgO->Branch("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0,		"fBkgTreeVariableisPrimaryV0/I");		
  tBkgO->Branch("fTreeVariableRapAssoc"                  ,&fBkgTreeVariableRapAssoc,		       	"fBkgTreeVariableRapAssoc/D");		       	
  tBkgO->Branch("fTreeVariableDcaXiDaughters"            ,&fBkgTreeVariableDcaXiDaughters,	       	"fBkgTreeVariableDcaXiDaughters/D");	       	
  tBkgO->Branch("fTreeVariableV0CosineOfPointingAngle"   ,&fBkgTreeVariableV0CosineOfPointingAngle,	"fBkgTreeVariableV0CosineOfPointingAngle/D");	        
  tBkgO->Branch("fTreeVariableXiCosineOfPointingAngle"   ,&fBkgTreeVariableXiCosineOfPointingAngle,	"fBkgTreeVariableXiCosineOfPointingAngle/D");	        
  tBkgO->Branch("fTreeVariablectau"                      ,&fBkgTreeVariablectau,			"fBkgTreeVariablectau/D");			       
  tBkgO->Branch("fTreeVariablePtV0"                      ,&fBkgTreeVariablePtV0,			"fBkgTreeVariablePtV0/D");			       
  tBkgO->Branch("fTreeVariableInvMassXi"                 ,&fBkgTreeVariableInvMassXi,		       	"fBkgTreeVariableInvMassXi/D");		       	
  tBkgO->Branch("fTreeVariableInvMassOmega"              ,&fBkgTreeVariableInvMassOmega,		"fBkgTreeVariableInvMassOmega/D");		       
  tBkgO->Branch("fTreeVariableInvMassLambda"             ,&fBkgTreeVariableInvMassLambda,		"fBkgTreeVariableInvMassLambda/D");		       
  tBkgO->Branch("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0,			"fBkgTreeVariableEtaV0/D");			       
  tBkgO->Branch("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0,			"fBkgTreeVariablePhiV0/D");			       
  tBkgO->Branch("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta,		       	"fBkgTreeVariableDeltaEta/D");		       	
  tBkgO->Branch("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi,		       	"fBkgTreeVariableDeltaPhi/D");		       	
  tBkgO->Branch("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta,		       	"fBkgTreeVariableDeltaTheta/D");		       	
  tBkgO->Branch("fTreeVariableMultiplicity"              ,&fBkgTreeVariableMultiplicity,		"fBkgTreeVariableMultiplicity/D");		       
  tBkgO->Branch("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex,			"fBkgTreeVariableZvertex/D");			
  tBkgO->Branch("fTreeVariablePDGCodeTrigger"            ,&fBkgTreeVariablePDGCodeTrigger,            	"fBkgTreeVariablePDGCodeTrigger/I");            	
  tBkgO->Branch("fTreeVariablePDGCodeAssoc"              ,&fBkgTreeVariablePDGCodeAssoc,               "fBkgTreeVariablePDGCodeAssoc/I");                 
  tBkgO->Branch("fTreeVariableSkipAssoc"                 ,&fBkgTreeVariableSkipAssoc,                  "fBkgTreeVariableSkipAssoc/O");                        
  tBkgO->Branch("fTreeVariableIsCommonParton"            ,&fBkgTreeVariableIsCommonParton,             "fBkgTreeVariableIsCommonParton/O");             
  tBkgO->Branch("fTreeVariablePdgCommonParton"           ,&fBkgTreeVariablePdgCommonParton,            "fBkgTreeVariablePdgCommonParton/I");             


  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_BulkReg[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPPtInt_JetReg[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTruePtInt[nummolt+1][numzeta][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTruePtInt[nummolt+1][numzeta][numPtTrigger];

  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_BulkReg[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCP_JetReg[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_CPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_SEbins_DeltaPhiProj_NOCPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];

  TH1F * hMultiplicityDIstributionACEv = new TH1F ("hMultiplicityDIstributionACEv", "hMultiplicityDIstributionACEv", 100, 0,100);

  TH1F * HistoInfo = new TH1F("HistoInfo", "HistoInfo",20,0.5, 20.5);
  HistoInfo->GetXaxis()->SetBinLabel(1, "% Ev. NT>0");
  HistoInfo->GetXaxis()->SetBinLabel(2, "NV0/NInt7");
  HistoInfo->GetXaxis()->SetBinLabel(3, "<p_{T,Xi}>");
  HistoInfo->GetXaxis()->SetBinLabel(4,"<p_{T,Trig}>");
  HistoInfo->GetXaxis()->SetBinLabel(5,"SE pairs");
  HistoInfo->GetXaxis()->SetBinLabel(6,"ME pairs");
  for (Int_t m=0; m< nummolt; m++){
    HistoInfo->GetXaxis()->SetBinLabel(7+m,Form("SEpairs Mult distr m%i",m));
    HistoInfo->GetXaxis()->SetBinLabel(12+m, Form("NV0/ACev m%i",m));
  }

  //--------histograms to get some info from Xi candidate having a common origin to the trigger particle----------------
  TH1D *hSign_PtAssoc_CPTrue       = new TH1D("hSign_PtAssoc_CPTrue",       "hSign_PtAssoc_CPTrue",       300, 0, 30);
  hSign_PtAssoc_CPTrue->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssoc_CPTrue->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssoc_CPTrue->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssoc_CPTrue->GetYaxis()->SetLabelSize(0.05);

  TH1D *hSign_PtAssoc_NOCPTrue       = new TH1D("hSign_PtAssoc_NOCPTrue",       "hSign_PtAssoc_NOCPTrue",       300, 0, 30);
  hSign_PtAssoc_NOCPTrue->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssoc_NOCPTrue->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssoc_NOCPTrue->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssoc_NOCPTrue->GetYaxis()->SetLabelSize(0.05);

  //------------------Histograms os selected particles (V0) for future efficiency calculation ----------------
  TH3F*    fHistSelectedV0PtTMaxPhi=new TH3F(Form("fHistSelectedV0PtTMaxPhi_%i", sysV0), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
  fHistSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F*    fHistCPSelectedV0PtTMaxPhi=new TH3F(Form("fHistCPSelectedV0PtTMaxPhi_%i", sysV0), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
  fHistCPSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistCPSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F*    fHistNOCPSelectedV0PtTMaxPhi=new TH3F(Form("fHistNOCPSelectedV0PtTMaxPhi_%i", sysV0), "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
  fHistNOCPSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistNOCPSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F *    fHistSelectedV0PtTMaxEta=new TH3F(Form("fHistSelectedV0PtTMaxEta_%i", sysV0), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles (Casc, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistCPSelectedV0PtTMaxEta=new TH3F(Form("fHistCPSelectedV0PtTMaxEta_%i", sysV0), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles with common parton (Casc, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistCPSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistCPSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistNOCPSelectedV0PtTMaxEta=new TH3F(Form("fHistNOCPSelectedV0PtTMaxEta_%i", sysV0), "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles with no common parton(Casc, primary, events w T>0)", 120, -30, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistNOCPSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistNOCPSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistSelectedV0PtPtTMax=new TH3F(Form("fHistSelectedV0PtPtTMax_%i",sysV0), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles (Casc, primary, events w T>0)", 300, 0, 30, 120, -30,30,  100, 0, 100 );
  fHistSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F *    fHistCPSelectedV0PtPtTMax=new TH3F(Form("fHistCPSelectedV0PtPtTMax_%i",sysV0), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles with common parton (Casc, primary, events w T>0)", 300, 0, 30, 120, -30,30,  100, 0, 100 );
  fHistCPSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistCPSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F *    fHistNOCPSelectedV0PtPtTMax=new TH3F(Form("fHistNOCPSelectedV0PtPtTMax_%i",sysV0), "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles with no common parton (Casc, primary, events w T>0)", 300, 0, 30, 120, -30,30,  100, 0, 100 );
  fHistNOCPSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistNOCPSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F*    fHistPrimaryV0[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPrimaryV0[j]=new TH3F(Form("fHistPrimaryV0_%i_cut%i",j, sysV0), "V0 MC (Casc, selected)", 4, 0.5, 4.5, 160, 0, 16, 60, 0, 30);
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s");
    fHistPrimaryV0[j]->GetYaxis()->SetTitle("p_{T}");
    fHistPrimaryV0[j]->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }

  //------------------Histograms os selected particles (trigger) for future efficiency calculation ----------------
  TH3F*      fHistSelectedTriggerPtPhi=new TH3F(Form("fHistSelectedTriggerPtPhi_%i",sysTrigger), "p_{T} and #phi distribution of selected trigger particles (primary)", 300, 0, 30, 400,0   , 2*TMath::Pi() ,  100, 0, 100);
  fHistSelectedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  TH3F*     fHistSelectedTriggerPtEta=new TH3F(Form("fHistSelectedTriggerPtEta_%i",sysTrigger), "p_{T} and #eta distribution of selected trigger particles (primary)", 300, 0, 30, 400,   1.2, 1.2,  100, 0, 100);
  fHistSelectedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtEta->GetYaxis()->SetTitle("#eta");


  TH2F*    fHistPrimaryTrigger[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPrimaryTrigger[j]=new TH2F(Form("fHistPrimaryTrigger_%i_cut%i", j, sysTrigger), "Trigger MC (selected)", 4, 0.5, 4.5, 100, 0,30 );
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(1,"Primary selected triggers");
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected triggers");
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected triggers");
    fHistPrimaryTrigger[j]->GetYaxis()->SetTitle("p_{T}");
  }
 
  /*-----------------------Pt Trigger  --------------------------- */
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
  /*-----------------------Pt Assoc  --------------------------- */
  TH1D *hSign_PtAssoc       = new TH1D("hSign_PtAssoc",       "hSign_PtAssoc",       300, 0, 30);
  hSign_PtAssoc->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssoc->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssoc->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssoc->GetYaxis()->SetLabelSize(0.05);

  TH1D *hSign_PtAssocTrue       = new TH1D("hSign_PtAssocTrue",       "hSign_PtAssocTrue",       300, 0, 30);
  hSign_PtAssocTrue->GetXaxis()->SetTitle("p_{T} (Gev/c)");
  hSign_PtAssocTrue->GetXaxis()->SetTitleSize(0.05);
  hSign_PtAssocTrue->GetXaxis()->SetLabelSize(0.05);
  hSign_PtAssocTrue->GetYaxis()->SetLabelSize(0.05);


  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hMassvsPt_SEbins[nummolt+1][numzeta][numPtTrigger];
  TH2D *hMassvsPt_SEbins_true[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCP[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCPTrue[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[nummolt+1][numzeta][numPtTrigger];

  //Form("hMassvsPt_"+tipo[type]+"_%i",molt
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	hMassvsPt_SEbins[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i", m), 100, LimInfMass[type], LimSupMass[type], 100, 0, 10);
	hMassvsPt_SEbins_true[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m), 100, LimInfMass[type], LimSupMass[type], 100, 0, 10);
	hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_CPTrue_PtInt", "SE_m"+Smolt[m]+"_CPTrue_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi()); 
	hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_NOCPTrue_PtInt", "SE_m"+Smolt[m]+"_NOCPTrue_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi()); 
	hDeltaEtaDeltaPhi_SEbins_CPPtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_CP_PtInt", "SE_m"+Smolt[m]+"_CP_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi()); 
	hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]= new TH2D("SE_m"+Smolt[m]+"_NOCP_PtInt", "SE_m"+Smolt[m]+"_NOCP_PtInt", 56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi()); 

	for(Int_t v=1; v<numPtV0; v++){
	  nameSE[m][z][v][tr]="SE_";
	  namemassSE[m][z][v][tr]="InvMassSE_";
	  //nameSE[m][z][v][tr]+=Form("m%i_z%i_v%i_tr%i",m,z,v,tr);
	  nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_SEbins[m][z][v][tr]= new TH1D(namemassSE[m][z][v][tr], namemassSE[m][z][v][tr],100, LimInfMass[type], LimSupMass[type]);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr], nameSE[m][z][v][tr],   56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi()); 
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

	  hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]= (TH2D*)	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_CPTrue");
	  hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]->SetTitle("AC for Xi (true) cand. with a common parton");
	  hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]= (TH2D*)	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_CP");
	  hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]->SetTitle("AC for Xi cand. with a common parton");
	  hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]= (TH2D*)	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_NOCPTrue");
	  hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]->SetTitle("AC for Xi (true) cand. without a common parton");
	  hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]= (TH2D*)	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Clone(nameSE[m][z][v][tr]+"_NOCP");
	  hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->SetTitle("AC for Xi cand. without a common parton");

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
	hMassvsPt_MEbins[m][z][tr]= new TH2D(Form("ME_hMassvsPt_"+tipo[type]+"_%i", m),Form("ME_hMassvsPt_"+tipo[type]+"_%i", m),100, LimInfMass[type], LimSupMass[type],100, 0, 10);
	hMassvsPt_MEbins_true[m][z][tr]= new TH2D(Form("ME_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("ME_hMassvsPt_"+tipo[type]+"_%i_true", m),100, LimInfMass[type], LimSupMass[type],100, 0, 10);
	for(Int_t v=1; v<numPtV0; v++){
	  nameME[m][z][v][tr]="ME_";
	  namemassME[m][z][v][tr]="InvMassME_";
	  nameME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassME[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_MEbins[m][z][v][tr]= new TH1D(namemassME[m][z][v][tr], namemassME[m][z][v][tr],   100, LimInfMass[type], LimSupMass[type]);
	  hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]= new TH2D(nameME[m][z][v][tr], nameME[m][z][v][tr],  56, -1.5, 1.5, 104,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
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
  EntriesBkg  =  tBkg ->GetEntries();

  Int_t InJet[nummolt+1]={0};
  Int_t OutJet[nummolt+1]={0};

    
  Bool_t MoltSel=kFALSE; 
  Float_t     fSignTreeVariableInvMassCasc= 0;
  Bool_t isCascTrue=kFALSE;

  // cout << "  entries Sign: " << EntriesSign<<endl;
  Int_t l=0;
  for(Int_t k = 0; k<EntriesSign; k++){
    //    if (k>100000) continue;
    tSign->GetEntry(k);  
    if (k==10000*l){
      l++;     
      cout << "SE ----" << (Float_t)k/EntriesSign<< endl;
    }
    fSignTreeVariableDeltaEta=fSignTreeVariableEtaV0-fSignTreeVariableEtaTrigger;

    //charge selection
    if ((type==0 || type ==2) && fSignTreeVariableChargeAssoc==1) continue;
    else if ((type==1 || type ==3) && fSignTreeVariableChargeAssoc==-1) continue;

    //inv mass definition
    fSignTreeVariableInvMassCasc= 0;
    if (type==0 || type==1 || type==4)      fSignTreeVariableInvMassCasc= fSignTreeVariableInvMassXi;
    else if (type==2 || type==3 || type==5)      fSignTreeVariableInvMassCasc= fSignTreeVariableInvMassOmega;

    if (CommonParton){
      MassInPeakCasc=0;
      if ( fSignTreeVariableInvMassCasc> 1.315 &&  fSignTreeVariableInvMassCasc<1.33)  MassInPeakCasc=1;
    }
    //rapidity selection
    if (israp==0 && TMath::Abs(fSignTreeVariableEtaV0)>0.8)continue; 
    else if (israp==1 && TMath::Abs(fSignTreeVariableRapAssoc)>0.5)continue;

    //definition of true cascade
    if (type<=3) isCascTrue= (fSignTreeVariablePDGCodeAssoc==PDGCode[type] );
    else if (type==4) isCascTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[0]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[1])  );
    else if (type==5) isCascTrue= ((fSignTreeVariablePDGCodeAssoc==PDGCode[2]) ||(fSignTreeVariablePDGCodeAssoc==PDGCode[3])  );

    if (SkipAssoc){    if (fSignTreeVariableSkipAssoc==1) continue;}

    if(isMC==0 || (isMC==1 && isEfficiency==1)){

      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fSignTreeVariablePtTrigger)<PtTrigMin) continue;
      if(TMath::Abs(fSignTreeVariablePtTrigger)>ptjmax) continue;

      //cuts on DCAz trigger*******************
      if(TMath::Abs(fSignTreeVariableDCAz)>DCAzTrigger) continue;
      CounterSignPairsIntermediate++;  

      //******************* some other cuts for sys studies**************************
      if (fSignTreeVariableDcaXiDaughters> DCAXiDaughters) continue;
      if (fSignTreeVariableV0CosineOfPointingAngle<CosinePAngleV0ToXi)continue;
      if (fSignTreeVariableXiCosineOfPointingAngle < CosinePAngleXiToPV)continue;
      if (TMath::Abs(fSignTreeVariableInvMassLambda -massLambda) > InvMassLambdaWindow)continue;
      if (fSignTreeVariablectau  > ctau)continue;
	
    }
    if (isCascTrue)      TrueCounterSignPairsAfterPtMinCut++;  
    CounterSignPairsAfterPtMinCut++;  

    tSignO->Fill();
    //**********************************************************************************************

    fSignTreeVariableDeltaPhi= -fSignTreeVariableDeltaPhi;
    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();


    if (isMC && isEfficiency) {
      if(fSignTreeVariableisPrimaryTrigger==1){
	fHistSelectedTriggerPtPhi->Fill(fSignTreeVariablePtTrigger,fSignTreeVariablePhiTrigger, fSignTreeVariableMultiplicity);
	fHistSelectedTriggerPtEta->Fill(fSignTreeVariablePtTrigger,fSignTreeVariableEtaTrigger, fSignTreeVariableMultiplicity);
      }

      for (Int_t m =0; m<5;m++){
	if(fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]){
	  for(Int_t p=1; p<=4; p++){
	    if (fSignTreeVariableisPrimaryTrigger==p){
	      fHistPrimaryTrigger[m]->Fill(p,fSignTreeVariablePtTrigger );
	    }
	  }
	}
      }
      for(Int_t p=1; p<=4; p++){
	if (fSignTreeVariableisPrimaryTrigger==p){
	  fHistPrimaryTrigger[5]->Fill(p,fSignTreeVariablePtTrigger );
	}
      }

    
      if(isCascTrue &&  fSignTreeVariableisPrimaryV0==1){
	if (CommonParton){
	  if ( fSignTreeVariableIsCommonParton)      {
	    fHistCPSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
	    fHistCPSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc , fSignTreeVariableMultiplicity);
	    fHistCPSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);

	  }
	  else     {
	    fHistNOCPSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
	    fHistNOCPSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc , fSignTreeVariableMultiplicity);
	    fHistNOCPSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
	  }
	}
	fHistSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
	fHistSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
	fHistSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger*fSignTreeVariableChargeAssoc , fSignTreeVariableMultiplicity);
      }
    

      for (Int_t m =0; m<5;m++){
	if(fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]){
	  for(Int_t p=1; p<=4; p++){
	    if (fSignTreeVariableisPrimaryV0==p){
	      fHistPrimaryV0[m]->Fill(p,fSignTreeVariablePtV0, fSignTreeVariablePtTrigger );
	    }
	  }
	}
      }
      for(Int_t p=1; p<=4; p++){
	if (fSignTreeVariableisPrimaryV0==p){
	  fHistPrimaryV0[5]->Fill(p,fSignTreeVariablePtV0, fSignTreeVariablePtTrigger);
	}
      }
    }

    hMultiplicityDIstributionACEv->Fill(fSignTreeVariableMultiplicity);
    hSign_PtAssoc->Fill(fSignTreeVariablePtV0);
    hSign_PtTrigger->Fill(fSignTreeVariablePtTrigger);
    if (isCascTrue) {
      hSign_PtAssocTrue->Fill(fSignTreeVariablePtV0);
    }
    if ((isMC && isCascTrue) || !isMC){
      if (CommonParton){
	if (fSignTreeVariableIsCommonParton)      hSign_PtAssoc_CPTrue->Fill(fSignTreeVariablePtV0);
	else      hSign_PtAssoc_NOCPTrue->Fill(fSignTreeVariablePtV0);
      }
    }


    for(Int_t m=0; m<nummolt+1; m++){
      if (m < nummolt){
	MoltSel = (fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1]);
      }
      else if(m==nummolt) MoltSel=kTRUE;
      // cout << m << "  " << nummolt << "  " << MoltSel<< endl;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  if(MoltSel){
	    CounterSignPairsAfterPtMinCutMult[m]++;  
	    if (isCascTrue)      TrueCounterSignPairsAfterPtMinCutMult[m]++;  
	    hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMassCasc, fSignTreeVariablePtV0); 
	    if(isCascTrue) hMassvsPt_SEbins_true[m][z][tr]->Fill(fSignTreeVariableInvMassCasc, fSignTreeVariablePtV0); 
	    if (CommonParton){
	      if ((isMC && isCascTrue)){
		if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CPTruePtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		else   hDeltaEtaDeltaPhi_SEbins_NOCPTruePtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }
	      if ( MassInPeakCasc){
	      if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CPPtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      else   hDeltaEtaDeltaPhi_SEbins_NOCPPtInt[m][z][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }
	      if (fSignTreeVariableIsCommonParton ){
		if (TMath::Abs(fSignTreeVariableDeltaPhi)<1) InJet[m]++;
		else OutJet[m]++;
	      }
	    }
	    
	  }

	  for(Int_t v=1; v<numPtV0; v++){
	    if(MoltSel && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      hInvMassK0Short_SEbins[m][z][v][tr]->Fill(fSignTreeVariableInvMassCasc);
	      if (CommonParton){
		if ((isMC && isCascTrue) || !isMC){
		  if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CPTrue[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		  else   hDeltaEtaDeltaPhi_SEbins_NOCPTrue[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		}
		if ( MassInPeakCasc){
		if (fSignTreeVariableIsCommonParton) hDeltaEtaDeltaPhi_SEbins_CP[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		else   hDeltaEtaDeltaPhi_SEbins_NOCP[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
		}

	      }

	    }
	  }
	}
      }
    }
  }


  Float_t     fBkgTreeVariableInvMassCasc= 0;
  l=0;
  for(Int_t k = 0; k<EntriesBkg; k++){
    // for(Int_t k = 0; k<1; k++){
    //    if (k>100000) continue;
    tBkg->GetEntry(k);
    if (k==10000*l){
      l++;     
      cout << "ME ----" << (Float_t)k/EntriesBkg<< endl;
    }
    fBkgTreeVariableDeltaEta=fBkgTreeVariableEtaV0-fBkgTreeVariableEtaTrigger;

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

    if (SkipAssoc){    if (fBkgTreeVariableSkipAssoc==1) continue;}

    if(isMC==0 || (isMC==1 && isEfficiency==1)){
      //************cuts on pT trigger min*********************************
      if(TMath::Abs(fBkgTreeVariablePtTrigger)<PtTrigMin) continue;
      if(TMath::Abs(fBkgTreeVariablePtTrigger)>ptjmax) continue;

      //cuts on DCAz trigger*******************
      if(TMath::Abs(fBkgTreeVariableDCAz)>DCAzTrigger) continue;

      CounterBkgPairsIntermediate++;  
      //******************* some other cuts for sys studies **************************

      if (fBkgTreeVariableDcaXiDaughters> DCAXiDaughters) continue;
      if (fBkgTreeVariableV0CosineOfPointingAngle<CosinePAngleV0ToXi)continue;
      if (fBkgTreeVariableXiCosineOfPointingAngle < CosinePAngleXiToPV)continue;
      if (TMath::Abs(fBkgTreeVariableInvMassLambda -massLambda) > InvMassLambdaWindow)continue;
      if (fBkgTreeVariablectau  > ctau)continue;

    }
    if (isCascTrue)      TrueCounterBkgPairsAfterPtMinCut++;  
    CounterBkgPairsAfterPtMinCut++;  

    tBkgO->Fill();
    //**********************************************************************************************

    fBkgTreeVariableDeltaPhi= -fBkgTreeVariableDeltaPhi;
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
	    CounterBkgPairsAfterPtMinCutMult[m]++;  
	    if (isCascTrue)      TrueCounterBkgPairsAfterPtMinCutMult[m]++;  
	    hMassvsPt_MEbins[m][z][tr]->Fill(fBkgTreeVariableInvMassCasc, fBkgTreeVariablePtV0);
	    if(isCascTrue) 	    hMassvsPt_MEbins_true[m][z][tr]->Fill(fBkgTreeVariableInvMassCasc, fBkgTreeVariablePtV0);
	  }
	  for(Int_t v=1; v<numPtV0; v++){
	    if(MoltSel && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0<NPtV0[v+1]){
	      hInvMassK0Short_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableInvMassCasc);
	      hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	    }
	  }
	}
      }
    }
  }

  cout << " ratio between Xi with CommonParton (CP) and Xi without in the different pt bins " << endl;
  Float_t      CPTrue[numPtV0] ={0};
  Float_t      NOCPTrue[numPtV0]={0};

  for (Int_t v=0; v< numPtV0; v++){
    for(Int_t b= hSign_PtAssoc_CPTrue->GetXaxis()->FindBin(NPtV0[v]+0.001); b<= hSign_PtAssoc_CPTrue->GetXaxis()->FindBin(NPtV0[v+1]-0.001); b++ ){
      CPTrue[v] +=  hSign_PtAssoc_CPTrue->GetBinContent(b);
      NOCPTrue[v] +=  hSign_PtAssoc_NOCPTrue->GetBinContent(b);
    }
    cout << "v " << v << " " <<  CPTrue[v]/NOCPTrue[v]<< endl;
  }

   

  
  cout << "Pt Min delle particelle trigger " << PtTrigMin<< endl;

  cout << "\nsignal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterSignPairsAfterPtMinCut <<" true: " <<   TrueCounterSignPairsAfterPtMinCut << " -> all entries in sign Tree were : " <<   EntriesSign << " " <<(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign << endl;
  cout << "bkg pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterBkgPairsAfterPtMinCut << " true: " <<   TrueCounterBkgPairsAfterPtMinCut << " -> all entries in sign Tree were : " <<   EntriesBkg << " " <<(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg <<  endl;

  for (Int_t m=0; m< nummolt+1; m++){
    cout << m << endl;  
    cout << "signal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterSignPairsAfterPtMinCutMult[m] << ", " << (Float_t)CounterSignPairsAfterPtMinCutMult[m]/CounterSignPairsAfterPtMinCut<<endl;
    cout << "bkgal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterBkgPairsAfterPtMinCutMult[m] << ", " << (Float_t)CounterBkgPairsAfterPtMinCutMult[m]/CounterBkgPairsAfterPtMinCut<<endl;
  }

  cout << "\nsignal pairs trigger-associated (after all selections, included Pt min cut) " << 	CounterSignPairsAfterPtMinCut <<" percentage of INT7 events with at least one selected V0 (calculated assuming 1V0 per event in which there is a V0) "<<(Float_t)CounterSignPairsAfterPtMinCut/TotEvtINT7  <<  endl;

  cout << "entries of V0 selected histogram (all Casc primary true ) "<<  fHistSelectedV0PtTMaxPhi->GetEntries()<< endl;
  cout << "entries of V0 primary histogram (all Casc true ) "<<  fHistPrimaryV0[5]->GetEntries()<< endl;

  cout << "\n Other useful information; " << endl;
  cout << "average pT cascade " <<   hSign_PtAssoc->GetMean()<< " entries: " <<  endl;
  cout << "average pT trigger particles in AC events selected " <<     hSign_PtTrigger->GetMean()<< endl;
  //  cout << "average pT trigger particles in events with NT>0  " << << endl;
  cout << "\n\npartendo dal file " << PathIn << " ho creato il file " << PathOut<< endl;

  cout <<"\n\nINT7     " << " Ev. NT>0    " << "NV0/INT7  " << "<pT,Xi>  " << "<pT,Trig>   " << "SE pairs   " << "ME pairs   " << "Mult. distr                   " << "NV0/ev mult " << endl;  
  cout << std::setprecision(2);
  cout << "    " << TotEvtINT7;
  cout << "    " <<   (Float_t)hMultiplicityBefAll->GetEntries()/TotEvtINT7;
  cout << "    " <<(Float_t)CounterSignPairsAfterPtMinCut/TotEvtINT7;
  cout << std::setprecision(3);
  cout << "    " << hSign_PtAssoc->GetMean();
  cout << "    " <<hSign_PtTrigger->GetMean();
  cout << std::setprecision(2);
  cout << "      " <<(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign;
  cout << "      " <<(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg;
  cout << "      " ;
  for (Int_t m=0; m<nummolt; m++){
    if (m<nummolt-1)    cout << ACcounter[m]/ACcounter[5]<< "-";
    else     cout << ACcounter[m]/ACcounter[5]<< "      ";
  }
  cout << std::setprecision(3);
  for (Int_t m=0; m<nummolt; m++){
    if (m<nummolt-1)    cout <<  hMultvsNumberAssoc_Proj[m]->GetMean() <<"-";
    else    cout <<  hMultvsNumberAssoc_Proj[m]->GetMean() <<endl;
  }

  HistoInfo->SetBinContent(1,(Float_t)hMultiplicityBefAll->GetEntries()/TotEvtINT7);
  HistoInfo->SetBinContent(2,(Float_t)CounterSignPairsAfterPtMinCut/TotEvtINT7);
  HistoInfo->SetBinContent(3,hSign_PtAssoc->GetMean());
  HistoInfo->SetBinContent(4,hSign_PtTrigger->GetMean());
  HistoInfo->SetBinContent(5,(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign);
  HistoInfo->SetBinContent(6,(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg);
  for (Int_t m=0; m< nummolt; m++){
    HistoInfo->SetBinContent(7+m,(Float_t)CounterSignPairsAfterPtMinCutMult[m]/CounterSignPairsAfterPtMinCut);
    HistoInfo->SetBinContent(12+m, hMultvsNumberAssoc_Proj[m]->GetMean());
  }

  for (Int_t m=0; m< nummolt+1; m++){
    cout << "Ratio between |#Delta#phi|<1 and all the rest for Xi cand. with origin in common with trigger  " << (Float_t)InJet[m]/OutJet[m] << endl;
  }


  cout <<"check if you have to change histogram binning for selected trigger particles (double the number of bins) " << endl;

  TCanvas *DeltaPhiProj;
  TCanvas *DeltaPhiProjMult[nummolt+1];


  TCanvas *canvas = new TCanvas ("canvas", "canvas",1300,800);
  TCanvas *canvasAC = new TCanvas ("canvasAC", "canvasAC",1300,800);
  canvas->Divide(4,2);
  canvasAC->Divide(4,2);
  TH1F*  hDeltaEtaDeltaPhi_SEbinsProj[nummolt+1][numzeta][numPtV0][numPtTrigger];

    for(Int_t m=0; m<nummolt+1; m++){
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for (Int_t v=1; v< numPtV0; v++){
	    //	    cout << " v " << v << endl;
	    canvasAC->cd(v+1);
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Draw("colz ");
	    canvas->cd(v+1);
	    hDeltaEtaDeltaPhi_SEbinsProj[m][z][v][tr] = (TH1F*)  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->ProjectionY( "SE_m"+Smolt[m]+"_v"+SPtV0[v]+"_Proj",0, -1, "E");

	    hDeltaEtaDeltaPhi_SEbinsProj[m][z][v][tr]->Rebin(2);
	    hDeltaEtaDeltaPhi_SEbinsProj[m][z][v][tr]->GetYaxis()->SetRangeUser(0,1.2*	    hDeltaEtaDeltaPhi_SEbinsProj[m][z][v][tr]->GetMaximum());
	    hDeltaEtaDeltaPhi_SEbinsProj[m][z][v][tr]->Draw("");
	  }
	}
      }
    }

    fout->WriteTObject(canvas);
    fout->WriteTObject(canvasAC);

    fout->Write();

    cout << "cut values " << endl;
    cout << "DcaXiDaughters:" << DCAXiDaughters  <<endl;
    cout << "CosinePAngleXiToPV   :"  <<CosinePAngleXiToPV     <<endl;
    cout << "CosinePAngleV0ToXi   :"  <<CosinePAngleV0ToXi    <<endl;
    cout << "InvMassLambdaWindow  :" << InvMassLambdaWindow << endl;
    cout << "ctau          :"  <<ctau            <<endl;
    cout << "DCAzTrigger   :"  <<DCAzTrigger            <<endl;

    cout << std::setprecision(8);
    cout << "signal pairs after selection/All signal pairs in initial tree      " <<(Float_t)CounterSignPairsAfterPtMinCut/EntriesSign<<endl;
    cout << "signal pairs after selection/All signal pairs in initial tree      " <<(Float_t)CounterSignPairsAfterPtMinCut << " over " << EntriesSign<<endl;
    cout << "intermediate counter for sign pairs " <<       CounterSignPairsIntermediate<< endl;
    cout << "bkg pairs after selection/All bkg pairs in initial tree      " <<(Float_t)CounterBkgPairsAfterPtMinCut/EntriesBkg<<endl;
    cout << "bkg pairs after selection/All bkg pairs in initial tree      " <<(Float_t)CounterBkgPairsAfterPtMinCut << " over " << EntriesBkg<<endl;
    cout << "intermediate counter for bkg pairs " <<       CounterBkgPairsIntermediate<< endl;

}

