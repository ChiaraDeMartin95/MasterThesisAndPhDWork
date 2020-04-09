
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
void readTreePLChiara_first(Bool_t ishhCorr=0, Float_t PtTrigMin=3.0, Float_t ptjmax=15, Int_t sysV0=0, bool isMC = 1,Bool_t isEfficiency=1,Int_t sysTrigger=0,
			    //TString PathIn ="./AnalysisResults2018d8_allrunswMult_MC.root", Int_t type=0 //type = 0 per K0short
			    TString year="2018f1_extra_onlyTriggerWithHighestPt", TString year0="2016", TString Path1 ="", Int_t type=0 //type = 0 per K0short 
			    
			    )
{

  if (ishhCorr) {
    cout << "This macro should not be run is hh correlation is studied; go directly to readTreePLChiara_second " << endl;
    return;
  }

  //2016k_onlyTriggerWithHighestPt
  if (sysV0>6) return;
  if (sysV0>2 && ishhCorr) return;
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
 
  PathIn+=Path1;
  PathIn+=".root";
  PathOut+=Path1;
  if (ishhCorr) PathOut +="_hhCorr";
  PathOut +=Form("_MassDistr_SysT%i_SysV0%i_PtMin%.1f",sysTrigger, sysV0, PtTrigMin); 
  PathOut+= ".root";

  TFile *fin = new TFile(PathIn);
  TDirectoryFile *d = (TDirectoryFile*)fin->Get("MyTask");


  TTree *tSign = (TTree *)d->Get("fSignalTree");
  TTree *tBkg  = (TTree *)d->Get("fBkgTree");
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Double_t ctauK0s = 2.6844;

  Double_t     fSignTreeVariablePtTrigger;		       
  Double_t     fSignTreeVariableChargeTrigger;		       
  Double_t     fSignTreeVariableEtaTrigger;		       
  Double_t     fSignTreeVariablePhiTrigger;		       
  Double_t     fSignTreeVariableDCAz;			       
  Double_t     fSignTreeVariableDCAxy;			  
  Double_t     fSignTreeVariableChargeAssoc;		       
  Double_t     fSignTreeVariableAssocDCAz;			       
  Double_t     fSignTreeVariableAssocDCAxy;			  
  Double_t     fSignTreeVariableisPrimaryTrigger;			  
  Double_t     fSignTreeVariableisPrimaryV0;			  
  Double_t     fSignTreeVariableRapK0Short;		       
  Double_t     fSignTreeVariableDcaV0ToPrimVertex;	       
  Double_t     fSignTreeVariableDcaPosToPrimVertex;	       
  Double_t     fSignTreeVariableDcaNegToPrimVertex;	       
  Double_t     fSignTreeVariableV0CosineOfPointingAngle;	       
  Double_t     fSignTreeVariablectau;			       
  Double_t     fSignTreeVariablePtV0;			       
  Double_t     fSignTreeVariableInvMassK0s;
  Double_t     fSignTreeVariableInvMassLambda;
  Double_t     fSignTreeVariableInvMassAntiLambda;		      		       
  Double_t     fSignTreeVariableEtaV0;			       
  Double_t     fSignTreeVariablePhiV0;			       
  Double_t     fSignTreeVariablePtArmenteros;		       
  Double_t     fSignTreeVariableAlpha;			       
  Double_t     fSignTreeVariableDeltaEta;		       
  Double_t     fSignTreeVariableDeltaPhi;		       
  Double_t     fSignTreeVariableDeltaTheta;		       
  Double_t     fSignTreeVariableMultiplicity;		       
  Double_t     fSignTreeVariableZvertex;                        
  Double_t     fSignTreeVariablePDGCode;                        


  Double_t        fBkgTreeVariablePtTrigger;
  Double_t        fBkgTreeVariableChargeTrigger;
  Double_t        fBkgTreeVariableEtaTrigger;
  Double_t        fBkgTreeVariablePhiTrigger;
  Double_t        fBkgTreeVariableDCAz;
  Double_t        fBkgTreeVariableDCAxy;
  Double_t        fBkgTreeVariableChargeAssoc;		       
  Double_t        fBkgTreeVariableAssocDCAz;			       
  Double_t        fBkgTreeVariableAssocDCAxy;			  
  Double_t        fBkgTreeVariableisPrimaryTrigger;			  
  Double_t        fBkgTreeVariableisPrimaryV0;			  
  Double_t        fBkgTreeVariableRapK0Short;
  Double_t        fBkgTreeVariableDcaV0ToPrimVertex;
  Double_t        fBkgTreeVariableDcaPosToPrimVertex;
  Double_t        fBkgTreeVariableDcaNegToPrimVertex;
  Double_t        fBkgTreeVariableV0CosineOfPointingAngle;
  Double_t        fBkgTreeVariablectau;
  Double_t        fBkgTreeVariablePtV0;
  Double_t        fBkgTreeVariableInvMassK0s;
  Double_t        fBkgTreeVariableInvMassLambda;
  Double_t        fBkgTreeVariableInvMassAntiLambda;		      		       
  Double_t        fBkgTreeVariableEtaV0;
  Double_t        fBkgTreeVariablePhiV0;
  Double_t        fBkgTreeVariablePtArmenteros;
  Double_t        fBkgTreeVariableAlpha;
  Double_t        fBkgTreeVariableDeltaEta;
  Double_t        fBkgTreeVariableDeltaPhi;
  Double_t        fBkgTreeVariableDeltaTheta;
  Double_t        fBkgTreeVariableMultiplicity;
  Double_t        fBkgTreeVariableZvertex;
  Double_t        fBkgTreeVariablePDGCode;

  Int_t CounterSignPairsAfterPtMinCut=0; 
  Int_t CounterBkgPairsAfterPtMinCut=0; 

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=7;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  TString tipo[numtipo]={"kK0s", "bo"};


  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "all"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  //TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-1.5", "1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,1.5,2,2.5,3,4,8};
  //   TString SPtV0[numPtV0]={"0-1", "1-2", "2,3", "3-4", "4-8"};
  //   Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  Double_t NPtTrigger[numPtTrigger+1]={PtTrigMin,ptjmax};

  
  //what is the fraction of AC events in each multiplicity class?
  TList *d1 = (TList*)d->Get("MyOutputContainer");
  if (!d1) return;
  TH1F* hMultiplicity=(TH1F*)  d1->FindObject("fHist_multiplicity");
  TH2F* hMultvsNumberAssoc=(TH2F*)  d1->FindObject("fHistMultvsV0");
  TH1F*  hMultvsNumberAssoc_Proj[nummolt];
  if (!hMultiplicity) cout << " no info about multiplicity distribution of AC events available " << endl;
  if (!hMultvsNumberAssoc) cout << " no info about multiplicity distribution of AC events available " << endl;
  Float_t ACcounter[nummolt+1];
  if (hMultiplicity && hMultvsNumberAssoc){
    ACcounter[5]= hMultiplicity->GetEntries();
    cout <<   hMultiplicity->GetEntries() << endl;
    for (Int_t m=0; m< nummolt; m++){ 
      hMultvsNumberAssoc_Proj[m] = (TH1F*)       hMultvsNumberAssoc->ProjectionX(Form("hMultvsNumberAssoc_%i", m), hMultvsNumberAssoc->GetYaxis()->FindBin(Nmolt[m]+0.001), hMultvsNumberAssoc->GetYaxis()->FindBin(Nmolt[m+1]-0.001));

    cout << " m " << m << endl;
    ACcounter[m] =0;
    for (Int_t b= hMultiplicity->GetXaxis()->FindBin(Nmolt[m]+0.001); b <=  hMultiplicity->GetXaxis()->FindBin(Nmolt[m+1]-0.001); b++){
      ACcounter[m] +=  hMultiplicity->GetBinContent(b);
    }
    cout << "fraction of events in mult bin " << Smolt[m] << ": " << ACcounter[m]/ACcounter[5] <<  " ~average V0 number " <<     hMultvsNumberAssoc_Proj[m]->GetMean() <<endl;
    }
    
  }

  
  //Signal
  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);		       
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);		       
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);		       
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);		       
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);			       
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);			       
  // tSign->SetBranchAddress("fTreeVariableChargeAssoc"             ,&fSignTreeVariableChargeAssoc);		       	   
  // tSign->SetBranchAddress("fTreeVariableAssocDCAz"                      ,&fSignTreeVariableAssocDCAz);			       
  // tSign->SetBranchAddress("fTreeVariableAssocDCAxy"                     ,&fSignTreeVariableAssocDCAxy);	

  tSign->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fSignTreeVariableisPrimaryTrigger);			       
  tSign->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fSignTreeVariableisPrimaryV0);			       
  tSign->SetBranchAddress("fTreeVariableRapK0Short"                ,&fSignTreeVariableRapK0Short);		       
  tSign->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex"         ,&fSignTreeVariableDcaV0ToPrimVertex);	       
  tSign->SetBranchAddress("fTreeVariableDcaPosToPrimVertex"        ,&fSignTreeVariableDcaPosToPrimVertex);	       
  tSign->SetBranchAddress("fTreeVariableDcaNegToPrimVertex"        ,&fSignTreeVariableDcaNegToPrimVertex);	       
  tSign->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fSignTreeVariableV0CosineOfPointingAngle);	       
  tSign->SetBranchAddress("fTreeVariablectau"                      ,&fSignTreeVariablectau);			       
  tSign->SetBranchAddress("fTreeVariablePtV0"                      ,&fSignTreeVariablePtV0);			       
  tSign->SetBranchAddress("fTreeVariableInvMassK0s"                ,&fSignTreeVariableInvMassK0s);		       
  tSign->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fSignTreeVariableInvMassLambda);		       
  tSign->SetBranchAddress("fTreeVariableInvMassAntiLambda"         ,&fSignTreeVariableInvMassAntiLambda);		       
  tSign->SetBranchAddress("fTreeVariableEtaV0"                     ,&fSignTreeVariableEtaV0);			       
  tSign->SetBranchAddress("fTreeVariablePhiV0"                     ,&fSignTreeVariablePhiV0);			       
  tSign->SetBranchAddress("fTreeVariablePtArmenteros"              ,&fSignTreeVariablePtArmenteros);		       
  tSign->SetBranchAddress("fTreeVariableAlpha"                     ,&fSignTreeVariableAlpha);			       
  tSign->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fSignTreeVariableDeltaEta);		       
  tSign->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fSignTreeVariableDeltaPhi);		       
  tSign->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fSignTreeVariableDeltaTheta);		       
  tSign->SetBranchAddress("fTreeVariableMultiplicity"              ,&fSignTreeVariableMultiplicity);		       
  tSign->SetBranchAddress("fTreeVariableZvertex"                   ,&fSignTreeVariableZvertex);                        
  tSign->SetBranchAddress("fTreeVariablePDGCode"                   ,&fSignTreeVariablePDGCode);                         
  /*  if(isMC){
      tSign->SetBranchAddress("tMCtruepair"       ,&tSignMCtruepair  );
      tSign->SetBranchAddress("tMCSameMother"     ,&tSignMCSameMother);
      tSign->SetBranchAddress("tMCMotherP1"       ,&tSignMCMotherP1  );
      tSign->SetBranchAddress("tMCMotherP2"       ,&tSignMCMotherP2  );
      tSign->SetBranchAddress("tMCptcTypeP1"      ,&tSignMCptcTypeP1 ); 
      tSign->SetBranchAddress("tMCptcTypeP2"      ,&tSignMCptcTypeP2 ); 
      tSign->SetBranchAddress("tMCSameGM"         ,&tSignMCSameGM    );
      tSign->SetBranchAddress("tMotherPDG"        ,&tSignMotherPDG   );
      tSign->SetBranchAddress("tMotherPDGP2"      ,&tSignMotherPDGP2 );
      tSign->SetBranchAddress("tpdgcodeP1"        ,&tSignpdgcodeP1   ); 
      tSign->SetBranchAddress("tpdgcodeP2"        ,&tSignpdgcodeP2   );
      tSign->SetBranchAddress("tKstarGen"         ,&tSignKstarGen    );
      }
  */
  //BackGround
  tBkg->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fBkgTreeVariablePtTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fBkgTreeVariableChargeTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fBkgTreeVariableEtaTrigger);		       
  tBkg->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fBkgTreeVariablePhiTrigger);		       
  tBkg->SetBranchAddress("fTreeVariableDCAz"                      ,&fBkgTreeVariableDCAz);			       
  tBkg->SetBranchAddress("fTreeVariableDCAxy"                     ,&fBkgTreeVariableDCAxy);	
  // tBkg->SetBranchAddress("fTreeVariableChargeAssoc"             ,&fBkgTreeVariableChargeAssoc);		       
  // tBkg->SetBranchAddress("fTreeVariableAssocDCAz"                      ,&fBkgTreeVariableAssocDCAz);			       
  // tBkg->SetBranchAddress("fTreeVariableAssocDCAxy"                     ,&fBkgTreeVariableAssocDCAxy);			    
  
  tBkg->SetBranchAddress("fTreeVariableisPrimaryTrigger"          ,&fBkgTreeVariableisPrimaryTrigger);			       
  tBkg->SetBranchAddress("fTreeVariableisPrimaryV0"               ,&fBkgTreeVariableisPrimaryV0);
  tBkg->SetBranchAddress("fTreeVariableRapK0Short"                ,&fBkgTreeVariableRapK0Short);		       
  tBkg->SetBranchAddress("fTreeVariableDcaV0ToPrimVertex"         ,&fBkgTreeVariableDcaV0ToPrimVertex);	       
  tBkg->SetBranchAddress("fTreeVariableDcaPosToPrimVertex"        ,&fBkgTreeVariableDcaPosToPrimVertex);	       
  tBkg->SetBranchAddress("fTreeVariableDcaNegToPrimVertex"        ,&fBkgTreeVariableDcaNegToPrimVertex);	       
  tBkg->SetBranchAddress("fTreeVariableV0CosineOfPointingAngle"   ,&fBkgTreeVariableV0CosineOfPointingAngle);	        
  tBkg->SetBranchAddress("fTreeVariablectau"                      ,&fBkgTreeVariablectau);			       
  tBkg->SetBranchAddress("fTreeVariablePtV0"                      ,&fBkgTreeVariablePtV0);			       
  tBkg->SetBranchAddress("fTreeVariableInvMassK0s"                ,&fBkgTreeVariableInvMassK0s);		       
  tBkg->SetBranchAddress("fTreeVariableInvMassLambda"             ,&fBkgTreeVariableInvMassLambda);		       
  tBkg->SetBranchAddress("fTreeVariableInvMassAntiLambda"         ,&fBkgTreeVariableInvMassAntiLambda);		       
  tBkg->SetBranchAddress("fTreeVariableEtaV0"                     ,&fBkgTreeVariableEtaV0);			       
  tBkg->SetBranchAddress("fTreeVariablePhiV0"                     ,&fBkgTreeVariablePhiV0);			       
  tBkg->SetBranchAddress("fTreeVariablePtArmenteros"              ,&fBkgTreeVariablePtArmenteros);		       
  tBkg->SetBranchAddress("fTreeVariableAlpha"                     ,&fBkgTreeVariableAlpha);			       
  tBkg->SetBranchAddress("fTreeVariableDeltaEta"                  ,&fBkgTreeVariableDeltaEta);		       
  tBkg->SetBranchAddress("fTreeVariableDeltaPhi"                  ,&fBkgTreeVariableDeltaPhi);		       
  tBkg->SetBranchAddress("fTreeVariableDeltaTheta"                ,&fBkgTreeVariableDeltaTheta);		       
  tBkg->SetBranchAddress("fTreeVariableMultiplicity"              ,&fBkgTreeVariableMultiplicity);		       
  tBkg->SetBranchAddress("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex);                          tBkg->SetBranchAddress("fTreeVariablePDGCode"                   ,&fBkgTreeVariablePDGCode);                        

  /* if(isMC){
     tBkg->SetBranchAddress("tMCtruepair"       ,&tBkgMCtruepair   );
     tBkg->SetBranchAddress("tMCSameMother"     ,&tBkgMCSameMother );
     tBkg->SetBranchAddress("tMCMotherP1"       ,&tBkgMCMotherP1   );
     tBkg->SetBranchAddress("tMCMotherP2"       ,&tBkgMCMotherP2   );
     tBkg->SetBranchAddress("tMCptcTypeP1"      ,&tBkgMCptcTypeP1  ); 
     tBkg->SetBranchAddress("tMCptcTypeP2"      ,&tBkgMCptcTypeP2  ); 
     tBkg->SetBranchAddress("tMCSameGM"         ,&tBkgMCSameGM     );
     tBkg->SetBranchAddress("tMotherPDG"        ,&tBkgMotherPDG    );
     tBkg->SetBranchAddress("tMotherPDGP2"      ,&tBkgMotherPDGP2  );
     tBkg->SetBranchAddress("tpdgcodeP1"        ,&tBkgpdgcodeP1    ); 
     tBkg->SetBranchAddress("tpdgcodeP2"        ,&tBkgpdgcodeP2    );
     tBkg->SetBranchAddress("tKstarGen"         ,&tBkgKstarGen     );
     }
  */

  Int_t EntriesSign = 0; 
  Int_t EntriesBkg  = 0; 
  
  TFile *fout = new TFile(PathOut,"RECREATE");
  TDirectory  *dirSign= fout->mkdir("SE");
  TDirectory  *dirBkg= fout->mkdir("ME");


  //------------------Histograms os selected particles (V0) for future efficiency calculation ----------------
  TH3F*    fHistSelectedV0PtTMaxPhi=new TH3F("fHistSelectedV0PtTMaxPhi", "p^{Trigg, Max}_{T} and #phi distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
  fHistSelectedV0PtTMaxPhi->GetXaxis()->SetTitle("p^{Trigg, Max}_{T}");
  fHistSelectedV0PtTMaxPhi->GetYaxis()->SetTitle("#phi");

  TH3F *    fHistSelectedV0PtTMaxEta=new TH3F("fHistSelectedV0PtTMaxEta", "p^{Trigg, Max}_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistSelectedV0PtTMaxEta->GetXaxis()->SetTitle("p^{Trigg, max}_{T}");
  fHistSelectedV0PtTMaxEta->GetYaxis()->SetTitle("#eta");

  TH3F *    fHistSelectedV0PtPtTMax=new TH3F("fHistSelectedV0PtPtTMax", "p_{T} and p^{Trigg, Max}_{T} distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 60, 0,30,  100, 0, 100 );
  fHistSelectedV0PtPtTMax->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtPtTMax->GetYaxis()->SetTitle("p^{Trigg, Max}_{T}");

  TH3F*    fHistPrimaryV0[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPrimaryV0[j]=new TH3F(Form("fHistPrimaryV0_%i",j), "V0 MC (K0s, selected)", 4, 0.5, 4.5, 160, 0, 16, 60, 0, 30);
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s");
    fHistPrimaryV0[j]->GetYaxis()->SetTitle("p_{T}");
    fHistPrimaryV0[j]->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");
  }

  //------------------Histograms os selected particles (trigger) for future efficiency calculation ----------------
  TH3F*      fHistSelectedTriggerPtPhi=new TH3F("fHistSelectedTriggerPtPhi", "p_{T} and #phi distribution of selected trigger particles (primary)", 300, 0, 30, 400,0   , 2*TMath::Pi() ,  100, 0, 100);
  fHistSelectedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  TH3F*     fHistSelectedTriggerPtEta=new TH3F("fHistSelectedTriggerPtEta_%i", "p_{T} and #eta distribution of selected trigger particles (primary)", 300, 0, 30, 400,   1.2, 1.2,  100, 0, 100);
  fHistSelectedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtEta->GetYaxis()->SetTitle("#eta");


  TH2F*    fHistPrimaryTrigger[nummolt+1];
  for(Int_t j=0; j<nummolt+1; j++){
    fHistPrimaryTrigger[j]=new TH2F(Form("fHistPrimaryTrigger_%i", j), "Trigger MC (selected)", 4, 0.5, 4.5, 100, 0,30 );
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


  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hInvMassK0Short_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hMassvsPt_SEbins[nummolt+1][numzeta][numPtTrigger];
  TH2D *hMassvsPt_SEbins_true[nummolt+1][numzeta][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];

  //Form("hMassvsPt_"+tipo[type]+"_%i",molt
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	hMassvsPt_SEbins[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i", m), 100, 0.4, 0.6, 100, 0, 10);
	hMassvsPt_SEbins_true[m][z][tr]= new TH2D(Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m),Form("SE_hMassvsPt_"+tipo[type]+"_%i_true", m), 100, 0.4, 0.6, 100, 0, 10);
	for(Int_t v=0; v<numPtV0; v++){
	  nameSE[m][z][v][tr]="SE_";
	  namemassSE[m][z][v][tr]="InvMassSE_";
	  //nameSE[m][z][v][tr]+=Form("m%i_z%i_v%i_tr%i",m,z,v,tr);
	  nameSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  namemassSE[m][z][v][tr]+="m"+ Smolt[m]+"_v"+SPtV0[v];
	  hInvMassK0Short_SEbins[m][z][v][tr]= new TH1D(namemassSE[m][z][v][tr], namemassSE[m][z][v][tr],100, 0.4, 0.6);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]= new TH2D(nameSE[m][z][v][tr], nameSE[m][z][v][tr],   50, -1.5, 1.5, 100,  -0.5*TMath::Pi(), 1.5*TMath::Pi());
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitle("#Delta #phi (rad)");
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetTitleOffset(1.5);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetLabelSize(0.05);
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetYaxis()->SetLabelSize(0.05);

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
	  hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #eta");
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
    
  Bool_t MoltSel=kFALSE; 

  dirSign->cd();
  // cout << "  entries Sign: " << EntriesSign<<endl;
  for(Int_t k = 0; k<EntriesSign; k++){
    tSign->GetEntry(k);  
    if(isMC==0 || (isMC==1 && isEfficiency==1)){

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

      //******************* some other cuts for sys studies**************************
      if (!ishhCorr){
	if(sysV0!=5){
	  if(TMath::Abs((fSignTreeVariableInvMassLambda - massLambda))<= 0.005) continue;
	  if(TMath::Abs((fSignTreeVariableInvMassAntiLambda - massLambda))<= 0.005) continue;
	  if(sysV0==1){   
	    if(fSignTreeVariableV0CosineOfPointingAngle <= 0.997) continue;
	  }
	  if(sysV0==2){   
	    if(fSignTreeVariablectau >=8 )continue;
	  }
	  if(sysV0==3){   
	    if(fSignTreeVariableRapK0Short >= 0.5)continue;
	  }
	  if(sysV0==4){   
	    if(TMath::Abs((fSignTreeVariableInvMassLambda - massLambda))<= 0.010) continue;
	    if(TMath::Abs((fSignTreeVariableInvMassAntiLambda - massLambda))<= 0.010) continue;
	  }
	  if(sysV0==6){   
	    if(fSignTreeVariableDcaV0ToPrimVertex>=0.3 )continue;
	  }
	}
      }

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

    }
  
    //**********************************************************************************************

    // if(TMath::Abs((fSignTreeVariableInvMassK0s - massK0s))> 0.01){
    //   hSign_PtTrigger->Fill(fSignTreeVariablePtTrigger);
    // }
   
    //    cofSignTreeVariableisPut << "*********************"<< endl;
    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();

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

    
    if(fSignTreeVariablePDGCode==310 &&  fSignTreeVariableisPrimaryV0==1){
    fHistSelectedV0PtTMaxPhi->Fill(fSignTreeVariablePtTrigger, fSignTreeVariablePhiV0, fSignTreeVariableMultiplicity);
    fHistSelectedV0PtTMaxEta->Fill(fSignTreeVariablePtTrigger, fSignTreeVariableEtaV0, fSignTreeVariableMultiplicity);
    fHistSelectedV0PtPtTMax->Fill(fSignTreeVariablePtV0,fSignTreeVariablePtTrigger , fSignTreeVariableMultiplicity);
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
	fHistPrimaryV0[5]->Fill(p,fSignTreeVariablePtV0, fSignTreeVariablePtTrigger );
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
	    hMassvsPt_SEbins[m][z][tr]->Fill(fSignTreeVariableInvMassK0s, fSignTreeVariablePtV0); 
	    if(fSignTreeVariablePDGCode==310) hMassvsPt_SEbins_true[m][z][tr]->Fill(fSignTreeVariableInvMassK0s, fSignTreeVariablePtV0); 
	  }
	  for(Int_t v=0; v<numPtV0; v++){
	    if(MoltSel && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0<NPtV0[v+1]){
	      hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      hInvMassK0Short_SEbins[m][z][v][tr]->Fill(fSignTreeVariableInvMassK0s);
	    }
	  }
	}
      }
    }
  }

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
	}
      }
    }
  }

   
  fout->Write();
  cout << "Pt Min delle particelle trigger " << PtTrigMin<< endl;
  cout << "signal pairs trigger-associated after Pt min cut " << 	CounterSignPairsAfterPtMinCut << endl;
  cout << "bkg pairs trigger-associated after Pt min cut " << 	CounterBkgPairsAfterPtMinCut << endl;
  cout << "partendo dal file " << PathIn << " ho creato il file " << PathOut<< endl;
}

