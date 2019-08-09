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
void readTreePLChiara_second(bool isMC = kFALSE,Bool_t isEfficiency=1,Int_t sysTrigger=0, Int_t sysV0=0,
			     TString year="2016l",  TString Path1 ="_15runs_6thtry", Int_t type=0,  Double_t ptjmin=3, Double_t ptjmax =30, Double_t nsigmamin=4, Double_t nsigmamax=10){

  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  const Float_t ctauK0s=2.6844;

  TString file;
  const Int_t numtipo=2;  
  TString tipo[numtipo]={"kK0s", "bo"};

  TString PathIn="./AnalysisResults";
  TString PathOut="./histo/AngularCorrelation";
  TString PathInMass= "invmass_distribution_thesis/";
  TString PathInMassDef;
  PathIn+=year;
  PathOut+=year;  
  PathInMass+=year;
  
  if(isMC && isEfficiency){ 
    PathIn+="_MCEff";
    PathOut+="_MCEff";
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
  PathOut +=Form("_SysT%i_SysV0%i",sysTrigger, sysV0); 
  PathOut+= ".root";

  TFile *fin = new TFile(PathIn);
  TFile *fileMassSigma;
  TDirectoryFile *d = (TDirectoryFile*)fin->Get("MyTask");
  TTree *tSign = (TTree *)d->Get("fSignalTree");
  TTree *tBkg  = (TTree *)d->Get("fBkgTree");

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  
 
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0, 5, 10, 30, 50, 100};
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  //  TString SPtTrigger[numPtTrigger]={"2-10"};
  Double_t NPtTrigger[numPtTrigger+1]={ptjmin,ptjmax};


  Double_t sigma[numtipo][nummolt+1][numPtV0];
  Double_t mass[numtipo][nummolt+1][numPtV0];
  Double_t meansigma;
  Double_t meanmass;
  TH1F* histoSigma;
  TH1F* histoMean;


  if(isMC==0 || (isMC==1 && isEfficiency==1)){
    for(Int_t m=0; m<nummolt+1; m++){
      PathInMassDef=PathInMass+"/invmass_distribution_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i.root", m, sysTrigger, sysV0);
      fileMassSigma= new TFile(PathInMassDef);
      histoSigma=(TH1F*)fileMassSigma->Get("histo_sigma");
      histoMean=(TH1F*)fileMassSigma->Get("histo_mean");
      for(Int_t v=0; v<numPtV0; v++){
	mass[type][m][v]=histoMean->GetBinContent(v+1);
	sigma[type][m][v]=histoSigma->GetBinContent(v+1);
	cout <<"mult interval " <<  m << " PtV0 interval " << v << " mean " << mass[type][m][v] << " sigma "<< sigma[type][m][v] << endl;
      }
    }
  }

  Double_t     fSignTreeVariablePtTrigger;		       
  Double_t     fSignTreeVariableChargeTrigger;		       
  Double_t     fSignTreeVariableEtaTrigger;		       
  Double_t     fSignTreeVariablePhiTrigger;		       
  Double_t     fSignTreeVariableDCAz;			       
  Double_t     fSignTreeVariableDCAxy;			  
  Double_t     fSignTreeVariableisPrimaryTrigger;			  
  Double_t     fSignTreeVariableisPrimaryV0;			  
  Double_t     fSignTreeVariableRapK0Short;		       
  Double_t     fSignTreeVariableDcaV0ToPrimVertex;	       
  Double_t     fSignTreeVariableDcaPosToPrimVertex;	       
  Double_t     fSignTreeVariableDcaNegToPrimVertex;	       
  Double_t     fSignTreeVariableV0CosineOfPointingAngle;	       
  Double_t     fSignTreeVariablePtV0;			       
  Double_t     fSignTreeVariablectau;			       
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


  Double_t        fBkgTreeVariablePtTrigger;
  Double_t        fBkgTreeVariableChargeTrigger;
  Double_t        fBkgTreeVariableEtaTrigger;
  Double_t        fBkgTreeVariablePhiTrigger;
  Double_t        fBkgTreeVariableDCAz;
  Double_t        fBkgTreeVariableDCAxy;
  Double_t        fBkgTreeVariableisPrimaryTrigger;			  
  Double_t        fBkgTreeVariableisPrimaryV0;			  
  Double_t        fBkgTreeVariableRapK0Short;
  Double_t        fBkgTreeVariableDcaV0ToPrimVertex;
  Double_t        fBkgTreeVariableDcaPosToPrimVertex;
  Double_t        fBkgTreeVariableDcaNegToPrimVertex;
  Double_t        fBkgTreeVariableV0CosineOfPointingAngle;
  Double_t        fBkgTreeVariablePtV0;
  Double_t        fBkgTreeVariablectau;			       
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

 
  //Signal
  tSign->SetBranchAddress("fTreeVariablePtTrigger"                 ,&fSignTreeVariablePtTrigger);		       
  tSign->SetBranchAddress("fTreeVariableChargeTrigger"             ,&fSignTreeVariableChargeTrigger);		       
  tSign->SetBranchAddress("fTreeVariableEtaTrigger"                ,&fSignTreeVariableEtaTrigger);		       
  tSign->SetBranchAddress("fTreeVariablePhiTrigger"                ,&fSignTreeVariablePhiTrigger);		       
  tSign->SetBranchAddress("fTreeVariableDCAz"                      ,&fSignTreeVariableDCAz);			       
  tSign->SetBranchAddress("fTreeVariableDCAxy"                     ,&fSignTreeVariableDCAxy);			       
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
  tBkg->SetBranchAddress("fTreeVariableZvertex"                   ,&fBkgTreeVariableZvertex);                        

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

  /*-----------------------DeltaEtaDeltaPhi in bin di molteplicita/ZVertex/pTV)/pTTrigger------------- */
  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TString namemassSE[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH2D *hDeltaEtaDeltaPhi_SEbins_sidebands[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hSign_PtTrigger[nummolt+1][numzeta];
  TH1D *hSign_PtV0[nummolt+1][numzeta];

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
	for(Int_t v=0; v<numPtV0; v++){
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
	for(Int_t v=0; v<numPtV0; v++){
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

  dirSign->cd();
  for(Int_t k = 0; k<EntriesSign; k++){
    tSign->GetEntry(k);
    if(isMC==0 || (isMC==1 && isEfficiency==1)){
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
    if(sysV0==1){   
      if(fSignTreeVariableV0CosineOfPointingAngle < 0.997) continue;
    }
    if(sysV0==2){   
      if(fSignTreeVariablectau >3*ctauK0s )continue;
    }
    if(sysV0==4){   
      if(TMath::Abs((fSignTreeVariableInvMassLambda - massLambda))< 0.010) continue;
      if(TMath::Abs((fSignTreeVariableInvMassAntiLambda - massLambda))< 0.010) continue;
    }
    if(sysV0==5){   
      if(TMath::Abs((fSignTreeVariableInvMassLambda - massLambda))< 0.005) continue;
      if(TMath::Abs((fSignTreeVariableInvMassAntiLambda - massLambda))< 0.005) continue;
    }
    if(sysV0==6){   
      if(fSignTreeVariableDcaV0ToPrimVertex>0.3 )continue;
    }
    }
 if(sysV0==3){   
      if(fSignTreeVariableRapK0Short > 0.5)continue;
    }
   
    //**********************************************************************************************

    if (fSignTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fSignTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fSignTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fSignTreeVariableDeltaPhi += 2.0*TMath::Pi();
    for(Int_t m=0; m<nummolt+1; m++){
      if(m< nummolt) BoolVar ==  fSignTreeVariableMultiplicity>=Nmolt[m] && fSignTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=0; v<numPtV0; v++){
	    if(isMC && !isEfficiency) BoolMC = kTRUE;
	    else BoolMC =TMath::Abs((fSignTreeVariableInvMassK0s - mass[type][m][v]))<3*sigma[type][m][v]; 
	    if(BoolMC){
	      hSign_PtTrigger[m][z]->Fill(fSignTreeVariablePtTrigger);
	      hSign_PtV0[m][z]->Fill(fSignTreeVariablePtV0);
	    }	  
	    if(BoolVar && fSignTreeVariablePtTrigger>=NPtTrigger[tr] && fSignTreeVariablePtTrigger<NPtTrigger[tr+1] && fSignTreeVariablePtV0>=NPtV0[v]&& fSignTreeVariablePtV0>=NPtV0[v]){
	      if(BoolMC){
		hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->Fill(fSignTreeVariableDeltaEta, fSignTreeVariableDeltaPhi);
	      }
	      if((!isMC || (isMC &&isEfficiency))&& TMath::Abs((fSignTreeVariableInvMassK0s - mass[type][m][v]))>nsigmamin*sigma[type][m][v] && TMath::Abs((fSignTreeVariableInvMassK0s - mass[type][m][v]))<nsigmamax*sigma[type][m][v]){
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
 
  cout << "ciao " << endl;
  dirBkg->cd();
  for(Int_t k = 0; k<EntriesBkg; k++){
    tBkg->GetEntry(k);
    if(isMC==0 || (isMC==1 && isEfficiency==1)){
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
    if(sysV0==1){   
      if(fBkgTreeVariableV0CosineOfPointingAngle < 0.997) continue;
    }
    if(sysV0==2){   
      if(fBkgTreeVariablectau >3*ctauK0s )continue;
    }
    if(sysV0==4){   
      if(TMath::Abs((fBkgTreeVariableInvMassLambda - massLambda))< 0.010) continue;
      if(TMath::Abs((fBkgTreeVariableInvMassAntiLambda - massLambda))< 0.010) continue;
    }
 if(sysV0==5){   
      if(TMath::Abs((fBkgTreeVariableInvMassLambda - massLambda))< 0.005) continue;
      if(TMath::Abs((fBkgTreeVariableInvMassAntiLambda - massLambda))< 0.005) continue;
    }
      
 if(sysV0==6){   
      if(fBkgTreeVariableDcaV0ToPrimVertex>0.3 )continue;
    }
    }
    if(sysV0==3){   
      if(fBkgTreeVariableRapK0Short > 0.5)continue;
    }

    //**********************************************************************************************

    if (fBkgTreeVariableDeltaPhi >  (1.5*TMath::Pi())) fBkgTreeVariableDeltaPhi -= 2.0*TMath::Pi();
    if (fBkgTreeVariableDeltaPhi < (-0.5*TMath::Pi())) fBkgTreeVariableDeltaPhi += 2.0*TMath::Pi();
    for(Int_t m=0; m<nummolt+1; m++){
      if(m< nummolt) BoolVar ==  fBkgTreeVariableMultiplicity>=Nmolt[m] && fBkgTreeVariableMultiplicity<Nmolt[m+1];
      else BoolVar=kTRUE;
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=0; v<numPtV0; v++){
	    if(isMC && !isEfficiency) BoolMC = kTRUE;
	    else BoolMC =TMath::Abs((fBkgTreeVariableInvMassK0s - mass[type][m][v]))<3*sigma[type][m][v]; 
	    if(BoolMC){
	      hBkg_PtTrigger[m][z]->Fill(fBkgTreeVariablePtTrigger);
	      hBkg_PtV0[m][z]->Fill(fBkgTreeVariablePtV0);
	    }	  
	    if(BoolVar && fBkgTreeVariablePtTrigger>=NPtTrigger[tr] && fBkgTreeVariablePtTrigger<NPtTrigger[tr+1] && fBkgTreeVariablePtV0>=NPtV0[v]&& fBkgTreeVariablePtV0>=NPtV0[v]){
	      if(BoolMC){
		hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->Fill(fBkgTreeVariableDeltaEta, fBkgTreeVariableDeltaPhi);
	      }
	      if((!isMC || (isMC &&isEfficiency)) && TMath::Abs((fBkgTreeVariableInvMassK0s - mass[type][m][v]))>nsigmamin*sigma[type][m][v] && TMath::Abs((fBkgTreeVariableInvMassK0s - mass[type][m][v]))<nsigmamax*sigma[type][m][v]){
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

  cout << "partendo dal file " << PathIn << " ho creato il file " << PathOut<< endl;
}
