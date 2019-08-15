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
#include <TLegend.h>



void AngularCorrelation(Int_t type=0, Bool_t isMC=0, Bool_t isEfficiency=0,   TString year="2016lk", TString yearMC="2018d8",  TString Path1 ="", TString Path2 ="_15runs_6thtry",Float_t ptjmin=3,  Float_t ptjmax=30, Int_t sysTrigger=0, Int_t sysV0=5, Int_t rebin=1,  Int_t rebinx=2,  Int_t rebiny=2, Float_t BulkLowValue=0.5, Float_t BulkUpValue=1.2,Float_t JetValue=0.4 ){ 

  TString file = year+Path1;
  if(isMC && isEfficiency) file = yearMC + "_MCEff" + Path1;
  if(isMC && !isEfficiency) file = yearMC + "_MCTruth" + Path1;
  TString fileMC = yearMC+"_MCEff" +Path2;
  Int_t MC=0;
  if (isMC) MC=1;

  TString PathIn= "histo/AngularCorrelation" + file + Form("_SysT%i_SysV0%i", sysTrigger, sysV0)+".root";
  TString PathInBis =  "AnalysisResults" + file  +".root";
  TString PathInEfficiency =  "histo/AngularCorrelation" +fileMC  +Form("_MC_Efficiency_SysT%i_SysV0%i", sysTrigger, sysV0)+".root";
  TString PathOut1="histo/AngularCorrelation" + file  + Form("_SysT%i_SysV0%i", sysTrigger, sysV0)+"_Output.root";
  TString PathInMass= "invmass_distribution_thesis/" + year + Path1 ;
  TString PathInMassDef;

  TFile *filepurezza;
  TFile *filein = new TFile(PathIn);
  TFile *fileinbis = new TFile(PathInBis);
  TFile *fileinEfficiency = new TFile(PathInEfficiency);

  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask");
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TList *list2 = (TList*)dir->Get("MyOutputContainer3");


  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "REMEMBER TO RUN READTREEPLCHIARA_first.C and READTREEPLCHIARA_second.C FIRST! "  << endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSB=2;
  const Int_t numeta=3;
  const Int_t numetabis=2;
  TString tipo[numtipo]={"kK0s", "bo"};
  
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Ssideband[numSB]={"", "_SB"};
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  TString SPtTrigger[numPtTrigger]={"3-30"};
  Double_t NPtTrigger[numPtTrigger+1]={ptjmin,ptjmax};
  TString Seta[numetabis]={"_eta<0.5", "_allEta"};
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Int_t numbin=2; //questo numero deve essere pari

  Float_t binwx;
  Float_t binwy;
 
  TH2D *fHistTriggervsMult         = (TH2D*)list->FindObject("fHistPtvsMultBefAll");
  TH1D *fHistTriggervsMult_MultProj= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj");
  TH2D *fHistTriggervsMultV0       = (TH2D*)list->FindObject("fHistPtvsMult");
  TH1D *fHistTriggervsMultV0_MultProj= (TH1D*)fHistTriggervsMultV0->ProjectionY("fHistTriggervsMultV0_MultProj");
  TH1D *HistoTriggerEfficiency     = (TH1D*)fileinEfficiency->Get("HistoTriggerEfficiency"); //trigger efficiency vs mult
  TH1D *HistContTriggerMolt        = (TH1D*)fileinEfficiency->Get("HistContTriggerMolt"); //trigger contamination factor vs mult

  cout << "\n \n *************************************************************************************" << endl;
  cout << "**********Sto effettuando Angular correlation in range di PtV0 e di molteplicita'********" << endl;

  Int_t Marker[nummolt+1]={7,4,20,22,29,28};
  Int_t Color[nummolt+1]={2,3,4,6,7, 9};

  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins_rapMolt[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MErap[nummolt+1][numzeta][numPtV0][numPtTrigger];
  // Float_t ssbFactor[nummolt+1][numPtV0];    
  // TH1F * histoSSB[nummolt+1];
  TH1F * histo_Bcentral[nummolt+1];
  TH1F * histo_Bside[nummolt+1];
  TH1F*  fHistPtvsMult_mult[nummolt+1];
  TH1D* HistBI[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D* HistBII[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];

  TH2D *hDeltaEtaDeltaPhi_ME_normbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB][numeta]; //SE/ME norm proiettato in delta Phi
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[nummolt+1][numzeta][numPtV0][numPtTrigger][numeta]; //sottrazione bkg V0 effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata  + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EtaAll[nummolt+1][numzeta][numPtV0][numPtTrigger];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis]; //sottrazione bkg V0 effettuata + corr efficienze + sottrazione distribuzione bulk effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis]; //istogramma finale su cui effettuo fit

  Float_t norm_MEbins[nummolt+1][numPtV0][numSB]={0}; //valore ME isto in (0,0)
  Float_t norm_MEbins_norm[nummolt+1][numPtV0][numSB]={0}; //valore ME isto normalizzato in (0,0)
  Int_t   NTrigger[nummolt+1]={0}; //total number of trigger particles 
  Int_t   NTriggerV0[nummolt+1]={0}; //total number of trigger particles only in events with V > 0
  Float_t NSpectrum[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles
  Float_t NSpectrumError[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles
  Float_t NSpectrumV0[nummolt+1][numPtV0][numetabis]={0}; //total number of trigger particles in events with V> 0
  TH1D*  fHistSpectrum[nummolt+1][numetabis];
  TH1D*  fHistSpectrumV0[nummolt+1][numetabis];
  TH1D*  fHistSpectrumV0NotEffCorr[nummolt+1][numetabis];

  TH1D*  fHistV0EfficiencyPtBins[nummolt+1]; //efficienza selezione V0 in PtV0 bins (vedi SPtV0)
  TH1D*  HistContV0PtBins[nummolt+1]; //contamination factor V0 in PtV0 bins (vedi SPtV0)

  TF1* 	  gauss[nummolt+1][numzeta][numPtV0][numPtTrigger][numetabis];

  Bool_t isProjectionPhi[nummolt+1][numzeta][numPtV0][numPtTrigger];



  //********************************************************************* 
  //**************calcolo numero particelle di trigger*******************
  //********************************************************************* 
  for(Int_t m=0; m<nummolt+1; m++){
    fHistV0EfficiencyPtBins[m]= (TH1D*)fileinEfficiency->Get("fHistV0EfficiencyPtBins_" + Smolt[m]);
    HistContV0PtBins[m]= (TH1D*)fileinEfficiency->Get("HistContV0PtBins_" + Smolt[m]);

    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m]+0.001); j<fHistTriggervsMult_MultProj->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
	NTriggerV0[m]+= fHistTriggervsMultV0_MultProj->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<fHistTriggervsMult_MultProj->GetNbinsX(); j++ ){
	NTrigger[m]+= fHistTriggervsMult_MultProj->GetBinContent(j);
	NTriggerV0[m]+= fHistTriggervsMultV0_MultProj->GetBinContent(j);
      }
    }
    cout << "n trigger in mult range (all triggers)    " << m << "  " <<  NTrigger[m] <<   endl;
    cout << "n trigger in mult range (only events V>0) " << m << "  " <<  NTriggerV0[m] <<   endl;
  }
    



  for(Int_t sb=0; sb< numSB; sb++){
    for(Int_t m=0; m<nummolt+1; m++){
      PathInMassDef=PathInMass+"/invmass_distribution_"+tipo[type]+Form("_molt%i_sysT%i_sysV0%i.root", m, sysTrigger, sysV0);
      filepurezza= new TFile(PathInMassDef);
      //    histoSSB[m]=(TH1F*)filepurezza->Get("histo_SSB");
      histo_Bcentral[m]=(TH1F*)filepurezza->Get("histo_Bcentral");
      histo_Bside[m]=(TH1F*)filepurezza->Get("histo_Bside");
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=0; v<numPtV0; v++){

	    cout << "\n\n *********************************************" << endl;
	    cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v]<< " per sideband(1) or not(0) = " << sb << endl;
	  
	    //	    cout << "ciao " << endl;
	    nameME[m][z][v][tr][sb]="ME_";
	    nameME[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];
	    nameSE[m][z][v][tr][sb]="SE_";
	    nameSE[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];

	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameSE[m][z][v][tr][sb]);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameME[m][z][v][tr][sb]);
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Rebin2D(rebinx,rebiny);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Rebin2D(rebinx,rebiny);

	    //	    cout << "\n ho preso istogrammi Me e Se" << endl;

	    /////////// primo modo per normalizzare ME 
	    binwx= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->GetBinWidth(1);
	    binwy= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->GetBinWidth(1);

	    cout << "\nbin x (delta phi) " << binwx <<  " bin y (delta eta) " << binwy << endl;
	    Int_t cont=0;
	    for(Int_t i=0; i < numbin ; i++ ){
	      for(Int_t j=0; j < numbin ; j++ ){
		norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001+(i-numbin/2)*binwx), hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001+(j-numbin/2)*binwy));
		cont++;
	      }
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/cont;  	    
  
	    cout << "\nvalore histo in DeltaEta=0 e DeltaPhi=0 (quadrato)" <<  norm_MEbins[m][v][sb] << endl;

	    /////////// secondo modo per normalizzare ME 
	    norm_MEbins[m][v][sb]=0;
	    cont=0;
	    for(Int_t i=0; i <hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetNbinsY()  ; i++ ){
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001),i+1);
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-0.000001),i+1);
	      cont++;
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/2/cont;
	    
	    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 (media in delta phi, 2 righe) " <<  norm_MEbins[m][v][sb] << endl;

	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_norm");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC");
	    if (norm_MEbins[m][v][sb]==0){
	      cout << "normalizzazione ME non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	      continue;  
	    }

	    //	    cout << "\n divido SE distribution per ME disatribution normalized" << endl;
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]-> Scale(1./norm_MEbins[m][v][sb]);
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]);

	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->SetTitle("ME distribution normalized at 1 in (0,0) in mult. " +Smolt[m] + " and Pt V0 " +  SPtV0[v]+ Ssideband[sb]);
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->SetTitle("SE distribution normalized for pair acceptance in mult " +Smolt[m] + " and  Pt V0 " + SPtV0[v]+ Ssideband[sb]);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta #Eta");
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta #Phi");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta #Eta");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta #Phi");
	    norm_MEbins_norm[m][v][sb]= hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001), hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001));
	    cout << "valore histo  in DeltaEta=0 e DeltaPhi=0, histo normalizzato  " <<  norm_MEbins_norm[m][v][sb] <<"\n"<< endl;

	    //*****************************************************************************************************************
	    //proietto in Delta Phi
	    //****************************************************************************************************************

	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][0]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaJet", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-JetValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(JetValue), "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][2]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaAll",0,-1, "E"));
	    // HistBI->Sumw2();
	    // HistBII->Sumw2();

	    HistBI[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaBulkI", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-BulkUpValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-BulkLowValue), "E"));
	    HistBII[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_etaBulkII", hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(BulkLowValue), hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetXaxis()->FindBin(BulkUpValue), "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]= (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][0]->Clone(nameME[m][z][v][tr][sb]+"_AC_phi_etaBI");

	    for(Int_t j=1; j < HistBI[m][z][v][tr][sb]->GetNbinsX(); j++){
	      hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->SetBinContent(j, HistBI[m][z][v][tr][sb]->GetBinContent(j)+ HistBII[m][z][v][tr][sb]->GetBinContent(j));
	    }

	    for(Int_t eta=0; eta<numeta; eta++){
	      //	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetName(nameME[m][z][v][tr][sb]+"_AC_phi");
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->SetTitle("#Delta #Phi projected angular correlation in mult. " +Smolt[m] + " and Pt V0 " + SPtV0[v]+ Ssideband[sb]);
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][eta]->GetXaxis()->SetTitle("#Delta #Phi");
	    Bool_t NotScale=kFALSE;
	    if (sb==1){
	      if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) NotScale=kTRUE;
	      if(!NotScale){
		hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->Scale((Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001))/histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001)));
	      }
	      cout << "Bcentral " << (Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001)) << " Bside " << histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001)) << " B central/Bside " << (Float_t)histo_Bcentral[m]->GetBinContent(histo_Bcentral[m]->FindBin(NPtV0[v]+0.0001))/histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))<< endl;
	    }
	    }//chiusura ciclo eta


	  
	  } //chiusura ciclo su pt della V0
	}//chiusura ciclo su Pt Trigger
      }//chiusura ciclo su z vertice
	
    }//chiusura ciclo su molteplcità
  }//chiusura ciclo su sideband or not


  cout << " \n\n ********** second part of the program *******************" << endl;
  for(Int_t m=0; m<nummolt+1; m++){
    for(Int_t etabis=0; etabis< numetabis; etabis++){
      fHistSpectrum[m][etabis]=new TH1D ("fHistSpectrum_"+Smolt[m]+Seta[etabis], "fHistSpectrum_"+Smolt[m]+Seta[etabis], numPtV0, NPtV0) ;
    fHistSpectrumV0[m][etabis]=new TH1D ("fHistSpectrumV0_"+Smolt[m]+Seta[etabis], "fHistSpectrumV0_"+Smolt[m]+Seta[etabis], numPtV0, NPtV0) ;
    }
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=0; v<numPtV0; v++){
    for(Int_t etabis=0; etabis< numetabis; etabis++){
	    fHistSpectrum[m][etabis]->SetBinContent(v+1, 0);
	  fHistSpectrumV0[m][etabis]->SetBinContent(v+1,0);
    }
	  cout << "\n\n *********************************************" << endl;
	  cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v] << endl;
	  hDeltaEtaDeltaPhi_MErap[m][z][v][tr]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][0]->Clone("ME_m"+Smolt[m]+"_v"+SPtV0[v]+"_rap");
	  hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Divide(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][1], hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][0]);
	  hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->SetTitle("Rapporto ME sideband/ME 3sigma region");
	  if (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0){
	      cout << "normalizzazione ME non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	      continue;  
	    }
	  //***************************************************************
	  //sottraggo fondo data da "finte" K0s 
	  //***************************************************************
	  for(Int_t eta=0; eta< numeta; eta++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta] =  (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->Clone(nameME[m][z][v][tr][0]+Form("%i_AC_phi_V0Sub", eta));
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->SetName(nameME[m][z][v][tr][0]+Form("_eta%i_AC_phi_V0Sub",eta));
	    for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->SetBinContent(i,  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0][eta]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1][eta]->GetBinContent(i));
	  }
	  //}
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetXaxis()->SetTitle("#Delta #Phi");
	  //cout << " \n \n error of phi projection raw and Vo bkg subtracted " << endl;
	  for(Int_t i=1; i<hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][eta]->GetNbinsX(); i++ ){
	    //the erros below are exactly the same
	    //   cout << "error of bin " << i << ": " << 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinError(i)<< endl;
	    //   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->SetBinError(i, sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0]->GetBinError(i),2) +pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1]->GetBinError(i),2)));
	    //   cout <<"error of bin " << i << " my calculation: "<<     hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinError(i)<< endl;
	  }
	  }
	  //***************************************************************
	  //sottraggo distribuzione del bulk 
	  //***************************************************************
	  Float_t ScaleFactor=(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue)))/2./(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)));
	  cout << (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(JetValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(-JetValue))) << "   " <<(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue)) - hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue)))<<"  " << (hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinUpEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkUpValue))) << "  " << hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->GetBinLowEdge(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][0]->GetXaxis()->FindBin(BulkLowValue));
	  cout << "\n" <<ScaleFactor << endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Scale(ScaleFactor);
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EtaAll[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EtaAll");
	  for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetNbinsX(); i++){
	 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->SetBinContent(i,  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->GetBinContent(i));
	  }
	 

	  //****************************************************************************************************************
	  //divido angular correlation proietatta in deltaphi per efficienza di selezione Trigger e V0 e per contamination factors e per SSB
	  //****************************************************************************************************************
	  cout << "\n \n Trigger selection efficiency: " <<HistoTriggerEfficiency->GetBinContent(m+1) << "\n V0 selection efficiency:  " << fHistV0EfficiencyPtBins[m]->GetBinContent(v+1) << "\n Trigger contamination factor:  " << (1-HistContTriggerMolt->GetBinContent(m+1)) << "\n V0 contamination factor:  " << (1-HistContV0PtBins[m]->GetBinContent(v+1))<< endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][0]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_BulkSub_EffCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][1]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_EtaAll[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_EtaAll_EffCorr");
	  if((HistoTriggerEfficiency->GetBinContent(m+1))!=0 && (fHistV0EfficiencyPtBins[m]->GetBinContent(v+1))!=0 ){
	    for(Int_t etabis=0; etabis<2; etabis++){
	      // hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->Scale(1./HistoTriggerEfficiency->GetBinContent(m+1)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContTriggerMolt->GetBinContent(m+1))*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	      	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->Scale(1./fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
	    }
	    isProjectionPhi[m][z][v][tr]=kTRUE;
	  }
	  else {
	    cout << "\n ****** efficienza pari a zero********* " << endl;
	    isProjectionPhi[m][z][v][tr]=kFALSE;
	    continue;
	  }

	  //****************************************************************************************************************
	  //bkg subtraction
	  //****************************************************************************************************************
	  for(Int_t etabis=0; etabis<2; etabis++){
	  Float_t DeltaPhiBkg=0;
	  Float_t counter=0;
	  Float_t DeltaPhiBkgError=0;

	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-1.1); i++){
	    counter++;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i);
	  }
	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.84); i++){
	    counter++;
	    //	    cout <<" content " <<  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i) << endl;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i);
	  }
	  DeltaPhiBkg=DeltaPhiBkg/counter;

	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(-1.1); i++){
	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i)-DeltaPhiBkg,2);
	  }
	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->FindBin(1.84); i++){//simmetrico attorno a pi/mezzi

	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i)-DeltaPhiBkg,2);
	  }
	  DeltaPhiBkgError=sqrt(DeltaPhiBkgError/counter/(counter+1));


	  cout << "numero di bin su cui si media per sottrarre fondo: " << counter << "\n valore del fondo di coppie non correlate: " << DeltaPhiBkg << " +- " << DeltaPhiBkgError<<  endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]= (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_BkgSub"+Seta[etabis]);
	  //hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Sumw2();

	  if(etabis==1){
	  for(Int_t i =1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetNbinsX(); i++){
	      hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->SetBinContent(i, (Float_t)(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinContent(i))-DeltaPhiBkg);
	    //the errors below are exactly the same
//	    cout << "\nbin n: " << i << endl;
//	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->GetBinError(i)<< endl;
//	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->GetBinError(i)<< endl;
	  }//
	  }

	  gauss[m][z][v][tr][etabis] = new TF1("gaus_m"+ Smolt[m]+"_v"+SPtV0[v]+Form("_eta%i",etabis),"gaus",-0.5*TMath::Pi(), 1.5*TMath::Pi());
	  gauss[m][z][v][tr][etabis]->SetLineColor(kRed);   
	  gauss[m][z][v][tr][etabis]->SetParameter(1,0);
	  gauss[m][z][v][tr][etabis]->SetParName(0, "norm");
	  gauss[m][z][v][tr][etabis]->SetParName(1, "mean");
	  gauss[m][z][v][tr][etabis]->SetParName(2, "sigma");
	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Rebin(rebin);
	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]-> Fit(gauss[m][z][v][tr][etabis], "", "", -1, 1); //non ho ben capito l;a differenza tra 0 e N
	  NSpectrum[m][v][etabis]=gauss[m][z][v][tr][etabis]->Integral(-0.5*TMath::Pi(), 1.5*TMath::Pi())/hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->GetBinWidth(1);
	  NSpectrumError[m][v][etabis]=gauss[m][z][v][tr][etabis]->IntegralError(-0.5*TMath::Pi(), 1.5*TMath::Pi())/hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->GetBinWidth(1);
	 
	  if (NTriggerV0[m]==0)	    cout << "\n \n ERROR: number of trigger particles is zero! " << endl;
	  if (NTriggerV0[m]!=0){ 
	    NSpectrumV0[m][v][etabis]= NSpectrum[m][v][etabis]/NTriggerV0[m];
	  }
	  
	  if (NTrigger[m]==0)	    cout << "\n \n ERROR: number of trigger particles is zero! " << endl;
	  if (NTrigger[m]!=0){ 
	    NSpectrum[m][v][etabis]= NSpectrum[m][v][etabis]/NTrigger[m];
	    NSpectrumError[m][v][etabis]= NSpectrumError[m][v][etabis]/NTrigger[m];
	  }
	  cout << "integrale/bin width/Ntrigger " << NSpectrum[m][v][etabis] << " n trigger " << NTrigger[m]<< endl;
	  cout << "integrale/bin width/Ntrigger" << NSpectrumV0[m][v][etabis] << " n trigger " << NTriggerV0[m]<< endl;
	  fHistSpectrum[m][etabis]->SetBinContent(v+1, NSpectrum[m][v][etabis]);
	  fHistSpectrum[m][etabis]->SetBinError(v+1, NSpectrumError[m][v][etabis]);
	  fHistSpectrumV0[m][etabis]->SetBinContent(v+1, NSpectrumV0[m][v][etabis]);
	  }
	  
	}
      }
    }
  }



  for(Int_t sb=0; sb< numSB; sb++){
    for(Int_t m=0; m<nummolt+1; m++){
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=0; v<numPtV0; v++){

    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Clone(nameME[m][z][v][tr][sb] + "_rapMolt");
    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->SetName(nameME[m][z][v][tr][sb] + "_rapMolt");
    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->SetTitle("Rapporto ME molteplicità indicata / ME molteplicità integrata");
       hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->Divide( hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_MEbins[5][z][v][tr][sb]);
  }
	}
      }
    }
  }

  auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->SetHeader("Selezioni applicate");     

  cout << "\n \n sto per fare la canvas dello spettro " << endl;
  TCanvas *canvasSpectrum1= new TCanvas ("canvasspectrum1", "spettro ottenuto da |deltaEta| < 0.5");
  canvasSpectrum1->cd();
  for(Int_t m=0; m< nummolt+1; m++){
    //    cout << m << endl;   
    //fHistSpectrum[m]->GetYaxis()->SetRangeUser(0,60);
    fHistSpectrum[m][0]->SetMarkerStyle(Marker[m]);
    fHistSpectrum[m][0]->SetLineColor(Color[m]);
    fHistSpectrum[m][0]->SetMarkerColor(Color[m]);
    legend->AddEntry(fHistSpectrum[m][0],Smolt[m],"pel");   
    fHistSpectrum[m][0]->Draw("samee");
    if(m==nummolt) legend->Draw();
    }

  TCanvas *canvasSpectrum2= new TCanvas ("canvasspectrum2", "spettro ottenuto da EtaAll");
  canvasSpectrum2->cd();
  for(Int_t m=0; m< nummolt+1; m++){
    //    cout << m << endl;   
    //fHistSpectrum[m]->GetYaxis()->SetRangeUser(0,60);
    fHistSpectrum[m][1]->SetMarkerStyle(Marker[m]);
    fHistSpectrum[m][1]->SetLineColor(Color[m]);
    fHistSpectrum[m][1]->SetMarkerColor(Color[m]);
    fHistSpectrum[m][1]->Draw("samee");
    if(m==nummolt) legend->Draw();
    }


  cout << "\n sto per scrivere su file " << endl;
  // salvo solo istogrammi con norm != 0 sia per sideband che per regione centrale, con efficienze diverse da zero e con integrale siiideband non =0 (ossia solo histo con significato fisico) 
  TFile *fileout = new TFile(PathOut1, "RECREATE");
  for(Int_t m=0; m<nummolt+1; m++){
    //    cout << "\n" <<m << endl;
    if (NTrigger[m]!=0) fHistSpectrum[m][0]->Write();
    if (NTrigger[m]!=0) fHistSpectrum[m][1]->Write();
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=0; v<numPtV0; v++){
	  if (histo_Bside[m]->GetBinContent(histo_Bside[m]->FindBin(NPtV0[v]+0.0001))==0) continue;
	  hDeltaEtaDeltaPhi_MErap[m][z][v][tr]->Write();
	    if (norm_MEbins[m][v][0]==0 || norm_MEbins[m][v][1]==0) continue;
	    if(!isProjectionPhi[m][z][v][tr]) continue;
	   
	  for(Int_t sb=0; sb<numSB; sb++){
		
	    hDeltaEtaDeltaPhi_MEbins_rapMolt[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][0]->Write();
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][1]->Write();
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb][2]->Write();
	    //   cout<<m << endl;
	  }
	   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][0]->Write();
	   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][1]->Write();
	   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr][2]->Write();
	   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub[m][z][v][tr]->Write();
	 
	   for(Int_t etabis=0; etabis < numetabis; etabis++){

	   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_BulkSub_EffCorr[m][z][v][tr][etabis]->Write();
	   hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr][etabis]->Write();
	   }
	}

      }
    }
  }


  fileout->Close(); 

  
   cout << "******************************************************************"<< endl;
   cout << "partendo dai file " << PathIn << " e "<< PathInBis << " e "<< PathInEfficiency <<  " ho creato: "<< endl;
   cout << "il file " << PathOut1 << endl;

}
