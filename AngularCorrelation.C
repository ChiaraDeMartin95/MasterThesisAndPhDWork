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

void AngularCorrelation(Int_t type=0, Int_t MC=0, Float_t ptjmin=3,  Float_t ptjmax=30, Int_t sysTrigger=0, Int_t sysV0=0){ //MC = 1 se si tratta del MC
  cout << "Remember to check the input file name " << endl;  

  TString file = "2016l_fifth";
  //  TString fileMC="2018d8_MC_try1";
  TString fileMC="2018d8_MC_try7_10runs";
  if(MC==1) file=fileMC;


  TString isMC[2]= {"", "_MC"}; 
  TString PathIn= "histo/AngularCorrelation" + file + isMC[MC]+ ".root";
  TString PathInBis =  "AnalysisResults" + file +isMC[MC] +".root";
  TString PathInEfficiency =  "histo/AngularCorrelation" + fileMC +"_MC_Efficiency.root";
  TString PathOut1="histo/AngularCorrelation" + file + isMC[MC] + "_Output.root";

  TFile *filepurezza;

  TFile *filein = new TFile(PathIn);
  TFile *fileinbis = new TFile(PathInBis);
  TFile *fileinEfficiency = new TFile(PathInEfficiency);

  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask");
  TList *list = (TList*)dir->Get("MyOutputContainer");


  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "REMEMBER TO RUN READTREEPLCHIARA.C FIRST! "  << endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;
  cout << "********************************************"<< endl;

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSB=2;
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
  Double_t massK0s = 0.497611;
  Double_t massLambda = 1.115683;
  Int_t numbin=4; //questo numero deve essere pari

  Float_t binwx;
  Float_t binwy;
  Int_t cont=0;

  TH1D *fHistTriggervsMult         = (TH1D*)list->FindObject("fHistTriggervsMult");
  TH1D *HistoTriggerEfficiency     = (TH1D*)fileinEfficiency->Get("HistoTriggerEfficiency"); //trigger efficiency vs muly
  TH1D *HistContTriggerMolt        = (TH1D*)fileinEfficiency->Get("HistContTriggerMolt"); //trigger contamination factor vs mult

  /*
    cout << endl;
    cout << "*************************************************************************************" << endl;
    cout << "**********Sto effettuando Angular correlation nell'intero range di PtV0 e di molteplicita'*******" << endl;
    TH2D *hDeltaEtaDeltaPhi_SE       = (TH2D*)filein->Get("hDeltaEtaDeltaPhi_SE");
    TH2D *hDeltaEtaDeltaPhi_ME       = (TH2D*)filein->Get("hDeltaEtaDeltaPhi_ME");

    TH2D *fHistPtvsMult              = (TH2D*)list->FindObject("fHistPtvsMult");//to get number of trigger particles in specific pt bin
    TH1D *HistoTriggerEfficiency     = (TH1D*)fileinEfficiency->Get("HistoTriggerEfficiency"); //trigger efficiency vs muly
    TH1D *HistContTriggerMolt        = (TH1D*)fileinEfficiency->Get("HistContTriggerMolt"); //trigger contamination factor vs mult

    //ecco come prendere content di un bin di un TH2
    /////////// primo modo per normalizzare ME 
    Float_t binwx= hDeltaEtaDeltaPhi_ME->GetXaxis()->GetBinWidth(1);
    Float_t binwy= hDeltaEtaDeltaPhi_ME->GetYaxis()->GetBinWidth(1);
 
    Float_t norm_ME=0;
    
    cout << "bin x (delta phi) " << binwx <<  "bin y (delta eta) " << binwy << endl;
    for(Int_t i=0; i < numbin ; i++ ){
    for(Int_t j=0; j < numbin ; j++ ){
    norm_ME+= hDeltaEtaDeltaPhi_ME->GetBinContent(hDeltaEtaDeltaPhi_ME->GetXaxis()->FindBin(0.000001+(i-numbin/2)*binwx), hDeltaEtaDeltaPhi_ME->GetYaxis()->FindBin(0.000001+(j-numbin/2)*binwy));
    cont++;
    }
    }
    norm_ME=norm_ME/cont;
    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 " <<  norm_ME << endl;
  

    /////////// secondo modo per normalizzare ME 
    norm_ME=0;
    cont=0;
    for(Int_t i=0; i <hDeltaEtaDeltaPhi_ME->GetNbinsY()  ; i++ ){
    norm_ME+= hDeltaEtaDeltaPhi_ME->GetBinContent(hDeltaEtaDeltaPhi_ME->GetXaxis()->FindBin(0.000001),i+1);
    cont++;
    }
    norm_ME=norm_ME/cont;
    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 " <<  norm_ME << endl;

    TH2D* hDeltaEtaDeltaPhi_ME_norm=(TH2D*)hDeltaEtaDeltaPhi_ME-> Clone("hDeltaEtaDeltaPhi_ME_norm"); //normalized ME distribution
    TH2D* hDeltaEtaDeltaPhi_AC=(TH2D*)hDeltaEtaDeltaPhi_ME-> Clone("hDeltaEtaDeltaPhi_AC"); // angular correlation distribution
    hDeltaEtaDeltaPhi_ME_norm->SetTitle("ME distribution normalized at 1 in (0,0)");
    hDeltaEtaDeltaPhi_AC->SetTitle("SE distribution normalized for pair acceptance");
    hDeltaEtaDeltaPhi_AC->GetXaxis()->SetTitle("#Delta #eta");
    hDeltaEtaDeltaPhi_AC->GetYaxis()->SetTitle("#Delta #phi (rad)");
    hDeltaEtaDeltaPhi_AC->GetXaxis()->SetTitleOffset(1.);
    hDeltaEtaDeltaPhi_AC->GetYaxis()->SetTitleOffset(1.);

    Int_t NTriggerTotal=0;  
    //Int_t NTriggerTotal=fHistPtvsMult->GetEntries();
  
    for(Int_t j=0; j<fHistTriggervsMult->GetNbinsX(); j++ ){
    NTriggerTotal+= fHistTriggervsMult->GetBinContent(j);
    }
  
  
    cout << "total number of trigger particles " << NTriggerTotal << endl;

  
    hDeltaEtaDeltaPhi_SE->Sumw2();
    hDeltaEtaDeltaPhi_ME->Sumw2();
    hDeltaEtaDeltaPhi_ME_norm->Sumw2();
    if(norm_ME=!0)  hDeltaEtaDeltaPhi_ME_norm-> Scale(1./norm_ME);

    //error check...

    else cout << "normalizzazione ME non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;

    hDeltaEtaDeltaPhi_AC->Divide(  hDeltaEtaDeltaPhi_SE,   hDeltaEtaDeltaPhi_ME);
						  
    Int_t norm_ME_bis= hDeltaEtaDeltaPhi_ME_norm->GetBinContent(hDeltaEtaDeltaPhi_ME_norm->GetXaxis()->FindBin(0.000001), hDeltaEtaDeltaPhi_ME_norm->GetYaxis()->FindBin(0.000001));
    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 , histo normalizzato " <<  norm_ME_bis << endl;

  */
  cout << endl;
  cout << endl;
  cout << "*************************************************************************************" << endl;
  cout << "**********Sto effettuando Angular correlation in range di PtV0 e di molteplicita'********" << endl;

  Int_t Marker[nummolt+1]={7,4,20,22,29,28};
  Int_t Color[nummolt+1]={2,3,4,6,8, 9};

  TString nameSE[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_SEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];

  Float_t ssbFactor[nummolt+1][numPtV0];    
  TH1F * histoSSB[nummolt+1];
  TH1F*  fHistPtvsMult_mult[nummolt+1];

  TString nameME[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_MEbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];

  TH2D *hDeltaEtaDeltaPhi_ME_normbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH2D *hDeltaEtaDeltaPhi_ACbins_NotCorr[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB];
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_NotCorr[nummolt+1][numzeta][numPtV0][numPtTrigger][numSB]; //non corretto per efficienze e contaminazioni
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg coppie non correlate effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_NotCorr[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg coppie non correlate effettuata
  TH1D *hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[nummolt+1][numzeta][numPtV0][numPtTrigger]; //sottrazione bkg V0 effettuata

  Float_t norm_MEbins[nummolt+1][numPtV0][numSB]={0}; //valore ME isto in (0,0)
  Float_t norm_MEbins_norm[nummolt+1][numPtV0][numSB]={0}; //valore ME isto normalizzato in (0,0)
  Int_t NTrigger[nummolt+1]={0}; //total number of trigger particles 
  Float_t NSpectrum[nummolt+1][numPtV0]={0}; //total number of trigger particles

  TH1D*  fHistV0EfficiencyPtBins[nummolt+1]; //efficienza selezione V0 in PtV0 bins (vedi SPtV0)
  TH1D*  HistContV0PtBins[nummolt+1]; //contamination factor V0 in PtV0 bins (vedi SPtV0)
  TH1D*  fHistSpectrum[nummolt+1];


  TF1* 	  gauss[nummolt+1][numzeta][numPtV0][numPtTrigger];

  Bool_t isProjectionPhi[nummolt+1][numzeta][numPtV0][numPtTrigger];
 
  for(Int_t m=0; m<nummolt; m++){
    fHistV0EfficiencyPtBins[m]= (TH1D*)fileinEfficiency->Get("fHistV0EfficiencyPtBins_" + Smolt[m]);
    HistContV0PtBins[m]= (TH1D*)fileinEfficiency->Get("HistContV0PtBins_" + Smolt[m]);
    if(m<nummolt){
      for(Int_t j=fHistTriggervsMult->GetXaxis()->FindBin(Nmolt[m]+0.001); j<fHistTriggervsMult->GetXaxis()->FindBin(Nmolt[m+1]-0.001); j++ ){
	NTrigger[m]+= fHistTriggervsMult->GetBinContent(j);
      }
    }
    else {
      for(Int_t j=1; j<fHistTriggervsMult->GetNbinsX(); j++ ){
	NTrigger[m]+= fHistTriggervsMult->GetBinContent(j);
      }
    }
    /*
      for(Int_t j=      fHistPtvsMult_mult[m]->FindBin(ptjmin); j =fHistPtvsMult_mult[m]->FindBin(ptjmax); j++ ){
      NTrigger[m] +=      fHistPtvsMult_mult[m]->GetBinContent(j);
      }
    */  
    cout << "n trigger in mult range " << m << "  " <<  NTrigger[m] <<   endl;
  }
    
  //per ora non integrato in molt
  for(Int_t sb=0; sb< numSB; sb++){
    for(Int_t m=0; m<nummolt; m++){
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=0; v<numPtV0; v++){

	    cout << "\n\n *********************************************" << endl;
	    cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v]<< " per sb = " << sb << endl;
	  
	    //	    cout << "ciao " << endl;
	    nameME[m][z][v][tr][sb]="ME_";
	    nameME[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];
	    nameSE[m][z][v][tr][sb]="SE_";
	    nameSE[m][z][v][tr][sb]+="m"+ Smolt[m]+"_v"+SPtV0[v]+Ssideband[sb];

	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameSE[m][z][v][tr][sb]);
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]= (TH2D*)filein-> Get(nameME[m][z][v][tr][sb]);

	    cout << "\n ho preso istogrammi Me e Se" << endl;

	    /////////// primo modo per normalizzare ME 
	    binwx= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->GetBinWidth(1);
	    binwy= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->GetBinWidth(1);

  
	    cout << "bin x (delta phi) " << binwx <<  "bin y (delta eta) " << binwy << endl;
	    cont=0;
	    for(Int_t i=0; i < numbin ; i++ ){
	      for(Int_t j=0; j < numbin ; j++ ){
		norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001+(i-numbin/2)*binwx), hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001+(j-numbin/2)*binwy));
		cont++;
	      }
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/cont;  
	    
  
	    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 " <<  norm_MEbins[m][v][sb] << endl;

	    /////////// secondo modo per normalizzare ME 
	    norm_MEbins[m][v][sb]=0;
	    cont=0;
	    for(Int_t i=0; i <hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetNbinsY()  ; i++ ){
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001),i+1);
	      norm_MEbins[m][v][sb]+= hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->GetXaxis()->FindBin(-0.000001),i+1);
	      cont++;
	    }
	    norm_MEbins[m][v][sb]=norm_MEbins[m][v][sb]/2/cont;
	    
	    cout << "valore histo in DeltaEta=0 e DeltaPhi=0 " <<  norm_MEbins[m][v][sb] << endl;

	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_norm");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC");
	    if (norm_MEbins[m][v][sb]==0){
	      cout << "normalizzazione ME non effettuata perche' denominatore =0; salvo histo non normalizzato " << endl;
	      continue;  
	    }

	    hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_MEbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->Sumw2();
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]-> Scale(1./norm_MEbins[m][v][sb]);
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Divide(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr][sb], hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]);
	    //	    cout << "ciao " << endl;
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->SetTitle("ME distribution normalized at 1 in (0,0) in mult. " +Smolt[m] + " and Pt V0 " +  SPtV0[v]+ Ssideband[sb]);
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->SetTitle("SE distribution normalized for pair acceptance in mult " +Smolt[m] + " and  Pt V0 " + SPtV0[v]+ Ssideband[sb]);
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta #Eta");
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta #Phi");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta #Eta");
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->GetYaxis()->SetTitle("#Delta #Phi");
	    //	    cout << "ciao " << endl;
	    hDeltaEtaDeltaPhi_ACbins_NotCorr[m][z][v][tr][sb]=(TH2D*)hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]-> Clone(nameME[m][z][v][tr][sb]+"_AC_NotCorr");//isto non diviso per efficienze e contamination factor
	    //	    cout << "ciao " << endl;
	    //hDeltaEtaDeltaPhi_ACbins_phi_NotCorr[m][z][v][tr][sb]->Sumw2();
	    //hDeltaEtaDeltaPhi_ACbins_phi_NotCorr[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi_NotCorr", 0, -1, "E"));
	    //	    cout << "ciao " << endl;
	    norm_MEbins_norm[m][v][sb]= hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetBinContent(hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetXaxis()->FindBin(0.000001), hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->GetYaxis()->FindBin(0.000001));
	    cout << "valore histo  in DeltaEta=0 e DeltaPhi=0, histo normalizzato  " <<  norm_MEbins_norm[m][v][sb] <<"\n"<< endl;

	    //*****************************************************************************************************************
	    //proietto in Delta Phi
	    //****************************************************************************************************************

	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb]= (TH1D*)(hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->ProjectionY(nameME[m][z][v][tr][sb]+"_AC_phi", 0, -1, "E"));
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb]->SetName(nameME[m][z][v][tr][sb]+"_AC_phi");
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb]->SetTitle("#Delta #Phi projected angular correlation in mult. " +Smolt[m] + " and Pt V0 " + SPtV0[v]+ Ssideband[sb]);
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb]->GetXaxis()->SetTitle("#Delta #Phi");
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb]->Sumw2();

	  

	  } //chiusura ciclo su pt della V0
	}//chiusura ciclo su Pt Trigger
      }//chiusura ciclo su z vertice
	
    }//chiusura ciclo su molteplcità
  }//chiusura ciclo su sideband or not

  for(Int_t m=0; m<nummolt; m++){
    filepurezza= new TFile(Form("invmass_distribution_thesis/invmass_distribution_"+ tipo[type]+"_%i.root",m), "");
    histoSSB[m]=(TH1F*)filepurezza->Get("histo_SSB");
    fHistSpectrum[m]=new TH1D ("fHistSpectrum_"+Smolt[m], "fHistSpectrum_"+Smolt[m], numPtV0, NPtV0) ;
    for(Int_t z=0; z<numzeta; z++){
      for(Int_t tr=0; tr<numPtTrigger; tr++){
	for(Int_t v=0; v<numPtV0; v++){
	  cout << "\n\n *********************************************" << endl;
	  cout << "analisi nell'intervallo di molt " << Smolt[m] <<" e nell\'intervallo di PtV0 " << SPtV0[v] << endl;
	  
	  //***************************************************************
	  //sottraggo fondo data da "finte" K0s 
	  //***************************************************************
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr] =  (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->SetName(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->Sumw2();
	  for(Int_t i=1; i< 	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->SetBinContent(i,  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0]->GetBinContent(i) -  hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1]->GetBinContent(i));
	  }
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetXaxis()->SetTitle("#Delta #Phi");
	  cout << " \n \n error of phi projection raw and Vo bkg subtracted " << endl;
	  for(Int_t i=1; i<hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetNbinsX(); i++ ){
	    //the erros below are exactly the same
	    //   cout << "error of bin " << i << ": " << 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinError(i)<< endl;
	    //   hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->SetBinError(i, sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][0]->GetBinError(i),2) +pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][1]->GetBinError(i),2)));
	    //   cout <<"error of bin " << i << " my calculation: "<<     hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinError(i)<< endl;
	  }
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_NotCorr[m][z][v][tr]=(TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_NotCorr");
	  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_NotCorr[m][z][v][tr]->SetName(nameME[m][z][v][tr][0]+"_AC_phi_V0Sub_NotCorr");
	  //****************************************************************************************************************
	  //divido angular correlation proietatta in deltaphi per efficienza di selezione Trigger e V0 e per contamination factors e per SSB
	  //****************************************************************************************************************
	  cout << "\n \n Trigger selection efficiency: " <<HistoTriggerEfficiency->GetBinContent(m+1) << "\n V0 selection efficiency:  " << fHistV0EfficiencyPtBins[m]->GetBinContent(v+1) << "\n Trigger contamination factor:  " << (1-HistContTriggerMolt->GetBinContent(m+1)) << "\n V0 contamination factor:  " << (1-HistContV0PtBins[m]->GetBinContent(v+1))<< endl;
	  if((HistoTriggerEfficiency->GetBinContent(m+1))!=0 && (fHistV0EfficiencyPtBins[m]->GetBinContent(v+1))!=0 ){
	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->Scale(1./HistoTriggerEfficiency->GetBinContent(m+1)/fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)*(1-HistContTriggerMolt->GetBinContent(m+1))*(1-HistContV0PtBins[m]->GetBinContent(v+1)));
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
	  Float_t DeltaPhiBkg=0;
	  Float_t counter=0;
	  Float_t DeltaPhiBkgError=0;

	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(-1.1); i++){
	    counter++;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinContent(i);
	  }
	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(1.84); i++){
	    counter++;
	    cout <<" content " <<  hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinContent(i) << endl;
	    DeltaPhiBkg+=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinContent(i);
	  }
	  DeltaPhiBkg=DeltaPhiBkg/counter;

	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(-0.5*TMath::Pi()); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(-1.1); i++){
	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinContent(i)-DeltaPhiBkg,2);
	  }
	  for(Int_t i=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(1.3); i < hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->FindBin(1.84); i++){//simmetrico attorno a pi/mezzi

	    DeltaPhiBkgError+=pow(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinContent(i)-DeltaPhiBkg,2);
	  }
	  DeltaPhiBkgError=sqrt(DeltaPhiBkgError/counter/(counter+1));


	  cout << "numero di bin su cui si media per sottrarre fondo: " << counter << "\n valore del fondo di coppie non correlate: " << DeltaPhiBkg << " +- " << DeltaPhiBkgError<<  endl;
	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]= (TH1D*)hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->Clone(nameME[m][z][v][tr][0]+"_AC_phi_BkgSub");
	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->Sumw2();

	  for(Int_t i =1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetNbinsX(); i++){
	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->SetBinContent(i, (Float_t)(hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinContent(i))-DeltaPhiBkg);
	    //the errors below are exactly the same
//	    cout << "\nbin n: " << i << endl;
//	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetBinError(i)<< endl;
//	    cout << "errore V0 bkg "<< 	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->GetBinError(i)<< endl;
	  }//
	  /*	  
	  Int_t sigmaAC=0;
	  for(Int_t i =1; i<=hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->GetNbinsX(); i++){
	   	  
	    //  cout << "\n \n \n sto analizzando il bin " << i << "-esimo della distribuzione proiettata in phi" << endl;
	    // for(Int_t eta=1; eta<= hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr]->GetNbinsY() ; eta++){
	    //	      sigmaAC+=pow(norm_MEbins_norm[m][v]/hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(eta, i),2)*hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(eta, i) + pow(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(eta, i)/hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(eta, i),2)*pow(sigmanorm[m][v],2) + pow(hDeltaEtaDeltaPhi_SEbins[m][z][v][tr]->GetBinContent(eta, i),2)/pow(hDeltaEtaDeltaPhi_MEbins[m][z][v][tr]->GetBinContent(eta, i),3)*pow(sigmanorm[m][v],2);

	    //}
	    sigmaAC=pow(hDeltaEtaDeltaPhi_ACbins_phi_NotCorr[m][z][v][tr]->GetBinError(i),2);
	    sigmaAC=sqrt(sigmaAC);///1
	    sigmaAC=(sigmaAC/hDeltaEtaDeltaPhi_ACbins_phi_NotCorr[m][z][v][tr]->GetBinContent(i) + HistoTriggerEfficiency->GetBinContent(m+1)/HistoTriggerEfficiency->GetBinError(m+1) + fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)/fHistV0EfficiencyPtBins[m]->GetBinError(v+1) + HistContTriggerMolt->GetBinContent(m+1)/HistContTriggerMolt->GetBinError(m+1) + HistContV0PtBins[m]->GetBinContent(v+1)/HistContV0PtBins[m]->GetBinError(v+1))*hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr]->GetBinContent(i) ;//2
	    //cout << " \n\n contribution from different terms in the sum: \n  ACbins: " <<  sigmaAC/hDeltaEtaDeltaPhi_ACbins_phi_NotCorr[m][z][v][tr]->GetBinContent(i) << "\n trigger efficiency " <<  HistoTriggerEfficiency->GetBinContent(m+1)/HistoTriggerEfficiency->GetBinError(m+1) << "\n V0 efficiency " << fHistV0EfficiencyPtBins[m]->GetBinContent(v+1)/fHistV0EfficiencyPtBins[m]->GetBinError(v+1) << "\n contamination trigger " << HistContTriggerMolt->GetBinContent(m+1)/HistContTriggerMolt->GetBinError(m+1) << "\n contaminatio V0 " << HistContV0PtBins[m]->GetBinContent(v+1)/HistContV0PtBins[m]->GetBinError(v+1);
	    
	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->SetBinContent(i, (Float_t)(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr]->GetBinContent(i))-DeltaPhiBkg);
	     hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->SetBinError(i, sqrt(pow(hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr]->GetBinError(i),2)+ pow(DeltaPhiBkgError,2))); //3
	    // cout << "\n\n valore bin phi: "<< hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->GetBinContent(i) << " +- " <<  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->GetBinError(i)<< "\n errore phi distribution senza sottrazione bkg " << sigmaAC << "\n errore bkg subtraction " << DeltaPhiBkgError << endl;
	  }
	  */

	  gauss[m][z][v][tr] = new TF1("gaus_m"+ Smolt[m]+"_v"+SPtV0[v],"gaus",-0.5*TMath::Pi(), 1.5*TMath::Pi());
	  gauss[m][z][v][tr]->SetLineColor(kRed);   
	  gauss[m][z][v][tr]->SetParameter(1,0);
	  gauss[m][z][v][tr]->SetParName(0, "norm");
	  gauss[m][z][v][tr]->SetParName(1, "mean");
	  gauss[m][z][v][tr]->SetParName(2, "sigma");
	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]-> Fit(gauss[m][z][v][tr], "", "", -1, 1); //non ho ben capito l;a differenza tra 0 e N
	  NSpectrum[m][v]=gauss[m][z][v][tr]->Integral(-0.5*TMath::Pi(), 1.5*TMath::Pi())/	  hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->GetBinWidth(1);
	  cout << NTrigger[m] << endl;
	  if (NTrigger[m]==0)	    cout << "\n \n ERROR: number of trigger particles is zero! " << endl;
	  if (NTrigger[m]!=0){ 
	    NSpectrum[m][v]= NSpectrum[m][v]/NTrigger[m];
	  }
	  cout << "integrale " << gauss[m][z][v][tr]->Integral(-0.5*TMath::Pi(), 1.5*TMath::Pi()) << " n trigger " << NTrigger[m]<< endl;
	  fHistSpectrum[m]->SetBinContent(v+1, NSpectrum[m][v]);
	}
      }
    }
  }


  cout << "\n \n sto per fare la canvas dello spettro " << endl;
  TCanvas *canvasSpectrum= new TCanvas ("canvasspectrum", "canvas");
  for(Int_t m=0; m< nummolt; m++){
    //    canvasSpectrum->cd();
    //fHistSpectrum[m]->GetYaxis()->SetRangeUser(0,60);
    fHistSpectrum[m]->SetMarkerStyle(Marker[m]);
    fHistSpectrum[m]->SetLineColor(Color[m]);
    fHistSpectrum[m]->SetMarkerColor(Color[m]);
    fHistSpectrum[m]->Draw("same");
    }


  cout << "\n sto per scrivere su file " << endl;

  TFile *fileout = new TFile(PathOut1, "RECREATE");
  for(Int_t m=0; m<nummolt+1; m++){
    if (NTrigger[m]!=0) fHistSpectrum[m]->Write();
      for(Int_t z=0; z<numzeta; z++){
	for(Int_t tr=0; tr<numPtTrigger; tr++){
	  for(Int_t v=0; v<numPtV0; v++){
	    for(Int_t sb=0; sb<numSB; sb++){
	    if (norm_MEbins[m][v][sb]==0) continue;
	    if(!isProjectionPhi[m][z][v][tr]) continue;
	    hDeltaEtaDeltaPhi_ME_normbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins[m][z][v][tr][sb]->Write();
	    hDeltaEtaDeltaPhi_ACbins_phi[m][z][v][tr][sb]->Write();
	    if(sb==0)	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub_NotCorr[m][z][v][tr]->Write();
	    if(sb==0)	    hDeltaEtaDeltaPhi_ACbins_phi_V0Sub[m][z][v][tr]->Write();
	    if(sb==0)	    hDeltaEtaDeltaPhi_ACbins_phi_BkgSub[m][z][v][tr]->Write();
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
