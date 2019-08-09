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

Double_t SetEfficiencyError(Int_t k, Int_t n){
  return ((Float_t)k+1)*((Float_t)k+2)/(n+2)/(n+3) - pow((Float_t)(k+1),2)/pow(n+2,2);
}

void Efficiency(Int_t type=0,   Float_t ptjmin=3,  Float_t ptjmax=30, Int_t sysTrigger=0, Int_t sysV0=0){
  //TString file = "AngularCorrelation2018d8_allrunswMult";
  //TString file = "AngularCorrelation2018d8_MC_runNew";
  //TString file = "AngularCorrelation2018d8_MC_try7_10runs";
  //TString file = "AngularCorrelation2017d20a2_extra_MCEffEvt_3runs";
  // TString file = "AngularCorrelation2018d8_MC_try2fullrun";
  TString file = "AngularCorrelation2018d8_MCEff_15runs_6thtry";
  // TString PathInBis =  "AnalysisResults2018d8_allrunswMult_MC.root";
  // TString PathInBis =  "AnalysisResults2018d8_MC_runNew.root";
  // TString PathInBis =  "AnalysisResults2017d20a2_extra_MCEffEvt_3runs.root";
 // TString PathInBis =  "AnalysisResults2018d8_MC_try2fullrun.root";
  TString PathInBis =  "AnalysisResults2018d8_MCEff_15runs_6thtry.root";
  //  TString PathInBis =  "AnalysisResults.root";
  TString PathOut2="histo/" + file + Form("_MC_Efficiency_SysT%i_SysV0%i.root", sysTrigger, sysV0);

  TFile *fileinbis = new TFile(PathInBis);
  TDirectoryFile *dir = (TDirectoryFile*)fileinbis->Get("MyTask");
  TList *list = (TList*)dir->Get("MyOutputContainer");
  TList *list2 = (TList*)dir->Get("MyOutputContainer3");

  const Int_t nummolt=5;
  const Int_t numzeta=1;
  const Int_t numPtV0=5;
  const Int_t numPtTrigger=1;
  const Int_t numtipo=2;
  const Int_t numSelTrigger=3;
  const Int_t numSelV0=6;
  const Int_t numSysTrigger = 3;
  const Int_t numSysV0 = 7;
  Int_t sysV0Gen=0;
  if(sysV0==3) sysV0Gen=1;

  TString tipo[numtipo]={"kK0s", "bo"};
  
  TString Smolt[nummolt+1]={"0-5", "5-10", "10-30", "30-50", "50-100", "_all"};
  Double_t Nmolt[nummolt+1]={0,5,10,30,50,100}; 
  TString Szeta[numzeta]={""};
  //  Double_t Nzeta[numzeta+1]={};
  TString SPtV0[numPtV0]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  Double_t NPtV0[numPtV0+1]={0,1,2,3,4,8};
  // TString SPtTrigger[numPtTrigger]={"2-10"};
  // Double_t NPtTrigger[numPtTrigger+1]={2,10};

  Int_t Marker[nummolt+1]={7,4,20,22,29, 35};
  Int_t Color[nummolt+1]={2,3,4,6,8,9};
  Int_t ColorSysTrigger[numSysTrigger]={2,3,4};
  Int_t ColorSysV0[numSysV0]={2,3,4, 6,8,9};

  TCanvas *canvasEff=new TCanvas ("canvasEff", "canvasEff", 1500, 800);
  TCanvas *canvasRes=new TCanvas ("canvasRes", "canvasRes", 1500, 800);
  TCanvas *canvasCont=new TCanvas ("canvasCont", "canvasCont", 1000, 700);
  TCanvas *canvasUsed=new TCanvas ("canvasUsed", "canvasUsed", 1000, 700);

  canvasEff->Divide(3,2);
  canvasRes->Divide(3,2);
  canvasCont->Divide(2,1);
  canvasUsed->Divide(2,2);


  //***************istogrammi delle efficienze**********************************************************
  TH3D*   fHistGeneratedTriggerPtPhi= (TH3D*)list2->FindObject("fHistGeneratedTriggerPtPhi");
  //  TH3D*   fHistGeneratedTriggerPtPhi= (TH3D*)list2->FindObject(Form("fHistGeneratedTriggerPtPhi_%i",sysV0Gen));
  TH3D*   fHistSelectedTriggerPtPhi=  (TH3D*)list2->FindObject(Form("fHistSelectedTriggerPtPhi_%i", sysTrigger));
  TH3D*   fHistGeneratedV0PtPhi=      (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtPhi_%i", sysV0Gen));
  TH3D*   fHistSelectedV0PtPhi=       (TH3D*)list2->FindObject(Form("fHistSelectedV0PtPhi_%i",sysV0));
  TH3D*   fHistGeneratedTriggerPtEta= (TH3D*)list2->FindObject("fHistGeneratedTriggerPtEta");
  TH3D*   fHistSelectedTriggerPtEta=  (TH3D*)list2->FindObject(Form("fHistSelectedTriggerPtEta_%i", sysTrigger));
  TH3D*   fHistGeneratedV0PtEta=      (TH3D*)list2->FindObject(Form("fHistGeneratedV0PtEta_%i", sysV0Gen));
  TH3D*   fHistSelectedV0PtEta=       (TH3D*)list2->FindObject(Form("fHistSelectedV0PtEta_%i", sysV0));
  TH3D*   fHistReconstructedV0PtMass= (TH3D*)list->FindObject("fHistReconstructedV0PtMass");
  TH3D*   fHistSelectedV0PtMass=      (TH3D*)list->FindObject("fHistSelectedV0PtMass");

  //***************istogrammi delle risoluzioni**********************************************************
  TH2D*   fHistResolution[3][2];
  TString nameRes[3][2];
  TH1D*   fHistResolution_1D[nummolt+1][3][2];
  TString nameRes_1D[nummolt+1][3][2]; 
  TString TorV[2]={"Trigger", "V0"};
  TString Var[3]={"Pt", "Phi", "Eta"};

  for(Int_t m=0; m< 3; m++){
    for(Int_t t=0; t< 2; t++){
      nameRes[m][t]="fHistResolution" + TorV[t] + Var[m];
      fHistResolution[m][t]=   (TH2D*)list->FindObject(nameRes[m][t]);
    }
  }

  auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->SetHeader("Intervalli di molteplicita'");     

  //distribuzioni  Pt, Phi, Eta in 2D delle selezionate e generate
  TH2D*  fHistSelected_2D[nummolt+1];
  TH2D*  fHistGenerated_2D[nummolt+1];
  
  TH2D*  fHistSelected_2D_TriggerPtPhi[nummolt+1];
  TH2D*  fHistGenerated_2D_TriggerPtPhi[nummolt+1];
  TH2D*  fHistSelected_2D_TriggerPtEta[nummolt+1];
  TH2D*  fHistGenerated_2D_TriggerPtEta[nummolt+1];

  TH2D*  fHistSelected_2D_V0PtPhi[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtPhi[nummolt+1];
  TH2D*  fHistSelected_2D_V0PtEta[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtEta[nummolt+1];
  TH2D*  fHistSelected_2D_V0PtPhi_clone[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtPhi_clone[nummolt+1];
  TH2D*  fHistSelected_2D_V0PtEta_clone[nummolt+1];
  TH2D*  fHistGenerated_2D_V0PtEta_clone[nummolt+1];

  TH2D*  fHistSelectedMass_2D[nummolt+1];
  TH2D*  fHistRecoMass_2D[nummolt+1];

  //distribuzioni Pt, Phi, Eta in 1D delle selezionate e generate
  TH1D*  fHistSelected_1D_TriggerPt[nummolt+1];
  TH1D*  fHistGenerated_1D_TriggerPt[nummolt+1];
  TH1D*  fHistSelected_1D_V0Pt[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Pt[nummolt+1];
  TH1D*  fHistSelected_1D_TriggerPhi[nummolt+1];
  TH1D*  fHistGenerated_1D_TriggerPhi[nummolt+1];
  TH1D*  fHistSelected_1D_V0Phi[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Phi[nummolt+1];
  TH1D*  fHistSelected_1D_TriggerEta[nummolt+1];
  TH1D*  fHistGenerated_1D_TriggerEta[nummolt+1];
  TH1D*  fHistSelected_1D_V0Eta[nummolt+1];
  TH1D*  fHistGenerated_1D_V0Eta[nummolt+1];
  TH1D*  fHistSelectedMass_1D[nummolt+1];
  TH1D*  fHistRecoMass_1D[nummolt+1];
  
  TH2D*  fHistTriggerEfficiencyPtPhi[nummolt+1];                                 
  TH2D*  fHistTriggerEfficiencyPtEta[nummolt+1]; 
  TH2D*  fHistV0EfficiencyPtPhi[nummolt+1];
  TH2D*  fHistV0EfficiencyPtEta[nummolt+1];
  TH2D*  fHistEfficiencyV0Selection[nummolt+1];

  TH1D*  fHistTriggerEfficiencyPt[nummolt+1];
  TH1D*  fHistTriggerEfficiencyPhi[nummolt+1];
  TH1D*  fHistTriggerEfficiencyEta[nummolt+1];
  TH1D*  fHistV0EfficiencyPt[nummolt+1];
  TH1D*  fHistV0EfficiencyPhi[nummolt+1];
  TH1D*  fHistV0EfficiencyEta[nummolt+1];
  
  TH1D*  fHistV0EfficiencyPtBins[nummolt+1];
  TH1D*  HistoTriggerEfficiency= new TH1D("HistoTriggerEfficiency", "Trigger selection efficiency vs centrality", nummolt,Nmolt );

  //  TH1D*  fHistV0MeanEfficiencyPt[nummolt+1];

  TH2D* fHistV0EfficiencyReco[nummolt+1];
  TH1D* fHistV0EfficiencyRecoPt[nummolt+1];
  // TH2D* fHistV0EfficiencyRecoMolt;
  // TH1D* fHistV0EfficiencyRecoPtMolt;
  
  //TH1D* fHistV0MeanEfficiencyRecoPt[nummolt+1];

  //istogrammi per contaminazioni
  TH2D * HistContTrigger[nummolt+1];
  TH2D * HistContV0[nummolt+1];
  TH1D * HistContTriggerPt[nummolt+1];
  TH1D * HistContV0Pt[nummolt+1];
  TH1D * HistContV0PtBins[nummolt+1];
  TH1D * HistInt[nummolt+1];
  TH1D * HistContTriggerMolt = new TH1D("HistContTriggerMolt", "Trigger contamination factor vs multiplicity", nummolt, Nmolt);

  Float_t ContTrigger[nummolt+1];
  Float_t ContV0[nummolt+1];
  Float_t ContTriggerInt;
  Float_t ContV0Int;
  Float_t ContTriggerIntError;
  Float_t ContV0IntError;

  Float_t TriggerEfficiency[nummolt+1];

  cout << "*************************************************************************************" << endl;
  cout << "Sto effettuando istogrammi efficienze in range di molteplicita'" << endl;
    
  for(Int_t molt=0; molt < nummolt+1; molt++){
    cout << "\n\n I'm analyzing multiplcity interval n. " << molt<< endl;
  
    cout << "Trigger 2D projection in Phi and Pt " << endl;
    if(molt < nummolt){
    fHistSelectedTriggerPtPhi->GetZaxis()->SetRange(fHistSelectedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    fHistGeneratedTriggerPtPhi->GetZaxis()->SetRange(fHistGeneratedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedTriggerPtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    fHistSelected_2D_TriggerPtPhi[molt] = (TH2D*)fHistSelectedTriggerPtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_TriggerPtPhi[molt] = (TH2D*)fHistGeneratedTriggerPtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axi
    fHistSelected_2D_TriggerPtPhi[molt] ->SetName("fHistSelected_2D_TriggerPtPhi_"+ Smolt[molt] );
    fHistGenerated_2D_TriggerPtPhi[molt]->SetName("fHistGenerated_2D_TriggerPtPhi_"+ Smolt[molt]);
    fHistTriggerEfficiencyPtPhi[molt]= new TH2D("fHistTriggerEfficiencyPtPhi_"+ Smolt[molt],"fHistTriggerEfficiencyPtPhi_"+ Smolt[molt],fHistSelectedTriggerPtPhi->GetNbinsX(),fHistSelectedTriggerPtPhi->GetXaxis()->GetXmin(), fHistSelectedTriggerPtPhi->GetXaxis()->GetXmax(),fHistSelectedTriggerPtPhi->GetNbinsY(),fHistSelectedTriggerPtPhi->GetYaxis()->GetBinLowEdge(1), fHistSelectedTriggerPtPhi->GetYaxis()->GetBinUpEdge(fHistSelectedTriggerPtPhi->GetNbinsY()) );
    fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->SetTitle("p_{T}");      
    fHistTriggerEfficiencyPtPhi[molt]->GetYaxis()->SetTitle("#phi");   

    fHistSelected_2D_TriggerPtPhi[molt] ->RebinX(1);   
    fHistSelected_2D_TriggerPtPhi[molt] ->RebinY(5);   
    
    fHistGenerated_2D_TriggerPtPhi[molt]->RebinX(1);    
    fHistGenerated_2D_TriggerPtPhi[molt]->RebinY(5);    
    
    fHistTriggerEfficiencyPtPhi[molt]->RebinX(1);    
    fHistTriggerEfficiencyPtPhi[molt]->RebinY(5);    
   
    fHistTriggerEfficiencyPtPhi[molt]->Divide(fHistSelected_2D_TriggerPtPhi[molt], fHistGenerated_2D_TriggerPtPhi[molt]); 
  
    cout << "Trigger 1D projection in Phi and Pt " << endl;
    fHistTriggerEfficiencyPt[molt]= new TH1D("fHistTriggerEfficiencyPt_"+ Smolt[molt] , "fHistTriggerEfficiencyPt_"+ Smolt[molt] ,  fHistTriggerEfficiencyPtPhi[molt]->GetNbinsX(), fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->GetXmin(), fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->GetXmax() );
    fHistTriggerEfficiencyPt[molt]->GetXaxis()->SetTitle("p_{T}");      
    fHistTriggerEfficiencyPt[molt]->SetTitle("fHistTriggerEfficiencyPt_"+ Smolt[molt]);
    fHistSelected_1D_TriggerPt[molt]=(TH1D*)fHistSelected_2D_TriggerPtPhi[molt]->ProjectionX("fHistSelected_1D_TriggerPt_"+ Smolt[molt]) ;
    fHistGenerated_1D_TriggerPt[molt]=(TH1D*)fHistGenerated_2D_TriggerPtPhi[molt]->ProjectionX("fHistGenerated_1D_TriggerPt_"+ Smolt[molt]) ;
    fHistTriggerEfficiencyPt[molt] ->Divide(  fHistSelected_1D_TriggerPt[molt],   fHistGenerated_1D_TriggerPt[molt] );


    
    canvasEff->cd(1);
    fHistTriggerEfficiencyPt[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyPt[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyPt[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyPt[molt]->SetMarkerColor(Color[molt]);
    legend->AddEntry(fHistTriggerEfficiencyPt[molt],Smolt[molt],"pel");   

    fHistTriggerEfficiencyPt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();

    fHistTriggerEfficiencyPhi[molt]= new TH1D("fHistTriggerEfficiencyPhi_"+ Smolt[molt] , "fHistTriggerEfficiencyPhi_"+ Smolt[molt] ,  fHistTriggerEfficiencyPtPhi[molt]->GetNbinsY(), fHistTriggerEfficiencyPtPhi[molt]->GetYaxis()->GetBinLowEdge(1), fHistTriggerEfficiencyPtPhi[molt]->GetYaxis()->GetBinUpEdge(fHistTriggerEfficiencyPtPhi[molt]->GetNbinsY()) );
    fHistTriggerEfficiencyPhi[molt]->GetXaxis()->SetTitle("#phi");      
    fHistTriggerEfficiencyPhi[molt]->SetTitle("fHistTriggerEfficiencyPhi_"+ Smolt[molt] );  
    fHistSelected_1D_TriggerPhi[molt]=(TH1D*)fHistSelected_2D_TriggerPtPhi[molt]->ProjectionY("fHistSelected_1D_TriggerPhi_"+ Smolt[molt],fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->FindBin(ptjmax -0.0001) );
    fHistGenerated_1D_TriggerPhi[molt]=(TH1D*)fHistGenerated_2D_TriggerPtPhi[molt]->ProjectionY("fHistGenerated_1D_TriggerPhi_"+ Smolt[molt],fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistTriggerEfficiencyPtPhi[molt]->GetXaxis()->FindBin(ptjmax -0.0001));
    fHistTriggerEfficiencyPhi[molt] ->Divide (  fHistSelected_1D_TriggerPhi[molt], fHistGenerated_1D_TriggerPhi[molt]);

    canvasEff->cd(2);
    fHistTriggerEfficiencyPhi[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyPhi[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyPhi[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyPhi[molt]->SetMarkerColor(Color[molt]);
    fHistTriggerEfficiencyPhi[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();
  
    //Trigger efficiency integrated in all variables except multiplicity (ptmin < pt < ptmax)
    Double_t SelEntries=0;
    Double_t GenEntries=0;
    for(Int_t j=fHistSelected_1D_TriggerPt[molt]->FindBin(ptjmin); j<fHistSelected_1D_TriggerPt[molt]->FindBin(ptjmax); j++ ){
      SelEntries+=fHistSelected_1D_TriggerPt[molt]->GetBinContent(j);
    }
    for(Int_t j=fHistGenerated_1D_TriggerPt[molt]->FindBin(ptjmin); j<fHistGenerated_1D_TriggerPt[molt]->FindBin(ptjmax); j++ ){
      GenEntries+=fHistGenerated_1D_TriggerPt[molt]->GetBinContent(j);
    }
    TriggerEfficiency[molt]=(SelEntries/GenEntries);
    HistoTriggerEfficiency->SetBinContent(molt+1, TriggerEfficiency[molt]);
    HistoTriggerEfficiency->SetBinError(molt+1, SetEfficiencyError(SelEntries, GenEntries));
    canvasUsed->cd(1);
    HistoTriggerEfficiency->GetYaxis()->SetRangeUser(0,0.15);
    HistoTriggerEfficiency->SetMarkerStyle(ColorSysTrigger[sysTrigger]);
    HistoTriggerEfficiency->SetLineColor(ColorSysTrigger[sysTrigger]);
    HistoTriggerEfficiency->SetMarkerColor(ColorSysTrigger[sysTrigger]);
    HistoTriggerEfficiency->Draw("e");
      
    //per la particella di trigger - PtEta
    if(molt < nummolt){
    fHistSelectedTriggerPtEta->GetZaxis()->SetRange(fHistSelectedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    fHistGeneratedTriggerPtEta->GetZaxis()->SetRange(fHistGeneratedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedTriggerPtEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    fHistSelected_2D_TriggerPtEta[molt] = (TH2D*)fHistSelectedTriggerPtEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_TriggerPtEta[molt] = (TH2D*)fHistGeneratedTriggerPtEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axi
    fHistSelected_2D_TriggerPtEta[molt] ->SetName("fHistSelected_2D_TriggerPtEta_"+ Smolt[molt] );
    fHistGenerated_2D_TriggerPtEta[molt]->SetName("fHistGenerated_2D_TriggerPtEta_"+ Smolt[molt]);
    fHistTriggerEfficiencyPtEta[molt]= new TH2D("fHistTriggerEfficiencyPtEta_"+ Smolt[molt],"fHistTriggerEfficiencyPtEta_"+ Smolt[molt],fHistSelectedTriggerPtEta->GetNbinsX(),fHistSelectedTriggerPtEta->GetXaxis()->GetXmin(), fHistSelectedTriggerPtEta->GetXaxis()->GetXmax(),fHistSelectedTriggerPtEta->GetNbinsY(),fHistSelectedTriggerPtEta->GetYaxis()->GetBinLowEdge(1), fHistSelectedTriggerPtEta->GetYaxis()->GetBinUpEdge(fHistSelectedTriggerPtEta->GetNbinsY()) );
    fHistTriggerEfficiencyPtEta[molt]->GetXaxis()->SetTitle("p_{T}");      
    fHistTriggerEfficiencyPtEta[molt]->GetYaxis()->SetTitle("#phi");      
    
    fHistSelected_2D_TriggerPtEta[molt] ->RebinX(1);   
    fHistSelected_2D_TriggerPtEta[molt] ->RebinY(10);   
    
    fHistGenerated_2D_TriggerPtEta[molt]->RebinX(1);    
    fHistGenerated_2D_TriggerPtEta[molt]->RebinY(10);    
    
    fHistTriggerEfficiencyPtEta[molt]->RebinX(1);    
    fHistTriggerEfficiencyPtEta[molt]->RebinY(10);    
    
    fHistTriggerEfficiencyPtEta[molt]->Divide(fHistSelected_2D_TriggerPtEta[molt], fHistGenerated_2D_TriggerPtEta[molt]); 
    
    fHistTriggerEfficiencyEta[molt]= new TH1D("fHistTriggerEfficiencyEta_"+ Smolt[molt] , "fHistTriggerEfficiencyEta_"+ Smolt[molt] ,  fHistTriggerEfficiencyPtEta[molt]->GetNbinsY(), fHistTriggerEfficiencyPtEta[molt]->GetYaxis()->GetBinLowEdge(1), fHistTriggerEfficiencyPtEta[molt]->GetYaxis()->GetBinUpEdge(fHistTriggerEfficiencyPtEta[molt]->GetNbinsY()) );
    fHistTriggerEfficiencyEta[molt]->GetXaxis()->SetTitle("#eta");      
    fHistTriggerEfficiencyEta[molt]->SetTitle("fHistTriggerEfficiencyEta_"+ Smolt[molt] );  
    fHistSelected_1D_TriggerEta[molt]=(TH1D*)fHistSelected_2D_TriggerPtEta[molt]->ProjectionY("fHistSelected_1D_TriggerEta_"+ Smolt[molt],fHistTriggerEfficiencyPtEta[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistTriggerEfficiencyPtEta[molt]->GetXaxis()->FindBin(ptjmax -0.0001));
    fHistGenerated_1D_TriggerEta[molt]=(TH1D*)fHistGenerated_2D_TriggerPtEta[molt]->ProjectionY("fHistGenerated_1D_TriggerEta_"+ Smolt[molt],fHistTriggerEfficiencyPtEta[molt]->GetXaxis()->FindBin(ptjmin +0.0001), fHistTriggerEfficiencyPtEta[molt]->GetXaxis()->FindBin(ptjmax -0.0001));
    
    fHistTriggerEfficiencyEta[molt] ->Divide (fHistSelected_1D_TriggerEta[molt], fHistGenerated_1D_TriggerEta[molt]);
    
    canvasEff->cd(3);
    fHistTriggerEfficiencyEta[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistTriggerEfficiencyEta[molt]->SetMarkerStyle(Marker[molt]);
    fHistTriggerEfficiencyEta[molt]->SetLineColor(Color[molt]);
    fHistTriggerEfficiencyEta[molt]->SetMarkerColor(Color[molt]);
    fHistTriggerEfficiencyEta[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();
  
    

    /* per calcolo errore efficienza
    for(Int_t i=1 ; i< fHistTriggerEfficiencyEta[molt]->GetNbinsX(); i++){
      //cfHistTriggerEfficiencyEta[molt]->SetBinContent(i,(Float_t)((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i)/((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i));
      if( ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i) != 0){
	//ccout <<  ((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i) << " " <<  ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i) << " " <<SetEfficiencyError(((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i), ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i))<< endl;
	fHistTriggerEfficiencyEta[molt]->SetBinError(i, SetEfficiencyError(((TH1D*)fHistSelected_2D->ProjectionY())->GetBinContent(i), ((TH1D*)fHistGenerated_2D->ProjectionY())->GetBinContent(i)));
      }
    }
    */

    cout << "V0 2D projection in Phi and Pt " << endl;
    if(molt < nummolt){
    fHistSelectedV0PtPhi->GetZaxis()->SetRange(fHistSelectedV0PtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    fHistGeneratedV0PtPhi->GetZaxis()->SetRange(fHistGeneratedV0PtPhi->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedV0PtPhi->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
    fHistSelected_2D_V0PtPhi[molt] = (TH2D*)fHistSelectedV0PtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
    fHistGenerated_2D_V0PtPhi[molt] = (TH2D*)fHistGeneratedV0PtPhi->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
  fHistSelected_2D_V0PtPhi[molt] ->SetName("fHistSelected_2D_V0PtPhi_"+ Smolt[molt]);
  fHistGenerated_2D_V0PtPhi[molt]->SetName("fHistGenerated_2D_V0PtPhi_"+ Smolt[molt]);
  fHistSelected_2D_V0PtPhi_clone[molt]  =  (TH2D*)fHistSelected_2D_V0PtPhi[molt]->Clone("fHistSelected_2D_V0PtPhi_clone_"+ Smolt[molt]);
  fHistGenerated_2D_V0PtPhi_clone[molt] =  (TH2D*)fHistGenerated_2D_V0PtPhi[molt]->Clone("fHistGenerated_2D_V0PtPhi_clone_"+ Smolt[molt]);

  fHistV0EfficiencyPtPhi[molt]= new TH2D("fHistV0EfficiencyPtPhi_"+ Smolt[molt],"fHistV0EfficiencyPtPhi_"+ Smolt[molt],fHistSelectedV0PtPhi->GetNbinsX(),fHistSelectedV0PtPhi->GetXaxis()->GetXmin(), fHistSelectedV0PtPhi->GetXaxis()->GetXmax(),fHistSelectedV0PtPhi->GetNbinsY(),fHistSelectedV0PtPhi->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtPhi->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtPhi->GetNbinsY()) );
  fHistV0EfficiencyPtPhi[molt]->GetXaxis()->SetTitle("p_{T}");      
  fHistV0EfficiencyPtPhi[molt]->GetYaxis()->SetTitle("#phi");

  fHistSelected_2D_V0PtPhi[molt] ->RebinX(1);   
  fHistSelected_2D_V0PtPhi[molt] ->RebinY(5);   

  fHistGenerated_2D_V0PtPhi[molt]->RebinX(1);    
  fHistGenerated_2D_V0PtPhi[molt]->RebinY(5);    

  fHistV0EfficiencyPtPhi[molt]->RebinX(1);    
  fHistV0EfficiencyPtPhi[molt]->RebinY(5);    

  fHistV0EfficiencyPtPhi[molt]->Divide(fHistSelected_2D_V0PtPhi[molt], fHistGenerated_2D_V0PtPhi[molt]); 

  cout << "V0 1D projection in Phi and Pt " << endl;
  fHistV0EfficiencyPt[molt]= new TH1D("fHistV0EfficiencyPt_"+ Smolt[molt] , "fHistV0EfficiencyPt_"+ Smolt[molt] ,  fHistSelected_2D_V0PtPhi[molt]->GetNbinsX(), fHistSelected_2D_V0PtPhi[molt]->GetXaxis()->GetXmin(), fHistSelected_2D_V0PtPhi[molt]->GetXaxis()->GetXmax() );
  fHistV0EfficiencyPt[molt]->GetXaxis()->SetTitle("p_{T}");      
  fHistV0EfficiencyPt[molt]->SetTitle("fHistV0EfficiencyPt_"+ Smolt[molt] );

  fHistSelected_1D_V0Pt[molt]=(TH1D*)fHistSelected_2D_V0PtPhi[molt]->ProjectionX("fHistSelected_1D_V0Pt_"+ Smolt[molt]);
  fHistGenerated_1D_V0Pt[molt]=(TH1D*)fHistGenerated_2D_V0PtPhi[molt]->ProjectionX("fHistGenerated_1D_V0Pt_"+ Smolt[molt]); 
  fHistV0EfficiencyPt[molt]->Divide(fHistSelected_1D_V0Pt[molt],  fHistGenerated_1D_V0Pt[molt]);

  canvasEff->cd(4);
  fHistV0EfficiencyPt[molt]->GetYaxis()->SetRangeUser(0,1);
  fHistV0EfficiencyPt[molt]->SetMarkerStyle(Marker[molt]);
  fHistV0EfficiencyPt[molt]->SetLineColor(Color[molt]);
  fHistV0EfficiencyPt[molt]->SetMarkerColor(Color[molt]);
  fHistV0EfficiencyPt[molt]->Draw("same");
  if(molt ==nummolt)     legend->Draw();
  
  cout << "V0 efficiency in Pt bins used in the analysis " << endl;
  fHistV0EfficiencyPtBins[molt]= new TH1D("fHistV0EfficiencyPtBins_" + Smolt[molt], "fHistV0EfficiencyPtBins_" + Smolt[molt],numPtV0, NPtV0 );
  fHistV0EfficiencyPtBins[molt]->GetXaxis()->SetTitle("p_{T}");      
  Float_t NumberOfSelected;
  Float_t NumberOfGenerated;
  cout << "\n V0 efficiency in Pt bins " << endl;
  for(Int_t j=0; j<numPtV0; j++){
    NumberOfSelected=0;
    NumberOfGenerated=0;
    for(Int_t i=fHistSelected_1D_V0Pt[molt]->GetXaxis()->FindBin(NPtV0[j]+0.00001); i<=fHistSelected_1D_V0Pt[molt]->GetXaxis()->FindBin(NPtV0[j+1]-0.00001); i++){
      NumberOfGenerated+=fHistGenerated_1D_V0Pt[molt]->GetBinContent(i);
      NumberOfSelected+=fHistSelected_1D_V0Pt[molt]->GetBinContent(i);
    }
    fHistV0EfficiencyPtBins[molt]->SetBinContent(j+1, NumberOfSelected/NumberOfGenerated);
    fHistV0EfficiencyPtBins[molt]->SetBinError(j+1, SetEfficiencyError(NumberOfSelected,NumberOfGenerated));
    cout << 	fHistV0EfficiencyPtBins[molt]->GetBinContent(j+1) << "+-" << 	fHistV0EfficiencyPtBins[molt]->GetBinError(j+1)<< endl;
  }
    canvasUsed->cd(2);
    fHistV0EfficiencyPtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    fHistV0EfficiencyPtBins[molt]->SetMarkerStyle(Marker[molt]);
    fHistV0EfficiencyPtBins[molt]->SetLineColor(Color[molt]);
    fHistV0EfficiencyPtBins[molt]->SetMarkerColor(Color[molt]);
    fHistV0EfficiencyPtBins[molt]->Draw("samee");
    if(molt ==nummolt)     legend->Draw();


  fHistV0EfficiencyPhi[molt]= new TH1D("fHistV0EfficiencyPhi_"+ Smolt[molt] , "fHistV0EfficiencyPhi_"+ Smolt[molt] , fHistSelected_2D_V0PtPhi[molt]->GetNbinsY(), fHistSelected_2D_V0PtPhi[molt]->GetYaxis()->GetBinLowEdge(1), fHistSelected_2D_V0PtPhi[molt]->GetYaxis()->GetBinUpEdge(fHistSelected_2D_V0PtPhi[molt]->GetNbinsY()) );
  fHistV0EfficiencyPhi[molt]->GetXaxis()->SetTitle("#phi");      
  fHistV0EfficiencyPhi[molt]->SetTitle("fHistV0EfficiencyPhi_"+ Smolt[molt] );

  fHistSelected_1D_V0Phi[molt]=(TH1D*)fHistSelected_2D_V0PtPhi[molt]->ProjectionY("fHistSelected_1D_V0Phi_"+ Smolt[molt]);
  fHistGenerated_1D_V0Phi[molt]=(TH1D*)fHistGenerated_2D_V0PtPhi[molt]->ProjectionY("fHistGenerated_1D_V0Phi_"+ Smolt[molt]);  
  fHistV0EfficiencyPhi[molt] ->Divide(  fHistSelected_1D_V0Phi[molt],  fHistGenerated_1D_V0Phi[molt]);
  
  canvasEff->cd(5);
  fHistV0EfficiencyPhi[molt]->GetYaxis()->SetRangeUser(0,1);
  fHistV0EfficiencyPhi[molt]->SetMarkerStyle(Marker[molt]);
  fHistV0EfficiencyPhi[molt]->SetLineColor(Color[molt]);
  fHistV0EfficiencyPhi[molt]->SetMarkerColor(Color[molt]);
  fHistV0EfficiencyPhi[molt]->Draw("same");
  if(molt ==nummolt)     legend->Draw();
  
  cout << "V0 2D projection in Eta and Pt " << endl;
    if(molt < nummolt){
  fHistSelectedV0PtEta->GetZaxis()->SetRange(fHistSelectedV0PtEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
  fHistGeneratedV0PtEta->GetZaxis()->SetRange(fHistGeneratedV0PtEta->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistGeneratedV0PtEta->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
  fHistSelected_2D_V0PtEta[molt]  = (TH2D*)fHistSelectedV0PtEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
  fHistGenerated_2D_V0PtEta[molt] = (TH2D*)fHistGeneratedV0PtEta->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
  fHistSelected_2D_V0PtEta[molt] ->SetName("fHistSelected_2D_V0PtEta_"+ Smolt[molt]);
  fHistGenerated_2D_V0PtEta[molt]->SetName("fHistGenerated_2D_V0PtEta_"+ Smolt[molt]);
  fHistSelected_2D_V0PtEta_clone[molt]  =  (TH2D*)fHistSelected_2D_V0PtEta[molt]->Clone("fHistSelected_2D_V0PtEta_clone_"+ Smolt[molt]);
  fHistGenerated_2D_V0PtEta_clone[molt] =  (TH2D*)fHistGenerated_2D_V0PtEta[molt]->Clone("fHistGenerated_2D_V0PtEta_clone_"+ Smolt[molt]);
  fHistV0EfficiencyPtEta[molt]= new TH2D("fHistV0EfficiencyPtEta_"+ Smolt[molt],"fHistV0EfficiencyPtEta_"+ Smolt[molt],fHistSelectedV0PtEta->GetNbinsX(),fHistSelectedV0PtEta->GetXaxis()->GetXmin(), fHistSelectedV0PtEta->GetXaxis()->GetXmax(),fHistSelectedV0PtEta->GetNbinsY(),fHistSelectedV0PtEta->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtEta->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtEta->GetNbinsY()) );
  fHistV0EfficiencyPtEta[molt]->GetXaxis()->SetTitle("p_{T}");      
  fHistV0EfficiencyPtEta[molt]->GetYaxis()->SetTitle("#eta");      

  fHistSelected_2D_V0PtEta[molt] ->RebinX(1);   
  fHistSelected_2D_V0PtEta[molt] ->RebinY(10);   

  fHistGenerated_2D_V0PtEta[molt]->RebinX(1);    
  fHistGenerated_2D_V0PtEta[molt]->RebinY(10);    

  fHistV0EfficiencyPtEta[molt]->RebinX(1);    
  fHistV0EfficiencyPtEta[molt]->RebinY(10);    

  fHistV0EfficiencyPtEta[molt]->Divide(fHistSelected_2D_V0PtEta[molt], fHistGenerated_2D_V0PtEta[molt]); 

  cout << "V0 1D projection in Eta " << endl;
  fHistV0EfficiencyEta[molt]= new TH1D("fHistV0EfficiencyEta_"+ Smolt[molt] , "fHistV0EfficiencyEta_"+ Smolt[molt] ,    fHistGenerated_2D_V0PtEta[molt]->GetNbinsY(),   fHistGenerated_2D_V0PtEta[molt]->GetYaxis()->GetBinLowEdge(1),   fHistGenerated_2D_V0PtEta[molt]->GetYaxis()->GetBinUpEdge(  fHistGenerated_2D_V0PtEta[molt]->GetNbinsY()) );
   fHistV0EfficiencyEta[molt]->GetXaxis()->SetTitle("#eta");      
  fHistV0EfficiencyEta[molt]->SetTitle( "fHistV0EfficiencyEta_"+ Smolt[molt] );
   fHistSelected_1D_V0Eta[molt]=(TH1D*)fHistSelected_2D_V0PtEta[molt]->ProjectionY("fHistSelected_1D_V0Eta_"+ Smolt[molt]);
  fHistGenerated_1D_V0Eta[molt]=(TH1D*)fHistGenerated_2D_V0PtEta[molt]->ProjectionY("fHistGenerated_1D_V0Eta_"+ Smolt[molt]); 
  fHistV0EfficiencyEta[molt] ->Divide(  fHistSelected_1D_V0Eta[molt],  fHistGenerated_1D_V0Eta[molt]);
   canvasEff->cd(6);
  fHistV0EfficiencyEta[molt]->GetYaxis()->SetRangeUser(0,1);
  fHistV0EfficiencyEta[molt]->SetMarkerStyle(Marker[molt]);
  fHistV0EfficiencyEta[molt]->SetLineColor(Color[molt]);
  fHistV0EfficiencyEta[molt]->SetMarkerColor(Color[molt]);
  fHistV0EfficiencyEta[molt]->Draw("same");
  if(molt ==nummolt)     legend->Draw();
  
 
  //per efficienza selezioni V0
    if(molt < nummolt){
  fHistSelectedV0PtMass->GetZaxis()->SetRange(fHistSelectedV0PtMass->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistSelectedV0PtMass->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
 fHistReconstructedV0PtMass->GetZaxis()->SetRange(fHistReconstructedV0PtMass->GetZaxis()->FindBin(Nmolt[molt]+0.001),fHistReconstructedV0PtMass->GetZaxis()->FindBin(Nmolt[molt+1]-0.001) );
    }
 
  fHistSelectedMass_2D[molt] = (TH2D*)fHistSelectedV0PtMass->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
  fHistRecoMass_2D[molt] = (TH2D*)fHistReconstructedV0PtMass->Project3D("yxo"); //y is on the vertical axis, x along horizontal axis
  fHistV0EfficiencyReco[molt]= new TH2D("fHistV0EfficiencyReco_"+ Smolt[molt],"fHistV0EfficiencyReco_"+ Smolt[molt],fHistSelectedV0PtMass->GetNbinsX(),fHistSelectedV0PtMass->GetXaxis()->GetXmin(), fHistSelectedV0PtMass->GetXaxis()->GetXmax(),fHistSelectedV0PtMass->GetNbinsY(),fHistSelectedV0PtMass->GetYaxis()->GetBinLowEdge(1), fHistSelectedV0PtMass->GetYaxis()->GetBinUpEdge(fHistSelectedV0PtMass->GetNbinsY()) );
  fHistV0EfficiencyReco[molt]->GetXaxis()->SetTitle("Mass");      
  fHistV0EfficiencyReco[molt]->GetYaxis()->SetTitle("p_{T}");      
  fHistV0EfficiencyReco[molt]->Divide(fHistSelectedMass_2D[molt], fHistRecoMass_2D[molt]);

 
  fHistV0EfficiencyRecoPt[molt]= new TH1D("fHistV0EfficiencyRecoPt_"+ Smolt[molt] , "fHistV0EfficiencyRecoPt_"+ Smolt[molt],  fHistV0EfficiencyReco[molt]->GetNbinsY(), fHistV0EfficiencyReco[molt]->GetYaxis()->GetBinLowEdge(1), fHistV0EfficiencyReco[molt]->GetYaxis()->GetBinUpEdge(fHistV0EfficiencyReco[molt]->GetNbinsY()));
  fHistSelectedMass_1D[molt]= (TH1D*)fHistSelectedMass_2D[molt]->ProjectionY("fHistSelectedMass_1D_"+ Smolt[molt]);
  fHistRecoMass_1D[molt]= (TH1D*)fHistRecoMass_2D[molt]->ProjectionY("fHistRecoMass_1D_"+ Smolt[molt]);  
  fHistV0EfficiencyRecoPt[molt] ->Divide(fHistSelectedMass_1D[molt], fHistRecoMass_1D[molt]);

  cout << "resolution histos" << endl;
    //**********************resolution histograms*****************************
    for(Int_t m=0; m< 3; m++){
    for(Int_t t=0; t< 2; t++){
      nameRes_1D[molt][m][t]="fHistResolution" + TorV[t] + Var[m] + Smolt[molt];
      if(molt< nummolt){
      fHistResolution_1D[molt][m][t]=(TH1D*)fHistResolution[m][t]->ProjectionX(nameRes_1D[molt][m][t], fHistResolution[m][t]->GetYaxis()->FindBin(Nmolt[molt]+0.0001), fHistResolution[m][t]->GetYaxis()->FindBin(Nmolt[molt+1]-0.0001));
      }
      else {
      fHistResolution_1D[molt][m][t]=(TH1D*)fHistResolution[m][t]->ProjectionX(nameRes_1D[molt][m][t]);
      }
      Int_t normal = fHistResolution_1D[molt][m][t]->GetEntries();
      fHistResolution_1D[molt][m][t]->Scale(1./normal);
      if (t==0)      canvasRes->cd(m+1);
      else  canvasRes->cd(m+4);
      fHistResolution_1D[molt][m][t]->GetYaxis()->SetRangeUser(0, 1);
      fHistResolution_1D[molt][m][t]->GetXaxis()->SetRangeUser(-0.05, 0.05);
      fHistResolution_1D[molt][m][t]->SetMarkerStyle(Marker[molt]);
      fHistResolution_1D[molt][m][t]->SetLineColor(Color[molt]);
      fHistResolution_1D[molt][m][t]->SetMarkerColor(Color[molt]);
      fHistResolution_1D[molt][m][t]->Draw("same");
      if(molt ==nummolt)     legend->Draw();
  
      cout << "molt: " << molt << " sto riempiendo histo " <<  nameRes_1D[molt][m][t]<< endl;
    }
    }

    //***************contamination factor for trigger and V0 in multiplicity intervals***************

    HistContTrigger[molt]=(TH2D*)list2->FindObject(Form("fHistPrimaryTrigger_%i_cut%i", molt, sysTrigger)); //histo tipo di trigger (primario, secondario) vs pT trigger
    cout << "deficit " << endl;
    HistContV0[molt]=(TH2D*)list2->FindObject(Form("fHistPrimaryV0_%i_cut%i", molt, 0));//histo tipo di V0 (primario, secondario) vs pT V0

    HistInt[molt]=HistContTrigger[molt]->ProjectionX("", HistContTrigger[molt]->GetYaxis()->FindBin(ptjmin+0.00001),HistContTrigger[molt]->GetYaxis()->FindBin(ptjmax-0.0001));
    if((HistInt[molt]->GetBinContent(1) + HistInt[molt]->GetBinContent(2) + HistInt[molt]->GetBinContent(3)+ HistInt[molt]->GetBinContent(4)) != 0){
      ContTrigger[molt]=1. - (Double_t)(HistInt[molt]->GetBinContent(1))/(HistInt[molt]->GetBinContent(1) + HistInt[molt]->GetBinContent(2) + HistInt[molt]->GetBinContent(3)+ HistInt[molt]->GetBinContent(4));
      HistContTriggerMolt->SetBinContent(molt+1, ContTrigger[molt]);
    }
    else     HistContTriggerMolt->SetBinContent(molt+1, 0);

    canvasUsed->cd(3);
    HistContTriggerMolt->GetYaxis()->SetRangeUser(0,0.05);
    HistContTriggerMolt->SetMarkerStyle(ColorSysTrigger[sysTrigger]);
    HistContTriggerMolt->SetLineColor(ColorSysTrigger[sysTrigger]);
    HistContTriggerMolt->SetMarkerColor(ColorSysTrigger[sysTrigger]);
    HistContTriggerMolt->Draw("e");
    HistContTriggerMolt->Draw();
    
    HistInt[molt]=      HistContV0[molt]->ProjectionX();
    ContV0[molt]=1. - (Double_t)(HistInt[molt]->GetBinContent(1))/(HistInt[molt]->GetBinContent(1) + HistInt[molt]->GetBinContent(2) + HistInt[molt]->GetBinContent(3)+ HistInt[molt]->GetBinContent(4));

    cout << " cont trigger " << ContTrigger[molt] << " cont v0 " << ContV0[molt] << endl;
      
    HistContTriggerPt[molt]=new TH1D("HistContTriggerPt_"+Smolt[molt], Form("HistContTriggerPt_%i", molt),HistContTrigger[molt]->GetNbinsY(),HistContTrigger[molt]->GetYaxis()->GetBinLowEdge(1), HistContTrigger[molt]->GetYaxis()->GetBinUpEdge(HistContTrigger[molt]->GetNbinsY()) );
    HistContTriggerPt[molt]->SetTitle("HistContTriggerPt_"+Smolt[molt]);
    HistContTriggerPt[molt]->GetXaxis()->SetTitle("p_T");
    HistContTriggerPt[molt]->GetYaxis()->SetTitle("C factor");

    HistContV0Pt[molt]=new TH1D("HistContV0Pt_"+Smolt[molt], Form("HistContV0Pt_%i", molt),HistContV0[molt]->GetNbinsY(),HistContV0[molt]->GetYaxis()->GetBinLowEdge(1), HistContV0[molt]->GetYaxis()->GetBinUpEdge(HistContV0[molt]->GetNbinsY()) );
    HistContV0Pt[molt]->SetTitle("HistContV0Pt_"+Smolt[molt]);
    HistContV0Pt[molt]->GetXaxis()->SetTitle("p_T");
    HistContV0Pt[molt]->GetYaxis()->SetTitle("C factor");

    HistContV0PtBins[molt]=new TH1D("HistContV0PtBins_"+Smolt[molt], Form("HistContV0PtBins_%i", molt),numPtV0, NPtV0 );
    HistContV0PtBins[molt]->SetTitle("HistContV0PtBins_"+Smolt[molt]);
    HistContV0PtBins[molt]->GetXaxis()->SetTitle("p_T");
    HistContV0PtBins[molt]->GetYaxis()->SetTitle("C factor");


    Double_t denom=0;
    Double_t num=0;

    //histo contamination factor for trigger particles in pt bins of trigger particle
    for(Int_t pt=0; pt<HistContTrigger[molt]->GetNbinsY(); pt++ ){
      num= (Double_t)(HistContTrigger[molt]->GetBinContent(1, pt+1)); //number of primary trigger particles
      denom= (HistContTrigger[molt]->GetBinContent(1, pt+1) + HistContTrigger[molt]->GetBinContent(2, pt+1) + HistContTrigger[molt]->GetBinContent(3, pt+1)+ HistContTrigger[molt]->GetBinContent(4, pt+1)); //total number of trigger particles
      if (denom!=0){
	ContTriggerInt=1. - num/denom;
	//	ContTriggerIntError=sqrt(pow(ContTriggerInt,2)*(denom-num)+pow((num-denom)/pow(denom,2),2)*num);
	//ContTriggerIntError=sqrt(pow((num-denom)/pow(denom,2),2)*num);
	HistContTriggerPt[molt]->SetBinContent(pt+1, ContTriggerInt);
	//	HistContTriggerPt[molt]->SetBinError(pt+1, ContTriggerIntError);
	cout << "pt bin " << pt << " cont factor trigger " << ContTriggerInt <<" +- " << ContTriggerIntError <<endl;
      }
    }


    //histo contamination factor for V0 particles in pt bins of V0 particle
    for(Int_t pt=0; pt<HistContV0[molt]->GetNbinsY(); pt++ ){
      num= (Double_t)(HistContV0[molt]->GetBinContent(1, pt+1)); //number of primary V0 
      denom= (HistContV0[molt]->GetBinContent(1, pt+1) + HistContV0[molt]->GetBinContent(2, pt+1) + HistContV0[molt]->GetBinContent(3, pt+1)+ HistContV0[molt]->GetBinContent(4, pt+1)); //total number of V0
      if(denom!=0){
	ContV0Int =1. - num/denom;
	//ContV0IntError=sqrt(pow(ContV0Int,2)*(denom-num)+pow((num-denom)/pow(denom,2),2)*num);
	HistContV0Pt[molt]->SetBinContent(pt+1, ContV0Int);
	//HistContV0Pt[molt]->SetBinContent(pt+1, ContV0IntError);
	cout << "pt bin " << pt << " cont factor V0 " << ContV0Int <<" +- " << ContV0IntError << endl;
      }
      else{
	HistContV0Pt[molt]->SetBinContent(pt+1, 0);
	HistContV0Pt[molt]->SetBinContent(pt+1, 0);
      } 
    }

    cout << "\n\n contamination factor for V0 particles in Pt bins " << endl;
    //histo contamination factor for V0 particles in pt bins of V0 particle (bins used in analysis)
    for(Int_t binpt=0; binpt<numPtV0; binpt++){
      num=0;
      denom=0;
      for(Int_t pt=HistContV0[molt]->GetYaxis()->FindBin(NPtV0[binpt]+0.0001); pt<HistContV0[molt]->GetYaxis()->FindBin(NPtV0[binpt+1]-0.0001); pt++){
	num+= (Double_t)(HistContV0[molt]->GetBinContent(1, pt+1)); //number of primary V0 
	denom+= (HistContV0[molt]->GetBinContent(1, pt+1) + HistContV0[molt]->GetBinContent(2, pt+1) + HistContV0[molt]->GetBinContent(3, pt+1)+ HistContV0[molt]->GetBinContent(4, pt+1)); //total number of V0
      }
      if(denom!=0){
	ContV0Int =1. - num/denom;
	//	ContV0IntError=sqrt(pow(ContV0Int,2)*(denom-num)+pow((num-denom)/pow(denom,2),2)*num);
	HistContV0PtBins[molt]->SetBinContent(binpt+1, ContV0Int);
	//HistContV0PtBins[molt]->SetBinContent(binpt+1, ContV0IntError);
	cout << "pt bin " << binpt << " cont factor V0 " << ContV0Int <<" +- " << ContV0IntError << endl;
      }
      else{
	HistContV0PtBins[molt]->SetBinContent(binpt+1, 0);
	
      }
    }


    canvasUsed->cd(4);
    HistContV0PtBins[molt]->GetYaxis()->SetRangeUser(0,1);
    HistContV0PtBins[molt]->SetMarkerStyle(Marker[molt]);
    HistContV0PtBins[molt]->SetLineColor(Color[molt]);
    HistContV0PtBins[molt]->SetMarkerColor(Color[molt]);
    HistContV0PtBins[molt]->Draw("samee");
    if(molt ==nummolt)     legend->Draw();


    canvasCont->cd(1);
    HistContTriggerPt[molt]->Rebin(2);
    HistContTriggerPt[molt]->GetYaxis()->SetRangeUser(0, 0.5);
    HistContTriggerPt[molt]->SetMarkerStyle(Marker[molt]);
    HistContTriggerPt[molt]->SetLineColor(Color[molt]);
    HistContTriggerPt[molt]->SetMarkerColor(Color[molt]);
    HistContTriggerPt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();


    canvasCont->cd(2);
    HistContV0Pt[molt]->Rebin(2);
    HistContV0Pt[molt]->GetYaxis()->SetRangeUser(0, 0.5);
    HistContV0Pt[molt]->SetMarkerStyle(Marker[molt]);
    HistContV0Pt[molt]->SetLineColor(Color[molt]);
    HistContV0Pt[molt]->SetMarkerColor(Color[molt]);
    HistContV0Pt[molt]->Draw("same");
    if(molt ==nummolt)     legend->Draw();

  
  }


 
  cout << "\n \n sto per scriver sul file efficienza " << endl;
  
  TFile *fileoutbis;  
  fileoutbis = new TFile(PathOut2, "RECREATE");
  // fHistSelected_2D_TriggerPtPhi->Write();
  // fHistGenerated_2D_TriggerPtPhi->Write();

/*
  fHistTriggerEfficiencyPtPhiMolt->Write();
  fHistTriggerEfficiencyPtEtaMolt->Write();
  fHistSelected_1D_TriggerPt->Write();
  fHistGenerated_1D_TriggerPt->Write();
  fHistTriggerEfficiencyPtMolt->Write();
  fHistSelected_1D_TriggerPhi->Write();
  fHistGenerated_1D_TriggerPhi->Write();
  fHistTriggerEfficiencyPhiMolt->Write();
  fHistSelected_1D_TriggerEta->Write();
  fHistGenerated_1D_TriggerEta->Write();
  fHistTriggerEfficiencyEtaMolt->Write();
    
  fHistV0EfficiencyPtPhiMolt->Write();
  fHistV0EfficiencyPtEtaMolt->Write();
  fHistSelected_1D_V0Pt->Write();
  fHistGenerated_1D_V0Pt->Write();
  fHistV0EfficiencyPtMolt ->Write();
  fHistSelected_1D_V0Phi->Write();
  fHistGenerated_1D_V0Phi->Write();
  fHistV0EfficiencyPhiMolt->Write();
  fHistSelected_1D_V0Eta->Write();
  fHistGenerated_1D_V0Eta->Write();
  fHistV0EfficiencyEtaMolt->Write();

  fHistSelectedMass_1D->Write();
  fHistRecoMass_1D->Write();
  fHistV0EfficiencyRecoMolt->Write();
  fHistV0EfficiencyRecoPtMolt->Write();
*/
  HistoTriggerEfficiency->Write();
  HistContTriggerMolt->Write();   


  for(Int_t m=0; m<nummolt+1; m++){

    fHistTriggerEfficiencyPtPhi[m]->Write();
    fHistTriggerEfficiencyPtEta[m]->Write();
    fHistSelected_1D_TriggerPt[m]->Write();
    fHistGenerated_1D_TriggerPt[m]->Write();
    fHistTriggerEfficiencyPt[m]->Write();
    fHistSelected_1D_TriggerPhi[m]->Write();
    fHistGenerated_1D_TriggerPhi[m]->Write();
    fHistTriggerEfficiencyPhi[m]->Write();
    fHistSelected_1D_TriggerEta[m]->Write();
    fHistGenerated_1D_TriggerEta[m]->Write();
    fHistTriggerEfficiencyEta[m]->Write();
      
    fHistV0EfficiencyPtPhi[m]->Write();
    fHistV0EfficiencyPtEta[m]->Write();
    fHistSelected_1D_V0Pt[m]->Write();
    fHistGenerated_1D_V0Pt[m]->Write();
    fHistV0EfficiencyPt[m]->Write();
    fHistSelected_1D_V0Phi[m]->Write();
    fHistGenerated_1D_V0Phi[m]->Write();
    fHistV0EfficiencyPhi[m]->Write();
    fHistSelected_1D_V0Eta[m]->Write();
    fHistGenerated_1D_V0Eta[m]->Write();
    fHistV0EfficiencyEta[m]->Write();
    
    fHistV0EfficiencyPtBins[m] ->Write();
      
    fHistV0EfficiencyReco[m]->Write();
    fHistV0EfficiencyRecoPt[m]->Write();
     
  
    HistContV0[m]->Write();
    HistContTrigger[m]->Write();
    HistContV0Pt[m]->Write();
    HistContV0PtBins[m]->Write();
    HistContTriggerPt[m]->Write();
  

    cout << m << endl;
    for(Int_t y=0; y< 3; y++){
      for(Int_t t=0; t< 2; t++){
	fHistResolution_1D[m][y][t]->Write();
      }
    }

  }
  fileoutbis->Close();
    
  
  
  cout << "******************************************************************"<< endl;
  cout << "partendo dai file "  << PathInBis << " ho creato: "<< endl;
  cout << "il file " << PathOut2 << endl;
}
