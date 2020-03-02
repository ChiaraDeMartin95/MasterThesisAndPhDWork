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
  return sqrt(((Float_t)k+1)*((Float_t)k+2)/(n+2)/(n+3) - pow((Float_t)(k+1),2)/pow(n+2,2));
}

void CfrEfficiencyCascades(TString path1="15runsNOFB_Try14"){
  const Int_t numSel= 22;
  const Int_t numPtV0=25;//36;

  TH2F *fHistEventXiTrue[4];

  TH1F *fHistEventXiTrue1D[numSel+1][4];
  TH1F *fHistEventXiTrue1DPtBins[numSel+1][4];
  TH3F *hGeneratedMultvsPt[4];
  TH2F* hGeneratedMultvsPt2D[4];
  TH1F *hGeneratedPt1DPtBins[4]; 
  TH1F *hGeneratedPt1D[4]; 
  TH1F *histo_Efficiency[numSel+1][4]; 
  TH1F *histo_EfficiencyFinal[4]; 

  
  Int_t Color[numSel+1]={1,401,801,628, 909,881, 860, 868, 841, 418, 1,401,801,628, 909,881, 860, 868, 841, 418,1, 401, 629};
  Int_t Marker[3]={22, 32, 30};
  //  Double_t binl[numPtV0+1]={0,0.4, 1,1.5, 2,2.5,3,4,8};
  Double_t binl[numPtV0+1]={0, 0.4,0.6, 0.8, 1,1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 3.8, 4.2, 4.6, 5.0,5.4, 5.8, 6.2, 6.6, 7.0, 8.0};// 5.4, 5.8, 6.2, 6.6, 7., 7.4, 8};

  auto legend = new TLegend(0.6, 0.1, 0.9, 0.4);
  legend->SetHeader("Progressive selections");     

  TString stringin= "FinalOutput/AnalysisResults2018f1_extra_MCEff_Cascades_"+ path1 + ".root";
  TString nome_TDir="MyTask";
  TFile *filein= new TFile(stringin, "");
  TDirectoryFile *dirinput = (TDirectoryFile*)filein->Get(nome_TDir);
  TList *listinput = (TList*)dirinput->Get("MyOutputContainer");
  TList *listinputGen = (TList*)dirinput->Get("MyOutputContainer3");
  TString stringout = "CfrEfficiencyCascadesOutput_"+path1+".root";
  TFile *fileout= new TFile(stringout, "RECREATE");
  Float_t 	nEventXiTrue1D=0;

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1500,1000);
  canvas->Divide(2,2);

  //Prendo histo cascades generate*********************************************************
  TString filewithfinaleff[4];
  TFile *  filefinaleff[4];
  filewithfinaleff[0] = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_MCEff_Cascades_"+path1+"_2018f1_extra_XiNeg_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0_Rap0.root";
  filewithfinaleff[1] = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_MCEff_Cascades_"+path1+"_2018f1_extra_XiPos_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0_Rap0.root";
  filewithfinaleff[2] = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_MCEff_Cascades_"+path1+"_2018f1_extra_XiNeg_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0_Rap1.root";
  filewithfinaleff[3] = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_MCEff_Cascades_"+path1+"_2018f1_extra_XiPos_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0_Rap1.root";
  hGeneratedMultvsPt[0] = (TH3F*)listinputGen->FindObject("fHistGeneratedXiPt_0");
  hGeneratedMultvsPt[1] = (TH3F*)listinputGen->FindObject("fHistGeneratedXiPt_0");
  hGeneratedMultvsPt[2] = (TH3F*)listinputGen->FindObject("fHistGeneratedXiPt_0");
  hGeneratedMultvsPt[3] = (TH3F*)listinputGen->FindObject("fHistGeneratedXiPt_0");
  fHistEventXiTrue[0] = (TH2F*)listinput ->FindObject("fHistEventXiTrueNeg");
  fHistEventXiTrue[1] = (TH2F*)listinput ->FindObject("fHistEventXiTruePos");
  fHistEventXiTrue[2] = (TH2F*)listinput ->FindObject("fHistEventXiTrueNegRapSel");
  fHistEventXiTrue[3] = (TH2F*)listinput ->FindObject("fHistEventXiTruePosRapSel");
  for (Int_t c=0; c<4; c++){
    if (!hGeneratedMultvsPt[c]) {cout << "the generated histogram " << c << " does not exist" << endl; return;}
    if (!fHistEventXiTrue[c]) {cout << "the selected histogram " << c << " does not exist" << endl; return;}
    filefinaleff[c] = new TFile(filewithfinaleff[c], "");
    histo_EfficiencyFinal[c] = (TH1F*)filefinaleff[c]->Get("histo_Efficiency");
    if (    histo_EfficiencyFinal[c])        histo_EfficiencyFinal[c]->SetName(Form("histo_Efficiency_c%i", c));
    if (    histo_EfficiencyFinal[c])    histo_EfficiencyFinal[c] ->SetLineColor(1);
  }

  

  for (Int_t c=0; c<8; c++){

    // if (c==0)    cout <<"Let's start with XiNeg with no rapidity selection " << endl;
    // else if (c==1)    cout <<"Let's start with XiPos with no rapidity selection " << endl;
    // else if (c==2)    cout <<"Let's start with XiNeg with |y|<0.5 " << endl;
    // else   cout <<"Let's start with XiPos with |y|<0.5" << endl;

    if (c==0 || c==2 || c==4 || c==6){
      hGeneratedMultvsPt[c]->GetXaxis()->SetRangeUser( ,-0.00001);
    }
    else {
      hGeneratedMultvsPt[c]->GetXaxis()->SetRangeUser(0.00001, );
    }

    hGeneratedPt1DPtBins[c] = new TH1F (Form("hGeneratedPt1DPtBins_%i",c),"Generated vs Pt", numPtV0, binl);
    hGeneratedMultvsPt2D[c] = (TH2F*)hGeneratedMultvsPt[c]->Project3D("yxo");
    hGeneratedPt1D[c] = (TH1F*)hGeneratedMultvsPt2D[c]->ProjectionX(Form("hGeneratedPt1D_%i", c),hGeneratedMultvsPt2D[c]->GetYaxis()->FindBin(0+0.001),hGeneratedMultvsPt2D[c]->GetYaxis()->FindBin(100-0.001)); 

    cout << " I get the generated histogram " << endl;
    //creo istogramma generati in bin di pt**********************************
    Float_t  NumberOfGenerated=0;
    for(Int_t ptv0=0; ptv0<numPtV0; ptv0++){
      NumberOfGenerated=0;
      for(Int_t i=hGeneratedPt1D[c]->GetXaxis()->FindBin(binl[ptv0]+0.00001); i<=hGeneratedPt1D[c]->GetXaxis()->FindBin(binl[ptv0+1]-0.00001); i++){
	NumberOfGenerated+=hGeneratedPt1D[c]->GetBinContent(i);
      }
      hGeneratedPt1DPtBins[c]->SetBinContent(ptv0+1, NumberOfGenerated);
      if (ptv0==0) {
	hGeneratedPt1DPtBins[c]->SetBinContent(ptv0+1, 0);
      }
    }
 
    for (Int_t l=1; l<=numSel; l++){
      if (l<=3) continue;
      cout << " Now I get the selected histograms for the different selections------------> selection n. " << l << endl;
      cout <<   fHistEventXiTrue[0]->GetXaxis()->GetBinLabel(l)<< endl;
      histo_Efficiency[l][c]= new TH1F (Form("histo_Efficiency_%i_sel%i", c, l),"Efficiency vs Pt", numPtV0, binl);
      histo_Efficiency[l][c]->SetTitle(fHistEventXiTrue[0]->GetXaxis()->GetBinLabel(l));
      if (c==0)      legend->AddEntry(      histo_Efficiency[l][c],fHistEventXiTrue[0]->GetXaxis()->GetBinLabel(l), "pl");
      //creo istogramma selezionati in bin di pt**********************************
      fHistEventXiTrue1D[l][c]= (TH1F*) fHistEventXiTrue[c]->ProjectionY(Form("fHistEventXiTrue%i_sel%i", c,l), l, l);
      if (!      fHistEventXiTrue1D[l][c]) cout << "projected  histogram has some problems " << endl;
      fHistEventXiTrue1DPtBins[l][c] = new TH1F (Form("fHistEventXiTrue1DPtBins_%i_sel%i",c, l),"Selected vs Pt", numPtV0, binl);
      fHistEventXiTrue1DPtBins[l][c]->SetTitle(fHistEventXiTrue[0]->GetXaxis()->GetBinLabel(l));
      for (Int_t ptv0=0; ptv0< numPtV0; ptv0++){
	nEventXiTrue1D=0;
	 for (Int_t pt=fHistEventXiTrue1D[l][c]->GetXaxis()->FindBin(binl[ptv0]+0.0001); pt <=fHistEventXiTrue1D[l][c]->GetXaxis()->FindBin(binl[ptv0+1]-0.0001); pt++){
	   nEventXiTrue1D+=fHistEventXiTrue1D[l][c]->GetBinContent(pt);
	 }
	fHistEventXiTrue1DPtBins[l][c]->SetBinContent(ptv0+1, nEventXiTrue1D);
	if (ptv0==0) {
	  fHistEventXiTrue1DPtBins[l][c] ->SetBinContent(ptv0+1, 0);
	}
	//	cout << "bin content of bin " << binl[ptv0]<< "-" << binl[ptv0+1] << " : " <<	  fHistEventXiTrue1DPtBins[l][c]->GetBinContent(ptv0+1)<<endl;
      }

    //creo istogramma efficienze in bin di pt************************************************
      cout << " Now I compute the efficiency for the different selections------------> selection n. " << l << endl;
      for(Int_t i=0; i <histo_Efficiency[l][c]->GetNbinsX(); i++ ){
	if (hGeneratedPt1DPtBins[c]->GetBinContent(i+1) !=0){
	  histo_Efficiency[l][c]->SetBinContent(i+1,((Float_t)fHistEventXiTrue1DPtBins[l][c]->GetBinContent(i+1)/ hGeneratedPt1DPtBins[c]->GetBinContent(i+1)));
	  histo_Efficiency[l][c]->SetBinError(i+1, SetEfficiencyError(fHistEventXiTrue1DPtBins[l][c]->GetBinContent(i+1), hGeneratedPt1DPtBins[c]->GetBinContent(i+1)));

	}
	else {
	  histo_Efficiency[l][c]->SetBinContent(i+1,0);
	  histo_Efficiency[l][c]->SetBinError(i+1, 0);
	}
	if (histo_Efficiency[l][c]->GetBinContent(i+1)==0) 	  histo_Efficiency[l][c]->SetBinError(i+1,0);
	//	cout << "bin content of bin " << binl[i]<< "-" << binl[i+1] << " : " <<histo_Efficiency[l][c]->GetBinContent(i+1)<<endl;
      }

      if (l<10)    histo_Efficiency[l][c]->SetMarkerStyle(33);
      else if (l<20)   histo_Efficiency[l][c]->SetMarkerStyle(23);
      else    histo_Efficiency[l][c]->SetMarkerStyle(30);
       histo_Efficiency[l][c]->SetMarkerSize(1.2);
       histo_Efficiency[l][c]->SetLineColor(Color[l]);
       histo_Efficiency[l][c]->SetMarkerColor(Color[l]);
      fHistEventXiTrue1DPtBins[l][c]->Write();
      histo_Efficiency[l][c]->Write();

    } //chiudo loop diverse selezioni
    hGeneratedPt1DPtBins[c]->Write();

    canvas->cd(c+1);
    if (histo_EfficiencyFinal[c]) {
     histo_EfficiencyFinal[c]->Draw("same pl");
       histo_EfficiencyFinal[c]->SetMarkerSize(3);
      histo_EfficiencyFinal[c]->GetYaxis()->SetRangeUser(0,1.3);
      if (c==0)      histo_EfficiencyFinal[c]->SetTitle("Efficiency vs Pt Xi Neg No selection on y");
      else       if (c==1)      histo_EfficiencyFinal[c]->SetTitle("Efficiency vs Pt Xi Pos No selection on y");
      else if (c==2)      histo_EfficiencyFinal[c]->SetTitle("Efficiency vs Pt Xi Neg |y|<0.5");
      else      histo_EfficiencyFinal[c]->SetTitle("Efficiency vs Pt Xi Pos |y|<0.5");

    }
    for (Int_t l=numSel; l>=4; l--){
      if (c==0)      histo_Efficiency[l][c]->SetTitle("Efficiency vs Pt Xi Neg No selection on y");
      else       if (c==1)      histo_Efficiency[l][c]->SetTitle("Efficiency vs Pt Xi Pos No selection on y");
      else if (c==2)      histo_Efficiency[l][c]->SetTitle("Efficiency vs Pt Xi Neg |y|<0.5");
      else      histo_Efficiency[l][c]->SetTitle("Efficiency vs Pt Xi Pos |y|<0.5");
      cout << l << endl;
      cout <<       histo_Efficiency[l][c]->GetBinContent(2)<< endl;
      histo_Efficiency[l][c]->Draw("same pe");
      histo_Efficiency[l][c]->GetYaxis()->SetRangeUser(0,1.2);
      legend->Draw("");
      //      histo_Efficiency[4][0]->Draw("same");
    }
  }//chiudo loop su XiPos/XiNeg e diverse rapidità

  fileout->cd();
  fileout->WriteTObject(canvas);
  fileout->Close();
  cout << "\n\npartendo dal file " << stringin << " ho creato il file " << stringout << endl;

  /*
  TCanvas *canvasbis= new TCanvas("canvasbis", "canvasbis", 1000, 800);
  canvasbis->Divide(4,1);
  for (Int_t c=0; c<4; c++){
     cout << " c " << c << endl;
     canvasbis->cd(c+1);
    //    for (Int_t l=1; l<=numSel; l++){
    for (Int_t l=1; l<=4; l++){
      if (l<=3) continue;
      cout << " l " << l << endl;
      // histo_Efficiency[l][c]->Draw("same pe");
           cout <<       histo_Efficiency[l][c]->GetBinContent(2)<< endl;
    }
  }
  */
 
}


