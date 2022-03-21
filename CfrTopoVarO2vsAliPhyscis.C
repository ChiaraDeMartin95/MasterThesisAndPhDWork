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
#include <Macros/ErrRatioCorr.C>
#include <Macros/BarlowVariable.C>

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Float_t mSize, Int_t linestyle){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetLineStyle(linestyle);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  //  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  //  histo->GetYaxis()->SetTitleOffset(yOffset); //1.2                                          
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
  histo->SetLineWidth(1);
}

const Int_t numTopoVar = 15;

void CfrTopoVarO2vsAliPhyscis(Bool_t isScale=0, TString PathIn ="AnalysisResultsESDMultiStrangeQA.root"/*"APQA.root"/*"outputAliPhysics_LHC20i2b_288909.root"*/, TString PathInO2 = "O2output_20i2bv4.root"/*"O2output_18ppbig.root"/*"O2output_LHC20i2b_Ter.root" /*"TriggerOutput_loose.root"/*"O2QA.root"/*"outputsO2Physics_LHC20i2b_288909.root"*/, Int_t AliPhysicsOrO2=0){

  //AliPhysicsOrO2 == 0 -> O2 vs AliPhysics comparison
  //AliPhysicsOrO2 == 1 -> only O2 plots
  //AliPhysicsOrO2 == 2 -> only AliPhysics plots

  TFile * filein = new TFile (PathIn, "");
  if (!filein) {cout << "file AliPhysics not present " << endl; return;}
  TDirectoryFile * dir = (TDirectoryFile*) filein->Get("PWGLFStrangeness.outputCheckCascade");
  if (!dir)  {cout << "dir AliPhysics not present " << endl; return;}
  TList * list = (TList*)dir->Get("fListHistMultistrangev2QA");
  if (!list)  {cout << "list AliPhysics not present " << endl; return;}

  TFile * fileinO2 = new TFile (PathInO2, "");
  if (!fileinO2) {cout << "file O2 not present " << endl; return;}
  TDirectoryFile * dirO2 = (TDirectoryFile*) fileinO2->Get("cascade-qa");
  if (!dirO2)  {cout << "dir O2 not present " << endl; return;}
  //  TDirectoryFile * listO2 = (TDirectoryFile*) dirO2->Get("QAHistos");
  //  if (!listO2)  {cout << "list O2 not present " << endl; return;}

  TString STopoVar[numTopoVar] = { "hV0Radius", "hCascRadius", "hV0CosPA", "hCascCosPA", "hDCAPosToPV", "hDCANegToPV", "hDCABachToPV", "hDCAV0ToPV", "hDCAV0Dau", "hDCACascDau", "hLambdaMass", "hMassXiMinus", "hMassOmegaMinus", "hMassXiPlus", "hMassOmegaPlus"};
  TString titleX[numTopoVar] = {"hV0Radius", "hCascRadius", "hV0CosPA", "hCascCosPA", "hDCAPosToPV", "hDCANegToPV", "hDCABachToPV", "hDCAV0ToPV", "hDCAV0Dau", "hDCACascDau", "hLambdaMass", "hMassXiMinus", "hMassOmegaMinus", "hMassXiPlus", "hMassOmegaPlus",};
  TH1F* hTopoVarAliPhysics[numTopoVar]; 
  TH1F* hTopoVarO2[numTopoVar]; 
  TH1F*     hLambdaMassFromXi;
  TH1F*     hLambdaMassFromOmega;

  gStyle->SetOptStat(0);
  TCanvas * canvas[4];
  for (Int_t i=0; i<4; i++){
    canvas[i]  = new TCanvas (Form("canvas%i", i), Form("canvas%i",i), 1600, 1000);
    canvas[i]->Divide(2,2);
  }

  for (Int_t i=0; i<numTopoVar; i++){
    cout << STopoVar[i]<< endl;
   hTopoVarAliPhysics[i] = (TH1F*)list->FindObject(STopoVar[i]);
   if (!hTopoVarAliPhysics[i]) {cout << "Histogram " << STopoVar[i] << " is not present!" << endl; return;}
   if (STopoVar[i] == "hLambdaMass") {
     hTopoVarAliPhysics[i]->Rebin(1);
   }
   else if (STopoVar[i].Index("Mass") !=-1)  hTopoVarAliPhysics[i]->Rebin(2);  
   else    hTopoVarAliPhysics[i]->Rebin(8); 
   if (isScale)   hTopoVarAliPhysics[i]->Scale(1./hTopoVarAliPhysics[i]->GetEntries());
   if (STopoVar[i] == "hMassXiMinus" || STopoVar[i] == "hMassXiPlus")   hTopoVarAliPhysics[i]->GetXaxis()->SetRangeUser(1.2, 1.4);
   if (STopoVar[i] == "hMassOmegaMinus" || STopoVar[i] == "hMassOmegaPlus")   hTopoVarAliPhysics[i]->GetXaxis()->SetRangeUser(1.5, 1.8);
   //   if (STopoVar[i] == "hCascCosPA" || STopoVar[i] == "hV0CosPA")   hTopoVarAliPhysics[i]->GetXaxis()->SetRangeUser(0.8,1);
   if (STopoVar[i] == "hLambdaMass")   hTopoVarAliPhysics[i]->GetXaxis()->SetRangeUser(1,1.2);
   cout << hTopoVarAliPhysics[i]->GetNbinsX() << endl;
   StyleHistoYield(hTopoVarAliPhysics[i], 10e-4, 1.3*hTopoVarAliPhysics[i]->GetMaximum(), 628, 1, titleX[i], "", STopoVar[i], 2, 1);

   hTopoVarO2[i] = (TH1F*)dirO2->Get(STopoVar[i]);
   if (!hTopoVarO2[i]) {cout << "Histogram " << STopoVar[i] << " is not present!" << endl; return;}
   //   if (STopoVar[i] == "hCascCosPA" || STopoVar[i] == "hV0CosPA")   hTopoVarO2[i]->Rebin(20);
   if (STopoVar[i] == "hLambdaMass")  {
     hTopoVarO2[i]->Rebin(1); 
     hLambdaMassFromXi = (TH1F*)dirO2->Get("hLambdaMassFromXi");
     hLambdaMassFromOmega = (TH1F*)dirO2->Get("hLambdaMassFromOmega");
     hLambdaMassFromXi->Scale(3);
     hLambdaMassFromOmega->Scale(3);
     hLambdaMassFromXi->Rebin(4);
     hLambdaMassFromOmega->Rebin(4);
   }
   else if (STopoVar[i].Index("Mass") !=-1)  hTopoVarO2[i]->Rebin(2);  
   else   hTopoVarO2[i]->Rebin(8);
   if (isScale)   hTopoVarO2[i]->Scale(1./hTopoVarO2[i]->GetEntries());
   if (STopoVar[i] == "hMassXiMinus" || STopoVar[i] == "hMassXiPlus")   hTopoVarO2[i]->GetXaxis()->SetRangeUser(1.2, 1.4);
   if (STopoVar[i] == "hMassOmegaMinus" || STopoVar[i] == "hMassOmegaPlus")   hTopoVarO2[i]->GetXaxis()->SetRangeUser(1.5, 1.8);
   //   if (STopoVar[i] == "hCascCosPA" || STopoVar[i] == "hV0CosPA")   hTopoVarO2[i]->GetXaxis()->SetRangeUser(0.8,1);
   if (STopoVar[i] == "hLambdaMass")   hTopoVarO2[i]->GetXaxis()->SetRangeUser(1,1.2);
   cout << hTopoVarO2[i]->GetNbinsX() << endl;
   StyleHistoYield(hTopoVarO2[i], 10e-4, 1.3*hTopoVarO2[i]->GetMaximum(), 9, 1, titleX[i], "", STopoVar[i], 2, 1);

   TLegend * legend = new TLegend (0.6, 0.6, 0.9, 0.9);
   if (isScale){
     legend->AddEntry(hTopoVarAliPhysics[i], Form("AliPhysics, entries: %i", (int)hTopoVarAliPhysics[i]->GetEntries()), "l");
     legend->AddEntry(hTopoVarO2[i], Form("O2, entries: %i", (int)hTopoVarO2[i]->GetEntries()), "l");
   }
   else {
     legend->AddEntry(hTopoVarAliPhysics[i], Form("AliPhysics, entries: %i - integral: %i", (int)hTopoVarAliPhysics[i]->GetEntries(), (int)hTopoVarAliPhysics[i]->Integral()), "l");
     legend->AddEntry(hTopoVarO2[i], Form("O2, entries: %i - integral: %i", (int)hTopoVarO2[i]->GetEntries(), (int)hTopoVarO2[i]->Integral()), "l");
   }

   if (i<4) {
     canvas[0]->cd(i+1);
   }
   else if (i<8) {
     canvas[1]->cd(i-4+1);
   }
   else if (i<11) {
     canvas[2]->cd(i-8+1);
   }
   else {
     canvas[3]->cd(i-11+1);
   }
   cout << hTopoVarAliPhysics[i]->GetMaximum() << " " <<hTopoVarO2[i]->GetMaximum()<<endl;
   if (STopoVar[i] == "hCascCosPA" || STopoVar[i] == "hV0CosPA") gPad->SetLogy();

   if (AliPhysicsOrO2==1){
     hTopoVarO2[i]->Draw("hist");
   }
   else if (AliPhysicsOrO2==2){
     hTopoVarAliPhysics[i]->Draw("hist");
   }
   else {
     if (hTopoVarAliPhysics[i]->GetMaximum()>hTopoVarO2[i]->GetMaximum()){
       hTopoVarAliPhysics[i]->Draw("hist");
       hTopoVarO2[i]->Draw("hist same");
     }
     else {
       hTopoVarO2[i]->Draw("hist");
       hTopoVarAliPhysics[i]->Draw("hist same");
     }
   }
   if (hLambdaMassFromXi && hLambdaMassFromOmega && STopoVar[i] == "hLambdaMass"){
     hLambdaMassFromXi->SetLineColor(881);
     hLambdaMassFromOmega->SetLineColor(807);
     hLambdaMassFromXi->Draw("hist same");
     hLambdaMassFromOmega->Draw("hist same");
   }
   legend->Draw("");
  }

  TString outputName = "CfrO2vsAliPhysics_Casc";
  if (isScale) outputName += "_Norm";
  if (AliPhysicsOrO2==1) outputName += "_O2Only";
  else if (AliPhysicsOrO2==2) outputName += "_AliPhysicsOnly";
  outputName += ".pdf";
  canvas[0]->SaveAs(outputName+"(");
  canvas[1]->SaveAs(outputName);
  canvas[2]->SaveAs(outputName);
  canvas[3]->SaveAs(outputName+")");

  cout << "PathIn: " << PathIn << endl;
  cout << "PathInO2: " << PathInO2 << endl;
}
