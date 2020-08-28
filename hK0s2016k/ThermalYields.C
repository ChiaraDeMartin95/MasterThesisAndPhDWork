#include "Riostream.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLine.h"
//#include "TRatioPlot.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TSystem.h"
#include "Math/IFunction.h"
//#include "SpecFuncMathMore.h"
//#include "/Applications/root6/root/math/mathmore/src/SpecFuncMathMore.cxx"
//#include "/Applications/root6/root/math/mathmore/inc/Math/SpecFuncMathMore.h"
#include "Math/SpecFuncMathMore.h"


TH1F * CreateRatio(TH1F * hNum, TH1F * hDenom, Float_t ErrDenom[], TString NameRatio){ //ratio between hNum and hDenom; num has no error
  TH1F * hRatio = (TH1F*)hNum->Clone(NameRatio);
  hRatio->Divide(hDenom);
  for (Int_t b=1; b<=hRatio->GetNbinsX(); b++){
    hRatio->SetBinError(b, ErrDenom[b-1]*hNum->GetBinContent(b)/pow(hDenom->GetBinContent(b),2));
  }
  return hRatio;
}

TH1F* HistoStyle(TH1F* hdNdyThermal,  Int_t numPart, Float_t dNdyThermal[], TString Part[], Int_t Color, TLegend *legendParVar, TString VarString){
  hdNdyThermal->SetTitle("Yields predicted by the thermal model in 0-10% PbPb collisions at #sqrt{s}=2.76 TeV");
  for (Int_t b=1; b<=hdNdyThermal->GetNbinsX(); b++){
    hdNdyThermal->SetBinContent(b, dNdyThermal[b-1]);
    hdNdyThermal->GetXaxis()->SetBinLabel(b, Part[b-1]);
  }
  hdNdyThermal->SetMarkerColor(Color);
  hdNdyThermal->SetLineColor(Color);
  hdNdyThermal->SetMarkerStyle(33);
  legendParVar->AddEntry(  hdNdyThermal, VarString, "pl");
  return   hdNdyThermal;
}

float YieldComponent(Int_t k, Float_t T, Float_t V, Float_t mu, Float_t mass, Float_t g, Bool_t IsFermion){
  gSystem->Load("libMathMore");
  //I use the modified Bessel function of the second kind (also called irregular modified cyilindrical Bessel function) [BesselK[nu,z] in the Wolfram Language]
  Float_t MeanN =0;
  TF1 * bessel = new TF1 ("bessel", "ROOT::Math::cyl_bessel_k([0], x)", 0, 10);
  bessel->SetParameter(0,2);
  //    TF1 * bessel = new TF1 ("bessel", "sin(x)", 0, 10);
  cout << " bessel " << bessel->Eval(k*mass/T)<< endl;
  Float_t Planckhbarc= 197.326980;
  if (IsFermion) MeanN = V*T*g/pow(Planckhbarc,3)/(2*pow(TMath::Pi(),2)) * pow(-1,k-1)/k*pow(exp(mu/T),k) * pow(mass,2) * bessel->Eval(k*mass/T);
  else MeanN = V*T*g/pow(Planckhbarc,3)/(2*pow(TMath::Pi(),2)) * pow(+1,k-1)/k*pow(exp(mu/T),k) * pow(mass,2) * bessel->Eval(k*mass/T);
  //  cout << V*T*g/pow(Planckhbarc,3) << " " << 1./(2*pow(TMath::Pi(),2)) << " " <<  pow(-1,k) << " " << 1./k*pow(exp(mu/T),k)<< " " << pow(mass,2)<< endl;
  return MeanN;
}

void ThermalYields(){

  //the publsihed yields are relative to PbPb collisions 0-10% and might be obsolete (dating back to 2015?)
  const Int_t numPart = 15;
  const Int_t HighestK=30;
  const Int_t NGraph= 8;
  const Int_t numParVariations=6;

  Float_t  VVar[numParVariations] = {5400, 4000, 5000, 6000, 5400, 5400};
  Float_t  TVar[numParVariations] = {156, 156, 156, 156, 165, 175};

  TString Part[numPart]={"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "K^{0}_{S}", "K^{*0}", "#phi", "p", "#bar{p}", "#Lambda", "#Xi^{+}", "#Xi^{-}", "#Omega^{+}", "#Omega^{-}", "d"};
  TString PartSmall[NGraph]={ "#pi", "K", "#phi", "p", "#Lambda","#Xi", "#Omega",  "d"};
  Float_t x[numPart] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

  Float_t dNdyPub[numPart]={669.5, 669, 100, 99.5, 98, 19.99, 12.75616, 31, 30.5, 24.1, 3.28, 3.34, 0.60, 0.58, 0.117};
  Float_t ErrdNdyPub[numPart]={48, 47, 8, 8.5, 8, 4.23, 1.49947, 2.5, 2.5, 2.6, 0.28, 0.25, 0.11, 0.10, 0.012};
  
  Float_t T = 156; //MeV
  Float_t V = 5400; //fm^3
  Float_t n[numPart] = {0};
  Float_t mass[numPart] = {139.57, 139.57, 493.677, 493.677, 497.611, 895.81, 1019.46, 938.27, 938.27, 1115.683, 1321.71, 1321.71, 1672.45, 1672.45, 1865.613};
  Float_t g[numPart] = {0}; //fattore di degenrazione di spin
  Float_t gSmall[NGraph] = {0}; //fattore di degenrazione di spin
  Float_t spin[numPart] = {0, 0, 0,0,0,1,1,0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1}; 
  Bool_t  IsFermion[numPart] = {0};
  Float_t dNdyThermal[numPart]={0};
  Float_t dNdyThermalDef[numPart]={0};
  Float_t dNdyThermalV4000[numPart]={0};
  Float_t dNdyThermalV5000[numPart]={0};
  Float_t dNdyThermalV6000[numPart]={0};
  Float_t dNdyThermalT165[numPart]={0};
  Float_t dNdyThermalT175[numPart]={0};
  Float_t dNdyThermalSmall[numPart]={0};

  //gSystem->Load("libMathMore");
  //  cout << " defining bessel function " << endl;
  //  cout << ROOT::Math::cyl_bessel_k(2, 1) << endl;
  //TF1 * bessel = new TF1 ("bessel", "ROOT::Math::cyl_bessel_k([0], x)", 0, 10);
  //TF1 * bessel = new TF1 ("bessel", "ROOT::Math::cyl_bessel_k(2, x)", 0, 10);
  //  cout << " setting param" << endl;
  //  bessel->SetParameter(0,2);

 
  for (Int_t part=0; part<numPart; part++){
    cout << "\n************+\n Analyzing particle: " << Part[part] << endl;
    //    cout << "dN/dy thermal for part " <<Part[part]<< endl;
    g[part] = 2*spin[part]+1;
    if (spin[part] ==0 || spin[part]==1) IsFermion[part]=0;
    else  IsFermion[part]=1;
    IsFermion[part]=1;
    //cout << "initial yield" <<  dNdyThermal[part]<< endl;
    for (Int_t par=0; par <numParVariations; par++){
      V = VVar[par];
      T = TVar[par];
      dNdyThermal[part] =0;
      Float_t YieldComponentK =0;
      for (Int_t k=1; k<= HighestK; k++){
	cout << "\n k " << k << endl;  
	YieldComponentK = YieldComponent(k, T, V, 0, mass[part], g[part], IsFermion[part]);
	if (TMath::Abs(YieldComponentK) < 0.001*TMath::Abs(dNdyThermal[part])) break;
	dNdyThermal[part] += YieldComponentK;
	cout << Part[part]<< ": " <<  " k component of yield: " << YieldComponentK << " total yield " << dNdyThermal[part] << "  " <<endl;
	//     if (part==numPart) cout << endl;
      }

      if (par==0)    dNdyThermalDef[part]=  dNdyThermal[part];
      if (par==1)    dNdyThermalV4000[part]=  dNdyThermal[part];
      if (par==2)    dNdyThermalV5000[part]=  dNdyThermal[part];
      if (par==3)    dNdyThermalV6000[part]=  dNdyThermal[part];
      if (par==4)    dNdyThermalT165[part]=  dNdyThermal[part];
      if (par==5)    dNdyThermalT175[part]=  dNdyThermal[part];
    }

    dNdyThermal[part]=  dNdyThermalDef[part];

    if (Part[part] == "#pi^{+}")  {
      dNdyThermalSmall[0] =       dNdyThermal[part];
      gSmall[0] = g[part];
    }
    if (Part[part] == "K^{+}") {
      dNdyThermalSmall[1] =       dNdyThermal[part];  
      gSmall[1] = g[part];
    }
    if (Part[part] == "#phi"){
      dNdyThermalSmall[2] =       dNdyThermal[part];  
      gSmall[2] = g[part];
    }
    if (Part[part] == "p"){
      dNdyThermalSmall[3] =       dNdyThermal[part];  
      gSmall[3] = g[part];
    }
    if (Part[part] == "#Lambda"){ 
      dNdyThermalSmall[4] =       dNdyThermal[part];
      gSmall[4] = g[part];
    }
    if (Part[part] == "#Xi^{+}"){
      dNdyThermalSmall[5] =       dNdyThermal[part]; 
      gSmall[5] = g[part];
    }
    if (Part[part] == "#Omega^{+}"){
      dNdyThermalSmall[6] =       dNdyThermal[part]; 
      gSmall[6] = g[part];
    }
    if (Part[part] == "d")          {
      dNdyThermalSmall[7] =       dNdyThermal[part]; 
      gSmall[7] = g[part];
    }
  }
  
  TLegend * legendParVar = new TLegend(0.7, 0.7, 0.9, 0.9);

  TCanvas *canvas = new TCanvas("canvas", "canvas",800,800);
  canvas->cd();
  gPad->SetLogy();

  TH1F * hdNdyPub = new TH1F("hdNdyPub", "hdNdyPub", numPart, 0.5, 15.5);
  hdNdyPub->SetTitle("Measured yields in 0-10% PbPb collisions at #sqrt{s}=2.76 TeV");
  cout << hdNdyPub->GetNbinsX() << endl;
  for (Int_t b=1; b<=hdNdyPub->GetNbinsX(); b++){
    hdNdyPub->SetBinContent(b, dNdyPub[b-1]);
    hdNdyPub->SetBinError(b, ErrdNdyPub[b-1]);
    hdNdyPub->GetXaxis()->SetBinLabel(b, Part[b-1]);
  }
  //  hdNdyPub->Sumw2();
  hdNdyPub->SetMarkerColor(628);
  hdNdyPub->SetLineColor(628);
  hdNdyPub->SetMarkerStyle(21);
  hdNdyPub->Draw("p");

  TH1F * hdNdyThermal = new TH1F("hdNdyThermal", "hdNdyThermal", numPart, 0.5, 15.5);
  HistoStyle(hdNdyThermal,numPart, dNdyThermal, Part, 862, legendParVar, "V5600 - T156");
  hdNdyThermal->Draw("p same");

  TH1F * hdNdyThermalT165 = new TH1F("hdNdyThermalT165", "hdNdyThermalT165", numPart, 0.5, 15.5);
  HistoStyle(hdNdyThermalT165,numPart, dNdyThermalT165, Part, 797, legendParVar, "V5600 - T165");
  hdNdyThermalT165->Draw("p same");

  TH1F * hdNdyThermalT175 = new TH1F("hdNdyThermalT175", "hdNdyThermalT175", numPart, 0.5, 15.5);
  HistoStyle(hdNdyThermalT175,numPart, dNdyThermalT175, Part, 812, legendParVar, "V5600 - T175");
  hdNdyThermalT175->Draw("p same");

  TH1F * hdNdyThermalV4000 = new TH1F("hdNdyThermalV4000", "hdNdyThermalV400", numPart, 0.5, 15.5);
  HistoStyle(hdNdyThermalV4000,numPart, dNdyThermalV4000, Part, 881, legendParVar, "V4000 - T156");
  hdNdyThermalV4000->Draw("p same");

  TH1F * hdNdyThermalV5000 = new TH1F("hdNdyThermalV5000", "hdNdyThermalV400", numPart, 0.5, 15.5);
  HistoStyle(hdNdyThermalV5000,numPart, dNdyThermalV5000, Part, 909, legendParVar, "V5000 - T156");
  hdNdyThermalV5000->Draw("p same");

  TH1F * hdNdyThermalV6000 = new TH1F("hdNdyThermalV6000", "hdNdyThermalV400", numPart, 0.5, 15.5);
  HistoStyle(hdNdyThermalV6000,numPart, dNdyThermalV6000, Part, 603, legendParVar, "V6000 - T156");
  hdNdyThermalV6000->Draw("p same");

  legendParVar->Draw("");
  /*
  TGraphErrors *GdNdyPub = new TGraphErrors(15, x, dNdyPub, 0, ErrdNdyPub);
  GdNdyPub->SetTitle("Yields in 0-10% PbPb collisions at #sqrt{s}=2.76 TeV");
  GdNdyPub->SetMarkerColor(628);
  GdNdyPub->SetLineColor(628);
  GdNdyPub->SetMarkerStyle(21);
  GdNdyPub->GetYaxis()->SetRangeUser(0.01, 1000);
  GdNdyPub->GetYaxis()->SetTickLength(0);
  */
  //  GdNdyPub->Draw("APsame"); //A seems to be necessary to make the graph appear!

  /*
    TText *t = new TText();
    t->SetTextFont(72);
    t->SetTextSize(0.02);
    t->SetTextAlign(372);
    for (Int_t i=0; i<numPart; i++){
    t->DrawText( x[i],-2,  Part[i]);
    }
  */

  TCanvas *canvasRatio = new TCanvas("canvasRatio", "canvasRatio",800,800);
  canvasRatio->cd();

  TH1F * hdNdyDataModelRatio =  CreateRatio(hdNdyThermal ,hdNdyPub , ErrdNdyPub, "hdNdyDataModelRatio");
  hdNdyDataModelRatio->Draw("p");

  TH1F * hdNdyDataModelRatioV4000 =  CreateRatio(hdNdyThermalV4000 ,hdNdyPub , ErrdNdyPub, "hdNdyDataModelRatioV4000");
  hdNdyDataModelRatioV4000->Draw("same p");

  TH1F * hdNdyDataModelRatioV5000 =  CreateRatio(hdNdyThermalV5000 ,hdNdyPub , ErrdNdyPub, "hdNdyDataModelRatioV5000");
  hdNdyDataModelRatioV5000->Draw("same p");

  TH1F * hdNdyDataModelRatioV6000 =  CreateRatio(hdNdyThermalV6000 ,hdNdyPub , ErrdNdyPub, "hdNdyDataModelRatioV6000");
  hdNdyDataModelRatioV6000->Draw("same p");

  TH1F * hdNdyDataModelRatioT165 =  CreateRatio(hdNdyThermalT165 ,hdNdyPub , ErrdNdyPub, "hdNdyDataModelRatioT165");
  hdNdyDataModelRatioT165->Draw("same p");

  TH1F * hdNdyDataModelRatioT175 =  CreateRatio(hdNdyThermalT175 ,hdNdyPub , ErrdNdyPub, "hdNdyDataModelRatioT175");
  hdNdyDataModelRatioT175->Draw("same p");

  TF1 * lineat1 = new TF1("pol0", "pol0", 0,15.5);
  lineat1->FixParameter(0,1);
  lineat1->SetLineColor(kBlack);
  lineat1->Draw("same");

  legendParVar->Draw("");

  Double_t MassFromPlot[NGraph] = {139.57,  493.677, 1019.46, 938.27, 1115.683, 1321.71,  1672.45, 1865.613};
  Double_t ThermalYieldFromPlot[NGraph]= {190.1977, 51.47391,4.123178, 6.142303,2.4777,0.8189,0.11669,0.03943};
  for(Int_t i=0; i<NGraph; i++) {
    ThermalYieldFromPlot[i] = ThermalYieldFromPlot[i]*gSmall[i];
  }
  Int_t first, second;

  /*  
ifstream fileYieldvsMass;
  fileYieldvsMass.open("DataYieldsvsMassBis.dat");
  Int_t npt=0;
  cout << " here mass and yield extracted with plot digitizer " << endl;
  while(1){
    fileYieldvsMass >> MassFromPlot[npt] >> ThermalYieldFromPlot[npt];
    //    fileYieldvsMass >> first >> second;
    if (!fileYieldvsMass.good()) break;
    cout << " \n" << npt << endl;
    //cout << first << " " << second << endl;
    npt++;
    cout <<  MassFromPlot[npt] << " " <<  ThermalYieldFromPlot[npt]<< endl;;
  }
  fileYieldvsMass.close();
  */
  TCanvas * cYvsMass = new TCanvas ("cYvsMass", "cYvsMass", 800, 500);
  gPad->SetLogy();

  TLegend *  legend = new TLegend (0.7, 0.7, 0.9, 0.9);

  TH1F * DummyHisto = new TH1F("DummyHisto", "ThermalYields", 2000, 0, 2000);
  DummyHisto->SetTitle("Thermal yields");
  DummyHisto->GetYaxis()->SetRangeUser(5e-2, 3e+3);
  DummyHisto->GetXaxis()->SetTitle("mass (MeV/c^{2})");
  DummyHisto->Draw("");

 
  TGraph *GdNdyThermalPlot = new TGraph(NGraph,MassFromPlot,ThermalYieldFromPlot);
  GdNdyThermalPlot->SetMarkerColor(628);
  GdNdyThermalPlot->SetLineColor(628);
  GdNdyThermalPlot->SetMarkerStyle(21);
  legend->AddEntry(  GdNdyThermalPlot, "ThYieldsAndronic", "pl");
  GdNdyThermalPlot->Draw("same p");
  //  GdNdyPub->GetYaxis()->SetRangeUser(0.01, 1000);

  TGraph *GdNdyMyResults = new TGraph(numPart,mass,dNdyThermal);
  GdNdyMyResults->SetMarkerColor(kBlue);
  GdNdyMyResults->SetLineColor(kBlue);
  GdNdyMyResults->SetMarkerStyle(21);
  legend->AddEntry(  GdNdyMyResults, "ThYieldsChiara", "pl");

  TGraph *GdNdyPub = new TGraph(numPart,mass,dNdyPub);
  GdNdyPub->SetMarkerColor(kGreen);
  GdNdyPub->SetLineColor(kGreen);
  GdNdyPub->SetMarkerStyle(21);
  legend->AddEntry(  GdNdyPub, "Measured Yields PbPb 0-10% 2.76TeV", "pl");
  GdNdyPub->Draw("same p");
 
  TLatex *latex[numPart];
  for (Int_t i=0; i <numPart; i++){  
    latex[i] = new TLatex(GdNdyMyResults->GetX()[i], GdNdyMyResults->GetY()[i],Part[i]); 
    if (Part[i] == "K^{0}_{S}") continue;
    latex[i]->SetTextSize(0.07);
    GdNdyMyResults->GetListOfFunctions()->Add(latex[i]); 
  }

  GdNdyMyResults->Draw("p same");


  legend->Draw("");


  TCanvas * cYvsMassRatio = new TCanvas ("cYvsMassRatio", "cYvsMassRatio", 800, 500);

  Double_t ThermalYieldRatio[NGraph]={0};
  for (Int_t i=0; i <NGraph; i++){
    ThermalYieldRatio[i] = dNdyThermalSmall[i]/ThermalYieldFromPlot[i];
    cout << MassFromPlot[i] << " " << dNdyThermalSmall[i] << " " << ThermalYieldFromPlot[i] << " " << ThermalYieldRatio[i] << endl;
  }

  TH1F * DummyHistoBis = new TH1F("DummyHistoBis", "ThermalYields", 2000, 0, 2000);
  DummyHistoBis->SetTitle("Ratio Chiara/Andronic thermal tields");
  DummyHistoBis->GetYaxis()->SetRangeUser(1,1.5);
  DummyHistoBis->GetXaxis()->SetTitle("mass (MeV/c^{2})");
  DummyHistoBis->Draw("");

  TGraph *GdNdyRatio = new TGraph(NGraph,MassFromPlot,ThermalYieldRatio);
  GdNdyRatio->SetMarkerColor(628);
  GdNdyRatio->SetLineColor(628);
  GdNdyRatio->SetMarkerStyle(21);

  TLatex *     latexBis[NGraph];
  for (Int_t i=0; i <NGraph; i++){  
    latexBis[i] = new TLatex(GdNdyRatio->GetX()[i], GdNdyRatio->GetY()[i],PartSmall[i]); 
    latexBis[i]->SetTextSize(0.07);
    GdNdyRatio->GetListOfFunctions()->Add(latexBis[i]); 
  }

  GdNdyRatio->Draw("same p");
}
