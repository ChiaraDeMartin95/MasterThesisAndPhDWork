#include <string>
#include <iostream>
#include <vector>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TKey.h>

using namespace std;

void SetStyle(TGraph *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath, const int nChips, string skipruns, int IBorOB);
TString SIBorOB[2]={"IB", "OB"};

//
// MAIN
//
void AnalyzeErrorsFHR(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*w_error_and_trig* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  int IBorOB;
  int IBorOBmin=0;
  int IBorOBmax=0;
  //IBorOB = 0 if I want to check all IB layers                                                              
  //IBorOB = 1 if I want to check all OB layers                                                              
  //IBorOB = 2 if I want to check all IB + OB layers or if I want to check a single layer                    

  if(fpath.find("IB")!=string::npos){
    IBorOB = 0;
    // IBorOBmin=0;
    // IBorOBmax=0;
  }
  else if (fpath.find("OB")!=string::npos){
    IBorOB = 1;
    // IBorOBmin=1;
    // IBorOBmax=1;
  }
  else if (fpath.find("all")!=string::npos){
    IBorOB = 2;
    // IBorOBmin=0;
    // IBorOBmax=1;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    IBorOB = 2;
    // IBorOBmin=0;
    // IBorOBmax=0;
  }

  //Choose whether to skip runs
  string skipans, skipruns;
  cout<<endl;
  cout<<"Would you like to skip some run(s)? [y/n] ";
  cin>>skipans;
  if(skipans=="y" || skipans=="Y"){
    cout<<endl;
    cout<<"Specify run number(s) separated by comma (no white spaces!):";
    cin>>skipruns;
    cout<<endl;
  }
  else
    skipruns=" ";


  //Call
  DoAnalysis(fpath, nchips, skipruns, IBorOB);

}

//
//Set Style
//
void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, string skipruns, int IBorOB){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> herr;
  std::vector<string> timestamps, runnumbers;
  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  TKey *key;
  TObject *obj;
  TIter next(list);
  TH2 *h2;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("err")==string::npos) continue;
    h2 = (TH2*)obj;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    herr.push_back(h2);
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
  }

  int nRuns = (int)runnumbers.size();

  /*
    Int_t iplotmin=1;
    Int_t iplotmax=(int)herr.size();
    Int_t indexSOBorOB=0;
  for (Int_t IBorOBindex = IBorOBmin; IBorOBindex<= IBorOBmax; IBorOBindex++)
  if (IBorOB==2){ //exclude case when only 1 layer is analysed
  if (IBorOBindex==0) {
  iplotmin =1; 
  iplotmax = (int)herr.size()/2;
  }
  else {
  iplotmin = (int)herr.size()/2; 
  iplotmax = (int)herr.size();
  }
  }
  SindexIBorOB = IBorOBindex; --> if analysing just one layer: if (layernum < =2) = 0, else =1
  */

  //sum all the histos in a single histogram (for summary plot) for each layer
  TH2D *hSummary =  (TH2D*)herr[0]->Clone("hSummary");
  for(int iplot=1; iplot<(int)herr.size(); iplot++){
    hSummary->Add(herr[iplot]);
  }

  //Make plots with Error IDs vs Run for each layer
  TGraph *trend[hSummary->GetNbinsY()][2];

  int ir = 0;
  double max = -1.;
  for(int iplot=0; iplot<(int)herr.size(); iplot++){
    for (Int_t i=0; i<=1; i++){
      if (i==0) {ibinMin =1; ibinMax = 144;} //IB
      else  {ibinMin =144; ibinMax = herr[iplot]->GetNbinsX();} //OB
      TH1D *hproj = (TH1D*)herr[iplot]->ProjectionY(Form("herr_%d",iplot), ibinMin, ibinMax);
      for(int ibin=1; ibin<=hproj->GetNbinsX(); ibin++){
	if(ir==0){
	  trend[ibin-1][i] = new TGraph();
	  trend[ibin-1][i]->SetName(Form("gr_errID%d_%s",ibin, SIBorOB[i].c_str()));
	  SetStyle(trend[ibin-1][i], col[ibin<=10?ibin-1:ibin<=20?ibin-11:ibin-21], ibin<=10?24:ibin<=20?25:26);
	}
	trend[ibin-1][i]->SetPoint(ir,ir, hproj->GetBinContent(ibin));
	if(hproj->GetBinContent(ibin)>max)
	  max=hproj->GetBinContent(ibin);
      }
    }
    delete hproj;
    ir++;
  }

  //Draw summary plot
  TCanvas canvas;
  canvas.cd();
  canvas.SetTickx();
  canvas.SetTicky();
  canvas.SetLogz();
  canvas.SetMargin(0.0988,0.2,0.194,0.0993);
  canvas.SetRightMargin(0.15);
  hSummary->SetTitle(Form("Errors , %s", filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  hSummary->Draw("colz");
  //hSummary[ilay]->GetXaxis()->SetNdivisions(530);
  //hSummary[ilay]->GetYaxis()->SetNdivisions(516);
  hSummary->GetXaxis()->SetLabelSize(0.045);
  hSummary->GetYaxis()->SetLabelSize(0.045);
  hSummary->GetZaxis()->SetLabelSize(0.045);
  hSummary->GetXaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleOffset(0.7);
  hSummary->GetZaxis()->SetTitleSize(0.05);
  hSummary->GetZaxis()->SetTitleOffset(0.9);

  canvas.SaveAs(Form("../Plots/FHRErrorPlotSummary_%s.pdf[", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas.SaveAs(Form("../Plots/FHRErrorPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

  //Draw trends
  int npoints = trend[0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; # Errors", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<(int)runnumbers.size(); ir++)
    hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[runnumbers.size()-1-ir])));

  TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
  leg->SetHeader("Error IDs");
  leg->SetNColumns(2);
  for(int iid=1; iid<=hSummary->GetNbinsY();iid++)
    leg->AddEntry(trend[iid-1], Form("%d",iid), "p");

  TCanvas canvas2;
  canvas2.cd();
  canvas2.SetTickx();
  canvas2.SetTicky();
  canvas2.SetLogy();
  canvas2.SetMargin(0.0988,0.1,0.194,0.0993);

  hfake->GetXaxis()->SetTitleOffset(2.8);
  hfake->SetTitle(Form("IB, Error trends %s", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  hfake->GetYaxis()->SetRangeUser(1, 10*max[0]);
  hfake->GetXaxis()->SetTitleOffset(2.8);
  hfake->DrawClone();
  for(int iid=1; iid<=hSummary->GetNbinsY();iid++){
    trend[iid-1][0]->Draw("P same");
  }
  leg->Draw("same");

  canvas2.SaveAs(Form("../Plots/FHRErrorPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

  TCanvas canvas3 = (TCanvas)canvas2->Clone();
  hfake->SetTitle(Form("OB, Error trends %s", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  hfake->GetYaxis()->SetRangeUser(1, 10*max[1]);
  hfake->Draw();
  for(int iid=1; iid<=hSummary->GetNbinsY();iid++){
    trend[iid-1][1]->Draw("P same");
  }
  leg->Draw("same");

  canvas3.SaveAs(Form("../Plots/FHRErrorPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas3.SaveAs(Form("../Plots/FHRErrorPlotSummary_%s.pdf]", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

}
