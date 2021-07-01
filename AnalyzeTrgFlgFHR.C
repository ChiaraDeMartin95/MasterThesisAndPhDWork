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
void DoAnalysis(string filepath, const int nChips, string skipruns);
TString SIBorOB[2]={"IB", "OB"};

//
// MAIN
//
void AnalyzeTrgFlgFHR(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*w_error_and_trig* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

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
  DoAnalysis(fpath, nchips, skipruns);

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
void DoAnalysis(string filepath, const int nChips, string skipruns){

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
    if(objname.find("trg")==string::npos) continue;
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

  //sum all the histos in a single histogram (for summary plot) for each layer
  TH2D *hSummary = (TH2D*)herr[0]->Clone("hSummary");
  for(int iplot=1; iplot<(int)herr.size(); iplot++){
    hSummary->Add(herr[iplot]);
  }

  //Make plots with Trg IDs vs Run for each layer
  TGraph *trend[hSummary->GetNbinsY()][2];

  int ibinMin=1;
  int ibinMax=1;
  int IBorOBindex = 0;
  int IBorOB = 0;
  if (hSummary->GetNbinsX() == 144) {  //for backward compatibility (plots for IB only)                        
    ibinMax = hSummary->GetNbinsX();
    IBorOB = 0;
    IBorOBindex = 0;
  }
  else if (hSummary->GetNbinsX() == 288) { //for backward compatibility (plots for OB only)                    
    IBorOB = 1;
    IBorOBindex = 1;
    ibinMax = hSummary->GetNbinsX();
  }
  else IBorOB=2; //plots for IB + OB     

  int ir = 0;
  double max[2] = {-1., -1.};
  TH1D *hproj;
  for(int iplot=0; iplot<(int)herr.size(); iplot++){
    for (Int_t i=0; i<=1; i++){
      if (IBorOB==2){
        IBorOBindex = i;
	if (i==0) {ibinMin =1; ibinMax = 144;} //IB                                                        
	else  {ibinMin =144; ibinMax = herr[iplot]->GetNbinsX();} //OB                               
      }   
      hproj = (TH1D*)herr[iplot]->ProjectionY(Form("herr_%d",iplot), ibinMin, ibinMax);
      for(int ibin=1; ibin<=hSummary->GetNbinsY(); ibin++){
	if(ir==0){
          trend[ibin-1][i] = new TGraph();
          trend[ibin-1][i]->SetName(Form("gr_trgID%d_%s",ibin, SIBorOB[IBorOBindex].Data()));
	  SetStyle(trend[ibin-1][i], col[ibin<=10?ibin-1:ibin-11], ibin<=10?24:25);
        }
        trend[ibin-1][i]->SetPoint(ir,ir, hproj->GetBinContent(ibin));
        if(hproj->GetBinContent(ibin)>max[i])
          max[i]=hproj->GetBinContent(ibin);
      }
    }
    delete hproj;
    ir++;
  }

  if (IBorOB==2) IBorOBindex=0;

  //Draw summary plot
  TCanvas canvas;
  canvas.cd();
  canvas.SetTickx();
  canvas.SetTicky();
  canvas.SetLogz();
  canvas.SetMargin(0.18,0.2,0.194,0.0993);
  canvas.SetRightMargin(0.15);
  hSummary->SetTitle(Form("Trigger & Flags, %s",filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  hSummary->Draw("colz");
  //hSummary[ilay]->GetXaxis()->SetNdivisions(530);
  //hSummary[ilay]->GetYaxis()->SetNdivisions(516);
  hSummary->GetXaxis()->SetLabelSize(0.045);
  hSummary->GetYaxis()->SetLabelSize(0.045);
  hSummary->GetYaxis()->SetTitleOffset(1.3);
  hSummary->GetZaxis()->SetLabelSize(0.045);
  hSummary->GetXaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleOffset(0.7);
  hSummary->GetZaxis()->SetTitleSize(0.05);
  hSummary->GetZaxis()->SetTitleOffset(0.9);

  canvas.SaveAs(Form("../Plots/FHRTrgFlgPlotSummary_%s.pdf[", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas.SaveAs(Form("../Plots/FHRTrgFlgPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

  //Draw trends
  int npoints = trend[0][0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; # Errors", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<(int)runnumbers.size(); ir++)
      hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[runnumbers.size()-1-ir])));

  TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
  leg->SetHeader("IDs");
  for(int iid=1; iid<=hSummary->GetNbinsY();iid++)
    leg->AddEntry(trend[iid-1][0], Form("%d", iid), "p");

  TCanvas canvas2;
  canvas2.cd();
  canvas2.SetTickx();
  canvas2.SetTicky();
  canvas2.SetLogy();
  canvas2.SetMargin(0.0988,0.1,0.194,0.0993);

  hfake->GetXaxis()->SetTitleOffset(2.8);
  hfake->SetTitle(SIBorOB[IBorOBindex] + Form(", Trigger & Flag trends %s", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  hfake->GetYaxis()->SetRangeUser(1, 10*max[0]);
  hfake->GetXaxis()->SetTitleOffset(2.8);
  hfake->Draw();
  for(int iid=1; iid<=hSummary->GetNbinsY();iid++){
    trend[iid-1][0]->Draw("P same");
  }
  leg->Draw("same");

  canvas2.SaveAs(Form("../Plots/FHRTrgFlgPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  if (IBorOB!=2)    canvas2.SaveAs(Form("../Plots/FHRTrgFlgPlotSummary_%s.pdf]", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

  TCanvas canvas3;
  canvas3.cd();
  canvas3.SetTickx();
  canvas3.SetTicky();
  canvas3.SetLogy();
  canvas3.SetMargin(0.0988,0.1,0.194,0.0993);

  hfake->SetTitle(Form("OB, Trigger & Flag trends %s", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  hfake->GetYaxis()->SetRangeUser(1, 10*max[1]);
  hfake->Draw();
  for(int iid=1; iid<=hSummary->GetNbinsY();iid++){
    trend[iid-1][1]->Draw("P same");
  }
  leg->Draw("same");

  if (IBorOB ==2){
    canvas3.SaveAs(Form("../Plots/FHRTrgFlgPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas3.SaveAs(Form("../Plots/FHRTrgFlgPlotSummary_%s.pdf]", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

}
