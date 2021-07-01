#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TText.h>
#include <TSystem.h>
#include <TKey.h>
#include <THnSparse.h>

#include "inc/constants.h"

using namespace std;

//Functions
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, string skipruns, long int refrun, int IBorOB);
int nStavesInLay[7]= {0};

void CompareNoisyPixelsInRuns(){
  string fpath;
  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  //IBorOB = 0 if I want to check all IB layers                                                              
  //IBorOB = 1 if I want to check all OB layers                                                              
  //IBorOB = 2 if I want to check all IB + OB layers or if I want to check a single layer                    

  int IBorOB=0;
  if(fpath.find("IB")!=string::npos){
    IBorOB = 0;
    nStavesInLay[0] = nStavesInLayAll[0];
    nStavesInLay[1] = nStavesInLayAll[1];
    nStavesInLay[2] = nStavesInLayAll[2];
  }
  else if (fpath.find("OB")!=string::npos){
    IBorOB = 1;
    nStavesInLay[0] = nStavesInLayAll[3];
    nStavesInLay[1] = nStavesInLayAll[4];
    nStavesInLay[2] = nStavesInLayAll[5];
    nStavesInLay[3] = nStavesInLayAll[6];
  }
  else if (fpath.find("all")!=string::npos){
    IBorOB = 2;
    for (Int_t i =0; i < 7; i ++){
      nStavesInLay[i] = nStavesInLayAll[i];
    }
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    IBorOB = 2;
    nStavesInLay[0] = nStavesInLayAll[stoi(layernum)];
  }

  cout << "nStavesInLay[0] " << nStavesInLay[0] << endl;

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

  cout<<"Available runs in your file:\n"<<endl;
  TFile *infile=new TFile(fpath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  vector<string> stavenums;
  while((key = ((TKey*)next()))){
    obj=key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
         && (!obj->InheritsFrom("THnSparse"))) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user
    if(stvnum.find("_")!=string::npos)
      stvnum = objname.substr(objname.find("Stv")+3,1);
    stavenums.push_back(stvnum);
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
    if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2])
      break;

    cout<<runnum<<endl;
  }

  long int refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;

  DoAnalysis(fpath, nchips, skipruns, refrun, IBorOB);
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, string skipruns, long int refrun, int IBorOB){

  gStyle->SetOptStat(0000);

  std::vector<THnSparse*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  vector<int> posrefrun;
  int nLayers=1, nRuns=1;

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  THnSparse *hsparse;
  while((key=((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    hsparse = (THnSparse*)obj;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(hsparse);

    if(stol(runnum)==refrun) //position of refence run
      posrefrun.push_back((int)hmaps.size()-1);

    if((int)laynums.size()>1 && laynum!=laynums[laynums.size()-1]){
      nLayers++;
    }

    if((int)stavenums.size()>1 && stvnum==stavenums[stavenums.size()-1])
      nRuns++;
    else nRuns=1;

    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);
    stavenums.push_back(stvnum);
  }

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  vector<string> runlabel;
  int istave = 0;
  for(int ilay=0; ilay<nLayers; ilay++){
    istave+=nStavesInLay[ilay];
  }
  istave--;
  int ilayer = nLayers-1;
  long int first[nLayers][nRuns], second[nLayers][nRuns], both[nLayers][nRuns];
  int irun=0;
  vector<array<long int,5>> noisypix;
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int i=0; i<nRuns; i++){
      first[ilay][i]=0; second[ilay][i]=0; both[ilay][i]=0;
    }
  cout << "hmaps size " << (int)hmaps.size() << endl;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    cout <<"ihist " << ihist << endl;
    if(runnumbers[ihist].find(std::to_string(refrun))!=string::npos){
      if(ihist>0){// in case ref run is the first into the list of runs
        if(stavenums[ihist-1]!=stavenums[ihist]){
          istave--;
          irun=0;
        }
        if(laynums[ihist-1]!=laynums[ihist]){
          ilayer--;
        }
      }
      continue;
    }
    //if(!hmaps[ihist]->GetEntries()) continue; // do not compare ref run with empty run (= empty maps)

    noisypix.push_back(CompareTwoRuns(hmaps[posrefrun[istave]], hmaps[ihist]));
    noisypix[noisypix.size()-1][3] = stol(stavenums[ihist]);
    noisypix[noisypix.size()-1][4] = stol(laynums[ihist]);
    first[ilayer][irun]+=noisypix[noisypix.size()-1][0];
    second[ilayer][irun]+=noisypix[noisypix.size()-1][1];
    both[ilayer][irun]+=noisypix[noisypix.size()-1][2];
    irun++;

    //cout<<noisypix[noisypix.size()-1][0]<<"  "<<noisypix[noisypix.size()-1][1]<<"  "<<noisypix[noisypix.size()-1][2]<<"  "<<noisypix[noisypix.size()-1][3]<<endl;
    if(ilayer==nLayers-1 && stavenums[ihist]=="0")
      runlabel.push_back(runnumbers[ihist]);

    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
        istave--;
        irun=0;
      }
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
    }
  }//end loop on histograms

  //Make plot for each layer and for each stave in the root file
  TGraphErrors *ge_nref[nLayers];
  TGraphErrors *ge_n2[nLayers];
  TGraphErrors *ge_ncom1[nLayers];
  TGraphErrors *ge_ncom2[nLayers];
  TGraphErrors *ge_nref_stave[nLayers][100];//100 is an abitrary large number of staves
  TGraphErrors *ge_n2_stave[nLayers][100];
  TGraphErrors *ge_ncom1_stave[nLayers][100];
  TGraphErrors *ge_ncom2_stave[nLayers][100];
  double xshift = 3.;
  double max[nLayers];
  double min[nLayers];

  cout << "making plots " << endl;
  for(int ilay=0; ilay<nLayers; ilay++){
    ge_nref[ilay] = new TGraphErrors();
    ge_n2[ilay] = new TGraphErrors();
    ge_ncom1[ilay] = new TGraphErrors();
    ge_ncom2[ilay] = new TGraphErrors();
    max[ilay] = -1.;
    min[ilay] = 1e35;

    for(int ir=0; ir<nRuns-1; ir++){//first the older data and last the most recent
      //first couple of bar on the left
      //int ipoint = (int)noisypix.size()-icomp-1;
      if(!ir) xshift=1.;
      else xshift = 3.;
      ge_nref[ilay]->SetPoint(ir, ir*xshift, (double)both[ilay][ir]/2.+(double)first[ilay][ir]/2.);
      ge_nref[ilay]->SetPointError(ir, 0.5, (double)first[ilay][ir]/2.);
      ge_ncom1[ilay]->SetPoint(ir, ir*xshift, 0.);
      ge_ncom1[ilay]->SetPointError(ir, 0.5, (double)both[ilay][ir]/2.);
      if((double)both[ilay][ir]/2.+(double)first[ilay][ir] > max[ilay]) max[ilay] = (double)both[ilay][ir]/2.+(double)first[ilay][ir];

      //second couple of bar on the right
      ge_n2[ilay]->SetPoint(ir, ir*xshift+1, -(double)both[ilay][ir]/2.-(double)second[ilay][ir]/2.);
      ge_n2[ilay]->SetPointError(ir, 0.5, (double)second[ilay][ir]/2.);
      ge_ncom2[ilay]->SetPoint(ir, ir*xshift+1, 0.);
      ge_ncom2[ilay]->SetPointError(ir, 0.5, (double)both[ilay][ir]/2.);
      if(-(double)both[ilay][ir]/2.-(double)second[ilay][ir] < min[ilay]) min[ilay] = -(double)both[ilay][ir]/2.-(double)second[ilay][ir];
    }//end first loop on runs

    //Style
    SetStyle(ge_nref[ilay], kBlue);
    SetStyle(ge_ncom1[ilay], kBlack);
    SetStyle(ge_ncom2[ilay], kBlack);
    SetStyle(ge_n2[ilay], kRed+2);
  }//end loop on layers

  //FILL PLOTS FOR EACH STAVE
  //fill the plots for each stave
  int cnt = 0;
  double maxs[nLayers][100];
  double mins[nLayers][100];
  for(int ilay=nLayers-1; ilay>=0; ilay--){
    cout << "\n nStavesInLay[ilay]-1 = " << nStavesInLay[ilay]-1 << endl;
    for(int is=nStavesInLay[ilay]-1; is>=0; is--){
      maxs[ilay][is] = -1.;
      mins[ilay][is] = 1e35;
      ge_nref_stave[ilay][is] = new TGraphErrors();
      ge_ncom1_stave[ilay][is] = new TGraphErrors();
      ge_n2_stave[ilay][is] = new TGraphErrors();
      ge_ncom2_stave[ilay][is] = new TGraphErrors();
      cout << "nRuns " << nRuns << endl;
      for(int ir=0; ir<nRuns-1; ir++){//first the older data and last the most recent
	cout << "ir " << ir << endl;
        if(!ir) xshift=1.;
        else xshift = 3.;
        ge_nref_stave[ilay][is]->SetPoint(ir, ir*xshift, (double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0]/2.);
        ge_nref_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][0]/2.);
        ge_ncom1_stave[ilay][is]->SetPoint(ir, ir*xshift, 0.);
        ge_ncom1_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
        if((double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0] > maxs[ilay][is]) maxs[ilay][is] = (double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0];

        //second couple of bar on the right
        ge_n2_stave[ilay][is]->SetPoint(ir, ir*xshift+1, -(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1]/2.);
        ge_n2_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][1]/2.);
        ge_ncom2_stave[ilay][is]->SetPoint(ir, ir*xshift+1, 0.);
        ge_ncom2_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
        if(-(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1] < mins[ilay][is]) mins[ilay][is] = -(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1];
        cnt++;
      }
      SetStyle(ge_nref_stave[ilay][is], kBlue);
      SetStyle(ge_ncom1_stave[ilay][is], kBlack);
      SetStyle(ge_ncom2_stave[ilay][is], kBlack);
      SetStyle(ge_n2_stave[ilay][is], kRed+2);
    }
  }

  //Legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref[0], "#splitline{#noisy pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2[0], "#splitline{#noisy pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1[0], "#splitline{#noisy pix}{both}");

  //Draw plot for each layer
  Int_t indexEff = 0;
  Int_t indexLaynums = 0;
  for(int ilay=0; ilay<nLayers; ilay++){
    if (ilay==0) indexLaynums = 0;
    else {
      indexEff += nStavesInLay[ilay-1];
      indexLaynums = nRuns*indexEff;
    }
    TCanvas *canvas = new TCanvas(Form("mycanvas_%d",ilay), Form("mycanvas_%d",ilay), 1300, 800);
    canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
    canvas->cd();

    //fake histo (just for the axes)
    double x2,y2;
    cout << "should be the number of run or sth like that " << ge_ncom2[ilay]->GetN() << endl;
    ge_ncom2[ilay]->GetPoint(ge_ncom2[ilay]->GetN()-1, x2,y2);
    TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
    //draw labels on x axis
    int counter = 0;
    for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
      hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
      counter++;
    }
    hfake->Draw();
    //canvas->SetLogy();
    cout << "fake region " << endl;
    cout << indexLaynums<< endl;
    cout << laynums[indexLaynums] << endl;
    //    hfake->SetTitle(Form("Layer-%s - %s%06ld compared to all",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(), filepath.find("run")==string::npos? "":"run",refrun));
    hfake->SetTitle(Form("Layer-%s - %s%06ld compared to all",laynums[indexLaynums].c_str(), filepath.find("run")==string::npos? "":"run",refrun));
    ge_nref[ilay]->Draw("P E2 same");
    ge_ncom1[ilay]->Draw("E2 same");
    ge_ncom2[ilay]->Draw("E2 same");
    ge_n2[ilay]->Draw("E2 same");
    hfake->GetYaxis()->SetRangeUser(min[ilay]+0.1*min[ilay], max[ilay]+0.1*max[ilay]);
    //hfake->GetYaxis()->SetLabelColor(kWhite);
    hfake->GetYaxis()->SetTickLength(0.005);
    hfake->GetYaxis()->SetMaxDigits(4);
    TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
    lineref->SetLineColor(kGray-1);
    lineref->SetLineStyle(2);
    lineref->Draw("same");

    //draw legend
    leg->Draw("same");

    canvas->SaveAs(Form("../Plots/Layer%s_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf", laynums[indexLaynums].c_str(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_NoisyPixComparison_%s%ld_compared_to_run_%s.root", laynums[indexLaynums].c_str(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

    delete canvas;
    delete hfake;
    delete lineref;
  }//end loop on layers

  //Draw plot for each stave
  //Draw
  indexEff = 0;
  indexLaynums = 0;
  for(int ilay=0; ilay<nLayers; ilay++){
    if (ilay==0) indexLaynums = 0;
    else {
      indexEff += nStavesInLay[ilay-1];
      indexLaynums = nRuns*indexEff;
    }
    for(int is=0; is<nStavesInLay[ilay]; is++){
      TCanvas *canvas = new TCanvas(Form("mycanvas_%d_%d",ilay,is), Form("mycanvas_%d_%d",ilay,is), 1300, 800);
      canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
      canvas->cd();

      //fake histo (just for the axes)
      double x2,y2;
      ge_ncom2_stave[ilay][is]->GetPoint(ge_ncom2_stave[ilay][is]->GetN()-1, x2,y2);
      TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
      //draw labels on x axis
      int counter = 0;
      for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
        hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
        counter++;
      }
      hfake->Draw();
      //canvas->SetLogy();
      hfake->SetTitle(Form("Layer-%s - Stave-%d - %s%06ld compared to all",laynums[indexLaynums].c_str(), is, filepath.find("run")==string::npos? "":"run",refrun));
      ge_nref_stave[ilay][is]->Draw("P E2 same");
      ge_ncom1_stave[ilay][is]->Draw("E2 same");
      ge_ncom2_stave[ilay][is]->Draw("E2 same");
      ge_n2_stave[ilay][is]->Draw("E2 same");
      hfake->GetYaxis()->SetRangeUser(mins[ilay][is]+0.1*mins[ilay][is], maxs[ilay][is]+0.1*maxs[ilay][is]);
      //hfake->GetYaxis()->SetLabelColor(kWhite);
      hfake->GetYaxis()->SetTickLength(0.005);
      hfake->GetYaxis()->SetMaxDigits(4);
      TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
      lineref->SetLineColor(kGray-1);
      lineref->SetLineStyle(2);
      lineref->Draw("same");
      leg->Draw("same");

      TString LayerTitle="";
      if (IBorOB==0) LayerTitle = "AllLayersIB_";
      else if (IBorOB==1) LayerTitle = "AllLayersOB_";
      else {
	if (nLayers==1) LayerTitle = Form("Layer%s_",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str());
	else LayerTitle = "AllLayers_";
      }
      //      nLayers==1 ? Form("Layer%s_",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str()) : "AllLayers_"
      if(!ilay && !is) canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf[", LayerTitle.Data(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf", LayerTitle.Data() ,filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      if(ilay==nLayers-1 && is==nStavesInLay[ilay]-1) canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf]", LayerTitle.Data(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

      delete canvas;
      delete hfake;
      delete lineref;
    }
  }

}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2){

  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int ibinref=0; ibinref<href->GetNbins(); ibinref++){
    int coordref[2];
    double bincref = href->GetBinContent(ibinref, coordref);
    bool isfound = false;
    for(int ibin2=0; ibin2<h2->GetNbins(); ibin2++){
      int coord2[2];
      double binc2 = h2->GetBinContent(ibin2,coord2);

      if((coordref[0]==coord2[0] && coordref[1]==coord2[1]) && (bincref>2 && binc2>2)){//noisy in both runs
        noisypix[2]++;
        isfound=true;
        break;
      }

    }
    if(!isfound && bincref>2)
      noisypix[0]++; //noisy only in ref run
  }

  for(int ibin2=0; ibin2<h2->GetNbins(); ibin2++){
    int coord2[2];
    double binc2 = h2->GetBinContent(ibin2, coord2);
    bool isfound = false;
    for(int ibinref=0; ibinref<href->GetNbins(); ibinref++){
      int coordref[2];
      double bincref = href->GetBinContent(ibinref,coordref);

      if((coordref[0]==coord2[0] && coordref[1]==coord2[1]) && (bincref>2 && binc2>2)){//noisy in both runs
        isfound=true;
        break;
      }

    }
    if(!isfound && binc2>2)
      noisypix[1]++; //noisy only in second run
  }

  /*cout<<"Stave/run: "<<h2->GetName()<<endl;
  cout<<"Ref run:  "<<noisypix[0]<<endl;
  cout<<"Sec run:  "<<noisypix[1]<<endl;
  cout<<"Both run: "<<noisypix[2]<<endl;*/

  return noisypix;
}

//
// Style
//
void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}
