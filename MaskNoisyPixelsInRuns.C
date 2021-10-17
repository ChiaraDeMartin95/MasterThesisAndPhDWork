#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TH3.h>
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
#include <TLatex.h>
#include <THnSparse.h>
#include <TStopwatch.h>
#include <TClass.h> //new
//#include "inc/constants.h"

using namespace std;
const int nMasked = 10; //5000
const int nStavesInLayAll[7] = {12, 16, 20, 24, 30, 42, 48};

//Functions
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap, bool HalfStave, bool IB,  bool isHotPixelMapDrawn);
int GetNchipsActive(THnSparse *hmap, int maxchip, int MaxRange, bool HalfStave, bool IB);
int GetNrunsWOhits(TH2 *hFhrStv);
void SetStyle(TH1 *h, int col, Style_t mkr);
void DoAnalysis(string filepath_hit, string skipruns, int IBorOB, bool isHotPixelMapDrawn);

void MaskNoisyPixelsInRuns(){
  string fpath;

  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;
  int IBorOB;

  //IBorOB = 0 if I want to check all IB layers
  //IBorOB = 1 if I want to check all OB layers
  //IBorOB = 2 if I want to check all IB + OB layers or if I want to check a single layer

  if(fpath.find("IB")!=string::npos){
    IBorOB = 0;
  }
  else if (fpath.find("OB")!=string::npos){
    IBorOB = 1;
  }
  else if (fpath.find("all")!=string::npos){
    IBorOB = 2;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    IBorOB = 2;
  }

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

  string drawpixelmap;
  bool isHotPixelMapDrawn =0;
  cout << "Would you like to save the hot pixel map? [y/n] ";
  cin >> drawpixelmap;
  if (drawpixelmap=="y" || drawpixelmap=="Y"){
    cout << endl;
    isHotPixelMapDrawn =1;
  }

  DoAnalysis(fpath, skipruns, IBorOB, isHotPixelMapDrawn);
}

//
// Analysis
//
void DoAnalysis(string filepath_hit, string skipruns, int IBorOB, bool isHotPixelMapDrawn){

  //TStopwatch t;
  //TStopwatch t1;
  //t.Start();
  //t1.Start();
  gStyle->SetOptStat(0000);
  int col[] = {810, 807, 797, 827, 417, 841, 868, 867, 860, 602, 921, 874};

  std::vector<THnSparse*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums, laynumsBis;
  vector<int> posrefrun;
  int nLayersInput=1;
  int nTimes=0;
  int nRunsB[7]={-1};
  for (int ilay=0; ilay < 7; ilay++){
    nRunsB[ilay] =-1;
  }

  //Read the file and the list of plots with entries (hitmaps!)
  TFile *infile=new TFile(filepath_hit.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  while((key=(TKey*)next())){
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
    THnSparse *h2 = (THnSparse*)obj;
    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }
    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    hmaps.push_back(h2);

    if((int)laynums.size()>1 && laynum!=laynums[laynums.size()-1]){
      nLayersInput++;
    }

    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);
    stavenums.push_back(stvnum);
  }

  //Read file with fhr maps for each layer --> ONLY TO EXTRACT NUMBER OF TRIGGERS!!
  TList *list2 = (TList*)infile->GetListOfKeys();
  TIter next2(list2);
  TH2 *h2_2[100][7];
  TKey *key2;
  TObject *obj2;
  int cntrun=0;
  while((key2=((TKey*)next2()))){
    obj2 = key2->ReadObj();
    if ((strcmp(obj2->IsA()->GetName(),"TProfile")!=0)
	&& (!obj2->InheritsFrom("TH2"))
	&& (!obj2->InheritsFrom("TH1")
	    && (!obj2->InheritsFrom("THnSparse")))
	) {
      cout<<"<W> Object "<<obj2->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
    }
    string objname2 = (string)obj2->GetName();
    if(objname2.find("Stv")!=string::npos) break;
    string laynum = objname2.substr(objname2.find("L")+1,1);
    string runnum =  objname2.find("run")==string::npos ? "norun":objname2.substr(objname2.find("run")+3, 6);
    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    laynumsBis.push_back(laynum);
    nTimes++;
    if (nRunsB[stoi(laynum)]==-1) nRunsB[stoi(laynum)]=0;
    if(nTimes>1){
      if (laynum==laynumsBis[laynumsBis.size()-2]) {
	nRunsB[stoi(laynum)]++;
      }
    }
    h2_2[nRunsB[stoi(laynum)]][stoi(laynum)] = (TH2*)obj2;
    cout<<"... Reading "<<obj2->GetName()<< " Run: " << nRunsB[stoi(laynum)] << " layer: " << laynum << endl;
  }
  
  if (nTimes ==0) {cout << "The input file contains no plots\n"; return; }

  int nLayers;
  nLayers=nLayersInput;

  int ilayEff=0;
  for (int ilay=0; ilay < nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + stoi(laynums[0]) ;
    else if (IBorOB==2) ilayEff = ilay + stoi(laynums[0]) ;
    else ilayEff = ilay;
    if (nRunsB[ilayEff] ==-1) continue;
    cout << "L" << ilayEff << " nRuns: " << nRunsB[ilayEff] << endl;
  }

  //Calculate number of triggers for each run
  cout<<endl;
  cout<<"... Extract the number of triggers for each run (for each layer)"<<endl;
  int MaxRange  = 0; 
  int MaxRangeY = 0;
  int NChipsPerHIC   = 1;  
  int ChipRowsPerHIC = 1;
  int HICsPerHS = 0; //Number of HICs in each half stave
  int HalfStaveFound = 0;
  vector<double> ntrig;
  int nEntrieshmaps =0; 
  int nRunsTotAllLayers =0; 

  cout << "nLayers " << nLayers << endl;
  for (int ilay=0; ilay<nLayers; ilay++){
    cout << "ilay " << ilay << endl;
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + stoi(laynums[0]) ;
    else if (IBorOB==2) ilayEff = ilay + stoi(laynums[0]) ;
    else ilayEff = ilay;
    if (nRunsB[ilayEff] ==-1) continue;
    cout << "ilay " << ilayEff << endl;
    if (ilay > 0) nEntrieshmaps += (nRunsB[ilayEff-1]+1) * nStavesInLayAll[ilayEff-1];
    nRunsTotAllLayers += (nRunsB[ilayEff]+1);

    for(int ir=0; ir<=nRunsB[ilayEff]; ir++){
      double fhr_run = h2_2[ir][ilayEff]->GetBinContent(1,1);
      int stavefound = 0;
      int chipfound = 0;
      int chipfoundeff = 0;
      if(fhr_run<1e-15){
	for(int ibinx=1; ibinx<=h2_2[ir][ilayEff]->GetNbinsX(); ibinx++){
	  for(int ibiny=h2_2[ir][ilayEff]->GetNbinsY(); ibiny>=1; ibiny--){
	    fhr_run = h2_2[ir][ilayEff]->GetBinContent(ibinx,ibiny);
	    if(fhr_run>1e-15) {
	      stavefound = ibiny-1;
	      chipfound = ibinx-1;
	      break;
	    }
	  }
	  if(fhr_run>1e-15) break;
	}
      }
      if(fhr_run<1e-15){
	ntrig.push_back(-1);
	cout<<"Run "<<runnumbers[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir]<<" for Layer " << ilayEff << " has "<<ntrig[ntrig.size()-1]<<" triggers (FHR < 10-15 for all staves and chips/HICs) (ignored in the calculation of the average fhr)"<<endl;
	continue;
	//cout<<"INVALID FHR... setting it to -1"<<endl;
      }

      MaxRange =  hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(0)->GetXmax();
      MaxRangeY =  hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(1)->GetXmax();

      chipfoundeff = chipfound;
      if (ilayEff > 2) {
	NChipsPerHIC=7;  
	ChipRowsPerHIC=2;
	if (ilayEff < 5) HICsPerHS=4;
	else  HICsPerHS=7;
	if (chipfound < HICsPerHS) HalfStaveFound=0;
	else {
	  HalfStaveFound =1; 
	  chipfoundeff = chipfound-HICsPerHS;
	}
      }

      hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(0) ->SetRange(1+1024*chipfoundeff*NChipsPerHIC, NChipsPerHIC*1024*(chipfoundeff+1));
      hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(1) ->SetRange(1+512*HalfStaveFound*ChipRowsPerHIC, ChipRowsPerHIC*512*(HalfStaveFound+1));

      TH2F *hprojsparse = (TH2F*)hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir]->Projection(1,0);
      hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(0) ->SetRange(1, MaxRange);//reset the range
      hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(1) ->SetRange(1, MaxRangeY);//reset the range

      double hits_chip = hprojsparse->Integral(1,1024*NChipsPerHIC,1,512*ChipRowsPerHIC);
      delete hprojsparse;

      /*
      cout << "\nlayer " << ilay << " nchips/HIC " << NChipsPerHIC << " chipRows/HIC " << ChipRowsPerHIC << endl;
      cout << "nEntrieshmaps " << nEntrieshmaps << endl;
      cout << "stavefound " << stavefound << " chipfound " << chipfound << " half stave found " << HalfStaveFound << endl;
      cout << "run " << runnumbers[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] << endl;
      cout << "stave " <<       stavenums[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] << endl;
      cout << "hits_chip " << hits_chip << endl;
      cout << "X range " << hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(0) ->GetXmax() << endl;
      cout << "Y range " << hmaps[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir] -> GetAxis(1) ->GetXmax() << endl;
      */
      if(hits_chip/(ChipRowsPerHIC*512.*NChipsPerHIC*1024.*fhr_run) < 1e-15){//to avoid bugs due to bad runs
	ntrig.push_back(-1);
	cout<<"Run "<<runnumbers[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir]<<" Layer " << ilayEff << " has "<<ntrig[ntrig.size()-1]<<" triggers (ignored in the calculation of the average fhr)"<<endl;
      }
      else{
	ntrig.push_back(hits_chip/(ChipRowsPerHIC*512.*NChipsPerHIC*1024.*fhr_run));
	cout<<"Run "<<runnumbers[nEntrieshmaps+stavefound*(nRunsB[ilayEff]+1)+ir]<<" Layer " << ilayEff << " has "<<ntrig[ntrig.size()-1]<<" trigger" << endl;
      }
    }
  }

  //Start masking hottest pixels for each stave in each run, Fill also the histo with the hot pixel maps for each stave
  cout<<endl;
  cout<<"... Analysing FHR with (hot) pixel masking"<<endl;
  vector<array<float,nMasked+1>> fhrall;
  vector<array<float,nMasked+1>> fhrall1;

  TH2F *hHotMap[nLayers][48];
  TH2F *hHotMapCloned[nLayers][48];
  int index =0;
  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + stoi(laynums[0]) ;
    else if (IBorOB==2) ilayEff = ilay + stoi(laynums[0]) ;
    else ilayEff = ilay;
    if (nRunsB[ilayEff] == -1) continue;
    //    cout <<     hmaps[index] ->GetAxis(0)->GetXmax() << " " << hmaps[index] ->GetAxis(1)->GetXmax() << endl;
    for(int istave=0; istave<48; istave++){
      if (isHotPixelMapDrawn) hHotMap[ilay][istave] = new TH2F(Form("hHotMap_L%i_Stv%d",ilayEff, istave), "; ; ", (int)hmaps[index] ->GetAxis(0)->GetXmax(),0,hmaps[index] ->GetAxis(0)->GetXmax(), (int)hmaps[index] ->GetAxis(1)->GetXmax(),0,hmaps[index] ->GetAxis(1)->GetXmax());
    }
    index += (nRunsB[ilayEff]+1)*nStavesInLayAll[ilayEff];
  }

  cout << "\n\n" << endl;
  for (int i=0; i<nRunsTotAllLayers; i++){
    cout << "i " << i << "ntrig " << ntrig[i] << endl;
  }
  cout << "\n\n" << endl;

  int numStavePart = 2;
  int irun = nRunsTotAllLayers-1; 
  int nChips =0;
  bool IB = 1;
  int  ilayIndex =0;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    //    if (ihist < (int)hmaps.size()-20) continue;
    cout << "\nihist " << ihist << " laynum " << laynums[ihist] << " stave num " << stavenums[ihist] << " runnumber " << runnumbers[ihist] << " irun " << irun << " ntrigger " << ntrig[irun] << endl;
    if (stoi(laynums[ihist]) < 3)  { //IB layers                                                   
      nChips=9; 
      MaxRange = 9216; // 1024 x nChips=9
      IB = 1;
      numStavePart = 1;
    }
    else if (stoi(laynums[ihist]) ==3 || stoi(laynums[ihist]) ==4) { //L3 and L4
      nChips=2*28.;  //#chips in a half stave: 7 chips in one HIC row x 2 chip rows in one HIC x 4 HICs in a half stave
      MaxRange = 28672; //1024 x 7 chips in one HIC row x 4 HICs in a half stave
      IB = 0;
      numStavePart =2; //HS Lower and HS Upper
    }
    else if (stoi(laynums[ihist]) ==5 || stoi(laynums[ihist]) ==6)  { //L5 and L6
      nChips=2*49.;  //#chips in a half stave: 7 chips in one HIC row x 2 chip rows in one HIC x 7 HICs in a half stave
      MaxRange = 50176; //1024 x 7 chips in one HIC row x 7 HICs in a half stave
      IB = 0;
      numStavePart =2; //HS Lower and HS Upper
    }
    cout << "numStavePart " << numStavePart << endl;
    //cout<<"\n*************Layer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<<" Run: "<<runnumbers[ihist]<< endl;
    for (int StavePart=0; StavePart< numStavePart; StavePart++){ //loop over the two Half Staves for OB      
      int nchipsactive = GetNchipsActive(hmaps[ihist],nChips, MaxRange, StavePart, IB);
      if (numStavePart==1){
      if(nchipsactive<nChips) cout<<"\nLayer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<<" Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nchipsactive<<endl;
      }
      else {
	if(nchipsactive<nChips) cout<<"\nLayer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<< " Half Stave " << StavePart << " Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nchipsactive<<endl;
      }
      ilayIndex = stoi(laynums[ihist]);
      if (IBorOB == 1)       ilayIndex = stoi(laynums[ihist]) -3;
      cout << "\nilayIndex " << ilayIndex << endl;
      if (StavePart==0)   {
	fhrall.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : ilayIndex][stoi(stavenums[ihist])], StavePart, IB, isHotPixelMapDrawn));
	if (numStavePart==1) fhrall1.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : ilayIndex][stoi(stavenums[ihist])], StavePart, IB,  isHotPixelMapDrawn));
      }
      else  {
	fhrall1.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : ilayIndex][stoi(stavenums[ihist])], StavePart, IB,  isHotPixelMapDrawn));
      }
    }
    irun-- ;

    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
        irun=nRunsTotAllLayers-1;
      }
      if(laynums[ihist-1]!=laynums[ihist]){
        nRunsTotAllLayers -= (nRunsB[stoi(laynums[ihist])]+1);
      }
    }
    cout << "\nfhr all (HS Upper) " << fhrall[fhrall.size()-1][0] << endl;
    cout << "fhr all (HS Lower) " << fhrall1[fhrall1.size()-1][0] << endl;
  }

  //special binning
  double binstart = 0.4;
  double binsmasked[nMasked+2];
  binsmasked[0] = 0.1;
  binsmasked[1] = 0.5;
  for(int i=2; i<=nMasked+1; i++){
    binsmasked[i] = binsmasked[i-1]+1.;
  }

  TH2F *hFhrStv[nLayers][100][2];

  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + stoi(laynums[0]) ;
    else if (IBorOB==2) ilayEff = ilay + stoi(laynums[0]) ;
    else ilayEff = ilay;
    if (nRunsB[ilayEff] == -1) continue;
    //    cout << "ilay " << ilay << " ilayEff " << ilayEff << endl;
    for(int is=0; is<nStavesInLayAll[ilayEff]; is++){
      for (int i=0; i<2; i++){
	hFhrStv[ilay][is][i] = new TH2F(Form("h2FhrStv_%i_%d_HS%i", ilayEff,is, i), Form("Layer-%i - Stave-%d; # Hot Pixel Masked;Run", ilayEff,is),nMasked+1, binsmasked, nRunsB[ilayEff]+1, 0.5, nRunsB[ilayEff]+1.5);
      }
    }
  }

  //Fill histogram
  int ilayer = nLayers-1;
  int istave = 0;
  irun=0;
  //  cout << "\n\nfhr size " << (int)fhrall.size() << " = should correspond to nEntrieshmaps " << endl;
  //  cout << "fhr size (nMasked) " << (int)fhrall[0].size() <<" = should correspond to nMasked+1 " << endl;
  //  cout << "hmpas size\n " << hmaps.size() << endl;
  for(int i=0; i<(int)fhrall.size(); i++){
    istave = stoi(stavenums[(int)hmaps.size()-1-i]);
    //    cout << "i " << i << " ilayer " << ilayer << " istave " << istave << " irun " << irun << endl;
    //    cout << " ilayer " << laynums[(int)hmaps.size()-1-i] << " istave " << stavenums[(int)hmaps.size()-1-i] << " irun " << runnumbers[(int)hmaps.size()-1-i] << endl;
    for(int ifhr=0; ifhr<(int)fhrall[i].size(); ifhr++){ //loop over number of pixels masked
      if (ifhr ==0) cout << "fhr " << fhrall[i][ifhr] << endl; //they match fhr from loop above
      for (int StavePart=0; StavePart< 2; StavePart++){ //loop over the two Half Staves for OB  
	if (StavePart==0) hFhrStv[ilayer][istave][StavePart]->SetBinContent(ifhr+1, irun+1, fhrall[i][ifhr]);
	else hFhrStv[ilayer][istave][StavePart]->SetBinContent(ifhr+1, irun+1, fhrall1[i][ifhr]);
      }
    }

    irun++;
    if(i<(int)fhrall.size()-1){// in case ref run is the first into the list of runs
      if(stavenums[fhrall.size()-i-2]!=stavenums[fhrall.size()-i-1]){
        istave--;
        irun=0;
      }
      if(laynums[fhrall.size()-i-2]!=laynums[fhrall.size()-i-1]){
        ilayer--;
        istave = nStavesInLayAll[ilayer]-1;
      }
    }
  }

  
  cout << "\nFHR vs #masked pixels " << endl;
  //Make FHR (averaged on all runs) vs #masked pix for all staves in a layer
  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + stoi(laynums[0]) ;
    else if (IBorOB==2) ilayEff = ilay + stoi(laynums[0]) ;
    else ilayEff = ilay;
    if (nRunsB[ilayEff] == -1) continue;
    cout << "ilay " << ilay << " ilayEff " << ilayEff << endl;

    //legend
    TLegend *leg = new TLegend(0.904, 0.127,0.997,0.898);
    if (ilayEff>=3) leg->SetNColumns(2);

    if (IBorOB==0) numStavePart=1;
    else if (IBorOB==1) numStavePart=2;
    else {
      if (ilayEff < 3) numStavePart=1;
      else numStavePart=2;
    }

    TString SStavePart[2] = {"HS Lower", "HS Upper"};
    if (numStavePart==1) {
      SStavePart[0] = "";
      SStavePart[1] = "";
    }

    for (int StavePart=0; StavePart< numStavePart; StavePart++){ //loop over the two Half Staves for OB      
      TCanvas cnv(Form("cnv_%d_HS%i",ilayEff, StavePart), Form("cnv_%d_HS%i",ilayEff, StavePart));
      cnv.cd();
      cnv.SetLogy();
      cnv.SetLogx();
      cnv.SetTickx();
      cnv.SetTicky();
      cnv.SetMargin(0.0988,0.1,0.1,0.0993);

      TH1F *hframe = cnv.DrawFrame(0.1,7e-15,3*(nMasked),1e-3, Form("Layer %i %s - Average FHR %s; # Hot Pixels masked ; FHR (/event/pixel)",ilayEff, SStavePart[StavePart].Data(), filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
      hframe->SetBins(nMasked+1,binsmasked);

      for(int is=0; is<nStavesInLayAll[ilayEff];is++){
	TH1F *proj = (TH1F*)hFhrStv[ilay][is][StavePart]->ProjectionX(Form("proj_%d%d_HS%i",ilayEff,is, StavePart));
	int runswohits = GetNrunsWOhits(hFhrStv[ilay][is][StavePart]);
	proj->Scale(1./(nRunsB[ilayEff]+1-runswohits)); //Divide by the number of runs minus the ones without hits

	if (ilayEff<3){
	  SetStyle(proj, col[is<nStavesInLayAll[ilayEff]/2 ? is : is-nStavesInLayAll[ilayEff]/2],is<nStavesInLayAll[ilayEff]/2 ? 24:26);
	}
	else if (ilayEff >= 3 && ilayEff<5){
	  if((is)<nStavesInLayAll[ilayEff]/3){
	    SetStyle(proj, col[is], 24);
	  }
	  else if ((is)<nStavesInLayAll[ilayEff]*2/3){
	    SetStyle(proj, col[is-nStavesInLayAll[ilayEff]/3], 26);
	  }
	  else{
	    SetStyle(proj, col[is-nStavesInLayAll[ilayEff]*2/3], 25);
	  }
	}
	else {
	  if((is)<int(nStavesInLayAll[ilayEff]/4))
	    SetStyle(proj, col[is], 24);
	  else if ((is)<2*int(nStavesInLayAll[ilayEff]/4))
	    SetStyle(proj, col[is-int(nStavesInLayAll[ilayEff]/4)], 26);
	  else if ((is)<3*nStavesInLayAll[ilayEff]/4)
	    SetStyle(proj, col[is-2*int(nStavesInLayAll[ilayEff]/4)], 25);
	  else
	    SetStyle(proj, col[is-int(nStavesInLayAll[ilayEff]*3/4)], 30);
	}
	proj->Draw("PL same");

	if (StavePart==0) leg->AddEntry(proj, Form("Stv%d",is),"p");
      }
      leg->Draw("same");
      TString NameCnv = "";

      if (StavePart==0){
	if (numStavePart==1) NameCnv = Form("../Plots/Layer%i_FHRpixmask_%s", ilayEff,filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
	else NameCnv = Form("../Plots/Layer%i_HSLower_FHRpixmask_%s", ilayEff, filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
      }
      else NameCnv = Form("../Plots/Layer%i_HSUpper_FHRpixmask_%s", ilayEff, filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
      cnv.SaveAs(NameCnv + ".pdf");
      cnv.SaveAs(NameCnv + ".root");
    }
  }

  //  cout << "Running time up to here: " << endl;
  //  t1.Stop();
  //  t1.Print();

  if (isHotPixelMapDrawn){
    cout << "\nDrawing hot pixel map for each layer " << endl;
    for(int ilay=0; ilay<nLayers; ilay++){
      if (nLayers==1) ilayEff = stoi(laynums[0]);
      else if (IBorOB==1) ilayEff = ilay + stoi(laynums[0]) ;
      else if (IBorOB==2) ilayEff = ilay + stoi(laynums[0]) ;
      else ilayEff = ilay;
      if (nRunsB[ilayEff] == -1) continue;
      cout << "ilay " << ilay << " ilayEff " << ilayEff << endl;

      TCanvas cnvT(Form("cnvT_%d",ilayEff), Form("cnv_%d",ilayEff),800,1200);
      TCanvas cnvB(Form("cnvB_%d",ilayEff), Form("cnvBOT_%d",ilayEff),800,1200);
      cnvT.SetTopMargin(0.4);
      if (ilayEff<3) cnvT.Divide(1,nStavesInLayAll[ilayEff],0,0);
      else cnvT.Divide(1,nStavesInLayAll[ilayEff]/2,0,0);
      cnvB.SetTopMargin(0.4);
      cnvB.Divide(1,nStavesInLayAll[ilayEff]/2,0,0);
      if (ilayEff >= 3) cnvT.SetTitle(Form("cnvTOP_%d",ilayEff));

      int istaveeff=0;
      for(int istave=0; istave<nStavesInLayAll[ilayEff]; istave++){
	istaveeff = istave;
	hHotMap[ilay][istave]->SetMarkerStyle(20);
	hHotMap[ilay][istave]->SetMarkerSize(0.6);
	hHotMap[ilay][istave]->SetMarkerColor(kRed);
	hHotMap[ilay][istave]->SetLineColor(kRed);

	if (ilayEff<3 || istave < nStavesInLayAll[ilayEff]/2) {
	  cnvT.cd(istave+1);
	  cnvT.GetPad(istave+1)->SetTickx();
	  cnvT.GetPad(istave+1)->SetTicky();
	  cnvT.GetPad(istave+1)->SetRightMargin(0.01);
	  if(!istave) cnvT.GetPad(istave+1)->SetTopMargin(0.3);
	}
	else {
	  istaveeff = istave - nStavesInLayAll[ilayEff]/2;
	  cnvB.cd(istaveeff+1);
	  cnvB.GetPad(istaveeff+1)->SetTickx();
	  cnvB.GetPad(istaveeff+1)->SetTicky();
	  cnvB.GetPad(istaveeff+1)->SetRightMargin(0.01);
	  if(!istaveeff) cnvB.GetPad(istaveeff+1)->SetTopMargin(0.3);
	}
	hHotMap[ilay][istave]->Rebin2D(4);
	hHotMapCloned[ilay][istave]= (TH2F*) hHotMap[ilay][istave]->Clone(Form("%s_", hHotMap[ilay][istave]->GetName()));
	hHotMap[ilay][istave]->Draw("P X+");
	hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.005);
	hHotMap[ilay][istave]->GetYaxis()->SetTickLength(0.005);
	hHotMap[ilay][istave]->GetYaxis()->SetLabelSize(0.13);
	hHotMap[ilay][istave]->GetXaxis()->SetLabelSize(0.13);
	if ( (numStavePart==1 && istave ==0) || (numStavePart ==2 && (istave == 0 || istave == nStavesInLayAll[ilayEff]/2)) ){
	  hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(0.003);
	  hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
	  hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
	}
	else{
	  hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(999);
	  hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
	  hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
	}
	if (numStavePart == 2)   hHotMap[ilay][istave]->GetYaxis()->SetNdivisions(8);

	TLatex lat;
	lat.SetTextAngle(90);
	lat.SetNDC();
	if (numStavePart==1) {
	  lat.SetTextSize(0.15);
	  lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
	}
	else {
	  if (istave < nStavesInLayAll[ilayEff]/2){
	    lat.SetTextSize(0.25);
	    lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
	  }
	  else {
	    lat.SetTextSize(0.25);
	    lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
	  }
	}
      }

      cnvT.cd();
      TLatex latT;
      latT.SetNDC();
      latT.SetTextSize(0.03);
      if (numStavePart==1)    latT.DrawLatex(0.01,0.98,Form("L%i",ilayEff));
      else     latT.DrawLatex(0.01,0.98,Form("L%i - TOP", ilayEff));

      if (numStavePart==2){
	cnvB.cd();
	TLatex latB;
	latB.SetNDC();
	latB.SetTextSize(0.03);
	latB.DrawLatex(0.01,0.98,Form("L%i - BOTTOM", ilayEff));
      }

      TFile * fileOutputT;
      TFile * fileOutputB;
      TString NameCnvT="";
      TString NameCnvB="";
      if (ilayEff<3) {
	NameCnvT = Form("../Plots/Layer%i_Hotpixmap_%s", ilayEff,filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
      }
      else {
	NameCnvT = Form("../Plots/Layer%i-TOP_Hotpixmap_%s", ilayEff,filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
	NameCnvB = Form("../Plots/Layer%i-BOT_Hotpixmap_%s", ilayEff,filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
	cnvB.SaveAs(NameCnvB + ".pdf");
	fileOutputB = new TFile(NameCnvB+ ".root", "RECREATE");
      }
      cnvT.SaveAs(NameCnvT + ".pdf");
      fileOutputT = new TFile(NameCnvT + ".root", "RECREATE");

      for(int istave=0; istave<nStavesInLayAll[ilayEff]; istave++){
	if (ilayEff < 3 || istave < nStavesInLayAll[ilayEff]/2) {
	  fileOutputT->WriteTObject(hHotMapCloned[ilay][istave]);
	}
	else {
	  fileOutputB->WriteTObject(hHotMapCloned[ilay][istave]);
	}
      }
      fileOutputT->Close();
      if (ilayEff>=3)    fileOutputB->Close();
      cout << "The following root files have been created:\n" << NameCnvT << ".root" << endl;
      if (ilayEff>=3) cout << NameCnvB << ".root"<< endl;
    }
  }
 
  //t.Stop();
  //cout << "\nRunning time: " << endl;
  //t.Print();
 
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap, bool HS , bool IB, bool isHotPixelMapDrawn){

  array<float,nMasked+1> fhrstave;

  THnSparse *hmapclone = (THnSparse*)hmap->Clone(Form("%s_clone",hmap->GetName()));
  //cout<<"\n" <<hmap->GetName()<<" --> "<<hmapclone->GetNbins()<<" noisy pixels"<<endl;
  int iyMin =1;
  int iyMax =hmapclone->GetAxis(1)->GetNbins();
  if (IB==0) {
    if (HS==0) iyMax = hmap->GetAxis(1)->GetNbins()/2;
    else iyMin = hmap->GetAxis(1)->GetNbins()/2 +1;
  }
  hmapclone->GetAxis(1)->SetRange(iyMin, iyMax);
  for(int iter=0; iter<nMasked+1; iter++){
    TH1F *hproj = (TH1F*)hmapclone->Projection(1);
    long int totalhits = hproj->Integral();
    float fhr = nchips==0 ? 0. : (float)totalhits / (512.*1024.*nchips*ntrig);

    if(ntrig<0) fhr=0.;
    fhrstave[iter] = fhr;

    int coord[2];
    double max = -1.;
    int x=0,y=0;
    long int binwithmax = 0;
    for(int ibin=0; ibin<hmapclone->GetNbins(); ibin++){
      double bincontent = hmapclone->GetBinContent(ibin, coord);
      if (IB==0){
	if (coord[1] < iyMin || coord[1] > iyMax) continue;
      }
      if(bincontent>max){
        max=bincontent;
        binwithmax = ibin;
        x=coord[0];
        y=coord[1];
      }
    }

    if(nchips) hmapclone->SetBinContent(binwithmax,0.);
    if (isHotPixelMapDrawn){
      if(totalhits!=0 && iter < 100) {
	hhotmap->SetBinContent(x, y, 1); // to avoid a marker in 0,0 for empty histos
      }
    }
    delete hproj;
  }

  delete hmapclone;

  return fhrstave;

}

//
//Function to return the number of active (i.e. enabled) chips in a stave (or in a half stave when dealing with OB)
//
int GetNchipsActive(THnSparse *hmap, int maxchip, int MaxRange, bool HS, bool IB){
  int ix1=1;
  int ix2=1024;
  int iy1=1;
  int iy2=512;
  int iy2Max = 0;
  int iy2Min = 1;
  int activechips = maxchip;
  iy2Max =hmap->GetAxis(1)->GetNbins();
  if (IB==0) {
    if (HS==0) iy2Max = hmap->GetAxis(1)->GetNbins()/2;
    else {
      iy1 = 1025;
      iy2 = 1024 + 512;
      iy2Min = hmap->GetAxis(1)->GetNbins()/2;
    }
  }
  while (iy2 <= iy2Max && iy1 >= iy2Min){
    ix1 = 1;
    ix2=1024;
    while(ix2<=hmap->GetAxis(0)->GetNbins()){
      hmap->GetAxis(0)->SetRange(ix1,ix2);
      hmap->GetAxis(1)->SetRange(iy1,iy2);
      TH2F *hproj = (TH2F*)hmap->Projection(1,0);
      if((hproj->Integral(1,1024,1,512)<1e-15)){
	activechips--;
      }
      ix1=ix2+1;
      ix2+=1024;
      delete hproj;
    }
    iy1=iy2+1;
    iy2+=512;
  }
  hmap->GetAxis(0)->SetRange(1,MaxRange);//reset range
  return activechips;
}

//
//Function to return the number of runs in which a stave has no hits
//
int GetNrunsWOhits(TH2 *hFhrStv){
  int runswohits = 0;
  for(int iy=1; iy<=hFhrStv->GetNbinsY(); iy++){
    if(hFhrStv->GetBinContent(1,iy)<1e-15)
      runswohits++;
  }
  return runswohits;
}

//
//Set Style
//
void SetStyle(TH1 *h, int col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}
