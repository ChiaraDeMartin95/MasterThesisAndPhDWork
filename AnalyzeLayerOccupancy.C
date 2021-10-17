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

void SetStyle(TGraph *h, int col, Style_t mkr);
void DoAnalysis(string filepath, const int nChips, string skipruns, int IBorOB);

//
// MAIN
//
void AnalyzeLayerOccupancy(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
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


  //Call
  DoAnalysis(fpath, nchips, skipruns, IBorOB);

}

//
//Set Style
//
void SetStyle(TGraph *h, int col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, int nChips, string skipruns, int IBorOB){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, laynums;
  string HS[2] = {"Upper", "Lower"};
  int nTimes=0;
  int nRunsB[7]={-1};
  for (int ilay=0; ilay < 7; ilay++){
    nRunsB[ilay] =-1;
  }
  int nRunsTot=0;
  int nRunsTotFixed=0;
  int nLayersInput =1;
  int col[] = {810, 807, 797, 827, 417, 841, 868, 867, 860, 602, 921, 874};

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
    if(objname.find("Stv")!=string::npos) break;
    h2 = (TH2*)obj;
    if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    if(skipruns.find(runnum)!=string::npos) continue;// eventually skip runs if specified
    hmaps.push_back(h2);
    //cout<<runnum<<"  "<<timestamp<<endl;
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);

    laynums.push_back(laynum);
    //    cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
    nTimes++;
    if (nRunsB[stoi(laynum)]==-1) nRunsB[stoi(laynum)]=0;
    if(nTimes>1){
      if (laynum==laynums[laynums.size()-2]) {
        nRunsB[stoi(laynum)]++;
      }
      else nLayersInput++;
    }
  }

  if (nTimes==0) {cout <<"\nInput files contains no plots\n" << endl; return;}

  /*
  for (int i=0; i<7; i++){
    cout <<"nRunsB " << nRunsB[i] << endl;
  }
  */

  int nLayers;
  if (nLayersInput==1) nLayers =1;
  else {
    if (IBorOB==0) nLayers = 3;
    else if (IBorOB==1) nLayers = 4;
    else nLayers = 7;
  }

  TGraph *trend[nLayers][100][2];
  for (int ilay =0; ilay<7; ilay++){
    trend[ilay][0][0] = new TGraph(); //to get npoints=0 for layers with no data
  }
  int ilayer=0;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){
    if (IBorOB==1)  ilayer = stoi(laynums[ihist])-3;
    else ilayer = stoi(laynums[ihist]);
    if (nLayersInput==1) ilayer = 0;
    //    cout << "ilayer " <<ilayer << endl;
    for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){
      for (int b=0; b<=1; b++){
	trend[ilayer][ibiny-1][b] = new TGraph();
      }
    }
  }

  TH1F *hproj = new TH1F();
  string histname = hmaps[0]->GetName();
  int irun=0;

  int ChipMin =1;
  int ChipMax =1;
  int numStavePart =0;

  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    if (IBorOB==1)  ilayer = stoi(laynums[ihist])-3;
    else ilayer = stoi(laynums[ihist]);
    if (nLayersInput==1) ilayer = 0;

    if (stoi(laynums[ihist]) < 3)   nChips=9; //IB layers
    else if (stoi(laynums[ihist]) ==3 || stoi(laynums[ihist]) ==4)    nChips=8./2;  //L3 and L4, Half stave
    else if (stoi(laynums[ihist]) ==5 || stoi(laynums[ihist]) ==6)    nChips=14./2; //L5 and L6, Half Stave

    if (stoi(laynums[ihist]) >=3) numStavePart = 1;
    else numStavePart=0;

    for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){//loop on y bins (stave)s
      //cout << "\n stave number: " << ibiny-1 << endl;

      TH1D *hproj = hmaps[ihist]->ProjectionX("proj",ibiny,ibiny); //single stave
      int deadchips = 0;

      for (int StavePart=0; StavePart<= numStavePart; StavePart++){ //loop over the two Half Staves for OB
	trend[ilayer][ibiny-1][StavePart]->SetName(Form("gr_L%s_stave%d_HS%i",laynums[ihist].c_str(),ibiny-1, StavePart));

	ChipMin =1;
	ChipMax =hmaps[ihist]->GetNbinsX() ;
	if (stoi(laynums[ihist]) >= 3){
	  if (StavePart ==0) { //HS Lower
	    ChipMin =1;
	    ChipMax = hmaps[ihist]->GetNbinsX()/2;
	  }
	  else { //HS Upper
	    ChipMin =hmaps[ihist]->GetNbinsX()/2 +1 ;
	    ChipMax =hmaps[ihist]->GetNbinsX() ;
	  }
	}

	deadchips = 0;
	for(int ibinx=ChipMin; ibinx<=ChipMax; ibinx++){//evaluate the number of disabled chips
	  if(hmaps[ihist]->GetBinContent(ibinx,ibiny)<1e-15)
	    deadchips++;
	}
	if(deadchips>0){
	  if (stoi(laynums[ihist]) < 3) {
	    cout<<"Layer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" Number of chips: " << nChips<<" --> Chips active:"<<nChips-deadchips<<endl;
	  } else {
	    cout<<"Layer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" Number of HICs: " << nChips << " --> HICs active:"<<nChips-deadchips<<endl;
	  }
	}

	double StaveOccupancy=0;
	if(deadchips!=nChips){
	  for (int b=ChipMin; b<= ChipMax; b++){
	    StaveOccupancy+= hproj->GetBinContent(b);
	  }
	  //	  cout << "full int " <<hproj->Integral() << " StaveOccupancy " << StaveOccupancy<<  " average stave occupancy per chip " <<  StaveOccupancy/(nChips-deadchips) << endl;
	  trend[ilayer][ibiny-1][StavePart]->SetPoint(irun, irun, StaveOccupancy/(nChips-deadchips));
	  //	  cout << "layer " << laynums[ihist] << " run " <<runnumbers[ihist] << " stave " << ibiny-1 << " " <<  StaveOccupancy/(nChips-deadchips)<< endl;
	}
	else
	  trend[ilayer][ibiny-1][StavePart]->SetPoint(irun, irun, 0.);

	double x=0;
	double y=0;

	if (stoi(laynums[ihist]) < 3){
	  if((ibiny-1)<hmaps[ihist]->GetNbinsY()/2)
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1], 24);
	  else
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1-hmaps[ihist]->GetNbinsY()/2], 26);
	}
	else if (stoi(laynums[ihist]) ==3 || stoi(laynums[ihist])==4){
	  if((ibiny-1)<hmaps[ihist]->GetNbinsY()/3){
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1], 24);
	  }
	  else if ((ibiny-1)<hmaps[ihist]->GetNbinsY()*2/3){
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1-hmaps[ihist]->GetNbinsY()/3], 26);
	  }
	  else{
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1-hmaps[ihist]->GetNbinsY()*2/3], 25);
	  }
	}
	else if (stoi(laynums[ihist]) ==5 || stoi(laynums[ihist])==6){
	  if((ibiny-1)<int(hmaps[ihist]->GetNbinsY()/4))
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1], 24);
	  else if ((ibiny-1)<2*int(hmaps[ihist]->GetNbinsY()/4))
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1-int(hmaps[ihist]->GetNbinsY()/4)], 26);
	  else if ((ibiny-1)<3*hmaps[ihist]->GetNbinsY()/4)
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1-2*int(hmaps[ihist]->GetNbinsY()/4)], 25);
	  else
	    SetStyle(trend[ilayer][ibiny-1][StavePart], col[ibiny-1-int(hmaps[ihist]->GetNbinsY()*3/4)], 30);
	}
	//	cout << "trend values for run "<< irun<< " "  << trend[ilayer][ibiny-1][StavePart]->GetPointY(irun) << endl;

      }//Stavepart end
    }
    irun++;
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
	irun=0;
      }
  }

  //  cout << "nLayers " <<  nLayers << endl;
  int ilayEff=0;
  int npoints=0;
  TH1F *hfake[7];
  nRunsTot =0;
  for(int ilay=0; ilay<nLayers; ilay++){    
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    if (nRunsB[ilayEff]!=-1) nRunsTot += (nRunsB[ilayEff]+1);   
    npoints    = trend[ilay][0][0]->GetN();
    //    cout << "ilay " << ilay << " ilayEff " << ilayEff << " nRunsTot " << nRunsTot<< endl;
    if (npoints==0) continue;
    hfake[ilay]= new TH1F(Form("hfake_L%i", ilay), "; Run; Fake-hit Rate (/event/pixel)", npoints, -0.5, (double)npoints-0.5);
    for(int ir=0; ir<=nRunsB[ilayEff]; ir++) {
      hfake[ilay]->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[nRunsTot-ir-1])));
    }
  }
  nRunsTotFixed = nRunsTot;
  //  cout << "Total number of runs analysed " << nRunsTotFixed<< endl;

  TString PathOut = Form("../Plots/FHR_%s.root",filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str());
  TFile * fileout = new TFile(PathOut, "RECREATE");

  //Draw
  nRunsTot =0;
  for(int ilay=0; ilay<nLayers; ilay++){
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    if (ilay>0 && nRunsB[ilayEff-1]!=-1) nRunsTot+= (nRunsB[ilayEff-1]+1);
    npoints    = trend[ilay][0][0]->GetN();
    if (npoints==0) continue;
    TCanvas *canvas = new TCanvas("canvas", "canvas");
    TCanvas *Secondcanvas = new TCanvas("Secondcanvas", "Secondcanvas");
    if (ilayEff>=3) {
      Secondcanvas->SetName("canvas_HSUpper");
      canvas->SetName("canvas_HSLower");
    }
    canvas->SetLogy();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.0988,0.15,0.194,0.0993);
    Secondcanvas->SetLogy();
    Secondcanvas->SetTickx();
    Secondcanvas->SetTicky();
    Secondcanvas->SetMargin(0.0988,0.15,0.194,0.0993);

    Bool_t IsTwoCanvas=0;
    if (ilayEff>=3) IsTwoCanvas = 1;

    TLegend *leg = new TLegend(0.857, 0.197,0.997,0.898);
    if (ilayEff>=3) leg->SetNColumns(2);  

    for(int istave=0; istave<hmaps[nRunsTot]->GetNbinsY(); istave++)
      leg->AddEntry(trend[ilay][istave][0], Form("Stv%d",istave), "p");
    hfake[ilay]->GetYaxis()->SetRangeUser(1e-14, 1e-2);
    hfake[ilay]->GetXaxis()->SetTitleOffset(2);
    hfake[ilay]->SetTitle(Form("Layer-%s, %s",laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->cd();
    if (IsTwoCanvas) hfake[ilay]->SetTitle(Form("Layer-%s, HS Lower, %s",laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    hfake[ilay]->Draw();
    for(int istave=0; istave<hmaps[nRunsTot]->GetNbinsY(); istave++){
      trend[ilay][istave][0]->Draw("P same");
      fileout->WriteTObject(trend[ilay][istave][0]);
    }
    leg->Draw("same");
    if (!IsTwoCanvas){
      canvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s.pdf", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      canvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s.root", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    }
    else {
      canvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s_HSLower.pdf", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      canvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s_HSLower.root", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    }

    if (IsTwoCanvas){
      Secondcanvas->cd();
      hfake[ilay]->SetTitle(Form("Layer-%s, HS Upper, %s",laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      hfake[ilay]->Draw();
      for(int istave=0; istave<hmaps[nRunsTot]->GetNbinsY(); istave++){
	trend[ilay][istave][1]->Draw("P same");      
	fileout->WriteTObject(trend[ilay][istave][1]);
      }

      leg->Draw("same");

      Secondcanvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s_HSUpper.pdf", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      Secondcanvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s_HSUpper.root", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    }

    delete canvas;
    delete Secondcanvas;
    delete leg;
  }
  fileout->Close();
  
  //Make GIF with TH2 for each run and for each layer
  nRunsTot =0;
  irun=1;
  gStyle->SetPalette(1);
  
  for(int ilay=0; ilay<nLayers; ilay++) {//remove images if they exist already
    if (nLayers==1) ilayEff = stoi(laynums[0]);
    else if (IBorOB==1) ilayEff = ilay + 3 ;
    else ilayEff = ilay;
    npoints    = trend[ilay][0][0]->GetN();
    if (ilay>0 && nRunsB[ilayEff-1]!=-1) nRunsTot+= (nRunsB[ilayEff-1]+1);
    if (npoints==0) continue;
    gSystem->Unlink(Form("../Plots/Layer%s_fakehitratemap_%s.gif", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }


  nRunsTot = nRunsTotFixed-1; 
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    if (ihist!=((int)hmaps.size()-1) && (stoi(laynums[ihist])!=stoi(laynums[ihist+1]))) {
      if (nRunsB[stoi(laynums[ihist])]!=-1)      nRunsTot-= (nRunsB[stoi(laynums[ihist])]+1);
    }
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogz();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.15);
    hmaps[ihist]->Draw("colz");
    hmaps[ihist]->SetMinimum(1e-14);
    hmaps[ihist]->SetMaximum(1e-3);
    hmaps[ihist]->GetZaxis()->SetTitle("Fake-hit Rate (/event/pixel)");
    hmaps[ihist]->SetTitle(Form("Layer-%s, Run %06d (%d/%d)", laynums[nRunsTot].c_str(), stoi(runnumbers[ihist]), irun, nRunsB[stoi(laynums[ihist])]+1));
    canvas->Print(Form("../Plots/Layer%s_fakehitratemap_%s.gif+40", laynums[nRunsTot].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    irun++;
    if(!ihist && nLayers==1){
      break;
    }
    if(nLayers>1){
      if(ihist>0){
        if(laynums[ihist-1]!=laynums[ihist]){
          irun=1;
        }
      }
    }
    delete canvas;
  }

  cout << "\nROOT file " << PathOut << " has been created" << endl;
}
