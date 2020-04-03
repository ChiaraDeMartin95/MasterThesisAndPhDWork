#include <Riostream.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TCutG.h>
#include "TFitResult.h"
#include "TLegend.h"

Bool_t reject;
Double_t fparab(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 0.474 && x[0] < 0.520) {
      TF1::RejectPoint();
      return 0;
   }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
    //return par[0] + par[1]*x[0] + par[2]*0;
}


void histos_exercise_thesis( Float_t PtTrigMin=3.0, Int_t type=0, Int_t sysTrigger=0, Int_t sysV0=0, Int_t syst=0, Double_t nsigmamax=9, TString year0="2016", TString year="2018f1_extra_onlyTriggerWithHighestPt", TString Path1="", Bool_t isMC=1, Bool_t isEfficiency=1){ //molt va da 0 a 4, type=0 (=K0s), type=1 (=Lambda), type=2 (=AntiLambda), type=3 (=Lambda + Antilambda) , cut = 0 per tagli definitivi (quelli riportati nella presentazione per il tirocinio), cut=1 per tagli di default (quelli riportati nell'analysis note); opzione valida per i soli K0s, per altre particelle usare cut = 0


  // //lista degli effetti  sistematici studiati in questa macro
 //sys=1 nsigmamin=5 (def:4)
  //sys=2 sigmacentral =4 (def:3)

  if (syst!=0 && syst!=1 && syst!=2) {
    cout << "syst should be changed " << endl;
    return;
  }
  if((sysV0!=0||sysTrigger!=0) && syst!=0){
    cout << "syst conflicting " << endl;
    return;
  }

  //sys=1 nsigmamin=5 (def:4)
  //sys=2 sigmacentral =4 (def:3)
  Double_t sigmacentral=3;
  Double_t nsigmamin=4;
  
  if (syst==1) {
    nsigmamin=5;
  }
  if (syst==2) {
    sigmacentral=4;
  }

  cout << "*****************************************************************************"<< endl;
  cout << "Attenzione:  I limiti superiori di errsigma e errmean dipendono dalla tipologia di taglio implementato, in modo da garantire che i vari istogrammi siano rimepiti solo se il fit viene \"bene\" " << endl;

  if (isMC && !isEfficiency) cout <<"**********************ERRORE****************"<< endl;
  const Int_t mult=5; //numero intervalli molteplicita'
  const Int_t num_tipo=4; //numero tipologie di particelle
  const Int_t num_histo=7; //numero intervalli di Pt
  const Int_t num_tagli=6; //numero di diversi tagli applicati alle K0s
  const Int_t NsysTrigger=3; 
  const Int_t NsysV0=7; 

 
  Float_t massLambda = 1.115683;
  Float_t massK0= 0.497611;
  TString tipo[num_tipo]={"kK0s", "Lambda", "AntiLambda","LambdaAntiLambda" };
  TString tipo1[num_tipo]={"K0", "lambda", "antilambda","lambda+antilambda" };

  TString nome_TDir="MyTask";
 
  TString invmass[num_tipo]={"#pi^{+} #pi^{-}", "p #pi^{-}", "overline{p} #pi^{+}", "p #pi^{-} + overline{p} #pi^{+}"};
  //  TString int_pt[num_histo]={"0-0.5", "0.5-1","1-1.5","1.5-2", "2-3","3-4","4-5","5-6","6-7","7-8","8-10", "10-12", "12-16"};
  //  TString int_pt[num_histo]={"0-1", "1-2", "2-3", "3-4", "4-8"};
  TString int_pt[num_histo]={"0-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString multiplicity[mult+1]={"0-5","5-10", "10-30", "30-50", "50-100", "all"};
  TString invmass_title[num_tipo] = {"(#pi^{+}, #pi^{-}) invariant mass (GeV/c^{2})", "massa invariante (p, #pi^{-}) (GeV/c^{2})", "massa invariante (overline{p}, #pi^{+})", "massa invariante [(p, #pi^{-}) + (overline{p}, #pi^{+})] (GeV/c^{2})"};
  Float_t min_range_signal[num_tipo]={0.47,1.105,1.105,1.105}; //scelto ad occhio
  Float_t max_range_signal[num_tipo]={0.53,1.125,1.125,1.125}; //scelto ad occhio
  Int_t entries_range[num_histo]={0};
  Float_t bin_content0[num_histo]={0.};
  Float_t bin_content1[num_histo]={0.};
  Float_t bin_content2[num_histo]={0.};
  Float_t bin_contents3[num_histo]={0.};
  Float_t bin_contentSSB1[num_histo]={0.};
  Float_t bin_contentSSB2[num_histo]={0.};
  Float_t sigmas1[num_histo]={0};
  Float_t sigmab1[num_histo]={0};
  Float_t errs3[num_histo]={0};
  Float_t err0[num_histo]={0};
  Float_t err1[num_histo]={0};
  Float_t err2[num_histo]={0};
  Float_t err3[num_histo]={0};
  Float_t err4[num_histo]={0};
  Float_t errSSB1[num_histo]={0};
  Float_t errSSB2[num_histo]={0};
  Float_t s1[num_histo]={0};
  Float_t s2[num_histo]={0};
  Float_t st[num_histo]={0};
  Float_t tot[num_histo]={0};
  Float_t b[num_histo]={0};
  Float_t bside[num_histo]={0};
 Float_t errb[num_histo]={0};
 Float_t errbside[num_histo]={0};
  
 Float_t sigma[num_histo]={0};
 Float_t errsigma[num_histo]={0};
  Float_t mean[num_histo]={0};
  Float_t errmean[num_histo]={0};


  TFitResultPtr fFitResultPtr0[num_histo]; //serve per poter calcolare in seguito l'errore di un integrale 
  TFitResultPtr fFitResultPtr1[num_histo];

  // Float_t min_histo[num_tipo]={0.45,1.09,1.09,1.09};  //estremi del range degli istogrammi
  // Float_t max_histo[num_tipo]={0.55,1.14,1.14,1.14 };
  Float_t min_histo[num_tipo]={0.45,1.09,1.09,1.09};  //estremi del range degli istogrammi
  Float_t max_histo[num_tipo]={0.55,1.14,1.14,1.14 };
  
  Float_t lim_inf_mean[num_tipo]={0.495,1.1153,1.1153,1.1153 }; //come l'ho scelto?
  Float_t lim_sup_mean[num_tipo]={0.500,1.1168,1.1168,1.1168 };
  Float_t lim_inf_sigma[num_tipo]={0,0,0, 0};
  Float_t lim_sup_sigma[num_tipo]={0.008,0.002,0.002,0.002 }; //first one must be 0.008 if we want to display the 7-15 mult values between 7 and 8 pT
  Float_t lim_inf_errmean[num_tipo]={0,0,0, 0}; 
  //  Float_t lim_sup_errmean[num_tipo]={0.001,0.0006,0.00035,0.0006 }; //per Arm e tagli standard
  // Float_t lim_sup_errmean[num_tipo]={0.002,0.0006,0.00035,0.0006 };//per Lrejection
  Float_t lim_sup_errmean[num_tipo]={10,0.0006,0.00035,0.0006 };//loooooose
  Float_t lim_inf_errsigma[num_tipo]={0,0,0, 0}; 
  //  Float_t lim_sup_errsigma[num_tipo]={0.001,0.0004,0.0015,0.0004 }; //per Arm e tagli standard
  //  Float_t lim_sup_errsigma[num_tipo]={0.002,0.0004,0.0015,0.0004 }; //per Lrejection
   Float_t lim_sup_errsigma[num_tipo]={10,0.0004,0.0015,0.0004 }; //loose

  Double_t multip_intervals[mult+1]={0, 5, 10, 30, 50,100};
  //  Double_t binl[num_histo+1]={0, 1,1.5, 2,2.5,3,4,8}
  Double_t binl[num_histo+1]={0, 1,1.5, 2,2.5,3,4,8};

  Int_t rebin[num_tipo][mult+1][num_histo]={1};
  /*
    {{1,1,1,1,1,1,1,2,2,2,2,2,2},{1,1,1,1,1,1,1,1,2,2,2,2,2},{1,1,1,1,1,1,1,2,2,2,2,2,2},{1,1,1,1,1,1,1,2,2,2,2,2,2},{1,1,1,1,1,1,2,2,2,2,2,2,2} },
    {{2,2,2,2,2,4,4,4,4,4,4,4,4},{2,2,2,2,2,2,4,4,4,4,4,4,4},{2,2,2,2,2,2,4,4,4,6,6,6,6},{4,2,2,2,2,4,4,4,5,5,5,5,5},{4,4,2,2,2,4,4,4,4,4,4,4,4}},
    {{2,2,2,2,2,4,4,4,4,4,4,4,4},{2,2,2,2,2,2,4,4,4,4,4,4,4},{2,2,2,2,2,2,4,4,4,6,6,6,6},{4,2,2,2,2,4,4,4,5,5,5,5,5},{4,4,4,4,4,4,4,4,4,4,4,4,4}},
    {{2,2,2,2,2,4,4,6,6,6,6,6,6},{2,2,2,2,2,2,4,4,4,4,4,4,4},{2,2,2,2,2,2,4,4,4,6,6,6,6},{4,2,2,2,2,4,4,4,6,6,6,6,6},{4,4,4,4,4,4,4,4,4,4,4,4,4}}
  };
  /*
  Double_t range_bkg[num_tipo][mult][num_histo]={
    {{5,5,5,5,7,7,7,7,7,7,7,7,7},{5,5,5,5,5,7,7,7,7,7,7,7,7},{5,5,5,5,5,7,7,7,5,5,7,7,7},{5,5,5,5,5,7,7,7,7,7,7,7,7},{5,5,5,5,5,7,7,7,7,7,7,7,7}}, //{5,5,5,5,5,7,7,7,7,7,7,7,7} in nolabel era questo (terza graffa)
    {{7,7,7,7,7,7,7,7,7,7,7,7,7},{6,6,6,6,6,6,6,6,6,6,6,6,6},{7,5,5,5,5,5,7,7,7,7,7,7,7},{7,7,7,7,7,7,6,7,7,7,7,7,7},{7,7,7,7,7,7,7,7,7,7,7,7,7}},
    {{7,7,7,7,7,7,7,7,7,7,7,7,7},{6,6,6,6,6,6,7,6,6,6,6,6,6},{7,5,5,5,7,5,7,7,7,7,7,7,7},{7,7,7,7,7,7,7,7,7,7,7,7,7},{7,7,8,8,8,8,7,7,7,7,7,7,7}},
    {{7,5,5,5,5,7,7,7,7,7,7,7,7},{6,5,5,5,5,6,7,6,6,6,6,6,6},{7,5,5,5,7,5,7,7,7,7,7,7,7},{7,7,7,7,7,7,7,7,7,7,7,7,7},{7,7,8,8,8,8,7,7,7,7,7,7,7}}
  };
  */
  Double_t range_bkg[num_tipo][mult+1][num_histo];
  for(Int_t i=0; i< num_tipo; i++){
    for(Int_t j=0; j< mult+1; j++){
      for(Int_t l=0; l< num_histo; l++){
	range_bkg[i][j][l]=nsigmamax;
      }    
    }  
  }
  cout<< " ciao " << endl;
  TFile *myfile; 
  TFile *myfileAnalysis; 
  TCanvas *canv[mult+1];

  for(Int_t molt=0; molt<mult+1; molt++){
  TString nome_file_1 ="FinalOutput/DATA"+year0+"/histo/AngularCorrelation"+year;
  TString nome_file_output[NsysTrigger][NsysV0];
  TString nome_file_analysis;
  nome_file_analysis="FinalOutput/AnalysisResults"+year;
  nome_file_analysis+=Path1+".root";
  if (isMC && isEfficiency) {
    nome_file_analysis="FinalOutput/AnalysisResults"+year+"_MCEff";
    nome_file_analysis+=Path1+".root";
  }
  nome_file_output[sysTrigger][sysV0] ="FinalOutput/DATA"+year0+"/invmass_distribution_thesis/invmass_distribution";

 if(isMC && isEfficiency){
    nome_file_1+="_MCEff";
    nome_file_output[sysTrigger][sysV0]+="_MCEff";
  }
 nome_file_1+=Path1;
 nome_file_1 +=Form("_MassDistr_SysT%i_SysV0%i_PtMin%.1f.root",sysTrigger, sysV0, PtTrigMin);
 nome_file_output[sysTrigger][sysV0]+=Path1;
 nome_file_output[sysTrigger][sysV0] +=Form("_"+year+"_"+tipo[type]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", molt, sysTrigger, sysV0, syst,PtTrigMin);

 myfile = new TFile(nome_file_1, ""); 
 myfileAnalysis = new TFile(nome_file_analysis, ""); 
 TDirectoryFile *dirinput = (TDirectoryFile*)myfileAnalysis->Get(nome_TDir);
 TList *listinput = (TList*)dirinput->Get("MyOutputContainer");

  canv[molt] = new TCanvas(Form("canv_molt%i", molt),"Distribuzione di massa invariante di " +invmass[type] + " dopo applicazione tagli",1600,1000);
  canv[molt]->Divide(4,2);

  Double_t events=1;
  

 TFile *f = new TFile(nome_file_output[sysTrigger][sysV0],"RECREATE");
  cout << "ciao " << endl;
  //  TH2F *hMassvsPt = (TH2F*)listinput->Get(Form("SE_hMassvsPt_" + tipo[type]+"_%i",molt));
  TH2F *hMassvsPt_tagli = (TH2F*)myfile->Get(Form("SE_hMassvsPt_"+tipo[type]+"_%i",molt));

  TH1F *isto[num_histo];
  TH1F *isto_tagli[num_histo];
  cout<< " ciao " << endl;

  Double_t par[num_histo][9];  //se cambio funzione di bkg, devo cambiare il secondo numero
  cout<< " ciao " << endl;

  TF1 **functions1 = new TF1*[num_histo];  //perche due star
  TF1 **functions2 = new TF1*[num_histo];
  cout<< " ciao " << endl;

  TLegend *legend[num_histo];
  TCanvas* canvas[num_histo];
  TF1 **bkg1 = new TF1*[num_histo];
  TF1 **bkg2 = new TF1*[num_histo];
  TF1 **bkg3 = new TF1*[num_histo];
  TF1 **total= new TF1*[num_histo];
  TF1 **totalbis =new TF1*[num_histo]; //mi serve solo per calcolare errore del bkg
  TF1 **total_one= new TF1*[num_histo];
  TF1 **totalbis_one =new TF1*[num_histo]; //mi serve solo per calcolare errore del bkg
  TH1F *histo_sigma = new TH1F ("histo_sigma","Sigma vs Pt", num_histo, binl);
  TH1F *histo_mean = new TH1F ("histo_mean","Mean vs Pt", num_histo, binl);
  TH1F *histo_chis = new TH1F ("histo_chis","ChiSquare/dof vs Pt", num_histo, binl);
  TH1F *histo_SB = new TH1F ("histo_SB","S/B entro 3sigma vs Pt", num_histo, binl);
  TH1F *histo_SSB = new TH1F ("histo_SSB","S/(S+B) entro 3sigma vs Pt", num_histo, binl);
  TH1F *histo_S = new TH1F ("histo_S","S entro 3sigma vs Pt", num_histo, binl);
  TH1F *histo_Bcentral = new TH1F ("histo_Bcentral","B entro nsigma vs Pt", num_histo, binl);
  TH1F *histo_Bside = new TH1F ("histo_Bside","B in sideband vs Pt", num_histo, binl);
  //  TH1F *histo_stat = new TH1F ("histo_stat","Numero totale entries histo non tagliati", num_histo, binl);
  TH1F *histo_signal_int = new TH1F ("histo_signal_int","Integrale segnale dopo applicazione tagli / numero eventi per classe di molteplicita'", num_histo, binl);
  TH1F *histo_signal_int_pure = new TH1F ("histo_signal_int_pure","Integrale segnale dopo applicazione tagli (numero di V0)", num_histo, binl);
  //histo_mean->GetYaxis()->SetRangeUser(lim_inf[type], lim_sup[type]);
  histo_sigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_sigma->GetYaxis()->SetTitle("#sigma (GeV/c^{2})");
  histo_mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_mean->GetYaxis()->SetTitle("#mu (GeV/c^{2})");
  histo_chis->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_chis->GetYaxis()->SetTitle("#chi^{2}/dof");
  histo_SB->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_SB->GetYaxis()->SetTitle("");
  histo_SSB->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_SSB->GetYaxis()->SetTitle("");
  histo_S->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_S->GetYaxis()->SetTitle("");
  histo_Bcentral->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_Bcentral->GetYaxis()->SetTitle("");
  histo_Bside->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_Bside->GetYaxis()->SetTitle("");

  // histo_stat->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  // histo_stat->GetYaxis()->SetTitle("Counts");
  histo_signal_int->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_signal_int->GetYaxis()->SetTitle("Signal integral/#events (GeV/c^{2})");
  histo_signal_int_pure->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_signal_int_pure->GetYaxis()->SetTitle("Signal integral (GeV/c^{2})");
  histo_sigma->GetYaxis()->SetTitleOffset(1.2);
  histo_mean->GetYaxis()->SetTitleOffset(1.2);
  histo_chis->GetYaxis()->SetTitleOffset(1.2);
  histo_SB->GetYaxis()->SetTitleOffset(1.2);
  histo_SSB->GetYaxis()->SetTitleOffset(1.2);
  //  histo_stat->GetYaxis()->SetTitleOffset(1.2);
  histo_signal_int->GetYaxis()->SetTitleOffset(1.2);
  histo_signal_int_pure->GetYaxis()->SetTitleOffset(1.2);


  cout << "ciao" << endl;
  for(Int_t j=0;j<num_histo; j++){
    //isto[j]= (TH1F*)hMassvsPt->ProjectionX("isto_"+tipo[type]+"_" +int_pt[j],hMassvsPt->GetYaxis()->FindBin(binl[j]+0.001),hMassvsPt->GetYaxis()->FindBin(binl[j+1]-0.001));  
    isto_tagli[j]= (TH1F*)hMassvsPt_tagli->ProjectionX("isto_tagli_"+tipo[type]+"_" +int_pt[j],hMassvsPt_tagli->GetYaxis()->FindBin(binl[j]+0.001),hMassvsPt_tagli->GetYaxis()->FindBin(binl[j+1]-0.001));  
    isto_tagli[j]->Rebin(rebin[type][molt][j]);
    //if(j<=4) isto_tagli[j]->Rebin(2);
    //else  isto_tagli[j]->Rebin(4);
    // isto[j]->SetTitle("Distribuzione di massa invariante nell'intervallo di p_{T} ["+ int_pt[j] + "] GeV/c e di molteplicita ["+ multiplicity[molt]+") ("+tipo[type]+")");
    // isto_tagli[j]->SetTitle("Invariant mass distribution nell'intervallo di p_{T} ["+ int_pt[j] + "] GeV/c e di molteplicità ["+ multiplicity[molt]+") ("+tipo[type]+")");
    //isto[j]->GetXaxis()->SetTitle(invmass_title[type]);
    //isto[j]->GetYaxis()->SetTitle("Counts");
    isto_tagli[j]->GetXaxis()->SetTitle(invmass_title[type]);
    isto_tagli[j]->GetYaxis()->SetTitle("Counts");
    isto_tagli[j]->GetYaxis()->SetTitleOffset(1.6);
  
    char fname1[num_histo]; 
    char fname2[num_histo]; 
    sprintf(fname1,"1f_%d",j);
    sprintf(fname2,"2f_%d",j);
    functions1[j] = new TF1(fname1,"gaus",min_range_signal[type],max_range_signal[type]);
    functions1[j]->SetLineColor(kRed);   
    functions1[j]->SetParameter(1, massK0);
    functions1[j]->SetParName(0, "norm");
    functions1[j]->SetParName(1, "mean");
    functions1[j]->SetParName(2, "sigma");
    functions1[j]->SetParLimits(2, 0.001,100);
    functions1[j]->SetParLimits(0, 0,1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));
    functions2[j] = new TF1(fname2,"gaus",min_range_signal[type],max_range_signal[type]);
    functions2[j]->SetLineColor(kMagenta);   
    functions2[j]->SetParameter(1, massK0);
    functions2[j]->SetParName(0, "norm");
    functions2[j]->SetParName(1, "mean");
    functions2[j]->SetParName(2, "sigma");
    functions2[j]->SetParLimits(2, 0.001,100);
    functions2[j]->SetParLimits(0,0, 1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));
        
    char fnamebkg1[num_histo];
    sprintf(fnamebkg1,"f%d",j);
    bkg1[j] = new TF1(fnamebkg1,"pol1",functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
    bkg1[j]->SetLineColor(kBlue);   
    bkg1[j]->SetParameter(1, 0);
    char fnamebkg2[num_histo]; 
    sprintf(fnamebkg2,"f%d",j); 
    //    bkg2[j] = new TF1(fnamebkg2,"pol2",functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
    bkg2[j] = new TF1(fnamebkg2,"pol2",0.45, 0.55);
    bkg2[j]->SetLineColor(kBlue);
    char fnamebkg3[num_histo]; 
    sprintf(fnamebkg3,"f%d",j); 
    //  bkg3[j] = new TF1(fnamebkg3,fparab,functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2),3);
    bkg3[j] = new TF1(fnamebkg3,fparab,0.45, 0.55,3);
    bkg3[j]->SetLineColor(kBlue);
    //bkg3[j]->FixParameter(3,functions1[j]->GetParameter(1)-4*functions1[j]->GetParameter(2) );
    //bkg3[j]->FixParameter(4,functions1[j]->GetParameter(1)+4*functions1[j]->GetParameter(2) );
  
    gStyle->SetOptFit(0111);
    //histo_sigma->Sumw2();
    //histo_mean->Sumw2();

    cout<<endl;
    cout << "********* fit nel range  " << binl[j]<<"-"<<binl[j+1]<<"****************"<<endl;
    cout<<endl;
    cout<<endl;

    //cut ==0:     if((j<6 && molt<2) || (j<9 && molt ==2) || (j<8 && molt ==3) || (j<5 && j>0 && molt ==4)){
    //cut==1   (Arm)
    // if((j<6 && molt==0) || (j<=7 && molt==1)||(j<9 && molt ==2) || (j<9 && molt ==3) || (j<3 && j>0 && molt ==4)){
    //    cut ==1 Lrejection 

   Bool_t UseTwoGauss;
   //if (sysV0!=4) UseTwoGauss=((j<3 && molt==0) ||(j<4 && molt==1) || (j<4 && molt ==2) || (j<4 && molt ==3) || (j<4 && j>0 && molt ==4) || (j<0 && j>0 && molt ==5)); //no cut
   //2016lk if (sysV0!=4) UseTwoGauss=((j<2 && molt==0) ||(j<4 && molt==1) || (j<3 && molt ==2) || (j<1 && molt ==3) || (j<3 && molt ==4)|| (j<3 && molt ==5)); //no cut
   //2017kh 
   //if (sysV0!=4) UseTwoGauss=((j<5 && molt==0) ||(j<4 && molt==1) || (j<5 && molt ==2) || (j<4 && molt ==3) || (j<1 && molt ==4)|| (j<5 && molt ==5)); //no cut
 //2018d8_MCEff
   //if (sysV0!=4) UseTwoGauss=(((j<3 && j>0) && molt==0) ||(j==1 && molt==1) || (j==1 && molt ==2) || (j<0 && molt ==3) || (j<0 && molt ==4)|| (j>0 &&j<3  && molt ==5)); //no cut
    //if (sysV0!=4) UseTwoGauss=((j<0 && molt==0) ||(j<0 && molt==1) || (j<0&& molt ==2) || (j<0 && molt ==3) || (j<0 && molt ==4)|| (j<0 && molt ==5)); //no cut

   // else if (sysV0==4) UseTwoGauss=((j<6 && molt==0) || (j<=7 && molt==1)||(j<9 && molt ==2) || (j<9 && molt ==3) || (j<3 && j>0 && molt ==4) || (j<0 && j>0 && molt ==5)); //Arm
   //   else if (cut==2) UseTwoGauss=((j<6 && molt==0) || (j<=8 && molt==1)||(j<9 && molt ==2) || (j<8 && molt ==3) || (j<3 && j>1 && molt ==4)); //Lrej
   //2018b
   //  UseTwoGauss=((j<4 && molt==0) ||(j<2 && molt==1) || (j<5 && molt ==2) || (j<3 && molt ==3) || (j<1 && molt ==4)|| (j>=0 &&j<5  && molt ==5)); //no cut
   //2018mb
   //    UseTwoGauss=((j<5 && molt==0) ||(j<3 && molt==1) || (j<4 && molt ==2) || (j<2 && molt ==3) || (j<4 && molt ==4)|| (j>=0 &&j<5  && molt ==5)); //no cut
 //2018mb2016l2017hkm
  
   // UseTwoGauss=((j<5 && molt==0) ||(j<4 && j>0 && molt==1) || (j<5 && molt ==2) || (j<4 && molt ==3) || (j<4 && molt ==4)|| (j>=0 &&j<5  && molt ==5)); //no cut
   //2016k
   // UseTwoGauss=((j<4 && molt==0) ||(j<4 && j>=0 && molt==1) || (j<4&& molt ==2) || (j<4&& molt ==3) || (j<4 && j!=1 && molt ==4)|| (j>=0 &&j<5  && molt ==5)); //no cut
   //2018d8_DCACorrFinal
   // UseTwoGauss=((j<4 && molt==0) ||(j<3 && j>=0 && molt==1) || (j==1&& molt ==2) || (j==1&& molt ==3) || (j<0  && molt ==4)|| (j>=0 &&j<5  && molt ==5)); //no cut
 
   //2016k con 7 bin pT
   if(year=="2016k"){
   UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2&& molt ==3) || (j<6 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
   if(sysV0==1)   UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2&& molt ==3) || (j<5 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==3)     UseTwoGauss=((j<5 && molt==0) ||(j<6 && j>=0 &&j!=3 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2&& molt ==3) || (j<4 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==4)     UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2&& molt ==3) || (j<5 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==5)     UseTwoGauss=((j<6 && j!=3 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& j!=2 &&molt ==2) || (j<6 &&j!=2&&j!=3&& molt ==3) || (j<5 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==6)   UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2&&j!=4&& molt ==3) || (j<6 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
   }


   if(year=="2016k_onlyTriggerWithHighestPt"){
   UseTwoGauss=((j<7 && molt==0) ||(j<7 && molt==1) || (j<7&& molt ==2) || (j<7 && molt ==3) || (j<5 && j!=3 && molt ==4)|| (j<6  && molt ==5));
   if(sysV0==1)   UseTwoGauss=((j<7 && molt==0) ||(j<7 &&j!=5 && molt==1) || (j<7&& molt ==2) || (j<7 && molt ==3) || (j!=4 &&j!=3 && molt ==4)|| (j<7  && molt ==5));
   if (sysV0==2) UseTwoGauss=((j<7 && molt==0) ||(j<7 && molt==1) || (j<7&& molt ==2) || (j<7 && molt ==3) || (j<5 &&j!=3&& molt ==4)|| (j<6  && molt ==5));
   if(sysV0==3)     UseTwoGauss=((j<7&& molt==0) ||(j<7 && molt==1) || (j<7&& molt ==2) || (j<7 &&molt ==3) || (j<6 &&j!=3 && molt ==4)|| (j<7  && molt ==5));
  if(sysV0==4)     UseTwoGauss=((j<7 && molt==0) ||(j<7 && molt==1) || (j<7&& molt ==2) || (j<7 &&j!=2&& molt ==3) || (j<5 &&j!=1 &&j!=3&& molt ==4)|| (j<7 && molt ==5));
  if(sysV0==5)     UseTwoGauss=((j!=5 && j!=3 &&j!=2&& molt==0) ||(j<7 && j!=2 && molt==1) || (j<7&& j!=5 &&molt ==2) || (j<7 &&j!=2&& molt ==3) || ((j==5 ||j ==0 ||j==2) && molt ==4)|| (j!=2 && j!=4 &&j<7  && molt ==5));
  if(sysV0==6)   UseTwoGauss=((j<7 && molt==0) ||(j<7 && j!=3 && molt==1) || (j<7&& molt ==2) || (j<7 && molt ==3) || (j<3 && j!=1 && molt ==4)|| (j<7  && molt ==5));
  //UseTwoGauss=kTRUE;
   }

  if(year=="2018f1_extra_onlyTriggerWithHighestPt"){
   UseTwoGauss=((j<7 && molt==0) ||(j<7 && molt==1) || (j<6&& molt ==2) || (j<6&&j!=4 && molt ==3) || (j<5 && j!=2 && molt ==4)|| (j<7  && molt ==5));
   if(sysV0==1)   UseTwoGauss=((j<7 && molt==0) ||(j<7  && molt==1) || (j<7&& molt ==2) || (j<6 && molt ==3) || (j<5 &&j!=2 && molt ==4)|| (j<7  && molt ==5));
   if (sysV0==2) UseTwoGauss=((j<7 && molt==0) ||(j<7 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=4 && molt ==3) || (j<5 &&j!=2&& molt ==4)|| (j<7  && molt ==5));
   if(sysV0==3)     UseTwoGauss=((j<7&& molt==0) ||(j<7 && molt==1) || (j<6&& molt ==2) || (j<4 &&molt ==3) || (j<6 &&j!=4 && molt ==4)|| (j<7  && molt ==5));
  if(sysV0==4)     UseTwoGauss=((j<7 && molt==0) ||(j<7 && molt==1) || (j<6&& molt ==2) || (j<4 && molt ==3) || (j<5 && molt ==4)|| (j<7 && molt ==5));
  if(sysV0==5)     UseTwoGauss=((j<6 && molt==0) ||(j!=4 && j!=3 && molt==1) || (j<7&& j!=6 &&molt ==2) || (j!=3 &&j!=4&&j!=6 && molt ==3) || (j!=6 && j!=2 && molt ==4)|| (j<7  && molt ==5));
  if(sysV0==6)   UseTwoGauss=((j<7 && molt==0) ||(j<7  && molt==1) || (j<6&& molt ==2) || (j<7 && molt ==3) || (j<6 && j!=1 && molt ==4)|| (j<7  && molt ==5));
  //UseTwoGauss=kTRUE;
   }

   if(year=="2016k" && PtTrigMin==4){
     UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || ((j==0 || j==5)&&molt ==3) || (j<3 &&j!=1 && j!=3 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
   if(sysV0==1)   UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=4 &&j!=2&& molt ==3) || (j<4 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==3)     UseTwoGauss=((j<5 && molt==0) ||(j<6 && j>=0 &&j!=3 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=3 && j!=4 &&j!=2&& molt ==3) || (j<4 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==4)     UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=1 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2 && j!=4 && molt ==3) || (j<5 &&j!=4 && j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==5)     UseTwoGauss=((j<6 && j!=3 && j!=2 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& j!=2 &&molt ==2) || (j<6 &&j!=2&&j!=3&&j!=1 && molt ==3) || (j<3 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
  if(sysV0==6)   UseTwoGauss=((j<6 && molt==0) ||(j<6 && j>=0 && molt==1) || (j<6&& molt ==2) || (j<6 &&j!=2&&j!=4&& molt ==3) || (j<6 &&j!=1 && molt ==4)|| (j>=0 &&j<7  && molt ==5));
   }


   if(year=="2017h"){
     UseTwoGauss=kTRUE;
     if(sysV0==0) UseTwoGauss=((j<6 && molt==0) || (j<7 && molt==1) || (j<6&& molt ==2) || (j<5 && molt ==3) || (j<6 &&j!=0 && molt ==4) ||  (j<7 && molt==5));
     if(sysV0==1) UseTwoGauss=((j<6 && molt==0) || (j<6 && molt==1) || (j<6&& molt ==2) || (j<5 && molt ==3) || (j<6  && molt ==4) ||  (j<7 && molt==5));
     if(sysV0==2) UseTwoGauss=((j<6 && molt==0) || (j<7 && molt==1) || (j<6&& molt ==2) || (j<5 && molt ==3) || (j<6 && j!=0  && molt ==4) ||  (j<7 && j!=1 && j!=2 && molt==5));
     if(sysV0==3) UseTwoGauss=((j<7 && molt==0) || (j<6 && molt==1) || (j<6&& molt ==2) || (j!=5 && molt ==3) || (j<6 && j!=0 && j!=4  && molt ==4) ||  (j<7  && molt==5));
     if(sysV0==4) UseTwoGauss=((j<6 && molt==0) || (j<6 && molt==1) || (j<7&& molt ==2) || (j!=5 && molt ==3) || (j<6  && molt ==4) ||  (j<7  && molt==5));
     if(sysV0==5) UseTwoGauss=((j<7 && molt==0) || (j<7 && j!=2 && molt==1) || (j<6&& molt ==2) || (j!=5 && j!=4 && molt ==3) || (j<6 && j!=0  && molt ==4) ||  (j<6 && j!=4  && molt==5));
     if(sysV0==6) UseTwoGauss=((j<6 && molt==0) || (j<7 && molt==1) || (j<6&& molt ==2) || (j<5 && molt ==3) || (j<6 && j!=0  && molt ==4) ||  (j<7 && molt==5));   
   }

   if(year=="2018f1_extra"){
     UseTwoGauss=kTRUE;
     if(sysV0==0) UseTwoGauss=((j<7 && molt==0) || (j<7 && molt==1) || (j<6&& molt ==2) || (j<6 && j!=4 && molt ==3) || (j<4 && molt ==4) ||  (j<7 && molt==5));
   if(sysV0==1) UseTwoGauss=((j<7 && molt==0) || (j<6 && molt==1) || (j<6&& molt ==2) || (j<6 && molt ==3) || (j<7  && molt ==4) ||  (j<7 && molt==5));
     if(sysV0==2) UseTwoGauss=((j<6 && molt==0) || (j<5 && molt==1) || (j<6&& molt ==2) || (j<6 && molt ==3) || (j<7 && molt ==4) ||  (j<7 && molt==5));
     if(sysV0==3) UseTwoGauss=((j<6 && molt==0) || (j<6 && molt==1) || (j<6&&j!=1&& molt ==2) || (j!=4 && j<6&& molt ==3) || (j<4  && molt ==4) ||  (j<7  && molt==5));
     if(sysV0==4) UseTwoGauss=((j<7 && molt==0) || (j<7 && molt==1) || (j<6&& molt ==2) || (j<6&& molt ==3) || (j<4  && molt ==4) ||  (j<7  && molt==5));
     if(sysV0==5) UseTwoGauss=((j<7 && molt==0) || (j<6 && molt==1) || (j<6&& molt ==2) || (j<6 && molt ==3) || (j<4 && molt ==4) ||  (j!=2  && molt==5));
     if(sysV0==6) UseTwoGauss=((j<7 && molt==0) || (j<7 && molt==1) || (j<6&& molt ==2) || (j<6 && molt ==3) || (j<4 && molt ==4) ||  (j<7 && molt==5));   
   }

    canvas[j]= new TCanvas(Form("canvas_v%i",j),Form("canvas_v%i",j), 800,600);
    legend[j]=new TLegend(0.6,0.6,0.9,0.9);
    if(UseTwoGauss){
      char fnametotal[num_histo]; 
      sprintf(fnametotal,"f%d",j); 
      //    total[j] = new TF1(fnametotal,"gaus(0)+gaus(3)+pol2(6)",functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));  //va cambiato se cambio funzione di bkg
      total[j] = new TF1(fnametotal,"gaus(0)+gaus(3)+pol2(6)",0.45, 0.55);
    total[j]->SetLineColor(7); //7   
    total[j]->SetParName(0, "norm");
    total[j]->SetParName(1, "mean");
    total[j]->SetParName(2, "sigma");
    total[j]->SetParName(3, "norm2");
    total[j]->SetParName(4, "mean2");
    total[j]->SetParName(5, "sigma2");

    //***********limits to obtain a good fit*******************

    total[j]->SetParLimits(0, 0.08*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
    //total[j]->SetParLimits(5, 0.001,10);
    // total[j]->SetParLimits(4, total[j]->GetParameter(1)-1*total[j]->GetParameter(2), total[j]->GetParameter(1)+1*total[j]->GetParameter(2) );
    total[j]->SetParLimits(1, 0.494, 0.501);
    total[j]->SetParLimits(2, 0.002,10); //it was 0.001

    total[j]->SetParLimits(3, 0.08*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  //maximum was wothout 0.3
    total[j]->SetParLimits(4, 0.494, 0.501);
    //     total[j]->SetParLimits(5, 0.003, 1);
   total[j]->SetParLimits(5, total[j]->GetParameter(2), 1); //old one
    //    total[j]->SetParLimits(5, 0.002, 1); //new
     //    if (molt==1)  total[j]->SetParLimits(5, 0.001, 1);

    //****************************************

    isto_tagli[j]-> Fit(functions1[j], "R0"); //non ho ben capito l;a differenza tra 0 e N
    isto_tagli[j]-> Fit(functions2[j], "R0");     
    //   bkg2[j]->SetRange(functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
      bkg2[j]->SetRange(0.45, 0.55);
      //bkg3[j]->SetRange(functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
      bkg3[j]->SetRange(0.45, 0.55);
      //    total[j]->SetRange(functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
      total[j]->SetRange(0.45, 0.55);
    isto_tagli[j]-> Fit(bkg3[j], "RB0");
    functions1[j]->GetParameters(&par[j][0]);
    functions2[j]->GetParameters(&par[j][3]);
    bkg3[j]->GetParameters(&par[j][6]);
    total[j]->SetParameters(par[j]);
    fFitResultPtr0[j] =  isto_tagli[j]->Fit(total[j],"SR");//per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
    //isto_tagli[j]-> Fit(total[j], "R");
    //qui faccio cose non ben chiare ma mi pare che funzioni
    totalbis[j] = (TF1*)total[j]->Clone();
    //fFitResultPtr1[j] =  isto_tagli[j]->Fit(totalbis[j],"SR");
    fFitResultPtr1[j] = fFitResultPtr0[j] ;
    functions1[j]->FixParameter(0,total[j]->GetParameter(0));
    functions1[j]->FixParameter(1,total[j]->GetParameter(1));
    functions1[j]->FixParameter(2,total[j]->GetParameter(2));
    functions2[j]->FixParameter(0,total[j]->GetParameter(3));
    functions2[j]->FixParameter(1,total[j]->GetParameter(4));
    functions2[j]->FixParameter(2,total[j]->GetParameter(5));
    bkg2[j]->FixParameter(0,total[j]->GetParameter(6));
    bkg2[j]->FixParameter(1,total[j]->GetParameter(7));
    bkg2[j]->FixParameter(2,total[j]->GetParameter(8));
    bkg3[j]->FixParameter(0,total[j]->GetParameter(6));
    bkg3[j]->FixParameter(1,total[j]->GetParameter(7));
    bkg3[j]->FixParameter(2,total[j]->GetParameter(8));
    //    isto_tagli[j]-> Fit(functions1[j], "RB+");
    isto_tagli[j]-> Fit(functions1[j], "RB+");
    isto_tagli[j]-> Fit(functions2[j], "RB+");
    //isto_tagli[j]-> Fit(bkg3[j], "RB+");
    //    isto_tagli[j]-> Fit(bkg2[j], "B+", "", functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
    isto_tagli[j]-> Fit(bkg2[j], "BR+");

    //isto[j]-> Write();
    //    isto_tagli[j]->GetXaxis()->SetRangeUser(min_histo[type],max_histo[type]);
    isto_tagli[j]->GetXaxis()->SetRangeUser(0.45, 0.55);
    cout << "qua ok" << endl;
    canvas[j]->cd();
    isto_tagli[j]->Draw("");
      legend[j]->AddEntry(total[j],"Total", "l");
      legend[j]->AddEntry(functions1[j],"Gaussian", "l");
      legend[j]->AddEntry(functions2[j],"Gaussian", "l");
      legend[j]->AddEntry(bkg2[j],"2nd degree polynomial","l");
    legend[j]->Draw("same");
    cout << "ho disegnato legenda" << endl;
    f->WriteTObject(canvas[j]);
    canvas[j]->Close();
    cout << " ho salvato canvas " << endl;
    isto_tagli[j]-> Write();

    canv[molt]->cd(j+1);
    isto_tagli[j]->DrawCopy();

    cout << "\n\n\n\n ho salvato canvas con istogramma per m "<< molt <<" e pt " << j << " \n\n\n" << endl;
    isto_tagli[j]-> DrawCopy();
    TMatrixDSym cov =   fFitResultPtr0[j]->GetCovarianceMatrix();
    Double_t cov_mean = cov[1][4];
    Double_t cov_sigma = cov[2][5];
    cout << pow(total[j]->GetParError(1),2) << " " << pow(total[j]->GetParError(4),2) << " " << cov_mean << endl;
    cout << pow(total[j]->GetParError(2),2) << " " << pow(total[j]->GetParError(5),2) << " " << cov_mean << endl;
    cout << cov_mean << " " << cov_sigma << endl;
    total[j]->FixParameter(6,0);
    total[j]->FixParameter(7,0);
    total[j]->FixParameter(8,0);
    totalbis[j]->FixParameter(0,0);
    totalbis[j]->FixParameter(1,0);
    totalbis[j]->FixParameter(2,0);
    totalbis[j]->FixParameter(3,0);
    totalbis[j]->FixParameter(4,0);
    totalbis[j]->FixParameter(5,0);

    mean[j]=(functions1[j]->GetParameter(1) + functions2[j]->GetParameter(1))/2;
    errmean[j]=sqrt(pow(total[j]->GetParError(1),2) +pow(total[j]->GetParError(4),2)+2*cov_mean)/2;
    sigma[j]=(functions1[j]->GetParameter(2) + functions2[j]->GetParameter(2))/2;
    errsigma[j]=sqrt(pow(total[j]->GetParError(2),2) +pow(total[j]->GetParError(5),2)+2*cov_sigma)/2;
  
  }

   else{
  char fnametotal[num_histo]; 
    sprintf(fnametotal,"f%d",j); 
    total[j] = new TF1(fnametotal,"gaus(0)+pol2(3)",functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));  //va cambiato se cambio funzione di bkg
    total[j]->SetLineColor(7); //7   
    total[j]->SetParName(0, "norm");
    total[j]->SetParName(1, "mean");
    total[j]->SetParName(2, "sigma");
    total[j]->SetParLimits(2, 0.001,10);
    total[j]->SetParLimits(0, 0,1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
    isto_tagli[j]-> Fit(functions1[j], "R0"); //non ho ben capito l;a differenza tra 0 e N
    //        bkg2[j]->SetRange(functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
           bkg2[j]->SetRange(0.45, 0.55);
	   //  bkg3[j]->SetRange(functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
	   bkg3[j]->SetRange(0.45, 0.55);
    total[j]->SetRange(functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2));
    isto_tagli[j]-> Fit(bkg3[j], "RB0");
    functions1[j]->GetParameters(&par[j][0]);
     bkg3[j]->GetParameters(&par[j][3]);
    total[j]->SetParameters(par[j]);
    fFitResultPtr0[j] =  isto_tagli[j]->Fit(total[j],"SR");//per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
    //isto_tagli[j]-> Fit(total[j], "R");
    //qui faccio cose non ben chiare ma mi pare che funzioni
    totalbis[j] = (TF1*)total[j]->Clone();
    //fFitResultPtr1[j] =  isto_tagli[j]->Fit(totalbis[j],"SR");
    fFitResultPtr1[j] = fFitResultPtr0[j] ;

    functions1[j]->FixParameter(0,total[j]->GetParameter(0));
    functions1[j]->FixParameter(1,total[j]->GetParameter(1));
    functions1[j]->FixParameter(2,total[j]->GetParameter(2));
    bkg2[j]->FixParameter(0,total[j]->GetParameter(3));
    bkg2[j]->FixParameter(1,total[j]->GetParameter(4));
    bkg2[j]->FixParameter(2,total[j]->GetParameter(5));
    bkg3[j]->FixParameter(0,total[j]->GetParameter(3));
    bkg3[j]->FixParameter(1,total[j]->GetParameter(4));
    bkg3[j]->FixParameter(2,total[j]->GetParameter(5));
    isto_tagli[j]-> Fit(functions1[j], "RB+");
    //isto_tagli[j]-> Fit(bkg3[j], "RB+");
    //    isto_tagli[j]-> Fit(bkg2[j], "B+", "", functions1[j]->GetParameter(1)-range_bkg[type][molt][j]*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+range_bkg[type][molt][j]*functions1[j]->GetParameter(2)); 
    isto_tagli[j]-> Fit(bkg2[j], "RB+"); 
    //isto[j]-> Write();
    isto_tagli[j]->GetXaxis()->SetRangeUser(min_histo[type],max_histo[type]);
    canvas[j]->cd();
    isto_tagli[j]->Draw("");
      legend[j]->AddEntry(total[j],"Total", "l");
      legend[j]->AddEntry(functions1[j],"Gaussian", "l");
      legend[j]->AddEntry(functions2[j],"Gaussian", "l");
      legend[j]->AddEntry(bkg2[j],"2nd degree polynomial","l");
    legend[j]->Draw("same");
    cout << "ho disegnato legenda" << endl;
    f->WriteTObject(canvas[j]);
    canvas[j]->Close();
    isto_tagli[j]-> Write();

    canv[molt]->cd(j+1);
    isto_tagli[j]->DrawCopy();

    total[j]->FixParameter(3,0);
    total[j]->FixParameter(4,0);
    total[j]->FixParameter(5,0);
    totalbis[j]->FixParameter(0,0);
    totalbis[j]->FixParameter(1,0);
    totalbis[j]->FixParameter(2,0);
    mean[j]=functions1[j]->GetParameter(1);
    errmean[j]=total[j]->GetParError(1);
    sigma[j]=functions1[j]->GetParameter(2);
    errsigma[j]=total[j]->GetParError(2);

    }
 
    //cout << "parametri del fit" << fFitResultPtr0[j]->GetParams() << endl;
    cout << "hello" << endl;
    for(Int_t l=isto_tagli[j]->GetXaxis()->FindBin(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2)); l<=isto_tagli[j]->GetXaxis()->FindBin(functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2)); l++){
      entries_range[j]+=isto_tagli[j]->GetBinContent(l);
    }

    ////////////////////////////////////
    s1[j]=functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2));
    s2[j]=functions2[j]->Integral(functions2[j]->GetParameter(1)-sigmacentral*functions2[j]->GetParameter(2),functions2[j]->GetParameter(1)+sigmacentral*functions2[j]->GetParameter(2));
    st[j]=total[j]->Integral(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j]);
    // cout <<"con totale " <<  st[j]<< endl;
    //   st[j]=s1[j]+s2[j];
    // cout <<"con somma " <<  st[j]<< endl;

    tot[j]=entries_range[j]*isto_tagli[j]->GetBinWidth(1);

    b[j]=bkg2[j]->Integral(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j]);
    errb[j]=bkg2[j]->IntegralError(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j]);
    //    bside[j]=bkg2[j]->Integral(mean[j]-nsigmamax*sigma[j],mean[j]-nsigmamin*sigma[j]) + bkg2[j]->Integral(mean[j]+nsigmamin*sigma[j], mean[j]+nsigmamax*sigma[j]);
    bside[j]=bkg2[j]->Integral(0.45,mean[j]-nsigmamin*sigma[j]) + bkg2[j]->Integral(mean[j]+nsigmamin*sigma[j], 0.55);
    //    errbside[j]=bkg2[j]->IntegralError(mean[j]-nsigmamax*sigma[j], mean[j]-nsigmamin*sigma[j]) + bkg2[j]->IntegralError(mean[j]+nsigmamin*sigma[j], mean[j]+nsigmamax*sigma[j]);
    errbside[j]=bkg2[j]->IntegralError(0.45, mean[j]-nsigmamin*sigma[j]) + bkg2[j]->IntegralError(mean[j]+nsigmamin*sigma[j], 0.55);
    bin_content0[j]=st[j]; //err0 e' l'errore associato
    bin_contents3[j]= tot[j]-b[j];//errs3
    bin_contentSSB1[j]= st[j]/(st[j]+b[j]);//errSSB1
    bin_contentSSB2[j]= (tot[j]-b[j])/tot[j];//errSSB2

    if (b[j]>0)  {
      bin_content1[j]=st[j]/b[j];//err1
      bin_content2[j]= (tot[j]-b[j])/b[j];//err3
    }
    else{
      bin_content1[j]=0;
      bin_content2[j]=0;
    }

    sigmas1[j]=total[j]->IntegralError(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j],fFitResultPtr0[j] ->GetParams(),(fFitResultPtr0[j]->GetCovarianceMatrix()).GetMatrixArray());
    sigmab1[j]=totalbis[j]->IntegralError(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j],fFitResultPtr1[j] ->GetParams(),(fFitResultPtr1[j]->GetCovarianceMatrix()).GetMatrixArray());
    err0[j]=sigmas1[j];
    err1[j]=   sqrt(pow(sigmas1[j],2)/pow(b[j],2) + pow(sigmab1[j],2)*pow(st[j],2)/pow(b[j],4));
    //fcredo sbagliato    err2[j]=sqrt( (pow(isto_tagli[j]->GetBinWidth(1),2)*entries_range[j] + pow(sigmab1[j],2) )/pow(b[j],2) +  pow(sigmab1[j],2)* pow((tot[j]-b[j]),2)/pow(b[j],4)) ;
    err2[j]=sqrt( (pow(isto_tagli[j]->GetBinWidth(1),2)*entries_range[j] )/pow(b[j],2) + pow(sigmab1[j],2)*pow(tot[j],2)/pow(b[j],4)) ;
    errs3[j]=sqrt( pow(isto_tagli[j]->GetBinWidth(1),2)*entries_range[j] + pow(sigmab1[j],2));
    //errori per histo_signal_int
    err3[j]=sqrt(pow(sigmas1[j],2)/pow(isto_tagli[j]->GetBinWidth(1),2)/pow(events,2)   +  pow(st[j],2)/pow(isto_tagli[j]->GetBinWidth(1),2)/pow(events,3) );
    err4[j]=sqrt(pow(errs3[j],2)/pow(isto_tagli[j]->GetBinWidth(1),2)/pow(events,2)   +  pow(tot[j]-b[j],2)/pow(isto_tagli[j]->GetBinWidth(1),2)/pow(events,3) );
    errSSB1[j]=(pow(st[j],2)*pow(sigmab1[j],2)+pow(b[j],2)*pow(sigmas1[j],2))/pow(st[j]+b[j],4);
    errSSB2[j]=pow(b[j],2)/pow(tot[j],3)-pow(sigmab1[j],2)/pow(tot[j],2);

    cout << "hello" << endl;

    Bool_t ShowInHisto=kFALSE;
    if(sysV0!=4) ShowInHisto=((molt==0 && j<9) || (molt==1 && j<11)|| (molt==2 && j<11)|| (molt==3 && j<11)||(molt==4 && j<7) ||(molt==5) ); //nocut da rivedere
    else if(sysV0==4) ShowInHisto=((molt==0 && j<10) || (molt==1 && j<10)|| (molt==2 && j<11)|| (molt==3 && j<9)||(molt==4 && j<7) || (molt==5)); //Arm
    //    else if(cut==2) ShowInHisto=((molt==0 && j<9) || (molt==1 && j<9)|| (molt==2 && j<11)|| (molt==3 && j<9)||(molt==4 && j<7)); //Lrej

    if(ShowInHisto){
    if(mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
      histo_sigma->SetBinContent(j+1,sigma[j]);
      histo_sigma->SetBinError(j+1,errsigma[j]);
      histo_mean->SetBinContent(j+1,mean[j]);
      histo_mean->SetBinError(j+1,errmean[j]);

    } else {
      histo_sigma->SetBinContent(j+1,0);
      histo_sigma->SetBinError(j+1,0);
      //histo_mean->SetBinContent(j+1,0);
      //histo_mean->SetBinError(j+1,0);
 
    }
    if(total[j]->GetNDF() != 0){
      histo_chis->SetBinContent(j+1,total[j]->GetChisquare()/total[j]->GetNDF());
    } else {
      histo_chis->SetBinContent(j+1,0);
    }
    // histo_stat->SetBinContent(j+1,isto[j]->GetEntries());
    cout << "hello" << endl;

    for(Int_t l=isto_tagli[j]->GetXaxis()->FindBin(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2)); l<=isto_tagli[j]->GetXaxis()->FindBin(functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2)); l++){
      entries_range[j]+=isto_tagli[j]->GetBinContent(l);
    }


    cout << "hello" << endl;
    if((total[j]->GetChisquare()/total[j]->GetNDF())<=3. && mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && bin_content1[j]>0 && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
      histo_SB->SetBinContent(j+1,bin_content1[j]);
      histo_SB->SetBinError(j+1, err1[j]);
      histo_SSB->SetBinContent(j+1,bin_contentSSB1[j]);
      histo_SSB->SetBinError(j+1, errSSB1[j]);
      histo_S->SetBinContent(j+1,bin_content0[j]/isto_tagli[j]->GetBinWidth(1));
      histo_S->SetBinError(j+1, err0[j]/isto_tagli[j]->GetBinWidth(1));
      histo_Bcentral->SetBinContent(j+1,b[j]);
      histo_Bcentral->SetBinError(j+1, errb[j]);
      histo_Bside->SetBinContent(j+1,bside[j]);
      histo_Bside->SetBinError(j+1, errbside[j]);
      histo_signal_int->SetBinContent(j+1,st[j]/isto_tagli[j]->GetBinWidth(1)/events);
      histo_signal_int->SetBinError(j+1,err3[j]);
      histo_signal_int_pure->SetBinContent(j+1,st[j]/isto_tagli[j]->GetBinWidth(1));
      histo_signal_int_pure->SetBinError(j+1,sigmas1[j]/isto_tagli[j]->GetBinWidth(1));

      cout << "caso A"<< endl;

    } else if((total[j]->GetChisquare()/total[j]->GetNDF())>3. && mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && bin_content2[j]>0 && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
      histo_SB->SetBinContent(j+1,bin_content2[j]);
      histo_SB->SetBinError(j+1,err2[j]);
      histo_SSB->SetBinContent(j+1,bin_contentSSB1[j]);
      histo_SSB->SetBinError(j+1, errSSB1[j]);
      histo_S->SetBinContent(j+1,bin_contents3[j]/isto_tagli[j]->GetBinWidth(1));
      histo_S->SetBinError(j+1,errs3[j]/isto_tagli[j]->GetBinWidth(1));
      histo_Bcentral->SetBinContent(j+1,b[j]);
      histo_Bcentral->SetBinError(j+1, errb[j]);
      histo_Bside->SetBinContent(j+1,bside[j]);
      histo_Bside->SetBinError(j+1, errbside[j]);
      histo_signal_int->SetBinContent(j+1,(tot[j]-b[j])/isto_tagli[j]->GetBinWidth(1)/events);
      histo_signal_int->SetBinError(j+1,err4[j]);
      histo_signal_int_pure->SetBinContent(j+1,(tot[j]-b[j])/isto_tagli[j]->GetBinWidth(1));
      histo_signal_int_pure->SetBinError(j+1,errs3[j]/isto_tagli[j]->GetBinWidth(1));

      cout << "caso B"<< endl;
      cout <<total[j]->GetParameter(1) << endl; 
      cout << functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))<< endl; 
      cout << bin_content2[j] << endl;
      cout << entries_range[j]*isto_tagli[j]->GetBinWidth(1)/(-functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))+entries_range[j]*isto_tagli[j]->GetBinWidth(1)) << endl;

    }  else {
      cout << " S su B varra' zero perche  condizioni non sono soddisfatte " << endl; 
      cout << "reduced chi " << total[j]->GetChisquare()/total[j]->GetNDF()<< "\nmean " << mean[j] << " > " << lim_inf_mean[type] << " < " <<lim_sup_mean[type] <<"\n S/B " <<bin_content2[j]<< "\n sigma " << sigma[j] << " < " << lim_sup_sigma[type] << "\n err sigma " << errsigma[j] << " < " << lim_sup_errsigma[type] << "\n err mean " << errmean[j] << " < " << lim_sup_errmean[type] << endl;
      histo_SB->SetBinContent(j+1,0);
      histo_SB->SetBinError(j+1,0);
      histo_SSB->SetBinContent(j+1,0);
      histo_SSB->SetBinError(j+1,0);
      histo_S->SetBinContent(j+1,0);
      histo_S->SetBinError(j+1,0);
      histo_Bcentral->SetBinContent(j+1,0);
      histo_Bcentral->SetBinError(j+1, 0);
      histo_Bside->SetBinContent(j+1,0);
      histo_Bside->SetBinError(j+1, 0);
      histo_signal_int->SetBinContent(j+1,0);
      histo_signal_int->SetBinError(j+1,0);
      histo_signal_int_pure->SetBinContent(j+1,0);
      histo_signal_int_pure->SetBinError(j+1,0);

      cout << "caso C " << endl;
    }
    
    cout << "****************"<<endl;
    cout << total[j]->GetChisquare()<<endl;
    cout << total[j]->GetNDF()<<endl;
    cout << total[j]->GetChisquare()/total[j]->GetNDF() << endl;
    cout <<functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2)) << endl;
    cout << bkg2[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))<< endl;
    
    /*
      cout<< "-----------------" << endl;
      cout << entries_range[j]*isto_tagli[j]->GetBinWidth(1) << endl;
      cout << (-functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))+entries_range[j]*isto_tagli[j]->GetBinWidth(1)) << endl;
      cout << bin_content2[j] << endl;
      cout << entries_range[j]*isto_tagli[j]->GetBinWidth(1)/(-functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))+entries_range[j]*isto_tagli[j]->GetBinWidth(1)) << endl;
    */
    cout<<"---------------"<<j<<endl;
    }  
  }

  for(Int_t j=0;j<num_histo; j++){
    if(bin_content1[j]<0 || bin_content2[j]<0){
      cout << "*************** ho ottenuto alcuni valori di S/B negativi**************" << endl;
    }
    //cout << "  entries nel range del segnale  " << entries_range[j] << endl;
    //cout << "  integrale nel range del segnale/ampiezza bin   "<< (functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2)))/isto_tagli[j]->GetBinWidth(1) << endl;
  }

  //  event_multiplicity->Write();
  histo_sigma->Write();
  histo_mean->Write();
  histo_chis->Write();
  histo_SB->Write();
  histo_SSB->Write();
  histo_S->Write();
  histo_Bcentral->Write();
  histo_Bside->Write();

  //  histo_stat->Write();
  // histo_signal_int->Write();
  // histo_signal_int_pure->Write();
  f->Close();
 
  
  cout << endl;
  cout << "*******************" << endl;
  cout << "Pt Trigger minimo " << PtTrigMin<< endl;
  cout << "sto analizzando l'intervallo di molteplicita' " << mult<< " per la particella " << tipo[type]<< endl;
  cout << "a partire dal file " << nome_file_1<< endl;
  cout <<"gli istogrammi sono nel file " <<nome_file_output[sysTrigger][sysV0] << endl;
  //  cout << "numero totale eventi "<< events_total << endl;
  cout << "*******************" << endl;
  cout <<   "run this for syst =0,1,2 (sysV0=0) and sysV0 =0,..,6" << endl;
  cout << endl;
  }
	
}
