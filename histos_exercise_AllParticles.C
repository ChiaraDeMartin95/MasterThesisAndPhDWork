#include <Riostream.h>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TLine.h>
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
  Float_t LimInf=0;
  Float_t LimSup=0;
    if (par[3]==0) {LimInf=0.474; LimSup=0.520;}
  else if (par[3]==4 || par[3]==5 || par[3]==8) {LimInf=1.310; LimSup=1.335;}
  if (reject && x[0] > LimInf && x[0] < LimSup) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t fretta(Double_t *x, Double_t *par)
{
  Float_t LimInf=0;
  Float_t LimSup=0;
  //if (par[2]==0) {LimInf=0.474; LimSup=0.520;}
  if (par[2]==0) {LimInf=0.47; LimSup=0.530;}
  else if (par[2]==4 || par[2]==5 || par[2]==8) {LimInf=1.310; LimSup=1.335;}
  if (reject && x[0] > LimInf && x[0] < LimSup) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0];
}

Double_t SetEfficiencyError(Int_t k, Int_t n){
  return sqrt(((Float_t)k+1)*((Float_t)k+2)/(n+2)/(n+3) - pow((Float_t)(k+1),2)/pow(n+2,2));
}


void histos_exercise_AllParticles(Bool_t isSignalFromIntegral=0, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Float_t PtTrigMin=3,Float_t PtTrigMax=15,  Int_t rap=0,Int_t type=0, Int_t sysTrigger=0, Int_t sysV0=0, Int_t syst=0, Double_t nsigmamax=9, TString year0="2016", TString year="2018f1_extra_onlyTriggerWithHighestPt"/*"2018f1_extra_hXi_65runs"*/, Bool_t isMC=1, Bool_t isEfficiency=1, TString path1=""){

  //isSignalFromIntegral permette di scegliere tra l'utilizzo di S = integral fit function o S = entries- integral bkg function
  //IsBkgPParab = fit con pl2 per fondo if kTRUE, fit con pol1 per fodno if kFALSE

  if (kTRUE) {
    //    cout << "controllare che la funzione di fit totale e il numero di parametri par[] siano adeguati alla scelta della funzioen di bkg (pol1 o pol2)" << endl;
    cout << "cotrollare range fretta e fparab siano adeguate a tipo di particella scelta " << endl;

    //    return
  }
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

  Double_t sigmacentral=4; //it was 3 for kaons
  Double_t nsigmamin=5; //it was four for kaons
  
  if (syst==1) {
    nsigmamin=6;
  }
  if (syst==2) {
    sigmacentral=5;
  }

  cout << "*****************************************************************************"<< endl;
  cout << "Attenzione:  I limiti superiori di errsigma e errmean dipendono dalla tipologia di taglio implementato, in modo da garantire che i vari istogrammi siano rimepiti solo se il fit viene \"bene\" " << endl;

  if (isMC && !isEfficiency) cout <<"**********************ERRORE****************"<< endl;
  const Int_t mult=5; //numero intervalli molteplicita'
  const Int_t num_tipo=10; //numero tipologie di particelle
  const Int_t num_histo=8; //numero intervalli di Pt
  const Int_t num_tagli=6; //numero di diversi tagli applicati alle K0s
  const Int_t NsysTrigger=3; 
  const Int_t NsysV0=9; 

  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};
  TString int_pt[num_histo]={"", "0.5-1", "1-1.5","1.5-2", "2-2.5","2.5-3", "3-4", "4-8"};
  TString multiplicity[mult+1]={"0-5","5-10", "10-30", "30-50", "50-100", "all"};
  Float_t YLimit[2] = {1.0, 0.5};
 
  Float_t massParticle[num_tipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString tipo[num_tipo]={"kK0s", "Lambda", "AntiLambda","LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPlus", "Xi", "Omega"};
  TString Globaltipo[num_tipo]={"kK0s", "Lambda", "AntiLambda","LambdaAntiLambda", "Xi", "Xi", "Omega", "Omega", "Xi", "Omega"};

  cout << "la particella scelta è" << tipo[type]<<endl;
  
  TString nome_TDir="MyTask";
  if (type!=0 ) nome_TDir +=Globaltipo[type];
  TString invmass[num_tipo]={"#pi^{+} #pi^{-}", "p #pi^{-}", "overline{p} #pi^{+}", "p #pi^{-} + overline{p} #pi^{+}", "#pi^{-} #Lambda", "#pi^{+} overline{#Lambda}", "K^{-} #Lambda", "K^{+} overline{#Lambda}", "#pi^{+} overline{#Lambda} + #pi^{-} #Lambda",  "K^{+} overline{#Lambda} + K^{-} #Lambda"};
  TString invmass_title[num_tipo] = {"(#pi^{+}, #pi^{-}) invariant mass (GeV/c^{2})", "massa invariante (p, #pi^{-}) (GeV/c^{2})", "massa invariante (overline{p}, #pi^{+})", "massa invariante [(p, #pi^{-}) + (overline{p}, #pi^{+})] (GeV/c^{2})", "(#pi^{-}, #Lambda) invariant mass (GeV/c^{2})", "(#pi^{+}, overline{#Lambda}) invariant mass (GeV/c^{2})", "(K^{-}, #Lambda) invariant mass (GeV/c^{2})", "(K^{+}, overline{#Lambda}) invariant mass (GeV/c^{2})", "(#pi^{+}, #Lambda) invariant mass (GeV/c^{2})", "(K^{+},#Lambda) invariant mass (GeV/c^{2})"};

  Float_t min_range_signal[num_tipo]={0.47,1.105,1.105,1.105, 1.31, 1.31, 1.66, 1.66, 1.31, 1.66}; //estremi region fit segnale (gaussiane)
  Float_t max_range_signal[num_tipo]={0.53,1.125,1.125,1.125, 1.334, 1.334, 1.68, 1.68,1.334, 1.681}; 
  Float_t min_histo[num_tipo]={0.45,1.09,1.09,1.09, 1.30, 1.30, 1.65, 1.65, 1.30, 1.65};  //estremi del range degli istogrammi
  Float_t max_histo[num_tipo]={0.55,1.14,1.14,1.14, 1.342, 1.342, 1.69, 1.69 , 1.342, 1.69};
  Float_t liminf[num_tipo]={0.495,1.1153,1.1153,1.1153, 1.30, 1.30, 1.66, 1.66, 1.30, 1.66}; //estremi regione fit del bkg e total
  Float_t limsup[num_tipo]={0.500,1.1168,1.1168,1.1168, 1.342, 1.342, 1.68, 1.68, 1.342, 1.68};  
  //era 1.30 e 1.342
  Float_t lim_inf_mean[num_tipo]={0.495,1.1153,1.1153,1.1153, 1.31, 1.31, 1.66, 1.66, 1.31, 1.66}; //come l'ho scelto?
  Float_t lim_sup_mean[num_tipo]={0.500,1.1168,1.1168,1.1168, 1.33, 1.33, 1.68, 1.68, 1.33, 1.68};
  Float_t lim_inf_sigma[num_tipo]={0};
  Float_t lim_sup_sigma[num_tipo]={0.008,0.002,0.002,0.002 , 0.008, 0.008, 0.008, 0.008, 0.008, 0.008}; //first one must be 0.008 if we want to display the 7-15 mult values between 7 and 8 pT
  Float_t lim_inf_errmean[num_tipo]={0}; 
  Float_t lim_sup_errmean[num_tipo]={10,0.0006,0.00035,0.0006, 10, 10, 10, 10 ,10,10};//loooooose
  Float_t lim_inf_errsigma[num_tipo]={0}; 
  Float_t lim_sup_errsigma[num_tipo]={10,0.0004,0.0015,0.0004 , 10, 10, 10, 10, 10, 10}; //loose

  Int_t entries_range[num_histo]={0};
  Int_t entries_range_true[num_histo]={0};
  Int_t entries_range_false[num_histo]={0};
  Int_t entries_sideband[num_histo]={0};
  Int_t entries_sideband_true[num_histo]={0};
  Int_t entries_sideband_false[num_histo]={0};
  Float_t bin_contentS1[num_histo]={0.};
  Float_t bin_contentSB1[num_histo]={0.};
  Float_t bin_contentSB2[num_histo]={0.};
  Float_t bin_contentS2[num_histo]={0.};
  Float_t bin_contentSSB1[num_histo]={0.};
  Float_t bin_contentSSB2[num_histo]={0.};
  Float_t sigmas1[num_histo]={0};
  Float_t sigmab1[num_histo]={0};
  Float_t errS2[num_histo]={0};
  Float_t errS1[num_histo]={0};
  Float_t errSB1[num_histo]={0};
  Float_t errSB2[num_histo]={0};
  Float_t errSSB1[num_histo]={0};
  Float_t errSSB2[num_histo]={0};
  Float_t s1[num_histo]={0};
  Float_t s2[num_histo]={0};
  Float_t st[num_histo]={0};
  Float_t IntegralSignalAllRange[num_histo]={0};
  Float_t tot[num_histo]={0};
  Float_t b[num_histo]={0};
  Float_t bside[num_histo]={0};
  Float_t errb[num_histo]={0};
  Float_t errbside[num_histo]={0};
  Float_t err_range_false[num_histo]={0};
  Float_t err_sideband_false[num_histo]={0};
  
  Float_t sigma[num_histo]={0};
  Float_t errsigma[num_histo]={0};
  Float_t mean[num_histo]={0};
  Float_t errmean[num_histo]={0};


  TFitResultPtr fFitResultPtr0[num_histo]; //serve per poter calcolare in seguito l'errore di un integrale 
  TFitResultPtr fFitResultPtr1[num_histo];

  Double_t multip_intervals[mult+1]={0, 5, 10, 30, 50,100};
  //  Double_t binl[num_histo+1]={0, 1,1.5, 2,2.5,3,4,8}
  Double_t binl[num_histo+1]={0,0.5, 1,1.5, 2,2.5,3,4,8};

  Int_t rebin[num_tipo][mult+1][num_histo]={1};
  Double_t range_bkg[num_tipo][mult+1][num_histo];
  for(Int_t i=0; i< num_tipo; i++){
    for(Int_t j=0; j< mult+1; j++){
      for(Int_t l=0; l< num_histo; l++){
	range_bkg[i][j][l]=nsigmamax;
      }    
    }  
  }
  //  cout<< " ciao " << endl;
  TFile *myfile; 
  TFile *myfileAnalysis; 
  TCanvas *canv[mult+1];
    Int_t     NEvents[mult+1];
    TString Srap[2] = {"_Eta0.8", "_y0.5"};
    for(Int_t molt=0; molt<mult+1; molt++){
            if (molt < mult) continue;
//      TString nome_file_1 ="FinalOutput/DATA"+year0+"/histo/AngularCorrelation"+year;
      TString nome_file_1 ="AngularCorrelation"+year; //local
      TString nome_file_output[NsysTrigger][NsysV0];
      TString nome_file_analysis;
//      nome_file_analysis="FinalOutput/AnalysisResults"+year+".root";
      nome_file_analysis="AnalysisResultsFolder/AnalysisResults"+year+".root"; //local
      //      if (isMC && isEfficiency) nome_file_analysis="FinalOutput/AnalysisResults"+year+"_MCEff.root";
      if (isMC && isEfficiency) nome_file_analysis="AnalysisResultsFolder/AnalysisResults"+year+"_MCEff.root"; //local
//      nome_file_output[sysTrigger][sysV0] ="FinalOutput/DATA"+year0+"/invmass_distribution_thesis/invmass_distribution";
      nome_file_output[sysTrigger][sysV0] ="invmass_distribution"; //local

      if(isMC && isEfficiency){
	nome_file_1+="_MCEff";
	nome_file_output[sysTrigger][sysV0]+="_MCEff";
      }

      nome_file_output[sysTrigger][sysV0]+=path1;
      nome_file_output[sysTrigger][sysV0]+="_";
      if (type!=0){ 
	nome_file_1+="_" +tipo[type];
      }

      nome_file_1+=path1;
      if (type!=0){
	nome_file_1 +=Srap[rap];
      }
      nome_file_1 +=Form("_MassDistr_SysT%i_SysV0%i_PtMin%.1f.root",sysTrigger, sysV0, PtTrigMin);
      nome_file_output[sysTrigger][sysV0] +=Form(year+"_"+tipo[type]+Srap[rap]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", molt, sysTrigger, sysV0, syst,PtTrigMin);
      
      cout << "questo è il nome del file prodotto dal task (nome_file_analysis) " << nome_file_analysis << endl;
      cout << "questo è il nome del file " << nome_file_1 << endl;
      cout << "questo è il nome del file di output " << nome_file_output[sysTrigger][sysV0] << endl;

      myfile = new TFile(nome_file_1, ""); 
      if (!myfile) { cout << "nome_file_1 is missing " << endl; return;}
      myfileAnalysis = new TFile(nome_file_analysis, ""); 
      if (!myfileAnalysis) { cout << "nome_file_analysis is missing " << endl; return;}

      TDirectoryFile *dirinput = (TDirectoryFile*)myfileAnalysis->Get(nome_TDir);
      if (!dirinput) {
	cout << " dir in nome_file_analysis is missing " << endl;
	return;
      }
      TList *listinputForNEvents = (TList*)dirinput->Get("MyOutputContainer");
      if (!listinputForNEvents) return;
      TH2D *fHistTriggervsMult         = (TH2D*)listinputForNEvents->FindObject("fHistPtMaxvsMultBefAll");
      TH1D *fHistMul= (TH1D*)fHistTriggervsMult->ProjectionY("fHistTriggervsMult_MultProj", fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMin+0.0001),fHistTriggervsMult->GetXaxis()->FindBin(PtTrigMax-0.00001) );

      //      TH1F* fHistMul = (TH1F*)listinputForNEvents->FindObject("fHist_multiplicity_EvwTrigger");
     
      if (!fHistMul) return;

      NEvents[molt] =0;
      if (molt!=mult){
      for (Int_t l=fHistMul->FindBin(multip_intervals[molt]+0.0001);l<fHistMul->FindBin(multip_intervals[molt+1]-0.0001);l++){
	NEvents[molt]+=fHistMul->GetBinContent(l);
      }
      }
      else {
      for (Int_t l=fHistMul->FindBin(0+0.0001);l<fHistMul->FindBin(100-0.0001);l++){
	NEvents[molt]+=fHistMul->GetBinContent(l);
      }
      }
      cout << "n events " << NEvents[molt]<< endl;

    canv[molt] = new TCanvas(Form("canv_molt%i", molt),"Distribuzione di massa invariante di " +invmass[type] + " dopo applicazione tagli",1600,1000);
    canv[molt]->Divide(4,2);

    Double_t events=1;
  
    TFile *f = new TFile(nome_file_output[sysTrigger][sysV0],"RECREATE");
    //  cout << "ciao " << endl;
    TH2F *hMassvsPt_tagli = (TH2F*)myfile->Get(Form("SE_hMassvsPt_"+tipo[type]+"_%i",molt));
    TH2F *hMassvsPt_tagli_true = (TH2F*)myfile->Get(Form("SE_hMassvsPt_"+tipo[type]+"_%i_true",molt)); //new for Casc

    TH1F *isto[num_histo];
    TH1F *isto_tagli[num_histo];
    TH1F *isto_tagli_true[num_histo];
    TH1F *isto_tagli_false[num_histo];
    //  cout<< " ciao " << endl;

    Double_t parTwoGaussParab[num_histo][9];  
    Double_t parTwoGaussRetta[num_histo][8];  
    Double_t parOneGaussParab[num_histo][6];  
    Double_t parOneGaussRetta[num_histo][5];  
    //  cout<< " ciao " << endl;

    TF1 **functionsFirst = new TF1*[num_histo];  
    TF1 **functionsSecond = new TF1*[num_histo]; 
    TF1 **functions1 = new TF1*[num_histo];  
    TF1 **functions2 = new TF1*[num_histo];
    //  cout<< " ciao " << endl;

    TLine *lineCentralSX[num_histo];
    TLine *lineCentralDX[num_histo];
    TLine *lineSidebandsSX[num_histo];
    TLine *lineSidebandsDX[num_histo];
    TLegend *legend[num_histo];
    TCanvas* canvas[num_histo];
    TF1 **bkg1 = new TF1*[num_histo];
    TF1 **bkg2 = new TF1*[num_histo];
    TF1 **bkg2retta = new TF1*[num_histo];
    TF1 **bkgretta = new TF1*[num_histo];
    TF1 **bkgparab = new TF1*[num_histo];
    TF1 **total= new TF1*[num_histo];
    TF1 **totalbis =new TF1*[num_histo]; //mi serve solo per calcolare errore del bkg
    TH1F *histo_sigma = new TH1F ("histo_sigma","Sigma vs Pt", num_histo, binl);
    TH1F *histo_mean = new TH1F ("histo_mean","Mean vs Pt", num_histo, binl);
    TH1F *histo_chis = new TH1F ("histo_chis","ChiSquare/dof vs Pt", num_histo, binl);
    TH1F *histo_SB = new TH1F ("histo_SB","S/B entro 3sigma vs Pt", num_histo, binl);
    TH1F *histo_SSB = new TH1F ("histo_SSB","S/(S+B) entro 3sigma vs Pt", num_histo, binl);
    TH1F *histo_S = new TH1F ("histo_S","S entro 3sigma vs Pt", num_histo, binl);
    TH1F *histo_Bcentral = new TH1F ("histo_Bcentral","B entro nsigma vs Pt", num_histo, binl);
    TH1F *histo_BsideFalse = new TH1F ("histo_BsideFalse","B in sideband from MC vs Pt", num_histo, binl);
    TH1F *histo_BcentralFalse = new TH1F ("histo_BcentralFalse","B entro nsigma from MC vs Pt", num_histo, binl);
    TH1F *histo_Bside = new TH1F ("histo_Bside","B in sideband vs Pt", num_histo, binl);
    TH1F *histo_FracStrangePeak = new TH1F ("histo_FracStrangePeak","Frazione picco bkg su totale vs Pt", num_histo, binl);  
  TH1F *histo_BRatioCentral = new TH1F ("histo_BRatioCentral","Rapporto B(integral)/B(verità MC) nsigma vs Pt", num_histo, binl);
    TH1F *histo_SRatioCentral = new TH1F ("histo_SRatioCentral","Rapporto S(integral)/S(verità MC) nsigma vs Pt", num_histo, binl);
    TH1F *histo_BRatioSide = new TH1F ("histo_BRatioSide","Rapporto B(integral)/B(verità MC) in sideband vs Pt", num_histo, binl);
    TH1F *histo_BDoubleRatio = new TH1F ("histo_BDoubleRatio","Rapporto Bside/Bcentral(integral/MCTruth) vs Pt", num_histo, binl);
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
    histo_BcentralFalse->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    histo_BcentralFalse->GetYaxis()->SetTitle("");
    histo_BsideFalse->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    histo_BsideFalse->GetYaxis()->SetTitle("");
    histo_signal_int->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    histo_signal_int->GetYaxis()->SetTitle("Signal integral/#events (GeV/c^{2})");
    histo_signal_int_pure->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    histo_signal_int_pure->GetYaxis()->SetTitle("Signal integral (GeV/c^{2})");
    histo_sigma->GetYaxis()->SetTitleOffset(1.2);
    histo_mean->GetYaxis()->SetTitleOffset(1.2);
    histo_chis->GetYaxis()->SetTitleOffset(1.2);
    histo_SB->GetYaxis()->SetTitleOffset(1.2);
    histo_SSB->GetYaxis()->SetTitleOffset(1.2);
    histo_signal_int->GetYaxis()->SetTitleOffset(1.2);
    histo_signal_int_pure->GetYaxis()->SetTitleOffset(1.2);

    //    cout << "ciao" << endl;
    for(Int_t j=1;j<num_histo; j++){
      if (type==4 || type==5 || type==8){
	liminf[type]=1.30;
	limsup[type]=1.342;
	min_histo[type]=1.30;
	max_histo[type]=1.342;
    }
      if (j >4)  {
	if (type==4 || type==5 || type==8){
	liminf[type]=1.29;
	limsup[type]=1.352;
	min_histo[type]=1.29;
	max_histo[type]=1.352;
	}
      }
      if (type==0){
	liminf[type]=0.46;//where the bkg and total fit is performed (the 2 gaussian fit is performed in the min_range_histo range
	limsup[type]=0.54;
	min_histo[type]=0.455; //where final bkg function is displayed
	max_histo[type]=0.54;
    }

      //    cout << "\n\n " << j << endl;
      //isto[j]= (TH1F*)hMassvsPt->ProjectionX("isto_"+tipo[type]+"_" +int_pt[j],hMassvsPt->GetYaxis()->FindBin(binl[j]+0.001),hMassvsPt->GetYaxis()->FindBin(binl[j+1]-0.001));  
      isto_tagli[j]= (TH1F*)hMassvsPt_tagli->ProjectionX("isto_tagli_"+tipo[type]+"_" +int_pt[j],hMassvsPt_tagli->GetYaxis()->FindBin(binl[j]+0.001),hMassvsPt_tagli->GetYaxis()->FindBin(binl[j+1]-0.001));  
      isto_tagli_true[j]= (TH1F*)hMassvsPt_tagli_true->ProjectionX("isto_tagli_true"+tipo[type]+"_" +int_pt[j],hMassvsPt_tagli->GetYaxis()->FindBin(binl[j]+0.001),hMassvsPt_tagli->GetYaxis()->FindBin(binl[j+1]-0.001));  
      isto_tagli_false[j]= (TH1F*)      isto_tagli[j]->Clone("isto_tagli_false"+tipo[type]+"_" +int_pt[j]);
      //    cout << "qui ok" << endl;
      isto_tagli[j]->Rebin(rebin[type][molt][j]);
      isto_tagli_true[j]->Rebin(rebin[type][molt][j]);
      isto_tagli_false[j]->Rebin(rebin[type][molt][j]);
      if ((type ==4 || type==5 || type ==8) && j>=5){
      isto_tagli[j]->Rebin(2);
      isto_tagli_true[j]->Rebin(2);
      isto_tagli_false[j]->Rebin(2);
      }
      //if(j<=4) isto_tagli[j]->Rebin(2);
      //else  isto_tagli[j]->Rebin(4);
      // isto[j]->SetTitle("Distribuzione di massa invariante nell'intervallo di p_{T} ["+ int_pt[j] + "] GeV/c e di molteplicita ["+ multiplicity[molt]+") ("+tipo[type]+")");
      // isto_tagli[j]->SetTitle("Invariant mass distribution nell'intervallo di p_{T} ["+ int_pt[j] + "] GeV/c e di molteplicità ["+ multiplicity[molt]+") ("+tipo[type]+")");
      //isto[j]->GetXaxis()->SetTitle(invmass_title[type]);
      //isto[j]->GetYaxis()->SetTitle("Counts");
      isto_tagli[j]->GetXaxis()->SetTitle(invmass_title[type]);
      isto_tagli[j]->GetYaxis()->SetTitle("Counts");
      isto_tagli[j]->GetYaxis()->SetTitleOffset(1.6);
      isto_tagli_true[j]->GetXaxis()->SetTitle(invmass_title[type]);
      isto_tagli_true[j]->GetYaxis()->SetTitle("Counts");
      isto_tagli_true[j]->GetYaxis()->SetTitleOffset(1.6);
      isto_tagli_true[j]->SetLineColor(kRed);
      isto_tagli_false[j]->SetLineColor(kBlue);
  
      char fname1[num_histo]; 
      char fname2[num_histo]; 
      sprintf(fname1,"1f_%d",j);
      sprintf(fname2,"2f_%d",j);
      functionsFirst[j] = new TF1(fname1,"gaus",min_range_signal[type],max_range_signal[type]);
      functionsFirst[j]->SetLineColor(881);   
      functionsFirst[j]->SetParameter(1, massParticle[type]);
      functionsFirst[j]->SetParName(0, "norm");
      functionsFirst[j]->SetParName(1, "mean");
      functionsFirst[j]->SetParName(2, "sigma");
      functionsFirst[j]->SetParLimits(1,min_range_signal[type] ,max_range_signal[type] );
      functionsFirst[j]->SetParLimits(2, 0.001,0.01);
      functionsFirst[j]->SetParLimits(0, 0,1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));
      functionsSecond[j] = new TF1(fname1,"gaus",min_range_signal[type],max_range_signal[type]);
      functionsSecond[j]->SetLineColor(868);   
      functionsSecond[j]->SetParameter(1, massParticle[type]);
      functionsSecond[j]->SetParName(0, "norm");
      functionsSecond[j]->SetParName(1, "mean");
      functionsSecond[j]->SetParName(2, "sigma");
      functionsSecond[j]->SetParLimits(1,min_range_signal[type] ,max_range_signal[type] );
      functionsSecond[j]->SetParLimits(2, 0.001,0.01);
      functionsSecond[j]->SetParLimits(0, 0,1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));
      functions1[j] = new TF1(fname1,"gaus",min_range_signal[type],max_range_signal[type]);
      functions1[j]->SetLineColor(kRed);   
      functions1[j]->SetParName(0, "norm");
      functions1[j]->SetParName(1, "mean");
      functions1[j]->SetParName(2, "sigma");
      functions2[j] = new TF1(fname2,"gaus",min_range_signal[type],max_range_signal[type]);
      functions2[j]->SetLineColor(kMagenta);   
      functions2[j]->SetParName(0, "norm");
      functions2[j]->SetParName(1, "mean");
      functions2[j]->SetParName(2, "sigma");
      //    cout << "qui ok" << endl;

      char fnamebkg2[num_histo]; 
      sprintf(fnamebkg2,"f%d",j); 
      //      bkg2[j] = new TF1(fnamebkg2,"pol2",liminf[type], limsup[type]);
      bkg2[j] = new TF1(fnamebkg2,"pol2",min_histo[type], max_histo[type]);
      bkg2[j]->SetLineColor(kBlue);

      char fnamebkg1[num_histo]; 
      sprintf(fnamebkg1,"f%d",j); 
      //      bkg1[j] = new TF1(fnamebkg1,"pol1",liminf[type], limsup[type]);
      bkg1[j] = new TF1(fnamebkg1,"pol1",min_histo[type], max_histo[type]);
      bkg1[j]->SetLineColor(412);

      char fnamebkgparab[num_histo]; 
      sprintf(fnamebkgparab,"f%d",j); 
      bkgparab[j] = new TF1(fnamebkgparab,fparab,liminf[type], limsup[type],4);
      bkgparab[j]->SetLineColor(kBlue);
      bkgparab[j]->FixParameter(3, type);

      char fnameretta[num_histo]; 
      sprintf(fnameretta,"f%d",j); 
      bkgretta[j] = new TF1(fnameretta,fretta,liminf[type], limsup[type],3);
      bkgretta[j]->SetLineColor(418);
      bkgretta[j]->FixParameter(2, type);
      
      cout << "qui ok" << endl;

      gStyle->SetOptFit(0111);
      cout << " ****ho creato funzioni di fit " << endl;
      cout << "**\n\n******* fit nel range  " << binl[j]<<"-"<<binl[j+1]<<"****************\n"<<endl;

      Bool_t UseTwoGauss=kTRUE;
      Bool_t IsOneGauss=kFALSE;

      //      if (tipo[type]=="XiNeg" || tipo[type]=="XiPos"){
      //	if(year=="2018f1_extra" || year=="2016k"){
	  UseTwoGauss=kTRUE;
	  //}
	  //      }

      canvas[j]= new TCanvas(Form("canvas_v%i",j),Form("canvas_v%i",j), 800,600);
      legend[j]=new TLegend(0.6,0.6,0.9,0.9);

      if(UseTwoGauss){
	char fnametotal[num_histo]; 
	sprintf(fnametotal,"f%d",j); 
	if (isBkgParab)      total[j] = new TF1(fnametotal,"gaus(0)+gaus(3)+pol2(6)",liminf[type], limsup[type]);
	else       total[j] = new TF1(fnametotal,"gaus(0)+gaus(3)+pol1(6)",liminf[type], limsup[type]);
	total[j]->SetLineColor(7); 
	total[j]->SetParName(0, "norm");
	total[j]->SetParName(1, "mean");
	total[j]->SetParName(2, "sigma");
	total[j]->SetParName(3, "norm2");
	total[j]->SetParName(4, "mean2");
	total[j]->SetParName(5, "sigma2");

	//***********limits to obtain a good fit for K0s*******************

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

	//****************************************limits to obtain a good fit for Xi (Cascades)

	if (tipo[type]=="XiNeg" || tipo[type]=="XiPos" || tipo[type]=="Xi"){
	  total[j]->SetParLimits(0, 0.08*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	  //total[j]->SetParLimits(5, 0.001,10);
	  // total[j]->SetParLimits(4, total[j]->GetParameter(1)-1*total[j]->GetParameter(2), total[j]->GetParameter(1)+1*total[j]->GetParameter(2) );
	  total[j]->SetParLimits(1, 1.318, 1.326);
	  total[j]->SetParLimits(2, 0.0012,0.010); //it was 0.001

	  total[j]->SetParLimits(3, 0.08*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  //maximum was wothout 0.3
	  total[j]->SetParLimits(4, 1.318, 1.326);
	  //     total[j]->SetParLimits(5, 0.003, 1);

	  total[j]->SetParLimits(5, 0.0012, 0.007); 
	  if(isMeanFixedPDG){
	  total[j]->FixParameter(1, massParticle[type]);
	  total[j]->FixParameter(4, massParticle[type]);
	  }
	 
	  //    if (molt==1)  total[j]->SetParLimits(5, 0.001, 1);

	}
	if (tipo[type]=="kK0s"){

	  total[j]->SetParLimits(0, 0.08*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	  //total[j]->SetParLimits(5, 0.001,10);
	  // total[j]->SetParLimits(4, total[j]->GetParameter(1)-1*total[j]->GetParameter(2), total[j]->GetParameter(1)+1*total[j]->GetParameter(2) );
	  total[j]->SetParLimits(1, 0.494, 0.501);
	  total[j]->SetParLimits(2, 0.0012,0.010); //it was 0.001

	  total[j]->SetParLimits(3, 0.08*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  //maximum was wothout 0.3
	  total[j]->SetParLimits(4, 0.494, 0.501);
	  //     total[j]->SetParLimits(5, 0.003, 1);

	  total[j]->SetParLimits(5, 0.0012, 0.007);
	  bkgretta[j]->SetParameter(1,0);
	  if(isMeanFixedPDG){
	  total[j]->FixParameter(1, massParticle[type]);
	  total[j]->FixParameter(4, massParticle[type]);
	  }
	  //    if (molt==1)  total[j]->SetParLimits(5, 0.001, 1);

	}

	cout << "\n\n fit gauss1 " << endl;
	isto_tagli[j]-> Fit(functionsFirst[j], "RB"); 
	cout << "\n\n fit gauss2 " << endl;
	// ho visto che peggiora il fit	functionsSecond[j]->SetParLimits(2,functionsSecond[j]->GetParameter(2), 100);
	isto_tagli[j]-> Fit(functionsSecond[j], "RB"); 

	// bkg1[j]->SetRange(liminf[type], limsup[type]);
	// bkg2[j]->SetRange(liminf[type], limsup[type]);
	bkg1[j]->SetRange(min_histo[type], max_histo[type]);
	bkg2[j]->SetRange(min_histo[type], max_histo[type]);
	bkgparab[j]->SetRange(liminf[type], limsup[type]);
	bkgretta[j]->SetRange(liminf[type], limsup[type]);
	total[j]->SetRange(liminf[type], limsup[type]);

	cout << "\n\n fit bkg " << endl;
	if (isBkgParab)    isto_tagli[j]-> Fit(bkgparab[j], "RB0");
	else   isto_tagli[j]-> Fit(bkgretta[j], "RB0");
	functionsFirst[j]->GetParameters(&parTwoGaussParab[j][0]);
	functionsFirst[j]->GetParameters(&parTwoGaussRetta[j][0]);
	functionsSecond[j]->GetParameters(&parTwoGaussParab[j][3]);
	functionsSecond[j]->GetParameters(&parTwoGaussRetta[j][3]);
	if (isBkgParab)  {
	  bkgparab[j]->GetParameters(&parTwoGaussParab[j][6]);
	  total[j]->SetParameters(parTwoGaussParab[j]);
	  if (isMeanFixedPDG){
	  // total[j]->FixParameter(1, massParticle[type]);
	  // total[j]->FixParameter(4, massParticle[type]);
	  }
	}
	else{
	  bkgretta[j]->GetParameters(&parTwoGaussRetta[j][6]);
	  total[j]->SetParameters(parTwoGaussRetta[j]);
	  if (isMeanFixedPDG){
	  // total[j]->FixParameter(1, massParticle[type]);
	  // total[j]->FixParameter(4, massParticle[type]);
	  }
	}


	//	total[j]->SetParLimits(5, functions1[j]->GetParameter(2), 1); 
	// non serve if perché due gaussiane sono uguali
	/*
	if (functions1[j]->GetParameter(2) > functions2[j]->GetParameter(2)){
	  total[j]->SetParLimits(5, functions2[j]->GetParameter(2), 1); 
	}
	else{
	  total[j]->SetParLimits(5, functions1[j]->GetParameter(2), 1); 
	}
	*/

	cout << "\n\n fit total " << endl;
	fFitResultPtr0[j] =  isto_tagli[j]->Fit(total[j],"SRB");//per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
	//la gaussiana più larga deve esserte quella più bassa
	if (total[j]->GetParameter(2) > total[j]->GetParameter(5)){
	  if (total[j]->GetParameter(0) > total[j]->GetParameter(3)) IsOneGauss=kTRUE;
	}
	else {
	  if (total[j]->GetParameter(0) < total[j]->GetParameter(3)) IsOneGauss=kTRUE;
	}
       

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
	if (isBkgParab){
	  bkg2[j]->FixParameter(0,total[j]->GetParameter(6));
	  bkg2[j]->FixParameter(1,total[j]->GetParameter(7));
	  bkg2[j]->FixParameter(2,total[j]->GetParameter(8));
	  bkgparab[j]->FixParameter(0,total[j]->GetParameter(6));
	  bkgparab[j]->FixParameter(1,total[j]->GetParameter(7));
	  bkgparab[j]->FixParameter(2,total[j]->GetParameter(8));
	}
	else{
	  bkg1[j]->FixParameter(0,total[j]->GetParameter(6));
	  bkg1[j]->FixParameter(1,total[j]->GetParameter(7));
	  bkgretta[j]->FixParameter(0,total[j]->GetParameter(6));
	  bkgretta[j]->FixParameter(1,total[j]->GetParameter(7));
	}

	functions1[j]->SetLineColor(kRed);   
	isto_tagli[j]-> Fit(functions1[j], "RB+");
	isto_tagli[j]-> Fit(functions2[j], "RB+");
	if (isBkgParab)    isto_tagli[j]-> Fit(bkg2[j], "BR+");
	else   isto_tagli[j]-> Fit(bkg1[j], "BR+");
	// isto_tagli[j]-> Fit(functionsFirst[j], "RB+"); 
	// isto_tagli[j]-> Fit(functionsSecond[j], "RB+");  
	isto_tagli[j]->GetYaxis()->SetRangeUser(0, 1.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));
	isto_tagli[j]->GetXaxis()->SetRangeUser(min_histo[type],max_histo[type]);
	isto_tagli_true[j]->GetXaxis()->SetRangeUser(min_histo[type],max_histo[type]);
	isto_tagli_false[j]->GetXaxis()->SetRangeUser(min_histo[type],max_histo[type]);
	//    cout << "qua ok" << endl;

	for (Int_t i=1; i<=isto_tagli[j]->GetNbinsX(); i++){
	  isto_tagli_false[j]->SetBinContent(i, 	isto_tagli[j]->GetBinContent(i)-isto_tagli_true[j]->GetBinContent(i));
	}

	canvas[j]->cd();
	gPad->SetLogy();
	isto_tagli[j]->Draw("");
	legend[j]->AddEntry(total[j],"Total", "l");
	legend[j]->AddEntry(functions1[j],"Gaussian", "l");
	legend[j]->AddEntry(functions2[j],"Gaussian", "l");
	if (isBkgParab)        legend[j]->AddEntry(bkg2[j],"2nd degree polynomial","l");
	else        legend[j]->AddEntry(bkg1[j],"1st degree polynomial","l");
	legend[j]->Draw("same");
	//    cout << "ho disegnato legenda" << endl;
	f->WriteTObject(canvas[j]);
	canvas[j]->Close();
	//    cout << " ho salvato canvas " << endl;
	isto_tagli[j]-> Write();
	isto_tagli_true[j]-> Write();
	isto_tagli_false[j]-> Write();

	canv[molt]->cd(j+1);
	//	gPad->SetLogy();
		isto_tagli[j]->GetYaxis()->SetRangeUser(0,0.02*isto_tagli[j]->GetMaximum());
	isto_tagli[j]->DrawCopy();
	cout << "\n\n\n\n ho salvato canvas con istogramma per m "<< molt <<" e pt " << j << " \n\n\n" << endl;
	TMatrixDSym cov =   fFitResultPtr0[j]->GetCovarianceMatrix();
	Double_t cov_mean = cov[1][4];
	Double_t cov_sigma = cov[2][5];
	// cout << pow(total[j]->GetParError(1),2) << " " << pow(total[j]->GetParError(4),2) << " " << cov_mean << endl;
	// cout << pow(total[j]->GetParError(2),2) << " " << pow(total[j]->GetParError(5),2) << " " << cov_mean << endl;
	// cout << cov_mean << " " << cov_sigma << endl;
	total[j]->FixParameter(6,0);
	total[j]->FixParameter(7,0);
	if (isBkgParab)    total[j]->FixParameter(8,0);
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

	total[j]->SetLineColor(881);
	total[j]->Draw("same");

	lineCentralSX[j]= new TLine(mean[j]-sigmacentral*sigma[j], 0, mean[j]-sigmacentral*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineCentralDX[j]= new TLine(mean[j]+sigmacentral*sigma[j], 0, mean[j]+sigmacentral*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineCentralSX[j]->Draw("same");
	lineCentralDX[j]->Draw("same");

	lineSidebandsSX[j]= new TLine(mean[j]-nsigmamin*sigma[j], 0, mean[j]-nsigmamin*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineSidebandsDX[j]= new TLine(mean[j]+nsigmamin*sigma[j], 0, mean[j]+nsigmamin*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineSidebandsSX[j]->Draw("same");
	lineSidebandsDX[j]->Draw("same");
  
      }


      if (IsOneGauss){
      //      else{
	char fnametotal[num_histo]; 
	sprintf(fnametotal,"f%d",j); 
	if (isBkgParab)      total[j] = new TF1(fnametotal,"gaus(0)+pol2(3)",liminf[type], limsup[type]);
	else       total[j] = new TF1(fnametotal,"gaus(0)+pol1(3)",liminf[type], limsup[type]);

	total[j]->SetLineColor(7);
	total[j]->SetParName(0, "norm");
	total[j]->SetParName(1, "mean");
	total[j]->SetParName(2, "sigma");
	total[j]->SetParLimits(2, 0.001,10);
	total[j]->SetParLimits(0,0,1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  

	if (tipo[type]=="XiNeg" || tipo[type]=="XiPos" || tipo[type]=="Xi") {
	  total[j]->SetParLimits(1, 1.318, 1.326);
	  total[j]->SetParLimits(0, 0.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	  total[j]->SetParLimits(2, 0.001,0.004);
	  if (isMeanFixedPDG){
	  total[j]->FixParameter(1, massParticle[type]);
	  }
	}

	if (tipo[type]=="kK0s"){
	  total[j]->SetParLimits(1, 0.494, 0.501);
	  total[j]->SetParLimits(0, 0.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()),1.1*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	  total[j]->SetParLimits(2, 0.001,0.004);
	  if (isMeanFixedPDG){
	  total[j]->FixParameter(1, massParticle[type]);
	  }
	}

	cout << "\n\n fit gauss1 " << endl;
	functions1[j]->SetLineColor(kRed);   
	isto_tagli[j]-> Fit(functionsFirst[j], "R0"); 
	// bkg2[j]->SetRange(liminf[type], limsup[type]);
	// bkg1[j]->SetRange(liminf[type], limsup[type]);
	bkg1[j]->SetRange(min_histo[type], max_histo[type]);
	bkg2[j]->SetRange(min_histo[type], max_histo[type]);
	bkgparab[j]->SetRange(liminf[type], limsup[type]);
	bkgretta[j]->SetRange(liminf[type], limsup[type]);
	total[j]->SetRange(liminf[type], limsup[type]);
	cout << "\n\n fit bkg " << endl;
	if (isBkgParab)    isto_tagli[j]-> Fit(bkgparab[j], "RB0");
	else     isto_tagli[j]-> Fit(bkgretta[j], "RB0");
	functionsFirst[j]->GetParameters(&parOneGaussParab[j][0]);
	functionsFirst[j]->GetParameters(&parOneGaussRetta[j][0]);
	if (isBkgParab)	{
	  bkgparab[j]->GetParameters(&parOneGaussParab[j][3]);
	  total[j]->SetParameters(parOneGaussParab[j]);
	  if (isMeanFixedPDG){
	    //	  total[j]->FixParameter(1, massParticle[type]);
	  }
	}
	else{
	bkgretta[j]->GetParameters(&parOneGaussRetta[j][3]);
	total[j]->SetParameters(parOneGaussRetta[j]);
	  if (isMeanFixedPDG){
	    //	  total[j]->FixParameter(1, massParticle[type]);
	  }
	}

	cout << "\n\n fit total " << endl;
	fFitResultPtr0[j] =  isto_tagli[j]->Fit(total[j],"SR");//per errore gaussiana, S indica che il risultato del fit e' accessibile da fFitResultPtr0
	//isto_tagli[j]-> Fit(total[j], "R");
	//qui faccio cose non ben chiare ma mi pare che funzioni
	totalbis[j] = (TF1*)total[j]->Clone();
	fFitResultPtr1[j] = fFitResultPtr0[j] ;

	functions1[j]->FixParameter(0,total[j]->GetParameter(0));
	functions1[j]->FixParameter(1,total[j]->GetParameter(1));
	functions1[j]->FixParameter(2,total[j]->GetParameter(2));
	if (isBkgParab){
	  bkg2[j]->FixParameter(0,total[j]->GetParameter(3));
	  bkg2[j]->FixParameter(1,total[j]->GetParameter(4));
	  bkg2[j]->FixParameter(2,total[j]->GetParameter(5));
	  bkgparab[j]->FixParameter(0,total[j]->GetParameter(3));
	  bkgparab[j]->FixParameter(1,total[j]->GetParameter(4));
	  bkgparab[j]->FixParameter(2,total[j]->GetParameter(5));
	}
	else{
	  bkg1[j]->FixParameter(0,total[j]->GetParameter(3));
	  bkg1[j]->FixParameter(1,total[j]->GetParameter(4));
	  bkgretta[j]->FixParameter(0,total[j]->GetParameter(3));
	  bkgretta[j]->FixParameter(1,total[j]->GetParameter(4));
	}

	if (isBkgParab)  isto_tagli[j]-> Fit(bkg2[j], "RB+"); 
	else  isto_tagli[j]-> Fit(bkg1[j], "RB+"); 
	isto_tagli[j]-> Fit(functions1[j], "RB+");
	//isto[j]-> Write();
	isto_tagli[j]->GetXaxis()->SetRangeUser(min_histo[type],max_histo[type]);
	isto_tagli[j]->GetYaxis()->SetRangeUser(0,0.004*isto_tagli[j]->GetMaximum());
	canvas[j]->cd();
	isto_tagli[j]->Draw("");
	legend[j]->AddEntry(total[j],"Total", "l");
	legend[j]->AddEntry(functions1[j],"Gaussian", "l");
	if (isBkgParab)  legend[j]->AddEntry(bkg2[j],"2nd degree polynomial","l");
	else  legend[j]->AddEntry(bkg1[j],"1st degree polynomial","l");
	legend[j]->Draw("same");
	//	cout << "ho disegnato legenda" << endl;
	f->WriteTObject(canvas[j]);
	canvas[j]->Close();
	isto_tagli[j]-> Write();

	canv[molt]->cd(j+1);
	isto_tagli[j]->DrawCopy();

	total[j]->FixParameter(3,0);
	total[j]->FixParameter(4,0);
	if (isBkgParab)    total[j]->FixParameter(5,0);
	totalbis[j]->FixParameter(0,0);
	totalbis[j]->FixParameter(1,0);
	totalbis[j]->FixParameter(2,0);
	mean[j]=functions1[j]->GetParameter(1);
	errmean[j]=total[j]->GetParError(1);
	sigma[j]=functions1[j]->GetParameter(2);
	errsigma[j]=total[j]->GetParError(2);

	total[j]->SetLineColor(881);
	total[j]->Draw("same");

	lineCentralSX[j]= new TLine(mean[j]-sigmacentral*sigma[j], 0, mean[j]-sigmacentral*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineCentralDX[j]= new TLine(mean[j]+sigmacentral*sigma[j], 0, mean[j]+sigmacentral*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineCentralSX[j]->Draw("same");
	lineCentralDX[j]->Draw("same");

	lineSidebandsSX[j]= new TLine(mean[j]-nsigmamin*sigma[j], 0, mean[j]-nsigmamin*sigma[j],0.2*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineSidebandsDX[j]= new TLine(mean[j]+nsigmamin*sigma[j], 0, mean[j]+nsigmamin*sigma[j],0.22*isto_tagli[j]->GetBinContent(isto_tagli[j]->GetMaximumBin()));  
	lineSidebandsSX[j]->Draw("same");
	lineSidebandsDX[j]->Draw("same");

      }

      //******Rimepio isto con sigma, mean, purezza,...**********************************************************************

      entries_range_false[j]=0;
      entries_range_true[j]=0;
      entries_range[j]=0;
      entries_sideband[j]=0;
      entries_sideband_false[j]=0;      
      entries_sideband_true[j]=0;

      for(Int_t l=isto_tagli[j]->GetXaxis()->FindBin(mean[j]-sigmacentral*sigma[j]); l<=isto_tagli[j]->GetXaxis()->FindBin(mean[j]+sigmacentral*sigma[j]); l++){
	entries_range[j]+=isto_tagli[j]->GetBinContent(l);
	entries_range_true[j]+=isto_tagli_true[j]->GetBinContent(l);
	entries_range_false[j]+=isto_tagli_false[j]->GetBinContent(l);
      }

      //      entries_range_false[j]=entries_range[j]-entries_range_true[j];
      err_range_false[j]=sqrt(entries_range[j]+entries_range_true[j]);
      for(Int_t l=isto_tagli[j]->GetXaxis()->FindBin(min_histo[type]); l<=isto_tagli[j]->GetXaxis()->FindBin(mean[j]-nsigmamin*sigma[j]); l++){
	entries_sideband[j]+=isto_tagli[j]->GetBinContent(l);
	entries_sideband_true[j]+=isto_tagli_true[j]->GetBinContent(l);
      }
      for(Int_t l=isto_tagli[j]->GetXaxis()->FindBin(mean[j]+nsigmamin*sigma[j]); l<=isto_tagli[j]->GetXaxis()->FindBin(max_histo[type]); l++){
	entries_sideband[j]+=isto_tagli[j]->GetBinContent(l);
	entries_sideband_true[j]+=isto_tagli_true[j]->GetBinContent(l);
      }
      entries_sideband_false[j]=entries_sideband[j]-entries_sideband_true[j];
      err_sideband_false[j]=sqrt(entries_sideband[j]+entries_sideband_true[j]);
      ////////////////////////////////////
      //never used    s1[j]=functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2));
      //never used    s2[j]=functions2[j]->Integral(functions2[j]->GetParameter(1)-sigmacentral*functions2[j]->GetParameter(2),functions2[j]->GetParameter(1)+sigmacentral*functions2[j]->GetParameter(2));
      st[j]=total[j]->Integral(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j]); //la funzione total ha parametri del bkg settati a zero!
      IntegralSignalAllRange[j]=total[j]->Integral(min_histo[type], max_histo[type]); //la funzione total ha parametri del bkg settati a zero!
      // cout <<"con totale " <<  st[j]<< endl;
      //   st[j]=s1[j]+s2[j];
      // cout <<"con somma " <<  st[j]<< endl;

      tot[j]=entries_range[j]*isto_tagli[j]->GetBinWidth(1);

      if (isBkgParab)  {
	b[j]=bkg2[j]->Integral(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j]);
	bside[j]=bkg2[j]->Integral(min_histo[type],mean[j]-nsigmamin*sigma[j]) + bkg2[j]->Integral(mean[j]+nsigmamin*sigma[j], max_histo[type]);
      }
      else{
	b[j]=bkg1[j]->Integral(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j]);
	bside[j]=bkg1[j]->Integral(min_histo[type],mean[j]-nsigmamin*sigma[j]) + bkg1[j]->Integral(mean[j]+nsigmamin*sigma[j], max_histo[type]);
      }

      //    errbside[j]=bkg2[j]->IntegralError(mean[j]-nsigmamax*sigma[j], mean[j]-nsigmamin*sigma[j]) + bkg2[j]->IntegralError(mean[j]+nsigmamin*sigma[j], mean[j]+nsigmamax*sigma[j]);
      bin_contentS1[j]=st[j]; //errS1 e' l'errore associato
      bin_contentS2[j]= tot[j]-b[j];//errS2
     
      bin_contentSSB1[j]= st[j]/(st[j]+b[j]);//errSSB1
      bin_contentSSB2[j]= (tot[j]-b[j])/tot[j];//errSSB2

      if (b[j]>0)  {
	bin_contentSB1[j]=st[j]/b[j];//errSB1
	bin_contentSB2[j]= (tot[j]-b[j])/b[j];//errSB2
      }
      else{
	bin_contentSB1[j]=0;
	bin_contentSB2[j]=0;
      }


      sigmas1[j]=total[j]->IntegralError(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j],fFitResultPtr0[j] ->GetParams(),(fFitResultPtr0[j]->GetCovarianceMatrix()).GetMatrixArray());//ok
      sigmab1[j]=totalbis[j]->IntegralError(mean[j]-sigmacentral*sigma[j],mean[j]+sigmacentral*sigma[j],fFitResultPtr1[j] ->GetParams(),(fFitResultPtr1[j]->GetCovarianceMatrix()).GetMatrixArray()); //ok
      errb[j]=sigmab1[j]; //ok
      errbside[j]=totalbis[j]->IntegralError(min_histo[type], mean[j]-nsigmamin*sigma[j],fFitResultPtr1[j] ->GetParams(),(fFitResultPtr1[j]->GetCovarianceMatrix()).GetMatrixArray())+ totalbis[j]->IntegralError(mean[j]+nsigmamin*sigma[j], limsup[type],fFitResultPtr1[j] ->GetParams(),(fFitResultPtr1[j]->GetCovarianceMatrix()).GetMatrixArray()); //ok
      errS1[j]=sigmas1[j]; //ok
      errS2[j]=sqrt( pow(isto_tagli[j]->GetBinWidth(1),2)*entries_range[j] + pow(sigmab1[j],2));    //ok
      errSB1[j]=   sqrt(pow(sigmas1[j],2)/pow(b[j],2) + pow(sigmab1[j],2)*pow(st[j],2)/pow(b[j],4));//ok
      errSB2[j]=sqrt( (pow(isto_tagli[j]->GetBinWidth(1),2)*entries_range[j] )/pow(b[j],2) + pow(sigmab1[j],2)*pow(tot[j],2)/pow(b[j],4)) ; //ok
      errSSB1[j]=sqrt((pow(st[j],2)*pow(sigmab1[j],2)+pow(b[j],2)*pow(sigmas1[j],2))/pow(st[j]+b[j],4)); //ok
      //credo errato      errSSB2[j]=pow(b[j],2)/pow(tot[j],3)-pow(sigmab1[j],2)/pow(tot[j],2);
      errSSB2[j]=sqrt(pow(b[j],2)/pow(tot[j],4)*pow(isto_tagli[j]->GetBinWidth(1),2)*entries_range[j]+pow(sigmab1[j],2)/pow(tot[j],2));//ok

      Bool_t ShowInHisto=kTRUE;
      if(ShowInHisto){
	if(mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
	  histo_sigma->SetBinContent(j+1,sigma[j]);
	  histo_sigma->SetBinError(j+1,errsigma[j]);
	  histo_mean->SetBinContent(j+1,mean[j]);
	  histo_mean->SetBinError(j+1,errmean[j]);

	} else {
	  histo_sigma->SetBinContent(j+1,0);
	  histo_sigma->SetBinError(j+1,0);
	 histo_mean->SetBinContent(j+1,0);
	 histo_mean->SetBinError(j+1,0);
 
	}
	if(total[j]->GetNDF() != 0){
	  histo_chis->SetBinContent(j+1,total[j]->GetChisquare()/total[j]->GetNDF());
	  histo_chis->SetBinError(j+1,0);
	} else {
	  histo_chis->SetBinContent(j+1,0);
	  histo_chis->SetBinError(j+1,0);
	}
	// histo_stat->SetBinContent(j+1,isto[j]->GetEntries());
	
	//	if((total[j]->GetChisquare()/total[j]->GetNDF())<=3. && mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && bin_contentSB1[j]>0 && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
	if( isSignalFromIntegral && mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && bin_contentSB1[j]>0 && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
	  histo_SB->SetBinContent(j+1,bin_contentSB1[j]);
	  histo_SB->SetBinError(j+1, errSB1[j]);
	  histo_SSB->SetBinContent(j+1,bin_contentSSB1[j]);
	  histo_SSB->SetBinError(j+1, errSSB1[j]);
	  histo_S->SetBinContent(j+1,bin_contentS1[j]/isto_tagli[j]->GetBinWidth(1)/NEvents[molt]/(binl[j+1]-binl[j]));
	  histo_S->SetBinError(j+1, errS1[j]/isto_tagli[j]->GetBinWidth(1)/NEvents[molt]/(binl[j+1]-binl[j]));
	  histo_Bcentral->SetBinContent(j+1,b[j]);
	  histo_Bcentral->SetBinError(j+1, errb[j]);
	  histo_Bside->SetBinContent(j+1,bside[j]);
	  histo_Bside->SetBinError(j+1, errbside[j]);
	  histo_BcentralFalse->SetBinContent(j+1,entries_range_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_BcentralFalse->SetBinError(j+1, err_range_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_BsideFalse->SetBinContent(j+1,entries_sideband_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_BsideFalse->SetBinError(j+1, err_sideband_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_FracStrangePeak->SetBinContent(j+1,(entries_range_false[j]*isto_tagli[j]->GetBinWidth(1)-b[j])/entries_range[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_FracStrangePeak->SetBinError(j+1,0);
	  histo_BRatioCentral->SetBinContent(j+1,b[j]/entries_range_false[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_SRatioCentral->SetBinContent(j+1,st[j]/entries_range_true[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_BRatioSide->SetBinContent(j+1,bside[j]/entries_sideband_false[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_BDoubleRatio->SetBinContent(j+1,(b[j]/bside[j])/((Float_t)entries_range_false[j]/entries_sideband_false[j]));

	  histo_BRatioCentral->SetBinError(j+1,0);
	  histo_BRatioSide->SetBinError(j+1,0);
	  histo_BDoubleRatio->SetBinError(j+1,0);
	  histo_SRatioCentral->SetBinError(j+1, sqrt(pow(sigmas1[j]/st[j],2) + 1./entries_range_true[j])*st[j]/entries_range_true[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_signal_int->SetBinContent(j+1,st[j]/isto_tagli[j]->GetBinWidth(1)/events);
	  histo_signal_int->SetBinError(j+1,sigmas1[j]/isto_tagli[j]->GetBinWidth(1)/events);
	  histo_signal_int_pure->SetBinContent(j+1,st[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_signal_int_pure->SetBinError(j+1,sigmas1[j]/isto_tagli[j]->GetBinWidth(1));

	  cout << "caso A"<< endl;

	  //	} else if((total[j]->GetChisquare()/total[j]->GetNDF())>3. && mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && bin_contentSB2[j]>0 && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
	} 
	else if(!isSignalFromIntegral && mean[j]>lim_inf_mean[type] && mean[j]<lim_sup_mean[type] && bin_contentSB2[j]>0 && sigma[j]<lim_sup_sigma[type] && errsigma[j]<lim_sup_errsigma[type] && errmean[j]<lim_sup_errmean[type]){
	  histo_SB->SetBinContent(j+1,bin_contentSB2[j]);
	  histo_SB->SetBinError(j+1,errSB2[j]);
	  histo_SSB->SetBinContent(j+1,bin_contentSSB2[j]);
	  histo_SSB->SetBinError(j+1, errSSB2[j]);
	  histo_S->SetBinContent(j+1,bin_contentS2[j]/isto_tagli[j]->GetBinWidth(1)/NEvents[molt]/(binl[j+1]-binl[j]));
	  histo_S->SetBinError(j+1,errS2[j]/isto_tagli[j]->GetBinWidth(1)/NEvents[molt]/(binl[j+1]-binl[j]));
	  histo_Bcentral->SetBinContent(j+1,b[j]);
	  histo_Bcentral->SetBinError(j+1, errb[j]);
	  histo_Bside->SetBinContent(j+1,bside[j]);
	  histo_Bside->SetBinError(j+1, errbside[j]);
	  histo_BcentralFalse->SetBinContent(j+1,entries_range_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_BcentralFalse->SetBinError(j+1, err_range_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_BsideFalse->SetBinContent(j+1,entries_sideband_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_BsideFalse->SetBinError(j+1, err_sideband_false[j]*isto_tagli[j]->GetBinWidth(1));
	  histo_FracStrangePeak->SetBinContent(j+1,(entries_range_false[j]*isto_tagli[j]->GetBinWidth(1)-b[j])/entries_range[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_BRatioCentral->SetBinContent(j+1,b[j]/entries_range_false[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_SRatioCentral->SetBinContent(j+1,st[j]/entries_range_true[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_BRatioSide->SetBinContent(j+1,bside[j]/entries_sideband_false[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_BDoubleRatio->SetBinContent(j+1,(b[j]/bside[j])/((Float_t)entries_range_false[j]/entries_sideband_false[j]));

	  histo_BRatioCentral->SetBinError(j+1,0);
	  histo_BRatioSide->SetBinError(j+1,0);
	  histo_BDoubleRatio->SetBinError(j+1,0);
	  histo_SRatioCentral->SetBinError(j+1, sqrt(pow(sigmas1[j]/st[j],2) + 1./entries_range_true[j])*st[j]/entries_range_true[j]/isto_tagli[j]->GetBinWidth(1));
	  histo_FracStrangePeak->SetBinError(j+1,0);
	  histo_signal_int->SetBinContent(j+1,(tot[j]-b[j])/isto_tagli[j]->GetBinWidth(1)/events);
	  histo_signal_int->SetBinError(j+1,errS2[j]/isto_tagli[j]->GetBinWidth(1)/events);
	  histo_signal_int_pure->SetBinContent(j+1,(tot[j]-b[j])/isto_tagli[j]->GetBinWidth(1));
	  histo_signal_int_pure->SetBinError(j+1,errS2[j]/isto_tagli[j]->GetBinWidth(1));

	  cout << "caso B"<< endl;
	  cout <<"mean1 of total " << total[j]->GetParameter(1) << endl; 
	  cout <<"mean2 of total " << total[j]->GetParameter(4) << endl; 
	  cout <<"mean of function1 " <<  functions1[j]->GetParameter(1)<< endl;	  cout << "mean put in histo " << mean[j]<< endl;
	  /*
cout << functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))<< endl; 
	  cout << bin_contentSB2[j] << endl;
	  cout << entries_range[j]*isto_tagli[j]->GetBinWidth(1)/(-functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))+entries_range[j]*isto_tagli[j]->GetBinWidth(1)) << endl;
	  */

	}  else {
	  cout << "caso C " << endl;
	  cout << "S su B varra' zero perche  condizioni non sono soddisfatte " << endl; 
	  cout << "reduced chi " << total[j]->GetChisquare()/total[j]->GetNDF()<< "\nmean " << mean[j] << " >? " << lim_inf_mean[type] << " <? " <<lim_sup_mean[type] <<"\ntotale in central region from bin counting " << tot[j] << "\nB in central region from integral " << b[j] << " \n S/B " <<bin_contentSB2[j]<< "\n sigma " << sigma[j] << " <? " << lim_sup_sigma[type] << "\nerr sigma " << errsigma[j] << " <? " << lim_sup_errsigma[type] << "\nerr mean " << errmean[j] << " <? " << lim_sup_errmean[type] << endl;
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
	  histo_BcentralFalse->SetBinContent(j+1,0);
	  histo_BcentralFalse->SetBinError(j+1, 0);
	  histo_BsideFalse->SetBinContent(j+1,0);
	  histo_BsideFalse->SetBinError(j+1, 0);
	  histo_FracStrangePeak->SetBinContent(j+1,0);
	  histo_BRatioCentral->SetBinContent(j+1,0);
	  histo_BRatioSide->SetBinContent(j+1,0);
	  histo_BDoubleRatio->SetBinContent(j+1,0);
	  histo_SRatioCentral->SetBinContent(j+1,0);
	  histo_BRatioCentral->SetBinError(j+1,0);
	  histo_BRatioSide->SetBinError(j+1,0);
	  histo_BDoubleRatio->SetBinError(j+1,0);
	  histo_SRatioCentral->SetBinError(j+1,0);
	  histo_FracStrangePeak->SetBinError(j+1,0);
	  histo_signal_int->SetBinContent(j+1,0);
	  histo_signal_int->SetBinError(j+1,0);
	  histo_signal_int_pure->SetBinContent(j+1,0);
	  histo_signal_int_pure->SetBinError(j+1,0);

	}

	cout << "\n more results: " << endl;
	cout << "mu - nsigma " << mean[j]-sigmacentral*sigma[j] << " mu + nsigmna " << mean[j]+sigmacentral*sigma[j]<< endl;
	cout << "min_histo "<< min_histo[type] << " mu - nsigmnamin *sigma" << mean[j]-nsigmamin*sigma[j]<< " mu + nsigmamin*sigma " <<  mean[j]+nsigmamin*sigma[j] <<  " limsup " << max_histo[type]<< endl;

	cout << "bside from integral " << bside[j]<< ", from bin counting "<< entries_sideband_false[j]*isto_tagli[j]->GetBinWidth(1)<< endl;
	cout << "b in central region from integral " << b[j]<< ", from bin counting "<< entries_range_false[j]*isto_tagli[j]->GetBinWidth(1)<< endl;

	cout << "****************"<<endl;
	// cout << total[j]->GetChisquare()<<endl;
	// cout << total[j]->GetNDF()<<endl;
	// cout << total[j]->GetChisquare()/total[j]->GetNDF() << endl;
	// cout <<functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2)) << endl;
	// cout << bkg2[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))<< endl;
    
	/*
	  cout<< "-----------------" << endl;
	  cout << entries_range[j]*isto_tagli[j]->GetBinWidth(1) << endl;
	  cout << (-functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))+entries_range[j]*isto_tagli[j]->GetBinWidth(1)) << endl;
	  cout << bin_contentSB2[j] << endl;
	  cout << entries_range[j]*isto_tagli[j]->GetBinWidth(1)/(-functions1[j]->Integral(functions1[j]->GetParameter(1)-sigmacentral*functions1[j]->GetParameter(2),functions1[j]->GetParameter(1)+sigmacentral*functions1[j]->GetParameter(2))+entries_range[j]*isto_tagli[j]->GetBinWidth(1)) << endl;
	*/
	cout<<"---------------"<<j<<endl;
      }
    }

    for(Int_t j=1;j<num_histo; j++){
      cout << "segnale central "  << st[j]<<" segnale all range " <<IntegralSignalAllRange[j]<< endl;
	cout << "Frazione dell'integrale della funzione di segnale compresa tra (mu - sigma*nsigmacentral, mu + sigma*nsigmacentral) " << st[j]/IntegralSignalAllRange[j]<< endl; 
     if(bin_contentSB1[j]<0 || bin_contentSB2[j]<0){
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
    histo_BcentralFalse->Write();
    histo_BsideFalse->Write();
    histo_FracStrangePeak->Write();
    histo_BRatioCentral->Write();
    histo_SRatioCentral->Write();
    histo_BRatioSide->Write();
    histo_BDoubleRatio->Write();
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

    for (Int_t molt=0; molt < mult+1; molt++){
    cout <<"n trigger particles for m " << molt << " " <<  NEvents[molt] << endl;
    }
}
