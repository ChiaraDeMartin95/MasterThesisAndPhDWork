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
#include <TLegend.h>
#include <THStack.h>
#include <TSpline.h>

void confronto_histo_thesis(Bool_t isSignalFromIntegral=0, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Float_t PtTrigMin=3, Int_t type=0, Int_t sysTrigger=0, Int_t sysV0=0, Int_t sys=0, TString yearMC="17pq_hXi"/*"17pq_hK0s"/*"2019h11c_extra_HM_hK0s"/"1617MC_hK0s"/*"2015g3b1"/"AllMC_hXi"/*"2018f1_extra_hK0s"*/, TString yeardati="17pq_hXi"/*"17pq_hK0s"/*"2016k_HM_hK0s"/"1617_hK0s"/*"2015f"/"2016kehjl_hK0s"/"Run2DataRed_MECorr_hXi"*/, TString year="17pq_hXi"/*"17pq_hK0s"/*"2016k_HM_hK0s"/"1617_hK0s"/*"2016kehjl_18f1_extra_hK0s"/*"2015f"*"Run2DataRed_MECorr_AllMC_hXi"*/, TString year0="2016", Bool_t isParticleStudy=1, Bool_t isSystV0Analysis=0, Int_t rap=0, TString Path1="", TString Path1MC = "",Bool_t SkipAssoc=1, Int_t PtBinning = 0, Bool_t isHM=0, Bool_t isINEL=1){ 

  if (isINEL){
    SkipAssoc =0;
    PtTrigMin=0.;
    rap =1;
    isHM =1;
    if (type==0){
      year = "161718_HM_hK0s_INELgt0";
      yeardati = "161718_HM_hK0s_INELgt0";
      yearMC ="";
    }
    PtBinning = 0;
  }
  TString  SisINEL[2] = {"", "_INEL"};

  //isSystV0Analysis=1 if I want to compare different selections used to identify the particle for a given multiplcity (to be chosen!); =0 if I want to compare diffrerent multiplcity ranges (but same selections)
  Int_t chosenmolt=5; //choose the multiplicity you desire
  const Int_t numsysV0=21; 
  const Int_t mult=5; //numero intervalli molteplicita'
  const Int_t num_tipo=12; //numero tipologie di particelle
  const Int_t num_isto=12; //numero di istogrammi di tipologia diversa per dati veri
  const Int_t num_isto_MC=13; //numero di istogrammi di tipologia diversa per MC
  const Int_t num_isto_cfr=2; //numero di istogrammi di confronto dati/MC
  const Int_t num_isto_lambda=1; //numero di istogrammi di confronto lambda/antilambda
  const Int_t num_tagli=3;

  TString Srap[2] = {"_Eta0.8", "_y0.5"};
  TString    SSkipAssoc[2]={"_AllAssoc", ""};

  TString isDataorMC[2]={"Data", "MC"};
  TString StringSystV0Analysis[2]={"", "_SystV0Analysis"};
  if (isParticleStudy==0 && numsysV0!=1) return;
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};
  Int_t num_isto_finale=num_isto;
  if (isParticleStudy) num_isto_finale=num_isto+1;
  TString tipo2[num_tipo]={"K^{0}_{s}", "#Lambda","#bar{#Lambda}", "#Lambda + #bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}", "#Xi", "#Omega"};
  TString tipoNotSign[num_tipo]={"K^{0}_{s}", "#Lambda","#bar{#Lambda}", "#Lambda + #bar{#Lambda}", "#Xi", "#Xi", "#Omega", "#Omega", "#Xi", "#Omega"};
  TString tipo[num_tipo]={"K0s", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega"};
  TString tipobis[num_tipo]={/*"kK0s"*/"K0s", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega"};
  TString SParticleType[num_tipo]={"", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos", "Xi", "Omega"};
  TString tipo_histo[num_isto+1]={"sigma", "mean", "chis", "SB", "S","SEffCorr", "SSB","SRatioCentral", "BRatioCentral", "BRatioSide", "BDoubleRatio", "Efficiency", "EfficiencyRelError"};
  TString nomi[num_isto+1]={"Invariant mass resolution", 
			    tipo2[type]+" mass",
			    "#chi^{2}/dof", 
			    "S/B within 3#sigma", 
			    "S within 3#sigma",
			    "dN/dp_{T} 1/N_{ev}",
			    "S/(S+B) within 3#sigma", 
			    //			    "Estimate of gaussian bkg/all selected within 3#sigma",
			    "S (bin counting - integral bkg) / S (MCTruth)",
			    "B within 3#sigma (bkg integral/MC truth)", 
			    "B in sidebands  (bkg integral/MCtruth)",
			    "B double ratio (sidebands/central)_{int} /  (sidebands/central)_{MCTruth}",
			    "Efficiency",
			    "Relative errors of efficiency"};
  TString tipo_histo_MC[num_isto_MC]={"sigma", "mean", "chis", "SB","SSB", "S", "signal_int","eff","eff_bis", "eff_tagli", "eff_tagli_bis","eff_tris","bontafit"};
  TString nomi_MC[num_isto_MC]={
    "Risoluzione della massa invariante",
    "Massa del " +tipo2[type],
    "#chi^{2}/dof", 
    "S/B entro 3#sigma", 
    "S/(S+B) entro 3#sigma", 
    "S entro 3#sigma", 
    "Numero di " + tipo2[type] +" per evento",  //ossia "Integrale segnale dopo applicazione tagli / numero eventi per classe di molteplicita' (ossia numero di V0 per evento in data classe di molteplicita')"
    "Efficienza x acc. x branching ratio (calcolo errore tipo Grazia)", 
    "Efficienza x accettanza x branching ratio, errore Bayes", 
    "Efficienza delle selezioni applicate (calcolo errore tipo Grazia)", //Efficienza dei tagli applicati (rapporto sig MC dopo e prima tagli)"
    "Frazione di " +tipo2[type]+  " identificati con le selezioni topologiche",
    "Efficienza dei tagli, divisione tra istogrammi",
    "Rapporto tra integrale segnale MC e segnale MC dopo selezioni"};
  TString nomi_completi[num_isto+1];
  TString nomi_MC_completi[num_isto_MC];
  TString nomi_cfr[num_isto_cfr] = {"Numero "+ tipo2[type]+ " dati veri/numero "+ tipo2[type]+ " MC (dopo selezioni)", "Numero di "+ tipo2[type]+ " per evento corretto per efficienza"}; //il primo integrale e':"Integrale dati veri per evento/integrale segnale MC per evento, entrambi dopo tagli";  il secondo istogramma e': Integrale segnale dopo applicazione tagli (ossia numero di V0 per evento) corretto per efficienza"
  TString nomi_lambda[num_isto_lambda] = {"Rapporto numero "+ tipo2[1]+ " e "+ tipo2[2]+ " dopo applicazione selezioni"}; //"Integrale segnale Lambda/integrale segnale AntiLambda dopo applicazione tagli (rapporto tra numero lambda e antilambda dopo applicazione tagli, cioe' non corretto per efficienza totale)"
  TString molteplicit[mult+1]={"0-5 %", "5-10 %", "10-30 %","30-50 %","50-100 %", "0-100 %" };//"[40-70)"

  if (isHM){
    molteplicit[0] = "0-0.001 %";
    molteplicit[1] = "0.001-0.005 %";
    molteplicit[2] = "0.005-0.01 %";
    molteplicit[3] = "0.01-0.05 %";
    molteplicit[4] = "0.05-0.1 %";
  }

  TString systematic[numsysV0]={"Default", "CosPointingAngleV0", "CosPointingAngleCasc", "DCAV0ToPV", "DCAMesonToPV", "DCABaryonToPV", "DCABachToPV", "DCAV0Daught", "CascDaught", "MassLambda", "Mixed1", "+CascPAngle>0.995", "+V0PAngle>0.99", "+V0PAngle>0.995", "+Both PAngle > 0.995", "+CascPAngle>0.992", "+CascV0<0.6", "+InvMassLambda< 0.005", "+BothTwoBefore", "Before+CascPAngle>0.995)", "Before+CosV0Angle>0.995"};
  // TString systematic[numsysV0]={"Default", "CosPointingAngleV0", "CosPointingAngleCasc"};
  TString titleX="p^{Assoc}_{T} (GeV/c)";
  TString titleY[num_isto+1]={"#sigma_{G} (GeV/c^{2})","#mu (GeV/c^{2})", "#chi^{2}/dof","S/B", "S", "dN/dp_{T} 1/N_{ev}", "S/(S+B)", "","B (integral/conteggi)", "B (integral/conteggi)" , "B integral (side/central)","Efficiency", "Relative error"};
  TString titleY_MC[num_isto_MC]={"#sigma (GeV/c^{2})","#mu (GeV/c^{2})", "#chi^{2}/dof","S/B","S", "S/(S+B)",  "Conteggi", "Numero di "+ tipo2[type]+ "/evento" , "Efficienza x B.R.","Efficienza x accettanza x B.R.", "Efficienza", "F", "Efficienza"}; //r e' rapporto integrale MC/segnale MC
  TString titleY_cfr[num_isto_cfr]={"Numero "+ tipo2[type]+ " dati veri/numero "+ tipo2[type]+ " dati MC", "Numero di "+ tipo2[type]+ "/evento/(efficienza x accettanza x B.R.)"};
  TString titleY_lambda[num_isto_lambda]={"Numero #Lambda/Numero #bar{#Lambda}"};
  TFile *myfile[20];
  TFile *myfileMC[mult+1];
  TString innameMC="FinalOutput/DATA" +year0+"/invmass_distribution_thesis/invmass_distribution";
  //  TString innamedati="invmass_distribution/invmass_distribution_"+tipo[type];
  TString innamedati="FinalOutput/DATA" +year0+"/invmass_distribution_thesis/invmass_distribution";

  TString innamecompleto[20];

  innameMC+="_MCEff";
  if (PtBinning==1){
    innameMC += "_PtBinning1";
    innamedati += "_PtBinning1";
  }
  innameMC+=Path1MC;
  innamedati+=Path1;

  innamedati+="_"+yeardati; 
  innameMC+="_"+yearMC; 
  TString innamedataorMC[2]={innamedati, innameMC};
  cout << "input file name (dati)" << innamedati << endl;
  cout << "input file name (MC)  " << innameMC << endl;
  TString outname="FinalOutput/DATA"+year0 + "/invmass_distribution_thesis/confronto_histogrammi_"+year +"_" +tipobis[type]+Path1 + Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+StringSystV0Analysis[isSystV0Analysis]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f" + SisINEL[isINEL],sysTrigger, sysV0, sys, PtTrigMin);
  TString outnamepdf = outname;
  outname += ".root";

  TFile *outfile = new TFile(outname, "RECREATE");
  cout << "output file name " << outname << endl;

  Int_t Color[21]={1, 401, 801, 628, 909, 881, 860, 868, 841, 418, 628, 909, 881, 867, 921, 401, 841, 862, 866, 865, 864};
  TCanvas *c[num_isto+1]; 
  TCanvas *c_MC[num_isto_MC-num_isto+1];
  TCanvas *c_cfr[num_isto_cfr];
  TCanvas *c_lambda[num_isto_lambda];
  TH1F *histo[num_isto+1][30];
  TH1F *histo_master[num_isto+1];
  TH1F *histo_ratio[num_isto+1][30];
  TH1F *histo_MC[num_isto_MC][mult+1];
  TH1F *histo_cfr[num_isto_cfr][mult+1];
  TH1F *histo_lambda[num_isto_lambda][mult+1];
  TH1F *histo1[num_isto_lambda][mult+1];
  TH1F *histo2[num_isto_lambda][mult+1];
  //TH1F *bontafit2; //rapporto tra integrale del segnale dati veri dopo applicazione tagli/integrale segnale MC dopo applicazione tagli, correttamente normalizzato
  //TH1F *histo_signal_int_corretto; // integrale del segnale dati veri dopo applicazione tagli/numero eventi per classe di molteplicita'/efficienza (totale)
  Float_t X1_MC[num_isto_MC]={0.7, 0.1,0.7, 0.7, 0.7,0.7,  0.7, 0.3, 0.3,0.1,0.7,0.1,  0.1 };
  Float_t Y1_MC[num_isto_MC]={0.7, 0.7,0.7, 0.7, 0.7,0.7,  0.7, 0.1, 0.1,0.7,0.1,0.7,  0.1 };
  Float_t X2_MC[num_isto_MC]={0.9, 0.3,0.9, 0.9, 0.9,0.9,  0.9, 0.5, 0.5,0.3,0.9,0.3,  0.3 };
  Float_t Y2_MC[num_isto_MC]={0.9, 0.9,0.9, 0.9, 0.9,0.9,  0.9, 0.3, 0.3,0.9,0.3,0.9,  0.3 };
  Float_t X1_cfr[num_isto_cfr]={0.7, 0.7 };
  Float_t Y1_cfr[num_isto_cfr]={0.7, 0.7 };
  Float_t X2_cfr[num_isto_cfr]={0.9, 0.9 };
  Float_t Y2_cfr[num_isto_cfr]={0.9, 0.9 };
  Float_t X1_lambda[num_isto_lambda]={0.7 };
  Float_t Y1_lambda[num_isto_lambda]={0.7 };
  Float_t X2_lambda[num_isto_lambda]={0.9 };
  Float_t Y2_lambda[num_isto_lambda]={0.9 };
  Int_t marker[mult+1]={7,20,20,22,29,25};

  THStack *hs[num_isto+1];//("hs","test stacked histograms");  
  THStack *hs_ratio[num_isto+1];//("hs","test stacked histograms");
  THStack *hs_MC[num_isto_MC];
  THStack *hs_cfr[num_isto_cfr];
  THStack *hs_lambda[num_isto_lambda];

  //for first 4 brackets (K0s and Lambda) use this!
  //  Float_t max_range[num_tipo][num_isto+1]= {{0.012, 0.505,8, 500,1,250000,100,100,1},{0.004, 1.119,45, 40,1.2,400000,0.04,60000,1}, {0.004, 1.119,45, 40,1.2,400000,0.04, 60000,1}, {0.003, 1.119,30, 20,1.2, 800000,0.09,60000,1}, {0.01, 1.325,3,20, 0.0001,0.004, 1,1,5,5,0.2}, {0.01, 1.325,3, 20, 0.0001,0.004,1,1,5,5,0.2}, {0.006, 1.678,8, 20, 20,500,1,250000,100,100,1}, {0.012, 1.678,8, 20,20, 500,1,250000,100,100,1} };
  //  Float_t min_range[num_tipo][num_isto+1]= {{0., 0.491,0, 0,0.92,0,0.000001,0,0},{0., 1.112,0, 0,0,0,0.000001,0,0}, {0., 1.112,0, 0,0,0,0.000001,0,0}, {0., 1.112,0, 0,0,0,0.00001,0,0},  {0., 1.32,0, 0,0,0,0.00001,0,0},  {0., 1.32,0, 0,0,0,0.00001,0,0},  {0., 1.665,0, 0,0,0,0.00001,0,0},  {0., 1.665,0, 0,0,0,0.00001,0,0} };
  //for first 4 brackets (K0s and Lambda) use this!

  Float_t min_range[num_tipo][num_isto+1]= {{0., 0.491,0, 0,0,0,0,0.92,0,0.000001,0,0,-0.001},{0., 1.112,0, 0,0,0,0,0,0,0.000001,0,0,-0.001}, {0., 1.112,0, 0,0,0,0,0,0,0.000001,0,0,-0.001}, {0., 1.112,0, 0,0,0,0,0,0.00001,0,0,-0.001},  {0., 1.32,0, 0,0,0,0.75,0.8,0,0.00001,0,0,-0.011},  {0., 1.32,0, 0,0,0,0.75,0.8,0,0.00001,0,0,-0.01},  {0., 1.665,0, 0,0,0,0,0,0,0.00001,0,0,-0.01},  {0., 1.665,0, 0,0,0,0,0,0,0.00001,0,0,-0.01 }, {0., 1.32,0, 0,0,0,0.8,0.8,0,0.00001,0,0,-0.011},  {0., 1.665,0, 0,0,0,0,0,0,0.00001,0,0,-0.01 }};
  Float_t max_range[num_tipo][num_isto+1]= {{0.012, 0.505,8, 500,1,5,1,250000,100,100,2,1, 0.001},{0.01, 1.119,45, 40,367,367,1.2,400000,0.04,60000,2,1,0.01}, {0.01, 1.119,45, 40,367,367,1.2,400000,0.04, 60000,2,1,0.01}, {0.003, 1.119,30, 20,367,367,1.2, 800000,0.09,60000,2,1,0.01}, {0.005, 1.325,3,20, 0.003,0.04, 1,1.5,5,5,2,0.3,0.1}, {0.005, 1.325,3, 20, 0.003,0.04,1,1.5,5,5,2,0.3,0.1}, {0.006, 1.678,8, 20, 20,500,1,250000,100,100,1,1,0.01}, {0.012, 1.678,8, 20,20, 500,1,250000,100,100,1,1,0.01}, {0.005, 1.325,3, 20, 0.003,0.04,1,1.5,5,5,2,0.3,0.1} ,{0.012, 1.678   ,8, 20,20, 500,1,250000,100,100,1,1,0.01}};

  Float_t canvas_size1[num_isto+1]={800, 800, 1400, 1400, 800,800, 1400, 800, 800, 800};
  Float_t canvas_size2[num_isto+1]={600, 600, 500, 500, 600, 1000, 500, 1000, 1000, 1000};

  Float_t massParticle[num_tipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245, 1.32171, 1.67245};
  TString Title[2];

  TF1 *retta = new TF1("retta", "pol0", 0,10);
  retta->SetLineColor(1);
  retta->SetLineWidth(1);
  retta->SetParameter(0,massParticle[type]);
  TF1 *rettaUno = new TF1("rettaUno", "pol0", 0,10);
  rettaUno->SetLineColor(1);
  rettaUno->SetLineWidth(1);
  rettaUno->SetParameter(0,1);


  //confronto con dati pubblicati*******
  cout << "\nprendo histo per confronto con dati pubblicati " << endl;
  TString PathDatiPubblicati ="";
  if (type==0) PathDatiPubblicati = "HEPData-ins1748157-v1-Table_1.root";
  else if (type==8) PathDatiPubblicati = "HEPData-1583750454-v1-Table_3.root";
  TFile *filedatipubblicati = new TFile(PathDatiPubblicati, "");
  if (!filedatipubblicati) {cout << "file dati pubblicati not there " << endl; return;}
  TString STable = "Table 3";
  if (type==0) STable = "Table 1";
  TDirectoryFile *dir = (TDirectoryFile*)filedatipubblicati->Get(STable);
  if (!dir)  {cout << "directory dati pubblicati not there " << endl; return;}

  TH1F* hspectrum[11];
  TH1F* hspectrum1[11];
  TH1F* hspectrum2[11];
  TH1F* hspectrum3[11];
  TH1F* hspectrumetot[11];
  TH1F* hspectrumCfr[6];
  TSpline3 *splineFio[6];

  for (Int_t i=0; i<11; i++){
    hspectrum[i] = (TH1F*)dir->Get(Form("Hist1D_y%i", i+1));
    hspectrum1[i] = (TH1F*)dir->Get(Form("Hist1D_y%i_e1", i+1));
    hspectrum2[i] = (TH1F*)dir->Get(Form("Hist1D_y%i_e2", i+1));
    hspectrum3[i] = (TH1F*)dir->Get(Form("Hist1D_y%i_e3", i+1));

    if (!hspectrum[i] ||     !hspectrum1[i] || !hspectrum2[i]|| !hspectrum3[i] ) { cout << "histo is missing " << endl; return;}
    hspectrumetot[i]= (TH1F*)    hspectrum3[i]->Clone(Form("Hist1D_y%i_etot", i+1));

    //    cout << " I clone a histo " << endl;
    if (i<=mult)    hspectrumCfr[i] = (TH1F*)hspectrum[i]->Clone(Form("Spectrum_mult%i",i+1));
    //    cout << " I clone a histo " << endl;
    for (Int_t b=1; b<= hspectrumCfr[0]->GetNbinsX();b++){
      //      hspectrumetot[i]->SetBinContent(b,   sqrt(pow(hspectrum1[i]->GetBinContent(b),2) + pow(hspectrum2[i]->GetBinContent(b),2) + pow(hspectrum3[i]->GetBinContent(b),2)) );
      hspectrumetot[i]->SetBinContent(b,   sqrt(pow(hspectrum1[i]->GetBinContent(b),2) + pow(hspectrum2[i]->GetBinContent(b),2) ) );//hspectrum3 contains the uncorrelated part of the syst error, hspectrum2 the total one
    }
  }

  //  cout << " I got all histos " << endl;

  for (Int_t b=1; b<= hspectrumCfr[0]->GetNbinsX();b++){
    hspectrumCfr[0]->SetBinContent(b,1./5 * (hspectrum[0]->GetBinContent(b)+ hspectrum[1]->GetBinContent(b)*4));
    hspectrumCfr[0]->SetBinError(b, sqrt(pow(hspectrumetot[0]->GetBinContent(b) * 1./5,2) + pow(hspectrumetot[1]->GetBinContent(b)*4./5,2)));

    hspectrumCfr[1]->SetBinContent(b,hspectrum[2]->GetBinContent(b));
    hspectrumCfr[1]->SetBinError(b, hspectrumetot[2]->GetBinContent(b));
    
    hspectrumCfr[2]->SetBinContent(b,1./20 * (hspectrum[3]->GetBinContent(b)*5+ hspectrum[4]->GetBinContent(b)*5 + hspectrum[5]->GetBinContent(b)*10));
    hspectrumCfr[2]->SetBinError(b, sqrt(pow(hspectrumetot[3]->GetBinContent(b)*5./20,2) + pow(hspectrumetot[4]->GetBinContent(b)*5./20,2) + pow(hspectrumetot[5]->GetBinContent(b)*10./20,2)));

    hspectrumCfr[3]->SetBinContent(b,1./20 * (hspectrum[6]->GetBinContent(b)*10+ hspectrum[7]->GetBinContent(b)*10));
    hspectrumCfr[3]->SetBinError(b, sqrt(pow(hspectrumetot[6]->GetBinContent(b) * 10./20,2) + pow(hspectrumetot[7]->GetBinContent(b)*10./20,2)));

    hspectrumCfr[4]->SetBinContent(b,1./50 * (hspectrum[8]->GetBinContent(b)*20+ hspectrum[9]->GetBinContent(b)*30));
    hspectrumCfr[4]->SetBinError(b, sqrt(pow(hspectrumetot[8]->GetBinContent(b) * 20./50,2) + pow(hspectrumetot[9]->GetBinContent(b)*30./50,2)));

    hspectrumCfr[5]->SetBinContent(b,hspectrum[10]->GetBinContent(b));
    hspectrumCfr[5]->SetBinError(b, hspectrumetot[10]->GetBinContent(b));

  }
  cout << " \n\n*************\nspline of spectra obtained from Fiorella's ones and evaluation of their mean relative error " << endl;
  Float_t  MeanRelErr[6]={0};
  Float_t  SigmaRelErr[6]={0};
  for (Int_t i=0; i<6; i++){
    if(type==4 || type==5)    hspectrumCfr[i]->Scale(1./2.);
    splineFio[i] = new TSpline3(hspectrumCfr[i],Form("splineFio_%i",i));
    MeanRelErr[i]=0;
    SigmaRelErr[i]=0;

    for(Int_t b=1; b <= hspectrumCfr[i]->GetNbinsX(); b++){
      MeanRelErr[i]+=(hspectrumCfr[i]->GetBinError(b)/hspectrumCfr[i]->GetBinContent(b));
    }
    MeanRelErr[i]= MeanRelErr[i]/hspectrumCfr[i]->GetNbinsX();
    for(Int_t b=1; b <= hspectrumCfr[i]->GetNbinsX(); b++){
      SigmaRelErr[i]+=pow(MeanRelErr[i]-(hspectrumCfr[i]->GetBinError(b)/hspectrumCfr[i]->GetBinContent(b)),2);
    }
    SigmaRelErr[i]=sqrt( SigmaRelErr[i]/(hspectrumCfr[i]->GetNbinsX()-1));
    cout << "\n\nmean relative error for spectrum at mult " << molteplicit[i] << ": " << MeanRelErr[i]<< " sigma of distribution " << SigmaRelErr[i] <<  endl;
  }
  //************************************

  for (Int_t j =0 ; j< num_isto_finale; j++){ //ciclo sui diversi istogrammi relativi ai dati veri

    //    cout << "tipo isto " << j << endl;
    nomi_completi[j] = nomi[j];
    //    if (isSystV0Analysis && (tipo_histo[j]!="Efficiency" && tipo_histo[j]!="SSB")) continue;
    //       if (!isSystV0Analysis && tipo_histo[j]!="SEffCorr") continue;
    //  if (j!=num_isto_finale-1) continue;
    //    if (j==3 || j==10) continue;
    cout << "\n\ntipo isto " << titleY[j]<< endl;
    //    auto legend = new TLegend(X1[j],Y1[j],X2[j],Y2[j]);
    //    auto legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    auto legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    if (!isSystV0Analysis)    legend->SetHeader("Multiplicity classes"); 
    else     legend->SetHeader("Selections"); 
    auto legendMass = new TLegend(0.6, 0.6, 0.9, 0.9);
    if (!isSystV0Analysis)    legendMass->SetHeader("Multiplicity classes"); 
    else     legendMass->SetHeader("Selections"); 
    //    c[j] = new TCanvas (Form("c%i",j),nomi[j], 1800,1000);//1400, 500 per averle in orizzontale oppure 800, 1000 per averle in verticale
    c[j] = new TCanvas (nomi[j],nomi[j], 1800,1000);
    //    c[j]->Divide(4,2); for comparison between particle and antiparticle
    c[j]->Divide(2,2);
    Int_t counterParticleType=0;
    for (Int_t ParticleType=type; ParticleType<=type/*+1*/; ParticleType++){
      cout <<"Particle Type " <<  tipo[ParticleType] << endl; 
      counterParticleType++; 
     
      if (isINEL){
	if (type==0){
	  if (tipo_histo[j]=="S")	  max_range[ParticleType][j] = 0.3;
	  if (tipo_histo[j]=="SSB")	 { min_range[ParticleType][j] = 0.85;  max_range[ParticleType][j] = 1.15; }
	}
      }
    
      for(Int_t dataorMC=0; dataorMC<2; dataorMC++){
	if (isINEL && dataorMC ==1) continue;
	cout << "data or MC " << dataorMC << endl;
	cout << num_isto_finale << endl;
	if (dataorMC==0 && j>=num_isto_finale-6 && tipo_histo[j]!="Efficiency") {cout << " I continue " << endl; continue;}
	if (dataorMC==1 && tipo_histo[j] == "EfficiencyRelError") continue;
	if (isINEL){
	if (dataorMC==0 && tipo_histo[j]=="SEffCorr") continue;
	if (dataorMC==0 && tipo_histo[j]=="Efficiency") continue;
	}
	hs[j] = new THStack(Form("hs%i",j),"");
	hs_ratio[j] = new THStack(Form("hs_ratio%i",j),"");
	Title[ParticleType-type]=isDataorMC[dataorMC]+"  "+ tipobis[ParticleType];
	Int_t ForLoopLimit=0;
	Int_t ForLoopLimInf=0;
	//	cout << " before loop " << endl;
	if (isSystV0Analysis){ ForLoopLimit=numsysV0; ForLoopLimInf=0;}
	if (!isSystV0Analysis) {ForLoopLimit=0;ForLoopLimInf=mult;}
	Int_t 	Icounter=-1;
	for (Int_t i =ForLoopLimInf ; i>= ForLoopLimit; i--){//ok for mult loop
	  cout << "mult loop " << i << endl;
        //for (Int_t i =ForLoopLimInf ; i< ForLoopLimit; i++){//ok for systloop
	  //cout << "\n\nsystematic: " << systematic[i] << endl;
	     //	     if (isSystV0Analysis && i>=3 && i<=6 ) continue;
	     //	     if (isSystV0Analysis && (i!=0 && i!=19)) continue;
	    //	    if (isSystV0Analysis && (i<19 || i>20) && i!=0) continue;
	  //	  if (isSystV0Analysis && i!=10 && i!=11 && i!=19 && i!=20 && i!=0) continue;
	    //	    if (isSystV0Analysis && (i<10 || i>20)) continue;
	    //	    if (isSystV0Analysis && systematic[i]=="Mixed3") continue;
	    //	    if (isSystV0Analysis && systematic[i]=="Mixed4") continue;
	  Icounter++;
	  //if multiplicity dependence is analyzed*****************
	  if (!isSystV0Analysis){
	    innamecompleto[i]= Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType] +Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f" +SisINEL[isINEL] +".root", i, sysTrigger, sysV0, sys, PtTrigMin);

	  }
	  cout <<"fileinput completo " <<  innamecompleto[i]<< endl;
	  //*******************************************************

	  //if particle selection dependence is analyzed*****************
	  if (isSystV0Analysis){
	    innamecompleto[i]= Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+Srap[rap]+SSkipAssoc[SkipAssoc]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", chosenmolt, sysTrigger, i, sys, PtTrigMin);

	  }
	 
	  //*******************************************************
	  myfile[i] = new TFile(innamecompleto[i], ""); 
	  histo[j][i] = (TH1F*)myfile[i]->Get("histo_"+tipo_histo[j]);
	  if (!histo[j][i]) {cout << " histo not here " << j << " " << tipo_histo[j] << endl; return;}
	  if (i<=25) 	  histo[j][i]->SetMarkerStyle(33);
	  else 	  histo[j][i]->SetMarkerStyle(20);
	  histo[j][i]->SetMarkerSize(1.2);
	  histo[j][i]->SetLineColor(Color[i]);
	  histo[j][i]->SetMarkerColor(Color[i]);
	  if (tipo_histo[j]=="SEffCorr" && !isSystV0Analysis){
	    hspectrumCfr[i]->SetMarkerStyle(27);
	    hspectrumCfr[i]->SetLineColor(Color[i]);	   
	    hspectrumCfr[i]->SetMarkerColor(Color[i]);	   
	    splineFio[i] ->SetLineColor(Color[i]);
	    hs[j]->Add(hspectrumCfr[i], "");
	    //hs[j]->Add(    splineFio[i], "");
	  }
	  if (tipo_histo[j]=="SEffCorr" && isSystV0Analysis && i==0){
	    hspectrumCfr[5]->SetMarkerStyle(27);
	    hspectrumCfr[5]->SetLineColor(Color[i]);	   
	    hspectrumCfr[5]->SetMarkerColor(Color[i]);	   
	      hs[j]->Add(hspectrumCfr[5], "");
	  }

	  hs[j]->Add(histo[j][i], "p");
	  hs[j]->Add(histo[j][i], "same");
	  
	  if (tipo_histo[j]=="EfficiencyRelError"){
	    for (Int_t l=0; l<histo[j][i]->GetNbinsX() ; l++){
	      //	      cout << "bin " << l<< "  " <<histo[j][i]->GetBinContent(l+1) << endl;
	    }
	  }
	  if (Icounter==0) 	{
	    histo_master[j] = (TH1F*)histo[j][i] ->Clone("histo_"+tipo_histo[j]+"_master");
	    histo_master[j]->Sumw2();
	  }
	  histo_ratio[j][i] = (TH1F*)	histo[j][i] ->Clone("histo_"+tipo_histo[j]+"_ratio");
	  histo_ratio[j][i]->Sumw2();
	  
	  Bool_t SkipIsto=kFALSE;
	  for (Int_t IPt=1; IPt <histo[j][i]->GetNbinsX(); IPt++){ //parto da IPt = 1 perché nelle Xi il primo bin è vuoto (< 0.3 GeV/c)
	    //	    cout << IPt << " " << SkipIsto << endl;
	    if (histo_master[j]->GetBinContent(IPt+1)==0) SkipIsto=kTRUE;
	  }
	  
	  //cout << "SkipIsto " << SkipIsto << endl;
	 	 
	  if (!SkipIsto)	{
	    histo_ratio[j][i] ->Divide( histo_master[j] );
	  }

	    if (!isSystV0Analysis){
	      //	      cout << "\n\n\n ****************************" << endl;
	      //	      cout << " multiplicity: " << i << endl;
	      if (tipo_histo[j]=="SEffCorr"){
		for(Int_t b=1; b <=  histo[j][i]->GetNbinsX(); b++){
		  histo_ratio[j][i]->SetBinContent(b, histo[j][i]->GetBinContent(b)/splineFio[i]->Eval(histo[j][i]->GetBinCenter(b)));
		  histo_ratio[j][i]->SetBinError(b, sqrt(pow(histo[j][i]->GetBinError(b),2) +pow(MeanRelErr[i]* histo[j][i]->GetBinContent(b)*histo_ratio[j][i]->GetBinContent(b),2)) /splineFio[i]->Eval(histo[j][i]->GetBinCenter(b)));

		  cout << "\n bin center " << histo[j][i]->GetBinCenter(b) << " number of sigmas from 1 "  << TMath::Abs(histo_ratio[j][i]->GetBinContent(b)-1)/(histo_ratio[j][i]->GetBinError(b)) << endl;
		  cout << "bin content " << histo_ratio[j][i]->GetBinContent(b) << " bin error " << histo_ratio[j][i]->GetBinError(b) << endl;
		  //  cout << "\n bin center " << histo[j][i]->GetBinCenter(b) << " histo value at bin center " << histo[j][i]->GetBinContent(b) << " spline value " << splineFio[i]->Eval(histo[j][i]->GetBinCenter(b)) << " ratio value 1 " <<histo[j][i]->GetBinContent(b)/splineFio[i]->Eval(histo[j][i]->GetBinCenter(b))<< " , 2: " <<  histo_ratio[j][i]->GetBinContent(b) <<endl;
		}
	      }
	    }

	  //I set at zero the first bin [0-0.3] GeV/c since not used
	    if (type==8){
	  histo_ratio[j][i]->SetBinContent(1,0);
	  histo[j][i]->SetBinContent(1,0);
	    }
	  //I set the correct errors to the ratios according to Barlow prescription; I don't do this for SEffCorr since in that case the ratios are mycorrectedspectra/Fiorella'spectra spline. The error for that ratio is correctly computed just above.
	  if (!(tipo_histo[j]=="SEffCorr")){
	  for (Int_t l=1; l<histo_ratio[j][i]->GetNbinsX(); l++){
	    histo_ratio[j][i]->SetBinError(l+1, sqrt(TMath::Abs(pow(histo_master[j]->GetBinError(l+1),2) -pow(histo[j][i]->GetBinError(l+1),2)))/histo_master[j]->GetBinContent(l+1));
	  }
	  }

	  if (i<=25)	  histo_ratio[j][i]->SetMarkerStyle(33);
	  else   histo_ratio[j][i]->SetMarkerStyle(20);
	  histo_ratio[j][i]->SetMarkerSize(1.2);
	  histo_ratio[j][i]->SetLineColor(Color[i]);
	  histo_ratio[j][i]->SetMarkerColor(Color[i]);
	  
	  if (Icounter!=0)	hs_ratio[j]->Add(histo_ratio[j][i], "P");
	  if (tipo_histo[j]=="SEffCorr") hs_ratio[j]->Add(histo_ratio[j][i], "P");
	  if (Icounter!=0 && tipo_histo[j]!="EfficiencyRelError")	  hs_ratio[j]->Add(histo_ratio[j][i]);
	  //	  cout << "ho aggiunto histo ratio" << endl;
	  if (dataorMC==0 && !isSystV0Analysis && counterParticleType==1) legend->AddEntry(histo[j][i],molteplicit[i],"pel");
	  if (dataorMC==0 && !isSystV0Analysis && counterParticleType==1) legendMass->AddEntry(histo[j][i],molteplicit[i],"pel");
	  //if particle selection dependence is analyzed*****************
	  if (dataorMC==0 && isSystV0Analysis && counterParticleType==1 && i!=0) legend->AddEntry(histo[j][i],systematic[i],"pel");
	  if (dataorMC==0 && isSystV0Analysis && counterParticleType==1 && i!=0) legendMass->AddEntry(histo[j][i],systematic[i],"pel");
	  if(i==0 && j==1 && dataorMC==0 && counterParticleType==1) legendMass->AddEntry(retta, tipoNotSign[ParticleType]+" PDG mass");
	  outfile->cd();
	  retta->SetLineColor(Color[i]);
	  //	  histo_ratio[j][i]->Fit(retta, "B+");
	  rettaUno->FixParameter(0,1);
	  //	  histo_ratio[j][i]->Fit(rettaUno, "B+");
	  // histo[j][i]->Write();
	  if (tipo_histo[j]=="SEffCorr") histo_ratio[j][i]->Write();
	  //	  cout << "I'm triying to close the file" << Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+ "_" +MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", chosenmolt, sysTrigger, i, sys, PtTrigMin) << endl;
	  //       	myfile[i]->Close();      
	  //	  cout << "file closed" << endl;
	} //molteplicità o sistematico
      
	if (ParticleType==type)      c[j]->cd(dataorMC+1);
	if (ParticleType==type+1)      c[j]->cd(dataorMC+3);
	//     if(j == 6) {gPad->SetLogy();}
	gPad->SetLeftMargin(0.2);
	hs[j]->Draw("nostack"); //cosi' me li disegna in maniera corretta
	hs[j]->SetTitle(Title[ParticleType-type]);
	hs[j]->GetXaxis()->SetTitle(titleX); //il set degli assi va fatto dopo l'opzione Draw!!
	hs[j]->GetXaxis()->SetTitleOffset(0.9);
	hs[j]->GetXaxis()->SetTitleSize(0.05);
	hs[j]->GetYaxis()->SetTitle(titleY[j]);
	hs[j]->GetYaxis()->SetTitleOffset(2);
	hs[j]->GetYaxis()->SetTitleSize(0.046); 
	hs[j]->SetMinimum(min_range[ParticleType][j]);
	hs[j]->SetMaximum(max_range[ParticleType][j]); //l'unico modo funzionante per settare estremi

	if (j==1)      legendMass->Draw();
	else legend->Draw();
	if (j==1)    retta->Draw("same");
	//	cout << "ok 2"<< endl;
	if (ParticleType==type)      c[j]->cd(dataorMC+1+2);
	if (ParticleType==type+1)      c[j]->cd(dataorMC+3+2);
	//     if(j == 6) {gPad->SetLogy();}
	//	cout << "ok 3"<< endl;
	gPad->SetLeftMargin(0.2);
	cout <<  " draw ratio stack for data/MC" << dataorMC << " and type " <<  tipo_histo[j] << endl; 
	hs_ratio[j]->Draw("nostack"); //cosi' me li disegna in maniera piu' corretta
	//	cout << "ok 4"<< endl;
	hs_ratio[j]->SetTitle(Title[ParticleType-type]);
	//	cout << "ok 4"<< endl;
	hs_ratio[j]->GetXaxis()->SetTitle(titleX); //il set degli assi va fatto dopo l'opzione Draw!!
	hs_ratio[j]->GetXaxis()->SetTitleOffset(0.9);
	hs_ratio[j]->GetXaxis()->SetTitleSize(0.05);
	hs_ratio[j]->GetYaxis()->SetTitle(titleY[j]);
	hs_ratio[j]->GetYaxis()->SetTitleOffset(1.2);
	hs_ratio[j]->GetYaxis()->SetTitleSize(0.046); 
	hs_ratio[j]->SetMinimum(0);
	hs_ratio[j]->SetMaximum(2); //l'unico modo funzionante per settare estremi
	if (j==4 || j ==5 || j == num_isto_finale -1){
	  hs_ratio[j]->SetMaximum(13);
	}
	if (j ==5){
	  hs_ratio[j]->SetMaximum(3);
	  hs_ratio[j]->SetMinimum(0);
	}
	//	if (isSystV0Analysis && tipo_histo[j]=="SSB"){
	if (tipo_histo[j]=="SSB"){
	  hs_ratio[j]->SetMinimum(0.95);
	  hs_ratio[j]->SetMaximum(1.05); //1.8
	}
	if (isSystV0Analysis && tipo_histo[j]=="Efficiency"){
	  hs_ratio[j]->SetMinimum(0.6);
	  hs_ratio[j]->SetMaximum(1);
	}
	//	cout << "ok 5"<< endl;
	if (j==1){
	 hs_ratio[j]->SetMinimum(0.9985);
	 hs_ratio[j]->SetMaximum(1.0015);
	}
	//cout << "ok 6"<< endl;
	//	legend->Draw();
	rettaUno->Draw("same");
	//cout << "ok 7"<< endl;
      } //data or MC
    } //particle type
  } //tipo isto

  for (Int_t j =0 ; j< num_isto_finale; j++){
    if (j==0) c[j]->SaveAs(outnamepdf+".pdf(");
    else if (j==num_isto_finale-1) c[j]->SaveAs(outnamepdf+".pdf)");
    else c[j]->SaveAs(outnamepdf+".pdf");
    //    if (isSystV0Analysis && (tipo_histo[j]!="Efficiency" && tipo_histo[j]!="SSB")) continue;
    //    if (!isSystV0Analysis && tipo_histo[j]!="SEffCorr" && tipo_histo[j]!="mean") continue;
       outfile->WriteTObject(c[j]);
  }

  /*
    for (Int_t j =0 ; j<= num_isto_MC-num_isto; j++){ 
    c_MC[j]->SaveAs("isto_confronto/confronto_"+ tipo_histo_MC[num_isto+j-1]+ def[cut]+"_"+tipo[type]+ "_mixed.pdf");
    }
  */
  /*
    for (Int_t j =0 ; j<num_isto_cfr; j++){ //ciclo sugli istogrammi di confronto dati/MC
    hs_cfr[j] = new THStack(Form("hs%i_cfr",j), nomi_cfr[j]);
    auto legend = new TLegend(X1_cfr[j],Y1_cfr[j],X2_cfr[j],Y2_cfr[j]);
    legend->SetHeader("Intervalli di molteplicita'"); 
    c_cfr[j] = new TCanvas (Form("c%i_cfr",j),nomi_cfr[j], 800, 600);
    if(j == 1) {c_cfr[j]->SetLogy();}
    for (Int_t i =0 ; i< mult+1; i++){
    myfile[i] = new TFile(Form(innamedati+"_molt%i_SysT%i_SysV0%i_Sys%i.root", i, syst, sysV0, sys), ""); 
    histo[j][i] = (TH1F*)myfile[i]->Get("histo_signal_int");
    myfileMC[i] = new TFile(Form(innameMC+"_%i" +def[cut]+".root", i), ""); 
    histo_MC[j][i] = (TH1F*)myfileMC[i]->Get("histo_"+tipo_histo_MC[j+6]);
    histo[j][i] ->Sumw2();
    histo_MC[j][i]->Sumw2();
    histo[j][i] -> Divide(histo_MC[j][i]);
    histo_cfr[j][i]=(TH1F*)histo[j][i]->Clone("");
    histo_cfr[j][i]->SetMarkerStyle(marker[i]);
    if (i+1 != 5) {histo_cfr[j][i]->SetLineColor(i+1); histo_cfr[j][i]->SetMarkerColor(i+1);}
    if (i+1 == 5) {histo_cfr[j][i]->SetLineColor(6);   histo_cfr[j][i]->SetMarkerColor(6);}
    if (i+1 == 3) {histo_cfr[j][i]->SetLineColor(8);   histo_cfr[j][i]->SetMarkerColor(8);}
    hs_cfr[j]->Add(histo_cfr[j][i]);
    legend->AddEntry(histo_cfr[j][i],molteplicit[i],"pel");
    outfile->cd();
    histo_cfr[j][i]->Write();
    }
    c_cfr[j]->cd();
    hs_cfr[j]->Draw("nostack"); //cosi' me li disegna in maniera piu' corretta
    hs_cfr[j]->GetXaxis()->SetTitle(titleX); //il set degli assi va fatto dopo l'opzione Draw!!
    hs_cfr[j]->GetXaxis()->SetTitleOffset(1.2);
    hs_cfr[j]->GetYaxis()->SetTitle(titleY_cfr[j]);
    hs_cfr[j]->GetYaxis()->SetTitleOffset(1.3);
    legend->Draw();
    }
  
    for (Int_t j =0 ; j<num_isto_cfr; j++){
    c_cfr[j]->SaveAs(Form("isto_confronto/confronto_cfr_%i"  +def[cut]+"_"+tipo[type]+".pdf", j));
    }
  */
  /*
    if (type==1){
    for (Int_t j =0 ; j<num_isto_lambda; j++){ //confronto tra lambda e antilambda
    hs_lambda[j] = new THStack(Form("hs%i_lambda",j), nomi_lambda[j]);
    auto legend = new TLegend(X1_lambda[j],Y1_lambda[j],X2_lambda[j],Y2_lambda[j]);
    legend->SetHeader("Intervalli di molteplicita'"); 
    c_lambda[j] = new TCanvas (Form("c%i_lambda",j),nomi_lambda[j], 800, 600);
    for (Int_t i =0 ; i< mult+1; i++){
    myfile[i] = new TFile(Form("invmass_distribution/invmass_distribution_"+tipo[type]+"_%i.root", i), "");
    histo1[j][i] = (TH1F*)myfile[i]->Get("histo_"+tipo_histo[j+5]);
    myfile[i] = new TFile(Form("invmass_distribution/invmass_distribution_"+tipo[type+1]+"_%i.root", i), "");
    histo2[j][i] = (TH1F*)myfile[i]->Get("histo_"+tipo_histo[j+5]);
    histo1[j][i] ->Sumw2();
    histo2[j][i] ->Sumw2();
    histo1[j][i] -> Divide(histo2[j][i]);
    histo_lambda[j][i]=(TH1F*)histo1[j][i]->Clone("");
    cout << " bin error"<<histo_lambda[j][i]->GetBinError(1) << "sqrt of bin content " <<sqrt(histo_lambda[j][i]->GetBinContent(1))<< endl;;
    histo_lambda[j][i]->SetMarkerStyle(marker[i]);
    if (i+1 != 5) {histo_lambda[j][i]->SetLineColor(i+1); histo_lambda[j][i]->SetMarkerColor(i+1);}
    if (i+1 == 5) {histo_lambda[j][i]->SetLineColor(6);   histo_lambda[j][i]->SetMarkerColor(6);}
    if (i+1 == 3) {histo_lambda[j][i]->SetLineColor(8);   histo_lambda[j][i]->SetMarkerColor(8);}
    hs_lambda[j]->Add(histo_lambda[j][i]);
    cout << "ciao" << endl;
    legend->AddEntry(histo_lambda[j][i],molteplicit[i],"pel");
    outfile->cd();
    cout << "ciao" << endl;
    histo_lambda[j][i]->Write();
    cout << "ciao " << endl;
    }
    cout << "ciao" << endl;
    c_lambda[j]->cd();
    hs_lambda[j]->Draw("nostack"); //cosi' me li disegna in maniera piu' corretta
    hs_lambda[j]->GetXaxis()->SetTitle(titleX); //il set degli assi va fatto dopo l'opzione Draw!!
    hs_lambda[j]->GetXaxis()->SetTitleOffset(1.2);
    hs_lambda[j]->GetYaxis()->SetTitle(titleY_lambda[j]);
    hs_lambda[j]->GetYaxis()->SetTitleOffset(1.3);
    legend->Draw();
    //c_lambda[j]->Write();
    c_lambda[j]->SaveAs("isto_confronto/confronto_lambda_con_antilambda.pdf");
    }
    // c_lambda[0]->SaveAs("isto_confronto/confronto_lambda_con_antilambda.pdf");
    } 
  */
  outfile->Close();
	Int_t ForLoopLimit=0;
	Int_t ForLoopLimInf=0;
	//	cout << " before loop " << endl;
	if (isSystV0Analysis){ ForLoopLimit=numsysV0; ForLoopLimInf=0;}
	if (!isSystV0Analysis) {ForLoopLimit=0;ForLoopLimInf=mult;}

  for (Int_t i =ForLoopLimInf ; i>= ForLoopLimit; i--){//ok for mult loop
  cout <<"fileinput completo " <<  innamecompleto[i]<< endl;
  }
  cout << "\n\n\n il file di output si chiama " << outname << " (anche .pdf)" << endl;
  cout << " errors on the ratios calculated according to Barlow presciption " << endl;


}
