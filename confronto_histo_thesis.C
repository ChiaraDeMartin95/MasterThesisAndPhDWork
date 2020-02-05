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
void confronto_histo_thesis(Bool_t isSignalFromIntegral=0, Bool_t isBkgParab=0, Bool_t isMeanFixedPDG=1, Float_t PtTrigMin=3.0, Int_t type=4, Int_t sysTrigger=0, Int_t sysV0=0, Int_t sys=0, TString yearMC="2018f1_extra", TString yeardati="2016k", TString year="2016k18f1_extra", TString year0="2016", Bool_t isParticleStudy=1, Bool_t isSystV0Analysis=0, Int_t rap=0){ 


  //isSystV0Analysis=1 if I want to compare different selections used to identify the particle for a given multiplcity (to be chosen!); =0 if I want to compare diffrerent multiplcity ranges (but same selections)
  Int_t chosenmolt=5; //choose the multiplicity you desire
  const Int_t numsysV0=10; 
  const Int_t mult=5; //numero intervalli molteplicita'
  const Int_t num_tipo=8; //numero tipologie di particelle
  const Int_t num_isto=10; //numero di istogrammi di tipologia diversa per dati veri
  const Int_t num_isto_MC=13; //numero di istogrammi di tipologia diversa per MC
  const Int_t num_isto_cfr=2; //numero di istogrammi di confronto dati/MC
  const Int_t num_isto_lambda=1; //numero di istogrammi di confronto lambda/antilambda
  const Int_t num_tagli=3;

  TString isDataorMC[2]={"Data", "MC"};
  TString StringSystV0Analysis[2]={"", "_SystV0Analysis"};
  if (isParticleStudy==0 && numsysV0!=1) return;
  TString BkgType[2]={"BkgRetta", "BkgParab"};
  TString MassFixedPDG[2]={"", "isMeanFixedPDG_"};
  Int_t num_isto_finale=num_isto;
  if (isParticleStudy) num_isto_finale=num_isto+1;
  TString tipo2[num_tipo]={"K^{0}_{s}", "#Lambda","#bar{#Lambda}", "#Lambda + #bar{#Lambda}", "#Xi^{-}", "#Xi^{+}", "#Omega^{-}", "#Omega^{+}"};
  TString tipoNotSign[num_tipo]={"K^{0}_{s}", "#Lambda","#bar{#Lambda}", "#Lambda + #bar{#Lambda}", "#Xi", "#Xi", "#Omega", "#Omega"};
  TString tipo[num_tipo]={"K0s", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos"};
  TString tipobis[num_tipo]={"kK0s", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos"};
  TString ParticleType[num_tipo]={"", "Lambda", "AntiLambda", "LambdaAntiLambda", "XiNeg", "XiPos", "OmegaNeg", "OmegaPos"};
  TString tipo_histo[num_isto+1]={"sigma", "mean", "chis", "SB", "S","SEffCorr", "SSB","FracStrangePeak", "BRatioCentral", "BRatioSide", "Efficiency"};
  TString nomi[num_isto+1]={"Invariant mass resolution", 
			    tipo2[type]+" mass",
			    "#chi^{2}/dof", 
			    "S/B within 3#sigma", 
			    "S within 3#sigma",
			    "dN/dp_{T} 1/N_{ev}",
			    "S/(S+B) within 3#sigma", 
			    "Estimate of gaussian bkg/all selected within 3#sigma",
			    "B within 3#sigma (bkg integral/MC truth)", 
			    "B in sidebands  (bkg integral/MCtruth)",
			    "Efficiency"};
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
  TString systematic[numsysV0]={"Default", "CosPointingAngleV0", "CosPointingAngleCasc", "DCAV0ToPV", "DCAMesonToPV", "DCABaryonToPV", "DCABachToPV", "DCAV0Daught", "CascDaught", "MassLambda"};
  // TString systematic[numsysV0]={"Default", "CosPointingAngleV0", "CosPointingAngleCasc"};
  TString titleX="p^{Assoc}_{T} (GeV/c)";
  TString titleY[num_isto+1]={"#sigma_{G} (GeV/c^{2})","#mu (GeV/c^{2})", "#chi^{2}/dof","S/B", "S", "dN/dp_{T} 1/N_{ev}", "S/(S+B)", "","B (integral/conteggi)", "B (integral/conteggi)" , "Efficiency"};
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
  if (type>=4){
    innameMC+="_Cascades";
    innamedati+="_Cascades";
  }
  innamedati+="_"+yeardati; 
  innameMC+="_"+yearMC; 
  TString innamedataorMC[2]={innamedati, innameMC};
  // cout << "input file name (dati)" << innamedati << endl;
  // cout << "input file name (MC)  " << innameMC << endl;
  TString outname="FinalOutput/DATA"+year0 + "/invmass_distribution_thesis/confronto_histogrammi_"+year +"_" +ParticleType[type]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+StringSystV0Analysis[isSystV0Analysis]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f.root",sysTrigger, sysV0, sys, PtTrigMin);
  if (type>=4) outname="FinalOutput/DATA"+year0 + "/invmass_distribution_thesis/confronto_histogrammi_"+year +"_" +ParticleType[type]+"_"+MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+StringSystV0Analysis[isSystV0Analysis]+Form("_SysT%i_SysV0%i_Sys%i_PtMin%.1f_Rap%i.root",sysTrigger, sysV0, sys, PtTrigMin, rap);
  TFile *outfile = new TFile(outname, "RECREATE");

  Int_t Color[10]={1, 401, 801, 628, 909, 881, 860, 868, 841, 418};
  TCanvas *c[num_isto+1]; 
  TCanvas *c_MC[num_isto_MC-num_isto+1];
  TCanvas *c_cfr[num_isto_cfr];
  TCanvas *c_lambda[num_isto_lambda];
  TH1F *histo[num_isto+1][10];
  TH1F *histo_master[num_isto+1];
  TH1F *histo_ratio[num_isto+1];
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

  Float_t min_range[num_tipo][num_isto+1]= {{0., 0.491,0, 0,0,0,0.92,0,0.000001,0,0},{0., 1.112,0, 0,0,0,0,0,0.000001,0,0}, {0., 1.112,0, 0,0,0,0,0,0.000001,0,0}, {0., 1.112,0, 0,0,0,0,0,0.00001,0,0},  {0., 1.32,0, 0,0,0,0,0,0.00001,0,0},  {0., 1.32,0, 0,0,0,0,0,0.00001,0,0},  {0., 1.665,0, 0,0,0,0,0,0.00001,0,0},  {0., 1.665,0, 0,0,0,0,0,0.00001,0,0} };
  Float_t max_range[num_tipo][num_isto+1]= {{0.012, 0.505,8, 500,367,367,1,250000,100,100,1},{0.004, 1.119,45, 40,367,367,1.2,400000,0.04,60000,1}, {0.004, 1.119,45, 40,367,367,1.2,400000,0.04, 60000,1}, {0.003, 1.119,30, 20,367,367,1.2, 800000,0.09,60000,1}, {0.01, 1.325,3,20, 0.0001,0.004, 1,1,5,5,0.2}, {0.01, 1.325,3, 20, 0.0001,0.004,1,1,5,5,0.2}, {0.006, 1.678,8, 20, 20,500,1,250000,100,100,1}, {0.012, 1.678,8, 20,20, 500,1,250000,100,100,1} };

  Float_t canvas_size1[num_isto+1]={800, 800, 1400, 1400, 800,800, 1400, 800, 800};
  Float_t canvas_size2[num_isto+1]={600, 600, 500, 500, 600, 1000, 500, 1000, 1000};

  Float_t massParticle[num_tipo]= {0.497611, 1.115683, 1.115683, 1.115683, 1.32171, 1.32171, 1.67245, 1.67245};
  TString Title[2];

  TF1 *retta = new TF1("retta", "pol0", 0,10);
  retta->SetLineColor(1);
  retta->SetLineWidth(1);
  retta->SetParameter(0,massParticle[type]);
  TF1 *rettaUno = new TF1("rettaUno", "pol0", 0,10);
  rettaUno->SetLineColor(1);
  rettaUno->SetLineWidth(1);
  rettaUno->SetParameter(0,1);

  for (Int_t j =0 ; j< num_isto_finale; j++){ //ciclo sui diversi istogrammi relativi ai dati veri
    //    cout << "tipo isto " << j << endl;
    nomi_completi[j] = nomi[j];
    //    if (j==3 || j==10) continue;
    cout << "\n\ntipo isto " << titleY[j]<< endl;
    //    auto legend = new TLegend(X1[j],Y1[j],X2[j],Y2[j]);
    auto legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    //    auto legend = new TLegend(0.6, 0.1, 0.9, 0.4);
    if (!isSystV0Analysis)    legend->SetHeader("Multiplicity classes"); 
    else     legend->SetHeader("Selections"); 
    auto legendMass = new TLegend(0.6, 0.6, 0.9, 0.9);
    if (!isSystV0Analysis)    legendMass->SetHeader("Multiplicity classes"); 
    else     legendMass->SetHeader("Selections"); 
    c[j] = new TCanvas (Form("c%i",j),nomi[j], 1800,1000);//1400, 500 per averle in orizzontale oppure 800, 1000 per averle in verticale
    c[j]->Divide(4,2);
    Int_t counterParticleType=0;
    for (Int_t ParticleType=type; ParticleType<=type+1; ParticleType++){
      cout <<"Particle Type " <<  tipo[ParticleType] << endl; 
      counterParticleType++;
     
      for(Int_t dataorMC=0; dataorMC<2; dataorMC++){
	cout << "data or MC " << dataorMC << endl;
	if (dataorMC==0 && j>=num_isto_finale-4) continue;
	hs[j] = new THStack(Form("hs%i",j),"");
	hs_ratio[j] = new THStack(Form("hs_ratio%i",j),"");
	Title[ParticleType-type]=isDataorMC[dataorMC]+"  "+ tipobis[ParticleType];
	Int_t ForLoopLimit=0;
	Int_t ForLoopLimInf=0;
	if (isSystV0Analysis){ ForLoopLimit=numsysV0; ForLoopLimInf=0;}
	if (!isSystV0Analysis) {ForLoopLimit=0;ForLoopLimInf=mult;}
	Int_t 	Icounter=-1;
	for (Int_t i =ForLoopLimInf ; i>= ForLoopLimit; i--){//ok for mult loop
	//for (Int_t i =ForLoopLimInf ; i< ForLoopLimit; i++){//ok for systloop

	  Icounter++;
	  //if multiplicity dependence is analyzed*****************
	  if (!isSystV0Analysis){
	    innamecompleto[i]= Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+ "_" +MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", i, sysTrigger, sysV0, sys, PtTrigMin);
	    if (type>=4) innamecompleto[i]= Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+ "_" +MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f_Rap%i.root", i, sysTrigger, sysV0, sys, PtTrigMin, rap);

	  }
	  //*******************************************************

	  //if particle selection dependence is analyzed*****************
	  if (isSystV0Analysis){
	    innamecompleto[i] = Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+ "_" +MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", chosenmolt, sysTrigger, i, sys, PtTrigMin); 
	    if (type>=4)		  innamecompleto[i] = Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+ "_" +MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f_Rap%i.root", chosenmolt, sysTrigger, i, sys, PtTrigMin, rap); 

	  }
	  cout <<"fileinput completo " <<  innamecompleto[i]<< endl;
	  //*******************************************************
	  myfile[i] = new TFile(innamecompleto[i], ""); 
	  histo[j][i] = (TH1F*)myfile[i]->Get("histo_"+tipo_histo[j]);
	  hs[j]->Add(histo[j][i], "p");
	  hs[j]->Add(histo[j][i]);
	  histo[j][i]->SetMarkerStyle(33);
	  histo[j][i]->SetMarkerSize(1.2);
	  histo[j][i]->SetLineColor(Color[i]);
	  histo[j][i]->SetMarkerColor(Color[i]);

	  if (Icounter==0) 	histo_master[j] = (TH1F*)histo[j][i] ->Clone("histo_"+tipo_histo[j]+"_master");
	  histo_ratio[j] = (TH1F*)	histo[j][i] ->Clone("histo_"+tipo_histo[j]+"_ratio");

	  
	  Bool_t SkipIsto=kFALSE;
	  for (Int_t IPt=1; IPt <histo[j][i]->GetNbinsX(); IPt++){ //parto da IPt = 1 perché nelle Xi il primo bin è vuoto (< 0.3 GeV/c)
	    cout << IPt << " " << SkipIsto << endl;
	    if (histo_master[j]->GetBinContent(IPt+1)==0) SkipIsto=kTRUE;
	  }
	  
	  //cout << "SkipIsto " << SkipIsto << endl;
	 	 
	  if (!SkipIsto)	  histo_ratio[j] ->Divide( histo_master[j] );
	  histo_ratio[j]->SetMarkerStyle(33);
	  histo_ratio[j]->SetMarkerSize(1.2);
	  histo_ratio[j]->SetLineColor(Color[i]);
	  histo_ratio[j]->SetMarkerColor(Color[i]);
	  
	  if (Icounter!=0)	  hs_ratio[j]->Add(histo_ratio[j], "p");
	  if (Icounter!=0)	  hs_ratio[j]->Add(histo_ratio[j]);
	  //	  cout << "ho aggiunto histo ratio" << endl;
	  if (dataorMC==1 && !isSystV0Analysis && counterParticleType==1) legend->AddEntry(histo[j][i],molteplicit[i],"pel");
	  if (dataorMC==1 && !isSystV0Analysis && counterParticleType==1) legendMass->AddEntry(histo[j][i],molteplicit[i],"pel");
	  //if particle selection dependence is analyzed*****************
	  if (dataorMC==1 && isSystV0Analysis && counterParticleType==1) legend->AddEntry(histo[j][i],systematic[i],"pel");
	  if (dataorMC==1 && isSystV0Analysis && counterParticleType==1) legendMass->AddEntry(histo[j][i],systematic[i],"pel");
	  if(i==0 && j==1 && dataorMC==0 && counterParticleType==1) legendMass->AddEntry(retta, tipoNotSign[ParticleType]+" PDG mass");
	  outfile->cd();
	  // histo[j][i]->Write();
	  // histo_ratio[j]->Write();
	  cout << "I'm triying to close the file" << Form(innamedataorMC[dataorMC]+"_" +tipobis[ParticleType]+ "_" +MassFixedPDG[isMeanFixedPDG]+ BkgType[isBkgParab]+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", chosenmolt, sysTrigger, i, sys, PtTrigMin) << endl;
	  //       	myfile[i]->Close();      
	  cout << "file closed" << endl;
	}
      
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
	if (ParticleType==type)      c[j]->cd(dataorMC+1+4);
	if (ParticleType==type+1)      c[j]->cd(dataorMC+3+4);
	//     if(j == 6) {gPad->SetLogy();}
	//	cout << "ok 3"<< endl;
	gPad->SetLeftMargin(0.2);
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
	//	cout << "ok 5"<< endl;
	if (j==1){
	 hs_ratio[j]->SetMinimum(0.9985);
	 hs_ratio[j]->SetMaximum(1.0015);
	}
	//cout << "ok 6"<< endl;
	legend->Draw();
	rettaUno->Draw("same");
	//cout << "ok 7"<< endl;
      } //data or MC
    } //particle type
  }

  /*
    for (Int_t j =0 ; j< num_isto; j++){ //ciclo sugli istogrammi relativi al MC
    if (j <6 && j!=1) {
    nomi_MC_completi[j] = nomi_MC[j]+" (MC)";
    }
    else if (j ==6 || j==1) {
    nomi_MC_completi[j] = nomi_MC[j]+" (MC)";
    }	
    else {
    nomi_MC_completi[j] = nomi_MC[j];
    }  
    hs_MC[j] = new THStack(Form("hs%i_MC",j), nomi_MC_completi[j]);
    auto legend = new TLegend(X1_MC[j],Y1_MC[j],X2_MC[j],Y2_MC[j]);
    legend->SetHeader("Intervalli di molteplicita'"); 
    //c_MC[j] = new TCanvas (Form("c%i_MC",j),nomi_MC[j]+ "_MC",800, 600 );
    for (Int_t i =0 ; i< mult+1; i++){
    myfile[i] = new TFile(Form(innameMC+"_molt%i_sysT%i_sysV0%i_Sys%i_PtMin%.1f.root", i, syst, sysV0, sys, PtTrigMin), ""); 
    histo_MC[j][i] = (TH1F*)myfile[i]->Get("histo_"+tipo_histo[j]);
    histo_MC[j][i]->SetMarkerStyle(marker[i]);
    if (i+1 != 5) {histo_MC[j][i]->SetLineColor(i+1); histo_MC[j][i]->SetMarkerColor(i+1);}
    if (i+1 == 5) {histo_MC[j][i]->SetLineColor(6);   histo_MC[j][i]->SetMarkerColor(6);}
    if (i+1 == 3) {histo_MC[j][i]->SetLineColor(8);   histo_MC[j][i]->SetMarkerColor(8);}
    hs_MC[j]->Add(histo_MC[j][i]);
    legend->AddEntry(histo_MC[j][i],molteplicit[i],"pel");
    outfile->cd();
    histo_MC[j][i]->Write();
    }
   
    if (j<num_isto){
    c[j]->cd(2);
    //      if(j == 6) {gPad->SetLogy();}
    }
    /*
    if (j>=num_isto-1){ 
    c_MC[j-num_isto+1] = new TCanvas (Form("c%i_MC",j-num_isto+1),nomi_MC[j]+ "_MC",800, 600 );
    c_MC[j-num_isto+1]->cd();
    cout << j - num_isto+1<< endl;
    // hs_MC[j]->SetTitle(Form("hs%i_MC",j), nomi_MC[j]+ " (" + tipo2[type] +")");
    }
   
    hs_MC[j]->Draw("nostack"); //cosi' me li disegna in maniera piu' corretta
    hs_MC[j]->GetXaxis()->SetTitle(titleX); //il set degli assi va fatto dopo l'opzione Draw!!
    hs_MC[j]->GetXaxis()->SetTitleOffset(1.2);
    hs_MC[j]->GetYaxis()->SetTitle(titleY_MC[j]);
    hs_MC[j]->GetYaxis()->SetTitleOffset(1.5); 
    hs_MC[j]->SetMinimum(min_range_MC[type][j]);
    hs_MC[j]->SetMaximum(max_range_MC[type][j]);
    legend->Draw();
    }
  */
  for (Int_t j =0 ; j< num_isto_finale; j++){
    //    c[j]->SaveAs("FinalOutput/DATA2016/IstoConfronto"+year+"_"+ tipo_histo[j]+ def[cut]+"_"+tipo[type]+Form("_PtMin%.1f.pdf", PtTrigMin));
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
  cout << "\n\n\n il file di output si chiama " << outname << endl;



}