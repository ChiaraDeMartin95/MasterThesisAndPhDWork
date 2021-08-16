#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLegend.h>
#include <Macros/ErrRatioCorr.C>
#include <Macros/BarlowVariable.C>

Double_t SetEfficiencyError(Int_t k, Int_t n){
  return sqrt(((Double_t)k+1)*((Double_t)k+2)/(n+2)/(n+3) - pow((Double_t)(k+1),2)/pow(n+2,2));
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString titleX, TString titleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset){
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(titleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(titleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); //1.2                                                                                       
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void MacroRatioHistos(Int_t RunVar=4, Int_t multDef = 5 /*only used when a given multiplicity has to be manually selected*/, Int_t type =0, Int_t AnalysisType =0, TString  CommonFileName = ""/*name common to all files */,  TString OutputName ="", Bool_t isBarlow=1, Bool_t ispp5TeV=0){

  //RunVar should be increased when you want to do a new comparison; also, OutputName And CommonFileName should be updated below!!
  TString TypeAnalysis[3] = {"Jet", "Bulk", "All"};
  TString TypeAnalysisBis[3] = {"Jet", "Bulk", "Inclusive"};
  TString TypeAnalysisTer[3] = {"near-side jet", "out of jet", "full"};
  TString tipo[2] = {"K0s", "Xi"};
  TString SCorrelation[3] = {"No corr.", "Full", "Partial"};
  Int_t CorrelationBtwHistos=0;
  Int_t numFiles=0;
  Int_t TypeOfComparison = 1; 

  Int_t sys=0;
  Int_t sysPhi=0;
  Int_t NSign=3;
  Int_t NSigma=2;
  //TypeOfComparison:
  //1) Comparison with respect to one plot (e.g., with respect to plot in file number 0, which might be a specific mult class or a specific pt,trig min). The number of the file with respect to which the ratio is computed is defined below
  //2) Comparison between a variable calculated in one way and the same variable calculated in another way -- done in CLASSES (e.g. multiplicity classes, pt,trig min classes)
  //3) Like 1, but for two different selections (e.g. I compare spectra calculated in one way (way A) to 0-100% calculated in way A, and spectra calculated in another way (way B) to 0-100% (way B) )

  if (RunVar==0){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelationRun2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3";
    numFiles=8;
  }
  else if (RunVar==1){
    CommonFileName= "FinalOutput/DATA2016/histo/AngularCorrelation2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName="hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOB";
    numFiles=4;
  }
  else if (RunVar==2){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOB_Comparison";
    numFiles=8;
    TypeOfComparison = 2; 
  }
  else if (RunVar==3){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius_Comparison";
    numFiles=12;
    TypeOfComparison = 2; 
  }
  else if (RunVar==4){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_SysT0_SysV00_Sys0_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius";
    numFiles=6;
  }
  else if (RunVar==5){
    CommonFileName = "FinalOutput/DATA2016/histo/AngularCorrelation2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_TOFOOBRadius_SkipAllAssocComp";
    numFiles=8;
    TypeOfComparison = 2; 
  }
  else if (RunVar==6){
    //comparison between Xi spectra obtained with default selections, only-TOF OOB pileup rejection and +RV0 < 34 cm
    CommonFileName = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_";
    OutputName = "hXiEstimateTriggerRun3/hXiEstimateForTriggerRun3_SpectraComparison";
    numFiles=3;
  }
  else if (RunVar==7){ 
    //*****TYPE 1: comparison bertween efficiency calculated in different ways (same multiplicity class)*****
    //comparison between Xi efficiency obtained with old AODs and with new (refiltered) AODs (AOD235)
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = "FinalOutput/DATA2016/Efficiency/EfficiencyComparisonOldNewAODsXi";
    numFiles=12; //2 x 6 mult classes
    TypeOfComparison = 2; 
  }
  else if (RunVar==8){
    //comparison norm factor for K0s obtained with 200M INT7 events to the one obtained with 300M INT7 events
    CommonFileName = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid";
    OutputName = "FinalOutput/DATA2016/MCClosureCheckComparisonNormFactor";
    numFiles=12; //2 x 6 mult classes
    TypeOfComparison = 2; 
  }
  else if (RunVar==9){
    //multiplicity comparison norm factor for K0s obtained with 200M INT7 events to the one obtained with 300M INT7 events
    CommonFileName = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid";
    OutputName = "FinalOutput/DATA2016/MCClosureCheckMultComparisonNormFactor";
    numFiles=6; //6 mult classes
  }
  else if (RunVar==10){
    //comparison norm factor for K0s obtained with 200M INT7 events to the one obtained with 300M INT7 events
    CommonFileName = "CompareYieldDifferentCollisions";
    OutputName = "CompareYieldHMK0s_All_16k_vs_AllBut16k";
    numFiles=2; //
  }
  else if (RunVar==11){
    //comparison between K0s efficiency from 18g4_extra to K0s efficiency from 1617MC
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = "CompareK0sEfficiency_18g4_extra_vs_1617MC";
    numFiles=12; //2+6 mult
    TypeOfComparison = 2; 
  }
  else if (RunVar==12){
    //comparison MB yield obtained with 18g4_extra eff (pt) to those obtained with 18g4_extra eff (pt, eta)
    CommonFileName = "CompareYieldDifferentCollisions";
    OutputName = "CompareYieldMB_EtaEff_vs_NoEtaEff_18g4extra_"+ TypeAnalysis[AnalysisType];
    numFiles=2; //
  }
  else if (RunVar==13){
    //comparison HM yields vs mult obtained from all but 16k data; comparison is between results corrected by eff 2019h11c_extra and those corrected by efficiency from 2019h11a,b,c_extra
    CommonFileName = "FinalOutput/DATA2016/SystematicAnalysisAllhK0sHM_RedNo16k_PtBinning1_K0s_Eta0.8_"+  TypeAnalysis[AnalysisType]+"Data_PtMin3.0_IsEtaEff";
    OutputName = "CompareYieldHMK0s_AllBut16k_Eff2019h11c_vs_2019h11abc_"+ TypeAnalysis[AnalysisType];
    numFiles=2; //
  }
  else if (RunVar==14){
    //comparison between 5 TeV and 13 TeV K0s efficiency
    CommonFileName ="FinalOutput/DATA2016/Efficiency/Efficiency"; 
    OutputName = "CompareK0sEfficiencyAcrossEnergies";
    numFiles=12; //2*6 mult
    TypeOfComparison = 2; 
  }
  else if (RunVar==15){
    //comparison K0s or Xi at 13 TeV vs K0s at 5 TeV (yields)
    CommonFileName = "FinalOutput/DATA2016/";
    //    OutputName = "CompareYieldK0s_5TeV_vs_13TeV_"+ TypeAnalysis[AnalysisType];
    OutputName = "CompareYield" + tipo[type] + "_5TeV_vs_13TeVNewAODs_AllNormCorr"+ TypeAnalysis[AnalysisType];
    numFiles=2; //
  }
  else if (RunVar==16){
    //comparison K0s at 13 TeV vs K0s at 5 TeV (spectra)
    CommonFileName = "FinalOutput/DATA2016/";
    //oldAODs 13 TeV    OutputName = "CompareSpectraK0s_5TeV_vs_13TeV_"+ TypeAnalysis[AnalysisType];
    OutputName = "CompareSpectraK0s_5TeV_vs_NewAODs13TeV_"+ TypeAnalysis[AnalysisType];
    numFiles=12; //
    TypeOfComparison=3;
  }
  else if (RunVar==17){
    //comparison K0s purity at 13 TeV vs K0s at 5 TeV
    CommonFileName = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_";
    OutputName = "ComparePurityK0s_5TeV_vs_13TeV";
    numFiles=12; //
    TypeOfComparison=3;
  }
  else if (RunVar==18){
    //comparison K0s efficiency at 5 TeV CENTwoSDD vs FAST
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hK0s_";
    OutputName = "CompareFASTCENTwoSDDeff_K0s_5TeV_vs_13TeV";
    numFiles=12; //
    TypeOfComparison = 2; 
  }
  else if (RunVar==19){
    //COMPARE EFFICIENCIES ACROSS MULTIPLICITIES
    //    CommonFileName ="FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_MCEff_Efficiency_SysT0_SysV00_PtMin3.0.root";
    //    OutputName = "CompareEfficiencyAcrossMult_2018f1_extra_hK0s_PtMin3_MasterThesis";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_hK0s_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //    OutputName = "CompareEfficiencyAcrossMult_2018f1_extra_hK0s_y0.5_AllAssoc_PtMin0.15";

    //CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency1617_GP_AOD235_With18c12b_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    //OutputName = "CompareEfficiencyAcrossMult_1617_GP_AOD235_extra_hK0s_Eta0.8_PtMin3";

    //CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_AOD235_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    //OutputName = "CompareEfficiencyAcrossMult_2018f1_extra_AOD235_hK0s_Eta0.8_PtMin3";

    //CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hK0s_pttrig0.15_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //OutputName = "CompareEfficiencyAcrossMult_LHC17pq_pp5TeV_hK0s_y0.5_PtMin0.15";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    //OutputName = "CompareEfficiencyAcrossMult_LHC17pq_pp5TeV_hK0s_PtMin3.0";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/EfficiencyAllMC_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    //    OutputName = "CompareEfficiencyAcrossMult_AllMC_hXi_PtMin3.0";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2016kl_hXi_FioPtBins_Xi_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //OutputName = "CompareEfficiencyAcrossMult_AllMC_hXi_y0.5_AllAssoc_PtMin0.15";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_15runs_NoTrackLengthSelectionFull_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //    OutputName = "CompareEfficiencyAcrossMult_NoTrackLengthSelection_hK0s_y0.5_AllAssoc_PtMin0.15";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_15runs_NoOOBPileUp_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //    OutputName = "CompareEfficiencyAcrossMult_NoOOBPileUp_hK0s_y0.5_AllAssoc_PtMin0.15";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_15runs_NSigma5_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //    OutputName = "CompareEfficiencyAcrossMult_NSigma5_hK0s_y0.5_AllAssoc_PtMin0.15";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_15runs_V0RadiusSelection_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
    //OutputName = "CompareEfficiencyAcrossMult_V0RadiusSelection_hK0s_y0.5_AllAssoc_PtMin0.15";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_15runs_Trig0Pt_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.0.root";
    //    OutputName = "CompareEfficiencyAcrossMult_Trig0Pt_hK0s_y0.5_AllAssoc_PtMin0.0";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency15g3c3_hK0s_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2_LooseCosinePAngle.root";
    //    OutputName = "CompareEfficiencyAcrossMult_15g3c3_hK0s_y0.5_AllAssoc_PtMin0.15_LooseCosinePAngle";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency17pq_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0_MultBinning3.root";
    //    OutputName = "CompareEfficiencyAcrossMult_17pq_hXi_MultBinning3";

    //    CommonFileName = "FinalOutput/DATA2016/Efficiency/EfficiencyLHC18_GP_AOD235_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    //    OutputName = "CompareEfficiencyAcrossMult_LHC18_AOD235_hXi";

    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency161718_AOD235_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    OutputName = "CompareEfficiencyAcrossMult_LHC161718_AOD235_hXi";

    numFiles=6; 
  }
  else if (RunVar==20){
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_";
    OutputName = "CompareEfficiencyAcrossMult_And_OldNewAODs";
    numFiles = 12;
    TypeOfComparison=3;
  }
  else if (RunVar==21){
    //Comparison of efficiencies vs pt
    if (type==0)    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_";
    else if (type==1) CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    //OutputName = "CompareEfficiency_OldNewAODs_" +tipo[type];
    OutputName = "CompareEfficiency_TwoDiffDatasets_" +tipo[type];
    numFiles = 12;
    TypeOfComparison=2;
  }
  else if (RunVar==22 || RunVar==23){
    //Comparison of efficiencies vs eta
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency2018f1_extra_";
    if (RunVar==22)    OutputName = "CompareEfficiencyEtaAcrossMult_And_OldNewAODs";
    if (RunVar==23)    OutputName = "CompareEfficiencyEta_OldNewAODs";
    numFiles = 12;
    if (RunVar==22)    TypeOfComparison=3;
    else  TypeOfComparison=2;
  }
  else if (RunVar==24 || RunVar==25){
    //RunVar==24: compare 2018f1_extra K0s eff with 2017e5_extra K0s eff; comparison done for both standard efficiency and efficiency after mass selection
    //RunVar==25: compare K0s afficiency obtained with and without mass selection; comparison done for both 2018f1_extra and 2017e5_extra
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    if (RunVar==24)    OutputName = "CompareEfficiencyAcrossPeriod_And_wAndwoMassSel_m50100";
    if (RunVar==25)    OutputName = "CompareEfficiency_wAndwoMassSel_m50100";
    numFiles = 4;
    if (RunVar==24)    TypeOfComparison=3;
    else  TypeOfComparison=2;
  }
  else if (RunVar==26 || RunVar==27){
    //Compare data to MC spectra in jet, out-of-jet and full; K0s (26) and Xi (27)
    CommonFileName = "";
    if (RunVar==26)    OutputName = "CompareDataMCK0s_" + TypeAnalysis[AnalysisType];
    else if (RunVar==27)    OutputName = "CompareDataMCXi_"  + TypeAnalysis[AnalysisType];
    if (RunVar==26) type=0;
    else if (RunVar==27) type=1;
    numFiles = 2*6; //2* mult
    TypeOfComparison=2;
  }
  else if (RunVar==28 || RunVar==29){
    //Compare spectra in jet, out-of-jet and full to 0-100%; done both for data and MC; K0s (28) and Xi (29)
    CommonFileName = "";
    if (RunVar==28)    OutputName = "ComparePtSpectraMultK0s_" + TypeAnalysis[AnalysisType];
    else if (RunVar==29)    OutputName = "ComparePtSpectraMultXi_"  + TypeAnalysis[AnalysisType];
    if (RunVar==28) type=0;
    else if (RunVar==29) type=1;
    numFiles = 2*6; //2* mult
    TypeOfComparison=3;
  }
  else if (RunVar==30){
    //Compare data vs MC yields vs mult in jet, out-of-jet and full to 0-100%; 
    CommonFileName = "";
    OutputName = "CompareYieldvsMultDatavsMC_" + TypeAnalysis[AnalysisType] + "_" + tipo[type];
    numFiles = 2; 
    TypeOfComparison=1;
  }
  else if (RunVar==31){
    //Compare K0s spectra of different periods (corrected by corresponding efficiency)
    CommonFileName = "FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_MCEff_PtBinning1";
    OutputName = "CompareK0sYieldDiffPeriods";
    numFiles = 5; 
    TypeOfComparison=1;
  }
  else if (RunVar==32){
    //Compare K0s efficiency (1617_GP) with and without period 18c12
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency1617_GP_AOD235";
    OutputName = "CompareK0sEfficiency1617GP";
    numFiles = 2*6; 
    TypeOfComparison=2;
  }
  else if (RunVar==33){
    //Compare old AOD spectra to new AOD spectra (no correction by etaEff and norm factor) in jet, out-of-jet and full; 
    CommonFileName = "";
    OutputName = "CompareSpectraOldNewAODs_" + tipo[type]+"_"  + TypeAnalysis[AnalysisType];
    numFiles = 2*6; //2* mult
    TypeOfComparison=2;
  }
  else if (RunVar==34){
    //Compare new AOD spectra with and without eta dependent efficiency in jet, out-of-jet and full; 
    CommonFileName = "";
    OutputName = "CompareSpectraWithWoEtaEff_" + tipo[type]+"_"  + TypeAnalysis[AnalysisType];
    numFiles = 2*6; //2* mult
    TypeOfComparison=2;
  }
  else if (RunVar==35){
    //Compare K0s efficiency vs multiplicity for INT7 and HM events
    CommonFileName = "";
    OutputName = "CompareSpectraWithWoEtaEff_" + tipo[type]+"_"  + TypeAnalysis[AnalysisType];
    numFiles = 2*6; //2* mult
    TypeOfComparison=2;
  }
  else if (RunVar==36){
    //compare K0s raw yields obtained with loosest and tightest selections with default raw yields. Study performed for systematic uncertainty determinaton.
    CommonFileName = " FinalOutput/DATA2016/invmass_distribution_thesis/invmass_distribution_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV0";
    OutputName = "CompareSSBK0sDifferentSelections_1617_AOD234";
    numFiles =3;
    TypeOfComparison = 1;
  }
  else if (RunVar==37){
    //comparison K0s efficiency 2018f1_extra and 2015g3c3 
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = "CompareSpectraK0sEfficiency_2018f1_2015g3c3_vsMultiplicity";
    numFiles=12; //
    TypeOfComparison=3;
  }
  else if (RunVar==38){
    //comparison HM K0s efficiency AOD230 vs new AODs
    CommonFileName = "";
    OutputName = "CompareHMK0sEfficiency_AOD230_vs_newAODs";
    numFiles=6; 
    TypeOfComparison=2;
  }
  else if (RunVar==39){
    //comparison of two methods to subtract OOJ (K0s) 
    CommonFileName = "FinalOutput/DATA2016/SystematicAnalysis1617_hK0s_PtBinning1_";
    OutputName = "CompareOOJMethodsK0s_Default_vs_NewOOJScaledto1dphi2";
    //OutputName = "CompareOOJMethodsK0s_Default_vs_NewOOJScaledto1dphi2_MultCorr";
    //OutputName = "CompareOOJMethodsK0s_Default_vs_NewOOJScaledtoDefBulk";
    //OutputName = "CompareOOJMethodsK0s_Default_vs_NewOOJScaledtoDefBulk_MultCorr";
    numFiles=12; 
    TypeOfComparison=2;
  }
  else if (RunVar==40){
    //ratio between Xi and K0s yields
    CommonFileName = "FinalOutput/DATA2016/";
    OutputName = "CompareXiK0sYields_13TeVvs5TeV";
    numFiles=12; 
    TypeOfComparison=2;
  }
  else if (RunVar==41){
    //ratio between Xi efficiency of 2016, 2017 and 2018
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = Form("CompareXiEfficiencyOfDifferentYears_m%i", multDef);
    numFiles=3; //16, 17, 18. multiplicity should be chosen 
    TypeOfComparison=1;
  }
  else if (RunVar==42){
    //ratio between eff (eta) for low pt, trig (0.15 < pt < 2.5) and the one for pt,trig > 3
    CommonFileName = "FinalOutput/DATA2016/Efficiency/Efficiency";
    OutputName = "Compare" + tipo[type] + "EfficiencyDiffPtTrig";
    numFiles=12; 
    TypeOfComparison=2;
  }
  else if (RunVar==43){
    //Comparison of norm factor obtained using different DeltaPhi and DeltaEta choices
    CommonFileName = "FinalOutput/DATA2016/MCClosureCheck_HybridtoTrue1617GP_hK0s_Hybrid_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0";
    OutputName = "FinalOutput/DATA2016/SystNormFactor" + tipo[type] + "_"+  TypeAnalysis[AnalysisType];
    numFiles=12; //4 x DeltaEta choice x 3 x DeltaPhi choice 
    TypeOfComparison=1;
    NSign = 4;
    NSigma=2;
  }

  //this macros superimpose in a same canvas pad same histos found in different files. A loop over the files is done. The ratio of the histos to the histo found in the first file is performed. You can choose if histos have to be considered fully correlated, uncorrelated, or if the histos are obtaiend from a subsample of the data used to obtain the histo found in the first file (i.e. partial correlation)

  TString VarName[numFiles] = {""}; 
  TString NameHisto[numFiles] = {""}; 
  TString histoName= "fHistQA6";
  if (RunVar==6) histoName = "histo_S";
  else   if (RunVar==7 || RunVar==11 || RunVar==14 || RunVar==18 || RunVar==19 || RunVar==20 || RunVar==21 || RunVar==24 || RunVar==25 || RunVar==32 || RunVar==37 || RunVar==38 || RunVar==41) histoName = "fHistV0EfficiencyPtBins";
  else if (RunVar==22 || RunVar==23 || RunVar==42) histoName = "fHistV0EfficiencyEta";
  else   if (RunVar==8 || RunVar==9 || RunVar==43) histoName = "SpectrumRatio" +  TypeAnalysis[AnalysisType]+"_m"; //All -> Jet, Bulk
  else if (RunVar==10 || RunVar==12) histoName = "histoYieldComparison";
  else if (RunVar==13) histoName = "fHistYieldvsErrSoloStat";
  else if (RunVar==15 || RunVar==16) histoName = "";
  else if (RunVar==17) histoName = "histo_SSB";
  else if (RunVar>=26 && RunVar<=29) histoName = "";
  else if (RunVar==30) histoName = "";
  else if (RunVar==31) histoName = "histo_SEffCorr";
  else if (RunVar==33 || RunVar==34) histoName = "";
  else if (RunVar==36) histoName = "histo_SSB";
  else if (RunVar==39) histoName = "fHistSpectrumPart_";
  else if (RunVar==40) histoName = "";

  Int_t numDef=3;
  if (RunVar>=7) numDef=0;
  Float_t num=0;
  Int_t numEff=0;

  Float_t LimSupPol0=8;
  Float_t LimInfPol0 = 0;
  if (RunVar==10 || RunVar==13) LimSupPol0 = 45;
  else if (RunVar==15 || RunVar==30) LimSupPol0 = 30;
  if (RunVar==13) LimInfPol0 = 25;
  if (RunVar==22 || RunVar==23) {LimInfPol0 = -0.8; LimSupPol0 = 0.8;}
  TF1 * pol0At1 = new TF1("pol0", "pol0", LimInfPol0, LimSupPol0);
  pol0At1->SetParameter(0,1);
  pol0At1->SetLineColor(1);
  pol0At1->SetLineWidth(0.5);
  TF1 * pol1[numFiles]; 
  TF1* 	lineatB0 = new TF1("pol0B", "pol0", 0, 2.5);
  lineatB0->SetParameter(0,0);
  lineatB0->SetLineColor(1);

  //histo style selections
  Float_t Low=10e-8;
  Float_t Up=0.003;
  Float_t LowRatio=10e-8;
  Float_t UpRatio=0.6;
  Float_t LowBarlow=0.5;
  Float_t UpBarlow=1.5;
  Float_t LimSupSys=1;
  Float_t UpError=0.1;
  //  Int_t color[10]={1,402 , 628, /*905, 881,*/601, 867, 418, 905, 881};
  Int_t color[13]={1, 401,801,628, 909, 881,860,868,841,  418, 881, 7, 1};
  Int_t style =27;
  TString titleX = "Multiplicity class";
  TString titleY="#hXi events / #INT7 events";
  TString title="Fraction of INT7 events containing trigger particle and Xi";
  TString titleRatio="Ratio to #it{p}_{T}^{trigg} > 3 GeV/#it{c}";
  TString Smolt[6] = {"_0-5", "_5-10", "_10-30", "_30-50", "_50-100", "__all"};  
  TString Smolt5TeV[6] = {"_0-10", "_10-100", "_100-100", "_100-100", "_100-100", "__all"};  
  TString SmoltHM[6] = {"_0-0.001", "_0.001-0.005", "_0.005-0.01", "_0.01-0.05", "_0.05-0.1", "__all"};  
  TString SmoltBis[6] = {" 0-5%", " 5-10%", " 10-30%", " 30-50%", " 50-100%", " 0-100%"};  
  TString SmoltBis5TeV[6] = {" 0-10%", " 10-100%", " 100-100%", " 100-100%", " 100-100%", " 0-100%"};  

  if (ispp5TeV){
    for (Int_t i=0; i<6; i++) { Smolt[i] = Smolt5TeV[i];  SmoltBis[i] = SmoltBis5TeV[i];}
  }

  if (RunVar==2){
    titleRatio = "Ratio of only TOF to (SPD+TOF) pileup rej";
    LowRatio = 0.5;
    UpRatio = 1;
  }

  else  if (RunVar==3){
    titleRatio = "Ratio of Run3 to Run 2 selections";
    LowRatio = 0.5;
    UpRatio = 1;
  }

  else  if (RunVar==5){
    titleRatio = "Ratio of SkipAssoc to AllAssoc with Run3Sel";
    LowRatio = 0.8;
    UpRatio = 1.3;
  }
  else  if (RunVar==6){
    titleRatio = "Ratio to default Run2Sel";
    Low = 10e-6;
    Up = 5*10e-3;
    LowRatio = 0;
    UpRatio = 1.5;
    titleY = "1/N_{trigg} dN^{#Xi}/dp_{T}";
    titleX = "p_{T} (GeV/c)";
    title = "Raw Xi spectrum 0-100%, p_{T}^{trigg} > 3 GeV/c";
  }
  else  if (RunVar==7){
    titleRatio = "Ratio to new AODs";
    Low = 0;
    Up = 0.1;//0.3
    LowRatio = 0.8;
    UpRatio = 1;
    titleY = "Xi efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "Xi efficiency";
  }
  else  if (RunVar==8){
    titleRatio = "Ratio to 300M events (new) norm factor";
    Low = 0.9;
    Up = 1.2;
    LowRatio = 0.95;
    UpRatio = 1.1;
    titleY = "K0s norm factor";
    titleX = "p_{T} (GeV/c)";
    title = "K0s norm factor";
  }
  else  if (RunVar==9){
    titleRatio = "Ratio to 0-100% norm factor";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s norm factor";
    titleX = "p_{T} (GeV/c)";
    title = "K0s norm factor";
  }
  else  if (RunVar==10){
    titleRatio = "Ratio to 16k only";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==11){
    titleRatio = "Ratio to 1617MC";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s efficiency";
    titleX = "p_{T}";
    title = "K0s efficiency";
  }
  else  if (RunVar==12){
    titleRatio = "Ratio to non-eta-dependent eff";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==13){
    titleRatio = "Ratio to 2019h11c_extra eff";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s yield";
    titleX = "dN/d#eta";
    title = "K0s yield";
  }
  else  if (RunVar==14){
    titleRatio = "Ratio to 13 TeV efficiency";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "K0s efficiency";
    titleX = "p_{T}";
    title = "K0s efficiency";
  }
  else  if (RunVar==15){
    titleRatio = "Ratio to 13 TeV";
    if (AnalysisType==0){
      Low = 0.02;
      Up = 0.04;
      if (type==1){
	Up = 0.004;
      }
    }
    else {
      Low = 0;
      Up = 0.3;
      if (type==1){
	Up = 0.02;
      }
    }
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = tipo[type]+" yield";
    titleX = "dN/d#eta";
    title = tipo[type] + " yield";
  }
  else  if (RunVar==16){
    titleRatio = "Ratio to 0-100%";
    if (AnalysisType==0){
      Low = 0.0005;
      Up = 0.02;
    }
    else {
      Low = 0.00005;
      Up = 0.3;
    }
    LowRatio = 0;
    UpRatio = 2;
    titleY = "1/N_{trigg} dN/dp_{T} (1/(GeV/c))";
    titleX = "p_{T} (GeV/c)";
    title = "K0s yield";
  }
  else  if (RunVar==17){
    titleRatio = "Ratio to 0-100%";
    Low = 0.95;
    Up = 1;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "Purity";
    titleX = "p_{T} (GeV/c)";
    title = "K0s purity";
  }
  else  if (RunVar==18){
    titleRatio = "Ratio to CENTwoSDD";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.95;
    UpRatio = 1.05;
    titleY = "Efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "K0s efficiency";
  }
  else  if (RunVar==19 || RunVar==20 || RunVar==21 || RunVar==24 || RunVar==25 || RunVar==37){
    titleRatio = "Ratio to 0-100%";
    if (RunVar==21)     titleRatio = "Ratio to old AODs";
    else if (RunVar==24)     titleRatio = "Ratio to 2018f1_extra";
    else if (RunVar==25)     titleRatio = "Ratio to standard efficiency";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.85;
    UpRatio = 1.15;
    titleY = tipo[type] +" efficiency";
    titleX = "p_{T} (GeV/c)";
    title = tipo[type] + " efficiency";
  }
  else  if (RunVar==22 || RunVar==23){
    titleRatio = "Ratio to 0-100%";
    if (RunVar==23)     titleRatio = "Ratio to old AODs";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.8;
    UpRatio = 1.2;
    titleY = "K0s efficiency";
    titleX = "#eta";
    title = "K0s efficiency";
  }
  else  if (RunVar>=26 && RunVar<=29){
    if (RunVar==26 || RunVar ==27)    titleRatio = "Ratio MC/DATA";
    else titleRatio = "Ratio to 0-100%";
    if (AnalysisType==0){
      Low = 0.0005;
      Up = 0.04;
      if (type==1){ //Xi
	Low = 10-7;
	Up = 0.008;
      }
    }
    else {
      Low = 0.00005;
      Up = 0.3;
      if (type==1){ //Xi
	Low = 10-6;
	Up = 0.01;
      }
    }
    LowRatio = 0;
    UpRatio = 2;
    titleY = "1/N_{trigg} dN/dp_{T} (1/(GeV/c))";
    titleX = "p_{T} (GeV/c)";
    title = tipo[type] + " yield "+ TypeAnalysisTer[AnalysisType];
  }
  else  if (RunVar==30){
    titleRatio = "Ratio MC/DATA";
    Low = 0;
    Up = 0.3;
    if (AnalysisType==0) Up = 0.1;
    if (type==1){
      if (AnalysisType==0) Up = 0.004;
      else Up = 0.03;
    }
    LowRatio = 0.4;
    UpRatio = 1.6;
    titleY = "N/N_{trigg} 1/#Delta#eta #Delta#varphi";
    titleX = "#LTdN_{ch}/d#eta#GT_{|#eta|<0.5}";
    title = tipo[type] + " yield "+ TypeAnalysisTer[AnalysisType];
  }
  else  if (RunVar==31){
    titleRatio = "Ratio to 16k";
    Low = 0.0005;
    Up = 1;
    LowRatio = 0;
    UpRatio = 2;
    titleY = "1/N_{evt} dN/dp_{T} (1/(GeV/c))";
    titleX = "p_{T} (GeV/c)";
    title = "K0s inclusive yield ";
  }
  else  if (RunVar==32){
    titleRatio = "Ratio to 1617_GP_AOD235";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.9;
    UpRatio = 1.1;
    titleY = "K0s efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "K0s efficiency";
  }
  else  if (RunVar==33 || RunVar==34){
    if (RunVar==33)    titleRatio = "Ratio to old AODs";
    else    titleRatio = "Ratio to non-eta dependent eff";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.9;
    UpRatio = 1.1;
    titleY = tipo[type] + "1/N_{trigg} dN/dp_{T} 1/#Delta#eta#Delta#varphi";
    titleX = "p_{T} (GeV/c)";
    title = tipo[type] + " yield " + TypeAnalysis[AnalysisType];

    if (AnalysisType==0){
      Low = 0.0005;
      Up = 0.04;
      if (type==1){ //Xi
	Low = 10-7;
	Up = 0.008;
      }
    }
    else {
      Low = 0.00005;
      Up = 0.3;
      if (type==1){ //Xi
	Low = 10-6;
	Up = 0.01;
      }
    }
  }
  else  if (RunVar==36){
    titleRatio = "Ratio to default selections";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.9;
    UpRatio = 1.1;
    titleY = "K0s raw yield";
    titleX = "p_{T} (GeV/c)";
    title = "K0s raw yield";
  }
  else  if (RunVar==38){
    titleRatio = "Ratio to old AODs";
    Low = 0;
    Up = 0.5;
    LowRatio = 0.85;
    UpRatio = 1.15;
    titleY = "K0s efficiency";
    //titleY = "Xi efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "K0s efficiency";
    //title = "Xi efficiency";
  }
  else  if (RunVar==39){
    titleRatio = "Ratio to default subtraction";
    Low = 0;
    Up = 0.02;
    LowRatio = 0.85;
    UpRatio = 1.15;
    LowBarlow=-7;
    UpBarlow=7;
    UpError=0.1;
    titleY = "dN/dp_{T}";
    //titleY = "Xi efficiency";
    titleX = "p_{T} (GeV/c)";
    title = "K0s yield";
    //title = "Xi efficiency";
  }
  else  if (RunVar==40){
    titleRatio = "Xi/K0s yield";
    Low = 0;
    Up = 0.02;
    LowRatio = 0.02;
    UpRatio = 0.13;
    titleY = "";
    titleX = "dN/d{#eta}";
    title = "yield";
  }
  else  if (RunVar==41){
    titleRatio = "Xi efficiency";
    Low = 0;
    Up = 0.3;
    LowRatio = 0.9;
    UpRatio = 1.1;
    titleY = "";
    titleX = "p_{T}";
    title = "efficiency";
  }
  else  if (RunVar==42){
    titleRatio = "Ratio to efficiency in events with 0.15 < pt,trig < 2.5 GeV/c";
    Low = 0;
    Up = 0.3;
    LowRatio = 0;
    UpRatio = 20;
    titleY = "";
    titleX = "#eta";
    title = "Efficiency vs #eta";
  }
  else  if (RunVar==43){
    titleRatio = "Ratio to default norm factor";
    Low = 0.95;
    Up = 1.05;
    LowRatio = 0.95;
    UpRatio = 1.05;
    LowBarlow=-7;
    UpBarlow=7;
    LimSupSys = 0.005;
    titleY = "K0s norm factor";
    titleX = "p_{T} (GeV/c)";
    title = "K0s norm factor";
  }

  TLegend * legend = new TLegend (0.6, 0.7, 0.9, 0.9);
  if (RunVar<6) legend->SetHeader("#it{p}_{T}^{trigg} > ");
  TString LegendName[numFiles]={""};

  TH1F * histo[numFiles];
  TH1F * hdummy[numFiles];
  TH1F * histoRatio[numFiles];
  TH1F * histoBarlow[numFiles];
  Bool_t IsBarlowSign=1;
  TH1F * histoSysError[numFiles];

  TString InputName="";
  TFile * InputFile;
  TString OutputNameRoot= OutputName +".root";
  TString OutputNamepdf= OutputName +".pdf";
  TFile * OutputFile= new TFile (OutputNameRoot, "RECREATE");

  gStyle->SetOptStat(0);
  TCanvas * canvas = new TCanvas("canvas", "canvas", 1300, 800);
  canvas->Divide(2,1);
  TCanvas * canvasB = new TCanvas("canvasB", "canvasB", 1300, 800);
  canvasB->Divide(2,1);

  cout << "TypeOfComparison: " <<     TypeOfComparison<< endl;

  for (Int_t i=0; i<numFiles; i++){
    cout << "\n\e[35mLoop n. " << i << " \e[39m\nfile name: " << InputName << endl;
    NameHisto[i] = "";
    num = numDef+i;
    numEff= num;

    if (RunVar==9 || RunVar==19) numEff = numFiles-num-1; 

    if (TypeOfComparison==2){
      if (i<numFiles/2)      numEff = num;
      else if (i>=numFiles/2) numEff=num-numFiles/2;
    }
    else if (TypeOfComparison==3){
      if (i<numFiles/2)      numEff = numFiles/2-num-1;
      else if (i>=numFiles/2) numEff=numFiles-num-1;
    }

    if (RunVar==0 || RunVar==1 || RunVar==4)    VarName[i] = Form("PtMin%.1f", num);
    else if (RunVar==2){
      if (i>=0 && i<numFiles/2) VarName[i] = Form("Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
      else if (i>=numFiles/2) VarName[i] =  Form("2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==3){
      if (i>=0 && i<numFiles/2) VarName[i] = Form("Run2DataRed_MECorr_hXi_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
      else if (i>=numFiles/2) VarName[i] =  Form("2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==5){
      if (i<=numFiles/2) VarName[i]+= "_AllAssoc";
      VarName[i] += Form("_SysT0_SysV00_Sys0_PtMin%.1f", (float)numEff);
    }
    else if (RunVar==6){
      if (i==0) VarName[i] = "Run2DataRed_MECorr_hXi_Xi_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      if (i==1) VarName[i] = "2016k_pass2_TOFOOBPileUp_Xi_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      if (i==2) VarName[i] = "2016k_TOFOOBPileUp_XiV0Rad34_AOD234_Try2_Xi_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
    }
    else if (RunVar==7){
      if (i>=0 && i<numFiles/2) VarName[i] = "161718_MD_EtaEff_LowPtTrig_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      //"161718_MD_EtaEff_PtTrig3_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      //"161718_MD_EtaEff_LowPtTrig_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      else if (i>=numFiles/2) VarName[i] = "AllMC_hXi_PtTrigMax2.5_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      //"AllMC_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      //"AllMC_hXi_PtTrigMax2.5_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==8){
      if (i>=0 && i<numFiles/2) VarName[i] = "_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
      else if (i>=numFiles/2) VarName[i] = "_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
      NameHisto[i] = Form("%i",numEff);
    }
    else if (RunVar==9){
      VarName[i] = "_New_vs_1617GP_hK0s_PtBinning1_K0s_Eta0.8_PtMin3.0.root";
      NameHisto[i] = Form("%i",numEff);
    }
    else if (RunVar==10){
      if (i==0)   VarName[i] = "_HM16kAll.root";
      else    VarName[i] = "All.root";
      NameHisto[i] = "";
    }
    else if (RunVar==11){
      if (i>=0 && i<numFiles/2) VarName[i] = "1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      else if (i>=numFiles/2) VarName[i] = "2018g4_extra_EtaEff_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==12){
      if (i==0)   VarName[i] = "_MBEff18g4extra_NoEtaEff" +  TypeAnalysis[AnalysisType]+".root";
      else    VarName[i] =  TypeAnalysis[AnalysisType]+".root";
      NameHisto[i] = "";
    }
    else if (RunVar==13){
      if (i==0)   VarName[i] = ".root";
      else    VarName[i] = "_AllEff2019.root";
      NameHisto[i] = "";
    }
    else if (RunVar==14){
      if (i>=0 && i<numFiles/2)   VarName[i] = "1617MC_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      else if (i>=numFiles/2)   VarName[i] = "17pq_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==15){
      if (i==0)  {
	//oldAODs	VarName[i] = "PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_" +  TypeAnalysisBis[AnalysisType]+ ".root";
	//old AODs      NameHisto[i] = "fHistYieldStat";

	//	VarName[i]+="SystematicAnalysis1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	//	NameHisto[i] += "fHistYieldvsErrSoloStat"; 

	if (type==0) {
	  VarName[i] += "PtSpectraBis_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_" + TypeAnalysisBis[AnalysisType]+"_isNormCorrFullyComputed.root";
	  NameHisto[i]+= "fHistYieldStat";
	}
	else if(type==1){
	  //	  VarName[i] += "PtSpectraBis_Xi_Eta0.8_PtMin3.0_" + TypeAnalysisBis[AnalysisType]+".root";
	  VarName[i] += "SystematicAnalysisRun2DataRed_MECorr_hXi_Jet0.75_Xi_Eta0.8_"+ TypeAnalysis[AnalysisType] +"Data_PtMin3.0.root"; 
	  //NameHisto[i]+= "fHistYieldStat";
	  NameHisto[i] = "fHistYieldvsErrSoloStat";
	}
      }
      else {
	if (type==0){
	  //	VarName[i] = "SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_" +  TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff_isNormCorr.root";
	  VarName[i] = "SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_" +  TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff_MultBinning3_isNormCorr.root";
	  NameHisto[i] = "fHistYieldvsErrSoloStat";
	}
	else if (type==1){
	  VarName[i] = "SystematicAnalysis17pq_hXi_Xi_Eta0.8_" + TypeAnalysis[AnalysisType]+ "Data_PtMin3.0_IsEtaEff_MultBinning3.root";
	  if (AnalysisType==0)  	  VarName[i] = "SystematicAnalysis17pq_hXi_OOJNoTriggerSmoothed_Xi_Eta0.8_" + TypeAnalysis[AnalysisType]+ "Data_PtMin3.0_IsEtaEff_MultBinning3.root";
	  NameHisto[i] = "fHistYieldvsErrSoloStat";
	}
      }
    }
    else if (RunVar==16){
      if (i<numFiles/2)  {
	//old AODs	VarName[i] = "PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_" +  TypeAnalysisBis[AnalysisType]+ ".root";
	//old AODs	NameHisto[i] = "fHistSpectrum"+Smolt[numEff];
	VarName[i]+="SystematicAnalysis1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	NameHisto[i] = Form("fHistSpectrumPart_m%i_syst0",numEff );
      }
      else {
	VarName[i] = "SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_" +  TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	NameHisto[i] = Form("fHistSpectrumPart_m%i_syst0",numEff );
      }
      cout << " molt " << Smolt[numEff] << endl;
    }
    else if (RunVar==17){
      if (i<numFiles/2)  {
	VarName[i] = Form("PtBinning1_1617_hK0s_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt%i_sysT0_sysV00_Sys0_PtMin3.0.root", numEff);
	NameHisto[i] = "";
      }
      else {
	VarName[i] = Form("PtBinning1_17pq_hK0s_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt%i_sysT0_sysV00_Sys0_PtMin3.0.root", numEff);
	NameHisto[i] = "";
      }
    }
    else if (RunVar==18){
      if (i<numFiles/2)  {
	VarName[i] = "CENTwoSDD_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
      else {
	VarName[i] = "FAST_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==19){
      VarName[i] = "";
      NameHisto[i] = Smolt[numEff];
    }
    else if (RunVar==20 || RunVar==21 || RunVar==22 || RunVar==23){
      if (i<numFiles/2)  {
	if (type==0)	VarName[i] = "hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	//	else if (type==1) VarName[i] = "AllMC_hXi_EtaEff_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	else if (type==1) VarName[i] = "161718_AOD235_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	//	VarName[i] = "AllMC_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	NameHisto[i] = Smolt[numEff];
      }
      else {
	if (type==0)	VarName[i] = "AOD235_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	else if (type==1) 	VarName[i] = "161718Full_AOD235_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	NameHisto[i] = Smolt[numEff];
      }
    }
    else if (RunVar==24 || RunVar==25){
      if (i<numFiles/2)  {
	if (i==0) VarName[i] = "2018f1";
	else if (i==1) VarName[i] = "2017e5";
	VarName[i] += "_extra_AOD235_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	NameHisto[i] = Smolt[4];
      }
      else {
	if (i==numFiles/2) VarName[i] = "2018f1";
	else if (i==numFiles/2 +1) VarName[i] = "2017e5";
	VarName[i] += "_extra_AOD235_hK0s_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0_isEfficiencyMassSel.root";
	NameHisto[i] = Smolt[4];
      }
    }
    else if (RunVar==26 || RunVar==28){ //K0s
      if (i<numFiles/2)  {
	VarName[i] += "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_"+ TypeAnalysisBis[AnalysisType] +"_isNormCorrFullyComputed.root"; //DATA
	NameHisto[i] = "fHistSpectrum" + Smolt[numEff];
      }
      else {
	VarName[i] += "FinalOutput/DATA2016/SystematicAnalysis1617GP_hK0s_PtBinning1_K0s_Eta0.8_AllAssoc_"+TypeAnalysis[AnalysisType] +"_MCTruth_PtMin3.0.root"; //MC
	NameHisto[i] = Form("fHistSpectrumPart_m%i_syst0", numEff);
      }
    }
    else if (RunVar==27 || RunVar==29){ //Xi
      if (i<numFiles/2)  {
	VarName[i] += "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_"+TypeAnalysisBis[AnalysisType] + "_isNormCorrFullyComputed.root"; //DATA
	NameHisto[i] = "fHistSpectrum" + Smolt[numEff];
      }
      else {
	VarName[i] += "FinalOutput/DATA2016/SystematicAnalysis161718_hXi_Xi_Eta0.8_AllAssoc_"+TypeAnalysis[AnalysisType] +"_MCTruth_PtMin3.0.root"; //MC
	NameHisto[i] = Form("fHistSpectrumPart_m%i_syst0", numEff);
      }
    }
    else if (RunVar==30){ 
      if (i<numFiles/2)  {
	if (type==0)	VarName[i] += "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_K0s_Eta0.8_PtMin3.0_"+ TypeAnalysisBis[AnalysisType] +"_isNormCorrFullyComputed.root"; //DATA 
	else VarName[i] += "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_"+TypeAnalysisBis[AnalysisType] + "_isNormCorrFullyComputed.root"; //DATA
	NameHisto[i] = "fHistYieldStat";
      }
      else {
	if (type==0) 	VarName[i] += "FinalOutput/DATA2016/SystematicAnalysis1617GP_hK0s_PtBinning1_K0s_Eta0.8_AllAssoc_"+TypeAnalysis[AnalysisType] +"_MCTruth_PtMin3.0.root"; //MC
	else VarName[i] += "FinalOutput/DATA2016/SystematicAnalysis161718_hXi_Xi_Eta0.8_AllAssoc_"+TypeAnalysis[AnalysisType] +"_MCTruth_PtMin3.0.root"; //MC
	NameHisto[i] = "fHistYieldvsErrSoloStat";
      }
    }
    else if (RunVar==31){ 
      if (i==0) VarName[i]+="_2018f1_extra_AOD235_hK0s_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      else if (i==1) VarName[i]+="_LHC16_GP_AOD235_17f9_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      else if (i==2) VarName[i]+="_LHC16_GP_AOD235_17e5_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      else if (i==3) VarName[i]+="_LHC16_GP_AOD235_17d16_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
      else if (i==4) VarName[i]+="_LHC17_GP_AOD235_18c12_K0s_Eta0.8_isMeanFixedPDG_BkgRetta_molt5_sysT0_sysV00_Sys0_PtMin3.0.root";
    }
    else if (RunVar==32){ 
      if (i>=numFiles/2) VarName[i]+="_AllBut18c12";
      VarName[i]+="_PtBinning1_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      NameHisto[i] += Smolt[numEff];
    }
    else if (RunVar==33){ 
      if (i>=numFiles/2) {
	if (type == 0)	VarName[i]+="FinalOutput/DATA2016/SystematicAnalysis1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0.root";
	else if (type==1) 	VarName[i]+="FinalOutput/DATA2016/SystematicAnalysis161718Full_AOD234_hXi_Xi_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0.root";
	NameHisto[i] +=	Form("fHistSpectrumPart_m%i_syst0", numEff);

      }
      else {
	if (type==0)	VarName[i]+= "FinalOutput/DATA2016/PtSpectraBis_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_" + TypeAnalysisBis[AnalysisType]+ ".root";
	else VarName[i] += "FinalOutput/DATA2016/PtSpectraBis_Xi_Eta0.8_PtMin3.0_" + TypeAnalysisBis[AnalysisType] +".root";
	NameHisto[i] += "fHistSpectrum"+ Smolt[numEff];
      }
    }
    else if (RunVar==34){ 
      if (i>=numFiles/2) {
	if (type==0)	VarName[i]+="FinalOutput/DATA2016/SystematicAnalysis1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	else 	VarName[i]+="FinalOutput/DATA2016/SystematicAnalysis161718Full_AOD234_hXi_Xi_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0_IsEtaEff.root";
	NameHisto[i] +=	Form("fHistSpectrumPart_m%i_syst0", numEff);
      }
      else {
	if (type==0)VarName[i]+="FinalOutput/DATA2016/SystematicAnalysis1617_AOD234_hK0s_PtBinning1_K0s_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0.root";
	else 	VarName[i]+="FinalOutput/DATA2016/SystematicAnalysis161718Full_AOD234_hXi_Xi_Eta0.8_" + TypeAnalysis[AnalysisType] +"Data_PtMin3.0.root";
	NameHisto[i] +=	Form("fHistSpectrumPart_m%i_syst0", numEff);
      }
    }
    else if (RunVar==36){ 
      if (i==0) VarName[i]+="Default_Sys0_PtMin3.0.root";
      else if (i==1) VarName[i] += "Tightest_Sys0_PtMin3.0.root";
      else if (i==2) VarName[i] += "Loosest_Sys0_PtMin3.0.root";
      else if (i==3) VarName[i] += "0_Sys0_PtMin3.0.root";
    }
    else if (RunVar==37){
      if (i<numFiles/2)  {
	VarName[i] = "2018f1_extra_hK0s_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
	NameHisto[i] = Smolt[numEff];
      }
      else {
	VarName[i] = "15g3c3_hK0s_PtBinning1_K0s_y0.5_AllAssoc_SysT0_SysV00_PtMin0.2.root";
	NameHisto[i] = Smolt[numEff];
      }
    }
    else if (RunVar==38){
      if (i<numFiles/2)  {
	VarName[i] = Form("hK0sHMMCOld/AOD230/Efficiency/Efficiency2019h11_extra_Child%i_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root", i+1);
	NameHisto[i] = "_0-0.1";
      }
      else {
	VarName[i] = Form("hK0sHMMCCorrect/Efficiency/EfficiencyLHC19h11_Child%i_K0s_Eta0.8_SysT0_SysV00_PtMin3.0.root", i+1-3);
	NameHisto[i] = "_0-0.1";
      }
    }
    else if (RunVar==39){
      NameHisto[i] = Form("m%i_syst0", numEff);
      if (i<numFiles/2)  {
	VarName[i] = "K0s_Eta0.8_JetData_PtMin3.0.root";
      }
      else {
	VarName[i] = "OOJNoTriggerSmoothed_K0s_Eta0.8_JetData_PtMin3.0.root";
	//	VarName[i] = "OOJNoTriggerSmoothedCorrMult_K0s_Eta0.8_JetData_PtMin3.0.root";
	//	VarName[i] = "OOJSmoothedBis_K0s_Eta0.8_JetData_PtMin3.0.root";
	//	VarName[i] = "OOJSmoothedBisCorrMult_K0s_Eta0.8_JetData_PtMin3.0.root";
      }
    }
    else if (RunVar==40){
      if (i<numFiles/2)  { //K0s
	if (i<numFiles/4) { //13 TeV
	  VarName[i] += "PtSpectraBis_PtBinning1_1617_AOD234_hK0s_K0s_Eta0.8_PtMin3.0_" + TypeAnalysisBis[i]+"_isNormCorrFullyComputed.root";
	  NameHisto[i]+= "fHistYieldStat";
	}
	else { //5 TeV
	  VarName[i] = "SystematicAnalysis17pq_hK0s_PtBinning1_K0s_Eta0.8_" +  TypeAnalysis[i-numFiles/4] +"Data_PtMin3.0_IsEtaEff_MultBinning3_isNormCorr.root";
	  NameHisto[i] = "fHistYieldvsErrSoloStat";
	}
      }
      else { //Xi
	if (i<3*numFiles/4){
	  VarName[i] += "PtSpectraBis_Xi_Eta0.8_PtMin3.0_" + TypeAnalysisBis[i-numFiles/2]+".root";
	  NameHisto[i] = "fHistYieldStat";
	  //	  VarName[i] += "SystematicAnalysisRun2DataRed_MECorr_hXi_Jet0.75_Xi_Eta0.8_"+ TypeAnalysis[i-numFiles/2] +"Data_PtMin3.0.root"; 
	  //	  NameHisto[i] = "fHistYieldvsErrSoloStat";
	}
	else { //5 TeV
	  VarName[i] = "SystematicAnalysis17pq_hXi_Xi_Eta0.8_" + TypeAnalysis[i-numFiles/4*3]+ "Data_PtMin3.0_IsEtaEff_MultBinning3.root";
	  if (i==3*numFiles/4) VarName[i] = "SystematicAnalysis17pq_hXi_OOJNoTriggerSmoothed_Xi_Eta0.8_" + TypeAnalysis[AnalysisType]+ "Data_PtMin3.0_IsEtaEff_MultBinning3.root";
	  NameHisto[i] = "fHistYieldvsErrSoloStat";
	}
      }
    }
    else if (RunVar==41){
      NameHisto[i] = Smolt[multDef];
      if (i == 0)  {
	VarName[i] = "LHC16_GP_AOD235";
      }
      else if (i ==1) {
	VarName[i] = "LHC17_GP_AOD235";
      }
      else {
	VarName[i] = "LHC18_GP_AOD235";
      }
      VarName[i] += "_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
    }
    else if (RunVar==42){
      NameHisto[i] = Smolt[numEff];
      if (i<numFiles/2)  { //low pt,trig
	if (type==1)	VarName[i]  = "161718_MD_EtaEff_LowPtTrig_hXi_PtTrigMax2.5_Xi_Eta0.8_SysT0_SysV00_PtMin0.2.root";

      }
      else { //pt,trig > 3
	//	if (type==1)	VarName[i]  = "161718Full_AOD235_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
	VarName[i] = "161718_MD_EtaEff_PtTrig3_hXi_Xi_Eta0.8_SysT0_SysV00_PtMin3.0.root";
      }
    }
    else if (RunVar==43){
      if (i==0)  {sys=0; sysPhi=0;}
      if (i==1)  {sys=0; sysPhi=1;}
      if (i==2)  {sys=0; sysPhi=2;}
      if (i==3)  {sys=4; sysPhi=0;}
      if (i==4)  {sys=4; sysPhi=1;}
      if (i==5)  {sys=4; sysPhi=2;}
      if (i==6)  {sys=5; sysPhi=0;}
      if (i==7)  {sys=5; sysPhi=1;}
      if (i==8)  {sys=5; sysPhi=2;}
      if (i==9)  {sys=6; sysPhi=0;}
      if (i==10) {sys=6; sysPhi=1;}
      if (i==11) {sys=6; sysPhi=2;}
      if (sys!=0 && sysPhi!=0)      VarName[i] = Form("_sys%i_sysPhi%i.root", sys, sysPhi);
      else if (sys!=0)      VarName[i] = Form("_sys%i.root", sys);
      else if (sysPhi!=0)   VarName[i] = Form("_sysPhi%i.root", sysPhi);
      else   VarName[i] = ".root";
      NameHisto[i] = Form("%i", multDef);
    }

    if (RunVar==43 && sys==5 && AnalysisType==0) continue;

    InputName = CommonFileName + VarName[i];
    //    if (RunVar==4) InputName += "_IsOnlypiKpemu";
    if (RunVar<6) InputName+="_IsEstimateRun3.root";
    InputFile = new TFile (InputName, "");
    if (!InputFile) return;
    histo[i] = (TH1F*)InputFile->Get(histoName + NameHisto[i]);
    cout << "histo name " << histoName << NameHisto[i]<< endl;
    if (!histo[i]) {cout << "histogram is not there: " << histoName << NameHisto[i] << endl; return;}
    histo[i] ->SetName(Form("histoName%i", i));
    histo[i]->Sumw2();

    canvas->cd(1);
    if (RunVar<7 || RunVar==16 || (RunVar>=26 && RunVar<=29) || RunVar==33 || RunVar==34)    gPad->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    if (TypeOfComparison==2 && i>=numFiles/2) style = 33;
    if (TypeOfComparison==3 && i>=numFiles/2) style = 33;
    if (TypeOfComparison==1) style = 33;
    if (RunVar==19) style = 33;
    if (RunVar==10 || RunVar==12 || RunVar==13 || RunVar==15) color[numEff-numDef]=1;
    if (RunVar==30 || RunVar==15) {
      if (i==0) color[numEff-numDef]= kRed+2;
      else  color[numEff-numDef]= kBlue-3;
    }
    Float_t msize = 1;
    if (RunVar==30) msize = 2;
    if (RunVar==40) msize = 2;
    StyleHistoYield(histo[i], Low, Up, color[numEff-numDef], style, titleX, titleY, title, msize, 1.2, 1.5);
    if ((RunVar==26 || RunVar==27) && AnalysisType==0){
      //      histo[i]->GetYaxis()->SetLabelSize(0.03);
      //      histo[i]->GetYaxis()->SetTitleSize(0.04);
    }
    if (RunVar==30){
     histo[i]->GetYaxis()->SetLabelSize(0.03);
     histo[i]->GetYaxis()->SetTitleSize(0.04);
    }
    if (RunVar<=6) LegendName[i] = Form("%.1f GeV/#it{c}", num);
    if (RunVar==2){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}, only TOF PU rej", (float)numEff);
    }
    else if (RunVar==3){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}, Run3 sel", (float)numEff);
    }
    else if (RunVar==5){
      if (i<numFiles/2)    LegendName[i] = Form("%.1f GeV/#it{c}, AllAsso", (float)numEff);
      else    LegendName[i] = Form("%.1f GeV/#it{c}", (float)numEff);
    }
    else if (RunVar==6){
      if (i==0)  LegendName[i] = " Run2 sel";
      if (i==1)  LegendName[i] = " only-TOF OOBpu rej";
      if (i==2)  LegendName[i] = " Run3 sel";
    }
    else if (RunVar==7){
      if (i<numFiles/2)    LegendName[i] = "New AODs " + SmoltBis[numEff];
      else    LegendName[i] = "Old AODs " + SmoltBis[numEff];
    }
    else if (RunVar==8){
      if (i<numFiles/2)    LegendName[i] = "New " + SmoltBis[numEff];
      else    LegendName[i] = "Old " + SmoltBis[numEff];
    }
    else if (RunVar==9){
      LegendName[i] = "New " + SmoltBis[numEff];
    }
    else if (RunVar==10){
      if (i==0)   LegendName[i] = "16k ";
      else if (i==1)    LegendName[i] = "AllBut16k ";
    }
    else if (RunVar==11){
      if (i<numFiles/2)    LegendName[i] = "1617MC " + SmoltBis[numEff];
      else    LegendName[i] = "18g4_extra " + SmoltBis[numEff];
    }
    else if (RunVar==12){
      if (i==0)   LegendName[i] = "18g4_extra "; //denominator
      else if (i==1)    LegendName[i] = "18g4_extra -- eta-dep eff"; //numerator
    }
    else if (RunVar==13){
      if (i==0)   LegendName[i] = "Eff: 2019h11c_extra "; //denominator
      else if (i==1)    LegendName[i] = "Eff: 2019h11abc_extra"; //numerator
    }
    else if (RunVar==14){
      if (i<numFiles/2)   LegendName[i] = "Eff: 13 TeV "+ SmoltBis[numEff]; //denominator
      else   LegendName[i] = "Eff: 5 TeV" + SmoltBis[numEff]; //numerator
    }
    else if (RunVar==15){
      if (i==0)   LegendName[i] = "13 TeV "; //denominator
      else if (i==1)    LegendName[i] = "5 TeV "; //numerator
    }
    else if (RunVar==16){
      if (i<numFiles/2)    LegendName[i] = "13 TeV" + SmoltBis[numEff];
      else    LegendName[i] = "5 TeV " + SmoltBis[numEff];
    }
    else if (RunVar==17){
      if (i<numFiles/2)    LegendName[i] = "13 TeV" + SmoltBis[numEff];
      else    LegendName[i] = "5 TeV " + SmoltBis[numEff];
    }
    else if (RunVar==18){
      if (i<numFiles/2)    LegendName[i] = "CENTwoSDD " + SmoltBis[numEff];
      else    LegendName[i] = "FAST " + SmoltBis[numEff];
    }
    else if (RunVar==19){
      LegendName[i] = SmoltBis[numEff];
    }
    else if (RunVar==20 || RunVar==21 || RunVar==22 || RunVar==23){
      if (i<numFiles/2)    LegendName[i] = "Old AODs " + SmoltBis[numEff];
      else    LegendName[i] = "AOD 235 " + SmoltBis[numEff];
    }
    else if (RunVar==24 || RunVar==25){
      if (i==0)    LegendName[i] = "2018f1_extra";
      else if (i==1) LegendName[i] = "2017e5_extra";
      else if (i==2) LegendName[i] = "2018f1_extra MassSel";
      else if (i==3) LegendName[i] = "2017e5_extra MassSel";
    }
    else if (RunVar>=26 && RunVar <=29){
      if (i<numFiles/2)    LegendName[i] = tipo[type] + " DATA " + SmoltBis[numEff];
      else    LegendName[i] = tipo[type] + " MC " + SmoltBis[numEff];
    }
    else if (RunVar==30){
      if (i<numFiles/2)    LegendName[i] = tipo[type] + " DATA ";
      else    LegendName[i] = tipo[type] + " MC ";
    }
    else if (RunVar==31){
      if (i==0)    LegendName[i] = "18f1";
      else if (i==1) LegendName[i] = "17f9";
      else if (i==2) LegendName[i] = "17e5";
      else if (i==3) LegendName[i] = "17d16";
      else if (i==4) LegendName[i] = "18c12";
    }
    else if (RunVar==32){
      if (i<numFiles/2)  LegendName[i] = "1617_GP_AOD235" + SmoltBis[numEff];
      else LegendName[i] = "1617_GP_AOD235_AllBut18c12" + SmoltBis[numEff];
    }
    else if (RunVar==33){
      if (i<numFiles/2)  LegendName[i] = "Old AOD " + SmoltBis[numEff];
      else LegendName[i] = "AOD234 " + SmoltBis[numEff];
    }
    else if (RunVar==34){
      if (i<numFiles/2)  LegendName[i] = " " + SmoltBis[numEff];
      else LegendName[i] = "#eta dependent eff " + SmoltBis[numEff];
    }
    else if (RunVar==36){
      if (i==0)  LegendName[i] = "Default selections";
      else if (i==1)  LegendName[i] = "Tightest selections";
      else if (i==2)  LegendName[i] = "Loosest selections";
      else if (i==3)  LegendName[i] = "Default2 selections";
    }
    else if (RunVar==37){
      if (i<numFiles/2)    LegendName[i] = "2018f1_extra " + SmoltBis[numEff];
      else    LegendName[i] = "2015g3c3 " + SmoltBis[numEff];
    }
    else if (RunVar==38){
      if (i==0)   LegendName[i] = "AOD230 19h11a_extra";
      else if (i==1)   LegendName[i] = "AOD230 19h11b_extra";
      else if (i==2)   LegendName[i] = "AOD230 19h11c_extra";
      if (i==3)   LegendName[i] = "newAOD 19h11a_extra";
      else if (i==4)   LegendName[i] = "newAOD 19h11b_extra";
      else if (i==5)   LegendName[i] = "newAOD 19h11c_extra";
    }
    else if (RunVar==39){
      if (i<numFiles/2)   LegendName[i] = "Default OOJ sub " + SmoltBis[numEff];
      else LegendName[i] = "New OOJ scaled to 1<dphi<2" + SmoltBis[numEff];
    }
    else if (RunVar==40){
      if (i<numFiles/4)   LegendName[i] = "K0s yield 13 TeV " + TypeAnalysis[i];
      else if (i<numFiles/4*2)   LegendName[i] = "K0s yield 5 TeV " + TypeAnalysis[i-3];
      else if (i<numFiles/4*3)   LegendName[i] = "Xi yield 13 TeV " + TypeAnalysis[i-6];
      else  LegendName[i] = "Xi yield 5 TeV " + TypeAnalysis[i-9];
    }
    else if (RunVar==41){
      if (i == 0)   LegendName[i] = "2016 " + SmoltBis[multDef];
      else if (i == 1)   LegendName[i] = "2017 " + SmoltBis[multDef];
      else   LegendName[i] = "2018 " + SmoltBis[multDef];
    }
    else if (RunVar==42){
      if (i<numFiles/2)   LegendName[i] = "low pt,trig " + SmoltBis[multDef];
      else    LegendName[i] = " pt,trig > 3 " + SmoltBis[multDef];
    }
    else if (RunVar==43){
      if (i==0) LegendName[i] = "Default  " + SmoltBis[multDef];
    }
    legend->AddEntry(histo[i], LegendName[i], "pl");
    //    cout << " n bins " <<     histo[i]->GetNbinsX() << endl;
    pol1[i] = new TF1(Form("pol1_%i", i), "pol1", 0,30);
    if (RunVar==22 || RunVar==23){
      for (Int_t b=1; b<= histo[i]->GetNbinsX(); b++){
	histo[i]->SetBinError(b, 0);
      }
    }

    if (TypeOfComparison==3){
      if (i<numFiles/2){
	//	if (Smolt[numEff] == "_10-30"|| Smolt[numEff] == "__all" || Smolt[numEff] == "_30-50" )     histo[i]->Draw("same");
	//	if (Smolt[numEff] == "_10-30"|| Smolt[numEff] == "__all" )     histo[i]->Draw("same");
	//	if (Smolt[numEff] == "_50-100" || Smolt[numEff] == "_30-50"  )     histo[i]->Draw("same");
      }
      else {
	//	if (Smolt[numEff] == "_5-10"|| Smolt[numEff] == "__all" || Smolt[numEff] == "_10-30" )     histo[i]->Draw("same");
	//	if (Smolt[numEff] == "_50-100" || Smolt[numEff] == "_30-50" || Smolt[numEff] == "_10-30")     histo[i]->Draw("same");
      }
      histo[i]->Draw("same ep");
    }
    else histo[i]->Draw("same pe");

    if (RunVar==15 && AnalysisType!=0) {
      if (i==1) pol1[i]->SetLineStyle(8);
      pol1[i]->SetLineColor(1);
      pol1[i]->SetLineWidth(0.3);
      histo[i]->Fit(pol1[i], "R+");
    }
    if (i==numFiles-1) legend->Draw("");

    //    cout << " first canvas ok " << endl;
    canvas->cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    if (RunVar == 22 || RunVar == 23) {
      histo[i]->GetXaxis()->SetRangeUser(-0.8, 0.8);
      // histo[i] ->Rebin(2);
    }
    histoRatio[i] = (TH1F*) histo[i]->Clone(histoName + "_Ratio");

    if (TypeOfComparison==2){
      if (i>=numFiles/2){
	histoRatio[i]->Divide(histo[i-numFiles/2]);
	if (RunVar!=22 && RunVar!=23 && RunVar!=26 && RunVar!=27){
	  ErrRatioCorr(histo[i], histo[i-numFiles/2], histoRatio[i], 0);
	  CorrelationBtwHistos = 2;
	}
	if (RunVar==39 || RunVar==34) {
	  ErrRatioCorr(histo[i], histo[i-numFiles/2], histoRatio[i], 1);
	  CorrelationBtwHistos = 1;
	}
	if (RunVar==22 && RunVar==23) {
	  for (Int_t b=1; b<= histoRatio[i]->GetNbinsX(); b++){
	    histoRatio[i]->SetBinError(b, 0);
	  }
	}
      }
    }
    else if (TypeOfComparison==3){
      if (i<numFiles/2){
	cout << " I divide " <<   histo[i]->GetName() << " by " << histo[0]->GetName() << endl;
	histoRatio[i]->Divide(histo[0]);
	ErrRatioCorr(histo[i], histo[0], histoRatio[i], 0);
	CorrelationBtwHistos = 2;
      }
      else {
	cout << " I divide " <<   histo[i]->GetName() << " by " << histo[numFiles/2]->GetName() << endl;
	histoRatio[i]->Divide(histo[numFiles/2]);
	//	histoRatio[i]->Divide(histo[0]);
	ErrRatioCorr(histo[i], histo[numFiles/2], histoRatio[i], 0);
	CorrelationBtwHistos = 2;
      }
      for (Int_t b=1; b<= histoRatio[i]->GetNbinsX() ; b++){
	cout << "num " << histo[i]->GetBinContent(b) << " +- " <<  histo[i]->GetBinError(b)<< " (" << histo[i]->GetBinError(b)/histo[i]->GetBinContent(b)<<  ")" <<endl;
	if (i<numFiles/2) cout << "denom " << histo[0]->GetBinContent(b) << " +- " <<  histo[0]->GetBinError(b)<< " (" <<histo[0]->GetBinError(b)/histo[0]->GetBinContent(b)  << ")" << endl;
	else cout << "denom " << histo[numFiles/2]->GetBinContent(b) << " +- " <<  histo[numFiles/2]->GetBinError(b)<< " (" <<histo[numFiles/2]->GetBinError(b)/histo[0]->GetBinContent(b)  << ")" << endl;
	cout << "ratio " <<histoRatio[i]->GetBinContent(b) << " +- " <<  histoRatio[i]->GetBinError(b)<<endl;
      }
    }
    else {
      if (i!=0){
	histoRatio[i]->Divide(histo[0]); //no corr
	CorrelationBtwHistos = 0;
	if (RunVar==12 || RunVar==13)	{
	  ErrRatioCorr(histo[i], histo[0], histoRatio[i], 1); //full corr
	  CorrelationBtwHistos = 1;
	}
	else if (RunVar!=15 && RunVar!=30 && RunVar!=31 && RunVar!=41) {
	  ErrRatioCorr(histo[i], histo[0], histoRatio[i], 0); //partial corr
	  CorrelationBtwHistos = 2;
	}
      }
    }

    StyleHistoYield(histoRatio[i], LowRatio, UpRatio, color[numEff-numDef], style, titleX, "Ratio", titleRatio, msize, 1.2, 1.4);
    if (RunVar==39) histoRatio[i]->GetXaxis()->SetRangeUser(0, 2.5);
    if (RunVar==40) histoRatio[i]->GetXaxis()->SetRangeUser(0, 25);

    //old    if (RunVar==2 || RunVar==3 || RunVar==5 || RunVar==7 || RunVar==8 || RunVar==9 || RunVar==11 || RunVar==14 || RunVar==18 || RunVar==19) {
    if (TypeOfComparison==2) {
      if (i>=numFiles/2) {
	histoRatio[i]->Draw("same pe");
	//	if (RunVar==8 || RunVar==14 || RunVar==18 || RunVar==22) pol0At1->Draw("same");
	pol0At1->Draw("same");
      }
    }
    else if (TypeOfComparison==3){
      if (numEff!=5 && numEff!=11){
	//	if (i<numFiles/2) histoRatio[i]->Draw("same pe");
	histoRatio[i]->Draw("same pe");
      }
      if (i==numFiles-1) {
	//	legend->Draw("");
	pol0At1->Draw("same");
      }
    }
    else {
      if (i!=0) histoRatio[i]->Draw("same");
      //      if (i==numFiles-1) legend->Draw("");
      if (i!=0 && (RunVar==10 || RunVar==13 || RunVar==15 || RunVar==19 || RunVar==30 || RunVar==31 || RunVar==41)) pol0At1->Draw("same");
    }

    hdummy[i]= new TH1F (Form("hdummy%i", i), Form("hdummy%i", i), 1000, 0, 30);
    if (RunVar==15 && AnalysisType!=0) {
      hdummy[i]->Add(pol1[i]);
      if (i==1)  {
	hdummy[1]->Divide(hdummy[0]);
	hdummy[1]->Draw("same");
      }
    }
  
    IsBarlowSign=0;
    if (isBarlow){
      canvasB->cd(1);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      histoBarlow[i] = (TH1F*) histo[i]->Clone(histoName + "_Barlow");
      histoSysError[i] = (TH1F*) histo[i]->Clone(histoName + "_Syst");
      StyleHistoYield(histoBarlow[i], LowBarlow, UpBarlow, color[numEff-numDef], style, titleX,"Barlow variable" , "Barlow variable", 2, 1.2, 1.5);
      StyleHistoYield(histoSysError[i], 0, UpError, color[numEff-numDef], style, titleX,"Relative syst. uncertainty (%)" , "Syst. uncertainty", msize, 1.2, 1.5);
      if (TypeOfComparison==2 && i>=numFiles/2) {
	BarlowVariable(histo[i], histo[i-numFiles/2], histoBarlow[i], histoSysError[i], 2 ,2, IsBarlowSign);
	histoBarlow[i]->GetXaxis()->SetRangeUser(0,2.5);
	histoBarlow[i]->Draw("same p");
	lineatB0->Draw("same");
	if (i==numFiles-1) legend->Draw("");
      }
      else if (TypeOfComparison ==1 && i>0){
	BarlowVariable(histo[i], histo[0], histoBarlow[i], histoSysError[i], NSigma ,NSign, IsBarlowSign);
	cout << "Variation n. "<< i << " Is the variation Barlow significant? "<< IsBarlowSign << endl;
	histoBarlow[i]->GetXaxis()->SetRangeUser(0,2.5);
	histoBarlow[i]->Draw("same p");
	lineatB0->Draw("same");
	if (i==numFiles-1) legend->Draw("");
      }

      canvasB->cd(2);
      gPad->SetLeftMargin(0.15);
      gPad->SetBottomMargin(0.15);
      if ((TypeOfComparison==2 && i>=numFiles/2) || (TypeOfComparison ==1 && i!=0)) {
	//	histoSysError[i]->GetXaxis()->SetRangeUser(0,2.5);
	histoSysError[i]->GetYaxis()->SetRangeUser(0,LimSupSys);
	histoSysError[i]->Draw("same");
      }
    }
  }
  cout << "\n\n"<< endl;
  canvas->SaveAs(OutputNamepdf);
  OutputFile->WriteTObject(canvas);
  OutputFile->WriteTObject(canvasB);
  OutputFile->Close();

  cout << "I produced the output file " << OutputNameRoot << " and " << OutputNamepdf << endl;

  cout << "\nCorrelation between numerator and denominator histos:" << endl;
  cout << SCorrelation[CorrelationBtwHistos] << endl;

  cout << "\nTypeOfComparison: " <<     TypeOfComparison<< endl;
}
