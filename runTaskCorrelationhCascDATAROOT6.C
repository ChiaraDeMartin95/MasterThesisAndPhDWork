void runTaskCorrelationhCascDATAROOT6(
				  const char* runtype = "grid", // local, proof or grid
				  const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
				  const char *AssocParticle = "Xi",
				  const bool isMC = 0, // 1 = MC, 0 = DATA
				  const bool isEff=0, //1 = MC per efficienza, 0 = MC truth
				  const bool isHybridMCTr=0,
				  const bool isInclusiveAnalysis = 0,
				  const bool isHM =0,
				  const float PercentileMin = 0.,
				  const float PercentileMax = 100.,
				  const Float_t EtaTrigger=0.8,
				  const Float_t EtahAssoc=0.8,
				  const Float_t EtaV0Assoc=0.8,
				  const Int_t FilterBitValue=128,
				  const Int_t EvtToMix = 5,
				  const Int_t year=2016,
				  const Long64_t nentries = 1234567890,//2000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
				  const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
				  const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
				  const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
				  const char *taskname = "Xijetpp_LHC16k" //LHC18f1_extra_15runsNOFB_150MeV"*/   // sets name of grid generated macros
		   
)
{

  Bool_t islocal= kFALSE;
  Bool_t isTestMode = kFALSE;
  if (strcmp(runtype, "local")==0)  islocal = kTRUE; 
  if (strcmp(gridmode, "test")==0)  isTestMode = kTRUE;

  Printf("%s analysis chosen",runtype);
  if (isMC==1 && isEff==1) Printf("MC analysis chosen: MC analyzed as data and efficiency computed");
  if (isMC==1 && isEff==0) Printf("MC analysis chosen: MC truth analyzed");
  if (isMC==1 && isHybridMCTr==1) Printf("MC analysis chosen: MC truth analyzed in hybrid mode (reco trigger particles, gen associated)");
  if (isMC==0) Printf("data analysis chosen");

  // load libraries
  gSystem->Load("libCore");        
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWGPP");
  gSystem->Load("libANALYSISaliceBase");
  gSystem->Load("libCORRFW");
  gSystem->Load("libOADB");
  //  gSystem->Load("libAliPythia6"), this makes the code break

    
  // since we will compile a class, tell root where to look for headers  
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/STEER");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/ANALYSIS");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/ESD");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/STEER/ESD");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
  gInterpreter->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER");
  gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS");
  gROOT->ProcessLine(".include $ALICE_ROOT/ESD");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER/ESD");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
#endif
     
  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  /*
  //PhysicsSelection Configuration
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* ps =  reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%i)", isMC)));
  if (!ps) { Printf("no physics selection task"); return;}
  */

  //MultSelection
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* ms =  reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
  //ms->SetAddInfo(kTRUE);
  if(!ms) { Printf("no multiplictyTask"); return; }
  
  //PID Response
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  const bool isPIDResponseMCTruth = isMC;
  cout << "\nisPIDResponseMCTruth " << isPIDResponseMCTruth<< endl;
  AliAnalysisTaskPIDResponse* pid = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro(Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%i)", isPIDResponseMCTruth)));
  if(!pid ) { Printf("no PID response task"); return; }
  cout << "PID task loaded \n \n" << endl;
  
  // compile the class and load the add task macro
  // here we have to differentiate between using the just-in-time compiler
  // from root6, or the interpreter of root5
  AliAnalysisTaskCorrelationhCascDATA *task;
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskCorrelationhCascDATA.cxx++g");
  gInterpreter->LoadMacro("AliAnalysisCorrelationEventCollection.cxx++g");
  task = reinterpret_cast<AliAnalysisTaskCorrelationhCascDATA*>(gInterpreter->ExecuteMacro(Form("AddTaskCorrelationhCascDATA.C(%.f, %d, %i, %i, %i, %i)", 0.15, 15, islocal, isMC, isEff, isHybridMCTr)));  
#else
  gROOT->LoadMacro("AliAnalysisTaskCorrelationhCascDATA.cxx++g");
  gROOT->LoadMacro("AddTaskCorrelationhCascDATA.C");
  gROOT->LoadMacro("AliAnalysisCorrelationEventCollection.cxx++g");
  task = AddTaskCorrelationhCascDATA();
#endif

  task->SetMinPt(0.15);
  task->SetMaxPt(15);
  task->SetMC(isMC);
  task->SetEff(isEff);
  task->SetHybridMCTruth(isHybridMCTr);
  task->SetEvtToMix(EvtToMix);
  task->SetEtaTrigger(EtaTrigger);
  task->SetEtahAssoc(EtahAssoc);
  task->SetEtaV0Assoc(EtaV0Assoc);
  task->SetFilterBit(FilterBitValue);
  task->SetYear(year);
  task->SetInclusive(isInclusiveAnalysis);
  task->SetHM(isHM);
  task->SetMinimumMultPercentile(PercentileMin);
  task->SetMaximumMultPercentile(PercentileMax);
  task->SetAssocParticle(AssocParticle);

  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  // mgr->SetUseProgressBar(1, 25);

  if(islocal) {
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    //    TChain* chain = new TChain("TE");
    // add a few files to the chain (change this so that your local files are added)
    chain->Add("AliAOD_MCOld.root");
    //    chain->Add("galice.root");
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
  } else {
    // if we want to run on grid, we create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("libAliPythia6 libpythia6_4_21.so AliAnalysisCorrelationEventCollection.h AliAnalysisTaskCorrelationhCascDATA.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskCorrelationhCascDATA.cxx AliAnalysisCorrelationEventCollection.cxx");
    // select the aliphysics version. all other packages
    // are LOADED AUTOMATICALLY!
    alienHandler->SetAliPhysicsVersion("vAN-20220501_ROOT6-1");
    // set the Alien API version
    alienHandler->SetAPIVersion("V1.1x");
    // select the input data
    alienHandler->SetGridDataDir("/alice/data/2016/LHC16k/");
    alienHandler->SetDataPattern("/pass2/AOD234/*/AliAOD.root");
    // MC has no prefix, data has prefix 000
    alienHandler->SetRunPrefix("000");
    // runnumber
    Int_t nrun[10]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387};
    for(Int_t i = 0; i<10; i++){
      alienHandler->AddRunNumber(nrun[i]);
    }
    // number of files per subjob
    alienHandler->SetSplitMaxInputFileNumber(40);
    alienHandler->SetExecutable(Form("%s_2d.sh",taskname));
    //set number of runs in subjob
    alienHandler->SetNrunsPerMaster(15);
    // specify how many seconds your job may take
    alienHandler->SetTTL(10000);
    alienHandler->SetJDLName(Form("%s_2d.jdl",taskname));

    alienHandler->SetOutputToRunNo(kTRUE);
    alienHandler->SetKeepLogs(kTRUE);
    // merging: run with kTRUE to merge on grid
    // after re-running the jobs in SetRunMode("terminate") 
    // (see below) mode, set SetMergeViaJDL(kFALSE) 
    // to collect final results
    alienHandler->SetMaxMergeStages(1);
    alienHandler->SetMergeViaJDL(kTRUE);

    // define the output folders
    alienHandler->SetGridWorkingDir("myWorkingDir");
    alienHandler->SetGridOutputDir("myOutputDir");

    // connect the alien plugin to the manager
    mgr->SetGridHandler(alienHandler);
    if(isTestMode){
      // speficy on how many files you want to run
      alienHandler->SetNtestFiles(1);
      // and launch the analysis
      alienHandler->SetRunMode("test");
      mgr->StartAnalysis("grid");
    } else {
      // else launch the full grid analysis
      alienHandler->SetRunMode("full");
      mgr->StartAnalysis("grid");
    }
  }
}
