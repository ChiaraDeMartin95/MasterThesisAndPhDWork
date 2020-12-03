class AliAnalysisGrid;

//______________________________________________________________________________
void runTaskCorrelationhCasc(
	     const char* runtype = "grid", // local, proof or grid
	     const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
	     const bool isMC = 1, // 1 = MC, 0 = DATA 
	     const bool isEff=1, //1 = MC per efficienza, 0 = MC truth
	     const bool isHybridMCTr=0,
	     const Float_t EtaTrigger=0.8,
	     const Float_t EtahAssoc=0.8,
	     const Float_t EtaV0Assoc=0.8,
	     const Int_t FilterBitValue=128,
	     const Int_t year=2016,
	     const Long64_t nentries = 1234567890,//2000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
	     const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
	     const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
	     const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
	     const char *taskname = "Xijetpp_LHC18f1_extra_MCTruth"   // sets name of grid generated macros
	     ){
  // check run type
  if(runtype != "local" && runtype != "proof" && runtype != "grid"){
    Printf("\n\tIncorrect run option, check first argument of run macro");
    Printf("\tint runtype = local, proof or grid\n");
    return;
  }
  Printf("%s analysis chosen",runtype);
  if (isMC==1 && isEff==1) Printf("MC analysis chosen: MC analyzed as data and efficiency computed");
  if (isMC==1 && isEff==0) Printf("MC analysis chosen: MC truth analyzed");
  if (isMC==0) Printf("data analysis chosen");
  if (isMC==1 && isHybridMCTr==1) Printf("MC analysis chosen: MC truth analyzed in hybrid mode (reco trigger particles, gen associated)");
  const bool bMCtruth = isMC;             // 1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
  const bool bMCphyssel = isMC;           // 1 = looking at MC truth or reconstructed, 0 = looking at real data
  const Bool_t isPIDResponseMCtruth = isMC; // 1 = MC, 0 = DATA 
   
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
  cout<<"2"<<endl;
  
  // add aliroot indlude path
  
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER");
  gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS");
  gROOT->ProcessLine(".include $ALICE_ROOT/ESD");
  gROOT->ProcessLine(".include $ALICE_ROOT/STEER/ESD");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  //gROOT->ProcessLine(".include $ALICE_PHYSICS/STEER");
  //gROOT->ProcessLine(".include $ALICE_PHYSICS/ANALYSIS");
  
  gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  gROOT->SetStyle("Plain");
  
  cout<<"3"<<endl;
  
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
  
  // create the alien handler and attach it to the manager
  AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset); 
  mgr->SetGridHandler(plugin);

  //  Input

  AliAODInputHandler* iH = new AliAODInputHandler("handler","handler for my analisys");
  mgr->SetInputEventHandler(iH);
  
  // mc event handler
  if(bMCtruth) {
    AliMCEventHandler* mchandler = new AliMCEventHandler();
    // Not reading track references
    mchandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mchandler);
  }   
  
  
  // === Physics Selection Task ===
  //
  // In SelectCollisionCandidate(), default is kMB, so the task UserExec() 
  // function is only called for these events.
  // Options are:
  //    kMB             Minimum Bias trigger
  //    kMBNoTRD        Minimum bias trigger where the TRD is not read out
  //    kMUON           Muon trigger
  //    kHighMult       High-Multiplicity Trigger
  //    kUserDefined    For manually defined trigger selection
  //
  // Multiple options possible with the standard AND/OR operators && and ||
  // These all have the usual offline SPD or V0 selections performed.
  //
  // With a pointer to the physics selection object using physSelTask->GetPhysicsSelection(),
  // one can manually set the selected and background classes using:
  //    AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL")
  //    AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL");
  //
  // One can also specify multiple classes at once, or require a class to NOT
  // trigger, for e.g.
  //    AddBGTriggerClass("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL");
  //
  // NOTE that manually setting the physics selection overrides the standard
  // selection, so it must be done in completeness.
  //
  // ALTERNATIVELY, one can make the physics selection inside the task
  // UserExec().
  // For this case, comment out the task->SelectCol.... line, 
  // and see AliBasicTask.cxx UserExec() function for details on this.
  
  /*
  //========================================
  // OCBD Connect
  //========================================

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C");
  AliTaskConfigOCDB *CDBconnect = AddTaskConfigOCDB("raw://");
  */


  /*
  //=========================================
  // PHYSICS SELECTION
  //=========================================
  
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
  if(!physSelTask) { Printf("no physSelTask"); return; }
  */
  //==========================================
  // PID Task
  //==========================================
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  //  AliAnalysisTask *PIDTask = AddTaskPIDResponse(kFALSE,kTRUE);
  //  AliAnalysisTask *PIDTask = AddTaskPIDResponse(isPIDResponseMCTruth); //dipende da Dati o MC
  AliAnalysisTask *PIDTask = AddTaskPIDResponse(0); //dipende da Dati o MC
  if(!PIDTask ) { Printf("no PID response task"); return; }
  cout << "PID task loaded \n \n" << endl;

  
  //==========================================
  // MULTIPLICITY Task
  //==========================================
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* multSelectionTask = AddTaskMultSelection(kFALSE);
  if(!multSelectionTask) { Printf("no multiplictyTask"); return; }
  //cout << "MULT loaded \n \n" << endl;
  

  // //==========================================
  // // PID QA Task
  // //==========================================
 
  // gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  // AliAnalysisTask *PIDQATask = AddTaskPIDqa();
     
  //==========================================
  // create task
  //==========================================

  gROOT->LoadMacro("AliAnalysisCorrelationEventCollection.cxx++g");
  gROOT->LoadMacro("AliAnalysisTaskCorrelationhCasc.cxx++g");
  gROOT->LoadMacro("AddTaskCorrelationhCasc.C");
  

  Bool_t isLocal=kFALSE;
  if (runtype == "local") isLocal = kTRUE;

  AliAnalysisTaskCorrelationhCasc *task = AddTaskCorrelationhCasc("name",3,30, isLocal, isMC, isEff,isHybridMCTr, 10, EtaTrigger, EtahAssoc,EtaV0Assoc, FilterBitValue, year, "Xi");

  // enable debug printouts
  mgr->SetDebugLevel(2);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  // start analysis
  Printf("Starting Analysis....");
  mgr->StartAnalysis(runtype,nentries,firstentry);

}

AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const char *proofcluster, const char *proofdataset)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(gridmode);
  
 
   // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-30-alice10-21");
  plugin->SetAliROOTVersion("v5-09-47c-1");
  plugin->SetAliPhysicsVersion("vAN-20190630-1");
   

  TString includes_str = "-Wno-deprecated -I$. -I$CGAL_DIR/include -I$FASTJET/include -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include";
  plugin->AddIncludePath(includes_str.Data()); // for grid running

  //con questa versione non funziona
  // plugin->SetAPIVersion("V1.1x");
  // plugin->SetROOTVersion("v5-34-30-alice-8");
  // plugin->SetAliROOTVersion("v5-09-47b-1");
  // plugin->SetAliPhysicsVersion("vAN-20190429-1");
    
  // are LOADED AUTOMATICALLY!
  // plugin->SetAliPhysicsVersion("vAN-20190429-1");
    
  cout<<"Q2"<<endl;

  // Declare input data to be processed.

  //for data    (16k)
  //    plugin->SetGridDataDir("/alice/data/2016/LHC16k/"); 
  //    plugin->SetDataPattern("/pass2/AOD208/*/AliAOD.root"); 
  //
  //for data (15f)
  //plugin->SetGridDataDir("/alice/data/2015/LHC15f/"); 
  //plugin->SetDataPattern("/pass2/AOD208/*/AliAOD.root"); 


  //for MC 18f1_extra
    plugin->SetGridDataDir("/alice/sim/2018/LHC18f1_extra"); 
   plugin->SetDataPattern("/AOD209/*/AliAOD.root"); 

  //for MC 15g3b1
   //  plugin->SetGridDataDir("/alice/sim/2015/LHC15g3b1"); 
   //  plugin->SetDataPattern("/AOD/*/AliAOD.root"); 


  plugin->SetRunPrefix(""); //put 000 if data

  //for LHC15f pass2 (the ones Fiorella used for her analysis)
  
  //  Int_t nrun[56] = {225000, 225011, 225016, 225026, 225031, 225035, 225037, 225041, 225043, 225050, 225051, 225052, 225106, 225305, 225307, 225309, 225310, 225313, 225314, 225315, 225322, 225576, 225578, 225579, 225580, 225582, 225586, 225587, 225705, 225707, 225708, 225709, 225710, 225716, 225717, 225719, 225753, 225757, 225762, 225763, 225766, 225768, 226062, 226170, 226220, 226225, 226444, 226445, 226452, 226466, 226468, 226472, 226476, 226483, 226495, 226500};
  //  Int_t nrun[51]={226500, 226495, 226483, 226476, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225705, 225587, 225586, 225579, 225578, 225576, 225322, 225315, 225314, 225313, 225310, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026};


  //utilizzati per 2016k
    Int_t nrun[161]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941};
  //  Int_t nrun[131]= {258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941};
  //  Int_t nrun[30]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202};
  //    Int_t nrun[15]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306};
  //  Int_t nrun[61]= { 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941 };
  //  Int_t nrun[50]= {258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936};
  //  Int_t nrun[25]= {258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063};
  //  Innt_t nrun[96] = { 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941};
  //Int_t nrun[10]={258537, 258499, 258477, 258456, 258454,  258452, 258426, 258393, 258391, 258387};

  
  for(Int_t i = 0; i<161;i++){
    plugin->AddRunNumber(nrun[i]);
  }

  plugin->SetNrunsPerMaster(5);
  plugin->SetOutputToRunNo();

  // plugin->SetUseSubmitPolicy(kTRUE);
  // comment out the next line when using the "terminate" option, unless
  // you want separate merged files for each run
  //   plugin->SetMergeViaJDL();

  // Method 2: Declare existing data files (raw collections, xml collections, root file)
  // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
  // XML collections added via this method can be combined with the first method if
  // the content is compatible (using or not tags)
  //   plugin->AddDataFile("tag.xml");
  //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(taskname);

  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/taskname/out
    
  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  plugin->SetAnalysisSource("AliAnalysisCorrelationEventCollection.cxx  AliAnalysisTaskCorrelationhCasc.cxx");

  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.

  //   plugin->SetAdditionalLibs("AliAnalysisTaskOverProton.h AliAnalysisTaskOverProton.cxx libGui.so libProof.so libMinuit.so libXMLParser.so libSTAT.so libCDB.so libRAWDatabase.so libRAWDatarec.so libSTEER.so libANALYSIS.so libANALYSISalice.so libANALYSIScalib.so libTENDER.so libCORRFW.so libPWGUDbase.so libTPCbase.so libTPCrec.so libTPCcalib.so libTRDbase.so libTRDrec.so libITSbase.so libITSrec.so libHMPIDbase.so libPWGPP.so libCore.so libPhysics.so libGeom.so libVMC.so libNet.so libTree.so libSTEERBase.so libESD.so libTOFbase.so libTOFsim.so libTOFcalib.so libTOFrec.so libOADB.so libJETAN.so libThread.so libVZERObase.so libVZEROrec.so libEMCALraw.so libEMCALUtils.so libEMCALbase.so libEMCALrec.so libT0base.so libT0rec.so");

  plugin->SetAdditionalLibs("libOADB.so AliAnalysisCorrelationEventCollection.h AliAnalysisCorrelationEventCollection.cxx AliAnalysisTaskCorrelationhCasc.h AliAnalysisTaskCorrelationhCasc.cxx");

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("AnalysisResultsLocalhCasc.root");
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s_2d.C",taskname));

  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(100);

  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s_2d.sh",taskname));

  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(1);

  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);

  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);

  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");

  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName(Form("%s_2d.jdl",taskname));

  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      

  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
    
  plugin->SetFastReadOption();
  cout<<"Q3"<<endl;

  //----------------------------------------------------------
  //---      PROOF MODE SPECIFIC SETTINGS         ------------
  //---------------------------------------------------------- 
  // Proof cluster
  plugin->SetProofCluster(proofcluster);
  // Dataset to be used   
  plugin->SetProofDataSet(proofdataset);
  // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
  plugin->SetProofReset(0);
  // May limit number of workers
  plugin->SetNproofWorkers(0);
  // May limit the number of workers per slave
  plugin->SetNproofWorkersPerSlave(1);    
  // May use a specific version of root installed in proof
  plugin->SetRootVersionForProof("current");
  // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
  plugin->SetAliRootMode("default"); // Loads AF libs by default
  // May request ClearPackages (individual ClearPackage not supported)
  plugin->SetClearPackages(kFALSE);
  // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
  plugin->SetFileForTestMode("fileMC.txt"); // file should contain path name to a local directory containg *ESDs.root etc
  // Request connection to alien upon connection to grid
  plugin->SetProofConnectGrid(kFALSE);

  return plugin;
}
  

