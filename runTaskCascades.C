class AliAnalysisGrid;

//______________________________________________________________________________
void runTaskCascades(	     //local should be matched with test
	     const char* corrtype = "hh", //hh for hadron-hadron correlation
	     const char* runtype = "grid", // local, proof or grid
	     const char *gridmode = "full", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
	     const bool isMC = 0, // 1 = MC, 0 = DATA 
	     const bool isEff=0, //1 = MC per efficienza, 0 = MC truth
	     const Float_t EtaTrigger=0.8,
	     const Float_t EtahAssoc=0.8,
	     const Long64_t nentries = 1234567890,//2000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
	     const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
	     const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
	     const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
	     const char *taskname = "K0sjetpp_LHC18f1_extra_Cascades_10run_Bis"   // sets name of grid generated macros
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

    const bool bMCtruth = isMC;             // 1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
  const bool bMCphyssel = isMC;           // 1 = looking at MC truth or reconstructed, 0 = looking at real data
  const bool isPIDResponseMCtruth = isMC; // 1 = MC, 0 = DATA 
   
  /* 
  TString file;
  if (isMC==0) file = "fileMC.txt";
  else file = "file.txt";
  cout<<"1"<<endl;
  */
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
  //// if(runtype!="local"){
  AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset); 
  mgr->SetGridHandler(plugin);
  //  }
  // AliVEventHandler* esdH = new AliESDInputHandler();
  // mgr->SetInputEventHandler(esdH);

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
  
  cout<<"4"<<endl;
  
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
  AliAnalysisTask *PIDTask = AddTaskPIDResponse(kFALSE); //dipende da Dati o MC
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

  gROOT->LoadMacro("AliAnalysisTaskCascades.cxx++g");
  gROOT->LoadMacro("AddMyTaskCascades.C");
  
  cout<<"------------------------------------------------------------------- QUI"<<endl;

  // AliAnalysisTaskSE *task = AddTaskNucleiv2SP("Deuteronv2SP",1,kFALSE, 10, "CL1");

  // gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Hypernuclei/AliAnalysisTaskNucleiv2SP.cxx");

  Bool_t ishhCorr=kFALSE;
  if (corrtype == "hh") ishhCorr = kTRUE;
  Bool_t isLocal=kFALSE;
  if (runtype == "local") isLocal = kTRUE;

  AliAnalysisTaskCascades *task = AddMyTaskCascades("name",3,30,ishhCorr, isLocal, isMC, isEff, 5, EtaTrigger, EtahAssoc);
  // mgr->AddTask(task,2,10);
  // mgr->AddTask(task);
   
  cout << "  5" << endl;
  // enable debug printouts
  mgr->SetDebugLevel(2);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  cout << " end of debug printouts \n \n" << endl;
  // start analysis
  Printf("Starting Analysis....");
  mgr->StartAnalysis(runtype,nentries,firstentry);

  //______________________________________________________________________________
  //  Printf("ciao Chiara");
  /*
  if(runtype == "local"){
    cout << "Running locally \n \n \n" << endl;
    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    // add a few files to the chain (change this so that your local files are added
    TString filetorunon="Data/AliAOD.root";
    cout << "running on "<< filetorunon << endl;
    chain->Add(filetorunon);
    //cout << "qua " << endl;
    //	 chain->Add("Data/AliAOD_MC.root");
    // start the analysis locally, reading the events from the tchain
    mgr->StartAnalysis("local", chain);
    //cout<<  " oi " << endl;
   
  } else {
    cout << "Not Running locally \n \n \n" << endl;
    mgr->StartAnalysis(runtype,nentries,firstentry);
    }*/
}

AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const char *proofcluster, const char *proofdataset)
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(gridmode);
  

  cout<<"Q1"<<endl;

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

  // Method 1: Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  //    plugin->SetGridDataDir("/alice/data/2011/LHC11h_2"); //LHC11h
    
  // plugin->SetGridDataDir("/alice/data/2015/LHC15f");
  // plugin->SetDataPattern("pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
  // plugin->SetRunPrefix("000");   // real data

  // ...then add run numbers to be considered
  // to be added

  // Pb-Pb
  // plugin->SetGridDataDir("/alice/data/2011/LHC11h_2/"); //sim
  // plugin->SetDataPattern("ESDs/pass2/AOD145/*AliAOD.root"); // sim
  // plugin->SetRunPrefix("000"); 

  // plugin->AddRunNumber(170593);     
  // plugin->SetNrunsPerMaster(1);
  // plugin->SetOutputToRunNo();

  // // pp
    
  // //LHC10b
    
  // plugin->SetGridDataDir("/alice/data/2010/LHC10b/"); //sim
  // plugin->SetDataPattern("ESDs/pass2/AOD137/*AliAOD.root"); // sim
  // plugin->SetRunPrefix("000"); 
    
  // Int_t nrun[22]={117222, 117220, 117116, 117112, 117099, 117092, 117063, 117060, 117059, 117053, 117052, 117050, 117048, 116645, 116643, 116574, 116571, 116562, 116403, 116402, 116288, 116102};
  // for(Int_t i = 0; i<1;i++){
  //   //    for(Int_t i = 0; i<1;i++){
  //   plugin->AddRunNumber(nrun[i]);
  // }
       
  // LHC16q

  //for data    
  //    plugin->SetGridDataDir("/alice/data/2016/LHC16k/"); 
  //    plugin->SetDataPattern("/pass2/AOD208/*/AliAOD.root"); 
  // //
  //for MC
  //     plugin->SetGridDataDir("/alice/sim/2018/LHC18j4_extra"); 
     plugin->SetGridDataDir("/alice/sim/2018/LHC18f1_extra"); 
  //  plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a2_extra/"); 
     plugin->SetDataPattern("/AOD209/*/AliAOD.root"); 

//  plugin->SetDataPattern("/pass1/AOD/*/AliAOD.root"); 
  plugin->SetRunPrefix(""); //put 000 if data

  //2017m
  //  Int_t nrun[108]={280140, 280135, 280134, 280131, 280126, 280118, 280114, 280111, 280108, 280107, 280066, 280052, 280051, 279879, 279855, 279854, 279853, 279830, 279827, 279826, 279773, 279749, 279747, 279719, 279718, 279715, 279689, 279688, 279687, 279684, 279683, 279682, 279679, 279677, 279676, 279642, 279641, 279632, 279630, 279559, 279550, 279491, 279488, 279487, 279483, 279441, 279439, 279435, 279410, 279391, 279355, 279354, 279349, 279348, 279344, 279342, 279312, 279310, 279309, 279274, 279273, 279270, 279268, 279267, 279265, 279264, 279242, 279238, 279235, 279234, 279232, 279208, 279207, 279201, 279199, 279157, 279155, 279130, 279123, 279122, 279118, 279117, 279107, 279106, 279075, 279074, 279073, 279069, 279068, 279044, 279043, 279041, 279036, 279035, 279008, 279007, 279005, 279000, 278999, 278964, 278963, 278960, 278959, 278941, 278939, 278936, 278915, 278914};


  //2018b data
  //  Int_t nrun[25]={  285396, 285365, 285364, 285347, 285328, 285327, 285224, 285222, 285203, 285202, 285200, 285165, 285127, 285125, 285108, 285106, 285066, 285065, 285064, 285015, 285014, 285013, 285012, 285011, 285009};

  //2018d8_extra MC (anchored to 2016l)
  //  Int_t nrun[79]={  258919, 258920, 258921, 258923, 258926, 258931, 258962, 258964, 259086, 259088, 259090, 259091, 259096, 259099, 259117, 259118, 259162, 259164, 259204, 259257, 259261, 259263, 259264, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259303, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259473, 259477, 259649, 259650, 259668, 259697, 259700, 259703, 259704, 259705, 259711, 259713, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259789, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888, 260010, 260011, 260014};


  //2016l
  //   Int_t nrun[58]={259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962};
  // Int_t nrun[1]={259888};
  // Int_t nrun[5]={258919, 258920, 258921, 258923, 258926};
  //2017h
  //Int_t nrun[88]={273103, 273100, 273099, 273077, 273010, 273009, 272985, 272983, 272976, 272949, 272947, 272939, 272935, 272934, 272933, 272932, 272905, 272903, 272880, 272873, 272871, 272870, 272836, 272834, 272833, 272829, 272828, 272784, 272783, 272782, 272764, 272763, 272760, 272749, 272747, 272712, 272691, 272690, 272620, 272610, 272608, 272607, 272585, 272577, 272575, 272574, 272521, 272468, 272466, 272463, 272462, 272461, 272413, 272411, 272400, 272399, 272395, 272394, 272389, 272388, 272360, 272359, 272340, 272335, 272194, 272156, 272155, 272154, 272153, 272152, 272151, 272123, 272101, 272100, 272076, 272042, 272040, 272039, 272038, 272036,   272020, 272018, 271886, 271880, 271874, 271873, 271871, 271870};

  //2017k
  // Int_t nrun[105]={276508, 276507, 276506, 276462, 276439, 276438, 276437, 276435, 276351, 276348, 276302, 276297, 276294, 276292, 276290, 276259, 276257, 276230, 276205, 276178, 276177, 276170, 276169, 276166, 276145, 276140, 276135, 276104, 276102, 276099, 276098, 276097, 275847, 275664, 275661, 275650, 275648, 275647, 275624, 275623, 275622, 275621, 275617, 275612, 275559, 275558, 275515, 275472, 275471, 275467, 275459, 275457, 275456, 275453, 275452, 275448, 275443, 275406, 275404, 275401, 275372, 275369, 275361, 275360, 275333, 275332, 275328, 275326, 275324, 275322, 275314, 275283, 275247, 275246, 275245, 275239, 275188, 275184, 275180, 275177, 275174, 275173, 275151, 275150, 275149, 275076, 275075, 275073, 275068, 275067, 274979, 274978, 274886, 274882, 274878, 274877, 274822, 274821, 274815, 274806, 274803, 274802, 274801, 274708, 274690};


  //2018m
  //  Int_t nrun[242]= {292839, 292836, 292834, 292832, 292831, 292811, 292810, 292809, 292804, 292803, 292752, 292750, 292748, 292747, 292744, 292739, 292737, 292704, 292701, 292698, 292696, 292695, 292693, 292586, 292584, 292563, 292560, 292559, 292557, 292554, 292553, 292526, 292524, 292523, 292521, 292500, 292497, 292496, 292495, 292461, 292460, 292457, 292456, 292434, 292432, 292430, 292429, 292428, 292406, 292405, 292398, 292397, 292298, 292273, 292265, 292242, 292241, 292240, 292218, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114, 292109, 292108, 292107, 292106, 292081, 292080, 292077, 292075, 292067, 292062, 292061, 292060, 292040, 292012, 291982, 291977, 291976, 291953, 291948, 291946, 291945, 291944, 291943, 291942, 291803, 291796, 291795, 291769, 291768, 291766, 291762, 291760, 291756, 291755, 291729, 291706, 291698, 291697, 291690, 291665, 291661, 291657, 291626, 291624, 291622, 291618, 291615, 291614, 291590, 291485, 291484, 291482, 291481, 291457, 291456, 291453, 291451, 291447, 291424, 291420, 291417, 291416, 291402, 291400, 291399, 291397, 291377, 291375, 291363, 291362, 291361, 291360, 291286, 291285, 291284, 291282, 291266, 291265, 291263, 291262, 291257, 291240, 291209, 291188, 291143, 291116, 291111, 291110, 291101, 291100, 291093, 291069, 291066, 291065, 291041, 291037, 291035, 291006, 291005, 291004, 291003, 291002, 290980, 290979, 290976, 290975, 290974, 290948, 290944, 290943, 290941, 290935, 290932, 290895, 290894, 290888, 290887, 290886, 290862, 290860, 290853, 290848, 290846, 290843, 290841, 290790, 290787, 290766, 290689, 290687, 290665, 290660, 290645, 290632, 290627, 290615, 290614, 290613, 290612, 290590, 290588, 290553, 290550, 290549, 290544, 290540, 290539, 290538, 290501, 290500, 290499, 290469, 290467, 290459, 290458, 290456, 290427, 290426, 290425, 290423, 290412, 290411, 290404, 290401, 290399, 290376, 290375, 290374, 290350, 290327, 290323};
  //Int_t nrun[118]= {292839, 292836, 292834, 292832, 292831, 292811, 292810, 292809, 292804, 292803, 292752, 292750, 292748, 292747, 292744, 292739, 292737, 292704, 292701, 292698, 292696, 292695, 292693, 292586, 292584, 292563, 292560, 292559, 292557, 292554, 292553, 292526, 292524, 292523, 292521, 292500, 292497, 292496, 292495, 292461, 292460, 292457, 292456, 292434, 292432, 292430, 292429, 292428, 292406, 292405, 292398, 292397, 292298, 292273, 292265, 292242, 292241, 292240, 292218, 292192, 292168, 292167, 292166, 292164, 292163, 292162, 292161, 292160, 292140, 292115, 292114, 292109, 292108, 292107, 292106, 292081, 292080, 292077, 292075, 292067, 292062, 292061, 292060, 292040, 292012, 291982, 291977, 291976, 291953, 291948, 291946, 291945, 291944, 291943, 291942, 291803, 291796, 291795, 291769, 291768, 291766, 291762, 291760, 291756, 291755, 291729, 291706, 291698, 291697, 291690, 291665, 291661, 291657, 291626, 291624, 291622, 291618, 291615};


//utilizzati per MC 2018d8
  ////utilizzati per 2016l_fifth e anche per MC associato (2018d8)
  //Int_t nrun[15]={259888,  259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788,  259781, 259756, 259752, 259751, 259750}; 
  // Int_t nrun[20]={ 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303};
  // Int_t nrun[58]={259888,  259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962};
//Int_t nrun[23]={ 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962};

//utilizzati per 2016k
  //Int_t nrun[161]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941};
//Int_t nrun[15]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306};
//  Int_t nrun[61]= { 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257012, 257011, 256944, 256942, 256941 };
 
  Int_t nrun[10]={258537, 258499, 258477, 258456, 258454,  258452, 258426, 258393, 258391, 258387};


  //2016o
  //Int_t nrun[71]= {264035, 264033, 263985, 263984, 263981, 263978, 263977, 263923, 263920, 263917, 263916, 263905, 263866, 263863, 263810, 263803, 263793, 263792, 263790, 263787, 263786, 263785, 263784, 263744, 263743, 263741, 263739, 263738, 263737, 263691, 263690, 263682, 263663, 263662, 263657, 263654, 263652, 263647, 263529, 263497, 263496, 263490, 263487, 263332, 263331, 262858, 262855, 262853, 262849, 262847, 262844, 262842, 262841, 262778, 262777, 262776, 262768, 262760, 262727, 262725, 262723, 262719, 262717, 262713, 262708, 262706, 262705, 262428, 262426, 262425, 262424 };

  //utilizzati per MC associato (2017d20a2_extra)
  //  Int_t nrun[1]={  258919};
 //Int_t nrun[74]={  258919, 258920, 258921, 258923, 258962, 258964, 259086, 259088, 259090, 259091, 259096, 259099, 259117, 259118, 259162, 259164, 259204, 259257, 259261, 259263, 259264, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259303, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259473, 259477, 259649, 259650, 259668, 259697, 259700, 259703, 259704, 259705, 259711, 259713, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259789, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888}
 
   //utilizzato per tentativo //Int_t nrun[1]={259888};

  //  Int_t nrun[80]={260014, 260011, 260010, 259979, 259961, 259954, 259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259792, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259713, 259711, 259705, 259704, 259703, 259700, 259697, 259668, 259650, 259649, 259477, 259473, 259471, 259470, 259469, 259396, 259395, 259394, 259389, 259388, 259382, 259381, 259379, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 259086, 258964, 258962}


  //per LHC2016l secondo giro
  // Int_t nrun[47]={260014, 260011, 260010, 259979, 259961, 259954, 259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259792, 259789, 259788, 259781, 259713, 259711, 259705, 259704, 259703, 259700, 259697, 259668, 259650, 259649, 259477, 259473, 259471, 259470, 259469, 259396, 259395, 259394, 259389, 259388, 259382, 259381, 259379, 259378, 259342, 259341, 259340, 259339, 259336}
  
  //Utilizzati per LHC2016l primo giro
   // Int_t nrun[27]={259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 259086, 258964, 258962};
  // Utilizzati per LHC2017o
  // Int_t nrun[22]={281961, 281956, 281953, 281946, 281940, 281939, 281932, 281931, 281928, 281920, 281919, 281918, 281916, 281915, 281895, 281894, 281893, 281892, 281756, 281755, 281754, 281753};
  //  for(Int_t i = 0; i<31;i++){
  
  for(Int_t i = 0; i<10;i++){
    plugin->AddRunNumber(nrun[i]);
  }
  //   plugin->AddRunNumber(265521)
  
  // //MC LHC17d
  // plugin->SetGridDataDir("/alice/sim/2017/LHC17f2a_fast_fix/"); //sim
  // plugin->SetDataPattern("/AOD202/*AliAOD.root"); // sim
  // // plugin->AddRunNumber(265309);
  // Int_t nrun[36]={265309, 265332, 265334, 265335, 265336, 265338, 265339, 265342, 265343, 265344, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525, 267163, 267164, 267165, 267166}
  // for(Int_t i = 0; i<36;i++){
  //   //  for(Int_t i = 0; i<1;i++){
  //   plugin->AddRunNumber(nrun[i]);
  // }

  //    plugin->SetNrunsPerMaster(1);
  //    plugin->SetNrunsPerMaster(5);
  plugin->SetNrunsPerMaster(5);
  //plugin->SetNrunsPerMaster(1);
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
  plugin->SetAnalysisSource(" AliAnalysisTaskCascades.cxx");

  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.

  //   plugin->SetAdditionalLibs("AliAnalysisTaskOverProton.h AliAnalysisTaskOverProton.cxx libGui.so libProof.so libMinuit.so libXMLParser.so libSTAT.so libCDB.so libRAWDatabase.so libRAWDatarec.so libSTEER.so libANALYSIS.so libANALYSISalice.so libANALYSIScalib.so libTENDER.so libCORRFW.so libPWGUDbase.so libTPCbase.so libTPCrec.so libTPCcalib.so libTRDbase.so libTRDrec.so libITSbase.so libITSrec.so libHMPIDbase.so libPWGPP.so libCore.so libPhysics.so libGeom.so libVMC.so libNet.so libTree.so libSTEERBase.so libESD.so libTOFbase.so libTOFsim.so libTOFcalib.so libTOFrec.so libOADB.so libJETAN.so libThread.so libVZERObase.so libVZEROrec.so libEMCALraw.so libEMCALUtils.so libEMCALbase.so libEMCALrec.so libT0base.so libT0rec.so");

   plugin->SetAdditionalLibs("libOADB.so AliAnalysisTaskCascades.h AliAnalysisTaskCascades.cxx");

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("AnalysisResultsLocal.root");
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
  

