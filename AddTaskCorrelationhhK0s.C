//used to instantiate an object of the task,define input and output and connect it to manager

AliAnalysisTaskCorrelationhhK0s* AddTaskCorrelationhhK0s(TString name = "name", Float_t minpt=3, Float_t maxpt=15, Bool_t ishhCorr,bool isLocal=kTRUE, bool isMC=kTRUE, bool isEff=kTRUE, Int_t EvtToMix=50, Float_t EtaTrigger=0.8, Float_t EtahAssoc=0.8,Float_t EtaV0Assoc=0.8, Int_t FilterBitValue=128, Int_t year=2010){

  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the fun2ction
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler, again via a static method. 
  // this handler is part of the managing system and feeds events
  // to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  // by default, a file is open for writing. here, we get the filename
  TString fileName;  
  if(isLocal) {
    fileName = "AnalysisResultsLocal.root";
  }
  else {
    fileName = AliAnalysisManager::GetCommonFileName();
  }
  fileName += ":MyTask";      // create a subfolder in the file

  // now we create an instance of your task
  AliAnalysisTaskCorrelationhhK0s* task = new AliAnalysisTaskCorrelationhhK0s(name.Data());   
  if(!task) return 0x0;
  //task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  // add your task to the manager
   task->SetMinPt(minpt); //setto i valori di min pt, modificando i valori di default, settati nel momento in cui viene creato l'oggetto task
   task->SetMaxPt(maxpt);
   task->SetCorr(ishhCorr);
   task->SetMC(isMC);
   task->SetEff(isEff);
   task->SetEvtToMix(EvtToMix);
   task->SetEtaTrigger(EtaTrigger);
   task->SetEtahAssoc(EtahAssoc);
   task->SetEtaV0Assoc(EtaV0Assoc);
   task->SetFilterBit(FilterBitValue);
   task->SetYear(year);

  mgr->AddTask(task);
  // your task needs input: here we connect the manager to your task
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer("MyOutputContainer1", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,3,mgr->CreateContainer("MyOutputContainer2", TTree::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,4,mgr->CreateContainer("MyOutputContainer3", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,5,mgr->CreateContainer("MyOutputContainer4", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  mgr->ConnectOutput(task,6,mgr->CreateContainer("Risoluzione", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
  // in the end, this macro returns a pointer to your task. this will be convenient later on
  // when you will run your analysis in an analysis train on grid
  return task;
}
