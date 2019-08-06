/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 ************************************v**************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TMath.h"
#include "Riostream.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskMyTask.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
//#include "AliMultSelectionBase.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
//#include "AliAnalysisKPEventCollectionChiara.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliVVertex.h"


class AliAnalysisTaskMyTask;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMyTask) // classimp: necessary for root

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() :AliAnalysisTaskSE(), 
  fAnalysisType("AOD"), 
  fCollidingSystem("pp"), 
  fAOD(0), 
  fPIDResponse(0),
//fMultSelection(0),
  fEventCuts(0), 			  			
  fOutputList(0), 
  fSignalTree(0), 
  fBkgTree(0), 
  fMCEvent(0), 
  fReadMCTruth(0),
// fEventColl(0x0), 
//fEvt(0x0), 
  fzVertexBins(20), 
  fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(50),
  fnEventsToMix(50),
  fHistPt(0), 
  fHistPtvsMult(0), 
  fHistZvertex(0),  
  fHist_eta_phi(0),  
  fHist_multiplicity(0),
						fHistEventMult(0), 
						fHistEventV0(0), 
						fHistTrack(0), 
  fHistPDG(0), 
  fMassV0(0), 
  fHistMultvsV0(0),
  fHistMassvsPt(0),
  fHistMassvsPt_tagli(0x0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistPtArmvsAlpha(0),
  fHistPtArmvsAlphaAfterSelection(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistMultiplicityOfMixedEvent(0),
  fHistGeneratedTriggerPtPhi(0),
  fHistSelectedTriggerPtPhi(0),
  fHistGeneratedV0PtPhi(0),
  fHistSelectedV0PtPhi(0),
  fHistGeneratedTriggerPtEta(0),
  fHistSelectedTriggerPtEta(0),
  fHistGeneratedV0PtEta(0),
  fHistSelectedV0PtEta(0),
  fHistReconstructedV0PtMass(0),
  fHistSelectedV0PtMass(0),
  fHistResolutionTriggerPt(0),
  fHistResolutionTriggerPhi(0),
  fHistResolutionTriggerEta(0),
  fHistResolutionV0Pt(0),
  fHistResolutionV0Phi(0),
  fHistResolutionV0Eta(0),
  fHistPrimaryTrigger(0),
  fHistPrimaryV0(0),
  fminPtj(2), 
  fmaxPtj(10), 
  fV0("kK0s"),  
  fminPtV0(0), 
  fmaxPtV0(30),  
  Evcounter(0), 
  Evcounterczero(0),
  fmolt(5),
  farrGT(0), 
  fTrackBufferSize(20200),
  fTreeVariablePtTrigger(0),		      
  fTreeVariableChargeTrigger(0),		      
  fTreeVariableEtaTrigger(0), 		      
  fTreeVariablePhiTrigger(0),		      
  fTreeVariableDCAz(0),			      
  fTreeVariableDCAxy(0),			      
  fTreeVariableisPrimaryTrigger(0),  
  fTreeVariableisPrimaryV0(0),  
  fTreeVariableRapK0Short(0),		      	      
  fTreeVariableDcaV0ToPrimVertex (0),	      	      
  fTreeVariableDcaPosToPrimVertex(0),	      	      
  fTreeVariableDcaNegToPrimVertex(0),	      	      
  fTreeVariableV0CosineOfPointingAngle(0),      	      
  fTreeVariablePtV0(0),			      
  fTreeVariableInvMassK0s(0),		      
  fTreeVariableInvMassLambda(0),		      
  fTreeVariableInvMassAntiLambda(0),		      
  fTreeVariableEtaV0(0),			      
  fTreeVariablePhiV0(0),			      
  fTreeVariablePtArmenteros(0),                   
  fTreeVariableAlpha(0),
  fTreeVariableDeltaEta(0),			       
  fTreeVariableDeltaPhi(0),
  fTreeVariableDeltaTheta(0),			       
  fTreeVariableMultiplicity(0),                   
  fTreeVariableZvertex(0),
  fTreeVariablePDGCode(0),
  FifoShiftok(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purpos, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
								 fAnalysisType("AOD"), 
								 fCollidingSystem("pp"), 
								 fAOD(0), 
								 fPIDResponse(0),
								 //fMultSelection(0), 			  			
								 fEventCuts(0),								 
								 fOutputList(0), 
								 fSignalTree(0), 
								 fBkgTree(0), 
								 fMCEvent(0), 
								 fReadMCTruth(0),
								 //fEventColl(0x0), 
								 //fEvt(0x0), 
								 fzVertexBins(20), 
								 fnMultBins(20),	 
								 fMaxFirstMult(50),
								 fMaxSecondMult(50),
								 fnEventsToMix(50),
  fHistPt(0), 
  fHistPtvsMult(0), 
  fHistZvertex(0),  
  fHist_eta_phi(0),  
  fHist_multiplicity(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistTrack(0), 
  fHistPDG(0), 
  fMassV0(0), 
  fHistMultvsV0(0),
  fHistMassvsPt(0),
  fHistMassvsPt_tagli(0x0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistPtArmvsAlpha(0),
  fHistPtArmvsAlphaAfterSelection(0),
  fHistTrigger(0),
  fHistTriggerMCTructh(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistMultiplicityOfMixedEvent(0),
  fHistGeneratedTriggerPtPhi(0),
  fHistSelectedTriggerPtPhi(0),
  fHistGeneratedV0PtPhi(0),
  fHistSelectedV0PtPhi(0),
  fHistGeneratedTriggerPtEta(0),
  fHistSelectedTriggerPtEta(0),
  fHistGeneratedV0PtEta(0),
  fHistSelectedV0PtEta(0),
  fHistReconstructedV0PtMass(0),
  fHistSelectedV0PtMass(0),
  fHistResolutionTriggerPt(0),
  fHistResolutionTriggerPhi(0),
  fHistResolutionTriggerEta(0),
  fHistResolutionV0Pt(0),
  fHistResolutionV0Phi(0),
  fHistResolutionV0Eta(0),
  fHistPrimaryTrigger(0),
  fHistPrimaryV0(0),
  fminPtj(2), 
  fmaxPtj(10), 
  fV0("kK0s"),  
  fminPtV0(0), 
  fmaxPtV0(30),  
  Evcounter(0), 
  Evcounterczero(0),
  fmolt(5),
  farrGT(0), 
  fTrackBufferSize(20200),
  fTreeVariablePtTrigger(0),		      
  fTreeVariableChargeTrigger(0),		      
  fTreeVariableEtaTrigger(0), 		      
  fTreeVariablePhiTrigger(0),		      
  fTreeVariableDCAz(0),			      
  fTreeVariableDCAxy(0),			      
  fTreeVariableisPrimaryTrigger(0),
  fTreeVariableisPrimaryV0(0),
  fTreeVariableRapK0Short(0),		      	      
  fTreeVariableDcaV0ToPrimVertex (0),	      	      
  fTreeVariableDcaPosToPrimVertex(0),	      	      
  fTreeVariableDcaNegToPrimVertex(0),	      	      
  fTreeVariableV0CosineOfPointingAngle(0),      	      
  fTreeVariablePtV0(0),			      
  fTreeVariableInvMassK0s(0),		      
  fTreeVariableInvMassLambda(0),		      
  fTreeVariableInvMassAntiLambda(0),		      
  fTreeVariableEtaV0(0),			      
  fTreeVariablePhiV0(0),			      
  fTreeVariablePtArmenteros(0),                   
  fTreeVariableAlpha(0),
  fTreeVariableDeltaEta(0),			       
  fTreeVariableDeltaPhi(0),
  fTreeVariableDeltaTheta(0),			       
  fTreeVariableMultiplicity(0),                   
  fTreeVariableZvertex(0),
  fTreeVariablePDGCode(0),
  FifoShiftok(kFALSE)
{
                    
  
  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it, 
  // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
  // you can add more output objects by calling DefineOutput(2, classname::Class())
  // if you add more output objects, make sure to call PostData for all of them, and to
  // make changes to your AddTask macro!

  DefineOutput(2, TTree::Class());  
  DefineOutput(3, TTree::Class());  
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
  if(fSignalTree) {
    delete fSignalTree;     // at the end of your task, it is deleted from memory by calling this function
  }
  if(fBkgTree) {
    delete fBkgTree;     // at the end of your task, it is deleted from memory by calling this function
  }
  if (farrGT)
    delete[] farrGT;
  farrGT=0;

  // if (fHistMassvsPt)
  // delete[] fHistMassvsPt;

  // if (fHistMassvsPt_tagli)
  // delete[] fHistMassvsPt_tagli ;

  /*
  for(unsigned short i=0; i < fzVertexBins; i++){
    for(unsigned short j=0; j < fnMultBins; j++){
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }
  delete[] fEventColl;
  */
}
  

void AliAnalysisTaskMyTask::ProcessMCParticles(Bool_t Generated, AliAODTrack *track, Int_t& labelPrimOrSec, Float_t lPercentiles, Bool_t isV0)
{
  // process MC particles
  //  TList *list=fAOD->GetList();
  Float_t moltep[6]={0,5,10,30,50,100};  //valori associati a centralita'  
  TClonesArray* AODMCTrackArraybis =0x0;  
  AODMCTrackArraybis = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArraybis == NULL){
    return;
    Printf("ERROR: stack not available");
  }
  
  if(Generated){
    // Loop over all primary MC particle
    for(Long_t i = 0; i < AODMCTrackArraybis->GetEntriesFast(); i++) {
      
      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(i));
      if (!particle) continue;
      
      if (isV0==kFALSE) fHistPDG->Fill(particle->GetPdgCode());      

      if((particle->Charge())!=0){
	if(TMath::Abs(particle->Eta())>0.8)continue; //I need to select particles within this eta range!
	if (!(particle->IsPhysicalPrimary()))continue; 
	fHistGeneratedTriggerPtPhi->Fill(particle->Pt(), particle->Phi(), lPercentiles);
	fHistGeneratedTriggerPtEta->Fill(particle->Pt(), particle->Eta(), lPercentiles);
      }
      if(isV0==kTRUE){ 
	if ((particle->GetPdgCode())!=310) continue;
	//cout << "K0s mass two methods " << particle->M() << "  " << particle->GetCalcMass()<< endl;
	if(TMath::Abs(particle->Eta())>0.8)continue;
	if (!(particle->IsPhysicalPrimary()))continue;
	cout << "these are all K0s generated passing selection criteria: label K0s " << particle->Label()<<endl; 
	fHistGeneratedV0PtPhi->Fill(particle->Pt(),particle->Phi(), lPercentiles );
	fHistGeneratedV0PtEta->Fill(particle->Pt(),particle->Eta(), lPercentiles );
      }
      
    }
  }
  else {
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(TMath::Abs(track->GetLabel())));
    // AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(label));
    //   if(particle->IsPhysicalPrimary() && TMath::Abs(particle->Eta())<=0.8 && particle->Charge()!=0){
   if(particle->IsPhysicalPrimary()){
      cout <<"selected trigger Pt " <<  track->Pt() << endl;
      fHistResolutionTriggerPt->Fill(track->Pt()- particle->Pt(), lPercentiles);
      fHistResolutionTriggerPhi->Fill(track->Phi()- particle->Phi(), lPercentiles);
      fHistResolutionTriggerEta->Fill(track->Eta()- particle->Eta(), lPercentiles);
      fHistSelectedTriggerPtPhi->Fill(track->Pt(), track->Phi(), lPercentiles);
      fHistSelectedTriggerPtEta->Fill(track->Pt(), track->Eta(), lPercentiles);    
      labelPrimOrSec=1;
    }
    else if(particle->IsSecondaryFromWeakDecay())      labelPrimOrSec=2;
    else if(particle->IsSecondaryFromMaterial())      labelPrimOrSec=3;
    else labelPrimOrSec=4;
    cout << "label is " << labelPrimOrSec<< endl;
    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	for(Int_t p=1; p<=4; p++){
	  if (labelPrimOrSec==p) fHistPrimaryTrigger[m]->Fill(p,particle->Pt() );
	}
      }
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use 
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output event
  //

  /*
  //file collection
  //OpenFile();
  fEventColl = new AliAnalysisKPEventCollectionChiara **[fzVertexBins]; 
  
  for (unsigned short i=0; i<fzVertexBins; i++) {
    fEventColl[i] = new AliAnalysisKPEventCollectionChiara *[fnMultBins];
    for (unsigned short j=0; j<fnMultBins; j++) {
      fEventColl[i][j] = new AliAnalysisKPEventCollectionChiara(fnEventsToMix+1, fMaxFirstMult, fMaxSecondMult);
    }
  }
  // end event collection
  */

  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];
  
  fOutputList = new TList();         
  fOutputList->SetOwner(kTRUE);     

  fSignalTree= new TTree("fSignalTree","fSignalTree");
  fSignalTree->Branch("fTreeVariablePtTrigger",          &fTreeVariablePtTrigger   , "fTreeVariablePtTrigger/D");
  fSignalTree->Branch("fTreeVariableChargeTrigger",      &fTreeVariableChargeTrigger, "fTreeVariableChargeTrigger/D");
  fSignalTree->Branch("fTreeVariableEtaTrigger",         &fTreeVariableEtaTrigger  , "fTreeVariableEtaTrigger/D");
  fSignalTree->Branch("fTreeVariablePhiTrigger",         &fTreeVariablePhiTrigger, "fTreeVariablePhiTrigger/D");
  fSignalTree->Branch("fTreeVariableDCAz",               &fTreeVariableDCAz  , "fTreeVariableDCAz/D");
  fSignalTree->Branch("fTreeVariableDCAxy",              &fTreeVariableDCAxy  , "fTreeVariableDCAxy/D");
  fSignalTree->Branch("fTreeVariableisPrimaryTrigger",              &fTreeVariableisPrimaryTrigger  , "fTreeVariableisPrimaryTrigger/D");
  fSignalTree->Branch("fTreeVariableisPrimaryV0",              &fTreeVariableisPrimaryV0  , "fTreeVariableisPrimaryV0/D");
  fSignalTree->Branch("fTreeVariableRapK0Short",         &fTreeVariableRapK0Short               ,"fTreeVariableRapK0Short/D");
  fSignalTree->Branch("fTreeVariableDcaV0ToPrimVertex",  &fTreeVariableDcaV0ToPrimVertex 	, "fTreeVariableDcaV0ToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaPosToPrimVertex", &fTreeVariableDcaPosToPrimVertex	, "fTreeVariableDcaPosToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaNegToPrimVertex", &fTreeVariableDcaNegToPrimVertex	, "fTreeVariableDcaNegToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableV0CosineOfPointingAngle", &fTreeVariableV0CosineOfPointingAngle	, "fTreeVariableV0CosineOfPointingAngle/D");
  fSignalTree->Branch("fTreeVariablePtV0",               &fTreeVariablePtV0   , "fTreeVariablePtV0/D");
  fSignalTree->Branch("fTreeVariableInvMassK0s",         &fTreeVariableInvMassK0s, "fTreeVariableInvMassK0s/D");
  fSignalTree->Branch("fTreeVariableInvMassLambda",      &fTreeVariableInvMassLambda, "fTreeVariableInvMassLambda/D");
  fSignalTree->Branch("fTreeVariableInvMassAntiLambda",  &fTreeVariableInvMassAntiLambda, "fTreeVariableInvMassAntiLambda/D");
  fSignalTree->Branch("fTreeVariableEtaV0",              &fTreeVariableEtaV0  , "fTreeVariableEtaV0/D");
  fSignalTree->Branch("fTreeVariablePhiV0",              &fTreeVariablePhiV0, "fTreeVariablePhiV0/D");
  fSignalTree->Branch("fTreeVariablePtArmenteros",       &fTreeVariablePtArmenteros  , "fTreeVariablePtArmenteros/D");
  fSignalTree->Branch("fTreeVariableAlpha",              &fTreeVariableAlpha  , "fTreeVariableAlpha/D");
  fSignalTree->Branch("fTreeVariableDeltaEta",              &fTreeVariableDeltaEta  , "fTreeVariableDeltaEta/D");
  fSignalTree->Branch("fTreeVariableDeltaPhi",              &fTreeVariableDeltaPhi, "fTreeVariableDeltaPhi/D");
  fSignalTree->Branch("fTreeVariableDeltaTheta",              &fTreeVariableDeltaTheta, "fTreeVariableDeltaTheta/D");
  fSignalTree->Branch("fTreeVariableMultiplicity",       &fTreeVariableMultiplicity , "fTreeVariableMultiplicity/D");
  fSignalTree->Branch("fTreeVariableZvertex",              &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fSignalTree->Branch("fTreeVariablePDGCode",              &fTreeVariablePDGCode  , "fTreeVariablePDGCode/D");

  fBkgTree= new TTree("fBkgTree","fBkgTree");
  fBkgTree->Branch("fTreeVariablePtTrigger",          &fTreeVariablePtTrigger   , "fTreeVariablePtTrigger/D");
  fBkgTree->Branch("fTreeVariableChargeTrigger",      &fTreeVariableChargeTrigger, "fTreeVariableChargeTrigger/D");
  fBkgTree->Branch("fTreeVariableEtaTrigger",         &fTreeVariableEtaTrigger  , "fTreeVariableEtaTrigger/D");
  fBkgTree->Branch("fTreeVariablePhiTrigger",         &fTreeVariablePhiTrigger, "fTreeVariablePhiTrigger/D");
  fBkgTree->Branch("fTreeVariableDCAz",               &fTreeVariableDCAz  , "fTreeVariableDCAz/D");
  fBkgTree->Branch("fTreeVariableDCAxy",              &fTreeVariableDCAxy  , "fTreeVariableDCAxy/D");
  fBkgTree->Branch("fTreeVariableisPrimaryTrigger",              &fTreeVariableisPrimaryTrigger  , "fTreeVariableisPrimaryTrigger/D");  
  fBkgTree->Branch("fTreeVariableisPrimaryV0",              &fTreeVariableisPrimaryV0  , "fTreeVariableisPrimaryV0/D");  
  fBkgTree->Branch("fTreeVariableRapK0Short",         &fTreeVariableRapK0Short               ,"fTreeVariableRapK0Short/D");
  fBkgTree->Branch("fTreeVariableDcaV0ToPrimVertex",  &fTreeVariableDcaV0ToPrimVertex 	, "fTreeVariableDcaV0ToPrimVertex/D");
  fBkgTree->Branch("fTreeVariableDcaPosToPrimVertex", &fTreeVariableDcaPosToPrimVertex	, "fTreeVariableDcaPosToPrimVertex/D");
  fBkgTree->Branch("fTreeVariableDcaNegToPrimVertex", &fTreeVariableDcaNegToPrimVertex	, "fTreeVariableDcaNegToPrimVertex/D");
  fBkgTree->Branch("fTreeVariableV0CosineOfPointingAngle", &fTreeVariableV0CosineOfPointingAngle	, "fTreeVariableV0CosineOfPointingAngle/D");
  fBkgTree->Branch("fTreeVariablePtV0",               &fTreeVariablePtV0   , "fTreeVariablePtV0/D");
  fBkgTree->Branch("fTreeVariableInvMassK0s",         &fTreeVariableInvMassK0s, "fTreeVariableInvMassK0s/D");
  fBkgTree->Branch("fTreeVariableInvMassLambda",      &fTreeVariableInvMassLambda, "fTreeVariableInvMassLambda/D");
  fBkgTree->Branch("fTreeVariableInvMassAntiLambda",  &fTreeVariableInvMassAntiLambda, "fTreeVariableInvMassAntiLambda/D");
  fBkgTree->Branch("fTreeVariableEtaV0",              &fTreeVariableEtaV0  , "fTreeVariableEtaV0/D");
  fBkgTree->Branch("fTreeVariablePhiV0",              &fTreeVariablePhiV0, "fTreeVariablePhiV0/D");
  fBkgTree->Branch("fTreeVariablePtArmenteros",       &fTreeVariablePtArmenteros  , "fTreeVariablePtArmenteros/D");
  fBkgTree->Branch("fTreeVariableAlpha",              &fTreeVariableAlpha  , "fTreeVariableAlpha/D");
  fBkgTree->Branch("fTreeVariableDeltaEta",              &fTreeVariableDeltaEta  , "fTreeVariableDeltaEta/D");
  fBkgTree->Branch("fTreeVariableDeltaPhi",              &fTreeVariableDeltaPhi  , "fTreeVariableDeltaPhi/D");
  fBkgTree->Branch("fTreeVariableDeltaTheta",              &fTreeVariableDeltaTheta, "fTreeVariableDeltaTheta/D");
  fBkgTree->Branch("fTreeVariableMultiplicity",       &fTreeVariableMultiplicity , "fTreeVariableMultiplicity/D");
  fBkgTree->Branch("fTreeVariableZvertex",              &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fBkgTree->Branch("fTreeVariablePDGCode",              &fTreeVariablePDGCode  , "fTreeVariablePDGCode/D");

  fHistPt = new TH1F("fHistPt", "p_{T} distribution of selected charged tracks in events with NT>0 and NV0>0", 300, 0, 30); 
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtvsMult= new TH2F("fHistPtvsMult", "p_{T} and centrality distribution of selected charged tracks in events with NT>0 and NV0>0", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMult->GetYaxis()->SetTitle("Centrality");

  fHistZvertex= new TH1F("fHistZvertex", "Z vertex distribution of selected events (NT>0 and NV0>0)", 40,-20,20);

  fHist_eta_phi= new TH2F("fHist_eta_phi", "Distribution of selected charged tracks in events with NT>0 and NV0>0", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_multiplicity=new TH1F("fHist_multiplicity", "fHist_multiplicity", 100, 0, 100); 
  fHist_multiplicity->SetTitle("Multiplicity distribution of selected events in events with NT>0 and NV0>0");

  fHistPDG=new TH1F("fHistPDG", "fHistPDG",3200, -3200, 3200);
  
  fHistEventMult=new TH1F("fHistEventMult", "fHistEventMult", 14, 0.5, 14.5);
  fHistEventMult->GetXaxis()->SetBinLabel(1,"All events");
  fHistEventMult->GetXaxis()->SetBinLabel(2,"Events w/PV and AOD");
  fHistEventMult->GetXaxis()->SetBinLabel(3,"Events w/|Vx|<10 cm"); 
  fHistEventMult->GetXaxis()->SetBinLabel(4,"Events w/ PID"); 
  fHistEventMult->GetXaxis()->SetBinLabel(5,"centrality <= 199"); 
  fHistEventMult->GetXaxis()->SetBinLabel(6,"NO PILE UP"); 
  fHistEventMult->GetXaxis()->SetBinLabel(7,"INT7"); 
  fHistEventMult->GetXaxis()->SetBinLabel(8,"ANY"); 
  fHistEventMult->GetXaxis()->SetBinLabel(9,"other type"); 
  fHistEventMult->GetXaxis()->SetBinLabel(10,"Ntrigger>0"); 
  fHistEventMult->GetXaxis()->SetBinLabel(11,"Ntrigger>0 (MC)"); 
  fHistEventMult->GetXaxis()->SetBinLabel(12,"NTrigger>1");
  fHistEventMult->GetXaxis()->SetBinLabel(13,"NTrigger>1 (MC)");
  fHistEventMult->GetXaxis()->SetBinLabel(14,"Ntrigger>0 && NV0>0");//tutti gli eventi in cui ho almeno una particella di trigger ricostruita e almeno una V0 (MC truth)

  fHistEventV0=new TH1F("fHistEventV0", "fHistEventV0",12, 0.5, 12.5);
  fHistEventV0->SetTitle("Number of V0 which progressively pass the listed selections");
  fHistEventV0->GetXaxis()->SetBinLabel(1,"All V0s");
  fHistEventV0->GetXaxis()->SetBinLabel(2,"V0s ok");
  fHistEventV0->GetXaxis()->SetBinLabel(3,"Filterbit daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(4,"Chis daughter tracks"); 
  fHistEventV0->GetXaxis()->SetBinLabel(5,"PID daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(6,"|eta daughters|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(7,"TPC track quality"); 
  fHistEventV0->GetXaxis()->SetBinLabel(8,"|eta_K0s|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(9,"All other cuts except pT"); 
  fHistEventV0->GetXaxis()->SetBinLabel(10,"pT cuts"); 
  fHistEventV0->GetXaxis()->SetBinLabel(11,"pT < pT trigger"); //numbersecondparticle (reco V0)
  fHistEventV0->GetXaxis()->SetBinLabel(12,"pT < pT trigger (MC)"); //numbersecondparticle (true V0)


  fHistTrack=new TH1F("fHistTrack", "fHistTrack", 11, 0.5, 11.5);
  fHistTrack->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistTrack->GetXaxis()->SetBinLabel(2,"Tracks after filterbit");
  fHistTrack->GetXaxis()->SetBinLabel(3,"Tracks with |eta| < 0.8"); 
  fHistTrack->GetXaxis()->SetBinLabel(4,"Track quality");
  fHistTrack->GetXaxis()->SetBinLabel(5,"TPCCrossedRows>70");
  fHistTrack->GetXaxis()->SetBinLabel(6,"Crossed rows/findable >0.8");
  fHistTrack->GetXaxis()->SetBinLabel(7,"Charged tracks");
  fHistTrack->GetXaxis()->SetBinLabel(8,"N.trigger"); //NumberPrimary in all slected events 
  fHistTrack->GetXaxis()->SetBinLabel(9,"N.trigger MC");
  fHistTrack->GetXaxis()->SetBinLabel(10,"N.trigger (NV0>0)"); //NumberPrimary in events with at least one V0 (one reco V0 for data, one true V0 for MC)
  fHistTrack->GetXaxis()->SetBinLabel(11,"N.trigger (NV0>0) MC");

  fMassV0= new TH1F("fMassV0", "Invariant mass of V0 candidates", 100, 0.45, 0.55);
  fMassV0->GetXaxis()->SetTitle("M_{#pi^+ #pi^-}");

  fHistMassvsPt = new TH2F *[5];
  fHistMassvsPt_tagli = new TH2F *[5];

  for(Int_t j=0; j<5; j++){
    fHistMassvsPt[j] = new TH2F(Form("fHistMassvsPt_" +fV0+ "_%i",j),Form("fHistMassvsPt_"+fV0 + "_%i"+" (cuts on PID, eta and track quality of daughters + eta V0) ",j),400,0.3,0.7,160,0,16); 
    fHistMassvsPt_tagli[j] = new TH2F(Form("fHistMassvsPt_" +fV0+ "_tagli_%i",j),Form("fHistMassvsPt_" +fV0+ "_tagli_%i" + " (all selections on V0 applied)",j),400,0.3,0.7,160,0,16);
    fHistMassvsPt[j]->GetXaxis()->SetTitle("Invariant mass of V0 candidate");
    fHistMassvsPt[j]->GetYaxis()->SetTitle("p_{T} of V0 candidate");   
    fHistMassvsPt_tagli[j]->GetXaxis()->SetTitle("Invariant mass of V0 candidate");
    fHistMassvsPt_tagli[j]->GetYaxis()->SetTitle("p_{T} of V0 candidate");   
  }
  
  fHistMultvsTrigger=new TH2F("fHistMultvsTrigger", "Multiplicity of selected events (T>0, V>0) vs number of trigger particles", 20, 0, 20, 100, 0, 100);
  fHistMultvsTrigger->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTrigger->GetYaxis()->SetTitle("Multiplicity");
  

  fHistMultvsTriggerMCTruth=new TH2F("fHistMultvsTriggerMCTruth", "Multiplicity of selected events (T>0, V>0) vs number of trigger particles, MC Truth", 20, 0, 20, 100, 0, 100);
  fHistMultvsTriggerMCTruth->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruth->GetYaxis()->SetTitle("Multiplicity");
  

  fHistMultvsV0=new TH2F("fHistMultvsV0", "Multiplicity of selected events (T>0, V0>0) vs number of V0s",20, 0, 20,100, 0, 100 );
  fHistMultvsV0->GetXaxis()->SetTitle("Number of V0 particles");
  fHistMultvsV0->GetYaxis()->SetTitle("Multiplicity");
  
  fHistPtArmvsAlpha=new TH2F("fHistPtArmvsAlpha", "Distribution of V0 candidates before cuts on Armenteros variables",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlpha->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlpha->GetYaxis()->SetTitle("Pt Armenteros");

  fHistPtArmvsAlphaAfterSelection=new TH2F("fHistPtArmvsAlphaAfterSelection", "Distribution of V0 candidates after applying cuts on Armenteros", 80, -1, 1,80, 0, 0.3);

  fHistTrigger=new TH1F("fHistTrigger", "Number of reco trigger particle distribution for selected events (also for events with NV0=0 (reco if data, true if MC))", 20, -0.5, 19.5); // each entry is an event

  fHistTriggerwV0=new TH1F("fHistTriggerwV0", "Number of true trigger particle distribution for selected events (for events with NV0>0 (reco if data, true if MC()", 20, -0.5, 19.5); // each entry is an event

  fHistTriggerMCTruth=new TH1F("fHistTriggerMCTruth", "Number of trigger particle distribution for selected events, MC Truth (also for events with NV0=0 (reco if data, true if MC))", 20, -0.5, 19.5); // each entry is an event

  fHistTriggerwV0MCTruth=new TH1F("fHistTriggerwV0MCTruth", "Number of trigger particle distribution for selected events, MC Truth (for events with NV0>0 (reco if data, true if MC))", 20, -0.5, 19.5); // each entry is an event

  fHistMultiplicityVsVertexZ=new TH2F("fHistMultiplicityVsVertexZ", "Multiplicity vs Z vertex of selected events with NT>0 and NV0>0 ",  20, -10, 10,100, 0, 100);
      
  fHistTriggervsMult=new TH1F("fHistTriggervsMult", "Numero di particelle di trigger nei vari intervalli di molteplicita'", 100, 0, 100);
  fHistTriggervsMult->GetXaxis()->SetTitle("Centrality");

  fHistGeneratedTriggerPtPhi=new TH3F("fHistGeneratedTriggerPtPhi", "p_{T} and #phi distribution of generated trigger particles (charged, primary)", 300, 0, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
  fHistGeneratedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  fHistGeneratedTriggerPtEta=new TH3F("fHistGeneratedTriggerPtEta", "p_{T} and #eta distribution of generated trigger particles(primary, charged)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistGeneratedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtEta->GetYaxis()->SetTitle("#eta");

  fHistSelectedTriggerPtPhi=new TH3F("fHistSelectedTriggerPtPhi", "p_{T} and #phi distribution of selected trigger particles (primary)", 300, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
  fHistSelectedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  fHistSelectedTriggerPtEta=new TH3F("fHistSelectedTriggerPtEta", "p_{T} and #eta distribution of selected trigger particles (primary)", 300, 0, 30, 400,-1.2, 1.2,  100, 0, 100);
  fHistGeneratedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtEta->GetYaxis()->SetTitle("#eta");

  fHistGeneratedV0PtPhi=new TH3F("fHistGeneratedV0PtPhi", "p_{T} and #phi distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
  fHistGeneratedV0PtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedV0PtPhi->GetYaxis()->SetTitle("#phi");

  fHistSelectedV0PtPhi=new TH3F("fHistSelectedV0PtPhi", "p_{T} and #phi distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
  fHistSelectedV0PtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtPhi->GetYaxis()->SetTitle("#phi");

  fHistGeneratedV0PtEta=new TH3F("fHistGeneratedV0PtEta", "p_{T} and #eta distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistGeneratedV0PtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedV0PtEta->GetYaxis()->SetTitle("#eta");

  fHistSelectedV0PtEta=new TH3F("fHistSelectedV0PtEta", "p_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistSelectedV0PtEta->GetXaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtEta->GetYaxis()->SetTitle("#eta");

  fHistReconstructedV0PtMass=new TH3F("fHistReconstructedV0PtMass", "p_{T} and mass distribution of reconstructed V0 particles(K0s, primary, event w T>0)", 100, 0.45, 0.55, 160, 0, 16,  100, 0, 100);
  fHistReconstructedV0PtMass->GetYaxis()->SetTitle("p_{T}");
  fHistReconstructedV0PtMass->GetXaxis()->SetTitle("M_{pi^{+} #pi^{-}}");

  fHistSelectedV0PtMass=new TH3F("fHistSelectedV0PtMass", "p_{T} and mass distribution of selected V0 particles (K0s, primary, event w T>0)", 100, 0.45, 0.55,  160, 0, 16, 100, 0, 100);
  fHistSelectedV0PtMass->GetYaxis()->SetTitle("p_{T}");
  fHistSelectedV0PtMass->GetXaxis()->SetTitle("M_{pi^{+} #pi^{-}}");

  fHistResolutionTriggerPt=new TH2F("fHistResolutionTriggerPt", "p_{T} resolution of selected trigger particles (primary)", 500, -0.5, 0.5, 100, 0, 100);
  fHistResolutionTriggerPt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionTriggerPt->GetYaxis()->SetTitle("Centrality");

  fHistResolutionTriggerPhi=new TH2F("fHistResolutionTriggerPhi", "#Phi resolution of selected trigger particles (primary)", 500, -0.5, 0.5, 100, 0, 100);
  fHistResolutionTriggerPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhi->GetYaxis()->SetTitle("Centrality");

  fHistResolutionTriggerEta=new TH2F("fHistResolutionTriggerEta", "#Eta resolution of selected trigger particles (primary)", 500, -0.5, 0.5, 100, 0, 100);
  fHistResolutionTriggerEta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionTriggerEta->GetYaxis()->SetTitle("Centrality");

  fHistResolutionV0Pt=new TH2F("fHistResolutionV0Pt", "p_{T} resolution of selected V0 particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 100, 0, 100);
  fHistResolutionV0Pt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionV0Pt->GetYaxis()->SetTitle("Centrality");

  fHistResolutionV0Phi=new TH2F("fHistResolutionV0Phi", "#Phi resolution of selected V0 particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 100, 0, 100);
  fHistResolutionV0Phi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionV0Phi->GetYaxis()->SetTitle("Centrality");
  
  fHistResolutionV0Eta=new TH2F("fHistResolutionV0Eta", "#Eta resolution of selected V0 particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 100, 0, 100);
  fHistResolutionV0Eta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionV0Eta->GetYaxis()->SetTitle("Centrality");

  fHistPrimaryTrigger= new TH2F*[5];
  for(Int_t j=0; j<5; j++){
    fHistPrimaryTrigger[j]=new TH2F(Form("fHistPrimaryTrigger_%i", j), "Trigger MC (selected)", 4, 0.5, 4.5, 100, 0,30 );
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(1,"Primary selected triggers");
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected triggers"); 
    fHistPrimaryTrigger[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected triggers"); 
    fHistPrimaryTrigger[j]->GetYaxis()->SetTitle("p_{T}"); 
  }

  fHistPrimaryV0= new TH2F*[5];
  for(Int_t j=0; j<5; j++){
    fHistPrimaryV0[j]=new TH2F(Form("fHistPrimaryV0_%i",j), "V0 MC (K0s, selected)", 4, 0.5, 4.5, 160, 0, 16);
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s"); 
    fHistPrimaryV0[j]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s"); 
    fHistPrimaryV0[j]->GetYaxis()->SetTitle("p_{T}"); 
  }

  //  TString molteplicit[6]={"0-7","7-15","15-25","25-40","40-70",">70"};
  fHistMultiplicityOfMixedEvent=new TH1F("fHistMultiplicityOfMixedEvent", "Distribution of number of events used for the mixing", 10, 0, 10);

  fEventCuts.AddQAplotsToList(fOutputList);

  fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
  fOutputList->Add(fHistPtvsMult);       
  fOutputList->Add(fHist_eta_phi);// your histogram in the output file, add it to the list!
  fOutputList->Add(fHistZvertex);
  fOutputList->Add(fHist_multiplicity);
  fOutputList->Add(fHistPDG);
  fOutputList->Add(fHistEventMult);
  fOutputList->Add(fHistEventV0);
  fOutputList->Add(fHistTrack); 
  fOutputList->Add(fMassV0);
  fOutputList->Add(fHistMultvsV0);

  for(Int_t j=0; j < 5; j++){
    fOutputList->Add(fHistMassvsPt[j]);
    fOutputList->Add(fHistMassvsPt_tagli[j]);  
    fOutputList->Add(fHistPrimaryTrigger[j]); 
    fOutputList->Add(fHistPrimaryV0[j]);      

  }
  
  fOutputList->Add(fHistPtArmvsAlpha);
  fOutputList->Add(fHistPtArmvsAlphaAfterSelection);  
  fOutputList->Add(fHistMultvsTrigger);
  fOutputList->Add(fHistMultvsTriggerMCTruth);
  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistTriggerwV0);
  fOutputList->Add(fHistTriggerMCTructh);
  fOutputList->Add(fHistTriggerwV0MCTruth);
  fOutputList->Add(fHistMultiplicityVsVertexZ);
  fOutputList->Add(fHistTriggervsMult);
  fOutputList->Add(fHistMultiplicityOfMixedEvent);
  fOutputList->Add(fHistGeneratedTriggerPtPhi);
  fOutputList->Add(fHistSelectedTriggerPtPhi);
  fOutputList->Add(fHistGeneratedTriggerPtEta);
  fOutputList->Add(fHistSelectedTriggerPtEta);
  fOutputList->Add(fHistGeneratedV0PtPhi); 
  fOutputList->Add(fHistSelectedV0PtPhi);
  fOutputList->Add(fHistGeneratedV0PtEta);
  fOutputList->Add(fHistSelectedV0PtEta);
  fOutputList->Add(fHistReconstructedV0PtMass);
  fOutputList->Add(fHistSelectedV0PtMass);
  fOutputList->Add(fHistResolutionTriggerPt);
  fOutputList->Add(fHistResolutionTriggerPhi);
  fOutputList->Add(fHistResolutionTriggerEta);
  fOutputList->Add(fHistResolutionV0Pt);
  fOutputList->Add(fHistResolutionV0Phi);
  fOutputList->Add(fHistResolutionV0Eta);
  
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
  // fOutputList object. the manager will in the end take care of writing your output to file
  // so it needs to know what's in the output
  PostData(2, fSignalTree);       
  PostData(3, fBkgTree);       
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
  //cout<<"enter user exec"<<endl;
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you 
  // have access to the current event. 
  // once you return from the UserExec function, the manager will retrieve the next event from the chain

  //  Float_t moltep[5]={0,7,15,25,40,70}; //valori associati a molteplicita'
  Float_t moltep[6]={0,5,10,30,50,100};  //valori associati a centralita'

  fHistEventMult->Fill(1);

  // AliMCEvent   *lMCevent  = 0x0;
  // AliStack     *lMCstack  = 0x0;
  // TClonesArray *arrayMC = 0x0;
   
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());  
  if(!fAOD) {
    AliWarning("Error: AOD event not available \n");
    PostData(1, fOutputList);       
    PostData(2, fSignalTree);       
    PostData(3,fBkgTree); return;
  }        
  
 
 /// Use the event cut class to apply the required selections
 if (!fEventCuts.AcceptEvent(fAOD)) {   
 PostData(1, fOutputList);
 PostData(2, fSignalTree );
 PostData(3,fBkgTree); return;
 }
 
  Int_t iTracks(fAOD->GetNumberOfTracks());         
  Int_t V0Tracks(fAOD->GetNumberOfV0s());           
  Evcounter++;  
  cout << "\n \n \n ********************************************************* "<< endl;
  cout << "numero dell'evento "<<  Evcounter << endl; 
  cout << "number of tracks before any cut " << iTracks << endl;
  cout << "number of V0 before any cut " << V0Tracks << endl;
  
  //VERTEX SELECTION AND TRIGGER
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  const AliAODVertex *lPrimaryBestAODVtx = fAOD->GetPrimaryVertex();
  if (!lPrimaryBestAODVtx){
    AliWarning("No prim. vertex in AOD... return!");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); return;
  }
  fHistEventMult->Fill(2);
  AliVVertex *vertexmain =0x0;
  vertexmain = (AliVVertex*) lPrimaryBestAODVtx;
  lPrimaryBestAODVtx->GetXYZ(lBestPrimaryVtxPos);


  if (TMath::Abs(lBestPrimaryVtxPos[2])>10.){
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    cout << "z vertex selection failed " << endl;
    PostData(3,fBkgTree); return;
  }
  fHistEventMult->Fill(3);

  /* ho implementato in modo diverso
     if (fReadMCTruth) {
     //      Printf("Reading MC truth!!! \n");
     arrayMC = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    
     if (!arrayMC) AliFatal("Error: MC particles branch not found!\n");
    
     }
  */

  //PID: retrieve AliPIDResponse object from manager + connect our pointer to this object
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  //  if (man) {
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
  UInt_t mask = inputHandler->IsEventSelected();
  //}

  if (!fPIDResponse){
    AliWarning("cannot get pid response");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); return;
  }
  fHistEventMult->Fill(4);


  //CENTRALITY SELECTIONS for Pb-Pb //preso da AliAnalysisTaskKPFemto
  //lcentrality=centrality->GetCentralityPercentile("V0A");
  
  Float_t lPercentiles = 0;
  
  
  //This will work for both ESDs and AODs
  AliMultSelection *MultSelection = (AliMultSelection*) fAOD -> FindListObject("MultSelection");

  if ( MultSelection ){
    cout << "mult sel ok" << endl;
    lPercentiles= MultSelection->GetMultiplicityPercentile("V0M");
  }else{
    AliInfo("Didn't find MultSelection!"); 
  }
   
  cout << "centralita' dell evento " << lPercentiles << endl;
  //Embedded in framework:
  // ---> if lPercentiles are 200 or above: event didn't pass event selections
  // ---> if lPercentiles are 0-100       : event passed, go ahead
 

  //// tentativo 
  //   Float_t lPercentiles = 0;
  //lPercentiles= fEventCuts.GetCentrality(0);
  //cout << "centralita' dell evento " << lPercentiles << endl;
  
   
  if ( lPercentiles > 199 ){
    PostData(1,fOutputList );
    PostData(2, fSignalTree );
    cout << "lPercentiles >= 200" << endl;
    PostData(3,fBkgTree); return;   
  }

  /*da non utilizzare ora

    Int_t lcentrality=((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08();//CONTROLLA!!
    cout << "lcentrality  " << lcentrality << endl;
    //fMultSelection =(AliMultSelection*)(fAOD->FindListObject("MultSelection"));
    //if(fMultSelection) centrality = fMultSelection->GetMultiplicityPercentile("V0M");
    if (lcentrality==-1){
    PostData(1,fOutputList );
    PostData(2, fSignalTree );
    cout << "return: lcentrality= -1 " << endl;
    PostData(3,fBkgTree); return;
    }
  */
  fHistEventMult->Fill(5);

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fAOD->IsPileupFromSPD();
  if(isPileUpSpd){ 
    PostData(1,fOutputList );
    PostData(2, fSignalTree );
    cout << "return: event is pile up " << endl;
    PostData(3,fBkgTree); return;
  }
  fHistEventMult->Fill(6);

  Bool_t isSelectedInt7        = kFALSE;
  Bool_t isSelectedAny         = kFALSE;
  Bool_t isSelected            = kFALSE;
 
  if(fCollidingSystem == "pp"){
    isSelectedInt7        = (mask & AliVEvent::kINT7);
    isSelectedAny         = (mask & AliVEvent::kAnyINT);
    
    if(isSelectedInt7 )
      isSelected = kTRUE;
  }
  
  else if(fCollidingSystem == "pPb"){
    isSelectedInt7        = (mask & AliVEvent::kINT7);
    isSelectedAny         = (mask & AliVEvent::kAny);

    if(isSelectedInt7 )
      isSelected = kTRUE;
  }
  
  if(isSelectedInt7)
    fHistEventMult->Fill(7);
  else if(isSelectedAny)
    fHistEventMult->Fill(8);
  else
    fHistEventMult->Fill(9) ; 
  
  // cout<<"Trigger mask: "<<fAODevent->GetTriggerMask()<<endl;
  // cout<<"Event type  : "<<fAODevent->GetEventType()<<endl;
  // FIXME : event selection to be added.. DONE

    
  if(!isSelected){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    cout << "event does not fulfil centrality selection criteria " << endl;     
    PostData(3,fBkgTree);
    return;
  }
  
  cout << "event has passed selection criteria.... first and second particles to be analyzed ...."<< endl;


  //MC
  Bool_t Generated=kTRUE; //TRUE se devo produrre istogrammi delle generate
  Int_t labelPrimOrSec=0; //mi dice se la particella e' primaria, secondaria,...
  Int_t label = 0; //il label della traccia
  AliAODTrack* track=0x0;

  Bool_t isV0=kFALSE;
  if(fReadMCTruth){
    fMCEvent= MCEvent();
    if (fMCEvent){
      ProcessMCParticles(Generated, track, labelPrimOrSec, lPercentiles, isV0);
    }
  }
  
  
  //  cout << "dopo centrality selection " << endl;
  const Float_t bfield = (InputEvent())->GetMagneticField();
  int fieldsign;
  if (bfield >=0.) fieldsign = 1;
  else fieldsign = -1;
 
  //Store event in the buffer to do mixing
  //find vertex...
  int zBin=0;
  int centralityBin=0;
    
  double zStep=2*10/double(fzVertexBins), zStart=-10.;
    
  for (int i=0; i<fzVertexBins; i++) {
    if ((lBestPrimaryVtxPos[2] > zStart+i*zStep) && (lBestPrimaryVtxPos[2] < zStart+(i+1)*zStep)) {
      zBin=i;
      break;
    }
  }
  //cout << " zbin " << zBin<<endl;
  //...and centrality
  if(lPercentiles < 0.01) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
  else if(lPercentiles < 0.1) centralityBin=18;
  else if(lPercentiles < 0.5) centralityBin=17;
  else if(lPercentiles < 1.0) centralityBin=16;
  else if(lPercentiles < 5.0) centralityBin=15;
  else if(lPercentiles < 10.) centralityBin=14;
  else if(lPercentiles < 20.) centralityBin=13;
  else if(lPercentiles < 30.) centralityBin=12;
  else if(lPercentiles < 40.) centralityBin=11;
  else if(lPercentiles < 50.) centralityBin=10;
  else if(lPercentiles < 70.) centralityBin=9; 
  else if(lPercentiles <= 100.) centralityBin=8;

  /*
    if(lPercentiles < 5.) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
    else if(lPercentiles < 10.) centralityBin=18;
    else if(lPercentiles < 15.) centralityBin=17;
    else if(lPercentiles < 20.) centralityBin=16;
    else if(lPercentiles < 25.) centralityBin=15;
    else if(lPercentiles < 30.) centralityBin=14;
    else if(lPercentiles < 35.) centralityBin=13;
    else if(lPercentiles < 40.) centralityBin=12;
    else if(lPercentiles < 45.) centralityBin=11;
    else if(lPercentiles < 50.) centralityBin=10;
    else if(lPercentiles < 55.) centralityBin=9; 
    else if(lPercentiles < 60.) centralityBin=8;
    else if(lPercentiles < 65.) centralityBin=7;
    else if(lPercentiles < 70.) centralityBin=6;
    else if(lPercentiles < 75.) centralityBin=5;
    else if(lPercentiles < 80.) centralityBin=4;
    else if(lPercentiles < 85.) centralityBin=3;
    else if(lPercentiles < 90.) centralityBin=2;   
    else if(lPercentiles < 95.) centralityBin=1;
    else if(lPercentiles <= 100.) centralityBin=0;
  */
  //cout << "6 " << endl;
  cout << "  " << zBin <<" centrbin " << centralityBin<< endl;
  if (((centralityBin+1) >fnMultBins) || ((zBin+1) > fzVertexBins)){ 
    cout<<" ##################  WARNING: I'm going to break bacause of dimensional issues ########################"<<endl;
  }
  // //se decommento, ho seg viol
  
  //try4 if(Evcounter==1) FifoShiftok=kTRUE;
    
  //try 4 if (FifoShiftok==kTRUE)  fEventColl[zBin][centralityBin]->FifoShift();
  //try6  fEventColl[zBin][centralityBin]->FifoShift();
  //cout << "6 " << endl;
  FifoShiftok=kFALSE;
  // // //se decommento, ho seg viol 
  //try6 fEvt = fEventColl[zBin][centralityBin]->fEvt;
  // // cout << "7 " << endl;
  //-----------------------------------LOOP OVER THE TRACKS

  Float_t nTPCCrossedRows=0.;
  float rationCrnFind=0;
  Float_t nsigmaTOFj=3;
  Float_t nsigmaTPCj=3;
  Int_t charge=0;
  Int_t NumberFirstParticle=0;
  Int_t NumberFirstParticleMC=0;
  Int_t NumberFirstParticle_finale=0;
  Int_t NumberSecondParticle=0;
  Int_t NumberSecondParticleMC=0;
  Int_t NumberSecondParticleAll=0;
  Double_t Ptintermediate=0;
  Double_t selectedtrackID=0;
  Int_t pos0or1=0;
  Int_t neg0or1=0;
  Int_t CharegFirstParticle=0;
  Double_t dz[2] = {-999.,-999.};

  //LOOP FOR FIRST PARTICLE
  //cout << "ciao " << endl;
 Float_t ptTriggerMinimoDati=10000;
  for(Int_t i=0; i < iTracks; i++) {
    //cout << " ciao, I enter the traxck loop " << endl;   
    track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));        
    fHistTrack->Fill(1);
    if(!track) continue;
    //if(track->GetID()<0)continue;
    //cout << "T1 " << endl;
    //  if(!track->IsOn(AliAODTrack::kTPCrefit)) continue;
    if(!track->TestFilterBit(128)) continue;
    //cout << "T2 " << endl;          
    fHistTrack->Fill(2);
    if(TMath::Abs(track->Eta())>0.8)  continue;
    fHistTrack->Fill(3);
    //cout << "T3 " << endl;
    
    if(track->Chi2perNDF()>4.)continue;
    fHistTrack->Fill(4);
    //cout << "T4 " << endl;
    
    nTPCCrossedRows=track->GetTPCNCrossedRows();
    if(nTPCCrossedRows<70) continue;
    fHistTrack->Fill(5);
    //cout << "T5 " << endl;
    
    rationCrnFind=nTPCCrossedRows/track->GetTPCNclsF();
    if(rationCrnFind<0.8)  continue;
    fHistTrack->Fill(6);
    //cout << "T6 " << endl;
    //cout << "pT della traccia " <<track->Pt()<< endl;

    if((track->Charge())==0) continue;
    fHistTrack->Fill(7);

    
    dz[0] = track->DCA();    // the TPC one should be applied the other biases the CF --> from Maciejs note --> FIXME to be checked 
    dz[1] = track->ZAtDCA(); // for those two lines check AliAODTrack.h // FIXME these two lines produce shifted distributions, known problem, asked Mac and Marian. 


    //to determine efficiency of trigger particle
    label=track->GetLabel();
    Generated = kFALSE;     
    if(fReadMCTruth){
      if(fMCEvent){
	ProcessMCParticles(Generated, track, labelPrimOrSec, lPercentiles, isV0);
      }
    }

    if(track->Pt()> fminPtj && track->Pt()<fmaxPtj){
      if(track->Pt()< ptTriggerMinimoDati) ptTriggerMinimoDati=track->Pt(); 
      if (track->Pt()>Ptintermediate){
	//cout << "T8 " << endl;
	Ptintermediate=track->Pt();
	//cout <<"Pt intermediate " <<  Ptintermediate << endl;
	selectedtrackID= track->GetID();
	//cout << "selectedtrackID " << selectedtrackID << endl;
	NumberFirstParticle_finale= NumberFirstParticle;     
	
      }
      
      if(!(fReadMCTruth)){
      //save first particle information (leading particle)
      // fEvt->fReconstructedFirst[NumberFirstParticle].fCharge       = track->Charge();
      // fEvt->fReconstructedFirst[NumberFirstParticle].fPt           = track->Pt();
      // fEvt->fReconstructedFirst[NumberFirstParticle].fEta          = track->Eta();
      // fEvt->fReconstructedFirst[NumberFirstParticle].fPhi          = track->Phi();
      // fEvt->fReconstructedFirst[NumberFirstParticle].fTheta        = track->Theta();
      // fEvt->fReconstructedFirst[NumberFirstParticle].fDCAz         = dz[1];
      // fEvt->fReconstructedFirst[NumberFirstParticle].fDCAxy        = dz[0];
      // fEvt->fReconstructedFirst[NumberFirstParticle].fMultiplicity = lPercentiles;
      // fEvt->fReconstructedFirst[NumberFirstParticle].fZvertex      = lBestPrimaryVtxPos[2];
      // fEvt->fReconstructedFirst[NumberFirstParticle].isP           = labelPrimOrSec;
      

      NumberFirstParticle++;
      }
        
    }

    //to get the contamination factor: ratio between 3+4 and 1: check if it changes placing these lines after selection of events with V0>0 and T>0 
    // cout << labelPrimOrSec << endl;
  }//end loop on tracks for first particl


 TClonesArray* AODMCTrackArray =0x0;  
 Float_t ptTriggerMinimoMC=10000;
 if(fReadMCTruth){
      if (fMCEvent){
	AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}
	for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
	 
	  AliAODMCParticle* trParticle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
	  
	  if (!trParticle) continue;
	  
	  if((trParticle->Charge())==0)continue;
	  if(TMath::Abs(trParticle->Eta())>0.8)continue; //I need to select particles within this eta range!
	  if (!(trParticle->IsPhysicalPrimary()))continue; 
	  if(trParticle->Pt()<= fminPtj && trParticle->Pt()>=fmaxPtj)continue;
	  if(trParticle->Pt()< ptTriggerMinimoMC) ptTriggerMinimoMC=trParticle->Pt();	  
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fCharge       = trParticle->Charge();
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fPt           = trParticle->Pt();
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fEta          = trParticle->Eta();
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fPhi          = trParticle->Phi();
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fTheta        = trParticle->Theta();
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fDCAz         = 0;
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fDCAxy        = 0;
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fMultiplicity = lPercentiles;
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].fZvertex      = lBestPrimaryVtxPos[2];
	  // fEvt->fReconstructedFirst[NumberFirstParticleMC].isP           = labelPrimOrSec;
	  NumberFirstParticleMC++;

	  	  
	}
      }
    }
 

  fHistTrigger->Fill(NumberFirstParticle);
  fHistTriggerMCTruth->Fill(NumberFirstParticleMC);

  fHistTrack->AddBinContent(8, NumberFirstParticle);
  fHistTrack->AddBinContent(9, NumberFirstParticleMC);
  if (NumberFirstParticle>0) fHistEventMult->Fill(10);   
  if (NumberFirstParticleMC>0) fHistEventMult->Fill(11);   
  if(NumberFirstParticle>1)fHistEventMult->Fill(13); 
  if(NumberFirstParticleMC>1)fHistEventMult->Fill(14); 

  if(NumberFirstParticle==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    cout  << "event does not have Trigger particles " << endl;     
    PostData(3,fBkgTree);
    return;
  }


  cout <<"candidati first particle nell'evento analizzato" << "   " <<  NumberFirstParticle << endl;        
  cout <<"candidati first particle nell'evento analizzato, MC truth" << "   " <<  NumberFirstParticleMC << endl;        
  cout << "\n \n \n Pt della particella con Pt maggiore " << Ptintermediate << endl;
   
  //LOOP FOR SECOND PARTICLE
  Float_t ycut[2]={0.5, 0.5};
  Float_t PosEtaCut[2]={0.8, 0.8};
  Float_t NegEtaCut[2]={0.8, 0.8};
  Float_t nsigmaTPC=3;
  Float_t nsigmaTOF=3;
  Float_t kMaxTPCSigmaPion=3;
  Double_t rapidityV0[2]={0};
  Double_t EV0[2]={0};
  bool goodPiPlus=kFALSE;
  bool goodPiMinus=kFALSE;
  bool goodPiPlusTPC=kFALSE;
  bool goodPiMinusTPC=kFALSE;
  bool goodPiPlusTOF=kFALSE;
  bool goodPiMinusTOF=kFALSE;
  
  Float_t kMaxTOFSigmaPion=3;
  Float_t kTOFLow[2]={1000, 1000}; //per ora e alto
  Float_t kMaxDCA[2]={0.5, 1000};
  Float_t kMinCosAngle[2]={0.997, 0.997};
  Float_t kMaxDCADaughters[2]={1, 1};
  Float_t kMinDCAPrimary[2]={0.06, 0.06};
  Float_t kMinDL[2]={0.5, 0.5};
  Float_t kctauval[2]={20, 30};
  Float_t kctau[2]={0, 0};
  Float_t kEtacut[2]={0.8,0.8};
  Float_t Mass[2]={0.497611, 1.115683};
  Int_t ParticleType=0;
  AliPIDResponse::EDetPidStatus statusTOFPos;
  AliPIDResponse::EDetPidStatus statusTOFNeg;

  Int_t labelPos=0;
  Int_t labelNeg=0;
  AliAODMCParticle* particlePos;
  AliAODMCParticle* particleNeg;
  Int_t PdgPos=0;
  Int_t PdgNeg=0;
  Int_t labelMotherPos=0;
  Int_t labelMotherNeg=0;
  AliAODMCParticle* MotherPos;
  AliAODMCParticle* MotherNeg;
  Int_t PdgMotherPos=0;
  Int_t PdgMotherNeg=0;
  Int_t V0PDGCode=0;
  

  if(fV0=="kK0s") ParticleType =0;
  if(fV0=="kLambda") ParticleType=1;

  //MC generated V0
  isV0=kTRUE;
  Generated=kTRUE;

  if(fReadMCTruth){
    if (fMCEvent){
      ProcessMCParticles(Generated, track, labelPrimOrSec, lPercentiles, isV0);
    }
  }
  
  cout << "\n \n  here I start the loop on v0s " << endl;
  for(Int_t i(0); i < V0Tracks; i++) {       
    cout << "i "<< i << endl;
    fHistEventV0->Fill(1);
    rapidityV0[2]={0};    
    EV0[2]={0};
    kctau[2]={0};
    goodPiPlus=kFALSE;
    goodPiMinus=kFALSE;
    AliAODv0* v0 = fAOD->GetV0(i);
    if(!v0) continue;       
    fHistEventV0->Fill(2);
    AliAODTrack* tempTrack=(AliAODTrack*)v0->GetDaughter(0);
    if(tempTrack->Charge()>0) {pos0or1=0; neg0or1=1;}
    else {pos0or1=1; neg0or1=0;}
    AliAODTrack* prongTrackPos=(AliAODTrack*)v0->GetDaughter(pos0or1);
    AliAODTrack* prongTrackNeg=(AliAODTrack*)v0->GetDaughter(neg0or1);

    // TClonesArray* AODMCTrackArray =0x0;  
    if (fReadMCTruth){
      cout << "\n sono in MC event: I analyze all v0 reconstructed by ALICE" << endl;
      if (fMCEvent){
	AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}
	labelPos=prongTrackPos->GetLabel();
	labelNeg=prongTrackNeg->GetLabel();
	cout <<	"label tracce figlie (pos e neg) "<< labelPos<< endl;
	cout <<	"label tracce figlie (pos e neg) "<< labelNeg<< endl;
	particlePos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelPos)));
	particleNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelNeg)));
	if(labelPos>=0)	PdgPos = particlePos->GetPdgCode();
	if(labelNeg>=0)	PdgNeg = particleNeg->GetPdgCode();
	// if(labelPos >=0) cout << "pdg code for label>=0 " <<	PdgPos << endl;
	// if(labelPos <0) cout << "pdg code for label<0 " <<	PdgPos << endl;
	// if(labelNeg >=0) cout << "pdg code for label>=0 " <<	PdgNeg << endl;
	// if(labelNeg <0) cout << "pdg code for label<0 " <<	PdgNeg << endl;

	labelMotherPos=particlePos->GetMother();
	labelMotherNeg=particleNeg->GetMother();
	cout << "label tracce madri (pos e neg) " << 	labelMotherPos<< endl;
	cout << "label tracce madri (pos e neg) " << 	labelMotherNeg<< endl;

	MotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelMotherPos));
	MotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelMotherNeg));
	//	cout << "sono qua " << endl;
	if (labelMotherPos>=0) PdgMotherPos = MotherPos->GetPdgCode();
	if (labelMotherNeg>=0)	PdgMotherNeg = MotherNeg->GetPdgCode();
	//	cout << "e ora sono qua "<< endl;
      }
    }

    //    cout << "sono oltre il primo readMC " << endl;
    if(ParticleType==0) EV0[ParticleType]= TMath::Sqrt(pow(v0->MassK0Short(),2)+ pow( v0->Px(),2)+ pow( v0->Py(),2)+ pow( v0->Pz(),2));
    if(ParticleType==1) EV0[ParticleType]= TMath::Sqrt(pow(v0->MassLambda(),2)+ pow( v0->Px(),2)+ pow( v0->Py(),2)+ pow( v0->Pz(),2));

    rapidityV0[ParticleType]  = 0.5*TMath::Log(( EV0[ParticleType] + v0->Pz()) / (EV0[ParticleType]   - v0->Pz()) );
    
    //attenzione! era sbagliato
    fTreeVariablePtV0=TMath::Sqrt( pow( v0->Px(),2)+ pow( v0->Py(),2));
    fTreeVariableInvMassK0s=v0->MassK0Short();    
    fTreeVariableInvMassLambda=v0->MassLambda();    
    fTreeVariableInvMassAntiLambda=v0->MassAntiLambda();    
    kctau[ParticleType]=Mass[ParticleType]*v0->DecayLengthV0(lBestPrimaryVtxPos)/TMath::Sqrt( pow( v0->Px(),2)+ pow( v0->Py(),2) + pow(v0->Pz(),2));
    //-------------------------daughter cuts---------------

    Float_t nTPCCrossedRowspos=prongTrackPos ->GetTPCNCrossedRows();
    Float_t nTPCCrossedRowsneg=prongTrackNeg ->GetTPCNCrossedRows();
    Float_t rationCrnFindpos=nTPCCrossedRowspos/prongTrackPos->GetTPCNclsF();
    Float_t rationCrnFindneg=nTPCCrossedRowsneg/prongTrackNeg->GetTPCNclsF();

    if(!prongTrackPos->TestFilterBit(1)) continue;
    if(!prongTrackNeg->TestFilterBit(1)) continue;
    fHistEventV0->Fill(3);
    if(prongTrackPos->Chi2perNDF()>4.)continue;
    if(prongTrackNeg->Chi2perNDF()>4.)continue;
    fHistEventV0->Fill(4);

    //TPC PID
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( prongTrackPos, (AliPID::EParticleType)2))< 3.) goodPiPlusTPC=kTRUE;
    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC( prongTrackNeg, (AliPID::EParticleType)2))< 3.) goodPiMinusTPC=kTRUE;
 
    //se non il TOF ingnora -> usa solo TPC
    bool isTOFPIDok = kFALSE;
    statusTOFPos = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,prongTrackPos);
    statusTOFNeg = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,prongTrackNeg);

    if ( (statusTOFPos ==  AliPIDResponse::kDetPidOk) && (statusTOFNeg ==  AliPIDResponse::kDetPidOk) ) {  
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF( prongTrackPos, (AliPID::EParticleType)2))< 3.) goodPiPlusTOF=kTRUE; //controlla
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF( prongTrackNeg, (AliPID::EParticleType)2))< 3.) goodPiMinusTOF=kTRUE;
      isTOFPIDok=kTRUE;
    }
    //cout << "numero della V0 analizzata  " << i << " is TOF PID ok  " << isTOFPIDok << endl;
    if(isTOFPIDok){
      if(goodPiPlusTPC && goodPiPlusTOF){
	goodPiPlus=kTRUE;
      }      
      if(goodPiMinusTPC && goodPiMinusTOF){
	goodPiMinus=kTRUE;
      }      
    } else {
      if(goodPiPlusTPC){
	goodPiPlus=kTRUE;
      }      
      if(goodPiMinusTPC){
	goodPiMinus=kTRUE;
      }
    }

    //cout << "selezioni su PID figli" << endl;
    if(!goodPiMinus || !goodPiPlus )             continue;
    //cout << goodPiPlus<<endl;
    //cout << goodPiMinus<<endl;
    //cout << TMath::Abs(fPIDResponse->NumberOfSigmasTPC( prongTrackPos, AliPID::kPion)) << endl;
    //cout << TMath::Abs(fPIDResponse->NumberOfSigmasTPC( prongTrackNeg, AliPID::kPion)) << endl;    
    fHistEventV0->Fill(5);

    if(TMath::Abs(v0->EtaProng(pos0or1))>0.8) continue;
    if(TMath::Abs(v0->EtaProng(neg0or1))>0.8) continue;
    fHistEventV0->Fill(6);
    
    //DAUGHTER TRACK QUALITY 
    if(nTPCCrossedRowspos<70) continue;
    if(nTPCCrossedRowsneg<70) continue;
    if(rationCrnFindpos<0.8)  continue;
    if(rationCrnFindneg<0.8)  continue;
    fHistEventV0->Fill(7);    
    
    //V0 cuts
    // if(TMath::Abs(rapidityV0[ParticleType])>ycut[ParticleType])             continue;
    //fHistEventV0->Fill(8);    
    if(TMath::Abs(v0->Eta()) > kEtacut[ParticleType])	   	continue;
    fHistEventV0->Fill(8);    
    
    //    if(v0->PtProng(pos0or1) < .15) continue;
    //if(v0->PtProng(neg0or1) < .15) continue;
    /*
    //load status for PID
    statusPos=prongTrackPos->GetStatus();
    if((statusPos&AliESDtrack::kTPCrefit)==0) continue;
    prongTrackPos->SetAODEvent(fAOD);
    statusNeg=prongTrackNeg->GetStatus();
    if((statusNeg&AliESDtrack::kTPCrefit)==0) continue;
    prongTrackNeg->SetAODEvent(fAOD);
    */
 
   
    //    cout << lPercentiles << endl;
    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	//	cout<< "ho riempito isto masse no tagli " << endl;
	fHistMassvsPt[m]->Fill(fTreeVariableInvMassK0s,fTreeVariablePtV0);
      }
    }
    if(fReadMCTruth){
      if (fMCEvent){
	cout << "\n this particle has passe dall but pt cuts: let's fill the mass Pt histo for true reco K0s "<< endl;
	if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg  && MotherPos->IsPhysicalPrimary()){
	  fHistReconstructedV0PtMass->Fill(fTreeVariableInvMassK0s,v0->Pt(), lPercentiles);
	
	}
      }
    }
    
    if(v0->CosPointingAngle(lBestPrimaryVtxPos) < kMinCosAngle[ParticleType]) 	continue;
    //if(v0->MassK0Short() < .2 || v0->MassK0Short() > .8)   	        continue;
    //    if(v0->DecayLengthV0(lBestPrimaryVtxPos) > kMaxDL)          	continue;
    if(v0->DecayLengthV0(lBestPrimaryVtxPos) < kMinDL[ParticleType])	   	continue;

    double v0Dca = v0->DcaV0ToPrimVertex();
    //if(!fCutCheck){
    if(v0->DcaNegToPrimVertex() < kMinDCAPrimary[ParticleType])      continue;
    if(v0->DcaPosToPrimVertex() < kMinDCAPrimary[ParticleType])      continue;
    if(v0->DcaV0Daughters() > kMaxDCADaughters[ParticleType])        continue;
    if(v0Dca > kMaxDCA[ParticleType]) 	        		     continue;
    if(kctau[ParticleType]>kctauval[ParticleType])                           continue;
    //}
    fHistPtArmvsAlpha->Fill(v0->AlphaV0(),v0->PtArmV0());    
    //if(v0->PtArmV0()< 0.2*TMath::Abs(v0->AlphaV0()))                    continue;
    //fHistPtArmvsAlphaAfterSelection->Fill(v0->AlphaV0(), v0->PtArmV0());    

    fHistEventV0->Fill(9);    
    if(v0->MassK0Short()> 0.55 || v0->MassK0Short()< 0.45) continue;
    bool skipV0=kFALSE;


    Int_t labelPrimOrSecV0=0;
    if(fReadMCTruth){
      if (fMCEvent){
	if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg){
	  if(MotherPos->IsPhysicalPrimary()){
	    cout << "\n this particle has passed selections: label mother pos selected" << labelMotherPos<< endl;
	    fHistSelectedV0PtMass->Fill(fTreeVariableInvMassK0s,fTreeVariablePtV0, lPercentiles);
	  }	    
	}
      }
    }
    
    
    if(!(v0->Pt()> fminPtV0 && v0->Pt()<fmaxPtV0) )continue;
    fHistEventV0->Fill(10);     
    //try6  for(Int_t j=0; j < NumberFirstParticle; j++){
    //try6    if (v0->Pt() >= fEvt->fReconstructedFirst[j].fPt){
    //try6  	skipV0=kTRUE;
    //try6  	//	  continue;
    //try6    }
    //try6  }
    //try6
    if (v0->Pt()>ptTriggerMinimoDati) skipV0=kTRUE;
    if (skipV0){
      cout << " pT V0 > Pt trigger " <<endl; 
      continue;
    }
    
    fHistEventV0->Fill(11);    
    if(fReadMCTruth){
      if (fMCEvent){
	V0PDGCode=0;
	if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg){
	    V0PDGCode=PdgMotherPos;
	  if(MotherPos->IsPhysicalPrimary()){
	    cout <<"selected v0 Pt " <<  v0->Pt() << endl;
	    fHistResolutionV0Pt->Fill(v0->Pt()- MotherPos->Pt(), lPercentiles);
	    fHistResolutionV0Phi->Fill(v0->Phi()- MotherPos->Phi(), lPercentiles);
	    fHistResolutionV0Eta->Fill(v0->Eta()- MotherPos->Eta(), lPercentiles);
	    fHistSelectedV0PtPhi->Fill(v0->Pt(), v0->Phi(), lPercentiles);
	    fHistSelectedV0PtEta->Fill(v0->Pt(), v0->Eta(), lPercentiles);    
	    labelPrimOrSecV0=1;
	  }
	  else if(MotherPos->IsSecondaryFromWeakDecay())      labelPrimOrSecV0=2;
	  else if(MotherPos->IsSecondaryFromMaterial())      labelPrimOrSecV0=3;
	  else labelPrimOrSecV0=4;
	    
	  for (Int_t m =0; m<5;m++){
	    if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	      for(Int_t p=1; p<=4; p++){
		if (labelPrimOrSecV0==p) fHistPrimaryV0[m]->Fill(p, v0->Pt());
	      }
	    }
	  }
	}
      }
    }

    
    //save second particle information (V0)
    //cout << "save second particle information (V0) "<< endl;
    if(!fReadMCTruth){
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sDcaPosV0     = v0->DcaPosToPrimVertex();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sDcaNegV0     = v0->DcaNegToPrimVertex();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sPtArmV0      = v0->PtArmV0();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sAlphaV0      = v0->AlphaV0();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sInvMassK0s   = v0->MassK0Short();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sInvMassLambda   = v0->MassLambda();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sInvMassAntiLambda = v0->MassAntiLambda();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sCosPointingAngle  = v0->CosPointingAngle(lBestPrimaryVtxPos);
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sDcaV0ToPV    = v0->DcaV0ToPrimVertex();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sRap          = rapidityV0[0];
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sPt           = v0->Pt();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sEta          = v0->Eta();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sTheta        = v0->Theta();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sPhi          = v0->Phi();
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].isP           = labelPrimOrSecV0;
    //try6 fEvt->fReconstructedSecond[NumberSecondParticle].sPDGcode      = V0PDGCode;
    

    NumberSecondParticle++;
    }
    fMassV0->Fill(v0->MassK0Short());    

    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	fHistMassvsPt_tagli[m]->Fill(fTreeVariableInvMassK0s,fTreeVariablePtV0);      
      }
    }	   
  } //end loop on tracks for second particle
  
  if(fReadMCTruth){
    if (fMCEvent){
      AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
      if (AODMCTrackArray == NULL){
	return;
	Printf("ERROR: stack not available");
      }
      for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
	Bool_t skipV0_MC=kFALSE;
	AliAODMCParticle* particleV0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
	if (!particleV0) continue;
	if((particleV0->GetPdgCode())!=310) continue;
	if(TMath::Abs(particleV0->Eta())>0.8)continue; //I need to select particles within this eta range!
	if (!(particleV0->IsPhysicalPrimary()))continue; 
	if(!(particleV0->Pt()> fminPtV0 && particleV0->Pt()<fmaxPtV0) )continue;

//try6	for(Int_t j=0; j < NumberFirstParticle; j++){
//try6	  if (particleV0->Pt() >= fEvt->fReconstructedFirst[j].fPt){
//try6	    skipV0_MC=kTRUE;
//try6	  }
//try6	}
	  if (particleV0->Pt()>ptTriggerMinimoMC) skipV0_MC=kTRUE;
	if (skipV0_MC)      continue;
   	 
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sDcaPosV0     = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sDcaNegV0     = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sPtArmV0      = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sAlphaV0      = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sInvMassK0s   = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sInvMassLambda   = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sInvMassAntiLambda = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sCosPointingAngle  = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sDcaV0ToPV    = 0;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sRap          = particleV0->Y();
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sPt           = particleV0->Pt();
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sEta          = particleV0->Eta();
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sTheta        = particleV0->Theta();
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sPhi          = particleV0->Phi();
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].isP           = 1;
	//try6 fEvt->fReconstructedSecond[NumberSecondParticleMC].sPDGcode      = particleV0->GetPdgCode();
	//try6   
	NumberSecondParticleMC++;
      }
    }
  }

  cout <<"candidati second particle nell'evento analizzato " <<  NumberSecondParticle << endl;      
  cout <<"candidati second particle nell'evento analizzato, MC truth " <<  NumberSecondParticleMC << endl;      
  fHistEventV0->AddBinContent(12, NumberSecondParticleMC);    
  if (!fReadMCTRuth) NumberSecondParticleAll=NumberSecondParticle;
  if (fReadMCTRuth) NumberSecondParticleAll=NumberSecondParticleMC;
  if(NumberSecondParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    cout << "event has trigger particle but no V0 " << endl;
    PostData(3,fBkgTree);
    return;
  }

  FifoShiftok=kTRUE;
  fHistEventMult->Fill(14);  
  if(!fReadMCTruth){
//try6  fEvt->fNumberCandidateFirst = NumberFirstParticle;
//try6  fEvt->fNumberCandidateSecond = NumberSecondParticle;
  }
  else{
//try6  fEvt->fNumberCandidateFirst = NumberFirstParticleMC;
//try6  fEvt->fNumberCandidateSecond = NumberSecondParticleMC;
  }

  //Fill histos about selected events (wt >0 Trigger part. and >0 V0)  
  fHistTriggerwV0->Fill(NumberFirstParticle);
  fHistTriggerwV0MCTruth->Fill(NumberFirstParticleMC);
  fHistMultvsV0->Fill(NumberSecondParticleAll,lPercentiles);
  fHist_multiplicity->Fill(lPercentiles);
  fHistZvertex->Fill(lBestPrimaryVtxPos[2]);
  fHistTrack->AddBinContent(10, NumberFirstParticle);
  fHistTrack->AddBinContent(11, NumberFirstParticleMC);
  fHistMultvsTrigger->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruth->Fill(NumberFirstParticleMC, lPercentiles);
  fHistMultiplicityVsVertexZ->Fill(lBestPrimaryVtxPos[2], lPercentiles);


  for(Int_t i=0; i< 100; i++){
    if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMult->AddBinContent(i+1,NumberFirstParticle);  
  }


  for(Int_t l=0; l< NumberFirstParticle; l++){
//try6     fHistPt->Fill(fEvt->fReconstructedFirst[l].fPt);
//try6     fHistPtvsMult->Fill(fEvt->fReconstructedFirst[l].fPt, lPercentiles);
//try6     fHist_eta_phi->Fill(fEvt->fReconstructedFirst[l].fPhi, fEvt->fReconstructedFirst[l].fEta);
  }
  
  /*
    if(fFirstpart == fSecondpart){ // particella carica sempre diversa da K0Short
    DoPairshh(lPercentiles, fieldsign);  
    else{
    //Remove candidates that are at the same time a ptc1 and ptc2
    for (int i=0; i < fEvt->fNumberCandidateFirst; i++) {
    for (int j=0; j<fEvt->fNumberCandidateSecond; j++) {
    if (fEvt->fReconstructedFirst[i].index == fEvt->fReconstructedSecond[j].index) {
    //cout<<"the track can be both tracks!"<<endl;
    fEvt->fReconstructedFirst[i].doSkipOver = kTRUE;
    fEvt->fReconstructedSecond[j].doSkipOver = kTRUE;
    }
    }
    }
  */
  //--------------------------------------------------------------
  //  DoPairsh1h2((Int_t)lPercentiles, fieldsign, lBestPrimaryVtxPos[2]);  
     
  PostData(1, fOutputList);     
  PostData(2,fSignalTree);
  PostData(3,fBkgTree);
}




//----------------------------------------------------------------------------------------------------
/*
void AliAnalysisTaskMyTask::DoPairsh1h2 ( const Float_t lPercentiles, int fieldsign, Double_t lBestPrimaryVtxPos)  {

  // UInt_t dimsparse;
  // if(!fReadMCTruth) dimsparse=16;
  // else dimsparse=25;
  // Double_t xsparseSignal[dimsparse];
  // for(UInt_t counter = 0; counter < dimsparse ; counter++) xsparseSignal[counter]=-999;//Default position for THnSparse 
  // Double_t xsparseBkg[dimsparse];
  // for(UInt_t counter = 0; counter < dimsparse ; counter++) xsparseBkg[counter]=-999;//Default position for THnSparse 

  //-----------
  double DCAxyP1  = -999. ;
  double DCAzP1   = -999. ;  
  double DCAxyP2  = -999. ; 
  double DCAzP2   = -999. ;  

  double  ptP1 = -999.;
  double  ptP2 = -999.;

  // Short_t chargeP1 = -999.;
  // Short_t chargeP2 = -999.;

  bool isP1  = kFALSE;
  bool isaP1 = kFALSE;
  bool isP2  = kFALSE;
  bool isaP2 = kFALSE;

  Int_t  SignP1 = -999;
  Int_t  SignP2 = -999;

  double dphi  = -999.;
  double dphis = -999.;
  double deta  = -999.;
  //  double detas = -999.;
  double dtheta = -999.;
  
  bool isMC1 = kFALSE;
  bool isMC2 = kFALSE;

  bool isMCvector = kFALSE;
  bool sameMother = kFALSE;
  bool sameGrandMother = kFALSE;
  
  Int_t mcMotherLabelP1 = -999;
  Int_t mcMotherLabelP2 = -999;

  Int_t mcGrandMotherLabelP1 = -999;
  Int_t mcGrandMotherLabelP2 = -999;
 
  Int_t typeP1 = -999;
  Int_t typeP2 = -999;
 
  Int_t mcPDGMotherP1 = 0;
  Int_t mcPDGMotherP2 = 0;
  //  Int_t mcMotherBin = 0;
 
  Int_t mcPDGGrandMother = 0;
  Int_t mcGrandMotherBin = 0;

  Int_t mcPDGcodeP1 = 0;
  Int_t mcPDGcodeP2 = 0;
  // Int_t mcPDG1Bin = 0;
  // Int_t mcPDG2Bin = 0;

  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
 
  double pairKstar   = 0.;
  double pairKstarMC = 0.;
  double pairMass  = 0.;
  double pairMassE = 0.;
  double pairKt    = 0.;

  cout << "fNumberCandidateFirst " << fEvt->fNumberCandidateFirst << "  fNumberCandidateSecond " <<  fEvt->fNumberCandidateSecond << endl; 
  for (int i=0; i<fEvt->fNumberCandidateFirst; i++) {
    
    //ch    if (fEvt->fReconstructedFirst[i].doSkipOver) continue;
    
    //c DCAxyP1         = fEvt->fReconstructedFirst[i].fDCAxy;
    //c DCAzP1          = fEvt->fReconstructedFirst[i].fDCAz;
    //c ptP1            = fEvt->fReconstructedFirst[i].fPt;
    //c if(ptP1 < fminPtj)continue; //superfluo?
    //c //    chargeP1        = fEvt->fReconstructedFirst[i].fCharge;
    //c     isMC1           = fEvt->fReconstructedFirst[i].isMCptc;
    //c     mcMotherLabelP1 = fEvt->fReconstructedFirst[i].fMCmumIdx;
    //c     typeP1          = fEvt->fReconstructedFirst[i].fMCcode;
    //c     //  if(typeP1 < 0) continue;
    //c  mcGrandMotherLabelP1 = fEvt->fReconstructedFirst[i].fMCgrandmumIdx;
    //c  mcPDGcodeP1          = fEvt->fReconstructedFirst[i].fPDGcode;

    */ //for try7********************************************************************************************+

    /*
      if(mcPDGcodeP1 == 11 || mcPDGcodeP1 == -11)                    mcPDG1Bin=1;
      else if(mcPDGcodeP1 == 13 || mcPDGcodeP1 == -13)               mcPDG1Bin=2;
      else if(mcPDGcodeP1 == 211 || mcPDGcodeP1 == -211)             mcPDG1Bin=3;
      else if(mcPDGcodeP1 == 321 || mcPDGcodeP1 == -321)             mcPDG1Bin=4;
      else if(mcPDGcodeP1 == 2212 || mcPDGcodeP1 == -2212)           mcPDG1Bin=5;
      else if(mcPDGcodeP1 == 1000010020|| mcPDGcodeP1 == -1000010020)mcPDG1Bin=6;
    */
    
    //c isP1  = fEvt->fReconstructedFirst[i].isP;
    //c isaP1 = fEvt->fReconstructedFirst[i].isaP;
    //c 
    //c if(isP1) SignP1 = 1;
    //c else if (isaP1) SignP1 = -1;

    //c mcPDGMotherP1    = fEvt->fReconstructedFirst[i].fMCmumPDG;
    //c mcPDGGrandMother = fEvt->fReconstructedFirst[i].fMCgrandmumPDG;
        
/* //for try7**************************************************************************************************
for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      // For same event pairs
      //cout << "numero dell'evento con cui mixo" <<  eventNumber << " "<<((fEvt+eventNumber)->fNumberCandidateSecond)<< endl;
      //if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateFirst)!=0.) evmultmixed++; 
      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateSecond)!=0.) {
	cout << "con questo evento faccio mixing " << endl;
	cout << evmultmixed << endl;
	evmultmixed++; 
      }
      for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateSecond; j++) {
	//	cout << " j " << j <<  endl;
	//        cout<<" event number "<<eventNumber<<endl;
	//	cout<<"eventNumber: "<<eventNumber<<" j: "<<j<<" i: "<<i<<" fEvt+eventNumber: "<<(fEvt+eventNumber)<<endl;
	
	//c 	if ((fEvt+eventNumber)->fReconstructedSecond[j].doSkipOver) continue;
         
	//c DCAxyP2         = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAxy;
	//c DCAzP2          = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAz;
	//c ptP2            = (fEvt+eventNumber)->fReconstructedSecond[j].sPt;
	//c if(ptP2< fminPtV0)continue; //superfluo?
	//c //chargeP2        = (fEvt+eventNumber)->fReconstructedSecond[j].sCharge;
	//c isMC2           = (fEvt+eventNumber)->fReconstructedSecond[j].isMCptc;
	//c mcMotherLabelP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumIdx;
	//c typeP2          = (fEvt+eventNumber)->fReconstructedSecond[j].sMCcode;
	//	if(typeP2 < 0) continue;
	//c mcGrandMotherLabelP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCgrandmumIdx;
	//c mcPDGcodeP2     = (fEvt+eventNumber)->fReconstructedSecond[j].sPDGcode;

	/*
	  if(mcPDGcodeP2 == 11 || mcPDGcodeP2 == -11)                    mcPDG2Bin=1;
	  else if(mcPDGcodeP2 == 13 || mcPDGcodeP2 == -13)               mcPDG2Bin=2;
	  else if(mcPDGcodeP2 == 211 || mcPDGcodeP2 == -211)             mcPDG2Bin=3;
	  else if(mcPDGcodeP2 == 321 || mcPDGcodeP2 == -321)             mcPDG2Bin=4;
	  else if(mcPDGcodeP2 == 2212 || mcPDGcodeP2 == -2212)           mcPDG2Bin=5;
	  else if(mcPDGcodeP2 == 1000010020|| mcPDGcodeP2 == -1000010020)mcPDG2Bin=6;
	*/

	//c	isP2  = (fEvt+eventNumber)->fReconstructedSecond[j].isP;
	//c        isaP2 = (fEvt+eventNumber)->fReconstructedSecond[j].isaP;
	//c	mcPDGMotherP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumPDG;

	//c if(isP2) SignP2 = 1;
	//c else if (isaP2) SignP2 = -1;
	//c 
	//c if(isMC1 && isMC2) isMCvector = kTRUE;
	//c else isMCvector = kFALSE;
	//c 
	//c if(mcMotherLabelP1 == mcMotherLabelP2 && mcMotherLabelP1!=-999)sameMother = kTRUE;
	//c else sameMother = kFALSE;
	
	// if(eventNumber == 0){
	//   cout<<"-- OUTSIDE"<<endl;
	//   cout<<"---------------------------------------------------> Mum PDG1 "<<fEvt->fReconstructedFirst[i].fMCmumPDG<<" PDG2 "<<(fEvt+eventNumber)->fReconstructedSecond[j].sMCmumPDG<<endl;
	//   cout<<"                                                     MC label 1: "<<mcMotherLabelP1<<" 2: "<<mcMotherLabelP2<<endl;
	// }
	
	//c if(isMCvector && sameMother){
	  
	//c if((fFirstpart == 3 && fSecondpart == 4) || (fFirstpart == 4 && fSecondpart == 3)){
	    
	// // // cout<<"---------------------------------------------------------------------------------------------> Mum PDG: "<<mcPDGMother<<endl;
	// if(eventNumber == 0){
	//   cout<<"-- OUTSIDE"<<endl;
	//   cout<<"---------------------------------------------------> Mum PDG1 "<<fEvt->fReconstructedFirst[i].fMCmumPDG<<" PDG2 "<<(fEvt+eventNumber)->fReconstructedSecond[j].sMCmumPDG<<endl;
	//   cout<<"                                                     MC label 1: "<<mcMotherLabelP1<<" 2: "<<mcMotherLabelP2<<endl;
	// }
	    
	//   if(mcPDGMother == 3124 || mcPDGMother == -3124 )  mcMotherBin = 1;  //1520
	//   /*
	//     else if(mcPDGMother == 3126 || mcPDGMother == -3126 )  mcMotherBin = 2;  //1820
	//     else if(mcPDGMother == 3128 || mcPDGMother == -3128 )  mcMotherBin = 3;  //2100
	//     else if(mcPDGMother == 13122|| mcPDGMother == -13122)  mcMotherBin = 4;  //1405
	//     else if(mcPDGMother == 13124|| mcPDGMother == -13124)  mcMotherBin = 5;  //1690
	//     else if(mcPDGMother == 13126|| mcPDGMother == -13124)  mcMotherBin = 6;  //1830
	//     else if(mcPDGMother == 23122|| mcPDGMother == -23122)  mcMotherBin = 7;  //1600
	//     else if(mcPDGMother == 23124|| mcPDGMother == -23124)  mcMotherBin = 8;  //1890
	//     else if(mcPDGMother == 23126|| mcPDGMother == -23126)  mcMotherBin = 9;  //2110
	//     else if(mcPDGMother == 33122|| mcPDGMother == -33122)  mcMotherBin = 10; //1670
	//     else if(mcPDGMother == 43122|| mcPDGMother == -43122)  mcMotherBin = 11; //1800
	//     else if(mcPDGMother == 53553|| mcPDGMother == -53553)  mcMotherBin = 12; //1810
	//
	//
	//   */
	//   
	//   else if(TMath::Abs(mcPDGMother)==1)  mcMotherBin = 2;  //quark d
	//   else if(TMath::Abs(mcPDGMother)==2)  mcMotherBin = 3;  //quark u
	//   else if(TMath::Abs(mcPDGMother)==3)  mcMotherBin = 4;  //quark s
	//   else if(TMath::Abs(mcPDGMother)==4)  mcMotherBin = 5;  //quark b
	//   else if(TMath::Abs(mcPDGMother)==5)  mcMotherBin = 6;  //quark c
	//   else if(TMath::Abs(mcPDGMother)==6)  mcMotherBin = 7;  //quark t
	//   else if(mcPDGMother == 21 || mcPDGMother == -21 )  mcMotherBin = 8;  //gluoni
	//   else if(mcPDGMother == 130|| mcPDGMother == -130)  mcMotherBin = 9;  //KoL1405
	//   else if(mcPDGMother == 313|| mcPDGMother == -313)  mcMotherBin = 10;  //K*892
	//
	//   else if(mcPDGMother == 4122 || mcPDGMother == -4122)   mcMotherBin = 13; //Lambda_c
	//   else if(mcPDGMother!=0) {
	//     // cout<<"--------------------------------------------------------------------------------> "<<mcPDGMother<<endl;
	//     mcMotherBin = 14; 
	//   }//
	//   //	    hMotherID->Fill(mcMotherBin);
	//   

	//c }
	//c }
	//c else if(isMCvector && ((fFirstpart == 3 && fSecondpart == 4) || (fFirstpart == 4 && fSecondpart == 3))){
	    
	//c if(mcGrandMotherLabelP1 == mcGrandMotherLabelP2){sameGrandMother = kTRUE;}
	  
	//  cout<<"GM-------------------------------- chargeP1: "<<chargeP1<<" ----------- chargeP2: "<<chargeP2<<" ---------------------> "<<mcPDGGrandMother<<endl;
	  
	//c f(TMath::Abs(mcPDGGrandMother)>= 1 && TMath::Abs(mcPDGGrandMother)<= 6)  mcGrandMotherBin = 1;  //quark
	//c lse if(TMath::Abs(mcPDGGrandMother)==2212)  mcGrandMotherBin = 2;  //p
	//c lse if(TMath::Abs(mcPDGGrandMother)==21)  mcGrandMotherBin = 3;    //g
	//c lse if(TMath::Abs(mcPDGGrandMother)> 400 && TMath::Abs(mcPDGGrandMother)< 500 )  mcGrandMotherBin = 5;  //D meson
	//c else if(mcPDGGrandMother!=0) {
	// cout<<"--------------------------------------------------------------------------------> "<<mcPDGMother<<endl;
	//c mcGrandMotherBin  = 6; 
	//c }//
	//}
	// else if(!sameMother)
	//   mcMotherBin = 15; 
	
	//cout<<"----------------- outside: "<<mcGrandMotherBin<<endl;

	// Pair histogramming
/* //for try7***********************************************************
	deta   = CalculateDeltaEta(fEvt->fReconstructedFirst[i].fEta, (fEvt+eventNumber)->fReconstructedSecond[j].sEta);
	dtheta = CalculateDeltaTheta(fEvt->fReconstructedFirst[i].fTheta, (fEvt+eventNumber)->fReconstructedSecond[j].sTheta);
	//dphi   = CalculateDeltaPhi(fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi);
	dphi = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,0); // 2 - 1
	
	//dphis = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,fRadius); // 2 - 1

	
	// // Apply two-track cuts //RAMONA
	//c if (fkApplyTtc) {
	//	  if (TMath::Abs(dphis)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
	//c if (TMath::Abs(dphi)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
	//c }

	

	if (eventNumber==0) {//Same event pair histogramming

	  fTreeVariablePtTrigger              = fEvt->fReconstructedFirst[i].fPt;		       
	  fTreeVariableChargeTrigger          = fEvt->fReconstructedFirst[i].fCharge;		       
	  fTreeVariableEtaTrigger             = fEvt->fReconstructedFirst[i].fEta; 		       
	  fTreeVariablePhiTrigger	      =	fEvt->fReconstructedFirst[i].fPhi;             
	  fTreeVariableDCAz		      =	fEvt->fReconstructedFirst[i].fDCAz;             
	  fTreeVariableDCAxy		      =	fEvt->fReconstructedFirst[i].fDCAxy;              
	  fTreeVariableRapK0Short	      =	fEvt->fReconstructedSecond[j].sRap;     	      
	  fTreeVariableDcaV0ToPrimVertex      = fEvt->fReconstructedSecond[j].sDcaV0ToPV;	       	      
	  fTreeVariableDcaPosToPrimVertex     = fEvt->fReconstructedSecond[j].sDcaPosV0; 	       	      
	  fTreeVariableDcaNegToPrimVertex     = fEvt->fReconstructedSecond[j].sDcaNegV0;  	       	      
	  fTreeVariableV0CosineOfPointingAngle= fEvt->fReconstructedSecond[j].sCosPointingAngle;      
	  fTreeVariablePtV0		      = fEvt->fReconstructedSecond[j].sPt;     
	  fTreeVariableInvMassK0s	      =	fEvt->fReconstructedSecond[j].sInvMassK0s;   
	  fTreeVariableInvMassLambda	      =	fEvt->fReconstructedSecond[j].sInvMassLambda;   
	  fTreeVariableInvMassAntiLambda      =	fEvt->fReconstructedSecond[j].sInvMassAntiLambda;   
	  fTreeVariableEtaV0		      =	fEvt->fReconstructedSecond[j].sEta;    
	  fTreeVariablePhiV0		      =	fEvt->fReconstructedSecond[j].sPhi;   
	  fTreeVariablePtArmenteros           = fEvt->fReconstructedSecond[j].sPtArmV0;     
	  fTreeVariableAlpha	              = fEvt->fReconstructedSecond[j].sAlphaV0;  
	  fTreeVariableDeltaEta	      	      = deta;  
	  fTreeVariableDeltaPhi		      = dphi;
	  fTreeVariableDeltaTheta             = dtheta;
	  fTreeVariableMultiplicity	      = lPercentiles;
	  fTreeVariableZvertex                = lBestPrimaryVtxPos;
	  fTreeVariablePDGCode =   fEvt->fReconstructedSecond[j].sPDGcode;
	  fTreeVariableisPrimaryTrigger       =  fEvt->fReconstructedFirst[i].isP ;
	  fTreeVariableisPrimaryV0            =  fEvt->fReconstructedSecond[j].isP ;
	  /*
	    if(fReadMCTruth == kTRUE){
	    tMCtruepair    =  isMCvector;		
	    tMCSameMother  =  sameMother;		
	    tMCMotherP1    =  mcPDGMotherP1;//mcMotherBin;		
	    tMCMotherP2    =  mcPDGMotherP2;//mcMotherBin;		
	    tMCptcTypeP1   =  typeP1 ;		
	    tMCptcTypeP2   =  typeP2 ;		
	    tMCSameGM	   =  sameGrandMother;	
	    tMotherPDG	   =  mcGrandMotherBin;	
	    tpdgcodeP1 	   =  mcPDGcodeP1;//mcPDG1Bin;		
	    tpdgcodeP2	   =  mcPDGcodeP2;//mcPDG2Bin;		 
	    tKstarGen      =  pairKstarMC;            
	  */

	  // // if((TMath::Abs(mcPDGcodeP1 == 321) && TMath::Abs(mcPDGcodeP2 == 2212))){
	  // cout<<"-- INSIDE"<<endl;
	  // cout<<"---------------------------------------------------> Mum PDG1 "<< mcPDGMother <<" PDG2 - "<<endl;
	  // cout<<"                                                     MC label 1: "<<mcMotherLabelP1<<" 2: "<<mcMotherLabelP2<<endl;
	  // //}
	  //c}
	  /* //for try7***********************************************************
	  fSignalTree->Fill();  
	  //fHistSparseSignal->Fill(xsparseSignal);
	}

	else {//Mixed-event pair histogramming
	  //  cout << "ciao!" << endl;


	  fTreeVariablePtTrigger              = fEvt->fReconstructedFirst[i].fPt;		       
	  fTreeVariableChargeTrigger          = fEvt->fReconstructedFirst[i].fCharge;		       
	  fTreeVariableEtaTrigger             = fEvt->fReconstructedFirst[i].fEta; 		       
	  fTreeVariablePhiTrigger	      =	fEvt->fReconstructedFirst[i].fPhi;             
	  fTreeVariableDCAz		      =	fEvt->fReconstructedFirst[i].fDCAz;             
	  fTreeVariableDCAxy		      =	fEvt->fReconstructedFirst[i].fDCAxy;              
	  fTreeVariableRapK0Short	      =	fEvt->fReconstructedSecond[j].sRap;     	      
	  fTreeVariableDcaV0ToPrimVertex      = fEvt->fReconstructedSecond[j].sDcaV0ToPV;	       	      
	  fTreeVariableDcaPosToPrimVertex     = fEvt->fReconstructedSecond[j].sDcaPosV0; 	       	      
	  fTreeVariableDcaNegToPrimVertex     = fEvt->fReconstructedSecond[j].sDcaNegV0;  	       	      
	  fTreeVariableV0CosineOfPointingAngle= fEvt->fReconstructedSecond[j].sCosPointingAngle;      
	  fTreeVariablePtV0		      = fEvt->fReconstructedSecond[j].sPt;     
	  fTreeVariableInvMassK0s	      =	fEvt->fReconstructedSecond[j].sInvMassK0s;   
	  fTreeVariableInvMassAntiLambda      =	fEvt->fReconstructedSecond[j].sInvMassAntiLambda;   
	  fTreeVariableInvMassLambda	      =	fEvt->fReconstructedSecond[j].sInvMassLambda;   
	  fTreeVariableEtaV0		      =	fEvt->fReconstructedSecond[j].sEta;    
	  fTreeVariablePhiV0		      =	fEvt->fReconstructedSecond[j].sPhi;   
	  fTreeVariablePtArmenteros           = fEvt->fReconstructedSecond[j].sPtArmV0;     
	  fTreeVariableAlpha	              = fEvt->fReconstructedSecond[j].sAlphaV0;  
	  fTreeVariableDeltaEta	       	      =deta;  
	  fTreeVariableDeltaPhi		      =dphi;
	  fTreeVariableDeltaTheta             =dtheta;      
	  fTreeVariableMultiplicity	      = lPercentiles;
	  fTreeVariableZvertex                = lBestPrimaryVtxPos;
	  fTreeVariablePDGCode =   fEvt->fReconstructedSecond[j].sPDGcode;
	  fTreeVariableisPrimaryTrigger       =  fEvt->fReconstructedFirst[i].isP;
	  fTreeVariableisPrimaryV0            =  fEvt->fReconstructedSecond[j].isP;

	  /*
	    if(fReadMCTruth == kTRUE){
	    tMCtruepair    =  isMCvector;		
	    tMCSameMother  =  sameMother;		
	    tMCMotherP1    =  mcPDGMotherP1;//mcMotherBin;		
	    tMCMotherP2    =  mcPDGMotherP2;//mcMotherBin;		
	    tMCptcTypeP1   =  typeP1 ;		
	    tMCptcTypeP2   =  typeP2 ;		
	    tMCSameGM	   =  sameGrandMother;	
	    tMotherPDG	   =  mcGrandMotherBin;	
	    tpdgcodeP1 	   =  mcPDGcodeP1;//mcPDG1Bin;		
	    tpdgcodeP2	   =  mcPDGcodeP2;//mcPDG2Bin;		 
	    tKstarGen      =  pairKstarMC;            
	    }
	  */ 
	  /* //for try7***********************************************************
	  fBkgTree->Fill();  
	  
	} //mixed
	
      } // second part

    }//end event loop

    if (evmultmixed!=0) multmixedcounted = kTRUE;
    
  } // first part
  
  if  (multmixedcounted) 
    fHistMultiplicityOfMixedEvent->Fill(evmultmixed); //mi dice con quanti eventi ho mixato 
  
}

/*
//----------------------------------------------------------------------------------------------
//void AliAnalysisTaskKPFemto::DoPairshh (const Float_t lPercentiles, int fieldsign) {
void AliAnalysisTaskKPFemto::DoPairshh (const Int_t lPercentiles, int fieldsign) {
return;
}

//-----------------------------------------------------------------------------------------------

double AliAnalysisTaskMyTask::CalculateMass(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

// Calculate Invariant Mass 
TLorentzVector  vP1,vP2,vSum;
  
vP1.SetXYZM(momentum1[0],momentum1[1],momentum1[2],mass1); 
vP2.SetXYZM(momentum2[0],momentum2[1],momentum2[2],mass2);       
vSum=vP1+vP2;

double mass = vSum.M();
return mass;
}
*/


//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskMyTask::CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad) { //AliFemtoUser/AliFemtoPairCutDetaDphi.h

  const Double_t unit_factor = 0.299792458 / 2.0;
  const Double_t b_field = 0.5006670488586 * magSign;

  Double_t  shift1 = TMath::ASin(unit_factor * chg1 * b_field * rad / ptv1);
  Double_t  shift2 = TMath::ASin(unit_factor * chg2 * b_field * rad / ptv2);

  double dps = (phi1 + shift1) - (phi2 + shift2);
  
  //  dps = TVector2::Phi_mpi_pi(dps); //to be checked

  return dps; //deltaphi;
  
}
//_______________________________________________________________


//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskMyTask::CalculateDeltaEta( Double_t eta1, Double_t eta2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double deta = eta2 - eta1;
  return deta;
}
//_______________________________________________________________
Double_t AliAnalysisTaskMyTask::CalculateDeltaTheta( Double_t theta1, Double_t theta2 )  {  
  const double dtheta = theta2 - theta1;
  return dtheta;
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskMyTask::CalculateDeltaPhi( Double_t phi1, Double_t phi2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double dphi = phi2 - phi1;
  return dphi;
}

/*
//----------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskMyTask::ThetaS( Double_t posSftR125[3] ) const { // Hans B
// Returns the longitudinal angle of the particles propagated
// position at R=1.25m. See
// https://edms.cern.ch/file/406391/2/ALICE-INT-2003-038.pdf
// for the ALICE coordinate system. Theta is zero at positive z,
// pi/2 at z = 0 aka the xy plane and pi at negative z 

// R^    ^  
//  |   /
//  |?'/
//  | / ?
//  |/____>z
// 
// Let's compute ?' and ? = pi/2 - ?'
// where ?' can even be and should 
// sometimes be negative
// tan(?') = z/R
// ?' = arctan(z/R)
// ? = pi/2 - ?'
//   = pi/2 - arctan(z/R)
// Note that in the doc above theta
// is calculated as arccos(z/sqrt(x^2+y^2+z^2))

// Array of positions is 85,105,125,..cm,
// we take the z position at R=1.25m
// return TMath::Pi()/2. - TMath::ATan(fXshifted[2][2]/125.);
return TMath::Pi()/2. - TMath::ATan(posSftR125[2]/125.); // ok here R is really there --> transverse plane 
}
//_______________________________________________________________
Double_t AliAnalysisTaskMyTask::EtaS( Double_t posSftR125[3] ) const {  // Hans B
// Returns the corresponding eta of a pri. part. 
// with this particles pos at R=1.25m

// http://en.wikipedia.org/wiki/Pseudorapidity
// ? = -ln[ tan(?/2)]
// printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
// ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));
return -TMath::Log( TMath::Tan(ThetaS(posSftR125 )/2.) );
}
*/
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
