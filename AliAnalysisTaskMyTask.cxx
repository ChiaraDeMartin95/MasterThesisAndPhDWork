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
#include "AliAnalysisKPEventCollectionChiara.h"
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
  fMultSelection(0),
fEventCuts(0), 			  			
  fOutputList(0), 
  fSignalTree(0), 
  fBkgTree(0), 
  fOutputList2(0),
  fOutputList3(0),  
  fMCEvent(0), 
  fReadMCTruth(0),
  isEfficiency(0),
 fEventColl(0x0), 
fEvt(0x0), 
  fzVertexBins(10), 
  fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(50),
  fnEventsToMix(50),
  fHistPt(0), 
  fHistDCAxym1(0),
  fHistDCAzm1(0),
  fHistDCAxym2(0),
  fHistDCAzm2(0),
  fHistPtV0(0), 
  fHistPtTMin(0),
  fHistPtTMinMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistZvertex(0),  
  fHist_eta_phi(0),  
						fHist_multiplicity(0),
						fHistEventMult(0), 
  fHistEventV0(0), 
  fHistTrack(0), 
  fHistPDG(0), 
  fHistSecondParticleAll(0),
  fHistSecondParticleTruthAll(0),
  fHistSecondParticle(0),
  fHistSecondParticleTruth(0),
  fMassV0(0), 
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0),
  fHistMultvsV0(0),
  fHistMultvsV0Truth(0),
  fHistMultvsV0MC(0),
  fHistMassvsPt(0),
  fHistMassvsPt_tagli(0x0),
  fHistMultvsTriggerBefAll(0),
  fHistMultvsTriggerMCTruthBefAll(0),
  fHistMultvsTriggerAll(0),
  fHistMultvsTriggerMCTruthAll(0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistMassPhoton(0),
  fHistMass2Photon(0),
  fHistPtArmvsAlpha(0),
  fHistPtArmvsAlphaAfterSelection(0),
  fHistPtArmvsAlphaAfterPhotonSelection(0),
  fHistPtArmvsAlphaAfterLambdaRejectionSelection(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistTriggervsMultMC(0),
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
  fTreeVariablectau(0),			      
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
								 								 fMultSelection(0), 			  		
							 						 fEventCuts(0),								 
								 fOutputList(0), 
								 fSignalTree(0), 
								 fBkgTree(0), 
								 fOutputList2(0), 
								 fOutputList3(0), 
								 fMCEvent(0), 
								 fReadMCTruth(0),
								 isEfficiency(0),
								 fEventColl(0x0), 
								 fEvt(0x0), 
								 fzVertexBins(10), 
								 fnMultBins(20),	 
								 fMaxFirstMult(50),
								 fMaxSecondMult(50),
								 fnEventsToMix(50),
								 fHistPt(0), 
  fHistDCAxym1(0),
  fHistDCAzm1(0),
  fHistDCAxym2(0),
  fHistDCAzm2(0),
  fHistPtV0(0), 
  fHistPtTMin(0),
  fHistPtTMinMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistZvertex(0),  
  fHist_eta_phi(0),  
  fHist_multiplicity(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistTrack(0), 
  fHistPDG(0), 
  fHistSecondParticleAll(0),
  fHistSecondParticleTruthAll(0),
  fHistSecondParticle(0),
  fHistSecondParticleTruth(0),
  fMassV0(0),
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0), 
  fHistMultvsV0(0),
  fHistMultvsV0Truth(0),
  fHistMultvsV0MC(0),
  fHistMassvsPt(0),
  fHistMassvsPt_tagli(0x0),
  fHistMultvsTriggerBefAll(0),
  fHistMultvsTriggerMCTruthBefAll(0),
  fHistMultvsTriggerAll(0),
  fHistMultvsTriggerMCTruthAll(0),
  fHistMultvsTrigger(0),
  fHistMultvsTriggerMCTruth(0),
  fHistMassPhoton(0),
  fHistMass2Photon(0),
  fHistPtArmvsAlpha(0),
  fHistPtArmvsAlphaAfterSelection(0),
  fHistPtArmvsAlphaAfterPhotonSelection(0),
  fHistPtArmvsAlphaAfterLambdaRejectionSelection(0),
  fHistTrigger(0),
  fHistTriggerMCTruth(0),
  fHistTriggerwV0(0),
  fHistTriggerwV0MCTruth(0),
  fHistMultiplicityVsVertexZ(0),
  fHistTriggervsMult(0),
  fHistTriggervsMultMC(0),
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
  fTreeVariablectau(0),			      
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
  DefineOutput(4, TList::Class());  
  DefineOutput(5, TList::Class());  
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
  if(fOutputList2) {
    delete fOutputList2;     // at the end of your task, it is deleted from memory by calling this function
  }
  if(fOutputList3) {
    delete fOutputList3;     // at the end of your task, it is deleted from memory by calling this function
  }
  if (farrGT)
    delete[] farrGT;
  farrGT=0;

  // if (fHistMassvsPt)
  // delete[] fHistMassvsPt;

  // if (fHistMassvsPt_tagli)
  // delete[] fHistMassvsPt_tagli ;

  
  for(unsigned short i=0; i < fzVertexBins; i++){
    for(unsigned short j=0; j < fnMultBins; j++){
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }
  delete[] fEventColl;
  
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
      
      if(isV0==kFALSE){
	fHistPDG->Fill(particle->GetPdgCode());      
	if((particle->Charge())==0) continue;	
	if(TMath::Abs(particle->Eta())>0.8)continue; //I need to select particles within this eta range!
	if (!(particle->IsPhysicalPrimary()))continue; 
	fHistGeneratedTriggerPtPhi->Fill(particle->Pt(), particle->Phi(), lPercentiles);
	fHistGeneratedTriggerPtEta->Fill(particle->Pt(), particle->Eta(), lPercentiles);
      }
      if(isV0==kTRUE){ 
	if ((particle->GetPdgCode())!=310) continue;
	if(TMath::Abs(particle->Eta())>0.8)continue;
	if (!(particle->IsPhysicalPrimary()))continue;
	//	cout << "these are all K0s generated passing selection criteria: label K0s " << particle->Label()<<endl; 
	fHistGeneratedV0PtPhi[0]->Fill(particle->Pt(),particle->Phi(), lPercentiles );
	fHistGeneratedV0PtEta[0]->Fill(particle->Pt(),particle->Eta(), lPercentiles );
	if (TMath::Abs(particle->Y())>0.5 )continue;
	//	cout << "these are all K0s generated passing selection criteria: label K0s " << particle->Label()<<endl; 
	fHistGeneratedV0PtPhi[1]->Fill(particle->Pt(),particle->Phi(), lPercentiles );
	fHistGeneratedV0PtEta[1]->Fill(particle->Pt(),particle->Eta(), lPercentiles );
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
      // fHistSelectedTriggerPtPhi[0]->Fill(track->Pt(), track->Phi(), lPercentiles);
      // fHistSelectedTriggerPtEta[0]->Fill(track->Pt(), track->Eta(), lPercentiles);    
      if(  (TMath::Abs(track->ZAtDCA()) < 2.)) {
	fHistSelectedTriggerPtPhi[1]->Fill(track->Pt(), track->Phi(), lPercentiles);
	fHistSelectedTriggerPtEta[1]->Fill(track->Pt(), track->Eta(), lPercentiles);    
      }
      if(  (TMath::Abs(track->ZAtDCA()) < 1.)) {
	fHistSelectedTriggerPtPhi[0]->Fill(track->Pt(), track->Phi(), lPercentiles);
	fHistSelectedTriggerPtEta[0]->Fill(track->Pt(), track->Eta(), lPercentiles);    
      }
      if( (TMath::Abs(track->ZAtDCA()) < 0.5)) {
	fHistSelectedTriggerPtPhi[2]->Fill(track->Pt(), track->Phi(), lPercentiles);
	fHistSelectedTriggerPtEta[2]->Fill(track->Pt(), track->Eta(), lPercentiles);    
      }
      labelPrimOrSec=1;
    }
    else if(particle->IsSecondaryFromWeakDecay())      labelPrimOrSec=2;
    else if(particle->IsSecondaryFromMaterial())      labelPrimOrSec=3;
    else labelPrimOrSec=4;
    //    cout << "label is " << labelPrimOrSec<< endl;
    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	for(Int_t p=1; p<=4; p++){
	  if (labelPrimOrSec==p){
	    if( (TMath::Abs(track->ZAtDCA()) < 1))   fHistPrimaryTrigger[m][0]->Fill(p,particle->Pt() );
	    if( (TMath::Abs(track->ZAtDCA()) < 2))   fHistPrimaryTrigger[m][1]->Fill(p,particle->Pt() );
	    if( (TMath::Abs(track->ZAtDCA()) < 0.5))   fHistPrimaryTrigger[m][2]->Fill(p,particle->Pt() );
	  }
	}
      }
    }
    for(Int_t p=1; p<=4; p++){
      if (labelPrimOrSec==p){
	if((TMath::Abs(track->ZAtDCA()) < 1))   fHistPrimaryTrigger[5][0]->Fill(p,particle->Pt() );
	if( (TMath::Abs(track->ZAtDCA()) < 2))   fHistPrimaryTrigger[5][1]->Fill(p,particle->Pt() );
	if( (TMath::Abs(track->ZAtDCA()) < 0.5))   fHistPrimaryTrigger[5][2]->Fill(p,particle->Pt() );
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
 

  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];

  fOutputList = new TList();         
  fOutputList->SetOwner(kTRUE);     
  fOutputList2 = new TList();         
  fOutputList2->SetOwner(kTRUE);     
  fOutputList3 = new TList();         
  fOutputList3->SetOwner(kTRUE);     
  
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
  fSignalTree->Branch("fTreeVariablectau",               &fTreeVariablectau   , "fTreeVariablectau/D");
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
  fBkgTree->Branch("fTreeVariablectau",               &fTreeVariablectau   , "fTreeVariablectau/D");
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

  
  fHistPt = new TH1F("fHistPt", "p_{T} distribution of selected charged tracks in events used for AC", 300, 0, 30); 
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 
  fHistDCAxym1 = new TH1F("fHistDCAxym1", "DCAxy method 1 before DCA cuts", -10, 0, 10); 
  fHistDCAxym1->GetXaxis()->SetTitle("DCAxy method 1 (cm)");
  fHistDCAxym2 = new TH1F("fHistDCAxym2", "DCAxy method 2 before DCA cuts", -10, 0, 10); 
  fHistDCAxym2->GetXaxis()->SetTitle("DCAxy method 2 (cm)");
  fHistDCAzm1 = new TH1F("fHistDCAzm1", "DCAz method 1 before DCA cuts", -10, 0, 10); 
  fHistDCAzm1->GetXaxis()->SetTitle("DCAz method 1 (cm)");
  fHistDCAzm2 = new TH1F("fHistDCAzm2", "DCAz method 2 before DCA cuts", -10, 0, 10); 
  fHistDCAzm2->GetXaxis()->SetTitle("DCAz method 2 (cm)");
 
  fHistPtV0 = new TH1F("fHistPtV0", "p_{T} distribution of selected V0 in events used for AC", 300, 0, 30); 
  fHistPtV0->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMin = new TH1F("fHistPtTMin", "p_{T} distribution of reco trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMin->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMinMC = new TH1F("fHistPtTMinMC", "p_{T} distribution of true trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtvsMultBefAll= new TH2F("fHistPtvsMultBefAll", "p_{T} and centrality distribution of charged tracks in events w T>0", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtvsMult= new TH2F("fHistPtvsMult", "p_{T} and centrality distribution of charged tracks in events used for AC", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMult->GetYaxis()->SetTitle("Centrality");


  fHistZvertex= new TH1F("fHistZvertex", "Z vertex distribution of selected events used for AC", 40,-20,20);

  fHist_eta_phi= new TH2F("fHist_eta_phi", "Distribution of charged tracks in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_multiplicity=new TH1F("fHist_multiplicity", "fHist_multiplicity", 100, 0, 100); 
  fHist_multiplicity->SetTitle("Centrality distribution of events used for AC");

  fHistPDG=new TH1F("fHistPDG", "fHistPDG",3200, -3200, 3200);
  
  fHistEventMult=new TH1F("fHistEventMult", "fHistEventMult", 22, 0.5, 22.5);
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
  fHistEventMult->GetXaxis()->SetBinLabel(14,"NFirstPartMC < NFirstReco");
  fHistEventMult->GetXaxis()->SetBinLabel(15,"NFirstPartMC=0 && NFirstPart=1");
  fHistEventMult->GetXaxis()->SetBinLabel(16,"NFirstPartMC!=0 && NFirstPart!=0");
  fHistEventMult->GetXaxis()->SetBinLabel(17,"NSecondPartMC==0 && NSecondPart!=0 (NT>0)");
  fHistEventMult->GetXaxis()->SetBinLabel(18,"NSecondPartMC<NSecondRecoTrue (NT>0)");
  fHistEventMult->GetXaxis()->SetBinLabel(19,"NTrigger>0 && DoubleCounted");
  fHistEventMult->GetXaxis()->SetBinLabel(20,"Common daughterPos");
  fHistEventMult->GetXaxis()->SetBinLabel(21,"Common daughterNeg");
  fHistEventMult->GetXaxis()->SetBinLabel(22,"SelEv (ACEvents)");//tutti gli eventi usati per correlazione angolare

  fHistEventV0=new TH1F("fHistEventV0", "fHistEventV0",17, 0.5, 17.5);
  fHistEventV0->SetTitle("Number of V0 which progressively pass the listed selections");
  fHistEventV0->GetXaxis()->SetBinLabel(1,"All V0s");
  fHistEventV0->GetXaxis()->SetBinLabel(2,"V0s ok");
  fHistEventV0->GetXaxis()->SetBinLabel(3,"Filterbit daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(4,"Chis daughter tracks"); 
  fHistEventV0->GetXaxis()->SetBinLabel(5,"PID daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(6,"|eta daughters|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(7,"TPC track quality"); 
  fHistEventV0->GetXaxis()->SetBinLabel(8,"|eta_K0s|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(9,"more cuts + 0.45 < Mass < 0.55"); 
  fHistEventV0->GetXaxis()->SetBinLabel(10,"Lrejection"); 
  fHistEventV0->GetXaxis()->SetBinLabel(11,"0<pT<pTV0max (reco up to here)"); 
  fHistEventV0->GetXaxis()->SetBinLabel(12,"NV0(reco) in ev wNT>0"); 
  fHistEventV0->GetXaxis()->SetBinLabel(13,"NV0(reco true) in ev wNT>0");
  fHistEventV0->GetXaxis()->SetBinLabel(14,"NV0(MC) in ev wNT>0");
  fHistEventV0->GetXaxis()->SetBinLabel(15,"NV0(reco) in SelEv");
  fHistEventV0->GetXaxis()->SetBinLabel(16,"NV0(reco true) in SelEv");
  fHistEventV0->GetXaxis()->SetBinLabel(17,"NV0(MC) in SelEv"); 

  fHistTrack=new TH1F("fHistTrack", "fHistTrack", 15, 0.5, 15.5);
  fHistTrack->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistTrack->GetXaxis()->SetBinLabel(2,"Tracks after filterbit");
  fHistTrack->GetXaxis()->SetBinLabel(3,"Tracks with |eta| < 0.8"); 
  fHistTrack->GetXaxis()->SetBinLabel(4,"Track quality");
  fHistTrack->GetXaxis()->SetBinLabel(5,"TPCCrossedRows>70");
  fHistTrack->GetXaxis()->SetBinLabel(6,"Crossed rows/findable >0.8");
  fHistTrack->GetXaxis()->SetBinLabel(7,"Charged tracks");
  fHistTrack->GetXaxis()->SetBinLabel(8,"DCAxy < 0.010+0.035/pt**1.1");
  fHistTrack->GetXaxis()->SetBinLabel(9,"DCAz <2");
  fHistTrack->GetXaxis()->SetBinLabel(10,"N.trigger"); //NumberPrimary in all slected events 
  fHistTrack->GetXaxis()->SetBinLabel(11,"N.trigger MC");
  fHistTrack->GetXaxis()->SetBinLabel(12,"N.trigger (NT>0)"); //NumberPrimary in all slected events 
  fHistTrack->GetXaxis()->SetBinLabel(13,"N.trigger MC (NT>0)");
  fHistTrack->GetXaxis()->SetBinLabel(14,"N.trigger (NV0>0)"); //NumberPrimary in events with at least one V0 (one reco V0 for data, one true V0 for MC)
  fHistTrack->GetXaxis()->SetBinLabel(15,"N.trigger MC (NV0>0)");

  fMassV0= new TH1F("fMassV0", "Invariant mass of V0 candidates", 100, 0.45, 0.55);
  fMassV0->GetXaxis()->SetTitle("M_{#pi^+ #pi^-}");

  fHistSecondParticleAll= new TH2F("fHistSecondParticleAll", "Number of V0 MCTrue vs number V0 reco (T>0) ", 10,-0.5,9.5,10,-0.5,9.5);
  fHistSecondParticleAll->GetXaxis()->SetTitle("Number (reco)");
  fHistSecondParticleAll->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticleTruthAll= new TH2F("fHistSecondParticleTruthAll", "Number of V0 MCTrue vs number V0 reco (true) (T>0)", 10,-0.5,9.5,10,-0.5,9.5);
  fHistSecondParticleTruthAll->GetXaxis()->SetTitle("Number (reco true)");
  fHistSecondParticleTruthAll->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticle= new TH2F("fHistSecondParticle", "Number of V0 MCTrue vs number V0 reco (T>0, V>0)", 10,-0.5,9.5,10,-0.5,9.5);
  fHistSecondParticle->GetXaxis()->SetTitle("Number (reco)");
  fHistSecondParticle->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticleTruth= new TH2F("fHistSecondParticleTruth", "Number of V0 MCTrue vs number V0 reco (true) (T>0, V>0)", 10,-0.5,9.5,10,-0.5,9.5);
  fHistSecondParticleTruth->GetXaxis()->SetTitle("Number (reco true)");
  fHistSecondParticleTruth->GetYaxis()->SetTitle("Number (MC)");

  fHistMassvsPt = new TH2F *[6];
  fHistMassvsPt_tagli = new TH2F *[6];

  for(Int_t j=0; j<6; j++){
    fHistMassvsPt[j] = new TH2F(Form("fHistMassvsPt_" +fV0+ "_%i",j),Form("fHistMassvsPt_"+fV0 + "_%i"+" (cuts on PID, eta and track quality of daughters + eta V0) ",j),400,0.3,0.7,160,0,16); 
    fHistMassvsPt_tagli[j] = new TH2F(Form("fHistMassvsPt_" +fV0+ "_tagli_%i",j),Form("fHistMassvsPt_" +fV0+ "_tagli_%i" + " (all selections on V0 applied)",j),400,0.3,0.7,160,0,16);
    fHistMassvsPt[j]->GetXaxis()->SetTitle("Invariant mass of V0 candidate");
    fHistMassvsPt[j]->GetYaxis()->SetTitle("p_{T} of V0 candidate");   
    fHistMassvsPt_tagli[j]->GetXaxis()->SetTitle("Invariant mass of V0 candidate");
    fHistMassvsPt_tagli[j]->GetYaxis()->SetTitle("p_{T} of V0 candidate");   
  }
  
  fHistMultvsTrigger=new TH2F("fHistMultvsTrigger", "Centrality of selected events (T>0, V>0) vs number of trigger particles", 20, -0.5, 19.5, 100, 0, 100);
  fHistMultvsTrigger->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTrigger->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerMCTruth=new TH2F("fHistMultvsTriggerMCTruth", "Centrality of selected events (T>0, V>0) vs number of trigger particles, MC Truth", 20, -0.5, 19.5, 100, 0, 100);
  fHistMultvsTriggerMCTruth->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruth->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerAll=new TH2F("fHistMultvsTriggerAll", "Centrality of events w T>0 vs number of trigger particles", 20, -0.5, 19.5, 100, 0, 100);
  fHistMultvsTriggerAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerAll->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsTriggerMCTruthAll=new TH2F("fHistMultvsTriggerMCTruthAll", "Centrality of events w T>0 vs number of trigger particles, MC Truth", 20, -0.5, 19.5, 100, 0, 100);
  fHistMultvsTriggerMCTruthAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruthAll->GetYaxis()->SetTitle("Centrality");

  fHistMultvsTriggerBefAll=new TH2F("fHistMultvsTriggerBefAll", "Centrality of events vs number of trigger particles", 20, -0.5, 19.5, 100, 0, 100);
  fHistMultvsTriggerBefAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerBefAll->GetYaxis()->SetTitle("Centrality");
  
 
  fHistMultvsTriggerMCTruthBefAll=new TH2F("fHistMultvsTriggerMCTruthBefAll", "Centrality of events vs number of trigger particles, MC Truth", 20, -0.5, 19.5, 100, 0, 100);
  fHistMultvsTriggerMCTruthBefAll->GetXaxis()->SetTitle("Number of Trigger particles");
  fHistMultvsTriggerMCTruthBefAll->GetYaxis()->SetTitle("Centrality");
  

  fHistMultvsV0=new TH2F("fHistMultvsV0", "Centrality of selected events (T>0, V0>0) vs number of reco V0s",20, -0.5, 19.5,100, 0, 100 );
  fHistMultvsV0->GetXaxis()->SetTitle("Number of reco V0 particles");
  fHistMultvsV0->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0Truth=new TH2F("fHistMultvsV0Truth", "Centrality of selected events (T>0, V0>0) vs number of reco true V0s",20, -0.5, 19.5,100, 0, 100 );
  fHistMultvsV0Truth->GetXaxis()->SetTitle("Number of reco true V0 particles");
  fHistMultvsV0Truth->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0MC=new TH2F("fHistMultvsV0MC", "Centrality of selected events (T>0, V0>0) vs number of true V0s",20, -0.5, 19.5,100, 0, 100 );
  fHistMultvsV0MC->GetXaxis()->SetTitle("Number of V0 true particles");
  fHistMultvsV0MC->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0All=new TH2F("fHistMultvsV0All", "Centrality of events w T>0 vs number of reco V0s",20, -0.5, 19.5,100, 0, 100 );
  fHistMultvsV0All->GetXaxis()->SetTitle("Number of V0 particles");
  fHistMultvsV0All->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0AllTruth=new TH2F("fHistMultvsV0AllTruth", "Centrality of events w T>0 vs number of reco true V0s",20, -0.5, 19.5,100, 0, 100 );
  fHistMultvsV0AllTruth->GetXaxis()->SetTitle("Number of V0 reco particles");
  fHistMultvsV0AllTruth->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0MCAll=new TH2F("fHistMultvsV0MCAll", "Centrality of events w T>0 vs number of true V0s",20, -0.5, 19.5,100, 0, 100 );
  fHistMultvsV0MCAll->GetXaxis()->SetTitle("Number of V0 true particles");
  fHistMultvsV0MCAll->GetYaxis()->SetTitle("Centrality");

  fHistMassPhoton=new TH1F("fHistMassPhoton", "Inv Mass of two V0 daughters as electrons, after V0 mass selection",  300, 0, 1.5);
  fHistMass2Photon=new TH1F("fHistMass2Photons", "Inv MassSquared of two daughters as electrons",  100, -3, 3);
  
  fHistPtArmvsAlpha=new TH2F("fHistPtArmvsAlpha", "Distribution of V0 candidates before cuts on Arm/Lrej/mass",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlpha->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlpha->GetYaxis()->SetTitle("Pt Armenteros");

  fHistPtArmvsAlphaAfterSelection=new TH2F("fHistPtArmvsAlphaAfterSelection", "Distribution of V0 candidates after applying mass cuts", 80, -1, 1,80, 0, 0.3);

  fHistPtArmvsAlphaAfterPhotonSelection=new TH2F("fHistPtArmvsAlphaAfterPhotonSelection", "Distribution of V0 candidates after photon cuts",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlphaAfterPhotonSelection->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlphaAfterPhotonSelection->GetYaxis()->SetTitle("Pt Armenteros");

  fHistPtArmvsAlphaAfterLambdaRejectionSelection=new TH2F("fHistPtArmvsAlphaAfterLambdaRejectionSelection", "Distribution of V0 candidates after Lrejection",  80, -1, 1,80, 0, 0.3);
  fHistPtArmvsAlphaAfterLambdaRejectionSelection->GetXaxis()->SetTitle("Alpha");
  fHistPtArmvsAlphaAfterLambdaRejectionSelection->GetYaxis()->SetTitle("Pt Armenteros");

  fHistTrigger=new TH1F("fHistTrigger", "Number of reco trigger particle distribution for selected events (also T=0)", 30, -0.5, 29.5); // each entry is an event

  fHistTriggerwV0=new TH1F("fHistTriggerwV0", "Number of reco trigger particle distribution for events used for AC", 30, -0.5, 29.5); // each entry is an event

  fHistTriggerMCTruth=new TH1F("fHistTriggerMCTruth", "Number of true trigger particle distribution for selected events (also T=0)", 30, -0.5, 29.5); // each entry is an event

  fHistTriggerwV0MCTruth=new TH1F("fHistTriggerwV0MCTruth", "Number of true trigger particle distribution for events used for AC", 30, -0.5, 29.5); // each entry is an event

  fHistMultiplicityVsVertexZ=new TH2F("fHistMultiplicityVsVertexZ", "Centrality vs Z vertex of selected events with NT>0 and NV0>0 ",  20, -10, 10,100, 0, 100);
      
  fHistTriggervsMult=new TH1F("fHistTriggervsMult", "Numero di particelle di trigger nei vari intervalli di centralita'", 100, 0, 100);
  fHistTriggervsMult->GetXaxis()->SetTitle("Centrality");

  fHistTriggervsMultMC=new TH1F("fHistTriggervsMultMC", "Numero di particelle di trigger (MCtruth) nei vari intervalli di centralita'", 100, 0, 100);
  fHistTriggervsMultMC->GetXaxis()->SetTitle("Centrality");

  fHistGeneratedTriggerPtPhi=new TH3F("fHistGeneratedTriggerPtPhi", "p_{T} and #phi distribution of generated trigger particles (charged, primary)", 300, 0, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
  fHistGeneratedTriggerPtPhi->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtPhi->GetYaxis()->SetTitle("#phi");

  fHistGeneratedTriggerPtEta=new TH3F("fHistGeneratedTriggerPtEta", "p_{T} and #eta distribution of generated trigger particles (primary, charged)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistGeneratedTriggerPtEta->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedTriggerPtEta->GetYaxis()->SetTitle("#eta");

  fHistSelectedTriggerPtPhi= new TH3F*[3];
  for(Int_t j=0; j<3; j++){
    fHistSelectedTriggerPtPhi[j]=new TH3F(Form("fHistSelectedTriggerPtPhi_%i",j), "p_{T} and #phi distribution of selected trigger particles (primary)", 300, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
    fHistSelectedTriggerPtPhi[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedTriggerPtPhi[j]->GetYaxis()->SetTitle("#phi");
  }

  fHistSelectedTriggerPtEta= new TH3F*[3];
  for(Int_t j=0; j<3; j++){
    fHistSelectedTriggerPtEta[j]=new TH3F(Form("fHistSelectedTriggerPtEta_%i",j), "p_{T} and #eta distribution of selected trigger particles (primary)", 300, 0, 30, 400,-1.2, 1.2,  100, 0, 100);
    fHistSelectedTriggerPtEta[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedTriggerPtEta[j]->GetYaxis()->SetTitle("#eta");
  }
  
    fHistGeneratedV0PtPhi=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
  fHistGeneratedV0PtPhi[j]=new TH3F(Form("fHistGeneratedV0PtPhi_%i",j), "p_{T} and #phi distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,0, 2*TMath::Pi(),  100, 0, 100 );
  fHistGeneratedV0PtPhi[j]->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedV0PtPhi[j]->GetYaxis()->SetTitle("#phi");
  }
  
  fHistSelectedV0PtPhi=new TH3F*[7];
  for(Int_t j=0; j<7; j++){
    fHistSelectedV0PtPhi[j]=new TH3F(Form("fHistSelectedV0PtPhi_%i",j), "p_{T} and #phi distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,0, 2*TMath::Pi() ,  100, 0, 100);
    fHistSelectedV0PtPhi[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedV0PtPhi[j]->GetYaxis()->SetTitle("#phi");
  }
  
  fHistGeneratedV0PtEta=new TH3F*[2];
 for(Int_t j=0; j<2; j++){
   fHistGeneratedV0PtEta[j]=new TH3F(Form("fHistGeneratedV0PtEta_%i",j), "p_{T} and #eta distribution of generated V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
  fHistGeneratedV0PtEta[j]->GetXaxis()->SetTitle("p_{T}");
  fHistGeneratedV0PtEta[j]->GetYaxis()->SetTitle("#eta");
 }
  
  fHistSelectedV0PtEta=new TH3F*[7];
  for(Int_t j=0; j<7; j++){
    fHistSelectedV0PtEta[j]=new TH3F(Form("fHistSelectedV0PtEta_%i",j), "p_{T} and #eta distribution of selected V0 particles (K0s, primary, events w T>0)", 300, 0, 30, 400,-1.2,1.2,  100, 0, 100 );
    fHistSelectedV0PtEta[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedV0PtEta[j]->GetYaxis()->SetTitle("#eta");
  }

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


  fHistPrimaryTrigger= new TH2F **[6];
  for(Int_t j=0; j<6; j++){
    fHistPrimaryTrigger[j]=new TH2F*[3];
    for(Int_t i=0; i<3; i++){
      fHistPrimaryTrigger[j][i]=new TH2F(Form("fHistPrimaryTrigger_%i_cut%i", j,i), "Trigger MC (selected)", 4, 0.5, 4.5, 100, 0,30 );
      fHistPrimaryTrigger[j][i]->GetXaxis()->SetBinLabel(1,"Primary selected triggers");
      fHistPrimaryTrigger[j][i]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected triggers"); 
      fHistPrimaryTrigger[j][i]->GetXaxis()->SetBinLabel(3,"Secondary from material selected triggers"); 
      fHistPrimaryTrigger[j][i]->GetYaxis()->SetTitle("p_{T}"); 
    }
  }


  fHistPrimaryV0= new TH2F**[6];
  for(Int_t j=0; j<6; j++){
    fHistPrimaryV0[j]=new TH2F*[7];
    for(Int_t i=0; i<7; i++){
      fHistPrimaryV0[j][i]=new TH2F(Form("fHistPrimaryV0_%i_cut%i",j,i), "V0 MC (K0s, selected)", 4, 0.5, 4.5, 160, 0, 16);
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s"); 
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s"); 
      fHistPrimaryV0[j][i]->GetYaxis()->SetTitle("p_{T}");
    } 
  }

  //  TString molteplicit[6]={"0-7","7-15","15-25","25-40","40-70",">70"};
  fHistMultiplicityOfMixedEvent=new TH1F("fHistMultiplicityOfMixedEvent", "Distribution of number of events used for the mixing", 20, 0.5, 20.5);

  fEventCuts.AddQAplotsToList(fOutputList);
  
  fOutputList->Add(fHistPDG);
  fOutputList->Add(fHistEventMult);
  fOutputList->Add(fHistEventV0);
  fOutputList->Add(fHistTrack); 
  
  //istogrammi riempiti ad ogni evento selezionato (ossia utilizzato per AC) con una entry
  fOutputList->Add(fHistZvertex);
  fOutputList->Add(fHist_multiplicity); 
  fOutputList->Add(fHistMultiplicityVsVertexZ);
  fOutputList->Add(fHistMultvsTriggerBefAll);
  fOutputList->Add(fHistMultvsTriggerMCTruthBefAll);
  fOutputList->Add(fHistMultvsTriggerAll);
  fOutputList->Add(fHistMultvsTriggerMCTruthAll);
  fOutputList->Add(fHistMultvsTrigger);
  fOutputList->Add(fHistMultvsTriggerMCTruth);
  fOutputList->Add(fHistMultvsV0All);
  fOutputList->Add(fHistMultvsV0AllTruth);
  fOutputList->Add(fHistMultvsV0MCAll);
  fOutputList->Add(fHistMultvsV0);
  fOutputList->Add(fHistMultvsV0Truth);
  fOutputList->Add(fHistMultvsV0MC);
  fOutputList->Add(fHistSecondParticleAll);
  fOutputList->Add(fHistSecondParticleTruthAll);
  fOutputList->Add(fHistSecondParticle);  
  fOutputList->Add(fHistSecondParticleTruth);
  fOutputList->Add(fHistPt);     
  fOutputList->Add(fHistDCAxym1);     
  fOutputList->Add(fHistDCAxym2);    
  fOutputList->Add(fHistDCAzm1);      
  fOutputList->Add(fHistDCAzm2);     
  fOutputList->Add(fHistPtV0);     
  fOutputList->Add(fHistPtTMin);     
  fOutputList->Add(fHistPtTMinMC);     
  fOutputList->Add(fHistPtvsMult);       
  fOutputList->Add(fHistPtvsMultBefAll);       
  fOutputList->Add(fHist_eta_phi);
  fOutputList->Add(fHistTriggervsMult);
  fOutputList->Add(fHistTriggervsMultMC);
  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistTriggerwV0);
  fOutputList->Add(fHistTriggerMCTruth);
  fOutputList->Add(fHistTriggerwV0MCTruth);

  fOutputList->Add(fMassV0);

  for(Int_t j=0; j < 6; j++){
    fOutputList->Add(fHistMassvsPt[j]);
    fOutputList->Add(fHistMassvsPt_tagli[j]);
    
    for(Int_t i=0; i<3; i++){  
      fOutputList3->Add(fHistPrimaryTrigger[j][i]); 
    }
   
    for(Int_t l=0; l<7; l++){
      fOutputList2->Add(fHistPrimaryV0[j][l]);      
    }
    
  }
  
  
  fOutputList->Add(fHistMass2Photon);
  fOutputList->Add(fHistMassPhoton);
  fOutputList->Add(fHistPtArmvsAlpha);
  fOutputList->Add(fHistPtArmvsAlphaAfterSelection);  
  fOutputList->Add(fHistPtArmvsAlphaAfterPhotonSelection);  
  fOutputList->Add(fHistPtArmvsAlphaAfterLambdaRejectionSelection);  
  fOutputList->Add(fHistMultiplicityOfMixedEvent);
  
  
  fOutputList3->Add(fHistGeneratedTriggerPtPhi);
  fOutputList3->Add(fHistGeneratedTriggerPtEta);
  
  for(Int_t j=0; j < 3; j++){
    fOutputList3->Add(fHistSelectedTriggerPtPhi[j]);
     fOutputList3->Add(fHistSelectedTriggerPtEta[j]);
  }
 
  for(Int_t j=0; j < 2; j++){
   fOutputList2->Add(fHistGeneratedV0PtPhi[j]); 
   fOutputList2->Add(fHistGeneratedV0PtEta[j]); 
  }
  for(Int_t j=0; j < 7; j++){
   fOutputList2->Add(fHistSelectedV0PtPhi[j]);
   fOutputList2->Add(fHistSelectedV0PtEta[j]);
  }
 
  fOutputList->Add(fHistReconstructedV0PtMass);
  fOutputList->Add(fHistSelectedV0PtMass);
  fOutputList->Add(fHistResolutionTriggerPt);
  fOutputList->Add(fHistResolutionTriggerPhi);
  fOutputList->Add(fHistResolutionTriggerEta);
  fOutputList->Add(fHistResolutionV0Pt);
  fOutputList->Add(fHistResolutionV0Phi);
  fOutputList->Add(fHistResolutionV0Eta);
 
  PostData(1, fOutputList);  
  PostData(2, fSignalTree);       
  PostData(3, fBkgTree); 
  PostData(4, fOutputList2);     
  PostData(5, fOutputList3);     
     
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
  Float_t LastzBin;
  Float_t LastcentralityBin;

  fHistEventMult->Fill(1);

  // AliMCEvent   *lMCevent  = 0x0;
  // AliStack     *lMCstack  = 0x0;
  // TClonesArray *arrayMC = 0x0;
   
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());  
  if(!fAOD) {
    AliWarning("Error: AOD event not available \n");
    PostData(1, fOutputList);       

    PostData(2, fSignalTree);       
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    return;
  }        
  
  
 /// Use the event cut class to apply the required selections
 if (!fEventCuts.AcceptEvent(fAOD)) {   
 PostData(1, fOutputList);
 PostData(2, fSignalTree );
 PostData(3,fBkgTree);
  PostData(4, fOutputList2);     
  PostData(5, fOutputList3);      return;
 }
  
  
  Int_t iTracks(fAOD->GetNumberOfTracks());         
  Int_t V0Tracks(fAOD->GetNumberOfV0s());           
  Evcounter++;  
  //c cout << "\n \n \n ********************************************************* "<< endl;
  //c cout << "numero dell'evento "<<  Evcounter << endl; 
  //c cout << "number of tracks before any cut " << iTracks << endl;
  //c cout << "number of V0 before any cut " << V0Tracks << endl;
  
  //VERTEX SELECTION AND TRIGGER
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  const AliAODVertex *lPrimaryBestAODVtx = fAOD->GetPrimaryVertex();
  if (!lPrimaryBestAODVtx){
    AliWarning("No prim. vertex in AOD... return!");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     return;
  }
  fHistEventMult->Fill(2);
  AliVVertex *vertexmain =0x0;
  vertexmain = (AliVVertex*) lPrimaryBestAODVtx;
  lPrimaryBestAODVtx->GetXYZ(lBestPrimaryVtxPos);


  if (TMath::Abs(lBestPrimaryVtxPos[2])>10.){
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    //c cout << "z vertex selection failed " << endl;
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
  PostData(5, fOutputList3);     
  return;
  }
  fHistEventMult->Fill(3);


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
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
  PostData(5, fOutputList3);     
 return;
  }
  fHistEventMult->Fill(4);


  //CENTRALITY SELECTIONS for Pb-Pb //preso da AliAnalysisTaskKPFemto
  //lcentrality=centrality->GetCentralityPercentile("V0A");
  
  Float_t lPercentiles = 0;
  
 
  //This will work for both ESDs and AODs
  AliMultSelection *MultSelection = (AliMultSelection*) fAOD -> FindListObject("MultSelection");

  if ( MultSelection ){
  //c cout << "mult sel ok" << endl;
  lPercentiles= MultSelection->GetMultiplicityPercentile("V0M");
  }else{
  AliInfo("Didn't find MultSelection!"); 
  }
 
  //c cout << "centralita' dell evento " << lPercentiles << endl;
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
    //c cout << "lPercentiles >= 200" << endl;
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
    PostData(5, fOutputList3);      return;   
  }

  fHistEventMult->Fill(5);

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fAOD->IsPileupFromSPD();
  if(isPileUpSpd){ 
    PostData(1,fOutputList );
    PostData(2, fSignalTree );
    //c cout << "return: event is pile up " << endl;
    PostData(3,fBkgTree); 
    PostData(4, fOutputList2); 
  PostData(5, fOutputList3);     
 return;
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
    //c cout << "event does not fulfil centrality selection criteria " << endl;     
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     
    return;
  }
  
  //c cout << "event has passed selection criteria.... first and second particles to be analyzed ...."<< endl;


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
  
  
  //  //c cout << "dopo centrality selection " << endl;
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
  if(lPercentiles < 5.0) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
  else if(lPercentiles < 10.0) centralityBin=18;
  else if(lPercentiles < 30.0) centralityBin=17;
  else if(lPercentiles < 50.) centralityBin=16;
  else if(lPercentiles <= 100.) centralityBin=15;

  ////c cout << "6 " << endl;
  //c cout << "  " << zBin <<" centrbin " << centralityBin<< endl;
  if (((centralityBin+1) >fnMultBins) || ((zBin+1) > fzVertexBins)){ 
    //c cout<<" ##################  WARNING: I'm going to break bacause of dimensional issues ########################"<<endl;
  }
  // //se decommento, ho seg viol
  
  //  if(Evcounter==1)  fEventColl[zBin][centralityBin]->FifoShift();
//unuseful  if (Evcounter==1){
//unuseful    LastzBin=zBin;
//unuseful    LastcentralityBin=centralityBin;
//unuseful  }
//unuseful  if (FifoShiftok==kTRUE && LastzBin==zBin && LastcentralityBin==centralityBin)  fEventColl[zBin][centralityBin]->FifoShift();
//unuseful  else  fEventColl[zBin][centralityBin]->FifoClear();
//unuseful  //working  fEventColl[zBin][centralityBin]->FifoShift();
  //cout << "6 " << endl;
  //unuseful FifoShiftok=kFALSE;
  // // //se decommento, ho seg viol 

  fEventColl[zBin][centralityBin]->FifoClear();
  fEvt = fEventColl[zBin][centralityBin]->fEvt;
  // // cout << "7 " << endl;
  //-----------------------------------LOOP OVER THE TRACKS

  Float_t nTPCCrossedRows=0.;
  float rationCrnFind=0;
  Float_t nsigmaTOFj=3;
  Float_t nsigmaTPCj=3;
  Int_t charge=0;
  Int_t NumberFirstParticleAll=0;
  Int_t NumberFirstParticle=0;
  Int_t NumberFirstParticleMC=0;
  Int_t NumberFirstParticle_finale=0;
  Int_t NumberSecondParticleRecoTrue=0;
  Int_t NumberSecondParticle=0;
  Int_t NumberSecondParticleMC=0;
  Int_t NumberSecondParticleAll=0;
  Double_t Ptintermediate=0;
  Double_t selectedtrackID=0;
  Int_t pos0or1=0;
  Int_t neg0or1=0;
  Int_t CharegFirstParticle=0;
  Double_t dz[2] = {-999.,-999.};
  Double_t d[2] = {-999.,-999.};

  //LOOP FOR FIRST PARTICLE
 
  Float_t ptTriggerMinimoDati=10000;
  
  for(Int_t i=0; i < iTracks; i++) {
    track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));        
    fHistTrack->Fill(1);
    if(!track) continue;
    //if(track->GetID()<0)continue;
    //  if(!track->IsOn(AliAODTrack::kTPCrefit)) continue;
    if(!track->TestFilterBit(128)) continue;
    fHistTrack->Fill(2);
    
    if(TMath::Abs(track->Eta())>0.8)  continue;
    fHistTrack->Fill(3);
    
    if(track->Chi2perNDF()>4.)continue;
    fHistTrack->Fill(4);
    
    nTPCCrossedRows=track->GetTPCNCrossedRows();
    if(nTPCCrossedRows<70) continue;
    fHistTrack->Fill(5);
        
    rationCrnFind=nTPCCrossedRows/track->GetTPCNclsF();
    if(rationCrnFind<0.8)  continue;
    fHistTrack->Fill(6);
    
    if((track->Charge())==0) continue;
    fHistTrack->Fill(7);

    dz[0] = track->DCA();    // the TPC one should be applied the other biases the CF --> from Maciejs note --> FIXME to be checked 
    dz[1] = track->ZAtDCA(); // for those two lines check AliAODTrack.h // FIXME these two lines produce shifted distributions, known problem, asked Mac and Marian. 

    fHistDCAxym1->Fill(track->DCA());
    fHistDCAzm1->Fill(track->ZAtDCA());

    Double_t  covd[3];
    //    if (track->DCA()==-999. || track->ZAtDCA()==-999.){
    AliAODTrack* track_clone=(AliAODTrack*)track->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
    Bool_t isDCA = track_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),9999.,d,covd);
    delete track_clone;

    fHistDCAxym2->Fill(d[0]);
    fHistDCAzm2 ->Fill(d[1]);


    //   if(TMath::Abs(track->DCA())> (0.0105 + 0.0350/pow(track->Pt(),1.1))) continue;
    if(TMath::Abs(d[0])> (0.0105 + 0.0350/pow(track->Pt(),1.1))) continue;
    fHistTrack->Fill(8);
    if(TMath::Abs(d[1])> 2.) continue;
    //    if(TMath::Abs(track->ZAtDCA())> 2.) continue;
    fHistTrack->Fill(9);

    //to determine efficiency of trigger particle
    label=track->GetLabel();
    Generated = kFALSE;     
    if(fReadMCTruth){
      if(fMCEvent){
	ProcessMCParticles(Generated, track, labelPrimOrSec, lPercentiles, isV0);
      }
    }

    if(track->Pt()> fminPtj && track->Pt()<fmaxPtj){
      NumberFirstParticle++;
      if(track->Pt()< ptTriggerMinimoDati) ptTriggerMinimoDati=track->Pt(); 
      if (track->Pt()>Ptintermediate){
	Ptintermediate=track->Pt();
	selectedtrackID= track->GetID();
	NumberFirstParticle_finale= NumberFirstParticle;     
      }
      
      if((!fReadMCTruth)|| (fReadMCTruth && isEfficiency)){
	//save first particle information (leading particle)
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fCharge       = track->Charge();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fPt           = track->Pt();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fEta          = track->Eta();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fPhi          = track->Phi();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fTheta        = track->Theta();
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fDCAz         = d[1];
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fDCAxy        = d[0];
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fMultiplicity = lPercentiles;
	fEvt->fReconstructedFirst[NumberFirstParticle-1].fZvertex      = lBestPrimaryVtxPos[2];
	fEvt->fReconstructedFirst[NumberFirstParticle-1].isP           = labelPrimOrSec;
      }
      fHistPtvsMultBefAll->Fill(track->Pt(), lPercentiles);
    }
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
	if(trParticle->Pt()<= fminPtj || trParticle->Pt()>=fmaxPtj)continue;
	if(trParticle->Pt()< ptTriggerMinimoMC) ptTriggerMinimoMC=trParticle->Pt();
	NumberFirstParticleMC++;
	if (isEfficiency) continue;	  
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fCharge       = trParticle->Charge();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPt           = trParticle->Pt();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fEta          = trParticle->Eta();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fPhi          = trParticle->Phi();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fTheta        = trParticle->Theta();
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fDCAz         = 0;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fDCAxy        = 0;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fMultiplicity = lPercentiles;
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].fZvertex      = lBestPrimaryVtxPos[2];
	fEvt->fReconstructedFirst[NumberFirstParticleMC-1].isP           = 1;

      }
    }
  }
  
  fHistPtTMin->Fill(ptTriggerMinimoDati);
  fHistPtTMinMC->Fill(ptTriggerMinimoMC);
  fHistTrigger->Fill(NumberFirstParticle);
  fHistTriggerMCTruth->Fill(NumberFirstParticleMC);
  fHistMultvsTriggerBefAll->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruthBefAll->Fill(NumberFirstParticleMC, lPercentiles);

  fHistTrack->AddBinContent(10, NumberFirstParticle);
  fHistTrack->AddBinContent(11, NumberFirstParticleMC);
  if (NumberFirstParticle>0) fHistEventMult->Fill(10);   
  if (NumberFirstParticleMC>0) fHistEventMult->Fill(11);   
  if(NumberFirstParticle>1)fHistEventMult->Fill(12); 
  if(NumberFirstParticleMC>1)fHistEventMult->Fill(13); 
  if(NumberFirstParticleMC<NumberFirstParticle) fHistEventMult->Fill(14); 
  if(NumberFirstParticleMC==0 && NumberFirstParticle==1) fHistEventMult->Fill(15); 
  if(NumberFirstParticleMC!=0 && NumberFirstParticle!=0) fHistEventMult->Fill(16); 

  
  if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)) NumberFirstParticleAll=NumberFirstParticle;
  else if (fReadMCTruth && !isEfficiency) NumberFirstParticleAll=NumberFirstParticleMC;
  if(NumberFirstParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    //c cout  << "event does not have Trigger particles " << endl;     
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     
    return;
  }
 
  fHistTrack->AddBinContent(12, NumberFirstParticle);
  fHistTrack->AddBinContent(13, NumberFirstParticleMC);
 

  //c cout <<"candidati first particle nell'evento analizzato" << "   " <<  NumberFirstParticle << endl;        
  //c cout <<"candidati first particle nell'evento analizzato, MC truth" << "   " <<  NumberFirstParticleMC << endl;        
  //c cout << "\n \n \n Pt della particella con Pt maggiore " << Ptintermediate << endl;
   
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
  Float_t kMinCosAngle[2]={0.995, 0.997};
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
  TLorentzVector vPos;
  TLorentzVector vNeg;
  TLorentzVector vPhoton;
  Double_t MassPhoton;
  Double_t MassPhoton2;
  Double_t MassElectron= 0.0005109989461;
  Double_t massLambda = 1.115683;
  Int_t isaK0s=0;
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
  
  //c cout << "\n \n  here I start the loop on v0s " << endl;
  for(Int_t i(0); i < V0Tracks; i++) {       
    //c cout << "i "<< i << endl;
    isaK0s=0; //it will be put to 1 for true K0s in MC
    fHistEventV0->Fill(1);
    rapidityV0[2]={0};    
    EV0[2]={0};
    kctau[2]={0};
    goodPiPlus=kFALSE;
    goodPiMinus=kFALSE;
    AliAODv0* v0 = fAOD->GetV0(i);
    if(!v0) continue;       
    if(v0->GetOnFlyStatus()) continue;
    fHistEventV0->Fill(2);
    AliAODTrack* tempTrack=(AliAODTrack*)v0->GetDaughter(0);
    if(tempTrack->Charge()>0) {pos0or1=0; neg0or1=1;}
    else {pos0or1=1; neg0or1=0;}
    AliAODTrack* prongTrackPos=(AliAODTrack*)v0->GetDaughter(pos0or1);
    AliAODTrack* prongTrackNeg=(AliAODTrack*)v0->GetDaughter(neg0or1);

    if (!prongTrackPos || !prongTrackNeg) {
      Printf("ERROR: Could not retreive one of the daughter track");
      continue;
    }

    // Filter like-sign V0 (next: add counter and distribution)
    if ( prongTrackPos->GetSign() == prongTrackNeg->GetSign()) {
      continue;
    }

    labelPos=prongTrackPos->GetLabel();
    labelNeg=prongTrackNeg->GetLabel();
    //c cout <<	"label tracce figlie (pos e neg) "<< labelPos<< endl;
    //c cout <<	"label tracce figlie (pos e neg) "<< labelNeg<< endl;

    // TClonesArray* AODMCTrackArray =0x0;  
    if (fReadMCTruth){
      //      cout << "\n sono in MC event: I analyze all v0 reconstructed by ALICE" << endl;
      if (fMCEvent){
	AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}
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
	//c cout << "label tracce madri (pos e neg) " << 	labelMotherPos<< endl;
	//c cout << "label tracce madri (pos e neg) " << 	labelMotherNeg<< endl;

	MotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelMotherPos));
	MotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(labelMotherNeg));
	
	if (labelMotherPos>=0) PdgMotherPos = MotherPos->GetPdgCode();
	if (labelMotherNeg>=0)	PdgMotherNeg = MotherNeg->GetPdgCode();
	
      }
    }

    if(ParticleType==0) EV0[ParticleType]= TMath::Sqrt(pow(v0->MassK0Short(),2)+ pow( v0->Px(),2)+ pow( v0->Py(),2)+ pow( v0->Pz(),2));
    if(ParticleType==1) EV0[ParticleType]= TMath::Sqrt(pow(v0->MassLambda(),2)+ pow( v0->Px(),2)+ pow( v0->Py(),2)+ pow( v0->Pz(),2));

    rapidityV0[ParticleType]  = 0.5*TMath::Log(( EV0[ParticleType] + v0->Pz()) / (EV0[ParticleType]   - v0->Pz()) );
    
    fTreeVariablePtV0=TMath::Sqrt( pow( v0->Px(),2)+ pow( v0->Py(),2));
    fTreeVariableInvMassK0s=v0->MassK0Short();    
    fTreeVariableInvMassLambda=v0->MassLambda();    
    fTreeVariableInvMassAntiLambda=v0->MassAntiLambda();    
    kctau[ParticleType]=Mass[ParticleType]*v0->DecayLengthV0(lBestPrimaryVtxPos)/TMath::Sqrt( pow( v0->Px(),2)+ pow( v0->Py(),2) + pow(v0->Pz(),2));
    vPos.SetPxPyPzE(prongTrackPos->Px(),prongTrackPos->Py(),prongTrackPos->Pz(),sqrt(pow(prongTrackPos->P(),2) + pow(MassElectron,2)));
    vNeg.SetPxPyPzE(prongTrackNeg->Px(),prongTrackNeg->Py(),prongTrackNeg->Pz(),sqrt(pow(prongTrackNeg->P(),2) + pow(MassElectron,2)));
    vPhoton = vPos+vNeg;
    MassPhoton2=vPhoton.M2();
    MassPhoton=vPhoton.M();
    fHistMass2Photon->Fill(MassPhoton2);
   

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
    
    //V0 cuts>
    // if(TMath::Abs(rapidityV0[ParticleType])>ycut[ParticleType])             continue;
    //fHistEventV0->Fill(8);    
    if(TMath::Abs(v0->Eta()) > kEtacut[ParticleType])	   	continue;
    fHistEventV0->Fill(8);    
    
    //    if(v0->PtProng(pos0or1) < .15) continue;
    //if(v0->PtProng(neg0or1) < .15) continue;
   
    //    cout << lPercentiles << endl;
    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	//	//c cout<< "ho riempito isto masse no tagli " << endl;
	fHistMassvsPt[m]->Fill(v0->MassK0Short(),v0->Pt());
      }
    }
    fHistMassvsPt[5]->Fill(v0->MassK0Short(),v0->Pt());

    if(fReadMCTruth){
      if (fMCEvent){
	cout << "\n this particle has passe dall but pt cuts: let's fill the mass Pt histo for true reco K0s "<< endl;
	if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg  && MotherPos->IsPhysicalPrimary()){
	  fHistReconstructedV0PtMass->Fill(v0->MassK0Short(),v0->Pt(), lPercentiles);
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
    if(kctau[ParticleType]>kctauval[ParticleType])                   continue;
    //}
   
    //if(v0->PtArmV0()< 0.2*TMath::Abs(v0->AlphaV0()))                    continue;
    //fHistPtArmvsAlphaAfterSelection->Fill(v0->AlphaV0(), v0->PtArmV0());    
    fHistPtArmvsAlpha->Fill(v0->AlphaV0(),v0->PtArmV0());       

    if(v0->MassK0Short()< 0.55 && v0->MassK0Short()> 0.45) {
      fHistPtArmvsAlphaAfterSelection->Fill(v0->AlphaV0(),v0->PtArmV0());    
    }

    fHistMassPhoton->Fill(MassPhoton);   
    if(MassPhoton>0.1){
      fHistPtArmvsAlphaAfterPhotonSelection->Fill(v0->AlphaV0(),v0->PtArmV0());    
    }

    if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
      fHistPtArmvsAlphaAfterLambdaRejectionSelection->Fill(v0->AlphaV0(),v0->PtArmV0());    
    }

    if(v0->MassK0Short()> 0.55 || v0->MassK0Short()< 0.45) continue;
    fHistEventV0->Fill(9);    
    // if(TMath::Abs((v0->MassLambda() - massLambda))< 0.005) continue;
    // if(TMath::Abs((v0->MassAntiLambda() - massLambda))< 0.005) continue;
    fHistEventV0->Fill(10);    

    bool skipV0=kFALSE;
    Int_t labelPrimOrSecV0=0;
    if(fReadMCTruth){
      if (fMCEvent){
	if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg){
	  if(MotherPos->IsPhysicalPrimary()){
	    //c cout << "\n this particle has passed selections: label mother pos selected" << labelMotherPos<< endl;
	    fHistSelectedV0PtMass->Fill(fTreeVariableInvMassK0s,fTreeVariablePtV0, lPercentiles);
	  }	    
	}
      }
    }
    
    
    if(!(v0->Pt()> fminPtV0 && v0->Pt()<fmaxPtV0) )continue;
    fHistEventV0->Fill(11);     
    //try6  for(Int_t j=0; j < NumberFirstParticle; j++){
    //try6    if (v0->Pt() >= fEvt->fReconstructedFirst[j].fPt){
    //try6  	skipV0=kTRUE;
    //try6  	//	  continue;
    //try6    }
    //try6  }
    //try6
    if (v0->Pt()>=ptTriggerMinimoDati) skipV0=kTRUE;
    if (skipV0){
      //c cout << " pT V0 > Pt trigger " <<endl; 
      continue;
    }
    
    fMassV0->Fill(v0->MassK0Short());    

    for (Int_t m =0; m<5;m++){
      if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	fHistMassvsPt_tagli[m]->Fill(v0->MassK0Short(),v0->Pt());      
      }
    }
    fHistMassvsPt_tagli[5]->Fill(v0->MassK0Short(),v0->Pt());      

    NumberSecondParticle++;

    if(fReadMCTruth){
      if (fMCEvent){
	V0PDGCode=0;
	if(PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 &&  PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg){
	  V0PDGCode=PdgMotherPos;
	  if(MotherPos->IsPhysicalPrimary()){
	    isaK0s=1;
	    //c cout <<"selected v0 Pt " <<  v0->Pt() << endl;
	    fHistResolutionV0Pt->Fill(v0->Pt()- MotherPos->Pt(), lPercentiles);
	    fHistResolutionV0Phi->Fill(v0->Phi()- MotherPos->Phi(), lPercentiles);
	    fHistResolutionV0Eta->Fill(v0->Eta()- MotherPos->Eta(), lPercentiles);
	    fHistSelectedV0PtPhi[5]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
	    fHistSelectedV0PtEta[5]->Fill(v0->Pt(), v0->Eta(), lPercentiles);


	    if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
	      fHistSelectedV0PtPhi[0]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
	      fHistSelectedV0PtEta[0]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      
	      if(v0->CosPointingAngle(lBestPrimaryVtxPos) > 0.997){
		fHistSelectedV0PtPhi[1]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
		fHistSelectedV0PtEta[1]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      } 
	      if(kctau[ParticleType]<0.4*kctauval[ParticleType]){
		fHistSelectedV0PtPhi[2]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
		fHistSelectedV0PtEta[2]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      } 
	      if(TMath::Abs(rapidityV0[0])<0.5){
		fHistSelectedV0PtPhi[3]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
		fHistSelectedV0PtEta[3]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      } 
	      if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.010 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.010) {
		fHistSelectedV0PtPhi[4]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
		fHistSelectedV0PtEta[4]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      } 
	      // if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
	      //   fHistSelectedV0PtPhi[5]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
	      //   fHistSelectedV0PtEta[5]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      // } 
	      if(v0Dca< 0.3){
		fHistSelectedV0PtPhi[6]->Fill(v0->Pt(), v0->Phi(), lPercentiles);
		fHistSelectedV0PtEta[6]->Fill(v0->Pt(), v0->Eta(), lPercentiles);
	      } 
	    }
	       
	    labelPrimOrSecV0=1;
	    NumberSecondParticleRecoTrue++;
	  }
	  else if(MotherPos->IsSecondaryFromWeakDecay())      labelPrimOrSecV0=2;
	  else if(MotherPos->IsSecondaryFromMaterial())      labelPrimOrSecV0=3;
	  else labelPrimOrSecV0=4;
	    
	  for (Int_t m =0; m<5;m++){
	    if(lPercentiles>=moltep[m] && lPercentiles<moltep[m+1]){
	      for(Int_t p=1; p<=4; p++){
		if (labelPrimOrSecV0==p) {
		  fHistPrimaryV0[m][5]->Fill(p, v0->Pt());   
		  if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
		    fHistPrimaryV0[m][0]->Fill(p, v0->Pt());   
		    if(v0->CosPointingAngle(lBestPrimaryVtxPos) > 0.997){
		      fHistPrimaryV0[m][1]->Fill(p, v0->Pt());
		    } 
		    if(kctau[ParticleType]<0.4*kctauval[ParticleType]){
		      fHistPrimaryV0[m][2]->Fill(p, v0->Pt());
		    } 
		    if(TMath::Abs(rapidityV0[0])<0.5){
		      fHistPrimaryV0[m][3]->Fill(p, v0->Pt());
		    } 
		    if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.010 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.010) {
		      fHistPrimaryV0[m][4]->Fill(p, v0->Pt());
		    } 
		    // if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
		    //   fHistPrimaryV0[m][5]->Fill(p, v0->Pt());
		    // } 
		    if(v0Dca< 0.3){
		      fHistPrimaryV0[m][6]->Fill(p, v0->Pt());
		    } 
		  }
		}
	      }
	    }
	  }
	  for(Int_t p=1; p<=4; p++){
	    if (labelPrimOrSecV0==p){
	      fHistPrimaryV0[5][5]->Fill(p, v0->Pt());
	      if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
		fHistPrimaryV0[5][0]->Fill(p, v0->Pt());
		if(v0->CosPointingAngle(lBestPrimaryVtxPos) > 0.997){
		  fHistPrimaryV0[5][1]->Fill(p, v0->Pt());
		} 
		if(kctau[ParticleType]<0.4*kctauval[ParticleType]){
		  fHistPrimaryV0[5][2]->Fill(p, v0->Pt());
		} 
		if(TMath::Abs(rapidityV0[0])<0.5){
		  fHistPrimaryV0[5][3]->Fill(p, v0->Pt());
		} 
		if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.010 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.010) {
		  fHistPrimaryV0[5][4]->Fill(p, v0->Pt());
		} 
		// if(TMath::Abs((v0->MassLambda() - massLambda))>= 0.005 && TMath::Abs((v0->MassAntiLambda() - massLambda))>= 0.005) {
		// 	fHistPrimaryV0[5][5]->Fill(p, v0->Pt());
		// } 
		if(v0Dca< 0.3){
		  fHistPrimaryV0[5][6]->Fill(p, v0->Pt());
		} 
	      }
	    }
	  }
	}
      }
    }

    //save second particle information (V0)
    //cout << "save second particle information (V0) "<< endl;
    if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)){
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sLabelMotherPos = labelMotherPos;
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sLabelPos = labelPos;
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sLabelNeg = labelNeg;
      fEvt->fReconstructedSecond[NumberSecondParticle-1].isMCptc       = isaK0s;
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sDcaPosV0     = v0->DcaPosToPrimVertex();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sDcaNegV0     = v0->DcaNegToPrimVertex();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sPtArmV0      = v0->PtArmV0();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sAlphaV0      = v0->AlphaV0();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sInvMassK0s   = v0->MassK0Short();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sInvMassLambda   = v0->MassLambda();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sInvMassAntiLambda = v0->MassAntiLambda();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sCosPointingAngle  = v0->CosPointingAngle(lBestPrimaryVtxPos);
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sDcaV0ToPV    = v0->DcaV0ToPrimVertex();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sRap          = rapidityV0[0];
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sPt           = v0->Pt();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sctau         = kctau[ParticleType];
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sEta          = v0->Eta();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sTheta        = v0->Theta();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sPhi          = v0->Phi();
      fEvt->fReconstructedSecond[NumberSecondParticle-1].isP           = labelPrimOrSecV0;
      fEvt->fReconstructedSecond[NumberSecondParticle-1].sPDGcode      = V0PDGCode;
      
      fHistPtV0->Fill(v0->Pt());
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
	if ((particleV0->Pt())>ptTriggerMinimoMC) skipV0_MC=kTRUE;
	if (skipV0_MC)      continue;
	fHistPtV0->Fill(particleV0->Pt());
	NumberSecondParticleMC++;
	if(isEfficiency) continue;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDcaPosV0     = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDcaNegV0     = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPtArmV0      = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sAlphaV0      = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sInvMassK0s   = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sInvMassLambda   = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sInvMassAntiLambda = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sCosPointingAngle  = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sDcaV0ToPV    = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sRap          = particleV0->Y();
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPt           = particleV0->Pt();
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sctau         = 0;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sEta          = particleV0->Eta();
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sTheta        = particleV0->Theta();
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPhi          = particleV0->Phi();
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].isP           = 1;
	fEvt->fReconstructedSecond[NumberSecondParticleMC-1].sPDGcode      = particleV0->GetPdgCode();
      }
    }
  }//end loop for second MC particle

  //c cout <<"candidati second particle nell'evento analizzato " <<  NumberSecondParticle << endl;      
  //c cout <<"candidati second particle nell'evento analizzato, MC truth " <<  NumberSecondParticleMC << endl;      

  fHistEventV0->AddBinContent(12, NumberSecondParticle);    
  fHistEventV0->AddBinContent(13, NumberSecondParticleRecoTrue);    
  fHistEventV0->AddBinContent(14, NumberSecondParticleMC);    
  fHistMultvsV0All->Fill(NumberSecondParticle,lPercentiles);
  fHistMultvsV0AllTruth->Fill(NumberSecondParticleRecoTrue,lPercentiles);
  fHistMultvsV0MCAll->Fill(NumberSecondParticleMC,lPercentiles);
  fHistMultvsTriggerAll->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruthAll->Fill(NumberFirstParticleMC, lPercentiles);
  fHistSecondParticleAll->Fill(NumberSecondParticle, NumberSecondParticleMC);
  fHistSecondParticleTruthAll->Fill(NumberSecondParticleRecoTrue, NumberSecondParticleMC);
  if(NumberSecondParticleMC==0 && NumberSecondParticle!=0) fHistEventMult->Fill(17); 
  if(NumberSecondParticleMC < NumberSecondParticleRecoTrue) fHistEventMult->Fill(18); 

  Bool_t DoubleCounted=kFALSE;
  Bool_t CommonDaughtPos=kFALSE;
  Bool_t CommonDaughtNeg=kFALSE;
  if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)){
    for(Int_t j=0; j < NumberSecondParticle; j++){
      if(  fEvt->fReconstructedSecond[j].isP           ==1){
	// fHistEventMult->Fill(20);
	for(Int_t i=0; i< NumberSecondParticle; i++){
	  if(  fEvt->fReconstructedSecond[i].isP           ==1  && i!=j){
	    if(    fEvt->fReconstructedSecond[j].sLabelMotherPos ==     fEvt->fReconstructedSecond[i].sLabelMotherPos) {
	      DoubleCounted=kTRUE;

	      if(    fEvt->fReconstructedSecond[j].sLabelPos ==     fEvt->fReconstructedSecond[i].sLabelPos) CommonDaughtPos=kTRUE;
	      if(    fEvt->fReconstructedSecond[j].sLabelNeg ==     fEvt->fReconstructedSecond[i].sLabelNeg) CommonDaughtNeg=kTRUE;
	    }
	  }
	}
      }
    }
  }
  if (DoubleCounted)  fHistEventMult->Fill(19);
  if (CommonDaughtPos)  fHistEventMult->Fill(20);
  if (CommonDaughtNeg)  fHistEventMult->Fill(21);
  if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)) NumberSecondParticleAll=NumberSecondParticle;
  else if (fReadMCTruth && !isEfficiency) NumberSecondParticleAll=NumberSecondParticleMC;

  if(NumberSecondParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    //c cout << "event has trigger particle but no V0 " << endl;
    PostData(3,fBkgTree);
    PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     
    return;
  }




  fHistSecondParticle->Fill(NumberSecondParticle, NumberSecondParticleMC);
  fHistSecondParticleTruth->Fill(NumberSecondParticleRecoTrue, NumberSecondParticleMC);

  fHistEventV0->AddBinContent(15, NumberSecondParticle);    
  fHistEventV0->AddBinContent(16, NumberSecondParticleRecoTrue);    
  fHistEventV0->AddBinContent(17, NumberSecondParticleMC);    
  // FifoShiftok=kTRUE;
  // LastzBin=zBin;
  // LastcentralityBin=centralityBin;
  
  fHistEventMult->Fill(22);  
  if(!fReadMCTruth || (fReadMCTruth && isEfficiency)){
    fEvt->fNumberCandidateFirst = NumberFirstParticle;
    fEvt->fNumberCandidateSecond = NumberSecondParticle;
  }
  else{
    fEvt->fNumberCandidateFirst = NumberFirstParticleMC;
    fEvt->fNumberCandidateSecond = NumberSecondParticleMC;
  }

  //Fill histos about selected events (events with which I do se and ME) (wt >0 Trigger part. (reco) and >0 V0)  
  fHistTriggerwV0->Fill(NumberFirstParticle);
  fHistTriggerwV0MCTruth->Fill(NumberFirstParticleMC);
  fHistMultvsV0->Fill(NumberSecondParticle,lPercentiles);
  fHistMultvsV0Truth->Fill(NumberSecondParticleRecoTrue,lPercentiles);
  fHistMultvsV0MC->Fill(NumberSecondParticleMC,lPercentiles);
  fHist_multiplicity->Fill(lPercentiles);
  fHistZvertex->Fill(lBestPrimaryVtxPos[2]);
  fHistTrack->AddBinContent(14, NumberFirstParticle);
  fHistTrack->AddBinContent(15, NumberFirstParticleMC);
  fHistMultvsTrigger->Fill(NumberFirstParticle, lPercentiles);
  fHistMultvsTriggerMCTruth->Fill(NumberFirstParticleMC, lPercentiles);
  fHistMultiplicityVsVertexZ->Fill(lBestPrimaryVtxPos[2], lPercentiles);
 
  for(Int_t i=0; i< 100; i++){
    if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMult->AddBinContent(i+1,NumberFirstParticle);  
  }
  for(Int_t i=0; i< 100; i++){
    if(( lPercentiles <i+1) && (lPercentiles >= i) ) fHistTriggervsMultMC->AddBinContent(i+1,NumberFirstParticleMC);  
  }


  for(Int_t l=0; l< NumberFirstParticleAll; l++){
    fHistPt->Fill(fEvt->fReconstructedFirst[l].fPt);
    fHistPtvsMult->Fill(fEvt->fReconstructedFirst[l].fPt, lPercentiles);
    fHist_eta_phi->Fill(fEvt->fReconstructedFirst[l].fPhi, fEvt->fReconstructedFirst[l].fEta);
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
  DoPairsh1h2((Int_t)lPercentiles, fieldsign, lBestPrimaryVtxPos[2]);  

  PostData(1, fOutputList);     
  PostData(2,fSignalTree);
  PostData(3,fBkgTree);
  PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     

  fEventColl[zBin][centralityBin]->FifoShift();       
}




//----------------------------------------------------------------------------------------------------

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

  //c cout << "fNumberCandidateFirst " << fEvt->fNumberCandidateFirst << "  fNumberCandidateSecond " <<  fEvt->fNumberCandidateSecond << endl; 
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
        

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      // For same event pairs
      //cout << "numero dell'evento con cui mixo" <<  eventNumber << " "<<((fEvt+eventNumber)->fNumberCandidateSecond)<< endl;
      //if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateFirst)!=0.) evmultmixed++; 
      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateSecond)!=0.) {
	//c cout << "con questo evento faccio mixing " << endl;
	//c cout << evmultmixed << endl;
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
	  fTreeVariablectau		      = fEvt->fReconstructedSecond[j].sctau;     
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
	  fTreeVariablectau		      = fEvt->fReconstructedSecond[j].sctau;     
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
