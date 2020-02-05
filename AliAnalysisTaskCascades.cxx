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
#include "AliAnalysisTaskCascades.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliAODcascade.h"
//#include "AliMultSelectionBase.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliAODv0.h"
#include "AliVVertex.h"

//class AliESDTrack;
class AliAnalysisTaskCascades;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCascades) // classimp: necessary for root

AliAnalysisTaskCascades::AliAnalysisTaskCascades() :AliAnalysisTaskSE(), 
  fAnalysisType("AOD"), 
  fCollidingSystem("pp"), 
  fAOD(0), 
  fPIDResponse(0),
//fMultSelection(0),
  fEventCuts(0), 			  			
  fOutputList(0), 
  fSignalTree(0), 
  fOutputList2(0),
  fOutputList3(0),  
  fMCEvent(0), 
  fReadMCTruth(0),
  fIshhCorr(0),
  isEfficiency(0),
  fzVertexBins(10), 
  fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(150),
  fnEventsToMix(50),
  fEtaTrigger(0.8),
  fEtahAssoc(0.8),
  fHistPt(0), 
						    fHistPtTriggerParticle(0), 
						    fHistDCAxym1(0),
  fHistDCAzm1(0),
  fHistDCAxym2(0),
  fHistDCAzm2(0),
  fHistPtV0(0), 
  fHistPthAssoc(0), 
  fHistPtTMaxBefAllCfrDataMC(0),
  fHistPtTMinBefAll(0),
  fHistPtTMinBefAllMC(0),
  fHistPtTMaxBefAll(0),
  fHistPtTMaxBefAllBis(0),
  fHistPtTMaxBefAllMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistPtMaxvsMult(0), 
  fHistPtMaxvsMultBefAll(0), 

  fHistZvertex(0),  
  fHistNumberChargedAllEvents(0),
  fHistNumberChargedNoTrigger(0),
  fHistNumberChargedTrigger(0),
  fHist_eta_phi(0),  
  fHist_eta_phi_PtMax(0),  
  fHist_multiplicity(0),
  fHist_multiplicity_EvwTrigger(0),
  fHistEventMult(0), 
  fRunNumber(0),
  fBunchCrossNumber(0),
  fHistEventV0(0), 
  fHistEventXiTrueNeg(0), 
  fHistEventXiTruePos(0), 
  fHistTrack(0), 
  fHistTriggerComposition(0), 
  fHistTriggerCompositionMCTruth(0), 
  fHistAssocComposition(0), 
  fHistAssocCompositionMCTruth(0), 
  fHistTrackAssoc(0), 
  fHistPDG(0),
  fHistPDGLambda(0),
  fHistPDGBachMom(0),
  fHistTheta(0),
  fHistEta(0),
  fHistPhi(0),
  fHistTrackBufferOverflow(0), 
  fHistSecondParticleAll(0),
  fHistSecondParticleTruthAll(0),
  fHistSecondParticle(0),
  fHistSecondParticleTruth(0),
  fMassXiPlus(0), 
  fMassXiMinus(0),
  fV0Lifetime(0), 
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0),
  fHistMultvsV0(0),
  fHistMultvsV0Truth(0),
  fHistMultvsV0MC(0),
  fHistTriggerNotLeading(0),
  fHistTriggerNotLeadingMC(0),
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
  fHistGeneratedXiPt(0),
  fHistSelectedXiPt(0),
  fHistGeneratedOmegaPt(0),
  fHistSelectedOmegaPt(0),
  fHistReconstructedV0PtMass(0),
  fHistSelectedV0PtMass(0),
  fHistResolutionTriggerPt(0),
  fHistResolutionTriggerPhi(0),
  fHistResolutionTriggerEta(0),
  fHistResolutionXiPt(0),
  fHistResolutionXiPhi(0),
  fHistResolutionXiEta(0),
  fHistResolutionOmegaPt(0),
  fHistResolutionOmegaPhi(0),
  fHistResolutionOmegaEta(0),
  fHistResolutionTriggerPhiPt(0),
  fHistResolutionTriggerPhiPdgCode(0),
  fHistPrimaryTrigger(0),
  fHistPrimaryV0(0),
  fminPtj(2), 
  fmaxPtj(10), 
  fV0("kK0s"),  
  fminPtV0(0), 
  fmaxPtV0(30),  
  fminPthAssoc(0), 
  fmaxPthAssoc(30),    
  Evcounter(0), 
  Evcounterczero(0),
  fmolt(5),
  farrGT(0), 
  fTrackBufferSize(20200),
  fTreeVariableMultiplicity(0),       
  fTreeVariableZvertex(0),             
  fTreeVariablePDGCode(0),             
  fTreeVariablePDGCodeBach(0),             
  fTreeVariablePDGCodeNeg(0),             
  fTreeVariablePDGCodePos(0),             
  fTreeVariablePDGCodeLambda(0),             
  fTreeVariablePDGCodeMotherLambda(0),             
  fTreeVariableRunNumber(0),           
  fTreeVariableBunchCrossNumber(0),    
  fTreeVariableNegNSigmaPion(0),
  fTreeVariableNegNSigmaProton(0),
  fTreeVariablePosNSigmaPion(0),
  fTreeVariablePosNSigmaProton(0),
  fTreeVariableBachNSigmaPion(0),
  fTreeVariableBachNSigmaKaon(0),
  fTreeVariablePosTrackLength(0),
  fTreeVariableNegTrackLength(0),
  fTreeVariableBachTrackLength(0),
  fTreeVariableDcaXiToPrimVertex(0),
  fTreeVariableXYDcaXiToPrimVertex(0),
  fTreeVariableZDcaXiToPrimVertex(0),
  fTreeVariableDcaV0ToPrimVertex(0),
  fTreeVariableDcaPosToPrimVertex(0), 
  fTreeVariableDcaNegToPrimVertex(0), 
  fTreeVariableDcaV0Daughters(0), 
  fTreeVariableDcaCascDaughters(0), 
  fTreeVariableDcaBachToPrimVertex(0), 
  fTreeVariableV0CosineOfPointingAngle(0),
  fTreeVariableV0CosineOfPointingAngleSpecial(0),
  fTreeVariableCascCosineOfPointingAngle(0),
  fTreeVariablePtCasc(0),              
  fTreeVariableChargeCasc(0),
  fTreeVariableEtaCasc(0),
  fTreeVariablePhiCasc(0),
  fTreeVariableThetaCasc(0),                    
  fTreeVariablectau(0),               
  fTreeVariableInvMassXi(0),
  fTreeVariableInvMassOmega(0),      
  fTreeVariableInvMassLambda(0),  
  fTreeVariableInvMassK0Short(0),
  fTreeVariableRapXi(0),               
  fTreeVariableRapOmega(0),            
  fTreeVariableCascRadius(0),     
  fTreeVariableV0Radius(0),     
  fTreeVariableLeastNbrClusters(0),    
  fTreeVariableV0Lifetime(0),
  fTreeVariableIsPrimaryXi(0),
  fTreeVariableIsPrimaryOmega(0),
  FifoShiftok(kFALSE)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purpos, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskCascades::AliAnalysisTaskCascades(const char* name) : AliAnalysisTaskSE(name),
								     fAnalysisType("AOD"), 
								     fCollidingSystem("pp"), 
								     fAOD(0), 
								     fPIDResponse(0),
								     //								 								 fMultSelection(0), 			  		
								     fEventCuts(0),								 
								     fOutputList(0), 
								     fSignalTree(0), 
								     fOutputList2(0), 
								     fOutputList3(0), 
								     fMCEvent(0), 
								     fReadMCTruth(0),
								     fIshhCorr(0),
								     isEfficiency(0),
								     fzVertexBins(10), 
								     fnMultBins(20),	 
  fMaxFirstMult(50),
  fMaxSecondMult(150),
  fnEventsToMix(50),
  fEtaTrigger(0.8),
  fEtahAssoc(0.8),
  fHistPt(0), 
  fHistPtTriggerParticle(0), 
  fHistDCAxym1(0),
  fHistDCAzm1(0),
  fHistDCAxym2(0),
  fHistDCAzm2(0),
  fHistPtV0(0), 
  fHistPthAssoc(0), 
  fHistPtTMaxBefAllCfrDataMC(0), 
  fHistPtTMinBefAll(0),
  fHistPtTMinBefAllMC(0),
  fHistPtTMaxBefAll(0),
  fHistPtTMaxBefAllBis(0),
  fHistPtTMaxBefAllMC(0),
  fHistPtvsMult(0), 
  fHistPtvsMultBefAll(0), 
  fHistPtMaxvsMult(0), 
  fHistPtMaxvsMultBefAll(0), 
  fHistZvertex(0),  
  fHistNumberChargedAllEvents(0),
  fHistNumberChargedNoTrigger(0),
  fHistNumberChargedTrigger(0),
  fHist_eta_phi(0),  
  fHist_eta_phi_PtMax(0),  
  fHist_multiplicity(0),
  fHist_multiplicity_EvwTrigger(0),
  fRunNumber(0),
  fBunchCrossNumber(0),
  fHistEventMult(0), 
  fHistEventV0(0), 
  fHistEventXiTrueNeg(0), 
  fHistEventXiTruePos(0), 
  fHistTrack(0), 
  fHistTriggerComposition(0), 
  fHistTriggerCompositionMCTruth(0), 
  fHistAssocComposition(0), 
  fHistAssocCompositionMCTruth(0), 
  fHistTrackAssoc(0), 
  fHistPDG(0), 
  fHistPDGLambda(0),
  fHistPDGBachMom(0),
  fHistTheta(0),
  fHistEta(0),
  fHistPhi(0),
  fHistTrackBufferOverflow(0), 
  fHistSecondParticleAll(0),
  fHistSecondParticleTruthAll(0),
  fHistSecondParticle(0),
  fHistSecondParticleTruth(0),
  fMassXiPlus(0), 
  fMassXiMinus(0),
  fV0Lifetime(0), 
  fHistMultvsV0All(0),
  fHistMultvsV0AllTruth(0),
  fHistMultvsV0MCAll(0), 
  fHistMultvsV0(0),
  fHistMultvsV0Truth(0),
  fHistMultvsV0MC(0),
  fHistTriggerNotLeading(0),
  fHistTriggerNotLeadingMC(0),
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
  fHistGeneratedXiPt(0),
  fHistSelectedXiPt(0),
  fHistGeneratedOmegaPt(0),
  fHistSelectedOmegaPt(0),
  fHistReconstructedV0PtMass(0),
  fHistSelectedV0PtMass(0),
  fHistResolutionTriggerPt(0),
  fHistResolutionTriggerPhi(0),
  fHistResolutionTriggerEta(0),
  fHistResolutionXiPt(0),
  fHistResolutionXiPhi(0),
  fHistResolutionXiEta(0),
  fHistResolutionOmegaPt(0),
  fHistResolutionOmegaPhi(0),
  fHistResolutionOmegaEta(0),
  fHistResolutionTriggerPhiPt(0),
  fHistResolutionTriggerPhiPdgCode(0),
  fHistPrimaryTrigger(0),
  fHistPrimaryV0(0),
  fminPtj(2), 
  fmaxPtj(10), 
  fV0("kK0s"),  
  fminPtV0(0), 
  fmaxPtV0(30),  
  fminPthAssoc(0), 
  fmaxPthAssoc(30),    
  Evcounter(0), 
  Evcounterczero(0),
  fmolt(5),
  farrGT(0), 
  fTrackBufferSize(20200),
  fTreeVariableMultiplicity(0),       
  fTreeVariableZvertex(0),             
  fTreeVariablePDGCode(0),             
  fTreeVariablePDGCodeBach(0),             
  fTreeVariablePDGCodeNeg(0),             
  fTreeVariablePDGCodePos(0),             
  fTreeVariablePDGCodeLambda(0),             
  fTreeVariablePDGCodeMotherLambda(0),             
  fTreeVariableRunNumber(0),           
  fTreeVariableBunchCrossNumber(0),    
  fTreeVariableNegNSigmaPion(0),
  fTreeVariableNegNSigmaProton(0),
  fTreeVariablePosNSigmaPion(0),
  fTreeVariablePosNSigmaProton(0),
  fTreeVariableBachNSigmaPion(0),
  fTreeVariableBachNSigmaKaon(0),
  fTreeVariablePosTrackLength(0),
  fTreeVariableNegTrackLength(0),
  fTreeVariableBachTrackLength(0),
  fTreeVariableDcaXiToPrimVertex(0),
  fTreeVariableXYDcaXiToPrimVertex(0),
  fTreeVariableZDcaXiToPrimVertex(0),
  fTreeVariableDcaV0ToPrimVertex(0),
  fTreeVariableDcaPosToPrimVertex(0), 
  fTreeVariableDcaNegToPrimVertex(0), 
  fTreeVariableDcaV0Daughters(0), 
  fTreeVariableDcaCascDaughters(0), 
  fTreeVariableDcaBachToPrimVertex(0), 
  fTreeVariableV0CosineOfPointingAngle(0),
  fTreeVariableV0CosineOfPointingAngleSpecial(0),
  fTreeVariableCascCosineOfPointingAngle(0),
  fTreeVariablePtCasc(0),              
  fTreeVariableChargeCasc(0),
  fTreeVariableEtaCasc(0),
  fTreeVariablePhiCasc(0),          
  fTreeVariableThetaCasc(0),          
  fTreeVariablectau(0),               
  fTreeVariableInvMassXi(0),        
  fTreeVariableInvMassOmega(0),      
  fTreeVariableInvMassLambda(0),  
  fTreeVariableInvMassK0Short(0),
  fTreeVariableRapXi(0),               
  fTreeVariableRapOmega(0),            
  fTreeVariableCascRadius(0),     
  fTreeVariableV0Radius(0),     
  fTreeVariableLeastNbrClusters(0),    
  fTreeVariableV0Lifetime(0),
  fTreeVariableIsPrimaryXi(0),
  fTreeVariableIsPrimaryOmega(0),
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
AliAnalysisTaskCascades::~AliAnalysisTaskCascades()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
  if(fSignalTree) {
    delete fSignalTree;     // at the end of your task, it is deleted from memory by calling this function
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

}
  

void AliAnalysisTaskCascades::ProcessMCParticles(Bool_t Generated, Float_t lPercentiles, Bool_t isV0,  Bool_t ishhCorr)
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
    cout << "loop for all generated" << endl;
    // Loop over all primary MC particle
    for(Long_t i = 0; i < AODMCTrackArraybis->GetEntriesFast(); i++) {
      
      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(i));
      if (!particle) continue;
      Float_t isXi=0;
      if(isV0==kFALSE){
	fHistPDG->Fill(particle->GetPdgCode());      
	if((particle->Charge())==0) continue;	
	if(TMath::Abs(particle->Eta())>fEtaTrigger)continue; //I need to select particles within this eta range!
	if (!(particle->IsPhysicalPrimary()))continue; 
      }
      if(isV0==kTRUE){
	if(!ishhCorr){ 
	  if(TMath::Abs(particle->Eta())>0.8) continue;
	  if (!(particle->IsPhysicalPrimary()))continue;
	  if ( ((particle->GetPdgCode())==-3312) || ((particle->GetPdgCode())==3312) ){
	    if (particle->GetPdgCode() == -3312) isXi=0.5;
	    if (particle->GetPdgCode() == 3312) isXi=-0.5;
	    fHistGeneratedXiPt[0]->Fill(particle->Pt(), lPercentiles,isXi );
	    if (TMath::Abs(particle->Y())<=0.5 ) fHistGeneratedXiPt[1]->Fill(particle->Pt(), lPercentiles, isXi);
	  }
	  else if ( ((particle->GetPdgCode())==-3334) || ((particle->GetPdgCode())==3334) ){
	    if (particle->GetPdgCode() == -3334) isXi=0.5;
	    if (particle->GetPdgCode() == 3334) isXi=-0.5;
	    fHistGeneratedOmegaPt[0]->Fill(particle->Pt(), lPercentiles,isXi );
	    if (TMath::Abs(particle->Y())<=0.5 ) fHistGeneratedOmegaPt[1]->Fill(particle->Pt(), lPercentiles, isXi);
	  }
	}
      }      
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskCascades::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use 
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output event
  //


  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];

  fOutputList = new TList();         
  fOutputList->SetOwner(kTRUE);     
  fOutputList2 = new TList();         
  fOutputList2->SetOwner(kTRUE);     
  fOutputList3 = new TList();         
  fOutputList3->SetOwner(kTRUE);     

  fSignalTree= new TTree("fSignalTree","fSignalTree");
  fSignalTree->Branch("fTreeVariableMultiplicity",       &fTreeVariableMultiplicity , "fTreeVariableMultiplicity/D");
  fSignalTree->Branch("fTreeVariableZvertex",              &fTreeVariableZvertex  , "fTreeVariableZvertex/D");
  fSignalTree->Branch("fTreeVariablePDGCode",              &fTreeVariablePDGCode  , "fTreeVariablePDGCode/D");
  fSignalTree->Branch("fTreeVariablePDGCodeBach",              &fTreeVariablePDGCodeBach  , "fTreeVariablePDGCodeBach/D");
  fSignalTree->Branch("fTreeVariablePDGCodePos",              &fTreeVariablePDGCodePos  , "fTreeVariablePDGCodePos/D");
  fSignalTree->Branch("fTreeVariablePDGCodeNeg",              &fTreeVariablePDGCodeNeg  , "fTreeVariablePDGCodeNeg/D");
  fSignalTree->Branch("fTreeVariablePDGCodeLambda",              &fTreeVariablePDGCodeLambda  , "fTreeVariablePDGCodeLambda/D");
  fSignalTree->Branch("fTreeVariablePDGCodeMotherLambda",              &fTreeVariablePDGCodeMotherLambda  , "fTreeVariablePDGCodeMotherLambda/D");
  fSignalTree->Branch("fTreeVariableRunNumber",              &fTreeVariableRunNumber  , "fTreeVariableRunNumber/D");
  fSignalTree->Branch("fTreeVariableBunchCrossNumber",              &fTreeVariableBunchCrossNumber  , "fTreeVariableBunchCrossNumber/D");		
  fSignalTree->Branch("fTreeVariableNegNSigmaPion",&fTreeVariableNegNSigmaPion  ,"fTreeVariableNegNSigmaPion/D");
  fSignalTree->Branch("fTreeVariableNegNSigmaProton",&fTreeVariableNegNSigmaProton,"fTreeVariableNegNSigmaProton/D");
  fSignalTree->Branch("fTreeVariablePosNSigmaPion",&fTreeVariablePosNSigmaPion  ,"fTreeVariablePosNSigmaPion/D");
  fSignalTree->Branch("fTreeVariablePosNSigmaProton",&fTreeVariablePosNSigmaProton,"fTreeVariablePosNSigmaProton/D");
  fSignalTree->Branch("fTreeVariableBachNSigmaPion",&fTreeVariableBachNSigmaPion ,"fTreeVariableBachNSigmaPion/D");
  fSignalTree->Branch("fTreeVariableBachNSigmaKaon",&fTreeVariableBachNSigmaKaon ,"fTreeVariableBachNSigmaKaon/D");
  fSignalTree->Branch("fTreeVariablePosTrackLength", &fTreeVariablePosTrackLength, "fTreeVariablePosTrackLength/D");
  fSignalTree->Branch("fTreeVariableNegTrackLength", &fTreeVariableNegTrackLength, "fTreeVariableNegTrackLength/D");
  fSignalTree->Branch("fTreeVariableBachTrackLength", &fTreeVariableBachTrackLength, "fTreeVariableBachTrackLength/D");
  fSignalTree->Branch("fTreeVariableDcaXiToPrimVertex",  &fTreeVariableDcaXiToPrimVertex 	, "fTreeVariableDcaXiToPrimVertex/D");  
  fSignalTree->Branch("fTreeVariableXYDcaXiToPrimVertex",  &fTreeVariableXYDcaXiToPrimVertex 	, "fTreeVariableXYDcaXiToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableZDcaXiToPrimVertex",  &fTreeVariableZDcaXiToPrimVertex 	, "fTreeVariableZDcaXiToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaV0ToPrimVertex",  &fTreeVariableDcaV0ToPrimVertex 	, "fTreeVariableDcaV0ToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaPosToPrimVertex", &fTreeVariableDcaPosToPrimVertex	, "fTreeVariableDcaPosToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaNegToPrimVertex", &fTreeVariableDcaNegToPrimVertex	, "fTreeVariableDcaNegToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableDcaV0Daughters", &fTreeVariableDcaV0Daughters, "fTreeVariableDcaV0Daughters/D");
  fSignalTree->Branch("fTreeVariableDcaCascDaughters", &fTreeVariableDcaCascDaughters, "fTreeVariableDcaCascDaughters/D");
  fSignalTree->Branch("fTreeVariableDcaBachToPrimVertex", &fTreeVariableDcaBachToPrimVertex, "fTreeVariableDcaBachToPrimVertex/D");
  fSignalTree->Branch("fTreeVariableV0CosineOfPointingAngle", &fTreeVariableV0CosineOfPointingAngle	, "fTreeVariableV0CosineOfPointingAngle/D");
  fSignalTree->Branch("fTreeVariableV0CosineOfPointingAngleSpecial", &fTreeVariableV0CosineOfPointingAngleSpecial	, "fTreeVariableV0CosineOfPointingAngleSpecial/D");
  fSignalTree->Branch("fTreeVariableCascCosineOfPointingAngle", &fTreeVariableCascCosineOfPointingAngle	, "fTreeVariableCascCosineOfPointingAngle/D");
  fSignalTree->Branch("fTreeVariablePtCasc",               &fTreeVariablePtCasc   , "fTreeVariablePtCasc/D");
  fSignalTree->Branch("fTreeVariableChargeCasc",               &fTreeVariableChargeCasc   , "fTreeVariableChargeCasc/D");
  fSignalTree->Branch("fTreeVariableEtaCasc",               &fTreeVariableEtaCasc   , "fTreeVariableEtaCasc/D");
  fSignalTree->Branch("fTreeVariablePhiCasc",               &fTreeVariablePhiCasc   , "fTreeVariablePhiCasc/D");
  fSignalTree->Branch("fTreeVariableThetaCasc",               &fTreeVariableThetaCasc   , "fTreeVariableThetaCasc/D");
  fSignalTree->Branch("fTreeVariablectau",               &fTreeVariablectau   , "fTreeVariablectau/D");
  fSignalTree->Branch("fTreeVariableInvMassXi",         &fTreeVariableInvMassXi, "fTreeVariableInvMassXi/D");
  fSignalTree->Branch("fTreeVariableInvMassOmega",      &fTreeVariableInvMassOmega, "fTreeVariableInvMassOmega/D");
  fSignalTree->Branch("fTreeVariableInvMassLambda",  &fTreeVariableInvMassLambda, "fTreeVariableInvMassLambda/D");
  fSignalTree->Branch("fTreeVariableInvMassK0Short",  &fTreeVariableInvMassK0Short, "fTreeVariableInvMassK0Short/D");
  fSignalTree->Branch("fTreeVariableRapXi",               &fTreeVariableRapXi   , "fTreeVariableRapXi/D");
  fSignalTree->Branch("fTreeVariableRapOmega",               &fTreeVariableRapOmega   , "fTreeVariableRapOmega/D");
  fSignalTree->Branch("fTreeVariableCascRadius",         &fTreeVariableCascRadius, "fTreeVariableCascRadius/D");
  fSignalTree->Branch("fTreeVariableV0Radius",      &fTreeVariableV0Radius, "fTreeVariableV0Radius/D");
  fSignalTree->Branch("fTreeVariableLeastNbrClusters",      &fTreeVariableLeastNbrClusters, "fTreeVariableLeastNbrClusters/D");
  fSignalTree->Branch("fTreeVariableV0Lifetime", &fTreeVariableV0Lifetime, "fTreeVariableV0Lifetime/D");
  fSignalTree->Branch("fTreeVariableIsPrimaryXi", &fTreeVariableIsPrimaryXi, "fTreeVariableIsPrimaryXi/D");
  fSignalTree->Branch("fTreeVariableIsPrimaryOmega", &fTreeVariableIsPrimaryOmega, "fTreeVariableIsPrimaryOmega/D");

  fHistPt = new TH1F("fHistPt", "p_{T} distribution of selected charged tracks in events used for AC", 300, 0, 30); 
  fHistPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtTriggerParticle = new TH1F("fHistPtTriggerParticle", "Particle-p_{T} distribution of selected charged tracks in events used for AC", 300, 0, 30); 
  fHistPtTriggerParticle->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 
  fHistDCAxym1 = new TH1F("fHistDCAxym1", "DCAxy method 1 before DCA cuts", 100, -10, 10); 
  fHistDCAxym1->GetXaxis()->SetTitle("DCAxy method 1 (cm)");
  fHistDCAxym2 = new TH1F("fHistDCAxym2", "DCAxy method 2 before DCA cuts", 100, -10, 10); 
  fHistDCAxym2->GetXaxis()->SetTitle("DCAxy method 2 (cm)");
  fHistDCAzm1 = new TH1F("fHistDCAzm1", "DCAz method 1 before DCA cuts", 100, -10, 10); 
  fHistDCAzm1->GetXaxis()->SetTitle("DCAz method 1 (cm)");
  fHistDCAzm2 = new TH1F("fHistDCAzm2", "DCAz method 2 before DCA cuts", 100, -10, 10); 
  fHistDCAzm2->GetXaxis()->SetTitle("DCAz method 2 (cm)");
 
  fHistPtV0 = new TH1F("fHistPtV0", "p_{T} distribution of selected V0 in events used for AC", 300, 0, 30); 
  fHistPtV0->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPthAssoc = new TH1F("fHistPthAssoc", "p_{T} distribution of selected associated charged particles in events used for AC", 300, 0, 30); 
  fHistPthAssoc->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMaxBefAllCfrDataMC= new TH2F("fHistPtTMaxBefAllCfrDataMC", "p_{T} distribution of trigger particle with Pt maximum in the event", 300, 0, 30, 300, 0,30); 
  fHistPtTMaxBefAllCfrDataMC->GetXaxis()->SetTitle("p^{Trigger, Max}_{T} (reco)(GeV/c)");
  fHistPtTMaxBefAllCfrDataMC->GetYaxis()->SetTitle("p^{Trigger, Max}_{T} (MC)(GeV/c)");

  fHistPtTMinBefAll = new TH1F("fHistPtTMinBefAll", "p_{T} distribution of reco trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMinBefAllMC = new TH1F("fHistPtTMinBefAllMC", "p_{T} distribution of true trigger particle with Pt minimum in the event", 300, 0, 30); 
  fHistPtTMinBefAllMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMaxBefAll = new TH1F("fHistPtTMaxBefAll", "p_{T} distribution of reco trigger particle with Pt maximum in the event", 300, 0, 30); 
  fHistPtTMaxBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMaxBefAllBis = new TH1F("fHistPtTMaxBefAllBis", "p_{T} distribution of reco trigger particle with Pt maximum in the event", 300, 0, 30); 
  fHistPtTMaxBefAllBis->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtTMaxBefAllMC = new TH1F("fHistPtTMaxBefAllMC", "p_{T} distribution of true trigger particle with Pt maximum in the event", 300, 0, 30); 
  fHistPtTMaxBefAllMC->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fHistPtvsMultBefAll= new TH2F("fHistPtvsMultBefAll", "p_{T} and centrality distribution of charged tracks in events w T>0", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtvsMult= new TH2F("fHistPtvsMult", "p_{T} and centrality distribution of charged tracks in events used for AC", 300, 0, 30, 100, 0, 100); 
  fHistPtvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtvsMult->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMultBefAll= new TH2F("fHistPtMaxvsMultBefAll", "p_{T} and centrality distribution of charged tracks with maxiumum pt in events w T>0", 300, 0, 30, 100, 0, 100); 
  fHistPtMaxvsMultBefAll->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMultBefAll->GetYaxis()->SetTitle("Centrality");

  fHistPtMaxvsMult= new TH2F("fHistPtMaxvsMult", "p_{T} and centrality distribution of charged tracks with maximum pT in events used for AC)", 300, 0, 30, 100, 0, 100); 
  fHistPtMaxvsMult->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fHistPtMaxvsMult->GetYaxis()->SetTitle("Centrality");


  fHistZvertex= new TH1F("fHistZvertex", "Z vertex distribution of selected events (where I look for Cascades)", 40,-20,20);

  fHistNumberChargedAllEvents=new TH3F("fHistNumberChargedAllEvents", "fHistNumberChargedAllEvents", 100,0,100, 100,0,100, 60, 0,30);
  fHistNumberChargedAllEvents->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedAllEvents->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedAllEvents->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHistNumberChargedNoTrigger=new TH3F("fHistNumberChargedNoTrigger", "fHistNumberChargedNoTrigger", 100,0,100, 100,0,100, 60, 0,30);
  fHistNumberChargedNoTrigger->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedNoTrigger->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedNoTrigger->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHistNumberChargedTrigger=new TH3F("fHistNumberChargedTrigger", "fHistNumberChargedTrigger", 100,0,100, 100,0,100, 60, 0,30);
  fHistNumberChargedTrigger->GetXaxis()->SetTitle("Multiplicity class");
  fHistNumberChargedTrigger->GetYaxis()->SetTitle("Number of charged primary particles");
  fHistNumberChargedTrigger->GetZaxis()->SetTitle("p^{Trigg, Max}_{T} (GeV/c)");

  fHist_eta_phi= new TH2F("fHist_eta_phi", "Distribution of charged tracks in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_eta_phi_PtMax= new TH2F("fHist_eta_phi_PtMax", "Distribution of charged tracks with maximum Pt in events used for AC", 400, 0,2*TMath::Pi(), 400,-0.8, 0.8);
  fHist_eta_phi_PtMax->GetYaxis()->SetTitle("Eta");
  fHist_eta_phi_PtMax->GetXaxis()->SetTitle("Phi (radians)"); 

  fHist_multiplicity=new TH1F("fHist_multiplicity", "fHist_multiplicity", 100, 0, 100); 
  fHist_multiplicity->SetTitle("Centrality distribution of selected events (where I look for Cascades)");
  fHist_multiplicity_EvwTrigger= new TH1F("fHist_multiplicity_EvwTrigger", "fHist_multiplicity_EvwTrigger", 100, 0, 100); 
  fHist_multiplicity_EvwTrigger->SetTitle("Centrality distribution of events with NT>0");

  fHistPDG=new TH1F("fHistPDG", "fHistPDG",6800, -3400, 3400);
  fHistPDGLambda=new TH1F("fHistPDGLambda", "fHistPDGLambda",6800, -3400, 3400);
  fHistPDGBachMom=new TH1F("fHistPDGBachMom", "fHistPDGBachMom",6800, -3400, 3400);
  fHistTheta=new TH1F("fHistTheta", "fHistTheta", 100,0, TMath::Pi());
  fHistEta=new TH1F("fHistEta", "fHistEta", 100,-5, 5);
  fHistPhi=new TH1F("fHistPhi", "fHistPhi", 100,0, 2*TMath::Pi());

  fHistTrackBufferOverflow = new TH1F("fHistTrackBufferOverflow","",2,0,2);
  

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
  fHistEventV0->GetXaxis()->SetBinLabel(3,"All daugthers available");
  fHistEventV0->GetXaxis()->SetBinLabel(4,"Filterbit daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(5,"Chis daughter tracks"); 
  fHistEventV0->GetXaxis()->SetBinLabel(6,"TPC refit"); 
  fHistEventV0->GetXaxis()->SetBinLabel(7,"PID daughters"); 
  fHistEventV0->GetXaxis()->SetBinLabel(8,"|eta daughters|<0.8"); 
  fHistEventV0->GetXaxis()->SetBinLabel(9,"V0 lifetime"); 
  fHistEventV0->GetXaxis()->SetBinLabel(10,"Xi lifetime"); 
  //  fHistEventV0->GetXaxis()->SetBinLabel(9,"more cuts + 0.45 < Mass < 0.55"); 
  //  fHistEventV0->GetXaxis()->SetBinLabel(10,"Lrejection"); 
  fHistEventV0->GetXaxis()->SetBinLabel(11,"Mass selected"); 
  fHistEventV0->GetXaxis()->SetBinLabel(12,"FB4"); 
  fHistEventV0->GetXaxis()->SetBinLabel(13,"FB1"); 
  fHistEventV0->GetXaxis()->SetBinLabel(14,"FB128"); 
  fHistEventV0->GetXaxis()->SetBinLabel(15,"NumberSecondParticle"); 
  fHistEventV0->GetXaxis()->SetBinLabel(16,"NumberTrueXi"); 
  fHistEventV0->GetXaxis()->SetBinLabel(17,"NumberTrueOmega"); 

  fHistEventXiTrueNeg=new TH2F("fHistEventXiTrueNeg", "fHistEventXiTrueNeg",14, 0.5, 14.5, 150,0,30);
  fHistEventXiTrueNeg->SetTitle("Number of V0 which progressively pass the listed selections");
  fHistEventXiTrueNeg->GetYaxis()->SetTitle("p_{T,reco}");
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(1,"All V0s");
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(2,"V0s ok");
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(3,"All daugthers available");
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(4,"Filterbit daughters"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(5,"Chis daughter tracks"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(6,"TPC refit"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(7,"PID daughters"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(8,"|eta daughters|<0.8"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(9,"V0 lifetime"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(10,"Xi lifetime"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(11,"Mass selected"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(12,"FB4"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(13,"FB1"); 
  fHistEventXiTrueNeg->GetXaxis()->SetBinLabel(14,"FB128"); 

  fHistEventXiTruePos=new TH2F("fHistEventXiTruePos", "fHistEventXiTruePos",14, 0.5, 14.5, 150, 0, 30);
  fHistEventXiTruePos->SetTitle("Number of V0 which progressively pass the listed selections");
  fHistEventXiTruePos->GetYaxis()->SetTitle("p_{T,reco}");
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(1,"All V0s");
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(2,"V0s ok");
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(3,"All daugthers available");
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(4,"Filterbit daughters"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(5,"Chis daughter tracks"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(6,"TPC refit"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(7,"PID daughters"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(8,"|eta daughters|<0.8"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(9,"V0 lifetime"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(10,"Xi lifetime"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(11,"Mass selected"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(12,"FB4"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(13,"FB1"); 
  fHistEventXiTruePos->GetXaxis()->SetBinLabel(14,"FB128"); 


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

  fHistTrackAssoc=new TH1F("fHistTrackAssoc", "fHistTrackAssoc", 16, 0.5, 16.5);
  fHistTrackAssoc->GetXaxis()->SetBinLabel(1,"All tracks");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(2,"TrackAssocs after filterbit");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(3,"TrackAssocs with |eta| < 0.8"); 
  fHistTrackAssoc->GetXaxis()->SetBinLabel(4,"TrackAssoc quality");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(5,"TPCCrossedRows>70");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(6,"Crossed rows/findable >0.8");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(7,"Charged tracks");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(8,"DCAxy < 0.010+0.035/pt**1.1");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(9,"DCAz <2");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(10,"0<pT<pTV0max (reco up to here)"); 
  fHistTrackAssoc->GetXaxis()->SetBinLabel(11,"NAssoc(reco) in ev wNT>0"); 
  fHistTrackAssoc->GetXaxis()->SetBinLabel(12,"NAssoc(reco true) in ev wNT>0");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(13,"NAssoc(MC) in ev wNT>0");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(14,"NAssoc(reco) in SelEv");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(15,"NAssoc(reco true) in SelEv");
  fHistTrackAssoc->GetXaxis()->SetBinLabel(16,"NAssoc(MC) in SelEv"); 

  fHistTriggerComposition=new TH2F("fHistTriggerComposition", "fHistTriggerComposition",10000 , -5000, 5000, 2, 0,2);
  fHistTriggerComposition->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");

  fHistTriggerCompositionMCTruth=new TH2F("fHistTriggerCompositionMCTruth", "fHistTriggerCompositionMCTruth",10000 , -5000, 5000, 2, 0,2);
  fHistTriggerCompositionMCTruth->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");

  fHistAssocComposition=new TH2F("fHistAssocComposition", "fHistAssocComposition",10000 , -5000, 5000, 2, 0,2);
  fHistAssocComposition->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");

  fHistAssocCompositionMCTruth=new TH2F("fHistAssocCompositionMCTruth", "fHistAssocCompositionMCTruth",10000 , -5000, 5000, 2, 0,2);
  fHistAssocCompositionMCTruth->GetYaxis()->SetTitle("0=NotPrim, 1=Primary");

  fMassXiPlus= new TH1F("fMassXiPlus", "Invariant mass of XiPlus candidates", 100,1.25, 1.40);
  fMassXiPlus->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}}");

  fMassXiMinus= new TH1F("fMassXiMinus", "Invariant mass of XiMinus candidates", 100,1.25, 1.40);
  fMassXiMinus->GetXaxis()->SetTitle("M_{#Lambda #pi^{-}}");

  fV0Lifetime= new TH1F("fV0Lifetime", "ctau of Lambda daughter of Xi candidates", 100,0, 100);
  fV0Lifetime->GetXaxis()->SetTitle("cm");

  fHistSecondParticleAll= new TH2F("fHistSecondParticleAll", "Number of V0 MCTrue vs number V0 reco (T>0) ", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticleAll->GetXaxis()->SetTitle("Number (reco)");
  fHistSecondParticleAll->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticleTruthAll= new TH2F("fHistSecondParticleTruthAll", "Number of V0 MCTrue vs number V0 reco (true) (T>0)", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticleTruthAll->GetXaxis()->SetTitle("Number (reco true)");
  fHistSecondParticleTruthAll->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticle= new TH2F("fHistSecondParticle", "Number of V0 MCTrue vs number V0 reco (T>0, V>0)", 60,-0.5,59.5,60,-0.5,59.5);
  fHistSecondParticle->GetXaxis()->SetTitle("Number (reco)");
  fHistSecondParticle->GetYaxis()->SetTitle("Number (MC)");

  fHistSecondParticleTruth= new TH2F("fHistSecondParticleTruth", "Number of V0 MCTrue vs number V0 reco (true) (T>0, V>0)", 60,-0.5,59.5,60,-0.5,59.5);
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
  

  fHistMultvsV0=new TH2F("fHistMultvsV0", "Centrality of selected events (where I look for Cascades) vs number of reco XI",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0->GetXaxis()->SetTitle("Number of reco V0 particles");
  fHistMultvsV0->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0Truth=new TH2F("fHistMultvsV0Truth", "Centrality of selected events (T>0, V0>0) vs number of reco true V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0Truth->GetXaxis()->SetTitle("Number of reco true V0 particles");
  fHistMultvsV0Truth->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0MC=new TH2F("fHistMultvsV0MC", "Centrality of selected events (T>0, V0>0) vs number of true V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0MC->GetXaxis()->SetTitle("Number of V0 true particles");
  fHistMultvsV0MC->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0All=new TH2F("fHistMultvsV0All", "Centrality of events w T>0 vs number of reco V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0All->GetXaxis()->SetTitle("Number of V0 particles");
  fHistMultvsV0All->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0AllTruth=new TH2F("fHistMultvsV0AllTruth", "Centrality of events w T>0 vs number of reco true V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0AllTruth->GetXaxis()->SetTitle("Number of V0 reco particles");
  fHistMultvsV0AllTruth->GetYaxis()->SetTitle("Centrality");

  fHistMultvsV0MCAll=new TH2F("fHistMultvsV0MCAll", "Centrality of events w T>0 vs number of true V0s",60, -0.5, 59.5,100, 0, 100 );
  fHistMultvsV0MCAll->GetXaxis()->SetTitle("Number of V0 true particles");
  fHistMultvsV0MCAll->GetYaxis()->SetTitle("Centrality");

  fHistTriggerNotLeading=new TH3F("fHistTriggerNotLeading", "Events with trigger not leading in all events with NT>0",60, -0.5, 59.5,100, 0, 100, 60, 0, 30 );
  fHistTriggerNotLeading->GetXaxis()->SetTitle("Number of V0 with p_{T}> p_{T} Trigger");
  fHistTriggerNotLeading->GetYaxis()->SetTitle("Multiplicity class");
  fHistTriggerNotLeading->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");

  fHistTriggerNotLeadingMC=new TH3F("fHistTriggerNotLeadingMC", "Events with trigger not leading in all events with NT>0 (MC Truth)",60, -0.5, 59.5,100, 0, 100, 60, 0, 30 );
  fHistTriggerNotLeadingMC->GetXaxis()->SetTitle("Number of V0 with p_{T}> p_{T} Trigger");
  fHistTriggerNotLeadingMC->GetYaxis()->SetTitle("Multiplicity class");
  fHistTriggerNotLeadingMC->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");

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

  fHistMultiplicityVsVertexZ=new TH2F("fHistMultiplicityVsVertexZ", "Centrality vs Z vertex of selected events (where I look for Cascades)",  20, -10, 10,100, 0, 100);
      
  fHistTriggervsMult=new TH1F("fHistTriggervsMult", "Numero di particelle di trigger nei vari intervalli di centralita'", 100, 0, 100);
  fHistTriggervsMult->GetXaxis()->SetTitle("Centrality");

  fHistTriggervsMultMC=new TH1F("fHistTriggervsMultMC", "Numero di particelle di trigger (MCtruth) nei vari intervalli di centralita'", 100, 0, 100);
  fHistTriggervsMultMC->GetXaxis()->SetTitle("Centrality");


  fHistGeneratedXiPt=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistGeneratedXiPt[j]=new TH3F(Form("fHistGeneratedXiPt_%i",j), "p_{T} distribution of generated Xi particles (primary)", 300, 0, 30,  100, 0, 100 , 2, -1, 1);
    fHistGeneratedXiPt[j]->GetXaxis()->SetTitle("p_{T}");
    fHistGeneratedXiPt[j]->GetYaxis()->SetTitle("Multiplicity class");
    fHistGeneratedXiPt[j]->GetZaxis()->SetTitle("Xi Minus or Plus");
  }

  fHistSelectedXiPt=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistSelectedXiPt[j]=new TH3F(Form("fHistSelectedXiPt_%i",j), "p_{T} distribution of selected Xi particles (primary)", 300, 0, 30,  100, 0, 100 , 2, -1, 1);
    fHistSelectedXiPt[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedXiPt[j]->GetYaxis()->SetTitle("Multiplicity class");
    fHistSelectedXiPt[j]->GetZaxis()->SetTitle("Xi Minus or Plus");
  }

  fHistGeneratedOmegaPt=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistGeneratedOmegaPt[j]=new TH3F(Form("fHistGeneratedOmegaPt_%i",j), "p_{T} distribution of generated Omega particles (primary)", 300, 0, 30,  100, 0, 100 , 2, -1, 1);
    fHistGeneratedOmegaPt[j]->GetXaxis()->SetTitle("p_{T}");
    fHistGeneratedOmegaPt[j]->GetYaxis()->SetTitle("Multiplicity class");
    fHistGeneratedOmegaPt[j]->GetZaxis()->SetTitle("Omega Minus or Plus");
  }

  fHistSelectedOmegaPt=new TH3F*[2];
  for(Int_t j=0; j<2; j++){
    fHistSelectedOmegaPt[j]=new TH3F(Form("fHistSelectedOmegaPt_%i",j), "p_{T} distribution of selected Omega particles (primary)", 300, 0, 30,  100, 0, 100 , 2, -1, 1);
    fHistSelectedOmegaPt[j]->GetXaxis()->SetTitle("p_{T}");
    fHistSelectedOmegaPt[j]->GetYaxis()->SetTitle("Multiplicity class");
    fHistSelectedOmegaPt[j]->GetZaxis()->SetTitle("Omega Minus or Plus");
  }


  fHistResolutionTriggerPt=new TH3F("fHistResolutionTriggerPt", "p_{T} resolution of selected trigger particles (primary)", 500, -0.5, 0.5, 100, 0, 100, 60,0,30);
  fHistResolutionTriggerPt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionTriggerPt->GetYaxis()->SetTitle("Centrality");
  fHistResolutionTriggerPt->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerPhi=new TH3F("fHistResolutionTriggerPhi", "#Phi resolution of selected trigger particles (primary)", 500, -0.5, 0.5, 100, 0, 100, 60,0,30);
  fHistResolutionTriggerPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhi->GetYaxis()->SetTitle("Centrality");
  fHistResolutionTriggerPhi->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerEta=new TH3F("fHistResolutionTriggerEta", "#Eta resolution of selected trigger particles (primary)", 500, -0.5, 0.5, 100, 0, 100, 60,0,30);
  fHistResolutionTriggerEta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionTriggerEta->GetYaxis()->SetTitle("Centrality");
  fHistResolutionTriggerEta->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionXiPt=new TH2F("fHistResolutionXiPt", "p_{T} resolution of selected Xi particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 60,0,30);
  fHistResolutionXiPt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionXiPt->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionXiPhi=new TH2F("fHistResolutionXiPhi", "#Phi resolution of selected Xi particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 60,0,30);
  fHistResolutionXiPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionXiPhi->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");
  
  fHistResolutionXiEta=new TH2F("fHistResolutionXiEta", "#Eta resolution of selected Xi particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 60,0,30);
  fHistResolutionXiEta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionXiEta->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");


  fHistResolutionOmegaPt=new TH2F("fHistResolutionOmegaPt", "p_{T} resolution of selected Omega particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 60,0,30);
  fHistResolutionOmegaPt->GetXaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionOmegaPt->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionOmegaPhi=new TH2F("fHistResolutionOmegaPhi", "#Phi resolution of selected Omega particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 60,0,30);
  fHistResolutionOmegaPhi->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionOmegaPhi->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");
  
  fHistResolutionOmegaEta=new TH2F("fHistResolutionOmegaEta", "#Eta resolution of selected Omega particles (K0s, primary, event w T>0)", 500, -0.5, 0.5, 60,0,30);
  fHistResolutionOmegaEta->GetXaxis()->SetTitle("Eta_{DATA}-Eta_{MC}");
  fHistResolutionOmegaEta->GetYaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerPhiPt=new TH3F("fHistResolutionTriggerPhiPt", "#Phi and Pt resolution of anomalous selected trigger particles (primary)", 500, -0.5, 0.5, 500, -0.5, 0.5, 60,0,30);
  fHistResolutionTriggerPhiPt->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhiPt->GetYaxis()->SetTitle("p_{T, DATA}-p_{T, MC}");
  fHistResolutionTriggerPhiPt->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

  fHistResolutionTriggerPhiPdgCode=new TH3F("fHistResolutionTriggerPhiPdgCode", "#Phi resolution and pdg code of anomalous selected trigger particles (primary)", 500, -0.5, 0.5,6400, -3200, 3200, 60,0,30);
  fHistResolutionTriggerPhiPdgCode->GetXaxis()->SetTitle("Phi_{DATA}-Phi_{MC}");
  fHistResolutionTriggerPhiPdgCode->GetYaxis()->SetTitle("PdgCode");
  fHistResolutionTriggerPhiPdgCode->GetZaxis()->SetTitle("p^{Trigg, max}_{T}");

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


  fHistPrimaryV0= new TH3F**[6];
  for(Int_t j=0; j<6; j++){
    fHistPrimaryV0[j]=new TH3F*[7];
    for(Int_t i=0; i<7; i++){
      fHistPrimaryV0[j][i]=new TH3F(Form("fHistPrimaryV0_%i_cut%i",j,i), "V0 MC (K0s, selected)", 4, 0.5, 4.5, 160, 0, 16, 60, 0, 30);
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(1,"Primary selected V0s");
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(2,"Secondary from w-decay selected V0s"); 
      fHistPrimaryV0[j][i]->GetXaxis()->SetBinLabel(3,"Secondary from material selected V0s"); 
      fHistPrimaryV0[j][i]->GetYaxis()->SetTitle("p_{T}");
      fHistPrimaryV0[j][i]->GetZaxis()->SetTitle("p^{Trigg, Max}_{T}");
    } 
  }

  //  TString molteplicit[6]={"0-7","7-15","15-25","25-40","40-70",">70"};
  fHistMultiplicityOfMixedEvent=new TH1F("fHistMultiplicityOfMixedEvent", "Distribution of number of events used for the mixing", 20, 0.5, 20.5);

  fEventCuts.AddQAplotsToList(fOutputList);
  
  fOutputList->Add(fHistEventMult);
  fOutputList->Add(fHistEventV0);
  fOutputList->Add(fHistEventXiTrueNeg);
  fOutputList->Add(fHistEventXiTruePos);
  
  //istogrammi riempiti ad ogni evento selezionato (ossia utilizzato per AC) con una entry
  fOutputList->Add(fHistZvertex);
  fOutputList->Add(fHist_multiplicity); 
  fOutputList->Add(fHistMultiplicityVsVertexZ);
  fOutputList->Add(fHistMultvsV0);
  fOutputList->Add(fMassXiPlus);
  fOutputList->Add(fMassXiMinus);
  fOutputList->Add(fV0Lifetime);
  fOutputList->Add(fHistPDG);
  fOutputList->Add(fHistPDGLambda);
  fOutputList->Add(fHistPDGBachMom);
  fOutputList->Add(fHistTheta);
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistPhi);
    
    
  for(Int_t j=0; j < 2; j++){
    fOutputList2->Add(fHistGeneratedXiPt[j]);
    fOutputList2->Add(fHistSelectedXiPt[j]); 
    fOutputList2->Add(fHistGeneratedOmegaPt[j]);
    fOutputList2->Add(fHistSelectedOmegaPt[j]); 
  }
 
  fOutputList->Add(fHistReconstructedV0PtMass);
  fOutputList->Add(fHistSelectedV0PtMass);
  fOutputList->Add(fHistResolutionTriggerPt);
  fOutputList->Add(fHistResolutionTriggerPhi);
  fOutputList->Add(fHistResolutionTriggerEta);
  fOutputList->Add(fHistResolutionXiPt);
  fOutputList->Add(fHistResolutionXiPhi);
  fOutputList->Add(fHistResolutionXiEta);
  fOutputList->Add(fHistResolutionOmegaPt);
  fOutputList->Add(fHistResolutionOmegaPhi);
  fOutputList->Add(fHistResolutionOmegaEta);
  fOutputList->Add(fHistResolutionTriggerPhiPt);
  fOutputList->Add(fHistResolutionTriggerPhiPdgCode);

  PostData(1, fOutputList);  
  PostData(2, fSignalTree);       
  PostData(4, fOutputList2);     
  PostData(5, fOutputList3);     
     
}
//_____________________________________________________________________________
void AliAnalysisTaskCascades::UserExec(Option_t *)
{
  //cout<<"enter user exec"<<endl;
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you 
  // have access to the current event. 
  // once you return from the UserExec function, the manager will retrieve the next event from the chain
  
  //  Float_t moltep[5]={0,7,15,25,40,70}; //valori associati a molteplicita'
  //  cout << " ciao, sto ora iniziando a eseguire lo user exec " << endl; 
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
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    return;
  }        
  
  
  /// Use the event cut class to apply the required selections
  if (!fEventCuts.AcceptEvent(fAOD)) {   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    PostData(4, fOutputList2);     
    PostData(5, fOutputList3);      return;
  }
  
  
  Int_t iTracks(fAOD->GetNumberOfTracks());         
  Int_t V0Tracks(fAOD->GetNumberOfV0s());           
  Int_t iCascades(fAOD->GetNumberOfCascades());           
  Evcounter++;  
  cout << "\n \n \n ********************************************************* "<< endl;
  
  //VERTEX SELECTION AND TRIGGER
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};
  const AliAODVertex *lPrimaryBestAODVtx = fAOD->GetPrimaryVertex();
  if (!lPrimaryBestAODVtx){
    AliWarning("No prim. vertex in AOD... return!");
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
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

  fRunNumber = fAOD->GetRunNumber();
  fBunchCrossNumber = fAOD->GetBunchCrossNumber();
  Double_t lMagneticField = -10;
  lMagneticField = fAOD->GetMagneticField( );

    
  if(!isSelected){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    //c cout << "event does not fulfil centrality selection criteria " << endl;     
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    return;
  }
  

  Int_t NumberSecondParticleRecoTrue=0;
  Int_t NumberSecondParticleRecoTrueXiPos=0;
  Int_t NumberSecondParticleRecoTrueXiNeg=0;
  Int_t NumberSecondParticle=0;
  Int_t NumberTruthXi=0;
  Int_t NumberTruthOmega=0;
  Int_t NumberSecondParticleMC=0;
  Int_t NumberSecondParticleAll=0;

  //*******************************************************************************************
  //********************************Loop for cascades*****************************************
  cout << " loop for cascades to begin " << endl;
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

  AliAODMCParticle* particlePos;
  AliAODMCParticle* particleNeg;
  AliAODMCParticle* particleBach;
  Int_t PdgPos=0;
  Int_t PdgNeg=0;
  Int_t PdgBach=0;
  Int_t labelMotherPos=0;
  Int_t labelMotherNeg=0;
  Int_t labelGMotherPos=0;
  Int_t labelGMotherNeg=0;
  Int_t labelMotherBach=0;
  AliAODMCParticle* MotherPos;
  AliAODMCParticle* MotherNeg;
  AliAODMCParticle* GMotherPos;
  AliAODMCParticle* GMotherNeg;
  AliAODMCParticle* MotherBach;
  Int_t PdgMotherPos=0;
  Int_t PdgMotherNeg=0;
  Int_t PdgMotherLambda=0;
  Int_t PdgGMotherPos=0;
  Int_t PdgGMotherNeg=0;
  Int_t PdgMotherBach=0;
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

  
  //MC generated cascades
  Bool_t  isV0=kTRUE;
  Bool_t  Generated=kTRUE;

  if(fReadMCTruth){
    fMCEvent=MCEvent(); 
    if (fMCEvent){
      ProcessMCParticles(Generated, lPercentiles, isV0, fIshhCorr);
    }
  }
 

  //------------------------------------------------
  // MAIN CASCADE LOOP STARTS HERE
  //------------------------------------------------
  // Code Credit: Antonin Maire (thanks^100)
  // ---> This is an adaptation

  for (Int_t iXi = 0; iXi < iCascades; iXi++) {
    fHistEventV0->Fill(1);  
    //------------------------------------------------
    // Initializations
    //------------------------------------------------
    //Double_t lTrkgPrimaryVtxRadius3D = -500.0;
    //Double_t lBestPrimaryVtxRadius3D = -500.0;

    // - 1st part of initialisation : variables needed to store AliESDCascade data members
    Double_t lEffMassXi      = 0. ;
    //Double_t lChi2Xi         = -1. ;
    Double_t lDcaXiDaughters = -1. ;
    Double_t lXiCosineOfPointingAngle = -1. ;
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
    Double_t lXiRadius = -1000. ;
    Double_t lXiDecayLength = -1000. ;
    Double_t kctauXi = -1000. ;

    // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
    Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
    Int_t    lNegTPCClusters    = -1; // For ESD only ...
    Int_t    lBachTPCClusters   = -1; // For ESD only ...

    // - 3rd part of initialisation : about V0 part in cascades
    Double_t lInvMassLambdaAsCascDghter = 0.;
    Double_t lInvMassK0sAsCascDghter = 0.;
    //Double_t lV0Chi2Xi         = -1. ;
    Double_t lDcaV0DaughtersXi = -1.;

    Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
    Double_t lDcaPosToPrimVertexXi  = -1.;
    Double_t lDcaNegToPrimVertexXi  = -1.;
    Double_t lDcaXiToPrimVertex= -1.;
    Double_t lV0CosineOfPointingAngleXi = -1. ;
    Double_t lV0CosineOfPointingAngleXiSpecial = -1. ;
    Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
    Double_t lV0RadiusXi = -1000.0;
    Double_t lV0quality  = 0.;

    // - 4th part of initialisation : Effective masses
    Double_t lInvMassXiMinus    = 0.;
    Double_t lInvMassXiPlus     = 0.;
    Double_t lInvMassOmegaMinus = 0.;
    Double_t lInvMassOmegaPlus  = 0.;

    // - 6th part of initialisation : extra info for QA
    Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
    Double_t lXiTransvMom  = 0. ;
    //Double_t lXiTransvMomMC= 0. ;
    Double_t lXiTotMom     = 0. ;

    Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
    //Double_t lBachTransvMom  = 0.;
    //Double_t lBachTotMom     = 0.;

    fTreeVariableNegNSigmaPion   = -100;
    fTreeVariableNegNSigmaProton = -100;
    fTreeVariablePosNSigmaPion   = -100;
    fTreeVariablePosNSigmaProton = -100;
    fTreeVariableBachNSigmaPion  = -100;
    fTreeVariableBachNSigmaKaon  = -100;

    Short_t  lChargeXi = -2;
    Double_t  lEtaXi = -999;
    Double_t  lThetaXi = -999;
    Double_t  lPhiXi = -999;
    //Double_t lV0toXiCosineOfPointingAngle = 0. ;

    Double_t lRapXi   = -20.0, lRapOmega = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
    //Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
    Float_t NominalMassXi =1.32171;
    // -------------------------------------
      
    AliAODcascade *xi = fAOD->GetCascade(iXi);
    if (!xi) continue;
    fHistEventV0->Fill(2);  

    //ChiSquare implementation
    /*
      fTreeCascVarChiSquareV0      = xi->Chi2V0();
      fTreeCascVarChiSquareCascade = xi->Chi2Xi();
    */
    //Xi decay vertex 
    lPosXi[0] = xi->DecayVertexXiX();
    lPosXi[1] = xi->DecayVertexXiY();
    lPosXi[2] = xi->DecayVertexXiZ();
    
    //Xi decay length and radius
    lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] ); //calculated wrt 0 since PVx/PVy approx = 0
    
    lXiDecayLength = TMath::Sqrt(
				 TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
				 TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
				 TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
				 );
    lChargeXi = xi->ChargeXi();
    lThetaXi = xi->Theta();
    lEtaXi = -TMath::Log(TMath::Tan(lThetaXi/2));
    lPhiXi = xi->Phi();
    
    //parte valida per esd ma non aod? 
    /*
    //V0 and Bach
    UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
    UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
    UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
    // Care track label can be negative in MC production (linked with the track quality)
    // However = normally, not the case for track index ...

    // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
    if(lBachIdx == lIdxNegXi) {
    AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");
    continue;
    }
    if(lBachIdx == lIdxPosXi) {
    AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");
    continue;
    }
    */

    AliAODTrack *pTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
    AliAODTrack *nTrackXi    = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
    AliAODTrack *bachTrackXi = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );
    if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
      AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
      continue;
    }
    fHistEventV0->Fill(3);   
    
    Int_t labelPos    = pTrackXi->GetLabel();
    Int_t labelNeg    = nTrackXi->GetLabel();
    Int_t labelBach = bachTrackXi->GetLabel();


    //-------------------------------------------------------
    //---------MC information--------------------------------
    //-------------------------------------------------------

    TClonesArray* AODMCTrackArray =0x0;
    if(fReadMCTruth){
      fMCEvent= MCEvent();
      cout << "hey there I'm getting pdg info! " << endl;
      if (fMCEvent){
	cout << "hey there I'm getting pdg info! (2)" << endl;
	AODMCTrackArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArray == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}
	particlePos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelPos)));
	particleNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelNeg)));
	particleBach = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelBach)));
	cout << "label 3 daughters (pos, neg, bach) " << labelPos << "  " << labelNeg << "  " << labelBach << endl;
	// if(labelPos>=0)	PdgPos = particlePos->GetPdgCode();
	// if(labelNeg>=0)	PdgNeg = particleNeg->GetPdgCode();
	// if(labelBach>=0)	PdgBach = particleBach->GetPdgCode();
	PdgPos = particlePos->GetPdgCode();
	PdgNeg = particleNeg->GetPdgCode();
	PdgBach = particleBach->GetPdgCode();

	labelMotherPos=particlePos->GetMother();
	labelMotherNeg=particleNeg->GetMother();
	labelMotherBach=particleBach->GetMother();
	//c cout << "label tracce madri (pos e neg) " << 	labelMotherPos<< endl;
	//c cout << "label tracce madri (pos e neg) " << 	labelMotherNeg<< endl;
	    
	MotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(labelMotherPos)));
	MotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(labelMotherNeg)));
	MotherBach = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelMotherBach)));
	cout << "label 3 mothers (pos, neg, bach) " << labelMotherPos << "  " << labelMotherNeg << "  " << labelMotherBach << endl;	    
	// if (labelMotherPos>=0) PdgMotherPos = MotherPos->GetPdgCode();
	// if (labelMotherNeg>=0)	PdgMotherNeg = MotherNeg->GetPdgCode();
	// if (labelMotherBach>=0)	PdgMotherBach = MotherBach->GetPdgCode();
	PdgMotherPos = MotherPos->GetPdgCode();
	PdgMotherNeg = MotherNeg->GetPdgCode();
	PdgMotherBach = MotherBach->GetPdgCode();
	    
	labelGMotherPos=MotherPos->GetMother();
	labelGMotherNeg=MotherNeg->GetMother();

	GMotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelGMotherPos)));
	GMotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelGMotherNeg)));
	cout << "label 2GMothers (pos, neg, bach) " << labelGMotherPos << "  " << labelGMotherNeg << endl;	   
	// if (labelGMotherPos>=0) PdgGMotherPos = GMotherPos->GetPdgCode();
	// if (labelGMotherNeg>=0)	PdgGMotherNeg = GMotherNeg->GetPdgCode();	
	PdgGMotherPos = GMotherPos->GetPdgCode();
	PdgGMotherNeg = GMotherNeg->GetPdgCode();	
	 	 
      }
    }

    Bool_t isXiNeg=kFALSE;
    Bool_t isXiPos=kFALSE;
    Bool_t isOmegaNeg=kFALSE;
    Bool_t isOmegaPos=kFALSE;
    Float_t PtMotherBach=-999;
    Float_t isXi=-999;
    if(fReadMCTruth){
      if (fMCEvent){

	isXiNeg = (PdgPos==2212 && PdgNeg==-211 && PdgMotherPos == 3122 && PdgMotherNeg == 3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==-211 && PdgMotherBach==3312 && PdgGMotherPos==3312 &&PdgGMotherNeg==3312 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);
	isXiPos = (PdgPos==211 && PdgNeg==-2212 && PdgMotherPos == -3122 &&  PdgMotherNeg == -3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==211 && PdgMotherBach==-3312 && PdgGMotherPos==-3312 && PdgGMotherNeg==-3312 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);
        isOmegaNeg = (PdgPos==2212 && PdgNeg==-211 && PdgMotherPos == 3122 && PdgMotherNeg == 3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==-321 && PdgMotherBach==3334 && PdgGMotherPos==3334 &&PdgGMotherNeg==3334 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);
	isOmegaPos = (PdgPos==211 && PdgNeg==-2212 && PdgMotherPos == -3122 &&  PdgMotherNeg == -3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==321 && PdgMotherBach==-3334 && PdgGMotherPos==-3334 && PdgGMotherNeg==-3334 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);

	PtMotherBach=MotherBach->Pt();
	Bool_t isV0Particle=(labelMotherPos==labelMotherNeg);
	//	    cout << "isXiPos " << isXiPos << "isXiNeg " << isXiNeg << endl;
	if(isXiPos || isXiNeg){
	  if (GMotherPos->IsPrimary()) {fTreeVariableIsPrimaryXi=1;}
	  else {fTreeVariableIsPrimaryXi=-1;}
	  fTreeVariableIsPrimaryOmega=0;
	  //	      fTreeVariablePDGCode=PdgGMotherPos;
	  fTreeVariablePDGCode=PdgMotherBach;
	  fTreeVariablePDGCodeBach=0;
	  fTreeVariablePDGCodeLambda=0;
	  fTreeVariablePDGCodePos=0;
	  fTreeVariablePDGCodeNeg=0;
	  fTreeVariablePDGCodeMotherLambda=0;
	}
	else if(isOmegaPos || isOmegaNeg){
	  if (GMotherPos->IsPrimary()) {fTreeVariableIsPrimaryOmega=1;}
	  else {fTreeVariableIsPrimaryOmega=-1;}
	  fTreeVariableIsPrimaryXi=0;
	  //	      fTreeVariablePDGCode=PdgGMotherPos;
	  fTreeVariablePDGCode=PdgMotherBach;
	  fTreeVariablePDGCodeBach=0;
	  fTreeVariablePDGCodeLambda=0;
	  fTreeVariablePDGCodePos=0;
	  fTreeVariablePDGCodeNeg=0;
	  fTreeVariablePDGCodeMotherLambda=0;
	}

	else {
	  if (isV0Particle) {
	    fTreeVariablePDGCodeLambda=PdgMotherPos; 	  
	    fTreeVariablePDGCodeMotherLambda=PdgGMotherPos;
	  }
	  else {
	    fTreeVariablePDGCodeLambda=0;
	    fTreeVariablePDGCodeMotherLambda=0;
	  }
	  if (isV0Particle) fTreeVariablePDGCodeBach=PdgBach;
	  else  fTreeVariablePDGCodeBach=0;
	  fTreeVariablePDGCodePos=PdgPos;
	  fTreeVariablePDGCodeNeg=PdgNeg;
	  fTreeVariablePDGCode=PdgMotherBach;
	  fTreeVariableIsPrimaryXi=0;
	  fTreeVariableIsPrimaryOmega=0;

	}
      }
    }

    //FilterBit for daugther tracks----------
    // if(!pTrackXi->TestFilterBit(4)) continue;
    // if(!nTrackXi->TestFilterBit(4)) continue;
    // if(!bachTrackXi->TestFilterBit(4)) continue;
    //----------------------------------------
    fHistEventV0->Fill(4);   
    if(fReadMCTruth){
      if (fMCEvent){
	if (isXiPos)    fHistEventXiTruePos->Fill(4,PtMotherBach);  
	if (isXiNeg)    fHistEventXiTrueNeg->Fill(4,PtMotherBach);  
      }
    }
    if (pTrackXi->TestFilterBit(4) && nTrackXi->TestFilterBit(4) &&bachTrackXi->TestFilterBit(4)){
      if(fReadMCTruth){
	if (fMCEvent){

	  if (isXiPos)    fHistEventXiTruePos->Fill(12,PtMotherBach);  
	  if (isXiNeg)    fHistEventXiTrueNeg->Fill(12,PtMotherBach);  
	}
      }
      fHistEventV0->Fill(12);   
    }
    if (pTrackXi->TestFilterBit(1) && nTrackXi->TestFilterBit(1) &&bachTrackXi->TestFilterBit(1)){
      if(fReadMCTruth){
	if (fMCEvent){

	  if (isXiPos)    fHistEventXiTruePos->Fill(13,PtMotherBach);  
	  if (isXiNeg)    fHistEventXiTrueNeg->Fill(13,PtMotherBach);  
	}
      }
      fHistEventV0->Fill(13);   
    }
    if (pTrackXi->TestFilterBit(128) && nTrackXi->TestFilterBit(128) &&bachTrackXi->TestFilterBit(128)){
      if(fReadMCTruth){
	if (fMCEvent){
  
	  if (isXiPos)  fHistEventXiTruePos->Fill(14,PtMotherBach);  
	  if (isXiNeg)    fHistEventXiTrueNeg->Fill(14,PtMotherBach);  
	}
      }
      fHistEventV0->Fill(14);   
    }

    //daughter track quality cuts------------
    if(pTrackXi->Chi2perNDF()>4.)continue;
    if(nTrackXi->Chi2perNDF()>4.)continue;
    if(bachTrackXi->Chi2perNDF()>4.)continue;
    //--------------------------------------
    fHistEventV0->Fill(5);  
    if(fReadMCTruth){
      if (fMCEvent){

	if (isXiPos)    fHistEventXiTruePos->Fill(5,PtMotherBach);  
	if (isXiNeg)    fHistEventXiTrueNeg->Fill(5,PtMotherBach);  
      }
    }

    Double_t lBMom[3], lNMom[3], lPMom[3];
    pTrackXi->GetPxPyPz( lBMom);
    nTrackXi->GetPxPyPz( lPMom);
    bachTrackXi->GetPxPyPz( lNMom);

    Float_t lBachTransMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] );
    Float_t lPosTransMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
    Float_t lNegTransMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );

    //------------------------------------------------
    // TPC dEdx information
    //------------------------------------------------
    fTreeVariableNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion   );
    fTreeVariableNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
    fTreeVariablePosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
    fTreeVariablePosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
    fTreeVariableBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
    fTreeVariableBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );

    //------------------------------------------------
    // TPC Number of clusters info
    // --- modified to save the smallest number
    // --- of TPC clusters for the 3 tracks
    //------------------------------------------------

    lPosTPCClusters   = pTrackXi->GetTPCNcls();
    lNegTPCClusters   = nTrackXi->GetTPCNcls();
    lBachTPCClusters  = bachTrackXi->GetTPCNcls();

    // 1 - Poor quality related to TPCrefit
    ULong_t pStatus    = pTrackXi->GetStatus();
    ULong_t nStatus    = nTrackXi->GetStatus();
    ULong_t bachStatus = bachTrackXi->GetStatus();
    //fTreeVariablekITSRefitBachelor = kTRUE;
    //fTreeVariablekITSRefitNegative = kTRUE;
    //fTreeVariablekITSRefitPositive = kTRUE;

    
    if ((pStatus&AliAODTrack::kTPCrefit)    == 0) {
      AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!");
      continue;
    }
    if ((nStatus&AliAODTrack::kTPCrefit)    == 0) {
      AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!");
      continue;
    }
    if ((bachStatus&AliAODTrack::kTPCrefit) == 0) {
      AliWarning("Pb / Bach.   track has no TPCrefit ... continue!");
      continue;
    }
    fHistEventV0->Fill(6);  
    if(fReadMCTruth){
      if (fMCEvent){
	if (isXiPos)    fHistEventXiTruePos->Fill(6,PtMotherBach);  
	if (isXiNeg)    fHistEventXiTrueNeg->Fill(6,PtMotherBach);  
      }
    }
    //***********************************************
    //Sandbox mode information, please
    //   fTreeCascVarMagneticField = fAOD->GetMagneticField();
    // fTreeCascVarPosDCAz = GetDCAz(pTrackXi);
    // fTreeCascVarNegDCAz = GetDCAz(nTrackXi);
    // fTreeCascVarBachDCAz = GetDCAz(bachTrackXi);
        
    Float_t lPosChi2PerCluster = pTrackXi->GetTPCchi2() / ((Float_t) lPosTPCClusters);
    Float_t lNegChi2PerCluster = nTrackXi->GetTPCchi2() / ((Float_t) lNegTPCClusters);
    Float_t lBachChi2PerCluster = bachTrackXi->GetTPCchi2() / ((Float_t) lBachTPCClusters);
        
     
    //********************************
      
    Int_t leastnumberofclusters = 1000;
    if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
    if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
    if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;

    //Extra track quality: min track length--------------------------------------------
    Float_t lSmallestTrackLength = 1000;
    Float_t lPosTrackLength = -1;
    Float_t lNegTrackLength = -1;
    Float_t lBachTrackLength = -1;
    lPosTrackLength = GetLengthInActiveZone( pTrackXi, /*1,*/ 2.0, 220.0, fAOD->GetMagneticField());
    lNegTrackLength = GetLengthInActiveZone( nTrackXi, /*1,*/ 2.0, 220.0, fAOD->GetMagneticField());
    lBachTrackLength = GetLengthInActiveZone( bachTrackXi, /*1,*/ 2.0, 220.0, fAOD->GetMagneticField());
        
    if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
    if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
    if ( lBachTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lBachTrackLength;

    //Calculate V0 lifetime for adaptive decay radius cut
    xi->GetXYZ( lPosV0Xi );
    lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
        
    Float_t lLeastNcrOverLength = 200;
    Float_t lPosTrackNcrOverLength = pTrackXi->GetTPCClusterInfo(2,1)/(lPosTrackLength-TMath::Max(lV0RadiusXi-85.,0.));
    Float_t lNegTrackNcrOverLength = nTrackXi->GetTPCClusterInfo(2,1)/(lNegTrackLength-TMath::Max(lV0RadiusXi-85.,0.));
    Float_t lBachTrackNcrOverLength = bachTrackXi->GetTPCClusterInfo(2,1)/(lBachTrackLength-TMath::Max(lXiRadius-85.,0.));
        
    lLeastNcrOverLength = (Float_t) lPosTrackNcrOverLength;
    if( lNegTrackNcrOverLength < lLeastNcrOverLength )
      lLeastNcrOverLength = (Float_t) lNegTrackNcrOverLength;
    if( lBachTrackNcrOverLength < lLeastNcrOverLength )
      lLeastNcrOverLength = (Float_t) lBachTrackNcrOverLength;
        

    //        fTreeCascVarMinTrackLength = lSmallestTrackLength;
    //---------------------------------------------------------------------------------------

    // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
    //    if(lPosTPCClusters  < 70 && lSmallestTrackLength<80) {
    if(lPosTPCClusters  < 70) {
      AliWarning("Pb / V0 Pos. track has less than 70 TPC clusters ... continue!");
      continue;
    }
    //    if(lNegTPCClusters  < 70  && lSmallestTrackLength<80) {
    if(lNegTPCClusters  < 70) {
      AliWarning("Pb / V0 Neg. track has less than 70 TPC clusters ... continue!");
      continue;
    }
    //    if(lBachTPCClusters < 70  && lSmallestTrackLength<80) {
    if(lBachTPCClusters < 70) {
      AliWarning("Pb / Bach.   track has less than 70 TPC clusters ... continue!");
      continue;
    }
    fHistEventV0->Fill(7);  
    if(fReadMCTruth){
      if (fMCEvent){

	if (isXiPos)    fHistEventXiTruePos->Fill(7,PtMotherBach);  
	if (isXiNeg)    fHistEventXiTrueNeg->Fill(7,PtMotherBach);  
      }
    }

    //eta of the three daughters 
    Double_t pEta   = pTrackXi->Eta();
    Double_t nEta    = nTrackXi->Eta();
    Double_t bachEta = bachTrackXi->Eta();
    if (TMath::Abs(pEta)>0.8) continue;
    if (TMath::Abs(nEta)>0.8) continue;
    if (TMath::Abs(bachEta)>0.8) continue;
    fHistEventV0->Fill(8);  
    if(fReadMCTruth){
      if (fMCEvent){

	if (isXiPos)    fHistEventXiTruePos->Fill(8,PtMotherBach);  
	if (isXiNeg)    fHistEventXiTrueNeg->Fill(8,PtMotherBach);  
      }
    }

    cout << "lBachTransMom: " << lBachTransMom << "lNegTransMom: " << lNegTransMom <<  "lPosTransMom: " << lPosTransMom << endl;
    //alternative way to define charge    lChargeXi = bachTrackXi->Charge();
    if ( lChargeXi < 0)
      lInvMassLambdaAsCascDghter    = xi->MassLambda();
    else
      lInvMassLambdaAsCascDghter    = xi->MassAntiLambda();


    lInvMassK0sAsCascDghter    = xi->MassK0Short();

    //Different DCA    
    lDcaXiDaughters 	        = xi->DcaXiDaughters();
    lDcaV0DaughtersXi 		= xi->DcaV0Daughters();
    lDcaXiToPrimVertex		= xi->DcaXiToPrimVertex(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2]);
    lDcaV0ToPrimVertexXi        = xi->DcaV0ToPrimVertex();
    lDcaBachToPrimVertexXi      = xi->DcaBachToPrimVertex();
    lDcaPosToPrimVertexXi        = xi->DcaPosToPrimVertex();
    lDcaNegToPrimVertexXi        = xi->DcaNegToPrimVertex();

    //cosine of pointing angle
    lXiCosineOfPointingAngle   = xi->CosPointingAngleXi( lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
    lV0CosineOfPointingAngleXi = xi->CosPointingAngle(lBestPrimaryVtxPos);
    //Modification: V0 CosPA wrt to Cascade decay vertex
    lV0CosineOfPointingAngleXiSpecial = xi->CosPointingAngle(xi->GetDecayVertexXi());
    
    ///========================================================================================
    
    //3D Distance travelled by the V0 in the cascade
    Float_t lV0DistanceTrav =  TMath::Sqrt(  TMath::Power( lPosV0Xi[0]-lPosXi[0] , 2)
					     + TMath::Power( lPosV0Xi[1]-lPosXi[1] , 2)
					     + TMath::Power( lPosV0Xi[2]-lPosXi[2] , 2) );
        
    //Total V0 momentum
    Float_t lV0TotMomentum = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
					   + TMath::Power( lNMom[1]+lPMom[1] , 2)
					   + TMath::Power( lNMom[2]+lPMom[2] , 2) );
        
    //V0 transverse momentum
    Float_t lV0Pt = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
				  + TMath::Power( lNMom[1]+lPMom[1] , 2) );
        
    //Calculate V0 lifetime: mL/p
    //	cout << "V0 total momentum " << lV0TotMomentum << endl;
    Double_t lV0Lifetime=-1;
    if( TMath::Abs(lV0TotMomentum)>1e-5 ){
      lV0Lifetime = 1.115683*lV0DistanceTrav / lV0TotMomentum;
    }else{
      lV0Lifetime = -1;
    }
    if (lV0Lifetime>100) continue;
    fHistEventV0->Fill(9);  
    if(fReadMCTruth){
      if (fMCEvent){

	if (isXiPos)	fHistEventXiTruePos->Fill(9,PtMotherBach);  
	if (isXiNeg)	fHistEventXiTrueNeg->Fill(9,PtMotherBach);  
      }
    }

    /*
    //========================================================================================

    //Cowboy/sailor info regarding V0 inside cascade
    //Calculate vec prod with momenta projected to xy plane
    //Provisions for cowboy/sailor check
    Double_t lModp1 = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
    Double_t lModp2 = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );
        
    Double_t lVecProd = (lPMom[0]*lNMom[1] - lPMom[1]*lNMom[0]) / (lModp1*lModp2);
        
    if ( lMagneticField < 0 ) lVecProd *= -1; //invert sign
        
    fTreeCascVarIsCowboy = kFALSE;
    if (lVecProd < 0) fTreeCascVarIsCowboy = kTRUE;
        
    fTreeCascVarCowboyness = lVecProd;
        
    Double_t lBachMod = TMath::Sqrt(lBMom[0]*lBMom[0]+lBMom[1]*lBMom[1]);
    Double_t lV0px = lPMom[0] + lNMom[0];
    Double_t lV0py = lPMom[1] + lNMom[1];
    Double_t lVecProdXi = (lV0px*lBMom[1] - lV0py*lBMom[0]) / (lV0Pt*lBachMod);
        
    if ( lMagneticField < 0 ) lVecProdXi *= -1; //invert sign
        
    fTreeCascVarIsCascadeCowboy = kFALSE;
    if (lVecProdXi < 0) fTreeCascVarIsCascadeCowboy = kTRUE;
        
    fTreeCascVarCascadeCowboyness = lVecProdXi;
    //========================================================================================
    */
        
    if ( lChargeXi < 0 )        lInvMassXiMinus     = xi->MassXi();
    if ( lChargeXi > 0 )        lInvMassXiPlus         = xi->MassXi();
    if ( lChargeXi < 0 )        lInvMassOmegaMinus     = xi->MassOmega();
    if ( lChargeXi > 0 )        lInvMassOmegaPlus     = xi->MassOmega();

		

    Double_t lXiMom[3];
	
    //	xi->GetPxPyPz(lXiMom);
    lXiMom[0]=xi->Px();
    lXiMom[1]=xi->Py();
    lXiMom[2]=xi->Pz();
    lXiMomX=lXiMom[0];
    lXiMomY=lXiMom[1];
    lXiMomZ=lXiMom[2];
    lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );

    lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );

    lRapXi    = xi->RapXi();
    lRapOmega = xi->RapOmega();

    if (	TMath::Sqrt( pow( lXiMomX,2)+ pow(lXiMomY,2) + pow(lXiMomZ,2)) > 1.e-5) {
      kctauXi=NominalMassXi* lXiDecayLength/TMath::Sqrt( pow( lXiMomX,2)+ pow(lXiMomY,2) + pow(lXiMomZ,2));
    }
    else {
      kctauXi=-999.;
    }
    if (kctauXi > 100) continue;
    fHistEventV0->Fill(10);  
    if(fReadMCTruth){
      if (fMCEvent){

	if (isXiPos)    fHistEventXiTruePos->Fill(10,PtMotherBach);  
	if (isXiNeg)    fHistEventXiTrueNeg->Fill(10,PtMotherBach);  
      }    
    }
    //	cout << "Pt Cascade " << lXiTransvMom << endl;
    //	cout << "P cascade component x " << lXiMom[0]<< endl;
	
    //----------------------------------------
    // Calculate Cascade DCA to PV, please
    //----------------------------------------
        
    Int_t lChargeCascade = 0;
        
    //cascade properties to get started
    Double_t xyzCascade[3], pxpypzCascade[3], cvCascade[21];
    for(Int_t ii=0;ii<21;ii++) cvCascade[ii]=0.0; //something small
        
    xyzCascade[0] = xi->DecayVertexXiX();
    xyzCascade[1] = xi->DecayVertexXiY();
    xyzCascade[2] = xi->DecayVertexXiZ();

    AliExternalTrackParam lCascTrajObject(xyzCascade,lXiMom,cvCascade,lChargeCascade), *hCascTraj = &lCascTrajObject;
        
    Double_t lCascDCAtoPVxy = TMath::Abs(hCascTraj->GetD(lBestPrimaryVtxPos[0],
							 lBestPrimaryVtxPos[1],
							 lMagneticField) );
    Float_t dzcascade[2];
    hCascTraj->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dzcascade );
    Double_t lCascDCAtoPVz = dzcascade[1];


    //----------------------------------------------------

    cout << "\n\nfReadMCTruth" << fReadMCTruth << endl;
    cout << "pdg of daughter pos " << PdgPos<< endl;
    cout << "pdg of daughter neg " << PdgNeg<< endl;
    cout << "pdg of daughter bach " << PdgBach<< endl;
    cout << "pdg of gmother pos (Xi candidate) " << PdgGMotherPos<< endl;
    cout << "V0 lifetime " <<             lV0Lifetime<< endl;
    cout << "Xi ctau " <<             kctauXi<< endl;
    cout << "Pt Cascade " << lXiTransvMom << endl;
    cout << "P cascade component x " << lXiMom[0]<< endl;
    cout << "P cascade component y " << lXiMom[1]<< endl;
    cout << "P cascade component z " << lXiMom[2]<< endl;
    cout << "DCA cascade to PV xy " <<         lCascDCAtoPVxy;
    cout << " DCA cascade to PV z " <<         lCascDCAtoPVz;
    cout << " DCA cascade to PV " <<         lDcaXiToPrimVertex << endl;
    //------------------------------------------------
    // Set Variables for adding to tree
    //------------------------------------------------

    //Copy Multiplicity information
    // fTreeVariableVarRefMultEta8 = fRefMultEta8;
    // fTreeVariableVarRefMultEta5 = fRefMultEta5;
    fTreeVariablePosTrackLength= lPosTrackLength;
    fTreeVariableNegTrackLength= lNegTrackLength;
    fTreeVariableBachTrackLength= lBachTrackLength;
    fTreeVariableRunNumber = fRunNumber;
    fTreeVariableBunchCrossNumber = fBunchCrossNumber;
    fTreeVariableMultiplicity	      = lPercentiles;
    fTreeVariableZvertex                = lBestPrimaryVtxPos[2];

    fTreeVariableChargeCasc	= lChargeXi;
    fTreeVariableEtaCasc	= lEtaXi;
    fTreeVariablePhiCasc	= lPhiXi;
    fTreeVariableThetaCasc	= lThetaXi;
    fTreeVariableRapXi = lRapXi ;
    fTreeVariableRapOmega = lRapOmega ;
    fTreeVariableDcaXiToPrimVertex      = 	  lDcaXiToPrimVertex;
    fTreeVariableXYDcaXiToPrimVertex      = 	  lCascDCAtoPVxy;
    fTreeVariableZDcaXiToPrimVertex      = 	  lCascDCAtoPVz;
    fTreeVariableDcaV0ToPrimVertex      = 	  lDcaV0ToPrimVertexXi;
    fTreeVariableDcaPosToPrimVertex     = 	    lDcaPosToPrimVertexXi;    	      
    fTreeVariableDcaNegToPrimVertex     = 	       	       lDcaNegToPrimVertexXi;
    fTreeVariableDcaV0Daughters = lDcaV0DaughtersXi;
    fTreeVariableDcaCascDaughters = lDcaXiDaughters;
    fTreeVariableDcaBachToPrimVertex= lDcaBachToPrimVertexXi;

    fTreeVariableV0CosineOfPointingAngleSpecial=  lV0CosineOfPointingAngleXiSpecial; 
    fTreeVariableV0CosineOfPointingAngle=  lV0CosineOfPointingAngleXi;
    fTreeVariableCascCosineOfPointingAngle = lXiCosineOfPointingAngle;
    fTreeVariablePtCasc		      =lXiTransvMom;
    fTreeVariablectau		      = kctauXi;
    fTreeVariableV0Lifetime = lV0Lifetime;
    if(lInvMassXiMinus!=0)    fTreeVariableInvMassXi = lInvMassXiMinus;
    if(lInvMassXiPlus!=0)     fTreeVariableInvMassXi = lInvMassXiPlus;
    if(lInvMassOmegaMinus!=0) fTreeVariableInvMassOmega = lInvMassOmegaMinus;
    if(lInvMassOmegaPlus!=0)  fTreeVariableInvMassOmega = lInvMassOmegaPlus;
    fTreeVariableInvMassLambda = lInvMassLambdaAsCascDghter;
    fTreeVariableInvMassK0Short = lInvMassK0sAsCascDghter;

    fTreeVariableCascRadius = lXiRadius;
    fTreeVariableV0Radius = lV0RadiusXi;
    fTreeVariableLeastNbrClusters = leastnumberofclusters;

    if( !( (fTreeVariableInvMassXi<1.32+0.075&&fTreeVariableInvMassXi>1.32-0.075) ||(fTreeVariableInvMassOmega<1.68+0.075&&fTreeVariableInvMassOmega>1.68-0.075) ) ) continue;
    //        for(Int_t i=0; i<20; i++) fTreeVariableVarRefMultDiffEta[i] = fRefMultDiffEta[i];

    //        fTreeVariableVarDistOverTotMom = TMath::Sqrt(
    //                                         TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
    //                               TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
    //                               TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
    //                           );
    //        fTreeVariableVarDistOverTotMom /= (lXiTotMom+1e-13);

    //All vars not specified here: specified elsewhere!

    //------------------------------------------------
    // Fill Tree!
    //------------------------------------------------

    // The conditional is meant to decrease excessive
    // memory usage! Be careful when loosening the
    // cut!

    //Xi    Mass window: 150MeV wide
    //Omega mass window: 150MeV wide


    fSignalTree->Fill();

    if(fReadMCTruth){
      if (fMCEvent){
	if (isXiPos)  	fHistEventXiTruePos->Fill(11,PtMotherBach);  
	if (isXiNeg) 	fHistEventXiTrueNeg->Fill(11,PtMotherBach);  
      }
    }
    NumberSecondParticle++;
    if (isXiNeg)       NumberSecondParticleRecoTrueXiNeg++;
    if (isXiPos)       NumberSecondParticleRecoTrueXiPos++;


    if(fReadMCTruth){
      if (fMCEvent){
	if(isXiPos || isXiNeg){
	  if (isXiPos) isXi=0.5;
	  else if (isXiNeg) isXi=-0.5;
	  fHistSelectedXiPt[0]->Fill(lXiTransvMom, lPercentiles,isXi );
	  if (TMath::Abs(lRapXi)<=0.5 ) {
	    //	cout << "these are all K0s generated passing selection criteria: label K0s " << particle->Label()<<endl; 
	    fHistSelectedXiPt[1]->Fill(lXiTransvMom, lPercentiles, isXi);
	  }
	  fHistResolutionXiPt->Fill(lXiTransvMom-MotherBach->Pt(),lXiTransvMom);
	  fHistResolutionXiPhi->Fill(xi->Phi()-MotherBach->Phi(),lXiTransvMom);
	  fHistResolutionXiEta->Fill(lEtaXi-MotherBach->Eta(), lXiTransvMom);
	  NumberTruthXi++;
	}
	else if(isOmegaPos || isOmegaNeg){
	  if (isOmegaPos) isXi=0.5;
	  else if (isOmegaNeg) isXi=-0.5;
	  fHistSelectedOmegaPt[0]->Fill(lXiTransvMom, lPercentiles,isXi);
	  if (TMath::Abs(lRapOmega)<=0.5 ) {
	    //	cout << "these are all K0s generated passing selection criteria: label K0s " << particle->Label()<<endl; 
	    fHistSelectedOmegaPt[1]->Fill(lXiTransvMom, lPercentiles, isXi);
	  }
	  fHistResolutionOmegaPt->Fill(lXiTransvMom-MotherBach->Pt(),lXiTransvMom);
	  fHistResolutionOmegaPhi->Fill(xi->Phi()-MotherBach->Phi(), lXiTransvMom);
	  fHistResolutionOmegaEta->Fill(lEtaXi-MotherBach->Eta(),lXiTransvMom);
	  NumberTruthOmega++;
	}

      }
    }
    //------------------------------------------------
    // Fill tree over.
    //------------------------------------------------
    fMassXiMinus->Fill(lInvMassXiMinus);    
    fMassXiPlus->Fill(lInvMassXiPlus);    
    fV0Lifetime->Fill(lV0Lifetime);    
    fHistPDG-> Fill(PdgGMotherPos);
    fHistPDGLambda-> Fill(PdgMotherPos);
    fHistPDGBachMom-> Fill(PdgMotherBach);
    fHistTheta-> Fill(lThetaXi);
    fHistEta-> Fill(lEtaXi);
    fHistPhi-> Fill(lPhiXi);
  }// end of the Cascade loop (AOD)


  if((fReadMCTruth && isEfficiency) || (!fReadMCTruth)) NumberSecondParticleAll=NumberSecondParticle;
  else if (fReadMCTruth && !isEfficiency) NumberSecondParticleAll=NumberSecondParticleMC;

  //Fill histos about selected events 
  fHistMultvsV0->Fill(NumberSecondParticle,lPercentiles);
  fHist_multiplicity->Fill(lPercentiles);
  fHistZvertex->Fill(lBestPrimaryVtxPos[2]);
  fHistMultiplicityVsVertexZ->Fill(lBestPrimaryVtxPos[2], lPercentiles);


  if(NumberSecondParticleAll==0){   
    PostData(1, fOutputList);
    PostData(2, fSignalTree );
    cout << "event has trigger particle but no cascade " << endl;
    PostData(4, fOutputList2);  
    PostData(5, fOutputList3);     
    return;
  }



  cout << "*******************************************************************************************************************************" << endl;
 
  fHistEventV0->AddBinContent(15, NumberSecondParticle);    
  fHistEventV0->AddBinContent(16, NumberTruthXi);    
  fHistEventV0->AddBinContent(17, NumberTruthOmega);    
  fHistEventMult->Fill(22);  

 
  PostData(1, fOutputList);     
  PostData(2,fSignalTree);
  PostData(4, fOutputList2);  
  PostData(5, fOutputList3);     

}

//_____________________________________________________________________________
void AliAnalysisTaskCascades::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}

//________________________________________________________________________
Float_t AliAnalysisTaskCascades::GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b ){
  // Input parameters:
  //   deltaY - user defined "dead region" in cm
  //   deltaZ - user defined "active region" in cm (250 cm drift lenght - 14 cm L1 delay
  //   b     - magnetic field 
  AliESDtrack esdTrack( gt );
  esdTrack.SetESDEvent((AliESDEvent*) gt->GetEvent() );
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(gt);
  esdTrack.ResetTrackParamIp(&etp);
  return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}
//_____________________________________________________________________________
