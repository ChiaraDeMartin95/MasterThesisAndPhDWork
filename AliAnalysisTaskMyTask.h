/*Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMyTask_H
#define AliAnalysisTaskMyTask_H
class AliPIDResponse;
class AliMultSelection;
class AliAODMCParticle;
class AliCentrality;
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisKPEventCollectionChiara.h"
#include "AliEventCuts.h"

class AliAnalysisTaskMyTask : public AliAnalysisTaskSE  
{
 public:
  AliAnalysisTaskMyTask();
  AliAnalysisTaskMyTask(const char *name);
  virtual                 ~AliAnalysisTaskMyTask();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);
  
  void SetMinPt(Float_t ptmin) {fminPtj = ptmin;}
  void SetMaxPt(Float_t ptmax) {fmaxPtj = ptmax;}
  void SetMC(Bool_t isMC){fReadMCTruth = isMC;}
  void SetEff(Bool_t isEff){isEfficiency = isEff;}
  void SetEvtToMix(Int_t EvtToMix){fnEventsToMix = EvtToMix;}

  void ProcessMCParticles(Bool_t Generated, AliAODTrack* track, Int_t& labelPrimOrSec, Float_t lPercentiles, Bool_t isV0);
  // double CalculateMass(double momentum1[3], double momentum2[3], double mass1, double mass2); */
  Double_t CalculateDeltaTheta( Double_t theta1, Double_t theta2 ); 
  Double_t CalculateDeltaPhi( Double_t phi1, Double_t phi2 ) ; 
  Double_t CalculateDeltaEta( Double_t eta1, Double_t eta2 ) ; 
  double CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad);

  //  Double_t ThetaS( Double_t posSftR125[3] )const; 
  //Double_t EtaS( Double_t posSftR125[3] ) const ; 

   void DoPairsh1h2 ( const Float_t lcentrality, int fieldsignf, Double_t lBestPrimaryVtxPos);
  //void DoPairshh ( const Float_t lcentrality, int fieldsignf);

 private:
  TString                 fAnalysisType;                  // "ESD" or "AOD" analysis type
  TString                 fCollidingSystem;               // "pp", "pPb", "PbPb" 
  AliAODEvent*            fAOD;             //! input event
  AliPIDResponse *        fPIDResponse;     //!PID response object 
  //  AliMultSelection *      fMultSelection;   //!
  //  AliEventCuts            fEventCuts; //!
  
  TList*                  fOutputList;      //! output list
  TTree*                  fSignalTree;      //! output tree
  TTree*                  fBkgTree;      //! output tree
  
  AliMCEvent *            fMCEvent;         //!
  Bool_t                  fReadMCTruth;
  Bool_t                  isEfficiency;
  AliAnalysisKPEventCollectionChiara ***fEventColl;  //!
  AliAnalysisKPEventChiara *    fEvt;                //!

  Int_t                    fzVertexBins; 
  Int_t                    fnMultBins;	 
  Int_t                    fMaxFirstMult;
  Int_t                    fMaxSecondMult;
  Int_t                    fnEventsToMix;

   TH1F*                   fHistPt;          //! 
   TH1F*                   fHistPtV0;          //! 
   TH1F*                   fHistPtTMin;          //! 
   TH1F*                   fHistPtTMinMC;          //! 
   TH2F*                   fHistPtvsMult;          //! 
   TH2F*                   fHistPtvsMultBefAll;          //! 
  TH1F *                  fHistZvertex;     //!
  TH2F *                  fHist_eta_phi ;   //!
  TH1F*                   fHist_multiplicity; //!
  TH1F*                   fHistEventMult;   //!
  TH1F*                   fHistEventV0;   //!
  TH1F*                   fHistTrack;       //!
  TH1F*                   fHistPDG;         //!	
  TH2F*                   fHistSecondParticleAll; //!
  TH2F*                   fHistSecondParticleTruthAll; //!
  TH2F*                   fHistSecondParticle; //!
  TH2F*                   fHistSecondParticleTruth; //!
  TH1F*                   fMassV0;          //!
  TH2F *                  fHistMultvsV0All; //!
  TH2F *                  fHistMultvsV0AllTruth; //!
  TH2F *                  fHistMultvsV0MCAll; //!
  TH2F *                  fHistMultvsV0; //!
  TH2F *                  fHistMultvsV0Truth; //!
  TH2F *                  fHistMultvsV0MC; //!
  TH2F**                   fHistMassvsPt;                               //!
  TH2F**                   fHistMassvsPt_tagli;                         //!
  TH2F*                  fHistMultvsTriggerBefAll; //!
  TH2F*                  fHistMultvsTriggerMCTruthBefAll; //!
  TH2F*                  fHistMultvsTriggerAll; //!
  TH2F*                  fHistMultvsTriggerMCTruthAll; //!
  TH2F*                  fHistMultvsTrigger; //!
  TH2F*                  fHistMultvsTriggerMCTruth; //!
  TH1F*                  fHistMassPhoton;  //!
  TH1F*                  fHistMass2Photon;  //!
  TH2F*                  fHistPtArmvsAlpha;  //!
  TH2F*                  fHistPtArmvsAlphaAfterSelection;  //!
  TH2F*                  fHistPtArmvsAlphaAfterPhotonSelection;  //!
  TH2F*                  fHistPtArmvsAlphaAfterLambdaRejectionSelection;  //!
  TH1F*                  fHistTrigger;//!
  TH1F*                  fHistTriggerMCTruth;//!
  TH1F*                  fHistTriggerwV0;//!
  TH1F*                  fHistTriggerwV0MCTruth;//!
  TH2F*                  fHistMultiplicityVsVertexZ; //!
  TH1F*                  fHistTriggervsMult; //!
  TH1F*                  fHistTriggervsMultMC; //!
  TH1F*                  fHistMultiplicityOfMixedEvent; //!
  TH3F *  fHistGeneratedTriggerPtPhi; //!
  TH3F **  fHistSelectedTriggerPtPhi; //!
  TH3F **  fHistGeneratedV0PtPhi; //!
  TH3F **  fHistSelectedV0PtPhi; //!
  TH3F *  fHistGeneratedTriggerPtEta; //!
  TH3F **  fHistSelectedTriggerPtEta; //!
   TH3F **  fHistGeneratedV0PtEta; //!
  TH3F **  fHistSelectedV0PtEta; //!
  TH3F *  fHistReconstructedV0PtMass; //!
  TH3F *  fHistSelectedV0PtMass; //!  

TH2F*  fHistResolutionTriggerPt; //!
TH2F*  fHistResolutionTriggerPhi; //!
TH2F*  fHistResolutionTriggerEta; //!
TH2F*  fHistResolutionV0Pt; //!
TH2F*  fHistResolutionV0Phi; //!
TH2F*  fHistResolutionV0Eta; //!

  TH2F *** fHistPrimaryTrigger; //!
  TH2F **  fHistPrimaryV0; //!  
 
  Float_t                 fminPtj;
  Float_t                 fmaxPtj;
  TString                 fV0;
  Float_t                 fminPtV0;
  Float_t                 fmaxPtV0;
  Int_t                   Evcounter;
  Int_t                   Evcounterczero;
  Int_t                   fmolt;
  
  Int_t *                 farrGT;                //!
  UShort_t                fTrackBufferSize;      // Size fo the above array, ~12000 for PbPb

  //tree leaf
  Double_t fTreeVariablePtTrigger;		       
  Double_t fTreeVariableChargeTrigger;		       
  Double_t fTreeVariableEtaTrigger; 		       
  Double_t fTreeVariablePhiTrigger;		       
  Double_t fTreeVariableDCAz;			       
  Double_t fTreeVariableDCAxy;			       
  Double_t fTreeVariableRapK0Short;		       	      
  Double_t fTreeVariableisPrimaryTrigger;
  Double_t fTreeVariableisPrimaryV0;
  Double_t fTreeVariableDcaV0ToPrimVertex ;	       	      
  Double_t fTreeVariableDcaPosToPrimVertex;	       	      
  Double_t fTreeVariableDcaNegToPrimVertex;	       	      
  Double_t fTreeVariableV0CosineOfPointingAngle;       	      
  Double_t fTreeVariablePtV0;			       
  Double_t fTreeVariablectau;			       
  Double_t fTreeVariableInvMassK0s;		       
  Double_t fTreeVariableInvMassLambda;		       
  Double_t fTreeVariableInvMassAntiLambda;		       
  Double_t fTreeVariableEtaV0;			       
  Double_t fTreeVariablePhiV0;			       
  Double_t fTreeVariablePtArmenteros;                   
  Double_t fTreeVariableAlpha;	   
  Double_t fTreeVariableDeltaEta;			       
  Double_t fTreeVariableDeltaPhi;			       
  Double_t fTreeVariableDeltaTheta;

  Double_t fTreeVariableMultiplicity;                   
  Double_t fTreeVariableZvertex;
  Double_t fTreeVariablePDGCode;

  bool FifoShiftok;	                        		      

  /* Double_t fTreeVariableDcaV0Daughters;			       */
  /* Double_t fTreeVariableV0Radius;			       */
  /* Double_t fTreeVariableLeastNbrCrossedRowsPos;		       */
  /* Double_t fTreeVariableLeastRatioCrossedRowsOverFindablePos;     */
  /* Double_t fTreeVariableLeastNbrCrossedRowsNeg;		       */
  /* Double_t fTreeVariableLeastRatioCrossedRowsOverFindableNeg;      */
  /* Double_t fTreeVariableMultiplicity;       */
  

  AliAnalysisTaskMyTask(const AliAnalysisTaskMyTask&); // not implemented
  AliAnalysisTaskMyTask& operator=(const AliAnalysisTaskMyTask&); // not implemented

  ClassDef(AliAnalysisTaskMyTask, 1);
};

#endif
