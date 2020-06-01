/*Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCascades_H
#define AliAnalysisTaskCascades_H
class AliPIDResponse;
class AliMultSelection;
class AliAODMCParticle;
class AliCentrality;
#include "AliAnalysisTaskSE.h"
#include "AliAODcascade.h"
#include "AliEventCuts.h"

class AliAnalysisTaskCascades : public AliAnalysisTaskSE  
{
 public:
  AliAnalysisTaskCascades();
  AliAnalysisTaskCascades(const char *name);
  virtual                 ~AliAnalysisTaskCascades();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);
  
  void SetMinPt(Float_t ptmin) {fminPtj = ptmin;}
  void SetMaxPt(Float_t ptmax) {fmaxPtj = ptmax;}
  void SetCorr(Bool_t ishhCorr){fIshhCorr = ishhCorr;}
  void SetMC(Bool_t isMC){fReadMCTruth = isMC;}
  void SetEff(Bool_t isEff){isEfficiency = isEff;}
  void SetEvtToMix(Int_t EvtToMix){fnEventsToMix = EvtToMix;}
  void SetEtaTrigger(Float_t EtaTrigger){fEtaTrigger = EtaTrigger;}
  void SetEtahAssoc(Float_t EtahAssoc){fEtahAssoc = EtahAssoc;}

  void ProcessMCParticles(Bool_t Generated, Float_t lPercentiles, Bool_t isV0, Bool_t fIshhCorr);

  Float_t GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b );

 private:
  TString                 fAnalysisType;                  // "ESD" or "AOD" analysis type
  TString                 fCollidingSystem;               // "pp", "pPb", "PbPb" 
  AliAODEvent*            fAOD;             //! input event
  AliPIDResponse *        fPIDResponse;     //!PID response object 
  //  AliMultSelection *      fMultSelection;   //! 
  AliEventCuts            fEventCuts; //! 
  
  TList*                  fOutputList;      //! output list
  TTree*                  fSignalTree;      //! output tree
  TList*                  fOutputList2;      //! output list
  TList*                  fOutputList3;      //! output list
  
  AliMCEvent *            fMCEvent;         //!
  Bool_t                  fReadMCTruth;
  Bool_t                  fIshhCorr;
  Bool_t                  isEfficiency;

  Int_t                    fzVertexBins; 
  Int_t                    fnMultBins;	 
  Int_t                    fMaxFirstMult;
  Int_t                    fMaxSecondMult;
  Int_t                    fnEventsToMix;
  Float_t                    fEtaTrigger;
  Float_t                    fEtahAssoc;

   TH1F*                   fHistPt;          //! 
   TH1F*                   fHistPtTriggerParticle;          //! 
   TH1F*                   fHistDCAxym1;          //! 
   TH1F*                   fHistDCAxym2;          //! 
   TH1F*                   fHistDCAzm1;          //! 
   TH1F*                   fHistDCAzm2;          //! 
   TH1F*                   fHistPtV0;          //! 
   TH1F*                   fHistPthAssoc;          //! 
   TH2F*                   fHistPtTMaxBefAllCfrDataMC; //!
   TH1F*                   fHistPtTMinBefAll;          //! 
   TH1F*                   fHistPtTMinBefAllMC;          //! 
   TH1F*                   fHistPtTMaxBefAll;          //! 
   TH1F*                   fHistPtTMaxBefAllBis;          //! 
   TH1F*                   fHistPtTMaxBefAllMC;          //! 
   TH2F*                   fHistPtvsMult;          //! 
   TH2F*                   fHistPtvsMultBefAll;          //! 
   TH2F*                   fHistPtMaxvsMult;          //! 
   TH2F*                   fHistPtMaxvsMultBefAll;          //! 
  TH1F *                  fHistZvertex;     //!
  TH3F* fHistNumberChargedAllEvents; //!
  TH3F* fHistNumberChargedNoTrigger; //!
  TH3F* fHistNumberChargedTrigger; //!
  TH2F *                  fHist_eta_phi;   //!
  TH2F *                  fHist_eta_phi_PtMax;   //!
  TH1F*                   fHist_multiplicity; //!
  TH1F*                   fHist_multiplicity_EvwTrigger; //!
  TH1F*                   fHistEventMult;   //!
Double_t   fRunNumber; //!
Double_t   fBunchCrossNumber; //!
  TH1F*                   fHistEventV0;   //!
  TH3F*                   fHistEventXiTrue;   //!
  TH3F*                   fHistEventOmegaTrue;   //!

TH2F *fHistDCApTrackXi;//!
TH2F *fHistDCAnTrackXi;//!
TH2F *fHistDCAbachTrackXi; //!
TH2F*  fHistLengthvsCrossedRowsPos; //!
TH2F*  fHistLengthvsCrossedRowsNeg; //!
TH2F* fHistLengthvsCrossedRowsBach; //!
TH2F*  fHistLengthvsCrossedRowsAfterSelPos; //!
TH2F*  fHistLengthvsCrossedRowsAfterSelNeg; //!
TH2F* fHistLengthvsCrossedRowsAfterSelBach; //!

 TH2F*   fHistCfrDiffDefXiPt; //!
 TH2F*   fHistCfrDiffDefXiP; //!

  TH1F*                   fHistTrack;       //!
  TH2F* fHistTriggerComposition; //! 
  TH2F* fHistTriggerCompositionMCTruth; //! 
  TH2F* fHistAssocComposition; //! 
  TH2F* fHistAssocCompositionMCTruth; //! 
  TH1F*                   fHistTrackAssoc;       //!
  TH1F*                   fHistPDG;         //!	
  TH1F*                   fHistPDGLambda;         //!	
  TH1F*                   fHistPDGBachMom;         //!	
  TH1F*                   fHistTheta; //!
  TH1F*                   fHistEta; //!
TH1F*  fHistRapGenXi; //!
TH1F*  fHistRapGenOmega; //!
TH1F*  fHistRapSelXi; //!
TH1F*  fHistRapSelOmega; //!

  TH1F*                   fHistPhi; //!
  TH1F*                   fHistTrackBufferOverflow;         //!	
  TH2F*                   fHistSecondParticleAll; //!
  TH2F*                   fHistSecondParticleTruthAll; //!
  TH2F*                   fHistSecondParticle; //!
  TH2F*                   fHistSecondParticleTruth; //!
  TH1F*                   fMassXiPlus;          //!
  TH1F*                   fMassXiMinus;          //!
  TH1F*                   fV0Lifetime;          //!
TH1F *  fV0DistanceTrav; //! 
TH1F *  fV0TotMomentum; //! 
  TH2F *                  fHistMultvsV0All; //!
  TH2F *                  fHistMultvsV0AllTruth; //!
  TH2F *                  fHistMultvsV0MCAll; //!
  TH2F *                  fHistMultvsV0; //!
  TH2F *                  fHistMultvsV0Truth; //!
  TH2F *                  fHistMultvsV0MC; //!
  TH3F*  fHistTriggerNotLeading; //!
    TH3F*  fHistTriggerNotLeadingMC; //!
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
  TH3F *  fHistGeneratedXiPt; //!
  TH3F *  fHistSelectedXiPt; //!
  TH3F *  fHistGeneratedOmegaPt; //!
  TH3F *  fHistSelectedOmegaPt; //!
  TH3F *  fHistReconstructedV0PtMass; //!
  TH3F *  fHistSelectedV0PtMass; //!  

TH2F*  fHistTriggerPtRecovsPtGenCasc; //!
TH2F*  fHistTriggerPtRecovsPtGenPos; //!
TH2F*  fHistTriggerPtRecovsPtGenNeg; //!
TH2F*  fHistTriggerPtRecovsPtGenBach; //!
TH2F*  fHistTriggerYRecovsYGenCasc; //!

TH3F*  fHistResolutionTriggerPt; //!
TH3F*  fHistResolutionTriggerPhi; //!
TH3F*  fHistResolutionTriggerEta; //!
TH2F*  fHistResolutionXiPt; //!
TH2F*  fHistResolutionXiPhi; //!
TH2F*  fHistResolutionXiEta; //!
TH2F *  fHistResolutionWithCrossedRowsOverFindXiPt; //!
TH2F *  fHistResolutionWithCrossedRowsOverFindXiPhi; //!
TH2F *  fHistResolutionWithCrossedRowsOverFindXiEta; //!

TH2F*  fHistResolutionOmegaPt; //!
TH2F*  fHistResolutionOmegaPhi; //!
TH2F*  fHistResolutionOmegaEta; //!

TH3F*  fHistResolutionTriggerPhiPt; //!
TH3F*  fHistResolutionTriggerPhiPdgCode; //!

TH1F*  fHistCrossedRowsOverFindablePos; //!
TH1F*  fHistCrossedRowsOverFindableNeg; //!
TH1F*  fHistCrossedRowsOverFindableBach; //!


  TH2F *** fHistPrimaryTrigger; //!
  TH3F ***  fHistPrimaryV0; //!  
 
  Float_t                 fminPtj;
  Float_t                 fmaxPtj;
  TString                 fV0;
  Float_t                 fminPtV0;
  Float_t                 fmaxPtV0;
  Float_t                 fminPthAssoc;
  Float_t                 fmaxPthAssoc;
  Int_t                   Evcounter;
  Int_t                   Evcounterczero;
  Int_t                   fmolt;
  
  Int_t *                 farrGT;                //!
  UShort_t                fTrackBufferSize;      // Size fo the above array, ~12000 for PbPb

  //tree leaf
Double_t  fTreeVariableMultiplicity;       
Double_t  fTreeVariableZvertex;             
Double_t  fTreeVariablePDGCode;             
Double_t  fTreeVariablePDGCodeBach;             
Double_t  fTreeVariablePDGCodeNeg;             
Double_t  fTreeVariablePDGCodePos;             
Double_t  fTreeVariablePDGCodeLambda;             
Double_t  fTreeVariablePDGCodeMotherLambda;             
Double_t  fTreeVariableRunNumber;           
Double_t  fTreeVariableBunchCrossNumber;    
Double_t  fTreeVariableNegNSigmaPion;
Double_t  fTreeVariableNegNSigmaProton;
Double_t  fTreeVariablePosNSigmaPion;
Double_t  fTreeVariablePosNSigmaProton;
Double_t  fTreeVariableBachNSigmaPion;
Double_t  fTreeVariableBachNSigmaKaon;
 Double_t	fTreeVariablePosTrackLength;
 Double_t	fTreeVariableNegTrackLength;
 Double_t	fTreeVariableBachTrackLength;
Double_t  fTreeVariableDcaXiToPrimVertex;
Double_t  fTreeVariableXYDcaXiToPrimVertex;
Double_t  fTreeVariableZDcaXiToPrimVertex;
Double_t  fTreeVariableDcaV0ToPrimVertex;
Double_t  fTreeVariableDcaPosToPrimVertex; 
Double_t  fTreeVariableDcaNegToPrimVertex; 
Double_t  fTreeVariableDcaV0Daughters; 
Double_t  fTreeVariableDcaCascDaughters; 
Double_t  fTreeVariableDcaBachToPrimVertex; 
Double_t  fTreeVariableV0CosineOfPointingAngle;
Double_t  fTreeVariableV0CosineOfPointingAngleSpecial;
Double_t  fTreeVariableCascCosineOfPointingAngle;
Double_t  fTreeVariablePtCasc;              
Double_t  fTreeVariableChargeCasc;          
Double_t  fTreeVariableEtaCasc;          
Double_t  fTreeVariablePhiCasc;          
Double_t  fTreeVariableThetaCasc;          
Double_t  fTreeVariablectau;               
Double_t  fTreeVariableInvMassXi;        
Double_t  fTreeVariableInvMassOmega;      
Double_t  fTreeVariableInvMassLambda;  
Double_t  fTreeVariableInvMassK0Short;  
Double_t  fTreeVariableRapXi;               
Double_t  fTreeVariableRapOmega;            
Double_t  fTreeVariableCascRadius;     
Double_t  fTreeVariableV0Radius;     
Double_t  fTreeVariableLeastNbrClusters;    
Double_t  fTreeVariableV0Lifetime;    
 Double_t fTreeVariableIsPrimaryXi;
 Double_t fTreeVariableIsPrimaryOmega;

  bool FifoShiftok;	                        		      

  /* Double_t fTreeVariableDcaV0Daughters;			       */
  /* Double_t fTreeVariableV0Radius;			       */
  /* Double_t fTreeVariableLeastNbrCrossedRowsPos;		       */
  /* Double_t fTreeVariableLeastRatioCrossedRowsOverFindablePos;     */
  /* Double_t fTreeVariableLeastNbrCrossedRowsNeg;		       */
  /* Double_t fTreeVariableLeastRatioCrossedRowsOverFindableNeg;      */
  /* Double_t fTreeVariableMultiplicity;       */
  

  AliAnalysisTaskCascades(const AliAnalysisTaskCascades&); // not implemented
  AliAnalysisTaskCascades& operator=(const AliAnalysisTaskCascades&); // not implemented

  ClassDef(AliAnalysisTaskCascades, 1);
};

#endif
