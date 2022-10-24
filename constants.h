Float_t massK0s = 0.497614; //GeV/c^2
Float_t ctauK0s= 2.6844; //cm 
Float_t massXi = 1.32171; //GeV/c^2
Float_t ctauXi = 4.91; //cm

//MCpredictions
const Int_t nummolt = 10;
TString SmoltGenOnTheFly[nummolt+1]={"0-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-300", "0-300"};
TString SBismoltGenOnTheFly[nummolt+1]={"0-30", "30-45", "45-54", "54-63", "63-72", "72-81", "81-90", "90-105", "105-120", "120-300", "0-300"};
Double_t NmoltGenOnTheFly[nummolt+1]={0, 30, 45, 54, 63, 72, 81, 90, 105, 120, 300};

//dNdeta final values computed using BK task (for DATA)

//DATA 13 TeV
Float_t dNdEtaFinal13TeV[nummolt] = {7.396, 10.703, 14.799, 18.930, 24.039, 32.047, 34.063, 37.6, 0, 0}; 
Float_t dNdEtaFinal13TeV_ErrorL[nummolt] = {0.143, 0.184, 0.21, 0.262, 0.221, 0.57,  0.62, 0.70, 0, 0}; 
Float_t dNdEtaFinal13TeV_ErrorR[nummolt] = {0.165, 0.175, 0.18, 0.229, 0.274, 0.48, 0.54, 0.65, 0, 0}; 

//DATA 13 TeV
Float_t dNdEtaFinal5TeV[nummolt] = {10.234, 16.955, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t dNdEtaFinal5TeV_ErrorL[nummolt] = {0.172, 0.292, 0, 0, 0, 0, 0, 0, 0, 0};
Float_t dNdEtaFinal5TeV_ErrorR[nummolt] = {0.139, 0.198, 0, 0, 0, 0, 0, 0, 0, 0};
