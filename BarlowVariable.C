#include "Riostream.h"
#include "TH1F.h"

void BarlowVariable(TH1F* hVar, TH1F* hDefault, TH1F* hBarlowVar, TH1F* hSyst, Float_t NSigma, Int_t NSign, Bool_t IsBarlowSign){

  Float_t BarlowVar=0;
  Int_t NBarlowSign=0;

  //  if (hVar->GetNbinsX() != hDefault->GetNbinsX()) return;
  //  if (hVar->GetNbinsX() != hBarlowVar->GetNbinsX()) return;
  //  if (hVar->GetNbinsX() != hSyst->GetNbinsX()) return;

  for (Int_t b=1; b<=hVar->GetNbinsX();b++){
    if (hVar->GetBinContent(b) ==0 || hDefault->GetBinContent(b) ==0 || hVar->GetBinError(b) ==0 || hDefault->GetBinError(b) ==0 || (hVar->GetBinError(b) ==  hDefault->GetBinError(b)) ){
      hBarlowVar->SetBinContent(b,0);
      hBarlowVar->SetBinError(b,0);      
      continue;
    }
    BarlowVar = (hVar->GetBinContent(b) - hDefault->GetBinContent(b)) / sqrt(TMath::Abs(pow(hVar->GetBinError(b), 2)  - pow( hDefault->GetBinError(b),2) ));
    hBarlowVar->SetBinContent(b,BarlowVar);
    hBarlowVar->SetBinError(b,0);
    //    cout << hVar->GetBinContent(b) << " " << hDefault->GetBinContent(b) << " " << hVar->GetBinError(b) << " " << hDefault->GetBinError(b) << endl;
    //    cout << "b " << b << " " << BarlowVar << endl;
    if (TMath::Abs(BarlowVar) > NSigma) NBarlowSign++;
  }
  if (NBarlowSign >= NSign) IsBarlowSign=kTRUE;
  if (IsBarlowSign)  hBarlowVar->SetMarkerStyle(20);
  else hBarlowVar->SetMarkerStyle(33);

  //---set syst uncertainty--------
  for (Int_t b=1; b<=hVar->GetNbinsX();b++){
    if (hVar->GetBinContent(b) ==0 || hDefault->GetBinContent(b) ==0 ){
      hSyst->SetBinContent(b,0);
      hSyst->SetBinError(b,0);
      continue;
    }
    hSyst-> SetBinContent(b, (TMath::Abs(hVar->GetBinContent(b) - hDefault->GetBinContent(b))/2) / hDefault->GetBinContent(b));
    hSyst->SetBinError(b,0);
  }
}
