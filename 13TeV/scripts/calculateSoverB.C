#include <stdio.h>
#include <iostream>
#include <TMath.h>

using namespace std;
using namespace ROOT;

void calculateSoverB()
{
 TString prefix = "ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID";

 TFile* fBkgMC[3];
 TH1D * h_tempMC[3];

 fBkgMC[0] = new TFile("rootFiles/GJ_"+prefix+".root","READ");
 h_tempMC[0] = (TH1D*)fBkgMC[0]->Get("h_mass_VarBin_MassCut");
 TH1D * h_Bkg = (TH1D*)h_tempMC[0]->Clone("h_Bkg");

 fBkgMC[1] = new TFile("rootFiles/DiJet_"+prefix+".root","READ");
 h_tempMC[1] = (TH1D*)fBkgMC[1]->Get("h_mass_VarBin_MassCut");
 h_Bkg->Add(h_tempMC[1]);
 
 fBkgMC[2] = new TFile("rootFiles/EWK_"+prefix+".root","READ");
 h_tempMC[2] = (TH1D*)fBkgMC[2]->Get("h_mass_VarBin_MassCut");
 h_Bkg->Add(h_tempMC[2]);


 TFile * fSig = new TFile("rootFiles/QstarToGJ_M2000_f1p0_1.root", "READ");
 TH1D * h_Sig = (TH1D*)fSig->Get("h_mass_VarBin_MassCut");

 Int_t lowBin  = h_Sig->FindBin(1609 + 0.5);
 Int_t highBin = h_Sig->FindBin(2436 - 0.5); 

 std::cout<<" Low Bin = "<<h_Sig->GetBinLowEdge(lowBin)<<std::endl;
 std::cout<<" High Bin = "<<h_Sig->GetBinLowEdge(highBin) + h_Sig->GetBinWidth(highBin)<<std::endl<<std::endl;
 
 std::cout<< (h_Sig->Integral(lowBin, highBin)/h_Sig->Integral())*100.0 << " % of signal, S/sqrt(B): "
   <<(h_Sig->Integral(lowBin, highBin)/sqrt(h_Bkg->Integral(lowBin, highBin)))<<std::endl;

}
