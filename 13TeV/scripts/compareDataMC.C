void compareDataMC(){

  const Int_t nFiles = 4; 

  TString FileName[nFiles] = {
    "rootFiles/GJets_allHTBins.root",
    "rootFiles/QCD_allPtBins.root",
    "rootFiles/Run2015D_PromptReco.root",
    "rootFiles/QstarToGJ_M2000_f1p0_1.root"
  };

  const Float_t scaleMC   = 1.0;


  TFile* fileBkg1 = new TFile(FileName[0],"READ");
  TFile* fileBkg2 = new TFile(FileName[1],"READ");
  TFile* fileData = new TFile(FileName[2],"READ");
  TFile* fileSig = new TFile(FileName[3],"READ");

//  TString var[1]  = {"h_mass_VarBin_MassCut"};
  TString var[1]  = {"h_CutExpFlow"};

  TH1F *hSig  = (TH1F*)fileSig->Get(var[0]);   //Signal 2 TeV
  TH1F *hBkg1 = (TH1F*)fileBkg1->Get(var[0]);   //GJ MC
  TH1F *hBkg2 = (TH1F*)fileBkg2->Get(var[0]);   //DiJet MC
  TH1F *hData = (TH1F*)fileData->Get(var[0]);   //Data

  hSig->SetStats(0);
  hSig->SetLineColor(kMagenta);
  hSig->SetLineWidth(2);

  hBkg1->SetStats(0);
  hBkg1->SetLineColor(kOrange+9);
  hBkg1->SetFillColor(kOrange+6);
  hBkg1->SetLineWidth(2);

  hBkg2->SetStats(0);
  hBkg2->SetLineColor(kAzure+4);
  hBkg2->SetLineWidth(2);
  hBkg2->SetFillColor(kAzure+2);

  hData->SetStats(0); 
  hData->SetMarkerStyle(8); 
  hData->SetMarkerColor(kBlack); 
  hData->SetMarkerSize(0.8); 
  hData->SetLineColor(kBlack);


  const Int_t nleg = 4 ;
  TString legendName[nleg] = { "q* = 2 TeV", "GJ (Bkg)", "DiJet (Bkg)", "Data" };
  TString legendFill[nleg] = { "L", "F", "F", "Lp" };

  TLegend *leg = new TLegend( 0.68, 0.7, 0.85, 0.85);
  leg->SetFillColor(0); leg->SetTextFont(42); leg->SetBorderSize(0);
  leg->AddEntry(hSig,legendName[0],legendFill[0]);
  leg->AddEntry(hBkg1,legendName[1],legendFill[1]);
  leg->AddEntry(hBkg2,legendName[2],legendFill[2]);
  leg->AddEntry(hData,legendName[3],legendFill[3]);
	   
  TPaveText *cms = new TPaveText(0.72,0.9,0.83,0.95,"blNDC");
  cms->AddText("594 pb^{-1} (13 TeV)");
  cms->SetBorderSize(0);
  cms->SetFillColor(0);
  cms->SetFillStyle(0);
  cms->SetTextFont(42);
  cms->SetTextSize(0.035);
  TPaveText *lumi = new TPaveText(0.65,0.82,0.85,0.92,"blNDC");
  lumi->AddText("CMS Preliminary");
  lumi->SetBorderSize(0);
  lumi->SetFillColor(0);
  lumi->SetFillStyle(0);
  lumi->SetTextFont(42);
  lumi->SetTextSize(0.035);

  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
 
  THStack* hs = new THStack("hs","");
  hs->SetMinimum(0.01);
  hs->Add(hBkg2);
  hs->Add(hBkg1); 
  hs->SetHistogram(hBkg1);

  TH1F *err_bkg = hBkg1->Clone();
  err_bkg->Add(hBkg2);
  err_bkg->SetFillColor(1); err_bkg->SetFillStyle(3004);

  hs->Draw("hist");
  hSig->Draw("histSAME");
  hData->Draw("SAMEP E1");
  err_bkg->Draw("SAME E2");
  leg->Draw();
  cms->Draw();
  lumi->Draw();
  c1->cd();

}
