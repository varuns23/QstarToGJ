compare8TeVand13TeVQstarShape(){
  //TFile* file1 = new TFile("rootFiles/QstarToGJ_M1000_f1p0_1.root","READ");
  //TFile* file1 = new TFile("rootFiles/QstarToGJ_M2000_f1p0_1.root","READ");
  TFile* file1 = new TFile("rootFiles/QstarToGJ_M3000_f1p0_1.root","READ");

  //TFile* file2 = new TFile("rootFiles/8TeV/QstarToGJ_M_1000.root","READ");
  //TFile* file2 = new TFile("rootFiles/8TeV/QstarToGJ_M_2000.root","READ");
  TFile* file2 = new TFile("rootFiles/8TeV/QstarToGJ_M_3000.root","READ");

  TH1F *hist1  = (TH1F*)file1->Get("h_mass_VarBin_noMassCut");
  hist1->SetStats(0);
  hist1->SetLineColor(kRed);
  hist1->SetLineWidth(2);
  
  TH1F *hist2  = (TH1F*)file2->Get("h_mass_VarBin_finalnoPU");
  hist2->SetStats(0);
  hist2->SetLineColor(kBlue);
  hist2->SetLineWidth(2);


  TLegend *leg = new TLegend( 0.68, 0.7, 0.85, 0.85);
  leg->SetFillColor(0); leg->SetTextFont(42); leg->SetBorderSize(0);
  leg->AddEntry(hist1,"13 TeV","F");
  leg->AddEntry(hist2,"8 TeV","F");


  hist2->DrawNormalized("hist");
  hist1->DrawNormalized("histsame");
  leg->Draw();

}
