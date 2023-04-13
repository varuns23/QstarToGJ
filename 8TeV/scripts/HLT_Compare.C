HLT_Compare(){
  compare();

}

compare()
{
  
  //TFile* fileIn = new TFile("Data/Run2015D_PromptReco_v4_1.root","READ");
//  TFile* fileIn = new TFile("Data/Run2015D.root","READ");
  TFile* fileIn = new TFile("Data/Run2015D_PromptReco.root","READ");

  TH1F *h_Num = (TH1F*)fileIn->Get("h_ptTrigPhoton_num");
  TH1F *h_Den = (TH1F*)fileIn->Get("h_ptTrigPhoton_deno");
//  TH1F *h_Num = (TH1F*)fileIn->Get("h_GJmassTrig_num");
//  TH1F *h_Den = (TH1F*)fileIn->Get("h_GJmassTrig_deno");

  TGraphAsymmErrors* triggerTurnOn_PT = new TGraphAsymmErrors();

  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
  c1->SetFillColor(0);
  c1->SetGridy();
  c1->SetFrameBorderMode(0) ;

  //  TLegend *leg = new TLegend(0.82,0.86,0.99,0.99);
  //  leg->SetFillColor(0);
  //  leg->SetTextFont(40);
  //  leg->SetBorderSize(1);

  //  leg->AddEntry(SepEff_d,"Data","p");

  //////////////////////////////
  //  triggerTurnOn_PT->BayesDivide(h_125,h_90_125);
  triggerTurnOn_PT->Divide(h_Num,h_Den);
//  triggerTurnOn_PT->GetXaxis()->SetRangeUser(80,300);
  triggerTurnOn_PT->SetMarkerStyle(8);
  triggerTurnOn_PT->SetMarkerColor(2);
  triggerTurnOn_PT->SetMarkerSize(.6);
  triggerTurnOn_PT->Draw("AP");
  ////////////////////////////

  //  leg->Draw();
  //  c1->SaveAs("HLT_Efficiency.gif");

}
