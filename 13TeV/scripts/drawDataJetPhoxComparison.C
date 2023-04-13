{
  gStyle->SetOptTitle(0);

  Double_t xMin = 500.0;
  Double_t xMax = 4000.0;

  TFile * fileData = new TFile("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root","READ");
  TH1F * h_data = (TH1F*)fileData->Get("h_mass_VarBin_MassCut");

  TFile * fileJetPhox = new TFile("rootFiles/OutputFileCombineSumIso_Histo_JetPhox_15March2016_v2.root","READ");
  TH1F * h_temp  = (TH1F*)fileJetPhox->Get("hmass");


  TH1F * h_jetphox  = (TH1F*)h_data->Clone("h_jetphox");
  h_jetphox->Reset();

  for(int ibin = 1; ibin <= h_data->GetNbinsX(); ++ibin)
    h_jetphox->SetBinContent( h_jetphox->FindBin( h_temp->GetBinCenter(ibin) ), h_temp->GetBinContent(ibin));

  h_data->SetStats(0);
  h_data->SetLineColor(kBlack);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.8);
  h_data->GetXaxis()->SetRangeUser(xMin, xMax);
  h_data->GetYaxis()->SetRangeUser(5e-5, 1.);
  h_data->Scale(1.0/ h_data->Integral());

  h_jetphox->SetStats(0);
  h_jetphox->SetLineColor(kRed);
  h_jetphox->SetFillColor(kRed);
  h_jetphox->SetFillStyle(3002);
  h_jetphox->Scale(1.0/ h_jetphox->Integral());


  TLegend *leg = new TLegend(0.6,0.7,0.87,0.87);
  leg->SetFillColor(0);  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->AddEntry(h_data, "Data", "LP");
  leg->AddEntry(h_jetphox, "JETPHOX", "F");

  TH1F * h_ratio = (TH1F*)h_data->Clone("h_ratio");
  h_ratio->Add(h_jetphox, -1);
  h_ratio->Divide(h_data);
  h_ratio->SetMarkerStyle(4);
  h_ratio->SetMarkerColor(kRed+2);
  h_ratio->SetLineColor(kRed+2);
  h_ratio->GetXaxis()->SetNdivisions(508);
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetTitleSize(0.13);
  h_ratio->GetYaxis()->SetNdivisions(508);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetRangeUser(-1.,1.);

  TLine * ratioLine = new TLine (xMin, 0.0,  xMax, 0.0);
  ratioLine->SetLineStyle(2);
  ratioLine->SetLineColor(38);


  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  pad1->Draw(); pad1->cd();                  
  pad1->SetLogy();                         
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetTitle(0);                         
  //pad1->SetRightMargin(padRmargin);          
  //pad1->SetLeftMargin(padLmargin);           
  //pad1->SetTopMargin(padTmargin);            
  pad1->SetBottomMargin(0);     
  h_data->DrawNormalized("PE");
  h_jetphox->DrawNormalized("sameHIST");
  h_data->DrawNormalized("samePE");
  leg->Draw();
  c1->cd();
  ///////////////////////////////////////////           
  // Second PAD for making ratio plot         
  TPad *pad2 = new TPad("pad2","",0,0.0,1,0.25);
  pad2->Draw(); pad2->cd();                             
  pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
  pad2->SetTitle(0);                          
  //pad2->SetRightMargin(padRmargin);           
  //pad2->SetLeftMargin(padLmargin);            
  pad2->SetBottomMargin(0.3);                 
  pad2->SetTopMargin(0);                                

  h_ratio->Draw("EP");                        
  ratioLine->Draw("SAME");                    
  h_ratio->Draw("sameEP");                    
  c1->cd();                                   

  c1->SaveAs("Data_JetPhox_comparison.pdf");
  delete c1;                    



}


