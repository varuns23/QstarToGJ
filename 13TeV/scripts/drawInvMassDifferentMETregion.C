{
 TFile * fData[4];
 fData[0] = new TFile("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root","READ");
 fData[1] = new TFile("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_belowMET100.root","READ");
 fData[2] = new TFile("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_betweenMET100and200.root","READ");
 fData[3] = new TFile("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_betweenMET200and300.root","READ");
 fData[4] = new TFile("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_aboveMET300.root","READ");

 TH1F * h_Data[5];
 for(int ihist = 0; ihist < 5; ++ihist)
   h_Data[ihist] = (TH1F*)fData[ihist]->Get("h_mass_VarBin_MassCut");

 TLegend *leg = new TLegend(0.4445, 0.6, 0.9, 0.9);
 leg->SetBorderSize(0);
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetTextSize(0.03);

 h_Data[0]->SetStats(0);
 h_Data[0]->SetMarkerColor(kBlack);
 h_Data[0]->SetMarkerStyle(20);
 h_Data[0]->SetMarkerSize(0.8);
 h_Data[0]->SetLineColor(kBlack);
 h_Data[0]->GetXaxis()->SetRangeUser(500.0, 4000.0);
 leg->AddEntry(h_Data[0], "Data (no MET cut)", "LP");

 h_Data[1]->SetStats(0);
 h_Data[1]->SetLineColor(kTeal);
 h_Data[1]->SetFillColor(kTeal);
 leg->AddEntry(h_Data[1], "Data (MET < 100 GeV)", "F");

 h_Data[2]->SetStats(0);
 h_Data[2]->SetLineColor(kOrange);
 h_Data[2]->SetFillColor(kOrange);
 leg->AddEntry(h_Data[2], "Data (100 <= MET < 200 GeV)", "F");

 h_Data[3]->SetStats(0);
 h_Data[3]->SetLineColor(kRed);
 h_Data[3]->SetFillColor(kRed);
 leg->AddEntry(h_Data[3], "Data (200 <= MET < 300 GeV)", "F");

 h_Data[4]->SetStats(0);
 h_Data[4]->SetLineColor(kMagenta);
 h_Data[4]->SetFillColor(kMagenta);
 leg->AddEntry(h_Data[4], "Data (MET > 300 GeV)", "F");

 THStack *hS_ = new THStack("hS_", "");
 hS_->Add(h_Data[4]);
 hS_->Add(h_Data[3]);
 hS_->Add(h_Data[2]);
 hS_->Add(h_Data[1]);

//
 TCanvas * c1 = new TCanvas("c1", "",80, 20, 700, 700);
 c1->cd();
 c1->SetLogy();
 c1->SetFillColor(0);
 c1->SetFrameBorderMode(0);
 c1->SetTitle(0);
 c1->SetRightMargin(0.05);

 h_Data[0]->Draw("P");
 hS_->Draw("sameHIST");
 h_Data[0]->Draw("sameP");

 leg->Draw();

c1->SaveAs("Mass_DifferentMETregion_Data.pdf");



}
