getCutFlowTable(){

  const Int_t nFiles = 4;
  TString FileName[nFiles] = {
    "rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root",
    "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root",
    "rootFiles/DiJet_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root",
    "rootFiles/QstarToGJ_M2000_f1p0_1.root"
  }

  TFile* fileData     = new TFile(FileName[0],"READ");
  TFile* fileBkgGJ    = new TFile(FileName[1],"READ");
  TFile* fileBkgDiJet = new TFile(FileName[2],"READ");
  TFile* fileSignal   = new TFile(FileName[3],"READ");

  TH1F * h_data   = (TH1F*)fileData     ->Get("h_CutExpFlow");
  TH1F * h_GJ     = (TH1F*)fileBkgGJ    ->Get("h_CutExpFlow");
  TH1F * h_DiJet  = (TH1F*)fileBkgDiJet ->Get("h_CutExpFlow");
  TH1F * h_Signal = (TH1F*)fileSignal   ->Get("h_CutExpFlow");

  Int_t nbins = h_Signal->GetNbinsX();

  std::vector<Double_t> ExpEntries_Sig, ExpEntries_GJ, ExpEntries_DiJet, ExpEntries_MC, ExpEntries_Data, ExpEntries_SbyB;
  std::vector<TString> Labels;

  ofstream fout("CutFlowTable_Pt190jetEta2p4dphi2p5dEta2p0M560LID.tex");
  fout<<"Cuts     &  Signal  &  PhotonJet   &   DiJet   &   Total Bkg   &    Data  &  \$S/\\sqrt{B}\$   \\\\"<<endl<<"\\hline"<<endl;

  for(Int_t k = 1; k <= nbins ; ++k){
  ExpEntries_Sig   .push_back(  h_Signal ->GetBinContent(k) );
  ExpEntries_GJ    .push_back(  h_GJ     ->GetBinContent(k) );
  ExpEntries_DiJet .push_back(  h_DiJet  ->GetBinContent(k) );
  ExpEntries_MC    .push_back(  h_GJ->GetBinContent(k) + h_DiJet->GetBinContent(k) );
  ExpEntries_Data  .push_back(  h_data->GetBinContent(k)    );
  ExpEntries_SbyB  .push_back(  (ExpEntries_Sig[k-1])/(pow(ExpEntries_MC[k-1], 0.5)) );
  Labels           .push_back(  h_data->GetXaxis()->GetBinLabel(k));
 
  fout<<Labels[k-1]<<" &  "<<ExpEntries_Sig[k-1]<<" &  "<<ExpEntries_GJ[k-1]<<" &  "<<ExpEntries_DiJet[k-1]<<" &  "<<ExpEntries_MC[k-1]<<" &  "<<ExpEntries_Data[k-1]<<" &  "<<ExpEntries_SbyB[k-1]<<" \\"<<"\\"<<endl;

  }
  fout.close();


}
