compilePlotsForLimit(){


  const Int_t nVar = 2;
  TString plotVar[nVar] = {
    "h_mass_VarBin_MassCut",  // Signal -0
    "h_mass_bin1_MassCut"     // Data   -1
  };

  TString OutFile = "Info_Pt190jetEta2p4dphi2p0M560.root";
  TFile* OutputFile = new TFile(OutFile,"RECREATE");
  

  TString InDataFile = "rootFiles/Run2015D_PromptReco_metFilters.root";

  const Int_t nSigFiles = 9;
  TString InSignalPath = "rootFiles/";
  TString InSignalFiles[nSigFiles] = {
    "QstarToGJ_M1000_f1p0_1.root",
    "QstarToGJ_M2000_f1p0_1.root",
    "QstarToGJ_M3000_f1p0_1.root",
    "QstarToGJ_M4000_f1p0_1.root",
    "QstarToGJ_M5000_f1p0_1.root",
    "QstarToGJ_M6000_f1p0_1.root",
    "QstarToGJ_M7000_f1p0_1.root",
    "QstarToGJ_M8000_f1p0_1.root",
    "QstarToGJ_M9000_f1p0_1.root"
  };

  TString SignalHistName[nSigFiles] = {"Qstar1000", "Qstar2000", "Qstar3000", "Qstar4000", "Qstar5000", "Qstar6000", "Qstar7000", "Qstar8000", "Qstar9000"};
  TString SignalCoupling =  "_f1p0";  // _f0p5, _f0p1

  //Data ----
  TFile* dataFile = new TFile(InDataFile,"READ");
  TH1F * h_Data_massBin1 = (TH1F*)dataFile->Get(plotVar[1]);

  OutputFile->cd();
  h_Data_massBin1 -> Write("h_Data_massBin1");
  h_Data_massBin1 -> Scale(1.0/h_Data_massBin1->Integral());
  h_Data_massBin1 -> Write("h_Data_massBin1");


  //Signal ----
  for(Int_t i = 0 ; i < nSigFiles ; ++i){
    TFile* signalFile = new TFile(InSignalPath+InSignalFiles[i],"READ");
    TH1F * h_Qstar_massVarBin = (TH1F*)signalFile->Get(plotVar[0]);
    
    OutputFile->cd();
    h_Qstar_massVarBin->Write("h_"+SignalHistName[i]+SignalCoupling+"_massVarBin");
    h_Qstar_massVarBin->Scale(1.0/h_Qstar_massVarBin->Integral());
    h_Qstar_massVarBin->Write("h_Prob"+SignalHistName[i]+SignalCoupling+"_massVarBin");

  
  }

}
