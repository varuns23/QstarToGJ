#include "makePlots.h"

Int_t main(){

  inputFilesPath = "rootFiles";
  prefix = "ptPhotJet190etaJet2p4dPhi1p5dEta1p8M695LID_76X"; // for Input root Files and also for output files (if any)

  plotsLocation  = "Plots";
  tablesLocation = "Tables";

  /// -- Input Files --- Nothing to change here -------------
  // Order matter for following lines -- do not mess with it.
  getInputFileNames(inputFilesPath, prefix); //calls the input files
  TFile * Inputfiles[fileNames.size()];
  TFile *sigFile_f1p0[sigName_f1p0.size()];
  TFile *sigFile_f0p5[sigName_f0p5.size()];
  TFile *sigFile_f0p1[sigName_f0p1.size()];
  // All scale factor etc. in Init function  
  getInitialized(Inputfiles, sigFile_f1p0, sigFile_f0p5, sigFile_f0p1);
  ////---------------------------------------------------------------------

///  // Make 1-D Histos:
///  draw1DHistos(Inputfiles, "ptPhoton_MassCut", 170., 2000.0, 0.01, 1e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "ptJet_MassCut", 170., 2000.0, 0.01, 1e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "mass_VarBin_MassCut", 300., 4000.0, 0.01, 1e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "etaPhoton_MassCut", -1.7, 1.7, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "etaJet_MassCut", -2.5, 2.5, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "phiPhoton_MassCut", -3.2, 3.2, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "phiJet_MassCut", -3.2, 3.2, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "PFMet_MassCut", 0.0, 1000.0, 0.01, 1e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "MetBySumMET_MassCut", 0.0, 0.5, 0.01, 1e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "DR_PhotonJet_MassCut", 2.4, 4.4, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "cosThetaStar_MassCut", 0.0, 0.8, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "dEta_MassCut", 0.0, 2.2, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "dphi_MassCut", 1.45, 3.2, 0.01, 1e5, "bottomCenter", true, 0);
///
///  draw1DHistos(Inputfiles, "Photon_SigmaIetaIeta_MassCut", 0.0, 0.018, 0.01, 1e4, "topRight", false, 0);
///  draw1DHistos(Inputfiles, "HoE_MassCut", 0.0, 0.06, 0.01, 1e5, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "CorrPFiso_Charged_MassCut", 0.0, 3.5, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "CorrPFiso_Neutral_MassCut", 0.0, 50.0, 0.01, 5e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "CorrPFiso_Photon_MassCut", 0.0, 8.0, 0.01, 5e4, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "PFiso_Electronveto_MassCut", 0.0, 3.0, 0.01, 1e5, "topRight", true, 0);
///
///  draw1DHistos(Inputfiles, "jet_NEF_MassCut", 0.0, 1.2, 0.01, 1e5, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "jet_NHF_MassCut", 0.0, 1.2, 0.01, 1e5, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "jet_CEF_MassCut", 0.0, 1.2, 0.01, 1e5, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "jet_CHF_MassCut", 0.0, 1.2, 0.01, 1e5, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "jet_NConstituents_MassCut", 0.0, 100.0, 0.01, 1e4, "bottomCenter", true, 0);
///  draw1DHistos(Inputfiles, "jet_ChargeMultiplicity_MassCut", 0.0, 100.0, 0.01, 1e4, "bottomCenter", true, 0);
///
///  draw1DHistos(Inputfiles, "nPhoton_MassCut", 0.0, 5.0, 0.01, 1e5, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "nIsoPhoton_MassCut", 0.0, 5.0, 0.01, 1e5, "topRight", true, 0);
///  draw1DHistos(Inputfiles, "nJet_MassCut", 0.0, 15.0, 0.01, 1e5, "topRight", true, 0);
///
///
/////  // Make 2-D Histos
///  for(Int_t ifile=0; ifile < fileNames.size(); ++ifile){
/////    draw2DHistos(Inputfiles[ifile], "PtPhotJet_MassCut",    150., 2000., 150., 2000., labelFiles[ifile]);
/////    draw2DHistos(Inputfiles[ifile], "etaPhotJet_MassCut",   -1.5,   1.5, -2.5,   2.5, labelFiles[ifile]);
///    draw2DHistos(Inputfiles[ifile], "pfMetgjMass_MassCut",  500., 4000.,   0., 1000., labelFiles[ifile]);
/////    draw2DHistos(Inputfiles[ifile], "nIsoPhotonPt_MassCut", 150., 2000.,   0.,   10., labelFiles[ifile]);
/////    draw2DHistos(Inputfiles[ifile], "nPhotonPt_MassCut",    150., 2000.,   0.,   10., labelFiles[ifile]);
///    draw2DHistos(Inputfiles[ifile], "nJetPt_MassCut",       150., 2000.,   0.,   20., labelFiles[ifile]);
///  }
//
//  // Compare two histos
//  compareHistos( "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root", "h_mass_VarBin_noPU_MassCut", "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root", "h_mass_VarBin_MassCut", "mass_VarBin_puComparison", 500.0, 4000.0, true, 1);
//  compareHistos( "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root", "h_ptPhoton_noPU_MassCut", "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root", "h_ptPhoton_MassCut", "ptPhoton_puComparison", 150.0, 2000.0, true, 1);
//  compareHistos( "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root", "h_ptJet_noPU_MassCut", "rootFiles/GJ_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID.root", "h_ptJet_MassCut", "ptJet_puComparison", 150.0, 2000.0, true, 1);
//
//
//  //Compare N-histos
//  std::vector<TString> fileName;
//  fileName.push_back("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560LID_74X.root");
//  fileName.push_back("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560MID_74X.root");
//  fileName.push_back("rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta2p0M560TID_74X.root");
//
//  std::vector<TString> histvar;
//  histvar.push_back("h_mass_VarBin_MassCut");
//  histvar.push_back("h_mass_VarBin_MassCut");
//  histvar.push_back("h_mass_VarBin_MassCut");
//
//  std::vector<TString> legName;
//  legName.push_back("74X - LID");
//  legName.push_back("74X - MID");
//  legName.push_back("74X - TID");
//  //compareNhistos( fileName, histvar, "74X-LMT_Comparison", legName,  500.0, 4000.0, true);
//
//
//  getPileUpHist(Inputfiles, "Vertices");
//  getPileUpHist(Inputfiles, "Vertices_MCut");

//  getHLTturnON(Inputfiles[0]);
//  getCutFlowTable(Inputfiles);
//  getSignalEfficiencies(sigName_f1p0, sigFile_f1p0,"f1p0");
//  getSignalEfficiencies(sigName_f0p5, sigFile_f0p5,"f0p5");
//  getSignalEfficiencies(sigName_f0p1, sigFile_f0p1,"f0p1");

} //-main
