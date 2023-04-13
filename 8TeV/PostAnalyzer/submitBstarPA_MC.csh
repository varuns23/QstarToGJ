#!/bin/tcsh

#####################################################
## This is a postanalyzer for b* -> gamma + b study
## author : Varun Sharma
## Email : varun.sharma@cern.ch
#####################################################

setenv pwd $PWD

## 0=Bstarf1p0,  1=Bstarf0p5,  2=Bstarf0p1,  3=Bkg GJ,   4=Bkg Dijet
#set runCase = 1 
foreach runCase ( 0 1 2  3 4)
set sampleIndex = 0 


setenv OutEos /eos/uscms/store/user/lpcqstar/Varun/Bstar8TeVResults/PostAnalyzer
##if( -d !${OutEos}) then
##mkdir -p ${OutEos}
##endif


###----------------------------------------------------------------------------------------------------
if( ${runCase} == 0 ) then
echo "******** Running for Signal -- Bstar, f=1.0 ***********"
setenv tmp BstarToGJ_
foreach i (${tmp}M500_f1p0 ${tmp}M700_f1p0 ${tmp}M1000_f1p0 ${tmp}M1200_f1p0 ${tmp}M1500_f1p0 ${tmp}M1700_f1p0 ${tmp}M2000_f1p0 ${tmp}M2500_f1p0 ${tmp}M3000_f1p0 ${tmp}M4000_f1p0)
set XS  = (1.348e0         2.096e-1        2.291e-2         6.52e-3          1.214e-3         4.321e-4         1.013e-4         1.072e-5         1.34e-6          3.707e-8)   
set totalEvents = (60000   60000           60000            120000           120000           120000           118402           160000           160000           198922 )    ## Confirm with Rocky if any Signal MC jobs failed
setenv sourceDir /eos/uscms/store/user/lpcqstar/NTuples/MC/BstarSignal/${i}/
setenv OutPath ${OutEos}/Bstarf1p0
set filesPerJob = 100 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 1 ) then
echo "******** Running for Signal -- Bstar, f=0.5 ***********"
setenv tmp BstarToGJ_
foreach i (${tmp}M500_f0p5 ${tmp}M700_f0p5 ${tmp}M1000_f0p5 ${tmp}M1200_f0p5 ${tmp}M1500_f0p5 ${tmp}M1700_f0p5 ${tmp}M2000_f0p5 ${tmp}M2500_f0p5 ${tmp}M3000_f0p5 ${tmp}M4000_f0p5)
set XS  = (3.415e-1        5.305e-2        5.841e-3         1.667e-3         3.057e-4         1.087e-4         2.505e-5         2.543e-6         2.91e-7          5.162e-9 )
set totalEvents = (60000   60000           120000           119626           120000           120000           160000           160000           198200           197291 )    ## Confirm with Rocky if any Signal MC jobs failed
setenv sourceDir /eos/uscms/store/user/lpcqstar/NTuples/MC/BstarSignal/${i}/
setenv OutPath ${OutEos}/Bstarf0p5
set filesPerJob = 100 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 2 ) then
echo "******** Running for Signal -- Bstar, f=0.1 ***********"
setenv tmp BstarToGJ_
foreach i (${tmp}M500_f0p1 ${tmp}M700_f0p1 ${tmp}M1000_f0p1 ${tmp}M1200_f0p1 ${tmp}M1500_f0p1 ${tmp}M1700_f0p1 ${tmp}M2000_f0p1 ${tmp}M2500_f0p1 ${tmp}M3000_f0p1 ${tmp}M4000_f0p1)
set XS  = (1.377e-2        2.143e-3        2.366e-4         6.691e-5         1.231e-5         4.339e-6         1.003e-6         1.001e-7         1.103e-8         1.531e-10 )
set totalEvents = (60000   120000          120000           120000           120000           160000           160000           198980           197872           196962 )    ## Confirm with Rocky if any Signal MC jobs failed
setenv sourceDir /eos/uscms/store/user/lpcqstar/NTuples/MC/BstarSignal/${i}/
setenv OutPath ${OutEos}/Bstarf0p1
set filesPerJob = 100 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 3 ) then
echo "******** Running for Bkg -- Photon+Jet Pythia ***********"
foreach i (G_Pt_120to170  G_Pt_170to300  G_Pt_300to470  G_Pt_470to800  G_Pt_800to1400  G_Pt_1400to1800  G_Pt_1800toInf )
set XS  = (108.0068       30.12207       2.138632       0.2119244      0.007077847     4.510327E-5      1.867141E-6)   
set totalEvents = (2000043.0 2000069.0   2000130.0      1975231.0      1973504.0       1993890.0        1939122.0)   ## Confirm with Rocky if any GJ MC jobs failed 
setenv sourceDir /eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/${i}/
setenv OutPath ${OutEos}/GJBkg
set filesPerJob = 100 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------
if( ${runCase} == 4 ) then
echo "******** Running for Bkg -- Dijet Pythia ***********"
setenv tmp QCDDiJet_Pt_
foreach i ( ${tmp}120to170  ${tmp}170to300  ${tmp}300to470  ${tmp}470to600  ${tmp}600to800  ${tmp}800to1000  ${tmp}1000to1400  ${tmp}1400to1800  ${tmp}1800toInf )
set XS  = ( 156293.3        34138.15        1759.549        113.8791        26.9921         3.550036         0.737844          0.03352235        0.001829005 )
set totalEvents = (5985732.0  5814398.0     5978500.0       3994848.0       3996864.0       3998563.0        1964088.0         2000062.0         977586.0)   ## Confirm with Rocky if any Dijet MC jobs failed 
setenv sourceDir /eos/uscms/store/user/lpcqstar/NTuples/MC/QCDDiJet/${i}/
setenv OutPath ${OutEos}/DiJetBkg
set filesPerJob = 100 ## Total files : max
endif
##
###----------------------------------------------------------------------------------------------------


setenv FileNameTag ${i}

@ sampleIndex = ${sampleIndex} + 1

echo " sampleIndex = ${sampleIndex} " 

end

##################MyEvent_MC.C#####################
cat>PostAnalyzerMC.C<<EOF
#define MyEvent_MC_cxx
#include "MyEvent_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MyEvent_MC::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
EOF


################################################################
cat>PostAnalyzerMC.h<<EOF
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 23 07:24:02 2015 by ROOT version 5.32/00
// from TTree myEvent/A Tree with Histograms
// found on file: /eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/G_Pt_120to170/AODSIM_G_Pt_120to170_110_1_dHz.root
//////////////////////////////////////////////////////////

#ifndef MyEvent_MC_h
#define MyEvent_MC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TTree.h>
#include <TBranch.h>
#include "TH1.h"
#include "TH2.h"
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <Riostream.h> 
#include <TGraphAsymmErrors.h>
#include <map>
#include <vector>

using namespace std;
using namespace ROOT;

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxhasMatchedGenPhoton = 1;
const Int_t kMaxhasMatchedGenParton_ToPFJet = 1;
const Int_t kMaxhasMatchedGenJet_ToPFJet = 1;
const Int_t kMaxIsScrapingEvent = 1;

class MyEvent_MC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          LumiNumber;
   UInt_t          BXNumber;
   UInt_t          totalIntensityBeam1;
   UInt_t          totalIntensityBeam2;
   Float_t         avgInsDelLumi;
   Float_t         avgInsDelLumiErr;
   Float_t         avgInsRecLumi;
   Float_t         avgInsRecLumiErr;
   Int_t           npuVertices;
   Int_t           npuVerticesm1;
   Int_t           npuVerticesp1;
   Int_t           npuVerticespm2;
   Float_t         trueNumofInteractions;
   Int_t           Photon_n;
   vector<double>  *Photon_E;
   vector<double>  *Photon_et;
   vector<double>  *Photon_pt;
   vector<double>  *Photon_eta;
   vector<double>  *Photon_phi;
   vector<double>  *Photon_theta;
   vector<double>  *Photon_px;
   vector<double>  *Photon_py;
   vector<double>  *Photon_pz;
   vector<double>  *Photon_vx;
   vector<double>  *Photon_vy;
   vector<double>  *Photon_vz;
   vector<float>   *Photon_r9;
   vector<float>   *Photon_maxEnergyXtal;
   vector<float>   *Photon_e1x5;
   vector<float>   *Photon_e2x5;
   vector<float>   *Photon_e3x3;
   vector<float>   *Photon_e5x5;
   vector<float>   *Photon_r1x5;
   vector<float>   *Photon_r2x5;
   vector<float>   *Photon_SigmaEtaEta;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_SigmaEtaPhi;
   vector<float>   *Photon_SigmaIEtaIPhi;
   vector<float>   *Photon_SigmaPhiPhi;
   vector<float>   *Photon_SigmaIPhiIPhi;
   vector<float>   *Photon_roundness;
   vector<float>   *Photon_angle;
   vector<float>   *Photon_swissCross;
   vector<float>   *Photon_s9;
   vector<float>   *Photon_e4Overe1;
   vector<float>   *Photon_e6Overe2;
   vector<float>   *Photon_e2Overe9;
   vector<float>   *Photon_rookFraction;
   vector<bool>    *Photon_isEB;
   vector<bool>    *Photon_isEE;
   vector<bool>    *Photon_isEBGap;
   vector<bool>    *Photon_isEEGap;
   vector<bool>    *Photon_isEBEEGap;
   vector<float>   *Photon_ecalRecHitSumEtConeDR03;
   vector<float>   *Photon_hcalTowerSumEtConeDR03;
   vector<float>   *Photon_hcalDepth1TowerSumEtConeDR03;
   vector<float>   *Photon_hcalDepth2TowerSumEtConeDR03;
   vector<float>   *Photon_trkSumPtSolidConeDR03;
   vector<float>   *Photon_trkSumPtHollowConeDR03;
   vector<int>     *Photon_nTrkSolidConeDR03;
   vector<int>     *Photon_nTrkHollowConeDR03;
   vector<float>   *Photon_ecalRecHitSumEtConeDR04;
   vector<float>   *Photon_hcalTowerSumEtConeDR04;
   vector<float>   *Photon_hcalDepth1TowerSumEtConeDR04;
   vector<float>   *Photon_hcalDepth2TowerSumEtConeDR04;
   vector<float>   *Photon_trkSumPtSolidConeDR04;
   vector<float>   *Photon_trkSumPtHollowConeDR04;
   vector<int>     *Photon_nTrkSolidConeDR04;
   vector<int>     *Photon_nTrkHollowConeDR04;
   vector<float>   *Photon_HoE;
   vector<float>   *Photon_SingleTowerHoE;
   vector<bool>    *Photon_hasConvTrk;
   vector<bool>    *Photon_hasPixelSeed;
   vector<bool>    *passedConvSafeElectronVeto;
   vector<int>     *Photon_SC_nOfBasicClusters;
   vector<double>  *Photon_SC_rawEnergy;
   vector<double>  *Photon_SC_preShowerEnergy;
   vector<double>  *Photon_SC_energy;
   vector<double>  *Photon_SC_eta;
   vector<double>  *Photon_SC_phi;
   vector<double>  *Photon_SC_x;
   vector<double>  *Photon_SC_y;
   vector<double>  *Photon_SC_z;
   vector<double>  *Photon_SC_etaWidth;
   vector<double>  *Photon_SC_phiWidth;
   vector<float>   *Photon_mipChi2;
   vector<float>   *Photon_mipTotEnergy;
   vector<float>   *Photon_mipSlope;
   vector<float>   *Photon_mipIntercept;
   vector<int>     *Photon_mipNhitCone;
   vector<bool>    *Photon_mipIsHalo;
   vector<unsigned int> *Photon_nConvTracks;
   vector<bool>    *Photon_isConverted;
   vector<float>   *Photon_pairInvariantMass;
   vector<float>   *Photon_pairCotThetaSeparation;
   vector<float>   *Photon_pairMomentum_x;
   vector<float>   *Photon_pairMomentum_y;
   vector<float>   *Photon_pairMomentum_z;
   vector<float>   *Photon_EoverP;
   vector<float>   *Photon_conv_vx;
   vector<float>   *Photon_conv_vy;
   vector<float>   *Photon_conv_vz;
   vector<float>   *Photon_zOfPrimaryVtxFromTrks;
   vector<float>   *Photon_distOfMinimumApproach;
   vector<float>   *Photon_dPhiTracksAtVtx;
   vector<float>   *Photon_dPhiTracksAtEcal;
   vector<float>   *Photon_dEtaTracksAtEcal;
   vector<int>     *Photon_nCrystals;
   vector<vector<float> > *Photon_xtal_timing;
   vector<vector<float> > *Photon_xtal_timeErr;
   vector<float>   *Photon_avgTimeAllxtals;
   vector<vector<float> > *Photon_xtal_energy;
   vector<vector<int> > *Photon_xtal_EBieta;
   vector<vector<int> > *Photon_xtal_EBiphi;
   vector<vector<int> > *Photon_xtal_EBrecoFlag;
   vector<double>  *PFIsoPhoton03;
   vector<double>  *PFIsoNeutral03;
   vector<double>  *PFIsoCharged03;
   vector<double>  *PFIsoSum03;
   vector<double>  *PFIsoChargedWorstvtx03;
   Int_t           GenPhoton_n;
   vector<int>     *GenPhoton_status;
   vector<double>  *GenPhoton_E;
   vector<double>  *GenPhoton_et;
   vector<double>  *GenPhoton_pt;
   vector<double>  *GenPhoton_eta;
   vector<double>  *GenPhoton_phi;
   vector<double>  *GenPhoton_px;
   vector<double>  *GenPhoton_py;
   vector<double>  *GenPhoton_pz;
   vector<double>  *GenPhoton_vx;
   vector<double>  *GenPhoton_vy;
   vector<double>  *GenPhoton_vz;
   vector<int>     *GenPhoton_mother_status;
   vector<int>     *GenPhoton_mother_pdgid;
   vector<double>  *GenPhoton_mother_E;
   vector<double>  *GenPhoton_mother_et;
   vector<double>  *GenPhoton_mother_pt;
   vector<double>  *GenPhoton_mother_eta;
   vector<double>  *GenPhoton_mother_phi;
   vector<double>  *GenPhoton_mother_px;
   vector<double>  *GenPhoton_mother_py;
   vector<double>  *GenPhoton_mother_pz;
   vector<double>  *GenPhoton_mother_vx;
   vector<double>  *GenPhoton_mother_vy;
   vector<double>  *GenPhoton_mother_vz;
   vector<int>     *GenPhoton_grandmother_status;
   vector<int>     *GenPhoton_grandmother_pdgid;
   vector<double>  *GenPhoton_grandmother_E;
   vector<double>  *GenPhoton_grandmother_et;
   vector<double>  *GenPhoton_grandmother_pt;
   vector<double>  *GenPhoton_grandmother_eta;
   vector<double>  *GenPhoton_grandmother_phi;
   vector<double>  *GenPhoton_grandmother_px;
   vector<double>  *GenPhoton_grandmother_py;
   vector<double>  *GenPhoton_grandmother_pz;
   vector<double>  *GenPhoton_grandmother_vx;
   vector<double>  *GenPhoton_grandmother_vy;
   vector<double>  *GenPhoton_grandmother_vz;
   vector<bool>    *hasMatchedGenPhoton_;
   vector<int>     *MatchedGenPhoton_status;
   vector<int>     *MatchedGenPhoton_pdgid;
   vector<double>  *MatchedGenPhoton_E;
   vector<double>  *MatchedGenPhoton_et;
   vector<double>  *MatchedGenPhoton_pt;
   vector<double>  *MatchedGenPhoton_eta;
   vector<double>  *MatchedGenPhoton_phi;
   vector<double>  *MatchedGenPhoton_px;
   vector<double>  *MatchedGenPhoton_py;
   vector<double>  *MatchedGenPhoton_pz;
   vector<double>  *MatchedGenPhoton_vx;
   vector<double>  *MatchedGenPhoton_vy;
   vector<double>  *MatchedGenPhoton_vz;
   vector<int>     *MatchedGenPhoton_nMothers;
   vector<vector<int> > *MatchedGenPhoton_mothersPdgid;
   vector<vector<int> > *MatchedGenPhoton_mothersStatus;
   Int_t           PFPatJet_n;
   vector<double>  *PFPatJet_E;
   vector<double>  *PFPatJet_et;
   vector<double>  *PFPatJet_pt;
   vector<double>  *PFPatJet_eta;
   vector<double>  *PFPatJet_phi;
   vector<double>  *PFPatJet_theta;
   vector<double>  *PFPatJet_px;
   vector<double>  *PFPatJet_py;
   vector<double>  *PFPatJet_pz;
   vector<double>  *PFPatJet_vx;
   vector<double>  *PFPatJet_vy;
   vector<double>  *PFPatJet_vz;
   vector<float>   *PFPatJet_ChargedEmEnergy;
   vector<float>   *PFPatJet_ChargedEmEnergyFrac;
   vector<float>   *PFPatJet_ChargedHadEnergy;
   vector<float>   *PFPatJet_ChargedHadEnergyFrac;
   vector<int>     *PFPatJet_ChargedHadMult;
   vector<int>     *PFPatJet_ChargedMult;
   vector<int>     *PFPatJet_NConstituents;
   vector<float>   *PFPatJet_HFEMEnergy;
   vector<float>   *PFPatJet_HFEMEnergyFrac;
   vector<int>     *PFPatJet_HFEMMult;
   vector<float>   *PFPatJet_HFHadEnergy;
   vector<float>   *PFPatJet_HFHadEnergyFrac;
   vector<int>     *PFPatJet_HFHadMult;
   vector<float>   *PFPatJet_NeutralEmEnergy;
   vector<float>   *PFPatJet_NeutralEmEnergyFrac;
   vector<float>   *PFPatJet_NeutralHadEnergy;
   vector<float>   *PFPatJet_NeutralHadEnergyFrac;
   vector<int>     *PFPatJet_NeutralHadMult;
   vector<int>     *PFPatJet_NeutralMult;
   vector<double>  *PFPatJet_jecUncertainity;
   vector<float>   *PFPatJet_puJetIdCutBased_MVA;
   vector<float>   *PFPatJet_puJetIdSimple_MVA;
   vector<float>   *PFPatJet_puJetIdFull_MVA;
   vector<bool>    *PFPatJet_PassPUJetIdCutBased_loose;
   vector<bool>    *PFPatJet_PassPUJetIdCutBased_medium;
   vector<bool>    *PFPatJet_PassPUJetIdCutBased_tight;
   vector<bool>    *PFPatJet_PassPUJetIdSimple_loose;
   vector<bool>    *PFPatJet_PassPUJetIdSimple_medium;
   vector<bool>    *PFPatJet_PassPUJetIdSimple_tight;
   vector<bool>    *PFPatJet_PassPUJetIdFull_loose;
   vector<bool>    *PFPatJet_PassPUJetIdFull_medium;
   vector<bool>    *PFPatJet_PassPUJetIdFull_tight;
   vector<float>   *PFPatJet_BJetDiscrByTrackCountingHighEff;
   vector<float>   *PFPatJet_BJetDiscrByTrackCountingHighPur;
   vector<float>   *PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff;
   vector<float>   *PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur;
   vector<float>   *PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff;
   vector<float>   *PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA;
   vector<float>   *PFPatJet_BJetDiscrByjetProbabilityBJetTags;
   vector<float>   *PFPatJet_BJetDiscrByjetBProbabilityBJetTags;
   vector<float>   *PFPatJet_BJetDiscrBySoftElectronBJetTags;
   vector<float>   *PFPatJet_BJetDiscrBySoftMuonBJetTags;
   vector<bool>    *hasMatchedGenParton_ToPFJet_;
   vector<bool>    *hasMatchedGenJet_ToPFJet_;
   vector<int>     *MatchedGenPartonToPFJet_pdgId;
   vector<int>     *MatchedGenPartonToPFJet_status;
   vector<double>  *MatchedGenJetToPFJet_E;
   vector<double>  *MatchedGenJetToPFJet_pt;
   vector<double>  *MatchedGenJetToPFJet_eta;
   vector<double>  *MatchedGenJetToPFJet_phi;
   vector<int>     *MatchedGenJetToPFJet_nDaughters;
   vector<int>     *PFPatJet_JetPartonFlavor;
   Int_t           PFRecoJet_n;
   vector<double>  *PFRecoJet_E;
   vector<double>  *PFRecoJet_et;
   vector<double>  *PFRecoJet_pt;
   vector<double>  *PFRecoJet_eta;
   vector<double>  *PFRecoJet_phi;
   vector<double>  *PFRecoJet_theta;
   vector<double>  *PFRecoJet_px;
   vector<double>  *PFRecoJet_py;
   vector<double>  *PFRecoJet_pz;
   vector<double>  *PFRecoJet_vx;
   vector<double>  *PFRecoJet_vy;
   vector<double>  *PFRecoJet_vz;
   vector<float>   *PFRecoJet_ChargedEmEnergy;
   vector<float>   *PFRecoJet_ChargedEmEnergyFrac;
   vector<float>   *PFRecoJet_ChargedHadEnergy;
   vector<float>   *PFRecoJet_ChargedHadEnergyFrac;
   vector<int>     *PFRecoJet_ChargedHadMult;
   vector<int>     *PFRecoJet_ChargedMult;
   vector<int>     *PFRecoJet_NConstituents;
   vector<float>   *PFRecoJet_HFEMEnergy;
   vector<float>   *PFRecoJet_HFEMEnergyFrac;
   vector<int>     *PFRecoJet_HFEMMult;
   vector<float>   *PFRecoJet_HFHadEnergy;
   vector<float>   *PFRecoJet_HFHadEnergyFrac;
   vector<int>     *PFRecoJet_HFHadMult;
   vector<float>   *PFRecoJet_NeutralEmEnergy;
   vector<float>   *PFRecoJet_NeutralEmEnergyFrac;
   vector<float>   *PFRecoJet_NeutralHadEnergy;
   vector<float>   *PFRecoJet_NeutralHadEnergyFrac;
   vector<int>     *PFRecoJet_NeutralHadMult;
   vector<int>     *PFRecoJet_NeutralMult;
   vector<double>  *PFRecoJet_jecUncertainity;
   vector<float>   *PFRecoJet_BJetDiscrByTrackCountingHighEff;
   vector<float>   *PFRecoJet_BJetDiscrByTrackCountingHighPur;
   vector<float>   *PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff;
   vector<float>   *PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur;
   vector<float>   *PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff;
   vector<float>   *PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA;
   vector<float>   *PFRecoJet_BJetDiscrByjetProbabilityBJetTags;
   vector<float>   *PFRecoJet_BJetDiscrByjetBProbabilityBJetTags;
   vector<float>   *PFRecoJet_BJetDiscrBySoftElectronBJetTags;
   vector<float>   *PFRecoJet_BJetDiscrBySoftMuonBJetTags;
   Int_t           GenJet_n;
   vector<double>  *GenJet_E;
   vector<double>  *GenJet_et;
   vector<double>  *GenJet_pt;
   vector<double>  *GenJet_eta;
   vector<double>  *GenJet_phi;
   vector<double>  *GenJet_px;
   vector<double>  *GenJet_py;
   vector<double>  *GenJet_pz;
   vector<double>  *GenJet_vx;
   vector<double>  *GenJet_vy;
   vector<double>  *GenJet_vz;
   vector<int>     *GenJet_nDaughters;
   vector<vector<int> > *GenJet_daughters_status;
   vector<vector<int> > *GenJet_daughters_pdgid;
   vector<vector<double> > *GenJet_daughters_E;
   vector<vector<double> > *GenJet_daughters_pt;
   vector<vector<double> > *GenJet_daughters_eta;
   vector<vector<double> > *GenJet_daughters_phi;
   Int_t           Vertex_n;
   vector<double>  *Vertex_x;
   vector<double>  *Vertex_y;
   vector<double>  *Vertex_z;
   vector<double>  *Vertex_chi2;
   vector<double>  *Vertex_nchi2;
   vector<double>  *Vertex_ndof;
   vector<int>     *Vertex_tracksSize;
   vector<bool>    *Vertex_isFake;
   vector<bool>    *Vertex_isValid;
   vector<double>  *Vertex_d0;
   Int_t           Tracks_n;
   vector<double>  *Track_pt;
   vector<double>  *Track_px;
   vector<double>  *Track_py;
   vector<double>  *Track_pz;
   vector<double>  *Track_vx;
   vector<double>  *Track_vy;
   vector<double>  *Track_vz;
   vector<double>  *Track_eta;
   vector<double>  *Track_phi;
   vector<double>  *Track_theta;
   vector<double>  *Track_chi2;
   Bool_t          IsScrapingEvent_;
   Float_t         Scraping_FractionOfGoodTracks;
   vector<string>  *HLT_Photon_triggers;
   vector<int>     *HLT_Photon_trig_prescales;
   vector<bool>    *HLT_Photon_ifTriggerPassed;
   Int_t           HLT_Photon_nTriggers;
   vector<int>     *HLT_Photon_triggerIndex;
   vector<int>     *HLT_Photon_nFilters;
   vector<string>  *HLT_Photon_FilterNames;
   vector<int>     *HLT_Photon_trigger_FilterStartPosition;
   vector<int>     *HLT_Photon_trigger_FilterEndPosition;
   vector<vector<double> > *HLT_Photon_FilterObjects_pt;
   vector<vector<double> > *HLT_Photon_FilterObjects_eta;
   vector<vector<double> > *HLT_Photon_FilterObjects_phi;
   Float_t         rho;
   Float_t         sigma;
   Float_t         rho25;
   Float_t         sigma25;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_LumiNumber;   //!
   TBranch        *b_BXNumber;   //!
   TBranch        *b_totalIntensityBeam1;   //!
   TBranch        *b_totalIntensityBeam2;   //!
   TBranch        *b_avgInsDelLumi;   //!
   TBranch        *b_avgInsDelLumiErr;   //!
   TBranch        *b_avgInsRecLumi;   //!
   TBranch        *b_avgInsRecLumiErr;   //!
   TBranch        *b_npuVertices;   //!
   TBranch        *b_npuVerticesm1;   //!
   TBranch        *b_npuVerticesp1;   //!
   TBranch        *b_npuVerticespm2;   //!
   TBranch        *b_trueNumofInteractions;   //!
   TBranch        *b_Photon_n;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_theta;   //!
   TBranch        *b_Photon_px;   //!
   TBranch        *b_Photon_py;   //!
   TBranch        *b_Photon_pz;   //!
   TBranch        *b_Photon_vx;   //!
   TBranch        *b_Photon_vy;   //!
   TBranch        *b_Photon_vz;   //!
   TBranch        *b_Photon_r9;   //!
   TBranch        *b_Photon_maxEnergyXtal;   //!
   TBranch        *b_Photon_e1x5;   //!
   TBranch        *b_Photon_e2x5;   //!
   TBranch        *b_Photon_e3x3;   //!
   TBranch        *b_Photon_e5x5;   //!
   TBranch        *b_Photon_r1x5;   //!
   TBranch        *b_Photon_r2x5;   //!
   TBranch        *b_Photon_SigmaEtaEta;   //!
   TBranch        *b_Photon_SigmaIEtaIEta;   //!
   TBranch        *b_Photon_SigmaEtaPhi;   //!
   TBranch        *b_Photon_SigmaIEtaIPhi;   //!
   TBranch        *b_Photon_SigmaPhiPhi;   //!
   TBranch        *b_Photon_SigmaIPhiIPhi;   //!
   TBranch        *b_Photon_roundness;   //!
   TBranch        *b_Photon_angle;   //!
   TBranch        *b_Photon_swissCross;   //!
   TBranch        *b_Photon_s9;   //!
   TBranch        *b_Photon_e4Overe1;   //!
   TBranch        *b_Photon_e6Overe2;   //!
   TBranch        *b_Photon_e2Overe9;   //!
   TBranch        *b_Photon_rookFraction;   //!
   TBranch        *b_Photon_isEB;   //!
   TBranch        *b_Photon_isEE;   //!
   TBranch        *b_Photon_isEBGap;   //!
   TBranch        *b_Photon_isEEGap;   //!
   TBranch        *b_Photon_isEBEEGap;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR03;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR03;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
   TBranch        *b_Photon_nTrkSolidConeDR03;   //!
   TBranch        *b_Photon_nTrkHollowConeDR03;   //!
   TBranch        *b_Photon_ecalRecHitSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
   TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
   TBranch        *b_Photon_nTrkSolidConeDR04;   //!
   TBranch        *b_Photon_nTrkHollowConeDR04;   //!
   TBranch        *b_Photon_HoE;   //!
   TBranch        *b_Photon_SingleTowerHoE;   //!
   TBranch        *b_Photon_hasConvTrk;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!
   TBranch        *b_passedConvSafeElectronVeto;   //!
   TBranch        *b_Photon_SC_nOfBasicClusters;   //!
   TBranch        *b_Photon_SC_rawEnergy;   //!
   TBranch        *b_Photon_SC_preShowerEnergy;   //!
   TBranch        *b_Photon_SC_energy;   //!
   TBranch        *b_Photon_SC_eta;   //!
   TBranch        *b_Photon_SC_phi;   //!
   TBranch        *b_Photon_SC_x;   //!
   TBranch        *b_Photon_SC_y;   //!
   TBranch        *b_Photon_SC_z;   //!
   TBranch        *b_Photon_SC_etaWidth;   //!
   TBranch        *b_Photon_SC_phiWidth;   //!
   TBranch        *b_Photon_mipChi2;   //!
   TBranch        *b_Photon_mipTotEnergy;   //!
   TBranch        *b_Photon_mipSlope;   //!
   TBranch        *b_Photon_mipIntercept;   //!
   TBranch        *b_Photon_mipNhitCone;   //!
   TBranch        *b_Photon_mipIsHalo;   //!
   TBranch        *b_Photon_nConvTracks;   //!
   TBranch        *b_Photon_isConverted;   //!
   TBranch        *b_Photon_pairInvariantMass;   //!
   TBranch        *b_Photon_pairCotThetaSeparation;   //!
   TBranch        *b_Photon_pairMomentum_x;   //!
   TBranch        *b_Photon_pairMomentum_y;   //!
   TBranch        *b_Photon_pairMomentum_z;   //!
   TBranch        *b_Photon_EoverP;   //!
   TBranch        *b_Photon_conv_vx;   //!
   TBranch        *b_Photon_conv_vy;   //!
   TBranch        *b_Photon_conv_vz;   //!
   TBranch        *b_Photon_zOfPrimaryVtxFromTrks;   //!
   TBranch        *b_Photon_distOfMinimumApproach;   //!
   TBranch        *b_Photon_dPhiTracksAtVtx;   //!
   TBranch        *b_Photon_dPhiTracksAtEcal;   //!
   TBranch        *b_Photon_dEtaTracksAtEcal;   //!
   TBranch        *b_Photon_nCrystals;   //!
   TBranch        *b_Photon_xtal_timing;   //!
   TBranch        *b_Photon_xtal_timeErr;   //!
   TBranch        *b_Photon_avgTimeAllxtals;   //!
   TBranch        *b_Photon_xtal_energy;   //!
   TBranch        *b_Photon_xtal_EBieta;   //!
   TBranch        *b_Photon_xtal_EBiphi;   //!
   TBranch        *b_Photon_xtal_EBrecoFlag;   //!
   TBranch        *b_PFIsoPhoton03;   //!
   TBranch        *b_PFIsoNeutral03;   //!
   TBranch        *b_PFIsoCharged03;   //!
   TBranch        *b_PFIsoSum03;   //!
   TBranch        *b_PFIsoChargedWorstvtx03;   //!
   TBranch        *b_GenPhoton_n;   //!
   TBranch        *b_GenPhoton_status;   //!
   TBranch        *b_GenPhoton_E;   //!
   TBranch        *b_GenPhoton_et;   //!
   TBranch        *b_GenPhoton_pt;   //!
   TBranch        *b_GenPhoton_eta;   //!
   TBranch        *b_GenPhoton_phi;   //!
   TBranch        *b_GenPhoton_px;   //!
   TBranch        *b_GenPhoton_py;   //!
   TBranch        *b_GenPhoton_pz;   //!
   TBranch        *b_GenPhoton_vx;   //!
   TBranch        *b_GenPhoton_vy;   //!
   TBranch        *b_GenPhoton_vz;   //!
   TBranch        *b_GenPhoton_mother_status;   //!
   TBranch        *b_GenPhoton_mother_pdgid;   //!
   TBranch        *b_GenPhoton_mother_E;   //!
   TBranch        *b_GenPhoton_mother_et;   //!
   TBranch        *b_GenPhoton_mother_pt;   //!
   TBranch        *b_GenPhoton_mother_eta;   //!
   TBranch        *b_GenPhoton_mother_phi;   //!
   TBranch        *b_GenPhoton_mother_px;   //!
   TBranch        *b_GenPhoton_mother_py;   //!
   TBranch        *b_GenPhoton_mother_pz;   //!
   TBranch        *b_GenPhoton_mother_vx;   //!
   TBranch        *b_GenPhoton_mother_vy;   //!
   TBranch        *b_GenPhoton_mother_vz;   //!
   TBranch        *b_GenPhoton_grandmother_status;   //!
   TBranch        *b_GenPhoton_grandmother_pdgid;   //!
   TBranch        *b_GenPhoton_grandmother_E;   //!
   TBranch        *b_GenPhoton_grandmother_et;   //!
   TBranch        *b_GenPhoton_grandmother_pt;   //!
   TBranch        *b_GenPhoton_grandmother_eta;   //!
   TBranch        *b_GenPhoton_grandmother_phi;   //!
   TBranch        *b_GenPhoton_grandmother_px;   //!
   TBranch        *b_GenPhoton_grandmother_py;   //!
   TBranch        *b_GenPhoton_grandmother_pz;   //!
   TBranch        *b_GenPhoton_grandmother_vx;   //!
   TBranch        *b_GenPhoton_grandmother_vy;   //!
   TBranch        *b_GenPhoton_grandmother_vz;   //!
   TBranch        *b_hasMatchedGenPhoton_;   //!
   TBranch        *b_MatchedGenPhoton_status;   //!
   TBranch        *b_MatchedGenPhoton_pdgid;   //!
   TBranch        *b_MatchedGenPhoton_E;   //!
   TBranch        *b_MatchedGenPhoton_et;   //!
   TBranch        *b_MatchedGenPhoton_pt;   //!
   TBranch        *b_MatchedGenPhoton_eta;   //!
   TBranch        *b_MatchedGenPhoton_phi;   //!
   TBranch        *b_MatchedGenPhoton_px;   //!
   TBranch        *b_MatchedGenPhoton_py;   //!
   TBranch        *b_MatchedGenPhoton_pz;   //!
   TBranch        *b_MatchedGenPhoton_vx;   //!
   TBranch        *b_MatchedGenPhoton_vy;   //!
   TBranch        *b_MatchedGenPhoton_vz;   //!
   TBranch        *b_MatchedGenPhoton_nMothers;   //!
   TBranch        *b_MatchedGenPhoton_mothersPdgid;   //!
   TBranch        *b_MatchedGenPhoton_mothersStatus;   //!
   TBranch        *b_PFPatJet_n;   //!
   TBranch        *b_PFPatJet_E;   //!
   TBranch        *b_PFPatJet_et;   //!
   TBranch        *b_PFPatJet_pt;   //!
   TBranch        *b_PFPatJet_eta;   //!
   TBranch        *b_PFPatJet_phi;   //!
   TBranch        *b_PFPatJet_theta;   //!
   TBranch        *b_PFPatJet_px;   //!
   TBranch        *b_PFPatJet_py;   //!
   TBranch        *b_PFPatJet_pz;   //!
   TBranch        *b_PFPatJet_vx;   //!
   TBranch        *b_PFPatJet_vy;   //!
   TBranch        *b_PFPatJet_vz;   //!
   TBranch        *b_PFPatJet_ChargedEmEnergy;   //!
   TBranch        *b_PFPatJet_ChargedEmEnergyFrac;   //!
   TBranch        *b_PFPatJet_ChargedHadEnergy;   //!
   TBranch        *b_PFPatJet_ChargedHadEnergyFrac;   //!
   TBranch        *b_PFPatJet_ChargedHadMult;   //!
   TBranch        *b_PFPatJet_ChargedMult;   //!
   TBranch        *b_PFPatJet_NConstituents;   //!
   TBranch        *b_PFPatJet_HFEMEnergy;   //!
   TBranch        *b_PFPatJet_HFEMEnergyFrac;   //!
   TBranch        *b_PFPatJet_HFEMMult;   //!
   TBranch        *b_PFPatJet_HFHadEnergy;   //!
   TBranch        *b_PFPatJet_HFHadEnergyFrac;   //!
   TBranch        *b_PFPatJet_HFHadMult;   //!
   TBranch        *b_PFPatJet_NeutralEmEnergy;   //!
   TBranch        *b_PFPatJet_NeutralEmEnergyFrac;   //!
   TBranch        *b_PFPatJet_NeutralHadEnergy;   //!
   TBranch        *b_PFPatJet_NeutralHadEnergyFrac;   //!
   TBranch        *b_PFPatJet_NeutralHadMult;   //!
   TBranch        *b_PFPatJet_NeutralMult;   //!
   TBranch        *b_PFPatJet_jecUncertainity;   //!
   TBranch        *b_PFPatJet_puJetIdCutBased_MVA;   //!
   TBranch        *b_PFPatJet_puJetIdSimple_MVA;   //!
   TBranch        *b_PFPatJet_puJetIdFull_MVA;   //!
   TBranch        *b_PFPatJet_PassPUJetIdCutBased_loose;   //!
   TBranch        *b_PFPatJet_PassPUJetIdCutBased_medium;   //!
   TBranch        *b_PFPatJet_PassPUJetIdCutBased_tight;   //!
   TBranch        *b_PFPatJet_PassPUJetIdSimple_loose;   //!
   TBranch        *b_PFPatJet_PassPUJetIdSimple_medium;   //!
   TBranch        *b_PFPatJet_PassPUJetIdSimple_tight;   //!
   TBranch        *b_PFPatJet_PassPUJetIdFull_loose;   //!
   TBranch        *b_PFPatJet_PassPUJetIdFull_medium;   //!
   TBranch        *b_PFPatJet_PassPUJetIdFull_tight;   //!
   TBranch        *b_PFPatJet_BJetDiscrByTrackCountingHighEff;   //!
   TBranch        *b_PFPatJet_BJetDiscrByTrackCountingHighPur;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur;   //!
   TBranch        *b_PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff;   //!
   TBranch        *b_PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA;   //!
   TBranch        *b_PFPatJet_BJetDiscrByjetProbabilityBJetTags;   //!
   TBranch        *b_PFPatJet_BJetDiscrByjetBProbabilityBJetTags;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySoftElectronBJetTags;   //!
   TBranch        *b_PFPatJet_BJetDiscrBySoftMuonBJetTags;   //!
   TBranch        *b_hasMatchedGenParton_ToPFJet_;   //!
   TBranch        *b_hasMatchedGenJet_ToPFJet_;   //!
   TBranch        *b_MatchedGenPartonToPFJet_pdgId;   //!
   TBranch        *b_MatchedGenPartonToPFJet_status;   //!
   TBranch        *b_MatchedGenJetToPFJet_E;   //!
   TBranch        *b_MatchedGenJetToPFJet_pt;   //!
   TBranch        *b_MatchedGenJetToPFJet_eta;   //!
   TBranch        *b_MatchedGenJetToPFJet_phi;   //!
   TBranch        *b_MatchedGenJetToPFJet_nDaughters;   //!
   TBranch        *b_PFPatJet_JetPartonFlavor;   //!
   TBranch        *b_PFRecoJet_n;   //!
   TBranch        *b_PFRecoJet_E;   //!
   TBranch        *b_PFRecoJet_et;   //!
   TBranch        *b_PFRecoJet_pt;   //!
   TBranch        *b_PFRecoJet_eta;   //!
   TBranch        *b_PFRecoJet_phi;   //!
   TBranch        *b_PFRecoJet_theta;   //!
   TBranch        *b_PFRecoJet_px;   //!
   TBranch        *b_PFRecoJet_py;   //!
   TBranch        *b_PFRecoJet_pz;   //!
   TBranch        *b_PFRecoJet_vx;   //!
   TBranch        *b_PFRecoJet_vy;   //!
   TBranch        *b_PFRecoJet_vz;   //!
   TBranch        *b_PFRecoJet_ChargedEmEnergy;   //!
   TBranch        *b_PFRecoJet_ChargedEmEnergyFrac;   //!
   TBranch        *b_PFRecoJet_ChargedHadEnergy;   //!
   TBranch        *b_PFRecoJet_ChargedHadEnergyFrac;   //!
   TBranch        *b_PFRecoJet_ChargedHadMult;   //!
   TBranch        *b_PFRecoJet_ChargedMult;   //!
   TBranch        *b_PFRecoJet_NConstituents;   //!
   TBranch        *b_PFRecoJet_HFEMEnergy;   //!
   TBranch        *b_PFRecoJet_HFEMEnergyFrac;   //!
   TBranch        *b_PFRecoJet_HFEMMult;   //!
   TBranch        *b_PFRecoJet_HFHadEnergy;   //!
   TBranch        *b_PFRecoJet_HFHadEnergyFrac;   //!
   TBranch        *b_PFRecoJet_HFHadMult;   //!
   TBranch        *b_PFRecoJet_NeutralEmEnergy;   //!
   TBranch        *b_PFRecoJet_NeutralEmEnergyFrac;   //!
   TBranch        *b_PFRecoJet_NeutralHadEnergy;   //!
   TBranch        *b_PFRecoJet_NeutralHadEnergyFrac;   //!
   TBranch        *b_PFRecoJet_NeutralHadMult;   //!
   TBranch        *b_PFRecoJet_NeutralMult;   //!
   TBranch        *b_PFRecoJet_jecUncertainity;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByTrackCountingHighEff;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByTrackCountingHighPur;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByjetProbabilityBJetTags;   //!
   TBranch        *b_PFRecoJet_BJetDiscrByjetBProbabilityBJetTags;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySoftElectronBJetTags;   //!
   TBranch        *b_PFRecoJet_BJetDiscrBySoftMuonBJetTags;   //!
   TBranch        *b_GenJet_n;   //!
   TBranch        *b_GenJet_E;   //!
   TBranch        *b_GenJet_et;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_px;   //!
   TBranch        *b_GenJet_py;   //!
   TBranch        *b_GenJet_pz;   //!
   TBranch        *b_GenJet_vx;   //!
   TBranch        *b_GenJet_vy;   //!
   TBranch        *b_GenJet_vz;   //!
   TBranch        *b_GenJet_nDaughters;   //!
   TBranch        *b_GenJet_daughters_status;   //!
   TBranch        *b_GenJet_daughters_pdgid;   //!
   TBranch        *b_GenJet_daughters_E;   //!
   TBranch        *b_GenJet_daughters_pt;   //!
   TBranch        *b_GenJet_daughters_eta;   //!
   TBranch        *b_GenJet_daughters_phi;   //!
   TBranch        *b_Vertex_n;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_nchi2;   //!
   TBranch        *b_Vertex_ndof;   //!
   TBranch        *b_Vertex_tracksSize;   //!
   TBranch        *b_Vertex_isFake;   //!
   TBranch        *b_Vertex_isValid;   //!
   TBranch        *b_Vertex_d0;   //!
   TBranch        *b_Tracks_n;   //!
   TBranch        *b_Track_pt;   //!
   TBranch        *b_Track_px;   //!
   TBranch        *b_Track_py;   //!
   TBranch        *b_Track_pz;   //!
   TBranch        *b_Track_vx;   //!
   TBranch        *b_Track_vy;   //!
   TBranch        *b_Track_vz;   //!
   TBranch        *b_Track_eta;   //!
   TBranch        *b_Track_phi;   //!
   TBranch        *b_Track_theta;   //!
   TBranch        *b_Track_chi2;   //!
   TBranch        *b_IsScrapingEvent_;   //!
   TBranch        *b_Scraping_FractionOfGoodTracks;   //!
   TBranch        *b_HLT_Photon_triggers;   //!
   TBranch        *b_HLT_Photon_trig_prescales;   //!
   TBranch        *b_HLT_Photon_ifTriggerPassed;   //!
   TBranch        *b_HLT_Photon_nTriggers;   //!
   TBranch        *b_HLT_Photon_triggerIndex;   //!
   TBranch        *b_HLT_Photon_nFilters;   //!
   TBranch        *b_HLT_Photon_FilterNames;   //!
   TBranch        *b_HLT_Photon_trigger_FilterStartPosition;   //!
   TBranch        *b_HLT_Photon_trigger_FilterEndPosition;   //!
   TBranch        *b_HLT_Photon_FilterObjects_pt;   //!
   TBranch        *b_HLT_Photon_FilterObjects_eta;   //!
   TBranch        *b_HLT_Photon_FilterObjects_phi;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_sigma;   //!
   TBranch        *b_rho25;   //!
   TBranch        *b_sigma25;   //!

   MyEvent_MC(TTree *tree=0);
   virtual ~MyEvent_MC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyEvent_MC_cxx
MyEvent_MC::MyEvent_MC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/G_Pt_120to170/AODSIM_G_Pt_120to170_110_1_dHz.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/uscms/store/user/lpcqstar/NTuples/MC/PhotonJet/G_Pt_120to170/AODSIM_G_Pt_120to170_110_1_dHz.root");
      }
      f->GetObject("myEvent",tree);

   }
   Init(tree);
}

MyEvent_MC::~MyEvent_MC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyEvent_MC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyEvent_MC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyEvent_MC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Photon_E = 0;
   Photon_et = 0;
   Photon_pt = 0;
   Photon_eta = 0;
   Photon_phi = 0;
   Photon_theta = 0;
   Photon_px = 0;
   Photon_py = 0;
   Photon_pz = 0;
   Photon_vx = 0;
   Photon_vy = 0;
   Photon_vz = 0;
   Photon_r9 = 0;
   Photon_maxEnergyXtal = 0;
   Photon_e1x5 = 0;
   Photon_e2x5 = 0;
   Photon_e3x3 = 0;
   Photon_e5x5 = 0;
   Photon_r1x5 = 0;
   Photon_r2x5 = 0;
   Photon_SigmaEtaEta = 0;
   Photon_SigmaIEtaIEta = 0;
   Photon_SigmaEtaPhi = 0;
   Photon_SigmaIEtaIPhi = 0;
   Photon_SigmaPhiPhi = 0;
   Photon_SigmaIPhiIPhi = 0;
   Photon_roundness = 0;
   Photon_angle = 0;
   Photon_swissCross = 0;
   Photon_s9 = 0;
   Photon_e4Overe1 = 0;
   Photon_e6Overe2 = 0;
   Photon_e2Overe9 = 0;
   Photon_rookFraction = 0;
   Photon_isEB = 0;
   Photon_isEE = 0;
   Photon_isEBGap = 0;
   Photon_isEEGap = 0;
   Photon_isEBEEGap = 0;
   Photon_ecalRecHitSumEtConeDR03 = 0;
   Photon_hcalTowerSumEtConeDR03 = 0;
   Photon_hcalDepth1TowerSumEtConeDR03 = 0;
   Photon_hcalDepth2TowerSumEtConeDR03 = 0;
   Photon_trkSumPtSolidConeDR03 = 0;
   Photon_trkSumPtHollowConeDR03 = 0;
   Photon_nTrkSolidConeDR03 = 0;
   Photon_nTrkHollowConeDR03 = 0;
   Photon_ecalRecHitSumEtConeDR04 = 0;
   Photon_hcalTowerSumEtConeDR04 = 0;
   Photon_hcalDepth1TowerSumEtConeDR04 = 0;
   Photon_hcalDepth2TowerSumEtConeDR04 = 0;
   Photon_trkSumPtSolidConeDR04 = 0;
   Photon_trkSumPtHollowConeDR04 = 0;
   Photon_nTrkSolidConeDR04 = 0;
   Photon_nTrkHollowConeDR04 = 0;
   Photon_HoE = 0;
   Photon_SingleTowerHoE = 0;
   Photon_hasConvTrk = 0;
   Photon_hasPixelSeed = 0;
   passedConvSafeElectronVeto = 0;
   Photon_SC_nOfBasicClusters = 0;
   Photon_SC_rawEnergy = 0;
   Photon_SC_preShowerEnergy = 0;
   Photon_SC_energy = 0;
   Photon_SC_eta = 0;
   Photon_SC_phi = 0;
   Photon_SC_x = 0;
   Photon_SC_y = 0;
   Photon_SC_z = 0;
   Photon_SC_etaWidth = 0;
   Photon_SC_phiWidth = 0;
   Photon_mipChi2 = 0;
   Photon_mipTotEnergy = 0;
   Photon_mipSlope = 0;
   Photon_mipIntercept = 0;
   Photon_mipNhitCone = 0;
   Photon_mipIsHalo = 0;
   Photon_nConvTracks = 0;
   Photon_isConverted = 0;
   Photon_pairInvariantMass = 0;
   Photon_pairCotThetaSeparation = 0;
   Photon_pairMomentum_x = 0;
   Photon_pairMomentum_y = 0;
   Photon_pairMomentum_z = 0;
   Photon_EoverP = 0;
   Photon_conv_vx = 0;
   Photon_conv_vy = 0;
   Photon_conv_vz = 0;
   Photon_zOfPrimaryVtxFromTrks = 0;
   Photon_distOfMinimumApproach = 0;
   Photon_dPhiTracksAtVtx = 0;
   Photon_dPhiTracksAtEcal = 0;
   Photon_dEtaTracksAtEcal = 0;
   Photon_nCrystals = 0;
   Photon_xtal_timing = 0;
   Photon_xtal_timeErr = 0;
   Photon_avgTimeAllxtals = 0;
   Photon_xtal_energy = 0;
   Photon_xtal_EBieta = 0;
   Photon_xtal_EBiphi = 0;
   Photon_xtal_EBrecoFlag = 0;
   PFIsoPhoton03 = 0;
   PFIsoNeutral03 = 0;
   PFIsoCharged03 = 0;
   PFIsoSum03 = 0;
   PFIsoChargedWorstvtx03 = 0;
   GenPhoton_status = 0;
   GenPhoton_E = 0;
   GenPhoton_et = 0;
   GenPhoton_pt = 0;
   GenPhoton_eta = 0;
   GenPhoton_phi = 0;
   GenPhoton_px = 0;
   GenPhoton_py = 0;
   GenPhoton_pz = 0;
   GenPhoton_vx = 0;
   GenPhoton_vy = 0;
   GenPhoton_vz = 0;
   GenPhoton_mother_status = 0;
   GenPhoton_mother_pdgid = 0;
   GenPhoton_mother_E = 0;
   GenPhoton_mother_et = 0;
   GenPhoton_mother_pt = 0;
   GenPhoton_mother_eta = 0;
   GenPhoton_mother_phi = 0;
   GenPhoton_mother_px = 0;
   GenPhoton_mother_py = 0;
   GenPhoton_mother_pz = 0;
   GenPhoton_mother_vx = 0;
   GenPhoton_mother_vy = 0;
   GenPhoton_mother_vz = 0;
   GenPhoton_grandmother_status = 0;
   GenPhoton_grandmother_pdgid = 0;
   GenPhoton_grandmother_E = 0;
   GenPhoton_grandmother_et = 0;
   GenPhoton_grandmother_pt = 0;
   GenPhoton_grandmother_eta = 0;
   GenPhoton_grandmother_phi = 0;
   GenPhoton_grandmother_px = 0;
   GenPhoton_grandmother_py = 0;
   GenPhoton_grandmother_pz = 0;
   GenPhoton_grandmother_vx = 0;
   GenPhoton_grandmother_vy = 0;
   GenPhoton_grandmother_vz = 0;
   hasMatchedGenPhoton_ = 0;
   MatchedGenPhoton_status = 0;
   MatchedGenPhoton_pdgid = 0;
   MatchedGenPhoton_E = 0;
   MatchedGenPhoton_et = 0;
   MatchedGenPhoton_pt = 0;
   MatchedGenPhoton_eta = 0;
   MatchedGenPhoton_phi = 0;
   MatchedGenPhoton_px = 0;
   MatchedGenPhoton_py = 0;
   MatchedGenPhoton_pz = 0;
   MatchedGenPhoton_vx = 0;
   MatchedGenPhoton_vy = 0;
   MatchedGenPhoton_vz = 0;
   MatchedGenPhoton_nMothers = 0;
   MatchedGenPhoton_mothersPdgid = 0;
   MatchedGenPhoton_mothersStatus = 0;
   PFPatJet_E = 0;
   PFPatJet_et = 0;
   PFPatJet_pt = 0;
   PFPatJet_eta = 0;
   PFPatJet_phi = 0;
   PFPatJet_theta = 0;
   PFPatJet_px = 0;
   PFPatJet_py = 0;
   PFPatJet_pz = 0;
   PFPatJet_vx = 0;
   PFPatJet_vy = 0;
   PFPatJet_vz = 0;
   PFPatJet_ChargedEmEnergy = 0;
   PFPatJet_ChargedEmEnergyFrac = 0;
   PFPatJet_ChargedHadEnergy = 0;
   PFPatJet_ChargedHadEnergyFrac = 0;
   PFPatJet_ChargedHadMult = 0;
   PFPatJet_ChargedMult = 0;
   PFPatJet_NConstituents = 0;
   PFPatJet_HFEMEnergy = 0;
   PFPatJet_HFEMEnergyFrac = 0;
   PFPatJet_HFEMMult = 0;
   PFPatJet_HFHadEnergy = 0;
   PFPatJet_HFHadEnergyFrac = 0;
   PFPatJet_HFHadMult = 0;
   PFPatJet_NeutralEmEnergy = 0;
   PFPatJet_NeutralEmEnergyFrac = 0;
   PFPatJet_NeutralHadEnergy = 0;
   PFPatJet_NeutralHadEnergyFrac = 0;
   PFPatJet_NeutralHadMult = 0;
   PFPatJet_NeutralMult = 0;
   PFPatJet_jecUncertainity = 0;
   PFPatJet_puJetIdCutBased_MVA = 0;
   PFPatJet_puJetIdSimple_MVA = 0;
   PFPatJet_puJetIdFull_MVA = 0;
   PFPatJet_PassPUJetIdCutBased_loose = 0;
   PFPatJet_PassPUJetIdCutBased_medium = 0;
   PFPatJet_PassPUJetIdCutBased_tight = 0;
   PFPatJet_PassPUJetIdSimple_loose = 0;
   PFPatJet_PassPUJetIdSimple_medium = 0;
   PFPatJet_PassPUJetIdSimple_tight = 0;
   PFPatJet_PassPUJetIdFull_loose = 0;
   PFPatJet_PassPUJetIdFull_medium = 0;
   PFPatJet_PassPUJetIdFull_tight = 0;
   PFPatJet_BJetDiscrByTrackCountingHighEff = 0;
   PFPatJet_BJetDiscrByTrackCountingHighPur = 0;
   PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff = 0;
   PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur = 0;
   PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff = 0;
   PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA = 0;
   PFPatJet_BJetDiscrByjetProbabilityBJetTags = 0;
   PFPatJet_BJetDiscrByjetBProbabilityBJetTags = 0;
   PFPatJet_BJetDiscrBySoftElectronBJetTags = 0;
   PFPatJet_BJetDiscrBySoftMuonBJetTags = 0;
   hasMatchedGenParton_ToPFJet_ = 0;
   hasMatchedGenJet_ToPFJet_ = 0;
   MatchedGenPartonToPFJet_pdgId = 0;
   MatchedGenPartonToPFJet_status = 0;
   MatchedGenJetToPFJet_E = 0;
   MatchedGenJetToPFJet_pt = 0;
   MatchedGenJetToPFJet_eta = 0;
   MatchedGenJetToPFJet_phi = 0;
   MatchedGenJetToPFJet_nDaughters = 0;
   PFPatJet_JetPartonFlavor = 0;
   PFRecoJet_E = 0;
   PFRecoJet_et = 0;
   PFRecoJet_pt = 0;
   PFRecoJet_eta = 0;
   PFRecoJet_phi = 0;
   PFRecoJet_theta = 0;
   PFRecoJet_px = 0;
   PFRecoJet_py = 0;
   PFRecoJet_pz = 0;
   PFRecoJet_vx = 0;
   PFRecoJet_vy = 0;
   PFRecoJet_vz = 0;
   PFRecoJet_ChargedEmEnergy = 0;
   PFRecoJet_ChargedEmEnergyFrac = 0;
   PFRecoJet_ChargedHadEnergy = 0;
   PFRecoJet_ChargedHadEnergyFrac = 0;
   PFRecoJet_ChargedHadMult = 0;
   PFRecoJet_ChargedMult = 0;
   PFRecoJet_NConstituents = 0;
   PFRecoJet_HFEMEnergy = 0;
   PFRecoJet_HFEMEnergyFrac = 0;
   PFRecoJet_HFEMMult = 0;
   PFRecoJet_HFHadEnergy = 0;
   PFRecoJet_HFHadEnergyFrac = 0;
   PFRecoJet_HFHadMult = 0;
   PFRecoJet_NeutralEmEnergy = 0;
   PFRecoJet_NeutralEmEnergyFrac = 0;
   PFRecoJet_NeutralHadEnergy = 0;
   PFRecoJet_NeutralHadEnergyFrac = 0;
   PFRecoJet_NeutralHadMult = 0;
   PFRecoJet_NeutralMult = 0;
   PFRecoJet_jecUncertainity = 0;
   PFRecoJet_BJetDiscrByTrackCountingHighEff = 0;
   PFRecoJet_BJetDiscrByTrackCountingHighPur = 0;
   PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff = 0;
   PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur = 0;
   PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff = 0;
   PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA = 0;
   PFRecoJet_BJetDiscrByjetProbabilityBJetTags = 0;
   PFRecoJet_BJetDiscrByjetBProbabilityBJetTags = 0;
   PFRecoJet_BJetDiscrBySoftElectronBJetTags = 0;
   PFRecoJet_BJetDiscrBySoftMuonBJetTags = 0;
   GenJet_E = 0;
   GenJet_et = 0;
   GenJet_pt = 0;
   GenJet_eta = 0;
   GenJet_phi = 0;
   GenJet_px = 0;
   GenJet_py = 0;
   GenJet_pz = 0;
   GenJet_vx = 0;
   GenJet_vy = 0;
   GenJet_vz = 0;
   GenJet_nDaughters = 0;
   GenJet_daughters_status = 0;
   GenJet_daughters_pdgid = 0;
   GenJet_daughters_E = 0;
   GenJet_daughters_pt = 0;
   GenJet_daughters_eta = 0;
   GenJet_daughters_phi = 0;
   Vertex_x = 0;
   Vertex_y = 0;
   Vertex_z = 0;
   Vertex_chi2 = 0;
   Vertex_nchi2 = 0;
   Vertex_ndof = 0;
   Vertex_tracksSize = 0;
   Vertex_isFake = 0;
   Vertex_isValid = 0;
   Vertex_d0 = 0;
   Track_pt = 0;
   Track_px = 0;
   Track_py = 0;
   Track_pz = 0;
   Track_vx = 0;
   Track_vy = 0;
   Track_vz = 0;
   Track_eta = 0;
   Track_phi = 0;
   Track_theta = 0;
   Track_chi2 = 0;
   HLT_Photon_triggers = 0;
   HLT_Photon_trig_prescales = 0;
   HLT_Photon_ifTriggerPassed = 0;
   HLT_Photon_triggerIndex = 0;
   HLT_Photon_nFilters = 0;
   HLT_Photon_FilterNames = 0;
   HLT_Photon_trigger_FilterStartPosition = 0;
   HLT_Photon_trigger_FilterEndPosition = 0;
   HLT_Photon_FilterObjects_pt = 0;
   HLT_Photon_FilterObjects_eta = 0;
   HLT_Photon_FilterObjects_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
   fChain->SetBranchAddress("BXNumber", &BXNumber, &b_BXNumber);
   fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1, &b_totalIntensityBeam1);
   fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2, &b_totalIntensityBeam2);
   fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi, &b_avgInsDelLumi);
   fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr, &b_avgInsDelLumiErr);
   fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_avgInsRecLumi);
   fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr, &b_avgInsRecLumiErr);
   fChain->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
   fChain->SetBranchAddress("npuVerticesm1", &npuVerticesm1, &b_npuVerticesm1);
   fChain->SetBranchAddress("npuVerticesp1", &npuVerticesp1, &b_npuVerticesp1);
   fChain->SetBranchAddress("npuVerticespm2", &npuVerticespm2, &b_npuVerticespm2);
   fChain->SetBranchAddress("trueNumofInteractions", &trueNumofInteractions, &b_trueNumofInteractions);
   fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
   fChain->SetBranchAddress("Photon_E", &Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon_et", &Photon_et, &b_Photon_et);
   fChain->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_theta", &Photon_theta, &b_Photon_theta);
   fChain->SetBranchAddress("Photon_px", &Photon_px, &b_Photon_px);
   fChain->SetBranchAddress("Photon_py", &Photon_py, &b_Photon_py);
   fChain->SetBranchAddress("Photon_pz", &Photon_pz, &b_Photon_pz);
   fChain->SetBranchAddress("Photon_vx", &Photon_vx, &b_Photon_vx);
   fChain->SetBranchAddress("Photon_vy", &Photon_vy, &b_Photon_vy);
   fChain->SetBranchAddress("Photon_vz", &Photon_vz, &b_Photon_vz);
   fChain->SetBranchAddress("Photon_r9", &Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_maxEnergyXtal", &Photon_maxEnergyXtal, &b_Photon_maxEnergyXtal);
   fChain->SetBranchAddress("Photon_e1x5", &Photon_e1x5, &b_Photon_e1x5);
   fChain->SetBranchAddress("Photon_e2x5", &Photon_e2x5, &b_Photon_e2x5);
   fChain->SetBranchAddress("Photon_e3x3", &Photon_e3x3, &b_Photon_e3x3);
   fChain->SetBranchAddress("Photon_e5x5", &Photon_e5x5, &b_Photon_e5x5);
   fChain->SetBranchAddress("Photon_r1x5", &Photon_r1x5, &b_Photon_r1x5);
   fChain->SetBranchAddress("Photon_r2x5", &Photon_r2x5, &b_Photon_r2x5);
   fChain->SetBranchAddress("Photon_SigmaEtaEta", &Photon_SigmaEtaEta, &b_Photon_SigmaEtaEta);
   fChain->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
   fChain->SetBranchAddress("Photon_SigmaEtaPhi", &Photon_SigmaEtaPhi, &b_Photon_SigmaEtaPhi);
   fChain->SetBranchAddress("Photon_SigmaIEtaIPhi", &Photon_SigmaIEtaIPhi, &b_Photon_SigmaIEtaIPhi);
   fChain->SetBranchAddress("Photon_SigmaPhiPhi", &Photon_SigmaPhiPhi, &b_Photon_SigmaPhiPhi);
   fChain->SetBranchAddress("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi, &b_Photon_SigmaIPhiIPhi);
   fChain->SetBranchAddress("Photon_roundness", &Photon_roundness, &b_Photon_roundness);
   fChain->SetBranchAddress("Photon_angle", &Photon_angle, &b_Photon_angle);
   fChain->SetBranchAddress("Photon_swissCross", &Photon_swissCross, &b_Photon_swissCross);
   fChain->SetBranchAddress("Photon_s9", &Photon_s9, &b_Photon_s9);
   fChain->SetBranchAddress("Photon_e4Overe1", &Photon_e4Overe1, &b_Photon_e4Overe1);
   fChain->SetBranchAddress("Photon_e6Overe2", &Photon_e6Overe2, &b_Photon_e6Overe2);
   fChain->SetBranchAddress("Photon_e2Overe9", &Photon_e2Overe9, &b_Photon_e2Overe9);
   fChain->SetBranchAddress("Photon_rookFraction", &Photon_rookFraction, &b_Photon_rookFraction);
   fChain->SetBranchAddress("Photon_isEB", &Photon_isEB, &b_Photon_isEB);
   fChain->SetBranchAddress("Photon_isEE", &Photon_isEE, &b_Photon_isEE);
   fChain->SetBranchAddress("Photon_isEBGap", &Photon_isEBGap, &b_Photon_isEBGap);
   fChain->SetBranchAddress("Photon_isEEGap", &Photon_isEEGap, &b_Photon_isEEGap);
   fChain->SetBranchAddress("Photon_isEBEEGap", &Photon_isEBEEGap, &b_Photon_isEBEEGap);
   fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", &Photon_ecalRecHitSumEtConeDR03, &b_Photon_ecalRecHitSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", &Photon_hcalTowerSumEtConeDR03, &b_Photon_hcalTowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", &Photon_hcalDepth1TowerSumEtConeDR03, &b_Photon_hcalDepth1TowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", &Photon_hcalDepth2TowerSumEtConeDR03, &b_Photon_hcalDepth2TowerSumEtConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", &Photon_trkSumPtSolidConeDR03, &b_Photon_trkSumPtSolidConeDR03);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", &Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", &Photon_nTrkSolidConeDR03, &b_Photon_nTrkSolidConeDR03);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", &Photon_nTrkHollowConeDR03, &b_Photon_nTrkHollowConeDR03);
   fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", &Photon_ecalRecHitSumEtConeDR04, &b_Photon_ecalRecHitSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", &Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", &Photon_hcalDepth1TowerSumEtConeDR04, &b_Photon_hcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", &Photon_hcalDepth2TowerSumEtConeDR04, &b_Photon_hcalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", &Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
   fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", &Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
   fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", &Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
   fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", &Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
   fChain->SetBranchAddress("Photon_HoE", &Photon_HoE, &b_Photon_HoE);
   fChain->SetBranchAddress("Photon_SingleTowerHoE", &Photon_SingleTowerHoE, &b_Photon_SingleTowerHoE);
   fChain->SetBranchAddress("Photon_hasConvTrk", &Photon_hasConvTrk, &b_Photon_hasConvTrk);
   fChain->SetBranchAddress("Photon_hasPixelSeed", &Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   fChain->SetBranchAddress("passedConvSafeElectronVeto", &passedConvSafeElectronVeto, &b_passedConvSafeElectronVeto);
   fChain->SetBranchAddress("Photon_SC_nOfBasicClusters", &Photon_SC_nOfBasicClusters, &b_Photon_SC_nOfBasicClusters);
   fChain->SetBranchAddress("Photon_SC_rawEnergy", &Photon_SC_rawEnergy, &b_Photon_SC_rawEnergy);
   fChain->SetBranchAddress("Photon_SC_preShowerEnergy", &Photon_SC_preShowerEnergy, &b_Photon_SC_preShowerEnergy);
   fChain->SetBranchAddress("Photon_SC_energy", &Photon_SC_energy, &b_Photon_SC_energy);
   fChain->SetBranchAddress("Photon_SC_eta", &Photon_SC_eta, &b_Photon_SC_eta);
   fChain->SetBranchAddress("Photon_SC_phi", &Photon_SC_phi, &b_Photon_SC_phi);
   fChain->SetBranchAddress("Photon_SC_x", &Photon_SC_x, &b_Photon_SC_x);
   fChain->SetBranchAddress("Photon_SC_y", &Photon_SC_y, &b_Photon_SC_y);
   fChain->SetBranchAddress("Photon_SC_z", &Photon_SC_z, &b_Photon_SC_z);
   fChain->SetBranchAddress("Photon_SC_etaWidth", &Photon_SC_etaWidth, &b_Photon_SC_etaWidth);
   fChain->SetBranchAddress("Photon_SC_phiWidth", &Photon_SC_phiWidth, &b_Photon_SC_phiWidth);
   fChain->SetBranchAddress("Photon_mipChi2", &Photon_mipChi2, &b_Photon_mipChi2);
   fChain->SetBranchAddress("Photon_mipTotEnergy", &Photon_mipTotEnergy, &b_Photon_mipTotEnergy);
   fChain->SetBranchAddress("Photon_mipSlope", &Photon_mipSlope, &b_Photon_mipSlope);
   fChain->SetBranchAddress("Photon_mipIntercept", &Photon_mipIntercept, &b_Photon_mipIntercept);
   fChain->SetBranchAddress("Photon_mipNhitCone", &Photon_mipNhitCone, &b_Photon_mipNhitCone);
   fChain->SetBranchAddress("Photon_mipIsHalo", &Photon_mipIsHalo, &b_Photon_mipIsHalo);
   fChain->SetBranchAddress("Photon_nConvTracks", &Photon_nConvTracks, &b_Photon_nConvTracks);
   fChain->SetBranchAddress("Photon_isConverted", &Photon_isConverted, &b_Photon_isConverted);
   fChain->SetBranchAddress("Photon_pairInvariantMass", &Photon_pairInvariantMass, &b_Photon_pairInvariantMass);
   fChain->SetBranchAddress("Photon_pairCotThetaSeparation", &Photon_pairCotThetaSeparation, &b_Photon_pairCotThetaSeparation);
   fChain->SetBranchAddress("Photon_pairMomentum_x", &Photon_pairMomentum_x, &b_Photon_pairMomentum_x);
   fChain->SetBranchAddress("Photon_pairMomentum_y", &Photon_pairMomentum_y, &b_Photon_pairMomentum_y);
   fChain->SetBranchAddress("Photon_pairMomentum_z", &Photon_pairMomentum_z, &b_Photon_pairMomentum_z);
   fChain->SetBranchAddress("Photon_EoverP", &Photon_EoverP, &b_Photon_EoverP);
   fChain->SetBranchAddress("Photon_conv_vx", &Photon_conv_vx, &b_Photon_conv_vx);
   fChain->SetBranchAddress("Photon_conv_vy", &Photon_conv_vy, &b_Photon_conv_vy);
   fChain->SetBranchAddress("Photon_conv_vz", &Photon_conv_vz, &b_Photon_conv_vz);
   fChain->SetBranchAddress("Photon_zOfPrimaryVtxFromTrks", &Photon_zOfPrimaryVtxFromTrks, &b_Photon_zOfPrimaryVtxFromTrks);
   fChain->SetBranchAddress("Photon_distOfMinimumApproach", &Photon_distOfMinimumApproach, &b_Photon_distOfMinimumApproach);
   fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", &Photon_dPhiTracksAtVtx, &b_Photon_dPhiTracksAtVtx);
   fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", &Photon_dPhiTracksAtEcal, &b_Photon_dPhiTracksAtEcal);
   fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", &Photon_dEtaTracksAtEcal, &b_Photon_dEtaTracksAtEcal);
   fChain->SetBranchAddress("Photon_nCrystals", &Photon_nCrystals, &b_Photon_nCrystals);
   fChain->SetBranchAddress("Photon_xtal_timing", &Photon_xtal_timing, &b_Photon_xtal_timing);
   fChain->SetBranchAddress("Photon_xtal_timeErr", &Photon_xtal_timeErr, &b_Photon_xtal_timeErr);
   fChain->SetBranchAddress("Photon_avgTimeAllxtals", &Photon_avgTimeAllxtals, &b_Photon_avgTimeAllxtals);
   fChain->SetBranchAddress("Photon_xtal_energy", &Photon_xtal_energy, &b_Photon_xtal_energy);
   fChain->SetBranchAddress("Photon_xtal_EBieta", &Photon_xtal_EBieta, &b_Photon_xtal_EBieta);
   fChain->SetBranchAddress("Photon_xtal_EBiphi", &Photon_xtal_EBiphi, &b_Photon_xtal_EBiphi);
   fChain->SetBranchAddress("Photon_xtal_EBrecoFlag", &Photon_xtal_EBrecoFlag, &b_Photon_xtal_EBrecoFlag);
   fChain->SetBranchAddress("PFIsoPhoton03", &PFIsoPhoton03, &b_PFIsoPhoton03);
   fChain->SetBranchAddress("PFIsoNeutral03", &PFIsoNeutral03, &b_PFIsoNeutral03);
   fChain->SetBranchAddress("PFIsoCharged03", &PFIsoCharged03, &b_PFIsoCharged03);
   fChain->SetBranchAddress("PFIsoSum03", &PFIsoSum03, &b_PFIsoSum03);
   fChain->SetBranchAddress("PFIsoChargedWorstvtx03", &PFIsoChargedWorstvtx03, &b_PFIsoChargedWorstvtx03);
   fChain->SetBranchAddress("GenPhoton_n", &GenPhoton_n, &b_GenPhoton_n);
   fChain->SetBranchAddress("GenPhoton_status", &GenPhoton_status, &b_GenPhoton_status);
   fChain->SetBranchAddress("GenPhoton_E", &GenPhoton_E, &b_GenPhoton_E);
   fChain->SetBranchAddress("GenPhoton_et", &GenPhoton_et, &b_GenPhoton_et);
   fChain->SetBranchAddress("GenPhoton_pt", &GenPhoton_pt, &b_GenPhoton_pt);
   fChain->SetBranchAddress("GenPhoton_eta", &GenPhoton_eta, &b_GenPhoton_eta);
   fChain->SetBranchAddress("GenPhoton_phi", &GenPhoton_phi, &b_GenPhoton_phi);
   fChain->SetBranchAddress("GenPhoton_px", &GenPhoton_px, &b_GenPhoton_px);
   fChain->SetBranchAddress("GenPhoton_py", &GenPhoton_py, &b_GenPhoton_py);
   fChain->SetBranchAddress("GenPhoton_pz", &GenPhoton_pz, &b_GenPhoton_pz);
   fChain->SetBranchAddress("GenPhoton_vx", &GenPhoton_vx, &b_GenPhoton_vx);
   fChain->SetBranchAddress("GenPhoton_vy", &GenPhoton_vy, &b_GenPhoton_vy);
   fChain->SetBranchAddress("GenPhoton_vz", &GenPhoton_vz, &b_GenPhoton_vz);
   fChain->SetBranchAddress("GenPhoton_mother_status", &GenPhoton_mother_status, &b_GenPhoton_mother_status);
   fChain->SetBranchAddress("GenPhoton_mother_pdgid", &GenPhoton_mother_pdgid, &b_GenPhoton_mother_pdgid);
   fChain->SetBranchAddress("GenPhoton_mother_E", &GenPhoton_mother_E, &b_GenPhoton_mother_E);
   fChain->SetBranchAddress("GenPhoton_mother_et", &GenPhoton_mother_et, &b_GenPhoton_mother_et);
   fChain->SetBranchAddress("GenPhoton_mother_pt", &GenPhoton_mother_pt, &b_GenPhoton_mother_pt);
   fChain->SetBranchAddress("GenPhoton_mother_eta", &GenPhoton_mother_eta, &b_GenPhoton_mother_eta);
   fChain->SetBranchAddress("GenPhoton_mother_phi", &GenPhoton_mother_phi, &b_GenPhoton_mother_phi);
   fChain->SetBranchAddress("GenPhoton_mother_px", &GenPhoton_mother_px, &b_GenPhoton_mother_px);
   fChain->SetBranchAddress("GenPhoton_mother_py", &GenPhoton_mother_py, &b_GenPhoton_mother_py);
   fChain->SetBranchAddress("GenPhoton_mother_pz", &GenPhoton_mother_pz, &b_GenPhoton_mother_pz);
   fChain->SetBranchAddress("GenPhoton_mother_vx", &GenPhoton_mother_vx, &b_GenPhoton_mother_vx);
   fChain->SetBranchAddress("GenPhoton_mother_vy", &GenPhoton_mother_vy, &b_GenPhoton_mother_vy);
   fChain->SetBranchAddress("GenPhoton_mother_vz", &GenPhoton_mother_vz, &b_GenPhoton_mother_vz);
   fChain->SetBranchAddress("GenPhoton_grandmother_status", &GenPhoton_grandmother_status, &b_GenPhoton_grandmother_status);
   fChain->SetBranchAddress("GenPhoton_grandmother_pdgid", &GenPhoton_grandmother_pdgid, &b_GenPhoton_grandmother_pdgid);
   fChain->SetBranchAddress("GenPhoton_grandmother_E", &GenPhoton_grandmother_E, &b_GenPhoton_grandmother_E);
   fChain->SetBranchAddress("GenPhoton_grandmother_et", &GenPhoton_grandmother_et, &b_GenPhoton_grandmother_et);
   fChain->SetBranchAddress("GenPhoton_grandmother_pt", &GenPhoton_grandmother_pt, &b_GenPhoton_grandmother_pt);
   fChain->SetBranchAddress("GenPhoton_grandmother_eta", &GenPhoton_grandmother_eta, &b_GenPhoton_grandmother_eta);
   fChain->SetBranchAddress("GenPhoton_grandmother_phi", &GenPhoton_grandmother_phi, &b_GenPhoton_grandmother_phi);
   fChain->SetBranchAddress("GenPhoton_grandmother_px", &GenPhoton_grandmother_px, &b_GenPhoton_grandmother_px);
   fChain->SetBranchAddress("GenPhoton_grandmother_py", &GenPhoton_grandmother_py, &b_GenPhoton_grandmother_py);
   fChain->SetBranchAddress("GenPhoton_grandmother_pz", &GenPhoton_grandmother_pz, &b_GenPhoton_grandmother_pz);
   fChain->SetBranchAddress("GenPhoton_grandmother_vx", &GenPhoton_grandmother_vx, &b_GenPhoton_grandmother_vx);
   fChain->SetBranchAddress("GenPhoton_grandmother_vy", &GenPhoton_grandmother_vy, &b_GenPhoton_grandmother_vy);
   fChain->SetBranchAddress("GenPhoton_grandmother_vz", &GenPhoton_grandmother_vz, &b_GenPhoton_grandmother_vz);
   fChain->SetBranchAddress("hasMatchedGenPhoton_", &hasMatchedGenPhoton_, &b_hasMatchedGenPhoton_);
   fChain->SetBranchAddress("MatchedGenPhoton_status", &MatchedGenPhoton_status, &b_MatchedGenPhoton_status);
   fChain->SetBranchAddress("MatchedGenPhoton_pdgid", &MatchedGenPhoton_pdgid, &b_MatchedGenPhoton_pdgid);
   fChain->SetBranchAddress("MatchedGenPhoton_E", &MatchedGenPhoton_E, &b_MatchedGenPhoton_E);
   fChain->SetBranchAddress("MatchedGenPhoton_et", &MatchedGenPhoton_et, &b_MatchedGenPhoton_et);
   fChain->SetBranchAddress("MatchedGenPhoton_pt", &MatchedGenPhoton_pt, &b_MatchedGenPhoton_pt);
   fChain->SetBranchAddress("MatchedGenPhoton_eta", &MatchedGenPhoton_eta, &b_MatchedGenPhoton_eta);
   fChain->SetBranchAddress("MatchedGenPhoton_phi", &MatchedGenPhoton_phi, &b_MatchedGenPhoton_phi);
   fChain->SetBranchAddress("MatchedGenPhoton_px", &MatchedGenPhoton_px, &b_MatchedGenPhoton_px);
   fChain->SetBranchAddress("MatchedGenPhoton_py", &MatchedGenPhoton_py, &b_MatchedGenPhoton_py);
   fChain->SetBranchAddress("MatchedGenPhoton_pz", &MatchedGenPhoton_pz, &b_MatchedGenPhoton_pz);
   fChain->SetBranchAddress("MatchedGenPhoton_vx", &MatchedGenPhoton_vx, &b_MatchedGenPhoton_vx);
   fChain->SetBranchAddress("MatchedGenPhoton_vy", &MatchedGenPhoton_vy, &b_MatchedGenPhoton_vy);
   fChain->SetBranchAddress("MatchedGenPhoton_vz", &MatchedGenPhoton_vz, &b_MatchedGenPhoton_vz);
   fChain->SetBranchAddress("MatchedGenPhoton_nMothers", &MatchedGenPhoton_nMothers, &b_MatchedGenPhoton_nMothers);
   fChain->SetBranchAddress("MatchedGenPhoton_mothersPdgid", &MatchedGenPhoton_mothersPdgid, &b_MatchedGenPhoton_mothersPdgid);
   fChain->SetBranchAddress("MatchedGenPhoton_mothersStatus", &MatchedGenPhoton_mothersStatus, &b_MatchedGenPhoton_mothersStatus);
   fChain->SetBranchAddress("PFPatJet_n", &PFPatJet_n, &b_PFPatJet_n);
   fChain->SetBranchAddress("PFPatJet_E", &PFPatJet_E, &b_PFPatJet_E);
   fChain->SetBranchAddress("PFPatJet_et", &PFPatJet_et, &b_PFPatJet_et);
   fChain->SetBranchAddress("PFPatJet_pt", &PFPatJet_pt, &b_PFPatJet_pt);
   fChain->SetBranchAddress("PFPatJet_eta", &PFPatJet_eta, &b_PFPatJet_eta);
   fChain->SetBranchAddress("PFPatJet_phi", &PFPatJet_phi, &b_PFPatJet_phi);
   fChain->SetBranchAddress("PFPatJet_theta", &PFPatJet_theta, &b_PFPatJet_theta);
   fChain->SetBranchAddress("PFPatJet_px", &PFPatJet_px, &b_PFPatJet_px);
   fChain->SetBranchAddress("PFPatJet_py", &PFPatJet_py, &b_PFPatJet_py);
   fChain->SetBranchAddress("PFPatJet_pz", &PFPatJet_pz, &b_PFPatJet_pz);
   fChain->SetBranchAddress("PFPatJet_vx", &PFPatJet_vx, &b_PFPatJet_vx);
   fChain->SetBranchAddress("PFPatJet_vy", &PFPatJet_vy, &b_PFPatJet_vy);
   fChain->SetBranchAddress("PFPatJet_vz", &PFPatJet_vz, &b_PFPatJet_vz);
   fChain->SetBranchAddress("PFPatJet_ChargedEmEnergy", &PFPatJet_ChargedEmEnergy, &b_PFPatJet_ChargedEmEnergy);
   fChain->SetBranchAddress("PFPatJet_ChargedEmEnergyFrac", &PFPatJet_ChargedEmEnergyFrac, &b_PFPatJet_ChargedEmEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_ChargedHadEnergy", &PFPatJet_ChargedHadEnergy, &b_PFPatJet_ChargedHadEnergy);
   fChain->SetBranchAddress("PFPatJet_ChargedHadEnergyFrac", &PFPatJet_ChargedHadEnergyFrac, &b_PFPatJet_ChargedHadEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_ChargedHadMult", &PFPatJet_ChargedHadMult, &b_PFPatJet_ChargedHadMult);
   fChain->SetBranchAddress("PFPatJet_ChargedMult", &PFPatJet_ChargedMult, &b_PFPatJet_ChargedMult);
   fChain->SetBranchAddress("PFPatJet_NConstituents", &PFPatJet_NConstituents, &b_PFPatJet_NConstituents);
   fChain->SetBranchAddress("PFPatJet_HFEMEnergy", &PFPatJet_HFEMEnergy, &b_PFPatJet_HFEMEnergy);
   fChain->SetBranchAddress("PFPatJet_HFEMEnergyFrac", &PFPatJet_HFEMEnergyFrac, &b_PFPatJet_HFEMEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_HFEMMult", &PFPatJet_HFEMMult, &b_PFPatJet_HFEMMult);
   fChain->SetBranchAddress("PFPatJet_HFHadEnergy", &PFPatJet_HFHadEnergy, &b_PFPatJet_HFHadEnergy);
   fChain->SetBranchAddress("PFPatJet_HFHadEnergyFrac", &PFPatJet_HFHadEnergyFrac, &b_PFPatJet_HFHadEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_HFHadMult", &PFPatJet_HFHadMult, &b_PFPatJet_HFHadMult);
   fChain->SetBranchAddress("PFPatJet_NeutralEmEnergy", &PFPatJet_NeutralEmEnergy, &b_PFPatJet_NeutralEmEnergy);
   fChain->SetBranchAddress("PFPatJet_NeutralEmEnergyFrac", &PFPatJet_NeutralEmEnergyFrac, &b_PFPatJet_NeutralEmEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_NeutralHadEnergy", &PFPatJet_NeutralHadEnergy, &b_PFPatJet_NeutralHadEnergy);
   fChain->SetBranchAddress("PFPatJet_NeutralHadEnergyFrac", &PFPatJet_NeutralHadEnergyFrac, &b_PFPatJet_NeutralHadEnergyFrac);
   fChain->SetBranchAddress("PFPatJet_NeutralHadMult", &PFPatJet_NeutralHadMult, &b_PFPatJet_NeutralHadMult);
   fChain->SetBranchAddress("PFPatJet_NeutralMult", &PFPatJet_NeutralMult, &b_PFPatJet_NeutralMult);
   fChain->SetBranchAddress("PFPatJet_jecUncertainity", &PFPatJet_jecUncertainity, &b_PFPatJet_jecUncertainity);
   fChain->SetBranchAddress("PFPatJet_puJetIdCutBased_MVA", &PFPatJet_puJetIdCutBased_MVA, &b_PFPatJet_puJetIdCutBased_MVA);
   fChain->SetBranchAddress("PFPatJet_puJetIdSimple_MVA", &PFPatJet_puJetIdSimple_MVA, &b_PFPatJet_puJetIdSimple_MVA);
   fChain->SetBranchAddress("PFPatJet_puJetIdFull_MVA", &PFPatJet_puJetIdFull_MVA, &b_PFPatJet_puJetIdFull_MVA);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdCutBased_loose", &PFPatJet_PassPUJetIdCutBased_loose, &b_PFPatJet_PassPUJetIdCutBased_loose);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdCutBased_medium", &PFPatJet_PassPUJetIdCutBased_medium, &b_PFPatJet_PassPUJetIdCutBased_medium);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdCutBased_tight", &PFPatJet_PassPUJetIdCutBased_tight, &b_PFPatJet_PassPUJetIdCutBased_tight);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdSimple_loose", &PFPatJet_PassPUJetIdSimple_loose, &b_PFPatJet_PassPUJetIdSimple_loose);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdSimple_medium", &PFPatJet_PassPUJetIdSimple_medium, &b_PFPatJet_PassPUJetIdSimple_medium);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdSimple_tight", &PFPatJet_PassPUJetIdSimple_tight, &b_PFPatJet_PassPUJetIdSimple_tight);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdFull_loose", &PFPatJet_PassPUJetIdFull_loose, &b_PFPatJet_PassPUJetIdFull_loose);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdFull_medium", &PFPatJet_PassPUJetIdFull_medium, &b_PFPatJet_PassPUJetIdFull_medium);
   fChain->SetBranchAddress("PFPatJet_PassPUJetIdFull_tight", &PFPatJet_PassPUJetIdFull_tight, &b_PFPatJet_PassPUJetIdFull_tight);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByTrackCountingHighEff", &PFPatJet_BJetDiscrByTrackCountingHighEff, &b_PFPatJet_BJetDiscrByTrackCountingHighEff);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByTrackCountingHighPur", &PFPatJet_BJetDiscrByTrackCountingHighPur, &b_PFPatJet_BJetDiscrByTrackCountingHighPur);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff", &PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff, &b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur", &PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur, &b_PFPatJet_BJetDiscrBySimpleSecondaryVertexHighPur);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff", &PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff, &b_PFPatJet_BJetDiscrByCombinedSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA", &PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA, &b_PFPatJet_BJetDiscrByCombinedSecondaryVertexMVA);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByjetProbabilityBJetTags", &PFPatJet_BJetDiscrByjetProbabilityBJetTags, &b_PFPatJet_BJetDiscrByjetProbabilityBJetTags);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrByjetBProbabilityBJetTags", &PFPatJet_BJetDiscrByjetBProbabilityBJetTags, &b_PFPatJet_BJetDiscrByjetBProbabilityBJetTags);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySoftElectronBJetTags", &PFPatJet_BJetDiscrBySoftElectronBJetTags, &b_PFPatJet_BJetDiscrBySoftElectronBJetTags);
   fChain->SetBranchAddress("PFPatJet_BJetDiscrBySoftMuonBJetTags", &PFPatJet_BJetDiscrBySoftMuonBJetTags, &b_PFPatJet_BJetDiscrBySoftMuonBJetTags);
   fChain->SetBranchAddress("hasMatchedGenParton_ToPFJet_", &hasMatchedGenParton_ToPFJet_, &b_hasMatchedGenParton_ToPFJet_);
   fChain->SetBranchAddress("hasMatchedGenJet_ToPFJet_", &hasMatchedGenJet_ToPFJet_, &b_hasMatchedGenJet_ToPFJet_);
   fChain->SetBranchAddress("MatchedGenPartonToPFJet_pdgId", &MatchedGenPartonToPFJet_pdgId, &b_MatchedGenPartonToPFJet_pdgId);
   fChain->SetBranchAddress("MatchedGenPartonToPFJet_status", &MatchedGenPartonToPFJet_status, &b_MatchedGenPartonToPFJet_status);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_E", &MatchedGenJetToPFJet_E, &b_MatchedGenJetToPFJet_E);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_pt", &MatchedGenJetToPFJet_pt, &b_MatchedGenJetToPFJet_pt);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_eta", &MatchedGenJetToPFJet_eta, &b_MatchedGenJetToPFJet_eta);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_phi", &MatchedGenJetToPFJet_phi, &b_MatchedGenJetToPFJet_phi);
   fChain->SetBranchAddress("MatchedGenJetToPFJet_nDaughters", &MatchedGenJetToPFJet_nDaughters, &b_MatchedGenJetToPFJet_nDaughters);
   fChain->SetBranchAddress("PFPatJet_JetPartonFlavor", &PFPatJet_JetPartonFlavor, &b_PFPatJet_JetPartonFlavor);
   fChain->SetBranchAddress("PFRecoJet_n", &PFRecoJet_n, &b_PFRecoJet_n);
   fChain->SetBranchAddress("PFRecoJet_E", &PFRecoJet_E, &b_PFRecoJet_E);
   fChain->SetBranchAddress("PFRecoJet_et", &PFRecoJet_et, &b_PFRecoJet_et);
   fChain->SetBranchAddress("PFRecoJet_pt", &PFRecoJet_pt, &b_PFRecoJet_pt);
   fChain->SetBranchAddress("PFRecoJet_eta", &PFRecoJet_eta, &b_PFRecoJet_eta);
   fChain->SetBranchAddress("PFRecoJet_phi", &PFRecoJet_phi, &b_PFRecoJet_phi);
   fChain->SetBranchAddress("PFRecoJet_theta", &PFRecoJet_theta, &b_PFRecoJet_theta);
   fChain->SetBranchAddress("PFRecoJet_px", &PFRecoJet_px, &b_PFRecoJet_px);
   fChain->SetBranchAddress("PFRecoJet_py", &PFRecoJet_py, &b_PFRecoJet_py);
   fChain->SetBranchAddress("PFRecoJet_pz", &PFRecoJet_pz, &b_PFRecoJet_pz);
   fChain->SetBranchAddress("PFRecoJet_vx", &PFRecoJet_vx, &b_PFRecoJet_vx);
   fChain->SetBranchAddress("PFRecoJet_vy", &PFRecoJet_vy, &b_PFRecoJet_vy);
   fChain->SetBranchAddress("PFRecoJet_vz", &PFRecoJet_vz, &b_PFRecoJet_vz);
   fChain->SetBranchAddress("PFRecoJet_ChargedEmEnergy", &PFRecoJet_ChargedEmEnergy, &b_PFRecoJet_ChargedEmEnergy);
   fChain->SetBranchAddress("PFRecoJet_ChargedEmEnergyFrac", &PFRecoJet_ChargedEmEnergyFrac, &b_PFRecoJet_ChargedEmEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_ChargedHadEnergy", &PFRecoJet_ChargedHadEnergy, &b_PFRecoJet_ChargedHadEnergy);
   fChain->SetBranchAddress("PFRecoJet_ChargedHadEnergyFrac", &PFRecoJet_ChargedHadEnergyFrac, &b_PFRecoJet_ChargedHadEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_ChargedHadMult", &PFRecoJet_ChargedHadMult, &b_PFRecoJet_ChargedHadMult);
   fChain->SetBranchAddress("PFRecoJet_ChargedMult", &PFRecoJet_ChargedMult, &b_PFRecoJet_ChargedMult);
   fChain->SetBranchAddress("PFRecoJet_NConstituents", &PFRecoJet_NConstituents, &b_PFRecoJet_NConstituents);
   fChain->SetBranchAddress("PFRecoJet_HFEMEnergy", &PFRecoJet_HFEMEnergy, &b_PFRecoJet_HFEMEnergy);
   fChain->SetBranchAddress("PFRecoJet_HFEMEnergyFrac", &PFRecoJet_HFEMEnergyFrac, &b_PFRecoJet_HFEMEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_HFEMMult", &PFRecoJet_HFEMMult, &b_PFRecoJet_HFEMMult);
   fChain->SetBranchAddress("PFRecoJet_HFHadEnergy", &PFRecoJet_HFHadEnergy, &b_PFRecoJet_HFHadEnergy);
   fChain->SetBranchAddress("PFRecoJet_HFHadEnergyFrac", &PFRecoJet_HFHadEnergyFrac, &b_PFRecoJet_HFHadEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_HFHadMult", &PFRecoJet_HFHadMult, &b_PFRecoJet_HFHadMult);
   fChain->SetBranchAddress("PFRecoJet_NeutralEmEnergy", &PFRecoJet_NeutralEmEnergy, &b_PFRecoJet_NeutralEmEnergy);
   fChain->SetBranchAddress("PFRecoJet_NeutralEmEnergyFrac", &PFRecoJet_NeutralEmEnergyFrac, &b_PFRecoJet_NeutralEmEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_NeutralHadEnergy", &PFRecoJet_NeutralHadEnergy, &b_PFRecoJet_NeutralHadEnergy);
   fChain->SetBranchAddress("PFRecoJet_NeutralHadEnergyFrac", &PFRecoJet_NeutralHadEnergyFrac, &b_PFRecoJet_NeutralHadEnergyFrac);
   fChain->SetBranchAddress("PFRecoJet_NeutralHadMult", &PFRecoJet_NeutralHadMult, &b_PFRecoJet_NeutralHadMult);
   fChain->SetBranchAddress("PFRecoJet_NeutralMult", &PFRecoJet_NeutralMult, &b_PFRecoJet_NeutralMult);
   fChain->SetBranchAddress("PFRecoJet_jecUncertainity", &PFRecoJet_jecUncertainity, &b_PFRecoJet_jecUncertainity);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByTrackCountingHighEff", &PFRecoJet_BJetDiscrByTrackCountingHighEff, &b_PFRecoJet_BJetDiscrByTrackCountingHighEff);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByTrackCountingHighPur", &PFRecoJet_BJetDiscrByTrackCountingHighPur, &b_PFRecoJet_BJetDiscrByTrackCountingHighPur);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff", &PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff, &b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur", &PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur, &b_PFRecoJet_BJetDiscrBySimpleSecondaryVertexHighPur);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff", &PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff, &b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexHighEff);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA", &PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA, &b_PFRecoJet_BJetDiscrByCombinedSecondaryVertexMVA);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByjetProbabilityBJetTags", &PFRecoJet_BJetDiscrByjetProbabilityBJetTags, &b_PFRecoJet_BJetDiscrByjetProbabilityBJetTags);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrByjetBProbabilityBJetTags", &PFRecoJet_BJetDiscrByjetBProbabilityBJetTags, &b_PFRecoJet_BJetDiscrByjetBProbabilityBJetTags);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySoftElectronBJetTags", &PFRecoJet_BJetDiscrBySoftElectronBJetTags, &b_PFRecoJet_BJetDiscrBySoftElectronBJetTags);
   fChain->SetBranchAddress("PFRecoJet_BJetDiscrBySoftMuonBJetTags", &PFRecoJet_BJetDiscrBySoftMuonBJetTags, &b_PFRecoJet_BJetDiscrBySoftMuonBJetTags);
   fChain->SetBranchAddress("GenJet_n", &GenJet_n, &b_GenJet_n);
   fChain->SetBranchAddress("GenJet_E", &GenJet_E, &b_GenJet_E);
   fChain->SetBranchAddress("GenJet_et", &GenJet_et, &b_GenJet_et);
   fChain->SetBranchAddress("GenJet_pt", &GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("GenJet_eta", &GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_phi", &GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_px", &GenJet_px, &b_GenJet_px);
   fChain->SetBranchAddress("GenJet_py", &GenJet_py, &b_GenJet_py);
   fChain->SetBranchAddress("GenJet_pz", &GenJet_pz, &b_GenJet_pz);
   fChain->SetBranchAddress("GenJet_vx", &GenJet_vx, &b_GenJet_vx);
   fChain->SetBranchAddress("GenJet_vy", &GenJet_vy, &b_GenJet_vy);
   fChain->SetBranchAddress("GenJet_vz", &GenJet_vz, &b_GenJet_vz);
   fChain->SetBranchAddress("GenJet_nDaughters", &GenJet_nDaughters, &b_GenJet_nDaughters);
   fChain->SetBranchAddress("GenJet_daughters_status", &GenJet_daughters_status, &b_GenJet_daughters_status);
   fChain->SetBranchAddress("GenJet_daughters_pdgid", &GenJet_daughters_pdgid, &b_GenJet_daughters_pdgid);
   fChain->SetBranchAddress("GenJet_daughters_E", &GenJet_daughters_E, &b_GenJet_daughters_E);
   fChain->SetBranchAddress("GenJet_daughters_pt", &GenJet_daughters_pt, &b_GenJet_daughters_pt);
   fChain->SetBranchAddress("GenJet_daughters_eta", &GenJet_daughters_eta, &b_GenJet_daughters_eta);
   fChain->SetBranchAddress("GenJet_daughters_phi", &GenJet_daughters_phi, &b_GenJet_daughters_phi);
   fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
   fChain->SetBranchAddress("Vertex_x", &Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex_y", &Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex_z", &Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex_chi2", &Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex_nchi2", &Vertex_nchi2, &b_Vertex_nchi2);
   fChain->SetBranchAddress("Vertex_ndof", &Vertex_ndof, &b_Vertex_ndof);
   fChain->SetBranchAddress("Vertex_tracksSize", &Vertex_tracksSize, &b_Vertex_tracksSize);
   fChain->SetBranchAddress("Vertex_isFake", &Vertex_isFake, &b_Vertex_isFake);
   fChain->SetBranchAddress("Vertex_isValid", &Vertex_isValid, &b_Vertex_isValid);
   fChain->SetBranchAddress("Vertex_d0", &Vertex_d0, &b_Vertex_d0);
   fChain->SetBranchAddress("Tracks_n", &Tracks_n, &b_Tracks_n);
   fChain->SetBranchAddress("Track_pt", &Track_pt, &b_Track_pt);
   fChain->SetBranchAddress("Track_px", &Track_px, &b_Track_px);
   fChain->SetBranchAddress("Track_py", &Track_py, &b_Track_py);
   fChain->SetBranchAddress("Track_pz", &Track_pz, &b_Track_pz);
   fChain->SetBranchAddress("Track_vx", &Track_vx, &b_Track_vx);
   fChain->SetBranchAddress("Track_vy", &Track_vy, &b_Track_vy);
   fChain->SetBranchAddress("Track_vz", &Track_vz, &b_Track_vz);
   fChain->SetBranchAddress("Track_eta", &Track_eta, &b_Track_eta);
   fChain->SetBranchAddress("Track_phi", &Track_phi, &b_Track_phi);
   fChain->SetBranchAddress("Track_theta", &Track_theta, &b_Track_theta);
   fChain->SetBranchAddress("Track_chi2", &Track_chi2, &b_Track_chi2);
   fChain->SetBranchAddress("IsScrapingEvent_", &IsScrapingEvent_, &b_IsScrapingEvent_);
   fChain->SetBranchAddress("Scraping_FractionOfGoodTracks", &Scraping_FractionOfGoodTracks, &b_Scraping_FractionOfGoodTracks);
   fChain->SetBranchAddress("HLT_Photon_triggers", &HLT_Photon_triggers, &b_HLT_Photon_triggers);
   fChain->SetBranchAddress("HLT_Photon_trig_prescales", &HLT_Photon_trig_prescales, &b_HLT_Photon_trig_prescales);
   fChain->SetBranchAddress("HLT_Photon_ifTriggerPassed", &HLT_Photon_ifTriggerPassed, &b_HLT_Photon_ifTriggerPassed);
   fChain->SetBranchAddress("HLT_Photon_nTriggers", &HLT_Photon_nTriggers, &b_HLT_Photon_nTriggers);
   fChain->SetBranchAddress("HLT_Photon_triggerIndex", &HLT_Photon_triggerIndex, &b_HLT_Photon_triggerIndex);
   fChain->SetBranchAddress("HLT_Photon_nFilters", &HLT_Photon_nFilters, &b_HLT_Photon_nFilters);
   fChain->SetBranchAddress("HLT_Photon_FilterNames", &HLT_Photon_FilterNames, &b_HLT_Photon_FilterNames);
   fChain->SetBranchAddress("HLT_Photon_trigger_FilterStartPosition", &HLT_Photon_trigger_FilterStartPosition, &b_HLT_Photon_trigger_FilterStartPosition);
   fChain->SetBranchAddress("HLT_Photon_trigger_FilterEndPosition", &HLT_Photon_trigger_FilterEndPosition, &b_HLT_Photon_trigger_FilterEndPosition);
   fChain->SetBranchAddress("HLT_Photon_FilterObjects_pt", &HLT_Photon_FilterObjects_pt, &b_HLT_Photon_FilterObjects_pt);
   fChain->SetBranchAddress("HLT_Photon_FilterObjects_eta", &HLT_Photon_FilterObjects_eta, &b_HLT_Photon_FilterObjects_eta);
   fChain->SetBranchAddress("HLT_Photon_FilterObjects_phi", &HLT_Photon_FilterObjects_phi, &b_HLT_Photon_FilterObjects_phi);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("sigma", &sigma, &b_sigma);
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("sigma25", &sigma25, &b_sigma25);
   Notify();
}

Bool_t MyEvent_MC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyEvent_MC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyEvent_MC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyEvent_MC_cxx
EOF


end  ## foreach i loop

