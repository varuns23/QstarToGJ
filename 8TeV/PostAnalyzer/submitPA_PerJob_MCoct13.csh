#!/bin/tcsh

### Here I have changed Variable binning for mass and pt, As I am changing both the cuts. 
### Also Pt bins start from 20.0 to put a cut on 170 geV.

setenv pwd $PWD
set sampleIndex = 0

setenv OutEos /eos/uscms/store/user/varun/work/QstarGJ_2012Analysis/PostAnalyzer/Preliminary/September2013/Pt170MassCut/dEta2p0M560
#setenv OutEos /eos/uscms/store/user/varun/work/QstarGJ_2012Analysis/PostAnalyzer/Preliminary/September2013/Pt170NoMassCut/dEta1p0

set M = 1 ## set 1 for MC
set D = 0 ## set 1 for Data

#--->Qstar to GJ 
#foreach i (QstarToGJ_M_700  QstarToGJ_M_1000  QstarToGJ_M_1200  QstarToGJ_M_1500  QstarToGJ_M_1700  QstarToGJ_M_2000  QstarToGJ_M_2500  QstarToGJ_M_3000  QstarToGJ_M_3500 QstarToGJ_M_4000 QstarToGJ_M_4500)
#set XS  = (24.85            4.18              1.552             0.4124            0.1843            0.05858           0.009768          0.001755          0.0003241        0.00006045       0.00001170) 
#set totalEvents = (60065.0  60192.0           60120.0           60138.0           60256.0           120105.0          120127.0          120032.0          160095.0         159703           306561)
#setenv sourceDir /pnfs/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/QstarSignal/${i}/
#setenv InputFilesPath /pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/QstarSignal/${i}/
#setenv OutPath ${OutEos}/Signal
#set filesPerJob = 100

foreach i (QstarToGJ_M_fhalf_700  QstarToGJ_M_700_CheckNewConfig )
set XS  = (6.249                  24.85 )
set totalEvents = (59013.0        60065.0 )

#---->Qstar to GJ half coupling
#foreach i (QstarToGJ_M_fhalf_700  QstarToGJ_M_fhalf_1000  QstarToGJ_M_fhalf_1500  QstarToGJ_M_fhalf_2000  QstarToGJ_M_fhalf_2500  QstarToGJ_M_fhalf_3000  QstarToGJ_M_fhalf_3500 QstarToGJ_M_fhalf_4000 QstarToGJ_M_fhalf_4500)
#set XS  = (6.249                  1.064                   0.105                   0.01484                 0.002425                4.304E-4                7.586E-5               1.258E-5               1.929E-6)
#set totalEvents = (59013.0        60196.0         60053.0                 120012.0                120054.0                120059.0                160008.0               159808                 160000)
setenv sourceDir /pnfs/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/QstarSignal/${i}/
setenv InputFilesPath /pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/QstarSignal/${i}/
setenv OutPath ${OutEos}/Signal
set filesPerJob = 100

#----> GammaJet Bkg ---- Pythia
#foreach i (G_Pt_120to170  G_Pt_170to300  G_Pt_300to470  G_Pt_470to800  G_Pt_800to1400  G_Pt_1400to1800  G_Pt_1800toInf  )
#set XS  = (108.0068       30.12207       2.138632       0.2119244      0.007077847     4.510327E-5      1.867141E-6)
#set totalEvents = (2000043.0 2000069.0   2000130.0      1975231.0      1973504.0       1983890.0        1939122.0) ## Reduced 10000 events in 1400to1800
#setenv sourceDir /pnfs/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/PhotonJet/${i}/
#setenv InputFilesPath /pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/PhotonJet/${i}/
#setenv OutPath ${OutEos}/GJBkg
#set filesPerJob = 25

#----> DiJet Bkg ----- Pythia
#foreach i (QCD_Pt120to170  QCD_Pt170to300  QCD_Pt300to470  QCD_Pt470to600  QCD_Pt600to800  QCD_Pt800to1000 QCD_Pt1000to1400 QCD_Pt1400to1800 QCD_Pt1800toInf )
#set XS  = (156293.3        34138.15        1759.549        113.8791        26.9921         3.550036        0.737844         0.03352235       0.001829005)
#set totalEvents = (5985732.0 5814398.0     5978500.0       3994848.0       3996864.0       3998563.0       1964088.0        2000062.0        977586.0)
#setenv sourceDir /pnfs/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/DiJet/${i}/
#setenv InputFilesPath /pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/DiJet/${i}/
#setenv OutPath ${OutEos}/DiJetBkg
#set filesPerJob = 25

#----------- EWK Background
#foreach i (WGToLNuG    WJetsToLNu    ZG_Inclusive)
#set XS  = (553.9       37510.0       172.1)
#set totalEvents = (4802358.0  18393090.0  6321549.0      )
#setenv sourceDir /pnfs/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/${i}/
#setenv InputFilesPath /pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/MC/Summer12/53X/${i}/
#setenv OutPath ${OutEos}/EWKBkg
#set filesPerJob = 20

#-----> Data ---- 22 Jan
#foreach i (Run2012D_22Jan  Run2012B_22Jan  Run2012C_22Jan  Run2012A_22Jan)
#set XS  = (1               1               1               1)
#set totalEvents = (1       1               1               1)
##set filesPerJob = (10      25              25              25)
#setenv sourceDir /pnfs/cms/WAX/11/store/user/varun/2012/Data22JanReReco/${i}/
#setenv InputFilesPath /pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/Data22JanReReco/${i}/
#setenv OutPath ${OutEos}/Data


setenv FileNameTag ${i}

set path_eos  = 0  ## if input files are in eos area
set path_pnfs = 1  ## if input files are in pnfs area

##### Check for inputDir---InCase of EOS area make change to ------++++++FullPathInputFile++++++++
#-------------------------------------------------------------------------------------------

@ sampleIndex = ${sampleIndex} + 1

setenv destination ${OutPath}

if( ! -d ${destination}) then
echo "Making Directory ${destination} "
mkdir ${destination}
chmod 775 ${destination}
endif

##_______________________________________________________________________________
## changes to read each file in dir and sub seperate jobs.
##-----------------------------------------------------------------

#Fill the range for input to code e.g 1-10 then r = 10
set r = ${filesPerJob}

set sf = 1          ## start file
set ef = ${r}       ## end file

ls -ltr --format=single-column ${sourceDir} > ${FileNameTag}_dataset.txt


##====================No CHANGE BELOW============================
#setenv datafile  ${pwd}/${FileNameTag}_dataset.txt
setenv datafile  ${FileNameTag}_dataset.txt
#==================================================================


set file_length=`wc -l $datafile | cut -c1-4`                  # Count Lines
set p = 0

set Tot = ${file_length}
#set Tot = 1

#run till the last line of the input files
while (${sf} <= ${Tot})

@ p = ${p} + 1
set DataName=`tail -n +$p ${datafile} | head -1`                  # Get name of dataset

#####More changes in constructor, in the outputfile and at the end of file--------------------------------


##################PostAnalyzeMC.C#####################
cat>PostAnalyzerMC.C<<EOF
#define PostAnalyzerMC_cxx
#include "PostAnalyzerMC.h"

    using namespace std;
    using namespace ROOT;

void PostAnalyzerMC::Loop()
{
    // Give values to all the parameters at the beginning 
    //Luminosity 
    Lumi           =  19740.0 ; //11361.7 // (pb^{-1})

    //Vertex Selection
    Cvertex_z      =  24.0 ; //(cm)
    Cvertex_ndof   =   4.  ;
    Cvertex_rho    =   2.0 ; //(cm)

    //Fiducial Cuts
    photon_pt_cut  = 170.0 ; //(GeV)   // Changing from 150-> 200
    photon_eta_cut = 1.4442;

    jet_pt_cut     = 170.0 ; //(GeV)   // Changing from 150-> 200
    jet_eta_cut    = 3.0   ; 

    mass_cut       = 560.0 ; //(GeV)

    PhotJetDPhi_cut    = 1.5   ;
    PhotJetDEta_cut    = 2.0   ;

    Bool_t passHLT;
    Bool_t nonScraping, primaryVtx;
    Bool_t tightJetID;

    MC   = ${M};
    DATA = ${D};

    //OUTPUT FILE
    TString OutputFile = "${FileNameTag}_${p}";
    TString OutputPath = "${destination}/";
    f1 = new TFile(OutputPath+OutputFile+".root","recreate");

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    const int nbins = 15;
    TString  CutFlowLabel[nbins] ={ "Total", "HLT", "Scraping", "PrimaryVtx", "PhotonID", "PhotonPt", "PhotonEta", "JetID", "JetPt ", "JetEta ", "Dphi1.5", "DEta1.0", "MassCut", "Dphi2.0", "Dphi2.5"};
    Double_t CutFlowNumber[nbins]={  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ;
    Double_t CutExpFlowNumber[nbins]={  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ;
    Double_t CutFlowBins[nbins+1]={  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    h_CutFlowTable    =new TH1F("h_CutFlowTable"," cut flow of selection",nbins,CutFlowBins);
    h_CutExpFlowTable =new TH1F("h_CutExpFlowTable","Expected cut flow of selection",nbins,CutFlowBins);

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    cout<<" Will Analyze = "<<nentries<<" events"<<endl;

    
    BookHistos();
    if(MC)LumiReWeighting();

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	nb = fChain->GetEntry(jentry);   nbytes += nb;

	Double_t pu_weight=1.0;  // weight from pileup---------------------------

	if(MC){
	    pu_weight = puweight(trueInteractions);
	    h_trueinteractions_MC->Fill(trueInteractions);
	    h_trueinteractions_MC_PU->Fill(trueInteractions,pu_weight);
	}

	if(MC)  OnlyEvtWeight = Lumi*($XS[${sampleIndex}]/$totalEvents[${sampleIndex}]);
	if(DATA)OnlyEvtWeight = 1.;
	EvtWeight = pu_weight*OnlyEvtWeight ;

	PC = -1;
	JC = -1;
	goodVertex = 0;

	nonScraping = NonScraping();
	primaryVtx  = PrimaryVertex(goodVertex);

	// Selecting photon candidate using tight photon ID.
	foundPhoton.clear();
	for(int ipho=0; ipho<Photon_n; ++ipho){
	    if( NoSpike(ipho)   &&  TightPhotonPFIso(ipho) ){
		foundPhoton.push_back(ipho);
	    }
	}
	if(foundPhoton.size() != 0) PC = foundPhoton[0];
	
	h_nIsoPhoton[0]->Fill(foundPhoton.size());
	h_nPhoton[0]->Fill(Photon_n);
	for(Int_t i_p = 0 ; i_p < foundPhoton.size() ; ++i_p){
	    h_nIsoPhotonPt[0]->Fill(Photon_pt[foundPhoton[i_p]],i_p+1);
	}


	// Selecting Jet candidate using tight Jet ID.
	foundJet.clear();
	if(PC >= 0){
	    for(int ijet=0; ijet<pfJet_n; ++ijet){
		Float_t drp = -1.0;
		drp = getDR(pfJet_eta[ijet],Photon_eta[PC],pfJet_phi[ijet], Photon_phi[PC]);
		if( drp > 0.5  &&  TightJetID(ijet) ){
		    foundJet.push_back(ijet);
		}
	    }
	}
	if(foundJet.size() != 0) JC = foundJet[0];
	
	h_nIsoJet[0]->Fill(foundJet.size());
	h_nJet[0]->Fill(pfJet_n);
	for(Int_t i_j = 0 ; i_j < foundJet.size() ; ++i_j){
	    h_nIsoJetPt[0]->Fill(pfJet_pt[foundJet[i_j]],i_j+1);
	}

	
	h_Vertices[0]->Fill(goodVertex);      
	h_Vertices[1]->Fill(goodVertex, pu_weight);

	//Kinematical Cuts  ----------------
	Bool_t PhoPt    = ( Photon_pt[PC] > photon_pt_cut );
	Bool_t PhoEtaEB = ( fabs(Photon_sc_eta[PC]) < photon_eta_cut && Photon_isEB[PC]);

	Bool_t JetPt    = ( pfJet_pt[JC] > jet_pt_cut );
	Bool_t JetEtaEB = ( fabs(pfJet_eta[JC]) < jet_eta_cut);

	Bool_t MassCut  = ( getMass(PC,JC) > mass_cut );

	Bool_t GJDeltaPhi = getDPhi(Photon_phi[PC], pfJet_phi[JC])    > PhotJetDPhi_cut;
	Bool_t GJDeltaEta = getDEta(Photon_sc_eta[PC], pfJet_eta[JC]) < PhotJetDEta_cut;

	Bool_t InnerRadius = ( fabs(Photon_sc_eta[PC]) < 0.7 && fabs(pfJet_eta[JC]) < 0.7) ;
	Bool_t OuterRadius = ( fabs(Photon_sc_eta[PC]) > 0.7 && fabs(Photon_sc_eta[PC]) < photon_eta_cut && fabs(pfJet_eta[JC]) > 0.7 && fabs(pfJet_eta[JC]) < jet_eta_cut ) ;

	if(MC){
	    passHLT     = true;
	    nonScraping = true;
	}

	//--- Clearing JES vectors
	mass_JESup = 0.0 ; mass_JESdown = 0.0 ;

	h_ptPhoton[0]->Fill(Photon_pt[0],EvtWeight);
	h_etaPhoton[0]->Fill(Photon_sc_eta[0],EvtWeight);
	h_ptJet[0]->Fill(pfJet_pt[0],EvtWeight);
	h_etaJet[0]->Fill(pfJet_eta[0],EvtWeight);
	h_mass_VarBin[0]->Fill(getMass(0,0),EvtWeight);
	h_mass_bin25[0]->Fill(getMass(0,0),EvtWeight);
	h_ptPFMet[0]->Fill(PFMetPt[0],EvtWeight);
	h_SumEtPFMet[0]->Fill(PFMetSumEt[0],EvtWeight);

	h_mass_bin1[0]->Fill(getMass(0,0),EvtWeight);
	h_mass_bin40[0]->Fill(getMass(0,0),EvtWeight);
	h_Photon_SigmaIetaIeta[0]->Fill(Photon_SigmaIetaIeta[0],EvtWeight);
	h_PtPhotJet[0]->Fill(Photon_pt[0],pfJet_pt[0],EvtWeight);
	h_dEtadPhi[0]->Fill(getDEta(Photon_sc_eta[0],pfJet_eta[0]),getDPhi(Photon_phi[0],pfJet_phi[0]),EvtWeight);
	h_DR_PhotonJet[0]->Fill(getDR(Photon_sc_eta[0], pfJet_eta[0],Photon_phi[0], pfJet_phi[0]),EvtWeight);
	h_YPhotJet[0]->Fill(getRapidity(Photon_E[0],Photon_pz[0]),getRapidity(pfJet_E[0],pfJet_pz[0]),EvtWeight);
	h_dY[0]->Fill(getRapidity(Photon_E[0],Photon_pz[0])-getRapidity(pfJet_E[0],pfJet_pz[0]),EvtWeight);
	h_etaPhotJet[0]->Fill(Photon_sc_eta[0],pfJet_eta[0],EvtWeight);
	h_dEta[0]->Fill(getDEta(Photon_sc_eta[0],pfJet_eta[0]),EvtWeight);
	h_phiPhotJet[0]->Fill(Photon_phi[0],pfJet_phi[0],EvtWeight);
	h_dphi[0]->Fill(getDPhi(Photon_phi[0],pfJet_phi[0]),EvtWeight);

	h_Photon_SigmaEtaEta[0]->Fill(Photon_SigmaEtaEta[0],EvtWeight);
	h_Photon_SigmaPhiPhi[0]->Fill(Photon_SigmaPhiPhi[0],EvtWeight);
	h_Photon_SigmaIphiIphi[0]->Fill(Photon_SigmaIphiIphi[0],EvtWeight);
	h_Photon_swissCross[0]->Fill(Photon_swissCross[0],EvtWeight);
	h_Photon_r9[0]->Fill(Photonr9[0],EvtWeight);

	h_PFiso_ChargedBarrel[0]->Fill(PFiso_Charged03[0],EvtWeight);
	h_PFiso_PhotonBarrel[0]->Fill(PFiso_Photon03[0],EvtWeight);
	h_PFiso_NeutralBarrel[0]->Fill(PFiso_Neutral03[0],EvtWeight);
	h_PFiso_SumBarrel[0]->Fill(PFiso_Sum03[0],EvtWeight);
	h_PFiso_ElectronvetoBarrel[0]->Fill(Photon_Electronveto[0],EvtWeight);
	h_HoEnewBarrel[0]->Fill(Photon_HoEnew[PC],EvtWeight);

	h_Photon_ecalRecHitSumEtConeDR04[0]->Fill(Photon_ecalRecHitSumEtConeDR04[0],EvtWeight);
	h_Photon_hcalTowerSumEtConeDR04[0]->Fill(Photon_hcalTowerSumEtConeDR04[0],EvtWeight);
	h_Photon_trkSumPtHollowConeDR04[0]->Fill(Photon_trkSumPtHollowConeDR04[0],EvtWeight);
	h_HoE[0]->Fill(Photon_HoE[0],EvtWeight);

	h_pfjet_NEF[0]->Fill(pfjet_NEF[0],EvtWeight);
	h_pfjet_NHF[0]->Fill(pfjet_NHF[0],EvtWeight);
	h_pfjet_CEF[0]->Fill(pfjet_CEF[0],EvtWeight);
	h_pfjet_CHF[0]->Fill(pfjet_CHF[0],EvtWeight);
	h_pfjet_NConstituents[0]->Fill(pfjet_NConstituents[0],EvtWeight);
	h_pfjet_ChargeMultiplicity[0]->Fill(pfjet_NCH[0],EvtWeight);

	h_ptSecondPhoton[0]->Fill(Photon_pt[1],EvtWeight);
	h_etaSecondPhoton[0]->Fill(Photon_sc_eta[1],EvtWeight);

	CutFlowNumber[0]++;
	CutExpFlowNumber[0] += EvtWeight;

	if(passHLT){
	    CutFlowNumber[1]++;
	    CutExpFlowNumber[1] += EvtWeight;

	    if(nonScraping){
		CutFlowNumber[2]++;
		CutExpFlowNumber[2] += EvtWeight;

		if(primaryVtx){
		    CutFlowNumber[3]++;
		    CutExpFlowNumber[3] += EvtWeight;

		    if(PC > -1){
			CutFlowNumber[4]++;
			CutExpFlowNumber[4] += EvtWeight;

			h_ptPhoton[1]->Fill(Photon_pt[PC],EvtWeight);
			h_etaPhoton[1]->Fill(Photon_sc_eta[PC],EvtWeight);
			h_ptJet[1]->Fill(pfJet_pt[0],EvtWeight);
			h_etaJet[1]->Fill(pfJet_eta[0],EvtWeight);
			h_mass_VarBin[1]->Fill(getMass(PC,0),EvtWeight);
			h_mass_bin25[1]->Fill(getMass(PC,0),EvtWeight);
			h_ptPFMet[1]->Fill(PFMetPt[0],EvtWeight);
			h_SumEtPFMet[1]->Fill(PFMetSumEt[0],EvtWeight);

			h_mass_bin1[1]->Fill(getMass(PC,0),EvtWeight);
			h_mass_bin40[1]->Fill(getMass(PC,0),EvtWeight);
			h_Photon_SigmaIetaIeta[1]->Fill(Photon_SigmaIetaIeta[PC],EvtWeight);
			h_PtPhotJet[1]->Fill(Photon_pt[PC],pfJet_pt[0],EvtWeight);
			h_dEtadPhi[1]->Fill(getDEta(Photon_sc_eta[PC],pfJet_eta[0]),getDPhi(Photon_phi[PC],pfJet_phi[0]),EvtWeight);
			h_DR_PhotonJet[1]->Fill(getDR(Photon_sc_eta[PC], pfJet_eta[0],Photon_phi[PC], pfJet_phi[0]),EvtWeight);
			h_YPhotJet[1]->Fill(getRapidity(Photon_E[PC],Photon_pz[PC]),getRapidity(pfJet_E[0],pfJet_pz[0]),EvtWeight);
			h_dY[1]->Fill(getRapidity(Photon_E[PC],Photon_pz[PC])-getRapidity(pfJet_E[0],pfJet_pz[0]),EvtWeight);
			h_etaPhotJet[1]->Fill(Photon_sc_eta[PC],pfJet_eta[0],EvtWeight);
			h_dEta[1]->Fill(getDEta(Photon_sc_eta[PC],pfJet_eta[0]),EvtWeight);
			h_phiPhotJet[1]->Fill(Photon_phi[PC],pfJet_phi[0],EvtWeight);
			h_dphi[1]->Fill(getDPhi(Photon_phi[PC],pfJet_phi[0]),EvtWeight);

			if(PhoPt){
			    CutFlowNumber[5]++;
			    CutExpFlowNumber[5] += EvtWeight;

			    h_ptPhoton[2]->Fill(Photon_pt[PC],EvtWeight);
			    h_etaPhoton[2]->Fill(Photon_sc_eta[PC],EvtWeight);
			    h_ptJet[2]->Fill(pfJet_pt[JC],EvtWeight);
			    h_etaJet[2]->Fill(pfJet_eta[JC],EvtWeight);
			    h_mass_VarBin[2]->Fill(getMass(PC,JC),EvtWeight);
			    h_mass_bin25[2]->Fill(getMass(PC,JC),EvtWeight);
			    h_ptPFMet[2]->Fill(PFMetPt[0],EvtWeight);
			    h_SumEtPFMet[2]->Fill(PFMetSumEt[0],EvtWeight);

			    if(PhoEtaEB){
				CutFlowNumber[6]++;
				CutExpFlowNumber[6] += EvtWeight;

				h_ptPhoton[3]->Fill(Photon_pt[PC],EvtWeight);
				h_etaPhoton[3]->Fill(Photon_sc_eta[PC],EvtWeight);
				h_ptJet[3]->Fill(pfJet_pt[JC],EvtWeight);
				h_etaJet[3]->Fill(pfJet_eta[JC],EvtWeight);
				h_mass_VarBin[3]->Fill(getMass(PC,JC),EvtWeight);
				h_mass_bin25[3]->Fill(getMass(PC,JC),EvtWeight);
				h_ptPFMet[3]->Fill(PFMetPt[0],EvtWeight);
				h_SumEtPFMet[3]->Fill(PFMetSumEt[0],EvtWeight);

				if(JC > -1){
				    CutFlowNumber[7]++;
				    CutExpFlowNumber[7] += EvtWeight;

				    h_ptPhoton[4]->Fill(Photon_pt[PC],EvtWeight);
				    h_etaPhoton[4]->Fill(Photon_sc_eta[PC],EvtWeight);
				    h_ptJet[4]->Fill(pfJet_pt[JC],EvtWeight);
				    h_etaJet[4]->Fill(pfJet_eta[JC],EvtWeight);
				    h_mass_VarBin[4]->Fill(getMass(PC,JC),EvtWeight);
				    h_mass_bin25[4]->Fill(getMass(PC,JC),EvtWeight);
				    h_ptPFMet[4]->Fill(PFMetPt[0],EvtWeight);
				    h_SumEtPFMet[4]->Fill(PFMetSumEt[0],EvtWeight);

				    if(JetPt){
					CutFlowNumber[8]++;
					CutExpFlowNumber[8] += EvtWeight;

					h_ptPhoton[5]->Fill(Photon_pt[PC],EvtWeight);
					h_etaPhoton[5]->Fill(Photon_sc_eta[PC],EvtWeight);
					h_ptJet[5]->Fill(pfJet_pt[JC],EvtWeight);
					h_etaJet[5]->Fill(pfJet_eta[JC],EvtWeight);
					h_mass_VarBin[5]->Fill(getMass(PC,JC),EvtWeight);
					h_mass_bin25[5]->Fill(getMass(PC,JC),EvtWeight);
					h_ptPFMet[5]->Fill(PFMetPt[0],EvtWeight);
					h_SumEtPFMet[5]->Fill(PFMetSumEt[0],EvtWeight);

					if(JetEtaEB){
					    CutFlowNumber[9]++;
					    CutExpFlowNumber[9] += EvtWeight;

					    h_ptPhoton[6]->Fill(Photon_pt[PC],EvtWeight);
					    h_etaPhoton[6]->Fill(Photon_sc_eta[PC],EvtWeight);
					    h_ptJet[6]->Fill(pfJet_pt[JC],EvtWeight);
					    h_etaJet[6]->Fill(pfJet_eta[JC],EvtWeight);
					    h_mass_VarBin[6]->Fill(getMass(PC,JC),EvtWeight);
					    h_mass_bin25[6]->Fill(getMass(PC,JC),EvtWeight);
					    h_ptPFMet[6]->Fill(PFMetPt[0],EvtWeight);
					    h_SumEtPFMet[6]->Fill(PFMetSumEt[0],EvtWeight);

					    if(GJDeltaPhi){
						CutFlowNumber[10]++;
						CutExpFlowNumber[10] += EvtWeight;

						if(GJDeltaEta){
						    CutFlowNumber[11]++;
						    CutExpFlowNumber[11] += EvtWeight;

						    if(MassCut){
							CutFlowNumber[12]++;
							CutExpFlowNumber[12] += EvtWeight;

							if(getDPhi(Photon_phi[PC], pfJet_phi[JC]) > 2.0){
							    CutFlowNumber[13]++;
							    CutExpFlowNumber[13] += EvtWeight;
							}
							if(getDPhi(Photon_phi[PC], pfJet_phi[JC]) > 2.5){
							    CutFlowNumber[14]++;
							    CutExpFlowNumber[14] += EvtWeight;
							}


							//---------------------- For JES
							pfjet4vec_JESup.SetPxPyPzE(pfJet_px[JC]*(1+pfJet_jecUncer[JC]), pfJet_py[JC]*(1+pfJet_jecUncer[JC]), pfJet_pz[JC]*(1+pfJet_jecUncer[JC]), pfJet_E[JC]*(1+pfJet_jecUncer[JC] ));
							pfjet4vec_JESdown.SetPxPyPzE(pfJet_px[JC]*(1-pfJet_jecUncer[JC]), pfJet_py[JC]*(1-pfJet_jecUncer[JC]), pfJet_pz[JC]*(1-pfJet_jecUncer[JC]), pfJet_E[JC]*(1-pfJet_jecUncer[JC] ));
							gj4vec_JESup.SetPxPyPzE(Photon_px[PC]+pfjet4vec_JESup[0], Photon_py[PC]+pfjet4vec_JESup[1], Photon_pz[PC]+pfjet4vec_JESup[2], Photon_E[PC]+pfjet4vec_JESup[3]);
							gj4vec_JESdown.SetPxPyPzE(Photon_px[PC]+pfjet4vec_JESdown[0], Photon_py[PC]+pfjet4vec_JESdown[1], Photon_pz[PC]+pfjet4vec_JESdown[2], Photon_E[PC]+pfjet4vec_JESdown[3]);
							mass_JESup = pow((gj4vec_JESup[3]*gj4vec_JESup[3] - gj4vec_JESup[0]*gj4vec_JESup[0] - gj4vec_JESup[1]*gj4vec_JESup[1] - gj4vec_JESup[2]*gj4vec_JESup[2] ),0.5); ;
							mass_JESdown = pow((gj4vec_JESdown[3]*gj4vec_JESdown[3] - gj4vec_JESdown[0]*gj4vec_JESdown[0] - gj4vec_JESdown[1]*gj4vec_JESdown[1] - gj4vec_JESdown[2]*gj4vec_JESdown[2] ),0.5); ;
						
							h_prof_JES[0]->Fill(getMass(PC,JC),((mass_JESup - getMass(PC,JC))/getMass(PC,JC)));
							h_prof_JES[1]->Fill(getMass(PC,JC),((mass_JESdown - getMass(PC,JC))/getMass(PC,JC)));
							h_mass_JES[0]->Fill(mass_JESup,EvtWeight);
							h_mass_JES[1]->Fill(mass_JESdown,EvtWeight);

							//-----------------------------------------------------------------------

							h_Vertices[2]->Fill(goodVertex, pu_weight);
							h_nIsoPhoton[1]->Fill(foundPhoton.size(),EvtWeight);
							h_nPhoton[1]->Fill(Photon_n,EvtWeight);
							for(Int_t i_p = 0 ; i_p < foundPhoton.size() ; ++i_p){
							    h_nIsoPhotonPt[1]->Fill(Photon_pt[foundPhoton[i_p]],i_p+1,EvtWeight);
							}
							h_nIsoJet[1]->Fill(foundJet.size(),EvtWeight);
							h_nJet[1]->Fill(pfJet_n,EvtWeight);
							for(Int_t i_j = 0 ; i_j < foundJet.size() ; ++i_j){
							    h_nIsoJetPt[1]->Fill(pfJet_pt[foundJet[i_j]],i_j+1,EvtWeight);
							}
							if(foundJet.size() > 1){
							    h_ptSecondPhoton[1]->Fill(Photon_pt[foundPhoton[1]],EvtWeight);
							    h_etaSecondPhoton[1]->Fill(Photon_sc_eta[foundPhoton[1]],EvtWeight);
							}
							h_PC->Fill(PC,EvtWeight);
							h_JC->Fill(JC,EvtWeight);

							// For Centrality Ratio
							if(InnerRadius)  h_InnerRadius_InvMass->Fill(getMass(PC,JC),EvtWeight);
							if(OuterRadius)  h_OuterRadius_InvMass->Fill(getMass(PC,JC),EvtWeight);

							// Basic Plots with PileUp
							h_ptPhoton[7]->Fill(Photon_pt[PC],EvtWeight);
							h_etaPhoton[7]->Fill(Photon_sc_eta[PC],EvtWeight);
							h_ptJet[7]->Fill(pfJet_pt[JC],EvtWeight);
							h_etaJet[7]->Fill(pfJet_eta[JC],EvtWeight);
							h_mass_VarBin[7]->Fill(getMass(PC,JC),EvtWeight);
							h_mass_bin25[7]->Fill(getMass(PC,JC),EvtWeight);
							h_ptPFMet[7]->Fill(PFMetPt[0],EvtWeight);
							h_SumEtPFMet[7]->Fill(PFMetSumEt[0],EvtWeight);

							h_ptPhoton_massVarBin->Fill(Photon_pt[PC],EvtWeight);
							h_ptPhoton_VarBin->Fill(Photon_pt[PC],EvtWeight);
							h_ptJet_massVarBin->Fill(pfJet_pt[JC],EvtWeight);
							h_ptJet_VarBin->Fill(pfJet_pt[JC],EvtWeight);

							h_mass_bin1[2]->Fill(getMass(PC,JC),EvtWeight);
							h_mass_bin40[2]->Fill(getMass(PC,JC),EvtWeight);
							h_Photon_SigmaIetaIeta[2]->Fill(Photon_SigmaIetaIeta[PC],EvtWeight);
							h_PtPhotJet[2]->Fill(Photon_pt[PC],pfJet_pt[JC],EvtWeight);
							h_dEtadPhi[2]->Fill(getDEta(Photon_sc_eta[PC],pfJet_eta[JC]),getDPhi(Photon_phi[PC],pfJet_phi[JC]),EvtWeight);
							h_DR_PhotonJet[2]->Fill(getDR(Photon_sc_eta[PC], pfJet_eta[JC],Photon_phi[PC], pfJet_phi[JC]),EvtWeight);
							h_YPhotJet[2]->Fill(getRapidity(Photon_E[PC],Photon_pz[PC]),getRapidity(pfJet_E[JC],pfJet_pz[JC]),EvtWeight);
							h_dY[2]->Fill(getRapidity(Photon_E[PC],Photon_pz[PC])-getRapidity(pfJet_E[JC],pfJet_pz[JC]),EvtWeight);
							h_etaPhotJet[2]->Fill(Photon_sc_eta[PC],pfJet_eta[JC],EvtWeight);
							h_dEta[2]->Fill(getDEta(Photon_sc_eta[PC],pfJet_eta[JC]),EvtWeight);
							h_phiPhotJet[2]->Fill(Photon_phi[PC],pfJet_phi[JC],EvtWeight);
							h_dphi[2]->Fill(getDPhi(Photon_phi[PC],pfJet_phi[JC]),EvtWeight);

							// Basic Plots withOUT PileUp
							h_ptPhoton[8]->Fill(Photon_pt[PC],OnlyEvtWeight);
							h_etaPhoton[8]->Fill(Photon_sc_eta[PC],OnlyEvtWeight);
							h_ptJet[8]->Fill(pfJet_pt[JC],OnlyEvtWeight);
							h_etaJet[8]->Fill(pfJet_eta[JC],OnlyEvtWeight);
							h_mass_VarBin[8]->Fill(getMass(PC,JC),OnlyEvtWeight);
							h_mass_bin25[8]->Fill(getMass(PC,JC),OnlyEvtWeight);
							h_ptPFMet[8]->Fill(PFMetPt[0],OnlyEvtWeight);
							h_SumEtPFMet[8]->Fill(PFMetSumEt[0],OnlyEvtWeight);

							//Photon Shape variables
							h_Photon_SigmaEtaEta[1]->Fill(Photon_SigmaEtaEta[PC],EvtWeight);
							h_Photon_SigmaPhiPhi[1]->Fill(Photon_SigmaPhiPhi[PC],EvtWeight);
							h_Photon_SigmaIphiIphi[1]->Fill(Photon_SigmaIphiIphi[PC],EvtWeight);
							//Photon distinguishing variables
							h_Photon_swissCross[1]->Fill(Photon_swissCross[PC],EvtWeight);
							h_Photon_r9[1]->Fill(Photonr9[PC],EvtWeight);

							//Photon PF Isolation
							h_PFiso_ChargedBarrel[1]->Fill(PFiso_Charged03[PC],EvtWeight);
							h_PFiso_PhotonBarrel[1]->Fill(PFiso_Photon03[PC],EvtWeight);
							h_PFiso_NeutralBarrel[1]->Fill(PFiso_Neutral03[PC],EvtWeight);
							h_PFiso_SumBarrel[1]->Fill(PFiso_Sum03[PC],EvtWeight);
							h_PFiso_ElectronvetoBarrel[1]->Fill(Photon_Electronveto[PC],EvtWeight);
							h_HoEnewBarrel[1]->Fill(Photon_HoEnew[PC],EvtWeight);

							h_CorrPFiso_Photon->Fill((TMath::Max(((PFiso_Photon03[PC])  - rho*EAElectronphoton(Photon_sc_eta[PC])),0.0)),EvtWeight);
							h_CorrPFiso_Charged->Fill((TMath::Max(((PFiso_Charged03[PC]) - rho*EAElectroncharged(Photon_sc_eta[PC])),0.0)),EvtWeight);
							h_CorrPFiso_Neutral->Fill((TMath::Max(((PFiso_Neutral03[PC]) - rho*EAElectronneutral(Photon_sc_eta[PC])),0.0)),EvtWeight);

							// Photon detector based Isolation variables
							h_Photon_ecalRecHitSumEtConeDR04[1]->Fill(Photon_ecalRecHitSumEtConeDR04[PC],EvtWeight);
							h_Photon_hcalTowerSumEtConeDR04[1]->Fill(Photon_hcalTowerSumEtConeDR04[PC],EvtWeight);
							h_Photon_trkSumPtHollowConeDR04[1]->Fill(Photon_trkSumPtHollowConeDR04[PC],EvtWeight);
							h_HoE[1]->Fill(Photon_HoE[PC],EvtWeight);

							//Jet ID
							h_pfjet_NEF[1]->Fill(pfjet_NEF[JC],EvtWeight);
							h_pfjet_NHF[1]->Fill(pfjet_NHF[JC],EvtWeight);
							h_pfjet_CEF[1]->Fill(pfjet_CEF[JC],EvtWeight);
							h_pfjet_CHF[1]->Fill(pfjet_CHF[JC],EvtWeight);
							h_pfjet_NConstituents[1]->Fill(pfjet_NConstituents[JC],EvtWeight);
							h_pfjet_ChargeMultiplicity[1]->Fill(pfjet_NCH[JC],EvtWeight);
						    
						    
						    
						    } //MassCut
						}//GJDeltaEta
					    }//GJDeltaPhi
					}//JetEtaEB
				    }//JetPt
				}//JC
			    }//PhoEtaEB
			}//PhoPt
		    }//PC
		}//PV
	    }//nonScraping
	}//passHLT

    }// Loop over Events


    //Fill Histograms
    for(Int_t ij = 1; ij <= nbins; ij++){
	h_CutFlowTable->SetBinContent(ij,CutFlowNumber[ij-1]);
	h_CutFlowTable->GetXaxis()->SetBinLabel(ij,CutFlowLabel[ij-1]);
	h_CutExpFlowTable->SetBinContent(ij,CutExpFlowNumber[ij-1]);
	h_CutExpFlowTable->GetXaxis()->SetBinLabel(ij,CutFlowLabel[ij-1]);
    }
}
EOF



################################################################

cat>PostAnalyzerMC.h<<EOF
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  4 04:59:31 2012 by ROOT version 5.32/00
// from TTree myEvent/a tree with histograms
// found on file: /eos/uscms/store/user/varun/Summer12/DiJet/QCD_120to170/AOD_Output_QCD_120to170_113_1_0OS.root
//////////////////////////////////////////////////////////

#ifndef PostAnalyzerMC_h
#define PostAnalyzerMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
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
#include <TProfile.h>
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
#ifdef __MAKECINT__ 
#pragma link C++ class vector<bool>+;
#endif 

using namespace std;
using namespace ROOT;

// Fixed size dimensions of array or collections stored in the TTree if any.

class PostAnalyzerMC {
    public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	// Variables defined by me

	TFile *f1;

	Bool_t MC, DATA;
	Double_t EvtWeight, OnlyEvtWeight;
	Double_t Lumi;

	std::vector<int> foundPhoton;
	std::vector<int> foundJet;

	Float_t Cvertex_z;
	Float_t Cvertex_ndof;
	Float_t Cvertex_rho;

	Double_t photon_pt_cut;
	Double_t photon_eta_cut;
	Double_t jet_pt_cut;
	Double_t jet_eta_cut;
	Double_t mass_cut;
	Double_t PhotJetDPhi_cut;
	Double_t PhotJetDEta_cut;

	TLorentzVector pfjet4vec_JESup, pfjet4vec_JESdown ;
	TLorentzVector gj4vec_JESup, gj4vec_JESdown ;
	Double_t mass_JESup, mass_JESdown ;

	Int_t goodVertex;
	Int_t PC;
	Int_t JC;

	// Histos Declarations
	TH1F *h_CutFlowTable, *h_CutExpFlowTable;
	TH1F *h_PC, *h_JC;
	
	TH1F *h_mass_JES[2];
	TProfile *h_prof_JES[2];

	TH1F *h_ptPhoton_massVarBin, *h_ptPhoton_VarBin;
	TH1F *h_ptJet_massVarBin, *h_ptJet_VarBin;

	TH1F *h_ptPhoton[9];
	TH1F *h_etaPhoton[9];
	TH1F *h_ptJet[9];
	TH1F *h_etaJet[9];
	TH1F *h_mass_VarBin[9];
	TH1F *h_mass_bin25[9];
	TH1F *h_ptPFMet[9];
	TH1F *h_SumEtPFMet[9];

	TH1F *h_mass_bin1[3];
	TH1F *h_mass_bin40[3];
	TH1F *h_Photon_SigmaIetaIeta[3];
	TH2F *h_PtPhotJet[3];
	TH2F *h_dEtadPhi[3];
	TH1F *h_DR_PhotonJet[3];
	TH2F *h_YPhotJet[3];
	TH1F *h_dY[3];
	TH2F *h_etaPhotJet[3];
	TH1F *h_dEta[3];
	TH2F *h_phiPhotJet[3];
	TH1F *h_dphi[3];
	
	TH1F *h_ptSecondPhoton[2];
	TH1F *h_etaSecondPhoton[2];

	TH1F *h_Photon_SigmaEtaEta[2];
	TH1F *h_Photon_SigmaPhiPhi[2];
	TH1F *h_Photon_SigmaIphiIphi[2];
	TH1F *h_Photon_swissCross[2];
	TH1F *h_Photon_r9[2];

	TH1F *h_PFiso_ChargedBarrel[2];
	TH1F *h_PFiso_PhotonBarrel[2];
	TH1F *h_PFiso_NeutralBarrel[2];
	TH1F *h_PFiso_SumBarrel[2];
	TH1F *h_PFiso_ElectronvetoBarrel[2];
	TH1F *h_HoEnewBarrel[2];

	TH1F *h_CorrPFiso_Photon;
	TH1F *h_CorrPFiso_Charged;
	TH1F *h_CorrPFiso_Neutral;

	TH1F *h_Photon_ecalRecHitSumEtConeDR04[2];
	TH1F *h_Photon_hcalTowerSumEtConeDR04[2];
	TH1F *h_Photon_trkSumPtHollowConeDR04[2];
	TH1F *h_HoE[2];

	TH1F *h_pfjet_NEF[2];
	TH1F *h_pfjet_NHF[2];
	TH1F *h_pfjet_CEF[2];
	TH1F *h_pfjet_CHF[2];
	TH1F *h_pfjet_NConstituents[2];
	TH1F *h_pfjet_ChargeMultiplicity[2];

	TH1F *h_mcPileUp60;
	TH1F *h_mcPileUp600;

	TH1F *h_nPhoton[2], *h_nIsoPhoton[2];
	TH1F *h_nJet[2], *h_nIsoJet[2];
	TH2F *h_nIsoPhotonPt[2], *h_nIsoJetPt[2];

	TH1F *h_Vertices[3];
	TH1F *h_trueinteractions_MC;
	TH1F *h_trueinteractions_MC_PU;

	TH1F *weights ;
	TH1F *MC_distr_ ;
	TH1F *Data_distr_ ;
	TH1F *den ;

	TH1F *h_InnerRadius_InvMass, *h_OuterRadius_InvMass;

	
	// Declaration of leaf types
	Int_t           nevents;
	UInt_t          run;
	UInt_t          event;
	UInt_t          luminosityBlock;
	UInt_t          beamCrossing;
	UInt_t          totalIntensityBeam1;
	UInt_t          totalIntensityBeam2;
	Float_t         avgInsDelLumi;
	Float_t         avgInsDelLumiErr;
	Float_t         avgInsRecLumi;
	Float_t         avgInsRecLumiErr;
	Int_t           ntriggers;
	vector<string>  *triggernames;
	vector<int>     *triggerprescales;
	vector<bool>    *ifTriggerpassed;
	vector<float>   *ObjectPt;
	vector<float>   *ObjectEta;
	vector<float>   *ObjectPhi;
	vector<string>  *FilterNames;
	vector<int>     *FilterStartPosition;
	vector<int>     *FilterEndPosition;
	vector<int>     *ObjectStartPosition;
	vector<int>     *ObjectEndPosition;

	// Vertex variables   
	Int_t           Vertex_n;
	Float_t         Vertex_x[200];   //[Vertex_n]
	Float_t         Vertex_y[200];   //[Vertex_n]
	Float_t         Vertex_z[200];   //[Vertex_n]
	Int_t           Vertex_tracksize[200];   //[Vertex_n]
	Int_t           Vertex_ndof[200];   //[Vertex_n]
	Float_t         Vertex_chi2[200];   //[Vertex_n]
	Float_t         Vertex_d0[200];   //[Vertex_n]
	Bool_t          Vertex_isFake[200];   //[Vertex_n]

	//scraping variables 
	Bool_t          Scraping_isScrapingEvent;
	Int_t           Scraping_numOfTracks;
	Float_t         Scraping_fractionOfGoodTracks;

	//pile-up information
	Int_t           npuVertices;
	Int_t           npuVerticesp1;
	Int_t           npuVerticesm1;
	Int_t           ootnpuVertices;
	Float_t         trueInteractions;

	// Track variables
	Int_t           Track_n;
	Float_t         Track_px[1000];   //[Track_n]
	Float_t         Track_py[1000];   //[Track_n]
	Float_t         Track_pz[1000];   //[Track_n]
	Float_t         Track_vx[1000];   //[Track_n]
	Float_t         Track_vy[1000];   //[Track_n]
	Float_t         Track_vz[1000];   //[Track_n]
	Float_t         Track_pt[1000];   //[Track_n]
	Float_t         Track_eta[1000];   //[Track_n]
	Float_t         Track_phi[1000];   //[Track_n]
/*
	//Jet Variables
	Int_t           Jet_n;
	Float_t         Jet_px[200];   //[Jet_n]
	Float_t         Jet_py[200];   //[Jet_n]
	Float_t         Jet_E[200];   //[Jet_n]
	Float_t         Jet_pz[200];   //[Jet_n]
	Float_t         Jet_vx[200];   //[Jet_n]
	Float_t         Jet_vy[200];   //[Jet_n]
	Float_t         Jet_vz[200];   //[Jet_n]
	Float_t         Jet_pt[200];   //[Jet_n]
	Float_t         Jet_eta[200];   //[Jet_n]
	Float_t         Jet_phi[200];   //[Jet_n]
	Float_t         Jet_emEnergyFraction[200];   //[Jet_n]
	Float_t         Jet_energyFractionHadronic[200];   //[Jet_n]
	Int_t           Jet_hitsInN90[200];   //[Jet_n]
	Int_t           Jet_n90Hits[200];   //[Jet_n]
	Int_t           Jet_nTowers[200];   //[Jet_n]
	Float_t         Jet_fHPD[200];   //[Jet_n]
	Float_t         Jet_fRBX[200];   //[Jet_n]
	Float_t         Jet_RHF[200];   //[Jet_n]
	Float_t         Jet_jecUncer[200];   //[Jet_n]
	Float_t         Jet_jecCorr[200];   //[Jet_n]

	// Uncorrected Jets
	Float_t         ucJet_px[200];   //[Jet_n]
	Float_t         ucJet_py[200];   //[Jet_n]
	Float_t         ucJet_E[200];   //[Jet_n]
	Float_t         ucJet_pz[200];   //[Jet_n]
	Float_t         ucJet_pt[200];   //[Jet_n]
	Float_t         ucJet_eta[200];   //[Jet_n]
	Float_t         ucJet_phi[200];   //[Jet_n]
*/
	//PFJets
	Int_t           pfJet_n;
	Float_t         pfJet_px[200];   //[pfJet_n]
	Float_t         pfJet_py[200];   //[pfJet_n]
	Float_t         pfJet_E[200];   //[pfJet_n]
	Float_t         pfJet_pz[200];   //[pfJet_n]
	Float_t         pfJet_vx[200];   //[pfJet_n]
	Float_t         pfJet_vy[200];   //[pfJet_n]
	Float_t         pfJet_vz[200];   //[pfJet_n]
	Float_t         pfJet_pt[200];   //[pfJet_n]
	Float_t         pfJet_eta[200];   //[pfJet_n]
	Float_t         pfJet_phi[200];   //[pfJet_n]
	Float_t         pfjet_CEF[200];   //[pfJet_n]
	Float_t         pfjet_CHF[200];   //[pfJet_n]
	Float_t         pfjet_NEF[200];   //[pfJet_n]
	Float_t         pfjet_NHF[200];   //[pfJet_n]
	Int_t           pfjet_NCH[200];   //[pfJet_n]
	Float_t         pfjet_HFHAE[200];   //[pfJet_n]
	Float_t         pfjet_HFEME[200];   //[pfJet_n]
	Int_t           pfjet_NConstituents[200];   //[pfJet_n]
	Int_t           pfJet_partonFlavor[200];   //[pfJet_n]
	Int_t           pfJet_partonStatus[200];   //[pfJet_n]

	//PU based Jet Id
	Float_t         pujetIdFull_mva[200];   //[pfJet_n]
	Float_t         pujetIdSimple_mva[200];   //[pfJet_n]
	Float_t         pujetIdCutBased_mva[200];   //[pfJet_n]

	Int_t           pujetIdFull_loose[200];   //[pfJet_n]
	Int_t           pujetIdFull_medium[200];   //[pfJet_n]
	Int_t           pujetIdFull_tight[200];   //[pfJet_n]

	Int_t           pujetIdSimple_loose[200];   //[pfJet_n]
	Int_t           pujetIdSimple_medium[200];   //[pfJet_n]
	Int_t           pujetIdSimple_tight[200];   //[pfJet_n]

	Int_t           pujetIdCutBased_loose[200];   //[pfJet_n]
	Int_t           pujetIdCutBased_medium[200];   //[pfJet_n]
	Int_t           pujetIdCutBased_tight[200];   //[pfJet_n]

	Float_t         pfjet_TrackCountHiEffBJetTags[200];   //[pfJet_n]
	Float_t         pfjet_TrackCountHiPurBJetTags[200];   //[pfJet_n]
	Float_t         pfjet_SimpleSVHiEffBJetTags[200];   //[pfJet_n]
	Float_t         pfjet_SimpleSVHiPurBJetTags[200];   //[pfJet_n]
	Float_t         pfJet_jecUncer[200];   //[pfJet_n]
	Float_t         pfJet_jecCorr[200];   //[pfJet_n]

	// Some uncorrectd jet information
	Float_t         ucpfJet_px[200];   //[pfJet_n]
	Float_t         ucpfJet_py[200];   //[pfJet_n]
	Float_t         ucpfJet_E[200];   //[pfJet_n]
	Float_t         ucpfJet_pz[200];   //[pfJet_n]
	Float_t         ucpfJet_pt[200];   //[pfJet_n]
	Float_t         ucpfJet_eta[200];   //[pfJet_n]
	Float_t         ucpfJet_phi[200];   //[pfJet_n]

	//Electron Information
	Int_t           Electron_n;
	Float_t         Electron_px[200];   //[Electron_n]
	Float_t         Electron_py[200];   //[Electron_n]
	Float_t         Electron_pz[200];   //[Electron_n]
	Float_t         Electron_vx[200];   //[Electron_n]
	Float_t         Electron_vy[200];   //[Electron_n]
	Float_t         Electron_vz[200];   //[Electron_n]
	Float_t         Electron_pt[200];   //[Electron_n]
	Float_t         Electron_eta[200];   //[Electron_n]
	Float_t         Electron_phi[200];   //[Electron_n]
	Float_t         Electron_energy[200];   //[Electron_n]
	Float_t         Electron_charge[200];   //[Electron_n]
	Float_t         Electron_trkIso[200];   //[Electron_n]
	Float_t         Electron_ecalIso[200];   //[Electron_n]
	Float_t         Electron_hcalIso[200];   //[Electron_n]
	Float_t         Electron_SigmaIetaIeta[200];   //[Electron_n]
	Float_t         Electron_dEtaIn[200];   //[Electron_n]
	Float_t         Electron_dPhiIn[200];   //[Electron_n]
	Float_t         Electron_HoE[200];   //[Electron_n]
	Float_t         Electron_sc_energy[200];   //[Electron_n]
	Float_t         Electron_sc_eta[200];   //[Electron_n]
	Float_t         Electron_sc_phi[200];   //[Electron_n]

	//Muon Information
	Int_t           Muon_n;
	Float_t         Muon_px[200];   //[Muon_n]
	Float_t         Muon_py[200];   //[Muon_n]
	Float_t         Muon_pz[200];   //[Muon_n]
	Float_t         Muon_vx[200];   //[Muon_n]
	Float_t         Muon_vy[200];   //[Muon_n]
	Float_t         Muon_vz[200];   //[Muon_n]
	Float_t         Muon_pt[200];   //[Muon_n]
	Float_t         Muon_eta[200];   //[Muon_n]
	Float_t         Muon_phi[200];   //[Muon_n]
	Float_t         Muon_energy[200];   //[Muon_n]
	Float_t         Muon_charge[200];   //[Muon_n]
	Bool_t          Muon_isGlobalMuon[200];   //[Muon_n]
	Bool_t          Muon_isTrackerMuon[200];   //[Muon_n]
	Bool_t          Muon_isStandAloneMuon[200];   //[Muon_n]
	Bool_t          Muon_InnerTrack_isNonnull[200];   //[Muon_n]
	Bool_t          Muon_OuterTrack_isNonnull[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_x[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_y[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_z[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_px[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_py[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_InnerPoint_pz[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_OuterPoint_x[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_OuterPoint_y[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_OuterPoint_z[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_OuterPoint_px[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_OuterPoint_py[200];   //[Muon_n]
	Float_t         Muon_OuterTrack_OuterPoint_pz[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_x[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_y[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_z[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_px[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_py[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_InnerPoint_pz[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_x[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_y[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_z[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_px[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_py[200];   //[Muon_n]
	Float_t         Muon_InnerTrack_OuterPoint_pz[200];   //[Muon_n]
	Float_t         Muon_trackIso[200];   //[Muon_n]
	Float_t         Muon_ecalIso[200];   //[Muon_n]
	Float_t         Muon_hcalIso[200];   //[Muon_n]
	Float_t         Muon_relIso[200];   //[Muon_n]
	Int_t           Muon_normChi2[200];   //[Muon_n]
	Int_t           Muon_validHits[200];   //[Muon_n]
	Int_t           Muon_tkHits[200];   //[Muon_n]
	Int_t           Muon_pixHits[200];   //[Muon_n]
	Int_t           Muon_numberOfMatches[200];   //[Muon_n]
	Float_t         Muon_OuterPoint_x[200];   //[Muon_n]
	Float_t         Muon_OuterPoint_y[200];   //[Muon_n]
	Float_t         Muon_OuterPoint_z[200];   //[Muon_n]
	Float_t         Muon_InnerPoint_x[200];   //[Muon_n]
	Float_t         Muon_InnerPoint_y[200];   //[Muon_n]
	Float_t         Muon_InnerPoint_z[200];   //[Muon_n]

/*	//AK5GenJet variables
	Int_t           genJet_n;
	Float_t         genJet_px[200];   //[genJet_n]
	Float_t         genJet_py[200];   //[genJet_n]
	Float_t         genJet_E[200];   //[genJet_n]
	Float_t         genJet_pz[200];   //[genJet_n]
	Float_t         genJet_pt[200];   //[genJet_n]
	Float_t         genJet_eta[200];   //[genJet_n]
	Float_t         genJet_phi[200];   //[genJet_n]
*/
	//Cosmic Muons
	Int_t           CosmicMuon_n;
	Float_t         CosmicMuon_px[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_py[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_pz[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_pt[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_eta[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_phi[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_energy[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_charge[200];   //[CosmicMuon_n]
	Bool_t          CosmicMuon_isGlobalMuon[200];   //[CosmicMuon_n]
	Bool_t          CosmicMuon_isTrackerMuon[200];   //[CosmicMuon_n]
	Bool_t          CosmicMuon_isStandAloneMuon[200];   //[CosmicMuon_n]
	Bool_t          CosmicMuon_InnerTrack_isNonnull[200];   //[CosmicMuon_n]
	Bool_t          CosmicMuon_OuterTrack_isNonnull[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_InnerPoint_x[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_InnerPoint_y[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_InnerPoint_z[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_InnerPoint_px[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_InnerPoint_py[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_InnerPoint_pz[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_OuterPoint_x[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_OuterPoint_y[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_OuterPoint_z[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_OuterPoint_px[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_OuterPoint_py[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterTrack_OuterPoint_pz[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_x[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_y[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_z[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_px[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_py[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_InnerPoint_pz[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_OuterPoint_x[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_OuterPoint_y[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_OuterPoint_z[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_OuterPoint_px[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_OuterPoint_py[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_InnerTrack_OuterPoint_pz[200];   //[CosmicMuon_n]

	//for  AOD only
	Float_t         CosmicMuon_OuterPoint_x[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterPoint_y[200];   //[CosmicMuon_n]
	Float_t         CosmicMuon_OuterPoint_z[200];   //[CosmicMuon_n]
/*
	//Gen Photon variables
	Float_t         gen_pthat;
	Int_t           ngenphotons;
	Float_t         gen_photonpt[1000];   //[ngenphotons]
	Float_t         gen_photoneta[1000];   //[ngenphotons]
	Float_t         gen_photonphi[1000];   //[ngenphotons]
	Float_t         gen_photonpx[1000];   //[ngenphotons]
	Float_t         gen_photonpy[1000];   //[ngenphotons]
	Float_t         gen_photonpz[1000];   //[ngenphotons]
	Float_t         gen_photonE[1000];   //[ngenphotons]
	Int_t           gen_photonstatus[1000];   //[ngenphotons]
	Int_t           gen_photonMotherID[1000];   //[ngenphotons]
	Float_t         gen_photonMotherPt[1000];   //[ngenphotons]
	Float_t         gen_photonMotherEta[1000];   //[ngenphotons]
	Float_t         gen_photonMotherPhi[1000];   //[ngenphotons]
	Int_t           gen_photonMotherStatus[1000];   //[ngenphotons]
	Int_t           gen_photonGrandmotherID[1000];   //[ngenphotons]
	Float_t         gen_photonGrandmotherPt[1000];   //[ngenphotons]
	Float_t         gen_photonGrandmotherEta[1000];   //[ngenphotons]
	Float_t         gen_photonGrandmotherPhi[1000];   //[ngenphotons]
	Int_t           gen_photonGrandmotherStatus[1000];   //[ngenphotons]

	Int_t           nhardphotons;
	Float_t         gen_hardphotonpt[2];   //[nhardphotons]
	Float_t         gen_hardphotoneta[2];   //[nhardphotons]
	Float_t         gen_hardphotonphi[2];   //[nhardphotons]
	Float_t         gen_hardphotonpx[2];   //[nhardphotons]
	Float_t         gen_hardphotonpy[2];   //[nhardphotons]
	Float_t         gen_hardphotonpz[2];   //[nhardphotons]
	Float_t         gen_hardphotonE[2];   //[nhardphotons]

	//Gen level Graviton info
	Float_t         gen_gravitonpt;
	Float_t         gen_gravitoneta;
	Float_t         gen_gravitonphi;
	Float_t         gen_gravitonpx;
	Float_t         gen_gravitonpy;
	Float_t         gen_gravitonpz;
	Float_t         gen_gravitonE;

	//Gen level Wdaughter info
	Float_t         gen_Wdaughterpt[2];
	Float_t         gen_Wdaughtereta[2];
	Float_t         gen_Wdaughterphi[2];
	Float_t         gen_Wdaughterpx[2];
	Float_t         gen_Wdaughterpy[2];
	Float_t         gen_Wdaughterpz[2];
	Float_t         gen_WdaughterE[2];
	Int_t           gen_Wdaughter_charge[2];
	Int_t           gen_WdaughterID[2];

	//Gen level variable for W
	Float_t         gen_Wbosonpt;
	Float_t         gen_Wbosoneta;
	Float_t         gen_Wbosonphi;
	Float_t         gen_Wbosonpx;
	Float_t         gen_Wbosonpy;
	Float_t         gen_Wbosonpz;
	Float_t         gen_WbosonE;
	Int_t           gen_Wbosoncharge;
	Int_t           gen_WbosonID;

	//Gen level variables for Zdaughter
	Float_t         gen_Zdaughterpt[2];
	Float_t         gen_Zdaughtereta[2];
	Float_t         gen_Zdaughterphi[2];
	Float_t         gen_Zdaughterpx[2];
	Float_t         gen_Zdaughterpy[2];
	Float_t         gen_Zdaughterpz[2];
	Float_t         gen_ZdaughterE[2];
	Int_t           gen_Zdaughter_charge[2];
	Int_t           gen_ZdaughterID[2];

	//Gen level variable for Z
	Float_t         gen_Zbosonpt;
	Float_t         gen_Zbosoneta;
	Float_t         gen_Zbosonphi;
	Float_t         gen_Zbosonpx;
	Float_t         gen_Zbosonpy;
	Float_t         gen_Zbosonpz;
	Float_t         gen_ZbosonE;
	Bool_t          is_signal_event;
	Bool_t          is_Z_event;
	Bool_t          is_W_event;
	Bool_t          is_Znunu_event;
	Bool_t          is_Zelec_event;
	Bool_t          is_Zmu_event;
	Bool_t          is_Ztau_event;
	Bool_t          is_Welec_event;
	Bool_t          is_Wmu_event;
	Bool_t          is_Wtau_event;
	Bool_t          is_SingleHardPhoton_event;
	Bool_t          is_diphoton_event;
	Bool_t          is_isr_photon_event;
	Int_t           n_signal_events;
	Int_t           n_Z_events;
	Int_t           n_W_events;
	Int_t           n_Znunu_events;
	Int_t           n_Zelec_events;
	Int_t           n_Zmu_events;
	Int_t           n_Ztau_events;
	Int_t           n_Welec_events;
	Int_t           n_Wmu_events;
	Int_t           n_Wtau_events;
	Int_t           n_SingleHardPhoton_events;
	Int_t           n_diphoton_events;


	Float_t         gen_MuonID[3];
	Float_t         gen_MuonStatus[3];
	Float_t         gen_MuonPt[3];
	Float_t         gen_MuonDaughterpt[3];
	Float_t         gen_MuonDaughtereta[3];
	Float_t         gen_MuonDaughterphi[3];
	Float_t         gen_MuonDaughterpx[3];
	Float_t         gen_MuonDaughterpy[3];
	Float_t         gen_MuonDaughterpz[3];
	Float_t         gen_MuonDaughterE[3];
	Int_t           gen_MuonDaughterCharge[3];
	Int_t           gen_MuonDaughterStatus[3];
	Int_t           gen_MuonDaughterID[3];
	Float_t         gen_tauID[3];
	Float_t         gen_tauStatus[3];
	Float_t         gen_tauPt[3];
	Float_t         gen_tauDaughterpt[3];
	Float_t         gen_tauDaughtereta[3];
	Float_t         gen_tauDaughterphi[3];
	Float_t         gen_tauDaughterpx[3];
	Float_t         gen_tauDaughterpy[3];
	Float_t         gen_tauDaughterpz[3];
	Float_t         gen_tauDaughterE[3];
	Int_t           gen_tauDaughterCharge[3];
	Int_t           gen_tauDaughterStatus[3];
	Int_t           gen_tauDaughterID[3];
*/
	// Photon Information
	Int_t           Photon_n;
	Float_t         Photon_E[200];   //[Photon_n]
	Float_t         Photon_pt[200];   //[Photon_n]
	Float_t         Photon_eta[200];   //[Photon_n]
	Float_t         Photon_phi[200];   //[Photon_n]
	Float_t         Photon_theta[200];   //[Photon_n]
	Float_t         Photon_et[200];   //[Photon_n]
	Float_t         Photon_swissCross[200];   //[Photon_n]
	Float_t         Photon_e6e2[200];   //[Photon_n]
	Float_t         Photon_e4e1[200];   //[Photon_n]
	Float_t         Photonr9[200];   //[Photon_n]
	Float_t         Photon_e1x5[200];   //[Photon_n]
	Float_t         Photon_e2x5[200];   //[Photon_n]
	Float_t         Photon_e3x3[200];   //[Photon_n]
	Float_t         Photon_e5x5[200];   //[Photon_n]
	Float_t         Photon_r1x5[200];   //[Photon_n]
	Float_t         Photon_r2x5[200];   //[Photon_n]
	Float_t         Photon_maxEnergyXtal[200];   //[Photon_n]
	Float_t         Photon_SigmaEtaEta[200];   //[Photon_n]
	Float_t         Photon_SigmaIetaIeta[200];   //[Photon_n]
	Float_t         Photon_SigmaEtaPhi[200];   //[Photon_n]
	Float_t         Photon_SigmaIetaIphi[200];   //[Photon_n]
	Float_t         Photon_SigmaPhiPhi[200];   //[Photon_n]
	Float_t         Photon_SigmaIphiIphi[200];   //[Photon_n]
	Float_t         Photon_Roundness[200];   //[Photon_n]
	Float_t         Photon_Angle[200];   //[Photon_n]

	Float_t         Photon_ecalRecHitSumEtConeDR03[200];   //[Photon_n]
	Float_t         Photon_hcalTowerSumEtConeDR03[200];   //[Photon_n]
	Float_t         Photon_trkSumPtSolidConeDR03[200];   //[Photon_n]
	Float_t         Photon_trkSumPtHollowConeDR03[200];   //[Photon_n]
	Int_t           Photon_nTrkSolidConeDR03[200];   //[Photon_n]
	Int_t           Photon_nTrkHollowConeDR03[200];   //[Photon_n]
	Float_t         Photon_hcalDepth1TowerSumEtConeDR03[200];   //[Photon_n]
	Float_t         Photon_hcalDepth2TowerSumEtConeDR03[200];   //[Photon_n]
	Float_t         Photon_ecalRecHitSumEtConeDR04[200];   //[Photon_n]
	Float_t         Photon_hcalTowerSumEtConeDR04[200];   //[Photon_n]
	Float_t         Photon_trkSumPtSolidConeDR04[200];   //[Photon_n]
	Float_t         Photon_trkSumPtHollowConeDR04[200];   //[Photon_n]
	Int_t           Photon_nTrkSolidConeDR04[200];   //[Photon_n]
	Int_t           Photon_nTrkHollowConeDR04[200];   //[Photon_n]
	Float_t         Photon_hcalDepth1TowerSumEtConeDR04[200];   //[Photon_n]
	Float_t         Photon_hcalDepth2TowerSumEtConeDR04[200];   //[Photon_n]

	Bool_t          Photon_hasPixelSeed[200];   //[Photon_n]
	Bool_t          Photon_isEB[200];   //[Photon_n]
	Bool_t          Photon_isEE[200];   //[Photon_n]
	Bool_t          Photon_isEBGap[200];   //[Photon_n]
	Bool_t          Photon_isEEGap[200];   //[Photon_n]
	Bool_t          Photon_isEBEEGap[200];   //[Photon_n]
	Float_t         Photon_e2e9[200];   //[Photon_n]
	Float_t         Photon_HoE[200];   //[Photon_n]
	Float_t         Photon_HoEnew[200];   //[Photon_n]
	Float_t         Photon_px[200];   //[Photon_n]
	Float_t         Photon_py[200];   //[Photon_n]
	Float_t         Photon_pz[200];   //[Photon_n]
	Float_t         Photon_vx[200];   //[Photon_n]
	Float_t         Photon_vy[200];   //[Photon_n]
	Float_t         Photon_vz[200];   //[Photon_n]
	Int_t           Photon_no_of_basic_clusters[200];   //[Photon_n]
	Float_t         Photon_sc_energy[200];   //[Photon_n]
	Float_t         Photon_sc_eta[200];   //[Photon_n]
	Float_t         Photon_sc_phi[200];   //[Photon_n]
	Float_t         Photon_sc_x[200];   //[Photon_n]
	Float_t         Photon_sc_y[200];   //[Photon_n]
	Float_t         Photon_sc_z[200];   //[Photon_n]
	Float_t         Photon_etaWidth[200];   //[Photon_n]
	Float_t         Photon_phiWidth[200];   //[Photon_n]
	Float_t         Photon_sc_et[200];   //[Photon_n]

	// Gen matched photon variables
	Float_t         matchphotonE[200];   //[Photon_n]
	Float_t         matchphotonpt[200];   //[Photon_n]
	Float_t         matchphotoneta[200];   //[Photon_n]
	Float_t         matchphotonphi[200];   //[Photon_n]
	Float_t         matchphotonpx[200];   //[Photon_n]
	Float_t         matchphotonpy[200];   //[Photon_n]
	Float_t         matchphotonpz[200];   //[Photon_n]
	Bool_t          ismatchedphoton[200];   //[Photon_n]

	//Converted photon variables
	Bool_t          Photon_hasConvTrk[200];   //[Photon_n]
	Int_t           Photon_ntracks[200];   //[Photon_n]
	Bool_t          Photon_isconverted[200];   //[Photon_n]
	Float_t         Photon_pairInvmass[200];   //[Photon_n]
	Float_t         Photon_pairCotThetaSeperation[200];   //[Photon_n]
	Float_t         Photon_pairmomentumX[200];   //[Photon_n]
	Float_t         Photon_pairmomentumY[200];   //[Photon_n]
	Float_t         Photon_pairmomentumZ[200];   //[Photon_n]
	Float_t         Photon_EoverP[200];   //[Photon_n]
	Float_t         Photon_ConvVx[200];   //[Photon_n]
	Float_t         Photon_ConvVy[200];   //[Photon_n]
	Float_t         Photon_ConvVz[200];   //[Photon_n]
	Float_t         Photon_ZOfPrimaryVertex[200];   //[Photon_n]
	Float_t         Photon_distOfMinimumApproach[200];   //[Photon_n]
	Float_t         Photon_dPhiTracksAtVtx[200];   //[Photon_n]
	Float_t         Photon_dPhiTracksAtEcal[200];   //[Photon_n]
	Float_t         Photon_dEtaTracksAtEcal[200];   //[Photon_n]

	Int_t           npho;
	Bool_t          Photon_Electronveto[200];   //[npho]
	Float_t         PFiso_Charged03[200];   //[npho]
	Float_t         PFiso_Photon03[200];   //[npho]
	Float_t         PFiso_Neutral03[200];   //[npho]
	Float_t         PFiso_Sum03[200];   //[npho]
	Float_t         PFWorstiso_Charged03[200];   //[npho]

	Int_t           Photon_ncrys[200];   //[Photon_n]
	Float_t         Photon_timing_xtal[200][100];   //[Photon_n]
	Float_t         Photon_timingavg_xtal[200];   //[Photon_n]
	Float_t         Photon_energy_xtal[200][100];   //[Photon_n]
	Int_t           Photon_ieta_xtalEB[200][100];   //[Photon_n]
	Int_t           Photon_iphi_xtalEB[200][100];   //[Photon_n]
	Int_t           Photon_recoFlag_xtalEB[200][100];   //[Photon_n]
	Float_t         Photon_timeError_xtal[200][100];   //[Photon_n]
	Float_t         Photon_s9[200];   //[Photon_n]
/*
	//HErechit information
	Int_t           HERecHit_subset_n;
	UInt_t          HERecHit_subset_detid[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_energy[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_time[10000];   //[HERecHit_subset_n]
	Int_t           HERecHit_subset_depth[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_phi[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_eta[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_x[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_y[10000];   //[HERecHit_subset_n]
	Float_t         HERecHit_subset_z[10000];   //[HERecHit_subset_n]
*/
	//MIP variable
	Float_t         Photon_mipChi2[200];   //[Photon_n]
	Float_t         Photon_mipTotEnergy[200];   //[Photon_n]
	Float_t         Photon_mipSlope[200];   //[Photon_n]
	Float_t         Photon_mipIntercept[200];   //[Photon_n]
	Int_t           Photon_mipNhitCone[200];   //[Photon_n]
	Bool_t          Photon_mipIsHalo[200];   //[Photon_n]

	// EBrechit variables
	Int_t           EBRecHit_size;
	Float_t         EBRecHit_eta[10000];   //[EBRecHit_size]
	Float_t         EBRecHit_phi[10000];   //[EBRecHit_size]
	Int_t           EBRecHit_ieta[10000];   //[EBRecHit_size]
	Int_t           EBRecHit_iphi[10000];   //[EBRecHit_size]
	Float_t         EBRecHit_e[10000];   //[EBRecHit_size]
	Float_t         EBRecHit_et[10000];   //[EBRecHit_size]
	Int_t           EBRecHit_flag[10000];   //[EBRecHit_size]
	Float_t         EBRecHit_time[10000];   //[EBRecHit_size]

	// EErechit variables
	Int_t           EERecHit_size;
	Float_t         EERecHit_eta[10000];   //[EERecHit_size]
	Float_t         EERecHit_phi[10000];   //[EERecHit_size]
	Int_t           EERecHit_ieta[10000];   //[EERecHit_size]
	Int_t           EERecHit_iphi[10000];   //[EERecHit_size]
	Float_t         EERecHit_e[10000];   //[EERecHit_size]
	Float_t         EERecHit_et[10000];   //[EERecHit_size]
	Int_t           EERecHit_flag[10000];   //[EERecHit_size]
	Float_t         EERecHit_time[10000];   //[EERecHit_size]
/*
	//Beam Halo variables
	Bool_t          isBeamHaloGlobalLoosePass;
	Bool_t          isBeamHaloGlobalTightPass;
	Bool_t          isBeamHaloHcalLoosePass;
	Bool_t          isBeamHaloHcalTightPass;
	Bool_t          isBeamHaloCSCLoosePass;
	Bool_t          isBeamHaloCSCTightPass;
	Bool_t          isBeamHaloEcalLoosePass;
	Bool_t          isBeamHaloEcalTightPass;
	Bool_t          isBeamHaloIDTightPass;
	Bool_t          isBeamHaloIDLoosePass;
	Bool_t          isSmellsLikeHalo_Tag;
	Bool_t          isLooseHalo_Tag;
	Bool_t          isTightHalo_Tag;
	Bool_t          isExtremeTightHalo_Tag;

	//Calo MET variables
	Float_t         CaloMetSigma;
	Float_t         CaloMetEz;
	Float_t         CaloEtFractionHadronic;
	Float_t         CaloEmEtFraction;
	Float_t         CaloHadEtInHB;
	Float_t         CaloHadEtInHE;
	Float_t         CaloHadEtInHO;
	Float_t         CaloHadEtInHF;
	Float_t         CaloEmEtInEB;
	Float_t         CaloEmEtInEE;
	Float_t         CaloEmEtInHF;
	Float_t         CaloMaxEtInEmTowers;
	Float_t         CaloMaxEtInHadTowers;
	Float_t         CaloMetPt[6];
	Float_t         CaloMetPx[6];
	Float_t         CaloMetPy[6];
	Float_t         CaloMetPhi[6];
	Float_t         CaloMetSumEt[6];
	Float_t         Delta_phi;
*/
	// PFMET variables
	Float_t         PFMetPt[6];
	Float_t         PFMetPx[6];
	Float_t         PFMetPy[6];
	Float_t         PFMetPhi[6];
	Float_t         PFMetSumEt[6];
	Float_t         Delta_phiPF;

/*
	Int_t           Photon_nummoth[200];   //[Photon_n]
	Int_t           Photon_mGenpdgId[200];   //[Photon_n]
	Int_t           Photon_mGenmompdgId[200][100];   //[Photon_n]
*/
	//CaloTowers
	Int_t           CaloTower_n;
	Float_t         CaloTower_eta[5000];   //[CaloTower_n]
	Float_t         CaloTower_phi[5000];   //[CaloTower_n]
	Float_t         CaloTower_E[5000];   //[CaloTower_n]
	Float_t         CaloTower_Et[5000];   //[CaloTower_n]
	Float_t         CaloTower_emEnergy[5000];   //[CaloTower_n]
	Float_t         CaloTower_hadEnergy[5000];   //[CaloTower_n]
	Float_t         CaloTower_p[5000];   //[CaloTower_n]
	Float_t         CaloTower_EMEt[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadEt[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadPhi[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadEta[5000];   //[CaloTower_n]
	Float_t         CaloTower_EMPhi[5000];   //[CaloTower_n]
	Float_t         CaloTower_EMEta[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadX[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadY[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadZ[5000];   //[CaloTower_n]
	Float_t         CaloTower_HE_E[5000];   //[CaloTower_n]
	Float_t         CaloTower_HB_E[5000];   //[CaloTower_n]
	Float_t         CaloTower_EMTime[5000];   //[CaloTower_n]
	Float_t         CaloTower_HadTime[5000];   //[CaloTower_n]

	Float_t         rho;
	Float_t         sigma;

	Float_t         rho25;
	Float_t         sigma25;

	// List of branches
	TBranch        *b_nevents;   //!
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
	TBranch        *b_ntriggers;   //!
	TBranch        *b_triggernames;   //!
	TBranch        *b_triggerprescales;   //!
	TBranch        *b_ifTriggerpassed;   //!
	TBranch        *b_ObjectPt;   //!
	TBranch        *b_ObjectEta;   //!
	TBranch        *b_ObjectPhi;   //!
	TBranch        *b_FilterNames;   //!
	TBranch        *b_FilterStartPosition;   //!
	TBranch        *b_FilterEndPosition;   //!
	TBranch        *b_ObjectStartPosition;   //!
	TBranch        *b_ObjectEndPosition;   //!
	TBranch        *b_Vertex_n;   //!
	TBranch        *b_Vertex_x;   //!
	TBranch        *b_Vertex_y;   //!
	TBranch        *b_Vertex_z;   //!
	TBranch        *b_Vertex_tracksize;   //!
	TBranch        *b_Vertex_ndof;   //!
	TBranch        *b_Vertex_chi2;   //!
	TBranch        *b_Vertex_d0;   //!
	TBranch        *b_Vertex_isFake;   //!
	TBranch        *b_Scraping_isScrapingEvent;   //!
	TBranch        *b_Scraping_numOfTracks;   //!
	TBranch        *b_Scraping_fractionOfGoodTracks;   //!
	TBranch        *b_npuVertices;   //!
	TBranch        *b_npuVerticesp1;   //!
	TBranch        *b_npuVerticesm1;   //!
	TBranch        *b_ootnpuVertices;   //!
	TBranch        *b_trueInteractions;   //!
	TBranch        *b_Track_n;   //!
	TBranch        *b_Track_px;   //!
	TBranch        *b_Track_py;   //!
	TBranch        *b_Track_pz;   //!
	TBranch        *b_Track_vx;   //!
	TBranch        *b_Track_vy;   //!
	TBranch        *b_Track_vz;   //!
	TBranch        *b_Track_pt;   //!
	TBranch        *b_Track_eta;   //!
	TBranch        *b_Track_phi;   //!
/*
        TBranch        *b_Jet_n;   //!
	TBranch        *b_Jet_px;   //!
	TBranch        *b_Jet_py;   //!
	TBranch        *b_Jet_E;   //!
	TBranch        *b_Jet_pz;   //!
	TBranch        *b_Jet_vx;   //!
	TBranch        *b_Jet_vy;   //!
	TBranch        *b_Jet_vz;   //!
	TBranch        *b_Jet_pt;   //!
	TBranch        *b_Jet_eta;   //!
	TBranch        *b_Jet_phi;   //!
	TBranch        *b_Jet_emEnergyFraction;   //!
	TBranch        *b_Jet_energyFractionHadronic;   //!
	TBranch        *b_Jet_hitsInN90;   //!
	TBranch        *b_Jet_n90Hits;   //!
	TBranch        *b_Jet_nTowers;   //!
	TBranch        *b_Jet_fHPD;   //!
	TBranch        *b_Jet_fRBX;   //!
	TBranch        *b_Jet_RHF;   //!
	TBranch        *b_Jet_jecUncer;   //!
	TBranch        *b_Jet_jecCorr;   //!
	TBranch        *b_ucJet_px;   //!
	TBranch        *b_ucJet_py;   //!
	TBranch        *b_ucJet_E;   //!
	TBranch        *b_ucJet_pz;   //!
	TBranch        *b_ucJet_pt;   //!
	TBranch        *b_ucJet_eta;   //!
	TBranch        *b_ucJet_phi;   //!
*/	
	TBranch        *b_pfJet_n;   //!
	TBranch        *b_pfJet_px;   //!
	TBranch        *b_pfJet_py;   //!
	TBranch        *b_pfJet_E;   //!
	TBranch        *b_pfJet_pz;   //!
	TBranch        *b_pfJet_vx;   //!
	TBranch        *b_pfJet_vy;   //!
	TBranch        *b_pfJet_vz;   //!
	TBranch        *b_pfJet_pt;   //!
	TBranch        *b_pfJet_eta;   //!
	TBranch        *b_pfJet_phi;   //!
	TBranch        *b_pfjet_CEF;   //!
	TBranch        *b_pfjet_CHF;   //!
	TBranch        *b_pfjet_NEF;   //!
	TBranch        *b_pfjet_NHF;   //!
	TBranch        *b_pfjet_NCH;   //!
	TBranch        *b_pfjet_HFHAE;   //!
	TBranch        *b_pfjet_HFEME;   //!
	TBranch        *b_pfjet_NConstituents;   //!
	TBranch        *b_pfJet_partonFlavor;   //!
	TBranch        *b_pfJet_partonStatus;   //!
	TBranch        *b_pujetIdFull_mva;   //!
	TBranch        *b_pujetIdSimple_mva;   //!
	TBranch        *b_pujetIdCutBased_mva;   //!
	TBranch        *b_pujetIdFull_loose;   //!
	TBranch        *b_pujetIdFull_medium;   //!
	TBranch        *b_pujetIdFull_tight;   //!
	TBranch        *b_pujetIdSimple_loose;   //!
	TBranch        *b_pujetIdSimple_medium;   //!
	TBranch        *b_pujetIdSimple_tight;   //!
	TBranch        *b_pujetIdCutBased_loose;   //!
	TBranch        *b_pujetIdCutBased_medium;   //!
	TBranch        *b_pujetIdCutBased_tight;   //!
	TBranch        *b_pfjet_TrackCountHiEffBJetTags;   //!
	TBranch        *b_pfjet_TrackCountHiPurBJetTags;   //!
	TBranch        *b_pfjet_SimpleSVHiEffBJetTags;   //!
	TBranch        *b_pfjet_SimpleSVHiPurBJetTags;   //!
	TBranch        *b_pfJet_jecUncer;   //!
	TBranch        *b_pfJet_jecCorr;   //!
	TBranch        *b_ucpfJet_px;   //!
	TBranch        *b_ucpfJet_py;   //!
	TBranch        *b_ucpfJet_E;   //!
	TBranch        *b_ucpfJet_pz;   //!
	TBranch        *b_ucpfJet_pt;   //!
	TBranch        *b_ucpfJet_eta;   //!
	TBranch        *b_ucpfJet_phi;   //!
	TBranch        *b_Electron_n;   //!
	TBranch        *b_Electron_px;   //!
	TBranch        *b_Electron_py;   //!
	TBranch        *b_Electron_pz;   //!
	TBranch        *b_Electron_vx;   //!
	TBranch        *b_Electron_vy;   //!
	TBranch        *b_Electron_vz;   //!
	TBranch        *b_Electron_pt;   //!
	TBranch        *b_Electron_eta;   //!
	TBranch        *b_Electron_phi;   //!
	TBranch        *b_Electron_energy;   //!
	TBranch        *b_Electron_charge;   //!
	TBranch        *b_Electron_trkIso;   //!
	TBranch        *b_Electron_ecalIso;   //!
	TBranch        *b_Electron_hcalIso;   //!
	TBranch        *b_Electron_SigmaIetaIeta;   //!
	TBranch        *b_Electron_dEtaIn;   //!
	TBranch        *b_Electron_dPhiIn;   //!
	TBranch        *b_Electron_HoE;   //!
	TBranch        *b_Electron_sc_energy;   //!
	TBranch        *b_Electron_sc_eta;   //!
	TBranch        *b_Electron_sc_phi;   //!
	TBranch        *b_Muon_n;   //!
	TBranch        *b_Muon_px;   //!
	TBranch        *b_Muon_py;   //!
	TBranch        *b_Muon_pz;   //!
	TBranch        *b_Muon_vx;   //!
	TBranch        *b_Muon_vy;   //!
	TBranch        *b_Muon_vz;   //!
	TBranch        *b_Muon_pt;   //!
	TBranch        *b_Muon_eta;   //!
	TBranch        *b_Muon_phi;   //!
	TBranch        *b_Muon_energy;   //!
	TBranch        *b_Muon_charge;   //!
	TBranch        *b_Muon_isGlobalMuon;   //!
	TBranch        *b_Muon_isTrackerMuon;   //!
	TBranch        *b_Muon_isStandAloneMuon;   //!
	TBranch        *b_Muon_InnerTrack_isNonnull;   //!
	TBranch        *b_Muon_OuterTrack_isNonnull;   //!
	TBranch        *b_Muon_OuterTrack_InnerPoint_x;   //!
	TBranch        *b_Muon_OuterTrack_InnerPoint_y;   //!
	TBranch        *b_Muon_OuterTrack_InnerPoint_z;   //!
	TBranch        *b_Muon_OuterTrack_InnerPoint_px;   //!
	TBranch        *b_Muon_OuterTrack_InnerPoint_py;   //!
	TBranch        *b_Muon_OuterTrack_InnerPoint_pz;   //!
	TBranch        *b_Muon_OuterTrack_OuterPoint_x;   //!
	TBranch        *b_Muon_OuterTrack_OuterPoint_y;   //!
	TBranch        *b_Muon_OuterTrack_OuterPoint_z;   //!
	TBranch        *b_Muon_OuterTrack_OuterPoint_px;   //!
	TBranch        *b_Muon_OuterTrack_OuterPoint_py;   //!
	TBranch        *b_Muon_OuterTrack_OuterPoint_pz;   //!
	TBranch        *b_Muon_InnerTrack_InnerPoint_x;   //!
	TBranch        *b_Muon_InnerTrack_InnerPoint_y;   //!
	TBranch        *b_Muon_InnerTrack_InnerPoint_z;   //!
	TBranch        *b_Muon_InnerTrack_InnerPoint_px;   //!
	TBranch        *b_Muon_InnerTrack_InnerPoint_py;   //!
	TBranch        *b_Muon_InnerTrack_InnerPoint_pz;   //!
	TBranch        *b_Muon_InnerTrack_OuterPoint_x;   //!
	TBranch        *b_Muon_InnerTrack_OuterPoint_y;   //!
	TBranch        *b_Muon_InnerTrack_OuterPoint_z;   //!
	TBranch        *b_Muon_InnerTrack_OuterPoint_px;   //!
	TBranch        *b_Muon_InnerTrack_OuterPoint_py;   //!
	TBranch        *b_Muon_InnerTrack_OuterPoint_pz;   //!
	TBranch        *b_Muon_trackIso;   //!
	TBranch        *b_Muon_ecalIso;   //!
	TBranch        *b_Muon_hcalIso;   //!
	TBranch        *b_Muon_relIso;   //!
	TBranch        *b_Muon_normChi2;   //!
	TBranch        *b_Muon_validHits;   //!
	TBranch        *b_Muon_tkHits;   //!
	TBranch        *b_Muon_pixHits;   //!
	TBranch        *b_Muon_numberOfMatches;   //!
	TBranch        *b_Muon_OuterPoint_x;   //!
	TBranch        *b_Muon_OuterPoint_y;   //!
	TBranch        *b_Muon_OuterPoint_z;   //!
	TBranch        *b_Muon_InnerPoint_x;   //!
	TBranch        *b_Muon_InnerPoint_y;   //!
	TBranch        *b_Muon_InnerPoint_z;   //!
/*
        TBranch        *b_genJet_n;   //!
	TBranch        *b_genJet_px;   //!
	TBranch        *b_genJet_py;   //!
	TBranch        *b_genJet_E;   //!
	TBranch        *b_genJet_pz;   //!
	TBranch        *b_genJet_pt;   //!
	TBranch        *b_genJet_eta;   //!
	TBranch        *b_genJet_phi;   //!
*/	
	TBranch        *b_CosmicMuon_n;   //!
	TBranch        *b_CosmicMuon_px;   //!
	TBranch        *b_CosmicMuon_py;   //!
	TBranch        *b_CosmicMuon_pz;   //!
	TBranch        *b_CosmicMuon_pt;   //!
	TBranch        *b_CosmicMuon_eta;   //!
	TBranch        *b_CosmicMuon_phi;   //!
	TBranch        *b_CosmicMuon_energy;   //!
	TBranch        *b_CosmicMuon_charge;   //!
	TBranch        *b_CosmicMuon_isGlobalMuon;   //!
	TBranch        *b_CosmicMuon_isTrackerMuon;   //!
	TBranch        *b_CosmicMuon_isStandAloneMuon;   //!
	TBranch        *b_CosmicMuon_InnerTrack_isNonnull;   //!
	TBranch        *b_CosmicMuon_OuterTrack_isNonnull;   //!
	TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_x;   //!
	TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_y;   //!
	TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_z;   //!
	TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_px;   //!
	TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_py;   //!
	TBranch        *b_CosmicMuon_OuterTrack_InnerPoint_pz;   //!
	TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_x;   //!
	TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_y;   //!
	TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_z;   //!
	TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_px;   //!
	TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_py;   //!
	TBranch        *b_CosmicMuon_OuterTrack_OuterPoint_pz;   //!
	TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_x;   //!
	TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_y;   //!
	TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_z;   //!
	TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_px;   //!
	TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_py;   //!
	TBranch        *b_CosmicMuon_InnerTrack_InnerPoint_pz;   //!
	TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_x;   //!
	TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_y;   //!
	TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_z;   //!
	TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_px;   //!
	TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_py;   //!
	TBranch        *b_CosmicMuon_InnerTrack_OuterPoint_pz;   //!
	TBranch        *b_CosmicMuon_OuterPoint_x;   //!
	TBranch        *b_CosmicMuon_OuterPoint_y;   //!
	TBranch        *b_CosmicMuon_OuterPoint_z;   //!
/*
        TBranch        *b_gen_pthat;   //!
	TBranch        *b_ngenphotons;   //!
	TBranch        *b_gen_photonpt;   //!
	TBranch        *b_gen_photoneta;   //!
	TBranch        *b_gen_photonphi;   //!
	TBranch        *b_gen_photonpx;   //!
	TBranch        *b_gen_photonpy;   //!
	TBranch        *b_gen_photonpz;   //!
	TBranch        *b_gen_photonE;   //!
	TBranch        *b_gen_photonstatus;   //!
	TBranch        *b_gen_photonMotherID;   //!
	TBranch        *b_gen_photonMotherPt;   //!
	TBranch        *b_gen_photonMotherEta;   //!
	TBranch        *b_gen_photonMotherPhi;   //!
	TBranch        *b_gen_photonMotherStatus;   //!
	TBranch        *b_gen_photonGrandmotherID;   //!
	TBranch        *b_gen_photonGrandmotherPt;   //!
	TBranch        *b_gen_photonGrandmotherEta;   //!
	TBranch        *b_gen_photonGrandmotherPhi;   //!
	TBranch        *b_gen_photonGrandmotherStatus;   //!
	TBranch        *b_nhardphotons;   //!
	TBranch        *b_gen_hardphotonpt;   //!
	TBranch        *b_gen_hardphotoneta;   //!
	TBranch        *b_gen_hardphotonphi;   //!
	TBranch        *b_gen_hardphotonpx;   //!
	TBranch        *b_gen_hardphotonpy;   //!
	TBranch        *b_gen_hardphotonpz;   //!
	TBranch        *b_gen_hardphotonE;   //!
	TBranch        *b_gen_graviton_pt;   //!
	TBranch        *b_gen_graviton_eta;   //!
	TBranch        *b_gen_graviton_phi;   //!
	TBranch        *b_gen_graviton_px;   //!
	TBranch        *b_gen_graviton_py;   //!
	TBranch        *b_gen_graviton_pz;   //!
	TBranch        *b_gen_graviton_E;   //!
	TBranch        *b_gen_Wdaughter_pt;   //!
	TBranch        *b_gen_Wdaughter_eta;   //!
	TBranch        *b_gen_Wdaughter_phi;   //!
	TBranch        *b_gen_Wdaughter_px;   //!
	TBranch        *b_gen_Wdaughter_py;   //!
	TBranch        *b_gen_Wdaughter_pz;   //!
	TBranch        *b_gen_Wdaughter_E;   //!
	TBranch        *b_gen_Wdaughter_charge;   //!
	TBranch        *b_gen_Wdaughter_ID;   //!
	TBranch        *b_gen_Wboson_pt;   //!
	TBranch        *b_gen_Wboson_eta;   //!
	TBranch        *b_gen_Wboson_phi;   //!
	TBranch        *b_gen_Wboson_px;   //!
	TBranch        *b_gen_Wboson_py;   //!
	TBranch        *b_gen_Wboson_pz;   //!
	TBranch        *b_gen_Wboson_E;   //!
	TBranch        *b_gen_Wboson_charge;   //!
	TBranch        *b_gen_Wboson_ID;   //!
	TBranch        *b_gen_Zdaughter_pt;   //!
	TBranch        *b_gen_Zdaughter_eta;   //!
	TBranch        *b_gen_Zdaughter_phi;   //!
	TBranch        *b_gen_Zdaughter_px;   //!
	TBranch        *b_gen_Zdaughter_py;   //!
	TBranch        *b_gen_Zdaughter_pz;   //!
	TBranch        *b_gen_Zdaughter_E;   //!
	TBranch        *b_gen_Zdaughter_charge;   //!
	TBranch        *b_gen_Zdaughter_ID;   //!
	TBranch        *b_gen_Zboson_pt;   //!
	TBranch        *b_gen_Zboson_eta;   //!
	TBranch        *b_gen_Zboson_phi;   //!
	TBranch        *b_gen_Zboson_px;   //!
	TBranch        *b_gen_Zboson_py;   //!
	TBranch        *b_gen_Zboson_pz;   //!
	TBranch        *b_gen_Zboson_E;   //!
	TBranch        *b_is_signal_event;   //!
	TBranch        *b_is_Z_event;   //!
	TBranch        *b_is_W_event;   //!
	TBranch        *b_is_Znunu_event;   //!
	TBranch        *b_is_Zelec_event;   //!
	TBranch        *b_is_Zmu_event;   //!
	TBranch        *b_is_Ztau_event;   //!
	TBranch        *b_is_Welec_event;   //!
	TBranch        *b_is_Wmu_event;   //!
	TBranch        *b_is_Wtau_event;   //!
	TBranch        *b_is_SingleHardPhoton_event;   //!
	TBranch        *b_is_diphoton_event;   //!
	TBranch        *b_is_isr_photon_event;   //!
	TBranch        *b_n_signal_events;   //!
	TBranch        *b_n_Z_events;   //!
	TBranch        *b_n_W_events;   //!
	TBranch        *b_n_Znunu_events;   //!
	TBranch        *b_n_Zelec_events;   //!
	TBranch        *b_n_Zmu_events;   //!
	TBranch        *b_n_Ztau_events;   //!
	TBranch        *b_n_Welec_events;   //!
	TBranch        *b_n_Wmu_events;   //!
	TBranch        *b_n_Wtau_events;   //!
	TBranch        *b_n_SingleHardPhoton_events;   //!
	TBranch        *b_n_diphoton_events;   //!
	TBranch        *b_gen_Muon_ID;   //!
	TBranch        *b_gen_Muon_Status;   //!
	TBranch        *b_gen_Muon_Pt;   //!
	TBranch        *b_gen_MuonDaughter_pt;   //!
	TBranch        *b_gen_MuonDaughter_eta;   //!
	TBranch        *b_gen_MuonDaughter_phi;   //!
	TBranch        *b_gen_MuonDaughter_px;   //!
	TBranch        *b_gen_MuonDaughter_py;   //!
	TBranch        *b_gen_MuonDaughter_pz;   //!
	TBranch        *b_gen_MuonDaughter_E;   //!
	TBranch        *b_gen_MuonDaughter_charge;   //!
	TBranch        *b_gen_MuonDaughter_status;   //!
	TBranch        *b_gen_MuonDaughter_ID;   //!
	TBranch        *b_gen_tau_ID;   //!
	TBranch        *b_gen_tau_Status;   //!
	TBranch        *b_gen_tau_Pt;   //!
	TBranch        *b_gen_tauDaughter_pt;   //!
	TBranch        *b_gen_tauDaughter_eta;   //!
	TBranch        *b_gen_tauDaughter_phi;   //!
	TBranch        *b_gen_tauDaughter_px;   //!
	TBranch        *b_gen_tauDaughter_py;   //!
	TBranch        *b_gen_tauDaughter_pz;   //!
	TBranch        *b_gen_tauDaughter_E;   //!
	TBranch        *b_gen_tauDaughter_charge;   //!
	TBranch        *b_gen_tauDaughter_status;   //!
	TBranch        *b_gen_tauDaughter_ID;   //!
*/	
	TBranch        *b_Photon_n;   //!
	TBranch        *b_Photon_E;   //!
	TBranch        *b_Photon_pt;   //!
	TBranch        *b_Photon_eta;   //!
	TBranch        *b_Photon_phi;   //!
	TBranch        *b_Photon_theta;   //!
	TBranch        *b_Photon_et;   //!
	TBranch        *b_Photon_swissCross;   //!
	TBranch        *b_Photon_e6e2;   //!
	TBranch        *b_Photon_e4e1;   //!
	TBranch        *b_Photonr9;   //!
	TBranch        *b_Photon_e1x5;   //!
	TBranch        *b_Photon_e2x5;   //!
	TBranch        *b_Photon_e3x3;   //!
	TBranch        *b_Photon_e5x5;   //!
	TBranch        *b_Photon_r1x5;   //!
	TBranch        *b_Photon_r2x5;   //!
	TBranch        *b_Photon_maxEnergyXtal;   //!
	TBranch        *b_Photon_SigmaEtaEta;   //!
	TBranch        *b_Photon_SigmaIetaIeta;   //!
	TBranch        *b_Photon_SigmaEtaPhi;   //!
	TBranch        *b_Photon_SigmaIetaIphi;   //!
	TBranch        *b_Photon_SigmaPhiPhi;   //!
	TBranch        *b_Photon_SigmaIphiIphi;   //!
	TBranch        *b_Photon_Roundness;   //!
	TBranch        *b_Photon_Angle;   //!
	TBranch        *b_Photon_ecalRecHitSumEtConeDR03;   //!
	TBranch        *b_Photon_hcalTowerSumEtConeDR03;   //!
	TBranch        *b_Photon_trkSumPtSolidConeDR03;   //!
	TBranch        *b_Photon_trkSumPtHollowConeDR03;   //!
	TBranch        *b_Photon_nTrkSolidConeDR03;   //!
	TBranch        *b_Photon_nTrkHollowConeDR03;   //!
	TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR03;   //!
	TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR03;   //!
	TBranch        *b_Photon_ecalRecHitSumEtConeDR04;   //!
	TBranch        *b_Photon_hcalTowerSumEtConeDR04;   //!
	TBranch        *b_Photon_trkSumPtSolidConeDR04;   //!
	TBranch        *b_Photon_trkSumPtHollowConeDR04;   //!
	TBranch        *b_Photon_nTrkSolidConeDR04;   //!
	TBranch        *b_Photon_nTrkHollowConeDR04;   //!
	TBranch        *b_Photon_hcalDepth1TowerSumEtConeDR04;   //!
	TBranch        *b_Photon_hcalDepth2TowerSumEtConeDR04;   //!
	TBranch        *b_Photon_hasPixelSeed;   //!
	TBranch        *b_Photon_isEB;   //!
	TBranch        *b_Photon_isEE;   //!
	TBranch        *b_Photon_isEBGap;   //!
	TBranch        *b_Photon_isEEGap;   //!
	TBranch        *b_Photon_isEBEEGap;   //!
	TBranch        *b_Photon_e2e9;   //!
	TBranch        *b_Photon_HoE;   //!
	TBranch        *b_Photon_HoEnew;   //!
	TBranch        *b_Photon_px;   //!
	TBranch        *b_Photon_py;   //!
	TBranch        *b_Photon_pz;   //!
	TBranch        *b_Photon_vx;   //!
	TBranch        *b_Photon_vy;   //!
	TBranch        *b_Photon_vz;   //!
	TBranch        *b_Photon_no_of_basic_clusters;   //!
	TBranch        *b_Photon_sc_energy;   //!
	TBranch        *b_Photon_sc_eta;   //!
	TBranch        *b_Photon_sc_phi;   //!
	TBranch        *b_Photon_sc_x;   //!
	TBranch        *b_Photon_sc_y;   //!
	TBranch        *b_Photon_sc_z;   //!
	TBranch        *b_Photon_etaWidth;   //!
	TBranch        *b_Photon_phiWidth;   //!
	TBranch        *b_Photon_sc_et;   //!
	TBranch        *b_matchphotonE;   //!
	TBranch        *b_matchphotonpt;   //!
	TBranch        *b_matchphotoneta;   //!
	TBranch        *b_matchphotonphi;   //!
	TBranch        *b_matchphotonpx;   //!
	TBranch        *b_matchphotonpy;   //!
	TBranch        *b_matchphotonpz;   //!
	TBranch        *b_ismatchedphoton;   //!
	TBranch        *b_Photon_hasConvTrk;   //!
	TBranch        *b_Photon_ntracks;   //!
	TBranch        *b_Photon_isconverted;   //!
	TBranch        *b_Photon_pairInvmass;   //!
	TBranch        *b_Photon_pairCotThetaSeperation;   //!
	TBranch        *b_Photon_pairmomentumX;   //!
	TBranch        *b_Photon_pairmomentumY;   //!
	TBranch        *b_Photon_pairmomentumZ;   //!
	TBranch        *b_Photon_EoverP;   //!
	TBranch        *b_Photon_ConvVx;   //!
	TBranch        *b_Photon_ConvVy;   //!
	TBranch        *b_Photon_ConvVz;   //!
	TBranch        *b_Photon_ZOfPrimaryVertex;   //!
	TBranch        *b_Photon_distOfMinimumApproach;   //!
	TBranch        *b_Photon_dPhiTracksAtVtx;   //!
	TBranch        *b_Photon_dPhiTracksAtEcal;   //!
	TBranch        *b_Photon_dEtaTracksAtEcal;   //!
	TBranch        *b_npho;   //!
	TBranch        *b_Photon_Electronveto;   //!
	TBranch        *b_PFiso_Charged03;   //!
	TBranch        *b_PFiso_Photon03;   //!
	TBranch        *b_PFiso_Neutral03;   //!
	TBranch        *b_PFiso_Sum03;   //!
	TBranch        *b_PFWorstiso_Charged03;   //!
	TBranch        *b_Photon_ncrys;   //!
	TBranch        *b_Photon_timing_xtal;   //!
	TBranch        *b_Photon_timingavg_xtal;   //!
	TBranch        *b_Photon_energy_xtal;   //!
	TBranch        *b_Photon_ieta_xtalEB;   //!
	TBranch        *b_Photon_iphi_xtalEB;   //!
	TBranch        *b_Photon_recoFlag_xtalEB;   //!
	TBranch        *b_Photon_timeError_xtal;   //!
	TBranch        *b_Photon_s9;   //!
/*	
	TBranch        *b_HERecHit_subset_n;   //!
	TBranch        *b_HERecHit_subset_detid;   //!
	TBranch        *b_HERecHit_subset_energy;   //!
	TBranch        *b_HERecHit_subset_time;   //!
	TBranch        *b_HERecHit_subset_depth;   //!
	TBranch        *b_HERecHit_subset_phi;   //!
	TBranch        *b_HERecHit_subset_eta;   //!
	TBranch        *b_HERecHit_subset_x;   //!
	TBranch        *b_HERecHit_subset_y;   //!
	TBranch        *b_HERecHit_subset_z;   //!
*/	
	TBranch        *b_Photon_mipChi2;   //!
	TBranch        *b_Photon_mipTotEnergy;   //!
	TBranch        *b_Photon_mipSlope;   //!
	TBranch        *b_Photon_mipIntercept;   //!
	TBranch        *b_Photon_mipNhitCone;   //!
	TBranch        *b_Photon_mipIsHalo;   //!
	TBranch        *b_EBRecHit_size;   //!
	TBranch        *b_EBRecHit_eta;   //!
	TBranch        *b_EBRecHit_phi;   //!
	TBranch        *b_EBRecHit_ieta;   //!
	TBranch        *b_EBRecHit_iphi;   //!
	TBranch        *b_EBRecHit_e;   //!
	TBranch        *b_EBRecHit_et;   //!
	TBranch        *b_EBRecHit_flag;   //!
	TBranch        *b_EBRecHit_time;   //!
	TBranch        *b_EERecHit_size;   //!
	TBranch        *b_EERecHit_eta;   //!
	TBranch        *b_EERecHit_phi;   //!
	TBranch        *b_EERecHit_ieta;   //!
	TBranch        *b_EERecHit_iphi;   //!
	TBranch        *b_EERecHit_e;   //!
	TBranch        *b_EERecHit_et;   //!
	TBranch        *b_EERecHit_flag;   //!
	TBranch        *b_EERecHit_time;   //!
/*	
	TBranch        *b_isBeamHaloGlobalLoosePass;   //!
	TBranch        *b_isBeamHaloGloablTightPass;   //!
	TBranch        *b_isBeamHaloHcalLoosePass;   //!
	TBranch        *b_isBeamHaloHcalTightPass;   //!
	TBranch        *b_isBeamHaloCSCLoosePass;   //!
	TBranch        *b_isBeamHaloCSCTightPass;   //!
	TBranch        *b_isBeamHaloEcalLoosePass;   //!
	TBranch        *b_isBeamHaloEcalTightPass;   //!
	TBranch        *b_isBeamHaloIDTightPass;   //!
	TBranch        *b_isBeamHaloIDLoosePass;   //!
	TBranch        *b_isSmellsLikeHalo_Tag;   //!
	TBranch        *b_isLooseHalo_Tag;   //!
	TBranch        *b_isTightHalo_Tag;   //!
	TBranch        *b_isExtremeTightHalo_Tag;   //!
	TBranch        *b_CaloMetSig;   //!
	TBranch        *b_CaloMetEz;   //!
	TBranch        *b_CaloEtFractionHadronic;   //!
	TBranch        *b_CaloEmEtFraction;   //!
	TBranch        *b_CaloHadEtInHB;   //!
	TBranch        *b_CaloHadEtInHE;   //!
	TBranch        *b_CaloHadEtInHO;   //!
	TBranch        *b_CaloHadEtInHF;   //!
	TBranch        *b_CaloEmEtInEB;   //!
	TBranch        *b_CaloEmEtInEE;   //!
	TBranch        *b_CaloEmEtInHF;   //!
	TBranch        *b_CaloMaxEtInEmTowers;   //!
	TBranch        *b_CaloMaxEtInHadTowers;   //!
	TBranch        *b_CaloMetPt;   //!
	TBranch        *b_CaloMetPx;   //!
	TBranch        *b_CaloMetPy;   //!
	TBranch        *b_CaloMetPhi;   //!
	TBranch        *b_CaloMetSumEt;   //!
	TBranch        *b_Delta_phi;   //!
*/	
	TBranch        *b_PFMetPt;   //!
	TBranch        *b_PFMetPx;   //!
	TBranch        *b_PFMetPy;   //!
	TBranch        *b_PFMetPhi;   //!
	TBranch        *b_PFMetSumEt;   //!
	TBranch        *b_Delta_phiPF;   //!
/*	
	TBranch        *b_Photon_nummoth;   //!
	TBranch        *b_Photon_mGenpdgId;   //!
	TBranch        *b_Photon_mGenmompdgId;   //!
*/	
	TBranch        *b_CaloTower_n;   //!
	TBranch        *b_CaloTower_eta;   //!
	TBranch        *b_CaloTower_phi;   //!
	TBranch        *b_CaloTower_E;   //!
	TBranch        *b_CaloTower_Et;   //!
	TBranch        *b_CaloTower_emEnergy;   //!
	TBranch        *b_CaloTower_hadEnergy;   //!
	TBranch        *b_CaloTower_p;   //!
	TBranch        *b_CaloTower_EMEt;   //!
	TBranch        *b_CaloTower_HadEt;   //!
	TBranch        *b_CaloTower_HadPhi;   //!
	TBranch        *b_CaloTower_HadEta;   //!
	TBranch        *b_CaloTower_EMPhi;   //!
	TBranch        *b_CaloTower_EMEta;   //!
	TBranch        *b_CaloTower_HadX;   //!
	TBranch        *b_CaloTower_HadY;   //!
	TBranch        *b_CaloTower_HadZ;   //!
	TBranch        *b_CaloTower_HE_E;   //!
	TBranch        *b_CaloTower_HB_E;   //!
	TBranch        *b_CaloTower_EMTime;   //!
	TBranch        *b_CaloTower_HadTime;   //!
	TBranch        *b_rho;   //!
	TBranch        *b_sigma;   //!
	TBranch        *b_rho25;   //!
	TBranch        *b_sigma25;   //!

	PostAnalyzerMC();
	virtual ~PostAnalyzerMC();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);
	virtual Bool_t   NonScraping();
	virtual Bool_t   NoSpike(Int_t ipho);
	virtual Bool_t   PrimaryVertex(Int_t &goodVertex);
	virtual Double_t getRapidity(Double_t r_E, Double_t r_Pz);
	virtual Double_t getDR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2);
	virtual Double_t getDEta(Double_t eta1, Double_t eta2);
	virtual Double_t getDPhi(Double_t phit1, Double_t phi2);
	virtual Double_t getMass(Int_t pho_i, Int_t jet_i);
	virtual Double_t getLICTD(Int_t i);
	virtual Bool_t   TightPhotonPFIso(Int_t ipho);
	virtual Bool_t   MediumPhotonPFIso(Int_t ipho);
	virtual Bool_t   LoosePhotonPFIso(Int_t ipho);
	virtual Bool_t   TightJetID( Int_t ijet);
	virtual Double_t EAElectroncharged(Double_t eta);
	virtual Double_t EAElectronneutral(Double_t eta);
	virtual Double_t EAElectronphoton(Double_t eta);
	virtual Double_t EXOtightPhoID(Int_t ipho);
	virtual void     LumiReWeighting();
	virtual Double_t puweight(Float_t npv);
	virtual void     BookHistos();
};

#endif

#ifdef PostAnalyzerMC_cxx
PostAnalyzerMC::PostAnalyzerMC(){

    TChain *chain = new TChain("myEvent");

    //add input files here
    ifstream datafile;
    datafile.open("${datafile}", ifstream::in);
    char datafilename[300];

    for(Int_t ii = 1; ii <= ${ef} && ii <= ${Tot}; ++ii){
	datafile >> datafilename;
	string fname(datafilename);
	string location        = "dcap://cmsgridftp.fnal.gov:24125";
	string main_path       = "${sourceDir}";                 
	string main_short_path = "${InputFilesPath}";
	Bool_t Path_eos = ${path_eos};
	Bool_t Path_pnfs = ${path_pnfs};

	if(ii >= ${sf} && ii <= ${Tot}){
	    if(Path_pnfs){ chain->Add((location+main_short_path+fname).c_str());
		cout<<((location+main_short_path+fname).c_str())<<endl;}
		if(Path_eos){ chain->Add((main_short_path+fname).c_str());
		    cout<<((main_short_path+fname).c_str())<<endl;}

		    cout<<chain->GetEntries()<<endl;
	}
    }

/*
    TString location        = "dcap://cmsgridftp.fnal.gov:24125";
    TString main_path       = "${sourceDir}";                 
    TString main_short_path = "${InputFilesPath}";
    Bool_t Path_eos = ${path_eos};
    Bool_t Path_pnfs = ${path_pnfs};


    TSystemDirectory sourceDir("hi",main_path);
    TList* fileList = sourceDir.GetListOfFiles();
    TIter next(fileList);
    TSystemFile* fileName;
    int fileNumber = 1;
    int maxFiles = -1; 

    while ((fileName = (TSystemFile*)next()) && fileNumber > maxFiles){
	if(TString(fileName->GetName()) == "." || TString(fileName->GetName()) == ".."  ){continue;}

	TString  FullPathInputFile;

	if(Path_pnfs) FullPathInputFile = (location+main_short_path+fileName->GetName());
	if(Path_eos)  FullPathInputFile = (main_short_path+fileName->GetName());
	cout<<FullPathInputFile<<endl; 
	chain->Add(FullPathInputFile); 
	fileNumber++;  

	//cout<<chain->GetEntries()<<endl;       
    }

    //    chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/varun/2012/QstarGJ/QstarToGJ_M_1000/AOD_Output_QstarToGJ_M_1000_4_1_602.root");
    //      chain->Add("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/lpcgg/sushil/Summer12/MC/52X_MC/PhotonJet/G_Pt_120to170/Ntuples_G_Pt_120to170_396_0_Hcz.root");
*/
    Init(chain);
}

PostAnalyzerMC::~PostAnalyzerMC(){

    if (!fChain) return;
    delete fChain->GetCurrentFile();

    f1->cd();
    f1->Write();
    //tree->Write();
    //   weights->Write();
    //   den->Write();
    MC_distr_->Write();
    Data_distr_->Write();
    f1->Close();

}

Int_t PostAnalyzerMC::GetEntry(Long64_t entry){
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t PostAnalyzerMC::LoadTree(Long64_t entry){
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

void PostAnalyzerMC::Init(TTree *tree){
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    triggernames = 0;
    triggerprescales = 0;
    ifTriggerpassed = 0;
    ObjectPt = 0;
    ObjectEta = 0;
    ObjectPhi = 0;
    FilterNames = 0;
    FilterStartPosition = 0;
    FilterEndPosition = 0;
    ObjectStartPosition = 0;
    ObjectEndPosition = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("nevents", &nevents, &b_nevents);
    fChain->SetBranchAddress("run", &run, &b_RunNumber);
    fChain->SetBranchAddress("event", &event, &b_EventNumber);
    fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_LumiNumber);
    fChain->SetBranchAddress("beamCrossing", &beamCrossing, &b_BXNumber);
    fChain->SetBranchAddress("totalIntensityBeam1", &totalIntensityBeam1, &b_totalIntensityBeam1);
    fChain->SetBranchAddress("totalIntensityBeam2", &totalIntensityBeam2, &b_totalIntensityBeam2);
    fChain->SetBranchAddress("avgInsDelLumi", &avgInsDelLumi, &b_avgInsDelLumi);
    fChain->SetBranchAddress("avgInsDelLumiErr", &avgInsDelLumiErr, &b_avgInsDelLumiErr);
    fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_avgInsRecLumi);
    fChain->SetBranchAddress("avgInsRecLumiErr", &avgInsRecLumiErr, &b_avgInsRecLumiErr);
    fChain->SetBranchAddress("ntriggers", &ntriggers, &b_ntriggers);
    fChain->SetBranchAddress("triggernames", &triggernames, &b_triggernames);
    fChain->SetBranchAddress("triggerprescales", &triggerprescales, &b_triggerprescales);
    fChain->SetBranchAddress("ifTriggerpassed", &ifTriggerpassed, &b_ifTriggerpassed);
    fChain->SetBranchAddress("ObjectPt", &ObjectPt, &b_ObjectPt);
    fChain->SetBranchAddress("ObjectEta", &ObjectEta, &b_ObjectEta);
    fChain->SetBranchAddress("ObjectPhi", &ObjectPhi, &b_ObjectPhi);
    fChain->SetBranchAddress("FilterNames", &FilterNames, &b_FilterNames);
    fChain->SetBranchAddress("FilterStartPosition", &FilterStartPosition, &b_FilterStartPosition);
    fChain->SetBranchAddress("FilterEndPosition", &FilterEndPosition, &b_FilterEndPosition);
    fChain->SetBranchAddress("ObjectStartPosition", &ObjectStartPosition, &b_ObjectStartPosition);
    fChain->SetBranchAddress("ObjectEndPosition", &ObjectEndPosition, &b_ObjectEndPosition);
    fChain->SetBranchAddress("Vertex_n", &Vertex_n, &b_Vertex_n);
    fChain->SetBranchAddress("Vertex_x", Vertex_x, &b_Vertex_x);
    fChain->SetBranchAddress("Vertex_y", Vertex_y, &b_Vertex_y);
    fChain->SetBranchAddress("Vertex_z", Vertex_z, &b_Vertex_z);
    fChain->SetBranchAddress("Vertex_tracksize", Vertex_tracksize, &b_Vertex_tracksize);
    fChain->SetBranchAddress("Vertex_ndof", Vertex_ndof, &b_Vertex_ndof);
    fChain->SetBranchAddress("Vertex_chi2", Vertex_chi2, &b_Vertex_chi2);
    fChain->SetBranchAddress("Vertex_d0", Vertex_d0, &b_Vertex_d0);
    fChain->SetBranchAddress("Vertex_isFake", Vertex_isFake, &b_Vertex_isFake);
    fChain->SetBranchAddress("Scraping_isScrapingEvent", &Scraping_isScrapingEvent, &b_Scraping_isScrapingEvent);
    fChain->SetBranchAddress("Scraping_numOfTracks", &Scraping_numOfTracks, &b_Scraping_numOfTracks);
    fChain->SetBranchAddress("Scraping_fractionOfGoodTracks", &Scraping_fractionOfGoodTracks, &b_Scraping_fractionOfGoodTracks);
    fChain->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
    fChain->SetBranchAddress("npuVerticesp1", &npuVerticesp1, &b_npuVerticesp1);
    fChain->SetBranchAddress("npuVerticesm1", &npuVerticesm1, &b_npuVerticesm1);
    fChain->SetBranchAddress("ootnpuVertices", &ootnpuVertices, &b_ootnpuVertices);
    fChain->SetBranchAddress("trueInteractions", &trueInteractions, &b_trueInteractions);
    fChain->SetBranchAddress("Track_n", &Track_n, &b_Track_n);
    fChain->SetBranchAddress("Track_px", Track_px, &b_Track_px);
    fChain->SetBranchAddress("Track_py", Track_py, &b_Track_py);
    fChain->SetBranchAddress("Track_pz", Track_pz, &b_Track_pz);
    fChain->SetBranchAddress("Track_vx", Track_vx, &b_Track_vx);
    fChain->SetBranchAddress("Track_vy", Track_vy, &b_Track_vy);
    fChain->SetBranchAddress("Track_vz", Track_vz, &b_Track_vz);
    fChain->SetBranchAddress("Track_pt", Track_pt, &b_Track_pt);
    fChain->SetBranchAddress("Track_eta", Track_eta, &b_Track_eta);
    fChain->SetBranchAddress("Track_phi", Track_phi, &b_Track_phi);
/*    
    fChain->SetBranchAddress("Jet_n", &Jet_n, &b_Jet_n);
    fChain->SetBranchAddress("Jet_px", Jet_px, &b_Jet_px);
    fChain->SetBranchAddress("Jet_py", Jet_py, &b_Jet_py);
    fChain->SetBranchAddress("Jet_E", Jet_E, &b_Jet_E);
    fChain->SetBranchAddress("Jet_pz", Jet_pz, &b_Jet_pz);
    fChain->SetBranchAddress("Jet_vx", Jet_vx, &b_Jet_vx);
    fChain->SetBranchAddress("Jet_vy", Jet_vy, &b_Jet_vy);
    fChain->SetBranchAddress("Jet_vz", Jet_vz, &b_Jet_vz);
    fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
    fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
    fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
    fChain->SetBranchAddress("Jet_emEnergyFraction", Jet_emEnergyFraction, &b_Jet_emEnergyFraction);
    fChain->SetBranchAddress("Jet_energyFractionHadronic", Jet_energyFractionHadronic, &b_Jet_energyFractionHadronic);
    fChain->SetBranchAddress("Jet_hitsInN90", Jet_hitsInN90, &b_Jet_hitsInN90);
    fChain->SetBranchAddress("Jet_n90Hits", Jet_n90Hits, &b_Jet_n90Hits);
    fChain->SetBranchAddress("Jet_nTowers", Jet_nTowers, &b_Jet_nTowers);
    fChain->SetBranchAddress("Jet_fHPD", Jet_fHPD, &b_Jet_fHPD);
    fChain->SetBranchAddress("Jet_fRBX", Jet_fRBX, &b_Jet_fRBX);
    fChain->SetBranchAddress("Jet_RHF", Jet_RHF, &b_Jet_RHF);
    fChain->SetBranchAddress("Jet_jecUncer", Jet_jecUncer, &b_Jet_jecUncer);
    fChain->SetBranchAddress("Jet_jecCorr", Jet_jecCorr, &b_Jet_jecCorr);
    fChain->SetBranchAddress("ucJet_px", ucJet_px, &b_ucJet_px);
    fChain->SetBranchAddress("ucJet_py", ucJet_py, &b_ucJet_py);
    fChain->SetBranchAddress("ucJet_E", ucJet_E, &b_ucJet_E);
    fChain->SetBranchAddress("ucJet_pz", ucJet_pz, &b_ucJet_pz);
    fChain->SetBranchAddress("ucJet_pt", ucJet_pt, &b_ucJet_pt);
    fChain->SetBranchAddress("ucJet_eta", ucJet_eta, &b_ucJet_eta);
    fChain->SetBranchAddress("ucJet_phi", ucJet_phi, &b_ucJet_phi);
*/    
    fChain->SetBranchAddress("pfJet_n", &pfJet_n, &b_pfJet_n);
    fChain->SetBranchAddress("pfJet_px", pfJet_px, &b_pfJet_px);
    fChain->SetBranchAddress("pfJet_py", pfJet_py, &b_pfJet_py);
    fChain->SetBranchAddress("pfJet_E", pfJet_E, &b_pfJet_E);
    fChain->SetBranchAddress("pfJet_pz", pfJet_pz, &b_pfJet_pz);
    fChain->SetBranchAddress("pfJet_vx", pfJet_vx, &b_pfJet_vx);
    fChain->SetBranchAddress("pfJet_vy", pfJet_vy, &b_pfJet_vy);
    fChain->SetBranchAddress("pfJet_vz", pfJet_vz, &b_pfJet_vz);
    fChain->SetBranchAddress("pfJet_pt", pfJet_pt, &b_pfJet_pt);
    fChain->SetBranchAddress("pfJet_eta", pfJet_eta, &b_pfJet_eta);
    fChain->SetBranchAddress("pfJet_phi", pfJet_phi, &b_pfJet_phi);
    fChain->SetBranchAddress("pfjet_CEF", pfjet_CEF, &b_pfjet_CEF);
    fChain->SetBranchAddress("pfjet_CHF", pfjet_CHF, &b_pfjet_CHF);
    fChain->SetBranchAddress("pfjet_NEF", pfjet_NEF, &b_pfjet_NEF);
    fChain->SetBranchAddress("pfjet_NHF", pfjet_NHF, &b_pfjet_NHF);
    fChain->SetBranchAddress("pfjet_NCH", pfjet_NCH, &b_pfjet_NCH);
    fChain->SetBranchAddress("pfjet_HFHAE", pfjet_HFHAE, &b_pfjet_HFHAE);
    fChain->SetBranchAddress("pfjet_HFEME", pfjet_HFEME, &b_pfjet_HFEME);
    fChain->SetBranchAddress("pfjet_NConstituents", pfjet_NConstituents, &b_pfjet_NConstituents);
    fChain->SetBranchAddress("pfJet_partonFlavor", pfJet_partonFlavor, &b_pfJet_partonFlavor);
    fChain->SetBranchAddress("pfJet_partonStatus", pfJet_partonStatus, &b_pfJet_partonStatus);
    fChain->SetBranchAddress("pujetIdFull_mva", pujetIdFull_mva, &b_pujetIdFull_mva);
    fChain->SetBranchAddress("pujetIdSimple_mva", pujetIdSimple_mva, &b_pujetIdSimple_mva);
    fChain->SetBranchAddress("pujetIdCutBased_mva", pujetIdCutBased_mva, &b_pujetIdCutBased_mva);
    fChain->SetBranchAddress("pujetIdFull_loose", pujetIdFull_loose, &b_pujetIdFull_loose);
    fChain->SetBranchAddress("pujetIdFull_medium", pujetIdFull_medium, &b_pujetIdFull_medium);
    fChain->SetBranchAddress("pujetIdFull_tight", pujetIdFull_tight, &b_pujetIdFull_tight);
    fChain->SetBranchAddress("pujetIdSimple_loose", pujetIdSimple_loose, &b_pujetIdSimple_loose);
    fChain->SetBranchAddress("pujetIdSimple_medium", pujetIdSimple_medium, &b_pujetIdSimple_medium);
    fChain->SetBranchAddress("pujetIdSimple_tight", pujetIdSimple_tight, &b_pujetIdSimple_tight);
    fChain->SetBranchAddress("pujetIdCutBased_loose", pujetIdCutBased_loose, &b_pujetIdCutBased_loose);
    fChain->SetBranchAddress("pujetIdCutBased_medium", pujetIdCutBased_medium, &b_pujetIdCutBased_medium);
    fChain->SetBranchAddress("pujetIdCutBased_tight", pujetIdCutBased_tight, &b_pujetIdCutBased_tight);
    fChain->SetBranchAddress("pfjet_TrackCountHiEffBJetTags", pfjet_TrackCountHiEffBJetTags, &b_pfjet_TrackCountHiEffBJetTags);
    fChain->SetBranchAddress("pfjet_TrackCountHiPurBJetTags", pfjet_TrackCountHiPurBJetTags, &b_pfjet_TrackCountHiPurBJetTags);
    fChain->SetBranchAddress("pfjet_SimpleSVHiEffBJetTags", pfjet_SimpleSVHiEffBJetTags, &b_pfjet_SimpleSVHiEffBJetTags);
    fChain->SetBranchAddress("pfjet_SimpleSVHiPurBJetTags", pfjet_SimpleSVHiPurBJetTags, &b_pfjet_SimpleSVHiPurBJetTags);
    fChain->SetBranchAddress("pfJet_jecUncer", pfJet_jecUncer, &b_pfJet_jecUncer);
    fChain->SetBranchAddress("pfJet_jecCorr", pfJet_jecCorr, &b_pfJet_jecCorr);
    fChain->SetBranchAddress("ucpfJet_px", ucpfJet_px, &b_ucpfJet_px);
    fChain->SetBranchAddress("ucpfJet_py", ucpfJet_py, &b_ucpfJet_py);
    fChain->SetBranchAddress("ucpfJet_E", ucpfJet_E, &b_ucpfJet_E);
    fChain->SetBranchAddress("ucpfJet_pz", ucpfJet_pz, &b_ucpfJet_pz);
    fChain->SetBranchAddress("ucpfJet_pt", ucpfJet_pt, &b_ucpfJet_pt);
    fChain->SetBranchAddress("ucpfJet_eta", ucpfJet_eta, &b_ucpfJet_eta);
    fChain->SetBranchAddress("ucpfJet_phi", ucpfJet_phi, &b_ucpfJet_phi);
    fChain->SetBranchAddress("Electron_n", &Electron_n, &b_Electron_n);
    fChain->SetBranchAddress("Electron_px", Electron_px, &b_Electron_px);
    fChain->SetBranchAddress("Electron_py", Electron_py, &b_Electron_py);
    fChain->SetBranchAddress("Electron_pz", Electron_pz, &b_Electron_pz);
    fChain->SetBranchAddress("Electron_vx", Electron_vx, &b_Electron_vx);
    fChain->SetBranchAddress("Electron_vy", Electron_vy, &b_Electron_vy);
    fChain->SetBranchAddress("Electron_vz", Electron_vz, &b_Electron_vz);
    fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
    fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
    fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
    fChain->SetBranchAddress("Electron_energy", Electron_energy, &b_Electron_energy);
    fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
    fChain->SetBranchAddress("Electron_trkIso", Electron_trkIso, &b_Electron_trkIso);
    fChain->SetBranchAddress("Electron_ecalIso", Electron_ecalIso, &b_Electron_ecalIso);
    fChain->SetBranchAddress("Electron_hcalIso", Electron_hcalIso, &b_Electron_hcalIso);
    fChain->SetBranchAddress("Electron_SigmaIetaIeta", Electron_SigmaIetaIeta, &b_Electron_SigmaIetaIeta);
    fChain->SetBranchAddress("Electron_dEtaIn", Electron_dEtaIn, &b_Electron_dEtaIn);
    fChain->SetBranchAddress("Electron_dPhiIn", Electron_dPhiIn, &b_Electron_dPhiIn);
    fChain->SetBranchAddress("Electron_HoE", Electron_HoE, &b_Electron_HoE);
    fChain->SetBranchAddress("Electron_sc_energy", Electron_sc_energy, &b_Electron_sc_energy);
    fChain->SetBranchAddress("Electron_sc_eta", Electron_sc_eta, &b_Electron_sc_eta);
    fChain->SetBranchAddress("Electron_sc_phi", Electron_sc_phi, &b_Electron_sc_phi);
    fChain->SetBranchAddress("Muon_n", &Muon_n, &b_Muon_n);
    fChain->SetBranchAddress("Muon_px", Muon_px, &b_Muon_px);
    fChain->SetBranchAddress("Muon_py", Muon_py, &b_Muon_py);
    fChain->SetBranchAddress("Muon_pz", Muon_pz, &b_Muon_pz);
    fChain->SetBranchAddress("Muon_vx", Muon_vx, &b_Muon_vx);
    fChain->SetBranchAddress("Muon_vy", Muon_vy, &b_Muon_vy);
    fChain->SetBranchAddress("Muon_vz", Muon_vz, &b_Muon_vz);
    fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
    fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
    fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
    fChain->SetBranchAddress("Muon_energy", Muon_energy, &b_Muon_energy);
    fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
    fChain->SetBranchAddress("Muon_isGlobalMuon", Muon_isGlobalMuon, &b_Muon_isGlobalMuon);
    fChain->SetBranchAddress("Muon_isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
    fChain->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
    fChain->SetBranchAddress("Muon_InnerTrack_isNonnull", Muon_InnerTrack_isNonnull, &b_Muon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_isNonnull", Muon_OuterTrack_isNonnull, &b_Muon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_x", Muon_OuterTrack_InnerPoint_x, &b_Muon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_y", Muon_OuterTrack_InnerPoint_y, &b_Muon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_z", Muon_OuterTrack_InnerPoint_z, &b_Muon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_px", Muon_OuterTrack_InnerPoint_px, &b_Muon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_py", Muon_OuterTrack_InnerPoint_py, &b_Muon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_InnerPoint_pz", Muon_OuterTrack_InnerPoint_pz, &b_Muon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_x", Muon_OuterTrack_OuterPoint_x, &b_Muon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_y", Muon_OuterTrack_OuterPoint_y, &b_Muon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_z", Muon_OuterTrack_OuterPoint_z, &b_Muon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_px", Muon_OuterTrack_OuterPoint_px, &b_Muon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_py", Muon_OuterTrack_OuterPoint_py, &b_Muon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_OuterTrack_OuterPoint_pz", Muon_OuterTrack_OuterPoint_pz, &b_Muon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_x", Muon_InnerTrack_InnerPoint_x, &b_Muon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_y", Muon_InnerTrack_InnerPoint_y, &b_Muon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_z", Muon_InnerTrack_InnerPoint_z, &b_Muon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_px", Muon_InnerTrack_InnerPoint_px, &b_Muon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_py", Muon_InnerTrack_InnerPoint_py, &b_Muon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_InnerPoint_pz", Muon_InnerTrack_InnerPoint_pz, &b_Muon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_x", Muon_InnerTrack_OuterPoint_x, &b_Muon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_y", Muon_InnerTrack_OuterPoint_y, &b_Muon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_z", Muon_InnerTrack_OuterPoint_z, &b_Muon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_px", Muon_InnerTrack_OuterPoint_px, &b_Muon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_py", Muon_InnerTrack_OuterPoint_py, &b_Muon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("Muon_InnerTrack_OuterPoint_pz", Muon_InnerTrack_OuterPoint_pz, &b_Muon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("Muon_trackIso", Muon_trackIso, &b_Muon_trackIso);
    fChain->SetBranchAddress("Muon_ecalIso", Muon_ecalIso, &b_Muon_ecalIso);
    fChain->SetBranchAddress("Muon_hcalIso", Muon_hcalIso, &b_Muon_hcalIso);
    fChain->SetBranchAddress("Muon_relIso", Muon_relIso, &b_Muon_relIso);
    fChain->SetBranchAddress("Muon_normChi2", Muon_normChi2, &b_Muon_normChi2);
    fChain->SetBranchAddress("Muon_validHits", Muon_validHits, &b_Muon_validHits);
    fChain->SetBranchAddress("Muon_tkHits", Muon_tkHits, &b_Muon_tkHits);
    fChain->SetBranchAddress("Muon_pixHits", Muon_pixHits, &b_Muon_pixHits);
    fChain->SetBranchAddress("Muon_numberOfMatches", Muon_numberOfMatches, &b_Muon_numberOfMatches);
    fChain->SetBranchAddress("Muon_OuterPoint_x", Muon_OuterPoint_x, &b_Muon_OuterPoint_x);
    fChain->SetBranchAddress("Muon_OuterPoint_y", Muon_OuterPoint_y, &b_Muon_OuterPoint_y);
    fChain->SetBranchAddress("Muon_OuterPoint_z", Muon_OuterPoint_z, &b_Muon_OuterPoint_z);
    fChain->SetBranchAddress("Muon_InnerPoint_x", Muon_InnerPoint_x, &b_Muon_InnerPoint_x);
    fChain->SetBranchAddress("Muon_InnerPoint_y", Muon_InnerPoint_y, &b_Muon_InnerPoint_y);
    fChain->SetBranchAddress("Muon_InnerPoint_z", Muon_InnerPoint_z, &b_Muon_InnerPoint_z);
/*    
    fChain->SetBranchAddress("genJet_n", &genJet_n, &b_genJet_n);
    fChain->SetBranchAddress("genJet_px", genJet_px, &b_genJet_px);
    fChain->SetBranchAddress("genJet_py", genJet_py, &b_genJet_py);
    fChain->SetBranchAddress("genJet_E", genJet_E, &b_genJet_E);
    fChain->SetBranchAddress("genJet_pz", genJet_pz, &b_genJet_pz);
    fChain->SetBranchAddress("genJet_pt", genJet_pt, &b_genJet_pt);
    fChain->SetBranchAddress("genJet_eta", genJet_eta, &b_genJet_eta);
    fChain->SetBranchAddress("genJet_phi", genJet_phi, &b_genJet_phi);
*/    
    fChain->SetBranchAddress("CosmicMuon_n", &CosmicMuon_n, &b_CosmicMuon_n);
    fChain->SetBranchAddress("CosmicMuon_px", CosmicMuon_px, &b_CosmicMuon_px);
    fChain->SetBranchAddress("CosmicMuon_py", CosmicMuon_py, &b_CosmicMuon_py);
    fChain->SetBranchAddress("CosmicMuon_pz", CosmicMuon_pz, &b_CosmicMuon_pz);
    fChain->SetBranchAddress("CosmicMuon_pt", CosmicMuon_pt, &b_CosmicMuon_pt);
    fChain->SetBranchAddress("CosmicMuon_eta", CosmicMuon_eta, &b_CosmicMuon_eta);
    fChain->SetBranchAddress("CosmicMuon_phi", CosmicMuon_phi, &b_CosmicMuon_phi);
    fChain->SetBranchAddress("CosmicMuon_energy", CosmicMuon_energy, &b_CosmicMuon_energy);
    fChain->SetBranchAddress("CosmicMuon_charge", CosmicMuon_charge, &b_CosmicMuon_charge);
    fChain->SetBranchAddress("CosmicMuon_isGlobalMuon", CosmicMuon_isGlobalMuon, &b_CosmicMuon_isGlobalMuon);
    fChain->SetBranchAddress("CosmicMuon_isTrackerMuon", CosmicMuon_isTrackerMuon, &b_CosmicMuon_isTrackerMuon);
    fChain->SetBranchAddress("CosmicMuon_isStandAloneMuon", CosmicMuon_isStandAloneMuon, &b_CosmicMuon_isStandAloneMuon);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_isNonnull", CosmicMuon_InnerTrack_isNonnull, &b_CosmicMuon_InnerTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_isNonnull", CosmicMuon_OuterTrack_isNonnull, &b_CosmicMuon_OuterTrack_isNonnull);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_x", CosmicMuon_OuterTrack_InnerPoint_x, &b_CosmicMuon_OuterTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_y", CosmicMuon_OuterTrack_InnerPoint_y, &b_CosmicMuon_OuterTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_z", CosmicMuon_OuterTrack_InnerPoint_z, &b_CosmicMuon_OuterTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_px", CosmicMuon_OuterTrack_InnerPoint_px, &b_CosmicMuon_OuterTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_py", CosmicMuon_OuterTrack_InnerPoint_py, &b_CosmicMuon_OuterTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_InnerPoint_pz", CosmicMuon_OuterTrack_InnerPoint_pz, &b_CosmicMuon_OuterTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_x", CosmicMuon_OuterTrack_OuterPoint_x, &b_CosmicMuon_OuterTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_y", CosmicMuon_OuterTrack_OuterPoint_y, &b_CosmicMuon_OuterTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_z", CosmicMuon_OuterTrack_OuterPoint_z, &b_CosmicMuon_OuterTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_px", CosmicMuon_OuterTrack_OuterPoint_px, &b_CosmicMuon_OuterTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_py", CosmicMuon_OuterTrack_OuterPoint_py, &b_CosmicMuon_OuterTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_OuterTrack_OuterPoint_pz", CosmicMuon_OuterTrack_OuterPoint_pz, &b_CosmicMuon_OuterTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_x", CosmicMuon_InnerTrack_InnerPoint_x, &b_CosmicMuon_InnerTrack_InnerPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_y", CosmicMuon_InnerTrack_InnerPoint_y, &b_CosmicMuon_InnerTrack_InnerPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_z", CosmicMuon_InnerTrack_InnerPoint_z, &b_CosmicMuon_InnerTrack_InnerPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_px", CosmicMuon_InnerTrack_InnerPoint_px, &b_CosmicMuon_InnerTrack_InnerPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_py", CosmicMuon_InnerTrack_InnerPoint_py, &b_CosmicMuon_InnerTrack_InnerPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_InnerPoint_pz", CosmicMuon_InnerTrack_InnerPoint_pz, &b_CosmicMuon_InnerTrack_InnerPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_x", CosmicMuon_InnerTrack_OuterPoint_x, &b_CosmicMuon_InnerTrack_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_y", CosmicMuon_InnerTrack_OuterPoint_y, &b_CosmicMuon_InnerTrack_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_z", CosmicMuon_InnerTrack_OuterPoint_z, &b_CosmicMuon_InnerTrack_OuterPoint_z);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_px", CosmicMuon_InnerTrack_OuterPoint_px, &b_CosmicMuon_InnerTrack_OuterPoint_px);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_py", CosmicMuon_InnerTrack_OuterPoint_py, &b_CosmicMuon_InnerTrack_OuterPoint_py);
    fChain->SetBranchAddress("CosmicMuon_InnerTrack_OuterPoint_pz", CosmicMuon_InnerTrack_OuterPoint_pz, &b_CosmicMuon_InnerTrack_OuterPoint_pz);
    fChain->SetBranchAddress("CosmicMuon_OuterPoint_x", CosmicMuon_OuterPoint_x, &b_CosmicMuon_OuterPoint_x);
    fChain->SetBranchAddress("CosmicMuon_OuterPoint_y", CosmicMuon_OuterPoint_y, &b_CosmicMuon_OuterPoint_y);
    fChain->SetBranchAddress("CosmicMuon_OuterPoint_z", CosmicMuon_OuterPoint_z, &b_CosmicMuon_OuterPoint_z);
/*    
    fChain->SetBranchAddress("gen_pthat", &gen_pthat, &b_gen_pthat);
    fChain->SetBranchAddress("ngenphotons", &ngenphotons, &b_ngenphotons);
    fChain->SetBranchAddress("gen_photonpt", gen_photonpt, &b_gen_photonpt);
    fChain->SetBranchAddress("gen_photoneta", gen_photoneta, &b_gen_photoneta);
    fChain->SetBranchAddress("gen_photonphi", gen_photonphi, &b_gen_photonphi);
    fChain->SetBranchAddress("gen_photonpx", gen_photonpx, &b_gen_photonpx);
    fChain->SetBranchAddress("gen_photonpy", gen_photonpy, &b_gen_photonpy);
    fChain->SetBranchAddress("gen_photonpz", gen_photonpz, &b_gen_photonpz);
    fChain->SetBranchAddress("gen_photonE", gen_photonE, &b_gen_photonE);
    fChain->SetBranchAddress("gen_photonstatus", gen_photonstatus, &b_gen_photonstatus);
    fChain->SetBranchAddress("gen_photonMotherID", gen_photonMotherID, &b_gen_photonMotherID);
    fChain->SetBranchAddress("gen_photonMotherPt", gen_photonMotherPt, &b_gen_photonMotherPt);
    fChain->SetBranchAddress("gen_photonMotherEta", gen_photonMotherEta, &b_gen_photonMotherEta);
    fChain->SetBranchAddress("gen_photonMotherPhi", gen_photonMotherPhi, &b_gen_photonMotherPhi);
    fChain->SetBranchAddress("gen_photonMotherStatus", gen_photonMotherStatus, &b_gen_photonMotherStatus);
    fChain->SetBranchAddress("gen_photonGrandmotherID", gen_photonGrandmotherID, &b_gen_photonGrandmotherID);
    fChain->SetBranchAddress("gen_photonGrandmotherPt", gen_photonGrandmotherPt, &b_gen_photonGrandmotherPt);
    fChain->SetBranchAddress("gen_photonGrandmotherEta", gen_photonGrandmotherEta, &b_gen_photonGrandmotherEta);
    fChain->SetBranchAddress("gen_photonGrandmotherPhi", gen_photonGrandmotherPhi, &b_gen_photonGrandmotherPhi);
    fChain->SetBranchAddress("gen_photonGrandmotherStatus", gen_photonGrandmotherStatus, &b_gen_photonGrandmotherStatus);
    fChain->SetBranchAddress("nhardphotons", &nhardphotons, &b_nhardphotons);
    fChain->SetBranchAddress("gen_hardphotonpt", &gen_hardphotonpt, &b_gen_hardphotonpt);
    fChain->SetBranchAddress("gen_hardphotoneta", &gen_hardphotoneta, &b_gen_hardphotoneta);
    fChain->SetBranchAddress("gen_hardphotonphi", &gen_hardphotonphi, &b_gen_hardphotonphi);
    fChain->SetBranchAddress("gen_hardphotonpx", &gen_hardphotonpx, &b_gen_hardphotonpx);
    fChain->SetBranchAddress("gen_hardphotonpy", &gen_hardphotonpy, &b_gen_hardphotonpy);
    fChain->SetBranchAddress("gen_hardphotonpz", &gen_hardphotonpz, &b_gen_hardphotonpz);
    fChain->SetBranchAddress("gen_hardphotonE", &gen_hardphotonE, &b_gen_hardphotonE);
    fChain->SetBranchAddress("gen_gravitonpt", &gen_gravitonpt, &b_gen_graviton_pt);
    fChain->SetBranchAddress("gen_gravitoneta", &gen_gravitoneta, &b_gen_graviton_eta);
    fChain->SetBranchAddress("gen_gravitonphi", &gen_gravitonphi, &b_gen_graviton_phi);
    fChain->SetBranchAddress("gen_gravitonpx", &gen_gravitonpx, &b_gen_graviton_px);
    fChain->SetBranchAddress("gen_gravitonpy", &gen_gravitonpy, &b_gen_graviton_py);
    fChain->SetBranchAddress("gen_gravitonpz", &gen_gravitonpz, &b_gen_graviton_pz);
    fChain->SetBranchAddress("gen_gravitonE", &gen_gravitonE, &b_gen_graviton_E);
    fChain->SetBranchAddress("gen_Wdaughterpt", gen_Wdaughterpt, &b_gen_Wdaughter_pt);
    fChain->SetBranchAddress("gen_Wdaughtereta", gen_Wdaughtereta, &b_gen_Wdaughter_eta);
    fChain->SetBranchAddress("gen_Wdaughterphi", gen_Wdaughterphi, &b_gen_Wdaughter_phi);
    fChain->SetBranchAddress("gen_Wdaughterpx", gen_Wdaughterpx, &b_gen_Wdaughter_px);
    fChain->SetBranchAddress("gen_Wdaughterpy", gen_Wdaughterpy, &b_gen_Wdaughter_py);
    fChain->SetBranchAddress("gen_Wdaughterpz", gen_Wdaughterpz, &b_gen_Wdaughter_pz);
    fChain->SetBranchAddress("gen_WdaughterE", gen_WdaughterE, &b_gen_Wdaughter_E);
    fChain->SetBranchAddress("gen_Wdaughter_charge", gen_Wdaughter_charge, &b_gen_Wdaughter_charge);
    fChain->SetBranchAddress("gen_WdaughterID", gen_WdaughterID, &b_gen_Wdaughter_ID);
    fChain->SetBranchAddress("gen_Wbosonpt", &gen_Wbosonpt, &b_gen_Wboson_pt);
    fChain->SetBranchAddress("gen_Wbosoneta", &gen_Wbosoneta, &b_gen_Wboson_eta);
    fChain->SetBranchAddress("gen_Wbosonphi", &gen_Wbosonphi, &b_gen_Wboson_phi);
    fChain->SetBranchAddress("gen_Wbosonpx", &gen_Wbosonpx, &b_gen_Wboson_px);
    fChain->SetBranchAddress("gen_Wbosonpy", &gen_Wbosonpy, &b_gen_Wboson_py);
    fChain->SetBranchAddress("gen_Wbosonpz", &gen_Wbosonpz, &b_gen_Wboson_pz);
    fChain->SetBranchAddress("gen_WbosonE", &gen_WbosonE, &b_gen_Wboson_E);
    fChain->SetBranchAddress("gen_Wbosoncharge", &gen_Wbosoncharge, &b_gen_Wboson_charge);
    fChain->SetBranchAddress("gen_WbosonID", &gen_WbosonID, &b_gen_Wboson_ID);
    fChain->SetBranchAddress("gen_Zdaughterpt", gen_Zdaughterpt, &b_gen_Zdaughter_pt);
    fChain->SetBranchAddress("gen_Zdaughtereta", gen_Zdaughtereta, &b_gen_Zdaughter_eta);
    fChain->SetBranchAddress("gen_Zdaughterphi", gen_Zdaughterphi, &b_gen_Zdaughter_phi);
    fChain->SetBranchAddress("gen_Zdaughterpx", gen_Zdaughterpx, &b_gen_Zdaughter_px);
    fChain->SetBranchAddress("gen_Zdaughterpy", gen_Zdaughterpy, &b_gen_Zdaughter_py);
    fChain->SetBranchAddress("gen_Zdaughterpz", gen_Zdaughterpz, &b_gen_Zdaughter_pz);
    fChain->SetBranchAddress("gen_ZdaughterE", gen_ZdaughterE, &b_gen_Zdaughter_E);
    fChain->SetBranchAddress("gen_Zdaughter_charge", gen_Zdaughter_charge, &b_gen_Zdaughter_charge);
    fChain->SetBranchAddress("gen_ZdaughterID", gen_ZdaughterID, &b_gen_Zdaughter_ID);
    fChain->SetBranchAddress("gen_Zbosonpt", &gen_Zbosonpt, &b_gen_Zboson_pt);
    fChain->SetBranchAddress("gen_Zbosoneta", &gen_Zbosoneta, &b_gen_Zboson_eta);
    fChain->SetBranchAddress("gen_Zbosonphi", &gen_Zbosonphi, &b_gen_Zboson_phi);
    fChain->SetBranchAddress("gen_Zbosonpx", &gen_Zbosonpx, &b_gen_Zboson_px);
    fChain->SetBranchAddress("gen_Zbosonpy", &gen_Zbosonpy, &b_gen_Zboson_py);
    fChain->SetBranchAddress("gen_Zbosonpz", &gen_Zbosonpz, &b_gen_Zboson_pz);
    fChain->SetBranchAddress("gen_ZbosonE", &gen_ZbosonE, &b_gen_Zboson_E);
    fChain->SetBranchAddress("is_signal_event", &is_signal_event, &b_is_signal_event);
    fChain->SetBranchAddress("is_Z_event", &is_Z_event, &b_is_Z_event);
    fChain->SetBranchAddress("is_W_event", &is_W_event, &b_is_W_event);
    fChain->SetBranchAddress("is_Znunu_event", &is_Znunu_event, &b_is_Znunu_event);
    fChain->SetBranchAddress("is_Zelec_event", &is_Zelec_event, &b_is_Zelec_event);
    fChain->SetBranchAddress("is_Zmu_event", &is_Zmu_event, &b_is_Zmu_event);
    fChain->SetBranchAddress("is_Ztau_event", &is_Ztau_event, &b_is_Ztau_event);
    fChain->SetBranchAddress("is_Welec_event", &is_Welec_event, &b_is_Welec_event);
    fChain->SetBranchAddress("is_Wmu_event", &is_Wmu_event, &b_is_Wmu_event);
    fChain->SetBranchAddress("is_Wtau_event", &is_Wtau_event, &b_is_Wtau_event);
    fChain->SetBranchAddress("is_SingleHardPhoton_event", &is_SingleHardPhoton_event, &b_is_SingleHardPhoton_event);
    fChain->SetBranchAddress("is_diphoton_event", &is_diphoton_event, &b_is_diphoton_event);
    fChain->SetBranchAddress("is_isr_photon_event", &is_isr_photon_event, &b_is_isr_photon_event);
    fChain->SetBranchAddress("n_signal_events", &n_signal_events, &b_n_signal_events);
    fChain->SetBranchAddress("n_Z_events", &n_Z_events, &b_n_Z_events);
    fChain->SetBranchAddress("n_W_events", &n_W_events, &b_n_W_events);
    fChain->SetBranchAddress("n_Znunu_events", &n_Znunu_events, &b_n_Znunu_events);
    fChain->SetBranchAddress("n_Zelec_events", &n_Zelec_events, &b_n_Zelec_events);
    fChain->SetBranchAddress("n_Zmu_events", &n_Zmu_events, &b_n_Zmu_events);
    fChain->SetBranchAddress("n_Ztau_events", &n_Ztau_events, &b_n_Ztau_events);
    fChain->SetBranchAddress("n_Welec_events", &n_Welec_events, &b_n_Welec_events);
    fChain->SetBranchAddress("n_Wmu_events", &n_Wmu_events, &b_n_Wmu_events);
    fChain->SetBranchAddress("n_Wtau_events", &n_Wtau_events, &b_n_Wtau_events);
    fChain->SetBranchAddress("n_SingleHardPhoton_events", &n_SingleHardPhoton_events, &b_n_SingleHardPhoton_events);
    fChain->SetBranchAddress("n_diphoton_events", &n_diphoton_events, &b_n_diphoton_events);
    fChain->SetBranchAddress("gen_MuonID", gen_MuonID, &b_gen_Muon_ID);
    fChain->SetBranchAddress("gen_MuonStatus", gen_MuonStatus, &b_gen_Muon_Status);
    fChain->SetBranchAddress("gen_MuonPt", gen_MuonPt, &b_gen_Muon_Pt);
    fChain->SetBranchAddress("gen_MuonDaughterpt", gen_MuonDaughterpt, &b_gen_MuonDaughter_pt);
    fChain->SetBranchAddress("gen_MuonDaughtereta", gen_MuonDaughtereta, &b_gen_MuonDaughter_eta);
    fChain->SetBranchAddress("gen_MuonDaughterphi", gen_MuonDaughterphi, &b_gen_MuonDaughter_phi);
    fChain->SetBranchAddress("gen_MuonDaughterpx", gen_MuonDaughterpx, &b_gen_MuonDaughter_px);
    fChain->SetBranchAddress("gen_MuonDaughterpy", gen_MuonDaughterpy, &b_gen_MuonDaughter_py);
    fChain->SetBranchAddress("gen_MuonDaughterpz", gen_MuonDaughterpz, &b_gen_MuonDaughter_pz);
    fChain->SetBranchAddress("gen_MuonDaughterE", gen_MuonDaughterE, &b_gen_MuonDaughter_E);
    fChain->SetBranchAddress("gen_MuonDaughterCharge", gen_MuonDaughterCharge, &b_gen_MuonDaughter_charge);
    fChain->SetBranchAddress("gen_MuonDaughterStatus", gen_MuonDaughterStatus, &b_gen_MuonDaughter_status);
    fChain->SetBranchAddress("gen_MuonDaughterID", gen_MuonDaughterID, &b_gen_MuonDaughter_ID);
    fChain->SetBranchAddress("gen_tauID", gen_tauID, &b_gen_tau_ID);
    fChain->SetBranchAddress("gen_tauStatus", gen_tauStatus, &b_gen_tau_Status);
    fChain->SetBranchAddress("gen_tauPt", gen_tauPt, &b_gen_tau_Pt);
    fChain->SetBranchAddress("gen_tauDaughterpt", gen_tauDaughterpt, &b_gen_tauDaughter_pt);
    fChain->SetBranchAddress("gen_tauDaughtereta", gen_tauDaughtereta, &b_gen_tauDaughter_eta);
    fChain->SetBranchAddress("gen_tauDaughterphi", gen_tauDaughterphi, &b_gen_tauDaughter_phi);
    fChain->SetBranchAddress("gen_tauDaughterpx", gen_tauDaughterpx, &b_gen_tauDaughter_px);
    fChain->SetBranchAddress("gen_tauDaughterpy", gen_tauDaughterpy, &b_gen_tauDaughter_py);
    fChain->SetBranchAddress("gen_tauDaughterpz", gen_tauDaughterpz, &b_gen_tauDaughter_pz);
    fChain->SetBranchAddress("gen_tauDaughterE", gen_tauDaughterE, &b_gen_tauDaughter_E);
    fChain->SetBranchAddress("gen_tauDaughterCharge", gen_tauDaughterCharge, &b_gen_tauDaughter_charge);
    fChain->SetBranchAddress("gen_tauDaughterStatus", gen_tauDaughterStatus, &b_gen_tauDaughter_status);
    fChain->SetBranchAddress("gen_tauDaughterID", gen_tauDaughterID, &b_gen_tauDaughter_ID);
*/    
    fChain->SetBranchAddress("Photon_n", &Photon_n, &b_Photon_n);
    fChain->SetBranchAddress("Photon_E", Photon_E, &b_Photon_E);
    fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
    fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
    fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
    fChain->SetBranchAddress("Photon_theta", Photon_theta, &b_Photon_theta);
    fChain->SetBranchAddress("Photon_et", Photon_et, &b_Photon_et);
    fChain->SetBranchAddress("Photon_swissCross", Photon_swissCross, &b_Photon_swissCross);
    fChain->SetBranchAddress("Photon_e6e2", Photon_e6e2, &b_Photon_e6e2);
    fChain->SetBranchAddress("Photon_e4e1", Photon_e4e1, &b_Photon_e4e1);
    fChain->SetBranchAddress("Photonr9", Photonr9, &b_Photonr9);
    fChain->SetBranchAddress("Photon_e1x5", Photon_e1x5, &b_Photon_e1x5);
    fChain->SetBranchAddress("Photon_e2x5", Photon_e2x5, &b_Photon_e2x5);
    fChain->SetBranchAddress("Photon_e3x3", Photon_e3x3, &b_Photon_e3x3);
    fChain->SetBranchAddress("Photon_e5x5", Photon_e5x5, &b_Photon_e5x5);
    fChain->SetBranchAddress("Photon_r1x5", Photon_r1x5, &b_Photon_r1x5);
    fChain->SetBranchAddress("Photon_r2x5", Photon_r2x5, &b_Photon_r2x5);
    fChain->SetBranchAddress("Photon_maxEnergyXtal", Photon_maxEnergyXtal, &b_Photon_maxEnergyXtal);
    fChain->SetBranchAddress("Photon_SigmaEtaEta", Photon_SigmaEtaEta, &b_Photon_SigmaEtaEta);
    fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta, &b_Photon_SigmaIetaIeta);
    fChain->SetBranchAddress("Photon_SigmaEtaPhi", Photon_SigmaEtaPhi, &b_Photon_SigmaEtaPhi);
    fChain->SetBranchAddress("Photon_SigmaIetaIphi", Photon_SigmaIetaIphi, &b_Photon_SigmaIetaIphi);
    fChain->SetBranchAddress("Photon_SigmaPhiPhi", Photon_SigmaPhiPhi, &b_Photon_SigmaPhiPhi);
    fChain->SetBranchAddress("Photon_SigmaIphiIphi", Photon_SigmaIphiIphi, &b_Photon_SigmaIphiIphi);
    fChain->SetBranchAddress("Photon_Roundness", Photon_Roundness, &b_Photon_Roundness);
    fChain->SetBranchAddress("Photon_Angle", Photon_Angle, &b_Photon_Angle);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR03", Photon_ecalRecHitSumEtConeDR03, &b_Photon_ecalRecHitSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR03", Photon_hcalTowerSumEtConeDR03, &b_Photon_hcalTowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR03", Photon_trkSumPtSolidConeDR03, &b_Photon_trkSumPtSolidConeDR03);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR03", Photon_trkSumPtHollowConeDR03, &b_Photon_trkSumPtHollowConeDR03);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR03", Photon_nTrkSolidConeDR03, &b_Photon_nTrkSolidConeDR03);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR03", Photon_nTrkHollowConeDR03, &b_Photon_nTrkHollowConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR03", Photon_hcalDepth1TowerSumEtConeDR03, &b_Photon_hcalDepth1TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR03", Photon_hcalDepth2TowerSumEtConeDR03, &b_Photon_hcalDepth2TowerSumEtConeDR03);
    fChain->SetBranchAddress("Photon_ecalRecHitSumEtConeDR04", Photon_ecalRecHitSumEtConeDR04, &b_Photon_ecalRecHitSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", Photon_hcalTowerSumEtConeDR04, &b_Photon_hcalTowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtSolidConeDR04", Photon_trkSumPtSolidConeDR04, &b_Photon_trkSumPtSolidConeDR04);
    fChain->SetBranchAddress("Photon_trkSumPtHollowConeDR04", Photon_trkSumPtHollowConeDR04, &b_Photon_trkSumPtHollowConeDR04);
    fChain->SetBranchAddress("Photon_nTrkSolidConeDR04", Photon_nTrkSolidConeDR04, &b_Photon_nTrkSolidConeDR04);
    fChain->SetBranchAddress("Photon_nTrkHollowConeDR04", Photon_nTrkHollowConeDR04, &b_Photon_nTrkHollowConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth1TowerSumEtConeDR04", Photon_hcalDepth1TowerSumEtConeDR04, &b_Photon_hcalDepth1TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hcalDepth2TowerSumEtConeDR04", Photon_hcalDepth2TowerSumEtConeDR04, &b_Photon_hcalDepth2TowerSumEtConeDR04);
    fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
    fChain->SetBranchAddress("Photon_isEB", Photon_isEB, &b_Photon_isEB);
    fChain->SetBranchAddress("Photon_isEE", Photon_isEE, &b_Photon_isEE);
    fChain->SetBranchAddress("Photon_isEBGap", Photon_isEBGap, &b_Photon_isEBGap);
    fChain->SetBranchAddress("Photon_isEEGap", Photon_isEEGap, &b_Photon_isEEGap);
    fChain->SetBranchAddress("Photon_isEBEEGap", Photon_isEBEEGap, &b_Photon_isEBEEGap);
    fChain->SetBranchAddress("Photon_e2e9", Photon_e2e9, &b_Photon_e2e9);
    fChain->SetBranchAddress("Photon_HoE", Photon_HoE, &b_Photon_HoE);
    fChain->SetBranchAddress("Photon_HoEnew", Photon_HoEnew, &b_Photon_HoEnew);
    fChain->SetBranchAddress("Photon_px", Photon_px, &b_Photon_px);
    fChain->SetBranchAddress("Photon_py", Photon_py, &b_Photon_py);
    fChain->SetBranchAddress("Photon_pz", Photon_pz, &b_Photon_pz);
    fChain->SetBranchAddress("Photon_vx", Photon_vx, &b_Photon_vx);
    fChain->SetBranchAddress("Photon_vy", Photon_vy, &b_Photon_vy);
    fChain->SetBranchAddress("Photon_vz", Photon_vz, &b_Photon_vz);
    fChain->SetBranchAddress("Photon_no_of_basic_clusters", Photon_no_of_basic_clusters, &b_Photon_no_of_basic_clusters);
    fChain->SetBranchAddress("Photon_sc_energy", Photon_sc_energy, &b_Photon_sc_energy);
    fChain->SetBranchAddress("Photon_sc_eta", Photon_sc_eta, &b_Photon_sc_eta);
    fChain->SetBranchAddress("Photon_sc_phi", Photon_sc_phi, &b_Photon_sc_phi);
    fChain->SetBranchAddress("Photon_sc_x", Photon_sc_x, &b_Photon_sc_x);
    fChain->SetBranchAddress("Photon_sc_y", Photon_sc_y, &b_Photon_sc_y);
    fChain->SetBranchAddress("Photon_sc_z", Photon_sc_z, &b_Photon_sc_z);
    fChain->SetBranchAddress("Photon_etaWidth", Photon_etaWidth, &b_Photon_etaWidth);
    fChain->SetBranchAddress("Photon_phiWidth", Photon_phiWidth, &b_Photon_phiWidth);
    fChain->SetBranchAddress("Photon_sc_et", Photon_sc_et, &b_Photon_sc_et);
    fChain->SetBranchAddress("matchphotonE", matchphotonE, &b_matchphotonE);
    fChain->SetBranchAddress("matchphotonpt", matchphotonpt, &b_matchphotonpt);
    fChain->SetBranchAddress("matchphotoneta", matchphotoneta, &b_matchphotoneta);
    fChain->SetBranchAddress("matchphotonphi", matchphotonphi, &b_matchphotonphi);
    fChain->SetBranchAddress("matchphotonpx", matchphotonpx, &b_matchphotonpx);
    fChain->SetBranchAddress("matchphotonpy", matchphotonpy, &b_matchphotonpy);
    fChain->SetBranchAddress("matchphotonpz", matchphotonpz, &b_matchphotonpz);
    fChain->SetBranchAddress("ismatchedphoton", ismatchedphoton, &b_ismatchedphoton);
    fChain->SetBranchAddress("Photon_hasConvTrk", Photon_hasConvTrk, &b_Photon_hasConvTrk);
    fChain->SetBranchAddress("Photon_ntracks", Photon_ntracks, &b_Photon_ntracks);
    fChain->SetBranchAddress("Photon_isconverted", Photon_isconverted, &b_Photon_isconverted);
    fChain->SetBranchAddress("Photon_pairInvmass", Photon_pairInvmass, &b_Photon_pairInvmass);
    fChain->SetBranchAddress("Photon_pairCotThetaSeperation", Photon_pairCotThetaSeperation, &b_Photon_pairCotThetaSeperation);
    fChain->SetBranchAddress("Photon_pairmomentumX", Photon_pairmomentumX, &b_Photon_pairmomentumX);
    fChain->SetBranchAddress("Photon_pairmomentumY", Photon_pairmomentumY, &b_Photon_pairmomentumY);
    fChain->SetBranchAddress("Photon_pairmomentumZ", Photon_pairmomentumZ, &b_Photon_pairmomentumZ);
    fChain->SetBranchAddress("Photon_EoverP", Photon_EoverP, &b_Photon_EoverP);
    fChain->SetBranchAddress("Photon_ConvVx", Photon_ConvVx, &b_Photon_ConvVx);
    fChain->SetBranchAddress("Photon_ConvVy", Photon_ConvVy, &b_Photon_ConvVy);
    fChain->SetBranchAddress("Photon_ConvVz", Photon_ConvVz, &b_Photon_ConvVz);
    fChain->SetBranchAddress("Photon_ZOfPrimaryVertex", Photon_ZOfPrimaryVertex, &b_Photon_ZOfPrimaryVertex);
    fChain->SetBranchAddress("Photon_distOfMinimumApproach", Photon_distOfMinimumApproach, &b_Photon_distOfMinimumApproach);
    fChain->SetBranchAddress("Photon_dPhiTracksAtVtx", Photon_dPhiTracksAtVtx, &b_Photon_dPhiTracksAtVtx);
    fChain->SetBranchAddress("Photon_dPhiTracksAtEcal", Photon_dPhiTracksAtEcal, &b_Photon_dPhiTracksAtEcal);
    fChain->SetBranchAddress("Photon_dEtaTracksAtEcal", Photon_dEtaTracksAtEcal, &b_Photon_dEtaTracksAtEcal);
    fChain->SetBranchAddress("npho", &npho, &b_npho);
    fChain->SetBranchAddress("Photon_Electronveto", Photon_Electronveto, &b_Photon_Electronveto);
    fChain->SetBranchAddress("PFiso_Charged03", PFiso_Charged03, &b_PFiso_Charged03);
    fChain->SetBranchAddress("PFiso_Photon03", PFiso_Photon03, &b_PFiso_Photon03);
    fChain->SetBranchAddress("PFiso_Neutral03", PFiso_Neutral03, &b_PFiso_Neutral03);
    fChain->SetBranchAddress("PFiso_Sum03", PFiso_Sum03, &b_PFiso_Sum03);
    fChain->SetBranchAddress("PFWorstiso_Charged03", PFWorstiso_Charged03, &b_PFWorstiso_Charged03);
    fChain->SetBranchAddress("Photon_ncrys", Photon_ncrys, &b_Photon_ncrys);
    fChain->SetBranchAddress("Photon_timing_xtal", Photon_timing_xtal, &b_Photon_timing_xtal);
    fChain->SetBranchAddress("Photon_timingavg_xtal", Photon_timingavg_xtal, &b_Photon_timingavg_xtal);
    fChain->SetBranchAddress("Photon_energy_xtal", Photon_energy_xtal, &b_Photon_energy_xtal);
    fChain->SetBranchAddress("Photon_ieta_xtalEB", Photon_ieta_xtalEB, &b_Photon_ieta_xtalEB);
    fChain->SetBranchAddress("Photon_iphi_xtalEB", Photon_iphi_xtalEB, &b_Photon_iphi_xtalEB);
    fChain->SetBranchAddress("Photon_recoFlag_xtalEB", Photon_recoFlag_xtalEB, &b_Photon_recoFlag_xtalEB);
    fChain->SetBranchAddress("Photon_timeError_xtal", Photon_timeError_xtal, &b_Photon_timeError_xtal);
    fChain->SetBranchAddress("Photon_s9", Photon_s9, &b_Photon_s9);
/*    
    fChain->SetBranchAddress("HERecHit_subset_n", &HERecHit_subset_n, &b_HERecHit_subset_n);
    fChain->SetBranchAddress("HERecHit_subset_detid", HERecHit_subset_detid, &b_HERecHit_subset_detid);
    fChain->SetBranchAddress("HERecHit_subset_energy", HERecHit_subset_energy, &b_HERecHit_subset_energy);
    fChain->SetBranchAddress("HERecHit_subset_time", HERecHit_subset_time, &b_HERecHit_subset_time);
    fChain->SetBranchAddress("HERecHit_subset_depth", HERecHit_subset_depth, &b_HERecHit_subset_depth);
    fChain->SetBranchAddress("HERecHit_subset_phi", HERecHit_subset_phi, &b_HERecHit_subset_phi);
    fChain->SetBranchAddress("HERecHit_subset_eta", HERecHit_subset_eta, &b_HERecHit_subset_eta);
    fChain->SetBranchAddress("HERecHit_subset_x", HERecHit_subset_x, &b_HERecHit_subset_x);
    fChain->SetBranchAddress("HERecHit_subset_y", HERecHit_subset_y, &b_HERecHit_subset_y);
    fChain->SetBranchAddress("HERecHit_subset_z", HERecHit_subset_z, &b_HERecHit_subset_z);
*/    
    fChain->SetBranchAddress("Photon_mipChi2", Photon_mipChi2, &b_Photon_mipChi2);
    fChain->SetBranchAddress("Photon_mipTotEnergy", Photon_mipTotEnergy, &b_Photon_mipTotEnergy);
    fChain->SetBranchAddress("Photon_mipSlope", Photon_mipSlope, &b_Photon_mipSlope);
    fChain->SetBranchAddress("Photon_mipIntercept", Photon_mipIntercept, &b_Photon_mipIntercept);
    fChain->SetBranchAddress("Photon_mipNhitCone", Photon_mipNhitCone, &b_Photon_mipNhitCone);
    fChain->SetBranchAddress("Photon_mipIsHalo", Photon_mipIsHalo, &b_Photon_mipIsHalo);
    fChain->SetBranchAddress("EBRecHit_size", &EBRecHit_size, &b_EBRecHit_size);
    fChain->SetBranchAddress("EBRecHit_eta", EBRecHit_eta, &b_EBRecHit_eta);
    fChain->SetBranchAddress("EBRecHit_phi", EBRecHit_phi, &b_EBRecHit_phi);
    fChain->SetBranchAddress("EBRecHit_ieta", EBRecHit_ieta, &b_EBRecHit_ieta);
    fChain->SetBranchAddress("EBRecHit_iphi", EBRecHit_iphi, &b_EBRecHit_iphi);
    fChain->SetBranchAddress("EBRecHit_e", EBRecHit_e, &b_EBRecHit_e);
    fChain->SetBranchAddress("EBRecHit_et", EBRecHit_et, &b_EBRecHit_et);
    fChain->SetBranchAddress("EBRecHit_flag", EBRecHit_flag, &b_EBRecHit_flag);
    fChain->SetBranchAddress("EBRecHit_time", EBRecHit_time, &b_EBRecHit_time);
    fChain->SetBranchAddress("EERecHit_size", &EERecHit_size, &b_EERecHit_size);
    fChain->SetBranchAddress("EERecHit_eta", EERecHit_eta, &b_EERecHit_eta);
    fChain->SetBranchAddress("EERecHit_phi", EERecHit_phi, &b_EERecHit_phi);
    fChain->SetBranchAddress("EERecHit_ieta", EERecHit_ieta, &b_EERecHit_ieta);
    fChain->SetBranchAddress("EERecHit_iphi", EERecHit_iphi, &b_EERecHit_iphi);
    fChain->SetBranchAddress("EERecHit_e", EERecHit_e, &b_EERecHit_e);
    fChain->SetBranchAddress("EERecHit_et", EERecHit_et, &b_EERecHit_et);
    fChain->SetBranchAddress("EERecHit_flag", EERecHit_flag, &b_EERecHit_flag);
    fChain->SetBranchAddress("EERecHit_time", EERecHit_time, &b_EERecHit_time);
/*    
    fChain->SetBranchAddress("isBeamHaloGlobalLoosePass", &isBeamHaloGlobalLoosePass, &b_isBeamHaloGlobalLoosePass);
    fChain->SetBranchAddress("isBeamHaloGlobalTightPass", &isBeamHaloGlobalTightPass, &b_isBeamHaloGloablTightPass);
    fChain->SetBranchAddress("isBeamHaloHcalLoosePass", &isBeamHaloHcalLoosePass, &b_isBeamHaloHcalLoosePass);
    fChain->SetBranchAddress("isBeamHaloHcalTightPass", &isBeamHaloHcalTightPass, &b_isBeamHaloHcalTightPass);
    fChain->SetBranchAddress("isBeamHaloCSCLoosePass", &isBeamHaloCSCLoosePass, &b_isBeamHaloCSCLoosePass);
    fChain->SetBranchAddress("isBeamHaloCSCTightPass", &isBeamHaloCSCTightPass, &b_isBeamHaloCSCTightPass);
    fChain->SetBranchAddress("isBeamHaloEcalLoosePass", &isBeamHaloEcalLoosePass, &b_isBeamHaloEcalLoosePass);
    fChain->SetBranchAddress("isBeamHaloEcalTightPass", &isBeamHaloEcalTightPass, &b_isBeamHaloEcalTightPass);
    fChain->SetBranchAddress("isBeamHaloIDTightPass", &isBeamHaloIDTightPass, &b_isBeamHaloIDTightPass);
    fChain->SetBranchAddress("isBeamHaloIDLoosePass", &isBeamHaloIDLoosePass, &b_isBeamHaloIDLoosePass);
    fChain->SetBranchAddress("isSmellsLikeHalo_Tag", &isSmellsLikeHalo_Tag, &b_isSmellsLikeHalo_Tag);
    fChain->SetBranchAddress("isLooseHalo_Tag", &isLooseHalo_Tag, &b_isLooseHalo_Tag);
    fChain->SetBranchAddress("isTightHalo_Tag", &isTightHalo_Tag, &b_isTightHalo_Tag);
    fChain->SetBranchAddress("isExtremeTightHalo_Tag", &isExtremeTightHalo_Tag, &b_isExtremeTightHalo_Tag);
    fChain->SetBranchAddress("CaloMetSigma", &CaloMetSigma, &b_CaloMetSig);
    fChain->SetBranchAddress("CaloMetEz", &CaloMetEz, &b_CaloMetEz);
    fChain->SetBranchAddress("CaloEtFractionHadronic", &CaloEtFractionHadronic, &b_CaloEtFractionHadronic);
    fChain->SetBranchAddress("CaloEmEtFraction", &CaloEmEtFraction, &b_CaloEmEtFraction);
    fChain->SetBranchAddress("CaloHadEtInHB", &CaloHadEtInHB, &b_CaloHadEtInHB);
    fChain->SetBranchAddress("CaloHadEtInHE", &CaloHadEtInHE, &b_CaloHadEtInHE);
    fChain->SetBranchAddress("CaloHadEtInHO", &CaloHadEtInHO, &b_CaloHadEtInHO);
    fChain->SetBranchAddress("CaloHadEtInHF", &CaloHadEtInHF, &b_CaloHadEtInHF);
    fChain->SetBranchAddress("CaloEmEtInEB", &CaloEmEtInEB, &b_CaloEmEtInEB);
    fChain->SetBranchAddress("CaloEmEtInEE", &CaloEmEtInEE, &b_CaloEmEtInEE);
    fChain->SetBranchAddress("CaloEmEtInHF", &CaloEmEtInHF, &b_CaloEmEtInHF);
    fChain->SetBranchAddress("CaloMaxEtInEmTowers", &CaloMaxEtInEmTowers, &b_CaloMaxEtInEmTowers);
    fChain->SetBranchAddress("CaloMaxEtInHadTowers", &CaloMaxEtInHadTowers, &b_CaloMaxEtInHadTowers);
    fChain->SetBranchAddress("CaloMetPt", CaloMetPt, &b_CaloMetPt);
    fChain->SetBranchAddress("CaloMetPx", CaloMetPx, &b_CaloMetPx);
    fChain->SetBranchAddress("CaloMetPy", CaloMetPy, &b_CaloMetPy);
    fChain->SetBranchAddress("CaloMetPhi", CaloMetPhi, &b_CaloMetPhi);
    fChain->SetBranchAddress("CaloMetSumEt", CaloMetSumEt, &b_CaloMetSumEt);
    fChain->SetBranchAddress("Delta_phi", &Delta_phi, &b_Delta_phi);
*/    
    fChain->SetBranchAddress("PFMetPt", PFMetPt, &b_PFMetPt);
    fChain->SetBranchAddress("PFMetPx", PFMetPx, &b_PFMetPx);
    fChain->SetBranchAddress("PFMetPy", PFMetPy, &b_PFMetPy);
    fChain->SetBranchAddress("PFMetPhi", PFMetPhi, &b_PFMetPhi);
    fChain->SetBranchAddress("PFMetSumEt", PFMetSumEt, &b_PFMetSumEt);
    fChain->SetBranchAddress("Delta_phiPF", &Delta_phiPF, &b_Delta_phiPF);
/*    
    fChain->SetBranchAddress("Photon_nummoth", Photon_nummoth, &b_Photon_nummoth);
    fChain->SetBranchAddress("Photon_mGenpdgId", Photon_mGenpdgId, &b_Photon_mGenpdgId);
    fChain->SetBranchAddress("Photon_mGenmompdgId", Photon_mGenmompdgId, &b_Photon_mGenmompdgId);
*/    
    fChain->SetBranchAddress("CaloTower_n", &CaloTower_n, &b_CaloTower_n);
    fChain->SetBranchAddress("CaloTower_eta", CaloTower_eta, &b_CaloTower_eta);
    fChain->SetBranchAddress("CaloTower_phi", CaloTower_phi, &b_CaloTower_phi);
    fChain->SetBranchAddress("CaloTower_E", CaloTower_E, &b_CaloTower_E);
    fChain->SetBranchAddress("CaloTower_Et", CaloTower_Et, &b_CaloTower_Et);
    fChain->SetBranchAddress("CaloTower_emEnergy", CaloTower_emEnergy, &b_CaloTower_emEnergy);
    fChain->SetBranchAddress("CaloTower_hadEnergy", CaloTower_hadEnergy, &b_CaloTower_hadEnergy);
    fChain->SetBranchAddress("CaloTower_p", CaloTower_p, &b_CaloTower_p);
    fChain->SetBranchAddress("CaloTower_EMEt", CaloTower_EMEt, &b_CaloTower_EMEt);
    fChain->SetBranchAddress("CaloTower_HadEt", CaloTower_HadEt, &b_CaloTower_HadEt);
    fChain->SetBranchAddress("CaloTower_HadPhi", CaloTower_HadPhi, &b_CaloTower_HadPhi);
    fChain->SetBranchAddress("CaloTower_HadEta", CaloTower_HadEta, &b_CaloTower_HadEta);
    fChain->SetBranchAddress("CaloTower_EMPhi", CaloTower_EMPhi, &b_CaloTower_EMPhi);
    fChain->SetBranchAddress("CaloTower_EMEta", CaloTower_EMEta, &b_CaloTower_EMEta);
    fChain->SetBranchAddress("CaloTower_HadX", CaloTower_HadX, &b_CaloTower_HadX);
    fChain->SetBranchAddress("CaloTower_HadY", CaloTower_HadY, &b_CaloTower_HadY);
    fChain->SetBranchAddress("CaloTower_HadZ", CaloTower_HadZ, &b_CaloTower_HadZ);
    fChain->SetBranchAddress("CaloTower_HE_E", CaloTower_HE_E, &b_CaloTower_HE_E);
    fChain->SetBranchAddress("CaloTower_HB_E", CaloTower_HB_E, &b_CaloTower_HB_E);
    fChain->SetBranchAddress("CaloTower_EMTime", CaloTower_EMTime, &b_CaloTower_EMTime);
    fChain->SetBranchAddress("CaloTower_HadTime", CaloTower_HadTime, &b_CaloTower_HadTime);
    fChain->SetBranchAddress("rho", &rho, &b_rho);
    fChain->SetBranchAddress("sigma", &sigma, &b_sigma);
    fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
    fChain->SetBranchAddress("sigma25", &sigma25, &b_sigma25);
    Notify();
}

Bool_t PostAnalyzerMC::Notify(){
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void PostAnalyzerMC::Show(Long64_t entry){
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t PostAnalyzerMC::Cut(Long64_t entry){
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}

void PostAnalyzerMC::BookHistos(){

    // Define Histograms here +++++++++++++++++++++++++++++
    f1->cd();

    char name[100];
    const Int_t nMassBins = 119;
    const Double_t MassBin[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 73, 86, 100, 115, 132, 150, 169, 189, 210, 232, 252, 273, 295, 318, 341, 365, 390, 416, 443, 471, 500, 530, 560, 593, 626, 660, 695, 731, 768, 806, 846, 887, 929, 972, 1017, 1063, 1110, 1159, 1209, 1261, 1315, 1370, 1427, 1486, 1547, 1609, 1673, 1739, 1807, 1877, 1950, 2025, 2102, 2182, 2264, 2349, 2436, 2526, 2619, 2714, 2812, 2913, 3018, 3126, 3237, 3352, 3470, 3592, 3718, 3847, 3980, 4117, 4259, 4405, 4556, 4711, 4871, 5036, 5206, 5381, 5562, 5748, 5940, 6138, 6342, 6552, 6769, 6993, 7223, 7461, 7706, 7959, 8219, 8487, 8764, 9049, 9343, 9646, 9958, 10280, 10612, 10954, 11307, 11671, 12046, 12432, 12830, 13241, 13664, 14000};

    const Int_t nPtMassBins = 110;
    const Double_t ptMassBins[nPtMassBins+1]  = {1 ,8 ,15 ,22 ,29 ,37 ,45 ,53 ,61 ,69 ,77 ,85 ,93 ,102 ,111 ,120 ,129 ,138 ,147 ,156 ,165 ,170 ,180 ,195 ,205 ,215 ,225 ,236 ,247 ,258 ,269 ,280 ,291 ,303 ,315 ,327 ,339 ,351 ,363 ,376 ,389 ,402 ,415 ,428 ,442 ,456 ,470 ,484 ,499 ,514 ,529 ,544 ,559 ,575 ,591 ,607 ,623 ,640 ,657 ,674 ,691 ,709 ,727 ,745 ,764 ,783 ,802 ,821 ,841 ,861 ,881 ,902 ,923 ,944 ,966 ,988 ,1010 ,1033 ,1056 ,1079 ,1103 ,1127 ,1151 ,1176 ,1201 ,1227 ,1253 ,1279 ,1306 ,1333 ,1361 ,1389 ,1418 ,1447 ,1476 ,1506 ,1536 ,1567 ,1598 ,1630 ,1662 ,1695 ,1728 ,1762 ,1796 ,1831 ,1866 ,1902 ,1938 ,1975 ,2013};

    const Int_t nPtBins = 28;
    const Double_t ptBins[nPtBins+1]  = {1 ,150 , 180 , 210, 240, 270, 300, 330, 360, 390, 420, 450, 480, 510, 540, 570, 610, 650, 690, 740, 790, 850, 910, 970, 1050, 1200, 1500, 2000, 2500 };

    std::string cut[2] = {"up", "down"};
    for( Int_t hi = 0 ; hi < 2 ; ++hi){
	
	sprintf(name, "h_mass_JES%s",cut[hi].c_str());
	h_mass_JES[hi] = new TH1F(name,"mass of photon+jet",nMassBins,MassBin);
	h_mass_JES[hi]->GetYaxis()->SetTitle("Events");    h_mass_JES[hi]->GetYaxis()->CenterTitle();
	h_mass_JES[hi]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)"); h_mass_JES[hi]->GetXaxis()->CenterTitle();
	h_mass_JES[hi]->Sumw2();

	sprintf(name, "h_prof_JES%s",cut[hi].c_str());
	h_prof_JES[hi] = new TProfile(name,"Profile for JES vs mass", nMassBins,MassBin);

    }

    h_ptPhoton_massVarBin  = new TH1F("h_ptPhoton_massVarBin","pt of photon",nPtMassBins,ptMassBins);
    h_ptPhoton_massVarBin->GetYaxis()->SetTitle("Events");                  h_ptPhoton_massVarBin->GetYaxis()->CenterTitle();
    h_ptPhoton_massVarBin->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_ptPhoton_massVarBin->GetXaxis()->CenterTitle();
    h_ptPhoton_massVarBin->Sumw2();

    h_ptPhoton_VarBin  = new TH1F("h_ptPhoton_VarBin","pt of photon",nPtBins,ptBins);
    h_ptPhoton_VarBin->GetYaxis()->SetTitle("Events");                  h_ptPhoton_VarBin->GetYaxis()->CenterTitle();
    h_ptPhoton_VarBin->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_ptPhoton_VarBin->GetXaxis()->CenterTitle();
    h_ptPhoton_VarBin->Sumw2();

    h_ptJet_massVarBin = new TH1F("h_ptJet_massVarBin","pt of jet ",nPtMassBins,ptMassBins);
    h_ptJet_massVarBin->GetYaxis()->SetTitle("Events");              h_ptJet_massVarBin->GetYaxis()->CenterTitle();
    h_ptJet_massVarBin->GetXaxis()->SetTitle("P_{T}^{jet} (GeV)");   h_ptJet_massVarBin->GetXaxis()->CenterTitle();
    h_ptJet_massVarBin->Sumw2();

    h_ptJet_VarBin = new TH1F("h_ptJet_VarBin","pt of jet ",nPtBins,ptBins);
    h_ptJet_VarBin->GetYaxis()->SetTitle("Events");              h_ptJet_VarBin->GetYaxis()->CenterTitle();
    h_ptJet_VarBin->GetXaxis()->SetTitle("P_{T}^{jet} (GeV)");   h_ptJet_VarBin->GetXaxis()->CenterTitle();
    h_ptJet_VarBin->Sumw2();

    std::string cut1[9] = {"NoCut", "PhotonID","PhoPt", "PhoEta", "JetID", "JetPt", "JetEta", "finalPU", "finalnoPU"};
    for( Int_t hi = 0 ; hi < 9 ; ++hi){

	sprintf(name, "h_ptPhoton_%s",cut1[hi].c_str());
	h_ptPhoton[hi]  = new TH1F(name,"pt of photon",100,20.0,2520.0);  // As we are putting a cut of 170.
	h_ptPhoton[hi]->GetYaxis()->SetTitle("Events/25 GeV");           h_ptPhoton[hi]->GetYaxis()->CenterTitle();
	h_ptPhoton[hi]->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");    h_ptPhoton[hi]->GetXaxis()->CenterTitle();
	h_ptPhoton[hi]->Sumw2();

	sprintf(name, "h_ptJet_%s",cut1[hi].c_str());
	h_ptJet[hi]     =new TH1F(name,"pt of jet      ",100,20.0,2520.0);  // As we are putting a cut of 170.
	h_ptJet[hi]->GetYaxis()->SetTitle("Events/25 GeV");       h_ptJet[hi]->GetYaxis()->CenterTitle();
	h_ptJet[hi]->GetXaxis()->SetTitle("P_{T}^{jet} (GeV)");   h_ptJet[hi]->GetXaxis()->CenterTitle();
	h_ptJet[hi]->Sumw2();

	sprintf(name, "h_mass_VarBin_%s",cut1[hi].c_str());
	h_mass_VarBin[hi]      =new TH1F(name,"mass of photon+jet",nMassBins,MassBin);
	h_mass_VarBin[hi]->GetYaxis()->SetTitle("Events");                      h_mass_VarBin[hi]->GetYaxis()->CenterTitle();
	h_mass_VarBin[hi]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_VarBin[hi]->GetXaxis()->CenterTitle();
	h_mass_VarBin[hi]->Sumw2(); 

	sprintf(name, "h_mass_bin25_%s",cut1[hi].c_str());
	h_mass_bin25[hi]      =new TH1F(name,"mass of photon+jet",160,0.0,4000.0);
	h_mass_bin25[hi]->GetYaxis()->SetTitle("Events/25 GeV");               h_mass_bin25[hi]->GetYaxis()->CenterTitle();
	h_mass_bin25[hi]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_bin25[hi]->GetXaxis()->CenterTitle();
	h_mass_bin25[hi]->Sumw2(); 

	sprintf(name, "h_etaPhoton_%s",cut1[hi].c_str());
	h_etaPhoton[hi] =new TH1F(name,"eta of photon",100,-2.5,2.5);
	h_etaPhoton[hi]->GetYaxis()->SetTitle("Events");        h_etaPhoton[hi]->GetYaxis()->CenterTitle();
	h_etaPhoton[hi]->GetXaxis()->SetTitle("#eta^{#gamma}"); h_etaPhoton[hi]->GetXaxis()->CenterTitle();
	h_etaPhoton[hi]->Sumw2();

	sprintf(name, "h_etaJet_%s",cut1[hi].c_str());
	h_etaJet[hi]    =new TH1F(name,"eta of jet      ",200,-5.0,5.0);
	h_etaJet[hi]->GetYaxis()->SetTitle("Events");      h_etaJet[hi]->GetYaxis()->CenterTitle();
	h_etaJet[hi]->GetXaxis()->SetTitle("#eta^{jet}");  h_etaJet[hi]->GetXaxis()->CenterTitle();
	h_etaJet[hi]->Sumw2();

	sprintf(name, "h_ptPFMet_%s",cut1[hi].c_str());
	h_ptPFMet[hi]  = new TH1F(name,"pt of PFMet",100,0.0,500.0);
	h_ptPFMet[hi]->GetYaxis()->SetTitle("Events/5 GeV");           h_ptPFMet[hi]->GetYaxis()->CenterTitle();
	h_ptPFMet[hi]->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");    h_ptPFMet[hi]->GetXaxis()->CenterTitle();
	h_ptPFMet[hi]->Sumw2();

	sprintf(name, "h_SumEtPFMet_%s",cut1[hi].c_str());
	h_SumEtPFMet[hi]  = new TH1F(name,"SumET PF Met",80,0.0,4000.0);
	h_SumEtPFMet[hi]->GetYaxis()->SetTitle("Events/5 GeV");               h_SumEtPFMet[hi]->GetYaxis()->CenterTitle();
	h_SumEtPFMet[hi]->GetXaxis()->SetTitle("#sum#slash{E}_{T} (GeV)");    h_SumEtPFMet[hi]->GetXaxis()->CenterTitle();
	h_SumEtPFMet[hi]->Sumw2();
    }

    std::string cut2[3] = {"NoCut", "PhotonID", "final"};
    for( Int_t hi = 0 ; hi < 3 ; ++hi){

	sprintf(name, "h_mass_bin1_%s",cut2[hi].c_str());
	h_mass_bin1[hi]      =new TH1F(name,"mass of photon+jet",14000,0.0,14000.0);
	h_mass_bin1[hi]->GetYaxis()->SetTitle("Events/1 GeV");               h_mass_bin1[hi]->GetYaxis()->CenterTitle();
	h_mass_bin1[hi]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_bin1[hi]->GetXaxis()->CenterTitle();
	h_mass_bin1[hi]->Sumw2(); 

	sprintf(name, "h_mass_bin40_%s",cut2[hi].c_str());
	h_mass_bin40[hi]      =new TH1F(name,"mass of photon+jet",100,0.0,4000.0);
	h_mass_bin40[hi]->GetYaxis()->SetTitle("Events/40 GeV");               h_mass_bin40[hi]->GetYaxis()->CenterTitle();
	h_mass_bin40[hi]->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");   h_mass_bin40[hi]->GetXaxis()->CenterTitle();
	h_mass_bin40[hi]->Sumw2(); 

	sprintf(name, "h_Photon_SigmaIetaIeta_%s",cut2[hi].c_str());
	h_Photon_SigmaIetaIeta[hi]  = new TH1F(name,"Photon SigmaIetaIeta",100,0.0,0.05);
	h_Photon_SigmaIetaIeta[hi]->GetYaxis()->SetTitle("Events");                 h_Photon_SigmaIetaIeta[hi]->GetYaxis()->CenterTitle();
	h_Photon_SigmaIetaIeta[hi]->GetXaxis()->SetTitle("#sigma_{i#eta i#eta}");   h_Photon_SigmaIetaIeta[hi]->GetXaxis()->CenterTitle(); 
	h_Photon_SigmaIetaIeta[hi]->Sumw2();   

	sprintf(name, "h_PtPhotJet_%s",cut2[hi].c_str());
	h_PtPhotJet[hi]    =new TH2F(name,"Pt of Photon vs Jet ",80,0.0,2000.0,80,0.0,2000.0);
	h_PtPhotJet[hi]->GetYaxis()->SetTitle("P_{T}^{Jet}");     h_PtPhotJet[hi]->GetYaxis()->CenterTitle();
	h_PtPhotJet[hi]->GetXaxis()->SetTitle("P_{T}^{#gamma}");  h_PtPhotJet[hi]->GetXaxis()->CenterTitle();
	h_PtPhotJet[hi]->Sumw2();

	sprintf(name, "h_dEtadPhi_%s",cut2[hi].c_str());
	h_dEtadPhi[hi]    =new TH2F(name,"dEta vs dPhi ",120,0,6,64,0,3.2);
	h_dEtadPhi[hi]->GetYaxis()->SetTitle("#Delta#phi");  h_dEtadPhi[hi]->GetYaxis()->CenterTitle();
	h_dEtadPhi[hi]->GetXaxis()->SetTitle("#Delta#eta");  h_dEtadPhi[hi]->GetXaxis()->CenterTitle();
	h_dEtadPhi[hi]->Sumw2();

	sprintf(name, "h_DR_PhotonJet_%s",cut2[hi].c_str());
	h_DR_PhotonJet[hi] = new TH1F(name,"DeltaR between photon n jet",100,0.0,10.0);
	h_DR_PhotonJet[hi]->GetYaxis()->SetTitle("Events");        h_DR_PhotonJet[hi]->GetYaxis()->CenterTitle();
	h_DR_PhotonJet[hi]->GetXaxis()->SetTitle("#Delta R");      h_DR_PhotonJet[hi]->GetXaxis()->CenterTitle();
	h_DR_PhotonJet[hi]->Sumw2();

	sprintf(name, "h_YPhotJet_%s",cut2[hi].c_str());
	h_YPhotJet[hi] = new TH2F(name,"Y of Photon vs Jet", 100, -2.5, 2.5, 100, -2.5, 2.5);
	h_YPhotJet[hi]->GetYaxis()->SetTitle("y_{jet}");      h_YPhotJet[hi]->GetYaxis()->CenterTitle();
	h_YPhotJet[hi]->GetXaxis()->SetTitle("y_{#gamma}");   h_YPhotJet[hi]->GetXaxis()->CenterTitle();
	h_YPhotJet[hi]->Sumw2();

	sprintf(name, "h_dY_%s",cut2[hi].c_str());
	h_dY[hi]      =new TH1F(name,"dy of photon+jet",120,0,6);
	h_dY[hi]->GetYaxis()->SetTitle("Events");     h_dY[hi]->GetYaxis()->CenterTitle();
	h_dY[hi]->GetXaxis()->SetTitle("#Delta y");   h_dY[hi]->GetXaxis()->CenterTitle();
	h_dY[hi]->Sumw2();

	sprintf(name, "h_etaPhotJet_%s",cut2[hi].c_str());
	h_etaPhotJet[hi]    =new TH2F(name,"eta of Photon vs Jet ",100,-2.5,2.5,100,-2.5,2.5);
	h_etaPhotJet[hi]->GetYaxis()->SetTitle("#eta^{Jet}");     h_etaPhotJet[hi]->GetYaxis()->CenterTitle();
	h_etaPhotJet[hi]->GetXaxis()->SetTitle("#eta^{#gamma}");  h_etaPhotJet[hi]->GetXaxis()->CenterTitle();
	h_etaPhotJet[hi]->Sumw2();

	sprintf(name, "h_dEta_%s",cut2[hi].c_str());
	h_dEta[hi]      =new TH1F(name,"dEta of photon+jet",120,0,6);
	h_dEta[hi]->GetYaxis()->SetTitle("Events");        h_dEta[hi]->GetYaxis()->CenterTitle();
	h_dEta[hi]->GetXaxis()->SetTitle("#Delta #eta");   h_dEta[hi]->GetXaxis()->CenterTitle();
	h_dEta[hi]->Sumw2();

	sprintf(name, "h_phiPhotJet_%s",cut2[hi].c_str());
	h_phiPhotJet[hi]    =new TH2F(name,"phi of Photon vs Jet ",128,0.0,6.4,128,0.0,6.4);
	h_phiPhotJet[hi]->GetYaxis()->SetTitle("#phi^{Jet}");     h_phiPhotJet[hi]->GetYaxis()->CenterTitle();
	h_phiPhotJet[hi]->GetXaxis()->SetTitle("#phi^{#gamma}");  h_phiPhotJet[hi]->GetXaxis()->CenterTitle();
	h_phiPhotJet[hi]->Sumw2();

	sprintf(name, "h_dphi_%s",cut2[hi].c_str());
	h_dphi[hi]      =new TH1F(name,"dphi of photon+jet",64,0,3.2);
	h_dphi[hi]->GetYaxis()->SetTitle("Events");        h_dphi[hi]->GetYaxis()->CenterTitle();
	h_dphi[hi]->GetXaxis()->SetTitle("#Delta #phi");   h_dphi[hi]->GetXaxis()->CenterTitle();
	h_dphi[hi]->Sumw2();
    }

    std::string cut3[2] = {"NoCut", "final"};
    for( Int_t hi = 0 ; hi < 2 ; ++hi){

	sprintf(name, "h_Photon_SigmaEtaEta_%s",cut3[hi].c_str());
	h_Photon_SigmaEtaEta[hi]  = new TH1F(name,"Photon SigmaEtaEta ",100,0.0,0.05);
	h_Photon_SigmaEtaEta[hi]->GetYaxis()->SetTitle("Events");              h_Photon_SigmaEtaEta[hi]->GetYaxis()->CenterTitle();
	h_Photon_SigmaEtaEta[hi]->GetXaxis()->SetTitle("#sigma_{#eta#eta}");   h_Photon_SigmaEtaEta[hi]->GetXaxis()->CenterTitle();
	h_Photon_SigmaEtaEta[hi]->Sumw2();

	sprintf(name, "h_Photon_SigmaPhiPhi_%s",cut3[hi].c_str());
	h_Photon_SigmaPhiPhi[hi]  = new TH1F(name,"Photon SigmaPhiPhi ",100,0.0,0.05);
	h_Photon_SigmaPhiPhi[hi]->GetYaxis()->SetTitle("Events");              h_Photon_SigmaPhiPhi[hi]->GetYaxis()->CenterTitle();
	h_Photon_SigmaPhiPhi[hi]->GetXaxis()->SetTitle("#sigma_{#phi#phi}");   h_Photon_SigmaPhiPhi[hi]->GetXaxis()->CenterTitle();
	h_Photon_SigmaPhiPhi[hi]->Sumw2();

	sprintf(name, "h_Photon_SigmaIphiIphi_%s",cut3[hi].c_str());
	h_Photon_SigmaIphiIphi[hi]  = new TH1F(name,"Photon SigmaIphiIphi ",100,0.0,0.05);
	h_Photon_SigmaIphiIphi[hi]->GetYaxis()->SetTitle("Events");                 h_Photon_SigmaIphiIphi[hi]->GetYaxis()->CenterTitle();
	h_Photon_SigmaIphiIphi[hi]->GetXaxis()->SetTitle("#sigma_{i#phi i#phi}");   h_Photon_SigmaIphiIphi[hi]->GetXaxis()->CenterTitle();
	h_Photon_SigmaIphiIphi[hi]->Sumw2();

	sprintf(name, "h_Photon_swissCross_%s",cut3[hi].c_str());
	h_Photon_swissCross[hi]  = new TH1F(name,"Photon Swiss Cross ",40,-2.0,2.0);
	h_Photon_swissCross[hi]->GetYaxis()->SetTitle("Events");              h_Photon_swissCross[hi]->GetYaxis()->CenterTitle();
	h_Photon_swissCross[hi]->GetXaxis()->SetTitle("swiss cross");         h_Photon_swissCross[hi]->GetXaxis()->CenterTitle();
	h_Photon_swissCross[hi]->Sumw2();

	sprintf(name, "h_Photon_r9_%s",cut3[hi].c_str());
	h_Photon_r9[hi]  = new TH1F(name,"Photon R9 ",55,0.0,1.1);
	h_Photon_r9[hi]->GetYaxis()->SetTitle("Events");              h_Photon_r9[hi]->GetYaxis()->CenterTitle();
	h_Photon_r9[hi]->GetXaxis()->SetTitle("photon r9");           h_Photon_r9[hi]->GetXaxis()->CenterTitle();
	h_Photon_r9[hi]->Sumw2();

	sprintf(name, "h_PFiso_ChargedBarrel_%s",cut3[hi].c_str());
	h_PFiso_ChargedBarrel[hi]  = new TH1F(name,"PFIso Charged barrel ",750,0,15);
	h_PFiso_ChargedBarrel[hi]->GetYaxis()->SetTitle("Events");              h_PFiso_ChargedBarrel[hi]->GetYaxis()->CenterTitle();
	h_PFiso_ChargedBarrel[hi]->GetXaxis()->SetTitle("PFIso Charged");       h_PFiso_ChargedBarrel[hi]->GetXaxis()->CenterTitle();
	h_PFiso_ChargedBarrel[hi]->Sumw2();

	sprintf(name, "h_PFiso_PhotonBarrel_%s",cut3[hi].c_str());
	h_PFiso_PhotonBarrel[hi]  = new TH1F(name,"PFIso Photon barrel",75,0,15);
	h_PFiso_PhotonBarrel[hi]->GetYaxis()->SetTitle("Events");              h_PFiso_PhotonBarrel[hi]->GetYaxis()->CenterTitle();
	h_PFiso_PhotonBarrel[hi]->GetXaxis()->SetTitle("PFIso Photon");        h_PFiso_PhotonBarrel[hi]->GetXaxis()->CenterTitle();
	h_PFiso_PhotonBarrel[hi]->Sumw2();

	sprintf(name, "h_PFiso_NeutralBarrel_%s",cut3[hi].c_str());
	h_PFiso_NeutralBarrel[hi]  = new TH1F(name,"PFIso Neutral barrel",75,0,15);
	h_PFiso_NeutralBarrel[hi]->GetYaxis()->SetTitle("Events");              h_PFiso_NeutralBarrel[hi]->GetYaxis()->CenterTitle();
	h_PFiso_NeutralBarrel[hi]->GetXaxis()->SetTitle("PFIso Neutral");       h_PFiso_NeutralBarrel[hi]->GetXaxis()->CenterTitle();
	h_PFiso_NeutralBarrel[hi]->Sumw2();

	sprintf(name, "h_PFiso_SumBarrel_%s",cut3[hi].c_str());
	h_PFiso_SumBarrel[hi]  = new TH1F(name,"PFIso Sum barrel ",75,0,15);
	h_PFiso_SumBarrel[hi]->GetYaxis()->SetTitle("Events");              h_PFiso_SumBarrel[hi]->GetYaxis()->CenterTitle();
	h_PFiso_SumBarrel[hi]->GetXaxis()->SetTitle("PFIso Sum");           h_PFiso_SumBarrel[hi]->GetXaxis()->CenterTitle();
	h_PFiso_SumBarrel[hi]->Sumw2();

	sprintf(name, "h_PFiso_ElectronvetoBarrel_%s",cut3[hi].c_str());
	h_PFiso_ElectronvetoBarrel[hi]  = new TH1F(name,"PFIso Electronveto barrel ",3,0,3);
	h_PFiso_ElectronvetoBarrel[hi]->GetYaxis()->SetTitle("Events");              h_PFiso_ElectronvetoBarrel[hi]->GetYaxis()->CenterTitle();
	h_PFiso_ElectronvetoBarrel[hi]->GetXaxis()->SetTitle("PFIso Electron veto"); h_PFiso_ElectronvetoBarrel[hi]->GetXaxis()->CenterTitle();
	h_PFiso_ElectronvetoBarrel[hi]->Sumw2();

	sprintf(name, "h_HoEnewBarrel_%s",cut3[hi].c_str());
	h_HoEnewBarrel[hi]  = new TH1F(name,"PFIso HoE Barrel ",50,0,0.1);
	h_HoEnewBarrel[hi]->GetYaxis()->SetTitle("Events");              h_HoEnewBarrel[hi]->GetYaxis()->CenterTitle();
	h_HoEnewBarrel[hi]->GetXaxis()->SetTitle("PFIso HoE");           h_HoEnewBarrel[hi]->GetXaxis()->CenterTitle();
	h_HoEnewBarrel[hi]->Sumw2();

	sprintf(name, "h_Photon_ecalRecHitSumEtConeDR04_%s",cut3[hi].c_str());
	h_Photon_ecalRecHitSumEtConeDR04[hi]  = new TH1F(name,"ECAL Iso ",60,0,30);
	h_Photon_ecalRecHitSumEtConeDR04[hi]->GetYaxis()->SetTitle("Events");                                 h_Photon_ecalRecHitSumEtConeDR04[hi]->GetYaxis()->CenterTitle();
	h_Photon_ecalRecHitSumEtConeDR04[hi]->GetXaxis()->SetTitle("Photon_ecalRecHitSumEtConeDR04");         h_Photon_ecalRecHitSumEtConeDR04[hi]->GetXaxis()->CenterTitle();
	h_Photon_ecalRecHitSumEtConeDR04[hi]->Sumw2();

	sprintf(name, "h_Photon_hcalTowerSumEtConeDR04_%s",cut3[hi].c_str());
	h_Photon_hcalTowerSumEtConeDR04[hi]  = new TH1F(name,"HCAL Iso ",60,0,30);
	h_Photon_hcalTowerSumEtConeDR04[hi]->GetYaxis()->SetTitle("Events");                                  h_Photon_hcalTowerSumEtConeDR04[hi]->GetYaxis()->CenterTitle();
	h_Photon_hcalTowerSumEtConeDR04[hi]->GetXaxis()->SetTitle("h_Photon_hcalTowerSumEtConeDR04");         h_Photon_hcalTowerSumEtConeDR04[hi]->GetXaxis()->CenterTitle();
	h_Photon_hcalTowerSumEtConeDR04[hi]->Sumw2();

	sprintf(name, "h_Photon_trkSumPtHollowConeDR04_%s",cut3[hi].c_str());
	h_Photon_trkSumPtHollowConeDR04[hi]  = new TH1F(name,"Track Iso ",60,0,30);
	h_Photon_trkSumPtHollowConeDR04[hi]->GetYaxis()->SetTitle("Events");                                  h_Photon_trkSumPtHollowConeDR04[hi]->GetYaxis()->CenterTitle();
	h_Photon_trkSumPtHollowConeDR04[hi]->GetXaxis()->SetTitle("h_Photon_trkSumPtHollowConeDR04");         h_Photon_trkSumPtHollowConeDR04[hi]->GetXaxis()->CenterTitle();
	h_Photon_trkSumPtHollowConeDR04[hi]->Sumw2();

	sprintf(name, "h_HoE_%s",cut3[hi].c_str());
	h_HoE[hi]  = new TH1F(name,"HoE",50,0,0.1);
	h_HoE[hi]->GetYaxis()->SetTitle("Events");              h_HoE[hi]->GetYaxis()->CenterTitle();
	h_HoE[hi]->GetXaxis()->SetTitle("HoE");                 h_HoE[hi]->GetXaxis()->CenterTitle();
	h_HoE[hi]->Sumw2();

   	sprintf(name, "h_pfjet_NEF_%s",cut3[hi].c_str());
	h_pfjet_NEF[hi] = new TH1F(name,"Neutral EM Fraction",25,0,1);
	h_pfjet_NEF[hi]->GetYaxis()->SetTitle("");                h_pfjet_NEF[hi]->GetYaxis()->CenterTitle();
	h_pfjet_NEF[hi]->GetXaxis()->SetTitle("pfjet_NEF");       h_pfjet_NEF[hi]->GetXaxis()->CenterTitle();
	h_pfjet_NEF[hi]->Sumw2();

	sprintf(name, "h_pfjet_NHF_%s",cut3[hi].c_str());
	h_pfjet_NHF[hi] = new TH1F(name,"Neutral hadron Fraction",25,0,1);
	h_pfjet_NHF[hi]->GetYaxis()->SetTitle("");                h_pfjet_NHF[hi]->GetYaxis()->CenterTitle();
	h_pfjet_NHF[hi]->GetXaxis()->SetTitle("pfjet_NHF");       h_pfjet_NHF[hi]->GetXaxis()->CenterTitle();
	h_pfjet_NHF[hi]->Sumw2();
	
	sprintf(name, "h_pfjet_CEF_%s",cut3[hi].c_str());
	h_pfjet_CEF[hi] = new TH1F(name,"Charged EM Fraction",25,0,1);
	h_pfjet_CEF[hi]->GetYaxis()->SetTitle("");                h_pfjet_CEF[hi]->GetYaxis()->CenterTitle();
	h_pfjet_CEF[hi]->GetXaxis()->SetTitle("pfjet_CEF");       h_pfjet_CEF[hi]->GetXaxis()->CenterTitle();
	h_pfjet_CEF[hi]->Sumw2();
	
	sprintf(name, "h_pfjet_CHF_%s",cut3[hi].c_str());
	h_pfjet_CHF[hi] = new TH1F(name,"Charged hadron Fraction",25,0,1);
	h_pfjet_CHF[hi]->GetYaxis()->SetTitle("");                h_pfjet_CHF[hi]->GetYaxis()->CenterTitle();
	h_pfjet_CHF[hi]->GetXaxis()->SetTitle("pfjet_CHF");       h_pfjet_CHF[hi]->GetXaxis()->CenterTitle();
	h_pfjet_CHF[hi]->Sumw2();
	
	sprintf(name, "h_pfjet_NConstituents_%s",cut3[hi].c_str());
	h_pfjet_NConstituents[hi] = new TH1F(name,"NConstituents",50,0,50);
	h_pfjet_NConstituents[hi]->GetYaxis()->SetTitle("");                     h_pfjet_NConstituents[hi]->GetYaxis()->CenterTitle();
	h_pfjet_NConstituents[hi]->GetXaxis()->SetTitle("pfjet_NConstituents");  h_pfjet_NConstituents[hi]->GetXaxis()->CenterTitle();
	h_pfjet_NConstituents[hi]->Sumw2();
	
	sprintf(name, "h_pfjet_ChargeMultiplicity_%s",cut3[hi].c_str());
	h_pfjet_ChargeMultiplicity[hi] = new TH1F(name,"ChargeMultiplicity",50,0,50);
	h_pfjet_ChargeMultiplicity[hi]->GetYaxis()->SetTitle("");                             h_pfjet_ChargeMultiplicity[hi]->GetYaxis()->CenterTitle();
	h_pfjet_ChargeMultiplicity[hi]->GetXaxis()->SetTitle("pfjet_ChargeMultiplicity");     h_pfjet_ChargeMultiplicity[hi]->GetXaxis()->CenterTitle();
	h_pfjet_ChargeMultiplicity[hi]->Sumw2(); 
    }

    h_CorrPFiso_Charged = new TH1F("h_CorrPFiso_Charged","Rho Corrected PFIso Charged  ",750,0,15);
    h_CorrPFiso_Charged->GetYaxis()->SetTitle("Events");                        h_CorrPFiso_Charged->GetYaxis()->CenterTitle();
    h_CorrPFiso_Charged->GetXaxis()->SetTitle("Corrected PFIso Charged");       h_CorrPFiso_Charged->GetXaxis()->CenterTitle();
    h_CorrPFiso_Charged->Sumw2();

    h_CorrPFiso_Photon = new TH1F("h_CorrPFiso_Photon","Rho Corrected PFIso Photon barrel",75,0,15);
    h_CorrPFiso_Photon->GetYaxis()->SetTitle("Events");                         h_CorrPFiso_Photon->GetYaxis()->CenterTitle();
    h_CorrPFiso_Photon->GetXaxis()->SetTitle("Corrected PFIso Photon");         h_CorrPFiso_Photon->GetXaxis()->CenterTitle();
    h_CorrPFiso_Photon->Sumw2();

    h_CorrPFiso_Neutral = new TH1F("h_CorrPFiso_Neutral","Rho Corrected PFIso Neutral barrel",75,0,15);
    h_CorrPFiso_Neutral->GetYaxis()->SetTitle("Events");                        h_CorrPFiso_Neutral->GetYaxis()->CenterTitle();
    h_CorrPFiso_Neutral->GetXaxis()->SetTitle("Corrected PFIso Neutral");       h_CorrPFiso_Neutral->GetXaxis()->CenterTitle();
    h_CorrPFiso_Neutral->Sumw2();

    std::string cut4[3] = {"NoPU", "PU", "final"};
    for( Int_t hi = 0 ; hi < 3 ; ++hi){
	sprintf(name, "h_Vertices_%s",cut4[hi].c_str());
	h_Vertices[hi] = new TH1F(name,"Vertices", 60, 0, 60);
    }

    // Histogram fornumber for photons and jets
    std::string cut6[2] ={ "noCut", "final"};
    for(Int_t hi = 0 ; hi < 2 ; hi++){
	sprintf(name, "h_nPhoton_%s",cut6[hi].c_str());
	h_nPhoton[hi] = new TH1F(name,"no of Photons",10,0,10);

	sprintf(name, "h_nIsoPhoton_%s",cut6[hi].c_str());
	h_nIsoPhoton[hi] = new TH1F(name,"no of Iso Photons",10,0,10);
	
	sprintf(name, "h_nIsoPhotonPt_%s",cut6[hi].c_str());
	h_nIsoPhotonPt[hi] = new TH2F(name,"no of Iso Photons vs Pt",80,0.0,2000.0,10,0.0,10.0);

	sprintf(name, "h_nJet_%s",cut6[hi].c_str());
	h_nJet[hi] = new TH1F(name,"no of Jets",20,0,20);

	sprintf(name, "h_nIsoJet_%s",cut6[hi].c_str());
	h_nIsoJet[hi] = new TH1F(name,"no of selected Jets",20,0,20);
	
	sprintf(name, "h_nIsoJetPt_%s",cut6[hi].c_str());
	h_nIsoJetPt[hi] = new TH2F(name,"no of selected Jets vs Pt",80,0.0,2000.0,20,0.0,20.0);
    }

    std::string cut7[2] = {"NoCut", "final"};
    for( Int_t hi = 0 ; hi < 2 ; ++hi){

	sprintf(name, "h_ptSecondPhoton_%s",cut7[hi].c_str());
	h_ptSecondPhoton[hi]  = new TH1F(name,"pt of Second photon",80,0.0,2000.0);
	h_ptSecondPhoton[hi]->GetYaxis()->SetTitle("Events/25 GeV");           h_ptSecondPhoton[hi]->GetYaxis()->CenterTitle();
	h_ptSecondPhoton[hi]->GetXaxis()->SetTitle("P_{T,2}^{#gamma} (GeV)");  h_ptSecondPhoton[hi]->GetXaxis()->CenterTitle();
	h_ptSecondPhoton[hi]->Sumw2();

	sprintf(name, "h_etaSecondPhoton_%s",cut7[hi].c_str());
	h_etaSecondPhoton[hi] =new TH1F(name,"eta of Secondphoton",200,-5.0,5.0);
	h_etaSecondPhoton[hi]->GetYaxis()->SetTitle("Events");            h_etaSecondPhoton[hi]->GetYaxis()->CenterTitle();
	h_etaSecondPhoton[hi]->GetXaxis()->SetTitle("#eta_{2}^{#gamma}"); h_etaSecondPhoton[hi]->GetXaxis()->CenterTitle();
	h_etaSecondPhoton[hi]->Sumw2();
}

    h_PC = new TH1F ("h_PC","Photon Candidate", 10, 0, 10);
    h_JC = new TH1F ("h_JC","Photon Candidate", 20, 0, 20);

    h_trueinteractions_MC = new TH1F ("h_trueinteractions_MC","trueinteractions", 60, 0, 60);
    h_trueinteractions_MC_PU = new TH1F ("h_trueinteractions_MC_PU","trueinteractions_MC_PU", 60, 0, 60);

    /// --- For the centrality ratio
    h_InnerRadius_InvMass  =new TH1F("h_InnerRadius_InvMass","mass of photon+jet",160,0.0,4000.0);
    h_OuterRadius_InvMass  =new TH1F("h_OuterRadius_InvMass","mass of photon+jet",160,0.0,4000.0);

} // end of bookhistos loop



//----------------
// Scraping variable
Bool_t PostAnalyzerMC::NonScraping(){
    Bool_t noScraping = !Scraping_isScrapingEvent ;
    return noScraping;
}

//--------------------
// Spike Cuts
Bool_t PostAnalyzerMC::NoSpike(Int_t ipho){
    Bool_t passSpike = false;

    if(  fabs(getLICTD(ipho))              < 5.0    && 
	 fabs(Photon_timing_xtal[ipho][0]) < 3.0    && 
	 Photon_SigmaIetaIeta[ipho]        > 0.001  && 
	 Photon_SigmaIphiIphi[ipho]        > 0.001  && 
         Photonr9[ipho]                    < 1.0){
	passSpike = true ;
    }

    return passSpike;

}

//----------------
//Primary Vertex
Bool_t PostAnalyzerMC::PrimaryVertex(Int_t &goodVertex){
    Bool_t VertexAccepted = false;
    goodVertex=0;

    for(Int_t i=0; i<Vertex_n && i<200;++i){
	if(  fabs(Vertex_z[i]) <= Cvertex_z    && 
	     Vertex_ndof[i]    >= Cvertex_ndof && 
 	     !Vertex_isFake[i]                 &&
	     fabs(Vertex_d0[i])<= Cvertex_rho)

	    goodVertex++ ;
    }
    if(goodVertex > 0)VertexAccepted=true;

    return VertexAccepted;
}

//-----------------
// Compute Rapidity
Double_t PostAnalyzerMC::getRapidity(Double_t r_E, Double_t r_Pz){
    Double_t rapdty  =  0.5*log( (r_E + r_Pz)/(r_E - r_Pz) );

    return rapdty ;
}

//-----------------
// Compute DeltaR
Double_t PostAnalyzerMC::getDR(Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){

    Double_t dPhi  = fabs(phi1 - phi2);
    Double_t twopi = 2.0*(TMath::Pi());
    Double_t dEta = fabs(eta1 - eta2);
    Double_t DR=0.0;

    if(dPhi < 0) dPhi = - dPhi; 
    if(dPhi >= (twopi - dPhi))dPhi = twopi - dPhi;

    DR = pow((dPhi*dPhi + dEta*dEta),0.5);
    return DR;
}

//-----------------
//Compute Delta Eta
Double_t PostAnalyzerMC::getDEta(Double_t eta1, Double_t eta2){

    Double_t dEta = fabs(eta1 - eta2);
    return dEta;
}

//------------------
//Compute Delta Phi
Double_t PostAnalyzerMC::getDPhi(Double_t phi1, Double_t phi2){

    Double_t dPhi  = fabs(phi1 - phi2);
    Double_t twopi = 2.0*(TMath::Pi());

    if(dPhi < 0) dPhi = - dPhi;
    if(dPhi >= (twopi - dPhi))dPhi = twopi - dPhi;

    return dPhi;
}

//----------------
//Compute Invariant mass for Gamma and Jet
Double_t PostAnalyzerMC::getMass(Int_t pho_i, Int_t jet_i){

    Double_t mass=0.0;

    Double_t E  = Photon_E[pho_i] + pfJet_E[jet_i];
    Double_t PX = Photon_px[pho_i] + pfJet_px[jet_i];
    Double_t PY = Photon_py[pho_i] + pfJet_py[jet_i];
    Double_t PZ = Photon_pz[pho_i] + pfJet_pz[jet_i];

    mass = pow((E*E - PX*PX - PY*PY - PZ*PZ),0.5);

    return mass;
}

//-----------------
// Compute LICTD LargeIntraClusterTimeDifference
Double_t PostAnalyzerMC::getLICTD(Int_t i){

    Double_t SeedTime = -999;
    Double_t SeedE    = -999;

    Int_t CrysIdx     = -1;

    for(Int_t k=0; k<Photon_ncrys[i]&&k<100;++k){
	Float_t CrysE = Photon_energy_xtal[i][k];
	if(CrysE > SeedE){
	    SeedE    = CrysE ;
	    SeedTime = Photon_timing_xtal[i][k];
	    CrysIdx  = k;
	}
    }

    Float_t LICTD = 99.0;

    if(fabs(SeedTime) < 3.0){
	LICTD = 0.0;
	Int_t CrysCrys   = -1;
	Int_t CrysThresh = 0;

	for(Int_t k=0;k<Photon_ncrys[i]&&k<100;++k){
	    if(CrysIdx == k)continue;
	    Float_t CrysE = Photon_energy_xtal[i][k];

	    if(CrysE > 1.0){
		CrysThresh++;
		Float_t timeDiff = Photon_timing_xtal[i][CrysIdx] - Photon_timing_xtal[i][k];
		if(fabs(timeDiff) > fabs(LICTD)){
		    LICTD    = timeDiff;
		    CrysCrys = k;
		}
	    }
	}
    }

    return LICTD;
}


//Tight Photon PF ISo for 2012 to select photon candidate
Bool_t PostAnalyzerMC::TightPhotonPFIso(Int_t ipho){

    Bool_t tightID = false;

    if(fabs(Photon_sc_eta[ipho]) <= 1.4442 ) {  // For Barrel
	tightID = (Photon_HoEnew[ipho]   < 0.05)          &&
	    (Photon_SigmaIetaIeta[ipho]  < 0.011)         &&
	    (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 0.7)                            &&
	    (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 0.4+0.04*Photon_pt[ipho])       &&
	    (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 0.5+0.005*Photon_pt[ipho])      &&
	    (Photon_Electronveto[ipho]  == 1);
    }

    if(fabs(Photon_sc_eta[ipho]) > 1.4442 )  {  // For EndCap
	tightID = (Photon_HoEnew[ipho]   < 0.05)          &&
	    (Photon_SigmaIetaIeta[ipho]  < 0.031)         &&
	    (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 0.5)                            &&
	    (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 1.5+0.04*Photon_pt[ipho])       &&
	    (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.0+0.005*Photon_pt[ipho])      &&
	    (Photon_Electronveto[ipho]  == 1);

    }
    return tightID;
}

Bool_t PostAnalyzerMC::MediumPhotonPFIso(Int_t ipho){

    Bool_t mediumID = false;

    if(fabs(Photon_sc_eta[ipho]) <= 1.4442 ) {  // For Barrel
	mediumID = (Photon_HoEnew[ipho]   < 0.05)          &&
	    (Photon_SigmaIetaIeta[ipho]  < 0.011)         &&
	    (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 1.5)                            &&
	    (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 1.0+0.04*Photon_pt[ipho])       &&
	    (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 0.7+0.005*Photon_pt[ipho])      &&
	    (Photon_Electronveto[ipho]  == 1);
    }

    if(fabs(Photon_sc_eta[ipho]) > 1.4442 )  {  // For EndCap
	mediumID = (Photon_HoEnew[ipho]   < 0.05)          &&
	    (Photon_SigmaIetaIeta[ipho]  < 0.033)         &&
	    (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 1.2)                            &&
	    (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 1.5+0.04*Photon_pt[ipho])       &&
	    (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.0+0.005*Photon_pt[ipho])      &&
	    (Photon_Electronveto[ipho]  == 1);

    }
    return mediumID;
}

Bool_t PostAnalyzerMC::LoosePhotonPFIso(Int_t ipho){

    Bool_t looseID = false;

    if(fabs(Photon_sc_eta[ipho]) <= 1.4442 ) {  // For Barrel
	looseID = (Photon_HoEnew[ipho]   < 0.05)          &&
	    (Photon_SigmaIetaIeta[ipho]  < 0.012)         &&
	    (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 2.6)                            &&
	    (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 3.5+0.04*Photon_pt[ipho])       &&
	    (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.3+0.005*Photon_pt[ipho])      &&
	    (Photon_Electronveto[ipho]  == 1);
    }

    if(fabs(Photon_sc_eta[ipho]) > 1.4442 )  {  // For EndCap
	looseID = (Photon_HoEnew[ipho]   < 0.05)          &&
	    (Photon_SigmaIetaIeta[ipho]  < 0.034)         &&
	    (TMath::Max(((PFiso_Charged03[ipho]) - rho*EAElectroncharged(Photon_sc_eta[ipho])),0.0) < 2.3)                            &&
	    (TMath::Max(((PFiso_Neutral03[ipho]) - rho*EAElectronneutral(Photon_sc_eta[ipho])),0.0) < 2.9+0.04*Photon_pt[ipho])       &&
//	    (TMath::Max(((PFiso_Photon03[ipho])  - rho*EAElectronphoton(Photon_sc_eta[ipho])),0.0)  < 1.0+0.005*Photon_pt[ipho])      &&
	    (Photon_Electronveto[ipho]  == 1);

    }
    return looseID;
}


// Tight Jet ID and jet seperation from photon to select an isolated jet.
Bool_t PostAnalyzerMC::TightJetID( Int_t ijet){

    Bool_t ID = false;
    if(fabs(pfJet_eta[ijet]) <= 2.4){
	ID = (pfjet_NEF[ijet]           < 0.90) &&  
	    (pfjet_NHF[ijet]            < 0.90) &&
	    (pfjet_NConstituents[ijet]  > 1)    &&
	    (pfjet_CEF[ijet]            < 0.99) &&
	    (pfjet_CHF[ijet]            > 0)    &&
	    (pfjet_NCH[ijet]            > 0);
    }

    if(fabs(pfJet_eta[ijet]) > 2.4){
	ID = (pfjet_NEF[ijet]            < 0.90) &&  
	    (pfjet_NHF[ijet]            < 0.90) &&
	    (pfjet_NConstituents[ijet]  > 1);
    }

    return ID;
}

// Effective area to be needed in PF Iso for photon ID
Double_t PostAnalyzerMC::EAElectroncharged(Double_t eta){
    Float_t EffectiveArea=0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.012;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.010;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.014;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.012;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.016;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.020;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.012;

    return EffectiveArea;
}

Double_t PostAnalyzerMC::EAElectronneutral(Double_t eta){
    Float_t EffectiveArea=0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.030;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.057;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.039;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.015;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.024;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.039;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.072;

    return EffectiveArea;
}


Double_t PostAnalyzerMC::EAElectronphoton(Double_t eta){
    Float_t EffectiveArea=0.;
    if(fabs(eta) >= 0.0   && fabs(eta) < 1.0   ) EffectiveArea = 0.148;
    if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.130;
    if(fabs(eta) >= 1.479 && fabs(eta) < 2.0   ) EffectiveArea = 0.112;
    if(fabs(eta) >= 2.0   && fabs(eta) < 2.2   ) EffectiveArea = 0.216;
    if(fabs(eta) >= 2.2   && fabs(eta) < 2.3   ) EffectiveArea = 0.262;
    if(fabs(eta) >= 2.3   && fabs(eta) < 2.4   ) EffectiveArea = 0.260;
    if(fabs(eta) >= 2.4                        ) EffectiveArea = 0.266;

    return EffectiveArea;
}


Double_t PostAnalyzerMC::EXOtightPhoID(Int_t ipho){

    Bool_t ID = false;

    if(fabs(Photon_sc_eta[ipho]) < 1.4442 ){
	ID = (  Photon_HoE[ipho] < 0.05                                                         &&
		Photon_ecalRecHitSumEtConeDR04[ipho] < 4.2 + 0.006*Photon_pt[ipho]  + 0.1830*rho25 &&
		Photon_hcalTowerSumEtConeDR04[ipho]  < 2.2 + 0.0025*Photon_pt[ipho] + 0.0620*rho25 &&
		Photon_trkSumPtHollowConeDR04[ipho]  < 2.0 + 0.001*Photon_pt[ipho]  + 0.0167*rho25 &&
		Photon_SigmaIetaIeta[ipho]           < 0.011                                    &&
		(!Photon_hasPixelSeed[ipho])
	     );
    }
    if(fabs(Photon_sc_eta[ipho]) > 1.4442 ){
	ID = (  Photon_HoE[ipho] < 0.05                                                      &&
		Photon_ecalRecHitSumEtConeDR04[ipho] < 4.2+0.006*Photon_pt[ipho]  + 0.090*rho25 &&
		Photon_hcalTowerSumEtConeDR04[ipho]  < 2.2+0.0025*Photon_pt[ipho] + 0.180*rho25 &&
		Photon_trkSumPtHollowConeDR04[ipho]  < 2.0+0.001*Photon_pt[ipho]  + 0.032*rho25 &&
		Photon_SigmaIetaIeta[ipho]           < 0.030                                 &&
		(!Photon_hasPixelSeed[ipho])
	     );
    }

    return ID;
}


void PostAnalyzerMC::LumiReWeighting(){

//    std::cout<<" Taking into account pileUp_${i}.root for MC & DataTrue2012PU60_hist.root for Data"<<std::endl;

//    TFile *fData = TFile::Open("/uscms_data/d3/varun/2012/QstarGJ_8TeV/PileUpHistograms/Data_19pb/pileUpTrue60/DataTrue2012PU60_hist.root");
    TFile *fData = TFile::Open("/uscms_data/d3/varun/2012/QstarGJ_8TeV/PileUpHistograms/Data13July/13JulPileUpTrue60/DataTrue60PU_13Jul_2012ABCD.root");
    TH1F* dataPU = (TH1F*)fData->Get("pileup");

    TFile *fMC = TFile::Open("/uscms_data/d3/varun/2012/QstarGJ_8TeV/PileUpHistograms/Summer12_MC53X/pileUp_QstarToGJ_M_700.root");
//    TFile *fMC = TFile::Open("/uscms_data/d3/varun/2012/QstarGJ_8TeV/PileUpHistograms/Summer12_MC53X/pileUp_${i}.root");
    TH1F* mcPU = (TH1F*)fMC->Get("h_mcPileUp60"); 

    std::vector<float> MCpileup;
    std::vector<float> datapileup;

    for(int i=0; i<60; i++){
	MCpileup.push_back(mcPU->GetBinContent(i+1));
	datapileup.push_back(dataPU->GetBinContent(i+1));
    } 

    if(MCpileup.size() != datapileup.size() ){
	std::cerr <<"ERROR: LumiReWeighting: input vectors have different sizes. Quitting... \n";
	return;
    }

    Int_t NBins = MCpileup.size();
    MC_distr_ = new TH1F("MC_distr","MC dist",NBins,0, 60);
    Data_distr_ = new TH1F("Data_distr","Data dist",NBins,0, 60);
    den = new TH1F("den","denominator",NBins,0, 60);
    weights = new TH1F("weights","weights",NBins,0, 60);

    for(int ibin = 1; ibin<NBins+1; ++ibin ){
	MC_distr_->SetBinContent(ibin, MCpileup[ibin-1]);
	den->SetBinContent(ibin, MCpileup[ibin-1]);
	Data_distr_->SetBinContent(ibin, datapileup[ibin-1]);
	weights->SetBinContent(ibin, datapileup[ibin-1]);
    }

    weights->Scale(1.0/weights->Integral());
    den->Scale(1.0/den->Integral());

    weights->Divide(den);

}

Double_t PostAnalyzerMC::puweight(Float_t npv){
    Int_t bin = weights->GetXaxis()->FindBin( npv );
    return weights->GetBinContent( bin );
}

#endif // #ifdef PostAnalyzerMC_cxx
EOF


cat>analysis_${FileNameTag}_${p}.C<<EOF
#include "PostAnalyzerMC.C"
#include "TROOT.h"

int main(){
    //PostAnalyzerMC *a = new PostAnalyzerMC();
    PostAnalyzerMC a ;
    a.Loop();
    return 0;
}
EOF


###Now compilation

g++ -Wno-deprecated analysis_${FileNameTag}_${p}.C -o ${FileNameTag}_${p}.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`


echo "--------------------Submitting Job for ${FileNameTag}_${p} Files"   
echo "------------------------Submitting Job #-> ${p}  for  ${sf} to ${ef} Files  -----Total = ${Tot}"

##change for next file
@ sf = ${sf} + ${r}
@ ef = ${ef} + ${r}
###Submit jobs    

./MakeCondorFiles_new.csh ${FileNameTag}_${p} ${FileNameTag}_dataset.txt


rm PostAnalyzerMC.C
rm PostAnalyzerMC.h 

end  ### for while loop


end  ### for foreach loop

