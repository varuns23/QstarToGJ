MyPlottingScript(){

    const Int_t nFiles = 5;
    const Int_t nSigFiles = 9;

    TString FileName[nFiles+1] = {
	"./rootFiles/DiJet_53X8T19pb_Pt170Mass560dEta20TID.root",        // DiJet Bkg
	"./rootFiles/GJ_53X8T19pb_Pt170Mass560dEta20TID.root",           // PhotonJet Bkg Pythia 
	"./rootFiles/EWK_53X8T19pb_Pt170Mass560dEta20TID.root",          // EWK Bkg
//	"./rootFiles/Data_53X8T19pb_Pt170Mass560dEta20TID_HLT150.root",
	"./rootFiles/Data_53X8T19pb_Pt170Mass560dEta20TID.root",         // Data
	"./rootFiles/QstarToGJ_M_1000.root",       // Qstar signal  1 TeV 
	"./rootFiles/QstarToGJ_M_3000.root"        // Qstar signal  3 TeV
    };

    TString FileSignals[nSigFiles] = { 
	"./rootFiles/QstarToGJ_M_700.root",
	"./rootFiles/QstarToGJ_M_1000.root",
	"./rootFiles/QstarToGJ_M_1200.root",
	"./rootFiles/QstarToGJ_M_1500.root",
	"./rootFiles/QstarToGJ_M_1700.root",
	"./rootFiles/QstarToGJ_M_2000.root",
	"./rootFiles/QstarToGJ_M_2500.root",
	"./rootFiles/QstarToGJ_M_3000.root",
	"./rootFiles/QstarToGJ_M_3500.root"
    };

    const Float_t scaleMC   = 0.98*0.964;
    const Float_t GJkfactor = 1.4;
    //    const Float_t GJkfactor = 1.44;
    TString Plotslocation  = "./Plots_TEMP/" ;
    Bool_t  SavePlots      = true ;

    //--- Give here the Data-Signal-Bkg-file names ++++++++++++++++++++++++
    TFile* fileBkg1 = new TFile(FileName[0],"READ");
    TFile* fileBkg2 = new TFile(FileName[1],"READ");
    TFile* fileBkg3 = new TFile(FileName[2],"READ");
    TFile* fileData = new TFile(FileName[3],"READ");
    TFile* fileSig1 = new TFile(FileName[4],"READ"); 
    TFile* fileSig2 = new TFile(FileName[5],"READ"); 
    //------------------------------------------------------------


//    Bool_t b_Pt=true, b_Eta=true, b_Vertex=true, b_PhoPFIso=true, b_PhoID=true, b_JetID=true, b_ShowerShape=true, b_2DHist=true, b_PC=true, b_CutFlow=true, b_TurnOn=true, b_2DAllSignal=true, b_1DAllSignal=true, b_1DAllSignalImpose=true, b_Cummulative=true, b_ReverseCummulative=true;
      Bool_t b_Pt=false, b_Eta=false, b_Vertex=false, b_PhoPFIso=false, b_PhoID=false, b_JetID=false, b_ShowerShape=false, b_2DHist=false, b_PC=false, b_CutFlow=false, b_TurnOn=true, b_2DAllSignal=false, b_1DAllSignal=false, b_1DAllSignalImpose=false, b_Cummulative=false, b_ReverseCummulative=false;



    if(b_Pt){ ///----->  1
	const Int_t nvar1 = 12;
	TString var1[nvar1] = {"h_ptPhoton_finalPU", "h_ptJet_finalPU", "h_ptPhoton_massVarBin", "h_ptPhoton_VarBin", "h_ptJet_massVarBin", "h_ptJet_VarBin", "h_mass_VarBin_finalPU", "h_mass_bin25_finalPU", "h_mass_bin40_final", "h_ptPFMet_finalPU", "h_SumEtPFMet_finalPU", "h_ptSecondPhoton_final"};
	TString PlotsName1[nvar1] = {"ptPhoton_final", "ptJet_final", "ptPhoton_massVarBin", "ptPhoton_VarBin", "ptJet_massVarBin", "ptJet_VarBin", "mass_VarBin_final", "mass_bin25_final", "mass_bin40_final", "ptPFMet_final", "SumEtPFMet_final", "ptSecondPhoton_final"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar1 ; ++i){

	  Double_t x_min , x_max;
	  if(i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5 ) x_min = 0; x_max = 2500;
	  if(i == 6 || i == 7 || i == 8 ) x_min = 0; x_max = 4000;
	  if(i == 9 || i == 10 || i == 11 ) x_min = 0; x_max = 4000;

	    Bool_t d_logy=true;

	    Double_t leg_xmin=0.7, leg_ymin=0.7, leg_xmax=0.9, leg_ymax=0.9;
	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var1,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName1, SavePlots, d_logy);
	}
    }

    if(b_Eta){ ///----->  2
	// Number of variables to be plotted with their names +++++++++++++++++
	//	const Int_t nvar2 = 1;
	//	TString var2[nvar2] = { "h_etaJet_PhotonID"};
	//	TString PlotsName2[nvar2] = {"etaJet_PhotonID"};
	const Int_t nvar2 = 7;
	TString var2[nvar2] = { "h_etaPhoton_finalPU", "h_etaJet_finalPU", "h_DR_PhotonJet_final", "h_dY_final", "h_dEta_final", "h_dphi_final", "h_etaSecondPhoton_final" };
	TString PlotsName2[nvar2] = {"etaPhoton_final", "etaJet_final", "DR_PhotonJet_final", "dY_final", "dEta_final", "dphi_final", "etaSecondPhoton_final"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar2 ; ++i){

	    //	    if(i == 0 || i == 1){ Double_t x_min = -5.0, x_max = 5.0;}
	    if(i == 0 ){ Double_t x_min = -2.5, x_max = 2.5;}
	    else if(i == 1){ Double_t x_min = -3.0, x_max = 3.0;}
	    else if(i == 2){ Double_t x_min = 1.2, x_max = 4.5;}
	    else if(i == 3 || i == 4){ Double_t x_min = 0.0, x_max = 1.5;}
	    else if(i == 5){ Double_t x_min = 1.4, x_max = 3.2;}
	    else if(i == 6){ Double_t x_min = -3.0, x_max = 3.0;}

	    Double_t leg_xmin=0.4, leg_ymin=0.05, leg_xmax=0.6, leg_ymax=0.3;
	    Bool_t d_logy=true;
	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var2,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName2, SavePlots, d_logy);
	    //		    scaleMC, Plotslocation, PlotsName2, false, d_logy);

	}
    }

    if(b_Vertex){ ///----->  3
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar3 = 2;
	//	TString var3[nvar3] = { "h_Vertices_NoPU", "h_Vertices_PU", "h_Vertices_final"};
	TString var3[nvar3] = { "h_Vertices_PU", "h_Vertices_NoPU"};
	//	TString PlotsName3[nvar3] = {"Vertices_NoPU", "Vertices_PU", "Vertices_final"};
	TString PlotsName3[nvar3] = {"Total Bkg after re-weighting", "Total Bkg before re-weighting"};
	TH1F *h3Bkg1[nvar3], *h3Bkg2[nvar3], *h3Bkg3[nvar3], *h3Data[nvar3];
	TLegend *leg3 = new TLegend(0.7,0.75,0.9,0.9);
	leg3->SetFillColor(0); leg3->SetTextFont(40); leg3->SetBorderSize(1);
//	for(Int_t j = 0 ; j < 2 ; ++j){	
	    for(Int_t i = 0 ; i < nvar3 ; ++i){

		h3Bkg1[i] = (TH1F*)fileBkg1->Get(var3[i]);   //DiJet MC
		h3Bkg2[i] = (TH1F*)fileBkg2->Get(var3[i]);   //GJ MC
		h3Bkg3[i] = (TH1F*)fileBkg3->Get(var3[i]);   //EWK
		h3Data[i] = (TH1F*)fileData->Get(var3[i]);   //Data


		h3Bkg1[i]->SetStats(0); 
		h3Bkg1[i]->Add(h3Bkg2[i],GJkfactor);
		h3Bkg1[i]->Add(h3Bkg3[i]);
		h3Bkg1[i]->Scale(1.0/h3Bkg1[i]->Integral());
		h3Data[i]->Scale(1.0/h3Data[i]->Integral());
		h3Data[i]->SetMarkerStyle(8);
		h3Data[i]->SetMarkerSize(1);
	    }
	    h3Bkg1[0]->SetFillColor(kYellow);
	    h3Bkg1[1]->SetLineColor(kBlue);
	    h3Bkg1[1]->SetFillColor(kBlue);
	    h3Bkg1[1]->SetFillStyle(3003);
	    h3Bkg1[1]->SetLineWidth(2);
	    
	    h3Bkg1[1]->Draw("HIST");
	    h3Bkg1[0]->Draw("sameHIST");
	    h3Data[1]->Draw("sameP");
	    
	    leg3->AddEntry(h3Data[0],"Data","P");
	    leg3->AddEntry(h3Bkg1[1],PlotsName3[1],"F");
	    leg3->AddEntry(h3Bkg1[0],PlotsName3[0],"F");
	    leg3->Draw();
	    //	    if( j == 0 ){c1->SaveAs(Plotslocation+"Vertices"+".eps");}else {
	    //		c1->SetLogy();
	    //		c1->SaveAs(Plotslocation+"VerticesLog"+".eps");
//	    }
//	}

    }

    if(b_PhoPFIso){ ///----->  4
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar4 = 8;
	TString var4[nvar4] = { "h_CorrPFiso_Charged", "h_CorrPFiso_Neutral", "h_CorrPFiso_Photon", "h_PFiso_ChargedBarrel_final", "h_PFiso_PhotonBarrel_final", "h_PFiso_NeutralBarrel_final", "h_PFiso_SumBarrel_final", "h_HoEnewBarrel_final"};
	TString PlotsName4[nvar4] = {"CorrPFiso_Charged", "CorrPFiso_Neutral", "CorrPFiso_Photon", "PFiso_ChargedBarrel_final", "PFiso_PhotonBarrel_final", "PFiso_NeutralBarrel_final", "PFiso_SumBarrel_final", "HoEnewBarrel_final"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar4 ; ++i){

	    Double_t x_min, x_max;
	    if(i == 0){x_min = 0.0; x_max = 0.8; }
	    else if(i == 1){x_min = 0.0; x_max = 20.0;}
	    else if(i == 2  ){x_min = 0.0; x_max = 10.0;}
	    else if(i == 3 ){x_min = 0.0; x_max = 1.4;}
	    else if(i == 4){x_min = 0.0; x_max = 10.0;}
	    else if(i == 5){x_min = 0.0; x_max = 20.0;}
	    else if(i == 6){x_min = 0.0; x_max = 20.0;}
	    else if(i == 7){x_min = 0.0; x_max = 0.06;}

	    Double_t leg_xmin=0.7, leg_ymin=0.7, leg_xmax=0.9, leg_ymax=0.9;
	    Bool_t d_logy=true;
	    FillHistoNoRatio(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var4,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName4, SavePlots, d_logy);
//	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var4,
//		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
//		    scaleMC, Plotslocation, PlotsName4, SavePlots, d_logy);

	}
    }

    if(b_PhoID){ ///----->  5
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar5 = 4;
	TString var5[nvar5] = { "h_Photon_ecalRecHitSumEtConeDR04_final", "h_Photon_hcalTowerSumEtConeDR04_final", "h_Photon_trkSumPtHollowConeDR04_final", "h_HoE_final"};
	TString PlotsName5[nvar5] = {"Photon_ecalRecHitSumEtConeDR04_final", "Photon_hcalTowerSumEtConeDR04_final", "Photon_trkSumPtHollowConeDR04_final", "HoE_final"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar5 ; ++i){

	    Double_t x_min, x_max;
	    if(i == 0 || i == 1 || i == 2){x_min = 0.0; x_max = 15.0; }
	    else{x_min = 0.0; x_max = 0.06;}

	    Double_t leg_xmin=0.7, leg_ymin=0.7, leg_xmax=0.9, leg_ymax=0.9;
	    Bool_t d_logy=true;
//	    FillHistoNoRatio(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var5,
//		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
//		    scaleMC, Plotslocation, PlotsName5, SavePlots, d_logy);
	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var5,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName5, SavePlots, d_logy);

	}
    }

    if(b_JetID){ ///----->  6
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar6 = 6;
	TString var6[nvar6] = { "h_pfjet_NEF_final", "h_pfjet_NHF_final", "h_pfjet_CEF_final", "h_pfjet_CHF_final", "h_pfjet_NConstituents_final", "h_pfjet_ChargeMultiplicity_final"};
	TString PlotsName6[nvar6] = {"pfjet_NEF_final", "pfjet_NHF_final", "pfjet_CEF_final", "pfjet_CHF_final", "pfjet_NConstituents_final", "pfjet_ChargeMultiplicity_final"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar6 ; ++i){

	    Double_t x_min, x_max;
	    if(i == 0 || i == 1 || i == 2 || i == 3){x_min = 0.0; x_max = 1.0; }
	    else if(i == 4 || i == 5 ){x_min = 0.0; x_max = 50.0;}

//	    Double_t leg_xmin=0.7, leg_ymin=0.7, leg_xmax=0.9, leg_ymax=0.9;
	    Double_t leg_xmin=0.4, leg_ymin=0.1, leg_xmax=0.6, leg_ymax=0.3;
	    Bool_t d_logy=true;
	    FillHistoNoRatio(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var6,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName6, SavePlots, d_logy);
//	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var6,
//		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
//		    scaleMC, Plotslocation, PlotsName6, SavePlots, d_logy);

	}
    }

    if(b_ShowerShape){ ///----->  7
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar7 = 6;
	TString var7[nvar7] = { "h_Photon_SigmaIetaIeta_final", "h_Photon_SigmaIphiIphi_final", "h_Photon_SigmaEtaEta_final", "h_Photon_SigmaPhiPhi_final", "h_Photon_swissCross_final", "h_Photon_r9_final"};
	TString PlotsName7[nvar7] = {"Photon_SigmaIetaIeta_final", "Photon_SigmaIphiIphi_final", "Photon_SigmaEtaEta_final", "Photon_SigmaPhiPhi_final", "Photon_swissCross_final", "Photon_r9_final"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar7 ; ++i){

	    Double_t x_min, x_max;
	    Bool_t d_logy;
	    if(i == 0 ){ x_min = 0.0; x_max = 0.02; d_logy=false;}
	    else if(i == 1 || i == 2 || i == 3 ){ x_min = 0.0; x_max = 0.02; d_logy=false;}
	    else if(i == 4 ){x_min = -1.2; x_max = 1.2; d_logy=true;}
	    else if(i == 5 ){x_min = 0.6; x_max = 1.0; d_logy=true;}

	    Double_t leg_xmin=0.7, leg_ymin=0.7, leg_xmax=0.9, leg_ymax=0.9;
	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var7,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName7, SavePlots, d_logy);
	}
    }

    if(b_2DHist){ ///----->  8
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar8 = 9;
	TString var8[nvar8] = { "h_PtPhotJet_final", "h_dEtadPhi_final", "h_YPhotJet_final", "h_etaPhotJet_final", "h_phiPhotJet_final", "h_nIsoPhotonPt_noCut", "h_nIsoPhotonPt_final", "h_nIsoJetPt_noCut", "h_nIsoJetPt_final" };
	TString PlotsName8[nvar8] = {"PtPhotJet_final", "dEtadPhi_final", "YPhotJet_final", "etaPhotJet_final", "phiPhotJet_final", "nIsoPhotonPt_noCut", "nIsoPhotonPt_final", "nIsoJetPt_noCut", "nIsoJetPt_final"};
	TString PlotsNameExt8[nFiles] = {"Sig", "DiJet", "PhotonJet", "EWK", "Data"} ;

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar8 ; ++i){
	    for(Int_t j = 0 ; j < nFiles ; ++j){

		TFile* file8In = new TFile(FileName[j],"READ");
		TH2F *h8hist = (TH2F*)file8In->Get(var8[i]);

		Double_t x_min, x_max, y_min, y_max;
		if(i == 0 ){ x_min = 180.0; x_max = 2000.0; y_min = 180.0; y_max = 2000.0; }
		else if(i == 1 ){ x_min = 0.0; x_max = 3.0; y_min = 1.5; y_max = 3.2; }
		else if(i == 2 || i == 3){ x_min = -1.6; x_max = 1.6; y_min = -1.6; y_max = 1.6; }
		else if(i == 4 ){ x_min = 0.0; x_max = 6.4; y_min = 0.0; y_max = 6.4; }
		else if(i == 5 || i == 6  ){ x_min = 0.0; x_max = 2000.0; y_min = 0.0; y_max = 5.0; }
		else if(i == 7 || i == 8  ){ x_min = 0.0; x_max = 2000.0; y_min = 0.0; y_max = 20.0; }

		Fill2DHisto(i,j,h8hist,x_min,x_max,y_min,y_max,scaleMC,Plotslocation,PlotsName8,PlotsNameExt8,SavePlots);
		//		Fill2DHisto(i,j,h9hist,x_min,x_max,y_min,y_max,scaleMC,Plotslocation,PlotsName9,PlotsNameExt9,SavePlots,PlotsHeading);
	    }
	}
    }

    if(b_PC){ ///----->  9
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar9 = 4;
	TString var9[nvar9] = { "h_nPhoton_final", "h_nJet_final", "h_PC", "h_JC"};
	TString PlotsName9[nvar9] = {"nPhoton_final", "nJet_final", "PC", "JC"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar9 ; ++i){

	    Double_t x_min, x_max;
	    if(i == 0 ){ x_min = 0.0; x_max = 5.0;}
	    else if(i == 1 ){x_min = 0.0; x_max = 20.0;}
	    else if(i == 2 || i == 3){x_min = 0.0; x_max = 20.0;}

	    Double_t leg_xmin=0.7, leg_ymin=0.7, leg_xmax=0.9, leg_ymax=0.9;
	    Bool_t d_logy=true;
	    FillHisto(i, fileSig1, fileSig2, fileBkg1, fileBkg2, fileBkg3, fileData, var9,
		    x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax,
		    scaleMC, Plotslocation, PlotsName9, SavePlots, d_logy);

	}
    }


    if(b_CutFlow){ ///----->  10

	const Int_t nvar10 = 2;
	TString var10[nvar10] = {"h_CutFlowTable", "h_CutExpFlowTable"};
	TString PlotsName10[nvar10] = {"CutFlowTable", "CutExpFlowTable"};

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar10 ; ++i){
	    TH1F *h10Sig  = (TH1F*)fileSig1->Get(var10[i]);   //Signal
	    TH1F *h10Bkg1 = (TH1F*)fileBkg1->Get(var10[i]);   //DiJet MC
	    TH1F *h10Bkg2 = (TH1F*)fileBkg2->Get(var10[i]);   //GJ MC
	    TH1F *h10Bkg3 = (TH1F*)fileBkg3->Get(var10[i]);   //EWK
	    TH1F *h10Data = (TH1F*)fileData->Get(var10[i]);   //Data

	    GetCutFlowTables(i, h10Sig, h10Bkg1, h10Bkg2, h10Bkg3, h10Data, scaleMC, PlotsName10, Plotslocation);

	}
    }


    if(b_TurnOn){ ///----->  11
	const Int_t nvar11 = 4;
	TString var11[nvar11] = {"h_HLT_Photon_EB", "h_HLT_Photon150_EB", "h_HLT_Photon_MassEB", "h_HLT_Photon150_MassEB" };
	TString denoVar[nvar11] = {"h_HLT_Photon90_EB_all", "h_HLT_Photon90_EB_all", "h_HLT_Photon90_MassEB", "h_HLT_Photon90_MassEB" };

	TString PlotsName11[nvar11] = {"HLT_Photon_EB", "HLT_Photon150_EB", "HLT_Photon_MassEB", "HLT_Photon150_MassEB" };

	// Script for plotting and saving 
//	for(Int_t i = 0 ; i < nvar11 ; ++i){
	for(Int_t i = 1 ; i < 2 ; ++i){
	    TH1F *h11Data = (TH1F*)fileData->Get(var11[i]);   //Data
	    TH1F *h11deno = (TH1F*)fileData->Get(denoVar[i]);   //Data

	    Double_t x_min, x_max, y_min, y_max;
	    //x_min = 125; x_max = 500.0 ; y_min = 0.85 ; y_max = 1.005; 
	    x_min = 125; x_max = 1500.0 ; y_min = 0.85 ; y_max = 1.005; 

//	    DrawTurnOn(i,h11Data,h11deno,x_min,x_max,y_min,y_max,Plotslocation,PlotsName11,SavePlots);
	    DrawTurnOn(i,h11Data,h11deno,x_min,x_max,y_min,y_max,Plotslocation,PlotsName11,false);
	}
    }

    if(b_2DAllSignal){ ///----->  12
	// Number of variables to be plotted with their names +++++++++++++++++
	const Int_t nvar12 = 5;
	TString var12[nvar12] = { "h_PtPhotJet_final", "h_dEtadPhi_final", "h_YPhotJet_final", "h_etaPhotJet_final", "h_phiPhotJet_final" };
	TString PlotsName12[nvar12] = {"PtPhotJet_final", "dEtadPhi_final", "YPhotJet_final", "etaPhotJet_final", "phiPhotJet_final"};
	TString PlotsNameExt12[nSigFiles] = { "Qstar700", "Qstar1000", "Qstar1200", "Qstar1500", "Qstar1700", "Qstar2000", "Qstar2500", "Qstar3000", "Qstar3500" };

	// Script for plotting and saving 
	for(Int_t i = 0 ; i < nvar12 ; ++i){
	    for(Int_t j = 0 ; j < nSigFiles ; ++j){

		TFile* file12In = new TFile(FileSignals[j],"READ");
		TH2F *h12hist = (TH2F*)file12In->Get(var12[i]);

		Double_t x_min, x_max, y_min, y_max;
		if(i == 0 ){ x_min = 180.0; x_max = 2000.0; y_min = 180.0; y_max = 2000.0; }
		else if(i == 1 ){ x_min = 0.0; x_max = 3.0; y_min = 1.5; y_max = 3.2; }
		else if(i == 2 || i == 3){ x_min = -1.6; x_max = 1.6; y_min = -1.6; y_max = 1.6; }
		else if(i == 4 ){ x_min = 0.0; x_max = 6.4; y_min = 0.0; y_max = 6.4; }

		Fill2DHisto(i,j,h12hist,x_min,x_max,y_min,y_max,scaleMC,Plotslocation,PlotsName12,PlotsNameExt12,SavePlots);
		//		Fill2DHisto(i,j,h9hist,x_min,x_max,y_min,y_max,scaleMC,Plotslocation,PlotsName9,PlotsNameExt9,SavePlots,PlotsHeading);
	    }
	}
    }

    if(b_1DAllSignal){  ///---->13
	const Int_t nvar13 = 2;
	TString var13[nvar13] = { "h_dEta_final", "h_dphi_final" };
	TString PlotsName13[nvar13] = {"dEta_final", "dPhi_final"};
	TString PlotsNameExt13[nSigFiles] = { "Qstar700", "Qstar1000", "Qstar1200", "Qstar1500", "Qstar1700", "Qstar2000", "Qstar2500", "Qstar3000", "Qstar3500" };

	for(Int_t i = 0 ; i < nvar13 ; ++i){
	    for(Int_t j = 0 ; j < nSigFiles ; ++j){

		TFile* file13In = new TFile(FileSignals[j],"READ");
		TH1F *h13hist = (TH1F*)file13In->Get(var13[i]);

		TPaveText *cms = new TPaveText(0.2,0.9,0.8,0.93,"blNDC");
		cms->AddText("CMS Preliminary          #sqrt{s} = 8 TeV         #int Ldt=19.7 fb^{-1}");
		cms->SetBorderSize(0);
		cms->SetFillColor(0);
		cms->SetFillStyle(0);
		cms->SetTextFont(42);
		cms->SetTextSize(0.035);
		TCanvas* c2 = new TCanvas("c2","",80,20,700,700);
		c2->SetLogy(); c2->SetFillColor(0); c2->SetFrameBorderMode(0);
		gStyle->SetOptTitle(0);

		h13hist->SetStats(0);
		if(i == 0)h13hist->GetXaxis()->SetRangeUser(0.0,3.0);
		if(i == 1)h13hist->GetXaxis()->SetRangeUser(1.5,3.2);
		h13hist->Draw();
		cms->Draw();
		c2->SaveAs(Plotslocation+PlotsName13[i]+"_"+PlotsNameExt13[j]+".pdf");
		c2->SaveAs(Plotslocation+PlotsName13[i]+"_"+PlotsNameExt13[j]+".eps");

		delete c2;

	    }
	}
    }

    if(b_1DAllSignalImpose){  ///---->14
	const Int_t nvar14 = 1;
	TString var14[nvar14] = { "h_dEta_final" };
	TString PlotsName14[nvar14] = {"dEta_AllSignal"};
	TString PlotsNameExt14[nSigFiles] = { "Qstar700", "Qstar1000", "Qstar1200", "Qstar1500", "Qstar1700", "Qstar2000", "Qstar2500", "Qstar3000", "Qstar3500" };

	for(Int_t i = 0 ; i < nvar14 ; ++i){
	    TH1F *h14hist[nSigFiles];
	    TCanvas* c2 = new TCanvas("c2","",80,20,700,700);
	    //c2->SetLogy(); c2->SetFillColor(0); c2->SetFrameBorderMode(0);
	    TLegend *leg = new TLegend(0.7,0.6,0.9,0.9);
	    leg->SetFillColor(0); leg->SetTextFont(40); leg->SetBorderSize(1);
	    TPaveText *cms = new TPaveText(0.2,0.9,0.8,0.93,"blNDC");
	    cms->AddText("CMS Preliminary          #sqrt{s} = 8 TeV         #int Ldt=19.7 fb^{-1}");
	    cms->SetBorderSize(0);
	    cms->SetFillColor(0);
	    cms->SetFillStyle(0);
	    cms->SetTextFont(42);
	    cms->SetTextSize(0.035);
	    gStyle->SetOptTitle(0);

	    for(Int_t j = 0 ; j < nSigFiles ; ++j){

		TFile* file14In = new TFile(FileSignals[j],"READ");
		h14hist[j] = (TH1F*)file14In->Get(var14[i]);

		leg->AddEntry(h14hist[j],PlotsNameExt14[j],"L");

		h14hist[j]->SetStats(0);
		h14hist[j]->SetLineColor(j+1);
		h14hist[j]->SetLineStyle(j);
		h14hist[j]->SetLineWidth(2);
		h14hist[j]->Scale(1.0/h14hist[j]->Integral());
		if(i == 0)h14hist[j]->GetXaxis()->SetRangeUser(0.0,3.0);

		if(j == 0 ){h14hist[j]->Draw("HIST");}
		else{h14hist[j]->Draw("sameHIST");}

	    }
	    leg->Draw();
	    cms->Draw();
	    c2->SaveAs(Plotslocation+PlotsName14[i]+".pdf");
	    c2->SaveAs(Plotslocation+PlotsName14[i]+".eps");
	    delete c2;
	}
    }


if(b_Cummulative){  // ----------> 15

	const Int_t nvar15 = 6;
	TString var15[nvar15] = {"h_ptPhoton_finalPU", "h_ptJet_finalPU", "h_ptPhoton_massVarBin", "h_ptPhoton_VarBin", "h_ptJet_massVarBin", "h_ptJet_VarBin"};
	TString PlotsName15[nvar15] = {"ptPhoton_Cum", "ptJet_Cum", "ptPhoton_Cum_MVB", "ptPhoton_Cum_VB", "ptJet_Cum_MVB", "ptJet_Cum_VB"};

    for(int i = 0 ; i < nvar15 ;  ++i){

//	if(i == 0 || i == 1)double x_min = 190, x_max = 1000;
//	if(i == 2)double x_min = 480, x_max = 2000;
	double x_min = 190, x_max = 2000;
	Bool_t d_logy=true;

	Double_t leg_xmin=0.4, leg_ymin=0.1, leg_xmax=0.6, leg_ymax=0.3;
	
	DrawCummulativePlots(i, nFiles, fileSig1, fileBkg1, fileBkg2, fileBkg3, fileData, var15, 
		x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax, 
		scaleMC, Plotslocation, PlotsName15, SavePlots, d_logy)   ;
    }

}

if(b_ReverseCummulative){  // ----------> 16

	const Int_t nvar16 = 6;
	TString var16[nvar16] = {"h_ptPhoton_finalPU", "h_ptJet_finalPU", "h_ptPhoton_massVarBin", "h_ptPhoton_VarBin", "h_ptJet_massVarBin", "h_ptJet_VarBin"};
	TString PlotsName16[nvar16] = {"ptPhoton_ReCum", "ptJet_ReCum", "ptPhoton_ReCum_MVB", "ptPhoton_ReCum_VB", "ptJet_ReCum_MVB", "ptJet_ReCum_VB"};

    for(int i = 0 ; i < nvar16 ;  ++i){

//	if(i == 0 || i == 1)double x_min = 190, x_max = 1000;
//	if(i == 2)double x_min = 480, x_max = 2000;
	double x_min = 190, x_max = 2000;
	Bool_t d_logy=true;

	Double_t leg_xmin=0.4, leg_ymin=0.1, leg_xmax=0.6, leg_ymax=0.3;
	
	DrawReverseCummulativePlots(i, nFiles, fileSig1, fileBkg1, fileBkg2, fileBkg3, fileData, var16, 
		x_min, x_max, leg_xmin, leg_ymin, leg_xmax, leg_ymax, 
		scaleMC, Plotslocation, PlotsName16, SavePlots, d_logy)   ;
    }
}
//    SignalEfficiencies();
//    fHalf_SignalEfficiencies();


}

//void FillHisto(Int_t i, TH1F* hSig, TH1F* hBkg1, TH1F* hBkg2, TH1F* hData, Double_t x_min, Double_t x_max, const Float_t scaleMC,TString Plotslocation,TString *PlotsName, Bool_t SavePlots){ 
void FillHisto(Int_t i, TFile* fileSig1, TFile* fileSig2, TFile* fileBkg1, TFile* fileBkg2, TFile* fileBkg3, TFile* fileData, TString *var,
	Double_t x_min, Double_t x_max, Double_t leg_xmin, Double_t leg_ymin, Double_t leg_xmax, Double_t leg_ymax,
	const Float_t scaleMC, TString Plotslocation,TString *PlotsName, Bool_t SavePlots, Bool_t D_logy){ 

    //    const Float_t GJkfactor = 1.44;
    const Float_t GJkfactor = 1.4;

    //---------------------
    TH1F *hSig1 = (TH1F*)fileSig1->Get(var[i]);   //Signal 1 TeV
    TH1F *hSig2 = (TH1F*)fileSig2->Get(var[i]);   //Signal 3 TeV
    TH1F *hBkg1 = (TH1F*)fileBkg1->Get(var[i]);   //DiJet MC
    TH1F *hBkg2 = (TH1F*)fileBkg2->Get(var[i]);   //GJ MC
    TH1F *hBkg3 = (TH1F*)fileBkg3->Get(var[i]);   //EWK
    TH1F *hData = (TH1F*)fileData->Get(var[i]);   //Data
    ///------------------------------
    // Legends+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const Int_t nleg = 6 ;
    TString legendName[nleg] = { "q* = 1 TeV", "q* = 3 TeV", "DiJet (Bkg)", "GJ (Bkg)", "EWK (Bkg) ", "Data" };
    TString legendFill[nleg] = { "L", "L", "F", "F", "F", "p" };
    //---------------------------------------------------------------------

    //    TLegend *leg = new TLegend(0.72,0.76,0.9,0.9);
    TLegend *leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax);
    leg->SetFillColor(0); leg->SetTextFont(40); leg->SetBorderSize(1);
    leg->AddEntry(hSig1,legendName[0],legendFill[0]);
    leg->AddEntry(hSig2,legendName[1],legendFill[1]);
    leg->AddEntry(hBkg1,legendName[2],legendFill[2]);
    leg->AddEntry(hBkg2,legendName[3],legendFill[3]);
    leg->AddEntry(hBkg3,legendName[4],legendFill[4]);
    leg->AddEntry(hData,legendName[5],legendFill[5]);

    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
    cms->AddText("CMS Preliminary                #sqrt{s} = 8 TeV                #int Ldt=19.7 fb^{-1}");
    //    cms->AddText("#sqrt{s} = 8 TeV #int Ldt=19.7 fb^{-1}");
    cms->SetBorderSize(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(0);
    cms->SetTextFont(42);
    cms->SetTextSize(0.04);
    TPaveText *lumi = new TPaveText(0.1,0.9,0.3,0.95,"blNDC");
    lumi->AddText("CMS Preliminary");
    lumi->SetBorderSize(0);
    lumi->SetFillColor(0);
    lumi->SetFillStyle(0);
    lumi->SetTextFont(42);
    lumi->SetTextSize(0.04);

    hSig1->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    hSig1->GetXaxis()->SetRangeUser(x_min,x_max);
    hSig1->SetLineColor(kBlue+1);
    hSig1->SetLineWidth(2);

    hSig2->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    hSig2->GetXaxis()->SetRangeUser(x_min,x_max);
    hSig2->SetLineColor(kMagenta);
    hSig2->SetLineWidth(2);

    hBkg1->SetStats(0);// hBkg1->SetLineColor(kMagenta-3); hBkg1->SetFillColor(kOrange+8);
    hBkg1->GetXaxis()->SetRangeUser(x_min,x_max);
    hBkg1->SetLineColor(kGreen+4); 
    hBkg1->SetLineWidth(2);
    hBkg1->SetLineStyle(1);
    hBkg1->SetFillColor(kGreen-3);
    hBkg1->SetFillStyle(3018);
    hBkg1->Scale(scaleMC);

    hBkg2->SetStats(0); //hBkg2->SetLineColor(kCyan+2); hBkg2->SetFillColor(kGray) ;
    hBkg2->GetXaxis()->SetRangeUser(x_min,x_max);
    hBkg2->SetLineColor(kRed); 
    hBkg2->SetLineWidth(2);
    hBkg2->SetLineStyle(1);
    hBkg2->SetFillColor(kGray);
    hBkg2->Scale(scaleMC*GJkfactor);

    hBkg3->SetStats(0); //hBkg2->SetLineColor(kCyan+2); hBkg2->SetFillColor(kGray) ;
    hBkg3->GetXaxis()->SetRangeUser(x_min,x_max);
    hBkg3->SetLineColor(kBlack); 
    hBkg3->SetLineWidth(1);
    hBkg3->SetLineStyle(1);
    hBkg3->SetFillColor(kSpring-9);
    hBkg3->Scale(scaleMC);

    hData->SetStats(0); 
    hData->SetMarkerStyle(8); 
    hData->SetMarkerColor(kBlack); 
    hData->SetMarkerSize(0.8); 
    hData->SetLineColor(kBlack);
    hData->GetXaxis()->SetRangeUser(x_min,x_max);


    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    //    c1->SetLogy(); c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->Draw(); pad1->cd();
    if(D_logy){pad1->SetLogy();}
    pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
    pad1->SetBottomMargin(0.);
    THStack* hs = new THStack("hs","");

    hs->SetMinimum(0.01);
    hs->Add(hBkg3); 
    hs->Add(hBkg1);
    hs->Add(hBkg2); 
    hs->SetHistogram(hBkg2);

    TH1F *err_bkg = hBkg2->Clone();
    err_bkg->Add(hBkg1);
    err_bkg->Add(hBkg3);
    err_bkg->SetFillColor(1); err_bkg->SetFillStyle(3004);
    //  gStyle->SetHatchesSpacing(2.0); gStyle->SetHatchesLineWidth(3);err_bkg->SetMarkerStyle(1);

    hs->Draw("hist");
    hSig1->Draw("histSAME");
    hSig2->Draw("histSAME");
    hData->Draw("SAMEP");
    err_bkg->Draw("SAME E2");
    leg->Draw();
    cms->Draw();
    //    lumi->Draw();
    c1->cd();

    ///////////////////////////////////////////
    // Second PAD for making ratio plot
    TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
    pad2->Draw(); pad2->cd();
    pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    //pad2->SetLabelSize(0.04, "XYZ");
    TH1F *hBkgRatio = hBkg2->Clone();
    hBkgRatio->Add(hBkg1);
    hBkgRatio->Add(hBkg3);
    TH1F *hRatio = hData->Clone();
    //    hRatio->Divide(hBkgRatio);
    hRatio->Add(hBkgRatio,-1);
    hRatio->Divide(hData);

    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitle("(Data-MC)/Data");
    hRatio->GetYaxis()->SetTitleOffset(.4);
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetYaxis()->SetRangeUser(-1.0,1.0);

    hRatio->GetXaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->Draw("ep");
    
    TLine *ratioLine = new TLine(x_min, 0.0, x_max, 0.0);
    ratioLine->Draw("same");
    c1->cd();

    if(SavePlots){
	 c1->SaveAs(Plotslocation+PlotsName[i]+".pdf");
	 c1->SaveAs(Plotslocation+PlotsName[i]+".eps");
	//	c1->SaveAs(Plotslocation+PlotsName[i]+".gif");
	delete c1;
    }
}


void FillHistoNoRatio(Int_t i, TFile* fileSig1, TFile* fileSig2, TFile* fileBkg1, TFile* fileBkg2, TFile* fileBkg3, TFile* fileData, TString *var,
	Double_t x_min, Double_t x_max, Double_t leg_xmin, Double_t leg_ymin, Double_t leg_xmax, Double_t leg_ymax,
	const Float_t scaleMC, TString Plotslocation,TString *PlotsName, Bool_t SavePlots, Bool_t D_logy){ 

    //    const Float_t GJkfactor = 1.44;
    const Float_t GJkfactor = 1.4;

    //---------------------
    TH1F *hSig1 = (TH1F*)fileSig1->Get(var[i]);   //Signal 1 TeV
    TH1F *hSig2 = (TH1F*)fileSig2->Get(var[i]);   //Signal 3 TeV
    TH1F *hBkg1 = (TH1F*)fileBkg1->Get(var[i]);   //DiJet MC
    TH1F *hBkg2 = (TH1F*)fileBkg2->Get(var[i]);   //GJ MC
    TH1F *hBkg3 = (TH1F*)fileBkg3->Get(var[i]);   //EWK
    TH1F *hData = (TH1F*)fileData->Get(var[i]);   //Data
    ///------------------------------
    // Legends+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const Int_t nleg = 6 ;
    TString legendName[nleg] = { "q* = 1 TeV", "q* = 3 TeV", "DiJet (Bkg)", "GJ (Bkg)", "EWK (Bkg) ", "Data" };
    TString legendFill[nleg] = { "L", "L", "F", "F", "F", "p" };
    //---------------------------------------------------------------------

    //    TLegend *leg = new TLegend(0.72,0.76,0.9,0.9);
    TLegend *leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax);
    leg->SetFillColor(0); leg->SetTextFont(40); leg->SetBorderSize(1);
    leg->AddEntry(hSig1,legendName[0],legendFill[0]);
    leg->AddEntry(hSig2,legendName[1],legendFill[1]);
    leg->AddEntry(hBkg1,legendName[2],legendFill[2]);
    leg->AddEntry(hBkg2,legendName[3],legendFill[3]);
    leg->AddEntry(hBkg3,legendName[4],legendFill[4]);
    leg->AddEntry(hData,legendName[5],legendFill[5]);

//    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
//    cms->AddText("CMS Preliminary                #sqrt{s} = 8 TeV                #int Ldt=19.7 fb^{-1}");
    cms->AddText("CMS Preliminary     #sqrt{s} = 8 TeV     #int Ldt=19.7 fb^{-1}");
    //    cms->AddText("#sqrt{s} = 8 TeV #int Ldt=19.7 fb^{-1}");
    cms->SetBorderSize(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(0);
    cms->SetTextFont(42);
    cms->SetTextSize(0.04);
    TPaveText *lumi = new TPaveText(0.1,0.9,0.3,0.95,"blNDC");
    lumi->AddText("CMS Preliminary");
    lumi->SetBorderSize(0);
    lumi->SetFillColor(0);
    lumi->SetFillStyle(0);
    lumi->SetTextFont(42);
    lumi->SetTextSize(0.04);

    hSig1->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    hSig1->GetXaxis()->SetRangeUser(x_min,x_max);
    hSig1->SetLineColor(kBlue+1);
    hSig1->SetLineWidth(2);

    hSig2->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    hSig2->GetXaxis()->SetRangeUser(x_min,x_max);
    hSig2->SetLineColor(kMagenta);
    hSig2->SetLineWidth(2);

    hBkg1->SetStats(0);// hBkg1->SetLineColor(kMagenta-3); hBkg1->SetFillColor(kOrange+8);
    hBkg1->GetXaxis()->SetRangeUser(x_min,x_max);
    hBkg1->SetLineColor(kGreen+4); 
    hBkg1->SetLineWidth(2);
    hBkg1->SetLineStyle(1);
    hBkg1->SetFillColor(kGreen-3);
    hBkg1->SetFillStyle(3018);
    hBkg1->Scale(scaleMC);

    hBkg2->SetStats(0); //hBkg2->SetLineColor(kCyan+2); hBkg2->SetFillColor(kGray) ;
    hBkg2->GetXaxis()->SetRangeUser(x_min,x_max);
    hBkg2->SetLineColor(kRed); 
    hBkg2->SetLineWidth(2);
    hBkg2->SetLineStyle(1);
    hBkg2->SetFillColor(kGray);
    hBkg2->Scale(scaleMC*GJkfactor);

    hBkg3->SetStats(0); //hBkg2->SetLineColor(kCyan+2); hBkg2->SetFillColor(kGray) ;
    hBkg3->GetXaxis()->SetRangeUser(x_min,x_max);
    hBkg3->SetLineColor(kBlack); 
    hBkg3->SetLineWidth(1);
    hBkg3->SetLineStyle(1);
    hBkg3->SetFillColor(kSpring-9);
    hBkg3->Scale(scaleMC);

    hData->SetStats(0); 
    hData->SetMarkerStyle(8); 
    hData->SetMarkerColor(kBlack); 
    hData->SetMarkerSize(0.8); 
    hData->SetLineColor(kBlack);
    hData->GetXaxis()->SetRangeUser(x_min,x_max);


    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    if(D_logy)c1->SetLogy();

    THStack* hs = new THStack("hs","");
    hs->SetMinimum(0.01);
    hs->Add(hBkg3); 
    hs->Add(hBkg1);
    hs->Add(hBkg2); 
    hs->SetHistogram(hBkg2);

    TH1F *err_bkg = hBkg2->Clone();
    err_bkg->Add(hBkg1);
    err_bkg->Add(hBkg3);
    err_bkg->SetFillColor(1); err_bkg->SetFillStyle(3004);
    //  gStyle->SetHatchesSpacing(2.0); gStyle->SetHatchesLineWidth(3);err_bkg->SetMarkerStyle(1);

    hs->Draw("hist");
    hSig1->Draw("histSAME");
    hSig2->Draw("histSAME");
    hData->Draw("SAMEP");
    err_bkg->Draw("SAME E2");
    leg->Draw();
    cms->Draw();
    //    lumi->Draw();

    if(SavePlots){
	c1->SaveAs(Plotslocation+PlotsName[i]+".pdf");
	c1->SaveAs(Plotslocation+PlotsName[i]+".eps");
	//	c1->SaveAs(Plotslocation+PlotsName[i]+".gif");
	delete c1;
    }
}




void Fill2DHisto(Int_t i,Int_t j, TH2F* hist, Double_t x_min, Double_t x_max, Double_t y_min, Double_t y_max, const Float_t scaleMC,TString Plotslocation,TString *PlotsName,TString *PlotsNameExt, Bool_t SavePlots){ 
    //void Fill2DHisto(Int_t i,Int_t j, TH2F* hist, Double_t x_min, Double_t x_max, Double_t y_min, Double_t y_max, const Float_t scaleMC,TString Plotslocation,TString *PlotsName,TString *PlotsNameExt, Bool_t SavePlots, TString *PlotsHeading){ 

    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
    cms->AddText("CMS Preliminary     #sqrt{s} = 8 TeV      #int Ldt=19.7 fb^{-1}");
    cms->SetBorderSize(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(0);
    cms->SetTextFont(42);
    cms->SetTextSize(0.04);
    TPaveText *name = new TPaveText(0.4,0.85,0.6,0.9,"blNDC");
    name->AddText(PlotsNameExt[j]);
    //    name->AddText(PlotsHeading[i]);
    name->SetBorderSize(0);
    name->SetFillColor(0);
    name->SetFillStyle(0);
    name->SetTextFont(42);
    name->SetTextSize(0.04);

    hist->SetStats(0);
    hist->GetXaxis()->SetRangeUser(x_min,x_max);
    hist->GetYaxis()->SetRangeUser(y_min,y_max);

    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    c1->SetLogz(); 
    c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    gStyle->SetOptTitle(0);
    hist->Draw("COLZ");
    cms->Draw();
    name->Draw();

    if(SavePlots){
	c1->SaveAs(Plotslocation+PlotsName[i]+"_"+PlotsNameExt[j]+".pdf");
	c1->SaveAs(Plotslocation+PlotsName[i]+"_"+PlotsNameExt[j]+".eps");
	//	c1->SaveAs(Plotslocation+PlotsName[i]+"_"+PlotsNameExt[j]+".gif");
	delete c1;
    }
}
void DrawCummulativePlots(Int_t i,const Int_t nFiles, TFile* fileSig1, TFile* fileBkg1, TFile* fileBkg2, TFile* fileBkg3, TFile* fileData, TString *var,
	Double_t x_min, Double_t x_max, Double_t leg_xmin, Double_t leg_ymin, Double_t leg_xmax, Double_t leg_ymax,
	const Float_t scaleMC, TString Plotslocation,TString *PlotsName, Bool_t SavePlots, Bool_t D_logy){ 

    const Float_t GJkfactor = 1.4;

    TH1F *h_hist[nFiles];
    TH1F *h_Cum[nFiles];

    Int_t nnbins[nFiles];
    Double_t total[nFiles] = { 0, 0, 0, 0, 0 };
 
    //---------------------
    h_hist[0] = (TH1F*)fileSig1->Get(var[i]);   //Signal 1 TeV
    h_hist[1] = (TH1F*)fileBkg1->Get(var[i]);   //DiJet MC
    h_hist[2] = (TH1F*)fileBkg2->Get(var[i]);   //GJ MC
    h_hist[3] = (TH1F*)fileBkg3->Get(var[i]);   //EWK
    h_hist[4] = (TH1F*)fileData->Get(var[i]);   //Data
    
    for(Int_t j = 0 ; j < nFiles ; ++j){
	h_Cum[j] = (TH1F*)h_hist[j]->Clone("h_Cum[j]");
	h_Cum[j]->Reset();

	nnbins[j] = h_hist[j]->GetNbinsX();
    }

    for(int j = 0 ; j < nFiles ; ++j){
	for(int k = 0 ; k < nnbins[j] ; ++k){
	    total[j] += h_hist[j]->GetBinContent(k);
	    h_Cum[j]->SetBinContent(k, total[j]) ;
	}
    }

    // Legends+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const Int_t nleg = 5 ;
    TString legendName[nleg] = { "q* = 1 TeV", "DiJet (Bkg)", "GJ (Bkg)", "EWK (Bkg) ", "Data" };
    TString legendFill[nleg] = { "L", "F", "F", "F", "p" };
    //---------------------------------------------------------------------

    //    TLegend *leg = new TLegend(0.72,0.76,0.9,0.9);
    TLegend *leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax);
    leg->SetFillColor(0); leg->SetTextFont(40); leg->SetBorderSize(1);
    leg->AddEntry(h_Cum[0],legendName[0],legendFill[0]);
    leg->AddEntry(h_Cum[2],legendName[2],legendFill[2]);
    leg->AddEntry(h_Cum[1],legendName[1],legendFill[1]);
    leg->AddEntry(h_Cum[3],legendName[3],legendFill[3]);
    leg->AddEntry(h_Cum[4],legendName[4],legendFill[4]);

    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
    cms->AddText("CMS Preliminary                #sqrt{s} = 8 TeV                #int Ldt=19.7 fb^{-1}");
    cms->SetBorderSize(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(0);
    cms->SetTextFont(42);
    cms->SetTextSize(0.04);

    h_Cum[0]->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    h_Cum[0]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[0]->SetLineColor(kBlue+1);
    h_Cum[0]->SetLineWidth(2);

    h_Cum[2]->SetStats(0);// h_Cum[0]->SetLineColor(kMagenta-3); h_Cum[0]->SetFillColor(kOrange+8);
    h_Cum[2]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[2]->SetLineColor(kGreen+4); 
    h_Cum[2]->SetLineWidth(2);
    h_Cum[2]->SetLineStyle(1);
    h_Cum[2]->SetFillColor(kGreen-3);
    h_Cum[2]->SetFillStyle(3018);
    h_Cum[2]->Scale(scaleMC*GJkfactor);

    h_Cum[1]->SetStats(0); //h_Cum[1]->SetLineColor(kCyan+2); h_Cum[1]->SetFillColor(kGray) ;
    h_Cum[1]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[1]->SetLineColor(kRed); 
    h_Cum[1]->SetLineWidth(2);
    h_Cum[1]->SetLineStyle(1);
    h_Cum[1]->SetFillColor(kGray);
    h_Cum[1]->Scale(scaleMC);

    h_Cum[3]->SetStats(0); //h_Cum[1]->SetLineColor(kCyan+2); h_Cum[1]->SetFillColor(kGray) ;
    h_Cum[3]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[3]->SetLineColor(kBlack); 
    h_Cum[3]->SetLineWidth(1);
    h_Cum[3]->SetLineStyle(1);
    h_Cum[3]->SetFillColor(kSpring-9);
    h_Cum[3]->Scale(scaleMC);

    h_Cum[4]->SetStats(0); 
    h_Cum[4]->SetMarkerStyle(8); 
    h_Cum[4]->SetMarkerColor(kBlack); 
    h_Cum[4]->SetMarkerSize(0.8); 
    h_Cum[4]->SetLineColor(kBlack);
    h_Cum[4]->GetXaxis()->SetRangeUser(x_min,x_max);

    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    //    c1->SetLogy(); c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->Draw(); pad1->cd();
    if(D_logy){pad1->SetLogy();}
    pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
    pad1->SetBottomMargin(0.);

    THStack* hs = new THStack("hs","");
    //    hs->SetMinimum(0.01);
    hs->Add(h_Cum[3]); 
    hs->Add(h_Cum[1]);
    hs->Add(h_Cum[2]); 
    hs->SetHistogram(h_Cum[2]);

    TH1F *err_bkg = h_Cum[2]->Clone();
    err_bkg->Add(h_Cum[1]);
    err_bkg->Add(h_Cum[3]);
    err_bkg->SetFillColor(1); err_bkg->SetFillStyle(3004);
    //  gStyle->SetHatchesSpacing(2.0); gStyle->SetHatchesLineWidth(3);err_bkg->SetMarkerStyle(1);

    hs->Draw("hist");
    h_Cum[0]->Draw("histSAME");
    h_Cum[4]->Draw("SAMEP");
//    err_bkg->Draw("SAME E2");
    leg->Draw();
    cms->Draw();
    //    lumi->Draw();
    c1->cd();

    ///////////////////////////////////////////
    // Second PAD for making ratio plot
    TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
    pad2->Draw(); pad2->cd();
    pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    //pad2->SetLabelSize(0.04, "XYZ");
    TH1F *hBkgRatio = h_Cum[2]->Clone();
    hBkgRatio->Add(h_Cum[1]);
    hBkgRatio->Add(h_Cum[3]);
    TH1F *hRatio = h_Cum[4]->Clone();
    //    hRatio->Divide(hBkgRatio);
    hRatio->Add(hBkgRatio,-1);
    hRatio->Divide(h_Cum[4]);

    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitle("(Data-MC)/Data");
    hRatio->GetYaxis()->SetTitleOffset(.4);
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetYaxis()->SetRangeUser(-1.0,1.0);

    hRatio->GetXaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->Draw("ep");
    c1->cd();

    if(SavePlots){
	c1->SaveAs(Plotslocation+PlotsName[i]+".pdf");
	c1->SaveAs(Plotslocation+PlotsName[i]+".eps");
	//	c1->SaveAs(Plotslocation+PlotsName[i]+".gif");
	delete c1;
    }
}


void DrawReverseCummulativePlots(Int_t i,const Int_t nFiles, TFile* fileSig1, TFile* fileBkg1, TFile* fileBkg2, TFile* fileBkg3, TFile* fileData, TString *var,
	Double_t x_min, Double_t x_max, Double_t leg_xmin, Double_t leg_ymin, Double_t leg_xmax, Double_t leg_ymax,
	const Float_t scaleMC, TString Plotslocation,TString *PlotsName, Bool_t SavePlots, Bool_t D_logy){ 

    const Float_t GJkfactor = 1.4;

    TH1F *h_hist[nFiles];
    TH1F *h_Cum[nFiles];

    Int_t nnbins[nFiles];
    Double_t total[nFiles] = { 0, 0, 0, 0, 0 };
 
    //---------------------
    h_hist[0] = (TH1F*)fileSig1->Get(var[i]);   //Signal 1 TeV
    h_hist[1] = (TH1F*)fileBkg1->Get(var[i]);   //DiJet MC
    h_hist[2] = (TH1F*)fileBkg2->Get(var[i]);   //GJ MC
    h_hist[3] = (TH1F*)fileBkg3->Get(var[i]);   //EWK
    h_hist[4] = (TH1F*)fileData->Get(var[i]);   //Data
    
    for(Int_t j = 0 ; j < nFiles ; ++j){
	h_Cum[j] = (TH1F*)h_hist[j]->Clone("h_Cum[j]");
	h_Cum[j]->Reset();

	nnbins[j] = h_hist[j]->GetNbinsX();
    }

    for(int j = 0 ; j < nFiles ; ++j){
	for(int k = nnbins[j]-1 ; k >= 0 ; --k){
	    total[j] += h_hist[j]->GetBinContent(k);
	    h_Cum[j]->SetBinContent(k, total[j]) ;
	}
    }

    // Legends+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const Int_t nleg = 5 ;
    TString legendName[nleg] = { "q* = 1 TeV", "DiJet (Bkg)", "GJ (Bkg)", "EWK (Bkg) ", "Data" };
    TString legendFill[nleg] = { "L", "F", "F", "F", "p" };
    //---------------------------------------------------------------------

    //    TLegend *leg = new TLegend(0.72,0.76,0.9,0.9);
    TLegend *leg = new TLegend(leg_xmin,leg_ymin,leg_xmax,leg_ymax);
    leg->SetFillColor(0); leg->SetTextFont(40); leg->SetBorderSize(1);
    leg->AddEntry(h_Cum[0],legendName[0],legendFill[0]);
    leg->AddEntry(h_Cum[2],legendName[2],legendFill[2]);
    leg->AddEntry(h_Cum[1],legendName[1],legendFill[1]);
    leg->AddEntry(h_Cum[3],legendName[3],legendFill[3]);
    leg->AddEntry(h_Cum[4],legendName[4],legendFill[4]);

    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
    cms->AddText("CMS Preliminary                #sqrt{s} = 8 TeV                #int Ldt=19.7 fb^{-1}");
    cms->SetBorderSize(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(0);
    cms->SetTextFont(42);
    cms->SetTextSize(0.04);

    h_Cum[0]->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    h_Cum[0]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[0]->SetLineColor(kBlue+1);
    h_Cum[0]->SetLineWidth(2);

    h_Cum[2]->SetStats(0);// h_Cum[0]->SetLineColor(kMagenta-3); h_Cum[0]->SetFillColor(kOrange+8);
    h_Cum[2]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[2]->SetLineColor(kGreen+4); 
    h_Cum[2]->SetLineWidth(2);
    h_Cum[2]->SetLineStyle(1);
    h_Cum[2]->SetFillColor(kGreen-3);
    h_Cum[2]->SetFillStyle(3018);
    h_Cum[2]->Scale(scaleMC*GJkfactor);

    h_Cum[1]->SetStats(0); //h_Cum[1]->SetLineColor(kCyan+2); h_Cum[1]->SetFillColor(kGray) ;
    h_Cum[1]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[1]->SetLineColor(kRed); 
    h_Cum[1]->SetLineWidth(2);
    h_Cum[1]->SetLineStyle(1);
    h_Cum[1]->SetFillColor(kGray);
    h_Cum[1]->Scale(scaleMC);

    h_Cum[3]->SetStats(0); //h_Cum[1]->SetLineColor(kCyan+2); h_Cum[1]->SetFillColor(kGray) ;
    h_Cum[3]->GetXaxis()->SetRangeUser(x_min,x_max);
    h_Cum[3]->SetLineColor(kBlack); 
    h_Cum[3]->SetLineWidth(1);
    h_Cum[3]->SetLineStyle(1);
    h_Cum[3]->SetFillColor(kSpring-9);
    h_Cum[3]->Scale(scaleMC);

    h_Cum[4]->SetStats(0); 
    h_Cum[4]->SetMarkerStyle(8); 
    h_Cum[4]->SetMarkerColor(kBlack); 
    h_Cum[4]->SetMarkerSize(0.8); 
    h_Cum[4]->SetLineColor(kBlack);
    h_Cum[4]->GetXaxis()->SetRangeUser(x_min,x_max);

    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    //    c1->SetLogy(); c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->Draw(); pad1->cd();
    if(D_logy){pad1->SetLogy();}
    pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
    pad1->SetBottomMargin(0.);

    THStack* hs = new THStack("hs","");
    //    hs->SetMinimum(0.01);
    hs->Add(h_Cum[3]); 
    hs->Add(h_Cum[1]);
    hs->Add(h_Cum[2]); 
    hs->SetHistogram(h_Cum[2]);

    TH1F *err_bkg = h_Cum[2]->Clone();
    err_bkg->Add(h_Cum[1]);
    err_bkg->Add(h_Cum[3]);
    err_bkg->SetFillColor(1); err_bkg->SetFillStyle(3004);
    //  gStyle->SetHatchesSpacing(2.0); gStyle->SetHatchesLineWidth(3);err_bkg->SetMarkerStyle(1);

    hs->Draw("hist");
    h_Cum[0]->Draw("histSAME");
    h_Cum[4]->Draw("SAMEP");
//    err_bkg->Draw("SAME E2");
    leg->Draw();
    cms->Draw();
    //    lumi->Draw();
    c1->cd();

    ///////////////////////////////////////////
    // Second PAD for making ratio plot
    TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.25);
    pad2->Draw(); pad2->cd();
    pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    //pad2->SetLabelSize(0.04, "XYZ");
    TH1F *hBkgRatio = h_Cum[2]->Clone();
    hBkgRatio->Add(h_Cum[1]);
    hBkgRatio->Add(h_Cum[3]);
    TH1F *hRatio = h_Cum[4]->Clone();
    //    hRatio->Divide(hBkgRatio);
    hRatio->Add(hBkgRatio,-1);
    hRatio->Divide(h_Cum[4]);

    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitle("(Data-MC)/Data");
    hRatio->GetYaxis()->SetTitleOffset(.4);
    hRatio->GetYaxis()->SetTitleSize(0.1);
    hRatio->GetYaxis()->SetRangeUser(-1.0,1.0);

    hRatio->GetXaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->Draw("ep");
    c1->cd();

    if(SavePlots){
	c1->SaveAs(Plotslocation+PlotsName[i]+".pdf");
	c1->SaveAs(Plotslocation+PlotsName[i]+".eps");
	//	c1->SaveAs(Plotslocation+PlotsName[i]+".gif");
	delete c1;
    }
}

void GetCutFlowTables(Int_t i, TH1F* hSig, TH1F* hBkg1, TH1F* hBkg2, TH1F* hBkg3, TH1F* hData, const Float_t scaleMC, TString *PlotsName, TString Plotslocation){

    Int_t nbins = hSig->GetNbinsX();
    //    const Float_t GJkfactor = 1.44;
    const Float_t GJkfactor = 1.4;

    std::vector<Double_t> ExpEntries_Sig, ExpEntries_Bkg1, ExpEntries_Bkg2, ExpEntries_Bkg3, ExpEntries_MC, ExpEntries_Data;
    std::vector<TString> Labels;

    ofstream fout(PlotsName[i]+=".tex");
    //    fout<<"\\rowcolor{dkgreen}\\color{white}{Cuts}     &  \\color{white}{Signal}  &   \\color{white}{DiJet}   &    \\color{white}{PhotonJet}  &   \\color{white}{EWK}  & \\color{white}{Total Bkg}   &    \\color{white}{Data}  & \\color{white}{S/$\sqrt$B} \\\\"<<endl<<"\\hline"<<endl;
    fout<<"Cuts     &  Signal  &   DiJetg   &    PhotonJet   &   EWK   &   Total Bkg   &    Data    \\\\"<<endl<<"\\hline"<<endl;

    for(Int_t k = 1; k <= nbins ; ++k){
	ExpEntries_Sig.push_back(hSig->GetBinContent(k));
	ExpEntries_Bkg1.push_back(hBkg1->GetBinContent(k));   // DiJet
	ExpEntries_Bkg2.push_back(hBkg2->GetBinContent(k));   // GJ
	ExpEntries_Bkg3.push_back(hBkg3->GetBinContent(k));   // EWK
//	ExpEntries_MC.push_back(hBkg1->GetBinContent(k)+(GJkfactor*(hBkg2->GetBinContent(k)))+hBkg3->GetBinContent(k));
	ExpEntries_MC.push_back((1.34*(hBkg1->GetBinContent(k)))+(1.33*(hBkg2->GetBinContent(k)))+hBkg3->GetBinContent(k));
	ExpEntries_Data.push_back(hData->GetBinContent(k));
	Labels.push_back(hData->GetXaxis()->GetBinLabel(k));

	//	fout<<"\\color{blue}"<<Labels[k-1]<<" &  "<<ExpEntries_Sig[k-1]<<" &  "<<scaleMC*ExpEntries_Bkg1[k-1]<<" &  "<<scaleMC*GJkfactor*ExpEntries_Bkg2[k-1]<<" &  "<<scaleMC*ExpEntries_Bkg3[k-1]<<" &  "<<scaleMC*ExpEntries_MC[k-1]<<" &  "<<ExpEntries_Data[k-1]<<" & "<<ExpEntries_Sig[k-1]/(pow(ExpEntries_MC[k-1],0.5))<<" \\"<<"\\"<<endl;
//	fout<<Labels[k-1]<<" &  "<<ExpEntries_Sig[k-1]<<" &  "<<scaleMC*ExpEntries_Bkg1[k-1]<<" &  "<<scaleMC*GJkfactor*ExpEntries_Bkg2[k-1]<<" &  "<<scaleMC*ExpEntries_Bkg3[k-1]<<" &  "<<scaleMC*ExpEntries_MC[k-1]<<" &  "<<ExpEntries_Data[k-1]<<" \\"<<"\\"<<endl;
	fout<<Labels[k-1]<<" &  "<<ExpEntries_Sig[k-1]<<" &  "<<scaleMC*1.34*ExpEntries_Bkg1[k-1]<<" &  "<<scaleMC*1.33*ExpEntries_Bkg2[k-1]<<" &  "<<scaleMC*ExpEntries_Bkg3[k-1]<<" &  "<<scaleMC*ExpEntries_MC[k-1]<<" &  "<<ExpEntries_Data[k-1]<<" \\"<<"\\"<<endl;

    }
    fout.close();


    TH1F *hTotalBkg = hBkg1->Clone();
    hTotalBkg->Add(hBkg2,GJkfactor);
    hTotalBkg->Add(hBkg3);
    hTotalBkg->Scale(scaleMC);
    //    hTotalBkg->Scale(1.0/hTotalBkg->Integral());
    //    hData->Scale(1.0/hData->Integral());
    //    hSig->Scale(1.0/hSig->Integral());

    hData->SetMarkerStyle(20);
    //    hData->SetFillColor(kGray+1);
    hData->SetFillColor(kBlue+1);
    hSig->SetMarkerStyle(33);
    hSig->SetMarkerSize(1.5);
    hSig->SetMarkerColor(kYellow);
    hTotalBkg->SetFillColor(kOrange+7);
    hTotalBkg->SetMarkerColor(kRed);

    ///----------------------------------
    TLegend *leg = new TLegend(0.75,0.7,0.9,0.9);
    leg->SetFillColor(0); leg->SetTextFont(40); leg->SetBorderSize(1);
    leg->AddEntry(hSig,"Signal (q*)","p");
    leg->AddEntry(hTotalBkg,"Total Bkg","F");
    leg->AddEntry(hData,"Data","F");

    hData->GetYaxis()->SetRangeUser(10e2,10e9);

    hData->SetBarWidth(0.4); hData->SetBarOffset(0.1); hData->SetStats(0);
    hData->Draw("b text90");
    hTotalBkg->SetBarWidth(0.4); hTotalBkg->SetBarOffset(0.5); hTotalBkg->SetStats(0);
    hTotalBkg->Draw("b text90 same ");
    hSig->Draw("sameP");
    leg->Draw();
    ////------------------------------
    /*
       hTotalBkg->SetStats(0);
       hTotalBkg->SetMinimum(100);
       hTotalBkg->Draw("HISTtext");
       hData->Draw("sameP");
       hSig->Draw("sameP");
     */
    c1->SetLogy();
    c1->SaveAs(Plotslocation+"Eff_"+PlotsName[i]+".pdf");
    c1->SaveAs(Plotslocation+"Eff_"+PlotsName[i]+".eps");

}

void DrawTurnOn(Int_t i,TH1F* hNum,TH1F* hDeno, Double_t x_min, Double_t x_max, Double_t y_min, Double_t y_max,TString Plotslocation,TString *PlotsName, Bool_t SavePlots){

    TCanvas * c1 = new TCanvas("c1","",80,20,700,700);
    c1->SetFillColor(0); c1->SetGridy(); c1->SetGridx(); c1->SetFrameBorderMode(0);

    TPaveText *cms = new TPaveText(0.1,0.9,0.9,0.95,"blNDC");
    cms->AddText("CMS Preliminary     #sqrt{s} = 8 TeV      #int Ldt=19.7 fb^{-1}");
    cms->SetBorderSize(0);
    cms->SetFillColor(0);
    cms->SetFillStyle(0);
    cms->SetTextFont(42);
    cms->SetTextSize(0.04);
    TPaveText *name = new TPaveText(0.4,0.85,0.6,0.9,"blNDC");
    name->AddText(PlotsName[i]);
    name->SetBorderSize(0);
    name->SetFillColor(0);
    name->SetFillStyle(0);
    name->SetTextFont(42);
    name->SetTextSize(0.04);

    TArrow *arrow = new TArrow(199.7307,0.9565306,199.7307,0.997449,0.02,">");
    arrow->SetFillColor(1);
    arrow->SetFillStyle(1001);
    arrow->SetLineColor(2);
    arrow->SetLineWidth(2);

    TGraphAsymmErrors* hlt_gr = new TGraphAsymmErrors();
    hlt_gr->Divide(hNum, hDeno);
    hlt_gr->SetMarkerStyle(8);
    hlt_gr->SetMarkerColor(2);
    hlt_gr->SetMarkerSize(.6);

    hlt_gr->GetXaxis()->SetLimits(x_min,x_max);
    hlt_gr->GetXaxis()->SetTitle("P_{T}^{#gamma} (GeV)");
    hlt_gr->GetXaxis()->CenterTitle();;
    hlt_gr->GetHistogram()->SetMinimum(y_min);
    hlt_gr->GetHistogram()->SetMaximum(y_max);
    hlt_gr->Draw("AP");
    cms->Draw();
    name->Draw();
    if(i == 2)arrow->Draw();


    if(SavePlots){
	c1->SaveAs(Plotslocation+PlotsName[i]+".pdf");
	c1->SaveAs(Plotslocation+PlotsName[i]+".eps");
	//	c1->SaveAs(Plotslocation+PlotsName[i]+".gif");
	delete c1;
    }

}


SignalEfficiencies(){

    const Int_t nSigFiles = 9;
    Int_t NoCut, Final;

    TString FileSignals[nSigFiles] = { 
	"./rootFiles/QstarToGJ_M_700.root",
	"./rootFiles/QstarToGJ_M_1000.root",
	"./rootFiles/QstarToGJ_M_1200.root",
	"./rootFiles/QstarToGJ_M_1500.root",
	"./rootFiles/QstarToGJ_M_1700.root",
	"./rootFiles/QstarToGJ_M_2000.root",
	"./rootFiles/QstarToGJ_M_2500.root",
	"./rootFiles/QstarToGJ_M_3000.root",
	"./rootFiles/QstarToGJ_M_3500.root"
    };


    const Int_t nvar = 2;
    TString var[nvar] = {"h_CutFlowTable", "h_CutExpFlowTable"};
    TString PlotsName[nvar] = {"SignalEffCutFlow", "SignalEffCutExpFlow"};
    TH1F *hSig[nSigFiles];

    // Script for plotting and saving 
    //    for(Int_t i = 1 ; i < nvar ; ++i){
    for(Int_t i = 0 ; i < 1 ; ++i){
	for(Int_t j = 0 ; j < nSigFiles ; ++j){
	    TFile* fileIn = new TFile(FileSignals[j],"READ");
	    hSig[j]  = (TH1F*)fileIn->Get(var[i]);   //Signal
	}
    }

    Int_t nbins = hSig[0]->GetNbinsX();

    std::vector<Double_t> ExpEntries_Sig700, ExpEntries_Sig1000, ExpEntries_Sig1200, ExpEntries_Sig1500, ExpEntries_Sig1700, ExpEntries_Sig2000, ExpEntries_Sig2500, ExpEntries_Sig3000, ExpEntries_Sig3500;
    std::vector<TString> Labels;

    ofstream fout(".tex");
    fout<<"Cuts     &  700  & 1000 & 1200 & 1500 & 1700 & 2000 & 2500 & 3000 & 3500  \\\\"<<endl<<"\\hline \\hline"<<endl;

    for(Int_t k = 1; k <= nbins ; ++k){ 
	ExpEntries_Sig700.push_back(hSig[0]->GetBinContent(k));
	ExpEntries_Sig1000.push_back(hSig[1]->GetBinContent(k));
	ExpEntries_Sig1200.push_back(hSig[2]->GetBinContent(k));
	ExpEntries_Sig1500.push_back(hSig[3]->GetBinContent(k));
	ExpEntries_Sig1700.push_back(hSig[4]->GetBinContent(k));
	ExpEntries_Sig2000.push_back(hSig[5]->GetBinContent(k));
	ExpEntries_Sig2500.push_back(hSig[6]->GetBinContent(k));
	ExpEntries_Sig3000.push_back(hSig[7]->GetBinContent(k));
	ExpEntries_Sig3500.push_back(hSig[8]->GetBinContent(k));

	Labels.push_back(hSig[0]->GetXaxis()->GetBinLabel(k));

//	fout<<Labels[k-1]<<" & "<<ExpEntries_Sig700[k-1]<<" &  "<<ExpEntries_Sig1000[k-1]<<" &  "<<ExpEntries_Sig1200[k-1]<<" &  "<<ExpEntries_Sig1500[k-1]<<" &  "<<ExpEntries_Sig1700[k-1]<<" &  "<<ExpEntries_Sig2000[k-1]<<" &  "<<ExpEntries_Sig2500[k-1]<<" &  "<<ExpEntries_Sig3000[k-1]<<" &  "<<ExpEntries_Sig3500[k-1]<<" \\"<<"\\"<<endl;
	fout<<Labels[k-1]<<" & "<<100*ExpEntries_Sig700[k-1]/ExpEntries_Sig700[0]<<" &  "<<100*ExpEntries_Sig1000[k-1]/ExpEntries_Sig1000[0]<<" &  "<<ExpEntries_Sig1200[k-1]<<" &  "<<ExpEntries_Sig1500[k-1]<<" &  "<<ExpEntries_Sig1700[k-1]<<" &  "<<ExpEntries_Sig2000[k-1]<<" &  "<<ExpEntries_Sig2500[k-1]<<" &  "<<ExpEntries_Sig3000[k-1]<<" &  "<<ExpEntries_Sig3500[k-1]<<" \\"<<"\\"<<endl;


	if(Labels[k-1] == "Total") NoCut = k-1;
	if(Labels[k-1] == "MassCut") Final = k-1;
    }

    //std::cout<<NoCut<<"  "<<Final<<std::endl;
    fout<<"\\hline"<<std::endl;
    fout<<"Efficiency" <<" & "<<ExpEntries_Sig700[Final]/ExpEntries_Sig700[NoCut]<<" &  "<<ExpEntries_Sig1000[Final]/ExpEntries_Sig1000[NoCut]<<" &  "<<ExpEntries_Sig1200[Final]/ExpEntries_Sig1200[NoCut]<<" &  "<<ExpEntries_Sig1500[Final]/ExpEntries_Sig1500[NoCut]<<" &  "<<ExpEntries_Sig1700[Final]/ExpEntries_Sig1700[NoCut]<<" &  "<<ExpEntries_Sig2000[Final]/ExpEntries_Sig2000[NoCut]<<" &  "<<ExpEntries_Sig2500[Final]/ExpEntries_Sig2500[NoCut]<<" &  "<<ExpEntries_Sig3000[Final]/ExpEntries_Sig3000[NoCut]<<" &  "<<ExpEntries_Sig3500[Final]/ExpEntries_Sig3500[NoCut]<<" \\"<<"\\"<<endl;
    fout<<"\\hline"<<std::endl;

    fout.close();
}


fHalf_SignalEfficiencies(){

    const Int_t nSigFiles = 6;

    std::vector<TString> FileSignals; 

    TString FileSignals[nSigFiles] = { 
	"./rootFiles/QstarToGJ_M_fhalf_1000.root",
	"./rootFiles/QstarToGJ_M_fhalf_1500.root",
	"./rootFiles/QstarToGJ_M_fhalf_2000.root",
	"./rootFiles/QstarToGJ_M_fhalf_2500.root",
	"./rootFiles/QstarToGJ_M_fhalf_3000.root",
	"./rootFiles/QstarToGJ_M_fhalf_3500.root"
    };

    const Int_t nvar = 2;
    TString var[nvar] = {"h_CutFlowTable", "h_CutExpFlowTable"};
    TString PlotsName[nvar] = {"fHalf_SignalEffCutFlow", "fHalf_SignalEffCutExpFlow"};
    TH1F *hSig[nSigFiles];

    // Script for plotting and saving 
    for(Int_t i = 0 ; i < nvar ; ++i){
	for(Int_t j = 0 ; j < nSigFiles ; ++j){
	    TFile* fileIn = new TFile(FileSignals[j],"READ");
	    hSig[j]  = (TH1F*)fileIn->Get(var[i]);   //Signal
	}
    }

    Int_t nbins = hSig[0]->GetNbinsX();

    std::vector<Double_t> ExpEntries_Sig1000, ExpEntries_Sig1500, ExpEntries_Sig2000, ExpEntries_Sig2500, ExpEntries_Sig3000, ExpEntries_Sig3500;
    std::vector<TString> Labels;

//    ofstream fout("fHalf_SignalEfficiencies.tex");
    ofstream fout(PlotsName[i]+=".tex");
    fout<<"Cuts     & 1000 & 1500 & 2000 & 2500 & 3000 & 3500  \\\\"<<endl<<"\\hline"<<endl;

    for(Int_t k = 1; k <= nbins ; ++k){ 
	ExpEntries_Sig1000.push_back(hSig[0]->GetBinContent(k));
	ExpEntries_Sig1500.push_back(hSig[1]->GetBinContent(k));
	ExpEntries_Sig2000.push_back(hSig[2]->GetBinContent(k));
	ExpEntries_Sig2500.push_back(hSig[3]->GetBinContent(k));
	ExpEntries_Sig3000.push_back(hSig[4]->GetBinContent(k));
	ExpEntries_Sig3500.push_back(hSig[5]->GetBinContent(k));

	Labels.push_back(hSig[0]->GetXaxis()->GetBinLabel(k));

	fout<<Labels[k-1]<<" & "<<ExpEntries_Sig1000[k-1]<<" &  "<<ExpEntries_Sig1500[k-1]<<" &  "<<ExpEntries_Sig2000[k-1]<<" &  "<<ExpEntries_Sig2500[k-1]<<" &  "<<ExpEntries_Sig3000[k-1]<<" &  "<<ExpEntries_Sig3500[k-1]<<" \\"<<"\\"<<endl;

    }
    fout.close();
}
