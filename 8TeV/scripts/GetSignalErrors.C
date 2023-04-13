GetSignalErrors(){



    const Int_t nSigFiles = 11;
    TString FileSignals[nSigFiles] = { 
	"./rootFiles/QstarToGJ_M_700.root",
	"./rootFiles/QstarToGJ_M_1000.root",
	"./rootFiles/QstarToGJ_M_1200.root",
	"./rootFiles/QstarToGJ_M_1500.root",
	"./rootFiles/QstarToGJ_M_1700.root",
	"./rootFiles/QstarToGJ_M_2000.root",
	"./rootFiles/QstarToGJ_M_2500.root",
	"./rootFiles/QstarToGJ_M_3000.root",
	"./rootFiles/QstarToGJ_M_3500.root",
	"./rootFiles/QstarToGJ_M_4000.root",
	"./rootFiles/QstarToGJ_M_4500.root"
    };

    TFile *InSigFile[nSigFiles];
    TH1F  *h_CutFlow[nSigFiles];

    Double_t InitialEvents[nSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Double_t InitialError[nSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Double_t masses[nSigFiles] = { 700, 1000, 1200, 1500, 1700, 2000, 2500, 3000, 3500, 4000, 4500 };

    Double_t FinalEvents[nSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Double_t FinalError[nSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    Double_t TotalError[nSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for( int i = 0 ; i < nSigFiles ; ++i){

	InSigFile[i]  = (TFile*)gROOT->GetListOfFiles()->FindObject(FileSignals[i]);
	if (!InSigFile[i]) InSigFile[i]  = new TFile(FileSignals[i]);
	InSigFile[i]->cd();

	h_CutFlow[i] = (TH1F*)InSigFile[i]->Get("h_CutFlowTable");

	InitialEvents[i] = h_CutFlow[i]->GetBinContent(1);
	FinalEvents[i] = h_CutFlow[i]->GetBinContent(13);

	TotalError[i] = sqrt(
		pow( (100*(FinalEvents[i]/InitialEvents[i])*sqrt( (1.0/InitialEvents[i]) + (1.0/FinalEvents[i]) )), 2) +     // statistical error a/b
		2.437*2.437 +       // lictd scale factor error
		2.337*2.337      // photonid scale factor error
		); 

	std::cout<<masses[i]<<"   "<<TotalError[i]<<std::endl;
    }

    cout<<endl;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const Int_t nhalfSigFiles = 8;
    TString FilehalfSignals[nhalfSigFiles] = { 
	"./rootFiles/QstarToGJ_M_fhalf_1000.root",
	"./rootFiles/QstarToGJ_M_fhalf_1500.root",
	"./rootFiles/QstarToGJ_M_fhalf_2000.root",
	"./rootFiles/QstarToGJ_M_fhalf_2500.root",
	"./rootFiles/QstarToGJ_M_fhalf_3000.root",
	"./rootFiles/QstarToGJ_M_fhalf_3500.root",
	"./rootFiles/QstarToGJ_M_fhalf_4000.root",
	"./rootFiles/QstarToGJ_M_fhalf_4500.root"
    }

    TFile *InhalfSigFile[nhalfSigFiles];
    TH1F  *h_halfCutFlow[nhalfSigFiles];

    Double_t InitialEventshalf[nhalfSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Double_t InitialErrorhalf[nhalfSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    TString masseshalf[nhalfSigFiles] = { "fhalf 1000", "fhalf 1500", "fhalf 2000", "fhalf 2500", "fhalf 3000", "fhalf 3500", "fhalf 4000", "fhalf 4500" };

    Double_t FinalEventshalf[nhalfSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    Double_t FinalErrorhalf[nhalfSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    Double_t TotalErrorhalf[nhalfSigFiles] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    for( int i = 0 ; i < nhalfSigFiles ; ++i){

	InhalfSigFile[i]  = (TFile*)gROOT->GetListOfFiles()->FindObject(FilehalfSignals[i]);
	if (!InhalfSigFile[i]) InhalfSigFile[i]  = new TFile(FilehalfSignals[i]);
	InhalfSigFile[i]->cd();

	h_halfCutFlow[i] = (TH1F*)InhalfSigFile[i]->Get("h_CutFlowTable");

	InitialEventshalf[i] = h_halfCutFlow[i]->GetBinContent(1);
	FinalEventshalf[i] = h_halfCutFlow[i]->GetBinContent(13);

	TotalErrorhalf[i] = sqrt(
		pow( (100*(FinalEventshalf[i]/InitialEventshalf[i])*sqrt( (1.0/InitialEventshalf[i]) + (1.0/FinalEventshalf[i]) )), 2) +     // statistical error a/b
		2.437*2.437 +       // lictd scale factor error
		2.337*2.337      // photonid scale factor error
		); 

	std::cout<<masseshalf[i]<<"   "<<TotalErrorhalf[i]<<std::endl;
    }


}

