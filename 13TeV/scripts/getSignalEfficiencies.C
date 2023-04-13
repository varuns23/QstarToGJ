getSignalEfficiencies(){

    const Int_t nSigFiles = 9;
    Int_t NoCut, Final;

    TString InSignalPath = "rootFiles/";
    TString OutPath = "";
    
    TString FileSignals[nSigFiles] = { 
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


    const Int_t nvar = 2;
    TString var[nvar] = {"h_CutFlow", "h_CutExpFlow"};
    TString PlotsName[nvar] = {"CutFlow", "CutExpFlow"};
    TH1F *hSig[nSigFiles];

    // Script for plotting and saving 
    for(Int_t i = 1 ; i < nvar ; ++i){
//    for(Int_t i = 0 ; i < 1 ; ++i){
	for(Int_t j = 0 ; j < nSigFiles ; ++j){
	    TFile* fileIn = new TFile(InSignalPath+FileSignals[j],"READ");
	    hSig[j]  = (TH1F*)fileIn->Get(var[i]);   //Signal
	}
    }

    Int_t nbins = hSig[0]->GetNbinsX();

    std::vector<Double_t> ExpEntries_Sig1000, ExpEntries_Sig2000, ExpEntries_Sig3000, ExpEntries_Sig4000, ExpEntries_Sig5000, ExpEntries_Sig6000, ExpEntries_Sig7000, ExpEntries_Sig8000, ExpEntries_Sig9000;
    std::vector<TString> Labels;

    ofstream fout(OutPath+"SignalEfficiencies.tex");
    fout<<"Cuts     &  1000  & 2000 & 3000 & 4000 & 5000 & 6000 & 7000 & 8000 & 9000 \\\\"<<endl<<"\\hline \\hline"<<endl;

    for(Int_t k = 1; k <= nbins ; ++k){ 
	ExpEntries_Sig1000.push_back(  hSig[0]->GetBinContent(k)  );
	ExpEntries_Sig2000.push_back(  hSig[1]->GetBinContent(k)  );
	ExpEntries_Sig3000.push_back(  hSig[2]->GetBinContent(k)  );
	ExpEntries_Sig4000.push_back(  hSig[3]->GetBinContent(k)  );
	ExpEntries_Sig5000.push_back(  hSig[4]->GetBinContent(k)  );
	ExpEntries_Sig6000.push_back(  hSig[5]->GetBinContent(k)  );
	ExpEntries_Sig7000.push_back(  hSig[6]->GetBinContent(k)  );
	ExpEntries_Sig8000.push_back(  hSig[7]->GetBinContent(k)  );
	ExpEntries_Sig9000.push_back(  hSig[8]->GetBinContent(k)  );

	Labels.push_back(hSig[0]->GetXaxis()->GetBinLabel(k));

//	fout<<Labels[k-1]<<" & "<<ExpEntries_Sig700[k-1]<<" &  "<<ExpEntries_Sig1000[k-1]<<" &  "<<ExpEntries_Sig1200[k-1]<<" &  "<<ExpEntries_Sig1500[k-1]<<" &  "<<ExpEntries_Sig1700[k-1]<<" &  "<<ExpEntries_Sig2000[k-1]<<" &  "<<ExpEntries_Sig2500[k-1]<<" &  "<<ExpEntries_Sig3000[k-1]<<" &  "<<ExpEntries_Sig3500[k-1]<<" &  "<<ExpEntries_Sig4000[k-1]<<" &  "<<ExpEntries_Sig4500[k-1]<<" \\"<<"\\"<<endl;

    fout<<Labels[k-1]<<" & "<<100*ExpEntries_Sig1000[k-1]/ExpEntries_Sig1000[0]<<" &  "<<100*ExpEntries_Sig2000[k-1]/ExpEntries_Sig2000[0]<<" &  "<<100*ExpEntries_Sig3000[k-1]/ExpEntries_Sig3000[0]<<" &  "<<100*ExpEntries_Sig4000[k-1]/ExpEntries_Sig4000[0]<<" &  "<<100*ExpEntries_Sig5000[k-1]/ExpEntries_Sig5000[0]<<" &  "<<100*ExpEntries_Sig6000[k-1]/ExpEntries_Sig6000[0]<<" &  "<<100*ExpEntries_Sig7000[k-1]/ExpEntries_Sig7000[0]<<" &  "<<100*ExpEntries_Sig8000[k-1]/ExpEntries_Sig8000[0]<<" &  "<<100*ExpEntries_Sig9000[k-1]/ExpEntries_Sig9000[0]<<" \\"<<"\\"<<endl;

	if(Labels[k-1] == "Total") NoCut = k-1;
	if(Labels[k-1] == "MassCut") Final = k-1;
    }

   //std::cout<<NoCut<<"  "<<Final<<std::endl;
    fout<<"\\hline"<<std::endl;
    fout<<"\\hline"<<std::endl;
    fout<<"Efficiency" <<" & "<<ExpEntries_Sig1000[Final]/ExpEntries_Sig1000[NoCut]<<" &  "<<ExpEntries_Sig2000[Final]/ExpEntries_Sig2000[NoCut]<<" &  "<<ExpEntries_Sig3000[Final]/ExpEntries_Sig3000[NoCut]<<" &  "<<ExpEntries_Sig4000[Final]/ExpEntries_Sig4000[NoCut]<<" &  "<<ExpEntries_Sig5000[Final]/ExpEntries_Sig5000[NoCut]<<" &  "<<ExpEntries_Sig6000[Final]/ExpEntries_Sig6000[NoCut]<<" &  "<<ExpEntries_Sig7000[Final]/ExpEntries_Sig7000[NoCut]<<" &  "<<ExpEntries_Sig8000[Final]/ExpEntries_Sig8000[NoCut]<<" &  "<<ExpEntries_Sig9000[Final]/ExpEntries_Sig9000[NoCut]<<" \\"<<"\\"<<endl;
 
    fout.close();
}
