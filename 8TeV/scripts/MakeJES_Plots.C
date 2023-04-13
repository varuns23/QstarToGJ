void MakeJES_Plots(){

    const Int_t nFiles = 4;

    TString FileName[nFiles] = {
	"./rootFiles/DiJet_53X8T19pb_Pt170Mass560dEta20TID.root",     // DiJet Bkg
	"./rootFiles/GJ_53X8T19pb_Pt170Mass560dEta20TID.root",        // PhotonJet Bkg Pythia 
	"./rootFiles/EWK_53X8T19pb_Pt170Mass560dEta20TID.root",       // EWK Bkg
	"./rootFiles/Data_53X8T19pb_Pt170Mass560dEta20TID.root"      // Data
    };

    //--- Give here the Data-Signal-Bkg-file names ++++++++++++++++++++++++
    TFile *InFile[nFiles];
    gStyle->SetOptStat(0);

//    for( int i = 0 ; i < nFiles ; ++i){
    for( int i = 1 ; i < 2 ; ++i){

	InFile[i]  = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName[i]);
	if (!InFile[i]) InFile[i]  = new TFile(FileName[i]);
	InFile[i]->cd();

//	h_prof_JES[hi] = new TProfile(name,"Profile for JES vs mass", nMassBins,MassBin);
	TProfile *Prof_JESup = (TProfile*)InFile[i]->Get("h_prof_JESup"); 
	TProfile *Prof_JESdown = (TProfile*)InFile[i]->Get("h_prof_JESdown");

	Prof_JESup->SetMinimum(-1.0);
	Prof_JESup->SetMaximum(1.0);
	Prof_JESup->GetXaxis()->SetRangeUser(300,6800);
	Prof_JESup->Scale(100.0);
	Prof_JESdown->Scale(100.0);

	Prof_JESup->SetFillColor(kSpring-9);
	Prof_JESup->SetLineColor(kBlack);
//	Prof_JESup->SetLineColor(kSpring-9);
	Prof_JESdown->SetFillColor(kSpring-9);
//	Prof_JESdown->SetLineColor(kSpring-9);
	Prof_JESdown->SetLineColor(kBlack);

	
	Prof_JESup->GetXaxis()->SetTitle("M_{#gamma j} (GeV)");
	Prof_JESup->GetXaxis()->CenterTitle();
	Prof_JESup->GetYaxis()->SetTitle("JES Uncertainity (%)");
	Prof_JESup->GetYaxis()->CenterTitle();

	Prof_JESup->Draw("histF");
	Prof_JESdown->Draw("same hist F");

	TLine *drawLineP1 =  new TLine(300, 0.75, 7000, 0.75);
	drawLineP1->SetLineColor(kBlack);
	drawLineP1->SetLineStyle(2);
	drawLineP1->SetLineWidth(1);
	drawLineP1->Draw("same");
	TLine *drawLineM1 =  new TLine(300, -0.75, 7000, -0.75);
	drawLineM1->SetLineColor(kBlack);
	drawLineM1->SetLineStyle(2);
	drawLineM1->SetLineWidth(1);
	drawLineM1->Draw("same");
	TLine *drawLine =  new TLine(300, 0.0, 7000, 0.0);
	drawLine->SetLineColor(kBlack);
	drawLine->SetLineStyle(2);
	drawLine->SetLineWidth(2);
	drawLine->Draw("same");
    }
}
