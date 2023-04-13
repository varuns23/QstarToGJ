DrawSigMassPoints(){

    const Int_t nSigFiles = 11;
    const Int_t nHalfSigFiles = 9;

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
    TString legSignals[nSigFiles] = { 
	"0.7 TeV ",
	"1.0 TeV",
	"1.2 TeV",
	"1.5 TeV",
	"1.7 TeV",
	"2.0 TeV",
	"2.5 TeV",
	"3.0 TeV",
	"3.5 TeV",
	"4.0 TeV",
	"4.5 TeV"
    };



    TString FileHalfSignals[nHalfSigFiles] = { 
	"./rootFiles/QstarToGJ_M_fhalf_700.root",
	"./rootFiles/QstarToGJ_M_fhalf_1000.root",
	"./rootFiles/QstarToGJ_M_fhalf_1500.root",
	"./rootFiles/QstarToGJ_M_fhalf_2000.root",
	"./rootFiles/QstarToGJ_M_fhalf_2500.root",
	"./rootFiles/QstarToGJ_M_fhalf_3000.root",
	"./rootFiles/QstarToGJ_M_fhalf_3500.root",
	"./rootFiles/QstarToGJ_M_fhalf_4000.root",
	"./rootFiles/QstarToGJ_M_fhalf_4500.root"
    };
    TString legHalfSignals[nHalfSigFiles] = { 
	"0.7 TeV",
	"1.0 TeV",
	"1.5 TeV",
	"2.0 TeV",
	"2.5 TeV",
	"3.0 TeV",
	"3.5 TeV",
	"4.0 TeV",
	"4.5 TeV"
    };


    TCanvas* c2 = new TCanvas("c2","Qstar f = 1",80,20,1000,700);
    c2->SetRightMargin(0.13);
    c2->SetLeftMargin(0.08);
    gStyle->SetOptTitle(0);

    TLegend *leg = new TLegend(0.88,0.1,1.,0.9);
    leg->SetFillColor(0); leg->SetTextFont(42); leg->SetBorderSize(0);
    leg->SetTextAlign(22);

    TH1F *h_Sighist[nSigFiles];
    for(Int_t j = 0 ; j < nSigFiles ; ++j){

	TFile* fileIn = new TFile(FileSignals[j],"READ");
	h_Sighist[j] = (TH1F*)fileIn->Get("h_mass_VarBin_finalPU");

	h_Sighist[j]->SetStats(0);
	if(j==9) h_Sighist[j]->SetLineColor(j+3); 
	else	h_Sighist[j]->SetLineColor(j+1);
	h_Sighist[j]->SetLineStyle(j);
	h_Sighist[j]->SetLineWidth(2);
	h_Sighist[j]->GetXaxis()->SetRangeUser(400.0,5000.0);
	h_Sighist[j]->Scale(1.0/h_Sighist[j]->Integral());

	if(j == 0) h_Sighist[j]->Draw("hist");
	else h_Sighist[j]->Draw("sameHist");
	
	leg->AddEntry(h_Sighist[j],legSignals[j],"F"); 
	leg->Draw();
    }

//    c2->SaveAs("QstarSignalShapes.pdf");
//    delete c2;
    
    TCanvas* c1 = new TCanvas("c1","Qstar f = 0.5",80,20,1000,700);
    c1->SetRightMargin(0.13);
    c1->SetLeftMargin(0.08);
    gStyle->SetOptTitle(0);

    TLegend *leg1 = new TLegend(0.88,0.1,1.,0.9);
    leg1->SetFillColor(0); leg1->SetTextFont(42); leg1->SetBorderSize(0);
    leg1->SetTextAlign(22);

    TH1F *h_HalfSighist[nHalfSigFiles];
    for(Int_t j = 0 ; j < nHalfSigFiles ; ++j){

	TFile* fileIn = new TFile(FileHalfSignals[j],"READ");
	h_HalfSighist[j] = (TH1F*)fileIn->Get("h_mass_VarBin_finalPU");

	h_HalfSighist[j]->SetStats(0);
	h_HalfSighist[j]->SetLineColor(j+1);
	h_HalfSighist[j]->SetLineStyle(j);
	h_HalfSighist[j]->SetLineWidth(2);
	h_HalfSighist[j]->GetXaxis()->SetRangeUser(400.0,5000.0);
	h_HalfSighist[j]->Scale(1.0/h_HalfSighist[j]->Integral());

	if(j == 0) h_HalfSighist[j]->Draw("hist");
	else h_HalfSighist[j]->Draw("sameHist");
	
	leg1->AddEntry(h_HalfSighist[j],legHalfSignals[j],"F"); 
	leg1->Draw();
    }

    //c1->SaveAs("QstarSignalShapes_fhalf.pdf");
}
