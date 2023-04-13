void MakeFinalPlots(){


    const Int_t nFiles = 5;
    const Int_t nSigFiles = 9;

    TString FileName[nFiles+1] = {
	"./rootFiles/GJ_53X8T19pb_Pt170Mass560dEta20TID.root",           // PhotonJet Bkg Pythia 
	"./rootFiles/DiJet_53X8T19pb_Pt170Mass560dEta20TID.root",        // DiJet Bkg
	"./rootFiles/EWK_53X8T19pb_Pt170Mass560dEta20TID.root",          // EWK Bkg
	"./rootFiles/Data_53X8T19pb_Pt170Mass560dEta20TID.root",         // Data
	"./rootFiles/QstarToGJ_M_1000.root",       // Qstar signal  1 TeV 
	"./rootFiles/QstarToGJ_M_2500.root"        // Qstar signal  3 TeV
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
    TString Plotslocation  = "../FinalPlots/" ;
    Bool_t  SavePlots      = false ;
  //  Bool_t  SavePlots      = true ;

//    TString var[1]  = {"h_mass_VarBin_finalPU"};
    TString var[1]  = {"h_etaJet_finalPU"};
    TString PlotsName[1] = {"mass_VarBin"};

//    double x_min = 150.0;
//    double x_max = 1500.0;

    double x_min = -3.2;
    double x_max = 3.2;

    //--- Give here the Data-Signal-Bkg-file names ++++++++++++++++++++++++
    TFile* fileBkg1 = new TFile(FileName[0],"READ");
    TFile* fileBkg2 = new TFile(FileName[1],"READ");
    TFile* fileBkg3 = new TFile(FileName[2],"READ");
    TFile* fileData = new TFile(FileName[3],"READ");
    TFile* fileSig1 = new TFile(FileName[4],"READ"); 
    TFile* fileSig2 = new TFile(FileName[5],"READ"); 
    //------------------------------------------------------------

    for( Int_t i = 0 ; i < 1 ; ++i){
    TH1F *hSig1 = (TH1F*)fileSig1->Get(var[i]);   //Signal 1 TeV
    TH1F *hSig2 = (TH1F*)fileSig2->Get(var[i]);   //Signal 3 TeV
    TH1F *hBkg1 = (TH1F*)fileBkg1->Get(var[i]);   //GJ MC
    TH1F *hBkg2 = (TH1F*)fileBkg2->Get(var[i]);   //DiJet MC
    TH1F *hBkg3 = (TH1F*)fileBkg3->Get(var[i]);   //EWK
    TH1F *hData = (TH1F*)fileData->Get(var[i]);   //Data

    hSig1->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    hSig1->SetLineColor(kBlue+1);
    hSig1->SetLineWidth(2);
    hSig1->Rebin(2);
    hSig1->GetXaxis()->SetRangeUser(x_min,x_max);

    hSig2->SetStats(0); // hSig->SetMarkerStyle(3); hSig->SetMarkerColor(kMagenta+4) ;
    hSig2->SetLineColor(kMagenta);
    hSig2->SetLineWidth(2);
    hSig2->Rebin(2);
    hSig2->GetXaxis()->SetRangeUser(x_min,x_max);

    hBkg1->SetStats(0); //hBkg2->SetLineColor(kCyan+2); hBkg2->SetFillColor(kGray) ;
    hBkg1->SetLineColor(kRed); 
    hBkg1->SetLineWidth(2);
    hBkg1->SetLineStyle(1);
    hBkg1->SetFillColor(kGray);
//    hBkg1->Scale(scaleMC*GJkfactor);
    hBkg1->Scale(scaleMC*1.33);
    hBkg1->Rebin(2);
    hBkg1->GetXaxis()->SetRangeUser(x_min,x_max);

    hBkg2->SetStats(0);// hBkg1->SetLineColor(kMagenta-3); hBkg1->SetFillColor(kOrange+8);
    hBkg2->SetLineColor(kGreen+4); 
    hBkg2->SetLineWidth(2);
    hBkg2->SetLineStyle(1);
    hBkg2->SetFillColor(kGreen-3);
    hBkg2->SetFillStyle(3018);
//    hBkg2->Scale(scaleMC);
    hBkg2->Scale(scaleMC*1.34);
    hBkg2->Rebin(2);
    hBkg2->GetXaxis()->SetRangeUser(x_min,x_max);

    hBkg3->SetStats(0); //hBkg2->SetLineColor(kCyan+2); hBkg2->SetFillColor(kGray) ;
    hBkg3->SetLineColor(kBlack); 
    hBkg3->SetLineWidth(1);
    hBkg3->SetLineStyle(1);
    hBkg3->SetFillColor(kSpring-9);
    hBkg3->Scale(scaleMC);
    hBkg3->Rebin(2);
    hBkg3->GetXaxis()->SetRangeUser(x_min,x_max);

    hData->SetStats(0); 
    hData->SetMarkerStyle(8); 
    hData->SetMarkerColor(kBlack); 
    hData->SetMarkerSize(0.8); 
    hData->SetLineColor(kBlack);
    hData->Rebin(2);
    hData->GetXaxis()->SetRangeUser(x_min,x_max);



    // Legends+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    const Int_t nleg = 6 ;
    TString legendName[nleg] = { "q* = 1 TeV", "q* = 2.5 TeV", "GJ (Bkg)", "DiJet (Bkg)", "EWK (Bkg) ", "Data" };
    TString legendFill[nleg] = { "L", "L", "F", "F", "F", "Lp" };
    //---------------------------------------------------------------------

    //    TLegend *leg = new TLegend(0.72,0.76,0.9,0.9);
//    TLegend *leg = new TLegend( 0.7, 0.7, 0.9, 0.9);
    TLegend *leg = new TLegend( 0.63, 0.55, 0.9, 0.88);
//    TLegend *leg = new TLegend( 0.4, 0.05, 0.6, 0.3);
    leg->SetFillColor(0); leg->SetTextFont(42); leg->SetBorderSize(0);
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


    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    //    c1->SetLogy(); c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->Draw(); pad1->cd();
    pad1->SetLogy();
    pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
    pad1->SetBottomMargin(0.);
    THStack* hs = new THStack("hs","");

    hs->SetMinimum(0.01);
    hs->Add(hBkg3); 
    hs->Add(hBkg2);
    hs->Add(hBkg1); 
    hs->SetHistogram(hBkg1);

    TH1F *err_bkg = hBkg1->Clone();
    err_bkg->Add(hBkg2);
    err_bkg->Add(hBkg3);
    err_bkg->SetFillColor(1); err_bkg->SetFillStyle(3004);
    //  gStyle->SetHatchesSpacing(2.0); gStyle->SetHatchesLineWidth(3);err_bkg->SetMarkerStyle(1);

    hs->Draw("hist");
    hSig1->Draw("histSAME");
    hSig2->Draw("histSAME");
    hData->Draw("SAMEP E1");
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
    TH1F *hBkgRatio = hBkg1->Clone();
    hBkgRatio->Add(hBkg2);
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

    hRatio->GetXaxis()->SetLabelSize(0.1);
    hRatio->GetXaxis()->SetTitleSize(0.1);
    hRatio->Draw("e1p");
    
    TLine *ratioLine = new TLine(x_min, 0.0, x_max, 0.0);
//    TLine *ratioLine = new TLine(x_min, 0.0, x_max+20, 0.0);
    ratioLine->SetLineStyle(2);
    ratioLine->SetLineWidth(2);
    ratioLine->SetLineColor(kRed);
    ratioLine->Draw("same");
    hRatio->Draw("psame ");

    c1->cd();

    if(SavePlots){
	 c1->SaveAs(Plotslocation+PlotsName[i]+".pdf");
	 c1->SaveAs(Plotslocation+PlotsName[i]+".eps");
	 c1->SaveAs(Plotslocation+PlotsName[i]+".gif");
	delete c1;
    }


    }

}
