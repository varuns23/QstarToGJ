void fittingInvMass_13TeV_onlyData(){


    // Set TDR Style
    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    
    float t_m = 0.08; //top margin
    float b_m = 0.4; //botton margin
    float l_m = 0.09; //left margin
    float r_m = 0.03; //right margin
    
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(t_m); 
    gStyle->SetPadBottomMargin(b_m);
    gStyle->SetPadLeftMargin(l_m);
    gStyle->SetPadRightMargin(r_m);
    
  
    Double_t minXFit = 600.0;
    Double_t maxXFit = 5000.0;

    // Input Files  for data
//    TFile *inputData = TFile::Open("./rootFiles/Data_ptPhotJet190etaJet2p4dPhi2p5dEta1p8M695LID_76X.root", "READ");
//    TFile *inputData = TFile::Open("./rootFiles/Data2016B_PromptReco_v2_ptPhotJet190etaJet2p4dPhi2p5dEta1p8M695LID_80X.root", "READ");
    TFile *inputData = TFile::Open("./rootFiles/Run2016_24fb.root", "READ");
    TH1F *h_GJMass_data = (TH1F*) inputData->Get("h_GJetInvtMass_VarBin_noBTag_MassCut");

    TH1F *hRatio = (TH1F*)h_GJMass_data->Clone("hRatio");
    hRatio->Reset();

    vector<double> func_Para, func_NDF, func_FCN;

    // Fit to data    
    //    TF1 *fit = new TF1("fit",fitFunc,minXFit,maxXFit,4); 
    TF1 *fit = new TF1("fit",fitFunc,minXFit,maxXFit,4); 
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,6.69761e-04);
    fit->SetParameter(1,4.77861e+00);
    fit->SetParameter(2,6.48175e+00);
    fit->SetParameter(3,4.07809e-01);

    fit->SetLineColor(kRed); //kAzure+7
    fit->SetLineWidth(3);
    fit->SetLineStyle(2);

    for(int a = 0; a < 10; ++a) 
    h_GJMass_data->Fit("fit","R");

    cout << "NDF = " << fit->GetNDF() << " FCN = " << fit->GetChisquare() << endl;

    func_NDF.push_back(fit->GetNDF());
    func_FCN.push_back(fit->GetChisquare());

    func_Para.push_back(fit->GetParameter(0));
    func_Para.push_back(fit->GetParameter(1));
    func_Para.push_back(fit->GetParameter(2));
    func_Para.push_back(fit->GetParameter(3));

    gStyle->SetOptFit(0000);

    ///////----------------------
    // Now making the data-fit ratio plot
    const double alpha = 0.3173;  // this following the CMS statistic commitee (https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars)
    int n, ibin; //float
    float dm, mass, xl, xh;
    float vx[1200],vy[1200],vexl[1200],vexh[1200],veyl[1200],veyh[1200];
    Double_t dataa,fitt, ratioo, errorr;

    for(ibin=0; ibin < h_GJMass_data->GetNbinsX(); ++ibin){
	n    = h_GJMass_data->GetBinContent(ibin+1);
	dm   = h_GJMass_data->GetBinWidth(ibin+1);
	mass = h_GJMass_data->GetBinCenter(ibin+1);
	vx[ibin]   = mass;
	vy[ibin]   = n;

	vexl[ibin] = dm/2.;
	vexh[ibin] = dm/2.;

	double LOW = (n == 0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,n,1.)); 
	double UP = ROOT::Math::gamma_quantile_c(alpha/2,n+1,1) ;

	veyl[ibin] = n - LOW;
	veyh[ibin] = UP - n;

	if(n !=0){
	    dataa = vy[ibin];
	    errorr= veyl[ibin];
	    fitt   = fit->Eval(vx[ibin],0,0);
	    //	    ratioo = (dataa-fitt)/dataa;
	    ratioo = (dataa-fitt)/errorr;
	    hRatio->SetBinContent(ibin+1,ratioo);

	}


    }

    // Graphs to Make Signal points ----------
    TGraphAsymmErrors *g = new TGraphAsymmErrors(ibin,vx,vy,vexl,vexh,veyl,veyh);

    //Gamjet Mass Cross Section with Fit	
    TCanvas* c1 = new TCanvas("c1","GamjetMass Cross Section with Fit", 900, 700);
    
    TPad * p1_1 = new TPad("p1_1","",0.0,0.24,1,1);
    p1_1->SetLogy();
    p1_1->SetBottomMargin(0);
    p1_1->Draw();
    p1_1->cd();

    h_GJMass_data->SetTitle("");
    h_GJMass_data->SetStats(0);
    h_GJMass_data->SetLineColor(1);
    h_GJMass_data->SetFillColor(1);
    h_GJMass_data->SetMarkerColor(1);
    h_GJMass_data->SetMarkerStyle(20);
    h_GJMass_data->GetYaxis()->SetTitle("Events per bin");
    h_GJMass_data->GetYaxis()->CenterTitle();
    h_GJMass_data->GetYaxis()->SetLabelSize(0.045);
    h_GJMass_data->GetYaxis()->SetLabelOffset(0.001);
    h_GJMass_data->GetYaxis()->SetTitleSize(0.06);
    h_GJMass_data->GetYaxis()->SetTitleOffset(0.7);    
    h_GJMass_data->GetYaxis()->SetRangeUser(0.3,1e5);
    h_GJMass_data->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    h_GJMass_data->GetXaxis()->SetTickLength(0.02);

    g->SetMarkerColor(1);
    g->SetMarkerSize(0.01);
    g->SetLineColor(1);


    h_GJMass_data->Draw("HIST P");
    //h_GJMass_jp->Draw("sameHist");
    g->Draw("samePE1");
    fit->Draw("same");
    //    h_err_bkg->Draw("same E2");
    //    h_GJMass_data->Draw("PZ E1same");

    //TLegend *leg = new TLegend(0.35,0.64,0.6,0.84);
    TLegend *leg = new TLegend(0.15,0.1,0.4,0.3);
    leg->SetTextFont(42);//
    leg->SetTextSize(0.048);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(h_GJMass_data,"Data","PLE"); 
    leg->AddEntry(fit,"Background-only fit to data","L");
    leg->Draw("same");

    float lumiTextSize = 0.6;
    float lumiTextOffset = 0.2;
    TLatex lumi;
    lumi.SetNDC();
    lumi.SetTextAngle(0);
    lumi.SetTextColor(kBlack);
    lumi.SetTextFont(42);
    lumi.SetTextAlign(31);
    lumi.SetTextSize(lumiTextSize*t_m);
    lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "24.45 fb^{-1} (13 TeV)");
     
    float cmsTextFont = 61;
    float cmsTextSize = 0.75;
    //    float posX_ = l_m     + 0.045 * (1 - l_m - r_m);  // Top left
    //    float posX_ = l_m     + 0.5   * (1 - l_m - r_m);  // Centre
    float posX_ = 1 - r_m + 0.045 * (1 - l_m - r_m);  // Top right
    float posY_ = 1 - t_m - 0.035 * (1 - t_m - b_m);
     
    TLatex cms; 
    cms.SetNDC();
    cms.SetTextFont(cmsTextFont);
    cms.SetTextSize(cmsTextSize * t_m);
    cms.SetTextAlign(33);  // 11-top left;  21-top centre;  31-top right
//    cms.DrawLatex(posX_, posY_, "CMS");
    cms.DrawLatex(0.9, 0.87, "CMS");
    cms.SetTextFont(52);
    cms.SetTextAlign(33);
    cms.SetTextSize(0.76*cmsTextSize*t_m);
    //cms.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t_m, extraText);
    cms.DrawLatex(0.92, 0.88 - 1.2*cmsTextSize*t_m, "Preliminary");   

    TLatex *fit_latex = new TLatex();
    fit_latex->SetTextAlign(13);
    fit_latex->SetTextFont(42);
    fit_latex->SetNDC();
    fit_latex->SetTextSize(0.055);
    fit_latex->DrawLatex(0.82,0.73,"q*#rightarrow q#gamma");

    p1_1->RedrawAxis();
    c1->cd();

    //---- Next PAD
    TPad *p1_2 = new TPad("p1_2","ratio",0.0,0.0,1,0.24);
    p1_2->SetTopMargin(0);
    p1_2->SetGridx();
    p1_2->SetGridy();
    p1_2->Draw();
    p1_2->cd();
    //c1_2->SetTickx(50);
  
    hRatio->SetXTitle("M_{#gamma jet} [GeV]");
    hRatio->SetStats(0);
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitleSize(0.18);
    hRatio->GetXaxis()->SetTitleOffset(1);
    hRatio->GetXaxis()->SetLabelSize(0.15); //0.14
    //hRatio->GetXaxis()->SetNdivisions(520);
    hRatio->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    
    hRatio->SetYTitle("#frac{(Data-Fit)}{#sigma_{Data}}");
    hRatio->GetYaxis()->CenterTitle();
    hRatio->GetYaxis()->SetRangeUser(-2.9,2.9);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.13);
    hRatio->GetYaxis()->SetTitleOffset(0.25);
    hRatio->GetYaxis()->SetLabelSize(0.12);
    hRatio->GetYaxis()->SetLabelOffset(0.007);

//    hRatio->SetFillColor(kYellow);
//    hRatio->SetLineColor(kGreen+2);
    hRatio->SetFillColor(kYellow);
    hRatio->SetLineColor(kRed);
    hRatio->SetLineWidth(2);
    hRatio->Draw("SAMEHIST");

//    c1->SaveAs("DataMC_Fit_PAS.pdf");
/*    c1->SaveAs("FinalFit_DataFitMC_pas.eps");
    
    TImage *img = TImage::Create();
    img->FromPad(c1);
    img->WriteImage("FinalFit_DataFitMC_pas.png");
    img->WriteImage("FinalFit_DataFitMC_pas.gif");
  */  

}

Double_t fitFunc( Double_t *m, Double_t *p) // function 1
{
    double x=m[0]/13000.;
    return p[0]*pow(1.0-x,p[1])/pow(x,p[2]+p[3]*log(x));
}
