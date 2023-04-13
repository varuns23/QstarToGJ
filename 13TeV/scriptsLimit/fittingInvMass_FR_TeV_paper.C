#include "vector.h"


void fittingInvMass_FR_TeV_paper(){


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

//    float t_m = 0.08; //top margin
//    float b_m = 0.12; //botton margin
//    float l_m = 0.12; //left margin
//    float r_m = 0.04; //right margin



    // Input Files  for data
    TFile *inputData = TFile::Open("./rootFiles/DataBkgSig_Mass_TeVScale.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_data");

    TH1F *hRatio = (TH1F*)hGJMass_data->Clone("hRatio");
    hRatio->Reset();

    // Input Files  for MC
    TFile *inputMC = TFile::Open("./rootFiles/DataBkgSig_Mass_TeVScale.root", "READ");
    TH1F *hGJMass_mc = (TH1F*)inputMC->Get("h_mass_MC");
    hGJMass_mc->SetLineWidth(2);

    TFile *inputSignal1 = TFile::Open("./rootFiles/DataBkgSig_Mass_TeVScale.root", "READ");
    TH1F *h_qstar_1000 = (TH1F*)inputSignal1->Get("h_mass_Sig1000");

    TFile *inputSignal2 = TFile::Open("./rootFiles/DataBkgSig_Mass_TeVScale.root", "READ");
    TH1F *h_qstar_1500 = (TH1F*)inputSignal2->Get("h_mass_Sig1500");

    TFile *inputSignal3 = TFile::Open("./rootFiles/DataBkgSig_Mass_TeVScale.root", "READ");
    TH1F *h_qstar_f2500 = (TH1F*)inputSignal3->Get("h_mass_Sig2500f");

    Double_t minXFit = 0.56 ;
    Double_t maxXFit = 3.1;

    vector<double> func_Para, func_NDF, func_FCN;

    // Fit to data    

    TF1 * fitFunc = new TF1("fitFunc","((2.74277e-03*pow(1.0-((x)/8000.0),3.33455))/(pow(((x)/8000.0), (8.06973 + 7.47310e-01*log((x)/8000.0)))))",560.0,3100.0);

    ///////----------------------
    // Now making the data-fit ratio plot
    const double alpha = 0.3173;  // this following the CMS statistic commitee (https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars)
    int n; //float
    float dm, mass, xl, xh;
    float vx[1200],vy[1200],vexl[1200],vexh[1200],veyl[1200],veyh[1200], vfx[1200], vfy[1200];
    int i,j;
    Double_t dataa,fitt, ratioo, errorr;
    std::vector<double> qstar1000_Y, qstar1000_X, qstar1000_Err;
    std::vector<double> qstar1500_Y, qstar1500_X, qstar1500_Err;
    std::vector<double> qstarf2500_Y, qstarf2500_X, qstarf2500_Err;


    //Draw fit function using graph
    for( j = 32; j <=73 ; ++j){
	vfx[j] = hGJMass_data->GetBinCenter(j);
	vfy[j] = fitFunc->Eval(vfx[j]*1000,0,0);
    }
   TGraph *grfit = new TGraph(j,vfx,vfy);
   grfit->SetLineColor(3);
   grfit->SetLineStyle(2);
   grfit->SetLineWidth(3);

    //Draw data with correct error bars, residuals
    for(i=0;i<=72;i++){

	n    = hGJMass_data->GetBinContent(i+1);
	dm   = hGJMass_data->GetBinWidth(i+1);
	mass = hGJMass_data->GetBinCenter(i+1);
	xl   = hGJMass_data->GetBinLowEdge(i+1);
	xh   = xl+dm; 
	vx[i]   = (xl+xh)/2.;
	vy[i]   = n;

	vexl[i] = dm/2.;
	vexh[i] = dm/2.;

	double LOW = (n == 0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,n,1.)); 
	double UP = ROOT::Math::gamma_quantile_c(alpha/2,n+1,1) ;

	veyl[i] = n - LOW;
	veyh[i] = UP - n;

	if(i == 68 || i == 70 || i == 72){
	    veyl[i] = 0.0;
	    veyh[i] = 1.8;

	    dataa = vy[i];
	    errorr= veyh[i];
	    fitt   = fitFunc->Eval(vx[i]*1000,0,0);
	    ratioo = (dataa-fitt)/errorr;
	    hRatio->SetBinContent(i+1,ratioo);
	}

	if(i > 72){
	    veyl[i] = 0.0;
	    veyh[i] = 0.0;
	}


	if(n !=0){
	    dataa = vy[i];
	    errorr= veyl[i];
	    fitt   = fitFunc->Eval(vx[i]*1000,0,0);
	    ratioo = (dataa-fitt)/errorr;
	    hRatio->SetBinContent(i+1,ratioo);
	}

	if(mass > 0.8*1.0 && mass < 1.3*1.0){
	    qstar1000_Y.push_back(h_qstar_1000->GetBinContent(i+1) + (hGJMass_mc->GetBinContent(i+1)*0.62));
	    qstar1000_X.push_back(mass);
	}

	if(mass > 0.85*1.5 && mass < 1.2*1.5){
	    qstar1500_Y.push_back(h_qstar_1500->GetBinContent(i+1) + (hGJMass_mc->GetBinContent(i+1)*0.62));
	    qstar1500_X.push_back(mass);
	}

	if(mass > 0.9*2.5 && mass < 1.12*2.5){
	    qstarf2500_Y.push_back(h_qstar_f2500->GetBinContent(i+1) + (hGJMass_mc->GetBinContent(i+1))); //0.9
	    qstarf2500_X.push_back(mass);
	}

    }

    // Graphs to Make Signal points ----------
    TGraphAsymmErrors *g = new TGraphAsymmErrors(i,vx,vy,vexl,vexh,veyl,veyh);

    TGraph* gr_qstar1000 = new TGraph(qstar1000_X.size() , &qstar1000_X[0], &qstar1000_Y[0]);
    gr_qstar1000->SetLineColor(kOrange+2);
    gr_qstar1000->SetLineStyle(9); //7
    gr_qstar1000->SetLineWidth(2);
    gr_qstar1000->SetTitle();


    TGraph* gr_qstar1500 = new TGraph(qstar1500_X.size() , &qstar1500_X[0], &qstar1500_Y[0]); 
    gr_qstar1500->SetLineColor(kCyan+2);
    gr_qstar1500->SetLineStyle(9); //7
    gr_qstar1500->SetLineWidth(2);
    gr_qstar1500->SetTitle();

    TGraph* gr_qstarf2500 = new TGraph(qstarf2500_X.size() , &qstarf2500_X[0], &qstarf2500_Y[0]); 
    gr_qstarf2500->SetLineColor(kMagenta+2);
    gr_qstarf2500->SetLineStyle(9); //7
    gr_qstarf2500->SetLineWidth(2);
    gr_qstarf2500->SetTitle();


    //Dijet Mass Cross Section with Fit	
    TCanvas* c1 = new TCanvas("c1","GammajetMass Cross Section with Fit",900,700);

    TPad *p1_1 = new TPad("p1_1","Photon + jet mass",0.0,0.24,1,1);
    p1_1->SetLogy();
    p1_1->SetBottomMargin(0);
    p1_1->Draw();
    p1_1->cd();

    hGJMass_data->SetTitle("");
    hGJMass_data->SetLineColor(1);
    hGJMass_data->SetFillColor(1);
    hGJMass_data->SetMarkerColor(1);
    hGJMass_data->SetMarkerStyle(20);

    hGJMass_mc->SetTitle("");
    hGJMass_mc->GetYaxis()->SetTitle("Events per bin");
    hGJMass_mc->GetYaxis()->CenterTitle();
    hGJMass_mc->GetYaxis()->SetLabelSize(0.045);
    hGJMass_mc->GetYaxis()->SetLabelOffset(0.001);
    hGJMass_mc->GetYaxis()->SetTitleSize(0.06);
    hGJMass_mc->GetYaxis()->SetTitleOffset(0.7);
    hGJMass_mc->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hGJMass_mc->GetXaxis()->SetTickLength(0.02);

    g->SetMarkerColor(1);
    g->SetLineColor(1);

    hGJMass_mc->Draw("Hist");
    g->Draw("same PE");
    grfit->Draw("same");
    gr_qstar1000->Draw("sameC");
    gr_qstar1500->Draw("sameC");
    gr_qstarf2500->Draw("sameC");

    TLatex *SignalName = new TLatex();
    SignalName->SetTextFont(62);
    SignalName->SetNDC();
    SignalName->SetTextSize(0.04);
    SignalName->SetTextColor(kOrange+2);
    SignalName->DrawLatex(0.19,0.87,"q* (1.0 TeV, f = 1.0)");
    SignalName->SetTextColor(kCyan+2);
    SignalName->DrawLatex(0.35,0.69,"q* (1.5 TeV, f = 1.0)");
    SignalName->SetTextColor(kMagenta+2);
    SignalName->DrawLatex(0.68,0.343,"q* (2.5 TeV, f = 0.5)");

//    TLegend *leg = new TLegend(0.58,.6,0.83,0.82);
    TLegend *leg = new TLegend(0.16,.09,0.41,0.31);
    leg->SetTextFont(42);//42
    leg->SetTextSize(0.048);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hGJMass_data,"Data","PLE"); 
    leg->AddEntry(hGJMass_mc,"Simulated background","L");
    leg->AddEntry(grfit,"Background-only fit to data","L");
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
    lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "19.7 fb^{-1} (8 TeV)");

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
    
    TLatex *fit_latex = new TLatex();
    fit_latex->SetTextAlign(13);
    fit_latex->SetTextFont(42);
    fit_latex->SetNDC();
    fit_latex->SetTextSize(0.055);
    fit_latex->DrawLatex(0.82,0.76,"q*#rightarrow q#gamma");


    p1_1->RedrawAxis();
    c1->cd();
    //p11_1->RedrawAxis();
    //p11_1->Update();

    //---- Next PAD
    TPad *p1_2 = new TPad("p1_2","ratio",0.0,0.0,1,0.24);
    p1_2->SetTopMargin(0);
    p1_2->SetGridx();
    p1_2->SetGridy();
    p1_2->Draw();
    p1_2->cd();
     
    hRatio->SetXTitle("M_{#gamma,jet} [TeV]");
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitleSize(0.18);  //18
    hRatio->GetXaxis()->SetTitleOffset(1); //1.1
    hRatio->GetXaxis()->SetLabelSize(0.15); //0.14
    hRatio->GetXaxis()->SetNdivisions(520);
    hRatio->GetXaxis()->SetRangeUser(minXFit,maxXFit);

    hRatio->SetYTitle("#frac{(Data-Fit)}{#sigma_{Data}}");
//    hRatio->SetYTitle("(Data-Fit)/#sigma_{Data}        ");
    hRatio->GetYaxis()->SetRangeUser(-3.5,3.5);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.13);
    hRatio->GetYaxis()->SetTitleOffset(0.25);
    hRatio->GetYaxis()->SetLabelSize(0.09);
    hRatio->GetYaxis()->SetLabelOffset(0.007);

    hRatio->SetLineWidth(0);
    hRatio->SetFillColor(2);
    hRatio->SetLineColor(1);
    hRatio->Draw("HIST");
    hRatio->Draw("sameHIST");

       c1->SaveAs("DataFitMC_TeV.pdf");
       
       //       c1->SaveAs("DataFitMC_FR_TeV.eps");

       TImage *img = TImage::Create();
       img->FromPad(c1);
       img->WriteImage("DataFitMC_TeV.png");
//       img->WriteImage("DataFitMC_FR_TeV.gif");
}

Double_t fitFunc1( Double_t *m, Double_t *p) // function 1
{
    double x=m[0]/8000.;
    return p[0]*pow(1.0-x,p[1])/pow(x,p[2]+p[3]*log(x));
}
