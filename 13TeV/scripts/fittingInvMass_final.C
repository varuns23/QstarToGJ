#include "vector.h"


void fittingInvMass_final(){


    // Set TDR Style	
    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // Input Files  for data
    TFile *inputData = TFile::Open("./rootFiles/Run2015D_PromptReco.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_VarBin_MassCut");

    TH1F *hRatio = (TH1F*)hGJMass_data->Clone("hRatio");
    hRatio->Reset();

    // Input Files  for MC
    TFile *inputMC1 = TFile::Open("./rootFiles/GJets_allHTBins.root", "READ");
    TH1F *hGJMass_mc1 = (TH1F*)inputMC1->Get("h_mass_VarBin_MassCut");
    TFile *inputMC2 = TFile::Open("./rootFiles/QCD_allPtBins.root", "READ");
    TH1F *hGJMass_mc2 = (TH1F*)inputMC2->Get("h_mass_VarBin_MassCut");

    TFile *inputSignal1 = TFile::Open("./rootFiles/QstarToGJ_M2000_f0p1_1.root", "READ");
    TH1F *h_qstar_2000 = (TH1F*)inputSignal1->Get("h_mass_VarBin_MassCut");

    TH1F *hGJMass_mc= (TH1F*)hGJMass_mc1->Clone("*hGJMass_mc") ;

    hGJMass_mc->Scale(1.0);
    hGJMass_mc->Add(hGJMass_mc2,1.0);
    hGJMass_mc->SetLineWidth(2);

    Double_t minXFit = 560.0 ;
    Double_t maxXFit = 3500.0;

    vector<double> func_Para, func_NDF, func_FCN;

    // Fit to data    
    //    TF1 *fit = new TF1("fit",fitFunc1,minXFit,maxXFit,4); 
    TF1 *fit = new TF1("fit",fitFunc1,minXFit,maxXFit,4); 
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,2.74281e-03);
    fit->SetParameter(1,3.33458e+00);
    fit->SetParameter(2,8.06973e+00);
    fit->SetParameter(3,7.47309e-01);

    /*    TF1 *fit = new TF1("fit",fitFunc2,minXFit,maxXFit,4); 
	  gStyle->SetOptFit(0000000); 
	  fit->SetParameter(0,3.49602e-02);
	  fit->SetParameter(1,-9.28015e+00);
	  fit->SetParameter(2,6.79711e+00);
	  fit->SetParameter(3,5.74076e-01);
     */

    //    fit->SetParameter(0,4.90863e-07);
    //    fit->SetParameter(1,9.99858e+00);
    //    fit->SetParameter(2,1.31627e+01);
    //    fit->SetParameter(3,1.57651e+00);

    fit->SetLineColor(3);
    fit->SetLineWidth(2);

    for(int a = 0; a < 10; ++a) hGJMass_data->Fit("fit","R");
    //    std::cout<<" Data : "<<hGJMass_data->Integral()<<"   Fit Integral : "<<fit->Integral(560,1000)<<std::endl;

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
    int n; //float
    float dm, mass, xl, xh;
    float vx[1200],vy[1200],vexl[1200],vexh[1200],veyl[1200],veyh[1200];
    int i;
    Double_t dataa,fitt, ratioo, errorr;

    for(i=0;i<hGJMass_data->GetNbinsX();i++){
    //for(i=0;i<72;i++){
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

	if(n !=0){
	    dataa = vy[i];
	    errorr= veyl[i];
	    fitt   = fit->Eval(vx[i],0,0);
	    //	    ratioo = (dataa-fitt)/dataa;
	    ratioo = (dataa-fitt)/errorr;
	    hRatio->SetBinContent(i+1,ratioo);
	}

    }

    // Graphs to Make Signal points ----------
    TGraphAsymmErrors *g = new TGraphAsymmErrors(i,vx,vy,vexl,vexh,veyl,veyh);


    //Gamjet Mass Cross Section with Fit	
    TCanvas* c1 = new TCanvas("c1","GamjetMass Cross Section with Fit");
    //    gStyle->SetOptFile(0);
    //    gStyle->SetOptStat(0);
    //    gStyle->SetPadBorderMode(0);
    c1->Divide(1,2,0,0,0);
    c1->cd(1);
    p11_1 = (TPad*)c1->GetPad(1);
    p11_1->SetPad(0.01,0.24,0.99,0.98);
    p11_1->SetLogy();
    p11_1->SetRightMargin(0.05);
    p11_1->SetTopMargin(0.06);

    TH1F *vFrame = p11_1->DrawFrame(minXFit,0.000000001,maxXFit,10.0);
    vFrame->SetTitle("");
    //    vFrame->SetXTitle("Mass M_{#gamma j} [GeV]");
    vFrame->GetXaxis()->SetTitleSize(0.06);
    vFrame->SetYTitle("");

    hGJMass_data->SetTitle("");
    hGJMass_data->SetLineColor(1);
    hGJMass_data->SetFillColor(1);
    hGJMass_data->SetMarkerColor(1);
    hGJMass_data->SetMarkerStyle(20);
    hGJMass_data->GetYaxis()->SetTitle("Events");
    hGJMass_data->GetYaxis()->SetLabelSize(0.045);
    hGJMass_data->GetYaxis()->SetTitleSize(0.05);
    hGJMass_data->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hGJMass_data->GetXaxis()->SetTickLength(0.02);

    hGJMass_mc->SetTitle("");
    hGJMass_mc->GetYaxis()->SetTitle("Events");
    hGJMass_mc->GetYaxis()->SetLabelSize(0.045);
    hGJMass_mc->GetYaxis()->SetTitleSize(0.05);
    hGJMass_mc->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hGJMass_mc->GetXaxis()->SetTickLength(0.02);

    h_qstar_2000->SetLineColor(kOrange+2);
    h_qstar_2000->SetLineStyle(7);
    h_qstar_2000->SetLineWidth(2);
    h_qstar_2000->SetTitle();

    g->SetMarkerColor(1);
    g->SetLineColor(1);


    //+++++++++++ for Error ++++++++++++++++++++++++++++++++++
    TH1F *h_err_bkg = hGJMass_mc->Clone("h_err_bkg");

    h_err_bkg->SetFillColor(kRed-2); 
    h_err_bkg->SetFillStyle(3244); 

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //    h_err_bkg->Draw("same E2");
    hGJMass_mc->Draw("Hist");
    g->Draw("samePE1");
    fit->Draw("same");
    //    hGJMass_data->Draw("PZ E1same");
    //   gr_qstar2000->Draw("sameC");
//    h_qstar_2000->Draw("same");
    //    hRatio->Draw("hist");

    //    TLegend *leg = new TLegend(0.2,0.1,0.5,0.4);
    TLegend *leg = new TLegend(0.62,.65,0.9,0.9);
    //    leg->SetTextSize(0.03146853);
    leg->SetTextFont(62);//42
    leg->SetTextSize(0.04);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hGJMass_data,"Data","PL"); 
    //  leg->AddEntry(hGJMass_data,Form("CMS  (%.3f fb^{-1})", 19.6),"PL"); 
    leg->AddEntry(fit,"Background from Fit","L");
    //   leg->AddEntry(h_tmp,"1#sigma Fit Error ","F");
    leg->AddEntry(hGJMass_mc,"Background from MC","L");
    //    leg->AddEntry(h_err_bkg,"Background Uncertainty ","F");
    leg->Draw("same");


    //    TPaveText *fit_para = new TPaveText(0.2,0.045,0.474,0.4,"NDC");
    TLatex *fit_latex = new TLatex();
    fit_latex->SetTextAlign(13);
    fit_latex->SetTextFont(62);
    fit_latex->SetNDC();
    fit_latex->SetTextSize(0.04);
    fit_latex->DrawLatex(0.202,0.36,"CMS Preliminary");
//    fit_latex->DrawLatex(0.202,0.36,"CMS");
    fit_latex->DrawLatex(0.2,0.30,"#sqrt{s} = 13 TeV");
    fit_latex->DrawLatex(0.2,0.22,"#int Ldt=594 pfb^{-1}");
    fit_latex->DrawLatex(0.2,0.14,"q* #rightarrow q#gamma");
    //it_latex->Draw();

    TPaveText *fit_para = new TPaveText(0.2,0.045,0.474,0.4,"NDC");
    fit_para->AddText("CMS Preliminary  ");
    fit_para->AddText("#sqrt{s} = 13 TeV");
    fit_para->AddText("#int Ldt=594 pb^{-1}");
    fit_para->AddText("q* #rightarrow q#gamma");
    // For Paper ---->

    TPaveText *fit_par1 = new TPaveText(0.2,0.15,0.474,0.4,"NDC");
    fit_par1->AddText("CMS");
    fit_par1->AddText("q* #rightarrow q#gamma");
    fit_par1->AddText("#sqrt{s} = 13 TeV");
    fit_par1->SetTextFont(62); //42
    fit_par1->SetTextAlign(12);
    fit_par1->SetTextSize(0.04);
    fit_par1->SetFillColor(0);
    fit_par1->SetLineColor(0);
    fit_par1->SetFillStyle(0);
    fit_par1->SetBorderSize(0);
//    fit_par1->Draw("same");

    TPaveText *fit_par2 = new TPaveText(0.2,0.045,0.474,0.1,"NDC");
    fit_par2->AddText("#int Ldt=594 pb^{-1}");
    fit_par2->SetTextFont(62); //42
    fit_par2->SetTextAlign(12);
    fit_par2->SetTextSize(0.04);
    fit_par2->SetFillColor(0);
    fit_par2->SetLineColor(0);
    fit_par2->SetFillStyle(0);
    fit_par2->SetBorderSize(0);
//    fit_par2->Draw("same");

    //    fit_para->AddEntry(hGJMass_data,"Data","PL"); 
    //  leg->AddEntry(hGJMass_data,Form("CMS  (%.3f fb^{-1})", 19.6),"PL"); 
    //    fit_para->AddEntry(fit,"Backgraound from Fit","L");
    //    fit_para->AddEntry(h_tmp,"1#sigma Fit Error ","F");
    //    fit_para->AddEntry(hGJMass_mc,"Backround from MC","L");
    //    fit_para->AddText(Form("param 1 : %2.2e",func_Para[0]));
    //    fit_para->AddText(Form("param 2 : %2.2e",func_Para[1]));
    //    fit_para->AddText(Form("param 3 : %2.2e",func_Para[2]));
    //    fit_para->AddText(Form("param 4 : %2.2e",func_Para[3]));
    //    fit_para->AddText(Form("#chi^{2}/ndf : %2.2f/%2.2f",func_FCN[0],func_NDF[0]));
    fit_para->SetTextFont(62); //42
    fit_para->SetTextAlign(13);
    fit_para->SetTextSize(0.04);
    fit_para->SetFillColor(0);
    fit_para->SetLineColor(0);
    fit_para->SetFillStyle(0);
    fit_para->SetBorderSize(0);
//    fit_para->Draw("same");
    //redraw axis

    TPaveText *fn = new TPaveText(0.3,0.2,0.4,0.2,"NDC");
    fn->AddText("Function 7");
    fn->SetFillColor(0);
    fn->SetLineColor(0);
    fn->SetFillStyle(0);
    fn->SetBorderSize(0);
    fn->SetTextFont(42);
    fn->SetTextSize(0.06);
    //    fn->Draw("same");


    TPaveText *pt_qstar1000 = new TPaveText(0.2084548,0.830821,0.4077533,0.9317112,"NDC");
    pt_qstar1000->SetFillColor(0);  
    pt_qstar1000->SetFillStyle(0);  
    pt_qstar1000->SetBorderSize(0); 
    pt_qstar1000->SetTextColor(kOrange+2);
    pt_qstar1000->AddText("q* (1.0 TeV, f = 1.0)");                        
//    pt_qstar1000->Draw();


    p11_1->RedrawAxis();
    p11_1->Update();

    //---- Next PAD
    c1->cd(2);
    p11_2 = (TPad*)c1->GetPad(2);
    p11_2->SetPad(0.01,0.01,0.99,0.24);
    p11_2->SetBottomMargin(0.35);
    p11_2->SetRightMargin(0.05);
    p11_2->SetGridx();
    p11_2->SetGridy();
    c1_2->SetTickx(50);


    TH1F *vFrame2 = p11_2->DrawFrame(p11_1->GetUxmin(), -3.5, p11_1->GetUxmax(), 3.5);

    vFrame2->SetTitle("");
    vFrame2->SetXTitle("M_{#gamma jet} [GeV]");
    vFrame2->GetXaxis()->CenterTitle();
    vFrame2->GetXaxis()->SetTitleSize(0.06);
    vFrame2->GetXaxis()->SetTitleOffset(0.90);
    vFrame2->GetXaxis()->SetTitleSize(0.18);
    vFrame2->GetXaxis()->SetLabelSize(0.14);
    vFrame2->SetYTitle("(Data-Fit)/#sigma_{Data}");
    vFrame2->GetYaxis()->SetTitleSize(0.12);
    vFrame2->GetYaxis()->SetTitleOffset(0.3);
    vFrame2->GetYaxis()->SetLabelSize(0.09);
    vFrame2->GetXaxis()->SetNdivisions(520);

    hRatio->SetLineWidth(0);
    hRatio->SetFillColor(2);
    hRatio->SetLineColor(1);
    hRatio->Draw("SAMEHIST");
    //    gr_Er_qstar1000->Draw("sameC");

/*    c1->SaveAs("FinalFit_DataFitMC_pas.pdf");
    c1->SaveAs("FinalFit_DataFitMC_pas.eps");
    
    TImage *img = TImage::Create();
    img->FromPad(c1);
    img->WriteImage("FinalFit_DataFitMC_pas.png");
    img->WriteImage("FinalFit_DataFitMC_pas.gif");
  */  
    //c1-SaveAs("./FittingPlots/Function1.pdf");

    //    TCanvas* c2 = new TCanvas("c2","Residuals ");
    //    hRatio->SetXTitle("M_{#gamma j} (GeV)");
    //    hRatio->GetXaxis()->CenterTitle();
    //    hRatio->GetXaxis()->SetTitleSize(0.06);
    //    hRatio->SetYTitle("(Data-Fit)/Error");
    //    hRatio->GetYaxis()->CenterTitle();

    //    hRatio->Draw("P E2");

}

Double_t fitFunc1( Double_t *m, Double_t *p) // function 1
{
    double x=m[0]/13000.;
    return p[0]*pow(1.0-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

Double_t fitFunc2( Double_t *m, Double_t *p)  // function 2
{
    double x=m[0]/13000.;
    return p[0]*pow(1.0+x,p[1])/pow(x,p[2]+p[3]*log(x));
}

Double_t fitFunc3( Double_t *m, Double_t *p)  // function 3
{
    double x=m[0]/13000.;
    return p[0]/pow((p[1]+x*p[2]+x*x),p[3]);
}

Double_t fitFunc4( Double_t *m, Double_t *p)  // function 4 
{
    double x=m[0]/13000.;
    return p[0]/pow((p[1]+x),p[2]);
}

/////////
Double_t fitFunc5( Double_t *m, Double_t *p)
{
    double x=m[0]/13000.;
    return p[0]*pow(1.-x+p[3]*x*x,p[1])/pow(m[0],p[2]);
}

// QCD fit function -- alternate 3 parameter fit function -- also used for QCD fit.
Double_t fitFunc6( Double_t *m, Double_t *p)
{
    double x=m[0]/13000.;
    return p[0]*pow(1.-x,p[1])/pow(m[0],p[2]);
}

Double_t fitFunc7( Double_t *m, Double_t *p)
{
    return p[0]/pow(m[0]+p[1],p[2]);
}

Double_t fitFuncForBand( Double_t m, Double_t dp0, Double_t dp1, Double_t dp2, Double_t dp3)  // function 2
{    
    double x = m/13000.;
    double p0 = 4.90863e-07 + 2.17946e-08*dp0;
    double p1 = 9.99858e+00 + 2.23551e-01*dp1;
    double p2 = 1.31627e+01 + 1.75293e-02*dp2;
    double p3 = 1.57651e+00 + 5.87351e-03*dp3;

    return p0*pow(1.0+x,p1)/pow(x,p2+p3*log(x));
}
