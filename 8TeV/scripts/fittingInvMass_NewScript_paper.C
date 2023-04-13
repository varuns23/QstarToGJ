#include "vector.h"


void fittingInvMass_NewScript_paper(){


    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    gStyle->SetOptTitle(0);
    gStyle->SetPadLeftMargin(0.10);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.45);


    // Input Files  for data
    TFile *inputData = TFile::Open("./rootFiles/Data_53X8T19pb_Pt170Mass560dEta20TID.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_VarBin_finalPU");

    TH1F *hRatio = (TH1F*)hGJMass_data->Clone("hRatio");
    hRatio->Reset();

    // Input Files  for MC
    TFile *inputMC1 = TFile::Open("./rootFiles/GJ_53X8T19pb_Pt170Mass560dEta20TID.root", "READ");
    TH1F *hGJMass_mc1 = (TH1F*)inputMC1->Get("h_mass_VarBin_finalPU");
    TFile *inputMC2 = TFile::Open("./rootFiles/DiJet_53X8T19pb_Pt170Mass560dEta20TID.root", "READ");
    TH1F *hGJMass_mc2 = (TH1F*)inputMC2->Get("h_mass_VarBin_finalPU");
    TFile *inputMC3 = TFile::Open("./rootFiles/EWK_53X8T19pb_Pt170Mass560dEta20TID.root", "READ");
    TH1F *hGJMass_mc3 = (TH1F*)inputMC3->Get("h_mass_VarBin_finalPU");

    TFile *inputSignal1 = TFile::Open("./rootFiles/QstarToGJ_M_1000.root", "READ");
    TH1F *h_qstar_1000 = (TH1F*)inputSignal1->Get("h_mass_VarBin_finalPU");

    TFile *inputSignal2 = TFile::Open("./rootFiles/QstarToGJ_M_1500.root", "READ");
    TH1F *h_qstar_1500 = (TH1F*)inputSignal2->Get("h_mass_VarBin_finalPU");

    TFile *inputSignal3 = TFile::Open("./rootFiles/QstarToGJ_M_fhalf_2500.root", "READ");
    TH1F *h_qstar_f2500 = (TH1F*)inputSignal3->Get("h_mass_VarBin_finalPU");

    TH1F *hGJMass_mc= (TH1F*)hGJMass_mc1->Clone("*hGJMass_mc") ;

    //    hGJMass_mc->Scale(1.42*0.98*0.95);
    hGJMass_mc->Scale(1.33*0.98*0.95);
    hGJMass_mc->Add(hGJMass_mc2,1.34*0.98*0.95);
    hGJMass_mc->Add(hGJMass_mc3,0.98*0.95);
    hGJMass_mc->SetLineWidth(2);

    Double_t minXFit = 560.0 ;
    Double_t maxXFit = 3100.0;

    vector<double> func_Para, func_NDF, func_FCN;

    // Fit to data    
    TF1 *fit = new TF1("fit",fitFunc1,minXFit,maxXFit,4); 
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,2.74281e-03);
    fit->SetParameter(1,3.33458e+00);
    fit->SetParameter(2,8.06973e+00);
    fit->SetParameter(3,7.47309e-01);

    fit->SetLineColor(3);
    fit->SetLineWidth(2);

    for(int a = 0; a < 10; ++a) hGJMass_data->Fit("fit","R");
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
    std::vector<double> qstar1000_Y, qstar1000_X, qstar1000_Err;
    std::vector<double> qstar1500_Y, qstar1500_X, qstar1500_Err;
    std::vector<double> qstarf2500_Y, qstarf2500_X, qstarf2500_Err;

    //    for(i=0;i<hGJMass_data->GetNbinsX();i++){
    for(i=0;i<72;i++){

	n    = hGJMass_data->GetBinContent(i+1);
	dm   = hGJMass_data->GetBinWidth(i+1);
	mass = hGJMass_data->GetBinCenter(i+1);
	xl   = hGJMass_data->GetBinLowEdge(i+1);
	xh   = xl+dm; 
	vx[i]   = (xl+xh)/2.;
	vy[i]   = n;

	vexl[i] = dm/2.;
	vexh[i] = dm/2.;

	//	veyl[i] = sqrt(n)/(lumi*dm);
	//      veyh[i] = sqrt(n)/(lumi*dm);
	//	veyl[i] = sqrt(n);
	//	veyh[i] = sqrt(n);
	double LOW = (n == 0) ? 0 : (ROOT::Math::gamma_quantile(alpha/2,n,1.)); 
	double UP = ROOT::Math::gamma_quantile_c(alpha/2,n+1,1) ;

	veyl[i] = n - LOW;
	veyh[i] = UP - n;

	if(i == 68 || i == 70){
	    veyl[i] = 0.0;
	    veyh[i] = 1.8;

	    dataa = vy[i];
	    errorr= veyh[i];
	    fitt   = fit->Eval(vx[i],0,0);
	    ratioo = (dataa-fitt)/errorr;
	    hRatio->SetBinContent(i+1,ratioo);
	}

	if(i >= 72){
	    veyl[i] = 0.0;
	    veyh[i] = 0.0;
	}


	if(n !=0){
	    dataa = vy[i];
	    errorr= veyl[i];
	    fitt   = fit->Eval(vx[i],0,0);
	    //	    ratioo = (dataa-fitt)/dataa;
	    ratioo = (dataa-fitt)/errorr;
	    hRatio->SetBinContent(i+1,ratioo);
	}

	if(mass > 0.8*1000 && mass < 1.3*1000){
	    qstar1000_Y.push_back(h_qstar_1000->GetBinContent(i+1) + (hGJMass_mc->GetBinContent(i+1)*0.62));
	    qstar1000_X.push_back(mass);
	}

	if(mass > 0.85*1500 && mass < 1.2*1500){
	    qstar1500_Y.push_back(h_qstar_1500->GetBinContent(i+1) + (hGJMass_mc->GetBinContent(i+1)*0.62));
	    qstar1500_X.push_back(mass);
	}

	if(mass > 0.9*2500 && mass < 1.12*2500){
	    qstarf2500_Y.push_back(h_qstar_f2500->GetBinContent(i+1) + (hGJMass_mc->GetBinContent(i+1)*0.9));
	    qstarf2500_X.push_back(mass);
	}

    }

    // Graphs to Make Signal points ----------
    TGraphAsymmErrors *g = new TGraphAsymmErrors(i,vx,vy,vexl,vexh,veyl,veyh);

    TGraph* gr_qstar1000 = new TGraph(qstar1000_X.size() , &qstar1000_X[0], &qstar1000_Y[0]);
    gr_qstar1000->SetLineColor(kOrange+2);
    gr_qstar1000->SetLineStyle(7);
    gr_qstar1000->SetLineWidth(2);
    gr_qstar1000->SetTitle();


    TGraph* gr_qstar1500 = new TGraph(qstar1500_X.size() , &qstar1500_X[0], &qstar1500_Y[0]); 
    gr_qstar1500->SetLineColor(kCyan+2);
    gr_qstar1500->SetLineStyle(7);
    gr_qstar1500->SetLineWidth(2);
    gr_qstar1500->SetTitle();

    TGraph* gr_qstarf2500 = new TGraph(qstarf2500_X.size() , &qstarf2500_X[0], &qstarf2500_Y[0]); 
    gr_qstarf2500->SetLineColor(kMagenta+2);
    gr_qstarf2500->SetLineStyle(7);
    gr_qstarf2500->SetLineWidth(2);
    gr_qstarf2500->SetTitle();


    //Dijet Mass Cross Section with Fit	
    TCanvas* c1 = new TCanvas("c1","GammajetMass Cross Section with Fit",800,700);

    TPad *p1_1 = new TPad("p1_1","Photon + jet mass",0.0,0.24,1,1);
    p1_1->SetLogy();
    p1_1->SetBottomMargin(0);
    p1_1->Draw();
    p1_1->cd();


    /*
       c1->Divide(1,2,0,0,0);
       c1->cd(1);
       p11_1 = (TPad*)c1->GetPad(1);
       p11_1->SetPad(0.01,0.24,0.99,0.98);
       p11_1->SetLogy();
       p11_1->SetRightMargin(0.05);
       p11_1->SetTopMargin(0.06);

       TH1F *vFrame = p11_1->DrawFrame(minXFit,0.000000001,maxXFit,10.0);
       vFrame->SetTitle("");
       vFrame->GetXaxis()->SetTitleSize(0.06);
       vFrame->SetYTitle("");
     */

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

    g->SetMarkerColor(1);
    g->SetLineColor(1);

    hGJMass_mc->Draw("Hist");
    g->Draw("samePE1");
    fit->Draw("same");
    gr_qstar1000->Draw("sameC");
    gr_qstar1500->Draw("sameC");
    gr_qstarf2500->Draw("sameC");

    TLegend *leg = new TLegend(0.62,.65,0.9,0.9);
    leg->SetTextFont(62);//42
    leg->SetTextSize(0.04);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hGJMass_data,"Data","PL"); 
    leg->AddEntry(fit,"Background from Fit","L");
    leg->AddEntry(hGJMass_mc,"Background from MC","L");
    leg->Draw("same");


    TLatex *fit_latex = new TLatex();
    fit_latex->SetTextAlign(13);
    fit_latex->SetTextFont(62);
    fit_latex->SetNDC();
    fit_latex->SetTextSize(0.04);
    fit_latex->DrawLatex(0.202,0.36,"CMS Preliminary");
    fit_latex->DrawLatex(0.2,0.30,"#sqrt{s} = 8 TeV");
    fit_latex->DrawLatex(0.2,0.22,"#int Ldt=19.7 fb^{-1}");
    fit_latex->DrawLatex(0.2,0.14,"q* #rightarrow q#gamma");

    TPaveText *pt_qstar1000 = new TPaveText(0.2084548,0.830821,0.4077533,0.9317112,"NDC");
    pt_qstar1000->SetFillColor(0);  
    pt_qstar1000->SetFillStyle(0);  
    pt_qstar1000->SetBorderSize(0); 
    pt_qstar1000->SetTextColor(kOrange+2);
    pt_qstar1000->AddText("q* (1.0 TeV, f = 1.0)");                        
    pt_qstar1000->Draw();

    TPaveText *pt_qstar1500 = new TPaveText(0.372449,0.6587141,0.5763028,0.7615826,"NDC");
    pt_qstar1500->SetFillColor(0);  
    pt_qstar1500->SetFillStyle(0);  
    pt_qstar1500->SetBorderSize(0); 
    pt_qstar1500->SetTextColor(kCyan+2);
    pt_qstar1500->AddText("q* (1.5 TeV, f = 1.0)");                        
    pt_qstar1500->Draw();

    TPaveText *pt_qstarf2500 = new TPaveText(0.7015762,0.3085658,0.9077077,0.3876954,"NDC");
    pt_qstarf2500->SetFillColor(0);  
    pt_qstarf2500->SetFillStyle(0);  
    pt_qstarf2500->SetBorderSize(0); 
    pt_qstarf2500->SetTextColor(kMagenta+2);
    pt_qstarf2500->AddText("q* (2.5 TeV, f = 0.5)");                        
    pt_qstarf2500->Draw();

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
    /*
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
     */
    hRatio->SetXTitle("M_{#gamma,jet} [GeV]");
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitleSize(0.06);
    hRatio->GetXaxis()->SetTitleOffset(0.90);
    hRatio->GetXaxis()->SetTitleSize(0.22);  //18
    hRatio->GetXaxis()->SetLabelSize(0.14);
    hRatio->GetXaxis()->SetNdivisions(520);
    hRatio->GetXaxis()->SetRangeUser(minXFit,maxXFit);

    hRatio->SetYTitle("#frac{(Data-Fit)}{#sigma_{Data}}");
    hRatio->GetYaxis()->SetRangeUser(-3.5,3.5);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.3);
    hRatio->GetYaxis()->SetLabelSize(0.09);

    hRatio->SetLineWidth(0);
    hRatio->SetFillColor(2);
    hRatio->SetLineColor(1);
    hRatio->Draw("HIST");
    //    gr_Er_qstar1000->Draw("sameC");

    /*
       c1->SaveAs("FinalFit_DataFitMC_pas.pdf");
       c1->SaveAs("FinalFit_DataFitMC_pas.eps");

       TImage *img = TImage::Create();
       img->FromPad(c1);
       img->WriteImage("FinalFit_DataFitMC_pas.png");
       img->WriteImage("FinalFit_DataFitMC_pas.gif");
     */
}

Double_t fitFunc1( Double_t *m, Double_t *p) // function 1
{
    double x=m[0]/8000.;
    return p[0]*pow(1.0-x,p[1])/pow(x,p[2]+p[3]*log(x));
}
