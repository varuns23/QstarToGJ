#include "vector.h"


void fittingInvMass_forSushil(){


    // Set TDR Style	
    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // Input Files  for data
//    TFile *inputData = TFile::Open("./rootFiles/Data_53X8T19pb_JetEta3TightID.root", "READ");
    TFile *inputData = TFile::Open("./rootFiles/Data_53X8T19pb_Pt170NoMassCutdEta25TID.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_VarBin_finalPU");
//    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_Data_massVarBin");
    
    TH1F *hRatio = (TH1F*)hGJMass_data->Clone("hRatio");
    hRatio->Reset();

    // Input Files  for MC
    TFile *inputMC1 = TFile::Open("./rootFiles/GJ_53X8T19pb_JetEta3TightID.root", "READ");
    TH1F *hGJMass_mc1 = (TH1F*)inputMC1->Get("h_mass_VarBin_finalPU");
    TFile *inputMC2 = TFile::Open("./rootFiles/DiJet_53X8T19pb_JetEta3TightID.root", "READ");
    TH1F *hGJMass_mc2 = (TH1F*)inputMC2->Get("h_mass_VarBin_finalPU");
    TFile *inputMC3 = TFile::Open("./rootFiles/EWK_53X8T19pb_JetEta3TightID.root", "READ");
    TH1F *hGJMass_mc3 = (TH1F*)inputMC3->Get("h_mass_VarBin_finalPU");

    TH1F *hGJMass_mc= (TH1F*)hGJMass_mc1->Clone("*hGJMass_mc") ;
    hGJMass_mc->Scale(1.3*0.98*0.95);
    hGJMass_mc->Add(hGJMass_mc2,0.98*0.95);
    hGJMass_mc->Add(hGJMass_mc3,0.98*0.95);
    hGJMass_mc->SetLineWidth(2);

    Double_t minXFit = 500.0 ;
    Double_t maxXFit = 3000.0;

    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc2,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
//    fit->SetParameter(0,1.67863e-19);
//    fit->SetParameter(1,5.80213e+01);
//    fit->SetParameter(2,2.93359e+01);
//    fit->SetParameter(3,4.12210e+00);
    fit->SetParameter(0,1.19582e-13);
    fit->SetParameter(1,2.24040e+01);
    fit->SetParameter(2,1.74540e+01);
    fit->SetParameter(3,2.23427e+00);
    fit->SetLineColor(3);

   fit->SetLineWidth(2);
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    hGJMass_data->Fit("fit_gj","R");	
    cout << "NDF = " << fit_gj->GetNDF() << " FCN = " << fit_gj->GetChisquare() << endl;

    Double_t NDF = fit_gj->GetNDF();
    Double_t FCN = fit_gj->GetChisquare();
    Double_t P1 = fit_gj->GetParameter(0);
    Double_t P2 = fit_gj->GetParameter(1);
    Double_t P3 = fit_gj->GetParameter(2);
    Double_t P4 = fit_gj->GetParameter(3);

    gStyle->SetOptFit(0000);

    ///////----------------------
    // Now making the data-fit ratio plot
    float n;
    float n, dm, mass, xl, xh;
    float vx[1200],vy[1200],vexl[1200],vexh[1200],veyl[1200],veyh[1200];
    int i;
    Double_t dataa,fitt, ratioo, errorr;

    for(i=0;i<hGJMass_data->GetNbinsX();i++){
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
	veyl[i] = sqrt(n);
        veyh[i] = sqrt(n);

	if(n !=0){
	    dataa = vy[i];
	    errorr= veyl[i];
	    fitt   = fit->Eval(vx[i],0,0);
//	    ratioo = (dataa-fitt)/dataa;
	    ratioo = (dataa-fitt)/errorr;

	    hRatio->SetBinContent(i+1,ratioo);
	}

    }


    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
    // TO implement the band in the fit function 
    Double_t Arry[16][4] = {
	{+1, +1, +1, +1},
	{+1, -1, -1, -1},
	{-1, +1, -1, -1},
	{-1, -1, +1, -1},
	{-1, -1, -1, +1},
	{+1, +1, -1, -1},
	{-1, -1, +1, +1},
	{+1, -1, -1, +1},
	{-1, +1, +1, -1},
	{+1, -1, +1, -1},
	{-1, +1, +1, -1},
	{+1, +1, +1, -1},
	{+1, +1, -1, +1},
	{+1, -1, +1, +1},
	{-1, +1, +1, +1},
	{-1, -1, -1, -1}
    };

    TH1F *h_Fitmin = new TH1F("h_Fitmin","",100,500,3000);
    TH1F *h_Fitmax = new TH1F("h_Fitmax","",100,500,3000);

    vector<double> X_Fitmin, X_Fitmax;

    for(Double_t Errmass = 512.5 ; Errmass < 3000 ; Errmass += 25 ){
	Double_t x_Fitmax = 0.0 ;
	Double_t x_Fitmin = 1000000.0 ;

	for(Int_t i = 0 ; i < 16 ; ++i){

	    Double_t Returnfunction = fitFuncForBand(Errmass,Arry[i][0],Arry[i][1],Arry[i][2],Arry[i][3]);

	    if(Returnfunction > x_Fitmax ) x_Fitmax = Returnfunction ;
	    if(Returnfunction < x_Fitmin ) x_Fitmin = Returnfunction ;
	}

	h_Fitmin->Fill(Errmass,x_Fitmin);
	h_Fitmax->Fill(Errmass,x_Fitmax);

    }

    THStack *h_fitErrorBand = new THStack("h_fitErrorBand","");
    TH1F *h_tmp = (TH1F*)h_Fitmax->Clone("h_tmp");
    h_tmp->Add(h_Fitmin,-1);
    h_fitErrorBand->Add(h_Fitmin);
    h_fitErrorBand->Add(h_tmp);

    h_Fitmin->SetLineColor(0);
    h_tmp->SetLineColor(5);
    h_tmp->SetFillColor(5);
    h_tmp->SetLineWidth(1);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=

    
    TCanvas* c1 = new TCanvas("c1","DijetMass Cross Section with Fit");
    c1->Divide(1,2,0,0,0);
    c1->cd(1);
    p11_1 = (TPad*)c1->GetPad(1);
    p11_1->SetPad(0.01,0.23,0.99,0.98);
    p11_1->SetLogy();
    p11_1->SetRightMargin(0.05);
    p11_1->SetTopMargin(0.05);

    TH1F *vFrame = p11_1->DrawFrame(minXFit,0.000000001,maxXFit,10.0);
    vFrame->SetTitle("");
    vFrame->SetXTitle("Mass M_{#gamma j} (GeV)");
    vFrame->GetXaxis()->SetTitleSize(0.06);
    vFrame->SetYTitle("");
    
    hGJMass_data->SetTitle("");
    hGJMass_data->SetLineColor(1);
    hGJMass_data->SetFillColor(1);
    hGJMass_data->SetMarkerColor(1);
    hGJMass_data->SetMarkerStyle(20);
    hGJMass_data->GetYaxis()->SetTitle("");
    hGJMass_data->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    
  
    //+++++++++++ for Error ++++++++++++++++++++++++++++++++++
    TH1F *h_err_bkg = hGJMass_mc->Clone("h_err_bkg");

    for( int k = 1 ; k <= h_err_bkg->GetNbinsX(); ++k){
	h_err_bkg->SetBinError(k, sqrt( pow(((0.03*h_err_bkg->GetBinContent(k))/0.93),2)+ pow(h_err_bkg->GetBinError(k),2) + 0.008*(h_err_bkg->GetBinContent(k)**2) + 
		    0.002*(h_err_bkg->GetBinContent(k)**2) + 0.026*(h_err_bkg->GetBinContent(k)**2) ));
    }  
    h_err_bkg->SetFillColor(kRed+2); h_err_bkg->SetFillStyle(3003); 
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
//    h_fitErrorBand->Draw("C");
//    h_err_bkg->Draw("same E2");
//    hGJMass_mc->Draw("sameHIST");
    hGJMass_data->Draw("PZ E1");
    fit->Draw("sameAZ ");

    //    hRatio->Draw("hist");
    
    TLegend *leg = new TLegend(0.2,0.1,0.5,0.3);
    leg->SetTextSize(0.04);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hGJMass_data,"Data","PL"); 

    leg->AddEntry(fit,"Fit","L");
    leg->AddEntry(h_tmp,"Fit Error ","F");
    leg->AddEntry(hGJMass_mc,"Bkg MC","L");
    leg->AddEntry(h_err_bkg,"Background Uncertainity ","F");
    leg->Draw("same");

    TPaveText *fit_para = new TPaveText(0.6,.5,0.95,0.92,"NDC");
//    fit_para->AddText("#sqrt{s} = 8 TeV");
    fit_para->AddText(Form("param 1 : %2.2e",P1));
    fit_para->AddText(Form("param 2 : %2.2e",P2));
    fit_para->AddText(Form("param 3 : %2.2e",P3));
    fit_para->AddText(Form("param 4 : %2.2e",P4));
    fit_para->AddText(Form("#chi^{2} : %2.2f",FCN));
    fit_para->AddText(Form("NDF : %2.2f",NDF));
    fit_para->SetFillColor(0);
    fit_para->SetLineColor(0);
    fit_para->SetFillStyle(0);
    fit_para->SetBorderSize(0);
    fit_para->SetTextFont(42);
    fit_para->SetTextSize(0.04);
    fit_para->Draw("same");
   //redraw axis
    TPaveText *fn = new TPaveText(0.3,0.2,0.4,0.2,"NDC");
    fn->AddText("Function 7");
    fn->SetFillColor(0);
    fn->SetLineColor(0);
    fn->SetFillStyle(0);
    fn->SetBorderSize(0);
    fn->SetTextFont(42);
    fn->SetTextSize(0.05);
//    fn->Draw("same");

    p11_1->RedrawAxis();
    p11_1->Update();

 
    //---- Next PAD
    c1->cd(2);
    p11_2 = (TPad*)c1->GetPad(2);
    p11_2->SetPad(0.01,0.02,0.99,0.24);
    p11_2->SetBottomMargin(0.35);
    p11_2->SetRightMargin(0.05);
    p11_2->SetGridx();
    p11_2->SetGridy();
    c1_2->SetTickx(50);


    TH1F *vFrame2 = p11_2->DrawFrame(p11_1->GetUxmin(), -3.0, p11_1->GetUxmax(), 3.0);
           
    vFrame2->SetTitle("");
    vFrame2->SetXTitle("M_{#gamma j} (GeV)");
    vFrame2->GetXaxis()->CenterTitle();
    vFrame2->GetXaxis()->SetTitleSize(0.06);
    vFrame2->SetYTitle("(Data-Fit)/Error");
    vFrame2->GetYaxis()->SetTitleSize(0.12);
    vFrame2->GetYaxis()->SetLabelSize(0.07);
    vFrame2->GetYaxis()->SetTitleOffset(0.50);
    vFrame2->GetXaxis()->SetTitleOffset(0.90);
    vFrame2->GetXaxis()->SetTitleSize(0.18);
    vFrame2->GetXaxis()->SetLabelSize(0.1);


    hRatio->SetLineWidth(0);
    hRatio->SetFillColor(2);
    hRatio->SetLineColor(1);
    hRatio->Draw("SAMEHIST");

}

Double_t fitFunc2( Double_t *m, Double_t *p)  // function 2
{
    double x=m[0]/8000.;
    return p[0]*pow(1.0+x,p[1])/pow(x,p[2]+p[3]*log(x));
}

Double_t fitFuncForBand( Double_t m, Double_t dp0, Double_t dp1, Double_t dp2, Double_t dp3)  // function 2
{    
    double x = m/8000.;
    double p0 = 2.77519e-10 + 1.39901e-11*dp0;
    double p1 = 2.08816e+01 + 2.73577e-01*dp1;
    double p2 = 1.69923e+01 + 1.91084e-02*dp2;
    double p3 = 2.16517e+00 + 6.16235e-03*dp3;

    return p0*pow(1.0+x,p1)/pow(x,p2+p3*log(x));
}
