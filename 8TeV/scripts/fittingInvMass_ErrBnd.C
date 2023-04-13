#include "vector.h"


void fittingInvMass_ErrBnd(){


    // Set TDR Style	
    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

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

    TH1F *hGJMass_mc= (TH1F*)hGJMass_mc1->Clone("*hGJMass_mc") ;
    hGJMass_mc->Scale(1.42*0.98*0.95);
    hGJMass_mc->Add(hGJMass_mc2,0.98*0.95);
    hGJMass_mc->Add(hGJMass_mc3,0.98*0.95);
    hGJMass_mc->SetLineWidth(2);

//    Double_t minXFit = 547.0 ;
    Double_t minXFit = 560.0 ;
//    Double_t maxXFit = 2965.0;
//    Double_t maxXFit = 2714.0;
    Double_t maxXFit = 3000.0;
/*
    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc1,minXFit,maxXFit,4); // function 1
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,1.14050e-07);
    fit->SetParameter(1,-6.57567e+00);
    fit->SetParameter(2,1.56204e+01);
    fit->SetParameter(3,2.32949e+00);
    fit->SetLineColor(2);
*/

    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc2,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,4.90863e-07);
    fit->SetParameter(1,9.99858e+00);
    fit->SetParameter(2,1.31627e+01);
    fit->SetParameter(3,1.57651e+00);
    fit->SetLineColor(3);

/*  
    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc3,minXFit,maxXFit,4); // function 3 Not converging
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,1.70138e+05);
    fit->SetParameter(1,7.61449e-01);
    fit->SetParameter(2,8.05594e+00);
    fit->SetParameter(3,9.68175e+00);
    fit->SetLineColor(4);
*/
/*    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc4,minXFit,maxXFit,3); // function 4 Not converging
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,2.92398e-03);
    fit->SetParameter(1,4.01131e-01);
    fit->SetParameter(2,2.04757e+01);
    fit->SetLineColor(5);
*/
/*    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc5,minXFit,maxXFit,4); // function 5
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,1.32682e+08);
    fit->SetParameter(1,-2.92762e+00);
    fit->SetParameter(2,1.18654e+00);
    fit->SetParameter(3,2.02310e+02);
    fit->SetLineColor(6);
*/
/*    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc6,minXFit,maxXFit,3); // function 6   Not Converging
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,5.29592e+08);
    fit->SetParameter(1,2.38179e+01);
    fit->SetParameter(2,1.38685e+00);
    fit->SetLineColor(7);
*/
/*    // Fit to data    
    TF1 *fit = new TF1("fit_gj",fitFunc7,minXFit,maxXFit,3); // function 7   Not Converging
    gStyle->SetOptFit(0000000); 
    fit->SetParameter(0,2.71844e+78);
    fit->SetParameter(1,3.26078e+03);
    fit->SetParameter(2,2.07350e+01);
    fit->SetLineColor(8);
*/
   fit->SetLineWidth(2);
    //    fit->SetLineColor(4);
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
    Int_t LowestMassXBin = 33;
    Int_t HighestMassXBin = 72;

    for(int xBin = LowestMassXBin; xBin <= HighestMassXBin; xBin++ ){


	double data_fromFit = fit->Eval((hGJMass_data->GetBinCenter(i)),0,0);
	double	x_Fitmax =  data_fromFit + sqrt(data_fromFit);
	double	x_Fitmin =  data_fromFit - sqrt(data_fromFit);

	h_Fitmin->Fill((hGJMass_data->GetBinCenter(i)),x_Fitmin);
	h_Fitmax->Fill((hGJMass_data->GetBinCenter(i)),x_Fitmax);

    }

/*    for(Double_t Errmass = 562.5 ; Errmass < 3000 ; Errmass += 25 ){
	Double_t x_Fitmax = 0.0 ;
	Double_t x_Fitmin = 1000000.0 ;

	for(Int_t i = 0 ; i < 16 ; ++i){

	    Double_t Returnfunction = fitFuncForBand(Errmass,Arry[i][0],Arry[i][1],Arry[i][2],Arry[i][3]);

	    if(Returnfunction > x_Fitmax ) x_Fitmax = Returnfunction ;
	    if(Returnfunction < x_Fitmin ) x_Fitmin = Returnfunction ;
	}
	
	h_Fitmin->Fill(Errmass,x_Fitmin);
	h_Fitmax->Fill(Errmass,x_Fitmax);
    }*/

    THStack *h_fitErrorBand = new THStack("h_fitErrorBand","");
    TH1F *h_tmp = (TH1F*)h_Fitmax->Clone("h_tmp");
    h_tmp->Add(h_Fitmin,-1);
    h_fitErrorBand->Add(h_Fitmin);
    h_fitErrorBand->Add(h_tmp);

    h_Fitmin->SetLineColor(0);
    h_tmp->SetLineColor(5);
    h_tmp->SetFillColor(5);
    h_tmp->SetLineWidth(1);
    h_fitErrorBand->SetMinimum(1);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=

    
    //Dijet Mass Cross Section with Fit	
    TCanvas* c1 = new TCanvas("c1","DijetMass Cross Section with Fit");
//    gStyle->SetOptFile(0);
//    gStyle->SetOptStat(0);
//    gStyle->SetPadBorderMode(0);
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
//    hGJMass_data->GetYaxis()->SetRangeUser(0.001,30000);
    
  
    //+++++++++++ for Error ++++++++++++++++++++++++++++++++++
    TH1F *h_err_bkg = hGJMass_mc->Clone("h_err_bkg");

    for( int k = 1 ; k <= h_err_bkg->GetNbinsX(); ++k){
//	std::cout<<"Bin Content : "<<h_err_bkg->GetBinContent(k)<<"      Err Before = "<<h_err_bkg->GetBinError(k)<<std::endl;
	h_err_bkg->SetBinError(k, sqrt( pow(((0.03*h_err_bkg->GetBinContent(k))/0.93),2)+ pow(h_err_bkg->GetBinError(k),2) + 0.008*(h_err_bkg->GetBinContent(k)**2) + 
		    0.002*(h_err_bkg->GetBinContent(k)**2) + 0.026*(h_err_bkg->GetBinContent(k)**2) ));
//	std::cout<<" +++++++++Err After = "<<h_err_bkg->GetBinError(k)<<std::endl;
    }  
//    h_err_bkg->SetFillColor(kRed+2); h_err_bkg->SetFillStyle(3002); 
    h_err_bkg->SetFillColor(kRed+2); h_err_bkg->SetFillStyle(3003); 

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    h_fitErrorBand->Draw("C");
    h_err_bkg->Draw("same E2");
    hGJMass_mc->Draw("sameHIST");
    fit->Draw("sameHiST");
    hGJMass_data->Draw("samePZ E1");

    //    hRatio->Draw("hist");
    
    TLegend *leg = new TLegend(0.2,0.1,0.5,0.3);
//    leg->SetTextSize(0.03146853);
    leg->SetTextSize(0.04);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hGJMass_data,"Data","PL"); 
  //  leg->AddEntry(hGJMass_data,Form("CMS  (%.3f fb^{-1})", 19.6),"PL"); 
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

    //c1-SaveAs("./FittingPlots/Function1.eps");
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
    double x=m[0]/8000.;
    return p[0]*pow(1.0-x,p[1])/pow(x,p[2]+p[3]*log(x));
}

Double_t fitFunc2( Double_t *m, Double_t *p)  // function 2
{
    double x=m[0]/8000.;
    return p[0]*pow(1.0+x,p[1])/pow(x,p[2]+p[3]*log(x));
}

Double_t fitFunc3( Double_t *m, Double_t *p)  // function 3
{
    double x=m[0]/8000.;
    return p[0]/pow((p[1]+x*p[2]+x*x),p[3]);
}

Double_t fitFunc4( Double_t *m, Double_t *p)  // function 4 
{
    double x=m[0]/8000.;
    return p[0]/pow((p[1]+x),p[2]);
}

/////////
Double_t fitFunc5( Double_t *m, Double_t *p)
{
        double x=m[0]/8000.;
	return p[0]*pow(1.-x+p[3]*x*x,p[1])/pow(m[0],p[2]);
}

// QCD fit function -- alternate 3 parameter fit function -- also used for QCD fit.
Double_t fitFunc6( Double_t *m, Double_t *p)
{
        double x=m[0]/8000.;
	    return p[0]*pow(1.-x,p[1])/pow(m[0],p[2]);
}
 
Double_t fitFunc7( Double_t *m, Double_t *p)
{
        return p[0]/pow(m[0]+p[1],p[2]);
}

Double_t fitFuncForBand( Double_t m, Double_t dp0, Double_t dp1, Double_t dp2, Double_t dp3)  // function 2
{    
    double x = m/8000.;
    double p0 = 4.90863e-07 + 2.17946e-08*dp0;
    double p1 = 9.99858e+00 + 2.23551e-01*dp1;
    double p2 = 1.31627e+01 + 1.75293e-02*dp2;
    double p3 = 1.57651e+00 + 5.87351e-03*dp3;

    return p0*pow(1.0+x,p1)/pow(x,p2+p3*log(x));
}
