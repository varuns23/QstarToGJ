#include "vector.h"

void fittingInvMass_DifferentRangeDiffPAD(){

    // Set TDR Style	
    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // Input Files  for data
    //    TFile *inputData = TFile::Open("./rootFiles/Data_53X8T19pb_Pt170Mass530dEta20TID.root", "READ");
    TFile *inputData = TFile::Open("/home/varuns/Research/work/LHC_CMS_13TeV/QstarToGJ/Analysis/postAnalyzerResults/Pt190jetEta2p4dphi2p0M560/rootFiles/Run2015D_PromptReco.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_VarBin_MassCut");

    
    TH1F *hRatio1 = (TH1F*)hGJMass_data->Clone("hRatio1");
    hRatio1->Reset();
    TH1F *hRatio2 = (TH1F*)hGJMass_data->Clone("hRatio2");
    hRatio2->Reset();
    TH1F *hRatio3 = (TH1F*)hGJMass_data->Clone("hRatio3");
    hRatio3->Reset();
    TH1F *hRatio4 = (TH1F*)hGJMass_data->Clone("hRatio4");
    hRatio4->Reset();
    TH1F *hRatio5 = (TH1F*)hGJMass_data->Clone("hRatio5");
    hRatio5->Reset();

    Double_t minXFit  = 560.0 ;
    Double_t maxXFit1 = 806.0;
    Double_t maxXFit2 = 1017.0;
    Double_t maxXFit3 = 1315.0;
    Double_t maxXFit4 = 1807.0;
    Double_t maxXFit5 = 3237.0;
    Double_t maxXFit = 3237.0;

    // Fit to data    
    TF1 *fit1 = new TF1("fit1",fitFunc,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    setfitFuncParameters(fit1);
    fit1->SetLineColor(2);
    fit1->SetLineWidth(2);

    for( int a = 0; a < 10; ++a){
	hGJMass_data->Fit("fit1","","",minXFit,maxXFit1);	
	//hGJMass_data->Fit("fit1","R"); // R is used to fit the range mentioned in function
    }
    Double_t fit1NDF = fit1->GetNDF();
    Double_t fit1FCN = fit1->GetChisquare();
    ///---------------------------------------------------------------------------------------
    TF1 *fit2 = new TF1("fit2",fitFunc,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    setfitFuncParameters(fit2);
    fit2->SetLineColor(3);
    fit2->SetLineWidth(2);

    for( int a = 0; a < 10; ++a){
	hGJMass_data->Fit("fit2","","",minXFit,maxXFit2);	
    }
    Double_t fit2NDF = fit2->GetNDF();
    Double_t fit2FCN = fit2->GetChisquare();
    ///---------------------------------------------------------------------------------------
    TF1 *fit3 = new TF1("fit3",fitFunc,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    setfitFuncParameters(fit3);
    fit3->SetLineColor(4);
    fit3->SetLineWidth(2);

    for( int a = 0; a < 10; ++a){
	hGJMass_data->Fit("fit3","","",minXFit,maxXFit3);	
    }
    Double_t fit3NDF = fit3->GetNDF();
    Double_t fit3FCN = fit3->GetChisquare();
    ///---------------------------------------------------------------------------------------
    TF1 *fit4 = new TF1("fit4",fitFunc,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    setfitFuncParameters(fit4);
    fit4->SetLineColor(5);
    fit4->SetLineWidth(2);

    for( int a = 0; a < 10; ++a){
	hGJMass_data->Fit("fit4","","",minXFit,maxXFit4);	
    }
    Double_t fit4NDF = fit4->GetNDF();
    Double_t fit4FCN = fit4->GetChisquare();
    ///---------------------------------------------------------------------------------------
    TF1 *fit5 = new TF1("fit5",fitFunc,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    setfitFuncParameters(fit5);
    fit5->SetLineColor(6);
    fit5->SetLineWidth(2);

    for( int a = 0; a < 10; ++a){
	hGJMass_data->Fit("fit5","","",minXFit,maxXFit5);	
    }
    Double_t fit5NDF = fit5->GetNDF();
    Double_t fit5FCN = fit5->GetChisquare();
    ///---------------------------------------------------------------------------------------

    ///////----------------------
    // Now making the data-fit ratio plot
    float n;
    float n, dm, mass, xl, xh;
    float vx[1200],vy[1200],veyl[1200];
    int i;
    Double_t dataa,fitt1,fitt2,fitt3,fitt4,fitt5, ratioo1,ratioo2,ratioo3,ratioo4,ratioo5, errorr;

    for(i=0;i<hGJMass_data->GetNbinsX();i++){
	n    = hGJMass_data->GetBinContent(i+1);
	dm   = hGJMass_data->GetBinWidth(i+1);
	mass = hGJMass_data->GetBinCenter(i+1);
	xl   = hGJMass_data->GetBinLowEdge(i+1);
	xh   = xl+dm; 
	vx[i]   = (xl+xh)/2.;
	vy[i]   = n;

	veyl[i] = sqrt(n);

	if(n !=0){
	    dataa = vy[i];
	    errorr= veyl[i];
	    fitt1   = fit1->Eval(vx[i],0,0);
	    fitt2   = fit2->Eval(vx[i],0,0);
	    fitt3   = fit3->Eval(vx[i],0,0);
	    fitt4   = fit4->Eval(vx[i],0,0);
	    fitt5   = fit5->Eval(vx[i],0,0);
	    ratioo1 = (dataa-fitt1)/errorr;
	    ratioo2 = (dataa-fitt2)/errorr;
	    ratioo3 = (dataa-fitt3)/errorr;
	    ratioo4 = (dataa-fitt4)/errorr;
	    ratioo5 = (dataa-fitt5)/errorr;

//	    if(vx[i] <= maxXFit1) ratioo1 = (dataa-fitt1)/errorr;
//	    if(vx[i] <= maxXFit2) ratioo2 = (dataa-fitt2)/errorr;
//	    if(vx[i] <= maxXFit3) ratioo3 = (dataa-fitt3)/errorr;
//	    if(vx[i] <= maxXFit4) ratioo4 = (dataa-fitt4)/errorr;
//	    if(vx[i] <= maxXFit5) ratioo5 = (dataa-fitt5)/errorr;

	    hRatio1->SetBinContent(i+1,ratioo1);
	    hRatio2->SetBinContent(i+1,ratioo2);
	    hRatio3->SetBinContent(i+1,ratioo3);
	    hRatio4->SetBinContent(i+1,ratioo4);
	    hRatio5->SetBinContent(i+1,ratioo5);
	}

    }


    //    TCanvas* c1 = new TCanvas("c1","PhotonJet Invariant Mass with Fit",100,100,700,700);
    //+++++++++++++ First Canvas starts ++++++++++++++++=

    TCanvas* c1 = new TCanvas("c1","PhotonJet Invariant Mass with Fit",600,650);
    c1->SetLogy();
    c1->SetPad(0.01,0.00,0.99,0.98);
    c1->SetLogy();
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.1);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.1);

    hGJMass_data->SetTitle("");
    hGJMass_data->SetLineColor(1);
    hGJMass_data->SetFillColor(1);
    hGJMass_data->SetMarkerColor(1);
    hGJMass_data->SetMarkerStyle(20);
    hGJMass_data->GetYaxis()->SetTitle("Events");
    hGJMass_data->GetYaxis()->SetLabelSize(0.04);
    hGJMass_data->GetYaxis()->SetTitleSize(0.05);
    hGJMass_data->GetYaxis()->SetTitleOffset(1);
    hGJMass_data->GetXaxis()->SetTitle("Mass M_{#gamma jet} (GeV)");
    hGJMass_data->GetXaxis()->SetLabelSize(0.04);
    hGJMass_data->GetXaxis()->SetTitleSize(0.04);
    hGJMass_data->GetXaxis()->SetTitleOffset(1.1);
    hGJMass_data->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    
    hGJMass_data->Draw("PZ E1");
    fit5->Draw("sameAZ ");
    fit4->Draw("sameAZ ");
    fit3->Draw("sameAZ ");
    fit2->Draw("sameAZ ");
    fit1->Draw("sameAZ ");
    
    TLegend *leg = new TLegend(0.4,0.7,0.85,0.92);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(hGJMass_data,"Data","PL"); 

    leg->AddEntry(fit1,Form("Fit (560 - %4.0f), #chi^{2}/NDF : %2.2f",maxXFit1,fit1FCN/fit1NDF),"L");
    leg->AddEntry(fit2,Form("Fit (560 - %4.0f), #chi^{2}/NDF : %2.2f",maxXFit2,fit2FCN/fit2NDF),"L");
    leg->AddEntry(fit3,Form("Fit (560 - %4.0f), #chi^{2}/NDF : %2.2f",maxXFit3,fit3FCN/fit3NDF),"L");
    leg->AddEntry(fit4,Form("Fit (560 - %4.0f), #chi^{2}/NDF : %2.2f",maxXFit4,fit4FCN/fit4NDF),"L");
    leg->AddEntry(fit5,Form("Fit (560 - %4.0f), #chi^{2}/NDF : %2.2f",maxXFit5,fit5FCN/fit5NDF),"L");
    leg->Draw("same");

    //+++++++++++++ First Canvas ends ++++++++++++++++++
    
    
    //+++++++++++++ Second Canvas starts here ++++++++++++++++++

    TCanvas* c2 = new TCanvas("c2","PhotonJet Invariant Mass with Fit",600,650);
    c2->GetWindowHeight();
    c2->GetWindowWidth();
//    c2->SetLogy();
    c2->Divide(1,5,0,0,0);

    //++++++ PAD 1
    c2->cd(1);
    p11_1 = (TPad*)c2->GetPad(1);
    p11_1->SetPad(0.01,0.78,0.99,0.95);
    p11_1->SetRightMargin(0.05);
    p11_1->SetLeftMargin(0.1);
    p11_1->SetTopMargin(0.05);
//    p11_1->SetBottomMargin(0.05);
    p11_1->SetGridx();
    p11_1->SetGridy();
    p11_1->SetTickx(50);

    hRatio1->SetTitle("");
    hRatio1->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hRatio1->GetYaxis()->SetRangeUser(-3.4, 3.3);
    hRatio1->GetYaxis()->SetLabelSize(0.14);
    hRatio1->GetYaxis()->SetLabelOffset(0.016);
    hRatio1->SetYTitle("");
    hRatio1->GetYaxis()->SetTitleSize(0.15);
    hRatio1->GetYaxis()->SetTitleOffset(0.24);

    hRatio1->SetLineWidth(0);
    hRatio1->SetFillColor(2);
    hRatio1->SetLineColor(1);
    hRatio1->Draw("HIST");

    //++++++ PAD 2
    c2->cd(2);
    p11_2 = (TPad*)c2->GetPad(2);
    p11_2->SetPad(0.01,0.61,0.99,0.78);
    p11_2->SetRightMargin(0.05);
    p11_2->SetLeftMargin(0.1);
//    p11_2->SetTopMargin(0.05);
    p11_2->SetGridx();
    p11_2->SetGridy();
    p11_2->SetTickx(50);

    hRatio2->SetTitle("");
    hRatio2->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hRatio2->GetYaxis()->SetRangeUser(-3.4, 3.3);
    hRatio2->GetYaxis()->SetLabelSize(0.14);
    hRatio2->GetYaxis()->SetLabelOffset(0.016);
    hRatio2->SetYTitle("");
    hRatio2->GetYaxis()->SetTitleSize(0.15);
    hRatio2->GetYaxis()->SetTitleOffset(0.24);

    hRatio2->SetLineWidth(0);
    hRatio2->SetFillColor(3);
    hRatio2->SetLineColor(1);
    hRatio2->Draw("HIST");

    //++++++ PAD 3
    c2->cd(3);
    p11_3 = (TPad*)c2->GetPad(3);
    p11_3->SetPad(0.01,0.44,0.99,0.61);
    p11_3->SetRightMargin(0.05);
    p11_3->SetLeftMargin(0.1);
//    p11_3->SetTopMargin(0.05);
    p11_3->SetGridx();
    p11_3->SetGridy();
    p11_3->SetTickx(50);

    hRatio3->SetTitle("");
    hRatio3->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hRatio3->GetYaxis()->SetRangeUser(-3.4, 3.3);
    hRatio3->GetYaxis()->SetLabelSize(0.14);
    hRatio3->GetYaxis()->SetLabelOffset(0.016);
    hRatio3->SetYTitle("(Data-Fit)/Error");
    hRatio3->GetYaxis()->SetTitleSize(0.15);
    hRatio3->GetYaxis()->SetTitleOffset(0.24);

    hRatio3->SetLineWidth(0);
    hRatio3->SetFillColor(4);
    hRatio3->SetLineColor(1);
    hRatio3->Draw("HIST");

    //++++++ PAD 4
    c2->cd(4);
    p11_4 = (TPad*)c2->GetPad(4);
    p11_4->SetPad(0.01,0.27,0.99,0.44);
    p11_4->SetRightMargin(0.05);
    p11_4->SetLeftMargin(0.1);
//    p11_4->SetTopMargin(0.05);
    p11_4->SetGridx();
    p11_4->SetGridy();
    p11_4->SetTickx(50);

    hRatio4->SetTitle("");
    hRatio4->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hRatio4->GetYaxis()->SetRangeUser(-3.4, 3.3);
    hRatio4->GetYaxis()->SetLabelSize(0.14);
    hRatio4->GetYaxis()->SetLabelOffset(0.016);
    hRatio4->SetYTitle("");
    hRatio4->GetYaxis()->SetTitleSize(0.15);
    hRatio4->GetYaxis()->SetTitleOffset(0.24);

    hRatio4->SetLineWidth(0);
    hRatio4->SetFillColor(5);
    hRatio4->SetLineColor(1);
    hRatio4->Draw("HIST");

    //++++++ PAD 5
    c2->cd(5);
    p11_5 = (TPad*)c2->GetPad(5);
    p11_5->SetPad(0.01,0.0,0.99,0.27);
    p11_5->SetRightMargin(0.05);
    p11_5->SetLeftMargin(0.1);
//    p11_5->SetTopMargin(0.05);
    p11_5->SetBottomMargin(0.4);
    p11_5->SetGridx();
    p11_5->SetGridy();
    p11_5->SetTickx(50);

    hRatio5->SetTitle("");
    hRatio5->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hRatio5->GetXaxis()->SetLabelSize(0.14);
    hRatio5->GetXaxis()->SetLabelOffset(0.016);
    hRatio5->GetXaxis()->SetTitleSize(0.15);
    hRatio5->GetXaxis()->SetTitleOffset(1.1);

    hRatio5->GetYaxis()->SetRangeUser(-3.2, 3.2);
    hRatio5->GetYaxis()->SetLabelSize(0.09);
    hRatio5->GetYaxis()->SetLabelOffset(0.016);
    hRatio5->SetYTitle("");
    hRatio5->GetYaxis()->SetTitleSize(0.1);
    hRatio5->GetYaxis()->SetTitleOffset(0.38);

    hRatio5->SetLineWidth(0);
    hRatio5->SetFillColor(6);
    hRatio5->SetLineColor(1);
    hRatio5->Draw("HIST");

}

Double_t fitFunc( Double_t *m, Double_t *p) // function 1
{
      double x=m[0]/13000.;
      return p[0]*pow(1.0-x,p[1])/pow(x,p[2]+p[3]*log(x));
}
Double_t setfitFuncParameters(TF1 *tempfit){
  tempfit->SetParameter(0,1.86812e-10); //-09
  tempfit->SetParameter(1,-1.34760e+01);
  tempfit->SetParameter(2,1.55225e+01);
  tempfit->SetParameter(3,2.11696e+00);
  tempfit->SetLineColor(6);

  return tempfit;
} 
//
//Double_t fitFunc( Double_t *m, Double_t *p)  // function 2
//{
//    double x=m[0]/13000.;
//    return p[0]*pow(1.0+x,p[1])/pow(x,p[2]+p[3]*log(x));
//}
//Double_t setfitFuncParameters(TF1 *tempfit){
//  tempfit->SetParameter(0,1.19582e-13);
//  tempfit->SetParameter(1,2.24040e+01);
//  tempfit->SetParameter(2,1.74540e+01);
//  tempfit->SetParameter(3,2.23427e+00);
//  tempfit->SetLineColor(6);
//
//  return tempfit;
//}
//
//Double_t fitFunc( Double_t *m, Double_t *p)  // function 3
//{
//      double x=m[0]/13000.;
//          return p[0]/pow((p[1]+x*p[2]+x*x),p[3]);
//}
//Double_t setfitFuncParameters(TF1 *tempfit){
//  tempfit->SetParameter(0,1.74081e+04);
//  tempfit->SetParameter(1,7.81222e-01);
//  tempfit->SetParameter(2,8.45393e+00);
//  tempfit->SetParameter(3,7.98450e+00);
//  tempfit->SetLineColor(6);
//
//  return tempfit;
//} 
//
//Double_t fitFunc( Double_t *m, Double_t *p) // function 4
//{
//      double x=m[0]/13000.;
//          return p[0]*pow(1.-x+p[3]*x*x,p[1])/pow(m[0],p[2]);
//}
//Double_t setfitFuncParameters(TF1 *tempfit){
//  tempfit->SetParameter(0,3.09891e+04);
//  tempfit->SetParameter(1,-2.88186e+00);
//  tempfit->SetParameter(2,1.56813e-01);
//  tempfit->SetParameter(3,2.41779e+02);
//  tempfit->SetLineColor(6);
//
//  return tempfit;
//} 
//
//Double_t fitFunc( Double_t *m, Double_t *p)  // function 5
//{
//      double x=m[0]/13000.;
//          return p[0]/pow((p[1]+x),p[2]);
//}
//Double_t setfitFuncParameters(TF1 *tempfit){
//  tempfit->SetParameter(0,4.37740e-04);
//  tempfit->SetParameter(1,9.89145e-02);
//  tempfit->SetParameter(2,8.40025e+00);
//  tempfit->SetLineColor(6);
//
//  return tempfit;
//}
