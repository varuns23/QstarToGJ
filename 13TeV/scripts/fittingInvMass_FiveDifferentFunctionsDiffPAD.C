#include "vector.h"

void fittingInvMass_FiveDifferentFunctionsDiffPAD(){

    // Set TDR Style	
    gROOT->ProcessLine(".L tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    // Input Files  for data
    TFile *inputData = TFile::Open("./rootFiles/Run2015D_PromptReco.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_VarBin_MassCut");

    Double_t minXFit = 560.0 ;
    Double_t maxXFit = 3500.0;

    ofstream fout("ErrorDueToDifferentFit.tex");

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

    vector<double> func_NDF, func_FCN;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 1
    TF1 *fit1 = new TF1("fit1",fitFunc1,minXFit,maxXFit,4); // function 1
    gStyle->SetOptFit(0000000); 
    fit1->SetParameter(0,4.24909e-10); // 1.26914e-05
    fit1->SetParameter(1,-1.52867e+01); // -2.89708e+00
    fit1->SetParameter(2,1.64354e+01); // 1.14282e+01
    fit1->SetParameter(3,2.26986e+00); // 1.32460e+00
    fit1->SetLineWidth(2);
    fit1->SetLineColor(2);

    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit1","","",minXFit,maxXFit);	
    std::cout<<"*********** DONE ******************"<<std::endl;
    //  cout << "Func1 NDF = " << fit1->GetNDF() << "Func1 FCN = " << fit1->GetChisquare() << endl;
    func_NDF.push_back(fit1->GetNDF()) ;
    func_FCN.push_back(fit1->GetChisquare()) ;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 2
    TF1 *fit2 = new TF1("fit2",fitFunc2,minXFit,maxXFit,4); // function 2
    gStyle->SetOptFit(0000000); 
    fit2->SetParameter(0,1.79102e-10); // 3.49595e-02
    fit2->SetParameter(1,2.33652e+01); // -9.28036e+00
    fit2->SetParameter(2,1.62451e+01); // 6.79711e+00
    fit2->SetParameter(3,2.14152e+00); // 5.74076e-01
    fit2->SetLineWidth(2);
    fit2->SetLineColor(3);

    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit2","","",minXFit,maxXFit);	
    std::cout<<"*********** DONE ******************"<<std::endl;

    //    cout << "Func2 NDF = " << fit2->GetNDF() << "Func2 FCN = " << fit2->GetChisquare() << endl;
    func_NDF.push_back(fit2->GetNDF()) ;
    func_FCN.push_back(fit2->GetChisquare()) ;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 3
    TF1 *fit3 = new TF1("fit3",fitFunc3,minXFit,maxXFit,4); 
    gStyle->SetOptFit(0000000); 
    fit3->SetParameter(0,1.74081e+04); // 1.56482e+05
    fit3->SetParameter(1,7.81222e-01); // 6.59912e-01
    fit3->SetParameter(2,8.45393e+00); // 8.46149e+00
    fit3->SetParameter(3,7.98450e+00); // 8.74513e+00
    fit3->SetLineWidth(2);
    fit3->SetLineColor(4);

    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit3","","",minXFit,maxXFit);	
    std::cout<<"*********** DONE ******************"<<std::endl;

    //    cout << "Func3 NDF = " << fit3->GetNDF() << "Func3 FCN = " << fit3->GetChisquare() << endl;
    func_NDF.push_back(fit3->GetNDF()) ;
    func_FCN.push_back(fit3->GetChisquare()) ;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 4
    TF1 *fit4 = new TF1("fit4",fitFunc4,minXFit,maxXFit,4); // function 4
    gStyle->SetOptFit(0000000); 
    fit4->SetParameter(0,3.09891e+04); // 4.39511e+11
    fit4->SetParameter(1,-2.88186e+00); // -2.30338e+00
    fit4->SetParameter(2,1.56813e-01); // 2.51384e+00
    fit4->SetParameter(3,2.41779e+02); // 1.17530e+02

    fit4->SetLineWidth(2);
    fit4->SetLineColor(5);
    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit4","","",minXFit,maxXFit);	
    std::cout<<"*********** DONE ******************"<<std::endl;

    //    cout << "Func4 NDF = " << fit4->GetNDF() << "Func4 FCN = " << fit4->GetChisquare() << endl;
    func_NDF.push_back(fit4->GetNDF()) ;
    func_FCN.push_back(fit4->GetChisquare()) ;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 5
    TF1 *fit5 = new TF1("fit5",fitFunc5,minXFit,maxXFit,3); // function 5
    gStyle->SetOptFit(0000000); 
    fit5->SetParameter(0,4.37805e-04); // 7.33129e-04
    fit5->SetParameter(1,9.89117e-02); // 8.31621e-02
    fit5->SetParameter(2,8.40008e+00); // 9.15528e+00
    fit5->SetLineWidth(2);
    fit5->SetLineColor(6);
    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit5","","",minXFit,maxXFit);	
    std::cout<<"*********** DONE ******************"<<std::endl;

    //    cout << "Func5 NDF = " << fit5->GetNDF() << "Func5 FCN = " << fit5->GetChisquare() << endl;
    func_NDF.push_back(fit5->GetNDF()) ;
    func_FCN.push_back(fit5->GetChisquare()) ;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 6
//    TF1 *fit6 = new TF1("fit6",fitFunc6,minXFit,maxXFit,3); // function 6
//    gStyle->SetOptFit(0000000); 
//    fit6->SetParameter(0,3.98042e+32);
//    fit6->SetParameter(1,6.65334e+02);
//    fit6->SetParameter(2,9.15552e+00);
//    fit6->SetLineWidth(2);
//    fit6->SetLineColor(7);
//    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit6","","",minXFit,maxXFit);	
//    std::cout<<"*********** DONE ******************"<<std::endl;
//
//    //    cout << "Func6 NDF = " << fit6->GetNDF() << "Func6 FCN = " << fit6->GetChisquare() << endl;
//    func_NDF.push_back(fit6->GetNDF()) ;
//    func_FCN.push_back(fit6->GetChisquare()) ;
//
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Fit to data with Function 7
    //    TF1 *fit7 = new TF1("fit7",fitFunc7,minXFit,maxXFit,3); // function 7   Not Converging
    //    gStyle->SetOptFit(0000000); 
    //    fit7->SetParameter(0,8.83045e+05);
    //    fit7->SetParameter(1,4.87766e+01);
    //    fit7->SetParameter(2,1.73136e-01);
    //    fit7->SetLineWidth(2);
    //    fit7->SetLineColor(8);
    //    for(Int_t a = 0 ; a < 10 ; ++a) hGJMass_data->Fit("fit7","","",minXFit,maxXFit);	

    //    cout << "Func7 NDF = " << fit7->GetNDF() << "Func7 FCN = " << fit7->GetChisquare() << endl;
    //    func_NDF.push_back(fit7->GetNDF()) ;
    //    func_FCN.push_back(fit7->GetChisquare()) ;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //    gStyle->SetOptFit(0000);

    float n;
    float n, dm, mass, xl, xh;
    float vx[1200],vy[1200],veyl[1200];
    int i;
    Double_t dataa, errorr;
    vector<double> fitt, ratioo;


fout<<"Mass (in GeV) & Function 1  & Function 2 & Function  & Function 3 & Function 4 & Function  5 \\\\ "<<endl;
    fout<<"700   &  "<<fit1->Eval(700,0,0)<<"  & "<<fit2->Eval(700,0,0)<<"  & "<<fit3->Eval(700,0,0)<<"  & "<<fit4->Eval(700,0,0)<<"  & "<<fit5->Eval(700,0,0)<<" \\"<<"\\"<<endl;
    fout<<"1000  &  "<<fit1->Eval(1000,0,0)<<"  & "<<fit2->Eval(1000,0,0)<<"  & "<<fit3->Eval(1000,0,0)<<"  & "<<fit4->Eval(1000,0,0)<<"  & "<<fit5->Eval(1000,0,0)<<" \\"<<"\\"<<endl;
    fout<<"1500  &  "<<fit1->Eval(1500,0,0)<<"  & "<<fit2->Eval(1500,0,0)<<"  & "<<fit3->Eval(1500,0,0)<<"  & "<<fit4->Eval(1500,0,0)<<"  & "<<fit5->Eval(1500,0,0)<<" \\"<<"\\"<<endl;
    fout<<"2000  &  "<<fit1->Eval(2000,0,0)<<"  & "<<fit2->Eval(2000,0,0)<<"  & "<<fit3->Eval(2000,0,0)<<"  & "<<fit4->Eval(2000,0,0)<<"  & "<<fit5->Eval(2000,0,0)<<" \\"<<"\\"<<endl;
    fout<<"2500  &  "<<fit1->Eval(2500,0,0)<<"  & "<<fit2->Eval(2500,0,0)<<"  & "<<fit3->Eval(2500,0,0)<<"  & "<<fit4->Eval(2500,0,0)<<"  & "<<fit5->Eval(2500,0,0)<<" \\"<<"\\"<<endl;
    fout<<"3000  &  "<<fit1->Eval(3000,0,0)<<"  & "<<fit2->Eval(3000,0,0)<<"  & "<<fit3->Eval(3000,0,0)<<"  & "<<fit4->Eval(3000,0,0)<<"  & "<<fit5->Eval(3000,0,0)<<" \\"<<"\\"<<endl;
fout.close();

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


	    fitt.push_back(fit1->Eval(vx[i],0,0));
	    fitt.push_back(fit2->Eval(vx[i],0,0));
	    fitt.push_back(fit3->Eval(vx[i],0,0));
	    fitt.push_back(fit4->Eval(vx[i],0,0));
	    fitt.push_back(fit5->Eval(vx[i],0,0));
	    
	    if(vx[i] >= minXFit){

		ratioo.push_back((dataa-fitt[0])/errorr);
		ratioo.push_back((dataa-fitt[1])/errorr);
		ratioo.push_back((dataa-fitt[2])/errorr);
		ratioo.push_back((dataa-fitt[3])/errorr);
		ratioo.push_back((dataa-fitt[4])/errorr);

	    }

	    hRatio1->SetBinContent(i+1,ratioo[0]);
	    hRatio2->SetBinContent(i+1,ratioo[1]);
	    hRatio3->SetBinContent(i+1,ratioo[2]);
	    hRatio4->SetBinContent(i+1,ratioo[3]);
	    hRatio5->SetBinContent(i+1,ratioo[4]);

	    fitt.clear();
	    ratioo.clear();
	}

    }          

    //Dijet Mass Cross Section with Fit	
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
    fit1->Draw("SAME");
    fit2->Draw("same");
    fit3->Draw("same");
    fit4->Draw("same");
    fit5->Draw("same");
    //    fit6->Draw("same");
    //    fit7->Draw("same");

    TLegend *leg = new TLegend(0.35,0.77,0.9,0.92);
    leg->SetTextSize(0.03);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->AddEntry(fit1,Form("Func1 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[0],func_NDF[0]),"L");
    leg->AddEntry(fit2,Form("Func2 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[1],func_NDF[1]),"L");
    leg->AddEntry(fit3,Form("Func3 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[2],func_NDF[2]),"L");
//    leg->AddEntry(fit4,Form("Func4 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[3],func_NDF[3]),"L");
    leg->AddEntry(fit4,"Func4 : (#chi2/NDF=20.84/34.00)","L");
    leg->AddEntry(fit5,Form("Func5 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[4],func_NDF[4]),"L");
    //    leg->AddEntry(fit6,Form("Func6 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[5],func_NDF[5]),"L");
    //    leg->AddEntry(fit7,Form("Func7 : (#chi2/NDF=%0.2f/%0.2f)",func_FCN[6],func_NDF[6]),"L");
    leg->Draw("same");
    hGJMass_data->Draw("sameP");
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
    p11_1->SetPad(0.01,0.77,0.99,0.95);
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
    p11_2->SetPad(0.01,0.59,0.99,0.77);
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
    p11_3->SetPad(0.01,0.41,0.99,0.59);
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
    p11_4->SetPad(0.01,0.23,0.99,0.41);
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
    p11_5->SetPad(0.01,0.0,0.99,0.23);
    p11_5->SetRightMargin(0.05);
    p11_5->SetLeftMargin(0.1);
    p11_5->SetBottomMargin(0.05);
//    p11_5->SetTopMargin(0.05);
    p11_5->SetGridx();
    p11_5->SetGridy();
    p11_5->SetTickx(50);

    hRatio5->SetTitle("");
    hRatio5->GetXaxis()->SetRangeUser(minXFit,maxXFit);
    hRatio5->GetYaxis()->SetRangeUser(-3.4, 3.3);
    hRatio5->GetYaxis()->SetLabelSize(0.14);
    hRatio5->GetYaxis()->SetLabelOffset(0.016);
    hRatio5->SetYTitle("");
    hRatio5->GetYaxis()->SetTitleSize(0.15);
    hRatio5->GetYaxis()->SetTitleOffset(0.24);

    hRatio5->SetLineWidth(0);
    hRatio5->SetFillColor(6);
    hRatio5->SetLineColor(1);
    hRatio5->Draw("HIST");

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

Double_t fitFunc4( Double_t *m, Double_t *p) // function 4
{
    double x=m[0]/13000.;
    return p[0]*pow(1.-x+p[3]*x*x,p[1])/pow(m[0],p[2]);
}

Double_t fitFunc5( Double_t *m, Double_t *p)  // function 5
{
    double x=m[0]/13000.;
    return p[0]/pow((p[1]+x),p[2]);
}

Double_t fitFunc6( Double_t *m, Double_t *p)
{
    return p[0]/pow(m[0]+p[1],p[2]);
}

// QCD fit function -- alternate 3 parameter fit function -- also used for QCD fit.
Double_t fitFunc7( Double_t *m, Double_t *p)
{
    double x=m[0]/13000.;
    return p[0]*pow(1.-x,p[1])/pow(m[0],p[2]);
}

