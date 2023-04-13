void CompareResonanceShapes(){

    const Int_t nFiles = 8;

    TString FileSignals[nFiles] = { 
	"./rootFiles/QstarToGJ_M_1000.root",
	"./rootFiles/QstarToGJ_M_1500.root",
	"./rootFiles/QstarToGJ_M_2000.root",
	"./rootFiles/QstarToGJ_M_2500.root",
	"./rootFiles/QstarToGJ_M_3000.root",
	"./rootFiles/QstarToGJ_M_3500.root",
	"./rootFiles/QstarToGJ_M_4000.root",
	"./rootFiles/QstarToGJ_M_4500.root"
    };


    TString FileHalfSignals[nFiles] = { 
	"./rootFiles/QstarToGJ_M_fhalf_1000.root",
	"./rootFiles/QstarToGJ_M_fhalf_1500.root",
	"./rootFiles/QstarToGJ_M_fhalf_2000.root",
	"./rootFiles/QstarToGJ_M_fhalf_2500.root",
	"./rootFiles/QstarToGJ_M_fhalf_3000.root",
	"./rootFiles/QstarToGJ_M_fhalf_3500.root",
	"./rootFiles/QstarToGJ_M_fhalf_4000.root",
	"./rootFiles/QstarToGJ_M_fhalf_4500.root"
    };

   Double_t mass[nFiles] = { 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500 };
   TString masses[nFiles] = { "1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500" };

   TString Plotslocation  = "ShapePlots/";
   TString PlotsName  = "ShapesComparison_";
   TString PlotsName1  = "grComparison_";

    TH1F * h_Sig_massVarBin[nFiles];
    TH1F * h_SigHalf_massVarBin[nFiles];
    TPaveText *HistFull_para[nFiles] ;
    TPaveText *HistHalf_para[nFiles] ;

    vector<double> Sig_Mean, Sig_RMS, SigHalf_Mean, SigHalf_RMS ;

    vector<double> qstarfull_X, qstarfull_Y, qstarhalf_X, qstarhalf_Y ;

    for(int i = 0; i < nFiles; ++i){
      TFile* InFull = new TFile(FileSignals[i],"READ");
      TFile* InHalf = new TFile(FileHalfSignals[i],"READ");

      Sig_Mean.clear(); Sig_RMS.clear();  SigHalf_Mean.clear(); SigHalf_RMS.clear();
      qstarfull_X.clear(); qstarfull_Y.clear(); qstarhalf_X.clear(); qstarhalf_Y.clear();
    
      h_Sig_massVarBin[i] = (TH1F*)InFull->Get("h_mass_VarBin_finalPU");
      h_Sig_massVarBin[i]->SetStats(0);
      h_Sig_massVarBin[i]->Scale(1.0/h_Sig_massVarBin[i]->Integral()); 
      h_Sig_massVarBin[i]->GetXaxis()->SetRangeUser(mass[i]-1000,mass[i]+1000);
      h_Sig_massVarBin[i]->GetYaxis()->SetRangeUser(0,0.45);
      
      Sig_Mean.push_back(h_Sig_massVarBin[i]->GetMean());
      Sig_RMS.push_back(h_Sig_massVarBin[i]->GetRMS());

      h_SigHalf_massVarBin[i] = (TH1F*)InHalf->Get("h_mass_VarBin_finalPU");
      h_SigHalf_massVarBin[i]->SetStats(0);
      h_SigHalf_massVarBin[i]->Scale(1.0/h_SigHalf_massVarBin[i]->Integral()); 
      h_SigHalf_massVarBin[i]->GetXaxis()->SetRangeUser(mass[i]-1000,mass[i]+1000);

      SigHalf_Mean.push_back(h_SigHalf_massVarBin[i]->GetMean());
      SigHalf_RMS.push_back(h_SigHalf_massVarBin[i]->GetRMS());

      h_Sig_massVarBin[i]->SetLineColor(kRed);
      h_SigHalf_massVarBin[i]->SetLineColor(kBlue);


     for( int j=0 ; j < h_Sig_massVarBin[i]->GetNbinsX() ; ++j){

	 if(h_Sig_massVarBin[i]->GetBinContent(j+1) != 0  && h_Sig_massVarBin[i]->GetBinCenter(j+1) >= mass[i]-1000 && h_Sig_massVarBin[i]->GetBinCenter(j+1) <= mass[i] +1000){
	     qstarfull_Y.push_back(h_Sig_massVarBin[i]->GetBinContent(j+1));
	     qstarfull_X.push_back(h_Sig_massVarBin[i]->GetBinCenter(j+1));

	     qstarhalf_Y.push_back(h_SigHalf_massVarBin[i]->GetBinContent(j+1));
	     qstarhalf_X.push_back(h_SigHalf_massVarBin[i]->GetBinCenter(j+1));
	 }
     }

     TGraph * gr_qstarfull = new TGraph(qstarfull_X.size() , &qstarfull_X[0], &qstarfull_Y[0]) ;
     gr_qstarfull->SetLineColor(kRed);
     gr_qstarfull->SetLineWidth(2);
     TGraph * gr_qstarhalf = new TGraph(qstarhalf_X.size() , &qstarhalf_X[0], &qstarhalf_Y[0]) ;
     gr_qstarhalf->SetLineColor(kBlue);
     gr_qstarhalf->SetLineWidth(2);


      TCanvas* c1 = new TCanvas("c1","Comparison of Full and Half couplings shapes");
      c1->SetFillColor(0); c1->SetFrameBorderMode(0); c1->SetBorderMode(0);

      HistFull_para[i] = new TPaveText(0.7,0.75,0.9,0.9,"NDC");
      HistFull_para[i]->SetTextColor(kRed);
      HistFull_para[i]->AddText("Coupling, f = 1.0");
      HistFull_para[i]->AddText(Form("Mean : %4.2f",Sig_Mean[i]));
      HistFull_para[i]->AddText(Form("RMS : %4.2f",Sig_RMS[i]));
//      HistFull_para[i]->AddText("");

      HistHalf_para[i] = new TPaveText(0.7,0.6,0.9,0.75,"NDC");
      HistHalf_para[i]->SetTextColor(kBlue);
      HistHalf_para[i]->AddText("Coupling, f = 0.5");
      HistHalf_para[i]->AddText(Form("Mean : %4.2f",SigHalf_Mean[i]));
      HistHalf_para[i]->AddText(Form("RMS : %4.2f",SigHalf_RMS[i]));


      HistFull_para[i]->SetTextFont(42);
      HistFull_para[i]->SetTextAlign(12);
      HistFull_para[i]->SetTextSize(0.04);
      HistFull_para[i]->SetFillColor(0);
      HistFull_para[i]->SetLineColor(0);
      HistFull_para[i]->SetFillStyle(0);
      HistFull_para[i]->SetBorderSize(0);

      HistHalf_para[i]->SetTextFont(42);
      HistHalf_para[i]->SetTextAlign(12);
      HistHalf_para[i]->SetTextSize(0.04);
      HistHalf_para[i]->SetFillColor(0);
      HistHalf_para[i]->SetLineColor(0);
      HistHalf_para[i]->SetFillStyle(0);
      HistHalf_para[i]->SetBorderSize(0);

      h_Sig_massVarBin[i]->Draw("hist");
      h_SigHalf_massVarBin[i]->Draw("histsame");
      HistFull_para[i]->Draw();
      HistHalf_para[i]->Draw();

//      c1->SaveAs(Plotslocation + PlotsName + masses[i] + ".pdf");
//      c1->SaveAs(Plotslocation + PlotsName + masses[i] + ".gif");

//      c1->SaveAs("ShapePlots/" + "ShapesComparison_" + masses[i] + ".pdf");
//      c1->SaveAs("ShapePlots/"+"ShapesComparison_"+masses[i]+".gif");
 //     delete c1;

      
      TCanvas* c2 = new TCanvas("c2","Comparison of Full and Half couplings shapes");
      gr_qstarhalf->Draw("ACP");
      gr_qstarfull->Draw("CP");
      HistFull_para[i]->Draw();
      HistHalf_para[i]->Draw();
//      c2->SaveAs(Plotslocation + PlotsName1 + masses[i] + ".pdf");
//      c2->SaveAs(Plotslocation + PlotsName1 + masses[i] + ".gif");
   //   delete c2;


    }

}
