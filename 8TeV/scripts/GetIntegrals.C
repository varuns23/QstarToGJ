#include "vector.h"

void GetIntegrals(){

    TFile *inputData = TFile::Open("./rootFiles/Data_53X8T19pb_Pt170Mass560dEta20TID.root", "READ");
    TH1F *hGJMass_data = (TH1F*) inputData->Get("h_mass_VarBin_finalPU");

    TF1 *fit = new TF1("fit","(0.00274277)*pow((1.0-x/8000.0),3.33455)/pow((x/8000.0),(8.06973 + 0.747310*log(x/8000.0)))");

    Double_t minX = 560.0;
//    Double_t maxX = 3018.0;
    Double_t maxX = 5036.0;

    Double_t integral_data = 0.0;


    for(int i=0;i<hGJMass_data->GetNbinsX();i++){

	if( hGJMass_data->GetBinLowEdge(i+1) >= minX && hGJMass_data->GetBinLowEdge(i+1) < maxX){

//	    std::cout<<i+1<<"   "<<hGJMass_data->GetBinContent(i+1)<<"   "<<hGJMass_data->GetBinWidth(i+1)<<std::endl;

	    integral_data = integral_data + (hGJMass_data->GetBinContent(i+1)*hGJMass_data->GetBinWidth(i+1)) ; 
	}
    }

   Double_t integral_fit =  fit->Integral(minX, maxX);

    std::cout<<"Data Integral : "<<integral_data<<"   Fit Integral : "<<integral_fit<<std::endl;


}



