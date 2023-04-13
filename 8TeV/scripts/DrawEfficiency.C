DrawEfficiency(){

    const Int_t n = 6;
//    Float_t Eff[n] = {   54.16,  55.73,  57.53,  57.59,  58.49,  58.69,  58.72,  58.51};
    // Float_t Eff[n] = { 48.8302 ,  58.861 ,  61.5333 ,  62.6118 ,  62.4607 ,  61.2001  }; // LID-f1p0
    // Float_t Eff[n] = { 0.439106 ,  0.535469 ,  0.561357 ,  0.573695 ,  0.575647 ,  0.565831 }; //MID-f1p0
    Float_t Eff[n] = { 0.386917 ,  0.474451 ,  0.498355 ,  0.510284 ,  0.510754 ,  0.498314 , }; //TID-f1p0
  //  Float_t mass[n] = { 1000.0, 1200.0, 1500.0, 1700.0, 2000.0, 2500.0, 3000.0, 3500.0};
    Float_t mass[n] = { 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0};

    const Int_t nhalf = 8;
    Float_t Eff_half[nhalf] = {   54.5,   57.4,   58.7,   58.9,   58.8,   58.7,   58.6,   53.3 };
    Float_t mass_half[nhalf] = {1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0 };
    
    TGraph *grhalfEff = new TGraph(nhalf,mass_half,Eff_half);
    
    TGraph *grEff = new TGraph(n,mass,Eff);

    grhalfEff->SetMarkerStyle(22);
    grEff->SetMarkerStyle(22);


    grEff->Draw("ALP");    // for full coupling
//    grhalfEff->Draw("AP");    // for half coupling

}
