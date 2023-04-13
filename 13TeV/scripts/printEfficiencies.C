printEfficiencies(){

  const Int_t nEff = 6;
  Double_t inputEff[nEff] = { 0.488178  ,    0.582729 ,    0.598829  ,    0.600764  ,    0.598127  , 0.583467};

  int count_lines=0;
  for( double imass = 1.0; imass <= 5.4; imass += 0.05){
    count_lines++;
    if( count_lines > 10 ){ 
      std::cout<<std::endl;
      count_lines = 1;
    }
    std::cout<<std::setw(8)<< imass << ", ";
  }

  std::cout<<std::endl<<std::endl;


  for( Int_t i = 0; i < nEff-1; ++i){
    Double_t diff = (inputEff[i+1] - inputEff[i])/20.0;

    for(Int_t j = 0; j < 20; ++j){
      if(j == 10)std::cout<<std::endl;
      std::cout<< std::setw(8) << (inputEff[i] + j * diff) <<", ";
    }

    std::cout<<std::endl;
  }
  std::cout<<inputEff[nEff-1]<<std::endl;


}
