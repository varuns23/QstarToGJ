#!/bin/bash

cat>printIntermediateEfficiencies.C<<EOF
#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>
 
#include <TGraph.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include <iostream>
#include <stdio.h>

using namespace std;

int main(int argc, char* argv[]){

  if(argc != 4) {
    std::cout << "Usage: printIntermediateEfficiencies <eff-Begin> <eff-End> <nSteps> " << endl;
    return 0;
  } 

  double eff_Begin =  atof(argv[1]);
  double eff_End   =  atof(argv[2]);
  int    nSteps    =  atoi(argv[3]);

  double printIntermediateEfficiencies(double eff_begin, double eff_end, int Nstep);
  printIntermediateEfficiencies(eff_Begin, eff_End, nSteps);

}

double printIntermediateEfficiencies(double eff_begin, double eff_end, int Nstep){

  double step_size = (eff_end - eff_begin)/(double) Nstep;

  for (int iStep = 0; iStep < Nstep; ++iStep)
    std::cout<<eff_begin + iStep*step_size<<",  ";

  std::cout<<std::endl;
}
EOF

g++ -Wno-deprecated printIntermediateEfficiencies.C -o runIntermediateEfficiencies.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`
./runIntermediateEfficiencies.exe ${1} ${2} ${3}

rm runIntermediateEfficiencies.exe
rm printIntermediateEfficiencies.C
