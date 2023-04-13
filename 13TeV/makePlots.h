#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMinuit.h>
#include <TRandom.h>
#include <string>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <stdio.h>
#include <TString.h>
#include "TH1.h"
#include "TH2.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH1I.h>
#include <TH2I.h>
#include <THStack.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TList.h>
#include <TGraphAsymmErrors.h>
#include <map>
#include "TMath.h"
#include <vector>
#include <TList.h>
#include <Riostream.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TArrow.h>


#ifdef __MAKECINT__                           
#pragma link C++ class map<string,TCanvas*>;
#pragma link C++ class map<string,TPad*>;
#pragma link C++ class map<string,TLegend*>;
#pragma link C++ class map<string,TGraphErrors*>;
#pragma link C++ class map<string,TLine*>;
#pragma link C++ class vector<bool>+;
#pragma link C++ class vector<TFile>+;
#pragma link C++ class vector<Double_t>+;
#pragma link C++ class vector<TString>+;
#pragma link C++ class map<TString,TH1D*>+;
#endif

using namespace std;        
using namespace ROOT;

TString inputFilesPath;
TString prefix;

std::vector < TString > fileNames;  
std::vector < TString > sigName_f1p0;
std::vector < TString > sigName_f0p5;
std::vector < TString > sigName_f0p1;

std::vector < Double_t > scaleFiles;
std::vector < TString > labelFiles;

Bool_t savePlots;
TString plotsLocation;
TString tablesLocation;

Double_t scaleMC;
Double_t scaleGJ;
Double_t scaleQCD;
Double_t scaleEWK;

TLegend *histLegend;

Double_t padLmargin;
Double_t padRmargin;
Double_t padTmargin;
Double_t padBmargin;

void getInitialized(TFile *Inputfiles[], TFile *sigFile_f1p0[], TFile *sigFile_f0p5[], TFile *sigFile_f0p1[]){
  scaleMC  = 1.0;
  scaleGJ  = 0.88;
  scaleQCD = 0.88;
  scaleEWK = 0.88;

  padLmargin = 0.11;
  padRmargin = 0.04;
  padTmargin = 0.06;
  padBmargin = 0.11;



  for(int ifile = 0; ifile < fileNames.size(); ++ifile)
    Inputfiles[ifile] = new TFile(fileNames.at(ifile), "READ");

  for(int ifile = 0; ifile < sigName_f1p0.size(); ++ifile){
    sigFile_f1p0[ifile] = new TFile(sigName_f1p0.at(ifile), "READ");
    sigFile_f0p5[ifile] = new TFile(sigName_f0p5.at(ifile), "READ");
    sigFile_f0p1[ifile] = new TFile(sigName_f0p1.at(ifile), "READ");
  }

  //WARNING::  make sure scaleFiles, labelFiles, and fileNames are all in same order
  scaleFiles.push_back(1.0); // Data
  scaleFiles.push_back(scaleMC*scaleGJ); // GJ
  scaleFiles.push_back(scaleMC*scaleQCD); // QCD
  scaleFiles.push_back(scaleMC*scaleEWK); // EWK
  scaleFiles.push_back(scaleMC); // Signal

  labelFiles.push_back("Data");
  labelFiles.push_back("GJ");
  labelFiles.push_back("DiJet");
  labelFiles.push_back("EWK");
  labelFiles.push_back("Sig");



} // getInitialized

void getInputFileNames(TString inputFilesPath=".", TString prefix=""){

  fileNames.push_back(inputFilesPath+"/Data_"+prefix+".root");
  fileNames.push_back(inputFilesPath+"/GJ_"+prefix+"_JetPhoxCorrected.root");
  //fileNames.push_back(inputFilesPath+"/GJ_"+prefix+".root");
  fileNames.push_back(inputFilesPath+"/DiJet_"+prefix+".root");
  fileNames.push_back(inputFilesPath+"/EWK_"+prefix+".root");
  fileNames.push_back(inputFilesPath+"/QstarToGJ_M2000_f1p0_1.root");

  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M1000_f1p0_1.root");
  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M2000_f1p0_1.root");
  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M3000_f1p0_1.root");
  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M4000_f1p0_1.root");
  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M5000_f1p0_1.root");
  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M6000_f1p0_1.root");
  sigName_f1p0.push_back(inputFilesPath+"/QstarToGJ_M7000_f1p0_1.root");

  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M1000_f0p5_1.root");
  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M2000_f0p5_1.root");
  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M3000_f0p5_1.root");
  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M4000_f0p5_1.root");
  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M5000_f0p5_1.root");
  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M6000_f0p5_1.root");
  sigName_f0p5.push_back(inputFilesPath+"/QstarToGJ_M7000_f0p5_1.root");

  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M1000_f0p1_1.root");
  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M2000_f0p1_1.root");
  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M3000_f0p1_1.root");
  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M4000_f0p1_1.root");
  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M5000_f0p1_1.root");
  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M6000_f0p1_1.root");
  sigName_f0p1.push_back(inputFilesPath+"/QstarToGJ_M7000_f0p1_1.root");
}

void draw1DHistos(TFile *inFile[], TString hist, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax, TString legPos, Bool_t logy, Bool_t b_ratio){

  TH1F * h_temp[fileNames.size()];

  for(Int_t i=0; i<fileNames.size(); ++i)                  
    h_temp[i] = (TH1F*)inFile[i]->Get("h_"+hist);

  h_temp[0]->SetStats(0); //0-Data
  h_temp[0]->SetMarkerStyle(8);
  h_temp[0]->SetMarkerColor(kBlack); 
  h_temp[0]->SetMarkerSize(0.8);    
  h_temp[0]->SetLineColor(kBlack);  
  h_temp[0]->GetXaxis()->SetRangeUser(xMin,xMax);
  h_temp[0]->Scale(scaleFiles[0]);

  h_temp[1]->SetStats(0); //1-GJ
  h_temp[1]->SetLineColor(kOrange+9);
  h_temp[1]->SetFillColor(kOrange+6);
  h_temp[1]->SetLineWidth(2);
  h_temp[1]->GetXaxis()->SetRangeUser(xMin,xMax);
  h_temp[1]->GetXaxis()->SetTitleOffset(1.3);
  h_temp[1]->GetYaxis()->SetTitleOffset(1.4);
  h_temp[1]->Scale(scaleFiles[1]);

  h_temp[2]->SetStats(0); //2-Dijet
  h_temp[2]->SetLineColor(kAzure+4);
  h_temp[2]->SetLineWidth(2);
  h_temp[2]->SetFillColor(kAzure+2);
  h_temp[2]->GetXaxis()->SetRangeUser(xMin,xMax);
  h_temp[2]->Scale(scaleFiles[2]);

  h_temp[3]->SetStats(0); //3-EWk
  h_temp[3]->SetLineColor(kYellow-7);
  h_temp[3]->SetLineWidth(2);
  h_temp[3]->SetFillColor(kYellow-9);
  h_temp[3]->GetXaxis()->SetRangeUser(xMin,xMax);
  h_temp[3]->Scale(scaleFiles[3]);;

  h_temp[4]->SetStats(0); //4-Signal
  h_temp[4]->SetLineColor(kMagenta);
  h_temp[4]->SetLineWidth(2);

  TH1F * h_err_bkg = (TH1F*)h_temp[1]->Clone("h_err_bkg");
  h_err_bkg->Add(h_temp[2]);
  h_err_bkg->Add(h_temp[3]);
  h_err_bkg->SetFillColor(1);
  h_err_bkg->SetFillStyle(3004);

  TH1F * h_Ratio = (TH1F*)h_temp[0]->Clone("h_Ratio");
  h_Ratio->Add(h_err_bkg, -1);
  h_Ratio->Divide(h_temp[0]);
  h_Ratio->SetMarkerColor(kRed+2);
  h_Ratio->SetLineColor(kRed+2);

  h_Ratio->GetYaxis()->SetLabelSize(0.08);
  h_Ratio->GetYaxis()->SetTitle("(Data-MC)/Data");
  h_Ratio->GetYaxis()->SetTitleOffset(.4);
  //  h_Ratio->GetYaxis()->SetTitleSize(0.4);
  h_Ratio->GetYaxis()->SetRangeUser(-1.0,1.0);
  h_Ratio->GetXaxis()->SetLabelSize(0.08);
  h_Ratio->GetXaxis()->SetTitleSize(0.1);

  TLine * ratioLine = new TLine (xMin, 0.0,  xMax, 0.0);
  ratioLine->SetLineStyle(2);
  ratioLine->SetLineColor(38);

  if(legPos == "topRight") histLegend = new TLegend(0.67, 0.58, 0.92, 0.81);
  if(legPos == "topCenter") histLegend = new TLegend(0.5, 0.7, 0.72, 0.92);
  if(legPos == "bottomCenter") histLegend = new TLegend(0.43, 0.17, 0.63, 0.36);
  histLegend->SetBorderSize(0);
  histLegend->SetFillColor(0);
  histLegend->SetFillStyle(0);
  histLegend->SetTextFont(42);
  histLegend->SetTextSize(0.03);
  histLegend->AddEntry(h_temp[0],        "Data", "Lp");
  histLegend->AddEntry(h_temp[1],    "GJ (Bkg)", "F");
  histLegend->AddEntry(h_temp[2], "DiJet (Bkg)", "F");
  histLegend->AddEntry(h_temp[3],   "EWK (Bkg)", "F");
  histLegend->AddEntry(h_temp[4], "q* = 2 TeV (f=1.0)", "L");


  THStack * hS_temp = new THStack("hS_temp","");
  hS_temp->SetMinimum(yMin);
  hS_temp->SetMaximum(yMax);
  hS_temp->Add(h_temp[3]);
  hS_temp->Add(h_temp[2]);
  hS_temp->Add(h_temp[1]);
  hS_temp->SetHistogram(h_temp[1]);

  if( b_ratio == false){
    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    c1->cd();
    c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    c1->SetTitle(0);
    c1->SetRightMargin(padRmargin);
    c1->SetLeftMargin(padLmargin);
    c1->SetTopMargin(padTmargin);
    c1->SetBottomMargin(padBmargin);

    hS_temp->Draw("hist"); // Total Bkg MC
    h_temp[4]->Draw("histSAME"); // Signal
    h_temp[0]->Draw("SAME PE"); // Data

    histLegend->Draw();

    TLatex *cmsText = new TLatex();
    cmsText->SetNDC();
    cmsText->SetTextAngle(0);
    cmsText->SetTextColor(kBlack);
    cmsText->SetTextFont(42);
    cmsText->SetTextAlign(31);
    cmsText->SetTextSize(0.6*padTmargin);
    cmsText->DrawLatex(1 - padRmargin, 1 - padTmargin +  0.2*padTmargin, "2.67 fb^{-1} (13 TeV)");

    cmsText->SetTextFont(61);
    cmsText->SetTextSize(0.65*padTmargin);
    cmsText->SetTextAlign(33); //11-top left;  21-top centre;  31-top right
    cmsText->DrawLatex(0.88, 0.91, "CMS");

    cmsText->SetTextFont(52);
    cmsText->SetTextAlign(33);
    cmsText->SetTextSize(0.76*0.65*padTmargin);
    cmsText->DrawLatex(0.91, 0.91 - 1.2*0.65*padTmargin, "Preliminary");


    if(logy) c1->SetLogy();
    c1->SaveAs(plotsLocation+"/"+hist+".pdf");
    delete c1;
  }// b_ratio-False
  if( b_ratio == true){

    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->Draw(); pad1->cd();
    if(logy) 
      pad1->SetLogy();
    pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
    pad1->SetTitle(0);
    pad1->SetRightMargin(padRmargin);
    pad1->SetLeftMargin(padLmargin);
    pad1->SetTopMargin(padTmargin);
    pad1->SetBottomMargin(0);

    hS_temp->Draw("hist"); // Total Bkg MC
    h_temp[4]->Draw("histSAME"); // Signal
    h_temp[0]->Draw("SAME PE"); // Data

    histLegend->Draw();

    TLatex *cmsText = new TLatex();
    cmsText->SetNDC();
    cmsText->SetTextAngle(0);
    cmsText->SetTextColor(kBlack);
    cmsText->SetTextFont(42);
    cmsText->SetTextAlign(31);
    cmsText->SetTextSize(0.6*padTmargin);
    cmsText->DrawLatex(1 - padRmargin, 1 - padTmargin +  0.2*padTmargin, "2.67 fb^{-1} (13 TeV)");

    cmsText->SetTextFont(61);
    cmsText->SetTextSize(0.65*padTmargin);
    cmsText->SetTextAlign(33); //11-top left;  21-top centre;  31-top right
    cmsText->DrawLatex(0.88, 0.91, "CMS");

    cmsText->SetTextFont(52);
    cmsText->SetTextAlign(33);
    cmsText->SetTextSize(0.76*0.65*padTmargin);
    cmsText->DrawLatex(0.91, 0.91 - 1.2*0.65*padTmargin, "Preliminary");

    c1->cd();

    ///////////////////////////////////////////           
    // Second PAD for making ratio plot 
    TPad *pad2 = new TPad("pad2","",0,0.0,1,0.25);
    pad2->Draw(); pad2->cd();                             
    pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
    pad2->SetTitle(0);
    pad2->SetRightMargin(padRmargin);
    pad2->SetLeftMargin(padLmargin);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetTopMargin(0);                                
    //    pad2->SetBottomMargin(0.3);

    h_Ratio->Draw("EP");
    ratioLine->Draw("SAME");
    h_Ratio->Draw("sameEP");
    c1->cd();

    c1->SaveAs(plotsLocation+"/"+hist+"_ratio.pdf");
    delete c1;

  }// b_ratio-True
} //-draw1DHistos

void draw2DHistos(TFile *inFile, TString hist, Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax, TString histName){

  TH2F * temp_hist = (TH2F*)inFile->Get("h_"+hist);

  temp_hist->SetStats(0);
  temp_hist->GetXaxis()->SetRangeUser(xMin,xMax);
  temp_hist->GetXaxis()->SetTitleOffset(1.1);
  temp_hist->GetXaxis()->SetTitleSize(0.035);
  temp_hist->GetYaxis()->SetRangeUser(yMin,yMax);
  temp_hist->GetYaxis()->SetTitleOffset(1.5);
  temp_hist->GetYaxis()->SetTitleSize(0.035);

  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
  c1->SetLogz(); 
  c1->SetFillColor(0); c1->SetFrameBorderMode(0);
  c1->SetTitle(0);
  c1->SetRightMargin(0.12);
  c1->SetLeftMargin(padLmargin);

  temp_hist->Draw("COLZ");

  c1->SaveAs(plotsLocation+"/"+hist+"_"+histName+".pdf");
  delete temp_hist;
  delete c1;
} //--draw2DHistos

void compareHistos(TString f1, TString var1, TString f2, TString var2, TString outName, Double_t xMin, Double_t xMax, Bool_t logy, Bool_t b_ratio){

  TFile * file1 = new TFile(f1,"READ");
  TH1F * hist1  = (TH1F*)file1->Get(var1);

  TFile * file2 = new TFile(f2,"READ");
  TH1F * hist2  = (TH1F*)file2->Get(var2);

  hist1->SetStats(0);
  hist1->SetLineColor(kRed);
  hist1->SetMarkerColor(kRed);
  hist1->GetXaxis()->SetRangeUser(xMin, xMax);

  hist2->SetStats(0);
  hist2->SetLineColor(kBlue);
  hist2->SetMarkerColor(kBlue);
  hist2->SetMarkerStyle(21);
  hist2->SetMarkerSize(0.5);
  hist2->GetXaxis()->SetRangeUser(xMin, xMax);

  TH1F * h_ratio = (TH1F*)hist1->Clone("h_ratio");
  h_ratio->Add(hist2, -1);
  h_ratio->Divide(hist1);
  h_ratio->SetMarkerStyle(4);
  h_ratio->SetMarkerColor(kRed+2);
  h_ratio->SetLineColor(kRed+2);

  h_ratio->GetYaxis()->SetLabelSize(0.08);
  h_ratio->GetYaxis()->SetTitle("(h1-h2)/h1");
  h_ratio->GetYaxis()->SetTitleOffset(.4);
  h_ratio->GetYaxis()->SetRangeUser(-0.5,0.5);
  h_ratio->GetXaxis()->SetLabelSize(0.08);
  h_ratio->GetXaxis()->SetTitleSize(0.1);

  TLine * ratioLine = new TLine (xMin, 0.0,  xMax, 0.0);
  ratioLine->SetLineStyle(2);
  ratioLine->SetLineColor(38);

  TLegend *leg = new TLegend(0.45,0.8,0.9,0.9);
  leg->SetFillColor(0);  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.025);
  leg->AddEntry(hist1, var1, "F");
  leg->AddEntry(hist2, var2, "LP");

  if(!b_ratio){
    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    c1->cd();
    c1->SetFillColor(0); c1->SetFrameBorderMode(0);
    c1->SetTitle(0);
    c1->SetRightMargin(padRmargin);
    c1->SetLeftMargin(padLmargin);
    c1->SetTopMargin(padTmargin);
    c1->SetBottomMargin(padBmargin);

    hist1->Draw("hist");
    hist2->Draw("SAME PE");
    leg->Draw();

    if(logy)c1->SetLogy();
    c1->SaveAs(plotsLocation+"/"+outName+".pdf");
    delete c1;

  }

  if(b_ratio){

    TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
    TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
    pad1->Draw(); pad1->cd();
    if(logy) 
      pad1->SetLogy();
    pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
    pad1->SetTitle(0);
    pad1->SetRightMargin(padRmargin);
    pad1->SetLeftMargin(padLmargin);
    pad1->SetTopMargin(padTmargin);
    pad1->SetBottomMargin(0);

    hist1->Draw("hist");
    hist2->Draw("SAME PE");
    leg->Draw();
    c1->cd();
    ///////////////////////////////////////////           
    // Second PAD for making ratio plot 
    TPad *pad2 = new TPad("pad2","",0,0.0,1,0.25);
    pad2->Draw(); pad2->cd();                             
    pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
    pad2->SetTitle(0);
    pad2->SetRightMargin(padRmargin);
    pad2->SetLeftMargin(padLmargin);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetTopMargin(0);                                
    //    pad2->SetBottomMargin(0.3);

    h_ratio->Draw("EP");
    ratioLine->Draw("SAME");
    h_ratio->Draw("sameEP");
    c1->cd();

    c1->SaveAs(plotsLocation+"/"+outName+"_ratio.pdf");
    delete c1;
 
  }
}

void compareNhistos(std::vector<TString> infile, std::vector<TString> var, TString outName, std::vector<TString> legName, Double_t xMin, Double_t xMax, Bool_t logy){

  const int nFiles = infile.size();

  TFile * inFile[nFiles];
  TH1F * hist[nFiles];

  padLmargin = 0.11;
  padRmargin = 0.04;
  padTmargin = 0.06;
  padBmargin = 0.11;

  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
  c1->cd();
  c1->SetFillColor(0); c1->SetFrameBorderMode(0);
  c1->SetTitle(0);
  c1->SetRightMargin(padRmargin);
  c1->SetLeftMargin(padLmargin);
  c1->SetTopMargin(padTmargin);
  c1->SetBottomMargin(padBmargin);

  TLegend *leg = new TLegend(0.65,0.8,0.9,0.9);
  leg->SetFillColor(0);  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.025);

  for(int i=0; i<nFiles; ++i){

    inFile[i] = new TFile(infile[i],"READ");
    hist[i]  = (TH1F*)inFile[i]->Get(var[i]);

    hist[i]->SetStats(0);
    hist[i]->SetTitle(0);
    hist[i]->SetLineColor(i+1);
    hist[i]->SetMarkerColor(i+1);

    leg->AddEntry(hist[i], legName[i], "F");
    
    if(i == 0){
      hist[i]->GetYaxis()->SetNdivisions(508);
      hist[i]->GetYaxis()->SetLabelSize(0.035);
      hist[i]->GetYaxis()->SetTitleOffset(1.2);
      hist[i]->GetYaxis()->SetTitleSize(0.035);
      //hist[i]->GetYaxis()->SetRangeUser(-1.,1.);
      hist[i]->GetXaxis()->SetLabelSize(0.035);
      hist[i]->GetXaxis()->SetTitleOffset(1.2);
      hist[i]->GetXaxis()->SetTitleSize(0.035);
      hist[i]->GetXaxis()->SetRangeUser(xMin, xMax);
      hist[i]->Draw("p");
    }else 
      hist[i]->Draw("same p");
  }
  leg->Draw();
    
  c1->cd();
  if(logy)c1->SetLogy();
  c1->SaveAs(plotsLocation+"/"+outName+".pdf");
  delete c1;
}

void getPileUpHist(TFile *inFile[], TString hist){

  TH1F * h_noPU[fileNames.size()-1];
  TH1F * h_PU[fileNames.size()-1];

  for(Int_t ifile=0; ifile<fileNames.size()-1; ++ifile){
    h_noPU[ifile] = (TH1F*)inFile[ifile]->Get("h_"+hist+"_noPU");
    h_PU[ifile] = (TH1F*)inFile[ifile]->Get("h_"+hist+"_PU");
  }

  TH1F * h_tmp_noPU = (TH1F*)h_noPU[1]->Clone("h_tmp_noPU");
  h_tmp_noPU->Scale(scaleFiles[1]);
  h_tmp_noPU->Add(h_noPU[2], scaleFiles[2]);
  h_tmp_noPU->Add(h_noPU[3], scaleFiles[3]);
  h_tmp_noPU->SetLineColor(kBlue);
  h_tmp_noPU->SetLineWidth(2);
  h_tmp_noPU->SetFillColor(kBlue);
  h_tmp_noPU->SetFillStyle(3003);
  h_tmp_noPU->SetStats(0);
  h_tmp_noPU->GetXaxis()->SetRangeUser(0,50);
  h_tmp_noPU->GetXaxis()->SetTitle("Number of Vertices");
 
  h_PU[0]->SetMarkerColor(kBlack);
  h_PU[0]->SetLineColor(kBlack);
  h_PU[0]->SetMarkerStyle(20);

  TH1F * h_tmp_PU = (TH1F*)h_PU[1]->Clone("h_tmp_PU");
  h_tmp_PU->Scale(scaleFiles[1]);
  h_tmp_PU->Add(h_PU[2], scaleFiles[2]);
  h_tmp_PU->Add(h_PU[3], scaleFiles[3]);
  h_tmp_PU->SetFillColor(kYellow);
  
  TLegend *leg = new TLegend(0.65,0.7,0.88,0.88);
  leg->SetFillColor(0); leg->SetTextFont(42); leg->SetBorderSize(0);
  leg->AddEntry(h_PU[0], "Data", "lP");
  leg->AddEntry(h_tmp_noPU, "Total Bkg before re-weighting", "F");
  leg->AddEntry(h_tmp_PU, "Total Bkg after re-weighting", "F");

  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
  c1->SetTitle(0);
  c1->SetFillColor(0);
  c1->SetFrameBorderMode(0) ;

  h_tmp_noPU->DrawNormalized("HIST");
  h_tmp_PU->DrawNormalized("sameHIST");
  h_PU[0]->DrawNormalized("sameP");
  leg->Draw();

  c1->SetLogy();
  c1->SaveAs(plotsLocation+"/"+hist+".pdf");
  delete c1;
}//--getPileUpHist

void getHLTturnON(TFile *fileIn){

  TH1F *h_Num = (TH1F*)fileIn->Get("h_ptTrigPhoton_num");
  TH1F *h_Den = (TH1F*)fileIn->Get("h_ptTrigPhoton_deno");

  TGraphAsymmErrors* triggerTurnOn = new TGraphAsymmErrors();

  TCanvas* c1 = new TCanvas("c1","",80,20,700,700);
  c1->SetFillColor(0);
  c1->SetGridy();
  c1->SetFrameBorderMode(0) ;

  triggerTurnOn->Divide(h_Num,h_Den);                                                                                                                                           
  triggerTurnOn->GetXaxis()->SetRangeUser(100,400);
  triggerTurnOn->GetXaxis()->SetTitle("P_{T}^{#gamma} [GeV]");
  triggerTurnOn->GetYaxis()->SetTitle("Relative Efficiency");
  triggerTurnOn->SetMarkerStyle(8);
  triggerTurnOn->SetMarkerColor(2);
  triggerTurnOn->SetMarkerSize(.7);
  triggerTurnOn->Draw("AP");

  TArrow *ar = new TArrow(190., 0.02, 190., 0.94, 0.03, ">");
  ar->SetLineWidth(2);
  ar->SetLineColor(2);
  ar->Draw();
  TLine *l = new TLine(190., 0.93, 190., 1.1);
  l->SetLineWidth(2);
  l->SetLineColor(2);
  l->SetLineStyle(3);
  l->Draw();

  c1->SaveAs(plotsLocation+"/ptTriggerTurnOn.pdf");
  //delete h_Num; delete h_Den; delete triggerTurnOn; delete c1;
}//--getHLTturnON

void getCutFlowTable(TFile *inFile[]){

  TH1D * h_temp[fileNames.size()];
  for(Int_t i=0; i<fileNames.size(); ++i)                  
    h_temp[i] = (TH1D*)inFile[i]->Get("h_CutExpFlow") ;

  ofstream fout(tablesLocation+"/Efficiencies_"+prefix+".tex");

  std::vector< map< TString, Double_t > > ExpEntries;
  std::vector< TString > Labels;

  fout<<std::setw(12)<<"Cuts"<<"  &  "<<std::setw(7)<<"Signal"<<"  &  "<<setw(12)<<"PhotonJet"<<"  &  "<<setw(12)<<"DiJet"<<"  &  "<<setw(10)<<"EWK"<<"  &  ";
  fout<<setw(12)<<"Total Bkg"<<"  &  "<<setw(12)<<"Data"<<"  &  "<<"$S/\\sqrt{B}$  \\\\"<<std::endl<<"\\hline"<<std::endl;

  for(int ibin=0; ibin < h_temp[0]->GetNbinsX(); ++ibin){
    std::map< TString, Double_t > map_temp;
    ExpEntries.push_back(map_temp);
    Labels.push_back(h_temp[0]->GetXaxis()->GetBinLabel(ibin+1));
    fout<<std::setw(12)<<Labels.at(ibin);

    for(int ifile = 0; ifile < fileNames.size(); ++ifile)
      ExpEntries.at(ibin)[ labelFiles[ifile] ] = h_temp[ifile]->GetBinContent(ibin+1) * scaleFiles[ifile];

    fout<<"  &  "<<std::setw(7)<<ExpEntries.at(ibin)["Sig"]<<"  &  "<<std::setw(12)<<ExpEntries.at(ibin)["GJ"]<<"  &  "<<std::setw(12)<<ExpEntries.at(ibin)["DiJet"];
    fout<<"  &  "<<setw(10)<<ExpEntries.at(ibin)["EWK"];
    fout<<"  &  "<<setw(12)<<ExpEntries.at(ibin)["GJ"] + ExpEntries.at(ibin)["DiJet"] + ExpEntries.at(ibin)["EWK"]<<"  &  "<<setw(12)<<ExpEntries.at(ibin)["Data"];
    fout<<"  &  "<<ExpEntries.at(ibin)["Sig"]/ pow( (ExpEntries.at(ibin)["GJ"] + ExpEntries.at(ibin)["DiJet"] + ExpEntries.at(ibin)["EWK"]) ,0.5)<<"  \\\\"<<std::endl;
  }
  fout.close();
  ExpEntries.clear();
} //--getCutFlowTable

void getSignalEfficiencies(vector<TString>& files, TFile *inFile[], TString outName){

  TH1D * h_temp[files.size()];
  for(Int_t i=0; i<files.size(); ++i)
    h_temp[i] = (TH1D*)inFile[i]->Get("h_CutExpFlow") ;

  Int_t eff = -1;

  std::vector< map< Int_t, Double_t > > ExpEntries;
  std::vector< TString > Labels;

  ofstream fout(tablesLocation+"/SignalEfficiencies_"+outName+"_"+prefix+".tex");

  fout<<"\\hline"<<std::endl<<std::setw(12)<<" Cuts ";
  for(Int_t j = 1; j <= files.size(); ++j)
    fout<<"  &  "<<std::setw(10)<<j*1000;
  fout<<"    \\\\"<<std::endl<<"\\hline\\hline"<<std::endl;

  for(Int_t i = 0; i <  h_temp[0]->GetNbinsX(); ++i){
    std::map< Int_t, Double_t > map_temp; 
    ExpEntries.push_back(map_temp);
    Labels.push_back(h_temp[0]->GetXaxis()->GetBinLabel(i+1));
    fout<<std::setw(12)<<Labels.at(i);

    for(Int_t j = 1; j <= files.size(); ++j){
      ExpEntries.at(i)[j*1000] = h_temp[j-1]->GetBinContent(i+1);
      fout<<"  &  "<<std::setw(10)<<100.0 * (ExpEntries.at(i)[j*1000]/ExpEntries.at(0)[j*1000]);
    }//j-files
    fout<<"    \\\\"<<std::endl;
    if(Labels[i] == "MassCut")  eff=i;
  } //i-nbins

  fout<<"\\hline\\hline"<<std::endl<<std::setw(12)<<"Efficiency";
  for(Int_t j = 1; j <= files.size(); ++j)
    fout<<"  &  "<<std::setw(10)<<ExpEntries.at(eff)[j*1000]/ExpEntries.at(0)[j*1000];
  fout<<"    \\\\"<<std::endl<<"\\hline";

  fout.close();
  ExpEntries.clear();
}//--getSignalEfficiencies
