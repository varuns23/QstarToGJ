#!/usr/bin/python

from ROOT import *
from array import array
import CMS_lumi, tdrstyle
import subprocess
import os
import imp
import multiprocessing
from itertools import repeat
import math

gStyle.SetOptFit(1111) 
#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_8TeV = "18.3 fb^{-1}"
CMS_lumi.lumi_13TeV = "803 pb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4
 
  
#fileNameSuffix = "test_range_1118_3704"
#fileNameSuffix = "test"
#fileNameSuffix = "JEC_L2L3Residuals"
#fileNameSuffix = "JEC_V4_firstbin1181"
fileNameSuffix = "DCSonly_JEC_Summer15_25nsV5_803pb-1"


#Fit functions
# 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" 
#
#1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#
#2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#   --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
#
#3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#   --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
#
#4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )" 
#   --> "log" extension wrt to DEFAULT     
#
#5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )" 
#   --> "log" extension wrt to DEFAULT     
#
#6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )" 
#   --> "log" extension wrt to DEFAULT     
#
number_of_variableWidth_bins = 103

massBins =[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430,10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000];

v_massBins = array("d",massBins)
lumi = 803.
#sf = 1.07
sf=1.
#minX_mass = 1000.
minX_mass = 1181.
#maxX_mass = 3019.
#maxX_mass = 5253.
#maxX_mass = 5663.
maxX_mass = 6564.
#================================================================================================================
  
def main():
  
  # data 
  #input_root_file = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_finalJSON_JEC_Summer15_50nsV2_L2L3Residuals/histo_data_mjj_fromTree.root"
  #input_root_file = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_finalJSON_25_07_15_JEC_Summer15_50nsV2/histo_data_mjj_fromTree.root"
  #input_root_file = "plots_data4T_Run2015B_plus_Run2015C_50ns_Cert_json_29Aug2015_xsecSpring15_withSF/histo_data_mjj_fromTree.root"
  #input_root_file = "plots_data4T_Run2015D_DCSonly_390pb-1_JEC_Summer15_25nsV3_withSF/histo_data_mjj_fromTree.root"
  input_root_file = "plots_data4T_Run2015D_DCSonly_JEC_Summer15_50nsV5_withSF/histo_data_mjj_fromTree.root"
  ###mc
  #input_root_file_mc = "/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_finalJSON_25_07_15_JEC_Summer15_50nsV2/histo_data_mjj_fromTree.root"
  #input_root_file_mc = "plots_data4T_Run2015B_plus_Run2015C_50ns_Cert_json_29Aug2015_xsecSpring15_withSF/histo_data_mjj_fromTree.root"
  input_root_file_mc = "plots_data4T_Run2015D_DCSonly_JEC_Summer15_50nsV5_withSF/histo_data_mjj_fromTree.root"
  ### input file and 1D histo
  
  file0 = TFile.Open( input_root_file )
  fileMC = TFile.Open( input_root_file_mc )
  input_1Dhistogram = "h_dat"
  input_1Dhistogram_mc = "hist_allCutsQCD"
  
  hist_mass_original = file0.Get(input_1Dhistogram)
  hist_binned = hist_mass_original.Rebin(number_of_variableWidth_bins,"hist_binned",v_massBins)
  hist_mass = TH1F("hist_mass","",number_of_variableWidth_bins,v_massBins)
  hist_mass_original.Scale(1/lumi)

  hist_mass_original_mc = fileMC.Get(input_1Dhistogram_mc)
  hist_mass_original_mc.Scale(sf)
  hist_binned_mc = hist_mass_original_mc.Rebin(number_of_variableWidth_bins,"hist_binned_MC",v_massBins)
  hist_mass_mc = TH1F("hist_mass_mc","",number_of_variableWidth_bins,v_massBins)
  
  for  i in range (1, number_of_variableWidth_bins):
    #data
    bincontent = hist_binned.GetBinContent(i)
    binwidth = hist_binned.GetBinWidth(i)
    binerror = hist_binned.GetBinError(i)
    hist_mass.SetBinContent(i,bincontent/(binwidth*lumi))   
    #mc
    bincontent_mc = hist_binned_mc.GetBinContent(i)
    binwidth_mc = hist_binned_mc.GetBinWidth(i)
    hist_mass_mc.SetBinContent(i,bincontent_mc/(binwidth_mc*lumi))   

  #hist_mass.Draw()
  #filetest = TFile("filetest.root","recreate")
  #filetest.cd()
  #hist_mass.Write()
  #filetest.Close()
  
  #######################################################
  #data in TGraph format (hist binned)
  alpha = 1 - 0.6827;
  
  x=[]
  y=[]
  exl=[]
  exh=[]
  eyl=[]
  eyh=[]
  x_mc=[]
  y_mc=[]
  
  for i in range(0,number_of_variableWidth_bins):
    n    = hist_binned.GetBinContent(i+1)
    dm   = hist_binned.GetBinWidth(i+1)
    mass = hist_binned.GetBinCenter(i+1)
    xl   = hist_binned.GetBinLowEdge(i+1)
    xh   = xl+dm
    x.append( (xl+xh)/2.)
    exl.append( dm/2.)
    exh.append( dm/2.)
    y.append( n / (dm*lumi))
    l = 0.5*TMath.ChisquareQuantile(alpha/2,2*n)
    h = 0.5*TMath.ChisquareQuantile(1-alpha/2,2*(n+1))
    eyl.append( (n-l)/(lumi*dm) )
    eyh.append( (h-n)/(lumi*dm) )
    #print "%f   %f    %f    %f    %f     %f" % (x[i],y[i],exl[i],exh[i],eyl[i],eyh[i])
    n_mc = hist_binned_mc.GetBinContent(i+1)
    if (i>=40 and i<103):
      x_mc.append((xl+xh)/2.)
      y_mc.append( n_mc / (dm*lumi)) 

  vx = array("f",x)
  vy = array("f",y)
  vexl = array("f",exl)
  vexh = array("f",exh)
  veyl = array("f",eyl)
  veyh = array("f",eyh)
  vx_mc = array("f",x_mc)
  vy_mc = array("f",y_mc)

  #data in TGraph format
  g = TGraphAsymmErrors(number_of_variableWidth_bins,vx,vy,vexl,vexh,veyl,veyh)
  g.SetName("g_data")
  g.Print()
  #mc
  g_mc = TGraph(number_of_variableWidth_bins-40,vx_mc,vy_mc)
  g_mc.SetName("g_mc")
  g_mc.Print()
  
  
  nBins_fit = hist_mass.FindBin(maxX_mass)- hist_mass.FindBin(minX_mass) 
  ##count zero bins
  zeroBins = 0
  for i in range(0,nBins_fit):
    if hist_mass.GetBinContent(hist_mass.FindBin(minX_mass)+i)==0: 
      zeroBins +=1
  
  nBins_fit = nBins_fit-zeroBins   
  FunctionTypes = [-2,-1,0,4,5,6]
  #FunctionTypes = [-1,0]#,4,5,6]
  list_RSS = []
  list_chi2 = []
  list_dof = []
  list_F = []
  list_CL = []
  list_pvalue_WaldTest = []
  list_CL_WaldTest = []
  i_f = 0
  for FunctionType in FunctionTypes:
    fitresult = doFitAndChi2(FunctionType,hist_mass,g,hist_mass_original)
    list_RSS.append(fitresult[0])
    list_chi2.append(fitresult[2])
    list_dof.append(fitresult[1])
    M1Bkg = fitresult[3]
    hist_fit_residual_vsMass = fitresult[4]
    nPar = nBins_fit - fitresult[1]# - 1
    result_WaldTest = WaldWolfowitzTest(hist_fit_residual_vsMass)
    DrawFit(g,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,fileNameSuffix)
    print "chi2 / dof for f%d = %f / %d" % (FunctionType,fitresult[2],fitresult[1])
    list_pvalue_WaldTest.append(result_WaldTest)
    list_CL_WaldTest.append(result_WaldTest)
    if (i_f > 0):
      result = FisherTest(list_RSS[i_f-1],list_RSS[i_f],list_dof[i_f-1],list_dof[i_f],nBins_fit)
      F = result[0]
      CL = result[1]
      list_F.append(F)
      list_CL.append(CL)
    i_f += 1
  
  result_mc = doChi2MC(g,hist_mass,hist_mass_mc)  
  chi2_mc = result_mc[0] 
  Ndof_mc = result_mc[1] 
  hist_mc_residual_vsMass = result_mc[2] 
  pvalue_WaldTest_mc = WaldWolfowitzTest(hist_mc_residual_vsMass)
  CL_WaldTest_mc =1- pvalue_WaldTest_mc
  DrawMC(g,g_mc,hist_mass,hist_mc_residual_vsMass,fileNameSuffix)

  print "nBins_Fit = "+str(nBins_fit)
  print "list_RSS = "+str(list_RSS)
  print "chi2 MC = "+str(chi2_mc)
  print "pvalue_WaldTest MC = "+str(pvalue_WaldTest_mc)
  print "CL_WaldTest MC = "+str(CL_WaldTest_mc)
  print "list_chi2 = "+str(list_chi2)
  print "list_dof ="+str(list_dof)
  print "list_F = "+str(list_F)
  print "list_CL = "+str(list_CL)
  print "list_pvalue_WaldTest = "+str(list_pvalue_WaldTest)
  print "list_CL_WaldTest = "+str(list_CL_WaldTest)

def doChi2MC(g,hist_mass,hist_mass_mc):
  hist_mc_residual_vsMass =  TH1D("hist_mc_residual_vsMass","hist_mc_residual_vsMass",number_of_variableWidth_bins,v_massBins)
  NumberOfObservations_VarBin = 0
  NumberOfVarBins = 0
  chi2_VarBin_zeroes = 0
  chi2_VarBin = 0
  for bin in range (1,number_of_variableWidth_bins):
    mc_residual = 0
    hist_mc_residual_vsMass.SetBinContent(bin,mc_residual)
  
    if( hist_mass_mc.GetXaxis().GetBinLowEdge(bin)>=minX_mass and hist_mass_mc.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
      NumberOfVarBins += 1
      #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
      data = hist_mass.GetBinContent(bin)
      err_data_low = g.GetErrorYlow(bin-1) 
      err_data_high= g.GetErrorYhigh(bin-1)
      mc = hist_mass_mc.GetBinContent(bin)
      #print "mc = %f" % mc 
      if(mc > data or mc == 0): err_tot = err_data_high
      else: err_tot = err_data_low
      # if err_tot==0: 
      #   print "!!!! ERR = 0 !!!"
      #   print "mc %f  data %f err_data_low %f    err_data_high %f " % (mc, data,err_data_low,err_data_high)
      mc_residual = (data - mc) / err_tot
      err_mc_residual = 1
      print "data %f   mc %f   error %f  residual %f" % (data,mc,err_tot,mc_residual)
      chi2_VarBin_zeroes += pow( (data - mc) , 2 ) / pow( err_tot , 2 )

      ##skip bin with zero entries
      if (hist_mass.GetBinContent(bin)>0): 
	NumberOfObservations_VarBin+=1
        chi2_VarBin += pow( (data - mc) , 2 ) / pow( err_tot , 2 )	 
      hist_mc_residual_vsMass.SetBinContent(bin,mc_residual)
    print "bin : %d   mc_residual : %f" % (bin,mc_residual) 
    
  ndf_VarBin = NumberOfObservations_VarBin #-1 
  ndf_VarBin_withzeroes = NumberOfVarBins #-1 
  print "============ MC ==============" 
  print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
  print "ndf_VarBin: %d" % ndf_VarBin 
  print "chi2_VarBin: %f" % chi2_VarBin
  print "ndf_VarBin with zeroes: %d" % ndf_VarBin_withzeroes 
  print "chi2_VarBin with zeroes: %f" % chi2_VarBin_zeroes
  print "============================"   
  return [chi2_VarBin,ndf_VarBin,hist_mc_residual_vsMass]
  

def doFitAndChi2(FunctionType,hist_mass,g,hist_mass_original):
  ### fit mass histogram with background function
  # -2: VARIATION-1 (3 par.) - " [0] / ( TMath::Power(x/13000,[2]) )" 
  if( FunctionType==-2 ):    
    nPar=2
    M1Bkg = TF1("M1Bkg"," [0] / ( TMath::Power(x/13000,[1]) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,2)
    M1Bkg.SetParLimits(0,0.,10);
    M1Bkg.SetParLimits(1,0.,20);
  
  
  # -1: VARIATION-1 (3 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]) )" 
  if( FunctionType==-1 ):    
    nPar=3
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.00001)
    M1Bkg.SetParameter(1,10)
    M1Bkg.SetParameter(2,5)
  
    #1 fb-1
    M1Bkg.SetParLimits(0,0.,1000);
    #10 fb-1
    #M1Bkg->SetParLimits(0,0,1.);
    #M1Bkg.SetParLimits(1,0.,100)
    #M1Bkg.SetParLimits(2,0.,20)
  
  # 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" 
  if( FunctionType==0 ):    
    nPar=4
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.0014)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
  
    M1Bkg.SetParLimits(0,0.,1000);
    #1 fb-1
    #M1Bkg.SetParLimits(0,0.,10);
    #10 fb-1
    #M1Bkg->SetParLimits(0,0,1.);
   # M1Bkg.SetParLimits(1,0.,10000)
   # M1Bkg.SetParLimits(2,1.,10000)
   # M1Bkg.SetParLimits(1,0,100)
   # M1Bkg.SetParLimits(2,0.,20)
   # M1Bkg.SetParLimits(3,-10.,10)
   # M1Bkg.FixParameter(3,0.)
     
  
  # 1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  if( FunctionType==1 ):    
    nPar=5
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.005)
    M1Bkg.SetParameter(1,9.3)
    M1Bkg.SetParameter(2,7.2)
    M1Bkg.SetParameter(3,0.4)
    M1Bkg.SetParameter(4,3.1)
  #   M1Bkg.SetParLimits(4,-5,5)
  
  # 2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    #  --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
  if( FunctionType==2 ):    
    nPar=6
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass);
    M1Bkg.SetParameter(0,0.005)
    M1Bkg.SetParameter(1,9.3)
    M1Bkg.SetParameter(2,7.2)
    M1Bkg.SetParameter(3,0.4)
    M1Bkg.SetParameter(4,3.1)
    M1Bkg.SetParameter(5,25.6)
    #M1Bkg.SetParLimits(5,10,50)      
      
  
  # 3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
  #    --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
  if( FunctionType==3 ) :   
    nPar=7
    M1Bkg =  TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" ,minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.005)
    M1Bkg.SetParameter(1,15.1)
    M1Bkg.SetParameter(2,7.2)
    M1Bkg.SetParameter(3,0.4)
    M1Bkg.SetParameter(4,13.0)
    M1Bkg.SetParameter(5,-4.0)
    M1Bkg.SetParameter(6,70.0)
    # M1Bkg.SetParLimits(4,-1,1)
  
  
  # 4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )" 
  #    --> "log" extension wrt to DEFAULT     
  if( FunctionType==4 ) :   
    nPar=5
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )",minX_mass,maxX_mass);
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
    M1Bkg.SetParameter(4,0.)
    #M1Bkg.SetParLimits(3,0,0.4)
    #M1Bkg.FixParameter(3,0.)
    #M1Bkg.FixParameter(4,0.)
    #M1Bkg.SetParLimits(1,0.,10000)
    #M1Bkg.SetParLimits(2,1.,10000)
    #M1Bkg.SetParLimits(0,0.,10);
    #M1Bkg.SetParLimits(1,0.,100)
    #M1Bkg.SetParLimits(2,0.,20)
    #M1Bkg.SetParLimits(3,-10,10)
    #M1Bkg.SetParLimits(4,-10,10)
   
  
  # 5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )" 
  #    --> "log" extension wrt to DEFAULT     
  if( FunctionType==5 ):    
    nPar=6
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
    M1Bkg.SetParameter(4,0.)
    M1Bkg.SetParameter(5,0.)
    #M1Bkg.SetParLimits(3,0,0.4)
    #M1Bkg.FixParameter(3,0.)
    #M1Bkg.FixParameter(4,0.)
    #M1Bkg.FixParameter(5,0.)
    #M1Bkg.SetParLimits(1,0.,10000)
    #M1Bkg.SetParLimits(2,1.,10000)
    #M1Bkg.SetParLimits(3,0.,10000)
    #M1Bkg.SetParLimits(0,0.,10);
    #M1Bkg.SetParLimits(1,1.,100)
    #M1Bkg.SetParLimits(2,1.,20)
    #M1Bkg.SetParLimits(3,-10,10)
    #M1Bkg.SetParLimits(4,-10,10)
    #M1Bkg.SetParLimits(5,-10,10)
   
  
  # 6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )" 
  #  //    --> "log" extension wrt to DEFAULT     
  if( FunctionType==6 ) :   
    nPar=7
    M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )",minX_mass,maxX_mass)
    M1Bkg.SetParameter(0,0.08)
    M1Bkg.SetParameter(1,12)
    M1Bkg.SetParameter(2,2)
    M1Bkg.SetParameter(3,-0.5)
    M1Bkg.SetParameter(4,0.)
    M1Bkg.SetParameter(5,0.)
    M1Bkg.SetParameter(6,0.)
    #M1Bkg.SetParLimits(3,0,0.4)
    #M1Bkg.FixParameter(3,0.)
    #M1Bkg.FixParameter(4,0.)
    #M1Bkg.FixParameter(5,0.)
    #M1Bkg.FixParameter(6,0.)
    #M1Bkg.SetParLimits(1,0.,10000)
    #M1Bkg.SetParLimits(2,1.,10000)
    #M1Bkg.SetParLimits(3,0.,10000)
    #M1Bkg.SetParLimits(0,0.,10);
    #M1Bkg.SetParLimits(1,0.,100)
    #M1Bkg.SetParLimits(2,0.,20)
    #M1Bkg.SetParLimits(3,-10,10)
    #M1Bkg.SetParLimits(4,-10,10)
    #M1Bkg.SetParLimits(5,-10,10)
    #M1Bkg.SetParLimits(6,-10,10)
  
  
  
  #TFitResultPtr r;
  stopProgram=1;
  for loop in range (0,10):
    r = hist_mass_original.Fit("M1Bkg","ELSR","",minX_mass,maxX_mass)      
    #r = hist_mass_original.Fit("M1Bkg","MSR","",minX_mass,maxX_mass)      
    fitStatus = int(r)
    print "fit status : %d" % fitStatus
    if(fitStatus==0):
      stopProgram=0
      r.Print("V")  
      break
   
  
  if(stopProgram==1):
    print "######################" 
    print"######################" 
    print "ERROR : Fit failed!!!!" 
    print "######################" 
    print "######################" 
  
  
  # fit residuals and chi2
  hist_fit_residual_vsMass =  TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,v_massBins)
  hist_fit_residual = TH1D("hist_fit_residual","hist_fit_residual",10,-5,5)
  NumberOfVarBins = 0
  NumberOfObservations_VarBin = 0
  chi2_VarBin = 0.
  chi2_VarBin_notNorm = 0.
  chi2_VarBin_zeroes = 0.

  for bin in range (1,number_of_variableWidth_bins):
  
    if( hist_mass.GetXaxis().GetBinLowEdge(bin)>=minX_mass and hist_mass.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
      NumberOfVarBins += 1
      #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
      data = hist_mass.GetBinContent(bin)
      err_data_low = g.GetErrorYlow(bin-1) 
      err_data_high= g.GetErrorYhigh(bin-1)
      fit = M1Bkg.Integral(hist_mass.GetXaxis().GetBinLowEdge(bin) , hist_mass.GetXaxis().GetBinUpEdge(bin) )
      fit = fit / ( hist_mass.GetBinWidth(bin) )
      if(fit > data): err_tot = err_data_high
      else: err_tot = err_data_low
      fit_residual = (data - fit) / err_tot
      err_fit_residual = 1
      ##skip bin with zero entries
      chi2_VarBin_zeroes += pow( (data - fit) , 2 ) / pow( err_tot , 2 )
      if (hist_mass.GetBinContent(bin)>0): 
	NumberOfObservations_VarBin+=1
        chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_tot , 2 )	 
        chi2_VarBin_notNorm += pow( (data - fit) , 2 ) 	 
  
  
      hist_fit_residual_vsMass.SetBinContent(bin,fit_residual)
      hist_fit_residual_vsMass.SetBinError(bin,err_fit_residual)
      hist_fit_residual.Fill(fit_residual)
    
  expected = M1Bkg.Integral(5058,10072) * lumi
  observed = hist_binned.Integral(hist_mass.FindBin(5058),hist_mass.FindBin(10072)) 
  ndf_VarBin = NumberOfObservations_VarBin - nPar# -1
  ndf_VarBin_withzeroes = NumberOfVarBins - nPar# -1
  print "============================" 
  print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
  print "ndf_VarBin: %d" % ndf_VarBin 
  print "ndf_VarBin with zeroes: %d" % ndf_VarBin_withzeroes 
  print "chi2_VarBin with zeroes: %f" % chi2_VarBin_zeroes
  print "chi2_VarBin: %f" % chi2_VarBin
  print "chi2_VarBin_notNorm: %f" % chi2_VarBin_notNorm
  print "expected events > 5 TeV (fit) : %f" % expected
  print "observed events > 5 TeV  : %f" % observed
  print "============================"   
  return [chi2_VarBin_notNorm,ndf_VarBin,chi2_VarBin,M1Bkg,hist_fit_residual_vsMass]


def FisherTest(RSS_1,RSS_2,dof_1,dof_2,N):
  RSS1 = RSS_1
  RSS2 = RSS_2
  n1 = N - dof_1 - 1
  n2 = N - dof_2 - 1
  print "n1 = %d    n2 = %d" % (n1,n2)
  F = ((RSS1-RSS2)/(n2-n1)) / (RSS2/(N-n2))
  #print "F = %f" % F
  F_dist = TF1("F_distr","TMath::Sqrt( (TMath::Power([0]*x,[0]) * TMath::Power([1],[1])) / (TMath::Power([0]*x + [1],[0]+[1])) ) / (x*TMath::Beta([0]/2,[1]/2))",0,1000)
  print "d1 = %d    d2 = %d" %(n2-n1,N-n2)
  F_dist.SetParameter(0, n2-n1)
  F_dist.SetParameter(1, N-n2)
  CL = 1 - F_distr.Integral(0.00000001,F)
  #c2 = TCanvas("c2","c2",600,600)
  #c2.cd()
  #c2.DrawFrame(0, 0, 5, 10, "Global Title;X Axis Title;Y Axis Title")
  #F_dist.Draw("same")
  return [F,CL]


def WaldWolfowitzTest(hist_fit_residual_vsMass):
  Nruns = 1
  Nplus = 0
  Nminus = 0
  N = hist_fit_residual_vsMass.FindBin(maxX_mass) - hist_fit_residual_vsMass.FindBin(minX_mass)
  for bin in range(hist_fit_residual_vsMass.FindBin(minX_mass)+1,hist_fit_residual_vsMass.FindBin(maxX_mass)):
    bincontent = hist_fit_residual_vsMass.GetBinContent(bin)
    previousbincontent = hist_fit_residual_vsMass.GetBinContent(bin-1)
    if (previousbincontent > 0): Nplus+= 1 
    if (previousbincontent < 0 ):  Nminus+= 1
    if( bincontent*previousbincontent < 0): Nruns += 1
  if (bincontent > 0): Nplus+= 1
  if (bincontent < 0 ):  Nminus+= 1

  print "N %d Nruns %d   Nplus %d   Nminus %d " %(N,Nruns,Nplus,Nminus)
  Pdf  = TF1("WaldWolfowitzProb","[0]*exp(-0.5*(x-[1])**2/[2]**2)",-10000,10000) 
  mu = float(2*Nplus*Nminus)/float(N) + 1
  sigma2 = float((mu-1)*(mu-2))/float(N-1)
  sigma = TMath.Sqrt(sigma2)
  norm = 1/(TMath.Sqrt(2*TMath.Pi())*sigma)
  print "mu %f   sigma %f" %(mu,sigma)
  Pdf.SetParameter(0, norm)
  Pdf.SetParameter(1,mu)
  Pdf.SetParameter(2,sigma)
  pvalue = Pdf.Integral(-10000.,Nruns)
  print "pvalue %f" %(pvalue)

  return pvalue

def DrawMC(g,g_mc,hist_mass,hist_mc_residual_vsMass,fileNameSuffix):
#  //### Draw plots
  W = 600
  H = 650
  H_ref = 650 
  W_ref = 600 
  T = 0.08*H_ref
  B = 0.12*H_ref
  L = 0.12*W_ref
  R = 0.04*W_ref
  
  c = TCanvas("c","DijetMass cross section with QCD MC",W,H)
  c.GetWindowHeight()
  c.GetWindowWidth()
  c.SetLogy()
  c.Divide(1,2,0,0,0)
  
  
  #------------ pad 1  ----------------
  c.cd(1)
  p11_1 = c.GetPad(1)
  p11_1.SetPad(0.01,0.23,0.99,0.98)
  p11_1.SetLogy()
  p11_1.SetRightMargin(0.05)
  p11_1.SetTopMargin(0.05)
  p11_1.SetFillColor(0)
  p11_1.SetBorderMode(0)
  p11_1.SetFrameFillStyle(0)
  p11_1.SetFrameBorderMode(0)
  
  #Pave text
  pave_fit = TPaveText(0.1558691,0.30735043,0.3750171,0.4070085,"NDC")
  pave_fit = TPaveText(0.2058691,0.20735043,0.4750171,0.3670085,"NDC")
    
  pave_fit.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3")
  pave_fit.AddText("M_{jj} > 1.2 TeV")
  pave_fit.AddText("Wide Jets")
  pave_fit.SetFillColor(0)
  pave_fit.SetLineColor(0)
  pave_fit.SetFillStyle(0)
  pave_fit.SetBorderSize(0)
  pave_fit.SetTextFont(42)
  pave_fit.SetTextSize(0.040)
  pave_fit.SetTextAlign(12) 
  
  
  vFrame = p11_1.DrawFrame(minX_mass,0.000005,maxX_mass,5.0)
  
  vFrame.SetTitle("")
  vFrame.SetXTitle("Dijet Mass (GeV)")
  vFrame.SetYTitle("d#sigma / dm_{jj}   (pb / GeV)")
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  vFrame.GetYaxis().SetTitleOffset(0.95)
  vFrame.GetYaxis().SetLabelSize(0.05)
  
  g_mc.SetLineWidth(2)
  g_mc.SetLineColor(kBlue)
  g_mc.Draw("c")
  g.SetMarkerSize(0.9)
  g.SetMarkerStyle(20)
  g.Draw("pe0 same")
    
  leg = TLegend(0.5564991,0.55,0.8903575,0.705812)
  #leg =  TLegend(0.5564991,0.55,0.8903575,0.80)
  leg.SetTextSize(0.03546853)
  leg.SetLineColor(0)
  leg.SetLineStyle(1)
  leg.SetLineWidth(1)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetMargin(0.35)
  leg.AddEntry(g,"data" ,"PL")
  leg.AddEntry(g_mc,"mc","L")
  leg.Draw("same")
  pave_fit.Draw("same")
  
  # writing the lumi information and the CMS "logo"
  #  CMS_lumi( p11_1, iPeriod, iPos );
  #redraw axis
  p11_1.RedrawAxis()
  p11_1.Update()
  p11_1.GetFrame().Draw()
  #draw the lumi text on the canvas
  CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)
  
  #--- Next PAD
  
  c.cd(2)
  p11_2 = c.GetPad(2)
  p11_2.SetPad(0.01,0.02,0.99,0.24)
  p11_2.SetBottomMargin(0.35)
  p11_2.SetRightMargin(0.05)
  p11_2.SetGridx()
  p11_2.SetGridy()
  
  vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -3., p11_1.GetUxmax(), 3.)
  
  vFrame2.SetTitle("")
  vFrame2.SetXTitle("Dijet Mass (GeV)")
  vFrame2.GetXaxis().SetTitleSize(0.06)
  vFrame2.SetYTitle("(Data-MC)/#sigma")
  vFrame2.GetYaxis().SetTitleSize(0.15)
  vFrame2.GetYaxis().SetTitleOffset(0.40)
  vFrame2.GetYaxis().SetLabelSize(0.09)
  vFrame2.GetXaxis().SetTitleSize(0.18)
  vFrame2.GetXaxis().SetTitleOffset(0.90)
  vFrame2.GetXaxis().SetLabelSize(0.15)
  
  hist_mc_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
  hist_mc_residual_vsMass.GetYaxis().SetRangeUser(-3.,3.)
  hist_mc_residual_vsMass.SetLineWidth(0)
  hist_mc_residual_vsMass.SetFillColor(kBlue)
  hist_mc_residual_vsMass.SetLineColor(1)
  hist_mc_residual_vsMass.Draw("SAMEHIST")
  
  line = TLine(minX_mass,0,maxX_mass,0)
  line.Draw("")
  p11_2.RedrawAxis()
  line2=TLine()
  line2.DrawLine(p11_2.GetUxmin(), p11_2.GetUymax(), p11_2.GetUxmax(), p11_2.GetUymax())
  line2.DrawLine(p11_2.GetUxmax(), p11_2.GetUymin(), p11_2.GetUxmax(), p11_2.GetUymax())
  	
  ### Output files
  
  output_root_file = "dijetFitResults_MC_%s.root" % (fileNameSuffix) 
  
  f_output = TFile(output_root_file,"RECREATE")
  f_output.cd()
  g.Write()
  #hist_mass_original.Write()
  #hist_binned.Write()
  hist_mass.Write()
  c.Write()
  #r_bin->Write()
  f_output.Close()
  c_fileName = "MCandResiduals_%s.png" %(fileNameSuffix)
  c.SaveAs(c_fileName)
  c_fileName = "MCandResiduals_%s.pdf" %(fileNameSuffix)
  c.SaveAs(c_fileName)


def DrawFit(g,M1Bkg,hist_fit_residual_vsMass,FunctionType,nPar,fileNameSuffix):
#  //### Draw plots
  W = 600
  H = 650
  H_ref = 650 
  W_ref = 600 
  T = 0.08*H_ref
  B = 0.12*H_ref
  L = 0.12*W_ref
  R = 0.04*W_ref
  
  c = TCanvas("c","DijetMass cross section with Fit and QCD MC",W,H)
  c.GetWindowHeight()
  c.GetWindowWidth()
  c.SetLogy()
  c.Divide(1,2,0,0,0)
  
  
  #------------ pad 1  ----------------
  c.cd(1)
  p11_1 = c.GetPad(1)
  p11_1.SetPad(0.01,0.23,0.99,0.98)
  p11_1.SetLogy()
  p11_1.SetRightMargin(0.05)
  p11_1.SetTopMargin(0.05)
  p11_1.SetFillColor(0)
  p11_1.SetBorderMode(0)
  p11_1.SetFrameFillStyle(0)
  p11_1.SetFrameBorderMode(0)
  
  #Pave text
  pave_fit = TPaveText(0.1558691,0.30735043,0.3750171,0.4070085,"NDC")
  pave_fit = TPaveText(0.2058691,0.20735043,0.4750171,0.3670085,"NDC")
    
  pave_fit.AddText("|#eta| < 2.5, |#Delta#eta| < 1.3")
  pave_fit.AddText("M_{jj} > 1.1 TeV")
  pave_fit.AddText("Wide Jets")
  pave_fit.SetFillColor(0)
  pave_fit.SetLineColor(0)
  pave_fit.SetFillStyle(0)
  pave_fit.SetBorderSize(0)
  pave_fit.SetTextFont(42)
  pave_fit.SetTextSize(0.040)
  pave_fit.SetTextAlign(12) 
  
  
  vFrame = p11_1.DrawFrame(minX_mass,0.000005,maxX_mass,5.0)
  
  vFrame.SetTitle("")
  vFrame.SetXTitle("Dijet Mass (GeV)")
  vFrame.SetYTitle("d#sigma / dm_{jj}   (pb / GeV)")
  vFrame.GetXaxis().SetTitleSize(0.06)
  vFrame.GetXaxis().SetTitleOffset(0.95)
  vFrame.GetXaxis().SetLabelSize(0.05)
  vFrame.GetYaxis().SetTitleSize(0.06)
  vFrame.GetYaxis().SetTitleOffset(0.95)
  vFrame.GetYaxis().SetLabelSize(0.05)
  
  g.SetMarkerSize(0.9)
  g.SetMarkerStyle(20)
  g.Draw("pe0")
  M1Bkg.SetLineWidth(2)
  M1Bkg.SetLineStyle(2)
  M1Bkg.SetLineColor(2)
  M1Bkg.Draw("same")
    
  leg = TLegend(0.5564991,0.55,0.8903575,0.705812)
  #leg =  TLegend(0.5564991,0.55,0.8903575,0.80)
  leg.SetTextSize(0.03546853)
  leg.SetLineColor(0)
  leg.SetLineStyle(1)
  leg.SetLineWidth(1)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  leg.SetMargin(0.35)
  leg.AddEntry(hist_mass,"data" ,"PL")
  leg.AddEntry(M1Bkg,"fit to data","L")
  leg.Draw("same")
  pave_fit.Draw("same")
  
  # writing the lumi information and the CMS "logo"
  #  CMS_lumi( p11_1, iPeriod, iPos );
  #redraw axis
  p11_1.RedrawAxis()
  p11_1.Update()
  p11_1.GetFrame().Draw()
  #draw the lumi text on the canvas
  CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)
  
  #--- Next PAD
  
  c.cd(2)
  p11_2 = c.GetPad(2)
  p11_2.SetPad(0.01,0.02,0.99,0.24)
  p11_2.SetBottomMargin(0.35)
  p11_2.SetRightMargin(0.05)
  p11_2.SetGridx()
  p11_2.SetGridy()
  
  vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -3., p11_1.GetUxmax(), 3.)
  
  vFrame2.SetTitle("")
  vFrame2.SetXTitle("Dijet Mass (GeV)")
  vFrame2.GetXaxis().SetTitleSize(0.06)
  vFrame2.SetYTitle("(Data-Fit)/#sigma")
  vFrame2.GetYaxis().SetTitleSize(0.15)
  vFrame2.GetYaxis().SetTitleOffset(0.40)
  vFrame2.GetYaxis().SetLabelSize(0.09)
  vFrame2.GetXaxis().SetTitleSize(0.18)
  vFrame2.GetXaxis().SetTitleOffset(0.90)
  vFrame2.GetXaxis().SetLabelSize(0.15)
  
  hist_fit_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
  hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-3.,3.)
  hist_fit_residual_vsMass.SetLineWidth(0)
  hist_fit_residual_vsMass.SetFillColor(2)
  hist_fit_residual_vsMass.SetLineColor(1)
  hist_fit_residual_vsMass.Draw("SAMEHIST")
  
  line = TLine(minX_mass,0,maxX_mass,0)
  line.Draw("")
  p11_2.RedrawAxis()
  line2=TLine()
  line2.DrawLine(p11_2.GetUxmin(), p11_2.GetUymax(), p11_2.GetUxmax(), p11_2.GetUymax())
  line2.DrawLine(p11_2.GetUxmax(), p11_2.GetUymin(), p11_2.GetUxmax(), p11_2.GetUymax())
  	
  ### Output files
  
  output_root_file = "dijetFitResults_FuncType%d_nParFit%d_%s.root" % (FunctionType,nPar,fileNameSuffix) 
  
  f_output = TFile(output_root_file,"RECREATE")
  f_output.cd()
  g.Write()
  #hist_mass_original.Write()
  #hist_binned.Write()
  hist_mass.Write()
  c.Write()
  #r_bin->Write()
  f_output.Close()
  c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.png" %(FunctionType,nPar,fileNameSuffix)
  c.SaveAs(c_fileName)
  c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.pdf" %(FunctionType,nPar,fileNameSuffix)
  c.SaveAs(c_fileName)

#if __name__ == "__main__":

#----- keep the GUI alive ------------
if __name__ == '__main__':
  main() 
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]



#  TCanvas *Canvas0 = new TCanvas("Canvas0","Canvas0",11,51,700,500);
#  Canvas0->cd();
#  Canvas0->SetLogy();
#  hist_mass->GetYaxis()->SetTitle("Events");
#  hist_mass->Draw();  
#  M1Bkg->SetLineColor(1);
#  M1Bkg->Draw("same");     
#
#  TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1",11,51,700,500);
#  Canvas1->cd();
#  Canvas1->SetLogy();
#  hist_mass_varbin->GetYaxis()->SetTitle("Events / bin width");
#  hist_mass_varbin->Draw();  
#  M1Bkg->SetLineColor(1);
#  M1Bkg->Draw("same");     
# 
#  TCanvas *Canvas2 = new TCanvas("Canvas2","Canvas2",11,51,700,500);
#  Canvas2->cd();
#  Canvas2->SetGridx();
#  Canvas2->SetGridy();
#  Canvas2->SetLogx();
#  hist_fit_residual_vsMass->GetYaxis()->SetLimits(-5,5);
#  hist_fit_residual_vsMass->GetYaxis()->SetRangeUser(-5,5);
#  hist_fit_residual_vsMass->GetYaxis()->SetTitle("(data - fit) / #sqrt{data}");
#  hist_fit_residual_vsMass->GetXaxis()->SetRangeUser(minX_mass,maxX_mass);
#  hist_fit_residual_vsMass->GetXaxis()->SetTitle("M_{jj} WideJets [GeV]");
#  hist_fit_residual_vsMass->Draw();
#
#  TCanvas *Canvas3 = new TCanvas("Canvas3","Canvas3",11,51,700,500);
#  Canvas3->cd();
#  hist_fit_residual->GetXaxis()->SetTitle("(data - fit) / #sqrt{data}");
#  hist_fit_residual->GetYaxis()->SetTitle("Number of bins");
#  hist_fit_residual->GetYaxis()->SetRangeUser(0,number_of_variableWidth_bins/3);
#  hist_fit_residual->Draw();
#  hist_fit_residual->Fit("gaus","L","",-3,3);
#
#  //### Output files
#  char output_root_file[500];
#  sprintf(output_root_file,"dijetFitResults_FuncType%d_nParFit%d_%s.root",FunctionType,nPar,fileNameSuffix); 
#
#  TFile f_output(output_root_file,"RECREATE");
#  f_output.cd();
#  Canvas0->Write();
#  Canvas1->Write();
#  Canvas2->Write();
#  Canvas3->Write();
#  f_output.Close();
#
#  //### Save figures from canvas
#  char c0_fileName[200];
#  sprintf(c0_fileName,"dijetmass_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
#  char c1_fileName[200];
#  sprintf(c1_fileName,"dijetmass_varbin_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
#  char c2_fileName[200];
#  sprintf(c2_fileName,"fitresiduals_vs_mass_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
#  char c3_fileName[200];
#  sprintf(c3_fileName,"fitresiduals_FuncType%d_nParFit%d_%s.png",FunctionType,nPar,fileNameSuffix);
#
#  Canvas0->SaveAs(c0_fileName);
#  Canvas1->SaveAs(c1_fileName);
#  Canvas2->SaveAs(c2_fileName);
#  Canvas3->SaveAs(c3_fileName);
#}
