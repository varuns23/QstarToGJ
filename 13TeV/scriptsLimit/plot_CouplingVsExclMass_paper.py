#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.045, "XYZ") #0.06
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.045, "XYZ") #0.06
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
t_m = 0.06  ##top margin 0.55
b_m = 0.11   ##botton margin 0.13
l_m = 0.12  ##left margin 0.13
r_m = 0.12  ##right margin 0.06
gStyle.SetPadTopMargin(t_m) ##0.06
gStyle.SetPadBottomMargin(b_m) ##0.15
gStyle.SetPadLeftMargin(l_m)  ###0.15
gStyle.SetPadRightMargin(r_m)  ###0.05
gROOT.ForceStyle()

Coupling = array('d', [ 0.06,  0.07,   0.1,  0.13,  0.15,  0.17,   0.2,   0.3,  0.35,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9,  1.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.2, 0.17, 0.15, 0.13, 0.1, 0.07, 0.06  ])

ExcludedMassBoth = array('d', [ 1.056, 1.114, 1.237, 1.472, 1.855, 1.923, 2.028, 2.306, 2.429, 2.811, 2.932, 3.009, 3.194, 3.315, 3.416, 3.50, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,  0.7, 0.7, 0.7, 0.7, 0.7, 0.7 ]) 

ExpectedCoupling = array('d', [ 0.05, 0.06, 0.07, 0.1, 0.13, 0.15, 0.17, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ExpectedExcludedMassBoth = array('d', [ 0.7, 0.879, 1.022, 1.339, 1.595, 1.726, 1.856, 2.042, 2.409, 2.538, 2.68, 2.867, 2.949, 3.054, 3.159, 3.257, 3.34 ])

grCoupling = TGraph(len(Coupling),ExcludedMassBoth,Coupling)
grCoupling.SetFillColor(kRed-4)
grCoupling.SetLineColor(kRed-4)

grCoupling.GetYaxis().SetTitle("Couplings [f = f_{s} = f']")
grCoupling.GetYaxis().CenterTitle()
grCoupling.GetYaxis().SetTitleSize(0.05)
grCoupling.GetYaxis().SetTitleOffset(1.2) ##1.1
grCoupling.GetYaxis().SetLabelSize(0.045) ##0.05
grCoupling.GetYaxis().SetLabelOffset(0.01) 
grCoupling.GetYaxis().SetRangeUser(0.0,1.0)


RightYaxis = TGaxis( 4., 0.0, 4., 1.0,0,1,510,"+L")
RightYaxis.SetLabelFont(42)
RightYaxis.SetTitle("M_{q*} / #Lambda")
RightYaxis.SetTitleFont(42)
RightYaxis.SetTitleSize(0.05)
RightYaxis.SetTitleOffset(1.1)
RightYaxis.CenterTitle()

grCoupling.GetXaxis().SetTitle("q* Mass [TeV]")
grCoupling.GetXaxis().CenterTitle()
grCoupling.GetXaxis().SetTitleSize(0.05)
grCoupling.GetXaxis().SetTitleOffset(1.1)
grCoupling.GetXaxis().SetLabelSize(0.04) ##0.05
grCoupling.GetXaxis().SetLabelOffset(0.008)
grCoupling.GetXaxis().SetLimits(0.7,4.)
grCoupling.GetXaxis().SetNdivisions(510)


grExpectedCoupling = TGraph(len(ExpectedCoupling),ExpectedExcludedMassBoth,ExpectedCoupling)
grExpectedCoupling.SetLineColor(kBlue+3)
grExpectedCoupling.SetLineStyle(9)
grExpectedCoupling.SetLineWidth(2)

c = TCanvas("c", "",800,800)
c.cd()

grCoupling.Draw("AFL")
grExpectedCoupling.Draw("L")
RightYaxis.Draw()
    
#legend = TLegend(0.37,0.14,0.68,0.22)
legend = TLegend(0.37,0.15,0.68,0.23)
#legend = TLegend(0.43,0.18,0.78,0.26)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)  ##35
legend.AddEntry(grExpectedCoupling,"Expected excluded q* reqion","L") 
legend.AddEntry(grCoupling,"Observed excluded q* region","F")
legend.Draw()

lumiTextSize = 0.6
lumiTextOffset = 0.2
lumi = TLatex()
lumi.SetNDC()
lumi.SetTextAngle(0)
lumi.SetTextColor(kBlack)
lumi.SetTextFont(42)
lumi.SetTextAlign(31)
lumi.SetTextSize(lumiTextSize*t_m)
lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "19.7 fb^{-1} (8 TeV)")

cmsTextFont = 61
cmsTextSize = 0.75
##posX_ = l_m     + 0.045 * (1 - l_m - r_m)  ##Top left
##posX_ = l_m     + 0.5   * (1 - l_m - r_m)   ## Centre
posX_ = 1 - r_m + 0.045 * (1 - l_m - r_m)  ## Top right
posY_ = 1 - t_m - 0.035 * (1 - t_m - b_m)

cms =  TLatex()
cms.SetNDC()
cms.SetTextFont(cmsTextFont)
cms.SetTextSize(cmsTextSize * t_m)
cms.SetTextAlign(33)  ### 11-top left;  21-top centre;  31-top right
##cms.DrawLatex(posX_, posY_, "CMS")
cms.DrawLatex(0.845, 0.88, "CMS")
#cms.DrawLatex(0.85, 0.9, "CMS")

l1 = TLatex()
l1.SetNDC()
l1.SetTextAlign(13)
l1.SetTextSize(0.04) ##0.035
l1.SetTextFont(42)
#l1.DrawLatex(0.75,0.85, "q* #rightarrow q#gamma")
l1.DrawLatex(0.735,0.79, "q*#rightarrow q#gamma")
#l1.DrawLatex(0.74,0.8, "q*#rightarrow q#gamma")

gPad.RedrawAxis();
c.SaveAs('CouplingvsMass__TEMP.pdf')
#c.SaveAs('CouplingvsMass.pdf')
#c.SaveAs('CouplingvsMass_paper.pdf')
#c.SaveAs('CouplingvsMass.png')
#c.SaveAs('CouplingvsMass_paper.eps')

