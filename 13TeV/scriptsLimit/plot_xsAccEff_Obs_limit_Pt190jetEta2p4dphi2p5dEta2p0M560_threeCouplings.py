#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array

gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.045, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.045, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
t_m = 0.06  ##top margin
b_m = 0.11   ##botton margin
l_m = 0.115  ##left margin
r_m = 0.04  ##right margin
gStyle.SetPadTopMargin(t_m)
gStyle.SetPadBottomMargin(b_m)
gStyle.SetPadLeftMargin(l_m)
gStyle.SetPadRightMargin(r_m)
gROOT.ForceStyle()


BR = 1.0


masses_gev = array('d', [1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0, 4400.0, 4600.0, 4800.0, 5000.0, 5200.0, 5400.0])
masses_tev = array('d', [ 0.001 * float(b) for b in masses_gev ])

xs_obs_limits_f1p0  = array('d', [0.0385836, 0.0328597, 0.0362624, 0.0307328, 0.0312984, 0.0233855, 0.011458, 0.00749401, 0.00391534, 0.00352713, 0.00381675, 0.00333843, 0.00265381, 0.00253488, 0.00230159, 0.00207015, 0.00181926, 0.0016805, 0.00161387, 0.00157083, 0.00156268, 0.00154644, 0.00154871])
xs_obs_limits_f1p0_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_f1p0 ])

xs_obs_limits_f0p5  = array('d', [0.0347955, 0.0293301, 0.0299501, 0.0225903, 0.0252556, 0.0205909, 0.00937511, 0.00691599, 0.00323465, 0.00313366, 0.00358346, 0.00311181, 0.00239662, 0.002439, 0.00223139, 0.00192447, 0.00171072, 0.00162075, 0.00156918, 0.00153427, 0.0016165, 0.00159358, 0.00149949])
xs_obs_limits_f0p5_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_f0p5 ])

xs_obs_limits_f0p1  = array('d', [0.0334031, 0.0281143, 0.0281072, 0.0202043, 0.023329, 0.0195144, 0.00868901, 0.00678801, 0.00307312, 0.00300843, 0.00348401, 0.00304116, 0.00228442, 0.00238282, 0.00221733, 0.00186647, 0.00170358, 0.00161163, 0.00155223, 0.00154747, 0.00152153, 0.0014992, 0.00147819])
xs_obs_limits_f0p1_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_f0p1 ])


masses_qstar      = array('d', [ 1.0,       1.1,       1.2,       1.3,       1.4,       1.5,       1.6,       1.7,       1.8,       1.9,
                                 2.0,       2.1,       2.2,       2.3,       2.4,       2.5,       2.6,       2.7,       2.8,       2.9,
				 3.0,       3.1,       3.2,       3.3,       3.4,       3.5,       3.6,       3.7,       3.8,       3.9,
				 4.0,       4.1,       4.2,       4.3,       4.4,       4.5,       4.6,       4.7,       4.8,       4.9,
				 5.0])

eff_qstar_MID    = array('d',  [ 0.439106,  0.4488,    0.4585,    0.4682,    0.4779,    0.4874,    0.4969,    0.50655,   0.5162,    0.525834,
                                 0.535469,  0.538234,  0.541,     0.5436,    0.5462,    0.54875,   0.5513,    0.5539,    0.5565,    0.55893,
				 0.561357,  0.562628,  0.5639,    0.56515,   0.5664,    0.567649,  0.5689,    0.5701499, 0.5714,    0.5725475,
				 0.573695,  0.573897,  0.5741,    0.57428,   0.57446,   0.574629,  0.5748,    0.57505,   0.5753,    0.575473,
				 0.575647])

xs_qstar_f1p0     = array('d', [ 1.635e+01, 1.068e+01, 7.036e+00, 4.827e+00, 3.363e+00, 2.373e+00, 1.709e+00, 1.243e+00, 9.292e-01, 6.885e-01,
                                 5.244e-01, 3.951e-01, 3.010e-01, 2.322e-01, 1.805e-01, 1.402e-01, 1.099e-01, 8.650e-02, 6.787e-02, 5.372e-02,
				 4.273e-02, 3.391e-02, 2.720e-02, 2.186e-02, 1.744e-02, 1.417e-02, 1.126e-02, 9.062e-03, 7.276e-03, 5.911e-03,
				 4.814e-03, 3.870e-03, 3.156e-03, 2.554e-03, 2.057e-03, 1.656e-03, 1.354e-03, 1.089e-03, 8.813e-04, 7.214e-04,
				 5.836e-04])
                                  
xs_qstar_f0p5     = array('d', [ 4.137e+00, 2.642e+00, 1.768e+00, 1.217e+00, 8.445e-01, 6.012e-01, 4.345e-01, 3.179e-01, 2.342e-01, 1.765e-01,
                                 1.328e-01, 1.005e-01, 7.712e-02, 5.922e-02, 4.583e-02, 3.601e-02, 2.799e-02, 2.206e-02, 1.746e-02, 1.378e-02,
				 1.096e-02, 8.642e-03, 7.002e-03, 5.531e-03, 4.407e-03, 3.554e-03, 2.860e-03, 2.302e-03, 1.851e-03, 1.488e-03,
				 1.211e-03, 9.753e-04, 7.847e-04, 6.374e-04, 5.156e-04, 4.187e-04, 3.360e-04, 2.728e-04, 2.189e-04, 1.770e-04,
				 1.437e-04])
                                  
xs_qstar_f0p1     = array('d', [ 1.655e-01, 1.057e-01, 7.134e-02, 4.932e-02, 3.421e-02, 2.440e-02, 1.750e-02, 1.284e-02, 9.433e-03, 7.075e-03,
                                 5.298e-03, 4.025e-03, 3.107e-03, 2.374e-03, 1.861e-03, 1.431e-03, 1.130e-03, 8.902e-04, 7.051e-04, 5.527e-04,
				 4.363e-04, 3.511e-04, 2.784e-04, 2.232e-04, 1.791e-04, 1.435e-04, 1.148e-04, 9.267e-05, 7.459e-05, 6.014e-05,
				 4.852e-05, 3.902e-05, 3.157e-05, 2.536e-05, 2.058e-05, 1.677e-05, 1.344e-05, 1.087e-05, 8.690e-06, 7.102e-06, 
				 5.739e-06])


xsEff_f1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_MID)])
xsEff_f1p0_fb = array('d', [1000.0 * float(b) for b in xsEff_f1p0])

xsEff_f0p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p5, eff_qstar_MID)])
xsEff_f0p5_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p5])

xsEff_f0p1 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p1, eff_qstar_MID)])
xsEff_f0p1_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p1])

graph_obs_f1p0 = TGraph(len(masses_tev),masses_tev,xs_obs_limits_f1p0_fb)
graph_obs_f1p0.GetXaxis().SetTitle("q* Mass [TeV]")
graph_obs_f1p0.GetYaxis().SetTitle("#sigma #times B #times A #times #epsilon [fb]")
graph_obs_f1p0.GetYaxis().SetRangeUser(5e-1,300)
graph_obs_f1p0.GetXaxis().SetNdivisions(510)
graph_obs_f1p0.GetXaxis().SetLimits(0.6,6)
graph_obs_f1p0.GetYaxis().CenterTitle()
graph_obs_f1p0.GetYaxis().SetLabelSize(0.04)
graph_obs_f1p0.GetYaxis().SetTitleOffset(1.1)
graph_obs_f1p0.GetXaxis().CenterTitle()
graph_obs_f1p0.GetXaxis().SetLabelSize(0.04)
graph_obs_f1p0.GetXaxis().SetTitleOffset(1.1)
graph_obs_f1p0.GetXaxis().CenterTitle()
graph_obs_f1p0.SetLineWidth(2)
graph_obs_f1p0.SetLineStyle(1)
graph_obs_f1p0.SetLineColor(2)
graph_obs_f1p0.SetMarkerStyle(21)
graph_obs_f1p0.SetMarkerColor(2)

graph_obs_f0p5 = TGraph(len(masses_tev),masses_tev,xs_obs_limits_f0p5_fb)
graph_obs_f0p5.SetLineWidth(2)
graph_obs_f0p5.SetLineStyle(1)
graph_obs_f0p5.SetLineColor(3)
graph_obs_f0p5.SetMarkerStyle(34)
graph_obs_f0p5.SetMarkerColor(3)

graph_obs_f0p1 = TGraph(len(masses_tev),masses_tev,xs_obs_limits_f0p1_fb)
graph_obs_f0p1.SetLineWidth(2)
graph_obs_f0p1.SetLineStyle(1)
graph_obs_f0p1.SetLineColor(4)
graph_obs_f0p1.SetMarkerStyle(29)
graph_obs_f0p1.SetMarkerColor(4)

graph_qstar_f1p0 = TGraph(len(masses_qstar),masses_qstar,xsEff_f1p0_fb)
graph_qstar_f1p0.SetLineWidth(2)                
graph_qstar_f1p0.SetLineColor(3)

graph_qstar_f0p5 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p5_fb)
graph_qstar_f0p5.SetLineWidth(2)                
graph_qstar_f0p5.SetLineColor(2)

graph_qstar_f0p1 = TGraph(len(masses_qstar),masses_qstar,xsEff_f0p1_fb)
graph_qstar_f0p1.SetLineWidth(2)                
graph_qstar_f0p1.SetLineColor(4)

c = TCanvas("c", "",800,800)
c.cd()

graph_obs_f1p0.Draw("ALP")
graph_obs_f0p5.Draw("LP")
graph_obs_f0p1.Draw("LP")

graph_qstar_f1p0.Draw("L")
graph_qstar_f0p5.Draw("L")
graph_qstar_f0p1.Draw("L")

##legend = TLegend(.55,.69,.85,.92)  ## for pas twiki
legend = TLegend(.52,.54,.85,.74)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetHeader('95% CL upper limits')
legend.AddEntry(graph_obs_f1p0,"Observed limit (f = 1.0)","lp")
legend.AddEntry(graph_obs_f0p5,"Observed limit (f = 0.5)","lp")
legend.AddEntry(graph_obs_f0p1,"Observed limit (f = 0.1)","lp")
legend.Draw()

legend1 = TLegend(.16,.16,.47,.34)
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetFillStyle(0)
legend1.SetTextFont(42)
legend1.SetTextSize(0.035)
legend1.AddEntry(graph_qstar_f1p0,"q* (f = 1.0)","l")
legend1.AddEntry(graph_qstar_f0p5,"q* (f = 0.5)","l")
legend1.AddEntry(graph_qstar_f0p1,"q* (f = 0.1)","l")
legend1.Draw()

lumiTextSize = 0.6
lumiTextOffset = 0.2
lumi = TLatex()
lumi.SetNDC()
lumi.SetTextAngle(0)
lumi.SetTextColor(kBlack)
lumi.SetTextFont(42)
lumi.SetTextAlign(31)
lumi.SetTextSize(lumiTextSize*t_m)
lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "2.2 fb^{-1} (13 TeV)")

cmsTextFont = 61
cmsTextSize = 0.75
extraTextFont = 52
extraOverCmsTextSize = 0.76
relExtraDY = 1.2
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
cms.DrawLatex(0.9, 0.9, "CMS")

### if extra text (Unpublished or Preliminary)
cms.SetTextFont(extraTextFont);
cms.SetTextAlign(33);
cms.SetTextSize(extraOverCmsTextSize*cmsTextSize*t_m);
##cms.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t_m, extraText);
cms.DrawLatex(0.9, 0.9 - relExtraDY*cmsTextSize*t_m, "Preliminary");

l1 = TLatex()
l1.SetNDC()
l1.SetTextAlign(13)
l1.SetTextFont(62)
l1.SetTextSize(0.04)
l1.SetTextFont(42)
l1.DrawLatex(0.785,0.805, "q*#rightarrow q#gamma") ##0.21,0.26


gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('Limit_xsAccEff_Obs_limit_Pt190jetEta2p4dphi2p5dEta2p0M560_CompareThreeCoup.pdf')
c.SaveAs('Limit_xsAccEff_Obs_limit_Pt190jetEta2p4dphi2p5dEta2p0M560_CompareThreeCoup.png')
