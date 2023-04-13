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

xs_obs_limits_All  = array('d', [0.0385836, 0.0328597, 0.0362624, 0.0307328, 0.0312984, 0.0233855, 0.011458, 0.00749401, 0.00391534, 0.00352713, 0.00381675, 0.00333843, 0.00265381, 0.00253488, 0.00230159, 0.00207015, 0.00184926, 0.0017205, 0.00164387, 0.00158083, 0.00158268, 0.00155644, 0.00154871])
xs_obs_limits_All_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_All ])

xs_obs_limits_Bkg  = array('d', [0.0373792, 0.0332222, 0.0339418, 0.02835, 0.0294415, 0.0218021, 0.0108877, 0.00732365, 0.00366479, 0.00343063, 0.00375058, 0.00330562, 0.0026211, 0.00252512, 0.00229473, 0.0020914, 0.00187781, 0.00171292, 0.00167785, 0.00162974, 0.00158326, 0.00159953, 0.00153226])
xs_obs_limits_Bkg_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_Bkg ])

xs_obs_limits_JES  = array('d', [0.0341481, 0.0274387, 0.0265453, 0.0212814, 0.0215731, 0.0179948, 0.00944963, 0.0066615, 0.00355723, 0.00329204, 0.00357586, 0.0031632, 0.00256294, 0.00245946, 0.00227873, 0.0020732, 0.00184783, 0.00171752, 0.00164106, 0.00157807, 0.0015813, 0.0015529, 0.00151673])
xs_obs_limits_JES_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_JES ])

xs_obs_limits_JER  = array('d', [0.03371, 0.0273352, 0.026055, 0.0207562, 0.021231, 0.0173484, 0.00927204, 0.00648455, 0.0034866, 0.0032714, 0.00357098, 0.00315633, 0.00252489, 0.00246574, 0.0022776, 0.00205882, 0.0018438, 0.00171686, 0.00164104, 0.00157779, 0.00158228, 0.00155076, 0.00151975])
xs_obs_limits_JER_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_JER ])

xs_obs_limits_Lumi  = array('d', [0.0339283, 0.0275311, 0.0261922, 0.0208042, 0.0213025, 0.0175492, 0.00929609, 0.00654876, 0.00348757, 0.0032717, 0.00360006, 0.00317669, 0.00252677, 0.00248203, 0.0022745, 0.00206811, 0.00185606, 0.00172721, 0.00165208, 0.00160263, 0.00160353, 0.00156549, 0.00152376])
xs_obs_limits_Lumi_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_Lumi ])

xs_obs_limits_Stat  = array('d', [0.0337138, 0.0273437, 0.0259661, 0.020629, 0.0210575, 0.017347, 0.0092202, 0.00649092, 0.00346004, 0.00324657, 0.0035681, 0.00315061, 0.00250857, 0.00246287, 0.00227806, 0.00205344, 0.0018422, 0.00171434, 0.00163972, 0.00157857, 0.00157927, 0.00155494, 0.00151404])
xs_obs_limits_Stat_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits_Stat ])


xs_obs_limits_All_Ratio = array('d', [float(b) / float(m) for b,m in zip( xs_obs_limits_All_fb,xs_obs_limits_Stat_fb) ])
xs_obs_limits_Bkg_Ratio = array('d', [float(b) / float(m) for b,m in zip( xs_obs_limits_Bkg_fb,xs_obs_limits_Stat_fb) ])
xs_obs_limits_JES_Ratio = array('d', [float(b) / float(m) for b,m in zip( xs_obs_limits_JES_fb,xs_obs_limits_Stat_fb) ])
xs_obs_limits_JER_Ratio = array('d', [float(b) / float(m) for b,m in zip( xs_obs_limits_JER_fb,xs_obs_limits_Stat_fb) ])
xs_obs_limits_Lumi_Ratio = array('d', [float(b) / float(m) for b,m in zip( xs_obs_limits_Lumi_fb,xs_obs_limits_Stat_fb) ])

graph_Ratio_All = TGraph(len(masses_tev),masses_tev,xs_obs_limits_All_Ratio)
graph_Ratio_All.SetFillColor(kYellow)
graph_Ratio_All.GetXaxis().SetTitle("q* Mass [TeV]")
graph_Ratio_All.GetYaxis().SetTitle("Limit Ratio")
graph_Ratio_All.GetYaxis().SetRangeUser(0.8,1.8)
graph_Ratio_All.GetXaxis().SetNdivisions(510)
graph_Ratio_All.GetXaxis().SetLimits(0.6,6)

graph_Ratio_All.GetYaxis().CenterTitle()
graph_Ratio_All.GetYaxis().SetLabelSize(0.04)
graph_Ratio_All.GetYaxis().SetTitleOffset(1.1)
graph_Ratio_All.GetXaxis().CenterTitle()
graph_Ratio_All.GetXaxis().SetLabelSize(0.04)
graph_Ratio_All.GetXaxis().SetTitleOffset(1.1)
graph_Ratio_All.GetXaxis().CenterTitle()
graph_Ratio_All.SetLineStyle(1)
graph_Ratio_All.SetLineColor(1)
graph_Ratio_All.SetMarkerStyle(8)
graph_Ratio_All.SetMarkerColor(1)

graph_Ratio_Bkg = TGraph(len(masses_tev),masses_tev,xs_obs_limits_Bkg_Ratio)
graph_Ratio_Bkg.SetLineStyle(7)
graph_Ratio_Bkg.SetLineColor(6)
graph_Ratio_Bkg.SetLineWidth(2)
graph_Ratio_Bkg.SetMarkerStyle(29)
graph_Ratio_Bkg.SetMarkerSize(1.2)
graph_Ratio_Bkg.SetMarkerColor(6)

graph_Ratio_JES = TGraph(len(masses_tev),masses_tev,xs_obs_limits_JES_Ratio)
graph_Ratio_JES.SetLineStyle(7)
graph_Ratio_JES.SetLineWidth(2)
graph_Ratio_JES.SetLineColor(4)
graph_Ratio_JES.SetMarkerStyle(25)
graph_Ratio_JES.SetMarkerSize(1.1)
graph_Ratio_JES.SetMarkerColor(4)

graph_Ratio_JER = TGraph(len(masses_tev),masses_tev,xs_obs_limits_JER_Ratio)
graph_Ratio_JER.SetLineStyle(2)
graph_Ratio_JER.SetLineWidth(2)
graph_Ratio_JER.SetLineColor(3)
graph_Ratio_JER.SetMarkerStyle(28)
graph_Ratio_JER.SetMarkerSize(1.2)
graph_Ratio_JER.SetMarkerColor(3)

graph_Ratio_Lumi = TGraph(len(masses_tev),masses_tev,xs_obs_limits_Lumi_Ratio)
graph_Ratio_Lumi.SetLineStyle(1)
graph_Ratio_Lumi.SetLineWidth(2)
graph_Ratio_Lumi.SetLineColor(2)
graph_Ratio_Lumi.SetMarkerStyle(32)
graph_Ratio_Lumi.SetMarkerSize(1.2)
graph_Ratio_Lumi.SetMarkerColor(2)



c = TCanvas("c", "",800,800)
c.cd()

graph_Ratio_All.Draw("ALP")
graph_Ratio_Bkg.Draw("LP")
graph_Ratio_JES.Draw("LP")
graph_Ratio_JER.Draw("LP")
graph_Ratio_Lumi.Draw("LP")

##legend = TLegend(.55,.69,.85,.92)  ## for pas twiki
legend = TLegend(.58,.52,.89,.72)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
#legend.SetHeader('95% CL upper limits')
legend.AddEntry(graph_Ratio_All,"All/Stat. only","lp")
legend.AddEntry(graph_Ratio_Bkg,"Bkg/Stat. only","lp")
legend.AddEntry(graph_Ratio_Lumi,"Lumi/Stat. only","lp")
legend.AddEntry(graph_Ratio_JES,"(JES + PES)/Stat. only","lp")
legend.AddEntry(graph_Ratio_JER,"(JER + PER)/Stat. only","lp")
legend.Draw()

legend1 = TLegend(.16,.18,.47,.24)
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetFillStyle(0)
legend1.SetTextFont(42)
legend1.SetTextSize(0.035)
#legend1.AddEntry(graph_qstar,"Excited quark (f = 1.0)","l")
#legend1.Draw()

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

#c.SetLogy()
c.SaveAs('Ratio_SysbyStat_xsAccEff_Pt190jetEta2p4dphi2p5dEta2p0M560MID.pdf')
c.SaveAs('Ratio_SysbyStat_xsAccEff_Pt190jetEta2p4dphi2p5dEta2p0M560MID.png')


