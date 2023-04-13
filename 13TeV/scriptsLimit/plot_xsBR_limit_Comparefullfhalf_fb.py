#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.055, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.055, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.1)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.1)
gROOT.ForceStyle()


masses_fhalf = array('d', [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0])

xs_obs_limits_fhalf_pb = array('d', [0.032172199999999998, 0.024260199999999999, 0.010855800000000001, 0.0065376799999999997, 0.0062754100000000004, 0.0082482000000000007, 0.0061187899999999998, 0.0042923099999999997, 0.0045939199999999996, 0.00316117, 0.00172926, 0.00128024, 0.0013677800000000001, 0.0012333800000000001, 0.0012960000000000001, 0.00143367, 0.0010339399999999999, 0.00092689699999999999, 0.00087992699999999997, 0.00062764299999999999, 0.00037498400000000003, 0.00031172799999999999, 0.000324381, 0.00032325600000000001, 0.00027191299999999998, 0.000216765, 0.00019314899999999999, 0.000184747, 0.00017984899999999999, 0.00017568799999999999, 0.00017252500000000001, 0.000169982, 0.00017652499999999999, 0.000175985, 0.00017442, 0.00017157500000000001, 0.000170388, 0.00016942899999999999])
xs_obs_limits_fhalf = array('d', [1000.0 * float(b) for b in xs_obs_limits_fhalf_pb ])

masses_full = array('d', [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0])
xs_obs_limits_full_pb = array('d', [0.043261300000000003, 0.0288829, 0.0127066, 0.0071527300000000004, 0.0073611700000000002, 0.0099175300000000008, 0.0079669600000000004, 0.0060442100000000004, 0.0060724100000000003, 0.0041595299999999998, 0.00233022, 0.00159549, 0.00165913, 0.0015969999999999999, 0.0016696899999999999, 0.0016410299999999999, 0.00134162, 0.0010888499999999999, 0.00096837799999999999, 0.00071037599999999998, 0.00043614800000000003, 0.00034915800000000002, 0.000340088, 0.00035707799999999999, 0.00029766599999999998, 0.000251463, 0.00022563499999999999, 0.00020038900000000001, 0.000201984, 0.00020212299999999999, 0.000196539, 0.000191881, 0.00018848000000000001, 0.000197457, 0.00019481500000000001, 0.00019233299999999999, 0.000190063, 0.000188073])
xs_obs_limits_full = array('d', [1000.0 * float(b) for b in xs_obs_limits_full_pb ]) 

## Atimes eff 
AEff_fhalf = array ('d', [0.464278, 0.489212, 0.510378, 0.549197, 0.553, 0.556752, 0.560435, 0.564032, 0.567526, 0.570899, 0.574136, 0.577218, 0.580128, 0.58285, 0.585367, 0.58766, 0.589714, 0.591511, 0.593034, 0.594266, 0.59519, 0.595788, 0.596045, 0.595941, 0.595461, 0.594588, 0.593304, 0.591592, 0.589435, 0.586816, 0.583718, 0.580124, 0.576016, 0.571379, 0.566193, 0.560444, 0.554112, 0.547182])
AEff_full = array ('d', [0.464278, 0.489212, 0.510378, 0.528148, 0.542879, 0.554907, 0.564549, 0.572106, 0.577858, 0.582068, 0.58498, 0.586821, 0.587797, 0.588096, 0.587891, 0.587332, 0.586553, 0.585669, 0.584776, 0.583953, 0.583258, 0.582735, 0.582404, 0.58227, 0.58232, 0.58252, 0.582819, 0.583148, 0.583418, 0.583523, 0.583339, 0.582721, 0.581508, 0.579519, 0.576556, 0.572402, 0.56682, 0.559556])

##Ge new obser, expe, 1sigma and 2 sigma
xs_obs_limits_n_fhalf = array('d', [float(b) / float(m) for b,m in zip(xs_obs_limits_fhalf, AEff_fhalf )])
xs_obs_limits_n_full = array('d', [float(b) / float(m) for b,m in zip(xs_obs_limits_full, AEff_full )])


masses_qstar_fhalf = array('d', [   1.5,     2.0,      2.5,      2.8,     2.85,      2.9,     2.95,         3.0,     3.05,      3.1,          3.5,       4.0,       4.5 ])
xs_qstar_fhalf_pb  = array('d', [ 0.105, 0.01484, 0.002425, 8.565e-4, 7.245e-4, 6.082e-4, 5.083e-4,   4.304e-04, 3.607e-4, 3.044e-4,    7.586e-05, 1.258e-05, 1.929e-06 ])
xs_qstar_fhalf     = array('d', [1000.0 * float(b) for b in xs_qstar_fhalf_pb ])
eff_qstar_5a_fhalf = array ('d', [ 1.0,  1.0, 1.0, 1.0,  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ])

xs_5a_fhalf = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_fhalf, eff_qstar_5a_fhalf)])

masses_qstar_full = array('d', [    1.5,    1.7,     2.0,      2.5,      3.0,      3.3,     3.35,      3.4,     3.45,       3.5,     3.55,      3.6,       4.0,       4.5 ])
xs_qstar_full_pb  = array('d', [ 0.4124, 0.1843, 0.05858, 0.009768, 0.001755, 6.366e-4, 5.428e-4, 4.662e-4, 3.947e-4, 3.385e-04, 2.821e-4, 2.388e-4, 6.241e-05, 1.270e-05 ])
xs_qstar_full     = array('d', [1000.0 * float(b) for b in xs_qstar_full_pb ])
eff_qstar_5a_full = array ('d', [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

xs_5a_full = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_full, eff_qstar_5a_full)])

result_fhalf= array('d',[0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4])
result_full = array('d',[0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4])

graph_obs_fhalf = TGraph(len(masses_fhalf),result_fhalf,xs_obs_limits_n_fhalf)
graph_obs_fhalf.SetLineColor(kRed)
#graph_obs.SetLineWidth(2)
graph_obs_fhalf.SetMarkerStyle(27)
graph_obs_fhalf.SetMarkerColor(kRed)
graph_obs_fhalf.GetXaxis().SetNdivisions(510)
graph_obs_fhalf.GetXaxis().SetTitle("q* Mass [TeV]")
graph_obs_fhalf.GetXaxis().CenterTitle()
graph_obs_fhalf.GetXaxis().SetTitleSize(0.05)
graph_obs_fhalf.GetXaxis().SetTitleOffset(0.94)
graph_obs_fhalf.GetXaxis().SetLabelSize(0.045)
graph_obs_fhalf.GetXaxis().SetLimits(0.4,4.7)

graph_obs_fhalf.GetYaxis().SetTitle("#sigma #times B [fb]")
graph_obs_fhalf.GetYaxis().CenterTitle()
graph_obs_fhalf.GetYaxis().SetTitleSize(0.05)
graph_obs_fhalf.GetYaxis().SetTitleOffset(0.94)
graph_obs_fhalf.GetYaxis().SetLabelSize(0.045)
graph_obs_fhalf.GetYaxis().SetRangeUser(1e-01,300)

graph_obs_full = TGraph(len(masses_full),result_full,xs_obs_limits_n_full)
graph_obs_full.SetLineColor(kBlue)
graph_obs_full.SetMarkerStyle(24)
graph_obs_full.SetMarkerColor(kBlue)
graph_obs_full.GetYaxis().SetTitleOffset(1.1)

graph_qstar_fhalf = TGraph(len(masses_qstar_fhalf),masses_qstar_fhalf,xs_5a_fhalf)
graph_qstar_fhalf.SetLineWidth(2)                
graph_qstar_fhalf.SetLineColor(kRed)
graph_qstar_fhalf.SetLineStyle(2)

graph_qstar_full = TGraph(len(masses_qstar_full),masses_qstar_full,xs_5a_full)
graph_qstar_full.SetLineWidth(2)                
graph_qstar_full.SetLineColor(kBlue)
graph_qstar_full.SetLineStyle(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_obs_fhalf.Draw("ALP")
graph_qstar_fhalf.Draw("L")
graph_obs_full.Draw("LP")
graph_qstar_full.Draw("L")

legend = TLegend(.53,.67,.83,.9)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
legend.SetHeader('95% CL Upper Limits')
legend.AddEntry(graph_obs_full,"Observed Limit (f = 1.0)","lp")
legend.AddEntry(graph_qstar_full,"Excited Quark (f = 1.0)","l")
legend.AddEntry(graph_obs_fhalf,"Observed Limit (f = 0.5)","lp")
legend.AddEntry(graph_qstar_fhalf,"Excited Quark (f = 0.5)","l")
legend.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.035)
l1.DrawLatex(0.16,0.37, "CMS Preliminary")
l1.DrawLatex(0.16,0.32, "#sqrt{s} = 8 TeV")
l1.DrawLatex(0.16,0.26, "#intLdt = 19.7 fb^{-1}")
l1.DrawLatex(0.17,0.19, "q* #rightarrow q#gamma")

gPad.RedrawAxis();

c.SetLogy()

#c.SaveAs('ExcitedQuarksToGJ_Comparef0p5f1p0_ObseExp_xs_Limits.pdf')
c.SaveAs('ExcitedQuarksToGJ_Comparef0p5f1p0_ObseExp_xs_Limits_TEMP.pdf')

c.Update() 
#img = TImage.Create()
#img.FromPad(c)   
#img.WriteImage('ExcitedQuarksToGJ_Comparef0p5f1p0_ObseExp_xs_Limits.png')
#img.WriteImage('ExcitedQuarksToGJ_Comparef0p5f1p0_ObseExp_xs_Limits.gif')

#c.SaveAs('ExcitedQuarksToGJ_Comparef0p5f1p0_ObseExp_xs_Limits.eps')
#c.SaveAs('ExcitedQuarksToGJ_Comparef0p5f1p0_ObseExp_xs_Limits.png')

