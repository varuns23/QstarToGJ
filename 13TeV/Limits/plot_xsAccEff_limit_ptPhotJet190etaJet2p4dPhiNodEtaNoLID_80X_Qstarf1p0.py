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


masses_gev = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0])
masses_tev = array('d', [ 0.001 * float(b) for b in masses_gev ])

xs_obs_limits  = array('d', [0.0232455, 0.0423358, 0.0363525, 0.02109, 0.0129538, 0.0103374, 0.0101107, 0.00758186, 0.00673586, 0.00568584, 0.00417352, 0.00379691, 0.00360912, 0.00263492, 0.00268502, 0.00321043, 0.00257726, 0.00194614, 0.00201415, 0.00210489, 0.00186416, 0.00150732, 0.00117112, 0.000962343, 0.00072086, 0.000600649, 0.00054003, 0.000543183, 0.0005525, 0.000544745, 0.000541047, 0.000528261, 0.000532373, 0.000523621, 0.000524023, 0.000512963, 0.000484502, 0.000448446, 0.000400122, 0.000348737, 0.000328568, 0.000282174, 0.000258208, 0.000237665, 0.000216714, 0.000202162, 0.000188485, 0.000180298, 0.000174827, 0.00016933])
xs_obs_limits_fb = array('d', [1000.0 * float(b) for b in xs_obs_limits ])

xs_exp_limits = array('d', [0.0319726, 0.0266186, 0.0206399, 0.0162276, 0.0140934, 0.0118763, 0.0100359, 0.0085482, 0.00726779, 0.00649961, 0.00520807, 0.00445246, 0.00395093, 0.00371607, 0.00305094, 0.00264514, 0.00242881, 0.00206417, 0.00187502, 0.0016213, 0.00140891, 0.00111219, 0.00101252, 0.000993473, 0.000914791, 0.000813974, 0.000738422, 0.000647738, 0.000531606, 0.000502216, 0.000486476, 0.000452818, 0.000386228, 0.000353286, 0.000312697, 0.000299492, 0.000281972, 0.00025838, 0.000257906, 0.000254771, 0.000247583, 0.000247025, 0.000234274, 0.000221434, 0.000210617, 0.000199807, 0.000189405, 0.000184108, 0.000177085, 0.000169924])
xs_exp_limits_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits ])

masses_exp_gev = array('d', [1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4500.0, 4600.0, 4700.0, 4800.0, 4900.0, 5000.0, 5100.0, 5200.0, 5300.0, 5400.0, 5500.0, 5600.0, 5700.0, 5800.0, 5900.0, 5900.0, 5800.0, 5700.0, 5600.0, 5500.0, 5400.0, 5300.0, 5200.0, 5100.0, 5000.0, 4900.0, 4800.0, 4700.0, 4600.0, 4500.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0, 1100.0, 1000.0])
masses_exp_tev = array('d', [ 0.001 * float(b) for b in masses_exp_gev ])

xs_exp_limits_1sigma  = array('d', [0.021114, 0.017909, 0.0145982, 0.0125503, 0.00954397, 0.00838801, 0.00691202, 0.00663239, 0.00524134, 0.00456755, 0.00386654, 0.00322311, 0.00268026, 0.00256502, 0.00226636, 0.00195592, 0.00174733, 0.00151518, 0.00131643, 0.00108516, 0.000947015, 0.000813341, 0.000783291, 0.000731446, 0.000620437, 0.000551813, 0.000539718, 0.000477782, 0.000393348, 0.000373275, 0.000354707, 0.000315487, 0.00029086, 0.000268545, 0.000248342, 0.000229733, 0.000212789, 0.000196311, 0.000201464, 0.000185132, 0.000189076, 0.000176244, 0.000170915, 0.000166633, 0.000161036, 0.000157422, 0.000153764, 0.000153098, 0.000151023, 0.000149815, 0.000234436, 0.000244153, 0.000253744, 0.000264092, 0.000269427, 0.000280742, 0.000291016, 0.000299285, 0.000310625, 0.000338126, 0.000330369, 0.000364114, 0.000371057, 0.000391823, 0.000424011, 0.000450613, 0.000508435, 0.00052978, 0.000614078, 0.000646427, 0.000715766, 0.000755623, 0.000892642, 0.00101408, 0.00115421, 0.00122065, 0.00135174, 0.00137919, 0.00167455, 0.00185252, 0.00231797, 0.00257573, 0.00303732, 0.00342313, 0.00367077, 0.00430615, 0.0049588, 0.00561712, 0.00628734, 0.00706287, 0.00934879, 0.00986279, 0.0109497, 0.0135751, 0.0173363, 0.0196777, 0.0215524, 0.028428, 0.0367653, 0.0448567])
xs_exp_limits_1sigma_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits_1sigma ])

xs_exp_limits_2sigma = array('d', [0.0169604, 0.0130364, 0.0106182, 0.00930662, 0.00747079, 0.00677251, 0.00511113, 0.00482341, 0.00388114, 0.003411, 0.0028603, 0.00225597, 0.00207889, 0.00208974, 0.0016246, 0.0014824, 0.00124427, 0.00105664, 0.00101012, 0.000865304, 0.000755533, 0.000660794, 0.000581475, 0.000522194, 0.000462986, 0.000404517, 0.000418123, 0.000367188, 0.000310569, 0.000306527, 0.00028297, 0.000244457, 0.000248232, 0.000201012, 0.000190297, 0.000175191, 0.000167834, 0.00016894, 0.000152226, 0.000157877, 0.000158333, 0.000154018, 0.000147269, 0.000149496, 0.000147356, 0.000146953, 0.000141391, 0.000141384, 0.000141037, 0.000136971, 0.000320816, 0.000330336, 0.000346268, 0.000363819, 0.000380632, 0.000385884, 0.000405001, 0.000446448, 0.000476452, 0.000469692, 0.000449247, 0.000437508, 0.000540654, 0.000578068, 0.000612357, 0.00066166, 0.000691368, 0.000781295, 0.000790427, 0.000888056, 0.000911034, 0.00108521, 0.00123295, 0.00128583, 0.00139228, 0.00180676, 0.00202874, 0.002175, 0.00219553, 0.00260135, 0.00322551, 0.00339132, 0.00430823, 0.00436464, 0.00514261, 0.00617126, 0.00669262, 0.00730707, 0.00836957, 0.0099852, 0.0134358, 0.0138605, 0.0160834, 0.0190717, 0.0220428, 0.0267119, 0.0296548, 0.038651, 0.0505649, 0.0634715])
xs_exp_limits_2sigma_fb = array('d', [1000.0 * float(b) for b in xs_exp_limits_2sigma ])


masses_qstar      = array('d', [ 1.0,       1.1,       1.2,       1.3,       1.4,       1.5,       1.6,       1.7,       1.8,       1.9,
                                 2.0,       2.1,       2.2,       2.3,       2.4,       2.5,       2.6,       2.7,       2.8,       2.9,
				 3.0,       3.1,       3.2,       3.3,       3.4,       3.5,       3.6,       3.7,       3.8,       3.9,
				 4.0,       4.1,       4.2,       4.3,       4.4,       4.5,       4.6,       4.7,       4.8,       4.9,
				 5.0,       5.1,       5.2,       5.3,       5.4,       5.5,       5.6,       5.7,       5.8,       5.9,
				 6.0])

eff_qstar_LID    = array('d',  [ .585588, .598179, .610769, .62336, .635951, .648542, .661132, .673723, .686314, .698904, 
                                 .711495, .714421, .717346, .720272, .723198, .726123, .729049, .731975, .734901, .737826, 
				 .740752, .741555, .742358, .743161, .743964, .744767, .74557, .746373, .747176, .747979, 
				 .748782, .748902, .749022, .749141, .749261, .749381, .749501, .749621, .74974, .74986, 
				 .74998, .748183, .746387, .74459, .742793, .740996, .7392, .737403, .735606, .73381, 
				 .73201])
                                
			
		
	

xs_qstar_f1p0     = array('d', [ 1.635e+01, 1.068e+01, 7.036e+00, 4.827e+00, 3.363e+00, 2.373e+00, 1.709e+00, 1.243e+00, 9.292e-01, 6.885e-01,
                                 5.244e-01, 3.951e-01, 3.010e-01, 2.322e-01, 1.805e-01, 1.402e-01, 1.099e-01, 8.650e-02, 6.787e-02, 5.372e-02,
				 4.273e-02, 3.391e-02, 2.720e-02, 2.186e-02, 1.744e-02, 1.417e-02, 1.126e-02, 9.062e-03, 7.276e-03, 5.911e-03,
				 4.814e-03, 3.870e-03, 3.156e-03, 2.554e-03, 2.057e-03, 1.656e-03, 1.354e-03, 1.089e-03, 8.813e-04, 7.214e-04,
				 5.836e-04, 4.734e-04, 3.807e-04, 3.108e-04, 2.517e-04, 2.051e-04, 1.650e-04, 1.339e-04, 1.072e-04, 8.685e-05,
                                 7.085e-05])
                                  
xs_qstar_f0p5     = array('d', [ 4.137e+00, 2.642e+00, 1.768e+00, 1.217e+00, 8.445e-01, 6.012e-01, 4.345e-01, 3.179e-01, 2.342e-01, 1.765e-01,
                                 1.328e-01, 1.005e-01, 7.712e-02, 5.922e-02, 4.583e-02, 3.601e-02, 2.799e-02, 2.206e-02, 1.746e-02, 1.378e-02,
				 1.096e-02, 8.642e-03, 7.002e-03, 5.531e-03, 4.407e-03, 3.554e-03, 2.860e-03, 2.302e-03, 1.851e-03, 1.488e-03,
				 1.211e-03, 9.753e-04, 7.847e-04, 6.374e-04, 5.156e-04, 4.187e-04, 3.360e-04, 2.728e-04, 2.189e-04, 1.770e-04,
				 1.437e-04])
                                  
xs_qstar_f0p1     = array('d',  [ 1.655e-01, 1.057e-01, 7.134e-02, 4.932e-02, 3.421e-02, 2.440e-02, 1.750e-02, 1.284e-02, 9.433e-03, 7.075e-03,
                                 5.298e-03, 4.025e-03, 3.107e-03, 2.374e-03, 1.861e-03, 1.431e-03, 1.130e-03, 8.902e-04, 7.051e-04, 5.527e-04,
				 4.363e-04, 3.511e-04, 2.784e-04, 2.232e-04, 1.791e-04, 1.435e-04, 1.148e-04, 9.267e-05, 7.459e-05, 6.014e-05,
				 4.852e-05, 3.902e-05, 3.157e-05, 2.536e-05, 2.058e-05, 1.677e-05, 1.344e-05, 1.087e-05, 8.690e-06, 7.102e-06, 
				 5.739e-06])


xsEff_f1p0 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f1p0, eff_qstar_LID)])
xsEff_f1p0_fb = array('d', [1000.0 * float(b) for b in xsEff_f1p0])

xsEff_f0p5 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p5, eff_qstar_LID)])
xsEff_f0p5_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p5])

xsEff_f0p1 = array('d', [float(b) * float(m) for b,m in zip(xs_qstar_f0p1, eff_qstar_LID)])
xsEff_f0p1_fb = array('d', [1000.0 * float(b) for b in xsEff_f0p1])


graph_exp_2sigma = TGraph(len(masses_exp_tev),masses_exp_tev,xs_exp_limits_2sigma_fb)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times B #times A #times #epsilon [fb]")
graph_exp_2sigma.GetYaxis().SetRangeUser(8e-2,300)
graph_exp_2sigma.GetXaxis().SetNdivisions(510)
graph_exp_2sigma.GetXaxis().SetLimits(0.6,6)

graph_exp_2sigma.GetYaxis().CenterTitle()
graph_exp_2sigma.GetYaxis().SetLabelSize(0.04)
graph_exp_2sigma.GetYaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetXaxis().CenterTitle()
graph_exp_2sigma.GetXaxis().SetLabelSize(0.04)
graph_exp_2sigma.GetXaxis().SetTitleOffset(1.1)
graph_exp_2sigma.GetXaxis().CenterTitle()

graph_exp_1sigma = TGraph(len(masses_exp_tev),masses_exp_tev,xs_exp_limits_1sigma_fb)
graph_exp_1sigma.SetFillColor(kGreen+1)

graph_exp = TGraph(len(masses_tev),masses_tev,xs_exp_limits_fb)
#graph_exp.SetMarkerStyle(24)
graph_exp.SetLineWidth(2)
graph_exp.SetLineStyle(2)
graph_exp.SetLineColor(4)

graph_obs = TGraph(len(masses_tev),masses_tev,xs_obs_limits_fb)
graph_obs.SetMarkerStyle(20)
graph_obs.SetLineWidth(2)
#graph_obs.SetLineStyle(1)
graph_obs.SetLineColor(1)

graph_qstar = TGraph(len(masses_qstar),masses_qstar,xsEff_f1p0_fb)
graph_qstar.SetLineWidth(2)                
graph_qstar.SetLineColor(2)

c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
graph_qstar.Draw("L")

##legend = TLegend(.55,.69,.85,.92)  ## for pas twiki
legend = TLegend(.58,.52,.89,.72)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetHeader('95% CL upper limits')
legend.AddEntry(graph_obs,"Observed limit","lp")
legend.AddEntry(graph_exp,"Expected limit","lp")
legend.AddEntry(graph_exp_1sigma,"Expected limit #pm 1#sigma","f")
legend.AddEntry(graph_exp_2sigma,"Expected limit #pm 2#sigma","f")
legend.Draw()

legend1 = TLegend(.16,.18,.47,.24)
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetFillStyle(0)
legend1.SetTextFont(42)
legend1.SetTextSize(0.035)
legend1.AddEntry(graph_qstar,"Excited quark (f = 1.0)","l")
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
lumi.DrawLatex(1 - r_m, 1 - t_m +  lumiTextOffset*t_m, "24.49 fb^{-1} (13 TeV)")

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
c.SaveAs('Limit_ObseExp_xsBRAccEff_ptPhotJet190etaJet2p4dPhiNodEtaNoLID_80X_Qstarf1p0.pdf')


