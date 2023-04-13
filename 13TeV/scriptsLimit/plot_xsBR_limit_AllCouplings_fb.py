#!/usr/bin/env python
import sys, os, subprocess, string, re
from ROOT import *
from array import array


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.06, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.11) #0.15 / 0.13
gStyle.SetPadRightMargin(0.03) #0.05 / 0.01
gStyle.SetPadTopMargin(0.06)  ##0.05
gStyle.SetPadBottomMargin(0.11) #0.15
gROOT.ForceStyle()


BR = 1.0



masses = array('d', [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0])
xs_obs_limits_pb  = array('d', [0.043261300000000003, 0.0288829, 0.0127066, 0.0071527300000000004, 0.0073611700000000002, 0.0099175300000000008, 0.0079669600000000004, 0.0060442100000000004, 0.0060724100000000003, 0.0041595299999999998, 0.00233022, 0.00159549, 0.00165913, 0.0015969999999999999, 0.0016696899999999999, 0.0016410299999999999, 0.00134162, 0.0010888499999999999, 0.00096837799999999999, 0.00071037599999999998, 0.00043614800000000003, 0.00034915800000000002, 0.000340088, 0.00035707799999999999, 0.00029766599999999998, 0.000251463, 0.00022563499999999999, 0.00020038900000000001, 0.000201984, 0.00020212299999999999, 0.000196539, 0.000191881, 0.00018848000000000001, 0.000197457, 0.00019481500000000001, 0.00019233299999999999, 0.000190063, 0.000188073])
xs_obs_limits = array('d', [1000.0 * float(b) for b in xs_obs_limits_pb])

masses_fhalf = array('d', [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0])
xs_fhalf_obs_limits_pb  = array('d', [0.032172199999999998, 0.024260199999999999, 0.010855800000000001, 0.0065376799999999997, 0.0062754100000000004, 0.0082482000000000007, 0.0061187899999999998, 0.0042923099999999997, 0.0045939199999999996, 0.00316117, 0.00172926, 0.00128024, 0.0013677800000000001, 0.0012333800000000001, 0.0012960000000000001, 0.00143367, 0.0010339399999999999, 0.00092689699999999999, 0.00087992699999999997, 0.00062764299999999999, 0.00037498400000000003, 0.00031172799999999999, 0.000324381, 0.00032325600000000001, 0.00027191299999999998, 0.000216765, 0.00019314899999999999, 0.000184747, 0.00017984899999999999, 0.00017568799999999999, 0.00017252500000000001, 0.000169982, 0.00017652499999999999, 0.000175985, 0.00017442, 0.00017157500000000001, 0.000170388, 0.00016942899999999999])
xs_fhalf_obs_limits = array('d', [1000.0 * float(b) for b in xs_fhalf_obs_limits_pb]) 

xs_exp_limits_pb = array('d', [0.035456799999999997, 0.024380499999999999, 0.017505199999999999, 0.0133283, 0.0095709200000000001, 0.0076610200000000002, 0.0062512799999999997, 0.0049046899999999997, 0.00376692, 0.0032564999999999998, 0.0027552599999999998, 0.0022690200000000001, 0.00179166, 0.0015299700000000001, 0.0012809500000000001, 0.0010582300000000001, 0.00091599699999999999, 0.00085092800000000004, 0.00077223199999999997, 0.00065015099999999996, 0.00056257400000000004, 0.00048872499999999997, 0.00044408300000000001, 0.00041765999999999998, 0.00040558199999999998, 0.00036096499999999998, 0.00033267300000000002, 0.00031854700000000003, 0.00031468099999999998, 0.00029590399999999999, 0.00027901500000000002, 0.000264391, 0.000251143, 0.00026234000000000001, 0.00025577500000000002, 0.00024911600000000002, 0.00024016, 0.00023427800000000001])
xs_exp_limits = array('d', [1000.0 * float(b) for b in xs_exp_limits_pb])

masses_exp = array('d', [700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0, 1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0, 2800.0, 2900.0, 3000.0, 3100.0, 3200.0, 3300.0, 3400.0, 3500.0, 3600.0, 3700.0, 3800.0, 3900.0, 4000.0, 4100.0, 4200.0, 4300.0, 4400.0, 4400.0, 4300.0, 4200.0, 4100.0, 4000.0, 3900.0, 3800.0, 3700.0, 3600.0, 3500.0, 3400.0, 3300.0, 3200.0, 3100.0, 3000.0, 2900.0, 2800.0, 2700.0, 2600.0, 2500.0, 2400.0, 2300.0, 2200.0, 2100.0, 2000.0, 1900.0, 1800.0, 1700.0, 1600.0, 1500.0, 1400.0, 1300.0, 1200.0, 1100.0, 1000.0, 900.0, 800.0, 700.0])
xs_exp_limits_1sigma_pb  = array('d', [0.023731700000000001, 0.0170703, 0.0129204, 0.0101826, 0.0069303899999999998, 0.0053394599999999999, 0.0040504900000000003, 0.0035441700000000001, 0.00267454, 0.00239916, 0.0019998799999999999, 0.00163528, 0.00139134, 0.0010822399999999999, 0.00088745000000000002, 0.00079694800000000001, 0.00068539700000000002, 0.00060026600000000002, 0.00051125900000000002, 0.00046553000000000002, 0.00042627099999999998, 0.000388878, 0.00035367200000000002, 0.000319816, 0.00030806099999999999, 0.00026832800000000002, 0.00025066200000000002, 0.000239371, 0.00023593400000000001, 0.00022542100000000001, 0.000212802, 0.00021074900000000001, 0.000201077, 0.000224336, 0.00021538699999999999, 0.00021091099999999999, 0.000207238, 0.000205362, 0.00028228800000000001, 0.00029176700000000001, 0.00030314600000000001, 0.00032448699999999999, 0.00034676200000000003, 0.00032883899999999999, 0.00036159700000000001, 0.00038773300000000002, 0.00042802600000000002, 0.00044344600000000001, 0.00045223000000000003, 0.000474105, 0.00050028900000000001, 0.00053848000000000004, 0.00059348699999999996, 0.00064777699999999999, 0.00071690999999999996, 0.00079286899999999997, 0.00089515799999999996, 0.00104039, 0.0011624700000000001, 0.00126694, 0.00151414, 0.00181631, 0.0021101200000000001, 0.0025153300000000001, 0.0031532800000000001, 0.0038069800000000002, 0.0044063699999999997, 0.0055911900000000002, 0.0070334799999999999, 0.0084836700000000004, 0.010704099999999999, 0.014936400000000001, 0.017494099999999999, 0.0246292, 0.032933700000000003, 0.049422199999999999])
xs_exp_limits_1sigma = array('d', [1000.0 * float(b) for b in xs_exp_limits_1sigma_pb])
xs_exp_limits_2sigma_pb = array('d', [0.016375600000000001, 0.012226900000000001, 0.0097508200000000003, 0.0074958899999999998, 0.00510371, 0.0039123099999999996, 0.0032265599999999998, 0.0025695599999999998, 0.0018817599999999999, 0.0016067900000000001, 0.0014405900000000001, 0.0011807499999999999, 0.00097551099999999998, 0.00088917100000000002, 0.00069930700000000003, 0.00062618699999999995, 0.00051909199999999997, 0.000462929, 0.00045401399999999998, 0.000385497, 0.00032846900000000002, 0.000294889, 0.000276218, 0.000256493, 0.00025067399999999999, 0.00023926000000000001, 0.00021855999999999999, 0.00020750200000000001, 0.00020844200000000001, 0.000198235, 0.00019420200000000001, 0.00018962800000000001, 0.00018410900000000001, 0.00020372399999999999, 0.000197388, 0.00019434399999999999, 0.00019121000000000001, 0.00018921000000000001, 0.00038368400000000002, 0.00040363699999999998, 0.000442516, 0.00047251599999999997, 0.00048144399999999998, 0.00045043299999999997, 0.00046281399999999998, 0.000528913, 0.00056158800000000002, 0.00061201500000000004, 0.00061670099999999999, 0.00062449100000000002, 0.00064168700000000005, 0.00073230400000000003, 0.00078713299999999997, 0.00087084399999999996, 0.00096422200000000004, 0.00114248, 0.00125212, 0.00138517, 0.00158241, 0.00186276, 0.0023061399999999999, 0.0024346599999999999, 0.00269345, 0.0034446899999999998, 0.0041143200000000003, 0.0050229000000000003, 0.0061420600000000004, 0.0072316500000000001, 0.0100077, 0.012219000000000001, 0.0143339, 0.020825300000000001, 0.023614699999999999, 0.032725600000000001, 0.044358099999999998, 0.069868299999999994])
xs_exp_limits_2sigma = array('d', [1000.0 * float(b) for b in xs_exp_limits_2sigma_pb])



## Atimes eff 
AEff = array ('d', [0.464278, 0.489212, 0.510378, 0.528148, 0.542879, 0.554907, 0.564549, 0.572106, 0.577858, 0.582068, 0.58498, 0.586821, 0.587797, 0.588096, 0.587891, 0.587332, 0.586553, 0.585669, 0.584776, 0.583953, 0.583258, 0.582735, 0.582404, 0.58227, 0.58232, 0.58252, 0.582819, 0.583148, 0.583418, 0.583523, 0.583339, 0.582721, 0.581508, 0.579519, 0.576556, 0.572402, 0.56682, 0.559556])

AEff_fhalf = array ('d', [0.464278, 0.489212, 0.510378, 0.549197, 0.553, 0.556752, 0.560435, 0.564032, 0.567526, 0.570899, 0.574136, 0.577218, 0.580128, 0.58285, 0.585367, 0.58766, 0.589714, 0.591511, 0.593034, 0.594266, 0.59519, 0.595788, 0.596045, 0.595941, 0.595461, 0.594588, 0.593304, 0.591592, 0.589435, 0.586816, 0.583718, 0.580124, 0.576016, 0.571379, 0.566193, 0.560444, 0.554112, 0.547182])
## for 1Sigma and 2 Sigma
AEff_Sigma = array ('d', [0.464278, 0.489212, 0.510378, 0.528148, 0.542879, 0.554907, 0.564549, 0.572106, 0.577858, 0.582068, 0.58498, 0.586821, 0.587797, 0.588096, 0.587891, 0.587332, 0.586553, 0.585669, 0.584776, 0.583953, 0.583258, 0.582735, 0.582404, 0.58227, 0.58232, 0.58252, 0.582819, 0.583148, 0.583418, 0.583523, 0.583339, 0.582721, 0.581508, 0.579519, 0.576556, 0.572402, 0.56682, 0.559556, 0.559556, 0.56682, 0.572402, 0.576556, 0.579519, 0.581508, 0.582721, 0.583339, 0.583523, 0.583418, 0.583148, 0.582819, 0.58252, 0.58232, 0.58227, 0.582404, 0.582735, 0.583258, 0.583953, 0.584776, 0.585669, 0.586553, 0.587332, 0.587891, 0.588096, 0.587797, 0.586821, 0.58498, 0.582068, 0.577858, 0.572106, 0.564549, 0.554907, 0.542879, 0.528148, 0.510378, 0.489212, 0.464278])
##Ge new obser, expe, 1sigma and 2 sigma
xs_obs_limits_n = array('d', [float(b) / float(m) for b,m in zip(xs_obs_limits, AEff )])
xs_exp_limits_n = array('d', [float(b) / float(m) for b,m in zip(xs_exp_limits, AEff )])
xs_exp_limits_2sigma_n =array('d', [float(b) / float(m) for b,m in zip(xs_exp_limits_2sigma, AEff_Sigma )])
xs_exp_limits_1sigma_n =array('d', [float(b) / float(m) for b,m in zip(xs_exp_limits_1sigma, AEff_Sigma )])

xs_fhalf_obs_limits_n = array('d', [float(b) / float(m) for b,m in zip(xs_fhalf_obs_limits, AEff_fhalf )])
#print xs_obs_limits_n

masses_qstar = array('d', [ 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
## cross secton
xs_qstar_pb = array ('d',[0.4124, 0.1843, 0.05858, 0.009768, 0.001755, 3.241e-04, 6.045e-05, 1.170e-05 ])
xs_qstar = array('d', [1000.0 * float(b) for b in xs_qstar_pb ])
eff_qstar_5a = array ('d', [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

xs_5a = array('d', [float(b) * float(m) for b,m in zip(xs_qstar, eff_qstar_5a)])

mass_qstar_f = array('d', [      1.0,      1.5,      2.0,      2.5,      3.0,      3.5,      4.0,      4.5 ] )

mass_qstar_fb = array('d', [ 0.7,   1.0,      1.5 ] )
xs_qstar_0p04_pb = array('d', [3.962e-2,  6.859e-3,   6.743e-4 ] )
xs_qstar_0p04 = array('d', [1000.0 * float(b) for b in xs_qstar_0p04_pb ])
xs_qstar_0p05_pb = array('d', [6.219e-2,  1.078e-2,   1.047e-3 ] )
xs_qstar_0p05 = array('d', [1000.0 * float(b) for b in xs_qstar_0p05_pb ])
xs_qstar_0p07_pb = array('d', [1.217e-1,  2.102e-2,   2.074e-3 ] )
xs_qstar_0p07 = array('d', [1000.0 * float(b) for b in xs_qstar_0p07_pb ])

xs_qstar_0p1_pb = array('d', [ 4.268e-2, 4.222e-3, 5.940e-4, 9.721e-5, 1.700e-5, 2.989e-6, 4.841e-7, 7.163e-8 ] )
xs_qstar_0p1 = array('d', [1000.0 * float(b) for b in xs_qstar_0p1_pb ])
xs_qstar_0p2_pb = array('d', [ 1.707e-1, 1.699e-2, 2.377e-3, 3.921e-4, 6.831e-5, 1.189e-5, 1.950e-6, 2.839e-7 ] )
xs_qstar_0p2 = array('d', [1000.0 * float(b) for b in xs_qstar_0p2_pb ])
xs_qstar_0p3_pb = array('d', [ 3.863e-1, 3.782e-2, 5.298e-3, 8.754e-4, 1.547e-4, 2.657e-5, 4.389e-6, 6.518e-7 ] )
xs_qstar_0p3 = array('d', [1000.0 * float(b) for b in xs_qstar_0p3_pb ])
xs_qstar_0p4_pb = array('d', [ 6.838e-1, 6.692e-2, 9.497e-3, 1.546e-3, 2.738e-4, 4.784e-5, 7.837e-6, 1.180e-6 ] )
xs_qstar_0p4 = array('d', [1000.0 * float(b) for b in xs_qstar_0p4_pb ])
xs_qstar_0p5_pb = array('d', [ 1.061e-0, 1.047e-1, 1.481e-2, 2.436e-3, 4.298e-4, 7.575e-5, 1.258e-5, 1.928e-6 ] )
xs_qstar_0p5 = array('d', [1000.0 * float(b) for b in xs_qstar_0p5_pb ])
xs_qstar_0p6_pb = array('d', [ 1.541e-0, 1.507e-1, 2.129e-2, 3.522e-3, 6.233e-4, 1.101e-4, 1.882e-5, 3.023e-6 ] )
xs_qstar_0p6 = array('d', [1000.0 * float(b) for b in xs_qstar_0p6_pb ])
xs_qstar_0p7_pb = array('d', [ 2.088e-0, 2.050e-1, 2.897e-2, 4.756e-3, 8.479e-4, 1.512e-4, 2.659e-5, 4.574e-6 ] )
xs_qstar_0p7 = array('d', [1000.0 * float(b) for b in xs_qstar_0p7_pb ])
xs_qstar_0p8_pb = array('d', [ 2.695e-0, 2.653e-1, 3.751e-2, 6.212e-3, 1.111e-3, 2.009e-4, 3.602e-5, 6.422e-6 ] )
xs_qstar_0p8 = array('d', [1000.0 * float(b) for b in xs_qstar_0p8_pb ])
xs_qstar_0p9_pb = array('d', [ 3.388e-0, 3.364e-1, 4.736e-2, 7.898e-3, 1.404e-3, 2.601e-4, 4.727e-5, 8.772e-6 ] )
xs_qstar_0p9 = array('d', [1000.0 * float(b) for b in xs_qstar_0p9_pb ])
xs_qstar_1p0_pb = array('d', [ 4.177e-0, 4.121e-1, 5.830e-2, 9.718e-3, 1.761e-3, 3.245e-4, 6.041e-5, 1.174e-5 ] )
xs_qstar_1p0 = array('d', [1000.0 * float(b) for b in xs_qstar_1p0_pb ])

result = array('d',[0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4])
result_fhalf = array('d',[0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4])

result_2sigma= array('d',[0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.4,4.3,4.2,4.1,4.0,3.9,3.8,3.7,3.6,3.5,3.4,3.3,3.2,3.1,3.0,2.9,2.8,2.7,2.6,2.5,2.4,2.3,2.2,2.1,2.0,1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7])

graph_exp_2sigma = TGraph(len(masses_exp),result_2sigma,xs_exp_limits_2sigma_n)
graph_exp_2sigma.SetFillColor(kYellow)
graph_exp_2sigma.GetYaxis().SetTitle("#sigma #times B [fb]")
graph_exp_2sigma.GetYaxis().SetTitleSize(0.05) #0.06
graph_exp_2sigma.GetYaxis().SetTitleOffset(0.95) #1.05
graph_exp_2sigma.GetYaxis().SetLabelSize(0.045) #0.05
graph_exp_2sigma.GetYaxis().SetRangeUser(1e-02,300)
graph_exp_2sigma.GetYaxis().CenterTitle()

graph_exp_2sigma.GetXaxis().SetTitle("q* Mass [TeV]")
graph_exp_2sigma.GetXaxis().SetNdivisions(510)
graph_exp_2sigma.GetXaxis().SetLimits(0.4,4.7)
graph_exp_2sigma.GetXaxis().SetTitleSize(0.05) #0.06
graph_exp_2sigma.GetXaxis().SetTitleOffset(0.9)
graph_exp_2sigma.GetXaxis().SetLabelSize(0.045) #0.05
graph_exp_2sigma.GetXaxis().CenterTitle()

graph_exp_1sigma = TGraph(len(masses_exp),result_2sigma,xs_exp_limits_1sigma_n)
graph_exp_1sigma.SetFillColor(kGreen+1)

graph_exp = TGraph(len(masses),result,xs_exp_limits_n)
#graph_exp.SetMarkerStyle(24)
graph_exp.SetLineWidth(2)
graph_exp.SetLineStyle(2)
graph_exp.SetLineColor(4)

graph_obs = TGraph(len(masses),result,xs_obs_limits_n)
graph_obs.SetMarkerStyle(20)
graph_obs.SetLineWidth(2)
#graph_obs.SetLineStyle(1)
graph_obs.SetLineColor(1)

graph_fhalf_obs = TGraph(len(masses_fhalf),result_fhalf,xs_fhalf_obs_limits_n)
graph_fhalf_obs.SetMarkerStyle(20)
graph_fhalf_obs.SetLineWidth(2)
#graph_fhalf_obs.SetLineStyle(1)
graph_fhalf_obs.SetLineColor(1)

graph_qstar = TGraph(len(masses_qstar),masses_qstar,xs_5a)
graph_qstar.SetLineWidth(2)                
graph_qstar.SetLineColor(2)

graph_qstar_0p1 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p1) 
graph_qstar_0p1.SetLineWidth(2)
graph_qstar_0p1.SetLineColor(49)

graph_qstar_0p2 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p2) 
graph_qstar_0p2.SetLineWidth(2)
graph_qstar_0p2.SetLineColor(46)

graph_qstar_0p3 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p3) 
graph_qstar_0p3.SetLineWidth(2)
graph_qstar_0p3.SetLineColor(41)

graph_qstar_0p4 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p4) 
graph_qstar_0p4.SetLineWidth(2)
graph_qstar_0p4.SetLineColor(9)

graph_qstar_0p5 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p5) 
graph_qstar_0p5.SetLineWidth(2)
graph_qstar_0p5.SetLineColor(7)

graph_qstar_0p6 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p6) 
graph_qstar_0p6.SetLineWidth(2)
graph_qstar_0p6.SetLineColor(6)

graph_qstar_0p7 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p7) 
graph_qstar_0p7.SetLineWidth(2)
graph_qstar_0p7.SetLineColor(4)

graph_qstar_0p8 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p8) 
graph_qstar_0p8.SetLineWidth(2)
graph_qstar_0p8.SetLineColor(2)

graph_qstar_0p9 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_0p9) 
graph_qstar_0p9.SetLineWidth(2)
graph_qstar_0p9.SetLineColor(13)

graph_qstar_1p0 = TGraph(len(mass_qstar_f), mass_qstar_f, xs_qstar_1p0) 
graph_qstar_1p0.SetLineWidth(2)
graph_qstar_1p0.SetLineColor(4)

graph_qstar_0p07 = TGraph(len(mass_qstar_fb), mass_qstar_fb, xs_qstar_0p07) 
graph_qstar_0p07.SetLineWidth(2)
graph_qstar_0p07.SetLineStyle(2)
graph_qstar_0p07.SetLineColor(41)

graph_qstar_0p05 = TGraph(len(mass_qstar_fb), mass_qstar_fb, xs_qstar_0p05) 
graph_qstar_0p05.SetLineWidth(2)
graph_qstar_0p05.SetLineStyle(2)
graph_qstar_0p05.SetLineColor(3)

graph_qstar_0p04 = TGraph(len(mass_qstar_fb), mass_qstar_fb, xs_qstar_0p04) 
graph_qstar_0p04.SetLineWidth(2)
graph_qstar_0p04.SetLineStyle(2)
graph_qstar_0p04.SetLineColor(2)


c = TCanvas("c", "",800,800)
c.cd()

graph_exp_2sigma.Draw("AF")
graph_exp_1sigma.Draw("F")
graph_exp.Draw("L")
graph_obs.Draw("LP")
#graph_fhalf_obs.Draw("LP")
#graph_qstar.Draw("L")

graph_qstar_0p04.Draw("L")
graph_qstar_0p05.Draw("L")
graph_qstar_0p07.Draw("L")
graph_qstar_0p1.Draw("L")
graph_qstar_0p2.Draw("L")
graph_qstar_0p3.Draw("L")
graph_qstar_0p4.Draw("L")
graph_qstar_0p5.Draw("L")
graph_qstar_0p6.Draw("L")
graph_qstar_0p7.Draw("L")
graph_qstar_0p8.Draw("L")
graph_qstar_0p9.Draw("L")
graph_qstar_1p0.Draw("L")

legend = TLegend(.59,.63,.91,.8)
legend.SetBorderSize(0)
legend.SetFillColor(0)
legend.SetFillStyle(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
legend.SetHeader('95% CL Upper Limits')
legend.AddEntry(graph_obs,"Observed Limit","lp")
legend.AddEntry(graph_exp,"Expected Limit","lp")
#legend.AddEntry(graph_qstar,"Excited Quark ( f = 1.0 )","l")
legend.AddEntry(graph_exp_1sigma,"Expected Limit #pm 1#sigma","f")
legend.AddEntry(graph_exp_2sigma,"Expected Limit #pm 2#sigma","f")
legend.Draw()

#leg = TLegend(.17,.14,.43,.46)
leg = TLegend(.14,.14,.4,.46)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.03)
#legend.AddText('Excited Quark with coupling')
leg.AddEntry(graph_qstar_0p04,"f = 0.04","l")
leg.AddEntry(graph_qstar_0p05,"f = 0.05","l")
leg.AddEntry(graph_qstar_0p1,"f = 0.1","l")
leg.AddEntry(graph_qstar_0p2,"f = 0.2","l")
leg.AddEntry(graph_qstar_0p3,"f = 0.3","l")
leg.AddEntry(graph_qstar_0p4,"f = 0.4","l")
leg.AddEntry(graph_qstar_0p5,"f = 0.5","l")
leg.AddEntry(graph_qstar_0p6,"f = 0.6","l")
leg.AddEntry(graph_qstar_0p7,"f = 0.7","l")
leg.AddEntry(graph_qstar_0p8,"f = 0.8","l")
leg.AddEntry(graph_qstar_0p9,"f = 0.9","l")
leg.AddEntry(graph_qstar_1p0,"f = 1.0","l")
leg.Draw()

l1 = TLatex()
l1.SetTextAlign(12)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)
l1.SetTextSize(0.03)
l1.DrawLatex(0.6,0.89, "CMS Preliminary")
l1.DrawLatex(0.6,0.83, "#intLdt = 19.7 fb^{-1}, #sqrt{s} = 8 TeV")
l1.DrawLatex(0.15,0.49, "q* with couplings")
#l1.DrawLatex(0.2,0.19, "q* with couplings")
l1.SetTextSize(0.04)

gPad.RedrawAxis();

c.SetLogy()
c.SaveAs('ExcitedQuarksToGJ_AllCouplings_ObseExp_xsAccEff_Limits_TEMP.pdf')
#c.SaveAs('ExcitedQuarksToGJ_AllCouplings_ObseExp_xsAccEff_Limits.eps')
#c.SaveAs('ExcitedQuarksToGJ_AllCouplings_ObseExp_xsAccEff_Limits.png')
#c.SaveAs('SigmaBR_ReReco_Pt170_M560GeV_dEta2p0_ObseExp_Limits_qstar_AllCouplings.pdf')
#c.SaveAs('SigmaBR_ReReco_Pt170_M560GeV_dEta2p0_ObseExp_Limits_qstar_AllCouplings.gif')


