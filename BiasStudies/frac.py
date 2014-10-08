#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import *
from rootutils import *
import math

#Dijet Mass Binning 
massBins = array('d',[220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509])
mass_low  = 354.
mass_high = 4500.
N = mass_high-mass_low

#Display
gStyle.SetOptStat(0)
canvas = TCanvas("c", "FractionFitter example", 700, 700)
canvas.Divide(2,2)

file_type1 = TFile.Open("../FisherStudies/dijetFitResults_FuncType5_nParFit6_Run2012BCD.root")
func_type1 = file_type1.FindObjectAny("M1Bkg").Clone("M1Bkg_func_type_1")
gen_type1 = TH1F("Generated_by_Type_1_Function","Generated_by_Type_1_Function", int(N), mass_low, mass_high)

background = TH1F("background","background", int(N), mass_low, mass_high)


#Generate by type 1

for i in xrange(0, int(N)):
	r = TRandom1()
	mjj_y = func_type1.Integral(mass_low+i, mass_low+i+1)
	gen_type1.SetBinContent(i+1, r.Poisson(mjj_y))
	gen_type1.SetBinError(i+1, TMath.Sqrt(gen_type1.GetBinContent(i+1)))
	print mjj_y, gen_type1.GetBinLowEdge(i+1), gen_type1.GetBinContent(i+1)
gen_type1.Rebin(1)

#Generate background
for i in xrange(0, int(N)):
	r = TRandom1()
	mjj_y = func_type1.Integral(mass_low+i, mass_low+i+1)
	background.SetBinContent(i+1, r.Poisson(mjj_y))
	background.SetBinError(i+1, TMath.Sqrt(background.GetBinContent(i+1)))
	print mjj_y, background.GetBinLowEdge(i+1), background.GetBinContent(i+1)
background.Rebin(1)

mc_file = TFile.Open("/Users/boraisildak/Dropbox/DataScouting/test/RSGravitonToGG_signals.root")
mc_hist = mc_file.FindObjectAny("hist_true400")
mc_hist.Rebin(2)

canvas.cd(1)
mc_hist.Draw()
canvas.cd(2)
background.Draw()

canvas.cd(1)
gen_type1.SetMarkerColor(ROOT.kBlue)
gen_type1.SetLineColor(ROOT.kBlue)
gen_type1.Draw('SAME')

#FractionFitter
mc = TObjArray(2)
mc.Add(background)
mc.Add(mc_hist) 

fit =  TFractionFitter(gen_type1, mc)
fit.Constrain(0, 0.0, 1.0)
fit.Constrain(1, 0.0, 1.0)                # constrain fraction 1 to be between 0 and 1
#fit.SetRangeX(3,13)
status = fit.Fit()	                   # perform the fit
print "fit status: ",  status
if status == 0:
	p0 = ROOT.Double()
	errP0 = ROOT.Double()
	p1 = ROOT.Double()
	errP1 = ROOT.Double()
	fit.GetResult( 0, p0, errP0)
	fit.GetResult( 1, p1, errP1)
	print p0, p1, p0+p1
	canvas.cd(3)
	result = fit.GetPlot()
	result.Draw()
	gen_type1.Draw("Epsame")
	mc_hist.Draw("Epsame")
keepGUIalive()

