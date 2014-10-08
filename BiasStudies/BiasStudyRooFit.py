#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import *
from rootutils import *
import math

#Dijet Mass Binning 
#massBins = array('d',[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509])
massBins = array('d', range(1, 4509))
mass_low  = 354.
mass_high = 4500.
N = mass_high-mass_low

file_type1 = TFile.Open("../FisherStudies/dijetFitResults_FuncType5_nParFit6_Run2012BCD.root")
func_type1 = file_type1.FindObjectAny("M1Bkg").Clone("M1Bkg_func_type_1")

gen_type1  = TH1F("Generated_by_Type_1_Function","Generated_by_Type_1_Function", len(massBins)-1, massBins)
gen_type1.Sumw2()



#Generate by type 1
for i in xrange(0, len(massBins)-1):
	
	r         = TRandom1()
	mjj_y     = func_type1.Integral(massBins[i], massBins[i+1])
	bin_width = massBins[i+1] - massBins[i]

	gen_type1.SetBinContent(i+1, r.Poisson(mjj_y))
	gen_type1.SetBinError(i+1, TMath.Sqrt(gen_type1.GetBinContent(i+1)))
	print mjj_y, gen_type1.GetBinLowEdge(i+1), gen_type1.GetBinContent(i+1)

gen_type1.Draw("SAME")
#gen_type1.Rebin(20)

#gen_type1.Rebin(87, "hist_mass_varbin", massBins)

can = TCanvas('can', 'can', 600, 600)

func_type1.Draw()

x = RooRealVar("mjj", "mjj", 354., 4500.)
#x.setBinning(RooUniformBinning(354., 4500., int((4500-354)/2)), "cache")

m = RooRealVar("mean", "mean", 400.)
s = RooRealVar("sigma", "sigma", 20.)
a = RooRealVar("alpha", "alpha", 1)
n = RooRealVar("n", "n", 1)

frac      = RooRealVar("frac" ,"Signal fraction" , 0.1 , 0.0, 1.0)
true_frac = RooRealVar("true_frac" ,"True Signal fraction" , 0.5)
hist1     = RooDataHist ("hist1","hist1", RooArgList(x), gen_type1)
#hist2 = RooDataHist ("hist2","hist2", RooArgList(x), hist_mass_varbin)

scale = 100.
p1 = RooRealVar('p1','p1', func_type1.GetParameter(1), func_type1.GetParameter(1)-scale*func_type1.GetParError(1), func_type1.GetParameter(1)+scale*func_type1.GetParError(1)) 
p2 = RooRealVar('p2','p2', func_type1.GetParameter(2), func_type1.GetParameter(2)-scale*func_type1.GetParError(2), func_type1.GetParameter(2)+scale*func_type1.GetParError(2))
p3 = RooRealVar('p3','p3', func_type1.GetParameter(3), func_type1.GetParameter(3)-scale*func_type1.GetParError(3), func_type1.GetParameter(3)+scale*func_type1.GetParError(3))
p4 = RooRealVar('p4','p4', func_type1.GetParameter(4), func_type1.GetParameter(4)-scale*func_type1.GetParError(4), func_type1.GetParameter(4)+scale*func_type1.GetParError(4))
p5 = RooRealVar('p5','p5', func_type1.GetParameter(5), func_type1.GetParameter(5)-scale*func_type1.GetParError(5), func_type1.GetParameter(5)+func_type1.GetParError(5))
p0 = RooRealVar('p0','p0', func_type1.GetParameter(0), func_type1.GetParameter(0)-scale*func_type1.GetParError(0), func_type1.GetParameter(0)+scale*func_type1.GetParError(0))

#bkg = RooGenericPdf('background', '@1*(pow(1-@0/8000,@2)*(1+@5*(@0/8000)+@6*pow(@0/8000,2))/pow(@0/8000,@3+@4*log(@0/8000)))', RooArgList(x, p0, p1, p2, p3, p4, p5))

bkg = RooGenericPdf("background", "@1*(pow(1-@0/8000,@2))/pow(@0/8000,(@3+@4*log(@0/8000)+@5*pow(log(@0/8000),2)+@6*pow(log(@0/8000),3)))", RooArgList(x, p0, p1, p2, p3, p4, p5) )

pdf1  = RooHistPdf ("pdf1", "pdf1", RooArgSet(x), hist1)

sig   = RooCBShape("CB Shape", "Crystal Ball", x, m, s, n, a)
bkg_plus_sig = RooAddPdf("model","Signal + Background", sig, bkg, true_frac)
data = bkg_plus_sig.generate(RooArgSet(x), 1e+4)
model = RooAddPdf("model","Signal + Background", sig, bkg, frac)
pull = 0
if pull == 1:
	can2 = TCanvas('can2', 'can2', 600, 900)
	can2.cd(1).SetBottomMargin(0.4)
else:
	can2 = TCanvas('can2', 'can2', 600, 600)

rllist = RooLinkedList()
rllist.Add(RooFit.Save())
rllist.Add(RooFit.Strategy(2))
can2.SetLogx(True)
can2.SetLogy(True)
#r_chi2_wgt = bkg.chi2FitTo(hist1, rllist)
model.fitTo(data)
frame = x.frame()
frame2 = x.frame()
data.plotOn(frame)
model.plotOn(frame, RooFit.Precision(1))
frame.GetXaxis().SetRangeUser(354., 4500.)
#frame.GetYaxis().SetRangeUser(1e-2, 1e+7)
frame.GetXaxis().SetTitle('m_{jj} (GeV)')
frame.Draw()
frame.GetXaxis().UnZoom()


if pull == 1:
	hpull = frame.pullHist()
	frame2.addPlotable(hpull,'p')
	pad = TPad('pad','pad',0.,0.,1.,1.)
	pad.SetTopMargin(0.6)
	pad.SetFillColor(0)
	pad.SetFillStyle(0)
	pad.Draw()
	pad.cd(0)
	frame2.SetMinimum(-10)
	frame2.SetMaximum(10)
	frame2.GetYaxis().SetNdivisions(505)
	frame2.GetXaxis().SetTitleOffset(0.9)
	frame2.GetYaxis().SetTitleOffset(0.8)
	frame2.GetYaxis().SetTickLength(0.06)
	frame2.GetYaxis().SetTitleSize(0.05)
	frame2.GetYaxis().SetLabelSize(0.03)
	frame2.GetYaxis().SetTitle('(Data-Fit)/Error')
	frame2.GetXaxis().SetTitle('m_{jj} (GeV)')
	frame2.Draw()

for i in xrange(0,6):
	print func_type1.GetParameter(i)
print frac.getVal(), frac.getAsymErrorLo()
keepGUIalive()

