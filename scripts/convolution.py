#!usr/bin/python
import sys, getopt
import ROOT
from ROOT import TFile, TH1F, TCanvas, TPad, TRandom3, gRandom
from ROOT import gROOT, gPad 
from ROOT import *
from rootutils import *
import math, sys, numpy as np
from bisect import bisect_left

print len(sys.argv)
if len(sys.argv) <= 2:
	sys.exit("Usage: python convolution2.py M(mass) N(number of pseudo experiments)")

massBins = array('d',[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250,7500,7750,8000])

res_bins_low  = np.arange(0., 1400., 200.)
res_bins_high = np.arange(1400., 4500., 500.)
res_bins      = np.concatenate([res_bins_low, res_bins_high])

mass = str(sys.argv[1])
N = int(float(sys.argv[2]))

resolution_file = TFile.Open('file:////Users/boraisildak/Desktop/tidy/DSComp/HLT_RECO_Smearing_Functions.root')

resolution_function = []

for i in xrange(len(res_bins)-1):
	fname = 'f_res_'+str(int(res_bins[i]))+'_'+str(int(res_bins[i+1]))
	func = resolution_file.FindObjectAny(fname)
	resolution_function.append(func)

file_to_write      = TFile('RSGravitonToQQbar_signals.root', 'UPDATE')
input_file         = TFile.Open('RSGravitonToQQbar_M_'+str(mass)+'.root')
input_histo        = input_file.FindObjectAny('h1_MjjWide_finalSel')
input_histo_varbin = input_file.FindObjectAny('h1_MjjWide_finalSel_varbin')
aux_histo          = input_histo.Clone('aux_histo')
aux_histo.Reset()

for i in xrange(int(input_histo.GetEntries())):
	aux_histo.Fill(input_histo.GetRandom())

can1 = TCanvas("mycanvas1","mycanvas1",600,600)
input_histo.Draw()
aux_histo.Draw('SAME')

hist_true = input_histo.Clone('hist_true')
hist_true.Reset()

hist_reco = input_histo.Clone('hist_reco')
hist_reco.Reset()

gRandom.SetSeed(0)
for i in xrange(N):
	progress = math.floor(50*i/N)
	progressbar(progress)
	m = input_histo.GetRandom()
	if m > 4500:
		continue
	index = bisect_left(res_bins, m)
	smearing = resolution_function[index-1]
	hist_true.Fill(m)
	smearing_factor = smearing.GetRandom()
	hist_reco.Fill(m*smearing_factor)

	#print m, index, resolution_function[index-1].GetName(), smearing_factor, m*smearing_factor

can2 = TCanvas("mycanvas2","mycanvas2",600,600)


hist_true.Scale(int(input_histo.GetEntries())/hist_true.Integral())
hist_true_varbin = hist_true.Rebin(len(massBins)-1, "hist_true_varbin", massBins)
hist_true_varbin.Draw()

hist_reco.Scale(int(input_histo.GetEntries())/hist_reco.Integral())
hist_reco_varbin = hist_reco.Rebin(len(massBins)-1, "hist_reco_varbin", massBins)
hist_reco_varbin.Draw("SAME")
hist_reco_varbin.SetLineColor(ROOT.kRed)

input_histo_varbin.Draw('SAME')
input_histo_varbin.SetLineColor(ROOT.kGreen)

can1.cd()
hist_reco_varbin.Draw("SAME")

file_to_write.cd()

input_histo.SetName("hist_true"+str(mass))
input_histo_varbin.SetName("hist_true_varbin_"+str(mass))

hist_reco.SetName("hist_smeared_"+str(mass))
hist_reco_varbin.SetName("hist_smeared_varbin_"+str(mass))

input_histo.Write()
input_histo_varbin.Write()

hist_reco.Write()
hist_reco_varbin.Write()

file_to_write.Close()

keepGUIalive()
