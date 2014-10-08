#!usr/bin/python
import ROOT
from ROOT import TFile, TF1, TH1F, TMath, TCanvas, TLegend, gROOT, gPad, TObject, TRandom1
from array import array

def update_progress(progress):
    print '\r[{0}] {1}%'.format('#'*(progress), progress)

gROOT.ProcessLine(".X ~/setTDRStyle.C")
progress = 0

#Dijet Mass Binning 
massBins = array('d',[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250,7500,7750,8000])

mass_low  = 354
mass_high = 5455
N = mass_high-mass_low
print N

#Get Fit Functions
file_type1 = TFile.Open("../FisherStudies/dijetFitResults_FuncType1_nParFit5_Run2012BCD.root")
func_type1 = file_type1.FindObjectAny("M1Bkg").Clone("M1Bkg_func_type_1")
func_type1.SetLineColor(1)
func_type1.SetNpx(5000)

file_type2 = TFile.Open("../FisherStudies/dijetFitResults_FuncType5_nParFit6_Run2012BCD.root")
func_type2 = file_type2.FindObjectAny("M1Bkg").Clone("M1Bkg_func_type_2")
func_type2.SetLineColor(2)
func_type2.SetNpx(5000)

#Generate PseudoExperiments
generated_spectra = TFile("generated_spectra.root", "RECREATE")

gen_type1 = TH1F("Generated_by_Type_1_Function","Generated_by_Type_1_Function", len(massBins)-1, massBins)
gen_type1 = TH1F("Generated_by_Type_1_Function","Generated_by_Type_1_Function", N, mass_low, mass_high)
gen_type1.Sumw2()

gen_type2 = TH1F("Generated_by_Type_2_Function","Generated_by_Type_2_Function", len(massBins)-1, massBins)
gen_type2 = TH1F("Generated_by_Type_2_Function","Generated_by_Type_2_Function", N, mass_low, mass_high)	
gen_type2.Sumw2()

#Generate by type 1
for i in range(0,N):
	progress = 100.0*i/(1.0*N)
	k = TMath.FloorNint(progress)
	update_progress(k)
	r = TRandom1()
	mjj_y = func_type1.Integral(mass_low+i, mass_low+i+1)
	gen_type1.SetBinContent(i+1, r.Poisson(mjj_y))
	gen_type1.SetBinError(i+1, TMath.Sqrt(gen_type1.GetBinContent(i+1)))
	print mjj_y, gen_type1.GetBinLowEdge(i+1), gen_type1.GetBinContent(i+1)
	mjj_y = func_type2.Integral(mass_low+i, mass_low+i+1)
	gen_type2.SetBinContent(i+1, r.Poisson(mjj_y))
	gen_type2.SetBinError(i+1, TMath.sqrt(gen_type2.GetBinContent(i+1)))

Canvas0 = TCanvas("Canvas0","Canvas0")
gen_type1.Draw("E")
func_type1.Draw("SAME")
gen_type1.GetXaxis().SetRangeUser(mass_low-100, mass_high+100)
gen_type1.GetXaxis().SetMoreLogLabels()
gen_type1.GetXaxis().SetNoExponent()
gPad.SetLogx()
gPad.SetLogy()

Canvas1 = TCanvas("Canvas1","Canvas1")
gen_type2.Draw("E")
func_type2.Draw("SAME")
gen_type2.GetXaxis().SetRangeUser(mass_low-100, mass_high+100)
gen_type2.GetXaxis().SetMoreLogLabels()
gen_type2.GetXaxis().SetNoExponent()
gPad.SetLogx()
gPad.SetLogy()

generated_spectra.cd()
gen_type1.Write()
gen_type2.Write()

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]
