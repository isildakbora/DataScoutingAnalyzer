#!usr/bin/python
import ROOT
from ROOT import TFile, TF1, TH1F, TMath, TCanvas, TLegend, gROOT, gPad, TObject, TRandom1, TPaveText
from array import array

gROOT.ProcessLine(".X ~/setTDRStyle.C")
progress = 0

#Dijet Mass Binning 

massBins = array('d',[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250,7500,7750,8000])

mass_low  = 354
mass_high = 5455

#Get Fit Functions
file_type1 = TFile.Open("../FisherStudies/dijetFitResults_FuncType1_nParFit5_Run2012BCD.root")
func_type1 = file_type1.FindObjectAny("M1Bkg").Clone("M1Bkg_func_type_1")
func_type1.SetRange(mass_low, mass_high)
func_type1.SetLineColor(1)
func_type1.SetNpx(10000)

file_type2 = TFile.Open("../FisherStudies/dijetFitResults_FuncType5_nParFit6_Run2012BCD.root")
func_type2 = file_type2.FindObjectAny("M1Bkg").Clone("M1Bkg_func_type_2")
func_type2.SetRange(mass_low, mass_high)
func_type2.SetLineColor(2)
func_type2.SetNpx(10000)

#Get Generated Spectra
file_gen_spectra =  TFile.Open("generated_spectra.root")
gen_type1        = file_gen_spectra.FindObjectAny("Generated_by_Type_1_Function")
gen_type2        = file_gen_spectra.FindObjectAny("Generated_by_Type_2_Function")


for fit_loop in range(0,5):
	gen_type1.Fit("M1Bkg_func_type_1", "LR")
	gen_type2.Fit("M1Bkg_func_type_2", "LR")

gen_type1_pull = TH1F("Generated_by_Type_1_Pull","Generated_by_Type_1_Pull", 40, -5, 5)
gen_type2_pull = TH1F("Generated_by_Type_2_Pull","Generated_by_Type_2_Pull", 40, -5, 5)

for i in range(0,gen_type1.GetNbinsX()):

	if (gen_type1.GetBinLowEdge(i+1) >= mass_low and (gen_type1.GetBinLowEdge(i+1) + gen_type1.GetBinWidth(i+1)) <= mass_high):
		pseudo_data = gen_type1.GetBinContent(i+1)
		fit = func_type2.Integral(gen_type1.GetBinLowEdge(i+1), gen_type1.GetBinLowEdge(i+1) + gen_type1.GetBinWidth(i+1))
		#print gen_type1.GetBinLowEdge(i+1), gen_type1.GetBinLowEdge(i+1) + gen_type1.GetBinWidth(i+1)
		fit = fit/ gen_type1.GetBinWidth(i+1)
		if(gen_type1.GetBinError(i+1)>0):
			gen_type1_pull.Fill((pseudo_data-fit)/gen_type1.GetBinError(i+1))

	if (gen_type2.GetBinLowEdge(i+1) >= mass_low and (gen_type2.GetBinLowEdge(i+1) + gen_type2.GetBinWidth(i+1)) <= mass_high):
		pseudo_data = gen_type2.GetBinContent(i+1)
		fit = func_type1.Integral(gen_type2.GetBinLowEdge(i+1), gen_type2.GetBinLowEdge(i+1) + gen_type2.GetBinWidth(i+1))
		#print gen_type2.GetBinLowEdge(i+1), gen_type2.GetBinLowEdge(i+1) + gen_type2.GetBinWidth(i+1)
		fit = fit/ gen_type2.GetBinWidth(i+1)
		if(gen_type2.GetBinError(i+1)>0):
			gen_type2_pull.Fill((pseudo_data-fit)/gen_type2.GetBinError(i+1))


Canvas0 = TCanvas("Canvas0","Canvas0")
gen_type1.Draw("E")
gen_type1.GetXaxis().SetRangeUser(mass_low-100, mass_high+100)
gen_type1.GetXaxis().SetMoreLogLabels()
gen_type1.GetXaxis().SetNoExponent()
gPad.SetLogx()
gPad.SetLogy()

Canvas1 = TCanvas("Canvas1","Canvas1")
gen_type1_pull.Fit("gaus","L","",-3,3)
gen_type1_pull.Draw()

leg1 = TPaveText(0.2,0.80,0.5,0.85, 'NDC')
leg1.InsertText('pull between generated')
leg1.InsertText('with type 1 fit by type 1')
leg1.SetFillColor(0)
leg1.SetLineColor(0)
leg1.SetShadowColor(0)
leg1.SetTextFont(42)
leg1.SetTextColor(1)
leg1.SetTextSize(0.03)
leg1.Draw()
Canvas1.SaveAs('pull_between_generated_with_type_1_fit_by_type_1.pdf')
Canvas1.SaveAs('pull_between_generated_with_type_1_fit_by_type_1.png')

Canvas2 = TCanvas("Canvas2","Canvas2")
gen_type2.Draw("E")
gen_type2.GetXaxis().SetRangeUser(mass_low-100, mass_high+100)
gen_type2.GetXaxis().SetMoreLogLabels()
gen_type2.GetXaxis().SetNoExponent()
gPad.SetLogx()
gPad.SetLogy()

Canvas3 = TCanvas("Canvas3","Canvas3")
gen_type2_pull.Fit("gaus","L","",-3,3)
gen_type2_pull.Draw()

leg2 = TPaveText(0.2,0.80,0.5,0.85, 'NDC')
leg2.InsertText('pull between generated')
leg2.InsertText('with type 2 fit by type 2')
leg2.SetFillColor(0)
leg2.SetLineColor(0)
leg2.SetShadowColor(0)
leg2.SetTextFont(42)
leg2.SetTextColor(1)
leg2.SetTextSize(0.03)
leg2.Draw()
Canvas3.SaveAs('pull_between_generated_with_type_2_fit_by_type_2.pdf')
Canvas3.SaveAs('pull_between_generated_with_type_2_fit_by_type_2.png')

#----- keep the GUI alive ------------
if __name__ == '__main__':
  rep = ''
  while not rep in ['q','Q']:
    rep = raw_input('enter "q" to quit: ')
    if 1 < len(rep):
      rep = rep[0]