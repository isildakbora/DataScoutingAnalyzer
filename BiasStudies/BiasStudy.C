void BiasStudy()
{
	/// Dijet Mass Binning ///
	int number_of_variableWidth_bins = 88 - 1;
	Double_t massBins[88] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250,7500,7750,8000}; 
	Double_t mass_low  = 354.;
	Double_t mass_high = 5455;

	/// Get Fit Functions ///
	TF1 *func_type1;
	TFile *file_type1 = TFile::Open("../FisherStudies/dijetFitResults_FuncType1_nParFit5_Run2012BCD.root");
	func_type1 = (TF1*)file_type1->FindObjectAny("M1Bkg")->Clone("M1Bkg_func_type_1");
	func_type1->SetRange(mass_low, mass_high);
	func_type1->SetLineColor(kRed);
	func_type1->SetNpx(10e+4);

	TF1 *func_type2;
	TFile *file_type2 = TFile::Open("../FisherStudies/dijetFitResults_FuncType5_nParFit6_Run2012BCD.root");
	func_type2 = (TF1*)file_type2->FindObjectAny("M1Bkg")->Clone("M1Bkg_func_type_2");
	func_type2->SetRange(mass_low, mass_high);
	func_type2->SetLineColor(kBlue);
	func_type2->SetNpx(1e+4);

	/// Get Generated Spectra ///
	TFile *file_gen_spectra =  TFile::Open("generated_spectra.root");
	TH1F *gen_type1 = (TH1F*)file_gen_spectra->FindObjectAny("Generated_by_Type_1_Function");
	TH1F *gen_type2 = (TH1F*)file_gen_spectra->FindObjectAny("Generated_by_Type_2_Function");

	for(int fit_loop=0; fit_loop < 20; ++ fit_loop)
	{
		gen_type1->Fit("M1Bkg_func_type_2", "LR");
	}
	TH1F* gen_type1_pull = new TH1F("Generated_by_Type_1_Pull","Generated_by_Type_1_Pull", 20, -5, 5);

	for (int i = 0; i < gen_type1->GetNbinsX(); ++i)
	{	
		if(gen_type1->GetBinLowEdge(i+1) >= mass_low && (gen_type1->GetBinLowEdge(i+1) + gen_type1->GetBinWidth(i+1)) <= mass_high)
		{
			Double_t pseudo_data = gen_type1->GetBinContent(i+1);
			Double_t fit = func_type2->Integral(gen_type1->GetBinLowEdge(i+1), gen_type1->GetBinLowEdge(i+1) + gen_type1->GetBinWidth(i+1));
			cout << gen_type1->GetBinLowEdge(i+1)<< " " << gen_type1->GetBinLowEdge(i+1) + gen_type1->GetBinWidth(i+1) << endl;
			fit = fit/ gen_type1->GetBinWidth(i+1);
			//cout << gen_type1->GetBinLowEdge(i+1)<< " " << pseudo_data << "  "<<fit << " " << gen_type1->GetBinError(i+1) << " "<<(pseudo_data-fit)/gen_type1->GetBinError(i+1) << endl;
			gen_type1_pull->Fill((pseudo_data-fit)/gen_type1->GetBinError(i+1));
		}
	}


	TCanvas *Canvas0 = new TCanvas("Canvas0","Canvas0");
	gen_type1->Draw("E");
	gen_type1->GetXaxis()->SetRangeUser(mass_low-100, mass_high+100);
	gen_type1->GetXaxis()->SetMoreLogLabels();
	gen_type1->GetXaxis()->SetNoExponent();
	gPad->SetLogx();
	gPad->SetLogy();

	TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1");
	gen_type1_pull->Fit("gaus","L","",-3,3);
	gen_type1_pull->Draw();
}