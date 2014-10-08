
void BiasStudy_Spectrum()
{
	/// Number of Events ///
	
	double progress = 0;

	/// Dijet Mass Binning ///
	int number_of_variableWidth_bins = 88 - 1;
	Double_t massBins[88] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649,  693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000, 7250,7500,7750,8000}; 
	Double_t mass_low  = 354.;
	Double_t mass_high = 5455.;
	int N = mass_high-mass_low;
	cout << N << endl;

	/// Get Fit Functions ///
	TF1 *func_type1;
	TFile *file_type1 = TFile::Open("../FisherStudies/dijetFitResults_FuncType1_nParFit5_Run2012BCD.root");
	func_type1 = (TF1*)file_type1->FindObjectAny("M1Bkg")->Clone("M1Bkg_func_type_1");
	func_type1->SetLineColor(kRed);
	func_type1->SetNpx(5000);

	TF1 *func_type2;
	TFile *file_type2 = TFile::Open("../FisherStudies/dijetFitResults_FuncType5_nParFit6_Run2012BCD.root");
	func_type2 = (TF1*)file_type2->FindObjectAny("M1Bkg")->Clone("M1Bkg_func_type_2");
	func_type2->SetLineColor(kBlue);
	func_type2->SetNpx(5000);

	/// Generate PseudoExperiments ///
	TFile* generated_spectra= new TFile("generated_spectra.root", "RECREATE");
	Double_t mjj;
	Double_t mjj_y;


	//TH1F *gen_type1 = new TH1F("Generated_by_Type_1_Function","Generated_by_Type_1_Function", number_of_variableWidth_bins, massBins);
	TH1F *gen_type1 = new TH1F("Generated_by_Type_1_Function","Generated_by_Type_1_Function", N, mass_low, mass_high);
	//gen_type1->Sumw2();

	//TH1F *gen_type2 = new TH1F("Generated_by_Type_2_Function","Generated_by_Type_2_Function", number_of_variableWidth_bins, massBins);
	TH1F *gen_type2 = new TH1F("Generated_by_Type_2_Function","Generated_by_Type_2_Function", N, mass_low, mass_high);	
	//gen_type2->Sumw2();

	/// Generate by type 1 ///
	for (int i = 0; i < N; ++i)
	{
		progress = 100.0*i/(1.0*N);
	    int k = TMath::FloorNint(progress);
	    std::cout << "\r" << "Generating by type 1: " << k << "% completed: " << "[" << std::string(k, '|') << std::string(99-k, ' ') << "]";
		std::cout.flush();
		TRandom1* r = new TRandom1();
		mjj_y = func_type1->Integral(mass_low+i, mass_low+i+1);
		//mjj_y = func_type1->Eval(mass_low+i+0.5);
		//mjj = GetRand(func_type1, mass_low, mass_high, r->Integer(1e+10));
		gen_type1->SetBinContent(i+1, r->Poisson(mjj_y));
		gen_type1->SetBinError(i+1, TMath::Sqrt(gen_type1->GetBinContent(i+1)));
		cout << mjj_y << "\t" <<gen_type1->GetBinLowEdge(i+1) << "\t"<< gen_type1->GetBinContent(i+1) << endl;
		//mjj = func_type1->GetRandom(2000., mass_high);
		//gen_type1->Fill(mjj);

		mjj_y = func_type2->Integral(mass_low+i, mass_low+i+1);
		//mjj = GetRand(func_type1, mass_low, mass_high, r->Integer(1e+10));
		gen_type2->SetBinContent(i+1, r->Poisson(mjj_y));
		gen_type2->SetBinError(i+1, TMath::Sqrt(gen_type2->GetBinContent(i+1)));
	}

	TCanvas *Canvas0 = new TCanvas("Canvas0","Canvas0");
	gen_type1->Draw("E");
	func_type1->Draw("SAME");
	gen_type1->GetXaxis()->SetRangeUser(mass_low-100, mass_high+100);
	gen_type1->GetXaxis()->SetMoreLogLabels();
	gen_type1->GetXaxis()->SetNoExponent();
	gPad->SetLogx();
	gPad->SetLogy();

	TCanvas *Canvas1 = new TCanvas("Canvas1","Canvas1");
	gen_type2->Draw("E");
	func_type2->Draw("SAME");
	gen_type2->GetXaxis()->SetRangeUser(mass_low-100, mass_high+100);
	gen_type2->GetXaxis()->SetMoreLogLabels();
	gen_type2->GetXaxis()->SetNoExponent();
	gPad->SetLogx();
	gPad->SetLogy();

	generated_spectra->cd();
	gen_type1->Write();
	gen_type2->Write();
}

Double_t GetRand(TF1 *func, Double_t X_MIN, Double_t X_MAX, int seed )
{
	Double_t Y_MAX;
	Y_MAX = func->GetMaximum();
	//std::cout << "Y_MAX=" << Y_MAX << std::endl;
	TRandom2* r = new TRandom3(seed);
	
	Double_t x =0;
	Double_t y   = r->Uniform();
	Double_t f_y = 3e+6;
	
	while(1)
	{
		x   = r->Uniform(X_MIN, X_MAX);
		//std::cout << "x=" << x << std::endl;
		y   = r->Uniform(0, Y_MAX);
		//std::cout << "y=" << y << std::endl;
		f_y = func->Eval(x);
		//std::cout << "f_y=" << f_y << std::endl;

		if(f_y > y)
			break;
	}
	//std::cout << "x budur=" << x << std::endl;
	return x;

}