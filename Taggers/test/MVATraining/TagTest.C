{

//Load effective sigma macro
gROOT->ProcessLine(".L effSigmaMacro.C");

//Names and values
const float signal_XS[4] = {43.92,3.748,2.2496,0.5608};
const float signal_BR[4] = {2.28e-3,2.28e-3,2.28e-3,2.28e-3};
const TString signal_Names[4] = {"ggH","VBF","VH","ttH"};

const float background_XS[2] = {0.0379,0.001776};
const float background_BR[2] = {17180.0,145400.0};
const TString background_Names[2] = {"GJp40","GJp20_40"};

const TString tag_Names[13] = {"untagged 0","untagged 1","untagged 2","untagged 3","untagged 4",
                                           "vbf 0","vbf 1","vbf 2","tthhad","tthlep","vhtight","vhloose","vhhad"}; 
//Canvas prep
	TCanvas c1("c1");
	gStyle.SetOptStat(111111);
	gStyle.SetOptFit(111);

//Output prep
	ofstream output;
	output.open("Output.csv");
	output << "Tag, Total Signal, ggH %, VBF %, VH %, ttH %, Sigma(Sim Fit), S/(S+B), Sigma(Sig Fit), Sigma(Effective)" << std::endl; 

	outputRoot = new TFile("Output.root","RECREATE");

//Read in .root files
	TFile *background_Files[2];
	background_Files[0] = new TFile(background_Names[0] + ".root","READ");
	background_Files[1] = new TFile(background_Names[1] + ".root","READ");

	TFile *signal_Files[4];
	signal_Files[0] = new TFile(signal_Names[0] + ".root","READ");
	signal_Files[1] = new TFile(signal_Names[1] + ".root","READ");
	signal_Files[2] = new TFile(signal_Names[2] + ".root","READ");
	signal_Files[3] = new TFile(signal_Names[3] + ".root","READ");

//Weightings
	float weights[6];
	
	//Backgrounds
	weights[0] = 1;

	int BG_Num1 = 20*background_XS[0]*background_BR[0]*1000;
	TVectorF *v = (TVectorF*)background_Files[1]->Get("properties");
	weights[1] = (BG_Num1/(background_XS[0]*background_BR[0]))/((*v)[0]/(background_XS[1]*background_BR[1]));

	//Signals
	float sig_Weight_Total = 0;
	for (int i=2;i<6;i++) {

		*v = (TVectorF*)signal_Files[i-2]->Get("properties");
		weights[i] = (BG_Num1/(background_XS[0]*background_BR[0]))/((*v)[0]/(signal_XS[i-2]*signal_BR[i-2]));
		sig_Weight_Total += weights[i];

	}

	

//Fractions and fits
int count;
int signal_Numbers[4];
for (int i=0;i<13;i++) {

	output << tag_Names[i] << ",";
	
	//Signal Fractions
	count = 0;
	for (int j=0;j<4;j++) { 	
	
		signal_Hist = (TH1F*)signal_Files[j].Get(tag_Names[i]);
		signal_Numbers[j] = signal_Hist->GetEntries()*(weights[j+2]/sig_Weight_Total)*100;
		count += signal_Hist->GetEntries()*(weights[j+2]/sig_Weight_Total);

	}
	
	if (count != 0) {

		output << count << ",";	
		output << signal_Numbers[0]/(float)count << ",";
		output << signal_Numbers[1]/(float)count << ",";
		output << signal_Numbers[2]/(float)count << ",";
		output << signal_Numbers[3]/(float)count << ",";

	

	//Signal plus background fits
	TH1F * background_Hists[2];	
	background_Hists[0] = (TH1F*)background_Files[0].Get(tag_Names[i]);
	background_Hists[1] = (TH1F*)background_Files[1].Get(tag_Names[i]);
	TH1F * full_BG_Hist = new TH1F(tag_Names[i] + " B",tag_Names[i] + " Background",50,100,150);
	
	background_Hists[0]->Scale(weights[0]);
	background_Hists[1]->Scale(weights[1]);

	full_BG_Hist->Add(background_Hists[0]);
	full_BG_Hist->Add(background_Hists[1]);

	TH1F * full_S_Hist = new TH1F(tag_Names[i] + " S", tag_Names[i] + " Signal",50,100,150);

	float sigma;
	float SB_Ratio;

	TH1F * signal_Hist;
	for (int j=0;j<4;j++) { 	

		signal_Hist = (TH1F*)signal_Files[j].Get(tag_Names[i]);	
		signal_Hist->Scale(weights[j+2]);
		full_S_Hist->Add(signal_Hist);
	
	}

	TH1F * SB_Hist = new TH1F(tag_Names[i] + " S+B",tag_Names[i] + " S+B",50,100,150);
	SB_Hist->Add(full_BG_Hist);
	SB_Hist->Add(full_S_Hist);
	
	//Fits
	//Signal Gaussian
	TFitResultPtr r = full_S_Hist->Fit("gaus","QS","",120,130);

	//Background Laurent
	TF1 * laurent = new TF1("laurent","[0]/pow(x,4) + [1]/pow(x,5)",0,4);
	laurent->SetParameter(0, 1e8);
	laurent->SetParameter(1, -1e8);
	full_BG_Hist->Fit("laurent","QS");
	TFitResultPtr s = full_BG_Hist->Fit("laurent","QS");



	if (!s->IsEmpty()) {
		if (fabs(s->Parameter(0)) < s->ParError(0) || fabs(s->Parameter(1)) < s->ParError(1)) {
			std::cout << tag_Names[i] +  ": poor background fit" << std::endl;
		}
	}



	//Simultaneous Laurent + Gaussian
	if((!r->IsEmpty()) && (!s->IsEmpty())) {
			
		TF1 * Sim_Fit = new TF1("Laurent_Gauss","[0]/pow(x,4) + [1]/pow(x,5) + [2]*exp(-pow((x-[3]),2)/(2*[4]))",0,4);
		Sim_Fit->SetParName(0,"p0");
		Sim_Fit->SetParameter(0, s->Parameter(0));
		Sim_Fit->SetParName(1,"p1");
		Sim_Fit->SetParameter(1, s->Parameter(1));
		Sim_Fit->SetParName(2,"Ampl.");
		Sim_Fit->SetParameter(2, r->Parameter(0));
		Sim_Fit->SetParName(3,"Mean");
		Sim_Fit->SetParameter(3, r->Parameter(1));
		Sim_Fit->SetParName(4,"Std Dev.");
		Sim_Fit->SetParameter(4, r->Parameter(2));
		TFitResultPtr t = SB_Hist->Fit("Laurent_Gauss","QS");
		
		
		//.pdf export
		TString pdf_Name = tag_Names[i] + ".pdf";

		SB_Hist->Draw();
		c1.Print(pdf_Name + "(");

		full_BG_Hist->SetLineColor(kRed);
		full_BG_Hist->Draw("same");
		c1.Print(pdf_Name);
			
		full_BG_Hist->SetLineColor(kBlue);
		full_BG_Hist->Draw();
		c1.Print(pdf_Name);

		full_S_Hist->Draw();
		c1.Print(pdf_Name + ")");
 
 		//Sigma and S/(S+B)
		int maxBin = signal_Hist->GetMaximumBin();
		int S = full_S_Hist->Integral(maxBin - r->Parameter(2)/1.25, maxBin + r->Parameter(2)/1.25);
		int B = full_BG_Hist->Integral(maxBin - r->Parameter(2)/1.25, maxBin + r->Parameter(2)/1.25);
			
		sigma = t->Parameter(4);	
		SB_Ratio = (float)S/(float)(S+B);

		if (!t->IsEmpty()) {
			if (t->Parameter(2) < t->ParError(2) || TMath::IsNaN(t->ParError(2)) == 1 ) {
				
				std::cout << tag_Names[i] << ": poor simultaneous fit" << std::endl;
			}
		}
		

	} else if ((!r->IsEmpty()) && (s->IsEmpty())) {
			
		TString pdf_Name = tag_Names[i] + ".pdf";
		full_S_Hist->Draw();
		c1.Print(pdf_Name);
		sigma = r->Parameter(2);
		SB_Ratio = 1;

	} else {

		std::cout << signal_Names[j] + " " + tag_Names[i] + " is empty" << std::endl;

	}

	output << sigma << "," << SB_Ratio << "," << r->Parameter(2) << "," << effSigmaMacro(full_S_Hist) << std::endl; 

	}else{

	output << ",,,,," <<  std::endl;

	}
	//.ROOT export
	outputRoot->cd();
	SB_Hist->Write();
	full_BG_Hist->Write();
	full_S_Hist->Write();
	
}

c1.Close();
outputRoot->Close();
background_Files[0]->Close();
background_Files[1]->Close();
signal_Files[0]->Close();
signal_Files[1]->Close();
signal_Files[2]->Close();
signal_Files[3]->Close();

}
