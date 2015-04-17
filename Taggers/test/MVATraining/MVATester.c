{

//Set values
	int events = 1000;

	TString events_string;
	events_string.Form("%d",events);
	TString Level  = "VBFDiPhoDiJet";
	TString path   = "${WORKSPACE}/test_diphodijet_training/";

	const float signal_XS[4] = {43.92,3.748,2.2496,0.5608};
	const float signal_BR[4] = {2.28e-3,2.28e-3,2.28e-3,2.28e-3};
	const TString signal_Names[4] = {"ggH","VBF","VH","ttH"};

	const float background_BR[3] = {1.0,0.0379,0.001776};
	const float background_XS[3] = {4746.0,17180.0,145400.0};

//Read in trees
	TFile *inputB1 = TFile::Open(path + "output_DYJetsToLL_M-50_13TeV-madgraph-pythia8_numEvent"+events_string+"_histos.root");
	TFile *inputB2 = TFile::Open(path + "output_GJet_Pt40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent"+events_string+"_histos.root");
	TFile *inputB3 = TFile::Open(path + "output_GJet_Pt20to40_doubleEMEnriched_TuneZ2star_13TeV-pythia6_numEvent"+events_string+"_histos.root");
	TFile *inputS1 = TFile::Open(path + "output_GluGluToHToGG_M-125_13TeV-powheg-pythia6_numEvent"+events_string+"_histos.root");
	TFile *inputS2 = TFile::Open(path + "output_VBF_HToGG_M-125_13TeV-powheg-pythia6_numEvent"+events_string+"_histos.root");	

	TTree *treeB1 = (TTree*)inputB1->Get(Level +"MVADumperNew/trees/dyJets_13TeV_All");
    	TTree *treeB2 = (TTree*)inputB2->Get(Level +"MVADumperNew/trees/gamJet_13TeV_All");
    	TTree *treeB3 = (TTree*)inputB3->Get(Level +"MVADumperNew/trees/gamJet_13TeV_All");
    	TTree *treeS1 = (TTree*)inputS1->Get(Level +"MVADumperNew/trees/ggh_m125_13TeV_All");
	TTree *treeS2 = (TTree*)inputS2->Get(Level +"MVADumperNew/trees/vbf_m125_13TeV_All");

//Set weights
	float weights_B[3];
	float weights_S[4];

     //Background weights
	weights_B[0] = 1;
	weights_B[1] = (events/(background_XS[0]*background_BR[0]))/((events/(background_XS[1]*background_BR[1])));
	weights_B[2] = (events/(background_XS[0]*background_BR[0]))/((events/(background_XS[2]*background_BR[2])));

     //Signal weights
	weights_S[0] = (events/(background_XS[0]*background_BR[0]))/((events/(signal_XS[0]*signal_BR[0])));
	weights_S[1] = (events/(background_XS[0]*background_BR[0]))/((events/(signal_XS[1]*signal_BR[1])));
	weights_S[2] = (events/(background_XS[0]*background_BR[0]))/((events/(signal_XS[2]*signal_BR[2])));
	weights_S[3] = (events/(background_XS[0]*background_BR[0]))/((events/(signal_XS[3]*signal_BR[3])));

//Set up histograms
	TH1F * histB1 = new TH1F("histB1","Drell Yan m_{gg}",50,100,150);
	TH1F * histB2 = new TH1F("histB2","GJets P_{t}:40 m_{gg}",50,100,150);
	TH1F * histB3 = new TH1F("histB3","GJets P_{t}:20 to 40 m_{gg}",50,100,150);
	TH1F * histS1 = new TH1F("histS1","ggH m_{gg}",50,100,150);
	TH1F * histS2 = new TH1F("histS2","VBF m_{gg}",50,100,150);

//Canvas prep
	TCanvas c1("c1");
	gStyle.SetOptStat(111111);
	gStyle.SetOptFit(111);

//Extract diphoton masses from trees
	float mgg;
	int numEntries;
	
	//B1: Drell-Yan
	numEntries = (int)treeB1->GetEntries();
	treeB1->SetBranchAddress("mgg",&mgg);
	for (int i=0;i<numEntries;i++) {
		treeB1->GetEntry(i);
		histB1->Fill(mgg);
	}
	//B2: GJ Pt40
	numEntries = (int)treeB2->GetEntries();
	treeB2->SetBranchAddress("mgg",&mgg);
	for (int i=0;i<numEntries;i++) {
		treeB2->GetEntry(i);
		histB2->Fill(mgg);
	}
	//B3: GJ Pt20-40
	numEntries = (int)treeB3->GetEntries();
	treeB3->SetBranchAddress("mgg",&mgg);
	for (int i=0;i<numEntries;i++) {
		treeB3->GetEntry(i);
		histB3->Fill(mgg);
	}
	//S1: ggH
	numEntries = (int)treeS1->GetEntries();
	treeS1->SetBranchAddress("mgg",&mgg);
	for (int i=0;i<numEntries;i++) {
		treeS1->GetEntry(i);
		histS1->Fill(mgg);
	}
	//S2: VBF
	numEntries = (int)treeS2->GetEntries();
	treeS2->SetBranchAddress("mgg",&mgg);
	for (int i=0;i<numEntries;i++) {
		treeS2->GetEntry(i);
		histS2->Fill(mgg);
	}

//Apply weightings
	histB1->Scale(weights_B[0]);
	histB2->Scale(weights_B[1]);
	histB3->Scale(weights_B[2]);
	histS1->Scale(weights_S[0]);
	histS2->Scale(weights_S[1]);

//Aggregate signal, background, and signal+background histograms
	TH1F * histB  = new TH1F("histB","Total Background",50,100,150);
	TH1F * histS  = new TH1F("histS","Total Signal",50,100,150);
	TH1F * histSB1 = new TH1F("histSB1","Signal + Background ggH",50,100,150); 
	TH1F * histSB2 = new TH1F("histSB2","Signal + Background VBF",50,100,150); 
	TH1F * histSBS = new TH1F("histSBS","Signal + Background ggH+VBF",50,100,150); 

	histB->Add(histB1);
	histB->Add(histB2);
	histB->Add(histB3);
	histS->Add(histS1);
	histS->Add(histS2);
	histSB1->Add(histB);
	histSB1->Add(histS1);
	histSB2->Add(histB);
	histSB2->Add(histS2);
	histSBS->Add(histB);
	histSBS->Add(histS);

//Fits
	float SBRatio1;
	float SBRatio2;
	float SBRatioS;

	//Background Laurent
	TF1 * laurent = new TF1("laurent","[0]/pow(x,4) + [1]/pow(x,5)",0,4);
	laurent->SetParameter(0, 1e8);
	laurent->SetParameter(1, -1e8);
	histB->Fit("laurent","QS");
	TFitResultPtr s = histB->Fit("laurent","QS");
	
	//Signal Gaussians
	TFitResultPtr r1 = histS1->Fit("gaus","QS","",120,130);
	TFitResultPtr r2 = histS2->Fit("gaus","QS","",120,130);
	TFitResultPtr rS = histS2->Fit("gaus","QS","",120,130);

	//Simultaneous G+L fits
	if((!r1->IsEmpty()) && (!s->IsEmpty()) && (!r2->IsEmpty())) {

		TF1 * Sim_Fit = new TF1("Laurent_Gauss","[0]/pow(x,4) + [1]/pow(x,5) + [2]*exp(-pow((x-[3]),2)/(2*[4]))",0,4);
		//ggH
		Sim_Fit->SetParName(0,"p0");
		Sim_Fit->SetParameter(0, s->Parameter(0));
		Sim_Fit->SetParName(1,"p1");
		Sim_Fit->SetParameter(1, s->Parameter(1));
		Sim_Fit->SetParName(2,"Ampl.");
		Sim_Fit->SetParameter(2, r1->Parameter(0));
		Sim_Fit->SetParName(3,"Mean");
		Sim_Fit->SetParameter(3, r1->Parameter(1));
		Sim_Fit->SetParName(4,"Std Dev.");
		Sim_Fit->SetParameter(4, r1->Parameter(2));
		TFitResultPtr t1 = histSB1->Fit("Laurent_Gauss","QS");
		//VBF
		Sim_Fit->SetParName(0,"p0");
		Sim_Fit->SetParameter(0, s->Parameter(0));
		Sim_Fit->SetParName(1,"p1");
		Sim_Fit->SetParameter(1, s->Parameter(1));
		Sim_Fit->SetParName(2,"Ampl.");
		Sim_Fit->SetParameter(2, r2->Parameter(0));
		Sim_Fit->SetParName(3,"Mean");
		Sim_Fit->SetParameter(3, r2->Parameter(1));
		Sim_Fit->SetParName(4,"Std Dev.");
		Sim_Fit->SetParameter(4, r2->Parameter(2));
		TFitResultPtr t2 = histSB2->Fit("Laurent_Gauss","QS");
		//ggH + VBF
		Sim_Fit->SetParName(0,"p0");
		Sim_Fit->SetParameter(0, s->Parameter(0));
		Sim_Fit->SetParName(1,"p1");
		Sim_Fit->SetParameter(1, s->Parameter(1));
		Sim_Fit->SetParName(2,"Ampl.");
		Sim_Fit->SetParameter(2, rS->Parameter(0));
		Sim_Fit->SetParName(3,"Mean");
		Sim_Fit->SetParameter(3, rS->Parameter(1));
		Sim_Fit->SetParName(4,"Std Dev.");
		Sim_Fit->SetParameter(4, rS->Parameter(2));
		TFitResultPtr tS = histSBS->Fit("Laurent_Gauss","QS");

		//S/(S+B)
		int maxBin1 = histS1->GetMaximumBin();
		int maxBin2 = histS2->GetMaximumBin();
		int maxBinS = histS->GetMaximumBin();
		
		float numS1 = histS1->Integral(maxBin1 - r1->Parameter(2)/1.25, maxBin1 + r1->Parameter(2)/1.25);
		float numS2 = histS2->Integral(maxBin2 - r2->Parameter(2)/1.25, maxBin2 + r2->Parameter(2)/1.25);
		float numS  = histS->Integral(maxBinS - rS->Parameter(2)/1.25, maxBinS + rS->Parameter(2)/1.25);
		float numB1 = histB->Integral(maxBin1 - r1->Parameter(2)/1.25, maxBin1 + r1->Parameter(2)/1.25);
		float numB2 = histB->Integral(maxBin2 - r2->Parameter(2)/1.25, maxBin2 + r2->Parameter(2)/1.25);
		float numBS = histB->Integral(maxBinS - rS->Parameter(2)/1.25, maxBinS + rS->Parameter(2)/1.25);
		
		SBRatio1 = numS1/sqrt(numS1+numB1);
		SBRatio2 = numS2/sqrt(numS2+numB2);
		SBRatioS = numS/sqrt(numS+numBS);
	
		std::cout << setw(12) << numS1;
		std::cout << setw(12) << numS2;
		std::cout << setw(12) << numS;
		std::cout << setw(12) << numB1;
		std::cout << setw(12) << numB2;
		std::cout << setw(12) << numBS;
		std::cout << std::endl;

		std::cout << "ggH: " << setw(12) << SBRatio1;
		std::cout << " VBF: " << setw(12) << SBRatio2;
		std::cout << " ggH + VBF: " << setw(12) << SBRatioS;
		std::cout << std::endl;

		//pdf Plotting output
		//ggH
		histSB1->Draw();
		c1.SaveAs("ggH.pdf(");
		histB->SetLineColor(kRed);
		histB->Draw("same");
		c1.SaveAs("ggH.pdf");
		histB->SetLineColor(kBlue);
		histB->Draw();
		c1.SaveAs("ggH.pdf");
		histS1->Draw();
		c1.SaveAs("ggH.pdf)");
		//VBF
		histSB2->Draw();
		c1.SaveAs("VBF.pdf(");
		histB->SetLineColor(kRed);
		histB->Draw("same");
		c1.SaveAs("VBF.pdf");
		histB->SetLineColor(kBlue);
		histB->Draw();
		c1.SaveAs("VBF.pdf");
		histS2->Draw();
		c1.SaveAs("VBF.pdf)");
		//ggH + VBF
		histSBS->Draw();
		c1.SaveAs("ggH_VBF.pdf(");
		histB->SetLineColor(kRed);
		histB->Draw("same");
		c1.SaveAs("ggH_VBF.pdf");
		histB->SetLineColor(kBlue);
		histB->Draw();
		c1.SaveAs("ggH_VBF.pdf");
		histS->Draw();
		c1.SaveAs("ggH_VBF.pdf)");

	}















}
