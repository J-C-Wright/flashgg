
void RadiationRecoveryPlots() {

    TFile *inputFile  = TFile::Open("output_numEvent10000.root");
    TFile *outputFile = new TFile("PlotOutput.root","RECREATE");

    TTree *tree = (TTree*)inputFile->Get("vbfTagDumper/trees/vbfh_13TeV_VBFDiJet");

    TObjArray *branchNames = tree->GetListOfBranches();

    std::vector<TString> variableNames;
    variableNames.push_back("J1J2_dEta"); variableNames.push_back("J1J3_dEta"); variableNames.push_back("J2J3_dEta");
    variableNames.push_back("dPhi_J1_J2");variableNames.push_back("dPhi_J1_J3");variableNames.push_back("dPhi_J2_J3");
    variableNames.push_back("J1J2_dR");   variableNames.push_back("J1J3_dR");   variableNames.push_back("J2J3_dR");
    variableNames.push_back("dR_min_J13J23");
    variableNames.push_back("J1J2_mjj"); variableNames.push_back("J1J3_mjj"); variableNames.push_back("J2J3_mjj");

    std::vector<TString> variableXLabels(variableNames.size());
    variableXLabels[0] = "#Delta#eta_{J1J2}";variableXLabels[1] = "#Delta#eta_{J1J3}";variableXLabels[2] = "#Delta#eta_{J2J3}";
    variableXLabels[3] = "|#Delta#phi_{#gamma#gammaJ1J2}|";variableXLabels[4] = "|#Delta#phi_{#gamma#gammaJ1J3}|";variableXLabels[5] = "|#Delta#phi_{#gamma#gammaJ2J3}|";
    variableXLabels[6] = "#DeltaR_{J1J2}";variableXLabels[7] = "#DeltaR_{J1J3}";variableXLabels[8] = "#DeltaR_{J2J3}";
    variableXLabels[9] = "dR_{min}";
    variableXLabels[10] = "m_{J1J2}";variableXLabels[11] = "m_{J1J3}";variableXLabels[12] = "m_{J2J3}";

    unsigned numBins(100);
    std::vector<std::pair<float,float>> rangeVector(variableNames.size());
    std::vector<std::vector<TH1F*>> histograms(variableNames.size());
    std::vector<TString> categoryLabels(4);
    categoryLabels[0] = "FFF"; categoryLabels[1] = "JFF"; categoryLabels[2] = "JJF"; categoryLabels[3] = "JJJ";
    for (unsigned variable(0);variable<variableNames.size();variable++) {

        rangeVector[variable].first  = 999.;
        rangeVector[variable].second = 0.0;
        
        //find range
        for (unsigned event(0);event<tree->GetEntries();event++) {    
            tree->GetEntry(event);
            float value = (float)tree->GetBranch(variableNames[variable])->GetLeaf(variableNames[variable])->GetValue();
            if (value < rangeVector[variable].first  && value > -990) rangeVector[variable].first = value;
            if (value > rangeVector[variable].second && value > -990) rangeVector[variable].second = value;
        }

        std::cout << "The min is " << rangeVector[variable].first;
        std::cout << " and the max is " << rangeVector[variable].second;
        std::cout << " for variable " << variableNames[variable] << std::endl;

        std::vector<TH1F*> catHists(categoryLabels.size());
        for (unsigned category(0);category<catHists.size();category++) {
            catHists[category] = new TH1F(variableNames[variable] + TString("_") + categoryLabels[category], 
                                          variableNames[variable],
                                          numBins,rangeVector[variable].first,rangeVector[variable].second);
            catHists[category]->GetXaxis()->SetTitle(variableXLabels[variable]);
        }
        histograms[variable] = catHists;
    }

    for (unsigned variable(0);variable<variableNames.size();variable++) {
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            unsigned category = (unsigned)tree->GetBranch("numberOfMatches")->GetLeaf("numberOfMatches")->GetValue();
            float value       = (float)tree->GetBranch(variableNames[variable])->GetLeaf(variableNames[variable])->GetValue();
            bool isTrijet     = (bool)tree->GetBranch("has3Jet")->GetLeaf("has3Jet")->GetValue();
            if (value > -9990.0 && isTrijet) {
                histograms[variable][category]->Fill(value);
            }
        }
    }

    TCanvas c1( "c1" );
    gStyle->SetOptStat( 0 );
    for (unsigned variable(0);variable<variableNames.size();variable++) {

        std::cout << variableNames[variable] << std::endl;
        float peak(0); unsigned maxCat(0);           
        std::vector<TH1F*> scaledHists(histograms[variable].size());

        for (unsigned category(0);category<categoryLabels.size();category++) {
            if (category == 4) {histograms[variable][category]->SetLineColor(6);}
            else if (category == 2) {histograms[variable][category]->SetLineColor(8);}
            else {histograms[variable][category]->SetLineColor(category+1);}
            histograms[variable][category]->SetLineWidth(2);

            scaledHists[category] = histograms[variable][category];
            if (scaledHists[category]->Integral() > 0.0) {
                scaledHists[category]->Scale(1/histograms[variable][category]->Integral());
            }
        }

        for (unsigned category(0);category<categoryLabels.size();category++) {
            if (scaledHists[category]->GetMaximum() > peak) {
                peak   = scaledHists[category]->GetMaximum();
                maxCat = category; 
            }
        }

        scaledHists[maxCat]->Draw();
        for (unsigned category(0);category<categoryLabels.size();category++) {
            if (category != maxCat) {scaledHists[category]->Draw("same");}
        }
        c1.Print("Plots/" + TString(variableNames[variable]) + TString(".pdf"));
    }

    //Make a separate legend
    c1.Clear();
    TLegend *legend = new TLegend(0.0,0.0,0.13,0.25);
    legend->SetTextFont(2);
    legend->SetTextSize(0.04);
    for (unsigned category(0);category<categoryLabels.size();category++) {
        legend->AddEntry(histograms[0][category],categoryLabels[category],"l");
    }
    legend->Draw();
    c1.Print("Plots/Legend.pdf");

    std::cout << std::endl;
    for (unsigned category(0);category<4;category++) {
        std::cout << setw(6) << categoryLabels[category];
        unsigned numInCategory(0);
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            unsigned numMatches = (unsigned)tree->GetBranch("numberOfMatches")->GetLeaf("numberOfMatches")->GetValue();
            bool isTrijet     = (bool)tree->GetBranch("has3Jet")->GetLeaf("has3Jet")->GetValue();
            if (numMatches == category && isTrijet) {
                numInCategory++;
            }
        }
        std::cout << setw(12) << numInCategory << std::endl;
    }

    //Experiment with inv. mass change
    std::vector<TH1F*> invMassComparison(4);
    for (unsigned category(0);category<4;category++) {  
        invMassComparison[category] = new TH1F(TString("InvMassComp") + categoryLabels[category],"Fractional change in m_{jj}",50,-0.5,1.5);
    }

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        unsigned category = (unsigned)tree->GetBranch("numberOfMatches")->GetLeaf("numberOfMatches")->GetValue();
        bool isTrijet     = (bool)tree->GetBranch("has3Jet")->GetLeaf("has3Jet")->GetValue();
        if (isTrijet) {
            float mjj  = (float)tree->GetBranch("J1J2_mjj")->GetLeaf("J1J2_mjj")->GetValue(); 
            float mjjj = (float)tree->GetBranch("Mjjj")->GetLeaf("Mjjj")->GetValue(); 
            invMassComparison[category]->Fill((mjjj-mjj)/mjj); 
        }
    }   

    outputFile->cd();
    for (unsigned category(0);category<4;category++) { invMassComparison[category]->Write(); }    
    outputFile->Close();

}
