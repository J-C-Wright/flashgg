
TGraph* getMergeVsDRCut(TTree *tree, unsigned colour) {
    
    Double_t mergeRate[17];
    Double_t dRCut[17];

    std::cout << setw(12) << "\% Merged" << setw(12) << "dR Cut" << std::endl;
    for (unsigned cut(0);cut<17;cut++) {
        unsigned count(0);
        dRCut[cut] = 0.4 + cut*0.1;
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event); 
            float value = (float)tree->GetBranch("jet3MinDR")->GetLeaf("jet3MinDR")->GetValue();
            if (value > 0 && value < dRCut[cut]) {
                count++; 
            }
        }
        std::cout << setw(12) << 100*count/(float)tree->GetEntries() << setw(12) << dRCut[cut] << std::endl;
        mergeRate[cut] = 100*count/(float)tree->GetEntries();
    }

    TGraph* graph = new TGraph(17,dRCut,mergeRate);
    graph->GetXaxis()->SetTitle("#DeltaR Cut");
    graph->GetYaxis()->SetTitle("Merger rate (%)");
    graph->SetLineColor(colour);
    graph->SetLineWidth(2);
    graph->SetTitle("Merger rate vs #DeltaR Cut");
    graph->SetName(Form("%d", colour));
    return graph;

};

void MergeCuts() {

    TString path = "Taggers/test/VBFProduction/output_vbfDumper_cfg_3rdjetmerge_dijet/";
    TString file[11];
    file[0]  = "output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0.root";
    file[1]  = "output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_0.root"; 
    file[2]  = "output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[3]  = "output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[4]  = "output_GluGluHToGG_M-125_13TeV_powheg_pythia8_0.root";
    file[5]  = "output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[6]  = "output_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[7]  = "output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[8]  = "output_VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_0.root";
    file[9]  =  "output_ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_0.root";
    file[10] = "output_VBFHToGG_M-125_13TeV_powheg_pythia8_0.root";

    TFile *inputFile[11];
    for (unsigned sample(0);sample<11;sample++){
        inputFile[sample] = TFile::Open(path + file[sample]);
        if (!inputFile[sample]->IsOpen()) {return;} else {std::cout << "\t" << path + file[sample] << " Opened successfully" << std::endl;}
    }
    std::cout << std::endl;

    TString treePath = "vbfTagDumper/trees/";
    TString treeName[11];
    treeName[0] = "dy_toll_m50_13TeV_VBFDiJet";
    treeName[1] = "gamgamjetbox_13TeV_VBFDiJet";
    treeName[2] = "gamJet_13TeV_VBFDiJet";
    treeName[3] = "gamJet_13TeV_VBFDiJet";
    treeName[4] = "ggf_m125_13TeV_VBFDiJet";
    treeName[5] = "qcd_13TeV_VBFDiJet";
    treeName[6] = "qcd_13TeV_VBFDiJet";
    treeName[7] = "qcd_13TeV_VBFDiJet";
    treeName[8] = "vh_m125_13TeV_VBFDiJet";
    treeName[9] = "tth_m125_13TeV_VBFDiJet";
    treeName[10] = "vbf_m125_13TeV_VBFDiJet";

    TTree* tree[11];
    for (unsigned sample(0);sample<11;sample++) {
        tree[sample] = (TTree*)inputFile[sample]->Get(treePath + treeName[sample]);
    }

    TGraph * mergeRate[11];
    for (unsigned sample(0);sample<11;sample++) {
        std::cout << "\n\tSample: " << file[sample] << std::endl;
        mergeRate[sample] = getMergeVsDRCut(tree[sample],sample+1);
        if (sample +1 == 10) {
            mergeRate[sample]->SetLineColor(1);
            mergeRate[sample]->SetLineStyle(2);
        }
    }

    TCanvas c1("c1");
    mergeRate[10]->Draw();
    for (unsigned sample(0);sample<10;sample++) {
        mergeRate[sample]->Draw("same");
    }

    leg = new TLegend(0.1,0.7,0.48,0.9);
    for (unsigned sample(0);sample<11;sample++) {
        leg->AddEntry(Form ("%d", sample+1),file[sample]);
    }
    leg->Draw("same");
    c1.Print("Plots/MergeRates.pdf");
    
}


