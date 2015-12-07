
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
        mergeRate[cut] = count/(float)tree->GetEntries();
    }

    TGraph* graph = new TGraph(17,dRCut,mergeRate);
    graph->SetLineColor(colour);
    graph->SetLineWidth(2);
    return graph;

};

void RR_Test() {

    TFile *inputFile  = TFile::Open("output.root");

    TTree *tree = (TTree*)inputFile->Get("vbfTagDumper/trees/vbfh_13TeV_VBFDiJet");
    std::cout << tree->GetEntries() << std::endl;

    TString sampleName = "Signal";
    TGraph * signal = getMergeVsDRCut(tree,1);
    signal->Draw();
}
