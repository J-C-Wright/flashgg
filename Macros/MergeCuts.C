
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

TH1F* getShiftDistribution(TTree *tree, TString mergedVar, TString unmergedVar, unsigned colour, float min, float max) {

    TH1F *hist = new TH1F(mergedVar + Form("%d", colour),mergedVar,50,min,max);

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        float merged   = (float)tree->GetBranch(mergedVar)->GetLeaf(mergedVar)->GetValue();
        float unmerged = (float)tree->GetBranch(unmergedVar)->GetLeaf(unmergedVar)->GetValue();
        float weight   = (float)tree->GetBranch("weight")->GetLeaf("weight")->GetValue();
        if (merged < -900 || unmerged < -900) continue;
        hist->Fill(merged-unmerged,weight);
    }

    hist->SetLineColor(colour);
    hist->GetXaxis()->SetTitle(mergedVar);
    hist->Scale(1/hist->Integral());
    return hist;
}; 

std::pair<float,float> getMaxMinShiftDistribution(TTree *tree, TString mergedVar, TString unmergedVar) {

    std::pair<float,float>  limits;
    limits.first  = 999.;
    limits.second = 0.0;

    std::cout << "\nFinding limits for " << mergedVar << " in tree " << tree->GetName() << std::endl;
    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        float merged   = (float)tree->GetBranch(mergedVar)->GetLeaf(mergedVar)->GetValue();
        float unmerged = (float)tree->GetBranch(unmergedVar)->GetLeaf(unmergedVar)->GetValue();
        if (merged < -900 || unmerged < -900) continue;
        if (merged-unmerged < limits.first) {
            limits.first = merged-unmerged;
            std::cout << "New min " << setw(12) << merged << setw(12) << unmerged << setw(12) << merged-unmerged << std::endl;
        }
        if (merged-unmerged > limits.second){ 
            limits.second = merged-unmerged;
            std::cout << "New max " << setw(12) << merged << setw(12) << unmerged << setw(12) << merged-unmerged << std::endl;
        }
    }

    return limits;
}






void MergeCuts() {

    TString path = "/vols/cms/jwright/RadRecTrees/dR_0_5/";
    TString file[9];
    file[0]  = "output_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0.root";
    file[1]  = "output_DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_0.root"; 
    file[2]  = "output_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[3]  = "output_GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[4]  = "output_GluGluHToGG_M-125_13TeV_powheg_pythia8_0.root";
    file[5]  = "output_QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[6]  = "output_QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[7]  = "output_QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8_0.root";
    file[8] = "output_VBFHToGG_M-125_13TeV_powheg_pythia8_0.root";

    TCanvas c1("c1");

    TFile *inputFile[9];
    for (unsigned sample(0);sample<9;sample++){
        inputFile[sample] = TFile::Open(path + file[sample]);
        if (!inputFile[sample]->IsOpen()) {return;} else {std::cout << "\t" << path + file[sample] << " Opened successfully" << std::endl;}
    }

    TString treePath = "vbfTagDumper/trees/";
    TString treeName[9];
    treeName[0] = "dy_toll_m50_13TeV_VBFDiJet";
    treeName[1] = "gamgamjetbox_13TeV_VBFDiJet";
    treeName[2] = "gamJet_13TeV_VBFDiJet";
    treeName[3] = "gamJet_13TeV_VBFDiJet";
    treeName[4] = "ggf_m125_13TeV_VBFDiJet";
    treeName[5] = "qcd_13TeV_VBFDiJet";
    treeName[6] = "qcd_13TeV_VBFDiJet";
    treeName[7] = "qcd_13TeV_VBFDiJet";
    treeName[8] = "vbf_m125_13TeV_VBFDiJet";

    TTree* tree[9];
    for (unsigned sample(0);sample<9;sample++) {
        tree[sample] = (TTree*)inputFile[sample]->Get(treePath + treeName[sample]);
        std::cout << "Tree " << treeName[sample] << " loaded" << std::endl;
    }
    
//Tree merging
    //QCD
    TList *qcdList = new TList;
    qcdList->Add(tree[5]);qcdList->Add(tree[6]);qcdList->Add(tree[7]);
    TTree *qcdTree = TTree::MergeTrees(qcdList);
    qcdTree->SetName("qcdTree");
    qcdTree->SetDirectory(0);
    //GJ
    TList *gjList = new TList;
    gjList->Add(tree[2]);gjList->Add(tree[3]);
    TTree *gjTree = TTree::MergeTrees(gjList);
    gjTree->SetName("gjTree");
    gjTree->SetDirectory(0);
    //Move elements of the array about
    tree[2] = gjTree;
    treeName[2] = gjTree->GetName();
    tree[3] = qcdTree;
    treeName[3] = qcdTree->GetName();
    tree[5] = tree[8];
    treeName[5] = treeName[8];

    TString labels[6];
    labels[0] = "Drell-Yan";
    labels[1] = "GGJetbox";
    labels[2] = "GJ";
    labels[3] = "QCD";
    labels[4] = "GGF";
    labels[5] = "VBF";

//Merge rates
    TGraph * mergeRate[6];
    for (unsigned sample(0);sample<6;sample++) {
        std::cout << "\n\tSample: " << treeName[sample] << std::endl;
        mergeRate[sample] = getMergeVsDRCut(tree[sample],sample+1);
    }

    mergeRate[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        mergeRate[sample]->Draw("same");
    }

    leg = new TLegend(0.1,0.7,0.48,0.9);
    for (unsigned sample(0);sample<6;sample++) {
        leg->AddEntry(Form ("%d", sample+1),labels[sample]);
    }
    leg->Draw("same");
    c1.Print("Plots/MergeRates.pdf");


//Test shift distribution
    //Invariant mass
    std::pair<float,float> mjjLimits;
    mjjLimits.first  = 999.;
    mjjLimits.second = 0.0;
    for (unsigned sample(0);sample<6;sample++) {
        std::pair<float,float> tempLimits = getMaxMinShiftDistribution(tree[sample], "dijet_Mjj", "J1J2_mjj"); 
        if (tempLimits.first < mjjLimits.first) mjjLimits.first = tempLimits.first;
        if (tempLimits.second > mjjLimits.second) mjjLimits.second = tempLimits.second;
    }
    std::cout << "Final mjjLimits" << std::endl;
    std::cout << setw(12) << mjjLimits.first << setw(12) << mjjLimits.second << std::endl;
        
    TH1F* mjjShift[6];
    for (unsigned sample(0);sample<6;sample++) {
        mjjShift[sample] = getShiftDistribution(tree[sample], "dijet_Mjj", "J1J2_mjj", sample+1, mjjLimits.first, mjjLimits.second);
    }
    mjjShift[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        mjjShift[sample]->Draw("HISTsame");
        if (sample +1 == 10) {
            mjjShift[sample]->SetLineColor(1);
            mjjShift[sample]->SetLineStyle(2);
        }
    }
    gPad->SetLogy();
    c1.Print("Plots/mjjShift.pdf");

    //Lead pT
    std::pair<float,float> leadPtLimits;
    leadPtLimits.first  = 999.;
    leadPtLimits.second = 0.0;
    for (unsigned sample(0);sample<6;sample++) {
        std::pair<float,float> tempLimits = getMaxMinShiftDistribution(tree[sample], "dijet_LeadJPt", "jet1_pt");
        if (tempLimits.first < leadPtLimits.first) leadPtLimits.first = tempLimits.first;
        if (tempLimits.second > leadPtLimits.second) leadPtLimits.second = tempLimits.second;
    }
    std::cout << "Final leadPtLimits" << std::endl;
    std::cout << setw(12) << leadPtLimits.first << setw(12) << leadPtLimits.second << std::endl;

    TH1F* leadPtShift[6];
    for (unsigned sample(0);sample<6;sample++) {
        leadPtShift[sample] = getShiftDistribution(tree[sample], "dijet_LeadJPt", "jet1_pt", sample+1, leadPtLimits.first, leadPtLimits.second);
    }
    leadPtShift[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        leadPtShift[sample]->Draw("HISTsame");
        if (sample +1 == 10) {
            leadPtShift[sample]->SetLineColor(1);
            leadPtShift[sample]->SetLineStyle(2);
        }
    }
    c1.Print("Plots/leadPtShift.pdf");

    //Sublead pT
    std::pair<float,float> subleadPtLimits;
    subleadPtLimits.first  = 999.;
    subleadPtLimits.second = 0.0;
    for (unsigned sample(0);sample<6;sample++) {
        std::pair<float,float> tempLimits = getMaxMinShiftDistribution(tree[sample], "dijet_SubJPt", "jet2_pt");
        if (tempLimits.first < subleadPtLimits.first) subleadPtLimits.first = tempLimits.first;
        if (tempLimits.second > subleadPtLimits.second) subleadPtLimits.second = tempLimits.second;
    }
    std::cout << "Final subleadPtLimits" << std::endl;
    std::cout << setw(12) << subleadPtLimits.first << setw(12) << subleadPtLimits.second << std::endl;

    TH1F* subleadPtShift[6];
    for (unsigned sample(0);sample<6;sample++) {
        subleadPtShift[sample] = getShiftDistribution(tree[sample], "dijet_SubJPt", "jet2_pt", sample+1, subleadPtLimits.first, subleadPtLimits.second);
    }
    subleadPtShift[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        subleadPtShift[sample]->Draw("HISTsame");
        if (sample +1 == 10) {
            subleadPtShift[sample]->SetLineColor(1);
            subleadPtShift[sample]->SetLineStyle(2);
        }
    }
    c1.Print("Plots/subleadPtShift.pdf");


    //dEta 
    std::pair<float,float> dEtaLimits;
    dEtaLimits.first  = 999.;
    dEtaLimits.second = 0.0;
    for (unsigned sample(0);sample<6;sample++) {
        std::pair<float,float> tempLimits = getMaxMinShiftDistribution(tree[sample], "dijet_abs_dEta", "J1J2_dEta");
        if (tempLimits.first < dEtaLimits.first) dEtaLimits.first = tempLimits.first;
        if (tempLimits.second > dEtaLimits.second) dEtaLimits.second = tempLimits.second;
    }
    std::cout << "Final dEtaLimits" << std::endl;
    std::cout << setw(12) << dEtaLimits.first << setw(12) << dEtaLimits.second << std::endl;

    TH1F* dEtaShift[6];
    for (unsigned sample(0);sample<6;sample++) {
        dEtaShift[sample] = getShiftDistribution(tree[sample], "dijet_abs_dEta", "J1J2_dEta", sample+1, dEtaLimits.first, dEtaLimits.second);
    }
    dEtaShift[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        dEtaShift[sample]->Draw("HISTsame");
        if (sample +1 == 10) {
            dEtaShift[sample]->SetLineColor(1);
            dEtaShift[sample]->SetLineStyle(2);
        }
    }
    c1.Print("Plots/dEtaShift.pdf");


    //dPhi
    std::pair<float,float> dPhiLimits;
    dPhiLimits.first  = 999.;
    dPhiLimits.second = 0.0;
    for (unsigned sample(0);sample<6;sample++) {
        std::pair<float,float> tempLimits = getMaxMinShiftDistribution(tree[sample], "dijet_dphi", "J1J2_dipho_dPhi");
        if (tempLimits.first < dPhiLimits.first) dPhiLimits.first = tempLimits.first;
        if (tempLimits.second > dPhiLimits.second) dPhiLimits.second = tempLimits.second;
    }
    std::cout << "Final dPhiLimits" << std::endl;
    std::cout << setw(12) << dPhiLimits.first << setw(12) << dPhiLimits.second << std::endl;

    TH1F* dPhiShift[6];
    for (unsigned sample(0);sample<6;sample++) {
        dPhiShift[sample] = getShiftDistribution(tree[sample], "dijet_dphi", "J1J2_dipho_dPhi", sample+1, dPhiLimits.first, dPhiLimits.second);
    }
    dPhiShift[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        dPhiShift[sample]->Draw("HISTsame");
        if (sample +1 == 10) {
            dPhiShift[sample]->SetLineColor(1);
            dPhiShift[sample]->SetLineStyle(2);
        }
    }
    c1.Print("Plots/dPhiShift.pdf");


    //Zep
    std::pair<float,float> ZepLimits;
    ZepLimits.first  = 999.;
    ZepLimits.second = 0.0;
    for (unsigned sample(0);sample<6;sample++) {
        std::pair<float,float> tempLimits = getMaxMinShiftDistribution(tree[sample], "dijet_Zep", "J1J2_Zep");
        if (tempLimits.first < ZepLimits.first) ZepLimits.first = tempLimits.first;
        if (tempLimits.second > ZepLimits.second) ZepLimits.second = tempLimits.second;
    }
    std::cout << "Final ZepLimits" << std::endl;
    std::cout << setw(12) << ZepLimits.first << setw(12) << ZepLimits.second << std::endl;

    TH1F* ZepShift[6];
    for (unsigned sample(0);sample<6;sample++) {
        ZepShift[sample] = getShiftDistribution(tree[sample], "dijet_Zep", "J1J2_Zep", sample+1, ZepLimits.first, ZepLimits.second);
    }
    ZepShift[5]->Draw();
    for (unsigned sample(0);sample<6;sample++) {
        ZepShift[sample]->Draw("HISTsame");
        if (sample +1 == 10) {
            ZepShift[sample]->SetLineColor(1);
            ZepShift[sample]->SetLineStyle(2);
        }
    }
    c1.Print("Plots/ZepShift.pdf");

}


