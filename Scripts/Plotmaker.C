
{
 struct MVAVarStruct {

    float leadingJetPt, subLeadingJetPt, subSubLeadingJetPt;
    float leadingJetEta, subLeadingJetEta, subSubLeadingJetEta;
    float leadingJetPhi, subLeadingJetPhi, subSubLeadingJetPhi;

    int   leadingJetHemisphere, subLeadingJetHemisphere, subSubLeadingJetHemisphere;   
    int   oppHemispheres_J1J2, oppHemispheres_J1J3, oppHemispheres_J2J3;
 
    float dR_12, dR_13, dR_23;
    float mjj_12, mjj_13, mjj_23, mjjj;
    float dEta_12, dEta_13, dEta_23;
    float zepjj_12, zepjj_13, zepjj_23, zepjjj;
    float dPhijj_12, dPhijj_13, dPhijj_23, dPhijjj;

    float dEta_J1J2J3, dEta_J2J3J1, dEta_J3J1J2;

    float mjj_d12_13plus23, mjj_d12_13, mjj_d12_23, mjj_d13_23;
    float dR_DP_12, dR_DP_13, dR_DP_23;
    float dR_Ph1_1,dR_Ph1_2,dR_Ph1_3,dR_Ph2_1,dR_Ph2_2,dR_Ph2_3;
    float dR_DP_123; 

    float leadingDR, subLeadingDR, subSubLeadingDR;
};
struct manualLimit {
    TString leafName;
    float   minValue;
    float   maxValue;
};



    std::cout << "Starting macro" << std::endl;

    TFile *inputFileVBF  = TFile::Open("MVA_Var_Trees_VBF.root");
    TFile *inputFileggH  = TFile::Open("MVA_Var_Trees_ggH.root");
    TFile *outputFile    = new TFile("Plots/MVA_Plots.root","RECREATE");

    MVAVarStruct recoLevel;
    MVAVarStruct genJetLevel;
    MVAVarStruct genParticleLevel;
    MVAVarStruct partonLevel;

    unsigned numTrees(9);
    std::vector<TString> treeNames(numTrees);
    treeNames[0] = "jjj";treeNames[1] = "jjf";treeNames[2] = "jff";treeNames[3] = "fff";
    treeNames[4] = "jj";treeNames[5] = "jf";treeNames[6] = "ff";
    treeNames[7] = "ggHjjj";treeNames[8] = "ggHjj";
    std::vector<TTree*> treeVector(numTrees);

    //VBF Trees
    for (unsigned tree(0);tree<numTrees-2;tree++) {
        treeVector[tree] = (TTree*)inputFileVBF->Get(treeNames[tree]);
        treeVector[tree]->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt);
        treeVector[tree]->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt);
        treeVector[tree]->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt);
        treeVector[tree]->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt);
    }        

    //ggH Trees
    TList *ggHjjj = new TList;
    for (unsigned tree(0);tree<4;tree++) {
        ggHjjj->Add((TTree*)inputFileggH->Get(treeNames[tree]));
    }        
    ggHjjjTree = TTree::MergeTrees(ggHjjj);
    ggHjjjTree->SetName("ggHjjj");
    ggHjjjTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt);
    ggHjjjTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt);
    ggHjjjTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt);
    ggHjjjTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt);

    TList *ggHjj = new TList;
    for (unsigned tree(4);tree<7;tree++) {
        ggHjj->Add((TTree*)inputFileggH->Get(treeNames[tree]));
    }        
    ggHjjTree = TTree::MergeTrees(ggHjj);
    ggHjjTree->SetName("ggHjj");
    ggHjjTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt);
    ggHjjTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt);
    ggHjjTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt);
    ggHjjTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt);

    //Would like the 3J ggH to be with the VBF jjj so it's easier to loop over them
    //Move tree 7 to position 5 and put jj ggH on the end, remove the spare entry
    treeVector.insert(treeVector.begin()+4,ggHjjjTree);
    treeVector[numTrees-1] = ggHjjTree;
    treeVector.pop_back();

    TObjArray *branchNames = treeVector[0]->GetListOfBranches();
    TObjArray *leafNames = treeVector[0]->GetBranch("recoLevel")->GetListOfLeaves(); 
    int numBranches = branchNames->GetEntries();
    int numLeaves   = leafNames->GetEntries();

    //Set manual ranges for certain variables - move this to an external file...
    std::vector<manualLimit> userLimits(11);
    userLimits[0].leafName = "leadingJetPt"; userLimits[0].minValue = 0; userLimits[0].maxValue = 400; 
    userLimits[1].leafName = "subLeadingJetPt"; userLimits[1].minValue = 0; userLimits[1].maxValue = 120; 
    userLimits[2].leafName = "subSubLeadingJetPt"; userLimits[2].minValue = 0; userLimits[2].maxValue = 70; 
    userLimits[3].leafName = "mjj_12"; userLimits[3].minValue = 0; userLimits[3].maxValue = 2000; 
    userLimits[4].leafName = "mjj_13"; userLimits[4].minValue = 0; userLimits[4].maxValue = 1000; 
    userLimits[5].leafName = "mjj_23"; userLimits[5].minValue = 0; userLimits[5].maxValue = 1000; 
    userLimits[6].leafName = "mjj_d12_13"; userLimits[6].minValue = -1000; userLimits[6].maxValue = 2000; 
    userLimits[7].leafName = "mjj_d12_23"; userLimits[7].minValue = -1000; userLimits[7].maxValue = 2000; 
    userLimits[8].leafName = "mjj_d13_23"; userLimits[8].minValue = -1000; userLimits[8].maxValue = 1000; 
    userLimits[9].leafName = "mjj_d12_13_plus23"; userLimits[9].minValue = -1000; userLimits[9].maxValue = 2000; 
    userLimits[10].leafName = "mjjj"; userLimits[10].minValue = 0; userLimits[10].maxValue = 3000; 



//Histogram construction and filling
    //Make vectors of hists. One for each tree, branch, and leaf
    std::vector<std::vector<std::vector<TH1F*>>> hists(numTrees);
    for (unsigned tree(0);tree<numTrees;tree++) {
        std::vector<std::vector<TH1F*>> branchHists(numBranches);
        for (unsigned branch(0);branch<numBranches;branch++) {
            std::vector<TH1F*> leafHists(numLeaves);
            branchHists[branch] = leafHists;
        }
        hists[tree] = branchHists;
    }

    //Build histograms
    std::vector<std::pair<float,float>> rangeVector(numLeaves); 
    unsigned numBins(100);
    for (unsigned branch(0);branch<numBranches;branch++) {
        for (unsigned leaf(0);leaf<numLeaves;leaf++) {

            rangeVector[leaf].first  = 999.;
            rangeVector[leaf].second = 0.0;
            //Is it user-defined?
            bool userDefinedLimits = false; unsigned userLimitIndex(0);
            for (unsigned var(0);var<userLimits.size();var++) {
                if (leafNames->At(leaf)->GetName() == userLimits[var].leafName) {
                    userDefinedLimits = true;
                    userLimitIndex = var;
                }
            }

            if (userDefinedLimits) {
                //Set predefined range
                std::cout << "User-defined range..." << std::endl;
                rangeVector[leaf].first  = userLimits[userLimitIndex].minValue;           
                rangeVector[leaf].second = userLimits[userLimitIndex].maxValue;
            }else{           
                //Find range
                for (unsigned tree(0);tree<numTrees;tree++) {
                    for (unsigned event(0);event<treeVector[tree]->GetEntries();event++) {    
                        treeVector[tree]->GetEntry(event);
                        float value = (float)treeVector[tree]->GetBranch(branchNames->At(branch)->GetName())->GetLeaf(leafNames->At(leaf)->GetName())->GetValue();
                        if (value < rangeVector[leaf].first  && value > -990) rangeVector[leaf].first = value;
                        if (value > rangeVector[leaf].second && value > -990) rangeVector[leaf].second = value;
                    }
                }

            }
            std::cout << "For branch " << branchNames->At(branch)->GetName() << " and leaf " << leafNames->At(leaf)->GetName();
            std::cout << " min is " << rangeVector[leaf].first << " and max is " << rangeVector[leaf].second << std::endl;  
            
            //Initialize histogram
            for (unsigned tree(0);tree<numTrees;tree++) {
                hists[tree][branch][leaf] = new TH1F(TString(branchNames->At(branch)->GetName()) + TString(leafNames->At(leaf)->GetName()) + TString("_") + treeNames[tree],
                                                     TString(branchNames->At(branch)->GetName()) + TString(" ") + TString(leafNames->At(leaf)->GetName()),
                                                     numBins,rangeVector[leaf].first,rangeVector[leaf].second);
            }
        }
    }
    //Histograms fill
    for (unsigned tree(0);tree<numTrees;tree++) {
        for (unsigned event(0);event<treeVector[tree]->GetEntries();event++) {
            treeVector[tree]->GetEntry(event);
            for (unsigned branch(0);branch<numBranches;branch++) {
                for (unsigned leaf(0);leaf<numLeaves;leaf++) {
                    if (-990 > (float)treeVector[tree]->GetBranch(branchNames->At(branch)->GetName())->GetLeaf(leafNames->At(leaf)->GetName())->GetValue()) continue;
                    hists[tree][branch][leaf]->Fill((float)treeVector[tree]->GetBranch(branchNames->At(branch)->GetName())->GetLeaf(leafNames->At(leaf)->GetName())->GetValue());
                } 
            }
        }
    }
    //Normalise histograms to area=1
    for (unsigned tree(0);tree<numTrees;tree++) {
        for (unsigned branch(0);branch<numBranches;branch++) {
            for (unsigned leaf(0);leaf<numLeaves;leaf++) {
                hists[tree][branch][leaf]->Scale(1/hists[tree][branch][leaf]->Integral());
            }
        }
    }

//Variable Comparison Plots
    TCanvas c1( "c1" );
    gStyle->SetOptStat( 0 );
    //Three Jet Plots
    unsigned firstDijetTree(5);
    for (unsigned branch(0);branch<numBranches;branch++) {
        for (unsigned leaf(0);leaf<numLeaves;leaf++) {
            //Find which tree has largest value
            float peak(0); unsigned maxTree(0);
            for (unsigned tree(0);tree<firstDijetTree;tree++) {
                if (hists[tree][branch][leaf]->GetMaximum() > peak) { 
                    peak  = hists[tree][branch][leaf]->GetMaximum();
                    maxTree = tree;
                }
                if (tree == 4) {hists[tree][branch][leaf]->SetLineColor(6);}else{hists[tree][branch][leaf]->SetLineColor(tree+1);}
            } 
            hists[maxTree][branch][leaf]->Draw();
            for (unsigned tree(0);tree<firstDijetTree;tree++) {
                if (tree != maxTree) {
                    hists[tree][branch][leaf]->Draw("same");
                }
            }            
            c1.Print("Plots/3J/" + TString(branchNames->At(branch)->GetName()) + "/"
                                 + TString(branchNames->At(branch)->GetName()) + "_" + TString(leafNames->At(leaf)->GetName()) + "_3J.pdf");
        }
    }
    //Two jet plots
    for (unsigned branch(0);branch<numBranches;branch++) {
        for (unsigned leaf(0);leaf<numLeaves;leaf++) {
            //Find which tree has largest value
            float peak(0); unsigned maxTree(0);
            for (unsigned tree(firstDijetTree);tree<numTrees;tree++) {
                if (hists[tree][branch][leaf]->GetMaximum() > peak) { 
                    peak  = hists[tree][branch][leaf]->GetMaximum();
                    maxTree = tree;
                }
                if (tree == 8) {hists[tree][branch][leaf]->SetLineColor(6);}else{hists[tree][branch][leaf]->SetLineColor(tree-3);}
            } 
            hists[maxTree][branch][leaf]->Draw();
            for (unsigned tree(firstDijetTree);tree<numTrees;tree++) {
                if (tree != maxTree) {
                    hists[tree][branch][leaf]->Draw("same");
                }
            }            
            c1.Print("Plots/2J/" + TString(branchNames->At(branch)->GetName()) + "/"
                                 + TString(branchNames->At(branch)->GetName()) + "_" + TString(leafNames->At(leaf)->GetName()) + "_2J.pdf");
        }
    }



//ROC Curves
    //Combination of hists into sigHist and bgrHist

    unsigned branch = 0;
    std::vector<TH1F*> sigHist(numLeaves);
    for(unsigned leaf(0);leaf<numLeaves;leaf++) {
        sigHist[leaf] = new TH1F("name",//TString(branchNames->At(branch)->GetName()) + TString(leafNames->At(leaf)->GetName()) + TString("_") + treeNames[tree] + TString("_S"),
                                 "",//TString(branchNames->At(branch)->GetName()) + TString(" ") + TString(leafNames->At(leaf)->GetName()) + TString("_S"),
                                 numBins,rangeVector[leaf].first,rangeVector[leaf].second);
        sigHist[leaf]->Add(hists[0][branch][leaf]);
        sigHist[leaf]->Add(hists[1][branch][leaf]);
        sigHist[leaf]->Scale(1/sigHist[leaf]->Integral());
    }
    std::vector<TH1F*> bgrHist(numLeaves);
    for(unsigned leaf(0);leaf<numLeaves;leaf++) {
        bgrHist[leaf] = new TH1F("name2",//TString(branchNames->At(branch)->GetName()) + TString(leafNames->At(leaf)->GetName()) + TString("_") + treeNames[tree] + TString("_BG"),
                                 "",//TString(branchNames->At(branch)->GetName()) + TString(" ") + TString(leafNames->At(leaf)->GetName()) + TString("_BG"),
                                 numBins,rangeVector[leaf].first,rangeVector[leaf].second);
        bgrHist[leaf]->Add(hists[2][branch][leaf]);
        bgrHist[leaf]->Add(hists[3][branch][leaf]);
        bgrHist[leaf]->Add(hists[4][branch][leaf]);
        bgrHist[leaf]->Scale(1/bgrHist[leaf]->Integral());
    }

    //ROC Curve construction
    Float_t FPR[101];
    Float_t TPR[101];

    for (unsigned leaf(0);leaf<numLeaves;leaf++) {
        for (unsigned i(0);i<numBins+1;i++) {
            FPR[i]  = bgrHist[leaf]->Integral(0,i);
            TPR[i]  = sigHist[leaf]->Integral(i,100);
        }
        TGraph* ROCCurve = new TGraph(101,FPR,TPR);
        ROCCurve->SetTitle(TString(leafNames->At(leaf)->GetName()) + TString("ROC Curve"));
        ROCCurve->SetLineWidth(3);
        ROCCurve->Draw();

        c1.Print("Plots/3J/" + TString(branchNames->At(branch)->GetName()) + "/ROCs/"
                             + TString(branchNames->At(branch)->GetName()) + "_" + TString(leafNames->At(leaf)->GetName()) + "_3J_ROC.pdf");
        c1.Clear();
    
    }

    //Scoring of the ROC Curves
    std::vector<std::pair<TString,float>> rankings;
    for (unsigned leaf(0);leaf<numLeaves;leaf++) {

        //Using trapezoid method to calculate area under ROC curve 
        float area(0.0);
        for (unsigned i(0);i<numBins-1;i++) {
            float FPR1,FPR2,TPR1,TPR2;
            FPR1 = bgrHist[leaf]->Integral(0,i);
            FPR2 = bgrHist[leaf]->Integral(0,i+1);
            TPR1 = sigHist[leaf]->Integral(i,100);
            TPR2 = sigHist[leaf]->Integral(i+1,100);
            area += 0.5*(FPR2-FPR1)*(TPR1+TPR2);
        }
        std::pair<TString,float> varPerf(TString(leafNames->At(leaf)->GetName()),area);
    
        //Insert in order
        unsigned int insertionIndex(0);
        for (unsigned int rank(0);rank<rankings.size();rank++) {
            if (varPerf.second < rankings[rank].second) { insertionIndex = rank + 1; }
        }
        rankings.insert( rankings.begin() + insertionIndex, varPerf);
    }

    std::cout << "The ranking of the variables by area under ROC curve is: " << std::endl;
    std::cout << setw(40) << "Variable Name" << setw(12) << "Area" << std::endl;
    for (unsigned int rank(0);rank<rankings.size();rank++) {
        if (rankings[rank].second < 1) {std::cout << setw(40) << rankings[rank].first << setw(12) << rankings[rank].second << std::endl;}
    }
        
}   


