
{
struct MVAVarStruct {

    float leadingJetPt, subLeadingJetPt, subSubLeadingJetPt;
    float leadingJetEta, subLeadingJetEta, subSubLeadingJetEta;
    float leadingJetPhi, subLeadingJetPhi, subSubLeadingJetPhi;

    int   leadingJetHemisphere, subLeadingJetHemisphere, subSubLeadingJetHemisphere;
   
    float dR_12, dR_13, dR_23;
    float mjj_12, mjj_13, mjj_23, mjjj;
    float dEta_12, dEta_13, dEta_23;
    float zepjj_12, zepjj_13, zepjj_23, zepjjj;
    float dPhijj_12, dPhijj_13, dPhijj_23, dPhijjj;

    float leadingDR, subLeadingDR, subSubLeadingDR;
};


    TFile *inputFile  = TFile::Open("MVA_Var_Trees_VBF.root");
    TFile *outputFile = new TFile("Plots/MVA_Plots.root","RECREATE");

    MVAVarStruct recoLevel;
    MVAVarStruct genJetLevel;
    MVAVarStruct genParticleLevel;
    MVAVarStruct partonLevel;

    unsigned numTrees(7);
    std::vector<TString> treeNames(numTrees);
    treeNames[0] = "jjj";treeNames[1] = "jjf";treeNames[2] = "jff";treeNames[3] = "fff";
    treeNames[4] = "jj";treeNames[5] = "jf";treeNames[6] = "ff";
    std::vector<TTree*> treeVector(numTrees);

    for (unsigned tree(0);tree<numTrees;tree++) {
        treeVector[tree] = (TTree*)inputFile->Get(treeNames[tree]);
        treeVector[tree]->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt);
        treeVector[tree]->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt);
        treeVector[tree]->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt);
        treeVector[tree]->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt);
    }        

    int numBranches = treeVector[0]->GetListOfBranches()->GetEntries();
    int numLeaves = treeVector[0]->GetBranch("recoLevel")->GetListOfLeaves()->GetEntries();
 
//Histograms
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
     
    //Initialize histograms
    TObjArray *branchNames = treeVector[0]->GetListOfBranches();
    TObjArray *leafNames = treeVector[0]->GetBranch("recoLevel")->GetListOfLeaves(); 
    for (unsigned tree(0);tree<numTrees;tree++) {
        for (unsigned branch(0);branch<numBranches;branch++) {
            for (unsigned leaf(0);leaf<numLeaves;leaf++) {
                //Find range
                float minVal(999.), maxVal(0.);
                for (unsigned event(0);event<treeVector[tree]->GetEntries();event++) {    
                    treeVector[tree]->GetEntry(event);
                    float value = (float)treeVector[tree]->GetBranch(branchNames->At(branch)->GetName())->GetLeaf(leafNames->At(leaf)->GetName())->GetValue();
                    if (value > maxVal) maxVal = value;
                    if (value < minVal) minVal = value;
                }
                hists[tree][branch][leaf] = new TH1F(TString(branchNames->At(branch)->GetName()) + TString(leafNames->At(leaf)->GetName()) + TString("_") + treeNames[tree],
                                                   "",
                                                   50,minVal,maxVal);
            }
        }
    }

    //Histograms fill
    for (unsigned tree(0);tree<numTrees;tree++) {
        for (unsigned event(0);event<treeVector[tree]->GetEntries();event++) {
            treeVector[tree]->GetEntry(event);
            for (unsigned branch(0);branch<numBranches;branch++) {
                for (unsigned leaf(0);leaf<numLeaves;leaf++) {
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

//Combine histograms
//Three Jet Plots
    TCanvas c1( "c1" );
    gStyle->SetOptStat( 0 );
    for (unsigned branch(0);branch<numBranches;branch++) {
        for (unsigned leaf(0);leaf<numLeaves;leaf++) {
            //Find which category has largest value
            float maxVal(0); unsigned maxTree(0);
            for (unsigned tree(0);tree<5;tree++) {
                if (hists[tree][branch][leaf]->GetMaximum() > maxVal) { 
                    maxVal  = hists[tree][branch][leaf]->GetMaximum();
                    maxTree = tree;
                }
                if (tree == 4) {hists[tree][branch][leaf]->SetLineColor(tree+2);}else{hists[tree][branch][leaf]->SetLineColor(tree);}
            } 
            hists[maxTree][branch][leaf]->Draw();
            for (unsigned tree(0);tree<5;tree++) {
                if (tree != maxTree) {
                    hists[tree][branch][leaf]->Draw("same");
                }
            }            
            c1.Print("Plots/" + TString(branchNames->At(branch)->GetName()) + TString(leafNames->At(leaf)->GetName()) + "_3J.pdf");
        }
    }
        

//Two jet plots
    for (unsigned branch(0);branch<numBranches;branch++) {
        for (unsigned leaf(0);leaf<numLeaves;leaf++) {
            //Find which category has largest value
            float maxVal(0); unsigned maxTree(0);
            for (unsigned tree(5);tree<numTrees;tree++) {
                if (hists[tree][branch][leaf]->GetMaximum() > maxVal) { 
                    maxVal  = hists[tree][branch][leaf]->GetMaximum();
                    maxTree = tree;
                }
                hists[tree][branch][leaf]->SetLineColor(tree-4);
            } 
            hists[maxTree][branch][leaf]->Draw();
            for (unsigned tree(5);tree<numTrees;tree++) {
                if (tree != maxTree) {
                    hists[tree][branch][leaf]->Draw("same");
                }
            }            
            c1.SetTitle(TString(leafNames->At(leaf)->GetName()) + TString(" at ") + TString(branchNames->At(branch)->GetName()));
            c1.Print("Plots/" + TString(branchNames->At(branch)->GetName()) + TString(leafNames->At(leaf)->GetName()) + "_2J.pdf");
        }
    }
}   




