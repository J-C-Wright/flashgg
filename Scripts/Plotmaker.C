//#include "flashgg/DataFormats/interface/VBFTagTruth.h"

{

    TFile *inputFile  = TFile::Open("VBF_Output.root");
    TFile *outputFile = new TFile("MVA_Plots.root","RECREATE");

    if (!TClassTable::GetDict("VBFTagTruth")) {
    gSystem->Load("flashgg/DataFormats/interface/VBFTagTruth.h");
    }

    TTree *jjjTree = (TTree*)inputFile->Get("jjj");

    VBFTagTruth truth = new VBFTagTruth;
    jjjTree->SetBranchAddress("jjj",&truth);   


}

