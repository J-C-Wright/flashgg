
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

    TTree *jjjTree = (TTree*)inputFile->Get("jjj");
    TTree *jjfTree = (TTree*)inputFile->Get("jjf");
    TTree *jffTree = (TTree*)inputFile->Get("jff");
    TTree *fffTree = (TTree*)inputFile->Get("fff");
    TTree *jjTree  = (TTree*)inputFile->Get("jj");
    TTree *jfTree  = (TTree*)inputFile->Get("jf");
    TTree *ffTree  = (TTree*)inputFile->Get("ff");

    jjjTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    jjjTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    jjjTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    jjjTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 
    jjfTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    jjfTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    jjfTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    jjfTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 
    jffTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    jffTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    jffTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    jffTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 
    fffTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    fffTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    fffTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    fffTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 
    jjTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    jjTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    jjTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    jjTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 
    jfTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    jfTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    jfTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    jfTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 
    ffTree->SetBranchAddress("recoLevel",&recoLevel.leadingJetPt); 
    ffTree->SetBranchAddress("genJetLevel",&genJetLevel.leadingJetPt); 
    ffTree->SetBranchAddress("genParticleLevel",&genParticleLevel.leadingJetPt); 
    ffTree->SetBranchAddress("partonLevel",&partonLevel.leadingJetPt); 

    



}

