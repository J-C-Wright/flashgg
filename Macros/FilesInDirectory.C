#include <vector>
#include <algorithm>

void FilesInDirectory() {

    gErrorIgnoreLevel = kInfo;

    const char *directoryName = "/vols/cms/jwright/RadRecTrees/dR_0_8";
    const char *extension = ".root";

    TFile *outputFile = new TFile("Output.root","RECREATE");
    outputFile->cd();

    TSystemDirectory dir(directoryName,directoryName);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fileName;
        TIter dirIter(files);

        std::vector<TString> channels;
        std::vector<TString> fileNames;

        while ((file=(TSystemFile*)dirIter())) {
            fileName = file->GetName();
            if (!file->IsDirectory() && fileName.EndsWith(extension)) {

                fileNames.push_back(fileName);

                TObjArray *array = fileName.Tokenize("_");
                TString temp = ((TObjString *)(array->At(1)))->String();

                unsigned count(0);
                for (unsigned channel(0);channel<channels.size();channel++) {
                    if (!temp.CompareTo(channels[channel])) { count++; }
                }
                if (count == 0) {channels.push_back(temp);}
            }
        }

        for (unsigned channel(0);channel<channels.size();channel++) {
            std::cout << channels[channel] << std::endl;
        }
        
        std::vector<TList*>  channelTreeLists(channels.size());
        for (unsigned channel(0);channel<channels.size();channel++) {
            channelTreeLists[channel] = new TList;
        }
        
        TString treePath = "vbfTagDumper/trees/";
        TString treeName[8];
        treeName[0] = "dy_toll_m50_13TeV_VBFDiJet";
        treeName[1] = "gamJet_13TeV_VBFDiJet";
        treeName[2] = "qcd_13TeV_VBFDiJet";
        treeName[3] = "gamgamjetbox_13TeV_VBFDiJet";
        treeName[4] = "vbf_m125_13TeV_VBFDiJet";
        treeName[5] = "ggf_m125_13TeV_VBFDiJet";
        treeName[6] = "tth_m125_13TeV_VBFDiJet";
        treeName[7] = "vh_m125_13TeV_VBFDiJet";

        for (unsigned file(0);file<10;file++) {
        //for (unsigned file(0);file<fileNames.size();file++) {

            TObjArray *array = fileNames[file].Tokenize("_");
            TString temp = ((TObjString *)(array->At(1)))->String();
            unsigned location(0);
            for(unsigned channel(0);channel<channels.size();channel++) {
                if (!temp.CompareTo(channels[channel])) location = channel;
            } 

            TFile *input = TFile::Open( TString(directoryName) + TString("/")+ fileNames[file] );
            TTree *treeTemp = (TTree*)input->Get(TString("vbfTagDumper/trees/") + treeName[location]);
            channelTreeLists[location]->Add(treeTemp);

        }


        std::vector<TTree*> mergedTrees(channels.size());
        for (unsigned channel(0);channel<channels.size();channel++) {
            mergedTrees[channel] = TTree::MergeTrees(channelTreeLists[channel]);
            mergedTrees[channel]->SetName(channels[channel]);
        }




    }
}













