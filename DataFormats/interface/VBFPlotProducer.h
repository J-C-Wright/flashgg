
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/VBFTagTruth.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace flashgg;

class VBFPlotProducer{

    private:
        VBFTagTruth truth_;         
        TString label_;
       
        TH1F *leadingJetPt_;
        TH1F *subLeadingJetPt_;
        TH1F *subSubLeadingJetPt_;
        TH1F *leadingJetEta_;
        TH1F *subLeadingJetEta_;
        TH1F *subSubLeadingJetEta_;
        TH1F *leadingJetPhi_;
        TH1F *subLeadingJetPhi_;
        TH1F *subSubLeadingJetPhi_;
       
    public:
        VBFPlotProducer(){}

        void setup(TString label) {

            label_ = label;

            leadingJetPt_ = new TH1F("leadingJetPt","leadingJetPt",50,0,500);
            subLeadingJetPt_ = new TH1F("subLeadingJetPt","subLeadingJetPt",50,0,500);
            subSubLeadingJetPt_ = new TH1F("subSubLeadingJetPt","subSubLeadingJetPt",50,0,500);
            leadingJetEta_ = new TH1F("leadingJetEta","leadingJetEta",50,-5,5);
            subLeadingJetEta_ = new TH1F("subLeadingJetEta","subLeadingJetEta",50,-5,5);
            subSubLeadingJetEta_ = new TH1F("subSubLeadingJetEta","subSubLeadingJetEta",50,-5,5);
            leadingJetPhi_ = new TH1F("leadingJetPhi","leadingJetPhi",50,-3.15,3.15);
            subLeadingJetPhi_ = new TH1F("subLeadingJetPhi","subLeadingJetPhi",50,-3.15,3.15);
            subSubLeadingJetPhi_ = new TH1F("subSubLeadingJetPhi","subSubLeadingJetPhi",50,-3.15,3.15);

        } 
    
        void fill(VBFTagTruth *truthPtr) {

            leadingJetPt_->Fill(truthPtr->leadingJet()->pt());
            subLeadingJetPt_->Fill(truthPtr->subLeadingJet()->pt());
            if (truthPtr->hasTrijet()) {subSubLeadingJetPt_->Fill(truthPtr->subSubLeadingJet()->pt());}
            leadingJetEta_->Fill(truthPtr->leadingJet()->eta());
            subLeadingJetEta_->Fill(truthPtr->subLeadingJet()->eta());
            if (truthPtr->hasTrijet()) {subSubLeadingJetEta_->Fill(truthPtr->subSubLeadingJet()->eta());}
            leadingJetPhi_->Fill(truthPtr->leadingJet()->phi());
            subLeadingJetPhi_->Fill(truthPtr->subLeadingJet()->phi());
            if (truthPtr->hasTrijet()) {subSubLeadingJetPhi_->Fill(truthPtr->subSubLeadingJet()->phi());}

        } 

        void write(TFile *file) {

            file->cd();
            TDirectory *cdlabel = file->mkdir( label_ );
            cdlabel->cd();

            leadingJetPt_->Write();
            subLeadingJetPt_->Write();
            subSubLeadingJetPt_->Write();
            leadingJetEta_->Write();
            subLeadingJetEta_->Write();
            subSubLeadingJetEta_->Write();
            leadingJetPhi_->Write();
            subLeadingJetPhi_->Write();
            subSubLeadingJetPhi_->Write();
        
        }


};
