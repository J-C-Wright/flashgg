
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

using namespace std;
using namespace edm;
using namespace flashgg;

class VBFPlotProducer{

    private:
        VBFTagTruth *truthPtr_;         
        TString label_;

        //Plots
        TH1F *leadJetPt_;
        TH1F *leadJetEta_;
        TH1F *leadJetPhi_;



    public:
        VBFPlotProducer(TString label) {
            label_ = label;
            //Testing
            std::cout << "Constructor test: the label is " << label_ << std::endl;

            //Set up plots


        } 
    
        void fill(VBFTagTruth *truthPtr) {
            std::cout << "Lead parton eta " << truthPtr->leadingParton()->eta() << std::endl;
        } 


};
