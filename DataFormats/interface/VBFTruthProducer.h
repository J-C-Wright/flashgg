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
#include "flashgg/DataFormats/interface/VBFTagTruth.h"

using namespace std;
using namespace edm;
using namespace flashgg;

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

class VBFTruthProducer {

    private:
        VBFTagTruth truth_;

    public:
        VBFTruthProducer(){};

        ~VBFTruthProducer(){};

        void produce ( unsigned int diPhotonIndex,
                                Handle<View<reco::GenParticle> > genParticles,
                                Handle<View<reco::GenJet> > genJets,
                                Handle<View<flashgg::DiPhotonCandidate> > diPhotonCollection,
                                std::vector<edm::Handle<edm::View<flashgg::Jet> > > jetCollections );

        VBFTagTruth truthObject() {return truth_;}

        MVAVarStruct recoLevelMVAVars();
        MVAVarStruct genJetLevelMVAVars();
        MVAVarStruct genParticleLevelMVAVars();
        MVAVarStruct partonLevelMVAVars();
                


};

