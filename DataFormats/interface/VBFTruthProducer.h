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

class VBFTruthProducer {

    private:

    public:
        VBFTruthProducer(){};

        ~VBFTruthProducer(){};

        VBFTagTruth produce ( unsigned int diPhotonIndex,
                                Handle<View<reco::GenParticle> > genParticles,
                                Handle<View<reco::GenJet> > genJets,
                                Handle<View<flashgg::DiPhotonCandidate> > diPhotonCollection,
                                std::vector<edm::Handle<edm::View<flashgg::Jet> > > jetCollections );



};
