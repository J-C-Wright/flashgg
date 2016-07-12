#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/EXOTag.h"
#include "flashgg/DataFormats/interface/Jet.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include <vector>
#include <algorithm>

using namespace std;
using namespace edm;

namespace flashgg {

    class EXOTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        EXOTagProducer( const ParameterSet & );

    private:
        void produce( Event &, const EventSetup & ) override;

        Int_t event_number;

        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > >  tokenJets_;
        std::vector<edm::InputTag>                           inputTagJets_;
        EDGetTokenT< edm::View<flashgg::Electron> > electronToken_;
        EDGetTokenT< edm::View<flashgg::Muon> > muonToken_;
        EDGetTokenT< edm::View<flashgg::DiPhotonCandidate> > diPhotonToken_;
        EDGetTokenT<double> rhoToken_;
        bool debug_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

    };
    

    EXOTagProducer::EXOTagProducer( const ParameterSet &iConfig ) :
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) ),
        electronToken_( consumes<View<flashgg::Electron>>(iConfig.getParameter<edm::InputTag>("ElectronTag"))),
        muonToken_( consumes<View<flashgg::Muon>>(iConfig.getParameter<edm::InputTag>("MuonTag"))),
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        rhoToken_( consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),       
        debug_( iConfig.getUntrackedParameter<bool>( "debug", false ) )
    {
        for (unsigned i = 0 ; i < inputTagJets_.size() ; i++) {
            auto token = consumes<View<flashgg::Jet> >(inputTagJets_[i]);
            tokenJets_.push_back(token);
        }

        produces<vector<EXOTag> >();

    }

    void EXOTagProducer::produce( Event &iEvent, const EventSetup & )
    {

        //DIPHOTON
        Handle<View<flashgg::DiPhotonCandidate> > diphotons;
        iEvent.getByToken( diPhotonToken_, diphotons );

        Handle<double> rhoHandle; 
        iEvent.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );

        //JETS
        JetCollectionVector Jets(inputTagJets_.size());
        for (size_t i=0;i<inputTagJets_.size();i++){
            iEvent.getByToken(tokenJets_[i],Jets[i]);
        }

        //ELECTRONS
        Handle<View<flashgg::Electron>> electrons;
        iEvent.getByToken(electronToken_,electrons);
        
        //MUONS
        Handle <View<flashgg::Muon>> muons; 
        iEvent.getByToken(muonToken_,muons); 

        std::auto_ptr<vector<EXOTag> > tags(new vector<EXOTag>);

        for( unsigned int candIndex = 0; candIndex < diphotons->size() ; candIndex++ ) {


            edm::Ptr<flashgg::DiPhotonCandidate> diphoton = diphotons->ptrAt( candIndex );
            unsigned jetCollectionIndex = diphoton->jetCollectionIndex();
            edm::Handle<edm::View<flashgg::Jet>> jets = Jets[jetCollectionIndex];

            EXOTag exoTag(diphoton,jets,electrons,muons,rhoFixedGrd,event_number);
            tags->push_back(exoTag);
        }
        iEvent.put( tags );

        event_number++;
    }
}

typedef flashgg::EXOTagProducer FlashggEXOTagProducer;
DEFINE_FWK_MODULE( FlashggEXOTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

