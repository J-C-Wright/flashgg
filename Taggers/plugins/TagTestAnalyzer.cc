// TagTestAnalyzer.cc by S. Zenz
//
// * Tests getting tags out of the event, sorting them, casting them
// * Dumps debugging output to the screen
// * Useful for quick tests of code changes, and should be kept up-to-date as tags are added/changed
// * Should NOT be included in productions
//
// Adapted from globelikePlotMakerWithTagSorter code by L. D. Corpe, which was
// Adapted from the flashggCommissioning plot maker code  by C. Favaro et al.

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

#include "flashgg/DataFormats/interface/VBFTag.h"
#include "flashgg/DataFormats/interface/UntaggedTag.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/TTHHadronicTag.h"
#include "flashgg/DataFormats/interface/TTHLeptonicTag.h"
#include "flashgg/DataFormats/interface/VHTightTag.h"
#include "flashgg/DataFormats/interface/VHEtTag.h"
#include "flashgg/DataFormats/interface/VHLooseTag.h"
#include "flashgg/DataFormats/interface/VHHadronicTag.h"
#include "flashgg/DataFormats/interface/VBFTagTruth.h"
#include "flashgg/DataFormats/interface/VBFTruthProducer.h"
#include "flashgg/DataFormats/interface/VBFPlotProducer.h"

#include "TTree.h"

using namespace std;
using namespace edm;

// **********************************************************************

namespace flashgg {

    class TagTestAnalyzer : public edm::EDAnalyzer
    {
    public:
        explicit TagTestAnalyzer( const edm::ParameterSet & );
        ~TagTestAnalyzer();

        static void fillDescriptions( edm::ConfigurationDescriptions &descriptions );


    private:

        edm::Service<TFileService> fs_;

        virtual void beginJob() override;
        virtual void analyze( const edm::Event &, const edm::EventSetup & ) override;
        virtual void endJob() override;

        edm::EDGetTokenT<edm::OwnVector<flashgg::DiPhotonTagBase> > TagSorterToken_;
        edm::EDGetTokenT<View<reco::GenParticle> > genPartToken_;
        edm::EDGetTokenT<View<reco::GenJet> > genJetToken_;
        EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        std::vector<edm::InputTag> inputTagJets_;

        typedef std::vector<edm::Handle<edm::View<flashgg::Jet> > > JetCollectionVector;

        TFile *outputFile_;
        TTree *jjjTree;
        TTree *jjfTree;
        TTree *jffTree;
        TTree *fffTree;
        VBFTagTruth truth;
        MVAVarStruct recoLevel;
        MVAVarStruct genJetLevel;
        MVAVarStruct genParticleLevel;
        MVAVarStruct partonLevel;

    };

// ******************************************************************************************
// ******************************************************************************************

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
    TagTestAnalyzer::TagTestAnalyzer( const edm::ParameterSet &iConfig ):
        TagSorterToken_( consumes<edm::OwnVector<flashgg::DiPhotonTagBase> >( iConfig.getParameter<InputTag> ( "TagSorter" ) ) ),
        genPartToken_( consumes<View<reco::GenParticle> >( iConfig.getUntrackedParameter<InputTag> ( "GenParticleTag", InputTag( "flashggPrunedGenParticles" ) ) ) ),
        genJetToken_( consumes<View<reco::GenJet> >( iConfig.getUntrackedParameter<InputTag> ( "GenJetTag", InputTag( "slimmedGenJets" ) ) ) ),
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "inputTagJets" ) )
    {
    }

    TagTestAnalyzer::~TagTestAnalyzer()
    {

    }

    void
    TagTestAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup )
    {

        bool debug = true;
//        float dRCut = 0.5;

        // ********************************************************************************
        // access edm objects

        Handle<edm::OwnVector<flashgg::DiPhotonTagBase> > TagSorter;
        iEvent.getByToken( TagSorterToken_, TagSorter );

        Handle<View<reco::GenParticle> > genParticles;
        iEvent.getByToken( genPartToken_, genParticles );

        Handle<View<reco::GenJet> > genJets;
        iEvent.getByToken( genJetToken_, genJets );

        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        iEvent.getByToken( diPhotonToken_, diPhotons );

        JetCollectionVector Jets( inputTagJets_.size() );
        for ( size_t j=0;j < inputTagJets_.size(); j++ ) {
            iEvent.getByLabel( inputTagJets_[j], Jets[j] );
        }

        if (diPhotons->size() == 0) {std::cout << "There are no preselected diphotons!" << std::endl; return;}
        if (genParticles->size() == 0) {std::cout << "There are no GenParticles" << std::endl; return; }        
        if (genJets->size() == 0) {std::cout << "There are no GenJets" << std::endl; return; }        

        unsigned candIndex(0);
        for (unsigned int dpIndex(0);dpIndex<diPhotons->size();dpIndex++) {
            if (diPhotons->ptrAt(dpIndex)->sumPt() > diPhotons->ptrAt(candIndex)->sumPt()) {candIndex = dpIndex;}
        }
        if (Jets[diPhotons->ptrAt(candIndex)->jetCollectionIndex()]->size() == 0) {std::cout << "There are no FLASHgg jets" << std::endl; return;}

        VBFTruthProducer truthProducer;
        truthProducer.produce(candIndex,genParticles,genJets,diPhotons,Jets);
        truth = truthProducer.truthObject();

        if (!truth.hasDijet()) {std::cout << "No dijet" << std::endl; return;}
        if (debug) { 
            truth = truthProducer.truthObject();
            std::cout << setw(36) << "Jets" << setw(36) << "Parton matches" << std::endl;
            std::cout << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi" << setw(12) << "Pt" << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
            std::cout << setw(12) << truth.leadingJet()->pt() << setw(12) << truth.leadingJet()->eta() << setw(12) << truth.leadingJet()->phi();
            std::cout << setw(12) << truth.closestPartonToLeadingJet()->eta() << setw(12) << truth.closestPartonToLeadingJet()->phi();
            std::cout << setw(12) << truth.dR_partonMatchingToJ1() << std::endl;
            std::cout << setw(12) << truth.subLeadingJet()->pt() << setw(12) << truth.subLeadingJet()->eta() << setw(12) << truth.subLeadingJet()->phi();
            std::cout << setw(12) << truth.closestPartonToSubLeadingJet()->eta() << setw(12) << truth.closestPartonToSubLeadingJet()->phi();
            std::cout << setw(12) << truth.dR_partonMatchingToJ2() << std::endl;
            if (truth.hasTrijet()) {
                std::cout << setw(12) << truth.subSubLeadingJet()->pt() << setw(12) << truth.subSubLeadingJet()->eta() << setw(12) << truth.subSubLeadingJet()->phi();
                std::cout << setw(12) << truth.closestPartonToSubSubLeadingJet()->eta() << setw(12) << truth.closestPartonToSubSubLeadingJet()->phi();
                std::cout << setw(12) << truth.dR_partonMatchingToJ3() << std::endl;
            }
        }

        std::cout << "Tree test: " << std::endl;
        recoLevel = truthProducer.recoLevelMVAVars();
        std::cout << "struct done..." << std::endl;
        std::cout << recoLevel.leadingJetPt << std::endl;
        jjjTree->Fill();

/*
        genJetLevel = truthProducer.genJetLevelMVAVars();
        genParticleLevel = truthProducer.genParticleLevelMVAVars();
        partonLevel = truthProducer.partonLevelMVAVars();
        unsigned matchesPostDRCut = truth.numberOfMatchesAfterDRCut(dRCut);
        //Look at trijet candidates, classify by matching, fill trees
        if (truth.hasTrijet()) {
            if (matchesPostDRCut == 3) {
                if (debug) {std::cout << "This is a JJJ event" << std::endl;}
                jjjTree->Fill();
            } else if (matchesPostDRCut == 2) {
                if (debug) {std::cout << "This is a JJF event" << std::endl;}
                jjfTree->Fill();
            } else if (matchesPostDRCut == 1) {
                if (debug) {std::cout << "This is a JFF event" << std::endl;}
                jffTree->Fill();
            } else if (matchesPostDRCut == 0) {
                if (debug) {std::cout << "This is a FFF event" << std::endl;}
                fffTree->Fill();
            }
        } 
*/




    } // analyze

    void
    TagTestAnalyzer::beginJob()
    {
        outputFile_ = new TFile( "VBF_Output.root", "RECREATE" );

        jjjTree = new TTree("jjj","ThreeTrueJets");
        jjjTree->Branch("recoLevel",&recoLevel.leadingJetPt,"leadingJetPt/F:subLeadingJetPt/F:subSubLeadingJetPt/F");
/*
        jjjTree->Branch("genJet",&genJetLevel);
        jjjTree->Branch("genPart",&genParticleLevel);
        jjjTree->Branch("parton",&partonLevel);

        jjfTree = new TTree("jjf","ThreeTrueJets");
        jjfTree->Branch("recoJet",&recoLevel);
        jjfTree->Branch("genJet",&genJetLevel);
        jjfTree->Branch("genPart",&genParticleLevel);
        jjfTree->Branch("parton",&partonLevel);

        jffTree = new TTree("jff","ThreeTrueJets");
        jffTree->Branch("recoJet",&recoLevel);
        jffTree->Branch("genJet",&genJetLevel);
        jffTree->Branch("genPart",&genParticleLevel);
        jffTree->Branch("parton",&partonLevel);

        fffTree = new TTree("fff","ThreeTrueJets");
        fffTree->Branch("recoJet",&recoLevel);
        fffTree->Branch("genJet",&genJetLevel);
        fffTree->Branch("genPart",&genParticleLevel);
        fffTree->Branch("parton",&partonLevel);
*/
    }


/*
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
*/


    void
    TagTestAnalyzer::endJob()
    {
        outputFile_->cd();
        jjjTree->Write();
//        jjfTree->Write();
//        jffTree->Write();
//        fffTree->Write();
        outputFile_->Close();
    }

    void
    TagTestAnalyzer::fillDescriptions( edm::ConfigurationDescriptions &descriptions )
    {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault( desc );
    }

} // namespace flashgg

typedef flashgg::TagTestAnalyzer FlashggTagTestAnalyzer;
DEFINE_FWK_MODULE( FlashggTagTestAnalyzer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

