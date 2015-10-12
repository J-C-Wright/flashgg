
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

#include "flashgg/DataFormats/interface/VBFTruthProducer.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "flashgg/DataFormats/interface/Jet.h"

using namespace std;
using namespace edm;
using namespace flashgg;

VBFTagTruth VBFTruthProducer::produce(  unsigned int diPhotonIndex,
                                        Handle<View<reco::GenParticle> > genParticles,
                                        Handle<View<reco::GenJet> > genJets,
                                        Handle<View<flashgg::DiPhotonCandidate> > diPhotonCollection,
                                        std::vector<edm::Handle<edm::View<flashgg::Jet> > > jetCollections ) {

    VBFTagTruth truthObject;

    //The diphoton
    edm::Ptr<flashgg::DiPhotonCandidate> diPhoton = diPhotonCollection->ptrAt(diPhotonIndex);
    truthObject.setDiPhoton(diPhoton);

//GenPartcle level
    //Find the partons
    std::vector<edm::Ptr<reco::GenParticle>> partons;
    for (unsigned int genLoop(0);genLoop < genParticles->size();genLoop++) {
        edm::Ptr<reco::GenParticle> gp = genParticles->ptrAt(genLoop);
        bool isGluon = abs( gp->pdgId() ) < 7 && gp->numberOfMothers() == 0;
        bool isQuark = gp->pdgId() == 21 && gp->numberOfMothers() == 0;
        if (isGluon || isQuark) {
            unsigned int insertionIndex(0);
            for (unsigned int parLoop(0);parLoop<partons.size();parLoop++) {
                if (gp->pt() < partons[parLoop]->pt()) { insertionIndex = parLoop + 1; }
            }
            partons.insert( partons.begin() + insertionIndex, gp);
        }
    }

    //Merge if close
/*
    Const issues need fixing  - leaving it as unmerged now
    std::pair<unsigned,unsigned> mergerIndices(0,0);
    for (unsigned i(0);i<partons.size();i++) {
        for (unsigned j(0); j < partons.size(); j++) {
            if (i <= j) continue;
            if (deltaR(partons[i]->eta(),partons[i]->phi(),partons[j]->eta(),partons[j]->phi()) < 0.4) {
                mergerIndices.first  = i;
                mergerIndices.second = j;
            }
        }
    }
    bool areDistinct = mergerIndices.first == 0 && mergerIndices.second == 0;
    if (!areDistinct) {
        partons[mergerIndices.first]->setP4(partons[mergerIndices.first]->p4() + partons[mergerIndices.second]->p4());
        partons[mergerIndices.first]->setPdgId(999);
        partons.erase(partons.begin()+mergerIndices.second);
        if (partons[0]->pt() < partons[1]->pt()) {std::swap(partons[0],partons[1]);}
    } 
*/

    //Add to truth object
    truthObject.setPtOrderedPartons(partons);
    if (partons.size() == 1) {truthObject.setLeadingParton(partons[0]);}
    if (partons.size() == 2) {truthObject.setLeadingParton(partons[0]);truthObject.setSubLeadingParton(partons[1]);}
    if (partons.size() == 3) {truthObject.setLeadingParton(partons[0]);truthObject.setSubLeadingParton(partons[1]);truthObject.setSubSubLeadingParton(partons[2]);}

//GenJet Level
    //Pt-ordered GenJets
    std::vector<edm::Ptr<reco::GenJet>> ptOrderedGenJets;
    for( unsigned int jetLoop( 0 ); jetLoop < genJets->size(); jetLoop++ ) {
        edm::Ptr<reco::GenJet> gj = genJets->ptrAt( jetLoop );
        unsigned insertionIndex( 0 );
        for( unsigned int i( 0 ); i < ptOrderedGenJets.size(); i++ ) {
            if( gj->pt() <= ptOrderedGenJets[i]->pt() && gj != ptOrderedGenJets[i] ) { insertionIndex = i + 1; }
        }
        //Remove photons        
        float dr_leadPhoton = deltaR( gj->eta(), gj->phi(),diPhoton->leadingPhoton()->eta(),diPhoton->leadingPhoton()->phi() ); 
        float dr_subLeadPhoton = deltaR( gj->eta(), gj->phi(),diPhoton->subLeadingPhoton()->eta(),diPhoton->subLeadingPhoton()->phi() ); 
        if( dr_leadPhoton > 0.1 && dr_subLeadPhoton > 0.1 ) {
            ptOrderedGenJets.insert( ptOrderedGenJets.begin() + insertionIndex, gj );
        }
    }

    //Add to truth object
    truthObject.setPtOrderedGenJets(ptOrderedGenJets);
    if (ptOrderedGenJets.size() == 1) {truthObject.setLeadingGenJet(ptOrderedGenJets[0]);}
    if (ptOrderedGenJets.size() == 2) {truthObject.setLeadingGenJet(ptOrderedGenJets[0]);truthObject.setSubLeadingGenJet(ptOrderedGenJets[1]);}
    if (ptOrderedGenJets.size() == 3) {
        truthObject.setLeadingGenJet(ptOrderedGenJets[0]);
        truthObject.setSubLeadingGenJet(ptOrderedGenJets[1]);
        truthObject.setSubSubLeadingGenJet(ptOrderedGenJets[2]);
    }

//FLASHgg Jet Level
    //Pt-ordered fggJets
    Handle<View<flashgg::Jet>> fggJets = jetCollections[diPhoton->jetCollectionIndex()]; 
    std::vector<edm::Ptr<flashgg::Jet>> ptOrderedFggJets;
    for (unsigned jetLoop(0);jetLoop<fggJets->size();jetLoop++) {
        edm::Ptr<flashgg::Jet> jet = fggJets->ptrAt(jetLoop);
        unsigned int insertionIndex( 0 );
        for ( unsigned int i(0);i<ptOrderedFggJets.size();i++) {
            if( jet->pt() <= ptOrderedFggJets[i]->pt() && jet != ptOrderedFggJets[i] ) { insertionIndex = i + 1; }
        }
        //Remove photons and pileup
        float dr_leadPhoton = deltaR( jet->eta(), jet->phi(),diPhoton->leadingPhoton()->eta(),diPhoton->leadingPhoton()->phi() ); 
        float dr_subLeadPhoton = deltaR( jet->eta(), jet->phi(),diPhoton->subLeadingPhoton()->eta(),diPhoton->subLeadingPhoton()->phi() ); 
        bool pileupRejection = true;
        if( dr_leadPhoton > 0.1 && dr_subLeadPhoton > 0.1 ) {
            if (pileupRejection && jet->passesPuJetId(diPhoton)) { 
                ptOrderedFggJets.insert( ptOrderedFggJets.begin() + insertionIndex, jet );
            } else if (!pileupRejection) {
                ptOrderedFggJets.insert( ptOrderedFggJets.begin() + insertionIndex, jet );
            }
        }
    } 
    //Add to truth object
    truthObject.setPtOrderedFggJets(ptOrderedFggJets);

//Closest matches to FLASHgg jets and the Diphoton
    //GenParticles
        //Lead
    if (ptOrderedFggJets.size() > 0) {
        float dr(999.0);
        unsigned gpIndex(0);
        for (unsigned partLoop(0);partLoop<genParticles->size();partLoop++) {
            edm::Ptr<reco::GenParticle> particle = genParticles->ptrAt(partLoop);
            float deltaR_temp = deltaR(ptOrderedFggJets[0]->eta(),ptOrderedFggJets[0]->phi(),particle->eta(),particle->phi());
            if (deltaR_temp < dr) {dr = deltaR_temp; gpIndex = partLoop;}
        }
        truthObject.setClosestParticleToLeadingJet(genParticles->ptrAt(gpIndex));
    } 
        //Sublead
    if (ptOrderedFggJets.size() > 1) {
        float dr(999.0);
        unsigned gpIndex(0);
        for (unsigned partLoop(0);partLoop<genParticles->size();partLoop++) {
            edm::Ptr<reco::GenParticle> particle = genParticles->ptrAt(partLoop);
            float deltaR_temp = deltaR(ptOrderedFggJets[1]->eta(),ptOrderedFggJets[1]->phi(),particle->eta(),particle->phi());
            if (deltaR_temp < dr) {dr = deltaR_temp; gpIndex = partLoop;}
        }
        truthObject.setClosestParticleToSubLeadingJet(genParticles->ptrAt(gpIndex));
    } 
        //Subsublead
    if (ptOrderedFggJets.size() > 2) {
        float dr(999.0);
        unsigned gpIndex(0);
        for (unsigned partLoop(0);partLoop<genParticles->size();partLoop++) {
            edm::Ptr<reco::GenParticle> particle = genParticles->ptrAt(partLoop);
            float deltaR_temp = deltaR(ptOrderedFggJets[2]->eta(),ptOrderedFggJets[2]->phi(),particle->eta(),particle->phi());
            if (deltaR_temp < dr) {dr = deltaR_temp; gpIndex = partLoop;}
        }
        truthObject.setClosestParticleToSubSubLeadingJet(genParticles->ptrAt(gpIndex));
    } 

    //GenJets
        //Lead
    if (ptOrderedFggJets.size() > 0) {
        float dr(999.0);
        unsigned gjIndex(0);
        for (unsigned jetLoop(0);jetLoop<ptOrderedGenJets.size();jetLoop++) {
            float deltaR_temp = deltaR(ptOrderedFggJets[0]->eta(),ptOrderedFggJets[0]->phi(),ptOrderedGenJets[jetLoop]->eta(),ptOrderedGenJets[jetLoop]->phi());
            if (deltaR_temp < dr) {dr = deltaR_temp; gjIndex = jetLoop;}
        }
        truthObject.setClosestGenJetToLeadingJet( ptOrderedGenJets[gjIndex] );
    }
        //Sublead
    if (ptOrderedFggJets.size() > 1) {
        float dr(999.0);
        unsigned gjIndex(0);
        for (unsigned jetLoop(0);jetLoop<ptOrderedGenJets.size();jetLoop++) {
            float deltaR_temp = deltaR(ptOrderedFggJets[1]->eta(),ptOrderedFggJets[1]->phi(),ptOrderedGenJets[jetLoop]->eta(),ptOrderedGenJets[jetLoop]->phi());
            if (deltaR_temp < dr) {dr = deltaR_temp; gjIndex = jetLoop;}
        }
        truthObject.setClosestGenJetToSubLeadingJet( ptOrderedGenJets[gjIndex] );
    }
        //Subsublead
    if (ptOrderedFggJets.size() > 2) {
        float dr(999.0);
        unsigned gjIndex(0);
        for (unsigned jetLoop(0);jetLoop<ptOrderedGenJets.size();jetLoop++) {
            float deltaR_temp = deltaR(ptOrderedFggJets[2]->eta(),ptOrderedFggJets[2]->phi(),ptOrderedGenJets[jetLoop]->eta(),ptOrderedGenJets[jetLoop]->phi());
            if (deltaR_temp < dr) {dr = deltaR_temp; gjIndex = jetLoop;}
        }
        truthObject.setClosestGenJetToSubSubLeadingJet( ptOrderedGenJets[gjIndex] );
    }

    //Diphoton-GenParticle Matching
    float dr_leadPhoton(999.),dr_subLeadPhoton(999.);
    unsigned gpIndex1(0),gpIndex2(0);
    for (unsigned partLoop(0);partLoop<genParticles->size();partLoop++) {
        edm::Ptr<reco::GenParticle> particle = genParticles->ptrAt(partLoop);
        float deltaR_temp = deltaR(diPhoton->leadingPhoton()->eta(),diPhoton->leadingPhoton()->phi(),particle->eta(),particle->phi());
        if (deltaR_temp < dr_leadPhoton) {dr_leadPhoton = deltaR_temp; gpIndex1 = partLoop;}
        deltaR_temp = deltaR(diPhoton->subLeadingPhoton()->eta(),diPhoton->subLeadingPhoton()->phi(),particle->eta(),particle->phi());
        if (deltaR_temp < dr_subLeadPhoton) {dr_subLeadPhoton = deltaR_temp; gpIndex2 = partLoop;}
    }
    truthObject.setClosestParticleToLeadingPhoton(genParticles->ptrAt(gpIndex1));
    truthObject.setClosestParticleToSubLeadingPhoton(genParticles->ptrAt(gpIndex2));

    return truthObject;
    
}
        



 



